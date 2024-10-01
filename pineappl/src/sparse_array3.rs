//! Module containing the `SparseArray3` struct.

use ndarray::{ArrayView3, Axis};
use serde::{Deserialize, Serialize};
use std::iter;
use std::mem;
use std::ops::{Index, IndexMut, Range};
use std::slice::{Iter, IterMut};

/// Struct for a sparse three-dimensional array, which is optimized for the sparsity of
/// interpolation grids.
#[derive(Clone, Deserialize, Serialize)]
pub struct SparseArray3<T> {
    entries: Vec<T>,
    indices: Vec<(usize, usize)>,
    start: usize,
    dimensions: (usize, usize, usize),
}

// TODO: write panic messages

impl<T> Index<[usize; 3]> for SparseArray3<T> {
    type Output = T;

    fn index(&self, mut index: [usize; 3]) -> &Self::Output {
        // index too small
        assert!(index[0] >= self.start);

        let dim1 = if self.dimensions.1 > self.dimensions.2 {
            index.swap(1, 2);
            self.dimensions.2
        } else {
            self.dimensions.1
        };

        // index too large
        assert!(index[0] < (self.start + (self.indices.len() - 1) / dim1));

        // index too large
        assert!(index[1] < dim1);

        let forward = dim1 * (index[0] - self.start) + index[1];
        let indices_a = &self.indices[forward];
        let indices_b = &self.indices[forward + 1];

        let zeros_left = indices_a.0;
        let offset = indices_a.1;
        let non_zeros = indices_b.1 - offset;

        // index too small
        assert!(index[2] >= zeros_left);

        // index too large
        assert!(index[2] < (non_zeros + zeros_left));

        &self.entries[offset + (index[2] - zeros_left)]
    }
}

impl<T: Clone + Default> IndexMut<[usize; 3]> for SparseArray3<T> {
    fn index_mut(&mut self, mut index: [usize; 3]) -> &mut Self::Output {
        let dim1 = if self.dimensions.1 > self.dimensions.2 {
            index.swap(1, 2);
            self.dimensions.2
        } else {
            self.dimensions.1
        };

        let max_index0 = self.start + (self.indices.len() - 1) / dim1;

        if index[0] < self.start {
            let elements = self.start - index[0];
            self.start = index[0];
            self.indices
                .splice(0..0, iter::repeat((0, 0)).take(elements * dim1));
        } else if index[0] >= self.dimensions.0 {
            panic!();
        } else if self.entries.is_empty() || (index[0] >= max_index0) {
            let elements = if self.entries.is_empty() {
                self.start = index[0];
                1
            } else {
                index[0] - max_index0 + 1
            };

            let insert = self.indices.len() - 1;
            self.indices.splice(
                insert..insert,
                iter::repeat((0, self.indices.last().unwrap().1)).take(elements * dim1),
            );
        }

        // index too large
        assert!(index[1] < dim1);

        let forward = dim1 * (index[0] - self.start) + index[1];
        let indices_a = &self.indices[forward];
        let indices_b = &self.indices[forward + 1];

        let zeros_left = indices_a.0;
        let offset = indices_a.1;
        let non_zeros = indices_b.1 - offset;

        let elements;
        let insert;

        if index[2] < zeros_left {
            elements = zeros_left - index[2];
            insert = offset;
            self.indices[forward].0 -= elements;
        } else if index[2] >= self.dimensions.2.max(self.dimensions.1) {
            panic!();
        } else if non_zeros == 0 {
            elements = 1;
            insert = offset;
            self.indices[forward].0 = index[2];
        } else if index[2] >= (zeros_left + non_zeros) {
            elements = index[2] - (zeros_left + non_zeros) + 1;
            insert = offset + non_zeros;
        } else {
            return &mut self.entries[offset + (index[2] - zeros_left)];
        }

        self.entries
            .splice(insert..insert, iter::repeat(T::default()).take(elements));
        self.indices
            .iter_mut()
            .skip(forward + 1)
            .for_each(|ix| ix.1 += elements);

        &mut self.entries[offset + (index[2] - self.indices[forward].0)]
    }
}

/// Immutable iterator over the elements of a `SparseArray3`.
pub struct IndexedIter<'a, T> {
    entry_iter: Iter<'a, T>,
    index_iter: Iter<'a, (usize, usize)>,
    offset_a: Option<&'a (usize, usize)>,
    offset_b: Option<&'a (usize, usize)>,
    tuple: (usize, usize, usize),
    dimensions: (usize, usize, usize),
}

impl<'a, T: Copy + Default + PartialEq> Iterator for IndexedIter<'a, T> {
    type Item = ((usize, usize, usize), T);

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(element) = self.entry_iter.next() {
            let offset_a = self.offset_a.unwrap();
            let offset_b = self.offset_b.unwrap();

            if self.dimensions.1 > self.dimensions.2 {
                self.tuple.1 = self.tuple.1.max(offset_a.0);

                if self.tuple.1 >= (offset_b.1 - offset_a.1 + offset_a.0) {
                    loop {
                        self.offset_a = self.offset_b;
                        self.offset_b = self.index_iter.next();

                        let offset_a = self.offset_a.unwrap();
                        let offset_b = self.offset_b?;

                        self.tuple.2 += 1;

                        if self.tuple.2 >= self.dimensions.2 {
                            self.tuple.0 += 1;
                            self.tuple.2 = 0;
                        }

                        if (offset_b.1 - offset_a.1) != 0 {
                            self.tuple.1 = offset_a.0;
                            break;
                        }
                    }
                }

                if *element == T::default() {
                    self.tuple.1 += 1;
                    self.next()
                } else {
                    let result = Some((self.tuple, *element));
                    self.tuple.1 += 1;
                    result
                }
            } else {
                self.tuple.2 = self.tuple.2.max(offset_a.0);

                if self.tuple.2 >= (offset_b.1 - offset_a.1 + offset_a.0) {
                    loop {
                        self.offset_a = self.offset_b;
                        self.offset_b = self.index_iter.next();

                        let offset_a = self.offset_a.unwrap();
                        let offset_b = self.offset_b?;

                        self.tuple.1 += 1;

                        if self.tuple.1 >= self.dimensions.1 {
                            self.tuple.0 += 1;
                            self.tuple.1 = 0;
                        }

                        if (offset_b.1 - offset_a.1) != 0 {
                            self.tuple.2 = offset_a.0;
                            break;
                        }
                    }
                }

                if *element == T::default() {
                    self.tuple.2 += 1;
                    self.next()
                } else {
                    let result = Some((self.tuple, *element));
                    self.tuple.2 += 1;
                    result
                }
            }
        } else {
            None
        }
    }
}

impl<T: Clone + Default + PartialEq> SparseArray3<T> {
    /// Constructs a new and empty `SparseArray3` with the specified dimensions `nx`, `ny` and
    /// `nz`.
    #[must_use]
    pub fn new(nx: usize, ny: usize, nz: usize) -> Self {
        Self {
            entries: vec![],
            indices: vec![(0, 0)],
            start: 0,
            dimensions: (nx, ny, nz),
        }
    }

    /// Converts `array` into a `SparseArray3`.
    #[must_use]
    pub fn from_ndarray(array: ArrayView3<T>, xstart: usize, xsize: usize) -> Self {
        let (_, ny, nz) = array.dim();
        let array = if ny > nz {
            let mut array = array;
            array.swap_axes(1, 2);
            array
        } else {
            array
        };

        let dimensions = (xsize, ny, nz);
        let mut entries = vec![];
        let mut indices = vec![];

        let mut offset = 0;

        for array2 in array.axis_iter(Axis(0)) {
            for array1 in array2.axis_iter(Axis(0)) {
                let start = array1.iter().position(|x| *x != T::default());

                if let Some(start) = start {
                    let end = array1.iter().enumerate().skip(start).fold(
                        start,
                        |last_non_zero, (index, x)| {
                            if *x == T::default() {
                                last_non_zero
                            } else {
                                index
                            }
                        },
                    ) + 1;
                    indices.push((start, offset));
                    offset += end - start;
                    entries.splice(
                        entries.len()..entries.len(),
                        array1.iter().skip(start).take(end - start).cloned(),
                    );
                } else {
                    indices.push((0, offset));
                }
            }
        }

        indices.push((0, offset));

        Self {
            entries,
            indices,
            start: xstart,
            dimensions,
        }
    }

    /// Clear the contents of the array.
    pub fn clear(&mut self) {
        self.entries.clear();
        self.indices.clear();
        self.indices.push((0, 0));
        self.start = 0;
    }

    /// Returns the dimensions of this array.
    #[must_use]
    pub const fn dimensions(&self) -> (usize, usize, usize) {
        self.dimensions
    }

    /// Returns the overhead for storing the explicitly zero and non-zero elements.
    #[must_use]
    pub fn overhead(&self) -> usize {
        (2 * self.indices.len() * mem::size_of::<usize>()) / mem::size_of::<f64>()
    }

    /// Returns the number of default (zero) elements in this array.
    #[must_use]
    pub fn zeros(&self) -> usize {
        self.entries.iter().filter(|x| **x == T::default()).count()
    }

    /// Returns the number of non-default (non-zero) elements in this array.
    #[must_use]
    pub fn len(&self) -> usize {
        self.entries.iter().filter(|x| **x != T::default()).count()
    }

    /// Returns `true` if the array contains no element.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Return an indexed `Iterator` over the non-zero elements of this array. The iterator element
    /// type is `((usize, usize, usize), T)`.
    #[must_use]
    pub fn indexed_iter(&self) -> IndexedIter<'_, T> {
        let mut result = IndexedIter {
            entry_iter: self.entries.iter(),
            index_iter: self.indices.iter(),
            offset_a: None,
            offset_b: None,
            tuple: (self.start, 0, 0),
            dimensions: self.dimensions,
        };

        result.offset_a = result.index_iter.next();
        result.offset_b = result.index_iter.next();

        result
    }

    /// Return an iterator over the elements, including zero elements.
    pub fn iter_mut(&mut self) -> IterMut<'_, T> {
        self.entries.iter_mut()
    }

    /// Return a half-open interval of indices that are filled for the first dimension.
    #[must_use]
    pub fn x_range(&self) -> Range<usize> {
        self.start
            ..(self.start + (self.indices.len() - 1) / self.dimensions.1.min(self.dimensions.2))
    }

    /// Increase the number of entries of the x-axis by one by inserting zeros at `x`.
    pub fn increase_x_at(&mut self, x: usize) {
        let dim1 = self.dimensions.1.min(self.dimensions.2);
        let nx = (self.indices.len() - 1) / dim1;

        if x <= self.start {
            self.start += 1;
        } else if x < self.start + nx {
            let at = (x - self.start) * dim1;
            let offset = self.indices[at].1;
            self.indices
                .splice(at..at, iter::repeat((0, offset)).take(dim1));
        } else if x <= self.dimensions.0 {
            // nothing to do here
        } else {
            self.dimensions.0 = x;
        }

        self.dimensions.0 += 1;
    }

    /// Removes all elements with the specified x coordinate.
    ///
    /// # Panics
    ///
    /// TODO
    pub fn remove_x(&mut self, x: usize) {
        let dim1 = self.dimensions.1.min(self.dimensions.2);
        let nx = (self.indices.len() - 1) / dim1;

        assert!((x >= self.start) && (x < self.start + nx));

        let index_a = (x - self.start) * dim1;
        let index_b = (x - self.start + 1) * dim1;
        let offset_a = self.indices[index_a].1;
        let offset_b = self.indices[index_b].1;

        self.entries.drain(offset_a..offset_b);
        self.indices
            .iter_mut()
            .skip(index_b)
            .for_each(|o| o.1 -= offset_b - offset_a);

        if (x != self.start) && (x != (self.start + nx - 1)) {
            self.indices
                .splice(index_a..index_b, iter::repeat((0, offset_a)).take(dim1));
        } else {
            if x == self.start {
                self.start += 1;
            }

            self.indices.drain(index_a..index_b);
        }

        if self.indices.last().unwrap().1 == 0 {
            self.clear();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array3;

    #[test]
    fn index_access() {
        let mut array = SparseArray3::new(40, 50, 50);

        // after creation the array must be empty
        assert_eq!(array.x_range(), 0..0);
        assert_eq!(array.overhead(), 2);
        assert!(array.is_empty());

        // insert the first element
        array[[5, 10, 10]] = 1.0;
        assert_eq!(array[[5, 10, 10]], 1.0);
        assert_eq!(array.len(), 1);
        assert_eq!(array.zeros(), 0);
        assert_eq!(array.x_range(), 5..6);
        assert_eq!(array.overhead(), 102);
        assert!(!array.is_empty());

        // insert an element after the first one
        array[[8, 10, 10]] = 2.0;
        assert_eq!(array[[8, 10, 10]], 2.0);
        assert_eq!(array[[5, 10, 10]], 1.0);
        assert_eq!(array.len(), 2);
        assert_eq!(array.zeros(), 0);
        assert_eq!(array.x_range(), 5..9);
        assert_eq!(array.overhead(), 402);
        assert!(!array.is_empty());

        // insert an element before the first one
        array[[1, 10, 10]] = 3.0;
        assert_eq!(array[[1, 10, 10]], 3.0);
        assert_eq!(array[[8, 10, 10]], 2.0);
        assert_eq!(array[[5, 10, 10]], 1.0);
        assert_eq!(array.len(), 3);
        assert_eq!(array.zeros(), 0);
        assert_eq!(array.x_range(), 1..9);
        assert_eq!(array.overhead(), 802);
        assert!(!array.is_empty());

        array[[1, 10, 11]] = 4.0;
        assert_eq!(array[[1, 10, 11]], 4.0);
        assert_eq!(array[[1, 10, 10]], 3.0);
        assert_eq!(array[[8, 10, 10]], 2.0);
        assert_eq!(array[[5, 10, 10]], 1.0);
        assert_eq!(array.len(), 4);
        assert_eq!(array.zeros(), 0);
        assert_eq!(array.x_range(), 1..9);
        assert_eq!(array.overhead(), 802);
        assert!(!array.is_empty());

        array[[1, 10, 9]] = 5.0;
        assert_eq!(array[[1, 10, 9]], 5.0);
        assert_eq!(array[[1, 10, 11]], 4.0);
        assert_eq!(array[[1, 10, 10]], 3.0);
        assert_eq!(array[[8, 10, 10]], 2.0);
        assert_eq!(array[[5, 10, 10]], 1.0);
        assert_eq!(array.len(), 5);
        assert_eq!(array.zeros(), 0);
        assert_eq!(array.x_range(), 1..9);
        assert_eq!(array.overhead(), 802);
        assert!(!array.is_empty());

        array[[1, 10, 0]] = 6.0;
        assert_eq!(array[[1, 10, 0]], 6.0);
        assert_eq!(array[[1, 10, 9]], 5.0);
        assert_eq!(array[[1, 10, 11]], 4.0);
        assert_eq!(array[[1, 10, 10]], 3.0);
        assert_eq!(array[[8, 10, 10]], 2.0);
        assert_eq!(array[[5, 10, 10]], 1.0);
        assert_eq!(array.len(), 6);
        assert_eq!(array.x_range(), 1..9);
        assert_eq!(array.overhead(), 802);
        assert!(!array.is_empty());

        // check zeros
        assert_eq!(array[[1, 10, 1]], 0.0);
        assert_eq!(array[[1, 10, 2]], 0.0);
        assert_eq!(array[[1, 10, 3]], 0.0);
        assert_eq!(array[[1, 10, 4]], 0.0);
        assert_eq!(array[[1, 10, 5]], 0.0);
        assert_eq!(array[[1, 10, 6]], 0.0);
        assert_eq!(array[[1, 10, 7]], 0.0);
        assert_eq!(array[[1, 10, 8]], 0.0);
        assert_eq!(array.zeros(), 8);

        // insert where previously a zero was
        array[[1, 10, 2]] = 7.0;
        assert_eq!(array[[1, 10, 2]], 7.0);
        assert_eq!(array[[1, 10, 0]], 6.0);
        assert_eq!(array[[1, 10, 9]], 5.0);
        assert_eq!(array[[1, 10, 11]], 4.0);
        assert_eq!(array[[1, 10, 10]], 3.0);
        assert_eq!(array[[8, 10, 10]], 2.0);
        assert_eq!(array[[5, 10, 10]], 1.0);
        assert_eq!(array.len(), 7);
        assert_eq!(array.x_range(), 1..9);
        assert_eq!(array.overhead(), 802);
        assert!(!array.is_empty());

        // check zeros
        assert_eq!(array[[1, 10, 1]], 0.0);
        assert_eq!(array[[1, 10, 3]], 0.0);
        assert_eq!(array[[1, 10, 4]], 0.0);
        assert_eq!(array[[1, 10, 5]], 0.0);
        assert_eq!(array[[1, 10, 6]], 0.0);
        assert_eq!(array[[1, 10, 7]], 0.0);
        assert_eq!(array[[1, 10, 8]], 0.0);
        assert_eq!(array.zeros(), 7);

        array[[1, 15, 2]] = 8.0;
        assert_eq!(array[[1, 15, 2]], 8.0);
        assert_eq!(array[[1, 10, 2]], 7.0);
        assert_eq!(array[[1, 10, 0]], 6.0);
        assert_eq!(array[[1, 10, 9]], 5.0);
        assert_eq!(array[[1, 10, 11]], 4.0);
        assert_eq!(array[[1, 10, 10]], 3.0);
        assert_eq!(array[[8, 10, 10]], 2.0);
        assert_eq!(array[[5, 10, 10]], 1.0);
        assert_eq!(array.len(), 8);
        assert_eq!(array.x_range(), 1..9);
        assert_eq!(array.overhead(), 802);
        assert!(!array.is_empty());

        // check zeros
        assert_eq!(array[[1, 10, 1]], 0.0);
        assert_eq!(array[[1, 10, 3]], 0.0);
        assert_eq!(array[[1, 10, 4]], 0.0);
        assert_eq!(array[[1, 10, 5]], 0.0);
        assert_eq!(array[[1, 10, 6]], 0.0);
        assert_eq!(array[[1, 10, 7]], 0.0);
        assert_eq!(array[[1, 10, 8]], 0.0);
        assert_eq!(array.zeros(), 7);

        array[[1, 15, 4]] = 9.0;
        assert_eq!(array[[1, 15, 4]], 9.0);
        assert_eq!(array[[1, 15, 2]], 8.0);
        assert_eq!(array[[1, 10, 2]], 7.0);
        assert_eq!(array[[1, 10, 0]], 6.0);
        assert_eq!(array[[1, 10, 9]], 5.0);
        assert_eq!(array[[1, 10, 11]], 4.0);
        assert_eq!(array[[1, 10, 10]], 3.0);
        assert_eq!(array[[8, 10, 10]], 2.0);
        assert_eq!(array[[5, 10, 10]], 1.0);
        assert_eq!(array.len(), 9);
        assert_eq!(array.x_range(), 1..9);
        assert_eq!(array.overhead(), 802);
        assert!(!array.is_empty());

        // check zeros
        assert_eq!(array[[1, 15, 3]], 0.0);
        assert_eq!(array[[1, 10, 1]], 0.0);
        assert_eq!(array[[1, 10, 3]], 0.0);
        assert_eq!(array[[1, 10, 4]], 0.0);
        assert_eq!(array[[1, 10, 5]], 0.0);
        assert_eq!(array[[1, 10, 6]], 0.0);
        assert_eq!(array[[1, 10, 7]], 0.0);
        assert_eq!(array[[1, 10, 8]], 0.0);
        assert_eq!(array.zeros(), 8);

        array[[1, 15, 0]] = 10.0;
        assert_eq!(array[[1, 15, 0]], 10.0);
        assert_eq!(array[[1, 15, 4]], 9.0);
        assert_eq!(array[[1, 15, 2]], 8.0);
        assert_eq!(array[[1, 10, 2]], 7.0);
        assert_eq!(array[[1, 10, 0]], 6.0);
        assert_eq!(array[[1, 10, 9]], 5.0);
        assert_eq!(array[[1, 10, 11]], 4.0);
        assert_eq!(array[[1, 10, 10]], 3.0);
        assert_eq!(array[[8, 10, 10]], 2.0);
        assert_eq!(array[[5, 10, 10]], 1.0);
        assert_eq!(array.len(), 10);
        assert_eq!(array.x_range(), 1..9);
        assert_eq!(array.overhead(), 802);
        assert!(!array.is_empty());

        // check zeros
        assert_eq!(array[[1, 15, 1]], 0.0);
        assert_eq!(array[[1, 15, 3]], 0.0);
        assert_eq!(array[[1, 10, 1]], 0.0);
        assert_eq!(array[[1, 10, 3]], 0.0);
        assert_eq!(array[[1, 10, 4]], 0.0);
        assert_eq!(array[[1, 10, 5]], 0.0);
        assert_eq!(array[[1, 10, 6]], 0.0);
        assert_eq!(array[[1, 10, 7]], 0.0);
        assert_eq!(array[[1, 10, 8]], 0.0);
        assert_eq!(array.zeros(), 9);
    }

    #[test]
    #[should_panic(expected = "explicit panic")]
    fn index_mut_panic_dim0() {
        let mut array = SparseArray3::new(40, 50, 50);

        array[[40, 0, 50]] = 1.0;
    }

    #[test]
    #[should_panic(expected = "assertion failed: index[1] < dim1")]
    fn index_mut_panic_dim1() {
        let mut array = SparseArray3::new(40, 50, 50);

        array[[0, 50, 0]] = 1.0;
    }

    #[test]
    #[should_panic(expected = "explicit panic")]
    fn index_mut_panic_dim2() {
        let mut array = SparseArray3::new(40, 50, 50);

        array[[0, 0, 50]] = 1.0;
    }

    #[test]
    #[should_panic(expected = "assertion failed: index[0] >= self.start")]
    fn index_panic_dim0_0() {
        let mut array = SparseArray3::new(40, 50, 50);

        array[[1, 0, 0]] = 1.0;

        assert_eq!(array[[0, 0, 0]], 0.0);
    }

    #[test]
    #[should_panic(
        expected = "assertion failed: index[0] < (self.start + (self.indices.len() - 1) / dim1)"
    )]
    fn index_panic_dim0_1() {
        let mut array = SparseArray3::new(40, 50, 50);

        array[[1, 0, 0]] = 1.0;

        assert_eq!(array[[2, 0, 0]], 0.0);
    }

    #[test]
    #[should_panic(expected = "assertion failed: index[1] < dim1")]
    fn index_panic_dim1() {
        let mut array = SparseArray3::new(40, 50, 50);

        array[[1, 0, 0]] = 1.0;

        assert_eq!(array[[1, 50, 0]], 0.0);
    }

    #[test]
    #[should_panic(expected = "assertion failed: index[2] >= zeros_left")]
    fn index_panic_dim2_0() {
        let mut array = SparseArray3::new(40, 50, 50);

        array[[0, 0, 1]] = 1.0;

        assert_eq!(array[[0, 0, 0]], 0.0);
    }

    #[test]
    #[should_panic(expected = "assertion failed: index[2] < (non_zeros + zeros_left)")]
    fn index_panic_dim2_1() {
        let mut array = SparseArray3::new(40, 50, 50);

        array[[0, 0, 1]] = 1.0;

        assert_eq!(array[[0, 0, 2]], 0.0);
    }

    #[test]
    fn indexed_iter() {
        let mut array = SparseArray3::new(40, 50, 50);

        // check empty iterator
        assert_eq!(array.indexed_iter().next(), None);

        // insert an element
        array[[2, 3, 4]] = 1.0;

        let mut iter = array.indexed_iter();

        // check iterator with one element
        assert_eq!(iter.next(), Some(((2, 3, 4), 1.0)));
        assert_eq!(iter.next(), None);

        // insert another element
        array[[2, 3, 6]] = 2.0;

        let mut iter = array.indexed_iter();

        assert_eq!(iter.next(), Some(((2, 3, 4), 1.0)));
        assert_eq!(iter.next(), Some(((2, 3, 6), 2.0)));
        assert_eq!(iter.next(), None);

        // insert yet another element
        array[[4, 5, 7]] = 3.0;

        let mut iter = array.indexed_iter();

        assert_eq!(iter.next(), Some(((2, 3, 4), 1.0)));
        assert_eq!(iter.next(), Some(((2, 3, 6), 2.0)));
        assert_eq!(iter.next(), Some(((4, 5, 7), 3.0)));
        assert_eq!(iter.next(), None);

        // insert at the very first position
        array[[2, 0, 0]] = 4.0;

        let mut iter = array.indexed_iter();

        assert_eq!(iter.next(), Some(((2, 0, 0), 4.0)));
        assert_eq!(iter.next(), Some(((2, 3, 4), 1.0)));
        assert_eq!(iter.next(), Some(((2, 3, 6), 2.0)));
        assert_eq!(iter.next(), Some(((4, 5, 7), 3.0)));
        assert_eq!(iter.next(), None);
    }

    #[test]
    fn iter_mut() {
        let mut array = SparseArray3::new(40, 50, 50);

        array[[3, 5, 1]] = 1.0;
        array[[7, 8, 9]] = 2.0;
        array[[7, 8, 13]] = 3.0;
        array[[9, 1, 4]] = 4.0;

        let mut iter = array.iter_mut();

        assert_eq!(iter.next(), Some(&mut 1.0));
        assert_eq!(iter.next(), Some(&mut 2.0));
        assert_eq!(iter.next(), Some(&mut 0.0));
        assert_eq!(iter.next(), Some(&mut 0.0));
        assert_eq!(iter.next(), Some(&mut 0.0));
        assert_eq!(iter.next(), Some(&mut 3.0));
        assert_eq!(iter.next(), Some(&mut 4.0));
        assert_eq!(iter.next(), None);
    }

    #[test]
    fn clear() {
        let mut array = SparseArray3::new(40, 50, 50);

        array[[3, 5, 1]] = 1.0;
        array[[7, 8, 9]] = 2.0;
        array[[9, 1, 4]] = 3.0;

        assert!(!array.is_empty());
        assert_eq!(array.len(), 3);
        assert_eq!(array.zeros(), 0);
        assert_eq!(array.x_range(), 3..10);

        array.clear();

        assert!(array.is_empty());
        assert_eq!(array.len(), 0);
        assert_eq!(array.zeros(), 0);
        assert_eq!(array.x_range(), 0..0);
    }

    #[test]
    fn remove_x() {
        let mut array = SparseArray3::new(40, 50, 50);

        array[[1, 5, 6]] = 1.0;
        array[[1, 6, 5]] = 2.0;
        array[[1, 2, 3]] = 3.0;
        array[[1, 9, 3]] = 4.0;
        array[[1, 8, 4]] = 5.0;
        array[[2, 0, 0]] = 6.0;
        array[[3, 4, 5]] = 7.0;
        array[[3, 4, 6]] = 8.0;
        array[[3, 4, 7]] = 9.0;
        array[[4, 0, 2]] = 10.0;
        array[[4, 0, 3]] = 11.0;
        array[[5, 0, 1]] = 12.0;
        array[[5, 0, 2]] = 13.0;

        assert_eq!(array.x_range(), 1..6);
        assert_eq!(array.len(), 13);
        assert_eq!(array.zeros(), 0);

        // remove the first five entries
        array.remove_x(1);

        assert_eq!(array.x_range(), 2..6);
        assert_eq!(array.len(), 8);
        assert_eq!(array.zeros(), 0);

        // remove the last two entries
        array.remove_x(5);

        assert_eq!(array.x_range(), 2..5);
        assert_eq!(array.len(), 6);
        assert_eq!(array.zeros(), 0);

        // remove the from the middle
        array.remove_x(3);

        assert_eq!(array.x_range(), 2..5);
        assert_eq!(array.len(), 3);
        assert_eq!(array.zeros(), 0);

        // remove also the rest
        array.remove_x(4);
        array.remove_x(2);

        assert_eq!(array.x_range(), 0..0);
        assert_eq!(array.len(), 0);
        assert_eq!(array.zeros(), 0);
    }

    #[test]
    #[should_panic(expected = "assertion failed: (x >= self.start) && (x < self.start + nx)")]
    fn remove_x_panic() {
        let mut array = SparseArray3::<f64>::new(40, 50, 50);

        array.remove_x(0);
    }

    #[test]
    fn increase_at_x() {
        let mut array = SparseArray3::new(1, 50, 50);

        array[[0, 0, 0]] = 1.0;
        array[[0, 2, 3]] = 2.0;
        array[[0, 2, 4]] = 3.0;
        array[[0, 2, 5]] = 4.0;
        array[[0, 3, 0]] = 5.0;
        array[[0, 49, 49]] = 6.0;

        assert_eq!(array.dimensions(), (1, 50, 50));
        assert_eq!(array[[0, 0, 0]], 1.0);
        assert_eq!(array[[0, 2, 3]], 2.0);
        assert_eq!(array[[0, 2, 4]], 3.0);
        assert_eq!(array[[0, 2, 5]], 4.0);
        assert_eq!(array[[0, 3, 0]], 5.0);
        assert_eq!(array[[0, 49, 49]], 6.0);

        // increase at the end
        array.increase_x_at(1);

        assert_eq!(array.dimensions(), (2, 50, 50));
        assert_eq!(array[[0, 0, 0]], 1.0);
        assert_eq!(array[[0, 2, 3]], 2.0);
        assert_eq!(array[[0, 2, 4]], 3.0);
        assert_eq!(array[[0, 2, 5]], 4.0);
        assert_eq!(array[[0, 3, 0]], 5.0);
        assert_eq!(array[[0, 49, 49]], 6.0);

        array[[1, 5, 0]] = 7.0;
        array[[1, 5, 5]] = 8.0;
        array[[1, 6, 3]] = 9.0;
        array[[1, 6, 0]] = 10.0;

        assert_eq!(array[[0, 0, 0]], 1.0);
        assert_eq!(array[[0, 2, 3]], 2.0);
        assert_eq!(array[[0, 2, 4]], 3.0);
        assert_eq!(array[[0, 2, 5]], 4.0);
        assert_eq!(array[[0, 3, 0]], 5.0);
        assert_eq!(array[[0, 49, 49]], 6.0);
        assert_eq!(array[[1, 5, 0]], 7.0);
        assert_eq!(array[[1, 5, 5]], 8.0);
        assert_eq!(array[[1, 6, 3]], 9.0);
        assert_eq!(array[[1, 6, 0]], 10.0);

        // increase at the start
        array.increase_x_at(0);

        assert_eq!(array.dimensions(), (3, 50, 50));
        assert_eq!(array[[1, 0, 0]], 1.0);
        assert_eq!(array[[1, 2, 3]], 2.0);
        assert_eq!(array[[1, 2, 4]], 3.0);
        assert_eq!(array[[1, 2, 5]], 4.0);
        assert_eq!(array[[1, 3, 0]], 5.0);
        assert_eq!(array[[1, 49, 49]], 6.0);
        assert_eq!(array[[2, 5, 0]], 7.0);
        assert_eq!(array[[2, 5, 5]], 8.0);
        assert_eq!(array[[2, 6, 3]], 9.0);
        assert_eq!(array[[2, 6, 0]], 10.0);

        // increase at the end
        array.increase_x_at(3);

        assert_eq!(array.dimensions(), (4, 50, 50));
        assert_eq!(array[[1, 0, 0]], 1.0);
        assert_eq!(array[[1, 2, 3]], 2.0);
        assert_eq!(array[[1, 2, 4]], 3.0);
        assert_eq!(array[[1, 2, 5]], 4.0);
        assert_eq!(array[[1, 3, 0]], 5.0);
        assert_eq!(array[[1, 49, 49]], 6.0);
        assert_eq!(array[[2, 5, 0]], 7.0);
        assert_eq!(array[[2, 5, 5]], 8.0);
        assert_eq!(array[[2, 6, 3]], 9.0);
        assert_eq!(array[[2, 6, 0]], 10.0);

        // increase after the end
        array.increase_x_at(5);

        assert_eq!(array.dimensions(), (6, 50, 50));
        assert_eq!(array[[1, 0, 0]], 1.0);
        assert_eq!(array[[1, 2, 3]], 2.0);
        assert_eq!(array[[1, 2, 4]], 3.0);
        assert_eq!(array[[1, 2, 5]], 4.0);
        assert_eq!(array[[1, 3, 0]], 5.0);
        assert_eq!(array[[1, 49, 49]], 6.0);
        assert_eq!(array[[2, 5, 0]], 7.0);
        assert_eq!(array[[2, 5, 5]], 8.0);
        assert_eq!(array[[2, 6, 3]], 9.0);
        assert_eq!(array[[2, 6, 0]], 10.0);

        // increase in the middle
        array.increase_x_at(2);

        assert_eq!(array.dimensions(), (7, 50, 50));
        assert_eq!(array[[1, 0, 0]], 1.0);
        assert_eq!(array[[1, 2, 3]], 2.0);
        assert_eq!(array[[1, 2, 4]], 3.0);
        assert_eq!(array[[1, 2, 5]], 4.0);
        assert_eq!(array[[1, 3, 0]], 5.0);
        assert_eq!(array[[1, 49, 49]], 6.0);
        assert_eq!(array[[3, 5, 0]], 7.0);
        assert_eq!(array[[3, 5, 5]], 8.0);
        assert_eq!(array[[3, 6, 3]], 9.0);
        assert_eq!(array[[3, 6, 0]], 10.0);
    }

    #[test]
    fn from_ndarray() {
        let mut ndarray = Array3::zeros((2, 50, 50));

        ndarray[[0, 4, 3]] = 1.0;
        ndarray[[0, 4, 4]] = 2.0;
        ndarray[[0, 4, 6]] = 3.0;
        ndarray[[0, 5, 1]] = 4.0;
        ndarray[[0, 5, 7]] = 5.0;
        ndarray[[1, 3, 9]] = 6.0;

        let array = SparseArray3::from_ndarray(ndarray.view(), 3, 40);

        assert_eq!(array[[3, 4, 3]], 1.0);
        assert_eq!(array[[3, 4, 4]], 2.0);
        assert_eq!(array[[3, 4, 5]], 0.0);
        assert_eq!(array[[3, 4, 6]], 3.0);
        assert_eq!(array[[3, 5, 1]], 4.0);
        assert_eq!(array[[3, 5, 2]], 0.0);
        assert_eq!(array[[3, 5, 3]], 0.0);
        assert_eq!(array[[3, 5, 4]], 0.0);
        assert_eq!(array[[3, 5, 5]], 0.0);
        assert_eq!(array[[3, 5, 6]], 0.0);
        assert_eq!(array[[3, 5, 7]], 5.0);
        assert_eq!(array[[4, 3, 9]], 6.0);

        assert_eq!(array.len(), 6);
        assert_eq!(array.zeros(), 6);
    }

    #[test]
    fn test_index_swap() {
        let mut array = SparseArray3::new(5, 50, 2);

        array[[0, 0, 0]] = 1.0;
        array[[0, 0, 1]] = 2.0;
        array[[1, 2, 1]] = 3.0;
        array[[1, 5, 1]] = 4.0;
        array[[1, 6, 0]] = 5.0;
        array[[1, 8, 0]] = 6.0;
        array[[1, 9, 0]] = 7.0;
        array[[2, 0, 0]] = 8.0;
        array[[3, 2, 1]] = 9.0;
        array[[3, 4, 0]] = 10.0;
        array[[3, 4, 1]] = 11.0;
        array[[4, 0, 0]] = 12.0;
        array[[4, 0, 1]] = 13.0;

        assert_eq!(array[[0, 0, 0]], 1.0);
        assert_eq!(array[[0, 0, 1]], 2.0);
        assert_eq!(array[[1, 2, 1]], 3.0);
        assert_eq!(array[[1, 5, 1]], 4.0);
        assert_eq!(array[[1, 6, 0]], 5.0);
        assert_eq!(array[[1, 8, 0]], 6.0);
        assert_eq!(array[[1, 9, 0]], 7.0);
        assert_eq!(array[[2, 0, 0]], 8.0);
        assert_eq!(array[[3, 2, 1]], 9.0);
        assert_eq!(array[[3, 4, 0]], 10.0);
        assert_eq!(array[[3, 4, 1]], 11.0);
        assert_eq!(array[[4, 0, 0]], 12.0);
        assert_eq!(array[[4, 0, 1]], 13.0);

        assert_eq!(array.x_range(), 0..5);

        let mut iter = array.indexed_iter();

        assert_eq!(iter.next(), Some(((0, 0, 0), 1.0)));
        assert_eq!(iter.next(), Some(((0, 0, 1), 2.0)));
        assert_eq!(iter.next(), Some(((1, 6, 0), 5.0)));
        assert_eq!(iter.next(), Some(((1, 8, 0), 6.0)));
        assert_eq!(iter.next(), Some(((1, 9, 0), 7.0)));
        assert_eq!(iter.next(), Some(((1, 2, 1), 3.0)));
        assert_eq!(iter.next(), Some(((1, 5, 1), 4.0)));
        assert_eq!(iter.next(), Some(((2, 0, 0), 8.0)));
        assert_eq!(iter.next(), Some(((3, 4, 0), 10.0)));
        assert_eq!(iter.next(), Some(((3, 2, 1), 9.0)));
        assert_eq!(iter.next(), Some(((3, 4, 1), 11.0)));
        assert_eq!(iter.next(), Some(((4, 0, 0), 12.0)));
        assert_eq!(iter.next(), Some(((4, 0, 1), 13.0)));
        assert_eq!(iter.next(), None);

        let mut ndarray = Array3::zeros((5, 50, 2));

        ndarray[[0, 0, 0]] = 1.0;
        ndarray[[0, 0, 1]] = 2.0;
        ndarray[[1, 2, 1]] = 3.0;
        ndarray[[1, 5, 1]] = 4.0;
        ndarray[[1, 6, 0]] = 5.0;
        ndarray[[1, 8, 0]] = 6.0;
        ndarray[[1, 9, 0]] = 7.0;
        ndarray[[2, 0, 0]] = 8.0;
        ndarray[[3, 2, 1]] = 9.0;
        ndarray[[3, 4, 0]] = 10.0;
        ndarray[[3, 4, 1]] = 11.0;
        ndarray[[4, 0, 0]] = 12.0;
        ndarray[[4, 0, 1]] = 13.0;

        let mut other = SparseArray3::from_ndarray(ndarray.view(), 0, 5);

        assert_eq!(other[[0, 0, 0]], 1.0);
        assert_eq!(other[[0, 0, 1]], 2.0);
        assert_eq!(other[[1, 2, 1]], 3.0);
        assert_eq!(other[[1, 5, 1]], 4.0);
        assert_eq!(other[[1, 6, 0]], 5.0);
        assert_eq!(other[[1, 8, 0]], 6.0);
        assert_eq!(other[[1, 9, 0]], 7.0);
        assert_eq!(other[[2, 0, 0]], 8.0);
        assert_eq!(other[[3, 2, 1]], 9.0);
        assert_eq!(other[[3, 4, 0]], 10.0);
        assert_eq!(other[[3, 4, 1]], 11.0);
        assert_eq!(other[[4, 0, 0]], 12.0);
        assert_eq!(other[[4, 0, 1]], 13.0);

        assert_eq!(other.x_range(), 0..5);

        other.remove_x(0);

        assert_eq!(other[[1, 2, 1]], 3.0);
        assert_eq!(other[[1, 5, 1]], 4.0);
        assert_eq!(other[[1, 6, 0]], 5.0);
        assert_eq!(other[[1, 8, 0]], 6.0);
        assert_eq!(other[[1, 9, 0]], 7.0);
        assert_eq!(other[[2, 0, 0]], 8.0);
        assert_eq!(other[[3, 2, 1]], 9.0);
        assert_eq!(other[[3, 4, 0]], 10.0);
        assert_eq!(other[[3, 4, 1]], 11.0);
        assert_eq!(other[[4, 0, 0]], 12.0);
        assert_eq!(other[[4, 0, 1]], 13.0);

        other.remove_x(3);

        assert_eq!(other[[1, 2, 1]], 3.0);
        assert_eq!(other[[1, 5, 1]], 4.0);
        assert_eq!(other[[1, 6, 0]], 5.0);
        assert_eq!(other[[1, 8, 0]], 6.0);
        assert_eq!(other[[1, 9, 0]], 7.0);
        assert_eq!(other[[2, 0, 0]], 8.0);
        assert_eq!(other[[4, 0, 0]], 12.0);
        assert_eq!(other[[4, 0, 1]], 13.0);

        other.remove_x(4);

        assert_eq!(other[[1, 2, 1]], 3.0);
        assert_eq!(other[[1, 5, 1]], 4.0);
        assert_eq!(other[[1, 6, 0]], 5.0);
        assert_eq!(other[[1, 8, 0]], 6.0);
        assert_eq!(other[[1, 9, 0]], 7.0);
        assert_eq!(other[[2, 0, 0]], 8.0);
    }

    // https://github.com/NNPDF/pineappl/issues/220
    #[test]
    fn regression_test_220() {
        let mut array = SparseArray3::new(1, 2, 4);

        array[[0, 0, 0]] = 1.0;

        assert_eq!(array[[0, 0, 0]], 1.0);

        assert_eq!(array.x_range(), 0..1);

        let mut iter = array.indexed_iter();

        assert_eq!(iter.next(), Some(((0, 0, 0), 1.0)));
        assert_eq!(iter.next(), None);

        array.increase_x_at(0);

        array[[0, 0, 0]] = 2.0;

        let mut iter = array.indexed_iter();

        assert_eq!(iter.next(), Some(((0, 0, 0), 2.0)));
        assert_eq!(iter.next(), Some(((1, 0, 0), 1.0)));
        assert_eq!(iter.next(), None);

        array.increase_x_at(1);

        array[[1, 0, 0]] = 3.0;

        let mut iter = array.indexed_iter();

        assert_eq!(iter.next(), Some(((0, 0, 0), 2.0)));
        assert_eq!(iter.next(), Some(((1, 0, 0), 3.0)));
        assert_eq!(iter.next(), Some(((2, 0, 0), 1.0)));
        assert_eq!(iter.next(), None);
    }
}
