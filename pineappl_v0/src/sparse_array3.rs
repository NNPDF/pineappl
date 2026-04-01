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
