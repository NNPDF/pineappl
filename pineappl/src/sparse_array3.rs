use std::iter;
use std::ops::{Index, IndexMut};
use std::slice::Iter;

#[derive(Clone)]
struct Offset {
    space: usize,
    offset: usize,
}

/// Struct for a sparse three-dimensional array, which is optimized for the sparsity of
/// interpolation grids.
pub struct SparseArray3<T> {
    entries: Vec<T>,
    indices: Vec<Offset>,
    start: usize,
    dimensions: (usize, usize, usize),
}

// TODO: write panic messages

impl<T> Index<[usize; 3]> for SparseArray3<T> {
    type Output = T;

    fn index(&self, index: [usize; 3]) -> &Self::Output {
        // index too small
        if index[0] < self.start {
            panic!();
        }

        // index too large
        if index[0] >= (self.start + (self.indices.len() - 1) / self.dimensions.1) {
            panic!();
        }

        // index too large
        if index[1] >= self.dimensions.1 {
            panic!();
        }

        let forward = self.dimensions.1 * (index[0] - self.start) + index[1];
        let indices_a = &self.indices[forward];
        let indices_b = &self.indices[forward + 1];

        let zeros_left = indices_a.space;
        let offset = indices_a.offset;
        let non_zeros = indices_b.offset - offset;

        // index too small
        if index[2] < zeros_left {
            panic!();
        }

        // index too large
        if index[2] >= (non_zeros + zeros_left) {
            panic!();
        }

        &self.entries[offset + (index[2] - zeros_left)]
    }
}

impl<T: Clone + Default> IndexMut<[usize; 3]> for SparseArray3<T> {
    fn index_mut(&mut self, index: [usize; 3]) -> &mut Self::Output {
        let max_index0 = self.start + (self.indices.len() - 1) / self.dimensions.1;

        if index[0] < self.start {
            let elements = self.start - index[0];
            self.start = index[0];
            self.indices.splice(
                0..0,
                iter::repeat(Offset {
                    space: 0,
                    offset: 0,
                })
                .take(elements * self.dimensions.1),
            );
        } else if index[0] >= self.dimensions.0 {
            panic!();
        } else if self.entries.is_empty() || (index[0] >= max_index0) {
            let elements;

            if self.entries.is_empty() {
                self.start = index[0];
                elements = 1;
            } else {
                elements = index[0] - max_index0 + 1;
            }

            let insert = self.indices.len() - 1;
            self.indices.splice(
                insert..insert,
                iter::repeat(Offset {
                    space: 0,
                    offset: self.indices.last().unwrap().offset,
                })
                .take(elements * self.dimensions.1),
            );
        }

        // index too large
        if index[1] >= self.dimensions.1 {
            panic!();
        }

        let forward = self.dimensions.1 * (index[0] - self.start) + index[1];
        let indices_a = &self.indices[forward];
        let indices_b = &self.indices[forward + 1];

        let zeros_left = indices_a.space;
        let offset = indices_a.offset;
        let non_zeros = indices_b.offset - offset;

        let elements;
        let insert;

        if index[2] < zeros_left {
            elements = zeros_left - index[2];
            insert = offset;
            self.indices[forward].space -= elements;
        } else if index[2] >= self.dimensions.2 {
            panic!();
        } else if non_zeros == 0 {
            elements = 1;
            insert = offset;
            self.indices[forward].space = index[2];
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
            .for_each(|ix| ix.offset += elements);

        &mut self.entries[offset + (index[2] - self.indices[forward].space)]
    }
}

pub struct SparseArray3Iter<'a, T> {
    entry_iter: Iter<'a, T>,
    index_iter: Iter<'a, Offset>,
    offset_a: Option<&'a Offset>,
    offset_b: Option<&'a Offset>,
    tuple: (usize, usize, usize),
    dimensions: (usize, usize, usize),
}

impl<'a, T: Default + PartialEq> Iterator for SparseArray3Iter<'a, T> {
    type Item = ((usize, usize, usize), &'a T);

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(element) = self.entry_iter.next() {
            self.tuple.2 += 1;

            let offset_a = self.offset_a.unwrap();
            let offset_b = self.offset_b.unwrap();

            if self.tuple.2 >= (offset_b.offset - offset_a.offset + offset_a.space) {
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

                    if (offset_b.offset - offset_a.offset) != 0 {
                        self.tuple.2 = offset_a.space;
                        break;
                    }
                }
            }

            if *element == T::default() {
                self.next()
            } else {
                Some((self.tuple, element))
            }
        } else {
            None
        }
    }
}

impl<'a, T> SparseArray3Iter<'a, T> {
    fn new(array: &'a SparseArray3<T>) -> Self {
        let mut result = Self {
            entry_iter: array.entries.iter(),
            index_iter: array.indices.iter(),
            offset_a: None,
            offset_b: None,
            tuple: (array.start, 0, 0),
            dimensions: array.dimensions,
        };

        result.offset_a = result.index_iter.next();
        result.offset_b = result.index_iter.next();

        result
    }
}

// TODO: implement conversion from Array3 to SparseArray3

impl<T: Default + PartialEq> SparseArray3<T> {
    /// Constructs a new and empty `SparseArray3` with the specified dimensions `nx`, `ny` and
    /// `nz`.
    pub fn new(nx: usize, ny: usize, nz: usize) -> Self {
        Self {
            entries: vec![],
            indices: vec![Offset {
                space: 0,
                offset: 0,
            }],
            start: 0,
            dimensions: (nx, ny, nz),
        }
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

    /// Return an `Iterator` over the elements of this array.
    pub fn iter<'a>(&'a self) -> SparseArray3Iter<'a, T> {
        SparseArray3Iter::new(&self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn index_access() {
        let mut array = SparseArray3::new(40, 50, 50);

        // after creation the array must be empty
        assert!(array.is_empty());

        // insert the first element
        array[[5, 10, 10]] = 1.0;
        assert_eq!(array[[5, 10, 10]], 1.0);
        assert_eq!(array.len(), 1);
        assert_eq!(array.zeros(), 0);
        assert!(!array.is_empty());

        // insert an element after the first one
        array[[8, 10, 10]] = 2.0;
        assert_eq!(array[[8, 10, 10]], 2.0);
        assert_eq!(array[[5, 10, 10]], 1.0);
        assert_eq!(array.len(), 2);
        assert_eq!(array.zeros(), 0);
        assert!(!array.is_empty());

        // insert an element before the first one
        array[[1, 10, 10]] = 3.0;
        assert_eq!(array[[1, 10, 10]], 3.0);
        assert_eq!(array[[8, 10, 10]], 2.0);
        assert_eq!(array[[5, 10, 10]], 1.0);
        assert_eq!(array.len(), 3);
        assert_eq!(array.zeros(), 0);
        assert!(!array.is_empty());

        array[[1, 10, 11]] = 4.0;
        assert_eq!(array[[1, 10, 11]], 4.0);
        assert_eq!(array[[1, 10, 10]], 3.0);
        assert_eq!(array[[8, 10, 10]], 2.0);
        assert_eq!(array[[5, 10, 10]], 1.0);
        assert_eq!(array.len(), 4);
        assert_eq!(array.zeros(), 0);
        assert!(!array.is_empty());

        array[[1, 10, 9]] = 5.0;
        assert_eq!(array[[1, 10, 9]], 5.0);
        assert_eq!(array[[1, 10, 11]], 4.0);
        assert_eq!(array[[1, 10, 10]], 3.0);
        assert_eq!(array[[8, 10, 10]], 2.0);
        assert_eq!(array[[5, 10, 10]], 1.0);
        assert_eq!(array.len(), 5);
        assert_eq!(array.zeros(), 0);
        assert!(!array.is_empty());

        array[[1, 10, 0]] = 6.0;
        assert_eq!(array[[1, 10, 0]], 6.0);
        assert_eq!(array[[1, 10, 9]], 5.0);
        assert_eq!(array[[1, 10, 11]], 4.0);
        assert_eq!(array[[1, 10, 10]], 3.0);
        assert_eq!(array[[8, 10, 10]], 2.0);
        assert_eq!(array[[5, 10, 10]], 1.0);
        assert_eq!(array.len(), 6);
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
    #[should_panic]
    fn index_mut_panic_dim0() {
        let mut array = SparseArray3::new(40, 50, 50);

        array[[40, 0, 50]] = 1.0;
    }

    #[test]
    #[should_panic]
    fn index_mut_panic_dim1() {
        let mut array = SparseArray3::new(40, 50, 50);

        array[[0, 50, 0]] = 1.0;
    }

    #[test]
    #[should_panic]
    fn index_mut_panic_dim2() {
        let mut array = SparseArray3::new(40, 50, 50);

        array[[0, 0, 50]] = 1.0;
    }

    #[test]
    #[should_panic]
    fn index_panic_dim0_0() {
        let mut array = SparseArray3::new(40, 50, 50);

        array[[1, 0, 0]] = 1.0;

        assert_eq!(array[[0, 0, 0]], 0.0);
    }

    #[test]
    #[should_panic]
    fn index_panic_dim0_1() {
        let mut array = SparseArray3::new(40, 50, 50);

        array[[1, 0, 0]] = 1.0;

        assert_eq!(array[[2, 0, 0]], 0.0);
    }

    #[test]
    #[should_panic]
    fn index_panic_dim1() {
        let mut array = SparseArray3::new(40, 50, 50);

        array[[1, 0, 0]] = 1.0;

        assert_eq!(array[[1, 50, 0]], 0.0);
    }

    #[test]
    #[should_panic]
    fn index_panic_dim2_0() {
        let mut array = SparseArray3::new(40, 50, 50);

        array[[0, 0, 1]] = 1.0;

        assert_eq!(array[[0, 0, 0]], 0.0);
    }

    #[test]
    #[should_panic]
    fn index_panic_dim2_1() {
        let mut array = SparseArray3::new(40, 50, 50);

        array[[0, 0, 1]] = 1.0;

        assert_eq!(array[[0, 0, 2]], 0.0);
    }

    #[test]
    fn iterator() {
        let mut array = SparseArray3::new(40, 50, 50);

        // check empty iterator
        assert_eq!(array.iter().next(), None);

        // insert an element
        array[[2, 3, 4]] = 1.0;

        let mut iter = array.iter();

        // check iterator with one element
        assert_eq!(iter.next(), Some(((2, 3, 4), &1.0)));
        assert_eq!(iter.next(), None);

        // insert another element
        array[[2, 3, 6]] = 2.0;

        let mut iter = array.iter();

        assert_eq!(iter.next(), Some(((2, 3, 4), &1.0)));
        assert_eq!(iter.next(), Some(((2, 3, 6), &2.0)));
        assert_eq!(iter.next(), None);

        // insert yet another element
        array[[4, 5, 7]] = 3.0;

        let mut iter = array.iter();

        assert_eq!(iter.next(), Some(((2, 3, 4), &1.0)));
        assert_eq!(iter.next(), Some(((2, 3, 6), &2.0)));
        assert_eq!(iter.next(), Some(((4, 5, 7), &3.0)));
        assert_eq!(iter.next(), None);
    }
}
