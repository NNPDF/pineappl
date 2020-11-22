use std::iter;
use std::ops::{Index, IndexMut};

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
            if self.entries.is_empty() {
                self.start = index[0];
            }

            let elements = index[0] - max_index0 + 1;
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
    array: &'a SparseArray3<T>,
}

impl<'a, T> Iterator for SparseArray3Iter<'a, T> {
    type Item = ((usize, usize, usize), &'a T);

    fn next(&mut self) -> Option<Self::Item> {
        todo!();
    }
}

impl<'a, T> SparseArray3Iter<'a, T> {
    fn new(array: &'a SparseArray3<T>) -> Self {
        todo!();
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

    //pub fn iter<'a>(&'a self) -> SparseArray3Iter<'a, T> {
    //    SparseArray3Iter { array: self }
    //}
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sparse_array3() {
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
}
