use std::iter;
use std::ops::{Index, IndexMut};

#[derive(Clone, Default)]
struct Offset {
    space: usize,
    offset: usize,
}

#[derive(Clone)]
struct SparseArray2<T> {
    entries: Vec<T>,
    indices: Vec<Offset>,
}

impl<T> SparseArray2<T> {
    fn new(ny: usize) -> Self {
        Self {
            entries: vec![],
            indices: vec![Offset::default(); ny + 1],
        }
    }
}

// TODO: merge SparseArray2 into SparseArray3

/// Struct for a sparse three-dimensional array, which is optimized for the sparsity of
/// interpolation grids.
pub struct SparseArray3<T> {
    entries: Vec<SparseArray2<T>>,
    start: usize,
    dimensions: (usize, usize, usize),
}

impl<T> Index<[usize; 3]> for SparseArray3<T> {
    type Output = T;

    fn index(&self, index: [usize; 3]) -> &Self::Output {
        let array2 = &self.entries[index[0] - self.start];
        let indices_a = &array2.indices[index[1]];
        let indices_b = &array2.indices[index[1] + 1];

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

        &array2.entries[(offset + index[2]) - zeros_left]
    }
}

impl<T: Clone + Default> IndexMut<[usize; 3]> for SparseArray3<T> {
    fn index_mut(&mut self, index: [usize; 3]) -> &mut Self::Output {
        if index[0] >= self.dimensions.0 {
            panic!();
        }

        if (index[0] < self.start) || (index[0] >= (self.start + self.entries.len())) {
            let mut count = 1;
            let mut range = 0..0;

            if self.entries.is_empty() {
                self.start = index[0];
            } else {
                let new_start = self.start.min(index[0]);

                if new_start == self.start {
                    range = self.entries.len()..self.entries.len();
                    count = (index[0] + 1) - (self.start + self.entries.len());
                } else {
                    count = self.start - new_start;
                }

                self.start = new_start;
            }

            self.entries.splice(
                range,
                iter::repeat(SparseArray2::new(self.dimensions.1)).take(count),
            );
        }

        let array2 = &mut self.entries[index[0] - self.start];
        let mut zeros_left = array2.indices[index[1]].space;
        let mut offset = array2.indices[index[1]].offset;
        let non_zeros = array2.indices[index[1] + 1].offset - offset;

        let elements;
        let mut splice = offset;

        if non_zeros == 0 {
            elements = 1;
            array2.indices[index[1]].space = index[2];
        } else if index[2] < zeros_left {
            elements = zeros_left - index[2];
            array2.indices[index[1]].space -= elements;
        } else if index[2] >= (non_zeros + zeros_left) {
            splice += non_zeros;
            elements = index[2] - zeros_left;
        } else {
            return &mut array2.entries[(offset + index[2]) - zeros_left];
        }

        array2
            .entries
            .splice(splice..splice, iter::repeat(T::default()).take(elements));
        array2
            .indices
            .iter_mut()
            .skip(index[1] + 1)
            .for_each(|i| i.offset += elements);

        zeros_left = array2.indices[index[1]].space;
        offset = array2.indices[index[1]].offset;

        &mut array2.entries[(offset + index[2]) - zeros_left]
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
    /// Returns a new array with the specified dimensions for the three dimensions.
    #[must_use]
    pub fn new(nx: usize, ny: usize, nz: usize) -> Self {
        Self {
            entries: vec![],
            start: 0,
            dimensions: (nx, ny, nz),
        }
    }

    /// Returns the number of default (zero) elements in this array.
    #[must_use]
    pub fn zeros(&self) -> usize {
        let mut zeros = 0;

        for array2 in &self.entries {
            zeros += array2
                .entries
                .iter()
                .filter(|x| **x == T::default())
                .count();
        }

        zeros
    }

    /// Returns the number of non-default (non-zero) elements in this array.
    #[must_use]
    pub fn len(&self) -> usize {
        let mut len = 0;

        for array2 in &self.entries {
            len += array2
                .entries
                .iter()
                .filter(|x| **x != T::default())
                .count();
        }

        len
    }

    /// Returns `true` if the array contains no element.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    pub fn iter<'a>(&'a self) -> SparseArray3Iter<'a, T> {
        SparseArray3Iter { array: self }
    }
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
