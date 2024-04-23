//! Provides the [`PackedArray`] struct.

use std::ops::{Index, IndexMut, MulAssign};
use std::{mem, vec};

use ndarray::ArrayView3;
use serde::{Deserialize, Serialize};

/// `D`-dimensional array similar to [`ndarray::ArrayBase`], except that `T::default()` is not
/// stored to save space. Instead, adjacent non-default elements are grouped together and the index
/// of their first element (`start_index`) and the length of the group (`lengths`) is stored.
#[derive(Clone, Deserialize, Serialize)]
pub struct PackedArray<T, const D: usize> {
    /// The actual values stored in the array. The length of `entries` is always the sum of the
    /// elements in `lengths`.
    entries: Vec<T>,
    /// The indices of the first elements in each group.
    start_indices: Vec<usize>,
    /// The length of each group.
    lengths: Vec<usize>,
    /// The shape (dimensions) of the array.
    shape: Vec<usize>,
}

impl<T: Copy + Default + PartialEq, const D: usize> PackedArray<T, D> {
    /// Constructs a new and empty `PackedArray` of shape `shape`.
    #[must_use]
    pub fn new(shape: [usize; D]) -> Self {
        Self {
            entries: vec![],
            start_indices: vec![],
            lengths: vec![],
            shape: shape.to_vec(),
        }
    }

    /// Returns `true` if the array contains no element.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Returns the shape of the array.
    #[must_use]
    pub fn shape(&self) -> &[usize] {
        &self.shape
    }

    /// Clears the contents of the array.
    pub fn clear(&mut self) {
        self.entries.clear();
        self.start_indices.clear();
        self.lengths.clear();
    }

    /// Returns the overhead of storing the `start_indices` and the `lengths` of the groups, in
    /// units of `f64`.
    #[must_use]
    pub fn overhead(&self) -> usize {
        ((self.start_indices.len() + self.lengths.len()) * mem::size_of::<usize>())
            / mem::size_of::<f64>()
    }

    /// Returns the number of default (zero) elements that are explicitly stored in `entries`. If
    /// there is one default element between adjacent groups, it is more economical to store the
    /// one default element explicitly and merge the two groups, than to store the `start_indices`
    /// and `lengths` of both groups.
    #[must_use]
    pub fn explicit_zeros(&self) -> usize {
        self.entries.iter().filter(|x| **x == T::default()).count()
    }

    /// Returns the number of non-default (non-zero) elements stored in the array.
    #[must_use]
    pub fn non_zeros(&self) -> usize {
        self.entries.iter().filter(|x| **x != T::default()).count()
    }

    /// Returns an `Iterator` over the non-default (non-zero) elements of this array. The type of
    /// an iterator element is `([usize; D], T)` where the first element of the tuple is the index
    /// and the second element is the value.
    pub fn indexed_iter(&self) -> impl Iterator<Item = ([usize; D], T)> + '_ {
        self.start_indices
            .iter()
            .zip(&self.lengths)
            .flat_map(|(&start_index, &length)| {
                (start_index..(start_index + length)).map(|i| unravel_index(i, &self.shape))
            })
            .zip(&self.entries)
            .filter(|&(_, entry)| *entry != Default::default())
            .map(|(indices, entry)| (indices, *entry))
    }
}

impl<T: Copy + MulAssign<T>, const D: usize> MulAssign<T> for PackedArray<T, D> {
    fn mul_assign(&mut self, rhs: T) {
        self.entries.iter_mut().for_each(|x| *x *= rhs);
    }
}

impl<T: Copy + Default + PartialEq> PackedArray<T, 3> {
    /// Converts `array` into a `PackedArray<T, 3>`.
    #[must_use]
    pub fn from_ndarray(array: ArrayView3<T>, xstart: usize, xsize: usize) -> Self {
        let shape = array.shape();

        let mut result = Self::new([xsize, shape[1], shape[2]]);

        for ((i, j, k), &entry) in array
            .indexed_iter()
            .filter(|(_, &entry)| entry != Default::default())
        {
            result[[i + xstart, j, k]] = entry;
        }

        result
    }
}

/// Converts a `multi_index` into a flat index.
fn ravel_multi_index<const D: usize>(multi_index: &[usize; D], dimensions: &[usize]) -> usize {
    assert_eq!(multi_index.len(), dimensions.len());

    multi_index
        .iter()
        .skip(1)
        .zip(dimensions.iter().skip(1))
        .fold(multi_index[0], |acc, (i, d)| acc * d + i)
}

/// Converts a flat `index` into a `multi_index`.
fn unravel_index<const D: usize>(index: usize, dimensions: &[usize]) -> [usize; D] {
    assert!(index < dimensions.iter().product());
    let mut indices = [0; D];
    indices
        .iter_mut()
        .rev()
        .zip(dimensions.iter().rev())
        .fold(index, |acc, (i, d)| {
            *i = acc % d;
            acc / d
        });
    indices
}

impl<T: Copy + Default + PartialEq, const D: usize> Index<[usize; D]> for PackedArray<T, D> {
    type Output = T;

    fn index(&self, index: [usize; D]) -> &Self::Output {
        debug_assert_eq!(index.len(), self.shape.len());
        assert!(
            index.iter().zip(self.shape.iter()).all(|(&i, &d)| i < d),
            "index {:?} is out of bounds for array of shape {:?}",
            index,
            self.shape
        );

        let raveled_index = ravel_multi_index(&index, &self.shape);
        let point = self.start_indices.partition_point(|&i| i <= raveled_index);

        assert!(
            point > 0,
            "entry at index {index:?} is implicitly set to the default value"
        );

        let start_index = self.start_indices[point - 1];
        let length = self.lengths[point - 1];

        let point_entries =
            self.lengths.iter().take(point - 1).sum::<usize>() + raveled_index - start_index;

        assert!(
            raveled_index < (start_index + length),
            "entry at index {index:?} is implicitly set to the default value"
        );

        &self.entries[point_entries]
    }
}

impl<T: Clone + Copy + Default + PartialEq, const D: usize> IndexMut<[usize; D]>
    for PackedArray<T, D>
{
    fn index_mut(&mut self, index: [usize; D]) -> &mut Self::Output {
        debug_assert_eq!(index.len(), self.shape.len());
        assert!(
            index.iter().zip(self.shape.iter()).all(|(&i, &d)| i < d),
            "index {:?} is out of bounds for array of shape {:?}",
            index,
            self.shape
        );

        let raveled_index = ravel_multi_index(&index, &self.shape);

        // the closest start_index that is greater than raveled_index
        let point = self.start_indices.partition_point(|&i| i <= raveled_index);

        // the index that is stored in `point`, translated to the indices of the entries
        let point_entries = self.lengths.iter().take(point).sum::<usize>();

        // maximum distance for merging regions
        let threshold_distance = 2;

        if point > 0 {
            let start_index = self.start_indices[point - 1];
            let length = self.lengths[point - 1];

            if raveled_index < start_index + length {
                return &mut self.entries[point_entries - length + raveled_index - start_index];
            } else if raveled_index < start_index + length + threshold_distance {
                let distance = raveled_index - (start_index + length) + 1;
                self.lengths[point - 1] += distance;

                for _ in 0..(distance) {
                    self.entries.insert(point_entries, Default::default());
                }

                if let Some(start_index_next) = self.start_indices.get(point) {
                    if raveled_index + threshold_distance >= *start_index_next {
                        // println!("test 4");
                        let distance_next = start_index_next - raveled_index;

                        self.lengths[point - 1] += distance_next - 1 + self.lengths[point];
                        self.lengths.remove(point);
                        self.start_indices.remove(point);
                        for _ in 0..(distance_next - 1) {
                            self.entries.insert(point_entries, Default::default());
                        }
                    }
                }

                return &mut self.entries[point_entries - 1 + distance];
            }
        }

        if let Some(start_index_next) = self.start_indices.get(point) {
            if raveled_index + threshold_distance >= *start_index_next {
                let distance = start_index_next - raveled_index;

                self.start_indices[point] = raveled_index;
                self.lengths[point] += distance;
                for _ in 0..distance {
                    self.entries.insert(point_entries, Default::default());
                }
                return &mut self.entries[point_entries];
            }
        }

        // in the most general case we have to insert a new start_index
        self.start_indices.insert(point, raveled_index);
        self.lengths.insert(point, 1);
        self.entries.insert(point_entries, Default::default());

        &mut self.entries[point_entries]
    }
}

#[cfg(test)]
mod tests {
    use ndarray::Array3;
    use std::mem;

    use super::*;

    #[test]
    fn unravel_index() {
        assert_eq!(super::unravel_index(0, &[3, 2]), [0, 0]);
        assert_eq!(super::unravel_index(1, &[3, 2]), [0, 1]);
        assert_eq!(super::unravel_index(2, &[3, 2]), [1, 0]);
        assert_eq!(super::unravel_index(3, &[3, 2]), [1, 1]);
        assert_eq!(super::unravel_index(4, &[3, 2]), [2, 0]);
        assert_eq!(super::unravel_index(5, &[3, 2]), [2, 1]);
    }

    #[test]
    fn ravel_multi_index() {
        assert_eq!(super::ravel_multi_index(&[0, 0], &[3, 2]), 0);
        assert_eq!(super::ravel_multi_index(&[0, 1], &[3, 2]), 1);
        assert_eq!(super::ravel_multi_index(&[1, 0], &[3, 2]), 2);
        assert_eq!(super::ravel_multi_index(&[1, 1], &[3, 2]), 3);
        assert_eq!(super::ravel_multi_index(&[2, 0], &[3, 2]), 4);
        assert_eq!(super::ravel_multi_index(&[2, 1], &[3, 2]), 5);
    }

    #[test]
    fn index() {
        let mut a = PackedArray::<f64, 2>::new([4, 2]);

        a[[0, 0]] = 1.0;
        assert_eq!(a[[0, 0]], 1.0);
        assert_eq!(a.entries, vec![1.0]);
        assert_eq!(a.start_indices, vec![0]);
        assert_eq!(a.lengths, vec![1]);

        a[[3, 0]] = 2.0;
        assert_eq!(a[[0, 0]], 1.0);
        assert_eq!(a[[3, 0]], 2.0);
        assert_eq!(a.entries, vec![1.0, 2.0]);
        assert_eq!(a.start_indices, vec![0, 6]);
        assert_eq!(a.lengths, vec![1, 1]);

        a[[3, 1]] = 3.0;
        assert_eq!(a[[0, 0]], 1.0);
        assert_eq!(a[[3, 0]], 2.0);
        assert_eq!(a[[3, 1]], 3.0);
        assert_eq!(a.entries, vec![1.0, 2.0, 3.0]);
        assert_eq!(a.start_indices, vec![0, 6]);
        assert_eq!(a.lengths, vec![1, 2]);

        a[[2, 0]] = 3.5;
        assert_eq!(a[[0, 0]], 1.0);
        assert_eq!(a[[3, 0]], 2.0);
        assert_eq!(a[[3, 1]], 3.0);
        assert_eq!(a[[2, 0]], 3.5);
        assert_eq!(a.entries, vec![1.0, 3.5, 0.0, 2.0, 3.0]);
        assert_eq!(a.start_indices, vec![0, 4]);
        assert_eq!(a.lengths, vec![1, 4]);

        a[[2, 0]] = 4.0;
        assert_eq!(a[[0, 0]], 1.0);
        assert_eq!(a[[3, 0]], 2.0);
        assert_eq!(a[[3, 1]], 3.0);
        assert_eq!(a[[2, 0]], 4.0);
        assert_eq!(a.entries, vec![1.0, 4.0, 0.0, 2.0, 3.0]);
        assert_eq!(a.start_indices, vec![0, 4]);
        assert_eq!(a.lengths, vec![1, 4]);

        a[[1, 0]] = 5.0;
        assert_eq!(a[[0, 0]], 1.0);
        assert_eq!(a[[3, 0]], 2.0);
        assert_eq!(a[[3, 1]], 3.0);
        assert_eq!(a[[2, 0]], 4.0);
        assert_eq!(a[[1, 0]], 5.0);
        assert_eq!(a.entries, vec![1.0, 0.0, 5.0, 0.0, 4.0, 0.0, 2.0, 3.0]);
        assert_eq!(a.start_indices, vec![0]);
        assert_eq!(a.lengths, vec![8]);
    }

    #[test]
    fn iter() {
        let mut a = PackedArray::<i32, 2>::new([6, 5]);
        a[[2, 2]] = 1;
        a[[2, 4]] = 2;
        a[[4, 1]] = 3;
        a[[4, 4]] = 4;
        a[[5, 0]] = 5;
        assert_eq!(
            a.indexed_iter().collect::<Vec<_>>(),
            &[
                ([2, 2], 1),
                ([2, 4], 2),
                ([4, 1], 3),
                ([4, 4], 4),
                ([5, 0], 5),
            ]
        );
    }

    #[test]
    fn index_access() {
        let mut array = PackedArray::new([40, 50, 50]);

        // after creation the array must be empty
        assert_eq!(array.overhead(), 0);
        assert!(array.is_empty());

        // insert the first element
        array[[5, 10, 10]] = 1.0;
        assert_eq!(array[[5, 10, 10]], 1.0);
        assert_eq!(array.non_zeros(), 1);
        assert_eq!(array.explicit_zeros(), 0);
        assert_eq!(array.overhead(), 2);
        assert!(!array.is_empty());

        // insert an element after the first one
        array[[8, 10, 10]] = 2.0;
        assert_eq!(array[[8, 10, 10]], 2.0);
        assert_eq!(array[[5, 10, 10]], 1.0);
        assert_eq!(array.non_zeros(), 2);
        assert_eq!(array.explicit_zeros(), 0);
        assert_eq!(array.overhead(), 4);
        assert!(!array.is_empty());

        // insert an element before the first one
        array[[1, 10, 10]] = 3.0;
        assert_eq!(array[[1, 10, 10]], 3.0);
        assert_eq!(array[[8, 10, 10]], 2.0);
        assert_eq!(array[[5, 10, 10]], 1.0);
        assert_eq!(array.non_zeros(), 3);
        assert_eq!(array.explicit_zeros(), 0);
        assert_eq!(array.overhead(), 6);
        assert!(!array.is_empty());

        array[[1, 10, 11]] = 4.0;
        assert_eq!(array[[1, 10, 11]], 4.0);
        assert_eq!(array[[1, 10, 10]], 3.0);
        assert_eq!(array[[8, 10, 10]], 2.0);
        assert_eq!(array[[5, 10, 10]], 1.0);
        assert_eq!(array.non_zeros(), 4);
        assert_eq!(array.explicit_zeros(), 0);
        assert_eq!(array.overhead(), 6);
        assert!(!array.is_empty());

        array[[1, 10, 9]] = 5.0;
        assert_eq!(array[[1, 10, 9]], 5.0);
        assert_eq!(array[[1, 10, 11]], 4.0);
        assert_eq!(array[[1, 10, 10]], 3.0);
        assert_eq!(array[[8, 10, 10]], 2.0);
        assert_eq!(array[[5, 10, 10]], 1.0);
        assert_eq!(array.non_zeros(), 5);
        assert_eq!(array.explicit_zeros(), 0);
        // dbg!(&array.start_indices);
        // dbg!(&array.lengths);
        assert_eq!(array.overhead(), 6);
        assert!(!array.is_empty());

        array[[1, 10, 0]] = 6.0;
        assert_eq!(array[[1, 10, 0]], 6.0);
        assert_eq!(array[[1, 10, 9]], 5.0);
        assert_eq!(array[[1, 10, 11]], 4.0);
        assert_eq!(array[[1, 10, 10]], 3.0);
        assert_eq!(array[[8, 10, 10]], 2.0);
        assert_eq!(array[[5, 10, 10]], 1.0);
        assert_eq!(array.non_zeros(), 6);
        assert_eq!(array.explicit_zeros(), 0);
        assert_eq!(array.overhead(), 8);
        assert!(!array.is_empty());

        array[[1, 10, 2]] = 7.0;
        assert_eq!(array[[1, 10, 2]], 7.0);
        assert_eq!(array[[1, 10, 0]], 6.0);
        assert_eq!(array[[1, 10, 9]], 5.0);
        assert_eq!(array[[1, 10, 11]], 4.0);
        assert_eq!(array[[1, 10, 10]], 3.0);
        assert_eq!(array[[8, 10, 10]], 2.0);
        assert_eq!(array[[5, 10, 10]], 1.0);
        assert_eq!(array.non_zeros(), 7);
        assert_eq!(array.overhead(), 8);
        assert!(!array.is_empty());

        // check zeros
        assert_eq!(array[[1, 10, 1]], 0.0);
        assert_eq!(array.explicit_zeros(), 1);

        array[[1, 15, 2]] = 8.0;
        assert_eq!(array[[1, 15, 2]], 8.0);
        assert_eq!(array[[1, 10, 2]], 7.0);
        assert_eq!(array[[1, 10, 0]], 6.0);
        assert_eq!(array[[1, 10, 9]], 5.0);
        assert_eq!(array[[1, 10, 11]], 4.0);
        assert_eq!(array[[1, 10, 10]], 3.0);
        assert_eq!(array[[8, 10, 10]], 2.0);
        assert_eq!(array[[5, 10, 10]], 1.0);
        assert_eq!(array.non_zeros(), 8);
        assert_eq!(array.overhead(), 10);
        assert!(!array.is_empty());

        // check zeros
        assert_eq!(array[[1, 10, 1]], 0.0);
        assert_eq!(array.explicit_zeros(), 1);

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
        assert_eq!(array.non_zeros(), 9);
        assert_eq!(array.overhead(), 10);
        assert!(!array.is_empty());

        // check zeros
        assert_eq!(array[[1, 15, 3]], 0.0);
        assert_eq!(array[[1, 10, 1]], 0.0);
        assert_eq!(array.explicit_zeros(), 2);

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
        assert_eq!(array.non_zeros(), 10);
        assert_eq!(array.overhead(), 10);
        assert!(!array.is_empty());

        // check zeros
        assert_eq!(array[[1, 15, 1]], 0.0);
        assert_eq!(array[[1, 15, 3]], 0.0);
        assert_eq!(array[[1, 10, 1]], 0.0);
        assert_eq!(array.explicit_zeros(), 3);
    }

    #[test]
    #[should_panic(expected = "index [40, 0, 50] is out of bounds for array of shape [40, 50, 50]")]
    fn index_mut_panic_dim0() {
        let mut array = PackedArray::new([40, 50, 50]);

        array[[40, 0, 50]] = 1.0;
    }

    #[test]
    #[should_panic(expected = "index [0, 50, 0] is out of bounds for array of shape [40, 50, 50]")]
    fn index_mut_panic_dim1() {
        let mut array = PackedArray::new([40, 50, 50]);

        array[[0, 50, 0]] = 1.0;
    }

    #[test]
    #[should_panic(expected = "index [0, 0, 50] is out of bounds for array of shape [40, 50, 50]")]
    fn index_mut_panic_dim2() {
        let mut array = PackedArray::new([40, 50, 50]);

        array[[0, 0, 50]] = 1.0;
    }

    #[test]
    #[should_panic(expected = "entry at index [0, 0, 0] is implicitly set to the default value")]
    fn index_panic_dim0_0() {
        let mut array = PackedArray::new([40, 50, 50]);

        array[[1, 0, 0]] = 1.0;

        assert_eq!(array[[0, 0, 0]], 0.0);
    }

    #[test]
    #[should_panic(expected = "entry at index [2, 0, 0] is implicitly set to the default value")]
    fn index_panic_dim0_1() {
        let mut array = PackedArray::new([40, 50, 50]);

        array[[1, 0, 0]] = 1.0;

        assert_eq!(array[[2, 0, 0]], 0.0);
    }

    #[test]
    #[should_panic(expected = "index [1, 50, 0] is out of bounds for array of shape [40, 50, 50]")]
    fn index_panic_dim1() {
        let mut array = PackedArray::new([40, 50, 50]);

        array[[1, 0, 0]] = 1.0;

        assert_eq!(array[[1, 50, 0]], 0.0);
    }

    #[test]
    #[should_panic(expected = "entry at index [0, 0, 0] is implicitly set to the default value")]
    fn index_panic_dim2_0() {
        let mut array = PackedArray::new([40, 50, 50]);

        array[[0, 0, 1]] = 1.0;

        assert_eq!(array[[0, 0, 0]], 0.0);
    }

    #[test]
    #[should_panic(expected = "entry at index [0, 0, 2] is implicitly set to the default value")]
    fn index_panic_dim2_1() {
        let mut array = PackedArray::new([40, 50, 50]);

        array[[0, 0, 1]] = 1.0;

        assert_eq!(array[[0, 0, 2]], 0.0);
    }

    #[test]
    fn indexed_iter() {
        let mut array = PackedArray::new([40, 50, 50]);

        // check empty iterator
        assert_eq!(array.indexed_iter().next(), None);

        // insert an element
        array[[2, 3, 4]] = 1.0;

        let mut iter = array.indexed_iter();

        // check iterator with one element
        assert_eq!(iter.next(), Some(([2, 3, 4], 1.0)));
        assert_eq!(iter.next(), None);

        mem::drop(iter);

        // insert another element
        array[[2, 3, 6]] = 2.0;

        let mut iter = array.indexed_iter();

        assert_eq!(iter.next(), Some(([2, 3, 4], 1.0)));
        assert_eq!(iter.next(), Some(([2, 3, 6], 2.0)));
        assert_eq!(iter.next(), None);

        mem::drop(iter);

        // insert yet another element
        array[[4, 5, 7]] = 3.0;

        let mut iter = array.indexed_iter();

        assert_eq!(iter.next(), Some(([2, 3, 4], 1.0)));
        assert_eq!(iter.next(), Some(([2, 3, 6], 2.0)));
        assert_eq!(iter.next(), Some(([4, 5, 7], 3.0)));
        assert_eq!(iter.next(), None);

        mem::drop(iter);

        // insert at the very first position
        array[[2, 0, 0]] = 4.0;

        let mut iter = array.indexed_iter();

        assert_eq!(iter.next(), Some(([2, 0, 0], 4.0)));
        assert_eq!(iter.next(), Some(([2, 3, 4], 1.0)));
        assert_eq!(iter.next(), Some(([2, 3, 6], 2.0)));
        assert_eq!(iter.next(), Some(([4, 5, 7], 3.0)));
        assert_eq!(iter.next(), None);
    }

    #[test]
    fn clear() {
        let mut array = PackedArray::new([40, 50, 50]);

        array[[3, 5, 1]] = 1.0;
        array[[7, 8, 9]] = 2.0;
        array[[9, 1, 4]] = 3.0;

        assert!(!array.is_empty());
        assert_eq!(array.non_zeros(), 3);
        assert_eq!(array.explicit_zeros(), 0);

        array.clear();

        assert!(array.is_empty());
        assert_eq!(array.non_zeros(), 0);
        assert_eq!(array.explicit_zeros(), 0);
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

        let array = PackedArray::from_ndarray(ndarray.view(), 3, 40);

        assert_eq!(array[[3, 4, 3]], 1.0);
        assert_eq!(array[[3, 4, 4]], 2.0);
        assert_eq!(array[[3, 4, 5]], 0.0);
        assert_eq!(array[[3, 4, 6]], 3.0);
        assert_eq!(array[[3, 5, 1]], 4.0);
        assert_eq!(array[[3, 5, 7]], 5.0);
        assert_eq!(array[[4, 3, 9]], 6.0);

        assert_eq!(array.explicit_zeros(), 1);
    }
}
