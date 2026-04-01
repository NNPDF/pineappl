//! Provides the [`PackedArray`] struct.

use ndarray::ArrayView3;
use serde::{Deserialize, Serialize};
use std::iter;
use std::mem;
use std::ops::{Index, IndexMut, MulAssign};

/// `D`-dimensional array similar to [`ndarray::ArrayBase`], except that `T::default()` is not
/// stored to save space. Instead, adjacent non-default elements are grouped together and the index
/// of their first element (`start_index`) and the length of the group (`lengths`) is stored.
#[derive(Clone, Deserialize, Serialize)]
pub struct PackedArray<T, const D: usize> {
    /// The actual values stored in the array. The length of `entries` is always the sum of the
    /// elements in `lengths`.
    entries: Vec<T>,
    /// The indices of the first elements in each group. `start_indices[i]` corresponds to the
    /// group with index `i`.
    start_indices: Vec<usize>,
    /// The length of each group. `lengths[i]` corresponds to the group with index `i`.
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
            .filter(|&(_, &entry)| entry != Default::default())
        {
            result[[i + xstart, j, k]] = entry;
        }

        result
    }
}

/// Converts a `multi_index` into a flat index.
fn ravel_multi_index<const D: usize>(multi_index: &[usize; D], shape: &[usize]) -> usize {
    assert_eq!(multi_index.len(), shape.len());

    multi_index
        .iter()
        .zip(shape)
        .fold(0, |acc, (i, d)| acc * d + i)
}

/// Converts a flat `index` into a `multi_index`.
fn unravel_index<const D: usize>(mut index: usize, shape: &[usize]) -> [usize; D] {
    assert!(index < shape.iter().product());
    let mut indices = [0; D];
    for (i, d) in indices.iter_mut().zip(shape).rev() {
        *i = index % d;
        index /= d;
    }
    indices
}

impl<T: Copy + Default + PartialEq, const D: usize> Index<[usize; D]> for PackedArray<T, D> {
    type Output = T;

    fn index(&self, index: [usize; D]) -> &Self::Output {
        assert_eq!(index.len(), self.shape.len());
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
        assert_eq!(index.len(), self.shape.len());

        // Panic if the index value for any dimension is greater or equal than the length of this
        // dimension.
        assert!(
            index.iter().zip(self.shape.iter()).all(|(&i, &d)| i < d),
            "index {:?} is out of bounds for array of shape {:?}",
            index,
            self.shape
        );

        // The insertion cases are:
        // 1. this array already stores an element at `index`:
        //    -> we just have to update this element
        // 2. this array does not store an element at `index`:
        //    a. the distance of the (raveled) `index` is `threshold_distance` away from the next
        //       or previous element that is already stored:
        //       -> we can merge the new element into already stored groups, potentially padding
        //          with `T::default()` elements
        //    b. the distance of the (raveled) `index` from the existing elements is greater than
        //       `threshold_distance`:
        //       -> we insert the element as a new group

        let raveled_index = ravel_multi_index(&index, &self.shape);

        // To determine which groups the new element is close to, `point` is the index of the
        // start_index of the first group after the new element. `point` is 0 if no elements before
        // the new element are stored, and point is `self.start_indices.len()` if no elements after
        // the new element are stored.
        let point = self.start_indices.partition_point(|&i| i <= raveled_index);

        // `point_entries` is the index of the first element of the next group, given in
        // `self.entries`, i.e. the element at index `self.start_indices[point]`.
        let point_entries = self.lengths.iter().take(point).sum::<usize>();

        // Maximum distance for merging groups. If the new element is within `threshold_distance`
        // of an existing group (i.e. there are `threshold_distance - 1` implicit elements
        // between them), we merge the new element into the existing group. We choose 2 as the
        // `threshold_distance` based on memory: in the case of `T` = `f64`, it is more economical
        // to store one zero explicitly than to store the start_index and length of a new group.
        let threshold_distance = 2;

        // If `point > 0`, there is at least one group preceding the new element. Thus, in the
        // following we determine if we can insert the new element into this group.
        if point > 0 {
            // start_index and length of the group before the new element, i.e. the group
            // (potentially) getting the new element
            let start_index = self.start_indices[point - 1];
            let length = self.lengths[point - 1];

            // Case 1: an element is already stored at this `index`
            if raveled_index < start_index + length {
                return &mut self.entries[point_entries - length + raveled_index - start_index];
            // Case 2a: the new element can be merged into the preceding group
            } else if raveled_index < start_index + length + threshold_distance {
                let distance = raveled_index - (start_index + length) + 1;
                // Merging happens by increasing the length of the group
                self.lengths[point - 1] += distance;
                // and inserting the necessary number of default elements.
                self.entries.splice(
                    point_entries..point_entries,
                    iter::repeat(Default::default()).take(distance),
                );

                // If the new element is within `threshold_distance` of the *next* group, we merge
                // the next group into this group.
                if let Some(start_index_next) = self.start_indices.get(point) {
                    if raveled_index + threshold_distance >= *start_index_next {
                        let distance_next = start_index_next - raveled_index;

                        // Increase the length of this group
                        self.lengths[point - 1] += distance_next - 1 + self.lengths[point];
                        // and remove the next group. we don't have to manipulate `self.entries`,
                        // since the grouping of the elements is handled only by
                        // `self.start_indices` and `self.lengths`
                        self.lengths.remove(point);
                        self.start_indices.remove(point);
                        // Insert the default elements between the groups.
                        self.entries.splice(
                            point_entries..point_entries,
                            iter::repeat(Default::default()).take(distance_next - 1),
                        );
                    }
                }

                return &mut self.entries[point_entries - 1 + distance];
            }
        }

        // Case 2a: the new element can be merged into the next group. No `self.lengths.remove` and
        // `self.start_indices.remove` here, since we are not merging two groups.
        if let Some(start_index_next) = self.start_indices.get(point) {
            if raveled_index + threshold_distance >= *start_index_next {
                let distance = start_index_next - raveled_index;

                self.start_indices[point] = raveled_index;
                self.lengths[point] += distance;
                self.entries.splice(
                    point_entries..point_entries,
                    iter::repeat(Default::default()).take(distance),
                );
                return &mut self.entries[point_entries];
            }
        }

        // Case 2b: we insert a new group of length 1
        self.start_indices.insert(point, raveled_index);
        self.lengths.insert(point, 1);
        self.entries.insert(point_entries, Default::default());

        &mut self.entries[point_entries]
    }
}
