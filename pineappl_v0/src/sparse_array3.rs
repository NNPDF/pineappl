//! Module containing the `SparseArray3` struct.

use serde::Deserialize;
use std::slice::Iter;

/// Struct for a sparse three-dimensional array, which is optimized for the sparsity of
/// interpolation grids.
#[derive(Deserialize)]
pub struct SparseArray3<T> {
    entries: Vec<T>,
    indices: Vec<(usize, usize)>,
    start: usize,
    dimensions: (usize, usize, usize),
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

impl<T: Copy + Default + PartialEq> Iterator for IndexedIter<'_, T> {
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
    /// Returns `true` if the array contains no element.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
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
}
