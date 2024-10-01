//! TODO

use super::interpolation::Interp;
use super::packed_array::PackedArray;
use super::subgrid::{Mu2, NodeValues, Stats, Subgrid, SubgridEnum, SubgridIndexedIter};
use itertools::izip;
use serde::{Deserialize, Serialize};
use std::mem;

/// TODO
#[derive(Clone, Deserialize, Serialize)]
pub struct PackedQ1X2SubgridV1 {
    array: PackedArray<f64>,
    node_values: Vec<NodeValues>,
}

impl PackedQ1X2SubgridV1 {
    /// Constructor.
    #[must_use]
    pub const fn new(array: PackedArray<f64>, node_values: Vec<NodeValues>) -> Self {
        Self { array, node_values }
    }
}

impl Subgrid for PackedQ1X2SubgridV1 {
    fn fill(&mut self, _: &[Interp], _: &[f64], _: f64) {
        panic!("PackedQ1X2SubgridV1 doesn't support the fill operation");
    }

    fn node_values(&self) -> Vec<NodeValues> {
        self.node_values.clone()
    }

    fn is_empty(&self) -> bool {
        self.array.is_empty()
    }

    fn merge(&mut self, other: &SubgridEnum, transpose: Option<(usize, usize)>) {
        let lhs_node_values = self.node_values();
        let mut rhs_node_values = other.node_values();
        let mut new_node_values = lhs_node_values.clone();
        if let Some((a, b)) = transpose {
            rhs_node_values.swap(a, b);
        }

        if self.node_values() != rhs_node_values {
            for (lhs, rhs) in new_node_values.iter_mut().zip(&rhs_node_values) {
                lhs.extend(rhs);
            }

            let mut array = PackedArray::new(new_node_values.iter().map(NodeValues::len).collect());

            for (indices, value) in self.array.indexed_iter3() {
                let target: Vec<_> = izip!(indices, &new_node_values, &lhs_node_values)
                    .map(|(index, new, lhs)| new.find(lhs.get(index)).unwrap())
                    .collect();

                array[target.as_slice()] = value;
            }

            self.array = array;
            self.node_values = new_node_values.clone();
        }

        for (mut indices, value) in other.indexed_iter() {
            if let Some((a, b)) = transpose {
                indices.swap(a, b);
            }

            let target: Vec<_> = izip!(indices, &new_node_values, &rhs_node_values)
                .map(|(index, new, rhs)| new.find(rhs.get(index)).unwrap())
                .collect();

            self.array[target.as_slice()] += value;
        }
    }

    fn scale(&mut self, factor: f64) {
        self.array *= factor;
    }

    fn symmetrize(&mut self, a: usize, b: usize) {
        let mut new_array = PackedArray::new(self.array.shape().to_vec());

        for (mut index, sigma) in self.array.indexed_iter3() {
            // TODO: why not the other way around?
            if index[b] < index[a] {
                index.swap(a, b);
            }

            new_array[index.as_slice()] += sigma;
        }

        self.array = new_array;
    }

    fn indexed_iter(&self) -> SubgridIndexedIter {
        Box::new(self.array.indexed_iter3())
    }

    fn stats(&self) -> Stats {
        Stats {
            total: self.array.shape().iter().product(),
            allocated: self.array.non_zeros() + self.array.explicit_zeros(),
            zeros: self.array.explicit_zeros(),
            overhead: self.array.overhead(),
            bytes_per_value: mem::size_of::<f64>(),
        }
    }

    fn static_scale(&self) -> Option<Mu2> {
        if let &[static_scale] = self.node_values()[0].values().as_slice() {
            Some(Mu2 {
                ren: static_scale,
                fac: static_scale,
                frg: -1.0,
            })
        } else {
            None
        }
    }
}

impl From<&SubgridEnum> for PackedQ1X2SubgridV1 {
    fn from(subgrid: &SubgridEnum) -> Self {
        // find smallest ranges
        let ranges: Vec<_> = subgrid.indexed_iter().fold(
            subgrid
                .node_values()
                .iter()
                .map(|values| values.values().len()..0)
                .collect(),
            |mut prev, (indices, _)| {
                for (i, index) in indices.iter().enumerate() {
                    prev[i].start = prev[i].start.min(*index);
                    prev[i].end = prev[i].end.max(*index + 1);
                }
                prev
            },
        );

        let mut new_node_values: Vec<_> = subgrid
            .node_values()
            .iter()
            .zip(&ranges)
            .map(|(values, range)| NodeValues::UseThese(values.values()[range.clone()].to_vec()))
            .collect();
        let static_scale = if let Some(Mu2 { ren, fac, frg }) = subgrid.static_scale() {
            assert_eq!(ren, fac);
            assert_eq!(frg, -1.0);
            new_node_values[0] = NodeValues::UseThese(vec![fac]);
            true
        } else {
            false
        };

        let mut array = PackedArray::new(
            new_node_values
                .iter()
                .map(|values| values.values().len())
                .collect(),
        );

        for (mut indices, value) in subgrid.indexed_iter() {
            for (idx, (index, range)) in indices.iter_mut().zip(&ranges).enumerate() {
                // TODO: generalize static scale detection
                *index = if static_scale && idx == 0 {
                    // if there's a static scale we want every value to be added to same grid point
                    0
                } else {
                    *index - range.start
                };
            }

            array[indices.as_slice()] += value;
        }

        Self::new(array, new_node_values)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::v0;

    #[test]
    #[should_panic(expected = "PackedQ1X2SubgridV1 doesn't support the fill operation")]
    fn fill_packed_q1x2_subgrid_v1() {
        let mut subgrid = PackedQ1X2SubgridV1::new(
            PackedArray::new(vec![0, 0, 0]),
            vec![NodeValues::UseThese(Vec::new()); 3],
        );
        subgrid.fill(&v0::default_interps(), &[0.0; 3], 0.0);
    }

    #[test]
    fn test_v1() {
        let x = vec![
            0.015625, 0.03125, 0.0625, 0.125, 0.1875, 0.25, 0.375, 0.5, 0.75, 1.0,
        ];
        let mut grid1: SubgridEnum = PackedQ1X2SubgridV1::new(
            PackedArray::new(vec![1, 10, 10]),
            vec![
                NodeValues::UseThese(vec![0.0]),
                NodeValues::UseThese(x.clone()),
                NodeValues::UseThese(x.clone()),
            ],
        )
        .into();

        assert_eq!(
            grid1.node_values(),
            vec![
                NodeValues::UseThese(vec![0.0]),
                NodeValues::UseThese(x.clone()),
                NodeValues::UseThese(x.clone())
            ]
        );

        assert!(grid1.is_empty());

        // only use exactly representable numbers here so that we can avoid using approx_eq
        if let SubgridEnum::PackedQ1X2SubgridV1(ref mut x) = grid1 {
            x.array[[0, 1, 2]] = 1.0;
            x.array[[0, 1, 3]] = 2.0;
            x.array[[0, 4, 3]] = 4.0;
            x.array[[0, 7, 1]] = 8.0;
        }

        assert!(!grid1.is_empty());

        assert_eq!(grid1.indexed_iter().next(), Some((vec![0, 1, 2], 1.0)));
        assert_eq!(grid1.indexed_iter().nth(1), Some((vec![0, 1, 3], 2.0)));
        assert_eq!(grid1.indexed_iter().nth(2), Some((vec![0, 4, 3], 4.0)));
        assert_eq!(grid1.indexed_iter().nth(3), Some((vec![0, 7, 1], 8.0)));

        // create grid with transposed entries, but different q2
        let mut grid2: SubgridEnum = PackedQ1X2SubgridV1::new(
            PackedArray::new(vec![1, 10, 10]),
            vec![
                NodeValues::UseThese(vec![1.0]),
                NodeValues::UseThese(x.clone()),
                NodeValues::UseThese(x.clone()),
            ],
        )
        .into();
        if let SubgridEnum::PackedQ1X2SubgridV1(ref mut x) = grid2 {
            x.array[[0, 2, 1]] = 1.0;
            x.array[[0, 3, 1]] = 2.0;
            x.array[[0, 3, 4]] = 4.0;
            x.array[[0, 1, 7]] = 8.0;
        }

        assert_eq!(grid2.indexed_iter().next(), Some((vec![0, 1, 7], 8.0)));
        assert_eq!(grid2.indexed_iter().nth(1), Some((vec![0, 2, 1], 1.0)));
        assert_eq!(grid2.indexed_iter().nth(2), Some((vec![0, 3, 1], 2.0)));
        assert_eq!(grid2.indexed_iter().nth(3), Some((vec![0, 3, 4], 4.0)));

        grid1.merge(&mut grid2, None);

        // the luminosity function is symmetric, so after symmetrization the result must be
        // unchanged
        grid1.symmetrize(1, 2);

        grid1.scale(2.0);

        assert_eq!(
            grid1.stats(),
            Stats {
                total: 200,
                allocated: 8,
                zeros: 0,
                overhead: 12,
                bytes_per_value: 8,
            }
        );
    }
}
