//! TODO

const EVOL_BASIS_IDS: [i32; 12] = [100, 103, 108, 115, 124, 135, 200, 203, 208, 215, 224, 235];

/// Translates IDs from the evolution basis into IDs using PDG Monte Carlo IDs.
#[must_use]
pub fn evol_to_pdg_mc_ids(id: i32) -> Vec<(i32, f64)> {
    match id {
        100 => vec![
            (2, 1.0),
            (-2, 1.0),
            (1, 1.0),
            (-1, 1.0),
            (3, 1.0),
            (-3, 1.0),
            (4, 1.0),
            (-4, 1.0),
            (5, 1.0),
            (-5, 1.0),
            (6, 1.0),
            (-6, 1.0),
        ],
        103 => vec![(2, 1.0), (-2, 1.0), (1, -1.0), (-1, -1.0)],
        108 => vec![
            (2, 1.0),
            (-2, 1.0),
            (1, 1.0),
            (-1, 1.0),
            (3, -2.0),
            (-3, -2.0),
        ],
        115 => vec![
            (2, 1.0),
            (-2, 1.0),
            (1, 1.0),
            (-1, 1.0),
            (3, 1.0),
            (-3, 1.0),
            (4, -3.0),
            (-4, -3.0),
        ],
        124 => vec![
            (2, 1.0),
            (-2, 1.0),
            (1, 1.0),
            (-1, 1.0),
            (3, 1.0),
            (-3, 1.0),
            (4, 1.0),
            (-4, 1.0),
            (5, -4.0),
            (-5, -4.0),
        ],
        135 => vec![
            (2, 1.0),
            (-2, 1.0),
            (1, 1.0),
            (-1, 1.0),
            (3, 1.0),
            (-3, 1.0),
            (4, 1.0),
            (-4, 1.0),
            (5, 1.0),
            (-5, 1.0),
            (6, -5.0),
            (-6, -5.0),
        ],
        200 => vec![
            (1, 1.0),
            (-1, -1.0),
            (2, 1.0),
            (-2, -1.0),
            (3, 1.0),
            (-3, -1.0),
            (4, 1.0),
            (-4, -1.0),
            (5, 1.0),
            (-5, -1.0),
            (6, 1.0),
            (-6, -1.0),
        ],
        203 => vec![(2, 1.0), (-2, -1.0), (1, -1.0), (-1, 1.0)],
        208 => vec![
            (2, 1.0),
            (-2, -1.0),
            (1, 1.0),
            (-1, -1.0),
            (3, -2.0),
            (-3, 2.0),
        ],
        215 => vec![
            (2, 1.0),
            (-2, -1.0),
            (1, 1.0),
            (-1, -1.0),
            (3, 1.0),
            (-3, -1.0),
            (4, -3.0),
            (-4, 3.0),
        ],
        224 => vec![
            (2, 1.0),
            (-2, -1.0),
            (1, 1.0),
            (-1, -1.0),
            (3, 1.0),
            (-3, -1.0),
            (4, 1.0),
            (-4, -1.0),
            (5, -4.0),
            (-5, 4.0),
        ],
        235 => vec![
            (2, 1.0),
            (-2, -1.0),
            (1, 1.0),
            (-1, -1.0),
            (3, 1.0),
            (-3, -1.0),
            (4, 1.0),
            (-4, -1.0),
            (5, 1.0),
            (-5, -1.0),
            (6, -5.0),
            (-6, 5.0),
        ],
        _ => vec![(id, 1.0)],
    }
}

/// Return the charge-conjugated PDG ID of `pid`.
#[must_use]
pub const fn charge_conjugate_pdg_pid(pid: i32) -> i32 {
    match pid {
        21 | 22 => pid,
        _ => -pid,
    }
}

/// Return the charge-conjugated particle ID of `pid` for the basis `lumi_id_types`. The returned
/// tuple contains a factor that possible arises during the carge conjugation.
///
/// # Panics
///
/// TODO
#[must_use]
pub fn charge_conjugate(lumi_id_types: &str, pid: i32) -> (i32, f64) {
    match (lumi_id_types, pid) {
        ("evol", 100 | 103 | 108 | 115 | 124 | 135) => (pid, 1.0),
        ("evol", 200 | 203 | 208 | 215 | 224 | 235) => (pid, -1.0),
        ("evol" | "pdg_mc_ids", _) => (charge_conjugate_pdg_pid(pid), 1.0),
        _ => todo!(),
    }
}

/// Given the particle IDs in `pids`, determine the right string for `lumi_id_types` stored in
/// `Grid`.
#[must_use]
pub fn determine_lumi_id_types(pids: &[i32]) -> String {
    // if we find more than 3 pids that are recognized to be from the evolution basis, declare
    // it to be the evolution basis (that's a heuristic), otherwise PDG MC IDs
    if pids
        .iter()
        .filter(|&pid| EVOL_BASIS_IDS.iter().any(|evol_pid| pid == evol_pid))
        .count()
        > 3
    {
        "evol".to_string()
    } else {
        "pdg_mc_ids".to_string()
    }
}

/// Given `tuples` represting a linear combination of PDG MC IDs, return a PID for the `evol`
/// basis. The order of each tuple in `tuples` is not relevant. This function inverts
/// [`evol_to_pdg_mc_ids`]. If the inversion is not possible, `None` is returned.
#[must_use]
pub fn pdg_mc_ids_to_evol(tuples: &[(i32, f64)]) -> Option<i32> {
    let mut tuples = tuples.to_vec();
    tuples.retain(|&(_, factor)| factor != 0.0);
    tuples.sort_by_key(|&(id, _)| id);
    let tuples = tuples;

    for &evol_pid in &EVOL_BASIS_IDS {
        let mut evol_vec = evol_to_pdg_mc_ids(evol_pid);
        evol_vec.sort_by_key(|&(id, _)| id);
        let evol_vec = evol_vec;

        if evol_vec == tuples {
            return Some(evol_pid);
        }
    }

    let non_zero: Vec<_> = tuples
        .into_iter()
        .filter(|&(_, factor)| factor != 0.0)
        .collect();

    if let &[(pid, factor)] = non_zero.as_slice() {
        if factor == 1.0 {
            return Some(pid);
        }
    }

    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test() {
        // check photon
        assert_eq!(evol_to_pdg_mc_ids(21), [(21, 1.0)]);

        // check gluon
        assert_eq!(evol_to_pdg_mc_ids(22), [(22, 1.0)]);

        // check singlet
        assert_eq!(
            evol_to_pdg_mc_ids(100),
            [
                (2, 1.0),
                (-2, 1.0),
                (1, 1.0),
                (-1, 1.0),
                (3, 1.0),
                (-3, 1.0),
                (4, 1.0),
                (-4, 1.0),
                (5, 1.0),
                (-5, 1.0),
                (6, 1.0),
                (-6, 1.0),
            ]
        );

        // check T3
        assert_eq!(
            evol_to_pdg_mc_ids(103),
            [(2, 1.0), (-2, 1.0), (1, -1.0), (-1, -1.0)]
        );

        // check T8
        assert_eq!(
            evol_to_pdg_mc_ids(108),
            [
                (2, 1.0),
                (-2, 1.0),
                (1, 1.0),
                (-1, 1.0),
                (3, -2.0),
                (-3, -2.0),
            ],
        );

        // check T15
        assert_eq!(
            evol_to_pdg_mc_ids(115),
            [
                (2, 1.0),
                (-2, 1.0),
                (1, 1.0),
                (-1, 1.0),
                (3, 1.0),
                (-3, 1.0),
                (4, -3.0),
                (-4, -3.0),
            ],
        );

        // check T24
        assert_eq!(
            evol_to_pdg_mc_ids(124),
            [
                (2, 1.0),
                (-2, 1.0),
                (1, 1.0),
                (-1, 1.0),
                (3, 1.0),
                (-3, 1.0),
                (4, 1.0),
                (-4, 1.0),
                (5, -4.0),
                (-5, -4.0),
            ],
        );

        // check T35
        assert_eq!(
            evol_to_pdg_mc_ids(135),
            [
                (2, 1.0),
                (-2, 1.0),
                (1, 1.0),
                (-1, 1.0),
                (3, 1.0),
                (-3, 1.0),
                (4, 1.0),
                (-4, 1.0),
                (5, 1.0),
                (-5, 1.0),
                (6, -5.0),
                (-6, -5.0),
            ],
        );

        // check valence
        assert_eq!(
            evol_to_pdg_mc_ids(200),
            [
                (1, 1.0),
                (-1, -1.0),
                (2, 1.0),
                (-2, -1.0),
                (3, 1.0),
                (-3, -1.0),
                (4, 1.0),
                (-4, -1.0),
                (5, 1.0),
                (-5, -1.0),
                (6, 1.0),
                (-6, -1.0),
            ],
        );

        // check V3
        assert_eq!(
            evol_to_pdg_mc_ids(203),
            [(2, 1.0), (-2, -1.0), (1, -1.0), (-1, 1.0)],
        );

        // check V8
        assert_eq!(
            evol_to_pdg_mc_ids(208),
            [
                (2, 1.0),
                (-2, -1.0),
                (1, 1.0),
                (-1, -1.0),
                (3, -2.0),
                (-3, 2.0),
            ],
        );

        // check V15
        assert_eq!(
            evol_to_pdg_mc_ids(215),
            [
                (2, 1.0),
                (-2, -1.0),
                (1, 1.0),
                (-1, -1.0),
                (3, 1.0),
                (-3, -1.0),
                (4, -3.0),
                (-4, 3.0),
            ],
        );

        // check V24
        assert_eq!(
            evol_to_pdg_mc_ids(224),
            [
                (2, 1.0),
                (-2, -1.0),
                (1, 1.0),
                (-1, -1.0),
                (3, 1.0),
                (-3, -1.0),
                (4, 1.0),
                (-4, -1.0),
                (5, -4.0),
                (-5, 4.0),
            ],
        );

        // check V35
        assert_eq!(
            evol_to_pdg_mc_ids(235),
            [
                (2, 1.0),
                (-2, -1.0),
                (1, 1.0),
                (-1, -1.0),
                (3, 1.0),
                (-3, -1.0),
                (4, 1.0),
                (-4, -1.0),
                (5, 1.0),
                (-5, -1.0),
                (6, -5.0),
                (-6, 5.0),
            ],
        );
    }

    #[test]
    fn test_pdg_mc_ids_to_evol() {
        assert_eq!(pdg_mc_ids_to_evol(&[]), None);

        // check photon
        assert_eq!(
            pdg_mc_ids_to_evol(&[
                (22, 1.0),
                (-6, 0.0),
                (-5, 0.0),
                (-4, 0.0),
                (-3, 0.0),
                (-2, 0.0),
                (-1, 0.0),
                (21, 0.0),
                (1, 0.0),
                (2, 0.0),
                (3, 0.0),
                (4, 0.0),
                (5, 0.0),
                (6, 0.0),
            ]),
            Some(22)
        );

        // check gluon
        assert_eq!(
            pdg_mc_ids_to_evol(&[
                (22, 0.0),
                (-6, 0.0),
                (-5, 0.0),
                (-4, 0.0),
                (-3, 0.0),
                (-2, 0.0),
                (-1, 0.0),
                (21, 1.0),
                (1, 0.0),
                (2, 0.0),
                (3, 0.0),
                (4, 0.0),
                (5, 0.0),
                (6, 0.0),
            ]),
            Some(21)
        );

        // check singlet
        assert_eq!(
            pdg_mc_ids_to_evol(&[
                (22, 0.0),
                (-6, 1.0),
                (-5, 1.0),
                (-4, 1.0),
                (-3, 1.0),
                (-2, 1.0),
                (-1, 1.0),
                (21, 0.0),
                (1, 1.0),
                (2, 1.0),
                (3, 1.0),
                (4, 1.0),
                (5, 1.0),
                (6, 1.0),
            ]),
            Some(100)
        );

        // check T3
        assert_eq!(
            pdg_mc_ids_to_evol(&[
                (22, 0.0),
                (-6, 0.0),
                (-5, 0.0),
                (-4, 0.0),
                (-3, 0.0),
                (-2, 1.0),
                (-1, -1.0),
                (21, 0.0),
                (1, -1.0),
                (2, 1.0),
                (3, 0.0),
                (4, 0.0),
                (5, 0.0),
                (6, 0.0),
            ]),
            Some(103)
        );

        // check T8
        assert_eq!(
            pdg_mc_ids_to_evol(&[
                (22, 0.0),
                (-6, 0.0),
                (-5, 0.0),
                (-4, 0.0),
                (-3, -2.0),
                (-2, 1.0),
                (-1, 1.0),
                (21, 0.0),
                (1, 1.0),
                (2, 1.0),
                (3, -2.0),
                (4, 0.0),
                (5, 0.0),
                (6, 0.0),
            ]),
            Some(108)
        );

        // check T15
        assert_eq!(
            pdg_mc_ids_to_evol(&[
                (22, 0.0),
                (-6, 0.0),
                (-5, 0.0),
                (-4, -3.0),
                (-3, 1.0),
                (-2, 1.0),
                (-1, 1.0),
                (21, 0.0),
                (1, 1.0),
                (2, 1.0),
                (3, 1.0),
                (4, -3.0),
                (5, 0.0),
                (6, 0.0),
            ]),
            Some(115)
        );

        // check T24
        assert_eq!(
            pdg_mc_ids_to_evol(&[
                (22, 0.0),
                (-6, 0.0),
                (-5, -4.0),
                (-4, 1.0),
                (-3, 1.0),
                (-2, 1.0),
                (-1, 1.0),
                (21, 0.0),
                (1, 1.0),
                (2, 1.0),
                (3, 1.0),
                (4, 1.0),
                (5, -4.0),
                (6, 0.0),
            ]),
            Some(124)
        );

        // check T35
        assert_eq!(
            pdg_mc_ids_to_evol(&[
                (22, 0.0),
                (-6, -5.0),
                (-5, 1.0),
                (-4, 1.0),
                (-3, 1.0),
                (-2, 1.0),
                (-1, 1.0),
                (21, 0.0),
                (1, 1.0),
                (2, 1.0),
                (3, 1.0),
                (4, 1.0),
                (5, 1.0),
                (6, -5.0),
            ]),
            Some(135)
        );

        // check valence
        assert_eq!(
            pdg_mc_ids_to_evol(&[
                (22, 0.0),
                (-6, -1.0),
                (-5, -1.0),
                (-4, -1.0),
                (-3, -1.0),
                (-2, -1.0),
                (-1, -1.0),
                (21, 0.0),
                (1, 1.0),
                (2, 1.0),
                (3, 1.0),
                (4, 1.0),
                (5, 1.0),
                (6, 1.0),
            ]),
            Some(200)
        );

        // check V3
        assert_eq!(
            pdg_mc_ids_to_evol(&[
                (22, 0.0),
                (-6, 0.0),
                (-5, 0.0),
                (-4, 0.0),
                (-3, 0.0),
                (-2, -1.0),
                (-1, 1.0),
                (21, 0.0),
                (1, -1.0),
                (2, 1.0),
                (3, 0.0),
                (4, 0.0),
                (5, 0.0),
                (6, 0.0),
            ]),
            Some(203)
        );

        // check V8
        assert_eq!(
            pdg_mc_ids_to_evol(&[
                (22, 0.0),
                (-6, 0.0),
                (-5, 0.0),
                (-4, 0.0),
                (-3, 2.0),
                (-2, -1.0),
                (-1, -1.0),
                (21, 0.0),
                (1, 1.0),
                (2, 1.0),
                (3, -2.0),
                (4, 0.0),
                (5, 0.0),
                (6, 0.0),
            ]),
            Some(208)
        );

        // check V15
        assert_eq!(
            pdg_mc_ids_to_evol(&[
                (22, 0.0),
                (-6, 0.0),
                (-5, 0.0),
                (-4, 3.0),
                (-3, -1.0),
                (-2, -1.0),
                (-1, -1.0),
                (21, 0.0),
                (1, 1.0),
                (2, 1.0),
                (3, 1.0),
                (4, -3.0),
                (5, 0.0),
                (6, 0.0),
            ]),
            Some(215)
        );

        // check V24
        assert_eq!(
            pdg_mc_ids_to_evol(&[
                (22, 0.0),
                (-6, 0.0),
                (-5, 4.0),
                (-4, -1.0),
                (-3, -1.0),
                (-2, -1.0),
                (-1, -1.0),
                (21, 0.0),
                (1, 1.0),
                (2, 1.0),
                (3, 1.0),
                (4, 1.0),
                (5, -4.0),
                (6, 0.0),
            ]),
            Some(224)
        );

        // check V35
        assert_eq!(
            pdg_mc_ids_to_evol(&[
                (22, 0.0),
                (-6, 5.0),
                (-5, -1.0),
                (-4, -1.0),
                (-3, -1.0),
                (-2, -1.0),
                (-1, -1.0),
                (21, 0.0),
                (1, 1.0),
                (2, 1.0),
                (3, 1.0),
                (4, 1.0),
                (5, 1.0),
                (6, -5.0),
            ]),
            Some(235)
        );
    }

    #[test]
    fn test_determine_lumi_id_types() {
        assert_eq!(
            determine_lumi_id_types(&[22, -6, -5, -4, -3, -2, -1, 21, 1, 2, 3, 4, 5, 6]),
            "pdg_mc_ids"
        );

        assert_eq!(
            determine_lumi_id_types(&[
                22, 100, 200, 21, 100, 103, 108, 115, 124, 135, 203, 208, 215, 224, 235
            ]),
            "evol"
        );
    }
}
