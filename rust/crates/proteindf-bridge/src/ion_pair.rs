// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::amino_acid::AminoAcid;
use crate::atom_group::AtomGroup;
use crate::position::Position;

type ResidueSites = Vec<(String, Vec<(Position, String)>)>;

/// Represents a detected ion-pair interaction between an anion and a cation site.
#[derive(Debug, Clone, PartialEq)]
pub struct IonPairRecord {
    pub anion_path: String,
    pub cation_path: String,
    pub anion_type: String,
    pub cation_type: String,
}

impl IonPairRecord {
    pub fn new(
        anion_path: impl Into<String>,
        cation_path: impl Into<String>,
        anion_type: impl Into<String>,
        cation_type: impl Into<String>,
    ) -> Self {
        Self {
            anion_path: anion_path.into(),
            cation_path: cation_path.into(),
            anion_type: anion_type.into(),
            cation_type: cation_type.into(),
        }
    }

    /// Converts this record into a tuple `(anion_path, cation_path, anion_type, cation_type)` matching Python.
    pub fn into_tuple(self) -> (String, String, String, String) {
        (
            self.anion_path,
            self.cation_path,
            self.anion_type,
            self.cation_type,
        )
    }
}

/// Detects salt bridges and ion-pair interactions in protein models corresponding to `proteindf_bridge.ionpair.IonPair`.
#[derive(Debug, Clone)]
pub struct IonPair {
    model: AtomGroup,
}

impl IonPair {
    /// Distance threshold for ion-pair interaction (4.0 Å).
    pub const MAX_DISTANCE: f64 = 4.0;

    /// Creates a new `IonPair` detector for the given model.
    pub fn new(model: &AtomGroup) -> Self {
        Self {
            model: model.clone(),
        }
    }

    /// Finds and returns all ion-pair interactions within `MAX_DISTANCE` (4.0 Å).
    pub fn get_ion_pairs(&self) -> Vec<IonPairRecord> {
        let mut ion_pairs = Vec::new();
        let (anion_list, cation_list) = self.get_ion_list();

        for (anion_path, anion_pos_array) in &anion_list {
            for (anion_pos, anion_type) in anion_pos_array {
                for (cation_path, cation_pos_array) in &cation_list {
                    for (cation_pos, cation_type) in cation_pos_array {
                        let d = anion_pos.distance_from(cation_pos);
                        if d < Self::MAX_DISTANCE {
                            ion_pairs.push(IonPairRecord::new(
                                anion_path,
                                cation_path,
                                anion_type,
                                cation_type,
                            ));
                        }
                    }
                }
            }
        }

        ion_pairs
    }

    /// Convenience static helper to detect ion pairs directly from an `AtomGroup`.
    pub fn find_ion_pairs(model: &AtomGroup) -> Vec<IonPairRecord> {
        let detector = Self::new(model);
        detector.get_ion_pairs()
    }

    /// Extracts anion and cation interaction sites from the structure.
    ///
    /// NOTE on divergence from Python:
    /// In the original Python implementation (`ionpair.py`), missing expected atoms
    /// (e.g., missing `OE2` in GLU or `CZ` in ARG) raises a `KeyError` during dictionary
    /// lookup, which crashes the entire ion-pair analysis.
    /// In this Rust implementation, `get_center_*` functions return `Option<Position>`.
    /// When any required atom is missing, the incomplete residue site is safely skipped
    /// (`None`), allowing the analysis of other intact residues in the model to proceed
    /// without crashing.
    fn get_ion_list(&self) -> (ResidueSites, ResidueSites) {
        let mut anion_list: ResidueSites = Vec::new();
        let mut cation_list: ResidueSites = Vec::new();

        for (_chain_key, chain) in self.model.groups() {
            for (_res_key, res) in chain.groups() {
                let name = res.name.trim();

                if name == "GLU" {
                    if let Some(pos) = Self::get_center_glu(res) {
                        anion_list.push((res.path().to_string(), vec![(pos, "GLU".to_string())]));
                    }
                } else if name == "ASP" {
                    if let Some(pos) = Self::get_center_asp(res) {
                        anion_list.push((res.path().to_string(), vec![(pos, "ASP".to_string())]));
                    }
                } else if name == "LYS" {
                    if let Some(pos) = Self::get_center_lys(res) {
                        cation_list.push((res.path().to_string(), vec![(pos, "LYS".to_string())]));
                    }
                } else if name == "ARG" {
                    let mut arg_sites = Vec::new();
                    if let Some(pos0) = Self::get_center_arg(res, 0) {
                        arg_sites.push((pos0, "ARG".to_string()));
                    }
                    if let Some(pos1) = Self::get_center_arg(res, 1) {
                        arg_sites.push((pos1, "ARG1".to_string()));
                    }
                    if let Some(pos2) = Self::get_center_arg(res, 2) {
                        arg_sites.push((pos2, "ARG2".to_string()));
                    }
                    if !arg_sites.is_empty() {
                        cation_list.push((res.path().to_string(), arg_sites));
                    }
                }

                if AminoAcid::is_aminoacid(res) {
                    if res.has_atom("H3") {
                        if let Some(pos) = Self::get_center_nterm(res) {
                            if let Some(entry) =
                                cation_list.iter_mut().find(|(p, _)| p == res.path())
                            {
                                entry.1.push((pos, "NTM".to_string()));
                            } else {
                                cation_list
                                    .push((res.path().to_string(), vec![(pos, "NTM".to_string())]));
                            }
                        }
                    }
                    if res.has_atom("OXT") {
                        if let Some(pos) = Self::get_center_cterm(res) {
                            if let Some(entry) =
                                anion_list.iter_mut().find(|(p, _)| p == res.path())
                            {
                                entry.1.push((pos, "CTM".to_string()));
                            } else {
                                anion_list
                                    .push((res.path().to_string(), vec![(pos, "CTM".to_string())]));
                            }
                        }
                    }
                }
            }
        }

        (anion_list, cation_list)
    }

    /// Returns the position of N if present; returns None if missing.
    fn get_center_nterm(res: &AtomGroup) -> Option<Position> {
        res.get_atom("N").map(|a| a.xyz)
    }

    /// Computes the geometric center of C, O, OXT. Returns None if any atom is missing.
    fn get_center_cterm(res: &AtomGroup) -> Option<Position> {
        let c = res.get_atom("C")?;
        let o = res.get_atom("O")?;
        let oxt = res.get_atom("OXT")?;
        Some((c.xyz + o.xyz + oxt.xyz) / 3.0)
    }

    /// Computes the geometric center of CD, OE1, OE2. Returns None if any atom is missing.
    fn get_center_glu(res: &AtomGroup) -> Option<Position> {
        let cd = res.get_atom("CD")?;
        let oe1 = res.get_atom("OE1")?;
        let oe2 = res.get_atom("OE2")?;
        Some((cd.xyz + oe1.xyz + oe2.xyz) / 3.0)
    }

    /// Computes the geometric center of CG, OD1, OD2. Returns None if any atom is missing.
    fn get_center_asp(res: &AtomGroup) -> Option<Position> {
        let cg = res.get_atom("CG")?;
        let od1 = res.get_atom("OD1")?;
        let od2 = res.get_atom("OD2")?;
        Some((cg.xyz + od1.xyz + od2.xyz) / 3.0)
    }

    /// Returns the position of NZ if present; returns None if missing.
    fn get_center_lys(res: &AtomGroup) -> Option<Position> {
        res.get_atom("NZ").map(|a| a.xyz)
    }

    /// Computes ARG site position (case 0: center of NH1, NH2, CZ; case 1: NH1; case 2: NH2).
    /// Returns None if any required atom is missing.
    fn get_center_arg(res: &AtomGroup, case: usize) -> Option<Position> {
        match case {
            0 => {
                let nh1 = res.get_atom("NH1")?;
                let nh2 = res.get_atom("NH2")?;
                let cz = res.get_atom("CZ")?;
                Some((nh1.xyz + nh2.xyz + cz.xyz) / 3.0)
            }
            1 => res.get_atom("NH1").map(|a| a.xyz),
            2 => res.get_atom("NH2").map(|a| a.xyz),
            _ => None,
        }
    }
}
