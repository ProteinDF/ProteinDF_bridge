// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use std::fs::File;
use std::io::{Read, Write};
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::atom_group::AtomGroup;
use crate::brd::sort_nicely;
use crate::ch_pi::{
    calc_ch_pi_interactions_with_thresholds, DEFAULT_MAX_ANGLE_DEG, DEFAULT_MAX_DISTANCE,
};
use crate::error::{BridgeError, Result};
use crate::hydrogen_bond::{calc_backbone_hbonds, calc_sidechain_hbonds};
use crate::ion_pair::IonPair;
use crate::ssbond::SSBond;

/// Kind of molecular interaction.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum InteractionKind {
    /// Disulfide bond between cysteine residues.
    Disulfide,
    /// Salt bridge / ion pair between acidic and basic residues.
    SaltBridge,
    /// Hydrogen bond (backbone or sidechain).
    HydrogenBond,
    /// CH-pi interaction between aliphatic/carbon atom and aromatic ring.
    ChPi,
}

/// Represents a detected non-covalent or disulfide interaction.
///
/// Matches the schema defined in `RUST_PORT_SPEC.md` §3.3:
/// `{ kind, atoms, distance, angle, donor_acceptor_role }`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Interaction {
    /// The classification of the interaction.
    pub kind: InteractionKind,
    /// Participating atom or residue paths (variable length, e.g. 2 for pairwise).
    pub atoms: Vec<String>,
    /// Interaction distance in Angstroms, if measured.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub distance: Option<f64>,
    /// Interaction angle in degrees, if measured.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub angle: Option<f64>,
    /// Specific role or description (e.g. "backbone", "sidechain", "GLU - LYS", "CH-pi (PHE)").
    #[serde(skip_serializing_if = "Option::is_none")]
    pub donor_acceptor_role: Option<String>,
}

impl Interaction {
    /// Creates a new `Interaction` with mandatory fields and `None` for optional metrics.
    pub fn new(kind: InteractionKind, atoms: Vec<String>) -> Self {
        Self {
            kind,
            atoms,
            distance: None,
            angle: None,
            donor_acceptor_role: None,
        }
    }

    /// Creates a new `Interaction` with all fields specified.
    pub fn with_metrics(
        kind: InteractionKind,
        atoms: Vec<String>,
        distance: Option<f64>,
        angle: Option<f64>,
        donor_acceptor_role: Option<String>,
    ) -> Self {
        Self {
            kind,
            atoms,
            distance,
            angle,
            donor_acceptor_role,
        }
    }
}

/// Aggregation set of detected interactions in a molecular structure.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct InteractionSet {
    /// List of all detected interactions.
    pub interactions: Vec<Interaction>,
}

fn get_atom_from_model<'a>(model: &'a AtomGroup, path: &str) -> Option<&'a crate::atom::Atom> {
    let model_path = model.path();
    let rel_path = if !model_path.is_empty() && path.starts_with(model_path) {
        &path[model_path.len()..]
    } else {
        path.trim_start_matches('/')
    };
    model.get_atom_by_path(rel_path)
}

impl InteractionSet {
    /// Creates an empty `InteractionSet`.
    pub fn new() -> Self {
        Self {
            interactions: Vec::new(),
        }
    }

    /// Creates an `InteractionSet` from a vector of interactions.
    pub fn with_interactions(interactions: Vec<Interaction>) -> Self {
        Self { interactions }
    }

    /// Number of interactions in the set.
    pub fn len(&self) -> usize {
        self.interactions.len()
    }

    /// Returns `true` if there are no interactions.
    pub fn is_empty(&self) -> bool {
        self.interactions.is_empty()
    }

    /// Iterator over interactions.
    pub fn iter(&self) -> std::slice::Iter<'_, Interaction> {
        self.interactions.iter()
    }

    /// Adds an interaction to the set.
    pub fn push(&mut self, interaction: Interaction) {
        self.interactions.push(interaction);
    }

    /// Filters interactions matching a specific `InteractionKind`.
    pub fn filter_by_kind(&self, kind: InteractionKind) -> Vec<&Interaction> {
        self.interactions
            .iter()
            .filter(|i| i.kind == kind)
            .collect()
    }

    /// Counts interactions of a given `InteractionKind`.
    pub fn count_by_kind(&self, kind: InteractionKind) -> usize {
        self.interactions.iter().filter(|i| i.kind == kind).count()
    }

    /// Detects all supported interactions in the given model:
    /// 1. Disulfide bonds (`SSBond::find_bonds`)
    /// 2. Salt bridges (`IonPair::find_ion_pairs`)
    /// 3. Hydrogen bonds (backbone `calc_backbone_hbonds` + sidechain `calc_sidechain_hbonds`)
    /// 4. CH-pi interactions (`calc_ch_pi_interactions_with_thresholds`)
    ///
    /// The CH-pi distance and angle thresholds can be customized via `ch_pi_max_distance`
    /// and `ch_pi_max_angle_deg`. If `None`, default thresholds (4.5 Å, 40.0°) are used.
    pub fn detect_all(
        model: &AtomGroup,
        ch_pi_max_distance: Option<f64>,
        ch_pi_max_angle_deg: Option<f64>,
    ) -> Self {
        let mut set = Self::new();

        // 1. Disulfide bonds
        let ssbonds = SSBond::find_bonds(model);
        for (path1, path2) in ssbonds {
            // Compute SG-SG distance if atoms are accessible
            let sg1_path = format!("{path1}SG");
            let sg2_path = format!("{path2}SG");
            let sg1_pos = get_atom_from_model(model, &sg1_path).map(|a| a.xyz);
            let sg2_pos = get_atom_from_model(model, &sg2_path).map(|a| a.xyz);
            let dist = match (sg1_pos, sg2_pos) {
                (Some(p1), Some(p2)) => Some(p1.distance_from(&p2)),
                _ => None,
            };

            set.push(Interaction::with_metrics(
                InteractionKind::Disulfide,
                vec![sg1_path, sg2_path],
                dist,
                None,
                Some("disulfide".to_string()),
            ));
        }

        // 2. Salt bridges
        let ion_pairs = IonPair::find_ion_pairs(model);
        for ip in ion_pairs {
            set.push(Interaction::with_metrics(
                InteractionKind::SaltBridge,
                vec![ip.anion_path.clone(), ip.cation_path.clone()],
                None,
                None,
                Some(format!("{} - {}", ip.anion_type, ip.cation_type)),
            ));
        }

        // 3. Backbone hydrogen bonds (per chain)
        let mut chain_keys = model.get_group_list();
        sort_nicely(&mut chain_keys);
        for c_key in &chain_keys {
            if let Some(chain) = model.get_group(c_key) {
                let bb_hbonds = calc_backbone_hbonds(chain);
                let chain_path = chain.path().to_string();
                for hb in bb_hbonds {
                    let d_path = format!("{}{}/N", chain_path, hb.donor_residue_key);
                    let a_path = format!("{}{}/O", chain_path, hb.acceptor_residue_key);

                    let d_pos = get_atom_from_model(model, &d_path).map(|a| a.xyz);
                    let a_pos = get_atom_from_model(model, &a_path).map(|a| a.xyz);
                    let dist = match (d_pos, a_pos) {
                        (Some(p1), Some(p2)) => Some(p1.distance_from(&p2)),
                        _ => None,
                    };

                    set.push(Interaction::with_metrics(
                        InteractionKind::HydrogenBond,
                        vec![d_path, a_path],
                        dist,
                        None,
                        Some(format!("backbone (E={:.3} kcal/mol)", hb.energy)),
                    ));
                }
            }
        }

        // 4. Sidechain hydrogen bonds
        let sc_hbonds = calc_sidechain_hbonds(model);
        for sc in sc_hbonds {
            set.push(Interaction::with_metrics(
                InteractionKind::HydrogenBond,
                vec![
                    format!("{}{}", sc.donor_path, sc.donor_atom),
                    format!("{}{}", sc.acceptor_path, sc.acceptor_atom),
                ],
                Some(sc.distance),
                sc.angle,
                Some("sidechain".to_string()),
            ));
        }

        // 5. CH-pi interactions
        let dist_thresh = ch_pi_max_distance.unwrap_or(DEFAULT_MAX_DISTANCE);
        let angle_thresh = ch_pi_max_angle_deg.unwrap_or(DEFAULT_MAX_ANGLE_DEG);
        let ch_pi_list = calc_ch_pi_interactions_with_thresholds(model, dist_thresh, angle_thresh);
        for ch in ch_pi_list {
            set.push(Interaction::with_metrics(
                InteractionKind::ChPi,
                vec![
                    format!("{}{}", ch.carbon_path, ch.carbon_atom),
                    ch.ring_path.clone(),
                ],
                Some(ch.distance),
                Some(ch.angle),
                Some(format!("CH-pi ({})", ch.ring_residue)),
            ));
        }

        set
    }

    /// Serializes `InteractionSet` to MessagePack byte vector.
    pub fn to_msgpack(&self) -> Result<Vec<u8>> {
        rmp_serde::to_vec_named(self).map_err(|e| BridgeError::MsgPack(e.to_string()))
    }

    /// Deserializes `InteractionSet` from MessagePack byte slice.
    pub fn from_msgpack(bytes: &[u8]) -> Result<Self> {
        rmp_serde::from_slice(bytes).map_err(|e| BridgeError::MsgPack(e.to_string()))
    }

    /// Saves `InteractionSet` as a MessagePack file.
    pub fn save_msgpack<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let bytes = self.to_msgpack()?;
        let mut file = File::create(path).map_err(|e| BridgeError::Io(e.to_string()))?;
        file.write_all(&bytes)
            .map_err(|e| BridgeError::Io(e.to_string()))?;
        Ok(())
    }

    /// Loads `InteractionSet` from a MessagePack file.
    pub fn load_msgpack<P: AsRef<Path>>(path: P) -> Result<Self> {
        let mut file = File::open(path).map_err(|e| BridgeError::Io(e.to_string()))?;
        let mut bytes = Vec::new();
        file.read_to_end(&mut bytes)
            .map_err(|e| BridgeError::Io(e.to_string()))?;
        Self::from_msgpack(&bytes)
    }

    /// Serializes `InteractionSet` to a YAML string.
    pub fn to_yaml(&self) -> Result<String> {
        serde_yaml_ng::to_string(self).map_err(|e| BridgeError::Yaml(e.to_string()))
    }

    /// Deserializes `InteractionSet` from a YAML string.
    pub fn from_yaml(yaml_str: &str) -> Result<Self> {
        serde_yaml_ng::from_str(yaml_str).map_err(|e| BridgeError::Yaml(e.to_string()))
    }

    /// Saves `InteractionSet` as a YAML file.
    pub fn save_yaml<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let yaml_str = self.to_yaml()?;
        let mut file = File::create(path).map_err(|e| BridgeError::Io(e.to_string()))?;
        file.write_all(yaml_str.as_bytes())
            .map_err(|e| BridgeError::Io(e.to_string()))?;
        Ok(())
    }

    /// Loads `InteractionSet` from a YAML file.
    pub fn load_yaml<P: AsRef<Path>>(path: P) -> Result<Self> {
        let mut file = File::open(path).map_err(|e| BridgeError::Io(e.to_string()))?;
        let mut yaml_str = String::new();
        file.read_to_string(&mut yaml_str)
            .map_err(|e| BridgeError::Io(e.to_string()))?;
        Self::from_yaml(&yaml_str)
    }
}
