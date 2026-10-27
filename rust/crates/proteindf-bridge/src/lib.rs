// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

pub mod amino_acid;
pub mod atom;
pub mod atom_group;
pub mod backbone_hydrogen;
pub mod bond;
pub mod brd;
pub mod ccd_templates;
pub mod ch_pi;
pub mod complete_atoms;
pub mod error;
pub mod format;
pub mod his_state;
pub mod hydrogen_bond;
pub mod hydrogenation;
pub mod ion_pair;
pub mod matrix;
pub mod missing_atoms;
pub mod modeling;
pub mod neutralize;
pub mod orchestrator;
pub mod periodic_table;
pub mod position;
pub mod ramachandran;
pub mod sasa;
pub mod secondary_structure;
pub mod selector;
pub mod spatial;
pub mod ssbond;
pub mod superposer;
pub mod superposer_quaternion;
pub mod vector;

pub mod qc_result;

pub use amino_acid::AminoAcid;
pub use atom::Atom;
pub use atom_group::{AtomGroup, BondRecord, SchemaViolation, Selector};
pub use backbone_hydrogen::{
    add_backbone_hydrogens_to_residue, add_backbone_hydrogens_to_residue_in_place,
    build_backbone_amide_hydrogen, build_nterm_hydrogens, STANDARD_AMIDE_NH_BOND_LENGTH,
    STANDARD_NTERM_NH_BOND_LENGTH,
};
pub use bond::Bond;
pub use brd::{
    load_atomgroup, load_brd_yui, load_msgpack, save_atomgroup, save_brd_yui, save_msgpack,
};
pub use ccd_templates::{CcdBondTemplate, CcdTemplateDb};
pub use ch_pi::{
    calc_ch_pi_interactions, calc_ch_pi_interactions_with_thresholds, calc_ring_center_and_normal,
    calc_ring_geometry, AromaticRing, AromaticRingDef, ChPiInteraction, AROMATIC_RINGS,
    DEFAULT_MAX_ANGLE_DEG, DEFAULT_MAX_DISTANCE,
};
pub use complete_atoms::{
    complete_missing_atoms, CompleteOptions, CompletionReport, CompletionStatus, ResidueCompletion,
};
pub use error::{BridgeError, Result};
pub use format::{
    AmberPrmtop, Cube, CubeAtom, CubeUnit, Format, MmcifStructureReport, MmcifWriteOptions, Pdb,
    SimpleGro, SimpleMmcif, SimpleMol2, StructConnPartnerUnresolved, StructConnRecord,
    UnresolvedStructConn, Xyz, BOHR_TO_ANGSTROM,
};
pub use his_state::{
    his_states, HisDecidedBy, HisReport, HisSiteReport, HisState, HisStateOptions,
    COORDINATING_METALS, DEFAULT_BURIAL_RADIUS, DEFAULT_METAL_DISTANCE, SALT_BRIDGE_DISTANCE,
};
pub use hydrogen_bond::{
    calc_backbone_hbonds, calc_kabsch_sander_energy, calc_pseudo_hydrogen, calc_sidechain_hbonds,
    calc_sidechain_hbonds_with_options, HydrogenBond, SidechainAtomType, SidechainHydrogenBond,
    SIDECHAIN_ATOM_TYPES,
};
pub use hydrogenation::{
    add_hydrogens_to_component, add_hydrogens_to_component_in_place,
    add_hydrogens_to_component_in_place_with_options, add_hydrogens_to_component_with_options,
    HydrogenationOptions, HydrogenationReport, MIN_SUPERPOSE_HEAVY_ATOMS,
};
pub use ion_pair::{IonPair, IonPairRecord};
pub use matrix::{identity_matrix, Matrix, SymmetricMatrix};
pub use missing_atoms::{
    check_missing_atoms, MissingAtomsOptions, MissingAtomsReport, ResidueAtomCheck, ResidueStatus,
};
pub use modeling::Modeling;
pub use neutralize::Neutralize;
pub use orchestrator::{
    hydrogenate_atomgroup, hydrogenate_atomgroup_with_his, OverallHydrogenationReport,
    DEFAULT_HIS_MODE,
};
pub use periodic_table::PeriodicTable;
pub use position::{dihedral_angle, Position};
pub use ramachandran::{calc_phi_psi, RamachandranAngle, PEPTIDE_BOND_MAX_DISTANCE};
pub use sasa::{AtomSasa, ResidueSasa, SasaOptions, SasaReport};
pub use secondary_structure::{
    apply_secondary_structure, calc_secondary_structure, SecondaryStructure, SsCode,
};
pub use selector::{
    SelectAtom, SelectAtomGroup, SelectName, SelectPath, SelectPathRegex, SelectPathSimple,
    SelectPathWildcard, SelectRange, SelectSymbol, Select_Atom, Select_AtomGroup, Select_Name,
    Select_Path, Select_PathRegex, Select_Path_simple, Select_Path_wildcard, Select_Range,
    Select_Symbol,
};
pub use spatial::CellList;
pub use ssbond::SSBond;
pub use superposer::Superposer;
pub use superposer_quaternion::{
    SuperposerQuaternion, SuperposerQuaternion as Superposer_quaternion,
};
pub use vector::Vector;

pub mod interaction_set;
pub use interaction_set::{Interaction, InteractionKind, InteractionSet};

pub use qc_result::{
    calc_orbital_mulliken_contributions, AtomicOrbital, DensityMatrix, FockMatrix, MoCoefficients,
    MullikenPopulation, OrbitalContribution, OrbitalEnergy, OrbitalInfo, OrbitalOccupation, QcAtom,
    QcBasisSet, QcBasisSetDefinition, QcCgto, QcMethod, QcMolecule, QcParameters, QcPgto, QcResult,
    RunType,
};

#[cfg(feature = "hdf5")]
pub use qc_result::PdfParamHdf5;
