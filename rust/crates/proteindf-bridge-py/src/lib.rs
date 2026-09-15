// Copyright (C) 2014 The ProteinDF development team.
// see also AUTHORS and README if provided.
//
// This file is a part of the ProteinDF software package.

use pyo3::prelude::*;

pub mod amino_acid;
pub mod atom;
pub mod atom_group;
pub mod bond;
pub mod error;
pub mod format;
pub mod ion_pair;
pub mod matrix;
pub mod periodic_table;
pub mod position;
pub mod selector;
pub mod ssbond;
pub mod superposer;
pub mod superposer_quaternion;
pub mod vector;

#[pymodule]
fn proteindf_bridge_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Exceptions
    m.add("BrError", m.py().get_type::<error::BrError>())?;
    m.add("BrInputError", m.py().get_type::<error::BrInputError>())?;
    m.add("BrValueError", m.py().get_type::<error::BrValueError>())?;

    // Foundation & Data Model (PR#14)
    m.add_class::<periodic_table::PyPeriodicTable>()?;
    m.add_class::<vector::PyVector>()?;
    m.add_class::<matrix::PyMatrix>()?;
    m.add_class::<matrix::PySymmetricMatrix>()?;
    m.add_class::<position::PyPosition>()?;
    m.add_class::<atom::PyAtom>()?;
    m.add_class::<bond::PyBond>()?;
    m.add_class::<atom_group::PyAtomGroup>()?;

    // Format I/O (PR#15)
    m.add_class::<format::PyFormat>()?;
    m.add_class::<format::PyXyz>()?;
    m.add_class::<format::PySimpleGro>()?;
    m.add_class::<format::PySimpleMol2>()?;
    m.add_class::<format::PyAmberPrmtop>()?;
    m.add_class::<format::PyPdb>()?;
    m.add_class::<format::PySimpleMmcif>()?;

    // Structural Operations (PR#16)
    m.add_class::<amino_acid::PyAminoAcid>()?;
    m.add_class::<ssbond::PySSBond>()?;
    m.add_class::<ion_pair::PyIonPair>()?;
    m.add_class::<superposer::PySuperposer>()?;
    m.add_class::<superposer_quaternion::PySuperposerQuaternion>()?;
    m.add("SuperposerQuaternion", m.getattr("Superposer_quaternion")?)?;

    // Selectors
    m.add_class::<selector::PySelectSymbol>()?;
    m.add("SelectSymbol", m.getattr("Select_Symbol")?)?;

    m.add_class::<selector::PySelectName>()?;
    m.add("SelectName", m.getattr("Select_Name")?)?;

    m.add_class::<selector::PySelectPathSimple>()?;
    m.add("SelectPathSimple", m.getattr("Select_Path_simple")?)?;

    m.add_class::<selector::PySelectPathWildcard>()?;
    m.add("SelectPathWildcard", m.getattr("Select_Path_wildcard")?)?;

    m.add_class::<selector::PySelectPathRegex>()?;
    m.add("SelectPathRegex", m.getattr("Select_PathRegex")?)?;

    m.add_class::<selector::PySelectPath>()?;
    m.add("SelectPath", m.getattr("Select_Path")?)?;

    m.add_class::<selector::PySelectRange>()?;
    m.add("SelectRange", m.getattr("Select_Range")?)?;

    m.add_class::<selector::PySelectAtom>()?;
    m.add("SelectAtom", m.getattr("Select_Atom")?)?;

    m.add_class::<selector::PySelectAtomGroup>()?;
    m.add("SelectAtomGroup", m.getattr("Select_AtomGroup")?)?;

    Ok(())
}
