// Copyright (C) 2014 The ProteinDF development team.
// see also AUTHORS and README if provided.
//
// This file is a part of the ProteinDF software package.

use pyo3::prelude::*;

pub mod atom;
pub mod atom_group;
pub mod bond;
pub mod error;
pub mod matrix;
pub mod periodic_table;
pub mod position;
pub mod vector;

#[pymodule]
fn proteindf_bridge_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Exceptions
    m.add("BrError", m.py().get_type::<error::BrError>())?;
    m.add("BrInputError", m.py().get_type::<error::BrInputError>())?;
    m.add("BrValueError", m.py().get_type::<error::BrValueError>())?;

    // Classes
    m.add_class::<periodic_table::PyPeriodicTable>()?;
    m.add_class::<vector::PyVector>()?;
    m.add_class::<matrix::PyMatrix>()?;
    m.add_class::<matrix::PySymmetricMatrix>()?;
    m.add_class::<position::PyPosition>()?;
    m.add_class::<atom::PyAtom>()?;
    m.add_class::<bond::PyBond>()?;
    m.add_class::<atom_group::PyAtomGroup>()?;

    Ok(())
}
