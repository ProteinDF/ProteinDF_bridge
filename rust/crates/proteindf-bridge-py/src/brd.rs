// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::atom_group::PyAtomGroup;
use crate::error::to_py_err;
use proteindf_bridge::brd as core_brd;
use pyo3::prelude::*;

/// Loads a bridge (.brd) MessagePack file directly into an AtomGroup.
#[pyfunction]
#[pyo3(signature = (brd_path))]
pub fn load_atomgroup(brd_path: &str) -> PyResult<PyAtomGroup> {
    let group = core_brd::load_atomgroup(brd_path).map_err(to_py_err)?;
    Ok(PyAtomGroup::from_core(group))
}

/// Saves an AtomGroup into a plain bridge (.brd) MessagePack file.
#[pyfunction]
#[pyo3(signature = (atomgroup, file_path))]
pub fn save_atomgroup(atomgroup: &PyAtomGroup, file_path: &str) -> PyResult<()> {
    core_brd::save_atomgroup(&atomgroup.inner, file_path).map_err(to_py_err)
}

/// Loads an AtomGroup from a file with a YUI-compatible header (Magic + Version + zstd/none).
#[pyfunction]
#[pyo3(signature = (path))]
pub fn load_brd_yui(path: &str) -> PyResult<PyAtomGroup> {
    let group = core_brd::load_brd_yui(path).map_err(to_py_err)?;
    Ok(PyAtomGroup::from_core(group))
}

/// Saves an AtomGroup with a YUI-compatible header (Magic + Version + zstd/none).
#[pyfunction]
#[pyo3(signature = (atomgroup, path, compress_zstd = false))]
pub fn save_brd_yui(atomgroup: &PyAtomGroup, path: &str, compress_zstd: bool) -> PyResult<()> {
    core_brd::save_brd_yui(&atomgroup.inner, path, compress_zstd).map_err(to_py_err)
}
