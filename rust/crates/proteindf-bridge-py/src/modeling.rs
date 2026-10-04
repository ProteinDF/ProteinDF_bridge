// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

#![allow(non_snake_case)]

use crate::atom::PyAtom;
use crate::atom_group::PyAtomGroup;
use crate::error::to_py_err;
use crate::matrix::PyMatrix;
use crate::position::extract_position;
use proteindf_bridge::modeling::Modeling as CoreModeling;
use pyo3::prelude::*;

/// Structural modeling and capping/neutralization utility corresponding to `proteindf_bridge.modeling.Modeling`.
#[pyclass(name = "Modeling", module = "proteindf_bridge.rs", skip_from_py_object)]
pub struct PyModeling {
    pub(crate) inner: CoreModeling,
}

#[pymethods]
impl PyModeling {
    /// Creates a new `Modeling` instance with embedded reference conformer structures.
    #[new]
    pub fn new() -> PyResult<Self> {
        let inner = CoreModeling::new().map_err(to_py_err)?;
        Ok(Self { inner })
    }

    /// Creates a new `Modeling` instance loading reference conformers from a directory.
    #[staticmethod]
    #[pyo3(signature = (data_dir))]
    pub fn from_data_dir(data_dir: &str) -> PyResult<Self> {
        let inner = CoreModeling::from_data_dir(data_dir).map_err(to_py_err)?;
        Ok(Self { inner })
    }

    /// Builds an ACE capping group by fitting against the optimal reference conformer.
    #[pyo3(signature = (res, next_aa = None))]
    pub fn get_ACE(
        &self,
        res: &PyAtomGroup,
        next_aa: Option<&PyAtomGroup>,
    ) -> PyResult<PyAtomGroup> {
        let next_inner = next_aa.map(|n| &n.inner);
        let ag = self
            .inner
            .get_ACE(&res.inner, next_inner)
            .map_err(to_py_err)?;
        Ok(PyAtomGroup::from_core(ag))
    }

    /// Builds an NME capping group by fitting against the optimal reference conformer.
    #[pyo3(signature = (res, next_aa = None))]
    pub fn get_NME(
        &self,
        res: &PyAtomGroup,
        next_aa: Option<&PyAtomGroup>,
    ) -> PyResult<PyAtomGroup> {
        let next_inner = next_aa.map(|n| &n.inner);
        let ag = self
            .inner
            .get_NME(&res.inner, next_inner)
            .map_err(to_py_err)?;
        Ok(PyAtomGroup::from_core(ag))
    }

    /// Turns the neighboring C-alpha position into an ACE methyl group.
    #[pyo3(signature = (next_aa))]
    pub fn get_ACE_simple(&self, next_aa: &PyAtomGroup) -> PyResult<PyAtomGroup> {
        let ag = self
            .inner
            .get_ACE_simple(&next_aa.inner)
            .map_err(to_py_err)?;
        Ok(PyAtomGroup::from_core(ag))
    }

    /// Turns the neighboring C-alpha position into an NME methyl group.
    #[pyo3(signature = (next_aa))]
    pub fn get_NME_simple(&self, next_aa: &PyAtomGroup) -> PyResult<PyAtomGroup> {
        let ag = self
            .inner
            .get_NME_simple(&next_aa.inner)
            .map_err(to_py_err)?;
        Ok(PyAtomGroup::from_core(ag))
    }

    /// Adds the hydrogens of -CH3 to C1, oriented away from C2.
    #[pyo3(signature = (C1, C2))]
    #[allow(non_snake_case)]
    pub fn add_methyl(&self, C1: &PyAtom, C2: &PyAtom) -> PyResult<PyAtomGroup> {
        let ag = self
            .inner
            .add_methyl(&C1.inner, &C2.inner)
            .map_err(to_py_err)?;
        Ok(PyAtomGroup::from_core(ag))
    }

    /// Generates an NH3 atom group with given angle and bond length.
    #[pyo3(signature = (angle = std::f64::consts::FRAC_PI_2, length = 1.0))]
    pub fn get_NH3(&self, angle: f64, length: f64) -> PyResult<PyAtomGroup> {
        let ag = self.inner.get_NH3(angle, length).map_err(to_py_err)?;
        Ok(PyAtomGroup::from_core(ag))
    }

    /// Returns consecutive amino acid residues from a chain group.
    #[pyo3(signature = (chain, from_resid, to_resid))]
    pub fn select_residues(
        &self,
        chain: &PyAtomGroup,
        from_resid: i64,
        to_resid: i64,
    ) -> PyAtomGroup {
        let ag = self
            .inner
            .select_residues(&chain.inner, from_resid, to_resid);
        PyAtomGroup::from_core(ag)
    }

    /// Computes the 3x3 rotation matrix `R` such that applying `R` to `in_b` yields
    /// a vector pointing in the direction of `in_a`.
    #[pyo3(signature = (in_a, in_b))]
    pub fn arbitary_rotate_matrix(
        &self,
        in_a: &Bound<'_, PyAny>,
        in_b: &Bound<'_, PyAny>,
    ) -> PyResult<PyMatrix> {
        let pos_a = extract_position(in_a)?;
        let pos_b = extract_position(in_b)?;
        let mat = self
            .inner
            .arbitary_rotate_matrix(pos_a, pos_b)
            .map_err(to_py_err)?;
        Ok(PyMatrix::from_core(mat))
    }

    /// Returns the maximum numerical index found in atom keys.
    #[pyo3(signature = (res))]
    pub fn get_last_index(&self, res: &PyAtomGroup) -> usize {
        self.inner.get_last_index(&res.inner)
    }

    /// Returns a Cl- group to neutralize the N-terminal residue.
    #[pyo3(signature = (res))]
    pub fn neutralize_Nterm(&self, res: &PyAtomGroup) -> PyResult<PyAtomGroup> {
        let ag = self.inner.neutralize_Nterm(&res.inner).map_err(to_py_err)?;
        Ok(PyAtomGroup::from_core(ag))
    }

    /// Returns an Na+ group to neutralize the C-terminal residue.
    #[pyo3(signature = (res))]
    pub fn neutralize_Cterm(&self, res: &PyAtomGroup) -> PyResult<PyAtomGroup> {
        let ag = self.inner.neutralize_Cterm(&res.inner).map_err(to_py_err)?;
        Ok(PyAtomGroup::from_core(ag))
    }

    /// Returns an Na+ group to neutralize a GLU side chain.
    #[pyo3(signature = (res))]
    pub fn neutralize_GLU(&self, res: &PyAtomGroup) -> PyResult<PyAtomGroup> {
        let ag = self.inner.neutralize_GLU(&res.inner).map_err(to_py_err)?;
        Ok(PyAtomGroup::from_core(ag))
    }

    /// Returns an Na+ group to neutralize an ASP side chain.
    #[pyo3(signature = (res))]
    pub fn neutralize_ASP(&self, res: &PyAtomGroup) -> PyResult<PyAtomGroup> {
        let ag = self.inner.neutralize_ASP(&res.inner).map_err(to_py_err)?;
        Ok(PyAtomGroup::from_core(ag))
    }

    /// Returns a Cl- group to neutralize a LYS side chain.
    #[pyo3(signature = (res))]
    pub fn neutralize_LYS(&self, res: &PyAtomGroup) -> PyResult<PyAtomGroup> {
        let ag = self.inner.neutralize_LYS(&res.inner).map_err(to_py_err)?;
        Ok(PyAtomGroup::from_core(ag))
    }

    /// Returns a Cl- group to neutralize an ARG side chain.
    /// case: 0 (center), 1 (NH1 side), 2 (NH2 side).
    #[pyo3(signature = (res, case = 0))]
    pub fn neutralize_ARG(&self, res: &PyAtomGroup, case: usize) -> PyResult<PyAtomGroup> {
        let ag = self
            .inner
            .neutralize_ARG(&res.inner, case)
            .map_err(to_py_err)?;
        Ok(PyAtomGroup::from_core(ag))
    }

    /// Neutralizes FAD phosphate groups by adding two Na+ ions.
    #[pyo3(signature = (ag))]
    pub fn neutralize_FAD(&self, ag: &PyAtomGroup) -> PyResult<PyAtomGroup> {
        let res_ag = self.inner.neutralize_FAD(&ag.inner).map_err(to_py_err)?;
        Ok(PyAtomGroup::from_core(res_ag))
    }

    fn __repr__(&self) -> String {
        "Modeling()".to_string()
    }
}
