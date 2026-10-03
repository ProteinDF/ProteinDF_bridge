// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::atom_group::PyAtomGroup;
use crate::error::to_py_err;
use proteindf_bridge::neutralize::Neutralize as CoreNeutralize;
use pyo3::prelude::*;

/// Neutralizes protein charges by adding counter-ions (Na+ / Cl-) to charged residues and termini.
/// Corresponding to `proteindf_bridge.neutralize.Neutralize`.
#[pyclass(
    name = "Neutralize",
    module = "proteindf_bridge_rs",
    skip_from_py_object
)]
pub struct PyNeutralize {
    neutralized_obj: Py<PyAtomGroup>,
}

#[pymethods]
impl PyNeutralize {
    /// Creates a new `Neutralize` instance and neutralizes the given protein model.
    #[new]
    #[pyo3(signature = (protein))]
    pub fn new(py: Python<'_>, protein: &PyAtomGroup) -> PyResult<Self> {
        let neut = CoreNeutralize::new(&protein.inner).map_err(to_py_err)?;
        let py_ag = PyAtomGroup::from_core(neut.into_neutralized());
        let neutralized_obj = Py::new(py, py_ag)?;
        Ok(Self { neutralized_obj })
    }

    /// Returns the neutralized `AtomGroup`.
    ///
    /// Per project convention, returns the same Python object on each access without re-copying.
    #[getter]
    pub fn neutralized(&self, py: Python<'_>) -> Py<PyAtomGroup> {
        self.neutralized_obj.clone_ref(py)
    }

    /// Computes the exempt list for a model using ion pair detection.
    ///
    /// Note: Each access returns a new list to prevent external mutation from affecting internal state.
    #[pyo3(signature = (model))]
    pub fn _exempt_list(&self, model: &PyAtomGroup) -> Vec<(String, String, String)> {
        CoreNeutralize::exempt_list(&model.inner)
    }

    /// Public alias for `_exempt_list`.
    #[pyo3(signature = (model))]
    pub fn exempt_list(&self, model: &PyAtomGroup) -> Vec<(String, String, String)> {
        self._exempt_list(model)
    }

    /// Divides an `AtomGroup` path into `(chain_name, res_name)`.
    #[pyo3(signature = (path))]
    pub fn _divide_path(&self, path: &str) -> (String, String) {
        CoreNeutralize::divide_path(path)
    }

    /// Public alias for `_divide_path`.
    #[pyo3(signature = (path))]
    pub fn divide_path(&self, path: &str) -> (String, String) {
        self._divide_path(path)
    }

    fn __repr__(&self, py: Python<'_>) -> String {
        format!(
            "Neutralize(atoms={})",
            self.neutralized_obj
                .borrow(py)
                .inner
                .get_number_of_all_atoms()
        )
    }
}
