// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::atom_group::PyAtomGroup;
use crate::error::to_py_err;
use proteindf_bridge::error::BridgeError;
use proteindf_bridge::format::mmcif::{
    MmcifStructureReport as CoreMmcifStructureReport, SimpleMmcif as CoreSimpleMmcif,
    StructConnPartnerUnresolved as CoreStructConnPartnerUnresolved,
    UnresolvedStructConn as CoreUnresolvedStructConn,
};
use proteindf_bridge::format::mmcif_writer::{validate_for_mmcif_write, MmcifWriteOptions};
use pyo3::prelude::*;

/// Describes why a partner atom in a `_struct_conn` record could not be resolved.
#[pyclass(
    name = "StructConnPartnerUnresolved",
    module = "proteindf_bridge.rs",
    skip_from_py_object
)]
#[derive(Clone)]
pub struct PyStructConnPartnerUnresolved {
    pub(crate) inner: CoreStructConnPartnerUnresolved,
}

impl From<CoreStructConnPartnerUnresolved> for PyStructConnPartnerUnresolved {
    fn from(inner: CoreStructConnPartnerUnresolved) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl PyStructConnPartnerUnresolved {
    /// Returns the reason category:
    /// - `"MissingSeqId"`
    /// - `"ChainNotFound"`
    /// - `"ResidueNotFound"`
    /// - `"AtomNotFound"`
    #[getter]
    pub fn reason(&self) -> String {
        match &self.inner {
            CoreStructConnPartnerUnresolved::MissingSeqId => "MissingSeqId".to_string(),
            CoreStructConnPartnerUnresolved::ChainNotFound { .. } => "ChainNotFound".to_string(),
            CoreStructConnPartnerUnresolved::ResidueNotFound { .. } => {
                "ResidueNotFound".to_string()
            }
            CoreStructConnPartnerUnresolved::AtomNotFound { .. } => "AtomNotFound".to_string(),
        }
    }

    #[getter]
    pub fn chain_id(&self) -> Option<String> {
        match &self.inner {
            CoreStructConnPartnerUnresolved::ChainNotFound { chain_id } => Some(chain_id.clone()),
            CoreStructConnPartnerUnresolved::ResidueNotFound { chain_id, .. } => {
                Some(chain_id.clone())
            }
            CoreStructConnPartnerUnresolved::AtomNotFound { chain_id, .. } => {
                Some(chain_id.clone())
            }
            _ => None,
        }
    }

    #[getter]
    pub fn res_key(&self) -> Option<String> {
        match &self.inner {
            CoreStructConnPartnerUnresolved::ResidueNotFound { res_key, .. } => {
                Some(res_key.clone())
            }
            CoreStructConnPartnerUnresolved::AtomNotFound { res_key, .. } => Some(res_key.clone()),
            _ => None,
        }
    }

    #[getter]
    pub fn atom_name(&self) -> Option<String> {
        match &self.inner {
            CoreStructConnPartnerUnresolved::AtomNotFound { atom_name, .. } => {
                Some(atom_name.clone())
            }
            _ => None,
        }
    }

    #[getter]
    pub fn message(&self) -> String {
        format!("{}", self.inner)
    }

    pub fn __str__(&self) -> String {
        format!("{}", self.inner)
    }

    pub fn __repr__(&self) -> String {
        format!(
            "StructConnPartnerUnresolved(reason='{}', message='{}')",
            self.reason(),
            self.inner
        )
    }
}

/// Information about a `_struct_conn` record that could not be resolved into a bond.
#[pyclass(
    name = "UnresolvedStructConn",
    module = "proteindf_bridge.rs",
    skip_from_py_object
)]
#[derive(Clone)]
pub struct PyUnresolvedStructConn {
    pub(crate) inner: CoreUnresolvedStructConn,
}

impl From<CoreUnresolvedStructConn> for PyUnresolvedStructConn {
    fn from(inner: CoreUnresolvedStructConn) -> Self {
        Self { inner }
    }
}

#[pymethods]
impl PyUnresolvedStructConn {
    #[getter]
    pub fn model_name(&self) -> String {
        self.inner.model_name.clone()
    }

    #[getter]
    pub fn conn_id(&self) -> String {
        self.inner.conn_id.clone()
    }

    #[getter]
    pub fn conn_type_id(&self) -> String {
        self.inner.conn_type_id.clone()
    }

    #[getter]
    pub fn ptnr1_unresolved(&self) -> Option<PyStructConnPartnerUnresolved> {
        self.inner
            .ptnr1_unresolved
            .clone()
            .map(PyStructConnPartnerUnresolved::from)
    }

    #[getter]
    pub fn ptnr2_unresolved(&self) -> Option<PyStructConnPartnerUnresolved> {
        self.inner
            .ptnr2_unresolved
            .clone()
            .map(PyStructConnPartnerUnresolved::from)
    }

    #[getter]
    pub fn message(&self) -> String {
        self.inner.message.clone()
    }

    pub fn __str__(&self) -> String {
        format!("{}", self.inner)
    }

    pub fn __repr__(&self) -> String {
        format!(
            "UnresolvedStructConn(conn_id='{}', conn_type_id='{}', model_name='{}', message='{}')",
            self.inner.conn_id, self.inner.conn_type_id, self.inner.model_name, self.inner.message
        )
    }
}

/// Result of parsing macromolecular structure from mmCIF, containing the atom hierarchy
/// and any `_struct_conn` records that could not be resolved.
#[pyclass(
    name = "MmcifStructureReport",
    module = "proteindf_bridge.rs",
    skip_from_py_object
)]
pub struct PyMmcifStructureReport {
    atomgroup: Py<PyAtomGroup>,
    unresolved_struct_conns: Vec<PyUnresolvedStructConn>,
    has_unresolved: bool,
}

impl PyMmcifStructureReport {
    pub fn new(py: Python<'_>, report: CoreMmcifStructureReport) -> PyResult<Self> {
        let has_unresolved = report.has_unresolved();
        let py_ag = Py::new(py, PyAtomGroup::from_core(report.atomgroup))?;
        let unresolved_struct_conns = report
            .unresolved_struct_conns
            .into_iter()
            .map(PyUnresolvedStructConn::from)
            .collect();
        Ok(Self {
            atomgroup: py_ag,
            unresolved_struct_conns,
            has_unresolved,
        })
    }
}

#[pymethods]
impl PyMmcifStructureReport {
    #[getter]
    pub fn atomgroup(&self, py: Python<'_>) -> Py<PyAtomGroup> {
        self.atomgroup.clone_ref(py)
    }

    #[getter]
    pub fn unresolved_struct_conns(&self) -> Vec<PyUnresolvedStructConn> {
        self.unresolved_struct_conns.clone()
    }

    pub fn has_unresolved(&self) -> bool {
        self.has_unresolved
    }

    pub fn __repr__(&self, py: Python<'_>) -> String {
        let atom_count = self.atomgroup.bind(py).borrow().get_number_of_all_atoms();
        format!(
            "MmcifStructureReport(atoms={}, unresolved_bonds={})",
            atom_count,
            self.unresolved_struct_conns.len()
        )
    }
}

#[pyclass(
    name = "SimpleMmcif",
    module = "proteindf_bridge.rs",
    skip_from_py_object
)]
#[derive(Clone)]
pub struct PySimpleMmcif {
    pub(crate) inner: CoreSimpleMmcif,
    atomgroup: Option<PyAtomGroup>,
    data_block_name: String,
    charge_to_b_factor: bool,
}

impl Default for PySimpleMmcif {
    fn default() -> Self {
        Self {
            inner: CoreSimpleMmcif::default(),
            atomgroup: None,
            data_block_name: "structure".to_string(),
            charge_to_b_factor: false,
        }
    }
}

#[pymethods]
impl PySimpleMmcif {
    #[new]
    #[pyo3(signature = (file_path=None))]
    pub fn new(file_path: Option<&str>) -> PyResult<Self> {
        let mut inner = CoreSimpleMmcif::new();
        if let Some(path) = file_path {
            inner.load(path).map_err(to_py_err)?;
        }
        Ok(Self {
            inner,
            atomgroup: None,
            data_block_name: "structure".to_string(),
            charge_to_b_factor: false,
        })
    }

    pub fn load(&mut self, file_path: &str) -> PyResult<()> {
        self.inner.load(file_path).map_err(to_py_err)
    }

    pub fn get_molecule_names(&self) -> Vec<String> {
        self.inner.get_molecule_names()
    }

    #[getter]
    pub fn data<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
        let root_dict = pyo3::types::PyDict::new(py);
        for (block_name, block) in self.inner.data() {
            let kv_dict = pyo3::types::PyDict::new(py);
            for (k, v) in &block.key_values {
                kv_dict.set_item(k, v)?;
            }
            let tables_list = pyo3::types::PyList::empty(py);
            for table in &block.tables {
                let table_list = pyo3::types::PyList::empty(py);
                for row in table {
                    let row_dict = pyo3::types::PyDict::new(py);
                    for (k, v) in row {
                        row_dict.set_item(k, v)?;
                    }
                    table_list.append(row_dict)?;
                }
                tables_list.append(table_list)?;
            }
            let block_list =
                pyo3::types::PyList::new(py, &[kv_dict.into_any(), tables_list.into_any()])?;
            root_dict.set_item(block_name, block_list)?;
        }
        Ok(root_dict)
    }

    #[staticmethod]
    pub fn from_data(_py: Python<'_>, data: &Bound<'_, pyo3::PyAny>) -> PyResult<Self> {
        let dict = data.cast::<pyo3::types::PyDict>().map_err(|_| {
            to_py_err(proteindf_bridge::error::BridgeError::input_error(
                "SimpleMmcif.from_data",
                "expected dict for data",
            ))
        })?;

        let mut cif_text = String::new();
        for (block_name_obj, block_data_obj) in dict.iter() {
            let block_name: String = block_name_obj.extract().map_err(|e| {
                to_py_err(proteindf_bridge::error::BridgeError::input_error(
                    "SimpleMmcif.from_data",
                    format!("invalid block name: {e}"),
                ))
            })?;
            if !block_name.starts_with("data_") {
                cif_text.push_str(&format!("data_{block_name}\n"));
            } else {
                cif_text.push_str(&format!("{block_name}\n"));
            }
            let block_list: Bound<'_, pyo3::types::PyList> =
                block_data_obj.extract().map_err(|e| {
                    to_py_err(proteindf_bridge::error::BridgeError::input_error(
                        "SimpleMmcif.from_data",
                        format!("invalid block data, expected list: {e}"),
                    ))
                })?;
            if block_list.len() < 2 {
                return Err(to_py_err(
                    proteindf_bridge::error::BridgeError::input_error(
                        "SimpleMmcif.from_data",
                        "block data must have at least 2 elements [dict, tables]",
                    ),
                ));
            }
            let item0 = block_list.get_item(0)?;
            let kv_dict: Bound<'_, pyo3::types::PyDict> = item0.extract().map_err(|e| {
                to_py_err(proteindf_bridge::error::BridgeError::input_error(
                    "SimpleMmcif.from_data",
                    format!("block[0] must be dict: {e}"),
                ))
            })?;
            for (k_obj, v_obj) in kv_dict.iter() {
                let k: String = k_obj.extract().map_err(|e| {
                    to_py_err(proteindf_bridge::error::BridgeError::input_error(
                        "SimpleMmcif.from_data",
                        format!("invalid key in block[0]: {e}"),
                    ))
                })?;
                let v: String = v_obj.extract().map_err(|e| {
                    to_py_err(proteindf_bridge::error::BridgeError::input_error(
                        "SimpleMmcif.from_data",
                        format!("invalid value in block[0]: {e}"),
                    ))
                })?;
                cif_text.push_str(&format!("{k} \"{v}\"\n"));
            }
            let item1 = block_list.get_item(1)?;
            let tables_list: Bound<'_, pyo3::types::PyList> = item1.extract().map_err(|e| {
                to_py_err(proteindf_bridge::error::BridgeError::input_error(
                    "SimpleMmcif.from_data",
                    format!("block[1] must be list of tables: {e}"),
                ))
            })?;
            for table_obj in tables_list.iter() {
                let table: Bound<'_, pyo3::types::PyList> = table_obj.extract().map_err(|e| {
                    to_py_err(proteindf_bridge::error::BridgeError::input_error(
                        "SimpleMmcif.from_data",
                        format!("table must be list of rows: {e}"),
                    ))
                })?;
                if table.is_empty() {
                    continue;
                }
                cif_text.push_str("loop_\n");
                let first_item = table.get_item(0)?;
                let first_row: Bound<'_, pyo3::types::PyDict> =
                    first_item.extract().map_err(|e| {
                        to_py_err(proteindf_bridge::error::BridgeError::input_error(
                            "SimpleMmcif.from_data",
                            format!("table row must be dict: {e}"),
                        ))
                    })?;
                let mut headers: Vec<String> = Vec::new();
                for (h_obj, _) in first_row.iter() {
                    let h: String = h_obj.extract().map_err(|e| {
                        to_py_err(proteindf_bridge::error::BridgeError::input_error(
                            "SimpleMmcif.from_data",
                            format!("invalid header in table row: {e}"),
                        ))
                    })?;
                    cif_text.push_str(&format!("{h}\n"));
                    headers.push(h);
                }
                for row_obj in table.iter() {
                    let row_dict: Bound<'_, pyo3::types::PyDict> =
                        row_obj.extract().map_err(|e| {
                            to_py_err(proteindf_bridge::error::BridgeError::input_error(
                                "SimpleMmcif.from_data",
                                format!("table row must be dict: {e}"),
                            ))
                        })?;
                    let mut line_vals = Vec::new();
                    for h in &headers {
                        if let Some(val_obj) = row_dict.get_item(h)? {
                            let v: String = val_obj.extract().map_err(|e| {
                                to_py_err(proteindf_bridge::error::BridgeError::input_error(
                                    "SimpleMmcif.from_data",
                                    format!("invalid value in table field: {e}"),
                                ))
                            })?;
                            if v.is_empty() {
                                line_vals.push("?".to_string());
                            } else if v.contains(' ') || v.contains('\t') || v.contains('\n') {
                                line_vals.push(format!("\"{v}\""));
                            } else {
                                line_vals.push(v);
                            }
                        } else {
                            line_vals.push("?".to_string());
                        }
                    }
                    cif_text.push_str(&format!("{}\n", line_vals.join(" ")));
                }
            }
        }
        let inner = CoreSimpleMmcif::from_str(&cif_text).map_err(to_py_err)?;
        Ok(Self {
            inner,
            atomgroup: None,
            data_block_name: "structure".to_string(),
            charge_to_b_factor: false,
        })
    }

    #[pyo3(signature = (name=None))]
    pub fn get_atomgroup(&self, name: Option<&str>) -> PyResult<PyAtomGroup> {
        let ag = match name {
            Some(n) => self.inner.get_atomgroup(n).map_err(to_py_err)?,
            None => {
                let names = self.inner.get_molecule_names();
                let first = names.first().ok_or_else(|| {
                    to_py_err(BridgeError::input_error(
                        "mmCIF",
                        "No data blocks found in file",
                    ))
                })?;
                self.inner.get_atomgroup(first).map_err(to_py_err)?
            }
        };
        Ok(PyAtomGroup::from_core(ag))
    }

    #[pyo3(signature = (select_model=None, select_altloc=Some("A")))]
    pub fn get_structure_atomgroup(
        &self,
        select_model: Option<usize>,
        select_altloc: Option<&str>,
    ) -> PyResult<PyAtomGroup> {
        let ag = self
            .inner
            .get_structure_atomgroup(select_model, select_altloc)
            .map_err(to_py_err)?;
        Ok(PyAtomGroup::from_core(ag))
    }

    #[pyo3(signature = (block_name, select_model=None, select_altloc=Some("A")))]
    pub fn get_structure_atomgroup_for_block(
        &self,
        block_name: &str,
        select_model: Option<usize>,
        select_altloc: Option<&str>,
    ) -> PyResult<PyAtomGroup> {
        let ag = self
            .inner
            .get_structure_atomgroup_for_block(block_name, select_model, select_altloc)
            .map_err(to_py_err)?;
        Ok(PyAtomGroup::from_core(ag))
    }

    #[pyo3(signature = (select_model=None, select_altloc=Some("A")))]
    pub fn get_structure_atomgroup_with_report(
        &self,
        py: Python<'_>,
        select_model: Option<usize>,
        select_altloc: Option<&str>,
    ) -> PyResult<PyMmcifStructureReport> {
        let report = self
            .inner
            .get_structure_atomgroup_with_report(select_model, select_altloc)
            .map_err(to_py_err)?;
        PyMmcifStructureReport::new(py, report)
    }

    #[pyo3(signature = (block_name, select_model=None, select_altloc=Some("A")))]
    pub fn get_structure_atomgroup_for_block_with_report(
        &self,
        py: Python<'_>,
        block_name: &str,
        select_model: Option<usize>,
        select_altloc: Option<&str>,
    ) -> PyResult<PyMmcifStructureReport> {
        let report = self
            .inner
            .get_structure_atomgroup_for_block_with_report(block_name, select_model, select_altloc)
            .map_err(to_py_err)?;
        PyMmcifStructureReport::new(py, report)
    }

    #[pyo3(signature = (
        atomgroup,
        data_block_name = "structure",
        charge_to_b_factor = false
    ))]
    pub fn set_by_atomgroup(
        &mut self,
        atomgroup: &PyAtomGroup,
        data_block_name: &str,
        charge_to_b_factor: bool,
    ) -> PyResult<()> {
        let opts = MmcifWriteOptions {
            data_block_name: data_block_name.to_string(),
            charge_to_b_factor,
        };
        validate_for_mmcif_write(&atomgroup.inner, &opts).map_err(to_py_err)?;
        self.atomgroup = Some(atomgroup.clone());
        self.data_block_name = data_block_name.to_string();
        self.charge_to_b_factor = charge_to_b_factor;
        Ok(())
    }

    #[pyo3(signature = (
        file_path,
        atomgroup = None,
        data_block_name = None,
        charge_to_b_factor = None
    ))]
    pub fn save(
        &self,
        file_path: &str,
        atomgroup: Option<&PyAtomGroup>,
        data_block_name: Option<&str>,
        charge_to_b_factor: Option<bool>,
    ) -> PyResult<()> {
        let ag = match (atomgroup, &self.atomgroup) {
            (Some(ag), _) => &ag.inner,
            (None, Some(ag)) => &ag.inner,
            (None, None) => {
                return Err(to_py_err(BridgeError::input_error(
                    "mmCIF",
                    "No AtomGroup specified for save(); call set_by_atomgroup() first or pass atomgroup",
                )));
            }
        };
        let c2b = charge_to_b_factor.unwrap_or(self.charge_to_b_factor);
        let block_name = data_block_name
            .map(|s| s.to_string())
            .unwrap_or_else(|| self.data_block_name.clone());
        let opts = MmcifWriteOptions {
            data_block_name: block_name,
            charge_to_b_factor: c2b,
        };
        CoreSimpleMmcif::save_structure(ag, file_path, &opts).map_err(to_py_err)
    }

    #[staticmethod]
    #[pyo3(signature = (
        atomgroup,
        file_path,
        data_block_name = "structure",
        charge_to_b_factor = false
    ))]
    pub fn save_structure(
        atomgroup: &PyAtomGroup,
        file_path: &str,
        data_block_name: &str,
        charge_to_b_factor: bool,
    ) -> PyResult<()> {
        let opts = MmcifWriteOptions {
            data_block_name: data_block_name.to_string(),
            charge_to_b_factor,
        };
        CoreSimpleMmcif::save_structure(&atomgroup.inner, file_path, &opts).map_err(to_py_err)
    }

    #[staticmethod]
    #[pyo3(signature = (
        atomgroup,
        data_block_name = "structure",
        charge_to_b_factor = false
    ))]
    pub fn write_structure(
        atomgroup: &PyAtomGroup,
        data_block_name: &str,
        charge_to_b_factor: bool,
    ) -> PyResult<String> {
        let opts = MmcifWriteOptions {
            data_block_name: data_block_name.to_string(),
            charge_to_b_factor,
        };
        let mut buf = Vec::new();
        CoreSimpleMmcif::write_structure(&atomgroup.inner, &mut buf, &opts).map_err(to_py_err)?;
        String::from_utf8(buf)
            .map_err(|e| to_py_err(BridgeError::input_error("write", e.to_string())))
    }

    pub fn get_text(&self) -> PyResult<String> {
        if let Some(ref ag) = self.atomgroup {
            let opts = MmcifWriteOptions {
                data_block_name: self.data_block_name.clone(),
                charge_to_b_factor: self.charge_to_b_factor,
            };
            let mut buf = Vec::new();
            CoreSimpleMmcif::write_structure(&ag.inner, &mut buf, &opts).map_err(to_py_err)?;
            String::from_utf8(buf)
                .map_err(|e| to_py_err(BridgeError::input_error("write", e.to_string())))
        } else {
            Err(to_py_err(BridgeError::input_error(
                "mmCIF",
                "No AtomGroup specified; call set_by_atomgroup() first",
            )))
        }
    }

    pub fn __str__(&self) -> PyResult<String> {
        self.get_text()
    }
}
