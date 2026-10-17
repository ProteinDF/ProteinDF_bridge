// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use proteindf_bridge::error::BridgeError;

use pyo3::create_exception;
use pyo3::exceptions::PyException;
use pyo3::prelude::*;

use std::ffi::CString;

create_exception!(
    proteindf_bridge.rs,
    BrError,
    PyException,
    "Base exception class for ProteinDF bridge errors."
);
create_exception!(
    proteindf_bridge.rs,
    BrInputError,
    BrError,
    "Exception raised for invalid input arguments or formats."
);
create_exception!(
    proteindf_bridge.rs,
    BrValueError,
    BrError,
    "Exception raised for invalid values or unsupported data."
);

pub fn register_exceptions(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = m.py();
    m.add("BrError", py.get_type::<BrError>())?;
    m.add("BrInputError", py.get_type::<BrInputError>())?;
    m.add("BrValueError", py.get_type::<BrValueError>())?;

    let py_code = CString::new(
        r#"
def _init_br_error(self, errmsg=""):
    self.errmsg = errmsg

def _init_br_input_error(self, expr="", msg=None):
    if msg is not None:
        self.errmsg = f"Input Error: {msg} ({expr})"
    else:
        self.errmsg = str(expr)

def _init_br_value_error(self, expr="", msg=None):
    if msg is not None:
        self.errmsg = f"Value Error: {msg} ({expr})"
    else:
        self.errmsg = str(expr)


def _str_br_error(self):
    return f"Bridge module error: {getattr(self, 'errmsg', '')}"

BrError.__init__ = _init_br_error
BrError.__str__ = _str_br_error
BrInputError.__init__ = _init_br_input_error
BrValueError.__init__ = _init_br_value_error
"#,
    )
    .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(e.to_string()))?;

    py.run(py_code.as_c_str(), Some(&m.dict()), None)?;

    let dict = m.dict();
    let _ = dict.del_item("_init_br_error");
    let _ = dict.del_item("_init_br_input_error");
    let _ = dict.del_item("_init_br_value_error");
    let _ = dict.del_item("_str_br_error");

    Ok(())
}

pub fn to_py_err(err: BridgeError) -> PyErr {
    match err {
        BridgeError::General(msg) => BrError::new_err(msg),
        BridgeError::InputError { expr, msg } => BrInputError::new_err((expr, msg)),
        BridgeError::ValueError { expr, msg } => BrValueError::new_err((expr, msg)),
        BridgeError::AtomicNumberNotFound(n) => BrValueError::new_err((format!(
            "PeriodicTable.get_symbol(): atomic number {} not found.",
            n
        ),)),
        BridgeError::SymbolNotFound(sym) => BrValueError::new_err((format!(
            "PeriodicTable.get_atomic_number(): symbol '{}' not found.",
            sym
        ),)),
        BridgeError::VdwRadiusNotFound(n) => {
            BrValueError::new_err((format!("PeriodicTable.vdw(): no VDW radius for atom {}", n),))
        }
        BridgeError::CovalentRadiusNotFound(n) => BrValueError::new_err((format!(
            "PeriodicTable.covalent_radius(): no covalent radius for atom {}",
            n
        ),)),
        BridgeError::AtomicWeightNotFound(n) => BrValueError::new_err((format!(
            "PeriodicTable.atomic_weight(): no atomic weight for atom {}",
            n
        ),)),
        BridgeError::Io(msg) => BrError::new_err(format!("I/O error: {}", msg)),
        BridgeError::MsgPack(msg) => {
            BrValueError::new_err((format!("MessagePack error: {}", msg),))
        }
        BridgeError::Yaml(msg) => BrValueError::new_err((format!("YAML error: {}", msg),)),
        BridgeError::Zstd(msg) => BrValueError::new_err((format!("Zstd error: {}", msg),)),
    }
}
