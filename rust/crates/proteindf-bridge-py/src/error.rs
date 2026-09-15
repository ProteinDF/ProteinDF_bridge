// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use proteindf_bridge::error::BridgeError;
use pyo3::create_exception;
use pyo3::exceptions::PyException;
use pyo3::prelude::*;

create_exception!(proteindf_bridge_rs, BrError, PyException);
create_exception!(proteindf_bridge_rs, BrInputError, BrError);
create_exception!(proteindf_bridge_rs, BrValueError, BrError);

pub fn to_py_err(err: BridgeError) -> PyErr {
    match err {
        BridgeError::General(msg) => BrError::new_err(msg),
        BridgeError::InputError { expr, msg } => {
            BrInputError::new_err(format!("Input Error: {} ({})", msg, expr))
        }
        BridgeError::ValueError { expr, msg } => {
            BrValueError::new_err(format!("Value Error: {} ({})", msg, expr))
        }
        BridgeError::AtomicNumberNotFound(n) => BrValueError::new_err(format!(
            "PeriodicTable.get_symbol(): atomic number {} not found.",
            n
        )),
        BridgeError::SymbolNotFound(sym) => BrValueError::new_err(format!(
            "PeriodicTable.get_atomic_number(): symbol '{}' not found.",
            sym
        )),
        BridgeError::VdwRadiusNotFound(n) => {
            BrValueError::new_err(format!("PeriodicTable.vdw(): no VDW radius for atom {}", n))
        }
        BridgeError::AtomicWeightNotFound(n) => BrValueError::new_err(format!(
            "PeriodicTable.atomic_weight(): no atomic weight for atom {}",
            n
        )),
        BridgeError::Io(msg) => BrError::new_err(format!("I/O error: {}", msg)),
        BridgeError::MsgPack(msg) => BrValueError::new_err(format!("MessagePack error: {}", msg)),
        BridgeError::Yaml(msg) => BrValueError::new_err(format!("YAML error: {}", msg)),
        BridgeError::Zstd(msg) => BrValueError::new_err(format!("Zstd error: {}", msg)),
    }
}
