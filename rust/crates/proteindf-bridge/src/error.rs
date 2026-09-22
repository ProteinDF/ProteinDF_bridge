// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use thiserror::Error;

/// Base error type for the Bridge module, corresponding to `BrError` and its subclasses.
#[derive(Error, Debug, Clone, PartialEq)]
pub enum BridgeError {
    /// General Bridge module error (`BrError`).
    #[error("Bridge module error: {0}")]
    General(String),

    /// Error raised for invalid input (`BrInputError`).
    #[error("Input Error: {msg} ({expr})")]
    InputError { expr: String, msg: String },

    /// Error raised for invalid value (`BrValueError`).
    #[error("Value Error: {msg} ({expr})")]
    ValueError { expr: String, msg: String },

    /// PeriodicTable error: atomic number not found.
    #[error("PeriodicTable.get_symbol(): atomic number {0} not found.")]
    AtomicNumberNotFound(usize),

    /// PeriodicTable error: element symbol not found.
    #[error("PeriodicTable.get_atomic_number(): symbol '{0}' not found.")]
    SymbolNotFound(String),

    /// PeriodicTable error: VDW radius not found.
    #[error("PeriodicTable.vdw(): no VDW radius for atom {0}")]
    VdwRadiusNotFound(usize),

    /// PeriodicTable error: covalent radius not found.
    #[error("PeriodicTable.covalent_radius(): no covalent radius for atom {0}")]
    CovalentRadiusNotFound(usize),

    /// PeriodicTable error: atomic weight not found.
    #[error("PeriodicTable.atomic_weight(): no atomic weight for atom {0}")]
    AtomicWeightNotFound(usize),

    /// I/O error during file operations.
    #[error("I/O error: {0}")]
    Io(String),

    /// MessagePack serialization or deserialization error.
    #[error("MessagePack error: {0}")]
    MsgPack(String),

    /// YAML serialization or deserialization error.
    #[error("YAML error: {0}")]
    Yaml(String),

    /// Zstd compression or decompression error.
    #[error("Zstd error: {0}")]
    Zstd(String),
}

impl From<std::io::Error> for BridgeError {
    fn from(err: std::io::Error) -> Self {
        BridgeError::Io(err.to_string())
    }
}

impl BridgeError {
    /// Helper to create a general `BrError`.
    pub fn general(msg: impl Into<String>) -> Self {
        Self::General(msg.into())
    }

    /// Helper to create a `BrInputError`.
    pub fn input_error(expr: impl Into<String>, msg: impl Into<String>) -> Self {
        Self::InputError {
            expr: expr.into(),
            msg: msg.into(),
        }
    }

    /// Helper to create a `BrValueError`.
    pub fn value_error(expr: impl Into<String>, msg: impl Into<String>) -> Self {
        Self::ValueError {
            expr: expr.into(),
            msg: msg.into(),
        }
    }
}

pub type Result<T> = std::result::Result<T, BridgeError>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_display() {
        let err = BridgeError::general("test message");
        assert_eq!(err.to_string(), "Bridge module error: test message");

        let err = BridgeError::input_error("foo", "invalid input");
        assert_eq!(err.to_string(), "Input Error: invalid input (foo)");

        let err = BridgeError::value_error("bar", "invalid value");
        assert_eq!(err.to_string(), "Value Error: invalid value (bar)");
    }
}
