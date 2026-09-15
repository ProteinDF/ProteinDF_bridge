// Copyright (C) 2014 The ProteinDF development team.
// see also AUTHORS and README if provided.
//
// This file is a part of the ProteinDF software package.
//
// The ProteinDF is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// The ProteinDF is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

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

    /// PeriodicTable error: atomic weight not found.
    #[error("PeriodicTable.atomic_weight(): no atomic weight for atom {0}")]
    AtomicWeightNotFound(usize),
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
