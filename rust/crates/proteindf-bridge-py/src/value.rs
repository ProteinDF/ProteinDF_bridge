// SPDX-FileCopyrightText: The ProteinDF development team
// SPDX-License-Identifier: GPL-3.0-or-later

use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
use pyo3::types::{PyBool, PyBytes, PyDict, PyList, PyTuple};
use rmpv::{Integer, Utf8String, Value};

/// Converts a Python object into a MessagePack `rmpv::Value`.
pub fn py_to_value(_py: Python<'_>, obj: &Bound<'_, PyAny>) -> PyResult<Value> {
    if obj.is_none() {
        return Ok(Value::Nil);
    }

    // Check bool before int, as PyBool inherits from PyLong (int) in Python.
    if let Ok(b) = obj.cast_exact::<PyBool>() {
        return Ok(Value::Boolean(b.is_true()));
    }

    if let Ok(i) = obj.extract::<i64>() {
        return Ok(Value::Integer(Integer::from(i)));
    }
    if let Ok(u) = obj.extract::<u64>() {
        return Ok(Value::Integer(Integer::from(u)));
    }

    if let Ok(f) = obj.extract::<f64>() {
        return Ok(Value::F64(f));
    }

    if let Ok(s) = obj.extract::<String>() {
        return Ok(Value::String(Utf8String::from(s)));
    }

    if let Ok(bytes) = obj.cast_exact::<PyBytes>() {
        return Ok(Value::Binary(bytes.as_bytes().to_vec()));
    }

    if let Ok(dict) = obj.cast::<PyDict>() {
        let mut map = Vec::with_capacity(dict.len());
        for (k, v) in dict.iter() {
            let k_val = py_to_value(_py, &k)?;
            let v_val = py_to_value(_py, &v)?;
            map.push((k_val, v_val));
        }
        return Ok(Value::Map(map));
    }

    if let Ok(list) = obj.cast::<PyList>() {
        let mut arr = Vec::with_capacity(list.len());
        for item in list.iter() {
            arr.push(py_to_value(_py, &item)?);
        }
        return Ok(Value::Array(arr));
    }

    if let Ok(tuple) = obj.cast::<PyTuple>() {
        let mut arr = Vec::with_capacity(tuple.len());
        for item in tuple.iter() {
            arr.push(py_to_value(_py, &item)?);
        }
        return Ok(Value::Array(arr));
    }

    Err(PyTypeError::new_err(format!(
        "Cannot convert Python type '{}' to MessagePack Value",
        obj.get_type()
    )))
}

/// Converts a MessagePack `rmpv::Value` into a Python object.
pub fn value_to_py<'py>(py: Python<'py>, val: &Value) -> PyResult<Bound<'py, PyAny>> {
    match val {
        Value::Nil => Ok(py.None().into_bound(py)),
        Value::Boolean(b) => Ok((*b).into_pyobject(py)?.to_owned().into_any()),
        Value::Integer(i) => {
            if let Some(v) = i.as_i64() {
                Ok(v.into_pyobject(py)?.into_any())
            } else if let Some(v) = i.as_u64() {
                Ok(v.into_pyobject(py)?.into_any())
            } else {
                Err(PyTypeError::new_err("Integer value out of range"))
            }
        }
        Value::F32(f) => Ok((*f as f64).into_pyobject(py)?.into_any()),
        Value::F64(f) => Ok((*f).into_pyobject(py)?.into_any()),
        Value::String(s) => {
            if let Some(utf8) = s.as_str() {
                Ok(utf8.into_pyobject(py)?.to_owned().into_any())
            } else {
                Ok(PyBytes::new(py, s.as_bytes()).into_any())
            }
        }
        Value::Binary(bytes) => Ok(PyBytes::new(py, bytes).into_any()),
        Value::Array(arr) => {
            let py_list = PyList::empty(py);
            for item in arr {
                py_list.append(value_to_py(py, item)?)?;
            }
            Ok(py_list.into_any())
        }
        Value::Map(map) => {
            let py_dict = PyDict::new(py);
            for (k, v) in map {
                let py_k = match k {
                    Value::String(s) => {
                        if let Some(utf8) = s.as_str() {
                            utf8.into_pyobject(py)?.to_owned().into_any()
                        } else {
                            PyBytes::new(py, s.as_bytes()).into_any()
                        }
                    }
                    Value::Binary(bytes) => {
                        if let Ok(utf8) = std::str::from_utf8(bytes) {
                            utf8.into_pyobject(py)?.to_owned().into_any()
                        } else {
                            PyBytes::new(py, bytes).into_any()
                        }
                    }
                    _ => value_to_py(py, k)?,
                };
                let py_v = value_to_py(py, v)?;
                py_dict.set_item(py_k, py_v)?;
            }
            Ok(py_dict.into_any())
        }
        Value::Ext(_, bytes) => Ok(PyBytes::new(py, bytes).into_any()),
    }
}
