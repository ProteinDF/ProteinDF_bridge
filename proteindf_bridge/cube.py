#!/usr/bin/env python
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: The ProteinDF development team
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Gaussian CUBE format reader wrapping Rust backend.
"""

from .rs import Cube

__all__ = ["Cube"]
