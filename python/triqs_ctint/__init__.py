# Copyright (c) 2017--present, The Simons Foundation
# This file is part of TRIQS/ctint and is licensed under the terms of GPLv3 or later.
# SPDX-License-Identifier: GPL-3.0-or-later
# See LICENSE in the root of this distribution for details.


import warnings
#warnings.simplefilter(action='ignore', category=FutureWarning)
#
#r"""
#DOC
#
#"""
from .solver import Solver
from .solver_core import SolverCore

__all__ = ['Solver','SolverCore']
