r"""
Families of Mixed Integer Linear Programs

Naming of methods follow :mod:`sage.numerical.interactive_simplex_method`
and https://github.com/mkoeppe/sage-numerical-interactive-mip/blob/master/sage_numerical_interactive_mip/interactive_milp_problem.py
where possible.

"""

# ****************************************************************************
#       Copyright (C) 2024 Matthias Koeppe
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.sets.family import LazyFamily


class LPProblem__max_cx__Ax_le_b__Family_b(LazyFamily):

    def __init__(self, A, b_set, c, *, objective_constant_term=0):
        self._A = A
        self._b_set = b_set
        self._c = c
        self._objective_constant_term = objective_constant_term

    # Aliases for the standard notation
    A = constraint_coefficients
    b_set = constant_terms_set
    c = objective_coefficients
    m = n_constraints
    n = n_variables
