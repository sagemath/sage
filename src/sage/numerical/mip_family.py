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
        super().__init__(b_set, )

    # Aliases for the standard notation
    A = constraint_coefficients
    b_set = constant_terms_set
    c = objective_coefficients
    m = n_constraints
    n = n_variables


    def restrict(self, smaller_b_set):
        r"""
        Return an :class:`LPProblem__max_cx__Ax_le_b__Family_b` for ..."""


    def optimal_solution_function(self):
        r"""

        """
        # Base case: value function is affine-linear in b


    def value_function(self):
        r"""
        Return a map from b_set to extended reals

        depend on ambient space of b_set

        RR^n, RR^['x1','x2','x3'] etc.
        """
        # Base case: value function is affine-linear in b

        # piecewise:
        #  - construct subfamilies, one for each piece
        #  - pick value_function of subfamily




    def activity_set(self, basis):
        r"""
        Compute (intersection of b_set with) polyhedron of b where basis is an optimal basis (if basis is dual feasible) = parameter set of primal feasibility
        """

    def wolsey(self, ....)
        r"""

        """

    def interactive_lp(self, b):
        r"""
        Specialize to b
        """
        return InteractiveLPProblem(....)


    def mixed_integer_linear_program(self, b, *, solver=None):
        return MixedIntegerLinearProgram()......
