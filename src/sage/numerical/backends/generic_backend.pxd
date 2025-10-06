##############################################################################
#       Copyright (C) 2010 Nathann Cohen <nathann.cohen@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

from sage.structure.sage_object cimport SageObject

# We inherit from SageObject to make some testing infrastructure available.

cdef class GenericBackend (SageObject):
    cpdef int add_variable(self, lower_bound=*, upper_bound=*, binary=*, continuous=*, integer=*, obj=*, name=*) except -1
    cpdef int add_variables(self, int, lower_bound=*, upper_bound=*, binary=*, continuous=*, integer=*, obj=*, names=*) except -1
    cpdef set_variable_type(self, int variable, int vtype)
    cpdef set_sense(self, int sense)
    cpdef objective_coefficient(self, int variable, coeff=*)
    cpdef objective_constant_term(self, d=*)
    cpdef set_objective(self, list coeff, d=*)
    cpdef set_verbosity(self, int level)
    cpdef add_linear_constraint(self, coefficients, lower_bound, upper_bound, name=*)
    cpdef add_linear_constraint_vector(self, degree, coefficients, lower_bound, upper_bound, name=*)
    cpdef remove_constraint(self, int)
    cpdef remove_constraints(self, constraints)
    cpdef add_col(self, indices, coeffs)
    cpdef add_linear_constraints(self, int number, lower_bound, upper_bound, names=*)
    cpdef int solve(self) except -1
    cpdef get_objective_value(self)
    cpdef best_known_objective_bound(self)
    cpdef get_relative_objective_gap(self)
    cpdef get_variable_value(self, int variable)
    cpdef bint is_maximization(self) noexcept
    cpdef write_lp(self, name)
    cpdef write_mps(self, name, int modern)
    cpdef row(self, int i)
    cpdef int ncols(self) noexcept
    cpdef int nrows(self) noexcept
    cpdef bint is_variable_binary(self, int) noexcept
    cpdef bint is_variable_integer(self, int) noexcept
    cpdef bint is_variable_continuous(self, int) noexcept
    cpdef problem_name(self, name=*)
    cpdef row_bounds(self, int index)
    cpdef col_bounds(self, int index)
    cpdef row_name(self, int index)
    cpdef col_name(self, int index)
    cpdef variable_upper_bound(self, int index, value=*)
    cpdef variable_lower_bound(self, int index, value=*)
    cpdef solver_parameter(self, name, value=*)
    cpdef zero(self)
    cpdef base_ring(self)
    cpdef __copy__(self)
    cpdef copy(self)
    cpdef bint is_variable_basic(self, int index) noexcept
    cpdef bint is_variable_nonbasic_at_lower_bound(self, int index) noexcept
    cpdef bint is_slack_variable_basic(self, int index) noexcept
    cpdef bint is_slack_variable_nonbasic_at_lower_bound(self, int index) noexcept

    cdef object obj_constant_term

cpdef GenericBackend get_solver(constraint_generation=?, solver=?, base_ring=?)
