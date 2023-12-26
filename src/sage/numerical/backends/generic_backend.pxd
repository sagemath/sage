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
    cpdef set_variable_type(self, int variable, int vtype) noexcept
    cpdef set_sense(self, int sense) noexcept
    cpdef objective_coefficient(self, int variable, coeff=*) noexcept
    cpdef objective_constant_term(self, d=*) noexcept
    cpdef set_objective(self, list coeff, d=*) noexcept
    cpdef set_verbosity(self, int level) noexcept
    cpdef add_linear_constraint(self, coefficients, lower_bound, upper_bound, name=*) noexcept
    cpdef add_linear_constraint_vector(self, degree, coefficients, lower_bound, upper_bound, name=*) noexcept
    cpdef remove_constraint(self, int) noexcept
    cpdef remove_constraints(self, constraints) noexcept
    cpdef add_col(self, indices, coeffs) noexcept
    cpdef add_linear_constraints(self, int number, lower_bound, upper_bound, names=*) noexcept
    cpdef int solve(self) except -1
    cpdef get_objective_value(self) noexcept
    cpdef best_known_objective_bound(self) noexcept
    cpdef get_relative_objective_gap(self) noexcept
    cpdef get_variable_value(self, int variable) noexcept
    cpdef bint is_maximization(self) noexcept
    cpdef write_lp(self, name) noexcept
    cpdef write_mps(self, name, int modern) noexcept
    cpdef row(self, int i) noexcept
    cpdef int ncols(self) noexcept
    cpdef int nrows(self) noexcept
    cpdef bint is_variable_binary(self, int) noexcept
    cpdef bint is_variable_integer(self, int) noexcept
    cpdef bint is_variable_continuous(self, int) noexcept
    cpdef problem_name(self, name = *) noexcept
    cpdef row_bounds(self, int index) noexcept
    cpdef col_bounds(self, int index) noexcept
    cpdef row_name(self, int index) noexcept
    cpdef col_name(self, int index) noexcept
    cpdef variable_upper_bound(self, int index, value = *) noexcept
    cpdef variable_lower_bound(self, int index, value = *) noexcept
    cpdef solver_parameter(self, name, value=*) noexcept
    cpdef zero(self) noexcept
    cpdef base_ring(self) noexcept
    cpdef __copy__(self) noexcept
    cpdef copy(self) noexcept
    cpdef bint is_variable_basic(self, int index) noexcept
    cpdef bint is_variable_nonbasic_at_lower_bound(self, int index) noexcept
    cpdef bint is_slack_variable_basic(self, int index) noexcept
    cpdef bint is_slack_variable_nonbasic_at_lower_bound(self, int index) noexcept

    cdef object obj_constant_term

cpdef GenericBackend get_solver(constraint_generation = ?, solver = ?, base_ring = ?) noexcept
