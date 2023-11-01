#*****************************************************************************
#       Copyright (C) 2014 Ingolfur Edvardsson <ingolfured@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
cdef class GenericSDPBackend:
    cpdef int add_variable(self, obj=*, name=*) except -1
    cpdef int add_variables(self, int, names=*) except -1
    cpdef set_sense(self, int sense) noexcept
    cpdef objective_coefficient(self, int variable, coeff=*) noexcept
    cpdef set_objective(self, list coeff, d=*) noexcept
    cpdef add_linear_constraint(self, constraints, name=*) noexcept
    cpdef add_linear_constraints(self, int number, names=*) noexcept
    cpdef int solve(self) except -1
    cpdef get_objective_value(self) noexcept
    cpdef get_variable_value(self, int variable) noexcept
    cpdef dual_variable(self, int variable, sparse=*) noexcept
    cpdef slack(self, int variable, sparse=*) noexcept
    cpdef bint is_maximization(self) noexcept
    cpdef row(self, int i) noexcept
    cpdef int ncols(self) noexcept
    cpdef int nrows(self) noexcept
    cpdef problem_name(self, name=*) noexcept
    cpdef row_name(self, int index) noexcept
    cpdef col_name(self, int index) noexcept
    cpdef solver_parameter(self, name, value=*) noexcept
    cpdef zero(self) noexcept
    cpdef base_ring(self) noexcept

    cdef obj_constant_term
    cdef dict matrices_dim

cpdef GenericSDPBackend get_solver(solver=?, base_ring=?) noexcept
