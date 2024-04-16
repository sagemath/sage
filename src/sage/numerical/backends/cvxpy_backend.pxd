# sage_setup: distribution = sagemath-polyhedra
##############################################################################
#       Copyright (C) 2010 Nathann Cohen <nathann.cohen@gmail.com>
#       Copyright (C) 2022 Matthias Koeppe <mkoeppe@math.ucdavis.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

from sage.numerical.backends.generic_backend cimport GenericBackend

cdef class CVXPYBackend(GenericBackend):

    cdef object variables
    cdef object problem
    cdef object prob_name
    cdef object constraint_names

    cdef object _cvxpy_solver
    cdef object _cvxpy_solver_args

    cdef list objective_coefficients
    cdef list Matrix

    cdef list row_lower_bound
    cdef list row_upper_bound
    cdef list col_lower_bound
    cdef list col_upper_bound

    cpdef int add_variable(self,
                           lower_bound=*,
                           upper_bound=*,
                           binary=*,
                           continuous=*,
                           integer=*,
                           obj=*,
                           name=*,
                           coefficients=*) \
                           except -1

    cpdef cvxpy_problem(self)
