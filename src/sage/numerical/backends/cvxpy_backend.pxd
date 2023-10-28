##############################################################################
#       Copyright (C) 2010 Nathann Cohen <nathann.cohen@gmail.com>
#                     2022 Matthias Koeppe <mkoeppe@math.ucdavis.edu>
#                     2023 Zhongling Xu
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.numerical.backends.matrix_backend cimport MatrixBackend

cdef class CVXPYBackend(MatrixBackend):

    cdef object variables
    cdef object problem

    cdef object _cvxpy_solver
    cdef object _cvxpy_solver_args

    cpdef cvxpy_problem(self)
