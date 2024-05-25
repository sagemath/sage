# sage_setup: distribution = sagemath-glpk
#*****************************************************************************
#       Copyright (C) 2016 Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.numerical.backends.glpk_backend cimport GLPKBackend

cdef class GLPKExactBackend(GLPKBackend):
    cpdef int add_variable(self, lower_bound=*, upper_bound=*, binary=*, continuous=*, integer=*, obj=*, name=*) except -1
    cpdef int add_variables(self, int, lower_bound=*, upper_bound=*, binary=*, continuous=*, integer=*, obj=*, names=*) except -1
    cpdef set_variable_type(self, int variable, int vtype)
