# sage_setup: distribution = sagemath-polyhedra
#*****************************************************************************
#       Copyright (C) 2017 Matthias Koeppe <mkoeppe@math.ucdavis.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.numerical.backends.generic_backend cimport GenericBackend

cdef class SCIPBackend(GenericBackend):

    cdef model
    cdef object variables
    cdef object constraints

    cpdef _get_model(self)
    cpdef get_row_prim(self, int i)
    cpdef write_cip(self, filename)
