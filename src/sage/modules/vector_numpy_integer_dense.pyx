# sage.doctest: optional - numpy
r"""
Dense integer vectors using a NumPy backend.

EXAMPLES::

    sage: from sage.modules.vector_numpy_integer_dense import Vector_numpy_integer_dense
    sage: v = Vector_numpy_integer_dense(FreeModule(ZZ, 3), [0, 0, 0]); v
    (0, 0, 0)
    sage: v[1] = 42
    sage: v
    (0, 42, 0)
    sage: v.numpy()
    array([ 0, 42,  0])               # 64-bit
    array([ 0, 42,  0], dtype=int64)  # 32-bit
"""

# ****************************************************************************
#       Copyright (C) 2021      Matthias Koeppe
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

cimport numpy
import numpy

from sage.rings.integer_ring import ZZ

# This is for the NumPy C API (the PyArray... functions) to work
numpy.import_array()


cdef class Vector_numpy_integer_dense(Vector_numpy_dense):

    def __cinit__(self, parent, entries, coerce=True, copy=True):
        self._numpy_dtype = numpy.dtype('int64')
        self._numpy_dtypeint = numpy.NPY_INT64
        self._python_dtype = int
        self._sage_dtype = ZZ
        self.__create_vector__()
        return
