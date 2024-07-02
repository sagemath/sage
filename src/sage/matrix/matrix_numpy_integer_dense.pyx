# sage.doctest: optional - numpy
r"""
Dense integer matrices using a NumPy backend

EXAMPLES::

    sage: from sage.matrix.matrix_numpy_integer_dense import Matrix_numpy_integer_dense
    sage: M = Matrix_numpy_integer_dense(MatrixSpace(ZZ, 2, 3)); M
    [0 0 0]
    [0 0 0]
    sage: M[1,2] = 47
    sage: M
    [ 0  0  0]
    [ 0  0 47]
    sage: M.numpy()
    array([[ 0,  0,  0], [ 0,  0, 47]])               # 64-bit
    array([[ 0,  0,  0], [ 0,  0, 47]], dtype=int64)  # 32-bit
"""

# ****************************************************************************
#       Copyright (C) 2021      Matthias Koeppe
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.integer_ring import ZZ

cimport numpy as cnumpy

numpy = None
scipy = None

# This is for the Numpy C API to work
cnumpy.import_array()


cdef class Matrix_numpy_integer_dense(Matrix_numpy_dense):
    r"""
    TESTS::

        sage: from sage.matrix.matrix_numpy_integer_dense import Matrix_numpy_integer_dense
        sage: M = Matrix_numpy_integer_dense(MatrixSpace(ZZ, 2, 3))
        sage: TestSuite(M).run()
    """

    def __cinit__(self):
        global numpy
        if numpy is None:
            import numpy
        self._numpy_dtype = numpy.dtype('int64')
        self._numpy_dtypeint = cnumpy.NPY_INT64
        self._python_dtype = int
        self._sage_dtype = ZZ
        self.__create_matrix__()
        return
