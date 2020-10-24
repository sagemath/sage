# distutils: language = c++
# distutils: libraries = CBLAS_LIBRARIES
# distutils: library_dirs = CBLAS_LIBDIR
# distutils: include_dirs = CBLAS_INCDIR
"""
Dense matrices over `\ZZ/n\ZZ` for `n < 2^{11}` using LinBox's ``Modular<float>``

AUTHORS:
- Burcin Erocal
- Martin Albrecht
"""
###############################################################################
#       Copyright (C) 2011 Burcin Erocal <burcin@erocal.org>
#       Copyright (C) 2011 Martin Albrecht <martinralbrecht@googlemail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.rings.finite_rings.stdint cimport *

from sage.libs.linbox.givaro cimport \
    Modular_float as ModField, \
    Poly1Dom, Dense

from sage.libs.linbox.linbox cimport \
    DenseMatrix_Modular_float as DenseMatrix, \
    reducedRowEchelonize

from sage.libs.linbox.fflas cimport \
    fgemm, pfgemm, fgemv, Det, pDet, Rank, pRank, ReducedRowEchelonForm, pReducedRowEchelonForm, applyP, \
    MinPoly, CharPoly, MinPoly, \
    ModFloatDensePolynomial as ModDensePoly

ctypedef Poly1Dom[ModField, Dense] ModDensePolyRing

# LinBox supports up to 2^11 using float but that's double dog slow,
# so we pick a smaller value for crossover
MAX_MODULUS = 2**8

include "matrix_modn_dense_template.pxi"

cdef class Matrix_modn_dense_float(Matrix_modn_dense_template):
    r"""
    Dense matrices over `\ZZ/n\ZZ` for `n < 2^{11}` using LinBox's ``Modular<float>``

    These are matrices with integer entries mod ``n`` represented as
    floating-point numbers in a 32-bit word for use with LinBox routines.
    This allows for ``n`` up to `2^{11}`.  The
    ``Matrix_modn_dense_double`` class is used for larger moduli.

    Routines here are for the most basic access, see the
    `matrix_modn_dense_template.pxi` file for higher-level routines.
    """
    def __cinit__(self):
        """
        The Cython constructor

        TESTS::

            sage: A = random_matrix(GF(7), 4, 4)
            sage: type(A[0,0])
            <type 'sage.rings.finite_rings.integer_mod.IntegerMod_int'>
        """
        self._get_template = self._base_ring.zero()

    cdef set_unsafe_int(self, Py_ssize_t i, Py_ssize_t j, int value):
        r"""
        Set the (i,j) entry of self to the int value.

        EXAMPLES::

            sage: A = random_matrix(GF(7), 4, 4); A
            [3 1 6 6]
            [4 4 2 2]
            [3 5 4 5]
            [6 2 2 1]
            sage: A[0,0] = 12; A
            [5 1 6 6]
            [4 4 2 2]
            [3 5 4 5]
            [6 2 2 1]

            sage: B = random_matrix(Integers(100), 4, 4); B
            [13 95  1 16]
            [18 33  7 31]
            [92 19 18 93]
            [82 42 15 38]
            sage: B[0,0] = 422; B
            [22 95  1 16]
            [18 33  7 31]
            [92 19 18 93]
            [82 42 15 38]
        """
        self._matrix[i][j] = <float>value

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, x):
        r"""
        Set the (i,j) entry with no bounds-checking, or any other checks.

        Assumes that `x` is in the base ring.

        EXAMPLES::

            sage: A = random_matrix(GF(13), 4, 4); A
            [ 0  0  2  9]
            [10  6 11  8]
            [10 12  8  8]
            [ 3  6  8  0]
            sage: K = A.base_ring()
            sage: x = K(27)
            sage: A[0,0] = x; A
            [ 1  0  2  9]
            [10  6 11  8]
            [10 12  8  8]
            [ 3  6  8  0]

            sage: B = random_matrix(Integers(200), 4, 4); B
            [ 13  95 101 116]
            [118 133   7 131]
            [192  19 118 193]
            [ 82 142 115  38]
            sage: R = B.base_ring()
            sage: x = R(311)
            sage: B[0,0] = x; B
            [111  95 101 116]
            [118 133   7 131]
            [192  19 118 193]
            [ 82 142 115  38]
        """
        self._matrix[i][j] = <float>(<IntegerMod_int>x).ivalue

    cdef IntegerMod_int get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        r"""
        Return the (i,j) entry with no bounds-checking.

        OUTPUT:

        A :class:`sage.rings.finite_rings.integer_mod.IntegerMod_int`
        object.

        EXAMPLES::

            sage: A = random_matrix(Integers(100), 4, 4); A
            [ 4 95 83 47]
            [44 57 91 53]
            [75 53 15 39]
            [26 25 10 74]
            sage: a = A[0,0]; a
            4
            sage: a in A.base_ring()
            True

            sage: B = random_matrix(Integers(100), 4, 4); B
            [13 95  1 16]
            [18 33  7 31]
            [92 19 18 93]
            [82 42 15 38]
            sage: b = B[0,0]; b
            13
            sage: b in B.base_ring()
            True
        """
        cdef float result = (<Matrix_modn_dense_template>self)._matrix[i][j]
        return (<Matrix_modn_dense_float>self)._get_template._new_c(<int_fast32_t>result)
