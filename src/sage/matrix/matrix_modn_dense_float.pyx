# sage_setup: distribution = sagemath-linbox
# distutils: language = c++
# distutils: libraries = CBLAS_LIBRARIES
# distutils: library_dirs = CBLAS_LIBDIR
# distutils: include_dirs = CBLAS_INCDIR
r"""
Dense matrices over `\ZZ/n\ZZ` for `n < 2^{8}` using LinBox's ``Modular<float>``

AUTHORS:
- Burcin Erocal
- Martin Albrecht
"""
# #############################################################################
#       Copyright (C) 2011 Burcin Erocal <burcin@erocal.org>
#       Copyright (C) 2011 Martin Albrecht <martinralbrecht@googlemail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.finite_rings.stdint cimport *

from sage.libs.linbox.givaro cimport \
    Modular_float as ModField, \
    Poly1Dom, Dense

from sage.libs.linbox.linbox cimport \
    DenseMatrix_Modular_float as DenseMatrix, \
    reducedRowEchelonize

from sage.libs.linbox.fflas cimport \
    fgemm, pfgemm, fgemv, Det, pDet, Rank, pRank, ReducedRowEchelonForm, pReducedRowEchelonForm, applyP, \
    MinPoly, CharPoly, \
    ModFloatDensePolynomial as ModDensePoly

ctypedef Poly1Dom[ModField, Dense] ModDensePolyRing

# LinBox supports up to 2^11 using float but that's double dog slow,
# so we pick a smaller value for crossover
MAX_MODULUS = 2**8

include "matrix_modn_dense_template.pxi"

cdef class Matrix_modn_dense_float(Matrix_modn_dense_template):
    r"""
    Dense matrices over `\ZZ/n\ZZ` for `n < 2^{8}` using LinBox's ``Modular<float>``.

    These are matrices with integer entries mod ``n`` represented as
    floating-point numbers in a 32-bit word for use with LinBox routines.
    This could allow for ``n`` up to `2^{11}`, but for performance reasons
    this is limited to ``n`` up to `2^{8}`, and for larger moduli the
    ``Matrix_modn_dense_double`` class is used.

    Routines here are for the most basic access, see the
    ``matrix_modn_dense_template.pxi`` file for higher-level routines.
    """
    def __cinit__(self):
        """
        The Cython constructor.

        TESTS::

            sage: A = random_matrix(GF(7), 4, 4)
            sage: type(A[0,0])
            <class 'sage.rings.finite_rings.integer_mod.IntegerMod_int'>
        """
        self._get_template = self._base_ring.zero()

    cdef void set_unsafe_int(self, Py_ssize_t i, Py_ssize_t j, int value) noexcept:
        r"""
        Set the (i,j) entry of ``self`` to the int value.

        EXAMPLES::

            sage: A = random_matrix(GF(7), 4, 4)
            sage: l = A.list()
            sage: A[0,0] = 12
            sage: A.list()[0] == 12
            True
            sage: l[1:] == A.list()[1:]
            True

            sage: B = random_matrix(Integers(100), 4, 4)
            sage: l = B.list()
            sage: B[0,0] = 422
            sage: B.list()[0] == 22
            True
            sage: l[1:] == B.list()[1:]
            True
        """
        self._matrix[i][j] = <float>value

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, x):
        r"""
        Set the (i,j) entry with no bounds-checking, or any other checks.

        Assumes that `x` is in the base ring.

        EXAMPLES::

            sage: A = random_matrix(GF(13), 4, 4)
            sage: l = A.list()
            sage: K = A.base_ring()
            sage: x = K(27)
            sage: A[0,0] = x
            sage: A[0,0] == x
            True
            sage: l[1:] == A.list()[1:]
            True

            sage: B = random_matrix(Integers(200), 4, 4)
            sage: l = B.list()
            sage: R = B.base_ring()
            sage: x = R(311)
            sage: B[0,0] = x
            sage: B.list()[0] == x
            True
            sage: l[1:] == B.list()[1:]
            True
        """
        self._matrix[i][j] = <float>(<IntegerMod_int>x).ivalue

    cdef IntegerMod_int get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        r"""
        Return the (i,j) entry with no bounds-checking.

        OUTPUT:

        A :class:`sage.rings.finite_rings.integer_mod.IntegerMod_int`
        object.

        EXAMPLES::

            sage: R = Integers(100)
            sage: l = [R.random_element() for _ in range(4*4)]
            sage: A = matrix(Integers(100), 4, 4, l)
            sage: a = A[0,0]
            sage: a == l[0]
            True
            sage: a in R
            True

            sage: l = [R.random_element() for _ in range(4*4)]
            sage: B = matrix(Integers(100), 4, 4, l)
            sage: b = B[0,0]
            sage: b == l[0]
            True
            sage: b in R
            True
        """
        cdef float result = (<Matrix_modn_dense_template>self)._matrix[i][j]
        return (<Matrix_modn_dense_float>self)._get_template._new_c(<int_fast32_t>result)
