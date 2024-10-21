# sage_setup: distribution = sagemath-linbox
# distutils: language = c++
# distutils: libraries = CBLAS_LIBRARIES
# distutils: library_dirs = CBLAS_LIBDIR
# distutils: include_dirs = CBLAS_INCDIR
# distutils: extra_compile_args = -D_XPG6
r"""
Dense matrices over `\ZZ/n\ZZ` for `n < 94906266` using LinBox's ``Modular<double>``

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
    Modular_double as ModField, \
    Poly1Dom, Dense

from sage.libs.linbox.linbox cimport \
    reducedRowEchelonize, \
    DenseMatrix_Modular_double as DenseMatrix

from sage.libs.linbox.fflas cimport \
    fgemm, pfgemm, fgemv, Det, pDet, Rank, pRank, ReducedRowEchelonForm, pReducedRowEchelonForm, applyP, \
    MinPoly, CharPoly, MinPoly, \
    ModDoubleDensePolynomial as ModDensePoly

ctypedef Poly1Dom[ModField, Dense] ModDensePolyRing

# Limit for LinBox Modular<double>
MAX_MODULUS = 94906266

from sage.rings.finite_rings.integer_mod cimport IntegerMod_int64

include "matrix_modn_dense_template.pxi"


cdef class Matrix_modn_dense_double(Matrix_modn_dense_template):
    r"""
    Dense matrices over `\ZZ/n\ZZ` for `n < 94906266` using LinBox's ``Modular<double>``.

    These are matrices with integer entries mod ``n`` represented as
    floating-point numbers in a 64-bit word for use with LinBox routines.
    This allows for ``n`` up to `94906266`. By default, the analogous
    ``Matrix_modn_dense_float`` class is used for smaller moduli, specifically
    for ``n`` up to `2^{8}`.

    Routines here are for the most basic access, see the
    ``matrix_modn_dense_template.pxi`` file for higher-level routines.
    """

    def __cinit__(self):
        """
        The Cython constructor.

        TESTS::

            sage: A = random_matrix(IntegerModRing(2^16), 4, 4)
            sage: type(A[0,0])
            <class 'sage.rings.finite_rings.integer_mod.IntegerMod_int64'>
        """
        self._get_template = self._base_ring.zero()
        # note that INTEGER_MOD_INT32_LIMIT is ceil(sqrt(2^31-1)) < 94906266
        self._fits_int32 = ((<Matrix_modn_dense_template>self).p <= INTEGER_MOD_INT32_LIMIT)

    cdef void set_unsafe_int(self, Py_ssize_t i, Py_ssize_t j, int value) noexcept:
        r"""
        Set the (i,j) entry of ``self`` to the int value.

        EXAMPLES::

            sage: A = random_matrix(GF(3016963), 4, 4)
            sage: l = A.list()
            sage: A[0,0] = 220082r
            sage: A.list()[1:] == l[1:]
            True
            sage: a = A[0,0]; a
            220082
            sage: ~a
            2859358

            sage: A = random_matrix(Integers(5099106), 4, 4)
            sage: l = A.list()
            sage: A[0,0] = 220082r
            sage: A.list()[1:] == l[1:]
            True
            sage: a = A[0,0]; a
            220082
            sage: a*a
            4777936
        """
        self._matrix[i][j] = <double>value

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, x):
        r"""
        Set the (i,j) entry with no bounds-checking, or any other checks.

        Assumes that `x` is in the base ring.

        EXAMPLES::

            sage: A = random_matrix(GF(3016963), 4, 4)
            sage: K = A.base_ring()
            sage: l = A.list()
            sage: A[0,0] = K(220082)
            sage: A.list()[1:] == l[1:]
            True
            sage: a = A[0,0]; a
            220082
            sage: ~a
            2859358

            sage: A = random_matrix(Integers(5099106), 4, 4)
            sage: K = A.base_ring()
            sage: l = A.list()
            sage: A[0,0] = K(220081)
            sage: A.list()[1:] == l[1:]
            True
            sage: a = A[0,0]; a
            220081
            sage: a*a
            4337773
        """
        if (<Matrix_modn_dense_double>self)._fits_int32:
            self._matrix[i][j] = <double>(<IntegerMod_int>x).ivalue
        else:
            self._matrix[i][j] = <double>(<IntegerMod_int64>x).ivalue

    cdef IntegerMod_abstract get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        r"""
        Return the (i,j) entry with no bounds-checking.

        OUTPUT:

        A :class:`sage.rings.finite_rings.integer_mod.IntegerMod_int`
        or
        :class:`sage.rings.finite_rings.integer_mod.IntegerMod_int64`
        object, depending on the modulus.

        EXAMPLES::

            sage: K = GF(3016963)
            sage: l = [K.random_element() for _ in range(4*4)]
            sage: A = matrix(K, 4, 4, l)
            sage: a = A[0,0]
            sage: a == l[0]
            True
            sage: a == 0 or ~a*a == 1
            True
            sage: a.parent() is K
            True

            sage: K = Integers(5099106)
            sage: l = [K.random_element() for _ in range(4*4)]
            sage: A = matrix(Integers(5099106), 4, 4, l)
            sage: a = A[0,1]
            sage: a == l[1]
            True
            sage: a*a == K(Integer(l[1]))^2
            True
        """
        cdef Matrix_modn_dense_double _self = <Matrix_modn_dense_double>self
        cdef double result = (<Matrix_modn_dense_template>self)._matrix[i][j]
        if _self._fits_int32:
            return (<IntegerMod_int>_self._get_template)._new_c(<int_fast32_t>result)
        else:
            return (<IntegerMod_int64>_self._get_template)._new_c(<int_fast64_t>result)
