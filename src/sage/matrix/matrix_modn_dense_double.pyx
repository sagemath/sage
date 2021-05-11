# distutils: language = c++
# distutils: libraries = CBLAS_LIBRARIES
# distutils: library_dirs = CBLAS_LIBDIR
# distutils: include_dirs = CBLAS_INCDIR
# distutils: extra_compile_args = -D_XPG6
"""
Dense matrices over `\ZZ/n\ZZ` for `n < 2^{23}` using LinBox's ``Modular<double>``

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
MAX_MODULUS = 2**23

from sage.rings.finite_rings.integer_mod cimport IntegerMod_int64

include "matrix_modn_dense_template.pxi"


cdef class Matrix_modn_dense_double(Matrix_modn_dense_template):
    r"""
    Dense matrices over `\ZZ/n\ZZ` for `n < 2^{23}` using LinBox's ``Modular<double>``

    These are matrices with integer entries mod ``n`` represented as
    floating-point numbers in a 64-bit word for use with LinBox routines.
    This allows for ``n`` up to `2^{23}`.  The analogous
    ``Matrix_modn_dense_float`` class is used for smaller moduli.

    Routines here are for the most basic access, see the
    ``matrix_modn_dense_template.pxi`` file for higher-level routines.
    """

    def __cinit__(self):
        """
        The Cython constructor

        TESTS::

            sage: A = random_matrix(IntegerModRing(2^16), 4, 4)
            sage: type(A[0,0])
            <type 'sage.rings.finite_rings.integer_mod.IntegerMod_int64'>
        """
        self._get_template = self._base_ring.zero()
        # note that INTEGER_MOD_INT32_LIMIT is ceil(sqrt(2^31-1)) < 2^23
        self._fits_int32 = ((<Matrix_modn_dense_template>self).p <= INTEGER_MOD_INT32_LIMIT)

    cdef set_unsafe_int(self, Py_ssize_t i, Py_ssize_t j, int value):
        r"""
        Set the (i,j) entry of self to the int value.

        EXAMPLES::

            sage: A = random_matrix(GF(3016963), 4, 4); A
            [ 220081 2824836  765701 2282256]
            [1795330  767112 2967421 1373921]
            [2757699 1142917 2720973 2877160]
            [1674049 1341486 2641133 2173280]
            sage: A[0,0] = 220082r; A
            [ 220082 2824836  765701 2282256]
            [1795330  767112 2967421 1373921]
            [2757699 1142917 2720973 2877160]
            [1674049 1341486 2641133 2173280]
            sage: a = A[0,0]; a
            220082
            sage: ~a
            2859358

            sage: A = random_matrix(Integers(5099106), 4, 4); A
            [2629491 1237101 2033003 3788106]
            [4649912 1157595 4928315 4382585]
            [4252686  978867 2601478 1759921]
            [1303120 1860486 3405811 2203284]
            sage: A[0,0] = 220082r; A
            [ 220082 1237101 2033003 3788106]
            [4649912 1157595 4928315 4382585]
            [4252686  978867 2601478 1759921]
            [1303120 1860486 3405811 2203284]
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

            sage: A = random_matrix(GF(3016963), 4, 4); A
            [ 220081 2824836  765701 2282256]
            [1795330  767112 2967421 1373921]
            [2757699 1142917 2720973 2877160]
            [1674049 1341486 2641133 2173280]
            sage: K = A.base_ring()
            sage: A[0,0] = K(220082); A
            [ 220082 2824836  765701 2282256]
            [1795330  767112 2967421 1373921]
            [2757699 1142917 2720973 2877160]
            [1674049 1341486 2641133 2173280]
            sage: a = A[0,0]; a
            220082
            sage: ~a
            2859358

            sage: A = random_matrix(Integers(5099106), 4, 4); A
            [2629491 1237101 2033003 3788106]
            [4649912 1157595 4928315 4382585]
            [4252686  978867 2601478 1759921]
            [1303120 1860486 3405811 2203284]
            sage: K = A.base_ring()
            sage: A[0,0] = K(220081)
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

            sage: A = random_matrix(GF(3016963), 4, 4); A
            [ 220081 2824836  765701 2282256]
            [1795330  767112 2967421 1373921]
            [2757699 1142917 2720973 2877160]
            [1674049 1341486 2641133 2173280]
            sage: a = A[0,0]; a
            220081
            sage: ~a
            697224
            sage: K = A.base_ring()
            sage: ~K(220081)
            697224

            sage: A = random_matrix(Integers(5099106), 4, 4); A
            [2629491 1237101 2033003 3788106]
            [4649912 1157595 4928315 4382585]
            [4252686  978867 2601478 1759921]
            [1303120 1860486 3405811 2203284]
            sage: a = A[0,1]; a
            1237101
            sage: a*a
            3803997
            sage: K = A.base_ring()
            sage: K(1237101)^2
            3803997
        """
        cdef Matrix_modn_dense_double _self = <Matrix_modn_dense_double>self
        cdef double result = (<Matrix_modn_dense_template>self)._matrix[i][j]
        if _self._fits_int32:
            return (<IntegerMod_int>_self._get_template)._new_c(<int_fast32_t>result)
        else:
            return (<IntegerMod_int64>_self._get_template)._new_c(<int_fast64_t>result)
