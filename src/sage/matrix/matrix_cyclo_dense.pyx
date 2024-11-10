# distutils: language = c++
# distutils: libraries = NTL_LIBRARIES
# distutils: extra_compile_args = NTL_CFLAGS
# distutils: include_dirs = NTL_INCDIR
# distutils: library_dirs = NTL_LIBDIR
# distutils: extra_link_args = NTL_LIBEXTRA
"""
Matrices over Cyclotomic Fields

The underlying matrix for a Matrix_cyclo_dense object is stored as
follows: given an n x m matrix over a cyclotomic field of degree d, we
store a d x (nm) matrix over QQ, each column of which corresponds to
an element of the original matrix. This can be retrieved via the
_rational_matrix method. Here is an example illustrating this:

EXAMPLES::

    sage: F.<zeta> = CyclotomicField(5)
    sage: M = Matrix(F, 2, 3, [zeta, 3, zeta**4+5, (zeta+1)**4, 0, 1])
    sage: M
    [                        zeta                            3  -zeta^3 - zeta^2 - zeta + 4]
    [3*zeta^3 + 5*zeta^2 + 3*zeta                            0                            1]

    sage: M._rational_matrix()
    [ 0  3  4  0  0  1]
    [ 1  0 -1  3  0  0]
    [ 0  0 -1  5  0  0]
    [ 0  0 -1  3  0  0]


AUTHORS:
   * William Stein
   * Craig Citro
"""

#*****************************************************************************
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from cysignals.signals cimport sig_on, sig_off

include "sage/libs/ntl/decl.pxi"

from sage.structure.element cimport Element
from sage.misc.randstate cimport randstate, current_randstate
from sage.libs.gmp.randomize cimport *

from sage.libs.flint.types cimport fmpz_t
from sage.libs.flint.fmpz cimport fmpz_init, fmpz_clear, fmpz_set_mpz, fmpz_one, fmpz_get_mpz, fmpz_add, fmpz_mul, fmpz_sub, fmpz_mul_si, fmpz_mul_si, fmpz_mul_si, fmpz_divexact, fmpz_lcm
from sage.libs.flint.fmpq cimport fmpq_is_zero, fmpq_set_mpq, fmpq_canonicalise
from sage.libs.flint.fmpq_mat cimport fmpq_mat_entry_num, fmpq_mat_entry_den, fmpq_mat_entry

from sage.matrix.args cimport MatrixArgs_init
from sage.matrix.constructor import matrix
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.matrix cimport Matrix
from sage.matrix import matrix_dense
from sage.matrix.matrix_integer_dense cimport _lift_crt
from sage.structure.element cimport Matrix as baseMatrix
from sage.matrix.misc_flint import matrix_integer_dense_rational_reconstruction

from sage.arith.misc import binomial, previous_prime
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
from sage.rings.real_mpfr import create_RealNumber as RealNumber
from sage.rings.integer cimport Integer
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.number_field.number_field_element cimport NumberFieldElement
from sage.rings.number_field.number_field_element_quadratic cimport NumberFieldElement_quadratic

from sage.structure.proof.proof import get_flag as get_proof_flag
from sage.misc.verbose import verbose

from sage.matrix.matrix_modn_dense_double import MAX_MODULUS as MAX_MODULUS_modn_dense_double
from sage.arith.multi_modular import MAX_MODULUS as MAX_MODULUS_multi_modular
MAX_MODULUS = min(MAX_MODULUS_modn_dense_double, MAX_MODULUS_multi_modular)

# parameters for tuning
echelon_primes_increment = 15
echelon_verbose_level = 1


cdef class Matrix_cyclo_dense(Matrix_dense):
    def __cinit__(self):
        """
        EXAMPLES::

            sage: from sage.matrix.matrix_cyclo_dense import Matrix_cyclo_dense
            sage: A = Matrix_cyclo_dense.__new__(Matrix_cyclo_dense, MatrixSpace(CyclotomicField(3),2), [0,1,2,3], True, True)
            sage: type(A)
            <class 'sage.matrix.matrix_cyclo_dense.Matrix_cyclo_dense'>

        Note that the entries of A haven't even been set yet above; that doesn't
        happen until ``__init__`` is called::

            sage: A[0,0]
            Traceback (most recent call last):
            ...
            ValueError: matrix entries not yet initialized
        """
        self._degree = self._base_ring.degree()
        self._n = int(self._base_ring._n())

    def __init__(self, parent, entries=None, copy=None, bint coerce=True):
        """
        Initialize a newly created cyclotomic matrix.

        INPUT:

        - ``parent`` -- a matrix space over a cyclotomic number field

        - ``entries`` -- see :func:`matrix`

        - ``copy`` -- ignored (for backwards compatibility)

        - ``coerce`` -- if ``False``, assume without checking that the
          entries lie in the base ring

        EXAMPLES:

        This function is called implicitly when you create new
        cyclotomic dense matrices::

            sage: W.<a> = CyclotomicField(100)
            sage: A = matrix(2, 3, [1, 1/a, 1-a,a, -2/3*a, a^19])
            sage: A
            [                        1 -a^39 + a^29 - a^19 + a^9                    -a + 1]
            [                        a                    -2/3*a                      a^19]
            sage: TestSuite(A).run()

        TESTS::

            sage: matrix(W, 2, 1, a)
            Traceback (most recent call last):
            ...
            TypeError: nonzero scalar matrix must be square

        We call __init__ explicitly below::

            sage: from sage.matrix.matrix_cyclo_dense import Matrix_cyclo_dense
            sage: A = Matrix_cyclo_dense.__new__(Matrix_cyclo_dense, MatrixSpace(CyclotomicField(3),2), [0,1,2,3], True, True)
            sage: A.__init__(MatrixSpace(CyclotomicField(3),2), [0,1,2,3], True, True)
            sage: A
            [0 1]
            [2 3]
        """
        ma = MatrixArgs_init(parent, entries)
        cdef list L = []
        for x in ma.iter(coerce):
            L += x.list()

        QQspace = MatrixSpace(QQ, self._nrows * self._ncols, self._degree)
        QQmat = Matrix_rational_dense(QQspace, L, False, False)
        self._matrix = QQmat.transpose()

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, value):
        """
        Set the ij-th entry of ``self``.

        WARNING: This function does no bounds checking whatsoever, as
        the name suggests. It also assumes certain facts about the
        internal representation of cyclotomic fields. This is intended
        for internal use only.

        EXAMPLES::

            sage: K.<z> = CyclotomicField(11) ; M = Matrix(K,2,range(4))
            sage: M[0,1] = z ; M
            [0 z]
            [2 3]

            sage: K.<z> = CyclotomicField(3) ; M = Matrix(K,2,range(4))
            sage: M[1,1] = z+1 ; M
            [    0     1]
            [    2 z + 1]

        TESTS:

        Since separate code exists for each quadratic field, we need
        doctests for each.::

            sage: K.<z> = CyclotomicField(4) ; M = Matrix(K,2,range(4))
            sage: M[1,1] = z+1 ; M
            [    0     1]
            [    2 z + 1]
            sage: K.<z> = CyclotomicField(6) ; M = Matrix(K,2,range(4))
            sage: M[1,1] = z+1 ; M
            [    0     1]
            [    2 z + 1]
        """
        # NEW FAST VERSION -- makes assumptions about how the
        # cyclotomic field is implemented.
        cdef Py_ssize_t k, c
        cdef NumberFieldElement v
        cdef NumberFieldElement_quadratic qv
        cdef mpz_t numer, denom
        cdef fmpz_t ftmp

        # The i,j entry is the (i * self._ncols + j)'th column.
        c = i * self._ncols + j

        if self._degree == 2:
            # Must be coded differently, since elements of
            # quadratic number fields are stored differently.
            qv = <NumberFieldElement_quadratic> value
            if self._n == 4:
                fmpz_set_mpz(fmpq_mat_entry_num(self._matrix._matrix, 0, c),
                             qv.a)
                fmpz_set_mpz(fmpq_mat_entry_den(self._matrix._matrix, 0, c),
                             qv.denom)
                fmpq_canonicalise(fmpq_mat_entry(self._matrix._matrix, 0, c))

                fmpz_set_mpz(fmpq_mat_entry_num(self._matrix._matrix, 1, c),
                             qv.b)
                fmpz_set_mpz(fmpq_mat_entry_den(self._matrix._matrix, 1, c),
                             qv.denom)
                fmpq_canonicalise(fmpq_mat_entry(self._matrix._matrix, 1, c))
            elif self._n == 3:
                # mat[0,c] = (x.a + x.b) / x.denom
                fmpz_set_mpz(fmpq_mat_entry_num(self._matrix._matrix, 0, c),
                        qv.a)

                # NOTE: it would be convenient here to have fmpz_add_mpz
                fmpz_init(ftmp)
                fmpz_set_mpz(ftmp, qv.b)
                fmpz_add(fmpq_mat_entry_num(self._matrix._matrix, 0, c),
                         fmpq_mat_entry_num(self._matrix._matrix, 0, c),
                         ftmp)
                fmpz_clear(ftmp)

                fmpz_set_mpz(fmpq_mat_entry_den(self._matrix._matrix, 0, c),
                             qv.denom)
                fmpq_canonicalise(fmpq_mat_entry(self._matrix._matrix, 0, c))

                # mat[1,c] = (2 x.b) / x.denom
                fmpz_set_mpz(fmpq_mat_entry_num(self._matrix._matrix, 1, c),
                             qv.b)
                fmpz_mul_si(fmpq_mat_entry_num(self._matrix._matrix, 1, c),
                            fmpq_mat_entry_num(self._matrix._matrix, 1, c), 2)
                fmpz_set_mpz(fmpq_mat_entry_den(self._matrix._matrix, 1, c),
                             qv.denom)
                fmpq_canonicalise(fmpq_mat_entry(self._matrix._matrix, 1, c))
            else: # self._n is 6
                fmpz_set_mpz(fmpq_mat_entry_num(self._matrix._matrix, 0, c),
                             qv.a)

                # NOTE: it would be convenient here to have fmpz_add_mpz
                fmpz_init(ftmp)
                fmpz_set_mpz(ftmp, qv.b)
                fmpz_sub(fmpq_mat_entry_num(self._matrix._matrix, 0, c),
                         fmpq_mat_entry_num(self._matrix._matrix, 0, c),
                         ftmp)
                fmpz_clear(ftmp)

                fmpz_set_mpz(fmpq_mat_entry_den(self._matrix._matrix, 0, c),
                             qv.denom)
                fmpq_canonicalise(fmpq_mat_entry(self._matrix._matrix, 0, c))

                fmpz_set_mpz(fmpq_mat_entry_num(self._matrix._matrix, 1, c),
                             qv.b)
                fmpz_mul_si(fmpq_mat_entry_num(self._matrix._matrix, 1, c),
                            fmpq_mat_entry_num(self._matrix._matrix, 1, c), 2)
                fmpz_set_mpz(fmpq_mat_entry_den(self._matrix._matrix, 1, c),
                             qv.denom)
                fmpq_canonicalise(fmpq_mat_entry(self._matrix._matrix, 1, c))
            return

        v = value

        mpz_init(numer)
        mpz_init(denom)

        v._ntl_denom_as_mpz(denom)
        for k in range(self._degree):
            v._ntl_coeff_as_mpz(numer, k)
            fmpz_set_mpz(fmpq_mat_entry_num(self._matrix._matrix, k, c), numer)
            fmpz_set_mpz(fmpq_mat_entry_den(self._matrix._matrix, k, c), denom)
            fmpq_canonicalise(fmpq_mat_entry(self._matrix._matrix, k, c))

        mpz_clear(numer)
        mpz_clear(denom)

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        """
        Get the ij-th of ``self``.

        WARNING: As the name suggests, expect segfaults if i,j are out
        of bounds!! This is for internal use only.

        EXAMPLES::

            sage: W.<a> = CyclotomicField(5)
            sage: A = matrix(2, 3, [9939208341, 1/a, 1-a,a, -2/3*a, a^19])

        This implicitly calls get_unsafe::

            sage: A[0,0]
            9939208341

        TESTS:

        Since separate code exists for each quadratic field, we need
        doctests for each.::

            sage: K.<z> = CyclotomicField(3) ; M = Matrix(K,2,range(4))
            sage: M[1,1] = z+1 ; M[1,1]
            z + 1
            sage: (M[1,1] - M[0,1])**3
            1
            sage: K.<z> = CyclotomicField(4) ; M = Matrix(K,2,range(4))
            sage: M[1,1] = z+1 ; M[1,1]
            z + 1
            sage: (M[1,1] - M[0,1])**4
            1
            sage: K.<z> = CyclotomicField(6) ; M = Matrix(K,2,range(4))
            sage: M[1,1] = z+1 ; M[1,1]
            z + 1
            sage: (M[1,1] - M[0,1])**6
            1
        """
        cdef Py_ssize_t k, c
        cdef NumberFieldElement x
        cdef NumberFieldElement_quadratic xq
        cdef mpz_t tmp
        cdef fmpz_t denom, ftmp
        cdef ZZ_c coeff

        if self._matrix is None:
            raise ValueError("matrix entries not yet initialized")

        c = i * self._ncols + j
        mpz_init(tmp)

        if self._degree == 2:
            fmpz_init(ftmp)
            xq = self._base_ring(0)
            if self._n == 4:
                fmpz_mul(ftmp, fmpq_mat_entry_num(self._matrix._matrix, 0, c),
                               fmpq_mat_entry_den(self._matrix._matrix, 1, c))
                fmpz_get_mpz(xq.a, ftmp)
                fmpz_mul(ftmp, fmpq_mat_entry_num(self._matrix._matrix, 1, c),
                               fmpq_mat_entry_den(self._matrix._matrix, 0, c))
                fmpz_get_mpz(xq.b, ftmp)
                fmpz_mul(ftmp, fmpq_mat_entry_den(self._matrix._matrix, 0, c),
                               fmpq_mat_entry_den(self._matrix._matrix, 1, c))
                fmpz_get_mpz(xq.denom, ftmp)

            else: # n is 3 or 6
                fmpz_mul(ftmp, fmpq_mat_entry_num(self._matrix._matrix, 0, c),
                               fmpq_mat_entry_den(self._matrix._matrix, 1, c))
                fmpz_mul_si(ftmp, ftmp, 2)
                fmpz_get_mpz(xq.a, ftmp)
                fmpz_mul(ftmp, fmpq_mat_entry_den(self._matrix._matrix, 0, c),
                               fmpq_mat_entry_num(self._matrix._matrix, 1, c))
                fmpz_get_mpz(tmp, ftmp)
                if self._n == 3:
                    mpz_sub(xq.a, xq.a, tmp)
                else: # n == 6
                    mpz_add(xq.a, xq.a, tmp)

                fmpz_mul(ftmp, fmpq_mat_entry_den(self._matrix._matrix, 0, c),
                               fmpq_mat_entry_num(self._matrix._matrix, 1, c))
                fmpz_get_mpz(xq.b, ftmp)

                fmpz_mul(ftmp, fmpq_mat_entry_den(self._matrix._matrix, 0, c),
                               fmpq_mat_entry_den(self._matrix._matrix, 1, c))
                fmpz_get_mpz(xq.denom, ftmp)
                mpz_mul_si(xq.denom, xq.denom, 2)

            xq._reduce_c_()
            mpz_clear(tmp)
            fmpz_clear(ftmp)
            return xq

        x = self._base_ring(0)
        fmpz_init(denom)
        fmpz_init(ftmp)
        fmpz_one(denom)

        # Get the least common multiple of the denominators in
        # this column.
        for k in range(self._degree):
            fmpz_lcm(denom, denom, fmpq_mat_entry_den(self._matrix._matrix, k, c))

        for k in range(self._degree):
            # set each entry of x to a*denom/b where a/b is the
            # k,c entry of _matrix.
            fmpz_mul(ftmp, fmpq_mat_entry_num(self._matrix._matrix, k, c), denom)
            fmpz_divexact(ftmp, ftmp, fmpq_mat_entry_den(self._matrix._matrix, k, c))
            # Now set k-th entry of x's numerator to tmp
            fmpz_get_mpz(tmp, ftmp)
            mpz_to_ZZ(&coeff, tmp)
            ZZX_SetCoeff(x._numerator, k, coeff)

        # Set the denominator of x to denom.
        fmpz_get_mpz(tmp, denom)
        mpz_to_ZZ(&x._denominator, tmp)
        fmpz_clear(denom)
        mpz_clear(tmp)
        fmpz_clear(ftmp)

        return x

    cdef bint get_is_zero_unsafe(self, Py_ssize_t i, Py_ssize_t j) except -1:
        r"""
        Return 1 if the entry ``(i, j)`` is zero, otherwise 0.

        EXAMPLES::

            sage: K.<z> = CyclotomicField(3)
            sage: A = matrix(K, 4, 3, [0, -z, -2, -2*z + 2, 2*z, z, z, 1-z, 2+3*z, z, 1+z, 0])
            sage: A.zero_pattern_matrix()  # indirect doctest
            [1 0 0]
            [0 0 0]
            [0 0 0]
            [0 0 1]
        """
        cdef int a
        for a in range(self._degree):
            if not self._matrix.get_is_zero_unsafe(a, j+i*self._ncols):
                return False
        return True

    def _pickle(self):
        """
        Used for pickling matrices. This function returns the
        underlying data and pickle version.

        OUTPUT:

        - data; output of pickle
        - version; integer

        EXAMPLES::

            sage: K.<z> = CyclotomicField(3)
            sage: w = matrix(K, 3, 3, [0, -z, -2, -2*z + 2, 2*z, z, z, 1-z, 2+3*z])
            sage: w._pickle()
            (('0 0 -2 2 0 0 0 1 2 0 -1 0 -2 2 1 1 -1 3', 0), 0)
        """
        data = self._matrix._pickle()
        return data, 0

    def _unpickle(self, data, int version):
        """
        Called when unpickling matrices.

        INPUT:

        - ``data`` -- string
        - ``version`` -- integer

        This modifies ``self``.

        EXAMPLES::

            sage: K.<z> = CyclotomicField(3)
            sage: w = matrix(K, 3, 3, [0, -z, -2, -2*z + 2, 2*z, z, z, 1-z, 2+3*z])
            sage: data, version = w._pickle()
            sage: k = w.parent()(0)
            sage: k._unpickle(data, version)
            sage: k == w
            True
        """
        #self.check_mutability()
        if version == 0:
            self._matrix = Matrix_rational_dense(MatrixSpace(QQ, self._degree, self._nrows*self._ncols), None, False, False)
            self._matrix._unpickle(*data)  # data is (data, matrix_QQ_version)
        else:
            raise RuntimeError("unknown matrix version (=%s)" % version)

    ########################################################################
    # LEVEL 2 functionality
    # x * cdef _add_
    # x * cdef _sub_
    #   * cdef _mul_
    # x * cdef _lmul_    -- scalar multiplication
    # x * cpdef _richcmp_
    # x * __neg__
    #   * __invert__
    # x * __copy__
    #   * _multiply_classical
    #   * _list -- list of underlying elements (need not be a copy)
    #   * _dict -- sparse dictionary of underlying elements (need not be a copy)
    ########################################################################

    cpdef _add_(self, right):
        """
        Return the sum of two dense cyclotomic matrices.

        INPUT:

        - ``self``, ``right`` -- dense cyclotomic matrices with the same
          parents

        OUTPUT: a dense cyclotomic matrix

        EXAMPLES::

            sage: W.<z> = CyclotomicField(5)
            sage: A = matrix(2, 3, [1,z,z^2,z^3,z^4,2/3*z]); B = matrix(2, 3, [-1,2*z,3*z^2,5*z+1,z^4,1/3*z])
            sage: A + B
            [                       0                      3*z                    4*z^2]
            [           z^3 + 5*z + 1 -2*z^3 - 2*z^2 - 2*z - 2                        z]

        Verify directly that the above output is correct::

            sage: (A+B).list() == [A.list()[i]+B.list()[i] for i in range(6)]
            True
        """
        cdef Matrix_cyclo_dense A = Matrix_cyclo_dense.__new__(Matrix_cyclo_dense, self.parent(), None, None, None)
        # Fix the redundancy here.
        A._matrix = self._matrix + (<Matrix_cyclo_dense>right)._matrix
        return A

    cpdef _sub_(self, right):
        """
        Return the difference of two dense cyclotomic matrices.

        INPUT:

        - ``self``, ``right`` -- dense cyclotomic matrices with the same
          parent

        OUTPUT: a dense cyclotomic matrix

        EXAMPLES::

            sage: W.<z> = CyclotomicField(5)
            sage: A = matrix(2, 3, [1,z,z^2,z^3,z^4,2/3*z]); B = matrix(2, 3, [-1,2*z,3*z^2,5*z+1,z^4,1/3*z])
            sage: A - B
            [            2            -z        -2*z^2]
            [z^3 - 5*z - 1             0         1/3*z]

        Verify directly that the above output is correct::

            sage: (A-B).list() == [A.list()[i]-B.list()[i] for i in range(6)]
            True
        """
        cdef Matrix_cyclo_dense A = Matrix_cyclo_dense.__new__(Matrix_cyclo_dense, self.parent(), None, None, None)
        A._matrix = self._matrix - (<Matrix_cyclo_dense>right)._matrix
        return A

    cpdef _lmul_(self, Element right):
        """
        Multiply a dense cyclotomic matrix by a scalar.

        INPUT:

        - ``right`` -- scalar in the base cyclotomic field

        EXAMPLES::

            sage: W.<z> = CyclotomicField(5)
            sage: A = matrix(2, 3, [1,z,z^2,z^3,z^4,2/3*z])
            sage: (1 + z/3)*A
            [                      1/3*z + 1                     1/3*z^2 + z                   1/3*z^3 + z^2]
            [2/3*z^3 - 1/3*z^2 - 1/3*z - 1/3            -z^3 - z^2 - z - 2/3                 2/9*z^2 + 2/3*z]

        Verify directly that the above output is correct::

            sage: ((1+z/3)*A).list() == [(1+z/3)*x for x in A.list()]
            True
        """
        if right.is_one():
            return self
        elif right.is_zero():
            return self.parent().zero()

        # Create a new matrix object but with the _matrix attribute not initialized:
        cdef Matrix_cyclo_dense A = Matrix_cyclo_dense.__new__(Matrix_cyclo_dense,
                                               self.parent(), None, None, None)

        if right.is_rational():
            A._matrix = self._matrix._lmul_(right._rational_())
        else:
            # Multiply by nontrivial element of the cyclotomic number field
            # We do this by finding the matrix of this element, then left
            # multiplying by it.  We have to tweak the matrix some to
            # get the right basis, etc.
            T = right.matrix().transpose()
            A._matrix = T * self._matrix
        return A

    cdef _matrix_times_matrix_(self, baseMatrix right):
        """
        Return the product of two cyclotomic dense matrices.

        INPUT:

        - ``self``, ``right`` -- cyclotomic dense matrices with compatible
          parents (same base ring, and compatible dimensions for matrix
          multiplication)

        OUTPUT: cyclotomic dense matrix

        ALGORITHM:

        Use a multimodular algorithm that involves multiplying the two matrices
        modulo split primes.

        EXAMPLES::

            sage: W.<z> = CyclotomicField(5)
            sage: A = matrix(3, 3, [1,z,z^2,z^3,z^4,2/3*z,-3*z,z,2+z]); B = matrix(3, 3, [-1,2*z,3*z^2,5*z+1,z^4,1/3*z,2-z,3-z,5-z])
            sage: A*B
            [        -z^3 + 7*z^2 + z - 1       -z^3 + 3*z^2 + 2*z + 1              -z^3 + 25/3*z^2]
            [-2*z^3 - 5/3*z^2 + 1/3*z + 4           -z^3 - 8/3*z^2 - 2     -2/3*z^2 + 10/3*z + 10/3]
            [             4*z^2 + 4*z + 4               -7*z^2 + z + 7  -9*z^3 - 2/3*z^2 + 3*z + 10]

        Verify that the answer above is consistent with what the
        generic sparse matrix multiply gives (which is a different
        implementation).::

            sage: A*B == A.sparse_matrix()*B.sparse_matrix()
            True

            sage: N1 = Matrix(CyclotomicField(6), 1, [1])
            sage: cf6 = CyclotomicField(6) ; z6 = cf6.0
            sage: N2 = Matrix(CyclotomicField(6), 1, 5, [0,1,z6,-z6,-z6+1])
            sage: N1*N2
            [         0          1      zeta6     -zeta6 -zeta6 + 1]
            sage: N1 = Matrix(CyclotomicField(6), 1, [-1])
            sage: N1*N2
            [        0        -1    -zeta6     zeta6 zeta6 - 1]

        Verify that a degenerate case bug reported at :issue:`5974` is fixed.

            sage: K.<zeta6>=CyclotomicField(6); matrix(K,1,2) * matrix(K,2,[0, 1, 0, -2*zeta6, 0, 0, 1, -2*zeta6 + 1])
            [0 0 0 0]

        TESTS:

        This is from :issue:`8666`::

            sage: K.<zeta4> = CyclotomicField(4)
            sage: m = matrix(K, [125])
            sage: n = matrix(K, [186])
            sage: m*n
            [23250]
            sage: (-m)*n
            [-23250]
        """
        A, denom_self = self._matrix._clear_denom()
        B, denom_right = (<Matrix_cyclo_dense>right)._matrix._clear_denom()

        # conservative but correct estimate: 2 is there to account for the
        # sign of the entries
        bound = 1 + 2 * A.height() * B.height() * self._ncols

        n = self._base_ring._n()
        p = previous_prime(MAX_MODULUS)
        prod = 1
        v = []
        while prod <= bound:
            while (n >= 2 and p % n != 1) or denom_self % p == 0 or denom_right % p == 0:
                if p == 2:
                    raise RuntimeError("we ran out of primes in matrix multiplication.")
                p = previous_prime(p)
            prod *= p
            Amodp, _ = self._reductions(p)
            Bmodp, _ = right._reductions(p)
            _, S = self._reduction_matrix(p)
            X = Amodp[0]._matrix_from_rows_of_matrices([Amodp[i] * Bmodp[i] for i in range(len(Amodp))])
            v.append(S*X)
            p = previous_prime(p)
        M = matrix(ZZ, self._base_ring.degree(), self._nrows*right.ncols())
        _lift_crt(M, v)
        d = denom_self * denom_right
        if d == 1:
            M = M.change_ring(QQ)
        else:
            M = (1/d)*M
        cdef Matrix_cyclo_dense C = Matrix_cyclo_dense.__new__(Matrix_cyclo_dense,
                    MatrixSpace(self._base_ring, self._nrows, right.ncols()),
                                                               None, None, None)
        C._matrix = M
        return C

    cdef long _hash_(self) except -1:
        """
        Return hash of an immutable matrix.

        This raises a :exc:`TypeError` if input matrix is mutable.

        EXAMPLES:

        This is called implicitly by the hash function.::

            sage: W.<z> = CyclotomicField(5)
            sage: A = matrix(W, 2, 2, [1,z,-z,1+z/2])

        The matrix must be immutable.::

            sage: hash(A)
            Traceback (most recent call last):
            ...
            TypeError: mutable matrices are unhashable
            sage: A.set_immutable()

        Yes, this works::

            sage: hash(A)  # random
            3107179158321342168

        ::

            sage: W.<z> = CyclotomicField(5)
            sage: A = matrix(W, 2, 2, [1,2/3*z+z^2,-z,1+z/2])
            sage: hash(A)
            Traceback (most recent call last):
            ...
            TypeError: mutable matrices are unhashable
            sage: A.set_immutable()
            sage: A.__hash__()  # random
            2347601038649299176
        """
        return hash(self._matrix)

    cpdef _richcmp_(self, right, int op):
        """
        Implement comparison of two cyclotomic matrices with
        identical parents.

        INPUT:

        - ``self``, ``right`` -- matrices with same parent

        OUTPUT: boolean

        EXAMPLES::

            sage: W.<z> = CyclotomicField(5)
            sage: A = matrix(W, 2, 2, [1,z,-z,1+z/2])

        These implicitly call richcmp::

            sage: A == 5
            False
            sage: A < 100
            True

        This function is called implicitly when comparisons with matrices
        are done::

            sage: W.<z> = CyclotomicField(5)
            sage: A = matrix(W, 2, 2, [1,2/3*z+z^2,-z,1+z/2])
            sage: A == A
            True
            sage: A < 2*A
            True
            sage: A >= 2*A
            False
        """
        return self._matrix._richcmp_((<Matrix_cyclo_dense>right)._matrix, op)

    def __copy__(self):
        """
        Make a copy of this matrix.

        EXAMPLES:

        We create a cyclotomic matrix::

            sage: W.<z> = CyclotomicField(5)
            sage: A = matrix(W, 2, 2, [1,2/3*z+z^2,-z,1+z/2])

        We make a copy of A::

            sage: C = A.__copy__()

        We make another reference to A::

            sage: B = A

        Changing this reference changes A itself::

            sage: B[0,0] = 10
            sage: A[0,0]
            10

        Changing the copy does not change A::

            sage: C[0,0] = 20
            sage: C[0,0]
            20
            sage: A[0,0]
            10
        """
        cdef Matrix_cyclo_dense A = Matrix_cyclo_dense.__new__(Matrix_cyclo_dense, self.parent(), None, None, None)
        A._matrix = self._matrix.__copy__()
        return A

    def __neg__(self):
        """
        Return the negative of this matrix.

        OUTPUT: matrix

        EXAMPLES::

            sage: W.<z> = CyclotomicField(5)
            sage: A = matrix(W, 2, 2, [1,2/3*z+z^2,-z,1+z/2])
            sage: -A
            [          -1 -z^2 - 2/3*z]
            [           z   -1/2*z - 1]
            sage: A.__neg__()
            [          -1 -z^2 - 2/3*z]
            [           z   -1/2*z - 1]
        """
        cdef Matrix_cyclo_dense A = Matrix_cyclo_dense.__new__(Matrix_cyclo_dense, self.parent(), None, None, None)
        A._matrix = -self._matrix
        return A

    ########################################################################
    # LEVEL 3 functionality (Optional)
    #    * __deepcopy__
    #    * __invert__
    #    * Matrix windows -- only if you need strassen for that base
    #    * Other functions (list them here):
    #    * Specialized echelon form
    #    * tensor product
    ########################################################################

    def set_immutable(self):
        """
        Change this matrix so that it is immutable.

        EXAMPLES::

            sage: W.<z> = CyclotomicField(5)
            sage: A = matrix(W, 2, 2, [1,2/3*z+z^2,-z,1+z/2])
            sage: A[0,0] = 10
            sage: A.set_immutable()
            sage: A[0,0] = 20
            Traceback (most recent call last):
            ...
            ValueError: matrix is immutable; please change a copy instead (i.e., use copy(M) to change a copy of M).

        Note that there is no function to set a matrix to be mutable
        again, since such a function would violate the whole point.
        Instead make a copy, which is always mutable by default.::

            sage: A.set_mutable()
            Traceback (most recent call last):
            ...
            AttributeError: 'sage.matrix.matrix_cyclo_dense.Matrix_cyclo_dense' object has no attribute 'set_mutable'...
            sage: B = A.__copy__()
            sage: B[0,0] = 20
            sage: B[0,0]
            20
        """
        self._matrix.set_immutable()
        matrix_dense.Matrix_dense.set_immutable(self)

    def _rational_matrix(self):
        """
        Return the underlying rational matrix corresponding to ``self``.

        EXAMPLES::

            sage: Matrix(CyclotomicField(7),4,4,range(16))._rational_matrix()
            [ 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
            sage: Matrix(CyclotomicField(7),4,4,[CyclotomicField(7).gen(0)**i for i in range(16)])._rational_matrix()
            [ 1  0  0  0  0  0 -1  1  0  0  0  0  0 -1  1  0]
            [ 0  1  0  0  0  0 -1  0  1  0  0  0  0 -1  0  1]
            [ 0  0  1  0  0  0 -1  0  0  1  0  0  0 -1  0  0]
            [ 0  0  0  1  0  0 -1  0  0  0  1  0  0 -1  0  0]
            [ 0  0  0  0  1  0 -1  0  0  0  0  1  0 -1  0  0]
            [ 0  0  0  0  0  1 -1  0  0  0  0  0  1 -1  0  0]
        """
        return self._matrix

    def denominator(self):
        """
        Return the denominator of the entries of this matrix.

        OUTPUT: integer; the smallest integer `d` so that ``d * self`` has
        entries in the ring of integers

        EXAMPLES::

            sage: W.<z> = CyclotomicField(5)
            sage: A = matrix(W, 2, 2, [-2/7,2/3*z+z^2,-z,1+z/19]); A
            [       -2/7 z^2 + 2/3*z]
            [         -z  1/19*z + 1]
            sage: d = A.denominator(); d
            399
        """
        return self._matrix.denominator()

    def coefficient_bound(self):
        r"""
        Return an upper bound for the (complex) absolute values of all
        entries of ``self`` with respect to all embeddings.

        Use ``self.height()`` for a sharper bound.

        This is computed using just the Cauchy-Schwarz inequality, i.e.,
        we use the fact that ::

             \left| \sum_i a_i\zeta^i \right| \leq \sum_i |a_i|,

        as `|\zeta| = 1`.

        EXAMPLES::

            sage: W.<z> = CyclotomicField(5)
            sage: A = matrix(W, 2, 2, [1+z, 0, 9*z+7, -3 + 4*z]); A
            [  z + 1       0]
            [9*z + 7 4*z - 3]
            sage: A.coefficient_bound()
            16

        The above bound is just `9 + 7`, coming from the lower left entry.
        A better bound would be the following::

            sage: (A[1,0]).abs()
            12.997543663...
        """
        cdef Py_ssize_t i, j

        bound = 0
        for i from 0 <= i < self._matrix._ncols:

            n = 0
            for j from 0 <= j < self._matrix._nrows:
                n += self._matrix[j, i].abs()
            if bound < n:
                bound = n

        return bound

    def height(self):
        r"""
        Return the height of ``self``.

        If we let `a_{ij}` be the `i,j` entry of self, then we define
        the height of ``self`` to be

            `\max_v \max_{i,j} |a_{ij}|_v`,

        where `v` runs over all complex embeddings of ``self.base_ring()``.

        EXAMPLES::

            sage: W.<z> = CyclotomicField(5)
            sage: A = matrix(W, 2, 2, [1+z, 0, 9*z+7, -3 + 4*z]); A
            [  z + 1       0]
            [9*z + 7 4*z - 3]
            sage: A.height()
            12.997543663...
            sage: (A[1,0]).abs()
            12.997543663...
        """
        cdef Py_ssize_t i, j

        emb = self._base_ring.complex_embeddings()

        ht = 0
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                t = max([ x.norm().sqrt() for x in [ f(self.get_unsafe(i,j)) for f in emb ] ])
                if t > ht:
                    ht = t

        return ht

    cdef _randomize_rational_column_unsafe(Matrix_cyclo_dense self,
        Py_ssize_t col, mpz_t nump1, mpz_t denp1, distribution=None):
        """
        Randomizes all entries in column ``col``.  This is a helper method
        used in the implementation of dense matrices over cyclotomic fields.

        INPUT:

        - ``col`` -- integer indicating the column; must be coercible to
          ``int``, and this must lie between 0 (inclusive) and
          ``self._ncols`` (exclusive), since no bounds-checking is performed
        - ``nump1`` -- integer; numerator bound plus one
        - ``denp1`` -- integer; denominator bound plus one
        - ``distribution`` -- ``None`` or '1/n' (default: ``None``); if '1/n'
          then ``num_bound``, ``den_bound`` are ignored and numbers are chosen
          using the GMP function ``mpq_randomize_entry_recip_uniform``
        - ``nonzero`` -- boolean (default: ``False``); whether the new entries
          are forced to be nonzero

        OUTPUT: none, the matrix is modified in-space

        WARNING:

        This method is quite unsafe.  It's called from the method
        ``randomize``, but probably shouldn't be called from another method
        without first carefully reading the source code!

        TESTS:

        The following doctests are all indirect::

            sage: MS = MatrixSpace(CyclotomicField(10), 4, 4)
            sage: A = MS.random_element(); A  # random
            [    -2*zeta10^3 + 2*zeta10^2 - zeta10    zeta10^3 + 2*zeta10^2 - zeta10 + 1                                     0 -2*zeta10^3 + zeta10^2 - 2*zeta10 + 2]
            [                                    0       -zeta10^3 + 2*zeta10^2 - zeta10                         -zeta10^3 + 1                     zeta10^3 + zeta10]
            [                         1/2*zeta10^2                       -2*zeta10^2 + 2      -1/2*zeta10^3 + 1/2*zeta10^2 + 2             2*zeta10^3 - zeta10^2 - 2]
            [                                    1                          zeta10^2 + 2                            2*zeta10^2                          2*zeta10 - 2]
            sage: A.parent() is MS
            True

        ::

            sage: B = MS.random_element(density=0.5)
            sage: all(a in (-2, -1, -1/2, 0, 1/2, 1, 2) for a in B._rational_matrix().list())
            True
            sage: while set(B._rational_matrix().list()) != set((-2, -1, -1/2, 0, 1/2, 1, 2)):
            ....:     B = MS.random_element(density=0.5)

        ::

            sage: C = MS.random_element(density=0.5, num_bound=20, den_bound=20)
            sage: all(abs(a.denominator()) <= 20 and abs(a.numerator()) <= 20 for a in C._rational_matrix().list())
            True
            sage: while not (any(abs(a.denominator()) == 20 for a in C._rational_matrix().list())
            ....:         and any(abs(a.numerator()) == 20 for a in C._rational_matrix().list())):
            ....:     C = MS.random_element(density=0.5, num_bound=20, den_bound=20)
        """
        cdef Py_ssize_t i
        cdef Matrix_rational_dense mat = self._matrix
        cdef mpq_t tmp

        sig_on()
        mpq_init(tmp)
        if distribution == "1/n":
            for i in range(mat._nrows):
                mpq_randomize_entry_recip_uniform(tmp)
                fmpq_set_mpq(fmpq_mat_entry(mat._matrix, i, col), tmp)
        elif mpz_cmp_si(denp1, 2):   # denom is > 1
            for i in range(mat._nrows):
                mpq_randomize_entry(tmp, nump1, denp1)
                fmpq_set_mpq(fmpq_mat_entry(mat._matrix, i, col), tmp)
        else:
            for i in range(mat._nrows):
                mpq_randomize_entry_as_int(tmp, nump1)
                fmpq_set_mpq(fmpq_mat_entry(mat._matrix, i, col), tmp)
        mpq_clear(tmp)
        sig_off()

    def randomize(self, density=1, num_bound=2, den_bound=2,
                  distribution=None, nonzero=False, *args, **kwds):
        r"""
        Randomize the entries of ``self``.

        Choose rational numbers according to ``distribution``, whose
        numerators are bounded by ``num_bound`` and whose denominators are
        bounded by ``den_bound``.

        EXAMPLES::

            sage: A = Matrix(CyclotomicField(5),2,2,range(4)) ; A
            [0 1]
            [2 3]
            sage: A.randomize()
            sage: A   # random output
            [       1/2*zeta5^2 + zeta5                        1/2]
            [        -zeta5^2 + 2*zeta5 -2*zeta5^3 + 2*zeta5^2 + 2]
        """
        # Problem 1:
        # We cannot simply call the ``randomize`` code in ``matrix2.pyx`` on
        # the underlying matrix, since this is a d x (mn) matrix, where d is
        # the degree of the field extension, which leads to an overly dense
        # matrix.
        #
        # Problem 2:
        # We cannot simply copy the code from ``matrix2.pyx``, since the
        # ``random_element`` method for cyclotomic fields does not support
        # the arguments ``num_bound`` and ``den_bound``, which are support by
        # the rational field.
        #
        # Proposed solution:
        # Randomly select a proportion of ``density`` of the elements in the
        # matrix over the cyclotomic field, that is, this many columns in the
        # underlying rational matrix.  Then, for each element in that column,
        # randomize it to a rational number, applying the arguments
        # ``num_bound`` and ``den_bound``.

        density = float(density)
        if density <= 0:
            return
        if density > 1:
            density = 1

        self.check_mutability()
        self.clear_cache()

        cdef Py_ssize_t col, i, k, num
        cdef randstate rstate = current_randstate()
        cdef Integer B, C
        cdef bint col_is_zero

        B = Integer(num_bound+1)
        C = Integer(den_bound+1)

        if nonzero:
            if density >= 1:
                for col in range(self._matrix._ncols):
                    col_is_zero = True
                    while col_is_zero:
                        self._randomize_rational_column_unsafe(col, B.value,
                            C.value, distribution)
                        # Check whether the new column is nonzero
                        for i in range(self._degree):
                            if not fmpq_is_zero(fmpq_mat_entry(self._matrix._matrix, i, col)):
                                col_is_zero = False
                                break
            else:
                num = int(self._nrows * self._ncols * density)
                for k in range(num):
                    col = rstate.c_random() % self._matrix._ncols
                    col_is_zero = True
                    while col_is_zero:
                        self._randomize_rational_column_unsafe(col, B.value,
                            C.value, distribution)
                        # Check whether the new column is nonzero
                        for i in range(self._degree):
                            if not fmpq_is_zero(fmpq_mat_entry(self._matrix._matrix, i, col)):
                                col_is_zero = False
                                break
        else:
            if density >= 1:
                for col in range(self._matrix._ncols):
                    self._randomize_rational_column_unsafe(col, B.value,
                        C.value, distribution)
            else:
                num = int(self._nrows * self._ncols * density)
                for k in range(num):
                    col = rstate.c_random() % self._matrix._ncols
                    self._randomize_rational_column_unsafe(col, B.value,
                        C.value, distribution)

    def _charpoly_bound(self):
        """
        Determine a bound for the coefficients of the characteristic
        polynomial of ``self``.

        We use the bound in Lemma 2.2 of:

          Dumas, J-G. "Bounds on the coefficients of characteristic
          and minimal polynomials." J. Inequal. Pure Appl. Math. 8
          (2007), no. 2.

        This bound only applies for ``self._nrows >= 4``, so in all
        smaller cases, we just use a naive bound.

        EXAMPLES::

            sage: A = Matrix(CyclotomicField(7),3,3,range(9))
            sage: A._charpoly_bound()
            2048
            sage: A.charpoly()
            x^3 - 12*x^2 - 18*x

        An example from the above paper, where the bound is sharp::

            sage: B = Matrix(CyclotomicField(7), 5,5, [1,1,1,1,1,1,1,-1,-1,-1,1,-1,1,-1,-1,1,-1,-1,1,-1,1,-1,-1,-1,1])
            sage: B._charpoly_bound()
            80
            sage: B.charpoly()
            x^5 - 5*x^4 + 40*x^2 - 80*x + 48
        """
        cdef Py_ssize_t i

        # should we even bother with this check, or just say in
        # the docstring that we assume it's square?
        if self._nrows != self._ncols:
            raise ArithmeticError("self must be a square matrix")

        if self.is_zero():
            return 1

        B = self.coefficient_bound()

        # TODO: should charpoly just hardcode the return value for
        # self.nrows() < 4?

        # this bound is only valid for n >= 4, use naive bounds
        # in other cases.
        if self._nrows <= 3:
            return max(1, 3*B, 6*B**2, 4*B**3)

        # This is 2*e^(1-(2(7\gamma-4))/(13(3-2\gamma))), where \gamma
        # is Euler's constant.
        delta = RealNumber('5.418236')
        # This is an approximation to 1/2. :)
        half = RealNumber('0.5')

        D = (((1+2*delta*self._nrows*(B**2)).sqrt()-1)/(delta*B**2)).ceil()

        # TODO: we don't check anything about overflows anywhere here;
        # should we?

        # i = 0 case
        M = ((self._nrows * B**2)**(self._nrows * half)).ceil()

        for i from 1 <= i < D:
            val = binomial(self._nrows, i) * \
                  (((self._nrows-i)*B**2)**((self._nrows-i)*half)).ceil()
            if val > M:
                M = val

        return M

    def charpoly(self, var='x', algorithm='multimodular', proof=None):
        r"""
        Return the characteristic polynomial of self, as a polynomial
        over the base ring.

        INPUT:

        - ``algorithm`` -- options:

            - ``'multimodular'`` (default): reduce modulo primes, compute
              charpoly mod p, and lift (very fast)
            - ``'pari'``: use pari (quite slow; comparable to Magma v2.14 though)
            - ``'hessenberg'``: put matrix in Hessenberg form (double dog slow)

        - ``proof`` -- boolean (default: ``None``); proof flag determined by
          global linalg proof

        OUTPUT: polynomial

        EXAMPLES::

            sage: K.<z> = CyclotomicField(5)
            sage: a = matrix(K, 3, [1,z,1+z^2, z/3,1,2,3,z^2,1-z])
            sage: f = a.charpoly(); f
            x^3 + (z - 3)*x^2 + (-16/3*z^2 - 2*z)*x - 2/3*z^3 + 16/3*z^2 - 5*z + 5/3
            sage: f(a)
            [0 0 0]
            [0 0 0]
            [0 0 0]
            sage: a.charpoly(algorithm='pari')
            x^3 + (z - 3)*x^2 + (-16/3*z^2 - 2*z)*x - 2/3*z^3 + 16/3*z^2 - 5*z + 5/3
            sage: a.charpoly(algorithm='hessenberg')
            x^3 + (z - 3)*x^2 + (-16/3*z^2 - 2*z)*x - 2/3*z^3 + 16/3*z^2 - 5*z + 5/3

            sage: Matrix(K, 1, [0]).charpoly()
            x
            sage: Matrix(K, 1, [5]).charpoly(var='y')
            y - 5

            sage: Matrix(CyclotomicField(13),3).charpoly()
            x^3
            sage: Matrix(CyclotomicField(13),3).charpoly()[2].parent()
            Cyclotomic Field of order 13 and degree 12

        TESTS::

            sage: Matrix(CyclotomicField(10),0).charpoly()
            1
        """
        key = 'charpoly-%s-%s'%(algorithm,proof)
        f = self.fetch(key)
        if f is not None:
            return f.change_variable_name(var)

        if self.nrows() != self.ncols():
            raise TypeError("self must be square")

        if self.is_zero():
            R = PolynomialRing(self.base_ring(), name=var)
            f = R.gen(0)**self.nrows()
            self.cache(key, f)
            return f

        if self.nrows() == 1:
            R = PolynomialRing(self.base_ring(), name=var)
            f = R.gen(0) - self[0,0]
            self.cache(key, f)
            return f

        if algorithm == 'multimodular':
            f = self._charpoly_multimodular(var, proof=proof)
        elif algorithm == 'pari':
            paripoly = self.__pari__().charpoly()
            f = self.base_ring()[var](paripoly)
        elif algorithm == 'hessenberg':
            f = self._charpoly_hessenberg(var)
        else:
            raise ValueError("unknown algorithm '%s'" % algorithm)
        self.cache(key, f)
        return f

    def _charpoly_mod(self, p):
        """
        Return the characteristic polynomial of self*denom modulo all
        primes over `p`.

        This is used internally by the multimodular charpoly algorithm.

        INPUT:

        - ``p`` -- a prime that splits completely

        OUTPUT: matrix over GF(p) whose columns correspond to the entries
        of all the characteristic polynomials of the reduction of ``self``
        modulo all the primes over `p`

        EXAMPLES::

            sage: W.<z> = CyclotomicField(5)
            sage: A = matrix(W, 2, 2, [1+z, 0, 9*z+7, -3 + 4*z]); A
            [  z + 1       0]
            [9*z + 7 4*z - 3]
            sage: A._charpoly_mod(11)
            [8 2 1]
            [1 6 0]
            [4 0 0]
            [0 0 0]
        """
        tm = verbose("Computing characteristic polynomial of cyclotomic matrix modulo %s." % p)
        # Reduce self modulo all primes over p
        R, _ = self._reductions(p)
        # Compute the characteristic polynomial of each reduced matrix
        F = [A.charpoly('x') for A in R]
        # Put the characteristic polynomials together as the rows of a mod-p matrix
        k = R[0].base_ring()
        S = matrix(k, len(F), self.nrows() + 1, [f.list() for f in F])
        # multiply by inverse of reduction matrix to lift
        _, L = self._reduction_matrix(p)
        X = L * S
        # Now the columns of the matrix X define the entries of the
        # charpoly modulo p.
        verbose("Finished computing charpoly mod %s." % p, tm)
        return X

    def _charpoly_multimodular(self, var='x', proof=None):
        """
        Compute the characteristic polynomial of ``self`` using a
        multimodular algorithm.

        INPUT:

        - ``proof`` -- boolean (default: global flag); if ``False``, compute
          using primes `p_i` until the lift modulo all primes up to `p_i` is
          the same as the lift modulo all primes up to `p_{i+3}` or the bound
          is reached

        EXAMPLES::

            sage: K.<z> = CyclotomicField(3)
            sage: A = matrix(3, [-z, 2*z + 1, 1/2*z + 2, 1, -1/2, 2*z + 2, -2*z - 2, -2*z - 2, 2*z - 1])
            sage: A._charpoly_multimodular()
            x^3 + (-z + 3/2)*x^2 + (17/2*z + 9/2)*x - 9/2*z - 23/2
            sage: A._charpoly_multimodular('T')
            T^3 + (-z + 3/2)*T^2 + (17/2*z + 9/2)*T - 9/2*z - 23/2
            sage: A._charpoly_multimodular('T', proof=False)
            T^3 + (-z + 3/2)*T^2 + (17/2*z + 9/2)*T - 9/2*z - 23/2

        TESTS:

        We test a degenerate case::

            sage: A = matrix(CyclotomicField(1),2,[1,2,3,4]); A.charpoly()
            x^2 - 5*x - 2
        """
        cdef Matrix_cyclo_dense A
        A = Matrix_cyclo_dense.__new__(Matrix_cyclo_dense, self.parent(),
                                       None, None, None)

        proof = get_proof_flag(proof, "linear_algebra")

        n = self._base_ring._n()
        p = previous_prime(MAX_MODULUS)
        prod = 1
        v = []
        #A, denom = self._matrix._clear_denom()
        # TODO: this might be stupidly slow
        denom = self._matrix.denominator()
        A._matrix = <Matrix_rational_dense>(denom*self._matrix)
        bound = A._charpoly_bound()
        L_last = 0
        while prod <= bound:
            while (n >= 2  and p % n != 1) or denom % p == 0:
                if p == 2:
                    raise RuntimeError("we ran out of primes in multimodular charpoly algorithm.")
                p = previous_prime(p)

            X = A._charpoly_mod(p)
            v.append(X)
            prod *= p
            p = previous_prime(p)

            # if we've used enough primes as determined by bound, or
            # if we've used 3 primes, we check to see if the result is
            # the same.
            if prod >= bound or (not proof and (len(v) % 3 == 0)):
                M = matrix(ZZ, self._base_ring.degree(), self._nrows+1)
                L = _lift_crt(M, v)
                if not proof and L == L_last:
                    break
                L_last = L

        # Now each column of L encodes a coefficient of the output polynomial,
        # with column 0 being the constant coefficient.
        K = self.base_ring()
        R = K[var]
        coeffs = [K(w.list()) for w in L.columns()]
        f = R(coeffs)

        # Rescale to account for denominator, if necessary
        if denom != 1:
            x = R.gen()
            f = f(x * denom) * (1 / (denom**f.degree()))

        return f

    def _reductions(self, p):
        """
        Compute the reductions modulo all primes over p of denom*self,
        where denom is the denominator of ``self``.

        INPUT:

        - ``p`` -- a prime that splits completely in the base cyclotomic field

        OUTPUT:

        - ``list`` -- of r distinct matrices modulo p, where r is
          the degree of the cyclotomic base field
        - ``denom`` -- integer

        EXAMPLES::

            sage: K.<z> = CyclotomicField(3)
            sage: w = matrix(K, 2, 3, [0, -z/5, -2/3, -2*z + 2, 2*z, z])
            sage: R, d = w._reductions(7)
            sage: R[0]
            [0 2 4]
            [1 1 4]
            sage: R[1]
            [0 1 4]
            [5 4 2]
            sage: d
            15
        """
        # Get matrix that defines the linear reduction maps modulo
        # each prime of the base ring over p.
        T, _ = self._reduction_matrix(p)
        # Clear denominator and get matrix over the integers suitable
        # for reduction.
        A, denom = self._matrix._clear_denom()
        # Actually reduce the matrix over the integers modulo the
        # prime p.
        B = A._mod_int(p)
        # Now multiply, which computes from B all the reductions of
        # self*denom modulo each of the primes over p.
        R = T * B
        # Finally compute the actual reductions by extracting them
        # from R (note that the rows of R define the reductions).
        ans = R._matrices_from_rows(self._nrows, self._ncols)
        return ans, denom

    def _reduction_matrix(self, p):
        """
        INPUT:

        - ``p`` -- a prime that splits completely in the base field

        OUTPUT:

        - Matrix over GF(p) whose action from the left gives the map from O_K
          to GF(p) x ... x GF(p) given by reducing modulo all the primes over p
        - inverse of this matrix

        EXAMPLES::

            sage: K.<z> = CyclotomicField(3)
            sage: w = matrix(K, 2, 3, [0, -z/5, -2/3, -2*z + 2, 2*z, z])
            sage: A, B = w._reduction_matrix(7)
            sage: A
            [1 4]
            [1 2]
            sage: B
            [6 2]
            [4 3]

        The reduction matrix is used to calculate the reductions mod primes
        above p. ::

            sage: K.<z> = CyclotomicField(5)
            sage: A = matrix(K, 2, 2, [1, z, z^2+1, 5*z^3]); A
            [      1       z]
            [z^2 + 1   5*z^3]
            sage: T, S = A._reduction_matrix(11)
            sage: T * A._rational_matrix().change_ring(GF(11))
            [ 1  9  5  4]
            [ 1  5  4  9]
            [ 1  4  6  1]
            [ 1  3 10  3]

        The rows of this product are the (flattened) matrices mod each prime above p::

            sage: roots = [r for r, e in K.defining_polynomial().change_ring(GF(11)).roots()]; roots
            [9, 5, 4, 3]
            sage: [r^2+1 for r in roots]
            [5, 4, 6, 10]
            sage: [5*r^3 for r in roots]
            [4, 9, 1, 3]

        The reduction matrix is cached::

            sage: w._reduction_matrix(7) is w._reduction_matrix(7)
            True
        """
        cache = self.fetch('reduction_matrices')
        if cache is None:
            cache = {}
            self.cache('reduction_matrices', cache)
        try:
            return cache[p]
        except KeyError:
            pass
        K = self.base_ring()
        phi = K.defining_polynomial()
        from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
        from sage.matrix.constructor import matrix
        F = GF(p)
        aa = [a for a, _ in phi.change_ring(F).roots()]
        n = K.degree()
        if len(aa) != n:
            raise ValueError("the prime p (=%s) must split completely but doesn't" % p)
        T = matrix(F, n)
        for i in range(n):
            a = aa[i]
            b = 1
            for j in range(n):
                T[i,j] = b
                b *= a
        T.set_immutable()
        ans = (T, T**(-1))
        cache[p] = ans
        return ans

    def echelon_form(self, algorithm='multimodular', height_guess=None):
        """
        Find the echelon form of self, using the specified algorithm.

        The result is cached for each algorithm separately.

        EXAMPLES::

            sage: W.<z> = CyclotomicField(3)
            sage: A = matrix(W, 2, 3, [1+z, 2/3, 9*z+7, -3 + 4*z, z, -7*z]); A
            [  z + 1     2/3 9*z + 7]
            [4*z - 3       z    -7*z]
            sage: A.echelon_form()
            [                  1                   0  -192/97*z - 361/97]
            [                  0                   1 1851/97*z + 1272/97]
            sage: A.echelon_form(algorithm='classical')
            [                  1                   0  -192/97*z - 361/97]
            [                  0                   1 1851/97*z + 1272/97]

        We verify that the result is cached and that the caches are separate::

            sage: A.echelon_form() is A.echelon_form()
            True
            sage: A.echelon_form() is A.echelon_form(algorithm='classical')
            False

        TESTS::

            sage: W.<z> = CyclotomicField(13)
            sage: A = Matrix(W, 2,3, [10^30*(1-z)^13, 1, 2, 3, 4, z])
            sage: B = Matrix(W, 2,3, [(1-z)^13, 1, 2, 3, 4, z])
            sage: A.echelon_form() == A.echelon_form('classical')  # long time (4s on sage.math, 2011)
            True
            sage: B.echelon_form() == B.echelon_form('classical')
            True

        A degenerate case with the degree 1 cyclotomic field::

            sage: A = matrix(CyclotomicField(1),2,3,[1,2,3,4,5,6])
            sage: A.echelon_form()
            [ 1  0 -1]
            [ 0  1  2]

        A case that checks the bug in :issue:`3500`::

            sage: cf4 = CyclotomicField(4) ; z4 = cf4.0
            sage: A = Matrix(cf4, 1, 2, [-z4, 1])
            sage: A.echelon_form()
            [    1 zeta4]

        Verify that the matrix on :issue:`10281` works::

            sage: K.<rho> = CyclotomicField(106)
           sage: coeffs = [(18603/107*rho^51 - 11583/107*rho^50 - 19907/107*rho^49 - 13588/107*rho^48 - 8722/107*rho^47 + 2857/107*rho^46 - 19279/107*rho^45 - 16666/107*rho^44 - 11327/107*rho^43 + 3802/107*rho^42 + 18998/107*rho^41 - 10798/107*rho^40 + 16210/107*rho^39 - 13768/107*rho^38 + 15063/107*rho^37 - 14433/107*rho^36 - 19434/107*rho^35 - 12606/107*rho^34 + 3786/107*rho^33 - 17996/107*rho^32 + 12341/107*rho^31 - 15656/107*rho^30 - 19092/107*rho^29 + 8382/107*rho^28 - 18147/107*rho^27 + 14024/107*rho^26 + 18751/107*rho^25 - 8301/107*rho^24 - 20112/107*rho^23 - 14483/107*rho^22 + 4715/107*rho^21 + 20065/107*rho^20 + 15293/107*rho^19 + 10072/107*rho^18 + 4775/107*rho^17 - 953/107*rho^16 - 19782/107*rho^15 - 16020/107*rho^14 + 5633/107*rho^13 - 17618/107*rho^12 - 18187/107*rho^11 + 7492/107*rho^10 + 19165/107*rho^9 - 9988/107*rho^8 - 20042/107*rho^7 + 10109/107*rho^6 - 17677/107*rho^5 - 17723/107*rho^4 - 12489/107*rho^3 - 6321/107*rho^2 - 4082/107*rho - 1378/107, 1, 4*rho + 1), (0, 1, rho + 4)]
            sage: m = matrix(2, coeffs)
            sage: a = m.echelon_form(algorithm='classical')
            sage: b = m.echelon_form(algorithm='multimodular')  # long time (5s on sage.math, 2012)
            sage: a == b  # long time (depends on previous)
            True
        """
        key = 'echelon_form-%s'%algorithm
        E = self.fetch(key)
        if E is not None:
            return E

        if self._nrows == 0:
            E = self.__copy__()
            self.cache(key, E)
            self.cache('pivots', ())
            return E

        if algorithm == 'multimodular':
            E = self._echelon_form_multimodular(height_guess=height_guess)
        elif algorithm == 'classical':
            E = (self*self.denominator())._echelon_classical()
        else:
            raise ValueError("unknown algorithm '%s'" % algorithm)

        self.cache(key, E)
        return E

    def _echelon_form_multimodular(self, num_primes=10, height_guess=None):
        """
        Use a multimodular algorithm to find the echelon form of ``self``.

        INPUT:

        - ``num_primes`` -- number of primes to work modulo

        - ``height_guess`` -- guess for the height of the echelon form of self

        OUTPUT: matrix in reduced row echelon form

        EXAMPLES::

            sage: W.<z> = CyclotomicField(3)
            sage: A = matrix(W, 2, 3, [1+z, 2/3, 9*z+7, -3 + 4*z, z, -7*z]); A
            [  z + 1     2/3 9*z + 7]
            [4*z - 3       z    -7*z]
            sage: A._echelon_form_multimodular(10)
            [                  1                   0  -192/97*z - 361/97]
            [                  0                   1 1851/97*z + 1272/97]

        TESTS:

        We test a degenerate case::

            sage: A = matrix(CyclotomicField(5),0); A
            []
            sage: A._echelon_form_multimodular(10)
            []
            sage: A.pivots()
            ()

            sage: A = matrix(CyclotomicField(13), 2, 3, [5, 1, 2, 46307, 46307*4, 46307])
            sage: A._echelon_form_multimodular()
            [   1    0 7/19]
            [   0    1 3/19]
        """
        cdef Matrix_cyclo_dense res
        cdef bint is_square

        verbose("entering _echelon_form_multimodular",
                level=echelon_verbose_level)

        denom = self._matrix.denominator()
        A = denom * self

        # This bound is chosen somewhat arbitrarily. Changing it affects the
        # runtime, not the correctness of the result.
        if height_guess is None:
            height_guess = (A.coefficient_bound() + 100) * 1000000

        # This is all setup to keep track of various data
        # in the loop below.
        p = previous_prime(MAX_MODULUS)
        found = 0
        prod = 1
        n = self._base_ring._n()
        height_bound = self._ncols * height_guess * A.coefficient_bound() + 1
        mod_p_ech_ls = []
        max_pivots = tuple()
        is_square = self._nrows == self._ncols

        verbose("using height bound %s" % height_bound,
                level=echelon_verbose_level)

        while True:
            # Generate primes to use, and find echelon form
            # modulo those primes.
            while found < num_primes or prod <= height_bound:
                if (n == 1) or p % n == 1:
                    try:
                        mod_p_ech, piv_ls = A._echelon_form_one_prime(p)
                    except ValueError:
                        # This means that we chose a prime which divides
                        # the denominator of the echelon form of self, so
                        # just skip it and continue
                        p = previous_prime(p)
                        continue
                    # if we have the identity, just return it, and
                    # we're done.
                    if is_square and len(piv_ls) == self._nrows:
                        self.cache('pivots', tuple(range(self._nrows)))
                        return self.parent().identity_matrix()
                    if piv_ls > max_pivots:
                        mod_p_ech_ls = [mod_p_ech]
                        max_pivots = piv_ls
                        # add this to the list of primes
                        prod = p
                        found = 1
                    elif piv_ls == max_pivots:
                        mod_p_ech_ls.append(mod_p_ech)
                        # add this to the list of primes
                        prod *= p
                        found += 1
                    else:
                        # this means that the rank profile mod this
                        # prime is worse than those that came before,
                        # so we just loop
                        p = previous_prime(p)
                        continue

                p = previous_prime(p)

            if found > num_primes:
                num_primes = found

            verbose("computed echelon form mod %s primes" % num_primes,
                    level=echelon_verbose_level)
            verbose("current product of primes used: %s" % prod,
                    level=echelon_verbose_level)

            # Use CRT to lift back to ZZ
            mat_over_ZZ = matrix(ZZ, self._base_ring.degree(),
                                 self._nrows * self._ncols)
            _lift_crt(mat_over_ZZ, mod_p_ech_ls)
            # note: saving the CRT intermediate MultiModularBasis does
            # not seem to affect the runtime at all

            # Attempt to use rational reconstruction to find
            # our echelon form
            try:
                verbose("attempting rational reconstruction ...",
                        level=echelon_verbose_level)
                res = Matrix_cyclo_dense.__new__(Matrix_cyclo_dense,
                                                 self.parent(),
                                                 None, None, None)
                res._matrix = <Matrix_rational_dense>matrix_integer_dense_rational_reconstruction(mat_over_ZZ, prod)

            except ValueError:
                # If a ValueError is raised here, it means that the
                # rational reconstruction failed. In this case, add
                # on a few more primes, and try again.

                num_primes += echelon_primes_increment
                verbose("rational reconstruction failed, trying with %s primes"%num_primes, level=echelon_verbose_level)
                continue

            verbose("rational reconstruction succeeded with %s primes!"%num_primes, level=echelon_verbose_level)

            if ((res * res.denominator()).coefficient_bound() *
                self.coefficient_bound() * self.ncols()) > prod:
                # In this case, we don't know the result to sufficient
                # "precision" (here precision is just the modulus,
                # prod) to guarantee its correctness, so loop.

                num_primes += echelon_primes_increment
                verbose("height not sufficient to determine echelon form",
                        level=echelon_verbose_level)
                continue

            verbose("found echelon form with %s primes, whose product is %s"%(num_primes, prod), level=echelon_verbose_level)
            self.cache('pivots', max_pivots)
            return res

    def _echelon_form_one_prime(self, p):
        """
        Find the echelon form of ``self`` mod the primes dividing p. Return
        the rational matrix representing this lift. If the pivots of the
        reductions mod the primes over p are different, then no such lift
        exists, and we raise a :exc:`ValueError`. If this happens, then the
        denominator of the echelon form of ``self`` is divisible by p. (Note
        that the converse need not be true.)

        INPUT:

        - ``p`` -- a prime that splits completely in the cyclotomic base field

        OUTPUT: tuple of

        - ``matrix`` -- Lift via CRT of the echelon forms of ``self`` modulo
          each of the primes over p.
        - ``tuple`` -- the tuple of pivots for the echelon form of ``self`` mod the
          primes dividing p

        EXAMPLES::

            sage: W.<z> = CyclotomicField(3)
            sage: A = matrix(W, 2, 3, [1+z, 2/3, 9*z+7, -3 + 4*z, z, -7*z]); A
            [  z + 1     2/3 9*z + 7]
            [4*z - 3       z    -7*z]
            sage: A._echelon_form_one_prime(7)
            (
            [1 0 4 0 1 2]
            [0 0 3 0 0 4], (0, 1)
            )
            sage: Matrix(W,2,3,[2*z+3,0,1,0,1,0])._echelon_form_one_prime(7)
            Traceback (most recent call last):
            ...
            ValueError: echelon form mod 7 not defined
        """
        cdef int i

        # Initialize variables
        ls, _ = self._reductions(p)

        # Find our first echelon form, and the associated list
        # of pivots
        ech_ls = [ls[0].echelon_form()]
        pivot_ls = ech_ls[0].pivots()
        # If we've found the identity matrix, we're all done.
        if self._nrows == self._ncols == len(pivot_ls):
            return (self.parent().identity_matrix(), range(self._nrows))

        # For each reduction of self (i.e. for each prime of
        # self.base_ring() over p), compute the echelon form, and
        # keep track of all reductions which have the largest
        # number of pivots seen so far.
        for i from 1 <= i < len(ls):
            ech = ls[i].echelon_form()

            # This should only occur when p divides the denominator
            # of the echelon form of self.
            if ech.pivots() != pivot_ls:
                raise ValueError("echelon form mod %s not defined" % p)

            ech_ls.append(ech)

        # Now, just lift back to ZZ and return it.

        # TODO: coercion going on here
        reduction = matrix(ZZ, len(ech_ls), self._nrows * self._ncols,
                           [ [y.lift() for y in E.list()] for E in ech_ls])

        # TODO: more coercion happening here
        _, Finv = self._reduction_matrix(p)

        lifted_matrix = Finv * reduction

        return (lifted_matrix, pivot_ls)

    def tensor_product(self, A, subdivide=True):
        r"""
        Return the tensor product of two matrices.

        INPUT:

        - ``A`` -- a matrix
        - ``subdivide`` -- boolean (default: ``True``); whether or not to return
          natural subdivisions with the matrix

        OUTPUT:

        Replace each element of ``self`` by a copy of ``A``, but first
        create a scalar multiple of ``A`` by the element it replaces.
        So if ``self`` is an `m\times n` matrix and ``A`` is a
        `p\times q` matrix, then the tensor product is an `mp\times nq`
        matrix.  By default, the matrix will be subdivided into
        submatrices of size `p\times q`.

        EXAMPLES::

            sage: C = CyclotomicField(12)
            sage: M = matrix.random(C, 3, 3)
            sage: N = matrix.random(C, 50, 50)
            sage: M.tensor_product(M) == super(type(M), M).tensor_product(M)
            True
            sage: N = matrix.random(C, 15, 20)
            sage: M.tensor_product(N) == super(type(M), M).tensor_product(N)
            True

        TESTS::

            sage: Mp = matrix.random(C, 2,3)
            sage: Np = matrix.random(C, 4,5)
            sage: subdiv = super(type(Mp),Mp).tensor_product(Np).subdivisions()
            sage: Mp.tensor_product(Np).subdivisions() == subdiv
            True

        Check that `m \times 0` and `0 \times m` matrices work
        (:issue:`22769`)::

            sage: m1 = matrix(C, 1, 0, [])
            sage: m2 = matrix(C, 2, 2, [1, 2, 3, 4])
            sage: m1.tensor_product(m2).dimensions()
            (2, 0)
            sage: m2.tensor_product(m1).dimensions()
            (2, 0)
            sage: m3 = matrix(C, 0, 3, [])
            sage: m3.tensor_product(m2).dimensions()
            (0, 6)
            sage: m2.tensor_product(m3).dimensions()
            (0, 6)
        """
        if not isinstance(A, Matrix):
            raise TypeError('tensor product requires a second matrix, not {0}'.format(A))

        if A.base_ring() is not self.base_ring():
            return super(Matrix_cyclo_dense, self).tensor_product(A, subdivide)

        cdef Matrix_cyclo_dense M
        l = []
        R = self.base_ring()
        X = R._generator_matrix()
        d = self._degree
        MS = MatrixSpace(QQ, d, d)
        for c in self._matrix.columns():
            v = c.list()
            for n in range(d-1):
                c = c * X
                v += c.list()
            rmul = MS([v[d*i+j] for j in range(d) for i in range(d)]) # We take the transpose
            l.append(rmul * A._rational_matrix())

        nr = self.nrows()
        nc = self.ncols()
        Anr = A.nrows()
        Anc = A.ncols()
        P = MatrixSpace(R, nr*Anr, nc*Anc)
        M = Matrix_cyclo_dense.__new__(Matrix_cyclo_dense, P,
                                       None, None, None)
        MS = MatrixSpace(QQ, d, P.nrows()*P.ncols())
        ret = [[l[mr*nc+mc][i,r*Anc+c] for mr in range(nr) for r in range(Anr)
                for mc in range(nc) for c in range(Anc)]
               for i in range(d)]
        M._matrix = MS(ret)
        if subdivide:
            M.subdivide([Anr*i for i in range(1,nr)], [Anc*i for i in range(1,nc)])
        return M
