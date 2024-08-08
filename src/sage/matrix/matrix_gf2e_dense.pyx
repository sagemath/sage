# distutils: libraries = m4rie M4RI_LIBRARIES m
# distutils: library_dirs = M4RI_LIBDIR
# distutils: include_dirs = M4RI_INCDIR
# distutils: extra_compile_args = M4RI_CFLAGS
r"""
Dense matrices over `\GF{2^e}` for `2 \leq e \leq 16` using the M4RIE library

The M4RIE library offers two matrix representations:

1) ``mzed_t``

  m x n matrices over `\GF{2^e}` are internally represented roughly as
  m x (en) matrices over `\GF{2}`. Several elements are packed into
  words such that each element is filled with zeroes until the next
  power of two. Thus, for example, elements of `\GF{2^3}` are
  represented as ``[0xxx|0xxx|0xxx|0xxx|...]``. This representation is
  wrapped as :class:`Matrix_gf2e_dense` in Sage.

  Multiplication and elimination both use "Newton-John" tables. These
  tables are simply all possible multiples of a given row in a matrix
  such that a scale+add operation is reduced to a table lookup +
  add. On top of Newton-John multiplication M4RIE implements
  asymptotically fast Strassen-Winograd multiplication. Elimination
  uses simple Gaussian elimination which requires `O(n^3)` additions
  but only `O(n^2 * 2^e)` multiplications.

2) ``mzd_slice_t``

  m x n matrices over `\GF{2^e}` are internally represented as slices
  of m x n matrices over `\GF{2}`. This representation allows for very
  fast matrix times matrix products using Karatsuba's polynomial
  multiplication for polynomials over matrices. However, it is not
  feature complete yet and hence not wrapped in Sage for now.

See http://m4ri.sagemath.org for more details on the M4RIE library.

EXAMPLES::

    sage: K.<a> = GF(2^8)
    sage: A = random_matrix(K, 3,4)
    sage: E = A.echelon_form()
    sage: A.row_space() == E.row_space()
    True
    sage: all(r[r.nonzero_positions()[0]] == 1 for r in E.rows() if r)
    True

AUTHOR:

* Martin Albrecht <martinralbrecht@googlemail.com>

TESTS::

    sage: TestSuite(sage.matrix.matrix_gf2e_dense.Matrix_gf2e_dense).run(verbose=True)
    running ._test_new() . . . pass
    running ._test_pickling() . . . pass

Test hashing::

    sage: K.<a> = GF(2^4)
    sage: A = random_matrix(K, 1000, 1000)
    sage: A.set_immutable()
    sage: {A:1}
    {1000 x 1000 dense matrix over Finite Field in a of size 2^4: 1}

.. TODO::

    Wrap ``mzd_slice_t``.

REFERENCES:

- [BB2009]_
"""

#*****************************************************************************
#       Copyright (C) 2010 Martin Albrecht <martinralbrecht@googlemail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from cysignals.signals cimport sig_on, sig_off

cimport sage.matrix.matrix_dense as matrix_dense
from sage.structure.element cimport Matrix
from sage.structure.element cimport Element
from sage.structure.richcmp cimport rich_to_bool
from sage.rings.finite_rings.element_base cimport Cache_base

from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
from sage.misc.randstate cimport randstate, current_randstate

from sage.matrix.matrix_mod2_dense cimport Matrix_mod2_dense
from sage.matrix.args cimport SparseEntry, MatrixArgs_init

from sage.libs.m4ri cimport m4ri_word, mzd_copy
from sage.libs.m4rie cimport *
from sage.libs.m4rie cimport mzed_t


# we must keep a copy of the internal finite field representation
# around to avoid re-creating it over and over again. Furthermore,
# M4RIE assumes pointer equivalence of identical fields.

_m4rie_finite_field_cache = {}

cdef class M4RIE_finite_field:
    """
    A thin wrapper around the M4RIE finite field class such that we
    can put it in a hash table. This class is not meant for public
    consumption.
    """
    cdef gf2e *ff

    def __dealloc__(self):
        """
        EXAMPLES::

            sage: from sage.matrix.matrix_gf2e_dense import M4RIE_finite_field
            sage: K = M4RIE_finite_field(); K
            <sage.matrix.matrix_gf2e_dense.M4RIE_finite_field object at 0x...>
            sage: del K
        """
        if self.ff:
            gf2e_free(self.ff)

cdef m4ri_word poly_to_word(f) noexcept:
    return f.to_integer()


cdef class Matrix_gf2e_dense(matrix_dense.Matrix_dense):
    def __cinit__(self, *args, bint alloc=True, **kwds):
        """
        INPUT:

        - ``alloc`` -- if ``True`` the matrix is allocated first (default: ``True``)

        EXAMPLES::

            sage: K.<a> = GF(2^4)
            sage: A = Matrix(K, 3, 4); A
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]

            sage: A.randomize()
            sage: TestSuite(A).run()

            sage: K.<a> = GF(2^3)
            sage: A = Matrix(K,3,4); A
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]

            sage: A.randomize()
            sage: TestSuite(A).run()
        """
        cdef M4RIE_finite_field FF

        R = self._base_ring
        f = R.polynomial()

        cdef long i
        cdef m4ri_word poly = sum(((<m4ri_word>c) << i) for (i, c) in enumerate(f))

        if alloc and self._nrows and self._ncols:
            if poly in _m4rie_finite_field_cache:
                self._entries = mzed_init((<M4RIE_finite_field>_m4rie_finite_field_cache[poly]).ff, self._nrows, self._ncols)
            else:
                FF = M4RIE_finite_field.__new__(M4RIE_finite_field)
                FF.ff = gf2e_init(poly)
                self._entries = mzed_init(FF.ff, self._nrows, self._ncols)
                _m4rie_finite_field_cache[poly] = FF

        # cache elements
        self._zero = self._base_ring(0)
        self._zero_word = poly_to_word(self._zero)
        self._one = self._base_ring(1)

    def __dealloc__(self):
        """
        TESTS::

            sage: K.<a> = GF(2^4)
            sage: A = Matrix(K, 1000, 1000)
            sage: del A
            sage: A = Matrix(K, 1000, 1000)
            sage: del A
        """
        if self._entries:
            mzed_free(self._entries)
            self._entries = NULL

    def __init__(self, parent, entries=None, copy=None, bint coerce=True):
        r"""
        Create new matrix over `GF(2^e)` for `2 \leq e \leq 16`.

        INPUT:

        - ``parent`` -- a matrix space over ``GF(2^e)``

        - ``entries`` -- see :func:`matrix`

        - ``copy`` -- ignored (for backwards compatibility)

        - ``coerce`` -- if ``False``, assume without checking that the
          entries lie in the base ring

        EXAMPLES::

            sage: K.<a> = GF(2^4)
            sage: l = [K.random_element() for _ in range(3*4)]

            sage: A = Matrix(K, 3, 4, l)
            sage: l == A.list()
            True

            sage: l[0] == A[0,0]
            True

            sage: A = Matrix(K, 3, 3, a); A
            [a 0 0]
            [0 a 0]
            [0 0 a]
        """
        ma = MatrixArgs_init(parent, entries)
        for t in ma.iter(coerce, True):
            se = <SparseEntry>t
            mzed_write_elem(self._entries, se.i, se.j, poly_to_word(se.entry))

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, value):
        """
        A[i,j] = value without bound checks.

        INPUT:

        - ``i`` -- row index
        - ``j`` -- column index
        - ``value`` -- a finite field element (not checked but assumed)

        EXAMPLES::

            sage: K.<a> = GF(2^4)
            sage: l = [K.random_element() for _ in range(3*4)]
            sage: A = Matrix(K, 3, 4, l)

            sage: i = randrange(3)
            sage: j = randrange(4)
            sage: A[i,j] = a                     # indirect doctest
            sage: A[i,j] == a == A.list()[j + 4*i]
            True
            sage: A.list()[:j + 4*i] == l[:j + 4*i]
            True
            sage: A.list()[j + 4*i + 1:] == l[j + 4*i + 1:]
            True
        """
        mzed_write_elem(self._entries, i, j, poly_to_word(value))

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        """
        Get A[i,j] without bound checks.

        INPUT:

        - ``i`` -- row index
        - ``j`` -- column index

        EXAMPLES::

            sage: K.<a> = GF(2^4)
            sage: l = [K.random_element() for _ in range(3*4)]
            sage: A = Matrix(K, 3, 4, l)
            sage: A[2,3] == l[3 + 4*2]           # indirect doctest
            True
            sage: K.<a> = GF(2^3)
            sage: l = [K.random_element() for _ in range(3*4)]
            sage: A = Matrix(K, 3, 4, l)
            sage: A.list() == l
            True
        """
        cdef int r = mzed_read_elem(self._entries, i, j)
        cdef Cache_base cache = <Cache_base> self._base_ring._cache
        return cache.fetch_int(r)

    cdef bint get_is_zero_unsafe(self, Py_ssize_t i, Py_ssize_t j) except -1:
        r"""
        Return 1 if the entry ``(i, j)`` is zero, otherwise 0.

        EXAMPLES::

            sage: K.<a> = GF(2^4)
            sage: A = Matrix(K, 2, 2, a)
            sage: A.zero_pattern_matrix()  # indirect doctest
            [0 1]
            [1 0]
        """
        return mzed_read_elem(self._entries, i, j) == self._zero_word

    cpdef _add_(self, right):
        r"""
        Return ``A+B``.

        INPUT:

        - ``right`` -- a matrix

        EXAMPLES::

            sage: K.<a> = GF(2^4)
            sage: A = random_matrix(K,3,4)
            sage: B = random_matrix(K,3,4)

            sage: C = A + B  # indirect doctest
            sage: all(C.list()[i] == A.list()[i] + B.list()[i] for i in range(12))
            True
        """
        cdef Matrix_gf2e_dense A
        A = Matrix_gf2e_dense.__new__(Matrix_gf2e_dense, self._parent, 0, 0, 0, alloc=False)
        if self._nrows == 0 or self._ncols == 0:
            return A
        A._entries = mzed_add(NULL, self._entries, (<Matrix_gf2e_dense>right)._entries)

        return A

    cpdef _sub_(self, right):
        """
        EXAMPLES::

            sage: from sage.matrix.matrix_gf2e_dense import Matrix_gf2e_dense
            sage: K.<a> = GF(2^4)
            sage: A = random_matrix(K,3,4)
            sage: B = random_matrix(K,3,4)
            sage: C = A - B  # indirect doctest
            sage: all(C.list()[i] == A.list()[i] - B.list()[i] for i in range(12))
            True
        """
        return self._add_(right)

    def _multiply_classical(self, Matrix right):
        """
        Classical cubic matrix multiplication.

        EXAMPLES::

            sage: K.<a> = GF(2^2)
            sage: A = random_matrix(K, 50, 50)
            sage: B = random_matrix(K, 50, 50)
            sage: A*B == A._multiply_classical(B)
            True

            sage: K.<a> = GF(2^4)
            sage: A = random_matrix(K, 50, 50)
            sage: B = random_matrix(K, 50, 50)
            sage: A*B == A._multiply_classical(B)
            True

            sage: K.<a> = GF(2^8)
            sage: A = random_matrix(K, 50, 50)
            sage: B = random_matrix(K, 50, 50)
            sage: A*B == A._multiply_classical(B)
            True

            sage: K.<a> = GF(2^10)
            sage: A = random_matrix(K, 50, 50)
            sage: B = random_matrix(K, 50, 50)
            sage: A*B == A._multiply_classical(B)
            True

        .. NOTE::

            This function is very slow. Use ``*`` operator instead.
        """
        if self._ncols != right._nrows:
            raise ArithmeticError("left ncols must match right nrows")

        cdef Matrix_gf2e_dense ans

        ans = self.new_matrix(nrows = self.nrows(), ncols = right.ncols())
        if self._nrows == 0 or self._ncols == 0 or right._ncols == 0:
            return ans
        sig_on()
        ans._entries = mzed_mul_naive(ans._entries, self._entries, (<Matrix_gf2e_dense>right)._entries)
        sig_off()
        return ans

    cdef _matrix_times_matrix_(self, Matrix right):
        r"""
        Return ``A*B``.

        Uses the M4RIE machinery to decide which function to call.

        INPUT:

        - ``right`` -- a matrix

        EXAMPLES::

            sage: K.<a> = GF(2^3)
            sage: A = random_matrix(K, 51, 50)
            sage: B = random_matrix(K, 50, 52)
            sage: A*B == A._multiply_newton_john(B) # indirect doctest
            True

            sage: K.<a> = GF(2^5)
            sage: A = random_matrix(K, 10, 50)
            sage: B = random_matrix(K, 50, 12)
            sage: A*B == A._multiply_newton_john(B)
            True

            sage: K.<a> = GF(2^7)
            sage: A = random_matrix(K,100, 50)
            sage: B = random_matrix(K, 50, 17)
            sage: A*B == A._multiply_classical(B)
            True
        """
        if self._ncols != right._nrows:
            raise ArithmeticError("left ncols must match right nrows")

        cdef Matrix_gf2e_dense ans

        ans = self.new_matrix(nrows = self.nrows(), ncols = right.ncols())
        if self._nrows == 0 or self._ncols == 0 or right._ncols == 0:
            return ans
        sig_on()
        ans._entries = mzed_mul(ans._entries, self._entries, (<Matrix_gf2e_dense>right)._entries)
        sig_off()
        return ans

    cpdef Matrix_gf2e_dense _multiply_newton_john(Matrix_gf2e_dense self, Matrix_gf2e_dense right):
        """
        Return A*B using Newton-John tables.

        We can write classical cubic multiplication (``C=A*B``) as::

        for i in range(A.ncols()):
           for j in range(A.nrows()):
             C[j] += A[j,i] * B[j]

        Hence, in the inner-most loop we compute multiples of ``B[j]``
        by the values ``A[j,i]``. If the matrix ``A`` is big and the
        finite field is small, there is a very high chance that
        ``e * B[j]`` is computed more than once for any ``e`` in the finite
        field. Instead, we compute all possible
        multiples of ``B[j]`` and re-use this data in the inner loop.
        This is what is called a "Newton-John" table in M4RIE.

        INPUT:

        - ``right`` -- a matrix

        EXAMPLES::

            sage: K.<a> = GF(2^2)
            sage: A = random_matrix(K, 50, 50)
            sage: B = random_matrix(K, 50, 50)
            sage: A._multiply_newton_john(B) == A._multiply_classical(B) # indirect doctest
            True

            sage: K.<a> = GF(2^4)
            sage: A = random_matrix(K, 50, 50)
            sage: B = random_matrix(K, 50, 50)
            sage: A._multiply_newton_john(B) == A._multiply_classical(B)
            True

            sage: K.<a> = GF(2^8)
            sage: A = random_matrix(K, 50, 50)
            sage: B = random_matrix(K, 50, 50)
            sage: A._multiply_newton_john(B) == A._multiply_classical(B)
            True

            sage: K.<a> = GF(2^10)
            sage: A = random_matrix(K, 50, 50)
            sage: B = random_matrix(K, 50, 50)
            sage: A._multiply_newton_john(B) == A._multiply_classical(B)
            True
        """
        if self._ncols != right._nrows:
            raise ArithmeticError("left ncols must match right nrows")

        cdef Matrix_gf2e_dense ans

        ans = self.new_matrix(nrows = self.nrows(), ncols = right.ncols())
        if self._nrows == 0 or self._ncols == 0 or right._ncols == 0:
            return ans

        sig_on()
        ans._entries = mzed_mul_newton_john(ans._entries, self._entries, (<Matrix_gf2e_dense>right)._entries)
        sig_off()
        return ans

    cpdef Matrix_gf2e_dense _multiply_karatsuba(Matrix_gf2e_dense self, Matrix_gf2e_dense right):
        r"""
        Matrix multiplication using Karatsuba over polynomials with
        matrix coefficients over GF(2).

        The idea behind Karatsuba multiplication for matrices over
        `\GF{p^n}` is to treat these matrices as polynomials with
        coefficients of matrices over `\GF{p}`. Then, Karatsuba-style
        formulas can be used to perform multiplication, cf. [BB2009]_.

        INPUT:

        - ``right`` -- a matrix

        EXAMPLES::

            sage: K.<a> = GF(2^2)
            sage: A = random_matrix(K, 50, 50)
            sage: B = random_matrix(K, 50, 50)
            sage: A._multiply_karatsuba(B) == A._multiply_classical(B) # indirect doctest
            True

            sage: K.<a> = GF(2^2)
            sage: A = random_matrix(K, 137, 11)
            sage: B = random_matrix(K, 11, 23)
            sage: A._multiply_karatsuba(B) == A._multiply_classical(B)
            True

            sage: K.<a> = GF(2^10)
            sage: A = random_matrix(K, 50, 50)
            sage: B = random_matrix(K, 50, 50)
            sage: A._multiply_karatsuba(B) == A._multiply_classical(B)
            True
        """
        if self._ncols != right._nrows:
            raise ArithmeticError("left ncols must match right nrows")

        cdef Matrix_gf2e_dense ans

        ans = self.new_matrix(nrows = self.nrows(), ncols = right.ncols())
        if self._nrows == 0 or self._ncols == 0 or right._ncols == 0:
            return ans

        sig_on()
        ans._entries = mzed_mul_karatsuba(ans._entries, self._entries, (<Matrix_gf2e_dense>right)._entries)
        sig_off()
        return ans

    cpdef Matrix_gf2e_dense _multiply_strassen(Matrix_gf2e_dense self, Matrix_gf2e_dense right, cutoff=0):
        """
        Winograd-Strassen matrix multiplication with Newton-John
        multiplication as base case.

        INPUT:

        - ``right`` -- a matrix
        - ``cutoff`` -- row or column dimension to switch over to
          Newton-John multiplication (default: 64)

        EXAMPLES::

            sage: K.<a> = GF(2^2)
            sage: A = random_matrix(K, 50, 50)
            sage: B = random_matrix(K, 50, 50)
            sage: A._multiply_strassen(B) == A._multiply_classical(B) # indirect doctest
            True

            sage: K.<a> = GF(2^4)
            sage: A = random_matrix(K, 50, 50)
            sage: B = random_matrix(K, 50, 50)
            sage: A._multiply_strassen(B) == A._multiply_classical(B)
            True

            sage: K.<a> = GF(2^8)
            sage: A = random_matrix(K, 50, 50)
            sage: B = random_matrix(K, 50, 50)
            sage: A._multiply_strassen(B) == A._multiply_classical(B)
            True

            sage: K.<a> = GF(2^10)
            sage: A = random_matrix(K, 50, 50)
            sage: B = random_matrix(K, 50, 50)
            sage: A._multiply_strassen(B) == A._multiply_classical(B)
            True
        """
        if self._ncols != right._nrows:
            raise ArithmeticError("left ncols must match right nrows")

        cdef Matrix_gf2e_dense ans

        ans = self.new_matrix(nrows = self.nrows(), ncols = right.ncols())
        if self._nrows == 0 or self._ncols == 0 or right._ncols == 0:
            return ans

        if cutoff == 0:
            cutoff = _mzed_strassen_cutoff(ans._entries, self._entries, (<Matrix_gf2e_dense>right)._entries)

        sig_on()
        ans._entries = mzed_mul_strassen(ans._entries, self._entries, (<Matrix_gf2e_dense>right)._entries, cutoff)
        sig_off()
        return ans

    cpdef _lmul_(self, Element right):
        """
        Return ``a*B`` for ``a`` an element of the base field.

        INPUT:

        - ``right`` -- an element of the base field

        EXAMPLES::

             sage: K.<a> = GF(4)
             sage: A = random_matrix(K,10,10)
             sage: B = a*A  # indirect doctest
             sage: all(B.list()[i] == a*A.list()[i] for i in range(100))
             True
        """
        cdef m4ri_word a = poly_to_word(right)
        cdef Matrix_gf2e_dense C = Matrix_gf2e_dense.__new__(Matrix_gf2e_dense, self._parent, 0, 0, 0)
        mzed_mul_scalar(C._entries, a, self._entries)
        return C

    def __neg__(self):
        """
        EXAMPLES::

            sage: K.<a> = GF(2^4)
            sage: A = random_matrix(K, 3, 4)
            sage: B = -A
            sage: all(B.list()[i] == -A.list()[i] for i in range(12))
            True
        """
        return self.__copy__()

    cpdef _richcmp_(self, right, int op):
        """
        EXAMPLES::

            sage: K.<a> = GF(2^4)
            sage: A = random_matrix(K,3,4)
            sage: A[0,0] = 0
            sage: B = copy(A)
            sage: A == B
            True
            sage: A[0,0] = a
            sage: A == B
            False
        """
        if self._nrows == 0 or self._ncols == 0:
            return rich_to_bool(op, 0)
        return rich_to_bool(op, mzed_cmp(self._entries,
                                         (<Matrix_gf2e_dense>right)._entries))

    def __copy__(self):
        """
        EXAMPLES::

            sage: K.<a> = GF(2^4)
            sage: m,n  = 3, 4
            sage: A = random_matrix(K,3,4)
            sage: A2 = copy(A)
            sage: A2.list() == A.list()
            True

            sage: A[0,0] = 1 if A[0,0] != 1 else 0
            sage: A2[0,0] == A[0,0]
            False
        """
        cdef Matrix_gf2e_dense A
        A = Matrix_gf2e_dense.__new__(Matrix_gf2e_dense, self._parent, 0, 0, 0)

        if self._nrows and self._ncols:
            mzed_copy(A._entries, <const_mzed_t *>self._entries)

        return A

    def __bool__(self):
        """
        Return if ``self`` is a zero matrix or not.

        EXAMPLES::

            sage: K.<a> = GF(2^4)
            sage: A = Matrix(K, 2, 2, a)
            sage: bool(A)
            True
            sage: zero = MatrixSpace(K, 3, 3).zero()
            sage: bool(zero)
            False
        """
        if self._nrows and self._ncols:
            return not mzed_is_zero(self._entries)
        return False

    def _list(self):
        """
        EXAMPLES::

            sage: K.<a> = GF(2^4)
            sage: l = [K.random_element() for _ in range(3*4)]
            sage: A = Matrix(K, 3, 4, l)
            sage: A.list() == l  # indirect doctest
            True
        """
        cdef int i,j
        l = []
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                l.append(self.get_unsafe(i,j))
        return l

    def randomize(self, density=1, nonzero=False, *args, **kwds):
        """
        Randomize ``density`` proportion of the entries of this matrix,
        leaving the rest unchanged.

        INPUT:

        - ``density`` -- float; proportion (roughly) to be considered for
          changes
        - ``nonzero`` -- boolean (default: ``False``); whether the new entries
          are forced to be nonzero

        OUTPUT: none, the matrix is modified in-place

        EXAMPLES::

            sage: K.<a> = GF(2^4)
            sage: total_count = 0
            sage: from collections import defaultdict
            sage: dic = defaultdict(Integer)
            sage: def add_samples():
            ....:     global dic, total_count
            ....:     for _ in range(100):
            ....:         A = Matrix(K,3,3)
            ....:         A.randomize()
            ....:         for a in A.list():
            ....:             dic[a] += 1
            ....:             total_count += 1.0
            sage: add_samples()
            sage: while not all(abs(dic[a]/total_count - 1/16) < 0.01 for a in dic):
            ....:     add_samples()

            sage: def add_sample(density):
            ....:     global density_sum, total_count
            ....:     total_count += 1.0
            ....:     density_sum += random_matrix(K, 1000, 1000, density=density).density()

            sage: density_sum = 0.0
            sage: total_count = 0.0
            sage: add_sample(0.1)
            sage: while abs(density_sum/total_count - 0.1) > 0.001:
            ....:     add_sample(0.1)

            sage: density_sum = 0.0
            sage: total_count = 0.0
            sage: add_sample(1.0)
            sage: while abs(density_sum/total_count - 1.0) > 0.001:
            ....:     add_sample(1.0)

            sage: density_sum = 0.0
            sage: total_count = 0.0
            sage: add_sample(0.5)
            sage: while abs(density_sum/total_count - 0.5) > 0.001:
            ....:     add_sample(0.5)

        Note, that the matrix is updated and not zero-ed out before
        being randomized::

            sage: def add_sample(density, nonzero):
            ....:     global density_sum, total_count
            ....:     total_count += 1.0
            ....:     A = matrix(K, 1000, 1000)
            ....:     A.randomize(nonzero=nonzero, density=density)
            ....:     A.randomize(nonzero=nonzero, density=density)
            ....:     density_sum += A.density()

            sage: density_sum = 0.0
            sage: total_count = 0.0
            sage: add_sample(0.1, True)
            sage: while abs(density_sum/total_count - (1 - 0.9^2)) > 0.001:
            ....:     add_sample(0.1, True)

            sage: density_sum = 0.0
            sage: total_count = 0.0
            sage: add_sample(0.1, False)
            sage: while abs(density_sum/total_count - (1 - 0.9^2)*15/16) > 0.001:
            ....:     add_sample(0.1, False)

            sage: density_sum = 0.0
            sage: total_count = 0.0
            sage: add_sample(0.05, True)
            sage: while abs(density_sum/total_count - (1 - 0.95^2)) > 0.001:
            ....:     add_sample(0.05, True)

            sage: density_sum = 0.0
            sage: total_count = 0.0
            sage: add_sample(0.5, True)
            sage: while abs(density_sum/total_count - (1 - 0.5^2)) > 0.001:
            ....:     add_sample(0.5, True)
        """
        if self._ncols == 0 or self._nrows == 0:
            return

        cdef Py_ssize_t i,j
        self.check_mutability()
        self.clear_cache()

        cdef m4ri_word mask = (1<<(self._parent.base_ring().degree())) - 1

        cdef randstate rstate = current_randstate()

        if self._ncols == 0 or self._nrows == 0:
            return

        cdef float _density = density
        if _density <= 0:
            return
        if _density > 1:
            _density = 1.0

        if _density == 1:
            if not nonzero:
                sig_on()
                for i in range(self._nrows):
                    for j in range(self._ncols):
                        tmp = rstate.c_random() & mask
                        mzed_write_elem(self._entries, i, j, tmp)
                sig_off()
            else:
                sig_on()
                for i in range(self._nrows):
                    for j in range(self._ncols):
                        tmp = rstate.c_random() & mask
                        while tmp == 0:
                            tmp = rstate.c_random() & mask
                        mzed_write_elem(self._entries, i, j, tmp)
                sig_off()
        else:
            if not nonzero:
                sig_on()
                for i in range(self._nrows):
                    for j in range(self._ncols):
                        if rstate.c_rand_double() <= _density:
                            tmp = rstate.c_random() & mask
                            mzed_write_elem(self._entries, i, j, tmp)
                sig_off()
            else:
                sig_on()
                for i in range(self._nrows):
                    for j in range(self._ncols):
                        if rstate.c_rand_double() <= _density:
                            tmp = rstate.c_random() & mask
                            while tmp == 0:
                                tmp = rstate.c_random() & mask
                            mzed_write_elem(self._entries, i, j, tmp)
                sig_off()

    def echelonize(self, algorithm='heuristic', reduced=True, **kwds):
        """
        Compute the row echelon form of ``self`` in place.

        INPUT:

        - ``algorithm`` -- one of the following
          - ``heuristic`` -- let M4RIE decide (default)
          - ``newton_john`` -- use newton_john table based algorithm
          - ``ple`` -- use PLE decomposition
          - ``naive`` -- use naive cubic Gaussian elimination (M4RIE implementation)
          - ``builtin`` -- use naive cubic Gaussian elimination (Sage implementation)
        - ``reduced`` -- if ``True`` return reduced echelon form. No
          guarantee is given that the matrix is *not* reduced if
          ``False`` (default: ``True``)

        EXAMPLES::

            sage: K.<a> = GF(2^4)
            sage: m,n  = 3, 5
            sage: A = random_matrix(K, 3, 5)
            sage: R = A.row_space()
            sage: A.echelonize()
            sage: all(r[r.nonzero_positions()[0]] == 1 for r in A.rows() if r)
            True
            sage: A.row_space() == R
            True

            sage: K.<a> = GF(2^3)
            sage: m,n  = 3, 5
            sage: MS = MatrixSpace(K,m,n)
            sage: A = random_matrix(K, 3, 5)
            sage: B = copy(A).echelon_form('newton_john')
            sage: C = copy(A).echelon_form('naive')
            sage: D = copy(A).echelon_form('builtin')
            sage: B == C == D
            True
            sage: all(r[r.nonzero_positions()[0]] == 1 for r in B.rows() if r)
            True
        """
        if self._nrows == 0 or self._ncols == 0:
            self.cache('in_echelon_form',True)
            self.cache('rank', 0)
            self.cache('pivots', [])
            return self

        cdef int full

        full = int(reduced)

        x = self.fetch('in_echelon_form')
        if x is not None:
            return  # already known to be in echelon form

        self.check_mutability()
        self.clear_cache()

        if algorithm == 'naive':
            sig_on()
            r =  mzed_echelonize_naive(self._entries, full)
            sig_off()

        elif algorithm == 'newton_john':
            sig_on()
            r =  mzed_echelonize_newton_john(self._entries, full)
            sig_off()

        elif algorithm == 'ple':
            sig_on()
            r =  mzed_echelonize_ple(self._entries, full)
            sig_off()

        elif algorithm == 'heuristic':
            sig_on()
            r =  mzed_echelonize(self._entries, full)
            sig_off()

        elif algorithm == 'builtin':
            self._echelon_in_place(algorithm='classical')

        else:
            raise ValueError("No algorithm '%s'."%algorithm)

        self.cache('in_echelon_form',True)
        self.cache('rank', r)
        self.cache('pivots', self._pivots())

    def _pivots(self):
        """
        EXAMPLES::

            sage: K.<a> = GF(2^8)
            sage: A = random_matrix(K, 15, 15)
            sage: while A.rank() != 15:
            ....:     A = random_matrix(K, 15, 15)
            sage: A.pivots()  # indirect doctest
            (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)
        """
        if not self.fetch('in_echelon_form'):
            raise ValueError("self must be in reduced row echelon form first.")
        pivots = []
        cdef Py_ssize_t i, j, nc
        nc = self._ncols
        i = 0
        while i < self._nrows:
            for j from i <= j < nc:
                if not self.get_is_zero_unsafe(i,j):
                    pivots.append(j)
                    i += 1
                    break
            if j == nc:
                break
        return pivots

    def __invert__(self):
        """
        EXAMPLES::

            sage: K.<a> = GF(2^3)
            sage: A = random_matrix(K, 3, 3)
            sage: while A.rank() != 3:
            ....:     A = random_matrix(K, 3, 3)
            sage: B = ~A
            sage: A*B
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        cdef Matrix_gf2e_dense A
        A = Matrix_gf2e_dense.__new__(Matrix_gf2e_dense, self._parent, 0, 0, 0)

        if self.rank() != self._nrows:
            raise ZeroDivisionError("Matrix does not have full rank.")

        if self._nrows:
            sig_on()
            mzed_invert_newton_john(A._entries, self._entries)
            sig_off()

        return A

    cdef rescale_row_c(self, Py_ssize_t row, multiple, Py_ssize_t start_col):
        """
        Return ``multiple * self[row][start_col:]``.

        INPUT:

        - ``row`` -- row index for row to rescale
        - ``multiple`` -- finite field element to scale by
        - ``start_col`` -- only start at this column index

        EXAMPLES::

            sage: K.<a> = GF(2^3)
            sage: A = random_matrix(K,3,3)
            sage: B = copy(A)
            sage: B.rescale_row(0, a, 0)
            sage: B[0] == a*A[0]
            True
            sage: B[1:] == A[1:]
            True

            sage: B = copy(A)
            sage: B.rescale_row(1, 0, 0)
            sage: B[0] == A[0]
            True
            sage: B[2] == A[2]
            True
            sage: B[1].is_zero()
            True

            sage: B = copy(A)
            sage: B.rescale_row(2, a^2, 2)
            sage: B.list()[:-1] == A.list()[:-1]
            True
            sage: B[2,2] == a^2*A[2,2]
            True
        """
        cdef m4ri_word x = poly_to_word(multiple)
        mzed_rescale_row(self._entries, row, start_col, x)

    cdef add_multiple_of_row_c(self, Py_ssize_t row_to, Py_ssize_t row_from, multiple, Py_ssize_t start_col):
        """
        Compute ``self[row_to][start_col:] += multiple * self[row_from][start_col:]``.

        INPUT:

        - ``row_to`` -- row index of source
        - ``row_from`` -- row index of destination
        - ``multiple`` -- finite field element
        - ``start_col`` -- only start at this column index

        EXAMPLES::

            sage: K.<a> = GF(2^3)
            sage: A = random_matrix(K,3,3)
            sage: B = copy(A)
            sage: B.add_multiple_of_row(0,1,a,0)
            sage: B[1:] == A[1:]
            True
            sage: B[0] == A[0] + a*A[1]
            True

            sage: B = copy(A)
            sage: B.add_multiple_of_row(2,1,a,2)
            sage: B.list()[:-1] == A.list()[:-1]
            True
            sage: B[2,2] == A[2,2] + a*A[1,2]
            True
        """

        cdef m4ri_word x = poly_to_word(multiple)
        mzed_add_multiple_of_row(self._entries, row_to, self._entries, row_from, x, start_col)

    cdef swap_rows_c(self, Py_ssize_t row1, Py_ssize_t row2):
        """
        Swap rows ``row1`` and ``row2``.

        INPUT:

        - ``row1`` -- row index
        - ``row2`` -- row index

        EXAMPLES::

            sage: K.<a> = GF(2^3)
            sage: A = random_matrix(K,3,3)
            sage: B = copy(A)
            sage: B.swap_rows(0,1)
            sage: B[0] == A[1]
            True
            sage: B[1] == A[0]
            True
            sage: B[2] == A[2]
            True
        """
        mzed_row_swap(self._entries, row1, row2)

    cdef swap_columns_c(self, Py_ssize_t col1, Py_ssize_t col2):
        """
        Swap columns ``col1`` and ``col2``.

        INPUT:

        - ``col1`` -- column index
        - ``col2`` -- column index

        EXAMPLES::

            sage: K.<a> = GF(2^3)
            sage: A = random_matrix(K,3,3)
            sage: B = copy(A)
            sage: B.swap_columns(0,1)
            sage: B.column(0) == A.column(1)
            True
            sage: B.column(1) == A.column(0)
            True
            sage: B.column(2) == A.column(2)
            True

            sage: A = random_matrix(K,4,16)
            sage: B = copy(A)
            sage: B.swap_columns(0,1)
            sage: B.swap_columns(0,1)
            sage: A == B
            True

            sage: A.swap_columns(0,15)
            sage: A.column(0) == B.column(15)
            True
            sage: A.swap_columns(14,15)
            sage: A.column(14) == B.column(0)
            True
        """
        mzed_col_swap(self._entries, col1, col2)

    def augment(self, Matrix_gf2e_dense right):
        """
        Augments ``self`` with ``right``.

        INPUT:

        - ``right`` -- a matrix

        EXAMPLES::

            sage: K.<a> = GF(2^4)
            sage: MS = MatrixSpace(K,3,3)
            sage: A = random_matrix(K,3,3)
            sage: B = A.augment(MS(1))
            sage: B.echelonize()
            sage: C = B.matrix_from_columns([3,4,5])
            sage: A.rank() < 3 or C == ~A
            True
            sage: A.rank() < 3 or C*A == MS(1)
            True

        TESTS::

            sage: K.<a> = GF(2^4)
            sage: A = random_matrix(K,2,3)
            sage: B = random_matrix(K,2,0)
            sage: A.augment(B) == A
            True

            sage: B.augment(A) == A
            True

            sage: M = Matrix(K, 0, 0, 0)
            sage: N = Matrix(K, 0, 19, 0)
            sage: W = M.augment(N)
            sage: W.ncols()
            19

            sage: M = Matrix(K, 0, 1, 0)
            sage: N = Matrix(K, 0, 1, 0)
            sage: M.augment(N)
            []
        """
        cdef Matrix_gf2e_dense A

        if self._nrows != right._nrows:
            raise TypeError("Both numbers of rows must match.")

        if self._ncols == 0:
            return right.__copy__()
        if right._ncols == 0:
            return self.__copy__()

        A = self.new_matrix(ncols = self._ncols + right._ncols)
        if self._nrows == 0:
            return A
        A._entries = mzed_concat(A._entries, self._entries, right._entries)
        return A

    cdef _stack_impl(self, bottom):
        r"""
        Stack ``self`` on top of ``bottom``.

        INPUT:

        - ``bottom`` -- a matrix with the same number of columns as ``self``

        EXAMPLES::

            sage: K.<a> = GF(2^4)
            sage: A = random_matrix(K,2,2)
            sage: B = random_matrix(K,2,2)
            sage: C = A.stack(B)
            sage: C[:2] == A
            True
            sage: C[2:] == B
            True
            sage: D = B.stack(A)
            sage: D[:2] == B
            True
            sage: D[2:] == A
            True

        TESTS::

            sage: A = random_matrix(K,0,3)
            sage: B = random_matrix(K,3,3)
            sage: A.stack(B) == B
            True
            sage: B.stack(A) == B
            True

            sage: M = Matrix(K, 0, 0, 0)
            sage: N = Matrix(K, 19, 0, 0)
            sage: W = M.stack(N)
            sage: W.nrows()
            19
            sage: M = Matrix(K, 1, 0, 0)
            sage: N = Matrix(K, 1, 0, 0)
            sage: M.stack(N)
            []

        Check that we can stack a vector (:issue:`31708`)::

            sage: R.<a> = GF(2^3)
            sage: M = matrix(R, [[1,1],[0,a+1]])
            sage: M.stack(vector(R, [a,0]))
            [    1     1]
            [    0 a + 1]
            [    a     0]
        """
        cdef Matrix_gf2e_dense other = <Matrix_gf2e_dense> bottom

        if self._nrows == 0:
            return other.__copy__()
        if other._nrows == 0:
            return self.__copy__()

        cdef Matrix_gf2e_dense A
        A = self.new_matrix(nrows=self._nrows + other._nrows)
        if self._ncols == 0:
            return A
        A._entries = mzed_stack(A._entries, self._entries, other._entries)
        return A

    def submatrix(self, Py_ssize_t row=0, Py_ssize_t col=0,
                        Py_ssize_t nrows=-1, Py_ssize_t ncols=-1):
        """
        Return submatrix from the index ``row,col`` (inclusive) with
        dimension ``nrows x ncols``.

        INPUT:

        - ``row`` -- index of start row
        - ``col`` -- index of start column
        - ``nrows`` -- number of rows of submatrix
        - ``ncols`` -- number of columns of submatrix

        EXAMPLES::

             sage: K.<a> = GF(2^10)
             sage: A = random_matrix(K,200,200)
             sage: A[0:2,0:2] == A.submatrix(0,0,2,2)
             True
             sage: A[0:100,0:100] == A.submatrix(0,0,100,100)
             True
             sage: A == A.submatrix(0,0,200,200)
             True

             sage: A[1:3,1:3] == A.submatrix(1,1,2,2)
             True
             sage: A[1:100,1:100] == A.submatrix(1,1,99,99)
             True
             sage: A[1:200,1:200] == A.submatrix(1,1,199,199)
             True

        TESTS for handling of default arguments (:issue:`18761`)::

             sage: A.submatrix(17,15) == A.submatrix(17,15,183,185)
             True
             sage: A.submatrix(row=100,col=37,nrows=1,ncols=3) == A.submatrix(100,37,1,3)
             True
        """
        if nrows < 0:
            nrows = self._nrows - row

        if ncols < 0:
            ncols = self._ncols - col

        cdef int highr = row + nrows
        cdef int highc = col + ncols

        if row < 0:
            raise TypeError("Expected row >= 0, but got %d instead."%row)

        if col < 0:
            raise TypeError("Expected col >= 0, but got %d instead."%col)

        if highc > self._entries.ncols:
            raise TypeError("Expected highc <= self.ncols(), but got %d > %d instead."%(highc, self._entries.ncols))

        if highr > self._entries.nrows:
            raise TypeError("Expected highr <= self.nrows(), but got %d > %d instead."%(highr, self._entries.nrows))

        cdef Matrix_gf2e_dense A = self.new_matrix(nrows = nrows, ncols = ncols)
        if ncols == 0 or nrows == 0:
            return A
        A._entries = mzed_submatrix(A._entries, self._entries, row, col, highr, highc)
        return A

    def rank(self):
        """
        Return the rank of this matrix (cached).

        EXAMPLES::

            sage: K.<a> = GF(2^4)
            sage: A = random_matrix(K, 10, 10, algorithm='unimodular')
            sage: A.rank()
            10
            sage: A = matrix(K, 10, 0)
            sage: A.rank()
            0
        """
        x = self.fetch('rank')
        if x is not None:
            return x
        if self._nrows == 0 or self._ncols == 0:
            return 0
        cdef mzed_t *A = mzed_copy(NULL, self._entries)

        cdef size_t r = mzed_echelonize(A, 0)
        mzed_free(A)
        self.cache('rank', r)
        return r

    def __reduce__(self):
        """
        EXAMPLES::

            sage: K.<a> = GF(2^8)
            sage: A = random_matrix(K,70,70)
            sage: f, s= A.__reduce__()
            sage: from sage.matrix.matrix_gf2e_dense import unpickle_matrix_gf2e_dense_v0
            sage: f == unpickle_matrix_gf2e_dense_v0
            True
            sage: f(*s) == A
            True

        See :issue:`21669`::

            sage: all(f(*s) == B
            ....:     for r,c in [(0,0),(0,1),(1,0)]
            ....:     for B in [Matrix(GF(4, 'a'), r,c)]
            ....:     for f,s in [B.__reduce__()])
            True
        """
        from sage.matrix.matrix_space import MatrixSpace

        cdef Matrix_mod2_dense A
        cdef int r,c

        r, c = self.nrows(), self.ncols()
        if r == 0 or c == 0:
            return unpickle_matrix_gf2e_dense_v0, (None, self.base_ring(), r, c)
        MS = MatrixSpace(GF(2), self._entries.x.nrows, self._entries.x.ncols)
        A = Matrix_mod2_dense.__new__(Matrix_mod2_dense, MS, 0, 0, 0, alloc = False)
        A._entries = mzd_copy( NULL, self._entries.x)
        return unpickle_matrix_gf2e_dense_v0, (A, self.base_ring(), self.nrows(), self.ncols())

    def slice(self):
        r"""
        Unpack this matrix into matrices over `\GF{2}`.

        Elements in `\GF{2^e}` can be represented as `\sum c_i a^i`
        where `a` is a root the minimal polynomial. This function
        returns a tuple of matrices `C` whose entry `C_i[x,y]` is the
        coefficient of `c_i` in `A[x,y]` if this matrix is `A`.

        EXAMPLES::

            sage: K.<a> = GF(2^2)
            sage: A = random_matrix(K, 5, 5)
            sage: A0, A1 = A.slice()
            sage: all(A.list()[i] == A0.list()[i] + a*A1.list()[i] for i in range(25))
            True

            sage: K.<a> = GF(2^3)
            sage: A = random_matrix(K, 5, 5)
            sage: A0, A1, A2 = A.slice()
            sage: all(A.list()[i] == A0.list()[i] + a*A1.list()[i] + a^2*A2.list()[i] for i in range(25))
            True

        Slicing and clinging are inverse operations::

            sage: B = matrix(K, 5, 5)
            sage: B.cling(A0, A1, A2)
            sage: B == A
            True
        """
        if self._entries.finite_field.degree > 4:
            raise NotImplementedError("Slicing is only implemented for degree <= 4.")

        from sage.matrix.matrix_space import MatrixSpace

        MS = MatrixSpace(GF(2), self._nrows, self._ncols)
        cdef mzd_slice_t *a = mzed_slice(NULL, self._entries)

        cdef Matrix_mod2_dense A0, A1, A2, A3
        A0 = Matrix_mod2_dense.__new__(Matrix_mod2_dense, MS, 0, 0, 0, alloc = True)
        A1 = Matrix_mod2_dense.__new__(Matrix_mod2_dense, MS, 0, 0, 0, alloc = True)
        mzd_copy(A0._entries, a.x[0])
        mzd_copy(A1._entries, a.x[1])
        if self._entries.finite_field.degree > 2:
            A2 = Matrix_mod2_dense.__new__(Matrix_mod2_dense, MS, 0, 0, 0, alloc = True)
            mzd_copy(A2._entries, a.x[2])
        if self._entries.finite_field.degree > 3:
            A3 = Matrix_mod2_dense.__new__(Matrix_mod2_dense, MS, 0, 0, 0, alloc = True)
            mzd_copy(A3._entries, a.x[3])

        mzd_slice_free(a)
        if self._entries.finite_field.degree == 2:
            return A0,A1
        elif self._entries.finite_field.degree == 3:
            return A0,A1,A2
        elif self._entries.finite_field.degree == 4:
            return A0,A1,A2,A3

    def cling(self, *C):
        r"""
        Pack the matrices over `\GF{2}` into this matrix over `\GF{2^e}`.

        Elements in `\GF{2^e}` can be represented as `\sum c_i a^i` where
        `a` is a root the minimal polynomial. If this matrix is `A`
        then this function writes `c_i a^i` to the entry `A[x,y]`
        where `c_i` is the entry `C_i[x,y]`.

        INPUT:

        - ``C`` -- list of matrices over GF(2)

        EXAMPLES::

            sage: K.<a> = GF(2^2)
            sage: A = matrix(K, 5, 5)
            sage: A0 = random_matrix(GF(2), 5, 5)
            sage: A1 = random_matrix(GF(2), 5, 5)
            sage: A.cling(A0, A1)
            sage: all(A.list()[i] == A0.list()[i] + a*A1.list()[i] for i in range(25))
            True

        Slicing and clinging are inverse operations::

            sage: B0, B1 = A.slice()
            sage: B0 == A0 and B1 == A1
            True

        TESTS::

            sage: K.<a> = GF(2^2)
            sage: A = matrix(K, 5, 5)
            sage: A0 = random_matrix(GF(2), 5, 5)
            sage: A1 = random_matrix(GF(2), 5, 5)
            sage: A.cling(A0, A1)
            sage: B = copy(A)
            sage: A.cling(A0, A1)
            sage: A == B
            True

            sage: A.cling(A0)
            Traceback (most recent call last):
            ...
            ValueError: The number of input matrices must be equal to the degree of the base field.

            sage: K.<a> = GF(2^5)
            sage: A = matrix(K, 5, 5)
            sage: A0 = random_matrix(GF(2), 5, 5)
            sage: A1 = random_matrix(GF(2), 5, 5)
            sage: A2 = random_matrix(GF(2), 5, 5)
            sage: A3 = random_matrix(GF(2), 5, 5)
            sage: A4 = random_matrix(GF(2), 5, 5)
            sage: A.cling(A0, A1, A2, A3, A4)
            Traceback (most recent call last):
            ...
            NotImplementedError: Cling is only implemented for degree <= 4.
        """
        cdef Py_ssize_t i

        if self._entries.finite_field.degree > 4:
            raise NotImplementedError("Cling is only implemented for degree <= 4.")

        if self._is_immutable:
            raise TypeError("Immutable matrices cannot be modified.")

        if len(C) != self._entries.finite_field.degree:
            raise ValueError("The number of input matrices must be equal to the degree of the base field.")

        cdef mzd_slice_t *v = mzd_slice_init(self._entries.finite_field, self._nrows, self._ncols)
        for i in range(self._entries.finite_field.degree):
            if not isinstance(C[i], Matrix_mod2_dense):
                mzd_slice_free(v)
                raise TypeError("All input matrices must be over GF(2).")
            mzd_copy(v.x[i], (<Matrix_mod2_dense>C[i])._entries)
        mzed_set_ui(self._entries, 0)
        mzed_cling(self._entries, v)
        mzd_slice_free(v)


def unpickle_matrix_gf2e_dense_v0(Matrix_mod2_dense a, base_ring, nrows, ncols):
    r"""
    EXAMPLES::

        sage: K.<a> = GF(2^2)
        sage: A = random_matrix(K,10,10)
        sage: f, s= A.__reduce__()
        sage: from sage.matrix.matrix_gf2e_dense import unpickle_matrix_gf2e_dense_v0
        sage: f == unpickle_matrix_gf2e_dense_v0
        True
        sage: f(*s) == A
        True

    We can still unpickle pickles from before :issue:`19240`::

        sage: old_pickle = b'x\x9c\x85RKo\xd3@\x10\xae\xdd$$\xdb&\xe5U\x1e-\x8f\xc2\xc9\x12RD#$\xce\xa0\xb4\x80\x07\xa2\xca\xc2\x07\x0e\xd5\xe2:\x1b\xdb\x8acg\x1c\xa7J\x85*!\xa4\x90\xe6\x07p\xe0\xc4\x01q\xe5\xc4\x19\xf5\xd0?\xc1\x81\xdf\x80\xb8q\x0b\xb3\x8eMS\xa1\x82V;;\xb3\xdf\xce\xf7\xcd\x8e\xe6\xb5j\xf7,GT;V\x1cy\x83\xf4\xe0\x9d\xb0Y\x13\xbc)\x82\x9e`\xfd\xa0\xeb\xd9m_\xf0\xbf1\xbe{\x97\xa1\xa2\x9d\xc6\xf0\x0f\x82,\x7f\x9d\xa1\xaa\x81\n\xb9m\x9c\xd7\xf4\xf1d2\x81-h\xc0#(\x03\x83\x15\xdas\xc9*\xc3\x13x\x0cu0\xd28\x97\x9e*(0\x9f\xfa\x1b\xd0\xd2\x7fH\x82\xb5\xf4\xa2@TO\xe19\x01I\xac\x136\x991\x9f\xa4\xf9&\xcd\x07i\xbe\xcb\xd4ib\t\xba\xa4\xf6\x02zIT\xd1\x8f2(u\x15\xfd\x9d<\xee@\x05V\xd3\x94E*\xb0\x0e\x0fH\xad\xa8\xbf\x97\xa0\r\x03\xfd\xf0\xb8\x1aU\xff\x92\x90\xe8?\xa5\xd6\x814_\xa5\xf9(\xcd\xafc\xe99\xe2\xd9\xa0\x06\xd4\xf5\xcf\xf2\xf2!\xbc\xd4\xdf\x90#\xc0\x8f\r\xccM\x1b\xdd\x8b\xa3\xbe\x1d\xf7#QmYv\x1cF{\xcc\x11\x81\x88<\x9b\xa71\xcf:\xce0\xaf\x9d\x96\xe3\x87a\xbb\xdf\xe5\x8e\x1f\xeeX>\xc3\x82\xb9\xb0\xe9\x05^,6=\xe17\xf1\xcc\xd0\xc0"u\xb0d\xe6wDl\xdd\x1fa)e\x8a\xbc\xc0\xe9U\xbd \x16\x8e\x88X\xc7j\x0b\x9e\x05\xc8L\xe5\x1e%.\x98\x8a5\xc4\xc5\xd9\xf7\xdd\xd0\xdf\x0b\xc2\x8eg\xf93.wZ\xb5\xc1\x94B\xf8\xa2#\x82\x98a\xf9\xffY\x12\xe3v\x18L\xff\x14Fl\xeb\x0ff\x10\xc4\xb0\xa2\xb9y\xcd-\xba%\xcd\xa5\x8ajT\xd1\x92\xa9\x0c\x86x\xb6a\xe6h\xf8\x02<g\xaa\xaf\xf6\xdd%\x89\xae\x13z\xfe \xc6\x0b\xfb1^4p\x99\x1e6\xc6\xd4\xebK\xdbx\xf9\xc4\x8f[Iw\xf8\x89\xef\xcbQf\xcfh\xe3\x95\x8c\xebj&\xb9\xe2.\x8f\x0c\\ui\x89\xf1x\xf4\xd6\xc0kf\xc1\xf1v\xad(\xc4\xeb\x89~\xfa\xf0\x06\xa8\xa4\x7f\x93\xf4\xd7\x0c\xbcE#\xad\x92\xfc\xed\xeao\xefX\\\x03'
        sage: loads(old_pickle)
        [    0     a]
        [a + 1     1]
    """
    from sage.matrix.matrix_space import MatrixSpace

    MS = MatrixSpace(base_ring, nrows, ncols)
    cdef Matrix_gf2e_dense A  = Matrix_gf2e_dense.__new__(Matrix_gf2e_dense, MS, 0, 0, 0)
    if nrows != 0 and ncols != 0:
        mzd_copy(A._entries.x, a._entries)
    return A
