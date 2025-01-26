# distutils: extra_compile_args = -D_XPG6 NTL_CFLAGS M4RI_CFLAGS
# distutils: extra_link_args = NTL_LIBEXTRA
# distutils: libraries = iml NTL_LIBRARIES m CBLAS_LIBRARIES
# distutils: library_dirs = NTL_LIBDIR CBLAS_LIBDIR
# distutils: include_dirs = NTL_INCDIR M4RI_INCDIR CBLAS_INCDIR
# distutils: language = c++
"""
Dense matrices over the rational field

EXAMPLES:

We create a 3x3 matrix with rational entries and do some
operations with it.

::

    sage: a = matrix(QQ, 3,3, [1,2/3, -4/5, 1,1,1, 8,2, -3/19]); a
    [    1   2/3  -4/5]
    [    1     1     1]
    [    8     2 -3/19]
    sage: a.det()
    2303/285
    sage: a.charpoly()
    x^3 - 35/19*x^2 + 1259/285*x - 2303/285
    sage: b = a^(-1); b
    [ -615/2303  -426/2303   418/2303]
    [ 2325/2303  1779/2303  -513/2303]
    [-1710/2303   950/2303    95/2303]
    sage: b.det()
    285/2303
    sage: a == b
    False
    sage: a < b
    False
    sage: b < a
    True
    sage: a > b
    True
    sage: a*b
    [1 0 0]
    [0 1 0]
    [0 0 1]

TESTS::

    sage: a = matrix(QQ, 2, range(4), sparse=False)
    sage: TestSuite(a).run()

Test hashing::

    sage: m = matrix(QQ, 2, [1/2, -1, 2, 3])
    sage: hash(m)
    Traceback (most recent call last):
    ...
    TypeError: mutable matrices are unhashable
    sage: m.set_immutable()
    sage: hash(m)
    2212268000387745777  # 64-bit
    1997752305           # 32-bit
"""

# ****************************************************************************
#       Copyright (C) 2004,2005,2006 William Stein <wstein@gmail.com>
#                     2017 Vincent Delecroix <20100.delecroix@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from libc.string cimport strcpy, strlen

from sage.categories.rings import Rings
from sage.cpython.string cimport char_to_str, str_to_bytes

from sage.modules.vector_rational_dense cimport Vector_rational_dense
from sage.ext.stdsage cimport PY_NEW
from sage.misc.randstate cimport randstate, current_randstate

from sage.modules.vector_rational_dense cimport Vector_rational_dense

from cysignals.signals cimport sig_on, sig_off
from cysignals.memory cimport sig_malloc, sig_free

from sage.arith.rational_reconstruction cimport mpq_rational_reconstruction

from sage.libs.gmp.types cimport mpz_t, mpq_t
from sage.libs.gmp.mpz cimport mpz_init, mpz_clear, mpz_cmp_si
from sage.libs.gmp.mpq cimport mpq_init, mpq_clear, mpq_set_si, mpq_mul, mpq_add, mpq_set
from sage.libs.gmp.randomize cimport (mpq_randomize_entry,
                                      mpq_randomize_entry_as_int,
                                      mpq_randomize_entry_recip_uniform,
                                      mpq_randomize_entry_nonzero,
                                      mpq_randomize_entry_as_int_nonzero,
                                      mpq_randomize_entry_recip_uniform_nonzero)

from sage.libs.flint.fmpz cimport *
from sage.libs.flint.fmpq cimport *
from sage.libs.flint.fmpz_mat cimport *
from sage.libs.flint.fmpq_mat cimport *

cimport sage.structure.element

from sage.structure.richcmp cimport rich_to_bool
from sage.rings.rational cimport Rational
from sage.matrix.matrix cimport Matrix
from sage.matrix.args cimport SparseEntry, MatrixArgs_init
from sage.matrix.matrix_integer_dense cimport Matrix_integer_dense, _lift_crt
from sage.structure.element cimport Element, Vector
from sage.rings.integer cimport Integer
from sage.rings.integer_ring import ZZ, IntegerRing_class
import sage.rings.abc
from sage.rings.rational_field import QQ

from sage.matrix.matrix2 import decomp_seq
from sage.misc.verbose import verbose

# ########################################################
# PARI C library
from sage.libs.pari.all import PariError
from sage.libs.pari.convert_gmp cimport INTFRAC_to_mpq
from sage.libs.pari.convert_flint cimport rational_matrix, _new_GEN_from_fmpq_mat_t
from cypari2.stack cimport clear_stack
from cypari2.paridecl cimport *
# ########################################################

cdef class Matrix_rational_dense(Matrix_dense):
    def __cinit__(self):
        """
        Create and allocate memory for the matrix.

        EXAMPLES::

            sage: from sage.matrix.matrix_rational_dense import Matrix_rational_dense
            sage: a = Matrix_rational_dense.__new__(Matrix_rational_dense, Mat(ZZ,3), 0,0,0)
            sage: type(a)
            <class 'sage.matrix.matrix_rational_dense.Matrix_rational_dense'>

        .. WARNING::

           This is for internal use only, or if you really know what
           you're doing.
        """
        sig_on()
        fmpq_mat_init(self._matrix, self._nrows, self._ncols)
        sig_off()

    cdef inline Matrix_rational_dense _new_matrix(self, Py_ssize_t nrows, Py_ssize_t ncols):
        if nrows == self._nrows and ncols == self._ncols:
            parent = self._parent
        else:
            parent = self.matrix_space(nrows, ncols)

        return Matrix_rational_dense.__new__(Matrix_rational_dense, parent, None, None, None)

    def __dealloc__(self):
        fmpq_mat_clear(self._matrix)

    def __init__(self, parent, entries=None, copy=None, bint coerce=True):
        r"""
        INPUT:

        - ``parent`` -- a matrix space over ``QQ``

        - ``entries`` -- see :func:`matrix`

        - ``copy`` -- ignored (for backwards compatibility)

        - ``coerce`` -- if ``False``, assume without checking that the
          entries are of type :class:`Rational`

        TESTS::

            sage: matrix(QQ, 2, 2, 1/4)
            [1/4   0]
            [  0 1/4]
            sage: matrix(QQ, 3, 1, [1/2, -3/4, 0])
            [ 1/2]
            [-3/4]
            [   0]
            sage: matrix(QQ, 2, 2, 0.5)
            [1/2   0]
            [  0 1/2]
        """
        ma = MatrixArgs_init(parent, entries)
        cdef Rational z
        for t in ma.iter(coerce, True):
            se = <SparseEntry>t
            z = <Rational>se.entry
            fmpq_set_mpq(fmpq_mat_entry(self._matrix, se.i, se.j), z.value)

    def matrix_from_columns(self, columns):
        """
        Return the matrix constructed from ``self`` using columns with indices
        in the columns list.

        EXAMPLES::

            sage: A = matrix(QQ, 3, range(9))
            sage: A
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: A.matrix_from_columns([2,1])
            [2 1]
            [5 4]
            [8 7]
            sage: A.matrix_from_columns((2,1,0,2))
            [2 1 0 2]
            [5 4 3 5]
            [8 7 6 8]
        """
        cdef Matrix_rational_dense A
        cdef Py_ssize_t k, r, col

        A = self._new_matrix(self._nrows, len(columns))
        k = 0
        for col in columns:
            if col < 0 or col >= self._ncols:
                raise IndexError("column out of range")
            for r in range(self._nrows):
                fmpq_set(fmpq_mat_entry(A._matrix, r, k), fmpq_mat_entry(self._matrix, r, col))
            k = k + 1
        return A

    def add_to_entry(self, Py_ssize_t i, Py_ssize_t j, elt):
        r"""
        Add ``elt`` to the entry at position ``(i,j)``.

        EXAMPLES::

            sage: m = matrix(QQ, 2, 2)
            sage: m.add_to_entry(0, 0, -1/3)
            sage: m
            [-1/3    0]
            [   0    0]
        """
        if not isinstance(elt, Rational):
            elt = Rational(elt)
        if i < 0:
            i += self._nrows
        if i < 0 or i >= self._nrows:
            raise IndexError("row index out of range")
        if j < 0:
            j += self._ncols
        if j < 0 or j >= self._ncols:
            raise IndexError("column index out of range")
        cdef fmpq_t tmp
        fmpq_init(tmp)
        fmpq_set_mpq(tmp, (<Rational>elt).value)
        fmpq_add(fmpq_mat_entry(self._matrix, i, j),
                 fmpq_mat_entry(self._matrix, i, j),
                 tmp)
        fmpq_clear(tmp)

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, value):
        fmpq_set_mpq(fmpq_mat_entry(self._matrix, i, j), (<Rational> value).value)

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        cdef Rational x
        x = Rational.__new__(Rational)
        fmpq_get_mpq(x.value, fmpq_mat_entry(self._matrix, i, j))
        return x

    cdef bint get_is_zero_unsafe(self, Py_ssize_t i, Py_ssize_t j) except -1:
        """
        Return 1 if the entry (i, j) is zero, otherwise 0.

        .. WARNING::

           This is very unsafe; it assumes i and j are in the right
           range.
        """
        return fmpq_is_zero(fmpq_mat_entry(self._matrix, i,j))

    cdef _add_ui_unsafe_assuming_int(self, Py_ssize_t i, Py_ssize_t j, unsigned long int n):
        # doesn't check immutability
        # doesn't do bounds checks.
        # assumes that self[i,j] is an integer.
        cdef fmpz * entry = fmpq_numref(fmpq_mat_entry(self._matrix, i, j))
        fmpz_add_ui(entry, entry, n)

    cdef _sub_ui_unsafe_assuming_int(self, Py_ssize_t i, Py_ssize_t j, unsigned long int n):
        # doesn't check immutability
        # doesn't do bounds checks.
        # assumes that self[i,j] is an integer.
        cdef fmpz * entry = fmpq_numref(fmpq_mat_entry(self._matrix, i, j))
        fmpz_sub_ui(entry, entry, n)

    def _pickle(self):
        return self._pickle_version0(), 0

    def _unpickle(self, data, int version):
        if version == 0:
            self._unpickle_version0(data)
        else:
            raise RuntimeError("unknown matrix version (=%s)" % version)

    cdef _pickle_version0(self):
        return self._export_as_string(32)

    cpdef _export_as_string(self, int base=10):
        """
        Return space separated string of the entries in this matrix, in the
        given base. This is optimized for speed.

        INPUT:

        - ``base`` -- integer (default: `10`)

        EXAMPLES::

            sage: m = matrix(QQ,2,3,[1,2/3,-3/4,1,-2/3,-45/17])
            sage: m._export_as_string(10)
            '1 2/3 -3/4 1 -2/3 -45/17'
            sage: m._export_as_string(16)
            '1 2/3 -3/4 1 -2/3 -2d/11'
        """
        cdef Py_ssize_t i, j, len_so_far, m, n
        cdef char *s
        cdef char *t
        cdef char *tmp

        if self._nrows == 0 or self._ncols == 0:
            data = ''
        else:
            n = self._nrows * self._ncols * 10
            s = <char*> sig_malloc(n * sizeof(char))
            t = s
            len_so_far = 0

            sig_on()
            for i in range(self._nrows):
                for j in range(self._ncols):
                    m = fmpz_sizeinbase (fmpq_mat_entry_num(self._matrix, i, j), base) + \
                        fmpz_sizeinbase (fmpq_mat_entry_den(self._matrix, i, j), base) + 3
                    if len_so_far + m + 1 >= n:
                        # copy to new string with double the size
                        n = 2*n + m + 1
                        tmp = <char*> sig_malloc(n * sizeof(char))
                        strcpy(tmp, s)
                        sig_free(s)
                        s = tmp
                        t = s + len_so_far
                    fmpq_get_str(t, base, fmpq_mat_entry(self._matrix, i, j))
                    m = strlen(t)
                    len_so_far = len_so_far + m + 1
                    t = t + m
                    t[0] = <char>32
                    t[1] = <char>0
                    t = t + 1
            sig_off()
            data = char_to_str(s)[:-1]
            sig_free(s)
        return data

    cdef _unpickle_version0(self, data):
        r"""
        TESTS::

            sage: a = random_matrix(QQ, 4, 3, num_bound=2**500, den_bound=2**500)
            sage: loads(dumps(a)) == a  # indirect doctest
            True
        """
        cdef Py_ssize_t i, j, k
        data = data.split()
        if len(data) != self._nrows * self._ncols:
            raise RuntimeError("invalid pickle data")
        k = 0
        for i in range(self._nrows):
            for j in range(self._ncols):
                s = data[k]
                k += 1
                if '/' in s:
                    num, den = [str_to_bytes(n) for n in s.split('/')]
                    if fmpz_set_str(fmpq_mat_entry_num(self._matrix, i, j), num, 32) or \
                       fmpz_set_str(fmpq_mat_entry_den(self._matrix, i, j), den, 32):
                        raise RuntimeError("invalid pickle data")
                else:
                    s = str_to_bytes(s)
                    if fmpz_set_str(fmpq_mat_entry_num(self._matrix, i, j), s, 32):
                        raise RuntimeError("invalid pickle data")
                    fmpz_one(fmpq_mat_entry_den(self._matrix, i, j))

    # #######################################################################
    # LEVEL 2 functionality
    # x * cdef _add_
    # x * cdef _mul_
    # x * cdef _vector_times_matrix_
    # x * cpdef _richcmp_
    # x * __neg__
    #   * __invert__
    # x * __copy__
    # x * _multiply_classical
    #   * _list -- list of underlying elements (need not be a copy)
    #   * _dict -- sparse dictionary of underlying elements (need not be a copy)
    # #######################################################################

    cpdef _lmul_(self, Element right):
        """
        EXAMPLES::

            sage: a = matrix(QQ, 2, range(6))
            sage: (3/4) * a
            [   0  3/4  3/2]
            [ 9/4    3 15/4]
        """
        cdef Matrix_rational_dense M
        cdef fmpq_t x
        fmpq_init(x)
        fmpq_set_mpq(x, (<Rational>right).value)
        M = Matrix_rational_dense.__new__(Matrix_rational_dense, self._parent, None, None, None)
        fmpq_mat_scalar_mul_fmpz(M._matrix, self._matrix, fmpq_numref(x))
        fmpq_mat_scalar_div_fmpz(M._matrix, M._matrix, fmpq_denref(x))
        fmpq_clear(x)
        return M

    cpdef _add_(self, right):
        """
        Add two dense matrices over QQ.

        EXAMPLES::

            sage: a = MatrixSpace(QQ,3)(range(9))
            sage: b = MatrixSpace(QQ,3)([1/n for n in range(1,10)])
            sage: a+b
            [   1  3/2  7/3]
            [13/4 21/5 31/6]
            [43/7 57/8 73/9]
            sage: b.swap_rows(1,2)
            sage: #a+b
        """
        cdef Matrix_rational_dense ans
        ans = Matrix_rational_dense.__new__(Matrix_rational_dense, self._parent, None, None, None)

        sig_on()
        fmpq_mat_add(ans._matrix, self._matrix, (<Matrix_rational_dense> right)._matrix)
        sig_off()
        return ans

    cpdef _sub_(self, right):
        """
        Subtract two dense matrices over QQ.

        EXAMPLES::

            sage: a = MatrixSpace(QQ,3)(range(9))
            sage: b = MatrixSpace(QQ,3)([1/n for n in range(1,10)])
            sage: a-b
            [  -1  1/2  5/3]
            [11/4 19/5 29/6]
            [41/7 55/8 71/9]
        """
        cdef Matrix_rational_dense ans
        ans = Matrix_rational_dense.__new__(Matrix_rational_dense, self._parent, None, None, None)

        sig_on()
        fmpq_mat_sub(ans._matrix, self._matrix, (<Matrix_rational_dense> right)._matrix)
        sig_off()
        return ans

    cpdef _richcmp_(self, right, int op):
        r"""
        TESTS::

            sage: M = MatrixSpace(QQ, 1)
            sage: M(1) < M(2)
            True
            sage: M(1/3) >= M(5/2)
            False
            sage: M(2) == M(2)
            True
            sage: M(3/4) != M(2)
            True

            sage: matrix(QQ, 2, 3) == matrix(QQ, 2, 3)
            True
            sage: matrix(QQ, 2, 2) == matrix(QQ, 2, 3)
            False
            sage: matrix(QQ, 2, 2) == matrix(QQ, 3, 2)
            False
            sage: matrix(QQ, 2, 3) == matrix(QQ, 3, 2)
            False

            sage: mats = [matrix(QQ, 2, 2, 1), matrix(QQ, 2, 2, -1), matrix(QQ, 2, 2, 0)]
            sage: mats.sort()
            sage: mats == [-1, 0, 1]
            True
        """
        cdef Py_ssize_t i, j
        cdef int k
        for i in range(self._nrows):
            for j in range(self._ncols):
                k = fmpq_cmp(fmpq_mat_entry(self._matrix, i, j),
                             fmpq_mat_entry((<Matrix_rational_dense> right)._matrix, i, j))
                if k:
                    if k > 0:
                        return rich_to_bool(op, 1)
                    else:
                        return rich_to_bool(op, -1)
        return rich_to_bool(op, 0)

    cdef _vector_times_matrix_(self, Vector v):
        r"""
        Return the vector times matrix product.

        INPUT:

        - ``v`` -- a free module element

        OUTPUT: the vector times matrix product v\*A

        EXAMPLES::

            sage: B = matrix(QQ,2, [1,2,3,4])
            sage: V = QQ^2
            sage: w = V([-1,5/2])
            sage: w * B
            (13/2, 8)
        """
        cdef Vector_rational_dense w, ans
        cdef Py_ssize_t i, j
        cdef mpq_t x, y, z

        M = self.row_ambient_module()
        w = <Vector_rational_dense> v
        ans = M.zero_vector()

        mpq_init(x)
        mpq_init(y)
        mpq_init(z)
        for i in range(self._ncols):
            mpq_set_si(x, 0, 1)
            for j in range(self._nrows):
                fmpq_get_mpq(z, fmpq_mat_entry(self._matrix, j, i))
                mpq_mul(y, w._entries[j], z)
                mpq_add(x, x, y)
            mpq_set(ans._entries[i], x)
        mpq_clear(x)
        mpq_clear(y)
        mpq_clear(z)
        return ans

    def __neg__(self):
        """
        Negate a matrix over QQ.

        EXAMPLES::

            sage: a = matrix(QQ, 3, [1/n for n in range(1,10)])
            sage: -a
            [  -1 -1/2 -1/3]
            [-1/4 -1/5 -1/6]
            [-1/7 -1/8 -1/9]
        """
        cdef Matrix_rational_dense ans
        ans = Matrix_rational_dense.__new__(Matrix_rational_dense, self._parent, None, None, None)
        fmpq_mat_neg(ans._matrix, self._matrix)
        return ans

    def __copy__(self):
        """
        Copy a matrix over QQ.

        TESTS::

            sage: a = matrix(QQ, 3, [1/n for n in range(1,10)])
            sage: b = a.__copy__()
            sage: a == b
            True
            sage: a is b
            False
            sage: b[0,0] = 5
            sage: a == b
            False

            sage: a.subdivide(2, 1)
            sage: b = a.__copy__()
            sage: b.subdivisions()
            ([2], [1])
            sage: a.subdivide(2, 2)
            sage: b.subdivisions()
            ([2], [1])
        """
        cdef Matrix_rational_dense ans
        ans = Matrix_rational_dense.__new__(Matrix_rational_dense, self._parent, None, None, None)
        fmpq_mat_set(ans._matrix, self._matrix)
        ans._subdivisions = self._subdivisions
        return ans

    # #######################################################################
    # LEVEL 3 functionality (Optional)
    # x * cdef _sub_
    #   * __deepcopy__
    #   * __invert__
    #   * Matrix windows -- only if you need strassen for that base
    #   * Other functions (list them here):
    # x * denom(self):
    # x * mpz_denom(self, mpz_t d):
    # x * _clear_denom(self):
    # o * echelon_modular(self, height_guess=None):
    # #######################################################################
    def __invert__(self):
        """
        EXAMPLES::

            sage: a = matrix(QQ,3,range(9))
            sage: a.inverse()
            Traceback (most recent call last):
            ...
            ZeroDivisionError: input matrix must be nonsingular
            sage: a = matrix(QQ, 2, [1, 5, 17, 3])
            sage: a.inverse()
            [-3/82  5/82]
            [17/82 -1/82]
            sage: ~matrix(QQ, 2, 3)
            Traceback (most recent call last):
            ...
            ArithmeticError: self must be a square matrix
        """
        return self.inverse()

    def _invert_flint(self):
        r"""
        TESTS::

            sage: matrix(QQ, 2, [1,2,3,4])._invert_flint()
            [  -2    1]
            [ 3/2 -1/2]
            sage: matrix(QQ, 1)._invert_flint()
            Traceback (most recent call last):
            ...
            ZeroDivisionError: input matrix must be nonsingular
        """
        cdef int ret
        cdef Matrix_rational_dense ans
        ans = Matrix_rational_dense.__new__(Matrix_rational_dense, self._parent, None, None, None)
        sig_on()
        ret = fmpq_mat_inv(ans._matrix, self._matrix)
        sig_off()
        if ret == 0:
            raise ZeroDivisionError("input matrix must be nonsingular")
        return ans

    def inverse(self, algorithm=None, check_invertible=True):
        """
        Return the inverse of this matrix.

        INPUT:

        - ``algorithm`` -- an optional specification of an algorithm. It can be one of

          - ``None``: (default) uses flint

          - ``'flint'``: uses flint library

          - ``'pari'``: uses PARI library

          - ``'iml'``: uses IML library

        - ``check_invertible`` -- only used when ``algorithm=iml``; whether to
          check that matrix is invertible

        EXAMPLES::

            sage: a = matrix(QQ,3,[1,2,5,3,2,1,1,1,1,])
            sage: a.inverse()
            [1/2 3/2  -4]
            [ -1  -2   7]
            [1/2 1/2  -2]

            sage: a = matrix(QQ, 2, [1, 5, 17, 3])
            sage: a.inverse(algorithm='flint')
            [-3/82  5/82]
            [17/82 -1/82]
            sage: a.inverse(algorithm='flint')  * a
            [1 0]
            [0 1]

            sage: a = matrix(QQ, 2, [-1, 5, 12, -3])
            sage: a.inverse(algorithm='iml')
            [1/19 5/57]
            [4/19 1/57]
            sage: a.inverse(algorithm='iml') * a
            [1 0]
            [0 1]

            sage: a = matrix(QQ, 4, primes_first_n(16))
            sage: a.inverse(algorithm='pari')
            [   3/11  -12/55    -1/5    2/11]
            [  -5/11   -2/55    3/10   -3/22]
            [ -13/22 307/440   -1/10   -9/88]
            [  15/22  -37/88       0    7/88]

        On singular matrices this method raises a :exc:`ZeroDivisionError`::

            sage: a = matrix(QQ, 2)
            sage: a.inverse(algorithm='flint')
            Traceback (most recent call last):
            ...
            ZeroDivisionError: input matrix must be nonsingular
            sage: a.inverse(algorithm='iml')
            Traceback (most recent call last):
            ...
            ZeroDivisionError: input matrix must be nonsingular
            sage: a.inverse(algorithm='pari')
            Traceback (most recent call last):
            ...
            ZeroDivisionError: input matrix must be nonsingular

        TESTS::

            sage: a = matrix(QQ, 2)
            sage: a.inverse(algorithm="IAmNotAnAlgorithm")
            Traceback (most recent call last):
            ...
            ValueError: unknown algorithm 'IAmNotAnAlgorithm'

            sage: for _ in range(30):
            ....:     dim = randint(1, 20)
            ....:     a = random_matrix(QQ, dim, num_bound=10, den_bound=10)
            ....:     while a.rank() != dim: a = random_matrix(QQ, dim)
            ....:     inv_flint = a.inverse(algorithm='flint')
            ....:     inv_pari = a.inverse(algorithm='pari')
            ....:     inv_iml = a.inverse(algorithm='iml')
            ....:     assert inv_flint == inv_pari == inv_iml
        """
        if self._nrows != self._ncols:
            raise ArithmeticError("self must be a square matrix")

        if self._nrows == 0:
            return self

        if algorithm is None:
            algorithm = "flint"

        if algorithm == "flint":
            return self._invert_flint()
        elif algorithm == "pari":
            try:
                return self._invert_pari()
            except PariError:
                raise ZeroDivisionError("input matrix must be nonsingular")
        elif algorithm == "iml":
            AZ, denom = self._clear_denom()
            B, d = AZ._invert_iml(check_invertible=check_invertible)
            return (denom/d)*B

        else:
            raise ValueError("unknown algorithm '%s'" % algorithm)

    def determinant(self, algorithm=None, proof=None):
        """
        Return the determinant of this matrix.

        INPUT:

        - ``algorithm`` -- an optional specification of an algorithm. It can be one of

          - ``None``: (default) uses flint

          - ``'flint'``: uses flint library

          - ``'pari'``: uses PARI library

          - ``'integer'``: removes denominator and call determinant on the corresponding
             integer matrix

          - ``'generic'``: calls the generic Sage implementation

        - ``proof`` -- boolean or ``None``; if ``None`` use
          proof.linear_algebra(); only relevant for the padic algorithm

        .. NOTE::

           It would be *VERY VERY* hard for det to fail even with
           proof=False.

        EXAMPLES::

            sage: m = matrix(QQ,3,[1,2/3,4/5, 2,2,2, 5,3,2/5])
            sage: m.determinant()
            -34/15
            sage: m.charpoly()
            x^3 - 17/5*x^2 - 122/15*x + 34/15

            sage: m = matrix(QQ, 3, [(1/i)**j for i in range(2,5) for j in range(3)])
            sage: m.determinant(algorithm='flint')
            -1/288

            sage: m = matrix(QQ, 4, [(-1)**n/n for n in range(1,17)])
            sage: m.determinant(algorithm='pari')
            2/70945875

            sage: m = matrix(QQ, 5, [1/(i+j+1) for i in range(5) for j in range(5)])
            sage: m.determinant(algorithm='integer')
            1/266716800000

        On non-square matrices, the method raises a :exc:`ValueError`::

            sage: matrix(QQ, 2, 3).determinant(algorithm='flint')
            Traceback (most recent call last):
            ...
            ValueError: non square matrix
            sage: matrix(QQ, 2, 3).determinant(algorithm='pari')
            Traceback (most recent call last):
            ...
            ValueError: non square matrix
            sage: matrix(QQ, 2, 3).determinant(algorithm='integer')
            Traceback (most recent call last):
            ...
            ValueError: non square matrix
            sage: matrix(QQ, 2, 3).determinant(algorithm='generic')
            Traceback (most recent call last):
            ...
            ValueError: non square matrix

        TESTS:

        Check that the four algorithms agree::

            sage: for _ in range(20):
            ....:     dim = randint(0, 30)
            ....:     m = random_matrix(QQ, dim, num_bound=10, den_bound=10)
            ....:     det_flint = m.determinant("flint"); m._clear_cache()
            ....:     det_pari = m.determinant("pari"); m._clear_cache()
            ....:     det_int = m.determinant("integer"); m._clear_cache()
            ....:     det_gen = m.determinant("generic")
            ....:     assert det_flint == det_pari == det_int == det_gen
        """
        if self._nrows != self._ncols:
            raise ValueError("non square matrix")

        det = self.fetch('det')
        if det is not None:
            return det

        if algorithm is None or algorithm == "flint":
            det = self._det_flint()
        elif algorithm == "pari":
            det = self._det_pari()
        elif algorithm == "integer":
            A, denom = self._clear_denom()
            det = Rational(A.determinant(proof=proof))
            if not denom.is_one():
                det = det / (denom ** self.nrows())
        elif algorithm == "generic":
            det = Matrix_dense.determinant(self)
        else:
            raise ValueError("unknown algorithm '%s'" % algorithm)

        self.cache('det', det)
        return det

    def _det_flint(self):
        r"""
        Return the determinant (computed using flint).

        EXAMPLES::

            sage: matrix(QQ, 2, [1/3, 2/5, 3/4, 7/8])._det_flint()
            -1/120
            sage: matrix(QQ, 0)._det_flint()
            1
            sage: matrix(QQ, 1, [0])._det_flint()
            0
        """
        cdef Rational d = Rational.__new__(Rational)
        cdef fmpq_t e
        fmpq_init(e)
        sig_on()
        fmpq_mat_det(e, self._matrix)
        fmpq_get_mpq(d.value, e)
        sig_off()
        return d

    def denominator(self):
        """
        Return the denominator of this matrix.

        OUTPUT: a Sage Integer

        EXAMPLES::

            sage: b = matrix(QQ,2,range(6)); b[0,0]=-5007/293; b
            [-5007/293         1         2]
            [        3         4         5]
            sage: b.denominator()
            293

            sage: matrix(QQ, 2, [1/2, 1/3, 1/4, 1/5]).denominator()
            60
        """
        cdef Integer z = Integer.__new__(Integer)
        cdef fmpz_t tmp
        fmpz_init(tmp)
        self.fmpz_denom(tmp)
        fmpz_get_mpz(z.value, tmp)
        fmpz_clear(tmp)
        return z

    cdef int fmpz_denom(self, fmpz_t d) except -1:
        cdef Py_ssize_t i, j
        sig_on()
        fmpz_one(d)
        for i in range(self._nrows):
            for j in range(self._ncols):
                fmpz_lcm(d, d, fmpq_mat_entry_den(self._matrix, i, j))
        sig_off()
        return 0

    def _clear_denom(self):
        r"""
        INPUT:

        - ``self`` -- a matrix

        OUTPUT: ``D*self, D``

        The product is a matrix over `\ZZ`.

        EXAMPLES::

            sage: a = matrix(QQ,2,[-1/6,-7,3,5/4]); a
            [-1/6   -7]
            [   3  5/4]
            sage: b, d = a._clear_denom()
            sage: b
            [ -2 -84]
            [ 36  15]
            sage: d
            12
            sage: b == d * a
            True
        """
        X = self.fetch('clear_denom')
        if X is not None:
            return X

        cdef Py_ssize_t i, j
        cdef Matrix_integer_dense A
        cdef fmpz * entry
        cdef fmpz_t denom
        fmpz_init(denom)
        self.fmpz_denom(denom)

        from sage.matrix.matrix_space import MatrixSpace
        MZ = MatrixSpace(ZZ, self._nrows, self._ncols, sparse=False)
        A = Matrix_integer_dense.__new__(Matrix_integer_dense, MZ, None, None, None)

        sig_on()
        for i in range(self._nrows):
            for j in range(self._ncols):
                entry = fmpz_mat_entry(A._matrix, i, j)
                fmpz_divexact(entry, denom, fmpq_mat_entry_den(self._matrix, i, j))
                fmpz_mul(entry, entry, fmpq_mat_entry_num(self._matrix, i, j))
        sig_off()

        cdef Integer D = PY_NEW(Integer)
        fmpz_get_mpz(D.value, denom)
        fmpz_clear(denom)
        X = (A, D)
        self.cache('clear_denom', X)
        return X

    def charpoly(self, var='x', algorithm=None):
        r"""
        Return the characteristic polynomial of this matrix.

        .. NOTE::

            The characteristic polynomial is defined as `\det(xI-A)`.

        INPUT:

        - ``var`` -- (optional) name of the variable as a string

        - ``algorithm`` -- an optional specification of an algorithm. It can be
          one of:

           - ``None``: (default) will use flint for small dimensions and linbox
             otherwise

           - ``'flint'``: uses flint library

           - ``'linbox'``: uses linbox library

           - ``'generic'``: uses Sage generic implementation

        OUTPUT: a polynomial over the rational numbers

        EXAMPLES::

            sage: a = matrix(QQ, 3, [4/3, 2/5, 1/5, 4, -3/2, 0, 0, -2/3, 3/4])
            sage: f = a.charpoly(); f
            x^3 - 7/12*x^2 - 149/40*x + 97/30
            sage: f(a)
            [0 0 0]
            [0 0 0]
            [0 0 0]

        TESTS:

        The cached polynomial should be independent of the ``var``
        argument (:issue:`12292`). We check (indirectly) that the
        second call uses the cached value by noting that its result is
        not cached::

            sage: M = MatrixSpace(QQ, 2)
            sage: A = M(range(0, 2^2))
            sage: type(A)
            <class 'sage.matrix.matrix_rational_dense.Matrix_rational_dense'>
            sage: A.charpoly('x')
            x^2 - 3*x - 2
            sage: A.charpoly('y')
            y^2 - 3*y - 2
            sage: A._cache['charpoly']
            x^2 - 3*x - 2

        Check consistency::

            sage: for _ in range(100):
            ....:     dim = randint(0, 10)
            ....:     m = random_matrix(QQ, dim, num_bound=8, den_bound=8)
            ....:     p_flint = m.charpoly(algorithm='flint'); m._clear_cache()
            ....:     p_linbox = m.charpoly(algorithm='linbox'); m._clear_cache()
            ....:     p_generic = m.charpoly(algorithm='generic')
            ....:     assert p_flint == p_linbox == p_generic
        """
        poly = self.fetch('charpoly')
        if poly is not None:
            return poly.change_variable_name(var)

        if algorithm is None:
            algorithm = 'flint' if self._nrows <= 40 else 'linbox'

        if algorithm == 'flint' or algorithm == 'linbox':
            A, denom = self._clear_denom()
            f = A.charpoly(var, algorithm=algorithm)
            x = f.parent().gen()
            g = f(x * denom) / denom ** f.degree()
        elif algorithm == 'generic':
            g = Matrix_dense.charpoly(self, var)
        else:
            raise ValueError("no algorithm '%s'" % algorithm)

        self.cache('charpoly', g)
        return g

    def minpoly(self, var='x', algorithm=None):
        """
        Return the minimal polynomial of this matrix.

        INPUT:

        - ``var`` -- (optional) the variable name as a string (default: ``'x'``)

        - ``algorithm`` -- an optional specification of an algorithm. It can
          be one of

           - ``None``: (default) will use linbox

           - ``'linbox'``: uses the linbox library

           - ``'generic'``: uses the generic Sage implementation

        OUTPUT: a polynomial over the rationals

        EXAMPLES::

            sage: a = matrix(QQ, 3, [4/3, 2/5, 1/5, 4, -3/2, 0, 0, -2/3, 3/4])
            sage: f = a.minpoly(); f
            x^3 - 7/12*x^2 - 149/40*x + 97/30
            sage: a = Mat(ZZ,4)(range(16))
            sage: f = a.minpoly(); f.factor()
            x * (x^2 - 30*x - 80)
            sage: f(a) == 0
            True

        ::

            sage: a = matrix(QQ, 4, [1..4^2])
            sage: factor(a.minpoly())
            x * (x^2 - 34*x - 80)
            sage: factor(a.minpoly('y'))
            y * (y^2 - 34*y - 80)
            sage: factor(a.charpoly())
            x^2 * (x^2 - 34*x - 80)
            sage: b = matrix(QQ, 4, [-1, 2, 2, 0, 0, 4, 2, 2, 0, 0, -1, -2, 0, -4, 0, 4])
            sage: a = matrix(QQ, 4, [1, 1, 0,0, 0,1,0,0, 0,0,5,0, 0,0,0,5])
            sage: c = b^(-1)*a*b
            sage: factor(c.minpoly())
            (x - 5) * (x - 1)^2
            sage: factor(c.charpoly())
            (x - 5)^2 * (x - 1)^2

        Check consistency::

            sage: for _ in range(100):
            ....:     dim = randint(0, 10)
            ....:     m = random_matrix(QQ, dim, num_bound=8, den_bound=8)
            ....:     p_linbox = m.charpoly(algorithm='linbox'); m._clear_cache()
            ....:     p_generic = m.charpoly(algorithm='generic')
            ....:     assert p_linbox == p_generic
        """
        poly = self.fetch('minpoly')
        if poly is not None:
            return poly.change_variable_name(var)

        if algorithm is None:
            algorithm = 'linbox'

        if algorithm == 'linbox':
            A, denom = self._clear_denom()
            f = A.minpoly(var, algorithm='linbox')
            x = f.parent().gen()
            g = f(x * denom) / denom**f.degree()
        elif algorithm == 'generic':
            g = Matrix_dense.minpoly(self, var)
        else:
            raise ValueError("no algorithm '%s'" % algorithm)

        self.cache('minpoly', g)
        return g

    cdef sage.structure.element.Matrix _matrix_times_matrix_(self, sage.structure.element.Matrix right):
        """
        EXAMPLES::

            sage: a = matrix(QQ, 3, range(9))/3
            sage: b = matrix(QQ, 3, range(1, 10))/5
            sage: a * b   # indirect doctest
            [ 6/5  7/5  8/5]
            [18/5 22/5 26/5]
            [   6 37/5 44/5]

            sage: matrix(QQ, 2, 3) * matrix(QQ, 4, 5)
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for *: 'Full MatrixSpace of 2 by 3 dense matrices over Rational Field' and 'Full MatrixSpace of 4 by 5 dense matrices over Rational Field'
        """
        return self._multiply_flint(right)

    def _multiply_flint(self, Matrix_rational_dense right):
        r"""
        Multiply this matrix by ``right`` using the flint library.

        EXAMPLES::

            sage: n = 3
            sage: a = matrix(QQ,n,range(n^2))/3
            sage: b = matrix(QQ,n,range(1, n^2 + 1))/5
            sage: a._multiply_flint(b)
            [ 6/5  7/5  8/5]
            [18/5 22/5 26/5]
            [   6 37/5 44/5]
        """
        if self._nrows == right._nrows:
            # self acts on the space of right
            parent = right.parent()
        if self._ncols == right._ncols:
            # right acts on the space of self
            parent = self.parent()
        else:
            parent = self.matrix_space(self._nrows, right._ncols)

        cdef Matrix_rational_dense ans
        ans = Matrix_rational_dense.__new__(Matrix_rational_dense, parent, None, None, None)

        sig_on()
        fmpq_mat_mul(ans._matrix, self._matrix, (<Matrix_rational_dense> right)._matrix)
        sig_off()
        return ans

    def _multiply_over_integers(self, Matrix_rational_dense right, algorithm='default'):
        r"""
        Multiply this matrix by right using a multimodular algorithm and
        return the result.

        INPUT:

        - ``self`` -- matrix over QQ

        - ``right`` -- matrix over QQ

        - ``algorithm``

           - 'default': use whatever is the default for A\*B when A, B
             are over ZZ.

           - 'multimodular': use a multimodular algorithm

        EXAMPLES::

            sage: a = MatrixSpace(QQ,10,5)(range(50))
            sage: b = MatrixSpace(QQ,5,12)([1/n for n in range(1,61)])
            sage: a._multiply_over_integers(b) == a._multiply_over_integers(b, algorithm='multimodular')
            True

        ::

            sage: a = MatrixSpace(QQ,3)(range(9))
            sage: b = MatrixSpace(QQ,3)([1/n for n in range(1,10)])
            sage: c = a._multiply_over_integers(b, algorithm = 'multimodular')
            sage: c
            [ 15/28   9/20   7/18]
            [  33/7 117/40   20/9]
            [249/28   27/5  73/18]
            sage: c == a._multiply_flint(b)
            True
        """
        cdef Matrix_integer_dense A, B, AB
        cdef Matrix_rational_dense res
        cdef Integer D
        sig_on()
        A, A_denom = self._clear_denom()
        B, B_denom = right._clear_denom()
        if algorithm == 'default' or algorithm == 'multimodular':
            AB = A*B
        else:
            sig_off()
            raise ValueError("unknown algorithm '%s'" % algorithm)
        D = A_denom * B_denom
        if self._nrows == right._nrows:
            # self acts on the space of right
            res = Matrix_rational_dense.__new__(Matrix_rational_dense, right.parent(), 0, 0, 0)
        if self._ncols == right._ncols:
            # right acts on the space of self
            res = Matrix_rational_dense.__new__(Matrix_rational_dense, self.parent(), 0, 0, 0)
        else:
            res = Matrix_rational_dense.__new__(Matrix_rational_dense, self.matrix_space(AB._nrows, AB._ncols), 0, 0, 0)
        for i in range(res._nrows):
            for j in range(res._ncols):
                fmpz_set(fmpq_mat_entry_num(res._matrix, i, j), fmpz_mat_entry(AB._matrix,i,j))
                fmpz_set_mpz(fmpq_mat_entry_den(res._matrix, i, j), D.value)
                fmpq_canonicalise(fmpq_mat_entry(res._matrix, i, j))
        sig_off()
        return res

    def height(self):
        """
        Return the height of this matrix, which is the maximum of the
        absolute values of all numerators and denominators of entries in
        this matrix.

        OUTPUT: integer

        EXAMPLES::

            sage: b = matrix(QQ,2,range(6)); b[0,0]=-5007/293; b
            [-5007/293         1         2]
            [        3         4         5]
            sage: b.height()
            5007
        """
        cdef Integer z
        cdef fmpz_t tmp
        fmpz_init(tmp)
        self.fmpz_height(tmp)
        z = PY_NEW(Integer)
        fmpz_get_mpz(z.value, tmp)
        fmpz_clear(tmp)
        return z

    cdef int fmpz_height(self, fmpz_t h) except -1:
        cdef fmpz_t x
        cdef int i, j
        sig_on()
        fmpz_init(x)
        fmpz_zero(h)
        for i in range(self._nrows):
            for j in range(self._ncols):
                fmpz_abs(x, fmpq_mat_entry_num(self._matrix, i, j))
                if fmpz_cmp(h, x) < 0:
                    fmpz_set(h, x)
                fmpz_abs(x, fmpq_mat_entry_den(self._matrix, i, j))
                if fmpz_cmp(h, x) < 0:
                    fmpz_set(h, x)
        fmpz_clear(x)
        sig_off()
        return 0

    def _adjugate(self):
        """
        Return the adjugate of this matrix.

        Assumes ``self`` is a square matrix (checked in adjugate).

        EXAMPLES::

            sage: m = matrix(QQ,3,[1..9])/9; m
            [1/9 2/9 1/3]
            [4/9 5/9 2/3]
            [7/9 8/9   1]
            sage: m.adjugate()
            [-1/27  2/27 -1/27]
            [ 2/27 -4/27  2/27]
            [-1/27  2/27 -1/27]
        """
        return self.parent()(self.__pari__().matadjoint().sage())

    def _magma_init_(self, magma):
        """
        EXAMPLES::

            sage: m = matrix(QQ,2,3,[1,2/3,-3/4,1,-2/3,-45/17])
            sage: m._magma_init_(magma)
            'Matrix(RationalField(),2,3,StringToIntegerSequence("204 136 -153 204 -136 -540"))/204'
            sage: magma(m)                                                # optional - magma
            [     1    2/3   -3/4]
            [     1   -2/3 -45/17]
        """
        X, d = self._clear_denom()
        s = X._magma_init_(magma).replace('IntegerRing','RationalField')
        if d != 1:
            s += '/%s' % d._magma_init_(magma)
        return s

    def prod_of_row_sums(self, cols):
        cdef Py_ssize_t i, c
        cdef fmpq_t s, pr
        fmpq_init(s)
        fmpq_init(pr)

        fmpq_one(pr)
        for i in range(self._nrows):
            fmpq_zero(s)
            for c in cols:
                if c < 0 or c >= self._ncols:
                    raise IndexError("matrix column index out of range")
                fmpq_add(s, s, fmpq_mat_entry(self._matrix, i, c))
            fmpq_mul(pr, pr, s)
        cdef Rational ans
        ans = Rational.__new__(Rational)
        fmpq_get_mpq(ans.value, pr)
        fmpq_clear(s)
        fmpq_clear(pr)
        return ans

    def _right_kernel_matrix(self, **kwds):
        r"""
        Return a pair that includes a matrix of basis vectors
        for the right kernel of ``self``.

        INPUT:

        - ``kwds`` -- these are provided for consistency with other versions
          of this method.  Here they are ignored as there is no optional
          behavior available.

        OUTPUT:

        Returns a pair.  First item is the string 'computed-iml-rational'
        that identifies the nature of the basis vectors.

        Second item is a matrix whose rows are a basis for the right kernel,
        over the rationals, as computed by the IML library.  Notice that the
        IML library returns a matrix that is in the 'pivot' format, once the
        whole matrix is multiplied by -1.  So the 'computed' format is very
        close to the 'pivot' format.

        EXAMPLES::

            sage: A = matrix(QQ, [
            ....:                 [1, 0, 1, -3, 1],
            ....:                 [-5, 1, 0, 7, -3],
            ....:                 [0, -1, -4, 6, -2],
            ....:                 [4, -1, 0, -6, 2]])
            sage: result = A._right_kernel_matrix()
            sage: result[0]
            'computed-iml-rational'
            sage: result[1]
            [-1  2 -2 -1  0]
            [ 1  2  0  0 -1]
            sage: X = result[1].transpose()
            sage: A*X == zero_matrix(QQ, 4, 2)
            True

        Computed result is the negative of the pivot basis, which
        is just slightly more efficient to compute. ::

            sage: A.right_kernel_matrix(basis='pivot') == -A.right_kernel_matrix(basis='computed')
            True

        TESTS:

        We test three trivial cases. ::

            sage: A = matrix(QQ, 0, 2)
            sage: A._right_kernel_matrix()[1]
            [1 0]
            [0 1]
            sage: A = matrix(QQ, 2, 0)
            sage: A._right_kernel_matrix()[1].parent()
            Full MatrixSpace of 0 by 0 dense matrices over Rational Field
            sage: A = zero_matrix(QQ, 4, 3)
            sage: A._right_kernel_matrix()[1]
            [1 0 0]
            [0 1 0]
            [0 0 1]
       """
        tm = verbose("computing right kernel matrix over the rationals for %sx%s matrix" % (self.nrows(), self.ncols()),level=1)
        # _rational_kernel_flint() gets the zero-row case wrong, fix it there
        if self.nrows()==0:
            from sage.matrix.constructor import identity_matrix
            K = identity_matrix(QQ, self.ncols())
        else:
            A, _ = self._clear_denom()
            K = A._rational_kernel_iml().transpose().change_ring(QQ)
        verbose("done computing right kernel matrix over the rationals for %sx%s matrix" % (self.nrows(), self.ncols()),level=1, t=tm)
        return 'computed-iml-rational', K

    # ###############################################
    # Change ring
    # ###############################################
    def change_ring(self, R):
        """
        Create the matrix over R with entries the entries of ``self`` coerced
        into R.

        EXAMPLES::

            sage: a = matrix(QQ,2,[1/2,-1,2,3])
            sage: a.change_ring(GF(3))
            [2 2]
            [2 0]
            sage: a.change_ring(ZZ)
            Traceback (most recent call last):
            ...
            TypeError: matrix has denominators so can...t change to ZZ
            sage: b = a.change_ring(QQ['x']); b
            [1/2  -1]
            [  2   3]
            sage: b.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Rational Field

        TESTS:

        Make sure that subdivisions are preserved when changing rings::

            sage: a = matrix(QQ, 3, range(9))
            sage: a.subdivide(2,1); a
            [0|1 2]
            [3|4 5]
            [-+---]
            [6|7 8]
            sage: a.change_ring(ZZ).change_ring(QQ)
            [0|1 2]
            [3|4 5]
            [-+---]
            [6|7 8]
            sage: a.change_ring(GF(3))
            [0|1 2]
            [0|1 2]
            [-+---]
            [0|1 2]
        """
        if R not in Rings():
            raise TypeError("R must be a ring")
        if R == self._base_ring:
            if self._is_immutable:
                return self
            return self.__copy__()
        if isinstance(R, IntegerRing_class):
            A, d = self._clear_denom()
            if not d.is_one():
                raise TypeError("matrix has denominators so can't change to ZZ")
            if self._subdivisions is not None:
                A.subdivide(self.subdivisions())
            return A

        from sage.matrix.matrix_modn_dense_double import MAX_MODULUS
        if isinstance(R, sage.rings.abc.IntegerModRing) and R.order() < MAX_MODULUS:
            b = R.order()
            A, d = self._clear_denom()
            if not b.gcd(d).is_one():
                raise TypeError("matrix denominator not coprime to modulus")
            B = A._mod_int(b)
            C = (1/(B.base_ring()(d))) * B
            if self._subdivisions is not None:
                C.subdivide(self.subdivisions())
            return C

        # fallback to the generic version
        return Matrix_dense.change_ring(self, R)

    # ###############################################
    # Echelon form
    # ###############################################
    def echelonize(self, algorithm=None,
                   height_guess=None, proof=None, **kwds):
        """
        Transform the matrix ``self`` into reduced row echelon form
        in place.

        INPUT:

        - ``algorithm`` -- an optional specification of an algorithm. One of

          - ``None``: (default) uses flint for small dimension and multimodular otherwise

          - ``'flint'``: use the flint library,

          - ``'padic'``: an algorithm based on the IML `p`-adic solver,

          - ``'multimodular'``: uses a multimodular algorithm the uses
            linbox modulo many primes (likely to be faster when coefficients
            are huge),

          - ``'classical'``: just clear each column using Gauss elimination.

        - ``height_guess``, ``**kwds`` -- all passed to the
          multimodular algorithm; ignored by other algorithms

        - ``proof`` -- boolean or ``None`` (default: None, see
          proof.linear_algebra or sage.structure.proof). Passed to the
          multimodular algorithm. Note that the Sage global default is
          ``proof=True``.

        EXAMPLES::

            sage: a = matrix(QQ, 4, range(16)); a[0,0] = 1/19; a[0,1] = 1/5; a
            [1/19  1/5    2    3]
            [   4    5    6    7]
            [   8    9   10   11]
            [  12   13   14   15]
            sage: a.echelonize()
            sage: a
            [      1       0       0 -76/157]
            [      0       1       0  -5/157]
            [      0       0       1 238/157]
            [      0       0       0       0]

        ::

            sage: a = matrix(QQ, 4, range(16)); a[0,0] = 1/19; a[0,1] = 1/5
            sage: a.echelonize(algorithm='multimodular')
            sage: a
            [      1       0       0 -76/157]
            [      0       1       0  -5/157]
            [      0       0       1 238/157]
            [      0       0       0       0]

        TESTS:

        Echelonizing a matrix in place throws away the cache of
        the old matrix (:issue:`14506`)::

            sage: for algo in ["flint", "padic", "multimodular", "classical"]:
            ....:      a = Matrix(QQ, [[1,2],[3,4]])
            ....:      _ = a.det()          # fills the cache
            ....:      _ = a._clear_denom() # fills the cache
            ....:      a.echelonize(algorithm=algo)
            ....:      assert sorted(a._cache.keys()) == ['echelon_form', 'in_echelon_form', 'pivots', 'rank'], (algo, a._cache.keys())
        """
        if self.fetch('in_echelon_form'):
            return  # already known to be in echelon form
        self.check_mutability()

        if algorithm is None:
            if self._nrows <= 25 or self._ncols <= 25:
                algorithm = 'flint'
            else:
                algorithm = 'multimodular'

        if algorithm == 'flint':
            pivots = self._echelonize_flint()
        elif algorithm == 'multimodular':
            pivots = self._echelonize_multimodular(height_guess, proof, **kwds)
        elif algorithm == 'classical':
            pivots = self._echelon_in_place_classical()
        elif algorithm == 'padic':
            pivots = self._echelonize_padic()
        else:
            raise ValueError("no algorithm '%s'" % algorithm)

        if type(pivots) is not tuple:
            raise RuntimeError("BUG: pivots must get set as a tuple. Got {} for algo {} with {}x{} matrix.".format(
                type(pivots), algorithm, self._nrows, self._ncols))

        self.cache('in_echelon_form', True)
        self.cache('echelon_form', self)
        self.cache('pivots', pivots)
        self.cache('rank', len(pivots))

    def echelon_form(self, algorithm=None,
                     height_guess=None, proof=None, **kwds):
        r"""
        Return the echelon form of this matrix.

        The (row) echelon form of a matrix, see :wikipedia:`Row_echelon_form`,
        is the matrix obtained by performing Gauss elimination on the rows
        of the matrix.

        INPUT: See :meth:`echelonize` for the options.

        EXAMPLES::

            sage: a = matrix(QQ, 4, range(16)); a[0,0] = 1/19; a[0,1] = 1/5; a
            [1/19  1/5    2    3]
            [   4    5    6    7]
            [   8    9   10   11]
            [  12   13   14   15]
            sage: a.echelon_form()
            [      1       0       0 -76/157]
            [      0       1       0  -5/157]
            [      0       0       1 238/157]
            [      0       0       0       0]
            sage: a.echelon_form(algorithm='multimodular')
            [      1       0       0 -76/157]
            [      0       1       0  -5/157]
            [      0       0       1 238/157]
            [      0       0       0       0]

        The result is an immutable matrix, so if you want to modify the result
        then you need to make a copy.  This checks that :issue:`10543` is
        fixed.::

            sage: A = matrix(QQ, 2, range(6))
            sage: E = A.echelon_form()
            sage: E.is_mutable()
            False
            sage: F = copy(E)
            sage: F[0,0] = 50
            sage: F
            [50  0 -1]
            [ 0  1  2]

        TESTS:

        Check consistency::

            sage: for _ in range(100):
            ....:     nrows = randint(0, 30)
            ....:     ncols = randint(0, 30)
            ....:     m = random_matrix(QQ, nrows, ncols, num_bound=10, den_bound=10)
            ....:     ech_flint = m.echelon_form('flint'); m._clear_cache()
            ....:     ech_padic = m.echelon_form('padic'); m._clear_cache()
            ....:     ech_multi = m.echelon_form('multimodular'); m._clear_cache()
            ....:     ech_class = m.echelon_form('classical')
            ....:     assert ech_flint == ech_padic == ech_multi == ech_class
        """
        x = self.fetch('echelon_form')
        if x is not None:
            return x
        if self.fetch('in_echelon_form'):
            raise RuntimeError('in_echelon_form set but not echelon_form')

        E = self.__copy__()
        E.echelonize(algorithm)
        E.set_immutable()
        self.cache('echelon_form', E)
        self.cache('pivots', E.pivots())
        self.cache('rank', len(E.pivots()))
        return E

    def _echelonize_flint(self):
        r"""
        EXAMPLES::

            sage: m = matrix(QQ, 4, range(16))
            sage: m._echelonize_flint()
            (0, 1)
            sage: m
            [ 1  0 -1 -2]
            [ 0  1  2  3]
            [ 0  0  0  0]
            [ 0  0  0  0]
            sage: m = matrix(QQ, 4, 6, [-1,0,0,-2,-1,-2,-1,0,0,-2,-1,0,3,3,-2,0,0,3,-2,-3,1,1,-2,3])
            sage: m._echelonize_flint()
            (0, 1, 2, 5)
            sage: m
            [   1    0    0    2    1    0]
            [   0    1    0 -4/3    1    0]
            [   0    0    1    1    3    0]
            [   0    0    0    0    0    1]
        """
        self.clear_cache()
        cdef long r

        sig_on()
        r = fmpq_mat_rref(self._matrix, self._matrix)
        sig_off()

        # compute pivots
        cdef long i, j, k
        cdef list p = []
        k = 0
        for i in range(r):
            for j in range(k, self._ncols):
                if not fmpq_is_zero(fmpq_mat_entry(self._matrix, i, j)):
                    p.append(j)
                    k = j+1  # so start at next position next time
                    break
            else:
                break
        return tuple(p)

    def _echelonize_padic(self):
        """
        Echelonize ``self`` using a `p`-adic nullspace algorithm.

        EXAMPLES::

            sage: m = matrix(QQ, 4, range(16))
            sage: m._echelonize_padic()
            (0, 1)
            sage: m
            [ 1  0 -1 -2]
            [ 0  1  2  3]
            [ 0  0  0  0]
            [ 0  0  0  0]

            sage: m = matrix(QQ, 4, 6, [-1,0,0,-2,-1,-2,-1,0,0,-2,-1,0,3,3,-2,0,0,3,-2,-3,1,1,-2,3])
            sage: m._echelonize_padic()
            (0, 1, 2, 5)
            sage: m
            [   1    0    0    2    1    0]
            [   0    1    0 -4/3    1    0]
            [   0    0    1    1    3    0]
            [   0    0    0    0    0    1]
        """
        cdef Matrix_integer_dense X
        cdef Integer d
        cdef fmpq * entry

        A, _ = self._clear_denom()
        pivots, nonpivots, X, d = A._rational_echelon_via_solve()
        self.clear_cache()

        # FIXME: we should always have X.nrows() == len(pivots)
        if X.nrows() != len(pivots):
            assert X.ncols() == len(pivots) == 0
            assert type(pivots) is list
            fmpq_mat_zero(self._matrix)
            return tuple(pivots)

        cdef Py_ssize_t i,j
        for i in range(X.nrows()):
            # 1 at pivot
            fmpq_one(fmpq_mat_entry(self._matrix, i, pivots[i]))

            # nonzero part
            for j in range(X.ncols()):
                entry = fmpq_mat_entry(self._matrix, i, nonpivots[j])
                fmpz_set(fmpq_numref(entry), fmpz_mat_entry(X._matrix, i, j))
                fmpz_set_mpz(fmpq_denref(entry), d.value)
                fmpq_canonicalise(entry)

            # zeros on the left of the pivot
            for j in range(pivots[i]):
                fmpq_zero(fmpq_mat_entry(self._matrix, i, j))

            # zeros on top of the other pivots
            for j in range(i):
                fmpq_zero(fmpq_mat_entry(self._matrix, j, pivots[i]))

        # Fill in the 0-rows at the bottom.
        for i in range(len(pivots), self._nrows):
            for j in range(self._ncols):
                fmpq_zero(fmpq_mat_entry(self._matrix, i, j))

        # FIXME: pivots should already be a tuple in all cases
        return tuple(pivots)

    def _echelonize_multimodular(self, height_guess=None, proof=None):
        """
        Echelonize ``self`` using multimodular recomposition.

        REFERENCE:

        - Chapter 7 of Stein's "Explicitly Computing Modular Forms".

        INPUT:

        - ``height_guess`` -- integer or ``None``

        - ``proof`` -- boolean (default: ``None``, see
          proof.linear_algebra or sage.structure.proof); Note that the Sage
          global default is ``proof=True``

        EXAMPLES::

            sage: m = matrix(QQ, 4, range(16))
            sage: m._echelonize_multimodular()
            (0, 1)
            sage: m
            [ 1  0 -1 -2]
            [ 0  1  2  3]
            [ 0  0  0  0]
            [ 0  0  0  0]
            sage: m = matrix(QQ, 4, 6, [-1,0,0,-2,-1,-2,-1,0,0,-2,-1,0,3,3,-2,0,0,3,-2,-3,1,1,-2,3])
            sage: m._echelonize_multimodular()
            (0, 1, 2, 5)
            sage: m
            [   1    0    0    2    1    0]
            [   0    1    0 -4/3    1    0]
            [   0    0    1    1    3    0]
            [   0    0    0    0    0    1]
        """
        from sage.matrix.misc import matrix_rational_echelon_form_multimodular
        E, pivots = matrix_rational_echelon_form_multimodular(self, height_guess, proof=proof)
        self.clear_cache()
        fmpq_mat_swap(self._matrix, (<Matrix_rational_dense>E)._matrix)
        return pivots

    cdef swap_rows_c(self, Py_ssize_t r1, Py_ssize_t r2):
        """
        EXAMPLES::

            sage: a = matrix(QQ,2,[1..6])
            sage: a.swap_rows(0,1)             # indirect doctest
            sage: a
            [4 5 6]
            [1 2 3]
        """
        # no bounds checking!
        cdef Py_ssize_t c
        for c in range(self._ncols):
            fmpq_swap(fmpq_mat_entry(self._matrix, r1, c),
                      fmpq_mat_entry(self._matrix, r2, c))

    cdef swap_columns_c(self, Py_ssize_t c1, Py_ssize_t c2):
        """
        EXAMPLES::

            sage: a = matrix(QQ,2,[1..6])
            sage: a.swap_columns(0,1)          # indirect doctest
            sage: a
            [2 1 3]
            [5 4 6]
        """
        # no bounds checking!
        for r in range(self._nrows):
            fmpq_swap(fmpq_mat_entry(self._matrix, r, c1),
                      fmpq_mat_entry(self._matrix, r, c2))

    def decomposition(self, is_diagonalizable=False, dual=False,
                      algorithm=None, height_guess=None, proof=None):
        """
        Return the decomposition of the free module on which this matrix A
        acts from the right (i.e., the action is x goes to x A), along with
        whether this matrix acts irreducibly on each factor. The factors
        are guaranteed to be sorted in the same way as the corresponding
        factors of the characteristic polynomial.

        Let A be the matrix acting from the on the vector space V of column
        vectors. Assume that A is square. This function computes maximal
        subspaces W_1, ..., W_n corresponding to Galois conjugacy classes
        of eigenvalues of A. More precisely, let f(X) be the characteristic
        polynomial of A. This function computes the subspace
        `W_i = ker(g_(A)^n)`, where g_i(X) is an irreducible
        factor of f(X) and g_i(X) exactly divides f(X). If the optional
        parameter is_diagonalizable is True, then we let W_i = ker(g(A)),
        since then we know that ker(g(A)) = `ker(g(A)^n)`.

        If dual is True, also returns the corresponding decomposition of V
        under the action of the transpose of A. The factors are guaranteed
        to correspond.

        INPUT:

        - ``is_diagonalizable`` -- ignored

        - ``dual`` -- whether to also return decompositions for
          the dual

        - ``algorithm`` -- an optional specification of an algorithm

           - ``None`` -- (default) use default algorithm for computing Echelon
             forms

           - 'multimodular': much better if the answers
             factors have small height

        - ``height_guess`` -- positive integer; only used by
          the multimodular algorithm

        - ``proof`` -- boolean or ``None`` (default: ``None``, see
          proof.linear_algebra or sage.structure.proof); only used by the
          multimodular algorithm. Note that the Sage global default is
          proof=True.


        .. NOTE::

           IMPORTANT: If you expect that the subspaces in the answer
           are spanned by vectors with small height coordinates, use
           algorithm='multimodular' and height_guess=1; this is
           potentially much faster than the default. If you know for a
           fact the answer will be very small, use
           algorithm='multimodular', height_guess=bound on height,
           proof=False.

        You can get very very fast decomposition with proof=False.

        EXAMPLES::

            sage: a = matrix(QQ,3,[1..9])
            sage: a.decomposition()
            [(Vector space of degree 3 and dimension 1 over Rational Field
              Basis matrix:
              [ 1 -2  1],
              True),
             (Vector space of degree 3 and dimension 2 over Rational Field
              Basis matrix:
              [ 1  0 -1]
              [ 0  1  2],
              True)]
        """
        X = self._decomposition_rational(is_diagonalizable=is_diagonalizable,
                                         echelon_algorithm=algorithm,
                                         height_guess=height_guess, proof=proof)
        if not dual:
            return X

        Y = self.transpose()._decomposition_rational(is_diagonalizable=is_diagonalizable,
            echelon_algorithm=algorithm, height_guess=height_guess, proof=proof)
        return X, Y

    def _decomposition_rational(self, is_diagonalizable = False,
                                echelon_algorithm=None,
                                kernel_algorithm='default',
                                **kwds):
        """
        Return the decomposition of the free module on which this matrix A
        acts from the right (i.e., the action is x goes to x A), along with
        whether this matrix acts irreducibly on each factor. The factors
        are guaranteed to be sorted in the same way as the corresponding
        factors of the characteristic polynomial.

        INPUT:

        - ``self`` -- a square matrix over the rational
          numbers

        - ``echelon_algorithm`` -- an optional algorithm to be passed to the
          method ``echelon_form``

        - ``'multimodular'`` -- use this if the answers have
          small height

        - ``**kwds`` -- passed on to echelon function

        .. NOTE::

           IMPORTANT: If you expect that the subspaces in the answer are
           spanned by vectors with small height coordinates, use
           algorithm='multimodular' and height_guess=1; this is potentially
           much faster than the default. If you know for a fact the answer
           will be very small, use algorithm='multimodular',
           height_guess=bound on height, proof=False

        OUTPUT:

        - ``Sequence`` -- list of tuples (V,t), where V is a
          vector spaces and t is ``True`` if and only if the charpoly of ``self`` on
          V is irreducible. The tuples are in order corresponding to the
          elements of the sorted list self.charpoly().factor().
        """
        cdef Py_ssize_t k

        if not self.is_square():
            raise ArithmeticError("self must be a square matrix")

        if self.nrows() == 0:
            return decomp_seq([])

        A, _ = self._clear_denom()

        f = A.charpoly('x')
        E = decomp_seq([])

        t = verbose('factoring the characteristic polynomial', level=2, caller_name='rational decomp')
        F = f.factor()
        verbose('done factoring', t=t, level=2, caller_name='rational decomp')

        if len(F) == 1:
            V = QQ**self.nrows()
            m = F[0][1]
            return decomp_seq([(V, m==1)])

        V = ZZ**self.nrows()
        v = V.random_element()
        num_iterates = max([0] + [f.degree() - g.degree() for g, _ in F if g.degree() > 1]) + 1

        S = []

        F.sort()
        for i in range(len(F)):
            g, m = F[i]

            if g.degree() == 1:
                # Just use kernel -- much easier.
                B = A.__copy__()
                for k from 0 <= k < A.nrows():
                    B[k,k] += g[0]
                if m > 1 and not is_diagonalizable:
                    B = B**m
                B = B.change_ring(QQ)
                W = B.kernel(algorithm = kernel_algorithm, **kwds)
                E.append((W, m==1))
                continue

            # General case, i.e., deg(g) > 1:
            W = None
            tries = m
            while True:

                # Compute the complementary factor of the charpoly.
                h = f // (g**m)
                v = h.list()

                while len(S) < tries:
                    t = verbose('%s-spinning %s-th random vector' % (num_iterates, len(S)),
                                level=2, caller_name='rational decomp')
                    S.append(A.iterates(V.random_element(x=-10,y=10), num_iterates))
                    verbose('done spinning', level=2, t=t, caller_name='rational decomp')

                for j in range(0 if W is None else W.nrows() // g.degree(), len(S)):
                    # Compute one element of the kernel of g(A)**m.
                    t = verbose('compute element of kernel of g(A), for g of degree %s' % g.degree(), level=2,
                                caller_name='rational decomp')
                    w = S[j].linear_combination_of_rows(h.list())
                    t = verbose('done computing element of kernel of g(A)', t=t,level=2, caller_name='rational decomp')

                    # Get the rest of the kernel.
                    t = verbose('fill out rest of kernel',level=2, caller_name='rational decomp')
                    if W is None:
                        W = A.iterates(w, g.degree())
                    else:
                        W = W.stack(A.iterates(w, g.degree()))
                    t = verbose('finished filling out more of kernel',level=2, t=t, caller_name='rational decomp')

                if W.rank() == m * g.degree():
                    W = W.change_ring(QQ)
                    t = verbose('now computing row space', level=2, caller_name='rational decomp')
                    W.echelonize(algorithm = echelon_algorithm, **kwds)
                    E.append((W.row_space(), m==1))
                    verbose('computed row space', level=2,t=t, caller_name='rational decomp')
                    break
                else:
                    verbose('we have not yet generated all the kernel (rank so far=%s, target rank=%s)' % (W.rank(), m*g.degree()),
                            level=2, caller_name='rational decomp')
                    tries += 1
                    if tries > 5*m:
                        raise RuntimeError("likely bug in decomposition")
                # end if
            # end while
        # end for
        return decomp_seq(E)

#      def simple_decomposition(self, echelon_algorithm='default', **kwds):
#          """
#          Returns the decomposition of the free module on which this
#          matrix A acts from the right (i.e., the action is x goes to x
#          A), as a direct sum of simple modules.

#          NOTE: self *must* be diagonalizable.

#          INPUT:
#              self -- a square matrix that is assumed to be diagonalizable
#              echelon_algorithm -- 'default'
#                                   'multimodular' -- use this if the answers
#                                   have small height
#              **kwds -- passed on to echelon function.

#          IMPORTANT NOTE:
#          If you expect that the subspaces in the answer are spanned by vectors
#          with small height coordinates, use algorithm='multimodular' and
#          height_guess=1; this is potentially much faster than the default.
#          If you know for a fact the answer will be very small, use
#             algorithm='multimodular', height_guess=bound on height, proof=False

#          OUTPUT:
#              Sequence -- list of tuples (V,g), where V is a subspace
#                          and an irreducible polynomial g, which is the
#                          charpoly (=minpoly) of self acting on V.
#          """
#          cdef Py_ssize_t k

#          if not self.is_square():
#              raise ArithmeticError("self must be a square matrix")

#          if self.nrows() == 0:
#              return decomp_seq([])

#          A, _ = self._clear_denom()

#          f = A.charpoly('x')
#          E = decomp_seq([])

#          t = verbose('factoring the characteristic polynomial', level=2, caller_name='simple decomp')
#          F = f.factor()
#          G = [g for g, _ in F]
#          minpoly = prod(G)
#          squarefree_degree = sum([g.degree() for g in G])
#          verbose('done factoring', t=t, level=2, caller_name='simple decomp')

#          V = ZZ**self.nrows()
#          v = V.random_element()
#          num_iterates = max([squarefree_degree - g.degree() for g in G]) + 1

#          S = [ ]

#          F.sort()
#          for i in range(len(F)):
#              g, m = F[i]

#              if g.degree() == 1:
#                  # Just use kernel -- much easier.
#                  B = A.__copy__()
#                  for k from 0 <= k < A.nrows():
#                      B[k,k] += g[0]
#                  if m > 1 and not is_diagonalizable:
#                      B = B**m
#                  W = B.change_ring(QQ).kernel()
#                  for b in W.basis():
#                      E.append((W.span(b), g))
#                  continue

#              # General case, i.e., deg(g) > 1:
#              W = None
#              while True:

#                  # Compute the complementary factor of the charpoly.
#                  h = minpoly // g
#                  v = h.list()

#                  while len(S) < m:
#                      t = verbose('%s-spinning %s-th random vector'%(num_iterates, len(S)),
#                                  level=2, caller_name='simple decomp')
#                      S.append(A.iterates(V.random_element(x=-10,y=10), num_iterates))
#                      verbose('done spinning', level=2, t=t, caller_name='simple decomp')

#                  for j in range(len(S)):
#                      # Compute one element of the kernel of g(A).
#                      t = verbose('compute element of kernel of g(A), for g of degree %s'%g.degree(),level=2,
#                              caller_name='simple decomp')
#                      w = S[j].linear_combination_of_rows(h.list())
#                      t = verbose('done computing element of kernel of g(A)', t=t,level=2, caller_name='simple decomp')

#                      # Get the rest of the kernel.
#                      t = verbose('fill out rest of kernel',level=2, caller_name='simple decomp')
#                      if W is None:
#                          W = A.iterates(w, g.degree())
#                      else:
#                          W = W.stack(A.iterates(w, g.degree()))
#                      t = verbose('finished filling out more of kernel',level=2, t=t, caller_name='simple decomp')

#                  if W.rank() == m * g.degree():
#                      W = W.change_ring(QQ)
#                      t = verbose('now computing row space', level=2, caller_name='simple decomp')
#                      W.echelonize(algorithm = echelon_algorithm, **kwds)
#                      E.append((W.row_space(), m==1))
#                      verbose('computed row space', level=2,t=t, caller_name='simple decomp')
#                      break
#                  else:
#                      verbose('we have not yet generated all the kernel (rank so far=%s, target rank=%s)'%(
#                          W.rank(), m*g.degree()), level=2, caller_name='simple decomp')
#                      j += 1
#                      if j > 3*m:
#                          raise RuntimeError("likely bug in decomposition")
#                  # end if
#              #end while
#          #end for
#          return E

    def _lift_crt_rr(self, res, mm):
        cdef Integer m
        cdef Matrix_integer_dense ZA
        cdef Matrix_rational_dense QA
        cdef Py_ssize_t i, j
        cdef mpz_t tmp
        cdef mpq_t tmp2
        mpz_init(tmp)
        mpq_init(tmp2)
        ZA = _lift_crt(res, mm)
        QA = Matrix_rational_dense.__new__(Matrix_rational_dense, self.parent(), None, None, None)
        m = mm.prod()
        for i in range(ZA._nrows):
            for j in range(ZA._ncols):
                fmpz_get_mpz(tmp, fmpz_mat_entry(ZA._matrix,i,j))
                mpq_rational_reconstruction(tmp2, tmp, m.value)
                fmpq_set_mpq(fmpq_mat_entry(QA._matrix, i, j), tmp2)
        mpz_clear(tmp)
        mpq_clear(tmp2)
        return QA

    def randomize(self, density=1, num_bound=2, den_bound=2,
                  distribution=None, nonzero=False):
        """
        Randomize ``density`` proportion of the entries of this matrix, leaving
        the rest unchanged.

        If ``x`` and ``y`` are given, randomized entries of this matrix have
        numerators and denominators bounded by ``x`` and ``y`` and have
        density 1.

        INPUT:

        - ``density`` -- number between 0 and 1 (default: 1)

        - ``num_bound`` -- numerator bound (default: 2)

        - ``den_bound`` -- denominator bound (default: 2)

        - ``distribution`` -- ``None`` or '1/n' (default: ``None``); if '1/n'
          then ``num_bound``, ``den_bound`` are ignored and numbers are chosen
          using the GMP function ``mpq_randomize_entry_recip_uniform``

        OUTPUT: none; the matrix is modified in-space

        EXAMPLES:

        The default distribution::

            sage: from collections import defaultdict
            sage: total_count = 0
            sage: dic = defaultdict(Integer)
            sage: def add_samples(distribution=None):
            ....:     global dic, total_count
            ....:     for _ in range(100):
            ....:         A = Matrix(QQ, 2, 4, 0)
            ....:         A.randomize(distribution=distribution)
            ....:         for a in A.list():
            ....:             dic[a] += 1
            ....:             total_count += 1.0

            sage: expected = {-2: 1/9, -1: 3/18, -1/2: 1/18, 0: 3/9,
            ....:             1/2: 1/18, 1: 3/18, 2: 1/9}
            sage: add_samples()
            sage: while not all(abs(dic[a]/total_count - expected[a]) < 0.001 for a in dic):
            ....:     add_samples()

        The distribution ``'1/n'``::

            sage: def mpq_randomize_entry_recip_uniform():
            ....:     r = 2*random() - 1
            ....:     if r == 0: r = 1
            ....:     num = int(4/(5*r))
            ....:     r = random()
            ....:     if r == 0: r = 1
            ....:     den = int(1/random())
            ....:     return Integer(num)/Integer(den)

            sage: total_count = 0
            sage: dic = defaultdict(Integer)
            sage: dic2 = defaultdict(Integer)
            sage: add_samples('1/n')
            sage: for _ in range(8):
            ....:     dic2[mpq_randomize_entry_recip_uniform()] += 1
            sage: while not all(abs(dic[a] - dic2[a])/total_count < 0.005 for a in dic):
            ....:     add_samples('1/n')
            ....:     for _ in range(800):
            ....:         dic2[mpq_randomize_entry_recip_uniform()] += 1

        The default can be used to obtain matrices of different rank::

            sage: ranks = [False]*11
            sage: while not all(ranks):
            ....:     for dens in (0.05, 0.1, 0.2, 0.5):
            ....:         A = Matrix(QQ, 10, 10, 0)
            ....:         A.randomize(dens)
            ....:         ranks[A.rank()] = True

        The default density is `6/9`::

            sage: def add_sample(density, num_rows, num_cols):
            ....:     global density_sum, total_count
            ....:     total_count += 1.0
            ....:     A = Matrix(QQ, num_rows, num_cols, 0)
            ....:     A.randomize(density)
            ....:     density_sum += float(A.density())

            sage: density_sum = 0.0
            sage: total_count = 0.0
            sage: expected_density = 6/9
            sage: add_sample(1.0, 100, 100)
            sage: while abs(density_sum/total_count - expected_density) > 0.001:
            ....:     add_sample(1.0, 100, 100)

        The modified density depends on the number of columns::

            sage: density_sum = 0.0
            sage: total_count = 0.0
            sage: expected_density = 6/9*0.5
            sage: add_sample(0.5, 100, 2)
            sage: while abs(density_sum/total_count - expected_density) > 0.001:
            ....:     add_sample(0.5, 100, 2)

            sage: density_sum = 0.0
            sage: total_count = 0.0
            sage: expected_density = 6/9*(1.0 - (99/100)^50)
            sage: expected_density
            0.263...

            sage: add_sample(0.5, 100, 100)
            sage: while abs(density_sum/total_count - expected_density) > 0.001:
            ....:     add_sample(0.5, 100, 100)

        Modifying the bounds for numerator and denominator::

            sage: num_dic = defaultdict(Integer)
            sage: den_dic = defaultdict(Integer)
            sage: while not (all(num_dic[i] for i in range(-200, 201))
            ....:            and all(den_dic[i] for i in range(1, 101))):
            ....:     a = matrix(QQ, 2, 4)
            ....:     a.randomize(num_bound=200, den_bound=100)
            ....:     for q in a.list():
            ....:         num_dic[q.numerator()] += 1
            ....:         den_dic[q.denominator()] += 1
            sage: len(num_dic)
            401
            sage: len(den_dic)
            100

        TESTS:

        Check that the option ``nonzero`` is meaningful (:issue:`22970`)::

            sage: a = matrix(QQ, 10, 10, 1)
            sage: b = a.__copy__()
            sage: b.randomize(nonzero=True)
            sage: a == b
            False
            sage: any(b[i,j].is_zero() for i in range(10) for j in range(10))
            False

        Check that :issue:`34103` is fixed::

            sage: a = matrix(QQ, 10, 10, 1)
            sage: a.randomize(nonzero=True, distribution='1/n')
            sage: bool(a)
            True
        """
        density = float(density)
        if density <= 0.0:
            return

        self.check_mutability()
        self.clear_cache()

        cdef Integer B, C
        cdef Py_ssize_t i, j, k, num_per_row
        cdef randstate rstate
        cdef mpq_t tmp

        B = Integer(num_bound + 1)
        C = Integer(den_bound + 1)

        mpq_init(tmp)

        if not nonzero:
            if density >= 1.0:
                if distribution == "1/n":
                    sig_on()
                    for i in range(self._nrows):
                        for j in range(self._ncols):
                            mpq_randomize_entry_recip_uniform(tmp)
                            fmpq_set_mpq(fmpq_mat_entry(self._matrix, i, j), tmp)
                    sig_off()
                elif mpz_cmp_si(C.value, 2):   # denom is > 1
                    sig_on()
                    for i in range(self._nrows):
                        for j in range(self._ncols):
                            mpq_randomize_entry(tmp, B.value, C.value)
                            fmpq_set_mpq(fmpq_mat_entry(self._matrix, i, j), tmp)
                    sig_off()
                else:
                    sig_on()
                    for i in range(self._nrows):
                        for j in range(self._ncols):
                            mpq_randomize_entry_as_int(tmp, B.value)
                            fmpq_set_mpq(fmpq_mat_entry(self._matrix, i, j), tmp)
                    sig_off()
            else:
                rstate = current_randstate()
                num_per_row = int(density * self._ncols)
                if distribution == "1/n":
                    sig_on()
                    for i in range(self._nrows):
                        for j in range(num_per_row):
                            k = rstate.c_random() % self._ncols
                            mpq_randomize_entry_recip_uniform(tmp)
                            fmpq_set_mpq(fmpq_mat_entry(self._matrix, i, k), tmp)
                    sig_off()
                elif mpz_cmp_si(C.value, 2):   # denom is > 1
                    sig_on()
                    for i in range(self._nrows):
                        for j in range(num_per_row):
                            k = rstate.c_random() % self._ncols
                            mpq_randomize_entry(tmp, B.value, C.value)
                            fmpq_set_mpq(fmpq_mat_entry(self._matrix, i, k), tmp)
                    sig_off()
                else:
                    sig_on()
                    for i in range(self._nrows):
                        for j in range(num_per_row):
                            k = rstate.c_random() % self._ncols
                            mpq_randomize_entry_as_int(tmp, B.value)
                            fmpq_set_mpq(fmpq_mat_entry(self._matrix, i, k), tmp)
                    sig_off()
        else:
            if density >= 1.0:
                if distribution == "1/n":
                    sig_on()
                    for i in range(self._nrows):
                        for j in range(self._ncols):
                            mpq_randomize_entry_recip_uniform_nonzero(tmp)
                            fmpq_set_mpq(fmpq_mat_entry(self._matrix, i, j), tmp)
                    sig_off()
                elif mpz_cmp_si(C.value, 2):   # denom is > 1
                    sig_on()
                    for i in range(self._nrows):
                        for j in range(self._ncols):
                            mpq_randomize_entry_nonzero(tmp, B.value, C.value)
                            fmpq_set_mpq(fmpq_mat_entry(self._matrix, i, j), tmp)
                    sig_off()
                else:
                    sig_on()
                    for i in range(self._nrows):
                        for j in range(self._ncols):
                            mpq_randomize_entry_as_int_nonzero(tmp, B.value)
                            fmpq_set_mpq(fmpq_mat_entry(self._matrix, i, j), tmp)
                    sig_off()
            else:
                rstate = current_randstate()
                num_per_row = int(density * self._ncols)
                if distribution == "1/n":
                    sig_on()
                    for i in range(self._nrows):
                        for j in range(num_per_row):
                            k = rstate.c_random() % self._ncols
                            mpq_randomize_entry_recip_uniform_nonzero(tmp)
                            fmpq_set_mpq(fmpq_mat_entry(self._matrix, i, k), tmp)
                    sig_off()
                elif mpz_cmp_si(C.value, 2):   # denom is > 1
                    sig_on()
                    for i in range(self._nrows):
                        for j in range(num_per_row):
                            k = rstate.c_random() % self._ncols
                            mpq_randomize_entry_nonzero(tmp, B.value, C.value)
                            fmpq_set_mpq(fmpq_mat_entry(self._matrix, i, k), tmp)
                    sig_off()
                else:
                    sig_on()
                    for i in range(self._nrows):
                        for j in range(num_per_row):
                            k = rstate.c_random() % self._ncols
                            mpq_randomize_entry_as_int_nonzero(tmp, B.value)
                            fmpq_set_mpq(fmpq_mat_entry(self._matrix, i, k), tmp)
                    sig_off()

        mpq_clear(tmp)

    def rank(self, algorithm=None):
        """
        Return the rank of this matrix.

        INPUT:

        - ``algorithm`` -- an optional specification of an algorithm. One of

          - ``None`` -- (default) will use flint

          - ``'flint'`` -- uses the flint library

          - ``'pari'`` -- uses the PARI library

          - ``'integer'`` -- eliminate denominators and calls the rank function
            on the corresponding integer matrix

        EXAMPLES::

            sage: matrix(QQ,3,[1..9]).rank()
            2
            sage: matrix(QQ,100,[1..100^2]).rank()
            2

        TESTS::

            sage: for _ in range(100):
            ....:     dim = randint(0, 30)
            ....:     m = random_matrix(QQ, dim, num_bound=2, density=0.5)
            ....:     r_pari = m.rank('pari'); m._clear_cache()
            ....:     r_flint = m.rank('flint'); m._clear_cache()
            ....:     r_int = m.rank('integer'); m._clear_cache()
            ....:     assert r_pari == r_flint == r_int
        """
        r = self.fetch('rank')
        if r is not None:
            return r

        if algorithm is None:
            algorithm = "flint"

        if algorithm == "flint":
            self.echelon_form(algorithm='flint')
            return self.fetch('rank')
        elif algorithm == "pari":
            r = self._rank_pari()
        elif algorithm == "integer":
            A, _ = self._clear_denom()
            r = A.rank()
        else:
            raise ValueError("unknown algorithm %s" % algorithm)

        self.cache('rank', r)
        return r

    def transpose(self):
        """
        Return the transpose of ``self``, without changing ``self``.

        EXAMPLES:

        We create a matrix, compute its transpose, and note that the
        original matrix is not changed.

        ::

            sage: A = matrix(QQ, 2, 3, range(6))
            sage: type(A)
            <class 'sage.matrix.matrix_rational_dense.Matrix_rational_dense'>
            sage: B = A.transpose()
            sage: print(B)
            [0 3]
            [1 4]
            [2 5]
            sage: print(A)
            [0 1 2]
            [3 4 5]

        ``.T`` is a convenient shortcut for the transpose::

            sage: print(A.T)
            [0 3]
            [1 4]
            [2 5]

        ::

            sage: A.subdivide(None, 1); A
            [0|1 2]
            [3|4 5]
            sage: A.transpose()
            [0 3]
            [---]
            [1 4]
            [2 5]
        """
        cdef Matrix_rational_dense ans
        if self._nrows == self._ncols:
            parent = self._parent
        else:
            parent = self._parent.matrix_space(self._ncols, self._nrows)
        ans = Matrix_rational_dense.__new__(Matrix_rational_dense, parent, None, None, None)
        sig_on()
        fmpq_mat_transpose(ans._matrix, self._matrix)
        sig_off()

        if self._subdivisions is not None:
            row_divs, col_divs = self.subdivisions()
            ans.subdivide(col_divs, row_divs)
        return ans

    def antitranspose(self):
        """
        Return the antitranspose of ``self``, without changing ``self``.

        EXAMPLES::

            sage: A = matrix(QQ,2,3,range(6))
            sage: type(A)
            <class 'sage.matrix.matrix_rational_dense.Matrix_rational_dense'>
            sage: A.antitranspose()
            [5 2]
            [4 1]
            [3 0]
            sage: A
            [0 1 2]
            [3 4 5]

            sage: A.subdivide(1,2); A
            [0 1|2]
            [---+-]
            [3 4|5]
            sage: A.antitranspose()
            [5|2]
            [-+-]
            [4|1]
            [3|0]
        """
        if self._nrows == self._ncols:
            parent = self._parent
        else:
            parent = self._parent.matrix_space(self._ncols, self._nrows)

        cdef Matrix_rational_dense ans
        ans = Matrix_rational_dense.__new__(Matrix_rational_dense, parent, None, None, None)

        cdef Py_ssize_t i,j
        cdef Py_ssize_t ri,rj  # reversed i and j
        sig_on()
        ri = self._nrows
        for i in range(self._nrows):
            rj = self._ncols
            ri = ri - 1
            for j in range(self._ncols):
                rj = rj - 1
                fmpq_set(fmpq_mat_entry(ans._matrix, rj, ri),
                         fmpq_mat_entry(self._matrix, i, j))
        sig_off()

        if self._subdivisions is not None:
            row_divs, col_divs = self.subdivisions()
            ans.subdivide([self._ncols - t for t in reversed(col_divs)],
                          [self._nrows - t for t in reversed(row_divs)])
        return ans

    def set_row_to_multiple_of_row(self, Py_ssize_t i, Py_ssize_t j, s):
        """
        Set row i equal to s times row j.

        EXAMPLES::

            sage: a = matrix(QQ,2,3,range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a.set_row_to_multiple_of_row(1,0,-3)
            sage: a
            [ 0  1  2]
            [ 0 -3 -6]
        """
        self.check_row_bounds_and_mutability(i, j)
        cdef Py_ssize_t k
        cdef fmpq_t ss
        fmpq_init(ss)
        fmpq_set_mpq(ss, (<Rational> Rational(s)).value)
        for k in range(self._ncols):
            fmpq_mul(fmpq_mat_entry(self._matrix, i, k),
                     fmpq_mat_entry(self._matrix, j, k),
                     ss)
        fmpq_clear(ss)

    def _set_row_to_negative_of_row_of_A_using_subset_of_columns(self, Py_ssize_t i, Matrix A,
                                                                 Py_ssize_t r, cols,
                                                                 cols_index=None):
        """
        Set row i of ``self`` to -(row r of A), but where we only take the
        given column positions in that row of A. We do not zero out the
        other entries of ``self``'s row i either.

        .. NOTE::

            This function exists just because it is useful for modular symbols presentations.

        INPUT:

        - ``i`` -- integer, index into the rows of self

        - ``A`` -- a matrix

        - ``r`` -- integer, index into rows of A

        - ``cols`` -- a *sorted* list of integers

        EXAMPLES::

            sage: a = matrix(QQ,2,3,range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a._set_row_to_negative_of_row_of_A_using_subset_of_columns(0,a,1,[1,2])
            sage: a
            [-4 -5  2]
            [ 3  4  5]
        """
        self.check_row_bounds_and_mutability(i,i)
        cdef Matrix_rational_dense _A
        if r < 0 or r >= A.nrows():
            raise IndexError("invalid row")
        cdef Py_ssize_t l = 0

        if not A.base_ring() == QQ:
            A = A.change_ring(QQ)
        if not A.is_dense():
            A = A.dense_matrix()

        _A = A
        for k in cols:
            entry = fmpq_mat_entry(self._matrix, i, l)
            fmpq_set(entry, fmpq_mat_entry(_A._matrix, r, k))
            fmpq_neg(entry, entry)
            l += 1

    def _add_col_j_of_A_to_col_i_of_self(self, Py_ssize_t i,
                                         Matrix_rational_dense A, Py_ssize_t j):
        """
        Unsafe technical function that very quickly adds the `j`-th column of
        A to the `i`-th column of ``self``.

        Does not check mutability.
        """
        if A._nrows != self._nrows:
            raise TypeError("nrows of self and A must be the same")
        cdef Py_ssize_t r
        for r in range(self._nrows):
            fmpq_add(fmpq_mat_entry(self._matrix, r, i),
                     fmpq_mat_entry(self._matrix, r, i),
                     fmpq_mat_entry(A._matrix, r, j))

    # ###############################################
    # Methods using PARI library
    # ###############################################

    def __pari__(self):
        """
        Return pari version of this matrix.

        EXAMPLES::

            sage: matrix(QQ,2,[1/5,-2/3,3/4,4/9]).__pari__()
            [1/5, -2/3; 3/4, 4/9]
        """
        return rational_matrix(self._matrix, False)

    def _det_pari(self, int flag=0):
        """
        Return the determinant of this matrix computed using pari.

        EXAMPLES::

            sage: matrix(QQ,3,[1..9])._det_pari()
            0
            sage: matrix(QQ,3,[1..9])._det_pari(1)
            0
            sage: matrix(QQ,3,[0]+[2..9])._det_pari()
            3
        """
        sig_on()
        cdef GEN d = det0(_new_GEN_from_fmpq_mat_t(self._matrix), flag)
        # now convert d to a Sage rational
        cdef Rational e = <Rational> Rational.__new__(Rational)
        INTFRAC_to_mpq(e.value, d)
        clear_stack()
        return e

    def _rank_pari(self):
        """
        Return the rank of this matrix computed using pari.

        EXAMPLES::

            sage: matrix(QQ,3,[1..9])._rank_pari()
            2
            sage: matrix(QQ, 0, 0)._rank_pari()
            0
        """
        sig_on()
        cdef long r = rank(_new_GEN_from_fmpq_mat_t(self._matrix))
        clear_stack()
        return r

    def _multiply_pari(self, Matrix_rational_dense right):
        """
        Return the product of ``self`` and ``right``, computed using PARI.

        EXAMPLES::

            sage: matrix(QQ,2,[1/5,-2/3,3/4,4/9])._multiply_pari(matrix(QQ,2,[1,2,3,4]))
            [  -9/5 -34/15]
            [ 25/12  59/18]

        We verify that 0 rows or columns works::

            sage: x = matrix(QQ,2,0); y = matrix(QQ,0,2); x*y
            [0 0]
            [0 0]
            sage: matrix(ZZ, 0, 0) * matrix(QQ, 0, 5)
            []
        """
        if self._ncols != right._nrows:
            raise ArithmeticError("self must be a square matrix")
        if not self._ncols*self._nrows or not right._ncols*right._nrows:
            # pari doesn't work in case of 0 rows or columns
            # This case is easy, since the answer must be the 0 matrix.
            return self.matrix_space(self._nrows, right._ncols).zero_matrix().__copy__()
        sig_on()
        cdef GEN M = gmul(_new_GEN_from_fmpq_mat_t(self._matrix),
                          _new_GEN_from_fmpq_mat_t(right._matrix))
        A = new_matrix_from_pari_GEN(self.matrix_space(self._nrows, right._ncols), M)
        clear_stack()
        return A

    def _invert_pari(self):
        """
        Return the inverse of this matrix computed using PARI.

        EXAMPLES::

            sage: matrix(QQ,2,[1,2,3,4])._invert_pari()
            [  -2    1]
            [ 3/2 -1/2]
            sage: matrix(QQ,2,[1,2,2,4])._invert_pari()
            Traceback (most recent call last):
            ...
            PariError: impossible inverse in ginv: [1, 2; 2, 4]
        """
        if self._nrows != self._ncols:
            raise ValueError("self must be a square matrix")
        cdef GEN M, d

        sig_on()
        M = _new_GEN_from_fmpq_mat_t(self._matrix)
        d = ginv(M)

        # Convert matrix back to Sage.
        A = new_matrix_from_pari_GEN(self._parent, d)
        clear_stack()
        return A

    def row(self, Py_ssize_t i, from_list=False):
        """
        Return the `i`-th row of this matrix as a dense vector.

        INPUT:

        - ``i`` -- integer

        - ``from_list`` -- ignored

        EXAMPLES::

            sage: m = matrix(QQ, 2, [1/5, -2/3, 3/4, 4/9])
            sage: m.row(0)
            (1/5, -2/3)
            sage: m.row(1)
            (3/4, 4/9)
            sage: m.row(1, from_list=True)
            (3/4, 4/9)
            sage: m.row(-2)
            (1/5, -2/3)

            sage: m.row(2)
            Traceback (most recent call last):
            ...
            IndexError: row index out of range
            sage: m.row(-3)
            Traceback (most recent call last):
            ...
            IndexError: row index out of range
        """
        if i < 0:
            i = i + self._nrows
        if i < 0 or i >= self._nrows:
            raise IndexError("row index out of range")

        cdef Py_ssize_t j
        parent = self.row_ambient_module()
        cdef Vector_rational_dense v = Vector_rational_dense.__new__(Vector_rational_dense)
        v._init(self._ncols, parent)
        for j in range(self._ncols):
            fmpq_get_mpq(v._entries[j], fmpq_mat_entry(self._matrix, i, j))
        return v

    def column(self, Py_ssize_t i, from_list=False):
        """
        Return the `i`-th column of this matrix as a dense vector.

        INPUT:

        - ``i`` -- integer

        - ``from_list`` -- ignored

        EXAMPLES::

            sage: m = matrix(QQ, 3, 2, [1/5,-2/3,3/4,4/9,-1,0])
            sage: m.column(1)
            (-2/3, 4/9, 0)
            sage: m.column(1,from_list=True)
            (-2/3, 4/9, 0)
            sage: m.column(-1)
            (-2/3, 4/9, 0)
            sage: m.column(-2)
            (1/5, 3/4, -1)

            sage: m.column(2)
            Traceback (most recent call last):
            ...
            IndexError: column index out of range
            sage: m.column(-3)
            Traceback (most recent call last):
            ...
            IndexError: column index out of range
        """
        if i < 0:
            i += self._ncols
        if i < 0 or i >= self._ncols:
            raise IndexError("column index out of range")

        cdef Py_ssize_t j
        parent = self.column_ambient_module()
        cdef Vector_rational_dense v = Vector_rational_dense.__new__(Vector_rational_dense)
        v._init(self._nrows, parent)
        for j in range(self._nrows):
            fmpq_get_mpq(v._entries[j], fmpq_mat_entry(self._matrix, j, i))
        return v

    # ###############################################
    # LLL
    # ###############################################

    def BKZ(self, *args, **kwargs):
        """
        Return the result of running Block Korkin-Zolotarev reduction on
        ``self`` interpreted as a lattice.

        The arguments ``*args`` and ``**kwargs`` are passed onto
        :meth:`sage.matrix.matrix_integer_dense.Matrix_integer_dense.BKZ`,
        see there for more details.

        EXAMPLES::

            sage: A = Matrix(QQ, 3, 3, [1/n for n in range(1, 10)])
            sage: A.BKZ()
            [ 1/28 -1/40 -1/18]
            [ 1/28 -1/40  1/18]
            [-1/14 -1/40     0]

            sage: A = random_matrix(QQ, 10, 10)
            sage: d = lcm(a.denom() for a in A.list())
            sage: A.BKZ() == (A * d).change_ring(ZZ).BKZ() / d
            True
        """
        A, d = self._clear_denom()
        return A.BKZ(*args, **kwargs) / d

    def LLL(self, *args, **kwargs):
        """
        Return an LLL reduced or approximated LLL reduced lattice for
        ``self`` interpreted as a lattice.

        The arguments ``*args`` and ``**kwargs`` are passed onto
        :meth:`sage.matrix.matrix_integer_dense.Matrix_integer_dense.LLL`,
        see there for more details.

        EXAMPLES::

            sage: A = Matrix(QQ, 3, 3, [1/n for n in range(1, 10)])
            sage: A.LLL()
            [ 1/28 -1/40 -1/18]
            [ 1/28 -1/40  1/18]
            [    0 -3/40     0]
            sage: L, U = A.LLL(transformation=True)
            sage: U * A == L
            True

            sage: A = random_matrix(QQ, 10, 10)
            sage: d = lcm(a.denom() for a in A.list())
            sage: A.LLL() == (A * d).change_ring(ZZ).LLL() / d
            True
        """
        A, d = self._clear_denom()
        if kwargs.get('transformation', False):
            L, U = A.LLL(*args, **kwargs)
            return L / d, U
        return A.LLL(*args, **kwargs) / d

    def is_LLL_reduced(self, delta=None, eta=None):
        r"""
        Return ``True`` if this lattice is `(\delta, \eta)`-LLL reduced.
        For a definition of LLL reduction, see
        :meth:`sage.matrix.matrix_integer_dense.Matrix_integer_dense.LLL`.

        EXAMPLES::

            sage: A = random_matrix(QQ, 10, 10)
            sage: L = A.LLL()
            sage: A.is_LLL_reduced()
            False
            sage: L.is_LLL_reduced()
            True
        """
        A, _ = self._clear_denom()
        return A.is_LLL_reduced(delta, eta)


cdef new_matrix_from_pari_GEN(parent, GEN d):
    """
    Given a PARI GEN with ``t_INT`` or ``t_FRAC entries, create a
    :class:`Matrix_rational_dense` from it.

    EXAMPLES::

        sage: matrix(QQ,2,[1..4])._multiply_pari(matrix(QQ,2,[2..5]))       # indirect doctest
        [10 13]
        [22 29]
    """
    cdef Py_ssize_t i, j
    cdef Matrix_rational_dense B = Matrix_rational_dense.__new__(
        Matrix_rational_dense, parent, None, None, None)
    cdef mpq_t tmp
    mpq_init(tmp)
    for i in range(B._nrows):
        for j in range(B._ncols):
            INTFRAC_to_mpq(tmp, gcoeff(d, i+1, j+1))
            fmpq_set_mpq(fmpq_mat_entry(B._matrix, i, j), tmp)
    mpq_clear(tmp)
    return B
