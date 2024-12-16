# distutils: extra_compile_args = NTL_CFLAGS M4RI_CFLAGS
# distutils: libraries = iml NTL_LIBRARIES gmp m CBLAS_LIBRARIES
# distutils: library_dirs = NTL_LIBDIR CBLAS_LIBDIR
# distutils: extra_link_args = NTL_LIBEXTRA
# distutils: include_dirs = NTL_INCDIR M4RI_INCDIR CBLAS_INCDIR
"""
Dense matrices over the integer ring

AUTHORS:

- William Stein

- Robert Bradshaw

- Marc Masdeu (August 2014). Implemented using FLINT, see :issue:`16803`.

- Jeroen Demeyer (October 2014): lots of fixes, see :issue:`17090` and
  :issue:`17094`.

- Vincent Delecroix (February 2015): make it faster, see :issue:`17822`.

- Vincent Delecroix (May 2017): removed duplication of entries and
  cleaner linbox interface

EXAMPLES::

    sage: a = matrix(ZZ, 3,3, range(9)); a
    [0 1 2]
    [3 4 5]
    [6 7 8]
    sage: a.det()
    0
    sage: a[0,0] = 10; a.det()
    -30
    sage: a.charpoly()
    x^3 - 22*x^2 + 102*x + 30
    sage: b = -3*a
    sage: a == b
    False
    sage: b < a
    True

TESTS::

    sage: a = matrix(ZZ,2,range(4), sparse=False)
    sage: TestSuite(a).run()
    sage: Matrix(ZZ,0,0).inverse()
    []
"""

#*****************************************************************************
#       Copyright (C) 2006,2007 William Stein
#       Copyright (C) 2014 Marc Masdeu
#       Copyright (C) 2014 Jeroen Demeyer
#       Copyright (C) 2015,2016,2017 Vincent Delecroix
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from libc.stdint cimport int64_t
from libc.string cimport strcpy, strlen

from sage.cpython.string cimport char_to_str, str_to_bytes
from sage.ext.stdsage cimport PY_NEW
from cysignals.signals cimport sig_check, sig_on, sig_str, sig_off
from cysignals.memory cimport sig_malloc, sig_free, check_allocarray

from sage.libs.gmp.mpz cimport *

from sage.modules.vector_integer_dense cimport Vector_integer_dense

from sage.misc.timing import cputime
from sage.misc.verbose import verbose, get_verbose

from sage.arith.misc import previous_prime
from sage.arith.long cimport integer_check_long_py
from sage.arith.power cimport generic_power
from sage.structure.element cimport Element
from sage.structure.proof.proof import get_flag as get_proof_flag
from sage.structure.richcmp cimport rich_to_bool
from sage.misc.randstate cimport randstate, current_randstate
from sage.misc.superseded import deprecated_function_alias
from sage.matrix.args cimport SparseEntry, MatrixArgs_init

#########################################################
# PARI C library
from cypari2.gen cimport Gen
from cypari2.stack cimport clear_stack, new_gen
from cypari2.paridecl cimport *
from sage.libs.pari import pari
from sage.libs.pari.convert_gmp cimport INT_to_mpz
from sage.libs.pari.convert_flint cimport (_new_GEN_from_fmpz_mat_t,
           _new_GEN_from_fmpz_mat_t_rotate90, integer_matrix)
from sage.libs.pari.convert_sage_matrix import gen_to_sage_matrix
#########################################################

from sage.arith.multi_modular cimport MultiModularBasis

from sage.libs.flint.fmpz cimport *
from sage.libs.flint.fmpz_mat cimport *

from sage.rings.integer cimport Integer
from sage.rings.rational_field import QQ
from sage.rings.real_double import RDF
from sage.rings.integer_ring import ZZ, IntegerRing_class
from sage.rings.integer_ring cimport IntegerRing_class
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.polynomial_integer_dense_flint cimport Polynomial_integer_dense_flint
from sage.structure.element cimport Element, Vector
from sage.structure.element import Vector

from sage.matrix.matrix_modn_dense_float cimport Matrix_modn_dense_template
from sage.matrix.matrix_modn_dense_float cimport Matrix_modn_dense_float
from sage.matrix.matrix_modn_dense_double cimport Matrix_modn_dense_double

from sage.matrix.matrix_mod2_dense import Matrix_mod2_dense
from sage.matrix.matrix_mod2_dense cimport Matrix_mod2_dense
from sage.rings.finite_rings.finite_field_constructor import GF


from sage.matrix.matrix2 import decomp_seq

from sage.matrix.matrix cimport Matrix

cimport sage.structure.element

import sage.matrix.matrix_space as matrix_space

################
# Used for modular HNF
from sage.rings.fast_arith cimport arith_int
cdef arith_int ai = arith_int()

######### linbox interface ##########
from sage.libs.linbox.linbox_flint_interface cimport *

########## iml -- integer matrix library ###########
from sage.libs.iml cimport *

fplll_fp_map = {None: None,
                'fp': 'double',
                'ld': 'long double',
                'xd': 'dpe',
                'rr': 'mpfr'}


cdef inline mpz_t * fmpz_mat_to_mpz_array(fmpz_mat_t m) except? NULL:
    cdef mpz_t * entries = <mpz_t *>check_allocarray(fmpz_mat_nrows(m), sizeof(mpz_t) * fmpz_mat_ncols(m))
    cdef size_t i, j
    cdef size_t k = 0
    sig_on()
    for i in range(fmpz_mat_nrows(m)):
        for j in range(fmpz_mat_ncols(m)):
            mpz_init(entries[k])
            fmpz_get_mpz(entries[k], fmpz_mat_entry(m, i, j))
            k += 1
    sig_off()
    return entries


cdef inline void mpz_array_clear(mpz_t * a, size_t length) noexcept:
    cdef size_t i
    for i in range(length):
        mpz_clear(a[i])
    sig_free(a)

cdef class Matrix_integer_dense(Matrix_dense):
    r"""
    Matrix over the integers, implemented using FLINT.

    On a 32-bit machine, they can have at most `2^{32}-1` rows or
    columns.  On a 64-bit machine, matrices can have at most
    `2^{64}-1` rows or columns.

    EXAMPLES::

        sage: a = MatrixSpace(ZZ,3)(2); a
        [2 0 0]
        [0 2 0]
        [0 0 2]
        sage: a = matrix(ZZ,1,3, [1,2,-3]); a
        [ 1  2 -3]
        sage: a = MatrixSpace(ZZ,2,4)(2); a
        Traceback (most recent call last):
        ...
        TypeError: nonzero scalar matrix must be square

    TESTS:

    Test hashing::

        sage: a = Matrix(ZZ, 2, [1,2,3,4])
        sage: hash(a)
        Traceback (most recent call last):
        ...
        TypeError: mutable matrices are unhashable
        sage: a.set_immutable()
        sage: hash(a)
        1846857684291126914  # 64-bit
        1591707266           # 32-bit
    """
    def __cinit__(self):
        """
        Create and allocate memory for the matrix. Does not actually
        initialize any of the memory.

        EXAMPLES::

            sage: from sage.matrix.matrix_integer_dense import Matrix_integer_dense
            sage: a = Matrix_integer_dense.__new__(Matrix_integer_dense, Mat(ZZ,3), 0,0,0)
            sage: type(a)
            <class 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'>

        TESTS::

            sage: Matrix(ZZ, sys.maxsize, sys.maxsize)
            Traceback (most recent call last):
            ...
            RuntimeError: FLINT exception
        """
        sig_str("FLINT exception")
        fmpz_mat_init(self._matrix, self._nrows, self._ncols)
        sig_off()

    def __dealloc__(self):
        """
        Frees all the memory allocated for this matrix.

        EXAMPLES::

            sage: a = Matrix(ZZ,2,[1,2,3,4])
            sage: del a
        """
        fmpz_mat_clear(self._matrix)

    def __init__(self, parent, entries=None, copy=None, bint coerce=True):
        r"""
        Initialize a dense matrix over the integers.

        INPUT:

        - ``parent`` -- a matrix space over ``ZZ``

        - ``entries`` -- see :func:`matrix`

        - ``copy`` -- ignored (for backwards compatibility)

        - ``coerce`` -- if ``False``, assume without checking that the
          entries are of type :class:`Integer`

        EXAMPLES:

        The __init__ function is called implicitly in each of the
        examples below to actually fill in the values of the matrix.

        We create a `2 \times 2` and a `1\times 4` matrix::

            sage: matrix(ZZ,2,2,range(4))
            [0 1]
            [2 3]
            sage: Matrix(ZZ,1,4,range(4))
            [0 1 2 3]

        If the number of columns isn't given, it is determined from the
        number of elements in the list.

        ::

            sage: matrix(ZZ,2,range(4))
            [0 1]
            [2 3]
            sage: matrix(ZZ,2,range(6))
            [0 1 2]
            [3 4 5]

        Another way to make a matrix is to create the space of matrices and
        coerce lists into it.

        ::

            sage: A = Mat(ZZ,2); A
            Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
            sage: A(range(4))
            [0 1]
            [2 3]

        Actually it is only necessary that the input can be coerced to a
        list, so the following also works::

            sage: v = reversed(range(4))
            sage: A(v)
            [3 2]
            [1 0]

        Matrices can have many rows or columns (in fact, on a 64-bit
        machine they could have up to `2^64-1` rows or columns)::

            sage: v = matrix(ZZ,1,10^5, range(10^5))
            sage: v.parent()
            Full MatrixSpace of 1 by 100000 dense matrices over Integer Ring
        """
        ma = MatrixArgs_init(parent, entries)
        cdef Integer z
        for t in ma.iter(coerce, True):
            se = <SparseEntry>t
            z = <Integer>se.entry
            fmpz_set_mpz(fmpz_mat_entry(self._matrix, se.i, se.j), z.value)

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, object x):
        """
        Set position i,j of this matrix to ``x``.

        The object ``x`` must be of type ``Integer``.

        INPUT:

        - ``i`` -- row

        - ``j`` -- column

        - ``x`` -- must be Integer! The value to set ``self[i,j]`` to.

        EXAMPLES::

            sage: a = matrix(ZZ,2,3, range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a[0,0] = 10
            sage: a
            [10  1  2]
            [ 3  4  5]
        """
        self.set_unsafe_mpz(i, j, (<Integer>x).value)

    cdef void set_unsafe_mpz(self, Py_ssize_t i, Py_ssize_t j, const mpz_t value) noexcept:
        """
        Set position i,j of this matrix to ``value``.

        INPUT:

        - ``i`` -- row

        - ``j`` -- column

        - ``value`` -- the value to set ``self[i,j]`` to; this will make a
          copy of ``value``

        EXAMPLES::

            sage: a = matrix(ZZ,2,3, range(6)); a
            [0 1 2]
            [3 4 5]
            sage: a[0,0] = 10
            sage: a
            [10  1  2]
            [ 3  4  5]
        """
        fmpz_set_mpz(fmpz_mat_entry(self._matrix,i,j), value)

    cdef void set_unsafe_int(self, Py_ssize_t i, Py_ssize_t j, int value) noexcept:
        """
        Set position i,j of this matrix to ``value``.
        """
        fmpz_set_si(fmpz_mat_entry(self._matrix,i,j), value)

    cdef void set_unsafe_double(self, Py_ssize_t i, Py_ssize_t j, double value) noexcept:
        """
        Set position i,j of this matrix to ``value``.
        """
        fmpz_set_d(fmpz_mat_entry(self._matrix,i,j), value)

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        """
        Return the (i, j) entry of ``self`` as a new Integer.

        .. WARNING::

           This is very unsafe; it assumes i and j are in the right
           range.

        EXAMPLES::

            sage: a = MatrixSpace(ZZ,3)(range(9)); a
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: a[1,2]
            5
            sage: a[4,7]
            Traceback (most recent call last):
            ...
            IndexError: matrix index out of range
            sage: a[-1,0]
            6
        """
        cdef Integer z = Integer.__new__(Integer)
        self.get_unsafe_mpz(i, j, z.value)
        return z

    cdef inline void get_unsafe_mpz(self, Py_ssize_t i, Py_ssize_t j, mpz_t value) noexcept:
        """
        Copy entry i,j of the matrix ``self`` to ``value``.

        .. WARNING::

           This is very unsafe; it assumes i and j are in the right
           range.

        EXAMPLES::

            sage: a = MatrixSpace(ZZ,3)(range(9)); a
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: a[1,2]
            5
            sage: a[4,7]
            Traceback (most recent call last):
            ...
            IndexError: matrix index out of range
            sage: a[-1,0]
            6
        """
        fmpz_get_mpz(value,fmpz_mat_entry(self._matrix, i, j))

    cdef inline int get_unsafe_int(self, Py_ssize_t i, Py_ssize_t j) noexcept:
        """
        Return the (i, j) entry of ``self`` as a new Integer.

        .. WARNING::

           This is very unsafe; it assumes i and j are in the right
           range.
        """
        return fmpz_get_si(fmpz_mat_entry(self._matrix, i, j))

    cdef inline double get_unsafe_double(self, Py_ssize_t i, Py_ssize_t j) noexcept:
        """
        Return the (i, j) entry of ``self`` as a new Integer.

        .. WARNING::

           This is very unsafe; it assumes i and j are in the right
           range.

        EXAMPLES::

            sage: a = MatrixSpace(ZZ,3)(range(9)); a
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: a[1,2]
            5
            sage: a[4,7]
            Traceback (most recent call last):
            ...
            IndexError: matrix index out of range
            sage: a[-1,0]
            6
        """
        return fmpz_get_d(fmpz_mat_entry(self._matrix, i, j))

    cdef bint get_is_zero_unsafe(self, Py_ssize_t i, Py_ssize_t j) except -1:
        """
        Return 1 if the entry (i, j) is zero, otherwise 0.

        .. WARNING::

           This is very unsafe; it assumes i and j are in the right
           range.
        """
        return fmpz_is_zero(fmpz_mat_entry(self._matrix, i,j))

    def _pickle(self):
        """
        EXAMPLES::

            sage: a = matrix(ZZ,2,3,[1,193,15,-2,3,0])
            sage: a._pickle() == (b'1 61 f -2 3 0', 0)
            True

            sage: S = ModularSymbols(250,4,sign=1).cuspidal_submodule().new_subspace().decomposition() # long time
            sage: S == loads(dumps(S)) # long time
            True
        """
        return self._pickle_version0(), 0

    cdef _pickle_version0(self):
        """
        EXAMPLES::

            sage: matrix(ZZ,1,3,[1,193,15])._pickle() == (b'1 61 f', 0)   # indirect doctest
            True
        """
        return str_to_bytes(self._export_as_string(32), 'ascii')

    cpdef _export_as_string(self, int base=10):
        """
        Return space separated string of the entries in this matrix, in the
        given base. This is optimized for speed.

        INPUT:

        - ``base`` -- integer <= 36; (default: 10)

        EXAMPLES::

            sage: m = matrix(ZZ,2,3,[1,2,-3,1,-2,-45])
            sage: m._export_as_string(10)
            '1 2 -3 1 -2 -45'
            sage: m._export_as_string(16)
            '1 2 -3 1 -2 -2d'
        """
        # TODO: *maybe* redo this to use mpz_import and mpz_export
        # from sec 5.14 of the GMP manual. ??
        cdef int i, j, len_so_far, m, n
        cdef char *s
        cdef char *t
        cdef char *tmp

        if self._nrows == 0 or self._ncols == 0:
            data = ''
        else:
            n = self._nrows*self._ncols*10
            s = <char*> sig_malloc(n * sizeof(char))
            t = s
            len_so_far = 0

            sig_on()
            for i from 0 <= i < self._nrows:
                for j from 0 <= j < self._ncols:
                    # mat_entry = fmpz_mat_entry(self._matrix,i,j)
                    m = fmpz_sizeinbase(fmpz_mat_entry(self._matrix,i,j), base)
                    if len_so_far + m + 2 >= n:
                        # copy to new string with double the size
                        n = 2*n + m + 1
                        tmp = <char*> sig_malloc(n * sizeof(char))
                        strcpy(tmp, s)
                        sig_free(s)
                        s = tmp
                        t = s + len_so_far
                    #endif
                    fmpz_get_str(t, base, fmpz_mat_entry(self._matrix,i,j))
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

    def _unpickle(self, data, int version):
        """
        TESTS::

            sage: b = matrix(ZZ,2,3, [0,0,0, 0, 0, 0])
            sage: s = b'1 61 f -2 3 0'
            sage: t = s.decode()
            sage: b._unpickle(s, 0) == b._unpickle(t, 0)
            True
            sage: b
            [  1 193  15]
            [ -2   3   0]
        """
        if version == 0:
            if isinstance(data, str):
                # old Py2 pickle: old "bytes" object reaches us as a
                # latin1-encoded string.
                data = data.encode('latin1')
            if isinstance(data, bytes):
                self._unpickle_version0(data)
            elif isinstance(data, list):
                self._unpickle_matrix_2x2_version0(data)
            else:
                raise RuntimeError("invalid pickle data")
        else:
            raise RuntimeError("unknown matrix version (=%s)"%version)

    cdef _unpickle_version0(self, data):
        cdef Py_ssize_t i, j, n, k
        data = data.split()
        n = self._nrows * self._ncols
        if len(data) != n:
            raise RuntimeError("invalid pickle data")
        k = 0
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                s = data[k]
                k += 1
                if fmpz_set_str(fmpz_mat_entry(self._matrix, i, j), s, 32):
                    raise RuntimeError("invalid pickle data")

    def _unpickle_matrix_2x2_version0(self, data):
        if len(data) != 4 or self._nrows != 2 or self._ncols != 2:
            raise RuntimeError("invalid pickle data")
        self.set_unsafe(0, 0, data[0])
        self.set_unsafe(0, 1, data[1])
        self.set_unsafe(1, 0, data[2])
        self.set_unsafe(1, 1, data[3])

    ########################################################################
    # LEVEL 1 helpers:
    #   These function support the implementation of the level 1 functionality.
    ########################################################################
    cdef Matrix_integer_dense _new(self, Py_ssize_t nrows, Py_ssize_t ncols):
        """
        Return a new matrix over the integers from given parent
        All memory is allocated for this matrix, but its
        entries have not yet been filled in.
        """
        if nrows == self._nrows and ncols == self._ncols:
            P = self._parent
        else:
            P = matrix_space.MatrixSpace(ZZ, nrows, ncols, sparse=False)
        cdef Matrix_integer_dense ans = Matrix_integer_dense.__new__(Matrix_integer_dense, P, None, None, None)
        return ans

    ########################################################################
    # LEVEL 2 functionality
    # x * cdef _add_
    # x * cdef _sub_
    # x * cdef _mul_
    # x * cpdef _richcmp_
    # x * __neg__
    # x * __invert__  -> SEE LEVEL 3 FUNCTIONALITIES
    # x * __copy__
    # x * _multiply_classical
    #   * _list -- list of underlying elements (need not be a copy)
    #   * _dict -- sparse dictionary of underlying elements (need not be a copy)
    ########################################################################

    def __copy__(self):
        r"""
        Return a new copy of this matrix.

        EXAMPLES::

            sage: a = matrix(ZZ,1,3, [1,2,-3]); a
            [ 1  2 -3]
            sage: b = a.__copy__(); b
            [ 1  2 -3]
            sage: b is a
            False
            sage: b == a
            True

            sage: M = MatrixSpace(ZZ,2,3)
            sage: m = M([1,2,3,3,2,1])
            sage: mc = m.__copy__()
            sage: mc == m and mc is not m
            True
        """
        cdef Matrix_integer_dense A
        A = self._new(self._nrows,self._ncols)

        sig_on()
        fmpz_mat_set(A._matrix,self._matrix)
        sig_off()
        if self._subdivisions is not None:
            A.subdivide(*self.subdivisions())
        return A

    def __bool__(self):
        r"""
        Test whether ``self`` is not the zero matrix.

        EXAMPLES::

            sage: a = MatrixSpace(ZZ, 2, 3)(range(6)); a
            [0 1 2]
            [3 4 5]
            sage: bool(a)
            True
            sage: bool(a - a)
            False

        ::

            sage: a = MatrixSpace(ZZ, 0, 3)()
            sage: bool(a)
            False
            sage: a = MatrixSpace(ZZ, 3, 0)()
            sage: bool(a)
            False
            sage: a = MatrixSpace(ZZ, 0, 0)()
            sage: bool(a)
            False
        """
        return not fmpz_mat_is_zero(self._matrix)

    def is_one(self):
        r"""
        Test whether ``self`` is the identity matrix.

        EXAMPLES::

            sage: matrix(2, [1,0,0,1]).is_one()
            True
            sage: matrix(2, [1,1,0,1]).is_one()
            False
            sage: matrix(2, 3, [1,0,0,0,1,0]).is_one()
            False
        """
        return self.is_square() and fmpz_mat_is_one(self._matrix)

    def _multiply_linbox(self, Matrix_integer_dense right):
        """
        Multiply matrices over ZZ using linbox.

        .. WARNING::

           This is very slow right now, i.e., linbox is very slow.

        EXAMPLES::

            sage: A = matrix(ZZ, 2, 3, range(6))
            sage: A * A.transpose()
            [ 5 14]
            [14 50]
            sage: A._multiply_linbox(A.transpose())
            [ 5 14]
            [14 50]

        TESTS:

        This fixes a bug found in :issue:`17094`::

            sage: A = identity_matrix(ZZ, 3)
            sage: A._multiply_linbox(A)
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        cdef Matrix_integer_dense ans
        cdef Matrix_integer_dense left = <Matrix_integer_dense> self

        ans = self._new(left._nrows, right._ncols)

        sig_on()
        linbox_fmpz_mat_mul(ans._matrix, left._matrix, right._matrix)
        sig_off()

        return ans

    def _multiply_classical(self, Matrix_integer_dense right):
        """
        EXAMPLES::

            sage: n = 3
            sage: a = MatrixSpace(ZZ,n,n)(range(n^2))
            sage: b = MatrixSpace(ZZ,n,n)(range(1, n^2 + 1))
            sage: a._multiply_classical(b)
            [ 18  21  24]
            [ 54  66  78]
            [ 90 111 132]

        TESTS::

            sage: for _ in range(100):
            ....:     nrows = randint(0, 10)
            ....:     nmid = randint(0, 10)
            ....:     ncols = randint(0, 10)
            ....:     m1 = random_matrix(ZZ, nrows, nmid)
            ....:     m2 = random_matrix(ZZ, nmid, ncols)
            ....:     ans_flint = m1 * m2
            ....:     ans_classical = m1._multiply_classical(m2)
            ....:     ans_linbox = m1._multiply_linbox(m2)
            ....:     if ans_flint != ans_classical:
            ....:         raise RuntimeError("ERROR\nm1=\n{}\nm2=\n{}\nans_flint=\n{}\nans_classical=\n{}".format(
            ....:                 m1.str(), m2.str(), ans_flint.str(), ans_classical.str()))
            ....:     if ans_flint != ans_linbox:
            ....:         raise RuntimeError("ERROR\nm1=\n{}\nm2=\n{}\nans_flint=\n{}\nans_linbox=\n{}".format(
            ....:                 m1.str(), m2.str(), ans_flint.str(), ans_linbox.str()))
        """
        if self._ncols != right._nrows:
            raise IndexError("Number of columns of self must equal number of rows of right.")

        cdef Py_ssize_t i, j, k, nr, nc, snc
        cdef object parent

        nr = self._nrows
        nc = right._ncols
        snc = self._ncols

        if self._nrows == right._nrows:
            # self acts on the space of right
            parent = right.parent()
        if self._ncols == right._ncols:
            # right acts on the space of self
            parent = self.parent()
        else:
            parent = self.matrix_space(nr, nc)

        cdef Matrix_integer_dense M, _right
        _right = right

        M = self._new(parent.nrows(),parent.ncols())

        cdef fmpz_t s
        fmpz_init(s)
        sig_on()
        for i from 0 <= i < nr:
            for j from 0 <= j < nc:
                fmpz_set_si(s,0)   # set s = 0
                for k from 0 <= k < snc:
                    fmpz_addmul(s, fmpz_mat_entry(self._matrix,i,k), fmpz_mat_entry(_right._matrix,k,j))
                fmpz_set(fmpz_mat_entry(M._matrix,i,j),s)
        sig_off()
        fmpz_clear(s)
        return M

    cdef sage.structure.element.Matrix _matrix_times_matrix_(self, sage.structure.element.Matrix right):
        cdef Matrix_integer_dense M

        if self._ncols != right._nrows:
            raise IndexError("Number of columns of self must equal number of rows of right.")

        M = self._new(self._nrows, right._ncols)

        sig_on()
        fmpz_mat_mul(M._matrix, self._matrix, (<Matrix_integer_dense>right)._matrix)
        sig_off()
        return M

    cpdef _lmul_(self, Element right):
        """
        EXAMPLES::

            sage: a = matrix(ZZ, 2, range(6))
            sage: 5 * a
            [ 0  5 10]
            [15 20 25]
        """
        cdef Integer x = Integer(right)
        cdef fmpz_t z
        cdef Matrix_integer_dense M = self._new(self._nrows, self._ncols)

        sig_on()
        fmpz_init_set_readonly(z, x.value)
        fmpz_mat_scalar_mul_fmpz(M._matrix, self._matrix, z)
        fmpz_clear_readonly(z)
        sig_off()
        return M

    cpdef _add_(self, right):
        """
        Add two dense matrices over ZZ.

        EXAMPLES::

            sage: a = MatrixSpace(ZZ,3)(range(9))
            sage: a+a
            [ 0  2  4]
            [ 6  8 10]
            [12 14 16]
            sage: b = MatrixSpace(ZZ,3)(range(9))
            sage: b.swap_rows(1,2)
            sage: a+b
            [ 0  2  4]
            [ 9 11 13]
            [ 9 11 13]
        """
        cdef Matrix_integer_dense M = self._new(self._nrows, self._ncols)

        sig_on()
        fmpz_mat_add(M._matrix,self._matrix,(<Matrix_integer_dense> right)._matrix)
        sig_off()
        return M

    cpdef _sub_(self, right):
        """
        Subtract two dense matrices over ZZ.

        EXAMPLES::

            sage: M = Mat(ZZ,3)
            sage: a = M(range(9)); b = M(reversed(range(9)))
            sage: a - b
            [-8 -6 -4]
            [-2  0  2]
            [ 4  6  8]
        """
        cdef Matrix_integer_dense M = self._new(self._nrows, self._ncols)

        sig_on()
        fmpz_mat_sub(M._matrix,self._matrix,(<Matrix_integer_dense> right)._matrix)
        sig_off()
        return M

    def __pow__(sself, n, dummy):
        r"""
        Return the ``n``-th power of this matrix.

        EXAMPLES::

            sage: M = MatrixSpace(ZZ,3)
            sage: m = M([1, 1, 1, 2, 1, 1, -3, -2, -1])
            sage: m ** 3
            [-3 -2 -1]
            [-3 -2  0]
            [ 2  1 -3]
            sage: m ** -2
            [ 2 -3 -1]
            [-4  4  1]
            [ 1  0  0]
            sage: M(range(9)) ** -1
            Traceback (most recent call last):
            ...
            ZeroDivisionError: matrix must be nonsingular

        TESTS::

            sage: m ** 3 == m ** 3r == (~m) ** (-3) == (~m) ** (-3r)
            True

        The following exponents do not fit in an unsigned long and the
        multiplication method fall back to the generic power implementation in
        :mod:`sage.structure.element`::

            sage: m = M.identity_matrix()
            sage: m ** (2**256)
            [1 0 0]
            [0 1 0]
            [0 0 1]
            sage: m ** (2r**256r)
            [1 0 0]
            [0 1 0]
            [0 0 1]

        In this case, the second argument to ``__pow__`` is a matrix,
        which should raise the correct error::

            sage: M = Matrix(2, 2, range(4))
            sage: None^M
            Traceback (most recent call last):
            ...
            TypeError: Cannot convert NoneType to sage.matrix.matrix_integer_dense.Matrix_integer_dense
            sage: M^M
            Traceback (most recent call last):
            ...
            NotImplementedError: the given exponent is not supported
        """
        cdef Matrix_integer_dense self = <Matrix_integer_dense?>sself

        if dummy is not None:
            raise ValueError
        if self._nrows != self._ncols:
            raise ArithmeticError("self must be a square matrix")

        cdef unsigned long e
        cdef long e_sgn
        cdef int err

        if integer_check_long_py(n, &e_sgn, &err) and not err:
            if e_sgn < 0:
                return (~self) ** (-n)
            e = <unsigned long>e_sgn
        else:
            if not isinstance(n, Integer):
                try:
                    n = Integer(n)
                except TypeError:
                    from sage.structure.element import Expression
                    if isinstance(n, Expression):
                        from sage.matrix.matrix2 import _matrix_power_symbolic
                        return _matrix_power_symbolic(self, n)
                    else:
                        raise NotImplementedError("the given exponent is not supported")
            if mpz_sgn((<Integer>n).value) < 0:
                return (~self) ** (-n)

            if mpz_fits_ulong_p((<Integer>n).value):
                e = mpz_get_ui((<Integer>n).value)
            else:
                # it is very likely that the following will never finish except
                # if self has only eigenvalues 0, 1 or -1.
                return generic_power(self, n)

        if e == 0:
            return self._parent.identity_matrix()
        if e == 1:
            return self

        cdef Matrix_integer_dense M = self._new(self._nrows, self._ncols)
        sig_on()
        fmpz_mat_pow(M._matrix, self._matrix, e)
        sig_off()
        return M

    def __neg__(self):
        r"""
        Return the negative of this matrix.

        TESTS::

            sage: a = matrix(ZZ,2,range(4))
            sage: a.__neg__()
            [ 0 -1]
            [-2 -3]
            sage: -a
            [ 0 -1]
            [-2 -3]
        """
        cdef Matrix_integer_dense M = self._new(self._nrows, self._ncols)
        sig_on()
        fmpz_mat_neg(M._matrix, self._matrix)
        sig_off()
        return M

    cpdef _richcmp_(self, right, int op):
        r"""
        Compare ``self`` with ``right``, examining entries in
        lexicographic (row major) ordering.

        EXAMPLES::

            sage: Matrix(ZZ, [[0, 10], [20, 30]]) == (Matrix(ZZ, [[0, 10], [20, 30]]))
            True
            sage: Matrix(ZZ, [[0, 10], [20, 30]]) < (Matrix(ZZ, [[0, 15], [20, 30]]))
            True
            sage: Matrix(ZZ, [[5, 10], [20, 30]]) > (Matrix(ZZ, [[0, 15], [20, 30]]))
            True
            sage: Matrix(ZZ, [[5, 10], [20, 30]]) <= (Matrix(ZZ, [[0, 10], [25, 30]]))
            False
        """
        cdef Py_ssize_t i, j
        cdef int k

        sig_on()
        for i in range(self._nrows):
            for j in range(self._ncols):
                k = fmpz_cmp(fmpz_mat_entry(self._matrix,i,j),
                             fmpz_mat_entry((<Matrix_integer_dense>right)._matrix,i,j))
                if k:
                    sig_off()
                    if k < 0:
                        return rich_to_bool(op, -1)
                    else:
                        return rich_to_bool(op, 1)
        sig_off()
        return rich_to_bool(op, 0)

    # TODO: Implement better
    cdef _vector_times_matrix_(self, Vector v):
        """
        Return the vector times matrix product.

        INPUT:

        - ``v`` -- a free module element

        OUTPUT: the vector times matrix product ``v*A``

        EXAMPLES::

            sage: B = matrix(ZZ,2, [1,2,3,4])
            sage: V = ZZ^2
            sage: w = V([-1,5])
            sage: w*B
            (14, 18)
        """
        cdef Vector_integer_dense w, ans
        cdef Py_ssize_t i, j
        cdef fmpz_t x
        cdef fmpz_t z

        M = self.row_ambient_module()
        w = <Vector_integer_dense> v
        ans = M.zero_vector()

        sig_on()
        fmpz_init(x)
        fmpz_init(z)
        for i from 0 <= i < self._ncols:
            fmpz_set_si(x, 0)
            for j from 0 <= j < self._nrows:
                fmpz_set_mpz(z,w._entries[j])
                fmpz_addmul(x, z, fmpz_mat_entry(self._matrix,j,i))
            fmpz_get_mpz(ans._entries[i], x)
        fmpz_clear(x)
        fmpz_clear(z)
        sig_off()
        return ans

    ########################################################################
    # LEVEL 3 functionality (Optional)
    #    * __deepcopy__
    #  x * __invert__
    #    * Matrix windows -- only if you need strassen for that base
    #    * Other functions (list them here):
    #    * Specialized echelon form
    ########################################################################

    def is_primitive(self):
        r"""
        Test whether the matrix is primitive.

        An integral matrix `A` is primitive if all its entries are nonnegative
        and for some positive integer `n` the matrix `A^n` has all its
        entries positive.

        EXAMPLES::

            sage: m = matrix(3, [1,1,0,0,0,1,1,0,0])
            sage: m.is_primitive()
            True
            sage: m**4
            [3 2 1]
            [1 1 1]
            [2 1 1]

            sage: m = matrix(4, [[1,1,0,0],[0,0,1,0],[0,0,0,1],[1,0,0,0]])
            sage: m.is_primitive()
            True
            sage: m**6
            [4 3 2 1]
            [1 1 1 1]
            [2 1 1 1]
            [3 2 1 1]

            sage: m = matrix(4, [[0,1,0,1],[1,0,1,0],[0,1,0,1],[1,0,1,0]])
            sage: m.is_primitive()
            False

        Testing extremal matrices::

            sage: def matrix1(d):
            ....:     m = matrix(d)
            ....:     m[0,0] = 1
            ....:     for i in range(d-1):
            ....:         m[i,i+1] = m[i+1,i] = 1
            ....:     return m
            sage: all(matrix1(d).is_primitive() for d in range(2,20))
            True

            sage: def matrix2(d):
            ....:     m = matrix(d)
            ....:     for i in range(d-1):
            ....:         m[i,i+1] = 1
            ....:     m[d-1,0] = m[d-1,1] = 1
            ....:     return m
            sage: all(matrix2(d).is_primitive() for d in range(2,20))
            True

        Non-primitive families::

            sage: def matrix3(d):
            ....:     m = matrix(d)
            ....:     m[0,0] = 1
            ....:     for i in range(d-1):
            ....:         m[i,i+1] = 1
            ....:     return m
            sage: any(matrix3(d).is_primitive() for d in range(2,20))
            False

        TESTS::

            sage: matrix(1, 2, [1,3]).is_primitive()
            Traceback (most recent call last):
            ...
            ValueError: not a square matrix
        """
        # NOTE: ideally, all computations should be done with 0/1 matrices
        # with the bitwise operations + = OR and * = AND. Instead we do
        # computation over integers and reset positive entries to 1 to
        # avoid coefficient blowup.

        cdef long dim = fmpz_mat_nrows(self._matrix)
        if dim != fmpz_mat_ncols(self._matrix):
            raise ValueError("not a square matrix")

        cdef int s, diag, zero
        cdef size_t i,j, N
        cdef fmpz_mat_t m

        try:
            fmpz_mat_init(m, dim, dim)

            # 1. check that self._matrix is nonnegative and set
            #    m as a 0/1 matrix
            zero = 0
            diag = 0
            for i in range(dim):
                for j in range(dim):
                    s = fmpz_sgn(fmpz_mat_entry(self._matrix, i, j))
                    if s == 0:
                        fmpz_zero(fmpz_mat_entry(m, i, j))
                        zero = 1
                    elif s == -1:
                        return False
                    elif s == 1:
                        fmpz_one(fmpz_mat_entry(m, i, j))
                        if i == j:
                            diag = 1
                    else:
                        raise RuntimeError

            # 2. is the matrix already positive?
            if not zero:
                return True

            # 3. try powers up to the given bounds
            #    2 * (dim - 1)  : if one diagonal entry is nonzero
            #    (dim - 1)^2 + 1: otherwise
            fmpz_mat_sqr(m, m)
            if diag:
                N = (dim - 1)
            else:
                N = ((dim - 1) * (dim - 1) + 1) >> 1

            while N:
                zero = False
                for i in range(dim):
                    for j in range(dim):
                        if fmpz_is_zero(fmpz_mat_entry(m, i, j)):
                            zero = True
                        else:
                            fmpz_one(fmpz_mat_entry(m, i, j))
                if zero:
                    fmpz_mat_sqr(m, m)
                    N >>= 1
                else:
                    return True

            for i in range(dim):
                for j in range(dim):
                    if fmpz_is_zero(fmpz_mat_entry(m, i, j)):
                        return False

            return True

        finally:
            fmpz_mat_clear(m)

    def _clear_denom(self):
        """
        INPUT:

        - ``self`` -- a matrix

        OUTPUT: ``self, 1``

        EXAMPLES::

            sage: a = matrix(ZZ,2,[1,2,3,4])
            sage: a._clear_denom()
            (
            [1 2]
            [3 4], 1
            )
        """
        return self, ZZ(1)

    def charpoly(self, var='x', algorithm=None):
        r"""
        .. NOTE::

            The characteristic polynomial is defined as `\det(xI-A)`.

        INPUT:

        - ``var`` -- a variable name

        - ``algorithm`` -- (default: ``'linbox'``) either ``'generic'``,
          ``'flint'`` or ``'linbox'``

        EXAMPLES::

            sage: A = matrix(ZZ,6, range(36))
            sage: f = A.charpoly(); f
            x^6 - 105*x^5 - 630*x^4
            sage: f(A) == 0
            True
            sage: g = A.charpoly(algorithm='flint')
            sage: f == g
            True
            sage: n=20; A = Mat(ZZ,n)(range(n^2))
            sage: A.charpoly()
            x^20 - 3990*x^19 - 266000*x^18
            sage: A.minpoly()
            x^3 - 3990*x^2 - 266000*x

        On non square matrices, this method raises an ArithmeticError::

            sage: matrix(ZZ, 2, 1).charpoly()
            Traceback (most recent call last):
            ...
            ArithmeticError: only valid for square matrix

        TESTS:

        The cached polynomial should be independent of the ``var``
        argument (:issue:`12292`). We check (indirectly) that the
        second call uses the cached value by noting that its result is
        not cached::

            sage: M = MatrixSpace(ZZ, 2)
            sage: A = M(range(0, 2^2))
            sage: type(A)
            <class 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'>
            sage: A.charpoly('x')
            x^2 - 3*x - 2
            sage: A.charpoly('y')
            y^2 - 3*y - 2
            sage: A._cache['charpoly_linbox']
            x^2 - 3*x - 2

        Test corner cases::

            sage: matrix([5]).charpoly('z', 'flint')
            z - 5
            sage: matrix([5]).charpoly('z', 'linbox')
            z - 5
            sage: matrix([5]).charpoly('z', 'generic')
            z - 5


            sage: matrix([]).charpoly('y', 'flint')
            1
            sage: matrix([]).charpoly('y', 'linbox')
            1
            sage: matrix([]).charpoly('y', 'generic')
            1

            sage: matrix([0]).charpoly('x', 'flint')
            x
            sage: matrix([0]).charpoly('x', 'linbox')
            x
            sage: matrix([0]).charpoly('x', 'generic')
            x

        Consistency on random inputs::

            sage: for _ in range(100):
            ....:     dim = randint(1, 20)
            ....:     m  = random_matrix(ZZ, dim)
            ....:     m._clear_cache(); ans_flint = m.charpoly(algorithm='flint')
            ....:     m._clear_cache(); ans_linbox = m.charpoly(algorithm='linbox')
            ....:     m._clear_cache(); ans_generic = m.charpoly(algorithm='generic')
            ....:     if ans_flint != ans_linbox or ans_flint != ans_generic:
            ....:         raise RuntimeError("ans_flint = {}, ans_linbox = {} and ans_generic = {} for\n{}".format(
            ....:                            ans_flint, ans_linbox, ans_generic, m.str()))
        """
        if self._nrows != self._ncols:
            raise ArithmeticError("only valid for square matrix")

        cdef Polynomial_integer_dense_flint g

        if algorithm is None:
            algorithm = 'linbox'

        cache_key = 'charpoly_%s' % algorithm
        g = self.fetch(cache_key)
        if g is not None:
            return g.change_variable_name(var)

        if algorithm == 'flint':
            g = (<Polynomial_integer_dense_flint> PolynomialRing(ZZ, names=var).gen())._new()
            sig_on()
            fmpz_mat_charpoly(g._poly, self._matrix)
            sig_off()
        elif algorithm == 'linbox':
            g = (<Polynomial_integer_dense_flint> PolynomialRing(ZZ, names=var).gen())._new()
            sig_on()
            linbox_fmpz_mat_charpoly(g._poly, self._matrix)
            sig_off()
        elif algorithm == 'generic':
            g = Matrix_dense.charpoly(self, var)
        else:
            raise ValueError("no algorithm '%s'"%algorithm)

        self.cache(cache_key, g)
        return g

    def minpoly(self, var='x', algorithm=None):
        r"""
        INPUT:

        - ``var`` -- a variable name

        - ``algorithm`` -- either ``'linbox'`` (default) or ``'generic'``

        EXAMPLES::

            sage: A = matrix(ZZ, 6, range(36))
            sage: A.minpoly()
            x^3 - 105*x^2 - 630*x

            sage: A = Mat(ZZ, 6)([k^2 for k in range(36)])
            sage: A.minpoly(algorithm='linbox')
            x^4 - 2695*x^3 - 257964*x^2 + 1693440*x
            sage: A.minpoly(algorithm='generic')
            x^4 - 2695*x^3 - 257964*x^2 + 1693440*x

        On non square matrices, this method raises an ArithmeticError::

            sage: matrix(ZZ, 2, 1).minpoly()
            Traceback (most recent call last):
            ...
            ArithmeticError: only valid for square matrix

        TESTS:

        Corner cases::

            sage: matrix([5]).minpoly('z', 'linbox')
            z - 5
            sage: matrix([5]).minpoly('z', 'generic')
            z - 5

            sage: matrix([]).minpoly('y', 'linbox')
            1
            sage: matrix([]).minpoly('y', 'generic')
            1

            sage: matrix(ZZ, 2).minpoly('x', 'linbox')
            x
            sage: matrix(ZZ, 2).minpoly('x', 'generic')
            x

        Consistency on random inputs::

            sage: for _ in range(100):
            ....:     dim = randint(1, 20)
            ....:     m  = random_matrix(ZZ, dim)
            ....:     m._clear_cache(); ans_generic = m.minpoly(algorithm='generic')
            ....:     m._clear_cache(); ans_linbox = m.minpoly(algorithm='linbox')
            ....:     if ans_generic != ans_linbox:
            ....:         raise RuntimeError("ans_generic = {} and ans_linbox = {} for\n{}".format(
            ....:                            ans_generic, ans_linbox, m.str()))
        """
        if self._nrows != self._ncols:
            raise ArithmeticError("only valid for square matrix")

        cdef Polynomial_integer_dense_flint g

        if algorithm is None:
            algorithm = 'linbox'

        key = 'minpoly_%s'%(algorithm)
        g = self.fetch(key)
        if g is not None:
            return g.change_variable_name(var)

        if algorithm == 'linbox':
            g = (<Polynomial_integer_dense_flint> PolynomialRing(ZZ, names=var).gen())._new()
            sig_on()
            linbox_fmpz_mat_minpoly(g._poly, self._matrix)
            sig_off()
        elif algorithm == 'generic':
            g = Matrix_dense.minpoly(self, var)
        else:
            raise ValueError("no algorithm '%s'"%algorithm)

        self.cache(key, g)
        return g

    def height(self):
        """
        Return the height of this matrix, i.e., the max absolute value of
        the entries of the matrix.

        OUTPUT: nonnegative integer

        EXAMPLES::

            sage: a = Mat(ZZ,3)(range(9))
            sage: a.height()
            8
            sage: a = Mat(ZZ,2,3)([-17,3,-389,15,-1,0]); a
            [ -17    3 -389]
            [  15   -1    0]
            sage: a.height()
            389
        """
        cdef Integer x = Integer.__new__(Integer)
        self.mpz_height(x.value)
        return x

    cdef int mpz_height(self, mpz_t height) except -1:
        """
        Used to compute the height of this matrix.

        INPUT:

        - ``height`` -- a GMP mpz_t which has been initialized

        OUTPUT: sets the value of height to the height of this matrix,
        i.e., the max absolute value of the entries of the matrix.
        """
        cdef fmpz_t x,h
        cdef Py_ssize_t i,j

        sig_on()
        fmpz_init(h)
        fmpz_init(x)
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                fmpz_abs(x, fmpz_mat_entry(self._matrix,i,j))
                if fmpz_cmp(h, x) < 0:
                    fmpz_set(h, x)
        fmpz_get_mpz(height,h)
        fmpz_clear(h)
        fmpz_clear(x)
        sig_off()
        return 0   # no error occurred.

    def _multiply_multi_modular(self, Matrix_integer_dense right):
        """
        Multiply this matrix by ``left`` using a multi modular algorithm.

        EXAMPLES::

            sage: M = Matrix(ZZ, 2, 3, range(5,11))
            sage: N = Matrix(ZZ, 3, 2, range(15,21))
            sage: M._multiply_multi_modular(N)
            [310 328]
            [463 490]
            sage: M._multiply_multi_modular(-N)
            [-310 -328]
            [-463 -490]
        """
        cdef Integer h
        cdef Matrix_integer_dense left = <Matrix_integer_dense>self
        cdef int i, k

        nr = left._nrows
        nc = right._ncols

        cdef Matrix_integer_dense result

        h = left.height() * right.height() * left.ncols()
        verbose('multiplying matrices of height %s and %s'%(left.height(),right.height()))
        mm = MultiModularBasis(h)
        res = left._reduce(mm)
        res_right = right._reduce(mm)
        k = len(mm)
        for i in range(k):  # yes, I could do this with zip, but to conserve memory...
            t = cputime()
            res[i] *= res_right[i]
            verbose('multiplied matrices modulo a prime (%s/%s)'%(i+1,k), t)
        result = left.new_matrix(nr,nc)
        _lift_crt(result, res, mm)  # changes result
        return result

    def _mod_int(self, modulus):
        """
        Reduce the integer matrix modulo a positive integer.

        EXAMPLES::

            sage: M = Matrix(ZZ, 2, [1,2,-2,3])
            sage: M._mod_int(2)
            [1 0]
            [0 1]
            sage: M._mod_int(1000000)
            [     1      2]
            [999998      3]
        """
        cdef mod_int c = modulus
        if int(c) != modulus:
            raise OverflowError
        else:
            return self._mod_int_c(modulus)

    cdef _mod_two(self):
        """
        TESTS:

        Check that bug discovered in :issue:`29839` is fixed::

            sage: M = Matrix(ZZ, [[0,1],[0,1]])
            sage: M._mod_int(2).transpose()
            [0 0]
            [1 1]
        """
        MS = matrix_space.MatrixSpace(GF(2), self._nrows, self._ncols)
        return Matrix_mod2_dense(MS, self, True, True)

    cdef _mod_int_c(self, mod_int p):
        from sage.matrix.matrix_modn_dense_float import MAX_MODULUS as MAX_MODULUS_FLOAT
        from sage.matrix.matrix_modn_dense_double import MAX_MODULUS as MAX_MODULUS_DOUBLE

        cdef Py_ssize_t i, j

        cdef float* res_row_f
        cdef Matrix_modn_dense_float res_f

        cdef double* res_row_d
        cdef Matrix_modn_dense_double res_d

        if p == 2:
            return self._mod_two()
        elif p < MAX_MODULUS_FLOAT:
            res_f = Matrix_modn_dense_float.__new__(Matrix_modn_dense_float,
                                                    matrix_space.MatrixSpace(IntegerModRing(p), self._nrows, self._ncols, sparse=False), None, None, None, zeroed_alloc=False)
            for i from 0 <= i < self._nrows:
                res_row_f = res_f._matrix[i]
                for j from 0 <= j < self._ncols:
                    res_row_f[j] = <float>fmpz_fdiv_ui(fmpz_mat_entry(self._matrix,i,j), p)
            return res_f

        elif p < MAX_MODULUS_DOUBLE:
            res_d = Matrix_modn_dense_double.__new__(Matrix_modn_dense_double,
                                                     matrix_space.MatrixSpace(IntegerModRing(p), self._nrows, self._ncols, sparse=False), None, None, None, zeroed_alloc=False)
            for i from 0 <= i < self._nrows:
                res_row_d = res_d._matrix[i]
                for j from 0 <= j < self._ncols:
                    res_row_d[j] = <double>fmpz_fdiv_ui(fmpz_mat_entry(self._matrix,i,j), p)
            return res_d
        else:
            raise ValueError("p to big.")

    def _reduce(self, moduli):
        from sage.matrix.matrix_modn_dense_float import MAX_MODULUS as MAX_MODULUS_FLOAT
        from sage.matrix.matrix_modn_dense_double import MAX_MODULUS as MAX_MODULUS_DOUBLE

        if isinstance(moduli, (int, Integer)):
            return self._mod_int(moduli)
        elif isinstance(moduli, list):
            moduli = MultiModularBasis(moduli)

        cdef MultiModularBasis mm
        mm = moduli

        res = []
        for p in mm:
            if p < MAX_MODULUS_FLOAT:
                res.append( Matrix_modn_dense_float.__new__(Matrix_modn_dense_float,
                                                            matrix_space.MatrixSpace(IntegerModRing(p), self._nrows, self._ncols, sparse=False),
                                                            None, None, None, zeroed_alloc=False) )
            elif p < MAX_MODULUS_DOUBLE:
                res.append( Matrix_modn_dense_double.__new__(Matrix_modn_dense_double,
                                                             matrix_space.MatrixSpace(IntegerModRing(p), self._nrows, self._ncols, sparse=False),
                                                             None, None, None, zeroed_alloc=False) )
            else:
                raise ValueError("p=%d too big."%p)

        cdef size_t i, k, n
        cdef Py_ssize_t nr, nc
        cdef mpz_t tmp
        mpz_init(tmp)
        n = len(mm)
        nr = self._nrows
        nc = self._ncols

        cdef mod_int *entry_list
        entry_list = <mod_int*>sig_malloc(sizeof(mod_int) * n)
        if entry_list == NULL:
            raise MemoryError("out of memory allocating multi-modular coefficient list")

        sig_on()
        for i from 0 <= i < nr:
            for j from 0 <= j < nc:
                self.get_unsafe_mpz(i,j,tmp)
                mm.mpz_reduce(tmp, entry_list)
                for k from 0 <= k < n:
                    if isinstance(res[k], Matrix_modn_dense_float):
                        (<Matrix_modn_dense_float>res[k])._matrix[i][j] = (<float>entry_list[k])%(<Matrix_modn_dense_float>res[k]).p
                    else:
                        (<Matrix_modn_dense_double>res[k])._matrix[i][j] = (<double>entry_list[k])%(<Matrix_modn_dense_double>res[k]).p
        sig_off()
        mpz_clear(tmp)
        sig_free(entry_list)
        return res

    cpdef _echelon_in_place(self, str algorithm):
        cdef Matrix_integer_dense E
        E = self.echelon_form()
        sig_on()
        fmpz_mat_set(self._matrix,E._matrix)
        sig_off()
        self.clear_cache()

    def _echelon_strassen(self):
        raise NotImplementedError

    def _magma_init_(self, magma):
        """
        EXAMPLES::

            sage: m = matrix(ZZ,2,3,[1,2,-3,1,-2,-45])
            sage: m._magma_init_(magma)
            'Matrix(IntegerRing(),2,3,StringToIntegerSequence("1 2 -3 1 -2 -45"))'
            sage: magma(m)                                               # optional - magma
            [  1   2  -3]
            [  1  -2 -45]
        """
        w = self._export_as_string(base=10)
        return 'Matrix(IntegerRing(),%s,%s,StringToIntegerSequence("%s"))'%(
            self.nrows(), self.ncols(), w)

    def symplectic_form(self):
        r"""
        Find a symplectic basis for ``self`` if ``self`` is an anti-symmetric,
        alternating matrix.

        Return a pair (F, C) such that the rows of C form a symplectic
        basis for ``self`` and ``F = C * self * C.transpose()``.

        Raise a :exc:`ValueError` if ``self`` is not anti-symmetric,
        or ``self`` is not alternating.

        Anti-symmetric means that `M = -M^t`. Alternating means
        that the diagonal of `M` is identically zero.

        A symplectic basis is a basis of the form
        `e_1, \ldots, e_j, f_1, \ldots f_j, z_1, \dots, z_k`
        such that

        -  `z_i M v^t = 0` for all vectors `v`

        -  `e_i M {e_j}^t = 0` for all `i, j`

        -  `f_i M {f_j}^t = 0` for all `i, j`

        -  `e_i M {f_i}^t = 1` for all `i`

        -  `e_i M {f_j}^t = 0` for all `i` not equal
            `j`.

        The ordering for the factors `d_{i} | d_{i+1}` and for
        the placement of zeroes was chosen to agree with the output of
        :meth:`smith_form`.

        See the example for a pictorial description of such a basis.

        EXAMPLES::

            sage: E = matrix(ZZ, 5, 5, [0, 14, 0, -8, -2, -14, 0, -3, -11, 4, 0, 3, 0, 0, 0, 8, 11, 0, 0, 8, 2, -4, 0, -8, 0]); E
            [  0  14   0  -8  -2]
            [-14   0  -3 -11   4]
            [  0   3   0   0   0]
            [  8  11   0   0   8]
            [  2  -4   0  -8   0]
            sage: F, C = E.symplectic_form()
            sage: F
            [ 0  0  1  0  0]
            [ 0  0  0  2  0]
            [-1  0  0  0  0]
            [ 0 -2  0  0  0]
            [ 0  0  0  0  0]
            sage: F == C * E * C.transpose()
            True
            sage: E.smith_form()[0]
            [1 0 0 0 0]
            [0 1 0 0 0]
            [0 0 2 0 0]
            [0 0 0 2 0]
            [0 0 0 0 0]
        """
        import sage.matrix.symplectic_basis
        return sage.matrix.symplectic_basis.symplectic_basis_over_ZZ(self)

    hermite_form = echelon_form

    def echelon_form(self, algorithm='default', proof=None, include_zero_rows=True,
                     transformation=False, D=None):
        r"""
        Return the echelon form of this matrix over the integers, also known
        as the hermite normal form (HNF).

        INPUT:

        - ``algorithm`` -- string; the algorithm to use. Valid options are:

          - ``'default'`` -- let Sage pick an algorithm (default).
            Up to 75 rows or columns with no transformation matrix,
            use pari with flag 0; otherwise, use flint.

          - ``'flint'`` -- use flint

          - ``'ntl'`` -- use NTL (only works for square matrices of
            full rank!)

          - ``'padic'`` -- an asymptotically fast `p`-adic modular
            algorithm, If your matrix has large coefficients and is
            small, you may also want to try this.

          - ``'pari'`` -- use PARI with flag 1

          - ``'pari0'`` -- use PARI with flag 0

          - ``'pari1'`` -- use PARI with flag 1

          - ``'pari4'`` -- use PARI with flag 4 (use heuristic LLL)

        - ``proof`` -- (default: ``True``) if proof=False certain
          determinants are computed using a randomized hybrid `p`-adic
          multimodular strategy until it stabilizes twice (instead of up to
          the Hadamard bound). It is *incredibly* unlikely that one would
          ever get an incorrect result with proof=False.

        - ``include_zero_rows`` -- boolean (default: ``True``); if ``False``,
          don't include zero rows

        - ``transformation`` -- if given, also compute
          transformation matrix; only valid for flint and padic algorithm

        - ``D`` -- (default: ``None``) if given and the algorithm
          is ``'ntl'``, then D must be a multiple of the determinant and this
          function will use that fact

        OUTPUT:

        The Hermite normal form (=echelon form over `\ZZ`) of ``self`` as
        an immutable matrix.

        EXAMPLES::

            sage: A = MatrixSpace(ZZ,2)([1,2,3,4])
            sage: A.echelon_form()
            [1 0]
            [0 2]
            sage: A = MatrixSpace(ZZ,5)(range(25))
            sage: A.echelon_form()
            [  5   0  -5 -10 -15]
            [  0   1   2   3   4]
            [  0   0   0   0   0]
            [  0   0   0   0   0]
            [  0   0   0   0   0]

        Getting a transformation matrix in the nonsquare case::

            sage: A = matrix(ZZ,5,3,[1..15])
            sage: H, U = A.hermite_form(transformation=True, include_zero_rows=False)
            sage: H
            [1 2 3]
            [0 3 6]
            sage: U
            [  0   0   0   4  -3]
            [  0   0   0  13 -10]
            sage: U*A == H
            True

        TESTS:

        Make sure the zero matrices are handled correctly::

            sage: m = matrix(ZZ,3,3,[0]*9)
            sage: m.echelon_form()
            [0 0 0]
            [0 0 0]
            [0 0 0]
            sage: m = matrix(ZZ,3,1,[0]*3)
            sage: m.echelon_form()
            [0]
            [0]
            [0]
            sage: m = matrix(ZZ,1,3,[0]*3)
            sage: m.echelon_form()
            [0 0 0]

        The ultimate border case!

        ::

            sage: m = matrix(ZZ,0,0,[])
            sage: m.echelon_form()
            []

        .. NOTE::

           If 'ntl' is chosen for a non square matrix this function
           raises a :exc:`ValueError`.

        Special cases: 0 or 1 rows::

            sage: a = matrix(ZZ, 1,2,[0,-1])
            sage: a.hermite_form()
            [0 1]
            sage: a.pivots()
            (1,)
            sage: a = matrix(ZZ, 1,2,[0,0])
            sage: a.hermite_form()
            [0 0]
            sage: a.pivots()
            ()
            sage: a = matrix(ZZ,1,3); a
            [0 0 0]
            sage: a.echelon_form(include_zero_rows=False)
            []
            sage: a.echelon_form(include_zero_rows=True)
            [0 0 0]

        Illustrate using various algorithms.::

            sage: matrix(ZZ,3,[1..9]).hermite_form(algorithm='pari')
            [1 2 3]
            [0 3 6]
            [0 0 0]
            sage: matrix(ZZ,3,[1..9]).hermite_form(algorithm='pari0')
            [1 2 3]
            [0 3 6]
            [0 0 0]
            sage: matrix(ZZ,3,[1..9]).hermite_form(algorithm='pari4')
            [1 2 3]
            [0 3 6]
            [0 0 0]
            sage: matrix(ZZ,3,[1..9]).hermite_form(algorithm='padic')
            [1 2 3]
            [0 3 6]
            [0 0 0]
            sage: matrix(ZZ,3,[1..9]).hermite_form(algorithm='default')
            [1 2 3]
            [0 3 6]
            [0 0 0]

        The 'ntl' algorithm doesn't work on matrices that do not have full rank.::

            sage: matrix(ZZ,3,[1..9]).hermite_form(algorithm='ntl')
            Traceback (most recent call last):
            ...
            ValueError: ntl only computes HNF for square matrices of full rank.
            sage: matrix(ZZ,3,[0] +[2..9]).hermite_form(algorithm='ntl')
            [1 0 0]
            [0 1 0]
            [0 0 3]

        TESTS:

        This example illustrated :issue:`2398`::

            sage: a = matrix([(0, 0, 3), (0, -2, 2), (0, 1, 2), (0, -2, 5)])
            sage: a.hermite_form()
            [0 1 2]
            [0 0 3]
            [0 0 0]
            [0 0 0]

        Check that :issue:`12280` is fixed::

            sage: m = matrix([(-2, 1, 9, 2, -8, 1, -3, -1, -4, -1),
            ....:             (5, -2, 0, 1, 0, 4, -1, 1, -2, 0),
            ....:             (-11, 3, 1, 0, -3, -2, -1, -11, 2, -2),
            ....:             (-1, 1, -1, -2, 1, -1, -1, -1, -1, 7),
            ....:             (-2, -1, -1, 1, 1, -2, 1, 0, 2, -4)]).stack(
            ....:             200 * identity_matrix(ZZ, 10))
            sage: matrix(ZZ,m).hermite_form(algorithm='pari', include_zero_rows=False)
            [  1   0   2   0  13   5   1 166  72  69]
            [  0   1   1   0  20   4  15 195  65 190]
            [  0   0   4   0  24   5  23  22  51 123]
            [  0   0   0   1  23   7  20 105  60 151]
            [  0   0   0   0  40   4   0  80  36  68]
            [  0   0   0   0   0  10   0 100 190 170]
            [  0   0   0   0   0   0  25   0 100 150]
            [  0   0   0   0   0   0   0 200   0   0]
            [  0   0   0   0   0   0   0   0 200   0]
            [  0   0   0   0   0   0   0   0   0 200]
            sage: matrix(ZZ,m).hermite_form(algorithm='padic', include_zero_rows=False)
            [  1   0   2   0  13   5   1 166  72  69]
            [  0   1   1   0  20   4  15 195  65 190]
            [  0   0   4   0  24   5  23  22  51 123]
            [  0   0   0   1  23   7  20 105  60 151]
            [  0   0   0   0  40   4   0  80  36  68]
            [  0   0   0   0   0  10   0 100 190 170]
            [  0   0   0   0   0   0  25   0 100 150]
            [  0   0   0   0   0   0   0 200   0   0]
            [  0   0   0   0   0   0   0   0 200   0]
            [  0   0   0   0   0   0   0   0   0 200]

        Check that the output is correct in corner cases, see :issue:`18613`::

            sage: m = matrix(2, 0)
            sage: m.parent()
            Full MatrixSpace of 2 by 0 dense matrices over Integer Ring
            sage: H, U = m.echelon_form(transformation=True)
            sage: H.parent()
            Full MatrixSpace of 2 by 0 dense matrices over Integer Ring
            sage: H.is_immutable()
            True
            sage: U
            [1 0]
            [0 1]
            sage: H == U * m
            True
            sage: H, U = m.echelon_form(transformation=True,
            ....:                       include_zero_rows=False)
            sage: H.parent()
            Full MatrixSpace of 0 by 0 dense matrices over Integer Ring
            sage: U.parent()
            Full MatrixSpace of 0 by 2 dense matrices over Integer Ring
            sage: H == U * m
            True
            sage: m = random_matrix(ZZ, 15, 15, x=-1000, y=1000, density=0.1)
            sage: m.parent()
            Full MatrixSpace of 15 by 15 dense matrices over Integer Ring
            sage: H, U = m.hermite_form(algorithm='flint',
            ....:                       transformation=True)
            sage: H == U*m
            True
        """
        key = 'hnf-%s-%s'%(include_zero_rows,transformation)
        ans = self.fetch(key)
        if ans is not None:
            return ans

        cdef Matrix_integer_dense H_m, U
        cdef Py_ssize_t nr, nc, n, i, j
        nr = self._nrows
        nc = self._ncols
        n = nr if nr >= nc else nc
        if algorithm == 'default':
            if transformation:
                algorithm = 'flint'
            else:
                if n < 75:
                    algorithm = 'pari0'
                else:
                    algorithm = 'flint'
        proof = get_proof_flag(proof, "linear_algebra")
        pivots = None

        if nr == 0 or nc == 0:
            pivots = ()
            if include_zero_rows:
                H_m = self.new_matrix()
                U = self.matrix_space(nr, nr).one()
            else:
                H_m = self.new_matrix(0, nc)
                U = self.new_matrix(0, nr)
        elif algorithm == "flint":
            H_m = self._new(nr, nc)

            if transformation:
                U = self._new(nr, nr)
                sig_on()
                fmpz_mat_hnf_transform(H_m._matrix, U._matrix, self._matrix)
                sig_off()
            else:
                sig_on()
                fmpz_mat_hnf(H_m._matrix, self._matrix)
                sig_off()
            if not include_zero_rows:
                r = H_m.rank()
                H_m = H_m[:r]
                if transformation:
                    U = U[:r]
        elif algorithm == "padic":
            from sage.matrix import matrix_integer_dense_hnf
            if transformation:
                H_m, U = matrix_integer_dense_hnf.hnf_with_transformation(self, proof=proof)
                if not include_zero_rows:
                    r = H_m.rank()
                    H_m = H_m[:r]
                    U = U[:r]
            else:
                H_m, pivots = matrix_integer_dense_hnf.hnf(self,
                                   include_zero_rows=include_zero_rows, proof=proof)
        elif transformation:
            raise ValueError("transformation matrix only available with p-adic algorithm")
        elif algorithm in ["pari", "pari0", "pari1", "pari4"]:
            flag = int(algorithm[-1]) if algorithm != "pari" else 1
            H_m = self._hnf_pari(flag, include_zero_rows=include_zero_rows)
        elif algorithm == 'ntl':
            if nr != nc:
                raise ValueError("ntl only computes HNF for square matrices of full rank.")

            import sage.libs.ntl.ntl_mat_ZZ
            v =  sage.libs.ntl.ntl_mat_ZZ.ntl_mat_ZZ(self._nrows,self._ncols)
            for i from 0 <= i < self._nrows:
                for j from 0 <= j < self._ncols:
                    v[i,j] = self.get_unsafe(nr-i-1,nc-j-1)

            try:
                w1 = v.HNF(D=D)
            except RuntimeError: # HNF may fail if a nxm matrix has rank < m
                raise ValueError("ntl only computes HNF for square matrices of full rank.")

            if include_zero_rows:
                H_m = self.new_matrix()
            else:
                H_m = self.new_matrix(nrows=w1.nrows())

            nr = w1.nrows()
            nc = w1.ncols()

            for i from 0 <= i < w1.nrows():
                for j from 0 <= j < w1.ncols():
                    H_m[i,j] = w1[nr-i-1,nc-j-1]

        else:
            raise ValueError("algorithm %r not understood" % algorithm)

        H_m.set_immutable()
        if pivots is None:
            from sage.matrix.matrix_integer_dense_hnf import pivots_of_hnf_matrix
            pivots = pivots_of_hnf_matrix(H_m)
        pivots = tuple(pivots)
        rank = len(pivots)
        H_m.cache('pivots', pivots)
        self.cache('pivots', pivots)
        H_m.cache('rank', rank)
        self.cache('rank',rank)
        H_m.cache('in_echelon_form', True)

        if transformation:
            U.set_immutable()
            ans = H_m, U
        else:
            ans = H_m
        self.cache(key, ans)
        return ans

    def saturation(self, p=0, proof=None, max_dets=5):
        r"""
        Return a saturation matrix of self, which is a matrix whose rows
        span the saturation of the row span of ``self``. This is not unique.

        The saturation of a `\ZZ` module `M`
        embedded in `\ZZ^n` is a module `S` that
        contains `M` with finite index such that
        `\ZZ^n/S` is torsion free. This function takes the
        row span `M` of self, and finds another matrix of full rank
        with row span the saturation of `M`.

        INPUT:

        - ``p`` -- (default: 0) if nonzero given, saturate
          only at the prime `p`, i.e., return a matrix whose row span
          is a `\ZZ`-module `S` that contains ``self`` and
          such that the index of `S` in its saturation is coprime to
          `p`. If `p` is None, return full saturation of ``self``.

        - ``proof`` -- (default: use proof.linear_algebra());
          if ``False``, the determinant calculations are done with
          ``proof=False``

        - ``max_dets`` -- (default: 5) technical parameter -
          max number of determinant to compute when bounding prime divisor of
          ``self`` in its saturation.

        OUTPUT: matrix over ZZ


        .. NOTE::

           The result is *not* cached.

        ALGORITHM: 1. Replace input by a matrix of full rank got from a
        subset of the rows. 2. Divide out any common factors from rows. 3.
        Check max_dets random dets of submatrices to see if their GCD
        (with p) is 1 - if so matrix is saturated and we're done. 4.
        Finally, use that if A is a matrix of full rank, then
        `hnf(transpose(A))^{-1}*A` is a saturation of A.

        EXAMPLES::

            sage: A = matrix(ZZ, 3, 5, [-51, -1509, -71, -109, -593, -19, -341, 4, 86, 98, 0, -246, -11, 65, 217])
            sage: A.echelon_form()
            [      1       5    2262   20364   56576]
            [      0       6   35653  320873  891313]
            [      0       0   42993  386937 1074825]
            sage: S = A.saturation(); S
            [  -51 -1509   -71  -109  -593]
            [  -19  -341     4    86    98]
            [   35   994    43    51   347]

        Notice that the saturation spans a different module than A.

        ::

            sage: S.echelon_form()
            [ 1  2  0  8 32]
            [ 0  3  0 -2 -6]
            [ 0  0  1  9 25]
            sage: V = A.row_space(); W = S.row_space()
            sage: V.is_submodule(W)
            True
            sage: V.index_in(W)
            85986
            sage: V.index_in_saturation()
            85986

        We illustrate each option::

            sage: S = A.saturation(p=2)
            sage: S = A.saturation(proof=False)
            sage: S = A.saturation(max_dets=2)
        """
        proof = get_proof_flag(proof, "linear_algebra")
        from sage.matrix.matrix_integer_dense_saturation import saturation
        return saturation(self, p=p, proof=proof, max_dets=max_dets)

    def index_in_saturation(self, proof=None):
        """
        Return the index of ``self`` in its saturation.

        INPUT:

        - ``proof`` -- (default: use proof.linear_algebra());
          if ``False``, the determinant calculations are done with
          ``proof=False``

        OUTPUT: positive integer; the index of the row span of
        this matrix in its saturation


        ALGORITHM:

        Use Hermite normal form twice to find an invertible matrix whose
        inverse transforms a matrix with the same row span as ``self`` to its
        saturation, then compute the determinant of that matrix.

        EXAMPLES::

            sage: A = matrix(ZZ, 2,3, [1..6]); A
            [1 2 3]
            [4 5 6]
            sage: A.index_in_saturation()
            3
            sage: A.saturation()
            [1 2 3]
            [1 1 1]
        """
        proof = get_proof_flag(proof, "linear_algebra")
        from sage.matrix.matrix_integer_dense_saturation import index_in_saturation
        return index_in_saturation(self, proof=proof)

    def pivots(self):
        """
        Return the pivot column positions of this matrix.

        OUTPUT: a tuple of Python integers: the position of the
        first nonzero entry in each row of the echelon form.

        EXAMPLES::

            sage: n = 3; A = matrix(ZZ,n,range(n^2)); A
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: A.pivots()
            (0, 1)
            sage: A.echelon_form()
            [ 3  0 -3]
            [ 0  1  2]
            [ 0  0  0]
        """
        p = self.fetch('pivots')
        if p is not None:
            return tuple(p)

        cdef Matrix_integer_dense E
        E = self.echelon_form()

        # Now we determine the pivots from the matrix E as quickly as we can.
        # For each row, we find the first nonzero position in that row -- it is the pivot.
        cdef Py_ssize_t i, j, k
        p = []
        k = 0
        for i from 0 <= i < E._nrows:
            for j from k <= j < E._ncols:
                if fmpz_cmp_si(fmpz_mat_entry(E._matrix,i,j), 0) != 0:  # nonzero position
                    p.append(j)
                    k = j+1  # so start at next position next time
                    break
        p = tuple(p)
        self.cache('pivots', p)
        return p

    #### Elementary divisors

    def elementary_divisors(self, algorithm='pari'):
        """
        Return the elementary divisors of self, in order.

        .. WARNING::

           This is MUCH faster than the :meth:`smith_form` function.

        The elementary divisors are the invariants of the finite abelian
        group that is the cokernel of *left* multiplication of this matrix.
        They are ordered in reverse by divisibility.

        INPUT:

        - ``self`` -- matrix

        - ``algorithm`` -- (default: ``'pari'``)

           - ``'pari'`` -- works robustly, but is slower.

           - ``'linbox'`` -- use linbox (currently off, broken)

        OUTPUT: list of integers


        .. NOTE::

           These are the invariants of the cokernel of *left* multiplication::

               sage: M = Matrix([[3,0,1],[0,1,0]])
               sage: M
               [3 0 1]
               [0 1 0]
               sage: M.elementary_divisors()
               [1, 1]
               sage: M.transpose().elementary_divisors()
               [1, 1, 0]

        EXAMPLES::

            sage: matrix(3, range(9)).elementary_divisors()
            [1, 3, 0]
            sage: matrix(3, range(9)).elementary_divisors(algorithm='pari')
            [1, 3, 0]
            sage: C = MatrixSpace(ZZ,4)([3,4,5,6,7,3,8,10,14,5,6,7,2,2,10,9])
            sage: C.elementary_divisors()
            [1, 1, 1, 687]

        ::

            sage: M = matrix(ZZ, 3, [1,5,7, 3,6,9, 0,1,2])
            sage: M.elementary_divisors()
            [1, 1, 6]

        This returns a copy, which is safe to change::

            sage: edivs = M.elementary_divisors()
            sage: edivs.pop()
            6
            sage: M.elementary_divisors()
            [1, 1, 6]

        .. SEEALSO::

           :meth:`smith_form`
        """
        d = self.fetch('elementary_divisors')
        if d is not None:
            return d[:]
        if self._nrows == 0 or self._ncols == 0:
            d = []
        elif algorithm == 'linbox':
            raise ValueError("linbox too broken -- currently Linbox SNF is disabled.")
        elif algorithm == 'pari':
            d = self.__pari__().matsnf(0).sage()
            i = d.count(0)
            d.sort()
            if i > 0:
                d = d[i:] + [d[0]]*i
        else:
            raise ValueError("algorithm (='%s') unknown"%algorithm)

        self.cache('elementary_divisors', d)
        return d[:]

    def smith_form(self, transformation=True, integral=None):
        r"""
        Return the smith normal form of this matrix, that is the diagonal
        matrix `S` with diagonal entries the ordered elementary divisors of
        this matrix.

        INPUT:

        - ``transformation`` -- boolean (default: ``True``); whether to
          return the transformation matrices `U` and `V` such that `S=U\cdot
          self\cdot V`

        - ``integral`` -- a subring of the base ring or ``True`` (default:
          ``None``); ignored for matrices with integer entries

        .. NOTE::

           The :meth:`elementary_divisors` function, which returns the diagonal
           entries of `S`, is VASTLY faster than this function.

        The elementary divisors are the invariants of the finite abelian
        group that is the cokernel of this matrix. They are ordered in
        reverse by divisibility.

        EXAMPLES::

            sage: A = MatrixSpace(IntegerRing(), 3)(range(9))
            sage: D, U, V = A.smith_form()
            sage: D
            [1 0 0]
            [0 3 0]
            [0 0 0]
            sage: U
            [ 0  2 -1]
            [ 0 -1  1]
            [ 1 -2  1]
            sage: V
            [ 0  0  1]
            [-1  2 -2]
            [ 1 -1  1]
            sage: U*A*V
            [1 0 0]
            [0 3 0]
            [0 0 0]

        It also makes sense for nonsquare matrices::

            sage: A = Matrix(ZZ,3,2,range(6))
            sage: D, U, V = A.smith_form()
            sage: D
            [1 0]
            [0 2]
            [0 0]
            sage: U
            [ 0  2 -1]
            [ 0 -1  1]
            [ 1 -2  1]
            sage: V
            [-1  1]
            [ 1  0]
            sage: U * A * V
            [1 0]
            [0 2]
            [0 0]

        Empty matrices are handled sensibly (see :issue:`3068`)::

            sage: m = MatrixSpace(ZZ, 2,0)(0); d,u,v = m.smith_form(); u*m*v == d
            True
            sage: m = MatrixSpace(ZZ, 0,2)(0); d,u,v = m.smith_form(); u*m*v == d
            True
            sage: m = MatrixSpace(ZZ, 0,0)(0); d,u,v = m.smith_form(); u*m*v == d
            True

        .. SEEALSO::

           :meth:`elementary_divisors`
        """
        X = self.matrix_space()([self[i,j] for i in range(self._nrows-1,-1,-1) for j in range(self._ncols-1,-1,-1)])
        v = X.__pari__().matsnf(1).sage()
        # need to reverse order of rows of U, columns of V, and both of D.
        D = self.matrix_space()([v[2][i,j] for i in range(self._nrows-1,-1,-1) for j in range(self._ncols-1,-1,-1)])

        if not transformation:
            return D

        if self._ncols == 0:
            v[0] = self.matrix_space(ncols = self._nrows)(1)
        if self._nrows == 0:
            v[1] = self.matrix_space(nrows = self._ncols)(1)

        if self._ncols == 0:
            # silly special cases for matrices with 0 columns (PARI has a unique empty matrix)
            U = self.matrix_space(ncols = self._nrows)(1)
        else:
            U = self.matrix_space(ncols = self._nrows)([v[0][i,j] for i in range(self._nrows-1,-1,-1) for j in range(self._nrows-1,-1,-1)])

        if self._nrows == 0:
            # silly special cases for matrices with 0 rows (PARI has a unique empty matrix)
            V = self.matrix_space(nrows = self._ncols)(1)
        else:
            V = self.matrix_space(nrows = self._ncols)([v[1][i,j] for i in range(self._ncols-1,-1,-1) for j in range(self._ncols-1,-1,-1)])

        return D, U, V

    def frobenius_form(self, flag=0, var='x'):
        """
        Return the Frobenius form (rational canonical form) of this
        matrix.

        INPUT:

        - ``flag`` -- 0 (default), 1 or 2 as follows:

            - ``0`` -- (default) return the Frobenius form of this
               matrix

            - ``1`` -- return only the elementary divisor
               polynomials, as polynomials in var

            - ``2`` -- return a two-components vector [F,B] where F
               is the Frobenius form and B is the basis change so that
               `M=B^{-1}FB`

        - ``var`` -- string (default: ``'x'``)

        ALGORITHM: uses PARI's :pari:`matfrobenius`

        EXAMPLES::

            sage: A = MatrixSpace(ZZ, 3)(range(9))
            sage: A.frobenius_form(0)
            [ 0  0  0]
            [ 1  0 18]
            [ 0  1 12]
            sage: A.frobenius_form(1)
            [x^3 - 12*x^2 - 18*x]
            sage: A.frobenius_form(1, var='y')
            [y^3 - 12*y^2 - 18*y]
            sage: F, B = A.frobenius_form(2)
            sage: A == B^(-1)*F*B
            True
            sage: a=matrix([])
            sage: a.frobenius_form(2)
            ([], [])
            sage: a.frobenius_form(0)
            []
            sage: a.frobenius_form(1)
            []
            sage: B = random_matrix(ZZ,2,3)
            sage: B.frobenius_form()
            Traceback (most recent call last):
            ...
            ArithmeticError: frobenius matrix of non-square matrix not defined.

        AUTHORS:

        - Martin Albrecht (2006-04-02)

        TODO: - move this to work for more general matrices than just over
        Z. This will require fixing how PARI polynomials are coerced to
        Sage polynomials.
        """
        if not self.is_square():
            raise ArithmeticError("frobenius matrix of non-square matrix not defined.")

        v = self.__pari__().matfrobenius(flag)
        if flag == 0:
            return self.matrix_space()(v.sage())
        elif flag == 1:
            r = PolynomialRing(self.base_ring(), names=var)
            retr = []
            for f in v:
                retr.append(eval(str(f).replace("^","**"),
                                 {'x': r.gen()}, r.gens_dict()))
            return retr
        elif flag == 2:
            F = matrix_space.MatrixSpace(QQ, self.nrows())(v[0].sage())
            B = matrix_space.MatrixSpace(QQ, self.nrows())(v[1].sage())
            return F, B

    frobenius = deprecated_function_alias(36396, frobenius_form)

    def _right_kernel_matrix(self, **kwds):
        r"""
        Return a pair that includes a matrix of basis vectors
        for the right kernel of ``self``.

        INPUT:

        - ``algorithm`` -- determines which algorithm to use, options are:

          - ``'flint'`` -- use the algorithm from the FLINT library
          - ``'pari'`` -- use the :pari:`matkerint` function from the PARI library
          - ``'padic'`` -- use the `p`-adic algorithm from the IML library
          - ``'default'`` -- use a heuristic to decide which of the three above
            routines is fastest.  This is the default value

        - ``proof`` -- this is passed to the `p`-adic IML algorithm;
          if not specified, the global flag for linear algebra will be used

        OUTPUT:

        Return a pair.  First item is the string is either
        'computed-flint-int', 'computed-pari-int', 'computed-flint-int', which identifies
        the nature of the basis vectors.

        Second item is a matrix whose rows are a basis for the right kernel,
        over the integers, as computed by either the FLINT, IML or PARI libraries.

        EXAMPLES::

            sage: A = matrix(ZZ, [[4, 7, 9, 7, 5, 0],
            ....:                 [1, 0, 5, 8, 9, 1],
            ....:                 [0, 1, 0, 1, 9, 7],
            ....:                 [4, 7, 6, 5, 1, 4]])

            sage: result = A._right_kernel_matrix(algorithm='pari')
            sage: result[0]
            'computed-pari-int'
            sage: X = result[1]; X
            [ 26 -31  30 -21  -2  10]
            [-47 -13  48 -14 -11  18]
            sage: A*X.transpose() == zero_matrix(ZZ, 4, 2)
            True

            sage: result = A._right_kernel_matrix(algorithm='padic')
            sage: result[0]
            'computed-iml-int'
            sage: X = result[1]; X
            [-469  214  -30  119  -37    0]
            [ 370 -165   18  -91   30   -2]
            sage: A*X.transpose() == zero_matrix(ZZ, 4, 2)
            True

            sage: result = A._right_kernel_matrix(algorithm='default')
            sage: result[0]
            'computed-flint-int'
            sage: result[1]
            [ 469 -214   30 -119   37    0]
            [-370  165  -18   91  -30    2]

            sage: result = A._right_kernel_matrix(algorithm='flint')
            sage: result[0]
            'computed-flint-int'
            sage: result[1]
            [ 469 -214   30 -119   37    0]
            [-370  165  -18   91  -30    2]

        With the 'default' given as the algorithm, several heuristics are
        used to determine if FLINT, PARI or IML ('padic') is used.  The code has
        exact details, but roughly speaking, relevant factors are: the
        absolute size of the matrix, or the relative dimensions, or the
        magnitude of the entries. ::

            sage: A = random_matrix(ZZ, 18, 11)
            sage: A._right_kernel_matrix(algorithm='default')[0]
            'computed-pari-int'
            sage: A = random_matrix(ZZ, 18, 11, x = 10^200)
            sage: A._right_kernel_matrix(algorithm='default')[0]
            'computed-iml-int'
            sage: A = random_matrix(ZZ, 60, 60)
            sage: A._right_kernel_matrix(algorithm='default')[0]
            'computed-iml-int'
            sage: A = random_matrix(ZZ, 60, 55)
            sage: A._right_kernel_matrix(algorithm='default')[0]
            'computed-pari-int'

        TESTS:

        We test three trivial cases. FLINT is used for small matrices,
        but we let the heuristic decide that.  ::

            sage: A = matrix(ZZ, 0, 2)
            sage: A._right_kernel_matrix()[1]
            []
            sage: A = matrix(ZZ, 2, 0)
            sage: A._right_kernel_matrix()[1].parent()
            Full MatrixSpace of 0 by 0 dense matrices over Integer Ring
            sage: A = zero_matrix(ZZ, 4, 3)
            sage: A._right_kernel_matrix()[1]
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        tm = verbose("computing right kernel matrix over the integers for %sx%s matrix" % (self.nrows(), self.ncols()),level=1)

        algorithm = kwds.pop('algorithm', None)
        if algorithm is None:
            algorithm = 'default'

        if algorithm == 'default':
            # The heuristic here could be auto-tuned, stored for
            # different architecture, etc.  What I've done below here
            # I just got by playing around with examples.  This is
            # *dramatically* better than doing absolutely nothing
            # (i.e., always choosing 'padic'), but is of course
            # far from optimal.   -- William Stein

            # I sometimes favor FLINT over PARI, but this should be better tuned. -- Marc Masdeu
            if max(self._nrows, self._ncols) <= 10:
                # Use FLINT for very small matrices, as long as entries aren't huge.
                algorithm = 'flint'
            elif max(self._nrows, self._ncols) <= 50:
                # when entries are huge, padic relatively good.
                h = self.height().ndigits()
                if h < 100:
                    algorithm = 'pari'
                else:
                    algorithm = 'padic'
            elif self._nrows <= self._ncols + 3:
                # the padic algorithm is much better for bigger
                # matrices if there are nearly more columns than rows
                # (that is its forte)
                algorithm = 'padic'
            else:
                algorithm = 'pari'

        if algorithm == 'flint':
            proof = kwds.pop('proof', None)
            proof = get_proof_flag(proof, "linear_algebra")
            K = self._rational_kernel_flint().transpose().saturation(proof=proof)
            format = 'computed-flint-int'
        elif algorithm == 'pari':
            K = self.__pari__().matkerint().mattranspose().sage()
            format = 'computed-pari-int'
        elif algorithm == 'padic':
            proof = kwds.pop('proof', None)
            proof = get_proof_flag(proof, "linear_algebra")
            K = self._rational_kernel_iml().transpose().saturation(proof=proof)
            format = 'computed-iml-int'
        else:
            raise ValueError('unknown algorithm: %s'%algorithm)
        tm = verbose("done computing right kernel matrix over the integers for %sx%s matrix" % (self.nrows(), self.ncols()),level=1, t=tm)
        return (format, K)

    # TODO: implement using flint function
    def _adjugate(self):
        """
        Return the adjugate (classical adjoint) of this matrix.

        Assumes ``self`` is a square matrix (checked in adjugate).

        EXAMPLES::

            sage: m = matrix(ZZ,3,[1..9])
            sage: m.adjugate()  # indirect doctest
            [ -3   6  -3]
            [  6 -12   6]
            [ -3   6  -3]
        """
        return self.parent()(self.__pari__().matadjoint().sage())

    def _ntl_(self):
        r"""
        ntl.mat_ZZ representation of ``self``.

        EXAMPLES::

            sage: a = MatrixSpace(ZZ,200).random_element(x=-2, y=2)    # -2 to 2
            sage: A = a._ntl_()

        .. NOTE::

           NTL only knows dense matrices, so if you provide a sparse
           matrix NTL will allocate memory for every zero entry.
        """
        import sage.libs.ntl.ntl_mat_ZZ
        return sage.libs.ntl.ntl_mat_ZZ.ntl_mat_ZZ(self._nrows,self._ncols, self.list())

    #######################################################################
    # LLL
    #######################################################################

    def BKZ(self, delta=None, algorithm='fpLLL', fp=None, block_size=10, prune=0,
            use_givens=False, precision=0, proof=None, **kwds):
        """
        Return the result of running Block Korkin-Zolotarev reduction on
        ``self`` interpreted as a lattice.

        INPUT:

        - ``delta`` -- (default: ``0.99``) LLL parameter

        - ``algorithm`` -- (default: ``'fpLLL'``) ``'fpLLL'`` or ``"NTL"``

        - ``fp`` -- floating point number implementation

            - ``None`` -- NTL's exact reduction or fpLLL's wrapper (default)

            - ``'fp'`` -- double precision: NTL's FP or fpLLL's double

            - ``'ld'`` -- long doubles (fpLLL only)

            - ``'qd'`` -- NTL's QP

            - ``'qd1'`` -- quad doubles: Uses ``quad_float`` precision
              to compute Gram-Schmidt, but uses double precision in
              the search phase of the block reduction algorithm. This
              seems adequate for most purposes, and is faster than
              ``'qd'``, which uses quad_float precision uniformly
              throughout (NTL only).

            - ``'xd'`` -- extended exponent: NTL's XD or fpLLL's dpe

            - ``'rr'`` -- arbitrary precision: NTL'RR or fpLLL's MPFR

        - ``block_size`` -- (default: ``10``) specifies the size
          of the blocks in the reduction.  High values yield
          shorter vectors, but the running time increases double
          exponentially with ``block_size``.  ``block_size``
          should be between 2 and the number of rows of ``self``.

        - ``proof`` -- (default: same as ``proof.linear_algebra()``)
          Insist on full BKZ reduction. If disabled and fplll is
          called, reduction is much faster but the result is not fully
          BKZ reduced.

        NTL SPECIFIC INPUT:

        - ``prune`` -- (default: ``0``) the optional parameter
          ``prune`` can be set to any positive number to invoke the
          Volume Heuristic from [SH1995]_. This can significantly reduce
          the running time, and hence allow much bigger block size,
          but the quality of the reduction is of course not as good in
          general. Higher values of ``prune`` mean better quality, and
          slower running time. When ``prune`` is ``0``, pruning is
          disabled. Recommended usage: for ``block_size==30``, set
          ``10 <= prune <=15``.

        - ``use_givens`` -- use Givens orthogonalization.  Only
          applies to approximate reduction using NTL.  This is a bit
          slower, but generally much more stable, and is really the
          preferred orthogonalization strategy.  For a nice
          description of this, see Chapter 5 of [GL1996]_.

        fpLLL SPECIFIC INPUT:

        - ``precision`` -- (default: ``0`` for automatic choice) bit
          precision to use if ``fp='rr'`` is set

        - ``**kwds`` -- keywords to be passed to :mod:`fpylll`; see
          :class:`fpylll.BKZ.Param` for details

        Also, if the verbose level is at least `2`, some output
        is printed during the computation.

        EXAMPLES::

            sage: A = Matrix(ZZ,3,3,range(1,10))
            sage: A.BKZ()
            [ 0  0  0]
            [ 2  1  0]
            [-1  1  3]

            sage: A = Matrix(ZZ,3,3,range(1,10))
            sage: A.BKZ(use_givens=True)
            [ 0  0  0]
            [ 2  1  0]
            [-1  1  3]

            sage: A = Matrix(ZZ,3,3,range(1,10))
            sage: A.BKZ(fp='fp')
            [ 0  0  0]
            [ 2  1  0]
            [-1  1  3]

        ALGORITHM:

        Calls either NTL or fpLLL.
        """
        if delta is None:
            delta = 0.99
        elif delta <= 0.25:
            raise TypeError("delta must be > 0.25")
        elif delta > 1:
            raise TypeError("delta must be <= 1")
        delta = float(delta)

        if prune < 0:
            raise TypeError("prune must be >= 0")
        prune = int(prune)

        verbose = get_verbose() >= 2

        block_size = int(block_size)

        if algorithm == "NTL":
            if fp is None:
                fp = "rr"

            if fp == "fp":
                algorithm = "BKZ_FP"
            elif fp == "qd":
                algorithm = "BKZ_QP"
            elif fp == "qd1":
                algorithm = "BKZ_QP1"
            elif fp == "xd":
                algorithm = "BKZ_XD"
            elif fp == "rr":
                algorithm = "BKZ_RR"
            else:
                raise TypeError("fp parameter not understood.")

            A = self._ntl_()

            if algorithm == "BKZ_FP":
                if not use_givens:
                    r = A.BKZ_FP(U=None, delta=delta, BlockSize=block_size,
                                 prune=prune, verbose=verbose)
                else:
                    r = A.G_BKZ_FP(U=None, delta=delta, BlockSize=block_size,
                                   prune=prune, verbose=verbose)

            elif algorithm == "BKZ_QP":
                if not use_givens:
                    r = A.BKZ_QP(U=None, delta=delta, BlockSize=block_size,
                                 prune=prune, verbose=verbose)
                else:
                    r = A.G_BKZ_QP(U=None, delta=delta, BlockSize=block_size,
                                   prune=prune, verbose=verbose)

            elif algorithm == "BKZ_QP1":
                if not use_givens:
                    r = A.BKZ_QP1(U=None, delta=delta, BlockSize=block_size,
                                  prune=prune, verbose=verbose)
                else:
                    r = A.G_BKZ_QP1(U=None, delta=delta, BlockSize=block_size,
                                    prune=prune, verbose=verbose)

            elif algorithm == "BKZ_XD":
                if not use_givens:
                    r = A.BKZ_XD(U=None, delta=delta, BlockSize=block_size,
                                 prune=prune, verbose=verbose)
                else:
                    r = A.G_BKZ_XD(U=None, delta=delta, BlockSize=block_size,
                                   prune=prune, verbose=verbose)

            elif algorithm == "BKZ_RR":
                if not use_givens:
                    r = A.BKZ_RR(U=None, delta=delta, BlockSize=block_size,
                                 prune=prune, verbose=verbose)
                else:
                    r = A.G_BKZ_RR(U=None, delta=delta, BlockSize=block_size,
                                   prune=prune, verbose=verbose)

            self.cache("rank",ZZ(r))
            R = <Matrix_integer_dense>self.new_matrix(
                    entries=[ZZ(z) for z in A.list()])

        elif algorithm == "fpLLL":
            from fpylll import BKZ, IntegerMatrix, load_strategies_json
            fp_ = fplll_fp_map[fp]
            if verbose:
                kwds["flags"] = kwds.get("flags", BKZ.DEFAULT) | BKZ.VERBOSE

            proof = get_proof_flag(proof, "linear_algebra")

            # enable performance improvements unless
            # 1. provable results are requested or
            # 2. the user has specified the relevant parameters already
            if "strategies" not in kwds:
                if proof is False:
                    import os
                    kwds["strategies"] = load_strategies_json(
                        os.path.normpath(os.path.join(
                            BKZ.DEFAULT_STRATEGY_PATH,
                            os.path.basename(BKZ.DEFAULT_STRATEGY))))

            if "auto_abort" not in kwds:
                if proof is False:
                    kwds["flags"] = kwds.get("flags", BKZ.DEFAULT) | BKZ.AUTO_ABORT
                    kwds["auto_abort"] = True

            A = IntegerMatrix.from_matrix(self)
            BKZ.reduction(A, BKZ.Param(block_size=block_size, delta=delta, **kwds),
                          float_type=fp_,
                          precision=precision)

            R = A.to_matrix(self.new_matrix())
        return R

    def LLL(self, delta=None, eta=None, algorithm='fpLLL:wrapper', fp=None, prec=0, early_red=False, use_givens=False, use_siegel=False, transformation=False, **kwds):
        r"""
        Return LLL-reduced or approximated LLL reduced matrix `R` of the lattice
        generated by the rows of ``self``.

        A set of vectors `(b_1, b_2, ..., b_d)` is `(\delta, \eta)`-LLL-reduced
        if the two following conditions hold:

        -  For any `i > j`, we have `\lvert \mu_{i,j} \rvert \leq \eta`.

        -  For any `i < d`, we have `\delta \lvert b_i^* \rvert^2 \leq
           \lvert b_{i + 1}^* + \mu_{i+1, i} b_i^* \rvert^2`,

        where `\mu_{i,j} = \langle b_i, b_j^* \rangle / \langle b_j^*, b_j^*
        \rangle` and `b_i^*` is the `i`-th vector of the Gram-Schmidt
        orthogonalisation of `(b_1, b_2, ..., b_d)`.

        The default reduction parameters are `\delta = 0.99` and `\eta = 0.501`.
        The parameters `\delta` and `\eta` must satisfy `0.25 < \delta
        \leq 1.0` and `0.5 \leq \eta < \sqrt{\delta}`. Polynomial time
        complexity is only guaranteed for `\delta < 1`. Not every algorithm
        admits the case `\delta = 1`.

        If the matrix has a nonzero kernel, the LLL-reduced matrix will contain
        zero rows, so that the output has the same dimensions as the input. The
        transformation matrix is always invertible over the integers.

        Also the rank of ``self`` is cached if it is computed during the
        reduction. Note that in general this only happens when
        ``self.rank() == self.ncols()`` and the exact algorithm is used.

        INPUT:

        - ``delta`` -- (default: ``0.99``) `\delta` parameter as described
          above, ignored by pari

        - ``eta`` -- (default: ``0.501``) `\eta` parameter as described above,
          ignored by NTL and pari

        - ``algorithm`` -- string; one of the algorithms listed below
          (default: ``'fpLLL:wrapper'``)

        - ``fp`` -- floating point number implementation, ignored by pari:

          - ``None`` -- NTL's exact reduction or fpLLL's wrapper
          - ``'fp'`` -- double precision: NTL's FP or fpLLL's double
          - ``'ld'`` -- long doubles (fpLLL only)
          - ``'qd'`` -- NTL's QP
          - ``'xd'`` -- extended exponent: NTL's XD or fpLLL's dpe
          - ``'rr'`` -- arbitrary precision: NTL's RR or fpLLL's MPFR

        - ``prec`` -- (default: auto choose) precision, ignored by NTL and pari

        - ``early_red`` -- boolean (default: ``False``); perform early reduction,
          ignored by NTL and pari

        - ``use_givens`` -- boolean (default: ``False``); use Givens
          orthogonalization.  Only applies to approximate reduction
          using NTL.  This is slower but generally more stable.

        - ``use_siegel`` -- boolean (default: ``False``); use Siegel's condition
          instead of Lovász's condition, ignored by NTL and pari

        - ``transformation`` -- boolean (default: ``False``); also return transformation
          matrix

        - ``**kwds`` -- keywords to be passed to :mod:`fpylll`; see
          :meth:`fpylll.LLL.reduction` for details

        Also, if the verbose level is at least `2`, some output
        is printed during the computation.

        AVAILABLE ALGORITHMS:

        - ``'NTL:LLL'`` -- NTL's LLL + choice of ``fp``

        - ``'fpLLL:heuristic'`` -- fpLLL's heuristic + choice of ``fp``

        - ``'fpLLL:fast'`` -- fpLLL's fast + choice of ``fp``

        - ``'fpLLL:proved'`` -- fpLLL's proved + choice of ``fp``

        - ``'fpLLL:wrapper'`` -- fpLLL's automatic choice (default)

        - ``'pari'`` -- pari's qflll

        OUTPUT: a matrix over the integers

        EXAMPLES::

            sage: A = Matrix(ZZ,3,3,range(1,10))
            sage: A.LLL()
            [ 0  0  0]
            [ 2  1  0]
            [-1  1  3]

        We compute the extended GCD of a list of integers using LLL, this
        example is from the Magma handbook::

            sage: Q = [ 67015143, 248934363018, 109210, 25590011055, 74631449,
            ....:       10230248, 709487, 68965012139, 972065, 864972271 ]
            sage: n = len(Q)
            sage: S = 100
            sage: X = Matrix(ZZ, n, n + 1)
            sage: for i in range(n):
            ....:     X[i, i + 1] = 1
            sage: for i in range(n):
            ....:     X[i, 0] = S * Q[i]
            sage: L = X.LLL()
            sage: M = L.row(n-1).list()[1:]
            sage: M
            [-3, -1, 13, -1, -4, 2, 3, 4, 5, -1]
            sage: add(Q[i]*M[i] for i in range(n))
            -1

        The case `\delta = 1` is not always supported::

            sage: L = X.LLL(delta=2)
            Traceback (most recent call last):
            ...
            TypeError: delta must be <= 1
            sage: L = X.LLL(delta=1)    # not tested, will eat lots of ram
            Traceback (most recent call last):
            ...
            RuntimeError: infinite loop in LLL
            sage: L = X.LLL(delta=1, algorithm='NTL:LLL')
            sage: L[-1]
            (-100, -3, -1, 13, -1, -4, 2, 3, 4, 5, -1)

        We return the transformation matrix::

            sage: A = random_matrix(ZZ, 10, 20)
            sage: R, U = A.LLL(transformation=True)
            sage: U * A == R
            True

            sage: R, U = A.LLL(algorithm='NTL:LLL', transformation=True)
            sage: U * A == R
            True

            sage: R, U = A.LLL(algorithm='pari', transformation=True)
            sage: U * A == R
            True

        Example with a nonzero kernel::

            sage: M = matrix(4,3,[1,2,3,2,4,6,7,0,1,-1,-2,-3])
            sage: M.LLL()[0:2]
            [0 0 0]
            [0 0 0]

            sage: M.LLL(algorithm="NTL:LLL")[0:2]
            [0 0 0]
            [0 0 0]

            sage: M.LLL(algorithm='pari')[0:2]
            [0 0 0]
            [0 0 0]

        TESTS::

            sage: matrix(ZZ, 0, 0).LLL()
            []
            sage: matrix(ZZ, 3, 0).LLL()
            []
            sage: matrix(ZZ, 0, 3).LLL()
            []

            sage: M = matrix(ZZ, [[1,2,3],[31,41,51],[101,201,301]])
            sage: A = M.LLL()
            sage: A
            [ 0  0  0]
            [-1  0  1]
            [ 1  1  1]
            sage: B = M.LLL(algorithm='NTL:LLL')
            sage: C = M.LLL(algorithm='NTL:LLL', fp=None)
            sage: D = M.LLL(algorithm='NTL:LLL', fp='fp')
            sage: F = M.LLL(algorithm='NTL:LLL', fp='xd')
            sage: G = M.LLL(algorithm='NTL:LLL', fp='rr')
            sage: A == B == C == D == F == G
            True
            sage: H = M.LLL(algorithm='NTL:LLL', fp='qd')
            Traceback (most recent call last):
            ...
            TypeError: algorithm NTL:LLL_QD not supported

            sage: A = random_matrix(ZZ, 0, 0)
            sage: R, U = A.LLL(transformation=True)

        Test rank caching::

            sage: M = matrix(4,3,[1,2,3,2,4,6,7,0,1,-1,-2,-3])
            sage: R = M.LLL(algorithm="NTL:LLL")
            sage: M._cache
            {'rank': 2}
            sage: M._clear_cache()
            sage: R = M.LLL(algorithm='pari')
            sage: M._cache
            {'rank': 2}

        Check that :issue:`37236` is fixed::

            sage: M = matrix(ZZ, 2, 2, [-1,1,1,1])
            sage: L = M.LLL(algorithm="NTL:LLL")
            sage: M.det() == L.det()
            True

        .. NOTE::

            See :mod:`sage.libs.ntl.ntl_mat_ZZ.ntl_mat_ZZ.LLL` and
            :mod:`fpylll.fplll.lll` for details on the algorithms used.

            Although LLL is a deterministic algorithm, the output for
            different implementations and CPUs (32-bit vs. 64-bit) may
            vary, while still being correct.
        """
        if self.ncols() == 0 or self.nrows() == 0:
            verbose("Trivial matrix, nothing to do")
            if transformation:
                return self, self.new_matrix(0,0)
            else:
                return self

        r = None  # rank
        cdef Matrix_integer_dense R = None  # LLL-reduced matrix
        cdef Matrix_integer_dense U = None  # transformation matrix

        tm = verbose("LLL of %sx%s matrix (algorithm %s)"%(self.nrows(), self.ncols(), algorithm))
        import sage.libs.ntl.all
        ntl_ZZ = sage.libs.ntl.all.ZZ

        verb = get_verbose() >= 2

        if prec < 0:
            raise TypeError("precision prec must be >= 0")
        prec = int(prec)

        if algorithm == 'NTL:LLL':
            if fp is None:
                algorithm = 'NTL:LLL'
            elif fp == 'fp':
                algorithm = 'NTL:LLL_FP'
            elif fp == 'qd':
                algorithm = 'NTL:LLL_QD'
            elif fp == 'xd':
                algorithm = 'NTL:LLL_XD'
            elif fp == 'rr':
                algorithm = 'NTL:LLL_RR'

            if delta is None:
                delta = ZZ(99)/ZZ(100)
            elif delta <= ZZ(1)/ZZ(4):
                raise TypeError("delta must be > 1/4")
            elif delta > 1:
                raise TypeError("delta must be <= 1")
            delta = QQ(delta)
            a = delta.numer()
            b = delta.denom()

        else:
            fp = fplll_fp_map[fp]
            if delta is None:
                delta = 0.99
            elif delta <= 0.25:
                raise TypeError("delta must be > 0.25")
            elif delta > 1:
                raise TypeError("delta must be <= 1")
            delta = float(delta)

            if eta is None:
                eta = 0.501
            elif eta < 0.5:
                raise TypeError("eta must be >= 0.5")

        if algorithm.startswith('NTL:'):

            A = sage.libs.ntl.all.mat_ZZ(self.nrows(),self.ncols(),
                    [ntl_ZZ(z) for z in self.list()])

            if algorithm == "NTL:LLL":
                if transformation:
                    r, _, UNTL = A.LLL(a, b, verbose=verb, return_U=True)
                else:
                    r, _ = A.LLL(a, b, verbose=verb)
            elif algorithm == "NTL:LLL_FP":
                if use_givens:
                    r = A.G_LLL_FP(delta, verbose=verb, return_U=transformation)
                else:
                    r = A.LLL_FP(delta, verbose=verb, return_U=transformation)
            elif algorithm == "NTL:LLL_QP":
                if use_givens:
                    r = A.G_LLL_QP(delta, verbose=verb, return_U=transformation)
                else:
                    r = A.LLL_QP(delta, verbose=verb, return_U=transformation)
            elif algorithm == "NTL:LLL_XD":
                if use_givens:
                    r = A.G_LLL_XD(delta, verbose=verb, return_U=transformation)
                else:
                    r = A.LLL_XD(delta, verbose=verb, return_U=transformation)
            elif algorithm == "NTL:LLL_RR":
                if use_givens:
                    r = A.G_LLL_RR(delta, verbose=verb, return_U=transformation)
                else:
                    r = A.LLL_XD(delta, verbose=verb, return_U=transformation)
            else:
                raise TypeError("algorithm %s not supported"%algorithm)

            if isinstance(r, tuple):
                r, UNTL = r

            if transformation:
                U = self.new_matrix(self.nrows(), self.nrows(),
                                    entries=[ZZ(z) for z in UNTL.list()])

            r = ZZ(r)

            R = <Matrix_integer_dense> self.new_matrix(
                    entries=[ZZ(z) for z in A.list()])

        elif algorithm.startswith('fpLLL:'):
            from fpylll import LLL, IntegerMatrix
            A = IntegerMatrix.from_matrix(self)
            Ufplll = None
            if transformation:
                Ufplll = IntegerMatrix(A.nrows, A.nrows)

            method = algorithm.replace("fpLLL:","")
            if verb:
                kwds["flags"] = kwds.get("flags", LLL.DEFAULT) | LLL.VERBOSE
            if use_siegel:
                kwds["flags"] = kwds.get("flags", LLL.DEFAULT) | LLL.SIEGEL
            if early_red:
                kwds["flags"] = kwds.get("flags", LLL.DEFAULT) | LLL.EARLY_RED

            LLL.reduction(A, delta=delta, eta=eta, method=method, float_type=fp, precision=prec, U=Ufplll)
            R = A.to_matrix(self.new_matrix())
            if transformation:
                U = Ufplll.to_matrix(self.new_matrix(Ufplll.nrows, Ufplll.ncols))

        elif algorithm == 'pari':
            # call pari with flag=4: kernel+image
            # pari uses column convention: need to transpose the matrices
            A = integer_matrix(self._matrix, 1)
            K, T = A.qflll(4)
            r = ZZ(T.length())
            # TODO: there is no optimized matrix converter pari -> sage
            U = gen_to_sage_matrix(pari.matconcat([K, T]).mattranspose())
            R = U * self

        else:
            raise TypeError("algorithm %s not supported"%algorithm)

        verbose("LLL finished", tm)
        if r is not None:
            self.cache("rank", r)
            R.cache("rank", r)
        if transformation:
            return R, U
        else:
            return R

    def is_LLL_reduced(self, delta=None, eta=None, algorithm='fpLLL'):
        r"""
        Return ``True`` if this lattice is `(\delta, \eta)`-LLL reduced.
        See ``self.LLL`` for a definition of LLL reduction.

        INPUT:

        - ``delta`` -- (default: `0.99`) parameter `\delta` as described above

        - ``eta`` -- (default: `0.501`) parameter `\eta` as described above

        - ``algorithm`` -- either ``'fpLLL'`` (default) or ``'sage'``

        EXAMPLES::

            sage: A = random_matrix(ZZ, 10, 10)
            sage: L = A.LLL()
            sage: A.is_LLL_reduced()
            False
            sage: L.is_LLL_reduced()
            True

        The ``'sage'`` algorithm currently does not work for matrices with
        linearly dependent rows::

            sage: A = matrix(ZZ, [[1, 2, 3], [2, 4, 6]])
            sage: A.is_LLL_reduced(algorithm='sage')
            Traceback (most recent call last):
            ...
            ValueError: linearly dependent input for module version of Gram-Schmidt
        """
        if eta is None:
            eta = 0.501
        if delta is None:
            delta = ZZ(99) / ZZ(100)

        if delta <= ZZ(1) / ZZ(4):
            raise TypeError("delta must be > 1/4")
        elif delta > 1:
            raise TypeError("delta must be <= 1")

        if eta < 0.5:
            raise TypeError("eta must be >= 0.5")

        if algorithm == 'fpLLL':
            from fpylll import LLL, IntegerMatrix
            A = IntegerMatrix.from_matrix(self)
            return LLL.is_reduced(A, delta=delta, eta=eta)
        elif algorithm == 'sage':
            # This is pretty slow
            import sage.modules.misc
            G, mu = sage.modules.misc.gram_schmidt(self.rows())
            # For any $i>j$, we have $|mu_{i, j}| <= \eta$
            for e in mu.list():
                if e.abs() > eta:
                    return False

            # For any $i<d$, we have $\delta |b_i^*|^2 <= |b_{i+1}^* + mu_{i+1, i} b_i^* |^2$
            norms = [G[i].norm()**2 for i in range(len(G))]
            for i in range(1,self.nrows()):
                if norms[i] < (delta - mu[i,i-1]**2) * norms[i-1]:
                    return False
            return True
        else:
            raise ValueError("algorithm must be one of 'fpLLL' or 'sage'")

    def prod_of_row_sums(self, cols):
        """
        Return the product of the sums of the entries in the submatrix of
        ``self`` with given columns.

        INPUT:

        - ``cols`` -- list (or set) of integers representing columns
          of ``self``

        OUTPUT: integer

        EXAMPLES::

            sage: a = matrix(ZZ,2,3,[1..6]); a
            [1 2 3]
            [4 5 6]
            sage: a.prod_of_row_sums([0,2])
            40
            sage: (1+3)*(4+6)
            40
            sage: a.prod_of_row_sums(set([0,2]))
            40
        """
        cdef Py_ssize_t c, row
        cdef fmpz_t s,pr
        fmpz_init(s)
        fmpz_init_set_ui(pr, 1)

        for row from 0 <= row < self._nrows:
            fmpz_set_si(s, 0)
            for c in cols:
                if c<0 or c >= self._ncols:
                    fmpz_clear(s)
                    fmpz_clear(pr)
                    raise IndexError("matrix column index out of range")
                fmpz_add(s, s, fmpz_mat_entry(self._matrix,row,c))
            fmpz_mul(pr, pr, s)
        cdef Integer z = Integer.__new__(Integer)
        fmpz_get_mpz(z.value, pr)
        fmpz_clear(s)
        fmpz_clear(pr)
        return z

    def rational_reconstruction(self, N):
        r"""
        Use rational reconstruction to lift ``self`` to a matrix over the
        rational numbers (if possible), where we view ``self`` as a matrix
        modulo N.

        INPUT:

        - ``N`` -- integer

        OUTPUT: matrix over `\QQ` or raise a :exc:`ValueError`

        EXAMPLES: We create a random 4x4 matrix over ZZ.

        ::

            sage: A = matrix(ZZ, 4, [4, -4, 7, 1, -1, 1, -1, -12, -1, -1, 1, -1, -3, 1, 5, -1])

        There isn't a unique rational reconstruction of it::

            sage: A.rational_reconstruction(11)
            Traceback (most recent call last):
            ...
            ValueError: rational reconstruction does not exist

        We throw in a denominator and reduce the matrix modulo 389 - it
        does rationally reconstruct.

        ::

            sage: B = (A/3 % 389).change_ring(ZZ)
            sage: B.rational_reconstruction(389) == A/3
            True

        TESTS:

        Check that :issue:`9345` is fixed::

            sage: A = random_matrix(ZZ, 3, 3)
            sage: A.rational_reconstruction(0)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: The modulus cannot be zero
        """
        from sage.matrix.misc_flint import matrix_integer_dense_rational_reconstruction
        return matrix_integer_dense_rational_reconstruction(self, N)

    def randomize(self, density=1, x=None, y=None, distribution=None,
                  nonzero=False):
        """
        Randomize ``density`` proportion of the entries of this matrix,
        leaving the rest unchanged.

        The parameters are the same as the ones for the integer ring's
        ``random_element`` function.

        If ``x`` and ``y`` are given, randomized entries of this matrix have
        to be between ``x`` and ``y`` and have density 1.

        INPUT:

        - ``self`` -- a mutable matrix over ZZ

        - ``density`` -- a float between 0 and 1

        - ``x``, ``y`` -- if not ``None``, these are passed to the
           ``ZZ.random_element`` function as the upper and lower endpoints in
           the  uniform distribution

        - ``distribution`` -- would also be passed into ``ZZ.random_element``
          if given

        - ``nonzero`` -- boolean (default: ``False``); whether the new entries
          are guaranteed to be zero

        OUTPUT: None, the matrix is modified in-place

        EXAMPLES::

            sage: A = matrix(ZZ, 2,3, [1..6])
            sage: ranks = [True, True, True]
            sage: while any(ranks):
            ....:    A.randomize()
            ....:    ranks[A.rank()] = False

            sage: mini = 0
            sage: maxi = 0
            sage: while mini != -30 and maxi != 30:
            ....:     A.randomize(x=-30, y=30)
            ....:     mini = min(min(A.list()), mini)
            ....:     maxi = min(min(A.list()), maxi)
        """
        density = float(density)
        if density <= 0:
            return
        if density > 1:
            density = float(1)

        self.check_mutability()
        self.clear_cache()

        cdef randstate rstate = current_randstate()

        cdef mpz_t tmp
        mpz_init(tmp)
        cdef Py_ssize_t i, j, k, nc, num_per_row
        global state, ZZ

        cdef IntegerRing_class the_integer_ring = ZZ

        if not nonzero:
            # Original code, before adding the ``nonzero`` option.
            sig_on()
            if density == 1:
                for i from 0 <= i < self._nrows:
                    for j from 0 <= j < self._ncols:
                        the_integer_ring._randomize_mpz(tmp, x, y,
                                                        distribution)
                        self.set_unsafe_mpz(i,j,tmp)
            else:
                nc = self._ncols
                num_per_row = int(density * nc)
                for i from 0 <= i < self._nrows:
                    for j from 0 <= j < num_per_row:
                        k = rstate.c_random()%nc
                        the_integer_ring._randomize_mpz(tmp,
                                                        x, y, distribution)
                        self.set_unsafe_mpz(i,k,tmp)
            sig_off()
        else:
            # New code, to implement the ``nonzero`` option.  Note that this
            # code is almost the same as above, the only difference being that
            # each entry is set until it's nonzero.
            sig_on()
            if density == 1:
                for i from 0 <= i < self._nrows:
                    for j from 0 <= j < self._ncols:
                        while fmpz_sgn(fmpz_mat_entry(self._matrix,i,j)) == 0:
                            the_integer_ring._randomize_mpz(tmp,
                                x, y, distribution)
                            self.set_unsafe_mpz(i,j,tmp)
            else:
                nc = self._ncols
                num_per_row = int(density * nc)
                for i from 0 <= i < self._nrows:
                    for j from 0 <= j < num_per_row:
                        k = rstate.c_random() % nc
                        while fmpz_sgn(fmpz_mat_entry(self._matrix,i,k)) == 0:
                            the_integer_ring._randomize_mpz(tmp,
                                                            x, y, distribution)
                            self.set_unsafe_mpz(i,k,tmp)

            sig_off()
        mpz_clear(tmp)

    #### Rank

    def rank(self, algorithm='modp'):
        """
        Return the rank of this matrix.

        INPUT:

        - ``algorithm`` -- either ``'modp'`` (default) or ``'flint'``
          or ``'linbox'``

        OUTPUT: nonnegative integer -- the rank

        .. NOTE::

            The rank is cached.

        ALGORITHM:

        If set to ``'modp'``, first check if the matrix has maximum
        possible rank by working modulo one random prime. If not, call
        LinBox's rank function.

        EXAMPLES::

            sage: a = matrix(ZZ,2,3,[1..6]); a
            [1 2 3]
            [4 5 6]
            sage: a.rank()
            2
            sage: a = matrix(ZZ,3,3,[1..9]); a
            [1 2 3]
            [4 5 6]
            [7 8 9]
            sage: a.rank()
            2

        Here is a bigger example - the rank is of course still 2::

            sage: a = matrix(ZZ,100,[1..100^2]); a.rank()
            2

        TESTS::

            sage: a.rank(algorithm='funky')
            Traceback (most recent call last):
            ...
            ValueError: algorithm must be one of 'modp', 'flint' or 'linbox'
        """
        if algorithm not in ['modp', 'flint', 'linbox']:
            raise ValueError("algorithm must be one of 'modp', 'flint' "
                             "or 'linbox'")

        r = self.fetch('rank')
        if r is not None:
            return r

        if algorithm == 'flint' or (self._nrows <= 6 and self._ncols <= 6
                                    and self.height().ndigits() <= 100):
            r = fmpz_mat_rank(self._matrix)
            self.cache('rank', r)
            return r
        elif algorithm == 'modp':
            # Can very quickly detect full rank by working modulo p.
            r = self._rank_modp()
            if r == self._nrows or r == self._ncols:
                self.cache('rank', r)
                return r
        # Algorithm is 'linbox' or detecting full rank didn't work --
        # use LinBox's general algorithm.
        r = self._rank_linbox()
        self.cache('rank', r)
        return r

    def _rank_linbox(self):
        """
        Compute the rank of this matrix using Linbox.

        TESTS::

            sage: matrix(ZZ, 4, 6, 0)._rank_linbox()
            0
            sage: matrix(ZZ, 3, 4, range(12))._rank_linbox()
            2
            sage: matrix(ZZ, 5, 10, [1+i+i^2 for i in range(50)])._rank_linbox()
            3
        """
        sig_on()
        cdef size_t r = linbox_fmpz_mat_rank(self._matrix)
        sig_off()
        return Integer(r)

    def _rank_modp(self, p=46337):
        A = self._mod_int_c(p)
        return A.rank()

    #### Determinant

    def determinant(self, algorithm='default', proof=None, stabilize=2):
        r"""
        Return the determinant of this matrix.

        INPUT:

        - ``algorithm``

          - ``'default'`` -- use ``flint``

          - ``'flint'`` -- let flint do the determinant

          - ``'padic'`` -- uses a `p`-adic / multimodular
            algorithm that relies on code in IML and linbox

          - ``'linbox'`` -- calls linbox det (you *must* set
            proof=False to use this!)

          - ``'ntl'`` -- calls NTL's det function

          - ``'pari'`` -- uses PARI

        - ``proof`` -- boolean or ``None``; if ``None`` use
          proof.linear_algebra(); only relevant for the padic algorithm

           .. NOTE::

              It would be *VERY VERY* hard for det to fail even with
              proof=False.

        - ``stabilize`` -- if proof is ``False``, require det to be the same
          for this many CRT primes in a row. Ignored if proof is ``True``.


        ALGORITHM: The `p`-adic algorithm works by first finding a random
        vector v, then solving `Ax = v` and taking the denominator
        `d`. This gives a divisor of the determinant. Then we
        compute `\det(A)/d` using a multimodular algorithm and the
        Hadamard bound, skipping primes that divide `d`.

        EXAMPLES::

            sage: A = matrix(ZZ,8,8,[3..66])
            sage: A.determinant()
            0

        ::

            sage: A = random_matrix(ZZ,20,20)
            sage: D1 = A.determinant()
            sage: A._clear_cache()
            sage: D2 = A.determinant(algorithm='ntl')
            sage: D1 == D2
            True

        We have a special-case algorithm for 4 x 4 determinants::

            sage: A = matrix(ZZ,4,[1,2,3,4,4,3,2,1,0,5,0,1,9,1,2,3])
            sage: A.determinant()
            270

        Next we try the Linbox det. Note that we must have proof=False.

        ::

            sage: A = matrix(ZZ,5,[1,2,3,4,5,4,6,3,2,1,7,9,7,5,2,1,4,6,7,8,3,2,4,6,7])
            sage: A.determinant(algorithm='linbox')
            Traceback (most recent call last):
            ...
            RuntimeError: you must pass the proof=False option to the determinant command to use LinBox's det algorithm
            sage: A.determinant(algorithm='linbox', proof=False)
            -21
            sage: A._clear_cache()
            sage: A.determinant()
            -21

        Try the other algorithms on the same example::

            sage: A._clear_cache(); A.determinant(algorithm='padic')
            -21
            sage: A._clear_cache(); A.determinant(algorithm='pari')
            -21
            sage: A._clear_cache(); A.determinant(algorithm='ntl')
            -21
            sage: A._clear_cache(); A.determinant(algorithm='padic')
            -21

        A bigger example::

            sage: A = random_matrix(ZZ,30)
            sage: d = A.determinant()
            sage: A._clear_cache()
            sage: A.determinant(algorithm='linbox',proof=False) == d
            True

        TESTS:

        This shows that we can compute determinants for all sizes up to
        80. The check that the determinant of a squared matrix is a
        square is a sanity check that the result is probably correct::

            sage: for s in [1..80]:  # long time
            ....:     M = random_matrix(ZZ, s)
            ....:     d = (M*M).determinant()
            ....:     assert d.is_square()

        Check consistency::

            sage: all(matrix(ZZ, 0).det(algorithm=algo).is_one() for algo in ['flint', 'padic', 'pari', 'ntl'])
            True
            sage: for _ in range(100):
            ....:     dim = randint(1, 10)
            ....:     m = random_matrix(ZZ, dim)
            ....:     det_flint = m.__copy__().det(algorithm='flint')
            ....:     det_padic = m.__copy__().det(algorithm='padic')
            ....:     det_pari = m.__copy__().det(algorithm='pari')
            ....:     det_ntl = m.__copy__().det(algorithm='ntl')
            ....:     if type(det_flint) is not Integer:
            ....:         raise RuntimeError("type(det_flint) = {}".format(type(det_flint)))
            ....:     if type(det_padic) is not Integer:
            ....:         raise RuntimeError("type(det_padic) = {}".format(type(det_padic)))
            ....:     if type(det_pari) is not Integer:
            ....:         raise RuntimeError("type(det_pari) = {}".format(type(det_pari)))
            ....:     if type(det_ntl) is not Integer:
            ....:         raise RuntimeError("type(det_ntl) = {}".format(type(det_ntl)))
            ....:     if det_flint != det_padic or det_flint != det_pari or det_flint != det_ntl:
            ....:         raise RuntimeError("ERROR\ndet_flint = {}\ndet_padic={}\ndet_pari={}\ndet_ntl={}\n{}".format(
            ....:                           det_flint, det_padic, det_pari, det_ntl, self.str()))
        """
        d = self.fetch('det')
        if d is not None:
            return d
        if not self.is_square():
            raise ValueError("self must be a square matrix")

        cdef Integer det = Integer()
        cdef fmpz_t e

        proof = get_proof_flag(proof, "linear_algebra")

        if algorithm == 'default':
            algorithm = 'flint'

        if algorithm == 'flint':
            fmpz_init(e)
            fmpz_mat_det(e, self._matrix)
            fmpz_get_mpz(det.value, e)
            fmpz_clear(e)
            d = det
        elif algorithm == 'padic':
            from sage.matrix import matrix_integer_dense_hnf
            d = matrix_integer_dense_hnf.det_padic(self, proof=proof, stabilize=stabilize)
        elif algorithm == 'linbox':
            if proof:
                raise RuntimeError("you must pass the proof=False option to the determinant command to use LinBox's det algorithm")
            d = self._det_linbox()
        elif algorithm == 'pari':
            d = self._det_pari()
        elif algorithm == 'ntl':
            d = self._det_ntl()
        else:
            raise TypeError("algorithm '%s' not known"%(algorithm))

        self.cache('det', d)
        return d

    def _det_linbox(self):
        """
        Compute the determinant of this matrix using Linbox.

        TESTS::

            sage: matrix(ZZ, 0)._det_linbox()
            1
        """
        if self._nrows != self._ncols:
            raise ArithmeticError("self must be a square matrix")
        if self._nrows == 0:
            return ZZ.one()

        cdef fmpz_t tmp
        fmpz_init(tmp)
        sig_on()
        linbox_fmpz_mat_det(tmp, self._matrix)
        sig_off()

        cdef Integer ans = PY_NEW(Integer)
        fmpz_get_mpz(ans.value, tmp)
        fmpz_clear(tmp)
        return ans

    def _det_pari(self, int flag=0):
        """
        Determinant of this matrix using Gauss-Bareiss. If (optional)
        flag is set to 1, use classical Gaussian elimination.

        For efficiency purposes, this det is computed entirely on the
        PARI stack then the PARI stack is cleared. This function is
        most useful for very small matrices.

        EXAMPLES::

            sage: matrix(ZZ, 0)._det_pari()
            1
            sage: matrix(ZZ, 0)._det_pari(1)
            1
            sage: matrix(ZZ, 3, [1..9])._det_pari()
            0
            sage: matrix(ZZ, 3, [1..9])._det_pari(1)
            0
        """
        sig_on()
        cdef GEN d = det0(pari_GEN(self), flag)
        # now convert d to a Sage integer e
        cdef Integer e = Integer.__new__(Integer)
        INT_to_mpz(e.value, d)
        clear_stack()
        return e

    def _det_ntl(self):
        """
        Compute the determinant of this matrix using NTL.

        EXAMPLES::

            sage: matrix(ZZ, 0)._det_ntl()
            1
            sage: matrix(ZZ, 3, [1..9])._det_ntl()
            0
            sage: matrix(ZZ, 3, [1,3,6,2,7,8,2,1,0])._det_ntl()
            -32
        """
        sig_on()
        d = self._ntl_().determinant()
        sig_off()
        return Integer(d)

    def _rational_kernel_iml(self):
        """
        Return the rational (left) kernel of this matrix.

        OUTPUT:

        A matrix ``K`` such that ``self * K = 0``, and the number of columns of
        K equals the nullity of ``self``.

        EXAMPLES::

            sage: m = matrix(ZZ, 5, 5, [1+i+i^2 for i in range(25)])
            sage: m._rational_kernel_iml()
            [ 1  3]
            [-3 -8]
            [ 3  6]
            [-1  0]
            [ 0 -1]

            sage: V1 = m._rational_kernel_iml().column_space().change_ring(QQ)
            sage: V2 = m._rational_kernel_flint().column_space().change_ring(QQ)
            sage: assert V1 == V2
        """
        if self._nrows == 0 or self._ncols == 0:
            return self.matrix_space(self._ncols, 0).zero_matrix()

        cdef long dim
        cdef unsigned long i,j,k
        cdef mpz_t *mp_N
        time = verbose('computing null space of %s x %s matrix using IML'%(self._nrows, self._ncols))
        cdef mpz_t * m = fmpz_mat_to_mpz_array(self._matrix)
        sig_on()
        dim = nullspaceMP(self._nrows, self._ncols, m, &mp_N)
        sig_off()
        # Now read the answer as a matrix.
        cdef Matrix_integer_dense M
        M = self._new(self._ncols, dim)
        k = 0
        for i in range(self._ncols):
            for j in range(dim):
                fmpz_set_mpz(fmpz_mat_entry(M._matrix, i, j), mp_N[k])
                k += 1
        mpz_array_clear(m, self._nrows * self._ncols)
        mpz_array_clear(mp_N, dim * self._ncols)
        verbose("finished computing null space", time)
        return M

    def _rational_kernel_flint(self):
        """
        Return the rational (left) kernel of this matrix.

        OUTPUT:

        A matrix ``K`` such that ``self * K = 0``, and the number of columns of
        K equals the nullity of ``self``.

        EXAMPLES::

            sage: m = matrix(ZZ, 4, 6, [i^2-2*i for i in range(24)])
            sage: m._rational_kernel_flint()
            [  1728   5184  10368]
            [ -5184 -13824 -25920]
            [  5184  10368  17280]
            [ -1728      0      0]
            [     0  -1728      0]
            [     0      0  -1728]

            sage: V1 = m._rational_kernel_iml().column_space().change_ring(QQ)
            sage: V2 = m._rational_kernel_flint().column_space().change_ring(QQ)
            sage: assert V1 == V2
        """
        if self._nrows == 0 or self._ncols == 0:
            return self.matrix_space(self._ncols, 0).zero_matrix()

        cdef long dim
        cdef fmpz_mat_t M0
        sig_on()
        fmpz_mat_init(M0, self._ncols, self._ncols)
        dim = fmpz_mat_nullspace(M0, self._matrix)
        sig_off()
        # Now read the answer as a matrix.
        cdef Matrix_integer_dense M = self._new(self._ncols, dim)
        cdef size_t i,j
        for i in range(self._ncols):
            for j in range(dim):
                fmpz_set(fmpz_mat_entry(M._matrix, i, j), fmpz_mat_entry(M0, i, j))
        fmpz_mat_clear(M0)
        return M

    def _invert_iml(self, use_nullspace=False, check_invertible=True):
        """
        Invert this matrix using IML. The output matrix is an integer
        matrix and a denominator.

        INPUT:

        - ``self`` -- an invertible matrix

        - ``use_nullspace`` -- boolean (default: ``False``); whether to
          use nullspace algorithm, which is slower, but doesn't require
          checking that the matrix is invertible as a precondition

        - ``check_invertible`` -- boolean (default: ``True``); whether to
          check that the matrix is invertible

        OUTPUT: `A`, `d` such that ``A*self == d``

        - ``A`` -- a matrix over ZZ

        - ``d`` -- integer

        ALGORITHM: Uses IML's `p`-adic nullspace function.

        EXAMPLES::

            sage: a = matrix(ZZ,3,[1,2,5, 3,7,8, 2,2,1])
            sage: b, d = a._invert_iml(); b,d
            (
            [  9  -8  19]
            [-13   9  -7]
            [  8  -2  -1], 23
            )
            sage: a*b
            [23  0  0]
            [ 0 23  0]
            [ 0  0 23]

        AUTHORS:

        - William Stein
        """
        if self._nrows != self._ncols:
            raise TypeError("self must be a square matrix.")

        P = self.parent()
        time = verbose('computing inverse of %s x %s matrix using IML'%(self._nrows, self._ncols))
        if use_nullspace:
            A = self.augment(P.identity_matrix())
            K = A._rational_kernel_iml()
            d = -K[self._nrows,0]
            if K.ncols() != self._ncols or d == 0:
                raise ZeroDivisionError("input matrix must be nonsingular")
            B = K[:self._nrows]
            verbose("finished computing inverse using IML", time)
            return B, d
        else:
            if check_invertible and self.rank() != self._nrows:
                raise ZeroDivisionError("input matrix must be nonsingular")
            return self._solve_iml(P.identity_matrix(), right=True)

    def _invert_flint(self):
        """
        Invert this matrix using FLINT. The output matrix is an integer
        matrix and a denominator.

        INPUT:

        - ``self`` -- an invertible matrix

        OUTPUT: `A`, `d` such that ``A*self == d``

        - ``A`` -- a matrix over ZZ

        - ``d`` -- integer

        EXAMPLES::

            sage: a = matrix(ZZ,3,[1,2,5, 3,7,8, 2,2,1])
            sage: b, d = a._invert_flint(); b,d
            (
            [  9  -8  19]
            [-13   9  -7]
            [  8  -2  -1], 23
            )
            sage: a*b
            [23  0  0]
            [ 0 23  0]
            [ 0  0 23]

        AUTHORS:

        - William Stein

        - Marc Masdeu -- (08/2014) Use FLINT
        """
        if self._nrows != self._ncols:
            raise ArithmeticError("self must be a square matrix.")

        cdef Matrix_integer_dense M
        cdef int res
        cdef Integer den = Integer(0)
        cdef fmpz_t fden
        fmpz_init(fden)
        M = self._new(self._nrows,self._ncols)
        verbose('computing inverse of %s x %s matrix using FLINT'%(self._nrows, self._ncols))
        sig_on()
        res = fmpz_mat_inv(M._matrix,fden,self._matrix)
        fmpz_get_mpz(den.value,fden)
        sig_off()
        fmpz_clear(fden)
        if res == 0:
            raise ZeroDivisionError('matrix must be nonsingular')
        if den < 0:
            return -M, -den
        else:
            return M, den

    def __invert__(self):
        r"""
        Return the inverse of ``self``.

        EXAMPLES::

            sage: M = MatrixSpace(ZZ,3)
            sage: m = M([1,2,3,3,4,5,1,2,-3])
            sage: ~m
            [-11/6     1  -1/6]
            [  7/6  -1/2   1/3]
            [  1/6     0  -1/6]
            sage: ~m * m == m * ~m == M.identity_matrix()
            True

        Note that inverse of determinant one integer matrices do not belong to
        the same parent::

            sage: (~M.identity_matrix()).parent()
            Full MatrixSpace of 3 by 3 dense matrices over Rational Field

        This is consistent with::

            sage: (~1).parent()
            Rational Field

        TESTS::

            sage: ~M.zero_matrix()
            Traceback (most recent call last):
            ...
            ZeroDivisionError: matrix must be nonsingular
        """
        A, d = self._invert_flint()
        return A / d

    def inverse_of_unit(self):
        r"""
        If ``self`` is a matrix with determinant `1` or `-1` return the inverse of
        ``self`` as a matrix over `ZZ`.

        EXAMPLES::

            sage: m = matrix(ZZ, 2, [2,1,1,1]).inverse_of_unit()
            sage: m
            [ 1 -1]
            [-1  2]
            sage: parent(m)
            Full MatrixSpace of 2 by 2 dense matrices over Integer Ring

            sage: matrix(2, [2,1,0,1]).inverse_of_unit()
            Traceback (most recent call last):
            ...
            ArithmeticError: non-invertible matrix
        """
        A, d = self._invert_flint()
        if not d.is_one():
            raise ArithmeticError("non-invertible matrix")
        return A

    def _solve_right_nonsingular_square(self, B, check_rank=True, algorithm='iml'):
        r"""
        If ``self`` is a matrix `A` of full rank, then this function
        returns a vector or matrix `X` such that `A X = B`.
        If `B` is a vector then `X` is a vector and if
        `B` is a matrix, then `X` is a matrix. The base
        ring of `X` is the integers unless a denominator is needed
        in which case the base ring is the rational numbers.

        .. NOTE::

           DEPRECATED. In Sage one can also write ``A \ B`` for
           ``A.solve_right(B)``, i.e., Sage implements the "the
           MATLAB/Octave backslash operator".

        .. NOTE::

           This is currently only implemented when A is square.

        INPUT:

        - ``B`` -- a matrix or vector

        - ``check_rank`` -- boolean (default: ``True``); if ``True``
          verify that in fact the rank is full

        - ``algorithm`` -- ``'iml'`` (default) or ``'flint'``

        OUTPUT: a matrix or vector over `\QQ`

        EXAMPLES::

            sage: a = matrix(ZZ, 2, [0, -1, 1, 0])
            sage: v = vector(ZZ, [2, 3])
            sage: a.solve_right(v)
            (3, -2)

        Note that the output vector or matrix is always over
        `\QQ`.

        ::

            sage: parent(a.solve_right(v))
            Vector space of dimension 2 over Rational Field

        We solve a bigger system where the answer is over the rationals.

        ::

            sage: a = matrix(ZZ, 3, 3, [1,2,3,4, 5, 6, 8, -2, 3])
            sage: v = vector(ZZ, [1,2,3])
            sage: w = a.solve_right(v); w
            (2/15, -4/15, 7/15)
            sage: parent(w)
            Vector space of dimension 3 over Rational Field
            sage: a * w
            (1, 2, 3)

        We solve a system where the right hand matrix has multiple
        columns.

        ::

            sage: a = matrix(ZZ, 3, 3, [1,2,3,4, 5, 6, 8, -2, 3])
            sage: b = matrix(ZZ, 3, 2, [1,5, 2, -3, 3, 0])
            sage: w = a.solve_right(b); w
            [ 2/15 -19/5]
            [-4/15 -27/5]
            [ 7/15 98/15]
            sage: a * w
            [ 1  5]
            [ 2 -3]
            [ 3  0]

        TESTS:

        We create a random 100x100 matrix and solve the
        corresponding system, then verify that the result is correct.
        (Note that this test is very risky without having a seeded
        random number generator!)

        ::

            sage: n = 100
            sage: a = random_matrix(ZZ,n)
            sage: v = vector(ZZ,n,range(n))
            sage: x = a.solve_right(v)
            sage: a * x == v
            True

            sage: n = 100
            sage: a = random_matrix(ZZ,n)
            sage: v = vector(ZZ,n,range(n))
            sage: x = a._solve_right_nonsingular_square(v,algorithm = 'flint')
            sage: a * x == v
            True
        """
        t = verbose('starting %s solve_right...' % algorithm)

        # It would probably be much better to rewrite linbox so it
        # throws an error instead of ** going into an infinite loop **
        # in the non-full rank case.  In any case, we do this for now,
        # since rank is very fast and infinite loops are evil.
        if check_rank and self.rank() < self.nrows():
            from sage.matrix.matrix2 import NotFullRankError
            raise NotFullRankError

        if not self.is_square():
            raise NotImplementedError("the input matrix must be square.")

        if isinstance(B, Vector):
            if self.nrows() != B.degree():
                raise ValueError("number of rows of self must equal degree of B.")
        elif self.nrows() != B.nrows():
                raise ValueError("number of rows of self must equal number of rows of B.")

        if self.nrows() == 0:
            return B

        matrix = True
        C = B
        if not isinstance(B, Matrix_integer_dense):
            if isinstance(B, Vector):
                matrix = False
                C = self.matrix_space(self.nrows(), 1)(B.list())
            else:
                raise NotImplementedError

        if C.ncols() >= 2*self.ncols():
            # likely much better to just invert then multiply
            X = self**(-1)*C
            verbose('finished solve_right (via inverse)', t)
            return X

        if algorithm == 'flint':
            X, d = self._solve_flint(C, right=True)
        elif algorithm == 'iml': # iml
            X, d = self._solve_iml(C, right = True)
        else:
            raise ValueError("Unknown algorithm '%s'"%algorithm)
        if d != 1:
            X = (1/d) * X
        if not matrix:
            # Convert back to a vector
            X = (X.base_ring() ** X.nrows())(X.list())
        verbose('finished solve_right via %s'%algorithm, t)
        return X

    def _solve_iml(self, Matrix_integer_dense B, right=True):
        """
        Let A equal ``self`` be a square matrix. Given B return an integer
        matrix C and an integer d such that ``self`` ``C*A == d*B`` if right is
        False or ``A*C == d*B`` if right is True.

        OUTPUT:

        - ``C`` -- integer matrix

        - ``d`` -- integer denominator

        EXAMPLES::

            sage: A = matrix(ZZ,4,4,[0, 1, -2, -1, -1, 1, 0, 2, 2, 2, 2, -1, 0, 2, 2, 1])
            sage: B = matrix(ZZ,3,4, [-1, 1, 1, 0, 2, 0, -2, -1, 0, -2, -2, -2])
            sage: C,d = A._solve_iml(B,right=False); C
            [  6 -18 -15  27]
            [  0  24  24 -36]
            [  4 -12  -6  -2]

        ::

            sage: d
            12

        ::

            sage: C*A == d*B
            True

        ::

            sage: A = matrix(ZZ,4,4,[0, 1, -2, -1, -1, 1, 0, 2, 2, 2, 2, -1, 0, 2, 2, 1])
            sage: B = matrix(ZZ,4,3, [-1, 1, 1, 0, 2, 0, -2, -1, 0, -2, -2, -2])
            sage: C,d = A._solve_iml(B)
            sage: C
            [ 12  40  28]
            [-12  -4  -4]
            [ -6 -25 -16]
            [ 12  34  16]

        ::

            sage: d
            12

        ::

            sage: A*C == d*B
            True

        Test wrong dimensions::

            sage: A = random_matrix(ZZ, 4, 4)
            sage: B = random_matrix(ZZ, 2, 3)
            sage: B._solve_iml(A)
            Traceback (most recent call last):
            ...
            ValueError: self must be a square matrix
            sage: A._solve_iml(B, right=False)
            Traceback (most recent call last):
            ...
            ArithmeticError: B's number of columns must match self's number of rows
            sage: A._solve_iml(B, right=True)
            Traceback (most recent call last):
            ...
            ArithmeticError: B's number of rows must match self's number of columns

        Check that this can be interrupted properly (:issue:`15453`)::

            sage: A = random_matrix(ZZ, 2000, 2000)
            sage: B = random_matrix(ZZ, 2000, 2000)
            sage: from sage.doctest.util import ensure_interruptible_after
            sage: with ensure_interruptible_after(2, max_wait_after_interrupt=8): A._solve_iml(B)

        ALGORITHM: Uses IML.

        AUTHORS:

        - Martin Albrecht
        """
        cdef unsigned long i, j, k
        cdef mpz_t *mp_N
        cdef mpz_t mp_D
        cdef Matrix_integer_dense M
        cdef Integer D

        if self._nrows != self._ncols:
            # This is *required* by the IML function we call below.
            raise ValueError("self must be a square matrix")

        if self._nrows == 1:
            return B, self[0,0]

        cdef SOLU_POS solu_pos

        if right:
            if self._ncols != B._nrows:
                raise ArithmeticError("B's number of rows must match self's number of columns")

            n = self._ncols
            m = B._ncols

            P = self.matrix_space(n, m)
            if self._nrows == 0 or self._ncols == 0:
                return P.zero_matrix(), Integer(1)

            if m == 0 or n == 0:
                return self.new_matrix(nrows = n, ncols = m), Integer(1)

            solu_pos = RightSolu

        else: # left
            if self._nrows != B._ncols:
                raise ArithmeticError("B's number of columns must match self's number of rows")

            n = self._ncols
            m = B._nrows

            P = self.matrix_space(m, n)
            if self._nrows == 0 or self._ncols == 0:
                return P.zero_matrix(), Integer(1)

            if m == 0 or n == 0:
                return self.new_matrix(nrows = m, ncols = n), Integer(1)

            solu_pos = LeftSolu

        sig_check()
        verbose("Initializing mp_N and mp_D")
        mp_N = <mpz_t *> sig_malloc( n * m * sizeof(mpz_t) )
        for i in range(n * m):
            mpz_init(mp_N[i])
        mpz_init(mp_D)
        verbose("Done with initializing mp_N and mp_D")
        cdef mpz_t * mself = fmpz_mat_to_mpz_array(self._matrix)
        cdef mpz_t * mB = fmpz_mat_to_mpz_array(B._matrix)
        try:
            verbose('Calling solver n = %s, m = %s'%(n,m))
            sig_on()
            nonsingSolvLlhsMM(solu_pos, n, m, mself, mB, mp_N, mp_D)
            sig_off()
            M = self._new(P.nrows(), P.ncols())
            k = 0
            for i from 0 <= i < M._nrows:
                for j from 0 <= j < M._ncols:
                    fmpz_set_mpz(fmpz_mat_entry(M._matrix,i,j), mp_N[k])
                    k += 1
            D = Integer.__new__(Integer)
            mpz_set(D.value, mp_D)
            return M, D
        finally:
            mpz_clear(mp_D)
            mpz_array_clear(mself, self.nrows() * self.ncols())
            mpz_array_clear(mB, B.nrows() * B.ncols())
            mpz_array_clear(mp_N, n*m)

    def _solve_flint(self, Matrix_integer_dense B, right=True):
        """
        Let A equal ``self`` be a square matrix. Given B return an integer
        matrix C and an integer d such that ``self`` ``C*A == d*B`` if right is
        ``False`` or ``A*C == d*B`` if right is ``True``.

        OUTPUT:

        - ``C`` -- integer matrix

        - ``d`` -- integer denominator

        EXAMPLES::

            sage: A = matrix(ZZ,4,4,[0, 1, -2, -1, -1, 1, 0, 2, 2, 2, 2, -1, 0, 2, 2, 1])
            sage: B = matrix(ZZ,3,4, [-1, 1, 1, 0, 2, 0, -2, -1, 0, -2, -2, -2])
            sage: C,d = A._solve_flint(B,right=False); C
            [  6 -18 -15  27]
            [  0  24  24 -36]
            [  4 -12  -6  -2]

        ::

            sage: d
            12

        ::

            sage: C*A == d*B
            True

        ::

            sage: A = matrix(ZZ,4,4,[0, 1, -2, -1, -1, 1, 0, 2, 2, 2, 2, -1, 0, 2, 2, 1])
            sage: B = matrix(ZZ,4,3, [-1, 1, 1, 0, 2, 0, -2, -1, 0, -2, -2, -2])
            sage: C,d = A._solve_flint(B)
            sage: C
            [ 12  40  28]
            [-12  -4  -4]
            [ -6 -25 -16]
            [ 12  34  16]

        ::

            sage: d
            12

        ::

            sage: A*C == d*B
            True

        Test wrong dimensions::

            sage: A = random_matrix(ZZ, 4, 4)
            sage: B = random_matrix(ZZ, 2, 3)
            sage: B._solve_flint(A)
            Traceback (most recent call last):
            ...
            ValueError: self must be a square matrix
            sage: A._solve_flint(B, right=False)
            Traceback (most recent call last):
            ...
            ArithmeticError: B's number of columns must match self's number of rows
            sage: A._solve_flint(B, right=True)
            Traceback (most recent call last):
            ...
            ArithmeticError: B's number of rows must match self's number of columns

        Check that this can be interrupted properly (:issue:`15453`)::

            sage: A = random_matrix(ZZ, 2000, 2000)
            sage: B = random_matrix(ZZ, 2000, 2000)
            sage: from sage.doctest.util import ensure_interruptible_after
            sage: with ensure_interruptible_after(2, max_wait_after_interrupt=8): A._solve_flint(B)

        AUTHORS:

        - Marc Masdeu (08/2014) following _solve_iml implementation of Martin Albrecht
        """
        cdef Matrix_integer_dense M
        cdef fmpz_t tmp
        cdef Integer den = Integer(0)
        if self._nrows != self._ncols:
            # This is *required* by the FLINT function we call below.
            raise ValueError("self must be a square matrix")

        if self.nrows() == 1:
            return B, self[0,0]

        if right:
            fmpz_init(tmp)
            if self._ncols != B._nrows:
                raise ArithmeticError("B's number of rows must match self's number of columns")

            n = self._ncols
            m = B._ncols

            if m == 0 or n == 0:
                return self.new_matrix(nrows = n, ncols = m), Integer(1)
            M = self._new(self._ncols,B._ncols)
            sig_on()
            fmpz_mat_solve(M._matrix,tmp,self._matrix,B._matrix)
            fmpz_get_mpz(den.value,tmp)
            fmpz_clear(tmp)
            if mpz_cmp_si(den.value,0) < 0:
                mpz_neg(den.value,den.value)
                fmpz_mat_neg(M._matrix,M._matrix)
            sig_off()
            return M,den

        else: # left
            if self._nrows != B._ncols:
                raise ArithmeticError("B's number of columns must match self's number of rows")

            n = self._nrows
            m = B._nrows

            if m == 0 or n == 0:
                return self.new_matrix(nrows = n, ncols = m), Integer(1)

            M,d = self.transpose()._solve_flint(B.transpose(), right=True)
            return M.transpose(),d

    def _rational_echelon_via_solve(self, solver='iml'):
        r"""
        Compute information that gives the reduced row echelon form (over
        QQ!) of a matrix with integer entries.

        INPUT:

        - ``self`` -- a matrix over the integers

        - ``solver`` -- either ``'iml'`` (default) or ``'flint'``

        OUTPUT:

        - ``pivots`` -- ordered list of integers that give the
          pivot column positions

        - ``nonpivots`` -- ordered list of the nonpivot column
          positions

        - ``X`` -- matrix with integer entries

        - ``d`` -- integer

        If you put standard basis vectors in order at the pivot columns,
        and put the matrix ``(1/d)*X`` everywhere else, then you get the
        reduced row echelon form of self, without zero rows at the bottom.

        .. NOTE::

           IML is the actual underlying `p`-adic solver that we
           use.

        AUTHORS:

        - William Stein

        ALGORITHM: I came up with this algorithm from scratch. As far as I
        know it is new. It's pretty simple, but it is ... (fast?!).

        Let A be the input matrix.


        #. Compute r = rank(A).

        #. Compute the pivot columns of the transpose `A^t` of
           `A`. One way to do this is to choose a random prime
           `p` and compute the row echelon form of `A^t`
           modulo `p` (if the number of pivot columns is not equal to
           `r`, pick another prime).

        #. Let `B` be the submatrix of `A` consisting of
           the rows corresponding to the pivot columns found in the previous
           step. Note that, aside from zero rows at the bottom, `B`
           and `A` have the same reduced row echelon form.

        #. Compute the pivot columns of `B`, again possibly by
           choosing a random prime `p` as in [2] and computing the
           Echelon form modulo `p`. If the number of pivot columns is
           not `r`, repeat with a different prime. Note - in this step
           it is possible that we mistakenly choose a bad prime `p`
           such that there are the right number of pivot columns modulo
           `p`, but they are at the wrong positions - e.g., imagine
           the augmented matrix `[pI|I]` - modulo `p` we would
           miss all the pivot columns. This is OK, since in the next step we
           would detect this, as the matrix we would obtain would not be in
           echelon form.

        #. Let `C` be the submatrix of `B` of pivot
           columns. Let `D` be the complementary submatrix of
           `B` of all non-pivot columns. Use a `p`-adic
           solver to find the matrix `X` and integer `d` such
           that `C (1/d) X=D`. I.e., solve a bunch of linear systems
           of the form `Cx = v`, where the columns of `X` are
           the solutions.

        #. Verify that we had chosen the correct pivot columns. Inspect the
           matrix `X` obtained in step 5. If when used to construct
           the echelon form of `B`, `X` indeed gives a matrix
           in reduced row echelon form, then we are done - output the pivot
           columns, `X`, and `d`. To tell if `X` is
           correct, do the following: For each column of `X`
           (corresponding to non-pivot column `i` of `B`),
           read up the column of `X` until finding the first nonzero
           position `j`; then verify that `j` is strictly less
           than the number of pivot columns of `B` that are strictly
           less than `i`. Otherwise, we got the pivots of `B`
           wrong - try again starting at step 4, but with a different random
           prime.
        """
        if self._nrows == 0:
            pivots = []
            nonpivots = list(range(self._ncols))
            X = self.__copy__()
            d = Integer(1)
            return pivots, nonpivots, X, d
        elif self._ncols == 0:
            pivots = []
            nonpivots = []
            X = self.__copy__()
            d = Integer(1)
            return pivots, nonpivots, X, d

        from sage.matrix.matrix_modn_dense_double import MAX_MODULUS
        A = self
        # Step 1: Compute the rank

        t = verbose('computing rank', level=2, caller_name='p-adic echelon')
        r = A.rank()
        verbose('done computing rank', level=2, t=t, caller_name='p-adic echelon')

        cdef randstate rstate = current_randstate()

        if r == self._nrows:
            # The input matrix already has full rank.
            B = A
        else:
            # Steps 2 and 3: Extract out a submatrix of full rank.
            i = 0
            while True:
                p = previous_prime(rstate.c_random() % (MAX_MODULUS-15000) + 10000)
                P = A._mod_int(p).transpose().pivots()
                if len(P) == r:
                    B = A.matrix_from_rows(P)
                    break
                else:
                    i += 1
                    if i == 50:
                        raise RuntimeError("Bug in _rational_echelon_via_solve in finding linearly independent rows.")

        _i = 0
        while True:
            _i += 1
            if _i == 50:
                raise RuntimeError("Bug in _rational_echelon_via_solve -- pivot columns keep being wrong.")

            # Step 4: Now we instead worry about computing the reduced row echelon form of B.
            i = 0
            while True:
                p = previous_prime(rstate.c_random() % (MAX_MODULUS-15000) + 10000)
                pivots = B._mod_int(p).pivots()
                if len(pivots) == r:
                    break
                else:
                    i += 1
                    if i == 50:
                        raise RuntimeError("Bug in _rational_echelon_via_solve in finding pivot columns.")

            # Step 5: Apply p-adic solver
            C = B.matrix_from_columns(pivots)
            pivots_ = set(pivots)
            non_pivots = [i for i in range(B.ncols()) if i not in pivots_]
            D = B.matrix_from_columns(non_pivots)
            t = verbose('calling %s solver'%solver, level=2, caller_name='p-adic echelon')
            if solver == 'iml':
                X, d = C._solve_iml(D, right=True)
            else:
                X, d = C._solve_flint(D, right=True)
            t = verbose('finished %s solver'%solver, level=2, caller_name='p-adic echelon', t=t)

            # Step 6: Verify that we had chosen the correct pivot columns.
            pivots_are_right = True
            for z in range(X.ncols()):
                if not pivots_are_right:
                    break
                i = non_pivots[z]
                np = len([k for k in pivots if k < i])
                for j in reversed(range(X.nrows())):
                    if X[j,z] != 0:
                        if j < np:
                            break # we're good -- go on to next column of X
                        else:
                            pivots_are_right = False
                            break
            if pivots_are_right:
                break

        # end while

        # Finally, return the answer.
        return pivots, non_pivots, X, d

    def decomposition(self, **kwds):
        """
        Return the decomposition of the free module on which this matrix A
        acts from the right (i.e., the action is x goes to x A), along with
        whether this matrix acts irreducibly on each factor. The factors
        are guaranteed to be sorted in the same way as the corresponding
        factors of the characteristic polynomial, and are saturated as ZZ
        modules.

        INPUT:

        - ``self`` -- a matrix over the integers

        - ``**kwds`` -- these are passed onto to the
          decomposition over QQ command

        EXAMPLES::

            sage: t = ModularSymbols(11,sign=1).hecke_matrix(2)
            sage: w = t.change_ring(ZZ)
            sage: w
            [ 3 -2]
            [ 0 -2]
            sage: w.charpoly().factor()
            (x - 3) * (x + 2)
            sage: w.decomposition()
            [
            (Free module of degree 2 and rank 1 over Integer Ring
            Echelon basis matrix:
            [ 5 -2], True),
            (Free module of degree 2 and rank 1 over Integer Ring
            Echelon basis matrix:
            [0 1], True)
            ]
        """
        F = self.charpoly().factor()
        if len(F) == 1:
            V = self.base_ring()**self.nrows()
            return decomp_seq([(V, F[0][1]==1)])

        A = self.change_ring(QQ)
        X = A.decomposition(**kwds)
        V = ZZ**self.nrows()
        if isinstance(X, tuple):
            D, E = X
            D = [(W.intersection(V), t) for W, t in D]
            E = [(W.intersection(V), t) for W, t in E]
            return decomp_seq(D), decomp_seq(E)
        else:
            return decomp_seq([(W.intersection(V), t) for W, t in X])

    def _add_row_and_maintain_echelon_form(self, row, pivots):
        r"""
        Assuming ``self`` is a full rank n x m matrix in reduced row Echelon
        form over `\ZZ` and row is a vector of degree m, this function creates
        a new matrix that is the echelon form of ``self`` with row appended to
        the bottom.

        .. WARNING::

           It is assumed that ``self`` is in echelon form.

        INPUT:

        - ``row`` -- a vector of degree m over `\ZZ`

        - ``pivots`` -- list of integers that are the pivot
          columns of ``self``

        OUTPUT: tuple of

        - ``matrix`` -- a matrix of in reduced row echelon form
          over `\ZZ`

        - ``pivots`` -- list of integers


        ALGORITHM: For each pivot column of self, we use the extended
        Euclidean algorithm to clear the column. The result is a new matrix
        B whose row span is the same as self.stack(row), and whose last row
        is 0 if and only if row is in the QQ-span of the rows of ``self``. If
        row is not in the QQ-span of the rows of self, then row is nonzero
        and suitable to be inserted into the top n rows of A to form a new
        matrix that is in reduced row echelon form. We then clear that
        corresponding new pivot column.

        EXAMPLES::

            sage: a = matrix(ZZ, 3, [1, 0, 110, 0, 3, 112, 0, 0, 221]); a
            [  1   0 110]
            [  0   3 112]
            [  0   0 221]
            sage: a._add_row_and_maintain_echelon_form(vector(ZZ,[1,2,3]),[0,1,2])
            (
            [1 0 0]
            [0 1 0]
            [0 0 1], [0, 1, 2]
            )
            sage: a._add_row_and_maintain_echelon_form(vector(ZZ,[0,0,0]),[0,1,2])
            (
            [  1   0 110]
            [  0   3 112]
            [  0   0 221], [0, 1, 2]
            )
            sage: a = matrix(ZZ, 2, [1, 0, 110, 0, 3, 112])
            sage: a._add_row_and_maintain_echelon_form(vector(ZZ,[1,2,3]),[0,1])
            (
            [  1   0 110]
            [  0   1 219]
            [  0   0 545], [0, 1, 2]
            )
        """
        cdef Py_ssize_t i, j, n = self._nrows, m = self._ncols

        from sage.matrix.constructor import matrix

        # 0. Base case
        if self.nrows() == 0:
            pos = row.nonzero_positions()
            if pos:
                pivots = [pos[0]]
                if row[pivots[0]] < 0:
                    row *= -1
            else:
                pivots = []
            return matrix([row]), pivots

        if row == 0:
            return self, pivots
        # 1. Create a new matrix that has row as the last row.
        row_mat = matrix(row)
        A = self.stack(row_mat)
        # 2. Working from the left, clear each column to put
        #    the resulting matrix back in echelon form.
        for i, p in enumerate(pivots):
            # p is the i-th pivot
            b = A[n,p]
            if not b:
                continue

            # (a). Take xgcd of pivot positions in last row and in ith
            # row.

            # TODO (optimize) -- change to use direct call to gmp and
            # no bounds checking!
            a = A[i,p]
            if not a:
                raise ZeroDivisionError("claimed pivot is not a pivot")
            if b % a == 0:
                # (b) Subtract a multiple of row i from row n.
                c = b // a
                if c:
                    for j in range(m):
                        A[n,j] -= c * A[i,j]
            else:
                # (b). More elaborate.
                #  Replace the ith row by s*A[i] + t*A[n], which will
                # have g in the i,p position, and replace the last row by
                # (b//g)*A[i] - (a//g)*A[n], which will have 0 in the i,p
                # position.
                g, s, t = a.xgcd(b)
                if not g:
                    raise ZeroDivisionError("claimed pivot is not a pivot (got a 0 gcd)")

                row_i = A.row(i)
                row_n = A.row(n)

                ag = a//g
                bg = b//g

                new_top = s*row_i  +  t*row_n
                new_bot = bg*row_i - ag*row_n

                # OK -- now we have to make sure the top part of the matrix
                # but with row i replaced by
                #     r = s*row_i[j]  +  t*row_n[j]
                # is put in rref.  We do this by recursively calling this
                # function with the top part of A (all but last row) and the
                # row r.

                zz = list(range(A.nrows() - 1))
                del zz[i]
                top_mat = A.matrix_from_rows(zz)
                new_pivots = list(pivots)
                del new_pivots[i]

                top_mat, pivots = top_mat._add_row_and_maintain_echelon_form(new_top, new_pivots)
                w = top_mat._add_row_and_maintain_echelon_form(new_bot, pivots)
                return w
        # 3. If it turns out that the last row is nonzero,
        #    insert last row in A sliding other rows down.
        v = A.row(n)
        new_pivots = list(pivots)
        if v != 0:
            i = v.nonzero_positions()[0]
            assert not (i in pivots), 'WARNING: bug in add_row -- i (=%s) should not be a pivot'%i

            # If pivot entry is negative negate this row.
            if v[i] < 0:
                v = -v

            # Determine where the last row should be inserted.
            new_pivots.append(i)
            new_pivots.sort()
            import bisect
            j = bisect.bisect(pivots, i)
            # The new row should go *before* row j, so it becomes row j
            A = A.insert_row(j, v)
        try:
            _clear_columns(A, new_pivots, A.nrows())
        except RuntimeError:
            raise ZeroDivisionError("mistake in claimed pivots")
        if A.row(A.nrows() - 1) == 0:
            A = A.matrix_from_rows(range(A.nrows()-1))
        return A, new_pivots

    #####################################################################################
    # Hermite form modulo D
    # This code below is by E. Burcin.  Thanks!
    #####################################################################################
    def _hnf_mod(self, D):
        """
        INPUT:

        - ``D`` -- small integer that is assumed to be a
          multiple of ``2*det(self)``

        OUTPUT: ``matrix`` -- the Hermite normal form of self

        EXAMPLES:

        A :exc:`ValueError` is raised if the matrix is not square,
        fixing :issue:`5548`::

            sage: random_matrix(ZZ,16,4)._hnf_mod(100)
            Traceback (most recent call last):
            ...
            ValueError: matrix is not square
        """
        t = verbose('hermite mod %s' % D, caller_name='matrix_integer_dense')
        if self._nrows != self._ncols:
            raise ValueError("matrix is not square")
        cdef Matrix_integer_dense res = self._new(self._nrows, self._ncols)
        self._hnf_modn(res, D)
        verbose('finished hnf mod', t, caller_name='matrix_integer_dense')
        return res

    cdef int _hnf_modn(Matrix_integer_dense self, Matrix_integer_dense res,
            unsigned int det) except -1:
        """
        Puts ``self`` into HNF form modulo det. Changes ``self`` in place.
        """
        cdef int* res_l
        cdef Py_ssize_t i,j,k
        res_l = self._hnf_modn_impl(det, self._nrows, self._ncols)
        k = 0
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                res.set_unsafe_int(i,j,res_l[k])
                k += 1
        sig_free(res_l)

    cdef int* _hnf_modn_impl(Matrix_integer_dense self, unsigned int det,
            Py_ssize_t nrows, Py_ssize_t ncols) except NULL:
        # NOTE: det should be at most 2^31-1, such that anything modulo
        # det fits in a 32-bit signed integer. To avoid overflow, we
        # need 64-bit arithmetic for some intermediate computations.
        cdef int* res
        cdef int* T_ent
        cdef int* *res_rows
        cdef int** T_rows
        cdef int* B
        cdef Py_ssize_t i, j, k
        cdef int R = det
        cdef int T_i_i, T_j_i, c1, c2, q, t
        cdef int u, v, d

        # allocate memory for result matrix
        res = <int*> sig_malloc(sizeof(int)*ncols*nrows)
        if res == NULL:
            raise MemoryError("out of memory allocating a matrix")
        res_rows = <int**> sig_malloc(sizeof(int*)*nrows)
        if res_rows == NULL:
            sig_free(res)
            raise MemoryError("out of memory allocating a matrix")

        # allocate memory for temporary matrix
        T_ent = <int*> sig_malloc(sizeof(int)*ncols*nrows)
        if T_ent == NULL:
            sig_free(res)
            sig_free(res_rows)
            raise MemoryError("out of memory allocating a matrix")
        T_rows = <int**> sig_malloc(sizeof(int*)*nrows)
        if T_rows == NULL:
            sig_free(res)
            sig_free(res_rows)
            sig_free(T_ent)
            raise MemoryError("out of memory allocating a matrix")

        # allocate memory for temporary row vector
        B = <int*>sig_malloc(sizeof(int)*nrows)
        if B == NULL:
            sig_free(res)
            sig_free(res_rows)
            sig_free(T_ent)
            sig_free(T_rows)
            raise MemoryError("out of memory allocating a matrix")

        # initialize the row pointers
        k = 0
        for i in range(nrows):
            res_rows[i] = res + k
            T_rows[i] = T_ent + k
            k += nrows

        cdef mpz_t tmp
        mpz_init(tmp)
        # copy entries from self to temporary matrix
        k = 0
        for i in range(nrows):
            for j in range(ncols):
                self.get_unsafe_mpz(i,j,tmp)
                mpz_mod_ui(tmp, tmp, det)
                T_ent[k] = mpz_get_si(tmp)
                k += 1
        mpz_clear(tmp)

        # initialize variables
        i = 0
        j = 0

        while True:
            if j == nrows-1:
                T_i_i = T_rows[i][i]
                d = ai.c_xgcd_int(T_i_i, R, &u, &v)
                for k from 0 <= k < i:
                    res_rows[i][k] = 0
                for k from i <= k < ncols:
                    t = ((<int64_t>u)*T_rows[i][k])%R
                    if t < 0:
                        t += R
                    res_rows[i][k] = t
                if res_rows[i][i] == 0:
                    res_rows[i][i] = R
                d = res_rows[i][i]
                for j from 0 <= j < i:
                    q = res_rows[j][i]/d
                    for k from i <= k < ncols:
                        u = (res_rows[j][k] - (<int64_t>q)*res_rows[i][k])%R
                        if u < 0:
                            u += R
                        res_rows[j][k] = u

                R = R/d
                i += 1
                j = i
                if i == nrows:
                    break # return res
                if T_rows[i][i] == 0:
                    T_rows[i][i] = R
                continue

            j += 1
            if T_rows[j][i] == 0:
                continue

            T_i_i = T_rows[i][i]
            T_j_i = T_rows[j][i]
            d = ai.c_xgcd_int(T_i_i, T_j_i, &u, &v)
            if d != T_i_i:
                for k from i <= k < ncols:
                    B[k] = ((<int64_t>u)*T_rows[i][k] + (<int64_t>v)*T_rows[j][k])%R
            c1 = T_i_i/d
            c2 = -T_j_i/d
            for k from i <= k < ncols:
                T_rows[j][k] = ((<int64_t>c1)*T_rows[j][k] + (<int64_t>c2)*T_rows[i][k])%R
            if d != T_i_i:
                for k from i <= k < ncols:
                    T_rows[i][k] = B[k]

        sig_free(B)
        sig_free(res_rows)
        sig_free(T_ent)
        sig_free(T_rows)
        return res

    #################################################################
    # operations with matrices
    #################################################################

    def row(self, Py_ssize_t i, from_list=False):
        """
        Return the `i`-th row of this matrix as a dense vector.

        INPUT:

        - ``i`` -- integer

        - ``from_list`` -- ignored

        EXAMPLES::

            sage: m = matrix(ZZ, 2, [1, -2, 3, 4])
            sage: m.row(0)
            (1, -2)
            sage: m.row(1)
            (3, 4)
            sage: m.row(1, from_list=True)
            (3, 4)
            sage: m.row(-2)
            (1, -2)

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
        cdef Vector_integer_dense v = parent.zero_vector()
        for j in range(self._ncols):
            fmpz_get_mpz(v._entries[j], fmpz_mat_entry(self._matrix, i, j))
        return v

    def column(self, Py_ssize_t i, from_list=False):
        """
        Return the `i`-th column of this matrix as a dense vector.

        INPUT:

        - ``i`` -- integer

        - ``from_list`` -- ignored

        EXAMPLES::

            sage: m = matrix(ZZ, 3, 2, [1, -2, 3, 4, -1, 0])
            sage: m.column(1)
            (-2, 4, 0)
            sage: m.column(1, from_list=True)
            (-2, 4, 0)
            sage: m.column(-1)
            (-2, 4, 0)
            sage: m.column(-2)
            (1, 3, -1)

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
        cdef Vector_integer_dense v = parent.zero_vector()
        for j in range(self._nrows):
            fmpz_get_mpz(v._entries[j], fmpz_mat_entry(self._matrix, j, i))
        return v

    cdef _stack_impl(self, bottom):
        r"""
        Return the matrix ``self`` on top of ``bottom``::

            [  self  ]
            [ bottom ]

        EXAMPLES::

            sage: M = Matrix(ZZ, 2, 3, range(6))
            sage: N = Matrix(ZZ, 1, 3, [10,11,12])
            sage: M.stack(N)
            [ 0  1  2]
            [ 3  4  5]
            [10 11 12]

        A vector may be stacked below a matrix. ::

            sage: A = matrix(ZZ, 2, 4, range(8))
            sage: v = vector(ZZ, 4, range(4))
            sage: A.stack(v)
            [0 1 2 3]
            [4 5 6 7]
            [0 1 2 3]

        The ``subdivide`` option will add a natural subdivision between
        ``self`` and ``bottom``.  For more details about how subdivisions
        are managed when stacking, see
        :meth:`sage.matrix.matrix1.Matrix.stack`.  ::

            sage: A = matrix(ZZ, 3, 4, range(12))
            sage: B = matrix(ZZ, 2, 4, range(8))
            sage: A.stack(B, subdivide=True)
            [ 0  1  2  3]
            [ 4  5  6  7]
            [ 8  9 10 11]
            [-----------]
            [ 0  1  2  3]
            [ 4  5  6  7]
        """
        cdef Matrix_integer_dense other = <Matrix_integer_dense>bottom
        cdef Matrix_integer_dense Z
        Z = self.new_matrix(nrows=self._nrows + other._nrows, ncols=self._ncols)

        cdef Py_ssize_t r, c
        cdef Py_ssize_t nr = self._nrows
        for r in range(self._nrows):
            for c in range(self._ncols):
                fmpz_set(fmpz_mat_entry(Z._matrix, r, c),fmpz_mat_entry(self._matrix, r, c))
        for r in range(other._nrows):
            for c in range(other._ncols):
                fmpz_set(fmpz_mat_entry(Z._matrix, r+nr, c),fmpz_mat_entry(other._matrix, r, c))

        return Z

    def augment(self, right, subdivide=False):
        r"""
        Return a new matrix formed by appending the matrix
        (or vector) ``right`` on the right side of ``self``.

        INPUT:

        - ``right`` -- a matrix, vector or free module element, whose
          dimensions are compatible with ``self``

        - ``subdivide`` -- (default: ``False``) request the resulting
          matrix to have a new subdivision, separating ``self`` from ``right``

        OUTPUT:

        A new matrix formed by appending ``right`` onto the right side of ``self``.
        If ``right`` is a vector (or free module element) then in this context
        it is appropriate to consider it as a column vector.  (The code first
        converts a vector to a 1-column matrix.)

        EXAMPLES::

            sage: A = matrix(ZZ, 4, 5, range(20))
            sage: B = matrix(ZZ, 4, 3, range(12))
            sage: A.augment(B)
            [ 0  1  2  3  4  0  1  2]
            [ 5  6  7  8  9  3  4  5]
            [10 11 12 13 14  6  7  8]
            [15 16 17 18 19  9 10 11]

        A vector may be augmented to a matrix. ::

            sage: A = matrix(ZZ, 3, 5, range(15))
            sage: v = vector(ZZ, 3, range(3))
            sage: A.augment(v)
            [ 0  1  2  3  4  0]
            [ 5  6  7  8  9  1]
            [10 11 12 13 14  2]

        The ``subdivide`` option will add a natural subdivision between
        ``self`` and ``right``.  For more details about how subdivisions
        are managed when augmenting, see
        :meth:`sage.matrix.matrix1.Matrix.augment`.  ::

            sage: A = matrix(ZZ, 3, 5, range(15))
            sage: B = matrix(ZZ, 3, 3, range(9))
            sage: A.augment(B, subdivide=True)
            [ 0  1  2  3  4| 0  1  2]
            [ 5  6  7  8  9| 3  4  5]
            [10 11 12 13 14| 6  7  8]

        Errors are raised if the sizes are incompatible. ::

            sage: A = matrix(ZZ, [[1, 2],[3, 4]])
            sage: B = matrix(ZZ, [[10, 20], [30, 40], [50, 60]])
            sage: A.augment(B)
            Traceback (most recent call last):
            ...
            TypeError: number of rows must be the same, not 2 != 3
        """
        if not isinstance(right, Matrix) and hasattr(right, '_vector_'):
            right = right.column()
        if self._nrows != right.nrows():
            raise TypeError('number of rows must be the same, not {0} != {1}'.format(self._nrows, right.nrows()))
        if not (self._base_ring is right.base_ring()):
            right = right.change_ring(self._base_ring)

        cdef Matrix_integer_dense other = right.dense_matrix()
        m = self._nrows
        ns, na = self._ncols, other._ncols
        n = ns + na

        cdef Matrix_integer_dense Z
        Z = self.new_matrix(nrows=m, ncols=n)
        cdef Py_ssize_t i, j
        for i from 0 <= i < m:
            for j from 0 <= j < ns:
                fmpz_set(fmpz_mat_entry(Z._matrix,i,j),
                         fmpz_mat_entry(self._matrix,i,j))
            for j from 0 <= j < na:
                fmpz_set(fmpz_mat_entry(Z._matrix,i,j + ns),
                         fmpz_mat_entry(other._matrix,i,j))
        if subdivide:
            Z._subdivide_on_augment(self, other)
        return Z

    def insert_row(self, Py_ssize_t index, row):
        """
        Create a new matrix from ``self`` with.

        INPUT:

        - ``index`` -- integer

        - ``row`` -- a vector

        EXAMPLES::

            sage: X = matrix(ZZ,3,range(9)); X
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: X.insert_row(1, [1,5,-10])
            [  0   1   2]
            [  1   5 -10]
            [  3   4   5]
            [  6   7   8]
            sage: X.insert_row(0, [1,5,-10])
            [  1   5 -10]
            [  0   1   2]
            [  3   4   5]
            [  6   7   8]
            sage: X.insert_row(3, [1,5,-10])
            [  0   1   2]
            [  3   4   5]
            [  6   7   8]
            [  1   5 -10]
        """
        cdef Matrix_integer_dense res = self._new(self._nrows + 1, self._ncols)
        cdef Py_ssize_t j
        cdef Integer z
        cdef fmpz_t zflint
        if index < 0:
            raise ValueError("index must be nonnegative")
        if index > self._nrows:
            raise ValueError("index must be less than number of rows")
        fmpz_init(zflint)

        for j from 0 <= j < self._ncols:
            for i from 0 <= i < index:
                fmpz_init_set(fmpz_mat_entry(res._matrix,i,j), fmpz_mat_entry(self._matrix,i,j))

            z = row[j]
            fmpz_set_mpz(zflint,z.value)
            fmpz_init_set(fmpz_mat_entry(res._matrix,index,j), zflint)

            for i from index <= i < self._nrows:
                fmpz_init_set(fmpz_mat_entry(res._matrix,i+1,j), fmpz_mat_entry(self._matrix,i,j))

        fmpz_clear(zflint)
        return res

    def _delete_zero_columns(self):
        """
        Return matrix obtained from ``self`` by deleting all zero columns along
        with the positions of those columns.

        OUTPUT: (matrix, list of integers)

        EXAMPLES::

            sage: a = matrix(ZZ, 2,3, [1,0,3,-1,0,5]); a
            [ 1  0  3]
            [-1  0  5]
            sage: a._delete_zero_columns()
            (
            [ 1  3]
            [-1  5], [1]
            )
        """
        zero_cols = [i for i, v in enumerate(self.columns()) if v.is_zero()]
        s = set(zero_cols)
        nonzero_cols = [i for i in range(self.ncols()) if i not in s]
        return self.matrix_from_columns(nonzero_cols), zero_cols

    def _insert_zero_columns(self, cols):
        """
        Return matrix obtained by ``self`` by inserting zero columns so that
        the columns with positions specified in cols are all 0.

        INPUT:

        - ``cols`` -- list of nonnegative integers

        OUTPUT: matrix

        EXAMPLES::

            sage: a = matrix(ZZ, 2,3, [1,0,3,-1,0,5]); a
            [ 1  0  3]
            [-1  0  5]
            sage: b, cols = a._delete_zero_columns()
            sage: b
            [ 1  3]
            [-1  5]
            sage: cols
            [1]
            sage: b._insert_zero_columns(cols)
            [ 1  0  3]
            [-1  0  5]
        """
        if not cols:
            return self
        cdef Py_ssize_t i, c, r, nc = max(self._ncols + len(cols), max(cols)+1)
        cdef Matrix_integer_dense A = self.new_matrix(self._nrows, nc)
        # Now fill in the entries of A that come from self.
        cols_set = set(cols)
        cols_ins = [j for j in range(nc) if j not in cols_set]
        for r from 0 <= r < self._nrows:
            i = 0
            for c in cols_ins:
                # The following does this quickly: A[r,c] = self[r,i]
                fmpz_set(fmpz_mat_entry(A._matrix,r,c),fmpz_mat_entry(self._matrix,r,i))
                i += 1
        return A

    def _factor_out_common_factors_from_each_row(self):
        """
        Very very quickly modify ``self`` so that the gcd of the entries in
        each row is 1 by dividing each row by the common gcd.

        EXAMPLES::

            sage: a = matrix(ZZ, 3, [-9, 3, -3, -36, 18, -5, -40, -5, -5, -20, -45, 15, 30, -15, 180])
            sage: a
            [ -9   3  -3 -36  18]
            [ -5 -40  -5  -5 -20]
            [-45  15  30 -15 180]
            sage: a._factor_out_common_factors_from_each_row()
            sage: a
            [ -3   1  -1 -12   6]
            [ -1  -8  -1  -1  -4]
            [ -3   1   2  -1  12]
        """
        self.check_mutability()

        cdef fmpz_t g
        cdef fmpz_t tmp
        fmpz_init(g)
        fmpz_init(tmp)
        cdef long i, j

        for i from 0 <= i < self._nrows:
            fmpz_set_ui(g, 0)
            for j from 0 <= j < self._ncols:
                fmpz_gcd(g, g, fmpz_mat_entry(self._matrix,i,j))
                if fmpz_cmp_ui(g, 1) == 0:
                    break
            if fmpz_cmp_ui(g, 1) != 0:
                # divide through row
                for j from 0 <= j < self._ncols:
                    fmpz_set(tmp,fmpz_mat_entry(self._matrix,i,j))
                    fmpz_divexact(tmp, tmp, g)
                    fmpz_set(fmpz_mat_entry(self._matrix,i,j),tmp)
        fmpz_clear(g)
        fmpz_clear(tmp)

    def gcd(self):
        """
        Return the gcd of all entries of self; very fast.

        EXAMPLES::

            sage: a = matrix(ZZ,2, [6,15,-6,150])
            sage: a.gcd()
            3
        """
        cdef Integer g = Integer(0)
        cdef Py_ssize_t i, j
        cdef fmpz_t tmpgcd
        fmpz_init(tmpgcd)
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                fmpz_gcd(tmpgcd,tmpgcd,fmpz_mat_entry(self._matrix,i,j))
                if fmpz_cmp_ui(tmpgcd, 1) == 0:
                    fmpz_get_mpz(g.value,tmpgcd)
                    return g
        fmpz_get_mpz(g.value,tmpgcd)
        fmpz_clear(tmpgcd)
        return g

    def _change_ring(self, ring):
        """
        Return the matrix obtained by coercing the entries of this matrix
        into the given ring.

        EXAMPLES::

            sage: a = matrix(ZZ,2,[1,-7,3,5])
            sage: a._change_ring(RDF)
            [ 1.0 -7.0]
            [ 3.0  5.0]
        """
        if ring == RDF:
            from sage.matrix.change_ring import integer_to_real_double_dense
            return integer_to_real_double_dense(self)
        else:
            raise NotImplementedError

    def _singular_(self, singular=None):
        r"""
        Return Singular representation of this integer matrix.

        INPUT:

        - ``singular`` -- Singular interface instance (default:
          ``None``)

        EXAMPLES::

            sage: A = matrix(ZZ, 3, 3, [-8, 2, 0, 0, 1, -1, 2, 1, -95])
            sage: As = singular(A); As
            -8     2     0
            0     1    -1
            2     1   -95
            sage: As.type()
            'intmat'
        """
        if singular is None:
            from sage.interfaces.singular import singular as singular_default
            singular = singular_default

        name = singular._next_var_name()
        values = str(self.list())[1:-1]
        singular.eval("intmat %s[%d][%d] = %s"%(name, self.nrows(), self.ncols(), values))

        from sage.interfaces.singular import SingularElement
        return SingularElement(singular, 'foobar', name, True)

    def transpose(self):
        """
        Return the transpose of self, without changing ``self``.

        EXAMPLES:

        We create a matrix, compute its transpose, and note that the
        original matrix is not changed.

        ::

            sage: A = matrix(ZZ, 2, 3, range(6))
            sage: type(A)
            <class 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'>
            sage: B = A.transpose()
            sage: print(B)
            [0 3]
            [1 4]
            [2 5]
            sage: print(A)
            [0 1 2]
            [3 4 5]

        ``.T`` is a convenient shortcut for the transpose::

            sage: A.T
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
        cdef Matrix_integer_dense A
        A = self._new(self._ncols,self._nrows)
        sig_on()
        fmpz_mat_transpose(A._matrix,self._matrix)
        sig_off()
        if self._subdivisions is not None:
            row_divs, col_divs = self.subdivisions()
            A.subdivide(col_divs, row_divs)
        return A

    def antitranspose(self):
        """
        Return the antitranspose of ``self``, without changing ``self``.

        EXAMPLES::

            sage: A = matrix(2,3,range(6))
            sage: type(A)
            <class 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'>
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
        nr, nc = self._nrows, self._ncols

        cdef Matrix_integer_dense A
        A = self._new(nc,nr)
        cdef Py_ssize_t i, j
        cdef Py_ssize_t ri, rj  # reversed i and j
        sig_on()
        ri = nr
        for i from 0 <= i < nr:
            rj = nc
            ri =  ri-1
            for j from 0 <= j < nc:
                rj = rj-1
                fmpz_init_set(fmpz_mat_entry(A._matrix, rj, ri),
                              fmpz_mat_entry(self._matrix, i, j))
        sig_off()

        if self._subdivisions is not None:
            row_divs, col_divs = self.subdivisions()
            A.subdivide([nc - t for t in reversed(col_divs)],
                        [nr - t for t in reversed(row_divs)])
        return A

    def __pari__(self):
        """
        Return PARI C-library version of this matrix.

        EXAMPLES::

            sage: a = matrix(ZZ,2,2,[1,2,3,4])
            sage: a.__pari__()
            [1, 2; 3, 4]
            sage: pari(a)
            [1, 2; 3, 4]
            sage: type(pari(a))
            <class 'cypari2.gen.Gen'>
        """
        return integer_matrix(self._matrix, 0)

    def _rank_pari(self):
        """
        Rank of this matrix, computed using PARI.  The computation is
        done entirely on the PARI stack, then the PARI stack is
        cleared.  This function is most useful for very small
        matrices.

        EXAMPLES::

            sage: matrix(ZZ,3,[1..9])._rank_pari()
            2
        """
        sig_on()
        cdef long r = rank(pari_GEN(self))
        clear_stack()
        return r

    def _hnf_pari(self, int flag=0, bint include_zero_rows=True):
        """
        Hermite normal form of this matrix, computed using PARI.

        INPUT:

        - ``flag`` -- 0 (default), 1, 3 or 4 (see docstring for
          ``pari.mathnf``)

        - ``include_zero_rows`` -- boolean; if ``False``, do not include
          any of the zero rows at the bottom of the matrix in the output

        .. NOTE::

            In no cases is the transformation matrix returned by this
            function.

        EXAMPLES::

            sage: matrix(ZZ,3,[1..9])._hnf_pari()
            [1 2 3]
            [0 3 6]
            [0 0 0]
            sage: matrix(ZZ,3,[1..9])._hnf_pari(1)
            [1 2 3]
            [0 3 6]
            [0 0 0]
            sage: matrix(ZZ,3,[1..9])._hnf_pari(3)
            [1 2 3]
            [0 3 6]
            [0 0 0]
            sage: matrix(ZZ,3,[1..9])._hnf_pari(4)
            [1 2 3]
            [0 3 6]
            [0 0 0]

        Check that ``include_zero_rows=False`` works correctly::

            sage: matrix(ZZ,3,[1..9])._hnf_pari(0, include_zero_rows=False)
            [1 2 3]
            [0 3 6]
            sage: matrix(ZZ,3,[1..9])._hnf_pari(1, include_zero_rows=False)
            [1 2 3]
            [0 3 6]
            sage: matrix(ZZ,3,[1..9])._hnf_pari(3, include_zero_rows=False)
            [1 2 3]
            [0 3 6]
            sage: matrix(ZZ,3,[1..9])._hnf_pari(4, include_zero_rows=False)
            [1 2 3]
            [0 3 6]

        Check that :issue:`12346` is fixed::

            sage: pari('mathnf(Mat([0,1]), 4)')
            [Mat(1), [1, 0; 0, 1]]
        """
        sig_on()
        A = _new_GEN_from_fmpz_mat_t_rotate90(self._matrix)
        H = mathnf0(A, flag)
        if typ(H) == t_VEC:
            H = gel(H, 1)
        GenH = new_gen(H)
        return extract_hnf_from_pari_matrix(self, GenH, include_zero_rows)

    def p_minimal_polynomials(self, p, s_max=None):
        r"""
        Compute `(p^s)`-minimal polynomials `\nu_s` of this matrix.

        For `s  \ge 0`, a `(p^s)`-minimal polynomial of
        a matrix `B` is a monic polynomial `f \in \ZZ[X]` of
        minimal degree such that all entries of `f(B)` are divisible
        by `p^s`.

        Compute a finite subset `\mathcal{S}` of the positive
        integers and `(p^s)`-minimal polynomials
        `\nu_s` for `s \in \mathcal{S}`.

        For `0 < t \le \max \mathcal{S}`, a `(p^t)`-minimal polynomial is
        given by `\nu_s` where
        `s = \min\{ r \in \mathcal{S} \mid r\ge t \}`.
        For `t > \max\mathcal{S}`, the minimal polynomial of `B` is
        also a `(p^t)`-minimal polynomial.

        INPUT:

        - ``p`` -- a prime in `\ZZ`

        - ``s_max`` -- positive integer (default: ``None``); if set, only
          `(p^s)`-minimal polynomials for ``s <= s_max`` are computed
          (see below for details)

        OUTPUT:

        A dictionary. Keys are the finite set `\mathcal{S}`, the
        values are the associated `(p^s)`-minimal polynomials `\nu_s`,
        `s\in\mathcal{S}`.

        Setting ``s_max`` only affects the output if ``s_max`` is at
        most `\max\mathcal{S}` where `\mathcal{S}` denotes the full
        set. In that case, only those `\nu_s` with ``s <= s_max`` are
        returned where ``s_max`` is always included even if it is not
        included in the full set `\mathcal{S}`.

        EXAMPLES::

            sage: B = matrix(ZZ, [[1, 0, 1], [1, -2, -1], [10, 0, 0]])
            sage: B.p_minimal_polynomials(2)
            {2: x^2 + 3*x + 2}

        .. SEEALSO::

            :mod:`~sage.matrix.compute_J_ideal`,
            :meth:`~sage.matrix.compute_J_ideal.ComputeMinimalPolynomials.p_minimal_polynomials`
        """
        from sage.matrix.compute_J_ideal import ComputeMinimalPolynomials
        return ComputeMinimalPolynomials(self).p_minimal_polynomials(p, s_max)

    def null_ideal(self, b=0):
        r"""
        Return the `(b)`-ideal of this matrix.

        Let `B` be a `n \times n` matrix. The *null ideal* modulo `b`,
        or `(b)`-ideal, is

        .. MATH::

            N_{(b)}(B) = \{f \in \ZZ[X] \mid f(B) \in M_n(b\ZZ)\}.

        INPUT:

        - ``b`` -- an element of `\ZZ` (default: 0)

        OUTPUT: an ideal in `\ZZ[X]`

        EXAMPLES::

            sage: B = matrix(ZZ, [[1, 0, 1], [1, -2, -1], [10, 0, 0]])
            sage: B.null_ideal()
            Principal ideal (x^3 + x^2 - 12*x - 20) of
                Univariate Polynomial Ring in x over Integer Ring
            sage: B.null_ideal(8)
            Ideal (8, x^3 + x^2 - 12*x - 20, 2*x^2 + 6*x + 4) of
                Univariate Polynomial Ring in x over Integer Ring
            sage: B.null_ideal(6)
            Ideal (6, 2*x^3 + 2*x^2 - 24*x - 40, 3*x^2 + 3*x) of
                Univariate Polynomial Ring in x over Integer Ring

        .. SEEALSO::

            :mod:`~sage.matrix.compute_J_ideal`,
            :meth:`~sage.matrix.compute_J_ideal.ComputeMinimalPolynomials.null_ideal`
        """
        from sage.matrix.compute_J_ideal import ComputeMinimalPolynomials
        return ComputeMinimalPolynomials(self).null_ideal(b)

    def integer_valued_polynomials_generators(self):
        r"""
        Determine the generators of the ring of integer valued polynomials on this
        matrix.

        OUTPUT:

        A pair ``(mu_B, P)`` where ``P`` is a list of polynomials in `\QQ[X]`
        such that

        .. MATH::

           \{f \in \QQ[X] \mid f(B) \in M_n(\ZZ)\}
               = \mu_B \QQ[X] + \sum_{g\in P} g \ZZ[X]

        where `B` is this matrix.

        EXAMPLES::

            sage: B = matrix(ZZ, [[1, 0, 1], [1, -2, -1], [10, 0, 0]])
            sage: B.integer_valued_polynomials_generators()
            (x^3 + x^2 - 12*x - 20, [1, 1/4*x^2 + 3/4*x + 1/2])

        .. SEEALSO::

            :mod:`~sage.matrix.compute_J_ideal`,
            :meth:`~sage.matrix.compute_J_ideal.ComputeMinimalPolynomials.integer_valued_polynomials_generators`
        """
        from sage.matrix.compute_J_ideal import ComputeMinimalPolynomials
        return ComputeMinimalPolynomials(self).integer_valued_polynomials_generators()


cdef inline GEN pari_GEN(Matrix_integer_dense B) noexcept:
    r"""
    Create the PARI GEN object on the stack defined by the integer
    matrix B. This is used internally by the function for conversion
    of matrices to PARI.

    For internal use only; this directly uses the PARI stack.
    One should call ``sig_on()`` before and ``sig_off()`` after.
    """
    cdef GEN A = _new_GEN_from_fmpz_mat_t(B._matrix)
    return A


cdef extract_hnf_from_pari_matrix(Matrix_integer_dense self, Gen H, bint include_zero_rows):
    cdef mpz_t tmp
    mpz_init(tmp)

    # Figure out how many columns we got back.
    cdef long H_nc = glength(H.g)  # number of columns
    # Now get the resulting Hermite form matrix back to Sage, suitably re-arranged.
    cdef Matrix_integer_dense B
    if include_zero_rows:
        B = self.new_matrix()
    else:
        B = self.new_matrix(nrows=H_nc)
    cdef long i, j
    for i in range(self._ncols):
        for j in range(H_nc):
            sig_check()
            INT_to_mpz(tmp, gcoeff(H.g, i+1, H_nc-j))
            fmpz_set_mpz(fmpz_mat_entry(B._matrix,j,self._ncols-i-1),tmp)
    mpz_clear(tmp)
    return B


cdef _clear_columns(Matrix_integer_dense A, pivots, Py_ssize_t n):
    # Clear all columns
    cdef Py_ssize_t i, k, p, l, m = A._ncols
    cdef fmpz_t c,t
    sig_on()
    fmpz_init(c)
    fmpz_init(t)
    for i from 0 <= i < len(pivots):
        p = pivots[i]
        for k from 0 <= k < n:
            if k != i:
                if fmpz_cmp_si(fmpz_mat_entry(A._matrix,k,p),0):
                    fmpz_fdiv_q(c, fmpz_mat_entry(A._matrix,k,p), fmpz_mat_entry(A._matrix,i,p))
                    # subtract off c*v from row k; resulting A[k,i] entry will be < b, hence in Echelon form.
                    for l from 0 <= l < m:
                        fmpz_mul(t, c, fmpz_mat_entry(A._matrix,i,l))
                        fmpz_sub(fmpz_mat_entry(A._matrix,k,l), fmpz_mat_entry(A._matrix,k,l), t)

    fmpz_clear(c)
    fmpz_clear(t)
    sig_off()


cpdef _lift_crt(Matrix_integer_dense M, residues, moduli=None):
    """
    INPUT:

    - ``M`` -- a ``Matrix_integer_dense``; will be modified to hold
      the output

    - ``residues`` -- list of ``Matrix_modn_dense_template``; the
      matrix to reconstruct modulo primes

    OUTPUT: the matrix whose reductions modulo primes are the input ``residues``

    TESTS::

        sage: from sage.matrix.matrix_integer_dense import _lift_crt
        sage: T1 = Matrix(Zmod(5), 4, 4, [1, 4, 4, 0, 2, 0, 1, 4, 2, 0, 4, 1, 1, 4, 0, 3])
        sage: T2 = Matrix(Zmod(7), 4, 4, [1, 4, 6, 0, 2, 0, 1, 2, 4, 0, 6, 6, 1, 6, 0, 5])
        sage: T3 = Matrix(Zmod(11), 4, 4, [1, 4, 10, 0, 2, 0, 1, 9, 8, 0, 10, 6, 1, 10, 0, 9])
        sage: _lift_crt(Matrix(ZZ, 4, 4), [T1, T2, T3])
        [ 1  4 -1  0]
        [ 2  0  1  9]
        [-3  0 -1  6]
        [ 1 -1  0 -2]

        sage: from sage.arith.multi_modular import MultiModularBasis
        sage: mm = MultiModularBasis([5,7,11])
        sage: _lift_crt(Matrix(ZZ, 4, 4), [T1, T2, T3], mm)
        [ 1  4 -1  0]
        [ 2  0  1  9]
        [-3  0 -1  6]
        [ 1 -1  0 -2]

    The modulus must be smaller than the maximum for the multi-modular
    reconstruction (using ``mod_int``) and also smaller than the limit
    for ``Matrix_modn_dense_double`` to be able to represent the
    ``residues`` ::

        sage: from sage.arith.multi_modular import MAX_MODULUS as MAX_multi_modular
        sage: from sage.matrix.matrix_modn_dense_double import MAX_MODULUS as MAX_modn_dense_double
        sage: MAX_MODULUS = min(MAX_multi_modular, MAX_modn_dense_double)
        sage: p0 = previous_prime(MAX_MODULUS)
        sage: p1 = previous_prime(p0)
        sage: mmod = [matrix(GF(p0), [[-1, 0, 1, 0, 0, 1, 1, 0, 0, 0, p0-1, p0-2]]),
        ....:         matrix(GF(p1), [[-1, 0, 1, 0, 0, 1, 1, 0, 0, 0, p1-1, p1-2]])]
        sage: _lift_crt(Matrix(ZZ, 1, 12), mmod)
        [-1  0  1  0  0  1  1  0  0  0 -1 -2]
    """

    cdef size_t i, j, k
    cdef Py_ssize_t nr, n
    cdef mpz_t *tmp = <mpz_t *>sig_malloc(sizeof(mpz_t) * M._ncols)
    n = len(residues)
    if n == 0:   # special case: obviously residues[0] wouldn't make sense here.
        return M
    nr = residues[0].nrows()
    nc = residues[0].ncols()

    if moduli is None:
        moduli = MultiModularBasis([m.base_ring().order() for m in residues])
    else:
        if len(residues) != len(moduli):
            raise IndexError("Number of residues (%s) does not match number of moduli (%s)"%(len(residues), len(moduli)))

    cdef MultiModularBasis mm
    mm = moduli

    for b in residues:
        if not (isinstance(b, Matrix_modn_dense_float) or
                isinstance(b, Matrix_modn_dense_double)):
            raise TypeError("Can only perform CRT on list of matrices mod n.")

    cdef mod_int **row_list
    row_list = <mod_int**>sig_malloc(sizeof(mod_int*) * n)
    if row_list == NULL:
        raise MemoryError("out of memory allocating multi-modular coefficient list")

    sig_on()
    for k in range(n):
        row_list[k] = <mod_int *>sig_malloc(sizeof(mod_int) * nc)
        if row_list[k] == NULL:
            raise MemoryError("out of memory allocating multi-modular coefficient list")

    for j in range(M._ncols):
        mpz_init(tmp[j])

    for i in range(nr):
        for k in range(n):
            (<Matrix_modn_dense_template>residues[k])._copy_row_to_mod_int_array(row_list[k],i)
        mm.mpz_crt_vec(tmp, row_list, nc)
        for j in range(nc):
            M.set_unsafe_mpz(i,j,tmp[j])

    for k in range(n):
        sig_free(row_list[k])
    for j in range(M._ncols):
        mpz_clear(tmp[j])
    sig_free(row_list)
    sig_free(tmp)
    sig_off()
    return M
