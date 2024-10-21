# sage_setup: distribution = sagemath-linbox
r"""
Sparse matrices over `\ZZ/n\ZZ` for `n` small

This is a compiled implementation of sparse matrices over
`\ZZ/n\ZZ` for `n` small.

.. TODO::

    move vectors into a Cython vector class - add _add_ and _mul_ methods.

EXAMPLES::

    sage: a = matrix(Integers(37),3,3,range(9),sparse=True); a
    [0 1 2]
    [3 4 5]
    [6 7 8]
    sage: type(a)
    <class 'sage.matrix.matrix_modn_sparse.Matrix_modn_sparse'>
    sage: parent(a)
    Full MatrixSpace of 3 by 3 sparse matrices over Ring of integers modulo 37
    sage: a^2
    [15 18 21]
    [ 5 17 29]
    [32 16  0]
    sage: a+a
    [ 0  2  4]
    [ 6  8 10]
    [12 14 16]
    sage: b = a.new_matrix(2,3,range(6)); b
    [0 1 2]
    [3 4 5]
    sage: a*b
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for *: 'Full MatrixSpace of 3 by 3 sparse matrices over Ring of integers modulo 37' and 'Full MatrixSpace of 2 by 3 sparse matrices over Ring of integers modulo 37'
    sage: b*a
    [15 18 21]
    [ 5 17 29]

::

    sage: TestSuite(a).run()
    sage: TestSuite(b).run()

::

    sage: a.echelonize(); a
    [ 1  0 36]
    [ 0  1  2]
    [ 0  0  0]
    sage: b.echelonize(); b
    [ 1  0 36]
    [ 0  1  2]
    sage: a.pivots()
    (0, 1)
    sage: b.pivots()
    (0, 1)
    sage: a.rank()
    2
    sage: b.rank()
    2
    sage: a[2,2] = 5
    sage: a.rank()
    3

TESTS::

    sage: matrix(Integers(37),0,0,sparse=True).inverse()
    []
"""

# ****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from libc.stdint cimport uint64_t
from libc.limits cimport UINT_MAX

from cysignals.memory cimport check_calloc, sig_free
from cysignals.signals cimport sig_on, sig_off

cimport sage.libs.linbox.givaro as givaro
cimport sage.libs.linbox.linbox as linbox

from sage.arith.misc import is_prime
from sage.data_structures.binary_search cimport *
from sage.ext.stdsage cimport PY_NEW
from sage.libs.flint.fmpz cimport fmpz_get_mpz, fmpz_set_mpz
from sage.libs.flint.fmpz_mat cimport fmpz_mat_entry
from sage.libs.gmp.mpz cimport mpz_set
from sage.libs.linbox.conversion cimport (get_method,
                                          METHOD_DEFAULT,
                                          METHOD_DENSE_ELIMINATION,
                                          METHOD_SPARSE_ELIMINATION,
                                          METHOD_BLACKBOX,
                                          METHOD_WIEDEMANN,
                                          new_linbox_matrix_modn_sparse,
                                          new_linbox_matrix_integer_sparse,
                                          new_linbox_vector_integer_dense,
                                          new_sage_vector_integer_dense)
from sage.matrix.args cimport SparseEntry, MatrixArgs_init
from sage.matrix.matrix2 import Matrix as Matrix2
from sage.matrix.matrix_dense cimport Matrix_dense
from sage.matrix.matrix_integer_dense cimport Matrix_integer_dense
from sage.matrix.matrix_sparse cimport Matrix_sparse
from sage.misc.verbose import verbose, get_verbose
from sage.modules.vector_integer_dense cimport Vector_integer_dense
from sage.modules.vector_integer_sparse cimport *
from sage.modules.vector_modn_sparse cimport *
from sage.rings.fast_arith cimport arith_int
from sage.rings.finite_rings.integer_mod cimport IntegerMod_int, IntegerMod_abstract
from sage.rings.integer cimport Integer
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.structure.element cimport Matrix

################
# TODO: change this to use extern cdef's methods.
cdef arith_int ai
ai = arith_int()
################

# The 46341 below is because the mod-n sparse code still uses
# int's, even on 64-bit computers.  Improving this is
# Github Issue #12679.
MAX_MODULUS = 46341

cdef class Matrix_modn_sparse(Matrix_sparse):
    def __cinit__(self):
        nr = self._nrows
        nc = self._ncols
        cdef int p = self._base_ring.order()
        self.p = p

        self.rows = <c_vector_modint*>check_calloc(nr, sizeof(c_vector_modint))

        cdef Py_ssize_t i
        for i in range(nr):
            init_c_vector_modint(&self.rows[i], p, nc, 0)

    def __dealloc__(self):
        cdef Py_ssize_t i
        if self.rows:
            for i in range(self._nrows):
                clear_c_vector_modint(&self.rows[i])
            sig_free(self.rows)

    def __init__(self, parent, entries=None, copy=None, bint coerce=True):
        r"""
        Create a sparse matrix over the integers modulo ``n``.

        INPUT:

        - ``parent`` -- a matrix space over the integers modulo ``n``

        - ``entries`` -- see :func:`matrix`

        - ``copy`` -- ignored (for backwards compatibility)

        - ``coerce`` -- if ``False``, assume without checking that the
          entries lie in the base ring
        """
        ma = MatrixArgs_init(parent, entries)
        for t in ma.iter(coerce, True):
            se = <SparseEntry>t
            z = se.entry
            if z:
                set_entry(&self.rows[se.i], se.j, z)

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, value):
        set_entry(&self.rows[i], j, (<IntegerMod_int> value).ivalue)

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        cdef IntegerMod_int n
        n =  IntegerMod_int.__new__(IntegerMod_int)
        IntegerMod_abstract.__init__(n, self._base_ring)
        n.ivalue = get_entry(&self.rows[i], j)
        return n

    cdef bint get_is_zero_unsafe(self, Py_ssize_t i, Py_ssize_t j) except -1:
        """
        Return 1 if the entry ``(i, j)`` is zero, otherwise 0.

        EXAMPLES::

            sage: M = matrix(GF(13), [[0,1,0],[0,0,0]], sparse=True)
            sage: M.zero_pattern_matrix()  # indirect doctest
            [1 0 1]
            [1 1 1]
        """
        return is_entry_zero_unsafe(&self.rows[i], j)

    def _dict(self):
        """
        Unsafe version of the dict method, mainly for internal use. This
        may return the dict of elements, but as an *unsafe* reference to
        the underlying dict of the object. It might be dangerous if you
        change entries of the returned dict.

        EXAMPLES::

            sage: MS = MatrixSpace(GF(13), 50, 50, sparse=True)
            sage: m = MS._random_nonzero_element(density=0.002)
            sage: d = m._dict()
            sage: for i in range(50):
            ....:     for j in range(50):
            ....:         if m[i, j] != 0:
            ....:             assert m[i, j] == d[i, j]
            ....:         else:
            ....:             assert (i, j) not in d

        TESTS::

            sage: [i, j] = list(d.keys())[0]
            sage: parent(m._dict()[i, j])
            Finite Field of size 13
        """
        d = self.fetch('dict')
        if d is not None:
            return d

        cdef Py_ssize_t i, j
        d = {}
        cdef IntegerMod_int n
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self.rows[i].num_nonzero:
                n = IntegerMod_int.__new__(IntegerMod_int)
                IntegerMod_abstract.__init__(n, self._base_ring)
                n.ivalue = self.rows[i].entries[j]
                d[(int(i),int(self.rows[i].positions[j]))] = n
        self.cache('dict', d)
        return d

    def _pickle(self):
        """
        TESTS::

            sage: M = Matrix( GF(2), [[1,1,1,1,0,0,0,0,0,0]], sparse=True )
            sage: loads(dumps(M))
            [1 1 1 1 0 0 0 0 0 0]
            sage: loads(dumps(M)) == M
            True
        """
        return self._dict(), 1

    def _unpickle(self, data, version):
        if version == 1:
            self.__init__(self.parent(), data, coerce=False)
        else:
            raise ValueError("unknown matrix format")

    cdef Matrix _matrix_times_matrix_(self, Matrix _right):
        """
        This code is implicitly called for multiplying ``self`` by another
        sparse matrix.

        EXAMPLES::

            sage: a = matrix(GF(43), 3, 3, range(9), sparse=True)
            sage: b = matrix(GF(43), 3, 3, range(10,19), sparse=True)
            sage: a*b
            [ 2  5  8]
            [33  2 14]
            [21 42 20]
            sage: a*a
            [15 18 21]
            [42 11 23]
            [26  4 25]
            sage: c = matrix(GF(43), 3, 8, range(24), sparse=True)
            sage: a*c
            [40  0  3  6  9 12 15 18]
            [26 38  7 19 31  0 12 24]
            [12 33 11 32 10 31  9 30]

        Even though sparse and dense matrices are represented
        differently, they still compare as equal if they have the
        same entries::

            sage: a*b == a._matrix_times_matrix_dense(b)
            True
            sage: d = matrix(GF(43), 3, 8, range(24))
            sage: a*c == a*d
            True

        TESTS:

        The following shows that :issue:`23669` has been addressed::

            sage: p = next_prime(2**15)
            sage: M = Matrix(GF(p), 1,3, lambda i,j: -1, sparse=True); M
            [32770 32770 32770]
            sage: M*M.transpose() # previously returned [32738]
            [3]
        """
        cdef Matrix_modn_sparse right, ans
        right = _right

        cdef c_vector_modint* v

        # Build a table that gives the nonzero positions in each column of right
        cdef list nonzero_positions_in_columns = [set() for _ in range(right._ncols)]
        cdef Py_ssize_t i, j, k
        for i in range(right._nrows):
            v = &(right.rows[i])
            for j in range(v.num_nonzero):
                (<set> nonzero_positions_in_columns[v.positions[j]]).add(i)

        ans = self.new_matrix(self._nrows, right._ncols)

        # Now do the multiplication, getting each row completely before filling it in.
        cdef int x, y, s
        cdef set c

        for i in range(self._nrows):
            v = &(self.rows[i])
            for j in range(right._ncols):
                s = 0
                c = <set> nonzero_positions_in_columns[j]
                for k in range(v.num_nonzero):
                    if v.positions[k] in c:
                        y = get_entry(&right.rows[v.positions[k]], j)
                        x = v.entries[k] * y
                        s = (s + x) % self.p
                set_entry(&ans.rows[i], j, s)
        return ans

    def _matrix_times_matrix_dense(self, Matrix _right):
        """
        Multiply ``self`` by the sparse matrix ``_right``, and return the
        result as a dense matrix.

        EXAMPLES::

            sage: a = matrix(GF(10007), 2, [1,2,3,4], sparse=True)
            sage: b = matrix(GF(10007), 2, 3, [1..6], sparse=True)
            sage: a * b
            [ 9 12 15]
            [19 26 33]
            sage: c = a._matrix_times_matrix_dense(b); c
            [ 9 12 15]
            [19 26 33]
            sage: type(c)
            <class 'sage.matrix.matrix_modn_dense_double.Matrix_modn_dense_double'>

            sage: a = matrix(GF(2), 20, 20, sparse=True)
            sage: a*a == a._matrix_times_matrix_dense(a)
            True
            sage: type(a._matrix_times_matrix_dense(a))
            <class 'sage.matrix.matrix_mod2_dense.Matrix_mod2_dense'>
        """
        cdef Matrix_modn_sparse right
        cdef Matrix_dense ans
        right = _right

        cdef c_vector_modint* v

        # Build a table that gives the nonzero positions in each column of right
        nonzero_positions_in_columns = [set([]) for _ in range(right._ncols)]
        cdef Py_ssize_t i, j, k
        for i from 0 <= i < right._nrows:
            v = &(right.rows[i])
            for j from 0 <= j < right.rows[i].num_nonzero:
                nonzero_positions_in_columns[v.positions[j]].add(i)

        ans = self.new_matrix(self._nrows, right._ncols, sparse=False)

        # Now do the multiplication, getting each row completely before filling it in.
        cdef int x, y, s

        for i from 0 <= i < self._nrows:
            v = &self.rows[i]
            for j from 0 <= j < right._ncols:
                s = 0
                c = nonzero_positions_in_columns[j]
                for k from 0 <= k < v.num_nonzero:
                    if v.positions[k] in c:
                        y = get_entry(&right.rows[v.positions[k]], j)
                        x = v.entries[k] * y
                        s = (s + x)%self.p
                ans.set_unsafe_int(i, j, s)
                #ans._matrix[i][j] = s
        return ans

    def swap_rows(self, r1, r2):
        self.check_bounds_and_mutability(r1,0)
        self.check_bounds_and_mutability(r2,0)
        self.swap_rows_c(r1, r2)

    cdef swap_rows_c(self, Py_ssize_t n1, Py_ssize_t n2):
        """
        Swap the rows in positions n1 and n2. No bounds checking.
        """
        cdef c_vector_modint tmp
        tmp = self.rows[n1]
        self.rows[n1] = self.rows[n2]
        self.rows[n2] = tmp

    cpdef _echelon_in_place(self, str algorithm):
        """
        Replace ``self`` by its reduction to reduced row echelon form.

        ALGORITHM: We use Gauss elimination, in a slightly intelligent way,
        in that we clear each column using a row with the minimum number of
        nonzero entries.

        TODO: Implement switching to a dense method when the matrix gets
        dense.
        """
        from sage.misc.verbose import verbose, get_verbose
        x = self.fetch('in_echelon_form')
        if x is not None and x:
            return  # already known to be in echelon form
        self.check_mutability()

        cdef Py_ssize_t i, r, c, min, min_row, start_row
        cdef int a_inverse, b, do_verb
        cdef c_vector_modint tmp
        start_row = 0
        pivots = []
        fifth = self._ncols / 10 + 1
        tm = verbose(caller_name = 'sparse_matrix_pyx matrix_modint echelon')
        do_verb = (get_verbose() >= 2)

        for c from 0 <= c < self._ncols:
            if do_verb and (c % fifth == 0 and c>0):
                tm = verbose('on column %s of %s'%(c, self._ncols),
                             level = 2,
                             caller_name = 'matrix_modn_sparse echelon')
            #end if
            min = self._ncols + 1
            min_row = -1
            for r from start_row <= r < self._nrows:
                if self.rows[r].num_nonzero > 0 and self.rows[r].num_nonzero < min:
                    # Since there is at least one nonzero entry, the first entry
                    # of the positions list is defined.  It is the first position
                    # of a nonzero entry, and it equals c precisely if row r
                    # is a row we could use to clear column c.
                    if self.rows[r].positions[0] == c:
                        min_row = r
                        min = self.rows[r].num_nonzero
                    #endif
                #endif
            #endfor
            if min_row != -1:
                r = min_row
                # print("min number of entries in a pivoting row = ", min)
                pivots.append(c)
                # Since we can use row r to clear column c, the
                # entry in position c in row r must be the first nonzero entry.
                a = self.rows[r].entries[0]
                if a != 1:
                    a_inverse = ai.c_inverse_mod_int(a, self.p)
                    scale_c_vector_modint(&self.rows[r], a_inverse)
                self.swap_rows_c(r, start_row)
                sig_on()
                for i from 0 <= i < self._nrows:
                    if i != start_row:
                        b = get_entry(&self.rows[i], c)
                        if b != 0:
                            add_c_vector_modint_init(&tmp, &self.rows[i],
                                                     &self.rows[start_row], self.p - b)
                            clear_c_vector_modint(&self.rows[i])
                            self.rows[i] = tmp
                sig_off()
                start_row = start_row + 1

        self.cache('pivots', tuple(pivots))
        self.cache('in_echelon_form', True)

    def _nonzero_positions_by_row(self, copy=True):
        """
        Return the list of pairs (i,j) such that self[i,j] != 0.

        It is safe to change the resulting list (unless you give the option copy=False).

        EXAMPLES::

            sage: M = Matrix(GF(7), [[0,0,0,1,0,0,0,0],[0,1,0,0,0,0,1,0]], sparse=True); M
            [0 0 0 1 0 0 0 0]
            [0 1 0 0 0 0 1 0]
            sage: M.nonzero_positions()
            [(0, 3), (1, 1), (1, 6)]
        """
        x = self.fetch('nonzero_positions')
        if x is not None:
            if copy:
                return list(x)
            return x
        nzp = []
        cdef Py_ssize_t i, j
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self.rows[i].num_nonzero:
                nzp.append((i,self.rows[i].positions[j]))
        self.cache('nonzero_positions', nzp)
        if copy:
            return list(nzp)
        return nzp

    def density(self):
        """
        Return the density of ``self``, i.e., the ratio of the number of
        nonzero entries of ``self`` to the total size of ``self``.

        EXAMPLES::

            sage: A = matrix(QQ,3,3,[0,1,2,3,0,0,6,7,8],sparse=True)
            sage: A.density()
            2/3

        Notice that the density parameter does not ensure the density
        of a matrix; it is only an upper bound.

        ::

            sage: A = random_matrix(GF(127), 200, 200, density=0.3, sparse=True)
            sage: density_sum = float(A.density())
            sage: total = 1
            sage: expected_density = 1.0 - (199/200)^60
            sage: expected_density
            0.2597...
            sage: while abs(density_sum/total - expected_density) > 0.001:
            ....:     A = random_matrix(GF(127), 200, 200, density=0.3, sparse=True)
            ....:     density_sum += float(A.density())
            ....:     total += 1
        """
        cdef Py_ssize_t i, nonzero_entries

        nonzero_entries = 0
        for i from 0 <= i < self._nrows:
            nonzero_entries += self.rows[i].num_nonzero

        return ZZ(nonzero_entries) / ZZ(self._nrows*self._ncols)

    def transpose(self):
        """
        Return the transpose of ``self``.

        EXAMPLES::

            sage: A = matrix(GF(127),3,3,[0,1,0,2,0,0,3,0,0],sparse=True)
            sage: A
            [0 1 0]
            [2 0 0]
            [3 0 0]
            sage: A.transpose()
            [0 2 3]
            [1 0 0]
            [0 0 0]

        ``.T`` is a convenient shortcut for the transpose::

            sage: A.T
            [0 2 3]
            [1 0 0]
            [0 0 0]
        """
        cdef int i, j
        cdef c_vector_modint row
        cdef Matrix_modn_sparse B

        B = self.new_matrix(nrows = self.ncols(), ncols = self.nrows())
        for i from 0 <= i < self._nrows:
            row = self.rows[i]
            for j from 0 <= j < row.num_nonzero:
                set_entry(&B.rows[row.positions[j]], i, row.entries[j])
        if self._subdivisions is not None:
            row_divs, col_divs = self.subdivisions()
            B.subdivide(col_divs, row_divs)
        return B

    def matrix_from_rows(self, rows):
        """
        Return the matrix constructed from ``self`` using rows with indices in
        the rows list.

        INPUT:

        - ``rows`` -- list or tuple of row indices

        EXAMPLES::

            sage: M = MatrixSpace(GF(127),3,3,sparse=True)
            sage: A = M(range(9)); A
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: A.matrix_from_rows([2,1])
            [6 7 8]
            [3 4 5]
        """
        cdef int i,k
        cdef Matrix_modn_sparse A
        cdef c_vector_modint row

        if not isinstance(rows, (list, tuple)):
            rows = list(rows)

        A = self.new_matrix(nrows = len(rows))

        k = 0
        for ii in rows:
            i = ii
            if i < 0 or i >= self.nrows():
                raise IndexError("row %s out of range" % i)

            row = self.rows[i]
            for j from 0 <= j < row.num_nonzero:
                set_entry(&A.rows[k], row.positions[j], row.entries[j])
            k += 1
        return A

    def matrix_from_columns(self, cols):
        """
        Return the matrix constructed from ``self`` using columns with indices
        in the columns list.

        EXAMPLES::

            sage: M = MatrixSpace(GF(127),3,3,sparse=True)
            sage: A = M(range(9)); A
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: A.matrix_from_columns([2,1])
            [2 1]
            [5 4]
            [8 7]
        """
        cdef int i,j
        cdef Matrix_modn_sparse A
        cdef c_vector_modint row

        if not isinstance(cols, (list, tuple)):
            cols = list(cols)

        A = self.new_matrix(ncols = len(cols))

        cols = dict(zip([int(e) for e in cols],range(len(cols))))

        for i from 0 <= i < self.nrows():
            row = self.rows[i]
            for j from 0 <= j < row.num_nonzero:
                if int(row.positions[j]) in cols:
                    set_entry(&A.rows[i], cols[int(row.positions[j])], row.entries[j])
        return A

    def _rank_det_linbox(self):
        """
        Return the rank and determinant using linbox.

        .. NOTE::

            This method does not perform any caching contrarily to
            :meth:`determinant` and :meth:`rank`.

        EXAMPLES::

            sage: m = matrix(Zmod(13), 1, sparse=True)
            sage: m[0,0] = 0
            sage: m._rank_det_linbox()
            (0, 0)
            sage: for i in range(1, 13):
            ....:     m[0,0] = i
            ....:     assert m._rank_det_linbox() == (1, i)

            sage: m = matrix(GF(5), 2, sparse=True)
            sage: m[0,0] = 1
            sage: m[0,1] = 2
            sage: m[1,0] = 1
            sage: m[1,1] = 3
            sage: m._rank_det_linbox()
            (2, 1)
            sage: m
            [1 2]
            [1 3]

        TESTS::

            sage: matrix(Zmod(3), 0, sparse=True)._rank_det_linbox()
            (0, 1)
        """
        if self._nrows == 0 or self._ncols == 0:
            # TODO: bug in linbox (gives segfault)
            return 0, self.base_ring().one()

        cdef size_t A_rank = 0
        cdef uint64_t A_det = 0

        if not is_prime(self.p):
            raise TypeError("only GF(p) supported via LinBox")

        cdef givaro.Modular_uint64 * F = new givaro.Modular_uint64(<uint64_t> self.p)
        cdef linbox.SparseMatrix_Modular_uint64 * A
        A = new_linbox_matrix_modn_sparse(F[0], self)

        cdef linbox.GaussDomain_Modular_uint64 * dom = new linbox.GaussDomain_Modular_uint64(F[0])

        # NOTE: see above for the reason of casting...
        if A.rowdim() >= <size_t> UINT_MAX or A.coldim() >= <size_t> UINT_MAX:
            raise ValueError("row/column size unsupported in LinBox")

        dom.InPlaceLinearPivoting(A_rank, A_det, A[0], A.rowdim(), A.coldim())

        del A
        del F
        del dom

        return <long> A_rank, self.base_ring()(A_det)

    def rank(self, algorithm=None):
        """
        Return the rank of this matrix.

        INPUT:

        - ``algorithm`` -- either ``'linbox'`` (only available for
          matrices over prime fields) or ``'generic'``

        EXAMPLES::

            sage: A = matrix(GF(127), 2, 2, sparse=True)
            sage: A[0,0] = 34
            sage: A[0,1] = 102
            sage: A[1,0] = 55
            sage: A[1,1] = 74
            sage: A.rank()
            2

            sage: A._clear_cache()
            sage: A.rank(algorithm='generic')
            2
            sage: A._clear_cache()
            sage: A.rank(algorithm='hey')
            Traceback (most recent call last):
            ...
            ValueError: no algorithm 'hey'

        TESTS::

            sage: matrix(GF(3), 0, sparse=True).rank(algorithm='generic')
            0
            sage: matrix(GF(3), 0, sparse=True).rank(algorithm='linbox')
            0

            sage: for _ in range(50):
            ....:     nrows = randint(0, 100)
            ....:     ncols = randint(0, 100)
            ....:     p = random_prime(10000)
            ....:     M = MatrixSpace(GF(p), nrows, ncols, sparse=True)
            ....:     m = M.random_element()
            ....:     rank_linbox = m.rank(algorithm='linbox')
            ....:     rank_generic = m.rank(algorithm='generic')
            ....:     if rank_linbox != rank_generic:
            ....:         print(m)
            ....:         raise RuntimeError

        REFERENCES:

        - Jean-Guillaume Dumas and Gilles Villars. 'Computing the Rank
          of Large Sparse Matrices over Finite
          Fields'. Proc. CASC'2002, The Fifth International Workshop
          on Computer Algebra in Scientific Computing, Big Yalta,
          Crimea, Ukraine, 22-27 sept. 2002, Springer-Verlag,
          http://perso.ens-lyon.fr/gilles.villard/BIBLIOGRAPHIE/POSTSCRIPT/rankjgd.ps

        .. NOTE::

           For very sparse matrices Gaussian elimination is faster
           because it barely has anything to do. If the fill in needs to
           be considered, 'Symbolic Reordering' is usually much faster.
        """
        if self._nrows == 0 or self._ncols == 0:
            return 0

        if not is_prime(self.p):
            raise ArithmeticError("rank not well defined for matrices over general ring")

        x = self.fetch('rank')
        if x is not None:
            return x

        if algorithm is None or algorithm == "linbox":
            rank, det = self._rank_det_linbox()
            self.cache("rank", rank)
            self.cache("det", det)
            return rank

        elif algorithm == "generic":
            return Matrix2.rank(self)

        else:
            raise ValueError("no algorithm '%s'"%algorithm)

    def determinant(self, algorithm=None):
        r"""
        Return the determinant of this matrix.

        INPUT:

        - ``algorithm`` -- either ``'linbox'`` (default) or ``'generic'``

        EXAMPLES::

            sage: A = matrix(GF(3), 4, range(16), sparse=True)
            sage: B = identity_matrix(GF(3), 4, sparse=True)
            sage: (A + B).det()
            2
            sage: (A + B).det(algorithm='linbox')
            2
            sage: (A + B).det(algorithm='generic')
            2
            sage: (A + B).det(algorithm='hey')
            Traceback (most recent call last):
            ...
            ValueError: no algorithm 'hey'

            sage: matrix(GF(11), 1, 2, sparse=True).det()
            Traceback (most recent call last):
            ...
            ValueError: self must be a square matrix

        TESTS::

            sage: matrix(GF(3), 0, sparse=True).det(algorithm='generic')
            1
            sage: matrix(GF(3), 0, sparse=True).det(algorithm='linbox')
            1

            sage: for _ in range(100):
            ....:     dim = randint(0, 50)
            ....:     p = random_prime(10000)
            ....:     M = MatrixSpace(GF(p), dim, sparse=True)
            ....:     m = M.random_element()
            ....:     det_linbox = m.det(algorithm='linbox')
            ....:     det_generic = m.det(algorithm='generic')
            ....:     assert parent(det_linbox) == m.base_ring()
            ....:     assert parent(det_generic) == m.base_ring()
            ....:     if det_linbox != det_generic:
            ....:         print(m)
            ....:         raise RuntimeError
        """
        if self._nrows != self._ncols:
            raise ValueError("self must be a square matrix")
        if self._nrows == 0:
            return self.base_ring().one()

        d = self.fetch('det')
        if d is not None:
            return d

        if algorithm is None or algorithm == "linbox":
            r, d = self._rank_det_linbox()
            self.cache('rank', r)
            self.cache('det', d)
            return d
        elif algorithm == 'generic':
            d = Matrix_sparse.determinant(self)
            self.cache('det', d)
            return d
        else:
            raise ValueError("no algorithm '%s'"%algorithm)
