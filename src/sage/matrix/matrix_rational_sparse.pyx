"""
Sparse rational matrices

AUTHORS:

- William Stein (2007-02-21)
- Soroosh Yazdani (2007-02-21)

TESTS::

    sage: a = matrix(QQ,2,range(4), sparse=True)
    sage: TestSuite(a).run()
    sage: matrix(QQ,0,0,sparse=True).inverse()
    []
"""

#*****************************************************************************
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cysignals.signals cimport sig_on, sig_off
from cysignals.memory cimport check_calloc, sig_free

from sage.data_structures.binary_search cimport *
from sage.modules.vector_integer_sparse cimport *
from sage.modules.vector_rational_sparse cimport *

from cpython.sequence cimport *

from sage.rings.rational cimport Rational
from sage.rings.integer  cimport Integer
from sage.matrix.matrix cimport Matrix
from sage.matrix.args cimport SparseEntry, MatrixArgs_init

from sage.libs.gmp.mpz cimport *
from sage.libs.gmp.mpq cimport *

from sage.libs.flint.fmpq cimport fmpq_set_mpq
from sage.libs.flint.fmpq_mat cimport fmpq_mat_entry

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

cimport sage.structure.element

import sage.matrix.matrix_space

from sage.matrix.matrix_integer_sparse cimport Matrix_integer_sparse
from sage.matrix.matrix_rational_dense cimport Matrix_rational_dense


cdef class Matrix_rational_sparse(Matrix_sparse):
    def __cinit__(self):
        self._matrix = <mpq_vector*>check_calloc(self._nrows, sizeof(mpq_vector))
        # initialize the rows
        cdef Py_ssize_t i
        for i in range(self._nrows):
            mpq_vector_init(&self._matrix[i], self._ncols, 0)

    def __dealloc__(self):
        cdef Py_ssize_t i
        if self._matrix is not NULL:
            for i in range(self._nrows):
                mpq_vector_clear(&self._matrix[i])
            sig_free(self._matrix)

    def __init__(self, parent, entries=None, copy=None, bint coerce=True):
        r"""
        Create a sparse matrix over the rational numbers.

        INPUT:

        - ``parent`` -- a matrix space over `\QQ`

        - ``entries`` -- see :func:`matrix`

        - ``copy`` -- ignored (for backwards compatibility)

        - ``coerce`` -- if ``False``, assume without checking that the
          entries are of type :class:`Rational`
        """
        ma = MatrixArgs_init(parent, entries)
        cdef Rational z
        for t in ma.iter(coerce, True):
            se = <SparseEntry>t
            z = <Rational>se.entry
            if z:
                mpq_vector_set_entry(&self._matrix[se.i], se.j, z.value)

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, x):
        mpq_vector_set_entry(&self._matrix[i], j, (<Rational> x).value)

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        cdef Rational x
        x = Rational()
        mpq_vector_get_entry(x.value, &self._matrix[i], j)
        return x

    cdef bint get_is_zero_unsafe(self, Py_ssize_t i, Py_ssize_t j) except -1:
        """
        Return 1 if the entry ``(i, j)`` is zero, otherwise 0.

        EXAMPLES::

            sage: M = matrix(QQ, [[0,1,0],[0,0,0]], sparse=True)
            sage: M.zero_pattern_matrix()  # indirect doctest
            [1 0 1]
            [1 1 1]
        """
        return mpq_vector_is_entry_zero_unsafe(&self._matrix[i], j)

    def add_to_entry(self, Py_ssize_t i, Py_ssize_t j, elt):
        r"""
        Add ``elt`` to the entry at position ``(i, j)``.

        EXAMPLES::

            sage: m = matrix(QQ, 2, 2, sparse=True)
            sage: m.add_to_entry(0, 0, -1/3)
            sage: m
            [-1/3    0]
            [   0    0]
            sage: m.add_to_entry(0, 0, 1/3)
            sage: m
            [0 0]
            [0 0]
            sage: m.nonzero_positions()
            []
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

        cdef mpq_t z
        mpq_init(z)
        mpq_vector_get_entry(z, &self._matrix[i], j)
        mpq_add(z, z, (<Rational>elt).value)
        mpq_vector_set_entry(&self._matrix[i], j, z)
        mpq_clear(z)

    ########################################################################
    # LEVEL 2 functionality
    #   * def _pickle
    #   * def _unpickle
    #   * cdef _add_
    #   * cdef _sub_
    #   * cdef _mul_
    #   * cpdef _richcmp_
    #   * __neg__
    #   * __invert__
    #   * __copy__
    #   * _multiply_classical
    #   * _matrix_times_matrix_
    #   * _list -- list of underlying elements (need not be a copy)
    #   * x _dict -- sparse dictionary of underlying elements (need not be a copy)

    cdef sage.structure.element.Matrix _matrix_times_matrix_(self, sage.structure.element.Matrix _right):
        cdef Matrix_rational_sparse right, ans
        right = _right

        cdef mpq_vector* v

        # Build a table that gives the nonzero positions in each column of right
        cdef list nonzero_positions_in_columns = [set() for _ in range(right._ncols)]
        cdef Py_ssize_t i, j, k
        for i in range(right._nrows):
            v = &(right._matrix[i])
            for j in range(v.num_nonzero):
                (<set> nonzero_positions_in_columns[v.positions[j]]).add(i)
        # pre-computes the list of nonzero columns of right
        cdef list right_indices
        right_indices = [j for j in range(right._ncols)
                         if nonzero_positions_in_columns[j]]

        ans = self.new_matrix(self._nrows, right._ncols)

        # Now do the multiplication, getting each row completely before filling it in.
        cdef set c
        cdef mpq_t x, y, s
        mpq_init(x)
        mpq_init(y)
        mpq_init(s)
        for i in range(self._nrows):
            v = &(self._matrix[i])
            if not v.num_nonzero:
                continue
            for j in right_indices:
                mpq_set_si(s, 0, 1)
                c = <set> nonzero_positions_in_columns[j]
                for k in range(v.num_nonzero):
                    if v.positions[k] in c:
                        mpq_vector_get_entry(y, &right._matrix[v.positions[k]], j)
                        mpq_mul(x, v.entries[k], y)
                        mpq_add(s, s, x)
                mpq_vector_set_entry(&ans._matrix[i], j, s)

        mpq_clear(x)
        mpq_clear(y)
        mpq_clear(s)
        return ans

    def _matrix_times_matrix_dense(self, sage.structure.element.Matrix _right):
        """
        Do the sparse matrix multiply, but return a dense matrix as the result.

        EXAMPLES::

            sage: a = matrix(QQ, 2, [1,2,3,4], sparse=True)
            sage: b = matrix(QQ, 2, 3, [1..6], sparse=True)
            sage: a * b
            [ 9 12 15]
            [19 26 33]
            sage: c = a._matrix_times_matrix_dense(b); c
            [ 9 12 15]
            [19 26 33]
            sage: type(c)
            <class 'sage.matrix.matrix_rational_dense.Matrix_rational_dense'>
        """
        cdef Matrix_rational_sparse right
        cdef Matrix_rational_dense ans
        right = _right

        cdef mpq_vector* v

        # Build a table that gives the nonzero positions in each column of right
        nonzero_positions_in_columns = [set([]) for _ in range(right._ncols)]
        cdef Py_ssize_t i, j, k
        for i in range(right._nrows):
            v = &(right._matrix[i])
            for j in range(right._matrix[i].num_nonzero):
                nonzero_positions_in_columns[v.positions[j]].add(i)

        ans = self.new_matrix(self._nrows, right._ncols, sparse=False)

        # Now do the multiplication, getting each row completely before filling it in.
        cdef mpq_t x, y, s
        mpq_init(x)
        mpq_init(y)
        mpq_init(s)
        for i in range(self._nrows):
            v = &self._matrix[i]
            for j in range(right._ncols):
                mpq_set_si(s, 0, 1)
                c = nonzero_positions_in_columns[j]
                for k in range(v.num_nonzero):
                    if v.positions[k] in c:
                        mpq_vector_get_entry(y, &right._matrix[v.positions[k]], j)
                        mpq_mul(x, v.entries[k], y)
                        mpq_add(s, s, x)
                fmpq_set_mpq(fmpq_mat_entry(ans._matrix, i, j), s)

        mpq_clear(x)
        mpq_clear(y)
        mpq_clear(s)
        return ans

    ########################################################################
    # def _pickle(self):
    # def _unpickle(self, data, int version):   # use version >= 0
    # cpdef _add_(self, right):
    # cdef _mul_(self, Matrix right):
    # cpdef _richcmp_(self, Matrix right, int op):
    # def __neg__(self):
    # def __invert__(self):
    # def __copy__(self):
    # def _multiply_classical(left, matrix.Matrix _right):
    # def _list(self):

    # TODO
    ##     cpdef _lmul_(self, Element right):
    ##         """
    ##         EXAMPLES::
    ##
    ##             sage: a = matrix(QQ,2,range(6))
    ##             sage: (3/4) * a
    ##             [   0  3/4  3/2]
    ##             [ 9/4    3 15/4]
    ##         """

    def _dict(self):
        """
        Unsafe version of the dict method, mainly for internal use.

        This may return the dict of elements, but as an *unsafe*
        reference to the underlying dict of the object.  It might
        be dangerous if you change entries of the returned dict.
        """
        d = self.fetch('dict')
        if d is not None:
            return d

        cdef Py_ssize_t i, j
        d = {}
        for i in range(self._nrows):
            for j in range(self._matrix[i].num_nonzero):
                x = Rational()
                mpq_set((<Rational>x).value, self._matrix[i].entries[j])
                d[(int(i), int(self._matrix[i].positions[j]))] = x
        self.cache('dict', d)
        return d

    ########################################################################
    # LEVEL 3 functionality (Optional)
    #    * cdef _sub_
    #    * __deepcopy__
    #    * __invert__
    #    * Matrix windows -- only if you need strassen for that base
    #    * Other functions (list them here):
    ########################################################################

    def _nonzero_positions_by_row(self, copy=True):
        """
        Return the list of pairs (i,j) such that self[i,j] != 0.

        It is safe to change the resulting list (unless you give the option copy=False).

        EXAMPLES::

            sage: M = Matrix(QQ, [[0,0,0,1,0,0,0,0],[0,1,0,0,0,0,1,0]], sparse=True); M
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
            for j from 0 <= j < self._matrix[i].num_nonzero:
                nzp.append((i,self._matrix[i].positions[j]))
        self.cache('nonzero_positions', nzp)
        if copy:
            return list(nzp)
        return nzp

    def height(self):
        """
        Return the height of this matrix, which is the least common
        multiple of all numerators and denominators of elements of
        this matrix.

        OUTPUT: integer

        EXAMPLES::

            sage: b = matrix(QQ,2,range(6), sparse=True); b[0,0]=-5007/293; b
            [-5007/293         1         2]
            [        3         4         5]
            sage: b.height()
            5007
        """
        cdef Integer z = Integer.__new__(Integer)
        self.mpz_height(z.value)
        return z

    cdef int mpz_height(self, mpz_t height) except -1:
        cdef mpz_t x, h
        mpz_init(x)
        mpz_init_set_si(h, 0)
        cdef int i, j
        sig_on()
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._matrix[i].num_nonzero:
                mpq_get_num(x, self._matrix[i].entries[j])
                mpz_abs(x, x)
                if mpz_cmp(h,x) < 0:
                    mpz_set(h,x)
                mpq_get_den(x, self._matrix[i].entries[j])
                mpz_abs(x, x)
                if mpz_cmp(h,x) < 0:
                    mpz_set(h,x)
        sig_off()
        mpz_set(height, h)
        mpz_clear(h)
        mpz_clear(x)
        return 0

    cdef int mpz_denom(self, mpz_t d) except -1:
        mpz_set_si(d, 1)
        cdef Py_ssize_t i, j

        sig_on()
        for i in range(self._nrows):
            for j in range(self._matrix[i].num_nonzero):
                mpz_lcm(d, d, mpq_denref(self._matrix[i].entries[j]))
        sig_off()
        return 0

    def denominator(self):
        """
        Return the denominator of this matrix.

        OUTPUT: Sage Integer

        EXAMPLES::

            sage: b = matrix(QQ,2,range(6)); b[0,0]=-5007/293; b
            [-5007/293         1         2]
            [        3         4         5]
            sage: b.denominator()
            293
        """
        cdef Integer z = Integer.__new__(Integer)
        self.mpz_denom(z.value)
        return z

    def _clear_denom(self):
        """
        OUTPUT: ``d*self, D``

        The product D*self is a matrix over ZZ

        EXAMPLES::

            sage: a = matrix(QQ,3,[-2/7, -1/4, -2, 0, 1/7, 1, 0, 1/2, 1/5],sparse=True)
            sage: a.denominator()
            140
            sage: a._clear_denom()
            (
            [ -40  -35 -280]
            [   0   20  140]
            [   0   70   28], 140
            )
        """
        cdef Integer D
        cdef Py_ssize_t i, j
        cdef Matrix_integer_sparse A
        cdef mpz_t t
        cdef mpq_vector* v

        D = Integer()
        self.mpz_denom(D.value)

        MZ = sage.matrix.matrix_space.MatrixSpace(ZZ, self._nrows, self._ncols, sparse=True)
        A = MZ.zero_matrix().__copy__()

        mpz_init(t)
        sig_on()
        for i from 0 <= i < self._nrows:
            v = &(self._matrix[i])
            for j from 0 <= j < v.num_nonzero:
                mpz_divexact(t, D.value, mpq_denref(v.entries[j]))
                mpz_mul(t, t, mpq_numref(v.entries[j]))
                mpz_vector_set_entry(&(A._matrix[i]), v.positions[j], t)
        sig_off()
        mpz_clear(t)
        return A, D

    ################################################
    # Echelon form
    ################################################
    def echelonize(self, height_guess=None, proof=True, **kwds):
        """
        Transform the matrix ``self`` into reduced row echelon form
        in place.

        INPUT:

        - ``height_guess``, ``proof``, ``**kwds`` -- all passed to the multimodular
          algorithm; ignored by the `p`-adic algorithm

        OUTPUT:

        Nothing. The matrix ``self`` is transformed into reduced row
        echelon form in place.

        ALGORITHM: a multimodular algorithm.

        EXAMPLES::

            sage: a = matrix(QQ, 4, range(16), sparse=True); a[0,0] = 1/19; a[0,1] = 1/5; a
            [1/19  1/5    2    3]
            [   4    5    6    7]
            [   8    9   10   11]
            [  12   13   14   15]
            sage: a.echelonize(); a
            [      1       0       0 -76/157]
            [      0       1       0  -5/157]
            [      0       0       1 238/157]
            [      0       0       0       0]

        :issue:`10319` has been fixed::

            sage: m = Matrix(QQ, [1], sparse=True); m.echelonize()
            sage: m = Matrix(QQ, [1], sparse=True); m.echelonize(); m
            [1]
        """

        x = self.fetch('in_echelon_form')
        if x is not None:
            return  # already known to be in echelon form
        self.check_mutability()

        pivots = self._echelonize_multimodular(height_guess, proof, **kwds)

        self.cache('in_echelon_form', True)
        self.cache('echelon_form', self)
        self.cache('pivots', pivots)
        self.cache('rank', len(pivots))

    def echelon_form(self, algorithm='default',
                     height_guess=None, proof=True, **kwds):
        """
        INPUT:

        - ``height_guess``, ``proof``, ``**kwds`` -- all passed to the multimodular
          algorithm; ignored by the `p`-adic algorithm

        OUTPUT: ``self`` is no in reduced row echelon form

        EXAMPLES::

            sage: a = matrix(QQ, 4, range(16), sparse=True); a[0,0] = 1/19; a[0,1] = 1/5; a
            [1/19  1/5    2    3]
            [   4    5    6    7]
            [   8    9   10   11]
            [  12   13   14   15]
            sage: a.echelon_form()
            [      1       0       0 -76/157]
            [      0       1       0  -5/157]
            [      0       0       1 238/157]
            [      0       0       0       0]
        """
        label = 'echelon_form_%s'%algorithm
        x = self.fetch(label)
        if x is not None:
            return x
        if self.fetch('in_echelon_form'):
            return self

        E, pivots = self._echelon_form_multimodular(height_guess, proof=proof)

        self.cache(label, E)
        self.cache('pivots', pivots)
        return E

    # Multimodular echelonization algorithms
    def _echelonize_multimodular(self, height_guess=None, proof=True, **kwds):
        cdef Matrix_rational_sparse E
        E, pivots = self._echelon_form_multimodular(height_guess, proof=proof, **kwds)
        self.clear_cache()
        # Swap the data of E and self (effectively moving E to self)
        self._matrix, E._matrix = E._matrix, self._matrix
        return pivots

    def _echelon_form_multimodular(self, height_guess=None, proof=True):
        """
        Return reduced row-echelon form using a multi-modular
        algorithm.  Does not change ``self``.

        INPUT:

        - ``height_guess`` -- integer or ``None``
        - ``proof`` -- boolean (default: ``True``)
        """
        from sage.matrix.misc import matrix_rational_echelon_form_multimodular
        cdef Matrix E
        E, pivots = matrix_rational_echelon_form_multimodular(self,
                                 height_guess=height_guess, proof=proof)
        E._parent = self._parent
        return E, pivots

    def set_row_to_multiple_of_row(self, i, j, s):
        """
        Set row i equal to s times row j.

        EXAMPLES::

            sage: a = matrix(QQ,2,3,range(6), sparse=True); a
            [0 1 2]
            [3 4 5]
            sage: a.set_row_to_multiple_of_row(1,0,-3)
            sage: a
            [ 0  1  2]
            [ 0 -3 -6]
        """
        self.check_row_bounds_and_mutability(i, j)
        cdef Rational _s
        _s = Rational(s)
        mpq_vector_scalar_multiply(&self._matrix[i], &self._matrix[j], _s.value)

    def dense_matrix(self):
        """
        Return dense version of this matrix.

        EXAMPLES::

            sage: a = matrix(QQ,2,[1..4],sparse=True); type(a)
            <class 'sage.matrix.matrix_rational_sparse.Matrix_rational_sparse'>
            sage: type(a.dense_matrix())
            <class 'sage.matrix.matrix_rational_dense.Matrix_rational_dense'>
            sage: a.dense_matrix()
            [1 2]
            [3 4]

        Check that subdivisions are preserved when converting between
        dense and sparse matrices::

            sage: a.subdivide([1,1], [2])
            sage: b = a.dense_matrix().sparse_matrix().dense_matrix()
            sage: b.subdivisions() == a.subdivisions()
            True
        """
        cdef Matrix_rational_dense B
        cdef mpq_vector* v

        B = self.matrix_space(sparse=False).zero_matrix().__copy__()
        for i from 0 <= i < self._nrows:
            v = &(self._matrix[i])
            for j from 0 <= j < v.num_nonzero:
                fmpq_set_mpq(fmpq_mat_entry(B._matrix, i, v.positions[j]), v.entries[j])
        B.subdivide(self.subdivisions())
        return B

##     def _set_row_to_negative_of_row_of_A_using_subset_of_columns(self, Py_ssize_t i, Matrix A,
##                                                                  Py_ssize_t r, cols):
##         B = self.__copy__()
##         B.x_set_row_to_negative_of_row_of_A_using_subset_of_columns(i, A, r, cols)
##         cdef Py_ssize_t l
##         l = 0
##         for z in range(self.ncols()):
##             self[i,z] = 0
##         for k in cols:
##             self.set_unsafe(i,l,-A.get_unsafe(r,k))               #self[i,l] = -A[r,k]
##             l += 1
##         if self != B:
##             print("correct =\n", self.str())
##             print("wrong = \n", B.str())
##             print("diff = \n", (self-B).str())

    def _set_row_to_negative_of_row_of_A_using_subset_of_columns(self, Py_ssize_t i, Matrix A,
                                                                 Py_ssize_t r, cols,
                                                                 cols_index=None):
        """
        Set row i of ``self`` to -(row r of A), but where we only take the
        given column positions in that row of A.  Note that we *DO*
        zero out the other entries of ``self``'s row i.

        INPUT:

        - ``i`` -- integer, index into the rows of self
        - ``A`` -- a sparse matrix
        - ``r`` -- integer, index into rows of A
        - ``cols`` -- a *sorted* list of integers
        - ``cols_index`` -- (optional) set it to this to vastly speed up
          calls to this function::

                dict([(cols[i], i) for i in range(len(cols))])

        EXAMPLES::

            sage: a = matrix(QQ,2,3,range(6), sparse=True); a
            [0 1 2]
            [3 4 5]

        Note that the row is zeroed out before being set in the sparse case. ::

            sage: a._set_row_to_negative_of_row_of_A_using_subset_of_columns(0,a,1,[1,2])
            sage: a
            [-4 -5  0]
            [ 3  4  5]
        """
        # this function exists just because it is useful for modular symbols presentations.
        self.check_row_bounds_and_mutability(i,i)
        if r < 0 or r >= A.nrows():
            raise IndexError("invalid row")

        if not A.is_sparse():
            A = A.sparse_matrix()

        if A.base_ring() != QQ:
            A = A.change_ring(QQ)

        cdef Matrix_rational_sparse _A
        _A = A

        cdef Py_ssize_t l, n

        cdef mpq_vector *v
        cdef mpq_vector *w
        v = &self._matrix[i]
        w = &_A._matrix[r]

        if cols_index is None:
            cols_index = dict([(cols[i], i) for i in range(len(cols))])

        _cols = set(cols)
        pos = [i for i from 0 <= i < w.num_nonzero if w.positions[i] in _cols]
        n = len(pos)

        mpq_vector_clear(v)
        allocate_mpq_vector(v, n)
        v.num_nonzero = n
        v.degree = self._ncols

        for l from 0 <= l < n:
            v.positions[l] = cols_index[w.positions[pos[l]]]
            mpq_mul(v.entries[l], w.entries[pos[l]], minus_one)

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
            ....:                 [4, -1, 0, -6, 2]],
            ....:             sparse=True)
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

            sage: A = matrix(QQ, 0, 2, sparse=True)
            sage: A._right_kernel_matrix()[1]
            [1 0]
            [0 1]
            sage: A = matrix(QQ, 2, 0, sparse=True)
            sage: A._right_kernel_matrix()[1].parent()
            Full MatrixSpace of 0 by 0 dense matrices over Rational Field
            sage: A = zero_matrix(QQ, 4, 3,  sparse=True)
            sage: A._right_kernel_matrix()[1]
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        return self.dense_matrix()._right_kernel_matrix()


#########################

# used for a function above
cdef mpq_t minus_one
mpq_init(minus_one)
mpq_set_si(minus_one, -1,1)
