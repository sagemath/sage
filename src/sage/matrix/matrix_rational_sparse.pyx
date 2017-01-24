"""
Sparse rational matrices.

AUTHORS:

- William Stein (2007-02-21)
- Soroosh Yazdani (2007-02-21)

TESTS::

    sage: a = matrix(QQ,2,range(4), sparse=True)
    sage: TestSuite(a).run()
    sage: matrix(QQ,0,0,sparse=True).inverse()
    []
"""

##############################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

from sage.data_structures.binary_search cimport *
from sage.modules.vector_integer_sparse cimport *
from sage.modules.vector_rational_sparse cimport *

include 'sage/ext/stdsage.pxi'
include "cysignals/signals.pxi"
from cpython.sequence cimport *

from sage.rings.rational cimport Rational
from sage.rings.integer  cimport Integer
from matrix cimport Matrix

from sage.libs.gmp.mpz cimport *
from sage.libs.gmp.mpq cimport *

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

cimport sage.structure.element

import sage.matrix.matrix_space

from matrix_integer_sparse cimport Matrix_integer_sparse
from matrix_rational_dense cimport Matrix_rational_dense

from sage.misc.misc import verbose

cdef class Matrix_rational_sparse(Matrix_sparse):

    ########################################################################
    # LEVEL 1 functionality
    #   * __cinit__
    #   * __dealloc__
    #   * __init__
    #   * set_unsafe
    #   * get_unsafe
    #   * __hash__       -- always simple
    ########################################################################
    def __cinit__(self, parent, entries, copy, coerce):
        # set the parent, nrows, ncols, etc.
        Matrix_sparse.__init__(self, parent)

        self._matrix = <mpq_vector*> sig_malloc(parent.nrows()*sizeof(mpq_vector))
        if self._matrix == NULL:
            raise MemoryError("error allocating sparse matrix")
        # initialize the rows
        for i from 0 <= i < parent.nrows():
            mpq_vector_init(&self._matrix[i], self._ncols, 0)

        # record that rows have been initialized
        self._initialized = True


    def __dealloc__(self):
        self._dealloc()

    cdef _dealloc(self):
        cdef Py_ssize_t i
        if self._initialized:
            for i from 0 <= i < self._nrows:
                mpq_vector_clear(&self._matrix[i])
        if self._matrix != NULL:
            sig_free(self._matrix)

    def __init__(self, parent, entries, copy, coerce):
        """
        Create a sparse matrix over the rational numbers.

        INPUT:

        - ``parent`` -- a matrix space

        - ``entries`` -- can be one of the following:

          * a Python dictionary whose items have the
            form ``(i, j): x``, where ``0 <= i < nrows``,
            ``0 <= j < ncols``, and ``x`` is coercible to
            a rational.  The ``i,j`` entry of ``self`` is
            set to ``x``.  The ``x``'s can be ``0``.
          * Alternatively, entries can be a list of *all*
            the entries of the sparse matrix, read
            row-by-row from top to bottom (so they would
            be mostly 0).

        - ``copy`` -- ignored

        - ``coerce`` -- ignored
        """
        cdef Py_ssize_t i, j, k
        cdef Rational z
        cdef PyObject** X

        if entries is None: return
        # fill in entries in the dict case
        if isinstance(entries, dict):
            R = self.base_ring()
            for ij, x in entries.iteritems():
                z = R(x)
                if z != 0:
                    i, j = ij  # nothing better to do since this is user input, which may be bogus.
                    if i < 0 or j < 0 or i >= self._nrows or j >= self._ncols:
                        raise IndexError("invalid entries list")
                    mpq_vector_set_entry(&self._matrix[i], j, z.value)
        elif isinstance(entries, list):
            # Dense input format -- fill in entries
            if len(entries) != self._nrows * self._ncols:
                raise TypeError("list of entries must be a dictionary of (i,j):x or a dense list of n * m elements")
            seq = PySequence_Fast(entries,"expected a list")
            X = PySequence_Fast_ITEMS(seq)
            k = 0
            R = self.base_ring()
            # Get fast access to the entries list.
            for i from 0 <= i < self._nrows:
                for  j from 0 <= j < self._ncols:
                    z = R(<object>X[k])
                    if z != 0:
                        mpq_vector_set_entry(&self._matrix[i], j, z.value)
                    k = k + 1
        else:
            # fill in entries in the scalar case
            z = Rational(entries)
            if z == 0:
                return
            if self._nrows != self._ncols:
                raise TypeError("matrix must be square to initialize with a scalar.")
            for i from 0 <= i < self._nrows:
                mpq_vector_set_entry(&self._matrix[i], i, z.value)


    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, x):
        mpq_vector_set_entry(&self._matrix[i], j, (<Rational> x).value)

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        cdef Rational x
        x = Rational()
        mpq_vector_get_entry(x.value, &self._matrix[i], j)
        return x

    def __hash__(self):
        return self._hash()


    ########################################################################
    # LEVEL 2 functionality
    #   * def _pickle
    #   * def _unpickle
    #   * cdef _add_
    #   * cdef _sub_
    #   * cdef _mul_
    #   * cpdef _cmp_
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
        nonzero_positions_in_columns = [set([]) for _ in range(right._ncols)]
        cdef Py_ssize_t i, j, k
        for i from 0 <= i < right._nrows:
            v = &(right._matrix[i])
            for j from 0 <= j < right._matrix[i].num_nonzero:
                nonzero_positions_in_columns[v.positions[j]].add(i)

        ans = self.new_matrix(self._nrows, right._ncols)

        # Now do the multiplication, getting each row completely before filling it in.
        cdef mpq_t x, y, s
        mpq_init(x)
        mpq_init(y)
        mpq_init(s)
        for i from 0 <= i < self._nrows:
            v = &self._matrix[i]
            for j from 0 <= j < right._ncols:
                mpq_set_si(s, 0, 1)
                c = nonzero_positions_in_columns[j]
                for k from 0 <= k < v.num_nonzero:
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
            <type 'sage.matrix.matrix_rational_dense.Matrix_rational_dense'>
        """
        cdef Matrix_rational_sparse right
        cdef Matrix_rational_dense ans
        right = _right

        cdef mpq_vector* v

        # Build a table that gives the nonzero positions in each column of right
        nonzero_positions_in_columns = [set([]) for _ in range(right._ncols)]
        cdef Py_ssize_t i, j, k
        for i from 0 <= i < right._nrows:
            v = &(right._matrix[i])
            for j from 0 <= j < right._matrix[i].num_nonzero:
                nonzero_positions_in_columns[v.positions[j]].add(i)

        ans = self.new_matrix(self._nrows, right._ncols, sparse=False)

        # Now do the multiplication, getting each row completely before filling it in.
        cdef mpq_t x, y, s
        mpq_init(x)
        mpq_init(y)
        mpq_init(s)
        for i from 0 <= i < self._nrows:
            v = &self._matrix[i]
            for j from 0 <= j < right._ncols:
                mpq_set_si(s, 0, 1)
                c = nonzero_positions_in_columns[j]
                for k from 0 <= k < v.num_nonzero:
                    if v.positions[k] in c:
                        mpq_vector_get_entry(y, &right._matrix[v.positions[k]], j)
                        mpq_mul(x, v.entries[k], y)
                        mpq_add(s, s, x)
                mpq_set(ans._matrix[i][j], s)

        mpq_clear(x)
        mpq_clear(y)
        mpq_clear(s)
        return ans


    ########################################################################
    # def _pickle(self):
    # def _unpickle(self, data, int version):   # use version >= 0
    # cpdef _add_(self, right):
    # cdef _mul_(self, Matrix right):
    # cpdef int _cmp_(self, Matrix right) except -2:
    # def __neg__(self):
    # def __invert__(self):
    # def __copy__(self):
    # def _multiply_classical(left, matrix.Matrix _right):
    # def _list(self):

# TODO
##     cpdef _lmul_(self, RingElement right):
##         """
##         EXAMPLES:
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
        if not d is None:
            return d

        cdef Py_ssize_t i, j, k
        d = {}
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._matrix[i].num_nonzero:
                x = Rational()
                mpq_set((<Rational>x).value, self._matrix[i].entries[j])
                d[(int(i),int(self._matrix[i].positions[j]))] = x
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
        Returns the list of pairs (i,j) such that self[i,j] != 0.

        It is safe to change the resulting list (unless you give the option copy=False).

        EXAMPLE::

            sage: M = Matrix(QQ, [[0,0,0,1,0,0,0,0],[0,1,0,0,0,0,1,0]], sparse=True); M
            [0 0 0 1 0 0 0 0]
            [0 1 0 0 0 0 1 0]
            sage: M.nonzero_positions()
            [(0, 3), (1, 1), (1, 6)]

        """
        x = self.fetch('nonzero_positions')
        if not x is None:
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

        OUTPUT:

            -- Integer

        EXAMPLES::

            sage: b = matrix(QQ,2,range(6), sparse=True); b[0,0]=-5007/293; b
            [-5007/293         1         2]
            [        3         4         5]
            sage: b.height()
            5007
        """
        cdef Integer z
        z = PY_NEW(Integer)
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
        mpz_set_si(d,1)
        cdef Py_ssize_t i, j
        cdef mpq_vector *v

        sig_on()
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._matrix[i].num_nonzero:
                mpz_lcm(d, d, mpq_denref(self._matrix[i].entries[j]))
        sig_off()
        return 0


    def denominator(self):
        """
        Return the denominator of this matrix.

        OUTPUT:

            -- Sage Integer

        EXAMPLES::

            sage: b = matrix(QQ,2,range(6)); b[0,0]=-5007/293; b
            [-5007/293         1         2]
            [        3         4         5]
            sage: b.denominator()
            293
        """
        cdef Integer z
        z = PY_NEW(Integer)
        self.mpz_denom(z.value)
        return z

    def _clear_denom(self):
        """
        INPUT:

        self -- a matrix

        OUTPUT:

        D*self, D

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
        A._initialized = 1
        return A, D


    ################################################
    # Echelon form
    ################################################
    def echelonize(self, height_guess=None, proof=True, **kwds):
        """
        Transform the matrix ``self`` into reduced row echelon form
        in place.

        INPUT:

        ``height_guess``, ``proof``, ``**kwds`` -- all passed to the multimodular
        algorithm; ignored by the p-adic algorithm.

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

        :trac:`10319` has been fixed::

            sage: m = Matrix(QQ, [1], sparse=True); m.echelonize()
            sage: m = Matrix(QQ, [1], sparse=True); m.echelonize(); m
            [1]
        """

        x = self.fetch('in_echelon_form')
        if not x is None: return  # already known to be in echelon form
        self.check_mutability()
        self.clear_cache()

        pivots = self._echelonize_multimodular(height_guess, proof, **kwds)

        self.cache('in_echelon_form', True)
        self.cache('pivots', tuple(pivots))


    def echelon_form(self, algorithm='default',
                     height_guess=None, proof=True, **kwds):
        """
        INPUT:

        ``height_guess``, ``proof``, ``**kwds`` -- all passed to the multimodular
        algorithm; ignored by the p-adic algorithm.

        OUTPUT:

            self is no in reduced row echelon form.

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
        if not x is None:
            return x
        if self.fetch('in_echelon_form'): return self

        E = self._echelon_form_multimodular(height_guess, proof=proof)

        self.cache(label, E)
        self.cache('pivots', E.pivots())
        return E

    # Multimodular echelonization algorithms
    def _echelonize_multimodular(self, height_guess=None, proof=True, **kwds):
        cdef Py_ssize_t i, j
        cdef Matrix_rational_sparse E
        cdef mpq_vector* v
        cdef mpq_vector* w
        E = self._echelon_form_multimodular(height_guess, proof=proof, **kwds)
        # Get rid of self's data
        self._dealloc()
        # Copy E's data to self's data.
        self._matrix = <mpq_vector*> sig_malloc(E._nrows * sizeof(mpq_vector))
        if self._matrix == NULL:
            raise MemoryError("error allocating sparse matrix")
        for i from 0 <= i < E._nrows:
            v = &self._matrix[i]
            w = &E._matrix[i]
            mpq_vector_init(v, E._ncols, w.num_nonzero)
            for j from 0 <= j < w.num_nonzero:
                mpq_set(v.entries[j], w.entries[j])
                v.positions[j] = w.positions[j]
        return E.pivots()


    def _echelon_form_multimodular(self, height_guess=None, proof=True):
        """
        Returns reduced row-echelon form using a multi-modular
        algorithm.  Does not change self.

        INPUT:

        - height_guess -- integer or None
        - proof -- boolean (default: True)
        """
        import misc
        cdef Matrix E = misc.matrix_rational_echelon_form_multimodular(self,
                                 height_guess=height_guess, proof=proof)
        E._parent = self._parent
        return E


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
            <type 'sage.matrix.matrix_rational_sparse.Matrix_rational_sparse'>
            sage: type(a.dense_matrix())
            <type 'sage.matrix.matrix_rational_dense.Matrix_rational_dense'>
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
                mpq_set(B._matrix[i][v.positions[j]], v.entries[j])
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
        Set row i of self to -(row r of A), but where we only take the
        given column positions in that row of A.  Note that we *DO*
        zero out the other entries of self's row i.

        INPUT:

        - i -- integer, index into the rows of self
        - A -- a sparse matrix
        - r -- integer, index into rows of A
        - cols -- a *sorted* list of integers.
        - cols_index -- (optional).  But set it to this to vastly speed up
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
        Returns a pair that includes a matrix of basis vectors
        for the right kernel of ``self``.

        INPUT:

        - ``kwds`` - these are provided for consistency with other versions
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

