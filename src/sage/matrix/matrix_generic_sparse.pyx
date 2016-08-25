r"""
Sparse Matrices over a general ring

EXAMPLES::

    sage: R.<x> = PolynomialRing(QQ)
    sage: M = MatrixSpace(QQ['x'],2,3,sparse=True); M
    Full MatrixSpace of 2 by 3 sparse matrices over Univariate Polynomial Ring in x over Rational Field
    sage: a = M(range(6)); a
    [0 1 2]
    [3 4 5]
    sage: b = M([x^n for n in range(6)]); b
    [  1   x x^2]
    [x^3 x^4 x^5]
    sage: a * b.transpose()
    [            2*x^2 + x           2*x^5 + x^4]
    [      5*x^2 + 4*x + 3 5*x^5 + 4*x^4 + 3*x^3]
    sage: pari(a)*pari(b.transpose())
    [2*x^2 + x, 2*x^5 + x^4; 5*x^2 + 4*x + 3, 5*x^5 + 4*x^4 + 3*x^3]
    sage: c = copy(b); c
    [  1   x x^2]
    [x^3 x^4 x^5]
    sage: c[0,0] = 5; c
    [  5   x x^2]
    [x^3 x^4 x^5]
    sage: b[0,0]
    1
    sage: c.dict()
    {(0, 0): 5, (0, 1): x, (0, 2): x^2, (1, 0): x^3, (1, 1): x^4, (1, 2): x^5}
    sage: c.list()
    [5, x, x^2, x^3, x^4, x^5]
    sage: c.rows()
    [(5, x, x^2), (x^3, x^4, x^5)]
    sage: TestSuite(c).run()
    sage: d = c.change_ring(CC['x']); d
    [5.00000000000000                x              x^2]
    [             x^3              x^4              x^5]
    sage: latex(c)
    \left(\begin{array}{rrr}
    5 & x & x^{2} \\
    x^{3} & x^{4} & x^{5}
    \end{array}\right)
    sage: c.sparse_rows()
    [(5, x, x^2), (x^3, x^4, x^5)]
    sage: d = c.dense_matrix(); d
    [  5   x x^2]
    [x^3 x^4 x^5]
    sage: parent(d)
    Full MatrixSpace of 2 by 3 dense matrices over Univariate Polynomial Ring in x over Rational Field
    sage: c.sparse_matrix() is c
    True
    sage: c.is_sparse()
    True
"""

cimport matrix
cimport matrix_sparse
cimport sage.structure.element
from sage.structure.element cimport ModuleElement

import sage.misc.misc as misc

cdef class Matrix_generic_sparse(matrix_sparse.Matrix_sparse):
    r"""
    Generic sparse matrix.

    The ``Matrix_generic_sparse`` class derives from
    :class:`~sage.matrix.matrix_sparse.Matrix_sparse`, and defines functionality
    for sparse matrices over any base ring. A generic sparse matrix is
    represented using a dictionary whose keys are pairs of integers `(i,j)` and
    values in the base ring. The values of the dictionary must never be zero.

    EXAMPLES::

        sage: R.<a,b> = PolynomialRing(ZZ,'a,b')
        sage: M = MatrixSpace(R,5,5,sparse=True)
        sage: M({(0,0):5*a+2*b, (3,4): -a})
        [5*a + 2*b         0         0         0         0]
        [        0         0         0         0         0]
        [        0         0         0         0         0]
        [        0         0         0         0        -a]
        [        0         0         0         0         0]
        sage: M(3)
        [3 0 0 0 0]
        [0 3 0 0 0]
        [0 0 3 0 0]
        [0 0 0 3 0]
        [0 0 0 0 3]
        sage: V = FreeModule(ZZ, 5,sparse=True)
        sage: m = M([V({0:3}), V({2:2, 4:-1}), V(0), V(0), V({1:2})])
        sage: m
        [ 3  0  0  0  0]
        [ 0  0  2  0 -1]
        [ 0  0  0  0  0]
        [ 0  0  0  0  0]
        [ 0  2  0  0  0]

    .. NOTE::

        The datastructure can potentially be optimized. Firstly, as noticed in
        :trac:`17663`, we lose time in using 2-tuples to store indices.
        Secondly, there is no fast way to access non-zero elements in a given
        row/column.
    """
    ########################################################################
    # LEVEL 1 functionality
    #   * __cinit__
    #   * __init__
    #   * __dealloc__
    #   * set_unsafe
    #   * get_unsafe
    #   * def _pickle
    #   * def _unpickle
    ########################################################################
    def __cinit__(self, parent, entries=0, coerce=True, copy=True):
        self._entries = {}  # crucial so that pickling works

    def __init__(self, parent, entries=None, coerce=True, copy=True):
        r"""
        Create a sparse matrix over the given base ring.

        INPUT:

        - ``parent`` -- a matrix space

        - ``entries`` -- can be one of the following:

          * a Python dictionary whose items have the
            form ``(i, j): x``, where ``0 <= i < nrows``,
            ``0 <= j < ncols``, and ``x`` is coercible to
            an element of the base ring.
            The ``i,j`` entry of ``self`` is
            set to ``x``.  The ``x``'s can be ``0``.
          * Alternatively, entries can be a list of *all*
            the entries of the sparse matrix, read
            row-by-row from top to bottom (so they would
            be mostly 0).

        - ``coerce`` (default: ``True``) -- whether the entries
          should be coerced into the base ring before being
          entered into the matrix

        - ``copy`` (default: ``True``) -- whether the list or
          dictionary ``entries`` (not the single entries
          themselves!) should be copied before being
          entered into the matrix

        TESTS::

            sage: R.<a> = PolynomialRing(ZZ,'a')
            sage: M = MatrixSpace(R,2,3,sparse=True)
            sage: m = M([4,1,0,0,0,2]); m
            [4 1 0]
            [0 0 2]
            sage: m2 = copy(m)
            sage: m2[0,0] = -1
            sage: m[0,0]
            4
            sage: loads(dumps(m)) == m
            True

            sage: R2.<a,b> = PolynomialRing(QQ,'a','b')
            sage: M2 = MatrixSpace(R2,2,3,sparse=True)
            sage: M2(m)
            [4 1 0]
            [0 0 2]
            sage: M2.has_coerce_map_from(M)
            True

            sage: M3 = M2.change_ring(R2.fraction_field())
            sage: M3.has_coerce_map_from(M2)
            True


        Check that it is not possible to use wrong indices::

            sage: M = MatrixSpace(R,2,2,sparse=True)
            sage: M({(3,0): 1})
            Traceback (most recent call last):
            ...
            IndexError: matrix indices (3, 0) out of range

            sage: M({(0,-3): 1})
            Traceback (most recent call last):
            ...
            IndexError: matrix indices (0, -3) out of range

        But negative indices are valid::

            sage: M({(-1,-1): 1})
            [0 0]
            [0 1]
        """
        matrix.Matrix.__init__(self, parent)
        R = self._base_ring
        self._zero = R.zero()

        if entries is None or not entries:
            # be careful here. We might get entries set to be an empty list
            # because of the code implemented in matrix_space.MatrixSpace
            # So the condtion
            #   if entries is None or not entries:
            #       ...
            # is valid. But
            #   if entries is None or entries == 0:
            #       ...
            # is not!
            return

        cdef Py_ssize_t i, j, k

        if not isinstance(entries, dict):
            # assume that entries is a scalar
            x = R(entries)
            entries = {}
            if self._nrows != self._ncols:
                raise TypeError("scalar matrix must be square")
            for i from 0 <= i < self._nrows:
                entries[(i,i)] = x

        if coerce:
            v = {}
            for key, x in entries.iteritems():
                i,j = key
                if i < 0: i += self._nrows
                if j < 0: j += self._ncols
                if (i < 0 or i >= self._nrows or j < 0 or j >= self._ncols):
                    raise IndexError("matrix indices {} out of range".format(key))
                w = R(x)
                if w:
                    v[(i,j)] = w
            entries = v
        else:
            # Here we do not pay attention to the indices. We just check that it
            # *converts* to a pair of Py_ssize_t. In particular it is possible
            # to do:
            #
            #    sage: R = QQ['a','b']
            #    sage: M = MatrixSpace(R, 3, 3, sparse=True)
            #    sage: m = M({(Zmod(3)(1), Zmod(6)(2)): R.one()}, coerce=False)
            #
            #  and this is bad since:
            #
            #    sage: map(type,m.dict().keys()[0])
            #    [<type 'sage.rings.finite_rings.integer_mod.IntegerMod_int'>,
            #     <type 'sage.rings.finite_rings.integer_mod.IntegerMod_int'>]
            #
            # But not that setting coerce=False is advanced usage and we assume
            # that in such case the user knows what he/she is doing.
            if copy:
                entries = entries.copy()
            for key in entries.keys():
                i,j = key
                if i < 0: i += self._nrows
                if j < 0: j += self._ncols
                if (i < 0 or i >= self._nrows or j < 0 or j >= self._ncols):
                    raise IndexError("matrix indices {} out of range".format(key))
                if not entries[key]:
                    del entries[key]

        self._entries = entries

    def __nonzero__(self):
        r"""
        Test wether this matrix is non-zero.

        TESTS::

            sage: R.<a,b> = Zmod(5)['a','b']
            sage: m = matrix(R,3,4,sparse=True)
            sage: bool(m)      # indirect doctest
            False
            sage: m[0,3] = 1
            sage: bool(m)      # indirect doctest
            True
            sage: m[0,3] = 0   # indirect doctest
            sage: bool(m)
            False
            sage: m.is_zero()  # indirect doctest
            True
        """
        return bool(self._entries)

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, value):
        if not value:
            try:
                del self._entries[(i,j)]
            except KeyError:
                pass
        else:
            self._entries[(i,j)] = value

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        return self._entries.get((i,j), self._zero)

    def _pickle(self):
        version = 0
        return self._entries, version

    def _unpickle(self, data, int version):
        """
        EXAMPLES::

            sage: a = matrix([[1,10],[3,4]],sparse=True); a
            [ 1 10]
            [ 3  4]
            sage: loads(dumps(a)) == a
            True
        """
        if version == 0:
            self._entries = data
            self._zero = self._base_ring(0)
        else:
            raise RuntimeError("unknown matrix version (=%s)"%version)

    def __hash__(self):
        return self._hash()

    ########################################################################
    # LEVEL 2 functionality
    # x  * cdef _add_
    #    * cdef _mul_
    #    * cpdef _cmp_
    #    * __neg__
    #    * __invert__
    # x  * __copy__
    #    * _multiply_classical
    # x  * _list -- copy of the list of underlying elements
    # x  * _dict -- copy of the sparse dictionary of underlying elements
    ########################################################################

    cpdef _add_(self, _other):
        """
        EXAMPLES::

            sage: a = matrix([[1,10],[3,4]],sparse=True); a
            [ 1 10]
            [ 3  4]
            sage: a+a
            [ 2 20]
            [ 6  8]

        ::

            sage: a = matrix([[1,10,-5/3],[2/8,3,4]],sparse=True); a
            [   1   10 -5/3]
            [ 1/4    3    4]
            sage: a+a
            [    2    20 -10/3]
            [  1/2     6     8]
        """
        # Compute the sum of two sparse matrices.
        # This is complicated because of how we represent sparse matrices.
        # Algorithm:
        #   1. Sort both entry coefficient lists.
        #   2. March through building a new list, adding when the two i,j are the same.
        cdef Py_ssize_t i, j, len_v, len_w
        cdef Matrix_generic_sparse other
        other = <Matrix_generic_sparse> _other
        cdef list v = self._entries.items()
        v.sort()
        cdef list w = other._entries.items()
        w.sort()
        s = {}
        i = 0  # pointer into self
        j = 0  # pointer into other
        len_v = len(v)
        len_w = len(w)
        while i < len_v and j < len_w:
            vi = v[i][0]
            wj = w[j][0]
            if vi < wj:
                s[vi] = v[i][1]
                i += 1
            elif vi > wj:
                s[wj] = w[j][1]
                j += 1
            else:  # equal
                sm = v[i][1] + w[j][1]
                if sm:
                    s[vi] = sm
                i += 1
                j += 1
        while i < len(v):
            s[v[i][0]] = v[i][1]
            i += 1
        while j < len(w):
            s[w[j][0]] = w[j][1]
            j += 1

        cdef Matrix_generic_sparse A
        A = Matrix_generic_sparse.__new__(Matrix_generic_sparse, self._parent, 0,0,0)
        matrix.Matrix.__init__(A, self._parent)
        A._entries = s
        A._zero = self._zero
        A._base_ring = self._base_ring
        return A

    def __copy__(self):
        A = self.__class__(self._parent, self._entries, copy = True, coerce=False)
        if self._subdivisions is not None:
            A.subdivide(*self.subdivisions())
        return A


    def _list(self):
        """
        Return all entries of self as a list of numbers of rows times
        number of columns entries.
        """
        cdef Py_ssize_t i,j
        cdef list v = self.fetch('list')
        if v is None:
            v = [self._zero]*(self._nrows * self._ncols)
            for (i,j), x in self._entries.iteritems():
                v[i*self._ncols + j] = x
            self.cache('list', v)
        return v

    def _dict(self):
        """
        Return the underlying dictionary of self.

        This is used in comparisons.

        TESTS::

            sage: R.<a,b> = Zmod(6)['a','b']
            sage: M = MatrixSpace(R, 4,3)
            sage: m = M({(0,3): a+3*b*a, (1,1): -b})
            sage: m == m    # indirect doctest
            True
            sage: M(0) == m # indirect doctest
            False
        """
        return self._entries

    ########################################################################
    # LEVEL 3 functionality -- matrix windows, etc.
    ########################################################################

    def _nonzero_positions_by_row(self, copy=True):
        r"""
        TESTS::

            sage: R.<a> = PolynomialRing(Zmod(8), 'a')
            sage: M = MatrixSpace(R,4,3,sparse=True)
            sage: m = M({(3,0): 1, (3,1): 2*a^2 + 1, (2,0): -1, (0,1): -2})
            sage: m._nonzero_positions_by_row()
            [(0, 1), (2, 0), (3, 0), (3, 1)]
        """
        cdef list v = self.fetch('nonzero_positions')
        if v is None:
            v = self._entries.keys()
            v.sort()
            self.cache('nonzero_positions', v)
        if copy:
            return v[:]
        return v

    def _nonzero_positions_by_column(self, copy=True):
        r"""
        TESTS::

            sage: R.<a> = PolynomialRing(Zmod(8), 'a')
            sage: M = MatrixSpace(R,4,3,sparse=True)
            sage: m = M({(3,0): 1, (3,1): 2*a^2 + 1, (2,0): -1, (0,1): -2})
            sage: m._nonzero_positions_by_column()
            [(2, 0), (3, 0), (0, 1), (3, 1)]
        """
        cdef list v = self.fetch('nonzero_positions_by_column')
        if v is None:
            v = self._entries.keys()
            v.sort(_cmp_backward)
            self.cache('nonzero_positions_by_column', v)
        if copy:
            return v[:]
        return v

    ######################
    # Echelon support
    ######################


##     def dense_matrix(self):
##         import sage.matrix.matrix
##         M = sage.matrix.matrix.MatrixSpace(self.base_ring(),
##                                            self._nrows, self._ncols,
##                                            sparse = False)(0)
##         for i, j, x in self._entries
##             M[i,j] = x
##         return M

##     def nonpivots(self):
##         # We convert the pivots to a set so we have fast
##         # inclusion testing
##         X = set(self.pivots())
##         # [j for j in xrange(self.ncols()) if not (j in X)]
##         np = []
##         for j in xrange(self.ncols()):
##             if not (j in X):
##                 np.append(j)
##         return np

##     def matrix_from_nonpivot_columns(self):
##         """
##         The sparse matrix got by deleted all pivot columns.
##         """
##         return self.matrix_from_columns(self.nonpivots())

##     def matrix_from_columns(self, columns):
##         """
##         Returns the sparse submatrix of self composed of the given
##         list of columns.

##         INPUT:
##             columns -- list of int's
##         OUTPUT:
##             a sparse matrix.
##         """
##         # Let A be this matrix and let B denote this matrix with
##         # the columns deleted.
##         # ALGORITHM:
##         # 1. Make a table that encodes the function
##         #          f : Z --> Z,
##         #    f(j) = column of B where the j-th column of A goes
##         # 2. Build new list of entries and return resulting matrix.
##         C = set(columns)
##         X = []
##         j = 0
##         for i in xrange(self.ncols()):
##             if i in C:
##                 X.append(j)
##                 j = j + 1
##             else:
##                 X.append(-1)    # column to be deleted.
##         entries2 = []
##         for i, j, x in self.entries():
##             if j in C:
##                 entries2.append((i,X[j],x))
##         return SparseMatrix(self.base_ring(), self.nrows(),
##                             len(C), entries2, coerce=False, sort=False)

##     def matrix_from_rows(self, rows):
##         """
##         Returns the sparse submatrix of self composed of the given
##         list of rows.

##         INPUT:
##             rows -- list of int's
##         OUTPUT:
##             a sparse matrix.
##         """
##         R = set(rows)
##         if not R.issubset(set(xrange(self.nrows()))):
##             raise ArithmeticError("Invalid rows.")
##         X = []
##         i = 0
##         for j in xrange(self.nrows()):
##             if j in R:
##                 X.append(i)
##                 i = i + 1
##             else:
##                 X.append(-1)    # row to be deleted.
##         entries2 = []
##         for i, j, x in self.entries():
##             if i in R:
##                 entries2.append((X[i],j,x))
##         return SparseMatrix(self.base_ring(), len(R),
##                             self.ncols(), entries2, coerce=False, sort=False)


##     def transpose(self):
##         """
##         Returns the transpose of self.
##         """
##         entries2 = [] # [(j,i,x) for i,j,x in self.entries()]
##         for i,j,x in self.entries():
##             entries2.append((j,i,x))
##         return SparseMatrix(self.base_ring(), self.ncols(),
##                             self.nrows(), entries2, coerce=False, sort=True)

##     def __rmul__(self, left):
##         return self.scalar_multiple(left)

##     def scalar_multiple(self, left):
##         R = self.base_ring()
##         left = R(left)
##         if left == R(1):
##             return self
##         if left == R(0):
##             return SparseMatrix(R, self.nrows(), self.ncols(), coerce=False, sort=False)

##         X = []
##         for i, j, x in self.list():
##             X.append((i,j,left*x))
##         return SparseMatrix(self.base_ring(), self.nrows(),
##                             self.ncols(), X, coerce=False, sort=False)


####################################################################################
# Various helper functions
####################################################################################
import matrix_space

def Matrix_sparse_from_rows(X):
    """
    INPUT:


    -  ``X`` - nonempty list of SparseVector rows


    OUTPUT: Sparse_matrix with those rows.

    EXAMPLES::

        sage: V = VectorSpace(QQ,20,sparse=True)
        sage: v = V(0)
        sage: v[9] = 4
        sage: from sage.matrix.matrix_generic_sparse import Matrix_sparse_from_rows
        sage: Matrix_sparse_from_rows([v])
        [0 0 0 0 0 0 0 0 0 4 0 0 0 0 0 0 0 0 0 0]
        sage: Matrix_sparse_from_rows([v, v, v, V(0)])
        [0 0 0 0 0 0 0 0 0 4 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 4 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 4 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
    """
    cdef Py_ssize_t i, j

    if not isinstance(X, (list, tuple)):
        raise TypeError("X (=%s) must be a list or tuple"%X)
    if len(X) == 0:
        raise ArithmeticError("X must be nonempty")
    entries = {}
    R = X[0].base_ring()
    ncols = X[0].degree()
    for i from 0 <= i < len(X):
        for j, x in X[i].iteritems():
            entries[(i,j)] = x
    M = matrix_space.MatrixSpace(R, len(X), ncols, sparse=True)
    return M(entries, coerce=False, copy=False)

def _cmp_backward(x, y):  # todo: speed up via Python/C API
    r"""
    TESTS::

        sage: from sage.matrix.matrix_generic_sparse import _cmp_backward
        sage: l0 = [(-1,-1), (0,0), (1,0), (-1,1), (0,1), (1,1), (-1,2)]
        sage: l = l0[:]
        sage: for _ in range(10):
        ....:   shuffle(l)
        ....:   l.sort(_cmp_backward)
        ....:   assert l0 == l
    """
    # compare two 2-tuples, but in reverse order, i.e., second entry than first
    cdef Py_ssize_t i,j
    i = x[1]
    j = y[1]
    if i < j:
        return -1
    elif i > j:
        return 1
    i = x[0]
    j = y[0]
    return i-j
