r"""
Sparse Matrices over a general ring

EXAMPLES::

    sage: R.<x> = PolynomialRing(QQ)
    sage: M = MatrixSpace(QQ['x'], 2, 3, sparse=True); M
    Full MatrixSpace of 2 by 3 sparse matrices over
     Univariate Polynomial Ring in x over Rational Field
    sage: a = M(range(6)); a
    [0 1 2]
    [3 4 5]
    sage: b = M([x^n for n in range(6)]); b
    [  1   x x^2]
    [x^3 x^4 x^5]
    sage: a * b.transpose()
    [            2*x^2 + x           2*x^5 + x^4]
    [      5*x^2 + 4*x + 3 5*x^5 + 4*x^4 + 3*x^3]
    sage: pari(a)*pari(b.transpose())                                                   # needs sage.libs.pari
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
    Full MatrixSpace of 2 by 3 dense matrices
     over Univariate Polynomial Ring in x over Rational Field
    sage: c.sparse_matrix() is c
    True
    sage: c.is_sparse()
    True
"""
cimport sage.matrix.matrix_sparse as matrix_sparse
cimport sage.structure.element
from sage.matrix.args cimport MatrixArgs_init


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
        :issue:`17663`, we lose time in using 2-tuples to store indices.
        Secondly, there is no fast way to access nonzero elements in a given
        row/column.
    """
    def __cinit__(self):
        self._entries = {}  # crucial so that pickling works
        self._zero = self._base_ring.zero()

    def __init__(self, parent, entries=None, copy=None, bint coerce=True):
        r"""
        Create a sparse matrix over the given base ring.

        INPUT:

        - ``parent`` -- a matrix space

        - ``entries`` -- see :func:`matrix`

        - ``copy`` -- ignored (for backwards compatibility)

        - ``coerce`` -- if ``False``, assume without checking that the
          entries lie in the base ring

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

            sage: R2.<a,b> = PolynomialRing(QQ)
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
            IndexError: invalid row index 3
            sage: M({(0,-3): 1})
            Traceback (most recent call last):
            ...
            IndexError: invalid column index -3

        But negative indices are valid and wrap around::

            sage: M({(-1,-1): 1})
            [0 0]
            [0 1]
        """
        ma = MatrixArgs_init(parent, entries)
        self._entries = ma.dict(coerce)

    def __bool__(self):
        r"""
        Test whether this matrix is nonzero.

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

    cdef bint get_is_zero_unsafe(self, Py_ssize_t i, Py_ssize_t j) except -1:
        """
        Return 1 if the entry ``(i, j)`` is zero, otherwise 0.

        EXAMPLES::

            sage: R.<a,b> = Zmod(5)['a','b']
            sage: m = matrix(R,2,4, {(1,3): a, (0,0):b}, sparse=True)
            sage: m.zero_pattern_matrix()  # indirect doctest
            [0 1 1 1]
            [1 1 1 0]
        """
        return (i,j) not in self._entries

    def _pickle(self):
        version = 0
        return self._entries, version

    def _unpickle(self, data, int version):
        """
        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: a = matrix(R, [[1,10],[3,4]],sparse=True); a
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

    ########################################################################
    # LEVEL 2 functionality
    # x  * cdef _add_
    #    * cdef _mul_
    #    * cpdef _richcmp_
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

            sage: R.<x> = QQ[]
            sage: a = matrix(R, [[1,10],[3,4]],sparse=True); a
            [ 1 10]
            [ 3  4]
            sage: a+a
            [ 2 20]
            [ 6  8]

        ::

            sage: a = matrix(R, [[1,10,-5/3],[2/8,3,4]], sparse=True); a
            [   1   10 -5/3]
            [ 1/4    3    4]
            sage: a+a
            [    2    20 -10/3]
            [  1/2     6     8]

        TESTS:

        Adding two subclasses should produce an instance of the subclass; not
        a generic matrix::

            sage: A = matrix(RDF, [[1]], sparse=True)
            sage: B = matrix(RDF, [[2]], sparse=True)
            sage: C = matrix(ZZ,  [[3]], sparse=True)
            sage: D = matrix(CDF, [[4]], sparse=True)
            sage: (A+B).__class__ == A.__class__
            True
            sage: (B+A).__class__ == A.__class__
            True
            sage: (A+C).__class__ == A.__class__
            True
            sage: (C+A).__class__ == A.__class__
            True
            sage: (A+D).__class__ == D.__class__
            True
        """
        # Compute the sum of two sparse matrices.
        # This is complicated because of how we represent sparse matrices.
        # Algorithm:
        #   1. Sort both entry coefficient lists.
        #   2. March through building a new list, adding when the two i,j are the same.
        cdef Py_ssize_t i, j, len_v, len_w
        cdef Matrix_generic_sparse other
        other = <Matrix_generic_sparse> _other
        cdef list v = sorted(self._entries.items())
        cdef list w = sorted(other._entries.items())
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

        # Use the parent of the left-hand summand as the parent of the
        # result. This allows you to (for example) add two sparse RDF
        # matrices together and get another one back, rather than
        # getting a generic sparse matrix back.
        MS = self._parent
        cdef type t = <type>type(self)
        A = <Matrix_generic_sparse>t.__new__(t, MS)
        A._entries = s
        A._zero = self._zero
        return A

    def __copy__(self):
        A = self.__class__(self._parent, self._entries, copy = True, coerce=False)
        if self._subdivisions is not None:
            A.subdivide(*self.subdivisions())
        return A

    def _list(self):
        """
        Return all entries of ``self`` as a list of numbers of rows times
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
        Return the underlying dictionary of ``self``.

        This is used in comparisons.

        TESTS::

            sage: R.<a,b> = Zmod(6)[]
            sage: M = MatrixSpace(R, 3, 4)
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
            v = sorted(self._entries)
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
            v = sorted(self._entries, key=lambda x: (x[1], x[0]))
            self.cache('nonzero_positions_by_column', v)
        if copy:
            return v[:]
        return v


####################################################################################
# Various helper functions
####################################################################################

def Matrix_sparse_from_rows(X):
    """
    INPUT:

    - ``X`` -- nonempty list of ``SparseVector`` rows

    OUTPUT: ``Sparse_matrix`` with those rows

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
        raise TypeError("X (=%s) must be a list or tuple" % X)
    if not X:
        raise ArithmeticError("X must be nonempty")

    from sage.matrix import matrix_space
    entries = {}
    R = X[0].base_ring()
    ncols = X[0].degree()
    for i from 0 <= i < len(X):
        for j, x in X[i].iteritems():
            entries[(i,j)] = x
    M = matrix_space.MatrixSpace(R, len(X), ncols, sparse=True)
    return M(entries, coerce=False, copy=False)
