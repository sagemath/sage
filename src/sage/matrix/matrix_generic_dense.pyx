# sage_setup: distribution = sagemath-modules
"""
Dense Matrices over a general ring
"""

cimport cython
from cpython.list cimport *
from cpython.number cimport *
from cpython.ref cimport *

cimport sage.matrix.matrix_dense as matrix_dense
from sage.matrix import matrix_dense
from sage.matrix.args cimport MatrixArgs_init

cimport sage.matrix.matrix as matrix


cdef class Matrix_generic_dense(matrix_dense.Matrix_dense):
    r"""
    The ``Matrix_generic_dense`` class derives from
    ``Matrix``, and defines functionality for dense
    matrices over any base ring. Matrices are represented by a list of
    elements in the base ring, and element access operations are
    implemented in this class.

    EXAMPLES::

        sage: A = random_matrix(Integers(25)['x'], 2)
        sage: type(A)
        <class 'sage.matrix.matrix_generic_dense.Matrix_generic_dense'>
        sage: TestSuite(A).run(skip='_test_minpoly')

    Test comparisons::

        sage: A = random_matrix(Integers(25)['x'], 2)
        sage: A == A
        True
        sage: A < A + 1 or A[0, 0].coefficients()[0] == 24
        True
        sage: A+1 < A and A[0, 0].coefficients()[0] != 24
        False

    Test hashing::

        sage: A = random_matrix(Integers(25)['x'], 2)
        sage: hash(A)
        Traceback (most recent call last):
        ...
        TypeError: mutable matrices are unhashable
        sage: A.set_immutable()
        sage: H = hash(A)
    """
    def __init__(self, parent, entries=None, copy=None, bint coerce=True):
        r"""
        Initialize a dense matrix.

        INPUT:

        - ``parent`` -- a matrix space

        - ``entries`` -- see :func:`matrix`

        - ``copy`` -- ignored (for backwards compatibility)

        - ``coerce`` -- if ``False``, assume without checking that the
          entries lie in the base ring

        TESTS:

        We check that the problem related to :issue:`9049` is not an issue any
        more::

            sage: # needs sage.rings.number_field
            sage: S.<t> = PolynomialRing(QQ)
            sage: F.<q> = QQ.extension(t^4 + 1)
            sage: R.<x,y> = PolynomialRing(F)
            sage: M = MatrixSpace(R, 1, 2)
            sage: from sage.matrix.matrix_generic_dense import Matrix_generic_dense
            sage: Matrix_generic_dense(M, (x, y), True, True)
            [x y]
        """
        ma = MatrixArgs_init(parent, entries)
        self._entries = ma.list(coerce)

    cdef Matrix_generic_dense _new(self, Py_ssize_t nrows, Py_ssize_t ncols):
        r"""
        Return a new dense matrix with no entries set.
        """
        if nrows == self._nrows and ncols == self._ncols:
            MS = self._parent
        else:
            MS = self.matrix_space(nrows, ncols)

        cdef type t = <type>type(self)
        return <Matrix_generic_dense>t.__new__(t, MS)

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, value):
        self._entries[i*self._ncols + j] = value

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        return self._entries[i*self._ncols + j]

    def _reverse_unsafe(self):
        r"""
        TESTS::

            sage: m = matrix(ZZ['x,y'], 2, 3, range(6))
            sage: m._reverse_unsafe()
            sage: m
            [5 4 3]
            [2 1 0]
        """
        self._entries.reverse()

    def _pickle(self):
        """
        EXAMPLES::

            sage: R.<x> = Integers(25)['x']; A = matrix(R, [1,x,x^3+1,2*x])
            sage: A._pickle()
            ([1, x, x^3 + 1, 2*x], 0)
        """
        return self._entries, 0

    def _unpickle(self, data, int version):
        """
        EXAMPLES::

            sage: R.<x> = Integers(25)['x']; A = matrix(R, [1,x,x^3+1,2*x]); B = A.parent()(0)
            sage: v = A._pickle()
            sage: B._unpickle(v[0], v[1])
            sage: B
            [      1       x x^3 + 1     2*x]
        """
        if version == 0:
            self._entries = data
        else:
            raise RuntimeError("unknown matrix version")

    ########################################################################
    # LEVEL 2 functionality
    # X  * cdef _add_
    #    * cdef _mul_
    #    * cpdef _richcmp_
    #    * __neg__
    #    * __invert__
    # x  * __copy__
    # x  * _multiply_classical
    # x  * _list -- copy of the list of underlying elements
    #    * _dict -- copy of the sparse dictionary of underlying elements
    ########################################################################

    def __copy__(self):
        """
        Create a copy of self, which may be changed without altering
        ``self``.

        EXAMPLES::

            sage: A = matrix(ZZ[['t']], 2,3,range(6)); A
            [0 1 2]
            [3 4 5]
            sage: A.subdivide(1,1); A
            [0|1 2]
            [-+---]
            [3|4 5]
            sage: B = A.__copy__(); B
            [0|1 2]
            [-+---]
            [3|4 5]
            sage: B == A
            True
            sage: B[0,0] = 100
            sage: B
            [100|  1   2]
            [---+-------]
            [  3|  4   5]
            sage: A
            [0|1 2]
            [-+---]
            [3|4 5]
            sage: R.<x> = QQ['x']
            sage: a = matrix(R,2,[x+1,2/3,  x^2/2, 1+x^3]); a
            [  x + 1     2/3]
            [1/2*x^2 x^3 + 1]
            sage: b = copy(a)
            sage: b[0,0] = 5
            sage: b
            [      5     2/3]
            [1/2*x^2 x^3 + 1]
            sage: a
            [  x + 1     2/3]
            [1/2*x^2 x^3 + 1]

        ::

            sage: b = copy(a)
            sage: f = b[0,0]; f[0] = 10
            Traceback (most recent call last):
            ...
            IndexError: polynomials are immutable
        """
        cdef Matrix_generic_dense A
        A = self._new(self._nrows, self._ncols)
        A._entries = self._entries[:]
        if self._subdivisions is not None:
            A.subdivide(*self.subdivisions())
        return A

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef _add_(self, right):
        """
        Add two generic dense matrices with the same parent.

        EXAMPLES::

            sage: # needs sage.combinat
            sage: R.<x,y> = FreeAlgebra(QQ, 2)
            sage: a = matrix(R, 2, 2, [1,2,x*y,y*x])
            sage: b = matrix(R, 2, 2, [1,2,y*x,y*x])
            sage: a._add_(b)
            [        2         4]
            [x*y + y*x     2*y*x]
        """
        cdef Py_ssize_t k
        cdef Matrix_generic_dense other = <Matrix_generic_dense> right
        cdef Matrix_generic_dense res = self._new(self._nrows, self._ncols)
        res._entries = [None]*(self._nrows*self._ncols)
        for k in range(self._nrows*self._ncols):
            res._entries[k] = self._entries[k] + other._entries[k]
        return res

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef _sub_(self, right):
        """
        Subtract two generic dense matrices with the same parent.

        EXAMPLES::

            sage: # needs sage.combinat
            sage: R.<x,y> = FreeAlgebra(QQ, 2)
            sage: a = matrix(R, 2, 2, [1,2,x*y,y*x])
            sage: b = matrix(R, 2, 2, [1,2,y*x,y*x])
            sage: a._sub_(b)
            [        0         0]
            [x*y - y*x         0]
        """
        cdef Py_ssize_t k
        cdef Matrix_generic_dense other = <Matrix_generic_dense> right
        cdef Matrix_generic_dense res = self._new(self._nrows, self._ncols)
        res._entries = [None]*(self._nrows*self._ncols)
        for k in range(self._nrows*self._ncols):
            res._entries[k] = self._entries[k] - other._entries[k]
        return res

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.overflowcheck(False)
    def _multiply_classical(left, matrix.Matrix _right):
        """
        Multiply the matrices left and right using the classical
        `O(n^3)` algorithm.

        EXAMPLES:

        We multiply two matrices over a fairly general ring::

            sage: R.<x,y> = Integers(8)['x,y']
            sage: a = matrix(R,2,[x,y,x^2,y^2]); a
            [  x   y]
            [x^2 y^2]
            sage: type(a)
            <class 'sage.matrix.matrix_generic_dense.Matrix_generic_dense'>
            sage: a*a
            [  x^2*y + x^2     y^3 + x*y]
            [x^2*y^2 + x^3   y^4 + x^2*y]
            sage: a.det()^2 == (a*a).det()
            True
            sage: a._multiply_classical(a)
            [  x^2*y + x^2     y^3 + x*y]
            [x^2*y^2 + x^3   y^4 + x^2*y]

            sage: A = matrix(QQ['x,y'], 2, [0,-1,2,-2])
            sage: B = matrix(QQ['x,y'], 2, [-1,-1,-2,-2])
            sage: A*B
            [2 2]
            [2 2]

        Sage fully supports degenerate matrices with 0 rows or 0 columns::

            sage: A = matrix(QQ['x,y'], 0, 4, []); A
            []
            sage: B = matrix(QQ['x,y'], 4,0, []); B
            []
            sage: A*B
            []
            sage: B*A
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]
        """
        cdef Py_ssize_t i, j, k, m, nr, nc, snc, p
        cdef Matrix_generic_dense right = _right

        if left._ncols != right._nrows:
            raise IndexError("Number of columns of left must equal number of rows of other.")

        nr = left._nrows
        nc = right._ncols
        snc = left._ncols

        R = left.base_ring()
        cdef list v = [None] * (left._nrows * right._ncols)
        zero = R.zero()
        p = 0
        for i in range(nr):
            for j in range(nc):
                z = zero
                m = i*snc
                for k in range(snc):
                    z += left._entries[m+k]._mul_(right._entries[k*nc+j])
                v[p] = z
                p += 1

        cdef Matrix_generic_dense A = left._new(nr, nc)
        A._entries = v
        return A

    def _list(self):
        """
        Return reference to list of entries of ``self``.  For internal use
        only, since this circumvents immutability.

        EXAMPLES::

            sage: A = random_matrix(Integers(25)['x'],2); A.set_immutable()
            sage: A._list()[0] = 0
            sage: A._list()[0]
            0
        """
        return self._entries

    ########################################################################
    # LEVEL 3 functionality (Optional)
    #    * cdef _sub_
    #    * __deepcopy__
    #    * __invert__
    #    * _multiply_classical
    #    * Matrix windows -- only if you need strassen for that base
    #    * Other functions (list them here):
    ########################################################################
