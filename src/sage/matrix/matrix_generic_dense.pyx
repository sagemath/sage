"""
Dense Matrices over a general ring
"""
cimport cython
from cpython.list cimport *
from cpython.number cimport *
from cpython.ref cimport *

cimport matrix_dense
import matrix_dense

cimport matrix

from sage.structure.element cimport parent_c

cdef class Matrix_generic_dense(matrix_dense.Matrix_dense):
    r"""
    The ``Matrix_generic_dense`` class derives from
    ``Matrix``, and defines functionality for dense
    matrices over any base ring. Matrices are represented by a list of
    elements in the base ring, and element access operations are
    implemented in this class.

    EXAMPLES::

        sage: A = random_matrix(Integers(25)['x'],2); A
        [    x^2 + 12*x + 2   4*x^2 + 13*x + 8]
        [ 22*x^2 + 2*x + 17 19*x^2 + 22*x + 14]
        sage: type(A)
        <type 'sage.matrix.matrix_generic_dense.Matrix_generic_dense'>
        sage: TestSuite(A).run()

    Test comparisons::

        sage: A = random_matrix(Integers(25)['x'],2)
        sage: cmp(A,A)
        0
        sage: cmp(A,A+1)
        -1
        sage: cmp(A+1,A)
        1
    """
    ########################################################################
    # LEVEL 1 functionality
    # 0 * __cinit__   (not needed)
    # x * __init__
    # 0 * __dealloc__   (not needed)
    # x * set_unsafe
    # x * get_unsafe
    # x * def _pickle
    # x * def _unpickle
    ########################################################################
    def __init__(self, parent, entries, copy, coerce):
        r"""
        See :class:`Matrix_generic_dense` for documentation.

        TESTS:

        We check that the problem related to :trac:`9049` is not an issue any
        more::

            sage: S.<t>=PolynomialRing(QQ)
            sage: F.<q>=QQ.extension(t^4+1)
            sage: R.<x,y>=PolynomialRing(F)
            sage: M = MatrixSpace(R, 1, 2)
            sage: from sage.matrix.matrix_generic_dense import Matrix_generic_dense
            sage: Matrix_generic_dense(M, (x, y), True, True)
            [x y]
        """
        matrix.Matrix.__init__(self, parent)

        cdef Py_ssize_t i,j
        cdef bint is_list

        R = parent.base_ring()
        zero = R.zero()

        # determine if entries is a list or a scalar
        if entries is None:
            entries = zero
            is_list = False
        elif parent_c(entries) is R:
            is_list = False
        elif type(entries) is list:
            # here we do a strong type checking as we potentially want to
            # assign entries to self._entries without copying it
            self._entries = entries
            is_list = True
        elif isinstance(entries, (list,tuple)):
            # it is needed to check for list here as for example Sequence
            # inherits from it but fails the strong type checking above
            self._entries = list(entries)
            is_list = True
            copy = False
        else:
            # not sure what entries is at this point... try scalar first
            try:
                entries = R(entries)
                is_list = False
            except TypeError:
                try:
                    self._entries = list(entries)
                    is_list = True
                    copy = False
                except TypeError:
                    raise TypeError("entries must be coercible to a list or the base ring")

        # now set self._entries
        if is_list:
            if len(self._entries) != self._nrows * self._ncols:
                raise TypeError("entries has the wrong length")
            if coerce:
                self._entries = [R(x) for x in self._entries]
            elif copy:
                self._entries = self._entries[:]
        elif self._nrows == self._ncols:
            self._entries = [zero]*(self._nrows*self._nrows)
            for i in range(self._nrows):
                self._entries[i+self._ncols*i]=entries
        elif entries == zero:
            self._entries = [zero]*(self._nrows*self._ncols)
        else:
            raise TypeError("nonzero scalar matrix must be square")

    cdef Matrix_generic_dense _new(self, Py_ssize_t nrows, Py_ssize_t ncols):
        r"""
        Return a new dense matrix with no entries set.
        """
        cdef Matrix_generic_dense res
        res = self.__class__.__new__(self.__class__, 0, 0, 0)

        if nrows == self._nrows and ncols == self._ncols:
            res._parent = self._parent
        else:
            res._parent = self.matrix_space(nrows, ncols)
        res._ncols  = ncols
        res._nrows  = nrows
        res._base_ring = self._base_ring
        return res

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, value):
        self._entries[i*self._ncols + j] = value

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        return self._entries[i*self._ncols + j]

    def _pickle(self):
        """
        EXAMPLES:
            sage: R.<x> = Integers(25)['x']; A = matrix(R, [1,x,x^3+1,2*x])
            sage: A._pickle()
            ([1, x, x^3 + 1, 2*x], 0)
        """
        return self._entries, 0

    def _unpickle(self, data, int version):
        """
        EXAMPLES:
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

    def __hash__(self):
        """
        EXAMPLES:
            sage: A = random_matrix(Integers(25)['x'],2)
            sage: hash(A)
            Traceback (most recent call last):
            ...
            TypeError: mutable matrices are unhashable
            sage: A.set_immutable()
            sage: hash(A)
            139665060168050560   # 64-bit
            -623270016           # 32-bit
        """
        return self._hash()

    ########################################################################
    # LEVEL 2 functionality
    #    * cdef _add_
    #    * cdef _mul_
    #    * cpdef _cmp_
    #    * __neg__
    #    * __invert__
    # x  * __copy__
    # x  * _multiply_classical
    # x  * _list -- copy of the list of underlying elements
    #    * _dict -- copy of the sparse dictionary of underlying elements
    ########################################################################

    def __copy__(self):
        """
        Creates a copy of self, which may be changed without altering
        self.

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
            <type 'sage.matrix.matrix_generic_dense.Matrix_generic_dense'>
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
        Return reference to list of entries of self.  For internal use
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
