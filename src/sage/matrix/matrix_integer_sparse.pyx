"""
Sparse integer matrices.

AUTHORS:

- William Stein (2007-02-21)
- Soroosh Yazdani (2007-02-21)

TESTS::

    sage: a = matrix(ZZ,2,range(4), sparse=True)
    sage: TestSuite(a).run()
    sage: Matrix(ZZ,0,0,sparse=True).inverse()
    []
"""

##############################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

include 'sage/modules/binary_search.pxi'
include 'sage/modules/vector_integer_sparse_h.pxi'
include 'sage/modules/vector_integer_sparse_c.pxi'
include 'sage/modules/vector_modn_sparse_h.pxi'
include 'sage/modules/vector_modn_sparse_c.pxi'
from cpython.sequence cimport *

include "cysignals/memory.pxi"

from sage.libs.gmp.mpz cimport *
from sage.rings.integer cimport Integer
from matrix cimport Matrix

from matrix_modn_sparse cimport Matrix_modn_sparse
from sage.structure.element cimport ModuleElement, RingElement, Element, Vector

import matrix_space

from sage.rings.integer_ring import ZZ
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing


cdef class Matrix_integer_sparse(Matrix_sparse):

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
        self._initialized = False
        # set the parent, nrows, ncols, etc.
        Matrix_sparse.__init__(self, parent)

        self._matrix = <mpz_vector*> sig_malloc(parent.nrows()*sizeof(mpz_vector))
        if self._matrix == NULL:
            raise MemoryError("error allocating sparse matrix")

        # initialize the rows
        for i from 0 <= i < parent.nrows():
            mpz_vector_init(&self._matrix[i], self._ncols, 0)
        # record that rows have been initialized
        self._initialized = True

    def __dealloc__(self):
        cdef Py_ssize_t i
        if self._initialized:
            for i from 0 <= i < self._nrows:
                mpz_vector_clear(&self._matrix[i])
        sig_free(self._matrix)

    def __init__(self, parent, entries, copy, coerce):
        """
        Create a sparse matrix over the integers.

        INPUT:

        - ``parent`` -- a matrix space

        - ``entries`` -- can be one of the following:

          * a Python dictionary whose items have the
            form ``(i, j): x``, where ``0 <= i < nrows``,
            ``0 <= j < ncols``, and ``x`` is coercible to
            an integer.  The ``i,j`` entry of ``self`` is
            set to ``x``.  The ``x``'s can be ``0``.
          * Alternatively, entries can be a list of *all*
            the entries of the sparse matrix, read
            row-by-row from top to bottom (so they would
            be mostly 0).

        - ``copy`` -- ignored

        - ``coerce`` -- ignored
        """
        cdef Py_ssize_t i, j, k
        cdef Integer z
        cdef PyObject** X

        # fill in entries in the dict case
        if entries is None: return
        if isinstance(entries, dict):
            R = self.base_ring()
            for ij, x in entries.iteritems():
                z = R(x)
                if z != 0:
                    i, j = ij  # nothing better to do since this is user input, which may be bogus.
                    if i < 0 or j < 0 or i >= self._nrows or j >= self._ncols:
                        raise IndexError("invalid entries list")
                    mpz_vector_set_entry(&self._matrix[i], j, z.value)

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
                        mpz_vector_set_entry(&self._matrix[i], j, z.value)
                    k = k + 1

        else:

            # fill in entries in the scalar case
            z = Integer(entries)
            if z == 0:
                return
            if self._nrows != self._ncols:
                raise TypeError("matrix must be square to initialize with a scalar.")
            for i from 0 <= i < self._nrows:
                mpz_vector_set_entry(&self._matrix[i], i, z.value)


    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, x):
        mpz_vector_set_entry(&self._matrix[i], j, (<Integer> x).value)

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        cdef Integer x
        x = Integer()
        mpz_vector_get_entry(x.value, &self._matrix[i], j)
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

    cpdef _lmul_(self, RingElement right):
        """
        EXAMPLES::

            sage: a = matrix(ZZ,2,range(6), sparse=True)
            sage: 3 * a
            [ 0  3  6]
            [ 9 12 15]
        """
        cdef Py_ssize_t i
        cdef mpz_vector *self_row
        cdef mpz_vector *M_row
        cdef Matrix_integer_sparse M
        cdef Integer _x
        _x = Integer(right)
        M = Matrix_integer_sparse.__new__(Matrix_integer_sparse, self._parent, None, None, None)
        for i from 0 <= i < self._nrows:
            self_row = &self._matrix[i]
            M_row = &M._matrix[i]
            mpz_vector_scalar_multiply(M_row, self_row, _x.value)
        return M

    cpdef _add_(self, right):
        cdef Py_ssize_t i, j
        cdef mpz_vector *self_row
        cdef mpz_vector *M_row
        cdef Matrix_integer_sparse M

        M = Matrix_integer_sparse.__new__(Matrix_integer_sparse, self._parent, None, None, None)
        cdef mpz_t mul
        mpz_init_set_si(mul,1)
        for i from 0 <= i < self._nrows:
            mpz_vector_clear(&M._matrix[i])
            add_mpz_vector_init(&M._matrix[i], &self._matrix[i], &(<Matrix_integer_sparse>right)._matrix[i], mul)
        mpz_clear(mul)
        return M

    cpdef _sub_(self, right):
        cdef Py_ssize_t i, j
        cdef mpz_vector *self_row
        cdef mpz_vector *M_row
        cdef Matrix_integer_sparse M

        M = Matrix_integer_sparse.__new__(Matrix_integer_sparse, self._parent, None, None, None)
        cdef mpz_t mul
        mpz_init_set_si(mul,-1)
        for i from 0 <= i < self._nrows:
            mpz_vector_clear(&M._matrix[i])
            add_mpz_vector_init(&M._matrix[i], &self._matrix[i], &(<Matrix_integer_sparse>right)._matrix[i], mul)
        mpz_clear(mul)
        return M

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
                x = Integer()
                mpz_set((<Integer>x).value, self._matrix[i].entries[j])
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

            sage: M = Matrix(ZZ, [[0,0,0,1,0,0,0,0],[0,1,0,0,0,0,1,0]], sparse=True); M
            [0 0 0 1 0 0 0 0]
            [0 1 0 0 0 0 1 0]
            sage: M._nonzero_positions_by_row()
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

    def _nonzero_positions_by_column(self, copy=True):
        """
        Returns the list of pairs (i,j) such that self[i,j] != 0, but
        sorted by columns, i.e., column j=0 entries occur first, then
        column j=1 entries, etc.

        It is safe to change the resulting list (unless you give the option copy=False).

        EXAMPLE::

            sage: M = Matrix(ZZ, [[0,0,0,1,0,0,0,0],[0,1,0,0,0,0,1,0]], sparse=True); M
            [0 0 0 1 0 0 0 0]
            [0 1 0 0 0 0 1 0]
            sage: M._nonzero_positions_by_column()
            [(1, 1), (0, 3), (1, 6)]

        """
        x = self.fetch('nonzero_positions_by_column')
        if not x is None:
            if copy:
                return list(x)
            return x
        nzc = [[] for _ in xrange(self._ncols)]
        cdef Py_ssize_t i, j
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._matrix[i].num_nonzero:
                p = self._matrix[i].positions[j]
                nzc[p].append((i,p))
        nzc = sum(nzc,[])
        self.cache('nonzero_positions_by_column', nzc)
        if copy:
            return list(nzc)
        return nzc

    def _mod_int(self, modulus):
        """
        Helper function in reducing matrices mod n.

        INPUT:

        - `modulus` - a number

        OUTPUT:

        This matrix, over `ZZ/nZZ`.

        TESTS::

            sage: M = Matrix(ZZ, sparse=True)
            sage: B = M._mod_int(7)
            sage: B.parent()
            Full MatrixSpace of 0 by 0 sparse matrices over Ring of integers modulo 7

        """
        return self._mod_int_c(modulus)

    cdef _mod_int_c(self, mod_int p):
        cdef Py_ssize_t i, j
        cdef Matrix_modn_sparse res
        cdef mpz_vector* self_row
        cdef c_vector_modint* res_row
        res = Matrix_modn_sparse.__new__(Matrix_modn_sparse, matrix_space.MatrixSpace(
            IntegerModRing(p), self._nrows, self._ncols, sparse=True), None, None, None)
        for i from 0 <= i < self._nrows:
            self_row = &(self._matrix[i])
            res_row = &(res.rows[i])
            for j from 0 <= j < self_row.num_nonzero:
                set_entry(res_row, self_row.positions[j], mpz_fdiv_ui(self_row.entries[j], p))
        return res


    def rational_reconstruction(self, N):
        """
        Use rational reconstruction to lift self to a matrix over the
        rational numbers (if possible), where we view self as a matrix
        modulo N.

        EXAMPLES::

            sage: A = matrix(ZZ, 3, 4, [(1/3)%500, 2, 3, (-4)%500, 7, 2, 2, 3, 4, 3, 4, (5/7)%500], sparse=True)
            sage: A.rational_reconstruction(500)
            [1/3   2   3  -4]
            [  7   2   2   3]
            [  4   3   4 5/7]

        TEST:

        Check that :trac:`9345` is fixed::

            sage: A = random_matrix(ZZ, 3, 3, sparse = True)
            sage: A.rational_reconstruction(0)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: The modulus cannot be zero
        """
        import misc
        return misc.matrix_integer_sparse_rational_reconstruction(self, N)

    def _right_kernel_matrix(self, **kwds):
        r"""
        Returns a pair that includes a matrix of basis vectors
        for the right kernel of ``self``.

        INPUT:

        - ``algorithm`` - determines which algorithm to use, options are:

          - 'pari' - use the ``matkerint()`` function from the PARI library
          - 'padic' - use the p-adic algorithm from the IML library
          - 'default' - use a heuristic to decide which of the two above
            routines is fastest.  This is the default value.

        - ``proof`` - this is passed to the p-adic IML algorithm.
          If not specified, the global flag for linear algebra will be used.

        OUTPUT:

        Returns a pair.  First item is the string is either
        'computed-pari-int', 'computed-iml-int' or 'computed-flint-int', which identifies
        the nature of the basis vectors.

        Second item is a matrix whose rows are a basis for the right kernel,
        over the integers, as computed by either the IML or PARI libraries.

        EXAMPLES::

            sage: A = matrix(ZZ, [[4, 7, 9, 7, 5, 0],
            ...                   [1, 0, 5, 8, 9, 1],
            ...                   [0, 1, 0, 1, 9, 7],
            ...                   [4, 7, 6, 5, 1, 4]],
            ...              sparse = True)

            sage: result = A._right_kernel_matrix(algorithm='pari')
            sage: result[0]
            'computed-pari-int'
            sage: X = result[1]; X
            [ 26 -31  30 -21  -2  10]
            [ 47  13 -48  14  11 -18]
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

            sage: result = A._right_kernel_matrix()
            sage: result[0]
            'computed-flint-int'
            sage: result[1]
            [ 469 -214   30 -119   37    0]
            [-370  165  -18   91  -30    2]

        With the 'default' given as the algorithm, several heuristics are
        used to determine if PARI or IML ('padic') is used.  The code has
        exact details, but roughly speaking, relevant factors are: the
        absolute size of the matrix, or the relative dimensions, or the
        magnitude of the entries. ::

            sage: A = random_matrix(ZZ, 18, 11, sparse=True)
            sage: A._right_kernel_matrix(algorithm='default')[0]
            'computed-pari-int'
            sage: A = random_matrix(ZZ, 18, 11, x = 10^200, sparse=True)
            sage: A._right_kernel_matrix(algorithm='default')[0]
            'computed-iml-int'
            sage: A = random_matrix(ZZ, 60, 60, sparse=True)
            sage: A._right_kernel_matrix(algorithm='default')[0]
            'computed-iml-int'
            sage: A = random_matrix(ZZ, 60, 55, sparse=True)
            sage: A._right_kernel_matrix(algorithm='default')[0]
            'computed-pari-int'

        TESTS:

        We test three trivial cases. PARI is used for small matrices,
        but we let the heuristic decide that. ::

            sage: A = matrix(ZZ, 0, 2, sparse=True)
            sage: A._right_kernel_matrix()[1]
            []
            sage: A = matrix(ZZ, 2, 0, sparse=True)
            sage: A._right_kernel_matrix()[1].parent()
            Full MatrixSpace of 0 by 0 dense matrices over Integer Ring
            sage: A = zero_matrix(ZZ, 4, 3, sparse=True)
            sage: A._right_kernel_matrix()[1]
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        return self.dense_matrix()._right_kernel_matrix(**kwds)

    hermite_form = Matrix.echelon_form

    def elementary_divisors(self, algorithm='pari'):
        """
        Return the elementary divisors of self, in order.

        The elementary divisors are the invariants of the finite
        abelian group that is the cokernel of *left* multiplication by
        this matrix.  They are ordered in reverse by divisibility.

        INPUT:

        - self -- matrix
        - algorithm -- (default: 'pari')

          * 'pari': works robustly, but is slower.
          * 'linbox' -- use linbox (currently off, broken)

        OUTPUT:

        list of integers

        EXAMPLES::

            sage: matrix(3, range(9),sparse=True).elementary_divisors()
            [1, 3, 0]
            sage: M = matrix(ZZ, 3, [1,5,7, 3,6,9, 0,1,2], sparse=True)
            sage: M.elementary_divisors()
            [1, 1, 6]

        This returns a copy, which is safe to change::

            sage: edivs = M.elementary_divisors()
            sage: edivs.pop()
            6
            sage: M.elementary_divisors()
            [1, 1, 6]

        ..SEEALSO:: :meth:`smith_form`
        """
        return self.dense_matrix().elementary_divisors(algorithm=algorithm)

    def smith_form(self):
        r"""
        Returns matrices S, U, and V such that S = U\*self\*V, and S is in
        Smith normal form. Thus S is diagonal with diagonal entries the
        ordered elementary divisors of S.

        This version is for sparse matrices and simply makes the matrix
        dense and calls the version for dense integer matrices.

        .. warning::

           The elementary_divisors function, which returns the
           diagonal entries of S, is VASTLY faster than this function.

        The elementary divisors are the invariants of the finite abelian
        group that is the cokernel of this matrix. They are ordered in
        reverse by divisibility.

        EXAMPLES::

            sage: A = MatrixSpace(IntegerRing(), 3, sparse=True)(range(9))
            sage: D, U, V = A.smith_form()
            sage: D
            [1 0 0]
            [0 3 0]
            [0 0 0]
            sage: U
            [ 0  1  0]
            [ 0 -1  1]
            [-1  2 -1]
            sage: V
            [-1  4  1]
            [ 1 -3 -2]
            [ 0  0  1]
            sage: U*A*V
            [1 0 0]
            [0 3 0]
            [0 0 0]

        It also makes sense for nonsquare matrices::

            sage: A = Matrix(ZZ,3,2,range(6), sparse=True)
            sage: D, U, V = A.smith_form()
            sage: D
            [1 0]
            [0 2]
            [0 0]
            sage: U
            [ 0  1  0]
            [ 0 -1  1]
            [-1  2 -1]
            sage: V
            [-1  3]
            [ 1 -2]
            sage: U * A * V
            [1 0]
            [0 2]
            [0 0]

        The examples above show that :trac:`10626` has been implemented.


        .. seealso::

           :meth:`elementary_divisors`
        """
        return self.dense_matrix().smith_form()
