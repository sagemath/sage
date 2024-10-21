# sage_setup: distribution = sagemath-linbox
r"""
Sparse integer matrices

AUTHORS:

- William Stein (2007-02-21)
- Soroosh Yazdani (2007-02-21)

TESTS::

    sage: a = matrix(ZZ,2,range(4), sparse=True)
    sage: TestSuite(a).run()
    sage: Matrix(ZZ,0,0,sparse=True).inverse()
    []
"""

# ****************************************************************************
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#       Copyright (C) 2018 Vincent Delecroix <20100.delecroix@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.memory cimport check_calloc, sig_free
from cysignals.signals cimport sig_on, sig_off

from cpython.long cimport PyLong_FromSize_t

from sage.ext.stdsage cimport PY_NEW
from sage.ext.mod_int cimport *

cimport sage.libs.linbox.givaro as givaro
cimport sage.libs.linbox.linbox as linbox
from sage.libs.linbox.conversion cimport (
    new_linbox_vector_integer_dense,
    new_sage_vector_integer_dense,
    new_linbox_matrix_integer_sparse,
    METHOD_DEFAULT, METHOD_DENSE_ELIMINATION,
    METHOD_SPARSE_ELIMINATION, METHOD_BLACKBOX,
    METHOD_WIEDEMANN, get_method)

from sage.data_structures.binary_search cimport *
from sage.modules.vector_integer_sparse cimport *
from sage.modules.vector_integer_dense cimport Vector_integer_dense
from sage.modules.vector_modn_sparse cimport *

from sage.libs.gmp.mpz cimport *

from sage.rings.integer cimport Integer
from sage.rings.polynomial.polynomial_integer_dense_flint cimport Polynomial_integer_dense_flint
from sage.matrix.matrix cimport Matrix

from sage.matrix.args cimport SparseEntry, MatrixArgs_init
from sage.matrix.matrix_integer_dense cimport Matrix_integer_dense
from sage.libs.flint.fmpz cimport fmpz_set_mpz, fmpz_get_mpz
from sage.libs.flint.fmpz_poly cimport fmpz_poly_fit_length, _fmpz_poly_set_length
from sage.libs.flint.fmpz_poly_sage cimport fmpz_poly_set_coeff_mpz
from sage.libs.flint.fmpz_mat cimport fmpz_mat_entry

from sage.matrix.matrix_modn_sparse cimport Matrix_modn_sparse
from sage.structure.element cimport Element

import sage.matrix.matrix_space as matrix_space

from sage.rings.integer_ring import ZZ
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing

cimport sage.structure.element

cdef class Matrix_integer_sparse(Matrix_sparse):
    def __cinit__(self):
        self._matrix = <mpz_vector*>check_calloc(self._nrows, sizeof(mpz_vector))
        # Initialize the rows
        cdef Py_ssize_t i
        for i in range(self._nrows):
            mpz_vector_init(&self._matrix[i], self._ncols, 0)

    def __dealloc__(self):
        cdef Py_ssize_t i
        if self._matrix is not NULL:
            for i in range(self._nrows):
                mpz_vector_clear(&self._matrix[i])
            sig_free(self._matrix)

    def __init__(self, parent, entries=None, copy=None, bint coerce=True):
        r"""
        Create a sparse matrix over the integers.

        INPUT:

        - ``parent`` -- a matrix space over ``ZZ``

        - ``entries`` -- see :func:`matrix`

        - ``copy`` -- ignored (for backwards compatibility)

        - ``coerce`` -- if ``False``, assume without checking that the
          entries are of type :class:`Integer`
        """
        ma = MatrixArgs_init(parent, entries)
        cdef Integer z
        for t in ma.iter(coerce, True):
            se = <SparseEntry>t
            z = <Integer>se.entry
            if z:
                mpz_vector_set_entry(&self._matrix[se.i], se.j, z.value)

    cdef set_unsafe(self, Py_ssize_t i, Py_ssize_t j, x):
        mpz_vector_set_entry(&self._matrix[i], j, (<Integer> x).value)

    cdef get_unsafe(self, Py_ssize_t i, Py_ssize_t j):
        cdef Integer x
        x = Integer()
        mpz_vector_get_entry(x.value, &self._matrix[i], j)
        return x

    cdef bint get_is_zero_unsafe(self, Py_ssize_t i, Py_ssize_t j) except -1:
        """
        Return 1 if the entry ``(i, j)`` is zero, otherwise 0.

        EXAMPLES::

            sage: M = matrix(ZZ, [[0,1,0],[0,0,0]], sparse=True)
            sage: M.zero_pattern_matrix()  # indirect doctest
            [1 0 1]
            [1 1 1]
        """
        return mpz_vector_is_entry_zero_unsafe(&self._matrix[i], j)

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

    cpdef _lmul_(self, Element right):
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
        cdef Py_ssize_t i
        cdef Matrix_integer_sparse M

        M = Matrix_integer_sparse.__new__(Matrix_integer_sparse, self._parent, None, None, None)
        cdef mpz_t mul
        mpz_init_set_si(mul,1)
        for i in range(self._nrows):
            mpz_vector_clear(&M._matrix[i])
            add_mpz_vector_init(&M._matrix[i], &self._matrix[i],
                                &(<Matrix_integer_sparse>right)._matrix[i], mul)
        mpz_clear(mul)
        return M

    cpdef _sub_(self, right):
        cdef Py_ssize_t i
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
        if d is not None:
            return d

        cdef Py_ssize_t i, j
        d = {}
        for i in range(self._nrows):
            for j in range(self._matrix[i].num_nonzero):
                x = Integer()
                mpz_set((<Integer>x).value, self._matrix[i].entries[j])
                d[(int(i), int(self._matrix[i].positions[j]))] = x
        self.cache('dict', d)
        return d

    cdef sage.structure.element.Matrix _matrix_times_matrix_(self, sage.structure.element.Matrix _right):
        """
        Return the product of the sparse integer matrices
        ``self`` and ``_right``.

       EXAMPLES::

            sage: a = matrix(ZZ, 2, [1,2,3,4], sparse=True)
            sage: b = matrix(ZZ, 2, 3, [1..6], sparse=True)
            sage: a * b
            [ 9 12 15]
            [19 26 33]
        """
        cdef Matrix_integer_sparse right, ans
        right = _right

        cdef mpz_vector* v

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
        cdef mpz_t x, y, s
        mpz_init(x)
        mpz_init(y)
        mpz_init(s)
        for i in range(self._nrows):
            v = &(self._matrix[i])
            if not v.num_nonzero:
                continue
            for j in right_indices:
                mpz_set_si(s, 0)
                c = <set> nonzero_positions_in_columns[j]
                for k in range(v.num_nonzero):
                    if v.positions[k] in c:
                        mpz_vector_get_entry(y, &right._matrix[v.positions[k]], j)
                        mpz_mul(x, v.entries[k], y)
                        mpz_add(s, s, x)
                mpz_vector_set_entry(&ans._matrix[i], j, s)

        mpz_clear(x)
        mpz_clear(y)
        mpz_clear(s)
        return ans

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
        Return the list of pairs (i,j) such that ``self[i,j] != 0``.

        It is safe to change the resulting list (unless you give the option copy=False).

        EXAMPLES::

            sage: M = Matrix(ZZ, [[0,0,0,1,0,0,0,0],[0,1,0,0,0,0,1,0]], sparse=True); M
            [0 0 0 1 0 0 0 0]
            [0 1 0 0 0 0 1 0]
            sage: M._nonzero_positions_by_row()
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

    def _nonzero_positions_by_column(self, copy=True):
        """
        Return the list of pairs (i,j) such that ``self[i,j] != 0``, but
        sorted by columns, i.e., column j=0 entries occur first, then
        column j=1 entries, etc.

        It is safe to change the resulting list (unless you give the option copy=False).

        EXAMPLES::

            sage: M = Matrix(ZZ, [[0,0,0,1,0,0,0,0],[0,1,0,0,0,0,1,0]], sparse=True); M
            [0 0 0 1 0 0 0 0]
            [0 1 0 0 0 0 1 0]
            sage: M._nonzero_positions_by_column()
            [(1, 1), (0, 3), (1, 6)]
        """
        x = self.fetch('nonzero_positions_by_column')
        if x is not None:
            if copy:
                return list(x)
            return x
        nzc = [[] for _ in range(self._ncols)]
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
        r"""
        Helper function in reducing matrices mod n.

        INPUT:

        - ``modulus`` -- a number

        OUTPUT: this matrix, over `\ZZ/n\ZZ`

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
        Use rational reconstruction to lift ``self`` to a matrix over the
        rational numbers (if possible), where we view ``self`` as a matrix
        modulo `N`.

        EXAMPLES::

            sage: A = matrix(ZZ, 3, 4, [(1/3)%500, 2, 3, (-4)%500,
            ....:                       7, 2, 2, 3,
            ....:                       4, 3, 4, (5/7)%500], sparse=True)
            sage: A.rational_reconstruction(500)
            [1/3   2   3  -4]
            [  7   2   2   3]
            [  4   3   4 5/7]

        TESTS:

        Check that :issue:`9345` is fixed::

            sage: A = random_matrix(ZZ, 3, 3, sparse=True)
            sage: A.rational_reconstruction(0)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: The modulus cannot be zero
        """
        from sage.matrix.misc import matrix_integer_sparse_rational_reconstruction
        return matrix_integer_sparse_rational_reconstruction(self, N)

    def _right_kernel_matrix(self, **kwds):
        r"""
        Return a pair that includes a matrix of basis vectors
        for the right kernel of ``self``.

        INPUT:

        - ``algorithm`` -- determines which algorithm to use, options are:

          - ``'pari'`` -- use the :pari:`matkerint` function from the PARI library
          - ``'padic'`` -- use the `p`-adic algorithm from the IML library
          - ``'default'`` -- use a heuristic to decide which of the two above
            routines is fastest.  This is the default value.

        - ``proof`` -- this is passed to the `p`-adic IML algorithm
          If not specified, the global flag for linear algebra will be used

        OUTPUT:

        Returns a pair.  First item is the string is either
        'computed-pari-int', 'computed-iml-int' or 'computed-flint-int', which identifies
        the nature of the basis vectors.

        Second item is a matrix whose rows are a basis for the right kernel,
        over the integers, as computed by either the IML or PARI libraries.

        EXAMPLES::

            sage: A = matrix(ZZ, [[4, 7, 9, 7, 5, 0],
            ....:                 [1, 0, 5, 8, 9, 1],
            ....:                 [0, 1, 0, 1, 9, 7],
            ....:                 [4, 7, 6, 5, 1, 4]],
            ....:            sparse = True)

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

        - ``self`` -- matrix
        - ``algorithm`` -- (default: ``'pari'``)

          * 'pari': works robustly, but is slower
          * 'linbox' -- use linbox (currently off, broken)

        OUTPUT: list of integers

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

        .. SEEALSO::

            :meth:`smith_form`
        """
        return self.dense_matrix().elementary_divisors(algorithm=algorithm)

    def smith_form(self, transformation=True, integral=None):
        r"""
        Return the smith normal form of this matrix, that is the diagonal
        matrix `S` with diagonal entries the ordered elementary divisors of
        this matrix.

        INPUT:

        - ``transformation`` -- boolean (default: ``True``); whether to
          return the transformation matrices `U` and `V` such that `S = U\cdot
          self\cdot V`

        - ``integral`` -- a subring of the base ring or ``True`` (default:
          ``None``); ignored for matrices with integer entries

        This version is for sparse matrices and simply makes the matrix
        dense and calls the version for dense integer matrices.

        .. NOTE::

           The :meth:`elementary_divisors` function, which returns the diagonal
           entries of S, is VASTLY faster than this function.

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

            sage: A = Matrix(ZZ,3,2,range(6), sparse=True)
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

        The examples above show that :issue:`10626` has been implemented.


        .. SEEALSO::

           :meth:`elementary_divisors`
        """
        return self.dense_matrix().smith_form(transformation=transformation)

    def rank(self, algorithm=None):
        r"""
        Compute the rank of this matrix.

        INPUT:

        - ``algorithm`` -- (optional) one of ``None``, ``'linbox'`` or
          ``'generic'``

        EXAMPLES::

            sage: M = MatrixSpace(ZZ, 3, 2, sparse=True)
            sage: m = M([1, 0, 2, 3, -1, 0])
            sage: m.rank()
            2
        """
        if self._nrows == 0 or self._ncols == 0:
            return 0

        x = self.fetch('rank')
        if x is not None:
            return x

        if algorithm is None or algorithm == "linbox":
            r = self._rank_linbox()
            self.cache("rank", r)
        elif algorithm == "generic":
            r = Matrix_sparse.rank(self)
            self.cache("rank", r)
        else:
            raise ValueError("no algorithm '%s'"%algorithm)

        return r

    def _rank_linbox(self):
        r"""
        Compute the rank using linbox.

        The result is not cached contrarily to the method ``rank``.

        EXAMPLES::

            sage: M = MatrixSpace(ZZ, 3, sparse=True)
            sage: m = M([1,0,1,0,2,0,2,0,2])
            sage: m._rank_linbox()
            2

        TESTS::

            sage: MatrixSpace(ZZ, 0, 0, sparse=True)()._rank_linbox()
            0
            sage: MatrixSpace(ZZ, 1, 0, sparse=True)()._rank_linbox()
            0
            sage: MatrixSpace(ZZ, 0, 1, sparse=True)()._rank_linbox()
            0
            sage: MatrixSpace(ZZ, 1, 1, sparse=True)()._rank_linbox()
            0
        """
        if self._nrows == 0 or self._ncols == 0:
            return 0

        cdef givaro.ZRing givZZ
        cdef linbox.SparseMatrix_integer * M = new_linbox_matrix_integer_sparse(givZZ, self)
        cdef size_t r = 0

        sig_on()
        linbox.rank(r, M[0])
        sig_off()

        del M

        return PyLong_FromSize_t(r)

    def _det_linbox(self):
        r"""
        Return the determinant computed with LinBox.

        .. NOTE::

            This method is much slower than converting to a dense matrix and
            computing the determinant there. There is not much point in making
            it available. See :issue:`28318`.

        EXAMPLES::

            sage: M = MatrixSpace(ZZ, 2, 2, sparse=True)
            sage: M([2,0,1,1])._det_linbox()
            2

        TESTS::

            sage: MatrixSpace(ZZ, 0, 0, sparse=True)()._det_linbox()
            1
            sage: MatrixSpace(ZZ, 1, 1, sparse=True)()._det_linbox()
            0

            sage: m = diagonal_matrix(ZZ, [2] * 46)
            sage: m._det_linbox() == 2**46
            True

            sage: m = diagonal_matrix(ZZ, [3] * 100)
            sage: m._det_linbox() == 3**100
            True
        """
        if self._nrows != self._ncols:
            raise ValueError("non square matrix")

        if self._nrows == 0:
            return Integer(1)

        cdef givaro.ZRing givZZ
        cdef linbox.SparseMatrix_integer * M = new_linbox_matrix_integer_sparse(givZZ, self)
        cdef givaro.Integer D

        sig_on()
        linbox.det(D, M[0])
        sig_off()

        cdef Integer d = PY_NEW(Integer)
        mpz_set(d.value, D.get_mpz_const())

        del M

        return d

    def charpoly(self, var='x', algorithm=None):
        r"""
        Return the characteristic polynomial of this matrix.

        INPUT:

        - ``var`` -- (default: ``'x'``) the name of the variable
          of the polynomial

        - ``algorithm`` -- (default: ``None``) one of ``None``,
          ``'linbox'``, or an algorithm accepted by
          :meth:`sage.matrix.matrix_sparse.Matrix_sparse.charpoly`

        EXAMPLES::

            sage: M = MatrixSpace(ZZ, 4, sparse=True)
            sage: m = M()
            sage: m[0,0] = m[1,2] = m[2,3] = m[3,3] = 1
            sage: m[0,2] = m[1,3] = m[2,0] = m[3,0] = -3
            sage: m[1,1] = 2
            sage: m
            [ 1  0 -3  0]
            [ 0  2  1 -3]
            [-3  0  0  1]
            [-3  0  0  1]
            sage: m.charpoly()
            x^4 - 4*x^3 - 4*x^2 + 16*x

        TESTS::

            sage: matrix(ZZ, 0, 0, sparse=True).charpoly(algorithm='linbox')
            1
            sage: matrix(ZZ, 0, 0, sparse=True).charpoly(algorithm='generic')
            1
        """
        if self._nrows != self._ncols:
            raise ArithmeticError("only valid for square matrix")

        cdef Polynomial_integer_dense_flint g

        if algorithm is None:
            algorithm = 'linbox'

        g = self.fetch('charpoly')
        if g is not None:
            return g.change_variable_name(var)

        if algorithm == 'linbox':
            g = self._charpoly_linbox()
        else:
            g = Matrix_sparse.charpoly(self, var, algorithm=algorithm)

        self.cache('charpoly', g)
        return g

    def _charpoly_linbox(self, var='x'):
        r"""
        Compute the charpoly using LinBox.

        EXAMPLES::

            sage: m = matrix(ZZ, 2, [2,1,1,1], sparse=True)
            sage: m._charpoly_linbox()
            x^2 - 3*x + 1

        TESTS::

            sage: matrix(ZZ, 0, 0, sparse=True)._charpoly_linbox()
            1
            sage: matrix(ZZ, 1, 1, sparse=True)._charpoly_linbox()
            x
        """
        cdef mpz_t tmp
        if self._nrows != self._ncols:
            raise ArithmeticError('only valid for square matrix')

        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        R = PolynomialRing(ZZ, names=var)

        # TODO: bug in LinBox (got SIGSEGV)
        if self._nrows == 0:
            return R.one()

        cdef givaro.ZRing givZZ
        cdef linbox.SparseMatrix_integer * M = new_linbox_matrix_integer_sparse(givZZ, self)
        cdef linbox.DensePolynomial_integer * p = new linbox.DensePolynomial_integer(givZZ, <size_t> self._nrows)
        cdef Polynomial_integer_dense_flint g = (<Polynomial_integer_dense_flint> R.gen())._new()

        sig_on()
        linbox.charpoly(p[0], M[0])
        sig_off()

        cdef size_t i
        fmpz_poly_fit_length(g._poly, p.size())
        for i in range(p.size()):
            tmp = p[0][i].get_mpz_const()
            fmpz_poly_set_coeff_mpz(g._poly, i, tmp)
        _fmpz_poly_set_length(g._poly, p.size())

        del M
        del p

        return g

    def minpoly(self, var='x', algorithm=None):
        r"""
        Return the minimal polynomial of this matrix.

        INPUT:

        - ``var`` -- (default: ``'x'``) the name of the variable
          of the polynomial

        - ``algorithm`` -- (default: ``None``) one of ``None``,
          ``'linbox'``, or an algorithm accepted by
          :meth:`sage.matrix.matrix_sparse.Matrix_sparse.minpoly`

        EXAMPLES::

            sage: M = MatrixSpace(ZZ, 4, sparse=True)
            sage: m = M({(0, 0):2, (1, 1):1, (2, 1):-8, (2, 2):2, (2, 3):-1, (3, 3):1})
            sage: m
            [ 2  0  0  0]
            [ 0  1  0  0]
            [ 0 -8  2 -1]
            [ 0  0  0  1]
            sage: m.minpoly()
            x^2 - 3*x + 2

        TESTS::

            sage: matrix(ZZ, 0, 0, sparse=True).minpoly(algorithm='linbox')
            1
            sage: matrix(ZZ, 0, 0, sparse=True).minpoly(algorithm='generic')
            1
        """
        if self._nrows != self._ncols:
            raise ArithmeticError("only valid for square matrix")

        cdef Polynomial_integer_dense_flint g

        if algorithm is None:
            algorithm = 'linbox'

        g = self.fetch('minpoly')
        if g is not None:
            return g.change_variable_name(var)

        if algorithm == 'linbox':
            g = self._minpoly_linbox()
        else:
            g = Matrix_sparse.minpoly(self, var, algorithm=algorithm)

        self.cache('minpoly', g)
        return g

    def _minpoly_linbox(self, var='x'):
        r"""
        Compute the minpoly using LinBox.

        EXAMPLES::

            sage: m = matrix(ZZ, 2, [2,1,1,1], sparse=True)
            sage: m._minpoly_linbox()
            x^2 - 3*x + 1

        TESTS::

            sage: matrix(ZZ, 0, 0, sparse=True)._minpoly_linbox()
            1
            sage: matrix(ZZ, 1, 1, sparse=True)._minpoly_linbox()
            x
        """
        if self._nrows != self._ncols:
            raise ArithmeticError('only valid for square matrix')

        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        R = PolynomialRing(ZZ, names=var)

        # TODO: bug in LinBox (got SIGSEGV)
        if self._nrows == 0:
            return R.one()

        cdef givaro.ZRing givZZ
        cdef linbox.SparseMatrix_integer * M = new_linbox_matrix_integer_sparse(givZZ, self)
        cdef linbox.DensePolynomial_integer * p = new linbox.DensePolynomial_integer(givZZ, <size_t> self._nrows)
        cdef Polynomial_integer_dense_flint g = (<Polynomial_integer_dense_flint> R.gen())._new()

        sig_on()
        linbox.minpoly(p[0], M[0])
        sig_off()

        cdef size_t i
        cdef mpz_t tmp
        fmpz_poly_fit_length(g._poly, p.size())
        for i in range(p.size()):
            tmp = p[0][i].get_mpz_const()
            fmpz_poly_set_coeff_mpz(g._poly, i, tmp)
        _fmpz_poly_set_length(g._poly, p.size())

        del M
        del p

        return g

    def _solve_right_nonsingular_square(self, B, algorithm=None, check_rank=False):
        r"""
        If ``self`` is a matrix `A`, then this function returns a
        vector or matrix `X` such that `A X = B`. If
        `B` is a vector then `X` is a vector and if
        `B` is a matrix, then `X` is a matrix.

        .. NOTE::

           DEPRECATED. In Sage one can also write ``A \ B`` for
           ``A.solve_right(B)``, i.e., Sage implements the "the
           MATLAB/Octave backslash operator".

        INPUT:

        - ``B`` -- a matrix or vector

        - ``algorithm`` -- one of the following:

            - ``'linbox'`` or ``'linbox_default'`` -- (default) use LinBox
              and let it chooses the appropriate algorithm

            - ``linbox_dense_elimination'`` -- use LinBox dense elimination

            - ``'linbox_sparse_elimination'`` -- use LinBox sparse elimination

            - ``'linbox_ blackbox'`` -- LinBox via a Blackbox algorithm

            - ``'linbox_wiedemann'`` -- use LinBox implementation of
              Wiedemann's algorithm

            - ``'generic'`` -- use the Sage generic implementation
              (via inversion)

        - ``check_rank`` -- whether to check that the rank is maximal

        OUTPUT: a matrix or vector

        EXAMPLES::

            sage: A = matrix(ZZ, 3, [1,2,3,-1,2,5,2,3,1], sparse=True)
            sage: b = vector(ZZ, [1,2,3])
            sage: x = A.solve_right(b)
            sage: x
            (-13/12, 23/12, -7/12)
            sage: A * x
            (1, 2, 3)

            sage: u = matrix(ZZ, 3, 2, [0,1,1,1,0,2])
            sage: x = A.solve_right(u)
            sage: x
            [-7/12  -1/6]
            [ 5/12   5/6]
            [-1/12  -1/6]
            sage: A * x
            [0 1]
            [1 1]
            [0 2]
        """
        if check_rank and self.rank() < self.nrows():
            from sage.matrix.matrix2 import NotFullRankError
            raise NotFullRankError

        if self.base_ring() != B.base_ring():
            B = B.change_ring(self.base_ring())
        if self.nrows() != B.nrows():
            raise ValueError("input matrices must have the same number of rows.")

        if algorithm == "generic":
            return Matrix_sparse.solve_right(self, B)
        else:
            if isinstance(B, sage.structure.element.Matrix):
                from sage.rings.rational_field import QQ
                from sage.matrix.special import diagonal_matrix
                m, d = self._solve_matrix_linbox(B, algorithm)
                return m  * diagonal_matrix([QQ((1,x)) for x in d])
            else:
                v, d = self._solve_vector_linbox(B, algorithm)
                return v / d

    def _solve_vector_linbox(self, v, algorithm=None):
        r"""
        Return a pair ``(a, d)`` so that ``d * b = m * a``.

        If there is no solution a :exc:`ValueError` is raised.

        INPUT:

        - ``b`` -- dense integer vector

        - ``algorithm`` -- (optional) either ``None``, ``'dense_elimination'``,
          ``'sparse_elimination'``, ``'wiedemann'`` or ``'blackbox'``

        OUTPUT: a pair ``(a, d)`` consisting of

        - ``a`` -- dense integer vector

        - ``d`` -- integer

        EXAMPLES::

            sage: m = matrix(ZZ, 4, sparse=True)
            sage: m[0,0] = m[1,2] = m[2,0] = m[3,3] = 2
            sage: m[0,2] = m[1,1] = -1
            sage: m[2,3] = m[3,0] = -3

            sage: b0 = vector((1,1,1,1))
            sage: m._solve_vector_linbox(b0)
            ((-1, -7, -3, -1), 1)
            sage: m._solve_vector_linbox(b0, 'dense_elimination')
            ((-1, -7, -3, -1), 1)
            sage: m._solve_vector_linbox(b0, 'sparse_elimination')
            ((-1, -7, -3, -1), 1)
            sage: m._solve_vector_linbox(b0, 'wiedemann')
            ((-1, -7, -3, -1), 1)
            sage: m._solve_vector_linbox(b0, 'blackbox')
            ((-1, -7, -3, -1), 1)

            sage: b1 = vector((1,2,3,4))
            sage: m._solve_vector_linbox(b1)
            ((-18, -92, -41, -17), 5)
            sage: m._solve_vector_linbox(b1, 'dense_elimination')
            ((-18, -92, -41, -17), 5)
            sage: m._solve_vector_linbox(b1, 'sparse_elimination')
            ((-18, -92, -41, -17), 5)
            sage: m._solve_vector_linbox(b1, 'wiedemann')
            ((-18, -92, -41, -17), 5)
            sage: m._solve_vector_linbox(b1, 'blackbox')
            ((-18, -92, -41, -17), 5)

            sage: a1, d1 = m._solve_vector_linbox(b1)
            sage: d1 * b1 == m * a1
            True

        TESTS::

            sage: algos = ["default", "dense_elimination", "sparse_elimination",
            ....:          "blackbox", "wiedemann"]
            sage: for i in range(20):
            ....:     dim = randint(1, 30)
            ....:     M = MatrixSpace(ZZ, dim, sparse=True)
            ....:     density = min(1, 4/dim)
            ....:     m = M.random_element(density=density)
            ....:     while m.rank() != dim:
            ....:         m = M.random_element(density=density)
            ....:     U = m.column_space().dense_module()
            ....:     for algo in algos:
            ....:         u, d = m._solve_vector_linbox(U.zero(), algorithm=algo)
            ....:         assert u.is_zero()
            ....:         b = U.random_element()
            ....:         x, d = m._solve_vector_linbox(b, algorithm=algo)
            ....:         assert m * x == d * b
        """
        Vin = self.column_ambient_module(base_ring=None, sparse=False)
        v = Vin(v)

        if self._nrows == 0 or self._ncols == 0:
            raise ValueError("not implemented for nrows=0 or ncols=0")

        # LinBox "solve" is mostly broken for nonsquare or singular matrices.
        # The conditions below could be removed once all LinBox issues has
        # been solved.
        if self._nrows != self._ncols or self.rank() != self._nrows:
            raise ValueError("only available for full rank square matrices")

        cdef givaro.ZRing givZZ
        cdef linbox.SparseMatrix_integer * A = new_linbox_matrix_integer_sparse(givZZ, self)
        cdef linbox.DenseVector_integer * b = new_linbox_vector_integer_dense(givZZ, v)
        cdef linbox.DenseVector_integer * res = new linbox.DenseVector_integer(givZZ, <size_t> self._ncols)
        cdef givaro.Integer D

        method = get_method(algorithm)

        sig_on()
        if method == METHOD_DEFAULT:
            linbox.solve(res[0], D, A[0], b[0])
        elif method == METHOD_WIEDEMANN:
            linbox.solve(res[0], D, A[0], b[0], linbox.Method.Wiedemann())
        elif method == METHOD_DENSE_ELIMINATION:
            linbox.solve(res[0], D, A[0], b[0], linbox.Method.DenseElimination())
        elif method == METHOD_SPARSE_ELIMINATION:
            linbox.solve(res[0], D, A[0], b[0], linbox.Method.SparseElimination())
        elif method == METHOD_BLACKBOX:
            linbox.solve(res[0], D, A[0], b[0], linbox.Method.Blackbox())
        sig_off()

        Vout = self.row_ambient_module(base_ring=None, sparse=False)
        res_sage = new_sage_vector_integer_dense(Vout, res[0])
        cdef Integer d = PY_NEW(Integer)
        mpz_set(d.value, D.get_mpz_const())

        del A
        del b
        del res

        return (res_sage, d)

    def _solve_matrix_linbox(self, mat, algorithm=None):
        r"""
        Solve the equation ``A x = mat`` where ``A`` is this matrix.

        EXAMPLES::

            sage: m = matrix(ZZ, [[1,2],[1,0]], sparse=True)
            sage: b = matrix(ZZ, 2, 4, [1,0,2,0,1,1,2,0], sparse=False)
            sage: u, d = m._solve_matrix_linbox(b)
            sage: u
            [ 1  2  2  0]
            [ 0 -1  0  0]
            sage: m * u == b * diagonal_matrix(d)
            True

            sage: u, d = m._solve_matrix_linbox([[1,3,4],[0,1,0]])
            sage: u
            [0 1 0]
            [1 1 2]
            sage: d
            (2, 1, 1)

        Test input::

            sage: m = matrix(ZZ, [[1,2],[1,0]], sparse=True)
            sage: b = matrix(ZZ, 3, 3, range(9))
            sage: m._solve_matrix_linbox(b)
            Traceback (most recent call last):
            ...
            ValueError: wrong matrix dimension

            sage: m._solve_matrix_linbox([[1,1],[2,3]], algorithm='hop')
            Traceback (most recent call last):
            ...
            ValueError: unknown algorithm

        TESTS::

            sage: algos = ["default", "dense_elimination", "sparse_elimination",
            ....:          "blackbox", "wiedemann"]

            sage: for _ in range(10):
            ....:     dim = randint(2, 10)
            ....:     M = MatrixSpace(ZZ, dim, sparse=True)
            ....:     m = M.random_element(density=min(1,10/dim))
            ....:     while m.rank() != dim:
            ....:         m = M.random_element(density=min(1,10/dim))
            ....:     b = random_matrix(ZZ, dim, 7)
            ....:     Mb = b.parent()
            ....:     for algo in algos:
            ....:         u, d = m._solve_matrix_linbox(b, algo)
            ....:         assert m * u == b * diagonal_matrix(d)
        """
        if self._nrows == 0 or self._ncols == 0:
            raise ValueError("not implemented for nrows=0 or ncols=0")

        from sage.matrix.constructor import matrix
        from sage.modules.free_module_element import vector

        cdef Matrix_integer_dense B
        if not isinstance(mat, Matrix):
            B = <Matrix_integer_dense?> matrix(ZZ, mat, sparse=False)
        else:
            B = <Matrix_integer_dense?> mat.change_ring(ZZ).dense_matrix()
        if B._nrows != self._nrows:
            raise ValueError("wrong matrix dimension")

        # LinBox "solve" is mostly broken for singular matrices. The
        # conditions below could be removed once all LinBox issues
        # have been solved.
        if self._nrows != self._ncols or self.rank() != self._nrows:
            raise ValueError("only available for full rank square matrices")

        cdef givaro.ZRing givZZ
        cdef linbox.SparseMatrix_integer * A = new_linbox_matrix_integer_sparse(givZZ, self)
        cdef linbox.DenseVector_integer * b = new linbox.DenseVector_integer(givZZ, <size_t> self._nrows)
        cdef linbox.DenseVector_integer * res = new linbox.DenseVector_integer(givZZ, <size_t> self._ncols)
        cdef givaro.Integer D

        cdef int algo = get_method(algorithm)

        cdef Matrix_integer_dense X = matrix(ZZ, A.coldim(), B.ncols(), sparse=False)  # solution
        cdef Vector_integer_dense d = vector(ZZ, X.ncols(), sparse=False)  # multipliers

        sig_on()
        cdef size_t i, j
        for i in range(X.ncols()):
            # set b to the i-th column of B
            for j in range(A.coldim()):
                fmpz_get_mpz(<mpz_t> b.getEntry(j).get_mpz(), fmpz_mat_entry(B._matrix, j, i))

            # solve the current row
            if algo == METHOD_DEFAULT:
                linbox.solve(res[0], D, A[0], b[0])
            elif algo == METHOD_DENSE_ELIMINATION:
                linbox.solve(res[0], D, A[0], b[0], linbox.Method.DenseElimination())
            elif algo == METHOD_SPARSE_ELIMINATION:
                linbox.solve(res[0], D, A[0], b[0], linbox.Method.SparseElimination())
            elif algo == METHOD_BLACKBOX:
                linbox.solve(res[0], D, A[0], b[0], linbox.Method.Blackbox())
            elif algo == METHOD_WIEDEMANN:
                linbox.solve(res[0], D, A[0], b[0], linbox.Method.Wiedemann())

            # set i-th column of X to be res
            for j in range(A.coldim()):
                fmpz_set_mpz(fmpz_mat_entry(X._matrix, j, i), res[0].getEntry(j).get_mpz())

            # compute common gcd
            mpz_set(d._entries[i], D.get_mpz_const())
        sig_off()

        del A
        del b
        del res

        return X, d
