r"""
Base class for sparse matrices
"""
# ****************************************************************************
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

cimport cython
from cysignals.memory cimport check_allocarray, sig_free
from cysignals.signals cimport sig_check

cimport sage.matrix.matrix as matrix
cimport sage.matrix.matrix0 as matrix0
from sage.categories.rings import Rings
from sage.structure.element cimport Element, Vector
from sage.structure.richcmp cimport richcmp_item, rich_to_bool

from cpython cimport *
from cpython.object cimport Py_EQ, Py_NE

import sage.matrix.matrix_space


cdef class Matrix_sparse(matrix.Matrix):

    cdef bint is_sparse_c(self) noexcept:
        return 1

    cdef bint is_dense_c(self) noexcept:
        return 0

    def change_ring(self, ring):
        """
        Return the matrix obtained by coercing the entries of this matrix
        into the given ring.

        Always returns a copy (unless ``self`` is immutable, in which case
        returns ``self``).

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: A = matrix(QQ['x,y'], 2, [0,-1,2*x,-2], sparse=True); A
            [  0  -1]
            [2*x  -2]
            sage: A.change_ring(QQ['x,y,z'])
            [  0  -1]
            [2*x  -2]

        Subdivisions are preserved when changing rings::

            sage: A.subdivide([2],[]); A
            [  0  -1]
            [2*x  -2]
            [-------]
            sage: A.change_ring(RR['x,y'])
            [                 0  -1.00000000000000]
            [2.00000000000000*x  -2.00000000000000]
            [-------------------------------------]
        """
        if ring not in Rings():
            raise TypeError("input must be a ring")
        if ring is self._base_ring:
            if self._is_immutable:
                return self
            return self.__copy__()

        M = sage.matrix.matrix_space.MatrixSpace(ring, self._nrows, self._ncols, sparse=self.is_sparse_c())
        mat = M(self.dict(), coerce=True, copy=False)
        if self._subdivisions is not None:
            mat.subdivide(self.subdivisions())
        return mat

    def __copy__(self):
        """
        Return a copy of this matrix. Changing the entries of the copy will
        not change the entries of this matrix.

        EXAMPLES::

            sage: A = matrix(QQ['x,y'], 2, [0,-1,2,-2], sparse=True); A
            [ 0 -1]
            [ 2 -2]
            sage: B = copy(A); B
            [ 0 -1]
            [ 2 -2]
            sage: B is A
            False
            sage: B[0,0]=10; B
            [10 -1]
            [ 2 -2]
            sage: A
            [ 0 -1]
            [ 2 -2]
        """
        A = self.new_matrix(entries=self.dict(), coerce=False, copy=False)
        if self._subdivisions is not None:
            A.subdivide(*self.subdivisions())
        return A

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef long _hash_(self) except -1:
        """
        Return the hash of this matrix.

        Equal matrices should have equal hashes, even if one is sparse
        and the other is dense. We also ensure that zero matrices hash
        to zero and that scalar matrices have the same hash as the
        scalar.

        EXAMPLES::

            sage: m = matrix(2, range(6), sparse=True)
            sage: m.set_immutable()
            sage: hash(m)
            -154991009345361003  # 64-bit
            -2003358827          # 32-bit

        The sparse and dense hashes should agree::

            sage: d = m.dense_matrix()
            sage: d.set_immutable()
            sage: hash(d) == hash(m)
            True

        ::

            sage: A = Matrix(ZZ[['t']], 2, 2, range(4), sparse=True)
            sage: hash(A)
            Traceback (most recent call last):
            ...
            TypeError: mutable matrices are unhashable
            sage: A.set_immutable()
            sage: B = A.__copy__(); B.set_immutable()
            sage: hash(A) == hash(B)
            True

        TESTS::

            sage: R.<x> = ZZ[]
            sage: M = matrix(R, 10, 20, sparse=True); M.set_immutable()
            sage: hash(M)
            0
            sage: M = matrix(R, 10, 10, x, sparse=True); M.set_immutable()
            sage: hash(M) == hash(x)
            True
        """
        cdef dict D = self._dict()
        cdef long C[5]
        self.get_hash_constants(C)

        cdef long h = 0, k, l
        cdef Py_ssize_t i, j
        for ij, x in D.iteritems():
            sig_check()
            i = (<tuple>ij)[0]
            j = (<tuple>ij)[1]
            k = C[0] if i == 0 else C[1] + C[2] * i
            l = C[3] * (i - j) * (i ^ j)
            h += (k ^ l) * hash(x)
        h *= C[4]

        if h == -1:
            return -2
        return h

    def _multiply_classical(Matrix_sparse left, Matrix_sparse right):
        """
        EXAMPLES::

            sage: A = matrix(QQ['x,y'], 2, [0,-1,2,-2], sparse=True)
            sage: type(A)
            <class 'sage.matrix.matrix_generic_sparse.Matrix_generic_sparse'>
            sage: B = matrix(QQ['x,y'], 2, [-1,-1,-2,-2], sparse=True)
            sage: A * B
            [2 2]
            [2 2]
            sage: B * A
            [-2  3]
            [-4  6]
        """
        cdef Py_ssize_t row, col, row_start, k1, k2, len_left, len_right, a, b
        cdef list left_nonzero = <list> left.nonzero_positions(copy=False, column_order=False)
        cdef list right_nonzero = <list> right.nonzero_positions(copy=False, column_order=True)
        len_left = len(left_nonzero)
        len_right = len(right_nonzero)

        cdef dict e = {}
        k1 = 0
        while k1 < len_left:
            row_start = k1
            row = get_ij(left_nonzero, row_start, 0)
            k2 = 0
            while k2 < len_right:
                sig_check()
                col = get_ij(right_nonzero, k2, 1)
                s = None
                k1 = row_start
                while (k1 < len_left and get_ij(left_nonzero,k1,0) == row
                       and k2 < len_right and get_ij(right_nonzero,k2,1) == col):
                    sig_check()
                    a = get_ij(left_nonzero, k1, 1)
                    b = get_ij(right_nonzero, k2, 0)
                    if a == b:
                        if s is None:
                            s = left.get_unsafe(row,a) * right.get_unsafe(a,col)
                        else:
                            s += left.get_unsafe(row,a) * right.get_unsafe(a,col)
                        k1 += 1
                        k2 += 1
                    elif a < b:
                        k1 += 1
                    else:
                        k2 += 1
                if s is not None:
                    e[row, col] = s
                while k2 < len_right and get_ij(right_nonzero, k2, 1) == col:
                    k2 += 1
            while k1 < len_left and get_ij(left_nonzero, k1, 0) == row:
                k1 += 1
        return left.new_matrix(left._nrows, right._ncols, entries=e, coerce=False, copy=False)

    def _multiply_classical_with_cache(Matrix_sparse left, Matrix_sparse right):
        """
        This function computes the locations of the end of the rows/columns
        in the nonzero entries list once O(rows+cols) time and space, then
        uses these values in the inner loops. For large matrices this can
        be a 2x or more speedup, but the matrices can no longer be
        arbitrarily large as the runtime and space requirements are no
        longer functions of the total number of entries only.

        EXAMPLES::

            sage: A = matrix(QQ['x,y'], 2, [0,-1,2,-2], sparse=True)
            sage: type(A)
            <class 'sage.matrix.matrix_generic_sparse.Matrix_generic_sparse'>
            sage: B = matrix(QQ['x,y'], 2, [-1,-1,-2,-2], sparse=True)
            sage: A._multiply_classical_with_cache(B)
            [2 2]
            [2 2]
        """
        cdef Py_ssize_t row, col, row_start, k1, k2, len_left, len_right, a, b, i
        cdef Py_ssize_t* next_row
        cdef Py_ssize_t* next_col
        left_nonzero = left.nonzero_positions(copy=False, column_order=False)
        right_nonzero = right.nonzero_positions(copy=False, column_order=True)
        len_left = len(left_nonzero)
        len_right = len(right_nonzero)
        next_row = <Py_ssize_t *>check_allocarray(left._nrows, sizeof(Py_ssize_t))
        next_col = <Py_ssize_t *>check_allocarray(right._ncols, sizeof(Py_ssize_t))

        i = len_left - 1
        for row from left._nrows > row >= 0:
            next_row[row] = i + 1
            while i >= 0 and get_ij(left_nonzero,i,0) == row:
                i = i - 1
        i = len_right - 1
        for col from right._ncols > col >= 0:
            next_col[col] = i + 1
            while i >= 0 and get_ij(right_nonzero,i,1) == col:
                i = i - 1

        e = {}
        k1 = 0
        while k1 < len_left:
            row_start = k1
            row = get_ij(left_nonzero, row_start, 0)
            k2 = 0
            while k2 < len_right:
                sig_check()
                col = get_ij(right_nonzero, k2, 1)
                sum = None
                k1 = row_start
                while k1 < next_row[row] and k2 < next_col[col]:
                    sig_check()
                    a = get_ij(left_nonzero, k1,1)
                    b = get_ij(right_nonzero,k2,0)
                    if a == b:
                        if sum is None:
                            sum = left.get_unsafe(row,a)*right.get_unsafe(a,col)
                        else:
                            sum = sum + left.get_unsafe(row,a)*right.get_unsafe(a,col)
                        k1 = k1 + 1
                        k2 = k2 + 1
                    elif a < b:
                        k1 = k1 + 1
                    else:
                        k2 = k2 + 1
                if sum is not None:
                    e[row, col] = sum
                k2 = next_col[col]
            k1 = next_row[row]

        sig_free(next_row)
        sig_free(next_col)

        return left.new_matrix(left._nrows, right._ncols, entries=e, coerce=False, copy=False)

    cpdef _lmul_(self, Element right):
        """
        Left scalar multiplication. Internal usage only.

        INPUT:

        - ``right`` -- a ring element which must already be in the base ring
          of ``self`` (no coercion done here)

        OUTPUT: the matrix ``self * right``

        EXAMPLES::

            sage: M = Matrix(QQ, 3, 6, range(18), sparse=true); M
            [ 0  1  2  3  4  5]
            [ 6  7  8  9 10 11]
            [12 13 14 15 16 17]
            sage: (2/3)*M
            [   0  2/3  4/3    2  8/3 10/3]
            [   4 14/3 16/3    6 20/3 22/3]
            [   8 26/3 28/3   10 32/3 34/3]
            sage: 7*M
            [  0   7  14  21  28  35]
            [ 42  49  56  63  70  77]
            [ 84  91  98 105 112 119]
            sage: (1/4)*M
            [   0  1/4  1/2  3/4    1  5/4]
            [ 3/2  7/4    2  9/4  5/2 11/4]
            [   3 13/4  7/2 15/4    4 17/4]

        Really Large Example you would not want to do with normal matrices::

            sage: M = MatrixSpace(QQ, 100000, 1000000, sparse=True)
            sage: m = M.random_element(density=1/100000000)
            sage: m == (97/42)*(42/97*m)
            True
        """
        cdef Py_ssize_t k, r, c
        cdef Matrix_sparse M
        nc, nr = self.ncols(), self.nrows()
        M = self.new_matrix(nr, nc, copy=False, coerce=False)
        nz = self.nonzero_positions(copy=False)
        for k from 0 <= k < len(nz):
            r = get_ij(nz, k, 0)
            c = get_ij(nz, k, 1)
            entry = self.get_unsafe(r,c)*right
            M.set_unsafe(r, c, entry)
        return M

    cdef bint _will_use_strassen(self, matrix0.Matrix right) except -2:
        # never use Strassen for sparse matrix multiply
        return 0

    def _pickle(self):
        version = -1
        data = self._dict()  # dict of all elements
        return data, version

    def _unpickle_generic(self, data, int version):
        cdef Py_ssize_t i, j
        if version == -1:
            for (i, j), x in data.iteritems():
                self.set_unsafe(i, j, x)
        else:
            raise RuntimeError("unknown matrix version (=%s)" % version)

    cpdef _richcmp_(self, right, int op):
        """
        Rich comparison.

        EXAMPLES::

            sage: M = matrix({(5,5): 2})
            sage: Mp = matrix({(5,5): 7, (3,1):-2})
            sage: M > Mp
            True
            sage: M == M.transpose()
            True
            sage: M != Mp
            True
        """
        other = <Matrix_sparse>right
        if op == Py_EQ:
            return self._dict() == other._dict()
        if op == Py_NE:
            return self._dict() != other._dict()
        cdef Py_ssize_t i, j
        # Parents are equal, so dimensions of self and other are equal
        for i in range(self._nrows):
            for j in range(self._ncols):
                lij = self.get_unsafe(i, j)
                rij = other.get_unsafe(i, j)
                r = richcmp_item(lij, rij, op)
                if r is not NotImplemented:
                    return bool(r)
        # Matrices are equal
        return rich_to_bool(op, 0)

    def transpose(self):
        """
        Return the transpose of ``self``, without changing ``self``.

        EXAMPLES: We create a matrix, compute its transpose, and note that
        the original matrix is not changed.

        ::

            sage: M = MatrixSpace(QQ, 2, sparse=True)
            sage: A = M([1,2,3,4]); A
            [1 2]
            [3 4]
            sage: B = A.transpose(); B
            [1 3]
            [2 4]

        ``.T`` is a convenient shortcut for the transpose::

           sage: A.T
           [1 3]
           [2 4]

        .. SEEALSO:: :meth:`antitranspose`
        """
        cdef Matrix_sparse A
        A = self.new_matrix(self._ncols, self._nrows)

        nz = self.nonzero_positions(copy=False)
        cdef Py_ssize_t i, j, k
        for k from 0 <= k < len(nz):
            i = get_ij(nz, k, 0)
            j = get_ij(nz, k, 1)
            A.set_unsafe(j,i,self.get_unsafe(i,j))
        if self._subdivisions is not None:
            row_divs, col_divs = self.subdivisions()
            A.subdivide(col_divs, row_divs)
        return A

    def antitranspose(self):
        """
        Return the antitranspose of ``self``, without changing ``self``.

        This is the mirror image along the other diagonal.

        EXAMPLES::

            sage: M = MatrixSpace(QQ, 2, sparse=True)
            sage: A = M([1,2,3,4]); A
            [1 2]
            [3 4]
            sage: A.antitranspose()
            [4 2]
            [3 1]

        .. SEEALSO:: :meth:`transpose`
        """
        cdef Matrix_sparse A
        A = self.new_matrix(self._ncols, self._nrows)

        nz = self.nonzero_positions(copy=False)
        cdef Py_ssize_t i, j, k
        for k from 0 <= k < len(nz):
            i = get_ij(nz, k, 0)
            j = get_ij(nz, k, 1)
            A.set_unsafe(self._ncols-j-1, self._nrows-i-1,self.get_unsafe(i,j))
        if self._subdivisions is not None:
            row_divs, col_divs = self.subdivisions()
            A.subdivide(list(reversed([self._ncols - t for t in col_divs])),
                            list(reversed([self._nrows - t for t in row_divs])))
        return A

    def _reverse_unsafe(self):
        r"""
        TESTS::

            sage: m = matrix(QQ, 3, 3, {(2,2): 1}, sparse=True)
            sage: m._reverse_unsafe()
            sage: m
            [1 0 0]
            [0 0 0]
            [0 0 0]
            sage: m = matrix(QQ, 3, 3, {(2,2): 1, (0,0):2}, sparse=True)
            sage: m._reverse_unsafe()
            sage: m
            [1 0 0]
            [0 0 0]
            [0 0 2]
            sage: m = matrix(QQ, 3, 3, {(1,2): 1}, sparse=True)
            sage: m._reverse_unsafe()
            sage: m
            [0 0 0]
            [1 0 0]
            [0 0 0]
            sage: m = matrix(QQ, 3, 3, {(1,1): 1}, sparse=True)
            sage: m._reverse_unsafe()
            sage: m
            [0 0 0]
            [0 1 0]
            [0 0 0]
        """
        cdef Py_ssize_t i, j, ii, jj
        for i,j in self.nonzero_positions(copy=False):
            ii = self._nrows - i - 1
            jj = self._ncols - j - 1
            if (i > ii or (i == ii and j >= jj)) and self.get_unsafe(ii, jj):
                # already swapped
                continue

            e1 = self.get_unsafe(i, j)
            e2 = self.get_unsafe(ii, jj)
            self.set_unsafe(i, j, e2)
            self.set_unsafe(ii, jj, e1)

    def charpoly(self, var='x', **kwds):
        """
        Return the characteristic polynomial of this matrix.

        .. NOTE::

            the generic sparse charpoly implementation in Sage is to
            just compute the charpoly of the corresponding dense
            matrix, so this could use a lot of memory. In particular,
            for this matrix, the charpoly will be computed using a
            dense algorithm.

        EXAMPLES::

            sage: A = matrix(ZZ, 4, range(16), sparse=True)
            sage: A.charpoly()
            x^4 - 30*x^3 - 80*x^2
            sage: A.charpoly('y')
            y^4 - 30*y^3 - 80*y^2
            sage: A.charpoly()
            x^4 - 30*x^3 - 80*x^2
        """
        f = self.fetch('charpoly')
        if f is not None:
            return f.change_variable_name(var)
        f = self.dense_matrix().charpoly(var=var, **kwds)
        self.cache('charpoly', f)
        return f

    def determinant(self, **kwds):
        """
        Return the determinant of this matrix.

        .. NOTE::

            the generic sparse determinant implementation in Sage is
            to just compute the determinant of the corresponding dense
            matrix, so this could use a lot of memory. In particular,
            for this matrix, the determinant will be computed using a
            dense algorithm.

        EXAMPLES::

            sage: A = matrix(ZZ, 4, range(16), sparse=True)
            sage: B = A + identity_matrix(ZZ, 4, sparse=True)
            sage: B.det()
            -49
        """
        d = self.fetch('det')
        if d is not None:
            return d
        d = self.dense_matrix().determinant(**kwds)
        self.cache('det', d)
        return d

    def _elementwise_product(self, right):
        r"""
        Return the elementwise product of two sparse
        matrices with identical base rings.

        This routine assumes that ``self`` and ``right``
        are both matrices, both sparse, with identical
        sizes and with identical base rings.  It is
        "unsafe" in the sense that these conditions
        are not checked and no sensible errors are
        raised.

        This routine is meant to be called from the
        :meth:`~sage.matrix.matrix2.Matrix.elementwise_product`
        method, which will ensure that this routine receives
        proper input.  More thorough documentation is provided
        there.

        EXAMPLES::

            sage: A = matrix(ZZ, 2, range(6), sparse=True)
            sage: B = matrix(ZZ, 2, [1,0,2,0,3,0], sparse=True)
            sage: A._elementwise_product(B)
            [ 0  0  4]
            [ 0 12  0]

        AUTHOR:

        - Rob Beezer (2009-07-14)
        """
        cdef Py_ssize_t k, r, c
        cdef Matrix_sparse other, prod

        nc, nr = self.ncols(), self.nrows()
        other = right
        prod = self.new_matrix(nr, nc, copy=False, coerce=False)
        nzleft = self.nonzero_positions(copy=False)
        nzright = other.nonzero_positions(copy=False)
        for k from 0 <= k < len(nzleft):
            r = get_ij(nzleft, k, 0)
            c = get_ij(nzleft, k, 1)
            if (r,c) in nzright:
                entry = self.get_unsafe(r,c)*other.get_unsafe(r,c)
                prod.set_unsafe(r,c,entry)
        return prod

    def apply_morphism(self, phi):
        """
        Apply the morphism ``phi`` to the coefficients of this sparse matrix.

        The resulting matrix is over the codomain of ``phi``.

        INPUT:

        - ``phi`` -- a morphism, so ``phi`` is callable and
           ``phi.domain()`` and ``phi.codomain()`` are defined. The
           codomain must be a ring.

        OUTPUT: a matrix over the codomain of ``phi``

        EXAMPLES::

            sage: m = matrix(ZZ, 3, range(9), sparse=True)
            sage: phi = ZZ.hom(GF(5))
            sage: m.apply_morphism(phi)
            [0 1 2]
            [3 4 0]
            [1 2 3]
            sage: m.apply_morphism(phi).parent()
            Full MatrixSpace of 3 by 3 sparse matrices
             over Finite Field of size 5
        """
        R = phi.codomain()
        M = sage.matrix.matrix_space.MatrixSpace(R, self._nrows,
                                                 self._ncols, sparse=True)
        return M({ij: phi(z) for ij, z in self.dict().iteritems()})

    def apply_map(self, phi, R=None, sparse=True):
        r"""
        Apply the given map ``phi`` (an arbitrary Python function or callable
        object) to this matrix.

        If ``R`` is not given, automatically determine the base ring
        of the resulting matrix.

        INPUT:

        - ``phi`` -- arbitrary Python function or callable object

        - ``R`` -- (optional) ring

        - ``sparse`` -- boolean (default: ``True``); whether to return
          a sparse or a dense matrix

        OUTPUT: a matrix over ``R``

        EXAMPLES::

            sage: m = matrix(ZZ, 10000, {(1,2): 17}, sparse=True)

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(9)
            sage: f = lambda x: k(x)
            sage: n = m.apply_map(f)
            sage: n.parent()
            Full MatrixSpace of 10000 by 10000 sparse matrices
             over Finite Field in a of size 3^2
            sage: n[1, 2]
            2

        An example where the codomain is explicitly specified.

        ::

            sage: n = m.apply_map(lambda x: x%3, GF(3))
            sage: n.parent()
            Full MatrixSpace of 10000 by 10000 sparse matrices
             over Finite Field of size 3
            sage: n[1, 2]
            2

        If we did not specify the codomain, the resulting matrix in the
        above case ends up over `\ZZ` again::

            sage: n = m.apply_map(lambda x: x%3)
            sage: n.parent()
            Full MatrixSpace of 10000 by 10000 sparse matrices over Integer Ring
            sage: n[1, 2]
            2

        If ``self`` is subdivided, the result will be as well::

            sage: m = matrix(2, 2, [0, 0, 3, 0])
            sage: m.subdivide(None, 1); m
            [0|0]
            [3|0]
            sage: m.apply_map(lambda x: x*x)
            [0|0]
            [9|0]

        If the map sends zero to a nonzero value, then it may be useful to
        get the result as a dense matrix.

        ::

            sage: m = matrix(ZZ, 3, 3, [0] * 7 + [1,2], sparse=True); m
            [0 0 0]
            [0 0 0]
            [0 1 2]
            sage: parent(m)
            Full MatrixSpace of 3 by 3 sparse matrices over Integer Ring
            sage: n = m.apply_map(lambda x: x+polygen(QQ), sparse=False); n
            [    x     x     x]
            [    x     x     x]
            [    x x + 1 x + 2]
            sage: parent(n)
            Full MatrixSpace of 3 by 3 dense matrices over Univariate Polynomial Ring in x over Rational Field

        TESTS::

            sage: m = matrix([], sparse=True)
            sage: m.apply_map(lambda x: x*x) == m
            True

            sage: m.apply_map(lambda x: x*x, sparse=False).parent()
            Full MatrixSpace of 0 by 0 dense matrices over Integer Ring

        Check that we do not unnecessarily apply phi to 0 in the sparse case::

            sage: m = matrix(QQ, 2, 2, range(1, 5), sparse=True)
            sage: m.apply_map(lambda x: 1/x)
            [  1 1/2]
            [1/3 1/4]

        Test subdivisions when phi maps 0 to nonzero::

            sage: m = matrix(2, 2, [0, 0, 3, 0])
            sage: m.subdivide(None, 1); m
            [0|0]
            [3|0]
            sage: m.apply_map(lambda x: x+1)
            [1|1]
            [4|1]

        When applying a map to a sparse zero matrix, the codomain is determined
        from the image of zero (:issue:`29214`)::

            sage: matrix(RR, 2, 2, sparse=True).apply_map(floor).base_ring() is ZZ
            True
        """
        if self._nrows==0 or self._ncols==0:
            if not sparse:
                return self.dense_matrix()
            else:
                return self.__copy__()
        self_dict = self._dict()
        if len(self_dict) < self._nrows * self._ncols:
            zero_res = phi(self.base_ring()(0))
        else:
            zero_res = None
        v = [(ij, phi(z)) for ij,z in self_dict.iteritems()]
        if R is None:
            w = [x for _, x in v]
            if zero_res is not None:
                w.append(zero_res)
            w = sage.structure.sequence.Sequence(w)
            R = w.universe()
            v = {v[i][0]: w[i] for i in range(len(v))}
        else:
            v = dict(v)
        if zero_res is not None and not zero_res.is_zero():
            M = sage.matrix.matrix_space.MatrixSpace(R, self._nrows,
                                                     self._ncols, sparse=sparse)
            m = M([zero_res] * (self._nrows * self._ncols))
            for i,n in v.items():
                m[i] = n
            if self._subdivisions is not None:
                m.subdivide(*self.subdivisions())
            return m

        M = sage.matrix.matrix_space.MatrixSpace(R, self._nrows,
                   self._ncols, sparse=sparse)
        m = M(v)
        if self._subdivisions is not None:
            m.subdivide(*self.subdivisions())
        return m

    def _derivative(self, var=None, R=None):
        """
        Differentiate with respect to ``var`` by differentiating each element
        with respect to ``var``.

        .. SEEALSO::

           :meth:`derivative`

        EXAMPLES::

            sage: m = matrix(2, [x^i for i in range(4)], sparse=True)                   # needs sage.symbolic
            sage: m._derivative(x)                                                      # needs sage.symbolic
            [    0     1]
            [  2*x 3*x^2]
        """
        # We would just use apply_map, except that Cython does not
        # allow lambda functions

        if self._nrows==0 or self._ncols==0:
            return self.__copy__()
        v = [(ij, z.derivative(var)) for ij, z in self.dict().iteritems()]
        if R is None:
            w = [x for _, x in v]
            w = sage.structure.sequence.Sequence(w)
            R = w.universe()
            v = {v[i][0]: w[i] for i in range(len(v))}
        else:
            v = dict(v)
        M = sage.matrix.matrix_space.MatrixSpace(R, self._nrows,
                   self._ncols, sparse=True)
        return M(v)

    def density(self):
        """
        Return the density of the matrix.

        By density we understand the ratio of the number of nonzero
        positions and the number ``self.nrows() * self.ncols()``,
        i.e. the number of possible nonzero positions.

        EXAMPLES::

            sage: a = matrix([[],[],[],[]], sparse=True); a.density()
            0
            sage: a = matrix(5000,5000,{(1,2): 1}); a.density()
            1/25000000
        """
        nr = self.nrows()
        nc = self.ncols()
        if nc == 0 or nr == 0:
            return 0
        from sage.rings.rational_field import QQ
        return QQ(len(self.nonzero_positions(copy=False))) / (nr * nc)

    def matrix_from_rows_and_columns(self, rows, columns):
        """
        Return the matrix constructed from ``self`` from the given rows and
        columns.

        EXAMPLES::

            sage: M = MatrixSpace(Integers(8),3,3, sparse=True)
            sage: A = M(range(9)); A
            [0 1 2]
            [3 4 5]
            [6 7 0]
            sage: A.matrix_from_rows_and_columns([1], [0,2])
            [3 5]
            sage: A.matrix_from_rows_and_columns([1,2], [1,2])
            [4 5]
            [7 0]

        Note that row and column indices can be reordered or repeated::

            sage: A.matrix_from_rows_and_columns([2,1], [2,1])
            [0 7]
            [5 4]

        For example here we take from row 1 columns 2 then 0 twice, and do
        this 3 times.

        ::

            sage: A.matrix_from_rows_and_columns([1,1,1],[2,0,0])
            [5 3 3]
            [5 3 3]
            [5 3 3]

        We can efficiently extract large submatrices::

            sage: A = random_matrix(ZZ, 100000, density=.00005, sparse=True)  # long time (4s on sage.math, 2012)
            sage: B = A[50000:,:50000]        # long time
            sage: count = 0
            sage: for i, j in A.nonzero_positions():  # long time
            ....:     if i >= 50000 and j < 50000:
            ....:         assert B[i-50000, j] == A[i, j]
            ....:         count += 1
            sage: count == sum(1 for _ in B.nonzero_positions())  # long time
            True

        We must pass in a list of indices::

            sage: A = random_matrix(ZZ,100,density=.02,sparse=True)
            sage: A.matrix_from_rows_and_columns(1,[2,3])
            Traceback (most recent call last):
            ...
            TypeError: 'sage.rings.integer.Integer' object is not iterable

            sage: A.matrix_from_rows_and_columns([1,2],3)
            Traceback (most recent call last):
            ...
            TypeError: 'sage.rings.integer.Integer' object is not iterable

        AUTHORS:

        - Jaap Spies (2006-02-18)

        - Didier Deshommes: some Pyrex speedups implemented

        - Jason Grout: sparse matrix optimizations
        """
        if not isinstance(rows, (list, tuple)):
            rows = list(rows)

        if not isinstance(columns, (list, tuple)):
            columns = list(columns)

        cdef Py_ssize_t nrows, ncols, k, i, j

        ncols = PyList_GET_SIZE(columns)
        nrows = PyList_GET_SIZE(rows)
        cdef Matrix_sparse A = self.new_matrix(nrows = nrows, ncols = ncols)

        tmp = [el for el in columns if 0 <= el < self._ncols]
        columns = tmp
        if ncols != PyList_GET_SIZE(columns):
            raise IndexError("column index out of range")

        tmp = [el for el in rows if 0 <= el < self._nrows]
        rows = tmp
        if nrows != PyList_GET_SIZE(rows):
            raise IndexError("row index out of range")

        row_map = {}
        for new_row, old_row in enumerate(rows):
            if old_row in row_map:
                row_map[old_row].append(new_row)
            else:
                row_map[old_row] = [new_row]

        col_map = {}
        for new_col, old_col in enumerate(columns):
            if old_col in col_map:
                col_map[old_col].append(new_col)
            else:
                col_map[old_col] = [new_col]

        nz = self.nonzero_positions(copy=False)
        for k in range(len(nz)):
            i = get_ij(nz, k, 0)
            j = get_ij(nz, k, 1)
            if i in row_map and j in col_map:
                entry = self.get_unsafe(i, j)
                for new_row in row_map[i]:
                    for new_col in col_map[j]:
                        A.set_unsafe(new_row, new_col, entry)
        return A

    cdef _stack_impl(self, bottom):
        r"""
        Stack ``self`` on top of ``bottom``::

            [ self  ]
            [ bottom ]

        EXAMPLES::

            sage: M = Matrix(QQ, 2, 3, range(6), sparse=True)
            sage: N = Matrix(QQ, 1, 3, [10,11,12], sparse=True)
            sage: M.stack(N)
            [ 0  1  2]
            [ 3  4  5]
            [10 11 12]

        A vector may be stacked below a matrix. ::

            sage: A = matrix(QQ, 2, 5, range(10), sparse=True)
            sage: v = vector(QQ, 5, range(5), sparse=True)
            sage: A.stack(v)
            [0 1 2 3 4]
            [5 6 7 8 9]
            [0 1 2 3 4]

        The ``subdivide`` option will add a natural subdivision between
        ``self`` and ``other``.  For more details about how subdivisions
        are managed when stacking, see
        :meth:`sage.matrix.matrix1.Matrix.stack`.  ::

            sage: A = matrix(ZZ, 3, 4, range(12), sparse=True)
            sage: B = matrix(ZZ, 2, 4, range(8), sparse=True)
            sage: A.stack(B, subdivide=True)
            [ 0  1  2  3]
            [ 4  5  6  7]
            [ 8  9 10 11]
            [-----------]
            [ 0  1  2  3]
            [ 4  5  6  7]

        TESTS:

        One can stack matrices over different rings (:issue:`16399`). ::

            sage: M = Matrix(ZZ, 2, 3, range(6), sparse=True)
            sage: N = Matrix(QQ, 1, 3, [10,11,12], sparse=True)
            sage: M.stack(N)
            [ 0  1  2]
            [ 3  4  5]
            [10 11 12]
            sage: N.stack(M)
            [10 11 12]
            [ 0  1  2]
            [ 3  4  5]
            sage: M2 = Matrix(ZZ['x'], 2, 3, range(6), sparse=True)
            sage: N.stack(M2)
            [10 11 12]
            [ 0  1  2]
            [ 3  4  5]
        """
        cdef Matrix_sparse other = <Matrix_sparse>bottom
        cdef Matrix_sparse Z
        Z = self.new_matrix(nrows=self._nrows + other._nrows, ncols=self._ncols)

        for i, j in self.nonzero_positions(copy=False):
            Z.set_unsafe(i, j, self.get_unsafe(i,j))
        for i, j in other.nonzero_positions(copy=False):
            Z.set_unsafe(i + self._nrows, j, other.get_unsafe(i,j))

        return Z

    def augment(self, right, subdivide=False):
        r"""
        Return the augmented matrix of the form::

            [self | right].

        EXAMPLES::

            sage: M = MatrixSpace(QQ, 2, 2, sparse=True)
            sage: A = M([1,2, 3,4])
            sage: A
            [1 2]
            [3 4]
            sage: N = MatrixSpace(QQ, 2, 1, sparse=True)
            sage: B = N([9,8])
            sage: B
            [9]
            [8]
            sage: A.augment(B)
            [1 2 9]
            [3 4 8]
            sage: B.augment(A)
            [9 1 2]
            [8 3 4]

        A vector may be augmented to a matrix. ::

            sage: A = matrix(QQ, 3, 4, range(12), sparse=True)
            sage: v = vector(QQ, 3, range(3), sparse=True)
            sage: A.augment(v)
            [ 0  1  2  3  0]
            [ 4  5  6  7  1]
            [ 8  9 10 11  2]

        The ``subdivide`` option will add a natural subdivision between
        ``self`` and ``right``.  For more details about how subdivisions
        are managed when augmenting, see
        :meth:`sage.matrix.matrix1.Matrix.augment`.  ::

            sage: A = matrix(QQ, 3, 5, range(15), sparse=True)
            sage: B = matrix(QQ, 3, 3, range(9), sparse=True)
            sage: A.augment(B, subdivide=True)
            [ 0  1  2  3  4| 0  1  2]
            [ 5  6  7  8  9| 3  4  5]
            [10 11 12 13 14| 6  7  8]

        TESTS:

        Verify that :issue:`12689` is fixed::

            sage: A = identity_matrix(QQ, 2, sparse=True)
            sage: B = identity_matrix(ZZ, 2, sparse=True)
            sage: A.augment(B)
            [1 0 1 0]
            [0 1 0 1]
        """
        if not isinstance(right, matrix.Matrix):
            if hasattr(right, '_vector_'):
                right = right.column()
            else:
                raise TypeError("right must be a matrix")

        if not (self._base_ring is right.base_ring()):
            right = right.change_ring(self._base_ring)

        cdef Matrix_sparse other = right.sparse_matrix()

        if self._nrows != other._nrows:
            raise TypeError("number of rows must be the same")

        cdef Matrix_sparse Z
        Z = self.new_matrix(ncols = self._ncols + other._ncols)
        for i, j in self.nonzero_positions(copy=False):
            Z.set_unsafe(i, j, self.get_unsafe(i,j))
        for i, j in other.nonzero_positions(copy=False):
            Z.set_unsafe(i, j + self._ncols, other.get_unsafe(i,j))
        if subdivide:
            Z._subdivide_on_augment(self, other)
        return Z

    cdef _vector_times_matrix_(self, Vector v):
        """
        Return the vector times matrix product.

        INPUT:

        - ``v`` -- a free module element

        OUTPUT: the vector times matrix product ``v*A``

        EXAMPLES::

            sage: v = FreeModule(ZZ, 3)([1, 2, 3])
            sage: m = matrix(QQ, 3, 4, range(12), sparse=True)
            sage: v * m
            (32, 38, 44, 50)

        TESTS::

            sage: (v * m).is_sparse()
            True
            sage: (v * m).parent() is m.row(0).parent()
            True
            """
        cdef int i, j
        if self._nrows != v._degree:
            raise ArithmeticError("number of rows of matrix must equal degree of vector")
        parent = self.row_ambient_module(base_ring=None, sparse=v.is_sparse_c())
        s = parent.zero_vector()
        for (i, j), a in self._dict().iteritems():
            s[j] += v[i] * a
        return s

    cdef _matrix_times_vector_(self, Vector v):
        """
        Return the matrix times vector product.

        INPUT:

        - ``v`` -- a free module element

        OUTPUT: the matrix times vector product ``A*v``

        EXAMPLES::

            sage: v = FreeModule(ZZ, 3)([1, 2, 3])
            sage: m = matrix(QQ, 4, 3, range(12), sparse=True)
            sage: m * v
            (8, 26, 44, 62)

        TESTS::

            sage: (m * v).is_sparse()
            True
            sage: (m * v).parent() is m.column(0).parent()
            True

        Check that the bug in :issue:`13854` has been fixed::

            sage: # needs sage.combinat sage.libs.singular
            sage: A.<x,y> = FreeAlgebra(QQ, 2)
            sage: P.<x,y> = A.g_algebra(relations={y*x: -x*y}, order='lex')
            sage: M = Matrix([[x]], sparse=True)
            sage: w = vector([y])
            doctest:...: UserWarning: You are constructing a free module
            over a noncommutative ring. Sage does not have a concept
            of left/right and both sided modules, so be careful.
            It's also not guaranteed that all multiplications are
            done from the right side.
            doctest:...: UserWarning: You are constructing a free module
            over a noncommutative ring. Sage does not have a concept
            of left/right and both sided modules, so be careful.
            It's also not guaranteed that all multiplications are
            done from the right side.
            sage: M*w
            (x*y)
        """
        cdef int i, j
        if self._ncols != v._degree:
            raise ArithmeticError("number of columns of matrix must equal degree of vector")
        parent = self.column_ambient_module(base_ring=None, sparse=v.is_sparse_c())
        s = parent.zero_vector()
        for (i, j), a in self._dict().iteritems():
            s[i] += a * v[j]
        return s


@cython.boundscheck(False)
@cython.wraparound(False)
# Return v[i][j] where v is a list of tuples.
# No checking is done, make sure you feed it valid input!
cdef inline Py_ssize_t get_ij(v, Py_ssize_t i, Py_ssize_t j) noexcept:
    t = (<list>v)[i]
    return (<tuple>t)[j]
