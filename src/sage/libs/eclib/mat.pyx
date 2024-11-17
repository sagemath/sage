"""
Cremona matrices
"""

from sage.libs.eclib cimport scalar, addscalar

from sage.matrix.matrix_space import MatrixSpace
from sage.rings.integer_ring import ZZ

from sage.matrix.matrix_integer_sparse cimport Matrix_integer_sparse
from sage.matrix.matrix_integer_dense cimport Matrix_integer_dense
from sage.rings.integer cimport Integer


cdef class Matrix:
    """
    A Cremona Matrix.

    EXAMPLES::

        sage: M = CremonaModularSymbols(225)
        sage: t = M.hecke_matrix(2)
        sage: type(t)
        <class 'sage.libs.eclib.mat.Matrix'>
        sage: t
        61 x 61 Cremona matrix over Rational Field

    TESTS::

        sage: t = CremonaModularSymbols(11).hecke_matrix(2); t
        3 x 3 Cremona matrix over Rational Field
        sage: type(t)
        <class 'sage.libs.eclib.mat.Matrix'>
    """
    def __repr__(self):
        """
        String representation of this matrix.  Use print(self.str()) to
        print out the matrix entries on the screen.

        EXAMPLES::

            sage: M = CremonaModularSymbols(23)
            sage: t = M.hecke_matrix(2); t
            5 x 5 Cremona matrix over Rational Field
            sage: print(t.str())
            [ 3  0  0  0  0]
            [-1 -1  0  0 -1]
            [ 1  1  0  1  1]
            [-1  1  1 -1  0]
            [ 0 -1  0  0  0]
        """
        return "%s x %s Cremona matrix over Rational Field" % (self.nrows(), self.ncols())

    def str(self):
        r"""
        Return full string representation of this matrix, never in compact form.

        EXAMPLES::

            sage: M = CremonaModularSymbols(22, sign=1)
            sage: t = M.hecke_matrix(13)
            sage: t.str()
            '[14  0  0  0  0]\n[-4 12  0  8  4]\n[ 0 -6  4 -6  0]\n[ 4  2  0  6 -4]\n[ 0  0  0  0 14]'
        """
        return self.sage_matrix_over_ZZ(sparse=False).str()

    def __dealloc__(self):
        del self.M

    def __getitem__(self, ij):
        """
        Return the (i,j) entry of this matrix.

        Here, ij is a 2-tuple (i,j) and the row and column indices start
        at 1 and not 0.

        EXAMPLES::

            sage: M = CremonaModularSymbols(19, sign=1)
            sage: t = M.hecke_matrix(13); t
            2 x 2 Cremona matrix over Rational Field
            sage: t.sage_matrix_over_ZZ()
            [ 28   0]
            [-12  -8]
            sage: [[t.__getitem__((i,j)) for j in [1,2]] for i in [1,2]]
            [[28, 0], [-12, -8]]
            sage: t.__getitem__((0,0))
            Traceback (most recent call last):
            ...
            IndexError: matrix indices out of range
        """
        cdef long i, j
        if self.M:
            i, j = ij
            if 0<i and i<=self.M[0].nrows() and 0<j and j<=self.M[0].ncols():
                return self.M.sub(i,j)
            raise IndexError("matrix indices out of range")
        raise IndexError("cannot index into an undefined matrix")

    def nrows(self):
        """
        Return the number of rows of this matrix.

        EXAMPLES::

            sage: M = CremonaModularSymbols(19, sign=1)
            sage: t = M.hecke_matrix(13); t
            2 x 2 Cremona matrix over Rational Field
            sage: t.nrows()
            2
        """
        return self.M[0].nrows()

    def ncols(self):
        """
        Return the number of columns of this matrix.

        EXAMPLES::

            sage: M = CremonaModularSymbols(1234, sign=1)
            sage: t = M.hecke_matrix(3); t.ncols()
            156
            sage: M.dimension()
            156
        """
        return self.M[0].ncols()

    # Commented out since it gives very weird
    # results when sign != 0.
##     def rank(self):
##         """
##         Return the rank of this matrix.

##         EXAMPLES::
##
##             sage: M = CremonaModularSymbols(389)
##             sage: t = M.hecke_matrix(2)
##             sage: t.rank()
##             65
##             sage: M = CremonaModularSymbols(389, cuspidal=True)
##             sage: t = M.hecke_matrix(2)
##             sage: t.rank()
##             64

##             sage: M = CremonaModularSymbols(389,sign=1)
##             sage: t = M.hecke_matrix(2)
##             sage: t.rank()   # known bug.
##             16
##         """
##         return rank(self.M[0])

    def add_scalar(self, scalar s):
        """
        Return new matrix obtained by adding `s` to each diagonal entry of ``self``.

        EXAMPLES::

            sage: M = CremonaModularSymbols(23, cuspidal=True, sign=1)
            sage: t = M.hecke_matrix(2); print(t.str())
            [ 0  1]
            [ 1 -1]
            sage: w = t.add_scalar(3); print(w.str())
            [3 1]
            [1 2]
        """
        return new_Matrix(addscalar(self.M[0], s))

    def charpoly(self, var='x'):
        """
        Return the characteristic polynomial of this matrix, viewed as
        as a matrix over the integers.

        ALGORITHM:

        Note that currently, this function converts this matrix into a
        dense matrix over the integers, then calls the charpoly
        algorithm on that, which I think is LinBox's.

        EXAMPLES::

            sage: M = CremonaModularSymbols(33, cuspidal=True, sign=1)
            sage: t = M.hecke_matrix(2)
            sage: t.charpoly()
            x^3 + 3*x^2 - 4
            sage: t.charpoly().factor()
            (x - 1) * (x + 2)^2
        """
        return self.sage_matrix_over_ZZ(sparse=False).charpoly(var)

    def sage_matrix_over_ZZ(self, sparse=True):
        """
        Return corresponding Sage matrix over the integers.

        INPUT:

        - ``sparse`` -- boolean (default: ``True``); whether the return matrix has
          a sparse representation

        EXAMPLES::

            sage: M = CremonaModularSymbols(23, cuspidal=True, sign=1)
            sage: t = M.hecke_matrix(2)
            sage: s = t.sage_matrix_over_ZZ(); s
            [ 0  1]
            [ 1 -1]
            sage: type(s)
            <class 'sage.matrix.matrix_integer_sparse.Matrix_integer_sparse'>
            sage: s = t.sage_matrix_over_ZZ(sparse=False); s
            [ 0  1]
            [ 1 -1]
            sage: type(s)
            <class 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'>
        """
        cdef long n = self.nrows()
        cdef long i, j, k
        cdef scalar* v = <scalar*> self.M.get_entries()   # coercion needed to deal with const

        cdef Matrix_integer_dense Td
        cdef Matrix_integer_sparse Ts

        # Ugly code...
        if sparse:
            Ts = MatrixSpace(ZZ, n, sparse=sparse).zero_matrix().__copy__()
            k = 0
            for i from 0 <= i < n:
                for j from 0 <= j < n:
                    if v[k]:
                        Ts.set_unsafe(i, j, Integer(v[k]))
                    k += 1
            return Ts
        else:
            Td = MatrixSpace(ZZ, n, sparse=sparse).zero_matrix().__copy__()
            k = 0
            for i from 0 <= i < n:
                for j from 0 <= j < n:
                    if v[k]:
                        Td.set_unsafe(i, j, Integer(v[k]))
                    k += 1
            return Td


cdef class MatrixFactory:
    cdef new_matrix(self, mat M):
        return new_Matrix(M)


cdef Matrix new_Matrix(mat M):
    cdef Matrix A = Matrix()
    A.M = new mat(M)
    return A
