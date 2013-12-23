r"""
Base class for matrices, part 2

For design documentation see matrix/docs.py.

AUTHORS:

- William Stein: initial version

- Miguel Marco (2010-06-19): modified eigenvalues and eigenvectors functions to
  allow the option extend=False

- Rob Beezer (2011-02-05): refactored all of the matrix kernel routines

TESTS::

    sage: m = matrix(ZZ['x'], 2, 3, [1..6])
    sage: TestSuite(m).run()
"""

#*****************************************************************************
#       Copyright (C) 2005, 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "sage/ext/stdsage.pxi"
include "sage/ext/python.pxi"
include "sage/ext/interrupt.pxi"

from sage.misc.randstate cimport randstate, current_randstate
from sage.structure.sequence import Sequence
from sage.structure.element import is_Vector
from sage.misc.misc import verbose, get_verbose
from sage.misc.temporary_file import graphics_filename
from sage.rings.number_field.all import is_NumberField
from sage.rings.integer_ring import ZZ, is_IntegerRing
from sage.rings.integer import Integer
from sage.rings.rational_field import QQ, is_RationalField
from sage.rings.real_double import RDF
from sage.rings.complex_double import CDF
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.misc.derivative import multi_derivative

import sage.modules.free_module
import matrix_space
import berlekamp_massey
from sage.modules.free_module_element import is_FreeModuleElement

cdef class Matrix(matrix1.Matrix):
    def _backslash_(self, B):
        r"""
        Used to compute `A \backslash B`, i.e., the backslash solver
        operator.

        EXAMPLES::

            sage: A = matrix(QQ, 3, [1,2,4,2,3,1,0,1,2])
            sage: B = matrix(QQ, 3, 2, [1,7,5, 2,1,3])
            sage: C = A._backslash_(B); C
            [  -1    1]
            [13/5 -3/5]
            [-4/5  9/5]
            sage: A*C == B
            True
            sage: A._backslash_(B) == A \ B
            True
            sage: A._backslash_(B) == A.solve_right(B)
            True
        """
        return self.solve_right(B)

    def subs(self, in_dict=None, **kwds):
        """
        EXAMPLES::

            sage: var('a,b,d,e')
            (a, b, d, e)
            sage: m = matrix([[a,b], [d,e]])
            sage: m.substitute(a=1)
            [1 b]
            [d e]
            sage: m.subs(a=b, b=d)
            [b d]
            [d e]
        """
        v = [a.subs(in_dict, **kwds) for a in self.list()]
        return self.new_matrix(self.nrows(), self.ncols(), v)

    def solve_left(self, B, check=True):
        """
        If self is a matrix `A`, then this function returns a
        vector or matrix `X` such that `X A = B`. If
        `B` is a vector then `X` is a vector and if
        `B` is a matrix, then `X` is a matrix.

        INPUT:


        -  ``B`` - a matrix

        -  ``check`` - bool (default: True) - if False and self
           is nonsquare, may not raise an error message even if there is no
           solution. This is faster but more dangerous.


        EXAMPLES::

            sage: A = matrix(QQ,4,2, [0, -1, 1, 0, -2, 2, 1, 0])
            sage: B = matrix(QQ,2,2, [1, 0, 1, -1])
            sage: X = A.solve_left(B)
            sage: X*A == B
            True

        TESTS::

            sage: A = matrix(QQ,4,2, [0, -1, 1, 0, -2, 2, 1, 0])
            sage: B = vector(QQ,2, [2,1])
            sage: X = A.solve_left(B)
            sage: X*A == B
            True
            sage: X
            (-1, 2, 0, 0)
            sage: A = Matrix(Zmod(128), 2, 3, [5, 29, 33, 64, 0, 7])
            sage: B = vector(Zmod(128), [31,39,56])
            sage: X = A.solve_left(B); X
            (19, 83)
            sage: X * A == B
            True

        """
        if is_Vector(B):
            return self.transpose().solve_right(B, check=check)
        else:
            return self.transpose().solve_right(B.transpose(), check=check).transpose()

    def solve_right(self, B, check=True):
        r"""
        If self is a matrix `A`, then this function returns a
        vector or matrix `X` such that `A X = B`. If
        `B` is a vector then `X` is a vector and if
        `B` is a matrix, then `X` is a matrix.

        .. note::

           In Sage one can also write ``A \backslash  B`` for
           ``A.solve_right(B)``, i.e., Sage implements the "the
           MATLAB/Octave backslash operator".

        INPUT:


        -  ``B`` - a matrix or vector

        -  ``check`` - bool (default: True) - if False and self
           is nonsquare, may not raise an error message even if there is no
           solution. This is faster but more dangerous.


        OUTPUT: a matrix or vector

        .. seealso::

           :meth:`solve_left`

        EXAMPLES::

            sage: A = matrix(QQ, 3, [1,2,3,-1,2,5,2,3,1])
            sage: b = vector(QQ,[1,2,3])
            sage: x = A \ b; x
            (-13/12, 23/12, -7/12)
            sage: A * x
            (1, 2, 3)

        We solve with A nonsquare::

            sage: A = matrix(QQ,2,4, [0, -1, 1, 0, -2, 2, 1, 0]); B = matrix(QQ,2,2, [1, 0, 1, -1])
            sage: X = A.solve_right(B); X
            [-3/2  1/2]
            [  -1    0]
            [   0    0]
            [   0    0]
            sage: A*X == B
            True

        Another nonsingular example::

            sage: A = matrix(QQ,2,3, [1,2,3,2,4,6]); v = vector([-1/2,-1])
            sage: x = A \ v; x
            (-1/2, 0, 0)
            sage: A*x == v
            True

        Same example but over `\ZZ`::

            sage: A = matrix(ZZ,2,3, [1,2,3,2,4,6]); v = vector([-1,-2])
            sage: A \ v
            (-1, 0, 0)

        An example in which there is no solution::

            sage: A = matrix(QQ,2,3, [1,2,3,2,4,6]); v = vector([1,1])
            sage: A \ v
            Traceback (most recent call last):
            ...
            ValueError: matrix equation has no solutions

        A ValueError is raised if the input is invalid::

            sage: A = matrix(QQ,4,2, [0, -1, 1, 0, -2, 2, 1, 0])
            sage: B = matrix(QQ,2,2, [1, 0, 1, -1])
            sage: X = A.solve_right(B)
            Traceback (most recent call last):
            ...
            ValueError: number of rows of self must equal number of rows of B

        We solve with A singular::

            sage: A = matrix(QQ,2,3, [1,2,3,2,4,6]); B = matrix(QQ,2,2, [6, -6, 12, -12])
            sage: X = A.solve_right(B); X
            [ 6 -6]
            [ 0  0]
            [ 0  0]
            sage: A*X == B
            True

        We illustrate left associativity, etc., of the backslash operator.

        ::

            sage: A = matrix(QQ, 2, [1,2,3,4])
            sage: A \ A
            [1 0]
            [0 1]
            sage: A \ A \ A
            [1 2]
            [3 4]
            sage: A.parent()(1) \ A
            [1 2]
            [3 4]
            sage: A \ (A \ A)
            [  -2    1]
            [ 3/2 -1/2]
            sage: X = A \ (A - 2); X
            [ 5 -2]
            [-3  2]
            sage: A * X
            [-1  2]
            [ 3  2]

        Solving over a polynomial ring::

            sage: x = polygen(QQ, 'x')
            sage: A = matrix(2, [x,2*x,-5*x^2+1,3])
            sage: v = vector([3,4*x - 2])
            sage: X = A \ v
            sage: X
            ((-8*x^2 + 4*x + 9)/(10*x^3 + x), (19*x^2 - 2*x - 3)/(10*x^3 + x))
            sage: A * X == v
            True

        Solving some systems over `\ZZ/n\ZZ`::

            sage: A = Matrix(Zmod(6), 3, 2, [1,2,3,4,5,6])
            sage: B = vector(Zmod(6), [1,1,1])
            sage: A.solve_right(B)
            (5, 1)
            sage: B = vector(Zmod(6), [5,1,1])
            sage: A.solve_right(B)
            Traceback (most recent call last):
            ...
            ValueError: matrix equation has no solutions
            sage: A = Matrix(Zmod(128), 2, 3, [23,11,22,4,1,0])
            sage: B = Matrix(Zmod(128), 2, 1, [1,0])
            sage: A.solve_right(B)
            [  1]
            [124]
            [  1]

        Solving a system over the p-adics::

            sage: k = Qp(5,4)
            sage: a = matrix(k, 3, [1,7,3,2,5,4,1,1,2]); a
            [    1 + O(5^4) 2 + 5 + O(5^4)     3 + O(5^4)]
            [    2 + O(5^4)     5 + O(5^5)     4 + O(5^4)]
            [    1 + O(5^4)     1 + O(5^4)     2 + O(5^4)]
            sage: v = vector(k, 3, [1,2,3])
            sage: x = a \ v; x
            (4 + 5 + 5^2 + 3*5^3 + O(5^4), 2 + 5 + 3*5^2 + 5^3 + O(5^4), 1 + 5 + O(5^4))
            sage: a * x == v
            True

        Solving a system of linear equation symbolically using symbolic matrices::

            sage: var('a,b,c,d,x,y')
            (a, b, c, d, x, y)
            sage: A=matrix(SR,2,[a,b,c,d]); A
            [a b]
            [c d]
            sage: result=vector(SR,[3,5]); result
            (3, 5)
            sage: soln=A.solve_right(result)
            sage: soln
            (-b*(3*c/a - 5)/(a*(b*c/a - d)) + 3/a, (3*c/a - 5)/(b*c/a - d))
            sage: (a*x+b*y).subs(x=soln[0],y=soln[1]).simplify_full()
            3
            sage: (c*x+d*y).subs(x=soln[0],y=soln[1]).simplify_full()
            5
            sage: (A*soln).apply_map(lambda x: x.simplify_full())
            (3, 5)
        """

        if is_Vector(B):
            if self.nrows() != B.degree():
                raise ValueError, "number of rows of self must equal degree of B"
        else:
            if self.nrows() != B.nrows():
                raise ValueError, "number of rows of self must equal number of rows of B"

        K = self.base_ring()
        if not K.is_integral_domain():
            from sage.rings.finite_rings.integer_mod_ring import is_IntegerModRing
            if is_IntegerModRing(K):
                from sage.libs.pari.all import pari
                A = pari(self.lift())
                b = pari([c.lift() for c in B]).Col()
                ret = A.matsolvemod(pari(K.cardinality()), b)
                if ret.type() == 't_INT':
                    raise ValueError("matrix equation has no solutions")
                ret = ret.Vec().sage()
                if is_Vector(B):
                    return (K ** self.ncols())(ret)
                return self.matrix_space(self.ncols(), 1)(ret)
            raise TypeError, "base ring must be an integral domain or a ring of integers mod n"
        if not K.is_field():
            K = K.fraction_field()
            self = self.change_ring(K)

        matrix = True
        if is_Vector(B):
            matrix = False
            C = self.matrix_space(self.nrows(), 1)(B.list())
        else:
            C = B

        if not self.is_square():
            X = self._solve_right_general(C, check=check)
            if not matrix:
                # Convert back to a vector
                return (X.base_ring() ** X.nrows())(X.list())
            else:
                return X

        if self.rank() != self.nrows():
            X = self._solve_right_general(C, check=check)
        else:
            X = self._solve_right_nonsingular_square(C, check_rank=False)

        if not matrix:
            # Convert back to a vector
            return X.column(0)
        else:
            return X

    def _solve_right_nonsingular_square(self, B, check_rank=True):
        r"""
        If self is a matrix `A` of full rank, then this function
        returns a matrix `X` such that `A X = B`.

        .. seealso::

           :meth:`solve_right` and :meth:`solve_left`

        INPUT:


        -  ``B`` - a matrix

        -  ``check_rank`` - bool (default: True)


        OUTPUT: matrix

        EXAMPLES::

            sage: A = matrix(QQ,3,[1,2,4,5,3,1,1,2,-1])
            sage: B = matrix(QQ,3,2,[1,5,1,2,1,5])
            sage: A._solve_right_nonsingular_square(B)
            [ -1/7 -11/7]
            [  4/7  23/7]
            [    0     0]
            sage: A._solve_right_nonsingular_square(B, check_rank=False)
            [ -1/7 -11/7]
            [  4/7  23/7]
            [    0     0]
            sage: X = A._solve_right_nonsingular_square(B, check_rank=False)
            sage: A*X == B
            True
        """
        D = self.augment(B).echelon_form()
        return D.matrix_from_columns(range(self.ncols(),D.ncols()))


    def pivot_rows(self):
        """
        Return the pivot row positions for this matrix, which are a topmost
        subset of the rows that span the row space and are linearly
        independent.

        OUTPUT: a tuple of integers

        EXAMPLES::

            sage: A = matrix(QQ,3,3, [0,0,0,1,2,3,2,4,6]); A
            [0 0 0]
            [1 2 3]
            [2 4 6]
            sage: A.pivot_rows()
            (1,)
            sage: A.pivot_rows() # testing cached value
            (1,)
        """
        v = self.fetch('pivot_rows')
        if v is not None:
            return tuple(v)
        v = self.transpose().pivots()
        self.cache('pivot_rows', v)
        return v

    def _solve_right_general(self, B, check=True):
        r"""
        This is used internally by the ``solve_right`` command
        to solve for self\*X = B when self is not square or not of full
        rank. It does some linear algebra, then solves a full-rank square
        system.

        INPUT:


        -  ``B`` - a matrix

        -  ``check`` - bool (default: True); if False, if there
           is no solution this function will not detect that fact.


        OUTPUT: matrix

        EXAMPLES::

            sage: A = matrix(QQ,2,3, [1,2,3,2,4,6]); B = matrix(QQ,2,2, [6, -6, 12, -12])
            sage: A._solve_right_general(B)
            [ 6 -6]
            [ 0  0]
            [ 0  0]
        """
        pivot_cols = self.pivots()
        A = self.matrix_from_columns(pivot_cols)
        pivot_rows = A.pivot_rows()
        A = A.matrix_from_rows(pivot_rows)
        X = A.solve_right(B.matrix_from_rows(pivot_rows), check=False)
        if len(pivot_cols) < self.ncols():
            # Now we have to put in zeros for the non-pivot ROWS, i.e.,
            # make a matrix from X with the ROWS of X interspersed with
            # 0 ROWS.
            Y = X.new_matrix(self.ncols(), X.ncols())
            # Put the columns of X into the matrix Y at the pivot_cols positions
            for i, c in enumerate(pivot_cols):
                Y.set_row(c, X.row(i))
            X = Y
        if check:
            # Have to check that we actually solved the equation.
            if self*X != B:
                raise ValueError, "matrix equation has no solutions"
        return X

    def prod_of_row_sums(self, cols):
        r"""
        Calculate the product of all row sums of a submatrix of `A`
        for a list of selected columns ``cols``.

        EXAMPLES::

            sage: a = matrix(QQ, 2,2, [1,2,3,2]); a
            [1 2]
            [3 2]
            sage: a.prod_of_row_sums([0,1])
            15

        Another example::

            sage: a = matrix(QQ, 2,3, [1,2,3,2,5,6]); a
            [1 2 3]
            [2 5 6]
            sage: a.prod_of_row_sums([1,2])
            55

        AUTHORS:

        - Jaap Spies (2006-02-18)
        """
        cdef Py_ssize_t c, row
        pr = 1
        for row from 0 <= row < self._nrows:
            tmp = []
            for c in cols:
#               if c<0 or c >= self._ncols:
#                   raise IndexError, "matrix column index out of range"
                tmp.append(self.get_unsafe(row, c))
            pr = pr * sum(tmp)
        return pr

    def elementwise_product(self, right):
        r"""
        Returns the elementwise product of two matrices
        of the same size (also known as the Hadamard product).

        INPUT:

        - ``right`` - the right operand of the product.  A matrix
          of the same size as ``self`` such that multiplication
          of elements of the base rings of ``self`` and ``right``
          is defined, once Sage's coercion model is applied.  If
          the matrices have different sizes, or if multiplication
          of individual entries cannot be achieved, a ``TypeError``
          will result.

        OUTPUT:

        A matrix of the same size as ``self`` and ``right``.  The
        entry in location `(i,j)` of the output is the product of
        the two entries in location `(i,j)` of ``self`` and
        ``right`` (in that order).

        The parent of the result is determined by Sage's coercion
        model.  If the base rings are identical, then the result
        is dense or sparse according to this property for
        the left operand.  If the base rings must be adjusted
        for one, or both, matrices then the result will be sparse
        only if both operands are sparse.  No subdivisions are
        present in the result.

        If the type of the result is not to your liking, or
        the ring could be "tighter," adjust the operands with
        :meth:`~sage.matrix.matrix0.Matrix.change_ring`.
        Adjust sparse versus dense inputs with the methods
        :meth:`~sage.matrix.matrix1.Matrix.sparse_matrix` and
        :meth:`~sage.matrix.matrix1.Matrix.dense_matrix`.

        EXAMPLES::

            sage: A = matrix(ZZ, 2, range(6))
            sage: B = matrix(QQ, 2, [5, 1/3, 2/7, 11/2, -3/2, 8])
            sage: C = A.elementwise_product(B)
            sage: C
            [   0  1/3  4/7]
            [33/2   -6   40]
            sage: C.parent()
            Full MatrixSpace of 2 by 3 dense matrices over Rational Field


        Notice the base ring of the results in the next two examples.  ::

            sage: D = matrix(ZZ[x],2,[1+x^2,2,3,4-x])
            sage: E = matrix(QQ,2,[1,2,3,4])
            sage: F = D.elementwise_product(E)
            sage: F
            [  x^2 + 1         4]
            [        9 -4*x + 16]
            sage: F.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Univariate Polynomial Ring in x over Rational Field

        ::

            sage: G = matrix(GF(3),2,[0,1,2,2])
            sage: H = matrix(ZZ,2,[1,2,3,4])
            sage: J = G.elementwise_product(H)
            sage: J
            [0 2]
            [0 2]
            sage: J.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Finite Field of size 3

        Non-commutative rings behave as expected.  These are the usual quaternions. ::

            sage: R.<i,j,k> = QuaternionAlgebra(-1, -1)
            sage: A = matrix(R, 2, [1,i,j,k])
            sage: B = matrix(R, 2, [i,i,i,i])
            sage: A.elementwise_product(B)
            [ i -1]
            [-k  j]
            sage: B.elementwise_product(A)
            [ i -1]
            [ k -j]

        Input that is not a matrix will raise an error.  ::

            sage: A = random_matrix(ZZ,5,10,x=20)
            sage: A.elementwise_product(vector(ZZ, [1,2,3,4]))
            Traceback (most recent call last):
            ...
            TypeError: operand must be a matrix, not an element of Ambient free module of rank 4 over the principal ideal domain Integer Ring

        Matrices of different sizes for operands will raise an error. ::

            sage: A = random_matrix(ZZ,5,10,x=20)
            sage: B = random_matrix(ZZ,10,5,x=40)
            sage: A.elementwise_product(B)
            Traceback (most recent call last):
            ...
            TypeError: incompatible sizes for matrices from: Full MatrixSpace of 5 by 10 dense matrices over Integer Ring and Full MatrixSpace of 10 by 5 dense matrices over Integer Ring

        Some pairs of rings do not have a common parent where
        multiplication makes sense.  This will raise an error. ::

            sage: A = matrix(QQ, 3, range(6))
            sage: B = matrix(GF(3), 3, [2]*6)
            sage: A.elementwise_product(B)
            Traceback (most recent call last):
            ...
            TypeError: no common canonical parent for objects with parents: 'Full MatrixSpace of 3 by 2 dense matrices over Rational Field' and 'Full MatrixSpace of 3 by 2 dense matrices over Finite Field of size 3'

        We illustrate various combinations of sparse and dense matrices.
        Notice how if base rings are unequal, both operands must be sparse
        to get a sparse result. ::

            sage: A = matrix(ZZ, 5, range(30), sparse=False)
            sage: B = matrix(ZZ, 5, range(30), sparse=True)
            sage: C = matrix(QQ, 5, range(30), sparse=True)
            sage: A.elementwise_product(C).is_dense()
            True
            sage: B.elementwise_product(C).is_sparse()
            True
            sage: A.elementwise_product(B).is_dense()
            True
            sage: B.elementwise_product(A).is_dense()
            True

        TESTS:

        Implementation for dense and sparse matrices are
        different, this will provide a trivial test that
        they are working identically. ::

            sage: A = random_matrix(ZZ, 10, x=1000, sparse=False)
            sage: B = random_matrix(ZZ, 10, x=1000, sparse=False)
            sage: C = A.sparse_matrix()
            sage: D = B.sparse_matrix()
            sage: E = A.elementwise_product(B)
            sage: F = C.elementwise_product(D)
            sage: E.is_dense() and F.is_sparse() and (E == F)
            True

        If the ring has zero divisors, the routines for setting
        entries of a sparse matrix should intercept zero results
        and not create an entry. ::

            sage: R = Integers(6)
            sage: A = matrix(R, 2, [3, 2, 0, 0], sparse=True)
            sage: B = matrix(R, 2, [2, 3, 1, 0], sparse=True)
            sage: C = A.elementwise_product(B)
            sage: len(C.nonzero_positions()) == 0
            True

        AUTHOR:

        - Rob Beezer (2009-07-13)
        """
        # Optimized routines for specialized classes would likely be faster
        # See the "pairwise_product" of vectors for some guidance on doing this
        from sage.structure.element import canonical_coercion
        if not PY_TYPE_CHECK(right, Matrix):
            raise TypeError('operand must be a matrix, not an element of %s' % right.parent())
        if (self.nrows() != right.nrows()) or (self.ncols() != right.ncols()):
            raise TypeError('incompatible sizes for matrices from: %s and %s'%(self.parent(), right.parent()))
        if self._parent is not (<Matrix>right)._parent:
            self, right = canonical_coercion(self, right)
        return self._elementwise_product(right)

    def permanent(self):
        r"""
        Calculate and return the permanent of the `m \times n`
        matrix ``self`` using Ryser's algorithm.

        Let `A = (a_{i,j})` be an `m \times n` matrix over
        any commutative ring, with `m \le n`. The permanent of
        `A` is

        .. MATH::

           \mathrm{per}(A)
           = \sum_\pi a_{1,\pi(1)} a_{2,\pi(2)} \cdots a_{m,\pi(m)}

        where the summation extends over all one-to-one functions
        `\pi` from `\{1, \ldots, m\}` to
        `\{1, \ldots, n\}`.

        The product
        `a_{1,\pi(1)} a_{2,\pi(2)} \cdots a_{m,\pi(m)}` is
        called diagonal product. So the permanent of an
        `m \times n` matrix `A` is the sum of all the
        diagonal products of `A`.

        INPUT:

        - ``A`` -- matrix of size `m \times n` with `m \leq n`

        OUTPUT:

        permanent of the matrix `A`

        ALGORITHM:

        Modification of theorem 7.1.1. from Brualdi and Ryser:
        Combinatorial Matrix Theory. Instead of deleting columns from
        `A`, we choose columns from `A` and calculate the
        product of the row sums of the selected submatrix.

        EXAMPLES::

            sage: M = MatrixSpace(ZZ,4,4)
            sage: A = M([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])
            sage: A.permanent()
            24

        ::

            sage: M = MatrixSpace(QQ,3,6)
            sage: A = M([1,1,1,1,0,0,0,1,1,1,1,0,0,0,1,1,1,1])
            sage: A.permanent()
            36

        ::

            sage: M = MatrixSpace(RR,3,6)
            sage: A = M([1.0,1.0,1.0,1.0,0,0,0,1.0,1.0,1.0,1.0,0,0,0,1.0,1.0,1.0,1.0])
            sage: A.permanent()
            36.0000000000000

        See Sloane's sequence OEIS A079908(3) = 36, "The Dancing School
        Problems"

        ::

            sage: oeis(79908)                           # optional -- internet
            A079908: Solution to the Dancing School Problem with 3 girls and n+3 boys: f(3,n).
            sage: _(3)                                  # optional -- internet
            36

        ::

            sage: M = MatrixSpace(ZZ,4,5)
            sage: A = M([1,1,0,1,1,0,1,1,1,1,1,0,1,0,1,1,1,0,1,0])
            sage: A.permanent()
            32

        See Minc: Permanents, Example 2.1, p. 5.

        ::

            sage: M = MatrixSpace(QQ,2,2)
            sage: A = M([1/5,2/7,3/2,4/5])
            sage: A.permanent()
            103/175

        ::

            sage: R.<a> = PolynomialRing(ZZ)
            sage: A = MatrixSpace(R,2)([[a,1], [a,a+1]])
            sage: A.permanent()
            a^2 + 2*a

        ::

            sage: R.<x,y> = PolynomialRing(ZZ,2)
            sage: A = MatrixSpace(R,2)([x, y, x^2, y^2])
            sage: A.permanent()
            x^2*y + x*y^2

        AUTHORS:

        - Jaap Spies (2006-02-16)

        - Jaap Spies (2006-02-21): added definition of permanent
        """
        cdef Py_ssize_t m, n, r
        cdef int sn

        perm = 0
        m = self._nrows
        n = self._ncols
        if not m <= n:
            raise ValueError, "must have m <= n, but m (=%s) and n (=%s)"%(m,n)

        for r from 1 <= r < m+1:
            lst = _choose(n, r)
            tmp = []
            for cols in lst:
                tmp.append(self.prod_of_row_sums(cols))
            s = sum(tmp)
            # sn = (-1)^(m-r)
            if (m - r) % 2 == 0:
                sn = 1
            else:
                sn = -1
            perm = perm + sn * _binomial(n-r, m-r) * s
        return perm


    def permanental_minor(self, Py_ssize_t k):
        r"""
        Return the permanental `k`-minor of an `m \times n` matrix.

        This is the sum of the permanents of all possible `k` by
        `k` submatrices of `A`.

        See Brualdi and Ryser: Combinatorial Matrix Theory, p. 203. Note
        the typo `p_0(A) = 0` in that reference! For applications
        see Theorem 7.2.1 and Theorem 7.2.4.

        Note that the permanental `m`-minor equals
        `\mathrm{per}(A)` if `m = n`.

        For a (0,1)-matrix `A` the permanental `k`-minor
        counts the number of different selections of `k` 1's of
        `A` with no two of the 1's on the same row and no two of the
        1's on the same column.

        INPUT:

        - ``self`` -- matrix of size `m \times n` with `m \leq n`

        OUTPUT:

        The permanental `k`-minor of the matrix ``self``.

        EXAMPLES::

            sage: M = MatrixSpace(ZZ,4,4)
            sage: A = M([1,0,1,0,1,0,1,0,1,0,10,10,1,0,1,1])
            sage: A.permanental_minor(2)
            114

        ::

            sage: M = MatrixSpace(ZZ,3,6)
            sage: A = M([1,1,1,1,0,0,0,1,1,1,1,0,0,0,1,1,1,1])
            sage: A.permanental_minor(0)
            1
            sage: A.permanental_minor(1)
            12
            sage: A.permanental_minor(2)
            40
            sage: A.permanental_minor(3)
            36

        Note that if `k = m = n`, the permanental `k`-minor equals
        `\mathrm{per}(A)`::

            sage: A.permanent()
            36

        ::

            sage: A.permanental_minor(5)
            0

        For `C` the "complement" of `A`::

            sage: M = MatrixSpace(ZZ,3,6)
            sage: C = M([0,0,0,0,1,1,1,0,0,0,0,1,1,1,0,0,0,0])
            sage: m, n = 3, 6
            sage: sum([(-1)^k * C.permanental_minor(k)*factorial(n-k)/factorial(n-m) for k in range(m+1)])
            36

        See Theorem 7.2.1 of Brualdi and Ryser: Combinatorial Matrix
        Theory: per(A)

        AUTHORS:

        - Jaap Spies (2006-02-19)
        """
        m = self._nrows
        n = self._ncols
        if not m <= n:
            raise ValueError("must have m <= n, but m (=%s) and n (=%s)"%(m,n))

        R = self._base_ring
        if k == 0:
            return R.one()
        if k > m:
            return R.zero()

        pm = 0
        for cols in _choose(n,k):
            for rows in _choose(m,k):
                pm = pm + self.matrix_from_rows_and_columns(rows, cols).permanent()
        return pm


    def rook_vector(self, check = False):
        r"""
        Return the rook vector of the matrix ``self``.

        Let `A` be an `m` by `n` (0,1)-matrix with `m \le n`. We identify
        `A` with a chessboard where rooks can be placed on the fields
        `(i, j)` with `a_{i,j} = 1`. The number
        `r_k = p_k(A)` (the permanental `k`-minor) counts the number of
        ways to place `k` rooks on this board so that no rook can attack
        another.

        The rook vector of the matrix `A` is the list consisting of
        `r_0, r_1, \ldots, r_m`.

        The rook polynomial is defined by
        `r(x) = \sum_{k=0}^m r_k x^k`.

        INPUT:

        - ``self`` -- an `m` by `n` (0,1)-matrix with `m \le n`

        - ``check`` -- Boolean (default: ``False``) determining whether
          to check that ``self`` is a (0,1)-matrix.

        OUTPUT:

        The rook vector of the matrix ``self``.

        EXAMPLES::

            sage: M = MatrixSpace(ZZ,3,6)
            sage: A = M([1,1,1,1,0,0,0,1,1,1,1,0,0,0,1,1,1,1])
            sage: A.rook_vector()
            [1, 12, 40, 36]

        ::

            sage: R.<x> = PolynomialRing(ZZ)
            sage: rv = A.rook_vector()
            sage: rook_polynomial = sum([rv[k] * x^k for k in range(len(rv))])
            sage: rook_polynomial
            36*x^3 + 40*x^2 + 12*x + 1

        AUTHORS:

        - Jaap Spies (2006-02-24)
        """
        m = self._nrows
        n = self._ncols
        if not m <= n:
            raise ValueError, "must have m <= n, but m (=%s) and n (=%s)"%(m,n)

        if check:
            # verify that self[i, j] in {0, 1}
            for i in range(m):
                for j in range(n):
                    x = self.get_unsafe(i, j)
                    if not (x == 0 or x == 1):
                        raise ValueError, "must have zero or one, but we have (=%s)"%x

        tmp = []
        for k in range(m+1):
            tmp.append(self.permanental_minor(k))
        return tmp

    def minors(self,k):
        r"""
        Return the list of all `k \times k` minors of self.

        Let `A` be an `m \times n` matrix and `k` an integer with
        `0 \leq k`, `k \leq m` and `k \leq n`.
        A `k \times k` minor of `A` is the determinant of a
        `k \times k` matrix obtained from `A` by deleting `m - k`
        rows and `n - k` columns.
        There are no `k \times k` minors of `A` if `k` is larger
        than either `m` or `n`.

        The returned list is sorted in lexicographical row major ordering,
        e.g., if A is a `3 \times 3` matrix then the minors returned are
        with these rows/columns: [ [0, 1]x[0, 1], [0, 1]x[0, 2], [0, 1]x[1, 2],
        [0, 2]x[0, 1], [0, 2]x[0, 2], [0, 2]x[1, 2], [1, 2]x[0, 1], [1,
        2]x[0, 2], [1, 2]x[1, 2] ].

        INPUT:

        - ``k`` -- integer

        EXAMPLES::

            sage: A = Matrix(ZZ,2,3,[1,2,3,4,5,6]); A
            [1 2 3]
            [4 5 6]
            sage: A.minors(2)
            [-3, -6, -3]
            sage: A.minors(1)
            [1, 2, 3, 4, 5, 6]
            sage: A.minors(0)
            [1]
            sage: A.minors(5)
            []

        ::

            sage: k = GF(37)
            sage: P.<x0,x1,x2> = PolynomialRing(k)
            sage: A = Matrix(P,2,3,[x0*x1, x0, x1, x2, x2 + 16, x2 + 5*x1 ])
            sage: A.minors(2)
            [x0*x1*x2 + 16*x0*x1 - x0*x2, 5*x0*x1^2 + x0*x1*x2 - x1*x2, 5*x0*x1 + x0*x2 - x1*x2 - 16*x1]
        """
        from sage.combinat.combination import Combinations
        all_rows = range(self.nrows())
        all_cols = range(self.ncols())
        m = []
        for rows in Combinations(all_rows,k):
            for cols in Combinations(all_cols,k):
                m.append(self.matrix_from_rows_and_columns(rows,cols).determinant())
        return m

    def determinant(self, algorithm=None):
        r"""
        Returns the determinant of self.

        ALGORITHM:

        For small matrices (n less than 4), this is computed using the naive
        formula. In the specific case of matrices over the integers modulo a
        non-prime, the determinant of a lift is computed over the integers.
        In general, the characteristic polynomial is computed either using
        the Hessenberg form (specified by ``"hessenberg"``) or the generic
        division-free algorithm (specified by ``"df"``).  When the base ring
        is an exact field, the default choice is ``"hessenberg"``, otherwise
        it is ``"df"``.  Note that for matrices over most rings, more
        sophisticated algorithms can be used. (Type ``A.determinant?`` to
        see what is done for a specific matrix ``A``.)

        INPUT:

        - ``algorithm`` - string:
            - ``"df"`` - Generic O(n^4) division-free algorithm
            - ``"hessenberg"`` - Use the Hessenberg form of the matrix

        EXAMPLES::

            sage: A = MatrixSpace(Integers(8),3)([1,7,3, 1,1,1, 3,4,5])
            sage: A.determinant()
            6
            sage: A.determinant() is A.determinant()
            True
            sage: A[0,0] = 10
            sage: A.determinant()
            7

        We compute the determinant of the arbitrary 3x3 matrix::

            sage: R = PolynomialRing(QQ,9,'x')
            sage: A = matrix(R,3,R.gens())
            sage: A
            [x0 x1 x2]
            [x3 x4 x5]
            [x6 x7 x8]
            sage: A.determinant()
            -x2*x4*x6 + x1*x5*x6 + x2*x3*x7 - x0*x5*x7 - x1*x3*x8 + x0*x4*x8

        We create a matrix over `\ZZ[x,y]` and compute its
        determinant.

        ::

            sage: R.<x,y> = PolynomialRing(IntegerRing(),2)
            sage: A = MatrixSpace(R,2)([x, y, x**2, y**2])
            sage: A.determinant()
            -x^2*y + x*y^2

        TESTS::

            sage: A = matrix(5, 5, [next_prime(i^2) for i in range(25)])
            sage: B = MatrixSpace(ZZ['x'], 5, 5)(A)
            sage: A.det() - B.det()
            0

        We verify that :trac:`5569` is resolved (otherwise the following
        would hang for hours)::

            sage: d = random_matrix(GF(next_prime(10^20)),50).det()
            sage: d = random_matrix(Integers(10^50),50).det()

        We verify that trac 7704 is resolved::

            sage: matrix(ZZ, {(0,0):1,(1,1):2,(2,2):3,(3,3):4}).det()
            24
            sage: matrix(QQ, {(0,0):1,(1,1):2,(2,2):3,(3,3):4}).det()
            24

        We verify that trac 10063 is resolved::

            sage: A = GF(2)['x,y,z']
            sage: A.inject_variables()
            Defining x, y, z
            sage: R = A.quotient(x^2 + 1).quotient(y^2 + 1).quotient(z^2 + 1)
            sage: R.inject_variables()
            Defining xbarbarbar, ybarbarbar, zbarbarbar
            sage: M = matrix([[1,1,1,1],[xbarbarbar,ybarbarbar,1,1],[0,1,zbarbarbar,1],[xbarbarbar,zbarbarbar,1,1]])
            sage: M.determinant()
            xbarbarbar*ybarbarbar*zbarbarbar + xbarbarbar*ybarbarbar + xbarbarbar*zbarbarbar + ybarbarbar*zbarbarbar + xbarbarbar + ybarbarbar + zbarbarbar + 1

        Check that the determinant is computed from a cached charpoly
        properly::

            sage: A = matrix(RR, [ [1, 0, 1/2],
            ...                    [0, 1, 0  ],
            ...                    [0, 0, -2 ] ])
            sage: B = copy(A)
            sage: _ = A.charpoly()
            sage: A.determinant() == B.determinant()
            True

        AUTHORS:

          - Unknown: No author specified in the file from 2009-06-25
          - Sebastian Pancratz (2009-06-25): Use the division-free
            algorithm for charpoly
          - Thierry Monteil (2010-10-05): Bugfix for trac ticket #10063,
            so that the determinant is computed even for rings for which
            the is_field method is not implemented.
        """

        from sage.rings.finite_rings.integer_mod_ring import is_IntegerModRing

        if self._nrows != self._ncols:
            raise ValueError, "self must be a square matrix"

        d = self.fetch('det')
        if not d is None:
            return d

        # If charpoly known, then det is easy.
        f = self.fetch('charpoly')
        if f is not None:
            c = f[0]
            if self._nrows % 2 != 0:
                c = -c
            d = self._coerce_element(c)
            self.cache('det', d)
            return d

        cdef Py_ssize_t n
        n = self._ncols
        R = self._base_ring

        # For small matrices, you can't beat the naive formula.
        if n <= 3:
            if n == 0:
                d = R(1)
            elif n == 1:
                d = self.get_unsafe(0,0)
            elif n == 2:
                d = self.get_unsafe(0,0)*self.get_unsafe(1,1) - self.get_unsafe(1,0)*self.get_unsafe(0,1)
            elif n == 3:
                d = self.get_unsafe(0,0) * (self.get_unsafe(1,1)*self.get_unsafe(2,2) - self.get_unsafe(1,2)*self.get_unsafe(2,1))    \
                    - self.get_unsafe(1,0) * (self.get_unsafe(0,1)*self.get_unsafe(2,2) - self.get_unsafe(0,2)*self.get_unsafe(2,1))  \
                    + self.get_unsafe(2,0) * (self.get_unsafe(0,1)*self.get_unsafe(1,2) - self.get_unsafe(0,2)*self.get_unsafe(1,1))
            self.cache('det', d)
            return d

        # Special case for Z/nZ or GF(p):
        if is_IntegerModRing(R) and self.is_dense():
            import sys
            # If the characteristic is prime and smaller than a machine
            # word, use PARI.
            ch = R.characteristic()
            if ch.is_prime() and ch < (2*sys.maxint):
                d = R(self._pari_().matdet())
            else:
                # Lift to ZZ and compute there.
                d = R(self.apply_map(lambda x : x.centerlift()).det())
            self.cache('det', d)
            return d

        # N.B.  The following comment should be obsolete now that the generic
        # code to compute the determinant has been replaced by generic
        # division-free code to compute the characteristic polynomial and read
        # off the determinant from this.
        #
        # If R is an exact integral domain, we could get the det by computing
        # charpoly.  The generic fraction field implementation is so slow that
        # the naive algorithm is much faster in practice despite
        # asymptotics.
        # TODO: Find a reasonable cutoff point.  (This is field specific, but
        # seems to be quite large for Q[x].)
        if algorithm is None:
            try:
                R_is_field_attempt = R.is_field()
            except NotImplementedError:
                R_is_field_attempt = False
        if (algorithm is None and R_is_field_attempt and R.is_exact()) or (algorithm == "hessenberg"):
            try:
                c = self.charpoly('x', algorithm="hessenberg")[0]
            except ValueError:
                # Hessenberg algorithm not supported, so we use whatever the default algorithm is.
                c = self.charpoly('x')[0]
            if self._nrows % 2:
                c = -c
            d = self._coerce_element(c)
            self.cache('det', d)
            return d

        # Generic division-free algorithm to compute the characteristic
        # polynomial.
        #
        # N.B.   The case of symbolic rings is quite specific.  It would be
        # nice to avoid hardcoding a reserved variable name as below, as this
        # is then assumed to not be a variable in the symbolic ring.  But this
        # resulted in further exceptions/ errors.
        from sage.symbolic.ring import is_SymbolicExpressionRing

        var = R('A0123456789') if is_SymbolicExpressionRing(R) else 'x'
        try:
            c = self.charpoly(var, algorithm="df")[0]
        except ValueError:
            # Division free algorithm not supported, so we use whatever the default algorithm is.
            c = self.charpoly(var)[0]

        if self._nrows % 2:
            c = -c
        d = self._coerce_element(c)
        self.cache('det', d)
        return d

    cdef _det_by_minors(self, Py_ssize_t level):
        """
        Compute the determinant of the upper-left level x level submatrix
        of self. Does not handle degenerate cases, level MUST be >= 2
        """
        cdef Py_ssize_t n, i
        if level == 2:
            return self.get_unsafe(0,0) * self.get_unsafe(1,1) - self.get_unsafe(0,1) * self.get_unsafe(1,0)
        else:
            level -= 1
            d = self.get_unsafe(level,level) * self._det_by_minors(level)
            # on each iteration, row i will be missing in the first (level) rows
            # swapping is much faster than taking submatrices
            for i from level > i >= 0:
                self.swap_rows(level, i)
                if (level - i) % 2:
                    d -= self.get_unsafe(level,level) * self._det_by_minors(level)
                else:
                     d += self.get_unsafe(level,level) * self._det_by_minors(level)
            # undo all our permutations to get us back to where we started
            for i from 0 <= i < level:
                self.swap_rows(level, i)
            return d

    def pfaffian(self, algorithm=None, check=True):
        r"""
        Return the Pfaffian of ``self``, assuming that ``self`` is an
        alternating matrix.

        INPUT:

        - ``algorithm`` -- string, the algorithm to use; currently the
          following algorithms have been implemented:

          * ``'definition'`` - using the definition given by perfect
            matchings

        - ``check`` (default: ``True``) -- Boolean determining whether to
          check ``self`` for alternatingness and squareness. This has to
          be set to ``False`` if ``self`` is defined over a non-discrete
          ring.

        The Pfaffian of an alternating matrix is defined as follows:

        Let `A` be an alternating `k \times k` matrix over a commutative
        ring. (Here, "alternating" means that `A^T = -A` and that the
        diagonal entries of `A` are zero.)
        If `k` is odd, then the Pfaffian of the matrix `A` is
        defined to be `0`. Let us now define it when `k` is even. In this
        case, set `n = k/2` (this is an integer). For every `i` and `j`,
        we denote the `(i, j)`-th entry of `A` by `a_{i, j}`. Let `M`
        denote the set of all perfect matchings of the set
        `\{ 1, 2, \ldots, 2n \}` (see
        :class:`sage.combinat.perfect_matching.PerfectMatchings` ).
        For every matching `m \in M`, define the sign `\mathrm{sign}(m)`
        of `m` by writing `m` as `\{ \{ i_1, j_1 \}, \{ i_2, j_2 \},
        \ldots, \{ i_n, j_n \} \}` with `i_k < j_k` for all `k`, and
        setting `\mathrm{sign}(m)` to be the sign of the permutation
        `( i_1, j_1, i_2, j_2, \ldots, i_n, j_n )` (written here in
        one-line notation). For every matching `m \in M`, define the
        weight `w(m)` of `m` by writing `m` as
        `\{ \{ i_1, j_1 \}, \{ i_2, j_2 \}, \ldots, \{ i_n, j_n \} \}`
        with `i_k < j_k` for all `k`, and setting
        `w(m) = a_{i_1, j_1} a_{i_2, j_2} \cdots a_{i_n, j_n}`. Now, the
        Pfaffian of the matrix `A` is defined to be the sum

        .. MATH::

            \sum_{m \in M} \mathrm{sign}(m) w(m).

        The Pfaffian of `A` is commonly denoted by `\mathrm{Pf}(A)`. It
        is well-known that `(\mathrm{Pf}(A))^2 = \det A` for every
        alternating matrix `A`, and that
        `\mathrm{Pf} (U^T A U) = \det U \cdot \mathrm{Pf}(A)` for any
        `n \times n` matrix `U` and any alternating `n \times n`
        matrix `A`.

        See [Kn95]_, [DW95]_ and [Rote2001]_, just to name three
        sources, for further properties of Pfaffians.

        ALGORITHM:

        The current implementation uses the definition given above.
        It checks alternatingness of the matrix ``self`` only if
        ``check`` is ``True`` (this is important because even if ``self``
        is alternating, a non-discrete base ring might prevent Sage
        from being able to check this).

        REFERENCES:

        .. [Kn95] Donald E. Knuth, *Overlapping Pfaffians*,
           :arxiv:`math/9503234v1`.

        .. [Rote2001] Gunter Rote,
           *Division-Free Algorithms for the Determinant and the
           Pfaffian: Algebraic and Combinatorial Approaches*,
           H. Alt (Ed.): Computational Discrete Mathematics, LNCS
           2122, pp. 119–135, 2001.
           http://page.mi.fu-berlin.de/rote/Papers/pdf/Division-free+algorithms.pdf

        .. [DW95] Andreas W.M. Dress, Walter Wenzel,
           *A Simple Proof of an Identity Concerning Pfaffians of
           Skew Symmetric Matrices*,
           Advances in Mathematics, volume 112, Issue 1, April
           1995, pp. 120-134.
           http://www.sciencedirect.com/science/article/pii/S0001870885710298

        .. TODO::

            Implement faster algorithms, including a division-free one.
            Does [Rote2001]_, section 3.3 give one?

            Check the implementation of the matchings used here for
            performance?

        EXAMPLES:

        A `3 \times 3` alternating matrix has Pfaffian 0 independently
        of its entries::

            sage: MSp = MatrixSpace(Integers(27), 3)
            sage: A = MSp([0, 2, -3,  -2, 0, 8,  3, -8, 0])
            sage: A.pfaffian()
            0
            sage: parent(A.pfaffian())
            Ring of integers modulo 27

        The Pfaffian of a `2 \times 2` alternating matrix is just its
        northeast entry::

            sage: MSp = MatrixSpace(QQ, 2)
            sage: A = MSp([0, 4,  -4, 0])
            sage: A.pfaffian()
            4
            sage: parent(A.pfaffian())
            Rational Field

        The Pfaffian of a `0 \times 0` alternating matrix is `1`::

            sage: MSp = MatrixSpace(ZZ, 0)
            sage: A = MSp([])
            sage: A.pfaffian()
            1
            sage: parent(A.pfaffian())
            Integer Ring

        Let us compute the Pfaffian of a generic `4 \times 4`
        alternating matrix::

            sage: R = PolynomialRing(QQ, 'x12,x13,x14,x23,x24,x34')
            sage: x12, x13, x14, x23, x24, x34 = R.gens()
            sage: A = matrix(R, [[   0,  x12,  x13,  x14],
            ....:                [-x12,    0,  x23,  x24],
            ....:                [-x13, -x23,    0,  x34],
            ....:                [-x14, -x24, -x34,    0]])
            sage: A.pfaffian()
            x14*x23 - x13*x24 + x12*x34
            sage: parent(A.pfaffian())
            Multivariate Polynomial Ring in x12, x13, x14, x23, x24, x34 over Rational Field

        The Pfaffian of an alternating matrix squares to its
        determinant::

            sage: A = [[0] * 6 for i in range(6)]
            sage: for i in range(6):
            ....:     for j in range(i):
            ....:         u = floor(random() * 10)
            ....:         A[i][j] = u
            ....:         A[j][i] = -u
            ....:     A[i][i] = 0
            sage: AA = Matrix(ZZ, A)
            sage: AA.pfaffian() ** 2 == AA.det()
            True

        AUTHORS:

        - Darij Grinberg (1 Oct 2013): first (slow)
          implementation.
        """
        k = self._nrows

        if check:
            if k != self._ncols:
                raise ValueError("self must be a square matrix")
            if not self.is_alternating():
                raise ValueError("self must be alternating, which includes the diagonal entries being 0")

        R = self.base_ring()

        if k % 2 == 1:
            return R.zero()

        if k == 0:
            return R.one()

        n = k // 2

        res = R.zero()

        from sage.combinat.perfect_matching import PerfectMatchings
        from sage.misc.flatten import flatten
        from sage.misc.misc_c import prod
        from sage.combinat.permutation import Permutation
        for m in PerfectMatchings(k):
            # We only need each edge of the matching to be sorted,
            # not the totality of edges.
            edges = [sorted(edge) for edge in list(m)]
            sgn = Permutation(flatten(edges)).sign()
            # Subtract 1 from everything for indexing purposes:
            edges2 = [[i-1 for i in edge] for edge in edges]
            # Product without base case since k == 0 case has
            # already been dealt with:
            res += sgn * prod([self.get_unsafe(edge[0], edge[1]) for edge in edges2])

        return res

    # shortcuts
    def det(self, *args, **kwds):
        """
        Synonym for self.determinant(...).

        EXAMPLES::

            sage: A = MatrixSpace(Integers(8),3)([1,7,3, 1,1,1, 3,4,5])
            sage: A.det()
            6
        """
        return self.determinant(*args, **kwds)

    def __abs__(self):
        """
        Synonym for self.determinant(...).

        EXAMPLES::

            sage: a = matrix(QQ, 2,2, [1,2,3,4]); a
            [1 2]
            [3 4]
            sage: abs(a)
            -2
        """
        return self.determinant()

    def characteristic_polynomial(self, *args, **kwds):
        """
        Synonym for self.charpoly(...).

        EXAMPLES::

            sage: a = matrix(QQ, 2,2, [1,2,3,4]); a
            [1 2]
            [3 4]
            sage: a.characteristic_polynomial('T')
            T^2 - 5*T - 2
        """
        return self.charpoly(*args, **kwds)

    def minimal_polynomial(self, var='x', **kwds):
        r"""
        This is a synonym for ``self.minpoly``

        EXAMPLES::

            sage: a = matrix(QQ, 4, range(16))
            sage: a.minimal_polynomial('z')
            z^3 - 30*z^2 - 80*z
            sage: a.minpoly()
            x^3 - 30*x^2 - 80*x
        """
        return self.minpoly(var, **kwds)

    def minpoly(self, var='x', **kwds):
        r"""
        Return the minimal polynomial of self.

        This uses a simplistic - and potentially very very slow - algorithm
        that involves computing kernels to determine the powers of the
        factors of the charpoly that divide the minpoly.

        EXAMPLES::

            sage: A = matrix(GF(9,'c'), 4, [1, 1, 0,0, 0,1,0,0, 0,0,5,0, 0,0,0,5])
            sage: factor(A.minpoly())
            (x + 1) * (x + 2)^2
            sage: A.minpoly()(A) == 0
            True
            sage: factor(A.charpoly())
            (x + 1)^2 * (x + 2)^2

        The default variable name is `x`, but you can specify
        another name::

            sage: factor(A.minpoly('y'))
            (y + 1) * (y + 2)^2

        We can take the minimal polynomial of symbolic matrices::

            sage: t = var('t')
            sage: m = matrix(2,[1,2,4,t])
            sage: m.minimal_polynomial()
            x^2 + (-t - 1)*x + t - 8
        """
        f = self.fetch('minpoly')
        if not f is None:
            return f.change_variable_name(var)
        f = self.charpoly(var=var, **kwds)
        if f.is_squarefree():  # is_squarefree for polys much faster than factor.
            # Then f must be the minpoly
            self.cache('minpoly', f)
            return f

        # Now we have to work harder.  We find the power of each
        # irreducible factor that divides the minpoly.
        mp = f.radical()
        for h, e in f.factor():
            if e > 1:
                # Find the power of B so that the dimension
                # of the kernel equals e*deg(h)
                B = h(self)   # this is the killer.
                C = B
                n = 1
                while C.kernel().dimension() < e*h.degree():
                    if n == e - 1:
                        n += 1
                        break
                    C *= B
                    n += 1
                mp *= h**(n-1)
        self.cache('minpoly', mp)
        return mp

    def charpoly(self, var = 'x', algorithm = None):
        r"""
        Returns the characteristic polynomial of self, as a polynomial over
        the base ring.

        ALGORITHM:

        In the generic case of matrices over a ring (commutative and with
        unity), there is a division-free algorithm, which can be accessed
        using ``"df"``, with complexity `O(n^4)`.  Alternatively, by
        specifying ``"hessenberg"``, this method computes the Hessenberg
        form of the matrix and then reads off the characteristic polynomial.
        Moreover, for matrices over number fields, this method can use
        PARI's charpoly implementation instead.

        The method's logic is as follows:  If no algorithm is specified,
        first check if the base ring is a number field (and then use PARI),
        otherwise check if the base ring is the ring of integers modulo n (in
        which case compute the characteristic polynomial of a lift of the
        matrix to the integers, and then coerce back to the base), next check
        if the base ring is an exact field (and then use the Hessenberg form),
        or otherwise, use the generic division-free algorithm.
        If an algorithm is specified explicitly, if
        ``algorithm == "hessenberg"``, use the Hessenberg form, or otherwise
        use the generic division-free algorithm.

        The result is cached.

        INPUT:

        - ``var`` - a variable name (default: 'x')
        - ``algorithm`` - string:
            - ``"df"`` - Generic `O(n^4)` division-free algorithm
            - ``"hessenberg"`` - Use the Hessenberg form of the matrix

        EXAMPLES:

        First a matrix over `\ZZ`::

            sage: A = MatrixSpace(ZZ,2)( [1,2,  3,4] )
            sage: f = A.charpoly('x')
            sage: f
            x^2 - 5*x - 2
            sage: f.parent()
            Univariate Polynomial Ring in x over Integer Ring
            sage: f(A)
            [0 0]
            [0 0]

        An example over `\QQ`::

            sage: A = MatrixSpace(QQ,3)(range(9))
            sage: A.charpoly('x')
            x^3 - 12*x^2 - 18*x
            sage: A.trace()
            12
            sage: A.determinant()
            0

        We compute the characteristic polynomial of a matrix over the
        polynomial ring `\ZZ[a]`::

            sage: R.<a> = PolynomialRing(ZZ)
            sage: M = MatrixSpace(R,2)([a,1,  a,a+1]); M
            [    a     1]
            [    a a + 1]
            sage: f = M.charpoly('x'); f
            x^2 + (-2*a - 1)*x + a^2
            sage: f.parent()
            Univariate Polynomial Ring in x over Univariate Polynomial Ring in a over Integer Ring
            sage: M.trace()
            2*a + 1
            sage: M.determinant()
            a^2

        We compute the characteristic polynomial of a matrix over the
        multi-variate polynomial ring `\ZZ[x,y]`::

            sage: R.<x,y> = PolynomialRing(ZZ,2)
            sage: A = MatrixSpace(R,2)([x, y, x^2, y^2])
            sage: f = A.charpoly('x'); f
            x^2 + (-y^2 - x)*x - x^2*y + x*y^2

        It's a little difficult to distinguish the variables. To fix this,
        we temporarily view the indeterminate as `Z`::

            sage: with localvars(f.parent(), 'Z'): print f
            Z^2 + (-y^2 - x)*Z - x^2*y + x*y^2

        We could also compute f in terms of Z from the start::

            sage: A.charpoly('Z')
            Z^2 + (-y^2 - x)*Z - x^2*y + x*y^2

        Here is an example over a number field::

            sage: x = QQ['x'].gen()
            sage: K.<a> = NumberField(x^2 - 2)
            sage: m = matrix(K, [[a-1, 2], [a, a+1]])
            sage: m.charpoly('Z')
            Z^2 - 2*a*Z - 2*a + 1
            sage: m.charpoly('a')(m) == 0
            True

        Here is an example over a general commutative ring, that is to say,
        as of version 4.0.2, SAGE does not even positively determine that
        ``S`` in the following example is an integral domain.  But the
        computation of the characteristic polynomial succeeds as follows::

            sage: R.<a,b> = QQ[]
            sage: S.<x,y> = R.quo((b^3))
            sage: A = matrix(S, [[x*y^2,2*x],[2,x^10*y]])
            sage: A
            [ x*y^2    2*x]
            [     2 x^10*y]
            sage: A.charpoly('T')
            T^2 + (-x^10*y - x*y^2)*T - 4*x

        TESTS::

            sage: P.<a,b,c> = PolynomialRing(Rationals())
            sage: u = MatrixSpace(P,3)([[0,0,a],[1,0,b],[0,1,c]])
            sage: Q.<x> = PolynomialRing(P)
            sage: u.charpoly('x')
            x^3 - c*x^2 - b*x - a

        A test case from :trac:`6442`. Prior to :trac:`12292`, the
        call to ``A.det()`` would attempt to use the cached charpoly,
        and crash if an empty dictionary was cached. We don't cache
        dictionaries anymore, but this test should still pass::

            sage: z = Zp(p=5)
            sage: A = matrix(z, [ [3 + O(5^1), 4 + O(5^1), 4 + O(5^1)],
            ...                   [2*5^2 + O(5^3), 2 + O(5^1), 1 + O(5^1)],
            ...                   [5 + O(5^2), 1 + O(5^1), 1 + O(5^1)] ])
            sage: A.charpoly(algorithm='hessenberg')
            Traceback (most recent call last):
            ...
            ValueError: element has negative valuation.
            sage: A.det()
            3 + O(5)

        The cached polynomial should be independent of the ``var``
        argument (:trac:`12292`). We check (indirectly) that the
        second call uses the cached value by noting that its result is
        not cached::

            sage: M = MatrixSpace(RR, 2)
            sage: A = M(range(0, 2^2))
            sage: type(A)
            <type 'sage.matrix.matrix_generic_dense.Matrix_generic_dense'>
            sage: A.charpoly('x')
            x^2 - 3.00000000000000*x - 2.00000000000000
            sage: A.charpoly('y')
            y^2 - 3.00000000000000*y - 2.00000000000000
            sage: A._cache['charpoly']
            x^2 - 3.00000000000000*x - 2.00000000000000

        AUTHORS:

        - Unknown: No author specified in the file from 2009-06-25
        - Sebastian Pancratz (2009-06-25): Include the division-free algorithm
        """

        f = self.fetch('charpoly')
        if f is not None:
            return f.change_variable_name(var)

        from sage.rings.finite_rings.integer_mod_ring import is_IntegerModRing

        if algorithm is None:
            R = self._base_ring
            if is_NumberField(R):
                f = self._charpoly_over_number_field(var)
            elif is_IntegerModRing(R):
                f = self.lift().charpoly(var).change_ring(R)
            elif R.is_field(proof = False) and R.is_exact():
                f = self._charpoly_hessenberg(var)
            else:
                f = self._charpoly_df(var)
        else:
            if algorithm == "hessenberg":
                f = self._charpoly_hessenberg(var)
            else:
                f = self._charpoly_df(var)

        # Cache the result, and return it.
        self.cache('charpoly', f)
        return f

    def _charpoly_df(self, var = 'x'):
        r"""
        Computes the characteristic polynomial of ``self`` without divisions.

        INPUT:

        - ``var`` - a variable name (default: ``'x'``)

        OUTPUT:

        - polynomial - the characteristic polynomial of ``self``

        EXAMPLES:

        Here is one easy example over the integers to illustrate this::

            sage: A = matrix(ZZ, [[1,24],[3,5]])
            sage: A
            [ 1 24]
            [ 3  5]
            sage: A._charpoly_df()
            x^2 - 6*x - 67

        The second example is a matrix over a univariate polynomial ring over the
        rationals::

            sage: R.<t> = QQ[]
            sage: A = matrix(R, [[7*t^2 - t - 9, -1/4*t - 1, -17*t^2 - t + 1], \
                                 [-t^2 + 1/4*t, t^2 + 5/7*t + 3, 1/5*t^2 +     \
                                  1662],                                       \
                                 [-2*t - 3, 2*t^2 + 6*t - 1/2, -1/6*t^2]])
            sage: A
            [    7*t^2 - t - 9        -1/4*t - 1   -17*t^2 - t + 1]
            [     -t^2 + 1/4*t   t^2 + 5/7*t + 3    1/5*t^2 + 1662]
            [         -2*t - 3 2*t^2 + 6*t - 1/2          -1/6*t^2]
            sage: A._charpoly_df()
            x^3 + (-47/6*t^2 + 2/7*t + 6)*x^2 + (79/15*t^4 - 13189/420*t^3 -
            1884709/560*t^2 - 279501/28*t + 807)*x - 901/30*t^6 - 423/8*t^5 +
            11218517/480*t^4 + 2797765/42*t^3 - 12987971/280*t^2 - 5235245/56*t + 2484

        Thirdly, an example over a ring which is not an integral domain::

            sage: A = matrix(ZZ.quo(12), 3, [5,8,0,10,2,1,8,7,9])
            sage: A
            [ 5  8  0]
            [10  2  1]
            [ 8  7  9]
            sage: A._charpoly_df()
            x^3 + 8*x^2 + 10*x + 1

        TESTS::

            sage: A = matrix(ZZ, 0, 0)
            sage: A
            []
            sage: A._charpoly_df()
            1

            sage: A = matrix(ZZ, 1, 1, [[23]])
            sage: A._charpoly_df()
            x - 23

        NOTES:

        The key feature of this implementation is that it is division-free.
        This means that it can be used as a generic implementation for any
        ring (commutative and with multiplicative identity).  The algorithm
        is described in full detail as Algorithm 3.1 in [Se02].

        Note that there is a missing minus sign in front of the last term in
        the penultimate line of Algorithm 3.1.

        REFERENCES:

        - [Se02] T. R. Seifullin, "Computation of determinants, adjoint
          matrices, and characteristic polynomials without division"

        AUTHORS:

        - Sebastian Pancratz (2009-06-12)
        """

        # Validate assertions
        #
        if not self.is_square():
            raise ValueError("self must be a square matrix")

        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

        # Extract parameters
        #
        cdef Matrix M  = <Matrix> self
        n  = M._ncols
        R  = M._base_ring
        S  = PolynomialRing(R, var)

        # Corner cases
        # N.B.  We already tested for M to be square, hence we do not need to
        # test for 0 x n or m x 0 matrices.
        #
        if n == 0:
            return S(1)

        # In the notation of Algorithm 3.1,
        #
        #   F  Matrix over R such that F[p, t] is the coefficient of X^{t-p}
        #      in the polynomial $F_t$;
        #   a  List of lists such that a[p, t] is a vector of length t;
        #   A  Matrix over R.
        #
        # But by looking at the algorithm, we see that in F, a and A we can
        # drop the second index t, reducing storage requirements.
        #
        # N.B.  The documentation is still 1-based, although the code, after
        # having been ported from Magma to SAGE, is 0-based.
        #
        from sage.matrix.constructor import matrix

        F = [R(0) for i in xrange(n)]
        cdef Matrix a = <Matrix> matrix(R, n-1, n)
        A = [R(0) for i in xrange(n)]

        F[0] = - M.get_unsafe(0, 0)
        for t in xrange(1,n):

            # Set a(1, t) to be M(<=t, t)
            #
            for i in xrange(t+1):
                a.set_unsafe(0, i, M.get_unsafe(i, t))

            # Set A[1, t] to be the (t)th entry in a[1, t]
            #
            A[0] = M.get_unsafe(t, t)

            for p in xrange(1, t):

                # Set a(p, t) to the product of M[<=t, <=t] * a(p-1, t)
                #
                for i in xrange(t+1):
                    s = R(0)
                    for j in xrange(t+1):
                        s = s + M.get_unsafe(i, j) * a.get_unsafe(p-1, j)
                    a.set_unsafe(p, i, s)

                # Set A[p, t] to be the (t)th entry in a[p, t]
                #
                A[p] = a.get_unsafe(p, t)

            # Set A[t, t] to be M[t, <=t] * a(p-1, t)
            #
            s = R(0)
            for j in xrange(t+1):
                s = s + M.get_unsafe(t, j) * a.get_unsafe(t-1, j)
            A[t] = s

            for p in xrange(t+1):
                s = F[p]
                for k in xrange(p):
                    s = s - A[k] * F[p-k-1]
                F[p] = s - A[p]

        X = S.gen(0)
        f = X ** n
        for p in xrange(n):
            f = f + F[p] * X ** (n-1-p)

        return f

    def _charpoly_over_number_field(self, var='x'):
        r"""
        Use PARI to compute the characteristic polynomial of self as a
        polynomial over the base ring.

        EXAMPLES::

            sage: x = QQ['x'].gen()
            sage: K.<a> = NumberField(x^2 - 2)
            sage: m = matrix(K, [[a-1, 2], [a, a+1]])
            sage: m._charpoly_over_number_field('Z')
            Z^2 - 2*a*Z - 2*a + 1
            sage: m._charpoly_over_number_field('a')(m) == 0
            True
            sage: m = matrix(K, [[0, a, 0], [-a, 0, 0], [0, 0, 0]])
            sage: m._charpoly_over_number_field('Z')
            Z^3 + 2*Z

        The remaining tests are indirect::

            sage: L.<b> = K.extension(x^3 - a)
            sage: m = matrix(L, [[b+a, 1], [a, b^2-2]])
            sage: m.charpoly('Z')
            Z^2 + (-b^2 - b - a + 2)*Z + a*b^2 - 2*b - 2*a
            sage: m.charpoly('a')
            a^2 + (-b^2 - b - a + 2)*a + a*b^2 - 2*b - 2*a
            sage: m.charpoly('a')(m) == 0
            True

        ::

            sage: M.<c> = L.extension(x^2 - a*x + b)
            sage: m = matrix(M, [[a+b+c, 0, b], [0, c, 1], [a-1, b^2+1, 2]])
            sage: f = m.charpoly('Z'); f
            Z^3 + (-2*c - b - a - 2)*Z^2 + ((b + 2*a + 4)*c - b^2 + (-a + 2)*b + 2*a - 1)*Z + (b^2 + (a - 3)*b - 4*a + 1)*c + a*b^2 + 3*b + 2*a
            sage: f(m) == 0
            True
            sage: f.base_ring() is M
            True
        """
        K = self.base_ring()
        if not is_NumberField(K):
            raise ValueError, "_charpoly_over_number_field called with base ring (%s) not a number field" % K

        paripoly = self._pari_().charpoly()
        return K[var](paripoly)

    def fcp(self, var='x'):
        """
        Return the factorization of the characteristic polynomial of self.

        INPUT:

        -  ``var`` - (default: 'x') name of variable of charpoly

        EXAMPLES::

            sage: M = MatrixSpace(QQ,3,3)
            sage: A = M([1,9,-7,4/5,4,3,6,4,3])
            sage: A.fcp()
            x^3 - 8*x^2 + 209/5*x - 286
            sage: A = M([3, 0, -2, 0, -2, 0, 0, 0, 0])
            sage: A.fcp('T')
            (T - 3) * T * (T + 2)
        """
        return self.charpoly(var).factor()

##     def minimal_polynomial(self, var, algorithm=''):
##         """
##         Synonym for self.charpoly(...).

##         EXAMPLES:
##             sage: ???
##         """
##         return self.minpoly(*args, **kwds)

##     def minpoly(self, *args, **kwds):
##         """
##         EXAMPLES:
##             sage: ???
##         """
##         raise NotImplementedError

    def denominator(self):
        r"""
        Return the least common multiple of the denominators of the
        elements of self.

        If there is no denominator function for the base field, or no LCM
        function for the denominators, raise a TypeError.

        EXAMPLES::

            sage: A = MatrixSpace(QQ,2)(['1/2', '1/3', '1/5', '1/7'])
            sage: A.denominator()
            210

        A trivial example::

            sage: A = matrix(QQ, 0,2)
            sage: A.denominator()
            1

        Denominators are not defined for real numbers::

            sage: A = MatrixSpace(RealField(),2)([1,2,3,4])
            sage: A.denominator()
            Traceback (most recent call last):
            ...
            TypeError: denominator not defined for elements of the base ring

        We can even compute the denominator of matrix over the fraction
        field of `\ZZ[x]`.

        ::

            sage: K.<x> = Frac(ZZ['x'])
            sage: A = MatrixSpace(K,2)([1/x, 2/(x+1), 1, 5/(x^3)])
            sage: A.denominator()
            x^4 + x^3

        Here's an example involving a cyclotomic field::

            sage: K.<z> = CyclotomicField(3)
            sage: M = MatrixSpace(K,3,sparse=True)
            sage: A = M([(1+z)/3,(2+z)/3,z/3,1,1+z,-2,1,5,-1+z])
            sage: print A
            [1/3*z + 1/3 1/3*z + 2/3       1/3*z]
            [          1       z + 1          -2]
            [          1           5       z - 1]
            sage: print A.denominator()
            3
        """
        if self.nrows() == 0 or self.ncols() == 0:
            return ZZ(1)
        R = self.base_ring()
        x = self.list()
        try:
            d = x[0].denominator()
        except AttributeError:
            raise TypeError, "denominator not defined for elements of the base ring"
        try:
            for y in x:
                d = d.lcm(y.denominator())
        except AttributeError:
            raise TypeError, "lcm function not defined for elements of the base ring"
        return d

    def diagonal(self):
      r"""
      Return the diagonal entries of ``self``.

      OUTPUT:

      A list containing the entries of the matrix that
      have equal row and column indices, in order of the
      indices.  Behavior is not limited to square matrices.

      EXAMPLES::

          sage: A = matrix([[2,5],[3,7]]); A
          [2 5]
          [3 7]
          sage: A.diagonal()
          [2, 7]

      Two rectangular matrices.  ::

          sage: B = matrix(3, 7, range(21)); B
          [ 0  1  2  3  4  5  6]
          [ 7  8  9 10 11 12 13]
          [14 15 16 17 18 19 20]
          sage: B.diagonal()
          [0, 8, 16]

          sage: C = matrix(3, 2, range(6)); C
          [0 1]
          [2 3]
          [4 5]
          sage: C.diagonal()
          [0, 3]

      Empty matrices behave properly. ::

          sage: E = matrix(0, 5, []); E
          []
          sage: E.diagonal()
          []
      """
      n = min(self.nrows(), self.ncols())
      return [self[i,i] for i in range(n)]

    def trace(self):
        """
        Return the trace of self, which is the sum of the diagonal entries
        of self.

        INPUT:


        -  ``self`` - a square matrix


        OUTPUT: element of the base ring of self

        EXAMPLES::

            sage: a = matrix(3,range(9)); a
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: a.trace()
            12
            sage: a = matrix({(1,1):10, (2,1):-3, (2,2):4/3}); a
            [  0   0   0]
            [  0  10   0]
            [  0  -3 4/3]
            sage: a.trace()
            34/3
        """
        if self._nrows != self._ncols:
            raise ValueError, "self must be a square matrix"
        R = self._base_ring
        cdef Py_ssize_t i
        cdef object s
        s = R(0)
        for i from 0 <= i < self._nrows:
            s = s + self.get_unsafe(i,i)
        return s

    def trace_of_product(self, Matrix other):
        """
        Returns the trace of self * other without computing the entire product.

        EXAMPLES::

            sage: M = random_matrix(ZZ, 10, 20)
            sage: N = random_matrix(ZZ, 20, 10)
            sage: M.trace_of_product(N)
            -1629
            sage: (M*N).trace()
            -1629
        """
        if self._nrows != other._ncols or other._nrows != self._ncols:
            raise ArithmeticError, "incompatible dimensions"
        s = self._base_ring(0)
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                s += self.get_unsafe(i, j) * other.get_unsafe(j, i)
        return s

    #####################################################################################
    # Generic Hessenberg Form and charpoly algorithm
    #####################################################################################
    def hessenberg_form(self):
        """
        Return Hessenberg form of self.

        If the base ring is merely an integral domain (and not a field),
        the Hessenberg form will (in general) only be defined over the
        fraction field of the base ring.

        EXAMPLES::

            sage: A = matrix(ZZ,4,[2, 1, 1, -2, 2, 2, -1, -1, -1,1,2,3,4,5,6,7])
            sage: h = A.hessenberg_form(); h
            [    2  -7/2 -19/5    -2]
            [    2   1/2 -17/5    -1]
            [    0  25/4  15/2   5/2]
            [    0     0  58/5     3]
            sage: parent(h)
            Full MatrixSpace of 4 by 4 dense matrices over Rational Field
            sage: A.hessenbergize()
            Traceback (most recent call last):
            ...
            TypeError: Hessenbergize only possible for matrices over a field
        """
        X = self.fetch('hessenberg_form')
        if not X is None:
            return X
        R = self._base_ring
        if not R.is_field():
            try:
                K = self._base_ring.fraction_field()
                H = self.change_ring(K)
                H.hessenbergize()
            except TypeError, msg:
                raise TypeError, "%s\nHessenberg form only possible for matrices over a field"%msg
        else:
            H = self.__copy__()
            H.hessenbergize()
        #end if
        self.cache('hessenberg_form', H)
        return H

    def hessenbergize(self):
        """
        Transform self to Hessenberg form.

        The hessenberg form of a matrix `A` is a matrix that is
        similar to `A`, so has the same characteristic polynomial
        as `A`, and is upper triangular except possible for entries
        right below the diagonal.

        ALGORITHM: See Henri Cohen's first book.

        EXAMPLES::

            sage: A = matrix(QQ,3, [2, 1, 1, -2, 2, 2, -1, -1, -1])
            sage: A.hessenbergize(); A
            [  2 3/2   1]
            [ -2   3   2]
            [  0  -3  -2]

        ::

            sage: A = matrix(QQ,4, [2, 1, 1, -2, 2, 2, -1, -1, -1,1,2,3,4,5,6,7])
            sage: A.hessenbergize(); A
            [    2  -7/2 -19/5    -2]
            [    2   1/2 -17/5    -1]
            [    0  25/4  15/2   5/2]
            [    0     0  58/5     3]

        You can't Hessenbergize an immutable matrix::

            sage: A = matrix(QQ, 3, [1..9])
            sage: A.set_immutable()
            sage: A.hessenbergize()
            Traceback (most recent call last):
            ...
            ValueError: matrix is immutable; please change a copy instead (i.e., use copy(M) to change a copy of M).
        """
        cdef Py_ssize_t i, j, m, n, r
        n = self._nrows

        tm = verbose("Computing Hessenberg Normal Form of %sx%s matrix"%(n,n))

        if not self.is_square():
            raise TypeError, "self must be square"

        if not self._base_ring.is_field():
            raise TypeError, "Hessenbergize only possible for matrices over a field"

        self.check_mutability()

        zero = self._base_ring(0)
        one = self._base_ring(1)
        for m from 1 <= m < n-1:
            # Search for a non-zero entry in column m-1
            i = -1
            for r from m+1 <= r < n:
                if self.get_unsafe(r, m-1) != zero:
                    i = r
                    break
            if i != -1:
                # Found a nonzero entry in column m-1 that is strictly below row m
                # Now set i to be the first nonzero position >= m in column m-1
                if self.get_unsafe(m,m-1) != zero:
                    i = m
                t = self.get_unsafe(i,m-1)
                t_inv = None
                if i > m:
                    self.swap_rows_c(i,m)
                    # We must do the corresponding column swap to
                    # maintain the characteristic polynomial (which is
                    # an invariant of Hessenberg form)
                    self.swap_columns_c(i,m)
                # Now the nonzero entry in position (m,m-1) is t.
                # Use t to clear the entries in column m-1 below m.
                for j from m+1 <= j < n:
                    x = self.get_unsafe(j, m-1)
                    if x != zero:
                        if t_inv is None:
                            t_inv = one / t
                        u = x * t_inv
                        self.add_multiple_of_row_c(j, m, -u, 0)
                        # To maintain charpoly, do the corresponding column operation,
                        # which doesn't mess up the matrix, since it only changes
                        # column m, and we're only worried about column m-1 right now.
                        # Add u*column_j to column_m.
                        self.add_multiple_of_column_c(m, j, u, 0)
        verbose("Finished Hessenberg Normal Form of %sx%s matrix"%(n,n),tm)



    def _charpoly_hessenberg(self, var):
        """
        Transforms self in place to its Hessenberg form then computes and
        returns the coefficients of the characteristic polynomial of this
        matrix.

        INPUT:

        -  ``var`` - name of the indeterminate of the charpoly

        The characteristic polynomial is represented as a vector of ints,
        where the constant term of the characteristic polynomial is the 0th
        coefficient of the vector.

        EXAMPLES::

            sage: matrix(QQ,3,range(9))._charpoly_hessenberg('Z')
            Z^3 - 12*Z^2 - 18*Z
            sage: matrix(ZZ,3,range(9))._charpoly_hessenberg('Z')
            Z^3 - 12*Z^2 - 18*Z
            sage: matrix(GF(7),3,range(9))._charpoly_hessenberg('Z')
            Z^3 + 2*Z^2 + 3*Z
            sage: matrix(QQ['x'],3,range(9))._charpoly_hessenberg('Z')
            Z^3 - 12*Z^2 - 18*Z
            sage: matrix(ZZ['ZZ'],3,range(9))._charpoly_hessenberg('Z')
            Z^3 - 12*Z^2 - 18*Z
        """
        if self._nrows != self._ncols:
            raise ArithmeticError, "charpoly not defined for non-square matrix."

        # Replace self by its Hessenberg form
        # (note the entries might now live in the fraction field)
        cdef Matrix H
        H = self.hessenberg_form()

        # We represent the intermediate polynomials that come up in
        # the calculations as rows of an (n+1)x(n+1) matrix, since
        # we've implemented basic arithmetic with such a matrix.
        # Please see the generic implementation of charpoly in
        # matrix.py to see more clearly how the following algorithm
        # actually works.  (The implementation is clearer (but slower)
        # if one uses polynomials to represent polynomials instead of
        # using the rows of a matrix.)  Also see Cohen's first GTM,
        # Algorithm 2.2.9.

        cdef Py_ssize_t i, m, n,
        n = self._nrows

        cdef Matrix c
        c = H.new_matrix(nrows=n+1,ncols=n+1)    # the 0 matrix
        one = H._coerce_element(1)
        c.set_unsafe(0,0,one)

        for m from 1 <= m <= n:
            # Set the m-th row of c to (x - H[m-1,m-1])*c[m-1] = x*c[m-1] - H[m-1,m-1]*c[m-1]
            # We do this by hand by setting the m-th row to c[m-1]
            # shifted to the right by one.  We then add
            # -H[m-1,m-1]*c[m-1] to the resulting m-th row.
            for i from 1 <= i <= n:
                c.set_unsafe(m, i, c.get_unsafe(m-1,i-1))
            c.add_multiple_of_row_c(m, m-1, -H.get_unsafe(m-1, m-1), 0)
            t = one
            for i from 1 <= i < m:
                t = t * H.get_unsafe(m-i,m-i-1)
                # Set the m-th row of c to c[m] - t*H[m-i-1,m-1]*c[m-i-1]
                c.add_multiple_of_row_c(m, m-i-1, - t*H.get_unsafe(m-i-1,m-1), 0)

        # The answer is now the n-th row of c.
        v = PyList_New(n+1)     # this is really sort of v = []..."
        for i from 0 <= i <= n:
            # Finally, set v[i] = c[n,i]
            o = c.get_unsafe(n,i)
            Py_INCREF(o); PyList_SET_ITEM(v, i, o)

        R = self._base_ring[var]    # polynomial ring over the base ring
        return R(v)

    #####################################################################################
    # Decomposition: kernel, image, decomposition
    #####################################################################################
    nullity = left_nullity

    def left_nullity(self):
        """
        Return the (left) nullity of this matrix, which is the dimension of
        the (left) kernel of this matrix acting from the right on row
        vectors.

        EXAMPLES::

            sage: M = Matrix(QQ,[[1,0,0,1],[0,1,1,0],[1,1,1,0]])
            sage: M.nullity()
            0
            sage: M.left_nullity()
            0

        ::

            sage: A = M.transpose()
            sage: A.nullity()
            1
            sage: A.left_nullity()
            1

        ::

            sage: M = M.change_ring(ZZ)
            sage: M.nullity()
            0
            sage: A = M.transpose()
            sage: A.nullity()
            1
        """
        # Use that rank + nullity = number of rows, since matrices act
        # from the right on row vectors.
        return self.nrows() - self.rank()

    def right_nullity(self):
        """
        Return the right nullity of this matrix, which is the dimension of
        the right kernel.

        EXAMPLES::

            sage: A = MatrixSpace(QQ,3,2)(range(6))
            sage: A.right_nullity()
            0

        ::

            sage: A = matrix(ZZ,3,range(9))
            sage: A.right_nullity()
            1
        """
        return self.ncols() - self.rank()

    ######################################
    # Kernel Helper Functions
    ######################################

    def _right_kernel_matrix_over_number_field(self):
        r"""
        Returns a pair that includes a matrix of basis vectors
        for the right kernel of ``self``.

        OUTPUT:

        Returns a pair.  First item is the string 'pivot-pari-numberfield'
        that identifies the nature of the basis vectors.

        Second item is a matrix whose rows are a basis for the right kernel,
        over the number field, as computed by PARI.

        EXAMPLES::

            sage: Q = QuadraticField(-7)
            sage: a = Q.gen(0)
            sage: A = matrix(Q, [[  2, 5-a, 15-a],
            ...                  [2+a,   a, -7 + 5*a]])
            sage: result = A._right_kernel_matrix_over_number_field()
            sage: result[0]
            'pivot-pari-numberfield'
            sage: P = result[1]; P
            [-a -3  1]
            sage: A*P.transpose() == zero_matrix(Q, 2, 1)
            True

        TESTS:

        We test some trivial cases. ::

            sage: Q = QuadraticField(-7)
            sage: A = matrix(Q, 0, 2)
            sage: A._right_kernel_matrix_over_number_field()[1]
            [1 0]
            [0 1]
            sage: A = matrix(Q, 2, 0)
            sage: A._right_kernel_matrix_over_number_field()[1].parent()
            Full MatrixSpace of 0 by 0 dense matrices over Number Field in a with defining polynomial x^2 + 7
            sage: A = zero_matrix(Q, 4, 3)
            sage: A._right_kernel_matrix_over_number_field()[1]
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        tm = verbose("computing right kernel matrix over a number field for %sx%s matrix" % (self.nrows(), self.ncols()),level=1)
        basis = self._pari_().matker()
        # Coerce PARI representations into the number field
        R = self.base_ring()
        basis = [[R(x) for x in row] for row in basis]
        verbose("done computing right kernel matrix over a number field for %sx%s matrix" % (self.nrows(), self.ncols()),level=1,t=tm)
        return 'pivot-pari-numberfield', matrix_space.MatrixSpace(R, len(basis), ncols=self._ncols)(basis)

    def _right_kernel_matrix_over_field(self, *args, **kwds):
        r"""
        Returns a pair that includes a matrix of basis vectors
        for the right kernel of ``self``.

        OUTPUT:

        Returns a pair.  First item is the string 'pivot-generic'
        that identifies the nature of the basis vectors.

        Second item is a matrix whose rows are a basis for the right kernel,
        over the field, as computed by general Python code.

        EXAMPLES::

            sage: C = CyclotomicField(14)
            sage: a = C.gen(0)
            sage: A = matrix(C, 3, 4, [[  1,    a,    1+a,  a^3+a^5],
            ...                        [  a,  a^4,  a+a^4,  a^4+a^8],
            ...                        [a^2, a^6, a^2+a^6, a^5+a^10]])
            sage: result = A._right_kernel_matrix_over_field()
            sage: result[0]
            'pivot-generic'
            sage: P = result[1]; P
            [       -1        -1         1         0]
            [-zeta14^3 -zeta14^4         0         1]
            sage: A*P.transpose() == zero_matrix(C, 3, 2)
            True

        TESTS:

        We test some trivial cases. ::

            sage: C = CyclotomicField(14)
            sage: A = matrix(C, 0, 2)
            sage: A._right_kernel_matrix_over_field()[1]
            [1 0]
            [0 1]
            sage: A = matrix(C, 2, 0)
            sage: A._right_kernel_matrix_over_field()[1].parent()
            Full MatrixSpace of 0 by 0 dense matrices over Cyclotomic Field of order 14 and degree 6
            sage: A = zero_matrix(C, 4, 3)
            sage: A._right_kernel_matrix_over_field()[1]
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        tm = verbose("computing right kernel matrix over an arbitrary field for %sx%s matrix" % (self.nrows(), self.ncols()),level=1)
        E = self.echelon_form(*args, **kwds)
        pivots = E.pivots()
        pivots_set = set(pivots)
        R = self.base_ring()
        zero = R(0)
        one = R(1)
        basis = []
        for i in xrange(self._ncols):
            if not (i in pivots_set):
                v = [zero]*self._ncols
                v[i] = one
                for r in range(len(pivots)):
                    v[pivots[r]] = -E[r,i]
                basis.append(v)
        tm = verbose("done computing right kernel matrix over an arbitrary field for %sx%s matrix" % (self.nrows(), self.ncols()),level=1,t=tm)
        return 'pivot-generic', matrix_space.MatrixSpace(R, len(basis), self._ncols)(basis)

    def _right_kernel_matrix_over_domain(self):
        r"""
        Returns a pair that includes a matrix of basis vectors
        for the right kernel of ``self``.

        OUTPUT:

        Returns a pair.  First item is the string 'computed-smith-form'
        that identifies the nature of the basis vectors.

        Second item is a matrix whose rows are a basis for the right kernel,
        over the field, as computed by general Python code.

        .. warning::

            This routine uses Smith normal form, which can fail
            if the domain is not a principal ideal domain.  Since we do
            not have a good test for PIDs, this is just a risk we take.
            See an example failure in the documentation for
            :meth:`right_kernel_matrix`.

        EXAMPLES:

        Univariate polynomials over a field form a PID.  ::

            sage: R.<y> = QQ[]
            sage: A = matrix(R, [[  1,   y, 1+y^2],
            ...                  [y^3, y^2, 2*y^3]])
            sage: result = A._right_kernel_matrix_over_domain()
            sage: result[0]
            'computed-smith-form'
            sage: P = result[1]; P
            [-1 -y  1]
            sage: A*P.transpose() == zero_matrix(R, 2, 1)
            True

        TESTS:

        We test some trivial cases. ::

            sage: R.<y> = QQ[]
            sage: A = matrix(R, 0, 2)
            sage: A._right_kernel_matrix_over_domain()[1]
            [1 0]
            [0 1]
            sage: A = matrix(R, 2, 0)
            sage: A._right_kernel_matrix_over_domain()[1].parent()
            Full MatrixSpace of 0 by 0 dense matrices over Univariate Polynomial Ring in y over Rational Field
            sage: A = zero_matrix(R, 4, 3)
            sage: A._right_kernel_matrix_over_domain()[1]
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        tm = verbose("computing right kernel matrix over a domain for %sx%s matrix" % (self.nrows(), self.ncols()),level=1)
        d, u, v = self.smith_form()
        basis = []
        for i in xrange(self.ncols()):
            if (i >= self.nrows()) or d[i][i] == 0:
                basis.append( v.column(i).list() )
        verbose("done computing right kernel matrix over a domain for %sx%s matrix" % (self.nrows(), self.ncols()),level=1,t=tm)
        return 'computed-smith-form', self.new_matrix(nrows=len(basis), ncols=self._ncols, entries=basis)

    def right_kernel_matrix(self, *args, **kwds):
        r"""
        Returns a matrix whose rows form a basis
        for the right kernel of ``self``.

        INPUT:

        - ``algorithm`` - default: 'default' - a keyword that selects the
          algorithm employed.  Allowable values are:

          - 'default' - allows the algorithm to be chosen automatically
          - 'generic' - naive algorithm usable for matrices over any field
          - 'pari' - PARI library code for matrices over number fields
            or the integers
          - 'padic' - padic algorithm from IML library for matrices
            over the rationals and integers
          - 'pluq' - PLUQ matrix factorization for matrices mod 2

        - ``basis`` - default: 'echelon' - a keyword that describes
          the format of the basis returned.  Allowable values are:

          - 'echelon': the basis matrix is in echelon form
          - 'pivot' : each basis vector is computed from the reduced
            row-echelon form of ``self`` by placing a single one in a
            non-pivot column and zeros in the remaining non-pivot columns.
            Only available for matrices over fields.
          - 'LLL': an LLL-reduced basis.  Only available for matrices
            over the integers.
          - 'computed': no work is done to transform the basis, it is
            returned exactly as provided by whichever routine actually
            computed the basis.  Request this for the least possible
            computation possible, but with no guarantees about the format
            of the basis.

        OUTPUT:

        A matrix ``X``  whose rows are an independent set spanning the
        right kernel of ``self``.  So ``self*X.transpose()`` is a zero matrix.

        The output varies depending on the choice of ``algorithm`` and the
        format chosen by ``basis``.

        The results of this routine are not cached, so you can call it again
        with different options to get possibly different output (like the basis
        format).  Conversely, repeated calls on the same matrix will always
        start from scratch.

        .. note::

            If you want to get the most basic description of a kernel, with a
            minimum of overhead, then ask for the right kernel matrix with
            the basis format requested as 'computed'.  You are then free to
            work with the output for whatever purpose.  For a left kernel,
            call this method on the transpose of your matrix.

            For greater convenience, plus cached results, request an actual
            vector space or free module with :meth:`right_kernel`
            or :meth:`left_kernel`.

        EXAMPLES:

        Over the Rational Numbers:

        Kernels are computed by the IML library in
        :meth:`~sage.matrix.matrix_rational_dense.Matrix_rational_dense._right_kernel_matrix`.
        Setting the `algorithm` keyword to 'default', 'padic' or unspecified
        will yield the same result, as there is no optional behavior.
        The 'computed' format of the basis vectors are exactly the negatives
        of the vectors in the 'pivot' format. ::

            sage: A = matrix(QQ, [[1, 0, 1, -3, 1],
            ...                   [-5, 1, 0, 7, -3],
            ...                   [0, -1, -4, 6, -2],
            ...                   [4, -1, 0, -6, 2]])
            sage: C = A.right_kernel_matrix(algorithm='default', basis='computed'); C
            [-1  2 -2 -1  0]
            [ 1  2  0  0 -1]
            sage: A*C.transpose() == zero_matrix(QQ, 4, 2)
            True
            sage: P = A.right_kernel_matrix(algorithm='padic', basis='pivot'); P
            [ 1 -2  2  1  0]
            [-1 -2  0  0  1]
            sage: A*P.transpose() == zero_matrix(QQ, 4, 2)
            True
            sage: C == -P
            True
            sage: E = A.right_kernel_matrix(algorithm='default', basis='echelon'); E
            [   1    0    1  1/2 -1/2]
            [   0    1 -1/2 -1/4 -1/4]
            sage: A*E.transpose() == zero_matrix(QQ, 4, 2)
            True

        Since the rationals are a field, we can call the general code
        available for any field by using the 'generic' keyword. ::

            sage: A = matrix(QQ, [[1, 0, 1, -3, 1],
            ...                   [-5, 1, 0, 7, -3],
            ...                   [0, -1, -4, 6, -2],
            ...                   [4, -1, 0, -6, 2]])
            sage: G = A.right_kernel_matrix(algorithm='generic', basis='echelon'); G
            [   1    0    1  1/2 -1/2]
            [   0    1 -1/2 -1/4 -1/4]
            sage: A*G.transpose() == zero_matrix(QQ, 4, 2)
            True

        We verify that the rational matrix code is called for both
        dense and sparse rational matrices, with equal result. ::

            sage: A = matrix(QQ, [[1, 0, 1, -3, 1],
            ...                   [-5, 1, 0, 7, -3],
            ...                   [0, -1, -4, 6, -2],
            ...                   [4, -1, 0, -6, 2]],
            ...              sparse=False)
            sage: B = copy(A).sparse_matrix()
            sage: set_verbose(1)
            sage: D = A.right_kernel(); D
            verbose ...
            verbose 1 (<module>) computing right kernel matrix over the rationals for 4x5 matrix
            ...
            verbose 1 (<module>) done computing right kernel matrix over the rationals for 4x5 matrix
            ...
            Vector space of degree 5 and dimension 2 over Rational Field
            Basis matrix:
            [   1    0    1  1/2 -1/2]
            [   0    1 -1/2 -1/4 -1/4]
            sage: S = B.right_kernel(); S
            verbose ...
            verbose 1 (<module>) computing right kernel matrix over the rationals for 4x5 matrix
            ...
            verbose 1 (<module>) done computing right kernel matrix over the rationals for 4x5 matrix
            ...
            Vector space of degree 5 and dimension 2 over Rational Field
            Basis matrix:
            [   1    0    1  1/2 -1/2]
            [   0    1 -1/2 -1/4 -1/4]
            sage: set_verbose(0)
            sage: D == S
            True

        Over Number Fields:

        Kernels are by default computed by PARI, (except for exceptions like
        the rationals themselves).  The raw results from PARI are a pivot
        basis, so the `basis` keywords 'computed' and 'pivot' will return
        the same results. ::

            sage: Q = QuadraticField(-7)
            sage: a = Q.gen(0)
            sage: A = matrix(Q, [[2, 5-a, 15-a, 16+4*a],
            ...                  [2+a, a, -7 + 5*a, -3+3*a]])
            sage: C = A.right_kernel_matrix(algorithm='default', basis='computed'); C
            [    -a     -3      1      0]
            [    -2 -a - 1      0      1]
            sage: A*C.transpose() == zero_matrix(Q, 2, 2)
            True
            sage: P = A.right_kernel_matrix(algorithm='pari', basis='pivot'); P
            [    -a     -3      1      0]
            [    -2 -a - 1      0      1]
            sage: A*P.transpose() == zero_matrix(Q, 2, 2)
            True
            sage: E = A.right_kernel_matrix(algorithm='default', basis='echelon'); E
            [                1                 0     7/88*a + 3/88 -3/176*a - 39/176]
            [                0                 1   -1/88*a - 13/88  13/176*a - 7/176]
            sage: A*E.transpose() == zero_matrix(Q, 2, 2)
            True

        We can bypass using PARI for number fields and use Sage's general
        code for matrices over any field.  The basis vectors as computed
        are in pivot format. ::

            sage: Q = QuadraticField(-7)
            sage: a = Q.gen(0)
            sage: A = matrix(Q, [[2, 5-a, 15-a, 16+4*a],[2+a, a, -7 + 5*a, -3+3*a]])
            sage: G = A.right_kernel_matrix(algorithm='generic', basis='computed'); G
            [    -a     -3      1      0]
            [    -2 -a - 1      0      1]
            sage: A*G.transpose() == zero_matrix(Q, 2, 2)
            True

        We check that number fields are handled by the right routine as part of
        typical right kernel computation. ::

            sage: Q = QuadraticField(-7)
            sage: a = Q.gen(0)
            sage: A = matrix(Q, [[2, 5-a, 15-a, 16+4*a],[2+a, a, -7 + 5*a, -3+3*a]])
            sage: set_verbose(1)
            sage: A.right_kernel(algorithm='default')
            verbose ...
            verbose 1 (<module>) computing right kernel matrix over a number field for 2x4 matrix
            verbose 1 (<module>) done computing right kernel matrix over a number field for 2x4 matrix
            ...
            Vector space of degree 4 and dimension 2 over Number Field in a with defining polynomial x^2 + 7
            Basis matrix:
            [                1                 0     7/88*a + 3/88 -3/176*a - 39/176]
            [                0                 1   -1/88*a - 13/88  13/176*a - 7/176]
            sage: set_verbose(0)

        Over the Finite Field of Order 2:

        Kernels are computed by the M4RI library using PLUQ matrix
        decomposition in the
        :meth:`~sage.matrix.matrix_mod2_dense.Matrix_mod2_dense._right_kernel_matrix`
        method. There are no options for the algorithm used.  ::

            sage: A = matrix(GF(2),[[0, 1, 1, 0, 0, 0],
            ...                     [1, 0, 0, 0, 1, 1,],
            ...                     [1, 0, 0, 0, 1, 1]])
            sage: E = A.right_kernel_matrix(algorithm='default', format='echelon'); E
            [1 0 0 0 0 1]
            [0 1 1 0 0 0]
            [0 0 0 1 0 0]
            [0 0 0 0 1 1]
            sage: A*E.transpose() == zero_matrix(GF(2), 3, 4)
            True

        Since GF(2) is a field we can route this computation to the generic
        code and obtain the 'pivot' form of the basis.  The ``algorithm``
        keywords, 'pluq', 'default' and unspecified, all have the
        same effect as there is no optional behavior. ::

            sage: A = matrix(GF(2),[[0, 1, 1, 0, 0, 0],
            ...                     [1, 0, 0, 0, 1, 1,],
            ...                     [1, 0, 0, 0, 1, 1]])
            sage: P = A.right_kernel_matrix(algorithm='generic', basis='pivot'); P
            [0 1 1 0 0 0]
            [0 0 0 1 0 0]
            [1 0 0 0 1 0]
            [1 0 0 0 0 1]
            sage: A*P.transpose() == zero_matrix(GF(2), 3, 4)
            True
            sage: DP = A.right_kernel_matrix(algorithm='default', basis='pivot'); DP
            [0 1 1 0 0 0]
            [0 0 0 1 0 0]
            [1 0 0 0 1 0]
            [1 0 0 0 0 1]
            sage: A*DP.transpose() == zero_matrix(GF(2), 3, 4)
            True
            sage: A.right_kernel_matrix(algorithm='pluq', basis='echelon')
            [1 0 0 0 0 1]
            [0 1 1 0 0 0]
            [0 0 0 1 0 0]
            [0 0 0 0 1 1]

        We test that the mod 2 code is called for matrices over GF(2). ::

            sage: A = matrix(GF(2),[[0, 1, 1, 0, 0, 0],
            ...                     [1, 0, 0, 0, 1, 1,],
            ...                     [1, 0, 0, 0, 1, 1]])
            sage: set_verbose(1)
            sage: A.right_kernel(algorithm='default')
            verbose ...
            verbose 1 (<module>) computing right kernel matrix over integers mod 2 for 3x6 matrix
            verbose 1 (<module>) done computing right kernel matrix over integers mod 2 for 3x6 matrix
            ...
            Vector space of degree 6 and dimension 4 over Finite Field of size 2
            Basis matrix:
            [1 0 0 0 0 1]
            [0 1 1 0 0 0]
            [0 0 0 1 0 0]
            [0 0 0 0 1 1]
            sage: set_verbose(0)

        Over Arbitrary Fields:

        For kernels over fields not listed above, totally general code
        will compute a set of basis vectors in the pivot format.
        These could be returned as a basis in echelon form.  ::

            sage: F.<a> = FiniteField(5^2)
            sage: A = matrix(F, 3, 4, [[  1,   a,     1+a,  a^3+a^5],
            ...                        [  a, a^4,   a+a^4,  a^4+a^8],
            ...                        [a^2, a^6, a^2+a^6, a^5+a^10]])
            sage: P = A.right_kernel_matrix(algorithm='default', basis='pivot'); P
            [      4       4       1       0]
            [  a + 2 3*a + 3       0       1]
            sage: A*P.transpose() == zero_matrix(F, 3, 2)
            True
            sage: E = A.right_kernel_matrix(algorithm='default', basis='echelon'); E
            [      1       0 3*a + 4 2*a + 2]
            [      0       1     2*a 3*a + 3]
            sage: A*E.transpose() == zero_matrix(F, 3, 2)
            True

        This general code can be requested for matrices over any field
        with the ``algorithm`` keyword 'generic'.  Normally, matrices
        over the rationals would be handled by specific routines from
        the IML library. The default format is an echelon basis, but a
        pivot basis may be requested, which is identical to the computed
        basis. ::

            sage: A = matrix(QQ, 3, 4, [[1,3,-2,4],
            ...                         [2,0,2,2],
            ...                         [-1,1,-2,0]])
            sage: G = A.right_kernel_matrix(algorithm='generic'); G
            [   1    0 -1/2 -1/2]
            [   0    1  1/2 -1/2]
            sage: A*G.transpose() == zero_matrix(QQ, 3, 2)
            True
            sage: C = A.right_kernel_matrix(algorithm='generic', basis='computed'); C
            [-1  1  1  0]
            [-1 -1  0  1]
            sage: A*C.transpose() == zero_matrix(QQ, 3, 2)
            True

        We test that the generic code is called for matrices over fields,
        lacking any more specific routine. ::

            sage: F.<a> = FiniteField(5^2)
            sage: A = matrix(F, 3, 4, [[  1,   a,     1+a,  a^3+a^5],
            ...                        [  a, a^4,   a+a^4,  a^4+a^8],
            ...                        [a^2, a^6, a^2+a^6, a^5+a^10]])
            sage: set_verbose(1)
            sage: A.right_kernel(algorithm='default')
            verbose ...
            verbose 1 (<module>) computing right kernel matrix over an arbitrary field for 3x4 matrix
            ...
            verbose 1 (<module>) done computing right kernel matrix over an arbitrary field for 3x4 matrix
            ...
            Vector space of degree 4 and dimension 2 over Finite Field in a of size 5^2
            Basis matrix:
            [      1       0 3*a + 4 2*a + 2]
            [      0       1     2*a 3*a + 3]
            sage: set_verbose(0)

        Over the Integers:

        Either the IML or PARI libraries are used to provide a set of basis
        vectors.  The ``algorithm`` keyword can be used to select either,
        or when set to 'default' a heuristic will choose between the two.
        Results can be returned in the 'compute' format, straight out of
        the libraries.  Unique to the integers, the basis vectors can be
        returned as an LLL basis.  Note the similarities and differences
        in the results. The 'pivot' format is not available, since the
        integers are not a field.  ::

            sage: A = matrix(ZZ, [[8, 0, 7, 1, 3, 4, 6],
            ...                   [4, 0, 3, 4, 2, 7, 7],
            ...                   [1, 4, 6, 1, 2, 8, 5],
            ...                   [0, 3, 1, 2, 3, 6, 2]])

            sage: X = A.right_kernel_matrix(algorithm='default', basis='echelon'); X
            [  1  12   3  14  -3 -10   1]
            [  0  35   0  25  -1 -31  17]
            [  0   0   7  12  -3  -1  -8]
            sage: A*X.transpose() == zero_matrix(ZZ, 4, 3)
            True

            sage: X = A.right_kernel_matrix(algorithm='padic', basis='LLL'); X
            [ -3  -1   5   7   2  -3  -2]
            [  3   1   2   5  -5   2  -6]
            [ -4 -13   2  -7   5   7  -3]
            sage: A*X.transpose() == zero_matrix(ZZ, 4, 3)
            True

            sage: X = A.right_kernel_matrix(algorithm='pari', basis='computed'); X
            [ 3  1 -5 -7 -2  3  2]
            [ 3  1  2  5 -5  2 -6]
            [ 4 13 -2  7 -5 -7  3]
            sage: A*X.transpose() == zero_matrix(ZZ, 4, 3)
            True

            sage: X = A.right_kernel_matrix(algorithm='padic', basis='computed'); X
            [ 265  345 -178   17 -297    0    0]
            [-242 -314  163  -14  271   -1    0]
            [ -36  -47   25   -1   40    0   -1]
            sage: A*X.transpose() == zero_matrix(ZZ, 4, 3)
            True

        We test that the code for integer matrices is called for matrices
        defined over the integers, both dense and sparse, with equal result. ::

            sage: A = matrix(ZZ, [[8, 0, 7, 1, 3, 4, 6],
            ...                   [4, 0, 3, 4, 2, 7, 7],
            ...                   [1, 4, 6, 1, 2, 8, 5],
            ...                   [0, 3, 1, 2, 3, 6, 2]],
            ...              sparse=False)
            sage: B = copy(A).sparse_matrix()
            sage: set_verbose(1)
            sage: D = A.right_kernel(); D
            verbose ...
            verbose 1 (<module>) computing right kernel matrix over the integers for 4x7 matrix
            verbose 1 (<module>) done computing right kernel matrix over the integers for 4x7 matrix
            ...
            Free module of degree 7 and rank 3 over Integer Ring
            Echelon basis matrix:
            [  1  12   3  14  -3 -10   1]
            [  0  35   0  25  -1 -31  17]
            [  0   0   7  12  -3  -1  -8]
            sage: S = B.right_kernel(); S
            verbose ...
            verbose 1 (<module>) computing right kernel matrix over the integers for 4x7 matrix
            verbose 1 (<module>) done computing right kernel matrix over the integers for 4x7 matrix
            ...
            Free module of degree 7 and rank 3 over Integer Ring
            Echelon basis matrix:
            [  1  12   3  14  -3 -10   1]
            [  0  35   0  25  -1 -31  17]
            [  0   0   7  12  -3  -1  -8]
            sage: set_verbose(0)
            sage: D == S
            True

        Over Principal Ideal Domains:

        Kernels can be computed using Smith normal form.  Only the default
        algorithm is available, and the 'pivot' basis format is
        not available. ::

            sage: R.<y> = QQ[]
            sage: A = matrix(R, [[  1,   y, 1+y^2],
            ...                  [y^3, y^2, 2*y^3]])
            sage: E = A.right_kernel_matrix(algorithm='default', basis='echelon'); E
            [-1 -y  1]
            sage: A*E.transpose() == zero_matrix(ZZ, 2, 1)
            True

        It can be computationally expensive to determine if an integral
        domain is a principal ideal domain.  The Smith normal form routine
        can fail for non-PIDs, as in this example. ::

            sage: D.<x> = ZZ[]
            sage: A = matrix(D, 2, 2, [[x^2 - x, -x + 5],
            ...                        [x^2 - 8, -x + 2]])
            sage: A.right_kernel_matrix()
            Traceback (most recent call last):
            ...
            ArithmeticError: Ideal Ideal (x^2 - x, x^2 - 8) of Univariate Polynomial Ring in x over Integer Ring not principal

        We test that the domain code is called for domains that lack any
        extra structure. ::

            sage: R.<y> = QQ[]
            sage: A = matrix(R, [[  1,   y, 1+y^2],
            ...                  [y^3, y^2, 2*y^3]])
            sage: set_verbose(1)
            sage: A.right_kernel(algorithm='default', basis='echelon')
            verbose ...
            verbose 1 (<module>) computing right kernel matrix over a domain for 2x3 matrix
            verbose 1 (<module>) done computing right kernel matrix over a domain for 2x3 matrix
            ...
            Free module of degree 3 and rank 1 over Univariate Polynomial Ring in y over Rational Field
            Echelon basis matrix:
            [-1 -y  1]
            sage: set_verbose(0)

        Trivial Cases:

        We test two trivial cases.  Any possible values for the
        keywords (``algorithm``, ``basis``) will return identical results. ::

            sage: A = matrix(ZZ, 0, 2)
            sage: A.right_kernel_matrix()
            [1 0]
            [0 1]
            sage: A = matrix(FiniteField(7), 2, 0)
            sage: A.right_kernel_matrix().parent()
            Full MatrixSpace of 0 by 0 dense matrices over Finite Field of size 7

        TESTS:

        The "usual" quaternions are a non-commutative ring and computations
        of kernels over these rings are not implemented. ::

            sage: Q.<i,j,k> = QuaternionAlgebra(-1,-1)
            sage: A = matrix(Q, 2, [i,j,-1,k])
            sage: A.right_kernel_matrix()
            Traceback (most recent call last):
            ...
            NotImplementedError: Cannot compute a matrix kernel over Quaternion Algebra (-1, -1) with base ring Rational Field

        We test error messages for improper choices of the 'algorithm'
        keyword. ::

            sage: matrix(ZZ, 2, 2).right_kernel_matrix(algorithm='junk')
            Traceback (most recent call last):
            ...
            ValueError: matrix kernel algorithm 'junk' not recognized
            sage: matrix(GF(2), 2, 2).right_kernel_matrix(algorithm='padic')
            Traceback (most recent call last):
            ...
            ValueError: 'padic' matrix kernel algorithm only available over the rationals and the integers, not over Finite Field of size 2
            sage: matrix(QQ, 2, 2).right_kernel_matrix(algorithm='pari')
            Traceback (most recent call last):
            ...
            ValueError: 'pari' matrix kernel algorithm only available over non-trivial number fields and the integers, not over Rational Field
            sage: matrix(Integers(6), 2, 2).right_kernel_matrix(algorithm='generic')
            Traceback (most recent call last):
            ...
            ValueError: 'generic' matrix kernel algorithm only available over a field, not over Ring of integers modulo 6
            sage: matrix(QQ, 2, 2).right_kernel_matrix(algorithm='pluq')
            Traceback (most recent call last):
            ...
            ValueError: 'pluq' matrix kernel algorithm only available over integers mod 2, not over Rational Field

        We test error messages for improper basis format requests. ::

            sage: matrix(ZZ, 2, 2).right_kernel_matrix(basis='junk')
            Traceback (most recent call last):
            ...
            ValueError: matrix kernel basis format 'junk' not recognized
            sage: matrix(ZZ, 2, 2).right_kernel_matrix(basis='pivot')
            Traceback (most recent call last):
            ...
            ValueError: pivot basis only available over a field, not over Integer Ring
            sage: matrix(QQ, 2, 2).right_kernel_matrix(basis='LLL')
            Traceback (most recent call last):
            ...
            ValueError: LLL-reduced basis only available over the integers, not over Rational Field

        Finally, error messages for the 'proof' keyword.  ::

            sage: matrix(ZZ, 2, 2).right_kernel_matrix(proof='junk')
            Traceback (most recent call last):
            ...
            ValueError: 'proof' must be one of True, False or None, not junk
            sage: matrix(QQ, 2, 2).right_kernel_matrix(proof=True)
            Traceback (most recent call last):
            ...
            ValueError: 'proof' flag only valid for matrices over the integers

        AUTHOR:

        - Rob Beezer (2011-02-05)
        """
        R = self.base_ring()

        # First: massage keywords
        # Determine algorithm to use for computation of kernel matrix
        algorithm = kwds.pop('algorithm', None)
        if algorithm is None:
            algorithm = 'default'
        elif not algorithm in ['default', 'generic', 'pari', 'padic', 'pluq']:
            raise ValueError("matrix kernel algorithm '%s' not recognized" % algorithm )
        elif algorithm == 'padic' and not (is_IntegerRing(R) or is_RationalField(R)):
            raise ValueError("'padic' matrix kernel algorithm only available over the rationals and the integers, not over %s" % R)
        elif algorithm == 'pari' and not (is_IntegerRing(R) or (is_NumberField(R) and not is_RationalField(R))):
            raise ValueError("'pari' matrix kernel algorithm only available over non-trivial number fields and the integers, not over %s" % R)
        elif algorithm == 'generic' and not R.is_field():
            raise ValueError("'generic' matrix kernel algorithm only available over a field, not over %s" % R)
        elif algorithm == 'pluq' and not isinstance(self, sage.matrix.matrix_mod2_dense.Matrix_mod2_dense):
            raise ValueError("'pluq' matrix kernel algorithm only available over integers mod 2, not over %s" % R)

        # Determine the basis format of independent spanning set to return
        basis = kwds.pop('basis', None)
        if basis is None:
            basis = 'echelon'
        elif not basis in ['computed', 'echelon', 'pivot', 'LLL']:
            raise ValueError("matrix kernel basis format '%s' not recognized" % basis )
        elif basis == 'pivot' and not R.is_field():
            raise ValueError('pivot basis only available over a field, not over %s' % R)
        elif basis == 'LLL' and not is_IntegerRing(R):
            raise ValueError('LLL-reduced basis only available over the integers, not over %s' % R)

        # Determine proof keyword for integer matrices
        proof = kwds.pop('proof', None)
        if not (proof in [None, True, False]):
            raise ValueError("'proof' must be one of True, False or None, not %s" % proof)
        if not (proof is None or is_IntegerRing(R)):
            raise ValueError("'proof' flag only valid for matrices over the integers")

        # We could sanitize/process remaining (un-popped) keywords here and
        # send them to the echelon form routine in the 'generic' algorithm case
        # Would need to handle 'algorithm' keyword properly

        # Second: Trivial cases over any integral domain
        # With zero columns the domain has dimension 0,
        #   so only an empty basis is possible
        # With zero rows the codomain is just the zero vector,
        #   thus the kernel is the whole domain, so return an identity matrix
        #   identity_matrix constructor will fail if ring does not have a one
        # For all keywords the results are identical
        if self._ncols == 0 and R.is_integral_domain():
            return self.new_matrix(nrows = 0, ncols = self._ncols)
        if self._nrows == 0 and R.is_integral_domain():
            import constructor
            return constructor.identity_matrix(R, self._ncols)

        # Third: generic first, if requested explicitly
        #   then try specialized class methods, and finally
        #   delegate to ad-hoc methods in greater generality
        M = None; format = ''

        if algorithm == 'generic':
            format, M = self._right_kernel_matrix_over_field()

        if M is None:
            try:
                format, M = self._right_kernel_matrix(algorithm=algorithm, proof=proof)
            except AttributeError:
                pass

        if M is None and is_NumberField(R):
            format, M = self._right_kernel_matrix_over_number_field()

        if M is None and R.is_field():
            format, M = self._right_kernel_matrix_over_field()

        if M is None and R.is_integral_domain():
            format, M = self._right_kernel_matrix_over_domain()

        if M is None:
            raise NotImplementedError("Cannot compute a matrix kernel over %s" % R)

        # Trivial kernels give empty matrices, which sometimes mistakenly have
        #   zero columns as well. (eg PARI?)  This could be fixed at the source
        #   with a careful study of the phenomenon.  Start by commenting out
        #   the following and running doctests in sage/matrix
        if M.nrows()==0 and M.ncols()!=self.ncols():
            M = M.new_matrix(nrows=0, ncols=self.ncols())

        # Convert basis to requested type and return the matrix
        #   format string leads with 'echelon', 'pivot' or 'LLL' if known
        #   to be of that format otherwise format string leads with
        #   'computed' if it needs work (ie raw results)
        if basis == 'computed':
            return M
        elif basis == 'echelon':
            if not format[:7] == 'echelon':
                return M.echelon_form()
            else:
                return M
        elif basis == 'pivot':
            # cannot get here unless over a field
            if not format[:5] == 'pivot':
                # convert non-pivot columns to identity matrix
                # this is the basis immediately obvious from echelon form
                # so C is always invertible (when working over a field)
                C = M.matrix_from_columns(self.nonpivots())
                return C.inverse()*M
            else:
                return M
        elif basis == 'LLL':
            # cannot get here unless over integers
            if not format[:3] == 'LLL':
                return M.LLL()
            else:
                return M

    def right_kernel(self, *args, **kwds):
        r"""
        Returns the right kernel of this matrix, as a vector space or
        free module. This is the set of vectors ``x`` such that ``self*x = 0``.

        .. note::

            For the left kernel, use :meth:`left_kernel`.  The method
            :meth:`kernel` is exactly equal to :meth:`left_kernel`.

        INPUT:

        - ``algorithm`` - default: 'default' - a keyword that selects the
          algorithm employed.  Allowable values are:

          - 'default' - allows the algorithm to be chosen automatically
          - 'generic' - naive algorithm usable for matrices over any field
          - 'pari' - PARI library code for matrices over number fields
            or the integers
          - 'padic' - padic algorithm from IML library for matrices
            over the rationals and integers
          - 'pluq' - PLUQ matrix factorization for matrices mod 2

        - ``basis`` - default: 'echelon' - a keyword that describes the
          format of the basis used to construct the left kernel.
          Allowable values are:

          - 'echelon': the basis matrix is in echelon form
          - 'pivot' : each basis vector is computed from the reduced
            row-echelon form of ``self`` by placing a single one in a
            non-pivot column and zeros in the remaining non-pivot columns.
            Only available for matrices over fields.
          - 'LLL': an LLL-reduced basis.  Only available for matrices
            over the integers.

        OUTPUT:

        A vector space or free module whose degree equals the number
        of columns in ``self`` and contains all the vectors ``x`` such
        that ``self*x = 0``.

        If ``self`` has 0 columns, the kernel has dimension 0, while if
        ``self`` has 0 rows the kernel is the entire ambient vector space.

        The result is cached.  Requesting the right kernel a second time,
        but with a different basis format, will return the cached result
        with the format from the first computation.

        .. note::

           For more detailed documentation on the selection of algorithms
           used and a more flexible method for computing a basis matrix
           for a right kernel (rather than computing a vector space), see
           :meth:`right_kernel_matrix`, which powers the computations for
           this method.

        EXAMPLES::

            sage: A = matrix(QQ, [[0, 0, 1, 2, 2, -5, 3],
            ...                   [-1, 5, 2, 2, 1, -7, 5],
            ...                   [0, 0, -2, -3, -3, 8, -5],
            ...                   [-1, 5, 0, -1, -2, 1, 0]])
            sage: K = A.right_kernel(); K
            Vector space of degree 7 and dimension 4 over Rational Field
            Basis matrix:
            [ 1  0  0  0 -1 -1 -1]
            [ 0  1  0  0  5  5  5]
            [ 0  0  1  0 -1 -2 -3]
            [ 0  0  0  1  0  1  1]
            sage: A*K.basis_matrix().transpose() == zero_matrix(QQ, 4, 4)
            True

        The default is basis vectors that form a matrix in echelon form.
        A "pivot basis" instead has a basis matrix where the columns of
        an identity matrix are in the locations of the non-pivot columns
        of the original matrix. This alternate format is available whenever
        the base ring is a field. ::

            sage: A = matrix(QQ, [[0, 0, 1, 2, 2, -5, 3],
            ...                   [-1, 5, 2, 2, 1, -7, 5],
            ...                   [0, 0, -2, -3, -3, 8, -5],
            ...                   [-1, 5, 0, -1, -2, 1, 0]])
            sage: A.rref()
            [ 1 -5  0  0  1  1 -1]
            [ 0  0  1  0  0 -1  1]
            [ 0  0  0  1  1 -2  1]
            [ 0  0  0  0  0  0  0]
            sage: A.nonpivots()
            (1, 4, 5, 6)
            sage: K = A.right_kernel(basis='pivot'); K
            Vector space of degree 7 and dimension 4 over Rational Field
            User basis matrix:
            [ 5  1  0  0  0  0  0]
            [-1  0  0 -1  1  0  0]
            [-1  0  1  2  0  1  0]
            [ 1  0 -1 -1  0  0  1]
            sage: A*K.basis_matrix().transpose() == zero_matrix(QQ, 4, 4)
            True

        Matrices may have any field as a base ring.  Number fields are
        computed by PARI library code, matrices over `GF(2)` are computed
        by the M4RI library, and matrices over the rationals are computed by
        the IML library.  For any of these specialized cases, general-purpose
        code can be called instead with the keyword setting
        ``algorithm='generic'``.

        Over an arbitrary field, with two basis formats.  Same vector space,
        different bases.  ::

            sage: F.<a> = FiniteField(5^2)
            sage: A = matrix(F, 3, 4, [[  1,   a,     1+a,  a^3+a^5],
            ...                        [  a, a^4,   a+a^4,  a^4+a^8],
            ...                        [a^2, a^6, a^2+a^6, a^5+a^10]])
            sage: K = A.right_kernel(); K
            Vector space of degree 4 and dimension 2 over Finite Field in a of size 5^2
            Basis matrix:
            [      1       0 3*a + 4 2*a + 2]
            [      0       1     2*a 3*a + 3]
            sage: A*K.basis_matrix().transpose() == zero_matrix(F, 3, 2)
            True
            sage: B = copy(A)
            sage: P = B.right_kernel(basis = 'pivot'); P
            Vector space of degree 4 and dimension 2 over Finite Field in a of size 5^2
            User basis matrix:
            [      4       4       1       0]
            [  a + 2 3*a + 3       0       1]
            sage: B*P.basis_matrix().transpose() == zero_matrix(F, 3, 2)
            True
            sage: K == P
            True

        Over number fields, PARI is used by default, but general-purpose code
        can be requested.  Same vector space, same bases, different code.::

            sage: Q = QuadraticField(-7)
            sage: a = Q.gen(0)
            sage: A = matrix(Q, [[  2, 5-a,     15-a, 16+4*a],
            ...                  [2+a,   a, -7 + 5*a, -3+3*a]])
            sage: K = A.right_kernel(algorithm='default'); K
            Vector space of degree 4 and dimension 2 over Number Field in a with defining polynomial x^2 + 7
            Basis matrix:
            [                1                 0     7/88*a + 3/88 -3/176*a - 39/176]
            [                0                 1   -1/88*a - 13/88  13/176*a - 7/176]
            sage: A*K.basis_matrix().transpose() == zero_matrix(Q, 2, 2)
            True
            sage: B = copy(A)
            sage: G = A.right_kernel(algorithm='generic'); G
            Vector space of degree 4 and dimension 2 over Number Field in a with defining polynomial x^2 + 7
            Basis matrix:
            [                1                 0     7/88*a + 3/88 -3/176*a - 39/176]
            [                0                 1   -1/88*a - 13/88  13/176*a - 7/176]
            sage: B*G.basis_matrix().transpose() == zero_matrix(Q, 2, 2)
            True
            sage: K == G
            True

        For matrices over the integers, several options are possible.
        The basis can be an LLL-reduced basis or an echelon basis.
        The pivot basis isnot available.  A heuristic will decide whether
        to use a p-adic algorithm from the IML library or an algorithm
        from the PARI library.  Note how specifying the algorithm can
        mildly influence the LLL basis. ::

            sage: A = matrix(ZZ, [[0, -1, -1, 2, 9, 4, -4],
            ...                   [-1, 1, 0, -2, -7, -1, 6],
            ...                   [2, 0, 1, 0, 1, -5, -2],
            ...                   [-1, -1, -1, 3, 10, 10, -9],
            ...                   [-1, 2, 0, -3, -7, 1, 6]])
            sage: A.right_kernel(basis='echelon')
            Free module of degree 7 and rank 2 over Integer Ring
            Echelon basis matrix:
            [  1   5  -8   3  -1  -1  -1]
            [  0  11 -19   5  -2  -3  -3]
            sage: B = copy(A)
            sage: B.right_kernel(basis='LLL')
            Free module of degree 7 and rank 2 over Integer Ring
            User basis matrix:
            [ 2 -1  3  1  0  1  1]
            [-5 -3  2 -5  1 -1 -1]
            sage: C = copy(A)
            sage: C.right_kernel(basis='pivot')
            Traceback (most recent call last):
            ...
            ValueError: pivot basis only available over a field, not over Integer Ring
            sage: D = copy(A)
            sage: D.right_kernel(algorithm='pari')
            Free module of degree 7 and rank 2 over Integer Ring
            Echelon basis matrix:
            [  1   5  -8   3  -1  -1  -1]
            [  0  11 -19   5  -2  -3  -3]
            sage: E = copy(A)
            sage: E.right_kernel(algorithm='padic', basis='LLL')
            Free module of degree 7 and rank 2 over Integer Ring
            User basis matrix:
            [-2  1 -3 -1  0 -1 -1]
            [ 5  3 -2  5 -1  1  1]

        Besides the integers, rings may be as general as principal ideal
        domains.  Results are then free modules.  ::

            sage: R.<y> = QQ[]
            sage: A = matrix(R, [[  1,   y, 1+y^2],
            ...                  [y^3, y^2, 2*y^3]])
            sage: K = A.right_kernel(algorithm='default', basis='echelon'); K
            Free module of degree 3 and rank 1 over Univariate Polynomial Ring in y over Rational Field
            Echelon basis matrix:
            [-1 -y  1]
            sage: A*K.basis_matrix().transpose() == zero_matrix(ZZ, 2, 1)
            True

        It is possible to compute a kernel for a matrix over an integral
        domain which is not a PID, but usually this will fail.  ::

            sage: D.<x> = ZZ[]
            sage: A = matrix(D, 2, 2, [[x^2 - x, -x + 5],
            ...                        [x^2 - 8, -x + 2]])
            sage: A.right_kernel()
            Traceback (most recent call last):
            ...
            ArithmeticError: Ideal Ideal (x^2 - x, x^2 - 8) of Univariate Polynomial Ring in x over Integer Ring not principal

        Matrices over non-commutative rings are not a good idea either.
        These are the "usual" quaternions.  ::

            sage: Q.<i,j,k> = QuaternionAlgebra(-1,-1)
            sage: A = matrix(Q, 2, [i,j,-1,k])
            sage: A.right_kernel()
            Traceback (most recent call last):
            ...
            NotImplementedError: Cannot compute a matrix kernel over Quaternion Algebra (-1, -1) with base ring Rational Field

        Sparse matrices, over the rationals and the integers,
        use the same routines as the dense versions. ::

            sage: A = matrix(ZZ, [[0, -1, 1, 1, 2],
            ...                   [1, -2, 0, 1, 3],
            ...                   [-1, 2, 0, -1, -3]],
            ...              sparse=True)
            sage: A.right_kernel()
            Free module of degree 5 and rank 3 over Integer Ring
            Echelon basis matrix:
            [ 1  0  0  2 -1]
            [ 0  1  0 -1  1]
            [ 0  0  1 -3  1]
            sage: B = A.change_ring(QQ)
            sage: B.is_sparse()
            True
            sage: B.right_kernel()
            Vector space of degree 5 and dimension 3 over Rational Field
            Basis matrix:
            [ 1  0  0  2 -1]
            [ 0  1  0 -1  1]
            [ 0  0  1 -3  1]

        With no columns, the kernel can only have dimension zero.
        With no rows, every possible vector is in the kernel.  ::

            sage: A = matrix(QQ, 2, 0)
            sage: A.right_kernel()
            Vector space of degree 0 and dimension 0 over Rational Field
            Basis matrix:
            []
            sage: A = matrix(QQ, 0, 2)
            sage: A.right_kernel()
            Vector space of degree 2 and dimension 2 over Rational Field
            Basis matrix:
            [1 0]
            [0 1]

        Every vector is in the kernel of a zero matrix, the
        dimension is the number of columns. ::

            sage: A = zero_matrix(QQ, 10, 20)
            sage: A.right_kernel()
            Vector space of degree 20 and dimension 20 over Rational Field
            Basis matrix:
            20 x 20 dense matrix over Rational Field

        Results are cached as the right kernel of the matrix.
        Subsequent requests for the right kernel will return
        the cached result, without regard for new values of the
        algorithm or format keyword.  Work with a copy if you
        need a new right kernel, or perhaps investigate the
        :meth:`right_kernel_matrix` method, which does not
        cache its results and is more flexible. ::

            sage: A = matrix(QQ, 3, range(9))
            sage: K1 = A.right_kernel(basis='echelon')
            sage: K1
            Vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [ 1 -2  1]
            sage: K2 = A.right_kernel(basis='pivot')
            sage: K2
            Vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [ 1 -2  1]
            sage: K1 is K2
            True
            sage: B = copy(A)
            sage: K3 = B.kernel(basis='pivot')
            sage: K3
            Vector space of degree 3 and dimension 1 over Rational Field
            User basis matrix:
            [ 1 -2  1]
            sage: K3 is K1
            False
            sage: K3 == K1
            True
        """
        K = self.fetch('right_kernel')
        if not K is None:
            verbose("retrieving cached right kernel for %sx%s matrix" % (self.nrows(), self.ncols()),level=1)
            return K

        R = self.base_ring()
        tm = verbose("computing a right kernel for %sx%s matrix over %s" % (self.nrows(), self.ncols(), R),level=1)

        # Sanitize basis format
        #   'computed' is OK in right_kernel_matrix(), but not here
        #   'echelon' is default here (and elsewhere)
        #   everything else gets checked in right_kernel_matrix()
        basis = kwds.pop('basis', None)
        if basis == 'computed':
            raise ValueError("kernel basis format 'computed' not supported for kernels (just right kernel matrices)")
        if basis is None:
            basis = 'echelon'
        kwds['basis'] = basis

        # Go get the kernel matrix, this is where it all happens
        M = self.right_kernel_matrix(*args, **kwds)

        ambient = R**self.ncols()
        if basis == 'echelon':
            K = ambient.submodule(M.rows(), already_echelonized=True, check=False)
        else:
            K = ambient.submodule_with_basis(M.rows(), already_echelonized=False, check=False)

        verbose("done computing a right kernel for %sx%s matrix over %s" % (self.nrows(), self.ncols(), R),level=1, t=tm)
        self.cache('right_kernel', K)
        return K

    def left_kernel(self, *args, **kwds):
        r"""
        Returns the left kernel of this matrix, as a vector space or free module.
        This is the set of vectors ``x`` such that ``x*self = 0``.

        .. note::

            For the right kernel, use :meth:`right_kernel`.  The method
            :meth:`kernel` is exactly equal to :meth:`left_kernel`.

        INPUT:

        - ``algorithm`` - default: 'default' - a keyword that selects the
          algorithm employed.  Allowable values are:

          - 'default' - allows the algorithm to be chosen automatically
          - 'generic' - naive algorithm usable for matrices over any field
          - 'pari' - PARI library code for matrices over number fields
            or the integers
          - 'padic' - padic algorithm from IML library for matrices
            over the rationals and integers
          - 'pluq' - PLUQ matrix factorization for matrices mod 2

        - ``basis`` - default: 'echelon' - a keyword that describes
          the format of the basis used to construct the left kernel.
          Allowable values are:

          - 'echelon': the basis matrix is in echelon form
          - 'pivot' : each basis vector is computed from the reduced
            row-echelon form of ``self`` by placing a single one in a
            non-pivot column and zeros in the remaining non-pivot columns.
            Only available for matrices over fields.
          - 'LLL': an LLL-reduced basis.  Only available for matrices
            over the integers.

        OUTPUT:

        A vector space or free module whose degree equals the number
        of rows in ``self`` and contains all the vectors ``x`` such
        that ``x*self = 0``.

        If ``self`` has 0 rows, the kernel has dimension 0, while if ``self``
        has 0 columns the kernel is the entire ambient vector space.

        The result is cached.  Requesting the left kernel a second time,
        but with a different basis format will return the cached result
        with the format from the first computation.

        .. note::

           For much more detailed documentation of the various options see
           :meth:`right_kernel`, since this method just computes
           the right kernel of the transpose of ``self``.

        EXAMPLES:

        Over the rationals with a basis matrix in echelon form. ::

            sage: A = matrix(QQ, [[1, 2, 4, -7, 4],
            ...                   [1, 1, 0, 2, -1],
            ...                   [1, 0, 3, -3, 1],
            ...                   [0, -1, -1, 3, -2],
            ...                   [0, 0, -1, 2, -1]])
            sage: A.left_kernel()
            Vector space of degree 5 and dimension 2 over Rational Field
            Basis matrix:
            [ 1  0 -1  2 -1]
            [ 0  1 -1  1 -4]

        Over a finite field, with a basis matrix in "pivot" format. ::

            sage: A = matrix(FiniteField(7), [[5, 0, 5, 2, 4],
            ...                               [1, 3, 2, 3, 6],
            ...                               [1, 1, 6, 5, 3],
            ...                               [2, 5, 6, 0, 0]])
            sage: A.kernel(basis='pivot')
            Vector space of degree 4 and dimension 2 over Finite Field of size 7
            User basis matrix:
            [5 2 1 0]
            [6 3 0 1]

        The left kernel of a zero matrix is the entire ambient vector
        space whose degree equals the number of rows of ``self``
        (i.e. everything).  ::

            sage: A = MatrixSpace(QQ, 3, 4)(0)
            sage: A.kernel()
            Vector space of degree 3 and dimension 3 over Rational Field
            Basis matrix:
            [1 0 0]
            [0 1 0]
            [0 0 1]

        We test matrices with no rows or columns. ::

            sage: A = matrix(QQ, 2, 0)
            sage: A.left_kernel()
            Vector space of degree 2 and dimension 2 over Rational Field
            Basis matrix:
            [1 0]
            [0 1]
            sage: A = matrix(QQ, 0, 2)
            sage: A.left_kernel()
            Vector space of degree 0 and dimension 0 over Rational Field
            Basis matrix:
            []

        The results are cached. Note that requesting a new format
        for the basis is ignored and the cached copy is returned.
        Work with a copy if you need a new left kernel, or perhaps
        investigate the :meth:`right_kernel_matrix` method on the
        transpose, which does not cache its results and is more
        flexible.  ::

            sage: A = matrix(QQ, [[1,1],[2,2]])
            sage: K1 = A.left_kernel()
            sage: K1
            Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [   1 -1/2]
            sage: K2 = A.left_kernel()
            sage: K1 is K2
            True
            sage: K3 = A.left_kernel(basis='pivot')
            sage: K3
            Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [   1 -1/2]
            sage: B = copy(A)
            sage: K3 = B.left_kernel(basis='pivot')
            sage: K3
            Vector space of degree 2 and dimension 1 over Rational Field
            User basis matrix:
            [-2  1]
            sage: K3 is K1
            False
            sage: K3 == K1
            True
        """
        K = self.fetch('left_kernel')
        if not K is None:
            verbose("retrieving cached left kernel for %sx%s matrix" % (self.nrows(), self.ncols()),level=1)
            return K

        tm = verbose("computing left kernel for %sx%s matrix" % (self.nrows(), self.ncols()),level=1)
        K = self.transpose().right_kernel(*args, **kwds)
        self.cache('left_kernel', K)
        verbose("done computing left kernel for %sx%s matrix" % (self.nrows(), self.ncols()),level=1,t=tm)
        return K

    # .kernel() is a an alias for .left_kernel()
    kernel = left_kernel

    def kernel_on(self, V, poly=None, check=True):
        """
        Return the kernel of self restricted to the invariant subspace V.
        The result is a vector subspace of V, which is also a subspace
        of the ambient space.

        INPUT:

        - ``V`` - vector subspace

        - ``check`` - (optional) default: True; whether to check that
          V is invariant under the action of self.

        - ``poly`` - (optional) default: None; if not None, compute instead
          the kernel of poly(self) on V.

        OUTPUT:

        - a subspace

        .. warning::

           This function does *not* check that V is in fact
           invariant under self if check is False.  With check False this
           function is much faster.

        EXAMPLES::

            sage: t = matrix(QQ, 4, [39, -10, 0, -12, 0, 2, 0, -1, 0, 1, -2, 0, 0, 2, 0, -2]); t
            [ 39 -10   0 -12]
            [  0   2   0  -1]
            [  0   1  -2   0]
            [  0   2   0  -2]
            sage: t.fcp()
            (x - 39) * (x + 2) * (x^2 - 2)
            sage: s = (t-39)*(t^2-2)
            sage: V = s.kernel(); V
            Vector space of degree 4 and dimension 3 over Rational Field
            Basis matrix:
            [1 0 0 0]
            [0 1 0 0]
            [0 0 0 1]
            sage: s.restrict(V)
            [0 0 0]
            [0 0 0]
            [0 0 0]
            sage: s.kernel_on(V)
            Vector space of degree 4 and dimension 3 over Rational Field
            Basis matrix:
            [1 0 0 0]
            [0 1 0 0]
            [0 0 0 1]
            sage: k = t-39
            sage: k.restrict(V)
            [  0 -10 -12]
            [  0 -37  -1]
            [  0   2 -41]
            sage: ker = k.kernel_on(V); ker
            Vector space of degree 4 and dimension 1 over Rational Field
            Basis matrix:
            [   1 -2/7    0 -2/7]
            sage: ker.0 * k
            (0, 0, 0, 0)

        Test that trac ticket #9425 is fixed.

        ::

            sage: V = span([[1/7,0,0] ,[0,1,0]], ZZ); V
            Free module of degree 3 and rank 2 over Integer Ring
            Echelon basis matrix:
            [1/7   0   0]
            [  0   1   0]
            sage: T = matrix(ZZ,3,[1,0,0,0,0,0,0,0,0]); T
            [1 0 0]
            [0 0 0]
            [0 0 0]
            sage: W = T.kernel_on(V); W.basis()
            [
            (0, 1, 0)
            ]
            sage: W.is_submodule(V)
            True
        """
        A = self.restrict(V, check=check)
        if not poly is None:
            A = poly(A)
        W = A.kernel()
        if V.is_ambient():
            return W
        else:
            A = V.basis_matrix()
            B = W.basis_matrix()
            C = B*A
            return C.row_module(base_ring=V.base_ring())


    def integer_kernel(self, ring=ZZ):
        """
        Return the kernel of this matrix over the given ring (which should be
        either the base ring, or a PID whose fraction field is the base ring).

        Assume that the base field of this matrix has a numerator and
        denominator functions for its elements, e.g., it is the rational
        numbers or a fraction field. This function computes a basis over
        the integers for the kernel of self.

        If the matrix is not coercible into QQ, then the PID itself should be
        given as a second argument, as in the third example below.

        EXAMPLES::

            sage: A = MatrixSpace(QQ, 4)(range(16))
            sage: A.integer_kernel()
            Free module of degree 4 and rank 2 over Integer Ring
            Echelon basis matrix:
            [ 1  0 -3  2]
            [ 0  1 -2  1]

        The integer kernel even makes sense for matrices with fractional
        entries::

            sage: A = MatrixSpace(QQ, 2)(['1/2',0,  0, 0])
            sage: A.integer_kernel()
            Free module of degree 2 and rank 1 over Integer Ring
            Echelon basis matrix:
            [0 1]

        An example over a bigger ring::

            sage: L.<w> = NumberField(x^2 - x + 2)
            sage: OL = L.ring_of_integers()
            sage: A = matrix(L, 2, [1, w/2])
            sage: A.integer_kernel(OL)
            Free module of degree 2 and rank 1 over Maximal Order in Number Field in w with defining polynomial x^2 - x + 2
            Echelon basis matrix:
            [    -1 -w + 1]

        """
        try:
            A, _ = self._clear_denom()
            return A.kernel()
        except AttributeError:
            d = self.denominator()
            A = self*d
            M = matrix_space.MatrixSpace(ring, self.nrows(), self.ncols())(A)
            return M.kernel()

    def image(self):
        """
        Return the image of the homomorphism on rows defined by this
        matrix.

        EXAMPLES::

            sage: MS1 = MatrixSpace(ZZ,4)
            sage: MS2 = MatrixSpace(QQ,6)
            sage: A = MS1.matrix([3,4,5,6,7,3,8,10,14,5,6,7,2,2,10,9])
            sage: B = MS2.random_element()

        ::

            sage: image(A)
            Free module of degree 4 and rank 4 over Integer Ring
            Echelon basis matrix:
            [  1   0   0 426]
            [  0   1   0 518]
            [  0   0   1 293]
            [  0   0   0 687]

        ::

            sage: image(B) == B.row_module()
            True
        """
        return self.row_module()

    def _row_ambient_module(self, base_ring=None):
        if base_ring is None:
            base_ring = self.base_ring()
        x = self.fetch('row_ambient_module_%s'%base_ring)
        if not x is None:
            return x
        x = sage.modules.free_module.FreeModule(base_ring, self.ncols(), sparse=self.is_sparse())
        self.cache('row_ambient_module',x)
        return x

    def row_module(self, base_ring=None):
        """
        Return the free module over the base ring spanned by the rows of
        self.

        EXAMPLES::

            sage: A = MatrixSpace(IntegerRing(), 2)([1,2,3,4])
            sage: A.row_module()
            Free module of degree 2 and rank 2 over Integer Ring
            Echelon basis matrix:
            [1 0]
            [0 2]
        """
        M = self._row_ambient_module(base_ring = base_ring)
        if (base_ring is None or base_ring == self.base_ring()) and self.fetch('in_echelon_form'):
            if self.rank() != self.nrows():
                rows = self.matrix_from_rows(range(self.rank())).rows()
            else:
                rows = self.rows()
            return M.span(rows, already_echelonized=True)
        else:
            return M.span(self.rows(), already_echelonized=False)

    def row_space(self, base_ring=None):
        """
        Return the row space of this matrix. (Synonym for
        self.row_module().)

        EXAMPLES::

            sage: t = matrix(QQ, 3, range(9)); t
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: t.row_space()
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [ 1  0 -1]
            [ 0  1  2]

        ::

            sage: m = Matrix(Integers(5),2,2,[2,2,2,2]);
            sage: m.row_space()
            Vector space of degree 2 and dimension 1 over Ring of integers modulo 5
            Basis matrix:
            [1 1]
        """
        return self.row_module(base_ring=base_ring)


    def _column_ambient_module(self):
        x = self.fetch('column_ambient_module')
        if not x is None:
            return x
        x = sage.modules.free_module.FreeModule(self.base_ring(), self.nrows(),
                                                sparse=self.is_sparse())
        self.cache('column_ambient_module',x)
        return x

    def column_module(self):
        """
        Return the free module over the base ring spanned by the columns of
        this matrix.

        EXAMPLES::

            sage: t = matrix(QQ, 3, range(9)); t
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: t.column_module()
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [ 1  0 -1]
            [ 0  1  2]
        """
        return self.transpose().row_module()

    def column_space(self):
        """
        Return the vector space over the base ring spanned by the columns
        of this matrix.

        EXAMPLES::

            sage: M = MatrixSpace(QQ,3,3)
            sage: A = M([1,9,-7,4/5,4,3,6,4,3])
            sage: A.column_space()
            Vector space of degree 3 and dimension 3 over Rational Field
            Basis matrix:
            [1 0 0]
            [0 1 0]
            [0 0 1]
            sage: W = MatrixSpace(CC,2,2)
            sage: B = W([1, 2+3*I,4+5*I,9]); B
            [                     1.00000000000000 2.00000000000000 + 3.00000000000000*I]
            [4.00000000000000 + 5.00000000000000*I                      9.00000000000000]
            sage: B.column_space()
            Vector space of degree 2 and dimension 2 over Complex Field with 53 bits of precision
            Basis matrix:
            [ 1.00000000000000 0.000000000000000]
            [0.000000000000000  1.00000000000000]
        """
        return self.column_module()



    def decomposition(self, algorithm='spin',
                      is_diagonalizable=False, dual=False):
        """
        Returns the decomposition of the free module on which this matrix A
        acts from the right (i.e., the action is x goes to x A), along with
        whether this matrix acts irreducibly on each factor. The factors
        are guaranteed to be sorted in the same way as the corresponding
        factors of the characteristic polynomial.

        Let A be the matrix acting from the on the vector space V of column
        vectors. Assume that A is square. This function computes maximal
        subspaces W_1, ..., W_n corresponding to Galois conjugacy classes
        of eigenvalues of A. More precisely, let `f(X)` be the characteristic
        polynomial of A. This function computes the subspace
        `W_i = ker(g_(A)^n)`, where `g_i(X)` is an irreducible
        factor of `f(X)` and `g_i(X)` exactly divides `f(X)`. If the optional
        parameter is_diagonalizable is True, then we let `W_i = ker(g(A))`,
        since then we know that `ker(g(A)) = ker(g(A)^n)`.

        INPUT:


        -  ``self`` - a matrix

        -  ``algorithm`` - 'spin' (default): algorithm involves
           iterating the action of self on a vector. 'kernel': naively just
           compute `ker(f_i(A))` for each factor `f_i`.

        -  ``dual`` - bool (default: False): If True, also
           returns the corresponding decomposition of V under the action of
           the transpose of A. The factors are guaranteed to correspond.

        -  ``is_diagonalizable`` - if the matrix is known to
           be diagonalizable, set this to True, which might speed up the
           algorithm in some cases.

        .. note::

           If the base ring is not a field, the kernel algorithm is
           used.


        OUTPUT:


        - ``Sequence`` - list of pairs (V,t), where V is a vector
          spaces and t is a bool, and t is True exactly when the
          charpoly of self on V is irreducible.


        - (optional) list - list of pairs (W,t), where W is a vector
          space and t is a bool, and t is True exactly when the
          charpoly of the transpose of self on W is irreducible.

        EXAMPLES::

            sage: A = matrix(ZZ, 4, [3,4,5,6,7,3,8,10,14,5,6,7,2,2,10,9])
            sage: B = matrix(QQ, 6, range(36))
            sage: B*11
            [  0  11  22  33  44  55]
            [ 66  77  88  99 110 121]
            [132 143 154 165 176 187]
            [198 209 220 231 242 253]
            [264 275 286 297 308 319]
            [330 341 352 363 374 385]
            sage: A.decomposition()
            [
            (Ambient free module of rank 4 over the principal ideal domain Integer Ring, True)
            ]
            sage: B.decomposition()
            [
            (Vector space of degree 6 and dimension 2 over Rational Field
            Basis matrix:
            [ 1  0 -1 -2 -3 -4]
            [ 0  1  2  3  4  5], True),
            (Vector space of degree 6 and dimension 4 over Rational Field
            Basis matrix:
            [ 1  0  0  0 -5  4]
            [ 0  1  0  0 -4  3]
            [ 0  0  1  0 -3  2]
            [ 0  0  0  1 -2  1], False)
            ]
        """
        if algorithm == 'kernel' or not self.base_ring().is_field():
            return self._decomposition_using_kernels(is_diagonalizable = is_diagonalizable, dual=dual)
        elif algorithm == 'spin':
            X = self._decomposition_spin_generic(is_diagonalizable = is_diagonalizable)
            if dual:
                Y = self.transpose()._decomposition_spin_generic(is_diagonalizable = is_diagonalizable)
                return X, Y
            return X
        else:
            raise ValueError, "no algorithm '%s'"%algorithm

    def _decomposition_spin_generic(self, is_diagonalizable=False):
        r"""
        Compute the decomposition of this matrix using the spin algorithm.

        INPUT:

        - ``self`` - a matrix with field entries

        OUTPUT: a list of reduced row echelon form basis

        AUTHORS:

        - William Stein
        """
        if not self.is_square():
            raise ValueError, "self must be a square matrix"

        if not self.base_ring().is_field():
            raise TypeError, "self must be over a field."

        if self.nrows() == 0:
            return decomp_seq([])

        f = self.charpoly('x')
        E = decomp_seq([])

        t = verbose('factoring the characteristic polynomial', level=2, caller_name='generic spin decomp')
        F = f.factor()
        verbose('done factoring', t=t, level=2, caller_name='generic spin decomp')

        if len(F) == 1:
            V = self.base_ring()**self.nrows()
            return decomp_seq([(V,F[0][1]==1)])

        V = self.base_ring()**self.nrows()
        v = V.random_element()
        num_iterates = max([0] + [f.degree() - g.degree() for g, _ in F if g.degree() > 1]) + 1

        S = [ ]

        F.sort()
        for i in range(len(F)):
            g, m = F[i]

            if g.degree() == 1:
                # Just use kernel -- much easier.
                B = self.__copy__()
                for k from 0 <= k < self.nrows():
                    B[k,k] += g[0]
                if m > 1 and not is_diagonalizable:
                    B = B**m
                W = B.kernel()
                E.append((W, m==1))
                continue

            # General case, i.e., deg(g) > 1:
            W = None
            tries = m
            while True:

                # Compute the complementary factor.
                h = f // (g**m)
                v = h.list()

                while len(S) < tries:
                    t = verbose('%s-spinning %s-th random vector'%(num_iterates, len(S)), level=2, caller_name='generic spin decomp')
                    S.append(self.iterates(V.random_element(), num_iterates))
                    verbose('done spinning', level=2, t=t, caller_name='generic spin decomp')

                for j in range(0 if W is None else W.nrows() // g.degree(), len(S)):
                    # Compute one element of the kernel of g(A)**m.
                    t = verbose('compute element of kernel of g(A), for g of degree %s'%g.degree(),level=2,
                                caller_name='generic spin decomp')
                    w = S[j].linear_combination_of_rows(h.list())
                    t = verbose('done computing element of kernel of g(A)', t=t,level=2, caller_name='generic spin decomp')

                    # Get the rest of the kernel.
                    t = verbose('fill out rest of kernel',level=2, caller_name='generic spin decomp')
                    if W is None:
                        W = self.iterates(w, g.degree())
                    else:
                        W = W.stack(self.iterates(w, g.degree()))
                    t = verbose('finished filling out more of kernel',level=2, t=t, caller_name='generic spin decomp')

                if W.rank() == m * g.degree():
                    t = verbose('now computing row space', level=2, caller_name='generic spin decomp')
                    W.echelonize()
                    E.append((W.row_space(), m==1))
                    verbose('computed row space', level=2,t=t, caller_name='generic spin decomp')
                    break
                else:
                    verbose('we have not yet generated all the kernel (rank so far=%s, target rank=%s)'%(
                        W.rank(), m*g.degree()), level=2, caller_name='generic spin decomp')
                    tries += 1
                    if tries > 1000*m:  # avoid an insanely long infinite loop
                        raise RuntimeError, "likely bug in decomposition"
                # end if
            #end while
        #end for
        return E

    def _decomposition_using_kernels(self, is_diagonalizable=False, dual=False):
        if not self.is_square():
            raise ValueError, "self must be a square matrix"

        if self.nrows() == 0:
            return decomp_seq([])

        f = self.charpoly('x')
        E = decomp_seq([])

        # Idea: For optimization, could compute powers of self
        #       up to max degree of any factor.  Then get g(self)
        #       by taking a linear combination.

        if dual:
            Edual = decomp_seq([])
        F = f.factor()
        if len(F) == 1:
            V = sage.modules.free_module.FreeModule(
                              self.base_ring(), self.nrows(), sparse=self.is_sparse())
            m = F[0][1]
            if dual:
                return decomp_seq([(V, m==1)]), decomp_seq([(V, m==1)])
            else:
                return decomp_seq([(V, m==1)])
        F.sort()
        for g, m in f.factor():
            t = verbose('decomposition -- Computing g(self) for an irreducible factor g of degree %s'%g.degree(),level=2)
            if is_diagonalizable:
                B = g(self)
            else:
                B = g(self)
                t2 = verbose('decomposition -- raising g(self) to the power %s'%m,level=2)
                B = B ** m
                verbose('done powering',t2)
            t = verbose('decomposition -- done computing g(self)', level=2, t=t)
            E.append((B.kernel(), m==1))
            t = verbose('decomposition -- time to compute kernel', level=2, t=t)
            if dual:
                Edual.append((B.transpose().kernel(), m==1))
                verbose('decomposition -- time to compute dual kernel', level=2, t=t)
        if dual:
            return E, Edual
        return E

    def decomposition_of_subspace(self, M, check_restrict = True, **kwds):
        """
        Suppose the right action of self on M leaves M invariant. Return
        the decomposition of M as a list of pairs (W, is_irred) where
        is_irred is True if the charpoly of self acting on the factor W is
        irreducible.

        Additional inputs besides M are passed onto the decomposition
        command.

        INPUT:

            - `M`                -- A subspace of the free module ``self`` acts on.
            - ``check_restrict`` -- A boolean (default: ``True``); Call restrict
                                    with or without check.
            - ``kwds``           -- Keywords that will be forwarded to :meth:`~.decomposition`.

        EXAMPLES::

            sage: t = matrix(QQ, 3, [3, 0, -2, 0, -2, 0, 0, 0, 0]); t
            [ 3  0 -2]
            [ 0 -2  0]
            [ 0  0  0]
            sage: t.fcp('X')   # factored charpoly
            (X - 3) * X * (X + 2)
            sage: v = kernel(t*(t+2)); v   # an invariant subspace
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [0 1 0]
            [0 0 1]
            sage: D = t.decomposition_of_subspace(v); D
            [
            (Vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [0 0 1], True),
            (Vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [0 1 0], True)
            ]
            sage: t.restrict(D[0][0])
            [0]
            sage: t.restrict(D[1][0])
            [-2]

        We do a decomposition over ZZ::

            sage: a = matrix(ZZ,6,[0, 0, -2, 0, 2, 0, 2, -4, -2, 0, 2, 0, 0, 0, -2, -2, 0, 0, 2, 0, -2, -4, 2, -2, 0, 2, 0, -2, -2, 0, 0, 2, 0, -2, 0, 0])
            sage: a.decomposition_of_subspace(ZZ^6)
            [
            (Free module of degree 6 and rank 2 over Integer Ring
            Echelon basis matrix:
            [ 1  0  1 -1  1 -1]
            [ 0  1  0 -1  2 -1], False),
            (Free module of degree 6 and rank 4 over Integer Ring
            Echelon basis matrix:
            [ 1  0 -1  0  1  0]
            [ 0  1  0  0  0  0]
            [ 0  0  0  1  0  0]
            [ 0  0  0  0  0  1], False)
            ]

        TESTS::

            sage: t = matrix(QQ, 3, [3, 0, -2, 0, -2, 0, 0, 0, 0]);
            sage: t.decomposition_of_subspace(v, check_restrict = False) == t.decomposition_of_subspace(v)
            True
        """
        if not sage.modules.free_module.is_FreeModule(M):
            raise TypeError, "M must be a free module."
        if not self.is_square():
            raise ArithmeticError, "self must be a square matrix"
        if M.base_ring() != self.base_ring():
            raise ArithmeticError, "base rings must be the same, but self is over %s and module is over %s"%(
                self.base_ring(), M.base_ring())
        if M.degree() != self.ncols():
            raise ArithmeticError, \
               "M must be a subspace of an %s-dimensional space"%self.ncols()

        time = verbose(t=0)

        # 1. Restrict
        B = self.restrict(M, check = check_restrict)
        time0 = verbose("decompose restriction -- ", time)

        # 2. Decompose restriction
        D = B.decomposition(**kwds)

        sum_dim = sum([A.dimension() for A,_ in D])
        assert sum_dim == M.dimension(), \
               "bug in decomposition; " + \
               "the sum of the dimensions (=%s) of the factors must equal the dimension (%s) of the acted on space:\nFactors found: %s\nSpace: %s"%(sum_dim, M.dimension(), D, M)

        # 3. Lift decomposition to subspaces of ambient vector space.
        # Each basis vector for an element of D defines a linear
        # combination of the basis of W, and these linear combinations
        # define the corresponding subspaces of the ambient space M.

        verbose("decomposition -- ", time0)
        C = M.basis_matrix()

        D = [((W.basis_matrix() * C).row_module(self.base_ring()), is_irred) for W, is_irred in D]
        D = decomp_seq(D)

        verbose(t=time)
        return D

    def restrict(self, V, check=True):
        """
        Returns the matrix that defines the action of self on the chosen
        basis for the invariant subspace V. If V is an ambient, returns
        self (not a copy of self).

        INPUT:


        -  ``V`` - vector subspace

        -  ``check`` - (optional) default: True; if False may
           not check that V is invariant (hence can be faster).


        OUTPUT: a matrix

        .. warning::

           This function returns an nxn matrix, where V has dimension
           n. It does *not* check that V is in fact invariant under
           self, unless check is True.

        EXAMPLES::

            sage: V = VectorSpace(QQ, 3)
            sage: M = MatrixSpace(QQ, 3)
            sage: A = M([1,2,0, 3,4,0, 0,0,0])
            sage: W = V.subspace([[1,0,0], [0,1,0]])
            sage: A.restrict(W)
            [1 2]
            [3 4]
            sage: A.restrict(W, check=True)
            [1 2]
            [3 4]

        We illustrate the warning about invariance not being checked by
        default, by giving a non-invariant subspace. With the default
        check=False this function returns the 'restriction' matrix, which
        is meaningless as check=True reveals.

        ::

            sage: W2 = V.subspace([[1,0,0], [0,1,1]])
            sage: A.restrict(W2, check=False)
            [1 2]
            [3 4]
            sage: A.restrict(W2, check=True)
            Traceback (most recent call last):
            ...
            ArithmeticError: subspace is not invariant under matrix
        """
        if not isinstance(V, sage.modules.free_module.FreeModule_generic):
            raise TypeError, "V must be a free module"
        #if V.base_ring() != self.base_ring():
        #     raise ValueError, "matrix and module must have the same base ring, but matrix is over %s and module is over %s"%(self.base_ring(), V.base_ring())
        if V.degree() != self.nrows():
            raise IndexError, "degree of V (=%s) must equal number of rows of self (=%s)"%(\
                V.degree(), self.nrows())
        if V.rank() == 0 or V.degree() == 0:
            return self.new_matrix(nrows=0, ncols=0)

        if not check and V.base_ring().is_field() and not V.has_user_basis():
            B = V.echelonized_basis_matrix()
            P = B.pivots()
            return B*self.matrix_from_columns(P)
        else:
            n = V.rank()
            try:
                # todo optimize so only involves matrix multiplies ?
                C = [V.coordinate_vector(b*self) for b in V.basis()]
            except ArithmeticError:
                raise ArithmeticError, "subspace is not invariant under matrix"
            return self.new_matrix(n, n, C, sparse=False)

    def restrict_domain(self, V):
        """
        Compute the matrix relative to the basis for V on the domain
        obtained by restricting self to V, but not changing the codomain of
        the matrix. This is the matrix whose rows are the images of the
        basis for V.

        INPUT:


        -  ``V`` - vector space (subspace of ambient space on
           which self acts)


        .. seealso::

           :meth:`restrict`

        EXAMPLES::

            sage: V = QQ^3
            sage: A = matrix(QQ,3,[1,2,0, 3,4,0, 0,0,0])
            sage: W = V.subspace([[1,0,0], [1,2,3]])
            sage: A.restrict_domain(W)
            [1 2 0]
            [3 4 0]
            sage: W2 = V.subspace_with_basis([[1,0,0], [1,2,3]])
            sage: A.restrict_domain(W2)
            [ 1  2  0]
            [ 7 10  0]
        """
        return V.basis_matrix() * self

    def restrict_codomain(self, V):
        r"""
        Suppose that self defines a linear map from some domain to a
        codomain that contains `V` and that the image of self is
        contained in `V`. This function returns a new matrix
        `A` that represents this linear map but as a map to
        `V`, in the sense that if `x` is in the domain,
        then `xA` is the linear combination of the elements of the
        basis of `V` that equals v\*self.

        INPUT:


        -  ``V`` - vector space (space of degree
           ``self.ncols()``) that contains the image of self.


        .. seealso::

           :meth:`restrict`, :meth:`restrict_domain`

        EXAMPLES::

            sage: A = matrix(QQ,3,[1..9])
            sage: V = (QQ^3).span([[1,2,3], [7,8,9]]); V
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [ 1  0 -1]
            [ 0  1  2]
            sage: z = vector(QQ,[1,2,5])
            sage: B = A.restrict_codomain(V); B
            [1 2]
            [4 5]
            [7 8]
            sage: z*B
            (44, 52)
            sage: z*A
            (44, 52, 60)
            sage: 44*V.0 + 52*V.1
            (44, 52, 60)
        """
        return V.basis_matrix().solve_left(self)

    def maxspin(self, v):
        """
        Computes the largest integer n such that the list of vectors
        `S=[v, v*A, ..., v * A^n]` are linearly independent, and
        returns that list.

        INPUT:


        -  ``self`` - Matrix

        -  ``v`` - Vector


        OUTPUT:


        -  ``list`` - list of Vectors


        ALGORITHM: The current implementation just adds vectors to a vector
        space until the dimension doesn't grow. This could be optimized by
        directly using matrices and doing an efficient Echelon form. Also,
        when the base is Q, maybe we could simultaneously keep track of
        what is going on in the reduction modulo p, which might make things
        much faster.

        EXAMPLES::

            sage: t = matrix(QQ, 3, range(9)); t
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: v = (QQ^3).0
            sage: t.maxspin(v)
            [(1, 0, 0), (0, 1, 2), (15, 18, 21)]
            sage: k = t.kernel(); k
            Vector space of degree 3 and dimension 1 over Rational Field
            Basis matrix:
            [ 1 -2  1]
            sage: t.maxspin(k.0)
            [(1, -2, 1)]
        """
        if v == 0:
            return []
        if not is_FreeModuleElement(v):
            raise TypeError, "v must be a FreeModuleElement"
        VS = v.parent()
        V = VS.span([v])
        w = v
        S = [v]
        while True:
            w = w*self
            W = V + VS.span([w])
            if W.dimension() == V.dimension():
                return S
            V = W
            S.append(w)


    def wiedemann(self, i, t=0):
        """
        Application of Wiedemann's algorithm to the i-th standard basis
        vector.

        INPUT:


        -  ``i`` - an integer

        -  ``t`` - an integer (default: 0) if t is nonzero, use
           only the first t linear recurrence relations.


        IMPLEMENTATION: This is a toy implementation.

        EXAMPLES::

            sage: t = matrix(QQ, 3, range(9)); t
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: t.wiedemann(0)
            x^2 - 12*x - 18
            sage: t.charpoly()
            x^3 - 12*x^2 - 18*x
        """
        i = int(i); t=int(t)
        if self.nrows() != self.ncols():
            raise ArithmeticError, "self must be a square matrix"
        n = self.nrows()
        v = sage.modules.free_module.VectorSpace(self.base_ring(), n).gen(i)
        tm = verbose('computing iterates...')
        cols = self.iterates(v, 2*n).columns()
        tm = verbose('computed iterates', tm)
        f = None
        # Compute the minimal polynomial of the linear recurrence
        # sequence corresponding to the 0-th entries of the iterates,
        # then the 1-th entries, etc.
        if t == 0:
            R = range(n)
        else:
            R = [t]
        for i in R:
            tm = verbose('applying berlekamp-massey')
            g = berlekamp_massey.berlekamp_massey(cols[i].list())
            verbose('berlekamp-massey done', tm)
            if f is None:
                f = g
            else:
                f = f.lcm(g)
            if f.degree() == n:
                break
        return f

    # Deprecated Aug 2008 Trac #3794
    # Removed July 2011
    # def eigenspaces(self, var='a', even_if_inexact=None):

    def _eigenspace_format(self, format):
        r"""
        Helper method to control output format for eigenspaces.

        INPUT:

        - ``format`` - ``None``, ``'all'`` or ``'galois'``

        OUTPUT:

        Any format except ``None`` is just passed through.  When the
        format is ``None`` a choice is made about the style of the output.
        If there is an algebraically closed field that will contain the
        possible eigenvalues, then 'all" of the eigenspaces are given.

        However if this is not the case, then only one eigenspace is output
        for each irreducible factor of the characteristic polynomial.

        EXAMPLES:

        Pass-through first.  ::

            sage: A = matrix(QQ, 2, range(4))
            sage: A._eigenspace_format('all') == 'all'
            True
            sage: A._eigenspace_format('galois') == 'galois'
            True

        The algebraic closure of the rationals, the field of algebraic numbers,
        (aka ``QQbar``) is implemented, while in general the algebraic closure
        of a finite field is not implemented.  ::

            sage: A = matrix(QQ, 2, range(4))
            sage: A._eigenspace_format(None) == 'all'
            True
            sage: B = matrix(GF(13), 2, range(4))
            sage: B._eigenspace_format(None) == 'galois'
            True

        Subrings are promoted to fraction fields and then checked for the
        existence of algebraic closures.  ::

            sage: A = matrix(ZZ, 2, range(4))
            sage: A._eigenspace_format(None) == 'all'
            True
        """
        if not format in [None, 'all', 'galois']:
            msg = "format keyword must be None, 'all' or 'galois', not {0}"
            raise ValueError(msg.format(format))

        # For subrings of implemented algebraically closed fields we
        #   default to all eigenspaces in the absence of a format keyword
        # Add newly implemented algebraically closed fields to list below
        #   and implement the determintion of the actual eigenvalues
        #   in the eigenspace_left() routine
        if format is None:
            R = self.base_ring()
            from sage.rings.qqbar import QQbar
            try:
                if R.fraction_field().algebraic_closure() in [QQbar]:
                    return 'all'
                else:
                    return 'galois'
            except (NotImplementedError, AttributeError):
                return 'galois'
        else:
            return format

    def eigenspaces_left(self, format='all', var='a', algebraic_multiplicity=False):
        r"""
        Compute the left eigenspaces of a matrix.

        Note that ``eigenspaces_left()`` and ``left_eigenspaces()``
        are identical methods.  Here "left" refers to the eigenvectors
        being placed to the left of the matrix.

        INPUT:

        - ``self`` - a square matrix over an exact field.  For inexact
          matrices consult the numerical or symbolic matrix classes.

        - ``format`` - default: ``None``

          - ``'all'`` - attempts to create every eigenspace.  This will
            always be possible for matrices with rational entries.
          - ``'galois'`` - for each irreducible factor of the characteristic
            polynomial, a single eigenspace will be output for a
            single root/eigenvalue for the irreducible factor.
          - ``None`` - Uses the 'all' format if the base ring is contained
            in an algebraically closed field which is implemented.
            Otherwise, uses the 'galois' format.

        - ``var`` - default: 'a' - variable name used to
          represent elements of the root field of each
          irreducible factor of the characteristic polynomial.
          If var='a', then the root fields will be in terms of
          a0, a1, a2, ...., where the numbering runs across all
          the irreducible factors of the characteristic polynomial,
          even for linear factors.

        - ``algebraic_multiplicity`` - default: False - whether or
          not to include the algebraic multiplicity of each eigenvalue
          in the output.  See the discussion below.

        OUTPUT:

        If algebraic_multiplicity=False, return a list of pairs (e, V)
        where e is an eigenvalue of the matrix, and V is the corresponding
        left eigenspace.  For Galois conjugates of eigenvalues, there
        may be just one representative eigenspace, depending on the
        ``format`` keyword.

        If algebraic_multiplicity=True, return a list of triples (e, V, n)
        where e and V are as above and n is the algebraic multiplicity of
        the eigenvalue.

        .. warning::

           Uses a somewhat naive algorithm (simply factors the
           characteristic polynomial and computes kernels directly
           over the extension field).

        EXAMPLES:

        We compute the left eigenspaces of a `3\times 3`
        rational matrix. First, we request `all` of the eigenvalues,
        so the results are in the field of algebraic numbers, `QQbar`.
        Then we request just one eigenspace per irreducible factor of
        the characteristic polynomial with the `galois` keyword.  ::

            sage: A = matrix(QQ,3,3,range(9)); A
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: es = A.eigenspaces_left(format='all'); es
            [
            (0, Vector space of degree 3 and dimension 1 over Rational Field
            User basis matrix:
            [ 1 -2  1]),
            (-1.348469228349535?, Vector space of degree 3 and dimension 1 over Algebraic Field
            User basis matrix:
            [                   1  0.3101020514433644? -0.3797958971132713?]),
            (13.34846922834954?, Vector space of degree 3 and dimension 1 over Algebraic Field
            User basis matrix:
            [                 1 1.289897948556636? 1.579795897113272?])
            ]

            sage: es = A.eigenspaces_left(format='galois'); es
            [
            (0, Vector space of degree 3 and dimension 1 over Rational Field
            User basis matrix:
            [ 1 -2  1]),
            (a1, Vector space of degree 3 and dimension 1 over Number Field in a1 with defining polynomial x^2 - 12*x - 18
            User basis matrix:
            [            1 1/15*a1 + 2/5 2/15*a1 - 1/5])
            ]
            sage: es = A.eigenspaces_left(format='galois', algebraic_multiplicity=True); es
            [
            (0, Vector space of degree 3 and dimension 1 over Rational Field
            User basis matrix:
            [ 1 -2  1], 1),
            (a1, Vector space of degree 3 and dimension 1 over Number Field in a1 with defining polynomial x^2 - 12*x - 18
            User basis matrix:
            [            1 1/15*a1 + 2/5 2/15*a1 - 1/5], 1)
            ]
            sage: e, v, n = es[0]; v = v.basis()[0]
            sage: delta = e*v - v*A
            sage: abs(abs(delta)) < 1e-10
            True

        The same computation, but with implicit base change to a field.  ::

            sage: A = matrix(ZZ,3,range(9)); A
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: A.eigenspaces_left(format='galois')
            [
            (0, Vector space of degree 3 and dimension 1 over Rational Field
            User basis matrix:
            [ 1 -2  1]),
            (a1, Vector space of degree 3 and dimension 1 over Number Field in a1 with defining polynomial x^2 - 12*x - 18
            User basis matrix:
            [            1 1/15*a1 + 2/5 2/15*a1 - 1/5])
            ]

        We compute the left eigenspaces of the matrix of the Hecke operator
        `T_2` on level 43 modular symbols, both with all eigenvalues (the default)
        and with one subspace per factor. ::

            sage: A = ModularSymbols(43).T(2).matrix(); A
            [ 3  0  0  0  0  0 -1]
            [ 0 -2  1  0  0  0  0]
            [ 0 -1  1  1  0 -1  0]
            [ 0 -1  0 -1  2 -1  1]
            [ 0 -1  0  1  1 -1  1]
            [ 0  0 -2  0  2 -2  1]
            [ 0  0 -1  0  1  0 -1]
            sage: A.base_ring()
            Rational Field
            sage: f = A.charpoly(); f
            x^7 + x^6 - 12*x^5 - 16*x^4 + 36*x^3 + 52*x^2 - 32*x - 48
            sage: factor(f)
            (x - 3) * (x + 2)^2 * (x^2 - 2)^2
            sage: A.eigenspaces_left(algebraic_multiplicity=True)
            [
            (3, Vector space of degree 7 and dimension 1 over Rational Field
            User basis matrix:
            [   1    0  1/7    0 -1/7    0 -2/7], 1),
            (-2, Vector space of degree 7 and dimension 2 over Rational Field
            User basis matrix:
            [ 0  1  0  1 -1  1 -1]
            [ 0  0  1  0 -1  2 -1], 2),
            (-1.414213562373095?, Vector space of degree 7 and dimension 2 over Algebraic Field
            User basis matrix:
            [                  0                   1                   0                  -1 0.4142135623730951?                   1                  -1]
            [                  0                   0                   1                   0                  -1                   0  2.414213562373095?], 2),
            (1.414213562373095?, Vector space of degree 7 and dimension 2 over Algebraic Field
            User basis matrix:
            [                   0                    1                    0                   -1  -2.414213562373095?                    1                   -1]
            [                   0                    0                    1                    0                   -1                    0 -0.4142135623730951?], 2)
            ]
            sage: A.eigenspaces_left(format='galois', algebraic_multiplicity=True)
            [
            (3, Vector space of degree 7 and dimension 1 over Rational Field
            User basis matrix:
            [   1    0  1/7    0 -1/7    0 -2/7], 1),
            (-2, Vector space of degree 7 and dimension 2 over Rational Field
            User basis matrix:
            [ 0  1  0  1 -1  1 -1]
            [ 0  0  1  0 -1  2 -1], 2),
            (a2, Vector space of degree 7 and dimension 2 over Number Field in a2 with defining polynomial x^2 - 2
            User basis matrix:
            [      0       1       0      -1 -a2 - 1       1      -1]
            [      0       0       1       0      -1       0 -a2 + 1], 2)
            ]

        Next we compute the left eigenspaces over the finite field of order 11. ::

            sage: A = ModularSymbols(43, base_ring=GF(11), sign=1).T(2).matrix(); A
            [ 3  9  0  0]
            [ 0  9  0  1]
            [ 0 10  9  2]
            [ 0  9  0  2]
            sage: A.base_ring()
            Finite Field of size 11
            sage: A.charpoly()
            x^4 + 10*x^3 + 3*x^2 + 2*x + 1
            sage: A.eigenspaces_left(format='galois', var = 'beta')
            [
            (9, Vector space of degree 4 and dimension 1 over Finite Field of size 11
            User basis matrix:
            [0 0 1 5]),
            (3, Vector space of degree 4 and dimension 1 over Finite Field of size 11
            User basis matrix:
            [1 6 0 6]),
            (beta2, Vector space of degree 4 and dimension 1 over Univariate Quotient Polynomial Ring in beta2 over Finite Field of size 11 with modulus x^2 + 9
            User basis matrix:
            [           0            1            0 5*beta2 + 10])
            ]

        This method is only applicable to exact matrices.
        The "eigenmatrix" routines for matrices with double-precision
        floating-point entries (``RDF``, ``CDF``) are the best
        alternative.  (Since some platforms return eigenvectors
        that are the negatives of those given here, this one example
        is not tested here.)  There are also "eigenmatrix" routines for
        matrices with symbolic entries.  ::

            sage: A = matrix(QQ, 3, 3, range(9))
            sage: A.change_ring(RR).eigenspaces_left()
            Traceback (most recent call last):
            ...
            NotImplementedError: eigenspaces cannot be computed reliably for inexact rings such as Real Field with 53 bits of precision,
            consult numerical or symbolic matrix classes for other options

            sage: em = A.change_ring(RDF).eigenmatrix_left()
            sage: eigenvalues = em[0]; eigenvalues.dense_matrix().zero_at(2e-15)
            [ 13.3484692283            0.0            0.0]
            [           0.0 -1.34846922835            0.0]
            [           0.0            0.0            0.0]
            sage: eigenvectors = em[1]; eigenvectors # not tested
            [ 0.440242867...  0.567868371...  0.695493875...]
            [ 0.897878732...  0.278434036... -0.341010658...]
            [ 0.408248290... -0.816496580...  0.408248290...]

            sage: x, y = var('x y')
            sage: S = matrix([[x, y], [y, 3*x^2]])
            sage: em = S.eigenmatrix_left()
            sage: eigenvalues = em[0]; eigenvalues
            [3/2*x^2 + 1/2*x - 1/2*sqrt(9*x^4 - 6*x^3 + x^2 + 4*y^2)                                                       0]
            [                                                      0 3/2*x^2 + 1/2*x + 1/2*sqrt(9*x^4 - 6*x^3 + x^2 + 4*y^2)]
            sage: eigenvectors = em[1]; eigenvectors
            [                                                    1 1/2*(3*x^2 - x - sqrt(9*x^4 - 6*x^3 + x^2 + 4*y^2))/y]
            [                                                    1 1/2*(3*x^2 - x + sqrt(9*x^4 - 6*x^3 + x^2 + 4*y^2))/y]


        A request for ``'all'`` the eigenvalues, when it is not
        possible, will raise an error.  Using the ``'galois'``
        format option is more likely to be successful.  ::

            sage: F.<b> = FiniteField(11^2)
            sage: A = matrix(F, [[b + 1, b + 1], [10*b + 4, 5*b + 4]])
            sage: A.eigenspaces_left(format='all')
            Traceback (most recent call last):
            ...
            NotImplementedError: unable to construct eigenspaces for eigenvalues outside the base field,
            try the keyword option: format='galois'

            sage: A.eigenspaces_left(format='galois')
            [
            (a0, Vector space of degree 2 and dimension 1 over Univariate Quotient Polynomial Ring in a0 over Finite Field in b of size 11^2 with modulus x^2 + (5*b + 6)*x + 8*b + 10
            User basis matrix:
            [               1 6*b*a0 + 3*b + 1])
            ]

        TESTS:

        We make sure that :trac:`13308` is fixed. ::

            sage: M = ModularSymbols(Gamma1(23), sign=1)
            sage: m = M.cuspidal_subspace().hecke_matrix(2)
            sage: [j*m==i[0]*j for i in m.eigenspaces_left(format='all') for j in i[1].basis()] # long time (4s)
            [True, True, True, True, True, True, True, True, True, True, True, True]

            sage: B = matrix(QQ, 2, 3, range(6))
            sage: B.eigenspaces_left()
            Traceback (most recent call last):
            ...
            TypeError: matrix must be square, not 2 x 3

            sage: B = matrix(QQ, 4, 4, range(16))
            sage: B.eigenspaces_left(format='junk')
            Traceback (most recent call last):
            ...
            ValueError: format keyword must be None, 'all' or 'galois', not junk

            sage: B.eigenspaces_left(algebraic_multiplicity='garbage')
            Traceback (most recent call last):
            ...
            ValueError: algebraic_multiplicity keyword must be True or False
        """
        if not algebraic_multiplicity in [True, False]:
            msg = 'algebraic_multiplicity keyword must be True or False'
            raise ValueError(msg.format(algebraic_multiplicity))
        if not self.is_square():
            msg = 'matrix must be square, not {0} x {1}'
            raise TypeError(msg.format(self.nrows(), self.ncols()))
        if not self.base_ring().is_exact():
            msg = ("eigenspaces cannot be computed reliably for inexact rings such as {0},\n",
                   "consult numerical or symbolic matrix classes for other options")
            raise NotImplementedError(''.join(msg).format(self.base_ring()))

        format = self._eigenspace_format(format)

        key = 'eigenspaces_left_' + format + '{0}'.format(var)
        x = self.fetch(key)
        if not x is None:
            if algebraic_multiplicity:
                return x
            else:
                return Sequence([(e[0],e[1]) for e in x], cr=True, check=False)

        # Possible improvements:
        # algorithm for dual_eigenvector in sage/modular/hecke/module.py
        # use minpoly when algebraic_multiplicity=False
        # as of 2007-03-25 minpoly is unreliable via linbox

        import sage.categories.homset
        import sage.rings.qqbar

        G = self.fcp()   # factored characteristic polynomial
        V = []
        i = -1 # variable name index, increments for each eigenvalue
        for h, e in G:
            i = i + 1
            if h.degree() == 1:
                alpha = -h[0]/h[1]
                F = alpha.parent()
                if F != self.base_ring():
                    self = self.change_ring(F)
                A = self - alpha
                W = A.kernel()
                V.append((alpha, W.ambient_module().span_of_basis(W.basis()), e))
            else:
                F = h.root_field('{0}{1}'.format(var,i))
                alpha = F.gen(0)
                A = self.change_ring(F) - alpha
                W = A.kernel()
                WB = W.basis()
                if format == 'galois':
                    V.append((alpha, W.ambient_module().span_of_basis(WB), e))
                elif format == 'all':
                    try:
                        alpha_conj = alpha.galois_conjugates(sage.rings.qqbar.QQbar)
                    except AttributeError:
                        msg = ("unable to construct eigenspaces for eigenvalues outside the base field,\n"
                               "try the keyword option: format='galois'")
                        raise NotImplementedError(''.join(msg))
                    for ev in alpha_conj:
                        m = sage.categories.homset.hom(alpha.parent(), ev.parent(), ev)
                        space = (ev.parent())**self.nrows()
                        evec_list = [(space)([m(v_j) for v_j in v]) for v in WB]
                        V.append((ev, space.span_of_basis(evec_list, already_echelonized=True), e))
        V = Sequence(V, cr=True, check=False)
        self.cache(key, V)
        if algebraic_multiplicity:
            return V
        else:
            return Sequence([(e[0],e[1]) for e in V], cr=True, check=False)

    left_eigenspaces = eigenspaces_left

    def eigenspaces_right(self, format='all', var='a', algebraic_multiplicity=False):
        r"""
        Compute the right eigenspaces of a matrix.

        Note that ``eigenspaces_right()`` and ``right_eigenspaces()``
        are identical methods.  Here "right" refers to the eigenvectors
        being placed to the right of the matrix.

        INPUT:

        - ``self`` - a square matrix over an exact field.  For inexact
          matrices consult the numerical or symbolic matrix classes.

        - ``format`` - default: ``None``

          - ``'all'`` - attempts to create every eigenspace.  This will
            always be possible for matrices with rational entries.
          - ``'galois'`` - for each irreducible factor of the characteristic
            polynomial, a single eigenspace will be output for a
            single root/eigenvalue for the irreducible factor.
          - ``None`` - Uses the 'all' format if the base ring is contained
            in an algebraically closed field which is implemented.
            Otherwise, uses the 'galois' format.

        - ``var`` - default: 'a' - variable name used to
          represent elements of the root field of each
          irreducible factor of the characteristic polynomial.
          If var='a', then the root fields will be in terms of
          a0, a1, a2, ...., where the numbering runs across all
          the irreducible factors of the characteristic polynomial,
          even for linear factors.

        - ``algebraic_multiplicity`` - default: False - whether or
          not to include the algebraic multiplicity of each eigenvalue
          in the output.  See the discussion below.

        OUTPUT:

        If algebraic_multiplicity=False, return a list of pairs (e, V)
        where e is an eigenvalue of the matrix, and V is the corresponding
        left eigenspace.  For Galois conjugates of eigenvalues, there
        may be just one representative eigenspace, depending on the
        ``format`` keyword.

        If algebraic_multiplicity=True, return a list of triples (e, V, n)
        where e and V are as above and n is the algebraic multiplicity of
        the eigenvalue.

        .. warning::

           Uses a somewhat naive algorithm (simply factors the
           characteristic polynomial and computes kernels directly
           over the extension field).

        EXAMPLES:

        Right eigenspaces are computed from the left eigenspaces of the
        transpose of the matrix.  As such, there is a greater collection
        of illustrative examples at the :meth:`eigenspaces_left`.

        We compute the right eigenspaces of a `3\times 3` rational matrix.  ::

            sage: A = matrix(QQ, 3 ,3, range(9)); A
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: A.eigenspaces_right()
            [
            (0, Vector space of degree 3 and dimension 1 over Rational Field
            User basis matrix:
            [ 1 -2  1]),
            (-1.348469228349535?, Vector space of degree 3 and dimension 1 over Algebraic Field
            User basis matrix:
            [                   1  0.1303061543300932? -0.7393876913398137?]),
            (13.34846922834954?, Vector space of degree 3 and dimension 1 over Algebraic Field
            User basis matrix:
            [                 1 3.069693845669907? 5.139387691339814?])
            ]
            sage: es = A.eigenspaces_right(format='galois'); es
            [
            (0, Vector space of degree 3 and dimension 1 over Rational Field
            User basis matrix:
            [ 1 -2  1]),
            (a1, Vector space of degree 3 and dimension 1 over Number Field in a1 with defining polynomial x^2 - 12*x - 18
            User basis matrix:
            [           1 1/5*a1 + 2/5 2/5*a1 - 1/5])
            ]
            sage: es = A.eigenspaces_right(format='galois', algebraic_multiplicity=True); es
            [
            (0, Vector space of degree 3 and dimension 1 over Rational Field
            User basis matrix:
            [ 1 -2  1], 1),
            (a1, Vector space of degree 3 and dimension 1 over Number Field in a1 with defining polynomial x^2 - 12*x - 18
            User basis matrix:
            [           1 1/5*a1 + 2/5 2/5*a1 - 1/5], 1)
            ]
            sage: e, v, n = es[0]; v = v.basis()[0]
            sage: delta = v*e - A*v
            sage: abs(abs(delta)) < 1e-10
            True

        The same computation, but with implicit base change to a field::

            sage: A = matrix(ZZ, 3, range(9)); A
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: A.eigenspaces_right(format='galois')
            [
            (0, Vector space of degree 3 and dimension 1 over Rational Field
            User basis matrix:
            [ 1 -2  1]),
            (a1, Vector space of degree 3 and dimension 1 over Number Field in a1 with defining polynomial x^2 - 12*x - 18
            User basis matrix:
            [           1 1/5*a1 + 2/5 2/5*a1 - 1/5])
            ]

        This method is only applicable to exact matrices.
        The "eigenmatrix" routines for matrices with double-precision
        floating-point entries (``RDF``, ``CDF``) are the best
        alternative.  (Since some platforms return eigenvectors
        that are the negatives of those given here, this one example
        is not tested here.)  There are also "eigenmatrix" routines for
        matrices with symbolic entries.  ::

            sage: B = matrix(RR, 3, 3, range(9))
            sage: B.eigenspaces_right()
            Traceback (most recent call last):
            ...
            NotImplementedError: eigenspaces cannot be computed reliably for inexact rings such as Real Field with 53 bits of precision,
            consult numerical or symbolic matrix classes for other options

            sage: em = B.change_ring(RDF).eigenmatrix_right()
            sage: eigenvalues = em[0]; eigenvalues.dense_matrix().zero_at(1e-15)
            [ 13.3484692283            0.0            0.0]
            [           0.0 -1.34846922835            0.0]
            [           0.0            0.0            0.0]
            sage: eigenvectors = em[1]; eigenvectors # not tested
            [ 0.164763817...  0.799699663...  0.408248290...]
            [ 0.505774475...  0.104205787... -0.816496580...]
            [ 0.846785134... -0.591288087...  0.408248290...]

            sage: x, y = var('x y')
            sage: S = matrix([[x, y], [y, 3*x^2]])
            sage: em = S.eigenmatrix_right()
            sage: eigenvalues = em[0]; eigenvalues
            [3/2*x^2 + 1/2*x - 1/2*sqrt(9*x^4 - 6*x^3 + x^2 + 4*y^2)                                                       0]
            [                                                      0 3/2*x^2 + 1/2*x + 1/2*sqrt(9*x^4 - 6*x^3 + x^2 + 4*y^2)]

            sage: eigenvectors = em[1]; eigenvectors
            [                                                    1                                                     1]
            [1/2*(3*x^2 - x - sqrt(9*x^4 - 6*x^3 + x^2 + 4*y^2))/y 1/2*(3*x^2 - x + sqrt(9*x^4 - 6*x^3 + x^2 + 4*y^2))/y]


        TESTS::

            sage: B = matrix(QQ, 2, 3, range(6))
            sage: B.eigenspaces_right()
            Traceback (most recent call last):
            ...
            TypeError: matrix must be square, not 2 x 3

            sage: B = matrix(QQ, 4, 4, range(16))
            sage: B.eigenspaces_right(format='junk')
            Traceback (most recent call last):
            ...
            ValueError: format keyword must be None, 'all' or 'galois', not junk

            sage: B.eigenspaces_right(algebraic_multiplicity='garbage')
            Traceback (most recent call last):
            ...
            ValueError: algebraic_multiplicity keyword must be True or False
        """
        if not algebraic_multiplicity in [True, False]:
            msg = 'algebraic_multiplicity keyword must be True or False'
            raise ValueError(msg.format(algebraic_multiplicity))
        if not self.is_square():
            msg = 'matrix must be square, not {0} x {1}'
            raise TypeError(msg.format(self.nrows(), self.ncols()))
        if not self.base_ring().is_exact():
            msg = ("eigenspaces cannot be computed reliably for inexact rings such as {0},\n",
                   "consult numerical or symbolic matrix classes for other options")
            raise NotImplementedError(''.join(msg).format(self.base_ring()))

        format = self._eigenspace_format(format)

        key = 'eigenspaces_right_' + format + '{0}'.format(var)
        x = self.fetch(key)
        if not x is None:
            if algebraic_multiplicity:
                return x
            else:
                return Sequence([(e[0],e[1]) for e in x], cr=True, check=False)

        V = self.transpose().eigenspaces_left(format=format, var=var, algebraic_multiplicity=True)

        self.cache(key, V)
        if algebraic_multiplicity:
            return V
        else:
            return Sequence([(e[0],e[1]) for e in V], cr=True, check=False)

    right_eigenspaces = eigenspaces_right

    def eigenvalues(self,extend=True):
        r"""
        Return a sequence of the eigenvalues of a matrix, with
        multiplicity. If the eigenvalues are roots of polynomials in QQ,
        then QQbar elements are returned that represent each separate
        root.

        If the option extend is set to False, only eigenvalues in the base
        ring are considered.

        EXAMPLES::

            sage: a = matrix(QQ, 4, range(16)); a
            [ 0  1  2  3]
            [ 4  5  6  7]
            [ 8  9 10 11]
            [12 13 14 15]
            sage: sorted(a.eigenvalues(), reverse=True)
            [32.46424919657298?, 0, 0, -2.464249196572981?]

        ::

            sage: a=matrix([(1, 9, -1, -1), (-2, 0, -10, 2), (-1, 0, 15, -2), (0, 1, 0, -1)])
            sage: a.eigenvalues()
            [-0.9386318578049146?, 15.50655435353258?, 0.2160387521361705? - 4.713151979747493?*I, 0.2160387521361705? + 4.713151979747493?*I]

        A symmetric matrix a+a.transpose() should have real eigenvalues

        ::

            sage: b=a+a.transpose()
            sage: ev = b.eigenvalues(); ev
            [-8.35066086057957?, -1.107247901349379?, 5.718651326708515?, 33.73925743522043?]

        The eigenvalues are elements of QQbar, so they really represent
        exact roots of polynomials, not just approximations.

        ::

            sage: e = ev[0]; e
            -8.35066086057957?
            sage: p = e.minpoly(); p
            x^4 - 30*x^3 - 171*x^2 + 1460*x + 1784
            sage: p(e) == 0
            True

        To perform computations on the eigenvalue as an element of a number
        field, you can always convert back to a number field element.

        ::

            sage: e.as_number_field_element()
            (Number Field in a with defining polynomial y^4 - 2*y^3 - 507*y^2 - 3972*y - 4264,
            a + 7,
            Ring morphism:
              From: Number Field in a with defining polynomial y^4 - 2*y^3 - 507*y^2 - 3972*y - 4264
              To:   Algebraic Real Field
              Defn: a |--> -15.35066086057957?)

        Notice the effect of the extend option.

        ::

            sage: M=matrix(QQ,[[0,-1,0],[1,0,0],[0,0,2]])
            sage: M.eigenvalues()
            [2, -1*I, 1*I]
            sage: M.eigenvalues(extend=False)
            [2]

        """
        x = self.fetch('eigenvalues')
        if not x is None:
            if not extend:
                x1=Sequence([])
                for i in x:
                    if i in self.base_ring():
                        x1.append(i)
                x=x1
            return x

        if not self.base_ring().is_exact():
            from warnings import warn
            warn("Using generic algorithm for an inexact ring, which will probably give incorrect results due to numerical precision issues.")

        from sage.rings.qqbar import QQbar
        G = self.fcp()   # factored charpoly of self.
        V = []
        i=0
        for h, e in G:
            if h.degree() == 1:
                alpha = [-h[0]/h[1]]
                V.extend(alpha*e)
            else:
                if extend:
                    F = h.root_field('%s%s'%('a',i))
                    try:
                        alpha = F.gen(0).galois_conjugates(QQbar)
                    except AttributeError:
                        raise NotImplementedError, "eigenvalues() is not implemented for matrices with eigenvalues that are not in the fraction field of the base ring or in QQbar"
                    V.extend(alpha*e)
            i+=1
        V = Sequence(V)
        if extend:
            self.cache('eigenvalues', V)
        if not extend:
            V1=Sequence([])
            for i in V:
                if i in self.base_ring():
                    V1.append(i)
            V=V1
        return V



    def eigenvectors_left(self,extend=True):
        r"""
        Compute the left eigenvectors of a matrix.

        For each distinct eigenvalue, returns a list of the form (e,V,n)
        where e is the eigenvalue, V is a list of eigenvectors forming a
        basis for the corresponding left eigenspace, and n is the algebraic
        multiplicity of the eigenvalue.

        If the option extend is set to False, then only the eigenvalues that
        live in the base ring are considered.

        EXAMPLES: We compute the left eigenvectors of a `3\times 3`
        rational matrix.

        ::

            sage: A = matrix(QQ,3,3,range(9)); A
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: es = A.eigenvectors_left(); es
            [(0, [
            (1, -2, 1)
            ], 1),
            (-1.348469228349535?, [(1, 0.3101020514433644?, -0.3797958971132713?)], 1),
            (13.34846922834954?, [(1, 1.289897948556636?, 1.579795897113272?)], 1)]
            sage: eval, [evec], mult = es[0]
            sage: delta = eval*evec - evec*A
            sage: abs(abs(delta)) < 1e-10
            True

        Notice the difference between considering ring extensions or not.

        ::

            sage: M=matrix(QQ,[[0,-1,0],[1,0,0],[0,0,2]])
            sage: M.eigenvectors_left()
            [(2, [
            (0, 0, 1)
            ], 1), (-1*I, [(1, -1*I, 0)], 1), (1*I, [(1, 1*I, 0)], 1)]
            sage: M.eigenvectors_left(extend=False)
            [(2, [
            (0, 0, 1)
            ], 1)]

        """
        x = self.fetch('eigenvectors_left')
        if not x is None:
            return x

        if not self.base_ring().is_exact():
            from warnings import warn
            warn("Using generic algorithm for an inexact ring, which may result in garbage from numerical precision issues.")

        V = []
        from sage.rings.qqbar import QQbar
        from sage.categories.homset import hom
        eigenspaces = self.eigenspaces_left(format='galois', algebraic_multiplicity=True)
        evec_list=[]
        n = self._nrows
        evec_eval_list = []
        F = self.base_ring().fraction_field()
        for ev in eigenspaces:
            eigval = ev[0]
            eigbasis = ev[1].basis()
            eigmult = ev[2]
            if eigval in self.base_ring() or extend:
                if eigval.parent().fraction_field() == F:
                    evec_eval_list.append((eigval, eigbasis, eigmult))
                else:
                    try:
                        eigval_conj = eigval.galois_conjugates(QQbar)
                    except AttributeError:
                        raise NotImplementedError, "eigenvectors are not implemented for matrices with eigenvalues that are not in the fraction field of the base ring or in QQbar"

                    for e in eigval_conj:
                        m = hom(eigval.parent(), e.parent(), e)
                        space = (e.parent())**n
                        evec_list = [(space)([m(i) for i in v]) for v in eigbasis]
                        evec_eval_list.append( (e, evec_list, eigmult))

        return evec_eval_list

    left_eigenvectors = eigenvectors_left

    def eigenvectors_right(self, extend=True):
        r"""
        Compute the right eigenvectors of a matrix.

        For each distinct eigenvalue, returns a list of the form (e,V,n)
        where e is the eigenvalue, V is a list of eigenvectors forming a
        basis for the corresponding right eigenspace, and n is the
        algebraic multiplicity of the eigenvalue. If ``extend = True``
        (the default), this will return eigenspaces over the algebraic
        closure of the base field where this is implemented; otherwise
        it will restrict to eigenvalues in the base field.

        EXAMPLES: We compute the right eigenvectors of a
        `3\times 3` rational matrix.

        ::

            sage: A = matrix(QQ,3,3,range(9)); A
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: es = A.eigenvectors_right(); es
            [(0, [
            (1, -2, 1)
            ], 1),
            (-1.348469228349535?, [(1, 0.1303061543300932?, -0.7393876913398137?)], 1),
            (13.34846922834954?, [(1, 3.069693845669907?, 5.139387691339814?)], 1)]
            sage: A.eigenvectors_right(extend=False)
            [(0, [
            (1, -2, 1)
            ], 1)]
            sage: eval, [evec], mult = es[0]
            sage: delta = eval*evec - A*evec
            sage: abs(abs(delta)) < 1e-10
            True
        """
        return self.transpose().eigenvectors_left(extend=extend)

    right_eigenvectors = eigenvectors_right

    def eigenmatrix_left(self):
        r"""
        Return matrices D and P, where D is a diagonal matrix of
        eigenvalues and P is the corresponding matrix where the rows are
        corresponding eigenvectors (or zero vectors) so that P\*self =
        D\*P.

        EXAMPLES::

            sage: A = matrix(QQ,3,3,range(9)); A
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: D, P = A.eigenmatrix_left()
            sage: D
            [                  0                   0                   0]
            [                  0 -1.348469228349535?                   0]
            [                  0                   0  13.34846922834954?]
            sage: P
            [                   1                   -2                    1]
            [                   1  0.3101020514433644? -0.3797958971132713?]
            [                   1   1.289897948556636?   1.579795897113272?]
            sage: P*A == D*P
            True

        Because P is invertible, A is diagonalizable.

        ::

            sage: A == (~P)*D*P
            True

        The matrix P may contain zero rows corresponding to eigenvalues for
        which the algebraic multiplicity is greater than the geometric
        multiplicity. In these cases, the matrix is not diagonalizable.

        ::

            sage: A = jordan_block(2,3); A
            [2 1 0]
            [0 2 1]
            [0 0 2]
            sage: A = jordan_block(2,3)
            sage: D, P = A.eigenmatrix_left()
            sage: D
            [2 0 0]
            [0 2 0]
            [0 0 2]
            sage: P
            [0 0 1]
            [0 0 0]
            [0 0 0]
            sage: P*A == D*P
            True

        TESTS:

        For matrices with floating point entries, some platforms will
        return eigenvectors that are negatives of those returned by the
        majority of platforms.  This test accounts for that possibility.
        Running this test independently, without adjusting the eigenvectors
        could indicate this situation on your hardware.  ::

            sage: A = matrix(QQ, 3, 3, range(9))
            sage: em = A.change_ring(RDF).eigenmatrix_left()
            sage: evalues = em[0]; evalues.dense_matrix().zero_at(2e-15)
            [ 13.3484692283            0.0            0.0]
            [           0.0 -1.34846922835            0.0]
            [           0.0            0.0            0.0]
            sage: evectors = em[1];
            sage: for i in range(3):
            ...       scale = evectors[i,0].sign()
            ...       evectors.rescale_row(i, scale)
            sage: evectors
            [ 0.440242867...  0.567868371...  0.695493875...]
            [ 0.897878732...  0.278434036... -0.341010658...]
            [ 0.408248290... -0.816496580...  0.408248290...]
        """
        from sage.misc.flatten import flatten
        evecs = self.eigenvectors_left()
        D = sage.matrix.constructor.diagonal_matrix(flatten([[e[0]]*e[2] for e in evecs]))
        rows = []
        for e in evecs:
            rows.extend(e[1]+[e[1][0].parent().zero_vector()]*(e[2]-len(e[1])))
        P = sage.matrix.constructor.matrix(rows)
        return D,P

    left_eigenmatrix = eigenmatrix_left

    def eigenmatrix_right(self):
        r"""
        Return matrices D and P, where D is a diagonal matrix of
        eigenvalues and P is the corresponding matrix where the columns are
        corresponding eigenvectors (or zero vectors) so that self\*P =
        P\*D.

        EXAMPLES::

            sage: A = matrix(QQ,3,3,range(9)); A
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: D, P = A.eigenmatrix_right()
            sage: D
            [                  0                   0                   0]
            [                  0 -1.348469228349535?                   0]
            [                  0                   0  13.34846922834954?]
            sage: P
            [                   1                    1                    1]
            [                  -2  0.1303061543300932?   3.069693845669907?]
            [                   1 -0.7393876913398137?   5.139387691339814?]
            sage: A*P == P*D
            True

        Because P is invertible, A is diagonalizable.

        ::

            sage: A == P*D*(~P)
            True

        The matrix P may contain zero columns corresponding to eigenvalues
        for which the algebraic multiplicity is greater than the geometric
        multiplicity. In these cases, the matrix is not diagonalizable.

        ::

            sage: A = jordan_block(2,3); A
            [2 1 0]
            [0 2 1]
            [0 0 2]
            sage: A = jordan_block(2,3)
            sage: D, P = A.eigenmatrix_right()
            sage: D
            [2 0 0]
            [0 2 0]
            [0 0 2]
            sage: P
            [1 0 0]
            [0 0 0]
            [0 0 0]
            sage: A*P == P*D
            True

        TESTS:

        For matrices with floating point entries, some platforms will
        return eigenvectors that are negatives of those returned by the
        majority of platforms.  This test accounts for that possibility.
        Running this test independently, without adjusting the eigenvectors
        could indicate this situation on your hardware.  ::

            sage: B = matrix(QQ, 3, 3, range(9))
            sage: em = B.change_ring(RDF).eigenmatrix_right()
            sage: evalues = em[0]; evalues.dense_matrix().zero_at(2e-15)
            [ 13.3484692283            0.0            0.0]
            [           0.0 -1.34846922835            0.0]
            [           0.0            0.0            0.0]
            sage: evectors = em[1];
            sage: for i in range(3):
            ...       scale = evectors[0,i].sign()
            ...       evectors.rescale_col(i, scale)
            sage: evectors
            [ 0.164763817...  0.799699663...  0.408248290...]
            [ 0.505774475...  0.104205787... -0.816496580...]
            [ 0.846785134... -0.591288087...  0.408248290...]
        """
        D,P=self.transpose().eigenmatrix_left()
        return D,P.transpose()

    right_eigenmatrix = eigenmatrix_right



    #####################################################################################
    # Generic Echelon Form
    ###################################################################################


    def rref(self, *args, **kwds):
        """
        Return the reduced row echelon form of the matrix, considered
        as a matrix over a field.

        If the matrix is over a ring, then an equivalent matrix is
        constructed over the fraction field, and then row reduced.

        All arguments are passed on to :meth:``echelon_form``.

        .. note::

            Because the matrix is viewed as a matrix over a field,
            every leading coefficient of the returned matrix will be
            one and will be the only nonzero entry in its column.

        EXAMPLES::

            sage: A=matrix(3,range(9)); A
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: A.rref()
            [ 1  0 -1]
            [ 0  1  2]
            [ 0  0  0]


        Note that there is a difference between :meth:`rref` and
        :meth:`echelon_form` when the matrix is not over a field (in
        this case, the integers instead of the rational numbers)::

            sage: A.base_ring()
            Integer Ring
            sage: A.echelon_form()
            [ 3  0 -3]
            [ 0  1  2]
            [ 0  0  0]

            sage: B=random_matrix(QQ,3,num_bound=10); B
            [ -4  -3   6]
            [  5  -5 9/2]
            [3/2  -4  -7]
            sage: B.rref()
            [1 0 0]
            [0 1 0]
            [0 0 1]

        In this case, since ``B`` is a matrix over a field (the
        rational numbers), :meth:`rref` and :meth:`echelon_form` are
        exactly the same::

            sage: B.echelon_form()
            [1 0 0]
            [0 1 0]
            [0 0 1]
            sage: B.echelon_form() is B.rref()
            True

        Since :meth:`echelon_form` is not implemented for every ring,
        sometimes behavior varies, as here::

            sage: R.<x>=ZZ[]
            sage: C = matrix(3,[2,x,x^2,x+1,3-x,-1,3,2,1])
            sage: C.rref()
            [1 0 0]
            [0 1 0]
            [0 0 1]
            sage: C.base_ring()
            Univariate Polynomial Ring in x over Integer Ring
            sage: C.echelon_form()
            Traceback (most recent call last):
            ...
            NotImplementedError: Ideal Ideal (2, x + 1) of Univariate Polynomial Ring in x over Integer Ring not principal
            Echelon form not implemented over 'Univariate Polynomial Ring in x over Integer Ring'.
            sage: C = matrix(3,[2,x,x^2,x+1,3-x,-1,3,2,1/2])
            sage: C.echelon_form()
            [                               2                                x                              x^2]
            [                               0                                1            15*x^2 - 3/2*x - 31/2]
            [                               0                                0 5/2*x^3 - 15/4*x^2 - 9/4*x + 7/2]
            sage: C.rref()
            [1 0 0]
            [0 1 0]
            [0 0 1]
            sage: C = matrix(3,[2,x,x^2,x+1,3-x,-1/x,3,2,1/2])
            sage: C.echelon_form()
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        R=self.base_ring()
        if R.is_field():
            return self.echelon_form()
        else:
            F=R.fraction_field()
            return self.change_ring(F).echelon_form()

    def _echelonize_ring(self, **kwds):
        r"""
        Echelonize self in place, where the base ring of self is assumed to
        be a ring (not a field).

        Right now this *only* works over ZZ and some principal ideal domains;
        otherwise a ``NotImplementedError`` is raised. In the special case of
        sparse matrices over ZZ it makes them dense, gets the echelon form of
        the dense matrix, then sets self equal to the result.

        EXAMPLES::

            sage: a = matrix(ZZ, 3, 4, [1..12], sparse=True); a
            [ 1  2  3  4]
            [ 5  6  7  8]
            [ 9 10 11 12]
            sage: a._echelonize_ring()
            sage: a
            [ 1  2  3  4]
            [ 0  4  8 12]
            [ 0  0  0  0]

            sage: L.<w> = NumberField(x^2 - x + 2)
            sage: OL = L.ring_of_integers()
            sage: m = matrix(OL, 2, 2, [1,2,3,4+w])
            sage: m.echelon_form()
            [    1     2]
            [    0 w - 2]
            sage: E, T = m.echelon_form(transformation=True); E,T
            (
            [    1     2]  [ 1  0]
            [    0 w - 2], [-3  1]
            )
            sage: E == T*m
            True

        TESTS:

        Check that http://trac.sagemath.org/sage_trac/ticket/11558 is fixed::

            sage: matrix(ZZ, [[1,2],[4,6]], sparse=False).echelon_form(transformation=True)
            (
            [1 0]  [-3  1]
            [0 2], [ 4 -1]
            )
            sage: matrix(ZZ, [[1,2],[4,6]], sparse=True).echelon_form(transformation=True)
            (
            [1 0]  [-3  1]
            [0 2], [ 4 -1]
            )
        """
        self.check_mutability()
        cdef Matrix d, a
        cdef Py_ssize_t r, c
        cdef bint transformation = kwds.has_key('transformation') and kwds['transformation']
        if self._base_ring == ZZ:
            if kwds.has_key('include_zero_rows') and not kwds['include_zero_rows']:
                raise ValueError, "cannot echelonize in place and delete zero rows"
            if transformation:
                d, a = self.dense_matrix().echelon_form(**kwds)
            else:
                d = self.dense_matrix().echelon_form(**kwds)
            for c from 0 <= c < self.ncols():
                for r from 0 <= r < self.nrows():
                    self.set_unsafe(r, c, d.get_unsafe(r,c))
            self.clear_cache()
            self.cache('pivots', d.pivots())
            self.cache('in_echelon_form', True)
        else:
            try:
                a, d, p = self._echelon_form_PID()
            except TypeError, msg:
                raise NotImplementedError, "%s\nechelon form over %s not yet implemented"%(msg, self.base_ring())

            for c from 0 <= c < self.ncols():
               for r from 0 <= r < self.nrows():
                    self.set_unsafe(r, c, d.get_unsafe(r,c))
            self.clear_cache()
            self.cache('pivots', tuple(p))
            self.cache('in_echelon_form', True)
        if transformation:
            return a
        else:
            return

    def echelonize(self, algorithm="default", cutoff=0, **kwds):
        r"""
        Transform ``self`` into a matrix in echelon form over the same
        base ring as self.

        .. note::

            This row reduction does not use division if the
            matrix is not over a field (e.g., if the matrix is over
            the integers).  If you want to calculate the echelon form
            using division, then use :meth:`rref`, which assumes that
            the matrix entries are in a field (specifically, the field
            of fractions of the base ring of the matrix).

        INPUT:

        - ``algorithm`` -- string. Which algorithm to use. Choices are

          - ``'default'``: Let Sage choose an algorithm (default).

          - ``'classical'``: Gauss elimination.

          - ``'strassen'``: use a Strassen divide and conquer
            algorithm (if available)

        - ``cutoff`` -- integer. Only used if the Strassen algorithm
          is selected.

        - ``transformation`` -- boolean. Whether to also return the
          transformation matrix. Some matrix backends do not provide
          this information, in which case this option is ignored.

        OUTPUT:

        The matrix ``self`` is put into echelon form. Nothing is
        returned unless the keyword option ``transformation=True`` is
        specified, in which case the transformation matrix is
        returned.

        EXAMPLES::

            sage: a = matrix(QQ,3,range(9)); a
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: a.echelonize()
            sage: a
            [ 1  0 -1]
            [ 0  1  2]
            [ 0  0  0]

        An immutable matrix cannot be transformed into echelon form. Use
        ``self.echelon_form()`` instead::

            sage: a = matrix(QQ,3,range(9)); a.set_immutable()
            sage: a.echelonize()
            Traceback (most recent call last):
            ...
            ValueError: matrix is immutable; please change a copy instead
            (i.e., use copy(M) to change a copy of M).
            sage: a.echelon_form()
            [ 1  0 -1]
            [ 0  1  2]
            [ 0  0  0]

        Echelon form over the integers is what is also classically often
        known as Hermite normal form::

            sage: a = matrix(ZZ,3,range(9))
            sage: a.echelonize(); a
            [ 3  0 -3]
            [ 0  1  2]
            [ 0  0  0]

        We compute an echelon form both over a domain and fraction field::

            sage: R.<x,y> = QQ[]
            sage: a = matrix(R, 2, [x,y,x,y])
            sage: a.echelon_form()               # not very useful? -- why two copies of the same row?
            [x y]
            [x y]

        ::

            sage: b = a.change_ring(R.fraction_field())
            sage: b.echelon_form()               # potentially useful
            [  1 y/x]
            [  0   0]

        Echelon form is not defined over arbitrary rings::

            sage: a = matrix(Integers(9),3,range(9))
            sage: a.echelon_form()
            Traceback (most recent call last):
            ...
            NotImplementedError: Echelon form not implemented over 'Ring of integers modulo 9'.

        Involving a sparse matrix::

            sage: m = matrix(3,[1, 1, 1, 1, 0, 2, 1, 2, 0], sparse=True); m
            [1 1 1]
            [1 0 2]
            [1 2 0]
            sage: m.echelon_form()
            [ 1  0  2]
            [ 0  1 -1]
            [ 0  0  0]
            sage: m.echelonize(); m
            [ 1  0  2]
            [ 0  1 -1]
            [ 0  0  0]

        The transformation matrix is optionally returned::

            sage: m_original = m
            sage: transformation_matrix = m.echelonize(transformation=True)
            sage: m == transformation_matrix * m_original
            True
        """
        self.check_mutability()

        if algorithm == 'default':
            if self._will_use_strassen_echelon():
                algorithm = 'strassen'
            else:
                algorithm = 'classical'
        try:
            if self.base_ring().is_field():
                if algorithm == 'classical':
                    self._echelon_in_place_classical()
                elif algorithm == 'strassen':
                    self._echelon_strassen(cutoff)
                else:
                    raise ValueError, "Unknown algorithm '%s'"%algorithm
            else:
                if not (algorithm in ['classical', 'strassen']):
                    kwds['algorithm'] = algorithm
                return self._echelonize_ring(**kwds)
        except ArithmeticError, msg:
            raise NotImplementedError, "%s\nEchelon form not implemented over '%s'."%(msg,self.base_ring())

    def echelon_form(self, algorithm="default", cutoff=0, **kwds):
        """
        Return the echelon form of self.

        .. note::

            This row reduction does not use division if the
            matrix is not over a field (e.g., if the matrix is over
            the integers).  If you want to calculate the echelon form
            using division, then use :meth:`rref`, which assumes that
            the matrix entries are in a field (specifically, the field
            of fractions of the base ring of the matrix).

        INPUT:

        - ``algorithm`` -- string. Which algorithm to use. Choices are

          - ``'default'``: Let Sage choose an algorithm (default).

          - ``'classical'``: Gauss elimination.

          - ``'strassen'``: use a Strassen divide and conquer
            algorithm (if available)

        - ``cutoff`` -- integer. Only used if the Strassen algorithm is selected.

        - ``transformation`` -- boolean. Whether to also return the
          transformation matrix. Some matrix backends do not provide
          this information, in which case this option is ignored.

        OUTPUT:

        The reduced row echelon form of ``self``, as an immutable
        matrix. Note that self is *not* changed by this command. Use
        :meth:`echelonize` to change ``self`` in place.

        If the optional parameter ``transformation=True`` is
        specified, the output consists of a pair `(E,T)` of matrices
        where `E` is the echelon form of ``self`` and `T` is the
        transformation matrix.

        EXAMPLES::

            sage: MS = MatrixSpace(GF(19),2,3)
            sage: C = MS.matrix([1,2,3,4,5,6])
            sage: C.rank()
            2
            sage: C.nullity()
            0
            sage: C.echelon_form()
            [ 1  0 18]
            [ 0  1  2]

        The matrix library used for `\ZZ/p`-matrices does not return
        the transformation matrix, so the ``transformation`` option is
        ignored::

            sage: C.echelon_form(transformation=True)
            [ 1  0 18]
            [ 0  1  2]

            sage: D = matrix(ZZ, 2, 3, [1,2,3,4,5,6])
            sage: D.echelon_form(transformation=True)
            (
            [1 2 3]  [ 1  0]
            [0 3 6], [ 4 -1]
            )
            sage: E, T = D.echelon_form(transformation=True)
            sage: T*D == E
            True
        """
        cdef bint transformation = (kwds.has_key('transformation') and kwds['transformation'])
        x = self.fetch('echelon_form')
        if x is not None:
            if not transformation:
                return x
            y = self.fetch('echelon_transformation')
            if y:
                return (x, y)

        E = self.__copy__()
        if algorithm == 'default':
            v = E.echelonize(cutoff=cutoff, **kwds)
        else:
            v = E.echelonize(algorithm = algorithm, cutoff=cutoff, **kwds)
        E.set_immutable()  # so we can cache the echelon form.
        self.cache('echelon_form', E)
        if v is not None:
            self.cache('echelon_transformation', v)
        self.cache('pivots', E.pivots())
        if transformation and v is not None:
            return (E, v)
        else:
            return E

    def _echelon_classical(self):
        """
        Return the echelon form of self.
        """
        E = self.fetch('echelon_classical')
        if not E is None:
            return E
        E = self.__copy__()
        E._echelon_in_place_classical()
        self.cache('echelon_classical', E)
        return E

    def _echelon_in_place_classical(self):
        """
        Transform self into echelon form and set the pivots of self.

        EXAMPLES::

            sage: t = matrix(QQ, 3, range(9)); t
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: E = t._echelon_in_place_classical(); t
            [ 1  0 -1]
            [ 0  1  2]
            [ 0  0  0]
        """
        tm = verbose('generic in-place Gauss elimination on %s x %s matrix'%(self._nrows, self._ncols))
        cdef Py_ssize_t start_row, c, r, nr, nc, i
        if self.fetch('in_echelon_form'):
            return

        self.check_mutability()
        cdef Matrix A

        nr = self._nrows
        nc = self._ncols
        A = self

        start_row = 0
        pivots = []

        for c from 0 <= c < nc:
            sig_check()
            for r from start_row <= r < nr:
                if A.get_unsafe(r, c):
                    pivots.append(c)
                    a_inverse = ~A.get_unsafe(r,c)
                    A.rescale_row(r, a_inverse, c)
                    A.swap_rows(r, start_row)
                    for i from 0 <= i < nr:
                        if i != start_row:
                            if A.get_unsafe(i,c):
                                minus_b = -A.get_unsafe(i, c)
                                A.add_multiple_of_row(i, start_row, minus_b, c)
                    start_row = start_row + 1
                    break
        self.cache('pivots', tuple(pivots))
        self.cache('in_echelon_form', True)
        self.cache('echelon_form', self)
        verbose('done with gauss echelon form', tm)

    def extended_echelon_form(self, subdivide=False, **kwds):
        r"""
        Returns the echelon form of ``self`` augmented with an identity matrix.

        INPUT:

        - ``subdivide`` - default: ``False`` - determines if the
          returned matrix is subdivided.  See the description of the
          (output) below for details.
        - ``kwds`` - additional keywords that can be passed to
          the method that computes the echelon form.

        OUTPUT:

        If `A` is an `m\times n` matrix, add the `m` columns of an
        `m\times m` identity matrix to the right of ``self``.  Then
        row-reduce this `m\times(n+m)` matrix.  This matrix is returned
        as an immutable matrix.

        If ``subdivide`` is ``True`` then the returned matrix has a single
        division among the columns and a single division among the rows.
        The column subdivision has `n` columns to the left and `m`
        columns to the right.  The row division separates the non-zero rows
        from the zero rows, when restricted to the first `n` columns.

        For a nonsingular matrix the final `m` columns of the extended
        echelon form are the inverse of ``self``.  For a matrix of
        any size, the final `m` columns provide a matrix that transforms
        ``self`` to echelon form when it multiplies ``self`` from the left.
        When the base ring is a field, the uniqueness of reduced row-echelon
        form implies that this transformation matrix can be taken as the
        coefficients giving a
        canonical set of linear combinations of the rows of ``self`` that
        yield reduced row-echelon form.

        When subdivided as described above, and again over a field, the
        parts of the subdivision in the upper-left corner and lower-right
        corner satisfy several interesting relationships with the row space,
        column space, left kernel and right kernel of ``self``.  See the
        examples below.

        .. note::

            This method returns an echelon form.  If the base ring
            is not a field, no atttempt is made to move to the fraction field.
            See an example below where the base ring is changed manually.

        EXAMPLES:

        The four relationships at the end of this example hold in general. ::

            sage: A = matrix(QQ, [[2, -1, 7, -1, 0, -3],
            ...                   [-1, 1, -5, 3, 4, 4],
            ...                   [2, -1, 7, 0, 2, -2],
            ...                   [2, 0, 4, 3, 6, 1],
            ...                   [2, -1, 7, 0, 2, -2]])
            sage: E = A.extended_echelon_form(subdivide=True); E
            [   1    0    2    0    0   -1|   0   -1    0    1   -1]
            [   0    1   -3    0   -2    0|   0   -2    0    2   -3]
            [   0    0    0    1    2    1|   0  2/3    0 -1/3  2/3]
            [-----------------------------+------------------------]
            [   0    0    0    0    0    0|   1  2/3    0 -1/3 -1/3]
            [   0    0    0    0    0    0|   0    0    1    0   -1]
            sage: J = E.matrix_from_columns(range(6,11)); J
            [   0   -1    0    1   -1]
            [   0   -2    0    2   -3]
            [   0  2/3    0 -1/3  2/3]
            [   1  2/3    0 -1/3 -1/3]
            [   0    0    1    0   -1]
            sage: J*A == A.rref()
            True
            sage: C = E.subdivision(0,0); C
            [ 1  0  2  0  0 -1]
            [ 0  1 -3  0 -2  0]
            [ 0  0  0  1  2  1]
            sage: L = E.subdivision(1,1); L
            [   1  2/3    0 -1/3 -1/3]
            [   0    0    1    0   -1]
            sage: A.right_kernel() == C.right_kernel()
            True
            sage: A.row_space() == C.row_space()
            True
            sage: A.column_space() == L.right_kernel()
            True
            sage: A.left_kernel() == L.row_space()
            True

        For a nonsingular matrix, the right half of the extended
        echelon form is the inverse matrix. ::

            sage: B = matrix(QQ, [[1,3,4], [1,4,4], [0,-2,-1]])
            sage: E = B.extended_echelon_form()
            sage: J = E.matrix_from_columns(range(3,6)); J
            [-4  5  4]
            [-1  1  0]
            [ 2 -2 -1]
            sage: J == B.inverse()
            True

        The result is in echelon form, so if the base ring is
        not a field, the leading entry of each row may not be 1.
        But you can easily change to the fraction field if necessary.  ::

            sage: A = matrix(ZZ, [[16, 20, 4, 5, -4, 13, 5],
            ...                   [10, 13, 3, -3, 7, 11, 6],
            ...                   [-12, -15, -3, -3, 2, -10, -4],
            ...                   [10, 13, 3, 3, -1, 9, 4],
            ...                   [4, 5, 1, 8, -10, 1, -1]])
            sage: E = A.extended_echelon_form(subdivide=True); E
            [ 2  0 -2  2 -9 -3 -4| 0  4 -3 -9  4]
            [ 0  1  1  2  0  1  1| 0  1  2  1  1]
            [ 0  0  0  3 -4 -1 -1| 0  3  1 -3  3]
            [--------------------+--------------]
            [ 0  0  0  0  0  0  0| 1  6  3 -6  5]
            [ 0  0  0  0  0  0  0| 0  7  2 -7  6]
            sage: J = E.matrix_from_columns(range(7,12)); J
            [ 0  4 -3 -9  4]
            [ 0  1  2  1  1]
            [ 0  3  1 -3  3]
            [ 1  6  3 -6  5]
            [ 0  7  2 -7  6]
            sage: J*A == A.echelon_form()
            True
            sage: B = A.change_ring(QQ)
            sage: B.extended_echelon_form(subdivide=True)
            [     1      0     -1      0  -19/6   -7/6   -5/3|     0      0 -89/42   -5/2    1/7]
            [     0      1      1      0    8/3    5/3    5/3|     0      0  34/21      2   -1/7]
            [     0      0      0      1   -4/3   -1/3   -1/3|     0      0   1/21      0    1/7]
            [------------------------------------------------+----------------------------------]
            [     0      0      0      0      0      0      0|     1      0    9/7      0   -1/7]
            [     0      0      0      0      0      0      0|     0      1    2/7     -1    6/7]

        Subdivided, or not, the result is immutable, so make a
        copy if you want to make changes.  ::

            sage: A = matrix(FiniteField(7), [[2,0,3], [5,5,3], [5,6,5]])
            sage: E = A.extended_echelon_form()
            sage: E.is_mutable()
            False
            sage: F = A.extended_echelon_form(subdivide=True)
            sage: F
            [1 0 0|0 4 6]
            [0 1 0|4 2 2]
            [0 0 1|5 2 3]
            [-----+-----]
            sage: F.is_mutable()
            False
            sage: G = copy(F)
            sage: G.subdivide([],[]); G
            [1 0 0 0 4 6]
            [0 1 0 4 2 2]
            [0 0 1 5 2 3]

        If you want to determine exactly which algorithm is
        used to compute the echelon form, you can add additional
        keywords to pass on to the ``echelon_form()`` routine
        employed on the augmented matrix.  Sending the flag
        ``include_zero_rows`` is a bit silly, since the
        extended echelon form will never have any zero rows. ::

            sage: A = matrix(ZZ, [[1,2], [5,0], [5,9]])
            sage: E = A.extended_echelon_form(algorithm='padic', include_zero_rows=False)
            sage: E
            [  1   0  36   1  -8]
            [  0   1   5   0  -1]
            [  0   0  45   1 -10]

        TESTS:

        The ``subdivide`` keyword is checked. ::

            sage: A = matrix(QQ, 2, range(4))
            sage: A.extended_echelon_form(subdivide='junk')
            Traceback (most recent call last):
            ...
            TypeError: subdivide must be True or False, not junk

        AUTHOR:

        - Rob Beezer (2011-02-02)
        """
        if not subdivide in [True, False]:
            raise TypeError("subdivide must be True or False, not %s" % subdivide)
        R = self.base_ring()
        import constructor
        ident = constructor.identity_matrix(R, self.nrows())
        E = self.augment(ident)
        extended = E.echelon_form(**kwds)
        if subdivide:
            from copy import copy
            extended = copy(extended)
            # pivots of E are cached from echelon form computation
            rank = len([c for c in E.pivots() if c < self.ncols()])
            extended.subdivide(rank, self.ncols())
            extended.set_immutable()
        return extended

    def weak_popov_form(self, ascend=True):
        """
        This function computes a weak Popov form of a matrix over a rational
        function field `k(x)`, for `k` a field.

        INPUT:

         - `ascend` - if True, rows of output matrix `W` are sorted so
           degree (= the maximum of the degrees of the elements in
           the row) increases monotonically, and otherwise degrees decrease.

        OUTPUT:

        A 3-tuple `(W,N,d)` consisting of:

        1. `W` - a matrix over `k(x)` giving a weak the Popov form of self
        2. `N` - a matrix over `k[x]` representing row operations used to
            transform `self` to `W`
        3. `d` - degree of respective columns of W; the degree of a column is
           the maximum of the degree of its elements

        `N` is invertible over `k(x)`. These matrices satisfy the relation
        `N*self = W`.

        EXAMPLES:

        The routine expects matrices over the rational function field, but
        other examples below show how one can provide matrices over the ring
        of polynomials (whose quotient field is the rational function field).

        ::

            sage: R.<t> = GF(3)['t']
            sage: K = FractionField(R)
            sage: M = matrix([[(t-1)^2/t],[(t-1)]])
            sage: M.weak_popov_form()
            (
            [          0]  [      t 2*t + 1]
            [(2*t + 1)/t], [      1       2], [-Infinity, 0]
            )

        If `self` is an `n x 1` matrix with at least one non-zero entry, `W` has
        a single non-zero entry and that entry is a scalar multiple of the
        greatest-common-divisor of the entries of `self`.

        ::

            sage: M1 = matrix([[t*(t-1)*(t+1)],[t*(t-2)*(t+2)],[t]])
            sage: output1 = M1.weak_popov_form()
            sage: output1
            (
            [0]  [        1         0 2*t^2 + 1]
            [0]  [        0         1 2*t^2 + 1]
            [t], [        0         0         1], [-Infinity, -Infinity, 1]
            )

        We check that the output is the same for a matrix `M` if its entries are
        rational functions intead of polynomials. We also check that the type of
        the output follows the documentation. See #9063

        ::

            sage: M2 = M1.change_ring(K)
            sage: output2 = M2.weak_popov_form()
            sage: output1 == output2
            True
            sage: output1[0].base_ring() is K
            True
            sage: output2[0].base_ring() is K
            True
            sage: output1[1].base_ring() is R
            True
            sage: output2[1].base_ring() is R
            True

        The following is the first half of example 5 in [H] *except* that we
        have transposed `self`; [H] uses column operations and we use row.

        ::

            sage: R.<t> = QQ['t']
            sage: M = matrix([[t^3 - t,t^2 - 2],[0,t]]).transpose()
            sage: M.weak_popov_form()
            (
            [      t    -t^2]  [ 1 -t]
            [t^2 - 2       t], [ 0  1], [2, 2]
            )

        The next example demonstrates what happens when `self` is a zero matrix.

        ::

            sage: R.<t> = GF(5)['t']
            sage: K = FractionField(R)
            sage: M = matrix([[K(0),K(0)],[K(0),K(0)]])
            sage: M.weak_popov_form()
            (
            [0 0]  [1 0]
            [0 0], [0 1], [-Infinity, -Infinity]
            )

        In the following example, `self` has more rows than columns.

        ::

            sage: R.<t> = QQ['t']
            sage: M = matrix([[t,t,t],[0,0,t]], ascend=False)
            sage: M.weak_popov_form()
            (
            [t t t]  [1 0]
            [0 0 t], [0 1], [1, 1]
            )

        The next example shows that M must be a matrix with
        coefficients in a rational function field `k(t)`.

        ::

            sage: M = matrix([[1,0],[1,1]])
            sage: M.weak_popov_form()
            Traceback (most recent call last):
            ...
            TypeError: the coefficients of M must lie in a univariate
            polynomial ring


        NOTES:

         - For consistency with LLL and other algorithms in sage, we have opted
           for row operations; however, references such as [H] transpose and use
           column operations.

         - There are multiple weak Popov forms of a matrix, so one may want to
           extend this code to produce a Popov form (see section 1.2 of [V]).  The
           latter is canonical, but more work to produce.

        REFERENCES:

        .. [H] F. Hess, "Computing Riemann-Roch spaces in algebraic function
          fields and related topics," J. Symbolic Comput. 33 (2002), no. 4,
          425--445.

        .. [MS] T. Mulders, A. Storjohann, "On lattice reduction for polynomial
          matrices," J. Symbolic Comput. 35 (2003), no. 4, 377--401

        .. [V] G. Villard, "Computing Popov and Hermite forms of polynomial
          matrices", ISSAC '96: Proceedings of the 1996 international symposium
          on Symbolic and algebraic computation, 250--258.

        """
        import sage.matrix.matrix_misc
        return sage.matrix.matrix_misc.weak_popov_form(self)

    #####################################################################################
    # Windowed Strassen Matrix Multiplication and Echelon
    # Precise algorithms invented and implemented by David Harvey and Robert Bradshaw
    # at William Stein's MSRI 2006 Summer Workshop on Modular Forms.
    #####################################################################################
    def _multiply_strassen(self, Matrix right, int cutoff=0):
        """
        Multiply self by the matrix right using a Strassen-based
        asymptotically fast arithmetic algorithm.

        ALGORITHM: Custom algorithm for arbitrary size matrices designed by
        David Harvey and Robert Bradshaw, based on Strassen's algorithm.

        INPUT:


        -  ``cutoff`` - integer (default: 0 - let class
           decide).


        EXAMPLES::

            sage: a = matrix(ZZ,4,range(16))
            sage: a._multiply_strassen(a,2)
            [ 56  62  68  74]
            [152 174 196 218]
            [248 286 324 362]
            [344 398 452 506]
        """
        if self._ncols != right._nrows:
            raise ArithmeticError, "Number of columns of self must equal number of rows of right."
        if not self._base_ring is right.base_ring():
            raise TypeError, "Base rings must be the same."

        if cutoff == 0:
            cutoff = self._strassen_default_cutoff(right)

        if cutoff <= 0:
            raise ValueError, "cutoff must be at least 1"

        output = self.new_matrix(self._nrows, right._ncols)
        # The following used to be a little faster, but meanwhile
        # the previous line is faster.
        #if self.is_sparse():
        #    output = self.matrix_space(self._nrows, right._ncols, sparse = True)(0)
        #else:
        #    output = self.matrix_space(self._nrows, right._ncols, sparse = False).zero_matrix().__copy__()

        self_window   = self.matrix_window()
        right_window  = right.matrix_window()
        output_window = output.matrix_window()


        import strassen
        strassen.strassen_window_multiply(output_window, self_window, right_window, cutoff)
        return output

    def _echelon_strassen(self, int cutoff=0):
        """
        In place Strassen echelon of self, and sets the pivots.

        ALGORITHM: Custom algorithm for arbitrary size matrices designed by
        David Harvey and Robert Bradshaw, based on Strassen's algorithm.

        EXAMPLES::

            sage: A = matrix(QQ, 4, range(16))
            sage: A._echelon_strassen(2)
            sage: A
            [ 1  0 -1 -2]
            [ 0  1  2  3]
            [ 0  0  0  0]
            [ 0  0  0  0]
        """
        tm = verbose('strassen echelon of %s x %s matrix'%(self._nrows, self._ncols))

        self.check_mutability()

        if not self._base_ring.is_field():
            raise ValueError, "Echelon form not defined over this base ring."

        if cutoff == 0:
            cutoff = self._strassen_default_echelon_cutoff()

        if cutoff < 1:
            raise ValueError, "cutoff must be at least 1"

        if self._nrows < cutoff or self._ncols < cutoff:
            self._echelon_in_place_classical()
            return

        import strassen
        pivots = strassen.strassen_echelon(self.matrix_window(), cutoff)
        self._set_pivots(pivots)
        verbose('done with strassen', tm)

    cpdef matrix_window(self, Py_ssize_t row=0, Py_ssize_t col=0,
                      Py_ssize_t nrows=-1, Py_ssize_t ncols=-1,
                      bint check=1):
        """
        Return the requested matrix window.

        EXAMPLES::

            sage: A = matrix(QQ, 3, range(9))
            sage: A.matrix_window(1,1, 2, 1)
            Matrix window of size 2 x 1 at (1,1):
            [0 1 2]
            [3 4 5]
            [6 7 8]

        We test the optional check flag.

        ::

            sage: matrix([1]).matrix_window(0,1,1,1, check=False)
            Matrix window of size 1 x 1 at (0,1):
            [1]
            sage: matrix([1]).matrix_window(0,1,1,1)
            Traceback (most recent call last):
            ...
            IndexError: matrix window index out of range

        Another test of bounds checking::

            sage: matrix([1]).matrix_window(1,1,1,1)
            Traceback (most recent call last):
            ...
            IndexError: matrix window index out of range
        """
        import matrix_window
        if nrows == -1:
            nrows = self._nrows - row
        if ncols == -1:
            ncols = self._ncols - col
        if check and (row < 0 or col < 0 or row + nrows > self._nrows or \
           col + ncols > self._ncols):
            raise IndexError, "matrix window index out of range"
        return matrix_window.MatrixWindow(self, row, col, nrows, ncols)

    def set_block(self, row, col, block):
        """
        Sets the sub-matrix of self, with upper left corner given by row,
        col to block.

        EXAMPLES::

            sage: A = matrix(QQ, 3, 3, range(9))/2
            sage: B = matrix(ZZ, 2, 1, [100,200])
            sage: A.set_block(0, 1, B)
            sage: A
            [  0 100   1]
            [3/2 200 5/2]
            [  3 7/2   4]

        We test that an exception is raised when the block is out of
        bounds::

            sage: matrix([1]).set_block(0,1,matrix([1]))
            Traceback (most recent call last):
            ...
            IndexError: matrix window index out of range
        """
        self.check_mutability()
        if block.base_ring() is not self.base_ring():
            block = block.change_ring(self.base_ring())
        window = self.matrix_window(row, col, block.nrows(), block.ncols(), check=True)
        window.set(block.matrix_window())

    def subdivide(self, row_lines=None, col_lines=None):
        """
        Divides self into logical submatrices which can then be queried and
        extracted. If a subdivision already exists, this method forgets the
        previous subdivision and flushes the cache.

        INPUT:


        -  ``row_lines`` - None, an integer, or a list of
           integers

        -  ``col_lines`` - None, an integer, or a list of
           integers


        OUTPUT: changes self

        .. note::

           One may also pass a tuple into the first argument which
           will be interpreted as (row_lines, col_lines)

        EXAMPLES::

            sage: M = matrix(5, 5, prime_range(100))
            sage: M.subdivide(2,3); M
            [ 2  3  5| 7 11]
            [13 17 19|23 29]
            [--------+-----]
            [31 37 41|43 47]
            [53 59 61|67 71]
            [73 79 83|89 97]
            sage: M.subdivision(0,0)
            [ 2  3  5]
            [13 17 19]
            sage: M.subdivision(1,0)
            [31 37 41]
            [53 59 61]
            [73 79 83]
            sage: M.subdivision_entry(1,0,0,0)
            31
            sage: M.subdivisions()
            ([2], [3])
            sage: M.subdivide(None, [1,3]); M
            [ 2| 3  5| 7 11]
            [13|17 19|23 29]
            [31|37 41|43 47]
            [53|59 61|67 71]
            [73|79 83|89 97]

        Degenerate cases work too.

        ::

            sage: M.subdivide([2,5], [0,1,3]); M
            [| 2| 3  5| 7 11]
            [|13|17 19|23 29]
            [+--+-----+-----]
            [|31|37 41|43 47]
            [|53|59 61|67 71]
            [|73|79 83|89 97]
            [+--+-----+-----]
            sage: M.subdivision(0,0)
            []
            sage: M.subdivision(0,1)
            [ 2]
            [13]
            sage: M.subdivide([2,2,3], [0,0,1,1]); M
            [|| 2|| 3  5  7 11]
            [||13||17 19 23 29]
            [++--++-----------]
            [++--++-----------]
            [||31||37 41 43 47]
            [++--++-----------]
            [||53||59 61 67 71]
            [||73||79 83 89 97]
            sage: M.subdivision(0,0)
            []
            sage: M.subdivision(2,4)
            [37 41 43 47]

        AUTHORS:

        - Robert Bradshaw (2007-06-14)
        """

        self.check_mutability()
        if col_lines is None and row_lines is not None and isinstance(row_lines, tuple):
            tmp = row_lines
            row_lines, col_lines = tmp
        if row_lines is None:
            row_lines = []
        elif not isinstance(row_lines, list):
            row_lines = [row_lines]
        if col_lines is None:
            col_lines = []
        elif not isinstance(col_lines, list):
            col_lines = [col_lines]
        row_lines = [0] + [int(ZZ(x)) for x in row_lines] + [self._nrows]
        col_lines = [0] + [int(ZZ(x)) for x in col_lines] + [self._ncols]
        if self._subdivisions is not None:
            self.clear_cache()
        self._subdivisions = (row_lines, col_lines)

    def subdivision(self, i, j):
        """
        Returns an immutable copy of the (i,j)th submatrix of self,
        according to a previously set subdivision.

        Before a subdivision is set, the only valid arguments are (0,0)
        which returns self.

        EXAMPLE::

            sage: M = matrix(3, 4, range(12))
            sage: M.subdivide(1,2); M
            [ 0  1| 2  3]
            [-----+-----]
            [ 4  5| 6  7]
            [ 8  9|10 11]
            sage: M.subdivision(0,0)
            [0 1]
            sage: M.subdivision(0,1)
            [2 3]
            sage: M.subdivision(1,0)
            [4 5]
            [8 9]

        It handles size-zero subdivisions as well.

        ::

            sage: M = matrix(3, 4, range(12))
            sage: M.subdivide([0],[0,2,2,4]); M
            [+-----++-----+]
            [| 0  1|| 2  3|]
            [| 4  5|| 6  7|]
            [| 8  9||10 11|]
            sage: M.subdivision(0,0)
            []
            sage: M.subdivision(1,1)
            [0 1]
            [4 5]
            [8 9]
            sage: M.subdivision(1,2)
            []
            sage: M.subdivision(1,0)
            []
            sage: M.subdivision(0,1)
            []
        """
        if self._subdivisions is None:
            self._subdivisions = ([0, self._nrows], [0, self._ncols])
        key = "subdivision %s %s"%(i,j)
        sd = self.fetch(key)
        if sd is None:
            sd = self[self._subdivisions[0][i]:self._subdivisions[0][i+1],
                      self._subdivisions[1][j]:self._subdivisions[1][j+1]]
            sd.set_immutable()
            self.cache(key, sd)
        return sd

    def subdivision_entry(self, i, j, x, y):
        """
        Returns the x,y entry of the i,j submatrix of self.

        EXAMPLES::

            sage: M = matrix(5, 5, range(25))
            sage: M.subdivide(3,3); M
            [ 0  1  2| 3  4]
            [ 5  6  7| 8  9]
            [10 11 12|13 14]
            [--------+-----]
            [15 16 17|18 19]
            [20 21 22|23 24]
            sage: M.subdivision_entry(0,0,1,2)
            7
            sage: M.subdivision(0,0)[1,2]
            7
            sage: M.subdivision_entry(0,1,0,0)
            3
            sage: M.subdivision_entry(1,0,0,0)
            15
            sage: M.subdivision_entry(1,1,1,1)
            24

        Even though this entry exists in the matrix, the index is invalid
        for the submatrix.

        ::

            sage: M.subdivision_entry(0,0,4,0)
            Traceback (most recent call last):
            ...
            IndexError: Submatrix 0,0 has no entry 4,0
        """
        if self._subdivisions is None:
            if not i and not j:
                return self[x,y]
            else:
                raise IndexError, "No such submatrix %s, %s"%(i,j)
        if x >= self._subdivisions[0][i+1]-self._subdivisions[0][i] or \
           y >= self._subdivisions[1][j+1]-self._subdivisions[1][j]:
            raise IndexError, "Submatrix %s,%s has no entry %s,%s"%(i,j, x, y)
        return self[self._subdivisions[0][i] + x , self._subdivisions[1][j] + y]

    def _subdivide_on_augment(self, left, right):
        r"""
        Helper method to manage subdivisions when augmenting a matrix.

        INPUT:

        - ``left``, ``right`` - two matrices, such that if ``left`` is
          augmented by placing ``right`` on the right side of ``left``,
          then the result is ``self``.  It is the responsibility of the
          calling routine to ensure this condition holds.

        OUTPUT:

        ``None``.  A new subdivision is created between ``left`` and
        ``right`` for ``self``.  If possible, row subdivisions are
        preserved in ``self``, but if the two sets of row subdivisions
        are incompatible, they are removed.

        EXAMPLE::

            sage: A = matrix(QQ, 3, 4, range(12))
            sage: B = matrix(QQ, 3, 6, range(18))
            sage: C = A.augment(B)
            sage: C._subdivide_on_augment(A, B)
            sage: C
            [ 0  1  2  3| 0  1  2  3  4  5]
            [ 4  5  6  7| 6  7  8  9 10 11]
            [ 8  9 10 11|12 13 14 15 16 17]

        More descriptive, but indirect, doctests are at
        :meth:`sage.matrix.matrix1.Matrix.augment`.
        """
        left_rows, left_cols = left.subdivisions()
        right_rows, right_cols = right.subdivisions()
        if left_rows == right_rows:
            self_rows = left_rows
        else:
            self_rows = None
        nc = left.ncols()
        self_cols = left_cols + [nc]
        for col in right_cols:
            self_cols.append(col+nc)
        self.subdivide(self_rows, self_cols)
        return None

    def _subdivide_on_stack(self, top, bottom):
        r"""
        Helper method to manage subdivisions when stacking a matrix.

        INPUT:

        - ``top``, ``bottom`` - two matrices, such that if ``top`` is
        stacked by placing ``top`` above ``bottom``, then the result
        is ``self``.  It is the responsibility of the calling routine
        to ensure this condition holds.

        OUTPUT:

        ``None``.  A new subdivision is created between ``top`` and
        ``bottom`` for ``self``.  If possible, column subdivisions are
        preserved in ``self``, but if the two sets of solumn subdivisions
        are incompatible, they are removed.

        EXAMPLE::

            sage: A = matrix(QQ, 3, 2, range(6))
            sage: B = matrix(QQ, 4, 2, range(8))
            sage: C = A.stack(B)
            sage: C._subdivide_on_stack(A, B)
            sage: C
            [0 1]
            [2 3]
            [4 5]
            [---]
            [0 1]
            [2 3]
            [4 5]
            [6 7]

        More descriptive, but indirect, doctests are at
        :meth:`sage.matrix.matrix1.Matrix.augment`.
        """
        top_rows, top_cols = top.subdivisions()
        bottom_rows, bottom_cols = bottom.subdivisions()
        if top_cols == bottom_cols:
            self_cols = top_cols
        else:
            self_cols = None
        nr = top.nrows()
        self_rows = top_rows + [nr]
        for row in bottom_rows:
            self_rows.append(row+nr)
        self.subdivide(self_rows, self_cols)
        return None

    def subdivisions(self):
        """
        Returns the current subdivision of self.

        EXAMPLES::

            sage: M = matrix(5, 5, range(25))
            sage: M.subdivisions()
            ([], [])
            sage: M.subdivide(2,3)
            sage: M.subdivisions()
            ([2], [3])
            sage: N = M.parent()(1)
            sage: N.subdivide(M.subdivisions()); N
            [1 0 0|0 0]
            [0 1 0|0 0]
            [-----+---]
            [0 0 1|0 0]
            [0 0 0|1 0]
            [0 0 0|0 1]
        """
        if self._subdivisions is None:
            return ([], [])
        else:
            return (self._subdivisions[0][1:-1], self._subdivisions[1][1:-1])

    # 'get_subdivisions' is kept for backwards compatibility: see #4983.
    get_subdivisions = subdivisions

    def tensor_product(self, A, subdivide=True):
        r"""
        Returns the tensor product of two matrices.

        INPUT:

        - ``A`` - a matrix
        - ``subdivide`` - default: True - whether or not to return
          natural subdivisions with the matrix

        OUTPUT:

        Replace each element of ``self`` by a copy of ``A``, but first
        create a scalar multiple of ``A`` by the element it replaces.
        So if ``self`` is an `m\times n` matrix and ``A`` is a
        `p\times q` matrix, then the tensor product is an `mp\times nq`
        matrix.  By default, the matrix will be subdivided into
        submatrices of size `p\times q`.

        EXAMPLES::

            sage: M1=Matrix(QQ,[[-1,0],[-1/2,-1]])
            sage: M2=Matrix(ZZ,[[1,-1,2],[-2,4,8]])
            sage: M1.tensor_product(M2)
            [  -1    1   -2|   0    0    0]
            [   2   -4   -8|   0    0    0]
            [--------------+--------------]
            [-1/2  1/2   -1|  -1    1   -2]
            [   1   -2   -4|   2   -4   -8]
            sage: M2.tensor_product(M1)
            [  -1    0|   1    0|  -2    0]
            [-1/2   -1| 1/2    1|  -1   -2]
            [---------+---------+---------]
            [   2    0|  -4    0|  -8    0]
            [   1    2|  -2   -4|  -4   -8]

        Subdivisions can be optionally suppressed.  ::

            sage: M1.tensor_product(M2, subdivide=False)
            [  -1    1   -2    0    0    0]
            [   2   -4   -8    0    0    0]
            [-1/2  1/2   -1   -1    1   -2]
            [   1   -2   -4    2   -4   -8]

        Different base rings are handled sensibly.  ::

            sage: A = matrix(ZZ, 2, 3, range(6))
            sage: B = matrix(FiniteField(23), 3, 4, range(12))
            sage: C = matrix(FiniteField(29), 4, 5, range(20))
            sage: D = A.tensor_product(B)
            sage: D.parent()
            Full MatrixSpace of 6 by 12 dense matrices over Finite Field of size 23
            sage: E = C.tensor_product(B)
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '*': 'Finite Field of size 29' and 'Full MatrixSpace of 3 by 4 dense matrices over Finite Field of size 23'

        The input is checked to be sure it is a matrix.  ::

            sage: A = matrix(QQ, 2, range(4))
            sage: A.tensor_product('junk')
            Traceback (most recent call last):
            ...
            TypeError: tensor product requires a second matrix, not junk
        """
        if not isinstance(A, Matrix):
            raise TypeError('tensor product requires a second matrix, not {0}'.format(A))
        return sage.matrix.constructor.block_matrix(self.nrows(),self.ncols(),[x*A for x in self.list()], subdivide=subdivide)

    def randomize(self, density=1, nonzero=False, *args, **kwds):
        """
        Randomize density proportion of the entries of this matrix, leaving
        the rest unchanged.

        .. note::

           We actually choose at random ``density`` proportion of entries of
           the matrix and set them to random elements. It's possible that the
           same position can be chosen multiple times, especially for a very
           small matrix.

        INPUT:

        -  ``density`` - ``float`` (default: 1); rough measure of the
           proportion of nonzero entries in the random matrix
        -  ``nonzero`` - Bool (default: ``False``); whether the new entries
           have to be non-zero
        -  ``*args, **kwds`` - Remaining parameters may be passed to the
           ``random_element`` function of the base ring

        EXAMPLES:

        We construct the zero matrix over a polynomial ring.

        ::

            sage: a = matrix(QQ['x'], 3); a
            [0 0 0]
            [0 0 0]
            [0 0 0]

        We then randomize roughly half the entries::

            sage: a.randomize(0.5)
            sage: a
            [      1/2*x^2 - x - 12 1/2*x^2 - 1/95*x - 1/2                      0]
            [-5/2*x^2 + 2/3*x - 1/4                      0                      0]
            [          -x^2 + 2/3*x                      0                      0]

        Now we randomize all the entries of the resulting matrix::

            sage: a.randomize()
            sage: a
            [     1/3*x^2 - x + 1             -x^2 + 1              x^2 - x]
            [ -1/14*x^2 - x - 1/4           -4*x - 1/5 -1/4*x^2 - 1/2*x + 4]
            [ 1/9*x^2 + 5/2*x - 3     -x^2 + 3/2*x + 1   -2/7*x^2 - x - 1/2]

        We create the zero matrix over the integers::

            sage: a = matrix(ZZ, 2); a
            [0 0]
            [0 0]

        Then we randomize it; the x and y parameters, which determine the
        size of the random elements, are passed onto the ZZ random_element
        method.

        ::

            sage: a.randomize(x=-2^64, y=2^64)
            sage: a
            [-12401200298100116246   1709403521783430739]
            [ -4417091203680573707  17094769731745295000]
        """
        randint = current_randstate().python_random().randint

        density = float(density)
        if density <= 0:
            return
        if density > 1:
            density = 1

        self.check_mutability()
        self.clear_cache()

        R = self.base_ring()

        cdef Py_ssize_t i, j, num

        if nonzero:
            if density >= 1:
                for i from 0 <= i < self._nrows:
                    for j from 0 <= j < self._ncols:
                        self.set_unsafe(i, j, R._random_nonzero_element(*args,\
                            **kwds))
            else:
                num = int(self._nrows * self._ncols * density)
                for i from 0 <= i < num:
                    self.set_unsafe(randint(0, self._nrows - 1),
                                    randint(0, self._ncols - 1),
                                    R._random_nonzero_element(*args, **kwds))
        else:
            if density >= 1:
                for i from 0 <= i < self._nrows:
                    for j from 0 <= j < self._ncols:
                        self.set_unsafe(i, j, R.random_element(*args, **kwds))
            else:
                num = int(self._nrows * self._ncols * density)
                for i from 0 <= i < num:
                    self.set_unsafe(randint(0, self._nrows - 1),
                                    randint(0, self._ncols - 1),
                                    R.random_element(*args, **kwds))

    def is_one(self):
        """
        Return True if this matrix is the identity matrix.

        EXAMPLES::

            sage: m = matrix(QQ,2,range(4))
            sage: m.is_one()
            False
            sage: m = matrix(QQ,2,[5,0,0,5])
            sage: m.is_one()
            False
            sage: m = matrix(QQ,2,[1,0,0,1])
            sage: m.is_one()
            True
            sage: m = matrix(QQ,2,[1,1,1,1])
            sage: m.is_one()
            False
        """
        return self.is_scalar(1)

    def is_scalar(self, a = None):
        """
        Return True if this matrix is a scalar matrix.

        INPUT

        - base_ring element a, which is chosen as self[0][0] if
          a = None

        OUTPUT

        - whether self is a scalar matrix (in fact the scalar matrix
          aI if a is input)

        EXAMPLES::

            sage: m = matrix(QQ,2,range(4))
            sage: m.is_scalar(5)
            False
            sage: m = matrix(QQ,2,[5,0,0,5])
            sage: m.is_scalar(5)
            True
            sage: m = matrix(QQ,2,[1,0,0,1])
            sage: m.is_scalar(1)
            True
            sage: m = matrix(QQ,2,[1,1,1,1])
            sage: m.is_scalar(1)
            False
        """
        if not self.is_square():
            return False
        cdef Py_ssize_t i, j
        if a is None:
            if self._nrows == 0:
                return True
            a = self.get_unsafe(0,0)
        else:
            a = self.base_ring()(a)
        zero = self.base_ring()(0)
        for i from 0 <= i < self._nrows:
            for j from 0 <= j < self._ncols:
                if i != j:
                    if self.get_unsafe(i,j) != zero:
                        return False
                else:
                    if self.get_unsafe(i, i) != a:
                        return False
        return True


    def is_unitary(self):
        r"""
        Returns ``True`` if the columns of the matrix are an orthonormal basis.

        For a matrix with real entries this determines if a matrix is
        "orthogonal" and for a matrix with complex entries this determines
        if the matrix is "unitary."

        OUTPUT:

        ``True`` if the matrix is square and its conjugate-transpose is
        its inverse, and ``False`` otherwise.  In other words, a matrix
        is orthogonal or unitary if the product of its conjugate-transpose
        times the matrix is the identity matrix.

        For numerical matrices a specialized routine available
        over ``RDF`` and ``CDF`` is a good choice.

        EXAMPLES::

            sage: A = matrix(QQbar, [[(1/sqrt(5))*(1+i), (1/sqrt(55))*(3+2*I), (1/sqrt(22))*(2+2*I)],
            ...                      [(1/sqrt(5))*(1-i), (1/sqrt(55))*(2+2*I),  (1/sqrt(22))*(-3+I)],
            ...                      [    (1/sqrt(5))*I, (1/sqrt(55))*(3-5*I),    (1/sqrt(22))*(-2)]])
            sage: A.is_unitary()
            True

        A permutation matrix is always orthogonal. ::

            sage: sigma = Permutation([1,3,4,5,2])
            sage: P = sigma.to_matrix(); P
            [1 0 0 0 0]
            [0 0 0 0 1]
            [0 1 0 0 0]
            [0 0 1 0 0]
            [0 0 0 1 0]
            sage: P.is_unitary()
            True
            sage: P.change_ring(GF(3)).is_unitary()
            True
            sage: P.change_ring(GF(3)).is_unitary()
            True

        A square matrix far from unitary. ::

            sage: A = matrix(QQ, 4, range(16))
            sage: A.is_unitary()
            False

        Rectangular matrices are never unitary.  ::

            sage: A = matrix(QQbar, 3, 4)
            sage: A.is_unitary()
            False
        """
        import sage.matrix.constructor
        if not self.is_square():
            return False
        if hasattr(self.base_ring().an_element(), 'conjugate'):
            P = self.conjugate().transpose()*self    # Unitary
        else:
            P = self.transpose()*self                # Orthogonal
        return P.is_scalar(1)

    def is_bistochastic(self, normalized = True):
        r"""
        Returns ``True`` if this matrix is bistochastic.

        A matrix is said to be bistochastic if both the sums of the
        entries of each row and the sum of the entries of each column
        are equal to 1.

        INPUT:

        - ``normalized`` -- if set to ``True`` (default), checks that
          the sums are equal to 1. When set to ``False``, checks that
          the row sums and column sums are all equal to some constant
          possibly different from 1.

        EXAMPLES:

        The identity matrix is clearly bistochastic::

            sage: Matrix(5,5,1).is_bistochastic()
            True

        The same matrix, multiplied by 2, is not bistochastic anymore,
        though is verifies the constraints of ``normalized == False``::

            sage: (2 * Matrix(5,5,1)).is_bistochastic()
            False
            sage: (2 * Matrix(5,5,1)).is_bistochastic(normalized = False)
            True
        """

        row_sums = map(sum, self.rows())
        col_sums = map(sum, self.columns())

        return self.is_square() and\
                col_sums[0] == row_sums[0] and\
                row_sums == col_sums and\
                row_sums == len(row_sums) * [col_sums[0]] and\
                ((not normalized) or col_sums[0] == self.base_ring()(1))

    def is_normal(self):
        r"""
        Returns ``True`` if the matrix commutes with its conjugate-transpose.

        OUTPUT:

        ``True`` if the matrix is square and commutes with its
        conjugate-transpose, and ``False`` otherwise.

        Normal matrices are precisely those that can be diagonalized
        by a unitary matrix.

        This routine is for matrices over exact rings and so may not
        work properly for matrices over ``RR`` or ``CC``.  For matrices with
        approximate entries, the rings of double-precision floating-point
        numbers, ``RDF`` and ``CDF``, are a better choice since the
        :meth:`sage.matrix.matrix_double_dense.Matrix_double_dense.is_normal`
        method has a tolerance parameter.  This provides control over
        allowing for minor discrepancies between entries when checking
        equality.

        The result is cached.

        EXAMPLES:

        Hermitian matrices are normal.  ::

            sage: A = matrix(QQ, 5, range(25)) + I*matrix(QQ, 5, range(0, 50, 2))
            sage: B = A*A.conjugate_transpose()
            sage: B.is_hermitian()
            True
            sage: B.is_normal()
            True

        Circulant matrices are normal.  ::

            sage: G = graphs.CirculantGraph(20, [3, 7])
            sage: D = digraphs.Circuit(20)
            sage: A = 3*D.adjacency_matrix() - 5*G.adjacency_matrix()
            sage: A.is_normal()
            True

        Skew-symmetric matrices are normal.  ::

            sage: A = matrix(QQ, 5, range(25))
            sage: B = A - A.transpose()
            sage: B.is_skew_symmetric()
            True
            sage: B.is_normal()
            True

        A small matrix that does not fit into any of the usual categories
        of normal matrices.  ::

            sage: A = matrix(ZZ, [[1, -1],
            ...                   [1,  1]])
            sage: A.is_normal()
            True
            sage: not A.is_hermitian() and not A.is_skew_symmetric()
            True

        Sage has several fields besides the entire complex numbers
        where conjugation is non-trivial. ::

            sage: F.<b> = QuadraticField(-7)
            sage: C = matrix(F, [[-2*b - 3,  7*b - 6, -b + 3],
            ...                  [-2*b - 3, -3*b + 2,   -2*b],
            ...                  [   b + 1,        0,     -2]])
            sage: C = C*C.conjugate_transpose()
            sage: C.is_normal()
            True

        A matrix that is nearly normal, but for a non-real
        diagonal entry. ::

            sage: A = matrix(QQbar, [[    2,   2-I, 1+4*I],
            ...                      [  2+I,   3+I, 2-6*I],
            ...                      [1-4*I, 2+6*I,     5]])
            sage: A.is_normal()
            False
            sage: A[1,1] = 132
            sage: A.is_normal()
            True

        Rectangular matrices are never normal.  ::

            sage: A = matrix(QQbar, 3, 4)
            sage: A.is_normal()
            False

        A square, empty matrix is trivially normal.  ::

            sage: A = matrix(QQ, 0, 0)
            sage: A.is_normal()
            True

        AUTHOR:

        - Rob Beezer (2011-03-31)
        """
        key = 'normal'
        n = self.fetch(key)
        if not n is None:
            return n
        if not self.is_square():
            self.cache(key, False)
            return False
        if self._nrows == 0:
            self.cache(key, True)
            return True

        CT = self.conjugate_transpose()
        cdef Matrix left = self*CT
        cdef Matrix right = CT*self

        cdef Py_ssize_t i,j
        normal = True
        # two products are Hermitian, need only check lower triangle
        for i from 0 <= i < self._nrows:
            for j from 0 <= j <= i:
                if left.get_unsafe(i,j) != right.get_unsafe(i,j):
                    normal = False
                    break
            if not normal:
                break
        self.cache(key, normal)
        return normal

    def as_sum_of_permutations(self):
        r"""
        Returns the current matrix as a sum of permutation matrices

        According to the Birkhoff-von Neumann Theorem, any bistochastic matrix
        can be written as a positive sum of permutation matrices, which also
        means that the polytope of bistochastic matrices is integer.

        As a non-bistochastic matrix can obviously not be written as a sum of
        permutations, this theorem is an equivalence.

        This function, given a bistochastic matrix, returns the corresponding
        decomposition.

        .. SEEALSO:

        - :meth:`bistochastic_as_sum_of_permutations <sage.combinat.permutation.bistochastic_as_sum_of_permutations>`
          -- for more information on this method.

        EXAMPLE:

        We create a bistochastic matrix from a convex sum of permutations, then
        try to deduce the decomposition from the matrix ::

            sage: L = []
            sage: L.append((9,Permutation([4, 1, 3, 5, 2])))
            sage: L.append((6,Permutation([5, 3, 4, 1, 2])))
            sage: L.append((3,Permutation([3, 1, 4, 2, 5])))
            sage: L.append((2,Permutation([1, 4, 2, 3, 5])))
            sage: M = sum([c * p.to_matrix() for (c,p) in L])
            sage: decomp = sage.combinat.permutation.bistochastic_as_sum_of_permutations(M)
            sage: print decomp
            2*B[[1, 4, 2, 3, 5]] + 3*B[[3, 1, 4, 2, 5]] + 9*B[[4, 1, 3, 5, 2]] + 6*B[[5, 3, 4, 1, 2]]

        An exception is raised when the matrix is not bistochastic::

            sage: M = Matrix([[2,3],[2,2]])
            sage: decomp = sage.combinat.permutation.bistochastic_as_sum_of_permutations(M)
            Traceback (most recent call last):
            ...
            ValueError: The matrix is not bistochastic
        """

        from sage.combinat.permutation import bistochastic_as_sum_of_permutations
        return bistochastic_as_sum_of_permutations(self)


    def visualize_structure(self, filename=None, maxsize=512):
        """
        Write a PNG image to 'filename' which visualizes self by putting
        black pixels in those positions which have nonzero entries.

        White pixels are put at positions with zero entries. If 'maxsize'
        is given, then the maximal dimension in either x or y direction is
        set to 'maxsize' depending on which is bigger. If the image is
        scaled, the darkness of the pixel reflects how many of the
        represented entries are nonzero. So if e.g. one image pixel
        actually represents a 2x2 submatrix, the dot is darker the more of
        the four values are nonzero.

        INPUT:


        -  ``filename`` - either a path or None in which case a
           filename in the current directory is chosen automatically
           (default:None)


        maxsize - maximal dimension in either x or y direction of the
        resulting image. If None or a maxsize larger than
        max(self.nrows(),self.ncols()) is given the image will have the
        same pixelsize as the matrix dimensions (default: 512)

        EXAMPLE::

            sage: M = random_matrix(CC, 4)
            sage: M.visualize_structure(os.path.join(SAGE_TMP, "matrix.png"))
        """
        import gd
        import os

        cdef int x, y, _x, _y, v, bi, bisq
        cdef int ir,ic
        cdef float b, fct

        mr, mc = self.nrows(), self.ncols()

        if maxsize is None:

            ir = mc
            ic = mr
            b = 1.0

        elif max(mr,mc) > maxsize:

            maxsize = float(maxsize)
            ir = int(mc * maxsize/max(mr,mc))
            ic = int(mr * maxsize/max(mr,mc))
            b = max(mr,mc)/maxsize

        else:

            ir = mc
            ic = mr
            b = 1.0

        bi = round(b)
        bisq = bi*bi
        fct = 255.0/bisq

        im = gd.image((ir,ic),1)
        white = im.colorExact((255,255,255))
        im.fill((0,0),white)

        # these speed things up a bit
        colorExact = im.colorExact
        setPixel = im.setPixel

        for x from 0 <= x < ic:
            for y from 0 <= y < ir:
                v = bisq
                for _x from 0 <= _x < bi:
                    for _y from 0 <= _y < bi:
                        if not self.get_unsafe(<int>(x*b + _x), <int>(y*b + _y)).is_zero():
                            v-=1 #increase darkness

                v =  round(v*fct)
                val = colorExact((v,v,v))
                setPixel((y,x), val)

        if filename is None:
            filename = graphics_filename()

        im.writePng(filename)

    def density(self):
        """
        Return the density of the matrix.

        By density we understand the ratio of the number of nonzero
        positions and the self.nrows() \* self.ncols(), i.e. the number of
        possible nonzero positions.

        EXAMPLE:

        First, note that the density parameter does not ensure the density
        of a matrix, it is only an upper bound.

        ::

            sage: A = random_matrix(GF(127),200,200,density=0.3)
            sage: A.density()
            5211/20000

        ::

            sage: A = matrix(QQ,3,3,[0,1,2,3,0,0,6,7,8])
            sage: A.density()
            2/3

        ::

            sage: a = matrix([[],[],[],[]])
            sage: a.density()
            0
        """
        cdef int x,y,k
        k = 0
        nr = self.nrows()
        nc = self.ncols()
        if nc == 0 or nr == 0:
            return 0
        for x from 0 <= x < nr:
            for y from 0 <= y < nc:
                if not self.get_unsafe(x,y).is_zero():
                    k+=1
        return QQ(k)/QQ(nr*nc)


    def inverse(self):
        """
        Returns the inverse of self, without changing self.

        Note that one can use the Python inverse operator to obtain the
        inverse as well.

        EXAMPLES::

            sage: m = matrix([[1,2],[3,4]])
            sage: m^(-1)
            [  -2    1]
            [ 3/2 -1/2]
            sage: m.inverse()
            [  -2    1]
            [ 3/2 -1/2]
            sage: ~m
            [  -2    1]
            [ 3/2 -1/2]

        ::

            sage: m = matrix([[1,2],[3,4]], sparse=True)
            sage: m^(-1)
            [  -2    1]
            [ 3/2 -1/2]
            sage: m.inverse()
            [  -2    1]
            [ 3/2 -1/2]
            sage: ~m
            [  -2    1]
            [ 3/2 -1/2]
            sage: m.I
            [  -2    1]
            [ 3/2 -1/2]

        TESTS::

            sage: matrix().inverse()
            []
       """
        return self.__invert__()

    def adjoint(self):
        """
        Returns the adjoint matrix of self (matrix of cofactors).

        OUTPUT:

        - ``N`` - the adjoint matrix, such that
          N \* M = M \* N = M.parent(M.det())

        ALGORITHM:

        Use PARI whenever the method ``self._adjoint`` is included to do so
        in an inheriting class.  Otherwise, use a generic division-free
        algorithm to compute the characteristic polynomial and hence the
        adjoint.

        The result is cached.

        EXAMPLES::

            sage: M = Matrix(ZZ,2,2,[5,2,3,4]) ; M
            [5 2]
            [3 4]
            sage: N = M.adjoint() ; N
            [ 4 -2]
            [-3  5]
            sage: M * N
            [14  0]
            [ 0 14]
            sage: N * M
            [14  0]
            [ 0 14]
            sage: M = Matrix(QQ,2,2,[5/3,2/56,33/13,41/10]) ; M
            [  5/3  1/28]
            [33/13 41/10]
            sage: N = M.adjoint() ; N
            [ 41/10  -1/28]
            [-33/13    5/3]
            sage: M * N
            [7363/1092         0]
            [        0 7363/1092]

        AUTHORS:

        - Unknown: No author specified in the file from 2009-06-25
        - Sebastian Pancratz (2009-06-25): Reflecting the change that
          ``_adjoint`` is now implemented in this class
        """

        if self._nrows != self._ncols:
            raise ValueError, "self must be a square matrix"

        X = self.fetch('adjoint')
        if not X is None:
            return X

        X = self._adjoint()
        self.cache('adjoint', X)
        return X

    def _adjoint(self):
        r"""
        Returns the adjoint of self.

        OUTPUT:

        - matrix - the adjoint of self

        EXAMPLES:

        Here is one example to illustrate this::

            sage: A = matrix(ZZ, [[1,24],[3,5]])
            sage: A
            [ 1 24]
            [ 3  5]
            sage: A._adjoint()
            [  5 -24]
            [ -3   1]

        Secondly, here is an example over a polynomial ring::

            sage: R.<t> = QQ[]
            sage: A = matrix(R, [[-2*t^2 + t + 3/2, 7*t^2 + 1/2*t - 1,      \
                                  -6*t^2 + t - 2/11],                       \
                                 [-7/3*t^2 - 1/2*t - 1/15, -2*t^2 + 19/8*t, \
                                  -10*t^2 + 2*t + 1/2],                     \
                                 [6*t^2 - 1/2, -1/7*t^2 + 9/4*t, -t^2 - 4*t \
                                  - 1/10]])
            sage: A
            [       -2*t^2 + t + 3/2       7*t^2 + 1/2*t - 1       -6*t^2 + t - 2/11]
            [-7/3*t^2 - 1/2*t - 1/15         -2*t^2 + 19/8*t     -10*t^2 + 2*t + 1/2]
            [            6*t^2 - 1/2        -1/7*t^2 + 9/4*t       -t^2 - 4*t - 1/10]
            sage: A._adjoint()
            [          4/7*t^4 + 1591/56*t^3 - 961/70*t^2 - 109/80*t 55/7*t^4 + 104/7*t^3 + 6123/1540*t^2 - 959/220*t - 1/10       -82*t^4 + 101/4*t^3 + 1035/88*t^2 - 29/22*t - 1/2]
            [   -187/3*t^4 + 13/6*t^3 + 57/10*t^2 - 79/60*t - 77/300            38*t^4 + t^3 - 793/110*t^2 - 28/5*t - 53/220 -6*t^4 + 44/3*t^3 + 4727/330*t^2 - 1147/330*t - 487/660]
            [          37/3*t^4 - 136/7*t^3 - 1777/840*t^2 + 83/80*t      292/7*t^4 + 107/14*t^3 - 323/28*t^2 - 29/8*t + 1/2   61/3*t^4 - 25/12*t^3 - 269/120*t^2 + 743/240*t - 1/15]

        Finally, an example over a general ring, that is to say, as of
        version 4.0.2, SAGE does not even determine that ``S`` in the following
        example is an integral domain::

            sage: R.<a,b> = QQ[]
            sage: S.<x,y> = R.quo((b^3))
            sage: A = matrix(S, [[x*y^2,2*x],[2,x^10*y]])
            sage: A
            [ x*y^2    2*x]
            [     2 x^10*y]
            sage: A.det()
            -4*x
            sage: A.charpoly('T')
            T^2 + (-x^10*y - x*y^2)*T - 4*x
            sage: A.adjoint()
            [x^10*y   -2*x]
            [    -2  x*y^2]
            sage: A.adjoint() * A
            [-4*x    0]
            [   0 -4*x]

        TESTS::

            sage: A = matrix(ZZ, 0, 0)
            sage: A
            []
            sage: A._adjoint()
            []
            sage: A = matrix(ZZ, [[2]])
            sage: A
            [2]
            sage: A._adjoint()
            [1]

        NOTES:

        The key feature of this implementation is that it is division-free.
        This means that it can be used as a generic implementation for any
        ring (commutative and with multiplicative identity).  The algorithm
        is described in full detail as Algorithm 3.1 in [Se02].

        Note that this method does not utilise a lookup if the adjoint has
        already been computed previously, and it does not cache the result.
        This is all left to the method `adjoint`.

        REFERENCES:

        - [Se02] T. R. Seifullin, "Computation of determinants, adjoint
          matrices, and characteristic polynomials without division"

        AUTHORS:

        - Sebastian Pancratz (2009-06-12): Initial version
        """

        # Validate assertions
        #
        if self._nrows != self._ncols:
            raise ValueError("self must be a square matrix")

        # Corner cases
        # N.B.  We already tested for the matrix  to be square, hence we do not
        # need to test for 0 x n or m x 0 matrices.
        #
        if self._ncols == 0:
            return self.copy()

        # Extract parameters
        #
        n  = self._ncols
        R  = self._base_ring
        MS = self._parent

        f = self.charpoly()

        # Let A denote the adjoint of M, which we want to compute, and
        # N denote a copy of M used to store powers of M.
        #
        A = f[1] * MS.identity_matrix()
        N = R(1) * MS.identity_matrix()
        for i in range(1, n):
            # Set N to be M^i
            #
            N = N * self
            A = A + f[i+1] * N
        if not (n % 2):
            A = - A

        return A

    def QR(self, full=True):
        r"""
        Returns a factorization of ``self`` as a unitary matrix
        and an upper-triangular matrix.

        INPUT:

        - ``full`` - default: ``True`` - if ``True`` then the
          returned matrices have dimensions as described below.
          If ``False`` the ``R`` matrix has no zero rows and the
          columns of ``Q`` are a basis for the column space of
          ``self``.

        OUTPUT:

        If ``self`` is an `m\times n` matrix and ``full=True`` then this
        method returns a pair of matrices: `Q` is an `m\times m` unitary
        matrix (meaning its inverse is its conjugate-transpose) and `R`
        is an `m\times n` upper-triangular matrix with non-negative entries
        on the diagonal.  For a matrix of full rank this factorization is
        unique (due to the restriction to positive entries on the diagonal).

        If ``full=False`` then `Q` has `m` rows and the columns form an
        orthonormal basis for the column space of ``self``.  So, in particular,
        the conjugate-transpose of `Q` times `Q` will be an identity matrix.
        The matrix `R` will still be upper-triangular but will also have full
        rank, in particular it will lack the zero rows present in a full
        factorization of a rank-deficient matrix.

        The results obtained when ``full=True`` are cached,
        hence `Q` and `R` are immutable matrices in this case.

        .. note::

            This is an exact computation, so limited to exact rings.
            Also the base ring needs to have a fraction field implemented
            in Sage and this field must contain square roots.  One example
            is the field of algebraic numbers, ``QQbar``, as used in the
            examples below.  If you need numerical results, convert the
            base ring to the field of complex double numbers, ``CDF``,
            which will use a faster routine that is careful about
            numerical subtleties.

        ALGORITHM:

        "Modified Gram-Schmidt," Algorithm 8.1 of [TREFETHEN-BAU]_.

        EXAMPLES:

        For a nonsingular matrix, the QR decomposition is unique. ::

            sage: A = matrix(QQbar, [[-2, 0, -4, -1, -1],
            ...                      [-2, 1, -6, -3, -1],
            ...                      [1, 1, 7, 4, 5],
            ...                      [3, 0, 8, 3, 3],
            ...                      [-1, 1, -6, -6, 5]])
            sage: Q, R = A.QR()
            sage: Q
            [ -0.4588314677411235?  -0.1260506983326509?   0.3812120831224489?   -0.394573711338418?      -0.687440062597?]
            [ -0.4588314677411235?   0.4726901187474409? -0.05198346588033394?    0.717294125164660?      -0.220962877263?]
            [  0.2294157338705618?   0.6617661662464172?   0.6619227988762521?   -0.180872093737548?      0.1964114464561?]
            [  0.6882472016116853?   0.1890760474989764?  -0.2044682991293135?    0.096630296654307?      -0.662888631790?]
            [ -0.2294157338705618?   0.5357154679137663?   -0.609939332995919?   -0.536422031427112?      0.0245514308070?]
            sage: R
            [  4.358898943540674? -0.4588314677411235?   13.07669683062202?   6.194224814505168?   2.982404540317303?]
            [                   0   1.670171752907625?  0.5987408170800917?  -1.292019657909672?   6.207996892883057?]
            [                   0                    0   5.444401659866974?   5.468660610611130?  -0.682716185228386?]
            [                   0                    0                    0   1.027626039419836?   -3.61930014968662?]
            [                   0                    0                    0                    0    0.02455143080702?]
            sage: Q.conjugate_transpose()*Q
            [1.000000000000000?            0.?e-18            0.?e-17            0.?e-15            0.?e-12]
            [           0.?e-18 1.000000000000000?            0.?e-16            0.?e-15            0.?e-12]
            [           0.?e-17            0.?e-16 1.000000000000000?            0.?e-15            0.?e-12]
            [           0.?e-15            0.?e-15            0.?e-15 1.000000000000000?            0.?e-12]
            [           0.?e-12            0.?e-12            0.?e-12            0.?e-12    1.000000000000?]
            sage: Q*R == A
            True


        An example with complex numbers in ``QQbar``, the field of algebraic
        numbers. ::

            sage: A = matrix(QQbar, [[-8, 4*I + 1, -I + 2, 2*I + 1],
            ...                      [1, -2*I - 1, -I + 3, -I + 1],
            ...                      [I + 7, 2*I + 1, -2*I + 7, -I + 1],
            ...                      [I + 2, 0, I + 12, -1]])
            sage: Q, R = A.QR()
            sage: Q
            [                          -0.7302967433402215?    0.2070566455055649? + 0.5383472783144687?*I    0.2463049809998642? - 0.0764456358723292?*I    0.2381617683194332? - 0.1036596032779695?*I]
            [                           0.0912870929175277?   -0.2070566455055649? - 0.3778783780476559?*I    0.3786559533863032? - 0.1952221495524667?*I      0.701244450214469? - 0.364371165098660?*I]
            [   0.6390096504226938? + 0.0912870929175277?*I    0.1708217325420910? + 0.6677576817554466?*I -0.03411475806452072? + 0.04090198741767143?*I    0.3140171085506763? - 0.0825191718705412?*I]
            [   0.1825741858350554? + 0.0912870929175277?*I  -0.03623491296347385? + 0.0724698259269477?*I   0.8632284069415110? + 0.06322839976356195?*I   -0.4499694867611521? - 0.0116119181208918?*I]
            sage: R
            [                          10.95445115010333?               0.?e-18 - 1.917028951268082?*I    5.385938482134133? - 2.190890230020665?*I  -0.2738612787525831? - 2.190890230020665?*I]
            [                                           0               4.829596256417300? + 0.?e-17*I   -0.869637911123373? - 5.864879483945125?*I   0.993871898426712? - 0.3054085521207082?*I]
            [                                           0                                            0               12.00160760935814? + 0.?e-16*I -0.2709533402297273? + 0.4420629644486323?*I]
            [                                           0                                            0                                            0               1.942963944258992? + 0.?e-16*I]
            sage: Q.conjugate_transpose()*Q
            [1.000000000000000? + 0.?e-19*I            0.?e-18 + 0.?e-17*I            0.?e-17 + 0.?e-17*I            0.?e-16 + 0.?e-16*I]
            [           0.?e-18 + 0.?e-17*I 1.000000000000000? + 0.?e-17*I            0.?e-17 + 0.?e-17*I            0.?e-16 + 0.?e-16*I]
            [           0.?e-17 + 0.?e-17*I            0.?e-17 + 0.?e-17*I 1.000000000000000? + 0.?e-16*I            0.?e-16 + 0.?e-16*I]
            [           0.?e-16 + 0.?e-16*I            0.?e-16 + 0.?e-16*I            0.?e-16 + 0.?e-16*I 1.000000000000000? + 0.?e-15*I]
            sage: Q*R - A
            [            0.?e-17 0.?e-17 + 0.?e-17*I 0.?e-16 + 0.?e-16*I 0.?e-16 + 0.?e-16*I]
            [            0.?e-18 0.?e-17 + 0.?e-17*I 0.?e-16 + 0.?e-16*I 0.?e-15 + 0.?e-15*I]
            [0.?e-17 + 0.?e-18*I 0.?e-17 + 0.?e-17*I 0.?e-16 + 0.?e-16*I 0.?e-16 + 0.?e-16*I]
            [0.?e-18 + 0.?e-18*I 0.?e-18 + 0.?e-17*I 0.?e-16 + 0.?e-16*I 0.?e-15 + 0.?e-16*I]

        A rank-deficient rectangular matrix, with both values of the ``full`` keyword.  ::

            sage: A = matrix(QQbar, [[2, -3, 3],
            ...                      [-1, 1, -1],
            ...                      [-1, 3, -3],
            ...                      [-5, 1, -1]])
            sage: Q, R = A.QR()
            sage: Q
            [  0.3592106040535498?  -0.5693261797050169?   0.7239227659930268?   0.1509015305256380?]
            [ -0.1796053020267749?   0.1445907757980996?                     0   0.9730546968377341?]
            [ -0.1796053020267749?   0.7048800320157352?    0.672213996993525?  -0.1378927778941174?]
            [ -0.8980265101338745?  -0.3976246334447737?   0.1551263069985058? -0.10667177157846818?]
            sage: R
            [ 5.567764362830022? -2.694079530401624?  2.694079530401624?]
            [                  0  3.569584777515583? -3.569584777515583?]
            [                  0                   0                   0]
            [                  0                   0                   0]
            sage: Q.conjugate_transpose()*Q
            [                 1            0.?e-18            0.?e-18            0.?e-18]
            [           0.?e-18                  1            0.?e-18            0.?e-18]
            [           0.?e-18            0.?e-18 1.000000000000000?            0.?e-18]
            [           0.?e-18            0.?e-18            0.?e-18 1.000000000000000?]

            sage: Q, R = A.QR(full=False)
            sage: Q
            [ 0.3592106040535498? -0.5693261797050169?]
            [-0.1796053020267749?  0.1445907757980996?]
            [-0.1796053020267749?  0.7048800320157352?]
            [-0.8980265101338745? -0.3976246334447737?]
            sage: R
            [ 5.567764362830022? -2.694079530401624?  2.694079530401624?]
            [                  0  3.569584777515583? -3.569584777515583?]
            sage: Q.conjugate_transpose()*Q
            [      1 0.?e-18]
            [0.?e-18       1]

        Another rank-deficient rectangular matrix, with complex entries,
        as a reduced decomposition. ::

            sage: A = matrix(QQbar, [[-3*I - 3, I - 3, -12*I + 1, -2],
            ...                      [-I - 1, -2, 5*I - 1, -I - 2],
            ...                      [-4*I - 4, I - 5, -7*I, -I - 4]])
            sage: Q, R = A.QR(full=False)
            sage: Q
            [ -0.4160251471689219? - 0.4160251471689219?*I   0.5370861555295747? + 0.1790287185098583?*I]
            [ -0.1386750490563073? - 0.1386750490563073?*I  -0.7519206177414046? - 0.2506402059138015?*I]
            [ -0.5547001962252291? - 0.5547001962252291?*I -0.2148344622118299? - 0.07161148740394329?*I]
            sage: R
            [                        7.211102550927979?  3.328201177351375? - 5.269651864139676?*I   7.904477796209515? + 8.45917799243475?*I  4.021576422632911? - 2.634825932069838?*I]
            [                                         0                         1.074172311059150?  -1.611258466588724? - 9.13046464400277?*I 1.611258466588724? + 0.5370861555295747?*I]
            sage: Q.conjugate_transpose()*Q
            [1.000000000000000? + 0.?e-18*I            0.?e-18 + 0.?e-18*I]
            [           0.?e-17 + 0.?e-17*I 1.000000000000000? + 0.?e-17*I]
            sage: Q*R-A
            [0.?e-18 + 0.?e-18*I 0.?e-18 + 0.?e-18*I 0.?e-17 + 0.?e-17*I 0.?e-18 + 0.?e-18*I]
            [0.?e-18 + 0.?e-18*I 0.?e-18 + 0.?e-18*I 0.?e-18 + 0.?e-17*I 0.?e-18 + 0.?e-18*I]
            [0.?e-18 + 0.?e-18*I 0.?e-17 + 0.?e-18*I 0.?e-17 + 0.?e-17*I 0.?e-18 + 0.?e-18*I]

        Results of full decompositions are cached and thus returned
        immutable.  ::

            sage: A = random_matrix(QQbar, 2, 2)
            sage: Q, R = A.QR()
            sage: Q.is_mutable()
            False
            sage: R.is_mutable()
            False

        Trivial cases return trivial results of the correct size,
        and we check `Q` itself in one case.  ::

            sage: A = zero_matrix(QQbar, 0, 10)
            sage: Q, R = A.QR()
            sage: Q.nrows(), Q.ncols()
            (0, 0)
            sage: R.nrows(), R.ncols()
            (0, 10)
            sage: A = zero_matrix(QQbar, 3, 0)
            sage: Q, R = A.QR()
            sage: Q.nrows(), Q.ncols()
            (3, 3)
            sage: R.nrows(), R.ncols()
            (3, 0)
            sage: Q
            [1 0 0]
            [0 1 0]
            [0 0 1]

        TESTS:

        Inexact rings are caught and ``CDF`` suggested.  ::

            sage: A = matrix(RealField(100), 2, range(4))
            sage: A.QR()
            Traceback (most recent call last):
            ...
            NotImplementedError: QR decomposition is implemented over exact rings, try CDF for numerical results, not Real Field with 100 bits of precision

        Without a fraction field, we cannot hope to run the algorithm. ::

            sage: A = matrix(Integers(6), 2, range(4))
            sage: A.QR()
            Traceback (most recent call last):
            ...
            ValueError: QR decomposition needs a fraction field of Ring of integers modulo 6

        The biggest obstacle is making unit vectors, thus requiring square
        roots, though some small cases pass through.  ::

            sage: A = matrix(ZZ, 3, range(9))
            sage: A.QR()
            Traceback (most recent call last):
            ...
            TypeError: QR decomposition unable to compute square roots in Rational Field

            sage: A = matrix(ZZ, 2, range(4))
            sage: Q, R = A.QR()
            sage: Q
            [0 1]
            [1 0]
            sage: R
            [2 3]
            [0 1]

        REFERENCES:

        .. [TREFETHEN-BAU] Trefethen, Lloyd N., Bau, David, III
           "Numerical Linear Algebra"
           SIAM, Philadelphia, 1997.

        AUTHOR:

        - Rob Beezer (2011-02-17)
        """
        from sage.modules.free_module_element import zero_vector
        from sage.matrix.constructor import zero_matrix, matrix
        from sage.functions.other import sqrt

        if full:
            QR = self.fetch('QR_factors')
            if QR is not None:
                return QR
        R = self.base_ring()
        if not R.is_exact():
            raise NotImplementedError('QR decomposition is implemented over exact rings, try CDF for numerical results, not %s' % R)
        try:
            F = R.fraction_field()
        except StandardError:
            raise ValueError("QR decomposition needs a fraction field of %s" % R)
        m = self.nrows()
        n = self.ncols()

        R = zero_matrix(F, m, n)
        V = self.columns(copy=True)
        Q = []
        row = 0  # row of R being filled
        for i in range(n):
            v = V[i]
            hip = v.hermitian_inner_product(v)
            if hip != 0:
                try:
                    scale = sqrt(hip)
                    q = (1/scale)*v
                    Q.append(q)
                    R[row,i] = scale
                    for j in range(i+1, n):
                        R[row,j] = q.hermitian_inner_product(V[j])
                        V[j] = V[j] - R[row,j]*q
                    row = row + 1
                except TypeError:
                    raise TypeError('QR decomposition unable to compute square roots in %s' % F)
        # complete to full orthonormal basis, or reduce to truncated R
        if full:
            Qt = matrix(Q) # as rows here
            if Qt.nrows() == 0:
                Qt = zero_matrix(F, 0, m)
            orthogonal = Qt.right_kernel().basis_matrix().transpose()
            Qperp, _ = orthogonal.QR(full=False)
            Q = Q + Qperp.columns()
        else:
            R = R[0:len(Q), 0:n]
        Q = matrix(Q).transpose()
        # Adjust rows of Q if empty
        if Q.ncols() == 0:
            Q = zero_matrix(F, m, 0)
        QR = (Q, R)
        if full:
            Q.set_immutable()
            R.set_immutable()
            self.cache('QR_factors', QR)
        return QR

    def _gram_schmidt_noscale(self):
        r"""
        Performs Gram-Schmidt orthogonalization, with no scaling to unit vectors.

        INPUT:

        - ``self`` - is a matrix whose columns are to be orthogonalized.
          The base ring of the matrix needs to have its fraction field
          implemented.

        OUTPUT:

        Two matrices, ``Q`` and ``R`` such that if ``A`` represents ``self``:

        - ``A = Q*R``
        - The columns of ``Q`` are an orthogonal set which spans the
          column space of ``A``.
        - The conjugate-transpose of ``Q`` times ``Q`` is a diagonal matrix.
        - ``R`` is a full-rank matrix, that has all entries below the
          main diagonal equal to zero.

        This is basically a "reduced" QR decomposition of ``self`` with
        the distinction that the orthogonal column vectors of ``Q`` have
        not been scaled to unit vectors, avoiding the need to take square
        roots.

        EXAMPLES:

        A rectangular matrix whose columns have full-rank.  Notice that the
        routine computes in the fraction field, requiring the column space
        check to step up to ``QQ``.  ::

            sage: A = matrix(ZZ, [[-1, -3,  0, -1],
            ...                   [ 1,  2, -1,  2],
            ...                   [-3, -6,  4, -7]])
            sage: Q,R = A._gram_schmidt_noscale()
            sage: Q
            [    -1 -10/11      0]
            [     1  -1/11   3/10]
            [    -3   3/11   1/10]
            sage: R
            [     1  23/11 -13/11  24/11]
            [     0      1  13/10 -13/10]
            [     0      0      1     -1]
            sage: Q*R == A
            True
            sage: Q.transpose()*Q
            [   11     0     0]
            [    0 10/11     0]
            [    0     0  1/10]
            sage: A.change_ring(QQ).column_space() == Q.column_space()
            True

        A matrix over a subfield of the complex numbers, with four
        columns but rank 3, so the orthogonal set has just 3 vectors
        as well.  Orthogonality comes from the Hermitian inner product
        so we need to check with the conjugate-transpose.  This
        example verifies that the bug on #10791 is fixed.  ::

            sage: F.<a> = QuadraticField(-5)
            sage: A = matrix(F, [[    1,   a - 3,   a - 2, a + 1],
            ...                  [    a, 2*a + 1, 3*a + 1,     1],
            ...                  [a + 1,   a - 6, 2*a - 5,     1],
            ...                  [  2*a,       a,     3*a,    -3],
            ...                  [    1,   a - 1,       a, a - 1]])
            sage: A.rank()
            3
            sage: Q, R = A._gram_schmidt_noscale()
            sage: Q
            [                      1         25/33*a - 38/11   641/1163*a + 453/1163]
            [                      a         17/11*a + 73/33  322/1163*a + 1566/1163]
            [                  a + 1        10/33*a - 173/33 -784/1163*a + 1614/1163]
            [                    2*a          1/11*a + 80/33  196/1163*a - 1234/1163]
            [                      1         25/33*a - 16/11  855/1163*a - 1717/1163]
            sage: R
            [                    1         8/33*a + 5/11        8/33*a + 16/11         2/11*a + 1/33]
            [                    0                     1                     1 -107/1163*a - 78/1163]
            [                    0                     0                     0                     1]
            sage: Q*R == A
            True
            sage: Q.transpose().conjugate()*Q
            [        33          0          0]
            [         0    2326/33          0]
            [         0          0 16532/1163]
            sage: Q.column_space() == A.column_space()
            True

        Some trivial cases.  ::

            sage: A = matrix(ZZ, 3, 0)
            sage: Q, R = A._gram_schmidt_noscale()
            sage: Q.parent()
            Full MatrixSpace of 3 by 0 dense matrices over Rational Field
            sage: R.parent()
            Full MatrixSpace of 0 by 0 dense matrices over Rational Field
            sage: Q*R == A
            True

            sage: A = matrix(ZZ, 0, 3)
            sage: Q, R = A._gram_schmidt_noscale()
            sage: Q.parent()
            Full MatrixSpace of 0 by 0 dense matrices over Rational Field
            sage: R.parent()
            Full MatrixSpace of 0 by 3 dense matrices over Rational Field
            sage: Q*R == A
            True

        TESTS:

        Without a fraction field, we cannot hope to proceed. ::

            sage: A = matrix(Integers(6), 2, range(4))
            sage: A._gram_schmidt_noscale()
            Traceback (most recent call last):
            ...
            TypeError: Gram-Schmidt orthogonalization requires a base ring with a fraction field, not Ring of integers modulo 6

        AUTHORS:

        - William Stein (2007-11-18)
        - Rob Beezer (2011-02-25)
        """
        import sage.matrix.constructor
        R = self.base_ring()
        try:
            F = R.fraction_field()
        except TypeError:
            raise TypeError("Gram-Schmidt orthogonalization requires a base ring with a fraction field, not %s" % R)
        n = self.ncols()
        B = self.columns()
        zero = F(0)
        Bstar = []
        R = sage.matrix.constructor.zero_matrix(F, n)
        nnz = 0 # number non-zero rows in R, or number of nonzero vectors in Bstar
        for i in range(n):
            ortho = B[i]
            for j in range(nnz):
                R[j,i] = Bstar[j].hermitian_inner_product(B[i])/Bstar[j].hermitian_inner_product(Bstar[j])
                ortho = ortho - R[j,i]*Bstar[j]
            if ortho.hermitian_inner_product(ortho) != zero:
                Bstar.append(ortho)
                R[nnz, i] = 1
                nnz = nnz + 1
        R = R[0:nnz]
        if Bstar == []:
            Q = sage.matrix.constructor.matrix(F, 0, self.nrows()).transpose()
        else:
            Q = sage.matrix.constructor.matrix(F, Bstar).transpose()
        return Q, R

    def gram_schmidt(self, orthonormal=False):
        r"""
        Performs Gram-Schmidt orthogonalization on the rows of the matrix,
        returning a new matrix and a matrix accomplishing the transformation.

        INPUT:

        - ``self`` - a matrix whose rows are to be orthogonalized.
        - ``orthonormal`` - default: ``False`` - if ``True`` the
          returned orthogonal vectors are unit vectors.  This keyword
          is ignored if the matrix is over ``RDF`` or ``CDF`` and the
          results are always orthonormal.

        OUTPUT:

        A pair of matrices, ``G`` and ``M`` such that if ``A``
        represents ``self``, where the parenthetical properties occur
        when ``orthonormal = True``:

        - ``A = M*G``
        - The rows of ``G`` are an orthogonal (resp. orthonormal)
          set of vectors.
        - ``G`` times the conjugate-transpose of ``G`` is a diagonal
          (resp. identity) matrix.
        - The row space of ``G`` equals the row space of ``A``.
        - ``M`` is a full-rank matrix with zeros above the diagonal.

        For exact rings, any zero vectors produced (when the original
        vectors are linearly dependent) are not output, thus the
        orthonormal set is linearly independent, and thus a basis for the
        row space of the original matrix.

        Any notion of a Gram-Schmidt procedure requires that the base
        ring of the matrix has a fraction field implemented.  In order
        to arrive at an orthonormal set, it must be possible to construct
        square roots of the elements of the base field.  In Sage, your
        best option is the field of algebraic numbers, ``QQbar``, which
        properly contains the rationals and number fields.

        If you have an approximate numerical matrix, then this routine
        requires that your base field be the real and complex
        double-precision floating point numbers, ``RDF`` and ``CDF``.
        In this case, the matrix is treated as having full rank, as no
        attempt is made to recognize linear dependence with approximate
        calculations.

        EXAMPLES:

        Inexact Rings, Numerical Matrices:

        First, the inexact rings, ``CDF`` and ``RDF``.  ::

            sage: A = matrix(CDF, [[ 0.6454 + 0.7491*I, -0.8662 + 0.1489*I,  0.7656 - 0.00344*I],
            ...                    [-0.2913 + 0.8057*I,  0.8321 + 0.8170*I, -0.6744 + 0.9248*I],
            ...                    [ 0.2554 + 0.3517*I, -0.4454 - 0.1715*I,  0.8325 - 0.6282*I]])
            sage: G, M = A.gram_schmidt()
            sage: G.round(6)  # random signs
            [-0.422243 - 0.490087*I  0.566698 - 0.097416*I -0.500882 + 0.002251*I]
            [-0.057002 - 0.495035*I  -0.35059 - 0.625323*I  0.255514 - 0.415284*I]
            [ 0.394105 - 0.421778*I -0.392266 - 0.039345*I  -0.352905 + 0.62195*I]
            sage: M.round(6)  # random
            [             -1.528503                    0.0                    0.0]
            [  0.459974 - 0.40061*I              -1.741233                    0.0]
            [-0.934304 + 0.148868*I   0.54833 + 0.073202*I              -0.550725]
            sage: (A - M*G).zero_at(10^-12)
            [0.0 0.0 0.0]
            [0.0 0.0 0.0]
            [0.0 0.0 0.0]
            sage: (G*G.conjugate_transpose()).zero_at(10^-12)
            [1.0 0.0 0.0]
            [0.0 1.0 0.0]
            [0.0 0.0 1.0]


        A rectangular matrix.  Note that the ``orthonormal`` keyword
        is ignored in these cases.  ::

            sage: A = matrix(RDF, [[-0.978325, -0.751994, 0.925305, -0.200512, 0.420458],
            ...                    [-0.474877, -0.983403, 0.089836,  0.132218, 0.672965]])
            sage: G, M = A.gram_schmidt(orthonormal=False)
            sage: G.round(6).zero_at(10^-6)
            [-0.607223 -0.466745  0.574315 -0.124453  0.260968]
            [ 0.123203 -0.617909 -0.530578  0.289773  0.487368]
            sage: M.round(6).zero_at(10^-6)
            [1.611147      0.0]
            [0.958116 0.867778]
            sage: (A-M*G).zero_at(10^-12)
            [0.0 0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0 0.0]
            sage: (G*G.transpose()).round(6).zero_at(10^-6)
            [1.0 0.0]
            [0.0 1.0]

        Even though a set of vectors may be linearly dependent, no effort
        is made to decide when a zero vector is really the result of a
        relation of linear dependence.  So in this regard, input matrices
        are treated as being of full rank.  Try one of the base rings that
        provide exact results if you need exact results.  ::

            sage: entries = [[1,1,2], [2,1,3], [3,1,4]]
            sage: A = matrix(QQ, entries)
            sage: A.rank()
            2
            sage: B = matrix(RDF, entries)
            sage: G, M = B.gram_schmidt()
            sage: G.round(6)  # random signs
            [-0.408248 -0.408248 -0.816497]
            [ 0.707107 -0.707107      -0.0]
            [ -0.57735  -0.57735   0.57735]
            sage: M.round(10)  # random
            [-2.4494897428           0.0           0.0]
            [-3.6742346142  0.7071067812           0.0]
            [-4.8989794856  1.4142135624           0.0]
            sage: (A - M*G).zero_at(1e-14)
            [0.0 0.0 0.0]
            [0.0 0.0 0.0]
            [0.0 0.0 0.0]
            sage: (G*G.transpose()).zero_at(1e-14)
            [1.0 0.0 0.0]
            [0.0 1.0 0.0]
            [0.0 0.0 1.0]

        Exact Rings, Orthonormalization:

        To scale a vector to unit length requires taking
        a square root, which often takes us outside the base ring.
        For the integers and the rationals, the field of algebraic numbers
        (``QQbar``) is big enough to contain what we need, but the price
        is that the computations are very slow, hence mostly of value
        for small cases or instruction. Now we need to use the
        ``orthonormal`` keyword.  ::

            sage: A = matrix(QQbar, [[6, -8,  1],
            ...                      [4,  1,  3],
            ...                      [6,  3,  3],
            ...                      [7,  1, -5],
            ...                      [7, -3,  5]])
            sage: G, M = A.gram_schmidt(orthonormal=True)
            sage: G
            [ 0.5970223141259934? -0.7960297521679913? 0.09950371902099891?]
            [ 0.6063218341690895?  0.5289635311888953?  0.5937772444966257?]
            [ 0.5252981913594170?  0.2941669871612735?  -0.798453250866314?]
            sage: M
            [ 10.04987562112089?                   0                   0]
            [ 1.890570661398980?  4.735582601355131?                   0]
            [ 1.492555785314984?  7.006153332071100?  1.638930357041381?]
            [ 2.885607851608969?  1.804330147889395?  7.963520581008761?]
            [ 7.064764050490923?  5.626248468100069? -1.197679876299471?]
            sage: M*G-A
            [0 0 0]
            [0 0 0]
            [0 0 0]
            [0 0 0]
            [0 0 0]
            sage: (G*G.transpose()-identity_matrix(3)).norm() < 10^-10
            True
            sage: G.row_space() == A.row_space()
            True

        After :trac:`14047`, the matrix can also be over the algebraic reals
        ``AA``::

            sage: A = matrix(AA, [[6, -8,  1],
            ...                   [4,  1,  3],
            ...                   [6,  3,  3],
            ...                   [7,  1, -5],
            ...                   [7, -3,  5]])
            sage: G, M = A.gram_schmidt(orthonormal=True)
            sage: G
            [ 0.5970223141259934? -0.7960297521679913? 0.09950371902099891?]
            [ 0.6063218341690895?  0.5289635311888953?  0.5937772444966257?]
            [ 0.5252981913594170?  0.2941669871612735?  -0.798453250866314?]
            sage: M
            [ 10.04987562112089?                   0                   0]
            [ 1.890570661398980?  4.735582601355131?                   0]
            [ 1.492555785314984?  7.006153332071100?  1.638930357041381?]
            [ 2.885607851608969?  1.804330147889395?  7.963520581008761?]
            [ 7.064764050490923?  5.626248468100069? -1.197679876299471?]

        Starting with complex numbers with rational real and imaginary parts.
        Note the use of the conjugate-transpose when checking the
        orthonormality. ::

            sage: A = matrix(QQbar, [[  -2,    -I - 1, 4*I + 2,       -1],
            ...                      [-4*I, -2*I + 17,       0,  9*I + 1],
            ...                      [   1,  -2*I - 6, -I + 11, -5*I + 1]])
            sage: G, M = A.gram_schmidt(orthonormal=True)
            sage: (M*G-A).norm() < 10^-10
            True
            sage: id3 = G*G.conjugate().transpose()
            sage: (id3 - identity_matrix(3)).norm() < 10^-10
            True
            sage: G.row_space() == A.row_space()  # long time
            True

        A square matrix with small rank.  The zero vectors produced as a
        result of linear dependence get eliminated, so the rows of ``G``
        are a basis for the row space of ``A``.  ::

            sage: A = matrix(QQbar, [[2, -6, 3, 8],
            ...                      [1, -3, 2, 5],
            ...                      [0,  0, 2, 4],
            ...                      [2, -6, 3, 8]])
            sage: A.change_ring(QQ).rank()
            2
            sage: G, M = A.gram_schmidt(orthonormal=True)
            sage: G
            [ 0.1881441736767195? -0.5644325210301583?  0.2822162605150792?  0.7525766947068779?]
            [-0.2502818123591464?   0.750845437077439?  0.3688363550555841?  0.4873908977520218?]
            sage: M
            [10.630145812734649?                   0]
            [ 6.208757731331742? 0.6718090752798139?]
            [ 3.574739299857670?  2.687236301119256?]
            [10.630145812734649?                   0]
            sage: M*G-A
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]
            sage: (G*G.transpose()-identity_matrix(2)).norm() < 10^-10
            True
            sage: G.row_space() == A.row_space()
            True

        Exact Rings, Orthogonalization:

        If we forego scaling orthogonal vectors to unit vectors, we
        can apply Gram-Schmidt to a much greater variety of rings.
        Use the ``orthonormal=False`` keyword (or assume it as the default).
        Note that now the orthogonality check creates a diagonal matrix
        whose diagonal entries are the squares of the lengths of the
        vectors.

        First, in the rationals, without involving ``QQbar``.  ::

            sage: A = matrix(QQ, [[-1,  3,  2,  2],
            ...                   [-1,  0, -1,  0],
            ...                   [-1, -2, -3, -1],
            ...                   [ 1,  1,  2,  0]])
            sage: A.rank()
            3
            sage: G, M = A.gram_schmidt()
            sage: G
            [    -1      3      2      2]
            [-19/18    1/6   -8/9    1/9]
            [  2/35  -4/35  -2/35   9/35]
            sage: M
            [     1      0      0]
            [ -1/18      1      0]
            [-13/18  59/35      1]
            [   1/3 -48/35     -2]
            sage: M*G-A
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]
            sage: G*G.transpose()
            [   18     0     0]
            [    0 35/18     0]
            [    0     0  3/35]
            sage: G.row_space() == A.row_space()
            True

        A complex subfield of the complex numbers.  ::

            sage: C.<z> = CyclotomicField(5)
            sage: A = matrix(C, [[              -z^3 - 2*z,             -z^3 - 1, 2*z^3 - 2*z^2 + 2*z,             1],
            ...                  [         z^3 - 2*z^2 + 1, -z^3 + 2*z^2 - z - 1,                  -1,       z^2 + z],
            ...                  [-1/2*z^3 - 2*z^2 + z + 1,         -z^3 + z - 2,    -2*z^3 + 1/2*z^2, 2*z^2 - z + 2]])
            sage: G, M = A.gram_schmidt(orthonormal=False)
            sage: G
            [                                                      -z^3 - 2*z                                                         -z^3 - 1                                              2*z^3 - 2*z^2 + 2*z                                                                1]
            [                   155/139*z^3 - 161/139*z^2 + 31/139*z + 13/139                 -175/139*z^3 + 180/139*z^2 - 125/139*z - 142/139                     230/139*z^3 + 124/139*z^2 + 6/139*z + 19/139                      -14/139*z^3 + 92/139*z^2 - 6/139*z - 95/139]
            [-10359/19841*z^3 - 36739/39682*z^2 + 24961/39682*z - 11879/39682  -28209/39682*z^3 - 3671/19841*z^2 + 51549/39682*z - 38613/39682    -42769/39682*z^3 - 615/39682*z^2 - 1252/19841*z - 14392/19841   4895/19841*z^3 + 57885/39682*z^2 - 46094/19841*z + 65747/39682]
            sage: M
            [                                                           1                                                            0                                                            0]
            [                14/139*z^3 + 47/139*z^2 + 145/139*z + 95/139                                                            1                                                            0]
            [              -7/278*z^3 + 199/278*z^2 + 183/139*z + 175/278 -3785/39682*z^3 + 3346/19841*z^2 - 3990/19841*z + 2039/19841                                                            1]
            sage: M*G - A
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]
            sage: G*G.conjugate().transpose()
            [                               15*z^3 + 15*z^2 + 28                                                   0                                                   0]
            [                                                  0                463/139*z^3 + 463/139*z^2 + 1971/139                                                   0]
            [                                                  0                                                   0 230983/19841*z^3 + 230983/19841*z^2 + 1003433/39682]
            sage: G.row_space() == A.row_space()
            True

        A slightly edited legacy example.  ::

            sage: A = matrix(ZZ, 3, [-1, 2, 5, -11, 1, 1, 1, -1, -3]); A
            [ -1   2   5]
            [-11   1   1]
            [  1  -1  -3]
            sage: G, mu = A.gram_schmidt()
            sage: G
            [     -1       2       5]
            [  -52/5    -1/5      -2]
            [  2/187  36/187 -14/187]
            sage: mu
            [     1      0      0]
            [   3/5      1      0]
            [  -3/5 -7/187      1]
            sage: G.row(0) * G.row(1)
            0
            sage: G.row(0) * G.row(2)
            0
            sage: G.row(1) * G.row(2)
            0

        The relation between mu and A is as follows.  ::

            sage: mu*G == A
            True
        """
        import sage.rings.real_double
        import sage.rings.complex_double
        R = self.base_ring()
        if R in [sage.rings.real_double.RDF, sage.rings.complex_double.CDF]:
            Q, R = self.transpose().QR()
            m = R.nrows(); n = R.ncols()
            if m > n:
                Q = Q[0:m, 0:n]
                R = R[0:n, 0:n]
        elif R.is_exact():
            if orthonormal:
                Q, R = self.transpose().QR(full=False)
            else:
                Q, R = self.transpose()._gram_schmidt_noscale()
        else:
            raise NotImplementedError("Gram-Schmidt orthogonalization not implemented for matrices over inexact rings, except for RDF and CDF")
        return Q.transpose(), R.transpose()

    def jordan_form(self, base_ring=None, sparse=False, subdivide=True, transformation=False, eigenvalues=None, check_input=True):
        r"""
        Compute the Jordan normal form of this square matrix `A`, if it exists.

        This computation is performed in a naive way using the ranks of powers
        of `A-xI`, where `x` is an eigenvalue of the matrix `A`.  If desired,
        a transformation matrix `P` can be returned, which is such that the
        Jordan canonical form is given by `P^{-1} A P`.

        INPUT:

        - ``base_ring`` - Ring in which to compute the Jordan form.

        - ``sparse`` - (default ``False``) If ``sparse=True``, return a sparse
          matrix.

        - ``subdivide`` - (default ``True``) If ``subdivide=True``, the
          subdivisions for the Jordan blocks in the matrix are shown.

        - ``transformation`` - (default ``False``) If ``transformation=True``,
          computes also the transformation matrix.

        - ``eigenvalues`` - (default ``None``) A complete set of roots, with
          multiplicity, of the characteristic polynomial of `A`, encoded as
          a list of pairs, each having the form `(r, m)` with `r` a root and
          `m` its multiplicity. If this is ``None``, then Sage computes this
          list itself, but this is only possible over base rings in whose
          quotient fields polynomial factorization is implemented. Over all
          other rings, providing this list manually is the only way to
          compute Jordan normal forms.

        - ``check_input`` - (default ``True``) A Boolean specifying whether
          the list ``eigenvalues`` (if provided) has to be checked for
          correctness. Set this to ``False`` for a speedup if the eigenvalues
          are known to be correct.

        NOTES:

        Currently, the Jordan normal form is not computed over inexact rings
        in any but the trivial cases when the matrix is either `0 \times 0`
        or `1 \times 1`.

        In the case of exact rings, this method does not compute any
        generalized form of the Jordan normal form, but is only able to
        compute the result if the characteristic polynomial of the matrix
        splits over the specific base ring.

        Note that the base ring must be a field or a ring with an implemented
        fraction field.

        EXAMPLES::

            sage: a = matrix(ZZ,4,[1, 0, 0, 0, 0, 1, 0, 0, 1, \
            -1, 1, 0, 1, -1, 1, 2]); a
            [ 1  0  0  0]
            [ 0  1  0  0]
            [ 1 -1  1  0]
            [ 1 -1  1  2]
            sage: a.jordan_form()
            [2|0 0|0]
            [-+---+-]
            [0|1 1|0]
            [0|0 1|0]
            [-+---+-]
            [0|0 0|1]
            sage: a.jordan_form(subdivide=False)
            [2 0 0 0]
            [0 1 1 0]
            [0 0 1 0]
            [0 0 0 1]
            sage: b = matrix(ZZ,3,range(9)); b
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: b.jordan_form()
            Traceback (most recent call last):
            ...
            RuntimeError: Some eigenvalue does not exist in Rational Field.
            sage: b.jordan_form(RealField(15))
            Traceback (most recent call last):
            ...
            ValueError: Jordan normal form not implemented over inexact rings.

        If you need the transformation matrix as well as the Jordan form of
        ``self``, then pass the option ``transformation=True``.

        ::

            sage: m = matrix([[5,4,2,1],[0,1,-1,-1],[-1,-1,3,0],[1,1,-1,2]]); m
            [ 5  4  2  1]
            [ 0  1 -1 -1]
            [-1 -1  3  0]
            [ 1  1 -1  2]
            sage: jf, p = m.jordan_form(transformation=True)
            sage: jf
            [2|0|0 0]
            [-+-+---]
            [0|1|0 0]
            [-+-+---]
            [0|0|4 1]
            [0|0|0 4]
            sage: ~p * m * p
            [2 0 0 0]
            [0 1 0 0]
            [0 0 4 1]
            [0 0 0 4]

        Note that for matrices over inexact rings and associated numerical
        stability problems, we do not attempt to compute the Jordan normal
        form.

        ::

            sage: b = matrix(ZZ,3,3,range(9))
            sage: jf, p = b.jordan_form(RealField(15), transformation=True)
            Traceback (most recent call last):
            ...
            ValueError: Jordan normal form not implemented over inexact rings.

        TESTS::

            sage: c = matrix(ZZ, 3, [1]*9); c
            [1 1 1]
            [1 1 1]
            [1 1 1]
            sage: c.jordan_form(subdivide=False)
            [3 0 0]
            [0 0 0]
            [0 0 0]

        ::

            sage: evals = [(i,i) for i in range(1,6)]
            sage: n = sum(range(1,6))
            sage: jf = block_diagonal_matrix([jordan_block(ev,size) for ev,size in evals])
            sage: p = random_matrix(ZZ,n,n)
            sage: while p.rank() != n: p = random_matrix(ZZ,n,n)
            sage: m = p * jf * ~p
            sage: mjf, mp = m.jordan_form(transformation=True)
            sage: mjf == jf
            True
            sage: m = diagonal_matrix([1,1,0,0])
            sage: jf,P = m.jordan_form(transformation=True)
            sage: jf == ~P*m*P
            True

        We verify that the bug from trac ticket #6942 is fixed::

            sage: M = Matrix(GF(2),[[1,0,1,0,0,0,1],[1,0,0,1,1,1,0],[1,1,0,1,1,1,1],[1,1,1,0,1,1,1],[1,1,1,0,0,1,0],[1,1,1,0,1,0,0],[1,1,1,1,1,1,0]])
            sage: J, T = M.jordan_form(transformation=True)
            sage: J
            [1 1|0 0|0 0|0]
            [0 1|0 0|0 0|0]
            [---+---+---+-]
            [0 0|1 1|0 0|0]
            [0 0|0 1|0 0|0]
            [---+---+---+-]
            [0 0|0 0|1 1|0]
            [0 0|0 0|0 1|0]
            [---+---+---+-]
            [0 0|0 0|0 0|1]
            sage: M * T == T * J
            True
            sage: T.rank()
            7
            sage: M.rank()
            7

        We verify that the bug from trac ticket #6932 is fixed::

            sage: M=Matrix(1,1,[1])
            sage: M.jordan_form(transformation=True)
            ([1], [1])

        We now go through three `10 \times 10` matrices to exhibit cases where
        there are multiple blocks of the same size::

            sage: A = matrix(QQ, [[15, 37/3, -16, -104/3, -29, -7/3, 0, 2/3, -29/3, -1/3], [2, 9, -1, -6, -6, 0, 0, 0, -2, 0], [24, 74/3, -41, -208/3, -58, -23/3, 0, 4/3, -58/3, -2/3], [-6, -19, 3, 21, 19, 0, 0, 0, 6, 0], [2, 6, 3, -6, -3, 1, 0, 0, -2, 0], [-96, -296/3, 176, 832/3, 232, 101/3, 0, -16/3, 232/3, 8/3], [-4, -2/3, 21, 16/3, 4, 14/3, 3, -1/3, 4/3, -25/3], [20, 26/3, -66, -199/3, -42, -41/3, 0, 13/3, -55/3, -2/3], [18, 57, -9, -54, -57, 0, 0, 0, -15, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 3]]); A
            [    15   37/3    -16 -104/3    -29   -7/3      0    2/3  -29/3   -1/3]
            [     2      9     -1     -6     -6      0      0      0     -2      0]
            [    24   74/3    -41 -208/3    -58  -23/3      0    4/3  -58/3   -2/3]
            [    -6    -19      3     21     19      0      0      0      6      0]
            [     2      6      3     -6     -3      1      0      0     -2      0]
            [   -96 -296/3    176  832/3    232  101/3      0  -16/3  232/3    8/3]
            [    -4   -2/3     21   16/3      4   14/3      3   -1/3    4/3  -25/3]
            [    20   26/3    -66 -199/3    -42  -41/3      0   13/3  -55/3   -2/3]
            [    18     57     -9    -54    -57      0      0      0    -15      0]
            [     0      0      0      0      0      0      0      0      0      3]
            sage: J, T = A.jordan_form(transformation=True); J
            [3 1 0|0 0 0|0 0 0|0]
            [0 3 1|0 0 0|0 0 0|0]
            [0 0 3|0 0 0|0 0 0|0]
            [-----+-----+-----+-]
            [0 0 0|3 1 0|0 0 0|0]
            [0 0 0|0 3 1|0 0 0|0]
            [0 0 0|0 0 3|0 0 0|0]
            [-----+-----+-----+-]
            [0 0 0|0 0 0|3 1 0|0]
            [0 0 0|0 0 0|0 3 1|0]
            [0 0 0|0 0 0|0 0 3|0]
            [-----+-----+-----+-]
            [0 0 0|0 0 0|0 0 0|3]
            sage: T * J * T**(-1) == A
            True
            sage: T.rank()
            10

        ::

            sage: A = matrix(QQ, [[15, 37/3, -16, -14/3, -29, -7/3, 0, 2/3, 1/3, 44/3], [2, 9, -1, 0, -6, 0, 0, 0, 0, 3], [24, 74/3, -41, -28/3, -58, -23/3, 0, 4/3, 2/3, 88/3], [-6, -19, 3, 3, 19, 0, 0, 0, 0, -9], [2, 6, 3, 0, -3, 1, 0, 0, 0, 3], [-96, -296/3, 176, 112/3, 232, 101/3, 0, -16/3, -8/3, -352/3], [-4, -2/3, 21, 16/3, 4, 14/3, 3, -1/3, 4/3, -25/3], [20, 26/3, -66, -28/3, -42, -41/3, 0, 13/3, 2/3, 82/3], [18, 57, -9, 0, -57, 0, 0, 0, 3, 28], [0, 0, 0, 0, 0, 0, 0, 0, 0, 3]]); A
            [    15   37/3    -16  -14/3    -29   -7/3      0    2/3    1/3   44/3]
            [     2      9     -1      0     -6      0      0      0      0      3]
            [    24   74/3    -41  -28/3    -58  -23/3      0    4/3    2/3   88/3]
            [    -6    -19      3      3     19      0      0      0      0     -9]
            [     2      6      3      0     -3      1      0      0      0      3]
            [   -96 -296/3    176  112/3    232  101/3      0  -16/3   -8/3 -352/3]
            [    -4   -2/3     21   16/3      4   14/3      3   -1/3    4/3  -25/3]
            [    20   26/3    -66  -28/3    -42  -41/3      0   13/3    2/3   82/3]
            [    18     57     -9      0    -57      0      0      0      3     28]
            [     0      0      0      0      0      0      0      0      0      3]
            sage: J, T = A.jordan_form(transformation=True); J
            [3 1 0|0 0 0|0 0|0 0]
            [0 3 1|0 0 0|0 0|0 0]
            [0 0 3|0 0 0|0 0|0 0]
            [-----+-----+---+---]
            [0 0 0|3 1 0|0 0|0 0]
            [0 0 0|0 3 1|0 0|0 0]
            [0 0 0|0 0 3|0 0|0 0]
            [-----+-----+---+---]
            [0 0 0|0 0 0|3 1|0 0]
            [0 0 0|0 0 0|0 3|0 0]
            [-----+-----+---+---]
            [0 0 0|0 0 0|0 0|3 1]
            [0 0 0|0 0 0|0 0|0 3]
            sage: T * J * T**(-1) == A
            True
            sage: T.rank()
            10

        ::

            sage: A = matrix(QQ, [[15, 37/3, -16, -104/3, -29, -7/3, 35, 2/3, -29/3, -1/3], [2, 9, -1, -6, -6, 0, 7, 0, -2, 0], [24, 74/3, -29, -208/3, -58, -14/3, 70, 4/3, -58/3, -2/3], [-6, -19, 3, 21, 19, 0, -21, 0, 6, 0], [2, 6, -1, -6, -3, 0, 7, 0, -2, 0], [-96, -296/3, 128, 832/3, 232, 65/3, -279, -16/3, 232/3, 8/3], [0, 0, 0, 0, 0, 0, 3, 0, 0, 0], [20, 26/3, -30, -199/3, -42, -14/3, 70, 13/3, -55/3, -2/3], [18, 57, -9, -54, -57, 0, 63, 0, -15, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 3]]); A
            [    15   37/3    -16 -104/3    -29   -7/3     35    2/3  -29/3   -1/3]
            [     2      9     -1     -6     -6      0      7      0     -2      0]
            [    24   74/3    -29 -208/3    -58  -14/3     70    4/3  -58/3   -2/3]
            [    -6    -19      3     21     19      0    -21      0      6      0]
            [     2      6     -1     -6     -3      0      7      0     -2      0]
            [   -96 -296/3    128  832/3    232   65/3   -279  -16/3  232/3    8/3]
            [     0      0      0      0      0      0      3      0      0      0]
            [    20   26/3    -30 -199/3    -42  -14/3     70   13/3  -55/3   -2/3]
            [    18     57     -9    -54    -57      0     63      0    -15      0]
            [     0      0      0      0      0      0      0      0      0      3]
            sage: J, T = A.jordan_form(transformation=True); J
            [3 1 0|0 0|0 0|0 0|0]
            [0 3 1|0 0|0 0|0 0|0]
            [0 0 3|0 0|0 0|0 0|0]
            [-----+---+---+---+-]
            [0 0 0|3 1|0 0|0 0|0]
            [0 0 0|0 3|0 0|0 0|0]
            [-----+---+---+---+-]
            [0 0 0|0 0|3 1|0 0|0]
            [0 0 0|0 0|0 3|0 0|0]
            [-----+---+---+---+-]
            [0 0 0|0 0|0 0|3 1|0]
            [0 0 0|0 0|0 0|0 3|0]
            [-----+---+---+---+-]
            [0 0 0|0 0|0 0|0 0|3]
            sage: T * J * T**(-1) == A
            True
            sage: T.rank()
            10

        Verify that we smoothly move to QQ from ZZ (:trac:`12693`), i.e.
        we work in the vector space over the field::

            sage: M = matrix(((2,2,2),(0,0,0),(-2,-2,-2)))
            sage: J, P = M.jordan_form(transformation=True)
            sage: J; P
            [0 1|0]
            [0 0|0]
            [---+-]
            [0 0|0]
            [ 2  1  0]
            [ 0  0  1]
            [-2  0 -1]
            sage: J - ~P * M * P
            [0 0 0]
            [0 0 0]
            [0 0 0]
            sage: parent(M)
            Full MatrixSpace of 3 by 3 dense matrices over Integer Ring
            sage: parent(J) == parent(P) == MatrixSpace(QQ, 3)
            True
            sage: M.jordan_form(transformation=True) == (M/1).jordan_form(transformation=True)
            True

        By providing eigenvalues ourselves, we can compute the Jordan form even
        lacking a polynomial factorization algorithm.  ::

            sage: Qx = PolynomialRing(QQ, 'x11, x12, x13, x21, x22, x23, x31, x32, x33')
            sage: x11, x12, x13, x21, x22, x23, x31, x32, x33 = Qx.gens()
            sage: M = matrix(Qx, [[0, 0, x31], [0, 0, x21], [0, 0, 0]])    # This is a nilpotent matrix.
            sage: M.jordan_form(eigenvalues=[(0, 3)])
            [0 1|0]
            [0 0|0]
            [---+-]
            [0 0|0]
            sage: M.jordan_form(eigenvalues=[(0, 2)])
            Traceback (most recent call last):
            ...
            ValueError: The provided list of eigenvalues is not correct.
            sage: M.jordan_form(transformation=True, eigenvalues=[(0, 3)])
            (
            [0 1|0]
            [0 0|0]  [x31   0   1]
            [---+-]  [x21   0   0]
            [0 0|0], [  0   1   0]
            )

        TESTS:

        The base ring for the matrix needs to have a fraction field
        and it needs to be implemented.  ::

            sage: A = matrix(Integers(6), 2, 2, range(4))
            sage: A.jordan_form()
            Traceback (most recent call last):
            ...
            ValueError: Matrix entries must be from a field, not Ring of integers modulo 6
        """
        from sage.matrix.constructor import block_diagonal_matrix, jordan_block, diagonal_matrix
        from sage.combinat.partition import Partition

        if self.ncols() != self.nrows():
            raise ValueError, "Jordan normal form not implemented for non-square matrices."

        # Set ``n`` to the number of rows and handle trivial cases, regardless
        # of the underlying ring.
        n = self.nrows()
        if n == 0:
            if not transformation:
                return self
            else:
                return self, self
        elif n == 1:
            if not transformation:
                return self
            else:
                return self, self.parent().identity_matrix()

        if (base_ring is None and not self.base_ring().is_exact()) or \
            (not base_ring is None and not base_ring.is_exact()):
            raise ValueError("Jordan normal form not implemented over inexact rings.")

        if base_ring is None:
            A = self
            base_ring = self.base_ring()

        # make sure we're working with a field..
        if not base_ring.is_field():
            try:
                base_field = base_ring.fraction_field()
            except (NotImplementedError, TypeError, AttributeError):
                raise ValueError("Matrix entries must be from a field, not {0}".
                                 format(base_ring))
            A = self.change_ring(base_field)

        # Compute the eigenvalues of the matrix, with multiplicities.  Here,
        # ``evals`` is a list of pairs, each first entry a root and each
        # second entry the corresponding multiplicity.
        if eigenvalues is not None:
            if check_input:    # Checking input for sanity.
                C1 = A.charpoly()
                Polyring = C1.parent()
                C2 = Polyring.one()
                x = Polyring.gens()[0]
                for z, i in eigenvalues:
                    C2 *= (x - z) ** i
                if C1 != C2:
                    raise ValueError("The provided list of eigenvalues is not correct.")
            evals = eigenvalues
        else:
            evals = A.charpoly().roots()
        if sum([mult for (_,mult) in evals]) < n:
            raise RuntimeError("Some eigenvalue does not exist in %s."  %(A.base_ring()))

        # Compute the block information.  Here, ``blocks`` is a list of pairs,
        # each first entry a root and each second entry the size of a block.
        # Note that in general there is more than one block per eigenvalue!
        blocks = []
        for eval, mult in evals:
            if mult == 1:
                blocks.append((eval,1))
            else:
                B = A - diagonal_matrix([eval]*n, sparse=sparse)
                C = B
                ranks = [n, C.rank()]
                i = 0
                while ranks[i] > ranks[i+1] and ranks[i+1] > n-mult:
                    C = B*C
                    ranks.append(C.rank())
                    i = i+1
                diagram = [ranks[i]-ranks[i+1] for i in xrange(len(ranks)-1)]
                blocks.extend([(eval, i) \
                    for i in Partition(diagram).conjugate()])

        # ``J`` is the matrix in Jordan canonical form.  Note that the blocks
        # are ordered firstly by the eigenvalues, in the same order as obeyed
        # by ``.roots()``, and secondly by size from greatest to smallest.
        J = block_diagonal_matrix([jordan_block(eval, size, sparse=sparse) \
            for (eval, size) in blocks], subdivide=subdivide)

        if transformation:
            from itertools import groupby

            # ``jordan_chains`` is a dictionary with keys the eigenvalues.
            # For every eigenvalue, we consider all Jordan blocks and find
            # a Jordan chain for each, adding the chain (a sequence of
            # vectors) to the entry for the eigenvalue (which is a list).
            jordan_chains = {}
            for eval,_ in evals:
                jordan_chains[eval] = []

                # Let B be the matrix `A - eval Id`.
                B = A - eval

                block_sizes = [size for e,size in blocks if e == eval]
                block_size_pairs = [(val,len(list(c))) \
                    for val,c in groupby(block_sizes)]

                # Y is a list of vectors, spanning everything that we have
                # covered by the Jordan chains we developed so far.
                Y = []

                for l,count in block_size_pairs:

                    # There are ``count`` Jordan blocks of size ``l``
                    # associated to this eigenvalue.

                    # We want to find elements in `\ker B^l - \ker B^{l-1}`.
                    Vlarge = (B**l).right_kernel().basis()
                    Vsmall = (B**(l-1)).right_kernel().basis()

                    for i in range(count):
                        # Let v be any vector in `\ker B^l` not in the kernel
                        # of `\ker B^{l-1}` which is also not in the span(Y),
                        # and start a chain from there.
                        v = _jordan_form_vector_in_difference(Vlarge, Vsmall+Y)
                        chain = [v]
                        for i in range(l-1):
                            chain.append(B*chain[-1])
                        chain.reverse()
                        Y.extend(chain)
                        jordan_chains[eval].append(chain)

            # Now ``jordan_chains`` has all the columns of the transformation
            # matrix; we just need to put them in the right order.
            jordan_basis = []
            for eval,size in blocks:
                # Find a block with the right size
                for index,chain in enumerate(jordan_chains[eval]):
                    if len(chain)==size:
                        jordan_basis += jordan_chains[eval].pop(index)
                        break

            transformation_matrix = (A.parent()(jordan_basis)).transpose()

        if transformation:
            return J, transformation_matrix
        else:
            return J

    def is_diagonalizable(self, base_field=None):
        r"""
        Determines if the matrix is similar to a diagonal matrix.

        INPUT:

        - ``base_field`` - a new field to use for entries
          of the matrix.

        OUTPUT:

        If ``self`` is the matrix `A`, then it is diagonalizable
        if there is an invertible matrix `S` and a diagonal matrix
        `D` such that

        .. math::

            S^{-1}AS = D

        This routine returns ``True`` if ``self`` is diagonalizable.
        The diagonal entries of the matrix `D` are the eigenvalues
        of `A`.  It may be necessary to "increase" the base field to
        contain all of the eigenvalues.  Over the rationals, the field
        of algebraic numbers, :mod:`sage.rings.qqbar` is a good choice.

        To obtain the matrices `S` and `D` use the :meth:`jordan_form`
        method with the ``transformation=True`` keyword.

        ALGORITHM:

        For each eigenvalue, this routine checks that the algebraic
        multiplicity (number of occurences as a root of the characteristic
        polynomial) is equal to the geometric multiplicity (dimension
        of the eigenspace), which is sufficient to ensure a basis of
        eigenvectors for the columns of `S`.

        EXAMPLES:

        A matrix that is diagonalizable over the rationals, as evidenced
        by its Jordan form.  ::

            sage: A = matrix(QQ, [[-7, 16, 12,  0,    6],
            ...                   [-9, 15,  0,  12, -27],
            ...                   [ 9, -8, 11, -12,  51],
            ...                   [ 3, -4,  0,  -1,   9],
            ...                   [-1,  0, -4,   4, -12]])
            sage: A.jordan_form(subdivide=False)
            [ 2  0  0  0  0]
            [ 0  3  0  0  0]
            [ 0  0  3  0  0]
            [ 0  0  0 -1  0]
            [ 0  0  0  0 -1]
            sage: A.is_diagonalizable()
            True

        A matrix that is not diagonalizable over the rationals, as evidenced
        by its Jordan form.  ::

            sage: A = matrix(QQ, [[-3, -14, 2, -1, 15],
            ...                   [4, 6, -2, 3, -8],
            ...                   [-2, -14, 0, 0, 10],
            ...                   [3, 13, -2, 0, -11],
            ...                   [-1, 6, 1, -3, 1]])
            sage: A.jordan_form(subdivide=False)
            [-1  1  0  0  0]
            [ 0 -1  0  0  0]
            [ 0  0  2  1  0]
            [ 0  0  0  2  1]
            [ 0  0  0  0  2]
            sage: A.is_diagonalizable()
            False

        If any eigenvalue of a matrix is outside the base ring, then
        this routine raises an error.  However, the ring can be
        "expanded" to contain the eigenvalues.  ::

            sage: A = matrix(QQ, [[1,  0,  1,  1, -1],
            ...                   [0,  1,  0,  4,  8],
            ...                   [2,  1,  3,  5,  1],
            ...                   [2, -1,  1,  0, -2],
            ...                   [0, -1, -1, -5, -8]])

            sage: [e in QQ for e in A.eigenvalues()]
            [False, False, False, False, False]
            sage: A.is_diagonalizable()
            Traceback (most recent call last):
            ...
            RuntimeError: an eigenvalue of the matrix is not contained in Rational Field

            sage: [e in QQbar for e in A.eigenvalues()]
            [True, True, True, True, True]
            sage: A.is_diagonalizable(base_field=QQbar)
            True

        Other exact fields may be employed, though it will not always
        be possible to expand their base fields to contain all
        the eigenvalues.  ::

            sage: F.<b> = FiniteField(5^2)
            sage: A = matrix(F, [[      4, 3*b + 2, 3*b + 1, 3*b + 4],
            ...                  [2*b + 1,     4*b,       0,       2],
            ...                  [    4*b,   b + 2, 2*b + 3,       3],
            ...                  [    2*b,     3*b, 4*b + 4, 3*b + 3]])
            sage: A.jordan_form()
            [      4       1|      0       0]
            [      0       4|      0       0]
            [---------------+---------------]
            [      0       0|2*b + 1       1]
            [      0       0|      0 2*b + 1]
            sage: A.is_diagonalizable()
            False

            sage: F.<c> = QuadraticField(-7)
            sage: A = matrix(F, [[   c + 3,   2*c - 2,   -2*c + 2,     c - 1],
            ...                  [2*c + 10, 13*c + 15, -13*c - 17, 11*c + 31],
            ...                  [2*c + 10, 14*c + 10, -14*c - 12, 12*c + 30],
            ...                  [       0,   2*c - 2,   -2*c + 2,   2*c + 2]])
            sage: A.jordan_form(subdivide=False)
            [    4     0     0     0]
            [    0    -2     0     0]
            [    0     0 c + 3     0]
            [    0     0     0 c + 3]
            sage: A.is_diagonalizable()
            True

        A trivial matrix is diagonalizable, trivially.  ::

            sage: A = matrix(QQ, 0, 0)
            sage: A.is_diagonalizable()
            True

        A matrix must be square to be diagonalizable. ::

            sage: A = matrix(QQ, 3, 4)
            sage: A.is_diagonalizable()
            False

        The matrix must have entries from a field,
        and it must be an exact field.  ::

            sage: A = matrix(ZZ, 4, range(16))
            sage: A.is_diagonalizable()
            Traceback (most recent call last):
            ...
            ValueError: matrix entries must be from a field, not Integer Ring

            sage: A = matrix(RDF, 4, range(16))
            sage: A.is_diagonalizable()
            Traceback (most recent call last):
            ...
            ValueError: base field must be exact, not Real Double Field

        AUTHOR:

        - Rob Beezer (2011-04-01)
        """
        if not self.is_square():
            return False
        if not base_field is None:
            self = self.change_ring(base_field)
        if not self.base_ring().is_exact():
            raise ValueError('base field must be exact, not {0}'.format(self.base_ring()))
        if not self.base_ring().is_field():
            raise ValueError('matrix entries must be from a field, not {0}'.format(self.base_ring()))

        evals = self.charpoly().roots()
        if sum([mult for (_,mult) in evals]) < self._nrows:
            raise RuntimeError('an eigenvalue of the matrix is not contained in {0}'.format(self.base_ring()))

        # Obtaining a generic minimal polynomial requires much more
        # computation with kernels and their dimensions than the following.
        # However, if a derived class has a fast minimal polynomial routine
        # then overriding this by checking for repeated factors might be faster.

        # check equality of algebraic multiplicity and geometric multiplicity
        for e, am in evals:
            gm = (self - e).right_kernel().dimension()
            if am != gm:
                return False
        return True

    def is_similar(self, other, transformation=False):
        r"""
        Returns ``True`` if ``self`` and ``other`` are similar,
        i.e. related by a change-of-basis matrix.

        INPUT:

        - ``other`` - a matrix, which should be square, and of the same size
          as ``self``, where the entries of the matrix have a fraction field
          equal to that of ``self``.  Inexact rings are not supported.

        - ``transformation`` - default: ``False`` - if ``True``, the output
          will include the change-of-basis matrix.  See below for an exact
          description.

        OUTPUT:

        Two matrices, $A$ and $B$ are similar if there is an invertible
        matrix $S$ such that $A=S^{-1}BS$.  $S$ can be interpreted as a
        change-of-basis matrix if $A$ and $B$ are viewed as matrix
        representations of the same linear transformation.

        When ``transformation=False`` this method will return ``True`` if
        such a matrix $S$ exists, otherwise it will return ``False``.  When
        ``transformation=True`` the method returns a pair.  The first part
        of the pair is ``True`` or ``False`` depending on if the matrices
        are similar and the second part is the change-of-basis matrix, or
        ``None`` should it not exist.

        When the transformation matrix is requested, it will satisfy
        ``self = S.inverse()*other*S``.

        If the base rings for any of the matrices is the integers, the
        rationals, or the field of algebraic numbers (``QQbar``), then the
        matrices are converted to have ``QQbar`` as their base ring prior
        to checking the equality of the base rings.

        It is possible for this routine to fail over most fields, even when
        the matrices are similar.  However, since the field of algebraic
        numbers is algebraically closed, the routine will always produce
        a result for matrices with rational entries.

        EXAMPLES:

        The two matrices in this example were constructed to be similar.
        The computations happen in the field of algebraic numbers, but we
        are able to convert the change-of-basis matrix back to the rationals
        (which may not always be possible). ::

            sage: A = matrix(ZZ, [[-5, 2, -11],
            ...                   [-6, 7, -42],
            ...                   [0, 1, -6]])
            sage: B = matrix(ZZ, [[ 1, 12,  3],
            ...                   [-1, -6, -1],
            ...                   [ 0,  6,  1]])
            sage: A.is_similar(B)
            True
            sage: _, T = A.is_similar(B, transformation=True)
            sage: T
            [ 1.0000000000000? + 0.?e-13*I           0.?e-13 + 0.?e-13*I           0.?e-13 + 0.?e-13*I]
            [-0.6666666666667? + 0.?e-13*I 0.16666666666667? + 0.?e-14*I -0.8333333333334? + 0.?e-13*I]
            [ 0.6666666666667? + 0.?e-13*I           0.?e-13 + 0.?e-13*I  -0.333333333334? + 0.?e-13*I]
            sage: T.change_ring(QQ)
            [   1    0    0]
            [-2/3  1/6 -5/6]
            [ 2/3    0 -1/3]
            sage: A == T.inverse()*B*T
            True

        Other exact fields are supported.  ::

            sage: F.<a> = FiniteField(7^2)
            sage: A = matrix(F,[[2*a + 5, 6*a + 6,   a + 3],
            ...                 [  a + 3, 2*a + 2, 4*a + 2],
            ...                 [2*a + 6, 5*a + 5,     3*a]])
            sage: B = matrix(F,[[5*a + 5, 6*a + 4,   a + 1],
            ...                 [  a + 5, 4*a + 3, 3*a + 3],
            ...                 [3*a + 5,   a + 4, 5*a + 6]])
            sage: A.is_similar(B)
            True
            sage: B.is_similar(A)
            True
            sage: _, T = A.is_similar(B, transformation=True)
            sage: T
            [      1       0       0]
            [6*a + 1 4*a + 3 4*a + 2]
            [6*a + 3 3*a + 5 3*a + 6]
            sage: A == T.inverse()*B*T
            True

        Two matrices with different sets of eigenvalues, so they
        cannot possibly be similar. ::

            sage: A = matrix(QQ, [[ 2,  3, -3, -6],
            ...                   [ 0,  1, -2, -8],
            ...                   [-3, -3,  4,  3],
            ...                   [-1, -2,  2,  6]])
            sage: B = matrix(QQ, [[ 1,  1,  2,  4],
            ...                   [-1,  2, -3, -7],
            ...                   [-2,  3, -4, -7],
            ...                   [ 0, -1,  0,  0]])
            sage: A.eigenvalues() == B.eigenvalues()
            False
            sage: A.is_similar(B, transformation=True)
            (False, None)

        Similarity is an equivalence relation, so this routine computes
        a representative of the equivalence class for each matrix, the
        Jordan form, as provided by :meth:`jordan_form`.  The matrices
        below have identical eigenvalues (as evidenced by equal
        characteristic polynomials), but slightly different Jordan forms,
        and hence are not similar.  ::

            sage: A = matrix(QQ, [[ 19, -7, -29],
            ...                   [-16, 11,  30],
            ...                   [ 15, -7, -25]])
            sage: B = matrix(QQ, [[-38, -63,  42],
            ...                   [ 14,  25, -14],
            ...                   [-14, -21,  18]])
            sage: A.charpoly() == B.charpoly()
            True
            sage: A.jordan_form()
            [-3| 0  0]
            [--+-----]
            [ 0| 4  1]
            [ 0| 0  4]
            sage: B.jordan_form()
            [-3| 0| 0]
            [--+--+--]
            [ 0| 4| 0]
            [--+--+--]
            [ 0| 0| 4]
            sage: A.is_similar(B)
            False

        Obtaining the Jordan form requires computing the eigenvalues of
        the matrix, which may not lie in the field used for entries of
        the matrix.  So the routine first checks the characteristic
        polynomials - if they are unequal, then the matrices cannot be
        similar. However, when the characteristic polynomials are equal,
        we must examine the Jordan form. In this case, the method may fail,
        EVEN when the matrices are similar.  This is not the case for
        matrices over the integers, rationals or algebraic numbers,
        since the computations are done in the algebraically closed
        field of algebraic numbers.

        Here is an example where the similarity is obvious, but the
        routine fails to compute a result.  ::

            sage: F.<a> = FiniteField(7^2)
            sage: C = matrix(F,[[  a + 2, 5*a + 4],
            ...                 [6*a + 6, 6*a + 4]])
            sage: S = matrix(ZZ, [[0, 1],
            ...                   [1, 0]])
            sage: D = S.inverse()*C*S
            sage: C.is_similar(D)
            Traceback (most recent call last):
            ...
            ValueError: unable to compute Jordan canonical form for a matrix
            sage: C.jordan_form()
            Traceback (most recent call last):
            ...
            RuntimeError: Some eigenvalue does not exist in Finite Field in a of size 7^2.

        Inexact rings and fields are also not supported.  ::

            sage: A = matrix(CDF, 2, 2, range(4))
            sage: B = copy(A)
            sage: A.is_similar(B)
            Traceback (most recent call last):
            ...
            ValueError: unable to compute Jordan canonical form for a matrix

        Rectangular matrices and mismatched sizes return quickly.  ::

            sage: A = matrix(3, 2, range(6))
            sage: B = copy(A)
            sage: A.is_similar(B)
            False
            sage: A = matrix(2, 2, range(4))
            sage: B = matrix(3, 3, range(9))
            sage: A.is_similar(B, transformation=True)
            (False, None)

        If the fraction fields of the entries are unequal, it is an error,
        except in the case when the rationals gets promoted to the
        algebraic numbers.  ::

            sage: A = matrix(ZZ, 2, 2, range(4))
            sage: B = matrix(GF(2), 2, 2, range(4))
            sage: A.is_similar(B, transformation=True)
            Traceback (most recent call last):
            ...
            TypeError: matrices need to have entries with identical fraction fields, not Algebraic Field and Finite Field of size 2
            sage: A = matrix(ZZ, 2, 2, range(4))
            sage: B = matrix(QQbar, 2, 2, range(4))
            sage: A.is_similar(B)
            True

        Inputs are checked.  ::

            sage: A = matrix(ZZ, 2, 2, range(4))
            sage: A.is_similar('garbage')
            Traceback (most recent call last):
            ...
            TypeError: similarity requires a matrix as an argument, not garbage
            sage: B = copy(A)
            sage: A.is_similar(B, transformation='junk')
            Traceback (most recent call last):
            ...
            ValueError: transformation keyword must be True or False, not junk
        """
        import sage.matrix.matrix
        import sage.rings.qqbar
        if not sage.matrix.matrix.is_Matrix(other):
            raise TypeError('similarity requires a matrix as an argument, not {0}'.format(other))
        if transformation not in [True, False]:
            raise ValueError('transformation keyword must be True or False, not {0}'.format(transformation))
        # easy false situations
        if not self.is_square() or not other.is_square():
            if transformation:
                return (False, None)
            else:
                return False
        if self.nrows() != other.nrows():
            if transformation:
                return (False, None)
            else:
                return False
        # convert to fraction fields for base rings
        A = self.matrix_over_field()
        B = other.matrix_over_field()
        # move rationals to algebraically closed algebraic numbers
        if A.base_ring() == QQ:
            A = A.change_ring(sage.rings.qqbar.QQbar)
        if B.base_ring() == QQ:
            B = B.change_ring(sage.rings.qqbar.QQbar)
        # require identical base fields
        if A.base_ring() != B.base_ring():
            raise TypeError('matrices need to have entries with identical fraction fields, not {0} and {1}'.format(A.base_ring(), B.base_ring()))
        # unequal characteristic polynomials implies not similar
        # and avoids any problems with eigenvalues not in the base field
        if A.charpoly() != B.charpoly():
            if transformation:
                return (False, None)
            else:
                return False
        # now more precisely compare Jordan form, and optionally get transformations
        try:
            if transformation:
                JA, SA = A.jordan_form(transformation=True)
            else:
                JA = A.jordan_form(transformation=False)
        except StandardError:
            raise ValueError('unable to compute Jordan canonical form for a matrix')
        try:
            if transformation:
                JB, SB = B.jordan_form(transformation=True)
            else:
                JB = B.jordan_form(transformation=False)
        except StandardError:
            raise ValueError('unable to compute Jordan canonical form for a matrix')
        similar = (JA == JB)
        transform = None
        if similar and transformation:
            transform = SB*SA.inverse()
        if transformation:
            return (similar, transform)
        else:
            return similar

    def symplectic_form(self):
        r"""
        Find a symplectic form for self if self is an anti-symmetric,
        alternating matrix defined over a field.

        Returns a pair (F, C) such that the rows of C form a symplectic
        basis for self and F = C \* self \* C.transpose().

        Raises a ValueError if not over a field, or self is not
        anti-symmetric, or self is not alternating.

        Anti-symmetric means that `M = -M^t`. Alternating means
        that the diagonal of `M` is identically zero.

        A symplectic basis is a basis of the form
        `e_1, \ldots, e_j, f_1, \ldots f_j, z_1, \dots, z_k`
        such that

        - `z_i M v^t` = 0 for all vectors `v`

        - `e_i M {e_j}^t = 0` for all `i, j`

        - `f_i M {f_j}^t = 0` for all `i, j`

        - `e_i M {f_i}^t = 1` for all `i`

        - `e_i M {f_j}^t = 0` for all `i` not equal
          `j`.

        See the example for a pictorial description of such a basis.

        EXAMPLES::

            sage: E = matrix(QQ, 8, 8, [0, -1/2, -2, 1/2, 2, 0, -2, 1, 1/2, 0, -1, -3, 0, 2, 5/2, -3, 2, 1, 0, 3/2, -1, 0, -1, -2, -1/2, 3, -3/2, 0, 1, 3/2, -1/2, -1/2, -2, 0, 1, -1, 0, 0, 1, -1, 0, -2, 0, -3/2, 0, 0, 1/2, -2, 2, -5/2, 1, 1/2, -1, -1/2, 0, -1, -1, 3, 2, 1/2, 1, 2, 1, 0]); E
            [   0 -1/2   -2  1/2    2    0   -2    1]
            [ 1/2    0   -1   -3    0    2  5/2   -3]
            [   2    1    0  3/2   -1    0   -1   -2]
            [-1/2    3 -3/2    0    1  3/2 -1/2 -1/2]
            [  -2    0    1   -1    0    0    1   -1]
            [   0   -2    0 -3/2    0    0  1/2   -2]
            [   2 -5/2    1  1/2   -1 -1/2    0   -1]
            [  -1    3    2  1/2    1    2    1    0]
            sage: F, C = E.symplectic_form(); F
            [ 0  0  0  0  1  0  0  0]
            [ 0  0  0  0  0  1  0  0]
            [ 0  0  0  0  0  0  1  0]
            [ 0  0  0  0  0  0  0  1]
            [-1  0  0  0  0  0  0  0]
            [ 0 -1  0  0  0  0  0  0]
            [ 0  0 -1  0  0  0  0  0]
            [ 0  0  0 -1  0  0  0  0]
            sage: F == C * E * C.transpose()
            True
            """
        import sage.matrix.symplectic_basis
        return sage.matrix.symplectic_basis.symplectic_basis_over_field(self)

    def _cyclic_subspace(self, v):
        r"""
        Helper function for computing with cyclic (Krylov) subspaces.

        For a square matrix `A` and a vector `v`, the cyclic subspace
        is spanned by the vectors

        .. math::

            \{v, Av, A^2v, A^3v, \dots \}

        INPUT:

        - ``self`` - a square matrix over a field.

        - ``v`` - a vector with a degree equal to the size of the matrix.

        There is no explicit error-checking, it is the responsibility of
        the calling routine to provide accurate input.

        OUTPUT:

        Four related items are output.  Principally this routine
        determines the dimension of a cyclic subspace, but also
        creates two bases for the subspace.  Let `k` be the smallest
        integer such that `A^kv` is a linear combination of the
        products with smaller powers of `A`, i.e. the dimension
        of the cyclic subspace.

        - A list of the vectors `v, Av, A^2v,\dots, A^{k-1}v`
          (the "iterates").  These vectors give one basis of
          the subspace.

        - A list of scalars giving a linear combination of
          `v, Av, A^2v,\dots, A^kv` that equals the zero vector.
          This is the unique such set of such scalars where the
          last one in the list is 1.  These can be used to form
          a monic polynomial in `A` that has `v` in its right kernel.
          the length of this list is `k+1`.

        - Form a matrix whose rows are the linearly independent iterates.
          Augment with a `k\times k` identity matrix.  Apply row operations,
          scaling and adding multiples of rows, but never swap rows.  Do
          this to create `k` pivot columns.  The third output is this
          augmented, nearly row-reduced, matrix.  The rows of the left
          portion will form a basis for the subspace, while the right
          portion will record linear combinations of the iterates that
          equal these basis vectors.

        - A list of length `k` with the location of the pivots
          in the augmented matrix.  Specifically, entry  ``i``  of this
          list is the column index of the pivot column containing its
          lone 1 in row ``i``.

        ALGORITHM:

        This could be called an "online echelon form" routine.  As each
        new power of the matrix is built, the iterate is added to the bottom
        of the augmented matrix and row operations are used to update
        the pivot columns.  Rows are never swapped, so this is not
        strictly reduced row-echelon form, but the running time will
        be similar.  The main difference is that it "discovers" the
        dimension of the subspace as quickly as possible.

        EXAMPLE::

            sage: A = matrix(QQ, [[5,4,2,1],[0,1,-1,-1],[-1,-1,3,0],[1,1,-1,2]])
            sage: v = vector(QQ, [0,1,0,0])
            sage: (QQ^4).span([v, A*v, A^2*v, A^3*v]).dimension()
            3

            sage: iterates, poly, augmented, pivots = A._cyclic_subspace(v)

            sage: iterates
            [(0, 1, 0, 0), (4, 1, -1, 1), (23, 1, -8, 8)]
            sage: poly
            [-16, 24, -9, 1]
            sage: lindep = iterates + [A^3*v]
            sage: sum(poly[i]*lindep[i] for i in range(4))
            (0, 0, 0, 0)
            sage: B = sum(poly[i]*A^i for i in range(4))
            sage: v in B.right_kernel()
            True

            sage: augmented
            [    0     1     0     0     1     0     0]
            [    1     0     0     0  -7/9   8/9  -1/9]
            [    0     0     1    -1 -19/9  23/9  -4/9]
            sage: pivots
            [1, 0, 2]
            sage: transform = augmented[:, 4:7]
            sage: transform*matrix(iterates) == augmented[:, 0:4]
            True
            sage: (QQ^4).span(iterates) == (QQ^4).span(augmented[:, 0:4].rows())
            True

        AUTHOR:

        - Rob Beezer (2011-05-20)
        """
        cdef Py_ssize_t n, i, j, k, pivcol
        cdef Matrix aug
        n = self.ncols()
        aug = self.new_matrix(nrows=n+1, ncols=n+(n+1))
        iterate = v.__copy__()
        iterates = []
        pivots = []
        for k in range(n+1):
            for j in range(n):
                aug[k, j] = iterate[j]
            # record keeping in augmented identity matrix
            aug[k, n+k] = 1
            # clear out pivot cols of row k, using pivots of previous rows
            for i in range(k):
                aug.add_multiple_of_row(k, i, -aug[k, pivots[i]])
            # identify new pivot
            # no new pivot is all zeros, ie linear dependence
            pivcol = -1
            for j in range(n):
                if aug[k, j] != 0:
                    pivcol = j
                    pivots.append(pivcol)
                    break
            # scale pivot, and clear its column
            if pivcol != -1:
                aug.rescale_row(k, 1/aug[k, pivcol])
                for i in range(k):
                    aug.add_multiple_of_row(i, k, -aug[i, pivcol])
                iterates.append(iterate)
                iterate = self*iterate
            else:
                break
        poly = []
        for j in range(n, n+k+1):
            poly.append(aug[k, j])
        return iterates, poly, aug.submatrix(0, 0, k, n+k), pivots

    def cyclic_subspace(self, v, var=None, basis='echelon'):
        r"""
        Create a cyclic subspace for a vector, and optionally,
        a minimal polynomial for the iterated powers.

        These subspaces are also known as Krylov subspaces.  They are
        spanned by the vectors

        .. math::

            \{v, Av, A^2v, A^3v, \dots \}

        INPUT:

        - ``self`` - a square matrix with entries from a field.

        - ``v`` - a vector with a degree equal to the size of the matrix
          and entries compatible with the entries of the matrix.

        - ``var`` - default: ``None`` - if specified as a string or
          a generator of a polynomial ring, then this will be used
          to construct a polynomial reflecting a relation of linear
          dependence on the powers `A^iv` *and* this will cause
          the polynomial to be returned along with the subspace.
          A generator must create polynomials with coefficients from
          the same field as the matrix entries.

        - ``basis`` - default: ``echelon`` - the basis for the
          subspace is "echelonized" by default, but the keyword
          'iterates' will return a subspace with a user basis
          equal to the largest linearly independent
          set `\{v, Av, A^2v, A^3v, \dots, A^{k-1}v \}`.

        OUTPUT:

        Suppose `k` is the smallest power such that
        `\{v, Av, A^2v, A^3v, \dots, A^{k}v \}` is linearly
        dependent.  Then the subspace returned will have
        dimension `k` and be spanned by the powers `0` through
        `k-1`.

        If a polynomial is requested through the use of the
        ``var`` keyword, then a pair is returned, with the
        polynomial first and the subspace second.  The polynomial
        is the unique monic polynomial whose coefficients provide
        a relation of linear dependence on the first `k` powers.

        For less convenient, but more flexible output, see the
        helper method "_cyclic_subspace" in this module.

        EXAMPLES::

            sage: A = matrix(QQ, [[5,4,2,1],[0,1,-1,-1],[-1,-1,3,0],[1,1,-1,2]])
            sage: v = vector(QQ, [0,1,0,0])
            sage: E = A.cyclic_subspace(v); E
            Vector space of degree 4 and dimension 3 over Rational Field
            Basis matrix:
            [ 1  0  0  0]
            [ 0  1  0  0]
            [ 0  0  1 -1]
            sage: F = A.cyclic_subspace(v, basis='iterates'); F
            Vector space of degree 4 and dimension 3 over Rational Field
            User basis matrix:
            [ 0  1  0  0]
            [ 4  1 -1  1]
            [23  1 -8  8]
            sage: E == F
            True
            sage: p, S = A.cyclic_subspace(v, var='T'); p
            T^3 - 9*T^2 + 24*T - 16
            sage: gen = polygen(QQ, 'z')
            sage: p, S = A.cyclic_subspace(v, var=gen); p
            z^3 - 9*z^2 + 24*z - 16
            sage: p.degree() == E.dimension()
            True

        The polynomial has coefficients that yield a non-trivial
        relation of linear dependence on the iterates.  Or,
        equivalently, evaluating the polynomial with the matrix
        will create a matrix that annihilates the vector.  ::

            sage: A = matrix(QQ, [[15, 37/3, -16, -104/3, -29, -7/3, 35, 2/3, -29/3, -1/3],
            ...                   [ 2, 9, -1, -6, -6, 0, 7, 0, -2, 0],
            ...                   [24, 74/3, -29, -208/3, -58, -14/3, 70, 4/3, -58/3, -2/3],
            ...                   [-6, -19, 3, 21, 19, 0, -21, 0, 6, 0],
            ...                   [2, 6, -1, -6, -3, 0, 7, 0, -2, 0],
            ...                   [-96, -296/3, 128, 832/3, 232, 65/3, -279, -16/3, 232/3, 8/3],
            ...                   [0, 0, 0, 0, 0, 0, 3, 0, 0, 0],
            ...                   [20, 26/3, -30, -199/3, -42, -14/3, 70, 13/3, -55/3, -2/3],
            ...                   [18, 57, -9, -54, -57, 0, 63, 0, -15, 0],
            ...                   [0, 0, 0, 0, 0, 0, 0, 0, 0, 3]])
            sage: u = zero_vector(QQ, 10); u[0] = 1
            sage: p, S = A.cyclic_subspace(u, var='t', basis='iterates')
            sage: S
            Vector space of degree 10 and dimension 3 over Rational Field
            User basis matrix:
            [   1    0    0    0    0    0    0    0    0    0]
            [  15    2   24   -6    2  -96    0   20   18    0]
            [  79   12  140  -36   12 -560    0  116  108    0]
            sage: p
            t^3 - 9*t^2 + 27*t - 27
            sage: k = p.degree()
            sage: coeffs = p.list()
            sage: iterates = S.basis() + [A^k*u]
            sage: sum(coeffs[i]*iterates[i] for i in range(k+1))
            (0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
            sage: u in p(A).right_kernel()
            True

        TESTS:

        A small case.  ::

            sage: A = matrix(QQ, 5, range(25))
            sage: u = zero_vector(QQ, 5)
            sage: A.cyclic_subspace(u)
            Vector space of degree 5 and dimension 0 over Rational Field
            Basis matrix:
            []

        Various problem inputs.  Notice the vector must have entries
        that coerce into the base ring of the matrix, and a polynomial
        ring generator must have a base ring that agrees with the
        base ring of the matrix.  ::

            sage: A = matrix(QQ, 4, range(16))
            sage: v = vector(QQ, 4, range(4))

            sage: A.cyclic_subspace('junk')
            Traceback (most recent call last):
            ...
            TypeError: first input should be a vector, not junk

            sage: A.cyclic_subspace(v, var=sin(x))
            Traceback (most recent call last):
            ...
            TypeError: polynomial variable must be a string or polynomial ring generator, not sin(x)

            sage: t = polygen(GF(7), 't')
            sage: A.cyclic_subspace(v, var=t)
            Traceback (most recent call last):
            ...
            TypeError: polynomial generator must be over the same ring as the matrix entries

            sage: A.cyclic_subspace(v, basis='garbage')
            Traceback (most recent call last):
            ...
            ValueError: basis format must be 'echelon' or 'iterates', not garbage

            sage: B = matrix(QQ, 4, 5, range(20))
            sage: B.cyclic_subspace(v)
            Traceback (most recent call last):
            ...
            TypeError: matrix must be square, not 4 x 5

            sage: C = matrix(QQ, 5, 5, range(25))
            sage: C.cyclic_subspace(v)
            Traceback (most recent call last):
            ...
            TypeError: vector must have degree equal to the size of the matrix, not 4

            sage: D = matrix(RDF, 4, 4, range(16))
            sage: D.cyclic_subspace(v)
            Traceback (most recent call last):
            ...
            TypeError: matrix entries must be from an exact ring, not Real Double Field

            sage: E = matrix(Integers(6), 4, 4, range(16))
            sage: E.cyclic_subspace(v)
            Traceback (most recent call last):
            ...
            TypeError: matrix entries must be from an exact field, not Ring of integers modulo 6

            sage: F.<a> = GF(2^4)
            sage: G = matrix(QQ, 4, range(16))
            sage: w = vector(F, 4, [1, a, a^2, a^3])
            sage: G.cyclic_subspace(w)
            Traceback (most recent call last):
            ...
            TypeError: unable to make vector entries compatible with matrix entries

        AUTHOR:

        - Rob Beezer (2011-05-20)
        """
        import sage.rings.polynomial.polynomial_ring
        n = self.ncols()
        R = self.base_ring()
        if not is_Vector(v):
            raise TypeError('first input should be a vector, not {0}'.format(v))
        if not (var is None  or isinstance(var, basestring)):
            generator = False
            try:
                generator = var.is_gen()
            except AttributeError:
                pass
            if not generator:
                raise TypeError('polynomial variable must be a string or polynomial ring generator, not {0}'.format(var))
            elif var.base_ring() != R:
                raise TypeError('polynomial generator must be over the same ring as the matrix entries')
        if not basis in ['echelon', 'iterates']:
            raise ValueError("basis format must be 'echelon' or 'iterates', not {0}".format(basis))
        if not self.is_square():
            raise TypeError('matrix must be square, not {0} x {1}'.format(self.nrows(), self.ncols()))
        if v.degree() != n:
            raise TypeError('vector must have degree equal to the size of the matrix, not {0}'.format(v.degree()))
        if not (R.is_field() and R.is_exact()):
            try:
                fraction_field = R.fraction_field()
            except TypeError:
                raise TypeError('matrix entries must be from an exact field, not {0}'.format(R))
            if fraction_field.is_exact():
                return self.change_ring(R.fraction_field()).cyclic_subspace(v, var, basis)
            raise TypeError('matrix entries must be from an exact ring, not {0}'.format(R))
        try:
            v = v.change_ring(R)
        except TypeError:
            raise TypeError('unable to make vector entries compatible with matrix entries')

        iterates, poly, augmented, pivots = self._cyclic_subspace(v)
        k = len(pivots)
        polynomial = not var is None
        if polynomial:
            x = sage.rings.polynomial.polynomial_ring.polygen(R, var)
            poly = sum([poly[i]*x**i for i in range(len(poly))])
        ambient = R**n
        if basis == 'echelon':
            echelon = []
            pivot_col_row = zip(pivots, range(k))
            pivot_col_row.sort()
            aug = augmented.submatrix(0, 0, k, n)
            for _, pivrow in pivot_col_row:
                echelon.append(aug.row(pivrow))
            subspace = ambient.subspace(echelon, check=False, already_echelonized=True)
        elif basis == 'iterates':
            subspace = ambient.subspace_with_basis(iterates, check=False)
        if polynomial:
            return poly, subspace
        else:
            return subspace

    def cholesky_decomposition(self):
        r"""
        Return the Cholesky decomposition of ``self``.

        .. WARNING::

            ``cholesky_decomposition()`` is deprecated,
            please use :meth:`cholesky` instead.

        The computed decomposition is cached and returned on
        subsequent calls. Methods such as :meth:`solve_left` may also
        take advantage of the cached decomposition depending on the
        exact implementation.

        INPUT:

        The input matrix must be:

        - real, symmetric, and positive definite; or

        - complex, Hermitian, and positive definite.

        If not, a ``ValueError`` exception will be raised.

        OUTPUT:

        An immutable lower triangular matrix `L` such that `L L^t` equals ``self``.

        ALGORITHM:

        Calls the method ``_cholesky_decomposition_``, which by
        default uses a standard recursion.

        .. warning::

            This implementation uses a standard recursion that is not known to
            be numerically stable.

        .. warning::

            It is potentially expensive to ensure that the input is
            positive definite.  Therefore this is not checked and it
            is possible that the output matrix is *not* a valid
            Cholesky decomposition of ``self``.  An example of this is
            given in the tests below.

        EXAMPLES:

        Here is an example over the real double field; internally, this uses SciPy::

            sage: r = matrix(RDF, 5, 5, [ 0,0,0,0,1, 1,1,1,1,1, 16,8,4,2,1, 81,27,9,3,1, 256,64,16,4,1 ])
            sage: m = r * r.transpose(); m
            [    1.0     1.0     1.0     1.0     1.0]
            [    1.0     5.0    31.0   121.0   341.0]
            [    1.0    31.0   341.0  1555.0  4681.0]
            [    1.0   121.0  1555.0  7381.0 22621.0]
            [    1.0   341.0  4681.0 22621.0 69905.0]
            sage: L = m.cholesky_decomposition(); L
            doctest:...: DeprecationWarning:
            cholesky_decomposition() is deprecated; please use cholesky() instead.
            See http://trac.sagemath.org/13045 for details.
            [          1.0           0.0           0.0           0.0           0.0]
            [          1.0           2.0           0.0           0.0           0.0]
            [          1.0          15.0 10.7238052948           0.0           0.0]
            [          1.0          60.0 60.9858144589 7.79297342371           0.0]
            [          1.0         170.0 198.623524155 39.3665667796 1.72309958068]
            sage: L.parent()
            Full MatrixSpace of 5 by 5 dense matrices over Real Double Field
            sage: L*L.transpose()
            [ 1.0     1.0     1.0     1.0     1.0]
            [ 1.0     5.0    31.0   121.0   341.0]
            [ 1.0    31.0   341.0  1555.0  4681.0]
            [ 1.0   121.0  1555.0  7381.0 22621.0]
            [ 1.0   341.0  4681.0 22621.0 69905.0]
            sage: ( L*L.transpose() - m ).norm(1) < 2^-30
            True

        The result is immutable::

            sage: L[0,0] = 0
            Traceback (most recent call last):
                ...
            ValueError: matrix is immutable; please change a copy instead (i.e., use copy(M) to change a copy of M).

        Here is an example over a higher precision real field::

            sage: r = matrix(RealField(100), 5, 5, [ 0,0,0,0,1, 1,1,1,1,1, 16,8,4,2,1, 81,27,9,3,1, 256,64,16,4,1 ])
            sage: m = r * r.transpose()
            sage: L = m.cholesky_decomposition()
            sage: L.parent()
            Full MatrixSpace of 5 by 5 dense matrices over Real Field with 100 bits of precision
            sage: ( L*L.transpose() - m ).norm(1) < 2^-50
            True

        Here is a Hermitian example::

            sage: r = matrix(CDF, 2, 2, [ 1, -2*I, 2*I, 6 ]); r
            [   1.0 -2.0*I]
            [ 2.0*I    6.0]
            sage: r.eigenvalues()
            [0.298437881284, 6.70156211872]
            sage: ( r - r.conjugate().transpose() ).norm(1) < 1e-30
            True
            sage: L = r.cholesky_decomposition(); L
            [          1.0           0.0]
            [        2.0*I 1.41421356237]
            sage: ( r - L*L.conjugate().transpose() ).norm(1) < 1e-30
            True
            sage: L.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Complex Double Field

        TESTS:

        The following examples are not positive definite::

            sage: m = -identity_matrix(3).change_ring(RDF)
            sage: m.cholesky_decomposition()
            Traceback (most recent call last):
            ...
            ValueError: The input matrix was not symmetric and positive definite

            sage: m = -identity_matrix(2).change_ring(RealField(100))
            sage: m.cholesky_decomposition()
            Traceback (most recent call last):
            ...
            ValueError: The input matrix was not symmetric and positive definite

        Here is a large example over a higher precision complex field::

            sage: r = MatrixSpace(ComplexField(100), 6, 6).random_element()
            sage: m = r * r.conjugate().transpose()
            sage: m.change_ring(CDF) # for display purposes
            [                     4.03491289478    1.65865229397 + 1.20395241554*I -0.275377464753 - 0.392579363912*I   0.646019176609 - 1.80427378747*I    1.15898675468 + 2.35202344518*I   -1.07920143474 + 1.37815737417*I]
            [   1.65865229397 - 1.20395241554*I                      5.19463366777  0.192784646673 + 0.217539084881*I   -1.24630239913 - 1.00510523556*I    1.65179716714 + 1.27031304403*I   1.17275462994 + 0.565615358757*I]
            [-0.275377464753 + 0.392579363912*I  0.192784646673 - 0.217539084881*I                      2.04647180997 -0.550558880479 + 0.379933796418*I   1.00862850855 + 0.945317139306*I  -0.740344951784 - 0.46578741292*I]
            [  0.646019176609 + 1.80427378747*I   -1.24630239913 + 1.00510523556*I -0.550558880479 - 0.379933796418*I                      4.07227967662   -2.81600845862 - 0.56060176804*I  -2.28695708255 - 0.360066613053*I]
            [   1.15898675468 - 2.35202344518*I    1.65179716714 - 1.27031304403*I   1.00862850855 - 0.945317139306*I   -2.81600845862 + 0.56060176804*I                      5.26892394669  0.964830717551 + 0.111780339251*I]
            [  -1.07920143474 - 1.37815737417*I   1.17275462994 - 0.565615358757*I  -0.740344951784 + 0.46578741292*I  -2.28695708255 + 0.360066613053*I  0.964830717551 - 0.111780339251*I                      3.49912480167]
            sage: eigs = m.change_ring(CDF).eigenvalues() # again for display purposes
            sage: all(abs(imag(e)) < 1.3e-15 for e in eigs)
            True
            sage: [real(e) for e in eigs]
            [12.1121838768, 5.17714373118, 0.183583821657, 0.798520682956, 2.03399202232, 3.81092266261]

            sage: ( m - m.conjugate().transpose() ).norm(1) < 1e-50
            True
            sage: L = m.cholesky_decomposition(); L.change_ring(CDF)
            [                      2.00870926089                                 0.0                                 0.0                                 0.0                                 0.0                                 0.0]
            [  0.825730396261 - 0.599366189511*I                       2.03802923221                                 0.0                                 0.0                                 0.0                                 0.0]
            [ -0.137091748475 + 0.195438618996*I    0.20761467212 - 0.145606613292*I                       1.38750721467                                 0.0                                 0.0                                 0.0]
            [  0.321609099528 + 0.898225453828*I -0.477666770113 + 0.0346666053769*I  -0.416429223553 - 0.094835914364*I                       1.65839194165                                 0.0                                 0.0]
            [   0.576980839012 - 1.17091282993*I   0.232362216253 - 0.318581071175*I   0.880672963687 - 0.692440838276*I  -0.920603548686 + 0.566479149373*I                      0.992988116915                                 0.0]
            [ -0.537261143636 - 0.686091014267*I   0.591339766401 + 0.158450627525*I  -0.561877938537 + 0.106470627954*I  -0.871217053358 + 0.176042897482*I  0.0516893015902 + 0.656402869037*I                      0.902427551681]
            sage: ( m - L*L.conjugate().transpose() ).norm(1) < 1e-20
            True
            sage: L.parent()
            Full MatrixSpace of 6 by 6 dense matrices over Complex Field with 100 bits of precision
            sage: L[0,0] = 0
            Traceback (most recent call last):
                ...
            ValueError: matrix is immutable; please change a copy instead (i.e., use copy(M) to change a copy of M).

        Here is an example that returns an incorrect answer, because the input is *not* positive definite::

            sage: r = matrix(CDF, 2, 2, [ 1, -2*I, 2*I, 0 ]); r
            [   1.0 -2.0*I]
            [ 2.0*I    0.0]
            sage: r.eigenvalues()
            [2.56155281281, -1.56155281281]
            sage: ( r - r.conjugate().transpose() ).norm(1) < 1e-30
            True
            sage: L = r.cholesky_decomposition(); L
            [  1.0   0.0]
            [2.0*I 2.0*I]
            sage: L*L.conjugate().transpose()
            [   1.0 -2.0*I]
            [ 2.0*I    8.0]

        This test verifies that the caching of the two variants
        of the Cholesky decomposition have been cleanly separated.
        It can be safely removed as part of removing this method
        at the end of the deprecation period.
        (From :trac:`13045`.)  ::

            sage: r = matrix(CDF, 2, 2, [ 0, -2*I, 2*I, 0 ]); r
            [   0.0 -2.0*I]
            [ 2.0*I    0.0]
            sage: r.cholesky_decomposition()
            [        0.0         0.0]
            [NaN + NaN*I NaN + NaN*I]
            sage: r.cholesky()
            Traceback (most recent call last):
            ...
            ValueError: matrix is not positive definite
            sage: r[0,0] = 0  # clears cache
            sage: r.cholesky()
            Traceback (most recent call last):
            ...
            ValueError: matrix is not positive definite
            sage: r.cholesky_decomposition()
            [        0.0         0.0]
            [NaN + NaN*I NaN + NaN*I]
        """
        from sage.misc.superseded import deprecation
        deprecation(13045, "cholesky_decomposition() is deprecated; please use cholesky() instead.")
        assert self._nrows == self._ncols, "Can only Cholesky decompose square matrices"
        if self._nrows == 0:
            return self.__copy__()
        return self._cholesky_decomposition_()

    def _cholesky_decomposition_(self):
        r"""
        Return the Cholesky decomposition of ``self``; see ``cholesky_decomposition``.

        This generic implementation uses a standard recursion.
        """
        L = self.fetch('cholesky_broken')
        if L is None:
            A = self.__copy__()
            L = A.parent()(0)
            n = self.nrows()
            for k in range(0, n-1 + 1):
                try:
                    L[k, k] = A[k, k].sqrt()
                except TypeError:
                    raise ValueError, "The input matrix was not symmetric and positive definite"

                for s in range(k+1, n):
                    L[s, k] = A[s, k] / L[k, k]
                for j in range(k+1, n):
                    for i in range(j, n):
                        A[i, j] -= L[i, k]*L[j, k].conjugate()
            L.set_immutable()
            self.cache('cholesky_broken', L)
        return L

    def cholesky(self):
        r"""
        Returns the Cholesky decomposition of a symmetric or Hermitian matrix.

        INPUT:

        A square matrix that is real, symmetric and positive definite.
        Or a square matrix that is complex, Hermitian and positive
        definite.  Generally, the base ring for the entries of the
        matrix needs to be a subfield of the algebraic numbers
        (``QQbar``).  Examples include the rational numbers (``QQ``),
        some number fields, and real algebraic numbers and the
        algebraic numbers themselves.

        OUTPUT:

        For a matrix `A` the routine returns a lower triangular
        matrix `L` such that,

        .. math::

            A = LL^\ast

        where `L^\ast` is the conjugate-transpose in the complex case,
        and just the transpose in the real case.  If the matrix
        fails to be positive definite (perhaps because it is not
        symmetric or Hermitian), then a ``ValueError`` results.

        ALGORITHM:

        Whether or not the matrix is positive definite is checked
        first in every case.  This is accomplished with an
        indefinite factorization (see :meth:`indefinite_factorization`)
        which caches its result.  This algorithm is of an order `n^3/3`.
        If the matrix is positive definite, this computation always
        succeeds, using just field operations.  The transistion to a
        Cholesky decomposition "only" requires computing square roots
        of the positive (real) entries of the diagonal matrix produced in
        the indefinite factorization.  Hence, there is no real penalty
        in the positive definite check (here, or prior to calling this
        routine), but a field extension with square roots may not be
        implemented in all reasonable cases.

        EXAMPLES:

        This simple example has a result with entries that remain
        in the field of rational numbers.  ::

            sage: A = matrix(QQ, [[ 4, -2,  4,  2],
            ...                   [-2, 10, -2, -7],
            ...                   [ 4, -2,  8,  4],
            ...                   [ 2, -7,  4,  7]])
            sage: A.is_symmetric()
            True
            sage: L = A.cholesky()
            sage: L
            [ 2  0  0  0]
            [-1  3  0  0]
            [ 2  0  2  0]
            [ 1 -2  1  1]
            sage: L.parent()
            Full MatrixSpace of 4 by 4 dense matrices over Rational Field
            sage: L*L.transpose() == A
            True

        This seemingly simple example requires first moving to
        the rational numbers for field operations, and then square
        roots necessitate that the result has entries in the field
        of algebraic numbers.  ::

            sage: A = matrix(ZZ, [[ 78, -30, -37,  -2],
            ...                   [-30, 102, 179, -18],
            ...                   [-37, 179, 326, -38],
            ...                   [ -2, -18, -38,  15]])
            sage: A.is_symmetric()
            True
            sage: L = A.cholesky()
            sage: L
            [   8.83176086632785?                    0                    0                    0]
            [ -3.396831102433787?    9.51112708681461?                    0                    0]
            [ -4.189425026335004?   17.32383862241232?   2.886751345948129?                    0]
            [-0.2264554068289192?  -1.973397116652010?  -1.649572197684645?   2.886751345948129?]
            sage: L.parent()
            Full MatrixSpace of 4 by 4 dense matrices over Algebraic Field
            sage: L*L.transpose() == A
            True

        Some subfields of the complex numbers, such as this number
        field of complex numbers with rational real and imaginary parts,
        allow for this computation.  ::

            sage: C.<I> = QuadraticField(-1)
            sage: A = matrix(C, [[        23,  17*I + 3,  24*I + 25,     21*I],
            ...                  [ -17*I + 3,        38, -69*I + 89, 7*I + 15],
            ...                  [-24*I + 25, 69*I + 89,        976, 24*I + 6],
            ...                  [     -21*I, -7*I + 15,  -24*I + 6,       28]])
            sage: A.is_hermitian()
            True
            sage: L = A.cholesky()
            sage: L
            [                4.79...?                         0                       0        0]
            [   0.62...? - 3.54...?*I                  5.00...?                       0        0]
            [   5.21...? - 5.00...?*I   13.58...? + 10.72...?*I               24.98...?        0]
            [             -4.37...?*I   -0.10...? -  0.85...?*I  -0.21...? + 0.37...?*I 2.81...?]
            sage: L.parent()
            Full MatrixSpace of 4 by 4 dense matrices over Algebraic Field
            sage: (L*L.conjugate_transpose() - A.change_ring(QQbar)).norm() < 10^-10
            True

        The field of algebraic numbers is an ideal setting for this
        computation.  ::

            sage: A = matrix(QQbar, [[        2,   4 + 2*I,   6 - 4*I],
            ...                      [ -2*I + 4,        11, 10 - 12*I],
            ...                      [  4*I + 6, 10 + 12*I,        37]])
            sage: A.is_hermitian()
            True
            sage: L = A.cholesky()
            sage: L
            [                       1.414213562373095?         0                   0]
            [2.828427124746190? - 1.414213562373095?*I         1                   0]
            [4.242640687119285? + 2.828427124746190?*I  -2*I + 2  1.732050807568878?]
            sage: L.parent()
            Full MatrixSpace of 3 by 3 dense matrices over Algebraic Field
            sage: (L*L.conjugate_transpose() - A.change_ring(QQbar)).norm() < 10^-10
            True


        Results are cached, hence immutable.  Use the ``copy`` function
        if you need to make a change.  ::

            sage: A = matrix(QQ, [[ 4, -2,  4,  2],
            ...                   [-2, 10, -2, -7],
            ...                   [ 4, -2,  8,  4],
            ...                   [ 2, -7,  4,  7]])
            sage: L = A.cholesky()
            sage: L.is_immutable()
            True

            sage: from copy import copy
            sage: LC = copy(L)
            sage: LC[0,0] = 1000
            sage: LC
            [1000    0    0    0]
            [  -1    3    0    0]
            [   2    0    2    0]
            [   1   -2    1    1]

        There are a variety of situations which will prevent the computation of a Cholesky decomposition.

        The base ring must be exact.  For numerical work, create a
        matrix with a base ring of ``RDF`` or ``CDF`` and use the
        :meth:`~sage.matrix.matrix_double_dense.Matrix_double_dense.cholesky`
        method for matrices of that type. ::

            sage: F = RealField(100)
            sage: A = matrix(F, [[1.0, 3.0], [3.0, -6.0]])
            sage: A.cholesky()
            Traceback (most recent call last):
            ...
            TypeError: base ring of the matrix must be exact, not Real Field with 100 bits of precision

        The base ring may not have a fraction field.  ::

            sage: A = matrix(Integers(6), [[2, 0], [0, 4]])
            sage: A.cholesky()
            Traceback (most recent call last):
            ...
            ValueError: unable to check positive definiteness because
            Unable to create the fraction field of Ring of integers modulo 6

        The base field may not have elements that are comparable to zero.  ::

            sage: F.<a> = FiniteField(5^4)
            sage: A = matrix(F, [[2+a^3, 3], [3, 3]])
            sage: A.cholesky()
            Traceback (most recent call last):
            ...
            ValueError: unable to check positive definiteness because
            cannot convert computations from Finite Field in a of size 5^4 into real numbers

        The algebraic closure of the fraction field of the base ring may not be implemented.  ::

            sage: F = Integers(7)
            sage: A = matrix(F, [[4, 0], [0, 3]])
            sage: A.cholesky()
            Traceback (most recent call last):
            ...
            TypeError: base field needs an algebraic closure with square roots,
            not Ring of integers modulo 7

        The matrix may not be positive definite.  ::

            sage: C.<I> = QuadraticField(-1)
            sage: B = matrix(C, [[      2, 4 - 2*I, 2 + 2*I],
            ...                  [4 + 2*I,       8,    10*I],
            ...                  [2 - 2*I,   -10*I,      -3]])
            sage: B.is_positive_definite()
            False
            sage: B.cholesky()
            Traceback (most recent call last):
            ...
            ValueError: matrix is not positive definite,
            so cannot compute Cholesky decomposition

        The matrix could be positive semi-definite, and thus
        lack a Cholesky decomposition.  ::

            sage: A = matrix(QQ, [[21, 15, 12, -3],
            ...                   [15, 12,  9,  12],
            ...                   [12,  9,  7,  3],
            ...                   [-3,  12,  3,  8]])
            sage: A.is_positive_definite()
            False
            sage: [A[:i,:i].determinant() for i in range(1,A.nrows()+1)]
            [21, 27, 0, 0]
            sage: A.cholesky()
            Traceback (most recent call last):
            ...
            ValueError: matrix is not positive definite,
            so cannot compute Cholesky decomposition

        Even in light of the above, you can sometimes get lucky
        and arrive at a situation where a particular matrix has
        a Cholesky decomposition when the general characteristics
        of the matrix suggest this routine would fail. In this
        example, the indefinite factorization produces a
        diagonal matrix whose elements from the finite field
        convert naturally to positive integers and are also
        perfect squares.  ::

            sage: F.<a> = FiniteField(5^3)
            sage: A = matrix(F, [[         4,       2*a^2 + 3,         4*a + 1],
            ...                  [ 2*a^2 + 3,         2*a + 2, 4*a^2 + 4*a + 4],
            ...                  [   4*a + 1, 4*a^2 + 4*a + 4,       a^2 + 4*a]])
            sage: A.is_symmetric()
            True
            sage: L = A.cholesky()
            sage: L*L.transpose() == A
            True

        TESTS:

        This verifies that :trac:`11274` is resolved.  ::

            sage: E = matrix(QQ, [[2, 1], [1, 1]])
            sage: E.is_symmetric()
            True
            sage: E.eigenvalues()
            [0.38...?, 2.61...?]
            sage: E.det()
            1
            sage: E.cholesky()
            [ 1.414213562373095?                   0]
            [0.7071067811865475? 0.7071067811865475?]

        AUTHOR:

        - Rob Beezer (2012-05-27)
        """
        from copy import copy
        C = self.fetch('cholesky')
        if C is None:
            if not self.is_square():
                msg = "matrix must be square, not {0} x {1}"
                raise ValueError(msg.format(self.nrows(), self.ncols()))
            if not self.base_ring().is_exact():
                msg = 'base ring of the matrix must be exact, not {0}'
                raise TypeError(msg.format(self.base_ring()))
            try:
                posdef = self.is_positive_definite()
            except (ValueError, TypeError) as e:
                msg = "unable to check positive definiteness because {0}"
                raise ValueError(msg.format(e))
            if not posdef:
                msg = "matrix is not positive definite, so cannot compute Cholesky decomposition"
                raise ValueError(msg)
            # the successful positive definite check will cache a Hermitian
            # or symmetric indefinite factorization, as appropriate
            factors = self.fetch('indefinite_factorization_hermitian')
            if factors is None:
                factors = self.fetch('indefinite_factorization_symmetric')
            L = factors[0]
            d = factors[1]
            F = L.base_ring()  # field really
            splits = []        # square roots of diagonal entries
            for x in d:
                try:
                    sqrt = F(x.sqrt())
                except (TypeError, ValueError):
                    try:
                        F = F.algebraic_closure()
                    except (NotImplementedError, AttributeError):
                        msg = "base field needs an algebraic closure with square roots, not {0}"
                        raise TypeError(msg.format(F))
                    # try again
                    sqrt = F(x.sqrt())
                splits.append(sqrt)
            # move square root of the diagonal matrix
            # into the lower triangular matrix
            # We need a copy, to break immutability
            # and the field may have changed as well
            C = copy(L)
            if F != C.base_ring():
                C = C.change_ring(F)
            for c in range(C.ncols()):
                C.rescale_col(c, splits[c])
            C.set_immutable()
            self.cache('cholesky', C)
        return C

    def LU(self, pivot=None, format='plu'):
        r"""
        Finds a decomposition into a lower-triangular matrix and
        an upper-triangular matrix.

        INPUT:

        - ``pivot`` - pivoting strategy

          - 'auto' (default) - see if the matrix entries are
            ordered (i.e. if they have an absolute value method),
            and if so, use a the partial pivoting strategy.
            Otherwise, fall back to the nonzero strategy.  This
            is the best choice for general routines that may
            call this for matrix entries of a variety of types.

          - 'partial' - each column is examined for
            the element with the largest absolute value and the
            row containing this element is swapped into place.

          - 'nonzero' - the first nonzero element in a column
            is located and the row with this element is used.

        - ``format`` - contents of output, see more discussion
          below about output.

          - 'plu' (default) - a triple; matrices P, L and U
            such that A = P*L*U.

          - 'compact' - a pair; row permutation as a tuple, and the
            matrices L and U combined into one matrix.

        OUTPUT:

        Suppose that `A` is an `m\times n` matrix, then an LU
        decomposition is a lower-triangular `m\times m` matrix
        `L` with every diagonal element equal to 1, and an
        upper-triangular `m\times n` matrix, `U` such that the
        product `LU`, after a permutation of the rows, is then
        equal to `A`.  For the 'plu' format the permutation is
        returned as an `m\times m` permutation matrix `P` such
        that

        .. math::

            A = PLU

        It is more common to place the permutation matrix just
        to the left of `A`.  If you desire this version, then
        use the inverse of `P` which is computed most efficiently
        as its transpose.

        If the 'partial' pivoting strategy is used, then the
        non-diagonal entries of `L` will be less than or equal
        to 1 in absolute value.  The 'nonzero' pivot strategy may
        be faster, but the growth of data structures for elements of
        the decomposition might counteract the advantage.

        By necessity, returned matrices have a base ring equal
        to the fraction field of the base ring of the original matrix.

        In the 'compact' format, the first returned value is a
        tuple that is a permutation of the rows of `LU` that yields
        `A`.  See the doctest for how you might employ this
        permutation.  Then the matrices `L` and `U` are merged
        into one matrix -- remove the diagonal of ones in `L`
        and the remaining nonzero entries can replace the
        entries of `U` beneath the diagonal.

        The results are cached, only in the compact format, separately
        for each pivot strategy called.  Repeated requests for the
        'plu' format will require just a small amount of overhead
        in each call to bust out the compact format to the three
        matrices.  Since only the compact format is cached, the
        components of the compact format are immutable, while the
        components of the 'plu' format are regenerated, and hence
        are mutable.

        Notice that while `U` is similar to row-echelon form and the
        rows of `U` span the row space of `A`, the rows of `U` are not
        generally linearly independent.  Nor are the pivot columns
        (or rank) immediately obvious.  However for rings without
        specialized echelon form routines, this method is about
        twice as fast as the generic echelon form routine since
        it only acts "below the diagonal", as would be predicted
        from a theoretical analysis of the algorithms.

        .. note::

            This is an exact computation, so limited to exact
            rings. If you need numerical results, convert the
            base ring to the field of real double numbers,
            ``RDF`` or the field of complex double numbers,
            ``CDF``, which will use a faster routine that
            is careful about numerical subtleties.

        ALGORITHM:

            "Gaussian Elimination with Partial Pivoting,"
            Algorithm 21.1 of [TREFETHEN-BAU]_.

        EXAMPLES:

        Notice the difference in the `L` matrix as a result of different
        pivoting strategies.  With partial pivoting, every entry of `L`
        has absolute value 1 or less.  ::

            sage: A = matrix(QQ, [[1, -1,  0,  2,  4,  7, -1],
            ...                   [2, -1,  0,  6,  4,  8, -2],
            ...                   [2,  0,  1,  4,  2,  6,  0],
            ...                   [1,  0, -1,  8, -1, -1, -3],
            ...                   [1,  1,  2, -2, -1,  1,  3]])
            sage: P, L, U = A.LU(pivot='partial')
            sage: P
            [0 0 0 0 1]
            [1 0 0 0 0]
            [0 0 0 1 0]
            [0 0 1 0 0]
            [0 1 0 0 0]
            sage: L
            [   1    0    0    0    0]
            [ 1/2    1    0    0    0]
            [ 1/2  1/3    1    0    0]
            [   1  2/3  1/5    1    0]
            [ 1/2 -1/3 -2/5    0    1]
            sage: U
            [    2    -1     0     6     4     8    -2]
            [    0   3/2     2    -5    -3    -3     4]
            [    0     0  -5/3  20/3    -2    -4 -10/3]
            [    0     0     0     0   2/5   4/5     0]
            [    0     0     0     0   1/5   2/5     0]
            sage: A == P*L*U
            True
            sage: P, L, U = A.LU(pivot='nonzero')
            sage: P
            [1 0 0 0 0]
            [0 1 0 0 0]
            [0 0 1 0 0]
            [0 0 0 1 0]
            [0 0 0 0 1]
            sage: L
            [ 1  0  0  0  0]
            [ 2  1  0  0  0]
            [ 2  2  1  0  0]
            [ 1  1 -1  1  0]
            [ 1  2  2  0  1]
            sage: U
            [ 1 -1  0  2  4  7 -1]
            [ 0  1  0  2 -4 -6  0]
            [ 0  0  1 -4  2  4  2]
            [ 0  0  0  0  1  2  0]
            [ 0  0  0  0 -1 -2  0]
            sage: A == P*L*U
            True

        An example of the compact format.  ::

            sage: B = matrix(QQ, [[ 1,  3,  5,  5],
            ...                   [ 1,  4,  7,  8],
            ...                   [-1, -4, -6, -6],
            ...                   [ 0, -2, -5, -8],
            ...                   [-2, -6, -6, -2]])
            sage: perm, M = B.LU(format='compact')
            sage: perm
            (4, 3, 0, 1, 2)
            sage: M
            [  -2   -6   -6   -2]
            [   0   -2   -5   -8]
            [-1/2    0    2    4]
            [-1/2 -1/2  3/4    0]
            [ 1/2  1/2 -1/4    0]

        We can easily illustrate the relationships between
        the two formats with a square matrix.  ::

            sage: C = matrix(QQ, [[-2,  3, -2, -5],
            ...                   [ 1, -2,  1,  3],
            ...                   [-4,  7, -3, -8],
            ...                   [-3,  8, -1, -5]])
            sage: P, L, U = C.LU(format='plu')
            sage: perm, M = C.LU(format='compact')
            sage: (L - identity_matrix(4)) + U == M
            True
            sage: p = [perm[i]+1 for i in range(len(perm))]
            sage: PP = Permutation(p).to_matrix()
            sage: PP == P
            True

        For a nonsingular matrix, and the 'nonzero' pivot
        strategy there is no need to permute rows, so the
        permutation matrix will be the identity.  Furthermore,
        it can be shown that then the `L` and `U` matrices
        are uniquely determined by requiring `L` to have ones
        on the diagonal.  ::

            sage: D = matrix(QQ, [[ 1,  0,  2,  0, -2, -1],
            ...                   [ 3, -2,  3, -1,  0,  6],
            ...                   [-4,  2, -3,  1, -1, -8],
            ...                   [-2,  2, -3,  2,  1,  0],
            ...                   [ 0, -1, -1,  0,  2,  5],
            ...                   [-1,  2, -4, -1,  5, -3]])
            sage: P, L, U = D.LU(pivot='nonzero')
            sage: P
            [1 0 0 0 0 0]
            [0 1 0 0 0 0]
            [0 0 1 0 0 0]
            [0 0 0 1 0 0]
            [0 0 0 0 1 0]
            [0 0 0 0 0 1]
            sage: L
            [   1    0    0    0    0    0]
            [   3    1    0    0    0    0]
            [  -4   -1    1    0    0    0]
            [  -2   -1   -1    1    0    0]
            [   0  1/2  1/4  1/2    1    0]
            [  -1   -1 -5/2   -2   -6    1]
            sage: U
            [   1    0    2    0   -2   -1]
            [   0   -2   -3   -1    6    9]
            [   0    0    2    0   -3   -3]
            [   0    0    0    1    0    4]
            [   0    0    0    0 -1/4 -3/4]
            [   0    0    0    0    0    1]
            sage: D == L*U
            True

        The base ring of the matrix may be any field, or a ring
        which has a fraction field implemented in Sage. The ring
        needs to be exact (there is a numerical LU decomposition
        for matrices over ``RDF`` and ``CDF``).  Matrices returned
        are over the original field, or the fraction field of the
        ring.  If the field is not ordered (i.e. the absolute value
        function is not implemented), then the pivot strategy needs
        to be 'nonzero'.  ::

            sage: A = matrix(RealField(100), 3, 3, range(9))
            sage: P, L, U = A.LU()
            Traceback (most recent call last):
            ...
            TypeError: base ring of the matrix must be exact, not Real Field with 100 bits of precision

            sage: A = matrix(Integers(6), 3, 2, range(6))
            sage: A.LU()
            Traceback (most recent call last):
            ...
            TypeError: base ring of the matrix needs a field of fractions, not Ring of integers modulo 6

            sage: R.<y> = PolynomialRing(QQ, 'y')
            sage: B = matrix(R, [[y+1, y^2+y], [y^2, y^3]])
            sage: P, L, U = B.LU(pivot='partial')
            Traceback (most recent call last):
            ...
            TypeError: cannot take absolute value of matrix entries, try 'pivot=nonzero'
            sage: P, L, U = B.LU(pivot='nonzero')
            sage: P
            [1 0]
            [0 1]
            sage: L
            [          1           0]
            [y^2/(y + 1)           1]
            sage: U
            [  y + 1 y^2 + y]
            [      0       0]
            sage: L.base_ring()
            Fraction Field of Univariate Polynomial Ring in y over Rational Field
            sage: B == P*L*U
            True

            sage: F.<a> = FiniteField(5^2)
            sage: C = matrix(F, [[a + 3, 4*a + 4, 2, 4*a + 2],
            ...                  [3, 2*a + 4, 2*a + 4, 2*a + 1],
            ...                  [3*a + 1, a + 3, 2*a + 4, 4*a + 3],
            ...                  [a, 3, 3*a + 1, a]])
            sage: P, L, U = C.LU(pivot='nonzero')
            sage: P
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]
            sage: L
            [      1       0       0       0]
            [3*a + 3       1       0       0]
            [    2*a 4*a + 2       1       0]
            [2*a + 3       2 2*a + 4       1]
            sage: U
            [  a + 3 4*a + 4       2 4*a + 2]
            [      0   a + 1   a + 3 2*a + 4]
            [      0       0       1 4*a + 2]
            [      0       0       0       0]
            sage: L.base_ring()
            Finite Field in a of size 5^2
            sage: C == P*L*U
            True

        With no pivoting strategy given (i.e. ``pivot=None``)
        the routine will try to use partial pivoting, but then
        fall back to the nonzero strategy. For the nonsingular
        matrix below, we see evidence of pivoting when viewed
        over the rationals, and no pivoting over the integers
        mod 29.  ::

            sage: entries = [3, 20, 11, 7, 16, 28, 5, 15, 21, 23, 22, 18, 8, 23, 15, 2]
            sage: A = matrix(Integers(29), 4, 4, entries)
            sage: perm, _ = A.LU(format='compact'); perm
            (0, 1, 2, 3)
            sage: B = matrix(QQ, 4, 4, entries)
            sage: perm, _ = B.LU(format='compact'); perm
            (2, 0, 1, 3)

        The `U` matrix is only guaranteed to be upper-triangular.
        The rows are not necessarily linearly independent, nor are
        the pivots columns or rank in evidence.  ::

            sage: A = matrix(QQ, [[ 1, -4,  1,  0, -2,  1, 3,  3,  2],
            ...                   [-1,  4,  0, -4,  0, -4, 5, -7, -7],
            ...                   [ 0,  0,  1, -4, -1, -3, 6, -5, -6],
            ...                   [-2,  8, -1, -4,  2, -4, 1, -8, -7],
            ...                   [ 1, -4,  2, -4, -3,  2, 5,  6,  4]])
            sage: P, L, U = A.LU()
            sage: U
            [   -2     8    -1    -4     2    -4     1    -8    -7]
            [    0     0   1/2    -2    -1    -2   9/2    -3  -7/2]
            [    0     0   3/2    -6    -2     0  11/2     2   1/2]
            [    0     0     0     0  -1/3    -1   5/3  -5/3  -5/3]
            [    0     0     0     0   1/3    -3   7/3 -19/3 -19/3]
            sage: A.rref()
            [ 1 -4  0  4  0  0 -1 -1 -1]
            [ 0  0  1 -4  0  0  1  0 -1]
            [ 0  0  0  0  1  0 -2 -1 -1]
            [ 0  0  0  0  0  1 -1  2  2]
            [ 0  0  0  0  0  0  0  0  0]
            sage: A.pivots()
            (0, 2, 4, 5)

        TESTS:

        Unknown keywords are caught.  ::

            sage: A = matrix(ZZ, 2, range(4))
            sage: A.LU(pivot='junk')
            Traceback (most recent call last):
            ...
            ValueError: pivot strategy must be None, 'partial' or 'nonzero', not junk
            sage: A.LU(format='garbage')
            Traceback (most recent call last):
            ...
            ValueError: format must be 'plu' or 'compact', not garbage

        Components of the 'compact' format are immutable, while
        components of the 'plu' format are not.  ::

            sage: A = matrix(ZZ, 2, range(4))
            sage: perm, M = A.LU(format='compact')
            sage: perm[0] = 25
            Traceback (most recent call last):
            ...
            TypeError: 'tuple' object does not support item assignment
            sage: M.is_immutable()
            True
            sage: P, L, U = A.LU(format='plu')
            sage: all([A.is_mutable() for A in [P, L, U]])
            True

        Partial pivoting is based on the absolute values of entries
        of a column.  Trac #12208 shows that the return value of the
        absolute value must be handled carefully.  This tests that
        situation in the case of cylotomic fields.  ::

            sage: C = SymmetricGroup(5).character_table()
            sage: C.base_ring()
            Cyclotomic Field of order 1 and degree 1
            sage: P, L, U = C.LU(pivot='partial')
            sage: C == P*L*U
            True

        AUTHOR:

        - Rob Beezer (2011-04-26)
        """
        if not pivot in [None, 'partial', 'nonzero']:
            msg = "pivot strategy must be None, 'partial' or 'nonzero', not {0}"
            raise ValueError(msg.format(pivot))
        if not format in ['compact', 'plu']:
            msg = "format must be 'plu' or 'compact', not {0}"
            raise ValueError(msg.format(format))

        # exact rings only, must have fraction field
        R = self.base_ring()
        if not R.is_exact():
            msg = 'base ring of the matrix must be exact, not {0}'
            raise TypeError(msg.format(R))
        if not R.is_field():
            try:
                F = R.fraction_field()
            except StandardError:
                msg = 'base ring of the matrix needs a field of fractions, not {0}'
                raise TypeError(msg.format(R))
        else:
            F = R

        # 'nonzero' strategy passes through untouched
        # 'partial' survives iff field has absolute value
        # None will use 'partial' iff possible, else fallback to nonzero
        if pivot in [None, 'partial']:
            try:
                abs(F.an_element())
                pivot = 'partial'
            except StandardError:
                if pivot == 'partial':
                    msg = "cannot take absolute value of matrix entries, try 'pivot=nonzero'"
                    raise TypeError(msg)
                pivot = 'nonzero'
        partial = (pivot == 'partial')

        cdef Py_ssize_t m, n, d, i, j, k, p, max_location
        cdef Matrix M

        # can now access cache, else compute
        #   the compact version of LU decomposition
        key = 'LU_' + pivot
        compact = self.fetch(key)
        if compact is None:
            if F == R:
                M = self.__copy__()
            else:
                M = self.change_ring(F)
            m, n = M._nrows, M._ncols
            d = min(m, n)
            perm = range(m)
            zero = F(0)
            for k in range(d):
                max_location = -1
                if partial:
                    # abs() necessary to convert zero to the
                    # correct type for comparisons (Trac #12208)
                    max_entry = abs(zero)
                    for i in range(k,m):
                        entry = abs(M.get_unsafe(i,k))
                        if entry > max_entry:
                            max_location = i
                            max_entry = entry
                else:
                    for i in range(k,m):
                        if M.get_unsafe(i,k) != zero:
                            max_location = i
                            break
                if max_location != -1:
                    perm[k], perm[max_location] = perm[max_location], perm[k]
                    M.swap_rows(k, max_location)
                    for j in range(k+1, m):
                        scale = -M.get_unsafe(j,k)/M.get_unsafe(k,k)
                        M.set_unsafe(j,k, -scale)
                        for p in range(k+1,n):
                            M.set_unsafe(j,p, M.get_unsafe(j,p) + scale*M.get_unsafe(k,p))
            perm = tuple(perm)
            M.set_immutable()
            compact = (perm, M)
            self.cache(key, compact)

        if format == 'compact':
            return compact
        elif format == 'plu':
            import sage.matrix.constructor
            import sage.combinat.permutation
            perm = compact[0]
            M = compact[1].__copy__()
            F = M.base_ring()
            m, n = M._nrows, M._ncols
            d = min(m, n)
            zero = F(0)
            perm = [perm[i]+1 for i in range(m)]
            P = sage.combinat.permutation.Permutation(perm).to_matrix()
            L = sage.matrix.constructor.identity_matrix(F, m)
            for i in range(1, m):
                for k in range(min(i,d)):
                    L[i,k] = M[i,k]
                    M[i,k] = zero
            return P, L, M

    def _indefinite_factorization(self, algorithm, check=True):
        r"""
        Utility function to decomposes a symmetric or
        Hermitian matrix into a lower triangular matrix
        and tuple of elements for the diagonal of a diagonal matrix.

        INPUT:

        - ``self`` - a matrix that is symmetric or Hermitian,
          over a ring that has a fraction field implemented.

        - ``algorithm`` - ``'symmetric'`` or ``'hermitian'``,
          according to the corresponding property of the matrix.

        - ``check`` - default: ``True`` - if ``True`` then
          performs the check that the matrix is consistent with the
          ``algorithm`` keyword.

        OUTPUT:

        Given a square matrix ``A``, the routine returns a
        pair: a matrix ``L`` and a list ``d``.

        ``L`` is a unit lower-triangular matrix.  ``d`` is
        the entries of a diagonal matrix.  Suppose this diagonal
        matrix is ``D``.  Then, for a symmetric matrix, these items
        are related as:

        .. math::

            A = LDL^T

        For a Hermitian matrix, the transpose can be replaced by
        the conjugate-transpose.

        If any leading principal submatrix is singular, then the
        computation cannot be performed and a ``ValueError`` results.

        Results are cached, and hence are immutable.  Caching
        eliminates redundant computations across
        :meth:`indefinite_factorization_`, :meth:`is_positive_definite`
        and :meth:`cholesky_decomposition`.

        EXAMPLES:

        A simple symmetric matrix. ::

            sage: A = matrix(QQ, [[ 4, -2,  4,  2],
            ...                   [-2, 10, -2, -7],
            ...                   [ 4, -2,  8,  4],
            ...                   [ 2, -7,  4,  7]])
            sage: A.is_symmetric()
            True
            sage: L, d = A._indefinite_factorization('symmetric')
            sage: L
            [   1    0    0    0]
            [-1/2    1    0    0]
            [   1    0    1    0]
            [ 1/2 -2/3  1/2    1]
            sage: d
            (4, 9, 4, 1)
            sage: A == L*diagonal_matrix(QQ, d)*L.transpose()
            True

        A Hermitian matrix. ::

            sage: x = var('x')
            sage: C.<I> = NumberField(x^2 + 1)
            sage: A = matrix(C, [[        23,  17*I + 3,  24*I + 25,     21*I],
            ...                  [ -17*I + 3,        38, -69*I + 89, 7*I + 15],
            ...                  [-24*I + 25, 69*I + 89,        976, 24*I + 6],
            ...                  [     -21*I, -7*I + 15,  -24*I + 6,       28]])
            sage: A.is_hermitian()
            True
            sage: L, d = A._indefinite_factorization('hermitian')
            sage: L
            [                1                    0                         0  0]
            [  -17/23*I + 3/23                    1                         0  0]
            [ -24/23*I + 25/23  617/288*I + 391/144                         1  0]
            [         -21/23*I     -49/288*I - 1/48  1336/89885*I - 773/89885  1]
            sage: d
            (23, 576/23, 89885/144, 142130/17977)
            sage: A == L*diagonal_matrix(C, d)*L.conjugate_transpose()
            True

        A matrix may have a singular submatrix in the upper-left
        corner ("a leading principal submatrix") which will unavoidably
        lead to division by zero.  This is the only situation when this
        algorithm fails.  ::

            sage: A = matrix(QQ, [[4, 6, 1], [6, 9, 5], [1, 5, 2]])
            sage: B = A[0:2, 0:2]; B.determinant()
            0
            sage: A._indefinite_factorization('symmetric')
            Traceback (most recent call last):
            ...
            ValueError: 2x2 leading principal submatrix is singular,
            so cannot create indefinite factorization

        TESTS:

        The matrix must be square.  ::

            sage: A = matrix(QQ, 3, 2, range(6))
            sage: A._indefinite_factorization('symmetric')
            Traceback (most recent call last):
            ...
            ValueError: matrix must be square, not 3 x 2

        The lone argument must describe the matrix.  ::

            sage: A = matrix(QQ, [[1, 5], [5, 8]])
            sage: A._indefinite_factorization('junk')
            Traceback (most recent call last):
            ...
            ValueError: 'algorithm' must be 'symmetric' or 'hermitian',
            not junk

        The matrix must contain entries from an exact ring.  ::

            sage: F = RealField(100)
            sage: A = matrix(F, [[1.0, 3.0], [3.0, -6.0]])
            sage: A._indefinite_factorization('symmetric')
            Traceback (most recent call last):
            ...
            TypeError: entries of the matrix must be in an exact ring,
            not Real Field with 100 bits of precision

        The base ring must have a fraction field.  ::

            sage: A = matrix(Integers(6), [[1, 5], [5, 8]])
            sage: A._indefinite_factorization('symmetric')
            Traceback (most recent call last):
            ...
            TypeError: Unable to create the fraction field of
            Ring of integers modulo 6

        When ``check`` is ``True`` (the default), the matrix is
        checked to see if it conforms with the ``algorithm``
        keyword.  ::

            sage: A = matrix(QQ, 4, 4, range(16))
            sage: A._indefinite_factorization('symmetric', check=True)
            Traceback (most recent call last):
            ...
            ValueError: matrix is not symmetric (maybe try the 'hermitian' keyword)

            sage: A = matrix([[3, 2+3*I], [5+6*I, 12]])
            sage: A._indefinite_factorization('hermitian', check=True)
            Traceback (most recent call last):
            ...
            ValueError: matrix is not hermitian

        Results are cached and hence immutable, according
        to the ``algorithm``.  ::

            sage: A = matrix(QQ, [[ 4, -2,  4,  2],
            ...                   [-2, 10, -2, -7],
            ...                   [ 4, -2,  8,  4],
            ...                   [ 2, -7,  4,  7]])
            sage: Ls, ds = A._indefinite_factorization('symmetric')
            sage: Lh, dh = A._indefinite_factorization('hermitian')
            sage: Ls.is_immutable(), Lh.is_immutable()
            (True, True)
            sage: isinstance(ds, tuple), isinstance(dh, tuple)
            (True, True)

        AUTHOR:

        - Rob Beezer (2012-05-24)
        """
        # Implementation note:  L begins as a copy of self.
        # Entries below the diagonal are replaced as the loops proceed.
        # Entries above the diagonal are used to store entries of D*L^T,
        #   which halves the number of multiplications needed otherwise.
        # The list d_inv holds reciprocals of the diagonal entries
        # So, below the diagonal, the main computation is:
        # A_ij = A_ij - (1/d_j)*sum(d_k*L_ik*L_jk, 0 <= k < j)

        cdef Py_ssize_t m, i, j, k
        cdef Matrix L

        if not algorithm in ['symmetric', 'hermitian']:
            msg = "'algorithm' must be 'symmetric' or 'hermitian', not {0}"
            raise ValueError(msg.format(algorithm))
        cache_string = 'indefinite_factorization_' + algorithm
        factors = self.fetch(cache_string)
        if factors is None:
            R = self.base_ring()
            if not self.is_square():
                msg = "matrix must be square, not {0} x {1}"
                raise ValueError(msg.format(self.nrows(), self.ncols()))
            if not algorithm in ['symmetric', 'hermitian']:
                msg = "'algorithm' must be 'symmetric' or 'hermitian', not {0}"
                raise ValueError(msg.format(algorithm))
            if not R.is_exact():
                msg = "entries of the matrix must be in an exact ring, not {0}"
                raise TypeError(msg.format(R))
            try:
                F = R.fraction_field()
            except (NotImplementedError, TypeError):
                msg = 'Unable to create the fraction field of {0}'
                raise TypeError(msg.format(R))
            if check and algorithm == 'symmetric':
                if not self.is_symmetric():
                    msg = "matrix is not symmetric (maybe try the 'hermitian' keyword)"
                    raise ValueError(msg)
            if check and algorithm == 'hermitian':
                if not self.is_hermitian():
                    raise ValueError('matrix is not hermitian')
            conjugate = (algorithm == 'hermitian')
            # we need a copy no matter what, so we
            # (potentially) change to fraction field at the same time
            L = self.change_ring(F)
            m = L._nrows
            zero = F(0)
            one = F(1)
            d = []
            d_inv = []
            for i in range(m):
                for j in range(i+1):
                    t = L.get_unsafe(i, j)
                    if conjugate:
                        for k in range(j):
                            t -= L.get_unsafe(k,i)*L.get_unsafe(j,k).conjugate()
                    else:
                        for k in range(j):
                            t -= L.get_unsafe(k,i)*L.get_unsafe(j,k)
                    if i == j:
                        if t == zero:
                            msg = "{0}x{0} leading principal submatrix is singular, so cannot create indefinite factorization"
                            raise ValueError(msg.format(i+1))
                        d.append(t)
                        d_inv.append(one/t)
                        L.set_unsafe(i, i, one)
                    else:
                        L.set_unsafe(j, i, t)
                        L.set_unsafe(i, j, (d_inv[j] * t))
            # Triangularize output matrix
            for i in range(m):
                for j in range(i+1, m):
                    L.set_unsafe(i, j, zero)
            L.set_immutable()
            d = tuple(d)
            factors = (L, d)
            self.cache(cache_string, factors)
        return factors

    def indefinite_factorization(self, algorithm='symmetric', check=True):
        r"""
        Decomposes a symmetric or Hermitian matrix into a
        lower triangular matrix and a diagonal matrix.

        INPUT:

        - ``self`` - a square matrix over a ring.  The base ring
          must have an implemented fraction field.

        - ``algorithm`` - default: ``'symmetric'``.  Either
          ``'symmetric'`` or ``'hermitian'``, according to if
          the input matrix is symmetric or hermitian.

        - ``check`` - default: ``True`` - if ``True`` then
          performs the check that the matrix is consistent with the
          ``algorithm`` keyword.

        OUTPUT:

        A lower triangular matrix `L` with each diagonal element
        equal to `1` and a vector of entries that form a
        diagonal matrix `D`.  The vector of diagonal entries
        can be easily used to form the matrix, as demonstrated
        below in the examples.

        For a symmetric matrix, `A`, these will be related by

        .. math::

            A = LDL^T

        If `A` is Hermitian matrix, then the transpose of `L`
        should be replaced by the conjugate-transpose of `L`.

        If any leading principal submatrix (a square submatrix
        in the upper-left corner) is singular then this method will
        fail with a ``ValueError``.

        ALGORITHM:

        The algorithm employed only uses field operations,
        but the compuation of each diagonal entry has the potential
        for division by zero.  The number of operations is of order
        `n^3/3`, which is half the count for an LU decomposition.
        This makes it an appropriate candidate for solving systems
        with symmetric (or Hermitian) coefficient matrices.

        EXAMPLES:

        There is no requirement that a matrix be positive definite, as
        indicated by the negative entries in the resulting diagonal
        matrix.  The default is that the input matrix is symmetric. ::

            sage: A = matrix(QQ, [[ 3,  -6,   9,   6,  -9],
            ...                   [-6,  11, -16, -11,  17],
            ...                   [ 9, -16,  28,  16, -40],
            ...                   [ 6, -11,  16,   9, -19],
            ...                   [-9,  17, -40, -19,  68]])
            sage: A.is_symmetric()
            True
            sage: L, d = A.indefinite_factorization()
            sage: D = diagonal_matrix(d)
            sage: L
            [ 1  0  0  0  0]
            [-2  1  0  0  0]
            [ 3 -2  1  0  0]
            [ 2 -1  0  1  0]
            [-3  1 -3  1  1]
            sage: D
            [ 3  0  0  0  0]
            [ 0 -1  0  0  0]
            [ 0  0  5  0  0]
            [ 0  0  0 -2  0]
            [ 0  0  0  0 -1]
            sage: A == L*D*L.transpose()
            True

        Optionally, Hermitian matrices can be factored
        and the result has a similar property (but not
        identical).  Here, the field is all complex numbers
        with rational real and imaginary parts.  As theory
        predicts, the diagonal entries will be real numbers.  ::

            sage: C.<I> = QuadraticField(-1)
            sage: B = matrix(C, [[      2, 4 - 2*I, 2 + 2*I],
            ...                  [4 + 2*I,       8,    10*I],
            ...                  [2 - 2*I,   -10*I,      -3]])
            sage: B.is_hermitian()
            True
            sage: L, d = B.indefinite_factorization(algorithm='hermitian')
            sage: D = diagonal_matrix(d)
            sage: L
            [      1       0       0]
            [  I + 2       1       0]
            [ -I + 1 2*I + 1       1]
            sage: D
            [ 2  0  0]
            [ 0 -2  0]
            [ 0  0  3]
            sage: B == L*D*L.conjugate_transpose()
            True

        If a leading principal submatrix has zero determinant, this
        algorithm will fail.  This will never happen with a positive
        definite matrix.  ::

            sage: A = matrix(QQ, [[21, 15, 12, -2],
            ...                   [15, 12,  9,  6],
            ...                   [12,  9,  7,  3],
            ...                   [-2,  6,  3,  8]])
            sage: A.is_symmetric()
            True
            sage: A[0:3,0:3].det() == 0
            True
            sage: A.indefinite_factorization()
            Traceback (most recent call last):
            ...
            ValueError: 3x3 leading principal submatrix is singular,
            so cannot create indefinite factorization

        This algorithm only depends on field operations, so
        outside of the singular submatrix situation, any matrix
        may be factored.  This provides a reasonable alternative
        to the Cholesky decomposition.  ::

            sage: F.<a> = FiniteField(5^3)
            sage: A = matrix(F,
            ...       [[      a^2 + 2*a, 4*a^2 + 3*a + 4,       3*a^2 + a, 2*a^2 + 2*a + 1],
            ...        [4*a^2 + 3*a + 4,       4*a^2 + 2,             3*a, 2*a^2 + 4*a + 2],
            ...        [      3*a^2 + a,             3*a,       3*a^2 + 2, 3*a^2 + 2*a + 3],
            ...        [2*a^2 + 2*a + 1, 2*a^2 + 4*a + 2, 3*a^2 + 2*a + 3, 3*a^2 + 2*a + 4]])
            sage: A.is_symmetric()
            True
            sage: L, d = A.indefinite_factorization()
            sage: D = diagonal_matrix(d)
            sage: L
            [              1               0               0               0]
            [4*a^2 + 4*a + 3               1               0               0]
            [              3   4*a^2 + a + 2               1               0]
            [      4*a^2 + 4 2*a^2 + 3*a + 3 2*a^2 + 3*a + 1               1]
            sage: D
            [      a^2 + 2*a               0               0               0]
            [              0 2*a^2 + 2*a + 4               0               0]
            [              0               0 3*a^2 + 4*a + 3               0]
            [              0               0               0       a^2 + 3*a]
            sage: A == L*D*L.transpose()
            True

        AUTHOR:

        - Rob Beezer (2012-05-24)
        """
        from sage.modules.free_module_element import vector
        L, d = self._indefinite_factorization(algorithm, check=check)
        return L, vector(L.base_ring(), d)

    def is_positive_definite(self):
        r"""
        Determines if a real or symmetric matrix is positive definite.

        A square matrix `A` is postive definite if it is
        symmetric with real entries or Hermitan with complex entries,
        and for every non-zero vector `\vec{x}`

        .. math::

            \vec{x}^\ast A\vec{x} > 0

        Here `\vec{x}^\ast` is the conjugate-transpose, which can be
        simplified to just the transpose in the real case.

        ALGORITHM:

        A matrix is positive definite if and only if the
        diagonal entries from the indefinite factorization
        are all positive (see :meth:`indefinite_factorization`).
        So this algorithm is of order ``n^3/3`` and may be applied
        to matrices with elements of any ring that has a fraction
        field contained within the reals or complexes.

        INPUT:

        Any square matrix.

        OUTPUT:

        This routine will return ``True`` if the matrix is square,
        symmetric or Hermitian, and meets the condition above
        for the quadratic form.

        The base ring for the elements of the matrix needs to
        have a fraction field implemented and the computations
        that result from the indefinite factorization must be
        convertable to real numbers that are comparable to zero.

        EXAMPLES:

        A real symmetric matrix that is positive definite,
        as evidenced by the positive entries for the diagonal
        matrix of the indefinite factorization and the postive
        determinants of the leading principal submatrices. ::

            sage: A = matrix(QQ, [[ 4, -2,  4,  2],
            ...                   [-2, 10, -2, -7],
            ...                   [ 4, -2,  8,  4],
            ...                   [ 2, -7,  4,  7]])
            sage: A.is_positive_definite()
            True
            sage: _, d = A.indefinite_factorization(algorithm='symmetric')
            sage: d
            (4, 9, 4, 1)
            sage: [A[:i,:i].determinant() for i in range(1,A.nrows()+1)]
            [4, 36, 144, 144]

        A real symmetric matrix which is not positive definite, along
        with a vector that makes the quadratic form negative. ::

            sage: A = matrix(QQ, [[ 3,  -6,   9,   6,  -9],
            ...                   [-6,  11, -16, -11,  17],
            ...                   [ 9, -16,  28,  16, -40],
            ...                   [ 6, -11,  16,   9, -19],
            ...                   [-9,  17, -40, -19,  68]])
            sage: A.is_positive_definite()
            False
            sage: _, d = A.indefinite_factorization(algorithm='symmetric')
            sage: d
            (3, -1, 5, -2, -1)
            sage: [A[:i,:i].determinant() for i in range(1,A.nrows()+1)]
            [3, -3, -15, 30, -30]
            sage: u = vector(QQ, [2, 2, 0, 1, 0])
            sage: u.row()*A*u
            (-3)

        A real symmetric matrix with a singular leading
        principal submatrix, that is therefore not positive definite.
        The vector ``u`` makes the quadratic form zero.  ::

            sage: A = matrix(QQ, [[21, 15, 12, -2],
            ...                   [15, 12,  9,  6],
            ...                   [12,  9,  7,  3],
            ...                   [-2,  6,  3,  8]])
            sage: A.is_positive_definite()
            False
            sage: [A[:i,:i].determinant() for i in range(1,A.nrows()+1)]
            [21, 27, 0, -75]
            sage: u = vector(QQ, [1,1,-3,0])
            sage: u.row()*A*u
            (0)

        An Hermitian matrix that is positive definite.  ::

            sage: C.<I> = NumberField(x^2 + 1)
            sage: A = matrix(C, [[        23,  17*I + 3,  24*I + 25,     21*I],
            ...                  [ -17*I + 3,        38, -69*I + 89, 7*I + 15],
            ...                  [-24*I + 25, 69*I + 89,        976, 24*I + 6],
            ...                  [     -21*I, -7*I + 15,  -24*I + 6,       28]])
            sage: A.is_positive_definite()
            True
            sage: _, d = A.indefinite_factorization(algorithm='hermitian')
            sage: d
            (23, 576/23, 89885/144, 142130/17977)
            sage: [A[:i,:i].determinant() for i in range(1,A.nrows()+1)]
            [23, 576, 359540, 2842600]

        An Hermitian matrix that is not positive definite.
        The vector ``u`` makes the quadratic form negative.  ::

            sage: C.<I> = QuadraticField(-1)
            sage: B = matrix(C, [[      2, 4 - 2*I, 2 + 2*I],
            ...                  [4 + 2*I,       8,    10*I],
            ...                  [2 - 2*I,   -10*I,      -3]])
            sage: B.is_positive_definite()
            False
            sage: _, d = B.indefinite_factorization(algorithm='hermitian')
            sage: d
            (2, -2, 3)
            sage: [B[:i,:i].determinant() for i in range(1,B.nrows()+1)]
            [2, -4, -12]
            sage: u = vector(C, [-5 + 10*I, 4 - 3*I, 0])
            sage: u.row().conjugate()*B*u
            (-50)

        A positive definite matrix over an algebraically closed field.  ::

            sage: A = matrix(QQbar, [[        2,   4 + 2*I,   6 - 4*I],
            ...                      [ -2*I + 4,        11, 10 - 12*I],
            ...                      [  4*I + 6, 10 + 12*I,        37]])
            sage: A.is_positive_definite()
            True
            sage: [A[:i,:i].determinant() for i in range(1,A.nrows()+1)]
            [2, 2, 6]

        TESTS:

        If the base ring lacks a ``conjugate`` method, it
        will be assumed to not be Hermitian and thus symmetric.
        If the base ring does not make sense as a subfield of
        the reals, then this routine will fail since comparison
        to zero is meaningless.  ::

            sage: F.<a> = FiniteField(5^3)
            sage: a.conjugate()
            Traceback (most recent call last):
            ...
            AttributeError: 'sage.rings.finite_rings.element_givaro.FiniteField_givaroElement'
            object has no attribute 'conjugate'
            sage: A = matrix(F,
            ...       [[      a^2 + 2*a, 4*a^2 + 3*a + 4,       3*a^2 + a, 2*a^2 + 2*a + 1],
            ...        [4*a^2 + 3*a + 4,       4*a^2 + 2,             3*a, 2*a^2 + 4*a + 2],
            ...        [      3*a^2 + a,             3*a,       3*a^2 + 2, 3*a^2 + 2*a + 3],
            ...        [2*a^2 + 2*a + 1, 2*a^2 + 4*a + 2, 3*a^2 + 2*a + 3, 3*a^2 + 2*a + 4]])
            sage: A.is_positive_definite()
            Traceback (most recent call last):
            ...
            TypeError: cannot convert computations from
            Finite Field in a of size 5^3 into real numbers

        AUTHOR:

        - Rob Beezer (2012-05-24)
        """
        from sage.rings.all import RR
        # check to see if the Hermitian routine is usable
        # otherwise we will assume the symmetric case
        imaginary = True
        a = self.base_ring().an_element()
        try:
            a.conjugate()
        except AttributeError:
            imaginary = False
        if imaginary:
            if not self.is_hermitian():
                return False
        else:
            if not self.is_symmetric():
                return False
        try:
            if imaginary:
                _, d = self._indefinite_factorization('hermitian', check=False)
            else:
                _, d = self._indefinite_factorization('symmetric', check=False)
        except ValueError as e:
            # a zero singular leading principal submatrix is one
            # indicator that the matrix is not positive definite
            if str(e).find('leading principal submatrix is singular') != -1:
                return False
            else:
                raise ValueError(e)
        # Now have diagonal entries (hopefully real) and so can
        # test with a generator (which will short-circuit)
        # positive definite iff all entries of d are positive
        try:
            posdef = all( RR(x) > 0 for x in d )
        except TypeError:
            universe = Sequence(d).universe()
            msg = "cannot convert computations from {0} into real numbers"
            raise TypeError(msg.format(universe))
        return posdef

    def hadamard_bound(self):
        r"""
        Return an int n such that the absolute value of the determinant of
        this matrix is at most `10^n`.

        This is got using both the row norms and the column norms.

        This function only makes sense when the base field can be coerced
        to the real double field RDF or the MPFR Real Field with 53-bits
        precision.

        EXAMPLES::

            sage: a = matrix(ZZ, 3, [1,2,5,7,-3,4,2,1,123])
            sage: a.hadamard_bound()
            4
            sage: a.det()
            -2014
            sage: 10^4
            10000

        In this example the Hadamard bound has to be computed
        (automatically) using MPFR instead of doubles, since doubles
        overflow::

            sage: a = matrix(ZZ, 2, [2^10000,3^10000,2^50,3^19292])
            sage: a.hadamard_bound()
            12215
            sage: len(str(a.det()))
            12215
        """
        from sage.rings.all import RDF, RealField
        try:
            A = self.change_ring(RDF)
            m1 = A._hadamard_row_bound()
            A = A.transpose()
            m2 = A._hadamard_row_bound()
            return min(m1, m2)
        except (OverflowError, TypeError):
            # Try using MPFR, which handles large numbers much better, but is slower.
            import misc
            R = RealField(53, rnd='RNDU')
            A = self.change_ring(R)
            m1 = misc.hadamard_row_bound_mpfr(A)
            A = A.transpose()
            m2 = misc.hadamard_row_bound_mpfr(A)
            return min(m1, m2)

    def find(self,f, indices=False):
        r"""
        Find elements in this matrix satisfying the constraints in the
        function `f`. The function is evaluated on each element of
        the matrix .

        INPUT:


        -  ``f`` - a function that is evaluated on each
           element of this matrix.

        -  ``indices`` - whether or not to return the indices
           and elements of this matrix that satisfy the function.


        OUTPUT: If ``indices`` is not specified, return a
        matrix with 1 where `f` is satisfied and 0 where it is not.
        If ``indices`` is specified, return a dictionary
        containing the elements of this matrix satisfying `f`.

        EXAMPLES::

            sage: M = matrix(4,3,[1, -1/2, -1, 1, -1, -1/2, -1, 0, 0, 2, 0, 1])
            sage: M.find(lambda entry:entry==0)
            [0 0 0]
            [0 0 0]
            [0 1 1]
            [0 1 0]

        ::

            sage: M.find(lambda u:u<0)
            [0 1 1]
            [0 1 1]
            [1 0 0]
            [0 0 0]

        ::

            sage: M = matrix(4,3,[1, -1/2, -1, 1, -1, -1/2, -1, 0, 0, 2, 0, 1])
            sage: len(M.find(lambda u:u<1 and u>-1,indices=True))
            5

        ::

            sage: M.find(lambda u:u!=1/2)
            [1 1 1]
            [1 1 1]
            [1 1 1]
            [1 1 1]

        ::

            sage: M.find(lambda u:u>1.2)
            [0 0 0]
            [0 0 0]
            [0 0 0]
            [1 0 0]

        ::

            sage: sorted(M.find(lambda u:u!=0,indices=True).keys()) == M.nonzero_positions()
            True
        """

        cdef Py_ssize_t size,i,j
        cdef object M

        if indices is False:
            L = self._list()
            size = PyList_GET_SIZE(L)
            M = PyList_New(0)

            for i from 0 <= i < size:
                PyList_Append(M,<object>f(<object>PyList_GET_ITEM(L,i)))

            return matrix_space.MatrixSpace(IntegerModRing(2),
                                            nrows=self._nrows,ncols=self._ncols).matrix(M)

        else:
            # return matrix along with indices in a dictionary
            d = {}
            for i from 0 <= i < self._nrows:
                for j from 0 <= j < self._ncols:
                    if f(self.get_unsafe(i,j)):
                        d[(i,j)] = self.get_unsafe(i,j)

            return d

    def conjugate(self):
        r"""
        Return the conjugate of self, i.e. the matrix whose entries are the
        conjugates of the entries of self.

        EXAMPLES::

            sage: A = matrix(CDF, [[1+I,1],[0,2*I]])
            sage: A.conjugate()
            [1.0 - 1.0*I         1.0]
            [        0.0      -2.0*I]

        A matrix over a not-totally-real number field::

            sage: K.<j> = NumberField(x^2+5)
            sage: M = matrix(K, [[1+j,1], [0,2*j]])
            sage: M.conjugate()
            [-j + 1      1]
            [     0   -2*j]

        There is a shortcut for the conjugate::

            sage: M.C
            [-j + 1      1]
            [     0   -2*j]

        There is also a shortcut for the conjugate transpose, or "Hermitian transpose"::

            sage: M.H
            [-j + 1      0]
            [     1   -2*j]

        Conjugates work (trivially) for matrices over rings that embed
        canonically into the real numbers::

            sage: M = random_matrix(ZZ, 2)
            sage: M == M.conjugate()
            True
            sage: M = random_matrix(QQ, 3)
            sage: M == M.conjugate()
            True
            sage: M = random_matrix(RR, 2)
            sage: M == M.conjugate()
            True
        """
        return self.new_matrix(self.nrows(), self.ncols(), [z.conjugate() for z in self.list()])

    def conjugate_transpose(self):
        r"""
        Returns the transpose of ``self`` after each entry has been converted to its complex conjugate.

        .. note::
            This function is sometimes known as the "adjoint" of a matrix,
            though there is substantial variation and some confusion with
            the use of that term.

        OUTPUT:

        A matrix formed by taking the complex conjugate of every entry
        of ``self`` and then transposing the resulting matrix.

        Complex conjugation is implemented for many subfields
        of the complex numbers.  See the examples below, or more
        at :meth:`~conjugate`.

        EXAMPLES::

            sage: M = matrix(SR, 2, 2, [[2-I, 3+4*I], [9-6*I, 5*I]])
            sage: M.base_ring()
            Symbolic Ring
            sage: M.conjugate_transpose()
             [   I + 2  6*I + 9]
            [-4*I + 3     -5*I]

            sage: P = matrix(CC, 3, 2, [0.95-0.63*I, 0.84+0.13*I, 0.94+0.23*I, 0.23+0.59*I, 0.52-0.41*I, -0.50+0.90*I])
            sage: P.base_ring()
            Complex Field with 53 bits of precision
            sage: P.conjugate_transpose()
            [ 0.950... + 0.630...*I  0.940... - 0.230...*I  0.520... + 0.410...*I]
            [ 0.840... - 0.130...*I  0.230... - 0.590...*I -0.500... - 0.900...*I]

        There is also a shortcut for the conjugate transpose, or "Hermitian transpose"::

            sage: M.H
             [   I + 2  6*I + 9]
            [-4*I + 3     -5*I]

        Matrices over base rings that can be embedded in the
        real numbers will behave as expected. ::

            sage: P = random_matrix(QQ, 3, 4)
            sage: P.conjugate_transpose() == P.transpose()
            True

        The conjugate of a matrix is formed by taking conjugates of
        all the entries.  Some specialized subfields of the complex numbers
        are implemented in Sage and complex conjugation can be applied.
        (Matrices over quadratic number fields are another
        class of examples.) ::

            sage: C = CyclotomicField(5)
            sage: a = C.gen(); a
            zeta5
            sage: CC(a)
            0.309016994374947 + 0.951056516295154*I
            sage: M = matrix(C, 1, 2, [a^2, a+a^3])
            sage: M.conjugate_transpose()
            [             zeta5^3]
            [-zeta5^3 - zeta5 - 1]

        Conjugation does not make sense over rings not containing complex numbers. ::

            sage: N = matrix(GF(5), 2, [0,1,2,3])
            sage: N.conjugate_transpose()
            Traceback (most recent call last):
            ...
            AttributeError: 'sage.rings.finite_rings.integer_mod.IntegerMod_int' object has no attribute 'conjugate'

        AUTHOR:

            Rob Beezer (2010-12-13)
        """
        # limited testing on a 1000 x 1000 matrix over CC:
        #   transpose is fast, conjugate is slow
        #   so perhaps direct speed improvements to the conjugate() method
        return self.conjugate().transpose()

    def norm(self, p=2):
        r"""
        Return the p-norm of this matrix, where `p` can be 1, 2,
        `\inf`, or the Frobenius norm.

        INPUT:


        -  ``self`` - a matrix whose entries are coercible into
           CDF

        -  ``p`` - one of the following options:

        -  ``1`` - the largest column-sum norm

        -  ``2 (default)`` - the Euclidean norm

        -  ``Infinity`` - the largest row-sum norm

        -  ``'frob'`` - the Frobenius (sum of squares) norm


        OUTPUT: RDF number

        .. SEEALSO::

            - :func:`sage.misc.functional.norm`

        EXAMPLES::

            sage: A = matrix(ZZ, [[1,2,4,3],[-1,0,3,-10]])
            sage: A.norm(1)
            13.0
            sage: A.norm(Infinity)
            14.0
            sage: B = random_matrix(QQ, 20, 21)
            sage: B.norm(Infinity) == (B.transpose()).norm(1)
            True

        ::

            sage: Id = identity_matrix(12)
            sage: Id.norm(2)
            1.0
            sage: A = matrix(RR, 2, 2, [13,-4,-4,7])
            sage: A.norm()
            15.0

        Norms of numerical matrices over high-precision reals are computed by this routine.
        Faster routines for double precision entries from `RDF` or `CDF` are provided by
        the :class:`~sage.matrix.matrix_double_dense.Matrix_double_dense` class.  ::

            sage: A = matrix(CC, 2, 3, [3*I,4,1-I,1,2,0])
            sage: A.norm('frob')
            5.65685424949
            sage: A.norm(2)
            5.47068444321
            sage: A.norm(1)
            6.0
            sage: A.norm(Infinity)
            8.41421356237
            sage: a = matrix([[],[],[],[]])
            sage: a.norm()
            0.0
            sage: a.norm(Infinity) == a.norm(1)
            True
        """

        if self._nrows == 0 or self._ncols == 0:
            return RDF(0)

        # 2-norm:
        if p == 2:
            A = self.change_ring(CDF)
            A = A.conjugate().transpose() * A
            U, S, V = A.SVD()
            return max(S.list()).real().sqrt()

        A = self.apply_map(abs).change_ring(RDF)

        # 1-norm: largest column-sum
        if p == 1:
            A = A.transpose()
            return max([sum(i) for i in list(A)])

        # Infinity norm: largest row-sum
        if p == sage.rings.infinity.Infinity:
            return max([sum(i) for i in list(A)])

        # Frobenius norm: square root of sum of squares of entries of self
        if p == 'frob':
            return sum([i**2 for i in A.list()]).sqrt()

    def _numerical_approx(self,prec=None,digits=None):
        r"""
        Return a numerical approximation of ``self`` as either
        a real or complex number with at least the requested number of bits
        or digits of precision.

        INPUT:


        -  ``prec`` - an integer: the number of bits of
           precision

        -  ``digits`` - an integer: digits of precision


        OUTPUT: A matrix coerced to a real or complex field with prec bits
        of precision.

        EXAMPLES::

            sage: d = matrix([[3, 0],[0,sqrt(2)]]) ;
            sage: b = matrix([[1, -1], [2, 2]]) ; e = b * d * b.inverse();e
            [ 1/2*sqrt(2) + 3/2 -1/4*sqrt(2) + 3/4]
            [      -sqrt(2) + 3  1/2*sqrt(2) + 3/2]

        ::

            sage: e.numerical_approx(53)
            [ 2.20710678118655 0.396446609406726]
            [ 1.58578643762690  2.20710678118655]

        ::

            sage: e.numerical_approx(20)
            [ 2.2071 0.39645]
            [ 1.5858  2.2071]

        ::

            sage: (e-I).numerical_approx(20)
            [2.2071 - 1.0000*I           0.39645]
            [           1.5858 2.2071 - 1.0000*I]

        ::

            sage: M=matrix(QQ,4,[i/(i+1) for i in range(12)]);M
            [    0   1/2   2/3]
            [  3/4   4/5   5/6]
            [  6/7   7/8   8/9]
            [ 9/10 10/11 11/12]

        ::

            sage: M.numerical_approx()
            [0.000000000000000 0.500000000000000 0.666666666666667]
            [0.750000000000000 0.800000000000000 0.833333333333333]
            [0.857142857142857 0.875000000000000 0.888888888888889]
            [0.900000000000000 0.909090909090909 0.916666666666667]

        ::

            sage: matrix(SR, 2, 2, range(4)).n()
            [0.000000000000000  1.00000000000000]
            [ 2.00000000000000  3.00000000000000]

        ::

            sage: numerical_approx(M)
            [0.000000000000000 0.500000000000000 0.666666666666667]
            [0.750000000000000 0.800000000000000 0.833333333333333]
            [0.857142857142857 0.875000000000000 0.888888888888889]
            [0.900000000000000 0.909090909090909 0.916666666666667]

        """

        if prec is None:
            if digits is None:
                prec = 53
            else:
                prec = int(digits * 3.4) + 2

        try:
            return self.change_ring(sage.rings.real_mpfr.RealField(prec))
        except TypeError:
            # try to return a complex result
            return self.change_ring(sage.rings.complex_field.ComplexField(prec))

    #This line is added so the doctest are visible
    numerical_approx=_numerical_approx
    n=_numerical_approx
    N=_numerical_approx

    def plot(self, *args, **kwds):
        """
        A plot of this matrix.

        Each (ith, jth) matrix element is given a different color value
        depending on its relative size compared to the other elements in
        the matrix.

        The tick marks drawn on the frame axes denote the (ith, jth)
        element of the matrix.

        This method just calls ``matrix_plot``.
        ``*args`` and ``**kwds`` are passed to
        ``matrix_plot``.

        EXAMPLES:

        A matrix over ZZ colored with different grey levels::

            sage: A = matrix([[1,3,5,1],[2,4,5,6],[1,3,5,7]])
            sage: A.plot()

        Here we make a random matrix over RR and use cmap='hsv' to color
        the matrix elements different RGB colors (see documentation for
        ``matrix_plot`` for more information on cmaps)::

            sage: A = random_matrix(RDF, 50)
            sage: plot(A, cmap='hsv')

        Another random plot, but over GF(389)::

            sage: A = random_matrix(GF(389), 10)
            sage: A.plot(cmap='Oranges')
        """
        from sage.plot.plot import matrix_plot
        return matrix_plot(self, *args, **kwds)

    #added this to make doctests visible to users
    numerical_approx=_numerical_approx

    def derivative(self, *args):
        """
        Derivative with respect to variables supplied in args.

        Multiple variables and iteration counts may be supplied; see
        documentation for the global derivative() function for more
        details.

        EXAMPLES::

            sage: v = vector([1,x,x^2])
            sage: v.derivative(x)
            (0, 1, 2*x)
            sage: type(v.derivative(x)) == type(v)
            True
            sage: v = vector([1,x,x^2], sparse=True)
            sage: v.derivative(x)
            (0, 1, 2*x)
            sage: type(v.derivative(x)) == type(v)
            True
            sage: v.derivative(x,x)
            (0, 0, 2)
        """
        return multi_derivative(self, args)

    def exp(self):
        r"""
        Calculate the exponential of this matrix X, which is the matrix

        .. math::

           e^X = \sum_{k=0}^{\infty} \frac{X^k}{k!}.

        This function depends on maxima's matrix exponentiation
        function, which does not deal well with floating point
        numbers.  If the matrix has floating point numbers, they will
        be rounded automatically to rational numbers during the
        computation.  If you want approximations to the exponential
        that are calculated numerically, you may get better results by
        first converting your matrix to RDF or CDF, as shown in the
        last example.

        EXAMPLES::

            sage: a=matrix([[1,2],[3,4]])
            sage: a.exp()
            [-1/22*(e^sqrt(33)*(sqrt(33) - 11) - sqrt(33) - 11)*e^(-1/2*sqrt(33) + 5/2)              2/33*(sqrt(33)*e^sqrt(33) - sqrt(33))*e^(-1/2*sqrt(33) + 5/2)]
            [             1/11*(sqrt(33)*e^sqrt(33) - sqrt(33))*e^(-1/2*sqrt(33) + 5/2)  1/22*(e^sqrt(33)*(sqrt(33) + 11) - sqrt(33) + 11)*e^(-1/2*sqrt(33) + 5/2)]

            sage: type(a.exp())
            <type 'sage.matrix.matrix_symbolic_dense.Matrix_symbolic_dense'>

            sage: a=matrix([[1/2,2/3],[3/4,4/5]])
            sage: a.exp()
            [-1/418*(e^(1/10*sqrt(209))*(3*sqrt(209) - 209) - 3*sqrt(209) - 209)*e^(-1/20*sqrt(209) + 13/20)                   20/627*(sqrt(209)*e^(1/10*sqrt(209)) - sqrt(209))*e^(-1/20*sqrt(209) + 13/20)]
            [                  15/418*(sqrt(209)*e^(1/10*sqrt(209)) - sqrt(209))*e^(-1/20*sqrt(209) + 13/20)  1/418*(e^(1/10*sqrt(209))*(3*sqrt(209) + 209) - 3*sqrt(209) + 209)*e^(-1/20*sqrt(209) + 13/20)]

            sage: a=matrix(RR,[[1,pi.n()],[1e2,1e-2]])
            sage: a.exp()
            [ 1/11882424341266*(e^(3/1275529100*sqrt(227345670387496707609))*(11*sqrt(227345670387496707609) + 5941212170633) - 11*sqrt(227345670387496707609) + 5941212170633)*e^(-3/2551058200*sqrt(227345670387496707609) + 101/200)                            445243650/75781890129165569203*(sqrt(227345670387496707609)*e^(3/1275529100*sqrt(227345670387496707609)) - sqrt(227345670387496707609))*e^(-3/2551058200*sqrt(227345670387496707609) + 101/200)]
            [                                     10000/53470909535697*(sqrt(227345670387496707609)*e^(3/1275529100*sqrt(227345670387496707609)) - sqrt(227345670387496707609))*e^(-3/2551058200*sqrt(227345670387496707609) + 101/200) -1/11882424341266*(e^(3/1275529100*sqrt(227345670387496707609))*(11*sqrt(227345670387496707609) - 5941212170633) - 11*sqrt(227345670387496707609) - 5941212170633)*e^(-3/2551058200*sqrt(227345670387496707609) + 101/200)]
            sage: a.change_ring(RDF).exp()
            [42748127.3153 7368259.24416]
            [234538976.138 40426191.4516]
        """
        from sage.symbolic.ring import SR
        return self.change_ring(SR).exp()

    def elementary_divisors(self):
        r"""
        If self is a matrix over a principal ideal domain R, return
        elements `d_i` for `1 \le i \le k = \min(r,s)`
        where `r` and `s` are the number of rows and
        columns of self, such that the cokernel of self is isomorphic to

        .. math::

           R/(d_1) \oplus R/(d_2) \oplus R/(d_k)

        with `d_i \mid d_{i+1}` for all `i`. These are
        the diagonal entries of the Smith form of self (see
        :meth:`smith_form()`).

        EXAMPLES::

            sage: OE = EquationOrder(x^2 - x + 2, 'w')
            sage: w = OE.ring_generators()[0]
            sage: m = Matrix([ [1, w],[w,7]])
            sage: m.elementary_divisors()
            [1, -w + 9]

        .. seealso::

           :meth:`smith_form`
        """
        d, u, v = self.smith_form()
        r = min(self.nrows(), self.ncols())
        return [d[i,i] for i in xrange(r)]

    def smith_form(self):
        r"""
        If self is a matrix over a principal ideal domain R, return
        matrices D, U, V over R such that D = U \* self \* V, U and V have
        unit determinant, and D is diagonal with diagonal entries the
        ordered elementary divisors of self, ordered so that
        `D_{i} \mid D_{i+1}`. Note that U and V are not uniquely
        defined in general, and D is defined only up to units.

        INPUT:


        -  ``self`` - a matrix over an integral domain. If the
           base ring is not a PID, the routine might work, or else it will
           fail having found an example of a non-principal ideal. Note that we
           do not call any methods to check whether or not the base ring is a
           PID, since this might be quite expensive (e.g. for rings of
           integers of number fields of large degree).


        ALGORITHM: Lifted wholesale from
        http://en.wikipedia.org/wiki/Smith_normal_form

        .. seealso::

           :meth:`elementary_divisors`

        AUTHORS:

        - David Loeffler (2008-12-05)

        EXAMPLES:

        An example over the ring of integers of a number field (of class
        number 1)::

            sage: OE = NumberField(x^2 - x + 2,'w').ring_of_integers()
            sage: w = OE.ring_generators()[0]
            sage: m = Matrix([ [1, w],[w,7]])
            sage: d, u, v = m.smith_form()
            sage: (d, u, v)
            (
            [     1      0]  [ 1  0]  [ 1 -w]
            [     0 -w + 9], [-w  1], [ 0  1]
            )
            sage: u * m * v == d
            True
            sage: u.base_ring() == v.base_ring() == d.base_ring() == OE
            True
            sage: u.det().is_unit() and v.det().is_unit()
            True

        An example over the polynomial ring QQ[x]::

            sage: R.<x> = QQ[]; m=x*matrix(R,2,2,1) - matrix(R, 2,2,[3,-4,1,-1]); m.smith_form()
            (
            [            1             0]  [    0    -1]  [    1 x + 1]
            [            0 x^2 - 2*x + 1], [    1 x - 3], [    0     1]
            )

        An example over a field::

            sage: m = matrix( GF(17), 3, 3, [11,5,1,3,6,8,1,16,0]); d,u,v = m.smith_form()
            sage: d
            [1 0 0]
            [0 1 0]
            [0 0 0]
            sage: u*m*v == d
            True

        Some examples over non-PID's work anyway::

            sage: R = EquationOrder(x^2 + 5, 's') # class number 2
            sage: s = R.ring_generators()[0]
            sage: A = matrix(R, 2, 2, [s-1,-s,-s,2*s+1])
            sage: D, U, V = A.smith_form()
            sage: D, U, V
            (
            [     1      0]  [    4 s + 4]  [       1 -5*s + 6]
            [     0 -s - 6], [    s s - 1], [       0        1]
            )
            sage: D == U*A*V
            True

        Others don't, but they fail quite constructively::

            sage: matrix(R,2,2,[s-1,-s-2,-2*s,-s-2]).smith_form()
            Traceback (most recent call last):
            ...
            ArithmeticError: Ideal Fractional ideal (2, s + 1) not principal

        Empty matrices are handled safely::

            sage: m = MatrixSpace(OE, 2,0)(0); d,u,v=m.smith_form(); u*m*v == d
            True
            sage: m = MatrixSpace(OE, 0,2)(0); d,u,v=m.smith_form(); u*m*v == d
            True
            sage: m = MatrixSpace(OE, 0,0)(0); d,u,v=m.smith_form(); u*m*v == d
            True

        Some pathological cases that crashed earlier versions::

            sage: m = Matrix(OE, [[2*w,2*w-1,-w+1],[2*w+2,-2*w-1,w-1],[-2*w-1,-2*w-2,2*w-1]]); d, u, v = m.smith_form(); u * m * v == d
            True
            sage: m = matrix(OE, 3, 3, [-5*w-1,-2*w-2,4*w-10,8*w,-w,w-1,-1,1,-8]); d,u,v = m.smith_form(); u*m*v == d
            True
        """
        R = self.base_ring()
        left_mat = self.new_matrix(self.nrows(), self.nrows(), 1)
        right_mat = self.new_matrix(self.ncols(), self.ncols(), 1)
        if self == 0 or (self.nrows() <= 1 and self.ncols() <= 1):
            return self.__copy__(), left_mat, right_mat

        # data type checks on R
        if not R.is_integral_domain() or not R.is_noetherian():
            raise TypeError, "Smith form only defined over Noetherian integral domains"
        if not R.is_exact():
            raise NotImplementedError, "Smith form over non-exact rings not implemented at present"

        # first clear the first row and column
        u,t,v = _smith_onestep(self)

        # now recurse: t now has a nonzero entry at 0,0 and zero entries in the rest
        # of the 0th row and column, so we apply smith_form to the smaller submatrix
        mm = t.submatrix(1,1)
        dd, uu, vv = mm.smith_form()
        mone = self.new_matrix(1, 1, [1])
        d = dd.new_matrix(1,1,[t[0,0]]).block_sum(dd)
        u = uu.new_matrix(1,1,[1]).block_sum(uu) * u
        v = v * vv.new_matrix(1,1,[1]).block_sum(vv)
        dp, up, vp = _smith_diag(d)
        return dp,up*u,v*vp

    def hermite_form(self, include_zero_rows=True, transformation=False):
        """
        Return the Hermite form of self, if it is defined.

        INPUT:

            - ``include_zero_rows`` -- bool (default: True); if False
              the zero rows in the output matrix are deleted.

            - ``transformation`` -- bool (default: False) a matrix U such that U*self == H.

        OUTPUT:

            - matrix H
            - (optional) transformation matrix U such that U*self == H, possibly with zero
              rows deleted...


        EXAMPLES::

            sage: M = FunctionField(GF(7),'x').maximal_order()
            sage: K.<x> = FunctionField(GF(7)); M = K.maximal_order()
            sage: A = matrix(M, 2, 3, [x, 1, 2*x, x, 1+x, 2])
            sage: A.hermite_form()
            [      x       1     2*x]
            [      0       x 5*x + 2]
            sage: A.hermite_form(transformation=True)
            (
            [      x       1     2*x]  [1 0]
            [      0       x 5*x + 2], [6 1]
            )
            sage: A = matrix(M, 2, 3, [x, 1, 2*x, 2*x, 2, 4*x])
            sage: A.hermite_form(transformation=True, include_zero_rows=False)
            ([  x   1 2*x], [1 0])
            sage: H, U = A.hermite_form(transformation=True, include_zero_rows=True); H, U
            (
            [  x   1 2*x]  [1 0]
            [  0   0   0], [5 1]
            )
            sage: U*A == H
            True
            sage: H, U = A.hermite_form(transformation=True, include_zero_rows=False)
            sage: U*A
            [  x   1 2*x]
            sage: U*A == H
            True
        """
        left, H, pivots = self._echelon_form_PID()
        if not include_zero_rows:
            i = H.nrows() - 1
            while H.row(i) == 0:
                i -= 1
            H = H[:i+1]
            if transformation:
                left = left[:i+1]
        if transformation:
            return H, left
        else:
            return H

    def _echelon_form_PID(self):
        r"""
        Return a triple (left, a, pivots) where left*self == a and a is row
        echelon form (which in this case we define to mean Hermite normal
        form).

        When ideals of the base ring have a "small_residue" method (as is the
        case for number field ideals), we use this to reduce entries above each
        column pivot.

        AUTHOR:

        - David Loeffler (2009-06-01)

        - Moritz Minzlaff (2011-03-17): corrected code for matrices of one row;
          this fixed trac 9053

        EXAMPLES::

            sage: L.<a> = NumberField(x^3 - 2)
            sage: OL = L.ring_of_integers()

        We check some degenerate cases::

            sage: m = matrix(OL, 0, 0, []); r,s,p = m._echelon_form_PID();
            sage: (r,s,p)
            ([], [], [])
            sage: r * m == s and r.det() == 1
            True
            sage: m = matrix(OL, 0, 1, []); r,s,p = m._echelon_form_PID();
            sage: (r,s,p)
            ([], [], [])
            sage: r * m == s and r.det() == 1
            True
            sage: m = matrix(OL, 1, 0, []); r,s,p = m._echelon_form_PID();
            sage: (r,s,p)
            ([1], [], [])
            sage: r * m == s and r.det() == 1
            True

        A 2x2 matrix::

            sage: m = matrix(OL, 2, 2, [1,0, a, 2]);
            sage: r,s,p = m._echelon_form_PID(); (r,s,p)
            (
            [ 1  0]  [1 0]
            [-a  1], [0 2], [0, 1]
            )
            sage: r * m == s and r.det() == 1
            True

        A larger example::

            sage: m = matrix(OL, 3, 5, [a^2 - 3*a - 1, a^2 - 3*a + 1, a^2 + 1,
            ...     -a^2 + 2, -3*a^2 - a - 1, -6*a - 1, a^2 - 3*a - 1,
            ...     2*a^2 + a + 5, -2*a^2 + 5*a + 1, -a^2 + 13*a - 3,
            ...     -2*a^2 + 4*a - 2, -2*a^2 + 1, 2*a, a^2 - 6, 3*a^2 - a ])
            sage: r,s,p = m._echelon_form_PID()
            sage: s[2]
            (0, 0, -3*a^2 - 18*a + 34, -68*a^2 + 134*a - 53, -111*a^2 + 275*a - 90)
            sage: r * m == s and r.det() == 1
            True

        We verify that trac 9053 is resolved::

            sage: R.<x> = GF(7)[]
            sage: A = R^3
            sage: L = A.span([x*A.0 + (x^3 + 1)*A.1, x*A.2])
            sage: M = A.span([x*L.0])
            sage: M.0 in L
            True

        """
        if self.ncols() == 0:
            return self.new_matrix(self.nrows(), self.nrows(), 1), self, []

        if self.nrows() == 0:
            return self.new_matrix(self.nrows(), self.nrows(), 1), self, []

        if self.nrows() == 1:
            if self.is_zero():
                return self.new_matrix(self.nrows(), self.nrows(), 1), self, []
            else:
                return self.new_matrix(self.nrows(), self.nrows(), 1), self, [
                    self.nonzero_positions_in_row(0)[0] ]

        R = self.base_ring()

        # data type checks on R
        if not R.is_integral_domain():
            raise TypeError, ( "Generic echelon form only defined over "
                "integral domains" )
        if not R.is_exact():
            raise NotImplementedError, ( "Echelon form over generic non-exact "
                "rings not implemented at present" )

        left_mat, a = _generic_clear_column(self)
        assert left_mat * self == a

        if a[0,0] != 0:
            aa = a.submatrix(1, 1)
            s, t, pivs = aa._echelon_form_PID()
            left_mat = s.new_matrix(1,1,[1]).block_sum(s) * left_mat
            a = left_mat * self
            pivs = [0] + [x + 1 for x in pivs]

        else:
            aa = a.submatrix(0, 1)
            s, t, pivs = aa._echelon_form_PID()
            left_mat = s * left_mat
            a = left_mat * self
            pivs = [x+1 for x in pivs]


        try:
            for i in xrange(1, len(pivs)):
                y = a[i][pivs[i]]
                I = R.ideal(y)
                s = a[0][pivs[i]]
                t = I.small_residue(s)
                v = R( (s-t) / y)

                left_mat.add_multiple_of_row(0, i, -v)
                a.add_multiple_of_row(0, i, -v)
                assert left_mat * self == a
        except AttributeError: # on I.small_residue
            pass

        return left_mat, a, pivs

    def _zigzag_form(self, basis=True):
        r"""
        Helper method for computation of ZigZag form.

        INPUT:

        - ``self`` - a square matrix over an exact field.

        - ``basis`` - default; ``True`` - controls whether or not to
          compute a change-of-basis matrix (also called a transformation
          matrix).

        OUTPUT:

        See the documentation for the :meth:`zigzag_form` method for a
        description of ZigZag form and notes on the algorithm employed.

        When ``basis=True`` four items are returned.  The first is a matrix
        whose columns form a basis so that a matrix representation of ``self``
        relative to this basis will be the ZigZag form.  More precisely, if
        `A` is the matrix, `U` is the change-of-basis matrix and `Z` is
        the ZigZag form, then

        .. math::

            U^{-1}*A*U = Z

        The second item returned is the ZigZag form of the matrix.  The
        third item is a list of lists, where each list is the coefficients
        of the polynomial associated with a companion matrix in the form.
        Finally, the fourth item is a list of the entries in the upper-left
        corner of each off-diagonal block, which will be a zero or a one.
        The length of the list of corner entries will be one less than the
        length of list of polynomials.

        When ``basis=False`` only three items are returned.  These are
        just as above, but without the change-of-basis matrix.

        The compuation of the change-of-basis matrix has not been optimized.
        As a helper method, no error checking is performed on the inputs -
        that should be performed by the calling method.

        ALGORITHM:

        ZigZag form, and its computation, are due to Arne Storjohann
        and are  described in [STORJOHANN-THESIS]_ and
        [STORJOHANN-ISACC98]_, where the former is more
        representative of the code here.

        EXAMPLE:

            sage: A = matrix(QQ, [[-68,   69, -27, -11, -65,   9, -181, -32],
            ...                   [-52,   52, -27,  -8, -52, -16, -133, -14],
            ...                   [ 92,  -97,  47,  14,  90,  32,  241,  18],
            ...                   [139, -144,  60,  18, 148, -10,  362,  77],
            ...                   [ 40,  -41,  12,   6,  45, -24,  105,  42],
            ...                   [-46,   48, -20,  -7, -47,   0, -122, -22],
            ...                   [-26,   27, -13,  -4, -29,  -6,  -66, -14],
            ...                   [-33,   34, -13,  -5, -35,   7,  -87, -23]])
            sage: Z, polys, corners = A._zigzag_form(basis=False)
            sage: Z
            [  0   0   0  40   1   0   0   0]
            [  1   0   0  52   0   0   0   0]
            [  0   1   0  18   0   0   0   0]
            [  0   0   1  -1   0   0   0   0]
            [  0   0   0   0   0   1   0   0]
            [  0   0   0   0 -25  10   0   0]
            [  0   0   0   0   1   0   0  -4]
            [  0   0   0   0   0   0   1  -4]
            sage: polys
            [[-40, -52, -18, 1, 1], [25, -10, 1], [4, 4, 1]]
            sage: corners
            [1, 1]
            sage: len(polys) == len(corners) + 1
            True

            sage: U, X, polys, corners = A._zigzag_form()
            sage: 38416*U
            [        38416      -2612288      -8797264     -47943168   83753660284  -16775074316    1574993574  -12141218334]
            [            0      -1997632      -7222208     -38108672   66323477868  -13362573684    3694412232  -12414630592]
            [            0       3534272      12523616      66920672 -116699384082   23442252424   -4176610942   18012693216]
            [            0       5339824      18631760     100842000 -175585304193   35191011843   -3619613634   22708786576]
            [            0       1536640       5301408      28812000  -50432310506   10071675864      79884994    5258177064]
            [            0      -1767136      -6261808     -33460336   58377759731  -11719115889    1849749314   -8725304756]
            [            0       -998816      -3611104     -19054336   33189768208   -6681286842    1651135654   -5815174372]
            [            0      -1267728      -4456256     -23933168   41831027339   -8380482791     785625330   -5536675718]
            sage: X == Z
            True
            sage: U.inverse()*A*U == X
            True

        AUTHOR:

        - Rob Beezer (2011-06-09)
        """
        import sage.matrix.constructor
        cdef Py_ssize_t n, s, c, i, j, k

        n = self._ncols
        R = self.base_ring()
        zero = R(0)
        one = R(1)
        cdef Matrix Z
        Z = self.__copy__()
        polys = []    # coefficients for polynomials of companion matrices
        corners = []  # zero or one in corner of off-diagonal blocks
        if basis:
            U = sage.matrix.constructor.identity_matrix(R, n) # transformation matrix
        # parity switch, True iff working on transpose
        # if False, mimic row operations only on U
        # if True,  mimic column operations only on U
        trans = False
        s = 0  # index of top row of current block
        c = 0  # index of current row of current block
        while s < n:
            zigging = True
            # check for totally zero column below diagonal
            while zigging:  # zigging means we are building a block
                nonzero = -1
                for i in range(c+1, n):
                    if Z[i, c] != 0:
                        nonzero = i
                        break
                zigging = (nonzero != -1)
                if zigging:
                    Z.swap_rows(nonzero, c+1)
                    Z.swap_columns(nonzero, c+1)
                    if basis:
                        if trans:
                            U.swap_columns(nonzero, c+1)
                        else:
                            U.swap_rows(nonzero, c+1)

                    # create a subdiagonal entry 1
                    scale = Z[c+1, c]
                    Z.rescale_row(c+1, 1/scale)
                    Z.rescale_col(c+1, scale)
                    if basis:
                        if trans:
                            U.rescale_col(c+1, scale)
                        else:
                            U.rescale_row(c+1, 1/scale)

                    # clear column throughout the block,and in all rows below
                    for i in range(s, n):
                        if i != c+1:
                            scale = Z[i, c]
                            Z.add_multiple_of_row(i, c+1, -scale)
                            Z.add_multiple_of_column(c+1, i, scale)
                            if basis:
                                if trans:
                                    U.add_multiple_of_column(c+1, i, scale)
                                else:
                                    U.add_multiple_of_row(i, c+1, -scale)
                    # move to next column
                    # column of coefficients will cause zero search to fail
                    # at the top of this loop, and then we leave the column alone
                    # having built a block
                    c = c + 1

            # now have a companion matrix between rows  s  and  c
            # (inclusive), use it to clear entries to the right
            # but first record polynomial for block just built
            # this is the full monic polynomial, with correct coefficients
            p = []
            for i in range(s, c+1):
                p.append(-Z[i, c])
            p.append(R(1))
            polys.append(p)

            # use new unit columns (i) to clear rows to right (j)
            # all but top row of block to the right will become zero
            # it is important that the loop on  i  goes in reverse
            for i in range(c-1, s-1, -1):
                for j in range(c+1, n):
                    scale = Z[i+1, j]
                    # Effectively: Z.add_multiple_of_column(j, i, -scale)
                    Z.set_unsafe(i+1, j, zero)
                    # row j leads with zeros, up to, and including column c
                    # Effectively: Z.add_multiple_of_row(i, j, scale)
                    for k in range(c+1, n):
                        # Z[i,k] = Z[i,k] + scale*Z[j,k]
                        Z.set_unsafe(i, k, Z.get_unsafe(i,k)+scale*Z.get_unsafe(j,k))
                    if basis:
                        if trans:
                            U.add_multiple_of_column(j, i, -scale)
                        else:
                            U.add_multiple_of_row(i, j, scale)

            # row  s  to the right of the block needs care
            # we will be left with a leading one, or totally zeros
            nonzero = -1
            # locate a nonzero entry in row  s
            for j in range(c+1, n):
                if Z[s, j] != 0:
                    nonzero = j
                    break
            if (nonzero != -1):
                # swap column wih nonzero entry just outside block
                if nonzero != c+1:
                    Z.swap_columns(c+1, nonzero)
                    Z.swap_rows(c+1, nonzero)
                    if basis:
                        if trans:
                            U.swap_columns(c+1, nonzero)
                        else:
                            U.swap_rows(c+1, nonzero)
                scale = Z[s, c+1]
                Z.rescale_col(c+1, 1/scale)
                Z.rescale_row(c+1, scale)
                if basis:
                    if trans:
                        U.rescale_col(c+1, 1/scale)
                    else:
                        U.rescale_row(c+1, scale)
                for j in range(c+2, n):
                    scale = Z[s, j]
                    # exploiting leading zeros does not seem too beneficial here
                    Z.add_multiple_of_column(j, c+1, -scale)
                    Z.add_multiple_of_row(c+1, j, scale)
                    if basis:
                        if trans:
                            U.add_multiple_of_column(j, c+1, -scale)
                        else:
                            U.add_multiple_of_row(c+1, j, scale)

            # maybe move this closer to where we know which it is
            if c < n-1:
                corners.append(Z[s, c+1])
            # reset indices for constructing next block
            # transpose all elements, this is a zig or a zag
            c = c+1
            s = c
            Z = Z.transpose()
            if basis:
                U = U.transpose()
            trans = not trans

        # maybe return to un-transposed state
        if trans:
            Z = Z.transpose()
            if basis:
                U = U.transpose()
        if basis:
            # Storjohann computes U, we return its inverse
            # after this step, the columns of U are basis for representing self
            # code could be re-arranged to compute U's inverse directly
            return U.inverse(), Z, polys, corners
        return Z, polys, corners

    def zigzag_form(self, subdivide=True, transformation=False):
        r"""
        Find a matrix in ZigZag form that is similar to ``self``.

        INPUT:

        - ``self`` - a square matrix with entries from an exact field.

        - ``transformation`` - default: False - if ``True`` return a
          change-of-basis matrix relating the matrix and its ZigZag form.

        - ``subdivide`` - default: ``True`` - if ``True`` the ZigZag
          form matrix is subdivided according to the companion matrices
          described in the output section below.

        OUTPUT:

        A matrix in ZigZag form has blocks on the main diagonal that are
        companion matrices.  The first companion matrix has ones just below
        the main diagonal.  The last column has the negatives of coefficients
        of a monic polynomial, but not the leading one.  Low degree monomials
        have their coefficients in the earlier rows.  The second companion
        matrix is like the first only transposed.  The third is like the first.
        The fourth is like the second.  And so on.

        These blocks on the main diagonal define blocks just off the diagonal.
        To the right of the first companion matrix, and above the second
        companion matrix is a block that is totally zero, except the entry
        of the first row and first column may be a one.  Below the second
        block and to the left of the third block is a block that is totally
        zero, except the entry of the first row and first column may be one.
        This alternating pattern continues. It may now be apparent how this
        form gets its name.  Any other entry of the matrix is zero.  So this
        form is reminiscent of rational canonical form and is a good precursor
        to that form.

        If ``transformation`` is ``True``, then the output is a pair of
        matrices.  The first is the form ``Z`` and the second is an invertible
        matrix ``U`` such that ``U.inverse()*self*U`` equals ``Z``.  In other
        words, the repsentation of ``self`` with respect to the columns
        of ``U`` will be ``Z``.

        If subdivide is ``True`` then the matrix returned as the form is
        partitioned according to the companion matrices and these may be
        manipulated by several different matrix methods.

        For output that may be more useful as input to other routines,
        see the helper method :meth:`_zigzag_form`.

        .. note::

            An efffort has been made to optimize computation of the form,
            but no such work has been done for the computation of the
            transformation matrix, so for fastest results do not request the
            transformation matrix.

        ALGORITHM:

        ZigZag form, and its computation, are due to Arne Storjohann
        and are  described in [STORJOHANN-THESIS]_ and
        [STORJOHANN-ISACC98]_, where the former is more
        representative of the code here.

        EXAMPLES:

        Two examples that illustrate ZigZag form well. Notice that this is
        *not* a canonical form.  The two matrices below are similar, since
        they have equal Jordan canonical forms, yet their ZigZag forms are
        quite different.  In other words, while the computation of the form
        is deterministic, the final result, when viewed as a property of a
        linear transformation, is dependent on the basis used for the
        matrix representation.  ::

            sage: A = matrix(QQ, [[-68,   69, -27, -11, -65,   9, -181, -32],
            ...                   [-52,   52, -27,  -8, -52, -16, -133, -14],
            ...                   [ 92,  -97,  47,  14,  90,  32,  241,  18],
            ...                   [139, -144,  60,  18, 148, -10,  362,  77],
            ...                   [ 40,  -41,  12,   6,  45, -24,  105,  42],
            ...                   [-46,   48, -20,  -7, -47,   0, -122, -22],
            ...                   [-26,   27, -13,  -4, -29,  -6,  -66, -14],
            ...                   [-33,   34, -13,  -5, -35,   7,  -87, -23]])
            sage: Z, U = A.zigzag_form(transformation=True)
            sage: Z
            [  0   0   0  40|  1   0|  0   0]
            [  1   0   0  52|  0   0|  0   0]
            [  0   1   0  18|  0   0|  0   0]
            [  0   0   1  -1|  0   0|  0   0]
            [---------------+-------+-------]
            [  0   0   0   0|  0   1|  0   0]
            [  0   0   0   0|-25  10|  0   0]
            [---------------+-------+-------]
            [  0   0   0   0|  1   0|  0  -4]
            [  0   0   0   0|  0   0|  1  -4]
            sage: U.inverse()*A*U == Z
            True

            sage: B = matrix(QQ, [[ 16,  69, -13,   2, -52,  143,   90,  -3],
            ...                   [ 26,  54,   6,  -5, -28,   73,   73, -48],
            ...                   [-16, -79,  12, -10,  64, -142, -115,  41],
            ...                   [ 27,  -7,  21, -33,  39,  -20,  -42,  43],
            ...                   [  8, -75,  34, -32,  86, -156, -130,  42],
            ...                   [  2, -17,   7,  -8,  20,  -33,  -31,  16],
            ...                   [-24, -80,   7,  -3,  56, -136, -112,  42],
            ...                   [ -6, -19,   0,  -1,  13,  -28,  -27,  15]])
            sage: Z, U = B.zigzag_form(transformation=True)
            sage: Z
            [   0    0    0    0    0 1000|   0|   0]
            [   1    0    0    0    0  900|   0|   0]
            [   0    1    0    0    0  -30|   0|   0]
            [   0    0    1    0    0 -153|   0|   0]
            [   0    0    0    1    0    3|   0|   0]
            [   0    0    0    0    1    9|   0|   0]
            [-----------------------------+----+----]
            [   0    0    0    0    0    0|  -2|   0]
            [-----------------------------+----+----]
            [   0    0    0    0    0    0|   1|  -2]
            sage: U.inverse()*B*U == Z
            True

            sage: A.jordan_form() == B.jordan_form()
            True

        Two more examples, illustrating the two extremes of the zig-zag
        nature of this form.  The first has a one in each of the
        off-diagonal blocks, the second has all zeros in each off-diagonal
        block.  Notice again that the two matrices are similar, since their
        Jordan canonical forms are equal. ::

            sage: C = matrix(QQ, [[2,  31, -10,  -9, -125,  13,  62, -12],
            ...                   [0,  48, -16, -16, -188,  20,  92, -16],
            ...                   [0,   9,  -1,   2,  -33,   5,  18,   0],
            ...                   [0,  15,  -5,   0,  -59,   7,  30,  -4],
            ...                   [0, -21,   7,   2,   84, -10, -42,   5],
            ...                   [0, -42,  14,   8,  167, -17, -84,  13],
            ...                   [0, -50,  17,  10,  199, -23, -98,  14],
            ...                   [0,  15,  -5,  -2,  -59,   7,  30, -2]])
            sage: Z, U = C.zigzag_form(transformation=True)
            sage: Z
            [2|1|0|0|0|0|0|0]
            [-+-+-+-+-+-+-+-]
            [0|2|0|0|0|0|0|0]
            [-+-+-+-+-+-+-+-]
            [0|1|2|1|0|0|0|0]
            [-+-+-+-+-+-+-+-]
            [0|0|0|2|0|0|0|0]
            [-+-+-+-+-+-+-+-]
            [0|0|0|1|2|1|0|0]
            [-+-+-+-+-+-+-+-]
            [0|0|0|0|0|2|0|0]
            [-+-+-+-+-+-+-+-]
            [0|0|0|0|0|1|2|1]
            [-+-+-+-+-+-+-+-]
            [0|0|0|0|0|0|0|2]
            sage: U.inverse()*C*U == Z
            True

            sage: D = matrix(QQ, [[ -4,   3,    7,   2,  -4,   5,    7,   -3],
            ...                   [ -6,   5,    7,   2,  -4,   5,    7,   -3],
            ...                   [ 21, -12,   89,  25,   8,  27,   98,  -95],
            ...                   [ -9,   5,  -44, -11,  -3, -13,  -48,   47],
            ...                   [ 23, -13,   74,  21,  12,  22,   85,  -84],
            ...                   [ 31, -18,  135,  38,  12,  47,  155, -147],
            ...                   [-33,  19, -138, -39, -13, -45, -156,  151],
            ...                   [ -7,   4,  -29,  -8,  -3, -10,  -34,  34]])
            sage: Z, U = D.zigzag_form(transformation=True)
            sage: Z
            [ 0 -4| 0  0| 0  0| 0  0]
            [ 1  4| 0  0| 0  0| 0  0]
            [-----+-----+-----+-----]
            [ 0  0| 0  1| 0  0| 0  0]
            [ 0  0|-4  4| 0  0| 0  0]
            [-----+-----+-----+-----]
            [ 0  0| 0  0| 0 -4| 0  0]
            [ 0  0| 0  0| 1  4| 0  0]
            [-----+-----+-----+-----]
            [ 0  0| 0  0| 0  0| 0  1]
            [ 0  0| 0  0| 0  0|-4  4]
            sage: U.inverse()*D*U == Z
            True

            sage: C.jordan_form() == D.jordan_form()
            True

        ZigZag form is achieved entirely with the operations of the field, so
        while the eigenvalues may lie outside the field, this does not impede
        the computation of the form.  ::

            sage: F.<a> = GF(5^4)
            sage: A = matrix(F, [[     a,      0,  0, a + 3],
            ...                  [     0,a^2 + 1,  0,     0],
            ...                  [     0,      0,a^3,     0],
            ...                  [a^2 +4 ,     0,   0,a + 2]])
            sage: A.zigzag_form()
            [                    0 a^3 + 2*a^2 + 2*a + 2|                    0|                    0]
            [                    1               2*a + 2|                    0|                    0]
            [-------------------------------------------+---------------------+---------------------]
            [                    0                     0|                  a^3|                    0]
            [-------------------------------------------+---------------------+---------------------]
            [                    0                     0|                    0|              a^2 + 1]
            sage: A.eigenvalues()
            Traceback (most recent call last):
            ...
            NotImplementedError: eigenvalues() is not implemented for matrices with eigenvalues that are not in the fraction field of the base ring or in QQbar

        Subdivisions are optional.  ::

            sage: F.<a> = GF(5^4)
            sage: A = matrix(F, [[     a,      0,  0, a + 3],
            ...                  [     0,a^2 + 1,  0,     0],
            ...                  [     0,      0,a^3,     0],
            ...                  [a^2 +4 ,     0,   0,a + 2]])
            sage: A.zigzag_form(subdivide=False)
            [                    0 a^3 + 2*a^2 + 2*a + 2                     0                     0]
            [                    1               2*a + 2                     0                     0]
            [                    0                     0                   a^3                     0]
            [                    0                     0                     0               a^2 + 1]

        TESTS::

            sage: A = matrix(QQ, 2, 3, range(6))
            sage: A.zigzag_form()
            Traceback (most recent call last):
            ...
            TypeError: matrix must be square, not 2 x 3

            sage: A = matrix(Integers(6), 2, 2, range(4))
            sage: A.zigzag_form()
            Traceback (most recent call last):
            ...
            TypeError: matrix entries must come from an exact field, not Ring of integers modulo 6

            sage: A = matrix(RDF, 2, 2, range(4))
            sage: A.zigzag_form()
            Traceback (most recent call last):
            ...
            TypeError: matrix entries must come from an exact field, not Real Double Field

            sage: A = matrix(QQ, 2, range(4))
            sage: A.zigzag_form(transformation='junk')
            Traceback (most recent call last):
            ...
            ValueError: 'transformation' keyword must be True or False, not junk

            sage: A = matrix(QQ, 2, range(4))
            sage: A.zigzag_form(subdivide='garbage')
            Traceback (most recent call last):
            ...
            ValueError: 'subdivide' keyword must be True or False, not garbage

        .. rubric:: Citations

        .. [STORJOHANN-THESIS] A. Storjohann, Algorithms
           for Matrix Canonical Forms. PhD Thesis. Department
           of Computer Science, Swiss Federal Institute of
           Technology -- ETH, 2000.

        .. [STORJOHANN-ISACC98] A. Storjohann, An O(n^3)
           algorithm for Frobenius normal form. Proceedings
           of the International Symposium on Symbolic and
           Algebraic Computation (ISSAC'98), ACM Press, 1998, pp. 101-104.

        AUTHOR:

        - Rob Beezer (2011-06-09)
        """
        R = self.base_ring()
        if not self.is_square():
            raise TypeError("matrix must be square, not {0} x {1}".format(self.nrows(), self.ncols()))
        if not (R.is_field() and R.is_exact()):
            raise TypeError("matrix entries must come from an exact field, not {0}".format(R))
        if transformation not in [True, False]:
            raise ValueError("'transformation' keyword must be True or False, not {0}".format(transformation))
        if subdivide not in [True, False]:
            raise ValueError("'subdivide' keyword must be True or False, not {0}".format(subdivide))
        if transformation:
            U, Z, polys, corners = self._zigzag_form(basis=True)
        else:
            Z, polys, corners = self._zigzag_form(basis=False)
        if subdivide:
            s = 0
            splits = []
            for p in polys[:-1]:
                s = s + len(p) - 1
                splits.append(s)
            Z.subdivide(splits, splits)
        if transformation:
            return Z, U
        else:
            return Z

    def rational_form(self, format='right', subdivide=True):
        r"""
        Returns the rational canonical form, also known as Frobenius form.

        INPUT:

        - ``self`` - a square matrix with entries from an exact field.

        - ``format`` - default: 'right' - one of 'right', 'bottom',
          'left', 'top' or 'invariants'.  The first four will cause a
          matrix to be returned with companion matrices dictated by the
          keyword.  The value 'invariants' will cause a list of lists to
          be returned, where each list contains coefficients of a
          polynomial associated with a companion matrix.

        - ``subdivide`` - default: 'True' - if 'True' and a matrix is
          returned, then it contains subdivisions delineating the
          companion matrices along the diagonal.

        OUTPUT:

        The rational form of a matrix is a similar matrix composed of
        submatrices ("blocks") placed on the main diagonal.  Each block
        is a companion matrix.  Associated with each companion matrix is a
        polynomial.  In rational form, the polynomial of one block will
        divide the polynomial of the next block (and thus, the polynomials
        of all subsequent blocks).

        Rational form, also known as Frobenius form, is a canonical form.
        In other words, two matrices are similar if and only if their
        rational canonical forms are equal.  The algorithm used does not
        provide the similarity transformation matrix (also known as the
        change-of-basis matrix).

        Companion matrices may be written in one of four styles, and any
        such style may be selected with the ``format`` keyword.  See the
        companion matrix constructor,
        :meth:`sage.matrix.constructor.companion_matrix`,
        for more information about companion matrices.

        If the 'invariants' value is used for the ``format`` keyword,
        then the return value is a list of lists, where each list is the
        coefficients of the polynomial associated with one of the companion
        matrices on the diagonal.  These coefficients include the leading
        one of the monic polynomial and are ready to be coerced into
        any polynomial ring over the same field (see examples of this below).
        This return value is intended to be the most compact representation
        and the easiest to use for testing equality of rational forms.

        Because the minimal and characteristic polynomials of a companion
        matrix are the associated polynomial, it is easy to see that the
        product of the polynomials of the blocks will be the characteristic
        polynomial and the final polynomial will be the minimal polynomial
        of the entire matrix.

        ALGORITHM:

        We begin with ZigZag form, which is due to Arne Storjohann
        and is documented at :meth:`zigzag_form`.  Then we eliminate
        ''corner'' entries enroute to rational form via an additional
        algorithm of Storjohann's [STORJOHANN-EMAIL]_.

        EXAMPLES:

        The lists of coefficients returned with the ``invariants`` keyword
        are designed to easily convert to the polynomials associated with
        the companion matrices.  This is illustrated by the construction
        below of the ``polys`` list.  Then we can test the divisibility
        condition on the list of polynomials. Also the minimal and
        characteristic polynomials are easy to determine from this list. ::

            sage: A = matrix(QQ, [[ 11,  14, -15,  -4, -38, -29,  1,  23,  14, -63,  17,  24,  36,  32],
            ...                   [ 18,   6, -17, -11, -31, -43, 12,  26,   0, -69,  11,  13,  17,  24],
            ...                   [ 11,  16, -22,  -8, -48, -34,  0,  31,  16, -82,  26,  31,  39,  37],
            ...                   [ -8, -18,  22,  10,  46,  33,  3, -27, -12,  70, -19, -20, -42, -31],
            ...                   [-13, -21,  16,  10,  52,  43,  4, -28, -25,  89, -37, -20, -53, -62],
            ...                   [ -2,  -6,   0,   0,   6,  10,  1,   1,  -7,  14, -11,  -3, -10, -18],
            ...                   [ -9, -19,  -3,   4,  23,  30,  8,  -3, -27,  55, -40,  -5, -40, -69],
            ...                   [  4,  -8,  -1,  -1,   5,  -4,  9,   5, -11,   4, -14,  -2, -13, -17],
            ...                   [  1,  -2,  16,  -1,  19,  -2, -1, -17,   2,  19,   5, -25,  -7,  14],
            ...                   [  7,   7, -13,  -4, -26,  -21, 3,  18,   5, -40,   7,  15,  20,  14],
            ...                   [ -6,  -7, -12,   4,  -1,  18,  3,   8, -11,  15, -18,  17, -15, -41],
            ...                   [  5,  11, -11,  -3, -26, -19, -1,  14,  10, -42,  14,  17,  25,  23],
            ...                   [-16, -15,   3,  10,  29,  45, -1, -13, -19,  71, -35,  -2, -35, -65],
            ...                   [  4,   2,   3,  -2,  -2, -10,  1,   0,   3, -11,   6,  -4,   6,  17]])
            sage: A.rational_form()
            [   0   -4|   0    0    0    0|   0    0    0    0    0    0    0    0]
            [   1    4|   0    0    0    0|   0    0    0    0    0    0    0    0]
            [---------+-------------------+---------------------------------------]
            [   0    0|   0    0    0   12|   0    0    0    0    0    0    0    0]
            [   0    0|   1    0    0   -4|   0    0    0    0    0    0    0    0]
            [   0    0|   0    1    0   -9|   0    0    0    0    0    0    0    0]
            [   0    0|   0    0    1    6|   0    0    0    0    0    0    0    0]
            [---------+-------------------+---------------------------------------]
            [   0    0|   0    0    0    0|   0    0    0    0    0    0    0 -216]
            [   0    0|   0    0    0    0|   1    0    0    0    0    0    0  108]
            [   0    0|   0    0    0    0|   0    1    0    0    0    0    0  306]
            [   0    0|   0    0    0    0|   0    0    1    0    0    0    0 -271]
            [   0    0|   0    0    0    0|   0    0    0    1    0    0    0  -41]
            [   0    0|   0    0    0    0|   0    0    0    0    1    0    0  134]
            [   0    0|   0    0    0    0|   0    0    0    0    0    1    0  -64]
            [   0    0|   0    0    0    0|   0    0    0    0    0    0    1   13]

            sage: R = PolynomialRing(QQ, 'x')
            sage: invariants = A.rational_form(format='invariants')
            sage: invariants
            [[4, -4, 1], [-12, 4, 9, -6, 1], [216, -108, -306, 271, 41, -134, 64, -13, 1]]
            sage: polys = [R(p) for p in invariants]
            sage: [p.factor() for p in polys]
            [(x - 2)^2, (x - 3) * (x + 1) * (x - 2)^2, (x + 1)^2 * (x - 3)^3 * (x - 2)^3]
            sage: all(polys[i].divides(polys[i+1]) for i in range(len(polys)-1))
            True
            sage: polys[-1] == A.minimal_polynomial(var='x')
            True
            sage: prod(polys) == A.characteristic_polynomial(var='x')
            True

        Rational form is a canonical form.  Any two matrices are similar
        if and only if their rational forms are equal.  By starting with
        Jordan canonical forms, the matrices ``C`` and ``D`` below were
        built as similar matrices, while ``E`` was built to be just
        slightly different.  All three matrices have equal characteristic
        polynomials though ``E``'s minimal polynomial differs. ::

            sage: C = matrix(QQ, [[2,  31, -10,  -9, -125,  13,  62, -12],
            ...                   [0,  48, -16, -16, -188,  20,  92, -16],
            ...                   [0,   9,  -1,   2,  -33,   5,  18,   0],
            ...                   [0,  15,  -5,   0,  -59,   7,  30,  -4],
            ...                   [0, -21,   7,   2,   84, -10, -42,   5],
            ...                   [0, -42,  14,   8,  167, -17, -84,  13],
            ...                   [0, -50,  17,  10,  199, -23, -98,  14],
            ...                   [0,  15,  -5,  -2,  -59,   7,  30, -2]])
            sage: C.minimal_polynomial().factor()
            (x - 2)^2
            sage: C.characteristic_polynomial().factor()
            (x - 2)^8
            sage: C.rational_form()
            [ 0 -4| 0  0| 0  0| 0  0]
            [ 1  4| 0  0| 0  0| 0  0]
            [-----+-----+-----+-----]
            [ 0  0| 0 -4| 0  0| 0  0]
            [ 0  0| 1  4| 0  0| 0  0]
            [-----+-----+-----+-----]
            [ 0  0| 0  0| 0 -4| 0  0]
            [ 0  0| 0  0| 1  4| 0  0]
            [-----+-----+-----+-----]
            [ 0  0| 0  0| 0  0| 0 -4]
            [ 0  0| 0  0| 0  0| 1  4]

            sage: D = matrix(QQ, [[ -4,   3,    7,   2,  -4,   5,    7,   -3],
            ...                   [ -6,   5,    7,   2,  -4,   5,    7,   -3],
            ...                   [ 21, -12,   89,  25,   8,  27,   98,  -95],
            ...                   [ -9,   5,  -44, -11,  -3, -13,  -48,   47],
            ...                   [ 23, -13,   74,  21,  12,  22,   85,  -84],
            ...                   [ 31, -18,  135,  38,  12,  47,  155, -147],
            ...                   [-33,  19, -138, -39, -13, -45, -156,  151],
            ...                   [ -7,   4,  -29,  -8,  -3, -10,  -34,  34]])
            sage: D.minimal_polynomial().factor()
            (x - 2)^2
            sage: D.characteristic_polynomial().factor()
            (x - 2)^8
            sage: D.rational_form()
            [ 0 -4| 0  0| 0  0| 0  0]
            [ 1  4| 0  0| 0  0| 0  0]
            [-----+-----+-----+-----]
            [ 0  0| 0 -4| 0  0| 0  0]
            [ 0  0| 1  4| 0  0| 0  0]
            [-----+-----+-----+-----]
            [ 0  0| 0  0| 0 -4| 0  0]
            [ 0  0| 0  0| 1  4| 0  0]
            [-----+-----+-----+-----]
            [ 0  0| 0  0| 0  0| 0 -4]
            [ 0  0| 0  0| 0  0| 1  4]

            sage: E = matrix(QQ, [[ 0, -8,   4, -6, -2,   5, -3,  11],
            ...                   [-2, -4,   2, -4, -2,   4, -2,   6],
            ...                   [ 5, 14,  -7, 12,  3,  -8,  6, -27],
            ...                   [-3, -8,   7, -5,  0,   2, -6,  17],
            ...                   [ 0,  5,   0,  2,  4,  -4,  1,   2],
            ...                   [-3, -7,   5, -6, -1,   5, -4,  14],
            ...                   [ 6, 18, -10, 14,  4, -10, 10, -28],
            ...                   [-2, -6,   4, -5, -1,   3,  -3, 13]])
            sage: E.minimal_polynomial().factor()
            (x - 2)^3
            sage: E.characteristic_polynomial().factor()
            (x - 2)^8
            sage: E.rational_form()
            [  2|  0   0|  0   0|  0   0   0]
            [---+-------+-------+-----------]
            [  0|  0  -4|  0   0|  0   0   0]
            [  0|  1   4|  0   0|  0   0   0]
            [---+-------+-------+-----------]
            [  0|  0   0|  0  -4|  0   0   0]
            [  0|  0   0|  1   4|  0   0   0]
            [---+-------+-------+-----------]
            [  0|  0   0|  0   0|  0   0   8]
            [  0|  0   0|  0   0|  1   0 -12]
            [  0|  0   0|  0   0|  0   1   6]


        The principal feature of rational canonical form is that it
        can be computed over any field using only field operations.  Other
        forms, such as Jordan canonical form, are complicated by the need to
        determine the eigenvalues of the matrix, which can lie outside
        the field.  The following matrix has all of its eigenvalues
        outside the rationals - some are irrational (`\pm\sqrt{2}`) and
        the rest are complex (`-1\pm 2i`).  ::

            sage: A = matrix(QQ,
            ...   [[-154,  -3,  -54,   44,   48, -244,  -19,   67, -326,   85,   355,   581],
            ...    [ 504,  25,  156, -145, -171,  793,   99, -213, 1036, -247, -1152, -1865],
            ...    [ 294,  -1,  112,  -89,  -90,  469,   36, -128,  634, -160,  -695, -1126],
            ...    [ -49,  -32,  25,    7,   37,  -64,  -58,   12,  -42,  -14,    72,   106],
            ...    [-261, -123,  65,   47,  169, -358, -254,   70, -309,  -29,   454,   673],
            ...    [-448, -123, -10,  109,  227, -668, -262,  163, -721,   95,   896,  1410],
            ...    [  38,    7,   8,  -14,  -17,   66,    6,  -23,   73,  -29,   -78,  -143],
            ...    [ -96,   10, -55,   37,   24, -168,   17,   56, -231,   88,   237,   412],
            ...    [ 310,   67,  31,  -81, -143,  473,  143, -122,  538,  -98,  -641, -1029],
            ...    [ 139,  -35,  99,  -49,  -18,  236,  -41,  -70,  370, -118,  -377,  -619],
            ...    [ 243,    9,  81,  -72,  -81,  386,   43, -105,  508, -124,  -564,  -911],
            ...    [-155,   -3, -55,   45,   50, -245,  -27,   65, -328,   77,   365,  583]])
            sage: A.characteristic_polynomial().factor()
            (x^2 - 2)^2 * (x^2 + 2*x + 5)^4
            sage: A.eigenvalues(extend=False)
            []
            sage: A.rational_form()
            [  0  -5|  0   0   0   0|  0   0   0   0   0   0]
            [  1  -2|  0   0   0   0|  0   0   0   0   0   0]
            [-------+---------------+-----------------------]
            [  0   0|  0   0   0  10|  0   0   0   0   0   0]
            [  0   0|  1   0   0   4|  0   0   0   0   0   0]
            [  0   0|  0   1   0  -3|  0   0   0   0   0   0]
            [  0   0|  0   0   1  -2|  0   0   0   0   0   0]
            [-------+---------------+-----------------------]
            [  0   0|  0   0   0   0|  0   0   0   0   0  50]
            [  0   0|  0   0   0   0|  1   0   0   0   0  40]
            [  0   0|  0   0   0   0|  0   1   0   0   0   3]
            [  0   0|  0   0   0   0|  0   0   1   0   0 -12]
            [  0   0|  0   0   0   0|  0   0   0   1   0 -12]
            [  0   0|  0   0   0   0|  0   0   0   0   1  -4]
            sage: F.<x> = QQ[]
            sage: polys = A.rational_form(format='invariants')
            sage: [F(p).factor() for p in polys]
            [x^2 + 2*x + 5, (x^2 - 2) * (x^2 + 2*x + 5), (x^2 - 2) * (x^2 + 2*x + 5)^2]

        Rational form may be computed over any field.  The matrix below is
        an example where the eigenvalues lie outside the field.  ::

            sage: F.<a> = FiniteField(7^2)
            sage: A = matrix(F,
            ...   [[5*a + 3, 4*a + 1, 6*a + 2, 2*a + 5,       6, 4*a + 5, 4*a + 5,       5,   a + 6,      5,  4*a + 4],
            ...    [6*a + 3, 2*a + 4,       0,       6, 5*a + 5,     2*a, 5*a + 1,       1, 5*a + 2,     4*a, 5*a + 6],
            ...    [3*a + 1, 6*a + 6,   a + 6,       2,       0, 3*a + 6, 5*a + 4, 5*a + 6, 5*a + 2,       3, 4*a + 2],
            ...    [    3*a,     6*a,     3*a,     4*a, 4*a + 4, 3*a + 6,     6*a,       4, 3*a + 4, 6*a + 2,     4*a],
            ...    [4*a + 5,   a + 1, 4*a + 3, 6*a + 5, 5*a + 2, 5*a + 2,     6*a, 4*a + 6, 6*a + 4, 5*a + 3, 3*a + 1],
            ...    [    3*a,     6*a, 4*a + 1, 6*a + 2, 2*a + 5, 4*a + 6,       2,   a + 5, 2*a + 4, 2*a + 1, 2*a + 1],
            ...    [4*a + 5, 3*a + 3,       6, 4*a + 1, 4*a + 3, 6*a + 3,       6, 3*a + 3,       3,   a + 3,       0],
            ...    [6*a + 6,   a + 4, 2*a + 6, 3*a + 5, 4*a + 3,       2,       a, 3*a + 4,     5*a, 2*a + 5, 4*a + 3],
            ...    [3*a + 5, 6*a + 2,     4*a,   a + 5,       0,     5*a, 6*a + 5, 2*a + 1, 3*a + 1, 3*a + 5, 4*a + 2],
            ...    [3*a + 2,   a + 3, 3*a + 6,       a, 3*a + 5, 5*a + 1, 3*a + 2,   a + 3,   a + 2, 6*a + 1, 3*a + 3],
            ...    [6*a + 6, 5*a + 1,     4*a,       2, 5*a + 5, 3*a + 5, 3*a + 1,     2*a,     2*a, 2*a + 4, 4*a + 2]])
            sage: A.rational_form()
            [  a + 2|      0       0       0|      0       0       0       0       0       0       0]
            [-------+-----------------------+-------------------------------------------------------]
            [      0|      0       0   a + 6|      0       0       0       0       0       0       0]
            [      0|      1       0 6*a + 4|      0       0       0       0       0       0       0]
            [      0|      0       1 6*a + 4|      0       0       0       0       0       0       0]
            [-------+-----------------------+-------------------------------------------------------]
            [      0|      0       0       0|      0       0       0       0       0       0     2*a]
            [      0|      0       0       0|      1       0       0       0       0       0 6*a + 3]
            [      0|      0       0       0|      0       1       0       0       0       0 6*a + 1]
            [      0|      0       0       0|      0       0       1       0       0       0   a + 2]
            [      0|      0       0       0|      0       0       0       1       0       0   a + 6]
            [      0|      0       0       0|      0       0       0       0       1       0 2*a + 1]
            [      0|      0       0       0|      0       0       0       0       0       1 2*a + 1]
            sage: invariants = A.rational_form(format='invariants')
            sage: invariants
            [[6*a + 5, 1], [6*a + 1, a + 3, a + 3, 1], [5*a, a + 4, a + 6, 6*a + 5, 6*a + 1, 5*a + 6, 5*a + 6, 1]]
            sage: R.<x> = F[]
            sage: polys = [R(p) for p in invariants]
            sage: [p.factor() for p in polys]
            [x + 6*a + 5, (x + 6*a + 5) * (x^2 + (2*a + 5)*x + 5*a), (x + 6*a + 5) * (x^2 + (2*a + 5)*x + 5*a)^3]
            sage: polys[-1] == A.minimal_polynomial()
            True
            sage: prod(polys) == A.characteristic_polynomial()
            True
            sage: A.eigenvalues()
            Traceback (most recent call last):
            ...
            NotImplementedError: eigenvalues() is not implemented for matrices with eigenvalues that are not in the fraction field of the base ring or in QQbar

        Companion matrices may be selected as any one of four different types.
        See the documentation for the companion matrix constructor,
        :meth:`sage.matrix.constructor.companion_matrix`, for more information. ::

            sage: A = matrix(QQ, [[35, -18, -2, -45],
            ...                   [22, -22, 12, -16],
            ...                   [ 5, -12, 12,   4],
            ...                   [16,  -6, -4, -23]])
            sage: A.rational_form(format='right')
            [ 2| 0  0  0]
            [--+--------]
            [ 0| 0  0 10]
            [ 0| 1  0 -1]
            [ 0| 0  1  0]
            sage: A.rational_form(format='bottom')
            [ 2| 0  0  0]
            [--+--------]
            [ 0| 0  1  0]
            [ 0| 0  0  1]
            [ 0|10 -1  0]
            sage: A.rational_form(format='left')
            [ 2| 0  0  0]
            [--+--------]
            [ 0| 0  1  0]
            [ 0|-1  0  1]
            [ 0|10  0  0]
            sage: A.rational_form(format='top')
            [ 2| 0  0  0]
            [--+--------]
            [ 0| 0 -1 10]
            [ 0| 1  0  0]
            [ 0| 0  1  0]

        TESTS::

            sage: A = matrix(QQ, 2, 3, range(6))
            sage: A.rational_form()
            Traceback (most recent call last):
            ...
            TypeError: matrix must be square, not 2 x 3

            sage: A = matrix(Integers(6), 2, 2, range(4))
            sage: A.rational_form()
            Traceback (most recent call last):
            ...
            TypeError: matrix entries must come from an exact field, not Ring of integers modulo 6

            sage: A = matrix(RDF, 2, 2, range(4))
            sage: A.rational_form()
            Traceback (most recent call last):
            ...
            TypeError: matrix entries must come from an exact field, not Real Double Field

            sage: A = matrix(QQ, 2, range(4))
            sage: A.rational_form(format='junk')
            Traceback (most recent call last):
            ...
            ValueError: 'format' keyword must be 'right', 'bottom', 'left', 'top' or 'invariants', not junk

            sage: A = matrix(QQ, 2, range(4))
            sage: A.rational_form(subdivide='garbage')
            Traceback (most recent call last):
            ...
            ValueError: 'subdivide' keyword must be True or False, not garbage

        .. rubric:: Citations

        .. [STORJOHANN-EMAIL] A. Storjohann, Email Communication. 30 May 2011.

        AUTHOR:

        - Rob Beezer (2011-06-09)
        """
        from sage.rings.arith import gcd   # remove if translated to object-oriented calls
        import sage.rings.polynomial.polynomial_ring_constructor
        import sage.matrix.constructor

        R = self.base_ring()
        if not self.is_square():
            raise TypeError("matrix must be square, not {0} x {1}".format(self.nrows(), self.ncols()))
        if not (R.is_field() and R.is_exact()):
            raise TypeError("matrix entries must come from an exact field, not {0}".format(R))
        if format not in ['right', 'bottom', 'left', 'top', 'invariants']:
            raise ValueError("'format' keyword must be 'right', 'bottom', 'left', 'top' or 'invariants', not {0}".format(format))
        if subdivide not in [True, False]:
            raise ValueError("'subdivide' keyword must be True or False, not {0}".format(subdivide))

        _, polys, corners = self._zigzag_form(basis=False)
        k = len(polys)
        F = sage.rings.polynomial.polynomial_ring_constructor.PolynomialRing(R, 'x')
        C = [F(p) for p in polys]
        B = [(b == 1) for b in corners]
        B.append(False)  # no last block, so no corner

        if B[0]:
            V = [F(1)]
        else:
            V = [F(0)]

        for j in range(1, k):
            V[j-1] = gcd([V[j-1],C[j],C[j-1]])
            for i in range(j-2, -1, -1):
                V[i] = gcd([V[i], V[i+1], C[i]])
            m = F(-1)
            for i in range(j):
                g = gcd(m*V[i], C[i])
                q, _ = C[i].quo_rem(g)
                C[i] = g
                if B[j]:
                    _, V[i] = m.quo_rem(C[i])
                else:
                    V[i] = F(0)
                m = m * q
            C[j] = m * C[j]
            if B[j]:
                V.append(m)
            else:
                V.append(F(0))

        # Leading constant polynomials in C are size zero blocks, so toss them
        # Massage remainder to have leading coefficient 1
        while (len(C) > 0) and C[0].degree() == 0:
            C.remove(C[0])
        for i in range(len(C)):
            unit = C[i].list()[-1]
            if unit != R(1):
                C[i] = (1/unit)*C[i]

        if format == 'invariants':
            inv = []
            for poly in C:
                inv.append(poly.list())
            return inv
        elif format in ['right', 'left', 'top', 'bottom']:
            companions = []
            for poly in C:
                companions.append(sage.matrix.constructor.companion_matrix(poly, format=format))
            return sage.matrix.constructor.block_diagonal_matrix(companions, subdivide=subdivide)

    # A limited number of access-only properties are provided for matrices
    property T:
        r"""
        Returns the transpose of a matrix.

        EXAMPLE::

            sage: A = matrix(QQ, 5, range(25))
            sage: A.T
            [ 0  5 10 15 20]
            [ 1  6 11 16 21]
            [ 2  7 12 17 22]
            [ 3  8 13 18 23]
            [ 4  9 14 19 24]
        """
        def __get__(self):
            return self.transpose()

    property C:
        r"""
        Returns the conjugate matrix.

        EXAMPLE::

            sage: A = matrix(QQbar, [[     -3,  5 - 3*I, 7 - 4*I],
            ...                      [7 + 3*I, -1 + 6*I, 3 + 5*I],
            ...                      [3 + 3*I, -3 + 6*I, 5 +   I]])
            sage: A.C
            [      -3  5 + 3*I  7 + 4*I]
            [ 7 - 3*I -1 - 6*I  3 - 5*I]
            [ 3 - 3*I -3 - 6*I  5 - 1*I]

        """
        def __get__(self):
            return self.conjugate()

    property H:
        r"""
        Returns the conjugate-transpose (Hermitian) matrix.

        EXAMPLE::

            sage: A = matrix(QQbar, [[     -3,  5 - 3*I, 7 - 4*I],
            ...                      [7 + 3*I, -1 + 6*I, 3 + 5*I],
            ...                      [3 + 3*I, -3 + 6*I, 5 +   I]])
            sage: A.H
            [      -3  7 - 3*I  3 - 3*I]
            [ 5 + 3*I -1 - 6*I -3 - 6*I]
            [ 7 + 4*I  3 - 5*I  5 - 1*I]
        """
        def __get__(self):
            return self.conjugate().transpose()

    property I:
        r"""
        Returns the inverse of the matrix, if it exists.

        EXAMPLES::

            sage: A = matrix(QQ, [[-5, -3, -1, -7],
            ...                   [ 1,  1,  1,  0],
            ...                   [-1, -2, -2,  0],
            ...                   [-2, -1,  0, -4]])
            sage: A.I
            [ 0  2  1  0]
            [-4 -8 -2  7]
            [ 4  7  1 -7]
            [ 1  1  0 -2]

            sage: B = matrix(QQ, [[-11, -5, 18,  -6],
            ...                   [  1,  2, -6,   8],
            ...                   [ -4, -2,  7,  -3],
            ...                   [  1, -2,  5, -11]])
            sage: B.I
            Traceback (most recent call last):
            ...
            ZeroDivisionError: input matrix must be nonsingular
        """
        def __get__(self):
            return self.inverse()

    # end of Matrix class methods

def _smith_diag(d):
    r"""
    For internal use by the smith_form routine. Given a diagonal matrix d
    over a ring r, return matrices d', a,b such that a\*d\*b = d' and
    d' is diagonal with each entry dividing the next.

    If any of the d's is a unit, it replaces it with 1 (but no other
    attempt is made to pick "good" representatives of ideals).

    EXAMPLE::

        sage: from sage.matrix.matrix2 import _smith_diag
        sage: OE = EquationOrder(x^2 - x + 2, 'w')
        sage: A = matrix(OE, 2, [2,0,0,3])
        sage: D,U,V = _smith_diag(A); D,U,V
        (
        [1 0]  [2 1]  [ 1 -3]
        [0 6], [3 2], [-1  4]
        )
        sage: D == U*A*V
        True
        sage: m = matrix(GF(7),2, [3,0,0,6]); d,u,v = _smith_diag(m); d
        [1 0]
        [0 1]
        sage: u*m*v == d
        True
    """

    dp = d.__copy__()
    n = min(d.nrows(), d.ncols())
    R = d.base_ring()
    left = d.new_matrix(d.nrows(), d.nrows(), 1)
    right = d.new_matrix(d.ncols(), d.ncols(), 1)
    for i in xrange(n):
        I = R.ideal(dp[i,i])

        if I == R.unit_ideal():
            if dp[i,i] != 1:
                left.add_multiple_of_row(i,i,R(R(1)/(dp[i,i])) - 1)
                dp[i,i] = R(1)
            continue

        for j in xrange(i+1,n):
            if dp[j,j] not in I:
                t = R.ideal([dp[i,i], dp[j,j]]).gens_reduced()
                if len(t) > 1: raise ArithmeticError
                t = t[0]
                # find lambda, mu such that lambda*d[i,i] + mu*d[j,j] = t
                lamb = R(dp[i,i]/t).inverse_mod( R.ideal(dp[j,j]/t))
                mu = R((t - lamb*dp[i,i]) / dp[j,j])

                newlmat = dp.new_matrix(dp.nrows(), dp.nrows(), 1)
                newlmat[i,i] = lamb
                newlmat[i,j] = 1
                newlmat[j,i] = R(-dp[j,j]*mu/t)
                newlmat[j,j] = R(dp[i,i]/t)
                newrmat = dp.new_matrix(dp.ncols(), dp.ncols(), 1)
                newrmat[i,i] = 1
                newrmat[i,j] = R(-dp[j,j]/t)
                newrmat[j,i] = mu
                newrmat[j,j] = R(lamb*dp[i,i] / t)

                left = newlmat*left
                right = right*newrmat
                dp = newlmat*dp*newrmat
    return dp, left, right

def _generic_clear_column(m):
    r"""
    Reduce the first column of m to canonical form -- that is, all entries
    below the first are nonzero -- by multiplying on the left by invertible
    matrices over the given base ring.  This assumes that the base ring is a
    PID. Returns a pair (left, a) where left*self = a and a has first column in
    canonical form.

    If the first column is zero, then this function doesn't do anything very
    exciting.

    Used by the smith_form and hermite_form methods.

    EXAMPLES::

        sage: L.<w> = NumberField(x^2 - x + 2)
        sage: OL = L.ring_of_integers(); w = OL(w)
        sage: m = matrix(OL, 8, 4, [2*w - 2, 2*w + 1, -2, w, 2, -2,-2*w - 2, -2*w + 2, -w + 2, 2*w + 1, -w + 2, -w - 2, -2*w, 2*w, -w+ 2, w - 1, -2*w + 2, 2*w + 2, 2*w - 1, -w, 2*w + 2, -w + 2, 2, 2*w -1, w - 4, -2*w - 2, 2*w - 1, 0, 6, 7, 2*w + 1, 14])
        sage: s,t = m.echelon_form(transformation=True); t*m == s # indirect doctest
        True
        sage: s[0]
        (w, 0, 0, 0)
    """
    if m.nrows() <= 1 or m.ncols() <= 0:
        return m.new_matrix(m.nrows(), m.nrows(), 1), m

    a = m.__copy__()
    left_mat = m.new_matrix(m.nrows(), m.nrows(), 1)
    R = m.base_ring()

    # case 1: if a_{0, 0} = 0 and a_{k, 0} != 0 for some k, swap rows 0 and k.
    if a[0, 0] == 0:
        k = 0
        while a[k, 0] == 0:
            k += 1
            if k == a.nrows(): # first column is zero
                return left_mat, a
        # k is now first row such that a[k, 0] is nonzero
        left_mat[0,0] = 0
        left_mat[k,k] = 0
        left_mat[0,k] = 1
        left_mat[k,0] = -1
        a = left_mat*a
        if left_mat * m != a:
            raise ArithmeticError, "Something went wrong"

    # case 2: if there is an entry at position (k,j) such that a_{0,j}
    # does not divide a_{k,j}, then we know that there are c,d in R such
    # that c.a_{0,j} - d.a_{k,j} = B where B = gcd(a_{0,j}, a_{k,j}) (which
    # is well-defined since R is a PID).
    # Then for appropriate e,f the matrix
    # [c,-d]
    # [e,f]
    # is invertible over R

    if a[0,0] != 0:
        I = R.ideal(a[0, 0]) # need to make sure we change this when a[0,0] changes
    else:
        I = R.zero_ideal()
    for k in xrange(1, a.nrows()):
        if a[k,0] not in I:
            try:
                v = R.ideal(a[0,0], a[k,0]).gens_reduced()
            except Exception, msg:
                raise ArithmeticError, "%s\nCan't create ideal on %s and %s" % (msg, a[0,0], a[k,0])
            if len(v) > 1:
                raise ArithmeticError, "Ideal %s not principal" %  R.ideal(a[0,0], a[k,0])
            B = v[0]

            # now we find c,d, using the fact that c * (a_{0,0}/B) - d *
            # (a_{k,0}/B) = 1, so c is the inverse of a_{0,0}/B modulo
            # a_{k,0}/B.
            # need to handle carefully the case when a_{k,0}/B is a unit, i.e. a_{k,0} divides
            # a_{0,0}.

            c = R(a[0,0] / B).inverse_mod(R.ideal(a[k,0] / B))
            d = R( (c*a[0,0] - B)/(a[k,0]) )

            # sanity check
            if c*a[0,0] - d*a[k,0] != B:
                raise ArithmeticError

            # now we find e,f such that e*d + c*f = 1 in the same way
            if c != 0:
                e = d.inverse_mod( R.ideal(c) )
                f = R((1 - d*e)/c)
            else:
                e = R(-a[k,0]/B) # here d is a unit and this is just 1/d
                f = R(1)

            if e*d + c*f != 1:
                raise ArithmeticError
            newlmat = left_mat.parent()(1)
            newlmat[0,0] = c
            newlmat[0,k] = -d
            newlmat[k,0] = e
            newlmat[k,k] = f
            if newlmat.det() != 1:
                raise ArithmeticError
            a = newlmat*a
            I = R.ideal(a[0,0])
            left_mat = newlmat*left_mat
            if left_mat * m != a:
                raise ArithmeticError

    # now everything in column 0 is divisible by the pivot
    for i in xrange(1,a.nrows()):
        s = R( a[i, 0]/a[0, 0])
        a.add_multiple_of_row(i, 0, -s )
        left_mat.add_multiple_of_row(i, 0, -s)
    if left_mat * m != a:
        raise ArithmeticError

    return left_mat, a

def _smith_onestep(m):
    r"""
    Carry out one step of Smith normal form for matrix m. Returns three matrices a,b,c over
    the same base ring as m, such that a \* m \* c = b, a and c have
    determinant 1, and the zeroth row and column of b have no nonzero
    entries except b[0,0].

    EXAMPLE::

        sage: from sage.matrix.matrix2 import _smith_onestep
        sage: OE = NumberField(x^2 - x + 2,'w').ring_of_integers()
        sage: w = OE.ring_generators()[0]
        sage: m = matrix(OE, 3,3,[1,0,7,2,w, w+17, 13+8*w, 0, 6])
        sage: a,b,c = _smith_onestep(m); b
        [         1          0          0]
        [         0          w      w + 3]
        [         0          0 -56*w - 85]
        sage: a * m * c == b
        True
    """

    a = m.__copy__()
    left_mat = m.new_matrix(m.nrows(), m.nrows(), 1)
    right_mat = m.new_matrix(m.ncols(), m.ncols(), 1)

    if m == 0 or (m.nrows() <= 1 and m.ncols() <= 1):
        return left_mat, m, right_mat

    # preparation: if column 0 is zero, swap it with the first nonzero column
    j = 0
    while a.column(j) == 0: j += 1
    if j > 0:
        right_mat[0,0] = right_mat[j,j] = 0
        right_mat[0,j] = 1
        right_mat[j,0] = -1
        a = a*right_mat
        if m * right_mat != a:
            raise ArithmeticError

    left_mat, a = _generic_clear_column(a)
    assert left_mat * m * right_mat == a

    # test if everything to the right of the pivot in row 0 is good as well
    isdone = True
    for jj in xrange(j+1, a.ncols()):
        if a[0,jj] != 0:
            isdone = False

    # if not we recurse -- algorithm must terminate if R is Noetherian.
    if isdone == False:
        s,t,u = _smith_onestep(a.transpose())
        left_mat = u.transpose() * left_mat
        a = t.transpose()
        right_mat = right_mat* s.transpose()

    return left_mat, a, right_mat

def _dim_cmp(x,y):
    """
    Used internally by matrix functions. Given 2-tuples (x,y), returns
    their comparison based on the first component.

    EXAMPLES::

        sage: from sage.matrix.matrix2 import _dim_cmp
        sage: V = [(QQ^3, 2), (QQ^2, 1)]
        sage: _dim_cmp(V[0], V[1])
        1
    """
    return cmp(x[0].dimension(), y[0].dimension())

def decomp_seq(v):
    """
    This function is used internally be the decomposition matrix
    method. It takes a list of tuples and produces a sequence that is
    correctly sorted and prints with carriage returns.

    EXAMPLES::

        sage: from sage.matrix.matrix2 import decomp_seq
        sage: V = [(QQ^3, 2), (QQ^2, 1)]
        sage: decomp_seq(V)
        [
        (Vector space of dimension 2 over Rational Field, 1),
        (Vector space of dimension 3 over Rational Field, 2)
        ]
    """
    list.sort(v, _dim_cmp)
    return Sequence(v, universe=tuple, check=False, cr=True)


def cmp_pivots(x,y):
    """
    Compare two sequences of pivot columns.

    - If x is shorter than y, return -1, i.e., x < y, "not as good".

    - If x is longer than y, x > y, "better".

    - If the length is the same then x is better, i.e., x > y if the
      entries of x are correspondingly >= those of y with one being
      greater.
    """
    if len(x) < len(y):
        return -1
    if len(x) > len(y):
        return 1
    if x < y:
        return 1
    elif x == y:
        return 0
    else:
        return -1


def _choose(Py_ssize_t n, Py_ssize_t t):
    """
    Returns all possible sublists of length t from range(n)

    Based on algorithm T from Knuth's taocp part 4: 7.2.1.3 p.5 This
    function replaces the one based on algorithm L because it is
    faster.

    EXAMPLES::

        sage: from sage.matrix.matrix2 import _choose
        sage: _choose(1,1)
        [[0]]
        sage: _choose(4,1)
        [[0], [1], [2], [3]]
        sage: _choose(4,4)
        [[0, 1, 2, 3]]

    AUTHORS:

    - Jaap Spies (2007-11-14)
    """
    cdef Py_ssize_t j, temp

    x = []               # initialize T1
    c = range(t)
    if t == n:
        x.append(c)
        return x
    c.append(n)
    c.append(0)
    j = t-1

    while True:
        x.append(c[:t])    # visit T2
        if j >= 0:
            c[j] = j+1
            j = j-1
            continue       # goto T2

        if c[0]+1 < c[1]:  # T3 easy case!
            c[0] = c[0]+1
            continue
        else:
            j = 1

        while True:
            c[j-1] = j-1      # T4 find j
            temp = c[j]+1
            if temp == c[j+1]:
                j = j+1
            else:
                break


        if j >= t:     # T5 stop?
            break

        c[j] = temp    # T6
        j = j-1

    return x


def _binomial(Py_ssize_t n, Py_ssize_t k):
    """
    Fast and unchecked implementation of binomial(n,k) This is only for
    internal use.

    EXAMPLES::

        sage: from sage.matrix.matrix2 import _binomial
        sage: _binomial(10,2)
        45
        sage: _binomial(10,5)
        252

    AUTHORS:

    - Jaap Spies (2007-10-26)
    """
    cdef Py_ssize_t i

    if k > (n/2):
        k = n-k
    if k == 0:
        return 1

    result = n
    n, k = n-1, k-1
    i = 2
    while k > 0:
        result = (result*n)/i
        i, n, k = i+1, n-1, k-1
    return result

def _jordan_form_vector_in_difference(V, W):
    r"""
    Given two lists of vectors ``V`` and ``W`` over the same base field,
    returns a vector in the difference ``V - W``.  If the difference is
    empty, returns ``None``.

    NOTES:

    This is meant to be a private helper method for the ``jordan_form`` method
    in the above class.

    TEST::

        sage: v = vector(ZZ, [1,0,0,0])
        sage: w = vector(ZZ, [0,1,0,0])
        sage: u = vector(ZZ, [1,1,0,0])
        sage: sage.matrix.matrix2._jordan_form_vector_in_difference([v,w], [u])
        (1, 0, 0, 0)
    """
    if len(V) == 0:
        return None
    if len(W) == 0:
        return V[0]
    W_space = sage.all.span(W)
    for v in V:
        if v not in W_space:
            return v
    return None

