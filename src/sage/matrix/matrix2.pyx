r"""
Base class for matrices, part 2

For design documentation see matrix/docs.py.

AUTHORS:

- William Stein: initial version

- Jaap Spies (2006-02-24): added ``prod_of_row_sums``, ``permanent``,
  ``permanental_minor``, ``rook_vector`` methods

- Robert Bradshaw (2007-06-14): added ``subdivide`` method

- Jaap Spies (2007-11-14): implemented ``_binomial``, ``_choose`` auxiliary functions

- William Stein (2007-11-18): added ``_gram_schmidt_noscale`` method

- David Loeffler (2008-12-05): added ``smith_form`` method

- David Loeffler (2009-06-01): added ``_echelon_form_PID`` method

- Sebastian Pancratz (2009-06-25): implemented ``adjoint`` and ``charpoly``
  methods; fixed ``adjoint`` reflecting the change that ``_adjoint`` is now
  implemented in :class:`Matrix`; used the division-free algorithm for
  ``charpoly``

- Rob Beezer (2009-07-13): added ``elementwise_product`` method

- Miguel Marco (2010-06-19): modified eigenvalues and eigenvectors functions to
  allow the option ``extend=False``

- Thierry Monteil (2010-10-05): bugfix for :issue:`10063`, so that the
  determinant is computed even for rings for which the ``is_field`` method is not
  implemented.

- Rob Beezer (2010-12-13): added ``conjugate_transpose`` method

- Rob Beezer (2011-02-05): refactored all of the matrix kernel routines; added
  ``extended_echelon_form``, ``right_kernel_matrix``, ``QR``,
  ``_gram_schmidt_noscale``, ``is_similar`` methods

- Moritz Minzlaff (2011-03-17): corrected ``_echelon_form_PID`` method for
  matrices of one row, fixed in :issue:`9053`

- Rob Beezer (2011-06-09): added ``is_normal``, ``is_diagonalizable``, ``LU``,
  ``cyclic_subspace``, ``zigzag_form``, ``rational_form`` methods

- Rob Beezer (2012-05-27): added ``indefinite_factorization``,
  ``is_positive_definite``, ``cholesky`` methods

- Darij Grinberg (2013-10-01): added first (slow) pfaffian implementation

- Mario Pernici (2014-07-01): modified ``rook_vector`` method

- Rob Beezer (2015-05-25): modified ``is_similar`` method

- Samuel Lelièvre (2020-09-18): improved method ``LLL_gram`` based on a patch
  by William Stein posted at :issue:`5178`, moving the method from its initial
  location in ``sage.matrix.integer_matrix_dense``

- Michael Jung (2020-10-02): added Bär-Faddeev-LeVerrier algorithm for the
  Pfaffian

- Moritz Firsching(2020-10-05): added ``quantum_determinant``

- Dima Pasechnik (2022-11-08): fixed ``echelonize`` for inexact matrices
"""

# ****************************************************************************
#       Copyright (C) 2005, 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from copy import copy

from cpython cimport *
from cysignals.signals cimport sig_check

from sage.categories.fields import Fields
from sage.categories.integral_domains import IntegralDomains
from sage.categories.rings import Rings
from sage.misc.lazy_string import lazy_string
from sage.misc.randstate cimport current_randstate
from sage.structure.coerce cimport py_scalar_parent
from sage.structure.sequence import Sequence
from sage.structure.coerce cimport coercion_model
from sage.structure.element import Vector
from sage.structure.element cimport have_same_parent
from sage.misc.verbose import verbose
from sage.rings.number_field.number_field_base import NumberField
from sage.rings.integer_ring import ZZ, IntegerRing_class
from sage.rings.rational_field import QQ, RationalField
import sage.rings.abc
from sage.arith.numerical_approx cimport digits_to_bits

import sage.modules.free_module
from sage.matrix import berlekamp_massey
from sage.modules.free_module_element import FreeModuleElement
from sage.matrix.matrix_misc import permanental_minor_polynomial

from sage.misc.misc_c import prod

# used to deprecate only adjoint method
from sage.misc.superseded import deprecation, deprecated_function_alias


# temporary hack to silence the warnings from #34806
def ideal_or_fractional(R, *args):
    if isinstance(R, sage.rings.abc.Order):
        R = R.number_field()
    return R.ideal(*args)


_Fields = Fields()


cdef class Matrix(Matrix1):
    """
    Base class for matrices, part 2.

    TESTS::

        sage: m = matrix(ZZ['x'], 2, 3, [1..6])
        sage: TestSuite(m).run()

    Check that a pair consisting of a matrix and its echelon form is
    pickled correctly (this used to give a wrong answer due to a Python
    bug, see :issue:`17527`)::

        sage: K.<x> = FractionField(QQ['x'])
        sage: m = Matrix([[1], [x]])
        sage: t = (m, m.echelon_form())
        sage: loads(dumps(t))
        (
        [1]  [1]
        [x], [0]
        )
    """
    def _backslash_(self, B):
        r"""
        Used to compute `A \backslash B`, i.e., the backslash solver
        operator.

        DEPRECATED

        EXAMPLES::

            sage: A = matrix(QQ, 3, [1,2,4,2,3,1,0,1,2])
            sage: B = matrix(QQ, 3, 2, [1,7,5, 2,1,3])
            sage: C = A._backslash_(B); C
            doctest:...: DeprecationWarning: the backslash operator has been deprecated; use A.solve_right(B) instead
            See https://github.com/sagemath/sage/issues/36394 for details.
            [  -1    1]
            [13/5 -3/5]
            [-4/5  9/5]
            sage: A*C == B
            True
            sage: A._backslash_(B) == A \ B
            doctest:...: DeprecationWarning: the backslash operator has been deprecated; use A.solve_right(B) instead
            See https://github.com/sagemath/sage/issues/36394 for details.
            True
            sage: A._backslash_(B) == A.solve_right(B)
            True
        """
        deprecation(36394, 'the backslash operator has been deprecated; use A.solve_right(B) instead')
        return self.solve_right(B)

    def subs(self, *args, **kwds):
        """
        Substitute values to the variables in that matrix.

        All the arguments are transmitted unchanged to the method ``subs`` of
        the coefficients.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: var('a,b,d,e')
            (a, b, d, e)
            sage: m = matrix([[a,b], [d,e]])
            sage: m.substitute(a=1)
            [1 b]
            [d e]
            sage: m.subs(a=b, b=d)
            [b d]
            [d e]
            sage: m.subs({a: 3, b: 2, d: 1, e: -1})
            [ 3  2]
            [ 1 -1]

        The parent of the newly created matrix might be different from the
        initial one. It depends on what the method ``.subs`` does on
        coefficients (see :issue:`19045`)::

            sage: x = polygen(ZZ)
            sage: m = matrix([[x]])
            sage: m2 = m.subs(x=2)
            sage: m2.parent()
            Full MatrixSpace of 1 by 1 dense matrices over Integer Ring
            sage: m1 = m.subs(x=RDF(1))
            sage: m1.parent()
            Full MatrixSpace of 1 by 1 dense matrices over Real Double Field

        However, sparse matrices remain sparse::

            sage: m = matrix({(3,2): -x, (59,38): x^2 + 2}, nrows=1000, ncols=1000)
            sage: m1 = m.subs(x=1)
            sage: m1.is_sparse()
            True
        """
        from sage.matrix.constructor import matrix
        if self.is_sparse():
            return matrix({ij: self[ij].subs(*args, **kwds) for ij in self.nonzero_positions()},
                          nrows=self._nrows, ncols=self._ncols, sparse=True)
        else:
            return matrix([a.subs(*args, **kwds) for a in self.list()],
                          nrows=self._nrows, ncols=self._ncols, sparse=False)

    def solve_left(self, B, check=True, *, extend=True):
        r"""
        Try to find a solution `X` to the equation `X A = B`.

        If ``self`` is a matrix `A`, then this function returns a
        vector or matrix `X` such that `X A = B`. If
        `B` is a vector then `X` is a vector and if
        `B` is a matrix, then `X` is a matrix.

        Over inexact rings, the output of this function may not be an
        exact solution. For example, over the real or complex double
        field, this method computes a least-squares solution if the
        system is not square.

        INPUT:

        - ``B`` -- a matrix or vector

        - ``extend`` -- boolean (default: ``True``); when set to ``True``,
          some solvers will return solutions over a larger ring than the
          base ring of the inputs (a typical case are rational solutions
          for integer linear systems). When set to ``False``, a solution
          over the base ring is returned, with a :exc:`ValueError`
          being raised if none exists.

        - ``check`` -- boolean (default: ``True``); verify the answer
          if the system is non-square or rank-deficient, and if its
          entries lie in an exact ring. Meaningless over most inexact
          rings, or when the system is square and of full rank.

        OUTPUT:

        If the system is square and has full rank, the unique solution
        is returned, and no check is done on the answer. Over inexact
        rings, you should expect this answer to be inexact.
        Moreover, due to the numerical issues involved, an error
        may be thrown in this case -- specifically if the system is
        singular but if SageMath fails to notice that.

        If the system is not square or does not have full rank, then a
        solution is attempted via other means. For example, over
        ``RDF`` or ``CDF`` a least-squares solution is returned, as
        with MATLAB's "backslash" operator. For most inexact rings, the
        ``check`` parameter is ignored because an approximate solution
        will be returned in any case. Over exact rings, on the other
        hand, setting the ``check`` parameter results in an additional
        test to determine whether or not the answer actually solves the
        system exactly. If a symbolic system involves only exact elements,
        its solution can still be checked.

        If `B` is a vector, the result is returned as a vector, as well,
        and as a matrix, otherwise.

        .. SEEALSO::

            :meth:`solve_right`

        EXAMPLES::

            sage: A = matrix(QQ, 4,2, [0, -1, 1, 0, -2, 2, 1, 0])
            sage: B = matrix(QQ, 2,2, [1, 0, 1, -1])
            sage: X = A.solve_left(B)
            sage: X*A == B
            True
            sage: X == B / A
            True

        ::

            sage: A = matrix([(3, -1, 0, 0), (1, 1, -2, 0), (0, 0, 0, -3)])
            sage: B = matrix(QQ, 3, 1, [0, 0, -1])
            sage: A.solve_left(B)
            Traceback (most recent call last):
            ...
            ValueError: number of columns of self must equal number of columns
            of right-hand side

        Over the reals::

            sage: A = matrix(RDF, 3,3, [1,2,5,7.6,2.3,1,1,2,-1]); A
            [ 1.0  2.0  5.0]
            [ 7.6  2.3  1.0]
            [ 1.0  2.0 -1.0]
            sage: b = vector(RDF,[1,2,3])
            sage: x = A.solve_left(b); x.zero_at(2e-17)  # fix noisy zeroes             # needs scipy
            (0.666666666..., 0.0, 0.333333333...)
            sage: x.parent()                                                            # needs scipy
            Vector space of dimension 3 over Real Double Field
            sage: x*A  # tol 1e-14                                                      # needs scipy
            (0.9999999999999999, 1.9999999999999998, 3.0)

        Over the complex numbers::

            sage: # needs sage.rings.complex_double sage.symbolic
            sage: A = matrix(CDF, [[      0, -1 + 2*I,  1 - 3*I,        I],
            ....:                  [2 + 4*I, -2 + 3*I, -1 + 2*I,   -1 - I],
            ....:                  [  2 + I,    1 - I,       -1,        5],
            ....:                  [    3*I,   -1 - I,   -1 + I,   -3 + I]])
            sage: b = vector(CDF, [2 -3*I, 3, -2 + 3*I, 8])
            sage: x = A.solve_left(b); x                                                # needs scipy
            (-1.55765124... - 0.644483985...*I, 0.183274021... + 0.286476868...*I,
             0.270818505... + 0.246619217...*I, -1.69003558... - 0.828113879...*I)
            sage: x.parent()                                                            # needs scipy
            Vector space of dimension 4 over Complex Double Field
            sage: abs(x*A - b) < 1e-14                                                  # needs scipy
            True

        If ``b`` is given as a matrix, the result will be a matrix, as well::

            sage: A = matrix(RDF, 3, 3, [2, 5, 0, 7, 7, -2, -4.3, 0, 1])
            sage: b = matrix(RDF, 2, 3, [2, -4, -5, 1, 1, 0.1])
            sage: A.solve_left(b)  # tol 1e-14                                          # needs scipy
            [  -6.495454545454545    4.068181818181818   3.1363636363636354]
            [  0.5277272727272727  -0.2340909090909091 -0.36818181818181817]

        If `A` is a non-square matrix, the result is a least-squares solution.
        For a tall matrix, this may give a solution with a least-squares error
        of almost zero::

            sage: A = matrix(RDF, 3, 2, [1, 3, 4, 2, 0, -3])
            sage: b = vector(RDF, [5, 6])
            sage: x = A.solve_left(b)                                                   # needs scipy
            sage: (x * A - b).norm() < 1e-14                                            # needs scipy
            True

        For a wide matrix `A`, the error is usually not small::

            sage: A = matrix(RDF, 2, 3, [1, 3, 4, 2, 0, -3])
            sage: b = vector(RDF, [5, 6, 1])
            sage: x = A.solve_left(b)                                                   # needs scipy
            sage: (x * A - b).norm()  # tol 1e-14                                       # needs scipy
            0.9723055853282466

        TESTS::

            sage: A = matrix(QQ, 4,2, [0, -1, 1, 0, -2, 2, 1, 0])
            sage: B = vector(QQ, 2, [2,1])
            sage: X = A.solve_left(B)
            sage: X*A == B
            True
            sage: X
            (-1, 2, 0, 0)

            sage: # needs sage.libs.pari
            sage: A = Matrix(Zmod(128), 2, 3, [5, 29, 33, 64, 0, 7])
            sage: B = vector(Zmod(128), [31,39,56])
            sage: X = A.solve_left(B); X
            (19, 83)
            sage: X * A == B
            True

            sage: M = matrix([(3,-1,0,0),(1,1,-2,0),(0,0,0,-3)])
            sage: B = matrix(QQ,3,1, [0,0,-1])
            sage: M.solve_left(B)
            Traceback (most recent call last):
            ...
            ValueError: number of columns of self must equal number of columns
            of right-hand side

        A degenerate case::

            sage: A = matrix(RDF, 0, 0, [])
            sage: A.solve_left(vector(RDF,[]))                                          # needs scipy
            ()

        Over an inexact ring like ``RDF``, the coefficient matrix of a
        square system must be nonsingular::

            sage: A = matrix(RDF, 5, range(25))
            sage: b = vector(RDF, [1,2,3,4,5])
            sage: A.solve_left(b)                                                       # needs scipy
            Traceback (most recent call last):
            ...
            LinAlgError: Matrix is singular.

        The vector of constants needs the correct degree::

            sage: A = matrix(RDF, 5, range(25))
            sage: b = vector(RDF, [1,2,3,4])
            sage: A.solve_left(b)                                                       # needs scipy
            Traceback (most recent call last):
            ...
            ValueError: number of columns of self must equal degree of
            right-hand side

        The vector of constants needs to be compatible with
        the base ring of the coefficient matrix::

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = FiniteField(27)
            sage: b = vector(F, [a,a,a,a,a])
            sage: A.solve_left(b)
            Traceback (most recent call last):
            ...
            TypeError: no common canonical parent for objects with parents: ...

        Check that coercions work correctly (:issue:`17405`)::

            sage: # needs sage.rings.complex_double sage.symbolic
            sage: A = matrix(RDF, 2, range(4))
            sage: b = vector(CDF, [1+I, 2])
            sage: A.solve_left(b)                                                       # needs scipy
            (0.5 - 1.5*I, 0.5 + 0.5*I)
            sage: b = vector(QQ[I], [1+I, 2])
            sage: x = A.solve_left(b)                                                   # needs scipy

        Over the inexact ring ``SR``, we can still verify the solution
        if all of the elements involved were exact to begin with; if
        any are inexact, however, the ``check`` is still skipped
        (:issue:`29729` and :issue:`33159`)::

            sage: A = matrix(SR, [[1, 1]])                                              # needs sage.symbolic
            sage: b = vector(SR, [2, 3])                                                # needs sage.symbolic
            sage: A.solve_left(b)                                                       # needs sage.symbolic
            Traceback (most recent call last):
            ...
            ValueError: matrix equation has no solutions

        In this case, turning off the ``check`` leads to a wrong result::

            sage: A.solve_left(b, check=False)                                          # needs sage.symbolic
            (2)
        """
        if isinstance(B, Vector):
            try:
                return self.transpose().solve_right(B, check=check, extend=extend)
            except ValueError as e:
                raise e.__class__(str(e).replace('row', 'column'))
        else:
            try:
                return self.transpose().solve_right(B.transpose(), check=check, extend=extend).transpose()
            except ValueError as e:
                raise e.__class__(str(e).replace('row', 'column'))

    def solve_right(self, B, check=True, *, extend=True):
        r"""
        Try to find a solution `X` to the equation `A X = B`.

        If ``self`` is a matrix `A`, then this function returns a
        vector or matrix `X` such that `A X = B`. If
        `B` is a vector then `X` is a vector and if
        `B` is a matrix, then `X` is a matrix.

        Over inexact rings, the output of this function may not be an
        exact solution. For example, over the real or complex double
        field, this method computes a least-squares solution if the
        system is not square.

        INPUT:

        - ``B`` -- a matrix or vector

        - ``extend`` -- boolean (default: ``True``); when set to ``True``,
          some solvers will return solutions over a larger ring than the
          base ring of the inputs (a typical case are rational solutions
          for integer linear systems). When set to ``False``, a solution
          over the base ring is returned, with a :exc:`ValueError`
          being raised if none exists.

        - ``check`` -- boolean (default: ``True``); verify the answer
          if the system is non-square or rank-deficient, and if its
          entries lie in an exact ring. Meaningless over most inexact
          rings, or when the system is square and of full rank.

        OUTPUT:

        If the system is square and has full rank, the unique solution
        is returned, and no check is done on the answer. Over inexact
        rings, you should expect this answer to be inexact.
        Moreover, due to the numerical issues involved, an error
        may be thrown in this case -- specifically if the system is
        singular but if SageMath fails to notice that.

        If the system is not square or does not have full rank, then a
        solution is attempted via other means. For example, over
        ``RDF`` or ``CDF`` a least-squares solution is returned, as
        with MATLAB's "backslash" operator. For most inexact rings, the
        ``check`` parameter is ignored because an approximate solution
        will be returned in any case. Over exact rings, on the other
        hand, setting the ``check`` parameter results in an additional
        test to determine whether or not the answer actually solves the
        system exactly. If a symbolic system involves only exact elements,
        its solution can still be checked.

        If `B` is a vector, the result is returned as a vector, as well,
        and as a matrix, otherwise.

        .. SEEALSO::

            :meth:`solve_left`

        EXAMPLES::

            sage: A = matrix(QQ, 3, [1,2,3,-1,2,5,2,3,1])
            sage: b = vector(QQ, [1,2,3])
            sage: x = A.solve_right(b); x
            (-13/12, 23/12, -7/12)
            sage: A * x
            (1, 2, 3)

        We solve with A nonsquare::

            sage: A = matrix(QQ, 2,4, [0, -1, 1, 0, -2, 2, 1, 0])
            sage: B = matrix(QQ, 2,2, [1, 0, 1, -1])
            sage: X = A.solve_right(B); X
            [-3/2  1/2]
            [  -1    0]
            [   0    0]
            [   0    0]
            sage: A*X == B
            True

        Another nonsingular example::

            sage: A = matrix(QQ,2,3, [1,2,3,2,4,6]); v = vector([-1/2,-1])
            sage: x = A.solve_right(v); x
            (-1/2, 0, 0)
            sage: A*x == v
            True

        Same example but over `\ZZ`::

            sage: A = matrix(ZZ,2,3, [1,2,3,2,4,6]); v = vector([-1,-2])
            sage: A.solve_right(v)
            (-1, 0, 0)

        An example in which there is no solution::

            sage: A = matrix(QQ,2,3, [1,2,3,2,4,6]); v = vector([1,1])
            sage: A.solve_right(v)
            Traceback (most recent call last):
            ...
            ValueError: matrix equation has no solutions

        A :exc:`ValueError` is raised if the input is invalid::

            sage: A = matrix(QQ, 4,2, [0, -1, 1, 0, -2, 2, 1, 0])
            sage: B = matrix(QQ, 2,2, [1, 0, 1, -1])
            sage: X = A.solve_right(B)
            Traceback (most recent call last):
            ...
            ValueError: number of rows of self must equal number of rows of
            right-hand side

        We solve with A singular::

            sage: A = matrix(QQ, 2,3, [1,2,3,2,4,6]); B = matrix(QQ, 2,2, [6, -6, 12, -12])
            sage: X = A.solve_right(B); X
            [ 6 -6]
            [ 0  0]
            [ 0  0]
            sage: A*X == B
            True

        We illustrate left associativity, etc., of the ``solve_right`` operator.

        ::

            sage: A = matrix(QQ, 2, [1,2,3,4])
            sage: A.solve_right(A)
            [1 0]
            [0 1]
            sage: (A.solve_right(A)).solve_right(A)
            [1 2]
            [3 4]
            sage: A.parent()(1).solve_right(A)
            [1 2]
            [3 4]
            sage: A.solve_right(A.solve_right(A))
            [  -2    1]
            [ 3/2 -1/2]
            sage: X = A.solve_right(A - 2); X
            [ 5 -2]
            [-3  2]
            sage: A * X
            [-1  2]
            [ 3  2]

        Solving over a polynomial ring::

            sage: x = polygen(QQ, 'x')
            sage: A = matrix(2, [x, 2*x, -5*x^2 + 1, 3])
            sage: v = vector([3, 4*x - 2])
            sage: X = A.solve_right(v)
            sage: X
            ((-4/5*x^2 + 2/5*x + 9/10)/(x^3 + 1/10*x),
             (19/10*x^2 - 1/5*x - 3/10)/(x^3 + 1/10*x))
            sage: A * X == v
            True

        Solving some systems over `\ZZ/n\ZZ`::

            sage: # needs sage.libs.pari
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
            [  5]
            [108]
            [127]
            sage: B = B.column(0)
            sage: A.solve_right(B)
            (5, 108, 127)
            sage: A = Matrix(Zmod(15), 3,4, range(12))
            sage: B = Matrix(Zmod(15), 3,3, range(3,12))
            sage: X = A.solve_right(B)
            sage: A*X == B
            True

        Solving a system over the `p`-adics::

            sage: # needs sage.rings.padics
            sage: k = Qp(5, 4)
            sage: a = matrix(k, 3, [1,7,3, 2,5,4, 1,1,2]); a
            [    1 + O(5^4) 2 + 5 + O(5^4)     3 + O(5^4)]
            [    2 + O(5^4)     5 + O(5^5)     4 + O(5^4)]
            [    1 + O(5^4)     1 + O(5^4)     2 + O(5^4)]
            sage: v = vector(k, 3, [1,2,3])
            sage: x = a.solve_right(v); x
            (4 + 5 + 5^2 + 3*5^3 + O(5^4), 2 + 5 + 3*5^2 + 5^3 + O(5^4), 1 + 5 + O(5^4))
            sage: a * x == v
            True

        Solving a system of linear equations symbolically using symbolic
        matrices::

            sage: # needs sage.symbolic
            sage: var('a,b,c,d,x,y')
            (a, b, c, d, x, y)
            sage: A = matrix(SR, 2, [a,b,c,d]); A
            [a b]
            [c d]
            sage: result = vector(SR, [3,5]); result
            (3, 5)
            sage: soln = A.solve_right(result); soln
            (-b*(3*c/a - 5)/(a*(b*c/a - d)) + 3/a, (3*c/a - 5)/(b*c/a - d))
            sage: (a*x + b*y).subs(x=soln[0], y=soln[1]).simplify_full()
            3
            sage: (c*x + d*y).subs(x=soln[0], y=soln[1]).simplify_full()
            5
            sage: (A*soln).apply_map(lambda x: x.simplify_full())
            (3, 5)

        Over inexact rings, the output of this function may not be an exact
        solution. For example, over the real or complex double field,
        this computes a least-squares solution::

            sage: A = matrix(RDF, 3, 2, [1, 3, 4, 2, 0, -3])
            sage: b = vector(RDF, [5, 6, 1])
            sage: A.solve_right(b)  # tol 1e-14                                         # needs scipy
            (1.4782608695652177, 0.35177865612648235)
            sage: ~(A.T * A) * A.T * b  # closed form solution, tol 1e-14               # needs scipy
            (1.4782608695652177, 0.35177865612648235)

        Over the reals::

            sage: A = matrix(RDF, 3,3, [1,2,5,7.6,2.3,1,1,2,-1]); A
            [ 1.0  2.0  5.0]
            [ 7.6  2.3  1.0]
            [ 1.0  2.0 -1.0]
            sage: b = vector(RDF, [1,2,3])
            sage: x = A.solve_right(b); x  # tol 1e-14                                  # needs scipy
            (-0.1136950904392765, 1.3901808785529717, -0.33333333333333337)
            sage: x.parent()                                                            # needs scipy
            Vector space of dimension 3 over Real Double Field
            sage: A*x  # tol 1e-14                                                      # needs scipy
            (1.0, 1.9999999999999996, 3.0000000000000004)

        Over the complex numbers::

            sage: # needs sage.rings.complex_double sage.symbolic
            sage: A = matrix(CDF, [[      0, -1 + 2*I,  1 - 3*I,        I],
            ....:                  [2 + 4*I, -2 + 3*I, -1 + 2*I,   -1 - I],
            ....:                  [  2 + I,    1 - I,       -1,        5],
            ....:                  [    3*I,   -1 - I,   -1 + I,   -3 + I]])
            sage: b = vector(CDF, [2 - 3*I, 3, -2 + 3*I, 8])
            sage: x = A.solve_right(b); x                                               # needs scipy
            (1.96841637... - 1.07606761...*I, -0.614323843... + 1.68416370...*I,
             0.0733985765... + 1.73487544...*I, -1.6018683... + 0.524021352...*I)
            sage: x.parent()                                                            # needs scipy
            Vector space of dimension 4 over Complex Double Field
            sage: abs(A*x - b) < 1e-14                                                  # needs scipy
            True

        If ``b`` is given as a matrix, the result will be a matrix, as well::

            sage: A = matrix(RDF, 3, 3, [1, 2, 2, 3, 4, 5, 2, 2, 2])
            sage: b = matrix(RDF, 3, 2, [3, 2, 3, 2, 3, 2])
            sage: A.solve_right(b)  # tol 1e-14                                         # needs scipy
            [ 0.0  0.0]
            [ 4.5  3.0]
            [-3.0 -2.0]

        If `A` is a non-square matrix, the result is a least-squares solution.
        For a wide matrix, this may give a solution with a least-squares error
        of almost zero::

            sage: A = matrix(RDF, 2, 3, [1, 3, 4, 2, 0, -3])
            sage: b = vector(RDF, [5, 6])
            sage: x = A.solve_right(b)                                                  # needs scipy
            sage: (A * x - b).norm() < 1e-14                                            # needs scipy
            True

        For a tall matrix `A`, the error is usually not small::

            sage: A = matrix(RDF, 3, 2, [1, 3, 4, 2, 0, -3])
            sage: b = vector(RDF, [5, 6, 1])
            sage: x = A.solve_right(b)                                                  # needs scipy
            sage: (A * x - b).norm()  # tol 1e-14                                       # needs scipy
            3.2692119900020438

        TESTS:

        Check that the arguments are coerced to a suitable parent
        (:issue:`12406`)::

            sage: A = matrix(QQ, 2, [1, 2, 3, 4])
            sage: b = vector(RDF, [pi, e])                                              # needs sage.symbolic
            sage: A.solve_right(b)  # tol 1e-15                                         # needs scipy sage.symbolic
            (-3.564903478720541, 3.353248066155167)

            sage: R.<t> = ZZ[]
            sage: b = vector(R, [1, t])
            sage: x = A.solve_right(b); x
            (t - 2, -1/2*t + 3/2)
            sage: A * x == b
            True
            sage: x.base_ring()
            Fraction Field of Univariate Polynomial Ring in t over Rational Field

        ::

            sage: A = Matrix(Zmod(6), 3, 2, [1,2,3,4,5,6])
            sage: b = vector(ZZ, [1,1,1])
            sage: A.solve_right(b).base_ring() is Zmod(6)                               # needs sage.libs.pari
            True

        Check that the coercion mechanism gives consistent results
        (:issue:`12406`)::

            sage: A = matrix(ZZ, [[1, 2, 3], [2, 0, 2], [3, 2, 5]])
            sage: b = vector(RDF, [1, 1, 1])
            sage: A.solve_right(b) == A.change_ring(RDF).solve_right(b)
            ...
            True

        A degenerate case::

            sage: A = matrix(RDF, 0, 0, [])
            sage: A.solve_right(vector(RDF,[]))
            ()

        Over an inexact ring like ``RDF``, the coefficient matrix of a
        square system must be nonsingular::

            sage: A = matrix(RDF, 5, range(25))
            sage: b = vector(RDF, [1,2,3,4,5])
            sage: A.solve_right(b)
            Traceback (most recent call last):
            ...
            LinAlgError: Matrix is singular.

        The vector of constants needs the correct degree.  ::

            sage: A = matrix(RDF, 5, range(25))
            sage: b = vector(RDF, [1,2,3,4])
            sage: A.solve_right(b)
            Traceback (most recent call last):
            ...
            ValueError: number of rows of self must equal degree of
            right-hand side

        The vector of constants needs to be compatible with
        the base ring of the coefficient matrix.  ::

            sage: F.<a> = FiniteField(27)                                               # needs sage.rings.finite_rings
            sage: b = vector(F, [a,a,a,a,a])                                            # needs sage.rings.finite_rings
            sage: A.solve_right(b)                                                      # needs sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            TypeError: no common canonical parent for objects with parents: ...

        Check that coercions work correctly (:issue:`17405`)::

            sage: # needs sage.rings.complex_double sage.symbolic
            sage: A = matrix(RDF, 2, range(4))
            sage: b = vector(CDF, [1+I, 2])
            sage: A.solve_right(b)
            (-0.5 - 1.5*I, 1.0 + 1.0*I)
            sage: b = vector(QQ[I], [1+I, 2])
            sage: x = A.solve_right(b)

        This method must be called with a vector or a matrix::

            sage: A = matrix(CDF, 2, [1 for i in range(4)])
            sage: x = A.solve_right([1]*2)
            Traceback (most recent call last):
            ...
            TypeError: the second argument must be a vector or a matrix

        Over inexact rings, the ``check`` parameter is ignored as the result is
        only an approximate solution (:issue:`13932`)::

            sage: RF = RealField(52)
            sage: B = matrix(RF, 2, 2, 1)
            sage: A = matrix(RF, [[0.24, 1, 0], [1, 0, 0]])
            sage: 0 <= (A * A.solve_right(B) - B).norm() < 1e-14
            True

        Over the inexact ring ``SR``, we can still verify the solution
        if all of the elements involved were exact to begin with; if
        any are inexact, however, the ``check`` is still skipped
        (:issue:`29729` and :issue:`33159`)::

            sage: m = matrix(SR, [0])                                                   # needs sage.symbolic
            sage: b = vector(SR, [1])                                                   # needs sage.symbolic
            sage: m.solve_right(b, check=True)                                          # needs sage.symbolic
            Traceback (most recent call last):
            ...
            ValueError: matrix equation has no solutions

        In this case, turning off the ``check`` leads to a wrong result::

            sage: m.solve_right(b, check=False)                                         # needs sage.symbolic
            (0)

        In the following, we have an inexact entry in the matrix, so
        the ``check`` is still skipped leading to a wrong result::

            sage: m = matrix(SR, [0.0])                                                 # needs sage.symbolic
            sage: m.solve_right(b, check=True)                                          # needs sage.symbolic
            (0)

        ::

            sage: # needs sage.symbolic
            sage: SC = SR.subring(no_variables=True)
            sage: m = matrix(SC, [0])
            sage: b = vector(SC, [1])
            sage: m.solve_right(b)
            Traceback (most recent call last):
            ...
            ValueError: matrix equation has no solutions

        Check that :issue:`28586` is fixed::

            sage: m = matrix(GF(3), 2, 2, [1,2,2,0], sparse=True)
            sage: v = vector(GF(3), [1,1])
            sage: m.solve_right(v)
            (2, 1)
        """
        try:
            L = B.base_ring()
        except AttributeError:
            raise TypeError("the second argument must be a vector or a matrix")
        b_is_vec = isinstance(B, Vector)
        if b_is_vec:
            if self.nrows() != B.degree():
                raise ValueError("number of rows of self must equal "
                                 "degree of right-hand side")
        else:
            if self.nrows() != B.nrows():
                raise ValueError("number of rows of self must equal "
                                 "number of rows of right-hand side")

        K = self.base_ring()
        L = B.base_ring()
        # first coerce both elements to parent over same base ring
        P = K if L is K else coercion_model.common_parent(K, L)
        if P not in _Fields and P.is_integral_domain() and extend:
            # the non-integral-domain case is handled separately below
            P = P.fraction_field()
        if L is not P:
            B = B.change_ring(P)
        if K is not P:
            K = P
            self = self.change_ring(P)

        # If our field is inexact, checking the answer is doomed in
        # most cases. But here we handle the special ones.
        if isinstance(K, sage.rings.abc.SymbolicRing):
            # Elements of SR "remember" whether or not they are exact.
            # If every element in the system is exact, we can probably
            # still check the solution over the inexact ring SR.
            check = (check and all(e.is_exact()
                                   for e in self.list() + B.list()))
        else:
            check = (check and K.is_exact())

        if not K.is_integral_domain():
            # The non-integral-domain case is handled almost entirely
            # separately.
            if isinstance(K, sage.rings.abc.IntegerModRing):
                from sage.libs.pari import pari
                A = pari(self.lift())
                b = pari(B).lift()
                if b.type() == "t_MAT":
                    X = []
                    for n in range(B.ncols()):
                        ret = A.matsolvemod(K.cardinality(), b[n])
                        if ret.type() == 't_INT':
                            raise ValueError("matrix equation has no solutions")
                        X.append(ret.sage())
                    X = self.matrix_space(B.ncols(), self.ncols())(X)
                    return X.T
                elif b.type() == "t_VEC":
                    b = b.Col()
                    ret = A.matsolvemod(K.cardinality(), b)
                    if ret.type() == 't_INT':
                        raise ValueError("matrix equation has no solutions")
                    ret = ret.Vec().sage()
                    return (K ** self.ncols())(ret)
            raise TypeError("base ring must be an integral domain or a ring of integers mod n")

        C = B.column() if b_is_vec else B

        if not extend:
            try:
                X = self._solve_right_hermite_form(C)
            except NotImplementedError:
                X = self._solve_right_smith_form(C)
            return X.column(0) if b_is_vec else X

        if not self.is_square():
            X = self._solve_right_general(C, check=check)
        else:
            try:
                X = self._solve_right_nonsingular_square(C, check_rank=True)
            except NotFullRankError:
                X = self._solve_right_general(C, check=check)

        return X.column(0) if b_is_vec else X

    def _solve_right_nonsingular_square(self, B, check_rank=True):
        r"""
        If ``self`` is a matrix `A` of full rank, then this function
        returns a matrix `X` such that `A X = B`.

        .. SEEALSO::

           :meth:`solve_right` and :meth:`solve_left`

        INPUT:

        - ``B`` -- a matrix

        - ``check_rank`` -- boolean (default: ``True``)

        OUTPUT: matrix

        EXAMPLES::

            sage: A = matrix(QQ, 3, [1,2,4,5,3,1,1,2,-1])
            sage: B = matrix(QQ, 3,2, [1,5,1,2,1,5])
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
        # this could probably be optimized so that the rank computation is
        # avoided
        if check_rank and self.rank() < self.nrows():
            raise NotFullRankError
        D = self.augment(B)
        D.echelonize()
        return D.matrix_from_columns(range(self.ncols(), D.ncols()))

    def pivot_rows(self):
        """
        Return the pivot row positions for this matrix, which are a topmost
        subset of the rows that span the row space and are linearly
        independent.

        OUTPUT: a tuple of integers

        EXAMPLES::

            sage: A = matrix(QQ, 3,3, [0,0,0,1,2,3,2,4,6]); A
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
        to solve for ``self*X = B`` when ``self`` is not square or not of full
        rank. It does some linear algebra, then solves a full-rank square
        system.

        INPUT:

        - ``B`` -- a matrix

        - ``check`` -- boolean (default: ``True``); if ``False``, if there
          is no solution this function will not detect that fact

        OUTPUT: matrix

        EXAMPLES::

            sage: A = matrix(QQ, 2,3, [1,2,3,2,4,6]); B = matrix(QQ, 2,2, [6, -6, 12, -12])
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
                raise ValueError("matrix equation has no solutions")
        return X

    def _solve_right_smith_form(self, B):
        r"""
        Solve a matrix equation over a PID using the Smith normal form.

        EXAMPLES::

            sage: A = matrix(ZZ, 5, 7, [(1 + x^4) % 55 for x in range(5*7)]); A
            [ 1  2 17 27 37 21 32]
            [37 27 17 46 12  2 17]
            [27 26 32 32 37 27  6]
            [ 2 12  2 17 16 37 32]
            [32 37 16 17  2 12  2]

            sage: y = vector(ZZ, [1, 2, 3, 4, 5])
            sage: x, = A._solve_right_smith_form(y.column()).columns()
            sage: x.parent()
            Ambient free module of rank 7 over the principal ideal domain Integer Ring
            sage: A * x == y
            True
            sage: x
            (10, 1530831087980480, -2969971929450215, -178745029498097, 2320752168397186, -806846536262381, -520939892126393)

            sage: z = vector(ZZ, [-4, -1, 1, 5, 14, 31, 4])
            sage: x, = A.transpose()._solve_right_smith_form(z.column()).columns()
            sage: x.parent()
            Ambient free module of rank 5 over the principal ideal domain Integer Ring
            sage: x * A == z
            True
            sage: x
            (-1, 0, 1, 1, -1)

            sage: X = A._solve_right_smith_form(identity_matrix(ZZ,5))
            sage: X.parent()
            Full MatrixSpace of 7 by 5 dense matrices over Integer Ring
            sage: X
            [              -1                0                2                0                1]
            [-156182670342972   -2199494166584  310625144000132          1293916  151907461896112]
            [ 303010665531453    4267248023008 -602645168043247         -2510332 -294716315371323]
            [  18236418267658     256820398262  -36269645268572          -151082  -17737230430447]
            [-236774176922867   -3334450741532  470910201757143          1961587  230292926737068]
            [  82318322106118    1159275026338 -163719448527234          -681977  -80065012022313]
            [  53148766104440     748485096017 -105705345467375          -440318  -51693918051894]
        """
        S,U,V = self.smith_form()

        n,m = self.dimensions()
        r = B.ncols()

        X_ = []
        for d, v in zip(S.diagonal(), (U*B).rows()):
            if d:
                X_.append(v / d)
            elif v:
                raise ValueError("matrix equation has no solutions")
            else:
                X_.append([0] * r)

        X_ += [[0] * r] * (m - n)  # arbitrary

        from sage.matrix.constructor import matrix
        try:
            X_ = matrix(self.base_ring(), m, r, X_)
        except TypeError:
            raise ValueError("matrix equation has no solutions")

        return V * X_

    def _solve_right_hermite_form(self, B):
        r"""
        Solve a matrix equation over a PID using the Hermite normal form.

        EXAMPLES::

            sage: A = matrix(ZZ, 5, 7, [(1 + x^4) % 55 for x in range(5*7)]); A
            [ 1  2 17 27 37 21 32]
            [37 27 17 46 12  2 17]
            [27 26 32 32 37 27  6]
            [ 2 12  2 17 16 37 32]
            [32 37 16 17  2 12  2]

            sage: y = vector(ZZ, [1, 2, 3, 4, 5])
            sage: x, = A._solve_right_hermite_form(y.column()).columns()
            sage: x.parent()
            Ambient free module of rank 7 over the principal ideal domain Integer Ring
            sage: A * x == y
            True
            sage: x
            (7903368738919, 1664769092987, -13561109876830, -5130794714802, 11219929764616, -2728570619520, 0)

            sage: z = vector(ZZ, [-4, -1, 1, 5, 14, 31, 4])
            sage: x, = A.transpose()._solve_right_hermite_form(z.column()).columns()
            sage: x.parent()
            Ambient free module of rank 5 over the principal ideal domain Integer Ring
            sage: x * A == z
            True
            sage: x
            (-1, 0, 1, 1, -1)

            sage: X = A._solve_right_hermite_form(identity_matrix(ZZ,5))
            sage: X.parent()
            Full MatrixSpace of 7 by 5 dense matrices over Integer Ring
            sage: X
            [  682282983347  1029536563053   550188901703   877862915448       -1147487]
            [  143716389918   216862037808   115892034032   184913433473        -241707]
            [-1170705152437 -1766545243552  -944049606627 -1506293815517        1968932]
            [ -442931873812  -668365722378  -357177603912  -569900577297         744938]
            [  968595469303  1461570161933   781069571508  1246248350502       -1629017]
            [ -235552378240  -355438713600  -189948023680  -303074680960         396160]
            [             0              0              0              0              0]
        """
        H,U = self.transpose().hermite_form(transformation=True)
        H = H.transpose()
        U = U.transpose()
#        assert self*U == H

        n,m = self.dimensions()
        r = B.ncols()

        from sage.matrix.constructor import matrix
        X_ = matrix(self.base_ring(), m, r)
        for i in range(min(n,m)):
            v = B[i,:]
            v -= H[i,:i] * X_[:i]
            d = H[i][i]
            try:
                X_[i] = v / d
            except (ZeroDivisionError, TypeError) as e:
                raise ValueError("matrix equation has no solution")
#        assert H*X_ == B

        return U * X_

    def prod_of_row_sums(self, cols):
        r"""
        Calculate the product of all row sums of a submatrix of `A`
        for a list of selected columns ``cols``.

        EXAMPLES::

            sage: a = matrix(QQ, 2, 2, [1,2,3,2]); a
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
        """
        cdef Py_ssize_t c, row
        pr = 1
        for row from 0 <= row < self._nrows:
            tmp = []
            for c in cols:
                tmp.append(self.get_unsafe(row, c))
            pr = pr * sum(tmp)
        return pr

    def elementwise_product(self, right):
        r"""
        Return the elementwise product of two matrices
        of the same size (also known as the Hadamard product).

        INPUT:

        - ``right`` -- the right operand of the product.  A matrix
          of the same size as ``self`` such that multiplication
          of elements of the base rings of ``self`` and ``right``
          is defined, once Sage's coercion model is applied.  If
          the matrices have different sizes, or if multiplication
          of individual entries cannot be achieved, a :exc:`TypeError`
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

            sage: A = matrix(ZZ, 2, 3, range(6))
            sage: B = matrix(QQ, 2, 3, [5, 1/3, 2/7, 11/2, -3/2, 8])
            sage: C = A.elementwise_product(B)
            sage: C
            [   0  1/3  4/7]
            [33/2   -6   40]
            sage: C.parent()
            Full MatrixSpace of 2 by 3 dense matrices over Rational Field


        Notice the base ring of the results in the next two examples.  ::

            sage: x = polygen(ZZ, 'x')
            sage: D = matrix(ZZ['x'],2,[1+x^2,2,3,4-x])
            sage: E = matrix(QQ,2,[1,2,3,4])
            sage: F = D.elementwise_product(E)
            sage: F
            [  x^2 + 1         4]
            [        9 -4*x + 16]
            sage: F.parent()
            Full MatrixSpace of 2 by 2 dense matrices
             over Univariate Polynomial Ring in x over Rational Field

        ::

            sage: G = matrix(GF(3), 2, [0, 1, 2, 2])
            sage: H = matrix(ZZ, 2, [1, 2, 3, 4])
            sage: J = G.elementwise_product(H)
            sage: J
            [0 2]
            [0 2]
            sage: J.parent()
            Full MatrixSpace of 2 by 2 dense matrices
             over Finite Field of size 3

        Non-commutative rings behave as expected.  These are the usual quaternions. ::

            sage: # needs sage.combinat
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

            sage: A = random_matrix(ZZ, 5, 10, x=20)
            sage: A.elementwise_product(vector(ZZ, [1,2,3,4]))
            Traceback (most recent call last):
            ...
            TypeError: no common canonical parent for objects with parents:
            'Full MatrixSpace of 5 by 10 dense matrices over Integer Ring' and
            'Ambient free module of rank 4 over the principal ideal domain Integer Ring'

            sage: A = matrix(2, 2, range(4))
            sage: A.elementwise_product(polygen(parent(A)))
            Traceback (most recent call last):
            ...
            TypeError: elementwise_product() argument should be a matrix
            or coercible to a matrix

        Matrices of different sizes for operands will raise an error. ::

            sage: A = random_matrix(ZZ, 5, 10, x=20)
            sage: B = random_matrix(ZZ, 10, 5, x=40)
            sage: A.elementwise_product(B)
            Traceback (most recent call last):
            ...
            TypeError: no common canonical parent for objects with parents:
            'Full MatrixSpace of 5 by 10 dense matrices over Integer Ring' and
            'Full MatrixSpace of 10 by 5 dense matrices over Integer Ring'

        Some pairs of rings do not have a common parent where
        multiplication makes sense.  This will raise an error. ::

            sage: A = matrix(QQ, 3, 2, range(6))
            sage: B = matrix(GF(3), 3, [2]*6)
            sage: A.elementwise_product(B)
            Traceback (most recent call last):
            ...
            TypeError: no common canonical parent for objects with parents:
             'Full MatrixSpace of 3 by 2 dense matrices over Rational Field' and
             'Full MatrixSpace of 3 by 2 dense matrices over Finite Field of size 3'

        We illustrate various combinations of sparse and dense matrices.
        The usual coercion rules apply::

            sage: A = matrix(ZZ, 5, 6, range(30), sparse=False)
            sage: B = matrix(ZZ, 5, 6, range(30), sparse=True)
            sage: C = matrix(QQ, 5, 6, range(30), sparse=True)
            sage: A.elementwise_product(C).is_sparse()
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
        """
        # Optimized routines for specialized classes would likely be faster
        # See the "pairwise_product" of vectors for some guidance on doing this
        if have_same_parent(self, right):
            return self._elementwise_product(right)
        A, B = coercion_model.canonical_coercion(self, right)
        if not isinstance(A, Matrix):
            # Canonical coercion is not a matrix?!
            raise TypeError("elementwise_product() argument should be a matrix or coercible to a matrix")
        return (<Matrix>A)._elementwise_product(B)

    def permanent(self, algorithm='Ryser'):
        r"""
        Return the permanent of this matrix.

        Let `A = (a_{i,j})` be an `m \times n` matrix over any
        commutative ring with `m \le n`. The permanent of `A` is

        .. MATH::

           \mathrm{per}(A)
           = \sum_\pi a_{1,\pi(1)} a_{2,\pi(2)} \cdots a_{m,\pi(m)}

        where the summation extends over all one-to-one functions
        `\pi` from `\{1, \ldots, m\}` to `\{1, \ldots, n\}`.

        The product
        `a_{1,\pi(1)} a_{2,\pi(2)} \cdots a_{m,\pi(m)}` is
        called *diagonal product*. So the permanent of an
        `m \times n` matrix `A` is the sum of all the
        diagonal products of `A`.

        By default, this method uses Ryser's algorithm, but setting
        ``algorithm`` to "ButeraPernici" you can use the algorithm of Butera and
        Pernici (which is well suited for band matrices, i.e. matrices whose
        entries are concentrated near the diagonal).

        INPUT:

        - ``A`` -- matrix of size `m \times n` with `m \leq n`

        - ``algorithm`` -- either "Ryser" (default) or "ButeraPernici". The
          Butera-Pernici algorithm takes advantage of presence of zeros and is
          very well suited for sparse matrices.

        ALGORITHM:

        The Ryser algorithm is implemented in the method
        :meth:`_permanent_ryser`. It is a modification of theorem 7.1.1. from
        Brualdi and Ryser: Combinatorial Matrix Theory. Instead of deleting
        columns from `A`, we choose columns from `A` and calculate the product
        of the row sums of the selected submatrix.

        The Butera-Pernici algorithm is implemented in the function
        :func:`~sage.matrix.matrix_misc.permanental_minor_polynomial`. It takes
        advantage of cancellations that may occur in the computations.

        EXAMPLES::

            sage: A = ones_matrix(4,4)
            sage: A.permanent()
            24

            sage: A = matrix(3,6, [1,1,1,1,0,0,0,1,1,1,1,0,0,0,1,1,1,1])
            sage: A.permanent()
            36
            sage: B = A.change_ring(RR)
            sage: B.permanent()
            36.0000000000000

        The permanent above is directed to the Sloane's sequence :oeis:`A079908`
        ("The Dancing School Problems") for which the third term is 36:

        ::

            sage: oeis(79908)                           # optional -- internet
            A079908: Solution to the Dancing School Problem with 3 girls and n+3 boys: f(3,n).
            sage: _(3)                                  # optional -- internet
            36

        ::

            sage: A = matrix(4,5, [1,1,0,1,1,0,1,1,1,1,1,0,1,0,1,1,1,0,1,0])
            sage: A.permanent()
            32

        A huge permanent that cannot be reasonably computed with the Ryser
        algorithm (a `50 \times 50` band matrix with width `5`)::

            sage: n, w = 50, 5
            sage: A = matrix(ZZ, n, n, lambda i,j: (i+j)%5 + 1 if abs(i-j) <= w else 0)
            sage: A.permanent(algorithm='ButeraPernici')
            57766972735511097036962481710892268404670105604676932908

        See Minc: Permanents, Example 2.1, p. 5.

        ::

            sage: A = matrix(QQ, 2,2, [1/5,2/7,3/2,4/5])
            sage: A.permanent()
            103/175

        ::

            sage: R.<a> = PolynomialRing(ZZ)
            sage: A = matrix(R, 2,2, [a,1,a,a+1])
            sage: A.permanent()
            a^2 + 2*a

        ::

            sage: R.<x,y> = PolynomialRing(ZZ, 2)
            sage: A = matrix(R, 2,2, [x, y, x^2, y^2])
            sage: A.permanent()
            x^2*y + x*y^2
        """
        if algorithm == "Ryser":
            return self._permanent_ryser()

        elif algorithm == "ButeraPernici":
            return permanental_minor_polynomial(self, True)

        else:
            raise ValueError("algorithm must be one of \"Ryser\" or \"ButeraPernici\".")

    def _permanent_ryser(self):
        r"""
        Return the permanent computed using Ryser algorithm.

        See :meth:`permanent` for the documentation.

        EXAMPLES::

            sage: m = matrix([[1,1],[1,1]])
            sage: m._permanent_ryser()
            2
        """
        cdef Py_ssize_t m, n, r
        cdef int sn

        perm = 0
        m = self._nrows
        n = self._ncols
        if not m <= n:
            raise ValueError("must have m <= n, but m (=%s) and n (=%s)"%(m, n))

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

    def permanental_minor(self, Py_ssize_t k, algorithm='Ryser'):
        r"""
        Return the permanental `k`-minor of this matrix.

        The *permanental* `k`-*minor* of a matrix `A` is the sum of the
        permanents of all possible `k` by `k` submatrices of `A`. Note that the
        maximal permanental minor is just the permanent.

        For a (0,1)-matrix `A` the permanental `k`-minor
        counts the number of different selections of `k` 1s of
        `A` with no two of the 1s on the same row and no two of the
        1s on the same column.

        See Brualdi and Ryser: Combinatorial Matrix Theory, p. 203. Note
        the typo `p_0(A) = 0` in that reference! For applications
        see Theorem 7.2.1 and Theorem 7.2.4.

        .. SEEALSO::

            The method :meth:`rook_vector` returns the list of all permanental
            minors.

        INPUT:

        - ``k`` -- the size of the minor

        - ``algorithm`` -- either "Ryser" (default) or "ButeraPernici"; the
          Butera-Pernici algorithm is well suited for band matrices

        EXAMPLES::

            sage: A = matrix(4, [1,0,1,0,1,0,1,0,1,0,10,10,1,0,1,1])
            sage: A.permanental_minor(2)
            114

        ::

            sage: A = matrix(3,6, [1,1,1,1,0,0,0,1,1,1,1,0,0,0,1,1,1,1])
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

        The permanental minors of the "complement" matrix of `A` is
        related to the permanent of `A`::

            sage: m, n = 3, 6
            sage: C = matrix(m, n, lambda i,j: 1 - A[i,j])
            sage: sum((-1)^k * C.permanental_minor(k)*factorial(n-k)/factorial(n-m)
            ....:     for k in range(m+1))
            36

        See Theorem 7.2.1 of Brualdi and Ryser: Combinatorial Matrix
        Theory: per(A)

        TESTS::

            sage: A.permanental_minor(5)
            0
        """
        if algorithm == "Ryser":
            return self._permanental_minor_ryser(k)

        elif algorithm == "ButeraPernici":
            p = permanental_minor_polynomial(self, prec=k+1)
            return p[k]

        else:
            raise ValueError("algorithm must be one of \"Ryser\" or \"ButeraPernici\".")

    def _permanental_minor_ryser(self, Py_ssize_t k):
        r"""
        Compute the `k`-th permanental minor using Ryser algorithm.

        See :meth:`permanental_minor` for the documentation.

        EXAMPLES::

            sage: m = matrix([[1,2,1],[3,4,3],[5,6,5]])
            sage: m._permanental_minor_ryser(1)
            30
            sage: m._permanental_minor_ryser(2)
            174
            sage: m._permanental_minor_ryser(3)
            136
        """
        m = self._nrows
        n = self._ncols

        R = self._base_ring
        if k == 0:
            return R.one()
        if k > m:
            return R.zero()

        pm = 0
        for cols in _choose(n, k):
            for rows in _choose(m, k):
                pm = pm + self.matrix_from_rows_and_columns(rows, cols).permanent()
        return pm

    def pseudoinverse(self, *, algorithm=None):
        """
        Return the Moore-Penrose pseudoinverse of this matrix.

        INPUT:

        - ``algorithm`` -- (default: guess) one of the following:

          - ``'numpy'`` -- use numpy's ``linalg.pinv()`` which is
            suitable over real or complex fields

          - ``'exact'`` -- use a simple algorithm which is not
            numerically stable but useful over exact fields. Assume that
            no conjugation is needed, that the conjugate transpose is
            just the transpose.

          - ``'exactconj'`` -- like ``exact`` but use the conjugate
            transpose

        OUTPUT: a matrix

        EXAMPLES::

            sage: # needs sage.rings.complex_double sage.symbolic
            sage: M = diagonal_matrix(CDF, [0, I, 1+I]); M
            [        0.0         0.0         0.0]
            [        0.0       1.0*I         0.0]
            [        0.0         0.0 1.0 + 1.0*I]
            sage: M.pseudoinverse()  # tol 1e-15
            [        0.0         0.0         0.0]
            [        0.0      -1.0*I         0.0]
            [        0.0         0.0 0.5 - 0.5*I]

        We check the properties of the pseudoinverse over an exact
        field::

            sage: M = random_matrix(QQ, 6, 3) * random_matrix(QQ, 3, 5)
            sage: Mx = M.pseudoinverse()
            sage: M * Mx * M == M
            True
            sage: Mx * M * Mx == Mx
            True
            sage: (M * Mx).is_symmetric()
            True
            sage: (Mx * M).is_symmetric()
            True

        Beware that the ``exact`` algorithm is not numerically stable,
        but the default ``numpy`` algorithm is::

            sage: M = matrix.hilbert(12, ring=RR)
            sage: (~M * M).norm()  # a considerable error                               # needs scipy
            1.3...
            sage: Mx = M.pseudoinverse(algorithm='exact')
            sage: (Mx * M).norm()  # huge error                                         # needs scipy
            11.5...
            sage: Mx = M.pseudoinverse(algorithm='numpy')                               # needs numpy
            sage: (Mx * M).norm()  # still OK
            1.00...

        When multiplying the given matrix with the pseudoinverse, the
        result is symmetric for the ``exact`` algorithm or hermitian
        for the ``exactconj`` algorithm::

            sage: # needs sage.rings.number_field sage.symbolic
            sage: M = matrix(QQbar, 2, 2, [1, sqrt(-3), -sqrt(-3), 3])
            sage: M * M.pseudoinverse()
            [   0.2500000000000000?  0.4330127018922193?*I]
            [-0.4330127018922193?*I     0.750000000000000?]
            sage: M * M.pseudoinverse(algorithm='exactconj')
            [                   1/4  0.4330127018922193?*I]
            [-0.4330127018922193?*I                    3/4]
            sage: M * M.pseudoinverse(algorithm='exact')
            [                -1/2 0.866025403784439?*I]
            [0.866025403784439?*I                  3/2]

        For an invertible matrix, the pseudoinverse is just the
        inverse::

            sage: M = matrix([[1,2], [3,4]])
            sage: ~M
            [  -2    1]
            [ 3/2 -1/2]
            sage: M.pseudoinverse()
            [  -2    1]
            [ 3/2 -1/2]

        Numpy gives a strange answer due to rounding errors::

            sage: M.pseudoinverse(algorithm='numpy')  # random                          # needs numpy
            [-1286742750677287/643371375338643 1000799917193445/1000799917193444]
            [  519646110850445/346430740566963  -300239975158034/600479950316067]

        Although it is not too far off::

            sage: (~M - M.pseudoinverse(algorithm='numpy')).norm() < 1e-14              # needs numpy
            True

        TESTS::

            sage: M.pseudoinverse(algorithm='exact')
            [  -2    1]
            [ 3/2 -1/2]
            sage: M.pseudoinverse(algorithm='whatever')
            Traceback (most recent call last):
            ...
            ValueError: unknown algorithm 'whatever', valid values are ('numpy', 'exact', 'exactconj')
            sage: M.change_ring(RealField(54)).pseudoinverse()
            Traceback (most recent call last):
            ...
            NotImplementedError: pseudoinverse for real/complex field is only implemented for <= 53 bits of precision
        """
        ring = self.base_ring()
        if algorithm is None:
            # Choose algorithm depending on base ring
            from sage.rings.real_mpfr import RealField
            from sage.rings.complex_mpfr import ComplexField
            is_complex = ComplexField(2).has_coerce_map_from(ring)
            if is_complex:
                if ring.is_exact():
                    is_real = RealField(2).has_coerce_map_from(ring)
                    algorithm = "exact" if is_real else "exactconj"
                else:
                    if ring.precision() <= 53:
                        algorithm = "numpy"
                    else:
                        raise NotImplementedError("pseudoinverse for real/complex field is only implemented for <= 53 bits of precision")
            else:
                algorithm = "exact"
        else:
            algos = ("numpy", "exact", "exactconj")
            if algorithm not in algos:
                raise ValueError("unknown algorithm {!r}, valid values are {}".format(algorithm, algos))

        if algorithm == "numpy":
            from numpy.linalg import pinv
            from sage.matrix.constructor import matrix
            ans = pinv(self.numpy())
            return matrix(ring.fraction_field(), ans)

        # We use a simple algorithm taken from
        # https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_pseudoinverse#Rank_decomposition
        # Write self as A * B such that A and B both have full rank.
        B = self.row_space().basis_matrix()
        A = B.solve_left(self)

        if algorithm.endswith("conj"):
            At = A.conjugate_transpose()
            Bt = B.conjugate_transpose()
        else:
            At = A.transpose()
            Bt = B.transpose()

        # Now the pseudoinverse is B^t (A^t A B B^t)^-1 A^t.
        Q = (At * A) * (B * Bt)
        return Bt * ~Q * At

    def rook_vector(self, algorithm=None, complement=False, use_complement=None):
        r"""
        Return the rook vector of this matrix.

        Let `A` be an `m` by `n` (0,1)-matrix. We identify `A` with a chessboard
        where rooks can be placed on the fields `(i, j)` with `A_{i,j} = 1`. The
        number `r_k = p_k(A)` (the permanental `k`-minor) counts the number of
        ways to place `k` rooks on this board so that no rook can attack
        another.

        The *rook vector* of the matrix `A` is the list consisting of `r_0,
        r_1, \ldots, r_h`, where `h = min(m,n)`. The *rook polynomial* is defined by
        `r(x) = \sum_{k=0}^h r_k x^k`.

        The rook vector can be generalized to matrices defined over any rings
        using permanental minors. Among the available algorithms, only "Godsil"
        needs the condition on the entries to be either `0` or `1`.

        See :wikipedia:`Rook_polynomial` for more information and also the
        method :meth:`permanental_minor` to compute individual permanental
        minor.

        See also ``sage.matrix.matrix2.permanental_minor_polynomial``
        and the graph method ``matching_polynomial``.

        INPUT:

        - ``self`` -- an `m` by `n` matrix

        - ``algorithm`` -- string which must be either ``'Ryser'`` or
          ``'ButeraPernici'`` (default) or ``'Godsil'``; Ryser one might be
          faster on simple and small instances. Godsil only accepts input in 0,1.

        - ``complement`` -- boolean (default: ``False``); whether we consider the
          rook vector of the complement matrix. If set to ``True`` then the
          matrix must have entries in {0, 1} and the complement matrix is the
          one for which the 0s are replaced by 1s and 1s by 0s.

        - ``use_complement`` -- boolean (default: ``None``); whether to compute the
          rook vector of a (0,1)-matrix from its complement. By default this is
          determined by the density of ones in the matrix.

        EXAMPLES:

        The standard chessboard is an `8` by `8` grid in which any positions is
        allowed. In that case one gets that the number of ways to position `4`
        non-attacking rooks is `117600` while for `8` rooks it is `40320`::

            sage: ones_matrix(8,8).rook_vector()
            [1, 64, 1568, 18816, 117600, 376320, 564480, 322560, 40320]

        These numbers are the coefficients of a modified Laguerre polynomial::

            sage: x = polygen(QQ)
            sage: factorial(8) * laguerre(8,-x)                                         # needs sage.symbolic
            x^8 + 64*x^7 + 1568*x^6 + 18816*x^5 + 117600*x^4 + 376320*x^3 +
            564480*x^2 + 322560*x + 40320

        The number of derangements of length `n` is the permanent
        of a matrix with 0 on the diagonal and 1 elsewhere;
        for `n=21` it is `18795307255050944540` (see :oeis:`A000166`):

           sage: A = identity_matrix(21)
           sage: A.rook_vector(complement=True)[-1]
           18795307255050944540
           sage: Derangements(21).cardinality()                                         # needs sage.combinat
           18795307255050944540

        An other example that we convert into a rook polynomial::

            sage: A = matrix(3,6, [1,1,1,1,0,0,0,1,1,1,1,0,0,0,1,1,1,1])
            sage: A
            [1 1 1 1 0 0]
            [0 1 1 1 1 0]
            [0 0 1 1 1 1]
            sage: A.rook_vector()
            [1, 12, 40, 36]

            sage: R = PolynomialRing(ZZ, 'x')
            sage: R(A.rook_vector())
            36*x^3 + 40*x^2 + 12*x + 1

        Different algorithms are available::

            sage: A = matrix([[1,0,0,1],[0,1,1,0],[0,1,1,0],[1,0,0,1]])
            sage: A.rook_vector(algorithm='ButeraPernici')
            [1, 8, 20, 16, 4]
            sage: A.rook_vector(algorithm='Ryser')
            [1, 8, 20, 16, 4]
            sage: A.rook_vector(algorithm='Godsil')                                     # needs sage.graphs sage.libs.flint
            [1, 8, 20, 16, 4]

        When the matrix `A` has more ones then zeroes it is usually faster
        to compute the rook polynomial of the complementary matrix, with
        zeroes and ones interchanged, and use the inclusion-exclusion theorem,
        giving for a `m \times n` matrix `A` with complementary matrix `B`

        .. MATH::

            r_k(A) = \sum_{j=0}^k (-1)^j \binom{m-j}{k-j} \binom{n-j}{k-j} (k-j)! r_j(B)

        see [Rio1958]_ or the introductory text [AS2011]_. This can be done
        setting the argument ``use_complement`` to ``True``.

        An example with an exotic matrix (for which only Butera-Pernici and
        Ryser algorithms are available)::

            sage: R.<x,y> = PolynomialRing(GF(5))
            sage: A = matrix(R, [[1, x, y], [x*y, x**2+y, 0]])
            sage: A.rook_vector(algorithm='ButeraPernici')
            [1, x^2 + x*y + x + 2*y + 1, 2*x^2*y + x*y^2 + x^2 + y^2 + y]
            sage: A.rook_vector(algorithm='Ryser')
            [1, x^2 + x*y + x + 2*y + 1, 2*x^2*y + x*y^2 + x^2 + y^2 + y]
            sage: A.rook_vector(algorithm='Godsil')
            Traceback (most recent call last):
            ...
            ValueError: coefficients must be zero or one, but we have 'x' in position (0,1).
            sage: B = A.transpose()
            sage: B.rook_vector(algorithm='ButeraPernici')
            [1, x^2 + x*y + x + 2*y + 1, 2*x^2*y + x*y^2 + x^2 + y^2 + y]
            sage: B.rook_vector(algorithm='Ryser')
            [1, x^2 + x*y + x + 2*y + 1, 2*x^2*y + x*y^2 + x^2 + y^2 + y]

        TESTS::

            sage: matrix([[0,0],[0,0]]).rook_vector(algorithm='ButeraPernici')
            [1, 0, 0]
            sage: matrix([[0,0],[0,0]]).rook_vector(algorithm='Ryser')
            [1, 0, 0]
            sage: matrix([[0,0],[0,0]]).rook_vector(algorithm='Godsil')                 # needs sage.graphs sage.libs.flint
            [1, 0, 0]
            sage: matrix.ones(4, 2).rook_vector("Ryser")
            [1, 8, 12]
            sage: matrix.ones(4, 2).rook_vector("Godsil")                               # needs sage.graphs sage.libs.flint
            [1, 8, 12]
            sage: m = matrix(ZZ,4,5)
            sage: m[:4,:4] = identity_matrix(4)
            sage: algos = ["Ryser", "ButeraPernici"]
            sage: algos += ["Godsil"]
            sage: for algorithm in algos:                                               # needs sage.graphs sage.libs.flint
            ....:     v = m.rook_vector(complement=True, use_complement=True, algorithm=algorithm)
            ....:     if v != [1, 16, 78, 128, 53]:
            ....:         print("ERROR with algorithm={} use_complement=True".format(algorithm))
            ....:     v = m.rook_vector(complement=True, use_complement=False, algorithm=algorithm)
            ....:     v = m.rook_vector(complement=True, use_complement=False)
            ....:     if v != [1, 16, 78, 128, 53]:
            ....:         print("ERROR with algorithm={} use_complement=False".format(algorithm))
        """
        cdef Py_ssize_t i, j
        cdef unsigned int num_ones
        cdef int m = self._nrows
        cdef int n = self._ncols
        cdef int mn = min(m, n)
        cdef Matrix B
        zero = self.base_ring().zero()
        one = self.base_ring().one()

        if algorithm is None:
            algorithm = "ButeraPernici"

        # we first run through the coefficients of the matrix to compute the
        # number of nonzero coefficients and to see whether or not it contains
        # only elements in {0,1}... but this is not always needed
        if complement or use_complement is None or algorithm == "Godsil":
            # z2 flag is True if all coefficients belong to {0,1}
            z2 = True
            num_ones = 1
            for i in range(m):
                for j in range(n):
                    x = self.get_unsafe(i, j)
                    if x != zero:
                        if x != one:
                            z2 = False
                            break
                        else:
                            num_ones += 1
                else:
                    continue
                break

            if not z2 and (complement or algorithm == "Godsil"):
                raise ValueError("coefficients must be zero or one, but we have '{}' in position ({},{}).".format(x, i, j))

            if use_complement is None:
                use_complement = z2 and num_ones > 0.55 * m * n

        if use_complement:
            B = self.new_matrix()
            for i in range(m):
                for j in range(n):
                    B.set_unsafe(i, j, one-self.get_unsafe(i, j))
            b = B.rook_vector(algorithm=algorithm, use_complement=False)
            complement = not complement

        elif algorithm == "Ryser":
            b = [self.permanental_minor(k, algorithm='Ryser')
                 for k in range(mn + 1)]

        elif algorithm == "ButeraPernici":
            p = permanental_minor_polynomial(self)
            b = [p[k] for k in range(mn + 1)]

        elif algorithm == "Godsil":
            from sage.graphs.bipartite_graph import BipartiteGraph
            p = BipartiteGraph(self).matching_polynomial()
            d = p.degree()
            b = [p[i] * (-1)**((d - i)//2) for i in range(d, d-2*mn-1, -2)]

        else:
            raise ValueError('algorithm must be one of "Ryser", "ButeraPernici" or "Godsil".')

        # now compute the permanental minor of the complement matrix if needed
        if complement:
            a = [one]
            c1 = 1
            for k in range(1, mn + 1):
                c1 = (c1 * (m-k+1) * (n-k+1)) // k
                c = c1
                s = c*b[0] + (-one)**k*b[k]
                for j in range(1, k):
                    c = -c * (k-j+1) // ((m-j+1) * (n-j+1))
                    s += c * b[j]
                a.append(s)
            return a
        else:
            return b

    def minors(self, k):
        r"""
        Return the list of all `k \times k` minors of ``self``.

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

            sage: A = Matrix(ZZ, 2,3, [1,2,3,4,5,6]); A
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
            sage: A = Matrix(P, 2, 3, [x0*x1, x0, x1, x2, x2 + 16, x2 + 5*x1])
            sage: A.minors(2)                                                           # needs sage.rings.finite_rings
            [x0*x1*x2 + 16*x0*x1 - x0*x2,
             5*x0*x1^2 + x0*x1*x2 - x1*x2,
             5*x0*x1 + x0*x2 - x1*x2 - 16*x1]

        This test addresses an issue raised at :issue:`20512`::

            sage: A.minors(0)[0].parent() == P
            True
        """
        from sage.combinat.combination import Combinations
        if k == 0:
            return [self.base_ring().one()]
        all_rows = range(self.nrows())
        all_cols = range(self.ncols())
        m = []
        for rows in Combinations(all_rows, k):
            for cols in Combinations(all_cols, k):
                m.append(self.matrix_from_rows_and_columns(rows, cols).determinant())
        return m

    def det(self, *args, **kwds):
        """
        Synonym for self.determinant(...).

        EXAMPLES::

            sage: A = MatrixSpace(Integers(8), 3)([1,7,3, 1,1,1, 3,4,5])
            sage: A.det()
            6
        """
        return self.determinant(*args, **kwds)

    def determinant(self, algorithm=None):
        r"""
        Return the determinant of ``self``.

        ALGORITHM:

        If the base ring has a method :meth:`_matrix_determinant`, we call it.

        Otherwise, for small matrices (n less than 4), this is computed using the
        naive formula. In the specific case of matrices over the integers modulo a
        non-prime, the determinant of a lift is computed over the integers.
        In general, the characteristic polynomial is computed either using
        the Hessenberg form (specified by ``'hessenberg'``) or the generic
        division-free algorithm (specified by ``'df'``).  When the base ring
        is an exact field, the default choice is ``'hessenberg'``, otherwise
        it is ``'df'``.  Note that for matrices over most rings, more
        sophisticated algorithms can be used. (Type ``A.determinant?`` to
        see what is done for a specific matrix ``A``.)

        INPUT:

        - ``algorithm`` -- string; one of
            - ``'df'`` -- generic O(n^4) division-free algorithm
            - ``'hessenberg'`` -- use the Hessenberg form of the matrix

        EXAMPLES::

            sage: A = MatrixSpace(Integers(8), 3)([1,7,3, 1,1,1, 3,4,5])
            sage: A.determinant()
            6
            sage: A.determinant() is A.determinant()
            True
            sage: A[0,0] = 10
            sage: A.determinant()
            7

        We compute the determinant of the arbitrary 3x3 matrix::

            sage: R = PolynomialRing(QQ, 9, 'x')
            sage: A = matrix(R, 3, R.gens())
            sage: A
            [x0 x1 x2]
            [x3 x4 x5]
            [x6 x7 x8]
            sage: A.determinant()
            -x2*x4*x6 + x1*x5*x6 + x2*x3*x7 - x0*x5*x7 - x1*x3*x8 + x0*x4*x8

        We create a matrix over `\ZZ[x,y]` and compute its
        determinant.

        ::

            sage: R.<x,y> = PolynomialRing(IntegerRing(), 2)
            sage: A = MatrixSpace(R,2)([x, y, x**2, y**2])
            sage: A.determinant()
            -x^2*y + x*y^2

        A matrix over a non-domain::

            sage: m = matrix(Integers(4), 2, [1,2,2,3])
            sage: m.determinant()
            3

        TESTS::

            sage: A = matrix(5, 5, [next_prime(i^2) for i in range(25)])                # needs sage.libs.pari
            sage: B = MatrixSpace(ZZ['x'], 5, 5)(A)                                     # needs sage.libs.pari
            sage: A.det() - B.det()                                                     # needs sage.libs.pari
            0

        We verify that :issue:`5569` is resolved (otherwise the following
        would hang for hours)::

            sage: d = random_matrix(GF(next_prime(10^20)), 50).det()                    # needs sage.rings.finite_rings
            sage: d = random_matrix(Integers(10^50), 50).det()                          # needs sage.rings.finite_rings

        We verify that :issue:`7704` is resolved::

            sage: matrix(ZZ, {(0,0): 1, (1,1): 2, (2,2): 3, (3,3): 4}).det()
            24
            sage: matrix(QQ, {(0,0): 1, (1,1): 2, (2,2): 3, (3,3): 4}).det()
            24

        We verify that :issue:`10063` is resolved::

            sage: # needs sage.libs.singular
            sage: A = GF(2)['x,y,z']
            sage: A.inject_variables()
            Defining x, y, z
            sage: R = A.quotient(x^2 + 1).quotient(y^2 + 1).quotient(z^2 + 1)
            sage: R.inject_variables()
            Defining xbarbarbar, ybarbarbar, zbarbarbar
            sage: M = matrix([[1, 1, 1, 1],
            ....:             [xbarbarbar, ybarbarbar, 1, 1],
            ....:             [0, 1, zbarbarbar, 1],
            ....:             [xbarbarbar, zbarbarbar, 1, 1]])
            sage: M.determinant()
            xbarbarbar*ybarbarbar*zbarbarbar + xbarbarbar*ybarbarbar
            + xbarbarbar*zbarbarbar + ybarbarbar*zbarbarbar + xbarbarbar
            + ybarbarbar + zbarbarbar + 1

        Check that the determinant is computed from a cached charpoly
        properly::

            sage: A = matrix(RR, [ [1, 0, 1/2],
            ....:                  [0, 1, 0  ],
            ....:                  [0, 0, -2 ] ])
            sage: B = copy(A)
            sage: _ = A.charpoly()
            sage: A.determinant() == B.determinant()
            True
        """
        cdef Py_ssize_t n
        n = self._ncols

        if self._nrows != n:
            raise ValueError("self must be a square matrix")

        d = self.fetch('det')
        if d is not None:
            return d

        # If charpoly known, then det is easy.
        f = self.fetch('charpoly')
        if f is not None:
            c = f[0]
            if n % 2:
                c = -c
            d = self._coerce_element(c)
            self.cache('det', d)
            return d

        # for base rings with their own specific determinant methods
        R = self._base_ring
        if hasattr(R, '_matrix_determinant'):
            d = R._matrix_determinant(self)
            self.cache('det', d)
            return d

        # For small matrices, you cannot beat the naive formula.
        if n <= 3:
            if n == 0:
                d = R.one()
            elif n == 1:
                d = self.get_unsafe(0, 0)
            elif n == 2:
                d = self.get_unsafe(0, 0)*self.get_unsafe(1, 1) - self.get_unsafe(1, 0)*self.get_unsafe(0, 1)
            elif n == 3:
                d = self.get_unsafe(0, 0) * (self.get_unsafe(1, 1)*self.get_unsafe(2, 2) - self.get_unsafe(1, 2)*self.get_unsafe(2, 1))    \
                    - self.get_unsafe(1, 0) * (self.get_unsafe(0, 1)*self.get_unsafe(2, 2) - self.get_unsafe(0, 2)*self.get_unsafe(2, 1))  \
                    + self.get_unsafe(2, 0) * (self.get_unsafe(0, 1)*self.get_unsafe(1, 2) - self.get_unsafe(0, 2)*self.get_unsafe(1, 1))
            self.cache('det', d)
            return d

        # Special case for Z/nZ or GF(p):
        if isinstance(R, sage.rings.abc.IntegerModRing) and self.is_dense():
            import sys
            # If the characteristic is prime and smaller than a machine
            # word, use PARI.
            ch = R.characteristic()
            if ch.is_prime() and ch < (2*sys.maxsize):
                try:
                    d = R(self.__pari__().matdet())
                except ImportError:
                    pass
                else:
                    self.cache('det', d)
                    return d
            # Lift to ZZ and compute there.
            d = R(self.apply_map(lambda x: x.lift_centered()).det())
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
        if (algorithm is None and R in _Fields and R.is_exact()) or (algorithm == "hessenberg"):
            try:
                charp = self.charpoly('x', algorithm='hessenberg')
            except ValueError:
                # Hessenberg algorithm not supported, so we use whatever the default algorithm is.
                charp = self.charpoly('x')
            c = charp[0]
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

        var = 'A0123456789' if isinstance(R, sage.rings.abc.SymbolicRing) else 'x'
        try:
            charp = self.charpoly(var, algorithm='df')
        except ValueError:
            # Division free algorithm not supported, so we use whatever the default algorithm is.
            charp = self.charpoly(var)
        c = charp[0]
        if self._nrows % 2:
            c = -c
        d = self._coerce_element(c)
        self.cache('det', d)
        return d

    cdef _det_by_minors(self, Py_ssize_t level):
        """
        Compute the determinant of the upper-left level x level submatrix
        of ``self``. Does not handle degenerate cases, level MUST be >= 2
        """
        cdef Py_ssize_t i
        if level == 2:
            return self.get_unsafe(0, 0) * self.get_unsafe(1, 1) - self.get_unsafe(0, 1) * self.get_unsafe(1, 0)
        else:
            level -= 1
            d = self.get_unsafe(level, level) * self._det_by_minors(level)
            # on each iteration, row i will be missing in the first (level) rows
            # swapping is much faster than taking submatrices
            for i from level > i >= 0:
                self.swap_rows(level, i)
                if (level - i) % 2:
                    d -= self.get_unsafe(level, level) * self._det_by_minors(level)
                else:
                    d += self.get_unsafe(level, level) * self._det_by_minors(level)
            # undo all our permutations to get us back to where we started
            for i from 0 <= i < level:
                self.swap_rows(level, i)
            return d

    def quantum_determinant(self, q=None):
        r"""
        Return the quantum deteminant of ``self``.

        The quantum determinant of a matrix `M = (m_{ij})_{i,j=1}^n`
        is defined by

        .. MATH::

            \det_q(M) =
            \sum_{\sigma \in S_n} (-q)^{\ell(\sigma)} M_{\sigma(i),j},

        where `S_n` is the symmetric group on `\{1, \ldots, n\}` and
        `\ell(\sigma)` denotes the length of `\sigma` written as simple
        transpositions (equivalently the number of inversions when
        written in one-line notation).

        INPUT:

        - ``q`` -- the parameter `q`; the default is `q \in F[q]`,
          where `F` is the base ring of ``self``

        EXAMPLES::

            sage: A = matrix(SR, 2, lambda i, j: f'a{i}{j}'); A                         # needs sage.symbolic
            [a00 a01]
            [a10 a11]
            sage: A.quantum_determinant()                                               # needs sage.symbolic
            -a01*a10*q + a00*a11

            sage: A = matrix(SR, 3, lambda i, j: f'a{i}{j}')                            # needs sage.symbolic
            sage: A.quantum_determinant()                                               # needs sage.symbolic
            -a02*a11*a20*q^3 + (a01*a12*a20 + a02*a10*a21)*q^2
             + (-a00*a12*a21 - a01*a10*a22)*q + a00*a11*a22

            sage: R.<q> = LaurentPolynomialRing(ZZ)
            sage: MS = MatrixSpace(Integers(8), 3)
            sage: A = MS([1,7,3, 1,1,1, 3,4,5])
            sage: A.det()
            6
            sage: A.quantum_determinant(q^-2)
            7*q^-6 + q^-4 + q^-2 + 5

            sage: S.<x,y> = PolynomialRing(GF(7))
            sage: R.<q> = LaurentPolynomialRing(S)
            sage: MS = MatrixSpace(S, 3, sparse=True)
            sage: A = MS([[x, y, 3], [4, 2+y, x^2], [0, 1-x, x+y]])
            sage: A.det()
            x^4 - x^3 + x^2*y + x*y^2 + 2*x^2 - 2*x*y + 3*y^2 + 2*x - 2
            sage: A.quantum_determinant()
            (2*x - 2)*q^2 + (x^4 - x^3 + 3*x*y + 3*y^2)*q + x^2*y + x*y^2 + 2*x^2 + 2*x*y
            sage: A.quantum_determinant(int(2))
            2*x^4 - 2*x^3 + x^2*y + x*y^2 + 2*x^2 + x*y - y^2 + x - 1
            sage: A.quantum_determinant(q*x + q^-1*y)
            (2*x*y^2 - 2*y^2)*q^-2 + (x^4*y - x^3*y + 3*x*y^2 + 3*y^3)*q^-1
             + (-2*x^2*y + x*y^2 + 2*x^2 - 2*x*y)
             + (x^5 - x^4 + 3*x^2*y + 3*x*y^2)*q + (2*x^3 - 2*x^2)*q^2
        """
        cdef Py_ssize_t n = self._ncols

        if self._nrows != n:
            raise ValueError("self must be a square matrix")

        if q is None:
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            q = PolynomialRing(self.base_ring(), 'q').gen()

        from sage.misc.misc_c import prod
        from sage.combinat.permutation import Permutations
        cdef Py_ssize_t i
        return sum((-q)**s.number_of_inversions()
                   * prod(self.get_unsafe(s[i] - 1, i) for i in range(n))
                   for s in Permutations(n))

    qdet = quantum_determinant

    def pfaffian(self, algorithm=None, check=True):
        r"""
        Return the Pfaffian of ``self``, assuming that ``self`` is an
        alternating matrix.

        The result is cached.

        INPUT:

        - ``algorithm`` -- string (default: ``None``); the algorithm to use.
          Currently the following algorithms have been implemented:

          * ``'bfl'`` -- using the Bär-Faddeev-LeVerrier algorithm
          * ``'definition'`` -- using the definition given by perfect
            matchings

        - ``check`` -- boolean (default: ``True``); whether to
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

        See [Knu1995]_, [DW1995]_ and [Rot2001]_, [Baer2020]_, just to name a
        few sources, for further properties of Pfaffians.

        ALGORITHM:

        If the matrix is small, namely up to size `4 \times 4`, the naive
        formulas are always used.

        The Bär-Faddeev-LeVerrier algorithm can be accessed using ``'bfl'``.
        It works over any `\QQ`-algebra or ring whose fraction field is an
        `\QQ`-algebra (see [Baer2020]_ for details). If that check fails,
        the implementation raises an error because correct results cannot be
        guaranteed.

        To access the algorithm using the above definition,
        use ``'definition'``.
        However, notice that this algorithm is usually very slow.

        By default, i.e. if no options are set, the implementation tries to
        apply the BFL algorithm first. If BFL is not applicable, it uses the
        definition by perfect matchings.

        The alternatingness of the matrix ``self`` is checked only if ``check``
        is ``True`` (this is important because even if ``self`` is alternating,
        a non-discrete base ring might prevent Sage from being able to check
        this).

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

        Let us compute the Pfaffian of a generic `4 \times 4` alternating
        matrix::

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

        In order to use the Bär-Faddeev-LeVerrier algorithm, the base ring
        must have characteristic zero::

            sage: A = matrix(GF(5), [(0, 3, 4, 1, 3, 4),
            ....:                    (2, 0, 2, 0, 1, 0),
            ....:                    (1, 3, 0, 4, 1, 0),
            ....:                    (4, 0, 1, 0, 2, 0),
            ....:                    (2, 4, 4, 3, 0, 0),
            ....:                    (1, 0, 0, 0, 0, 0)])
            sage: A.pfaffian(algorithm='bfl')
            Traceback (most recent call last):
            ...
            TypeError: Bär-Faddeev-LeVerrier algorithm not applicable,
             use another algorithm instead

        In that case, the definition by perfect matchings is used instead::

            sage: A.pfaffian()                                                          # needs sage.combinat
            2
        """

        pf = self.fetch('pfaffian')  # check out cache
        if pf is not None:
            return pf

        k = self._nrows

        if check:
            if k != self._ncols:
                raise ValueError("self must be a square matrix")
            if not self.is_alternating():
                raise ValueError("self must be alternating, which includes the diagonal entries being 0")

        # trivial cases:
        R = self.base_ring()
        if k % 2 == 1:
            pf = R.zero()
            # cache the result, and return it:
            self.cache('pfaffian', pf)
            return pf
        # For small matrices, you can't beat the naive formula:
        elif k <= 4:
            if k == 0:
                pf = R.one()
            elif k == 2:
                pf = self.get_unsafe(0, 1)
            elif k == 4:
                pf = self.get_unsafe(0, 1) * self.get_unsafe(2, 3) \
                     - self.get_unsafe(0, 2) * self.get_unsafe(1, 3) \
                     + self.get_unsafe(1, 2) * self.get_unsafe(0, 3)
            # cache the result, and return it:
            self.cache('pfaffian', pf)
            return pf

        # choose algorithm:
        if algorithm is None:
            if R in IntegralDomains():
                F = R.fraction_field()
            else:
                F = R
            if QQ.is_subring(F):
                temp = <Matrix> self.change_ring(F)
                pf = self._coerce_element(temp._pf_bfl())
            else:
                pf = self._pf_perfect_matchings()
        else:
            if algorithm == 'definition':
                pf = self._pf_perfect_matchings()
            elif algorithm == 'bfl':
                if R in IntegralDomains():
                    F = R.fraction_field()
                else:
                    F = R
                if not QQ.is_subring(F):
                    raise TypeError('Bär-Faddeev-LeVerrier algorithm not '
                                    'applicable, use another algorithm instead')
                temp = <Matrix> self.change_ring(F)
                pf = self._coerce_element(temp._pf_bfl())
            else:
                raise NotImplementedError("algorithm '%s' not recognized" % algorithm)
        # cache the result, and return it:
        self.cache('pfaffian', pf)
        return pf

    def _pf_perfect_matchings(self):
        r"""
        Compute the Pfaffian of ``self`` using the definition given by perfect
        matchings.

        OUTPUT: an element of the base ring of ``self`` representing the Pfaffian

        EXAMPLES::

            sage: A = matrix([(0, 2, -1/2, -2, 2, -1/2),
            ....:             (-2, 0, -1, 1, -1, 3/2),
            ....:             (1/2, 1, 0, 0, 3/2, 1),
            ....:             (2, -1, 0, 0, 1, 5/2),
            ....:             (-2, 1, -3/2, -1, 0, 1/2),
            ....:             (1/2, -3/2, -1, -5/2, -1/2, 0)])
            sage: A._pf_perfect_matchings()                                             # needs sage.combinat
            -1/2
        """
        R = self._base_ring
        k = self._nrows

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

    cdef _pf_bfl(self):
        r"""
        Compute the Pfaffian of ``self`` using the Baer-Faddeev-LeVerrier
        algorithm.

        .. WARNING::

            This method assumes that the base ring is an `\QQ`-algebra.

        OUTPUT:

        - an element (possibly coerced) originated from the base ring of
          ``self`` representing the Pfaffian

        EXAMPLES:

        Pfaffian of some matrix over the rationals using the
        Bär-Faddeev-LeVerrier algorithm::

            sage: A = matrix([(0, 2, -1/2, -2, 2, -1/2),
            ....:             (-2, 0, -1, 1, -1, 3/2),
            ....:             (1/2, 1, 0, 0, 3/2, 1),
            ....:             (2, -1, 0, 0, 1, 5/2),
            ....:             (-2, 1, -3/2, -1, 0, 1/2),
            ....:             (1/2, -3/2, -1, -5/2, -1/2, 0)])
            sage: A.pfaffian(algorithm='bfl')
            -1/2

        TESTS::

            sage: A = random_matrix(ZZ['x'], 6)
            sage: A = A - A.transpose()
            sage: A.pfaffian(algorithm='bfl') == A._pf_perfect_matchings()              # needs sage.combinat
            True
        """
        cdef Py_ssize_t n = self._ncols
        cdef Py_ssize_t q = n // 2
        cdef Py_ssize_t i, k

        # apply J:
        cdef Matrix A = <Matrix> copy(self)
        for i in range(0, n, 2):
            A.swap_columns_c(i, i+1)  # avoid checks
            for k in range(n):
                A.set_unsafe(k, i+1, -A.get_unsafe(k, i+1))

        cdef Matrix M = <Matrix> copy(A)

        # Baer-Faddeev-Leverrier algorithm:
        for k in range(1, q):
            c = -M.trace() / (2*k)
            # add c along the diagonal
            for i in range(n):
                M.set_unsafe(i, i, M.get_unsafe(i, i) + c)
            M = A * M
        c = -M.trace() / (2*q)
        return (-1)**q * c

    def apply_morphism(self, phi):
        """
        Apply the morphism phi to the coefficients of this dense matrix.

        The resulting matrix is over the codomain of phi.

        INPUT:

        - ``phi`` -- a morphism, so ``phi`` is callable and ``phi.domain()``
          and ``phi.codomain()`` are defined. The codomain must be a ring

        OUTPUT: a matrix over the codomain of phi

        EXAMPLES::

            sage: m = matrix(ZZ, 3, 3, range(9))
            sage: phi = ZZ.hom(GF(5))
            sage: m.apply_morphism(phi)
            [0 1 2]
            [3 4 0]
            [1 2 3]
            sage: parent(m.apply_morphism(phi))
            Full MatrixSpace of 3 by 3 dense matrices
             over Finite Field of size 5

        We apply a morphism to a matrix over a polynomial ring::

            sage: R.<x,y> = QQ[]
            sage: m = matrix(2, [x,x^2 + y, 2/3*y^2-x, x]); m
            [          x     x^2 + y]
            [2/3*y^2 - x           x]
            sage: phi = R.hom([y,x])
            sage: m.apply_morphism(phi)
            [          y     y^2 + x]
            [2/3*x^2 - y           y]
        """
        M = self.parent().change_ring(phi.codomain())
        if self.is_sparse():
            values = {(i, j): phi(z) for (i, j), z in self.dict()}
        else:
            values = [phi(z) for z in self.list()]
        image = M(values)
        if self._subdivisions is not None:
            image.subdivide(*self.subdivisions())
        return image

    def apply_map(self, phi, R=None, sparse=None):
        """
        Apply the given map phi (an arbitrary Python function or callable
        object) to this dense matrix. If R is not given, automatically
        determine the base ring of the resulting matrix.

        INPUT:

        - ``sparse`` -- boolean (default: ``False); ``True`` to make the output
          a sparse matrix

        - ``phi`` -- arbitrary Python function or callable object

        - ``R`` -- (optional) ring

        OUTPUT: a matrix over R

        EXAMPLES::

            sage: m = matrix(ZZ, 3, 3, range(9))
            sage: k.<a> = GF(9)                                                         # needs sage.rings.finite_rings
            sage: f = lambda x: k(x)
            sage: n = m.apply_map(f); n                                                 # needs sage.rings.finite_rings
            [0 1 2]
            [0 1 2]
            [0 1 2]
            sage: n.parent()                                                            # needs sage.rings.finite_rings
            Full MatrixSpace of 3 by 3 dense matrices
             over Finite Field in a of size 3^2

        In this example, we explicitly specify the codomain.

        ::

            sage: s = GF(3)
            sage: f = lambda x: s(x)
            sage: n = m.apply_map(f, k); n                                              # needs sage.rings.finite_rings
            [0 1 2]
            [0 1 2]
            [0 1 2]
            sage: n.parent()                                                            # needs sage.rings.finite_rings
            Full MatrixSpace of 3 by 3 dense matrices
             over Finite Field in a of size 3^2

        If ``self`` is subdivided, the result will be as well::

            sage: m = matrix(2, 2, srange(4))
            sage: m.subdivide(None, 1); m
            [0|1]
            [2|3]
            sage: m.apply_map(lambda x: x*x)
            [0|1]
            [4|9]

        If the matrix is sparse, the result will be as well::

            sage: m = matrix(ZZ,100,100,sparse=True)
            sage: m[18,32] = -6
            sage: m[1,83] = 19
            sage: n = m.apply_map(abs, R=ZZ)
            sage: n.dict()
            {(1, 83): 19, (18, 32): 6}
            sage: n.is_sparse()
            True

        If the map sends most of the matrix to zero, then it may be useful
        to get the result as a sparse matrix.

        ::

            sage: m = matrix(ZZ, 3, 3, range(1, 10))
            sage: n = m.apply_map(lambda x: 1//x, sparse=True); n
            [1 0 0]
            [0 0 0]
            [0 0 0]
            sage: n.parent()
            Full MatrixSpace of 3 by 3 sparse matrices over Integer Ring

        TESTS::

            sage: m = matrix([])
            sage: m.apply_map(lambda x: x*x) == m
            True

            sage: m.apply_map(lambda x: x*x, sparse=True).parent()
            Full MatrixSpace of 0 by 0 sparse matrices over Integer Ring

        Check that :issue:`19920` is fixed::

            sage: matrix.ones(2).apply_map(lambda x: int(-3))
            [-3 -3]
            [-3 -3]
        """
        if self._nrows == 0 or self._ncols == 0:
            if sparse is None or self.is_sparse() is sparse:
                return self.__copy__()
            elif sparse:
                return self.sparse_matrix()
            else:
                return self.dense_matrix()

        if self.is_sparse():
            values = {(i, j): phi(v) for (i, j), v in self.dict().iteritems()}
            if R is None:
                R = sage.structure.sequence.Sequence(values.values()).universe()
        else:
            values = [phi(v) for v in self.list()]
            if R is None:
                R = sage.structure.sequence.Sequence(values).universe()

        if isinstance(R, type):
            R = py_scalar_parent(R)
        if R not in Rings():
            raise TypeError("unable to find a common ring for all elements")

        if sparse is None or sparse is self.is_sparse():
            M = self.parent().change_ring(R)
        else:
            from sage.matrix.matrix_space import MatrixSpace
            M = MatrixSpace(R, self._nrows,
                            self._ncols, sparse=sparse)
        image = M(values)
        if self._subdivisions is not None:
            image.subdivide(*self.subdivisions())
        return image

    def characteristic_polynomial(self, *args, **kwds):
        """
        Synonym for self.charpoly(...).

        EXAMPLES::

            sage: a = matrix(QQ, 2,2, [1,2,3,4]); a
            [1 2]
            [3 4]
            sage: a.characteristic_polynomial('T')                                      # needs sage.libs.pari
            T^2 - 5*T - 2
        """
        return self.charpoly(*args, **kwds)

    def minimal_polynomial(self, var='x', **kwds):
        r"""
        This is a synonym for ``self.minpoly``.

        EXAMPLES::

            sage: a = matrix(QQ, 4, 4, range(16))
            sage: a.minimal_polynomial('z')                                             # needs sage.libs.pari
            z^3 - 30*z^2 - 80*z
            sage: a.minpoly()                                                           # needs sage.libs.pari
            x^3 - 30*x^2 - 80*x
        """
        return self.minpoly(var, **kwds)

    def minpoly(self, var='x', **kwds):
        r"""
        Return the minimal polynomial of ``self``.

        This uses a simplistic - and potentially very very slow - algorithm
        that involves computing kernels to determine the powers of the
        factors of the charpoly that divide the minpoly.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: A = matrix(GF(9, 'c'), 4, [1,1,0,0, 0,1,0,0, 0,0,5,0, 0,0,0,5])
            sage: factor(A.minpoly())
            (x + 1) * (x + 2)^2
            sage: A.minpoly()(A) == 0
            True
            sage: factor(A.charpoly())
            (x + 1)^2 * (x + 2)^2

        The default variable name is `x`, but you can specify
        another name::

            sage: factor(A.minpoly('y'))                                                # needs sage.rings.finite_rings
            (y + 1) * (y + 2)^2
        """
        f = self.fetch('minpoly')
        if f is not None:
            return f.change_variable_name(var)
        f = self.charpoly(var=var, **kwds)
        try:
            no_sq = f.is_squarefree()
            # is_squarefree for polys much faster than factor.
        except (TypeError, NotImplementedError):
            pass
        else:
            if no_sq:
                # Then f must be the minpoly
                self.cache('minpoly', f)
                return f

        # Now we have to work harder.  We find the power of each
        # irreducible factor that divides the minpoly.
        try:
            mp = f.radical()
        except (TypeError, NotImplementedError):
            raise NotImplementedError("minimal polynomial not implemented")
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

    def _test_minpoly(self, **options):
        """
        Check that :meth:`minpoly` works.

        EXAMPLES::

            sage: a = matrix([[1,2], [3,4]])
            sage: a._test_minpoly()                                                     # needs sage.libs.pari
        """
        if self.nrows() == self.ncols() and self.base_ring().is_exact():
            tester = self._tester(**options)
            # At least check that the minimal polynomial kills the matrix
            tester.assertTrue(self.minpoly().subs(x=self).is_zero())

    def charpoly(self, var='x', algorithm=None):
        r"""
        Return the characteristic polynomial of self, as a polynomial over
        the base ring.

        ALGORITHM:

        If the base ring has a method ``_matrix_charpoly``, we use it.

        In the generic case of matrices over a ring (commutative and with
        unity), there is a division-free algorithm, which can be accessed
        using ``'df'``, with complexity `O(n^4)`.  Alternatively, by
        specifying ``'hessenberg'``, this method computes the Hessenberg
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

        - ``var`` -- a variable name (default: ``'x'``)
        - ``algorithm`` -- string; one of
            - ``'df'`` -- generic `O(n^4)` division-free algorithm
            - ``'hessenberg'`` -- use the Hessenberg form of the matrix

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

            sage: A = MatrixSpace(QQ, 3)(range(9))
            sage: A.charpoly('x')                                                       # needs sage.libs.pari
            x^3 - 12*x^2 - 18*x
            sage: A.trace()
            12
            sage: A.determinant()
            0

        We compute the characteristic polynomial of a matrix over the
        polynomial ring `\ZZ[a]`::

            sage: R.<a> = PolynomialRing(ZZ)
            sage: M = MatrixSpace(R, 2)([a,1,  a,a+1]); M
            [    a     1]
            [    a a + 1]
            sage: f = M.charpoly('x'); f
            x^2 + (-2*a - 1)*x + a^2
            sage: f.parent()
            Univariate Polynomial Ring in x
             over Univariate Polynomial Ring in a over Integer Ring
            sage: M.trace()
            2*a + 1
            sage: M.determinant()
            a^2

        We compute the characteristic polynomial of a matrix over the
        multi-variate polynomial ring `\ZZ[x,y]`::

            sage: R.<x,y> = PolynomialRing(ZZ,2)
            sage: A = MatrixSpace(R, 2)([x, y, x^2, y^2])
            sage: f = A.charpoly('x'); f
            x^2 + (-y^2 - x)*x - x^2*y + x*y^2

        It's a little difficult to distinguish the variables. To fix this,
        we temporarily view the indeterminate as `Z`::

            sage: with localvars(f.parent(), 'Z'): print(f)
            Z^2 + (-y^2 - x)*Z - x^2*y + x*y^2

        We could also compute f in terms of Z from the start::

            sage: A.charpoly('Z')
            Z^2 + (-y^2 - x)*Z - x^2*y + x*y^2

        Here is an example over a number field::

            sage: # needs sage.rings.number_field
            sage: x = QQ['x'].gen()
            sage: K.<a> = NumberField(x^2 - 2)
            sage: m = matrix(K, [[a-1, 2], [a, a+1]])
            sage: m.charpoly('Z')
            Z^2 - 2*a*Z - 2*a + 1
            sage: m.charpoly('a')(m) == 0
            True

        Over integers modulo `n` with composite `n`::

            sage: A = Mat(Integers(6), 3, 3)(range(9))
            sage: A.charpoly()
            x^3

        Here is an example over a general commutative ring, that is to say,
        as of version 4.0.2, Sage does not even positively determine that
        ``S`` in the following example is an integral domain.  But the
        computation of the characteristic polynomial succeeds as follows::

            sage: # needs sage.libs.singular
            sage: R.<a,b> = QQ[]
            sage: S.<x,y> = R.quo((b^3))
            sage: A = matrix(S, [[x*y^2, 2*x], [2, x^10*y]]); A
            [ x*y^2    2*x]
            [     2 x^10*y]
            sage: A.charpoly('T')
            T^2 + (-x^10*y - x*y^2)*T - 4*x

        TESTS::

            sage: P.<a,b,c> = PolynomialRing(Rationals())
            sage: u = MatrixSpace(P, 3)([[0,0,a], [1,0,b], [0,1,c]])
            sage: Q.<x> = PolynomialRing(P)
            sage: u.charpoly('x')
            x^3 - c*x^2 - b*x - a

        A test case from :issue:`6442`. Prior to :issue:`12292`, the
        call to ``A.det()`` would attempt to use the cached charpoly,
        and crash if an empty dictionary was cached. We don't cache
        dictionaries anymore, but this test should still pass::

            sage: # needs sage.rings.padics
            sage: z = Zp(p=5)
            sage: A = matrix(z, [ [3 + O(5^1), 4 + O(5^1), 4 + O(5^1)],
            ....:                 [2*5^2 + O(5^3), 2 + O(5^1), 1 + O(5^1)],
            ....:                 [5 + O(5^2), 1 + O(5^1), 1 + O(5^1)] ])
            sage: A.charpoly(algorithm='hessenberg')
            Traceback (most recent call last):
            ...
            ValueError: negative valuation
            sage: A.det()
            3 + O(5)

        The cached polynomial should be independent of the ``var``
        argument (:issue:`12292`). We check (indirectly) that the
        second call uses the cached value by noting that its result is
        not cached::

            sage: M = MatrixSpace(RR, 2)
            sage: A = M(range(2^2))
            sage: type(A)
            <class 'sage.matrix.matrix_generic_dense.Matrix_generic_dense'>
            sage: A.charpoly('x')
            x^2 - 3.00000000000000*x - 2.00000000000000
            sage: A.charpoly('y')
            y^2 - 3.00000000000000*y - 2.00000000000000
            sage: A._cache['charpoly']
            x^2 - 3.00000000000000*x - 2.00000000000000
        """
        f = self.fetch('charpoly')
        if f is not None:
            return f.change_variable_name(var)

        R = self._base_ring

        if algorithm is None:
            if hasattr(R, '_matrix_charpoly'):
                f = R._matrix_charpoly(self, var)
            if f is None:
                if R in _Fields and R.is_exact():
                    f = self._charpoly_hessenberg(var)
                else:
                    f = self._charpoly_df(var)
        else:
            if algorithm == "hessenberg":
                f = self._charpoly_hessenberg(var)
            elif algorithm == "df":
                f = self._charpoly_df(var)
            else:
                raise ValueError('algorithm must be "hessenberg" or "df"')

        # Cache the result, and return it.
        self.cache('charpoly', f)
        return f

    def _charpoly_df(self, var='x'):
        r"""
        Compute the characteristic polynomial of ``self`` without divisions.

        INPUT:

        - ``var`` -- a variable name (default: ``'x'``)

        OUTPUT: polynomial; the characteristic polynomial of ``self``

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
            sage: A = matrix(R,[[7*t^2 - t - 9, -1/4*t - 1, -17*t^2 - t + 1],
            ....:               [-t^2 + 1/4*t, t^2 + 5/7*t + 3, 1/5*t^2 + 1662],
            ....:               [-2*t - 3, 2*t^2 + 6*t - 1/2, -1/6*t^2]])
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

        Test that :issue:`27937` is fixed::

            sage: R = FreeAbelianMonoid('u,v').algebra(QQ)                              # needs sage.combinat
            sage: matrix(4, 4, lambda i, j: R.an_element())._charpoly_df()              # needs sage.combinat
            B[1]*x^4 - 4*B[u]*x^3

        .. NOTE::

            The key feature of this implementation is that it is division-free.
            This means that it can be used as a generic implementation for any
            ring (commutative and with multiplicative identity).  The algorithm
            is described in full detail as Algorithm 3.1 in [Sei2002]_.

            Note that there is a missing minus sign in front of the last term in
            the penultimate line of Algorithm 3.1.
        """

        # Validate assertions
        if not self.is_square():
            raise ValueError("self must be a square matrix")

        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

        # Extract parameters
        cdef Matrix M  = <Matrix> self
        n  = M._ncols
        R  = M._base_ring
        S  = PolynomialRing(R, var)

        # Corner cases
        # N.B.  We already tested for M to be square, hence we do not need to
        # test for 0 x n or m x 0 matrices.
        if n == 0:
            return S.one()

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
        # having been ported from Magma to Sage, is 0-based.
        from sage.matrix.constructor import matrix

        F = [R.zero()] * n
        cdef Matrix a = <Matrix> matrix(R, n-1, n)
        A = [R.zero()] * n

        F[0] = - M.get_unsafe(0, 0)
        for t in range(1, n):

            # Set a(1, t) to be M(<=t, t)
            for i in range(t+1):
                a.set_unsafe(0, i, M.get_unsafe(i, t))

            # Set A[1, t] to be the (t)th entry in a[1, t]
            A[0] = M.get_unsafe(t, t)

            for p in range(1, t):

                # Set a(p, t) to the product of M[<=t, <=t] * a(p-1, t)
                for i in range(t+1):
                    s = R.zero()
                    for j in range(t+1):
                        s = s + M.get_unsafe(i, j) * a.get_unsafe(p-1, j)
                    a.set_unsafe(p, i, s)

                # Set A[p, t] to be the (t)th entry in a[p, t]
                A[p] = a.get_unsafe(p, t)

            # Set A[t, t] to be M[t, <=t] * a(p-1, t)
            s = R.zero()
            for j in range(t+1):
                s = s + M.get_unsafe(t, j) * a.get_unsafe(t-1, j)
            A[t] = s

            for p in range(t+1):
                s = F[p]
                for k in range(p):
                    s = s - A[k] * F[p-k-1]
                F[p] = s - A[p]

        X = S.gen(0)
        f = X**n + S(list(reversed(F)))

        return f

    def fcp(self, var='x'):
        """
        Return the factorization of the characteristic polynomial of ``self``.

        INPUT:

        - ``var`` -- (default: ``'x'``) name of variable of charpoly

        EXAMPLES::

            sage: M = MatrixSpace(QQ, 3, 3)
            sage: A = M([1,9,-7, 4/5,4,3, 6,4,3])
            sage: A.fcp()                                                               # needs sage.libs.pari
            x^3 - 8*x^2 + 209/5*x - 286
            sage: A = M([3, 0, -2, 0, -2, 0, 0, 0, 0])
            sage: A.fcp('T')                                                            # needs sage.libs.pari
            (T - 3) * T * (T + 2)
        """
        return self.charpoly(var).factor()

    def denominator(self):
        r"""
        Return the least common multiple of the denominators of the
        elements of ``self``.

        If there is no denominator function for the base field, or no LCM
        function for the denominators, raise a :exc:`TypeError`.

        EXAMPLES::

            sage: A = MatrixSpace(QQ, 2)([1/2, 1/3, 1/5, 1/7])
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

            sage: # needs sage.rings.number_field
            sage: K.<z> = CyclotomicField(3)
            sage: M = MatrixSpace(K, 3, sparse=True)
            sage: A = M([(1+z)/3, (2+z)/3, z/3, 1, 1+z, -2, 1, 5, -1+z])
            sage: print(A)
            [1/3*z + 1/3 1/3*z + 2/3       1/3*z]
            [          1       z + 1          -2]
            [          1           5       z - 1]
            sage: print(A.denominator())
            3
        """
        if self.nrows() == 0 or self.ncols() == 0:
            return ZZ.one()
        x = self.list()
        try:
            d = x[0].denominator()
        except AttributeError:
            raise TypeError("denominator not defined for elements of the base ring")
        try:
            for y in x:
                d = d.lcm(y.denominator())
        except AttributeError:
            raise TypeError("lcm function not defined for elements of the base ring")
        return d

    def diagonal(self):
        r"""
        Return the diagonal entries of ``self``.

        OUTPUT:

        A list containing the entries of the matrix that
        have equal row and column indices, in order of the
        indices.  Behavior is not limited to square matrices.

        EXAMPLES::

            sage: A = matrix([[2,5], [3,7]]); A
            [2 5]
            [3 7]
            sage: A.diagonal()
            [2, 7]

        Two rectangular matrices. ::

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
        return [self[i, i] for i in range(n)]

    def trace(self):
        """
        Return the trace of self, which is the sum of the diagonal entries
        of ``self``.

        INPUT:

        - ``self`` -- a square matrix

        OUTPUT: element of the base ring of self

        EXAMPLES::

            sage: a = matrix(3, 3, range(9)); a
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: a.trace()
            12
            sage: a = matrix({(1,1): 10, (2,1): -3, (2,2): 4/3}); a
            [  0   0   0]
            [  0  10   0]
            [  0  -3 4/3]
            sage: a.trace()
            34/3
        """
        if self._nrows != self._ncols:
            raise ValueError("self must be a square matrix")
        R = self._base_ring
        cdef Py_ssize_t i
        cdef object s
        s = R(0)
        for i from 0 <= i < self._nrows:
            s = s + self.get_unsafe(i, i)
        return s

    def trace_of_product(self, Matrix other):
        """
        Return the trace of ``self * other`` without computing the entire product.

        EXAMPLES::

            sage: M = random_matrix(ZZ, 10, 20)
            sage: N = random_matrix(ZZ, 20, 10)
            sage: M.trace_of_product(N) == (M*N).trace()
            True
        """
        if self._nrows != other._ncols or other._nrows != self._ncols:
            raise ArithmeticError("incompatible dimensions")
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
        Return Hessenberg form of ``self``.

        If the base ring is merely an integral domain (and not a field),
        the Hessenberg form will (in general) only be defined over the
        fraction field of the base ring.

        EXAMPLES::

            sage: A = matrix(ZZ, 4, [2, 1, 1, -2, 2, 2, -1, -1, -1,1,2,3,4,5,6,7])
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
        if X is not None:
            return X
        R = self._base_ring
        if R not in _Fields:
            try:
                K = self._base_ring.fraction_field()
                H = self.change_ring(K)
                H.hessenbergize()
            except TypeError as msg:
                raise TypeError("%s\nHessenberg form only possible for matrices over a field"%msg)
        else:
            H = self.__copy__()
            H.hessenbergize()
        # end if
        self.cache('hessenberg_form', H)
        return H

    def hessenbergize(self):
        """
        Transform ``self`` to Hessenberg form.

        The hessenberg form of a matrix `A` is a matrix that is
        similar to `A`, so has the same characteristic polynomial
        as `A`, and is upper triangular except possible for entries
        right below the diagonal.

        ALGORITHM: See Henri Cohen's first book.

        EXAMPLES::

            sage: A = matrix(QQ, 3, [2, 1, 1, -2, 2, 2, -1, -1, -1])
            sage: A.hessenbergize(); A
            [  2 3/2   1]
            [ -2   3   2]
            [  0  -3  -2]

        ::

            sage: A = matrix(QQ, 4, [2, 1, 1, -2, 2, 2, -1, -1, -1,1,2,3,4,5,6,7])
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

        tm = verbose(f"Computing Hessenberg Normal Form of {n}x{n} matrix",
                     level=2)

        if not self.is_square():
            raise TypeError("self must be square")

        self.check_mutability()

        base = self._base_ring
        if hasattr(base, '_matrix_hessenbergize'):
            base._matrix_hessenbergize(self)
            return

        if self._base_ring not in _Fields:
            raise TypeError("Hessenbergize only possible for matrices over a field")

        zero = self._base_ring(0)
        one = self._base_ring(1)
        for m from 1 <= m < n-1:
            # Search for a nonzero entry in column m-1
            i = -1
            for r from m+1 <= r < n:
                if not self.get_is_zero_unsafe(r, m-1):
                    i = r
                    break
            if i != -1:
                # Found a nonzero entry in column m-1 that is strictly below row m
                # Now set i to be the first nonzero position >= m in column m-1
                if not self.get_is_zero_unsafe(m, m-1):
                    i = m
                t = self.get_unsafe(i, m-1)
                t_inv = None
                if i > m:
                    self.swap_rows_c(i, m)
                    # We must do the corresponding column swap to
                    # maintain the characteristic polynomial (which is
                    # an invariant of Hessenberg form)
                    self.swap_columns_c(i, m)
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
        verbose(f"Finished Hessenberg Normal Form of {n}x{n} matrix",
                level=2, t=tm)

    def _charpoly_hessenberg(self, var):
        """
        Transform ``self`` in place to its Hessenberg form then computes and
        returns the coefficients of the characteristic polynomial of this
        matrix.

        INPUT:

        - ``var`` -- name of the indeterminate of the charpoly

        The characteristic polynomial is represented as a vector of ints,
        where the constant term of the characteristic polynomial is the `0`-th
        coefficient of the vector.

        EXAMPLES::

            sage: matrix(QQ, 3,3, range(9))._charpoly_hessenberg('Z')
            Z^3 - 12*Z^2 - 18*Z
            sage: matrix(ZZ, 3,3, range(9))._charpoly_hessenberg('Z')
            Z^3 - 12*Z^2 - 18*Z
            sage: matrix(GF(7), 3, 3, range(9))._charpoly_hessenberg('Z')
            Z^3 + 2*Z^2 + 3*Z
            sage: matrix(QQ['x'], 3,3, range(9))._charpoly_hessenberg('Z')
            Z^3 - 12*Z^2 - 18*Z
            sage: matrix(ZZ['ZZ'], 3,3, range(9))._charpoly_hessenberg('Z')
            Z^3 - 12*Z^2 - 18*Z
        """
        if self._nrows != self._ncols:
            raise ArithmeticError("charpoly not defined for non-square matrix.")

        # Replace ``self`` by its Hessenberg form
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

        cdef Py_ssize_t i, m, n
        n = self._nrows

        cdef Matrix c
        c = H.new_matrix(nrows=n+1, ncols=n+1)    # the 0 matrix
        one = H._coerce_element(1)
        c.set_unsafe(0, 0, one)

        for m from 1 <= m <= n:
            # Set the m-th row of c to (x - H[m-1,m-1])*c[m-1] = x*c[m-1] - H[m-1,m-1]*c[m-1]
            # We do this by hand by setting the m-th row to c[m-1]
            # shifted to the right by one.  We then add
            # -H[m-1,m-1]*c[m-1] to the resulting m-th row.
            for i from 1 <= i <= n:
                c.set_unsafe(m, i, c.get_unsafe(m-1, i-1))
            c.add_multiple_of_row_c(m, m-1, -H.get_unsafe(m-1, m-1), 0)
            t = one
            for i from 1 <= i < m:
                t = t * H.get_unsafe(m-i, m-i-1)
                # Set the m-th row of c to c[m] - t*H[m-i-1,m-1]*c[m-i-1]
                c.add_multiple_of_row_c(m, m-i-1, - t*H.get_unsafe(m-i-1, m-1), 0)

        # The answer is now the n-th row of c.
        v = PyList_New(n+1)     # this is really sort of v = []..."
        for i from 0 <= i <= n:
            # Finally, set v[i] = c[n,i]
            o = c.get_unsafe(n, i)
            Py_INCREF(o)
            PyList_SET_ITEM(v, i, o)

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

            sage: M = Matrix(QQ, [[1,0,0,1], [0,1,1,0], [1,1,1,0]])
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

            sage: A = matrix(ZZ,3,3,range(9))
            sage: A.right_nullity()
            1
        """
        return self.ncols() - self.rank()

    #####################################################################################
    # Kernel Helper Functions
    #####################################################################################

    def _right_kernel_matrix_over_number_field(self):
        r"""
        Return a pair that includes a matrix of basis vectors
        for the right kernel of ``self``.

        OUTPUT:

        Returns a pair.  First item is the string ``'pivot-pari-numberfield'``
        that identifies the nature of the basis vectors.

        Second item is a matrix whose rows are a basis for the right kernel,
        over the number field, as computed by PARI.

        EXAMPLES::

            sage: # needs sage.libs.pari sage.rings.number_field
            sage: Q = QuadraticField(-7)
            sage: a = Q.gen(0)
            sage: A = matrix(Q, [[    2, 5 - a,   15 - a],
            ....:                [2 + a,     a, -7 + 5*a]])
            sage: result = A._right_kernel_matrix_over_number_field()
            sage: result[0]
            'pivot-pari-numberfield'
            sage: P = result[1]; P
            [-a -3  1]
            sage: A * P.transpose() == zero_matrix(Q, 2, 1)
            True

        TESTS:

        We test some trivial cases. ::

            sage: # needs sage.libs.pari sage.rings.number_field
            sage: Q = QuadraticField(-7)
            sage: A = matrix(Q, 0, 2)
            sage: A._right_kernel_matrix_over_number_field()[1]
            [1 0]
            [0 1]
            sage: A = matrix(Q, 2, 0)
            sage: A._right_kernel_matrix_over_number_field()[1].parent()
            Full MatrixSpace of 0 by 0 dense matrices
             over Number Field in a with defining polynomial x^2 + 7
              with a = 2.645751311064591?*I
            sage: A = zero_matrix(Q, 4, 3)
            sage: A._right_kernel_matrix_over_number_field()[1]
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        from sage.matrix.matrix_space import MatrixSpace
        tm = verbose("computing right kernel matrix over a number field for %sx%s matrix" % (self.nrows(), self.ncols()), level=2)
        basis = self.__pari__().matker()
        # Coerce PARI representations into the number field
        R = self.base_ring()
        basis = [[R(x) for x in row] for row in basis]
        verbose("done computing right kernel matrix over a number field for %sx%s matrix" % (self.nrows(), self.ncols()), level=2, t=tm)
        return 'pivot-pari-numberfield', MatrixSpace(R, len(basis), ncols=self._ncols)(basis)

    def _right_kernel_matrix_over_field(self, *args, **kwds):
        r"""
        Return a pair that includes a matrix of basis vectors
        for the right kernel of ``self``.

        OUTPUT:

        Returns a pair.  First item is the string 'pivot-generic'
        that identifies the nature of the basis vectors.

        Second item is a matrix whose rows are a basis for the right kernel,
        over the field, as computed by general Python code.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: C = CyclotomicField(14)
            sage: a = C.gen(0)
            sage: A = matrix(C, 3, 4, [[  1,    a,    1+a,  a^3+a^5],
            ....:                      [  a,  a^4,  a+a^4,  a^4+a^8],
            ....:                      [a^2, a^6, a^2+a^6, a^5+a^10]])
            sage: result = A._right_kernel_matrix_over_field()
            sage: result[0]
            'pivot-generic'
            sage: P = result[1]; P
            [       -1        -1         1         0]
            [-zeta14^3 -zeta14^4         0         1]
            sage: A * P.transpose() == zero_matrix(C, 3, 2)
            True

        TESTS:

        We test some trivial cases. ::

            sage: # needs sage.rings.number_field
            sage: C = CyclotomicField(14)
            sage: A = matrix(C, 0, 2)
            sage: A._right_kernel_matrix_over_field()[1]
            [1 0]
            [0 1]
            sage: A = matrix(C, 2, 0)
            sage: A._right_kernel_matrix_over_field()[1].parent()
            Full MatrixSpace of 0 by 0 dense matrices
             over Cyclotomic Field of order 14 and degree 6
            sage: A = zero_matrix(C, 4, 3)
            sage: A._right_kernel_matrix_over_field()[1]
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        from sage.matrix.matrix_space import MatrixSpace
        tm = verbose("computing right kernel matrix over an arbitrary field for %sx%s matrix" % (self.nrows(), self.ncols()), level=2)
        E = self.echelon_form(*args, **kwds)
        pivots = E.pivots()
        pivots_set = set(pivots)
        zero = self._base_ring.zero()
        one = self._base_ring.one()
        MS = self.matrix_space(self._ncols-len(pivots), self._ncols)
        cdef Py_ssize_t i, r, cur_row
        if self.is_sparse():
            entries = {}
            cur_row = 0
            for i in range(self._ncols):
                if i not in pivots_set:
                    entries[cur_row, i] = one
                    for r, p in enumerate(pivots):
                        entries[cur_row, p] = -E[r, i]
                    cur_row += 1
            M = MS(entries, copy=False, coerce=False)
        else:
            basis = []
            for i in range(self._ncols):
                if i not in pivots_set:
                    v = [zero] * self._ncols
                    v[i] = one
                    for r, p in enumerate(pivots):
                        v[p] = -E[r, i]
                    basis.append(v)
            M = MS(basis, coerce=False)
        tm = verbose("done computing right kernel matrix over an arbitrary field for %sx%s matrix"
                     % (self.nrows(), self.ncols()), level=2, t=tm)
        return 'pivot-generic', M

    def _right_kernel_matrix_over_domain(self):
        r"""
        Return a pair that includes a matrix of basis vectors
        for the right kernel of ``self``.

        OUTPUT:

        Returns a pair.  First item is the string 'computed-smith-form'
        that identifies the nature of the basis vectors.

        Second item is a matrix whose rows are a basis for the right kernel,
        over the field, as computed by general Python code.

        .. WARNING::

            This routine uses Smith normal form, which can fail
            if the domain is not a principal ideal domain.  Since we do
            not have a good test for PIDs, this is just a risk we take.
            See an example failure in the documentation for
            :meth:`right_kernel_matrix`.

        EXAMPLES:

        Univariate polynomials over a field form a PID::

            sage: R.<y> = QQ[]
            sage: A = matrix(R, [[  1,   y, 1+y^2],
            ....:                [y^3, y^2, 2*y^3]])
            sage: result = A._right_kernel_matrix_over_domain()
            sage: result[0]
            'computed-smith-form'
            sage: P = result[1]; P
            [-1 -y  1]
            sage: A * P.transpose() == zero_matrix(R, 2, 1)
            True

        TESTS:

        We test some trivial cases::

            sage: R.<y> = QQ[]
            sage: A = matrix(R, 0, 2)
            sage: A._right_kernel_matrix_over_domain()[1]
            [1 0]
            [0 1]
            sage: A = matrix(R, 2, 0)
            sage: A._right_kernel_matrix_over_domain()[1].parent()
            Full MatrixSpace of 0 by 0 dense matrices
             over Univariate Polynomial Ring in y over Rational Field
            sage: A = zero_matrix(R, 4, 3)
            sage: A._right_kernel_matrix_over_domain()[1]
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        tm = verbose("computing right kernel matrix over a domain for %sx%s matrix"
                     % (self.nrows(), self.ncols()), level=2)
        d, _, v = self.smith_form()
        basis = []
        cdef Py_ssize_t i, nrows = self._nrows
        for i in range(self._ncols):
            if i >= nrows or d[i, i] == 0:
                basis.append(v.column(i))
        verbose("done computing right kernel matrix over a domain for %sx%s matrix"
                % (self.nrows(), self.ncols()), level=2, t=tm)
        return 'computed-smith-form', self.new_matrix(nrows=len(basis), ncols=self._ncols, entries=basis)

    def _right_kernel_matrix_over_integer_mod_ring(self):
        r"""
        Return a pair that includes a matrix of basis vectors
        for the right kernel of ``self``.

        OUTPUT:

        Returns a pair. First item is the string 'computed-pari-matkermod'
        that identifies the nature of the basis vectors.

        Second item is a matrix whose rows are a basis for the right kernel,
        over the integer modular ring, as computed by :pari:`matkermod`.

        EXAMPLES::

            sage: # needs sage.libs.pari
            sage: A = matrix(Zmod(24480), [[1,2,3,4,5], [7,7,7,7,7]])
            sage: result = A._right_kernel_matrix_over_integer_mod_ring()
            sage: result[0]
            'computed-pari-matkermod'
            sage: P = result[1]; P
            [    1 24478     1     0     0]
            [    2 24477     0     1     0]
            [    3 24476     0     0     1]
            sage: A * P.transpose() == 0
            True
        """
        R = self.base_ring()
        pariM = self.change_ring(ZZ).__pari__()
        basis = list(pariM.matkermod(R.characteristic()))
        from sage.matrix.matrix_space import MatrixSpace
        return 'computed-pari-matkermod', MatrixSpace(R, len(basis), ncols=self._ncols)(basis)

    def right_kernel_matrix(self, *args, **kwds):
        r"""
        Return a matrix whose rows form a basis
        for the right kernel of ``self``.

        INPUT:

        - ``algorithm`` -- (default: ``'default'``) a keyword that selects the
          algorithm employed.  Allowable values are:

          - ``'default'`` -- allows the algorithm to be chosen automatically
          - ``'generic'`` -- naive algorithm usable for matrices over any field
          - ``'flint'`` -- FLINT library code for matrices over the rationals
            or the integers
          - ``'pari'`` -- PARI library code for matrices over number fields
            or the integers
          - ``'padic'`` -- padic algorithm from IML library for matrices
            over the rationals and integers
          - ``'pluq'`` -- PLUQ matrix factorization for matrices mod 2

        - ``basis`` -- (default: ``'default'``) a keyword that describes
          the format of the basis returned.  Allowable values are:

          - 'default': uses 'echelon' over fields; 'computed' otherwise.
          - 'echelon': the basis matrix is returned in echelon form.
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
            of the basis. This option is recommended for inexact rings.

        OUTPUT:

        A matrix ``X``  whose rows are an independent set spanning the
        right kernel of ``self``.  So ``self*X.transpose()`` is a zero matrix.

        The output varies depending on the choice of ``algorithm`` and the
        format chosen by ``basis``.

        The results of this routine are not cached, so you can call it again
        with different options to get possibly different output (like the basis
        format).  Conversely, repeated calls on the same matrix will always
        start from scratch.

        .. NOTE::

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
            ....:                 [-5, 1, 0, 7, -3],
            ....:                 [0, -1, -4, 6, -2],
            ....:                 [4, -1, 0, -6, 2]])
            sage: C = A.right_kernel_matrix(algorithm='default', basis='computed'); C
            [-1  2 -2 -1  0]
            [ 1  2  0  0 -1]
            sage: A * C.transpose() == zero_matrix(QQ, 4, 2)
            True
            sage: P = A.right_kernel_matrix(algorithm='padic', basis='pivot'); P
            [ 1 -2  2  1  0]
            [-1 -2  0  0  1]
            sage: A * P.transpose() == zero_matrix(QQ, 4, 2)
            True
            sage: C == -P
            True
            sage: E = A.right_kernel_matrix(algorithm='default', basis='echelon'); E
            [   1    0    1  1/2 -1/2]
            [   0    1 -1/2 -1/4 -1/4]
            sage: A * E.transpose() == zero_matrix(QQ, 4, 2)
            True

        Since the rationals are a field, we can call the general code
        available for any field by using the 'generic' keyword. ::

            sage: A = matrix(QQ, [[1, 0, 1, -3, 1],
            ....:                 [-5, 1, 0, 7, -3],
            ....:                 [0, -1, -4, 6, -2],
            ....:                 [4, -1, 0, -6, 2]])
            sage: G = A.right_kernel_matrix(algorithm='generic', basis='echelon'); G
            [   1    0    1  1/2 -1/2]
            [   0    1 -1/2 -1/4 -1/4]
            sage: A * G.transpose() == zero_matrix(QQ, 4, 2)
            True

        We verify that the rational matrix code is called for both
        dense and sparse rational matrices, with equal result. ::

            sage: A = matrix(QQ, [[1, 0, 1, -3, 1],
            ....:                 [-5, 1, 0, 7, -3],
            ....:                 [0, -1, -4, 6, -2],
            ....:                 [4, -1, 0, -6, 2]],
            ....:            sparse=False)
            sage: B = copy(A).sparse_matrix()
            sage: from sage.misc.verbose import set_verbose
            sage: set_verbose(2)
            sage: D = A.right_kernel(); D
            verbose 2 (<module>) computing a right kernel for 4x5 matrix over Rational Field
            ...
            Vector space of degree 5 and dimension 2 over Rational Field
            Basis matrix:
            [   1    0    1  1/2 -1/2]
            [   0    1 -1/2 -1/4 -1/4]
            sage: S = B.right_kernel(); S
            verbose 2 (<module>) computing a right kernel for 4x5 matrix over Rational Field
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

            sage: # needs sage.rings.number_field
            sage: Q = QuadraticField(-7)
            sage: a = Q.gen(0)
            sage: A = matrix(Q, [[2, 5-a, 15-a, 16+4*a],
            ....:                [2+a, a, -7 + 5*a, -3+3*a]])
            sage: C = A.right_kernel_matrix(algorithm='default', basis='computed'); C
            [    -a     -3      1      0]
            [    -2 -a - 1      0      1]
            sage: A*C.transpose() == zero_matrix(Q, 2, 2)
            True
            sage: P = A.right_kernel_matrix(algorithm='pari', basis='pivot'); P         # needs sage.libs.pari
            [    -a     -3      1      0]
            [    -2 -a - 1      0      1]
            sage: A*P.transpose() == zero_matrix(Q, 2, 2)                               # needs sage.libs.pari
            True
            sage: E = A.right_kernel_matrix(algorithm='default', basis='echelon'); E
            [                1                 0     7/88*a + 3/88 -3/176*a - 39/176]
            [                0                 1   -1/88*a - 13/88  13/176*a - 7/176]
            sage: A*E.transpose() == zero_matrix(Q, 2, 2)
            True

        We can bypass using PARI for number fields and use Sage's general
        code for matrices over any field.  The basis vectors as computed
        are in pivot format. ::

            sage: # needs sage.rings.number_field
            sage: Q = QuadraticField(-7)
            sage: a = Q.gen(0)
            sage: A = matrix(Q, [[2, 5-a, 15-a, 16+4*a], [2+a, a, -7 + 5*a, -3+3*a]])
            sage: G = A.right_kernel_matrix(algorithm='generic', basis='computed'); G
            [    -a     -3      1      0]
            [    -2 -a - 1      0      1]
            sage: A*G.transpose() == zero_matrix(Q, 2, 2)
            True

        We check that number fields are handled by the right routine as part of
        typical right kernel computation. ::

            sage: # needs sage.rings.number_field
            sage: Q = QuadraticField(-7)
            sage: a = Q.gen(0)
            sage: A = matrix(Q, [[2, 5-a, 15-a, 16+4*a], [2+a, a, -7 + 5*a, -3+3*a]])
            sage: set_verbose(2)
            sage: A.right_kernel(algorithm='default')
            verbose ...
            verbose 2 (<module>) computing right kernel matrix over a number field for 2x4 matrix
            verbose 2 (<module>) done computing right kernel matrix over a number field for 2x4 matrix
            ...
            Vector space of degree 4 and dimension 2 over
             Number Field in a with defining polynomial x^2 + 7 with a = 2.645751311064591?*I
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
            ....:                   [1, 0, 0, 0, 1, 1,],
            ....:                   [1, 0, 0, 0, 1, 1]])
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

            sage: A = matrix(GF(2), [[0, 1, 1, 0, 0, 0],
            ....:                    [1, 0, 0, 0, 1, 1,],
            ....:                    [1, 0, 0, 0, 1, 1]])
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
            sage: A.right_kernel_matrix(algorithm='pluq', basis='echelon')              # needs sage.libs.m4ri
            [1 0 0 0 0 1]
            [0 1 1 0 0 0]
            [0 0 0 1 0 0]
            [0 0 0 0 1 1]

        We test that the mod 2 code is called for matrices over GF(2). ::

            sage: A = matrix(GF(2), [[0, 1, 1, 0, 0, 0],
            ....:                    [1, 0, 0, 0, 1, 1,],
            ....:                    [1, 0, 0, 0, 1, 1]])
            sage: set_verbose(2)                                                        # needs sage.rings.finite_rings
            sage: A.right_kernel(algorithm='default')
            verbose ...
            verbose ... (<module>) computing right kernel matrix over integers mod 2 for 3x6 matrix
            verbose ... (<module>) done computing right kernel matrix over integers mod 2 for 3x6 matrix
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

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = FiniteField(5^2)
            sage: A = matrix(F, 3, 4, [[  1,   a,     1+a,  a^3+a^5],
            ....:                      [  a, a^4,   a+a^4,  a^4+a^8],
            ....:                      [a^2, a^6, a^2+a^6, a^5+a^10]])
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
            ....:                       [2,0,2,2],
            ....:                       [-1,1,-2,0]])
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

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = FiniteField(5^2)
            sage: A = matrix(F, 3, 4, [[  1,   a,     1+a,  a^3+a^5],
            ....:                      [  a, a^4,   a+a^4,  a^4+a^8],
            ....:                      [a^2, a^6, a^2+a^6, a^5+a^10]])
            sage: set_verbose(2)
            sage: A.right_kernel(algorithm='default')
            verbose ...
            verbose 2 (<module>) computing right kernel matrix over an arbitrary field for 3x4 matrix
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
            ....:                 [4, 0, 3, 4, 2, 7, 7],
            ....:                 [1, 4, 6, 1, 2, 8, 5],
            ....:                 [0, 3, 1, 2, 3, 6, 2]])

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

            sage: # needs sage.libs.pari
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
            ....:                 [4, 0, 3, 4, 2, 7, 7],
            ....:                 [1, 4, 6, 1, 2, 8, 5],
            ....:                 [0, 3, 1, 2, 3, 6, 2]],
            ....:            sparse=False)
            sage: B = copy(A).sparse_matrix()
            sage: set_verbose(2)
            sage: D = A.right_kernel(); D
            verbose ... (<module>) computing a right kernel for 4x7 matrix over Integer Ring
            verbose ... (<module>) computing right kernel matrix over the integers for 4x7 matrix
            ...
            verbose ... (<module>) done computing right kernel matrix over the integers for 4x7 matrix
            ...
            Free module of degree 7 and rank 3 over Integer Ring
            Echelon basis matrix:
            [  1  12   3  14  -3 -10   1]
            [  0  35   0  25  -1 -31  17]
            [  0   0   7  12  -3  -1  -8]
            sage: S = B.right_kernel(); S
            verbose ... (<module>) computing a right kernel for 4x7 matrix over Integer Ring
            verbose ... (<module>) computing right kernel matrix over the integers for 4x7 matrix
            ...
            verbose ... (<module>) done computing right kernel matrix over the integers for 4x7 matrix
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
            ....:                [y^3, y^2, 2*y^3]])
            sage: E = A.right_kernel_matrix(algorithm='default', basis='echelon'); E
            [-1 -y  1]
            sage: A*E.transpose() == zero_matrix(ZZ, 2, 1)
            True

        It can be computationally expensive to determine if an integral
        domain is a principal ideal domain.  The Smith normal form routine
        can fail for non-PIDs, as in this example. ::

            sage: D.<x> = ZZ[]
            sage: A = matrix(D, 2, 2, [[x^2 - x, -x + 5],
            ....:                      [x^2 - 8, -x + 2]])
            sage: A.right_kernel_matrix()
            Traceback (most recent call last):
            ...
            ArithmeticError: Ideal Ideal (x^2 - x, x^2 - 8) of
            Univariate Polynomial Ring in x over Integer Ring not principal

        We test that the domain code is called for domains that lack any
        extra structure. ::

            sage: R.<y> = QQ[]
            sage: A = matrix(R, [[  1,   y, 1+y^2],
            ....:                [y^3, y^2, 2*y^3]])
            sage: set_verbose(2)
            sage: A.right_kernel(algorithm='default', basis='echelon')
            verbose ...
            verbose 2 (<module>) computing right kernel matrix over a domain for 2x3 matrix
            verbose 2 (<module>) done computing right kernel matrix over a domain for 2x3 matrix
            ...
            Free module of degree 3 and rank 1 over
             Univariate Polynomial Ring in y over Rational Field
            Echelon basis matrix:
            [-1 -y  1]
            sage: set_verbose(0)

        Over inexact rings:

        For inexact rings one should avoid echelonizing if possible::

            sage: A = Matrix(
            ....: [[          0.0,           0.5,  0.8090169944],
            ....:  [          0.0,           0.5, -0.8090169944],
            ....:  [          0.0,          -0.5,  0.8090169944],
            ....:  [          0.0,          -0.5, -0.8090169944],
            ....:  [          0.5,  0.8090169944,           0.0],
            ....:  [          0.5, -0.8090169944,           0.0],
            ....:  [         -0.5,  0.8090169944,           0.0],
            ....:  [         -0.5, -0.8090169944,           0.0],
            ....:  [ 0.8090169944,           0.0,           0.5],
            ....:  [-0.8090169944,           0.0,           0.5],
            ....:  [ 0.8090169944,           0.0,          -0.5],
            ....:  [-0.8090169944,           0.0,          -0.5]]).transpose()
            sage: (A * A.right_kernel_matrix().transpose()).norm() > 2
            True
            sage: (A * A.right_kernel_matrix(basis='computed').transpose()).norm() < 1e-15
            True

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
            NotImplementedError: Cannot compute a matrix kernel over
            Quaternion Algebra (-1, -1) with base ring Rational Field

        We test error messages for improper choices of the 'algorithm'
        keyword. ::

            sage: matrix(ZZ, 2, 2).right_kernel_matrix(algorithm='junk')
            Traceback (most recent call last):
            ...
            ValueError: matrix kernel algorithm 'junk' not recognized
            sage: matrix(GF(2), 2, 2).right_kernel_matrix(algorithm='padic')
            Traceback (most recent call last):
            ...
            ValueError: 'padic' matrix kernel algorithm only available over
            the rationals and the integers, not over Finite Field of size 2
            sage: matrix(QQ, 2, 2).right_kernel_matrix(algorithm='pari')
            Traceback (most recent call last):
            ...
            ValueError: 'pari' matrix kernel algorithm only available over
            non-trivial number fields and the integers, not over Rational Field
            sage: matrix(QQ, 2, 2).right_kernel_matrix(algorithm='pluq')                # needs sage.libs.m4ri
            Traceback (most recent call last):
            ...
            ValueError: 'pluq' matrix kernel algorithm only available over
            integers mod 2, not over Rational Field

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
            ValueError: LLL-reduced basis only available over the integers,
            not over Rational Field

        Finally, error messages for the 'proof' keyword.  ::

            sage: matrix(ZZ, 2, 2).right_kernel_matrix(proof='junk')
            Traceback (most recent call last):
            ...
            ValueError: 'proof' must be one of True, False or None, not junk
            sage: matrix(QQ, 2, 2).right_kernel_matrix(proof=True)
            Traceback (most recent call last):
            ...
            ValueError: 'proof' flag only valid for matrices over the integers
        """
        R = self.base_ring()

        # First: massage keywords
        # Determine algorithm to use for computation of kernel matrix
        algorithm = kwds.pop('algorithm', None)
        if algorithm is None:
            algorithm = 'default'
        elif algorithm not in ['default', 'generic', 'flint', 'pari', 'padic', 'pluq']:
            raise ValueError("matrix kernel algorithm '%s' not recognized" % algorithm)
        elif algorithm == 'padic' and not isinstance(R, (IntegerRing_class,
                                                         RationalField)):
            raise ValueError("'padic' matrix kernel algorithm only available over the rationals and the integers, not over %s" % R)
        elif algorithm == 'flint' and not isinstance(R, (IntegerRing_class,
                                                         RationalField)):
            raise ValueError("'flint' matrix kernel algorithm only available over the rationals and the integers, not over %s" % R)
        elif algorithm == 'pari' and not (isinstance(R, (IntegerRing_class, NumberField)) and not isinstance(R, RationalField)):
            raise ValueError("'pari' matrix kernel algorithm only available over non-trivial number fields and the integers, not over %s" % R)
        elif algorithm == 'generic' and R not in _Fields:
            raise ValueError("'generic' matrix kernel algorithm only available over a field, not over %s" % R)
        elif algorithm == 'pluq' and not isinstance(self, sage.matrix.matrix_mod2_dense.Matrix_mod2_dense):
            raise ValueError("'pluq' matrix kernel algorithm only available over integers mod 2, not over %s" % R)

        # Determine the basis format of independent spanning set to return
        basis = kwds.pop('basis', None)
        if basis is None:
            basis = 'default'
        elif basis not in ['default', 'computed', 'echelon', 'pivot', 'LLL']:
            raise ValueError("matrix kernel basis format '%s' not recognized" % basis)
        elif basis == 'pivot' and R not in _Fields:
            raise ValueError('pivot basis only available over a field, not over %s' % R)
        elif basis == 'LLL' and not isinstance(R, IntegerRing_class):
            raise ValueError('LLL-reduced basis only available over the integers, not over %s' % R)
        if basis == 'default':
            basis = 'echelon' if R in _Fields else 'computed'

        # Determine proof keyword for integer matrices
        proof = kwds.pop('proof', None)
        if proof not in [None, True, False]:
            raise ValueError("'proof' must be one of True, False or None, not %s" % proof)
        if not (proof is None or isinstance(R, IntegerRing_class)):
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
            return self.new_matrix(nrows=0, ncols=self._ncols)
        if self._nrows == 0 and R.is_integral_domain():
            return self.matrix_space(self._ncols, self._ncols).identity_matrix()

        # Third: generic first, if requested explicitly
        #   then try specialized class methods, and finally
        #   delegate to ad-hoc methods in greater generality
        M = None
        format = ''

        if algorithm == 'generic':
            format, M = self._right_kernel_matrix_over_field()

        if M is None:
            try:
                format, M = self._right_kernel_matrix(algorithm=algorithm, proof=proof)
            except AttributeError:
                pass

        if M is None and isinstance(R, NumberField):
            format, M = self._right_kernel_matrix_over_number_field()

        if M is None and R in _Fields:
            format, M = self._right_kernel_matrix_over_field()

        if M is None and R.is_integral_domain():
            format, M = self._right_kernel_matrix_over_domain()

        if M is None and isinstance(R, sage.rings.abc.IntegerModRing):
            format, M = self._right_kernel_matrix_over_integer_mod_ring()

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
                M.echelonize()
                return M
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

    def left_kernel_matrix(self, *args, **kwds):
        r"""
        Return a matrix whose rows form a basis for the left kernel
        of ``self``.

        This method is a thin wrapper around :meth:`right_kernel_matrix`.
        For supported parameters and input/output formats, see there.

        EXAMPLES::

            sage: M = matrix([[1,2],[3,4],[5,6]])
            sage: K = M.left_kernel_matrix(); K
            [ 1 -2  1]
            sage: K * M
            [0 0]
        """
        return self.transpose().right_kernel_matrix(*args, **kwds)

    def right_kernel(self, *args, **kwds):
        r"""
        Return the right kernel of this matrix, as a vector space or
        free module. This is the set of vectors ``x`` such that ``self*x = 0``.

        .. NOTE::

            For the left kernel, use :meth:`left_kernel`.  The method
            :meth:`kernel` is exactly equal to :meth:`left_kernel`.

            For inexact rings use :meth:`right_kernel_matrix` with
            ``basis='computed'`` to avoid echelonizing.

        INPUT:

        - ``algorithm`` -- (default: ``'default'``) a keyword that selects the
          algorithm employed.  Allowable values are:

          - ``'default'`` -- allows the algorithm to be chosen automatically
          - ``'generic'`` -- naive algorithm usable for matrices over any field
          - ``'flint'`` -- FLINT library code for matrices over the rationals
            or the integers
          - ``'pari'`` -- PARI library code for matrices over number fields
            or the integers
          - ``'padic'`` -- padic algorithm from IML library for matrices
            over the rationals and integers
          - ``'pluq'`` -- PLUQ matrix factorization for matrices mod 2

        - ``basis`` -- (default: ``'echelon'``) a keyword that describes the
          format of the basis used to construct the right kernel.
          Allowable values are:

          - 'echelon': the basis matrix is returned in echelon form
          - 'pivot' : each basis vector is computed from the reduced
            row-echelon form of ``self`` by placing a single one in a
            non-pivot column and zeros in the remaining non-pivot columns.
            Only available for matrices over fields.
          - 'LLL': an LLL-reduced basis.  Only available for matrices
            over the integers.

        OUTPUT:

        A vector space or free module whose degree equals the number
        of columns in ``self`` and which contains all the vectors ``x``
        such that ``self*x = 0``.

        If ``self`` has 0 columns, the kernel has dimension 0, while if
        ``self`` has 0 rows the kernel is the entire ambient vector space.

        The result is cached.  Requesting the right kernel a second time,
        but with a different basis format, will return the cached result
        with the format from the first computation.

        .. NOTE::

           For more detailed documentation on the selection of algorithms
           used and a more flexible method for computing a basis matrix
           for a right kernel (rather than computing a vector space), see
           :meth:`right_kernel_matrix`, which powers the computations for
           this method.

        EXAMPLES::

            sage: A = matrix(QQ, [[0, 0, 1, 2, 2, -5, 3],
            ....:                 [-1, 5, 2, 2, 1, -7, 5],
            ....:                 [0, 0, -2, -3, -3, 8, -5],
            ....:                 [-1, 5, 0, -1, -2, 1, 0]])
            sage: K = A.right_kernel(); K
            Vector space of degree 7 and dimension 4 over Rational Field
            Basis matrix:
            [ 1  0  0  0 -1 -1 -1]
            [ 0  1  0  0  5  5  5]
            [ 0  0  1  0 -1 -2 -3]
            [ 0  0  0  1  0  1  1]
            sage: A * K.basis_matrix().transpose() == zero_matrix(QQ, 4, 4)
            True

        The default is basis vectors that form a matrix in echelon form.
        A "pivot basis" instead has a basis matrix where the columns of
        an identity matrix are in the locations of the non-pivot columns
        of the original matrix. This alternate format is available whenever
        the base ring is a field. ::

            sage: A = matrix(QQ, [[0, 0, 1, 2, 2, -5, 3],
            ....:                 [-1, 5, 2, 2, 1, -7, 5],
            ....:                 [0, 0, -2, -3, -3, 8, -5],
            ....:                 [-1, 5, 0, -1, -2, 1, 0]])
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
            sage: A * K.basis_matrix().transpose() == zero_matrix(QQ, 4, 4)
            True

        Matrices may have any field as a base ring.  Number fields are
        computed by PARI library code, matrices over `GF(2)` are computed
        by the M4RI library, and matrices over the rationals are computed by
        the IML library.  For any of these specialized cases, general-purpose
        code can be called instead with the keyword setting
        ``algorithm='generic'``.

        Over an arbitrary field, with two basis formats.  Same vector space,
        different bases.  ::

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = FiniteField(5^2)
            sage: A = matrix(F, 3, 4, [[  1,   a,     1+a,  a^3+a^5],
            ....:                      [  a, a^4,   a+a^4,  a^4+a^8],
            ....:                      [a^2, a^6, a^2+a^6, a^5+a^10]])
            sage: K = A.right_kernel(); K
            Vector space of degree 4 and dimension 2
             over Finite Field in a of size 5^2
            Basis matrix:
            [      1       0 3*a + 4 2*a + 2]
            [      0       1     2*a 3*a + 3]
            sage: A * K.basis_matrix().transpose() == zero_matrix(F, 3, 2)
            True

        In the following test, we have to force usage of
        :class:`~sage.matrix.matrix_generic_dense.Matrix_generic_dense`,
        since the option ``basis = 'pivot'`` would simply yield the same
        result as the previous test, if the optional meataxe package is
        installed. ::

            sage: from sage.matrix.matrix_generic_dense import Matrix_generic_dense
            sage: B = Matrix_generic_dense(A.parent(), A.list(), False, False)
            sage: P = B.right_kernel(basis='pivot'); P                                  # needs sage.rings.finite_rings
            Vector space of degree 4 and dimension 2
             over Finite Field in a of size 5^2
            User basis matrix:
            [      4       4       1       0]
            [  a + 2 3*a + 3       0       1]

        If the optional meataxe package is installed, we again have to make sure
        to work with a copy of B that has the same type as ``P.basis_matrix()``::

            sage: (B.parent()(B.list()) * P.basis_matrix().transpose()                  # needs sage.rings.finite_rings
            ....:   == zero_matrix(F, 3, 2))
            True
            sage: K == P                                                                # needs sage.rings.finite_rings
            True

        Over number fields, PARI is used by default, but general-purpose code
        can be requested.  Same vector space, same bases, different code.::

            sage: # needs sage.rings.number_field
            sage: Q = QuadraticField(-7)
            sage: a = Q.gen(0)
            sage: A = matrix(Q, [[    2, 5 - a,   15 - a, 16 + 4*a],
            ....:                [2 + a,     a, -7 + 5*a, -3 + 3*a]])
            sage: K = A.right_kernel(algorithm='default'); K
            Vector space of degree 4 and dimension 2
             over Number Field in a with defining polynomial x^2 + 7 with a = 2.645751311064591?*I
            Basis matrix:
            [                1                 0     7/88*a + 3/88 -3/176*a - 39/176]
            [                0                 1   -1/88*a - 13/88  13/176*a - 7/176]
            sage: A * K.basis_matrix().transpose() == zero_matrix(Q, 2, 2)
            True
            sage: B = copy(A)
            sage: G = A.right_kernel(algorithm='generic'); G
            Vector space of degree 4 and dimension 2
             over Number Field in a with defining polynomial x^2 + 7 with a = 2.645751311064591?*I
            Basis matrix:
            [                1                 0     7/88*a + 3/88 -3/176*a - 39/176]
            [                0                 1   -1/88*a - 13/88  13/176*a - 7/176]
            sage: B * G.basis_matrix().transpose() == zero_matrix(Q, 2, 2)
            True
            sage: K == G
            True

        For matrices over the integers, several options are possible.
        The basis can be an LLL-reduced basis or an echelon basis.
        The pivot basis isnot available.  A heuristic will decide whether
        to use a `p`-adic algorithm from the IML library or an algorithm
        from the PARI library.  Note how specifying the algorithm can
        mildly influence the LLL basis. ::

            sage: A = matrix(ZZ, [[0, -1, -1, 2, 9, 4, -4],
            ....:                 [-1, 1, 0, -2, -7, -1, 6],
            ....:                 [2, 0, 1, 0, 1, -5, -2],
            ....:                 [-1, -1, -1, 3, 10, 10, -9],
            ....:                 [-1, 2, 0, -3, -7, 1, 6]])
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
            ....:                [y^3, y^2, 2*y^3]])
            sage: K = A.right_kernel(algorithm='default', basis='echelon'); K
            Free module of degree 3 and rank 1
             over Univariate Polynomial Ring in y over Rational Field
             Echelon basis matrix:
             [-1 -y  1]
            sage: A * K.basis_matrix().transpose() == zero_matrix(ZZ, 2, 1)
            True

        It is possible to compute a kernel for a matrix over an integral
        domain which is not a PID, but usually this will fail.  ::

            sage: D.<x> = ZZ[]
            sage: A = matrix(D, 2, 2, [[x^2 - x, -x + 5],
            ....:                      [x^2 - 8, -x + 2]])
            sage: A.right_kernel()
            Traceback (most recent call last):
            ...
            ArithmeticError: Ideal Ideal (x^2 - x, x^2 - 8) of Univariate
            Polynomial Ring in x over Integer Ring not principal

        Matrices over non-commutative rings are not a good idea either.
        These are the "usual" quaternions.  ::

            sage: Q.<i,j,k> = QuaternionAlgebra(-1,-1)
            sage: A = matrix(Q, 2, [i,j,-1,k])
            sage: A.right_kernel()
            Traceback (most recent call last):
            ...
            NotImplementedError: Cannot compute a matrix kernel over
            Quaternion Algebra (-1, -1) with base ring Rational Field

        Sparse matrices, over the rationals and the integers,
        use the same routines as the dense versions. ::

            sage: A = matrix(ZZ, [[0, -1, 1, 1, 2],
            ....:                 [1, -2, 0, 1, 3],
            ....:                 [-1, 2, 0, -1, -3]],
            ....:            sparse=True)
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

            sage: A = matrix(QQ, 3, 3, range(9))
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
        if K is not None:
            return K

        R = self.base_ring()
        tm = verbose(lazy_string("computing a right kernel for %sx%s matrix over %s", self.nrows(), self.ncols(), R), level=2)

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

        verbose(lazy_string("done computing a right kernel for %sx%s matrix over %s", self.nrows(), self.ncols(), R), level=2, t=tm)
        self.cache('right_kernel', K)
        return K

    def left_kernel(self, *args, **kwds):
        r"""
        Return the left kernel of this matrix, as a vector space or free module.

        This is the set of vectors ``x`` such that ``x*self = 0``.

        .. NOTE::

            For the right kernel, use :meth:`right_kernel`.  The method
            :meth:`kernel` is exactly equal to :meth:`left_kernel`.

            For inexact rings use :meth:`right_kernel_matrix` with
            ``basis='computed'`` (on the transpose of the matrix) to avoid echelonizing.

        INPUT:

        - ``algorithm`` -- (default: ``'default'``) a keyword that selects the
          algorithm employed.  Allowable values are:

          - ``'default'`` -- allows the algorithm to be chosen automatically
          - ``'generic'`` -- naive algorithm usable for matrices over any field
          - ``'flint'`` -- FLINT library code for matrices over the rationals
            or the integers
          - ``'pari'`` -- PARI library code for matrices over number fields
            or the integers
          - ``'padic'`` -- padic algorithm from IML library for matrices
            over the rationals and integers
          - ``'pluq'`` -- PLUQ matrix factorization for matrices mod 2

        - ``basis`` -- (default: ``'echelon'``) a keyword that describes
          the format of the basis used to construct the left kernel.
          Allowable values are:

          - 'echelon': the basis matrix is returned in echelon form
          - 'pivot' : each basis vector is computed from the reduced
            row-echelon form of ``self`` by placing a single one in a
            non-pivot column and zeros in the remaining non-pivot columns.
            Only available for matrices over fields.
          - 'LLL': an LLL-reduced basis.  Only available for matrices
            over the integers.

        OUTPUT:

        A vector space or free module whose degree equals the number
        of rows in ``self`` and which contains all the vectors ``x`` such
        that ``x*self = 0``.

        If ``self`` has 0 rows, the kernel has dimension 0, while if ``self``
        has 0 columns the kernel is the entire ambient vector space.

        The result is cached.  Requesting the left kernel a second time,
        but with a different basis format, will return the cached result
        with the format from the first computation.

        .. NOTE::

           For much more detailed documentation of the various options see
           :meth:`right_kernel`, since this method just computes
           the right kernel of the transpose of ``self``.

        EXAMPLES:

        Over the rationals with a basis matrix in echelon form. ::

            sage: A = matrix(QQ, [[1, 2, 4, -7, 4],
            ....:                 [1, 1, 0, 2, -1],
            ....:                 [1, 0, 3, -3, 1],
            ....:                 [0, -1, -1, 3, -2],
            ....:                 [0, 0, -1, 2, -1]])
            sage: A.left_kernel()
            Vector space of degree 5 and dimension 2 over Rational Field
            Basis matrix:
            [ 1  0 -1  2 -1]
            [ 0  1 -1  1 -4]

        Over a finite field, with a basis matrix in "pivot" format. ::

            sage: A = matrix(FiniteField(7), [[5, 0, 5, 2, 4],
            ....:                             [1, 3, 2, 3, 6],
            ....:                             [1, 1, 6, 5, 3],
            ....:                             [2, 5, 6, 0, 0]])
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
        if K is not None:
            return K

        tm = verbose("computing left kernel for %sx%s matrix" % (self.nrows(), self.ncols()), level=2)
        K = self.transpose().right_kernel(*args, **kwds)
        self.cache('left_kernel', K)
        verbose("done computing left kernel for %sx%s matrix" % (self.nrows(), self.ncols()), level=2, t=tm)
        return K

    kernel = left_kernel

    def kernel_on(self, V, poly=None, check=True):
        """
        Return the kernel of ``self`` restricted to the invariant subspace `V`.
        The result is a vector subspace of `V`, which is also a subspace
        of the ambient space.

        INPUT:

        - ``V`` -- vector subspace

        - ``check`` -- boolean (default: ``True``); whether to check that
          `V` is invariant under the action of ``self``

        - ``poly`` -- (default: ``None``) if not ``None``, compute instead
          the kernel of ``poly(self)`` on `V`

        OUTPUT: a subspace

        .. warning::

           This function does *not* check that `V` is in fact
           invariant under ``self`` if check is ``False``.  With check
           ``False`` this function is much faster.

        EXAMPLES::

            sage: t = matrix(QQ, 4, [39, -10, 0, -12, 0, 2, 0, -1, 0, 1, -2, 0, 0, 2, 0, -2]); t
            [ 39 -10   0 -12]
            [  0   2   0  -1]
            [  0   1  -2   0]
            [  0   2   0  -2]
            sage: t.fcp()                                                               # needs sage.libs.pari
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

        Test that :issue:`9425` is fixed.

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
            [(0, 1, 0)]
            sage: W.is_submodule(V)
            True
        """
        A = self.restrict(V, check=check)
        if poly is not None:
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
        the integers for the kernel of ``self``.

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

            sage: A = MatrixSpace(QQ, 2)([1/2, 0, 0, 0])
            sage: A.integer_kernel()
            Free module of degree 2 and rank 1 over Integer Ring
            Echelon basis matrix:
            [0 1]

        An example over a bigger ring::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: L.<w> = NumberField(x^2 - x + 2)
            sage: OL = L.ring_of_integers()
            sage: A = matrix(L, 2, [1, w/2])
            sage: A.integer_kernel(OL)
            Free module of degree 2 and rank 1 over
             Maximal Order generated by w in Number Field in w
             with defining polynomial x^2 - x + 2
             Echelon basis matrix:
             [    -1 -w + 1]
        """
        try:
            A, _ = self._clear_denom()
            return A.kernel()
        except AttributeError:
            from sage.matrix.matrix_space import MatrixSpace
            d = self.denominator()
            A = self * d
            M = MatrixSpace(ring, self.nrows(), self.ncols())(A)
            return M.kernel()

    def image(self):
        """
        Return the image of the homomorphism on rows defined by right
        multiplication by this matrix: that is, the row-space.

        EXAMPLES::

            sage: MS1 = MatrixSpace(ZZ, 4)
            sage: MS2 = MatrixSpace(QQ, 6)
            sage: A = MS1.matrix([3,4,5,6, 7,3,8,10, 14,5,6,7, 2,2,10,9])
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
            sage: image(B) == B.transpose().column_module()
            True

        .. SEEALSO::

            :meth:`row_module`, :meth:`column_module`
        """
        return self.row_module()

    def row_module(self, base_ring=None):
        """
        Return the free module over the base ring spanned by the rows of
        ``self``.

        EXAMPLES::

            sage: A = MatrixSpace(IntegerRing(), 2)([1,2,3,4])
            sage: A.row_module()
            Free module of degree 2 and rank 2 over Integer Ring
            Echelon basis matrix:
            [1 0]
            [0 2]
        """
        M = self.row_ambient_module(base_ring=base_ring)
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

            sage: t = matrix(QQ, 3, 3, range(9)); t
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: t.row_space()
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [ 1  0 -1]
            [ 0  1  2]

        ::

            sage: m = Matrix(Integers(5),2,2,[2,2,2,2])
            sage: m.row_space()
            Vector space of degree 2 and dimension 1 over Ring of integers modulo 5
            Basis matrix:
            [1 1]
        """
        return self.row_module(base_ring=base_ring)

    def column_module(self):
        """
        Return the free module over the base ring spanned by the columns of
        this matrix.

        EXAMPLES::

            sage: t = matrix(QQ, 3, 3, range(9)); t
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

            sage: M = MatrixSpace(QQ, 3, 3)
            sage: A = M([1,9,-7, 4/5,4,3, 6,4,3])
            sage: A.column_space()
            Vector space of degree 3 and dimension 3 over Rational Field
            Basis matrix:
            [1 0 0]
            [0 1 0]
            [0 0 1]

            sage: # needs sage.rings.real_mpfr sage.symbolic
            sage: W = MatrixSpace(CC, 2, 2)
            sage: B = W([1, 2 + 3*I, 4 + 5*I, 9]); B
            [                     1.00000000000000 2.00000000000000 + 3.00000000000000*I]
            [4.00000000000000 + 5.00000000000000*I                      9.00000000000000]
            sage: B.column_space()
            Vector space of degree 2 and dimension 2
             over Complex Field with 53 bits of precision
             Basis matrix:
             [ 1.00000000000000 0.000000000000000]
             [0.000000000000000  1.00000000000000]
        """
        return self.column_module()

    def decomposition(self, algorithm='spin',
                      is_diagonalizable=False, dual=False):
        """
        Return the decomposition of the free module on which this matrix A
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

        - ``self`` -- a matrix

        - ``algorithm`` -- string (default: ``'spin'``); ``'spin'``: involves
          iterating the action of ``self`` on a vector. ``'kernel'``: naively
          just compute `ker(f_i(A))` for each factor `f_i`.

        - ``dual`` -- boolean (default: ``False``); if True, also
          returns the corresponding decomposition of V under the action of
          the transpose of A. The factors are guaranteed to correspond.

        - ``is_diagonalizable`` -- if the matrix is known to
          be diagonalizable, set this to True, which might speed up the
          algorithm in some cases.

        .. NOTE::

           If the base ring is not a field, the kernel algorithm is
           used.

        OUTPUT:

        - ``Sequence`` -- list of pairs (V,t), where V is a vector
          spaces and t is a boolean, and t is ``True`` exactly when the
          charpoly of ``self`` on V is irreducible.

        - (optional) list -- list of pairs (W,t), where W is a vector
          space and t is a boolean, and t is ``True`` exactly when the
          charpoly of the transpose of ``self`` on W is irreducible.

        EXAMPLES::

            sage: A = matrix(ZZ, 4, [3,4,5,6, 7,3,8,10, 14,5,6,7, 2,2,10,9])
            sage: B = matrix(QQ, 6, 6, range(36))
            sage: B*11
            [  0  11  22  33  44  55]
            [ 66  77  88  99 110 121]
            [132 143 154 165 176 187]
            [198 209 220 231 242 253]
            [264 275 286 297 308 319]
            [330 341 352 363 374 385]
            sage: A.decomposition()                                                     # needs sage.libs.pari
            [(Ambient free module of rank 4 over the principal ideal domain Integer Ring,
              True)]
            sage: B.decomposition()                                                     # needs sage.libs.pari
            [(Vector space of degree 6 and dimension 2 over Rational Field
              Basis matrix:
              [ 1  0 -1 -2 -3 -4]
              [ 0  1  2  3  4  5],
              True),
             (Vector space of degree 6 and dimension 4 over Rational Field
              Basis matrix:
              [ 1  0  0  0 -5  4]
              [ 0  1  0  0 -4  3]
              [ 0  0  1  0 -3  2]
              [ 0  0  0  1 -2  1],
              False)]
        """
        if algorithm == 'kernel' or self.base_ring() not in _Fields:
            return self._decomposition_using_kernels(is_diagonalizable=is_diagonalizable, dual=dual)
        elif algorithm == 'spin':
            X = self._decomposition_spin_generic(is_diagonalizable=is_diagonalizable)
            if dual:
                Y = self.transpose()._decomposition_spin_generic(is_diagonalizable=is_diagonalizable)
                return X, Y
            return X
        else:
            raise ValueError("no algorithm '%s'"%algorithm)

    def _decomposition_spin_generic(self, is_diagonalizable=False):
        r"""
        Compute the decomposition of this matrix using the spin algorithm.

        INPUT:

        - ``self`` -- a matrix with field entries

        OUTPUT: list of reduced row echelon form basis
        """
        if not self.is_square():
            raise ValueError("self must be a square matrix")

        if self.base_ring() not in _Fields:
            raise TypeError("self must be over a field.")

        if self.nrows() == 0:
            return decomp_seq([])

        f = self.charpoly('x')
        E = decomp_seq([])

        t = verbose('factoring the characteristic polynomial', level=2, caller_name='generic spin decomp')
        F = f.factor()
        verbose('done factoring', t=t, level=2, caller_name='generic spin decomp')

        if len(F) == 1:
            V = self.base_ring()**self.nrows()
            return decomp_seq([(V, F[0][1]==1)])

        V = self.base_ring()**self.nrows()
        v = V.random_element()
        num_iterates = max([0] + [f.degree() - g.degree() for g, _ in F if g.degree() > 1]) + 1

        S = []

        F.sort()
        for i in range(len(F)):
            g, m = F[i]

            if g.degree() == 1:
                # Just use kernel -- much easier.
                B = self.__copy__()
                for k from 0 <= k < self.nrows():
                    B[k, k] += g[0]
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
                    t = verbose('compute element of kernel of g(A), for g of degree %s'%g.degree(), level=2,
                                caller_name='generic spin decomp')
                    w = S[j].linear_combination_of_rows(h.list())
                    t = verbose('done computing element of kernel of g(A)', t=t, level=2, caller_name='generic spin decomp')

                    # Get the rest of the kernel.
                    t = verbose('fill out rest of kernel', level=2, caller_name='generic spin decomp')
                    if W is None:
                        W = self.iterates(w, g.degree())
                    else:
                        W = W.stack(self.iterates(w, g.degree()))
                    t = verbose('finished filling out more of kernel', level=2, t=t, caller_name='generic spin decomp')

                if W.rank() == m * g.degree():
                    t = verbose('now computing row space', level=2, caller_name='generic spin decomp')
                    W.echelonize()
                    E.append((W.row_space(), m==1))
                    verbose('computed row space', level=2, t=t, caller_name='generic spin decomp')
                    break
                else:
                    verbose('we have not yet generated all the kernel (rank so far=%s, target rank=%s)'%(
                        W.rank(), m*g.degree()), level=2, caller_name='generic spin decomp')
                    tries += 1
                    if tries > 1000*m:  # avoid an insanely long infinite loop
                        raise RuntimeError("likely bug in decomposition")
                # end if
            # end while
        # end for
        return E

    def _decomposition_using_kernels(self, is_diagonalizable=False, dual=False):
        if not self.is_square():
            raise ValueError("self must be a square matrix")

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
            V = self.column_ambient_module()
            m = F[0][1]
            if dual:
                return decomp_seq([(V, m==1)]), decomp_seq([(V, m==1)])
            else:
                return decomp_seq([(V, m==1)])
        F.sort()
        for g, m in f.factor():
            t = verbose('decomposition -- Computing g(self) for an irreducible factor g of degree %s'%g.degree(), level=2)
            if is_diagonalizable:
                B = g(self)
            else:
                B = g(self)
                t2 = verbose('decomposition -- raising g(self) to the power %s'%m, level=2)
                B = B ** m
                verbose('done powering', level=2, t=t2)
            t = verbose('decomposition -- done computing g(self)', level=2, t=t)
            E.append((B.kernel(), m==1))
            t = verbose('decomposition -- time to compute kernel', level=2, t=t)
            if dual:
                Edual.append((B.transpose().kernel(), m==1))
                verbose('decomposition -- time to compute dual kernel', level=2, t=t)
        if dual:
            return E, Edual
        return E

    def decomposition_of_subspace(self, M, check_restrict=True, **kwds):
        """
        Suppose the right action of ``self`` on M leaves M invariant. Return
        the decomposition of M as a list of pairs (W, is_irred) where
        is_irred is ``True`` if the charpoly of ``self`` acting on the factor W is
        irreducible.

        Additional inputs besides M are passed onto the decomposition
        command.

        INPUT:

            - ``M`` -- a subspace of the free module ``self`` acts on
            - ``check_restrict`` -- boolean (default: ``True``); call restrict
              with or without check
            - ``kwds`` -- keywords that will be forwarded to :meth:`~.decomposition`

        EXAMPLES::

            sage: # needs sage.libs.pari
            sage: t = matrix(QQ, 3, [3, 0, -2, 0, -2, 0, 0, 0, 0]); t
            [ 3  0 -2]
            [ 0 -2  0]
            [ 0  0  0]
            sage: t.fcp('X')                # factored charpoly
            (X - 3) * X * (X + 2)
            sage: v = kernel(t*(t+2)); v    # an invariant subspace
            Vector space of degree 3 and dimension 2 over Rational Field
            Basis matrix:
            [0 1 0]
            [0 0 1]
            sage: D = t.decomposition_of_subspace(v); D
            [(Vector space of degree 3 and dimension 1 over Rational Field
              Basis matrix:
              [0 0 1],
              True),
             (Vector space of degree 3 and dimension 1 over Rational Field
              Basis matrix:
              [0 1 0],
              True)]
            sage: t.restrict(D[0][0])
            [0]
            sage: t.restrict(D[1][0])
            [-2]

        We do a decomposition over ZZ::

            sage: a = matrix(ZZ, 6, [0, 0, -2, 0, 2, 0,
            ....:                    2, -4, -2, 0, 2, 0,
            ....:                    0, 0, -2, -2, 0, 0,
            ....:                    2, 0, -2, -4, 2, -2,
            ....:                    0, 2, 0, -2, -2, 0,
            ....:                    0, 2, 0, -2, 0, 0])
            sage: a.decomposition_of_subspace(ZZ^6)                                     # needs sage.libs.pari
            [(Free module of degree 6 and rank 2 over Integer Ring
              Echelon basis matrix:
              [ 1  0  1 -1  1 -1]
              [ 0  1  0 -1  2 -1],
              False),
             (Free module of degree 6 and rank 4 over Integer Ring
              Echelon basis matrix:
              [ 1  0 -1  0  1  0]
              [ 0  1  0  0  0  0]
              [ 0  0  0  1  0  0]
              [ 0  0  0  0  0  1],
              False)]

        TESTS::

            sage: t = matrix(QQ, 3, [3, 0, -2, 0, -2, 0, 0, 0, 0])
            sage: t.decomposition_of_subspace(v, check_restrict=False) == t.decomposition_of_subspace(v)                # needs sage.libs.pari
            True
        """
        if not isinstance(M, sage.modules.free_module.FreeModule_generic):
            raise TypeError("M must be a free module.")
        if not self.is_square():
            raise ArithmeticError("self must be a square matrix")
        if M.base_ring() != self.base_ring():
            raise ArithmeticError("base rings must be the same, but self is over %s and module is over %s"%(
                self.base_ring(), M.base_ring()))
        if M.degree() != self.ncols():
            raise ArithmeticError("M must be a subspace of an %s-dimensional space" % self.ncols())

        time = verbose(level=2, t=0)

        # 1. Restrict
        B = self.restrict(M, check=check_restrict)
        time0 = verbose("decompose restriction -- ", level=2, t=time)

        # 2. Decompose restriction
        D = B.decomposition(**kwds)

        sum_dim = sum([A.dimension() for A, _ in D])
        assert sum_dim == M.dimension(), \
               "bug in decomposition; " + \
               "the sum of the dimensions (=%s) of the factors must equal the dimension (%s) of the acted on space:\nFactors found: %s\nSpace: %s"%(sum_dim, M.dimension(), D, M)

        # 3. Lift decomposition to subspaces of ambient vector space.
        # Each basis vector for an element of D defines a linear
        # combination of the basis of W, and these linear combinations
        # define the corresponding subspaces of the ambient space M.

        verbose("decomposition -- ", level=2, t=time0)
        C = M.basis_matrix()

        D = [((W.basis_matrix() * C).row_module(self.base_ring()), is_irred) for W, is_irred in D]
        D = decomp_seq(D)

        verbose(level=2, t=time)
        return D

    def restrict(self, V, check=True):
        """
        Return the matrix that defines the action of ``self`` on the chosen
        basis for the invariant subspace V. If V is an ambient, returns
        ``self`` (not a copy of self).

        INPUT:

        - ``V`` -- vector subspace

        - ``check`` -- boolean (default: ``True``); if ``False`` may
          not check that V is invariant (hence can be faster)

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
            raise TypeError("V must be a free module")
        if V.degree() != self.nrows():
            raise IndexError("degree of V (=%s) must equal number of rows of self (=%s)" % (V.degree(), self.nrows()))
        if V.rank() == 0 or V.degree() == 0:
            return self.new_matrix(nrows=0, ncols=0)

        if not check and V.base_ring() in _Fields and not V.has_user_basis():
            B = V.echelonized_basis_matrix()
            P = B.pivots()
            return B*self.matrix_from_columns(P)
        else:
            n = V.rank()
            try:
                # todo optimize so only involves matrix multiplies ?
                C = [V.coordinate_vector(b*self) for b in V.basis()]
            except ArithmeticError:
                raise ArithmeticError("subspace is not invariant under matrix")
            return self.new_matrix(n, n, C, sparse=False)

    def restrict_domain(self, V):
        """
        Compute the matrix relative to the basis for V on the domain
        obtained by restricting ``self`` to V, but not changing the codomain of
        the matrix. This is the matrix whose rows are the images of the
        basis for V.

        INPUT:

        - ``V`` -- vector space (subspace of ambient space on
          which ``self`` acts)


        .. SEEALSO::

           :meth:`restrict`

        EXAMPLES::

            sage: V = QQ^3
            sage: A = matrix(QQ, 3, [1,2,0, 3,4,0, 0,0,0])
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
        Suppose that ``self`` defines a linear map from some domain to a
        codomain that contains `V` and that the image of ``self`` is
        contained in `V`. This function returns a new matrix
        `A` that represents this linear map but as a map to
        `V`, in the sense that if `x` is in the domain,
        then `xA` is the linear combination of the elements of the
        basis of `V` that equals v\*self.

        INPUT:

        - ``V`` -- vector space (space of degree
          ``self.ncols()``) that contains the image of ``self``

        .. SEEALSO::

           :meth:`restrict`, :meth:`restrict_domain`

        EXAMPLES::

            sage: A = matrix(QQ, 3, [1..9])
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
        Compute the largest integer n such that the list of vectors
        `S=[v, v*A, ..., v * A^n]` are linearly independent, and
        returns that list.

        INPUT:

        - ``self`` -- matrix

        - ``v`` -- vector

        OUTPUT: list of Vectors

        ALGORITHM: The current implementation just adds vectors to a vector
        space until the dimension doesn't grow. This could be optimized by
        directly using matrices and doing an efficient Echelon form. Also,
        when the base is Q, maybe we could simultaneously keep track of
        what is going on in the reduction modulo p, which might make things
        much faster.

        EXAMPLES::

            sage: t = matrix(QQ, 3, 3, range(9)); t
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
        if not isinstance(v, FreeModuleElement):
            raise TypeError("v must be a FreeModuleElement")
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
        Application of Wiedemann's algorithm to the `i`-th standard basis
        vector.

        INPUT:

        - ``i`` -- integer

        - ``t`` -- integer (default: 0); if t is nonzero, use
          only the first t linear recurrence relations

        IMPLEMENTATION: This is a toy implementation.

        EXAMPLES::

            sage: t = matrix(QQ, 3, 3, range(9)); t
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: t.wiedemann(0)
            x^2 - 12*x - 18
            sage: t.charpoly()                                                          # needs sage.libs.pari
            x^3 - 12*x^2 - 18*x
        """
        i = int(i)
        t = int(t)
        if self.nrows() != self.ncols():
            raise ArithmeticError("self must be a square matrix")
        n = self.nrows()
        v = sage.modules.free_module.VectorSpace(self.base_ring(), n).gen(i)
        tm = verbose('computing iterates...', level=2)
        cols = self.iterates(v, 2*n).columns()
        tm = verbose('computed iterates', level=2, t=tm)
        f = None
        # Compute the minimal polynomial of the linear recurrence
        # sequence corresponding to the 0-th entries of the iterates,
        # then the 1-th entries, etc.
        if t == 0:
            R = list(range(n))
        else:
            R = [t]
        for i in R:
            tm = verbose('applying berlekamp-massey', level=2)
            g = berlekamp_massey.berlekamp_massey(cols[i].list())
            verbose('berlekamp-massey done', level=2, t=tm)
            if f is None:
                f = g
            else:
                f = f.lcm(g)
            if f.degree() == n:
                break
        return f

    def _eigenspace_format(self, format):
        r"""
        Helper method to control output format for eigenspaces.

        INPUT:

        - ``format`` -- ``None``, ``'all'`` or ``'galois'``

        OUTPUT:

        Any format except ``None`` is just passed through.  When the
        format is ``None``, a choice is made about the style of the output.
        If there is an algebraically closed field that will contain the
        possible eigenvalues, then ``'all'`` of the eigenspaces are given.

        However if this is not the case, then only one eigenspace is output
        for each irreducible factor of the characteristic polynomial.

        EXAMPLES:

        Pass-through first::

            sage: A = matrix(QQ, 2, range(4))
            sage: A._eigenspace_format('all') == 'all'
            True
            sage: A._eigenspace_format('galois') == 'galois'
            True

        The algebraic closure of the rationals (the field ``QQbar`` of
        algebraic numbers) is implemented, as are algebraic closures
        of finite fields::

            sage: A = matrix(QQ, 2, range(4))
            sage: A._eigenspace_format(None) == 'all'                                   # needs sage.rings.number_field
            True
            sage: B = matrix(GF(13), 2, range(4))
            sage: B._eigenspace_format(None)
            'all'

        Subrings are promoted to fraction fields and then checked for the
        existence of algebraic closures.  ::

            sage: A = matrix(ZZ, 2, range(4))
            sage: A._eigenspace_format(None) == 'all'                                   # needs sage.rings.number_field
            True
        """
        if format not in [None, 'all', 'galois']:
            msg = "format keyword must be None, 'all' or 'galois', not {0}"
            raise ValueError(msg.format(format))

        # In the absence of a format keyword, we default to 'all' for
        # subrings of fields of which an algebraic closure is implemented.
        if format is None:
            try:
                self.base_ring().fraction_field().algebraic_closure()
                return 'all'
            except (NotImplementedError, AttributeError):
                return 'galois'
        else:
            return format

    def eigenspaces_left(self, format='all', var='a', algebraic_multiplicity=False):
        r"""
        Compute the left eigenspaces of a matrix.

        Note that :meth:`eigenspaces_left` and :meth:`left_eigenspaces`
        are identical methods.  Here "left" refers to the eigenvectors
        being placed to the left of the matrix.

        INPUT:

        - ``self`` -- a square matrix over an exact field.  For inexact
          matrices consult the numerical or symbolic matrix classes

        - ``format`` -- one of:

          - ``'all'`` -- attempts to create every eigenspace.  This will
            always be possible for matrices with rational entries
          - ``'galois'`` -- for each irreducible factor of the characteristic
            polynomial, a single eigenspace will be output for a
            single root/eigenvalue for the irreducible factor
          - ``None`` -- default; uses the ``'all'`` format if the base ring is contained
            in an algebraically closed field which is implemented.
            Otherwise, uses the ``'galois'`` format.

        - ``var`` -- string (default: ``'a'``); variable name used to
          represent elements of the root field of each
          irreducible factor of the characteristic polynomial.
          If ``var='a'``, then the root fields will be in terms of
          ``a0, a1, a2, ...``, where the numbering runs across all
          the irreducible factors of the characteristic polynomial,
          even for linear factors.

        - ``algebraic_multiplicity`` -- boolean (default: ``False``);
          whether to include the algebraic multiplicity of each eigenvalue
          in the output.  See the discussion below.

        OUTPUT:

        If ``algebraic_multiplicity=False``, return a list of pairs `(e, V)`
        where `e` is an eigenvalue of the matrix, and `V` is the corresponding
        left eigenspace.  For Galois conjugates of eigenvalues, there
        may be just one representative eigenspace, depending on the
        ``format`` keyword.

        If ``algebraic_multiplicity=True``, return a list of triples `(e, V, n)`
        where `e` and `V` are as above and `n` is the algebraic multiplicity of
        the eigenvalue.

        .. warning::

           Uses a somewhat naive algorithm (simply factors the
           characteristic polynomial and computes kernels directly
           over the extension field).

        EXAMPLES:

        We compute the left eigenspaces of a `3\times 3`
        rational matrix. First, we request ``'all'`` of the eigenvalues,
        so the results are in the field of algebraic numbers, ``QQbar``.
        Then we request just one eigenspace per irreducible factor of
        the characteristic polynomial with ``format='galois'``.  ::

            sage: A = matrix(QQ, 3, 3, range(9)); A
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: es = A.eigenspaces_left(format='all'); es                             # needs sage.rings.number_field
            [(0,
              Vector space of degree 3 and dimension 1 over Rational Field
              User basis matrix:
              [ 1 -2  1]),
             (-1.348469228349535?,
              Vector space of degree 3 and dimension 1 over Algebraic Field
              User basis matrix:
              [                   1  0.3101020514433644? -0.3797958971132713?]),
             (13.34846922834954?,
              Vector space of degree 3 and dimension 1 over Algebraic Field
              User basis matrix:
              [                 1 1.289897948556636? 1.579795897113272?])]

            sage: # needs sage.rings.number_field
            sage: es = A.eigenspaces_left(format='galois'); es
            [(0,
              Vector space of degree 3 and dimension 1 over Rational Field
              User basis matrix:
              [ 1 -2  1]),
             (a1,
              Vector space of degree 3 and dimension 1 over Number Field in a1 with defining polynomial x^2 - 12*x - 18
              User basis matrix:
              [            1 1/15*a1 + 2/5 2/15*a1 - 1/5])]
            sage: es = A.eigenspaces_left(format='galois',
            ....:                         algebraic_multiplicity=True); es
            [(0,
              Vector space of degree 3 and dimension 1 over Rational Field
              User basis matrix:
              [ 1 -2  1],
              1),
             (a1,
              Vector space of degree 3 and dimension 1 over Number Field in a1 with defining polynomial x^2 - 12*x - 18
              User basis matrix:
              [            1 1/15*a1 + 2/5 2/15*a1 - 1/5],
              1)]
            sage: e, v, n = es[0]; v = v.basis()[0]
            sage: delta = e*v - v*A
            sage: abs(abs(delta)) < 1e-10
            True

        The same computation, but with implicit base change to a field.  ::

            sage: A = matrix(ZZ, 3, 3, range(9)); A
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: A.eigenspaces_left(format='galois')                                   # needs sage.rings.number_field
            [(0,
              Vector space of degree 3 and dimension 1 over Rational Field
              User basis matrix:
              [ 1 -2  1]),
             (a1,
              Vector space of degree 3 and dimension 1 over Number Field in a1 with defining polynomial x^2 - 12*x - 18
              User basis matrix:
              [            1 1/15*a1 + 2/5 2/15*a1 - 1/5])]

        We compute the left eigenspaces of the matrix of the Hecke operator
        `T_2` on level 43 modular symbols, both with all eigenvalues (the default)
        and with one subspace per factor. ::

            sage: # needs sage.modular
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
            [(3,
              Vector space of degree 7 and dimension 1 over Rational Field
              User basis matrix:
              [   1    0  1/7    0 -1/7    0 -2/7],
              1),
             (-2,
              Vector space of degree 7 and dimension 2 over Rational Field
              User basis matrix:
              [ 0  1  0  1 -1  1 -1]
              [ 0  0  1  0 -1  2 -1],
              2),
             (-1.414213562373095?,
              Vector space of degree 7 and dimension 2 over Algebraic Field
              User basis matrix:
              [                  0                   1                   0                  -1 0.4142135623730951?                   1                  -1]
              [                  0                   0                   1                   0                  -1                   0  2.414213562373095?],
              2),
             (1.414213562373095?,
              Vector space of degree 7 and dimension 2 over Algebraic Field
              User basis matrix:
              [                   0                    1                    0                   -1  -2.414213562373095?                    1                   -1]
              [                   0                    0                    1                    0                   -1                    0 -0.4142135623730951?],
              2)]
            sage: A.eigenspaces_left(format='galois', algebraic_multiplicity=True)
            [(3,
              Vector space of degree 7 and dimension 1 over Rational Field
              User basis matrix:
              [   1    0  1/7    0 -1/7    0 -2/7],
              1),
             (-2,
              Vector space of degree 7 and dimension 2 over Rational Field
              User basis matrix:
              [ 0  1  0  1 -1  1 -1]
              [ 0  0  1  0 -1  2 -1],
              2),
             (a2,
              Vector space of degree 7 and dimension 2 over Number Field in a2 with defining polynomial x^2 - 2
              User basis matrix:
              [      0       1       0      -1 -a2 - 1       1      -1]
              [      0       0       1       0      -1       0 -a2 + 1],
              2)]

        Next we compute the left eigenspaces over the finite field of order 11. ::

            sage: # needs sage.modular sage.rings.finite_rings
            sage: A = ModularSymbols(43, base_ring=GF(11), sign=1).T(2).matrix(); A
            [ 3  0  9  0]
            [ 0  9  0 10]
            [ 0  0 10  1]
            [ 0  0  1  1]
            sage: A.base_ring()
            Finite Field of size 11
            sage: A.charpoly()
            x^4 + 10*x^3 + 3*x^2 + 2*x + 1
            sage: A.eigenspaces_left(format='galois', var='beta')
            [(9,
              Vector space of degree 4 and dimension 1 over Finite Field of size 11
              User basis matrix:
              [0 1 5 6]),
             (3,
              Vector space of degree 4 and dimension 1 over Finite Field of size 11
              User basis matrix:
              [1 0 1 6]),
             (beta2,
              Vector space of degree 4 and dimension 1 over Univariate Quotient Polynomial Ring in beta2 over Finite Field of size 11 with modulus x^2 + 9
              User basis matrix:
              [        0         0         1 beta2 + 1])]

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
            NotImplementedError: eigenspaces cannot be computed reliably
            for inexact rings such as Real Field with 53 bits of precision,
            consult numerical or symbolic matrix classes for other options

            sage: # needs scipy
            sage: em = A.change_ring(RDF).eigenmatrix_left()
            sage: eigenvalues = em[0]; eigenvalues.dense_matrix()  # abs tol 1e-13
            [13.348469228349522                0.0                 0.0]
            [               0.0 -1.348469228349534                 0.0]
            [               0.0                0.0                 0.0]
            sage: eigenvectors = em[1]; eigenvectors # not tested
            [ 0.440242867...  0.567868371...  0.695493875...]
            [ 0.897878732...  0.278434036... -0.341010658...]
            [ 0.408248290... -0.816496580...  0.408248290...]

            sage: # needs sage.symbolic
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

            sage: # needs sage.rings.finite_rings
            sage: F.<b> = FiniteField(11^2)
            sage: A = matrix(F, [[b + 1, b + 1], [10*b + 4, 5*b + 4]])
            sage: A.eigenspaces_left(format='all')                                      # needs sage.rings.number_field
            Traceback (most recent call last):
            ...
            NotImplementedError: unable to construct eigenspaces for eigenvalues outside the base field,
            try the keyword option: format='galois'
            sage: A.eigenspaces_left(format='galois')                                   # needs sage.rings.number_field
            [(a0,
              Vector space of degree 2 and dimension 1 over Univariate Quotient Polynomial Ring in a0 over Finite Field in b of size 11^2 with modulus x^2 + (5*b + 6)*x + 8*b + 10
              User basis matrix:
              [               1 6*b*a0 + 3*b + 1])]

        TESTS:

        We make sure that :issue:`13308` is fixed. ::

            sage: M = ModularSymbols(Gamma1(23), sign=1)                                # needs sage.modular
            sage: m = M.cuspidal_subspace().hecke_matrix(2)                             # needs sage.modular
            sage: [j*m == i[0]*j                                # long time (4s)        # needs sage.modular
            ....:  for i in m.eigenspaces_left(format='all') for j in i[1].basis()]
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
        if algebraic_multiplicity not in [True, False]:
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
        if x is not None:
            if algebraic_multiplicity:
                return x
            else:
                return Sequence([(e[0], e[1]) for e in x], cr=True, check=False)

        # Possible improvements:
        # algorithm for dual_eigenvector in sage/modular/hecke/module.py
        # use minpoly when algebraic_multiplicity=False
        # as of 2007-03-25 minpoly is unreliable via linbox

        import sage.categories.homset
        import sage.rings.qqbar

        G = self.fcp()   # factored characteristic polynomial
        V = []
        i = -1  # variable name index, increments for each eigenvalue
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
                F = h.root_field('{0}{1}'.format(var, i))
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
            return Sequence([(e[0], e[1]) for e in V], cr=True, check=False)

    left_eigenspaces = eigenspaces_left

    def eigenspaces_right(self, format='all', var='a', algebraic_multiplicity=False):
        r"""
        Compute the right eigenspaces of a matrix.

        Note that ``eigenspaces_right()`` and ``right_eigenspaces()``
        are identical methods.  Here "right" refers to the eigenvectors
        being placed to the right of the matrix.

        INPUT:

        - ``self`` -- a square matrix over an exact field.  For inexact
          matrices consult the numerical or symbolic matrix classes

        - ``format`` -- (default: ``None``)

          - ``'all'`` -- attempts to create every eigenspace.  This will
            always be possible for matrices with rational entries
          - ``'galois'`` -- for each irreducible factor of the characteristic
            polynomial, a single eigenspace will be output for a
            single root/eigenvalue for the irreducible factor
          - ``None`` -- uses the 'all' format if the base ring is contained
            in an algebraically closed field which is implemented.
            Otherwise, uses the 'galois' format.

        - ``var`` -- (default: ``'a'``) variable name used to
          represent elements of the root field of each
          irreducible factor of the characteristic polynomial.
          If var='a', then the root fields will be in terms of
          a0, a1, a2, ...., where the numbering runs across all
          the irreducible factors of the characteristic polynomial,
          even for linear factors.

        - ``algebraic_multiplicity`` -- (default: ``False``) whether or
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

            sage: # needs sage.rings.number_field
            sage: A = matrix(QQ, 3, 3, range(9)); A
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: A.eigenspaces_right()
            [(0,
              Vector space of degree 3 and dimension 1 over Rational Field
              User basis matrix:
              [ 1 -2  1]),
             (-1.348469228349535?,
              Vector space of degree 3 and dimension 1 over Algebraic Field
              User basis matrix:
              [                   1  0.1303061543300932? -0.7393876913398137?]),
             (13.34846922834954?,
              Vector space of degree 3 and dimension 1 over Algebraic Field
              User basis matrix:
              [                 1 3.069693845669907? 5.139387691339814?])]
            sage: es = A.eigenspaces_right(format='galois'); es
            [(0,
              Vector space of degree 3 and dimension 1 over Rational Field
              User basis matrix:
              [ 1 -2  1]),
             (a1,
              Vector space of degree 3 and dimension 1 over Number Field in a1 with defining polynomial x^2 - 12*x - 18
              User basis matrix:
              [           1 1/5*a1 + 2/5 2/5*a1 - 1/5])]
            sage: es = A.eigenspaces_right(format='galois',
            ....:                          algebraic_multiplicity=True); es
            [(0,
              Vector space of degree 3 and dimension 1 over Rational Field
              User basis matrix:
              [ 1 -2  1],
              1),
             (a1,
              Vector space of degree 3 and dimension 1 over Number Field in a1 with defining polynomial x^2 - 12*x - 18
              User basis matrix:
              [           1 1/5*a1 + 2/5 2/5*a1 - 1/5],
              1)]
            sage: e, v, n = es[0]; v = v.basis()[0]
            sage: delta = v*e - A*v
            sage: abs(abs(delta)) < 1e-10
            True

        The same computation, but with implicit base change to a field::

            sage: A = matrix(ZZ, 3, range(9)); A
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: A.eigenspaces_right(format='galois')                                  # needs sage.rings.number_field
            [(0,
              Vector space of degree 3 and dimension 1 over Rational Field
              User basis matrix:
              [ 1 -2  1]),
             (a1,
              Vector space of degree 3 and dimension 1 over Number Field in a1 with defining polynomial x^2 - 12*x - 18
              User basis matrix:
              [           1 1/5*a1 + 2/5 2/5*a1 - 1/5])]

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
            NotImplementedError: eigenspaces cannot be computed reliably
            for inexact rings such as Real Field with 53 bits of precision,
            consult numerical or symbolic matrix classes for other options

            sage: # needs scipy
            sage: em = B.change_ring(RDF).eigenmatrix_right()
            sage: eigenvalues = em[0]; eigenvalues.dense_matrix()  # abs tol 1e-13
            [13.348469228349522                0.0                0.0]
            [               0.0 -1.348469228349534                0.0]
            [               0.0                0.0                0.0]
            sage: eigenvectors = em[1]; eigenvectors # not tested
            [ 0.164763817...  0.799699663...  0.408248290...]
            [ 0.505774475...  0.104205787... -0.816496580...]
            [ 0.846785134... -0.591288087...  0.408248290...]

            sage: # needs sage.symbolic
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
        if algebraic_multiplicity not in [True, False]:
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
        if x is not None:
            if algebraic_multiplicity:
                return x
            else:
                return Sequence([(e[0], e[1]) for e in x], cr=True, check=False)

        V = self.transpose().eigenspaces_left(format=format, var=var, algebraic_multiplicity=True)

        self.cache(key, V)
        if algebraic_multiplicity:
            return V
        else:
            return Sequence([(e[0], e[1]) for e in V], cr=True, check=False)

    right_eigenspaces = eigenspaces_right

    def eigenvalues(self, extend=True):
        r"""
        Return a sequence of the eigenvalues of a matrix, with
        multiplicity. If the eigenvalues are roots of polynomials in ``QQ``,
        then ``QQbar`` elements are returned that represent each separate
        root.

        If the option ``extend`` is set to ``False``, only eigenvalues in the base
        ring are considered.

        EXAMPLES::

            sage: a = matrix(ZZ, 4, range(16)); a
            [ 0  1  2  3]
            [ 4  5  6  7]
            [ 8  9 10 11]
            [12 13 14 15]
            sage: sorted(a.eigenvalues(), reverse=True)                                 # needs sage.rings.number_field
            [32.46424919657298?, 0, 0, -2.464249196572981?]

        ::

            sage: a = matrix([(1, 9, -1, -1),
            ....:             (-2, 0, -10, 2),
            ....:             (-1, 0, 15, -2),
            ....:             (0, 1, 0, -1)])
            sage: a.eigenvalues()                                                       # needs sage.rings.number_field
            [-0.9386318578049146?,
             15.50655435353258?,
             0.2160387521361705? - 4.713151979747493?*I,
             0.2160387521361705? + 4.713151979747493?*I]

        A symmetric matrix ``a + a.transpose()`` should have real eigenvalues

        ::

            sage: b = a + a.transpose()
            sage: ev = b.eigenvalues(); ev                                              # needs sage.rings.number_field
            [-8.35066086057957?, -1.107247901349379?,
             5.718651326708515?, 33.73925743522043?]

        The eigenvalues are elements of ``QQbar``, so they really represent
        exact roots of polynomials, not just approximations.

        ::

            sage: e = ev[0]; e                                                          # needs sage.rings.number_field
            -8.35066086057957?
            sage: p = e.minpoly(); p                                                    # needs sage.rings.number_field
            x^4 - 30*x^3 - 171*x^2 + 1460*x + 1784
            sage: p(e) == 0                                                             # needs sage.rings.number_field
            True

        To perform computations on the eigenvalue as an element of a number
        field, you can always convert back to a number field element.

        ::

            sage: e.as_number_field_element()                                           # needs sage.rings.number_field
            (Number Field in a
              with defining polynomial y^4 - 2*y^3 - 507*y^2 - 3972*y - 4264,
             a + 7,
             Ring morphism:
               From: Number Field in a with defining polynomial y^4 - 2*y^3 - 507*y^2 - 3972*y - 4264
               To:   Algebraic Real Field
               Defn: a |--> -15.35066086057957?)

        Notice the effect of the ``extend`` option.

        ::

            sage: M = matrix(QQ, [[0,-1,0], [1,0,0], [0,0,2]])
            sage: M.eigenvalues()                                                       # needs sage.rings.number_field
            [2, -1*I, 1*I]
            sage: M.eigenvalues(extend=False)                                           # needs sage.libs.pari
            [2]

        The method also works for matrices over finite fields::

            sage: M = matrix(GF(3), [[0,1,1], [1,2,0], [2,0,1]])
            sage: ev = sorted(M.eigenvalues()); ev                                      # needs sage.rings.finite_rings
            [2*z3, 2*z3 + 1, 2*z3 + 2]

        Similarly as in the case of ``QQbar``, the eigenvalues belong to some
        algebraic closure but they can be converted to elements of a finite
        field::

            sage: e = ev[0]                                                             # needs sage.rings.finite_rings
            sage: e.parent()                                                            # needs sage.rings.finite_rings
            Algebraic closure of Finite Field of size 3
            sage: e.as_finite_field_element()                                           # needs sage.rings.finite_rings
            (Finite Field in z3 of size 3^3,
             2*z3,
             Ring morphism:
               From: Finite Field in z3 of size 3^3
               To:   Algebraic closure of Finite Field of size 3
               Defn: z3 |--> z3)
        """
        x = self.fetch('eigenvalues')
        if x is not None:
            if not extend:
                x = Sequence(i for i in x if i in self.base_ring())
            return x

        if not self.base_ring().is_exact():
            from warnings import warn
            warn("Using generic algorithm for an inexact ring, which will probably give incorrect results due to numerical precision issues.")

        if not extend:
            return Sequence(r for r, m in self.charpoly().roots() for _ in range(m))

        # now we need to find a natural algebraic closure for the base ring
        K = self.base_ring()

        if K not in _Fields:
            if not K.is_integral_domain():
                raise NotImplementedError("eigenvalues() not implemented for non integral domains")
            K = K.fraction_field()

        try:
            A = K.algebraic_closure()
        except (AttributeError, ValueError):
            raise NotImplementedError("algebraic closure is not implemented for %s" % K)

        res = []
        for f, e in self.charpoly().change_ring(K).factor():
            if f.degree() == 1:
                res.extend([-f.constant_coefficient()]*e)
            else:
                for r, ee in f.change_ring(A).roots():
                    res.extend([r]*(e*ee))

        eigenvalues = Sequence(res)
        self.cache('eigenvalues', eigenvalues)
        return eigenvalues

    def eigenvectors_left(self, other=None, *, extend=True):
        r"""
        Compute the left eigenvectors of a matrix.

        INPUT:

        - ``other`` -- a square matrix `B` (default: ``None``) in a generalized
          eigenvalue problem; if ``None``, an ordinary eigenvalue problem is
          solved (currently supported only if the base ring of ``self`` is
          ``RDF`` or ``CDF``)

        - ``extend`` -- boolean (default: ``True``)

        OUTPUT:

        For each distinct eigenvalue, returns a list of the form (e,V,n)
        where e is the eigenvalue, V is a list of eigenvectors forming a
        basis for the corresponding left eigenspace, and n is the algebraic
        multiplicity of the eigenvalue.

        If the option extend is set to False, then only the eigenvalues that
        live in the base ring are considered.

        EXAMPLES:

        We compute the left eigenvectors of a `3\times 3` rational matrix.

        ::

            sage: # needs sage.rings.number_field
            sage: A = matrix(QQ, 3, 3, range(9)); A
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: es = A.eigenvectors_left(); es
            [(0, [(1, -2, 1)], 1),
             (-1.348469228349535?, [(1, 0.3101020514433644?, -0.3797958971132713?)], 1),
             (13.34846922834954?, [(1, 1.289897948556636?, 1.579795897113272?)], 1)]
            sage: eval, [evec], mult = es[0]
            sage: delta = eval*evec - evec*A
            sage: abs(abs(delta)) < 1e-10
            True

        Notice the difference between considering ring extensions or not.

        ::

            sage: M = matrix(QQ, [[0,-1,0], [1,0,0], [0,0,2]])
            sage: M.eigenvectors_left()                                                 # needs sage.rings.number_field
            [(2, [(0, 0, 1)], 1), (-1*I, [(1, -1*I, 0)], 1), (1*I, [(1, 1*I, 0)], 1)]
            sage: M.eigenvectors_left(extend=False)                                     # needs sage.rings.number_field
            [(2, [(0, 0, 1)], 1)]

        TESTS::

            sage: A = matrix(QQ, [[1, 2], [3, 4]])
            sage: B = matrix(QQ, [[1, 1], [0, 1]])
            sage: A.eigenvectors_left(B)
            Traceback (most recent call last):
            ...
            NotImplementedError: generalized eigenvector decomposition is
            implemented for RDF and CDF, but not for Rational Field

        Check the deprecation::

            sage: matrix(QQ, [[1, 2], [3, 4]]).eigenvectors_left(False)                 # needs sage.rings.number_field
            doctest:...: DeprecationWarning: "extend" should be used as keyword argument
            See https://github.com/sagemath/sage/issues/29243 for details.
            []

        Check :issue:`30518`::

            sage: # needs sage.rings.number_field
            sage: K.<i> = QuadraticField(-1)
            sage: m = matrix(K, 4, [2,4*i,-i,0, -4*i,2,-1,0, 2*i,-2,0,0, 4*i+4, 4*i-4,1-i,-2])
            sage: assert all(m*v == e*v for e, vs, _ in m.eigenvectors_right() for v in vs)
        """
        if other is not None:
            if isinstance(other, bool):
                # for backward compatibility
                from sage.misc.superseded import deprecation
                deprecation(29243,
                            '"extend" should be used as keyword argument')
                extend = other
                other = None
            else:
                raise NotImplementedError('generalized eigenvector '
                                          'decomposition is implemented '
                                          'for RDF and CDF, but not for %s'
                                          % self.base_ring())

        x = self.fetch('eigenvectors_left')
        if x is not None:
            return x

        if not self.base_ring().is_exact():
            from warnings import warn
            warn("Using generic algorithm for an inexact ring, which may result in garbage from numerical precision issues.")

        from sage.categories.homset import hom
        eigenspaces = self.eigenspaces_left(format='galois', algebraic_multiplicity=True)
        evec_list = []
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
                        from sage.rings.qqbar import QQbar
                        eigval_conj = eigval.galois_conjugates(QQbar)
                    except AttributeError:
                        raise NotImplementedError("eigenvectors are not implemented for matrices with eigenvalues that are not in the fraction field of the base ring or in QQbar")

                    for e in eigval_conj:
                        m = hom(eigval.parent(), e.parent(), e)
                        space = (e.parent())**n
                        evec_list = [(space)([m(i) for i in v]) for v in eigbasis]
                        evec_eval_list.append((e, evec_list, eigmult))

        return evec_eval_list

    left_eigenvectors = eigenvectors_left

    def eigenvectors_right(self, other=None, *, extend=True):
        r"""
        Compute the right eigenvectors of a matrix.

        INPUT:

        - ``other`` -- a square matrix `B` (default: ``None``) in a generalized
          eigenvalue problem; if ``None``, an ordinary eigenvalue problem is
          solved (currently supported only if the base ring of ``self`` is
          ``RDF`` or ``CDF``)

        - ``extend`` -- boolean (default: ``True``)

        OUTPUT:

        For each distinct eigenvalue, returns a list of the form (e,V,n)
        where e is the eigenvalue, V is a list of eigenvectors forming a
        basis for the corresponding right eigenspace, and n is the
        algebraic multiplicity of the eigenvalue. If ``extend = True``
        (the default), this will return eigenspaces over the algebraic
        closure of the base field where this is implemented; otherwise
        it will restrict to eigenvalues in the base field.

        EXAMPLES:

        We compute the right eigenvectors of a `3\times 3` rational matrix.

        ::

            sage: # needs sage.rings.number_field
            sage: A = matrix(QQ, 3, 3, range(9)); A
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: es = A.eigenvectors_right(); es
            [(0, [(1, -2, 1)], 1),
             (-1.348469228349535?, [(1, 0.1303061543300932?, -0.7393876913398137?)], 1),
             (13.34846922834954?, [(1, 3.069693845669907?, 5.139387691339814?)], 1)]
            sage: A.eigenvectors_right(extend=False)
            [(0, [(1, -2, 1)], 1)]
            sage: eval, [evec], mult = es[0]
            sage: delta = eval*evec - A*evec
            sage: abs(abs(delta)) < 1e-10
            True

        TESTS::

            sage: A = matrix(QQ, [[1, 2], [3, 4]])
            sage: B = matrix(QQ, [[1, 1], [0, 1]])
            sage: A.eigenvectors_right(B)
            Traceback (most recent call last):
            ...
            NotImplementedError: generalized eigenvector decomposition is
            implemented for RDF and CDF, but not for Rational Field
        """
        return self.transpose().eigenvectors_left(other=other, extend=extend)

    right_eigenvectors = eigenvectors_right

    def eigenmatrix_left(self, other=None):
        r"""
        Return matrices `D` and `P`, where `D` is a diagonal matrix of
        eigenvalues and the rows of `P` are corresponding eigenvectors
        (or zero vectors).

        INPUT:

        - ``other`` -- a square matrix `B` (default: ``None``) in a generalized
          eigenvalue problem; if ``None``, an ordinary eigenvalue problem is
          solved

        OUTPUT:

        If ``self`` is a square matrix `A`, then the output is a diagonal
        matrix `D` and a matrix `P` such that

        .. MATH::

            P A = D P,

        where the rows of `P` are eigenvectors of `A` and the diagonal entries
        of `D` are the corresponding eigenvalues.

        If a matrix `B` is passed as optional argument, the output is a
        solution to the generalized eigenvalue problem such that

        .. MATH::

            P A = D P B.

        The ordinary eigenvalue problem is equivalent to the generalized one if
        `B` is the identity matrix.

        The generalized eigenvector decomposition is currently only implemented
        for matrices over ``RDF`` and ``CDF``.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: A = matrix(QQ, 3, 3, range(9)); A
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

        Because `P` is invertible, `A` is diagonalizable.

        ::

            sage: A == (~P)*D*P                                                         # needs sage.rings.number_field
            True

        The matrix `P` may contain zero rows corresponding to eigenvalues for
        which the algebraic multiplicity is greater than the geometric
        multiplicity. In these cases, the matrix is not diagonalizable.

        ::

            sage: # needs sage.rings.number_field
            sage: A = jordan_block(2, 3); A
            [2 1 0]
            [0 2 1]
            [0 0 2]
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

        A generalized eigenvector decomposition::

            sage: # needs scipy
            sage: A = matrix(RDF, [[1, -2], [3, 4]])
            sage: B = matrix(RDF, [[0, 7], [2, -3]])
            sage: D, P = A.eigenmatrix_left(B)
            sage: (P * A - D * P * B).norm() < 1e-14
            True

        The matrix `B` in a generalized eigenvalue problem may be singular::

            sage: # needs scipy sage.rings.complex_double sage.symbolic
            sage: A = matrix.identity(CDF, 2)
            sage: B = matrix(CDF, [[2, 1+I], [4, 2+2*I]])
            sage: D, P = A.eigenmatrix_left(B)
            sage: D.diagonal()  # tol 1e-14
            [0.2 - 0.1*I, +infinity]

        In this case, we can still verify the eigenvector equation for the
        first eigenvalue and first eigenvector::

            sage: # needs scipy sage.rings.complex_double sage.symbolic
            sage: l = D[0, 0]
            sage: v = P[0, :]
            sage: (v * A - l * v * B).norm() < 1e-14
            True

        The second eigenvector is contained in the left kernel of `B`::

            sage: (P[1, :] * B).norm() < 1e-14                                          # needs scipy sage.rings.complex_double sage.symbolic
            True

        .. SEEALSO::

            :meth:`eigenvalues`,
            :meth:`eigenvectors_left`,
            :meth:`.Matrix_double_dense.eigenvectors_left`,
            :meth:`eigenmatrix_right`.

        TESTS:

        For matrices with floating point entries, some platforms will
        return eigenvectors that are negatives of those returned by the
        majority of platforms.  This test accounts for that possibility.
        Running this test independently, without adjusting the eigenvectors
        could indicate this situation on your hardware.  ::

            sage: A = matrix(QQ, 3, 3, range(9))
            sage: em = A.change_ring(RDF).eigenmatrix_left()
            sage: evalues = em[0]; evalues.dense_matrix()  # abs tol 1e-13
            [13.348469228349522                0.0                 0.0]
            [               0.0 -1.348469228349534                 0.0]
            [               0.0                0.0                 0.0]
            sage: evectors = em[1]
            sage: for i in range(3):
            ....:     scale = evectors[i,0].sign()
            ....:     evectors.rescale_row(i, scale)
            sage: evectors  # abs tol 1e-13
            [0.44024286723591904  0.5678683713143027  0.6954938753926869]
            [ 0.8978787322617111 0.27843403682172374 -0.3410106586182631]
            [ 0.4082482904638625 -0.8164965809277263 0.40824829046386324]

        The following example shows that :issue:`20439` has been resolved::

            sage: # needs sage.rings.complex_double sage.symbolic
            sage: A = matrix(CDF, [[-2.53634347567,  2.04801738686, -0.0, -62.166145304],
            ....:                  [ 0.7, -0.6, 0.0, 0.0],
            ....:                  [0.547271128842, 0.0, -0.3015, -21.7532081652],
            ....:                  [0.0, 0.0, 0.3, -0.4]])
            sage: D, P = A.eigenmatrix_left()
            sage: (P*A - D*P).norm() < 10^(-2)
            True

        The following example shows that the fix for :issue:`20439` (conjugating
        eigenvectors rather than eigenvalues) is the correct one::

            sage: # needs sage.rings.complex_double sage.symbolic
            sage: A = Matrix(CDF,[[I,0],[0,1]])
            sage: D, P = A.eigenmatrix_left()
            sage: (P*A - D*P).norm() < 10^(-2)
            True

        For some symbolic matrices, the Maxima backend fails to correctly
        compute some eigenvectors, returning either none or more vectors than
        the algebraic multiplicity. The following examples show that these
        cases are detected (:issue:`27842`)::

            sage: # needs sage.symbolic
            sage: A = matrix(SR, [(225/548, 0, -175/274*sqrt(193/1446)),
            ....:                 (0, 1/2, 0),
            ....:                 (-63/548*sqrt(723/386), 0, 49/548)])
            sage: A.eigenmatrix_left()
            Traceback (most recent call last):
            ...
            RuntimeError: failed to compute eigenvectors for eigenvalue ...,
            check eigenvectors_left() for partial results
            sage: B = matrix(SR, [(1/2, -7/2*sqrt(1/386), 0, 49/2*sqrt(1/279078)),
            ....:                 (-7/2*sqrt(1/386), 211/772, 0, -8425/772*sqrt(1/723)),
            ....:                 (0, 0, 1/2, 0),
            ....:                 (49/2*sqrt(1/279078), -8425/772*sqrt(1/723), 0, 561/772)])
            sage: B.eigenmatrix_left()          # long time (1.2 seconds)
            Traceback (most recent call last):
            ...
            RuntimeError: failed to compute eigenvectors for eigenvalue ...,
            check eigenvectors_left() for partial results

        The following example shows that :issue:`12595` has been resolved::

            sage: # needs sage.rings.complex_double
            sage: m = Matrix(CDF, 8, [[-1, -1, -1, -1, 1, -3, -1, -1],
            ....:                     [1, 1, 1, 1, -1, -1, 1, -3],
            ....:                     [-1, 3, -1, -1, 1, 1, -1, -1],
            ....:                     [-1, -1, -1, 3, 1, 1, -1, -1],
            ....:                     [1, 1, -3, 1, -1, -1, 1, 1],
            ....:                     [1, 1, 1, 1, -1, -1, -3, 1],
            ....:                     [3, -1, -1, -1, 1, 1, -1, -1],
            ....:                     [1, 1, 1, 1, 3, -1, 1, 1]])
            sage: d, p = m.eigenmatrix_left()
            sage: mm = p.inverse() * d * p  # This should be equal to m up to numerical noise
            sage: norm(mm - m) < 1e-12
            True
        """
        from sage.misc.flatten import flatten
        from sage.matrix.constructor import diagonal_matrix, matrix
        evecs = self.eigenvectors_left(other=other)
        D = diagonal_matrix(flatten([[e[0]]*e[2] for e in evecs]))
        rows = []
        for e in evecs:
            defect = e[2] - len(e[1])
            if e[1] and defect >= 0:
                rows.extend(e[1] + [e[1][0].parent().zero_vector()] * defect)
            else:
                # see Issue #27842
                raise RuntimeError(
                        "failed to compute eigenvectors for eigenvalue %s, "
                        "check eigenvectors_left() for partial results" % e[0])
        P = matrix(rows)
        return D, P

    left_eigenmatrix = eigenmatrix_left

    def eigenmatrix_right(self, other=None):
        r"""
        Return matrices `D` and `P`, where `D` is a diagonal matrix of
        eigenvalues and the columns of `P` are corresponding eigenvectors
        (or zero vectors).

        INPUT:

        - ``other`` -- a square matrix `B` (default: ``None``) in a generalized
          eigenvalue problem; if ``None``, an ordinary eigenvalue problem is
          solved

        OUTPUT:

        If ``self`` is a square matrix `A`, then the output is a diagonal
        matrix `D` and a matrix `P` such that

        .. MATH::

            A P = P D,

        where the columns of `P` are eigenvectors of `A` and the diagonal
        entries of `D` are the corresponding eigenvalues.

        If a matrix `B` is passed as optional argument, the output is a
        solution to the generalized eigenvalue problem such that

        .. MATH::

            A P = B P D.

        The ordinary eigenvalue problem is equivalent to the generalized one if
        `B` is the identity matrix.

        The generalized eigenvector decomposition is currently only implemented
        for matrices over ``RDF`` and ``CDF``.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: A = matrix(QQ, 3, 3, range(9)); A
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

        Because `P` is invertible, `A` is diagonalizable.

        ::

            sage: A == P*D*(~P)                                                         # needs sage.rings.number_field
            True

        The matrix `P` may contain zero columns corresponding to eigenvalues
        for which the algebraic multiplicity is greater than the geometric
        multiplicity. In these cases, the matrix is not diagonalizable.

        ::

            sage: # needs sage.rings.number_field
            sage: A = jordan_block(2, 3); A
            [2 1 0]
            [0 2 1]
            [0 0 2]
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

        A generalized eigenvector decomposition::

            sage: # needs scipy
            sage: A = matrix(RDF, [[1, -2], [3, 4]])
            sage: B = matrix(RDF, [[0, 7], [2, -3]])
            sage: D, P = A.eigenmatrix_right(B)
            sage: (A * P - B * P * D).norm() < 1e-14
            True

        The matrix `B` in a generalized eigenvalue problem may be singular::

            sage: # needs scipy
            sage: A = matrix.identity(RDF, 2)
            sage: B = matrix(RDF, [[3, 5], [6, 10]])
            sage: D, P = A.eigenmatrix_right(B); D   # tol 1e-14
            [0.07692307692307694                 0.0]
            [                0.0           +infinity]

        In this case, we can still verify the eigenvector equation for the
        first eigenvalue and first eigenvector::

            sage: # needs scipy
            sage: l = D[0, 0]
            sage: v = P[:, 0]
            sage: (A * v  - B * v * l).norm() < 1e-14
            True

        The second eigenvector is contained in the right kernel of `B`::

            sage: (B * P[:, 1]).norm() < 1e-14                                          # needs scipy
            True

        .. SEEALSO::

            :meth:`eigenvalues`,
            :meth:`eigenvectors_right`,
            :meth:`.Matrix_double_dense.eigenvectors_right`,
            :meth:`eigenmatrix_left`.

        TESTS:

        For matrices with floating point entries, some platforms will
        return eigenvectors that are negatives of those returned by the
        majority of platforms.  This test accounts for that possibility.
        Running this test independently, without adjusting the eigenvectors
        could indicate this situation on your hardware.  ::

            sage: # needs scipy
            sage: B = matrix(QQ, 3, 3, range(9))
            sage: em = B.change_ring(RDF).eigenmatrix_right()
            sage: evalues = em[0]; evalues.dense_matrix()  # abs tol 1e-13
            [13.348469228349522                0.0                0.0]
            [               0.0 -1.348469228349534                0.0]
            [               0.0                0.0                0.0]
            sage: evectors = em[1]
            sage: for i in range(3):
            ....:     scale = evectors[0,i].sign()
            ....:     evectors.rescale_col(i, scale)
            sage: evectors  # abs tol 1e-13
            [ 0.1647638172823028   0.799699663111865 0.40824829046386285]
            [ 0.5057744759005657 0.10420578771917821 -0.8164965809277261]
            [ 0.8467851345188293 -0.5912880876735089  0.4082482904638632]

        The following example shows that :issue:`20439` has been resolved::

            sage: # needs scipy sage.rings.complex_double
            sage: A = matrix(CDF, [[-2.53634347567,  2.04801738686, -0.0, -62.166145304],
            ....:                  [ 0.7, -0.6, 0.0, 0.0],
            ....:                  [0.547271128842, 0.0, -0.3015, -21.7532081652],
            ....:                  [0.0, 0.0, 0.3, -0.4]])
            sage: D, P = A.eigenmatrix_right()
            sage: (A*P - P*D).norm() < 10^(-2)
            True

        The following example shows that the fix for :issue:`20439` (conjugating
        eigenvectors rather than eigenvalues) is the correct one::

            sage: # needs scipy sage.rings.complex_double sage.symbolic
            sage: A = Matrix(CDF,[[I,0],[0,1]])
            sage: D, P = A.eigenmatrix_right()
            sage: (A*P - P*D).norm() < 10^(-2)
            True
        """
        D, P = self.transpose().eigenmatrix_left(None if other is None
                                                 else other.transpose())
        return D, P.transpose()

    right_eigenmatrix = eigenmatrix_right

    def eigenvalue_multiplicity(self, s):
        r"""
        Return the multiplicity of ``s`` as a generalized eigenvalue
        of the matrix.

        EXAMPLES::

            sage: M = Matrix(QQ, [[0,1],[0,0]])
            sage: M.eigenvalue_multiplicity(0)
            2
            sage: M.eigenvalue_multiplicity(1)
            0

            sage: M = posets.DiamondPoset(5).coxeter_transformation()                   # needs sage.graphs sage.libs.flint
            sage: [M.eigenvalue_multiplicity(x) for x in [-1, 1]]                       # needs sage.graphs sage.libs.flint
            [3, 2]

        TESTS::

            sage: M = Matrix(QQ, [[0,1,2],[0,0,0]])
            sage: M.eigenvalue_multiplicity(1)
            Traceback (most recent call last):
            ...
            TypeError: matrix must be square, not 2 x 3
        """
        if not self.is_square():
            msg = 'matrix must be square, not {0} x {1}'
            raise TypeError(msg.format(self.nrows(), self.ncols()))

        r = self.dimensions()[0]
        n = 0
        m1 = self - s
        while True:
            m1 *= m1
            m2 = m1.extended_echelon_form(subdivide=True)
            t = r - m2.subdivisions()[0][0]
            n += t
            if t == 0 or t == r:
                return n
            m3 = m2.subdivision(0, 0)
            m4 = m2.submatrix(0, r, r, r)
            m5 = m3 * m4.inverse()
            m1 = m5.submatrix(0, 0, r - t, r - t)
            r -= t

    #####################################################################################
    # Generic Echelon Form
    #####################################################################################

    def rref(self, *args, **kwds):
        """
        Return the reduced row echelon form of the matrix, considered
        as a matrix over a field.

        If the matrix is over a ring, then an equivalent matrix is
        constructed over the fraction field, and then row reduced.

        All arguments are passed on to :meth:`echelon_form`.

        .. NOTE::

            Because the matrix is viewed as a matrix over a field,
            every leading coefficient of the returned matrix will be
            one and will be the only nonzero entry in its column.

        EXAMPLES::

            sage: A = matrix(3,range(9)); A
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

            sage: B = random_matrix(QQ, 3, num_bound=10)
            sage: while B.rank() != 3:
            ....:     B = random_matrix(QQ, 3, num_bound=10)
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

            sage: R.<x> = ZZ[]
            sage: C = matrix(3, [2,x,x^2, x+1,3-x,-1, 3,2,1])
            sage: C.rref()
            [1 0 0]
            [0 1 0]
            [0 0 1]
            sage: C.base_ring()
            Univariate Polynomial Ring in x over Integer Ring
            sage: C.echelon_form()
            Traceback (most recent call last):
            ...
            NotImplementedError: Ideal Ideal (2, x + 1) of Univariate
            Polynomial Ring in x over Integer Ring not principal
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
            sage: C = matrix(3, [2,x,x^2, x+1,3-x,-1/x, 3,2,1/2])
            sage: C.echelon_form()
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        R = self.base_ring()
        if R in _Fields:
            return self.echelon_form()
        else:
            F = R.fraction_field()
            return self.change_ring(F).echelon_form()

    def _echelonize_ring(self, **kwds):
        r"""
        Echelonize ``self`` in place, where the base ring of ``self`` is assumed to
        be a ring (not a field).

        Right now this *only* works over ZZ and some principal ideal domains;
        otherwise a :exc:`NotImplementedError` is raised. In the special case of
        sparse matrices over ZZ it makes them dense, gets the echelon form of
        the dense matrix, then sets ``self`` equal to the result.

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

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: L.<w> = NumberField(x^2 - x + 2)
            sage: OL = L.ring_of_integers()
            sage: m = matrix(OL, 2, 2, [1,2,3,4+w])
            sage: m.echelon_form()
            [    1     2]
            [    0 w - 2]
            sage: E, T = m.echelon_form(transformation=True); E, T
            (
            [    1     2]  [ 1  0]
            [    0 w - 2], [-3  1]
            )
            sage: E == T*m
            True

        TESTS:

        Check that :issue:`11558` is fixed::

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
        cdef bint transformation = 'transformation' in kwds and kwds['transformation']
        if self._base_ring == ZZ:
            if 'include_zero_rows' in kwds and not kwds['include_zero_rows']:
                raise ValueError("cannot echelonize in place and delete zero rows")
            if transformation:
                d, a = self.dense_matrix().echelon_form(**kwds)
            else:
                d = self.dense_matrix().echelon_form(**kwds)
            for c from 0 <= c < self.ncols():
                for r from 0 <= r < self.nrows():
                    self.set_unsafe(r, c, d.get_unsafe(r, c))
            self.clear_cache()
            self.cache('pivots', d.pivots())
            self.cache('in_echelon_form', True)
        else:
            try:
                a, d, p = self._echelon_form_PID()
            except TypeError as msg:
                raise NotImplementedError("%s\nechelon form over %s not yet implemented"%(msg, self.base_ring()))

            for c from 0 <= c < self.ncols():
                for r from 0 <= r < self.nrows():
                    self.set_unsafe(r, c, d.get_unsafe(r, c))
            self.clear_cache()
            self.cache('pivots', tuple(p))
            self.cache('in_echelon_form', True)
        if transformation:
            return a
        else:
            return

    def echelonize(self, algorithm='default', cutoff=0, **kwds):
        r"""
        Transform ``self`` into a matrix in echelon form over the same
        base ring as ``self``.

        .. NOTE::

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

          - ``'partial_pivoting'``: Gauss elimination, using partial pivoting
            (if base ring has absolute value)

          - ``'scaled_partial_pivoting'``: Gauss elimination, using scaled
            partial pivoting (if base ring has absolute value)

          - ``'scaled_partial_pivoting_valuation'``: Gauss elimination, using
            scaled partial pivoting (if base ring has valuation)

          - ``'strassen'``: use a Strassen divide and conquer
            algorithm (if available)

        - ``cutoff`` -- integer; only used if the Strassen algorithm
          is selected

        - ``transformation`` -- boolean; whether to also return the
          transformation matrix. Some matrix backends do not provide
          this information, in which case this option is ignored.

        OUTPUT:

        The matrix ``self`` is put into echelon form. Nothing is
        returned unless the keyword option ``transformation=True`` is
        specified, in which case the transformation matrix is
        returned.

        EXAMPLES::

            sage: a = matrix(QQ, 3,3, range(9)); a
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

            sage: a = matrix(QQ, 3,3, range(9)); a.set_immutable()
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

            sage: a = matrix(ZZ, 3,3, range(9))
            sage: a.echelonize(); a
            [ 3  0 -3]
            [ 0  1  2]
            [ 0  0  0]

        We compute an echelon form both over a domain and fraction field::

            sage: R.<x,y> = QQ[]
            sage: a = matrix(R, 2, [x,y, x,y])
            sage: a.echelon_form()               # not very useful? -- why two copies of the same row?                  # needs sage.rings.function_field
            [x y]
            [x y]

        ::

            sage: b = a.change_ring(R.fraction_field())
            sage: b.echelon_form()               # potentially useful
            [  1 y/x]
            [  0   0]

        We check that the echelon form works for matrices over `p`-adics.
        See :issue:`17272`::

            sage: # needs sage.rings.padics
            sage: R = ZpCA(5,5,print_mode='val-unit')
            sage: A = matrix(R, 3,3, [250,2369,1147, 106,927,362, 90,398,2483])
            sage: A
            [5^3 * 2 + O(5^5)    2369 + O(5^5)    1147 + O(5^5)]
            [    106 + O(5^5)     927 + O(5^5)     362 + O(5^5)]
            [ 5 * 18 + O(5^5)     398 + O(5^5)    2483 + O(5^5)]
            sage: K = R.fraction_field()
            sage: A.change_ring(K).augment(identity_matrix(K,3)).echelon_form()
            [      1 + O(5^5)           O(5^5)           O(5^5) 5 * 212 + O(5^5)    3031 + O(5^5)    2201 + O(5^5)]
            [          O(5^5)       1 + O(5^5)           O(5^5)    1348 + O(5^5) 5 * 306 + O(5^5)    2648 + O(5^5)]
            [          O(5^5)           O(5^5)       1 + O(5^5)    1987 + O(5^5) 5 * 263 + O(5^5)     154 + O(5^5)]

        Echelon form is not defined over arbitrary rings::

            sage: a = matrix(Integers(9), 3,3, range(9))
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

        TESTS:

        Check that :issue:`34724` is fixed (indirect doctest)::

            sage: a = 6.12323399573677e-17
            sage: m = matrix(RR,[[-a, -1.72508242466029], [ 0.579682446302195, a]])
            sage: (~m*m).norm()
            1.0
        """
        self.check_mutability()
        basring = self.base_ring()

        if algorithm == 'default':
            from sage.categories.discrete_valuation import DiscreteValuationFields
            if self._will_use_strassen_echelon():
                algorithm = 'strassen'
            # Currently we only use scaled partial pivoting in discrete valuation fields
            # In general, we would like to do so in any rank one valuation ring,
            # but this should be done by introducing a category of general valuation rings and fields,
            # which we don't have at the moment
            elif basring in DiscreteValuationFields():
                try:
                    basring.one().abs()
                    algorithm = 'scaled_partial_pivoting'
                except (AttributeError, TypeError):
                    algorithm = 'scaled_partial_pivoting_valuation'
            elif basring.is_exact():
                algorithm = 'classical'
            else:
                try:
                    (basring(0.42)).abs()
                    algorithm = 'scaled_partial_pivoting'
                except (AttributeError, ArithmeticError, TypeError):
                    algorithm = 'classical'
        try:
            if basring in _Fields:
                if algorithm in ['classical', 'partial_pivoting', 'scaled_partial_pivoting', 'scaled_partial_pivoting_valuation']:
                    self._echelon_in_place(algorithm)
                elif algorithm == 'strassen':
                    self._echelon_strassen(cutoff)
                else:
                    raise ValueError("Unknown algorithm '%s'" % algorithm)
            else:
                if not (algorithm in ['classical', 'strassen']):
                    kwds['algorithm'] = algorithm
                return self._echelonize_ring(**kwds)
        except ArithmeticError as msg:
            raise NotImplementedError("%s\nEchelon form not implemented over '%s'."%(msg, basring))

    def echelon_form(self, algorithm='default', cutoff=0, **kwds):
        r"""
        Return the echelon form of ``self``.

        .. NOTE::

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

          - ``'partial_pivoting'``: Gauss elimination, using partial pivoting
            (if base ring has absolute value)

          - ``'scaled_partial_pivoting'`` -- Gauss elimination, using scaled
            partial pivoting (if base ring has absolute value)

          - ``'scaled_partial_pivoting_valuation'``: Gauss elimination, using
            scaled partial pivoting (if base ring has valuation)

          - ``'strassen'``: use a Strassen divide and conquer
            algorithm (if available)

        - ``cutoff`` -- integer; only used if the Strassen algorithm is selected

        - ``transformation`` -- boolean; whether to also return the
          transformation matrix. Some matrix backends do not provide
          this information, in which case this option is ignored.

        OUTPUT:

        The reduced row echelon form of ``self``, as an immutable
        matrix. Note that ``self`` is *not* changed by this command. Use
        :meth:`echelonize` to change ``self`` in place.

        If the optional parameter ``transformation=True`` is
        specified, the output consists of a pair `(E,T)` of matrices
        where `E` is the echelon form of ``self`` and `T` is the
        transformation matrix.

        EXAMPLES::

            sage: MS = MatrixSpace(GF(19), 2, 3)
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
        cdef bint transformation = ('transformation' in kwds and kwds['transformation'])
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
            v = E.echelonize(algorithm=algorithm, cutoff=cutoff, **kwds)
        E.set_immutable()  # so we can cache the echelon form.
        self.cache('echelon_form', E)
        if v is not None:
            self.cache('echelon_transformation', v)
        self.cache('pivots', E.pivots())
        if transformation and v is not None:
            return (E, v)
        else:
            return E

    cpdef _echelon(self, str algorithm):
        """
        Return the echelon form of ``self`` using ``algorithm``.

        EXAMPLES::

            sage: t = matrix(QQ, 3, 3, range(9)); t
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: t._echelon('classical')
            [ 1  0 -1]
            [ 0  1  2]
            [ 0  0  0]
            sage: a = matrix(QQ, 2, [1..6])
            sage: a._echelon('classical')
            [ 1  0 -1]
            [ 0  1  2]

            sage: # needs sage.rings.padics
            sage: R = ZpCA(5, 5, print_mode='val-unit')
            sage: A = matrix(R, 3, 3, [250,2369,1147, 106,927,362, 90,398,2483]); A
            [5^3 * 2 + O(5^5)    2369 + O(5^5)    1147 + O(5^5)]
            [    106 + O(5^5)     927 + O(5^5)     362 + O(5^5)]
            [ 5 * 18 + O(5^5)     398 + O(5^5)    2483 + O(5^5)]
            sage: A._echelon('partial_pivoting')
            [1 + O(5^5)     O(5^5)     O(5^5)]
            [    O(5^5) 1 + O(5^5)     O(5^5)]
            [    O(5^5)     O(5^5) 1 + O(5^5)]

        The following example is an example where partial pivoting fails,
        but scaled partial pivoting succeeds, taken from 'Numerical Analysis (9th edition)'
        by R.L. Burden and J.D. Faires (with minor adjustments)::

            sage: RR13 = RealField(prec=13)
            sage: A = Matrix(RR13, 2, 3, [30, 591400, 591700, 5.291, -6.130, 46.78])
            sage: A
            [   30.0 591000. 592000.]
            [   5.29   -6.13    46.8]
            sage: A._echelon('classical')
            [ 1.00 0.000  12.0]
            [0.000  1.00  1.00]
            sage: A._echelon('partial_pivoting')
            [ 1.00 0.000  12.0]
            [0.000  1.00  1.00]
            sage: A._echelon('scaled_partial_pivoting')
            [ 1.00 0.000  10.0]
            [0.000  1.00  1.00]
        """
        E = self.fetch('echelon_' + algorithm)
        if E is not None:
            return E
        E = self.__copy__()
        E._echelon_in_place(algorithm)
        self.cache('echelon_' + algorithm, E)
        return E

    # for backward compatibility

    def _echelon_classical(self):
        """
        Return the echelon form of ``self`` using the classical algorithm.

        EXAMPLES::

            sage: t = matrix(QQ, 3, 3, range(9)); t
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: t._echelon_classical()
            [ 1  0 -1]
            [ 0  1  2]
            [ 0  0  0]
            sage: a = matrix(QQ, 2, [1..6])
            sage: a._echelon_classical()
            [ 1  0 -1]
            [ 0  1  2]
        """
        return self._echelon('classical')

    cpdef _echelon_in_place(self, str algorithm):
        """
        Transform ``self`` into echelon form and return the pivots of ``self``.

        EXAMPLES::

            sage: t = matrix(QQ, 3, 3, range(9)); t
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: E = t._echelon_in_place('classical'); t
            [ 1  0 -1]
            [ 0  1  2]
            [ 0  0  0]
            sage: a = matrix(QQ,2,[1..6])
            sage: P = a._echelon_in_place('classical'); a
            [ 1  0 -1]
            [ 0  1  2]

            sage: # needs sage.rings.padics
            sage: R = ZpCA(5, 5, print_mode='val-unit')
            sage: A = matrix(R,3,3,[250,2369,1147,106,927,362,90,398,2483]); A
            [5^3 * 2 + O(5^5)    2369 + O(5^5)    1147 + O(5^5)]
            [    106 + O(5^5)     927 + O(5^5)     362 + O(5^5)]
            [ 5 * 18 + O(5^5)     398 + O(5^5)    2483 + O(5^5)]
            sage: P = A._echelon_in_place('partial_pivoting'); A
            [1 + O(5^5)     O(5^5)     O(5^5)]
            [    O(5^5) 1 + O(5^5)     O(5^5)]
            [    O(5^5)     O(5^5) 1 + O(5^5)]

        The following example is an example where partial pivoting fails,
        but scaled partial pivoting succeeds, taken from 'Numerical Analysis (9th edition)'
        by R.L. Burden and J.D. Faires (with minor adjustments)::

            sage: RR13 = RealField(prec=13)
            sage: A = Matrix(RR13, 2, 3, [30, 591400, 591700, 5.291, -6.130, 46.78])
            sage: A
            [   30.0 591000. 592000.]
            [   5.29   -6.13    46.8]
            sage: P = A._echelon_in_place('classical'); A
            [ 1.00 0.000  12.0]
            [0.000  1.00  1.00]
            sage: A = Matrix(RR13, 2, 3, [30, 591400, 591700, 5.291, -6.130, 46.78])
            sage: P = A._echelon_in_place('partial_pivoting'); A
            [ 1.00 0.000  12.0]
            [0.000  1.00  1.00]
            sage: A = Matrix(RR13, 2, 3, [30, 591400, 591700, 5.291, -6.130, 46.78])
            sage: P = A._echelon_in_place('scaled_partial_pivoting'); A
            [ 1.00 0.000  10.0]
            [0.000  1.00  1.00]
        """
        cdef Py_ssize_t start_row, c, r, nr, nc, i, best_r
        if self.fetch('in_echelon_form'):
            return self.fetch('pivots')

        _ = verbose('generic in-place Gauss elimination on %s x %s matrix using %s algorithm' % (self._nrows, self._ncols, algorithm), level=2)
        self.check_mutability()
        cdef Matrix A

        nr = self._nrows
        nc = self._ncols
        A = self

        start_row = 0
        pivots = []
        cdef list scale_factors = []

        if algorithm == "classical":
            for c in range(nc):
                sig_check()
                for r in range(start_row, nr):
                    if A.get_unsafe(r, c):
                        pivots.append(c)
                        a_inverse = ~A.get_unsafe(r, c)
                        A.rescale_row(r, a_inverse, c)
                        A.swap_rows(r, start_row)
                        for i in range(nr):
                            if i != start_row:
                                if A.get_unsafe(i, c):
                                    minus_b = -A.get_unsafe(i, c)
                                    A.add_multiple_of_row(i, start_row, minus_b, c)
                        start_row = start_row + 1
                        break
        else:
            if algorithm == 'scaled_partial_pivoting':
                for r in range(nr):
                    scale_factor = 0
                    for c in range(nc):
                        abs_val = A.get_unsafe(r, c).abs()
                        if abs_val > scale_factor:
                            scale_factor = abs_val
                    scale_factors.append(scale_factor)
            elif algorithm == 'scaled_partial_pivoting_valuation':
                for r in range(nr):
                    scale_factor = None
                    for c in range(nc):
                        valuation = A.get_unsafe(r, c).valuation()
                        if scale_factor is None or valuation < scale_factor:
                            scale_factor = valuation
                    scale_factors.append(scale_factor)

            for c in range(nc):
                sig_check()
                if algorithm == 'scaled_partial_pivoting_valuation':
                    max_abs_val = None
                else:
                    max_abs_val = 0
                best_r = start_row - 1

                if algorithm == 'partial_pivoting':
                    for r in range(start_row, nr):
                        abs_val = A.get_unsafe(r, c).abs()
                        if abs_val > max_abs_val:
                            max_abs_val = abs_val
                            best_r = r
                elif algorithm == 'scaled_partial_pivoting':
                    for r in range(start_row, nr):
                        if scale_factors[r]:
                            abs_val = A.get_unsafe(r, c).abs() / scale_factors[r]
                            if abs_val > max_abs_val:
                                max_abs_val = abs_val
                                best_r = r
                else:  # algorithm == 'scaled_partial_pivoting_valuation':
                    for r in range(start_row, nr):
                        if scale_factors[r] is not None:
                            abs_val = scale_factors[r] - A.get_unsafe(r, c).valuation()
                            if max_abs_val is None or abs_val > max_abs_val:
                                max_abs_val = abs_val
                                best_r = r

                if max_abs_val or (algorithm == 'scaled_partial_pivoting_valuation' and max_abs_val is not None):
                    pivots.append(c)
                    a_inverse = ~A.get_unsafe(best_r, c)
                    A.rescale_row(best_r, a_inverse, c)
                    A.swap_rows(best_r, start_row)

                    if algorithm == 'scaled_partial_pivoting' or algorithm == 'scaled_partial_pivoting_valuation':
                        scale_factors[best_r], scale_factors[start_row] = scale_factors[start_row], scale_factors[best_r]

                    for i in range(nr):
                        if i != start_row:
                            if A.get_unsafe(i, c):
                                minus_b = -A.get_unsafe(i, c)
                                A.add_multiple_of_row(i, start_row, minus_b, c)
                    start_row = start_row + 1

        pivots = tuple(pivots)
        self.cache('pivots', pivots)
        self.cache('in_echelon_form', True)
        self.cache('echelon_form', self)
        return pivots

    # for backward compatibility

    def _echelon_in_place_classical(self):
        """
        Transform ``self`` into echelon form and return the pivots of ``self``.

        EXAMPLES::

            sage: t = matrix(QQ, 3, 3, range(9)); t
            [0 1 2]
            [3 4 5]
            [6 7 8]
            sage: E = t._echelon_in_place_classical(); t
            [ 1  0 -1]
            [ 0  1  2]
            [ 0  0  0]
            sage: a = matrix(QQ, 2, [1..6])
            sage: P = a._echelon_in_place_classical(); a
            [ 1  0 -1]
            [ 0  1  2]
        """
        return self._echelon_in_place('classical')

    def extended_echelon_form(self, subdivide=False, **kwds):
        r"""
        Return the echelon form of ``self`` augmented with an identity matrix.

        INPUT:

        - ``subdivide`` -- boolean (default: ``False``); whether to
          subdivide the returned matrix.  See the description of the
          (output) below for details.
        - ``kwds`` -- additional keywords that can be passed to
          the method that computes the echelon form

        OUTPUT:

        If `A` is an `m\times n` matrix, add the `m` columns of an
        `m\times m` identity matrix to the right of ``self``.  Then
        row-reduce this `m\times(n+m)` matrix.  This matrix is returned
        as an immutable matrix.

        If ``subdivide`` is ``True`` then the returned matrix has a single
        division among the columns and a single division among the rows.
        The column subdivision has `n` columns to the left and `m`
        columns to the right.  The row division separates the nonzero rows
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

        .. NOTE::

            This method returns an echelon form.  If the base ring
            is not a field, no attempt is made to move to the fraction field.
            See an example below where the base ring is changed manually.

        EXAMPLES:

        The four relationships at the end of this example hold in general. ::

            sage: A = matrix(QQ, [[2, -1, 7, -1, 0, -3],
            ....:                 [-1, 1, -5, 3, 4, 4],
            ....:                 [2, -1, 7, 0, 2, -2],
            ....:                 [2, 0, 4, 3, 6, 1],
            ....:                 [2, -1, 7, 0, 2, -2]])
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
            ....:                 [10, 13, 3, -3, 7, 11, 6],
            ....:                 [-12, -15, -3, -3, 2, -10, -4],
            ....:                 [10, 13, 3, 3, -1, 9, 4],
            ....:                 [4, 5, 1, 8, -10, 1, -1]])
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
            sage: G.subdivide([], []); G
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
        """
        if subdivide not in [True, False]:
            raise TypeError("subdivide must be True or False, not %s" % subdivide)
        ident = self.matrix_space(self.nrows(), self.nrows()).one()
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

    #####################################################################################
    # Functions for symmetries of a matrix under row and column permutations
    #####################################################################################
    def as_bipartite_graph(self):
        r"""
        Construct a bipartite graph ``B`` representing
        the matrix uniquely.

        Vertices are labeled 1 to ``nrows``
        on the left and ``nrows`` + 1 to ``nrows`` + ``ncols`` on
        the right, representing rows and columns correspondingly.
        Each row is connected to each column with an edge
        weighted by the value of the corresponding matrix entry.

        This graph is a helper for calculating automorphisms of
        a matrix under row and column permutations. See
        :meth:`automorphisms_of_rows_and_columns`.

        OUTPUT: a bipartite graph

        EXAMPLES::

            sage: M = matrix(QQ, [[1/3, 7], [6, 1/4], [8, -5]])
            sage: M
            [1/3   7]
            [  6 1/4]
            [  8  -5]

            sage: # needs sage.graphs
            sage: B = M.as_bipartite_graph(); B
            Bipartite graph on 5 vertices
            sage: B.edges(sort=True)
            [(1, 4, 1/3), (1, 5, 7), (2, 4, 6), (2, 5, 1/4), (3, 4, 8), (3, 5, -5)]
            sage: len(B.left) == M.nrows()
            True
            sage: len(B.right) == M.ncols()
            True
        """
        from sage.graphs.bipartite_graph import BipartiteGraph
        nrows = self.nrows()
        ncols = self.ncols()
        B = BipartiteGraph()
        B.add_vertices(range(1, nrows + 1), left=True)
        B.add_vertices(range(nrows + 1, nrows + ncols + 1), right=True)
        for i in range(nrows):
            for j in range(ncols):
                B.add_edge(i + 1, nrows + 1 + j, self[i][j])
        return B

    def automorphisms_of_rows_and_columns(self):
        r"""
        Return the automorphisms of ``self`` under
        permutations of rows and columns as a list of
        pairs of ``PermutationGroupElement`` objects.

        EXAMPLES::

            sage: # needs sage.graphs sage.groups
            sage: M = matrix(ZZ,[[1,0],[1,0],[0,1]]); M
            [1 0]
            [1 0]
            [0 1]
            sage: A = M.automorphisms_of_rows_and_columns(); A
            [((), ()), ((1,2), ())]
            sage: M = matrix(ZZ, [[1,1,1,1],[1,1,1,1]])
            sage: A = M.automorphisms_of_rows_and_columns()
            sage: len(A)
            48

        One can now apply these automorphisms to ``M`` to show
        that it leaves it invariant::

            sage: all(M.with_permuted_rows_and_columns(*i) == M for i in A)             # needs sage.graphs sage.groups
            True

        Check that :issue:`25426` is fixed::

            sage: j = matrix([(3, 2, 1, 0, 0),
            ....:             (2, 2, 0, 1, 0),
            ....:             (1, 0, 3, 0, 2),
            ....:             (0, 1, 0, 2, 1),
            ....:             (0, 0, 2, 1, 2)])
            sage: j.automorphisms_of_rows_and_columns()                                 # needs sage.graphs sage.groups
            [((), ()), ((1,3)(2,5), (1,3)(2,5))]
        """
        from sage.groups.perm_gps.constructor import \
            PermutationGroupElement
        B = self.as_bipartite_graph()
        nrows = self.nrows()
        ncols = self.ncols()
        A = B.automorphism_group(edge_labels=True)

        # Convert to elements of Sym(self) from S_n
        permutations = []
        for p in A:
            if p(1) <= nrows:  # Check that rows are mapped to rows
                permutations.append(
                    (PermutationGroupElement([p(1 + i) for i in range(nrows)]),
                     PermutationGroupElement([p(1 + nrows + i) - nrows for i in range(ncols)])
                     ))
        return permutations

    def permutation_normal_form(self, check=False):
        r"""
        Take the set of matrices that are ``self``
        permuted by any row and column permutation,
        and return the maximal one of the set where
        matrices are ordered lexicographically
        going along each row.

        INPUT:

        - ``check`` -- boolean (default: ``False``); if ``True`` return a tuple of
          the maximal matrix and the permutations taking ``self``
          to the maximal matrix.
          If ``False``, return only the maximal matrix.

        OUTPUT: the maximal matrix

        EXAMPLES::

            sage: M = matrix(ZZ, [[0, 0, 1], [1, 0, 2], [0, 0, 0]])
            sage: M
            [0 0 1]
            [1 0 2]
            [0 0 0]

            sage: M.permutation_normal_form()
            [2 1 0]
            [1 0 0]
            [0 0 0]

            sage: M = matrix(ZZ, [[-1, 3], [-1, 5], [2, 4]])
            sage: M
            [-1  3]
            [-1  5]
            [ 2  4]

            sage: M.permutation_normal_form(check=True)                                 # needs sage.graphs sage.groups
            (
            [ 5 -1]
            [ 4  2]
            [ 3 -1],
            ((1,2,3), (1,2))
            )

        TESTS::

            sage: # needs sage.graphs
            sage: M = matrix(ZZ, [[3, 4, 5], [3, 4, 5], [3, 5, 4], [2, 0,1]])
            sage: M.permutation_normal_form()
            [5 4 3]
            [5 4 3]
            [4 5 3]
            [1 0 2]
            sage: M = matrix(ZZ, 0, 0, [])
            sage: M.permutation_normal_form()
            []
            sage: M.permutation_normal_form(check=True)
            ([], ((), ()))
        """
        nrows = self.nrows()
        ncols = self.ncols()
        if not nrows or not ncols:
            if not check:
                return self
            from sage.groups.perm_gps.permgroup_named import SymmetricGroup
            return (self, (SymmetricGroup(nrows).one(),
                           SymmetricGroup(ncols).one()))

        # A helper
        def new_sorted_matrix(m):
            return self.new_matrix(
                ncols, nrows,
                sorted(m.columns(), reverse=True)).transpose()

        # Let us sort each row:
        sorted_rows = [sorted([self[i][j] for j in range(ncols)], reverse=True)
                       for i in range(nrows)]
        # and find the maximal one:
        first_row = max(sorted_rows)
        first_rows = [j for j in range(nrows) if sorted_rows[j] == first_row]
        # We construct an array S, which will record the subsymmetries of the
        # columns, i.e. S records the automorphisms with respect to the column
        # swappings of the upper block already constructed. For example, if
        # allowed, and so on. S is in decreasing order and takes values between
        # S:=[a, a, ..., a (i-th), a-1, a-1, ...] then any swap between 1 and i is
        # me + nc and me + 1 such that no entry is present already in the matrix.
        S = [first_row[0] + ncols] + [None]*(ncols-1)
        for j in range(1, ncols):
            S[j] = S[j - 1] if first_row[j] == first_row[j - 1] else S[j - 1] - 1
        # If we want to sort the i-th row with respect to a symmetry determined
        # by S, then we will just sort the augmented row [[S[j], PM[i, j]] :
        # j in [1 .. nc]], S having lexicographic priority.

        # MS is a list of non-isomorphic matrices where one of the maximal rows
        # has been replaced by S
        MS = []
        aM = []
        if self.is_immutable():
            X = copy(self)
        else:
            X = self
        for i in first_rows:
            N = new_sorted_matrix(X.with_swapped_rows(0, i))
            aN = copy(N)
            aN.set_row(0, S)
            if not any(aN.is_permutation_of(j) for j in aM):
                MS.append(N)
                aM.append(aN)
        # We construct line l:
        for l in range(1, nrows - 1):
            if not S == list(range(first_row[0] + ncols, first_row[0], -1)):
                # Sort each row with respect to S for the first matrix in X = MS
                X = copy(MS)
                SM = [sorted([(S[j], X[0][k][j]) for j in range(ncols)], reverse=True)
                               for k in range(l, nrows)]
                SM = [[k[1] for k in s] for s in SM]

                # and pick the maximal row
                b = max(SM)
                # Find all rows equal to the maximal (potential new cases)
                m = [[j for j in range(nrows - l) if SM[j] == b]]
                w = 0  # keeps track of how many entries we have removed from MS
                # Let us find the maximal row in each of the entries in X = MS
                for i in range(1, len(X)):
                    SN = [sorted([(S[j], X[i][k][j]) for j in range(ncols)], reverse=True)
                          for k in range(l, nrows)]
                    SN = [[k[1] for k in s] for s in SN]
                    # Maximal row in this entry of X = MS
                    nb = max(SN)
                    # Number of occurrences.
                    n = [j for j in range(nrows - l) if SN[j] == nb]
                    # Now compare to our previous max
                    if b < nb:
                        # Bigger so save line
                        b = nb
                        m = [n]
                        # and delete all previous attempts
                        u = i - w
                        del MS[0:u]
                        w += u
                    elif b == nb:
                        # Same, so save symmetry
                        m.append(n)
                    else:
                        # Smaller, so forget about it!
                        MS.pop(i - w)
                        w += 1
                # Update symmetries
                test = [(S[i], b[i]) for i in range(ncols)]
                for j in range(1, ncols):
                    S[j] = S[j - 1] if (test[j] == test[j - 1]) else S[j - 1] - 1
                # For each case we check the isomorphism as previously, if
                # test fails we add a new non-isomorphic case to MS. We
                # pick our choice of maximal line l and sort the columns of
                # the matrix (this preserves symmetry automatically).
                n = len(MS)
                for i in range(n):
                    MS[i] = new_sorted_matrix(MS[i].with_swapped_rows(l, m[i][0] + l))
                    if len(m[i]) > 1:
                        aX = MS[i].submatrix(l, 0, nrows - l, ncols)
                        aX.set_row(0, S)
                        aX = [aX]
                        for j in m[i][1:]:
                            N = new_sorted_matrix(MS[i].with_swapped_rows(l, j + 1))
                            aN = N.submatrix(l, 0, nrows - l, ncols)
                            aN.set_row(0, S)
                            if not any(aN.is_permutation_of(k) for k in aX):
                                MS.append(N)
                                aX.append(aN)
            else:
                MS = [self.new_matrix(nrows, ncols, sorted(s.rows(), reverse=True)) for s in MS]
                break
        MS_max = max(MS)
        if check:
            return MS_max, self.is_permutation_of(MS_max, True)[1]
        else:
            return MS_max

    def is_permutation_of(self, N, check=False):
        r"""
        Return ``True`` if there exists a permutation of rows and
        columns sending ``self`` to ``N`` and ``False`` otherwise.

        INPUT:

        - ``N`` -- a matrix

        - ``check`` -- boolean (default: ``False``); if ``False``
          return Boolean indicating whether there exists a permutation of
          rows and columns sending ``self`` to ``N`` and ``False`` otherwise.
          If ``True`` return a tuple of a Boolean and a permutation mapping
          ``self`` to ``N`` if such a permutation exists, and
          (``False``, ``None``) if it does not.

        OUTPUT: a Boolean or a tuple of a Boolean and a permutation

        EXAMPLES::

            sage: M = matrix(ZZ, [[1,2,3], [3,5,3], [2,6,4]])
            sage: M
            [1 2 3]
            [3 5 3]
            [2 6 4]
            sage: N = matrix(ZZ, [[1,2,3], [2,6,4], [3,5,3]])
            sage: N
            [1 2 3]
            [2 6 4]
            [3 5 3]
            sage: M.is_permutation_of(N)                                                # needs sage.graphs
            True

        Some examples that are not permutations of each other::

            sage: N = matrix(ZZ, [[1,2,3], [4,5,6], [7,8,9]]); N
            [1 2 3]
            [4 5 6]
            [7 8 9]
            sage: M.is_permutation_of(N)                                                # needs sage.graphs
            False
            sage: N = matrix(ZZ, [[1,2], [3,4]]); N
            [1 2]
            [3 4]
            sage: M.is_permutation_of(N)
            False

        And for when ``check`` is True::

            sage: # needs sage.graphs sage.groups
            sage: N = matrix(ZZ, [[3,5,3], [2,6,4], [1,2,3]]); N
            [3 5 3]
            [2 6 4]
            [1 2 3]
            sage: r = M.is_permutation_of(N, check=True)
            sage: r
            (True, ((1,2,3), ()))
            sage: p = r[1]
            sage: M.with_permuted_rows_and_columns(*p) == N
            True
        """
        ncols = self.ncols()
        nrows = self.nrows()
        if N.ncols() != ncols or N.nrows() != nrows:
            if check:
                return (False, None)
            else:
                return False
        M_B = self.as_bipartite_graph()
        N_B = N.as_bipartite_graph()
        if check:
            truth, perm = N_B.is_isomorphic(M_B, certificate=True, edge_labels=True)
            from sage.groups.perm_gps.constructor import PermutationGroupElement
            if perm:
                s = sorted(perm.items(), key=lambda x:x[0])
                row_perms = [value for k, value in s if k <= nrows]
                col_perms = [value - nrows for k, value in s if k > nrows]
                perm = (PermutationGroupElement(row_perms), PermutationGroupElement(col_perms))
            return truth, perm
        else:
            return N_B.is_isomorphic(M_B, certificate=False, edge_labels=True)

    #####################################################################################
    # Windowed Strassen Matrix Multiplication and Echelon
    # Precise algorithms invented and implemented by David Harvey and Robert Bradshaw
    # at William Stein's MSRI 2006 Summer Workshop on Modular Forms.
    #####################################################################################
    def _multiply_strassen(self, Matrix right, int cutoff=0):
        """
        Multiply ``self`` by the matrix right using a Strassen-based
        asymptotically fast arithmetic algorithm.

        ALGORITHM: Custom algorithm for arbitrary size matrices designed by
        David Harvey and Robert Bradshaw, based on Strassen's algorithm.

        INPUT:

        - ``cutoff`` -- integer (default: 0 -- let class
          decide)

        EXAMPLES::

            sage: a = matrix(ZZ, 4,4, range(16))
            sage: a._multiply_strassen(a,2)
            [ 56  62  68  74]
            [152 174 196 218]
            [248 286 324 362]
            [344 398 452 506]
        """
        if self._ncols != right._nrows:
            raise ArithmeticError("Number of columns of self must equal number of rows of right.")
        if self._base_ring is not right.base_ring():
            raise TypeError("Base rings must be the same.")

        if cutoff == 0:
            cutoff = self._strassen_default_cutoff(right)

        if cutoff <= 0:
            raise ValueError("cutoff must be at least 1")

        output = self.new_matrix(self._nrows, right._ncols)
        # The following used to be a little faster, but meanwhile
        # the previous line is faster.
        # if self.is_sparse():
        #    output = self.matrix_space(self._nrows, right._ncols, sparse = True)(0)
        # else:
        #    output = self.matrix_space(self._nrows, right._ncols, sparse = False).zero_matrix().__copy__()

        self_window = self.matrix_window()
        right_window = right.matrix_window()
        output_window = output.matrix_window()

        from sage.matrix import strassen
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
        tm = verbose('strassen echelon of %s x %s matrix'%(self._nrows, self._ncols), level=2)

        self.check_mutability()

        if self._base_ring not in _Fields:
            raise ValueError("Echelon form not defined over this base ring.")

        if cutoff == 0:
            cutoff = self._strassen_default_echelon_cutoff()

        if cutoff < 1:
            raise ValueError("cutoff must be at least 1")

        if self._nrows < cutoff or self._ncols < cutoff:
            self._echelon_in_place_classical()
            return

        from sage.matrix import strassen
        pivots = strassen.strassen_echelon(self.matrix_window(), cutoff)
        self.cache('pivots', pivots)
        verbose('done with strassen', level=2, t=tm)

    cpdef matrix_window(self, Py_ssize_t row=0, Py_ssize_t col=0,
                        Py_ssize_t nrows=-1, Py_ssize_t ncols=-1,
                        bint check=1):
        """
        Return the requested matrix window.

        EXAMPLES::

            sage: A = matrix(QQ, 3, 3, range(9))
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
        from sage.matrix import matrix_window
        if nrows == -1:
            nrows = self._nrows - row
        if ncols == -1:
            ncols = self._ncols - col
        if check and (row < 0 or col < 0 or row + nrows > self._nrows or
                      col + ncols > self._ncols):
            raise IndexError("matrix window index out of range")
        return matrix_window.MatrixWindow(self, row, col, nrows, ncols)

    def set_block(self, row, col, block):
        """
        Set the sub-matrix of self, with upper left corner given by row,
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
        r"""
        Divide ``self`` into logical submatrices which can then be queried
        and extracted.

        If a subdivision already exists, this method forgets the
        previous subdivision and flushes the cache.

        INPUT:

        - ``row_lines`` -- ``None``, an integer, or a list of
          integers (lines at which ``self`` must be split)

        - ``col_lines`` -- ``None``, an integer, or a list of
          integers (columns at which ``self`` must be split)

        OUTPUT: none but changes ``self``

        .. NOTE::

           One may also pass a tuple into the first argument which
           will be interpreted as ``(row_lines, col_lines)``.

        EXAMPLES::

            sage: # needs sage.libs.pari
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

        Degenerate cases work too::

            sage: # needs sage.libs.pari
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

        Indices do not need to be in the right order (:issue:`14064`)::

            sage: M.subdivide([4, 2], [3, 1]); M                                        # needs sage.libs.pari
            [ 2| 3  5| 7 11]
            [13|17 19|23 29]
            [--+-----+-----]
            [31|37 41|43 47]
            [53|59 61|67 71]
            [--+-----+-----]
            [73|79 83|89 97]

        TESTS:

        Input such that the matrix has no subdivision results in
        the ``_subdivision`` attribute being set to ``None``::

            sage: A = matrix.identity(QQ, 4)
            sage: A._subdivisions is None
            True
            sage: A.subdivide()
            sage: A._subdivisions is None
            True
            sage: A.subdivide(2, 3)  # perform a subdivision
            sage: A._subdivisions is None
            False
            sage: A.subdivide(([], []))  # now reset
            sage: A._subdivisions is None
            True
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
        if self._subdivisions is not None:
            self.clear_cache()
        if (not row_lines) and (not col_lines):
            self._subdivisions = None
        else:
            l_row = sorted(row_lines)
            l_col = sorted(col_lines)
            l_row = [0] + [int(ZZ(x)) for x in l_row] + [self._nrows]
            l_col = [0] + [int(ZZ(x)) for x in l_col] + [self._ncols]
            self._subdivisions = (l_row, l_col)

    def subdivision(self, i, j):
        """
        Return an immutable copy of the (i,j)th submatrix of ``self``,
        according to a previously set subdivision.

        Before a subdivision is set, the only valid arguments are (0,0)
        which returns ``self``.

        EXAMPLES::

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
        key = "subdivision %s %s"%(i, j)
        sd = self.fetch(key)
        if sd is None:
            sd = self[self._subdivisions[0][i]:self._subdivisions[0][i+1],
                      self._subdivisions[1][j]:self._subdivisions[1][j+1]]
            sd.set_immutable()
            self.cache(key, sd)
        return sd

    def subdivision_entry(self, i, j, x, y):
        """
        Return the x,y entry of the i,j submatrix of ``self``.

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
                return self[x, y]
            else:
                raise IndexError("No such submatrix %s, %s"%(i, j))
        if x >= self._subdivisions[0][i+1]-self._subdivisions[0][i] or \
           y >= self._subdivisions[1][j+1]-self._subdivisions[1][j]:
            raise IndexError("Submatrix %s,%s has no entry %s,%s"%(i, j, x, y))
        return self[self._subdivisions[0][i] + x, self._subdivisions[1][j] + y]

    def _subdivide_on_augment(self, left, right):
        r"""
        Helper method to manage subdivisions when augmenting a matrix.

        INPUT:

        - ``left``, ``right`` -- two matrices, such that if ``left`` is
          augmented by placing ``right`` on the right side of ``left``,
          then the result is ``self``.  It is the responsibility of the
          calling routine to ensure this condition holds.

        OUTPUT:

        ``None``.  A new subdivision is created between ``left`` and
        ``right`` for ``self``.  If possible, row subdivisions are
        preserved in ``self``, but if the two sets of row subdivisions
        are incompatible, they are removed.

        EXAMPLES::

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

        - ``top``, ``bottom`` -- two matrices, such that if ``top`` is
        stacked by placing ``top`` above ``bottom``, then the result
        is ``self``.  It is the responsibility of the calling routine
        to ensure this condition holds.

        OUTPUT:

        ``None``.  A new subdivision is created between ``top`` and
        ``bottom`` for ``self``.  If possible, column subdivisions are
        preserved in ``self``, but if the two sets of column subdivisions
        are incompatible, they are removed.

        EXAMPLES::

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
        Return the current subdivision of ``self``.

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

    # for backwards compatibility: see #4983.
    get_subdivisions = subdivisions

    def tensor_product(self, A, subdivide=True):
        r"""
        Return the tensor product of two matrices.

        INPUT:

        - ``A`` -- a matrix
        - ``subdivide`` -- (default: ``True``) whether or not to return
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
            TypeError: unsupported operand parent(s) for *:
             'Finite Field of size 29' and
             'Full MatrixSpace of 3 by 4 dense matrices over Finite Field of size 23'

        The input is checked to be sure it is a matrix.  ::

            sage: A = matrix(QQ, 2, 2, range(4))
            sage: A.tensor_product('junk')
            Traceback (most recent call last):
            ...
            TypeError: tensor product requires a second matrix, not junk

        TESTS:

        Check that `m \times 0` and `0 \times m` matrices work
        (:issue:`22769`)::

            sage: m1 = matrix(QQ, 1, 0, [])
            sage: m2 = matrix(QQ, 2, 2, [1, 2, 3, 4])
            sage: m1.tensor_product(m2).dimensions()
            (2, 0)
            sage: m2.tensor_product(m1).dimensions()
            (2, 0)
            sage: m3 = matrix(QQ, 0, 3, [])
            sage: m3.tensor_product(m2).dimensions()
            (0, 6)
            sage: m2.tensor_product(m3).dimensions()
            (0, 6)

            sage: m1 = MatrixSpace(GF(5), 3, 2).an_element()
            sage: m2 = MatrixSpace(GF(5), 0, 4).an_element()
            sage: m1.tensor_product(m2).parent()
            Full MatrixSpace of 0 by 8 dense matrices over Finite Field of size 5
        """
        if not isinstance(A, Matrix):
            raise TypeError('tensor product requires a second matrix, not {0}'.format(A))
        from sage.matrix.constructor import block_matrix
        # Special case when one of the matrices is 0 \times m or m \times 0
        if self.nrows() == 0 or self.ncols() == 0 or A.nrows() == 0 or A.ncols() == 0:
            return self.matrix_space(self.nrows()*A.nrows(),
                                     self.ncols()*A.ncols()).zero_matrix().__copy__()
        return block_matrix(self.nrows(), self.ncols(),
                            [x * A for x in self.list()], subdivide=subdivide)

    def randomize(self, density=1, nonzero=False, *args, **kwds):
        """
        Replace a proportion of the entries of a matrix by random elements,
        leaving the remaining entries unchanged.

        .. NOTE::

           The locations of the entries of the matrix to change are
           determined randomly, with the total number of locations
           determined by the ``density`` keyword. These locations
           are not guaranteed to be distinct.  So it is possible
           that the same position can be chosen multiple times,
           especially for a very small matrix.  The exception is
           when ``density = 1``, in which case every entry of the
           matrix will be changed.

        INPUT:

        - ``density`` -- ``float`` (default: ``1``); upper bound for the
          proportion of entries that are changed
        - ``nonzero`` -- boolean (default: ``False``); if ``True``, then new
          entries will be nonzero
        - ``*args, **kwds`` -- remaining parameters may be passed to the
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
            sage: a.density() <= 0.5
            True

        Now we randomize all the entries of the resulting matrix::

            sage: while a.density() < 0.9:
            ....:     a = matrix(QQ['x'], 3)
            ....:     a.randomize()

        We create the zero matrix over the integers::

            sage: a = matrix(ZZ, 2); a
            [0 0]
            [0 0]

        Then we randomize it; the ``x`` and ``y`` keywords, which determine the
        size of the random elements, are passed on to the ``random_element``
        method for ``ZZ``.

        ::

            sage: a.randomize(x=-2^64, y=2^64)
            sage: while all(abs(b) < 2^63 for b in a.list()):
            ....:     a.randomize(x=-2^64, y=2^64)
            sage: all(abs(b) < 2^64 for b in a.list())
            True
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
                        self.set_unsafe(i, j, R._random_nonzero_element(*args, **kwds))
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
        Return ``True`` if this matrix is the identity matrix.

        EXAMPLES::

            sage: m = matrix(QQ, 2,2, range(4))
            sage: m.is_one()
            False
            sage: m = matrix(QQ, 2, [5,0,0,5])
            sage: m.is_one()
            False
            sage: m = matrix(QQ, 2, [1,0,0,1])
            sage: m.is_one()
            True
            sage: m = matrix(QQ, 2, [1,1,1,1])
            sage: m.is_one()
            False
        """
        return self.is_scalar(self.base_ring().one())

    def is_scalar(self, a=None):
        """
        Return ``True`` if this matrix is a scalar matrix.

        INPUT:

        - ``a`` -- base_ring element; chosen as ``self[0][0]`` if
          ``a==None``

        OUTPUT:

        Whether ``self`` is a scalar matrix (in fact the scalar matrix
        aI if a is input).

        EXAMPLES::

            sage: m = matrix(QQ, 2,2, range(4))
            sage: m.is_scalar(5)
            False
            sage: m = matrix(QQ, 2, [5,0,0,5])
            sage: m.is_scalar(5)
            True
            sage: m = matrix(QQ, 2, [1,0,0,1])
            sage: m.is_scalar(1)
            True
            sage: m = matrix(QQ, 2, [1,1,1,1])
            sage: m.is_scalar(1)
            False
        """
        if not self.is_square():
            return False
        cdef Py_ssize_t i, j
        if a is None:
            if self._nrows == 0:
                return True
            a = self.get_unsafe(0, 0)
        else:
            a = self.base_ring()(a)
        for i in range(self._nrows):
            for j in range(self._ncols):
                if i != j:
                    if not self.get_unsafe(i, j).is_zero():
                        return False
                else:
                    if self.get_unsafe(i, i) != a:
                        return False
        return True

    def is_diagonal(self) -> bool:
        """
        Return ``True`` if this matrix is a diagonal matrix.

        OUTPUT: boolean

        EXAMPLES::

            sage: m = matrix(QQ, 2,2, range(4))
            sage: m.is_diagonal()
            False
            sage: m = matrix(QQ, 2, [5,0,0,5])
            sage: m.is_diagonal()
            True
            sage: m = matrix(QQ, 2, [1,0,0,1])
            sage: m.is_diagonal()
            True
            sage: m = matrix(QQ, 2, [1,1,1,1])
            sage: m.is_diagonal()
            False
        """
        if not self.is_square():
            return False
        cdef Py_ssize_t i, j

        for i in range(self._nrows):
            for j in range(self._ncols):
                if i != j:
                    if not self.get_unsafe(i, j).is_zero():
                        return False
        return True

    def is_triangular(self, side='lower') -> bool:
        """
        Return ``True`` if this matrix is a triangular matrix.

        INPUT:

        - ``side`` -- either ``'lower'`` (default) or ``'upper'``

        OUTPUT: boolean

        EXAMPLES::

            sage: m = matrix(QQ, 2, 2, range(4))
            sage: m.is_triangular()
            False
            sage: m = matrix(QQ, 2, [5, 0, 0, 5])
            sage: m.is_triangular()
            True
            sage: m = matrix(QQ, 2, [1, 2, 0, 1])
            sage: m.is_triangular("upper")
            True
            sage: m.is_triangular("lower")
            False
        """
        if not self.is_square():
            return False
        cdef Py_ssize_t i, j

        if side == "upper":
            for i in range(1, self._nrows):
                for j in range(i):
                    if not self.get_unsafe(i, j).is_zero():
                        return False
        else:
            for i in range(self._nrows - 1):
                for j in range(i + 1, self._ncols):
                    if not self.get_unsafe(i, j).is_zero():
                        return False
        return True

    def is_unitary(self) -> bool:
        r"""
        Return ``True`` if the columns of the matrix are an orthonormal basis.

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

            sage: A = matrix(QQbar,                                                     # needs sage.rings.number_field sage.symbolic
            ....:            [[(1/sqrt(5))*(1+i), (1/sqrt(55))*(3+2*I), (1/sqrt(22))*(2+2*I)],
            ....:             [(1/sqrt(5))*(1-i), (1/sqrt(55))*(2+2*I),  (1/sqrt(22))*(-3+I)],
            ....:             [    (1/sqrt(5))*I, (1/sqrt(55))*(3-5*I),    (1/sqrt(22))*(-2)]])
            sage: A.is_unitary()                                                        # needs sage.rings.number_field sage.symbolic
            True

        A permutation matrix is always orthogonal. ::

            sage: # needs sage.combinat
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

            sage: A = matrix(QQbar, 3, 4)                                               # needs sage.rings.number_field
            sage: A.is_unitary()
            False
        """
        if not self.is_square():
            return False
        if hasattr(self.base_ring().an_element(), 'conjugate'):
            P = self.conjugate().transpose() * self    # Unitary
        else:
            P = self.transpose() * self                # Orthogonal
        return P.is_scalar(1)

    def is_bistochastic(self, normalized=True):
        r"""
        Return ``True`` if this matrix is bistochastic.

        A matrix is said to be bistochastic if both the sums of the
        entries of each row and the sum of the entries of each column
        are equal to 1 and all entries are nonnegative.

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
            sage: (2 * Matrix(5,5,1)).is_bistochastic(normalized=False)
            True

        Here is a matrix whose row and column sums is 1, but not all entries are
        nonnegative::

            sage: m = matrix([[-1,2],[2,-1]])
            sage: m.is_bistochastic()
            False
        """

        row_sums = [sum(r) for r in self.rows()]
        col_sums = [sum(c) for c in self.columns()]

        return self.is_square() and\
               col_sums[0] == row_sums[0] and\
               row_sums == col_sums and\
               row_sums == len(row_sums) * [col_sums[0]] and\
               ((not normalized) or col_sums[0] == self.base_ring()(1)) and\
               all(entry >= 0 for row in self for entry in row)

    def is_normal(self):
        r"""
        Return ``True`` if the matrix commutes with its conjugate-transpose.

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

            sage: # needs sage.symbolic
            sage: A = matrix(QQ, 5, 5, range(25)) + I*matrix(QQ, 5, 5, range(0, 50, 2))
            sage: B = A*A.conjugate_transpose()
            sage: B.is_hermitian()
            True
            sage: B.is_normal()
            True

        Circulant matrices are normal.  ::

            sage: # needs sage.graphs
            sage: G = graphs.CirculantGraph(20, [3, 7])
            sage: D = digraphs.Circuit(20)
            sage: A = 3*D.adjacency_matrix() - 5*G.adjacency_matrix()
            sage: A.is_normal()
            True

        Skew-symmetric matrices are normal.  ::

            sage: A = matrix(QQ, 5, 5, range(25))
            sage: B = A - A.transpose()
            sage: B.is_skew_symmetric()
            True
            sage: B.is_normal()
            True

        A small matrix that does not fit into any of the usual categories
        of normal matrices.  ::

            sage: A = matrix(ZZ, [[1, -1],
            ....:                 [1,  1]])
            sage: A.is_normal()
            True
            sage: not A.is_hermitian() and not A.is_skew_symmetric()
            True

        Sage has several fields besides the entire complex numbers
        where conjugation is non-trivial. ::

            sage: # needs sage.rings.number_field
            sage: F.<b> = QuadraticField(-7)
            sage: C = matrix(F, [[-2*b - 3,  7*b - 6, -b + 3],
            ....:                [-2*b - 3, -3*b + 2,   -2*b],
            ....:                [   b + 1,        0,     -2]])
            sage: C = C*C.conjugate_transpose()
            sage: C.is_normal()
            True

        A matrix that is nearly normal, but for a non-real
        diagonal entry. ::

            sage: # needs sage.rings.number_field
            sage: A = matrix(QQbar, [[    2,   2-I, 1+4*I],
            ....:                    [  2+I,   3+I, 2-6*I],
            ....:                    [1-4*I, 2+6*I,     5]])
            sage: A.is_normal()
            False
            sage: A[1,1] = 132
            sage: A.is_normal()
            True

        Rectangular matrices are never normal.  ::

            sage: A = matrix(QQbar, 3, 4)                                               # needs sage.rings.number_field
            sage: A.is_normal()                                                         # needs sage.rings.number_field
            False

        A square, empty matrix is trivially normal.  ::

            sage: A = matrix(QQ, 0, 0)
            sage: A.is_normal()
            True
        """
        key = 'normal'
        n = self.fetch(key)
        if n is not None:
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

        cdef Py_ssize_t i, j
        normal = True
        # two products are Hermitian, need only check lower triangle
        for i from 0 <= i < self._nrows:
            for j from 0 <= j <= i:
                if left.get_unsafe(i, j) != right.get_unsafe(i, j):
                    normal = False
                    break
            if not normal:
                break
        self.cache(key, normal)
        return normal

    def is_nilpotent(self) -> bool:
        r"""
        Return if ``self`` is a nilpotent matrix.

        A matrix `A` is *nilpotent* if there exists a positive integer
        `k` such that `A^k = 0`. We test this by using the Cayley-Hamilton
        theorem to see if the characteristic polynomial is `x^n = 0`.

        EXAMPLES::

            sage: A = matrix([[0,2,1,6], [0,0,1,2], [0,0,0,3], [0,0,0,0]])
            sage: A.is_nilpotent()
            True
            sage: B = matrix([[2,2,2,2,-4], [7,1,1,1,-5], [1,7,1,1,-5],
            ....:             [1,1,7,1,-5], [1,1,1,7,-5]])
            sage: B.is_nilpotent()
            True
            sage: C = matrix(GF(7), [[1, 2], [2, 6]])
            sage: C.is_nilpotent()
            False
            sage: D = matrix([[1,0],[0,0]])
            sage: D.is_nilpotent()
            False
            sage: Z = matrix.zero(QQ, 5)
            sage: Z.is_nilpotent()
            True

        TESTS:

        Check the corner case of the `0 \times 0` matrix::

            sage: Z = matrix.zero(QQ, 0); Z
            []
            sage: Z.charpoly()
            1
            sage: Z.is_nilpotent()
            True
        """
        if self.trace():
            return False
        chi = self.charpoly()
        return chi.is_monomial()

    def is_semisimple(self) -> bool:
        r"""
        Return if ``self`` is semisimple.

        A (square) matrix `A` is *semisimple* if the
        :meth:`minimal polynomial <minpoly>` of `A` is sqaure-free.

        If `A` represents a linear map from `F^n \to F^n` for some field `F`,
        then this is equivalent to every `A`-invariant subspace of `F^n`
        has a complementary `A`-invariant subspace. This is also equivalent
        to saying the matrix is diagonalizable over `\bar{F}`, the algebraic
        closure of `F`.

        EXAMPLES::

            sage: A = matrix([[0, -1], [1, 0]]); A
            [ 0 -1]
            [ 1  0]
            sage: A.is_semisimple()
            True
            sage: A.change_ring(QQ).is_diagonalizable()
            False
            sage: A.change_ring(CyclotomicField(4)).is_diagonalizable()
            True
        """
        return self.minpoly().is_squarefree()

    def as_sum_of_permutations(self):
        r"""
        Return the current matrix as a sum of permutation matrices.

        According to the Birkhoff-von Neumann Theorem, any bistochastic matrix
        can be written as a positive sum of permutation matrices, which also
        means that the polytope of bistochastic matrices is integer.

        As a non-bistochastic matrix can obviously not be written as a sum of
        permutations, this theorem is an equivalence.

        This function, given a bistochastic matrix, returns the corresponding
        decomposition.

        .. SEEALSO::

            - :func:`bistochastic_as_sum_of_permutations
              <sage.combinat.permutation.bistochastic_as_sum_of_permutations>`
              -- for more information on this method.

            - :meth:`~sage.geometry.polyhedron.library.Polytopes.Birkhoff_polytope`

        EXAMPLES:

        We create a bistochastic matrix from a convex sum of permutations, then
        try to deduce the decomposition from the matrix ::

            sage: L = []
            sage: L.append((9, Permutation([4, 1, 3, 5, 2])))
            sage: L.append((6, Permutation([5, 3, 4, 1, 2])))
            sage: L.append((3, Permutation([3, 1, 4, 2, 5])))
            sage: L.append((2, Permutation([1, 4, 2, 3, 5])))
            sage: M = sum([c * p.to_matrix() for c, p in L])
            sage: from sage.combinat.permutation import bistochastic_as_sum_of_permutations
            sage: decomp = bistochastic_as_sum_of_permutations(M)                       # needs sage.combinat sage.graphs
            sage: print(decomp)                                                         # needs sage.combinat sage.graphs
            2*B[[1, 4, 2, 3, 5]] + 3*B[[3, 1, 4, 2, 5]] + 9*B[[4, 1, 3, 5, 2]] + 6*B[[5, 3, 4, 1, 2]]

        An exception is raised when the matrix is not bistochastic::

            sage: M = Matrix([[2,3],[2,2]])
            sage: decomp = bistochastic_as_sum_of_permutations(M)                       # needs sage.graphs
            Traceback (most recent call last):
            ...
            ValueError: The matrix is not bistochastic
        """

        from sage.combinat.permutation import bistochastic_as_sum_of_permutations
        return bistochastic_as_sum_of_permutations(self)

    def visualize_structure(self, maxsize=512):
        r"""
        Visualize the nonzero entries.

        White pixels are put at positions with zero entries. If 'maxsize'
        is given, then the maximal dimension in either x or y direction is
        set to 'maxsize' depending on which is bigger. If the image is
        scaled, the darkness of the pixel reflects how many of the
        represented entries are nonzero. So if e.g. one image pixel
        actually represents a 2x2 submatrix, the dot is darker the more of
        the four values are nonzero.

        INPUT:

        - ``maxsize`` -- integer (default: `512`); maximal dimension
          in either x or y direction of the resulting image. If
          ``None`` or a maxsize larger than
          ``max(self.nrows(),self.ncols())`` is given the image will
          have the same pixelsize as the matrix dimensions.

        OUTPUT:

        Bitmap image as an instance of
        :class:`~sage.repl.image.Image`.

        EXAMPLES::

            sage: # needs sage.rings.real_mpfr
            sage: M = random_matrix(CC, 5, 7)
            sage: for i in range(5):  M[i,i] = 0
            sage: M[4, 0] = M[0, 6] = M[4, 6] = 0
            sage: img = M.visualize_structure();  img                                   # needs pillow
            7x5px 24-bit RGB image

        You can use :meth:`~sage.repl.image.Image.save` to save the
        resulting image::

            sage: # needs pillow sage.rings.real_mpfr
            sage: filename = tmp_filename(ext='.png')
            sage: img.save(filename)
            sage: with open(filename, 'rb') as fobj:
            ....:     fobj.read().startswith(b'\x89PNG')
            True

        TESTS:

        Test :issue:`17341`::

            sage: random_matrix(GF(2), 8, 586, sparse=True).visualize_structure()       # needs pillow
            512x6px 24-bit RGB image
        """
        cdef Py_ssize_t x, y, _x, _y, v, bi, bisq
        cdef Py_ssize_t ir, ic
        cdef double b, fct
        mr, mc = self.nrows(), self.ncols()
        if maxsize is None:
            ir = mc
            ic = mr
            b = 1.0
        elif max(mr, mc) > maxsize:
            maxsize = float(maxsize)
            ir = int(mc * maxsize/max(mr, mc))
            ic = int(mr * maxsize/max(mr, mc))
            b = max(mr, mc)/maxsize
        else:
            ir = mc
            ic = mr
            b = 1.0
        bi = int(round(b))
        bisq = bi*bi
        fct = 255.0/bisq
        from sage.repl.image import Image
        img = Image('RGB', (ir, ic))
        pixel = img.pixels()
        for x in range(ic):
            for y in range(ir):
                v = bisq
                for _x in range(bi):
                    for _y in range(bi):
                        if not self.get_unsafe(<Py_ssize_t>(x*b + _x), <Py_ssize_t>(y*b + _y)).is_zero():
                            v -= 1  # increase darkness
                v = <Py_ssize_t>(v * fct + 0.5)
                pixel[y, x] = (v, v, v)
        return img

    def density(self):
        r"""
        Return the density of the matrix.

        By density we understand the ratio of the number of nonzero
        positions and the self.nrows() \* self.ncols(), i.e. the number of
        possible nonzero positions.

        EXAMPLES:

        First, note that the density parameter does not ensure the density
        of a matrix, it is only an upper bound.

        ::

            sage: A = random_matrix(GF(127), 200, 200, density=0.3)
            sage: A.density() <= 0.3
            True

        ::

            sage: A = matrix(QQ, 3,3, [0,1,2,3,0,0,6,7,8])
            sage: A.density()
            2/3

        ::

            sage: a = matrix([[],[],[],[]])
            sage: a.density()
            0
        """
        cdef int x, y, k
        k = 0
        nr = self.nrows()
        nc = self.ncols()
        if nc == 0 or nr == 0:
            return 0
        for x from 0 <= x < nr:
            for y from 0 <= y < nc:
                if not self.get_unsafe(x, y).is_zero():
                    k+=1
        return QQ(k)/QQ(nr*nc)

    def inverse(self):
        """
        Return the inverse of self, without changing ``self``.

        Note that one can use the Python inverse operator to obtain the
        inverse as well.

        .. SEEALSO::

              :meth:`inverse_positive_definite`

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

        TESTS::

            sage: matrix().inverse()
            []

        Test :issue:`27473`::

            sage: F.<t> = LaurentSeriesRing(GF(2))
            sage: M = Matrix([[t,1], [0,t]])
            sage: ~M
            [t^-1 t^-2]
            [   0 t^-1]
        """
        return ~self

    def adjugate(self):
        r"""
        Return the adjugate matrix of ``self`` (that is, the transpose of the
        matrix of cofactors).

        Let `M` be an `n \times n`-matrix. The adjugate matrix of `M` is the `n
        \times n`-matrix `N` whose `(i, j)`-th entry is `(-1)^{i + j}
        \det(M_{j, i})`, where `M_{j,i}` is the matrix `M` with its `j`-th row
        and `i`-th column removed. It is known to satisfy `NM = MN = \det(M)I`.

        EXAMPLES::

            sage: M = Matrix(ZZ,2,2,[5,2,3,4]); M
            [5 2]
            [3 4]
            sage: N = M.adjugate(); N
            [ 4 -2]
            [-3  5]
            sage: M * N
            [14  0]
            [ 0 14]
            sage: N * M
            [14  0]
            [ 0 14]
            sage: M = Matrix(QQ, 2, 2, [5/3,2/56, 33/13,41/10]); M
            [  5/3  1/28]
            [33/13 41/10]
            sage: N = M.adjugate(); N                                                   # needs sage.libs.pari
            [ 41/10  -1/28]
            [-33/13    5/3]
            sage: M * N                                                                 # needs sage.libs.pari
            [7363/1092         0]
            [        0 7363/1092]

        An alias is :meth:`adjoint_classical`, which replaces the deprecated
        :meth:`adjoint` method::

            sage: M.adjoint()                                                           # needs sage.libs.pari
            ...: DeprecationWarning: adjoint is deprecated. Please use adjugate instead.
            See https://github.com/sagemath/sage/issues/10501 for details.
            [ 41/10  -1/28]
            [-33/13    5/3]
            sage: M.adjoint_classical()                                                 # needs sage.libs.pari
            [ 41/10  -1/28]
            [-33/13    5/3]

        ALGORITHM:

        Use PARI whenever the method ``self._adjugate`` is included to do so in
        an inheriting class. Otherwise, use a generic division-free algorithm
        that computes the adjugate matrix from the characteristic polynomial.

        The result is cached.
        """
        if self._nrows != self._ncols:
            raise ValueError("must be a square matrix")

        X = self.fetch('adjugate')
        if X is not None:
            return X

        X = self._adjugate()
        self.cache('adjugate', X)
        return X

    adjoint = deprecated_function_alias(10501, adjugate)

    adjoint_classical = adjugate

    def _adjugate(self):
        r"""
        Return the adjugate of this matrix.

        OUTPUT: matrix; the adjugate of the matrix

        EXAMPLES:

        Here is one example to illustrate this::

            sage: A = matrix(ZZ, [[1,24],[3,5]])
            sage: A
            [ 1 24]
            [ 3  5]
            sage: A._adjugate()
            [  5 -24]
            [ -3   1]

        Secondly, here is an example over a polynomial ring::

            sage: R.<t> = QQ[]
            sage: A = matrix(R, [[-2*t^2 + t + 3/2, 7*t^2 + 1/2*t - 1,
            ....:                 -6*t^2 + t - 2/11],
            ....:                [-7/3*t^2 - 1/2*t - 1/15, -2*t^2 + 19/8*t,
            ....:                 -10*t^2 + 2*t + 1/2],
            ....:                [6*t^2 - 1/2, -1/7*t^2 + 9/4*t, -t^2 - 4*t
            ....:                 - 1/10]])
            sage: A
            [       -2*t^2 + t + 3/2       7*t^2 + 1/2*t - 1       -6*t^2 + t - 2/11]
            [-7/3*t^2 - 1/2*t - 1/15         -2*t^2 + 19/8*t     -10*t^2 + 2*t + 1/2]
            [            6*t^2 - 1/2        -1/7*t^2 + 9/4*t       -t^2 - 4*t - 1/10]
            sage: A._adjugate()
            [          4/7*t^4 + 1591/56*t^3 - 961/70*t^2 - 109/80*t 55/7*t^4 + 104/7*t^3 + 6123/1540*t^2 - 959/220*t - 1/10       -82*t^4 + 101/4*t^3 + 1035/88*t^2 - 29/22*t - 1/2]
            [   -187/3*t^4 + 13/6*t^3 + 57/10*t^2 - 79/60*t - 77/300            38*t^4 + t^3 - 793/110*t^2 - 28/5*t - 53/220 -6*t^4 + 44/3*t^3 + 4727/330*t^2 - 1147/330*t - 487/660]
            [          37/3*t^4 - 136/7*t^3 - 1777/840*t^2 + 83/80*t      292/7*t^4 + 107/14*t^3 - 323/28*t^2 - 29/8*t + 1/2   61/3*t^4 - 25/12*t^3 - 269/120*t^2 + 743/240*t - 1/15]

        Finally, an example over a general ring ``S`` that is
        not an integral domain::

            sage: # needs sage.libs.singular
            sage: R.<a,b> = QQ[]
            sage: S.<x,y> = R.quo((b^3))
            sage: A = matrix(S, [[x*y^2, 2*x], [2, x^10*y]])
            sage: A
            [ x*y^2    2*x]
            [     2 x^10*y]
            sage: A.det()
            -4*x
            sage: A.charpoly('T')
            T^2 + (-x^10*y - x*y^2)*T - 4*x
            sage: A.adjugate()
            [x^10*y   -2*x]
            [    -2  x*y^2]
            sage: A.adjugate() * A
            [-4*x    0]
            [   0 -4*x]

        TESTS:

        Ensure correct handling of very small matrices::

            sage: A = matrix(ZZ, 0, 0)
            sage: A
            []
            sage: A._adjugate()
            []
            sage: A = matrix(ZZ, [[2]])
            sage: A
            [2]
            sage: A._adjugate()
            [1]

        Ensure proper computation of the adjugate matrix even in the
        presence of non-integral powers of the variable `x`
        (:issue:`14403`)::

            sage: x = var('x')                                                          # needs sage.symbolic
            sage: Matrix([[sqrt(x),x], [1,0]]).adjugate()                               # needs sage.symbolic
            [      0      -x]
            [     -1 sqrt(x)]

        .. NOTE::

            The key feature of this implementation is that it is division-free.
            This means that it can be used as a generic implementation for any
            ring (commutative and with multiplicative identity).  The algorithm
            is described in full detail as Algorithm 3.1 in [Sei2002]_.

            Note that this method does not utilise a lookup if the adjugate has
            already been computed previously, and it does not cache the result.
            This is all left to the method `adjugate`.
        """
        n = self._ncols

        if self._nrows != n:
            raise ValueError("self must be a square matrix")

        A = self.charpoly().shift(-1)(self)
        return A if n%2 else -A

    def QR(self, full=True):
        r"""
        Return a factorization of ``self`` as a unitary matrix
        and an upper-triangular matrix.

        INPUT:

        - ``full`` -- (default: ``True``) if ``True`` then the
          returned matrices have dimensions as described below.
          If ``False`` the ``R`` matrix has no zero rows and the
          columns of ``Q`` are a basis for the column space of
          ``self``.

        OUTPUT:

        If ``self`` is an `m\times n` matrix and ``full=True`` then this
        method returns a pair of matrices: `Q` is an `m\times m` unitary
        matrix (meaning its inverse is its conjugate-transpose) and `R`
        is an `m\times n` upper-triangular matrix with nonnegative entries
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

        .. NOTE::

            This is an exact computation, so limited to exact rings.
            Also the base ring needs to have a fraction field implemented
            in Sage and this field must contain square roots.  One example
            is the field of algebraic numbers, ``QQbar``, as used in the
            examples below.  If you need numerical results, convert the
            base ring to the field of complex double numbers, ``CDF``,
            which will use a faster routine that is careful about
            numerical subtleties.

        ALGORITHM:

        "Modified Gram-Schmidt," Algorithm 8.1 of [TB1997]_.

        EXAMPLES:

        For a nonsingular matrix, the QR decomposition is unique. ::

            sage: # needs sage.rings.number_field
            sage: A = matrix(QQbar, [[-2, 0, -4, -1, -1],
            ....:                    [-2, 1, -6, -3, -1],
            ....:                    [1, 1, 7, 4, 5],
            ....:                    [3, 0, 8, 3, 3],
            ....:                    [-1, 1, -6, -6, 5]])
            sage: Q, R = A.QR()
            sage: Q
            [ -0.4588314677411235?  -0.1260506983326509?   0.3812120831224489?   -0.394573711338418?     -0.6874400625964?]
            [ -0.4588314677411235?   0.4726901187474409? -0.05198346588033394?   0.7172941251646595?     -0.2209628772631?]
            [  0.2294157338705618?   0.6617661662464172?   0.6619227988762521?  -0.1808720937375480?      0.1964114464561?]
            [  0.6882472016116853?   0.1890760474989764?  -0.2044682991293135?   0.0966302966543065?     -0.6628886317894?]
            [ -0.2294157338705618?   0.5357154679137663?   -0.609939332995919?   -0.536422031427112?      0.0245514308070?]
            sage: R
            [  4.358898943540674? -0.4588314677411235?   13.07669683062202?   6.194224814505168?   2.982404540317303?]
            [                   0   1.670171752907625?  0.5987408170800917?  -1.292019657909672?   6.207996892883057?]
            [                   0                    0   5.444401659866974?   5.468660610611130? -0.6827161852283857?]
            [                   0                    0                    0   1.027626039419836?  -3.619300149686620?]
            [                   0                    0                    0                    0   0.024551430807012?]
            sage: Q.conjugate_transpose()*Q
            [1.000000000000000?            0.?e-18            0.?e-17            0.?e-16            0.?e-13]
            [           0.?e-18 1.000000000000000?            0.?e-17            0.?e-16            0.?e-13]
            [           0.?e-17            0.?e-17 1.000000000000000?            0.?e-16            0.?e-13]
            [           0.?e-16            0.?e-16            0.?e-16 1.000000000000000?            0.?e-13]
            [           0.?e-13            0.?e-13            0.?e-13            0.?e-13   1.0000000000000?]
            sage: Q * R == A
            True


        An example with complex numbers in ``QQbar``, the field of algebraic
        numbers. ::

            sage: # needs sage.rings.number_field
            sage: A = matrix(QQbar, [[-8, 4*I + 1, -I + 2, 2*I + 1],
            ....:                    [1, -2*I - 1, -I + 3, -I + 1],
            ....:                    [I + 7, 2*I + 1, -2*I + 7, -I + 1],
            ....:                    [I + 2, 0, I + 12, -1]])
            sage: Q, R = A.QR()
            sage: Q
            [                          -0.7302967433402215?    0.2070566455055649? + 0.5383472783144687?*I    0.2463049809998642? - 0.0764456358723292?*I    0.2381617683194332? - 0.1036596032779695?*I]
            [                           0.0912870929175277?   -0.2070566455055649? - 0.3778783780476559?*I    0.3786559533863033? - 0.1952221495524667?*I     0.701244450214469? - 0.3643711650986595?*I]
            [   0.6390096504226938? + 0.0912870929175277?*I    0.1708217325420910? + 0.6677576817554466?*I -0.03411475806452072? + 0.04090198741767143?*I    0.3140171085506764? - 0.0825191718705412?*I]
            [   0.1825741858350554? + 0.0912870929175277?*I  -0.03623491296347385? + 0.0724698259269477?*I   0.8632284069415110? + 0.06322839976356195?*I   -0.4499694867611521? - 0.0116119181208918?*I]
            sage: R
            [                          10.95445115010333?               0.?e-18 - 1.917028951268082?*I    5.385938482134133? - 2.190890230020665?*I  -0.2738612787525831? - 2.190890230020665?*I]
            [                                           0               4.829596256417300? + 0.?e-18*I   -0.869637911123373? - 5.864879483945125?*I   0.993871898426712? - 0.3054085521207082?*I]
            [                                           0                                            0               12.00160760935814? + 0.?e-16*I -0.2709533402297273? + 0.4420629644486323?*I]
            [                                           0                                            0                                            0               1.942963944258992? + 0.?e-16*I]
            sage: Q.conjugate_transpose()*Q
            [1.000000000000000? + 0.?e-19*I            0.?e-18 + 0.?e-17*I            0.?e-17 + 0.?e-17*I            0.?e-16 + 0.?e-16*I]
            [           0.?e-18 + 0.?e-17*I 1.000000000000000? + 0.?e-17*I            0.?e-17 + 0.?e-17*I            0.?e-16 + 0.?e-16*I]
            [           0.?e-17 + 0.?e-17*I            0.?e-17 + 0.?e-17*I 1.000000000000000? + 0.?e-17*I            0.?e-16 + 0.?e-16*I]
            [           0.?e-16 + 0.?e-16*I            0.?e-16 + 0.?e-16*I            0.?e-16 + 0.?e-16*I 1.000000000000000? + 0.?e-16*I]
            sage: Q*R - A
            [            0.?e-17 0.?e-17 + 0.?e-17*I 0.?e-16 + 0.?e-16*I 0.?e-16 + 0.?e-16*I]
            [            0.?e-18 0.?e-17 + 0.?e-17*I 0.?e-16 + 0.?e-16*I 0.?e-16 + 0.?e-16*I]
            [0.?e-17 + 0.?e-18*I 0.?e-17 + 0.?e-17*I 0.?e-16 + 0.?e-16*I 0.?e-16 + 0.?e-16*I]
            [0.?e-18 + 0.?e-18*I 0.?e-18 + 0.?e-18*I 0.?e-16 + 0.?e-16*I 0.?e-16 + 0.?e-16*I]

        A rank-deficient rectangular matrix, with both values of the ``full`` keyword.  ::

            sage: # needs sage.rings.number_field
            sage: A = matrix(QQbar, [[2, -3, 3],
            ....:                    [-1, 1, -1],
            ....:                    [-1, 3, -3],
            ....:                    [-5, 1, -1]])
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
            sage: Q.conjugate_transpose() * Q
            [                 1            0.?e-18            0.?e-18            0.?e-18]
            [           0.?e-18                  1            0.?e-18            0.?e-18]
            [           0.?e-18            0.?e-18 1.000000000000000?            0.?e-18]
            [           0.?e-18            0.?e-18            0.?e-18 1.000000000000000?]

            sage: # needs sage.rings.number_field
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

            sage: # needs sage.rings.number_field
            sage: A = matrix(QQbar, [[-3*I - 3, I - 3, -12*I + 1, -2],
            ....:                    [-I - 1, -2, 5*I - 1, -I - 2],
            ....:                    [-4*I - 4, I - 5, -7*I, -I - 4]])
            sage: Q, R = A.QR(full=False)
            sage: Q
            [ -0.4160251471689219? - 0.4160251471689219?*I   0.5370861555295747? + 0.1790287185098583?*I]
            [ -0.1386750490563073? - 0.1386750490563073?*I  -0.7519206177414046? - 0.2506402059138015?*I]
            [ -0.5547001962252291? - 0.5547001962252291?*I -0.2148344622118299? - 0.07161148740394329?*I]
            sage: R
            [                        7.211102550927979?  3.328201177351375? - 5.269651864139676?*I   7.904477796209515? + 8.45917799243475?*I  4.021576422632911? - 2.634825932069838?*I]
            [                                         0                         1.074172311059150?  -1.611258466588724? - 9.13046464400277?*I 1.611258466588724? + 0.5370861555295747?*I]
            sage: Q.conjugate_transpose()*Q
            [1 0]
            [0 1]
            sage: Q*R - A
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]

        Results of full decompositions are cached and thus returned
        immutable.  ::

            sage: # needs sage.rings.number_field
            sage: A = random_matrix(QQbar, 2, 2)
            sage: Q, R = A.QR()
            sage: Q.is_mutable()
            False
            sage: R.is_mutable()
            False

        Trivial cases return trivial results of the correct size,
        and we check `Q` itself in one case.  ::

            sage: # needs sage.rings.number_field
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

            sage: A = matrix(RealField(100), 2, 2, range(4))
            sage: A.QR()
            Traceback (most recent call last):
            ...
            NotImplementedError: QR decomposition is implemented over exact rings, try CDF for numerical results, not Real Field with 100 bits of precision

        Without a fraction field, we cannot hope to run the algorithm. ::

            sage: A = matrix(Integers(6), 2, 2, range(4))
            sage: A.QR()
            Traceback (most recent call last):
            ...
            ValueError: QR decomposition needs a fraction field of Ring of integers modulo 6

        The biggest obstacle is making unit vectors, thus requiring square
        roots, though some small cases pass through.  ::

            sage: A = matrix(ZZ, 3, 3, range(9))
            sage: A.QR()                                                                # needs sage.symbolic
            Traceback (most recent call last):
            ...
            TypeError: QR decomposition unable to compute square roots in Rational Field

            sage: A = matrix(ZZ, 2, 2, range(4))
            sage: Q, R = A.QR()
            sage: Q
            [0 1]
            [1 0]
            sage: R
            [2 3]
            [0 1]
        """
        from sage.matrix.constructor import zero_matrix, matrix
        from sage.misc.functional import sqrt

        if full:
            QR = self.fetch('QR_factors')
            if QR is not None:
                return QR
        R = self.base_ring()
        if not R.is_exact():
            raise NotImplementedError('QR decomposition is implemented over exact rings, try CDF for numerical results, not %s' % R)
        try:
            F = R.fraction_field()
        except Exception:
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
                    R[row, i] = scale
                    for j in range(i+1, n):
                        R[row, j] = q.hermitian_inner_product(V[j])
                        V[j] = V[j] - R[row, j]*q
                    row = row + 1
                except TypeError:
                    raise TypeError('QR decomposition unable to compute square roots in %s' % F)
        # complete to full orthonormal basis, or reduce to truncated R
        if full:
            Qt = matrix(Q)  # as rows here
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
        Perform Gram-Schmidt orthogonalization, with no scaling to unit vectors.

        INPUT:

        - ``self`` -- is a matrix whose columns are to be orthogonalized.
          The base ring of the matrix needs to have its fraction field
          implemented.

        OUTPUT: two matrices, ``Q`` and ``R`` such that if ``A`` represents ``self``:

        - ``A = Q*R``
        - The columns of ``Q`` are an orthogonal set which spans the
          column space of ``A``
        - The conjugate-transpose of ``Q`` times ``Q`` is a diagonal matrix
        - ``R`` -- a full-rank matrix, that has all entries below the
          main diagonal equal to zero

        This is basically a "reduced" QR decomposition of ``self`` with
        the distinction that the orthogonal column vectors of ``Q`` have
        not been scaled to unit vectors, avoiding the need to take square
        roots.

        EXAMPLES:

        A rectangular matrix whose columns have full-rank.  Notice that the
        routine computes in the fraction field, requiring the column space
        check to step up to ``QQ``.  ::

            sage: A = matrix(ZZ, [[-1, -3,  0, -1],
            ....:                 [ 1,  2, -1,  2],
            ....:                 [-3, -6,  4, -7]])
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
        example verifies that the bug on :issue:`10791` is fixed.  ::

            sage: # needs sage.rings.number_field
            sage: F.<a> = QuadraticField(-5)
            sage: A = matrix(F, [[    1,   a - 3,   a - 2, a + 1],
            ....:                [    a, 2*a + 1, 3*a + 1,     1],
            ....:                [a + 1,   a - 6, 2*a - 5,     1],
            ....:                [  2*a,       a,     3*a,    -3],
            ....:                [    1,   a - 1,       a, a - 1]])
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

            sage: A = matrix(Integers(6), 2, 2, range(4))
            sage: A._gram_schmidt_noscale()
            Traceback (most recent call last):
            ...
            TypeError: Gram-Schmidt orthogonalization requires a base ring
            with a fraction field, not Ring of integers modulo 6
        """
        from sage.matrix.constructor import matrix, zero_matrix
        R = self.base_ring()
        try:
            F = R.fraction_field()
        except TypeError:
            raise TypeError("Gram-Schmidt orthogonalization requires a base ring with a fraction field, not %s" % R)
        n = self.ncols()
        B = self.columns()
        zero = F(0)
        Bstar = []
        R = zero_matrix(F, n)
        nnz = 0  # number nonzero rows in R, or number of nonzero vectors in Bstar
        for i in range(n):
            ortho = B[i]
            for j in range(nnz):
                R[j, i] = Bstar[j].hermitian_inner_product(B[i])/Bstar[j].hermitian_inner_product(Bstar[j])
                ortho = ortho - R[j, i]*Bstar[j]
            if ortho.hermitian_inner_product(ortho) != zero:
                Bstar.append(ortho)
                R[nnz, i] = 1
                nnz = nnz + 1
        R = R[0:nnz]
        if not Bstar:
            Q = matrix(F, 0, self.nrows()).transpose()
        else:
            Q = matrix(F, Bstar).transpose()
        return Q, R

    def gram_schmidt(self, orthonormal=False):
        r"""
        Perform Gram-Schmidt orthogonalization on the rows of the matrix,
        returning a new matrix and a matrix accomplishing the transformation.

        INPUT:

        - ``self`` -- a matrix whose rows are to be orthogonalized
        - ``orthonormal`` -- (default: ``False``) if ``True`` the
          returned orthogonal vectors are unit vectors.  This keyword
          is ignored if the matrix is over ``RDF`` or ``CDF`` and the
          results are always orthonormal.

        OUTPUT:

        A pair of matrices, ``G`` and ``M`` such that if ``A``
        represents ``self``, where the parenthetical properties occur
        when ``orthonormal = True``:

        - ``A = M*G``
        - The rows of ``G`` are an orthogonal (resp. orthonormal)
          set of vectors
        - ``G`` times the conjugate-transpose of ``G`` is a diagonal
          (resp. identity) matrix
        - The row space of ``G`` equals the row space of ``A``
        - ``M`` is a full-rank matrix with zeros above the diagonal

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

            sage: # needs scipy sage.rings.complex_double sage.symbolic
            sage: A = matrix(CDF, [[ 0.6454 + 0.7491*I, -0.8662 + 0.1489*I,  0.7656 - 0.00344*I],
            ....:                  [-0.2913 + 0.8057*I,  0.8321 + 0.8170*I, -0.6744 + 0.9248*I],
            ....:                  [ 0.2554 + 0.3517*I, -0.4454 - 0.1715*I,  0.8325 - 0.6282*I]])
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
            sage: (G*G.conjugate_transpose())  # random
            [0.9999999999999999                0.0                0.0]
            [               0.0 0.9999999999999997                0.0]
            [               0.0                0.0                1.0]

        A rectangular matrix.  Note that the ``orthonormal`` keyword
        is ignored in these cases.  ::

            sage: # needs scipy
            sage: A = matrix(RDF, [[-0.978325, -0.751994, 0.925305, -0.200512, 0.420458],
            ....:                  [-0.474877, -0.983403, 0.089836,  0.132218, 0.672965]])
            sage: G, M = A.gram_schmidt(orthonormal=False)
            sage: G.round(6).zero_at(10^-6)
            [-0.607223 -0.466745  0.574315 -0.124453  0.260968]
            [ 0.123203 -0.617909 -0.530578  0.289773  0.487368]
            sage: M.round(6).zero_at(10^-6)
            [1.611147      0.0]
            [0.958116 0.867778]
            sage: (A - M*G).zero_at(10^-12)
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

            sage: # needs scipy
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
            sage: (G*G.transpose())  # abs tol 1e-14
            [0.9999999999999997                0.0                0.0]
            [               0.0 0.9999999999999998                0.0]
            [               0.0                0.0                1.0]

        Exact Rings, Orthonormalization:

        To scale a vector to unit length requires taking
        a square root, which often takes us outside the base ring.
        For the integers and the rationals, the field of algebraic numbers
        (``QQbar``) is big enough to contain what we need, but the price
        is that the computations are very slow, hence mostly of value
        for small cases or instruction. Now we need to use the
        ``orthonormal`` keyword.  ::

            sage: # needs sage.rings.number_field
            sage: A = matrix(QQbar, [[6, -8,  1],
            ....:                    [4,  1,  3],
            ....:                    [6,  3,  3],
            ....:                    [7,  1, -5],
            ....:                    [7, -3,  5]])
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
            sage: M*G - A
            [0 0 0]
            [0 0 0]
            [0 0 0]
            [0 0 0]
            [0 0 0]
            sage: (G*G.transpose() - identity_matrix(3)).norm() < 10^-10
            True
            sage: G.row_space() == A.row_space()
            True

        After :issue:`14047`, the matrix can also be over the algebraic reals
        ``AA``::

            sage: # needs sage.rings.number_field
            sage: A = matrix(AA, [[6, -8,  1],
            ....:                 [4,  1,  3],
            ....:                 [6,  3,  3],
            ....:                 [7,  1, -5],
            ....:                 [7, -3,  5]])
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

            sage: # needs sage.rings.number_field
            sage: A = matrix(QQbar, [[  -2,    -I - 1, 4*I + 2,       -1],
            ....:                    [-4*I, -2*I + 17,       0,  9*I + 1],
            ....:                    [   1,  -2*I - 6, -I + 11, -5*I + 1]])
            sage: G, M = A.gram_schmidt(orthonormal=True)
            sage: (M*G - A).norm() < 10^-10
            True
            sage: id3 = G*G.conjugate().transpose()
            sage: (id3 - identity_matrix(3)).norm() < 10^-10
            True
            sage: G.row_space() == A.row_space()        # long time
            True

        A square matrix with small rank.  The zero vectors produced as a
        result of linear dependence get eliminated, so the rows of ``G``
        are a basis for the row space of ``A``.  ::

            sage: # needs sage.rings.number_field
            sage: A = matrix(QQbar, [[2, -6, 3, 8],
            ....:                    [1, -3, 2, 5],
            ....:                    [0,  0, 2, 4],
            ....:                    [2, -6, 3, 8]])
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
            sage: M*G - A
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]
            sage: (G*G.transpose() - identity_matrix(2)).norm() < 10^-10
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
            ....:                 [-1,  0, -1,  0],
            ....:                 [-1, -2, -3, -1],
            ....:                 [ 1,  1,  2,  0]])
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

            sage: # needs sage.rings.number_field
            sage: C.<z> = CyclotomicField(5)
            sage: A = matrix(C,
            ....:       [[              -z^3 - 2*z,             -z^3 - 1, 2*z^3 - 2*z^2 + 2*z,             1],
            ....:        [         z^3 - 2*z^2 + 1, -z^3 + 2*z^2 - z - 1,                  -1,       z^2 + z],
            ....:        [-1/2*z^3 - 2*z^2 + z + 1,         -z^3 + z - 2,    -2*z^3 + 1/2*z^2, 2*z^2 - z + 2]])
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
        R = self.base_ring()
        if isinstance(R, (sage.rings.abc.RealDoubleField, sage.rings.abc.ComplexDoubleField)):
            Q, R = self.transpose().QR()
            m = R.nrows()
            n = R.ncols()
            if m > n:
                Q = Q[0:m, 0:n]
                R = R[0:n, 0:n]
        elif R.is_exact():
            if self == 1:
                # this special case occurs frequently enough to deserve a shortcut
                return (self, self)

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
        Jordan canonical form is given by `P^{-1} A P`; if the matrix is
        diagonalizable, this equals to *eigendecomposition* or *spectral
        decomposition*.

        INPUT:

        - ``base_ring`` -- ring in which to compute the Jordan form

        - ``sparse`` -- (default: ``False``) if ``sparse=True``, return a sparse
          matrix

        - ``subdivide`` -- (default: ``True``) if ``subdivide=True``, the
          subdivisions for the Jordan blocks in the matrix are shown

        - ``transformation`` -- (default: ``False``) if ``transformation=True``,
          computes also the transformation matrix

        - ``eigenvalues`` -- (default: ``None``) a complete set of roots, with
          multiplicity, of the characteristic polynomial of `A`, encoded as
          a list of pairs, each having the form `(r, m)` with `r` a root and
          `m` its multiplicity. If this is ``None``, then Sage computes this
          list itself, but this is only possible over base rings in whose
          quotient fields polynomial factorization is implemented. Over all
          other rings, providing this list manually is the only way to
          compute Jordan normal forms.

        - ``check_input`` -- boolean (default: ``True``); whether
          the list ``eigenvalues`` (if provided) has to be checked for
          correctness. Set this to ``False`` for a speedup if the eigenvalues
          are known to be correct.

        .. NOTE::

            Currently, the Jordan normal form is not computed over
            inexact rings in any but the trivial cases when the matrix
            is either `0 \times 0` or `1 \times 1`.

            In the case of exact rings, this method does not compute any
            generalized form of the Jordan normal form, but is only able to
            compute the result if the characteristic polynomial of the matrix
            splits over the specific base ring.

            Note that the base ring must be a field or a ring with an
            implemented fraction field.

        EXAMPLES::

            sage: # needs sage.combinat sage.libs.pari
            sage: a = matrix(ZZ,4,[1, 0, 0, 0, 0, 1, 0, 0,
            ....:                  1, -1, 1, 0, 1, -1, 1, 2]); a
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
            sage: b = matrix(ZZ,3,3,range(9)); b
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

        Here we need to specify a field, since the eigenvalues are not defined
        in the smallest ring containing the matrix entries (:issue:`14508`)::

            sage: c = matrix([[0,1,0], [0,0,1], [1,0,0]])
            sage: c.jordan_form(CyclotomicField(3))                                     # needs sage.combinat sage.rings.number_field
            [         1|         0|         0]
            [----------+----------+----------]
            [         0|     zeta3|         0]
            [----------+----------+----------]
            [         0|         0|-zeta3 - 1]

        If you need the transformation matrix as well as the Jordan form of
        ``self``, then pass the option ``transformation=True``. For example::

            sage: # needs sage.combinat sage.libs.pari
            sage: m = matrix([[5,4,2,1], [0,1,-1,-1], [-1,-1,3,0], [1,1,-1,2]]); m
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

        Note that for matrices over inexact rings, we do not attempt to
        compute the Jordan normal form, since it is not numerically
        stable::

            sage: b = matrix(ZZ, 3, 3, range(9))
            sage: jf, p = b.jordan_form(RealField(15), transformation=True)             # needs sage.combinat
            Traceback (most recent call last):
            ...
            ValueError: Jordan normal form not implemented over inexact rings.

        TESTS::

            sage: c = matrix(ZZ, 3, [1]*9); c
            [1 1 1]
            [1 1 1]
            [1 1 1]
            sage: c.jordan_form(subdivide=False)                                        # needs sage.combinat sage.libs.pari
            [3 0 0]
            [0 0 0]
            [0 0 0]

        ::

            sage: # needs sage.combinat sage.libs.pari
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
            sage: jf, P = m.jordan_form(transformation=True)
            sage: jf == ~P*m*P
            True

        We verify that the bug from :issue:`6942` is fixed::

            sage: # needs sage.combinat sage.libs.pari
            sage: M = Matrix(GF(2), [[1,0,1,0,0,0,1], [1,0,0,1,1,1,0], [1,1,0,1,1,1,1],
            ....:                    [1,1,1,0,1,1,1], [1,1,1,0,0,1,0], [1,1,1,0,1,0,0],
            ....:                    [1,1,1,1,1,1,0]])
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

        We verify that the bug from :issue:`6932` is fixed::

            sage: M = Matrix(1, 1, [1])
            sage: M.jordan_form(transformation=True)                                    # needs sage.combinat
            ([1], [1])

        We now go through three `10 \times 10` matrices to exhibit cases where
        there are multiple blocks of the same size::

            sage: # needs sage.combinat sage.libs.pari
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

            sage: # needs sage.combinat sage.libs.pari
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

            sage: # needs sage.combinat sage.libs.pari
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

        Verify that we smoothly move to QQ from ZZ (:issue:`12693`), i.e.
        we work in the vector space over the field::

            sage: # needs sage.combinat sage.libs.pari
            sage: M = matrix(((2,2,2), (0,0,0), (-2,-2,-2)))
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
            sage: M.jordan_form(eigenvalues=[(0, 3)])                                   # needs sage.combinat
            [0 1|0]
            [0 0|0]
            [---+-]
            [0 0|0]
            sage: M.jordan_form(eigenvalues=[(0, 2)])                                   # needs sage.combinat
            Traceback (most recent call last):
            ...
            ValueError: The provided list of eigenvalues is not correct.
            sage: M.jordan_form(transformation=True, eigenvalues=[(0, 3)])              # needs sage.combinat
            (
            [0 1|0]
            [0 0|0]  [x31   0   1]
            [---+-]  [x21   0   0]
            [0 0|0], [  0   1   0]
            )

        The base ring for the matrix needs to have a fraction field
        and it needs to be implemented.  ::

            sage: A = matrix(Integers(6), 2, 2, range(4))
            sage: A.jordan_form()                                                       # needs sage.combinat
            Traceback (most recent call last):
            ...
            ValueError: Matrix entries must be from a field, not Ring of integers modulo 6

        Test for :issue:`10563`::

            sage: R = FractionField(PolynomialRing(RationalField(), 'a'))
            sage: a = R.gen()
            sage: A = matrix(R, [[1,a], [a,1]])
            sage: A.jordan_form()                                                       # needs sage.combinat
            [ a + 1|     0]
            [------+------]
            [     0|-a + 1]
        """
        from sage.matrix.constructor import block_diagonal_matrix, jordan_block, diagonal_matrix
        from sage.combinat.partition import Partition

        if self.ncols() != self.nrows():
            raise ValueError("Jordan normal form not implemented for non-square matrices.")

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

        inferred_base_ring = base_ring

        if base_ring is None:
            inferred_base_ring = self.base_ring()

        if not inferred_base_ring.is_exact():
            raise ValueError("Jordan normal form not implemented over inexact rings.")

        # Make sure we're working with a field.
        if inferred_base_ring in _Fields:
            if base_ring is not None:
                A = self.change_ring(inferred_base_ring)
            else:
                A = self
        else:
            try:
                base_field = inferred_base_ring.fraction_field()
            except (NotImplementedError, TypeError, AttributeError):
                raise ValueError("Matrix entries must be from a field, not {0}".
                                 format(inferred_base_ring))
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
        if sum(mult for (_, mult) in evals) < n:
            raise RuntimeError("Some eigenvalue does not exist in %s." % (A.base_ring()))

        # Compute the block information.  Here, ``blocks`` is a list of pairs,
        # each first entry a root and each second entry the size of a block.
        # Note that in general there is more than one block per eigenvalue!
        blocks = []
        for eval, mult in evals:
            if mult == 1:
                blocks.append((eval, 1))
            else:
                B = A - diagonal_matrix([eval]*n, sparse=sparse)
                C = B
                ranks = [n, C.rank()]
                i = 0
                while ranks[i] > ranks[i+1] and ranks[i+1] > n-mult:
                    C = B*C
                    ranks.append(C.rank())
                    i += 1
                diagram = [ranks[i]-ranks[i+1] for i in range(len(ranks)-1)]
                blocks.extend([(eval, i)
                               for i in Partition(diagram).conjugate()])

        # ``J`` is the matrix in Jordan canonical form.  Note that the blocks
        # are ordered firstly by the eigenvalues, in the same order as obeyed
        # by ``.roots()``, and secondly by size from greatest to smallest.
        J = block_diagonal_matrix([jordan_block(eval, size, sparse=sparse)
                                   for (eval, size) in blocks],
                                  subdivide=subdivide)

        if transformation:
            from itertools import groupby

            # ``jordan_chains`` is a dictionary with keys the eigenvalues.
            # For every eigenvalue, we consider all Jordan blocks and find
            # a Jordan chain for each, adding the chain (a sequence of
            # vectors) to the entry for the eigenvalue (which is a list).
            jordan_chains = {}
            for eval, _ in evals:
                jordan_chains[eval] = []

                # Let B be the matrix `A - eval Id`.
                B = A - eval

                block_sizes = [size for e, size in blocks if e == eval]
                block_size_pairs = [(val, len(list(c)))
                                    for val, c in groupby(block_sizes)]

                # Y is a list of vectors, spanning everything that we have
                # covered by the Jordan chains we developed so far.
                Y = []

                for l, count in block_size_pairs:

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
            for eval, size in blocks:
                # Find a block with the right size
                for index, chain in enumerate(jordan_chains[eval]):
                    if len(chain) == size:
                        jordan_basis += jordan_chains[eval].pop(index)
                        break

            transformation_matrix = (A.parent()(jordan_basis)).transpose()

        if transformation:
            return J, transformation_matrix
        else:
            return J

    def jordan_decomposition(self):
        r"""
        Return the Jordan decomposition of ``self``.

        The Jordan decomposition of a matrix `A` is a pair of
        matrices `(S, N)` such that

        - `A = S + N`,
        - `S` is semisimple,
        - `N` is nilpotent.

        EXAMPLES::

            sage: A = matrix(QQ, 5, 5, {(0,1): -1, (1,0): 1, (2,3): -1}); A
            [ 0 -1  0  0  0]
            [ 1  0  0  0  0]
            [ 0  0  0 -1  0]
            [ 0  0  0  0  0]
            [ 0  0  0  0  0]
            sage: S, N = A.jordan_decomposition()
            sage: S
            [ 0 -1  0  0  0]
            [ 1  0  0  0  0]
            [ 0  0  0  0  0]
            [ 0  0  0  0  0]
            [ 0  0  0  0  0]
            sage: N
            [ 0  0  0  0  0]
            [ 0  0  0  0  0]
            [ 0  0  0 -1  0]
            [ 0  0  0  0  0]
            [ 0  0  0  0  0]
            sage: A == S + N
            True
            sage: S.is_semisimple()
            True
            sage: N.is_nilpotent()
            True
            sage: A.jordan_form()
            Traceback (most recent call last):
            ...
            RuntimeError: Some eigenvalue does not exist in Rational Field.

        TESTS::

            sage: X = random_matrix(QQ, 4)
            sage: S, N = X.jordan_decomposition()
            sage: X == S + N
            True
            sage: S is X.jordan_decomposition()[0]  # result is cached
            True
            sage: N is X.jordan_decomposition()[1]  # result is cached
            True
            sage: A = matrix(ZZ, 5, 5, {(0,1): -1, (1,0): 1, (2,3): -1})
            sage: A.jordan_decomposition()
            Traceback (most recent call last):
            ...
            ValueError: unable to compute Jordan decomposition
            sage: B = A.change_ring(RR)
            sage: B.jordan_decomposition()
            Traceback (most recent call last):
            ...
            NotImplementedError: Jordan decomposition not implemented over inexact rings
        """
        if not self.base_ring().is_exact():
            raise NotImplementedError("Jordan decomposition not implemented over inexact rings")
        JD = self.fetch('jordan_decomposition')
        if JD is not None:
            return JD
        f = self.minpoly()
        h = f // f.gcd(f.diff())
        o, p, q = h.xgcd(h.diff())
        if not o.is_one():
            raise ValueError("unable to compute Jordan decomposition")
        A = self
        hq = h * q
        # very bad bound, but requires no extra computation
        # a better bound is the maximum multiplicity in the minpoly,
        #   but this requires factoring the minpoly
        for _ in range(self.nrows()):
            if not h(A):
                ret = (A, self - A)
                ret[0].set_immutable()
                ret[1].set_immutable()
                self.cache('jordan_decomposition', ret)
                return ret
            A -= hq(A)
        raise ValueError("Jordan decomposition does not exist")

    def diagonalization(self, base_field=None):
        """
        Return a diagonal matrix similar to ``self`` along with the
        transformation matrix.

        INPUT:

        - ``base_field`` -- if given, ``self`` is regarded as a matrix over it

        OUTPUT: a diagonal matrix `D` and an invertible matrix `P` such
        that `P^{-1}AP=D`, if ``self`` is a diagonalizable matrix `A`.

        EXAMPLES::

            sage: # needs sage.libs.pari
            sage: A = matrix(QQ, 4, [-4, 6, 3, 3, -3, 5, 3, 3, 3, -6, -4, -3, -3, 6, 3, 2])
            sage: A
            [-4  6  3  3]
            [-3  5  3  3]
            [ 3 -6 -4 -3]
            [-3  6  3  2]
            sage: A.is_diagonalizable()
            True
            sage: A.diagonalization()
            (
            [ 2  0  0  0]  [ 1  1  0  0]
            [ 0 -1  0  0]  [ 1  0  1  0]
            [ 0  0 -1  0]  [-1  0  0  1]
            [ 0  0  0 -1], [ 1  1 -2 -1]
            )
            sage: D, P = A.diagonalization()
            sage: P^-1*A*P == D
            True

            sage: # needs sage.libs.pari
            sage: A = matrix(QQ, 2, [0, 2, 1, 0])
            sage: A.is_diagonalizable()
            False
            sage: A.is_diagonalizable(QQbar)                                            # needs sage.rings.number_field
            True
            sage: D, P = A.diagonalization(QQbar)                                       # needs sage.rings.number_field
            sage: P^-1*A*P == D                                                         # needs sage.rings.number_field
            True

        Matrices may fail to be diagonalizable for various reasons::

            sage: A = matrix(QQ, 2, [1,2,3, 4,5,6]); A
            [1 2 3]
            [4 5 6]
            sage: A.diagonalization()
            Traceback (most recent call last):
            ...
            TypeError: not a square matrix

            sage: B = matrix(ZZ, 2, [1, 2, 3, 4]); B
            [1 2]
            [3 4]
            sage: B.diagonalization()
            Traceback (most recent call last):
            ...
            ValueError: matrix entries must be from a field

            sage: C = matrix(RR, 2, [1., 2., 3., 4.]); C
            [1.00000000000000 2.00000000000000]
            [3.00000000000000 4.00000000000000]
            sage: C.diagonalization()
            Traceback (most recent call last):
            ...
            ValueError: base field must be exact,
            but Real Field with 53 bits of precision is not

            sage: D = matrix(QQ, 2, [0, 2, 1, 0]); D
            [0 2]
            [1 0]
            sage: D.diagonalization()                                                   # needs sage.libs.pari
            Traceback (most recent call last):
            ...
            ValueError: not diagonalizable over Rational Field

            sage: E = matrix(QQ, 2, [3, 1, 0, 3]); E
            [3 1]
            [0 3]
            sage: E.diagonalization()                                                   # needs sage.libs.pari
            Traceback (most recent call last):
            ...
            ValueError: not diagonalizable
            sage: E.jordan_form()                                                       # needs sage.combinat sage.libs.pari
            [3 1]
            [0 3]
        """
        if not self.is_square():
            raise TypeError('not a square matrix')

        if base_field is not None:
            A = self.change_ring(base_field)
        else:
            A = self

        if not A.base_ring() in _Fields:
            raise ValueError('matrix entries must be from a field')
        if not A.base_ring().is_exact():
            raise ValueError('base field must be exact, but {} is not'.format(A.base_ring()))

        from sage.matrix.constructor import diagonal_matrix, matrix

        # check if the sum of algebraic multiplicities equals the number of rows
        evals = A.charpoly().roots()
        if sum(mult for (_, mult) in evals) < self._nrows:
            raise ValueError('not diagonalizable over {}'.format(A.base_ring()))

        # compute diagonalization from the eigenspaces
        dia = []
        bas = []
        for e, am in evals:
            b = (A - e).right_kernel().basis()
            if len(b) != am:
                raise ValueError('not diagonalizable')
            dia += [e]*am
            bas += b

        return diagonal_matrix(dia), matrix(bas).transpose()

    def is_diagonalizable(self, base_field=None):
        r"""
        Determine if the matrix is similar to a diagonal matrix.

        INPUT:

        - ``base_field`` -- a new field to use for entries of the matrix

        OUTPUT:

        If ``self`` is the matrix `A`, then it is diagonalizable
        if there is an invertible matrix `P` and a diagonal matrix
        `D` such that

        .. MATH::

            P^{-1}AP = D

        This routine returns ``True`` if ``self`` is diagonalizable.
        The diagonal entries of the matrix `D` are the eigenvalues
        of `A`.

        A matrix not diagonalizable over the base field may become
        diagonalizable by extending the base field to contain all of the
        eigenvalues. Over the rationals, the field of algebraic numbers,
        :mod:`sage.rings.qqbar` is a good choice.

        To obtain the matrices `D` and `P`, use the :meth:`diagonalization`
        method.

        ALGORITHM:

        For each eigenvalue, this routine checks that the algebraic
        multiplicity (number of occurrences as a root of the characteristic
        polynomial) is equal to the geometric multiplicity (dimension
        of the eigenspace), which is sufficient to ensure a basis of
        eigenvectors for the columns of `P`.

        EXAMPLES:

        A matrix that is diagonalizable over the rationals::

            sage: A = matrix(QQ, [[-7, 16, 12,  0,    6],
            ....:                 [-9, 15,  0,  12, -27],
            ....:                 [ 9, -8, 11, -12,  51],
            ....:                 [ 3, -4,  0,  -1,   9],
            ....:                 [-1,  0, -4,   4, -12]])
            sage: A.is_diagonalizable()                                                 # needs sage.libs.pari
            True
            sage: A.diagonalization()                                                   # needs sage.libs.pari
            (
            [ 2  0  0  0  0]  [    1     1     0     1     0]
            [ 0  3  0  0  0]  [  1/2     0     1     0     1]
            [ 0  0  3  0  0]  [  1/6     1  -3/2   2/3 -14/9]
            [ 0  0  0 -1  0]  [ -1/6     0  -1/4     0  -1/3]
            [ 0  0  0  0 -1], [ -1/6  -1/3   1/3  -1/3   4/9]
            )

        This is a matrix not diagonalizable over the rationals, but you
        can still get its Jordan form. ::

            sage: A = matrix(QQ, [[-3, -14, 2, -1, 15],
            ....:                 [4, 6, -2, 3, -8],
            ....:                 [-2, -14, 0, 0, 10],
            ....:                 [3, 13, -2, 0, -11],
            ....:                 [-1, 6, 1, -3, 1]])
            sage: A.is_diagonalizable()                                                 # needs sage.libs.pari
            False
            sage: A.jordan_form(subdivide=False)                                        # needs sage.combinat sage.libs.pari
            [-1  1  0  0  0]
            [ 0 -1  0  0  0]
            [ 0  0  2  1  0]
            [ 0  0  0  2  1]
            [ 0  0  0  0  2]

        If any eigenvalue of a matrix is outside the base ring, then this
        routine raises an error. However, the ring can be extended to contain
        the eigenvalues. ::

            sage: A = matrix(QQ, [[1,  0,  1,  1, -1],
            ....:                 [0,  1,  0,  4,  8],
            ....:                 [2,  1,  3,  5,  1],
            ....:                 [2, -1,  1,  0, -2],
            ....:                 [0, -1, -1, -5, -8]])

            sage: [e in QQ for e in A.eigenvalues()]                                    # needs sage.rings.number_field
            [False, False, False, False, False]
            sage: A.is_diagonalizable()                                                 # needs sage.libs.pari
            False
            sage: A.diagonalization()                                                   # needs sage.libs.pari
            Traceback (most recent call last):
            ...
            ValueError: not diagonalizable over Rational Field

            sage: [e in QQbar for e in A.eigenvalues()]                                 # needs sage.rings.number_field
            [True, True, True, True, True]
            sage: A.is_diagonalizable(base_field=QQbar)                                 # needs sage.rings.number_field
            True

        Other exact fields may be employed, though it will not always
        be possible to extend their base fields to contain all
        the eigenvalues. ::

            sage: # needs sage.rings.finite_rings
            sage: F.<b> = FiniteField(5^2)
            sage: A = matrix(F, [[      4, 3*b + 2, 3*b + 1, 3*b + 4],
            ....:                [2*b + 1,     4*b,       0,       2],
            ....:                [    4*b,   b + 2, 2*b + 3,       3],
            ....:                [    2*b,     3*b, 4*b + 4, 3*b + 3]])
            sage: A.is_diagonalizable()
            False
            sage: A.jordan_form()                                                       # needs sage.combinat
            [      4       1|      0       0]
            [      0       4|      0       0]
            [---------------+---------------]
            [      0       0|2*b + 1       1]
            [      0       0|      0 2*b + 1]

            sage: # needs sage.rings.number_field
            sage: F.<c> = QuadraticField(-7)
            sage: A = matrix(F, [[   c + 3,   2*c - 2,   -2*c + 2,     c - 1],
            ....:                [2*c + 10, 13*c + 15, -13*c - 17, 11*c + 31],
            ....:                [2*c + 10, 14*c + 10, -14*c - 12, 12*c + 30],
            ....:                [       0,   2*c - 2,   -2*c + 2,   2*c + 2]])
            sage: A.is_diagonalizable()
            True
            sage: A.diagonalization()
            (
            [    4     0     0     0]  [   1    0    1    0]
            [    0    -2     0     0]  [   4    1    0    1]
            [    0     0 c + 3     0]  [   5    1 -2/9 10/9]
            [    0     0     0 c + 3], [   1    0 -4/9  2/9]
            )

        A trivial matrix is diagonalizable, trivially. ::

            sage: A = matrix(QQ, 0, 0)
            sage: A.is_diagonalizable()                                                 # needs sage.libs.pari
            True

        A matrix must be square to be diagonalizable. ::

            sage: A = matrix(QQ, 3, 4)
            sage: A.is_diagonalizable()
            Traceback (most recent call last):
            ...
            TypeError: not a square matrix

        The matrix must have entries from a field, and it must be an exact
        field. ::

            sage: A = matrix(ZZ, 4, range(16))
            sage: A.is_diagonalizable()
            Traceback (most recent call last):
            ...
            ValueError: matrix entries must be from a field

            sage: A = matrix(RDF, 4, range(16))
            sage: A.is_diagonalizable()
            Traceback (most recent call last):
            ...
            ValueError: base field must be exact, but Real Double Field is not
        """
        if not self.is_square():
            raise TypeError('not a square matrix')

        if base_field is not None:
            A = self.change_ring(base_field)
        else:
            A = self

        if not A.base_ring() in _Fields:
            raise ValueError('matrix entries must be from a field')
        if not A.base_ring().is_exact():
            raise ValueError('base field must be exact, but {} is not'.format(A.base_ring()))

        # check if the sum of algebraic multiplicities equals to the number of rows
        evals = A.charpoly().roots()
        if sum(mult for (_, mult) in evals) < self._nrows:
            return False

        # Obtaining a generic minimal polynomial requires much more
        # computation with kernels and their dimensions than the following.
        # However, if a derived class has a fast minimal polynomial routine
        # then overriding this by checking for repeated factors might be faster.

        # check equality of algebraic multiplicity and geometric multiplicity
        for e, am in evals:
            gm = (A - e).right_kernel().dimension()
            if am != gm:
                return False

        return True

    def is_similar(self, other, transformation=False):
        r"""
        Return ``True`` if ``self`` and ``other`` are similar,
        i.e. related by a change-of-basis matrix.

        INPUT:

        - ``other`` -- a matrix, which should be square, and of the same size
          as ``self``

        - ``transformation`` -- (default: ``False``) if ``True``, the output
          may include the change-of-basis matrix (also known as the similarity
          transformation). See below for an exact description.

        OUTPUT:

        Two matrices, `A` and `B` are similar if they are square
        matrices of the same size and there is an invertible matrix
        `S` such that `A=S^{-1}BS`. `S` can be interpreted as a
        change-of-basis matrix if `A` and `B` are viewed as matrix
        representations of the same linear transformation from a vector space
        to itself.

        When ``transformation=False`` this method will return ``True`` if
        such a matrix `S` exists, otherwise it will return ``False``.  When
        ``transformation=True`` the method returns a pair. The first part
        of the pair is ``True`` or ``False`` depending on if the matrices
        are similar. The second part of the pair is the change-of-basis
        matrix when the matrices are similar and ``None`` when the matrices
        are not similar.

        When a similarity transformation matrix ``S``  is requested,
        it will satisfy ``self = S.inverse()*other*S``.

        .. RUBRIC:: rings and coefficients

        Inexact rings are not supported. Only rings having a
        fraction field can be used as coefficients.

        The base rings for the matrices are promoted to a common
        field for the similarity check using rational form over this field.

        If the fraction fields of both matrices are the same, this
        field is used. Otherwise, if the fraction fields are only
        related by a canonical coercion, the common coercion field
        is used.

        In all cases, the result is about similarity over a common field.

        .. RUBRIC:: similarity transformation

        For computation of the similarity transformation, the
        matrices are first checked to be similar over their common
        field.

        In this case, a similarity transformation is then
        searched for over the common field. If this fails, the
        matrices are promoted to the algebraic closure of their
        common field (whenever it is available) and a similarity
        transformation is looked for over the algebraic closure.

        For example, matrices over the rationals
        may be promoted to the field of algebraic numbers (``QQbar``)
        for computation of the similarity transformation.

        .. WARNING::

            When the two matrices are similar, this routine may
            fail to find the similarity transformation.  A technical
            explanation follows.

        The similarity check is accomplished with rational form, which
        will be successful for any pair of matrices over the same field.
        However, the computation of rational form does not provide a
        transformation. So we instead compute Jordan form, which does
        provide a transformation.  But Jordan form will require that
        the eigenvalues of the matrix can be represented within Sage,
        requiring the existence of the appropriate extension field.
        When this is not possible, a :exc:`RuntimeError` is raised, as
        demonstrated in an example below.

        EXAMPLES:

        The two matrices in this example were constructed to be similar.
        The computations happen in the field of algebraic numbers, but we
        are able to convert the change-of-basis matrix back to the rationals
        (which may not always be possible). ::

            sage: A = matrix(ZZ, [[-5, 2, -11],
            ....:                 [-6, 7, -42],
            ....:                 [0, 1, -6]])
            sage: B = matrix(ZZ, [[ 1, 12,  3],
            ....:                 [-1, -6, -1],
            ....:                 [ 0,  6,  1]])
            sage: A.is_similar(B)
            True

            sage: # needs sage.combinat sage.libs.pari
            sage: _, T = A.is_similar(B, transformation=True)
            sage: T
            [ 1.00000000000000? + 0.?e-14*I            0.?e-14 + 0.?e-14*I            0.?e-14 + 0.?e-14*I]
            [-0.66666666666667? + 0.?e-15*I 0.166666666666667? + 0.?e-15*I -0.83333333333334? + 0.?e-14*I]
            [ 0.66666666666667? + 0.?e-14*I            0.?e-14 + 0.?e-14*I -0.33333333333333? + 0.?e-14*I]
            sage: T.change_ring(QQ)
            [   1    0    0]
            [-2/3  1/6 -5/6]
            [ 2/3    0 -1/3]
            sage: A == T.inverse()*B*T
            True

        Other exact fields are supported.  ::

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = FiniteField(7^2)
            sage: A = matrix(F, [[2*a + 5, 6*a + 6,   a + 3],
            ....:                [  a + 3, 2*a + 2, 4*a + 2],
            ....:                [2*a + 6, 5*a + 5,     3*a]])
            sage: B = matrix(F, [[5*a + 5, 6*a + 4,   a + 1],
            ....:                [  a + 5, 4*a + 3, 3*a + 3],
            ....:                [3*a + 5,   a + 4, 5*a + 6]])
            sage: A.is_similar(B)
            True
            sage: B.is_similar(A)
            True
            sage: _, T = A.is_similar(B, transformation=True)
            sage: T
            [      1       0       0]
            [6*a + 1 4*a + 3 4*a + 2]
            [6*a + 3 3*a + 5 3*a + 6]
            sage: A == T.inverse() * B * T
            True

        Two matrices with different sets of eigenvalues, so they
        cannot possibly be similar. ::

            sage: A = matrix(QQ, [[ 2,  3, -3, -6],
            ....:                 [ 0,  1, -2, -8],
            ....:                 [-3, -3,  4,  3],
            ....:                 [-1, -2,  2,  6]])
            sage: B = matrix(QQ, [[ 1,  1,  2,  4],
            ....:                 [-1,  2, -3, -7],
            ....:                 [-2,  3, -4, -7],
            ....:                 [ 0, -1,  0,  0]])
            sage: A.eigenvalues() == B.eigenvalues()                                    # needs sage.rings.number_field
            False
            sage: A.is_similar(B, transformation=True)
            (False, None)

        Similarity is an equivalence relation, so this routine computes
        a representative of the equivalence class for each matrix, the
        rational form, as provided by :meth:`rational_form`.  The matrices
        below have identical eigenvalues (as evidenced by equal
        characteristic polynomials), but slightly different rational forms,
        and hence are not similar.  ::

            sage: A = matrix(QQ, [[ 19, -7, -29],
            ....:                 [-16, 11,  30],
            ....:                 [ 15, -7, -25]])
            sage: B = matrix(QQ, [[-38, -63,  42],
            ....:                 [ 14,  25, -14],
            ....:                 [-14, -21,  18]])
            sage: A.charpoly() == B.charpoly()                                          # needs sage.libs.pari
            True
            sage: A.rational_form()
            [  0   0 -48]
            [  1   0   8]
            [  0   1   5]
            sage: B.rational_form()
            [ 4| 0  0]
            [--+-----]
            [ 0| 0 12]
            [ 0| 1  1]
            sage: A.is_similar(B)
            False

        Obtaining the transformation between two similar matrices
        requires the Jordan form, which requires computing the
        eigenvalues of the matrix, which may not lie in the field
        used for entries of the matrix.  In this unfortunate case,
        the computation of the transformation may fail with a
        :exc:`RuntimeError`, EVEN when the matrices are similar.  This
        is not the case for matrices over the integers, rationals
        or algebraic numbers, since the computations are done in
        the algebraically closed field of algebraic numbers.
        Here is an example where the similarity is obvious by
        design, but we are not able to resurrect a similarity
        transformation.  ::

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = FiniteField(7^2)
            sage: C = matrix(F, [[  a + 2, 5*a + 4],
            ....:                [6*a + 6, 6*a + 4]])
            sage: S = matrix(F, [[0, 1],
            ....:                [1, 0]])
            sage: D = S.inverse()*C*S
            sage: C.is_similar(D)
            True
            sage: C.is_similar(D, transformation=True)                                  # needs sage.combinat
            Traceback (most recent call last):
            ...
            RuntimeError: unable to compute transformation for similar matrices
            sage: C.jordan_form()                                                       # needs sage.combinat
            Traceback (most recent call last):
            ...
            RuntimeError: Some eigenvalue does not exist in
            Finite Field in a of size 7^2.

        An example over a finite field of prime order, which uses the
        algebraic closure of this field to find the change-of-basis
        matrix::

            sage: # needs sage.combinat sage.graphs sage.rings.finite_rings
            sage: cox = posets.TamariLattice(3).coxeter_transformation()
            sage: M = cox.change_ring(GF(3))
            sage: M.is_similar(M**3, True)      # long time
            (
                  [1 0 0 0 0]
                  [0 1 1 0 2]
                  [0 0 0 0 1]
                  [1 2 0 2 1]
            True, [0 0 1 0 0]
            )

        Inexact rings and fields are not supported.  ::

            sage: A = matrix(CDF, 2, 2, range(4))
            sage: B = copy(A)
            sage: A.is_similar(B)
            Traceback (most recent call last):
            ...
            TypeError: matrix entries must come from an exact field,
            not Complex Double Field

        Base rings for the matrices need to have a fraction field. So
        in particular, the ring needs to be at least an integral domain.  ::

            sage: Z6 = Integers(6)
            sage: A = matrix(Z6, 2, 2, range(4))
            sage: A.is_similar(A)
            Traceback (most recent call last):
            ...
            ValueError: base ring of a matrix needs a fraction field,
            maybe the ring is not an integral domain

        If the fraction fields of the entries are unequal and do not
        coerce in a common field, it is an error.  ::

            sage: A = matrix(GF(3), 2, 2, range(4))
            sage: B = matrix(GF(2), 2, 2, range(4))
            sage: A.is_similar(B, transformation=True)
            Traceback (most recent call last):
            ...
            TypeError: no common canonical parent for objects with parents:
            'Full MatrixSpace of 2 by 2 dense matrices over Finite Field
            of size 3' and
            'Full MatrixSpace of 2 by 2 dense matrices over Finite Field
            of size 2'

        A matrix over the integers and a matrix over the algebraic
        numbers will be compared over the algebraic numbers (by coercion
        of ``QQ`` in ``QQbar``).  ::

            sage: A = matrix(ZZ, 2, 2, range(4))
            sage: B = matrix(QQbar, 2, 2, range(4))                                     # needs sage.rings.number_field
            sage: A.is_similar(B)                                                       # needs sage.rings.number_field
            True

        TESTS:

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
            ValueError: transformation keyword must be True or False

        Mismatched sizes raise an error::

            sage: A = matrix(2, 2, range(4))
            sage: B = matrix(3, 3, range(9))
            sage: A.is_similar(B, transformation=True)
            Traceback (most recent call last):
            ...
            ValueError: matrices do not have the same size

        Rectangular matrices and mismatched sizes raise an error::

            sage: A = matrix(3, 2, range(6))
            sage: B = copy(A)
            sage: A.is_similar(B)
            Traceback (most recent call last):
            ...
            ValueError: similarity only makes sense for square matrices
        """
        from sage.structure.element import Matrix

        if not isinstance(other, Matrix):
            raise TypeError('similarity requires a matrix as an argument, not {0}'.format(other))
        if transformation not in [True, False]:
            raise ValueError('transformation keyword must be True or False')
        if not self.is_square() or not other.is_square():
            raise ValueError('similarity only makes sense for square matrices')
        if self.nrows() != other.nrows():
            raise ValueError('matrices do not have the same size')

        # convert to fraction fields for base rings
        try:
            A = self.matrix_over_field()
            B = other.matrix_over_field()
        except TypeError:
            mesg = "base ring of a matrix needs a fraction field, "
            mesg += "maybe the ring is not an integral domain"
            raise ValueError(mesg)
        if not have_same_parent(A, B):
            A, B = coercion_model.canonical_coercion(A, B)
        # base rings are equal now, via above check

        similar = (A.rational_form() == B.rational_form())
        if not transformation:
            return similar
        elif not similar:
            return (False, None)
        else:
            # rational form routine does not provide transformation
            # so if possible, get transformations to Jordan form

            # first try to look for Jordan forms over the fraction field
            try:
                _, SA = A.jordan_form(transformation=True)
                _, SB = B.jordan_form(transformation=True)
                return (True, SB * SA.inverse())
            except (ValueError, RuntimeError):
                pass

            # now move to the algebraic closure
            # so as to contain potential complex eigenvalues
            try:
                ring = A.base_ring()
                closure = ring.algebraic_closure()
                A = A.change_ring(closure)
                B = B.change_ring(closure)
                _, SA = A.jordan_form(transformation=True)
                _, SB = B.jordan_form(transformation=True)
                return (True, SB * SA.inverse())
            except (ValueError, RuntimeError, NotImplementedError, TypeError):
                raise RuntimeError('unable to compute transformation for similar matrices')

    def symplectic_form(self):
        r"""
        Find a symplectic form for ``self`` if ``self`` is an anti-symmetric,
        alternating matrix defined over a field.

        Returns a pair (F, C) such that the rows of C form a symplectic
        basis for ``self`` and ``F = C \* self \* C.transpose()``.

        Raises a :exc:`ValueError` if not over a field, or ``self`` is not
        anti-symmetric, or ``self`` is not alternating.

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

        .. MATH::

            \{v, Av, A^2v, A^3v, \dots \}

        INPUT:

        - ``self`` -- a square matrix over a field

        - ``v`` -- a vector with a degree equal to the size of the matrix

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

        EXAMPLES::

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

        .. MATH::

            \{v, Av, A^2v, A^3v, \dots \}

        INPUT:

        - ``self`` -- a square matrix with entries from a field

        - ``v`` -- a vector with a degree equal to the size of the matrix
          and entries compatible with the entries of the matrix

        - ``var`` -- (default: ``None``) if specified as a string or
          a generator of a polynomial ring, then this will be used
          to construct a polynomial reflecting a relation of linear
          dependence on the powers `A^iv` *and* this will cause
          the polynomial to be returned along with the subspace.
          A generator must create polynomials with coefficients from
          the same field as the matrix entries.

        - ``basis`` -- (default: ``echelon``) the basis for the
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
            sage: p.degree() == E.dimension()
            True

        The polynomial has coefficients that yield a non-trivial
        relation of linear dependence on the iterates.  Or,
        equivalently, evaluating the polynomial with the matrix
        will create a matrix that annihilates the vector.  ::

            sage: A = matrix(QQ, [[15, 37/3, -16, -104/3, -29, -7/3, 35, 2/3, -29/3, -1/3],
            ....:                 [ 2, 9, -1, -6, -6, 0, 7, 0, -2, 0],
            ....:                 [24, 74/3, -29, -208/3, -58, -14/3, 70, 4/3, -58/3, -2/3],
            ....:                 [-6, -19, 3, 21, 19, 0, -21, 0, 6, 0],
            ....:                 [2, 6, -1, -6, -3, 0, 7, 0, -2, 0],
            ....:                 [-96, -296/3, 128, 832/3, 232, 65/3, -279, -16/3, 232/3, 8/3],
            ....:                 [0, 0, 0, 0, 0, 0, 3, 0, 0, 0],
            ....:                 [20, 26/3, -30, -199/3, -42, -14/3, 70, 13/3, -55/3, -2/3],
            ....:                 [18, 57, -9, -54, -57, 0, 63, 0, -15, 0],
            ....:                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 3]])
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

            sage: A.cyclic_subspace(v, var=sin(x))                                      # needs sage.symbolic
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

            sage: F.<a> = GF(2^4)                                                       # needs sage.rings.finite_rings
            sage: G = matrix(QQ, 4, range(16))
            sage: w = vector(F, 4, [1, a, a^2, a^3])                                    # needs sage.rings.finite_rings
            sage: G.cyclic_subspace(w)                                                  # needs sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            TypeError: unable to make vector entries compatible with matrix entries
        """
        import sage.rings.polynomial.polynomial_ring
        n = self.ncols()
        R = self.base_ring()
        if not isinstance(v, Vector):
            raise TypeError('first input should be a vector, not {0}'.format(v))
        if not (var is None or isinstance(var, str)):
            generator = False
            try:
                generator = var.is_gen()
            except AttributeError:
                pass
            if not generator:
                raise TypeError('polynomial variable must be a string or polynomial ring generator, not {0}'.format(var))
            elif var.base_ring() != R:
                raise TypeError('polynomial generator must be over the same ring as the matrix entries')
        if basis not in ['echelon', 'iterates']:
            raise ValueError("basis format must be 'echelon' or 'iterates', not {0}".format(basis))
        if not self.is_square():
            raise TypeError('matrix must be square, not {0} x {1}'.format(self.nrows(), self.ncols()))
        if v.degree() != n:
            raise TypeError('vector must have degree equal to the size of the matrix, not {0}'.format(v.degree()))
        if not (R in _Fields and R.is_exact()):
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
        polynomial = (var is not None)
        if polynomial:
            x = sage.rings.polynomial.polynomial_ring.polygen(R, var)
            poly = sum([poly[i] * x**i for i in range(len(poly))])
        ambient = R**n
        if basis == 'echelon':
            echelon = []
            pivot_col_row = [(v, i) for i, v in enumerate(pivots)]
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

    def cholesky(self):
        r"""
        Return the Cholesky decomposition of a Hermitian matrix.

        Applies to a positive-definite matrix. Generally, the base ring for the
        entries of the matrix needs to be a subfield of the algebraic numbers
        (``QQbar``). Examples include the rational numbers (``QQ``), some
        number fields, and real algebraic numbers and the algebraic numbers
        themselves. Symbolic matrices can also occasionally be factored.

        OUTPUT:

        For a matrix `A` the routine returns a lower triangular
        matrix `L` such that,

        .. MATH::

            A = LL^\ast

        where `L^\ast` is the conjugate-transpose. If the matrix is
        not positive-definite (for example, if it is not Hermitian)
        then a :exc:`ValueError` results.

        If possible, the output matrix will be over the fraction field
        of the base ring of the input matrix. If that fraction field
        is missing the requisite square roots but if no imaginaries
        are encountered, then the algebraic-reals will be used.
        Otherwise, the algebraic closure of the fraction field
        (typically ``QQbar``) will be used.

        ALGORITHM:

        First we ensure that the matrix `A`
        :meth:`~.Matrix.is_hermitian`. Afterwards, we attempt to
        compute a classical :meth:`block_ldlt` factorization, `A =
        LDL^{*}`, of the matrix. If that fails, then the matrix was
        not positive-definite and an error is raised. Otherwise we
        take the entrywise square-root `\sqrt{D}` of the diagonal
        matrix `D` (whose entries are the positive eigenvalues of the
        original matrix) to obtain the Cholesky factorization `A =
        \left(L\sqrt{D}\right)\left(L\sqrt{D}\right)^{*}`. If the
        necessary square roots cannot be taken in the fraction field
        of original base ring, then we move to either its algebraic
        closure or the algebraic reals, depending on whether or not
        imaginary numbers are required.

        EXAMPLES:

        This simple example has a result with entries that remain
        in the field of rational numbers::

            sage: A = matrix(QQ, [[ 4, -2,  4,  2],
            ....:                 [-2, 10, -2, -7],
            ....:                 [ 4, -2,  8,  4],
            ....:                 [ 2, -7,  4,  7]])
            sage: A.is_symmetric()
            True
            sage: L = A.cholesky(); L
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
        of algebraic numbers::

            sage: A = matrix(ZZ, [[ 78, -30, -37,  -2],
            ....:                 [-30, 102, 179, -18],
            ....:                 [-37, 179, 326, -38],
            ....:                 [ -2, -18, -38,  15]])
            sage: A.is_symmetric()
            True
            sage: L = A.cholesky(); L                                                   # needs sage.rings.number_field
            [   8.83176086632785?                    0                    0                    0]
            [ -3.396831102433787?    9.51112708681461?                    0                    0]
            [ -4.189425026335004?   17.32383862241232?   2.886751345948129?                    0]
            [-0.2264554068289192?  -1.973397116652010?  -1.649572197684645?   2.886751345948129?]
            sage: L.parent()                                                            # needs sage.rings.number_field
            Full MatrixSpace of 4 by 4 dense matrices over Algebraic Real Field
            sage: L*L.transpose() == A                                                  # needs sage.rings.number_field
            True

        Some subfields of the complex numbers, such as this number
        field of complex numbers with rational real and imaginary parts,
        allow for this computation::

            sage: # needs sage.rings.number_field
            sage: C.<I> = QuadraticField(-1)
            sage: A = matrix(C, [[        23,  17*I + 3,  24*I + 25,     21*I],
            ....:                [ -17*I + 3,        38, -69*I + 89, 7*I + 15],
            ....:                [-24*I + 25, 69*I + 89,        976, 24*I + 6],
            ....:                [     -21*I, -7*I + 15,  -24*I + 6,       28]])
            sage: A.is_hermitian()
            True
            sage: L = A.cholesky(); L
            [                4.79...?                         0                       0        0]
            [   0.62...? - 3.54...?*I                  5.00...?                       0        0]
            [   5.21...? - 5.00...?*I   13.58...? + 10.72...?*I               24.98...?        0]
            [             -4.37...?*I   -0.10...? -  0.85...?*I  -0.21...? + 0.37...?*I 2.81...?]
            sage: L.parent()
            Full MatrixSpace of 4 by 4 dense matrices over Algebraic Field
            sage: (L*L.conjugate_transpose() - A.change_ring(QQbar)).norm() < 10^-10
            True

        The field of algebraic numbers is an ideal setting for this
        computation::

            sage: # needs sage.rings.number_field
            sage: A = matrix(QQbar, [[        2,   4 + 2*I,   6 - 4*I],
            ....:                    [ -2*I + 4,        11, 10 - 12*I],
            ....:                    [  4*I + 6, 10 + 12*I,        37]])
            sage: A.is_hermitian()
            True
            sage: L = A.cholesky()
            sage: L
            [                       1.414213562373095?          0                    0]
            [2.828427124746190? - 1.414213562373095?*I          1                    0]
            [4.242640687119285? + 2.828427124746190?*I   -2*I + 2   1.732050807568878?]
            sage: L.parent()
            Full MatrixSpace of 3 by 3 dense matrices over Algebraic Field
            sage: L*L.conjugate_transpose() == A
            True

        Results are cached, hence immutable.  Use the ``copy`` function
        if you need to make a change::

            sage: A = matrix(QQ, [[ 4, -2,  4,  2],
            ....:                 [-2, 10, -2, -7],
            ....:                 [ 4, -2,  8,  4],
            ....:                 [ 2, -7,  4,  7]])
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

        The base ring need not be exact, although you should expect
        the result to be inexact (correct only in the norm) as well
        in that case::

            sage: F = RealField(100)
            sage: A = A = matrix(F, [[1.0, 2.0], [2.0, 6.0]])
            sage: L = A.cholesky(); L
            [ 1.000... 0.000...]
            [ 2.000... 1.414...]
            sage: (L*L.transpose() - A).norm() < 1e-10                                  # needs scipy
            True

        Even symbolic matrices can sometimes be factored::

            sage: A = matrix(SR, [[pi,0], [0,pi]])                                      # needs sage.symbolic
            sage: A.cholesky()                                                          # needs sage.symbolic
            [sqrt(pi)        0]
            [       0 sqrt(pi)]

        There are a variety of situations which will prevent the
        computation of a Cholesky decomposition.

        The base ring may not be able to be viewed as a subset of the
        complex numbers, implying that "Hermitian" is meaningless::

            sage: A = matrix(Integers(6), [[2, 0], [0, 4]])
            sage: A.cholesky()
            Traceback (most recent call last):
            ...
            AttributeError: 'sage.rings.finite_rings.integer_mod.IntegerMod_int'
            object has no attribute 'conjugate'

        The matrix may not be Hermitian::

            sage: F.<a> = FiniteField(5^4)                                              # needs sage.rings.finite_rings
            sage: A = matrix(F, [[2+a^3, 3], [3, 3]])                                   # needs sage.rings.finite_rings
            sage: A.cholesky()                                                          # needs sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            ValueError: matrix is not Hermitian

        The matrix may not be positive-definite::

            sage: # needs sage.rings.number_field
            sage: C.<I> = QuadraticField(-1)
            sage: B = matrix(C, [[      2, 4 - 2*I, 2 + 2*I],
            ....:                [4 + 2*I,       8,    10*I],
            ....:                [2 - 2*I,   -10*I,      -3]])
            sage: B.is_positive_definite()
            False
            sage: B.cholesky()
            Traceback (most recent call last):
            ...
            ValueError: matrix is not positive definite

        ::

            sage: A = matrix(QQ, [[21, 15, 12, -3],
            ....:                 [15, 12,  9,  12],
            ....:                 [12,  9,  7,  3],
            ....:                 [-3,  12,  3,  8]])
            sage: A.is_positive_definite()
            False
            sage: A.cholesky()
            Traceback (most recent call last):
            ...
            ValueError: matrix is not positive definite

        TESTS:

        This verifies that :issue:`11274` is resolved::

            sage: E = matrix(QQ, [[2, 1], [1, 1]])
            sage: E.is_symmetric()
            True
            sage: E.eigenvalues()                                                       # needs sage.rings.number_field
            [0.38...?, 2.61...?]
            sage: E.det()
            1
            sage: E.cholesky()                                                          # needs sage.rings.number_field
            [ 1.414213562373095?                   0]
            [0.7071067811865475? 0.7071067811865475?]

        We check that if the input is a real matrix then the output is real as
        well (:issue:`18381`)::

            sage: E = matrix(QQ, [[4, 2], [2, 10/9]])
            sage: E.cholesky().base_ring()
            Rational Field

            sage: E = matrix(QQ, [[2, 1], [1, 1]])
            sage: E.cholesky().base_ring()                                              # needs sage.rings.number_field
            Algebraic Real Field

        Check that sparse floating-point matrices can be factored
        using a toy example reported as part of :issue:`13674`::

            sage: A = matrix(RDF, [[1, 1], [1, 2]], sparse=True)
            sage: A.cholesky()
            [1.0 0.0]
            [1.0 1.0]
            sage: A = matrix(CDF, [[1, I], [-I, 2]], sparse=True)                       # needs sage.rings.number_field
            sage: A.cholesky()
            [   1.0    0.0]
            [-1.0*I    1.0]

        Try the trivial case (:issue:`33107`)::

            sage: all( matrix(R,[]).cholesky() == matrix(R,[])                          # needs sage.rings.number_field
            ....:      for R in (RR,CC,RDF,CDF,ZZ,QQ,AA,QQbar) )
            True
            sage: all( matrix(R,[]).cholesky().is_immutable()                           # needs sage.rings.number_field
            ....:      for R in (RR,CC,RDF,CDF,ZZ,QQ,AA,QQbar) )
            True
        """
        cdef Matrix C  # output matrix
        C = self.fetch('cholesky')
        if C is not None:
            return C

        cdef Py_ssize_t n = self.nrows()

        if not self.is_hermitian():
            raise ValueError("matrix is not Hermitian")

        if n == 0:
            # The trivial matrix has a trivial cholesky decomposition.
            # We special-case this after is_hermitian() to ensure that
            # the matrix is square.
            C = self.__copy__()
            C.set_immutable()
            return C

        # Use classical=True to ensure that we don't get a permuted L.
        cdef Matrix L  # block_ldlt() results
        cdef list d    # block_ldlt() results
        try:
            _, L, d = self._block_ldlt(True)
        except ValueError:
            # If the matrix was positive-definite, that would
            # have worked.
            raise ValueError("matrix is not positive definite")

        F = L.base_ring()  # field really
        zero = F.zero()

        cdef list splits = []        # square roots of diagonal entries
        cdef bint extend = False
        for X in d:
            # The X are guaranteed to be one-by-one blocks.
            x = X[0, 0]

            if x <= zero:
                raise ValueError("matrix is not positive definite")

            if not extend and x.is_square():
                sqrt = x.sqrt()
            else:
                extend = True
                if hasattr(F, 'algebraic_closure'):
                    # This can fail if, say, F is the symbolic ring which does
                    # contain the requisite roots in its own special way but
                    # does not have an algebraic_closure() method.
                    F_ac = F.algebraic_closure()
                sqrt = F_ac(x).sqrt()
            splits.append(sqrt)
        # move square root of the diagonal matrix
        # into the lower triangular matrix
        # We need a copy, to break immutability
        # and the field may have changed as well
        if extend:
            # Try to use the algebraic reals but fall back to what is
            # probably QQbar if we find an entry in "L" that won't fit
            # in AA. This was Github issue #18381.
            from sage.rings.qqbar import AA
            try:
                C = L.change_ring(AA)
            except ValueError:  # cannot coerce...
                C = L.change_ring(F_ac)
        else:
            C = L.__copy__()

        # Overwrite the (strict) upper-triangular part of "C", since a
        # priori it contains junk after _block_ldlt().
        zero = C.base_ring().zero()
        cdef Py_ssize_t i, j  # loop indices
        for i in range(n):
            C.rescale_col_c(i, splits[i], 0)
            for j in range(i+1, n):
                C.set_unsafe(i, j, zero)
        C.set_immutable()
        self.cache('cholesky', C)
        return C

    def inverse_positive_definite(self):
        r"""
        Compute the inverse of a positive-definite matrix.

        In accord with :meth:`is_positive_definite`, only Hermitian
        matrices are considered positive-definite. Positive-definite
        matrices have several factorizations (Cholesky, LDLT, et
        cetera) that allow them to be inverted in a fast,
        numerically-stable way. This method uses an appropriate
        factorization, and is akin to the ``cholinv`` and ``chol2inv``
        functions available in R, Octave, and Stata.

        You should ensure that your matrix is positive-definite before
        using this method. When in doubt, use the generic
        :meth:`inverse` method instead.

        OUTPUT:

        If the given matrix is positive-definite, the return value is
        the same as that of the :meth:`inverse` method. If the matrix
        is not positive-definite, the behavior of this function is
        undefined.

        .. SEEALSO::

              :meth:`inverse`,
              :meth:`is_positive_definite`,
              :meth:`cholesky`,
              :meth:`indefinite_factorization`

        EXAMPLES:

        A simple two-by-two matrix with rational entries::

            sage: A = matrix(QQ, [[ 2, -1],
            ....:                 [-1,  2]])
            sage: A.is_positive_definite()
            True
            sage: A.inverse_positive_definite()
            [2/3 1/3]
            [1/3 2/3]
            sage: A.inverse_positive_definite() == A.inverse()
            True

        A matrix containing real roots::

            sage: # needs sage.rings.number_field sage.symbolic
            sage: A = matrix(AA, [ [1,       0,       sqrt(2)],
            ....:                  [0,       sqrt(3), 0      ],
            ....:                  [sqrt(2), 0,       sqrt(5)] ])
            sage: A.is_positive_definite()
            True
            sage: B = matrix(AA, [ [2*sqrt(5) + 5, 0, -sqrt(8*sqrt(5) + 18)],
            ....:                  [0,             sqrt(1/3),             0],
            ....:                  [-sqrt(8*sqrt(5) + 18), 0, sqrt(5) + 2] ])
            sage: A.inverse_positive_definite() == B
            True
            sage: A*B == A.matrix_space().identity_matrix()
            True

        A Hermitian (but not symmetric) matrix with complex entries::

            sage: # needs sage.rings.number_field sage.symbolic
            sage: A = matrix(QQbar, [ [ 1,  0,        I  ],
            ....:                     [ 0,  sqrt(5),  0  ],
            ....:                     [-I,  0,        3  ] ])
            sage: A.is_positive_definite()
            True
            sage: B = matrix(QQbar, [ [ 3/2, 0,        -I/2 ],
            ....:                     [ 0,   sqrt(1/5), 0   ],
            ....:                     [ I/2, 0,         1/2 ] ])
            sage: A.inverse_positive_definite() == B
            True
            sage: A*B == A.matrix_space().identity_matrix()
            True

        TESTS:

        Check that the naive inverse agrees with the fast one for a
        somewhat-random, positive-definite matrix with integer or
        rational entries::

            sage: from sage.misc.prandom import choice
            sage: n = ZZ.random_element(5)
            sage: ring = choice([ZZ, QQ])
            sage: A = matrix.random(ring, n)
            sage: I = matrix.identity(ring, n)
            sage: A = A*A.transpose() + I
            sage: A.is_positive_definite()
            True
            sage: actual = A.inverse_positive_definite()
            sage: expected = A.inverse()
            sage: actual == expected
            True

        Check that the naive inverse agrees with the fast one for a
        somewhat-random, possibly complex, positive-definite matrix
        with algebraic entries. This test is separate from the integer
        and rational one because inverting a matrix with algebraic
        entries is harder and requires smaller test cases::

            sage: # needs sage.rings.number_field
            sage: from sage.misc.prandom import choice
            sage: n = ZZ.random_element(2)
            sage: ring = choice([AA, QQbar])
            sage: A = matrix.random(ring, n)
            sage: I = matrix.identity(ring, n)
            sage: A = A*A.conjugate_transpose() + I
            sage: A.is_positive_definite()
            True
            sage: actual = A.inverse_positive_definite()
            sage: expected = A.inverse()
            sage: actual == expected
            True
        """
        P, L, D = self.block_ldlt()

        # The default "echelonize" inverse() method works just fine for
        # triangular matrices.
        L_inv = L.inverse()

        # Take A = PLDL^{*}P^{T} and simply invert.
        return P*L_inv.conjugate_transpose()*D.inverse()*L_inv*P.transpose()

    def LU(self, pivot=None, format='plu'):
        r"""
        Finds a decomposition into a lower-triangular matrix and
        an upper-triangular matrix.

        INPUT:

        - ``pivot`` -- pivoting strategy

          - 'auto' (default) -- see if the matrix entries are
            ordered (i.e. if they have an absolute value method),
            and if so, use a the partial pivoting strategy.
            Otherwise, fall back to the nonzero strategy.  This
            is the best choice for general routines that may
            call this for matrix entries of a variety of types.

          - ``'partial'`` -- each column is examined for
            the element with the largest absolute value and the
            row containing this element is swapped into place

          - ``'nonzero'`` -- the first nonzero element in a column
            is located and the row with this element is used

        - ``format`` -- contents of output, see more discussion
          below about output

          - 'plu' (default) -- a triple; matrices P, L and U
            such that A = P*L*U

          - ``'compact'`` -- a pair; row permutation as a tuple, and the
            matrices L and U combined into one matrix

        OUTPUT:

        Suppose that `A` is an `m\times n` matrix, then an LU
        decomposition is a lower-triangular `m\times m` matrix
        `L` with every diagonal element equal to 1, and an
        upper-triangular `m\times n` matrix, `U` such that the
        product `LU`, after a permutation of the rows, is then
        equal to `A`.  For the 'plu' format the permutation is
        returned as an `m\times m` permutation matrix `P` such
        that

        .. MATH::

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

        .. NOTE::

            This is an exact computation, so limited to exact
            rings. If you need numerical results, convert the
            base ring to the field of real double numbers,
            ``RDF`` or the field of complex double numbers,
            ``CDF``, which will use a faster routine that
            is careful about numerical subtleties.

        ALGORITHM:

            "Gaussian Elimination with Partial Pivoting,"
            Algorithm 21.1 of [TB1997]_.

        EXAMPLES:

        Notice the difference in the `L` matrix as a result of different
        pivoting strategies.  With partial pivoting, every entry of `L`
        has absolute value 1 or less.  ::

            sage: A = matrix(QQ, [[1, -1,  0,  2,  4,  7, -1],
            ....:                 [2, -1,  0,  6,  4,  8, -2],
            ....:                 [2,  0,  1,  4,  2,  6,  0],
            ....:                 [1,  0, -1,  8, -1, -1, -3],
            ....:                 [1,  1,  2, -2, -1,  1,  3]])
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
            ....:                 [ 1,  4,  7,  8],
            ....:                 [-1, -4, -6, -6],
            ....:                 [ 0, -2, -5, -8],
            ....:                 [-2, -6, -6, -2]])
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
            ....:                 [ 1, -2,  1,  3],
            ....:                 [-4,  7, -3, -8],
            ....:                 [-3,  8, -1, -5]])
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
            ....:                 [ 3, -2,  3, -1,  0,  6],
            ....:                 [-4,  2, -3,  1, -1, -8],
            ....:                 [-2,  2, -3,  2,  1,  0],
            ....:                 [ 0, -1, -1,  0,  2,  5],
            ....:                 [-1,  2, -4, -1,  5, -3]])
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
            TypeError: base ring of the matrix must be exact,
            not Real Field with 100 bits of precision

            sage: A = matrix(Integers(6), 3, 2, range(6))
            sage: A.LU()
            Traceback (most recent call last):
            ...
            TypeError: base ring of the matrix needs a field of fractions,
            not Ring of integers modulo 6

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

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = FiniteField(5^2)
            sage: C = matrix(F, [[a + 3, 4*a + 4, 2, 4*a + 2],
            ....:                [3, 2*a + 4, 2*a + 4, 2*a + 1],
            ....:                [3*a + 1, a + 3, 2*a + 4, 4*a + 3],
            ....:                [a, 3, 3*a + 1, a]])
            sage: P, L, U = C.LU(pivot='nonzero')
            sage: P
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]
            sage: L                                                                     # needs sage.combinat
            [      1       0       0       0]
            [3*a + 3       1       0       0]
            [    2*a 4*a + 2       1       0]
            [2*a + 3       2 2*a + 4       1]
            sage: U                                                                     # needs sage.combinat
            [  a + 3 4*a + 4       2 4*a + 2]
            [      0   a + 1   a + 3 2*a + 4]
            [      0       0       1 4*a + 2]
            [      0       0       0       0]
            sage: L.base_ring()                                                         # needs sage.combinat
            Finite Field in a of size 5^2
            sage: C == P*L*U                                                            # needs sage.combinat
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
            ....:                 [-1,  4,  0, -4,  0, -4, 5, -7, -7],
            ....:                 [ 0,  0,  1, -4, -1, -3, 6, -5, -6],
            ....:                 [-2,  8, -1, -4,  2, -4, 1, -8, -7],
            ....:                 [ 1, -4,  2, -4, -3,  2, 5,  6,  4]])
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
            sage: all(A.is_mutable() for A in [P, L, U])
            True

        Partial pivoting is based on the absolute values of entries
        of a column. :issue:`12208` shows that the return value of the
        absolute value must be handled carefully.  This tests that
        situation in the case of cyclotomic fields.  ::

            sage: # needs sage.groups sage.rings.number_field
            sage: C = SymmetricGroup(5).character_table()
            sage: C.base_ring()
            Cyclotomic Field of order 1 and degree 1
            sage: P, L, U = C.LU(pivot='partial')
            sage: C == P*L*U
            True

        Check that :issue:`32736` is solved::

            sage: M = Matrix(FiniteField(11), [[2,3],[4,5]])
            sage: P, L, U = M.LU()
            sage: P.base_ring()
            Finite Field of size 11
        """
        if pivot not in [None, 'partial', 'nonzero']:
            msg = "pivot strategy must be None, 'partial' or 'nonzero', not {0}"
            raise ValueError(msg.format(pivot))
        if format not in ['compact', 'plu']:
            msg = "format must be 'plu' or 'compact', not {0}"
            raise ValueError(msg.format(format))

        # exact rings only, must have fraction field
        R = self.base_ring()
        if not R.is_exact():
            msg = 'base ring of the matrix must be exact, not {0}'
            raise TypeError(msg.format(R))
        if R not in _Fields:
            try:
                F = R.fraction_field()
            except Exception:
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
            except Exception:
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
            perm = list(range(m))
            zero = F(0)
            for k in range(d):
                max_location = -1
                if partial:
                    # abs() necessary to convert zero to the
                    # correct type for comparisons (Issue #12208)
                    max_entry = abs(zero)
                    for i in range(k, m):
                        entry = abs(M.get_unsafe(i, k))
                        if entry > max_entry:
                            max_location = i
                            max_entry = entry
                else:
                    for i in range(k, m):
                        if M.get_unsafe(i, k) != zero:
                            max_location = i
                            break
                if max_location != -1:
                    perm[k], perm[max_location] = perm[max_location], perm[k]
                    M.swap_rows(k, max_location)
                    inv = M.get_unsafe(k, k).inverse()
                    for j in range(k+1, m):
                        scale = -M.get_unsafe(j, k) * inv
                        M.set_unsafe(j, k, -scale)
                        for p in range(k+1, n):
                            M.set_unsafe(j, p, M.get_unsafe(j, p) + scale*M.get_unsafe(k, p))
            perm = tuple(perm)
            M.set_immutable()
            compact = (perm, M)
            self.cache(key, compact)

        if format == 'compact':
            return compact
        elif format == 'plu':
            import sage.combinat.permutation
            perm = compact[0]
            M = compact[1].__copy__()
            F = M.base_ring()
            m, n = M._nrows, M._ncols
            d = min(m, n)
            zero = F(0)
            perm = [perm[i]+1 for i in range(m)]
            P = sage.combinat.permutation.Permutation(perm).to_matrix()
            P = P.change_ring(F)
            L = M.matrix_space(m, m).identity_matrix().__copy__()
            for i in range(1, m):
                for k in range(min(i, d)):
                    L[i, k] = M[i, k]
                    M[i, k] = zero
            return P, L, M

    def _indefinite_factorization(self, algorithm, check=True):
        r"""
        Utility function to decomposes a symmetric or
        Hermitian matrix into a lower triangular matrix
        and tuple of elements for the diagonal of a diagonal matrix.

        INPUT:

        - ``self`` -- a matrix that is symmetric or Hermitian,
          over a ring that has a fraction field implemented

        - ``algorithm`` -- ``'symmetric'`` or ``'hermitian'``,
          according to the corresponding property of the matrix

        - ``check`` -- (default: ``True``) if ``True`` then
          performs the check that the matrix is consistent with the
          ``algorithm`` keyword

        OUTPUT:

        Given a square matrix ``A``, the routine returns a
        pair: a matrix ``L`` and a list ``d``.

        ``L`` is a unit lower-triangular matrix.  ``d`` is
        the entries of a diagonal matrix.  Suppose this diagonal
        matrix is ``D``.  Then, for a symmetric matrix, these items
        are related as:

        .. MATH::

            A = LDL^T

        For a Hermitian matrix, the transpose can be replaced by
        the conjugate-transpose.

        If any leading principal submatrix is singular, then the
        computation cannot be performed and a :exc:`ValueError` results.

        Results are cached, and hence are immutable.  Caching
        eliminates redundant computations across
        :meth:`indefinite_factorization_`, :meth:`is_positive_definite`
        and :meth:`cholesky_decomposition`.

        EXAMPLES:

        A simple symmetric matrix. ::

            sage: A = matrix(QQ, [[ 4, -2,  4,  2],
            ....:                 [-2, 10, -2, -7],
            ....:                 [ 4, -2,  8,  4],
            ....:                 [ 2, -7,  4,  7]])
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

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ)
            sage: C.<I> = NumberField(x^2 + 1)
            sage: A = matrix(C, [[        23,  17*I + 3,  24*I + 25,     21*I],
            ....:                [ -17*I + 3,        38, -69*I + 89, 7*I + 15],
            ....:                [-24*I + 25, 69*I + 89,        976, 24*I + 6],
            ....:                [     -21*I, -7*I + 15,  -24*I + 6,       28]])
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
            sage: A.indefinite_factorization(algorithm='symmetric')
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

            sage: A = matrix([[3, 2+3*I], [5+6*I, 12]])                                 # needs sage.rings.number_field
            sage: A._indefinite_factorization('hermitian', check=True)
            Traceback (most recent call last):
            ...
            ValueError: matrix is not hermitian

        Results are cached and hence immutable, according
        to the ``algorithm``.  ::

            sage: A = matrix(QQ, [[ 4, -2,  4,  2],
            ....:                 [-2, 10, -2, -7],
            ....:                 [ 4, -2,  8,  4],
            ....:                 [ 2, -7,  4,  7]])
            sage: Ls, ds = A._indefinite_factorization('symmetric')
            sage: Lh, dh = A._indefinite_factorization('hermitian')
            sage: Ls.is_immutable(), Lh.is_immutable()
            (True, True)
            sage: isinstance(ds, tuple), isinstance(dh, tuple)
            (True, True)

        We check that :issue:`16633` is fixed::

            sage: A = matrix(QQ, [[ 4, -2,  4,  2],
            ....:                 [-2, 10, -2, -7],
            ....:                 [ 4, -2,  8,  4],
            ....:                 [ 2, -7,  4,  7]])
            sage: A.set_immutable()
            sage: L,d = A._indefinite_factorization('symmetric')
            sage: A
            [ 4 -2  4  2]
            [-2 10 -2 -7]
            [ 4 -2  8  4]
            [ 2 -7  4  7]
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

        if algorithm not in ['symmetric', 'hermitian']:
            msg = "'algorithm' must be 'symmetric' or 'hermitian', not {0}"
            raise ValueError(msg.format(algorithm))
        cache_string = 'indefinite_factorization_' + algorithm
        factors = self.fetch(cache_string)
        if factors is None:
            R = self.base_ring()
            if not self.is_square():
                msg = "matrix must be square, not {0} x {1}"
                raise ValueError(msg.format(self.nrows(), self.ncols()))
            if algorithm not in ['symmetric', 'hermitian']:
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
            if R is F:
                L = self.__copy__()
            else:
                L = self.change_ring(F)

            m = L._nrows
            zero = F.zero()
            one = F.one()
            d = []
            d_inv = []
            for i in range(m):
                for j in range(i+1):
                    t = L.get_unsafe(i, j)
                    if conjugate:
                        for k in range(j):
                            t -= L.get_unsafe(k, i)*L.get_unsafe(j, k).conjugate()
                    else:
                        for k in range(j):
                            t -= L.get_unsafe(k, i)*L.get_unsafe(j, k)
                    if i == j:
                        if not t:
                            self.cache(cache_string, (False, i+1))
                            return (False, i+1)
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
            factors = (L, tuple(d))
            self.cache(cache_string, factors)
        return factors

    def indefinite_factorization(self, algorithm='symmetric', check=True):
        r"""
        Decomposes a symmetric or Hermitian matrix into a
        lower triangular matrix and a diagonal matrix.

        INPUT:

        - ``self`` -- a square matrix over a ring; the base ring
          must have an implemented fraction field

        - ``algorithm`` -- (default: ``'symmetric'``) either
          ``'symmetric'`` or ``'hermitian'``, according to whether
          the input matrix is symmetric or hermitian

        - ``check`` -- (default: ``True``) if ``True`` then
          performs the check that the matrix is consistent with the
          ``algorithm`` keyword

        OUTPUT:

        A lower triangular matrix `L` with each diagonal element
        equal to `1` and a vector of entries that form a
        diagonal matrix `D`.  The vector of diagonal entries
        can be easily used to form the matrix, as demonstrated
        below in the examples.

        For a symmetric matrix, `A`, these will be related by

        .. MATH::

            A = LDL^T

        If `A` is Hermitian matrix, then the transpose of `L`
        should be replaced by the conjugate-transpose of `L`.

        If any leading principal submatrix (a square submatrix
        in the upper-left corner) is singular then this method will
        fail with a :exc:`ValueError`.

        ALGORITHM:

        The algorithm employed only uses field operations,
        but the computation of each diagonal entry has the potential
        for division by zero.  The number of operations is of order
        `n^3/3`, which is half the count for an LU decomposition.
        This makes it an appropriate candidate for solving systems
        with symmetric (or Hermitian) coefficient matrices.

        .. SEEALSO::

            :meth:`block_ldlt`

        EXAMPLES:

        There is no requirement that a matrix be positive definite, as
        indicated by the negative entries in the resulting diagonal
        matrix.  The default is that the input matrix is symmetric. ::

            sage: A = matrix(QQ, [[ 3,  -6,   9,   6,  -9],
            ....:                 [-6,  11, -16, -11,  17],
            ....:                 [ 9, -16,  28,  16, -40],
            ....:                 [ 6, -11,  16,   9, -19],
            ....:                 [-9,  17, -40, -19,  68]])
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

            sage: # needs sage.rings.number_field
            sage: C.<I> = QuadraticField(-1)
            sage: B = matrix(C, [[      2, 4 - 2*I, 2 + 2*I],
            ....:                [4 + 2*I,       8,    10*I],
            ....:                [2 - 2*I,   -10*I,      -3]])
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
            ....:                 [15, 12,  9,  6],
            ....:                 [12,  9,  7,  3],
            ....:                 [-2,  6,  3,  8]])
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

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = FiniteField(5^3)
            sage: A = matrix(F,
            ....:     [[      a^2 + 2*a, 4*a^2 + 3*a + 4,       3*a^2 + a, 2*a^2 + 2*a + 1],
            ....:      [4*a^2 + 3*a + 4,       4*a^2 + 2,             3*a, 2*a^2 + 4*a + 2],
            ....:      [      3*a^2 + a,             3*a,       3*a^2 + 2, 3*a^2 + 2*a + 3],
            ....:      [2*a^2 + 2*a + 1, 2*a^2 + 4*a + 2, 3*a^2 + 2*a + 3, 3*a^2 + 2*a + 4]])
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

        This works correctly for the 0x0 matrix::

            sage: Matrix(0).indefinite_factorization()
            ([], ())
        """
        from sage.modules.free_module_element import vector
        L, d = self._indefinite_factorization(algorithm, check=check)
        if L is False:
            msg = "{0}x{0} leading principal submatrix is singular, so cannot create indefinite factorization"
            raise ValueError(msg.format(d))
        return L, vector(L.base_ring(), d)

    cdef tuple _block_ldlt(self, bint classical):
        r"""
        Perform a user-unfriendly block-`LDL^{T}` factorization of the
        Hermitian matrix `A`

        This function is used internally to compute the factorization
        for the user-friendly :meth:`block_ldlt` method. Whereas that
        function returns three nice matrices, this one returns

          * An array ``p`` of the first `n` natural numbers, permuted
            in a way that represents the `n`-by-`n` permutation matrix
            `P`,
          * A matrix whose lower-triangular portion is ``L``, but whose
            (strict) upper-triangular portion is junk,
          * A list of the block-diagonal entries of ``D``.

        This is mainly useful to avoid having to "undo" the
        construction of the matrix `D` when we don't need it. For
        example, it's much easier to compute the inertia of a matrix
        from the list of blocks than it is from the block-diagonal
        matrix itself; given a block-diagonal matrix, you would
        first have to figure out where the blocks are!

        All of the real documentation, examples, and tests for this
        method can be found in the user-facing :meth:`block_ldlt`
        method.
        """
        cdef str cache_string = "_block_ldlt"
        if classical:
            cache_string += "_classical"
        cdef tuple result = self.fetch(cache_string)
        if result is not None:
            return result

        cdef Py_ssize_t i, j, k  # loop indices
        cdef Py_ssize_t r        # another row/column index

        # We need to construct 1x1 and 2x2 matrices to stick in d.
        from sage.matrix.constructor import matrix

        # We have to make at least one copy of the input matrix so
        # that we can change the base ring to its fraction field. Both
        # "L" and the intermediate Schur complements will potentially
        # have entries in the fraction field. However, we don't need
        # to make *two* copies.  We can't store the entries of "D" and
        # "L" in the same matrix if "D" will contain any 2x2 blocks;
        # but we can still store the entries of "L" in the copy of "A"
        # that we're going to make.  Contrast this with the non-block
        # LDL^T factorization where the entries of both "L" and "D"
        # overwrite the lower-left half of "A".
        #
        # This grants us an additional speedup, since we don't have to
        # permute the rows/columns of "L" *and* "A" at each iteration.
        #
        # Beware that the diagonals of "L" are all set to ``1`` only
        # at the end of the function, not as its columns are computed.
        ring = self.base_ring().fraction_field()

        cdef Matrix A  # A copy of the input matrix
        if self.base_ring() == ring:
            A = self.__copy__()
        else:
            # Changing the ring of a large matrix can take a loooong
            # time, compared with the short (but predictable) time we
            # might waste here checking if we need to do it.
            A = self.change_ring(ring)

        zero = ring.zero()
        one = ring.one()

        # The magic constant (1 + sqrt(17))/8 used by Bunch-Kaufman.
        # This is mainly useful for numerical stability, so we use its
        # numerical approximation to speed up the comparisons we're
        # going to make with it.
        cdef double alpha = 0.6403882032022076

        # Likewise, these two values are only ever used in comparisons
        # that determine which row/column swaps we make. It's quite
        # pointless to make long, slow comparisons when we happen to be
        # working in exact arithmetic where the process is stable anyway.
        # So, we define these constants to be C doubles, forcing any
        # comparisons to be made quickly.
        cdef double omega_1, omega_r = 0

        # Keep track of the permutations and diagonal blocks in a vector
        # rather than in a matrix, for efficiency.
        cdef Py_ssize_t n = A._nrows

        # Use a low-level array of unsigned integers for the permutation.
        from array import array
        p = array('I', range(n))

        # The list of diagonal blocks.
        cdef list d = []

        # And the parent of those diagonal blocks that are 1x1...
        one_by_one_space = A.matrix_space(1, 1)

        # The case n == 0 is *almost* handled by skipping the
        # forthcoming loop entirely. However, we must stick a trivial
        # matrix in "d" to let block_diagonal_matrix() know what its
        # base ring should be.
        if n == 0:
            d.append(A)

        k = 0
        while k < n:
            # At each step, we're considering the k-by-k submatrix
            # contained in the lower-right corner of "A", because that's
            # where we're storing the next iterate. So our indices are
            # always "k" greater than those of Higham or B&K.

            A_kk = A.get_unsafe(k, k)

            if k == (n-1):
                # Handle this trivial case manually, since otherwise the
                # algorithm's references to the e.g. "subdiagonal" are
                # meaningless. The corresponding entry of "L" will be
                # fixed later (since it's an on-diagonal element, it gets
                # set to one eventually).
                d.append(one_by_one_space(A_kk))
                k += 1
                continue

            if classical:
                try:
                    # This is a "step zero" that doesn't appear in the published algorithms.
                    # It's a back door that lets us escape with only the standard non-block
                    # non-pivoting LDL^T factorization. This allows us to implement e.g.
                    # indefinite_factorization() in terms of this method.
                    d.append(one_by_one_space(A_kk))
                    _block_ldlt_pivot1x1(A, k)
                    k += 1
                    continue
                except ZeroDivisionError:
                    raise ValueError("matrix has no classical LDL^T factorization")

            # Find the largest subdiagonal entry (in magnitude) in the
            # kth column. This occurs prior to Step (1) in Higham,
            # but is part of Step (1) in Bunch and Kaufman. We adopt
            # Higham's "omega" notation instead of B&K's "lambda"
            # because "lambda" can lead to some confusion.
            #
            # Note: omega_1 is defined as a C double, but the abs()
            # below would make a complex number approximate anyway.
            omega_1 = 0
            for i in range(k+1, n):
                a_ik_abs = A.get_unsafe(i, k).abs()
                if a_ik_abs > omega_1:
                    omega_1 = a_ik_abs
                    # We record the index "r" that corresponds to
                    # omega_1 for later. This is still part of Step
                    # (1) in B&K, but occurs later in the "else"
                    # branch of Higham's Step (1), separate from
                    # his computation of omega_1.
                    r = i

            if omega_1 == 0:
                # In this case, our matrix looks like
                #
                #   [ a 0 ]
                #   [ 0 B ]
                #
                # and we can simply skip to the next step after recording
                # the 1x1 pivot "a" in the top-left position. The entry "a"
                # will be adjusted to "1" later on to ensure that "L" is
                # (block) unit-lower-triangular.
                d.append(one_by_one_space(A_kk))
                k += 1
                continue

            if A_kk.abs() > alpha*omega_1:
                # This is the first case in Higham's Step (1), and B&K's
                # Step (2). Note that we have skipped the part of B&K's
                # Step (1) where we determine "r", since "r" is not yet
                # needed and we may waste some time computing it
                # otherwise. We are performing a 1x1 pivot, but the
                # rows/columns are already where we want them, so nothing
                # needs to be permuted.
                d.append(one_by_one_space(A_kk))
                _block_ldlt_pivot1x1(A, k)
                k += 1
                continue

            # Continuing the "else" branch of Higham's Step (1), and
            # onto B&K's Step (3) where we find the largest
            # off-diagonal entry (in magnitude) in column "r". Since
            # the matrix is Hermitian, we need only look at the
            # above-diagonal entries to find the off-diagonal of
            # maximal magnitude.
            #
            # Note: omega_r is defined as a C double, but the abs()
            # below would make a complex number approximate anyway.
            omega_r = 0
            for j in range(k, r):
                a_rj_abs = A.get_unsafe(r, j).abs()
                if a_rj_abs > omega_r:
                    omega_r = a_rj_abs

            if A_kk.abs()*omega_r >= alpha*(omega_1**2):
                # Step (2) in Higham or Step (4) in B&K.
                d.append(one_by_one_space(A_kk))
                _block_ldlt_pivot1x1(A, k)
                k += 1
                continue

            A_rr = A.get_unsafe(r, r)
            if A_rr.abs() > alpha*omega_r:
                # This is Step (3) in Higham or Step (5) in B&K. Still
                # a 1x1 pivot, but this time we need to swap
                # rows/columns k and r.
                d.append(one_by_one_space(A_rr))
                A.swap_columns_c(k, r)
                A.swap_rows_c(k, r)
                p_k = p[k]
                p[k] = p[r]
                p[r] = p_k
                _block_ldlt_pivot1x1(A, k)
                k += 1
                continue

            # If we've made it this far, we're at Step (4) in Higham
            # or Step (6) in B&K, where we perform a 2x2 pivot.  See
            # pivot1x1() for an explanation of why it's OK to permute
            # the entries of "L" here as well.
            A.swap_columns_c(k+1, r)
            A.swap_rows_c(k+1, r)
            p_k = p[k+1]
            p[k+1] = p[r]
            p[r] = p_k

            # The top-left 2x2 submatrix (starting at position k,k) is
            # now our pivot.
            E = A[k:k+2, k:k+2]
            d.append(E)

            C = A[k+2:n, k:k+2]
            B = A[k+2:, k+2:]

            # We don't actually need the inverse of E, what we really need
            # is C*E.inverse(), and that can be found by setting
            #
            #   X = C*E.inverse()   <====>   XE = C.
            #
            # Then "X" can be found easily by solving a system.  Note: I
            # do not actually know that sage solves the system more
            # intelligently, but this is still The Right Thing To Do.
            CE_inverse = E.solve_left(C)

            schur_complement = B - (CE_inverse*C.conjugate_transpose())

            # Compute the Schur complement that we'll work on during
            # the following iteration, and store it back in the lower-
            # right-hand corner of "A".
            for i in range(n-k-2):
                for j in range(i+1):
                    A.set_unsafe(k+2+i, k+2+j, schur_complement[i, j])
                    A.set_unsafe(k+2+j, k+2+i, schur_complement[j, i])

            # The on- and above-diagonal entries of "L" will be fixed
            # later, so we only need to worry about the lower-left entry
            # of the 2x2 identity matrix that belongs at the top of the
            # new column of "L".
            A.set_unsafe(k+1, k, zero)
            for i in range(n-k-2):
                for j in range(2):
                    # Store the new (k and (k+1)st) columns of "L" within
                    # the lower-left-hand corner of "A".
                    A.set_unsafe(k+i+2, k+j, CE_inverse[i, j])

            k += 2

        for i in range(n):
            # We skipped this during the main loop, but it's necessary for
            # correctness.
            A.set_unsafe(i, i, one)

        result = (p, A, d)
        self.cache(cache_string, result)
        return result

    def block_ldlt(self, classical=False):
        r"""
        Compute a block-`LDL^{T}` factorization of a Hermitian
        matrix.

        The standard `LDL^{T}` factorization of a positive-definite
        matrix `A` factors it as `A = LDL^{T}` where `L` is
        unit-lower-triangular and `D` is diagonal. If one allows
        row/column swaps via a permutation matrix `P`, then this
        factorization can be extended to many positive-semidefinite
        matrices `A` via the factorization `P^{T}AP = LDL^{T}` that
        places the zeros at the bottom of `D` to avoid division by
        zero. These factorizations extend easily to complex Hermitian
        matrices when one replaces the transpose by the
        conjugate-transpose.

        However, we can go one step further. If, in addition, we allow
        `D` to potentially contain `2 \times 2` blocks on its
        diagonal, then every real or complex Hermitian matrix `A` can
        be factored as `A = PLDL^{*}P^{T}`. When the row/column swaps
        are made intelligently, this process is numerically stable
        over inexact rings like ``RDF``.  Bunch and Kaufman describe
        such a "pivot" scheme that is suitable for the solution of
        Hermitian systems, and that is how we choose our row and
        column swaps.

        INPUT:

        - ``classical`` -- boolean (default: ``False``); whether or not to
          attempt a classical non-block `LDL^{T}` factorization with no
          row/column swaps.

        .. WARNING::

            Not all matrices have a classical `LDL^{T}` factorization.
            Set ``classical=True`` at your own risk, preferably after
            verifying that your matrix is positive-definite and (over
            inexact rings) not ill-conditioned.

        OUTPUT:

        If the input matrix is not Hermitian, the output from this
        function is undefined. Otherwise, we return a triple `(P,L,D)`
        such that `A = PLDL^{*}P^{T}` and

          * `P` is a permutation matrix,
          * `L` is unit lower-triangular,
          * `D` is a block-diagonal matrix whose blocks are of size
            one or two.

        With ``classical=True``, the permutation matrix `P` is always
        an identity matrix and the diagonal blocks are always
        one-by-one. A :exc:`ValueError` is raised if the matrix has no
        classical `LDL^{T}` factorization.

        ALGORITHM:

        We essentially follow "Algorithm A" in the paper by Bunch and
        Kaufman [BK1977]_ that describes the stable pivoting strategy.
        The same scheme is described by Higham [Hig2002]_.

        .. SEEALSO::

            :meth:`indefinite_factorization`

        REFERENCES:

        - [BK1977]_
        - [Hig2002]_

        EXAMPLES:

        This three-by-three real symmetric matrix has one positive, one
        negative, and one zero eigenvalue -- so it is not any flavor of
        (semi)definite, yet we can still factor it::

            sage: A =  matrix(QQ, [[0, 1, 0],
            ....:                  [1, 1, 2],
            ....:                  [0, 2, 0]])
            sage: P,L,D = A.block_ldlt()
            sage: P
            [0 0 1]
            [1 0 0]
            [0 1 0]
            sage: L
            [  1   0   0]
            [  2   1   0]
            [  1 1/2   1]
            sage: D
            [ 1| 0| 0]
            [--+--+--]
            [ 0|-4| 0]
            [--+--+--]
            [ 0| 0| 0]
            sage: P.transpose()*A*P == L*D*L.transpose()
            True

        This two-by-two matrix has no classical factorization, but it
        constitutes its own block-factorization::

            sage: A = matrix(QQ, [ [0,1],
            ....:                  [1,0] ])
            sage: A.block_ldlt(classical=True)
            Traceback (most recent call last):
            ...
            ValueError: matrix has no classical LDL^T factorization
            sage: A.block_ldlt()
            (
            [1 0]  [1 0]  [0 1]
            [0 1], [0 1], [1 0]
            )

        The same is true of the following complex Hermitian matrix::

            sage: # needs sage.rings.number_field
            sage: A = matrix(QQbar, [ [ 0,I],
            ....:                     [-I,0] ])
            sage: A.block_ldlt(classical=True)
            Traceback (most recent call last):
            ...
            ValueError: matrix has no classical LDL^T factorization
            sage: A.block_ldlt()
            (
            [1 0]  [1 0]  [ 0  I]
            [0 1], [0 1], [-I  0]
            )

        Complete diagonal pivoting could cause problems for the
        following matrix, since the diagonal entries are small
        compared to the off-diagonals that must be zeroed; however,
        the block algorithm refuses to factor it::

            sage: A = matrix(RDF, 2, 2, [ [1e-10, 1    ],
            ....:                         [1    , 2e-10] ])
            sage: _,L,D = A.block_ldlt(classical=True)
            sage: L*D*L.T
            [1e-10   1.0]
            [  1.0   0.0]
            sage: A.block_ldlt()                                                        # needs scipy
            (
            [1.0 0.0]  [1.0 0.0]  [1e-10   1.0]
            [0.0 1.0], [0.0 1.0], [  1.0 2e-10]
            )

        The factorization over an inexact ring is necessarily inexact,
        but `P^{T}AP` will ideally be close to `LDL^{*}` in the metric
        induced by the norm::

            sage: # needs scipy sage.rings.complex_double sage.symbolic
            sage: A = matrix(CDF, 2, 2, [ [-1.1933, -0.3185 - 1.3553*I],
            ....:                         [-0.3185 + 1.3553*I, 1.5729 ] ])
            sage: P,L,D = A.block_ldlt()
            sage: P.T*A*P == L*D*L.H
            False
            sage: (P.T*A*P - L*D*L.H).norm() < 1e-10
            True

        This matrix has a singular three-by-three leading principal
        submatrix, and therefore has no classical factorization::

            sage: A = matrix(QQ, [[21, 15, 12, -2],
            ....:                 [15, 12,  9,  6],
            ....:                 [12,  9,  7,  3],
            ....:                 [-2,  6,  3,  8]])
            sage: A[0:3,0:3].det() == 0
            True
            sage: A.block_ldlt(classical=True)
            Traceback (most recent call last):
            ...
            ValueError: matrix has no classical LDL^T factorization
            sage: A.block_ldlt()
            (
            [1 0 0 0]  [     1      0      0      0]
            [0 0 1 0]  [ -2/21      1      0      0]
            [0 0 0 1]  [   5/7  39/41      1      0]
            [0 1 0 0], [   4/7 87/164  48/79      1],
            <BLANKLINE>
            [     21|      0|      0|      0]
            [-------+-------+-------+-------]
            [      0| 164/21|      0|      0]
            [-------+-------+-------+-------]
            [      0|      0|-237/41|      0]
            [-------+-------+-------+-------]
            [      0|      0|      0| 25/316]
            )

        An indefinite symmetric matrix that happens to have a
        classical factorization::

            sage: A = matrix(QQ, [[ 3,  -6,   9,   6,  -9],
            ....:                 [-6,  11, -16, -11,  17],
            ....:                 [ 9, -16,  28,  16, -40],
            ....:                 [ 6, -11,  16,   9, -19],
            ....:                 [-9,  17, -40, -19,  68]])
            sage: A.block_ldlt(classical=True)[1:]
            (
                              [ 3| 0| 0| 0| 0]
                              [--+--+--+--+--]
                              [ 0|-1| 0| 0| 0]
                              [--+--+--+--+--]
            [ 1  0  0  0  0]  [ 0| 0| 5| 0| 0]
            [-2  1  0  0  0]  [--+--+--+--+--]
            [ 3 -2  1  0  0]  [ 0| 0| 0|-2| 0]
            [ 2 -1  0  1  0]  [--+--+--+--+--]
            [-3  1 -3  1  1], [ 0| 0| 0| 0|-1]
            )

        An indefinite Hermitian matrix that happens to have a
        classical factorization::

            sage: F.<I> = QuadraticField(-1)                                            # needs sage.rings.number_field
            sage: A = matrix(F, [[      2, 4 - 2*I, 2 + 2*I],                           # needs sage.rings.number_field
            ....:                [4 + 2*I,       8,    10*I],
            ....:                [2 - 2*I,   -10*I,      -3]])
            sage: A.block_ldlt(classical=True)[1:]                                      # needs sage.rings.number_field
            (
                                       [ 2| 0| 0]
                                       [--+--+--]
            [      1       0       0]  [ 0|-2| 0]
            [  I + 2       1       0]  [--+--+--]
            [ -I + 1 2*I + 1       1], [ 0| 0| 3]
            )

        TESTS:

        All three factors should be the identity when the input matrix is::

            sage: n = ZZ.random_element(6)
            sage: I = matrix.identity(QQ,n)
            sage: P,L,D = I.block_ldlt()
            sage: P == I and L == I and D == I
            True

        Ensure that a "random" real symmetric matrix is factored correctly::

            sage: n = ZZ.random_element(6)
            sage: A = matrix.random(QQ, n)
            sage: A = A + A.transpose()
            sage: P,L,D = A.block_ldlt()
            sage: A == P*L*D*L.transpose()*P.transpose()
            True

        Ensure that a "random" complex Hermitian matrix is factored
        correctly::

            sage: # needs sage.rings.number_field
            sage: n = ZZ.random_element(6)
            sage: F = QuadraticField(-1, 'I')
            sage: A = matrix.random(F, n)
            sage: A = A + A.conjugate_transpose()
            sage: P,L,D = A.block_ldlt()
            sage: A == P*L*D*L.conjugate_transpose()*P.conjugate_transpose()
            True

        Ensure that a "random" complex positive-semidefinite matrix is
        factored correctly and that the resulting block-diagonal matrix
        is in fact diagonal::

            sage: # needs sage.rings.number_field
            sage: n = ZZ.random_element(6)
            sage: F = QuadraticField(-1, 'I')
            sage: A = matrix.random(F, n)
            sage: A = A * A.conjugate_transpose()
            sage: P,L,D = A.block_ldlt()
            sage: A == P * L * D * L.conjugate_transpose() * P.conjugate_transpose()
            True
            sage: diagonal_matrix(D.diagonal()) == D
            True

        The factorization should be a no-op on diagonal matrices::

            sage: n = ZZ.random_element(6)
            sage: A = matrix.diagonal(random_vector(QQ, n))
            sage: I = matrix.identity(QQ, n)
            sage: P,L,D = A.block_ldlt()
            sage: P == I and L == I and A == D
            True

        All three factors have the same base ring, even when they're
        trivial::

            sage: A = matrix(QQ,0,[])
            sage: P,L,D = A.block_ldlt()
            sage: P.base_ring() == L.base_ring()
            True
            sage: L.base_ring() == D.base_ring()
            True

        Ensure that a "random" real positive-definite symmetric matrix
        has a classical factorization that agrees with
        :meth:`indefinite_factorization`::

            sage: n = ZZ.random_element(6)
            sage: A = matrix.random(QQ, n)
            sage: A = A*A.transpose() + matrix.identity(QQ, n)
            sage: _,L,D = A.block_ldlt(classical=True)
            sage: l,d = A.indefinite_factorization()
            sage: L == l and D == matrix.diagonal(d)
            True
        """
        cdef Py_ssize_t n     # size of the matrices
        cdef Py_ssize_t i, j  # loop indices
        cdef Matrix P, L, D   # output matrices

        p, L, d = self._block_ldlt(classical)
        MS = L.matrix_space()
        P = MS.matrix(lambda i, j: p[j] == i)

        # Warning: when n == 0, this works, but returns a matrix
        # whose (nonexistent) entries are in ZZ rather than in
        # the base ring of P and L. Problematic? Who knows.
        from sage.matrix.constructor import block_diagonal_matrix
        D = block_diagonal_matrix(d)

        # Overwrite the (strict) upper-triangular part of "L", since a
        # priori it contains the same entries as "A" did after _block_ldlt().
        n = L._nrows
        zero = MS.base_ring().zero()
        for i in range(n):
            for j in range(i+1, n):
                L.set_unsafe(i, j, zero)

        return (P, L, D)

    cdef bint _is_positive_definite_or_semidefinite(self, bint semi) except -1:
        """
        This is an internal wrapper that allows us to implement both
        :meth:`is_positive_definite` and
        :meth:`is_positive_semidefinite` with essentially the same
        code. The boolean ``semi`` argument exists only to change
        "greater than zero" into "greater than or equal to zero."
        """
        from sage.rings.real_lazy import RLF, CLF

        R = self.base_ring()

        if not (RLF.has_coerce_map_from(R) or
                R.has_coerce_map_from(RLF) or
                CLF.has_coerce_map_from(R) or
                R.has_coerce_map_from(CLF) or
                isinstance(R, sage.rings.abc.SymbolicRing)):
            # This is necessary to avoid "going through the motions"
            # with e.g. a one-by-one identity matrix over the finite
            # field of order 5^2, which might otherwise look positive-
            # definite.
            raise ValueError("Could not see {} as a subring of the "
                             "real or complex numbers".format(R))

        if not self.is_hermitian():
            return False

        if self._nrows == 0:
            return True  # vacuously

        cdef list d
        _, _, d = self._block_ldlt(False)

        # Check each 1x1 block for a nonpositive (negative) entry. If
        # we don't find any, the matrix is positive-(semi)definite. The
        # presence of any 2x2 blocks also indicates indefiniteness.
        import operator
        op = operator.gt
        if semi:
            op = operator.ge

        return all(d_i.nrows() == 1 and op(d_i[0, 0], 0) for d_i in d)

    def is_positive_semidefinite(self):
        r"""
        Return whether or not this matrix is positive-semidefinite.

        By SageMath convention, positive (semi)definite matrices must
        be either real symmetric or complex Hermitian.

        ALGORITHM:

        Bunch and Kaufman [BK1977]_ describe a fast,
        numerically-stable scheme for computing the "inertia" of a
        matrix by way Sylvester's inertia theorem and a
        block-`LDL^{T}` factorization. We perform this factorization,
        and read off the signs of the eigenvalues from the resulting
        diagonal blocks.

        REFERENCES:

        - [BK1977]_

        .. SEEALSO::

            :meth:`block_ldlt`, :meth:`is_positive_definite`

        EXAMPLES:

        A positive-definite matrix::

            sage: A = matrix(QQ, [ [2,1],
            ....:                  [1,2] ] )
            sage: A.eigenvalues()                                                       # needs sage.rings.number_field
            [3, 1]
            sage: A.is_positive_semidefinite()
            True

        A positive-semidefinite (but not positive-definite) matrix::

            sage: A = matrix(QQ, [ [1,1],
            ....:                  [1,1] ] )
            sage: A.eigenvalues()                                                       # needs sage.rings.number_field
            [2, 0]
            sage: A.is_positive_semidefinite()
            True

        And finally, an indefinite matrix::

            sage: A = matrix(QQ, [ [0,1],
            ....:                  [1,0] ] )
            sage: A.eigenvalues()                                                       # needs sage.rings.number_field
            [1, -1]
            sage: A.is_positive_semidefinite()
            False

        A non-Hermitian matrix cannot be positive-semidefinite,
        regardless of its eigenvalues::

            sage: A = matrix(QQ, [ [2,1],
            ....:                  [0,0] ])
            sage: A.eigenvalues()                                                       # needs sage.rings.number_field
            [2, 0]
            sage: A.is_positive_semidefinite()
            False

        Any of the preceding examples are valid over inexact rings and
        with complex numbers as well::

            sage: # needs sage.rings.complex_double sage.symbolic
            sage: A = matrix(CDF, [ [ 2, I],
            ....:                   [-I, 2] ] )
            sage: A.is_positive_semidefinite()
            True
            sage: A = matrix(CDF, [ [ 1, I],
            ....:                   [-I, 1] ] )
            sage: A.is_positive_semidefinite()
            True
            sage: A = matrix(CDF, [ [0,I],
            ....:                   [I,0] ] )
            sage: A.is_positive_semidefinite()
            False
            sage: A = matrix(CDF, [ [2,I],
            ....:                   [0,0] ])
            sage: A.is_positive_semidefinite()
            False

        TESTS:

        The trivial matrix is vacuously positive-semidefinite::

            sage: matrix(QQ, 0).is_positive_semidefinite()
            True
            sage: matrix(CDF, 0).is_positive_semidefinite()
            True

        Check that the naive and fast implementations are the same for
        a Hermitian matrix (for a non-Hermitian matrix, both "obviously"
        return ``False``)::

            sage: rings = [ZZ, QQ, RDF]
            sage: rings.append(CDF)                                                     # needs sage.rings.complex_double
            sage: rings.append(QuadraticField(-1, 'I'))                                 # needs sage.rings.number_field
            sage: from sage.misc.prandom import choice
            sage: ring = choice(rings)
            sage: A = matrix.random(ring, 10); A = A + A.conjugate_transpose()
            sage: def is_positive_semidefinite_naive(A):
            ....:     if A.nrows() == 0:
            ....:         return True
            ....:     return ( A.is_hermitian() and
            ....:              all(v >= 0 for v in A.eigenvalues()) )
            sage: expected = is_positive_semidefinite_naive(A)                          # needs scipy
            sage: actual = A.is_positive_semidefinite()                                 # needs scipy
            sage: actual == expected                                                    # needs scipy
            True

        We reject matrices whose base fields cannot be coerced to
        either real numbers, complex numbers, or symbolics; otherwise
        we risk returning nonsensical results::

            sage: F = FiniteField(5^2)                                                  # needs sage.rings.finite_rings
            sage: A = matrix.identity(F, 1)                                             # needs sage.rings.finite_rings
            sage: A.is_positive_semidefinite()                                          # needs sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            ValueError: Could not see Finite Field in z2 of size 5^2
            as a subring of the real or complex numbers
        """
        return self._is_positive_definite_or_semidefinite(True)

    def is_positive_definite(self, certificate=False):
        r"""
        Determine if a matrix is positive-definite.

        A matrix `A` is positive definite if it
        :meth:`~.Matrix.is_hermitian` and if, for every nonzero
        vector `x`,

        .. MATH::

            \left\langle Ax, x \right\rangle > 0.

        ALGORITHM:

        A Hermitian matrix is positive-definite if and only if the
        diagonal blocks in its :meth:`block_ldlt` factorization are
        all 1-by-1 and have positive entries. We first check that the
        matrix :meth:`~.Matrix.is_hermitian`, and then compute this
        factorization.

        INPUT:

        - ``self`` -- a matrix
        - ``certificate`` -- boolean (default: ``False``); return the
          lower-triangular and diagonal parts of the :meth:`block_ldlt`
          factorization when the matrix is positive-definite. Deprecated.

        OUTPUT:

        This routine will return ``True`` if the matrix is Hermitian
        and meets the condition above for the quadratic form.

        The base ring for the elements of the matrix must

        1. Have a fraction field implemented; and
        2. Be a subring of the real numbers, complex numbers,
           or symbolic ring.

        If ``certificate`` is ``True``, a triplet ``(b, L, d)`` will
        be returned instead, with ``b`` containing the result (true or
        false). If the matrix is positive-definite, then ``L`` and
        ``d`` will contain the lower-triangular and diagonal parts of
        the :meth:`block_ldlt` factorization, respectively. Or if the
        matrix is not positive-definite (that is, if ``b`` is
        ``False``), then both ``L`` and ``d`` will be ``None``.

        .. SEEALSO::

            :meth:`block_ldlt`, :meth:`~.Matrix.is_hermitian`,
            :meth:`is_positive_semidefinite`

        EXAMPLES:

        A real symmetric matrix that is positive-definite, as
        evidenced by the positive determinants of its leading
        principal submatrices::

            sage: A = matrix(QQ, [[ 4, -2,  4,  2],
            ....:                 [-2, 10, -2, -7],
            ....:                 [ 4, -2,  8,  4],
            ....:                 [ 2, -7,  4,  7]])
            sage: A.is_positive_definite()
            True
            sage: [A[:i,:i].determinant() for i in range(1, A.nrows() + 1)]
            [4, 36, 144, 144]

        A real symmetric matrix that is not positive-definite and a
        vector ``u`` that makes the corresponding quadratic form
        negative::

            sage: A = matrix(QQ, [[ 3,  -6,   9,   6,  -9],
            ....:                 [-6,  11, -16, -11,  17],
            ....:                 [ 9, -16,  28,  16, -40],
            ....:                 [ 6, -11,  16,   9, -19],
            ....:                 [-9,  17, -40, -19,  68]])
            sage: A.is_positive_definite()
            False
            sage: u = vector(QQ, [2, 2, 0, 1, 0])
            sage: (A*u).inner_product(u)
            -3

        Another real symmetric matrix that is not positive-definite
        and a vector ``u`` that makes the corresponding quadratic form
        zero::

            sage: A = matrix(QQ, [[21, 15, 12, -2],
            ....:                 [15, 12,  9,  6],
            ....:                 [12,  9,  7,  3],
            ....:                 [-2,  6,  3,  8]])
            sage: A.is_positive_definite()
            False
            sage: u = vector(QQ, [1,1,-3,0])
            sage: (A*u).inner_product(u)
            0

        A complex Hermitian matrix that is positive-definite,
        confirmed by the positive determinants of its leading
        principal submatrices::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: C.<I> = NumberField(x^2 + 1, embedding=CC(0,1))
            sage: A = matrix(C, [[        23,  17*I + 3,  24*I + 25,     21*I],
            ....:                [ -17*I + 3,        38, -69*I + 89, 7*I + 15],
            ....:                [-24*I + 25, 69*I + 89,        976, 24*I + 6],
            ....:                [     -21*I, -7*I + 15,  -24*I + 6,       28]])
            sage: A.is_positive_definite()
            True
            sage: [A[:i,:i].determinant() for i in range(1,A.nrows()+1)]
            [23, 576, 359540, 2842600]

        An Hermitian matrix that is not positive-definite and a vector
        ``u`` that makes the corresponding quadratic form negative::

            sage: # needs sage.rings.number_field
            sage: C.<I> = QuadraticField(-1)
            sage: B = matrix(C, [[      2, 4 - 2*I, 2 + 2*I],
            ....:                [4 + 2*I,       8,    10*I],
            ....:                [2 - 2*I,   -10*I,      -3]])
            sage: B.is_positive_definite()
            False
            sage: u = vector(C, [-5 + 10*I, 4 - 3*I, 0])
            sage: (B*u).hermitian_inner_product(u)
            -50

        A positive-definite matrix over an algebraically-closed field,
        confirmed by the positive determinants of its leading
        principal submatrices::

            sage: # needs sage.rings.number_field
            sage: A = matrix(QQbar, [[        2,   4 + 2*I,   6 - 4*I],
            ....:                    [ -2*I + 4,        11, 10 - 12*I],
            ....:                    [  4*I + 6, 10 + 12*I,        37]])
            sage: A.is_positive_definite()
            True
            sage: [A[:i,:i].determinant() for i in range(1, A.nrows() + 1)]
            [2, 2, 6]

        TESTS:

        If the base ring does not make sense as a subfield of the real
        numbers, complex numbers, or symbolic ring, then this routine
        will fail since comparison to zero is meaningless::

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = FiniteField(5^3)
            sage: a.conjugate()
            Traceback (most recent call last):
            ...
            TypeError: cardinality of the field must be a square number
            sage: A = matrix(F,
            ....:     [[      a^2 + 2*a, 4*a^2 + 3*a + 4,       3*a^2 + a, 2*a^2 + 2*a + 1],
            ....:      [4*a^2 + 3*a + 4,       4*a^2 + 2,             3*a, 2*a^2 + 4*a + 2],
            ....:      [      3*a^2 + a,             3*a,       3*a^2 + 2, 3*a^2 + 2*a + 3],
            ....:      [2*a^2 + 2*a + 1, 2*a^2 + 4*a + 2, 3*a^2 + 2*a + 3, 3*a^2 + 2*a + 4]])
            sage: A.is_positive_definite()
            Traceback (most recent call last):
            ...
            ValueError: Could not see Finite Field in a of size 5^3 as a subring
            of the real or complex numbers

        The 0x0 matrix is trivially positive-definite::

            sage: Matrix(0).is_positive_definite()
            True

        We can check positive-definiteness of matrices over
        approximate real/complex and symbolic rings::

            sage: # needs sage.rings.real_mpfr
            sage: matrix.identity(RR,4).is_positive_definite()
            True
            sage: matrix.identity(CC,4).is_positive_definite()
            True

            sage: matrix.identity(SR,4).is_positive_definite()                          # needs sage.symbolic
            True
        """
        result = self._is_positive_definite_or_semidefinite(False)
        if certificate:
            from sage.misc.superseded import deprecation
            msg = "the 'certificate' argument is deprecated; if you "
            msg += "need the corresponding factorization, you can "
            msg += "simply compute it yourself (the results are cached)"
            deprecation(31619, msg)
            L = None
            d = None
            if result:
                from sage.modules.free_module_element import vector
                _, L, D = self.block_ldlt()
                d = vector(D.base_ring(), D.diagonal())
            return (result, L, d)
        else:
            return result

    def principal_square_root(self, check_positivity=True):
        r"""
        Return the principal square root of a positive definite matrix.

        A positive definite matrix `A` has a unique positive definite
        matrix `M` such that `M^2 = A`.

        See :wikipedia:`Square_root_of_a_matrix`.

        EXAMPLES::

            sage: A = Matrix([[1,-1/2,0], [-1/2,1,-1/2], [0,-1/2,1]])
            sage: B = A.principal_square_root()
            sage: A == B^2
            True
        """
        from sage.matrix.special import diagonal_matrix
        from sage.misc.functional import sqrt

        if check_positivity and not self.is_positive_definite():
            return False

        d, L = self.eigenmatrix_left()
        return L.inverse() * diagonal_matrix([sqrt(a) for a in d.diagonal()]) * L

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

            sage: a = matrix(ZZ, 2, [2^10000, 3^10000, 2^50, 3^19292])
            sage: a.hadamard_bound()
            12215
            sage: len(str(a.det()))
            12215
        """
        from sage.rings.real_double import RDF
        try:
            A = self.change_ring(RDF)
            m1 = A._hadamard_row_bound()
            A = A.transpose()
            m2 = A._hadamard_row_bound()
            return min(m1, m2)
        except (OverflowError, TypeError):
            from sage.rings.real_mpfr import RealField
            # Try using MPFR, which handles large numbers much better, but is slower.
            from sage.matrix.misc_mpfr import hadamard_row_bound_mpfr
            R = RealField(53, rnd='RNDU')
            A = self.change_ring(R)
            m1 = hadamard_row_bound_mpfr(A)
            A = A.transpose()
            m2 = hadamard_row_bound_mpfr(A)
            return min(m1, m2)

    def find(self, f, indices=False):
        r"""
        Find elements in this matrix satisfying the constraints in the
        function `f`. The function is evaluated on each element of
        the matrix .

        INPUT:

        - ``f`` -- a function that is evaluated on each
          element of this matrix

        - ``indices`` -- whether or not to return the indices
          and elements of this matrix that satisfy the function

        OUTPUT: if ``indices`` is not specified, return a
        matrix with 1 where `f` is satisfied and 0 where it is not.
        If ``indices`` is specified, return a dictionary
        containing the elements of this matrix satisfying `f`.

        EXAMPLES::

            sage: M = matrix(4,3,[1, -1/2, -1, 1, -1, -1/2, -1, 0, 0, 2, 0, 1])
            sage: M.find(lambda entry: entry == 0)
            [0 0 0]
            [0 0 0]
            [0 1 1]
            [0 1 0]

        ::

            sage: M.find(lambda u: u < 0)
            [0 1 1]
            [0 1 1]
            [1 0 0]
            [0 0 0]

        ::

            sage: M = matrix(4,3,[1, -1/2, -1, 1, -1, -1/2, -1, 0, 0, 2, 0, 1])
            sage: len(M.find(lambda u:u<1 and u>-1,indices=True))
            5

        ::

            sage: M.find(lambda u: u != 1/2)
            [1 1 1]
            [1 1 1]
            [1 1 1]
            [1 1 1]

        ::

            sage: M.find(lambda u: u > 1.2)
            [0 0 0]
            [0 0 0]
            [0 0 0]
            [1 0 0]

        ::

            sage: sorted(M.find(lambda u: u != 0, indices=True).keys()) == M.nonzero_positions()
            True
        """
        from sage.matrix.matrix_space import MatrixSpace
        cdef Py_ssize_t size, i, j
        cdef object M

        if not indices:
            L = self._list()
            size = PyList_GET_SIZE(L)
            M = PyList_New(0)

            for i from 0 <= i < size:
                PyList_Append(M, <object>f(<object>PyList_GET_ITEM(L, i)))

            from sage.rings.finite_rings.integer_mod_ring import IntegerModRing

            return MatrixSpace(IntegerModRing(2),
                               nrows=self._nrows, ncols=self._ncols).matrix(M)

        else:
            # return matrix along with indices in a dictionary
            d = {}
            for i from 0 <= i < self._nrows:
                for j from 0 <= j < self._ncols:
                    if f(self.get_unsafe(i, j)):
                        d[(i, j)] = self.get_unsafe(i, j)

            return d

    def conjugate(self):
        r"""
        Return the conjugate of self, i.e. the matrix whose entries are the
        conjugates of the entries of ``self``.

        EXAMPLES::

            sage: # needs sage.rings.complex_double sage.symbolic
            sage: A = matrix(CDF, [[1+I,1],[0,2*I]])
            sage: A.conjugate()
            [1.0 - 1.0*I         1.0]
            [        0.0      -2.0*I]

        A matrix over a not-totally-real number field::

            sage: x = polygen(ZZ, 'x')
            sage: K.<j> = NumberField(x^2 + 5)                                          # needs sage.rings.number_field
            sage: M = matrix(K, [[1+j,1], [0,2*j]])                                     # needs sage.rings.number_field
            sage: M.conjugate()                                                         # needs sage.rings.number_field
            [-j + 1      1]
            [     0   -2*j]

        There is a shortcut for the conjugate::

            sage: M.C                                                                   # needs sage.rings.number_field
            [-j + 1      1]
            [     0   -2*j]

        There is also a shortcut for the conjugate transpose, or "Hermitian transpose"::

            sage: M.H                                                                   # needs sage.rings.number_field
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
        Return the transpose of ``self`` after each entry has been
        converted to its complex conjugate.

        .. NOTE::

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

            sage: M = matrix(SR, 2, 2, [[2-I, 3+4*I], [9-6*I, 5*I]])                    # needs sage.symbolic
            sage: M.base_ring()                                                         # needs sage.symbolic
            Symbolic Ring
            sage: M.conjugate_transpose()                                               # needs sage.symbolic
            [   I + 2  6*I + 9]
            [-4*I + 3     -5*I]

            sage: # needs sage.rings.real_mpfr sage.symbolic
            sage: P = matrix(CC, 3, 2, [0.95-0.63*I, 0.84+0.13*I,
            ....:                       0.94+0.23*I, 0.23+0.59*I,
            ....:                       0.52-0.41*I, -0.50+0.90*I])
            sage: P.base_ring()
            Complex Field with 53 bits of precision
            sage: P.conjugate_transpose()
            [ 0.950... + 0.630...*I  0.940... - 0.230...*I  0.520... + 0.410...*I]
            [ 0.840... - 0.130...*I  0.230... - 0.590...*I -0.500... - 0.900...*I]

        There is also a shortcut for the conjugate transpose, or "Hermitian transpose"::

            sage: M.H                                                                   # needs sage.symbolic
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

            sage: # needs sage.rings.number_field
            sage: C = CyclotomicField(5)
            sage: a = C.gen(); a
            zeta5
            sage: CC(a)
            0.309016994374947 + 0.951056516295154*I
            sage: M = matrix(C, 1, 2, [a^2, a+a^3])
            sage: M.conjugate_transpose()
            [             zeta5^3]
            [-zeta5^3 - zeta5 - 1]

        Furthermore, this method can be applied to matrices over
        quadratic extensions of finite fields::

            sage: F.<a> = GF(9,'a')                                                     # needs sage.rings.finite_rings
            sage: N = matrix(F, 2, [0,a,-a,1]); N                                       # needs sage.rings.finite_rings
            [  0   a]
            [2*a   1]
            sage: N.conjugate_transpose()                                               # needs sage.rings.finite_rings
            [      0   a + 2]
            [2*a + 1       1]

        Conjugation does not make sense over rings not containing complex
        numbers or finite fields which are not a quadratic extension::

            sage: N = matrix(GF(5), 2, [0,1,2,3])
            sage: N.conjugate_transpose()
            Traceback (most recent call last):
            ...
            AttributeError: 'sage.rings.finite_rings.integer_mod.IntegerMod_int' object
            has no attribute 'conjugate'...
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

        - ``self`` -- a matrix whose entries are coercible into ``CDF``

        - ``p`` -- one of the following options:

          - ``1`` -- the largest column-sum norm

          - ``2`` -- (default) the Euclidean norm

          - ``Infinity`` -- the largest row-sum norm

          - ``'frob'`` -- the Frobenius (sum of squares) norm

        OUTPUT: RDF number

        .. SEEALSO::

            - :func:`sage.misc.functional.norm`

        EXAMPLES::

            sage: A = matrix(ZZ, [[1,2,4,3], [-1,0,3,-10]])
            sage: A.norm(1)
            13.0
            sage: A.norm(Infinity)
            14.0
            sage: B = random_matrix(QQ, 20, 21)
            sage: B.norm(Infinity) == (B.transpose()).norm(1)
            True

        ::

            sage: Id = identity_matrix(12)
            sage: Id.norm(2)                                                            # needs scipy
            1.0

            sage: # needs scipy sage.rings.real_mpfr
            sage: A = matrix(RR, 2, 2, [13,-4,-4,7])
            sage: A.norm()  # rel tol 2e-16
            14.999999999999998

        Norms of numerical matrices over high-precision reals are computed by this routine.
        Faster routines for double precision entries from `RDF` or `CDF` are provided by
        the :class:`~sage.matrix.matrix_double_dense.Matrix_double_dense` class.  ::

            sage: # needs sage.rings.real_mpfr sage.symbolic
            sage: A = matrix(CC, 2, 3, [3*I,4,1-I,1,2,0])
            sage: A.norm('frob')
            5.656854249492381
            sage: A.norm(2)
            5.470684443210...
            sage: A.norm(1)
            6.0
            sage: A.norm(Infinity)
            8.414213562373096

            sage: a = matrix([[],[],[],[]])
            sage: a.norm()
            0.0
            sage: a.norm(Infinity) == a.norm(1)
            True

        TESTS:

        Check that a sparse zero matrix is handled (:issue:`29214`)::

            sage: matrix(CDF, 2, 2, sparse=True).norm(1)
            0.0
        """
        from sage.rings.real_double import RDF

        if self._nrows == 0 or self._ncols == 0:
            return RDF(0)

        # 2-norm:
        if p == 2:
            from sage.rings.complex_double import CDF

            A = self.change_ring(CDF)
            A = A.conjugate().transpose() * A
            S = A.SVD()[1]
            return max(S.list()).real().sqrt()

        A = self.apply_map(abs, R=RDF)

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

    def numerical_approx(self, prec=None, digits=None, algorithm=None):
        r"""
        Return a numerical approximation of ``self`` with ``prec`` bits
        (or decimal ``digits``) of precision.

        INPUT:

        - ``prec`` -- precision in bits

        - ``digits`` -- precision in decimal digits (only used if
          ``prec`` is not given)

        - ``algorithm`` -- ignored for matrices

        OUTPUT: a matrix converted to a real or complex field

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: d = matrix([[3, 0], [0, sqrt(2)]])
            sage: b = matrix([[1, -1], [2, 2]]); e = b * d * b.inverse(); e
            [ 1/2*sqrt(2) + 3/2 -1/4*sqrt(2) + 3/4]
            [      -sqrt(2) + 3  1/2*sqrt(2) + 3/2]
            sage: e.numerical_approx(53)
            [ 2.20710678118655 0.396446609406726]
            [ 1.58578643762690  2.20710678118655]
            sage: e.numerical_approx(20)
            [ 2.2071 0.39645]
            [ 1.5858  2.2071]
            sage: (e - I).numerical_approx(20)
            [2.2071 - 1.0000*I           0.39645]
            [           1.5858 2.2071 - 1.0000*I]

        ::

            sage: M = matrix(QQ, 4, [i/(i+1) for i in range(12)]); M
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

            sage: matrix(SR, 2, 2, range(4)).n()                                        # needs sage.symbolic
            [0.000000000000000  1.00000000000000]
            [ 2.00000000000000  3.00000000000000]

        ::

            sage: numerical_approx(M)
            [0.000000000000000 0.500000000000000 0.666666666666667]
            [0.750000000000000 0.800000000000000 0.833333333333333]
            [0.857142857142857 0.875000000000000 0.888888888888889]
            [0.900000000000000 0.909090909090909 0.916666666666667]

        We check that :issue:`29700` is fixed::

            sage: M = matrix(3, [1,1,1,1,0,0,0,1,0])
            sage: A, B = M.diagonalization(QQbar)                                       # needs sage.rings.number_field
            sage: _ = A.n()                                                             # needs sage.rings.number_field
        """
        from sage.rings.real_mpfr import RealField
        from sage.rings.complex_mpfr import ComplexField

        if prec is None:
            prec = digits_to_bits(digits)

        try:
            return self.change_ring(RealField(prec))
        except (TypeError, ValueError):
            # try to return a complex result
            return self.change_ring(ComplexField(prec))

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
            sage: A.plot()                                                              # needs sage.plot
            Graphics object consisting of 1 graphics primitive

        Here we make a random matrix over ``RR`` and use ``cmap='hsv'`` to color
        the matrix elements different RGB colors (see documentation for
        ``matrix_plot`` for more information on cmaps)::

            sage: A = random_matrix(RDF, 50)
            sage: plot(A, cmap='hsv')                                                   # needs sage.plot
            Graphics object consisting of 1 graphics primitive

        Another random plot, but over GF(389)::

            sage: A = random_matrix(GF(389), 10)
            sage: A.plot(cmap='Oranges')                                                # needs sage.plot
            Graphics object consisting of 1 graphics primitive
        """
        from sage.plot.matrix_plot import matrix_plot
        return matrix_plot(self, *args, **kwds)

    def derivative(self, *args):
        """
        Derivative with respect to variables supplied in args.

        Multiple variables and iteration counts may be supplied; see
        documentation for the global derivative() function for more
        details.

        EXAMPLES::

            sage: # needs sage.symbolic
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
        from sage.misc.derivative import multi_derivative
        return multi_derivative(self, args)

    def exp(self):
        r"""
        Calculate the exponential of this matrix X, which is the matrix.

        .. MATH::

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

            sage: # needs sage.symbolic
            sage: a = matrix([[1,2], [3,4]])
            sage: a.exp()
            [-1/22*((sqrt(33) - 11)*e^sqrt(33) - sqrt(33) - 11)*e^(-1/2*sqrt(33) + 5/2)              2/33*(sqrt(33)*e^sqrt(33) - sqrt(33))*e^(-1/2*sqrt(33) + 5/2)]
            [             1/11*(sqrt(33)*e^sqrt(33) - sqrt(33))*e^(-1/2*sqrt(33) + 5/2)  1/22*((sqrt(33) + 11)*e^sqrt(33) - sqrt(33) + 11)*e^(-1/2*sqrt(33) + 5/2)]

            sage: type(a.exp())                                                         # needs sage.symbolic
            <class 'sage.matrix.matrix_symbolic_dense.Matrix_symbolic_dense'>

            sage: a = matrix([[1/2,2/3], [3/4,4/5]])
            sage: a.exp()                                                               # needs sage.symbolic
            [-1/418*((3*sqrt(209) - 209)*e^(1/10*sqrt(209)) - 3*sqrt(209) - 209)*e^(-1/20*sqrt(209) + 13/20)                   20/627*(sqrt(209)*e^(1/10*sqrt(209)) - sqrt(209))*e^(-1/20*sqrt(209) + 13/20)]
            [                  15/418*(sqrt(209)*e^(1/10*sqrt(209)) - sqrt(209))*e^(-1/20*sqrt(209) + 13/20)  1/418*((3*sqrt(209) + 209)*e^(1/10*sqrt(209)) - 3*sqrt(209) + 209)*e^(-1/20*sqrt(209) + 13/20)]

            sage: a = matrix(RR, [[1,pi.n()], [1e2,1e-2]])                              # needs sage.symbolic
            sage: a.exp()                                                               # needs sage.symbolic
            [ 1/11882424341266*((11*sqrt(227345670387496707609) + 5941212170633)*e^(3/1275529100*sqrt(227345670387496707609)) - 11*sqrt(227345670387496707609) + 5941212170633)*e^(-3/2551058200*sqrt(227345670387496707609) + 101/200)                            445243650/75781890129165569203*(sqrt(227345670387496707609)*e^(3/1275529100*sqrt(227345670387496707609)) - sqrt(227345670387496707609))*e^(-3/2551058200*sqrt(227345670387496707609) + 101/200)]
            [                                     10000/53470909535697*(sqrt(227345670387496707609)*e^(3/1275529100*sqrt(227345670387496707609)) - sqrt(227345670387496707609))*e^(-3/2551058200*sqrt(227345670387496707609) + 101/200) -1/11882424341266*((11*sqrt(227345670387496707609) - 5941212170633)*e^(3/1275529100*sqrt(227345670387496707609)) - 11*sqrt(227345670387496707609) - 5941212170633)*e^(-3/2551058200*sqrt(227345670387496707609) + 101/200)]
            sage: a.change_ring(RDF).exp()  # rel tol 1e-13                             # needs sage.symbolic
            [42748127.31532951 7368259.244159399]
            [234538976.1381042 40426191.45156228]

        TESTS:

        Sparse matrices are handled correctly (:issue:`28935`), but may
        require a patched version of maxima (:issue:`32898`) for now::

            sage: matrix.diagonal([0], sparse=True).exp()       # not tested            # needs sage.symbolic
            [1]
            sage: matrix.zero(CBF, 2, sparse=True).exp()                                # needs sage.libs.flint sage.symbolic
            [1.000000000000000                 0]
            [                0 1.000000000000000]
        """
        if self.is_sparse():
            # exp is only implemented for dense matrices (:issue:`28935`)
            return self.dense_matrix().exp().sparse_matrix()
        from sage.symbolic.ring import SR
        return self.change_ring(SR).exp()

    def elementary_divisors(self, algorithm=None):
        r"""
        If ``self`` is a matrix over a principal ideal domain `R`, return
        elements `d_i` for `1 \le i \le k = \min(r,s)`
        where `r` and `s` are the number of rows and
        columns of self, such that the cokernel of ``self`` is isomorphic to

        .. MATH::

           R/(d_1) \oplus R/(d_2) \oplus R/(d_k)

        with `d_i \mid d_{i+1}` for all `i`. These are
        the diagonal entries of the Smith form of ``self`` (see
        :meth:`smith_form`).

        INPUT:

        - ``algorithm`` -- ignored

        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: OE.<w> = EquationOrder(x^2 - x + 2)                                   # needs sage.rings.number_field
            sage: m = Matrix([[1, w], [w, 7]])                                          # needs sage.rings.number_field
            sage: m.elementary_divisors()                                               # needs sage.rings.number_field
            [1, -w + 9]

        .. SEEALSO::

           :meth:`smith_form`
        """
        d = self.smith_form(transformation=False)
        r = min(self.nrows(), self.ncols())
        return [d[i, i] for i in range(r)]

    def smith_form(self, transformation=True, integral=None, exact=True):
        r"""
        Return a Smith normal form of this matrix.

        For a matrix `M`, a Smith normal form is a matrix `S = UMV` such that:

        * `U` and `V` are invertible matrices
        * the only non-vanishing entries of `S` are located on the diagonal
          (though `S` might not be a square matrix)
        * if `d_i` denotes the entry of `S` at `(i,i)`, then `d_i` divides
          `d_{i+1}` for all `i`, i.e., the `d_i` are the ordered
          :meth:`elementary_divisors` of `M`

        Note that the matrices `U` and `V` are not uniquely determined and the
        `d_i` are only uniquely determined up to units. For some base rings,
        such as local rings, the `d_i` might be further normalized, see
        ``LOCAL RINGS`` below.

        If the base ring is not a PID, the routine might work, or else it will
        fail having found an example of a non-principal ideal. Note that we do
        not call any methods to check whether or not the base ring is a PID,
        since this might be quite expensive (e.g. for rings of integers of
        number fields of large degree).

        INPUT:

        - ``transformation`` -- boolean (default: ``True``); whether the
          matrices `U` and `V` should be returned

        - ``integral`` -- a subring of the base ring, boolean or ``None``
          (default: ``None``); the entries of `U` and `V` are taken
          from this subring. If ``True``, the ring is taken to be the
          ring of integers of the base ring; if ``False`` the fraction field
          of the base ring; if ``None`` the base ring itself.
          When a subring is specified, multiplying
          by the denominator must map the entries into the subring; in this
          case the transformation matrices will have entries in this subring.

        - ``exact`` -- boolean (default: ``True``); only used for local rings/fields.
          See ``LOCAL RINGS`` for more details.

        OUTPUT:

        The matrices `S, U, V` or the matrix `S` depending on
        ``transformation``.

        ALGORITHM:

        If the base ring has a method ``_matrix_smith_form``, use it; note that
        ``_matrix_smith_form`` might choose to further normalize the output.

        Otherwise, use the algorithm from :wikipedia:`Smith_normal_form`.

        LOCAL RINGS:

        Over local rings, we normalize `S` to only contain powers of the uniformizer.

        In order to simplify the precision handling, we truncate the absolute precision
        of the input matrix to the minimum absolute precision of any of its entries.
        As long as all of the elementary divisors are nonzero modulo this precision,
        they can be determined exactly since they are defined to be powers of the
        uniformizer.  In this case, which is specified by the keyword ``exact=True``,
        one of the transformation matrices will be inexact: `U` in the case that
        the number of rows is at least the number of columns, and `V` otherwise.

        If ``exact=False``, we instead return an inexact Smith form.  Now the
        transformation matrices are exact and we can deal gracefully with
        elementary divisors that are zero modulo the working precision.  However,
        the off-diagonal entries of the smith form carry a precision that
        can affect the precision of future calculations.

        See ``_matrix_smith_form`` on the base ring for more detail.

        EXAMPLES:

        An example over the ring of integers of a number field (of class
        number 1)::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: OE.<w> = EquationOrder(x^2 - x + 2)
            sage: m = Matrix([[1, w], [w, 7]])
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

            sage: R.<x> = QQ[]; m = x*matrix(R, 2, 2, 1) - matrix(R, 2, 2, [3,-4,1,-1])
            sage: m.smith_form()
            (
            [            1             0]  [    0    -1]  [    1 x + 1]
            [            0 x^2 - 2*x + 1], [    1 x - 3], [    0     1]
            )

        An example over a field::

            sage: m = matrix(GF(17), 3, 3, [11,5,1, 3,6,8, 1,16,0])
            sage: d, u, v = m.smith_form()
            sage: d
            [1 0 0]
            [0 1 0]
            [0 0 0]
            sage: u * m * v == d
            True

        When the base ring has a ``ring_of_integers`` method and supports denominators,
        you can get an integral version of the smith form::

            sage: m = matrix(QQ, 2, 2, [17/6, 47/6, 25/6, 23/2])
            sage: m.smith_form()
            (
            [1 0]  [6/17    0]  [     1 -47/17]
            [0 1], [  75  -51], [     0      1]
            )
            sage: m.smith_form(integral=True)
            (
            [1/6   0]  [  3  -2]  [ 1  3]
            [  0 1/3], [-25  17], [ 0 -1]
            )

        Some examples over non-PID's work anyway::

            sage: # needs sage.rings.number_field
            sage: R.<s> = EquationOrder(x^2 + 5)  # class number 2
            sage: A = matrix(R, 2, 2, [s - 1, -s, -s, 2*s + 1])
            sage: D, U, V = A.smith_form()
            sage: D, U, V
            (
            [     1      0]  [    4 s + 4]  [       1 -5*s + 6]
            [     0 -s - 6], [    s s - 1], [       0        1]
            )
            sage: D == U * A * V
            True

        Others don't, but they fail quite constructively::

            sage: matrix(R, 2, 2, [s - 1, -s - 2, -2*s, -s - 2]).smith_form()           # needs sage.rings.number_field
            Traceback (most recent call last):
            ...
            ArithmeticError: Ideal Fractional ideal (2, s + 1) not principal

        Empty matrices are handled safely::

            sage: # needs sage.rings.number_field
            sage: m = MatrixSpace(OE, 2,0)(0)
            sage: d, u, v = m.smith_form(); u * m * v == d
            True
            sage: m = MatrixSpace(OE, 0,2)(0)
            sage: d, u, v = m.smith_form(); u * m * v == d
            True
            sage: m = MatrixSpace(OE, 0,0)(0)
            sage: d, u, v = m.smith_form(); u * m * v == d
            True

        Some pathological cases that crashed earlier versions::

            sage: # needs sage.rings.number_field
            sage: m = Matrix(OE, [[ 2*w,      2*w - 1,  -w + 1],
            ....:                 [ 2*w + 2, -2*w - 1,   w - 1],
            ....:                 [-2*w - 1, -2*w - 2, 2*w - 1]])
            sage: d, u, v = m.smith_form(); u * m * v == d
            True
            sage: m = matrix(OE, [[-5*w - 1, -2*w - 2, 4*w - 10],
            ....:                 [ 8*w,       -w,       w - 1],
            ....:                 [      -1,        1,      -8]])
            sage: d, u, v = m.smith_form(); u * m * v == d
            True

        Over local fields, we can request the transformation matrices to be integral:;

            sage: K = Qp(2, 5, print_mode='terse')                                      # needs sage.rings.padics
            sage: M = matrix(K, 2, 3, [1/2, 1, 2, 1/3, 1, 3])                           # needs sage.rings.padics
            sage: M.smith_form(integral=True)                                           # needs sage.rings.padics
            (
            [1/2 + O(2^4)            0            0]  [ 1 + O(2^5)           0]
            [           0   1 + O(2^5)            0], [42 + O(2^6)  1 + O(2^5)],
            <BLANKLINE>
            [ 1 + O(2^5) 26 + O(2^5)  6 + O(2^5)]
            [     O(2^4)  3 + O(2^4) 11 + O(2^4)]
            [          0           0  1 + O(2^5)]
            )
        """
        R = self.base_ring()
        if hasattr(R, '_matrix_smith_form'):
            return R._matrix_smith_form(self, transformation=transformation, integral=integral, exact=exact)
        if integral is True:
            integral = R.ring_of_integers()
        elif integral is R:
            integral = None
        if integral is False:
            self = self.change_ring(R.fraction_field())
        elif integral is not None:
            # Try to clear denominators
            den = self.denominator()
            self = (den * self).change_ring(integral)
        if transformation:
            left_mat = self.new_matrix(self.nrows(), self.nrows(), 1)
            right_mat = self.new_matrix(self.ncols(), self.ncols(), 1)
        if self == 0 or (self.nrows() <= 1 and self.ncols() <= 1):
            if transformation:
                return self.__copy__(), left_mat, right_mat
            else:
                return self.__copy__()

        # data type checks on R
        if not R.is_integral_domain() or not R.is_noetherian():
            raise TypeError("Smith form only defined over Noetherian integral domains")
        if not R.is_exact():
            raise NotImplementedError("Smith form over non-exact rings not implemented at present")

        # first clear the first row and column
        u, t, v = _smith_onestep(self)

        # now recurse: t now has a nonzero entry at 0,0 and zero entries in the rest
        # of the 0th row and column, so we apply smith_form to the smaller submatrix
        mm = t.submatrix(1, 1)
        if transformation:
            dd, uu, vv = mm.smith_form(transformation=True)
        else:
            dd = mm.smith_form(transformation=False)
        d = dd.new_matrix(1, 1, [t[0, 0]]).block_sum(dd)
        if transformation:
            u = uu.new_matrix(1, 1, [1]).block_sum(uu) * u
            v = v * vv.new_matrix(1, 1, [1]).block_sum(vv)
        dp, up, vp = _smith_diag(d, transformation=transformation)
        if integral is False:
            dp = dp.change_ring(R)
        elif integral is not None:
            dp = dp.change_ring(R) / den
        if transformation:
            return dp, up*u, v*vp
        else:
            return dp

    def fitting_ideal(self, i):
        r"""
        Return the `i`-th Fitting ideal of the matrix. This is the ideal generated
        by the `n - i` minors, where `n` is the number of columns.

        INPUT:

        - ``i`` -- integer

        OUTPUT: an ideal on the base ring

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: M = matrix(R, [[2*x-z, 0, y-z^2, 1], [0, z - y, z - x, 0],[z - y, x^2 - y, 0, 0]])
            sage: M
            [ 2*x - z        0 -z^2 + y        1]
            [       0   -y + z   -x + z        0]
            [  -y + z  x^2 - y        0        0]
            sage: [R.ideal(M.minors(i)) == M.fitting_ideal(4-i) for i in range(5)]
            [True, True, True, True, True]
            sage: M.fitting_ideal(0)
            Ideal (0) of Multivariate Polynomial Ring in x, y, z over Rational Field
            sage: M.fitting_ideal(1)
            Ideal (2*x^4 - 3*x^3*z + x^2*z^2 + y^2*z^2 - 2*y*z^3 + z^4 - 2*x^2*y - y^3 + 3*x*y*z + 2*y^2*z - 2*y*z^2, -x^3 + x^2*z + x*y - y*z, y^2 - 2*y*z + z^2, x*y - x*z - y*z + z^2) of Multivariate Polynomial Ring in x, y, z over Rational Field
            sage: M.fitting_ideal(3)
            Ideal (2*x - z, -z^2 + y, 1, -y + z, -x + z, -y + z, x^2 - y) of Multivariate Polynomial Ring in x, y, z over Rational Field
            sage: M.fitting_ideal(4)
            Ideal (1) of Multivariate Polynomial Ring in x, y, z over Rational Field


        If the base ring is a field, the Fitting ideals are zero under the corank::

            sage: M = matrix(QQ, [[2,1,3,5],[4,2,6,6],[0,3,2,0]])
            sage: M
            [2 1 3 5]
            [4 2 6 6]
            [0 3 2 0]
            sage: M.fitting_ideal(0)
            Principal ideal (0) of Rational Field
            sage: M.fitting_ideal(1)
            Principal ideal (1) of Rational Field
            sage: M.fitting_ideal(2)
            Principal ideal (1) of Rational Field
            sage: M.fitting_ideal(3)
            Principal ideal (1) of Rational Field
            sage: M.fitting_ideal(4)
            Principal ideal (1) of Rational Field


        In the case of principal ideal domains, it is given by the elementary
        divisors::

            sage: M = matrix([[2,1,3,5],[4,2,6,6],[0,3,2,0]])
            sage: M
            [2 1 3 5]
            [4 2 6 6]
            [0 3 2 0]
            sage: M.fitting_ideal(0)
            Principal ideal (0) of Integer Ring
            sage: M.fitting_ideal(1)
            Principal ideal (4) of Integer Ring
            sage: M.fitting_ideal(2)
            Principal ideal (1) of Integer Ring
            sage: M.fitting_ideal(3)
            Principal ideal (1) of Integer Ring
            sage: M.fitting_ideal(4)
            Principal ideal (1) of Integer Ring
            sage: M.elementary_divisors()
            [1, 1, 4]

        This is also true for univariate polynomials over a field::

            sage: R.<x> = QQ[]
            sage: M = matrix(R,[[x^2-2*x+1, x-1,x^2-1],[0,x+1,1]])
            sage: M.fitting_ideal(0)
            Principal ideal (0) of Univariate Polynomial Ring in x over Rational Field
            sage: M.fitting_ideal(1)
            Principal ideal (x - 1) of Univariate Polynomial Ring in x over Rational Field
            sage: M.fitting_ideal(2)
            Principal ideal (1) of Univariate Polynomial Ring in x over Rational Field
            sage: M.smith_form()[0]
            [    1     0     0]
            [    0 x - 1     0]
        """
        R = self.base_ring()
        if not R.is_exact():
            raise NotImplementedError("Fitting ideals over non-exact rings not implemented at present")
        n = self.ncols()
        rank_minors = n - i
        if rank_minors > self.nrows():
            return R.ideal([R.zero()])
        elif rank_minors <= 0:
            return R.ideal([R.one()])
        elif rank_minors == 1:
            return R.ideal(self.coefficients())
        if R in _Fields:
            if self.rank() >= rank_minors:
                return R.ideal([1])
            else:
                return R.ideal([0])
        try:
            elemdiv = self.elementary_divisors()
            if rank_minors > len(elemdiv):
                return R.ideal([0])
            return R.ideal(prod(elemdiv[:rank_minors]))
        except (TypeError, NotImplementedError, ArithmeticError):
            pass
        for (nr, r) in enumerate(self.rows()):
            nz = [e for e in enumerate(r) if e[1]]
            if len(nz) == 0:
                N = self.delete_rows([nr])
                return N.fitting_ideal(i)
            elif len(nz) == 1:
                N = self.delete_rows([nr])
                F1 = N.fitting_ideal(i)
                N = N.delete_columns([nz[0][0]])
                F2 = N.fitting_ideal(i)
                return F1 + nz[0][1]*F2
        for (nc, c) in enumerate(self.columns()):
            nz = [e for e in enumerate(c) if e[1]]
            if len(nz) == 0:
                N = self.delete_columns([nc])
                return N.fitting_ideal(i - 1)
            elif len(nz) == 1:
                N = self.delete_columns([nc])
                F1 = N.fitting_ideal(i-1)
                N = N.delete_rows([nz[0][0]])
                F2 = N.fitting_ideal(i)
                return F1 + nz[0][1]*F2
        if hasattr(self, '_fitting_ideal'):
            try:
                return self._fitting_ideal(i)
            except NotImplementedError:
                pass
        return R.ideal(self.minors(rank_minors))

    def _hermite_form_euclidean(self, transformation=False, normalization=None):
        """
        Transform the matrix in place to hermite normal form and optionally
        return the transformation matrix.

        The matrix is assumed to be over a Euclidean domain. In particular,
        ``xgcd()`` method should be available for the elements of the domain.

        INPUT:

        - ``transformation`` -- boolean (default: ``False``); if ``True``,
          return the transformation matrix

        - ``normalization`` -- function (default: ``None``); if given, the
          function is applied to each pivot to get a normalization coefficient,
          which is multiplied to the pivot.

        EXAMPLES::

            sage: B = matrix(ZZ, 3, [-1,-2,-3,4,5,6,7,8,9]); B
            [-1 -2 -3]
            [ 4  5  6]
            [ 7  8  9]
            sage: C = B.__copy__()
            sage: U = C._hermite_form_euclidean(transformation=True)
            sage: C
            [1 2 3]
            [0 3 6]
            [0 0 0]
            sage: U
            [-1  0  0]
            [-4 -1  0]
            [-1 -2  1]
            sage: U * B == C
            True

            sage: P.<x> = PolynomialRing(QQ)
            sage: A = matrix(P,3,[-(x-1)^((i-j) % 3) for i in range(3) for j in range(3)])
            sage: A
            [            -1 -x^2 + 2*x - 1         -x + 1]
            [        -x + 1             -1 -x^2 + 2*x - 1]
            [-x^2 + 2*x - 1         -x + 1             -1]
            sage: H = A.__copy__()
            sage: U = H._hermite_form_euclidean(transformation=True,
            ....:                               normalization=lambda p: ~p.lc())
            sage: H
            [                    1         x^2 - 2*x + 1                 x - 1]
            [                    0 x^3 - 3*x^2 + 3*x - 2                     0]
            [                    0                     0 x^3 - 3*x^2 + 3*x - 2]
            sage: U * A == H
            True
        """
        cdef Matrix A = self
        cdef Matrix U

        cdef Py_ssize_t m = A.nrows()
        cdef Py_ssize_t n = A.ncols()

        cdef Py_ssize_t i = 0
        cdef Py_ssize_t j = 0

        cdef Py_ssize_t k, l, c

        if transformation:
            from sage.matrix.constructor import identity_matrix
            U = identity_matrix(A.base_ring(), m)

        pivot_cols = []
        while j < n:
            k = i
            while k < m and A.get_unsafe(k, j).is_zero():  # first nonzero entry
                k += 1
            if k < m:
                l = k + 1
                while l < m:
                    while l < m and A.get_unsafe(l, j).is_zero():  # nonzero entry below
                        l += 1
                    if l >= m: break

                    a = A.get_unsafe(k, j)
                    b = A.get_unsafe(l, j)
                    d, p, q = a.xgcd(b)     # p * a + q * b = d = gcd(a,b)
                    e = a // d
                    f = b // d

                    for c in range(j, n):
                        Akc = A.get_unsafe(k, c)
                        Alc = A.get_unsafe(l, c)
                        A.set_unsafe(k, c, p * Akc + q * Alc)
                        A.set_unsafe(l, c, (-f) * Akc + e * Alc)
                    if transformation:
                        for c in range(m):
                            Ukc = U.get_unsafe(k, c)
                            Ulc = U.get_unsafe(l, c)
                            U.set_unsafe(k, c, p * Ukc + q * Ulc)
                            U.set_unsafe(l, c, (-f) * Ukc + e * Ulc)
                if i != k:
                    A.swap_rows_c(i, k)

                    if transformation:
                        U.swap_rows_c(i, k)
                pivot_cols.append(j)
                i += 1
            j += 1

        # reduce entries above pivots
        for i in range(len(pivot_cols)):
            j = pivot_cols[i]
            pivot = A.get_unsafe(i, j)

            # possibly normalize the pivot
            if normalization:
                coeff = normalization(pivot)
                for c in range(j, n):
                    A.set_unsafe(i, c, A.get_unsafe(i, c) * coeff)
                if transformation:
                    for c in range(m):
                        U.set_unsafe(i, c, U.get_unsafe(i, c) * coeff)

            pivot = A.get_unsafe(i, j)
            for k in range(i):
                q = - (A.get_unsafe(k, j) // pivot)
                if not q.is_zero():
                    for c in range(j, n):
                        A.set_unsafe(k, c, A.get_unsafe(k, c) + q * A.get_unsafe(i, c))
                    if transformation:
                        for c in range(m):
                            U.set_unsafe(k, c, U.get_unsafe(k, c) + q * U.get_unsafe(i, c))

        if transformation:
            return U

    def hermite_form(self, include_zero_rows=True, transformation=False):
        """
        Return the Hermite form of self, if it is defined.

        INPUT:

        - ``include_zero_rows`` -- boolean (default: ``True``); if ``False``
          the zero rows in the output matrix are deleted

        - ``transformation`` -- boolean (default: ``False``); a matrix `U` such
          that ``U*self == H``

        OUTPUT:

        - matrix H
        - (optional) transformation matrix `U` such that ``U*self == H``,
          possibly with zero rows deleted

        EXAMPLES::

            sage: M = FunctionField(GF(7), 'x').maximal_order()
            sage: K.<x> = FunctionField(GF(7)); M = K.maximal_order()
            sage: A = matrix(M, 2, 3, [x, 1, 2*x, x, 1 + x, 2])
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
            sage: H, U = A.hermite_form(transformation=True, include_zero_rows=True)
            sage: H, U
            (
            [  x   1 2*x]  [1 0]
            [  0   0   0], [5 1]
            )
            sage: U * A == H
            True
            sage: H, U = A.hermite_form(transformation=True, include_zero_rows=False)
            sage: U * A
            [  x   1 2*x]
            sage: U * A == H
            True
        """
        left, H, _ = self._echelon_form_PID()
        if not include_zero_rows:
            i = H.nrows() - 1
            while H.row(i) == 0:
                i -= 1
            H = H[:i + 1]
            if transformation:
                left = left[:i + 1]
        if transformation:
            return H, left
        return H

    def _echelon_form_PID(self):
        r"""
        Return a triple (left, a, pivots) where left*self == a and a is row
        echelon form (which in this case we define to mean Hermite normal
        form).

        When ideals of the base ring have a "small_residue" method (as is the
        case for number field ideals), we use this to reduce entries above each
        column pivot.
        EXAMPLES::

            sage: x = polygen(ZZ, 'x')
            sage: L.<a> = NumberField(x^3 - 2)                                          # needs sage.rings.number_field
            sage: OL = L.ring_of_integers()                                             # needs sage.rings.number_field

        We check some degenerate cases::

            sage: # needs sage.rings.number_field
            sage: m = matrix(OL, 0, 0, []); r,s,p = m._echelon_form_PID()
            sage: (r,s,p)
            ([], [], [])
            sage: r * m == s and r.det() == 1
            True
            sage: m = matrix(OL, 0, 1, []); r,s,p = m._echelon_form_PID()
            sage: (r,s,p)
            ([], [], [])
            sage: r * m == s and r.det() == 1
            True
            sage: m = matrix(OL, 1, 0, []); r,s,p = m._echelon_form_PID()
            sage: (r,s,p)
            ([1], [], [])
            sage: r * m == s and r.det() == 1
            True

        A 2x2 matrix::

            sage: m = matrix(OL, 2, 2, [1, 0, a, 2])                                    # needs sage.rings.number_field
            sage: r,s,p = m._echelon_form_PID(); (r,s,p)                                # needs sage.rings.number_field
            (
            [ 1  0]  [1 0]
            [-a  1], [0 2], [0, 1]
            )
            sage: r * m == s and r.det() == 1                                           # needs sage.rings.number_field
            True

        A larger example::

            sage: # needs sage.rings.number_field
            sage: m = matrix(OL, 3, 5, [a^2 - 3*a - 1, a^2 - 3*a + 1, a^2 + 1,
            ....:   -a^2 + 2, -3*a^2 - a - 1, -6*a - 1, a^2 - 3*a - 1,
            ....:   2*a^2 + a + 5, -2*a^2 + 5*a + 1, -a^2 + 13*a - 3,
            ....:   -2*a^2 + 4*a - 2, -2*a^2 + 1, 2*a, a^2 - 6, 3*a^2 - a ])
            sage: r,s,p = m._echelon_form_PID()
            sage: s[2]
            (0, 0, 3*a^2 + 18*a - 34, 68*a^2 - 134*a + 53, 111*a^2 - 275*a + 90)
            sage: r * m == s and r.det() == 1
            True

        We verify that :issue:`9053` is resolved::

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
                    self.nonzero_positions_in_row(0)[0]]

        R = self.base_ring()

        # data type checks on R
        if not R.is_integral_domain():
            raise TypeError("Generic echelon form only defined over "
                            "integral domains")
        if not R.is_exact():
            raise NotImplementedError("Echelon form over generic non-exact "
                                      "rings not implemented at present")

        left_mat, a = _generic_clear_column(self)
        assert left_mat * self == a

        if a[0, 0] != 0:
            aa = a.submatrix(1, 1)
            s, t, pivs = aa._echelon_form_PID()
            left_mat = s.new_matrix(1, 1, [1]).block_sum(s) * left_mat
            a = left_mat * self
            pivs = [0] + [x + 1 for x in pivs]

        else:
            aa = a.submatrix(0, 1)
            s, t, pivs = aa._echelon_form_PID()
            left_mat = s * left_mat
            a = left_mat * self
            pivs = [x+1 for x in pivs]

        try:
            for i in range(1, len(pivs)):
                y = a[i][pivs[i]]
                I = ideal_or_fractional(R, y)
                s = a[0][pivs[i]]
                t = I.small_residue(s)
                v = R((s-t) / y)

                left_mat.add_multiple_of_row(0, i, -v)
                a.add_multiple_of_row(0, i, -v)
                assert left_mat * self == a
        except AttributeError:  # on I.small_residue
            pass

        return left_mat, a, pivs

    def _zigzag_form(self, bint basis=True):
        r"""
        Helper method for computation of ZigZag form.

        INPUT:

        - ``self`` -- a square matrix over an exact field

        - ``basis`` -- (default; ``True``); controls whether or not to
          compute a change-of-basis matrix (also called a transformation
          matrix)

        OUTPUT:

        See the documentation for the :meth:`zigzag_form` method for a
        description of ZigZag form and notes on the algorithm employed.

        When ``basis=True`` four items are returned.  The first is a matrix
        whose columns form a basis so that a matrix representation of ``self``
        relative to this basis will be the ZigZag form.  More precisely, if
        `A` is the matrix, `U` is the change-of-basis matrix and `Z` is
        the ZigZag form, then

        .. MATH::

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

        The computation of the change-of-basis matrix has not been optimized.
        As a helper method, no error checking is performed on the inputs -
        that should be performed by the calling method.

        ALGORITHM:

        ZigZag form, and its computation, are due to Arne Storjohann
        and are  described in [Sto2000]_ and
        [Sto1998]_, where the former is more
        representative of the code here.

        EXAMPLES::

            sage: A = matrix(QQ, [[-68,   69, -27, -11, -65,   9, -181, -32],
            ....:                 [-52,   52, -27,  -8, -52, -16, -133, -14],
            ....:                 [ 92,  -97,  47,  14,  90,  32,  241,  18],
            ....:                 [139, -144,  60,  18, 148, -10,  362,  77],
            ....:                 [ 40,  -41,  12,   6,  45, -24,  105,  42],
            ....:                 [-46,   48, -20,  -7, -47,   0, -122, -22],
            ....:                 [-26,   27, -13,  -4, -29,  -6,  -66, -14],
            ....:                 [-33,   34, -13,  -5, -35,   7,  -87, -23]])
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
        """
        cdef Py_ssize_t n, s, c, i, j, k

        n = self._ncols
        R = self.base_ring()
        zero = R.zero()
        one = R.one()
        cdef Matrix Z = self.__copy__()
        cdef list polys = []    # coefficients for polynomials of companion matrices
        cdef list corners = []  # zero or one in corner of off-diagonal blocks
        if basis:
            from sage.matrix.constructor import identity_matrix
            U = identity_matrix(R, n)  # transformation matrix
        # parity switch, True iff working on transpose
        # if False, mimic row operations only on U
        # if True,  mimic column operations only on U
        cdef bint trans = False, zigging
        s = 0  # index of top row of current block
        c = 0  # index of current row of current block
        while s < n:
            zigging = True
            # check for totally zero column below diagonal
            while zigging:  # zigging means we are building a block
                nonzero = -1
                for i in range(c+1, n):
                    if Z.get_unsafe(i, c):
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
                    scale = Z.get_unsafe(c+1, c)
                    Z.rescale_row(c+1, 1/scale)
                    Z.rescale_col(c+1, scale)
                    if basis:
                        if trans:
                            U.rescale_col(c+1, scale)
                        else:
                            U.rescale_row(c+1, ~scale)

                    # clear column throughout the block,and in all rows below
                    for i in range(s, n):
                        if i != c+1:
                            scale = Z.get_unsafe(i, c)
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
            p = [-Z.get_unsafe(i, c) for i in range(s, c+1)]
            p.append(one)
            polys.append(p)

            # use new unit columns (i) to clear rows to right (j)
            # all but top row of block to the right will become zero
            # it is important that the loop on  i  goes in reverse
            for i in range(c-1, s-1, -1):
                for j in range(c+1, n):
                    scale = Z.get_unsafe(i+1, j)
                    # Effectively: Z.add_multiple_of_column(j, i, -scale)
                    Z.set_unsafe(i+1, j, zero)
                    # row j leads with zeros, up to, and including column c
                    # Effectively: Z.add_multiple_of_row(i, j, scale)
                    for k in range(c+1, n):
                        # Z[i,k] = Z[i,k] + scale*Z[j,k]
                        Z.set_unsafe(i, k, Z.get_unsafe(i, k)+scale*Z.get_unsafe(j, k))
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
                if Z.get_unsafe(s, j):
                    nonzero = j
                    break
            if (nonzero != -1):
                # swap column with nonzero entry just outside block
                if nonzero != c+1:
                    Z.swap_columns(c+1, nonzero)
                    Z.swap_rows(c+1, nonzero)
                    if basis:
                        if trans:
                            U.swap_columns(c+1, nonzero)
                        else:
                            U.swap_rows(c+1, nonzero)
                scale = Z.get_unsafe(s, c+1)
                Z.rescale_col(c+1, ~scale)
                Z.rescale_row(c+1, scale)
                if basis:
                    if trans:
                        U.rescale_col(c+1, ~scale)
                    else:
                        U.rescale_row(c+1, scale)
                for j in range(c+2, n):
                    scale = Z.get_unsafe(s, j)
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
                corners.append(Z.get_unsafe(s, c+1))
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

        - ``self`` -- a square matrix with entries from an exact field

        - ``transformation`` -- (default: ``False``) if ``True`` return a
          change-of-basis matrix relating the matrix and its ZigZag form

        - ``subdivide`` -- (default: ``True``) if ``True`` the ZigZag
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
        words, the representation of ``self`` with respect to the columns
        of ``U`` will be ``Z``.

        If subdivide is ``True`` then the matrix returned as the form is
        partitioned according to the companion matrices and these may be
        manipulated by several different matrix methods.

        For output that may be more useful as input to other routines,
        see the helper method :meth:`_zigzag_form`.

        .. NOTE::

            An effort has been made to optimize computation of the form,
            but no such work has been done for the computation of the
            transformation matrix, so for fastest results do not request the
            transformation matrix.

        ALGORITHM:

        ZigZag form, and its computation, are due to Arne Storjohann
        and are  described in [Sto2000]_ and
        [Sto1998]_, where the former is more
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
            ....:                 [-52,   52, -27,  -8, -52, -16, -133, -14],
            ....:                 [ 92,  -97,  47,  14,  90,  32,  241,  18],
            ....:                 [139, -144,  60,  18, 148, -10,  362,  77],
            ....:                 [ 40,  -41,  12,   6,  45, -24,  105,  42],
            ....:                 [-46,   48, -20,  -7, -47,   0, -122, -22],
            ....:                 [-26,   27, -13,  -4, -29,  -6,  -66, -14],
            ....:                 [-33,   34, -13,  -5, -35,   7,  -87, -23]])
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
            ....:                 [ 26,  54,   6,  -5, -28,   73,   73, -48],
            ....:                 [-16, -79,  12, -10,  64, -142, -115,  41],
            ....:                 [ 27,  -7,  21, -33,  39,  -20,  -42,  43],
            ....:                 [  8, -75,  34, -32,  86, -156, -130,  42],
            ....:                 [  2, -17,   7,  -8,  20,  -33,  -31,  16],
            ....:                 [-24, -80,   7,  -3,  56, -136, -112,  42],
            ....:                 [ -6, -19,   0,  -1,  13,  -28,  -27,  15]])
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

            sage: A.jordan_form() == B.jordan_form()                                    # needs sage.combinat sage.libs.pari
            True

        Two more examples, illustrating the two extremes of the zig-zag
        nature of this form.  The first has a one in each of the
        off-diagonal blocks, the second has all zeros in each off-diagonal
        block.  Notice again that the two matrices are similar, since their
        Jordan canonical forms are equal. ::

            sage: C = matrix(QQ, [[2,  31, -10,  -9, -125,  13,  62, -12],
            ....:                 [0,  48, -16, -16, -188,  20,  92, -16],
            ....:                 [0,   9,  -1,   2,  -33,   5,  18,   0],
            ....:                 [0,  15,  -5,   0,  -59,   7,  30,  -4],
            ....:                 [0, -21,   7,   2,   84, -10, -42,   5],
            ....:                 [0, -42,  14,   8,  167, -17, -84,  13],
            ....:                 [0, -50,  17,  10,  199, -23, -98,  14],
            ....:                 [0,  15,  -5,  -2,  -59,   7,  30, -2]])
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
            ....:                 [ -6,   5,    7,   2,  -4,   5,    7,   -3],
            ....:                 [ 21, -12,   89,  25,   8,  27,   98,  -95],
            ....:                 [ -9,   5,  -44, -11,  -3, -13,  -48,   47],
            ....:                 [ 23, -13,   74,  21,  12,  22,   85,  -84],
            ....:                 [ 31, -18,  135,  38,  12,  47,  155, -147],
            ....:                 [-33,  19, -138, -39, -13, -45, -156,  151],
            ....:                 [ -7,   4,  -29,  -8,  -3, -10,  -34,  34]])
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

            sage: C.jordan_form() == D.jordan_form()                                    # needs sage.combinat sage.libs.pari
            True

        ZigZag form is achieved entirely with the operations of the field, so
        while the eigenvalues may lie outside the field, this does not impede
        the computation of the form.  ::

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = GF(5^4)
            sage: A = matrix(F, [[     a,      0,  0, a + 3],
            ....:                [     0,a^2 + 1,  0,     0],
            ....:                [     0,      0,a^3,     0],
            ....:                [a^2 +4 ,     0,   0,a + 2]])
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
            TypeError: no canonical coercion from Finite Field in a of size 5^4 to Finite Field in z4 of size 5^4

        Subdivisions are optional.  ::

            sage: F.<a> = GF(5^4)                                                       # needs sage.rings.finite_rings
            sage: A = matrix(F, [[     a,      0,  0, a + 3],                           # needs sage.rings.finite_rings
            ....:                [     0,a^2 + 1,  0,     0],
            ....:                [     0,      0,a^3,     0],
            ....:                [a^2 +4 ,     0,   0,a + 2]])
            sage: A.zigzag_form(subdivide=False)                                        # needs sage.rings.finite_rings
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
        """
        R = self.base_ring()
        if not self.is_square():
            raise TypeError("matrix must be square, not {0} x {1}".format(self.nrows(), self.ncols()))
        if not (R in _Fields and R.is_exact()):
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
        Return the rational canonical form, also known as Frobenius form.

        INPUT:

        - ``self`` -- a square matrix with entries from an exact field

        - ``format`` -- (default: ``'right'``) one of ``'right'``, ``'bottom'``,
          ``'left'``, ``'top'`` or ``'invariants'``.  The first four will cause
          a matrix to be returned with companion matrices dictated by the
          keyword.  The value 'invariants' will cause a list of lists to
          be returned, where each list contains coefficients of a
          polynomial associated with a companion matrix.

        - ``subdivide`` -- (default: ``'True'``) if 'True' and a matrix is
          returned, then it contains subdivisions delineating the
          companion matrices along the diagonal

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
        algorithm of Storjohann's [Sto2011]_.

        EXAMPLES:

        The lists of coefficients returned with the ``invariants`` keyword
        are designed to easily convert to the polynomials associated with
        the companion matrices.  This is illustrated by the construction
        below of the ``polys`` list.  Then we can test the divisibility
        condition on the list of polynomials. Also the minimal and
        characteristic polynomials are easy to determine from this list. ::

            sage: A = matrix(QQ, [[ 11,  14, -15,  -4, -38, -29,  1,  23,  14, -63,  17,  24,  36,  32],
            ....:                 [ 18,   6, -17, -11, -31, -43, 12,  26,   0, -69,  11,  13,  17,  24],
            ....:                 [ 11,  16, -22,  -8, -48, -34,  0,  31,  16, -82,  26,  31,  39,  37],
            ....:                 [ -8, -18,  22,  10,  46,  33,  3, -27, -12,  70, -19, -20, -42, -31],
            ....:                 [-13, -21,  16,  10,  52,  43,  4, -28, -25,  89, -37, -20, -53, -62],
            ....:                 [ -2,  -6,   0,   0,   6,  10,  1,   1,  -7,  14, -11,  -3, -10, -18],
            ....:                 [ -9, -19,  -3,   4,  23,  30,  8,  -3, -27,  55, -40,  -5, -40, -69],
            ....:                 [  4,  -8,  -1,  -1,   5,  -4,  9,   5, -11,   4, -14,  -2, -13, -17],
            ....:                 [  1,  -2,  16,  -1,  19,  -2, -1, -17,   2,  19,   5, -25,  -7,  14],
            ....:                 [  7,   7, -13,  -4, -26,  -21, 3,  18,   5, -40,   7,  15,  20,  14],
            ....:                 [ -6,  -7, -12,   4,  -1,  18,  3,   8, -11,  15, -18,  17, -15, -41],
            ....:                 [  5,  11, -11,  -3, -26, -19, -1,  14,  10, -42,  14,  17,  25,  23],
            ....:                 [-16, -15,   3,  10,  29,  45, -1, -13, -19,  71, -35,  -2, -35, -65],
            ....:                 [  4,   2,   3,  -2,  -2, -10,  1,   0,   3, -11,   6,  -4,   6,  17]])
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
            sage: [p.factor() for p in polys]                                           # needs sage.rings.finite_rings
            [(x - 2)^2, (x - 3) * (x + 1) * (x - 2)^2, (x + 1)^2 * (x - 3)^3 * (x - 2)^3]
            sage: all(polys[i].divides(polys[i+1]) for i in range(len(polys)-1))
            True
            sage: polys[-1] == A.minimal_polynomial(var='x')                            # needs sage.libs.pari
            True
            sage: prod(polys) == A.characteristic_polynomial(var='x')                   # needs sage.libs.pari
            True

        Rational form is a canonical form.  Any two matrices are similar
        if and only if their rational forms are equal.  By starting with
        Jordan canonical forms, the matrices ``C`` and ``D`` below were
        built as similar matrices, while ``E`` was built to be just
        slightly different.  All three matrices have equal characteristic
        polynomials though ``E``'s minimal polynomial differs. ::

            sage: C = matrix(QQ, [[2,  31, -10,  -9, -125,  13,  62, -12],
            ....:                 [0,  48, -16, -16, -188,  20,  92, -16],
            ....:                 [0,   9,  -1,   2,  -33,   5,  18,   0],
            ....:                 [0,  15,  -5,   0,  -59,   7,  30,  -4],
            ....:                 [0, -21,   7,   2,   84, -10, -42,   5],
            ....:                 [0, -42,  14,   8,  167, -17, -84,  13],
            ....:                 [0, -50,  17,  10,  199, -23, -98,  14],
            ....:                 [0,  15,  -5,  -2,  -59,   7,  30, -2]])
            sage: C.minimal_polynomial().factor()                                       # needs sage.libs.pari
            (x - 2)^2
            sage: C.characteristic_polynomial().factor()                                # needs sage.libs.pari
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
            ....:                 [ -6,   5,    7,   2,  -4,   5,    7,   -3],
            ....:                 [ 21, -12,   89,  25,   8,  27,   98,  -95],
            ....:                 [ -9,   5,  -44, -11,  -3, -13,  -48,   47],
            ....:                 [ 23, -13,   74,  21,  12,  22,   85,  -84],
            ....:                 [ 31, -18,  135,  38,  12,  47,  155, -147],
            ....:                 [-33,  19, -138, -39, -13, -45, -156,  151],
            ....:                 [ -7,   4,  -29,  -8,  -3, -10,  -34,  34]])
            sage: D.minimal_polynomial().factor()                                       # needs sage.libs.pari
            (x - 2)^2
            sage: D.characteristic_polynomial().factor()                                # needs sage.libs.pari
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
            ....:                 [-2, -4,   2, -4, -2,   4, -2,   6],
            ....:                 [ 5, 14,  -7, 12,  3,  -8,  6, -27],
            ....:                 [-3, -8,   7, -5,  0,   2, -6,  17],
            ....:                 [ 0,  5,   0,  2,  4,  -4,  1,   2],
            ....:                 [-3, -7,   5, -6, -1,   5, -4,  14],
            ....:                 [ 6, 18, -10, 14,  4, -10, 10, -28],
            ....:                 [-2, -6,   4, -5, -1,   3,  -3, 13]])
            sage: E.minimal_polynomial().factor()                                       # needs sage.libs.pari
            (x - 2)^3
            sage: E.characteristic_polynomial().factor()                                # needs sage.libs.pari
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
            ....: [[-154,  -3,  -54,   44,   48, -244,  -19,   67, -326,   85,   355,   581],
            ....:  [ 504,  25,  156, -145, -171,  793,   99, -213, 1036, -247, -1152, -1865],
            ....:  [ 294,  -1,  112,  -89,  -90,  469,   36, -128,  634, -160,  -695, -1126],
            ....:  [ -49,  -32,  25,    7,   37,  -64,  -58,   12,  -42,  -14,    72,   106],
            ....:  [-261, -123,  65,   47,  169, -358, -254,   70, -309,  -29,   454,   673],
            ....:  [-448, -123, -10,  109,  227, -668, -262,  163, -721,   95,   896,  1410],
            ....:  [  38,    7,   8,  -14,  -17,   66,    6,  -23,   73,  -29,   -78,  -143],
            ....:  [ -96,   10, -55,   37,   24, -168,   17,   56, -231,   88,   237,   412],
            ....:  [ 310,   67,  31,  -81, -143,  473,  143, -122,  538,  -98,  -641, -1029],
            ....:  [ 139,  -35,  99,  -49,  -18,  236,  -41,  -70,  370, -118,  -377,  -619],
            ....:  [ 243,    9,  81,  -72,  -81,  386,   43, -105,  508, -124,  -564,  -911],
            ....:  [-155,   -3, -55,   45,   50, -245,  -27,   65, -328,   77,   365,  583]])
            sage: A.characteristic_polynomial().factor()                                # needs sage.libs.pari
            (x^2 - 2)^2 * (x^2 + 2*x + 5)^4
            sage: A.eigenvalues(extend=False)                                           # needs sage.libs.pari
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
            sage: [F(p).factor() for p in polys]                                        # needs sage.libs.pari
            [x^2 + 2*x + 5, (x^2 - 2) * (x^2 + 2*x + 5), (x^2 - 2) * (x^2 + 2*x + 5)^2]

        Rational form may be computed over any field.  The matrix below is
        an example where the eigenvalues lie outside the field.  ::

            sage: # needs sage.rings.finite_rings
            sage: F.<a> = FiniteField(7^2)
            sage: A = matrix(F,
            ....: [[5*a + 3, 4*a + 1, 6*a + 2, 2*a + 5,       6, 4*a + 5, 4*a + 5,       5,   a + 6,      5,  4*a + 4],
            ....:  [6*a + 3, 2*a + 4,       0,       6, 5*a + 5,     2*a, 5*a + 1,       1, 5*a + 2,     4*a, 5*a + 6],
            ....:  [3*a + 1, 6*a + 6,   a + 6,       2,       0, 3*a + 6, 5*a + 4, 5*a + 6, 5*a + 2,       3, 4*a + 2],
            ....:  [    3*a,     6*a,     3*a,     4*a, 4*a + 4, 3*a + 6,     6*a,       4, 3*a + 4, 6*a + 2,     4*a],
            ....:  [4*a + 5,   a + 1, 4*a + 3, 6*a + 5, 5*a + 2, 5*a + 2,     6*a, 4*a + 6, 6*a + 4, 5*a + 3, 3*a + 1],
            ....:  [    3*a,     6*a, 4*a + 1, 6*a + 2, 2*a + 5, 4*a + 6,       2,   a + 5, 2*a + 4, 2*a + 1, 2*a + 1],
            ....:  [4*a + 5, 3*a + 3,       6, 4*a + 1, 4*a + 3, 6*a + 3,       6, 3*a + 3,       3,   a + 3,       0],
            ....:  [6*a + 6,   a + 4, 2*a + 6, 3*a + 5, 4*a + 3,       2,       a, 3*a + 4,     5*a, 2*a + 5, 4*a + 3],
            ....:  [3*a + 5, 6*a + 2,     4*a,   a + 5,       0,     5*a, 6*a + 5, 2*a + 1, 3*a + 1, 3*a + 5, 4*a + 2],
            ....:  [3*a + 2,   a + 3, 3*a + 6,       a, 3*a + 5, 5*a + 1, 3*a + 2,   a + 3,   a + 2, 6*a + 1, 3*a + 3],
            ....:  [6*a + 6, 5*a + 1,     4*a,       2, 5*a + 5, 3*a + 5, 3*a + 1,     2*a,     2*a, 2*a + 4, 4*a + 2]])
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
            TypeError: no canonical coercion from Finite Field in a of size 7^2 to Finite Field in z2 of size 7^2

        Companion matrices may be selected as any one of four different types.
        See the documentation for the companion matrix constructor,
        :meth:`sage.matrix.constructor.companion_matrix`, for more information. ::

            sage: A = matrix(QQ, [[35, -18, -2, -45],
            ....:                 [22, -22, 12, -16],
            ....:                 [ 5, -12, 12,   4],
            ....:                 [16,  -6, -4, -23]])
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
        """
        from sage.arith.misc import GCD as gcd
        import sage.rings.polynomial.polynomial_ring_constructor
        from sage.matrix.constructor import (block_diagonal_matrix,
                                             companion_matrix)

        R = self.base_ring()
        if not self.is_square():
            raise TypeError("matrix must be square, not {0} x {1}".format(self.nrows(), self.ncols()))
        if not (R in _Fields and R.is_exact()):
            raise TypeError("matrix entries must come from an exact field, not {0}".format(R))
        if format not in ['right', 'bottom', 'left', 'top', 'invariants']:
            raise ValueError("'format' keyword must be 'right', 'bottom', 'left', 'top' or 'invariants', not {0}".format(format))
        if subdivide not in [True, False]:
            raise ValueError("'subdivide' keyword must be True or False, not {0}".format(subdivide))

        _, polys, corners = self._zigzag_form(basis=False)
        cdef Py_ssize_t k = len(polys), j, i
        F = sage.rings.polynomial.polynomial_ring_constructor.PolynomialRing(R, 'x')
        cdef list C = [F(p) for p in polys]
        cdef list B = [b.is_one() for b in corners]
        B.append(False)  # no last block, so no corner

        if B[0]:
            V = [F.one()]
        else:
            V = [F.zero()]

        for j in range(1, k):
            V[j-1] = gcd([V[j-1], C[j], C[j-1]])
            for i in range(j-2, -1, -1):
                V[i] = gcd([V[i], V[i+1], C[i]])
            m = -F.one()
            for i in range(j):
                g = gcd(m*V[i], C[i])
                q, _ = C[i].quo_rem(g)
                C[i] = g
                if B[j]:
                    _, V[i] = m.quo_rem(C[i])
                else:
                    V[i] = F.zero()
                m = m * q
            C[j] = m * C[j]
            if B[j]:
                V.append(m)
            else:
                V.append(F.zero())

        # Leading constant polynomials in C are size zero blocks, so toss them
        # Massage remainder to have leading coefficient 1
        while C and not C[0].degree():
            del C[0]
        for i in range(len(C)):
            unit = C[i].leading_coefficient()
            if not unit.is_one():
                C[i] = ~unit*C[i]

        if format == 'invariants':
            inv = []
            for poly in C:
                inv.append(poly.list())
            return inv
        elif format in ['right', 'left', 'top', 'bottom']:
            companions = []
            for poly in C:
                companions.append(companion_matrix(poly, format=format))
            return block_diagonal_matrix(companions, subdivide=subdivide)

    def is_positive_operator_on(self, K1, K2=None):
        r"""
        Determine if this matrix is a positive operator on a cone.

        A matrix is a positive operator on a cone if the image of the
        cone under the matrix is itself a subset of the cone. That
        concept can be extended to two cones: a matrix is a positive
        operator on a pair of cones if the image of the first cone is
        contained in the second cone.

        To reliably check whether or not this matrix is a positive
        operator, its base ring must be either exact (for example, the
        rationals) or the symbolic ring. An exact ring is more reliable,
        but in some cases a matrix whose entries contain symbolic
        constants like `e` and `\pi` will work. Performance is best
        for integer or rational matrices, for which we can check the
        "is a subset of the other cone" condition quickly.

        INPUT:

        - ``K1`` -- a polyhedral closed convex cone

        - ``K2`` -- (default: ``K1``) the codomain cone; this matrix is
          a positive operator if the image of ``K1`` is a subset of ``K2``

        OUTPUT:

        If the base ring of this matrix is exact, then ``True`` will be
        returned if and only if this matrix is a positive operator.

        If the base ring of this matrix is symbolic, then the situation is
        more complicated:

        - ``True`` will be returned if it can be proven that this matrix
          is a positive operator.
        - ``False`` will be returned if it can be proven that this matrix
          is not a positive operator.
        - ``False`` will also be returned if we can't decide; specifically
          if we arrive at a symbolic inequality that cannot be resolved.

        .. SEEALSO::

              :meth:`is_cross_positive_on`,
              :meth:`is_Z_operator_on`,
              :meth:`is_lyapunov_like_on`

        REFERENCES:

        A. Berman and P. Gaiha. A generalization of irreducible
        monotonicity. Linear Algebra and its Applications, 5:29-38,
        1972.

        A. Berman and R. J. Plemmons. Nonnegative Matrices in the
        Mathematical Sciences. SIAM, Philadelphia, 1994.

        EXAMPLES:

        Nonnegative matrices are positive operators on the nonnegative
        orthant::

            sage: K = Cone([(1,0,0), (0,1,0), (0,0,1)])                                 # needs sage.geometry.polyhedron
            sage: L = random_matrix(QQ, 3).apply_map(abs)
            sage: L.is_positive_operator_on(K)                                          # needs sage.geometry.polyhedron
            True

        Symbolic entries also work in some easy cases::

            sage: K = Cone([(1,0,0), (0,1,0), (0,0,1)])                                 # needs sage.geometry.polyhedron
            sage: L = matrix(SR, [ [0,       e, 0 ],                                    # needs sage.symbolic
            ....:                  [0,       2, pi],
            ....:                  [sqrt(2), 0, 0 ] ])
            sage: L.is_positive_operator_on(K)                                          # needs sage.geometry.polyhedron
            True

        Your matrix can be over any exact ring, for example the ring of
        univariate polynomials with rational coefficients::

            sage: # needs sage.geometry.polyhedron
            sage: K = Cone([(1,0), (-1,0), (0,1), (0,-1)])
            sage: K.is_full_space()
            True
            sage: x = polygen(ZZ, 'x')
            sage: L = matrix(QQ[x], [[x,0],[0,1]])
            sage: L.is_positive_operator_on(K)
            True

        TESTS:

        The identity matrix is always a positive operator::

            sage: # needs sage.geometry.polyhedron
            sage: K = random_cone(max_ambient_dim=8)
            sage: R = K.lattice().vector_space().base_ring()
            sage: L = identity_matrix(R, K.lattice_dim())
            sage: L.is_positive_operator_on(K)
            True

        The zero matrix is always a positive operator::

            sage: # needs sage.geometry.polyhedron
            sage: K = random_cone(max_ambient_dim=8)
            sage: R = K.lattice().vector_space().base_ring()
            sage: L = zero_matrix(R, K.lattice_dim())
            sage: L.is_positive_operator_on(K)
            True

        Everything in ``K1.positive_operators_gens(K2)`` should be
        positive on ``K1`` with respect to ``K2``, even if we make
        the underlying ring symbolic (the usual case is tested by
        the ``positive_operators_gens`` method)::

            sage: # needs sage.geometry.polyhedron sage.symbolic
            sage: K1 = random_cone(max_ambient_dim=5)
            sage: K2 = random_cone(max_ambient_dim=5)
            sage: results = ( L.change_ring(SR).is_positive_operator_on(K1, K2)
            ....:             for L in K1.positive_operators_gens(K2) )
            sage: all(results)                  # long time
            True

        Technically we could test this, but for now only closed convex cones
        are supported as our ``K1`` and ``K2`` arguments::

            sage: K = [ vector([1,2,3]), vector([5,-1,7]) ]
            sage: L = identity_matrix(3)
            sage: L.is_positive_operator_on(K)
            Traceback (most recent call last):
            ...
            TypeError: K1 and K2 must be cones.

        We can't give reliable answers over inexact rings::

            sage: K = Cone([(1,2,3), (4,5,6)])                                          # needs sage.geometry.polyhedron
            sage: L = identity_matrix(RR, 3)
            sage: L.is_positive_operator_on(K)                                          # needs sage.geometry.polyhedron
            Traceback (most recent call last):
            ...
            ValueError: The base ring of the matrix is neither symbolic nor
            exact.

        Symbolic subrings are fine::

            sage: # needs sage.geometry.polyhedron sage.symbolic
            sage: SCR = SR.subring(no_variables=True); SCR
            Symbolic Constants Subring
            sage: K = Cone([(1,2,3), (4,5,6)])
            sage: L = identity_matrix(SCR, 3)
            sage: L.is_positive_operator_on(K)
            True
        """
        import sage.geometry.abc

        if K2 is None:
            K2 = K1
        if not (isinstance(K1, sage.geometry.abc.ConvexRationalPolyhedralCone)
                and isinstance(K2, sage.geometry.abc.ConvexRationalPolyhedralCone)):
            raise TypeError('K1 and K2 must be cones.')
        if not self.base_ring().is_exact() and not isinstance(self.base_ring(), sage.rings.abc.SymbolicRing):
            msg = 'The base ring of the matrix is neither symbolic nor exact.'
            raise ValueError(msg)

        if self.base_ring() is ZZ or self.base_ring() is QQ:
            # This should be way faster than computing the dual and
            # checking a bunch of inequalities, but it doesn't work if
            # ``self*x`` is symbolic, polynomial, or otherwise
            # contains something unexpected. For example, ``e in
            # Cone([(1,)])`` is true, but returns ``False``.
            return all(self * x in K2 for x in K1)
        else:
            # Fall back to inequality-checking when the entries of
            # this matrix might be symbolic, polynomial, or something
            # else weird.
            return all(s * (self * x) >= 0 for x in K1 for s in K2.dual())

    def is_cross_positive_on(self, K):
        r"""
        Determine if this matrix is cross-positive on a cone.

        We say that a matrix `L` is cross-positive on a closed convex
        cone `K` if the inner product of `Lx` and `s` is nonnegative for
        all pairs of orthogonal vectors `x` in `K` and `s` in the dual
        of `K`. This property need only be checked for generators of `K`
        and its dual.

        To reliably check whether or not this matrix is cross-positive,
        its base ring must be either exact (for example, the rationals)
        or the symbolic ring. An exact ring is more reliable, but in
        some cases a matrix whose entries contain symbolic constants
        like `e` and `\pi` will work.

        INPUT:

        - ``K`` -- a polyhedral closed convex cone

        OUTPUT:

        If the base ring of this matrix is exact, then ``True`` will be
        returned if and only if this matrix is cross-positive on ``K``.

        If the base ring of this matrix is symbolic, then the situation
        is more complicated:

        - ``True`` will be returned if it can be proven that this matrix
          is cross-positive on ``K``.
        - ``False`` will be returned if it can be proven that this matrix
          is not cross-positive on ``K``.
        - ``False`` will also be returned if we can't decide; specifically
          if we arrive at a symbolic inequality that cannot be resolved.

        .. SEEALSO::

              :meth:`is_positive_operator_on`,
              :meth:`is_Z_operator_on`,
              :meth:`is_lyapunov_like_on`

        REFERENCES:

        H. Schneider and M. Vidyasagar. Cross-positive matrices. SIAM
        Journal on Numerical Analysis, 7:508-519, 1970.

        EXAMPLES:

        Negative Z-matrices are cross-positive operators on the
        nonnegative orthant::

            sage: K = Cone([(1,0,0), (0,1,0), (0,0,1)])                                 # needs sage.geometry.polyhedron
            sage: L = matrix(SR, [ [-1, 2, 0],                                          # needs sage.symbolic
            ....:                  [ 0, 2, 7],
            ....:                  [ 3, 0, 3] ])
            sage: L.is_cross_positive_on(K)                                             # needs sage.geometry.polyhedron sage.symbolic
            True

        Symbolic entries also work in some easy cases::

            sage: K = Cone([(1,0,0),(0,1,0),(0,0,1)])                                   # needs sage.geometry.polyhedron
            sage: L = matrix(SR, [ [-1,       e, 0 ],                                   # needs sage.symbolic
            ....:                  [ 0,       2, pi],
            ....:                  [ sqrt(2), 0, 3 ] ])
            sage: L.is_cross_positive_on(K)                                             # needs sage.geometry.polyhedron sage.symbolic
            True

        TESTS:

        The identity matrix is always cross-positive::

            sage: # needs sage.geometry.polyhedron
            sage: K = random_cone(max_ambient_dim=8)
            sage: R = K.lattice().vector_space().base_ring()
            sage: L = identity_matrix(R, K.lattice_dim())
            sage: L.is_cross_positive_on(K)
            True

        The zero matrix is always cross-positive::

            sage: # needs sage.geometry.polyhedron
            sage: K = random_cone(max_ambient_dim=8)
            sage: R = K.lattice().vector_space().base_ring()
            sage: L = zero_matrix(R, K.lattice_dim())
            sage: L.is_cross_positive_on(K)
            True

        Everything in ``K.cross_positive_operators_gens()`` should be
        cross-positive on ``K``, even if we make the underlying ring
        symbolic (the usual case is tested by the
        ``cross_positive_operators_gens`` method)::

            sage: # needs sage.geometry.polyhedron sage.symbolic
            sage: K = random_cone(max_ambient_dim=5)
            sage: results = ( L.change_ring(SR).is_cross_positive_on(K)
            ....:             for L in K.cross_positive_operators_gens() )
            sage: all(results)                  # long time
            True

        Technically we could test this, but for now only closed convex cones
        are supported as our ``K`` argument::

            sage: L = identity_matrix(3)
            sage: K = [ vector([8,2,-8]), vector([5,-5,7]) ]
            sage: L.is_cross_positive_on(K)
            Traceback (most recent call last):
            ...
            TypeError: K must be a cone.

        We can't give reliable answers over inexact rings::

            sage: K = Cone([(1,2,3), (4,5,6)])                                          # needs sage.geometry.polyhedron
            sage: L = identity_matrix(RR, 3)
            sage: L.is_cross_positive_on(K)                                             # needs sage.geometry.polyhedron
            Traceback (most recent call last):
            ...
            ValueError: The base ring of the matrix is neither symbolic nor
            exact.

        Symbolic subrings are fine::

            sage: # needs sage.geometry.polyhedron sage.symbolic
            sage: SCR = SR.subring(no_variables=True); SCR
            Symbolic Constants Subring
            sage: K = Cone([(1,2,3), (4,5,6)])
            sage: L = identity_matrix(SCR, 3)
            sage: L.is_cross_positive_on(K)
            True
        """
        import sage.geometry.abc

        if not isinstance(K, sage.geometry.abc.ConvexRationalPolyhedralCone):
            raise TypeError('K must be a cone.')
        if not self.base_ring().is_exact() and not isinstance(self.base_ring(), sage.rings.abc.SymbolicRing):
            msg = 'The base ring of the matrix is neither symbolic nor exact.'
            raise ValueError(msg)

        return all(s * (self * x) >= 0
                   for (x, s) in K.discrete_complementarity_set())

    def is_Z_operator_on(self, K):
        r"""
        Determine if this matrix is a Z-operator on a cone.

        We say that a matrix `L` is a Z-operator on a closed convex cone
        `K` if the inner product of `Lx` and `s` is nonpositive for all
        pairs of orthogonal vectors `x` in `K` and `s` in the dual of
        `K`. This property need only be checked for generators of `K`
        and its dual.

        A matrix is a Z-operator on `K` if and only if its negation is a
        cross-positive operator on `K`.

        To reliably check whether or not this matrix is a Z operator,
        its base ring must be either exact (for example, the rationals)
        or the symbolic ring. An exact ring is more reliable, but in
        some cases a matrix whose entries contain symbolic constants
        like `e` and `\pi` will work.

        INPUT:

        - ``K`` -- a polyhedral closed convex cone

        OUTPUT:

        If the base ring of this matrix is exact, then ``True`` will be
        returned if and only if this matrix is a Z-operator on ``K``.

        If the base ring of this matrix is symbolic, then the situation
        is more complicated:

        - ``True`` will be returned if it can be proven that this matrix
          is a Z-operator on ``K``.
        - ``False`` will be returned if it can be proven that this matrix
          is not a Z-operator on ``K``.
        - ``False`` will also be returned if we can't decide; specifically
          if we arrive at a symbolic inequality that cannot be resolved.

        .. SEEALSO::

              :meth:`is_positive_operator_on`,
              :meth:`is_cross_positive_on`,
              :meth:`is_lyapunov_like_on`

        REFERENCES:

        A. Berman and R. J. Plemmons. Nonnegative Matrices in the
        Mathematical Sciences. SIAM, Philadelphia, 1994.

        EXAMPLES:

        Z-matrices are Z-operators on the nonnegative orthant::

            sage: K = Cone([(1,0,0), (0,1,0), (0,0,1)])                                 # needs sage.geometry.polyhedron
            sage: L = matrix(SR, [ [-1, -2,  0],                                        # needs sage.symbolic
            ....:                  [ 0,  2, -7],
            ....:                  [-3,  0,  3] ])
            sage: L.is_Z_operator_on(K)                                                 # needs sage.geometry.polyhedron sage.symbolic
            True

        Symbolic entries also work in some easy cases::

            sage: K = Cone([(1,0,0), (0,1,0), (0,0,1)])                                 # needs sage.geometry.polyhedron
            sage: L = matrix(SR, [ [-1,      -e,  0 ],                                  # needs sage.symbolic
            ....:                  [ 0,       2, -pi],
            ....:                  [-sqrt(2), 0,  3 ] ])
            sage: L.is_Z_operator_on(K)                                                 # needs sage.geometry.polyhedron sage.symbolic
            True

        TESTS:

        The identity matrix is always a Z-operator::

            sage: # needs sage.geometry.polyhedron
            sage: K = random_cone(max_ambient_dim=8)
            sage: R = K.lattice().vector_space().base_ring()
            sage: L = identity_matrix(R, K.lattice_dim())
            sage: L.is_Z_operator_on(K)
            True

        The zero matrix is always a Z-operator::

            sage: # needs sage.geometry.polyhedron
            sage: K = random_cone(max_ambient_dim=8)
            sage: R = K.lattice().vector_space().base_ring()
            sage: L = zero_matrix(R, K.lattice_dim())
            sage: L.is_Z_operator_on(K)
            True

        Everything in ``K.Z_operators_gens()`` should be a Z-operator on
        ``K``, , even if we make the underlying ring symbolic (the usual
        case is tested by the ``Z_operators_gens`` method)::

            sage: K = random_cone(max_ambient_dim=5)                                    # needs sage.geometry.polyhedron
            sage: all(L.change_ring(SR).is_Z_operator_on(K)     # long time             # needs sage.geometry.polyhedron sage.symbolic
            ....:     for L in K.Z_operators_gens())
            True

        Technically we could test this, but for now only closed convex cones
        are supported as our ``K`` argument::

            sage: L = identity_matrix(3)
            sage: K = [ vector([-4,20,3]), vector([1,-5,2]) ]
            sage: L.is_Z_operator_on(K)
            Traceback (most recent call last):
            ...
            TypeError: K must be a cone.

        We can't give reliable answers over inexact rings::

            sage: K = Cone([(1,2,3), (4,5,6)])                                          # needs sage.geometry.polyhedron
            sage: L = identity_matrix(RR, 3)
            sage: L.is_Z_operator_on(K)                                                 # needs sage.geometry.polyhedron
            Traceback (most recent call last):
            ...
            ValueError: The base ring of the matrix is neither symbolic nor
            exact.
        """
        return (-self).is_cross_positive_on(K)

    def is_lyapunov_like_on(self, K):
        r"""
        Determine if this matrix is Lyapunov-like on a cone.

        We say that a matrix `L` is Lyapunov-like on a closed convex
        cone `K` if the inner product of `Lx` and `s` is zero for all
        pairs of orthogonal vectors `x` in `K` and `s` in the dual of
        `K`. This property need only be checked for generators of `K`
        and its dual.

        An operator is Lyapunov-like on `K` if and only if both the
        operator itself and its negation are cross-positive on `K`.

        To reliably check whether or not this matrix is Lyapunov-like,
        its base ring must be either exact (for example, the rationals)
        or the symbolic ring. An exact ring is more reliable, but in
        some cases a matrix whose entries contain symbolic constants
        like `e` and `\pi` will work.

        INPUT:

        - ``K`` -- a polyhedral closed convex cone

        OUTPUT:

        If the base ring of this matrix is exact, then ``True`` will be
        returned if and only if this matrix is Lyapunov-like on ``K``.

        If the base ring of this matrix is symbolic, then the situation
        is more complicated:

        - ``True`` will be returned if it can be proven that this matrix
          is Lyapunov-like on ``K``.
        - ``False`` will be returned if it can be proven that this matrix
          is not Lyapunov-like on ``K``.
        - ``False`` will also be returned if we can't decide; specifically
          if we arrive at a symbolic inequality that cannot be resolved.

        .. SEEALSO::

              :meth:`is_positive_operator_on`,
              :meth:`is_cross_positive_on`,
              :meth:`is_Z_operator_on`

        REFERENCES:

        - [Or2017]_

        EXAMPLES:

        Diagonal matrices are Lyapunov-like operators on the nonnegative
        orthant::

            sage: K = Cone([(1,0,0), (0,1,0), (0,0,1)])                                 # needs sage.geometry.polyhedron
            sage: L = diagonal_matrix(random_vector(QQ, 3))
            sage: L.is_lyapunov_like_on(K)                                              # needs sage.geometry.polyhedron
            True

        Symbolic entries also work in some easy cases::

            sage: K = Cone([(1,0,0), (0,1,0), (0,0,1)])                                 # needs sage.geometry.polyhedron
            sage: L = matrix(SR, [ [e, 0,  0      ],                                    # needs sage.symbolic
            ....:                  [0, pi, 0      ],
            ....:                  [0, 0,  sqrt(2)] ])
            sage: L.is_lyapunov_like_on(K)                                              # needs sage.geometry.polyhedron
            True

        TESTS:

        The identity matrix is always Lyapunov-like::

            sage: # needs sage.geometry.polyhedron
            sage: K = random_cone(max_ambient_dim=8)
            sage: R = K.lattice().vector_space().base_ring()
            sage: L = identity_matrix(R, K.lattice_dim())
            sage: L.is_lyapunov_like_on(K)
            True

        The zero matrix is always Lyapunov-like::

            sage: # needs sage.geometry.polyhedron
            sage: K = random_cone(max_ambient_dim=8)
            sage: R = K.lattice().vector_space().base_ring()
            sage: L = zero_matrix(R, K.lattice_dim())
            sage: L.is_lyapunov_like_on(K)
            True

        Everything in ``K.lyapunov_like_basis()`` should be
        Lyapunov-like on ``K``, even if we make the underlying ring
        symbolic (the usual case is tested by the
        ``lyapunov_like_basis`` method)::

            sage: K = random_cone(max_ambient_dim=5)                                    # needs sage.geometry.polyhedron
            sage: all(L.change_ring(SR).is_lyapunov_like_on(K)  # long time             # needs sage.geometry.polyhedron sage.symbolic
            ....:     for L in K.lyapunov_like_basis())
            True

        Technically we could test this, but for now only closed convex cones
        are supported as our ``K`` argument::

            sage: L = identity_matrix(3)
            sage: K = [ vector([2,2,-1]), vector([5,4,-3]) ]
            sage: L.is_lyapunov_like_on(K)
            Traceback (most recent call last):
            ...
            TypeError: K must be a cone.

        We can't give reliable answers over inexact rings::

            sage: K = Cone([(1,2,3), (4,5,6)])                                          # needs sage.geometry.polyhedron
            sage: L = identity_matrix(RR, 3)
            sage: L.is_lyapunov_like_on(K)                                              # needs sage.geometry.polyhedron
            Traceback (most recent call last):
            ...
            ValueError: The base ring of the matrix is neither symbolic nor
            exact.

        Symbolic subrings are fine::

            sage: # needs sage.geometry.polyhedron sage.symbolic
            sage: SCR = SR.subring(no_variables=True); SCR
            Symbolic Constants Subring
            sage: K = Cone([(1,2,3), (4,5,6)])
            sage: L = identity_matrix(SCR, 3)
            sage: L.is_lyapunov_like_on(K)
            True

        A matrix is Lyapunov-like on a cone if and only if both the
        matrix and its negation are cross-positive on the cone::

            sage: # needs sage.geometry.polyhedron
            sage: K = random_cone(max_ambient_dim=5)
            sage: R = K.lattice().vector_space().base_ring()
            sage: L = random_matrix(R, K.lattice_dim())
            sage: actual = L.is_lyapunov_like_on(K)             # long time
            sage: expected = (L.is_cross_positive_on(K) and     # long time
            ....:             (-L).is_cross_positive_on(K))
            sage: actual == expected                            # long time
            True
        """
        import sage.geometry.abc

        if not isinstance(K, sage.geometry.abc.ConvexRationalPolyhedralCone):
            raise TypeError('K must be a cone.')
        if not self.base_ring().is_exact() and not isinstance(self.base_ring(), sage.rings.abc.SymbolicRing):
            msg = 'The base ring of the matrix is neither symbolic nor exact.'
            raise ValueError(msg)

        # Even though ``discrete_complementarity_set`` is a cached
        # method of cones, this is faster than calling
        # :meth:`is_cross_positive_on` twice: doing so checks twice as
        # many inequalities as the number of equalities that we're
        # about to check.
        return all(s * (self * x) == 0
                   for (x, s) in K.discrete_complementarity_set())

    def LLL_gram(self, flag=0):
        """
        Return the LLL transformation matrix for this Gram matrix.

        That is, the transformation matrix U over ZZ of determinant 1
        that transforms the lattice with this matrix as Gram matrix
        to a lattice that is LLL-reduced.

        Always works when ``self`` is positive definite,
        might work in some semidefinite and indefinite cases.

        INPUT:

        - ``self`` -- the Gram matrix of a quadratic form or of
          a lattice equipped with a bilinear form

        - ``flag`` -- an optional flag passed to ``qflllgram``
          According  to :pari:`qflllgram`'s documentation the options are:

            - ``0`` -- (default), assume that ``self`` has either exact
              (integral or rational) or real floating point entries.
              The matrix is rescaled, converted to integers and the
              behavior is then as in ``flag=1``.

            - ``1`` -- assume that G is integral.
              Computations involving Gram-Schmidt vectors are
              approximate, with precision varying as needed.

        OUTPUT:

        A dense matrix ``U`` over the integers with determinant 1
        such that ``U.T * M * U`` is LLL-reduced.

        ALGORITHM:

        Calls PARI's :pari:`qflllgram`.

        EXAMPLES:

        Create a Gram matrix and LLL-reduce it::

            sage: M = Matrix(ZZ, 2, 2, [5, 3, 3, 2])
            sage: U = M.LLL_gram()                                                      # needs sage.libs.pari
            sage: MM = U.transpose() * M * U                                            # needs sage.libs.pari
            sage: M, U, MM                                                              # needs sage.libs.pari
            (
            [5 3]  [-1  1]  [1 0]
            [3 2], [ 1 -2], [0 1]
            )

        For a Gram matrix over RR with a length one first vector and
        a very short second vector, the LLL-reduced basis is obtained
        by swapping the two basis vectors (and changing sign to
        preserve orientation). ::

            sage: M = Matrix(RDF, 2, 2, [1, 0, 0, 1e-5])
            sage: M.LLL_gram()                                                          # needs sage.libs.pari
            [ 0 -1]
            [ 1  0]

        The algorithm might work for some semidefinite and indefinite forms::

            sage: Matrix(ZZ, 2, 2, [2, 6, 6, 3]).LLL_gram()                             # needs sage.libs.pari
            [-3 -1]
            [ 1  0]
            sage: Matrix(ZZ, 2, 2, [1, 0, 0, -1]).LLL_gram()                            # needs sage.libs.pari
            [ 0 -1]
            [ 1  0]

        However, it might fail for others, either raising a :exc:`ValueError`::

            sage: Matrix(ZZ, 1, 1, [0]).LLL_gram()                                      # needs sage.libs.pari
            Traceback (most recent call last):
            ...
            ValueError: qflllgram did not return a square matrix,
            perhaps the matrix is not positive definite

            sage: Matrix(ZZ, 2, 2, [0, 1, 1, 0]).LLL_gram()                             # needs sage.libs.pari
            Traceback (most recent call last):
            ...
            ValueError: qflllgram did not return a square matrix,
            perhaps the matrix is not positive definite

        or running forever::

            sage: Matrix(ZZ, 2, 2, [-5, -1, -1, -5]).LLL_gram()         # not tested, needs sage.libs.pari
            Traceback (most recent call last):
            ...
            RuntimeError: infinite loop while calling qflllgram

        Nonreal input leads to a value error::

           sage: Matrix(2, 2, [CDF(1, 1), 0, 0, 1]).LLL_gram()                          # needs sage.libs.pari
           Traceback (most recent call last):
           ...
           ValueError: qflllgram failed, perhaps the matrix is not positive definite
        """
        if self._nrows != self._ncols:
            raise ArithmeticError("self must be a square matrix")
        n = self.nrows()
        P = self.__pari__()
        try:
            if self.base_ring() == ZZ:
                U = P.lllgramint()
            else:
                U = P.qflllgram(flag)
        except (RuntimeError, ArithmeticError) as msg:
            raise ValueError("qflllgram failed, "
                             "perhaps the matrix is not positive definite")
        if U.matsize() != [n, n]:
            raise ValueError("qflllgram did not return a square matrix, "
                             "perhaps the matrix is not positive definite")
        from sage.matrix.matrix_space import MatrixSpace
        MS = MatrixSpace(ZZ, n)
        U = MS(U.sage())
        # Fix last column so that det = +1
        from sage.rings.finite_rings.finite_field_constructor import FiniteField
        if U.change_ring(FiniteField(3)).det() != 1:  # p = 3 is enough to decide
            U.rescale_col(n - 1, -1)
        return U

    # a limited number of access-only properties are provided for matrices
    @property
    def T(self):
        r"""
        Return the transpose of a matrix.

        EXAMPLES::

            sage: A = matrix(QQ, 5, range(25))
            sage: A.T
            [ 0  5 10 15 20]
            [ 1  6 11 16 21]
            [ 2  7 12 17 22]
            [ 3  8 13 18 23]
            [ 4  9 14 19 24]
        """
        return self.transpose()

    @property
    def C(self):
        r"""
        Return the conjugate matrix.

        EXAMPLES::

            sage: A = matrix(QQbar, [[     -3,  5 - 3*I, 7 - 4*I],                      # needs sage.rings.number_field
            ....:                    [7 + 3*I, -1 + 6*I, 3 + 5*I],
            ....:                    [3 + 3*I, -3 + 6*I, 5 +   I]])
            sage: A.C                                                                   # needs sage.rings.number_field
            [      -3  5 + 3*I  7 + 4*I]
            [ 7 - 3*I -1 - 6*I  3 - 5*I]
            [ 3 - 3*I -3 - 6*I  5 - 1*I]
        """
        return self.conjugate()

    @property
    def H(self):
        r"""
        Return the conjugate-transpose (Hermitian) matrix.

        EXAMPLES::

            sage: A = matrix(QQbar, [[     -3,  5 - 3*I, 7 - 4*I],                      # needs sage.rings.number_field
            ....:                    [7 + 3*I, -1 + 6*I, 3 + 5*I],
            ....:                    [3 + 3*I, -3 + 6*I, 5 +   I]])
            sage: A.H                                                                   # needs sage.rings.number_field
            [      -3  7 - 3*I  3 - 3*I]
            [ 5 + 3*I -1 - 6*I -3 - 6*I]
            [ 7 + 4*I  3 - 5*I  5 - 1*I]
        """
        return self.conjugate().transpose()


def _smith_diag(d, transformation=True):
    r"""
    For internal use by the smith_form routine. Given a diagonal matrix d
    over a ring r, return matrices d', a,b such that a\*d\*b = d' and
    d' is diagonal with each entry dividing the next.

    If any of the d's is a unit, it replaces it with 1 (but no other
    attempt is made to pick "good" representatives of ideals).

    EXAMPLES::

        sage: # needs sage.rings.number_field
        sage: from sage.matrix.matrix2 import _smith_diag
        sage: x = polygen(ZZ, 'x')
        sage: OE = EquationOrder(x^2 - x + 2, 'w')
        sage: A = matrix(OE, 2, [2, 0, 0, 3])
        sage: D,U,V = _smith_diag(A); D,U,V
        (
        [1 0]  [2 1]  [ 1 -3]
        [0 6], [3 2], [-1  4]
        )
        sage: D == U*A*V
        True
        sage: m = matrix(GF(7), 2, [3,0,0,6]); d,u,v = _smith_diag(m); d
        [1 0]
        [0 1]
        sage: u*m*v == d
        True
    """

    dp = d.__copy__()
    n = min(d.nrows(), d.ncols())
    R = d.base_ring()
    if transformation:
        left = d.new_matrix(d.nrows(), d.nrows(), 1)
        right = d.new_matrix(d.ncols(), d.ncols(), 1)
    else:
        left = right = None
    for i in range(n):
        I = ideal_or_fractional(R, dp[i, i])

        if I == ideal_or_fractional(R, 1):
            if dp[i, i] != 1:
                if transformation:
                    left.add_multiple_of_row(i, i, R(R(1)/(dp[i, i])) - 1)
                dp[i, i] = R(1)
            continue

        for j in range(i+1, n):
            if dp[j, j] not in I:
                t = ideal_or_fractional(R, [dp[i, i], dp[j, j]]).gens_reduced()
                if len(t) > 1:
                    raise ArithmeticError
                t = t[0]
                # find lambda, mu such that lambda*d[i,i] + mu*d[j,j] = t
                lamb = R(dp[i, i]/t).inverse_mod(ideal_or_fractional(R, dp[j, j]/t))
                mu = R((t - lamb*dp[i, i]) / dp[j, j])

                newlmat = dp.new_matrix(dp.nrows(), dp.nrows(), 1)
                newlmat[i, i] = lamb
                newlmat[i, j] = 1
                newlmat[j, i] = R(-dp[j, j]*mu/t)
                newlmat[j, j] = R(dp[i, i]/t)
                newrmat = dp.new_matrix(dp.ncols(), dp.ncols(), 1)
                newrmat[i, i] = 1
                newrmat[i, j] = R(-dp[j, j]/t)
                newrmat[j, i] = mu
                newrmat[j, j] = R(lamb*dp[i, i] / t)

                if transformation:
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

        sage: # needs sage.rings.number_field
        sage: x = polygen(ZZ, 'x')
        sage: L.<w> = NumberField(x^2 - x + 2)
        sage: OL = L.ring_of_integers(); w = OL(w)
        sage: m = matrix(OL, 8, 4, [2*w - 2, 2*w + 1, -2, w, 2, -2, -2*w - 2, -2*w + 2, -w + 2, 2*w + 1, -w + 2, -w - 2, -2*w,
        ....:     2*w, -w+ 2, w - 1, -2*w + 2, 2*w + 2, 2*w - 1, -w, 2*w + 2, -w + 2, 2, 2*w -1, w - 4, -2*w - 2, 2*w - 1, 0, 6, 7, 2*w + 1, 14])
        sage: s,t = m.echelon_form(transformation=True); t*m == s  # indirect doctest
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
            if k == a.nrows():  # first column is zero
                return left_mat, a
        # k is now first row such that a[k, 0] is nonzero
        left_mat[0, 0] = 0
        left_mat[k, k] = 0
        left_mat[0, k] = 1
        left_mat[k, 0] = -1
        a = left_mat*a
        if left_mat * m != a:
            raise ArithmeticError("Something went wrong")

    # case 2: if there is an entry at position (k,j) such that a_{0,j}
    # does not divide a_{k,j}, then we know that there are c,d in R such
    # that c.a_{0,j} - d.a_{k,j} = B where B = gcd(a_{0,j}, a_{k,j}) (which
    # is well-defined since R is a PID).
    # Then for appropriate e,f the matrix
    # [c,-d]
    # [e,f]
    # is invertible over R

    I = ideal_or_fractional(R, a[0, 0])  # need to make sure we change this when a[0,0] changes
    for k in range(1, a.nrows()):
        if a[k, 0] not in I:
            new_ideal = ideal_or_fractional(R, a[0, 0], a[k, 0])
            try:
                v = new_ideal.gens_reduced()
            except Exception as msg:
                raise ArithmeticError("%s\nCan't create ideal on %s and %s" % (msg, a[0, 0], a[k, 0]))
            if len(v) > 1:
                raise ArithmeticError("Ideal %s not principal" % new_ideal)
            B = v[0]

            # now we find c,d, using the fact that c * (a_{0,0}/B) - d *
            # (a_{k,0}/B) = 1, so c is the inverse of a_{0,0}/B modulo
            # a_{k,0}/B.
            # need to handle carefully the case when a_{k,0}/B is a unit, i.e. a_{k,0} divides
            # a_{0,0}.

            c = R(a[0, 0] / B).inverse_mod(ideal_or_fractional(R, a[k, 0] / B))
            d = R((c*a[0, 0] - B)/(a[k, 0]))

            # sanity check
            if c*a[0, 0] - d*a[k, 0] != B:
                raise ArithmeticError

            # now we find e,f such that e*d + c*f = 1 in the same way
            if c != 0:
                e = d.inverse_mod(ideal_or_fractional(R, c))
                f = R((1 - d*e)/c)
            else:
                e = R(-a[k, 0]/B)  # here d is a unit and this is just 1/d
                f = R(1)

            if e*d + c*f != 1:
                raise ArithmeticError
            newlmat = left_mat.parent()(1)
            newlmat[0, 0] = c
            newlmat[0, k] = -d
            newlmat[k, 0] = e
            newlmat[k, k] = f
            if newlmat.det() != 1:
                raise ArithmeticError
            a = newlmat*a
            I = ideal_or_fractional(R, a[0, 0])
            left_mat = newlmat*left_mat
            if left_mat * m != a:
                raise ArithmeticError

    # now everything in column 0 is divisible by the pivot
    for i in range(1, a.nrows()):
        s = R(a[i, 0]/a[0, 0])
        a.add_multiple_of_row(i, 0, -s)
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

    EXAMPLES::

        sage: # needs sage.rings.number_field
        sage: from sage.matrix.matrix2 import _smith_onestep
        sage: x = polygen(ZZ, 'x')
        sage: OE.<w> = EquationOrder(x^2 - x + 2)
        sage: m = matrix(OE, 3, 3, [1, 0, 7, 2, w, w+17, 13+8*w, 0, 6])
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
        right_mat[0, 0] = right_mat[j, j] = 0
        right_mat[0, j] = 1
        right_mat[j, 0] = -1
        a = a*right_mat
        if m * right_mat != a:
            raise ArithmeticError

    left_mat, a = _generic_clear_column(a)
    assert left_mat * m * right_mat == a

    # test if everything to the right of the pivot in row 0 is good as well
    isdone = True
    for jj in range(j+1, a.ncols()):
        if a[0, jj] != 0:
            isdone = False

    # if not we recurse -- algorithm must terminate if R is Noetherian.
    if not isdone:
        s, t, u = _smith_onestep(a.transpose())
        left_mat = u.transpose() * left_mat
        a = t.transpose()
        right_mat = right_mat* s.transpose()

    return left_mat, a, right_mat


def decomp_seq(v):
    """
    This function is used internally be the decomposition matrix
    method. It takes a list of tuples and produces a sequence that is
    correctly sorted and prints with carriage returns.

    EXAMPLES::

        sage: from sage.matrix.matrix2 import decomp_seq
        sage: V = [(QQ^3, 2), (QQ^2, 1)]
        sage: decomp_seq(V)
        [(Vector space of dimension 2 over Rational Field, 1),
         (Vector space of dimension 3 over Rational Field, 2)]
    """
    list.sort(v, key=lambda x: x[0].dimension())
    return Sequence(v, universe=tuple, check=False, cr=True)


def _choose(Py_ssize_t n, Py_ssize_t t):
    """
    Return all possible sublists of length `t` from ``range(n)``.

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
    """
    cdef Py_ssize_t j, temp

    x = []               # initialize T1
    c = list(range(t))
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
    """
    cdef Py_ssize_t i

    if 2 * k > n:
        k = n - k
    if k == 0:
        return 1

    result = n
    n, k = n - 1, k - 1
    i = 2
    while k > 0:
        result = (result * n) // i
        i, n, k = i + 1, n - 1, k - 1
    return result


def _jordan_form_vector_in_difference(V, W):
    r"""
    Given two lists of vectors ``V`` and ``W`` over the same base field,
    returns a vector in the difference ``V - W``.  If the difference is
    empty, returns ``None``.

    .. NOTE::

        This is meant to be a private helper method for the
        ``jordan_form`` method in the above class.

    TESTS::

        sage: v = vector(ZZ, [1,0,0,0])
        sage: w = vector(ZZ, [0,1,0,0])
        sage: u = vector(ZZ, [1,1,0,0])
        sage: sage.matrix.matrix2._jordan_form_vector_in_difference([v,w], [u])
        (1, 0, 0, 0)
    """
    from sage.modules.free_module import span

    if not V:
        return None
    if not W:
        return V[0]
    W_space = span(W)
    for v in V:
        if v not in W_space:
            return v
    return None


def _matrix_power_symbolic(A, n):
    r"""
    Return the symbolic `n`-th power `A^n` of the matrix `A`.

    This function implements the computation of `A^n` for symbolic `n`,
    relying on the Jordan normal form of `A`, available for exact rings
    as :meth:`jordan_form`. See [Hig2008]_, §1.2, for further details.

    INPUT:

    - ``A`` -- a square matrix over an exact field

    - ``n`` -- a symbolic exponent

    OUTPUT:

    The matrix `A^n` (with symbolic entries involving the exponent).

    EXAMPLES:

    General power of a two by two matrix::

        sage: # needs sage.symbolic
        sage: n = SR.var('n')
        sage: A = matrix(QQ, [[2, -1], [1,  0]])
        sage: B = A^n; B
        [ n + 1     -n]
        [     n -n + 1]
        sage: all(A^k == B.subs({n: k}) for k in range(8))
        True

    General power of a three by three matrix in Jordan form::

        sage: # needs sage.symbolic
        sage: n = SR.var('n')
        sage: A = matrix(QQ, 3, [[2, 1, 0], [0, 2, 0], [0, 0, 3]]); A
        [2 1 0]
        [0 2 0]
        [0 0 3]
        sage: B = A^n; B
        [        2^n 2^(n - 1)*n           0]
        [          0         2^n           0]
        [          0           0         3^n]
        sage: all(A^k == B.subs({n: k}) for k in range(8))
        True

    General power of a three by three matrix not in Jordan form::

        sage: # needs sage.symbolic
        sage: A = matrix([[4, 1, 2], [0, 2, -4], [0, 1, 6]]); A
        [ 4  1  2]
        [ 0  2 -4]
        [ 0  1  6]
        sage: B = A^n; B
        [                 4^n          4^(n - 1)*n        2*4^(n - 1)*n]
        [                   0 -2*4^(n - 1)*n + 4^n       -4*4^(n - 1)*n]
        [                   0          4^(n - 1)*n  2*4^(n - 1)*n + 4^n]
        sage: [B.subs({n: k}) for k in range(4)]
        [
        [1 0 0]  [ 4  1  2]  [ 16   8  16]  [  64   48   96]
        [0 1 0]  [ 0  2 -4]  [  0   0 -32]  [   0  -32 -192]
        [0 0 1], [ 0  1  6], [  0   8  32], [   0   48  160]
        ]
        sage: all(A^k == B.subs({n: k}) for k in range(8))
        True

    TESTS:

    Testing exponentiation in the symbolic ring::

        sage: n = var('n')                                                              # needs sage.symbolic
        sage: A = matrix([[pi, e],[0, -2*I]])                                           # needs sage.symbolic
        sage: (A^n).list()                                                              # needs sage.symbolic
        [pi^n,
         -(-2*I)^n/(pi*e^(-1) + 2*I*e^(-1)) + pi^n/(pi*e^(-1) + 2*I*e^(-1)),
         0,
         (-2*I)^n]

    If the base ring is inexact, the Jordan normal form is not available::

        sage: A = matrix(RDF, [[2, -1], [1,  0]])
        sage: A^n                                                                       # needs sage.symbolic
        Traceback (most recent call last):
        ...
        ValueError: Jordan normal form not implemented over inexact rings.

    Testing exponentiation in the integer ring::

        sage: A = matrix(ZZ, [[1,-1], [-1,1]])
        sage: A^(2*n+1)                                                                 # needs sage.symbolic
        [ 1/2*2^(2*n + 1) + 1/2*kronecker_delta(0, 2*n + 1) -1/2*2^(2*n + 1) + 1/2*kronecker_delta(0, 2*n + 1)]
        [-1/2*2^(2*n + 1) + 1/2*kronecker_delta(0, 2*n + 1)  1/2*2^(2*n + 1) + 1/2*kronecker_delta(0, 2*n + 1)]

    Check if :issue:`23215` is fixed::

        sage: a, b, k = var('a, b, k')                                                  # needs sage.symbolic
        sage: (matrix(2, [a, b, -b, a])^k).list()                                       # needs sage.symbolic
        [1/2*(a + I*b)^k + 1/2*(a - I*b)^k,
         -1/2*I*(a + I*b)^k + 1/2*I*(a - I*b)^k,
         1/2*I*(a + I*b)^k - 1/2*I*(a - I*b)^k,
         1/2*(a + I*b)^k + 1/2*(a - I*b)^k]

    Check if :issue:`36838` is fixed:Checking symbolic power of
    nilpotent matrix::

        sage: A = matrix([[0,1],[0,0]]); A
        [0 1]
        [0 0]
        sage: n = var('n'); n
        n
        sage: B = A^n; B
        [  kronecker_delta(0, n) n*kronecker_delta(1, n)]
        [                      0   kronecker_delta(0, n)]
    """
    from sage.rings.qqbar import AlgebraicNumber
    from sage.matrix.constructor import matrix
    from sage.functions.other import binomial
    from sage.symbolic.ring import SR
    from sage.rings.qqbar import QQbar
    from sage.functions.generalized import kronecker_delta

    got_SR = A.base_ring() == SR

    # Change to QQbar if possible
    try:
        A = A.change_ring(QQbar)
    except (TypeError, NotImplementedError):
        pass

    # Get Jordan matrix J and invertible matrix P such that A = P*J*~P
    # From that, we will compute M = J^n, and obtain A^n = P*J^n*~P
    J, P = A.jordan_form(transformation=True)

    # Where each Jordan block starts, and number of blocks
    block_start = [0] + J.subdivisions()[0]
    num_blocks = len(block_start)

    # Prepare matrix M to store `J^n`, computed by Jordan block
    M = matrix(SR, J.ncols())
    M.subdivide(J.subdivisions())

    for k in range(num_blocks):

        # Jordan block Jk, its dimension nk, the eigenvalue m
        Jk = J.subdivision(k, k)
        nk = Jk.ncols()
        mk = Jk[0, 0]

        # First row of block Mk; its entries are of the form
        # D^i(f) / i! with f = x^n and D = differentiation wrt x
        if hasattr(mk, 'radical_expression'):
            mk = mk.radical_expression()

        # When the variable "mk" is equal to zero, it is advisable to employ the Kronecker delta function
        # instead of utilizing the numerical value zero. This choice is made to encompass scenarios where
        # the power of zero is also equal to zero.
        if mk:
            vk = [(binomial(n, i) * mk._pow_(n-i)).simplify_full()
                  for i in range(nk)]
        else:
            vk = [(binomial(n, i).simplify_full() * kronecker_delta(n,i))
                  for i in range(nk)]

        # Form block Mk and insert it in M
        Mk = matrix(SR, [[SR.zero()]*i + vk[:nk-i] for i in range(nk)])
        M.set_block(block_start[k], block_start[k], Mk)

    # Change entries of P and P^-1 into symbolic expressions
    if not got_SR:
        Pinv = (~P).apply_map(AlgebraicNumber.radical_expression)
        P = P.apply_map(AlgebraicNumber.radical_expression)
    else:
        Pinv = ~P

    return P * M * Pinv


class NotFullRankError(ValueError):
    """
    An error that indicates that a matrix is not of full rank.

    The fact that a square system is rank-deficient sometimes only becomes
    apparent while attempting to solve it. The methods
    :meth:`.Matrix.solve_left` and :meth:`.Matrix.solve_right` defer to
    :meth:`.Matrix._solve_right_nonsingular_square` for square systems, and
    that method raises this error if the system turns out to be singular.
    """
    pass


cdef inline bint _block_ldlt_pivot1x1(Matrix A, Py_ssize_t k) except 1:
    r"""
    Update the `n`-by-`n` matrix `A` as part of a 1x1 pivot in the
    ``k,k`` position (whose value is ``pivot``). Relies on the fact
    that `A` is passed in by reference, since for performance reasons
    this routine should overwrite its argument.

    There is no return value from this function, as its intended
    effect is to update the matrix `A` in-place, but we allow it
    to return zero/one so that ``1`` can be used to indicate that
    a python exception occurred.
    """
    cdef Py_ssize_t i, j  # dumy loop indices
    cdef Py_ssize_t n = A._nrows
    pivot = A.get_unsafe(k, k)

    # Compute the Schur complement that we'll work on during
    # the following iteration, and store it back in the lower-
    # right-hand corner of "A".
    for i in range(n-k-1):
        for j in range(i+1):
            A.set_unsafe(k+1+i,
                         k+1+j,
                         (A.get_unsafe(k+1+i, k+1+j) -
                          A.get_unsafe(k+1+i, k)*A.get_unsafe(k, k+1+j)/pivot))
            A.set_unsafe(k+1+j,
                         k+1+i,
                         A.get_unsafe(k+1+i, k+1+j).conjugate())

    for i in range(n-k-1):
        # Store the new (kth) column of "L" within the lower-
        # left-hand corner of "A".
        A.set_unsafe(k+i+1,
                     k,
                     A.get_unsafe(k+i+1, k)/ pivot)

    return 0
