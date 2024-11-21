# sage.doctest: optional - numpy
"""
Dense matrices using a NumPy backend

This serves as a base class for dense matrices over
Real Double Field and Complex Double Field.

AUTHORS:

- Jason Grout, Sep 2008: switch to NumPy backend, factored out the Matrix_double_dense class

- Josh Kantor

- William Stein: many bug fixes and touch ups.

EXAMPLES::

    sage: b = Mat(RDF,2,3).basis()
    sage: b[0,0]
    [1.0 0.0 0.0]
    [0.0 0.0 0.0]

We deal with the case of zero rows or zero columns::

    sage: m = MatrixSpace(RDF,0,3)
    sage: m.zero_matrix()
    []

TESTS::

    sage: a = matrix(RDF,2,range(4), sparse=False)
    sage: TestSuite(a).run()
    sage: a = matrix(CDF,2,range(4), sparse=False)
    sage: TestSuite(a).run()
"""

# ****************************************************************************
#       Copyright (C) 2004-2006 Joshua Kantor <kantor.jm@gmail.com>
#       Copyright (C) 2008      Georg S. Weber
#       Copyright (C) 2008-2011 Mike Hansen
#       Copyright (C) 2008-2012 Jason Grout
#       Copyright (C) 2009      Dag Sverre Seljebotn
#       Copyright (C) 2009      Yann Laigle-Chapuy
#       Copyright (C) 2009-2010 Florent Hivert
#       Copyright (C) 2010-2012 Rob Beezer
#       Copyright (C) 2011      Martin Raum
#       Copyright (C) 2011-2012 J. H. Palmieri
#       Copyright (C) 2011-2014 André Apitzsch
#       Copyright (C) 2011-2018 Jeroen Demeyer
#       Copyright (C) 2012      Kenneth Smith
#       Copyright (C) 2016-2019 Frédéric Chapoton
#       Copyright (C) 2017      Kiran Kedlaya
#       Copyright (C) 2019      Chaman Agrawal
#       Copyright (C) 2019-2021 Markus Wageringel
#       Copyright (C) 2020      Michael Orlitzky
#       Copyright (C) 2020      Victor Santos
#       Copyright (C) 2021      Jonathan Kliem
#       Copyright (C) 2021      Travis Scrimshaw
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import math

import sage.rings.real_double
import sage.rings.complex_double

from sage.structure.element cimport Vector
from sage.matrix.constructor import matrix
cimport sage.structure.element

cimport numpy as cnumpy

numpy = None
scipy = None

# This is for the Numpy C API to work
cnumpy.import_array()


cdef class Matrix_double_dense(Matrix_numpy_dense):
    """
    Base class for matrices over the Real Double Field and the Complex
    Double Field.  These are supposed to be fast matrix operations
    using C doubles. Most operations are implemented using numpy which
    will call the underlying BLAS on the system.

    This class cannot be instantiated on its own.  The numpy matrix
    creation depends on several variables that are set in the
    subclasses.

    EXAMPLES::

        sage: m = Matrix(RDF, [[1,2],[3,4]])
        sage: m**2
        [ 7.0 10.0]
        [15.0 22.0]
        sage: m^(-1)        # rel tol 1e-15                                             # needs scipy
        [-1.9999999999999996  0.9999999999999998]
        [ 1.4999999999999998 -0.4999999999999999]

    TESTS:

    Test hashing::

        sage: A = matrix(RDF, 3, range(1,10))
        sage: hash(A)
        Traceback (most recent call last):
        ...
        TypeError: mutable matrices are unhashable
        sage: A.set_immutable()
        sage: hash(A)
        6694819972852100501  # 64-bit
        1829383573           # 32-bit
        sage: A = matrix(CDF, 3, range(1,10))
        sage: hash(A)
        Traceback (most recent call last):
        ...
        TypeError: mutable matrices are unhashable
        sage: A.set_immutable()
        sage: hash(A)
        6694819972852100501  # 64-bit
        1829383573           # 32-bit
    """

    def LU_valid(self):
        r"""
        Return ``True`` if the LU form of this matrix has
        already been computed.

        EXAMPLES::

            sage: A = random_matrix(RDF,3) ; A.LU_valid()
            False
            sage: P, L, U = A.LU()
            sage: A.LU_valid()
            True
        """
        return self.fetch('PLU_factors') is not None

    ########################################################################
    # LEVEL 2 functionality
    #   * def _pickle
    #   * def _unpickle
    cpdef _add_(self, right):
        """
        Add two matrices together.

        EXAMPLES::

            sage: A = matrix(RDF,3,range(1,10))
            sage: A+A
            [ 2.0  4.0  6.0]
            [ 8.0 10.0 12.0]
            [14.0 16.0 18.0]
        """
        if self._nrows == 0 or self._ncols == 0:
            return self.__copy__()

        cdef Matrix_double_dense M, _right, _left
        _right = right
        _left = self

        M = self._new()
        M._matrix_numpy = _left._matrix_numpy + _right._matrix_numpy
        return M

    cpdef _sub_(self, right):
        """
        Return ``self - right``.

        EXAMPLES::

            sage: A = matrix(RDF,3,range(1,10))
            sage: (A-A).is_zero()
            True
        """
        if self._nrows == 0 or self._ncols == 0:
            return self.__copy__()

        cdef Matrix_double_dense M,_right,_left
        _right = right
        _left = self

        M = self._new()
        M._matrix_numpy = _left._matrix_numpy - _right._matrix_numpy
        return M

    def __neg__(self):
        """
        Negate this matrix.

        EXAMPLES::

            sage: A = matrix(RDF,3,range(1,10))
            sage: -A
            [-1.0 -2.0 -3.0]
            [-4.0 -5.0 -6.0]
            [-7.0 -8.0 -9.0]
            sage: B = -A ; (A+B).is_zero()
            True
        """
        if self._nrows == 0 or self._ncols == 0:
            return self.__copy__()

        cdef Matrix_double_dense M
        M = self._new()
        M._matrix_numpy = -self._matrix_numpy
        return M

    # x * __copy__
    #   * _list -- list of underlying elements (need not be a copy)
    #   * _dict -- sparse dictionary of underlying elements (need not be a copy)
    ########################################################################
    # def _pickle(self):                        #unsure how to implement
    # def _unpickle(self, data, int version):   # use version >= 0 #unsure how to implement
    ######################################################################
    cdef sage.structure.element.Matrix _matrix_times_matrix_(self, sage.structure.element.Matrix right):
        r"""
        Multiply ``self * right`` as matrices.

        EXAMPLES::

            sage: A = matrix(RDF,3,range(1,10))
            sage: B = matrix(RDF,3,range(1,13))
            sage: A*B
            [ 38.0  44.0  50.0  56.0]
            [ 83.0  98.0 113.0 128.0]
            [128.0 152.0 176.0 200.0]

        TESTS:

        Check that :issue:`31234` is fixed::

            sage: matrix.identity(QQ, 4) * matrix(RDF, 4, 0)
            []

        Check that an empty matrix is initialized correctly; see :issue:`27366`:

            sage: A = matrix(RDF, 3, 0)
            sage: A*A.transpose()
            [0.0 0.0 0.0]
            [0.0 0.0 0.0]
            [0.0 0.0 0.0]
        """
        if self._ncols != right._nrows:
            raise IndexError("Number of columns of self must equal number of rows of right")

        cdef Matrix_double_dense M, _right, _left

        if self._nrows == 0 or self._ncols == 0 or right._nrows == 0 or right._ncols == 0:
            M = self._new(self._nrows, right._ncols)
            M._matrix_numpy.fill(0)
            return M

        M = self._new(self._nrows, right._ncols)
        _right = right
        _left = self
        global numpy
        if numpy is None:
            import numpy

        M._matrix_numpy = numpy.dot(_left._matrix_numpy, _right._matrix_numpy)
        return M

    def __invert__(self):
        """
        Invert this matrix.

        EXAMPLES::

            sage: A = Matrix(RDF, [[10, 0], [0, 100]])
            sage: (~A).det()
            0.001

            sage: A = matrix(RDF,3,[2,3,5,7,8,9,11,13,17]); A
            [ 2.0  3.0  5.0]
            [ 7.0  8.0  9.0]
            [11.0 13.0 17.0]
            sage: ~A  # tol 1e-14
            [-2.7142857142857184  -2.000000000000004  1.8571428571428603]
            [  2.857142857142863   3.000000000000006 -2.4285714285714333]
            [-0.4285714285714305  -1.000000000000002  0.7142857142857159]

        Note that if this matrix is (nearly) singular, finding
        its inverse will not help much and will give slightly different
        answers on similar platforms depending on the hardware
        and other factors::

            sage: A = matrix(RDF,3,range(1,10));A
            [1.0 2.0 3.0]
            [4.0 5.0 6.0]
            [7.0 8.0 9.0]

            sage: A.determinant() < 10e-12
            True

        TESTS::

            sage: ~Matrix(RDF, 0,0)
            []
            sage: ~Matrix(RDF, 0,3)
            Traceback (most recent call last):
            ...
            ArithmeticError: self must be a square matrix
        """
# see github issue 4502 --- there is an issue with the "#random" pragma that needs to be fixed
#                          as for the mathematical side, scipy v0.7 is expected to fix the invertibility failures
#
#            sage: A = Matrix(RDF, [[1, 0], [0, 0]])
#            sage: A.inverse().det()        # random - on some computers, this will be invertible due to numerical error.
#            Traceback (most recent call last):
#            ...
#            LinAlgError: singular matrix
#            sage: A = matrix(RDF,3,range(1,10));A
#            [1.0 2.0 3.0]
#            [4.0 5.0 6.0]
#            [7.0 8.0 9.0]
#
#            sage: A.determinant() < 10e-12
#            True
#            sage: ~A                       # random - on some computers, this will be invertible due to numerical error.
#            Traceback (most recent call last):
#            ...
#            ZeroDivisionError: singular matrix
#
        if self._nrows != self._ncols:
            raise ArithmeticError("self must be a square matrix")
        if self._nrows == 0 and self._ncols == 0:
            return self.__copy__()

        # Maybe we should cache the (P)LU decomposition and use scipy.lu_solve?
        cdef Matrix_double_dense M
        M = self._new()
        global scipy
        if scipy is None:
            import scipy
        import scipy.linalg
        from numpy.linalg import LinAlgError
        try:  # Standard error reporting for Sage.
            M._matrix_numpy = scipy.linalg.inv(self._matrix_numpy)
        except LinAlgError:
            raise ZeroDivisionError("input matrix must be nonsingular")
        return M

    # def _list(self):
    # def _dict(self):

    ########################################################################
    # LEVEL 3 functionality (Optional)
    #    * cdef _sub_
    #    * __deepcopy__
    #    * __invert__
    #    * Matrix windows -- only if you need strassen for that base
    #    * Other functions (list them here):
    #
    #    compute_LU(self)
    #
    ########################################################################

    def condition(self, p='frob'):
        r"""
        Return the condition number of a square nonsingular matrix.

        Roughly speaking, this is a measure of how sensitive
        the matrix is to round-off errors in numerical computations.
        The minimum possible value is 1.0, and larger numbers indicate
        greater sensitivity.

        INPUT:

        - ``p`` -- (default: ``'frob'``) controls which norm is used
          to compute the condition number, allowable values are
          'frob' (for the Frobenius norm), integers -2, -1, 1, 2,
          positive and negative infinity. See output discussion
          for specifics.

        OUTPUT:

        The condition number of a matrix is the product of a norm
        of the matrix times the norm of the inverse of the matrix.
        This requires that the matrix be square and invertible
        (nonsingular, full rank).

        Returned value is a double precision floating point value
        in ``RDF``, or ``Infinity``.  Row and column sums described below are
        sums of the absolute values of the entries, where the
        absolute value of the complex number `a+bi` is `\sqrt{a^2+b^2}`.
        Singular values are the "diagonal" entries of the "S" matrix in
        the singular value decomposition.

        - ``p = 'frob'``: the default norm employed in computing
          the condition number, the Frobenius norm, which for a
          matrix `A=(a_{ij})` computes

          .. MATH::

                \left(\sum_{i,j}\left\lvert{a_{i,j}}\right\rvert^2\right)^{1/2}

        - ``p = 'sv'``: the quotient of the maximal and minimal singular value.
        - ``p = Infinity`` or ``p = oo``: the maximum row sum.
        - ``p = -Infinity`` or ``p = -oo``: the minimum column sum.
        - ``p = 1``: the maximum column sum.
        - ``p = -1``: the minimum column sum.
        - ``p = 2``: the 2-norm, equal to the maximum singular value.
        - ``p = -2``: the minimum singular value.

        ALGORITHM:

        Computation is performed by the ``cond()`` function of
        the SciPy/NumPy library.

        EXAMPLES:

        First over the reals.  ::

            sage: A = matrix(RDF, 4, [(1/4)*x^3 for x in range(16)]); A
            [   0.0   0.25    2.0   6.75]
            [  16.0  31.25   54.0  85.75]
            [ 128.0 182.25  250.0 332.75]
            [ 432.0 549.25  686.0 843.75]
            sage: A.condition()
            9923.88955...
            sage: A.condition(p='frob')
            9923.88955...
            sage: A.condition(p=Infinity)  # tol 3e-14
            22738.50000000045
            sage: A.condition(p=-Infinity)  # tol 2e-14
            17.50000000000028
            sage: A.condition(p=1)
            12139.21...
            sage: A.condition(p=-1)  # tol 2e-14
            550.0000000000093
            sage: A.condition(p=2)
            9897.8088...
            sage: A.condition(p=-2)
            0.000101032462...

        And over the complex numbers.  ::

            sage: B = matrix(CDF, 3, [x + x^2*I for x in range(9)]); B
            [         0.0  1.0 + 1.0*I  2.0 + 4.0*I]
            [ 3.0 + 9.0*I 4.0 + 16.0*I 5.0 + 25.0*I]
            [6.0 + 36.0*I 7.0 + 49.0*I 8.0 + 64.0*I]
            sage: B.condition()
            203.851798...
            sage: B.condition(p='frob')
            203.851798...
            sage: B.condition(p=Infinity)
            369.55630...
            sage: B.condition(p=-Infinity)
            5.46112969...
            sage: B.condition(p=1)
            289.251481...
            sage: B.condition(p=-1)
            20.4566639...
            sage: B.condition(p=2)
            202.653543...
            sage: B.condition(p=-2)
            0.00493453005...

        Hilbert matrices are famously ill-conditioned, while
        an identity matrix can hit the minimum with the right norm.  ::

            sage: A = matrix(RDF, 10, [1/(i+j+1) for i in range(10) for j in range(10)])
            sage: A.condition()  # tol 2e-4
            16332197709146.014
            sage: id = identity_matrix(CDF, 10)
            sage: id.condition(p=1)
            1.0

        Return values are in `RDF`.  ::

            sage: A = matrix(CDF, 2, range(1,5))
            sage: A.condition() in RDF
            True

        Rectangular and singular matrices raise errors if p is not 'sv'.  ::

            sage: A = matrix(RDF, 2, 3, range(6))
            sage: A.condition()
            Traceback (most recent call last):
            ...
            TypeError: matrix must be square if p is not 'sv', not 2 x 3

            sage: A.condition('sv')
            7.34...

            sage: A = matrix(QQ, 5, range(25))
            sage: A.is_singular()
            True
            sage: B = A.change_ring(CDF)
            sage: B.condition()
            +Infinity

        Improper values of ``p`` are caught.  ::

            sage: A = matrix(CDF, 2, range(1,5))
            sage: A.condition(p='bogus')
            Traceback (most recent call last):
            ...
            ValueError: condition number 'p' must be +/- infinity, 'frob', 'sv' or an integer, not bogus
            sage: A.condition(p=632)
            Traceback (most recent call last):
            ...
            ValueError: condition number integer values of 'p' must be -2, -1, 1 or 2, not 632

        TESTS:

        Some condition numbers, first by the definition which also exercises
        :meth:`norm`, then by this method.  ::

            sage: A = matrix(CDF, [[1,2,4],[5,3,9],[7,8,6]])
            sage: c = A.norm(2)*A.inverse().norm(2)
            sage: d = A.condition(2)
            sage: abs(c-d) < 1.0e-12
            True
            sage: c = A.norm(1)*A.inverse().norm(1)
            sage: d = A.condition(1)
            sage: abs(c-d) < 1.0e-12
            True
        """
        if not self.is_square() and p != 'sv':
            raise TypeError("matrix must be square if p is not 'sv', not %s x %s" % (self.nrows(), self.ncols()))
        global numpy
        if numpy is None:
            import numpy
        import sage.rings.infinity
        import sage.rings.integer
        from sage.rings.real_double import RDF
        if p == sage.rings.infinity.Infinity:
            p = numpy.inf
        elif p == -sage.rings.infinity.Infinity:
            p = -numpy.inf
        elif p == 'frob':
            p = 'fro'
        elif p == 'sv':
            p = None
        else:
            try:
                p = sage.rings.integer.Integer(p)
            except TypeError:
                raise ValueError("condition number 'p' must be +/- infinity, 'frob', 'sv' or an integer, not %s" % p)
            if p not in [-2, -1, 1, 2]:
                raise ValueError("condition number integer values of 'p' must be -2, -1, 1 or 2, not %s" % p)
        # may raise a LinAlgError if matrix is singular
        c = numpy.linalg.cond(self._matrix_numpy, p=p)
        if c == numpy.inf:
            return sage.rings.infinity.Infinity
        else:
            return RDF(c.real if numpy.iscomplexobj(c) else c)

    def norm(self, p=2):
        r"""
        Return the norm of the matrix.

        INPUT:

        - ``p`` -- (default: 2) controls which norm is computed,
          allowable values are 'frob' (for the Frobenius norm),
          integers -2, -1, 1, 2, positive and negative infinity.  See
          output discussion for specifics.

        OUTPUT:

        Returned value is a double precision floating point value
        in ``RDF``.  Row and column sums described below are
        sums of the absolute values of the entries, where the
        absolute value of the complex number `a+bi` is `\sqrt{a^2+b^2}`.
        Singular values are the "diagonal" entries of the "S" matrix in
        the singular value decomposition.

        - ``p = 'frob'``: the Frobenius norm, which for
          a matrix `A=(a_{ij})` computes

          .. MATH::

                \left(\sum_{i,j}\left\lvert{a_{i,j}}\right\rvert^2\right)^{1/2}

        - ``p = Infinity`` or ``p = oo``: the maximum row sum.
        - ``p = -Infinity`` or ``p = -oo``: the minimum column sum.
        - ``p = 1``: the maximum column sum.
        - ``p = -1``: the minimum column sum.
        - ``p = 2``: the induced 2-norm, equal to the maximum singular value.
        - ``p = -2``: the minimum singular value.

        ALGORITHM:

        Computation is performed by the :func:`~scipy:scipy.linalg.norm`
        function of the SciPy/NumPy library.

        EXAMPLES:

        First over the reals.  ::

            sage: A = matrix(RDF, 3, range(-3, 6)); A
            [-3.0 -2.0 -1.0]
            [ 0.0  1.0  2.0]
            [ 3.0  4.0  5.0]
            sage: A.norm()
            7.99575670...
            sage: A.norm(p='frob')
            8.30662386...
            sage: A.norm(p=Infinity)
            12.0
            sage: A.norm(p=-Infinity)
            3.0
            sage: A.norm(p=1)
            8.0
            sage: A.norm(p=-1)
            6.0
            sage: A.norm(p=2)
            7.99575670...
            sage: A.norm(p=-2) < 10^-15
            True

        And over the complex numbers.  ::

            sage: # needs sage.symbolic
            sage: B = matrix(CDF, 2, [[1+I, 2+3*I],[3+4*I,3*I]]); B
            [1.0 + 1.0*I 2.0 + 3.0*I]
            [3.0 + 4.0*I       3.0*I]
            sage: B.norm()
            6.66189877...
            sage: B.norm(p='frob')
            7.0
            sage: B.norm(p=Infinity)
            8.0
            sage: B.norm(p=-Infinity)
            5.01976483...
            sage: B.norm(p=1)
            6.60555127...
            sage: B.norm(p=-1)
            6.41421356...
            sage: B.norm(p=2)
            6.66189877...
            sage: B.norm(p=-2)
            2.14921023...

        Since it is invariant under unitary multiplication, the
        Frobenius norm is equal to the square root of the sum of
        squares of the singular values.  ::

            sage: A = matrix(RDF, 5, range(1,26))
            sage: f = A.norm(p='frob')
            sage: U, S, V = A.SVD()
            sage: s = sqrt(sum([S[i,i]^2 for i in range(5)]))
            sage: abs(f-s) < 1.0e-12
            True

        Return values are in `RDF`. ::

            sage: A = matrix(CDF, 2, range(4))
            sage: A.norm() in RDF
            True

        Improper values of ``p`` are caught.  ::

            sage: A.norm(p='bogus')
            Traceback (most recent call last):
            ...
            ValueError: matrix norm 'p' must be +/- infinity, 'frob' or an integer, not bogus
            sage: A.norm(p=632)
            Traceback (most recent call last):
            ...
            ValueError: matrix norm integer values of 'p' must be -2, -1, 1 or 2, not 632
        """
        global numpy
        if numpy is None:
            import numpy

        import sage.rings.infinity
        import sage.rings.integer
        import sage.rings.real_double
        if p == sage.rings.infinity.Infinity:
            p = numpy.inf
        elif p == -sage.rings.infinity.Infinity:
            p = -numpy.inf
        elif p == 'frob':
            p = 'fro'
        else:
            try:
                p = sage.rings.integer.Integer(p)
            except TypeError:
                raise ValueError("matrix norm 'p' must be +/- infinity, 'frob' or an integer, not %s" % p)
            if p not in [-2, -1, 1, 2]:
                raise ValueError("matrix norm integer values of 'p' must be -2, -1, 1 or 2, not %s" % p)
        return sage.rings.real_double.RDF(numpy.linalg.norm(self._matrix_numpy, ord=p))

    def singular_values(self, eps=None):
        r"""
        Return a sorted list of the singular values of the matrix.

        INPUT:

        - ``eps`` -- (default: ``None``) the largest number which
          will be considered to be zero.  May also be set to the
          string 'auto'.  See the discussion below.

        OUTPUT:

        A sorted list of the singular values of the matrix, which are the
        diagonal entries of the "S" matrix in the SVD decomposition.  As such,
        the values are real and are returned as elements of ``RDF``.  The
        list is sorted with larger values first, and since theory predicts
        these values are always positive, for a rank-deficient matrix the
        list should end in zeros (but in practice may not).  The length of
        the list is the minimum of the row count and column count for the
        matrix.

        The number of nonzero singular values will be the rank of the
        matrix.  However, as a numerical matrix, it is impossible to
        control the difference between zero entries and very small
        nonzero entries.  As an informed consumer it is up to you
        to use the output responsibly.  We will do our best, and give
        you the tools to work with the output, but we cannot
        give you a guarantee.

        With ``eps`` set to ``None`` you will get the raw singular
        values and can manage them as you see fit.  You may also set
        ``eps`` to any positive floating point value you wish.  If you
        set ``eps`` to 'auto' this routine will compute a reasonable
        cutoff value, based on the size of the matrix, the largest
        singular value and the smallest nonzero value representable
        by the 53-bit precision values used.  See the discussion
        at page 268 of [Wat2010]_.

        See the examples for a way to use the "verbose" facility
        to easily watch the zero cutoffs in action.

        ALGORITHM:

        The singular values come from the SVD decomposition
        computed by SciPy/NumPy using :func:`scipy:scipy.linalg.svd`.

        EXAMPLES:

        Singular values close to zero have trailing digits that may vary
        on different hardware.  For exact matrices, the number of nonzero
        singular values will equal the rank of the matrix.  So for some of
        the doctests we round the small singular values that ideally would
        be zero, to control the variability across hardware.

        This matrix has a determinant of one.  A chain of two or
        three theorems implies the product of the singular values
        must also be one.  ::

            sage: A = matrix(QQ, [[ 1,  0,  0,  0,  0,  1,  3],
            ....:                 [-2,  1,  1, -2,  0, -4,  0],
            ....:                 [ 1,  0,  1, -4, -6, -3,  7],
            ....:                 [-2,  2,  1,  1,  7,  1, -1],
            ....:                 [-1,  0, -1,  5,  8,  4, -6],
            ....:                 [ 4, -2, -2,  1, -3,  0,  8],
            ....:                 [-2,  1,  0,  2,  7,  3, -4]])
            sage: A.determinant()
            1
            sage: B = A.change_ring(RDF)
            sage: sv = B.singular_values(); sv  # tol 1e-12
            [20.523980658874265, 8.486837028536643, 5.86168134845073, 2.4429165899286978, 0.5831970144724045, 0.26933287286576313, 0.0025524488076110402]
            sage: prod(sv)  # tol 1e-12
            0.9999999999999525

        An exact matrix that is obviously not of full rank, and then
        a computation of the singular values after conversion
        to an approximate matrix. ::

            sage: A = matrix(QQ, [[1/3, 2/3, 11/3],
            ....:                 [2/3, 1/3,  7/3],
            ....:                 [2/3, 5/3, 27/3]])
            sage: A.rank()
            2
            sage: B = A.change_ring(CDF)
            sage: sv = B.singular_values()
            sage: sv[0:2]
            [10.1973039..., 0.487045871...]
            sage: sv[2] < 1e-14
            True

        A matrix of rank 3 over the complex numbers.  ::

            sage: A = matrix(CDF, [[46*I - 28, -47*I - 50, 21*I + 51, -62*I - 782, 13*I + 22],
            ....:                  [35*I - 20, -32*I - 46, 18*I + 43, -57*I - 670, 7*I + 3],
            ....:                  [22*I - 13, -23*I - 23, 9*I + 24, -26*I - 347, 7*I + 13],
            ....:                  [-44*I + 23, 41*I + 57, -19*I - 54, 60*I + 757, -11*I - 9],
            ....:                  [30*I - 18, -30*I - 34, 14*I + 34, -42*I - 522, 8*I + 12]])
            sage: sv = A.singular_values()
            sage: sv[0:3]  # tol 1e-14
            [1440.7336659952966, 18.404403413369227, 6.839707797136151]
            sage: (sv[3] < 10^-13) or sv[3]
            True
            sage: (sv[4] < 10^-14) or sv[4]
            True

        A full-rank matrix that is ill-conditioned.  We use this to
        illustrate ways of using the various possibilities for ``eps``,
        including one that is ill-advised. Notice that the automatically
        computed cutoff gets this (difficult) example slightly wrong.
        This illustrates the impossibility of any automated process always
        getting this right.  Use with caution and judgement.  ::

            sage: entries = [1/(i+j+1) for i in range(12) for j in range(12)]
            sage: B = matrix(QQ, 12, 12, entries)
            sage: B.rank()
            12
            sage: A = B.change_ring(RDF)
            sage: A.condition() > 1.59e16 or A.condition()
            True

            sage: A.singular_values(eps=None)  # abs tol 7e-16
            [1.7953720595619975, 0.38027524595503703, 0.04473854875218107, 0.0037223122378911614, 0.0002330890890217751, 1.116335748323284e-05, 4.082376110397296e-07, 1.1228610675717613e-08, 2.2519645713496478e-10, 3.1113486853814003e-12, 2.6500422260778388e-14, 9.87312834948426e-17]
            sage: A.singular_values(eps='auto')  # abs tol 7e-16
            [1.7953720595619975, 0.38027524595503703, 0.04473854875218107, 0.0037223122378911614, 0.0002330890890217751, 1.116335748323284e-05, 4.082376110397296e-07, 1.1228610675717613e-08, 2.2519645713496478e-10, 3.1113486853814003e-12, 2.6500422260778388e-14, 0.0]
            sage: A.singular_values(eps=1e-4)  # abs tol 7e-16
            [1.7953720595619975, 0.38027524595503703, 0.04473854875218107, 0.0037223122378911614, 0.0002330890890217751, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

        With Sage's "verbose" facility, you can compactly see the cutoff
        at work.  In any application of this routine, or those that build upon
        it, it would be a good idea to conduct this exercise on samples.
        We also test here that all the  values are returned in `RDF` since
        singular values are always real. ::

            sage: A = matrix(CDF, 4, range(16))
            sage: from sage.misc.verbose import set_verbose
            sage: set_verbose(1)
            sage: sv = A.singular_values(eps='auto'); sv
            verbose 1 (<module>) singular values,
            smallest-non-zero:cutoff:largest-zero,
            2.2766...:6.2421...e-14:...
            [35.13996365902..., 2.27661020871472..., 0.0, 0.0]
            sage: set_verbose(0)

            sage: all(s in RDF for s in sv)
            True

        TESTS:

        Bogus values of the ``eps`` keyword will be caught::

            sage: A.singular_values(eps='junk')
            Traceback (most recent call last):
            ...
            ValueError: could not convert string to float: ...

        AUTHOR:

        - Rob Beezer - (2011-02-18)
        """
        from sage.misc.verbose import verbose
        from sage.rings.real_double import RDF
        global scipy
        # get SVD decomposition, which is a cached quantity
        _, S, _ = self.SVD()
        diag = min(self._nrows, self._ncols)
        sv = [RDF(S[i,i]) for i in range(diag)]
        # no cutoff, send raw data back
        if eps is None:
            verbose("singular values, no zero cutoff specified", level=1)
            return sv
        # set cutoff as RDF element
        if eps == 'auto':
            if scipy is None: import scipy
            eps = 2*max(self._nrows, self._ncols)*numpy.finfo(float).eps*sv[0]
        eps = RDF(eps)
        # locate nonzero entries
        rank = 0
        while rank < diag and sv[rank] > eps:
            rank = rank + 1
        # capture info for watching zero cutoff behavior at verbose level 1
        if rank == 0:
            small_nonzero = None
        else:
            small_nonzero = sv[rank-1]
        if rank < diag:
            large_zero = sv[rank]
        else:
            large_zero = None
        # convert small values to zero, then done
        for i in range(rank, diag):
            sv[i] = RDF(0)
        verbose("singular values, smallest-non-zero:cutoff:largest-zero, %s:%s:%s" % (small_nonzero, eps, large_zero), level=1)
        return sv

    def LU(self):
        r"""
        Return a decomposition of the (row-permuted) matrix as a product of
        a lower-triangular matrix ("L") and an upper-triangular matrix ("U").

        OUTPUT:

        For an `m\times n` matrix ``A`` this method returns a triple of
        immutable matrices ``P, L, U`` such that

        - ``A = P*L*U``
        - ``P`` is a square permutation matrix, of size `m\times m`,
          so is all zeroes, but with exactly a single one in each
          row and each column
        - ``L`` is lower-triangular, square of size `m\times m`,
          with every diagonal entry equal to one
        - ``U`` is upper-triangular with size `m\times n`, i.e.
          entries below the "diagonal" are all zero

        The computed decomposition is cached and returned on
        subsequent calls, thus requiring the results to be immutable.

        Effectively, ``P`` permutes the rows of ``A``.  Then ``L``
        can be viewed as a sequence of row operations on this matrix,
        where each operation is adding a multiple of a row to a
        subsequent row.  There is no scaling (thus 1s on the diagonal
        of ``L``) and no row-swapping (``P`` does that).  As a result
        ``U`` is close to being the result of Gaussian-elimination.
        However, round-off errors can make it hard to determine
        the zero entries of ``U``.

        .. NOTE::
            The behaviour of ``LU()`` has changed in Sage version 9.1.
            Earlier, ``LU()`` returned ``P,L,U`` such that ``P*A=L*U``,
            where ``P`` represents the permutation and is
            the matrix inverse of the ``P`` returned by this method.
            The computation of this matrix inverse can be accomplished
            quickly with just a transpose as the matrix is orthogonal/unitary.

            For details see :issue:`18365`.

        EXAMPLES::

            sage: m = matrix(RDF,4,range(16))
            sage: P,L,U = m.LU()
            sage: P*L*U # rel tol 2e-16
            [ 0.0  1.0  2.0  3.0]
            [ 4.0  5.0  6.0  7.0]
            [ 8.0  9.0 10.0 11.0]
            [12.0 13.0 14.0 15.0]

        Below example illustrates the change in behaviour of ``LU()``. ::

            sage: (m - P*L*U).norm() < 1e-14
            True
            sage: (P*m - L*U).norm() < 1e-14
            False

        :issue:`10839` made this routine available for rectangular matrices.  ::

            sage: A = matrix(RDF, 5, 6, range(30)); A
            [ 0.0  1.0  2.0  3.0  4.0  5.0]
            [ 6.0  7.0  8.0  9.0 10.0 11.0]
            [12.0 13.0 14.0 15.0 16.0 17.0]
            [18.0 19.0 20.0 21.0 22.0 23.0]
            [24.0 25.0 26.0 27.0 28.0 29.0]
            sage: P, L, U = A.LU()
            sage: P
            [0.0 1.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0 1.0]
            [0.0 0.0 1.0 0.0 0.0]
            [0.0 0.0 0.0 1.0 0.0]
            [1.0 0.0 0.0 0.0 0.0]
            sage: L.zero_at(0)   # Use zero_at(0) to get rid of signed zeros
            [ 1.0  0.0  0.0  0.0  0.0]
            [ 0.0  1.0  0.0  0.0  0.0]
            [ 0.5  0.5  1.0  0.0  0.0]
            [0.75 0.25  0.0  1.0  0.0]
            [0.25 0.75  0.0  0.0  1.0]
            sage: U.zero_at(0)   # Use zero_at(0) to get rid of signed zeros
            [24.0 25.0 26.0 27.0 28.0 29.0]
            [ 0.0  1.0  2.0  3.0  4.0  5.0]
            [ 0.0  0.0  0.0  0.0  0.0  0.0]
            [ 0.0  0.0  0.0  0.0  0.0  0.0]
            [ 0.0  0.0  0.0  0.0  0.0  0.0]
            sage: P.transpose()*A-L*U
            [0.0 0.0 0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0 0.0 0.0]
            sage: P*L*U
            [ 0.0  1.0  2.0  3.0  4.0  5.0]
            [ 6.0  7.0  8.0  9.0 10.0 11.0]
            [12.0 13.0 14.0 15.0 16.0 17.0]
            [18.0 19.0 20.0 21.0 22.0 23.0]
            [24.0 25.0 26.0 27.0 28.0 29.0]

        Trivial cases return matrices of the right size and
        characteristics.  ::

            sage: A = matrix(RDF, 5, 0)
            sage: P, L, U = A.LU()
            sage: P.parent()
            Full MatrixSpace of 5 by 5 dense matrices over Real Double Field
            sage: L.parent()
            Full MatrixSpace of 5 by 5 dense matrices over Real Double Field
            sage: U.parent()
            Full MatrixSpace of 5 by 0 dense matrices over Real Double Field
            sage: A-P*L*U
            []

        The results are immutable since they are cached.  ::

            sage: P, L, U = matrix(RDF, 2, 2, range(4)).LU()
            sage: L[0,0] = 0
            Traceback (most recent call last):
                ...
            ValueError: matrix is immutable; please change a copy instead (i.e., use copy(M) to change a copy of M).
            sage: P[0,0] = 0
            Traceback (most recent call last):
                ...
            ValueError: matrix is immutable; please change a copy instead (i.e., use copy(M) to change a copy of M).
            sage: U[0,0] = 0
            Traceback (most recent call last):
                ...
            ValueError: matrix is immutable; please change a copy instead (i.e., use copy(M) to change a copy of M).
        """
        global scipy, numpy
        cdef Matrix_double_dense P, L, U
        m = self._nrows
        n = self._ncols

        # scipy fails on trivial cases
        if m == 0 or n == 0:
            P = self._new(m, m)
            for i in range(m):
                P[i,i]=1
            P.set_immutable()
            L = P
            U = self._new(m,n)
            U.set_immutable()
            return P, L, U

        PLU = self.fetch('PLU_factors')
        if PLU is not None:
            return PLU
        if scipy is None:
            import scipy
        import scipy.linalg
        if numpy is None:
            import numpy
        PM, LM, UM = scipy.linalg.lu(self._matrix_numpy)
        # TODO: It's an awful waste to store a huge matrix for P, which
        # is just a simple permutation, really.
        P = self._new(m, m)
        L = self._new(m, m)
        U = self._new(m, n)
        P._matrix_numpy = numpy.ascontiguousarray(PM)
        L._matrix_numpy = numpy.ascontiguousarray(LM)
        U._matrix_numpy = numpy.ascontiguousarray(UM)
        PLU = (P, L, U)
        for M in PLU:
            M.set_immutable()
        self.cache('PLU_factors', PLU)
        return PLU

    def eigenvalues(self, other=None, algorithm='default', tol=None, *,
                    homogeneous=False):
        r"""
        Return a list of ordinary or generalized eigenvalues.

        INPUT:

        - ``self`` -- a square matrix

        - ``other`` -- a square matrix `B` (default: ``None``) in a generalized
          eigenvalue problem; if ``None``, an ordinary eigenvalue problem is
          solved; if ``algorithm`` is ``'symmetric'`` or ``'hermitian'``, `B`
          must be real symmetric or hermitian positive definite, respectively

        - ``algorithm`` -- (default: ``'default'``)

          - ``'default'`` -- applicable to any matrix
            with double-precision floating point entries.
            Uses the :func:`~scipy:scipy.linalg.eigvals` function from SciPy.

          - ``'symmetric'`` -- converts the matrix into a real matrix
            (i.e. with entries from :class:`~sage.rings.real_double.RDF`),
            then applies the algorithm for Hermitian matrices.  This
            algorithm can be significantly faster than the
            ``'default'`` algorithm.

          - ``'hermitian'`` -- uses the :func:`~scipy:scipy.linalg.eigh`
            function from SciPy, which applies only to real symmetric or
            complex Hermitian matrices.  Since Hermitian is defined as a matrix
            equaling its conjugate-transpose, for a matrix with real
            entries this property is equivalent to being symmetric.
            This algorithm can be significantly faster than the
            ``'default'`` algorithm.

        - ``'tol'`` -- (default: ``None``) if set to a value other than
          ``None``, this is interpreted as a small real number used to aid in
          grouping eigenvalues that are numerically similar, but is ignored
          when ``homogeneous`` is set.  See the output description for more
          information.

        - ``homogeneous`` -- boolean (default: ``False``); if ``True``, use
          homogeneous coordinates for the output
          (see :meth:`eigenvectors_right` for details)

        .. WARNING::

            When using the ``'symmetric'`` or ``'hermitian'`` algorithms,
            no check is made on the input matrix, and only the entries below,
            and on, the main diagonal are employed in the computation.

            Methods such as :meth:`is_symmetric` and :meth:`is_hermitian`
            could be used to verify this beforehand.

        OUTPUT:

        Default output for a square matrix of size `n` is a list of `n`
        eigenvalues from the complex double field,
        :class:`~sage.rings.complex_double.CDF`.  If the ``'symmetric'``
        or ``'hermitian'`` algorithms are chosen, the returned eigenvalues
        are from the real double field,
        :class:`~sage.rings.real_double.RDF`.

        If a tolerance is specified, an attempt is made to group eigenvalues
        that are numerically similar.  The return is then a list of pairs,
        where each pair is an eigenvalue followed by its multiplicity.
        The eigenvalue reported is the mean of the eigenvalues computed,
        and these eigenvalues are contained in an interval (or disk) whose
        radius is less than ``5*tol`` for `n < 10,000` in the worst case.

        More precisely, for an `n\times n` matrix, the diameter of the
        interval containing similar eigenvalues could be as large as sum
        of the reciprocals of the first `n` integers times ``tol``.

        .. WARNING::

            Use caution when using the  ``tol`` parameter to group
            eigenvalues.  See the examples below to see how this can go wrong.

        EXAMPLES::

            sage: m = matrix(RDF, 2, 2, [1,2,3,4])
            sage: ev = m.eigenvalues(); ev
            [-0.372281323..., 5.37228132...]
            sage: ev[0].parent()
            Complex Double Field

            sage: m = matrix(RDF, 2, 2, [0,1,-1,0])
            sage: m.eigenvalues(algorithm='default')
            [1.0*I, -1.0*I]

            sage: m = matrix(CDF, 2, 2, [I,1,-I,0])                                     # needs sage.symbolic
            sage: m.eigenvalues()                                                       # needs sage.symbolic
            [-0.624810533... + 1.30024259...*I, 0.624810533... - 0.30024259...*I]

        The adjacency matrix of a graph will be symmetric, and the
        eigenvalues will be real.  ::

            sage: # needs sage.graphs
            sage: A = graphs.PetersenGraph().adjacency_matrix()
            sage: A = A.change_ring(RDF)
            sage: ev = A.eigenvalues(algorithm='symmetric'); ev  # tol 1e-14
            [-2.0, -2.0, -2.0, -2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 3.0]
            sage: ev[0].parent()
            Real Double Field

        The matrix ``A`` is "random", but the construction of ``C``
        provides a positive-definite Hermitian matrix.  Note that
        the eigenvalues of a Hermitian matrix are real, and the
        eigenvalues of a positive-definite matrix will be positive.  ::

            sage: # needs sage.symbolic
            sage: A = matrix([[ 4*I + 5,  8*I + 1,  7*I + 5, 3*I + 5],
            ....:             [ 7*I - 2, -4*I + 7, -2*I + 4, 8*I + 8],
            ....:             [-2*I + 1,  6*I + 6,  5*I + 5,  -I - 4],
            ....:             [ 5*I + 1,  6*I + 2,    I - 4, -I + 3]])
            sage: C = (A*A.conjugate_transpose()).change_ring(CDF)
            sage: ev = C.eigenvalues(algorithm='hermitian'); ev
            [2.68144025..., 49.5167998..., 274.086188..., 390.71557...]
            sage: ev[0].parent()
            Real Double Field

        A tolerance can be given to aid in grouping eigenvalues that
        are similar numerically.  However, if the parameter is too small
        it might split too finely.  Too large, and it can go wrong very
        badly.  Use with care.  ::

            sage: # needs sage.graphs
            sage: G = graphs.PetersenGraph()
            sage: G.spectrum()
            [3, 1, 1, 1, 1, 1, -2, -2, -2, -2]
            sage: A = G.adjacency_matrix().change_ring(RDF)
            sage: A.eigenvalues(algorithm='symmetric', tol=1.0e-5)  # tol 1e-15
            [(-2.0, 4), (1.0, 5), (3.0, 1)]
            sage: A.eigenvalues(algorithm='symmetric', tol=2.5)  # tol 1e-15
            [(-2.0, 4), (1.3333333333333333, 6)]

        An (extreme) example of properly grouping similar eigenvalues.  ::

            sage: # needs sage.graphs
            sage: G = graphs.HigmanSimsGraph()
            sage: A = G.adjacency_matrix().change_ring(RDF)
            sage: A.eigenvalues(algorithm='symmetric', tol=1.0e-5)  # tol 2e-15
            [(-8.0, 22), (2.0, 77), (22.0, 1)]

        In this generalized eigenvalue problem, the homogeneous coordinates
        explain the output obtained for the eigenvalues::

            sage: A = matrix.identity(RDF, 2)
            sage: B = matrix(RDF, [[3, 5], [6, 10]])
            sage: A.eigenvalues(B)  # tol 1e-14
            [0.0769230769230769, +infinity]
            sage: E = A.eigenvalues(B, homogeneous=True); E  # random
            [(0.9999999999999999, 13.000000000000002), (0.9999999999999999, 0.0)]
            sage: [alpha/beta for alpha, beta in E]  # tol 1e-14
            [0.0769230769230769, NaN + NaN*I]

        .. SEEALSO::

            :meth:`eigenvectors_left`,
            :meth:`eigenvectors_right`,
            :meth:`.Matrix.eigenmatrix_left`,
            :meth:`.Matrix.eigenmatrix_right`.

        TESTS:

        Testing bad input.  ::

            sage: A = matrix(CDF, 2, range(4))
            sage: A.eigenvalues(algorithm='junk')
            Traceback (most recent call last):
            ...
            ValueError: algorithm must be 'default', 'symmetric', or 'hermitian', not junk

            sage: A = matrix(CDF, 2, 3, range(6))
            sage: A.eigenvalues()
            Traceback (most recent call last):
            ...
            ValueError: matrix must be square, not 2 x 3
            sage: matrix.identity(CDF, 2).eigenvalues(A)
            Traceback (most recent call last):
            ...
            ValueError: other matrix must be square, not 2 x 3

            sage: A = matrix(CDF, 2, [1, 2, 3, 4*I])
            sage: A.eigenvalues(algorithm='symmetric')
            Traceback (most recent call last):
            ...
            TypeError: cannot apply symmetric algorithm to matrix with complex entries

            sage: A = matrix(CDF, 2, 2, range(4))
            sage: A.eigenvalues(tol='junk')
            Traceback (most recent call last):
            ...
            TypeError: tolerance parameter must be a real number, not junk

            sage: A = matrix(CDF, 2, 2, range(4))
            sage: A.eigenvalues(tol=-0.01)
            Traceback (most recent call last):
            ...
            ValueError: tolerance parameter must be positive, not -0.01

        A very small matrix.  ::

            sage: matrix(CDF,0,0).eigenvalues()
            []

        Check that homogeneous coordinates work for hermitian positive definite
        input::

            sage: A = matrix.identity(CDF, 2)
            sage: B = matrix(CDF, [[2, 1+I], [1-I, 3]])
            sage: A.eigenvalues(B, algorithm='hermitian', homogeneous=True)  # tol 1e-14
            [(0.25, 1.0), (1.0, 1.0)]

        Test the deprecation::

            sage: # needs sage.graphs
            sage: A = graphs.PetersenGraph().adjacency_matrix().change_ring(RDF)
            sage: ev = A.eigenvalues('symmetric', 1e-13)
            doctest:...: DeprecationWarning: "algorithm" and "tol" should be used
            as keyword argument only
            See https://github.com/sagemath/sage/issues/29243 for details.
            sage: ev  # tol 1e-13
            [(-2.0, 4), (1.0, 5), (3.0, 1)]
            sage: A.eigenvalues('symmetric', 1e-13, tol=1e-12)
            Traceback (most recent call last):
            ...
            TypeError: eigenvalues() got multiple values for keyword argument 'tol'
            sage: A.eigenvalues('symmetric', algorithm='hermitian')
            Traceback (most recent call last):
            ...
            TypeError: eigenvalues() got multiple values for keyword argument 'algorithm'
        """
        from sage.rings.real_double import RDF
        from sage.rings.complex_double import CDF
        if isinstance(other, str):
            # for backward compatibility, allow algorithm to be passed as first
            # positional argument and tol as second positional argument
            from sage.misc.superseded import deprecation
            deprecation(29243, '"algorithm" and "tol" should be used as '
                               'keyword argument only')
            if algorithm != 'default':
                if isinstance(algorithm, str):
                    raise TypeError("eigenvalues() got multiple values for "
                                    "keyword argument 'algorithm'")
                if tol is not None:
                    raise TypeError("eigenvalues() got multiple values for "
                                    "keyword argument 'tol'")
                tol = algorithm
            algorithm = other
            other = None
        if algorithm not in ['default', 'symmetric', 'hermitian']:
            msg = "algorithm must be 'default', 'symmetric', or 'hermitian', not {0}"
            raise ValueError(msg.format(algorithm))
        if not self.is_square():
            raise ValueError('matrix must be square, not %s x %s'
                             % (self.nrows(), self.ncols()))
        if other is not None and not other.is_square():
            raise ValueError('other matrix must be square, not %s x %s'
                             % (other.nrows(), other.ncols()))
        if algorithm == 'symmetric':
            if self.base_ring() != RDF:
                try:
                    self = self.change_ring(RDF)  # check side effect
                except TypeError:
                    raise TypeError('cannot apply symmetric algorithm to matrix with complex entries')
            if other is not None and other.base_ring() != RDF:
                try:
                    other = other.change_ring(RDF)  # check side effect
                except TypeError:
                    raise TypeError('cannot apply symmetric algorithm to matrix with complex entries')
            algorithm = 'hermitian'
        if homogeneous:
            tol = None
        multiplicity = (tol is not None)
        if multiplicity:
            try:
                tol = float(tol)
            except (ValueError, TypeError):
                msg = 'tolerance parameter must be a real number, not {0}'
                raise TypeError(msg.format(tol))
            if tol < 0:
                msg = 'tolerance parameter must be positive, not {0}'
                raise ValueError(msg.format(tol))

        if self._nrows == 0:
            return []
        global scipy
        if scipy is None:
            import scipy
        import scipy.linalg
        global numpy
        if numpy is None:
            import numpy
        other_numpy = None if other is None else other.numpy()
        # generic eigenvalues, or real eigenvalues for Hermitian
        if algorithm == 'default':
            return_class = CDF
            evalues = scipy.linalg.eigvals(self._matrix_numpy, other_numpy,
                                           homogeneous_eigvals=homogeneous)
        elif algorithm == 'hermitian':
            return_class = RDF
            evalues = scipy.linalg.eigh(self._matrix_numpy, other_numpy,
                                        eigvals_only=True)
            if homogeneous:
                # eigh does not support homogeneous output
                evalues = evalues, [RDF.one()] * len(evalues)

        if homogeneous:
            return [(return_class(a), return_class(b))
                    for a, b in zip(*evalues)]
        elif not multiplicity:
            return [return_class(e) for e in evalues]
        else:
            # pairs in ev_group are
            #   slot 0: the sum of "equal" eigenvalues, "s"
            #   slot 1: number of eigenvalues in this sum, "m"
            #   slot 2: average of these eigenvalues, "avg"
            # we test if "new" eigenvalues are close to the group average
            ev_group = []
            for e in evalues:
                location = None
                best_fit = tol
                for i in range(len(ev_group)):
                    _, m, avg = ev_group[i]
                    d = numpy.abs(avg - e)
                    if d < best_fit:
                        best_fit = d
                        location = i
                if location is None:
                    ev_group.append([e, 1, e])
                else:
                    ev_group[location][0] += e
                    ev_group[location][1] += 1
                    ev_group[location][2] = ev_group[location][0]/ev_group[location][1]
            return [(return_class(avg), m) for _, m, avg in ev_group]

    def left_eigenvectors(self, other=None, *, homogeneous=False):
        r"""
        Compute the ordinary or generalized left eigenvectors of a matrix of
        double precision real or complex numbers (i.e. ``RDF`` or ``CDF``).

        INPUT:

        - ``other`` -- a square matrix `B` (default: ``None``) in a generalized
          eigenvalue problem; if ``None``, an ordinary eigenvalue problem is
          solved

        - ``homogeneous`` -- boolean (default: ``False``); if ``True``, use
          homogeneous coordinates for the eigenvalues in the output

        OUTPUT:

        A list of triples, each of the form ``(e,[v],1)``,
        where ``e`` is the eigenvalue, and ``v`` is an associated
        left eigenvector such that

        .. MATH::

            v A = e v.

        If the matrix `A` is of size `n`, then there are `n` triples.

        If a matrix `B` is passed as optional argument, the output is a
        solution to the generalized eigenvalue problem such that

        .. MATH::

            v A = e v B.

        If ``homogeneous`` is set, each eigenvalue is returned as a tuple
        `(\alpha, \beta)` of homogeneous coordinates such that

        .. MATH::

            \beta v A = \alpha v B.

        The format of the output is designed to match the format
        for exact results.  However, since matrices here have numerical
        entries, the resulting eigenvalues will also be numerical.  No
        attempt is made to determine if two eigenvalues are equal, or if
        eigenvalues might actually be zero.  So the algebraic multiplicity
        of each eigenvalue is reported as 1.  Decisions about equal
        eigenvalues or zero eigenvalues should be addressed in the
        calling routine.

        The SciPy routines used for these computations produce eigenvectors
        normalized to have length 1, but on different hardware they may vary
        by a complex sign. So for doctests we have normalized output by forcing
        their eigenvectors to have their first nonzero entry equal to one.

        ALGORITHM:

        Values are computed with the SciPy library using
        :func:`scipy:scipy.linalg.eig`.

        EXAMPLES::

            sage: m = matrix(RDF, [[-5, 3, 2, 8],[10, 2, 4, -2],[-1, -10, -10, -17],[-2, 7, 6, 13]])
            sage: m
            [ -5.0   3.0   2.0   8.0]
            [ 10.0   2.0   4.0  -2.0]
            [ -1.0 -10.0 -10.0 -17.0]
            [ -2.0   7.0   6.0  13.0]
            sage: spectrum = m.left_eigenvectors()
            sage: for i in range(len(spectrum)):
            ....:     spectrum[i][1][0] = matrix(RDF, spectrum[i][1]).echelon_form()[0]
            sage: spectrum[0]  # tol 1e-13
            (2.0, [(1.0, 1.0, 1.0, 1.0)], 1)
            sage: spectrum[1]  # tol 1e-13
            (1.0, [(1.0, 0.8, 0.8, 0.6)], 1)
            sage: spectrum[2]  # tol 1e-13
            (-2.0, [(1.0, 0.4, 0.6, 0.2)], 1)
            sage: spectrum[3]  # tol 1e-13
            (-1.0, [(1.0, 1.0, 2.0, 2.0)], 1)

        A generalized eigenvalue problem::

            sage: A = matrix(CDF, [[1+I, -2], [3, 4]])
            sage: B = matrix(CDF, [[0, 7-I], [2, -3]])
            sage: E = A.eigenvectors_left(B)
            sage: all((v * A - e * v * B).norm() < 1e-14 for e, [v], _ in E)
            True

        In a generalized eigenvalue problem with a singular matrix `B`, we can
        check the eigenvector property using homogeneous coordinates, even
        though the quotient `\alpha/\beta` is not always defined::

            sage: A = matrix.identity(CDF, 2)
            sage: B = matrix(CDF, [[2, 1+I], [4, 2+2*I]])
            sage: E = A.eigenvectors_left(B, homogeneous=True)
            sage: all((beta * v * A - alpha * v * B).norm() < 1e-14
            ....:     for (alpha, beta), [v], _ in E)
            True

        .. SEEALSO::

            :meth:`eigenvalues`,
            :meth:`eigenvectors_right`,
            :meth:`.Matrix.eigenmatrix_left`.

        TESTS:

        The following example shows that :issue:`20439` has been resolved::

            sage: A = matrix(CDF, [[-2.53634347567,  2.04801738686, -0.0, -62.166145304],
            ....:                  [ 0.7, -0.6, 0.0, 0.0],
            ....:                  [0.547271128842, 0.0, -0.3015, -21.7532081652],
            ....:                  [0.0, 0.0, 0.3, -0.4]])
            sage: spectrum = A.left_eigenvectors()
            sage: all((Matrix(spectrum[i][1])*(A - spectrum[i][0])).norm() < 10^(-2)
            ....:     for i in range(A.nrows()))
            True

        The following example shows that the fix for :issue:`20439` (conjugating
        eigenvectors rather than eigenvalues) is the correct one::

            sage: A = Matrix(CDF,[[I,0],[0,1]])
            sage: spectrum = A.left_eigenvectors()
            sage: for i in range(len(spectrum)):
            ....:   spectrum[i][1][0] = matrix(CDF, spectrum[i][1]).echelon_form()[0]
            sage: spectrum
            [(1.0*I, [(1.0, 0.0)], 1), (1.0, [(0.0, 1.0)], 1)]
        """
        if not self.is_square():
            raise ArithmeticError("self must be a square matrix")
        if other is not None and not other.is_square():
            raise ArithmeticError("other must be a square matrix")
        if self._nrows == 0:
            return [], self.__copy__()
        global scipy
        if scipy is None:
            import scipy
        import scipy.linalg
        v, eig = scipy.linalg.eig(self._matrix_numpy,
                                  None if other is None else other.numpy(),
                                  right=False, left=True,
                                  homogeneous_eigvals=homogeneous)
        # scipy puts eigenvectors in columns, we will extract from rows
        eig = matrix(eig.T)
        if other is not None:
            # scipy fails to normalize generalized left eigenvectors
            # (see https://github.com/scipy/scipy/issues/11550),
            # FIXME: remove this normalization step once that issue is resolved
            eig = [v.normalized() for v in eig]
        from sage.rings.complex_double import CDF
        if homogeneous:
            v = [(CDF(a), CDF(b)) for a, b in v.T]
        else:
            v = [CDF(e) for e in v]
        return [(v[i], [eig[i].conjugate()], 1) for i in range(len(v))]

    eigenvectors_left = left_eigenvectors

    def right_eigenvectors(self, other=None, *, homogeneous=False):
        r"""
        Compute the ordinary or generalized right eigenvectors of a matrix of
        double precision real or complex numbers (i.e. ``RDF`` or ``CDF``).

        INPUT:

        - ``other`` -- a square matrix `B` (default: ``None``) in a generalized
          eigenvalue problem; if ``None``, an ordinary eigenvalue problem is
          solved

        - ``homogeneous`` -- boolean (default: ``False``); if ``True``, use
          homogeneous coordinates for the eigenvalues in the output

        OUTPUT:

        A list of triples, each of the form ``(e,[v],1)``,
        where ``e`` is the eigenvalue, and ``v`` is an associated
        right eigenvector such that

        .. MATH::

            A v = e v.

        If the matrix `A` is of size `n`, then there are `n` triples.

        If a matrix `B` is passed as optional argument, the output is a
        solution to the generalized eigenvalue problem such that

        .. MATH::

            A v = e B v.

        If ``homogeneous`` is set, each eigenvalue is returned as a tuple
        `(\alpha, \beta)` of homogeneous coordinates such that

        .. MATH::

            \beta A v = \alpha B v.

        The format of the output is designed to match the format
        for exact results.  However, since matrices here have numerical
        entries, the resulting eigenvalues will also be numerical.  No
        attempt is made to determine if two eigenvalues are equal, or if
        eigenvalues might actually be zero.  So the algebraic multiplicity
        of each eigenvalue is reported as 1.  Decisions about equal
        eigenvalues or zero eigenvalues should be addressed in the
        calling routine.

        The SciPy routines used for these computations produce eigenvectors
        normalized to have length 1, but on different hardware they may vary
        by a complex sign. So for doctests we have normalized output by forcing
        their eigenvectors to have their first nonzero entry equal to one.

        ALGORITHM:

        Values are computed with the SciPy library using
        :func:`scipy:scipy.linalg.eig`.

        EXAMPLES::

            sage: m = matrix(RDF, [[-9, -14, 19, -74],[-1, 2, 4, -11],[-4, -12, 6, -32],[0, -2, -1, 1]])
            sage: m
            [ -9.0 -14.0  19.0 -74.0]
            [ -1.0   2.0   4.0 -11.0]
            [ -4.0 -12.0   6.0 -32.0]
            [  0.0  -2.0  -1.0   1.0]
            sage: spectrum = m.right_eigenvectors()
            sage: for i in range(len(spectrum)):
            ....:   spectrum[i][1][0] = matrix(RDF, spectrum[i][1]).echelon_form()[0]
            sage: spectrum[0]  # tol 1e-13
            (2.0, [(1.0, -2.0, 3.0, 1.0)], 1)
            sage: spectrum[1]  # tol 1e-13
            (1.0, [(1.0, -0.666666666666633, 1.333333333333286, 0.33333333333331555)], 1)
            sage: spectrum[2]  # tol 1e-13
            (-2.0, [(1.0, -0.2, 1.0, 0.2)], 1)
            sage: spectrum[3]  # tol 1e-12
            (-1.0, [(1.0, -0.5, 2.0, 0.5)], 1)

        A generalized eigenvalue problem::

            sage: A = matrix(CDF, [[1+I, -2], [3, 4]])
            sage: B = matrix(CDF, [[0, 7-I], [2, -3]])
            sage: E = A.eigenvectors_right(B)
            sage: all((A * v - e * B * v).norm() < 1e-14 for e, [v], _ in E)
            True

        In a generalized eigenvalue problem with a singular matrix `B`, we can
        check the eigenvector property using homogeneous coordinates, even
        though the quotient `\alpha/\beta` is not always defined::

            sage: A = matrix.identity(RDF, 2)
            sage: B = matrix(RDF, [[3, 5], [6, 10]])
            sage: E = A.eigenvectors_right(B, homogeneous=True)
            sage: all((beta * A * v - alpha * B * v).norm() < 1e-14
            ....:     for (alpha, beta), [v], _ in E)
            True

        .. SEEALSO::

            :meth:`eigenvalues`,
            :meth:`eigenvectors_left`,
            :meth:`.Matrix.eigenmatrix_right`.

        TESTS:

        The following example shows that :issue:`20439` has been resolved::

            sage: A = matrix(CDF, [[-2.53634347567,  2.04801738686, -0.0, -62.166145304],
            ....:                  [ 0.7, -0.6, 0.0, 0.0],
            ....:                  [0.547271128842, 0.0, -0.3015, -21.7532081652],
            ....:                  [0.0, 0.0, 0.3, -0.4]])
            sage: spectrum = A.right_eigenvectors()
            sage: all(((A - spectrum[i][0]) * Matrix(spectrum[i][1]).transpose()).norm() < 10^(-2)
            ....:     for i in range(A.nrows()))
            True

        The following example shows that the fix for :issue:`20439` (conjugating
        eigenvectors rather than eigenvalues) is the correct one::

            sage: A = Matrix(CDF,[[I,0],[0,1]])
            sage: spectrum = A.right_eigenvectors()
            sage: for i in range(len(spectrum)):
            ....:     spectrum[i][1][0] = matrix(CDF, spectrum[i][1]).echelon_form()[0]
            sage: spectrum
            [(1.0*I, [(1.0, 0.0)], 1), (1.0, [(0.0, 1.0)], 1)]
        """
        if not self.is_square():
            raise ArithmeticError("self must be a square matrix")
        if other is not None and not other.is_square():
            raise ArithmeticError("other must be a square matrix")
        if self._nrows == 0:
            return [], self.__copy__()
        global scipy
        if scipy is None:
            import scipy
        import scipy.linalg
        v, eig = scipy.linalg.eig(self._matrix_numpy,
                                  None if other is None else other.numpy(),
                                  right=True, left=False,
                                  homogeneous_eigvals=homogeneous)
        # scipy puts eigenvectors in columns, we will extract from rows
        eig = matrix(eig.T)
        from sage.rings.complex_double import CDF
        if homogeneous:
            v = [(CDF(a), CDF(b)) for a, b in v.T]
        else:
            v = [CDF(e) for e in v]
        return [(v[i], [eig[i]], 1) for i in range(len(v))]

    eigenvectors_right = right_eigenvectors

    def _solve_right_nonsingular_square(self, B, check_rank=False):
        """
        Find a solution `X` to the equation `A X = B` if ``self`` is a square
        matrix `A`.

        ALGORITHM:

        Uses the function :func:`scipy:scipy.linalg.solve` from SciPy.

        TESTS::

            sage: # needs sage.symbolic
            sage: A = matrix(CDF, [[1, 2], [3, 3+I]])
            sage: b = matrix(CDF, [[1, 0], [2, 1]])
            sage: x = A._solve_right_nonsingular_square(b)
            sage: (A * x - b).norm() < 1e-14
            True
        """
        global scipy
        if scipy is None:
            import scipy
        import scipy.linalg
        X = self._new(self._ncols, B.ncols())
        # may raise a LinAlgError for a singular matrix
        X._matrix_numpy = scipy.linalg.solve(self._matrix_numpy, B.numpy())
        return X

    def _solve_right_general(self, B, check=False):
        """
        Compute a least-squares solution `X` to the equation `A X = B` where
        ``self`` is the matrix `A`.

        ALGORITHM:

        Uses the function :func:`scipy:scipy.linalg.lstsq` from SciPy.

        TESTS::

            sage: A = matrix(RDF, 3, 2, [1, 3, 4, 2, 0, -3])
            sage: b = matrix(RDF, 3, 2, [5, 6, 1, 0, 0, 2])
            sage: x = A._solve_right_general(b)
            sage: y = ~(A.T * A) * A.T * b  # closed form solution
            sage: (x - y).norm() < 1e-14
            True
        """
        global scipy
        if scipy is None:
            import scipy
        import scipy.linalg
        X = self._new(self._ncols, B.ncols())
        arr = scipy.linalg.lstsq(self._matrix_numpy, B.numpy())[0]
        X._matrix_numpy = arr
        return X

    def determinant(self):
        """
        Return the determinant of ``self``.

        ALGORITHM:

        Uses :func:`scipy:scipy.linalg.det`.

        EXAMPLES::

            sage: m = matrix(RDF,2,range(4)); m.det()
            -2.0
            sage: m = matrix(RDF,0,[]); m.det()
            1.0
            sage: m = matrix(RDF, 2, range(6)); m.det()
            Traceback (most recent call last):
            ...
            ValueError: self must be a square matrix
        """
        if not self.is_square():
            raise ValueError("self must be a square matrix")
        if self._nrows == 0 or self._ncols == 0:
            return self._sage_dtype(1)
        global scipy
        if scipy is None:
            import scipy
        import scipy.linalg

        return self._sage_dtype(scipy.linalg.det(self._matrix_numpy))

    def log_determinant(self):
        """
        Compute the log of the absolute value of the determinant
        using LU decomposition.

        .. NOTE::

            This is useful if the usual determinant overflows.

        EXAMPLES::

            sage: m = matrix(RDF,2,2,range(4)); m
            [0.0 1.0]
            [2.0 3.0]
            sage: RDF(log(abs(m.determinant())))
            0.6931471805599453
            sage: m.log_determinant()
            0.6931471805599453
            sage: m = matrix(RDF,0,0,[]); m
            []
            sage: m.log_determinant()
            0.0
            sage: m = matrix(CDF,2,2,range(4)); m
            [0.0 1.0]
            [2.0 3.0]
            sage: RDF(log(abs(m.determinant())))
            0.6931471805599453
            sage: m.log_determinant()
            0.6931471805599453
            sage: m = matrix(CDF,0,0,[]); m
            []
            sage: m.log_determinant()
            0.0
        """
        global numpy
        cdef Matrix_double_dense U

        if self._nrows == 0 or self._ncols == 0:
            return sage.rings.real_double.RDF(0)

        if not self.is_square():
            raise ArithmeticError("self must be a square matrix")

        _, _, U = self.LU()
        if numpy is None:
            import numpy

        return sage.rings.real_double.RDF(sum(numpy.log(abs(numpy.diag(U._matrix_numpy)))))

    def conjugate(self):
        r"""
        Return the conjugate of this matrix, i.e. the matrix whose entries are
        the conjugates of the entries of ``self``.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: A = matrix(CDF, [[1+I, 3-I], [0, 2*I]])
            sage: A.conjugate()
            [1.0 - 1.0*I 3.0 + 1.0*I]
            [        0.0      -2.0*I]

        There is a shorthand notation::

            sage: A.conjugate() == A.C                                                  # needs sage.symbolic
            True

        Conjugates work (trivially) for real matrices::

            sage: B = matrix.random(RDF, 3)
            sage: B == B.conjugate()
            True

        TESTS::

            sage: matrix(CDF, 0).conjugate()
            []
        """
        cdef Matrix_double_dense A
        A = self._new(self._nrows, self._ncols)
        A._matrix_numpy = self._matrix_numpy.conjugate()
        if self._subdivisions is not None:
            A.subdivide(*self.subdivisions())
        return A

    def SVD(self):
        r"""
        Return the singular value decomposition of this matrix.

        The `U` and `V` matrices are not unique and may be returned with different
        values in the future or on different systems. The `S` matrix is unique
        and contains the singular values in descending order.

        The computed decomposition is cached and returned on subsequent calls.

        INPUT:

        - ``A`` -- a matrix

        OUTPUT:

        ``U, S, V`` -- immutable matrices such that ``A = U*S*V.conj().transpose()``
        where `U` and `V` are orthogonal and `S` is zero off of the diagonal

        Note that if ``self`` is m-by-n, then the dimensions of the
        matrices that this returns are (m,m), (m,n), and (n, n).

        .. NOTE::

            If all you need is the singular values of the matrix, see
            the more convenient :meth:`singular_values`.

        EXAMPLES::

            sage: m = matrix(RDF,4,range(1,17))
            sage: U,S,V = m.SVD()
            sage: U*S*V.transpose()  # tol 1e-14
            [0.9999999999999993 1.9999999999999987  3.000000000000001  4.000000000000002]
            [ 4.999999999999998  5.999999999999998  6.999999999999998                8.0]
            [ 8.999999999999998  9.999999999999996 10.999999999999998               12.0]
            [12.999999999999998               14.0               15.0               16.0]

        A non-square example::

            sage: m = matrix(RDF, 2, range(1,7)); m
            [1.0 2.0 3.0]
            [4.0 5.0 6.0]
            sage: U, S, V = m.SVD()
            sage: U*S*V.transpose()  # tol 1e-14
            [0.9999999999999994 1.9999999999999998  2.999999999999999]
            [ 4.000000000000001  5.000000000000002  6.000000000000001]

        S contains the singular values::

            sage: S.round(4)
            [ 9.508    0.0    0.0]
            [   0.0 0.7729    0.0]
            sage: [N(sqrt(abs(x)), digits=4) for x in (S*S.transpose()).eigenvalues()]
            [9.508, 0.7729]

        U and V are orthogonal matrices::

            sage: U # random, SVD is not unique
            [-0.386317703119 -0.922365780077]
            [-0.922365780077  0.386317703119]
            [-0.274721127897 -0.961523947641]
            [-0.961523947641  0.274721127897]
            sage: (U*U.transpose())  # tol 1e-15
            [               1.0                0.0]
            [               0.0 1.0000000000000004]
            sage: V # random, SVD is not unique
            [-0.428667133549  0.805963908589  0.408248290464]
            [-0.566306918848  0.112382414097 -0.816496580928]
            [-0.703946704147 -0.581199080396  0.408248290464]
            sage: (V*V.transpose())  # tol 1e-15
            [0.9999999999999999                0.0                0.0]
            [               0.0                1.0                0.0]
            [               0.0                0.0 0.9999999999999999]

        TESTS::

            sage: m = matrix(RDF,3,2,range(1, 7)); m
            [1.0 2.0]
            [3.0 4.0]
            [5.0 6.0]
            sage: U,S,V = m.SVD()
            sage: U*S*V.transpose()  # tol 1e-15
            [0.9999999999999996 1.9999999999999998]
            [               3.0 3.9999999999999996]
            [ 4.999999999999999  6.000000000000001]

            sage: m = matrix(RDF, 3, 0, []); m
            []
            sage: m.SVD()
            ([], [], [])
            sage: m = matrix(RDF, 0, 3, []); m
            []
            sage: m.SVD()
            ([], [], [])
            sage: def shape(x): return (x.nrows(), x.ncols())
            sage: m = matrix(RDF, 2, 3, range(6))
            sage: list(map(shape, m.SVD()))
            [(2, 2), (2, 3), (3, 3)]
            sage: for x in m.SVD(): x.is_immutable()
            True
            True
            True
        """
        global scipy, numpy
        cdef Py_ssize_t i
        cdef Matrix_double_dense U, S, V

        if self._nrows == 0 or self._ncols == 0:
            U_t = self.new_matrix(self._nrows, self._ncols)
            S_t = self.new_matrix(self._nrows, self._ncols)
            V_t = self.new_matrix(self._ncols, self._nrows)
            return U_t, S_t, V_t

        USV = self.fetch('SVD_factors')
        if USV is None:
            # TODO: More efficient representation of non-square diagonal matrix S
            if scipy is None:
                import scipy
            import scipy.linalg
            if numpy is None:
                import numpy
            U_mat, S_diagonal, V_mat = scipy.linalg.svd(self._matrix_numpy)

            U = self._new(self._nrows, self._nrows)
            S = self._new(self._nrows, self._ncols)
            V = self._new(self._ncols, self._ncols)

            S_mat = numpy.zeros((self._nrows, self._ncols), dtype=self._numpy_dtype)
            for i in range(S_diagonal.shape[0]):
                S_mat[i,i] = S_diagonal[i]

            U._matrix_numpy = numpy.ascontiguousarray(U_mat)
            S._matrix_numpy = S_mat
            V._matrix_numpy = numpy.ascontiguousarray(V_mat.conj().T)
            USV = U, S, V
            for M in USV: M.set_immutable()
            self.cache('SVD_factors', USV)

        return USV

    def QR(self):
        r"""
        Return a factorization into a unitary matrix and an
        upper-triangular matrix.

        Applies to any matrix over ``RDF`` or ``CDF``.

        OUTPUT:

        ``Q``, ``R`` -- a pair of matrices such that if `A`
        is the original matrix, then

        .. MATH::

            A = QR, \quad Q^\ast Q = I

        where `R` is upper-triangular.  `Q^\ast` is the
        conjugate-transpose in the complex case, and just
        the transpose in the real case. So `Q` is a unitary
        matrix (or rather, orthogonal, in the real case),
        or equivalently `Q` has orthogonal columns.  For a
        matrix of full rank this factorization is unique
        up to adjustments via multiples of rows and columns
        by multiples with scalars having modulus `1`.  So
        in the full-rank case, `R` is unique if the diagonal
        entries are required to be positive real numbers.

        The resulting decomposition is cached.

        ALGORITHM:

        Calls :func:`scipy:scipy.linalg.qr` from SciPy, which is in turn an
        interface to LAPACK routines.

        EXAMPLES:

        Over the reals, the inverse of ``Q`` is its transpose,
        since including a conjugate has no effect.  In the real
        case, we say ``Q`` is orthogonal. ::

            sage: A = matrix(RDF, [[-2, 0, -4, -1, -1],
            ....:                  [-2, 1, -6, -3, -1],
            ....:                  [1, 1, 7, 4, 5],
            ....:                  [3, 0, 8, 3, 3],
            ....:                  [-1, 1, -6, -6, 5]])
            sage: Q, R = A.QR()

        At this point, ``Q`` is only well-defined up to the
        signs of its columns, and similarly for ``R`` and its
        rows, so we normalize them::

            sage: Qnorm = Q._normalize_columns()
            sage: Rnorm = R._normalize_rows()
            sage: Qnorm.round(6).zero_at(10^-6)
            [ 0.458831  0.126051  0.381212  0.394574   0.68744]
            [ 0.458831  -0.47269 -0.051983 -0.717294  0.220963]
            [-0.229416 -0.661766  0.661923  0.180872 -0.196411]
            [-0.688247 -0.189076 -0.204468  -0.09663  0.662889]
            [ 0.229416 -0.535715 -0.609939  0.536422 -0.024551]
            sage: Rnorm.round(6).zero_at(10^-6)
            [ 4.358899 -0.458831 13.076697  6.194225  2.982405]
            [      0.0  1.670172  0.598741  -1.29202  6.207997]
            [      0.0       0.0  5.444402  5.468661 -0.682716]
            [      0.0       0.0       0.0  1.027626   -3.6193]
            [      0.0       0.0       0.0       0.0  0.024551]
            sage: (Q*Q.transpose())  # tol 1e-14
            [0.9999999999999994                0.0                0.0                0.0                0.0]
            [               0.0                1.0                0.0                0.0                0.0]
            [               0.0                0.0 0.9999999999999999                0.0                0.0]
            [               0.0                0.0                0.0 0.9999999999999998                0.0]
            [               0.0                0.0                0.0                0.0 1.0000000000000002]
            sage: (Q*R - A).zero_at(10^-14)
            [0.0 0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0 0.0]

        Now over the complex numbers, demonstrating that the SciPy libraries
        are (properly) using the Hermitian inner product, so that ``Q`` is
        a unitary matrix (its inverse is the conjugate-transpose).  ::

            sage: A = matrix(CDF, [[-8, 4*I + 1, -I + 2, 2*I + 1],
            ....:                  [1, -2*I - 1, -I + 3, -I + 1],
            ....:                  [I + 7, 2*I + 1, -2*I + 7, -I + 1],
            ....:                  [I + 2, 0, I + 12, -1]])
            sage: Q, R = A.QR()
            sage: Q._normalize_columns()  # tol 1e-6
            [                           0.7302967433402214    0.20705664550556482 + 0.5383472783144685*I   0.24630498099986423 - 0.07644563587232917*I   0.23816176831943323 - 0.10365960327796941*I]
            [                         -0.09128709291752768  -0.20705664550556482 - 0.37787837804765584*I   0.37865595338630315 - 0.19522214955246678*I    0.7012444502144682 - 0.36437116509865947*I]
            [  -0.6390096504226938 - 0.09128709291752768*I    0.17082173254209104 + 0.6677576817554466*I -0.03411475806452064 + 0.040901987417671426*I   0.31401710855067644 - 0.08251917187054114*I]
            [ -0.18257418583505536 - 0.09128709291752768*I  -0.03623491296347384 + 0.07246982592694771*I    0.8632284069415112 + 0.06322839976356195*I  -0.44996948676115206 - 0.01161191812089182*I]
            sage: R._normalize_rows().zero_at(1e-15)  # tol 1e-6
            [                        10.954451150103322                      -1.9170289512680814*I   5.385938482134133 - 2.1908902300206643*I -0.2738612787525829 - 2.1908902300206643*I]
            [                                       0.0                            4.8295962564173  -0.8696379111233719 - 5.864879483945123*I  0.993871898426711 - 0.30540855212070794*I]
            [                                       0.0                                        0.0                          12.00160760935814 -0.2709533402297273 + 0.4420629644486325*I]
            [                                       0.0                                        0.0                                        0.0                         1.9429639442589917]
            sage: (Q.conjugate().transpose()*Q).zero_at(1e-15)  # tol 1e-15
            [               1.0                0.0                0.0                0.0]
            [               0.0 0.9999999999999994                0.0                0.0]
            [               0.0                0.0 1.0000000000000002                0.0]
            [               0.0                0.0                0.0 1.0000000000000004]
            sage: (Q*R - A).zero_at(10^-14)
            [0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0]

        An example of a rectangular matrix that is also rank-deficient.
        If you run this example yourself, you may see a very small, nonzero
        entries in the third row, in the third column, even though the exact
        version of the matrix has rank 2.  The final two columns of ``Q``
        span the left kernel of ``A`` (as evidenced by the two zero rows of
        ``R``).  Different platforms will compute different bases for this
        left kernel, so we do not exhibit the actual matrix.  ::

            sage: Arat = matrix(QQ, [[2, -3, 3],
            ....:                    [-1, 1, -1],
            ....:                    [-1, 3, -3],
            ....:                    [-5, 1, -1]])
            sage: Arat.rank()
            2
            sage: A = Arat.change_ring(CDF)
            sage: Q, R = A.QR()
            sage: R._normalize_rows()  # abs tol 1e-14
            [     5.567764362830022    -2.6940795304016243     2.6940795304016243]
            [                   0.0     3.5695847775155825    -3.5695847775155825]
            [                   0.0                    0.0 2.4444034681064287e-16]
            [                   0.0                    0.0                    0.0]
            sage: (Q.conjugate_transpose()*Q)  # abs tol 1e-14
            [     1.0000000000000002  -5.185196889911925e-17 -4.1457180570414476e-17  -2.909388767229071e-17]
            [ -5.185196889911925e-17      1.0000000000000002  -9.286869233696149e-17 -1.1035822863186828e-16]
            [-4.1457180570414476e-17  -9.286869233696149e-17                     1.0  4.4159215672155694e-17]
            [ -2.909388767229071e-17 -1.1035822863186828e-16  4.4159215672155694e-17                     1.0]

        Results are cached, meaning they are immutable matrices.
        Make a copy if you need to manipulate a result. ::

            sage: A = random_matrix(CDF, 2, 2)
            sage: Q, R = A.QR()
            sage: Q.is_mutable()
            False
            sage: R.is_mutable()
            False
            sage: Q[0,0] = 0
            Traceback (most recent call last):
            ...
            ValueError: matrix is immutable; please change a copy instead (i.e., use copy(M) to change a copy of M).
            sage: Qcopy = copy(Q)
            sage: Qcopy[0,0] = 679
            sage: Qcopy[0,0]
            679.0

        TESTS:

        Trivial cases return trivial results of the correct size,
        and we check ``Q`` itself in one case, verifying a fix for
        :issue:`10795`.  ::

            sage: A = zero_matrix(RDF, 0, 10)
            sage: Q, R = A.QR()
            sage: Q.nrows(), Q.ncols()
            (0, 0)
            sage: R.nrows(), R.ncols()
            (0, 10)
            sage: A = zero_matrix(RDF, 3, 0)
            sage: Q, R = A.QR()
            sage: Q.nrows(), Q.ncols()
            (3, 3)
            sage: R.nrows(), R.ncols()
            (3, 0)
            sage: Q
            [1.0 0.0 0.0]
            [0.0 1.0 0.0]
            [0.0 0.0 1.0]
        """
        global scipy
        cdef Matrix_double_dense Q,R

        if self._nrows == 0 or self._ncols == 0:
            return self.new_matrix(self._nrows, self._nrows, entries=1), self.new_matrix()

        QR = self.fetch('QR_factors')
        if QR is None:
            Q = self._new(self._nrows, self._nrows)
            R = self._new(self._nrows, self._ncols)
            if scipy is None:
                import scipy
            import scipy.linalg
            Q._matrix_numpy, R._matrix_numpy = scipy.linalg.qr(self._matrix_numpy)
            Q.set_immutable()
            R.set_immutable()
            QR = (Q, R)
            self.cache('QR_factors', QR)
        return QR

    def is_unitary(self, tol=1e-12, algorithm='orthonormal'):
        r"""
        Return ``True`` if the columns of the matrix are an orthonormal basis.

        For a matrix with real entries this determines if a matrix is
        "orthogonal" and for a matrix with complex entries this determines
        if the matrix is "unitary."

        INPUT:

        - ``tol`` -- (default: ``1e-12``) the largest value of the
          absolute value of the difference between two matrix entries
          for which they will still be considered equal

        - ``algorithm`` -- (default: ``'orthonormal'``) set to
          ``'orthonormal'`` for a stable procedure and set to 'naive' for a
          fast procedure

        OUTPUT:

        ``True`` if the matrix is square and its conjugate-transpose is
        its inverse, and ``False`` otherwise.  In other words, a matrix
        is orthogonal or unitary if the product of its conjugate-transpose
        times the matrix is the identity matrix.

        The tolerance parameter is used to allow for numerical values
        to be equal if there is a slight difference due to round-off
        and other imprecisions.

        The result is cached, on a per-tolerance and per-algorithm basis.

        ALGORITHMS:

        The naive algorithm simply computes the product of the
        conjugate-transpose with the matrix and compares the entries
        to the identity matrix, with equality controlled by the
        tolerance parameter.

        The orthonormal algorithm first computes a Schur decomposition
        (via the :meth:`schur` method) and checks that the result is a
        diagonal matrix with entries of modulus 1, which is equivalent to
        being unitary.

        So the naive algorithm might finish fairly quickly for a matrix
        that is not unitary, once the product has been computed.
        However, the orthonormal algorithm will compute a Schur
        decomposition before going through a similar check of a
        matrix entry-by-entry.

        EXAMPLES:

        A matrix that is far from unitary. ::

            sage: A = matrix(RDF, 4, range(16))
            sage: A.conjugate().transpose()*A
            [224.0 248.0 272.0 296.0]
            [248.0 276.0 304.0 332.0]
            [272.0 304.0 336.0 368.0]
            [296.0 332.0 368.0 404.0]
            sage: A.is_unitary()
            False
            sage: A.is_unitary(algorithm='naive')
            False
            sage: A.is_unitary(algorithm='orthonormal')
            False

        The QR decomposition will produce a unitary matrix as Q and the
        SVD decomposition will create two unitary matrices, U and V. ::

            sage: # needs sage.symbolic
            sage: A = matrix(CDF, [[   1 - I,   -3*I,  -2 + I,        1, -2 + 3*I],
            ....:                  [   1 - I, -2 + I, 1 + 4*I,        0,    2 + I],
            ....:                  [      -1, -5 + I,  -2 + I,    1 + I, -5 - 4*I],
            ....:                  [-2 + 4*I,  2 - I, 8 - 4*I,  1 - 8*I,  3 - 2*I]])
            sage: Q, R = A.QR()
            sage: Q.is_unitary()
            True
            sage: U, S, V = A.SVD()
            sage: U.is_unitary(algorithm='naive')
            True
            sage: U.is_unitary(algorithm='orthonormal')
            True
            sage: V.is_unitary(algorithm='naive')
            True

        If we make the tolerance too strict we can get misleading results.  ::

            sage: A = matrix(RDF, 10, 10, [1/(i+j+1) for i in range(10) for j in range(10)])
            sage: Q, R = A.QR()
            sage: Q.is_unitary(algorithm='naive', tol=1e-16)
            False
            sage: Q.is_unitary(algorithm='orthonormal', tol=1e-17)
            False

        Rectangular matrices are not unitary/orthogonal, even if their
        columns form an orthonormal set.  ::

            sage: A = matrix(CDF, [[1,0], [0,0], [0,1]])
            sage: A.is_unitary()
            False

        The smallest cases::

            sage: P = matrix(CDF, 0, 0)
            sage: P.is_unitary(algorithm='naive')
            True

            sage: P = matrix(CDF, 1, 1, [1])
            sage: P.is_unitary(algorithm='orthonormal')
            True

            sage: P = matrix(CDF, 0, 0,)
            sage: P.is_unitary(algorithm='orthonormal')
            True

        TESTS::

            sage: P = matrix(CDF, 2, 2)
            sage: P.is_unitary(tol='junk')
            Traceback (most recent call last):
            ...
            TypeError: tolerance must be a real number, not junk

            sage: P.is_unitary(tol=-0.3)
            Traceback (most recent call last):
            ...
            ValueError: tolerance must be positive, not -0.3

            sage: P.is_unitary(algorithm='junk')
            Traceback (most recent call last):
            ...
            ValueError: algorithm must be 'naive' or 'orthonormal', not junk


        AUTHOR:

        - Rob Beezer (2011-05-04)
        """
        if self.dimensions() == (0,0):
            # The "orthonormal" algorithm would otherwise fail in this
            # corner case. Returning ``True`` is consistent with the
            # other implementations of this method.
            return True

        global numpy
        try:
            tol = float(tol)
        except Exception:
            raise TypeError('tolerance must be a real number, not {0}'.format(tol))
        if tol <= 0:
            raise ValueError('tolerance must be positive, not {0}'.format(tol))
        if algorithm not in ['naive', 'orthonormal']:
            raise ValueError("algorithm must be 'naive' or 'orthonormal', not {0}".format(algorithm))
        key = 'unitary_{0}_{1}'.format(algorithm, tol)
        b = self.fetch(key)
        if b is not None:
            return b
        if not self.is_square():
            self.cache(key, False)
            return False
        if numpy is None:
            import numpy
        cdef Py_ssize_t i, j
        cdef Matrix_double_dense T, P
        if algorithm == 'orthonormal':
            # Schur decomposition over CDF will be unitary
            # iff diagonal with unit modulus entries
            _, T = self.schur(base_ring=sage.rings.complex_double.CDF)
            unitary = T._is_lower_triangular(tol)
            if unitary:
                for 0 <= i < self._nrows:
                    if abs(abs(T.get_unsafe(i,i)) - 1) > tol:
                        unitary = False
                        break
        elif algorithm == 'naive':
            P = self.conjugate().transpose()*self
            unitary = True
            for i from 0 <= i < self._nrows:
                # off-diagonal, since P is Hermitian
                for j from 0 <= j < i:
                    if abs(P.get_unsafe(i,j)) > tol:
                        unitary = False
                        break
                # at diagonal
                if abs(P.get_unsafe(i,i) - 1) > tol:
                    unitary = False
                if not unitary:
                    break
        self.cache(key, unitary)
        return unitary

    def _is_hermitian_orthonormal(self, tol=1e-12, skew=False):
        r"""
        Return ``True`` if the matrix is (skew-)Hermitian.

        For internal purposes. This function is used in `is_hermitian`
        and `is_skew_hermitian` functions.

        INPUT:

        - ``tol`` -- (default: ``1e-12``) the largest value of the
          absolute value of the difference between two matrix entries
          for which they will still be considered equal

        - ``skew`` -- (default: ``False``) specifies the type of the
          test. Set to ``True`` to check whether the matrix is
          skew-Hermitian.

        OUTPUT:

        ``True`` if the matrix is square and (skew-)Hermitian, and
        ``False`` otherwise.


        Note that if conjugation has no effect on elements of the base
        ring (such as for integers), then the :meth:`is_(skew_)symmetric`
        method is equivalent and faster.

        The tolerance parameter is used to allow for numerical values
        to be equal if there is a slight difference due to round-off
        and other imprecisions.

        The result is cached, on a per-tolerance basis.

        ALGORITHMS:

        The orthonormal algorithm first computes a Schur decomposition
        (via the :meth:`schur` method) and checks that the result is a
        diagonal matrix with real entries.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: A = matrix(CDF, [[ 1 + I,  1 - 6*I, -1 - I],
            ....:                  [-3 - I,     -4*I,     -2],
            ....:                  [-1 + I, -2 - 8*I,  2 + I]])
            sage: A._is_hermitian_orthonormal()
            False
            sage: B = A*A.conjugate_transpose()
            sage: B._is_hermitian_orthonormal()
            True

        A matrix that is nearly Hermitian, but for one non-real
        diagonal entry::

            sage: # needs sage.symbolic
            sage: A = matrix(CDF, [[    2,   2-I, 1+4*I],
            ....:                  [  2+I,   3+I, 2-6*I],
            ....:                  [1-4*I, 2+6*I,     5]])
            sage: A._is_hermitian_orthonormal()
            False
            sage: A[1,1] = 132
            sage: A._is_hermitian_orthonormal()
            True

        A square, empty matrix is trivially Hermitian::

            sage: A = matrix(RDF, 0, 0)
            sage: A._is_hermitian_orthonormal()
            True

        Rectangular matrices are never Hermitian::

            sage: A = matrix(CDF, 3, 4)
            sage: A._is_hermitian_orthonormal()
            False

        A matrix that is skew-Hermitian::

            sage: # needs sage.symbolic
            sage: A = matrix(CDF, [[-I, 2.0+I], [-2.0+I, 0.0]])
            sage: A._is_hermitian_orthonormal()
            False
            sage: A._is_hermitian_orthonormal(skew=True)
            True

        AUTHOR:

        - Rob Beezer (2011-03-30)
        """
        import sage.rings.complex_double
        global numpy
        tol = float(tol)

        key = ("_is_hermitian_orthonormal", tol, skew)
        h = self.fetch(key)
        if h is not None:
            return h
        if not self.is_square():
            self.cache(key, False)
            return False
        if self._nrows == 0:
            self.cache(key, True)
            return True
        if numpy is None:
            import numpy
        cdef Py_ssize_t i
        cdef Matrix_double_dense T
        # A matrix M is skew-hermitian iff I*M is hermitian
        T = self.__mul__(1j) if skew else self.__copy__()

        # Schur decomposition over CDF will be diagonal and real iff Hermitian
        _, T = T.schur(base_ring=sage.rings.complex_double.CDF)
        hermitian = T._is_lower_triangular(tol)
        if hermitian:
            for i in range(T._nrows):
                if abs(T.get_unsafe(i, i).imag()) > tol:
                    hermitian = False
                    break
        self.cache(key, hermitian)
        return hermitian

    def is_hermitian(self, tol=1e-12, algorithm="naive"):
        r"""
        Return ``True`` if the matrix is equal to its conjugate-transpose.

        INPUT:

        - ``tol`` -- (default: ``1e-12``) the largest value of the
          absolute value of the difference between two matrix entries
          for which they will still be considered equal.

        - ``algorithm`` -- string (default: ``'naive'``); either ``'naive'``
          or ``'orthonormal'``

        OUTPUT:

        ``True`` if the matrix is square and equal to the transpose with
        every entry conjugated, and ``False`` otherwise.

        Note that if conjugation has no effect on elements of the base
        ring (such as for integers), then the :meth:`is_symmetric`
        method is equivalent and faster.

        The tolerance parameter is used to allow for numerical values
        to be equal if there is a slight difference due to round-off
        and other imprecisions.

        The result is cached, on a per-tolerance and per-algorithm basis.

        ALGORITHMS:

        The naive algorithm simply compares corresponding entries on either
        side of the diagonal (and on the diagonal itself) to see if they are
        conjugates, with equality controlled by the tolerance parameter.

        The orthonormal algorithm first computes a Schur decomposition
        (via the :meth:`schur` method) and checks that the result is a
        diagonal matrix with real entries.

        So the naive algorithm can finish quickly for a matrix that is not
        Hermitian, while the orthonormal algorithm will always compute a
        Schur decomposition before going through a similar check of the matrix
        entry-by-entry.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: A = matrix(CDF, [[ 1 + I,  1 - 6*I, -1 - I],
            ....:                  [-3 - I,     -4*I,     -2],
            ....:                  [-1 + I, -2 - 8*I,  2 + I]])
            sage: A.is_hermitian(algorithm='orthonormal')
            False
            sage: A.is_hermitian(algorithm='naive')
            False
            sage: B = A*A.conjugate_transpose()
            sage: B.is_hermitian(algorithm='orthonormal')
            True
            sage: B.is_hermitian(algorithm='naive')
            True

        A matrix that is nearly Hermitian, but for one non-real
        diagonal entry. ::

            sage: # needs sage.symbolic
            sage: A = matrix(CDF, [[    2,   2-I, 1+4*I],
            ....:                  [  2+I,   3+I, 2-6*I],
            ....:                  [1-4*I, 2+6*I,     5]])
            sage: A.is_hermitian(algorithm='orthonormal')
            False
            sage: A[1,1] = 132
            sage: A.is_hermitian(algorithm='orthonormal')
            True

        We get a unitary matrix from the SVD routine and use this
        numerical matrix to create a matrix that should be Hermitian
        (indeed it should be the identity matrix), but with some
        imprecision.  We use this to illustrate that if the tolerance
        is set too small, then we can be too strict about the equality
        of entries and may achieve the wrong result (depending on
        the system)::

            sage: # needs sage.symbolic
            sage: A = matrix(CDF, [[ 1 + I,  1 - 6*I, -1 - I],
            ....:                  [-3 - I,     -4*I,     -2],
            ....:                  [-1 + I, -2 - 8*I,  2 + I]])
            sage: U, _, _ = A.SVD()
            sage: B = U*U.conjugate_transpose()
            sage: B.is_hermitian(algorithm='naive')
            True
            sage: B.is_hermitian(algorithm='naive', tol=1.0e-17)  # random
            False
            sage: B.is_hermitian(algorithm='naive', tol=1.0e-15)
            True

        A square, empty matrix is trivially Hermitian. ::

            sage: A = matrix(RDF, 0, 0)
            sage: A.is_hermitian()
            True

        Rectangular matrices are never Hermitian, no matter which
        algorithm is requested. ::

            sage: A = matrix(CDF, 3, 4)
            sage: A.is_hermitian()
            False

        TESTS:

        The ``algorithm`` keyword gets checked. ::

            sage: A = matrix(RDF, 2, range(4))
            sage: A.is_hermitian(algorithm='junk')
            Traceback (most recent call last):
            ...
            ValueError: algorithm must be 'naive' or 'orthonormal', not junk

        AUTHOR:

        - Rob Beezer (2011-03-30)
        """
        if algorithm == "naive":
            return super()._is_hermitian(skew=False, tolerance=tol)
        elif algorithm == "orthonormal":
            return self._is_hermitian_orthonormal(tol=tol, skew=False)
        else:
            raise ValueError("algorithm must be 'naive' or 'orthonormal', not {0}".format(algorithm))

    def is_skew_hermitian(self, tol=1e-12, algorithm='orthonormal'):
        r"""
        Return ``True`` if the matrix is equal to the negative of its
        conjugate transpose.

        INPUT:

        - ``tol`` -- (default: ``1e-12``) the largest value of the
          absolute value of the difference between two matrix entries
          for which they will still be considered equal.

        - ``algorithm`` -- (default: ``'orthonormal'``) set to
          ``'orthonormal'`` for a stable procedure and set to ``'naive'`` for a
          fast procedure

        OUTPUT:

        ``True`` if the matrix is square and equal to the negative of
        its conjugate transpose, and ``False`` otherwise.

        Note that if conjugation has no effect on elements of the base
        ring (such as for integers), then the :meth:`is_skew_symmetric`
        method is equivalent and faster.

        The tolerance parameter is used to allow for numerical values
        to be equal if there is a slight difference due to round-off
        and other imprecisions.

        The result is cached, on a per-tolerance and per-algorithm basis.

        ALGORITHMS:

        The naive algorithm simply compares corresponding entries on either
        side of the diagonal (and on the diagonal itself) to see if they are
        conjugates, with equality controlled by the tolerance parameter.

        The orthonormal algorithm first computes a Schur decomposition
        (via the :meth:`schur` method) and checks that the result is a
        diagonal matrix with real entries.

        So the naive algorithm can finish quickly for a matrix that is not
        Hermitian, while the orthonormal algorithm will always compute a
        Schur decomposition before going through a similar check of the matrix
        entry-by-entry.

        EXAMPLES::

            sage: A = matrix(CDF, [[0, -1],
            ....:                  [1,  0]])
            sage: A.is_skew_hermitian(algorithm='orthonormal')
            True
            sage: A.is_skew_hermitian(algorithm='naive')
            True

        A matrix that is nearly skew-Hermitian, but for a non-real
        diagonal entry. ::

            sage: # needs sage.symbolic
            sage: A = matrix(CDF, [[  -I, -1, 1-I],
            ....:                  [   1,  1,  -1],
            ....:                  [-1-I,  1,  -I]])
            sage: A.is_skew_hermitian()
            False
            sage: A[1,1] = -I
            sage: A.is_skew_hermitian()
            True

        We get a unitary matrix from the SVD routine and use this
        numerical matrix to create a matrix that should be
        skew-Hermitian (indeed it should be the identity matrix
        multiplied by `I`), but with some imprecision.  We use this to
        illustrate that if the tolerance is set too small, then we can
        be too strict about the equality of entries and may achieve
        the wrong result (depending on the system)::

            sage: # needs sage.symbolic
            sage: A = matrix(CDF, [[ 1 + I,  1 - 6*I, -1 - I],
            ....:                  [-3 - I,     -4*I,     -2],
            ....:                  [-1 + I, -2 - 8*I,  2 + I]])
            sage: U, _, _ = A.SVD()
            sage: B = 1j*U*U.conjugate_transpose()
            sage: B.is_skew_hermitian(algorithm='naive')
            True
            sage: B.is_skew_hermitian(algorithm='naive', tol=1.0e-17)  # random
            False
            sage: B.is_skew_hermitian(algorithm='naive', tol=1.0e-15)
            True

        A square, empty matrix is trivially Hermitian.  ::

            sage: A = matrix(RDF, 0, 0)
            sage: A.is_skew_hermitian()
            True

        Rectangular matrices are never Hermitian, no matter which
        algorithm is requested.  ::

            sage: A = matrix(CDF, 3, 4)
            sage: A.is_skew_hermitian()
            False

        TESTS:

        The ``algorithm`` keyword gets checked.  ::

            sage: A = matrix(RDF, 2, range(4))
            sage: A.is_skew_hermitian(algorithm='junk')
            Traceback (most recent call last):
            ...
            ValueError: algorithm must be 'naive' or 'orthonormal', not junk

        AUTHOR:

        - Rob Beezer (2011-03-30)
        """
        if algorithm == "naive":
            return super()._is_hermitian(skew=True, tolerance=tol)
        elif algorithm == "orthonormal":
            return self._is_hermitian_orthonormal(tol=tol, skew=True)
        else:
            raise ValueError("algorithm must be 'naive' or 'orthonormal', not {0}".format(algorithm))

    def is_normal(self, tol=1e-12, algorithm='orthonormal'):
        r"""
        Return ``True`` if the matrix commutes with its conjugate-transpose.

        INPUT:

        - ``tol`` -- (default: ``1e-12``) the largest value of the
          absolute value of the difference between two matrix entries
          for which they will still be considered equal.

        - ``algorithm`` -- (default: ``'orthonormal'``) set to
          ``'orthonormal'`` for a stable procedure and set to ``'naive'`` for a
          fast procedure

        OUTPUT:

        ``True`` if the matrix is square and commutes with its
        conjugate-transpose, and ``False`` otherwise.

        Normal matrices are precisely those that can be diagonalized
        by a unitary matrix.

        The tolerance parameter is used to allow for numerical values
        to be equal if there is a slight difference due to round-off
        and other imprecisions.

        The result is cached, on a per-tolerance and per-algorithm basis.

        ALGORITHMS:

        The naive algorithm simply compares entries of the two possible
        products of the matrix with its conjugate-transpose, with equality
        controlled by the tolerance parameter.

        The orthonormal algorithm first computes a Schur decomposition
        (via the :meth:`schur` method) and checks that the result is a
        diagonal matrix.  An orthonormal diagonalization
        is equivalent to being normal.

        So the naive algorithm can finish fairly quickly for a matrix
        that is not normal, once the products have been computed.
        However, the orthonormal algorithm will compute a Schur
        decomposition before going through a similar check of a
        matrix entry-by-entry.

        EXAMPLES:

        First over the complexes.  ``B`` is Hermitian, hence normal.  ::

            sage: # needs sage.symbolic
            sage: A = matrix(CDF, [[ 1 + I,  1 - 6*I, -1 - I],
            ....:                  [-3 - I,     -4*I,     -2],
            ....:                  [-1 + I, -2 - 8*I,  2 + I]])
            sage: B = A*A.conjugate_transpose()
            sage: B.is_hermitian()
            True
            sage: B.is_normal(algorithm='orthonormal')
            True
            sage: B.is_normal(algorithm='naive')
            True
            sage: B[0,0] = I
            sage: B.is_normal(algorithm='orthonormal')
            False
            sage: B.is_normal(algorithm='naive')
            False

        Now over the reals.  Circulant matrices are normal.  ::

            sage: # needs sage.graphs
            sage: G = graphs.CirculantGraph(20, [3, 7])
            sage: D = digraphs.Circuit(20)
            sage: A = 3*D.adjacency_matrix() - 5*G.adjacency_matrix()
            sage: A = A.change_ring(RDF)
            sage: A.is_normal()
            True
            sage: A.is_normal(algorithm='naive')
            True
            sage: A[19,0] = 4.0
            sage: A.is_normal()
            False
            sage: A.is_normal(algorithm='naive')
            False

        Skew-Hermitian matrices are normal.  ::

            sage: # needs sage.symbolic
            sage: A = matrix(CDF, [[ 1 + I,  1 - 6*I, -1 - I],
            ....:                  [-3 - I,     -4*I,     -2],
            ....:                  [-1 + I, -2 - 8*I,  2 + I]])
            sage: B = A - A.conjugate_transpose()
            sage: B.is_hermitian()
            False
            sage: B.is_normal()
            True
            sage: B.is_normal(algorithm='naive')
            True

        A small matrix that does not fit into any of the usual categories
        of normal matrices.  ::

            sage: A = matrix(RDF, [[1, -1],
            ....:                  [1,  1]])
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

        A square, empty matrix is trivially normal.  ::

            sage: A = matrix(CDF, 0, 0)
            sage: A.is_normal()
            True

        Rectangular matrices are never normal, no matter which
        algorithm is requested.  ::

            sage: A = matrix(CDF, 3, 4)
            sage: A.is_normal()
            False

        TESTS:

        The tolerance must be strictly positive.  ::

            sage: A = matrix(RDF, 2, range(4))
            sage: A.is_normal(tol = -3.1)
            Traceback (most recent call last):
            ...
            ValueError: tolerance must be positive, not -3.1

        The ``algorithm`` keyword gets checked.  ::

            sage: A = matrix(RDF, 2, range(4))
            sage: A.is_normal(algorithm='junk')
            Traceback (most recent call last):
            ...
            ValueError: algorithm must be 'naive' or 'orthonormal', not junk

        AUTHOR:

         - Rob Beezer (2011-03-31)
        """
        import sage.rings.complex_double
        global numpy
        tol = float(tol)
        if tol <= 0:
            raise ValueError('tolerance must be positive, not {0}'.format(tol))
        if algorithm not in ['naive', 'orthonormal']:
            raise ValueError("algorithm must be 'naive' or 'orthonormal', not {0}".format(algorithm))

        key = 'normal_{0}_{1}'.format(algorithm, tol)
        b = self.fetch(key)
        if b is not None:
            return b
        if not self.is_square():
            self.cache(key, False)
            return False
        if self._nrows == 0:
            self.cache(key, True)
            return True
        cdef Py_ssize_t i, j
        cdef Matrix_double_dense T, left, right
        if algorithm == 'orthonormal':
            # Schur decomposition over CDF will be diagonal iff normal
            _, T = self.schur(base_ring=sage.rings.complex_double.CDF)
            normal = T._is_lower_triangular(tol)
        elif algorithm == 'naive':
            if numpy is None:
                import numpy
            CT = self.conjugate_transpose()
            left = self*CT
            right = CT*self
            normal = True
            # two products are Hermitian, need only check lower triangle
            for i in range(self._nrows):
                for j in range(i+1):
                    if abs(left.get_unsafe(i,j) - right.get_unsafe(i,j)) > tol:
                        normal = False
                        break
                if not normal:
                    break
        self.cache(key, normal)
        return normal

    def schur(self, base_ring=None):
        r"""
        Return the Schur decomposition of the matrix.

        INPUT:

        - ``base_ring`` -- defaults to the base ring of ``self``; use this to
          request the base ring of the returned matrices, which will affect the
          format of the results

        OUTPUT:

        A pair of immutable matrices.  The first is a unitary matrix `Q`.
        The second, `T`, is upper-triangular when returned over the complex
        numbers, while it is almost upper-triangular over the reals.  In the
        latter case, there can be some `2\times 2` blocks on the diagonal
        which represent a pair of conjugate complex eigenvalues of ``self``.

        If ``self`` is the matrix `A`, then

        .. MATH::

            A = QT({\overline Q})^t

        where the latter matrix is the conjugate-transpose of ``Q``, which
        is also the inverse of ``Q``, since ``Q`` is unitary.

        Note that in the case of a normal matrix (Hermitian, symmetric, and
        others), the upper-triangular matrix is  a diagonal matrix with
        eigenvalues of ``self`` on the diagonal, and the unitary matrix
        has columns that form an orthonormal basis composed of eigenvectors
        of ``self``.  This is known as "orthonormal diagonalization".

        .. WARNING::

            The Schur decomposition is not unique, as there may be numerous
            choices for the vectors of the orthonormal basis, and consequently
            different possibilities for the upper-triangular matrix.  However,
            the diagonal of the upper-triangular matrix will always contain the
            eigenvalues of the matrix (in the complex version), or `2\times 2`
            block matrices in the real version representing pairs of conjugate
            complex eigenvalues.

            In particular, results may vary across systems and processors.

        EXAMPLES:

        First over the complexes.  The similar matrix is always
        upper-triangular in this case.  ::

            sage: # needs sage.symbolic
            sage: A = matrix(CDF, 4, 4, range(16)) + matrix(CDF, 4, 4,
            ....:                                           [x^3*I for x in range(0, 16)])
            sage: Q, T = A.schur()
            sage: (Q*Q.conjugate().transpose()).zero_at(1e-12)  # tol 1e-12
            [ 0.999999999999999                0.0                0.0                0.0]
            [               0.0 0.9999999999999996                0.0                0.0]
            [               0.0                0.0 0.9999999999999992                0.0]
            [               0.0                0.0                0.0 0.9999999999999999]
            sage: all(T.zero_at(1.0e-12)[i,j] == 0 for i in range(4) for j in range(i))
            True
            sage: (Q*T*Q.conjugate().transpose() - A).zero_at(1.0e-11)
            [0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0]
            sage: eigenvalues = [T[i,i] for i in range(4)]; eigenvalues
            [30.733... + 4648.541...*I, -0.184... - 159.057...*I, -0.523... + 11.158...*I, -0.025... - 0.642...*I]
            sage: A.eigenvalues()
            [30.733... + 4648.541...*I, -0.184... - 159.057...*I, -0.523... + 11.158...*I, -0.025... - 0.642...*I]
            sage: abs(A.norm()-T.norm()) < 1e-10
            True

        We begin with a real matrix but ask for a decomposition over the
        complexes.  The result will yield an upper-triangular matrix over
        the complex numbers for ``T``. ::

            sage: A = matrix(RDF, 4, 4, [x^3 for x in range(16)])
            sage: Q, T = A.schur(base_ring=CDF)
            sage: (Q*Q.conjugate().transpose()).zero_at(1e-12)  # tol 1e-12
            [0.9999999999999987                0.0                0.0                0.0]
            [               0.0 0.9999999999999999                0.0                0.0]
            [               0.0                0.0 1.0000000000000013                0.0]
            [               0.0                0.0                0.0 1.0000000000000007]
            sage: T.parent()
            Full MatrixSpace of 4 by 4 dense matrices over Complex Double Field
            sage: all(T.zero_at(1.0e-12)[i,j] == 0 for i in range(4) for j in range(i))
            True
            sage: (Q*T*Q.conjugate().transpose() - A).zero_at(1.0e-11)
            [0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0]

        Now totally over the reals.  But with complex eigenvalues, the
        similar matrix may not be upper-triangular. But "at worst" there
        may be some `2\times 2` blocks on the diagonal which represent
        a pair of conjugate complex eigenvalues. These blocks will then
        just interrupt the zeros below the main diagonal.  This example
        has a pair of these of the blocks. ::

            sage: A = matrix(RDF, 4, 4, [[1, 0, -3, -1],
            ....:                        [4, -16, -7, 0],
            ....:                        [1, 21, 1, -2],
            ....:                        [26, -1, -2, 1]])
            sage: Q, T = A.schur()
            sage: (Q*Q.conjugate().transpose())  # tol 1e-12
            [0.9999999999999994                0.0                0.0                0.0]
            [               0.0 1.0000000000000013                0.0                0.0]
            [               0.0                0.0 1.0000000000000004                0.0]
            [               0.0                0.0                0.0 1.0000000000000016]
            sage: all(T.zero_at(1.0e-12)[i,j] == 0 for i in range(4) for j in range(i))
            False
            sage: all(T.zero_at(1.0e-12)[i,j] == 0 for i in range(4) for j in range(i-1))
            True
            sage: (Q*T*Q.conjugate().transpose() - A).zero_at(1.0e-11)
            [0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0]
            sage: sorted(T[0:2,0:2].eigenvalues() + T[2:4,2:4].eigenvalues())
            [-5.710... - 8.382...*I, -5.710... + 8.382...*I, -0.789... - 2.336...*I, -0.789... + 2.336...*I]
            sage: sorted(A.eigenvalues())
            [-5.710... - 8.382...*I, -5.710... + 8.382...*I, -0.789... - 2.336...*I, -0.789... + 2.336...*I]
            sage: abs(A.norm()-T.norm()) < 1e-12
            True

        Starting with complex numbers and requesting a result over the reals
        will never happen.  ::

            sage: # needs sage.symbolic
            sage: A = matrix(CDF, 2, 2, [[2+I, -1+3*I], [5-4*I, 2-7*I]])
            sage: A.schur(base_ring=RDF)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert input matrix over CDF to a matrix over RDF

        If theory predicts your matrix is real, but it contains some
        very small imaginary parts, you can specify the cutoff for "small"
        imaginary parts, then request the output as real matrices, and let
        the routine do the rest. ::

            sage: A = matrix(RDF, 2, 2, [1, 1, -1, 0]) + matrix(CDF, 2, 2, [1.0e-14*I]*4)
            sage: B = A.zero_at(1.0e-12)
            sage: B.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Complex Double Field
            sage: Q, T = B.schur(RDF)
            sage: Q.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Real Double Field
            sage: T.parent()
            Full MatrixSpace of 2 by 2 dense matrices over Real Double Field
            sage: Q.round(6)
            [ 0.707107  0.707107]
            [-0.707107  0.707107]
            sage: T.round(6)
            [ 0.5  1.5]
            [-0.5  0.5]
            sage: (Q*T*Q.conjugate().transpose() - B).zero_at(1.0e-11)
            [0.0 0.0]
            [0.0 0.0]

        A Hermitian matrix has real eigenvalues, so the similar matrix
        will be upper-triangular.  Furthermore, a Hermitian matrix is
        diagonalizable with respect to an orthonormal basis, composed
        of eigenvectors of the matrix.  Here that basis is the set of
        columns of the unitary matrix.  ::

            sage: # needs sage.symbolic
            sage: A = matrix(CDF, [[        52,   -9*I - 8,    6*I - 187,  -188*I + 2],
            ....:                  [   9*I - 8,         12,   -58*I + 59,   30*I + 42],
            ....:                  [-6*I - 187,  58*I + 59,         2677, 2264*I + 65],
            ....:                  [ 188*I + 2, -30*I + 42, -2264*I + 65,       2080]])
            sage: Q, T = A.schur()
            sage: T = T.zero_at(1.0e-12).change_ring(RDF)
            sage: T.round(6)
            [4680.13301        0.0        0.0        0.0]
            [       0.0 102.715967        0.0        0.0]
            [       0.0        0.0  35.039344        0.0]
            [       0.0        0.0        0.0    3.11168]
            sage: (Q*Q.conjugate().transpose()).zero_at(1e-12)  # tol 1e-12
            [1.0000000000000004                0.0                0.0                0.0]
            [               0.0 0.9999999999999989                0.0                0.0]
            [               0.0                0.0 1.0000000000000002                0.0]
            [               0.0                0.0                0.0 0.9999999999999992]
            sage: (Q*T*Q.conjugate().transpose() - A).zero_at(1.0e-11)
            [0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0]

        Similarly, a real symmetric matrix has only real eigenvalues,
        and there is an orthonormal basis composed of eigenvectors of
        the matrix.  ::

            sage: A = matrix(RDF, [[ 1, -2, 5, -3],
            ....:                  [-2,  9, 1,  5],
            ....:                  [ 5,  1, 3 , 7],
            ....:                  [-3,  5, 7, -8]])
            sage: Q, T = A.schur()
            sage: Q.round(4)
            [-0.3027  -0.751   0.576 -0.1121]
            [  0.139 -0.3892 -0.2648  0.8713]
            [ 0.4361   0.359  0.7599  0.3217]
            [ -0.836  0.3945  0.1438  0.3533]
            sage: T = T.zero_at(10^-12)
            sage: all(abs(e) < 10^-4
            ....:     for e in (T - diagonal_matrix(RDF, [-13.5698, -0.8508, 7.7664, 11.6542])).list())
            True
            sage: (Q*Q.transpose())  # tol 1e-12
            [0.9999999999999998                0.0                0.0                0.0]
            [               0.0                1.0                0.0                0.0]
            [               0.0                0.0 0.9999999999999998                0.0]
            [               0.0                0.0                0.0 0.9999999999999996]
            sage: (Q*T*Q.transpose() - A).zero_at(1.0e-11)
            [0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0]
            [0.0 0.0 0.0 0.0]

        The results are cached, both as a real factorization and also as a
        complex factorization.  This means the returned matrices are
        immutable.  ::

            sage: A = matrix(RDF, 2, 2, [[0, -1], [1, 0]])
            sage: Qr, Tr = A.schur(base_ring=RDF)
            sage: Qc, Tc = A.schur(base_ring=CDF)
            sage: all(M.is_immutable() for M in [Qr, Tr, Qc, Tc])
            True
            sage: Tr.round(6) != Tc.round(6)
            True

        TESTS:

        The Schur factorization is only defined for square matrices.  ::

            sage: A = matrix(RDF, 4, 5, range(20))
            sage: A.schur()
            Traceback (most recent call last):
            ...
            ValueError: Schur decomposition requires a square matrix, not a 4 x 5 matrix

        A base ring request is checked. ::

            sage: A = matrix(RDF, 3, range(9))
            sage: A.schur(base_ring=QQ)
            Traceback (most recent call last):
            ...
            ValueError: base ring of Schur decomposition matrices must be RDF or CDF, not Rational Field

        AUTHOR:

        - Rob Beezer (2011-03-31)
        """
        global scipy
        from sage.rings.real_double import RDF
        from sage.rings.complex_double import CDF

        cdef Matrix_double_dense Q, T

        if not self.is_square():
            raise ValueError('Schur decomposition requires a square matrix, not a {0} x {1} matrix'.format(self.nrows(), self.ncols()))
        if base_ring is None:
            base_ring = self.base_ring()
        if base_ring not in [RDF, CDF]:
            raise ValueError('base ring of Schur decomposition matrices must be RDF or CDF, not {0}'.format(base_ring))

        if self.base_ring() != base_ring:
            try:
                self = self.change_ring(base_ring)
            except TypeError:
                raise TypeError('unable to convert input matrix over CDF to a matrix over RDF')
        if base_ring == CDF:
            format = 'complex'
        else:
            format = 'real'

        schur = self.fetch('schur_factors_' + format)
        if schur is not None:
            return schur
        if scipy is None:
            import scipy
        import scipy.linalg
        Q = self._new(self._nrows, self._nrows)
        T = self._new(self._nrows, self._nrows)
        T._matrix_numpy, Q._matrix_numpy = scipy.linalg.schur(self._matrix_numpy, output=format)
        Q.set_immutable()
        T.set_immutable()
        # Our return order is the reverse of NumPy's
        schur = (Q, T)
        self.cache('schur_factors_' + format, schur)
        return schur

    def cholesky(self):
        r"""
        Return the Cholesky factorization of a matrix that
        is real symmetric, or complex Hermitian.

        INPUT:

        Any square matrix with entries from ``RDF`` that is symmetric, or
        with entries from ``CDF`` that is Hermitian.  The matrix must
        be positive definite for the Cholesky decomposition to exist.

        OUTPUT:

        For a matrix `A` the routine returns a lower triangular
        matrix `L` such that,

        .. MATH::

            A = LL^\ast

        where `L^\ast` is the conjugate-transpose in the complex case,
        and just the transpose in the real case.  If the matrix fails
        to be positive definite (perhaps because it is not symmetric
        or Hermitian), then this function raises a :exc:`ValueError`.

        IMPLEMENTATION:

        The existence of a Cholesky decomposition and the
        positive definite property are equivalent.  So this
        method and the :meth:`is_positive_definite` method compute and
        cache both the Cholesky decomposition and the
        positive-definiteness.  So the :meth:`is_positive_definite`
        method or catching a :exc:`ValueError` from the :meth:`cholesky`
        method are equally expensive computationally and if the
        decomposition exists, it is cached as a side-effect of either
        routine.

        EXAMPLES:

        A real matrix that is symmetric, Hermitian, and positive definite::

            sage: M = matrix(RDF,[[ 1,  1,    1,     1,     1],
            ....:                 [ 1,  5,   31,   121,   341],
            ....:                 [ 1, 31,  341,  1555,  4681],
            ....:                 [ 1,121, 1555,  7381, 22621],
            ....:                 [ 1,341, 4681, 22621, 69905]])
            sage: M.is_symmetric()
            True
            sage: M.is_hermitian()
            True
            sage: L = M.cholesky()
            sage: L.round(6).zero_at(10^-10)
            [   1.0    0.0         0.0        0.0     0.0]
            [   1.0    2.0         0.0        0.0     0.0]
            [   1.0   15.0   10.723805        0.0     0.0]
            [   1.0   60.0   60.985814   7.792973     0.0]
            [   1.0  170.0  198.623524  39.366567  1.7231]
            sage: (L*L.transpose()).round(6).zero_at(10^-10)
            [ 1.0     1.0     1.0     1.0     1.0]
            [ 1.0     5.0    31.0   121.0   341.0]
            [ 1.0    31.0   341.0  1555.0  4681.0]
            [ 1.0   121.0  1555.0  7381.0 22621.0]
            [ 1.0   341.0  4681.0 22621.0 69905.0]

        A complex matrix that is Hermitian and positive definite.  ::

            sage: # needs sage.symbolic
            sage: A = matrix(CDF, [[        23,  17*I + 3,  24*I + 25,     21*I],
            ....:                  [ -17*I + 3,        38, -69*I + 89, 7*I + 15],
            ....:                  [-24*I + 25, 69*I + 89,        976, 24*I + 6],
            ....:                  [     -21*I, -7*I + 15,  -24*I + 6,       28]])
            sage: A.is_hermitian()
            True
            sage: L = A.cholesky()
            sage: L.round(6).zero_at(10^-10)
            [               4.795832                     0.0                    0.0       0.0]
            [  0.625543 - 3.544745*I                5.004346                    0.0       0.0]
            [   5.21286 - 5.004346*I 13.588189 + 10.721116*I              24.984023       0.0]
            [            -4.378803*I  -0.104257 - 0.851434*I  -0.21486 + 0.371348*I  2.811799]
            sage: (L*L.conjugate_transpose()).round(6).zero_at(10^-10)
            [         23.0  3.0 + 17.0*I 25.0 + 24.0*I        21.0*I]
            [ 3.0 - 17.0*I          38.0 89.0 - 69.0*I  15.0 + 7.0*I]
            [25.0 - 24.0*I 89.0 + 69.0*I         976.0  6.0 + 24.0*I]
            [      -21.0*I  15.0 - 7.0*I  6.0 - 24.0*I          28.0]

        This routine will recognize when the input matrix is not
        positive definite.  The negative eigenvalues are an
        equivalent indicator.  (Eigenvalues of a Hermitian matrix
        must be real, so there is no loss in ignoring the imprecise
        imaginary parts).  ::

            sage: A = matrix(RDF, [[ 3,  -6,   9,   6,  -9],
            ....:                  [-6,  11, -16, -11,  17],
            ....:                  [ 9, -16,  28,  16, -40],
            ....:                  [ 6, -11,  16,   9, -19],
            ....:                  [-9,  17, -40, -19,  68]])
            sage: A.is_symmetric()
            True
            sage: A.eigenvalues()
            [108.07..., 13.02..., -0.02..., -0.70..., -1.37...]
            sage: A.cholesky()
            Traceback (most recent call last):
            ...
            ValueError: matrix is not positive definite

            sage: # needs sage.symbolic
            sage: B = matrix(CDF, [[      2, 4 - 2*I, 2 + 2*I],
            ....:                  [4 + 2*I,       8,    10*I],
            ....:                  [2 - 2*I,   -10*I,      -3]])
            sage: B.is_hermitian()
            True
            sage: [ev.real() for ev in B.eigenvalues()]
            [15.88..., 0.08..., -8.97...]
            sage: B.cholesky()
            Traceback (most recent call last):
            ...
            ValueError: matrix is not positive definite

        TESTS:

        A trivial case. ::

            sage: A = matrix(RDF, 0, [])
            sage: A.cholesky()
            []

        The Cholesky factorization is only defined for Hermitian (in
        particular, square) matrices::

            sage: A = matrix(RDF, 4, 5, range(20))
            sage: A.cholesky()
            Traceback (most recent call last):
            ...
            ValueError: matrix is not Hermitian

        ::

            sage: # needs sage.symbolic
            sage: A = matrix(CDF, [[1+I]])
            sage: A.cholesky()
            Traceback (most recent call last):
            ...
            ValueError: matrix is not Hermitian
        """
        from sage.rings.real_double import RDF
        from sage.rings.complex_double import CDF

        cdef Matrix_double_dense L
        cache_cholesky = 'cholesky'
        cache_posdef = 'positive_definite'
        L = self.fetch(cache_cholesky)
        if L is not None:
            return L

        if not self.is_hermitian():
            self.cache(cache_posdef, False)
            raise ValueError("matrix is not Hermitian")

        if self._nrows == 0:   # special case
            self.cache(cache_posdef, True)
            L = self.__copy__()
            L.set_immutable()
            return L

        L = self._new()
        from scipy.linalg import cholesky
        from numpy.linalg import LinAlgError
        try:
            L._matrix_numpy = cholesky(self._matrix_numpy, lower=1)
        except LinAlgError:
            self.cache(cache_posdef, False)
            raise ValueError("matrix is not positive definite")
        L.set_immutable()
        self.cache(cache_cholesky, L)
        self.cache(cache_posdef, True)

        return L

    def is_positive_definite(self):
        r"""
        Determines if a matrix is positive definite.

        A matrix `A` is positive definite if it is square,
        is Hermitian (which reduces to symmetric in the real case),
        and for every nonzero vector `\vec{x}`,

        .. MATH::

            \vec{x}^\ast A \vec{x} > 0

        where `\vec{x}^\ast` is the conjugate-transpose in the
        complex case and just the transpose in the real case.
        Equivalently, a positive definite matrix has only positive
        eigenvalues and only positive determinants of leading
        principal submatrices.

        Applies to any matrix over ``RDF`` or ``CDF``.

        OUTPUT:

        ``True`` if and only if the matrix is square, Hermitian,
        and meets the condition above on the quadratic form.
        The result is cached.

        IMPLEMENTATION:

        The existence of a Cholesky decomposition and the
        positive definite property are equivalent.  So this
        method and the :meth:`cholesky` method compute and
        cache both the Cholesky decomposition and the
        positive-definiteness.  So the :meth:`is_positive_definite`
        method or catching a :exc:`ValueError` from the :meth:`cholesky`
        method are equally expensive computationally and if the
        decomposition exists, it is cached as a side-effect of either
        routine.

        EXAMPLES:

        A matrix over ``RDF`` that is positive definite.  ::

            sage: M = matrix(RDF,[[ 1,  1,    1,     1,     1],
            ....:                 [ 1,  5,   31,   121,   341],
            ....:                 [ 1, 31,  341,  1555,  4681],
            ....:                 [ 1,121, 1555,  7381, 22621],
            ....:                 [ 1,341, 4681, 22621, 69905]])
            sage: M.is_symmetric()
            True
            sage: M.eigenvalues()
            [77547.66..., 82.44..., 2.41..., 0.46..., 0.011...]
            sage: [round(M[:i,:i].determinant()) for i in range(1, M.nrows()+1)]
            [1, 4, 460, 27936, 82944]
            sage: M.is_positive_definite()
            True

        A matrix over ``CDF`` that is positive definite.  ::

            sage: # needs sage.symbolic
            sage: C = matrix(CDF, [[        23,  17*I + 3,  24*I + 25,     21*I],
            ....:                  [ -17*I + 3,        38, -69*I + 89, 7*I + 15],
            ....:                  [-24*I + 25, 69*I + 89,        976, 24*I + 6],
            ....:                  [     -21*I, -7*I + 15,  -24*I + 6,       28]])
            sage: C.is_hermitian()
            True
            sage: [x.real() for x in C.eigenvalues()]
            [991.46..., 55.96..., 3.69..., 13.87...]
            sage: [round(C[:i,:i].determinant().real()) for i in range(1, C.nrows()+1)]
            [23, 576, 359540, 2842600]
            sage: C.is_positive_definite()
            True

        A matrix over ``RDF`` that is not positive definite.  ::

            sage: A = matrix(RDF, [[ 3,  -6,   9,   6,  -9],
            ....:                  [-6,  11, -16, -11,  17],
            ....:                  [ 9, -16,  28,  16, -40],
            ....:                  [ 6, -11,  16,   9, -19],
            ....:                  [-9,  17, -40, -19,  68]])
            sage: A.is_symmetric()
            True
            sage: A.eigenvalues()
            [108.07..., 13.02..., -0.02..., -0.70..., -1.37...]
            sage: [round(A[:i,:i].determinant()) for i in range(1, A.nrows()+1)]
            [3, -3, -15, 30, -30]
            sage: A.is_positive_definite()
            False

        A matrix over ``CDF`` that is not positive definite.  ::

            sage: # needs sage.symbolic
            sage: B = matrix(CDF, [[      2, 4 - 2*I, 2 + 2*I],
            ....:                  [4 + 2*I,       8,    10*I],
            ....:                  [2 - 2*I,   -10*I,      -3]])
            sage: B.is_hermitian()
            True
            sage: [ev.real() for ev in B.eigenvalues()]
            [15.88..., 0.08..., -8.97...]
            sage: [round(B[:i,:i].determinant().real()) for i in range(1, B.nrows()+1)]
            [2, -4, -12]
            sage: B.is_positive_definite()
            False

        A large random matrix that is guaranteed by theory to be
        positive definite. ::

            sage: R = random_matrix(CDF, 200)
            sage: H = R.conjugate_transpose()*R
            sage: H.is_positive_definite()
            True

        TESTS:

        A trivially small case.  ::

            sage: S = matrix(CDF, [])
            sage: S.nrows(), S.ncols()
            (0, 0)
            sage: S.is_positive_definite()
            True

        A rectangular matrix will never be positive definite.  ::

            sage: R = matrix(RDF, 2, 3, range(6))
            sage: R.is_positive_definite()
            False

        A non-Hermitian matrix will never be positive definite::

            sage: T = matrix(CDF, 8, 8, range(64))
            sage: T.is_positive_definite()
            False

        ::

            sage: # needs sage.symbolic
            sage: A = matrix(CDF, [[1+I]])
            sage: A.is_positive_definite()
            False

        AUTHOR:

        - Rob Beezer (2012-05-28)
        """
        cache_str = 'positive_definite'
        posdef = self.fetch(cache_str)
        if posdef is None:
            try:
                self.cholesky()
            except ValueError:
                pass
            posdef = self.fetch(cache_str)
        return posdef

    cdef _vector_times_matrix_(self, Vector v):
        if self._nrows == 0 or self._ncols == 0:
            return self.row_ambient_module().zero_vector()
        global numpy
        if numpy is None:
            import numpy

        v_numpy = numpy.array([self._python_dtype(i) for i in v])

        M = self.row_ambient_module()
        ans = numpy.dot(v_numpy,self._matrix_numpy)
        return M(ans)

    cdef _matrix_times_vector_(self, Vector v):
        if self._nrows == 0 or self._ncols == 0:
            return self.column_ambient_module().zero_vector()

        global numpy
        if numpy is None:
            import numpy

        v_numpy = numpy.array([self._python_dtype(i) for i in v], dtype=self._numpy_dtype)

        M = self.column_ambient_module()
        ans = numpy.dot(self._matrix_numpy, v_numpy)
        return M(ans)

    def _replace_self_with_numpy32(self, numpy_matrix):
        """

        EXAMPLES::

            sage: import numpy
            sage: a = numpy.array([[1,2],[3,4]], 'float32')
            sage: m = matrix(RDF,2,2,0)
            sage: m._replace_self_with_numpy32(a)
            sage: m
            [1.0 2.0]
            [3.0 4.0]
        """
        # TODO find where this is used and change it
        self._replace_self_with_numpy(numpy_matrix)

    def _hadamard_row_bound(self):
        r"""
        Return an integer n such that the absolute value of the
        determinant of this matrix is at most `10^n`.

        EXAMPLES::

            sage: a = matrix(RDF, 3, [1,2,5,7,-3,4,2,1,123])
            sage: a._hadamard_row_bound()
            4
            sage: a.det()
            -2014.0
            sage: 10^4
            10000
        """
        cdef double d = 0, s
        cdef Py_ssize_t i, j
        for i from 0 <= i < self._nrows:
            s = 0
            for j from 0 <= j < self._ncols:
                s += self.get_unsafe(i, j)**2
            d += math.log(s)
        d /= 2
        return int(math.ceil(d / math.log(10)))

    def exp(self):
        r"""
        Calculate the exponential of this matrix X, which is the matrix

        .. MATH::

            e^X = \sum_{k=0}^{\infty} \frac{X^k}{k!}.

        EXAMPLES::

            sage: A = matrix(RDF, 2, [1,2,3,4]); A
            [1.0 2.0]
            [3.0 4.0]
            sage: A.exp()  # tol 5e-14
            [51.968956198705044  74.73656456700327]
            [112.10484685050491 164.07380304920997]
            sage: A = matrix(CDF, 2, [1,2+I,3*I,4]); A                                  # needs sage.symbolic
            [        1.0 2.0 + 1.0*I]
            [      3.0*I         4.0]
            sage: A.exp()  # tol 3e-14                                                  # needs sage.symbolic
            [-19.614602953804912 + 12.517743846762578*I   3.7949636449582176 + 28.88379930658099*I]
            [ -32.383580980922254 + 21.88423595789845*I   2.269633004093535 + 44.901324827684824*I]

        TESTS::

            sage: A = matrix(RDF, 2, [1,2,3,4])
            sage: A.exp()   # tol 5e-14
            [51.968956198705044  74.73656456700327]
            [112.10484685050491 164.07380304920997]

            sage: A = matrix(CDF, 2, [1,2+I,3*I,4])                                     # needs sage.symbolic
            sage: A.exp()  # tol 3e-14                                                  # needs sage.symbolic
            [-19.614602953804923 + 12.51774384676257*I 3.7949636449582016 + 28.883799306580997*I]
            [-32.38358098092227 + 21.884235957898433*I  2.2696330040935084 + 44.90132482768484*I]
        """
        global scipy
        if scipy is None:
            import scipy
        import scipy.linalg

        cdef Matrix_double_dense M
        M = self._new()
        M._matrix_numpy = scipy.linalg.expm(self._matrix_numpy)
        return M

    def zero_at(self, eps):
        """
        Return a copy of the matrix where elements smaller than or
        equal to ``eps`` are replaced with zeroes. For complex matrices,
        the real and imaginary parts are considered individually.

        This is useful for modifying output from algorithms which have large
        relative errors when producing zero elements, e.g. to create reliable
        doctests.

        INPUT:

        - ``eps`` -- cutoff value

        OUTPUT: a modified copy of the matrix

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: a = matrix(CDF, [[1, 1e-4r, 1+1e-100jr], [1e-8+3j, 0, 1e-58r]])
            sage: a
            [           1.0         0.0001 1.0 + 1e-100*I]
            [ 1e-08 + 3.0*I            0.0          1e-58]
            sage: a.zero_at(1e-50)
            [          1.0        0.0001           1.0]
            [1e-08 + 3.0*I           0.0           0.0]
            sage: a.zero_at(1e-4)
            [  1.0   0.0   1.0]
            [3.0*I   0.0   0.0]
        """
        global numpy
        cdef Matrix_double_dense M
        if numpy is None:
            import numpy
        eps = float(eps)
        out = self._matrix_numpy.copy()
        if self._sage_dtype is sage.rings.complex_double.CDF:
            out.real[numpy.abs(out.real) <= eps] = 0
            out.imag[numpy.abs(out.imag) <= eps] = 0
        else:
            out[numpy.abs(out) <= eps] = 0
        M = self._new()
        M._matrix_numpy = out
        return M

    def round(self, ndigits=0):
        """
        Return a copy of the matrix where all entries have been rounded
        to a given precision in decimal digits (default: 0 digits).

        INPUT:

        - ``ndigits`` -- the precision in number of decimal digits

        OUTPUT: a modified copy of the matrix

        EXAMPLES::

            sage: M = matrix(CDF, [[10.234r + 34.2343jr, 34e10r]])
            sage: M
            [10.234 + 34.2343*I     340000000000.0]
            sage: M.round(2)
            [10.23 + 34.23*I  340000000000.0]
            sage: M.round()
            [ 10.0 + 34.0*I 340000000000.0]
        """
        global numpy
        cdef Matrix_double_dense M
        if numpy is None:
            import numpy
        ndigits = int(ndigits)
        M = self._new()
        M._matrix_numpy = numpy.round(self._matrix_numpy, ndigits)
        return M

    def _normalize_columns(self):
        """
        Return a copy of the matrix where each column has been
        multiplied by plus or minus 1, to guarantee that the real
        part of the leading entry of each nonzero column is positive.

        This is useful for modifying output from algorithms which
        produce matrices which are only well-defined up to signs of
        the columns, for example an algorithm which should produce an
        orthogonal matrix.

        OUTPUT: a modified copy of the matrix

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: a = matrix(CDF, [[1, -2+I, 0, -3*I], [2, 2, -2, 2], [-3, -3, -3, -2]])
            sage: a
            [         1.0 -2.0 + 1.0*I          0.0       -3.0*I]
            [         2.0          2.0         -2.0          2.0]
            [        -3.0         -3.0         -3.0         -2.0]
            sage: a._normalize_columns()
            [        1.0 2.0 - 1.0*I         0.0      -3.0*I]
            [        2.0        -2.0         2.0         2.0]
            [       -3.0         3.0         3.0        -2.0]
        """
        M = self.__copy__()
        cdef Py_ssize_t i, j
        for j from 0 <= j < M.ncols():
            for i from 0 <= i < M.column(j).degree():
                a = M.column(j)[i].real()
                if a != 0:
                    if a < 0:
                        M.rescale_col(j, -1)
                    break
        return M

    def _normalize_rows(self):
        """
        Return a copy of the matrix where each row has been
        multiplied by plus or minus 1, to guarantee that the real
        part of the leading entry of each nonzero row is positive.

        This is useful for modifying output from algorithms which
        produce matrices which are only well-defined up to signs of
        the rows, for example an algorithm which should produce an
        upper triangular matrix.

        OUTPUT: a modified copy of the matrix

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: a = matrix(CDF, [[1, 2, -3], [-2+I, 2, -3], [0, -2, -3], [-3*I, 2, -2]])
            sage: a
            [         1.0          2.0         -3.0]
            [-2.0 + 1.0*I          2.0         -3.0]
            [         0.0         -2.0         -3.0]
            [      -3.0*I          2.0         -2.0]
            sage: a._normalize_rows()
            [        1.0         2.0        -3.0]
            [2.0 - 1.0*I        -2.0         3.0]
            [        0.0         2.0         3.0]
            [     -3.0*I         2.0        -2.0]
        """
        return self.transpose()._normalize_columns().transpose()
