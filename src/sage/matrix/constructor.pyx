# sage_setup: distribution = sagemath-modules
# cython: binding=True
"""
General matrix Constructor and display options
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#       Copyright (C) 2016 Jeroen Demeyer <J.Demeyer@UGent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.matrix.args cimport MatrixArgs
from sage.structure.global_options import GlobalOptions


def matrix(*args, **kwds):
    """
    Create a matrix.

    This implements the ``matrix`` constructor::

        sage: matrix([[1,2],[3,4]])
        [1 2]
        [3 4]

    It also contains methods to create special types of matrices, see
    ``matrix.[tab]`` for more options. For example::

        sage: matrix.identity(2)
        [1 0]
        [0 1]

    INPUT:

    The :func:`matrix` command takes the entries of a matrix, optionally
    preceded by a ring and the dimensions of the matrix, and returns a
    matrix.

    The entries of a matrix can be specified as a flat list of
    elements, a list of lists (i.e., a list of rows), a list of Sage
    vectors, a callable object, or a dictionary having positions as
    keys and matrix entries as values (see the examples). If you pass
    in a callable object, then you must specify the number of rows and
    columns. You can create a matrix of zeros by passing an empty list
    or the integer zero for the entries.  To construct a multiple of
    the identity (`cI`), you can specify square dimensions and pass in
    `c`. Calling :func:`matrix` with a Sage object may return something that
    makes sense. Calling :func:`matrix` with a NumPy array will convert the
    array to a matrix.

    All arguments (even the positional ones) are optional.

    Positional and keyword arguments:

    - ``base_ring`` -- parent of the entries of the matrix (despite the
      name, this is not a priori required to be a ring). By default,
      determine this from the given entries, falling back to ``ZZ`` if
      no entries are given.

    - ``nrows`` -- the number of rows in the matrix, or a finite or
      enumerated family of arbitrary objects that index the rows of the matrix

    - ``ncols`` -- the number of columns in the matrix, or a finite or
      enumerated family of arbitrary objects that index the columns of the matrix

    - ``entries`` -- see examples below

    If any of ``nrows``, ``ncols``, ``row_keys``, ``column_keys`` is
    given as keyword argument, then none of these may be given as
    positional arguments.

    Keyword-only arguments:

    - ``sparse`` -- boolean; create a sparse matrix. This defaults to
      ``True`` when the entries are given as a dictionary, otherwise
      defaults to ``False``.

    - ``row_keys`` -- a finite or enumerated family of arbitrary objects
      that index the rows of the matrix

    - ``column_keys`` -- a finite or enumerated family of arbitrary objects
      that index the columns of the matrix

    - ``space`` -- matrix space which will be the parent of the output
      matrix. This determines ``base_ring``, ``nrows``, ``row_keys``,
      ``ncols``, ``column_keys``, and ``sparse``.

    - ``immutable`` -- boolean; make the matrix immutable. By default,
      the output matrix is mutable.

    OUTPUT: a matrix or, more generally, a homomorphism between free
    modules

    EXAMPLES::

        sage: m = matrix(2); m; m.parent()
        [0 0]
        [0 0]
        Full MatrixSpace of 2 by 2 dense matrices over Integer Ring

    ::

        sage: m = matrix(2, 3); m; m.parent()
        [0 0 0]
        [0 0 0]
        Full MatrixSpace of 2 by 3 dense matrices over Integer Ring

    ::

        sage: m = matrix(QQ, [[1,2,3], [4,5,6]]); m; m.parent()
        [1 2 3]
        [4 5 6]
        Full MatrixSpace of 2 by 3 dense matrices over Rational Field

    ::

        sage: m = matrix(QQ, 3, 3, lambda i, j: i + j); m
        [0 1 2]
        [1 2 3]
        [2 3 4]
        sage: m = matrix(3, lambda i, j: i - j); m
        [ 0 -1 -2]
        [ 1  0 -1]
        [ 2  1  0]

    ::

        sage: matrix(QQ, 2, 3, lambda x, y: x + y)
        [0 1 2]
        [1 2 3]
        sage: matrix(QQ, 5, 5, lambda x, y: (x+1) / (y+1))
        [  1 1/2 1/3 1/4 1/5]
        [  2   1 2/3 1/2 2/5]
        [  3 3/2   1 3/4 3/5]
        [  4   2 4/3   1 4/5]
        [  5 5/2 5/3 5/4   1]

    ::

        sage: v1 = vector((1,2,3))
        sage: v2 = vector((4,5,6))
        sage: m = matrix([v1,v2]); m; m.parent()
        [1 2 3]
        [4 5 6]
        Full MatrixSpace of 2 by 3 dense matrices over Integer Ring

    ::

        sage: m = matrix(QQ, 2, [1,2,3,4,5,6]); m; m.parent()
        [1 2 3]
        [4 5 6]
        Full MatrixSpace of 2 by 3 dense matrices over Rational Field

    ::

        sage: m = matrix(QQ, 2, 3, [1,2,3,4,5,6]); m; m.parent()
        [1 2 3]
        [4 5 6]
        Full MatrixSpace of 2 by 3 dense matrices over Rational Field

    ::

        sage: m = matrix({(0,1): 2, (1,1): 2/5}); m; m.parent()
        [  0   2]
        [  0 2/5]
        Full MatrixSpace of 2 by 2 sparse matrices over Rational Field

    ::

        sage: m = matrix(QQ, 2, 3, {(1,1): 2}); m; m.parent()
        [0 0 0]
        [0 2 0]
        Full MatrixSpace of 2 by 3 sparse matrices over Rational Field

    ::

        sage: import numpy                                                              # needs numpy
        sage: n = numpy.array([[1,2], [3,4]], float)                                    # needs numpy
        sage: m = matrix(n); m; m.parent()                                              # needs numpy
        [1.0 2.0]
        [3.0 4.0]
        Full MatrixSpace of 2 by 2 dense matrices over Real Double Field

    ::

        sage: v = vector(ZZ, [1, 10, 100])
        sage: m = matrix(v); m; m.parent()
        [  1  10 100]
        Full MatrixSpace of 1 by 3 dense matrices over Integer Ring
        sage: m = matrix(GF(7), v); m; m.parent()
        [1 3 2]
        Full MatrixSpace of 1 by 3 dense matrices over Finite Field of size 7
        sage: m = matrix(GF(7), 3, 1, v); m; m.parent()
        [1]
        [3]
        [2]
        Full MatrixSpace of 3 by 1 dense matrices over Finite Field of size 7

    ::

        sage: matrix(pari.mathilbert(3))                                                # needs sage.libs.pari
        [  1 1/2 1/3]
        [1/2 1/3 1/4]
        [1/3 1/4 1/5]

    ::

        sage: g = graphs.PetersenGraph()                                                # needs sage.graphs
        sage: m = matrix(g); m; m.parent()                                              # needs sage.graphs
        [0 1 0 0 1 1 0 0 0 0]
        [1 0 1 0 0 0 1 0 0 0]
        [0 1 0 1 0 0 0 1 0 0]
        [0 0 1 0 1 0 0 0 1 0]
        [1 0 0 1 0 0 0 0 0 1]
        [1 0 0 0 0 0 0 1 1 0]
        [0 1 0 0 0 0 0 0 1 1]
        [0 0 1 0 0 1 0 0 0 1]
        [0 0 0 1 0 1 1 0 0 0]
        [0 0 0 0 1 0 1 1 0 0]
        Full MatrixSpace of 10 by 10 dense matrices over Integer Ring

    ::

        sage: matrix(ZZ, 10, 10, range(100), sparse=True).parent()
        Full MatrixSpace of 10 by 10 sparse matrices over Integer Ring

    ::

        sage: R = PolynomialRing(QQ, 9, 'x')
        sage: A = matrix(R, 3, 3, R.gens()); A
        [x0 x1 x2]
        [x3 x4 x5]
        [x6 x7 x8]
        sage: det(A)
        -x2*x4*x6 + x1*x5*x6 + x2*x3*x7 - x0*x5*x7 - x1*x3*x8 + x0*x4*x8

    ::

        sage: M = Matrix([[1,2,3], [4,5,6], [7,8,9]], immutable=True)
        sage: M[0] = [9,9,9]
        Traceback (most recent call last):
        ...
        ValueError: matrix is immutable; please change a copy instead
        (i.e., use copy(M) to change a copy of M).

    Using ``row_keys`` and ``column_keys``::

        sage: M = matrix([[1,2,3], [4,5,6]],
        ....:            column_keys=['a','b','c'], row_keys=['u','v']); M
        Generic morphism:
          From: Free module generated by {'a', 'b', 'c'} over Integer Ring
          To:   Free module generated by {'u', 'v'} over Integer Ring
        sage: print(M._unicode_art_matrix())
          a b c
        u⎛1 2 3⎞
        v⎝4 5 6⎠

    It is allowed to specify dimensions redundantly::

        sage: M = matrix(2, 3, [[1,2,3], [4,5,6]],
        ....:            column_keys=['a','b','c'], row_keys=['u','v']); M
        Generic morphism:
        From: Free module generated by {'a', 'b', 'c'} over Integer Ring
        To:   Free module generated by {'u', 'v'} over Integer Ring

    TESTS:

    There are many ways to create an empty matrix::

        sage: m = matrix(); m; m.parent()
        []
        Full MatrixSpace of 0 by 0 dense matrices over Integer Ring
        sage: m = matrix(sparse=True); m; m.parent()
        []
        Full MatrixSpace of 0 by 0 sparse matrices over Integer Ring
        sage: m = matrix(QQ); m; m.parent()
        []
        Full MatrixSpace of 0 by 0 dense matrices over Rational Field
        sage: m = matrix(base_ring=QQ); m; m.parent()
        []
        Full MatrixSpace of 0 by 0 dense matrices over Rational Field
        sage: m = matrix(0); m; m.parent()
        []
        Full MatrixSpace of 0 by 0 dense matrices over Integer Ring
        sage: m = matrix(0, 0, base_ring=QQ); m; m.parent()
        []
        Full MatrixSpace of 0 by 0 dense matrices over Rational Field
        sage: m = matrix([]); m; m.parent()
        []
        Full MatrixSpace of 0 by 0 dense matrices over Integer Ring
        sage: m = matrix(QQ, []); m; m.parent()
        []
        Full MatrixSpace of 0 by 0 dense matrices over Rational Field
        sage: m = matrix(QQ, {}); m; m.parent()
        []
        Full MatrixSpace of 0 by 0 sparse matrices over Rational Field

    Only a ring and dimensions::

        sage: m = matrix(2); m; m.parent()
        [0 0]
        [0 0]
        Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
        sage: m = matrix(QQ, 2); m; m.parent()
        [0 0]
        [0 0]
        Full MatrixSpace of 2 by 2 dense matrices over Rational Field
        sage: m = matrix(QQ, 2, 3); m; m.parent()
        [0 0 0]
        [0 0 0]
        Full MatrixSpace of 2 by 3 dense matrices over Rational Field

    A ring, dimensions and a scalar::

        sage: m = matrix(2,2,1); m; m.parent()
        [1 0]
        [0 1]
        Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
        sage: m = matrix(QQ,2,2,5); m; m.parent()
        [5 0]
        [0 5]
        Full MatrixSpace of 2 by 2 dense matrices over Rational Field

    For non-square matrices, only zero works::

        sage: m = matrix(2,3,0); m; m.parent()
        [0 0 0]
        [0 0 0]
        Full MatrixSpace of 2 by 3 dense matrices over Integer Ring
        sage: m = matrix(QQ,2,3,0); m; m.parent()
        [0 0 0]
        [0 0 0]
        Full MatrixSpace of 2 by 3 dense matrices over Rational Field
        sage: matrix(QQ,2,3,1)
        Traceback (most recent call last):
        ...
        TypeError: nonzero scalar matrix must be square
        sage: matrix(QQ,2,3,5)
        Traceback (most recent call last):
        ...
        TypeError: nonzero scalar matrix must be square

    Matrices specified by a list of entries::

        sage: m = matrix([[1,2,3],[4,5,6]]); m; m.parent()
        [1 2 3]
        [4 5 6]
        Full MatrixSpace of 2 by 3 dense matrices over Integer Ring
        sage: m = matrix(QQ, 2, [[1,2,3],[4,5,6]]); m; m.parent()
        [1 2 3]
        [4 5 6]
        Full MatrixSpace of 2 by 3 dense matrices over Rational Field
        sage: m = matrix(QQ, 3, [[1,2,3],[4,5,6]]); m; m.parent()
        Traceback (most recent call last):
        ...
        ValueError: inconsistent number of rows: should be 3 but got 2
        sage: m = matrix(QQ, 2, 3, [[1,2,3],[4,5,6]]); m; m.parent()
        [1 2 3]
        [4 5 6]
        Full MatrixSpace of 2 by 3 dense matrices over Rational Field
        sage: m = matrix(QQ, 2, 4, [[1,2,3],[4,5,6]]); m; m.parent()
        Traceback (most recent call last):
        ...
        ValueError: sequence too short (expected length 4, got 3)
        sage: m = matrix([(1,2,3),(4,5,6)]); m; m.parent()
        [1 2 3]
        [4 5 6]
        Full MatrixSpace of 2 by 3 dense matrices over Integer Ring
        sage: m = matrix([1,2,3,4,5,6]); m; m.parent()
        [1 2 3 4 5 6]
        Full MatrixSpace of 1 by 6 dense matrices over Integer Ring
        sage: m = matrix((1,2,3,4,5,6)); m; m.parent()
        [1 2 3 4 5 6]
        Full MatrixSpace of 1 by 6 dense matrices over Integer Ring
        sage: m = matrix(QQ, [1,2,3,4,5,6]); m; m.parent()
        [1 2 3 4 5 6]
        Full MatrixSpace of 1 by 6 dense matrices over Rational Field
        sage: m = matrix(QQ, 3, 2, [1,2,3,4,5,6]); m; m.parent()
        [1 2]
        [3 4]
        [5 6]
        Full MatrixSpace of 3 by 2 dense matrices over Rational Field
        sage: m = matrix(QQ, 2, 4, [1,2,3,4,5,6]); m; m.parent()
        Traceback (most recent call last):
        ...
        ValueError: sequence too short (expected length 8, got 6)
        sage: m = matrix(QQ, 5, [1,2,3,4,5,6]); m; m.parent()
        Traceback (most recent call last):
        ...
        ValueError: sequence too long (expected length 5, got more)

    Matrices specified by a dict of entries::

        sage: m = matrix({(1,1): 2}); m; m.parent()
        [0 0]
        [0 2]
        Full MatrixSpace of 2 by 2 sparse matrices over Integer Ring
        sage: m = matrix({(1,1): 2}, sparse=False); m; m.parent()
        [0 0]
        [0 2]
        Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
        sage: m = matrix(QQ, {(1,1): 2}); m; m.parent()
        [0 0]
        [0 2]
        Full MatrixSpace of 2 by 2 sparse matrices over Rational Field
        sage: m = matrix(QQ, 3, {(1,1): 2}); m; m.parent()
        [0 0 0]
        [0 2 0]
        [0 0 0]
        Full MatrixSpace of 3 by 3 sparse matrices over Rational Field
        sage: m = matrix(QQ, 3, 4, {(1,1): 2}); m; m.parent()
        [0 0 0 0]
        [0 2 0 0]
        [0 0 0 0]
        Full MatrixSpace of 3 by 4 sparse matrices over Rational Field
        sage: m = matrix(QQ, 2, {(1,1): 2}); m; m.parent()
        [0 0]
        [0 2]
        Full MatrixSpace of 2 by 2 sparse matrices over Rational Field
        sage: m = matrix(QQ, 1, {(1,1): 2}); m; m.parent()
        Traceback (most recent call last):
        ...
        IndexError: invalid row index 1
        sage: m = matrix({}); m; m.parent()
        []
        Full MatrixSpace of 0 by 0 sparse matrices over Integer Ring
        sage: m = matrix(QQ, {}); m; m.parent()
        []
        Full MatrixSpace of 0 by 0 sparse matrices over Rational Field
        sage: m = matrix(QQ, 2, {}); m; m.parent()
        [0 0]
        [0 0]
        Full MatrixSpace of 2 by 2 sparse matrices over Rational Field
        sage: m = matrix(QQ, 2, 3, {}); m; m.parent()
        [0 0 0]
        [0 0 0]
        Full MatrixSpace of 2 by 3 sparse matrices over Rational Field
        sage: m = matrix(2, {}); m; m.parent()
        [0 0]
        [0 0]
        Full MatrixSpace of 2 by 2 sparse matrices over Integer Ring
        sage: m = matrix(2, 3, {}); m; m.parent()
        [0 0 0]
        [0 0 0]
        Full MatrixSpace of 2 by 3 sparse matrices over Integer Ring

    Matrices with zero rows or columns::

        sage: m = matrix(0,2); m; m.parent()
        []
        Full MatrixSpace of 0 by 2 dense matrices over Integer Ring
        sage: m = matrix(2,0); m; m.parent()
        []
        Full MatrixSpace of 2 by 0 dense matrices over Integer Ring
        sage: m = matrix(0, [1]); m; m.parent()
        Traceback (most recent call last):
        ...
        ValueError: sequence too long (expected length 0, got more)
        sage: m = matrix(1, 0, []); m; m.parent()
        []
        Full MatrixSpace of 1 by 0 dense matrices over Integer Ring
        sage: m = matrix(0, 1, []); m; m.parent()
        []
        Full MatrixSpace of 0 by 1 dense matrices over Integer Ring
        sage: m = matrix(0, []); m; m.parent()
        []
        Full MatrixSpace of 0 by 0 dense matrices over Integer Ring
        sage: m = matrix(0, {}); m; m.parent()
        []
        Full MatrixSpace of 0 by 0 sparse matrices over Integer Ring
        sage: m = matrix(0, {(1,1): 2}); m; m.parent()
        Traceback (most recent call last):
        ...
        IndexError: invalid row index 1
        sage: m = matrix(2, 0, {(1,1): 2}); m; m.parent()
        Traceback (most recent call last):
        ...
        IndexError: invalid column index 1

    Check conversion from numpy::

        sage: # needs numpy
        sage: import numpy
        sage: n = numpy.array([[complex(0,1),complex(0,2)], [3,4]], complex)
        sage: m = matrix(n); m; m.parent()
        [1.0*I 2.0*I]
        [  3.0   4.0]
        Full MatrixSpace of 2 by 2 dense matrices over Complex Double Field
        sage: n = numpy.array([[1,2], [3,4]], 'int32')
        sage: m = matrix(n); m; m.parent()
        [1 2]
        [3 4]
        Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
        sage: n = numpy.array([[1,2,3], [4,5,6], [7,8,9]], 'float32')
        sage: m = matrix(n); m; m.parent()
        [1.0 2.0 3.0]
        [4.0 5.0 6.0]
        [7.0 8.0 9.0]
        Full MatrixSpace of 3 by 3 dense matrices over Real Double Field
        sage: n = numpy.matrix([[1,2,3], [4,5,6], [7,8,9]], 'float64')
        sage: m = matrix(n); m; m.parent()
        [1.0 2.0 3.0]
        [4.0 5.0 6.0]
        [7.0 8.0 9.0]
        Full MatrixSpace of 3 by 3 dense matrices over Real Double Field
        sage: n = numpy.array([[1,2,3], [4,5,6], [7,8,9]], 'complex64')
        sage: m = matrix(n); m; m.parent()
        [1.0 2.0 3.0]
        [4.0 5.0 6.0]
        [7.0 8.0 9.0]
        Full MatrixSpace of 3 by 3 dense matrices over Complex Double Field
        sage: n = numpy.matrix([[1,2,3], [4,5,6], [7,8,9]], 'complex128')
        sage: m = matrix(n); m; m.parent()
        [1.0 2.0 3.0]
        [4.0 5.0 6.0]
        [7.0 8.0 9.0]
        Full MatrixSpace of 3 by 3 dense matrices over Complex Double Field
        sage: a = matrix([[1,2], [3,4]])
        sage: b = matrix(a.numpy()); b
        [1 2]
        [3 4]
        sage: a == b
        True
        sage: c = matrix(a.numpy('float32')); c
        [1.0 2.0]
        [3.0 4.0]
        sage: matrix(numpy.array([[5]]))
        [5]
        sage: matrix(numpy.matrix([[5]]))
        [5]

    A ring and a numpy array::

        sage: # needs numpy
        sage: n = numpy.array([[1,2,3], [4,5,6], [7,8,9]], 'float32')
        sage: m = matrix(ZZ, n); m; m.parent()
        [1 2 3]
        [4 5 6]
        [7 8 9]
        Full MatrixSpace of 3 by 3 dense matrices over Integer Ring
        sage: n = matrix(QQ, 2, 2, [1, 1/2, 1/3, 1/4]).numpy(); n
        array([[1.        , 0.5       ],
               [0.33333333, 0.25      ]])
        sage: matrix(QQ, n)
        [  1 1/2]
        [1/3 1/4]

    The dimensions of a matrix may be given as numpy types::

        sage: matrix(numpy.int32(2), numpy.int32(3))                                    # needs numpy
        [0 0 0]
        [0 0 0]

        sage: matrix(nrows=numpy.int32(2), ncols=numpy.int32(3))                        # needs numpy
        [0 0 0]
        [0 0 0]

    The dimensions of a matrix must have an integral type::

        sage: matrix(RR, 2.0, 2.0)
        Traceback (most recent call last):
        ...
        TypeError: too many arguments in matrix constructor

    More tests::

        sage: v = vector(ZZ, [1, 10, 100])
        sage: m = matrix(ZZ['x'], v); m; m.parent()
        [  1  10 100]
        Full MatrixSpace of 1 by 3 dense matrices
         over Univariate Polynomial Ring in x over Integer Ring
        sage: matrix(ZZ, 10, 10, range(100)).parent()
        Full MatrixSpace of 10 by 10 dense matrices over Integer Ring
        sage: m = matrix(GF(7), [[1/3,2/3,1/2], [3/4,4/5,7]]); m; m.parent()
        [5 3 4]
        [6 5 0]
        Full MatrixSpace of 2 by 3 dense matrices over Finite Field of size 7
        sage: m = matrix([[1,2,3], [RDF(2), CDF(1,2), 3]]); m; m.parent()
        [        1.0         2.0         3.0]
        [        2.0 1.0 + 2.0*I         3.0]
        Full MatrixSpace of 2 by 3 dense matrices over Complex Double Field
        sage: m = matrix(3, 3, 1/2); m; m.parent()
        [1/2   0   0]
        [  0 1/2   0]
        [  0   0 1/2]
        Full MatrixSpace of 3 by 3 dense matrices over Rational Field
        sage: matrix([[1],[2,3]])
        Traceback (most recent call last):
        ...
        ValueError: inconsistent number of columns: should be 1 but got 2
        sage: matrix([[1],2])
        Traceback (most recent call last):
        ...
        TypeError: 'sage.rings.integer.Integer' object is not iterable
        sage: matrix(vector(RR, [1,2,3])).parent()
        Full MatrixSpace of 1 by 3 dense matrices over Real Field with 53 bits of precision

    Check :issue:`10158`::

        sage: matrix(ZZ, [[0] for i in range(10^5)]).is_zero()
        True

    Check :issue:`24459`::

        sage: # needs sage.libs.linbox
        sage: Matrix(ZZ, sys.maxsize, sys.maxsize)
        Traceback (most recent call last):
        ...
        RuntimeError...

    Test a simple ``_matrix_`` method. Note that we are ignoring
    ``base`` which is inefficient but allowed::

        sage: class MatrixTest():
        ....:     def _matrix_(self, base=None):
        ....:         return matrix(ZZ, 2, 2, [1,2,3,4])
        sage: e = MatrixTest()
        sage: matrix(e)
        [1 2]
        [3 4]
        sage: S = MatrixSpace(ZZ, 2, 2)
        sage: M = S(e); M
        [1 2]
        [3 4]
        sage: parent(M) is S
        True
        sage: matrix(RDF, e)
        [1.0 2.0]
        [3.0 4.0]
        sage: S = MatrixSpace(RDF, 2, 2)
        sage: M = S(e); M
        [1.0 2.0]
        [3.0 4.0]
        sage: parent(M) is S
        True

    A redundant ``base_ring`` argument::

        sage: matrix(ZZ, 3, 3, base_ring=ZZ)
        Traceback (most recent call last):
        ...
        TypeError: too many arguments in matrix constructor

    Some calls using an iterator::

        sage: matrix(QQ, 3, 6, range(18), sparse=true)
        [ 0  1  2  3  4  5]
        [ 6  7  8  9 10 11]
        [12 13 14 15 16 17]
        sage: matrix(4, 4, range(16))
        [ 0  1  2  3]
        [ 4  5  6  7]
        [ 8  9 10 11]
        [12 13 14 15]

    AUTHORS:

    - William Stein: Initial implementation

    - Jason Grout (2008-03): almost a complete rewrite, with bits and
      pieces from the original implementation

    - Jeroen Demeyer (2016-02-05): major clean up, see :issue:`20015`
      and :issue:`20016`

    - Jeroen Demeyer (2018-02-20): completely rewritten using
      :class:`MatrixArgs`, see :issue:`24742`
    """
    immutable = kwds.pop('immutable', False)
    return MatrixArgs(*args, **kwds).element(immutable=immutable)


from sage.matrix.special import *


Matrix = matrix


@matrix_method
class options(GlobalOptions):
    r"""
    Global options for matrices.

    @OPTIONS@

    EXAMPLES::

        sage: matrix.options.max_cols = 6
        sage: matrix.options.max_rows = 3
        sage: matrix(ZZ, 3, 6)
        [0 0 0 0 0 0]
        [0 0 0 0 0 0]
        [0 0 0 0 0 0]
        sage: matrix(ZZ, 3, 7)
        3 x 7 dense matrix over Integer Ring...
        sage: matrix(ZZ, 4, 6)
        4 x 6 dense matrix over Integer Ring...
        sage: matrix.options._reset()

    The precision can also be set via the IPython magic::

        sage: from sage.repl.interpreter import get_test_shell
        sage: shell = get_test_shell()
        sage: shell.run_cell('%precision 5')
        '%.5f'
        sage: matrix.options.precision
        5
        sage: A = matrix(RR, [[200/3]]); A
        [66.667]

    The number format can be specified as well::

        sage: matrix.options.format_numeric = '{:.{prec}e}'
        sage: A
        [6.66667e+1]
        sage: matrix.options.format_numeric = '{:.{prec}f}'
        sage: A
        [66.66667]
        sage: matrix.options.format_numeric = '{:+.{prec}g}'
        sage: A
        [+66.667]
        sage: matrix.options._reset()
    """
    NAME = 'Matrix'
    max_cols = dict(default=49,
                    description='maximum number of columns to display',
                    checker=lambda val: val >= 0)
    max_rows = dict(default=19,
                    description='maximum number of rows to display',
                    checker=lambda val: val >= 0)
    precision = \
        dict(default=None,
             description='number of digits to display for floating point '
                         'entries; if ``None``, the exact representation is '
                         'used instead. This option is also set by the '
                         '`IPython magic <https://ipython.readthedocs.io/en/stable/interactive/magics.html#magic-precision>`_ '
                         '``%precision``.',
             checker=lambda val: val is None or val >= 0)
    format_numeric = \
        dict(default='{:.{prec}}',
             description='string used for formatting floating point numbers of'
                         ' an (optional) precision ``prec``; only supported '
                         'for entry types implementing ``__format__``',
             checker=lambda val: isinstance(val.format(3.1415, prec=3), str))
