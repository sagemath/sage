"""
General matrix Constructor
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#       Copyright (C) 2016 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
import types
from .matrix_space import MatrixSpace
from sage.rings.ring import is_Ring
from sage.rings.all import ZZ, RDF, CDF
from sage.arith.srange import srange
from sage.structure.coerce cimport is_numpy_type, py_scalar_parent
from sage.structure.element cimport Vector
from sage.structure.sequence import Sequence
from sage.misc.long cimport pyobject_to_long


class MatrixFactory(object):
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

    The matrix command takes the entries of a matrix, optionally
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
    `c`. Calling matrix() with a Sage object may return something that
    makes sense. Calling matrix() with a NumPy array will convert the
    array to a matrix.

    The ring, number of rows, and number of columns of the matrix can
    be specified by setting the ``ring``, ``nrows``, or ``ncols``
    keyword parameters or by passing them as the first arguments to the
    function in specified order. The ring defaults to ``ZZ`` if it is
    not specified and cannot be determined from the entries. If the
    number of rows and columns are not specified and cannot be
    determined, then an empty 0x0 matrix is returned.

    INPUT:

    -  ``ring`` -- the base ring for the entries of the
       matrix.

    -  ``nrows`` -- the number of rows in the matrix.

    -  ``ncols`` -- the number of columns in the matrix.

    -  ``sparse`` -- create a sparse matrix. This defaults to ``True``
       when the entries are given as a dictionary, otherwise defaults to
       ``False``.

    - ``entries`` -- see examples below.

    OUTPUT:

    a matrix

    EXAMPLES::

        sage: m = matrix(2); m; m.parent()
        [0 0]
        [0 0]
        Full MatrixSpace of 2 by 2 dense matrices over Integer Ring

    ::

        sage: m = matrix(2,3); m; m.parent()
        [0 0 0]
        [0 0 0]
        Full MatrixSpace of 2 by 3 dense matrices over Integer Ring

    ::

        sage: m = matrix(QQ,[[1,2,3],[4,5,6]]); m; m.parent()
        [1 2 3]
        [4 5 6]
        Full MatrixSpace of 2 by 3 dense matrices over Rational Field

    ::

        sage: m = matrix(QQ, 3, 3, lambda i, j: i+j); m
        [0 1 2]
        [1 2 3]
        [2 3 4]
        sage: m = matrix(3, lambda i,j: i-j); m
        [ 0 -1 -2]
        [ 1  0 -1]
        [ 2  1  0]

    ::

        sage: matrix(QQ, 2, 3, lambda x, y: x+y)
        [0 1 2]
        [1 2 3]
        sage: matrix(QQ, 5, 5, lambda x, y: (x+1) / (y+1))
        [  1 1/2 1/3 1/4 1/5]
        [  2   1 2/3 1/2 2/5]
        [  3 3/2   1 3/4 3/5]
        [  4   2 4/3   1 4/5]
        [  5 5/2 5/3 5/4   1]

    ::

        sage: v1=vector((1,2,3))
        sage: v2=vector((4,5,6))
        sage: m = matrix([v1,v2]); m; m.parent()
        [1 2 3]
        [4 5 6]
        Full MatrixSpace of 2 by 3 dense matrices over Integer Ring

    ::

        sage: m = matrix(QQ,2,[1,2,3,4,5,6]); m; m.parent()
        [1 2 3]
        [4 5 6]
        Full MatrixSpace of 2 by 3 dense matrices over Rational Field

    ::

        sage: m = matrix(QQ,2,3,[1,2,3,4,5,6]); m; m.parent()
        [1 2 3]
        [4 5 6]
        Full MatrixSpace of 2 by 3 dense matrices over Rational Field

    ::

        sage: m = matrix({(0,1): 2, (1,1):2/5}); m; m.parent()
        [  0   2]
        [  0 2/5]
        Full MatrixSpace of 2 by 2 sparse matrices over Rational Field

    ::

        sage: m = matrix(QQ,2,3,{(1,1): 2}); m; m.parent()
        [0 0 0]
        [0 2 0]
        Full MatrixSpace of 2 by 3 sparse matrices over Rational Field

    ::

        sage: import numpy
        sage: n = numpy.array([[1,2],[3,4]],float)
        sage: m = matrix(n); m; m.parent()
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

    ::

        sage: g = graphs.PetersenGraph()
        sage: m = matrix(g); m; m.parent()
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
        sage: m = matrix(ring=QQ); m; m.parent()
        []
        Full MatrixSpace of 0 by 0 dense matrices over Rational Field
        sage: m = matrix(0); m; m.parent()
        []
        Full MatrixSpace of 0 by 0 dense matrices over Integer Ring
        sage: m = matrix(0, 0, ring=QQ); m; m.parent()
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
        sage: m = matrix(QQ,2); m; m.parent()
        [0 0]
        [0 0]
        Full MatrixSpace of 2 by 2 dense matrices over Rational Field
        sage: m = matrix(QQ,2,3); m; m.parent()
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
        TypeError: identity matrix must be square
        sage: matrix(QQ,2,3,5)
        Traceback (most recent call last):
        ...
        TypeError: nonzero scalar matrix must be square

    Matrices specified by a list of entries::

        sage: m = matrix([[1,2,3],[4,5,6]]); m; m.parent()
        [1 2 3]
        [4 5 6]
        Full MatrixSpace of 2 by 3 dense matrices over Integer Ring
        sage: m = matrix(QQ,2,[[1,2,3],[4,5,6]]); m; m.parent()
        [1 2 3]
        [4 5 6]
        Full MatrixSpace of 2 by 3 dense matrices over Rational Field
        sage: m = matrix(QQ,3,[[1,2,3],[4,5,6]]); m; m.parent()
        Traceback (most recent call last):
        ...
        ValueError: number of rows does not match up with specified number
        sage: m = matrix(QQ,2,3,[[1,2,3],[4,5,6]]); m; m.parent()
        [1 2 3]
        [4 5 6]
        Full MatrixSpace of 2 by 3 dense matrices over Rational Field
        sage: m = matrix(QQ,2,4,[[1,2,3],[4,5,6]]); m; m.parent()
        Traceback (most recent call last):
        ...
        ValueError: number of columns does not match up with specified number
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
        sage: m = matrix(QQ,[1,2,3,4,5,6]); m; m.parent()
        [1 2 3 4 5 6]
        Full MatrixSpace of 1 by 6 dense matrices over Rational Field
        sage: m = matrix(QQ,3,2,[1,2,3,4,5,6]); m; m.parent()
        [1 2]
        [3 4]
        [5 6]
        Full MatrixSpace of 3 by 2 dense matrices over Rational Field
        sage: m = matrix(QQ,2,4,[1,2,3,4,5,6]); m; m.parent()
        Traceback (most recent call last):
        ...
        ValueError: entries has the wrong length
        sage: m = matrix(QQ,5,[1,2,3,4,5,6]); m; m.parent()
        Traceback (most recent call last):
        ...
        TypeError: cannot construct an element of
        Full MatrixSpace of 5 by 1 dense matrices over Rational Field
        from [1, 2, 3, 4, 5, 6]!

    Matrices specified by a dict of entries::

        sage: m = matrix({(1,1): 2}); m; m.parent()
        [0 0]
        [0 2]
        Full MatrixSpace of 2 by 2 sparse matrices over Integer Ring
        sage: m = matrix({(1,1): 2}, sparse=False); m; m.parent()
        [0 0]
        [0 2]
        Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
        sage: m = matrix(QQ,{(1,1): 2}); m; m.parent()
        [0 0]
        [0 2]
        Full MatrixSpace of 2 by 2 sparse matrices over Rational Field
        sage: m = matrix(QQ,3,{(1,1): 2}); m; m.parent()
        [0 0 0]
        [0 2 0]
        [0 0 0]
        Full MatrixSpace of 3 by 3 sparse matrices over Rational Field
        sage: m = matrix(QQ,3,4,{(1,1): 2}); m; m.parent()
        [0 0 0 0]
        [0 2 0 0]
        [0 0 0 0]
        Full MatrixSpace of 3 by 4 sparse matrices over Rational Field
        sage: m = matrix(QQ,2,{(1,1): 2}); m; m.parent()
        [0 0]
        [0 2]
        Full MatrixSpace of 2 by 2 sparse matrices over Rational Field
        sage: m = matrix(QQ,1,{(1,1): 2}); m; m.parent()
        Traceback (most recent call last):
        ...
        IndexError: invalid entries list
        sage: m = matrix({}); m; m.parent()
        []
        Full MatrixSpace of 0 by 0 sparse matrices over Integer Ring
        sage: m = matrix(QQ,{}); m; m.parent()
        []
        Full MatrixSpace of 0 by 0 sparse matrices over Rational Field
        sage: m = matrix(QQ,2,{}); m; m.parent()
        [0 0]
        [0 0]
        Full MatrixSpace of 2 by 2 sparse matrices over Rational Field
        sage: m = matrix(QQ,2,3,{}); m; m.parent()
        [0 0 0]
        [0 0 0]
        Full MatrixSpace of 2 by 3 sparse matrices over Rational Field
        sage: m = matrix(2,{}); m; m.parent()
        [0 0]
        [0 0]
        Full MatrixSpace of 2 by 2 sparse matrices over Integer Ring
        sage: m = matrix(2,3,{}); m; m.parent()
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
        sage: m = matrix(0,[1]); m; m.parent()
        Traceback (most recent call last):
        ...
        ValueError: entries has the wrong length
        sage: m = matrix(1,0,[]); m; m.parent()
        []
        Full MatrixSpace of 1 by 0 dense matrices over Integer Ring
        sage: m = matrix(0,1,[]); m; m.parent()
        []
        Full MatrixSpace of 0 by 1 dense matrices over Integer Ring
        sage: m = matrix(0,[]); m; m.parent()
        []
        Full MatrixSpace of 0 by 0 dense matrices over Integer Ring
        sage: m = matrix(0,{}); m; m.parent()
        []
        Full MatrixSpace of 0 by 0 sparse matrices over Integer Ring
        sage: m = matrix(0,{(1,1):2}); m; m.parent()
        Traceback (most recent call last):
        ...
        IndexError: invalid entries list
        sage: m = matrix(2,0,{(1,1):2}); m; m.parent()
        Traceback (most recent call last):
        ...
        IndexError: invalid entries list

    Check conversion from numpy::

        sage: import numpy
        sage: n = numpy.array([[numpy.complex(0,1),numpy.complex(0,2)],[3,4]],complex)
        sage: m = matrix(n); m; m.parent()
        [1.0*I 2.0*I]
        [  3.0   4.0]
        Full MatrixSpace of 2 by 2 dense matrices over Complex Double Field
        sage: n = numpy.array([[1,2],[3,4]],'int32')
        sage: m = matrix(n); m; m.parent()
        [1 2]
        [3 4]
        Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
        sage: n = numpy.array([[1,2,3],[4,5,6],[7,8,9]],'float32')
        sage: m = matrix(n); m; m.parent()
        [1.0 2.0 3.0]
        [4.0 5.0 6.0]
        [7.0 8.0 9.0]
        Full MatrixSpace of 3 by 3 dense matrices over Real Double Field
        sage: n = numpy.matrix([[1,2,3],[4,5,6],[7,8,9]],'float64')
        sage: m = matrix(n); m; m.parent()
        [1.0 2.0 3.0]
        [4.0 5.0 6.0]
        [7.0 8.0 9.0]
        Full MatrixSpace of 3 by 3 dense matrices over Real Double Field
        sage: n = numpy.array([[1,2,3],[4,5,6],[7,8,9]],'complex64')
        sage: m = matrix(n); m; m.parent()
        [1.0 2.0 3.0]
        [4.0 5.0 6.0]
        [7.0 8.0 9.0]
        Full MatrixSpace of 3 by 3 dense matrices over Complex Double Field
        sage: n = numpy.matrix([[1,2,3],[4,5,6],[7,8,9]],'complex128')
        sage: m = matrix(n); m; m.parent()
        [1.0 2.0 3.0]
        [4.0 5.0 6.0]
        [7.0 8.0 9.0]
        Full MatrixSpace of 3 by 3 dense matrices over Complex Double Field
        sage: a = matrix([[1,2],[3,4]])
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

        sage: n = numpy.array([[1,2,3],[4,5,6],[7,8,9]],'float32')
        sage: m = matrix(ZZ, n); m; m.parent()
        [1 2 3]
        [4 5 6]
        [7 8 9]
        Full MatrixSpace of 3 by 3 dense matrices over Integer Ring
        sage: n = matrix(QQ, 2, 2, [1, 1/2, 1/3, 1/4]).numpy(); n
        array([[ 1.        ,  0.5       ],
               [ 0.33333333,  0.25      ]])
        sage: matrix(QQ, n)
        [  1 1/2]
        [1/3 1/4]

    The dimensions of a matrix may be given as numpy types::

        sage: matrix(numpy.int32(2), ncols=numpy.int32(3))
        [0 0 0]
        [0 0 0]

    The dimensions of a matrix must have an integral type::

        sage: matrix(RR, 2.0, 2.0)
        Traceback (most recent call last):
        ...
        TypeError: invalid matrix constructor: type matrix? for help

    More tests::

        sage: v = vector(ZZ, [1, 10, 100])
        sage: m = matrix(ZZ['x'], v); m; m.parent()
        [  1  10 100]
        Full MatrixSpace of 1 by 3 dense matrices over Univariate Polynomial Ring in x over Integer Ring
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
        sage: m = matrix(3,3,1/2); m; m.parent()
        [1/2   0   0]
        [  0 1/2   0]
        [  0   0 1/2]
        Full MatrixSpace of 3 by 3 dense matrices over Rational Field
        sage: matrix([[1],[2,3]])
        Traceback (most recent call last):
        ...
        ValueError: list of rows is not valid (rows are wrong types or lengths)
        sage: matrix([[1],2])
        Traceback (most recent call last):
        ...
        ValueError: list of rows is not valid (rows are wrong types or lengths)
        sage: matrix(vector(RR,[1,2,3])).parent()
        Full MatrixSpace of 1 by 3 dense matrices over Real Field with 53 bits of precision

    Check :trac:`10158`::

        sage: matrix(ZZ, [[0] for i in range(10^5)]).is_zero()
        True

    Test conversion using a ``_matrix_`` method::

        sage: A = gap(MatrixSpace(QQ, 2, 2)(range(4)))
        sage: matrix(QQ, A)
        [0 1]
        [2 3]
        sage: matrix(A, ring=QQ)
        [0 1]
        [2 3]
        sage: matrix(A, QQ)
        doctest:...: DeprecationWarning: when constructing a matrix, the ring must be the first argument
        See http://trac.sagemath.org/20015 for details.
        [0 1]
        [2 3]

    A redundant ``ring`` argument::

        sage: matrix(ZZ, 3, 3, ring=ZZ)
        Traceback (most recent call last):
        ...
        TypeError: invalid matrix constructor: type matrix? for help

    TESTS:

    Some calls using an iterator (note that xrange is no longer available
    in Python 3)::

        sage: from six.moves import range
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

    - Jeroen Demeyer (2016-02-05): major clean up, see :trac:`20015`
      and :trac:`20016`
    """
    def __call__(self, *Args, ring=None, nrows=None, ncols=None, sparse=None):
        cdef list args = list(Args)

        # Check for deprecated (matrixable object, ring) argument
        if len(args) == 2 and hasattr(args[0], '_matrix_'):
            from sage.misc.superseded import deprecation
            deprecation(20015, "when constructing a matrix, the ring must be the first argument")
            args = [args[1], args[0]]

        # ring argument
        if ring is None and args and is_Ring(args[0]):
            ring = args.pop(0)

        # object with _matrix_ method
        if args:
            try:
                makematrix = args[0]._matrix_
            except AttributeError:
                pass
            else:
                if ring is None:
                    args.pop(0)
                else:
                    args[0] = ring
                if sparse is None:
                    return makematrix(*args)
                else:
                    return makematrix(*args, sparse=sparse)

        # nrows argument
        if nrows is None and args:
            arg = args[0]
            try:
                if is_numpy_type(type(arg)):
                    import numpy
                    if isinstance(arg, numpy.ndarray):
                        raise TypeError
                nrows = pyobject_to_long(arg)
            except TypeError:
                pass
            else:
                args.pop(0)

        # ncols argument
        if ncols is None and args:
            arg = args[0]
            try:
                if is_numpy_type(type(arg)):
                    import numpy
                    if isinstance(arg, numpy.ndarray):
                        raise TypeError
                ncols = pyobject_to_long(arg)
            except TypeError:
                pass
            else:
                args.pop(0)

        # Now we've taken care of initial ring, nrows, and ncols arguments.
        # We've also taken care of the Sage object case.

        # Now the rest of the arguments are a list of rows, a flat list of
        # entries, a callable, a dict, a numpy array, or a single value.
        entry_ring = ZZ
        if not args:
            # If no entries are specified, pass back a zero matrix
            entries = 0
        elif len(args) == 1:
            arg = args[0]
            if isinstance(arg, (types.FunctionType, types.LambdaType, types.MethodType)):
                if ncols is None and nrows is None:
                    raise TypeError("when passing in a callable, the dimensions of the matrix must be specified")
                if ncols is None:
                    ncols = nrows
                elif nrows is None:
                    nrows = ncols

                irange = srange(nrows)
                jrange = srange(ncols)
                arg = [[arg(i, j) for j in jrange] for i in irange]

            if isinstance(arg, xrange):
                arg = list(arg)
            if isinstance(arg, (list, tuple)):
                if not arg:
                    # no entries are specified, pass back the zero matrix
                    entries = 0
                elif isinstance(arg[0], (list, tuple)) or isinstance(arg[0], Vector):
                    # Ensure we have a list of lists, each inner list having the same number of elements
                    first_len = len(arg[0])
                    if not all( (isinstance(v, (list, tuple)) or isinstance(v, Vector)) and len(v) == first_len for v in arg):
                        raise ValueError("list of rows is not valid (rows are wrong types or lengths)")
                    # We have a list of rows or vectors
                    if nrows is None:
                        nrows = len(arg)
                    elif nrows != len(arg):
                        raise ValueError("number of rows does not match up with specified number")
                    if ncols is None:
                        ncols = len(arg[0])
                    elif ncols != len(arg[0]):
                        raise ValueError("number of columns does not match up with specified number")

                    entries = []
                    for v in arg:
                        entries.extend(v)

                else:
                    # We have a flat list; figure out nrows and ncols
                    if nrows is None:
                        nrows = 1

                    if nrows > 0:
                        if ncols is None:
                            ncols = len(arg) // nrows
                        elif ncols != len(arg) // nrows:
                            raise ValueError("entries has the wrong length")
                    elif len(arg) > 0:
                        raise ValueError("entries has the wrong length")

                    entries = arg

                if nrows > 0 and ncols > 0 and ring is None:
                    entries, ring = prepare(entries)

            elif isinstance(arg, dict):
                # We have a dictionary: default to sparse
                if sparse is None:
                    sparse = True
                if not arg:
                    # no entries are specified, pass back the zero matrix
                    entries = 0
                else:
                    entries, entry_ring = prepare_dict(arg)
                    if nrows is None:
                        nrows = nrows_from_dict(entries)
                        ncols = ncols_from_dict(entries)
                    # note that ncols can still be None if nrows is set --
                    # it will be assigned nrows down below.

                # See the construction after the numpy case below.
            else:
                if is_numpy_type(type(arg)):
                    import numpy
                    if isinstance(arg, numpy.ndarray):
                        # Convert to a numpy array if it was a matrix.
                        if type(arg) is not numpy.ndarray:
                            arg = numpy.array(arg)

                        str_dtype = str(arg.dtype)

                        if not (arg.flags.c_contiguous is True or arg.flags.f_contiguous is True):
                            raise TypeError('numpy matrix must be either c_contiguous or f_contiguous')

                        if str_dtype.count('float32') == 1:
                            m = matrix(RDF, arg.shape[0], arg.shape[1], 0)
                            m._replace_self_with_numpy32(arg)
                        elif str_dtype.count('float64') == 1:
                            m = matrix(RDF, arg.shape[0], arg.shape[1], 0)
                            m._replace_self_with_numpy(arg)
                        elif str_dtype.count('complex64') == 1:
                            m = matrix(CDF, arg.shape[0], arg.shape[1], 0)
                            m._replace_self_with_numpy32(arg)
                        elif str_dtype.count('complex128') == 1:
                            m = matrix(CDF, arg.shape[0], arg.shape[1], 0)
                            m._replace_self_with_numpy(arg)
                        elif str_dtype.count('int') == 1:
                            m = matrix(ZZ, [list(row) for row in list(arg)])
                        elif str_dtype.count('object') == 1:
                            # Get the raw nested list from the numpy array
                            # and feed it back into matrix
                            m = matrix([list(row) for row in list(arg)])
                        else:
                            raise TypeError("cannot convert NumPy matrix to Sage matrix")

                        if ring is not None and m.base_ring() is not ring:
                            m = m.change_ring(ring)

                        return m
                elif nrows is not None and ncols is not None:
                    # assume that we should just pass the thing into the
                    # MatrixSpace constructor and hope for the best
                    # This is not documented, but it is doctested
                    if ring is None:
                        entry_ring = arg.parent()
                    entries = arg
                else:
                    raise TypeError("invalid matrix constructor: type matrix? for help")
        else:  # len(args) >= 2
            raise TypeError("invalid matrix constructor: type matrix? for help")

        if ring is None:
            ring = entry_ring
        if nrows is None:
            nrows = 0
        if ncols is None:
            ncols = nrows

        return MatrixSpace(ring, nrows, ncols, sparse=sparse)(entries)

Matrix = matrix = MatrixFactory()


def prepare(w):
    """
    Given a list w of numbers, find a common ring that they all
    canonically map to, and return the list of images of the elements
    of w in that ring along with the ring.

    This is for internal use by the matrix function.

    INPUT:

    - ``w`` - list

    OUTPUT:

    list, ring

    EXAMPLES::

        sage: sage.matrix.constructor.prepare([-2, Mod(1,7)])
        ([5, 1], Ring of integers modulo 7)

    Notice that the elements must all canonically coerce to a common
    ring (since Sequence is called)::

        sage: sage.matrix.constructor.prepare([2/1, Mod(1,7)])
        Traceback (most recent call last):
        ...
        TypeError: unable to find a common ring for all elements

    TESTS:

    Check that :trac:`19920` is fixed::

        sage: import numpy
        sage: matrix([[numpy.int8(1)]])
        [1]
    """
    if not w:
        return Sequence([], ZZ), ZZ
    entries = Sequence(w)
    ring = entries.universe()
    if isinstance(ring,type):
        ring = py_scalar_parent(ring)
    if not is_Ring(ring):
        raise TypeError("unable to find a common ring for all elements")
    return entries, ring

def prepare_dict(w):
    """
    Given a dictionary w of numbers, find a common ring that they all
    canonically map to, and return the dictionary of images of the
    elements of w in that ring along with the ring.

    This is for internal use by the matrix function.

    INPUT:

    - ``w`` - dict

    OUTPUT:

    dict, ring

    EXAMPLES::

        sage: sage.matrix.constructor.prepare_dict({(0,1):2, (4,10):Mod(1,7)})
        ({(0, 1): 2, (4, 10): 1}, Ring of integers modulo 7)
    """
    Z = w.items()
    X = [x for _, x in Z]
    entries, ring = prepare(X)
    return dict([(Z[i][0],entries[i]) for i in range(len(entries))]), ring


def nrows_from_dict(d):
    """
    Given a dictionary that defines a sparse matrix, return the number
    of rows that matrix should have.

    This is for internal use by the matrix function.

    INPUT:

    - ``d`` - dict

    OUTPUT:

    integer

    EXAMPLES::

        sage: sage.matrix.constructor.nrows_from_dict({})
        0

    Here the answer is 301 not 300, since there is a 0-th row.

    ::

        sage: sage.matrix.constructor.nrows_from_dict({(300,4):10})
        301
    """
    if 0 == len(d):
        return 0
    return max([0] + [ij[0] for ij in d.keys()]) + 1

def ncols_from_dict(d):
    """
    Given a dictionary that defines a sparse matrix, return the number
    of columns that matrix should have.

    This is for internal use by the matrix function.

    INPUT:

    - ``d`` - dict

    OUTPUT:

    integer

    EXAMPLES::

        sage: sage.matrix.constructor.ncols_from_dict({})
        0

    Here the answer is 301 not 300, since there is a 0-th row.

    ::

        sage: sage.matrix.constructor.ncols_from_dict({(4,300):10})
        301
    """
    if 0 == len(d):
        return 0
    return max([0] + [ij[1] for ij in d.keys()]) + 1


from .special import *
