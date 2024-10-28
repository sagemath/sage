"""
Constructors for special matrices

This module gathers several constructors for special, commonly used or
interesting matrices. These can be reached through ``matrix.<tab>``.

For example, here is a circulant matrix of order five::

    sage: matrix.circulant(SR.var('a b c d e'))                                         # needs sage.symbolic
    [a b c d e]
    [e a b c d]
    [d e a b c]
    [c d e a b]
    [b c d e a]

The following constructions are available:

.. csv-table::
    :class: contentstable
    :widths: 30
    :delim: |

    :meth:`~sage.matrix.special.block_diagonal_matrix`
    :meth:`~sage.matrix.special.block_matrix`
    :meth:`~sage.matrix.special.circulant`
    :meth:`~sage.matrix.special.column_matrix`
    :meth:`~sage.matrix.special.companion_matrix`
    :meth:`~sage.matrix.special.diagonal_matrix`
    :meth:`~sage.matrix.special.elementary_matrix`
    :meth:`~sage.matrix.special.hankel`
    :meth:`~sage.matrix.special.hilbert`
    :meth:`~sage.matrix.special.identity_matrix`
    :meth:`~sage.matrix.special.ith_to_zero_rotation_matrix`
    :meth:`~sage.matrix.special.jordan_block`
    :meth:`~sage.matrix.special.lehmer`
    :meth:`~sage.matrix.special.ones_matrix`
    :meth:`~sage.matrix.special.random_matrix`
    :meth:`~sage.matrix.special.random_diagonalizable_matrix`
    :meth:`~sage.matrix.special.random_echelonizable_matrix`
    :meth:`~sage.matrix.special.random_rref_matrix`
    :meth:`~sage.matrix.special.random_subspaces_matrix`
    :meth:`~sage.matrix.special.random_unimodular_matrix`
    :meth:`~sage.matrix.special.toeplitz`
    :meth:`~sage.matrix.special.vandermonde`
    :meth:`~sage.matrix.special.vector_on_axis_rotation_matrix`
    :meth:`~sage.matrix.special.zero_matrix`

The Combinatorics module provides further matrix constructors, such as Hadamard
matrices and Latin squares. See:

    - :mod:`sage.combinat.matrices.hadamard_matrix`
    - :mod:`sage.combinat.matrices.latin`
"""

# ****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from copy import copy

import sage.matrix.matrix_space as matrix_space
from sage.categories.rings import Rings
from sage.modules.free_module_element import vector
from sage.structure.element import Matrix, parent, RingElement
from sage.structure.sequence import Sequence
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.integer import Integer
from sage.misc.misc_c import running_total
from sage.misc.prandom import randint, shuffle
from .constructor import matrix
import sage.categories.pushout


def matrix_method(func=None, name=None):
    """
    Allow a function to be tab-completed on the global matrix
    constructor object.

    INPUT:

    - ``*function`` -- a single argument; the function that is being
      decorated

    - ``**kwds`` -- a single optional keyword argument
      ``name=<string>``. The name of the corresponding method in the
      global matrix constructor object. If not given, it is derived
      from the function name.

    EXAMPLES::

        sage: from sage.matrix.constructor import matrix_method
        sage: def foo_matrix(n): return matrix.diagonal(range(n))
        sage: matrix_method(foo_matrix)
        <function foo_matrix at ...>
        sage: matrix.foo(5)
        [0 0 0 0 0]
        [0 1 0 0 0]
        [0 0 2 0 0]
        [0 0 0 3 0]
        [0 0 0 0 4]
        sage: matrix_method(foo_matrix, name='bar')
        <function foo_matrix at ...>
        sage: matrix.bar(3)
        [0 0 0]
        [0 1 0]
        [0 0 2]
    """
    if func is not None:
        if name is None:
            name = func.__name__.replace('matrix', '').strip('_')
        prefix = "    This function is available as %s(...) and matrix.%s(...)." % (
            func.__name__, name)
        func.__doc__ = "%s\n\n%s" % (prefix, func.__doc__)
        setattr(matrix, name, func)
        return func
    else:
        return lambda func: matrix_method(func, name=name)


@matrix_method
def column_matrix(*args, **kwds):
    r"""
    Construct a matrix, and then swap rows for columns and columns for rows.

    .. NOTE::

        Linear algebra in Sage favors rows over columns.  So,
        generally, when creating a matrix, input vectors and lists are
        treated as rows.  This function is a convenience that turns
        around this convention when creating a matrix.  If you are not
        familiar with the usual :func:`matrix`
        constructor, you might want to consider it first.

    INPUT:

    Inputs are almost exactly the same as for the :func:`matrix`
    constructor, which are documented there.  But see
    examples below for how dimensions are handled.

    OUTPUT:

    Output is exactly the transpose of what the :func:`matrix`
    constructor would return.  In other words, the
    ``matrix`` constructor builds a matrix and then this function
    exchanges rows for columns, and columns for rows.

    EXAMPLES:

    The most compelling use of this function is when you have a
    collection of lists or vectors that you would like to become the
    columns of a matrix. In almost any other situation, the
    :func:`matrix` constructor can probably do the
    job just as easily, or easier. ::

        sage: col_1 = [1,2,3]
        sage: col_2 = [4,5,6]
        sage: column_matrix([col_1, col_2])
        [1 4]
        [2 5]
        [3 6]

        sage: v1 = vector(QQ, [10, 20])
        sage: v2 = vector(QQ, [30, 40])
        sage: column_matrix(QQ, [v1, v2])
        [10 30]
        [20 40]

    If you only specify one dimension along with a flat list of entries,
    then it will be the number of columns in the result (which is different
    from the behavior of the ``matrix`` constructor).  ::

        sage: column_matrix(ZZ, 8, range(24))
        [ 0  3  6  9 12 15 18 21]
        [ 1  4  7 10 13 16 19 22]
        [ 2  5  8 11 14 17 20 23]

    And when you specify two dimensions, then they should be number of
    columns first, then the number of rows, which is the reverse of how
    they would be specified for the ``matrix`` constructor. ::

        sage: column_matrix(QQ, 5, 3, range(15))
        [ 0  3  6  9 12]
        [ 1  4  7 10 13]
        [ 2  5  8 11 14]

    And a few unproductive, but illustrative, examples. ::

        sage: A = matrix(ZZ, 3, 4, range(12))
        sage: B = column_matrix(ZZ, 3, 4, range(12))
        sage: A == B.transpose()
        True

        sage: A = matrix(QQ, 7, 12, range(84))
        sage: A == column_matrix(A.columns())
        True

        sage: A = column_matrix(QQ, matrix(ZZ, 3, 2, range(6)) )
        sage: A
        [0 2 4]
        [1 3 5]
        sage: A.parent()
        Full MatrixSpace of 2 by 3 dense matrices over Rational Field
    """
    return matrix(*args, **kwds).transpose()


@matrix_method
def random_matrix(ring, nrows, ncols=None, algorithm='randomize', implementation=None, *args, **kwds):
    r"""
    Return a random matrix with entries in a specified ring, and possibly with additional properties.

    INPUT:

    - ``ring`` -- base ring for entries of the matrix

    - ``nrows`` -- integer; number of rows

    - ``ncols`` -- (default: ``None``) number of columns. If ``None``
      defaults to ``nrows``.

    - ``algorithm`` -- (default: ``'randomize'``) determines what properties
      the matrix will have.  See examples below for possible additional
      arguments.

      - ``'randomize'`` -- create a matrix of random elements from the
        base ring, possibly controlling the density of nonzero entries

      - ``'echelon_form'`` -- creates a matrix in echelon form

      - ``'echelonizable'`` -- creates a matrix that has a predictable
        echelon form

      - ``'subspaces'`` -- creates a matrix whose four subspaces, when
        explored, have reasonably sized, integral valued, entries

      - ``'unimodular'`` -- creates a matrix of determinant 1

      - ``'diagonalizable'`` -- creates a diagonalizable matrix whose
        eigenvectors, if computed by hand, will have only integer entries

    - ``implementation`` -- (``None`` or string or a matrix class) a possible
      implementation. See the documentation of the constructor of
      :class:`~sage.matrix.matrix_space.MatrixSpace`.

    - ``*args, **kwds`` -- arguments and keywords to describe additional
      properties. See more detailed documentation below

    .. warning::

        Matrices generated are not uniformly distributed. For unimodular
        matrices over finite field this function does not even generate
        all of them: for example ``Matrix.random(GF(3), 2, algorithm='unimodular')``
        never generates ``[[2,0],[0,2]]``. This function is made for
        teaching purposes.

    .. warning::

        An upper bound on the absolute value of the entries may be set
        when the ``algorithm`` is ``echelonizable`` or ``unimodular``.
        In these cases it is possible for this constructor to fail with
        a :exc:`ValueError`.  If you *must* have this routine return
        successfully, do not set ``upper_bound``.  This behavior can
        be partially controlled by a ``max_tries`` keyword.

    .. NOTE::

        When constructing matrices with random entries and no
        additional properties (i.e. when ``algorithm='randomize'``),
        most of the randomness is controlled by the ``random_element``
        method for elements of the base ring of the matrix, so the
        documentation of that method may be relevant or useful.

    EXAMPLES:

    Random integer matrices.  With no arguments, the majority of the entries
    are zero, -1, and 1, and rarely "large." ::

        sage: from collections import defaultdict
        sage: total_count = 0
        sage: dic = defaultdict(Integer)
        sage: def add_samples(*args, **kwds):
        ....:     global dic, total_count
        ....:     for _ in range(100):
        ....:         A = random_matrix(*args, **kwds)
        ....:         for a in A.list():
        ....:             dic[a] += 1
        ....:             total_count += 1.0

        sage: expected = lambda n : 2 / (5*abs(n)*(abs(n) + 1)) if n != 0 else 1/5
        sage: expected(0)
        1/5
        sage: expected(0) == expected(1) == expected(-1)
        True
        sage: expected(100)
        1/25250
        sage: add_samples(ZZ, 5, 5)
        sage: while not all(abs(dic[a]/total_count - expected(a)) < 0.001 for a in dic):
        ....:     add_samples(ZZ, 5, 5)

    The ``distribution`` keyword  set to ``uniform`` will limit values
    between -2 and 2. ::

        sage: expected = lambda n : 1/5 if n in range(-2, 3) else 0
        sage: total_count = 0
        sage: dic = defaultdict(Integer)
        sage: add_samples(ZZ, 5, 5, distribution='uniform')
        sage: while not all(abs(dic[a]/total_count - expected(a)) < 0.001 for a in dic):
        ....:     add_samples(ZZ, 5, 5, distribution='uniform')

    The ``x`` and ``y`` keywords can be used to distribute entries uniformly.
    When both are used ``x`` is the minimum and ``y`` is one greater than
    the maximum. ::

        sage: expected = lambda n : 1/30 if n in range(70, 100) else 0
        sage: total_count = 0
        sage: dic = defaultdict(Integer)
        sage: add_samples(ZZ, 4, 8, x=70, y=100)
        sage: while not all(abs(dic[a]/total_count - expected(a)) < 0.001 for a in dic):
        ....:     add_samples(ZZ, 4, 8, x=70, y=100)

        sage: expected = lambda n : 1/10 if n in range(-5, 5) else 0
        sage: total_count = 0
        sage: dic = defaultdict(Integer)
        sage: add_samples(ZZ, 3, 7, x=-5, y=5)
        sage: while not all(abs(dic[a]/total_count - expected(a)) < 0.001 for a in dic):
        ....:     add_samples(ZZ, 3, 7, x=-5, y=5)

    If only ``x`` is given, then it is used as the upper bound of a range
    starting at 0. ::

        sage: expected = lambda n : 1/25 if n in range(25) else 0
        sage: total_count = 0
        sage: dic = defaultdict(Integer)
        sage: add_samples(ZZ, 5, 5, x=25)
        sage: while not all(abs(dic[a]/total_count - expected(a)) < 0.001 for a in dic):
        ....:     add_samples(ZZ, 5, 5, x=25)

    To control the number of nonzero entries, use the ``density`` keyword
    at a value strictly below the default of 1.0.  The ``density`` keyword
    is used to compute the number of entries per row that will be nonzero, but the
    same entry may be selected more than once.  So the value provided will
    be an upper bound for the density of the created matrix.  Note that for
    a square matrix it is only necessary to set a single dimension. ::

        sage: def add_sample(*args, **kwds):
        ....:     global density_sum, total_count
        ....:     total_count += 1.0
        ....:     A = random_matrix(*args, **kwds)
        ....:     density_sum += float(A.density())

        sage: # needs sage.libs.linbox (otherwise timeout)
        sage: density_sum = 0.0
        sage: total_count = 0.0
        sage: add_sample(ZZ, 5, x=-10, y=10, density=0.75)
        sage: expected_density = (1 - (4/5)^3)
        sage: float(expected_density)
        0.488
        sage: while abs(density_sum/total_count - expected_density) > 0.001:
        ....:     add_sample(ZZ, 5, x=-10, y=10, density=0.75)

        sage: # needs sage.libs.linbox (otherwise timeout)
        sage: density_sum = 0.0
        sage: total_count = 0.0
        sage: add_sample(ZZ, 5, x=20, y=30, density=0.75)
        sage: while abs(density_sum/total_count - expected_density) > 0.001:
        ....:     add_sample(ZZ, 5, x=20, y=30, density=0.75)

        sage: # needs sage.libs.linbox (otherwise timeout)
        sage: density_sum = 0.0
        sage: total_count = 0.0
        sage: add_sample(ZZ, 100, x=20, y=30, density=0.75)
        sage: expected_density = (1 - (99/100)^75)
        sage: float(expected_density)
        0.529...
        sage: while abs(density_sum/total_count - expected_density) > 0.001:
        ....:     add_sample(ZZ, 100, x=20, y=30, density=0.75)

    For a matrix with low density it may be advisable to insist on a sparse
    representation, as this representation is not selected automatically. ::

        sage: A = random_matrix(ZZ, 5, 5)
        sage: A.is_sparse()
        False
        sage: A = random_matrix(ZZ, 5, 5, sparse=True)
        sage: A.is_sparse()
        True

    For algorithm testing you might want to control the number of bits,
    say 10,000 entries, each limited to 16 bits.  ::

        sage: # needs sage.libs.linbox (otherwise timeout)
        sage: A = random_matrix(ZZ, 100, 100, x=2^16); A
        100 x 100 dense matrix over Integer Ring (use the '.str()' method to see the entries)

    One can prescribe a specific matrix implementation::

        sage: K.<a> = FiniteField(2^8)                                                  # needs sage.rings.finite_rings
        sage: type(random_matrix(K, 2, 5))                                              # needs sage.libs.m4ri sage.rings.finite_rings
        <class 'sage.matrix.matrix_gf2e_dense.Matrix_gf2e_dense'>
        sage: type(random_matrix(K, 2, 5, implementation='generic'))                    # needs sage.rings.finite_rings
        <class 'sage.matrix.matrix_generic_dense.Matrix_generic_dense'>

    Random rational matrices.  Now ``num_bound`` and ``den_bound`` control the
    generation of random elements, by specifying limits on the absolute value of
    numerators and denominators (respectively).  Entries will be positive and
    negative (map the absolute value function through the entries to get all
    positive values).  If either the numerator or denominator bound (or both)
    is not used, then the values default to ``2``::

        sage: A = random_matrix(QQ, 2, 8, num_bound=20, den_bound=4)
        sage: A.dimensions()
        (2, 8)
        sage: type(A)
        <class 'sage.matrix.matrix_rational_dense.Matrix_rational_dense'>
        sage: all(a.numerator() in range(-20, 21) and
        ....:     a.denominator() in range(1, 5)
        ....:     for a in A.list())
        True

        sage: A = random_matrix(QQ, 4, density=0.5, sparse=True)
        sage: type(A)
        <class 'sage.matrix.matrix_rational_sparse.Matrix_rational_sparse'>
        sage: A.density() <= 0.5
        True

        sage: A = random_matrix(QQ, 3, 10, num_bound = 99, den_bound = 99)
        sage: positives = list(map(abs, A.list()))
        sage: A1 = matrix(QQ, 3, 10, positives)
        sage: all(abs(A.list()[i]) == A1.list()[i] for i in range(30))
        True
        sage: all(a.numerator() in range(100) and
        ....:     a.denominator() in range(1, 100)
        ....:     for a in A1.list())
        True

        sage: A = random_matrix(QQ, 4, 10, den_bound = 10)
        sage: all(a.numerator() in range(-2, 3) and
        ....:     a.denominator() in range(1, 11)
        ....:     for a in A.list())
        True

        sage: A = random_matrix(QQ, 4, 10)
        sage: all(a.numerator() in range(-2, 3) and
        ....:     a.denominator() in range(1, 3)
        ....:     for a in A.list())
        True

    Random matrices over other rings.  Several classes of matrices have specialized
    ``randomize()`` methods.  You can locate these with the Sage command::

        search_def('randomize')

    The default implementation of :meth:`~sage.matrix.matrix2.randomize` relies
    on the ``random_element()`` method for the base ring.  The ``density`` and
    ``sparse`` keywords behave as described above. Since we have a different
    randomisation when using the optional meataxe package, we have to make sure
    that we use the default implementation in this test::

        sage: K.<a> = FiniteField(3^2)                                                  # needs sage.rings.finite_rings
        sage: A = random_matrix(K, 2, 5, implementation='generic')                      # needs sage.rings.finite_rings
        sage: type(A)
        <class 'sage.matrix.matrix_generic_dense.Matrix_generic_dense'>
        sage: A.base_ring() is K                                                        # needs sage.rings.finite_rings
        True
        sage: TestSuite(A).run()

        sage: A = random_matrix(RR, 3, 4, density=0.66)
        sage: type(A)
        <class 'sage.matrix.matrix_generic_dense.Matrix_generic_dense'>
        sage: A.base_ring() is RR
        True
        sage: TestSuite(A).run()

        sage: A = random_matrix(ComplexField(32), 3, density=0.8, sparse=True)
        sage: A.is_sparse()
        True
        sage: type(A)
        <class 'sage.matrix.matrix_generic_sparse.Matrix_generic_sparse'>
        sage: A.base_ring() is ComplexField(32)
        True
        sage: TestSuite(A).run()

    Random matrices in echelon form.  The ``algorithm='echelon_form'`` keyword,
    along with a requested number of nonzero rows (``num_pivots``) will return
    a random matrix in echelon form.  When the base ring is ``QQ`` the result has integer
    entries.  Other exact rings may be also specified. ::

        sage: A = random_matrix(QQ, 4, 8, algorithm='echelon_form', num_pivots=3)
        sage: A.base_ring()
        Rational Field
        sage: (A.nrows(), A.ncols())
        (4, 8)
        sage: A in sage.matrix.matrix_space.MatrixSpace(ZZ, 4, 8)
        True
        sage: A.rank()
        3
        sage: A == A.rref()
        True

    For more, see the documentation of the :func:`~sage.matrix.constructor.random_rref_matrix`
    function.  In the notebook or at the Sage command-line, first execute the following to make
    this further documentation available::

        from sage.matrix.constructor import random_rref_matrix

    Random matrices with predictable echelon forms.  The ``algorithm='echelonizable'``
    keyword, along with a requested rank (``rank``) and optional size control
    (``upper_bound``) will return a random matrix in echelon form.  When the
    base ring is ``ZZ`` or ``QQ`` the result has integer entries, whose magnitudes
    can be limited by the value of ``upper_bound``, and the echelon form of the
    matrix also has integer entries.  Other exact rings may be also
    specified, but there is no notion of controlling the size.  Square matrices
    of full rank generated by this function always have determinant one, and
    can be constructed with the ``unimodular`` keyword. ::

        sage: A = random_matrix(QQ, 4, 8, algorithm='echelonizable', rank=3, upper_bound=60)
        sage: A.base_ring()
        Rational Field
        sage: (A.nrows(), A.ncols())
        (4, 8)
        sage: A in sage.matrix.matrix_space.MatrixSpace(ZZ, 4, 8)
        True
        sage: A.rank()
        3
        sage: all(abs(x)<60 for x in A.list())
        True
        sage: A.rref() in sage.matrix.matrix_space.MatrixSpace(ZZ, 4, 8)
        True

    For more, see the documentation of the :func:`~sage.matrix.constructor.random_echelonizable_matrix`
    function.  In the notebook or at the Sage command-line, first execute the following to make
    this further documentation available::

        from sage.matrix.constructor import random_echelonizable_matrix

    Random diagonalizable matrices.  The ``algorithm='diagonalizable'`` keyword,
    along with a requested matrix size (``size``) and optional lists of
    eigenvalues (``eigenvalues``) and the corresponding eigenspace
    dimensions (``dimensions``) will return a random diagonalizable matrix.
    When the eigenvalues and dimensions are not specified the result will have
    randomly generated values for both that fit with the designated size. ::

        sage: A = random_matrix(QQ, 5, algorithm='diagonalizable',  # random
        ....:                   eigenvalues=[2,3,-1], dimensions=[1,2,2]); A
        sage: all(x in ZZ for x in (A - (2*identity_matrix(5))).rref().list())
        True
        sage: all(x in ZZ for x in (A - 3*identity_matrix(5)).rref().list())
        True
        sage: all(x in ZZ for x in (A - (-1)*identity_matrix(5)).rref().list())
        True
        sage: A.jordan_form()                                                           # needs sage.combinat sage.libs.pari
        [ 2| 0| 0| 0| 0]
        [--+--+--+--+--]
        [ 0| 3| 0| 0| 0]
        [--+--+--+--+--]
        [ 0| 0| 3| 0| 0]
        [--+--+--+--+--]
        [ 0| 0| 0|-1| 0]
        [--+--+--+--+--]
        [ 0| 0| 0| 0|-1]

    For more, see the documentation of the :func:`~sage.matrix.constructor.random_diagonalizable_matrix`
    function.  In the notebook or at the Sage command-line, first execute the following to make
    this further documentation available::

        from sage.matrix.constructor import random_diagonalizable_matrix

    Random matrices with predictable subspaces.  The ``algorithm='subspaces'``
    keyword, along with an optional rank (``rank``) will return
    a matrix whose natural basis vectors for its four fundamental subspaces, if computed as
    described in the documentation of the :func:`~sage.matrix.constructor.random_subspaces_matrix`
    contain only integer entries.  If ``rank``, is not set, the
    rank of the matrix will be generated randomly. ::

        sage: B = random_matrix(QQ, 5, 6, algorithm='subspaces', rank=3); B  # random
        sage: B_expanded = B.augment(identity_matrix(5)).rref()
        sage: (B.nrows(), B.ncols())
        (5, 6)
        sage: all(x in ZZ for x in B_expanded.list())
        True
        sage: C = B_expanded.submatrix(0, 0, B.nrows() - B.nullity(), B.ncols())
        sage: L = B_expanded.submatrix(B.nrows() - B.nullity(), B.ncols())
        sage: B.right_kernel() == C.right_kernel()
        True
        sage: B.row_space() == C.row_space()
        True
        sage: B.column_space() == L.right_kernel()
        True
        sage: B.left_kernel() == L.row_space()
        True

    For more, see the documentation of the :func:`~sage.matrix.constructor.random_subspaces_matrix`
    function.  In the notebook or at the Sage command-line, first execute the following to make
    this further documentation available::

        from sage.matrix.constructor import random_subspaces_matrix

    Random unimodular matrices.  The ``algorithm='unimodular'``
    keyword, along with an optional entry size control (``upper_bound``)
    will return a matrix of determinant 1. When the base ring is ``ZZ``
    or ``QQ`` the result has integer entries, whose magnitudes
    can be limited by the value of ``upper_bound``. ::

        sage: C = random_matrix(QQ, 5, algorithm='unimodular', upper_bound=70); C  # random
        sage: det(C)
        1
        sage: C.base_ring()
        Rational Field
        sage: (C.nrows(), C.ncols())
        (5, 5)
        sage: all(abs(x)<70 for x in C.list())
        True

    For more, see the documentation of the :func:`~sage.matrix.constructor.random_unimodular_matrix`
    function.  In the notebook or at the Sage command-line, first execute the following to make
    this further documentation available::

        from sage.matrix.constructor import random_unimodular_matrix

    TESTS:

    We return an error for a bogus value of ``algorithm``::

        sage: random_matrix(ZZ, 5, algorithm='bogus')
        Traceback (most recent call last):
        ...
        ValueError: random matrix algorithm "bogus" is not recognized

    AUTHOR:

    - William Stein (2007-02-06)

    - Rob Beezer (2010-08-25) Documentation, code to allow additional types of output
    """
    if ncols is None:
        ncols = nrows
    sparse = kwds.pop('sparse', False)
    # Construct the parent of the desired matrix
    parent = matrix_space.MatrixSpace(ring, nrows, ncols, sparse=sparse, implementation=implementation)
    if algorithm == 'randomize':
        density = kwds.pop('density', None)
        # zero matrix is immutable, copy is mutable
        A = copy(parent.zero_matrix())
        if density is None:
            A.randomize(density=float(1), nonzero=False, *args, **kwds)
        else:
            A.randomize(density=density, nonzero=True, *args, **kwds)
        return A
    elif algorithm == 'echelon_form':
        return random_rref_matrix(parent, *args, **kwds)
    elif algorithm == 'echelonizable':
        return random_echelonizable_matrix(parent, *args, **kwds)
    elif algorithm == 'diagonalizable':
        return random_diagonalizable_matrix(parent, *args, **kwds)
    elif algorithm == 'subspaces':
        return random_subspaces_matrix(parent, *args, **kwds)
    elif algorithm == 'unimodular':
        return random_unimodular_matrix(parent, *args, **kwds)
    else:
        raise ValueError('random matrix algorithm "%s" is not recognized' % algorithm)


@matrix_method
def diagonal_matrix(arg0=None, arg1=None, arg2=None, sparse=True):
    r"""
    Return a square matrix with specified diagonal entries, and zeros elsewhere.

    FORMATS:

      1. diagonal_matrix(entries)

      2. diagonal_matrix(nrows, entries)

      3. diagonal_matrix(ring, entries)

      4. diagonal_matrix(ring, nrows, entries)

    INPUT:

    - ``entries`` -- the values to place along the diagonal
      of the returned matrix.  This may be a flat list, a
      flat tuple, a vector or free module element, or
      a one-dimensional NumPy array.

    - ``nrows`` -- the size of the returned matrix, which
      will have an equal number of columns

    - ``ring`` -- the ring containing the entries of the
      diagonal entries.  This may not be specified in
      combination with a NumPy array.

    - ``sparse`` -- boolean (default: ``True``); whether or not
      the result has a sparse implementation

    OUTPUT:

    A square matrix over the given ``ring`` with a size
    given by ``nrows``.  If the ring is not given it
    is inferred from the given entries.  The values on
    the diagonal of the returned matrix come from ``entries``.
    If the number of entries is not enough to fill the whole
    diagonal, it is padded with zeros.

    EXAMPLES:

    We first demonstrate each of the input formats with various
    different ways to specify the entries.

    Format 1: a flat list of entries.  ::

        sage: A = diagonal_matrix([2, 1.3, 5]); A
        [ 2.00000000000000 0.000000000000000 0.000000000000000]
        [0.000000000000000  1.30000000000000 0.000000000000000]
        [0.000000000000000 0.000000000000000  5.00000000000000]
        sage: A.parent()
        Full MatrixSpace of 3 by 3 sparse matrices over Real Field with 53 bits of precision

    Format 2: size specified, a tuple with initial entries. Note that a short list of entries
    is effectively padded with zeros.  ::

        sage: A = diagonal_matrix(3, (4, 5)); A
        [4 0 0]
        [0 5 0]
        [0 0 0]
        sage: A.parent()
        Full MatrixSpace of 3 by 3 sparse matrices over Integer Ring

    Format 3: ring specified, a vector of entries. ::

        sage: A = diagonal_matrix(QQ, vector(ZZ, [1,2,3])); A
        [1 0 0]
        [0 2 0]
        [0 0 3]
        sage: A.parent()
        Full MatrixSpace of 3 by 3 sparse matrices over Rational Field

    Format 4: ring, size and list of entries. ::

        sage: A = diagonal_matrix(FiniteField(3), 3, [2, 16]); A
        [2 0 0]
        [0 1 0]
        [0 0 0]
        sage: A.parent()
        Full MatrixSpace of 3 by 3 sparse matrices over Finite Field of size 3

    NumPy arrays may be used as input. ::

        sage: # needs numpy
        sage: import numpy
        sage: entries = numpy.array([1.2, 5.6]); entries
        array([1.2, 5.6])
        sage: A = diagonal_matrix(3, entries); A
        [1.2 0.0 0.0]
        [0.0 5.6 0.0]
        [0.0 0.0 0.0]
        sage: A.parent()
        Full MatrixSpace of 3 by 3 sparse matrices over Real Double Field

        sage: # needs numpy
        sage: j = complex(0,1)
        sage: entries = numpy.array([2.0+j, 8.1, 3.4+2.6*j]); entries
        array([2. +1.j , 8.1+0.j , 3.4+2.6j])
        sage: A = diagonal_matrix(entries); A
        [2.0 + 1.0*I         0.0         0.0]
        [        0.0         8.1         0.0]
        [        0.0         0.0 3.4 + 2.6*I]
        sage: A.parent()
        Full MatrixSpace of 3 by 3 sparse matrices over Complex Double Field

        sage: # needs numpy
        sage: entries = numpy.array([4, 5, 6])
        sage: A = diagonal_matrix(entries); A
        [4 0 0]
        [0 5 0]
        [0 0 6]
        sage: A.parent()
        Full MatrixSpace of 3 by 3 sparse matrices over Integer Ring

        sage: entries = numpy.array([4.1, 5.2, 6.3])                                    # needs numpy
        sage: A = diagonal_matrix(ZZ, entries); A                                       # needs numpy
        Traceback (most recent call last):
        ...
        TypeError: unable to convert 4.1 to an element of Integer Ring

    By default returned matrices have a sparse implementation.  This can be changed
    when using any of the formats.  ::

        sage: A = diagonal_matrix([1,2,3], sparse=False)
        sage: A.parent()
        Full MatrixSpace of 3 by 3 dense matrices over Integer Ring

    An empty list and no ring specified defaults to the integers. ::

        sage: A = diagonal_matrix([])
        sage: A.parent()
        Full MatrixSpace of 0 by 0 sparse matrices over Integer Ring

    Giving the entries improperly may first complain about not being iterable::

        sage: diagonal_matrix(QQ, 5, 10)
        Traceback (most recent call last):
        ...
        TypeError: 'sage.rings.integer.Integer' object is not iterable

    Giving too many entries will raise an error. ::

        sage: diagonal_matrix(QQ, 3, [1,2,3,4])
        Traceback (most recent call last):
        ...
        ValueError: number of diagonal matrix entries (4) exceeds the requested matrix size (3)

    A negative size sometimes causes the error that there are too many elements. ::

        sage: diagonal_matrix(-2, [2])
        Traceback (most recent call last):
        ...
        ValueError: number of diagonal matrix entries (1) exceeds the requested matrix size (-2)

    Types for the entries need to be iterable (tuple, list, vector, NumPy array,
    etc)::

        sage: diagonal_matrix(x^2)                                                      # needs sage.symbolic
        Traceback (most recent call last):
        ...
        TypeError: 'sage.symbolic.expression.Expression' object is not iterable

    TESTS::

        sage: A = diagonal_matrix(reversed(range(4)))

    AUTHOR:

    - Rob Beezer (2011-01-11): total rewrite
    """
    # Roll arguments leftward
    #
    # Leads with a ring?
    # Formats 3, 4, else remains None
    ring = None
    if arg0 in Rings():
        ring = arg0
        arg0 = arg1
        arg1 = arg2
    # Size of matrix specified?
    # Formats 2, 4
    nrows = None
    if isinstance(arg0, (Integer, int)):
        nrows = arg0
        arg0 = arg1
    # Object holding entries
    # Formats 1, 2, 3, 4
    entries = arg0

    # sanity check for entries
    types = (list, tuple)
    try:
        from numpy import ndarray
    except ImportError:
        pass
    else:
        types += (ndarray,)
    if not isinstance(entries, types):
        entries = list(entries)

    # Reconcile matrix size and number of entries
    try:
        nentries = len(entries)
    except TypeError:
        raise TypeError('unable to determine number of entries for diagonal matrix construction')
    # sometimes catches a negative size
    if nrows is not None and nentries > nrows:
        raise ValueError('number of diagonal matrix entries (%s) exceeds the requested matrix size (%s)' % (nentries, nrows))
    if nrows is None:
        nrows = nentries

    # provide a default ring for an empty list
    if not nentries and ring is None:
        ring = ZZ

    # Convert entries to a list v over a common ring
    from sage.modules.free_module_element import prepare
    v, ring = prepare(entries, ring)

    # Create a "diagonal" dictionary for matrix constructor
    # If nentries < nrows, diagonal is effectively padded with zeros at end
    w = {(i, i): v[i] for i in range(len(v))}

    # Ship ring, matrix size, dictionary to matrix constructor
    if ring is None:
        return matrix(nrows, nrows, w, sparse=sparse)
    else:
        return matrix(ring, nrows, nrows, w, sparse=sparse)


@matrix_method
def identity_matrix(ring, n=0, sparse=False):
    r"""
    Return the `n \times n` identity matrix over the given
    ring.

    The default ring is the integers.

    EXAMPLES::

        sage: M = identity_matrix(QQ, 2); M
        [1 0]
        [0 1]
        sage: M.parent()
        Full MatrixSpace of 2 by 2 dense matrices over Rational Field
        sage: M = identity_matrix(2); M
        [1 0]
        [0 1]
        sage: M.parent()
        Full MatrixSpace of 2 by 2 dense matrices over Integer Ring
        sage: M.is_mutable()
        True
        sage: M = identity_matrix(3, sparse=True); M
        [1 0 0]
        [0 1 0]
        [0 0 1]
        sage: M.parent()
        Full MatrixSpace of 3 by 3 sparse matrices over Integer Ring
        sage: M.is_mutable()
        True
    """
    if isinstance(ring, (Integer, int)):
        n = ring
        ring = ZZ
    return matrix_space.MatrixSpace(ring, n, n, sparse)(1)


@matrix_method
def lehmer(ring, n=0):
    r"""
    Return the `n \times n` Lehmer matrix.

    The default ring is the rationals.

    Element `(i, j)` in the Lehmer matrix is
    `min(i, j)/max(i, j)`.

    See :wikipedia:`Lehmer_matrix`.

    EXAMPLES::

        sage: matrix.lehmer(3)
        [  1 1/2 1/3]
        [1/2   1 2/3]
        [1/3 2/3   1]
    """
    from sage.sets.integer_range import IntegerRange

    if isinstance(ring, (Integer, int)):
        n = ring
        ring = QQ
    return matrix_space.MatrixSpace(ring, n, n).matrix([[min(i, j)/max(i, j) for i in IntegerRange(1, n+1)] for j in IntegerRange(1, n+1)])


@matrix_method
def zero_matrix(ring, nrows=None, ncols=None, sparse=False):
    r"""
    Return the `nrows \times ncols` zero matrix over the given
    ring.

    The default ring is the integers.

    EXAMPLES::

        sage: M = zero_matrix(QQ, 2); M
        [0 0]
        [0 0]
        sage: M.parent()
        Full MatrixSpace of 2 by 2 dense matrices over Rational Field
        sage: M = zero_matrix(2, 3); M
        [0 0 0]
        [0 0 0]
        sage: M.parent()
        Full MatrixSpace of 2 by 3 dense matrices over Integer Ring
        sage: M.is_mutable()
        True
        sage: M = zero_matrix(3, 1, sparse=True); M
        [0]
        [0]
        [0]
        sage: M.parent()
        Full MatrixSpace of 3 by 1 sparse matrices over Integer Ring
        sage: M.is_mutable()
        True
        sage: matrix.zero(5)
        [0 0 0 0 0]
        [0 0 0 0 0]
        [0 0 0 0 0]
        [0 0 0 0 0]
        [0 0 0 0 0]

    TESTS:

    Check that :issue:`38221` is fixed::

        sage: # needs sage.groups
        sage: G = CyclicPermutationGroup(7)
        sage: R = GF(2)
        sage: A = G.algebra(R)
        sage: zero_matrix(A, 3, 3)
        [0 0 0]
        [0 0 0]
        [0 0 0]
    """
    if isinstance(ring, (Integer, int)):
        nrows, ncols = (ring, nrows)
        ring = ZZ
    return matrix_space.MatrixSpace(ring, nrows, ncols, sparse).matrix(None)


@matrix_method
def ones_matrix(ring, nrows=None, ncols=None, sparse=False):
    r"""
    Return a matrix with all entries equal to 1.

    CALL FORMATS:

    In each case, the optional keyword ``sparse`` can be used.

      1. ones_matrix(ring, nrows, ncols)
      2. ones_matrix(ring, nrows)
      3. ones_matrix(nrows, ncols)
      4. ones_matrix(nrows)

    INPUT:

    - ``ring`` -- (default: ``ZZ``) base ring for the matrix
    - ``nrows`` -- number of rows in the matrix
    - ``ncols`` -- number of columns in the matrix;
      if omitted, defaults to the number of rows, producing a square matrix
    - ``sparse`` -- (default: ``False``) if ``True`` creates a sparse representation

    OUTPUT:

    A matrix of size ``nrows`` by ``ncols`` over the ``ring`` with every
    entry equal to 1.  While the result is far from sparse, you may wish
    to choose a sparse representation when mixing this matrix with
    other sparse matrices.

    EXAMPLES:

    A call specifying the ring and the size.  ::

        sage: M= ones_matrix(QQ, 2, 5); M
        [1 1 1 1 1]
        [1 1 1 1 1]
        sage: M.parent()
        Full MatrixSpace of 2 by 5 dense matrices over Rational Field

    Without specifying the number of columns, the result is square. ::

        sage: M = ones_matrix(RR, 2); M
        [1.00000000000000 1.00000000000000]
        [1.00000000000000 1.00000000000000]
        sage: M.parent()
        Full MatrixSpace of 2 by 2 dense matrices over Real Field with 53 bits of precision

    The ring defaults to the integers if not given. ::

        sage: M = ones_matrix(2, 3); M
        [1 1 1]
        [1 1 1]
        sage: M.parent()
        Full MatrixSpace of 2 by 3 dense matrices over Integer Ring

    A lone integer input produces a square matrix over the integers. ::

        sage: M = ones_matrix(3); M
        [1 1 1]
        [1 1 1]
        [1 1 1]
        sage: M.parent()
        Full MatrixSpace of 3 by 3 dense matrices over Integer Ring

    The result can have a sparse implementation. ::

        sage: M = ones_matrix(3, 1, sparse=True); M
        [1]
        [1]
        [1]
        sage: M.parent()
        Full MatrixSpace of 3 by 1 sparse matrices over Integer Ring

    Giving just a ring will yield an error. ::

        sage: ones_matrix(CC)
        Traceback (most recent call last):
        ...
        ValueError: constructing an all ones matrix requires at least one dimension
    """
    if isinstance(ring, (Integer, int)):
        nrows, ncols = (ring, nrows)
        ring = ZZ
    if nrows is None:
        raise ValueError("constructing an all ones matrix requires at least one dimension")
    if ncols is None:
        nents = nrows**2
    else:
        nents = nrows*ncols
    one = ring(1)
    return matrix_space.MatrixSpace(ring, nrows, ncols, sparse).matrix([one]*nents)


@matrix_method
def elementary_matrix(arg0, arg1=None, **kwds):
    r"""
    Create a square matrix that corresponds to a row operation or a column operation.

    FORMATS:

    In each case, ``R`` is the base ring, and is optional. ``n`` is the size
    of the square matrix created.  Any call may include the ``sparse`` keyword
    to determine the representation used.  The default is ``False`` which
    leads to a dense representation.  We describe the matrices by their
    associated row operation, see the output description for more.

    - ``elementary_matrix(R, n, row1=i, row2=j)``

      The matrix which swaps rows ``i`` and ``j``.

    - ``elementary_matrix(R, n, row1=i, scale=s)``

      The matrix which multiplies row ``i`` by ``s``.

    - ``elementary_matrix(R, n, row1=i, row2=j, scale=s)``

      The matrix which multiplies row ``j`` by ``s``
      and adds it to row ``i``.

    Elementary matrices representing column operations are created
    in an entirely analogous way, replacing ``row1`` by ``col1`` and
    replacing ``row2`` by ``col2``.

    Specifying the ring for entries of the matrix is optional.  If it
    is not given, and a scale parameter is provided, then a ring containing
    the value of ``scale`` will be used.  Otherwise, the ring defaults
    to the integers.

    OUTPUT:

    An elementary matrix is a square matrix that is very close to being
    an identity matrix.  If ``E`` is an elementary matrix and ``A`` is any
    matrix with the same number of rows, then ``E*A`` is the result of
    applying a row operation to ``A``.  This is how the three types
    created by this function are described.  Similarly, an elementary matrix
    can be associated with a column operation, so if ``E`` has the same number
    of columns as ``A`` then ``A*E`` is the result of performing a column
    operation on ``A``.

    An elementary matrix representing a row operation is created if ``row1``
    is specified, while an elementary matrix representing a column operation
    is created if ``col1`` is specified.

    EXAMPLES:

    Over the integers, creating row operations. Recall that row
    and column numbering begins at zero.  ::

        sage: A = matrix(ZZ, 4, 10, range(40)); A
        [ 0  1  2  3  4  5  6  7  8  9]
        [10 11 12 13 14 15 16 17 18 19]
        [20 21 22 23 24 25 26 27 28 29]
        [30 31 32 33 34 35 36 37 38 39]

        sage: E = elementary_matrix(4, row1=1, row2=3); E
        [1 0 0 0]
        [0 0 0 1]
        [0 0 1 0]
        [0 1 0 0]
        sage: E*A
        [ 0  1  2  3  4  5  6  7  8  9]
        [30 31 32 33 34 35 36 37 38 39]
        [20 21 22 23 24 25 26 27 28 29]
        [10 11 12 13 14 15 16 17 18 19]

        sage: E = elementary_matrix(4, row1=2, scale=10); E
        [ 1  0  0  0]
        [ 0  1  0  0]
        [ 0  0 10  0]
        [ 0  0  0  1]
        sage: E*A
        [  0   1   2   3   4   5   6   7   8   9]
        [ 10  11  12  13  14  15  16  17  18  19]
        [200 210 220 230 240 250 260 270 280 290]
        [ 30  31  32  33  34  35  36  37  38  39]

        sage: E = elementary_matrix(4, row1=2, row2=1, scale=10); E
        [ 1  0  0  0]
        [ 0  1  0  0]
        [ 0 10  1  0]
        [ 0  0  0  1]
        sage: E*A
        [  0   1   2   3   4   5   6   7   8   9]
        [ 10  11  12  13  14  15  16  17  18  19]
        [120 131 142 153 164 175 186 197 208 219]
        [ 30  31  32  33  34  35  36  37  38  39]

    Over the rationals, now as column operations. Recall that row
    and column numbering begins at zero.  Checks now have the
    elementary matrix on the right.  ::

        sage: A = matrix(QQ, 5, 4, range(20)); A
        [ 0  1  2  3]
        [ 4  5  6  7]
        [ 8  9 10 11]
        [12 13 14 15]
        [16 17 18 19]

        sage: E = elementary_matrix(QQ, 4, col1=1, col2=3); E
        [1 0 0 0]
        [0 0 0 1]
        [0 0 1 0]
        [0 1 0 0]
        sage: A*E
        [ 0  3  2  1]
        [ 4  7  6  5]
        [ 8 11 10  9]
        [12 15 14 13]
        [16 19 18 17]

        sage: E = elementary_matrix(QQ, 4, col1=2, scale=1/2); E
        [  1   0   0   0]
        [  0   1   0   0]
        [  0   0 1/2   0]
        [  0   0   0   1]
        sage: A*E
        [ 0  1  1  3]
        [ 4  5  3  7]
        [ 8  9  5 11]
        [12 13  7 15]
        [16 17  9 19]

        sage: E = elementary_matrix(QQ, 4, col1=2, col2=1, scale=10); E
        [ 1  0  0  0]
        [ 0  1 10  0]
        [ 0  0  1  0]
        [ 0  0  0  1]
        sage: A*E
        [  0   1  12   3]
        [  4   5  56   7]
        [  8   9 100  11]
        [ 12  13 144  15]
        [ 16  17 188  19]

    An elementary matrix is always nonsingular.  Then repeated row
    operations can be represented by products of elementary matrices,
    and this product is again nonsingular.  If row operations are to
    preserve fundamental properties of a matrix (like rank), we do not
    allow scaling a row by zero.  Similarly, the corresponding elementary
    matrix is not constructed.  Also, we do not allow adding a multiple
    of a row to itself, since this could also lead to a new zero row.  ::

        sage: A = matrix(QQ, 4, 10, range(40)); A
        [ 0  1  2  3  4  5  6  7  8  9]
        [10 11 12 13 14 15 16 17 18 19]
        [20 21 22 23 24 25 26 27 28 29]
        [30 31 32 33 34 35 36 37 38 39]

        sage: E1 = elementary_matrix(QQ, 4, row1=0, row2=1)
        sage: E2 = elementary_matrix(QQ, 4, row1=3, row2=0, scale=100)
        sage: E = E2*E1
        sage: E.is_singular()
        False
        sage: E*A
        [  10   11   12   13   14   15   16   17   18   19]
        [   0    1    2    3    4    5    6    7    8    9]
        [  20   21   22   23   24   25   26   27   28   29]
        [1030 1131 1232 1333 1434 1535 1636 1737 1838 1939]

        sage: E3 = elementary_matrix(QQ, 4, row1=3, scale=0)
        Traceback (most recent call last):
        ...
        ValueError: scale parameter of row of elementary matrix must be nonzero

        sage: E4 = elementary_matrix(QQ, 4, row1=3, row2=3, scale=12)
        Traceback (most recent call last):
        ...
        ValueError: cannot add a multiple of a row to itself

    If the ring is not specified, and a scale parameter is given, the
    base ring for the matrix is chosen to contain the scale parameter.
    Otherwise, if no ring is given, the default is the integers. ::

        sage: E = elementary_matrix(4, row1=1, row2=3)
        sage: E.parent()
        Full MatrixSpace of 4 by 4 dense matrices over Integer Ring

        sage: E = elementary_matrix(4, row1=1, scale=4/3)
        sage: E.parent()
        Full MatrixSpace of 4 by 4 dense matrices over Rational Field

        sage: # needs sage.symbolic
        sage: E = elementary_matrix(4, row1=1, scale=I)
        sage: E.parent()
        Full MatrixSpace of 4 by 4 dense matrices over
         Number Field in I with defining polynomial x^2 + 1 with I = 1*I

        sage: # needs sage.rings.complex_double sage.symbolic
        sage: E = elementary_matrix(4, row1=1, scale=CDF(I))
        sage: E.parent()
        Full MatrixSpace of 4 by 4 dense matrices over Complex Double Field

        sage: # needs sage.rings.number_field sage.symbolic
        sage: E = elementary_matrix(4, row1=1, scale=QQbar(I))
        sage: E.parent()
        Full MatrixSpace of 4 by 4 dense matrices over Algebraic Field

    Returned matrices have a dense implementation by default,
    but a sparse implementation may be requested.  ::

        sage: E = elementary_matrix(4, row1=0, row2=1)
        sage: E.is_dense()
        True

        sage: E = elementary_matrix(4, row1=0, row2=1, sparse=True)
        sage: E.is_sparse()
        True

    And the ridiculously small cases.  The zero-row matrix cannot be built
    since then there are no rows to manipulate. ::

        sage: elementary_matrix(QQ, 1, row1=0, row2=0)
        [1]
        sage: elementary_matrix(QQ, 0, row1=0, row2=0)
        Traceback (most recent call last):
        ...
        ValueError: size of elementary matrix must be 1 or greater, not 0

    TESTS::

        sage: E = elementary_matrix('junk', 5, row1=3, row2=1, scale=12)
        Traceback (most recent call last):
        ...
        TypeError: optional first parameter must be a ring, not junk

        sage: E = elementary_matrix(5, row1=3, scale='junk')
        Traceback (most recent call last):
        ...
        TypeError: scale must be an element of some ring, not junk

        sage: E = elementary_matrix(ZZ, 5, row1=3, col2=3, scale=12)
        Traceback (most recent call last):
        ...
        ValueError: received an unexpected keyword: col2=3

        sage: E = elementary_matrix(QQ, row1=3, scale=12)
        Traceback (most recent call last):
        ...
        ValueError: size of elementary matrix must be given

        sage: E = elementary_matrix(ZZ, 4/3, row1=3, row2=1, scale=12)
        Traceback (most recent call last):
        ...
        TypeError: size of elementary matrix must be an integer, not 4/3

        sage: E = elementary_matrix(ZZ, -3, row1=3, row2=1, scale=12)
        Traceback (most recent call last):
        ...
        ValueError: size of elementary matrix must be 1 or greater, not -3

        sage: E = elementary_matrix(ZZ, 5, row2=1, scale=12)
        Traceback (most recent call last):
        ...
        ValueError: row1 or col1 must be specified

        sage: E = elementary_matrix(ZZ, 5, row1=3, col1=3, scale=12)
        Traceback (most recent call last):
        ...
        ValueError: cannot specify both row1 and col1

        sage: E = elementary_matrix(ZZ, 5, row1=4/3, row2=1, scale=12)
        Traceback (most recent call last):
        ...
        TypeError: row of elementary matrix must be an integer, not 4/3

        sage: E = elementary_matrix(ZZ, 5, col1=5, col2=1, scale=12)
        Traceback (most recent call last):
        ...
        ValueError: column of elementary matrix must be positive and smaller than 5, not 5

        sage: E = elementary_matrix(ZZ, 5, col1=3, col2=4/3, scale=12)
        Traceback (most recent call last):
        ...
        TypeError: column of elementary matrix must be an integer, not 4/3

        sage: E = elementary_matrix(ZZ, 5, row1=3, row2=-1, scale=12)
        Traceback (most recent call last):
        ...
        ValueError: row of elementary matrix must be positive and smaller than 5, not -1

        sage: E = elementary_matrix(ZZ, 5, row1=3, row2=1, scale=4/3)
        Traceback (most recent call last):
        ...
        TypeError: scale parameter of elementary matrix must an element of Integer Ring, not 4/3

        sage: E = elementary_matrix(ZZ, 5, row1=3)
        Traceback (most recent call last):
        ...
        ValueError: insufficient parameters provided to construct elementary matrix

        sage: E = elementary_matrix(ZZ, 5, row1=3, row2=3, scale=12)
        Traceback (most recent call last):
        ...
        ValueError: cannot add a multiple of a row to itself

        sage: E = elementary_matrix(ZZ, 5, col1=3, scale=0)
        Traceback (most recent call last):
        ...
        ValueError: scale parameter of column of elementary matrix must be nonzero

    AUTHOR:

    - Rob Beezer (2011-03-04)
    """
    # determine ring and matrix size
    if arg1 is not None and arg0 not in Rings():
        raise TypeError('optional first parameter must be a ring, not {0}'.format(arg0))
    scale = kwds.pop('scale', None)
    if arg0 in Rings():
        R = arg0
        arg0 = arg1
    elif scale is not None:
        if not isinstance(scale, RingElement):
            raise TypeError('scale must be an element of some ring, not {0}'.format(scale))
        R = scale.parent()
    else:
        R = ZZ
    if arg0 is None:
        raise ValueError('size of elementary matrix must be given')
    try:
        n = Integer(arg0)
    except TypeError:
        raise TypeError('size of elementary matrix must be an integer, not {0}'.format(arg0))
    if n <= 0:
        raise ValueError('size of elementary matrix must be 1 or greater, not {0}'.format(n))
    # row operations or column operations?
    # column operation matrix will be transpose of a row operation matrix
    row1 = kwds.pop('row1', None)
    col1 = kwds.pop('col1', None)
    if row1 is None and col1 is None:
        raise ValueError('row1 or col1 must be specified')
    if row1 is not None and col1 is not None:
        raise ValueError('cannot specify both row1 and col1')
    rowop = row1 is not None
    if rowop:
        opstring = "row"
        row2 = kwds.pop('row2', None)
    else:
        opstring = "column"
        row1 = col1
        row2 = kwds.pop('col2', None)
    sparse = kwds.pop('sparse', False)
    if kwds:
        extra = kwds.popitem()
        raise ValueError('received an unexpected keyword: {0}={1}'.format(extra[0], extra[1]))

    # analyze parameters to determine matrix type
    try:
        row1 = Integer(row1)
    except TypeError:
        raise TypeError('{0} of elementary matrix must be an integer, not {1}'.format(opstring, row1))
    if row1 < 0 or row1 >= n:
        raise ValueError('{0} of elementary matrix must be positive and smaller than {1}, not {2}'.format(opstring, n, row1))
    if row2 is not None:
        try:
            row2 = Integer(row2)
        except TypeError:
            raise TypeError('{0} of elementary matrix must be an integer, not {1}'.format(opstring, row2))
        if row2 < 0 or row2 >= n:
            raise ValueError('{0} of elementary matrix must be positive and smaller than {1}, not {2}'.format(opstring, n, row2))
    if scale is not None:
        try:
            scale = R(scale)
        except Exception:
            raise TypeError('scale parameter of elementary matrix must an element of {0}, not {1}'.format(R, scale))

    # determine type of matrix and adjust an identity matrix
    # return row operation matrix or the transpose as a column operation matrix
    elem = identity_matrix(R, n, sparse=sparse)
    if row2 is None and scale is None:
        raise ValueError('insufficient parameters provided to construct elementary matrix')
    elif row2 is not None and scale is not None:
        if row1 == row2:
            raise ValueError('cannot add a multiple of a {0} to itself'.format(opstring))
        elem[row1, row2] = scale
    elif row2 is not None and scale is None:
        elem[row1, row1] = 0
        elem[row2, row2] = 0
        elem[row1, row2] = 1
        elem[row2, row1] = 1
    elif row2 is None and scale is not None:
        if scale == 0:
            raise ValueError('scale parameter of {0} of elementary matrix must be nonzero'.format(opstring))
        elem[row1, row1] = scale
    if rowop:
        return elem
    else:
        return elem.transpose()


@matrix_method
def circulant(v, sparse=None):
    r"""
    Return the circulant matrix specified by its 1st row `v`.

    A circulant `n \times n` matrix specified by the 1st row `v=(v_0...v_{n-1})` is
    the matrix `(c_{ij})_{0 \leq i,j\leq n-1}`, where `c_{ij}=v_{j-i \mod b}`.

    INPUT:

    - ``v`` -- list or a vector of values

    - ``sparse`` -- ``None`` by default; if ``sparse`` is set to ``True``, the output
      will be sparse.  Respectively, setting it to ``False`` produces dense output.
      If ``sparse`` is not set, and if ``v`` is a vector, the output sparsity is determined
      by the sparsity of ``v``; else, the output will be dense.

    EXAMPLES::

        sage: v = [1,2,3,4,8]
        sage: matrix.circulant(v)
        [1 2 3 4 8]
        [8 1 2 3 4]
        [4 8 1 2 3]
        [3 4 8 1 2]
        [2 3 4 8 1]
        sage: m = matrix.circulant(vector(GF(3),[0,1,-1],sparse=True)); m
        [0 1 2]
        [2 0 1]
        [1 2 0]
        sage: m.is_sparse()
        True

    TESTS::

        sage: m = matrix.circulant(vector(GF(3),[0,1,-1],sparse=False))
        sage: m.is_sparse()
        False
        sage: matrix.circulant([0,1,-1]).is_sparse()
        False
        sage: matrix.circulant([0,1,-1], sparse=True).is_sparse()
        True
    """
    if sparse is None:
        try:
            sparse = v.is_sparse()
        except AttributeError:
            sparse = False
    n = len(v)
    return matrix(n, n, lambda i, j: v[(j - i) % n], sparse=sparse)


def _determine_block_matrix_grid(sub_matrices):
    r"""
    For internal use. This tries to determine the dimensions
    of rows/columns when assembling the matrices in sub_matrices in a
    rectangular grid. It returns a pair of lists containing
    respectively the sizes of rows and columns.

    sub_matrices must be a list of lists of matrices. All sublists
    are expected to be the same size.

    Non-zero scalars are considered to be square matrices of any size,
    and zeroes are considered to be zero matrices of any size.

    A :exc:`ValueError` is raised if there is insufficient or
    conflicting information.

    TESTS::

        sage: from sage.matrix.special import _determine_block_matrix_grid
        sage: A = matrix(QQ, 2, 2, [3,9,6,10])
        sage: _determine_block_matrix_grid([[A, A], [A, A]])
        ([2, 2], [2, 2])

        sage: B = matrix(QQ, 1, 1, [ 1 ] )
        sage: C = matrix(QQ, 2, 2, [ 2, 3, 4, 5 ] )
        sage: _determine_block_matrix_grid([[B, 0], [0, C]])
        ([1, 2], [1, 2])
    """

    nrows = len(sub_matrices)
    if nrows == 0:
        return ([], [])
    ncols = len(sub_matrices[0])

    if ncols == 0:
        return ([0] * nrows, [])

    row_heights = [None] * nrows
    col_widths = [None] * ncols

    changing = True
    while changing:
        changing = False
        for i in range(nrows):
            for j in range(ncols):
                M = sub_matrices[i][j]
                sub_width = None
                sub_height = None
                if isinstance(M, Matrix):
                    sub_width = M.ncols()
                    sub_height = M.nrows()
                elif M:  # nonzero scalar is interpreted as a square matrix
                    if row_heights[i] is None:
                        sub_width = col_widths[j]
                    else:
                        sub_width = row_heights[i]
                    sub_height = sub_width
                if sub_width is not None:
                    if col_widths[j] is None:
                        changing = True
                        col_widths[j] = sub_width
                    elif col_widths[j] != sub_width:
                        raise ValueError("incompatible submatrix widths")
                if sub_height is not None:
                    if row_heights[i] is None:
                        changing = True
                        row_heights[i] = sub_height
                    elif row_heights[i] != sub_height:
                        raise ValueError("incompatible submatrix heights")

    if None in row_heights or None in col_widths:
        if None in row_heights or None in col_widths:
            raise ValueError("insufficient information to determine dimensions.")

    return (row_heights, col_widths)


def _determine_block_matrix_rows(sub_matrices):
    """
    For internal use. This tests if the matrices in sub_matrices
    fit in a rectangular matrix when assembled a row at a time.

    sub_matrices must be a list of lists of matrices.

    It returns a pair (row_heights, zero_widths, width) where
    row_heights is the list of row heights, zero_widths is the
    total width filled up by zero matrices per row, and width
    is the total width of the resulting matrix.

    Non-zero scalars are considered to be square matrices of any size,
    and zeroes are considered to be zero matrices of any size.

    A :exc:`ValueError` is raised if there is insufficient or
    conflicting information.

    TESTS::

        sage: from sage.matrix.special import _determine_block_matrix_rows
        sage: A = Matrix(ZZ, 1, 4, [1, 2, 3, 4])
        sage: _determine_block_matrix_rows([ [1, 1], [ A ] ])
        ([2, 1], [0, 0], 4)

        sage: B = Matrix(ZZ, 2, 2, [1, 2, 3, 4])
        sage: _determine_block_matrix_rows([ [B, B], [B, 1] ])
        ([2, 2], [0, 0], 4)
    """
    total_width = None
    row_heights = [None] * len(sub_matrices)
    zero_widths = [0] * len(sub_matrices)

    # We first do a pass to see if we can determine the width
    unknowns = False
    for i in range(len(sub_matrices)):
        R = sub_matrices[i]
        height = None
        # We first do a pass to see if we can determine the height
        # of this row
        found_zeroes = False
        for M in R:
            if isinstance(M, Matrix):
                if height is None:
                    height = M.nrows()
                elif height != M.nrows():
                    raise ValueError("incompatible submatrix heights")
            elif not M:
                found_zeroes = True
        if not R:
            height = 0

        # If we have a height, then we know the dimensions of any
        # nonzero scalars, and can maybe compute the width
        if height is not None and not found_zeroes:
            width = 0
            for M in R:
                if isinstance(M, Matrix):
                    width += M.ncols()
                else:
                    # nonzero scalar
                    width += height
            if total_width is None:
                total_width = width
            elif total_width != width:
                raise ValueError("incompatible submatrix widths")
            row_heights[i] = height
        else:
            # We don't set height here even if we know it,
            # to signal this row hasn't been fit yet.
            unknowns = True

    if total_width is None:
        raise ValueError("insufficient information to determine submatrix widths")

    if unknowns:
        # Do a second pass and see if the remaining rows can be
        # determined now that we know the width of the matrix.
        for i in range(len(sub_matrices)):
            if row_heights[i] is not None:
                continue
            R = sub_matrices[i]
            zero_state = 0
            # 0: no zeroes found
            # 1: consecutive zeroes found
            # 2: consecutive zeroes followed by nonzero found
            # 3: non-consecutive zeroes found
            scalars = 0
            width = 0
            height = None
            for j in range(len(R)):
                M = R[j]
                if isinstance(M, Matrix):
                    height = M.nrows()
                    width += M.ncols()
                    if zero_state == 1:
                        zero_state = 2
                elif not M:
                    if zero_state == 0:
                        zero_state = 1
                    elif zero_state == 2:
                        zero_state = 3
                else:
                    scalars += 1

            remaining_width = total_width - width
            # This remaining width has to be split over the
            # zeroes and (nonzero) scalars

            if height is not None:
                remaining_width -= scalars * height
                if remaining_width < 0:
                    raise ValueError("incompatible submatrix widths")
                if remaining_width > 0 and zero_state == 3:
                    raise ValueError("insufficient information to determine submatrix widths")
                if remaining_width > 0 and zero_state == 0:
                    raise ValueError("incompatible submatrix widths")
                # otherwise, things fit
                row_heights[i] = height
                zero_widths[i] = remaining_width
            elif zero_state != 0:
                # if we don't know the height, and there are zeroes,
                # we can't determine the height
                raise ValueError("insufficient information to determine submatrix heights")
            elif total_width % len(R):
                raise ValueError("incompatible submatrix widths")
            else:
                height = int(total_width / len(R))
                row_heights[i] = height

    # If we got this far, then everything fits
    return (row_heights, zero_widths, total_width)


@matrix_method
def block_matrix(*args, **kwds):
    r"""
    Return a larger matrix made by concatenating submatrices
    (rows first, then columns). For example, the matrix

    ::

        [ A B ]
        [ C D ]

    is made up of submatrices A, B, C, and D.

    INPUT:

    The block_matrix command takes a list of submatrices to add
    as blocks, optionally preceded by a ring and the number of block rows
    and block columns, and returns a matrix.

    The submatrices can be specified as a list of matrices (using
    ``nrows`` and ``ncols`` to determine their layout), or a list
    of lists of matrices, where each list forms a row.

    - ``ring`` -- the base ring

    - ``nrows`` -- the number of block rows

    - ``ncols`` -- the number of block cols

    - ``sub_matrices`` -- matrices (see below for syntax)

    - ``subdivide`` -- boolean, whether or not to add
      subdivision information to the matrix

    - ``sparse`` -- boolean, whether to make the resulting matrix sparse

    EXAMPLES::

        sage: A = matrix(QQ, 2, 2, [3,9,6,10])
        sage: block_matrix([ [A, -A], [~A, 100*A] ])
        [    3     9|   -3    -9]
        [    6    10|   -6   -10]
        [-----------+-----------]
        [-5/12   3/8|  300   900]
        [  1/4  -1/8|  600  1000]

    If the number of submatrices in each row is the same,
    you can specify the submatrices as a single list too::

        sage: block_matrix(2, 2, [ A, A, A, A ])
        [ 3  9| 3  9]
        [ 6 10| 6 10]
        [-----+-----]
        [ 3  9| 3  9]
        [ 6 10| 6 10]

    One can use constant entries::

        sage: block_matrix([ [1, A], [0, 1] ])
        [ 1  0| 3  9]
        [ 0  1| 6 10]
        [-----+-----]
        [ 0  0| 1  0]
        [ 0  0| 0  1]

    A zero entry may represent any square or non-square zero matrix::

        sage: B = matrix(QQ, 1, 1, [ 1 ] )
        sage: C = matrix(QQ, 2, 2, [ 2, 3, 4, 5 ] )
        sage: block_matrix([ [B, 0], [0, C] ])
        [1|0 0]
        [-+---]
        [0|2 3]
        [0|4 5]

    One can specify the number of rows or columns as keywords too::

        sage: block_matrix([A, -A, ~A, 100*A], ncols=4)
        [    3     9|   -3    -9|-5/12   3/8|  300   900]
        [    6    10|   -6   -10|  1/4  -1/8|  600  1000]

        sage: block_matrix([A, -A, ~A, 100*A], nrows=1)
        [    3     9|   -3    -9|-5/12   3/8|  300   900]
        [    6    10|   -6   -10|  1/4  -1/8|  600  1000]

    It handles base rings nicely too::

        sage: R.<x> = ZZ['x']
        sage: block_matrix(2, 2, [1/2, A, 0, x-1])
        [  1/2     0|    3     9]
        [    0   1/2|    6    10]
        [-----------+-----------]
        [    0     0|x - 1     0]
        [    0     0|    0 x - 1]

        sage: block_matrix(2, 2, [1/2, A, 0, x-1]).parent()
        Full MatrixSpace of 4 by 4 dense matrices over Univariate Polynomial Ring in x over Rational Field

    Subdivisions are optional. If they are disabled, the columns need not line up::

        sage: B = matrix(QQ, 2, 3, range(6))
        sage: block_matrix([ [~A, B], [B, ~A] ], subdivide=False)
        [-5/12   3/8     0     1     2]
        [  1/4  -1/8     3     4     5]
        [    0     1     2 -5/12   3/8]
        [    3     4     5   1/4  -1/8]

    Without subdivisions it also deduces dimensions for scalars if possible::

        sage: C = matrix(ZZ, 1, 2, range(2))
        sage: block_matrix([ [ C, 0 ], [ 3, 4 ], [ 5, 6, C ] ], subdivide=False )
        [0 1 0 0]
        [3 0 4 0]
        [0 3 0 4]
        [5 6 0 1]

    If all submatrices are sparse (unless there are none at all), the result
    will be a sparse matrix. Otherwise it will be dense by default. The
    ``sparse`` keyword can be used to override this::

        sage: A = Matrix(ZZ, 2, 2, [0, 1, 0, 0], sparse=True)
        sage: block_matrix([ [ A ], [ A ] ]).parent()
        Full MatrixSpace of 4 by 2 sparse matrices over Integer Ring
        sage: block_matrix([ [ A ], [ A ] ], sparse=False).parent()
        Full MatrixSpace of 4 by 2 dense matrices over Integer Ring

    Consecutive zero submatrices are consolidated.  ::

        sage: B = matrix(2, range(4))
        sage: C = matrix(2, 8, range(16))
        sage: block_matrix(2, [[B,0,0,B],[C]], subdivide=False)
        [ 0  1  0  0  0  0  0  1]
        [ 2  3  0  0  0  0  2  3]
        [ 0  1  2  3  4  5  6  7]
        [ 8  9 10 11 12 13 14 15]

    Ambiguity is not tolerated.  ::

        sage: B = matrix(2, range(4))
        sage: C = matrix(2, 8, range(16))
        sage: block_matrix(2, [[B,0,B,0],[C]], subdivide=False)
        Traceback (most recent call last):
        ...
        ValueError: insufficient information to determine submatrix widths

    Giving only a flat list of submatrices does not work::

        sage: A = matrix(2, 3, range(6))
        sage: B = matrix(3, 3, range(9))
        sage: block_matrix([A, A, B, B])
        Traceback (most recent call last):
        ...
        ValueError: must specify either nrows or ncols

    TESTS::

        sage: A = matrix(ZZ, 2, 2, [3,5,8,13])
        sage: block_matrix(A)
        [ 3  5]
        [ 8 13]
        sage: block_matrix([[A, 0r], [1r, A]])
        [ 3  5| 0  0]
        [ 8 13| 0  0]
        [-----+-----]
        [ 1  0| 3  5]
        [ 0  1| 8 13]
    """
    args = list(args)
    sparse = kwds.get('sparse', None)

    if not args:
        if sparse is not None:
            return matrix_space.MatrixSpace(ZZ, 0, 0, sparse=sparse)([])
        else:
            return matrix_space.MatrixSpace(ZZ, 0, 0)([])

    if len(args) >= 1 and args[0] in Rings():
        # A ring is specified
        if kwds.get('ring', args[0]) != args[0]:
            raise ValueError("base ring specified twice and they are different")
        ring = args[0]
        args.pop(0)
    else:
        ring = kwds.get('ring', None)

    if len(args) >= 1:
        try:
            nrows = int(args[0])
            args.pop(0)
            if kwds.get('nrows', nrows) != nrows:
                raise ValueError("number of rows specified twice and they are different")
        except TypeError:
            nrows = kwds.get('nrows', None)
    else:
        nrows = kwds.get('nrows', None)

    if len(args) >= 1:
        # check to see if additionally, the number of columns is specified
        try:
            ncols = int(args[0])
            args.pop(0)
            if kwds.get('ncols', ncols) != ncols:
                raise ValueError("number of columns specified twice and they are different")
        except TypeError:
            ncols = kwds.get('ncols', None)
    else:
        ncols = kwds.get('ncols', None)

    # Now we've taken care of initial ring, nrows, and ncols arguments.

    # Now the rest of the arguments are a list of rows, a flat list of
    # matrices, or a single value.

    if not args:
        args = [[]]
    if len(args) > 1:
        print(args)
        raise TypeError("invalid block_matrix invocation")

    sub_matrices = args[0]

    if isinstance(sub_matrices, Matrix):
        M = sub_matrices
        # a single matrix (check nrows/ncols/ring)
        if (nrows is not None and nrows != 1) or \
           (ncols is not None and ncols != 1):
            raise ValueError("invalid nrows/ncols passed to block_matrix")
        if ring is not None:
            M = M.change_ring(ring)
        if sparse is not None and M.is_sparse() != sparse:
            M = M.sparse_matrix() if sparse else M.dense_matrix()
        return M

    if not isinstance(sub_matrices, (list, tuple)):
        raise TypeError("invalid block_matrix invocation")

    subdivide = kwds.get('subdivide', True)

    # Will we try to place the matrices in a rectangular grid?
    try_grid = True

    if not sub_matrices:
        if (nrows is not None and nrows != 0) or \
           (ncols is not None and ncols != 0):
            raise ValueError("invalid nrows/ncols passed to block_matrix")
    elif isinstance(sub_matrices[0], (list, tuple)):
        # A list of lists: verify all elements are lists, and if
        # ncols is set, the lengths match.
        if nrows is not None and len(sub_matrices) != nrows:
            raise ValueError("invalid nrows passed to block_matrix")
        first_len = len(sub_matrices[0])
        if ncols is not None and first_len != ncols:
            raise ValueError("invalid ncols passed to block_matrix")
        same_length = all(isinstance(v, (list, tuple)) and len(v) == first_len for v in sub_matrices)
        if subdivide and not same_length:
            raise ValueError("list of rows is not valid (rows are wrong types or lengths)")
        try_grid = same_length
    else:
        # A flat list
        # determine the block dimensions
        n = len(sub_matrices)
        if nrows is None:
            if ncols is None:
                raise ValueError("must specify either nrows or ncols")
            else:
                nrows = n // ncols
        elif ncols is None:
            ncols = n // nrows
        if nrows * ncols != n:
            raise ValueError("given number of rows (%s), columns (%s) incompatible with number of submatrices (%s)" % (nrows, ncols, n))
        # Now create a list of lists from this
        sub_matrices = [sub_matrices[i * ncols: (i + 1) * ncols]
                        for i in range(nrows)]

    # At this point sub_matrices is a list of lists

    # determine the base ring and sparsity
    if ring is None:
        ring = ZZ
        for row in sub_matrices:
            for M in row:
                R = M.base_ring() if isinstance(M, Matrix) else parent(M)
                if R is not ZZ:
                    ring = sage.categories.pushout.pushout(ring, R)

    if sparse is None:
        sparse = True
        for row in sub_matrices:
            for M in row:
                if sparse and isinstance(M, Matrix) and not M.is_sparse():
                    sparse = False

    row_heights = None
    col_widths = None
    zero_widths = None
    total_width = None

    # We first try to place the matrices in a rectangular grid
    if try_grid:
        try:
            (row_heights, col_widths) = _determine_block_matrix_grid(sub_matrices)
        except ValueError as e:
            if subdivide:
                raise ValueError(e)

    if col_widths is None:
        # Try placing the matrices in rows instead
        # (Only if subdivide is False)
        (row_heights, zero_widths, total_width) = _determine_block_matrix_rows(sub_matrices)

    # Success, so assemble the final matrix

    big = None
    for i in range(len(sub_matrices)):
        R = sub_matrices[i]
        row = None
        for j in range(len(R)):
            M = R[j]

            if isinstance(M, Matrix):
                if M.base_ring() is not ring:
                    M = M.change_ring(ring)
                if M.is_sparse() != sparse:
                    M = M.sparse_matrix() if sparse else M.dense_matrix()
            elif not M and zero_widths is not None:
                if zero_widths[i] > 0:
                    M = matrix(ring, row_heights[i], zero_widths[i], 0, sparse=sparse)
                    zero_widths[i] = 0
                else:
                    continue  # zero-width matrix
            else:
                if zero_widths is not None:
                    M = matrix(ring, row_heights[i], row_heights[i], M, sparse=sparse)
                else:
                    M = matrix(ring, row_heights[i], col_widths[j], M, sparse=sparse)

            # append M to this row
            if row is None:
                row = M
            else:
                row = row.augment(M)

        # append row to final matrix
        if big is None:
            big = row
        else:
            big = big.stack(row)

    if big is None:
        if ring is None:
            ring = ZZ
        big = matrix(ring, 0, 0)

    if subdivide:
        big.subdivide(running_total(row_heights[:-1]),
                      running_total(col_widths[:-1]))

    return big


@matrix_method
def block_diagonal_matrix(*sub_matrices, **kwds):
    """
    Create a block matrix whose diagonal block entries are given by
    sub_matrices, with zero elsewhere.

    See also :meth:`block_matrix`.

    EXAMPLES::

        sage: A = matrix(ZZ, 2, [1,2,3,4])
        sage: block_diagonal_matrix(A, A)
        [1 2|0 0]
        [3 4|0 0]
        [---+---]
        [0 0|1 2]
        [0 0|3 4]

    The sub-matrices need not be square::

        sage: B = matrix(QQ, 2, 3, range(6))
        sage: block_diagonal_matrix(~A, B)
        [  -2    1|   0    0    0]
        [ 3/2 -1/2|   0    0    0]
        [---------+--------------]
        [   0    0|   0    1    2]
        [   0    0|   3    4    5]
    """
    if isinstance(sub_matrices, (list, tuple)) and len(sub_matrices) == 1:
        sub_matrices = sub_matrices[0]
    n = len(sub_matrices)
    entries = [ZZ.zero()] * n**2
    for i in range(n):
        entries[n*i+i] = sub_matrices[i]
    return block_matrix(n, n, entries, **kwds)


@matrix_method
def jordan_block(eigenvalue, size, sparse=False):
    r"""
    Return the Jordan block for the given eigenvalue with given size.

    INPUT:

    - ``eigenvalue`` -- eigenvalue for the diagonal entries of the block
    - ``size`` -- size of the square matrix
    - ``sparse`` -- (default: ``False``) if ``True``, return a sparse matrix

    EXAMPLES::

        sage: jordan_block(5, 3)
        [5 1 0]
        [0 5 1]
        [0 0 5]

    TESTS::

        sage: jordan_block(6.2, 'junk')
        Traceback (most recent call last):
        ...
        TypeError: size of Jordan block needs to be an integer, not junk
        sage: jordan_block(6.2, -1)
        Traceback (most recent call last):
        ...
        ValueError: size of Jordan block must be nonnegative, not -1
    """
    try:
        size = ZZ(size)
    except TypeError:
        msg = "size of Jordan block needs to be an integer, not {0}"
        raise TypeError(msg.format(size))
    if size < 0:
        msg = "size of Jordan block must be nonnegative, not {0}"
        raise ValueError(msg.format(size))
    block = diagonal_matrix([eigenvalue] * size, sparse=sparse)
    for i in range(size - 1):
        block[i, i + 1] = 1
    return block


@matrix_method
def companion_matrix(poly, format='right'):
    r"""
    Create a companion matrix from a monic polynomial.

    INPUT:

    - ``poly`` -- a univariate polynomial, or an iterable containing
      the coefficients of a polynomial, with low-degree coefficients first.
      The polynomial (or the polynomial implied by the coefficients) must
      be monic.  In other words, the leading coefficient must be one.
      A symbolic expression that might also be a polynomial is not
      proper input, see examples below.

    - ``format`` -- (default: ``'right'``) specifies one of four
      variations of a companion matrix.  Allowable values are
      ``'right'``, ``'left'``, ``'top'`` and ``'bottom'``, which indicates
      which border of the matrix contains the negatives of the coefficients.

    OUTPUT:

    A square matrix with a size equal to the degree of the polynomial.
    The returned matrix has ones above, or below the diagonal, and the
    negatives of the coefficients along the indicated border of the
    matrix (excepting the leading one coefficient).
    See the first examples below for precise illustrations.

    EXAMPLES:

    Each of the four possibilities.  Notice that the coefficients are
    specified and their negatives become the entries of the matrix.  The
    leading one must be given, but is not used.  The permutation matrix
    ``P`` is the identity matrix, with the columns reversed.  The last three
    statements test the general relationships between the four variants.  ::

        sage: poly = [-2, -3, -4, -5, -6, 1]
        sage: R = companion_matrix(poly, format='right'); R
        [0 0 0 0 2]
        [1 0 0 0 3]
        [0 1 0 0 4]
        [0 0 1 0 5]
        [0 0 0 1 6]
        sage: L = companion_matrix(poly, format='left'); L
        [6 1 0 0 0]
        [5 0 1 0 0]
        [4 0 0 1 0]
        [3 0 0 0 1]
        [2 0 0 0 0]
        sage: B = companion_matrix(poly, format='bottom'); B
        [0 1 0 0 0]
        [0 0 1 0 0]
        [0 0 0 1 0]
        [0 0 0 0 1]
        [2 3 4 5 6]
        sage: T = companion_matrix(poly, format='top'); T
        [6 5 4 3 2]
        [1 0 0 0 0]
        [0 1 0 0 0]
        [0 0 1 0 0]
        [0 0 0 1 0]

        sage: perm = Permutation([5, 4, 3, 2, 1])
        sage: P = perm.to_matrix()
        sage: L == P*R*P
        True
        sage: B == R.transpose()
        True
        sage: T == P*R.transpose()*P
        True

    A polynomial may be used as input, however a symbolic expression,
    even if it looks like a polynomial, is not regarded as such when used
    as input to this routine.  Obtaining the list of coefficients from a
    symbolic polynomial is one route to the companion matrix. ::

        sage: x = polygen(QQ, 'x')
        sage: p = x^3 - 4*x^2 + 8*x - 12
        sage: companion_matrix(p)
        [ 0  0 12]
        [ 1  0 -8]
        [ 0  1  4]

        sage: # needs sage.symbolic
        sage: y = var('y')
        sage: q = y^3 - 2*y + 1
        sage: companion_matrix(q)
        Traceback (most recent call last):
        ...
        TypeError: input must be a polynomial (not a symbolic expression, see docstring),
        or other iterable, not y^3 - 2*y + 1
        sage: coeff_list = [q(y=0)] + [q.coefficient(y^k)
        ....:                          for k in range(1, q.degree(y) + 1)]
        sage: coeff_list
        [1, -2, 0, 1]
        sage: companion_matrix(coeff_list)
        [ 0  0 -1]
        [ 1  0  2]
        [ 0  1  0]

    The minimal polynomial of a companion matrix is equal to the
    polynomial used to create it.  Used in a block diagonal
    construction, they can be used to create matrices with
    any desired minimal polynomial, or characteristic polynomial.  ::

        sage: t = polygen(QQ, 't')
        sage: p = t^12 - 7*t^4 + 28*t^2 - 456
        sage: C = companion_matrix(p, format='top')
        sage: q = C.minpoly(var='t'); q                                                 # needs sage.libs.pari
        t^12 - 7*t^4 + 28*t^2 - 456
        sage: p == q                                                                    # needs sage.libs.pari
        True

        sage: p = t^3 + 3*t - 8
        sage: q = t^5 + t - 17
        sage: A = block_diagonal_matrix( companion_matrix(p),
        ....:                            companion_matrix(p^2),
        ....:                            companion_matrix(q),
        ....:                            companion_matrix(q) )
        sage: A.charpoly(var='t').factor()                                              # needs sage.libs.pari
        (t^3 + 3*t - 8)^3 * (t^5 + t - 17)^2
        sage: A.minpoly(var='t').factor()                                               # needs sage.libs.pari
        (t^3 + 3*t - 8)^2 * (t^5 + t - 17)

    TESTS::

        sage: companion_matrix([4, 5, 1], format='junk')
        Traceback (most recent call last):
        ...
        ValueError: format must be 'right', 'left', 'top' or 'bottom', not junk

        sage: companion_matrix(sin(x))                                                  # needs sage.symbolic
        Traceback (most recent call last):
        ...
        TypeError: input must be a polynomial (not a symbolic expression, see docstring),
        or other iterable, not sin(x)

        sage: companion_matrix([2, 3, 896])
        Traceback (most recent call last):
        ...
        ValueError: polynomial (or the polynomial implied by coefficients) must be monic,
        not a leading coefficient of 896

        sage: F.<a> = GF(2^2)                                                           # needs sage.rings.finite_rings
        sage: companion_matrix([4/3, a+1, 1])                                           # needs sage.rings.finite_rings
        Traceback (most recent call last):
        ...
        TypeError: unable to find common ring for coefficients from polynomial

        sage: A = companion_matrix([1])
        sage: A.nrows(); A.ncols()
        0
        0

        sage: A = companion_matrix([])
        Traceback (most recent call last):
        ...
        ValueError: polynomial cannot be specified by an empty list

        sage: companion_matrix([QQ.one()]).parent()
        Full MatrixSpace of 0 by 0 dense matrices over Rational Field

    AUTHOR:

    - Rob Beezer (2011-05-19)
    """
    import sage.matrix.constructor
    if format not in ['right', 'left', 'top', 'bottom']:
        raise ValueError("format must be 'right', 'left', 'top' or 'bottom', not {0}".format(format))
    try:
        poly = Sequence(poly)
    except TypeError:
        raise TypeError('input must be a polynomial (not a symbolic expression, see docstring), or other iterable, not {0}'.format(poly))
    n = len(poly) - 1
    if n == -1:
        raise ValueError('polynomial cannot be specified by an empty list')
    if not poly[n] == 1:
        raise ValueError('polynomial (or the polynomial implied by coefficients) must be monic, not a leading coefficient of {0}'.format(poly[n]))
    try:
        M = sage.matrix.constructor.matrix(poly.universe(), n, n)
    except TypeError:
        raise TypeError("unable to find common ring for coefficients from polynomial")
    # 1s below diagonal, or above diagonal
    if format in ['right', 'top']:
        for i in range(n - 1):
            M[i+1, i] = 1
    else:
        for i in range(n-1):
            M[i, i+1] = 1
    # right side, left side (reversed), bottom edge, top edge (reversed)
    if format == 'right':
        for i in range(n):
            M[i, n-1] = -poly[i]
    elif format == 'left':
        for i in range(n):
            M[n-1-i, 0] = -poly[i]
    elif format == 'bottom':
        for i in range(n):
            M[n-1, i] = -poly[i]
    elif format == 'top':
        for i in range(n):
            M[0, n-1-i] = -poly[i]
    return M


@matrix_method
def random_rref_matrix(parent, num_pivots):
    r"""
    Generate a matrix in reduced row-echelon form with a specified number of nonzero rows.

    INPUT:

    - ``parent`` -- a matrix space specifying the base ring, dimensions and
      representation (dense/sparse) for the result.  The base ring must be exact.

    - ``num_pivots`` -- the number of nonzero rows in the result, i.e. the rank

    OUTPUT:

    A matrix in reduced row echelon form with ``num_pivots`` nonzero rows. If the
    base ring is `ZZ` or `QQ` then the entries are all integers.

    .. NOTE::

        It is easiest to use this function via a call to the
        :func:`~sage.matrix.constructor.random_matrix`
        function with the ``algorithm='echelon_form'`` keyword.  We provide
        one example accessing this function directly, while the remainder will
        use this more general function.

    EXAMPLES:

    Matrices generated are in reduced row-echelon form with specified rank. If the
    base ring is `QQ` the result has only integer entries.  ::

        sage: from sage.matrix.constructor import random_rref_matrix
        sage: matrix_space = sage.matrix.matrix_space.MatrixSpace(QQ, 5, 6)
        sage: A = random_rref_matrix(matrix_space, num_pivots=4); A  # random
        [ 1  0  0 -6  0 -3]
        [ 0  1  0  2  0  3]
        [ 0  0  1 -4  0 -2]
        [ 0  0  0  0  1  3]
        [ 0  0  0  0  0  0]
        sage: A.base_ring()
        Rational Field
        sage: (A.nrows(), A.ncols())
        (5, 6)
        sage: A in sage.matrix.matrix_space.MatrixSpace(ZZ, 5, 6)
        True
        sage: A.rank()
        4
        sage: A == A.rref()
        True

    Matrices can be generated over other exact rings. ::

        sage: B = random_matrix(FiniteField(7), 4, 4,  # random
        ....:                   algorithm='echelon_form', num_pivots=3); B
        [1 0 0 0]
        [0 1 0 6]
        [0 0 1 1]
        [0 0 0 0]
        sage: B.rank() == 3
        True
        sage: B.base_ring()
        Finite Field of size 7
        sage: B == B.rref()
        True

    TESTS:

    Rank zero::

        sage: random_matrix(QQ, 1, 1, algorithm='echelon_form', num_pivots=0)
        [0]

    Rank of a matrix must be an integer. ::

        sage: random_matrix(QQ, 120, 56, algorithm='echelon_form', num_pivots=61/2)
        Traceback (most recent call last):
        ...
        TypeError: the number of pivots must be an integer

    Matrices must be generated over exact fields. ::

        sage: random_matrix(RR, 40, 88, algorithm='echelon_form', num_pivots=39)
        Traceback (most recent call last):
        ...
        TypeError: the base ring must be exact

    Matrices must have the number of pivot columns be less than or equal to the number of rows. ::

        sage: C = random_matrix(ZZ, 6,4, algorithm='echelon_form', num_pivots=7); C
        Traceback (most recent call last):
        ...
        ValueError: number of pivots cannot exceed the number of rows or columns

    Matrices must have the number of pivot columns be less than or equal to the number of columns. ::

        sage: D = random_matrix(QQ, 1,3, algorithm='echelon_form', num_pivots=5); D
        Traceback (most recent call last):
        ...
        ValueError: number of pivots cannot exceed the number of rows or columns

    Matrices must have the number of pivot columns be greater than zero. ::

        sage: random_matrix(QQ, 5, 4, algorithm='echelon_form', num_pivots=-1)
        Traceback (most recent call last):
        ...
        ValueError: the number of pivots must be zero or greater

    AUTHOR:

    Billy Wonderly (2010-07)
    """
    import sage.probability.probability_distribution as pd

    try:
        num_pivots = ZZ(num_pivots)
    except TypeError:
        raise TypeError("the number of pivots must be an integer")
    if num_pivots < 0:
        raise ValueError("the number of pivots must be zero or greater")
    ring = parent.base_ring()
    if not ring.is_exact():
        raise TypeError("the base ring must be exact")
    num_row = parent.nrows()
    num_col = parent.ncols()
    if num_pivots > num_row or num_pivots > num_col:
        raise ValueError("number of pivots cannot exceed the number of rows or columns")

    if num_pivots == 0:
        return parent.zero()

    one = ring.one()
    # Create a matrix of the desired size to be modified and then returned.
    return_matrix = copy(parent.zero_matrix())

    # No harm if no pivots at all.
    subset = list(range(1, num_col))
    shuffle(subset)
    pivots = [0] + sorted(subset[:num_pivots - 1])

    # Use the list of pivot columns to set the pivot entries of the return_matrix to leading ones.
    for pivot_row, pivot in enumerate(pivots):
        return_matrix[pivot_row, pivot] = one
    if ring is QQ or ring is ZZ:
        # Keep track of the non-pivot columns by using the pivot_index, start at the first column to
        # the right of the initial pivot column, go until the first column to the left of the next
        # pivot column.
        for pivot_index in range(num_pivots - 1):
            for non_pivot_column_index in range(pivots[pivot_index] + 1, pivots[pivot_index + 1]):
                entry_generator1 = pd.RealDistribution("beta", [6, 4])
                # Experimental distribution used to generate the values.
                for non_pivot_column_entry in range(pivot_index + 1):
                    sign1 = (2 * randint(0, 1) - 1)
                    return_matrix[non_pivot_column_entry, non_pivot_column_index] = sign1 * int(entry_generator1.get_random_element() * ((1 - non_pivot_column_entry / return_matrix.ncols()) * 7))
        # Use index to fill entries of the columns to the right of the last pivot column.
        for rest_non_pivot_column in range(pivots[num_pivots - 1] + 1, num_col):
            entry_generator2 = pd.RealDistribution("beta", [2.6, 4])
            # experimental distribution to generate small values.
            for rest_entries in range(num_pivots):
                sign2 = (2 * randint(0, 1) - 1)
                return_matrix[rest_entries, rest_non_pivot_column] = sign2 * int(entry_generator2.get_random_element() * 5)
    else:
        for pivot_index in range(num_pivots - 1):
            for non_pivot_column_index in range(pivots[pivot_index] + 1, pivots[pivot_index + 1]):
                for non_pivot_column_entry in range(pivot_index + 1):
                    return_matrix[non_pivot_column_entry, non_pivot_column_index] = ring.random_element()
        for rest_non_pivot_column in range(pivots[num_pivots - 1] + 1, num_col):
            for rest_entries in range(num_pivots):
                return_matrix[rest_entries, rest_non_pivot_column] = ring.random_element()
    return return_matrix


@matrix_method
def random_echelonizable_matrix(parent, rank, upper_bound=None, max_tries=100):
    r"""
    Generate a matrix of a desired size and rank, over a desired ring, whose reduced
    row-echelon form has only integral values.

    INPUT:

    - ``parent`` -- a matrix space specifying the base ring, dimensions and
      representation (dense/sparse) for the result.  The base ring must be exact.

    - ``rank`` -- rank of result, i.e the number of nonzero rows in the
      reduced row echelon form

    - ``upper_bound`` -- if designated, size control of the matrix entries is desired
      Set ``upper_bound`` to 1 more than the maximum value entries can achieve.
      If ``None``, no size control occurs. But see the warning below.  (default: ``None``)

    - ``max_tries`` -- if designated, number of tries used to generate each new random row;s
      only matters when upper_bound!=None. Used to prevent endless looping. (default: 100)

    OUTPUT: a matrix not in reduced row-echelon form with the desired
    dimensions and properties

    .. warning::

        When ``upper_bound`` is set, it is possible for this constructor to
        fail with a :exc:`ValueError`.  This may happen when the ``upper_bound``,
        ``rank`` and/or matrix dimensions are all so small that it becomes
        infeasible or unlikely to create the requested matrix.  If you *must*
        have this routine return successfully, do not set ``upper_bound``.

    .. NOTE::

        It is easiest to use this function via a call to the
        :func:`~sage.matrix.constructor.random_matrix`
        function with the ``algorithm='echelonizable'`` keyword.  We provide
        one example accessing this function directly, while the remainder will
        use this more general function.

    EXAMPLES:

    Generated matrices have the desired dimensions, rank and entry size. The
    matrix in reduced row-echelon form has only integer entries. ::

        sage: from sage.matrix.constructor import random_echelonizable_matrix
        sage: matrix_space = sage.matrix.matrix_space.MatrixSpace(QQ, 5, 6)
        sage: A = random_echelonizable_matrix(matrix_space, rank=4, upper_bound=40)
        sage: A.rank()
        4
        sage: max(map(abs,A.list())) < 40
        True
        sage: A.rref() == A.rref().change_ring(ZZ)
        True

    An example with default settings (i.e. no entry size control). ::

        sage: C = random_matrix(QQ, 6, 7, algorithm='echelonizable', rank=5)
        sage: C.rank()
        5
        sage: C.rref() == C.rref().change_ring(ZZ)
        True

    A matrix without size control may have very large entry sizes. ::

        sage: D = random_matrix(ZZ, 7, 8, algorithm='echelonizable', rank=6); D  # random
        [    1     2     8   -35  -178  -239  -284   778]
        [    4     9    37  -163  -827 -1111 -1324  3624]
        [    5     6    21   -88  -454  -607  -708  1951]
        [   -4    -5   -22    97   491   656   779 -2140]
        [    4     4    13   -55  -283  -377  -436  1206]
        [    4    11    43  -194  -982 -1319 -1576  4310]
        [   -1    -2   -13    59   294   394   481 -1312]

    Matrices can be generated over any exact ring. ::

        sage: # needs sage.rings.finite_rings
        sage: F.<a> = GF(2^3)
        sage: B = random_matrix(F, 4, 5, algorithm='echelonizable', rank=4,
        ....:                   upper_bound=None)
        sage: B.rank()
        4
        sage: B.base_ring() is F
        True

    Square matrices over ZZ or QQ with full rank are always unimodular. ::

        sage: E = random_matrix(QQ, 7, 7, algorithm='echelonizable', rank=7)
        sage: det(E)
        1
        sage: E = random_matrix(ZZ, 7, 7, algorithm='echelonizable', rank=7)
        sage: det(E)
        1

    TESTS:

    Matrices must have a rank zero or greater, and less than
    both the number of rows and the number of columns. ::

        sage: random_matrix(QQ, 3, 4, algorithm='echelonizable', rank=-1)
        Traceback (most recent call last):
        ...
        ValueError: matrices must have rank zero or greater.
        sage: random_matrix(QQ, 3, 8, algorithm='echelonizable', rank=4)
        Traceback (most recent call last):
        ...
        ValueError: matrices cannot have rank greater than min(ncols,nrows).
        sage: random_matrix(QQ, 8, 3, algorithm='echelonizable', rank=4)
        Traceback (most recent call last):
        ...
        ValueError: matrices cannot have rank greater than min(ncols,nrows).

    The base ring must be exact. ::

        sage: random_matrix(RR, 3, 3, algorithm='echelonizable', rank=2)
        Traceback (most recent call last):
        ...
        TypeError: the base ring must be exact

    Works for rank==1, too. ::

        sage: random_matrix(QQ, 3, 3, algorithm='echelonizable', rank=1).ncols()
        3


    AUTHOR:

    Billy Wonderly (2010-07)
    """
    from sage.misc.prandom import randint

    ring = parent.base_ring()
    rows = parent.nrows()
    if rank < 0:
        raise ValueError("matrices must have rank zero or greater.")
    if rank > min(rows,parent.ncols()):
        raise ValueError("matrices cannot have rank greater than min(ncols,nrows).")
    matrix = random_rref_matrix(parent, rank)

    # Entries of matrices over the ZZ or QQ can get large, entry size is regulated by finding the largest
    # entry of the resultant matrix after addition of scalar multiple of a row.
    if ring is QQ or ring is ZZ:
        # If upper_bound is not set, don't control entry size.
        if upper_bound is None:
            # If size control is not desired, the routine will run slightly faster, particularly with large matrices.
            for pivots in range(rank-1, -1, -1):
                row_index = 0
                while row_index < rows:
                    if pivots == row_index:
                        row_index += 1
                    if pivots != row_index and row_index != rows:
                        matrix.add_multiple_of_row(row_index,
                                                   matrix.pivot_rows()[pivots],
                                                   randint(-5, 5))
                        row_index += 1
            if rows > 1:
                matrix.add_multiple_of_row(0, randint(1,rows-1), randint(-3,3))
        else:
            if rank == 1:  # would be better just to have a special generator...
                tries = 0
                while max(abs(c) for c in matrix.list()) >= upper_bound:
                    matrix = random_rref_matrix(parent, rank)
                    tries += 1
                    if tries > max_tries:  # to prevent endless attempts
                        raise ValueError("tried "+str(max_tries)+" times to get a rank 1 random matrix. Try bigger upper_bound?")
                matrix_copy = matrix

            for pivots in range(len(matrix.pivots()) - 1, -1, -1):
                # keep track of the pivot column positions from the pivot column with the largest index to
                # the one with the smallest.
                row_index = 0
                tries = 0
                while row_index < rows:
                    # To each row in a pivot column add a scalar multiple of the pivot row.
                    # for full rank, square matrices, using only this row operation preserves the determinant of 1.
                    if pivots != row_index:
                        # To ensure a leading one is not removed by the addition of the pivot row by its
                        # additive inverse.
                        matrix_copy = matrix.with_added_multiple_of_row(row_index,matrix.pivot_rows()[pivots],randint(-5,5))
                        tries += 1
                        # Range for scalar multiples determined experimentally.
                    if max(map(abs, matrix_copy.list())) < upper_bound:
                        # Continue if the largest entry after a row operation is within the bound.
                        matrix = matrix_copy
                        row_index += 1
                        tries = 0
                    if tries > max_tries:  # to prevent endless unsuccessful row adding
                        raise ValueError("tried "+str(max_tries)+" times to get row number "+str(row_index)+". Try bigger upper_bound?")
            # The leading one in row one has not been altered, so add a scalar multiple of a random row
            # to row one.
            row1 = 0
            if rows > 1:
                while row1 < 1:
                    matrix_copy = matrix.with_added_multiple_of_row(0,randint(1,rows-1),randint(-3,3))
                    if max(map(abs,matrix_copy.list())) < upper_bound:
                        matrix = matrix_copy
                        row1 += 1
    # If the matrix generated over a different ring, random elements from the designated ring are used as and
    # the routine is run similarly to the size unchecked version for rationals and integers.
    else:
        for pivots in range(rank-1,-1,-1):
            row_index = 0
            while row_index < rows:
                if pivots == row_index:
                    row_index += 1
                if pivots != row_index and row_index != rows:
                    matrix.add_multiple_of_row(row_index,matrix.pivot_rows()[pivots],ring.random_element())
                    row_index += 1
        if rows > 1:
            matrix.add_multiple_of_row(0,randint(1,rows-1),ring.random_element())
    return matrix


@matrix_method
def random_subspaces_matrix(parent, rank=None):
    r"""
    Create a matrix of the designated size and rank whose right and
    left null spaces, column space, and row space have desirable
    properties that simplify the subspaces.

    INPUT:

    - ``parent`` -- a matrix space specifying the base ring, dimensions, and
      representation (dense/sparse) for the result.  The base ring must be exact.

    - ``rank`` -- the desired rank of the return matrix (default: ``None``)

    OUTPUT:

    A matrix whose natural basis vectors for its four subspaces, when
    computed, have reasonably sized, integral valued, entries.

    .. NOTE::

        It is easiest to use this function via a call to the
        :func:`~sage.matrix.constructor.random_matrix`
        function with the ``algorithm='subspaces'`` keyword.  We provide
        one example accessing this function directly, while the remainder will
        use this more general function.

    EXAMPLES:

    A 6x8 matrix with designated rank of 3.  The four subspaces are
    determined using one simple routine in which we augment the
    original matrix with the equal row dimension identity matrix.  The
    resulting matrix is then put in reduced row-echelon form and the
    subspaces can then be determined by analyzing subdivisions of this
    matrix. See the four subspaces routine in [Bee]_ for more. ::

        sage: from sage.matrix.constructor import random_subspaces_matrix
        sage: matrix_space = sage.matrix.matrix_space.MatrixSpace(QQ, 6, 8)
        sage: B = random_subspaces_matrix(matrix_space, rank=3)
        sage: B.rank()
        3
        sage: B.nullity()
        3
        sage: (B.nrows(), B.ncols())
        (6, 8)
        sage: all(x in ZZ for x in B.list())
        True
        sage: B_expanded = B.augment(identity_matrix(6)).rref()
        sage: all(x in ZZ for x in B_expanded.list())
        True

    Check that we fixed :issue:`10543` (echelon forms should be immutable)::

        sage: B_expanded.is_immutable()
        True

    We want to modify B_expanded, so replace it with a copy::

        sage: B_expanded = copy(B_expanded)
        sage: B_expanded.subdivide(B.nrows()-B.nullity(), B.ncols())
        sage: C = B_expanded.subdivision(0, 0)
        sage: L = B_expanded.subdivision(1, 1)
        sage: B.right_kernel() == C.right_kernel()
        True
        sage: B.row_space() == C.row_space()
        True
        sage: B.column_space() == L.right_kernel()
        True
        sage: B.left_kernel() == L.row_space()
        True

    A matrix to show that the null space of the L matrix is the column space of the starting matrix. ::

        sage: A = random_matrix(QQ, 5, 7, algorithm='subspaces', rank=None)
        sage: (A.nrows(), A.ncols())
        (5, 7)
        sage: all(x in ZZ for x in A.list())
        True
        sage: A_expanded = A.augment(identity_matrix(5)).rref()
        sage: all(x in ZZ for x in A_expanded.list())
        True
        sage: C = A_expanded.submatrix(0, 0, A.nrows() - A.nullity(), A.ncols())
        sage: L = A_expanded.submatrix(A.nrows() - A.nullity(), A.ncols())
        sage: A.right_kernel() == C.right_kernel()
        True
        sage: A.row_space() == C.row_space()
        True
        sage: A.column_space() == L.right_kernel()
        True
        sage: A.left_kernel() == L.row_space()
        True

    TESTS:

    The designated rank of the L matrix cannot be greater than the
    number of desired rows, nor can the rank be negative. ::

        sage: random_matrix(QQ, 19, 20, algorithm='subspaces', rank=21)
        Traceback (most recent call last):
        ...
        ValueError: rank cannot exceed the number of rows or columns.
        sage: random_matrix(QQ, 19, 20, algorithm='subspaces', rank=-1)
        Traceback (most recent call last):
        ...
        ValueError: matrices must have rank zero or greater.

    AUTHOR:

    Billy Wonderly (2010-07)
    """

    import sage.probability.probability_distribution as pd

    ring = parent.base_ring()
    rows = parent.nrows()
    columns = parent.ncols()

    # If rank is not designated, generate using probability distribution
    # skewing to smaller numbers, always at least 1.
    if rank is None:
        left_nullity_generator = pd.RealDistribution("beta", [1.4, 5.5])
        nullity = int(left_nullity_generator.get_random_element()*(rows-1) + 1)
        rank = rows - nullity
    if rank < 0:
        raise ValueError("matrices must have rank zero or greater.")
    if rank > rows or rank > columns:
        raise ValueError("rank cannot exceed the number of rows or columns.")
    nullity = rows - rank
    B = random_matrix(ring, rows, columns, algorithm='echelon_form',
            num_pivots=rank)

    # Create a nonsingular matrix whose columns will be used to stack a matrix
    # over the L matrix, forming a nonsingular matrix.
    K_nonzero_columns = random_matrix(ring, rank, rank,
            algorithm='echelonizable', rank=rank)
    K = matrix(QQ, rank, rows)
    L = random_matrix(ring, nullity, rows, algorithm='echelon_form',
            num_pivots=nullity)
    for column in range(len(L.nonpivots())):
        for entry in range(rank):
            K[entry, L.nonpivots()[column]] = K_nonzero_columns[entry, column]
    J = K.stack(L)

    # By multiplying the B matrix by J.inverse() we hide the B matrix of the
    # solution using row operations required to change the solution K matrix to
    # the identity matrix.
    return J.inverse() * B


@matrix_method
def random_unimodular_matrix(parent, upper_bound=None, max_tries=100):
    r"""
    Generate a random unimodular (determinant 1) matrix of a desired size over a desired ring.

    INPUT:

    - ``parent`` -- a matrix space specifying the base ring, dimensions
      and representation (dense/sparse) for the result.  The base ring
      must be exact.

    - ``upper_bound`` -- for large matrices over QQ or ZZ,
      ``upper_bound`` is the largest value matrix entries can achieve.  But
      see the warning below.

    - ``max_tries`` -- if designated, number of tries used to generate each new random row;
      only matters when upper_bound!=None. Used to prevent endless looping. (default: 100)

    A matrix not in reduced row-echelon form with the desired dimensions and properties.

    OUTPUT: an invertible matrix with the desired properties and determinant 1

    .. warning::

        When ``upper_bound`` is set, it is possible for this constructor to
        fail with a :exc:`ValueError`.  This may happen when the ``upper_bound``,
        ``rank`` and/or matrix dimensions are all so small that it becomes
        infeasible or unlikely to create the requested matrix.  If you *must*
        have this routine return successfully, do not set ``upper_bound``.

    .. NOTE::

        It is easiest to use this function via a call to the
        :func:`~sage.matrix.constructor.random_matrix`
        function with the ``algorithm='unimodular'`` keyword.  We provide
        one example accessing this function directly, while the remainder will
        use this more general function.

    EXAMPLES:

    A matrix size 5 over QQ. ::

        sage: from sage.matrix.constructor import random_unimodular_matrix
        sage: matrix_space = sage.matrix.matrix_space.MatrixSpace(QQ, 5)
        sage: A = random_unimodular_matrix(matrix_space)
        sage: det(A)
        1

    A matrix size 6 with entries no larger than 50. ::

        sage: B = random_matrix(ZZ, 7, algorithm='unimodular', upper_bound=50)
        sage: det(B)
        1
        sage: all(abs(b) < 50 for b in B.list())
        True

    A matrix over the number Field in `y` with defining polynomial `y^2-2y-2`. ::

        sage: # needs sage.rings.number_field
        sage: y = polygen(ZZ, 'y')
        sage: K = NumberField(y^2 - 2*y - 2, 'y')
        sage: C = random_matrix(K, 3, algorithm='unimodular')
        sage: det(C)
        1
        sage: C.base_ring() is K
        True

    TESTS:

    Unimodular matrices are square. ::

        sage: random_matrix(QQ, 5, 6, algorithm='unimodular')
        Traceback (most recent call last):
        ...
        TypeError: a unimodular matrix must be square.

    Only matrices over ZZ and QQ can have size control. ::

        sage: F.<a> = GF(5^7)                                                           # needs sage.rings.finite_rings
        sage: random_matrix(F, 5, algorithm='unimodular', upper_bound=20)               # needs sage.rings.finite_rings
        Traceback (most recent call last):
        ...
        TypeError: only matrices over ZZ or QQ can have size control.

    AUTHOR:

    Billy Wonderly (2010-07)
    """

    ring = parent.base_ring()
    size = parent.nrows()
    if parent.nrows() != parent.ncols():
        raise TypeError("a unimodular matrix must be square.")
    if upper_bound is not None and (ring != ZZ and ring != QQ):
        raise TypeError("only matrices over ZZ or QQ can have size control.")
    if upper_bound is None:
        # random_echelonizable_matrix() always returns a determinant one matrix if given full rank.
        return random_matrix(ring, size, algorithm='echelonizable', rank=size)
    elif upper_bound is not None and (ring == ZZ or ring == QQ):
        return random_matrix(ring, size,algorithm='echelonizable',rank=size, upper_bound=upper_bound, max_tries=max_tries)


@matrix_method
def random_diagonalizable_matrix(parent,eigenvalues=None,dimensions=None):
    """
    Create a random matrix that diagonalizes nicely.

    To be used as a teaching tool.  Return matrices have only real
    eigenvalues.

    INPUT:

    If eigenvalues and dimensions are not specified in a list,
    they will be assigned randomly.

    - ``parent`` -- the desired size of the square matrix

    - ``eigenvalues`` -- the list of desired eigenvalues (default=None)

    - ``dimensions`` -- the list of dimensions corresponding to each
      eigenspace (default=None)

    OUTPUT:

    A square, diagonalizable, matrix with only integer entries. The
    eigenspaces of this matrix, if computed by hand, give basis
    vectors with only integer entries.

    .. NOTE::

        It is easiest to use this function via a call to the
        :func:`~sage.matrix.constructor.random_matrix`
        function with the ``algorithm='diagonalizable'`` keyword.  We provide
        one example accessing this function directly, while the remainder will
        use this more general function.

    EXAMPLES:

    A diagonalizable matrix, size 5. ::

        sage: from sage.matrix.constructor import random_diagonalizable_matrix
        sage: matrix_space = sage.matrix.matrix_space.MatrixSpace(QQ, 5)
        sage: A = random_diagonalizable_matrix(matrix_space)

        sage: # needs sage.rings.number_field
        sage: eigenvalues = A.eigenvalues()
        sage: S = A.right_eigenmatrix()[1]
        sage: eigenvalues2 = (S.inverse()*A*S).diagonal()
        sage: sorted(eigenvalues) == sorted(eigenvalues2)
        True

    A diagonalizable matrix with eigenvalues and dimensions designated,
    with a check that if eigenvectors were calculated by hand
    entries would all be integers. ::

        sage: eigenvalues = [ZZ.random_element() for _ in range(3)]
        sage: B = random_matrix(QQ, 6, algorithm='diagonalizable',
        ....:                   eigenvalues=eigenvalues, dimensions=[2,3,1])
        sage: all(x in ZZ for x in (B-(-12*identity_matrix(6))).rref().list())
        True
        sage: all(x in ZZ for x in (B-(4*identity_matrix(6))).rref().list())
        True
        sage: all(x in ZZ for x in (B-(6*identity_matrix(6))).rref().list())
        True

        sage: # needs sage.rings.number_field
        sage: S = B.right_eigenmatrix()[1]
        sage: eigenvalues2 = (S.inverse()*B*S).diagonal()
        sage: all(e in eigenvalues for e in eigenvalues2)
        True

    TESTS:

    Eigenvalues must all be integers. ::

        sage: random_matrix(QQ, 3, algorithm='diagonalizable',                          # needs sage.symbolic
        ....:               eigenvalues=[2+I, 2-I, 2], dimensions=[1,1,1])
        Traceback (most recent call last):
        ...
        TypeError: eigenvalues must be integers.

    Diagonal matrices must be square. ::

        sage: random_matrix(QQ, 5, 7, algorithm='diagonalizable', eigenvalues=[-5,2,-3], dimensions=[1,1,3])
        Traceback (most recent call last):
        ...
        TypeError: a diagonalizable matrix must be square.

    A list of eigenvalues must be accompanied with a list of dimensions. ::

        sage: random_matrix(QQ,10,algorithm='diagonalizable',eigenvalues=[4,8])
        Traceback (most recent call last):
        ...
        ValueError: the list of eigenvalues must have a list of dimensions corresponding to each eigenvalue.

    A list of dimensions must be accompanied with a list of eigenvalues. ::

        sage: random_matrix(QQ, 10,algorithm='diagonalizable',dimensions=[2,2,4,2])
        Traceback (most recent call last):
        ...
        ValueError: the list of dimensions must have a list of corresponding eigenvalues.

    The sum of the eigenvalue dimensions must equal the size of the matrix. ::

        sage: random_matrix(QQ,12,algorithm='diagonalizable',eigenvalues=[4,2,6,-1],dimensions=[2,3,5,1])
        Traceback (most recent call last):
        ...
        ValueError: the size of the matrix must equal the sum of the dimensions.

    Each eigenspace dimension must be at least 1. ::

        sage: random_matrix(QQ,9,algorithm='diagonalizable',eigenvalues=[-15,22,8,-4,90,12],dimensions=[4,2,2,4,-3,0])
        Traceback (most recent call last):
        ...
        ValueError: eigenspaces must have a dimension of at least 1.

    Each eigenvalue must have a corresponding eigenspace dimension. ::

        sage: random_matrix(QQ,12,algorithm='diagonalizable',eigenvalues=[4,2,6,-1],dimensions=[4,3,5])
        Traceback (most recent call last):
        ...
        ValueError: each eigenvalue must have a corresponding dimension and each dimension a corresponding eigenvalue.

    Each dimension must have an eigenvalue paired to it. ::

        sage: random_matrix(QQ,12,algorithm='diagonalizable',eigenvalues=[4,2,6],dimensions=[2,3,5,2])
        Traceback (most recent call last):
        ...
        ValueError: each eigenvalue must have a corresponding dimension and each dimension a corresponding eigenvalue.

    .. TODO::

        Modify the routine to allow for complex eigenvalues.

    AUTHOR:

    Billy Wonderly (2010-07)
    """

    from sage.misc.prandom import randint

    size = parent.nrows()
    if parent.nrows() != parent.ncols():
        raise TypeError("a diagonalizable matrix must be square.")
    if eigenvalues is not None and dimensions is None:
        raise ValueError("the list of eigenvalues must have a list of dimensions corresponding to each eigenvalue.")
    if eigenvalues is None and dimensions is not None:
        raise ValueError("the list of dimensions must have a list of corresponding eigenvalues.")
    if eigenvalues is None and dimensions is None:
        values = []
        # create a list with "size" number of entries
        for eigen_index in range(size):
            eigenvalue = randint(-10, 10)
            values.append(eigenvalue)
        values.sort()
        dimensions = []
        eigenvalues = []
        # create a list with no duplicate values to be the eigenvalues
        for eigenvalue in range(size):
            if values[eigenvalue] not in eigenvalues:
                eigenvalues.append(values[eigenvalue])
        for dimension in range(len(eigenvalues)):
            # dimension is equal to how many times an eigenvalue was generated in the 'values' list
            dimensions.append(values.count(eigenvalues[dimension]))
    size_check = 0
    for check in range(len(dimensions)):
        size_check = size_check + dimensions[check]
    if not all(x in ZZ for x in eigenvalues):
        raise TypeError("eigenvalues must be integers.")
    if size != size_check:
        raise ValueError("the size of the matrix must equal the sum of the dimensions.")
    if min(dimensions) < 1:
        raise ValueError("eigenspaces must have a dimension of at least 1.")
    if len(eigenvalues) != len(dimensions):
        raise ValueError("each eigenvalue must have a corresponding dimension and each dimension a corresponding eigenvalue.")
    # sort the dimensions in order of increasing size, and sort the eigenvalues list in an identical fashion, to maintain corresponding values.
    dimensions_sort = sorted(zip(dimensions, eigenvalues))
    dimensions = [x[0] for x in dimensions_sort]
    eigenvalues = [x[1] for x in dimensions_sort]
    # Create the matrix of eigenvalues on the diagonal.  Use a lower limit and upper limit determined by the eigenvalue dimensions.
    diagonal_matrix = matrix(QQ, size)
    up_bound = 0
    low_bound = 0
    for row_index in range(len(dimensions)):
        up_bound = up_bound + dimensions[row_index]
        for entry in range(low_bound,up_bound):
            diagonal_matrix[entry, entry] = eigenvalues[row_index]
        low_bound = low_bound+dimensions[row_index]
    # Create a matrix to hold each of the eigenvectors as its columns, begin with an identity matrix so that after row and column
    # operations the resulting matrix will be unimodular.
    eigenvector_matrix = matrix(QQ, size, size, 1)
    upper_limit = 0
    lower_limit = 0
    # run the routine over the necessary number of columns corresponding eigenvalue dimension.
    for dimension_index in range(len(dimensions)-1):
        upper_limit = upper_limit+dimensions[dimension_index]
        lowest_index_row_with_one = size-dimensions[dimension_index]
        # assign a one to the row that is the eigenvalue dimension rows up from the bottom row then assign ones diagonally down to the right.
        for eigen_ones in range(lower_limit,upper_limit):
            eigenvector_matrix[lowest_index_row_with_one,eigen_ones] = 1
            lowest_index_row_with_one += 1
        lower_limit = lower_limit+dimensions[dimension_index]
    # Create a list to give the eigenvalue dimension corresponding to each column.
    dimension_check = []
    for i in range(len(dimensions)):
        for k in range(dimensions[i]):
            dimension_check.append(dimensions[i])
    # run routine over the rows that are in the range of the protected ones.  Use addition of column multiples to fill entries.
    for dimension_multiplicity in range(max(dimensions),min(dimensions),-1):
        highest_one_row = size-dimension_multiplicity
        highest_one_column = 0
        # find the column with the protected one in the lowest indexed row.
        while eigenvector_matrix[highest_one_row,highest_one_column] == 0:
            highest_one_column += 1
        # dimension_check determines if column has a low enough eigenvalue dimension to take a column multiple.
        for bottom_entry_filler in range(len(dimension_check)):
            if dimension_check[bottom_entry_filler] < dimension_multiplicity and eigenvector_matrix[highest_one_row,bottom_entry_filler] == 0:
                # randint range determined experimentally to keep entries manageable.
                eigenvector_matrix.add_multiple_of_column(bottom_entry_filler,highest_one_column,randint(-4,4))
    # Fill remaining rows using scalar row addition.
    for row in range(size-max(dimensions),size):
        for upper_row in range(size-max(dimensions)):
            # range of multiplier determined experimentally so that entries stay manageable for small matrices
            eigenvector_matrix.add_multiple_of_row(upper_row,row,randint(-4,4))
    return eigenvector_matrix*diagonal_matrix*(eigenvector_matrix.inverse())


@matrix_method
def vector_on_axis_rotation_matrix(v, i, ring=None):
    r"""
    Return a rotation matrix `M` such that `det(M)=1` sending the vector
    `v` on the `i`-th axis so that all other coordinates of `Mv` are zero.

    .. NOTE::

        Such a matrix is not uniquely determined. This function returns one
        such matrix.

    INPUT:

    - ``v`` -- vector
    - ``i`` -- integer
    - ``ring`` -- ring (default: ``None``) of the resulting matrix

    OUTPUT: a matrix

    EXAMPLES::

        sage: from sage.matrix.constructor import vector_on_axis_rotation_matrix
        sage: v = vector((1,2,3))
        sage: vector_on_axis_rotation_matrix(v, 2) * v                                  # needs sage.symbolic
        (0, 0, sqrt(14))
        sage: vector_on_axis_rotation_matrix(v, 1) * v                                  # needs sage.symbolic
        (0, sqrt(14), 0)
        sage: vector_on_axis_rotation_matrix(v, 0) * v                                  # needs sage.symbolic
        (sqrt(14), 0, 0)

    ::

        sage: # needs sage.symbolic
        sage: x,y = var('x,y')
        sage: v = vector((x,y))
        sage: vector_on_axis_rotation_matrix(v, 1)
        [ y/sqrt(x^2 + y^2) -x/sqrt(x^2 + y^2)]
        [ x/sqrt(x^2 + y^2)  y/sqrt(x^2 + y^2)]
        sage: vector_on_axis_rotation_matrix(v, 0)
        [ x/sqrt(x^2 + y^2)  y/sqrt(x^2 + y^2)]
        [-y/sqrt(x^2 + y^2)  x/sqrt(x^2 + y^2)]
        sage: vector_on_axis_rotation_matrix(v, 0) * v
        (x^2/sqrt(x^2 + y^2) + y^2/sqrt(x^2 + y^2), 0)
        sage: vector_on_axis_rotation_matrix(v, 1) * v
        (0, x^2/sqrt(x^2 + y^2) + y^2/sqrt(x^2 + y^2))

    ::

        sage: v = vector((1,2,3,4))
        sage: vector_on_axis_rotation_matrix(v, 0) * v                                  # needs sage.symbolic
        (sqrt(30), 0, 0, 0)
        sage: vector_on_axis_rotation_matrix(v, 0, ring=RealField(10))
        [ 0.18  0.37  0.55  0.73]
        [-0.98 0.068  0.10  0.14]
        [ 0.00 -0.93  0.22  0.30]
        [ 0.00  0.00 -0.80  0.60]
        sage: vector_on_axis_rotation_matrix(v, 0, ring=RealField(10)) * v
        (5.5, 0.00..., 0.00..., 0.00...)

    AUTHORS:

    Sébastien Labbé (April 2010)
    """
    dim = len(v)
    v = vector(v)
    m = identity_matrix(dim, sparse=True)
    L = list(range(i - 1, -1, -1)) + list(range(dim - 1, i, -1))
    for i in L:
        rot = ith_to_zero_rotation_matrix(v, i, ring=ring)
        v = rot * v
        m = rot * m
    return m


@matrix_method
def ith_to_zero_rotation_matrix(v, i, ring=None):
    r"""
    Return a rotation matrix that sends the `i`-th coordinates of the
    vector v to zero by doing a rotation with the `(i-1)`-th coordinate.

    INPUT:

    - ``v`` -- vector
    - ``i`` -- integer
    - ``ring`` -- ring (default: ``None``) of the resulting matrix

    OUTPUT: a matrix

    EXAMPLES::

        sage: from sage.matrix.constructor import ith_to_zero_rotation_matrix
        sage: v = vector((1,2,3))
        sage: ith_to_zero_rotation_matrix(v, 2)                                         # needs sage.symbolic
        [             1              0              0]
        [             0  2/13*sqrt(13)  3/13*sqrt(13)]
        [             0 -3/13*sqrt(13)  2/13*sqrt(13)]
        sage: ith_to_zero_rotation_matrix(v, 2) * v                                     # needs sage.symbolic
        (1, sqrt(13), 0)

    ::

        sage: ith_to_zero_rotation_matrix(v, 0)                                         # needs sage.symbolic
        [ 3/10*sqrt(10)              0 -1/10*sqrt(10)]
        [             0              1              0]
        [ 1/10*sqrt(10)              0  3/10*sqrt(10)]
        sage: ith_to_zero_rotation_matrix(v, 1)                                         # needs sage.symbolic
        [ 1/5*sqrt(5)  2/5*sqrt(5)            0]
        [-2/5*sqrt(5)  1/5*sqrt(5)            0]
        [           0            0            1]
        sage: ith_to_zero_rotation_matrix(v, 2)                                         # needs sage.symbolic
        [             1              0              0]
        [             0  2/13*sqrt(13)  3/13*sqrt(13)]
        [             0 -3/13*sqrt(13)  2/13*sqrt(13)]

    ::

        sage: ith_to_zero_rotation_matrix(v, 0) * v                                     # needs sage.symbolic
        (0, 2, sqrt(10))
        sage: ith_to_zero_rotation_matrix(v, 1) * v                                     # needs sage.symbolic
        (sqrt(5), 0, 3)
        sage: ith_to_zero_rotation_matrix(v, 2) * v                                     # needs sage.symbolic
        (1, sqrt(13), 0)

    Other ring::

        sage: ith_to_zero_rotation_matrix(v, 2, ring=RR)
        [  1.00000000000000  0.000000000000000  0.000000000000000]
        [ 0.000000000000000  0.554700196225229  0.832050294337844]
        [ 0.000000000000000 -0.832050294337844  0.554700196225229]
        sage: ith_to_zero_rotation_matrix(v, 2, ring=RDF)
        [                1.0                 0.0                 0.0]
        [                0.0  0.5547001962252291  0.8320502943378437]
        [                0.0 -0.8320502943378437  0.5547001962252291]

    On the symbolic ring::

        sage: # needs sage.symbolic
        sage: x,y,z = var('x,y,z')
        sage: v = vector((x,y,z))
        sage: ith_to_zero_rotation_matrix(v, 2)
        [                 1                  0                  0]
        [                 0  y/sqrt(y^2 + z^2)  z/sqrt(y^2 + z^2)]
        [                 0 -z/sqrt(y^2 + z^2)  y/sqrt(y^2 + z^2)]
        sage: ith_to_zero_rotation_matrix(v, 2) * v
        (x, y^2/sqrt(y^2 + z^2) + z^2/sqrt(y^2 + z^2), 0)

    TESTS::

        sage: ith_to_zero_rotation_matrix((1,0,0), 0)
        [ 0  0 -1]
        [ 0  1  0]
        [ 1  0  0]
        sage: ith_to_zero_rotation_matrix((1,0,0), 1)
        [1 0 0]
        [0 1 0]
        [0 0 1]
        sage: ith_to_zero_rotation_matrix((1,0,0), 2)
        [1 0 0]
        [0 1 0]
        [0 0 1]

    AUTHORS:

    Sébastien Labbé (April 2010)
    """
    if ring is not None:
        # coerce the vector so that computations
        # are done in that ring
        v = vector(ring, v)
    dim = len(v)
    i = i % dim
    j = (i - 1) % dim
    a, b = v[j], v[i]
    if b == 0:
        return identity_matrix(dim, sparse=True)
    from sage.misc.functional import sqrt
    norm = sqrt(a * a + b * b)
    aa = a / norm
    bb = b / norm
    entries = {(k, k): 1 for k in range(dim)}
    entries.update({(j, j): aa, (j, i): bb, (i, j): -bb, (i, i): aa})
    return matrix(entries, nrows=dim, base_ring=ring)


@matrix_method
def hilbert(dim, ring=QQ):
    r"""
    Return a Hilbert matrix of the given dimension.

    The `n` dimensional Hilbert matrix is a square matrix with entries being
    unit fractions,

    .. MATH::

        H_{ij} = \frac{1}{i+j-1},\qquad i, j = 1,\ldots, n.

    For more information see the :wikipedia:`Hilbert_matrix`.

    INPUT:

    - ``dim`` -- integer; the dimension of the Hilbert matrix

    - ``ring`` -- base ring (default: `\QQ`) of the resulting matrix

    EXAMPLES::

        sage: matrix.hilbert(5)
        [  1 1/2 1/3 1/4 1/5]
        [1/2 1/3 1/4 1/5 1/6]
        [1/3 1/4 1/5 1/6 1/7]
        [1/4 1/5 1/6 1/7 1/8]
        [1/5 1/6 1/7 1/8 1/9]
    """
    def entries(i, j):
        return ZZ.one() / (i + j + 1)
    return matrix(entries, nrows=dim, ncols=dim, base_ring=ring)


@matrix_method
def vandermonde(v, ring=None):
    r"""
    Return a Vandermonde matrix of the given vector.

    The `n` dimensional Vandermonde matrix is a square matrix with columns
    being the powers of a given vector `v`,

    .. MATH::

        V_{ij} = v_i^{j-1},\qquad i, j = 1,\ldots, n.

    For more information see the :wikipedia:`Vandermonde_matrix`.

    INPUT:

    - ``v`` -- vector, the second column of the Vandermonde matrix

    - ``ring`` -- base ring (default: ``None``) of the resulting matrix

    EXAMPLES:

    A Vandermonde matrix of order three over the symbolic ring::

        sage: matrix.vandermonde(SR.var(['x0', 'x1', 'x2']))                            # needs sage.symbolic
        [   1   x0 x0^2]
        [   1   x1 x1^2]
        [   1   x2 x2^2]
    """
    def entries(i, j):
        return v[i]**j
    return matrix(entries, nrows=len(v), ncols=len(v), base_ring=ring)


@matrix_method
def toeplitz(c, r, ring=None):
    r"""
    Return a Toeplitz matrix of given first column and first row.

    In a Toeplitz matrix, each descending diagonal from left to right is
    constant, such that:

    .. MATH:: T_{i,j} = T_{i+1, j+1}.

    For more information see the :wikipedia:`Toeplitz_matrix`.

    INPUT:

    - ``c`` -- vector, first column of the Toeplitz matrix

    - ``r`` -- vector, first row of the Toeplitz matrix, counting from the
      second column

    - ``ring`` -- base ring (default: ``None``) of the resulting matrix

    EXAMPLES:

    A rectangular Toeplitz matrix::

        sage: matrix.toeplitz([1..4], [5..6])
        [1 5 6]
        [2 1 5]
        [3 2 1]
        [4 3 2]

    The following `N\times N` Toeplitz matrix arises in the discretization of
    boundary value problems::

        sage: N = 4
        sage: matrix.toeplitz([-2, 1] + [0]*(N-2), [1] + [0]*(N-2))
        [-2  1  0  0]
        [ 1 -2  1  0]
        [ 0  1 -2  1]
        [ 0  0  1 -2]
    """
    def entries(i, j):
        return c[i - j] if i >= j else r[j - i - 1]
    return matrix(entries, nrows=len(c), ncols=len(r)+1, base_ring=ring)


@matrix_method
def hankel(c, r=None, ring=None):
    r"""
    Return a Hankel matrix of given first column and whose elements are zero
    below the first anti-diagonal.

    The Hankel matrix is symmetric and constant across the anti-diagonals,
    with elements

    .. MATH::

        H_{ij} = v_{i+j-1},\qquad i = 1,\ldots, m,~j = 1,\ldots, n,

    where the vector `v_i = c_i` for `i = 1,\ldots, m` and `v_{m+i} = r_i` for
    `i = 1, \ldots, n-1` completely determines the Hankel matrix. If the last
    row, `r`, is not given, the Hankel matrix is square by default and `r = 0`.
    For more information see the :wikipedia:`Hankel_matrix`.

    INPUT:

    - ``c`` -- vector, first column of the Hankel matrix

    - ``r`` -- vector (default: ``None``); last row of the Hankel matrix, from
      the second to the last column

    - ``ring`` -- base ring (default: ``None``) of the resulting matrix

    EXAMPLES:

    A Hankel matrix with symbolic entries::

        sage: matrix.hankel(SR.var('a, b, c, d, e'))                                    # needs sage.symbolic
        [a b c d e]
        [b c d e 0]
        [c d e 0 0]
        [d e 0 0 0]
        [e 0 0 0 0]

    We can also pass the elements of the last row, starting at the second column::

        sage: matrix.hankel(SR.var('a, b, c, d, e'), SR.var('f, g, h, i'))              # needs sage.symbolic
        [a b c d e]
        [b c d e f]
        [c d e f g]
        [d e f g h]
        [e f g h i]

    A third order Hankel matrix in the integers::

        sage: matrix.hankel([1, 2, 3])
        [1 2 3]
        [2 3 0]
        [3 0 0]

    The second argument allows to customize the last row::

        sage: matrix.hankel([1..3], [7..10])
        [ 1  2  3  7  8]
        [ 2  3  7  8  9]
        [ 3  7  8  9 10]
    """
    m = len(c)
    r = [0] * (m - 1) if r is None else list(r)
    n = len(r)

    def entries(i):
        return c[i] if i < m else r[i - m]
    return matrix(lambda i, j: entries(i + j), nrows=m, ncols=n + 1, base_ring=ring)
