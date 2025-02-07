# sage.doctest: needs sage.modules
r"""
Mutually Orthogonal Latin Squares (MOLS)

The main function of this module is :func:`mutually_orthogonal_latin_squares`
and can be can be used to generate MOLS (or check that they exist)::

    sage: MOLS = designs.mutually_orthogonal_latin_squares(4,8)                         # needs sage.schemes

For more information on MOLS, see the :wikipedia:`Wikipedia entry on MOLS
<Graeco-Latin_square#Mutually_orthogonal_Latin_squares>`. If you are only
interested by latin squares, see :mod:`~sage.combinat.matrices.latin`.

The functions defined here are

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`mutually_orthogonal_latin_squares` | Return `k` Mutually Orthogonal `n\times n` Latin Squares.
    :meth:`are_mutually_orthogonal_latin_squares` | Check that the list ``l`` of matrices in are MOLS.
    :meth:`latin_square_product` | Return the product of two (or more) latin squares.
    :meth:`MOLS_table` | Print the MOLS table.

**Table of MOLS**

Sage can produce a table of MOLS similar to the one from the Handbook of
Combinatorial Designs [DesignHandbook]_ (`available here
<http://books.google.fr/books?id=S9FA9rq1BgoC&dq=handbook%20combinatorial%20designs%20MOLS%2010000&pg=PA176>`_).

::

    sage: from sage.combinat.designs.latin_squares import MOLS_table
    sage: MOLS_table(600) # long time
           0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19
        ________________________________________________________________________________
      0| +oo +oo   1   2   3   4   1   6   7   8   2  10   5  12   4   4  15  16   5  18
     20|   4   5   3  22   7  24   4  26   5  28   4  30  31   5   4   5   8  36   4   5
     40|   7  40   5  42   5   6   4  46   8  48   6   5   5  52   5   6   7   7   5  58
     60|   5  60   5   6  63   7   5  66   5   6   6  70   7  72   5   7   6   6   6  78
     80|   9  80   8  82   6   6   6   6   7  88   6   7   6   6   6   6   7  96   6   8
    100|   8 100   6 102   7   7   6 106   6 108   6   6  13 112   6   7   6   8   6   6
    120|   7 120   6   6   6 124   6 126 127   7   6 130   6   7   6   7   7 136   6 138
    140|   6   7   6  10  10   7   6   7   6 148   6 150   7   8   8   7   6 156   7   6
    160|   9   7   6 162   6   7   6 166   7 168   6   8   6 172   6   6  14   9   6 178
    180|   6 180   6   6   7   9   6  10   6   8   6 190   7 192   6   7   6 196   6 198
    200|   7   7   6   7   6   8   6   8  14  11  10 210   6   7   6   7   7   8   6  10
    220|   6  12   6 222  13   8   6 226   6 228   6   7   7 232   6   7   6   7   6 238
    240|   7 240   6 242   6   7   6  12   7   7   6 250   6  12   9   7 255 256   6  12
    260|   6   8   8 262   7   8   7  10   7 268   7 270  15  16   6  13  10 276   6   9
    280|   7 280   6 282   6  12   6   7  15 288   6   6   6 292   6   6   7  10  10  12
    300|   7   7   7   7  15  15   6 306   7   7   7 310   7 312   7  10   7 316   7  10
    320|  15  15   6  16   8  12   6   7   7   9   6 330   7   8   7   6   7 336   6   7
    340|   6  10  10 342   7   7   6 346   6 348   8  12  18 352   6   9   7   9   6 358
    360|   7 360   6   7   7   7   6 366  15  15   7  15   7 372   7  15   7  13   7 378
    380|   7  12   7 382  15  15   7  15   7 388   7  16   7   7   7   7   8 396   7   7
    400|  15 400   7  15  11   8   7  15   8 408   7  13   8  12  10   9  18  15   7 418
    420|   7 420   7  15   7  16   6   7   7   7   6 430  15 432   6  15   6  18   7 438
    440|   7  15   7 442   7  13   7  11  15 448   7  15   7   7   7  15   7 456   7  16
    460|   7 460   7 462  15  15   7 466   8   8   7  15   7  15  10  18   7  15   6 478
    480|  15  15   6  15   8   7   6 486   7  15   6 490   6  16   6   7  15  15   6 498
    500|   7   8   9 502   7  15   6  15   7 508   6  15 511  18   7  15   8  12   8  15
    520|   8 520  10 522  12  15   8  16  15 528   7  15   8  12   7  15   8  15  10  15
    540|  12 540   7  15  18   7   7 546   7   8   7  18   7   7   7   7   7 556   7  12
    560|  15   7   7 562   7   7   6   7   7 568   6 570   7   7  15  22   8 576   7   7
    580|   7   8   7  10   7   8   7 586   7  18  17   7  15 592   8  15   7   7   8 598

Comparison with the results from the Handbook of Combinatorial Designs (2ed)
[DesignHandbook]_::

    sage: MOLS_table(600,compare=True) # long time
            0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19
        ________________________________________________________________________________
      0|                                                           +               +
     20|
     40|
     60|
     80|
    100|
    120|
    140|
    160|
    180|
    200|       -
    220|
    240|
    260|
    280|
    300|
    320|                                                                   -
    340|
    360|   -                   -
    380|                                                       -
    400|
    420|                                       -
    440|
    460|
    480|
    500|       -
    520|
    540|
    560|
    580|

.. TODO::

    Look at [ColDin01]_.

REFERENCES:

.. [Stinson2004] Douglas R. Stinson,
  *Combinatorial designs: construction and analysis*,
  Springer, 2004.

.. [ColDin01] Charles Colbourn, Jeffrey Dinitz,
  *Mutually orthogonal latin squares: a brief survey of constructions*,
  Volume 95, Issues 1-2, Pages 9-48,
  Journal of Statistical Planning and Inference,
  Springer, 1 May 2001.

Functions
---------
"""
from itertools import repeat
from sage.rings.integer import Integer
from sage.categories.sets_cat import EmptySetError
from sage.misc.unknown import Unknown
from sage.arith.misc import is_prime_power
from sage.rings.finite_rings.finite_field_constructor import GF


def are_mutually_orthogonal_latin_squares(l, verbose=False):
    r"""
    Check whether the list of matrices in ``l`` form mutually orthogonal latin
    squares.

    INPUT:

    - ``verbose`` -- if ``True`` then print why the list of matrices provided are
      not mutually orthogonal latin squares

    EXAMPLES::

        sage: from sage.combinat.designs.latin_squares import are_mutually_orthogonal_latin_squares
        sage: m1 = matrix([[0,1,2],[2,0,1],[1,2,0]])
        sage: m2 = matrix([[0,1,2],[1,2,0],[2,0,1]])
        sage: m3 = matrix([[0,1,2],[2,0,1],[1,2,0]])
        sage: are_mutually_orthogonal_latin_squares([m1,m2])
        True
        sage: are_mutually_orthogonal_latin_squares([m1,m3])
        False
        sage: are_mutually_orthogonal_latin_squares([m2,m3])
        True
        sage: are_mutually_orthogonal_latin_squares([m1,m2,m3], verbose=True)
        Squares 0 and 2 are not orthogonal
        False

        sage: m = designs.mutually_orthogonal_latin_squares(7,8)                        # needs sage.schemes
        sage: are_mutually_orthogonal_latin_squares(m)                                  # needs sage.schemes
        True

    TESTS:

    Not a latin square::

        sage: m1 = matrix([[0,1,0],[2,0,1],[1,2,0]])
        sage: m2 = matrix([[0,1,2],[1,2,0],[2,0,1]])
        sage: are_mutually_orthogonal_latin_squares([m1,m2], verbose=True)
        Matrix 0 is not row latin
        False
        sage: m1 = matrix([[0,1,2],[1,0,2],[1,2,0]])
        sage: are_mutually_orthogonal_latin_squares([m1,m2], verbose=True)
        Matrix 0 is not column latin
        False
        sage: m1 = matrix([[0,0,0],[1,1,1],[2,2,2]])
        sage: m2 = matrix([[0,1,2],[0,1,2],[0,1,2]])
        sage: are_mutually_orthogonal_latin_squares([m1,m2])
        False
    """

    if not l:
        raise ValueError("the list must be non empty")

    n = l[0].ncols()
    k = len(l)
    if any(M.ncols() != n or M.nrows() != n for M in l):
        if verbose:
            print("Not all matrices are square matrices of the same dimensions")
        return False

    # Check that all matrices are latin squares
    for i,M in enumerate(l):
        if any(len(set(R)) != n for R in M):
            if verbose:
                print("Matrix {} is not row latin".format(i))
            return False
        if any(len(set(R)) != n for R in zip(*M)):
            if verbose:
                print("Matrix {} is not column latin".format(i))
            return False

    from .designs_pyx import is_orthogonal_array
    return is_orthogonal_array(list(zip(*[[x for R in M for x in R] for M in l])),k,n, verbose=verbose, terminology='MOLS')


def mutually_orthogonal_latin_squares(k, n, partitions=False, check=True):
    r"""
    Return `k` Mutually Orthogonal `n\times n` Latin Squares (MOLS).

    For more information on Mutually Orthogonal Latin Squares, see
    :mod:`~sage.combinat.designs.latin_squares`.

    INPUT:

    - ``k`` -- integer; number of MOLS. If ``k`` is ``None`` it is set to the largest
      value available

    - ``n`` -- integer; size of the latin square

    - ``partitions`` -- boolean; a Latin Square can be seen as 3 partitions of
      the `n^2` cells of the array into `n` sets of size `n`, respectively:

      * The partition of rows
      * The partition of columns
      * The partition of number (cells numbered with 0, cells numbered with 1,
        ...)

      These partitions have the additional property that any two sets from
      different partitions intersect on exactly one element.

      When ``partitions`` is set to ``True``, this function returns a list of `k+2`
      partitions satisfying this intersection property instead of the `k+2` MOLS
      (though the data is exactly the same in both cases).

    - ``check`` -- boolean (default: ``True``); whether to check that output is
      correct before returning it. As this is expected to be useless, you may
      want to disable it whenever you want speed.

    EXAMPLES::

        sage: designs.mutually_orthogonal_latin_squares(4,5)                            # needs sage.schemes
        [
        [0 2 4 1 3]  [0 3 1 4 2]  [0 4 3 2 1]  [0 1 2 3 4]
        [4 1 3 0 2]  [3 1 4 2 0]  [2 1 0 4 3]  [4 0 1 2 3]
        [3 0 2 4 1]  [1 4 2 0 3]  [4 3 2 1 0]  [3 4 0 1 2]
        [2 4 1 3 0]  [4 2 0 3 1]  [1 0 4 3 2]  [2 3 4 0 1]
        [1 3 0 2 4], [2 0 3 1 4], [3 2 1 0 4], [1 2 3 4 0]
        ]

        sage: designs.mutually_orthogonal_latin_squares(3,7)                            # needs sage.schemes
        [
        [0 2 4 6 1 3 5]  [0 3 6 2 5 1 4]  [0 4 1 5 2 6 3]
        [6 1 3 5 0 2 4]  [5 1 4 0 3 6 2]  [4 1 5 2 6 3 0]
        [5 0 2 4 6 1 3]  [3 6 2 5 1 4 0]  [1 5 2 6 3 0 4]
        [4 6 1 3 5 0 2]  [1 4 0 3 6 2 5]  [5 2 6 3 0 4 1]
        [3 5 0 2 4 6 1]  [6 2 5 1 4 0 3]  [2 6 3 0 4 1 5]
        [2 4 6 1 3 5 0]  [4 0 3 6 2 5 1]  [6 3 0 4 1 5 2]
        [1 3 5 0 2 4 6], [2 5 1 4 0 3 6], [3 0 4 1 5 2 6]
        ]

        sage: designs.mutually_orthogonal_latin_squares(2,5,partitions=True)            # needs sage.schemes
        [[[0, 1, 2, 3, 4],
          [5, 6, 7, 8, 9],
          [10, 11, 12, 13, 14],
          [15, 16, 17, 18, 19],
          [20, 21, 22, 23, 24]],
         [[0, 5, 10, 15, 20],
          [1, 6, 11, 16, 21],
          [2, 7, 12, 17, 22],
          [3, 8, 13, 18, 23],
          [4, 9, 14, 19, 24]],
         [[0, 8, 11, 19, 22],
          [3, 6, 14, 17, 20],
          [1, 9, 12, 15, 23],
          [4, 7, 10, 18, 21],
          [2, 5, 13, 16, 24]],
         [[0, 9, 13, 17, 21],
          [2, 6, 10, 19, 23],
          [4, 8, 12, 16, 20],
          [1, 5, 14, 18, 22],
          [3, 7, 11, 15, 24]]]

    What is the maximum number of MOLS of size 8 that Sage knows how to build?::

        sage: designs.orthogonal_arrays.largest_available_k(8)-2                        # needs sage.schemes
        7

    If you only want to know if Sage is able to build a given set of
    MOLS, query the ``orthogonal_arrays.*`` functions::

        sage: designs.orthogonal_arrays.is_available(5+2, 5) # 5 MOLS of order 5
        False
        sage: designs.orthogonal_arrays.is_available(4+2,6)  # 4 MOLS of order 6        # needs sage.schemes
        False

    Sage, however, is not able to prove that the second MOLS do not exist::

        sage: designs.orthogonal_arrays.exists(4+2,6)  # 4 MOLS of order 6              # needs sage.schemes
        Unknown

    If you ask for such a MOLS then you will respectively get an informative
    :exc:`EmptySetError` or :exc:`NotImplementedError`::

        sage: designs.mutually_orthogonal_latin_squares(5, 5)
        Traceback (most recent call last):
        ...
        EmptySetError: there exist at most n-1 MOLS of size n if n>=2
        sage: designs.mutually_orthogonal_latin_squares(4,6)                            # needs sage.schemes
        Traceback (most recent call last):
        ...
        NotImplementedError: I don't know how to build 4 MOLS of order 6

    TESTS:

    The special case `n=1`::

        sage: designs.mutually_orthogonal_latin_squares(3, 1)
        [[0], [0], [0]]

    Wrong input for `k`::

        sage: designs.mutually_orthogonal_latin_squares(None, 1)
        Traceback (most recent call last):
        ...
        TypeError: k must be a positive integer

        sage: designs.mutually_orthogonal_latin_squares(-1, 1)
        Traceback (most recent call last):
        ...
        ValueError: k must be positive

        sage: designs.mutually_orthogonal_latin_squares(2,10)
        [
        [1 8 9 0 2 4 6 3 5 7]  [1 7 6 5 0 9 8 2 3 4]
        [7 2 8 9 0 3 5 4 6 1]  [8 2 1 7 6 0 9 3 4 5]
        [6 1 3 8 9 0 4 5 7 2]  [9 8 3 2 1 7 0 4 5 6]
        [5 7 2 4 8 9 0 6 1 3]  [0 9 8 4 3 2 1 5 6 7]
        [0 6 1 3 5 8 9 7 2 4]  [2 0 9 8 5 4 3 6 7 1]
        [9 0 7 2 4 6 8 1 3 5]  [4 3 0 9 8 6 5 7 1 2]
        [8 9 0 1 3 5 7 2 4 6]  [6 5 4 0 9 8 7 1 2 3]
        [2 3 4 5 6 7 1 8 9 0]  [3 4 5 6 7 1 2 8 0 9]
        [3 4 5 6 7 1 2 0 8 9]  [5 6 7 1 2 3 4 0 9 8]
        [4 5 6 7 1 2 3 9 0 8], [7 1 2 3 4 5 6 9 8 0]
        ]
    """
    from sage.combinat.designs.orthogonal_arrays import orthogonal_array
    from sage.matrix.constructor import Matrix
    from .database import MOLS_constructions

    if k is None:
        raise TypeError('k must be a positive integer')
    try:
        Integer(k)
    except TypeError:
        raise
    if k < 0:
        raise ValueError('k must be positive')

    if n == 1:
        matrices = [Matrix([[0]])] * k

    elif k >= n:
        raise EmptySetError("there exist at most n-1 MOLS of size n if n>=2")

    elif n in MOLS_constructions and k <= MOLS_constructions[n][0]:
        _, construction = MOLS_constructions[n]

        matrices = construction()[:k]

    elif is_prime_power(n):
        F = list(GF(n))

        # We need the first element of the list to be 0
        assert F[0] == 0

        # This dictionary is used to convert from field elements to integers
        conv = {F[i] : i for i in range(n)}
        
        # Make the matrices
        matrices = [Matrix([[conv[F[i] + F[r]*F[j]] for i in range(n)]
                 for j in range(n)]) for r in range(1, k+1)]

    elif orthogonal_array(k + 2, n, existence=True) is not Unknown:
        # Forwarding non-existence results
        if orthogonal_array(k + 2, n, existence=True):
            pass
        else:
            raise EmptySetError("there does not exist {} MOLS of order {}!".format(k, n))

        # make sure that the first two columns are "11, 12, ..., 1n, 21, 22, ..."
        OA = sorted(orthogonal_array(k + 2, n, check=False))

        # We first define matrices as lists of n^2 values
        matrices = [[] for _ in repeat(None, k)]
        for L in OA:
            for i in range(2, k + 2):
                matrices[i-2].append(L[i])

        # The real matrices
        matrices = [[M[i*n:(i+1)*n] for i in range(n)] for M in matrices]
        matrices = [Matrix(M) for M in matrices]

    else:
        raise NotImplementedError("I don't know how to build {} MOLS of order {}".format(k, n))

    if check:
        assert are_mutually_orthogonal_latin_squares(matrices)

    # partitions have been requested but have not been computed yet
    if partitions is True:
        partitions = [[[i*n+j for j in range(n)] for i in range(n)],
                      [[j*n+i for j in range(n)] for i in range(n)]]
        for m in matrices:
            partition = [[] for _ in repeat(None, n)]
            for i in range(n):
                for j in range(n):
                    partition[m[i,j]].append(i*n+j)
            partitions.append(partition)

    if partitions:
        return partitions
    else:
        return matrices


def latin_square_product(M, N, *others):
    r"""
    Return the product of two (or more) latin squares.

    Given two Latin Squares `M,N` of respective sizes `m,n`, the direct product
    `M\times N` of size `mn` is defined by `(M\times
    N)((i_1,i_2),(j_1,j_2))=(M(i_1,j_1),N(i_2,j_2))` where `i_1,j_1\in [m],
    i_2,j_2\in [n]`

    Each pair of values `(i,j)\in [m]\times [n]` is then relabeled to `in+j`.

    This is Lemma 6.25 of [Stinson2004]_.

    INPUT:

    - ``M``, ``N``, ``*others`` -- an arbitrary number of latin squares
      (greater than or equal to 2)

    EXAMPLES::

        sage: from sage.combinat.designs.latin_squares import latin_square_product
        sage: m=designs.mutually_orthogonal_latin_squares(3,4)[0]                       # needs sage.schemes
        sage: latin_square_product(m,m,m)                                               # needs sage.schemes
        64 x 64 sparse matrix over Integer Ring (use the '.str()' method to see the entries)
    """
    from sage.matrix.constructor import Matrix
    m = M.nrows()
    n = N.nrows()

    D = {((i,j),(ii,jj)):(M[i,ii],N[j,jj])
         for i in range(m)
         for ii in range(m)
         for j in range(n)
         for jj in range(n)}

    L = lambda i_j: i_j[0] * n + i_j[1]
    D = {(L(c[0]), L(c[1])): L(v) for c, v in D.items()}
    P = Matrix(D)

    if others:
        return latin_square_product(P, others[0], *others[1:])
    else:
        return P


def MOLS_table(start, stop=None, compare=False, width=None):
    r"""
    Print the MOLS table that Sage can produce.

    INPUT:

    - ``start``, ``stop`` -- integers; print the table of MOLS for value of
      `n` such that ``start<=n<stop``. If only one integer is given as input,
      it is interpreted as the value of ``stop`` with ``start=0`` (same
      behaviour as ``range``).

    - ``compare`` -- boolean; if sets to ``True`` the MOLS displays
      with `+` and `-` entries its difference with the table from the
      Handbook of Combinatorial Designs (2ed).

    - ``width`` -- integer; the width of each column of the table. By default,
      it is computed from range of values determined by the parameters ``start``
      and ``stop``.

    EXAMPLES::

        sage: # needs sage.schemes
        sage: from sage.combinat.designs.latin_squares import MOLS_table
        sage: MOLS_table(100)
               0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19
            ________________________________________________________________________________
          0| +oo +oo   1   2   3   4   1   6   7   8   2  10   5  12   4   4  15  16   5  18
         20|   4   5   3  22   7  24   4  26   5  28   4  30  31   5   4   5   8  36   4   5
         40|   7  40   5  42   5   6   4  46   8  48   6   5   5  52   5   6   7   7   5  58
         60|   5  60   5   6  63   7   5  66   5   6   6  70   7  72   5   7   6   6   6  78
         80|   9  80   8  82   6   6   6   6   7  88   6   7   6   6   6   6   7  96   6   8
        sage: MOLS_table(100, width=4)
                 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19
             ____________________________________________________________________________________________________
           0|  +oo  +oo    1    2    3    4    1    6    7    8    2   10    5   12    4    4   15   16    5   18
          20|    4    5    3   22    7   24    4   26    5   28    4   30   31    5    4    5    8   36    4    5
          40|    7   40    5   42    5    6    4   46    8   48    6    5    5   52    5    6    7    7    5   58
          60|    5   60    5    6   63    7    5   66    5    6    6   70    7   72    5    7    6    6    6   78
          80|    9   80    8   82    6    6    6    6    7   88    6    7    6    6    6    6    7   96    6    8
        sage: MOLS_table(100, compare=True)
               0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19
            ________________________________________________________________________________
          0|                                                           +               +
         20|
         40|
         60|
         80|
        sage: MOLS_table(50, 100, compare=True)
               0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19
            ________________________________________________________________________________
         40|
         60|
         80|
    """
    from .orthogonal_arrays import largest_available_k
    if stop is None:
        start,stop = 0,start
    # make start and stop be congruent to 0 mod 20
    start = start - (start % 20)
    stop = stop-1
    stop = stop + (20-(stop % 20))
    assert start % 20 == 0 and stop % 20 == 0
    if stop <= start:
        return

    # choose an appropriate width (needs to be >= 3 because "+oo" should fit)
    if width is None:
        width = max(3, Integer(stop-1).ndigits(10))

    print(" " * (width + 2) + " ".join("{i:>{width}}".format(i=i,width=width)
                                       for i in range(20)))
    print(" " * (width + 1) + "_" * ((width + 1) * 20), end="")
    for i in range(start,stop):
        if i % 20 == 0:
            print("\n{:>{width}}|".format(i, width=width), end="")
        k = largest_available_k(i)-2
        if compare:
            from . import MOLS_handbook_data
            lower_bound = MOLS_handbook_data.lower_bound(i)
            if i < 2 or lower_bound == k:
                c = ""
            elif lower_bound < k:
                c = "+"
            else:
                c = "-"
        else:
            if i < 2:
                c = "+oo"
            else:
                c = k
        print(' {:>{width}}'.format(c, width=width), end="")
