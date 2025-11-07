# sage.doctest: needs sage.modules sage.rings.finite_rings
r"""
Access functions to online databases for coding theory
"""
from sage.misc.lazy_import import lazy_import

# Import the following function so that it is available as
# - sage.coding.databases.self_dual_binary_codes
# - codes.databases.self_dual_binary_codes (which functions as a catalog).
lazy_import('sage.coding.self_dual_codes', 'self_dual_binary_codes')

del lazy_import

# Do not put any global imports here since this module is accessible as
# - sage.coding.databases.<tab>
# - codes.databases.<tab>


def best_linear_code_in_guava(n, k, F):
    r"""
    Return the linear code of length ``n``, dimension ``k`` over field ``F``
    with the maximal minimum distance which is known to the GAP package GUAVA.

    The function uses the tables described in :func:`bounds_on_minimum_distance_in_guava` to
    construct this code. This requires the optional GAP package GUAVA.

    INPUT:

    - ``n`` -- the length of the code to look up

    - ``k`` -- the dimension of the code to look up

    - ``F`` -- the base field of the code to look up

    OUTPUT: a :class:`LinearCode` which is a best linear code of the given parameters known to GUAVA

    EXAMPLES::

        sage: codes.databases.best_linear_code_in_guava(10,5,GF(2))    # long time; optional - gap_package_guava
        [10, 5] linear code over GF(2)
        sage: libgap.LoadPackage('guava')                              # long time; optional - gap_package_guava
        ...
        sage: libgap.BestKnownLinearCode(10,5,libgap.GF(2))            # long time; optional - gap_package_guava
        a linear [10,5,4]2..4 shortened code

    This means that the best possible binary linear code of length 10 and
    dimension 5 is a code with minimum distance 4 and covering radius s somewhere
    between 2 and 4. Use ``bounds_on_minimum_distance_in_guava(10,5,GF(2))``
    for further details.
    """
    from sage.features.gap import GapPackage
    from .linear_code import LinearCode
    GapPackage('guava', spkg='gap_packages').require()
    from sage.libs.gap.libgap import libgap
    libgap.load_package('guava')
    C = libgap.BestKnownLinearCode(n, k, F)
    return LinearCode(C.GeneratorMat()._matrix_(F))


def bounds_on_minimum_distance_in_guava(n, k, F):
    r"""
    Compute a lower and upper bound on the greatest minimum distance of a
    `[n,k]` linear code over the field ``F``.

    This function requires the optional GAP package GUAVA.

    The function returns a GAP record with the two bounds and an explanation for
    each bound. The method ``Display`` can be used to show the explanations.

    The values for the lower and upper bound are obtained from a table
    constructed by Cen Tjhai for GUAVA, derived from the table of
    Brouwer. See http://www.codetables.de/ for the most recent data.
    These tables contain lower and upper bounds for `q=2` (when ``n <= 257``),
    `q=3` (when ``n <= 243``), `q=4` (``n <= 256``). (Current as of
    11 May 2006.) For codes over other fields and for larger word lengths,
    trivial bounds are used.

    INPUT:

    - ``n`` -- the length of the code to look up

    - ``k`` -- the dimension of the code to look up

    - ``F`` -- the base field of the code to look up

    OUTPUT: a GAP record object. See below for an example

    EXAMPLES::

        sage: gap_rec = codes.databases.bounds_on_minimum_distance_in_guava(10,5,GF(2))  # optional - gap_package_guava
        sage: gap_rec.Display()                                                          # optional - gap_package_guava
        rec(
          construction := [ <Operation "ShortenedCode">,
            [ [ <Operation "UUVCode">,
              [ [ <Operation "DualCode">,
              [ [ <Operation "RepetitionCode">, [ 8, 2 ] ] ] ],
              [ <Operation "UUVCode">, [ [ <Operation "DualCode">,
                [ [ <Operation "RepetitionCode">, [ 4, 2 ] ] ] ],
                [ <Operation "RepetitionCode">, [ 4, 2 ] ] ] ] ] ],
            [ 1, 2, 3, 4, 5, 6 ] ] ],
          k := 5,
          lowerBound := 4,
          lowerBoundExplanation := ...
          n := 10,
          q := 2,
          references := rec(
               ),
          upperBound := 4,
          upperBoundExplanation := ... )
    """
    from sage.features.gap import GapPackage
    GapPackage('guava', spkg='gap_packages').require()
    from sage.libs.gap.libgap import libgap
    libgap.load_package('guava')
    return libgap.BoundsMinimumDistance(n, k, F)


def best_linear_code_in_codetables_dot_de(n, k, F, verbose=False):
    r"""
    Return the best linear code and its construction as per the web database
    http://www.codetables.de/

    INPUT:

    - ``n`` -- integer; the length of the code

    - ``k`` -- integer; the dimension of the code

    - ``F`` -- finite field, of order 2, 3, 4, 5, 7, 8, or 9

    - ``verbose`` -- boolean (default: ``False``)

    OUTPUT: an unparsed text explaining the construction of the code

    EXAMPLES::

        sage: L = codes.databases.best_linear_code_in_codetables_dot_de(72, 36, GF(2))    # optional - internet
        sage: print(L)                                                                    # optional - internet
        Construction of a linear code
        [72,36,15] over GF(2):
        [1]:  [73, 36, 16] Cyclic Linear Code over GF(2)
             CyclicCode of length 73 with generating polynomial x^37 + x^36 + x^34 +
        x^33 + x^32 + x^27 + x^25 + x^24 + x^22 + x^21 + x^19 + x^18 + x^15 + x^11 +
        x^10 + x^8 + x^7 + x^5 + x^3 + 1
        [2]:  [72, 36, 15] Linear Code over GF(2)
             Puncturing of [1] at 1
        <BLANKLINE>
        last modified: 2002-03-20

    This function raises an :exc:`IOError` if an error occurs downloading data or
    parsing it. It raises a :exc:`ValueError` if the ``q`` input is invalid.

    AUTHORS:

    - Steven Sivek (2005-11-14)
    - David Joyner (2008-03)
    """
    from urllib.request import urlopen
    from sage.cpython.string import bytes_to_str
    q = F.order()
    if q not in [2, 3, 4, 5, 7, 8, 9]:
        raise ValueError("q (=%s) must be in [2,3,4,5,7,8,9]" % q)
    n = int(n)
    k = int(k)

    param = ("?q=%s&n=%s&k=%s" % (q, n, k)).replace('L', '')

    url = "http://www.codetables.de/" + "BKLC/BKLC.php" + param
    if verbose:
        print("Looking up the bounds at %s" % url)
    with urlopen(url) as f:
        s = f.read()

    s = bytes_to_str(s)
    i = s.find("<PRE>")
    j = s.find("</PRE>")
    if i == -1 or j == -1:
        raise OSError("Error parsing data (missing pre tags).")
    return s[i+5:j].strip()


def self_orthogonal_binary_codes(n, k, b=2, parent=None, BC=None, equal=False,
    in_test=None):
    """
    Return a Python iterator which generates a complete set of
    representatives of all permutation equivalence classes of
    self-orthogonal binary linear codes of length in ``[1..n]`` and
    dimension in ``[1..k]``.

    INPUT:

    - ``n`` -- integer; maximal length

    - ``k`` -- integer; maximal dimension

    - ``b`` -- integer; requires that the generators all have weight divisible
      by ``b`` (if ``b=2``, all self-orthogonal codes are generated, and if
      ``b=4``, all doubly even codes are generated). Must be an even positive
      integer.

    - ``parent`` -- used in recursion (default: ``None``)

    - ``BC`` -- used in recursion (default: ``None``)

    - ``equal`` -- if ``True``, generates only [n, k] codes (default: ``False``)

    - ``in_test`` -- used in recursion (default: ``None``)

    EXAMPLES:

    Generate all self-orthogonal codes of length up to 7 and dimension up
    to 3::

        sage: # needs sage.groups
        sage: for B in codes.databases.self_orthogonal_binary_codes(7,3):
        ....:    print(B)
        [2, 1] linear code over GF(2)
        [4, 2] linear code over GF(2)
        [6, 3] linear code over GF(2)
        [4, 1] linear code over GF(2)
        [6, 2] linear code over GF(2)
        [6, 2] linear code over GF(2)
        [7, 3] linear code over GF(2)
        [6, 1] linear code over GF(2)

    Generate all doubly-even codes of length up to 7 and dimension up
    to 3::

        sage: # needs sage.groups
        sage: for B in codes.databases.self_orthogonal_binary_codes(7,3,4):
        ....:    print(B); print(B.generator_matrix())
        [4, 1] linear code over GF(2)
        [1 1 1 1]
        [6, 2] linear code over GF(2)
        [1 1 1 1 0 0]
        [0 1 0 1 1 1]
        [7, 3] linear code over GF(2)
        [1 0 1 1 0 1 0]
        [0 1 0 1 1 1 0]
        [0 0 1 0 1 1 1]

    Generate all doubly-even codes of length up to 7 and dimension up
    to 2::

        sage: # needs sage.groups
        sage: for B in codes.databases.self_orthogonal_binary_codes(7,2,4):
        ....:    print(B); print(B.generator_matrix())
        [4, 1] linear code over GF(2)
        [1 1 1 1]
        [6, 2] linear code over GF(2)
        [1 1 1 1 0 0]
        [0 1 0 1 1 1]

    Generate all self-orthogonal codes of length equal to 8 and
    dimension equal to 4::

        sage: # needs sage.groups
        sage: for B in codes.databases.self_orthogonal_binary_codes(8, 4, equal=True):
        ....:     print(B); print(B.generator_matrix())
        [8, 4] linear code over GF(2)
        [1 0 0 1 0 0 0 0]
        [0 1 0 0 1 0 0 0]
        [0 0 1 0 0 1 0 0]
        [0 0 0 0 0 0 1 1]
        [8, 4] linear code over GF(2)
        [1 0 0 1 1 0 1 0]
        [0 1 0 1 1 1 0 0]
        [0 0 1 0 1 1 1 0]
        [0 0 0 1 0 1 1 1]

    Since all the codes will be self-orthogonal, b must be divisible by
    2::

        sage: list(codes.databases.self_orthogonal_binary_codes(8, 4, 1, equal=True))
        Traceback (most recent call last):
        ...
        ValueError: b (1) must be a positive even integer.
    """
    from sage.rings.finite_rings.finite_field_constructor import FiniteField
    from sage.matrix.constructor import Matrix

    d = int(b)
    if d != b or d % 2 == 1 or d <= 0:
        raise ValueError("b (%s) must be a positive even integer." % b)
    from .linear_code import LinearCode
    from .binary_code import BinaryCode, BinaryCodeClassifier
    if k < 1 or n < 2:
        return
    if equal:
        in_test = lambda M: (M.ncols() - M.nrows()) <= (n-k)
        out_test = lambda C: (C.dimension() == k) and (C.length() == n)
    else:
        in_test = lambda M: True
        out_test = lambda C: True
    if BC is None:
        BC = BinaryCodeClassifier()
    if parent is None:
        for j in range(d, n+1, d):
            M = Matrix(FiniteField(2), [[1]*j])
            if in_test(M):
                for N in self_orthogonal_binary_codes(n, k, d, M, BC, in_test=in_test):
                    if out_test(N):
                        yield N
    else:
        C = LinearCode(parent)
        if out_test(C):
            yield C
        if k == parent.nrows():
            return
        for nn in range(parent.ncols()+1, n+1):
            if in_test(parent):
                for child in BC.generate_children(BinaryCode(parent), nn, d):
                    for N in self_orthogonal_binary_codes(n, k, d, child, BC, in_test=in_test):
                        if out_test(N):
                            yield N
