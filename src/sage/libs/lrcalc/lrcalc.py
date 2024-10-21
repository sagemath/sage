# sage_setup: distribution = sagemath-combinat
# sage.doctest: optional - lrcalc_python, needs sage.combinat
r"""
An interface to Anders Buch's Littlewood-Richardson Calculator ``lrcalc``

The "Littlewood-Richardson Calculator" is a C library for fast
computation of Littlewood-Richardson (LR) coefficients and products of
Schubert polynomials. It handles single LR coefficients, products of
and coproducts of Schur functions, skew Schur functions, and
fusion products. All of the above are achieved by counting LR
(skew)-tableaux (also called Yamanouchi (skew)-tableaux) of
appropriate shape and content by iterating through them.
Additionally, ``lrcalc`` handles products of Schubert polynomials.

The web page of ``lrcalc`` is
`<http://sites.math.rutgers.edu/~asbuch/lrcalc/>`_.

The following describes the Sage interface to this library.

EXAMPLES::

    sage: import sage.libs.lrcalc.lrcalc as lrcalc

Compute a single Littlewood-Richardson coefficient::

    sage: lrcalc.lrcoef([3,2,1],[2,1],[2,1])
    2

Compute a product of Schur functions; return the coefficients in the
Schur expansion::

    sage: lrcalc.mult([2,1], [2,1])
    {[2, 2, 1, 1]: 1,
     [2, 2, 2]: 1,
     [3, 1, 1, 1]: 1,
     [3, 2, 1]: 2,
     [3, 3]: 1,
     [4, 1, 1]: 1,
     [4, 2]: 1}

Same product, but include only partitions with at most 3 rows.  This
corresponds to computing in the representation ring of `\mathfrak{gl}(3)`::

    sage: lrcalc.mult([2,1], [2,1], 3)
    {[2, 2, 2]: 1, [3, 2, 1]: 2, [3, 3]: 1, [4, 1, 1]: 1, [4, 2]: 1}

We can also compute the fusion product, here for `\mathfrak{sl}(3)`
and level 2::

    sage: lrcalc.mult([3,2,1], [3,2,1], 3,2)
    {[4, 4, 4]: 1, [5, 4, 3]: 1}

Compute the expansion of a skew Schur function::

    sage: lrcalc.skew([3,2,1],[2,1])
    {[1, 1, 1]: 1, [2, 1]: 2, [3]: 1}

Compute the coproduct of a Schur function::

    sage: lrcalc.coprod([3,2,1])
    {([1, 1, 1], [2, 1]): 1,
     ([2, 1], [2, 1]): 2,
     ([2, 1], [3]): 1,
     ([2, 1, 1], [1, 1]): 1,
     ([2, 1, 1], [2]): 1,
     ([2, 2], [1, 1]): 1,
     ([2, 2], [2]): 1,
     ([2, 2, 1], [1]): 1,
     ([3, 1], [1, 1]): 1,
     ([3, 1], [2]): 1,
     ([3, 1, 1], [1]): 1,
     ([3, 2], [1]): 1,
     ([3, 2, 1], []): 1}

Multiply two Schubert polynomials::

    sage: lrcalc.mult_schubert([4,2,1,3], [1,4,2,5,3])
    {[4, 5, 1, 3, 2]: 1,
     [5, 3, 1, 4, 2]: 1,
     [5, 4, 1, 2, 3]: 1,
     [6, 2, 1, 4, 3, 5]: 1}

Same product, but include only permutations of 5 elements in the result.
This corresponds to computing in the cohomology ring of `Fl(5)`::

    sage: lrcalc.mult_schubert([4,2,1,3], [1,4,2,5,3], 5)
    {[4, 5, 1, 3, 2]: 1, [5, 3, 1, 4, 2]: 1, [5, 4, 1, 2, 3]: 1}

List all Littlewood-Richardson tableaux of skew shape `\mu/\nu`; in
this example `\mu=[3,2,1]` and `\nu=[2,1]`. Specifying a third entry
`M' = ``maxrows`` restricts the alphabet to `\{1,2,\ldots,M\}`::

    sage: list(lrcalc.lrskew([3,2,1],[2,1]))
    [[[None, None, 1], [None, 1], [1]], [[None, None, 1], [None, 1], [2]],
    [[None, None, 1], [None, 2], [1]], [[None, None, 1], [None, 2], [3]]]

    sage: list(lrcalc.lrskew([3,2,1],[2,1],maxrows=2))
    [[[None, None, 1], [None, 1], [1]], [[None, None, 1], [None, 1], [2]],
     [[None, None, 1], [None, 2], [1]]]

.. TODO::

    Use this library in the :class:`SymmetricFunctions` code, to
    make it easy to apply it to linear combinations of Schur functions.

.. SEEALSO::

    - :func:`lrcoef`
    - :func:`mult`
    - :func:`coprod`
    - :func:`skew`
    - :func:`lrskew`
    - :func:`mult_schubert`

.. RUBRIC:: Underlying algorithmic in lrcalc

Here is some additional information regarding the main low-level
C-functions in `lrcalc`. Given two partitions ``outer`` and ``inner``
with ``inner`` contained in ``outer``, the function::

    skewtab *st_new(vector *outer, vector *inner, vector *conts, int maxrows)

constructs and returns the (lexicographically) first LR skew tableau
of shape ``outer / inner``. Further restrictions can be imposed using
``conts`` and ``maxrows``.

Namely, the integer ``maxrows`` is a bound on the integers that can be
put in the tableau.  The name is chosen because this will limit the
partitions in the output of :func:`skew` or :func:`mult` to partitions
with at most this number of rows.

The vector ``conts`` is the content of an empty tableau(!!). More
precisely, this vector is added to the usual content of a tableau
whenever the content is needed.  This affects which tableaux are
considered LR tableaux (see :func:`mult` below).  ``conts`` may also
be the ``NULL`` pointer, in which case nothing is added.

The other function::

    int *st_next(skewtab *st)

computes in place the (lexicographically) next skew tableau with the
same constraints, or returns 0 if ``st`` is the last one.

For a first example, see the :func:`skew` function code in the
``lrcalc`` source code. We want to compute a skew Schur function, so
create a skew LR tableau of the appropriate shape with ``st_new``
(with ``conts = NULL``), then iterate through all the LR tableaux with
``st_next()``. For each skew tableau, we use that ``st->conts`` is the
content of the skew tableau, find this shape in the ``res`` hash table
and add one to the value.

For a second example, see ``mult(vector *sh1, vector *sh2, maxrows)``.
Here we call ``st_new()`` with the shape ``sh1 / (0)`` and use ``sh2``
as the ``conts`` argument.  The effect of using ``sh2`` in this way is
that ``st_next`` will iterate through semistandard tableaux `T` of
shape ``sh1`` such that the following tableau::

         111111
         22222    <--- minimal tableau of shape sh2
         333
    *****
    **T**
    ****
    **

is a LR skew tableau, and ``st->conts`` contains the content of the
combined tableaux.

More generally, ``st_new(outer, inner, conts, maxrows)`` and
``st_next`` can be used to compute the Schur expansion of the product
``S_{outer/inner} * S_conts``, restricted to partitions with at most
``maxrows`` rows.

AUTHORS:

- Mike Hansen (2010): core of the interface

- Anne Schilling, Nicolas M. Thiéry, and Anders Buch (2011): fusion
  product, iterating through LR tableaux, finalization, documentation
"""
# ****************************************************************************
#  Copyright (C) 2010 Mike Hansen <mhansen@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.combinat.partition import _Partitions
from sage.combinat.permutation import Permutation
from sage.combinat.skew_tableau import SemistandardSkewTableaux
from sage.combinat.skew_partition import SkewPartition
from sage.rings.integer import Integer
import lrcalc


def _lrcalc_dict_to_sage(result) -> dict:
    r"""
    Translate from lrcalc output format to Sage expected format.

    TESTS::

        sage: from sage.libs.lrcalc.lrcalc import mult
        sage: mult([2,1],[3,2,1],3) # indirect doctest
        {[3, 3, 3]: 1, [4, 3, 2]: 2, [4, 4, 1]: 1, [5, 2, 2]: 1, [5, 3, 1]: 1}
    """
    return {_Partitions.element_class(_Partitions, [Integer(p) for p in la]):
            Integer(k) for la, k in result.items()}


def lrcoef_unsafe(outer, inner1, inner2):
    r"""
    Compute a single Littlewood-Richardson coefficient.

    Return the coefficient of ``outer`` in the product of the Schur
    functions indexed by ``inner1`` and ``inner2``.

    INPUT:

    - ``outer`` -- a partition (weakly decreasing list of nonnegative integers)
    - ``inner1`` -- a partition
    - ``inner2`` -- a partition

    .. WARNING::

       This function does not do any check on its input.  If you want
       to use a safer version, use :func:`lrcoef`.

    EXAMPLES::

        sage: from sage.libs.lrcalc.lrcalc import lrcoef_unsafe
        sage: lrcoef_unsafe([3,2,1], [2,1], [2,1])
        2
        sage: lrcoef_unsafe([3,3], [2,1], [2,1])
        1
        sage: lrcoef_unsafe([2,1,1,1,1], [2,1], [2,1])
        0
    """
    return Integer(lrcalc.lrcoef(outer, inner1, inner2))


def lrcoef(outer, inner1, inner2):
    """
    Compute a single Littlewood-Richardson coefficient.

    Return the coefficient of ``outer`` in the product of the Schur
    functions indexed by ``inner1`` and ``inner2``.

    INPUT:

    - ``outer`` -- a partition (weakly decreasing list of nonnegative integers)
    - ``inner1`` -- a partition
    - ``inner2`` -- a partition

    .. NOTE::

       This function converts its inputs into :func:`Partition`'s.  If
       you don't need these checks and your inputs are valid, then you
       can use :func:`lrcoef_unsafe`.

    EXAMPLES::

        sage: from sage.libs.lrcalc.lrcalc import lrcoef
        sage: lrcoef([3,2,1], [2,1], [2,1])
        2
        sage: lrcoef([3,3], [2,1], [2,1])
        1
        sage: lrcoef([2,1,1,1,1], [2,1], [2,1])
        0
    """
    return lrcoef_unsafe(_Partitions(outer), _Partitions(inner1), _Partitions(inner2))


def mult(part1, part2, maxrows=None, level=None, quantum=None) -> dict:
    r"""
    Compute a product of two Schur functions.

    Return the product of the Schur functions indexed by the
    partitions ``part1`` and ``part2``.

    INPUT:

    - ``part1`` -- a partition
    - ``part2`` -- a partition
    - ``maxrows`` -- integer (optional)
    - ``level`` -- integer (optional)
    - ``quantum`` -- an element of a ring (optional)

    If ``maxrows`` is specified, then only partitions with at most
    this number of rows are included in the result.

    If both ``maxrows`` and ``level`` are specified, then the function
    calculates the fusion product for `\mathfrak{sl}(\mathrm{maxrows})`
    of the given level.

    If ``quantum`` is set, then this returns the product in the quantum
    cohomology ring of the Grassmannian. In particular, both ``maxrows``
    and ``level`` need to be specified.

    EXAMPLES::

        sage: from sage.libs.lrcalc.lrcalc import mult
        sage: mult([2],[])
        {[2]: 1}
        sage: sorted(mult([2],[2]).items())
        [([2, 2], 1), ([3, 1], 1), ([4], 1)]
        sage: sorted(mult([2,1],[2,1]).items())
        [([2, 2, 1, 1], 1), ([2, 2, 2], 1), ([3, 1, 1, 1], 1),
         ([3, 2, 1], 2), ([3, 3], 1), ([4, 1, 1], 1), ([4, 2], 1)]
        sage: sorted(mult([2,1],[2,1],maxrows=2).items())
        [([3, 3], 1), ([4, 2], 1)]
        sage: mult([2,1],[3,2,1],3)
        {[3, 3, 3]: 1, [4, 3, 2]: 2, [4, 4, 1]: 1, [5, 2, 2]: 1, [5, 3, 1]: 1}
        sage: mult([2,1],[2,1],3,3)
        {[2, 2, 2]: 1, [3, 2, 1]: 2, [3, 3]: 1, [4, 1, 1]: 1}
        sage: mult([2,1],[2,1],None,3)
        Traceback (most recent call last):
        ...
        ValueError: maxrows needs to be specified if you specify the level

     The quantum product::

        sage: q = polygen(QQ, 'q')
        sage: sorted(mult([1],[2,1], 2, 2, quantum=q).items())
        [([], q), ([2, 2], 1)]
        sage: sorted(mult([2,1],[2,1], 2, 2, quantum=q).items())
        [([1, 1], q), ([2], q)]

        sage: mult([2,1],[2,1], quantum=q)
        Traceback (most recent call last):
        ...
        ValueError: missing parameters maxrows or level
    """
    if maxrows is None and level is not None:
        raise ValueError('maxrows needs to be specified if you specify'
                         ' the level')
    if quantum is not None and (level is None or maxrows is None):
        raise ValueError('missing parameters maxrows or level')

    if quantum is None:
        if level is not None:
            return _lrcalc_dict_to_sage(lrcalc.mult_fusion(part1, part2, maxrows, level))
        if maxrows is None:
            maxrows = -1
        return _lrcalc_dict_to_sage(lrcalc.mult(part1, part2, maxrows))

    # Otherwise do quantum multiplication
    result = lrcalc.mult_quantum(part1, part2, maxrows, level, degrees=True)
    P = quantum.parent()
    output = {}
    for i, k in result.items():
        la = _Partitions(i[0])
        output[la] = output.get(la, P.zero()) + k * quantum**(i[1])
    return output


def skew(outer, inner, maxrows=-1) -> dict:
    """
    Compute the Schur expansion of a skew Schur function.

    Return a linear combination of partitions representing the Schur
    function of the skew Young diagram ``outer / inner``, consisting
    of boxes in the partition ``outer`` that are not in ``inner``.

    INPUT:

    - ``outer`` -- a partition
    - ``inner`` -- a partition
    - ``maxrows`` -- integer or ``None``

    If ``maxrows`` is specified, then only partitions with at most
    this number of rows are included in the result.

    EXAMPLES::

        sage: from sage.libs.lrcalc.lrcalc import skew
        sage: sorted(skew([2,1],[1]).items())
        [([1, 1], 1), ([2], 1)]
    """
    return _lrcalc_dict_to_sage(lrcalc.skew(outer, inner, maxrows))


def coprod(part, all=0) -> dict:
    """
    Compute the coproduct of a Schur function.

    Return a linear combination of pairs of partitions representing
    the coproduct of the Schur function given by the partition
    ``part``.

    INPUT:

    - ``part`` -- a partition
    - ``all`` -- integer

    If ``all`` is nonzero then all terms are included in the result.
    If ``all`` is zero, then only pairs of partitions ``(part1,
    part2)`` for which the weight of ``part1`` is greater than or
    equal to the weight of ``part2`` are included; the rest of the
    coefficients are redundant because Littlewood-Richardson
    coefficients are symmetric.

    EXAMPLES::

        sage: from sage.libs.lrcalc.lrcalc import coprod
        sage: sorted(coprod([2,1]).items())
        [(([1, 1], [1]), 1), (([2], [1]), 1), (([2, 1], []), 1)]
    """
    result = lrcalc.coprod(part, all)
    return {tuple([_Partitions.element_class(_Partitions,
                                             [Integer(p) for p in mu])
                   for mu in la]): Integer(k)
            for la, k in result.items()}


def mult_schubert(w1, w2, rank=0) -> dict:
    r"""
    Compute a product of two Schubert polynomials.

    Return a linear combination of permutations representing the
    product of the Schubert polynomials indexed by the permutations
    ``w1`` and ``w2``.

    INPUT:

    - ``w1`` -- a permutation
    - ``w2`` -- a permutation
    - ``rank`` -- integer

    If ``rank`` is nonzero, then only permutations from the symmetric
    group `S(\mathrm{rank})` are included in the result.

    EXAMPLES::

        sage: from sage.libs.lrcalc.lrcalc import mult_schubert
        sage: result = mult_schubert([3, 1, 5, 2, 4], [3, 5, 2, 1, 4])
        sage: sorted(result.items())
        [([5, 4, 6, 1, 2, 3], 1), ([5, 6, 3, 1, 2, 4], 1),
         ([5, 7, 2, 1, 3, 4, 6], 1), ([6, 3, 5, 1, 2, 4], 1),
         ([6, 4, 3, 1, 2, 5], 1), ([6, 5, 2, 1, 3, 4], 1),
         ([7, 3, 4, 1, 2, 5, 6], 1), ([7, 4, 2, 1, 3, 5, 6], 1)]
    """
    result = lrcalc.schubmult(w1, w2, rank)
    return {Permutation(list(la)): Integer(k) for la, k in result.items()}


def lrskew(outer, inner, weight=None, maxrows=-1):
    r"""
    Iterate over the skew LR tableaux of shape ``outer / inner``.

    INPUT:

    - ``outer`` -- a partition
    - ``inner`` -- a partition
    - ``weight`` -- a partition (optional)
    - ``maxrows`` -- positive integer (optional)

    OUTPUT: an iterator of :class:`SkewTableau`

    Specifying ``maxrows`` = `M` restricts the alphabet to `\{1,2,\ldots,M\}`.

    Specifying ``weight`` returns only those tableaux of given content/weight.

    EXAMPLES::

        sage: from sage.libs.lrcalc.lrcalc import lrskew
        sage: for st in lrskew([3,2,1],[2]):
        ....:     st.pp()
        .  .  1
        1  1
        2
        .  .  1
        1  2
        2
        .  .  1
        1  2
        3

        sage: for st in lrskew([3,2,1],[2], maxrows=2):
        ....:     st.pp()
        .  .  1
        1  1
        2
        .  .  1
        1  2
        2

        sage: list(lrskew([3,2,1],[2], weight=[3,1]))
        [[[None, None, 1], [1, 1], [2]]]

    TESTS::

        sage: from sage.libs.lrcalc.lrcalc import lrskew
        sage: list(lrskew([3,2,1],[2], weight=[]))
        []
        sage: list(lrskew([3,2,1],[2], weight=[0]))
        []
        sage: list(lrskew([3,2,1],[3,2,1], weight=[]))
        [[[None, None, None], [None, None], [None]]]
        sage: list(lrskew([3,2,1],[3,2,1], weight=[0]))
        [[[None, None, None], [None, None], [None]]]
        sage: list(lrskew([3,2,1],[3,2,1], weight=[1]))
        []
    """
    iterator = lrcalc.lr_iterator(outer, inner, maxrows)
    shape = SkewPartition([outer, inner])

    if weight is None:
        ST = SemistandardSkewTableaux(shape)
        for data in iterator:
            yield ST.from_shape_and_word(shape, [i + 1 for i in data])
    else:
        wt = _Partitions(weight)
        ST = SemistandardSkewTableaux(shape, wt)
        m = len(wt)
        for data in iterator:
            w = [0] * m
            for j in data:
                if j >= m:
                    # We know they are not equal, so make the check below quick
                    w = None
                    break
                w[j] += 1
            if w == wt:
                yield ST.from_shape_and_word(shape, [i + 1 for i in data])
