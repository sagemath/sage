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

The web page of ``lrcalc`` is `<http://sites.math.rutgers.edu/~asbuch/lrcalc/>`_.

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
corresponds to computing in the representation ring of gl(3)::

    sage: lrcalc.mult([2,1], [2,1], 3)
    {[2, 2, 2]: 1, [3, 2, 1]: 2, [3, 3]: 1, [4, 1, 1]: 1, [4, 2]: 1}

We can also compute the fusion product, here for sl(3) and level 2::

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
This corresponds to computing in the cohomology ring of Fl(5)::

    sage: lrcalc.mult_schubert([4,2,1,3], [1,4,2,5,3], 5)
    {[4, 5, 1, 3, 2]: 1, [5, 3, 1, 4, 2]: 1, [5, 4, 1, 2, 3]: 1}

List all Littlewood-Richardson tableaux of skew shape `\mu/\nu`; in
this example `\mu=[3,2,1]` and `\nu=[2,1]`. Specifying a third entry
`maxrows` restricts the alphabet to `\{1,2,\ldots,maxrows\}`::

    sage: list(lrcalc.lrskew([3,2,1],[2,1]))
    [[[None, None, 1], [None, 1], [1]], [[None, None, 1], [None, 1], [2]],
    [[None, None, 1], [None, 2], [1]], [[None, None, 1], [None, 2], [3]]]

    sage: list(lrcalc.lrskew([3,2,1],[2,1],maxrows=2))
    [[[None, None, 1], [None, 1], [1]], [[None, None, 1], [None, 1], [2]], [[None, None, 1], [None, 2], [1]]]

.. TODO::

    use this library in the :class:`SymmetricFunctions` code, to
    make it easy to apply it to linear combinations of Schur functions.

.. SEEALSO::

    - :func:`lrcoef`

    - :func:`mult`

    - :func:`coprod`

    - :func:`skew`

    - :func:`lrskew`

    - :func:`mult_schubert`

.. rubric:: Underlying algorithmic in lrcalc

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

from sage.rings.integer cimport Integer
from sage.structure.parent cimport Parent
from sage.combinat.partition import _Partitions
from sage.combinat.permutation import Permutation
from sage.combinat.skew_tableau import SkewTableau


cdef vector* iterable_to_vector(it):
    """
    Return an lrcalc vector (which is a list of integers) from a Python iterable.

    TESTS::

        sage: from sage.libs.lrcalc.lrcalc import test_iterable_to_vector
        sage: x = test_iterable_to_vector(Partition([3,2,1])); x   #indirect doctest
        [3, 2, 1]
    """
    cdef vector* v
    cdef list itr = list(it)
    cdef int n = len(itr)
    cdef int i
    v = v_new(n)
    for i in range(n):
        v.array[i] = int(itr[i])
    return v


cdef list vector_to_list(vector *v):
    """
    Converts a lrcalc vector to Python list.

    TESTS::

        sage: from sage.libs.lrcalc.lrcalc import test_iterable_to_vector
        sage: x = test_iterable_to_vector([]); x         #indirect doctest
        []
    """
    cdef int i, n
    n = v_length(v)
    cdef list result = [None]*n
    for i in range(n):
        result[i] = Integer(v_elem(v, i))
    return result


def test_iterable_to_vector(it):
    """
    A wrapper function for the cdef function ``iterable_to_vector``
    and ``vector_to_list``, to test that they are working correctly.

    EXAMPLES::

        sage: from sage.libs.lrcalc.lrcalc import test_iterable_to_vector
        sage: x = test_iterable_to_vector([3,2,1]); x
        [3, 2, 1]
    """
    cdef vector *v = iterable_to_vector(it)
    result = vector_to_list(v)
    v_free(v)
    return result


cdef skewtab_to_SkewTableau(skewtab *st):
    """
    A wrapper function which transforms the data set ``st`` used in
    ``lrcalc`` to a ``SkewTableau`` in Sage.

    TESTS::

        sage: from sage.libs.lrcalc.lrcalc import test_skewtab_to_SkewTableau
        sage: test_skewtab_to_SkewTableau([],[])
        []
    """
    inner = vector_to_list(st.inner)
    outer = vector_to_list(st.outer)
    return SkewTableau(expr=[[inner[y] for y in range(len(outer))],
                             [[st.matrix[x + y * st.cols] + 1
                                for x in range(inner[y], outer[y])]
                              for y in range(len(outer) - 1, -1, -1)]])


def test_skewtab_to_SkewTableau(outer, inner):
    """
    A wrapper function for the cdef function ``skewtab_to_SkewTableau``
    for testing purposes.

    It constructs the first LR skew tableau of shape ``outer/inner``
    as an ``lrcalc`` ``skewtab``, and converts it to a
    :class:`SkewTableau`.

    EXAMPLES::

        sage: from sage.libs.lrcalc.lrcalc import test_skewtab_to_SkewTableau
        sage: test_skewtab_to_SkewTableau([3,2,1],[])
        [[1, 1, 1], [2, 2], [3]]
        sage: test_skewtab_to_SkewTableau([4,3,2,1],[1,1]).pp()
        .  1  1  1
        .  2  2
        1  3
        2
    """
    cdef vector* o = iterable_to_vector(outer)
    cdef vector* i = iterable_to_vector(inner+[0]*(len(outer)-len(inner)))
    cdef skewtab* st = st_new(o, i, NULL, 0)
    return skewtab_to_SkewTableau(st)


cdef dict sf_hashtab_to_dict(hashtab *ht):
    """
    Return a dictionary representing a Schur function. The keys are
    partitions and the values are integers <class 'sage.rings.integer.Integer'>.

    EXAMPLES::

        sage: from sage.libs.lrcalc.lrcalc import mult
        sage: sorted(mult([1],[1]).items())        #indirect doctest
        [([1, 1], 1), ([2], 1)]
        sage: assert isinstance(mult([1],[1]),dict)#indirect doctest
    """
    cdef hash_itr itr
    cdef dict result = {}
    cdef list p
    hash_first(ht, itr)
    while hash_good(itr):
        p = vector_to_list(<vector*> hash_key(itr))
        result[_Partitions(p)] = Integer(hash_intvalue(itr))
        hash_next(itr)
    return result


cdef dict schubert_hashtab_to_dict(hashtab *ht):
    """
    Return a dictionary corresponding to a Schubert polynomial whose keys
    are permutations and whose values are integers <class 'sage.rings.integer.Integer'>.

    EXAMPLES::

        sage: from sage.libs.lrcalc.lrcalc import mult_schubert
        sage: mult_schubert([3,2,1], [1,2,3])      #indirect doctest
        {[3, 2, 1]: 1}
    """
    cdef hash_itr itr
    cdef dict result = {}
    hash_first(ht, itr)
    while hash_good(itr):
        p = vector_to_list(<vector*> hash_key(itr))
        result[Permutation(p)] = Integer(hash_intvalue(itr))
        hash_next(itr)
    return result


cdef dict vp_hashtab_to_dict(hashtab *ht):
    """
    Return a dictionary corresponding to the coproduct of a Schur function whose keys are
    pairs of partitions and whose values are integers <class 'sage.rings.integer.Integer'>.

    EXAMPLES::

        sage: from sage.libs.lrcalc.lrcalc import coprod
        sage: coprod([1])      #indirect doctest
        {([1], []): 1}
    """
    cdef hash_itr itr
    cdef vecpair* vp
    cdef dict result = {}
    hash_first(ht, itr)
    while hash_good(itr):
        vp = <vecpair*> hash_key(itr)
        p1 = _Partitions(vector_to_list(vp_first(vp)))
        p2 = _Partitions(vector_to_list(vp_second(vp)))
        result[(p1, p2)] = Integer(hash_intvalue(itr))
        hash_next(itr)
    return result


def lrcoef_unsafe(outer, inner1, inner2):
    r"""
    Compute a single Littlewood-Richardson coefficient.

    Return the coefficient of ``outer`` in the product of the Schur
    functions indexed by ``inner1`` and ``inner2``.

    INPUT:

    - ``outer`` -- a partition (weakly decreasing list of non-negative integers).

    - ``inner1`` -- a partition.

    - ``inner2`` -- a partition.

    .. warning::

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
    cdef long long result
    cdef vector *o
    cdef vector *i1
    cdef vector *i2
    o = iterable_to_vector(outer)
    i1 = iterable_to_vector(inner1)
    i2 = iterable_to_vector(inner2)
    result = lrcoef_c(o, i1, i2)
    v_free(o)
    v_free(i1)
    v_free(i2)
    return Integer(result)


def lrcoef(outer, inner1, inner2):
    """
    Compute a single Littlewood-Richardson coefficient.

    Return the coefficient of ``outer`` in the product of the Schur
    functions indexed by ``inner1`` and ``inner2``.

    INPUT:

    - ``outer`` -- a partition (weakly decreasing list of non-negative integers).

    - ``inner1`` -- a partition.

    - ``inner2`` -- a partition.

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


def mult(part1, part2, maxrows=None, level=None, quantum=None):
    r"""
    Compute a product of two Schur functions.

    Return the product of the Schur functions indexed by the
    partitions ``part1`` and ``part2``.

    INPUT:

    - ``part1`` -- a partition

    - ``part2`` -- a partition

    - ``maxrows`` -- (optional) an integer

    - ``level`` -- (optional) an integer

    - ``quantum`` -- (optional) an element of a ring

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
        [([2, 2, 1, 1], 1), ([2, 2, 2], 1), ([3, 1, 1, 1], 1), ([3, 2, 1], 2), ([3, 3], 1), ([4, 1, 1], 1), ([4, 2], 1)]
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

    cdef vector* v1 = iterable_to_vector(part1)
    cdef vector* v2 = iterable_to_vector(part2)
    if maxrows is None:
        maxrows = 0
    cdef hashtab* ht = mult_c(v1, v2, int(maxrows))
    cdef hashtab* tab
    cdef dict result

    if quantum is None:
        if level is not None:
            fusion_reduce_c(ht, int(maxrows), int(level), int(0))
        result = sf_hashtab_to_dict(ht)
        v_free(v1)
        v_free(v2)
        hash_free(ht)
        return result

    # Otherwise do quantum multiplication
    cdef _list *qlist
    cdef dict temp
    qlist = quantum_reduce_c(ht, int(maxrows), int(level))
    # The above call frees the memory associated with ht
    v_free(v1)
    v_free(v2)

    cdef Parent P = quantum.parent()
    result = {}
    for i in range(qlist.length):
        tab = <hashtab*>(qlist.array[i])
        temp = sf_hashtab_to_dict(tab)
        for k in temp:
            result[k] = result.get(k, P.zero()) + quantum**i * temp[k]
        hash_free(tab)
    l_free(qlist)
    return result


def skew(outer, inner, maxrows=0):
    """
    Compute the Schur expansion of a skew Schur function.

    Return a linear combination of partitions representing the Schur
    function of the skew Young diagram ``outer / inner``, consisting
    of boxes in the partition ``outer`` that are not in ``inner``.

    INPUT:

    - ``outer`` -- a partition.

    - ``inner`` -- a partition.

    - ``maxrows`` -- an integer or ``None``.

    If ``maxrows`` is specified, then only partitions with at most
    this number of rows are included in the result.

    EXAMPLES::

        sage: from sage.libs.lrcalc.lrcalc import skew
        sage: sorted(skew([2,1],[1]).items())
        [([1, 1], 1), ([2], 1)]
    """
    cdef vector* v1 = iterable_to_vector(outer)
    cdef vector* v2 = iterable_to_vector(inner)
    cdef hashtab* ht = skew_c(v1, v2, int(maxrows))
    result = sf_hashtab_to_dict(ht)
    v_free(v1)
    v_free(v2)
    hash_free(ht)
    return result


def coprod(part, all=0):
    """
    Compute the coproduct of a Schur function.

    Return a linear combination of pairs of partitions representing
    the coproduct of the Schur function given by the partition
    ``part``.

    INPUT:

    - ``part`` -- a partition.

    - ``all`` -- an integer.

    If ``all`` is non-zero then all terms are included in the result.
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
    cdef vector* v1 = iterable_to_vector(part)
    cdef hashtab* ht = coprod_c(v1, int(all))
    result = vp_hashtab_to_dict(ht)
    v_free(v1)
    hash_free(ht)
    return result


def mult_schubert(w1, w2, rank=0):
    r"""
    Compute a product of two Schubert polynomials.

    Return a linear combination of permutations representing the
    product of the Schubert polynomials indexed by the permutations
    ``w1`` and ``w2``.

    INPUT:

    - ``w1`` -- a permutation.

    - ``w2`` -- a permutation.

    - ``rank`` -- an integer.

    If ``rank`` is non-zero, then only permutations from the symmetric
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
    cdef vector* v1 = iterable_to_vector(w1)
    cdef vector* v2 = iterable_to_vector(w2)
    cdef hashtab* ht = mult_schubert_c(v1, v2, int(rank))
    result = schubert_hashtab_to_dict(ht)
    v_free(v1)
    v_free(v2)
    hash_free(ht)
    return result


def lrskew(outer, inner, weight=None, maxrows=0):
    r"""
    Iterate over the skew LR tableaux of shape ``outer / inner``.

    INPUT:

    - ``outer`` -- a partition

    - ``inner`` -- a partition

    - ``weight`` -- a partition (optional)

    - ``maxrows`` -- an integer (optional)

    OUTPUT: an iterator of :class:`SkewTableau`

    Specifying ``maxrows`` restricts the alphabet to `\{1,2,\ldots,maxrows\}`.

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
    """
    cdef vector* o = iterable_to_vector(outer)
    cdef vector* i = iterable_to_vector(inner + [0]*(len(outer) - len(inner)))
    cdef skewtab* st = st_new(o, i, NULL, int(maxrows))

    if weight is None:
        yield skewtab_to_SkewTableau(st)
        while st_next(st):
            yield skewtab_to_SkewTableau(st)
    else:
        wt = _Partitions(weight)
        r = skewtab_to_SkewTableau(st)
        if r.weight() == wt:
            yield r
        while st_next(st):
            r = skewtab_to_SkewTableau(st)
            if r.weight() == wt:
                yield r
    st_free(st)
