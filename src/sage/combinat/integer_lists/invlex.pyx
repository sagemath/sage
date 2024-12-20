r"""
Enumerated set of lists of integers with constraints, in inverse lexicographic order

- :class:`IntegerListsLex`: the enumerated set of lists of nonnegative
  integers with specified constraints, in inverse lexicographic order.

- :class:`IntegerListsBackend_invlex`: Cython back-end for
  :class:`IntegerListsLex`.

HISTORY:

This generic tool was originally written by Hivert and Thiery in
MuPAD-Combinat in 2000 and ported over to Sage by Mike Hansen in
2007. It was then completely rewritten in 2015 by Gillespie,
Schilling, and Thiery, with the help of many, to deal with
limitations and lack of robustness w.r.t. input.
"""

# ****************************************************************************
#       Copyright (C) 2015 Bryan Gillespie <Brg008@gmail.com>
#                          Nicolas M. Thiery <nthiery at users.sf.net>
#                          Anne Schilling <anne@math.ucdavis.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import builtins

from sage.misc.classcall_metaclass import ClasscallMetaclass, typecall
from sage.misc.cachefunc import cached_method
from sage.combinat.integer_lists.base cimport IntegerListsBackend
from sage.combinat.integer_lists.lists import IntegerLists
from sage.combinat.integer_lists.base import Infinity


class IntegerListsLex(IntegerLists, metaclass=ClasscallMetaclass):
    r"""
    Lists of nonnegative integers with constraints, in inverse
    lexicographic order.

    An *integer list* is a list `l` of nonnegative integers, its *parts*. The
    slope (at position `i`) is the difference ``l[i+1]-l[i]`` between two
    consecutive parts.

    This class allows to construct the set `S` of all integer lists
    `l` satisfying specified bounds on the sum, the length, the slope,
    and the individual parts, enumerated in *inverse* lexicographic
    order, that is from largest to smallest in lexicographic
    order. Note that, to admit such an enumeration, `S` is almost
    necessarily finite (see :ref:`IntegerListsLex_finiteness`).

    The main purpose is to provide a generic iteration engine for all the
    enumerated sets like :class:`Partitions`, :class:`Compositions`,
    :class:`IntegerVectors`. It can also be used to generate many other
    combinatorial objects like Dyck paths, Motzkin paths, etc. Mathematically
    speaking, this is a special case of set of integral points of a polytope (or
    union thereof, when the length is not fixed).

    INPUT:

    - ``min_sum`` -- nonnegative integer (default: 0);
      a lower bound on ``sum(l)``

    - ``max_sum`` -- nonnegative integer or `\infty` (default: `\infty`);
      an upper bound on ``sum(l)``

    - ``n`` -- nonnegative integer (optional); if specified, this
      overrides ``min_sum`` and ``max_sum``

    - ``min_length`` -- nonnegative integer (default: `0`); a lower
      bound on ``len(l)``

    - ``max_length`` -- nonnegative integer or `\infty` (default:
      `\infty`); an upper bound on ``len(l)``

    - ``length`` -- integer (optional); overrides ``min_length``
      and ``max_length`` if specified

    - ``min_part`` -- nonnegative integer; a lower bounds on all
      parts: ``min_part <= l[i]`` for ``0 <= i < len(l)``

    - ``floor`` -- list of nonnegative integers or a function; lower
      bounds on the individual parts `l[i]`

      If ``floor`` is a list of integers, then ``floor<=l[i]`` for ``0
      <= i < min(len(l), len(floor)``. Similarly, if ``floor`` is a
      function, then ``floor(i) <= l[i]`` for ``0 <= i < len(l)``.

    - ``max_part`` -- nonnegative integer or `\infty`; an upper
      bound on all parts: ``l[i] <= max_part`` for ``0 <= i < len(l)``

    - ``ceiling`` -- upper bounds on the individual parts ``l[i]``;
      this takes the same type of input as ``floor``, except that
      `\infty` is allowed in addition to integers, and the default
      value is `\infty`.

    - ``min_slope`` -- integer or `-\infty` (default: `-\infty`);
      a lower bound on the slope between consecutive parts:
      ``min_slope <= l[i+1]-l[i]`` for ``0 <= i < len(l)-1``

    - ``max_slope`` -- integer or `+\infty` (defaults: `+\infty`);
      an upper bound on the slope between consecutive parts:
      ``l[i+1]-l[i] <= max_slope`` for ``0 <= i < len(l)-1``

    - ``category`` -- a category (default: :class:`FiniteEnumeratedSets`)

    - ``check`` -- boolean (default: ``True``); whether to display the
      warnings raised when functions are given as input to ``floor``
      or ``ceiling`` and the errors raised when there is no proper
      enumeration.

    - ``name`` -- string or ``None`` (default: ``None``); if set,
      this will be passed down to :meth:`Parent.rename` to specify the
      name of ``self``. It is recommended to use rename method directly
      because this feature may become deprecated.

    - ``element_constructor`` -- a function (or callable) that creates
      elements of ``self`` from a list. See also :class:`Parent`

    - ``element_class`` -- a class for the elements of ``self``
      (default: `ClonableArray`). This merely sets the attribute
      ``self.Element``. See the examples for details.

    .. NOTE::

        When several lists satisfying the constraints differ only by
        trailing zeroes, only the shortest one is enumerated (and
        therefore counted). The others are still considered valid.
        See the examples below.

        This feature is questionable. It is recommended not to rely on
        it, as it may eventually be discontinued.

    EXAMPLES:

    We create the enumerated set of all lists of nonnegative integers
    of length `3` and sum `2`::

        sage: C = IntegerListsLex(2, length=3)
        sage: C
        Integer lists of sum 2 satisfying certain constraints
        sage: C.cardinality()
        6
        sage: [p for p in C]
        [[2, 0, 0], [1, 1, 0], [1, 0, 1], [0, 2, 0], [0, 1, 1], [0, 0, 2]]

        sage: [2, 0, 0] in C
        True
        sage: [2, 0, 1] in C
        False
        sage: "a" in C
        False
        sage: ["a"] in C
        False
        sage: C.first()
        [2, 0, 0]

    One can specify lower and upper bounds on each part::

        sage: list(IntegerListsLex(5, length=3, floor=[1,2,0], ceiling=[3,2,3]))
        [[3, 2, 0], [2, 2, 1], [1, 2, 2]]

    When the length is fixed as above, one can also use
    :class:`IntegerVectors`::

        sage: IntegerVectors(2,3).list()
        [[2, 0, 0], [1, 1, 0], [1, 0, 1], [0, 2, 0], [0, 1, 1], [0, 0, 2]]

    Using the slope condition, one can generate integer partitions
    (but see :class:`Partitions`)::

        sage: list(IntegerListsLex(4, max_slope=0))
        [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]

    The following is the list of all partitions of `7` with parts at least `2`::

        sage: list(IntegerListsLex(7, max_slope=0, min_part=2))
        [[7], [5, 2], [4, 3], [3, 2, 2]]


    .. RUBRIC:: floor and ceiling conditions

    Next we list all partitions of `5` of length at most `3` which are
    bounded below by ``[2,1,1]``::

        sage: list(IntegerListsLex(5, max_slope=0, max_length=3, floor=[2,1,1]))
        [[5], [4, 1], [3, 2], [3, 1, 1], [2, 2, 1]]

    Note that ``[5]`` is considered valid, because the floor
    constraints only apply to existing positions in the list. To
    obtain instead the partitions containing ``[2,1,1]``, one needs to
    use ``min_length`` or ``length``::

        sage: list(IntegerListsLex(5, max_slope=0, length=3, floor=[2,1,1]))
        [[3, 1, 1], [2, 2, 1]]

    Here is the list of all partitions of `5` which are contained in
    ``[3,2,2]``::

        sage: list(IntegerListsLex(5, max_slope=0, max_length=3, ceiling=[3,2,2]))
        [[3, 2], [3, 1, 1], [2, 2, 1]]

    This is the list of all compositions of `4` (but see :class:`Compositions`)::

        sage: list(IntegerListsLex(4, min_part=1))
        [[4], [3, 1], [2, 2], [2, 1, 1], [1, 3], [1, 2, 1], [1, 1, 2], [1, 1, 1, 1]]

    This is the list of all integer vectors of sum `4` and length `3`::

        sage: list(IntegerListsLex(4, length=3))
        [[4, 0, 0], [3, 1, 0], [3, 0, 1], [2, 2, 0], [2, 1, 1],
         [2, 0, 2], [1, 3, 0], [1, 2, 1], [1, 1, 2], [1, 0, 3],
         [0, 4, 0], [0, 3, 1], [0, 2, 2], [0, 1, 3], [0, 0, 4]]

    For whatever it is worth, the ``floor`` and ``min_part``
    constraints can be combined::

        sage: L = IntegerListsLex(5, floor=[2,0,2], min_part=1)
        sage: L.list()
        [[5], [4, 1], [3, 2], [2, 3], [2, 1, 2]]

    This is achieved by updating the floor upon constructing ``L``::

        sage: [L.floor(i) for i in range(5)]
        [2, 1, 2, 1, 1]

    Similarly, the ``ceiling`` and ``max_part`` constraints can be
    combined::

        sage: L = IntegerListsLex(4, ceiling=[2,3,1], max_part=2, length=3)
        sage: L.list()
        [[2, 2, 0], [2, 1, 1], [1, 2, 1]]
        sage: [L.ceiling(i) for i in range(5)]
        [2, 2, 1, 2, 2]


    This can be used to generate Motzkin words (see
    :wikipedia:`Motzkin_number`)::

        sage: def motzkin_words(n):
        ....:     return IntegerListsLex(length=n+1, min_slope=-1, max_slope=1,
        ....:                ceiling=[0]+[+oo for i in range(n-1)]+[0])
        sage: motzkin_words(4).list()
        [[0, 1, 2, 1, 0],
         [0, 1, 1, 1, 0],
         [0, 1, 1, 0, 0],
         [0, 1, 0, 1, 0],
         [0, 1, 0, 0, 0],
         [0, 0, 1, 1, 0],
         [0, 0, 1, 0, 0],
         [0, 0, 0, 1, 0],
         [0, 0, 0, 0, 0]]
        sage: [motzkin_words(n).cardinality() for n in range(8)]
        [1, 1, 2, 4, 9, 21, 51, 127]
        sage: oeis(_)                  # optional -- internet
        0: A001006: Motzkin numbers: number of ways of drawing any number
        of nonintersecting chords joining n (labeled) points on a circle.
        1: ...
        2: ...

    or Dyck words (see also :class:`DyckWords`), through the bijection
    with paths from `(0,0)` to `(n,n)` with left and up steps that remain
    below the diagonal::

        sage: def dyck_words(n):
        ....:     return IntegerListsLex(length=n, ceiling=list(range(n+1)), min_slope=0)
        sage: [dyck_words(n).cardinality() for n in range(8)]
        [1, 1, 2, 5, 14, 42, 132, 429]
        sage: dyck_words(3).list()
        [[0, 1, 2], [0, 1, 1], [0, 0, 2], [0, 0, 1], [0, 0, 0]]


    .. _IntegerListsLex_finiteness:

    .. RUBRIC:: On finiteness and inverse lexicographic enumeration

    The set of all lists of integers cannot be enumerated in inverse
    lexicographic order, since there is no largest list (take `[n]`
    for `n` as large as desired)::

        sage: IntegerListsLex().first()
        Traceback (most recent call last):
        ...
        ValueError: could not prove that the specified constraints yield a finite set

    Here is a variant which could be enumerated in lexicographic order
    but not in inverse lexicographic order::

        sage: L = IntegerListsLex(length=2, ceiling=[Infinity, 0], floor=[0,1])
        sage: for l in L: print(l)
        Traceback (most recent call last):
        ...
        ValueError: infinite upper bound for values of m

    Even when the sum is specified, it is not necessarily possible to
    enumerate *all* elements in inverse lexicographic order. In the
    following example, the list ``[1, 1, 1]`` will never appear in the
    enumeration::

        sage: IntegerListsLex(3).first()
        Traceback (most recent call last):
        ...
        ValueError: could not prove that the specified constraints yield a finite set

    If one wants to proceed anyway, one can sign a waiver by setting
    ``check=False`` (again, be warned that some valid lists may never appear)::

        sage: L = IntegerListsLex(3, check=False)
        sage: it = iter(L)
        sage: [next(it) for i in range(6)]
        [[3], [2, 1], [2, 0, 1], [2, 0, 0, 1], [2, 0, 0, 0, 1], [2, 0, 0, 0, 0, 1]]

    In fact, being inverse lexicographically enumerable is almost
    equivalent to being finite. The only infinity that can occur would
    be from a tail of numbers `0,1` as in the previous example, where
    the `1` moves further and further to the right. If there is any
    list that is inverse lexicographically smaller than such a
    configuration, the iterator would not reach it and hence would not
    be considered iterable. Given that the infinite cases are very
    specific, at this point only the finite cases are supported
    (without signing the waiver).

    The finiteness detection is not complete yet, so some finite cases
    may not be supported either, at least not without disabling the
    checks. Practical examples of such are welcome.

    .. RUBRIC:: On trailing zeroes, and their caveats

    As mentioned above, when several lists satisfying the constraints
    differ only by trailing zeroes, only the shortest one is listed::

        sage: L = IntegerListsLex(max_length=4, max_part=1)
        sage: L.list()
        [[1, 1, 1, 1],
         [1, 1, 1],
         [1, 1, 0, 1],
         [1, 1],
         [1, 0, 1, 1],
         [1, 0, 1],
         [1, 0, 0, 1],
         [1],
         [0, 1, 1, 1],
         [0, 1, 1],
         [0, 1, 0, 1],
         [0, 1],
         [0, 0, 1, 1],
         [0, 0, 1],
         [0, 0, 0, 1],
         []]

    and counted::

        sage: L.cardinality()
        16

    Still, the others are considered as elements of `L`::

        sage: L = IntegerListsLex(4,min_length=3,max_length=4)
        sage: L.list()
        [..., [2, 2, 0], ...]

        sage: [2, 2, 0] in L       # in L.list()
        True
        sage: [2, 2, 0, 0] in L    # not in L.list() !
        True
        sage: [2, 2, 0, 0, 0] in L
        False

    .. RUBRIC:: Specifying functions as input for the floor or ceiling

    We construct all lists of sum `4` and length `4` such that ``l[i] <= i``::

        sage: list(IntegerListsLex(4, length=4, ceiling=lambda i: i, check=False))
        [[0, 1, 2, 1], [0, 1, 1, 2], [0, 1, 0, 3], [0, 0, 2, 2], [0, 0, 1, 3]]

    .. WARNING::

        When passing a function as ``floor`` or ``ceiling``, it may
        become undecidable to detect improper inverse lexicographic
        enumeration. For example, the following example has a finite
        enumeration::

            sage: L = IntegerListsLex(3, floor=lambda i: 1 if i>=2 else 0, check=False)
            sage: L.list()
            [[3],
             [2, 1],
             [2, 0, 1],
             [1, 2],
             [1, 1, 1],
             [1, 0, 2],
             [1, 0, 1, 1],
             [0, 3],
             [0, 2, 1],
             [0, 1, 2],
             [0, 1, 1, 1],
             [0, 0, 3],
             [0, 0, 2, 1],
             [0, 0, 1, 2],
             [0, 0, 1, 1, 1]]

        but one cannot decide whether the following has an improper
        inverse lexicographic enumeration without computing the floor
        all the way to ``Infinity``::

            sage: L = IntegerListsLex(3, floor=lambda i: 0, check=False)
            sage: it = iter(L)
            sage: [next(it) for i in range(6)]
            [[3], [2, 1], [2, 0, 1], [2, 0, 0, 1], [2, 0, 0, 0, 1], [2, 0, 0, 0, 0, 1]]

        Hence a warning is raised when a function is specified as
        input, unless the waiver is signed by setting ``check=False``::

            sage: L = IntegerListsLex(3, floor=lambda i: 1 if i>=2 else 0)
            doctest:...
            A function has been given as input of the floor=[...] or ceiling=[...]
            arguments of IntegerListsLex. Please see the documentation for the caveats.
            If you know what you are doing, you can set check=False to skip this warning.

        Similarly, the algorithm may need to search forever for a
        solution when the ceiling is ultimately zero::

            sage: L = IntegerListsLex(2,ceiling=lambda i:0, check=False)
            sage: L.first()           # not tested: will hang forever
            sage: L = IntegerListsLex(2,ceiling=lambda i:0 if i<20 else 1, check=False)
            sage: it = iter(L)
            sage: next(it)
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1]
            sage: next(it)
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1]
            sage: next(it)
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1]


    .. RUBRIC:: Tip: using disjoint union enumerated sets for additional flexibility

    Sometimes, specifying a range for the sum or the length may be too
    restrictive. One would want instead to specify a list, or
    iterable `L`, of acceptable values. This is easy to achieve using
    a :class:`disjoint union of enumerated sets <DisjointUnionEnumeratedSets>`.
    Here we want to accept the values `n=0,2,3`::

        sage: C = DisjointUnionEnumeratedSets(Family([0,2,3],
        ....:         lambda n: IntegerListsLex(n, length=2)))
        sage: C
        Disjoint union of Finite family
        {0: Integer lists of sum 0 satisfying certain constraints,
         2: Integer lists of sum 2 satisfying certain constraints,
         3: Integer lists of sum 3 satisfying certain constraints}
        sage: C.list()
        [[0, 0],
         [2, 0], [1, 1], [0, 2],
         [3, 0], [2, 1], [1, 2], [0, 3]]

    The price to pay is that the enumeration order is now *graded
    lexicographic* instead of lexicographic: first choose the value
    according to the order specified by `L`, and use lexicographic
    order within each value. Here is we reverse `L`::

        sage: DisjointUnionEnumeratedSets(Family([3,2,0],
        ....:     lambda n: IntegerListsLex(n, length=2))).list()
        [[3, 0], [2, 1], [1, 2], [0, 3],
         [2, 0], [1, 1], [0, 2],
         [0, 0]]

    Note that if a given value appears several times, the
    corresponding elements will be enumerated several times, which
    may, or not, be what one wants::

        sage: DisjointUnionEnumeratedSets(Family([2,2],
        ....:     lambda n: IntegerListsLex(n, length=2))).list()
        [[2, 0], [1, 1], [0, 2], [2, 0], [1, 1], [0, 2]]

    Here is a variant where we specify acceptable values for the
    length::

        sage: DisjointUnionEnumeratedSets(Family([0,1,3],
        ....:     lambda l: IntegerListsLex(2, length=l))).list()
        [[2],
         [2, 0, 0], [1, 1, 0], [1, 0, 1], [0, 2, 0], [0, 1, 1], [0, 0, 2]]


    This technique can also be useful to obtain a proper enumeration
    on infinite sets by using a graded lexicographic enumeration::

        sage: C = DisjointUnionEnumeratedSets(Family(NN,
        ....:         lambda n: IntegerListsLex(n, length=2)))
        sage: C
        Disjoint union of Lazy family (<lambda>(i))_{i in Non negative integer semiring}
        sage: it = iter(C)
        sage: [next(it) for i in range(10)]
        [[0, 0],
         [1, 0], [0, 1],
         [2, 0], [1, 1], [0, 2],
         [3, 0], [2, 1], [1, 2], [0, 3]]


    .. RUBRIC:: Specifying how to construct elements

    This is the list of all monomials of degree `4` which divide the
    monomial `x^3y^1z^2` (a monomial being identified with its
    exponent vector)::

        sage: R.<x,y,z> = QQ[]
        sage: m = [3,1,2]
        sage: def term(exponents):
        ....:     return x^exponents[0] * y^exponents[1] * z^exponents[2]
        sage: list( IntegerListsLex(4, length=len(m), ceiling=m, element_constructor=term) )
        [x^3*y, x^3*z, x^2*y*z, x^2*z^2, x*y*z^2]

    Note the use of the ``element_constructor`` option to specify how
    to construct elements from a plain list.

    A variant is to specify a class for the elements. With the default
    element constructor, this class should take as input the parent
    ``self`` and a list.

    .. WARNING::

        The protocol for specifying the element class and constructor
        is subject to changes.

    ALGORITHM:

    The iteration algorithm uses a depth first search through the
    prefix tree of the list of integers (see also
    :ref:`section-generic-integerlistlex`). While doing so, it does
    some lookahead heuristics to attempt to cut dead branches.

    In most practical use cases, most dead branches are cut. Then,
    roughly speaking, the time needed to iterate through all the
    elements of `S` is proportional to the number of elements, where
    the proportion factor is controlled by the length `l` of the
    longest element of `S`. In addition, the memory usage is also
    controlled by `l`, which is to say negligible in practice.

    Still, there remains much room for efficiency improvements; see
    :issue:`18055`, :issue:`18056`.

    .. NOTE::

        The generation algorithm could in principle be extended to
        deal with non-constant slope constraints and with negative
        parts.

    TESTS:

    This example from the combinatorics tutorial used to fail before
    :issue:`17979` because the floor conditions did not satisfy the
    slope conditions::

        sage: I = IntegerListsLex(16, min_length=2, max_slope=-1, floor=[5,3,3])
        sage: I.list()
        [[13, 3], [12, 4], [11, 5], [10, 6], [9, 7], [9, 4, 3], [8, 5, 3], [8, 4, 3, 1],
         [7, 6, 3], [7, 5, 4], [7, 5, 3, 1], [7, 4, 3, 2], [6, 5, 4, 1], [6, 5, 3, 2],
         [6, 4, 3, 2, 1]]

    ::

        sage: Partitions(2, max_slope=-1, length=2).list()                              # needs sage.combinat
        []
        sage: list(IntegerListsLex(0, floor=ConstantFunction(1), min_slope=0))
        [[]]
        sage: list(IntegerListsLex(0, floor=ConstantFunction(1), min_slope=0, max_slope=0))
        [[]]
        sage: list(IntegerListsLex(0, max_length=0, floor=ConstantFunction(1), min_slope=0, max_slope=0))
        [[]]
        sage: list(IntegerListsLex(0, max_length=0, floor=ConstantFunction(0), min_slope=0, max_slope=0))
        [[]]
        sage: list(IntegerListsLex(0, min_part=1, min_slope=0))
        [[]]
        sage: list(IntegerListsLex(1, min_part=1, min_slope=0))
        [[1]]
        sage: list(IntegerListsLex(0, min_length=1, min_part=1, min_slope=0))
        []
        sage: list(IntegerListsLex(0, min_length=1, min_slope=0))
        [[0]]
        sage: list(IntegerListsLex(3, max_length=2))
        [[3], [2, 1], [1, 2], [0, 3]]
        sage: partitions = {"min_part": 1, "max_slope": 0}
        sage: partitions_min_2 = {"floor": ConstantFunction(2), "max_slope": 0}
        sage: compositions = {"min_part": 1}
        sage: integer_vectors = lambda l: {"length": l}
        sage: lower_monomials = lambda c: {"length": c, "floor": lambda i: c[i]}
        sage: upper_monomials = lambda c: {"length": c, "ceiling": lambda i: c[i]}
        sage: constraints = { "min_part":1, "min_slope": -1, "max_slope": 0}
        sage: list(IntegerListsLex(6, **partitions))
        [[6],
         [5, 1],
         [4, 2],
         [4, 1, 1],
         [3, 3],
         [3, 2, 1],
         [3, 1, 1, 1],
         [2, 2, 2],
         [2, 2, 1, 1],
         [2, 1, 1, 1, 1],
         [1, 1, 1, 1, 1, 1]]
        sage: list(IntegerListsLex(6, **constraints))
        [[6],
         [3, 3],
         [3, 2, 1],
         [2, 2, 2],
         [2, 2, 1, 1],
         [2, 1, 1, 1, 1],
         [1, 1, 1, 1, 1, 1]]
        sage: list(IntegerListsLex(1, **partitions_min_2))
        []
        sage: list(IntegerListsLex(2, **partitions_min_2))
        [[2]]
        sage: list(IntegerListsLex(3, **partitions_min_2))
        [[3]]
        sage: list(IntegerListsLex(4, **partitions_min_2))
        [[4], [2, 2]]
        sage: list(IntegerListsLex(5, **partitions_min_2))
        [[5], [3, 2]]
        sage: list(IntegerListsLex(6, **partitions_min_2))
        [[6], [4, 2], [3, 3], [2, 2, 2]]
        sage: list(IntegerListsLex(7, **partitions_min_2))
        [[7], [5, 2], [4, 3], [3, 2, 2]]
        sage: list(IntegerListsLex(9, **partitions_min_2))
        [[9], [7, 2], [6, 3], [5, 4], [5, 2, 2], [4, 3, 2], [3, 3, 3], [3, 2, 2, 2]]
        sage: list(IntegerListsLex(10, **partitions_min_2))
        [[10],
         [8, 2],
         [7, 3],
         [6, 4],
         [6, 2, 2],
         [5, 5],
         [5, 3, 2],
         [4, 4, 2],
         [4, 3, 3],
         [4, 2, 2, 2],
         [3, 3, 2, 2],
         [2, 2, 2, 2, 2]]
        sage: list(IntegerListsLex(4, **compositions))
        [[4], [3, 1], [2, 2], [2, 1, 1], [1, 3], [1, 2, 1], [1, 1, 2], [1, 1, 1, 1]]
        sage: list(IntegerListsLex(6, min_length=1, floor=[7]))
        []
        sage: L = IntegerListsLex(10**100,length=1)
        sage: L.list()
        [[10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000]]

    Noted on :issue:`17898`::

        sage: list(IntegerListsLex(4, min_part=1, length=3, min_slope=1))
        []
        sage: IntegerListsLex(6, ceiling=[4,2], floor=[3,3]).list()
        []
        sage: IntegerListsLex(6, min_part=1, max_part=3, max_slope=-4).list()
        []

    Noted in :issue:`17548`, which are now fixed::

        sage: IntegerListsLex(10, min_part=2, max_slope=-1).list()
        [[10], [8, 2], [7, 3], [6, 4], [5, 3, 2]]
        sage: IntegerListsLex(5, min_slope=1, floor=[2,1,1], max_part=2).list()
        []
        sage: IntegerListsLex(4, min_slope=0, max_slope=0).list()
        [[4], [2, 2], [1, 1, 1, 1]]
        sage: IntegerListsLex(6, min_slope=-1, max_slope=-1).list()
        [[6], [3, 2, 1]]
        sage: IntegerListsLex(6, min_length=3, max_length=2, min_part=1).list()
        []
        sage: I = IntegerListsLex(3, max_length=2, min_part=1)
        sage: I.list()
        [[3], [2, 1], [1, 2]]
        sage: [1,1,1] in I
        False
        sage: I = IntegerListsLex(10, ceiling=[4], max_length=1, min_part=1)
        sage: I.list()
        []
        sage: [4,6] in I
        False
        sage: I = IntegerListsLex(4, min_slope=1, min_part=1, max_part=2)
        sage: I.list()
        []
        sage: I = IntegerListsLex(7, min_slope=1, min_part=1, max_part=4)
        sage: I.list()
        [[3, 4], [1, 2, 4]]
        sage: I = IntegerListsLex(4, floor=[2,1], ceiling=[2,2], max_length=2, min_slope=0)
        sage: I.list()
        [[2, 2]]
        sage: I = IntegerListsLex(10, min_part=1, max_slope=-1)
        sage: I.list()
        [[10], [9, 1], [8, 2], [7, 3], [7, 2, 1], [6, 4], [6, 3, 1], [5, 4, 1],
         [5, 3, 2], [4, 3, 2, 1]]


    .. RUBRIC:: TESTS from comments on :issue:`17979`

    Comment 191::

        sage: list(IntegerListsLex(1, min_length=2, min_slope=0, max_slope=0))
        []

    Comment 240::

        sage: L = IntegerListsLex(min_length=2, max_part=0)
        sage: L.list()
        [[0, 0]]

    .. RUBRIC:: Tests on the element constructor feature and mutability

    Internally, the iterator works on a single list that is mutated
    along the way. Therefore, you need to make sure that the
    ``element_constructor`` actually **copies** its input. This example
    shows what can go wrong::

        sage: P = IntegerListsLex(n=3, max_slope=0, min_part=1, element_constructor=lambda x: x)
        sage: list(P)
        [[], [], []]

    However, specifying ``list()`` as constructor solves this problem::

        sage: P = IntegerListsLex(n=3, max_slope=0, min_part=1, element_constructor=list)
        sage: list(P)
        [[3], [2, 1], [1, 1, 1]]

    Same, step by step::

        sage: it = iter(P)
        sage: a = next(it); a
        [3]
        sage: b = next(it); b
        [2, 1]
        sage: a
        [3]
        sage: a is b
        False

    Tests from `MuPAD-Combinat <http://mupad-combinat.svn.sourceforge.net/viewvc/mupad-combinat/trunk/MuPAD-Combinat/lib/COMBINAT/TEST/MachineIntegerListsLex.tst>`_::

        sage: IntegerListsLex(7, min_length=2, max_length=6, floor=[0,0,2,0,0,1], ceiling=[3,2,3,2,1,2]).cardinality()
        83
        sage: IntegerListsLex(7, min_length=2, max_length=6, floor=[0,0,2,0,1,1], ceiling=[3,2,3,2,1,2]).cardinality()
        53
        sage: IntegerListsLex(5, min_length=2, max_length=6, floor=[0,0,2,0,0,0], ceiling=[2,2,2,2,2,2]).cardinality()
        30
        sage: IntegerListsLex(5, min_length=2, max_length=6, floor=[0,0,1,1,0,0], ceiling=[2,2,2,2,2,2]).cardinality()
        43

        sage: IntegerListsLex(0, min_length=0, max_length=7, floor=[1,1,0,0,1,0], ceiling=[4,3,2,3,2,2,1]).first()
        []

        sage: IntegerListsLex(0, min_length=1, max_length=7, floor=[0,1,0,0,1,0], ceiling=[4,3,2,3,2,2,1]).first()
        [0]
        sage: IntegerListsLex(0, min_length=1, max_length=7, floor=[1,1,0,0,1,0], ceiling=[4,3,2,3,2,2,1]).cardinality()
        0

        sage: IntegerListsLex(2, min_length=0, max_length=7, floor=[1,1,0,0,0,0], ceiling=[4,3,2,3,2,2,1]).first()  # Was [1,1], due to slightly different specs
        [2]
        sage: IntegerListsLex(1, min_length=1, max_length=7, floor=[1,1,0,0,0,0], ceiling=[4,3,2,3,2,2,1]).first()
        [1]
        sage: IntegerListsLex(1, min_length=2, max_length=7, floor=[1,1,0,0,0,0], ceiling=[4,3,2,3,2,2,1]).cardinality()
        0
        sage: IntegerListsLex(2, min_length=5, max_length=7, floor=[1,1,0,0,0,0], ceiling=[4,3,2,3,2,2,1]).first()
        [1, 1, 0, 0, 0]
        sage: IntegerListsLex(2, min_length=5, max_length=7, floor=[1,1,0,0,0,1], ceiling=[4,3,2,3,2,2,1]).first()
        [1, 1, 0, 0, 0]
        sage: IntegerListsLex(2, min_length=5, max_length=7, floor=[1,1,0,0,1,0], ceiling=[4,3,2,3,2,2,1]).cardinality()
        0

        sage: IntegerListsLex(4, min_length=3, max_length=6, floor=[2, 1, 2, 1, 1, 1], ceiling=[3, 1, 2, 3, 2, 2]).cardinality()
        0
        sage: IntegerListsLex(5, min_length=3, max_length=6, floor=[2, 1, 2, 1, 1, 1], ceiling=[3, 1, 2, 3, 2, 2]).first()
        [2, 1, 2]
        sage: IntegerListsLex(6, min_length=3, max_length=6, floor=[2, 1, 2, 1, 1, 1], ceiling=[3, 1, 2, 3, 2, 2]).first()
        [3, 1, 2]
        sage: IntegerListsLex(12, min_length=3, max_length=6, floor=[2, 1, 2, 1, 1, 1], ceiling=[3, 1, 2, 3, 2, 2]).first()
        [3, 1, 2, 3, 2, 1]
        sage: IntegerListsLex(13, min_length=3, max_length=6, floor=[2, 1, 2, 1, 1, 1], ceiling=[3, 1, 2, 3, 2, 2]).first()
        [3, 1, 2, 3, 2, 2]
        sage: IntegerListsLex(14, min_length=3, max_length=6, floor=[2, 1, 2, 1, 1, 1], ceiling=[3, 1, 2, 3, 2, 2]).cardinality()
        0

    This used to hang (see comment 389 and fix in :meth:`Envelope.__init__`)::

        sage: IntegerListsLex(7, max_part=0, ceiling=lambda i:i, check=False).list()
        []
    """
    backend_class = IntegerListsBackend_invlex

    @staticmethod
    def __classcall_private__(cls, n=None, **kwargs):
        r"""
        Specifying a list or iterable as argument was deprecated in
        :issue:`17979`. Please use ``DisjointUnionEnumeratedSets`` or
        the ``min_sum`` and ``max_sum`` arguments instead.
        """
        return typecall(cls, n=n, **kwargs)


cdef class IntegerListsBackend_invlex(IntegerListsBackend):
    """
    Cython back-end of a set of lists of integers with specified
    constraints enumerated in inverse lexicographic order.
    """
    def __init__(self, *args, check=True, **kwds):
        """
        Initialize ``self``.

        TESTS::

            sage: C = IntegerListsLex(2, length=3)
            sage: C == loads(dumps(C))
            True
            sage: C.cardinality().parent() is ZZ
            True
            sage: TestSuite(C).run()
            sage: IntegerListsLex(min_part=-1)
            Traceback (most recent call last):
            ...
            NotImplementedError: strictly negative min_part
        """
        IntegerListsBackend.__init__(self, *args, **kwds)

        self.check = check

        if self.min_part < 0:
            raise NotImplementedError("strictly negative min_part")

        if self.check and (
                self.floor.limit_start() == Infinity or
                self.ceiling.limit_start() == Infinity):
            from warnings import warn
            warn("""
A function has been given as input of the floor=[...] or ceiling=[...]
arguments of IntegerListsLex. Please see the documentation for the caveats.
If you know what you are doing, you can set check=False to skip this warning.""")

    @cached_method
    def _check_finiteness(self):
        """
        Check that the constraints define a finite set.

        As mentioned in the description of this class, being finite is
        almost equivalent to being inverse lexicographic iterable,
        which is what we really care about.

        This set is finite if and only if:

        #. For each `i` such that there exists a list of length at
           least `i+1` satisfying the constraints, there exists a
           direct or indirect upper bound on the `i`-th part, that
           is ``self.ceiling(i)`` is finite.

        #. There exists a global upper bound on the length.

        Failures for 1. are detected and reported later, during the
        iteration, namely the first time a prefix including the `i`-th
        part is explored.

        This method therefore focuses on 2., namely trying to prove
        the existence of an upper bound on the length. It may fail
        to do so even when the set is actually finite.

        OUTPUT:

        ``True`` if this method finds a proof that there
        exists an upper bound on the length. Otherwise a
        :exc:`ValueError` is raised.

        Note that :func:`cached_method` does not work with methods
        returning ``None``, so ``True`` is returned instead.

        EXAMPLES::

            sage: L = IntegerListsLex(4, max_length=4)
            sage: L._check_finiteness()

        The following example is infinite::

            sage: L = IntegerListsLex(4)
            sage: L._check_finiteness()
            Traceback (most recent call last):
            ...
            ValueError: could not prove that the specified constraints yield a finite set

        Indeed::

            sage: it = iter(IntegerListsLex(4, check=False))
            sage: for _ in range(10): print(next(it))
            [4]
            [3, 1]
            [3, 0, 1]
            [3, 0, 0, 1]
            [3, 0, 0, 0, 1]
            [3, 0, 0, 0, 0, 1]
            [3, 0, 0, 0, 0, 0, 1]
            [3, 0, 0, 0, 0, 0, 0, 1]
            [3, 0, 0, 0, 0, 0, 0, 0, 1]
            [3, 0, 0, 0, 0, 0, 0, 0, 0, 1]

        Unless ``check=False``, :meth:`_check_finiteness` is called as
        soon as an iteration is attempted::

            sage: iter(L)
            Traceback (most recent call last):
            ...
            ValueError: could not prove that the specified constraints yield a finite set

        Some other infinite examples::

            sage: L = IntegerListsLex(ceiling=[0], min_slope=1, max_slope=2)
            sage: L.list()
            Traceback (most recent call last):
            ...
            ValueError: could not prove that the specified constraints yield a finite set

            sage: L = IntegerListsLex(ceiling=[0], min_slope=1, max_slope=1)
            sage: L.list()
            Traceback (most recent call last):
            ...
            ValueError: could not prove that the specified constraints yield a finite set

            sage: IntegerListsLex(ceiling=[0], min_slope=1, max_slope=1).list()
            Traceback (most recent call last):
            ...
            ValueError: could not prove that the specified constraints yield a finite set

        The following example is actually finite, but not detected as such::

            sage: IntegerListsLex(7, floor=[4], max_part=4, min_slope=-1).list()
            Traceback (most recent call last):
            ...
            ValueError: could not prove that the specified constraints yield a finite set

        This is sad because the following equivalent example works just fine::

            sage: IntegerListsLex(7, floor=[4,3], max_part=4, min_slope=-1).list()
            [[4, 3]]

        Detecting this properly would require some deeper lookahead,
        and the difficulty is to decide how far this lookahead should
        search. Until this is fixed, one can disable the checks::

            sage: IntegerListsLex(7, floor=[4], max_part=4, min_slope=-1, check=False).list()
            [[4, 3]]

        If the ceiling or floor is a function, it is much more likely
        that a finite set will not be detected as such::

            sage: IntegerListsLex(ceiling=lambda i: max(3-i,0))._check_finiteness()
            Traceback (most recent call last):
            ...
            ValueError: could not prove that the specified constraints yield a finite set

            sage: IntegerListsLex(7, ceiling=lambda i:0).list()
            Traceback (most recent call last):
            ...
            ValueError: could not prove that the specified constraints yield a finite set

        The next example shows a case that is finite because we remove
        trailing zeroes::

            sage: list(IntegerListsLex(ceiling=[0], max_slope=0))
            [[]]
            sage: L = IntegerListsLex(ceiling=[1], min_slope=1, max_slope=1)
            sage: L.list()
            Traceback (most recent call last):
            ...
            ValueError: could not prove that the specified constraints yield a finite set

        In the next examples, there is either no solution, or the region
        is bounded::

            sage: IntegerListsLex(min_sum=10, max_sum=5).list()
            []
            sage: IntegerListsLex(max_part=1, min_slope=10).list()
            [[1], []]
            sage: IntegerListsLex(max_part=100, min_slope=10).first()
            [100]
            sage: IntegerListsLex(ceiling=[1,Infinity], max_part=2, min_slope=1).list()
            [[1, 2], [1], [0, 2], [0, 1, 2], [0, 1], []]
            sage: IntegerListsLex(min_sum=1, floor=[1,2], max_part=1).list()
            [[1]]

            sage: IntegerListsLex(min_length=2, max_length=1).list()
            []
            sage: IntegerListsLex(min_length=-2, max_length=-1).list()
            []
            sage: IntegerListsLex(min_length=-1, max_length=-2).list()
            []
            sage: IntegerListsLex(min_length=2, max_slope=0, min_slope=1).list()
            []
            sage: IntegerListsLex(min_part=2, max_part=1).list()
            [[]]

            sage: IntegerListsLex(floor=[0,2], ceiling=[3,1]).list()
            [[3], [2], [1], []]
            sage: IntegerListsLex(7, ceiling=[2], floor=[4]).list()
            []
            sage: IntegerListsLex(7, max_part=0).list()
            []
            sage: IntegerListsLex(5, max_part=0, min_slope=0).list()
            []
            sage: IntegerListsLex(max_part=0).list()
            [[]]
            sage: IntegerListsLex(max_sum=1, min_sum=4, min_slope=0).list()
            []
        """
        # Trivial cases
        if self.max_length < Infinity:
            return True
        if self.max_sum < self.min_sum:
            return True
        if self.min_slope > self.max_slope:
            return True
        if self.max_slope < 0:
            return True
        if self.ceiling.limit() < self.floor.limit():
            return True
        if self.ceiling.limit() == 0:
            # This assumes no trailing zeroes
            return True
        if self.min_slope > 0 and self.ceiling.limit() < Infinity:
            return True

        # Compute a lower bound on the sum of floor(i) for i=1 to infinity
        if self.floor.limit() > 0 or self.min_slope > 0:
            floor_sum_lower_bound = Infinity
        elif self.floor.limit_start() < Infinity:
            floor_sum_lower_bound = sum(self.floor(i) for i in range(self.floor.limit_start()))
        else:
            floor_sum_lower_bound = 0
        if floor_sum_lower_bound > 0 and self.min_slope >= 0:
            floor_sum_lower_bound = Infinity

        if self.max_sum < floor_sum_lower_bound:
            return True
        if self.max_sum == floor_sum_lower_bound and self.max_sum < Infinity:
            # This assumes no trailing zeroes
            return True

        # Variant on ceiling.limit() ==0 where we actually discover that the ceiling limit is 0
        if ( self.max_slope == 0 and
             (self.max_sum < Infinity or
              (self.ceiling.limit_start() < Infinity and
               any(self.ceiling(i) == 0 for i in range(self.ceiling.limit_start()+1)))
             ) ):
            return True

        limit_start = max(self.ceiling.limit_start(), self.floor.limit_start())
        if limit_start < Infinity:
            for i in range(limit_start+1):
                if self.ceiling(i) < self.floor(i):
                    return True

        raise ValueError("could not prove that the specified constraints yield a finite set")

    def _iter(self):
        """
        Return an iterator for the elements of ``self``.

        EXAMPLES::

            sage: C = IntegerListsLex(2, length=3)
            sage: list(C)     # indirect doctest
            [[2, 0, 0], [1, 1, 0], [1, 0, 1], [0, 2, 0], [0, 1, 1], [0, 0, 2]]
        """
        if self.check:
            self._check_finiteness()
        return IntegerListsLexIter(self)


# Constants for IntegerListsLexIter._next_state
LOOKAHEAD = 5
PUSH      = 4
ME        = 3
DECREASE  = 2
POP       = 1
STOP      = 0


class IntegerListsLexIter(builtins.object):
    r"""
    Iterator class for IntegerListsLex.

    Let ``T`` be the prefix tree of all lists of nonnegative
    integers that satisfy all constraints except possibly for
    ``min_length`` and ``min_sum``; let the children of a list
    be sorted decreasingly according to their last part.

    The iterator is based on a depth-first search exploration of a
    subtree of this tree, trying to cut branches that do not
    contain a valid list. Each call of ``next`` iterates through
    the nodes of this tree until it finds a valid list to return.

    Here are the attributes describing the current state of the
    iterator,  and their invariants:

    - ``backend`` -- the :class:`IntegerListsBackend` object this is
      iterating on;

    - ``_current_list`` -- the list corresponding to the current
      node of the tree;

    - ``_j`` -- the index of the last element of ``_current_list``:
      ``self._j == len(self._current_list) - 1``;

    - ``_current_sum`` -- the sum of the parts of ``_current_list``;

    - ``_search_ranges`` -- list of same length as
      ``_current_list``: the range for each part

    Furthermore, we assume that there is no obvious contradiction
    in the constraints:

    - ``self.backend.min_length <= self.backend.max_length``;
    - ``self.backend.min_slope <= self.backend.max_slope``
      unless ``self.backend.min_length <= 1``.

    Along this iteration, ``next`` switches between the following
    states:

    - LOOKAHEAD: determine whether the current list could be a
      prefix of a valid list;
    - PUSH: go deeper into the prefix tree by appending the
      largest possible part to the current list;
    - ME: check whether the current list is valid and if yes return it
    - DECREASE: decrease the last part;
    - POP: pop the last part of the current list;
    - STOP: the iteration is finished.

    The attribute ``_next_state`` contains the next state ``next``
    should enter in.
    """
    def __init__(self, backend):
        """
        TESTS::

            sage: from sage.combinat.integer_lists.invlex import IntegerListsLexIter
            sage: C = IntegerListsLex(2, length=3)
            sage: I = IntegerListsLexIter(C)
            sage: I._search_ranges
            []
            sage: I._current_list
            []
            sage: I._j
            -1
            sage: I._current_sum
            0
        """
        self.backend = backend

        self._search_ranges = []
        self._current_list = []
        self._j = -1     # index of last element of _current_list
        self._current_sum = 0     # sum of parts in _current_list

        # Make sure that some invariants are respected in the iterator
        if (backend.min_length <= backend.max_length and
           (backend.min_slope <= backend.max_slope or backend.min_length <= 1)):
            self._next_state = PUSH
        else:
            self._next_state = STOP

    def __iter__(self):
        """
        Return ``self`` as per the iterator protocol.

        EXAMPLES::

            sage: from sage.combinat.integer_lists.invlex import IntegerListsLexIter
            sage: C = IntegerListsLex(2, length=3)
            sage: it = IntegerListsLexIter(C)
            sage: it.__iter__() is it
            True
        """
        return self

    def _push_search(self):
        """
        Push search forward, resetting attributes.

        The push may fail if it is discovered that
        ``self._current_list`` cannot be extended in a valid way.

        OUTPUT: boolean; whether the push succeeded

        EXAMPLES::

            sage: C = IntegerListsLex(2, length=3)
            sage: I = C._iter()
            sage: I._j
            -1
            sage: I._search_ranges
            []
            sage: I._current_list
            []
            sage: I._current_sum
            0
            sage: I._push_search()
            True
            sage: I._j
            0
            sage: I._search_ranges
            [(0, 2)]
            sage: I._current_list
            [2]
            sage: I._current_sum
            2
            sage: I._push_search()
            True
            sage: I._j
            1
            sage: I._search_ranges
            [(0, 2), (0, 0)]
            sage: I._current_list
            [2, 0]
            sage: I._current_sum
            2
        """
        max_sum = self.backend.max_sum
        min_length = self.backend.min_length
        max_length = self.backend.max_length
        if self._j + 1 >= max_length:
            return False
        if self._j + 1 >= min_length and self._current_sum == max_sum:
            # Cannot add trailing zeroes
            return False

        if self._j >= 0:
            prev = self._current_list[self._j]
        else:
            prev = None
        interval = self._m_interval(self._j+1, self.backend.max_sum - self._current_sum, prev)
        if interval[0] > interval[1]:
            return False

        self._j += 1
        m = interval[1]
        self._search_ranges.append(interval)
        self._current_list.append(m)
        self._current_sum += m
        return True

    def _pop_search(self):
        """
        Go back in search tree. Resetting attributes.

        EXAMPLES::

            sage: C = IntegerListsLex(2, length=3)
            sage: I = C._iter()
            sage: I._push_search()
            True
            sage: I._j
            0
            sage: I._search_ranges
            [(0, 2)]
            sage: I._current_sum
            2
            sage: I._current_list
            [2]
            sage: I._pop_search()
            sage: I._j
            -1
            sage: I._search_ranges
            []
            sage: I._current_sum
            0
            sage: I._current_list
            []
        """
        if self._j >= 0:  # TODO: get rid of this condition
            self._j -= 1
            self._search_ranges.pop()
            self._current_sum -= self._current_list[-1]
            self._current_list.pop()

    def __next__(self):
        r"""
        Return the next element in the iteration.

        EXAMPLES::

            sage: from sage.combinat.integer_lists.invlex import IntegerListsLexIter
            sage: C = IntegerListsLex(2, length=3)
            sage: I = IntegerListsLexIter(C)
            sage: next(I)
            [2, 0, 0]
            sage: next(I)
            [1, 1, 0]
        """
        p = self.backend
        min_sum = p.min_sum
        max_length = p.max_length
        search_ranges = self._search_ranges

        while True:
            assert self._j == len(self._current_list) - 1
            assert self._j == len(self._search_ranges) - 1

            # LOOK AHEAD
            if self._next_state == LOOKAHEAD:
                if self._lookahead():
                    self._next_state = PUSH
                else:
                    # We should reuse information about the
                    # reasons for this failure, to avoid when
                    # possible retrying with smaller values.
                    # We just do a special case for now:
                    if self._j + 1 == max_length and self._current_sum < min_sum:
                        self._next_state = POP
                    else:
                        self._next_state = DECREASE

            # PUSH
            if self._next_state == PUSH:
                if self._push_search():
                    self._next_state = LOOKAHEAD
                    continue
                self._next_state = ME

            # ME
            if self._next_state == ME:
                if self._j == -1:
                    self._next_state = STOP
                else:
                    self._next_state = DECREASE
                if self._internal_list_valid():
                    return self._current_list

            # DECREASE
            if self._next_state == DECREASE:
                self._current_list[-1] -= 1
                self._current_sum -= 1
                if self._current_list[-1] >= search_ranges[self._j][0]:
                    self._next_state = LOOKAHEAD
                    continue
                self._next_state = POP

            # POP
            if self._next_state == POP:
                self._pop_search()
                self._next_state = ME
                continue

            # STOP
            if self._next_state == STOP:
                raise StopIteration()

            assert False

    def _internal_list_valid(self):
        """
        Return whether the current list in the iteration variable
        ``self._current_list`` is a valid list.

        This method checks whether the sum of the parts in
        ``self._current_list`` is in the right range, whether its
        length is in the required range, and whether there are trailing
        zeroes.  It does not check all of the necessary conditions to
        verify that an arbitrary list satisfies the constraints from
        the corresponding ``IntegerListsLex`` object, and should
        not be used except internally in the iterator class.

        EXAMPLES::

            sage: from sage.combinat.integer_lists.invlex import IntegerListsLexIter
            sage: C = IntegerListsLex(2, length=3)
            sage: I = IntegerListsLexIter(C)
            sage: I._current_list
            []
            sage: I._internal_list_valid()
            False
            sage: next(I)
            [2, 0, 0]
            sage: I._current_list
            [2, 0, 0]
            sage: I._internal_list_valid()
            True
        """
        p = self.backend
        mu = self._current_list
        nu = self._current_sum
        l = self._j + 1
        return (nu >= p.min_sum and nu <= p.max_sum # Good sum
                 and l >= p.min_length and l <= p.max_length # Good length
                 and (l <= max(p.min_length,0) or mu[-1] != 0)) # No trailing zeros

    def _m_interval(self, i, max_sum, prev=None):
        r"""
        Return coarse lower and upper bounds for the part ``m``
        at position ``i``.

        INPUT:

        - ``i`` -- nonnegative integer (position)

        - ``max_sum`` -- nonnegative integer or ``+oo``

        - ``prev`` -- nonnegative integer or ``None``

        Return coarse lower and upper bounds for the value ``m``
        of the part at position ``i`` so that there could exists
        some list suffix `v_i,\ldots,v_k` of sum bounded by
        ``max_sum`` and satisfying the floor and upper bound
        constraints. If ``prev`` is specified, then the slope
        conditions between ``v[i-1]=prev`` and ``v[i]=m`` should
        also be satisfied.

        Additionally, this raises an error if it can be detected
        that some part is neither directly nor indirectly bounded
        above, which implies that the constraints possibly do not
        allow for an inverse lexicographic iterator.

        OUTPUT:

        A tuple of two integers ``(lower_bound, upper_bound)``.

        EXAMPLES::

            sage: from sage.combinat.integer_lists.invlex import IntegerListsLexIter
            sage: C = IntegerListsLex(2, length=3)
            sage: I = IntegerListsLexIter(C)
            sage: I._m_interval(1,2)
            (0, 2)

        The second part is not bounded above, hence we cannot
        iterate lexicographically through all the elements::

            sage: IntegerListsLex(ceiling=[2,infinity,3], max_length=3).first()
            Traceback (most recent call last):
            ...
            ValueError: infinite upper bound for values of m

        Same here::

            sage: IntegerListsLex(ceiling=[2,infinity,2], max_length=3, min_slope=-1).cardinality()
            Traceback (most recent call last):
            ...
            ValueError: infinite upper bound for values of m

        In the following examples however, all parts are
        indirectly bounded above::

            sage: IntegerListsLex(ceiling=[2,infinity,2], length=3,     min_slope=-1).cardinality()
            24
            sage: IntegerListsLex(ceiling=[2,infinity,2], max_length=3, max_slope=1).cardinality()
            24

            sage: IntegerListsLex(max_part=2, max_length=3).cardinality()
            27
            sage: IntegerListsLex(3, max_length=3).cardinality()      # parts bounded by n
            10
            sage: IntegerListsLex(max_length=0, min_length=1).list()  # no part!
            []
            sage: IntegerListsLex(length=0).list()                    # no part!
            [[]]
        """
        p = self.backend

        lower_bound = max(0, p.floor(i))
        upper_bound = min(max_sum, p.ceiling(i))
        if prev is not None:
            lower_bound = max(lower_bound, prev + p.min_slope)
            upper_bound = min(upper_bound, prev + p.max_slope)

        ## check for infinite upper bound, in case max_sum is infinite
        if p.check and upper_bound == Infinity:
            # This assumes that there exists a valid list (which
            # is not yet always guaranteed). Then we just
            # discovered that part 'i' of this list can be made as
            # large as desired, which implies that `self.backend`
            # cannot be iterated in inverse lexicographic order
            raise ValueError("infinite upper bound for values of m")

        return (lower_bound, upper_bound)

    def _lookahead(self):
        r"""
        Return whether the current list can possibly be a prefix
        of a valid list.

        OUTPUT:

        ``False`` if it is guaranteed that the current list
        cannot be a prefix of a valid list and ``True`` otherwise.

        EXAMPLES::

            sage: it = IntegerListsLex(length=3, min_sum=2, max_sum=2)._iter()
            sage: it._current_list = [0,1]    # don't do this at home, kids
            sage: it._current_sum = 1
            sage: it._j = 1
            sage: it._lookahead()
            True

            sage: it = IntegerListsLex(length=3, min_sum=3, max_sum=2)._iter()
            sage: it._current_list = [0,1]
            sage: it._current_sum = 1
            sage: it._j = 1
            sage: it._lookahead()
            False

            sage: it = IntegerListsLex(min_length=2, max_part=0)._iter()
            sage: it._current_list = [0]
            sage: it._current_sum = 0
            sage: it._j = 0
            sage: it._lookahead()
            True
            sage: it._current_list = [0, 0]
            sage: it._j = 1
            sage: it._lookahead()
            True
            sage: it._current_list = [0, 0, 0]
            sage: it._j = 2
            sage: it._lookahead()
            False

            sage: n = 10**100
            sage: it = IntegerListsLex(n, length=1)._iter()
            sage: it._current_list = [n-1]
            sage: it._current_sum = n-1
            sage: it._j = 0
            sage: it._lookahead()
            False

            sage: it = IntegerListsLex(n=3, min_part=2, min_sum=3, max_sum=3)._iter()
            sage: it._current_list = [2]
            sage: it._current_sum = 2
            sage: it._j = 0
            sage: it._lookahead()
            False

        ALGORITHM:

        Let ``j = self._j`` be the position of the last part `m` of
        ``self._current_list``. The current algorithm computes,
        for `k = j, j+1, \ldots`, a lower bound `l_k` and an upper
        bound `u_k` for `v_0+\dots+v_k`, and stops if none of the
        intervals `[l_k, u_k]` intersect ``[min_sum, max_sum]``.

        The lower bound `l_k` is given by the area below
        `v_0,\dots,v_{j-1}` prolongated by the lower envelope
        between `j` and `k` and starting at `m`. The upper bound
        `u_k` is given similarly using the upper envelope.

        The complexity of this algorithm is bounded above by
        ``O(max_length)``. When ``max_length=oo``, the algorithm
        is guaranteed to terminate, unless ``floor`` is a function
        which is eventually constant with value `0`, or which
        reaches the value `0` while ``max_slope=0``.

        Indeed, the lower bound `l_k` is increasing with `k`; in
        fact it is strictly increasing, unless the local lower bound
        at `k` is `0`. Furthermore as soon as ``l_k >= min_sum``,
        we can conclude; we can also conclude if we know that the
        floor is eventually constant with value `0`, or there is a
        local lower bound at `k` is `0` and ``max_slope=0``.

        .. RUBRIC:: Room for improvement

        Improved prediction: the lower bound `l_k` does not take
        the slope conditions into account, except for those imposed
        by the value `m` at `j`. Similarly for `u_k`.

        Improved speed: given that `l_k` is increasing with `k`,
        possibly some dichotomy could be used to search for `k`,
        with appropriate caching / fast calculation of the partial
        sums. Also, some of the information gained at depth `j`
        could be reused at depth `j+1`.

        TESTS::

            sage: it = IntegerListsLex(1, min_length=2, min_slope=0, max_slope=0, min_sum=1, max_sum=1)._iter()
            sage: it._current_list = [0]
            sage: it._current_sum = 0
            sage: it._j = 0
            sage: it._lookahead()
            False
        """
        # Check code for various termination conditions.  Possible cases:
        # 0. interval [lower, upper] intersects interval [min_sum, max_sum] -- terminate True
        # 1. lower sum surpasses max_sum -- terminate False
        # 2. iteration surpasses max_length -- terminate False
        # 3. upper envelope is smaller than lower envelope -- terminate False
        # 4. max_slope <= 0 -- terminate False after upper passes 0
        # 5. ceiling_limit == 0 -- terminate False after reaching larger limit point
        # 6. (uncomputable) ceiling function == 0 for all but finitely many input values -- terminate False after reaching (unknown) limit point -- currently hangs

        m = self._current_list[-1]
        j = self._j
        min_sum = self.backend.min_sum - (self._current_sum-m)
        max_sum = self.backend.max_sum - (self._current_sum-m)

        if min_sum > max_sum:
            return False

        p = self.backend

        # Beware that without slope conditions, the functions below
        # currently forget about the value m at k!
        lower_envelope = self.backend.floor.adapt(m,j)
        upper_envelope = self.backend.ceiling.adapt(m,j)
        lower = 0    # The lower bound `l_k`
        upper = 0    # The upper bound `u_k`

        assert j >= 0
        # get to smallest valid number of parts
        for k in range(j, p.min_length-1):
            # We are looking at lists `v_j,...,v_k`
            lo = m if k == j else lower_envelope(k)
            up = m if k == j else upper_envelope(k)
            if lo > up:
                return False
            lower += lo
            upper += up

        if j < p.min_length and min_sum <= upper and lower <= max_sum:
            # There could exist a valid list `v_j,\dots,v_{min_length-1}`
            return True

        k = max(p.min_length-1,j)
        # Check if any of the intervals intersect the target interval
        while k < p.max_length:
            lo = m if k == j else lower_envelope(k)
            up = m if k == j else upper_envelope(k)
            if lo > up:
                # There exists no valid list of length >= k
                return False
            lower += lo
            upper += up
            assert lower <= upper

            if lower > max_sum:
                # There cannot exist a valid list `v_j,\dots,v_l` with l>=k
                return False

            if ( (p.max_slope <= 0 and up <= 0)
                 or (p.ceiling.limit() == 0 and k > p.ceiling.limit_start()) ):
                # This implies v_l=0 for l>=k: that is we would be generating
                # a list with trailing zeroes
                return False

            if min_sum <= upper and lower <= max_sum:
                # There could exist a valid list `v_j,\dots,v_k`
                return True

            k += 1

        return False
