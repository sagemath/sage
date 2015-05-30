r"""
Enumerated set of lists of integers with constraints, in inverse lexicographic order

- :class:`IntegerListsLex`: the enumerated set of lists of nonnegative
  integers with specified constraints, in inverse lexicographic order.

- :class:`Envelope`: a utility class for upper (lower) envelope of a
  function under constraints.

HISTORY:

This generic tool was originally written by Hivert and Thiery in
MuPAD-Combinat in 2000 and ported over to Sage by Mike Hansen in
2007. It was then completely rewritten in 2015 by Gillespie,
Schilling, and Thiery, with the help of many, to deal with
limitations and lack of robustness w.r.t. input.
"""
#*****************************************************************************
#       Copyright (C) 2015 Bryan Gillespie <Brg008@gmail.com>
#                          Nicolas M. Thiery <nthiery at users.sf.net>
#                          Anne Schilling <anne@math.ucdavis.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from inspect import ismethod
from sage.misc.classcall_metaclass import ClasscallMetaclass, typecall
from sage.misc.constant_function import ConstantFunction
from sage.misc.cachefunc import cached_method
from sage.categories.enumerated_sets import EnumeratedSets
from sage.structure.list_clone import ClonableArray
from sage.structure.parent import Parent
from sage.sets.family import Family
from sage.rings.integer_ring import ZZ

Infinity = float('+inf')

class IntegerListsLex(Parent):
    r"""
    Lists of nonnegative integers with constraints, in inverse lexicographic order.

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

    - ``min_sum`` -- a nonnegative integer (default: 0):
      a lower bound on ``sum(l)``.

    - ``max_sum`` -- a nonnegative integer or `\infty` (default: `\infty`):
      an upper bound on ``sum(l)``.

    - ``n`` -- a nonnegative integer (optional): if specified, this
      overrides ``min_sum`` and ``max_sum``.

    - ``min_length`` -- a nonnegative integer (default: `0`): a lower
      bound on ``len(l)``.

    - ``max_length`` -- a nonnegative integer or `\infty` (default:
      `\infty`): an upper bound on ``len(l)``.

    - ``length`` -- an integer (optional); overrides ``min_length``
      and ``max_length`` if specified;

    - ``min_part`` -- a nonnegative integer: a lower bounds on all
       parts: ``min_part <= l[i]`` for ``0 <= i < len(l)``.

    - ``floor`` -- a list of nonnegative integers or a function: lower
      bounds on the individual parts `l[i]`.

      If ``floor`` is a list of integers, then ``floor<=l[i]`` for ``0
      <= i < min(len(l), len(floor)``. Similarly, if ``floor`` is a
      function, then ``floor(i) <= l[i]`` for ``0 <= i < len(l)``.

    - ``max_part`` -- a nonnegative integer or `\infty`: an upper
      bound on all parts: ``l[i] <= max_part`` for ``0 <= i < len(l)``.

    - ``ceiling`` -- upper bounds on the individual parts ``l[i]``;
      this takes the same type of input as ``floor``, except that
      `\infty` is allowed in addition to integers, and the default
      value is `\infty`.

    - ``min_slope`` -- an integer or `-\infty` (default: `-\infty`):
      an lower bound on the slope between consecutive parts:
      ``min_slope <= l[i+1]-l[i]`` for ``0 <= i < len(l)-1``

    - ``max_slope`` -- an integer or `+\infty` (defaults: `+\infty`)
      an upper bound on the slope between consecutive parts:
      ``l[i+1]-l[i] <= max_slope`` for ``0 <= i < len(l)-1``

    - ``category`` -- a category (default: :class:`FiniteEnumeratedSets`)

    - ``check`` -- boolean (default: ``True``): whether to display the
      warnings raised when functions are given as input to ``floor``
      or ``ceiling`` and the errors raised when there is no proper
      enumeration.

    - ``name`` -- a string or ``None`` (default: ``None``) if set,
      this will be passed down to :meth:`Parent.rename` to specify the
      name of ``self``. It is recommented to use directly the rename
      method as this feature may become deprecated.

    - ``element_constructor`` -- a function (or callable) that creates
      elements of ``self`` from a list. See also :class:`Parent`.

    - ``element_class`` -- a class for the elements of ``self``
      (default: `ClonableArray`). This merely sets the attribute
      ``self.Element``. See the examples for details.

    - ``global_options`` -- a :class:`~sage.structure.global_options.GlobalOptions`
      object that will be assigned to the attribute
      ``_global_options``; for internal use only (subclasses, ...).


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

        sage: [L._floor(i) for i in range(5)]
        [2, 1, 2, 1, 1]

    Similarly, the ``ceiling`` and ``max_part`` constraints can be
    combined::

        sage: L = IntegerListsLex(4, ceiling=[2,3,1], max_part=2, length=3)
        sage: L.list()
        [[2, 2, 0], [2, 1, 1], [1, 2, 1]]
        sage: [L._ceiling(i) for i in range(5)]
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

    or Dyck words (see also :class:`DyckWords`), through the bijection
    with paths from `(0,0)` to `(n,n)` with left and up steps that remain
    below the diagonal::

        sage: def dyck_words(n):
        ....:     return IntegerListsLex(length=n, ceiling=range(n+1), min_slope=0)
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
        ValueError: Could not prove that the specified constraints yield a finite set

    Here is a variant which could be enumerated in lexicographic order
    but not in inverse lexicographic order::

        sage: L = IntegerListsLex(length=2, ceiling=[Infinity, 0], floor=[0,1])
        sage: for l in L: print l
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
        ValueError: Could not prove that the specified constraints yield a finite set

    If one wants to proceed anyway, one can sign a waiver by setting
    ``check=False`` (again, be warned that some valid lists may never appear)::

        sage: L = IntegerListsLex(3, check=False)
        sage: it = iter(L)
        sage: [it.next() for i in range(6)]
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
            sage: [it.next() for i in range(6)]
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
            sage: it.next()
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1]
            sage: it.next()
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1]
            sage: it.next()
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
        sage: [it.next() for i in range(10)]
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
    ``self`` and a list. Here we want the elements to be constructed
    in the class :class:`Partition`::

        sage: IntegerListsLex(3, max_slope=0, element_class=Partition, global_options=Partitions.global_options).list()
        [[3], [2, 1], [1, 1, 1]]

    Note that the :class:`Partition` further assumes the existence of
    an attribute ``_global_options`` in the parent, hence the use of the
    ``global_options`` parameter.

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
    :trac:`18055`, :trac:`18056`.

    .. NOTE::

        The generation algorithm could in principle be extended to
        deal with non-constant slope constraints and with negative
        parts.

    TESTS:

    This example from the combinatorics tutorial used to fail before
    :trac:`17979` because the floor conditions did not satisfy the
    slope conditions::

        sage: I = IntegerListsLex(16, min_length=2, max_slope=-1, floor=[5,3,3])
        sage: I.list()
        [[13, 3], [12, 4], [11, 5], [10, 6], [9, 7], [9, 4, 3], [8, 5, 3], [8, 4, 3, 1],
         [7, 6, 3], [7, 5, 4], [7, 5, 3, 1], [7, 4, 3, 2], [6, 5, 4, 1], [6, 5, 3, 2],
         [6, 4, 3, 2, 1]]

    ::

        sage: Partitions(2, max_slope=-1, length=2).list()
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

    Noted on :trac:`17898`::

        sage: list(IntegerListsLex(4, min_part=1, length=3, min_slope=1))
        []
        sage: IntegerListsLex(6, ceiling=[4,2], floor=[3,3]).list()
        []
        sage: IntegerListsLex(6, min_part=1, max_part=3, max_slope=-4).list()
        []

    Noted in :trac:`17548`, which are now fixed::

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
        sage: I=IntegerListsLex(10, ceiling=[4], max_length=1, min_part=1)
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


    .. RUBRIC:: TESTS from comments on :trac:`17979`

    Comment 191::

        sage: list(IntegerListsLex(1, min_length=2, min_slope=0, max_slope=0))
        []

    Comment 240::

        sage: L = IntegerListsLex(min_length=2, max_part=0)
        sage: L.list()
        [[0, 0]]

    .. RUBRIC:: Tests on the element constructor feature and mutability

    Internally, the iterator works on a single list that is mutated
    along the way. The following test makes sure that we actually make a copy of
    this list before passing it to ``element_constructor`` in order to
    avoid reference effects::

        sage: from sage.misc.c3_controlled import identity
        sage: P = IntegerListsLex(n=3, max_slope=0, min_part=1, element_constructor=identity)
        sage: list(P)
        [[3], [2, 1], [1, 1, 1]]

    Same, step by step::

        sage: it = iter(P)
        sage: a = it.next(); a
        [3]
        sage: b = it.next(); b
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
    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, n=None, **kwargs):
        r"""
        Return a disjoint union if ``n`` is a list or iterable.

        TESTS:

        Specifying a list or iterable as argument is deprecated::

            sage: IntegerListsLex([2,2], length=2).list()
            doctest:...: DeprecationWarning: Calling IntegerListsLex with n an iterable is deprecated. Please use DisjointUnionEnumeratedSets or the min_sum and max_sum arguments instead
            See http://trac.sagemath.org/17979 for details.
            [[2, 0], [1, 1], [0, 2], [2, 0], [1, 1], [0, 2]]
            sage: IntegerListsLex(NN, max_length=3)
            Disjoint union of Lazy family (<lambda>(i))_{i in Non negative integer semiring}
        """
        import collections
        if isinstance(n, collections.Iterable):
            from sage.misc.superseded import deprecation
            deprecation(17979, 'Calling IntegerListsLex with n an iterable is deprecated. Please use DisjointUnionEnumeratedSets or the min_sum and max_sum arguments instead')
            from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
            return DisjointUnionEnumeratedSets(Family(n, lambda i: IntegerListsLex(i, **kwargs)))
        else:
            return typecall(cls, n=n, **kwargs)

    def __init__(self,
                 n=None,
                 length=None, min_length=0, max_length=Infinity,
                 floor=None, ceiling=None,
                 min_part=0, max_part=Infinity,
                 min_slope=-Infinity, max_slope=Infinity,
                 min_sum=0, max_sum=Infinity,
                 name=None,
                 category=None,
                 element_constructor=None, element_class=None,
                 global_options=None,
                 check=True):
        """
        Initialize ``self``.

        TESTS::

            sage: C = IntegerListsLex(2, length=3)
            sage: C == loads(dumps(C))
            True
            sage: C == loads(dumps(C)) # this did fail at some point, really!
            True
            sage: C is loads(dumps(C)) # todo: not implemented
            True
            sage: C.cardinality().parent() is ZZ
            True
            sage: TestSuite(C).run()

            sage: IntegerListsLex(min_sum=Infinity).list()
            Traceback (most recent call last):
            ...
            TypeError: unable to coerce <class 'sage.rings.infinity.PlusInfinity'> to an integer
            sage: IntegerListsLex(min_sum=1.4).list()
            Traceback (most recent call last):
            ...
            TypeError: Attempt to coerce non-integral RealNumber to Integer
        """
        if category is None:
            category = EnumeratedSets().Finite()

        self._check = check

        if n is not None:
            min_sum = n
            max_sum = n
        self._min_sum = ZZ(min_sum)
        self._max_sum = ZZ(max_sum) if max_sum != Infinity else Infinity

        if length is not None:
            min_length = length
            max_length = length
        self._min_length = max(ZZ(min_length), 0)
        self._max_length = ZZ(max_length) if max_length != Infinity else Infinity

        self._min_slope = ZZ(min_slope) if min_slope != -Infinity else -Infinity
        self._max_slope = ZZ(max_slope) if max_slope !=  Infinity else Infinity

        self._min_part = ZZ(min_part)
        if self._min_part < 0:
            raise NotImplementedError("strictly negative min_part")
        self._max_part = ZZ(max_part) if max_part != Infinity else Infinity

        # self._floor_or_ceiling_is_function will be set to ``True``
        # if a function is given as input for floor or ceiling; in
        # this case a warning will be emitted, unless the user sets
        # check=False. See the documentation.
        self._floor_or_ceiling_is_function = False
        if floor is None:
            floor = 0
        elif isinstance(floor, (list, tuple)):
            floor = tuple(ZZ(i) for i in floor)
            if not all(i >= 0 for i in floor):
                raise NotImplementedError("negative parts in floor={}".format(floor))
        elif callable(floor):
            self._floor_or_ceiling_is_function = True
        else:
            raise TypeError("floor should be a list, tuple, or function")
        self._floor = Envelope(floor, sign=-1,
                               min_part=  self._min_part,  max_part= self._max_part,
                               min_slope= self._min_slope, max_slope=self._max_slope,
                               min_length=self._min_length)

        if ceiling is None:
            ceiling = Infinity
        elif isinstance(ceiling, (list, tuple)):
            ceiling = tuple(ZZ(i) if i != Infinity else Infinity
                            for i in ceiling)
            if not all(i >= 0 for i in ceiling):
                raise NotImplementedError("negative parts in ceiling={}".format(ceiling))
        elif callable(ceiling):
            self._floor_or_ceiling_is_function = True
        else:
            raise ValueError("Unable to parse value of parameter ceiling")
        self._ceiling = Envelope(ceiling, sign=1,
                               min_part=  self._min_part,  max_part= self._max_part,
                               min_slope= self._min_slope, max_slope=self._max_slope,
                               min_length=self._min_length)

        if name is not None:
            self.rename(name)

        if self._floor_or_ceiling_is_function and self._check:
            from warnings import warn
            warn("""
A function has been given as input of the floor=[...] or ceiling=[...]
arguments of IntegerListsLex. Please see the documentation for the caveats.
If you know what you are doing, you can set check=False to skip this warning.""")

        # Customization of the class and constructor for the elements

        # We set the following attribute to True if the element
        # constructor is known to be safe and does not claim ownership
        # on the input list. In this case, we can save a copy in Iter.next.
        self._element_constructor_is_copy_safe = False
        if element_class is not None:
            self.Element = element_class
        if element_constructor is not None:
            if element_constructor is list or element_constructor is tuple:
                self._element_constructor_is_copy_safe = True
        elif issubclass(self.Element, ClonableArray):
            # Not all element class support check=False
            element_constructor = self._element_constructor_nocheck
            self._element_constructor_is_copy_safe = True
        if global_options is not None:
            self.global_options = global_options

        Parent.__init__(self, element_constructor=element_constructor,
                        category=category)

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
           is ``self._ceiling(i)`` is finite.

        #. There exists a global upper bound on the length.

        Failures for 1. are detected and reported later, during the
        iteration, namely the first time a prefix including the `i`-th
        part is explored.

        This method therefore focuses on 2., namely trying to prove
        the existence of an upper bound on the length. It may fail
        to do so even when the set is actually finite.

        OUTPUT: ``None`` if this method finds a proof that there
        exists an upper bound on the length. Otherwise a
        ``ValueError`` is raised.

        EXAMPLES::

            sage: L = IntegerListsLex(4, max_length=4)
            sage: L._check_finiteness()

        The following example is infinite::

            sage: L = IntegerListsLex(4)
            sage: L._check_finiteness()
            Traceback (most recent call last):
            ...
            ValueError: Could not prove that the specified constraints yield a finite set

        Indeed::

            sage: it = iter(IntegerListsLex(4, check=False))
            sage: for _ in range(10): print next(it)
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

        Unless ``check=False`, :meth:`_check_finiteness` is called as
        soon as an iteration is attempted::

            sage: iter(L)
            Traceback (most recent call last):
            ...
            ValueError: Could not prove that the specified constraints yield a finite set

        Some other infinite examples::

            sage: L = IntegerListsLex(ceiling=[0], min_slope=1, max_slope=2)
            sage: L.list()
            Traceback (most recent call last):
            ...
            ValueError: Could not prove that the specified constraints yield a finite set

            sage: L = IntegerListsLex(ceiling=[0], min_slope=1, max_slope=1)
            sage: L.list()
            Traceback (most recent call last):
            ...
            ValueError: Could not prove that the specified constraints yield a finite set

            sage: IntegerListsLex(ceiling=[0], min_slope=1, max_slope=1).list()
            Traceback (most recent call last):
            ...
            ValueError: Could not prove that the specified constraints yield a finite set

        The following example is actually finite, but not detected as such::

            sage: IntegerListsLex(7, floor=[4], max_part=4, min_slope=-1).list()
            Traceback (most recent call last):
            ...
            ValueError: Could not prove that the specified constraints yield a finite set

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
            ValueError: Could not prove that the specified constraints yield a finite set

            sage: IntegerListsLex(7, ceiling=lambda i:0).list()
            Traceback (most recent call last):
            ...
            ValueError: Could not prove that the specified constraints yield a finite set

        The next example shows a case that is finite because we remove
        trailing zeroes::

            sage: list(IntegerListsLex(ceiling=[0], max_slope=0))
            [[]]
            sage: L = IntegerListsLex(ceiling=[1], min_slope=1, max_slope=1)
            sage: L.list()
            Traceback (most recent call last):
            ...
            ValueError: Could not prove that the specified constraints yield a finite set

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
        if self._max_length < Infinity:
            return
        if self._max_sum < self._min_sum:
            return
        if self._min_slope > self._max_slope:
            return
        if self._max_slope < 0:
            return
        if self._ceiling.limit() < self._floor.limit():
            return
        if self._ceiling.limit() == 0:
            # This assumes no trailing zeroes
            return
        if self._min_slope > 0 and self._ceiling.limit() < Infinity:
            return

        # Compute a lower bound on the sum of floor(i) for i=1 to infinity
        if self._floor.limit() > 0 or self._min_slope > 0:
            floor_sum_lower_bound = Infinity
        elif self._floor.limit_start() < Infinity:
            floor_sum_lower_bound = sum(self._floor(i) for i in range(self._floor.limit_start()))
        else:
            floor_sum_lower_bound = 0
        if floor_sum_lower_bound > 0 and self._min_slope >= 0:
            floor_sum_lower_bound = Infinity

        if self._max_sum < floor_sum_lower_bound:
            return
        if self._max_sum == floor_sum_lower_bound and self._max_sum < Infinity:
                # This assumes no trailing zeroes
                return

        # Variant on ceiling.limit() ==0 where we actually discover that the ceiling limit is 0
        if self._max_slope == 0 and \
           (self._max_sum < Infinity or
            (self._ceiling.limit_start() < Infinity and
             any(self._ceiling(i) == 0 for i in range(self._ceiling.limit_start()+1)))):
            return

        limit_start = max(self._ceiling.limit_start(), self._floor.limit_start())
        if limit_start < Infinity:
            for i in range(limit_start+1):
                if self._ceiling(i) < self._floor(i):
                    return

        raise ValueError("Could not prove that the specified constraints yield a finite set")


    @staticmethod
    def _list_function(l, default):
        r"""
        Generate a function on the nonnegative integers from input.

        This method generates a function on the nonnegative integers
        whose values are taken from ``l`` when the input is a valid index
        in the list ``l``, and has a default value ``default`` otherwise.

        INPUT:

        - ``l`` -- a list to use as a source of values

        - ``default`` -- a default value to use for indices outside of the list

        OUTPUT:

        A function on the nonnegative integers.

        EXAMPLES::

            sage: C = IntegerListsLex(2, length=3)
            sage: f = C._list_function([1,2], Infinity)
            sage: f(1)
            2
            sage: f(3)
            +Infinity
        """
        return lambda i: l[i] if i < len(l) else default

    def __eq__(self, other):
        r"""
        Return whether ``self == other``.

        EXAMPLES::

            sage: C = IntegerListsLex(2, length=3)
            sage: D = IntegerListsLex(2, length=3); L = D.list();
            sage: E = IntegerListsLex(2, min_length=3)
            sage: F = IntegerListsLex(2, length=3, element_constructor=list)
            sage: G = IntegerListsLex(4, length=3)
            sage: C == C
            True
            sage: C == D
            True
            sage: C == E
            False
            sage: C == F
            False
            sage: C == None
            False
            sage: C == G
            False

        This is a minimal implementation enabling pickling tests. It
        is safe, but one would want the two following objects to be
        detected as equal::

            sage: C = IntegerListsLex(2, ceiling=[1,1,1])
            sage: D = IntegerListsLex(2, ceiling=[1,1,1])
            sage: C == D
            False

        TESTS:

        This used to fail due to poor equality testing. See
        :trac:`17979`, comment 433::

            sage: DisjointUnionEnumeratedSets(Family([2,2],
            ....:     lambda n: IntegerListsLex(n, length=2))).list()
            [[2, 0], [1, 1], [0, 2], [2, 0], [1, 1], [0, 2]]
            sage: DisjointUnionEnumeratedSets(Family([2,2],
            ....:     lambda n: IntegerListsLex(n, length=1))).list()
            [[2], [2]]
        """
        if self.__class__ != other.__class__:
            return False
        for key in ["_min_length", "_max_length", "_floor", "_ceiling", "_min_part", "_max_part", "_min_sum", "_max_sum", "Element"]:
            if getattr(self, key) != getattr(other, key):
                return False
        a = self._element_constructor
        b = other._element_constructor
        if ismethod(a):
            a = a.im_func
        if ismethod(b):
            b = b.im_func
        return a == b

    def __ne__(self, other):
        r"""
        Return whether ``self != other``.

        EXAMPLES::

            sage: C = IntegerListsLex(2, length=3)
            sage: D = IntegerListsLex(2, length=3); L = D.list();
            sage: E = IntegerListsLex(2, max_length=3)
            sage: C != D
            False
            sage: C != E
            True
        """
        return not self == other

    def _repr_(self):
        """
        Return the name of this enumerated set.

        EXAMPLES::

            sage: C = IntegerListsLex(2, length=3)
            sage: C # indirect doctest
            Integer lists of sum 2 satisfying certain constraints

            sage: C = IntegerListsLex(2, length=3, name="A given name")
            sage: C
            A given name
        """
        if self._min_sum == self._max_sum:
            return "Integer lists of sum {} satisfying certain constraints".format(self._min_sum)
        elif self._max_sum == Infinity:
            if self._min_sum == 0:
                return "Integer lists with arbitrary sum satisfying certain constraints"
            else:
                return "Integer lists of sum at least {} satisfying certain constraints".format(self._min_sum)
        else:
            return "Integer lists of sum between {} and {} satisfying certain constraints".format(self._min_sum,self._max_sum)

    def __contains__(self, comp):
        """
        Return ``True`` if ``comp`` meets the constraints imposed by the arguments.

        EXAMPLES::

            sage: C = IntegerListsLex(n=2, max_length=3, min_slope=0)
            sage: all([l in C for l in C])
            True
        """
        if len(comp) < self._min_length or len(comp) > self._max_length:
            return False
        n = sum(comp)
        if n < self._min_sum or n > self._max_sum:
            return False
        for i in range(len(comp)):
            if comp[i] < self._floor(i):
                return False
            if comp[i] > self._ceiling(i):
                return False
        for i in range(len(comp)-1):
            slope = comp[i+1] - comp[i]
            if slope < self._min_slope or slope > self._max_slope:
                return False
        return True


    def __iter__(self):
        """
        Return an iterator for the elements of ``self``.

        EXAMPLES::

            sage: C = IntegerListsLex(2, length=3)
            sage: list(C)     # indirect doctest
            [[2, 0, 0], [1, 1, 0], [1, 0, 1], [0, 2, 0], [0, 1, 1], [0, 0, 2]]
        """
        if self._check:
            self._check_finiteness()
        return IntegerListsLexIter(self)

    def _element_constructor_nocheck(self, l):
        r"""
        A variant of the standard element constructor that passes
        ``check=False`` to the element class.

        EXAMPLES::

            sage: L = IntegerListsLex(4, max_slope=0)
            sage: L._element_constructor_nocheck([1,2,3])
            [1, 2, 3]

        When relevant, this is assigned to
        ``self._element_constructor`` by :meth:`__init__`, to avoid
        overhead when constructing elements from trusted data in the
        iterator::

            sage: L._element_constructor
            <bound method IntegerListsLex._element_constructor_nocheck of ...>
            sage: L._element_constructor([1,2,3])
            [1, 2, 3]
        """
        return self.element_class(self, l, check=False)

    class Element(ClonableArray):
        """
        Element class for :class:`IntegerListsLex`.
        """
        def check(self):
            """
            Check to make sure this is a valid element in its
            :class:`IntegerListsLex` parent.

            EXAMPLES::

                sage: C = IntegerListsLex(4)
                sage: C([4]).check()
                True
                sage: C([5]).check()   # not implemented
                False
            """
            return self.parent().__contains__(self)


# Constants for IntegerListsLexIter._next_state
LOOKAHEAD = 5
PUSH      = 4
ME        = 3
DECREASE  = 2
POP       = 1
STOP      = 0

class IntegerListsLexIter:
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

    - ``_parent`` -- the :class:`IntegerListsLex` object this is
      iterating on;

    - ``_current_list`` -- the list corresponding to the current
      node of the tree;

    - ``_j`` -- the index of the last element of ``_current_list``:
      ``self._j == len(self._current_list) - 1``;

    - ``_current_sum`` -- the sum of the parts of ``_current_list``;

    - ``_search_ranges`` -- a list of same length as
      ``_current_list``: the range for each part.

    Furthermore, we assume that there is no obvious contradiction
    in the contraints:

    - ``self._parent._min_length <= self._parent._max_length``;
    - ``self._parent._min_slope <= self._parent._max_slope``
      unless ``self._parent._min_length <= 1``.

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
    def __init__(self, parent):
        """
        TESTS::

            sage: from sage.combinat.integer_list import IntegerListsLexIter
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
        self._parent = parent

        self._search_ranges = []
        self._current_list = []
        self._j = -1     # index of last element of _current_list
        self._current_sum = 0     # sum of parts in _current_list

        # Make sure that some invariants are respected in the iterator
        if parent._min_length <= parent._max_length and \
           (parent._min_slope <= parent._max_slope or parent._min_length <= 1):
            self._next_state = PUSH
        else:
            self._next_state = STOP

    def __iter__(self):
        """
        Return ``self`` as per the iterator protocol.

        EXAMPLES::

            sage: from sage.combinat.integer_list import IntegerListsLexIter
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

        OUTPUT: a boolean: whether the push succeeded

        EXAMPLES::

            sage: C = IntegerListsLex(2, length=3)
            sage: I = C.__iter__()
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
        p = self._parent
        max_sum = p._max_sum
        min_length = p._min_length
        max_length = p._max_length
        if  self._j+1 >= max_length:
            return False
        if self._j+1 >= min_length and self._current_sum == max_sum:
            # Cannot add trailing zeroes
            return False

        if self._j >= 0:
            prev = self._current_list[self._j]
        else:
            prev = None
        interval = self._m_interval(self._j+1, self._parent._max_sum - self._current_sum, prev)
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
            sage: I = C.__iter__()
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

    def next(self):
        r"""
        Return the next element in the iteration.

        EXAMPLES::

            sage: from sage.combinat.integer_list import IntegerListsLexIter
            sage: C = IntegerListsLex(2, length=3)
            sage: I = IntegerListsLexIter(C)
            sage: I.next()
            [2, 0, 0]
            sage: I.next()
            [1, 1, 0]
        """
        p = self._parent
        min_sum = p._min_sum
        max_length = p._max_length
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
                    return p._element_constructor(
                        self._current_list
                        if p._element_constructor_is_copy_safe
                        else self._current_list[:])

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
        Return whether the current list in the iteration variable ``self._current_list`` is a valid list.

        This method checks whether the sum of the parts in ``self._current_list``
        is in the right range, whether its length is in the
        required range, and whether there are trailing zeroes.  It does not check all of the
        necessary conditions to verify that an arbitrary list satisfies the
        constraints from the corresponding ``IntegerListsLex`` object, and should
        not be used except internally in the iterator class.

        EXAMPLES::

            sage: from sage.combinat.integer_list import IntegerListsLexIter
            sage: C = IntegerListsLex(2, length=3)
            sage: I = IntegerListsLexIter(C)
            sage: I._current_list
            []
            sage: I._internal_list_valid()
            False
            sage: I.next()
            [2, 0, 0]
            sage: I._current_list
            [2, 0, 0]
            sage: I._internal_list_valid()
            True
        """
        p = self._parent
        mu = self._current_list
        nu = self._current_sum
        l = self._j + 1
        good_sum = (nu >= p._min_sum and nu <= p._max_sum)
        good_length = (l >= p._min_length and l <= p._max_length)
        no_trailing_zeros = (l <= max(p._min_length,0) or mu[-1] != 0)
        return good_sum and good_length and no_trailing_zeros

    def _m_interval(self, i, max_sum, prev=None):
        r"""
        Return coarse lower and upper bounds for the part ``m`` at position ``i``.

        INPUT:

        - ``i`` -- a nonnegative integer (position)

        - ``max_sum`` -- a nonnegative integer or ``+oo``

        - ``prev`` -- a nonnegative integer or ``None``

        Return coarse lower and upper bounds for the value ``m``
        of the part at position ``i`` so that there could exists
        some list suffix `v_i,\ldots,v_k` of sum bounded by
        ``max_sum`` and satisfying the floor and upper bound
        constraints. If ``prev`` is specified, then the slope
        conditions between ``v[i-1]=prev`` and ``v[i]=m`` should
        also be satisfied.

        Additionally, this raises an error if it can be detected
        that some part is neither directly nor indirectly bounded
        above, which implies that the constraints possibly do not allow for
        an inverse lexicographic iterator.

        OUTPUT:

        A tuple of two integers ``(lower_bound, upper_bound)``.

        EXAMPLES::

            sage: from sage.combinat.integer_list import IntegerListsLexIter
            sage: C = IntegerListsLex(2, length=3)
            sage: I = IntegerListsLexIter(C)
            sage: I._m_interval(1,2)
            (0, 2)

        The second part is not bounded above, hence we can not
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
        p = self._parent

        lower_bound = max(0, p._floor(i))
        upper_bound = min(max_sum, p._ceiling(i))
        if prev != None:
            lower_bound = max(lower_bound, prev + p._min_slope)
            upper_bound = min(upper_bound, prev + p._max_slope)

        ## check for infinite upper bound, in case max_sum is infinite
        if p._check and upper_bound == Infinity:
            # This assumes that there exists a valid list (which
            # is not yet always guaranteed). Then we just
            # discovered that part 'i' of this list can be made as
            # large as desired, which implies that `self._parent`
            # cannot be iterated in inverse lexicographic order
            raise ValueError("infinite upper bound for values of m")

        return (lower_bound, upper_bound)

    def _lookahead(self):
        r"""
        Return whether the current list can possibly be a prefix of a valid list.

        OUTPUT: ``False`` if it is guaranteed that the current
        list cannot be a prefix of a valid list and ``True``
        otherwise.

        EXAMPLES::

            sage: it = iter(IntegerListsLex(length=3, min_sum=2, max_sum=2))
            sage: it._current_list = [0,1]    # don't do this at home, kids
            sage: it._current_sum = 1
            sage: it._j = 1
            sage: it._lookahead()
            True

            sage: it = iter(IntegerListsLex(length=3, min_sum=3, max_sum=2))
            sage: it._current_list = [0,1]
            sage: it._current_sum = 1
            sage: it._j = 1
            sage: it._lookahead()
            False

            sage: it = iter(IntegerListsLex(min_length=2, max_part=0))
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
            sage: it = iter(IntegerListsLex(n, length=1))
            sage: it._current_list = [n-1]
            sage: it._current_sum = n-1
            sage: it._j = 0
            sage: it._lookahead()
            False

            sage: it = iter(IntegerListsLex(n=3, min_part=2, min_sum=3, max_sum=3))
            sage: it._current_list = [2]
            sage: it._current_sum = 2
            sage: it._j = 0
            sage: it._lookahead()
            False

        ALGORITHM:

        Let ``j=self._j`` be the position of the last part `m` of
        ``self._current_list``. The current algorithm computes,
        for `k=j,j+1,\ldots`, a lower bound `l_k` and an upper
        bound `u_k` for `v_0+\dots+v_k`, and stops if none of the
        invervals `[l_k, u_k]` intersect ``[min_sum, max_sum]``.

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

            sage: it = iter(IntegerListsLex(1, min_length=2, min_slope=0, max_slope=0, min_sum=1, max_sum=1))
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
        min_sum = self._parent._min_sum - (self._current_sum-m)
        max_sum = self._parent._max_sum - (self._current_sum-m)

        if min_sum > max_sum:
            return False

        p = self._parent

        # Beware that without slope conditions, the functions below
        # currently forget about the value m at k!
        lower_envelope = self._parent._floor.adapt(m,j)
        upper_envelope = self._parent._ceiling.adapt(m,j)
        lower = 0    # The lower bound `l_k`
        upper = 0    # The upper bound `u_k`

        assert j >= 0
        # get to smallest valid number of parts
        for k in range(j, p._min_length-1):
            # We are looking at lists `v_j,...,v_k`
            lo = m if k == j else lower_envelope(k)
            up = m if k == j else upper_envelope(k)
            if lo > up:
                return False
            lower += lo
            upper += up

        if j < p._min_length and min_sum <= upper and lower <= max_sum:
            # There could exist a valid list `v_j,\dots,v_{min_length-1}`
            return True

        k = max(p._min_length-1,j)
        # Check if any of the intervals intersect the target interval
        while k < p._max_length:
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

            if (p._max_slope <= 0 and up <= 0) or \
               (p._ceiling.limit() == 0 and k > p._ceiling.limit_start()):
                # This implies v_l=0 for l>=k: that is we would be generating
                # a list with trailing zeroes
                return False

            if min_sum <= upper and lower <= max_sum:
                # There could exist a valid list `v_j,\dots,v_k`
                return True

            k += 1

        return False


class Envelope(object):
    """
    The (currently approximated) upper (lower) envelope of a function
    under the specified constraints.

    INPUT:

    - ``f`` -- a function, list, or tuple; if ``f`` is a list, it is
      considered as the function ``f(i)=f[i]``, completed for larger
      `i` with ``f(i)=max_part``.

    - ``min_part``, ``max_part``, ``min_slope``, ``max_slope``, ...
      as for :class:`IntegerListsLex` (please consult for details).

    - ``sign`` -- (+1 or -1) multiply the input values with ``sign``
      and multiply the output with ``sign``. Setting this to `-1` can
      be used to implement a lower envelope.

    The *upper envelope* `U(f)` of `f` is the (pointwise) largest
    function which is bounded above by `f` and satisfies the
    ``max_part`` and ``max_slope`` conditions. Furthermore, for
    ``i,i+1<min_length``, the upper envelope also satisfies the
    ``min_slope`` condition.

    Upon computing `U(f)(i)`, all the previous values
    for `j\leq i` are computed and cached; in particular `f(i)` will
    be computed at most once for each `i`.

    .. TODO::

        - This class is a good candidate for Cythonization, especially
          to get the critical path in ``__call__`` super fast.

        - To get full envelopes, we would want both the ``min_slope``
          and ``max_slope`` conditions to always be satisfied. This is
          only properly defined for the restriction of `f` to a finite
          interval `0,..,k`, and depends on `k`.

        - This is the core "data structure" of
          ``IntegerListsLex``. Improving the lookahead there
          essentially depends on having functions with a good
          complexity to compute the area below an envelope; and in
          particular how it evolves when increasing the length.

    EXAMPLES::

        sage: from sage.combinat.integer_list import Envelope

    Trivial upper and lower envelopes::

        sage: f = Envelope([3,2,2])
        sage: [f(i) for i in range(10)]
        [3, 2, 2, inf, inf, inf, inf, inf, inf, inf]
        sage: f = Envelope([3,2,2], sign=-1)
        sage: [f(i) for i in range(10)]
        [3, 2, 2, 0, 0, 0, 0, 0, 0, 0]

    A more interesting lower envelope::

        sage: f = Envelope([4,1,5,3,5], sign=-1, min_part=2, min_slope=-1)
        sage: [f(i) for i in range(10)]
        [4, 3, 5, 4, 5, 4, 3, 2, 2, 2]

    Currently, adding ``max_slope`` has no effect::

        sage: f = Envelope([4,1,5,3,5], sign=-1, min_part=2, min_slope=-1, max_slope=0)
        sage: [f(i) for i in range(10)]
        [4, 3, 5, 4, 5, 4, 3, 2, 2, 2]

    unless ``min_length`` is large enough::

        sage: f = Envelope([4,1,5,3,5], sign=-1, min_part=2, min_slope=-1, max_slope=0, min_length=2)
        sage: [f(i) for i in range(10)]
        [4, 3, 5, 4, 5, 4, 3, 2, 2, 2]

        sage: f = Envelope([4,1,5,3,5], sign=-1, min_part=2, min_slope=-1, max_slope=0, min_length=4)
        sage: [f(i) for i in range(10)]
        [5, 5, 5, 4, 5, 4, 3, 2, 2, 2]

        sage: f = Envelope([4,1,5,3,5], sign=-1, min_part=2, min_slope=-1, max_slope=0, min_length=5)
        sage: [f(i) for i in range(10)]
        [5, 5, 5, 5, 5, 4, 3, 2, 2, 2]

    A non trivial upper envelope::

        sage: f = Envelope([9,1,5,4], max_part=7, max_slope=2)
        sage: [f(i) for i in range(10)]
        [7, 1, 3, 4, 6, 7, 7, 7, 7, 7]

    TESTS::

        sage: f = Envelope(3, min_slope=1)
        sage: [f(i) for i in range(10)]
        [3, 3, 3, 3, 3, 3, 3, 3, 3, 3]

        sage: f = Envelope(3, min_slope=1, min_length=5)
        sage: [f(i) for i in range(10)]
        [-1, 0, 1, 2, 3, 3, 3, 3, 3, 3]

        sage: f = Envelope(3, sign=-1, min_slope=1)
        sage: [f(i) for i in range(10)]
        [3, 4, 5, 6, 7, 8, 9, 10, 11, 12]

        sage: f = Envelope(3, sign=-1, max_slope=-1, min_length=4)
        sage: [f(i) for i in range(10)]
        [6, 5, 4, 3, 3, 3, 3, 3, 3, 3]

    .. automethod:: __init__
    """
    def __init__(self, f,
                 min_part=0, max_part=Infinity,
                 min_slope=-Infinity, max_slope=Infinity,
                 min_length=0, max_length=Infinity, sign=1):
        r"""
        Initialize this envelope.

        TESTS::

            sage: from sage.combinat.integer_list import Envelope
            sage: f = Envelope(3, sign=-1, max_slope=-1, min_length=4)
            sage: f.__dict__
            {'_f': The constant function (...) -> inf,
             '_f_limit_start': 0,
             '_max_part': -3,
             '_max_slope': inf,
             '_min_slope': 1,
             '_precomputed': [-6, -5, -4, -3],
             '_sign': -1}
            sage: TestSuite(f).run(skip="_test_pickling")
            sage: Envelope(3, sign=1/3, max_slope=-1, min_length=4)
            Traceback (most recent call last):
            ...
            TypeError: no conversion of this rational to integer
            sage: Envelope(3, sign=-2, max_slope=-1, min_length=4)
            Traceback (most recent call last):
            ...
            ValueError: sign should be +1 or -1
        """
        # self._sign = sign for the output values (the sign change for
        # f is handled here in __init__)
        self._sign = ZZ(sign)
        if self._sign == 1:
            self._max_part = max_part
            self._min_slope = min_slope
            self._max_slope = max_slope
            if max_part == 0:
                # This uses that all entries are nonnegative.
                # This is not for speed optimization but for
                # setting the limit start and avoid hangs.
                # See #17979: comment 389
                f = 0
        elif self._sign == -1:
            self._max_part = -min_part
            self._min_slope = -max_slope
            self._max_slope = -min_slope
        else:
            raise ValueError("sign should be +1 or -1")

        # Handle different types of f and multiply f with sign
        if f == Infinity or f == -Infinity or f in ZZ:
            limit_start = 0
            self._max_part = min(self._sign * f, self._max_part)
            f = ConstantFunction(Infinity)
        elif isinstance(f, (list, tuple)):
            limit_start = len(f)
            f_tab = [self._sign * i for i in f]
            f = lambda k: f_tab[k] if k < len(f_tab) else Infinity
        else:
            g = f
            f = lambda k: self._sign * g(k)
            # At this point, this is not really used
            limit_start = Infinity

        self._f = f
        # For i >= limit_start, f is constant
        # This does not necessarily means that self is constant!
        self._f_limit_start = limit_start
        self._precomputed = []

        if min_length > 0:
            self(min_length-1)
            for i in range(min_length-1,0,-1):
                self._precomputed[i-1] = min(self._precomputed[i-1], self._precomputed[i] - self._min_slope)

    def __eq__(self, other):
        r"""
        Return whether ``self == other``.

        This is a minimal implementation enabling pickling tests.

        EXAMPLES::

            sage: from sage.combinat.integer_list import Envelope
            sage: f = Envelope([3,2,2])
            sage: g = Envelope([3,2,2])
            sage: h = Envelope([3,2,2], min_part=2)
            sage: f == f, f == h, f == None
            (True, False, False)

        This would be desirable::

            sage: f == g          # todo: not implemented
            True
        """
        return self.__class__ == other.__class__ and self.__dict__ == other.__dict__

    def __ne__(self, other):
        r"""
        Return whether ``self != other``.

        This is a minimal implementation enabling pickling tests.

        EXAMPLES::

            sage: from sage.combinat.integer_list import Envelope
            sage: f = Envelope([3,2,2])
            sage: g = Envelope([3,2,2])
            sage: h = Envelope([3,2,2], min_part=2)
            sage: f != f, f != h, f != None
            (False, True, True)

        This would be desirable::

            sage: f != g           # todo: not implemented
            False
        """
        return not self == other

    def limit_start(self):
        """
        Return from which `i` on the bound returned by ``limit`` holds.

        .. SEEALSO:: :meth:`limit` for the precise specifications.

        EXAMPLES::

            sage: from sage.combinat.integer_list import Envelope
            sage: Envelope([4,1,5]).limit_start()
            3
            sage: Envelope([4,1,5], sign=-1).limit_start()
            3

            sage: Envelope([4,1,5], max_part=2).limit_start()
            3

            sage: Envelope(4).limit_start()
            0
            sage: Envelope(4, sign=-1).limit_start()
            0

            sage: Envelope(lambda x: 3).limit_start() == Infinity
            True
            sage: Envelope(lambda x: 3, max_part=2).limit_start() == Infinity
            True

            sage: Envelope(lambda x: 3, sign=-1, min_part=2).limit_start() == Infinity
            True

        """
        return self._f_limit_start

    def limit(self):
        """
        Return a bound on the limit of ``self``.

        OUTPUT: a nonnegative integer or `\infty`

        This returns some upper bound for the accumulation points of
        this upper envelope. For a lower envelope, a lower bound is
        returned instead.

        In particular this gives a bound for the value of ``self`` at
        `i` for `i` large enough. Special case: for a lower envelop,
        and when the limit is `\infty`, the envelope is guaranteed to
        tend to `\infty` instead.

        When ``s=self.limit_start()`` is finite, this bound is
        guaranteed to be valid for `i>=s`.

        Sometimes it's better to have a loose bound that starts early;
        sometimes the converse holds. At this point which specific
        bound and starting point is returned is not set in stone, in
        order to leave room for later optimizations.

        EXAMPLES::

            sage: from sage.combinat.integer_list import Envelope
            sage: Envelope([4,1,5]).limit()
            inf
            sage: Envelope([4,1,5], max_part=2).limit()
            2
            sage: Envelope([4,1,5], max_slope=0).limit()
            1
            sage: Envelope(lambda x: 3, max_part=2).limit()
            2

        Lower envelopes::

            sage: Envelope(lambda x: 3, min_part=2, sign=-1).limit()
            2
            sage: Envelope([4,1,5], min_slope=0, sign=-1).limit()
            5
            sage: Envelope([4,1,5], sign=-1).limit()
            0

        .. SEEALSO:: :meth:`limit_start`
        """
        if self.limit_start() < Infinity and self._max_slope <= 0:
            return self(self.limit_start())
        else:
            return self._max_part * self._sign

    def __call__(self, k):
        """
        Return the value of this envelope at `k`.

        EXAMPLES::

            sage: from sage.combinat.integer_list import Envelope
            sage: f = Envelope([4,1,5,3,5])
            sage: f.__call__(2)
            5
            sage: [f(i) for i in range(10)]
            [4, 1, 5, 3, 5, inf, inf, inf, inf, inf]

        .. NOTE::

            See the documentation of :class:`Envelope` for tests and
            examples.
        """
        if k >= len(self._precomputed):
            for i in range(len(self._precomputed), k+1):
                value = min(self._f(i), self._max_part)
                if i>0:
                    value = min(value, self._precomputed[i-1] + self._max_slope)
                self._precomputed.append(value)
        return self._precomputed[k] * self._sign

    def adapt(self, m, j):
        """
        Return this envelope adapted to an additional local constraint.

        INPUT:

        - ``m`` -- a nonnegative integer (starting value)

        - ``j`` -- a nonnegative integer (position)

        This method adapts this envelope to the additional local
        constraint imposed by having a part `m` at position `j`.
        Namely, this returns a function which computes, for any `i>j`,
        the minimum of the ceiling function and the value restriction
        given by the slope conditions.

        EXAMPLES::

            sage: from sage.combinat.integer_list import Envelope
            sage: f = Envelope(3)
            sage: g = f.adapt(1,1)
            sage: g is f
            True
            sage: [g(i) for i in range(10)]
            [3, 3, 3, 3, 3, 3, 3, 3, 3, 3]

            sage: f = Envelope(3, max_slope=1)
            sage: g = f.adapt(1,1)
            sage: [g(i) for i in range(10)]
            [0, 1, 2, 3, 3, 3, 3, 3, 3, 3]

        Note that, in both cases above, the adapted envelope is only
        guaranteed to be valid for `i>j`! This is to leave potential
        room in the future for sharing similar adapted envelopes::

            sage: g = f.adapt(0,0)
            sage: [g(i) for i in range(10)]
            [0, 1, 2, 3, 3, 3, 3, 3, 3, 3]

            sage: g = f.adapt(2,2)
            sage: [g(i) for i in range(10)]
            [0, 1, 2, 3, 3, 3, 3, 3, 3, 3]

            sage: g = f.adapt(3,3)
            sage: [g(i) for i in range(10)]
            [0, 1, 2, 3, 3, 3, 3, 3, 3, 3]

        Now with a lower envelope::

            sage: f = Envelope(1, sign=-1, min_slope=-1)
            sage: g = f.adapt(2,2)
            sage: [g(i) for i in range(10)]
            [4, 3, 2, 1, 1, 1, 1, 1, 1, 1]
            sage: g = f.adapt(1,3)
            sage: [g(i) for i in range(10)]
            [4, 3, 2, 1, 1, 1, 1, 1, 1, 1]
        """
        if self._max_slope == Infinity:
            return self
        m *= self._sign
        m = m - j * self._max_slope
        return lambda i: self._sign * min(m + i*self._max_slope, self._sign*self(i) )


def IntegerListsNN(**kwds):
    """
    Lists of nonnegative integers with constraints.

    This function returns the union of ``IntegerListsLex(n, **kwds)``
    where `n` ranges over all nonnegative integers.

    .. WARNING:: this function is likely to disappear in :trac:`17927`.

    EXAMPLES::

        sage: from sage.combinat.integer_list import IntegerListsNN
        sage: L = IntegerListsNN(max_length=3, max_slope=-1)
        sage: L
        Disjoint union of Lazy family (<lambda>(i))_{i in Non negative integer semiring}
        sage: it = iter(L)
        sage: for _ in range(20):
        ....:     print next(it)
        []
        [1]
        [2]
        [3]
        [2, 1]
        [4]
        [3, 1]
        [5]
        [4, 1]
        [3, 2]
        [6]
        [5, 1]
        [4, 2]
        [3, 2, 1]
        [7]
        [6, 1]
        [5, 2]
        [4, 3]
        [4, 2, 1]
        [8]
    """
    from sage.rings.semirings.non_negative_integer_semiring import NN
    from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
    return DisjointUnionEnumeratedSets(Family(NN, lambda i: IntegerListsLex(i, **kwds)))
