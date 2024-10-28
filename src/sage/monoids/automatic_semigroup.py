"""
Automatic Semigroups

Semigroups defined by generators living in an ambient semigroup and represented by an automaton.

AUTHORS:

- Nicolas M. Thiéry
- Aladin Virmaux
"""
# ****************************************************************************
#       Copyright (C) 2010-2015 Nicolas M. Thiéry
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import operator
from sage.misc.cachefunc import cached_method
from sage.categories.semigroups import Semigroups
from sage.categories.sets_cat import Sets
from sage.categories.monoids import Monoids
from sage.categories.groups import Groups
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element_wrapper import ElementWrapper
from sage.sets.family import Family
from sage.rings.integer import Integer
from sage.cpython.getattr import raw_getattr


class AutomaticSemigroup(UniqueRepresentation, Parent):
    r"""
    Semigroups defined by generators living in an ambient semigroup.

    This implementation lazily constructs all the elements of the
    semigroup, and the right Cayley graph relations between them, and
    uses the latter as an automaton.

    EXAMPLES::

        sage: from sage.monoids.automatic_semigroup import AutomaticSemigroup
        sage: R = IntegerModRing(12)
        sage: M = AutomaticSemigroup(Family({1: R(3), 2: R(5)}), one=R.one())
        sage: M in Monoids()
        True
        sage: M.one()
        1
        sage: M.one() in M
        True
        sage: g = M._generators; g
        Finite family {1: 3, 2: 5}
        sage: g[1]*g[2]
        3
        sage: M.some_elements()
        [1, 3, 5, 9]

        sage: M.list()
        [1, 3, 5, 9]

        sage: M.idempotents()
        [1, 9]

    As can be seen above, elements are represented by default the
    corresponding element in the ambient monoid. One can also represent
    the elements by their reduced word::

        sage: M.repr_element_method("reduced_word")
        sage: M.list()
        [[], [1], [2], [1, 1]]

    In case the reduced word has not yet been calculated, the element
    will be represented by the corresponding element in the ambient
    monoid::

        sage: R = IntegerModRing(13)
        sage: N = AutomaticSemigroup(Family({1: R(3), 2: R(5)}), one=R.one())
        sage: N.repr_element_method("reduced_word")
        sage: n = N.an_element()
        sage: n
        [1]
        sage: n*n
        9

    Calling :meth:`construct`, :meth:`cardinality`, or :meth:`list`,
    or iterating through the monoid will trigger its full construction
    and, as a side effect, compute all the reduced words. The order of
    the elements, and the induced choice of reduced word is currently
    length-lexicographic (i.e. the chosen reduced word is of minimal
    length, and then minimal lexicographically w.r.t. the order of the
    indices of the generators)::

        sage: M.cardinality()
        4
        sage: M.list()
        [[], [1], [2], [1, 1]]
        sage: g = M._generators

        sage: g[1]*g[2]
        [1]

        sage: g[1].transition(1)
        [1, 1]
        sage: g[1] * g[1]
        [1, 1]
        sage: g[1] * g[1] * g[1]
        [1]
        sage: g[1].transition(2)
        [1]
        sage: g[1] * g[2]
        [1]

        sage: [ x.lift() for x in M.list() ]
        [1, 3, 5, 9]

        sage: # needs sage.graphs
        sage: G = M.cayley_graph(side='twosided'); G
        Looped multi-digraph on 4 vertices
        sage: G.edges(sort=True, key=str)
        [([1, 1], [1, 1], (2, 'left')),
         ([1, 1], [1, 1], (2, 'right')),
         ([1, 1], [1], (1, 'left')),
         ([1, 1], [1], (1, 'right')),
         ([1], [1, 1], (1, 'left')),
         ([1], [1, 1], (1, 'right')),
         ([1], [1], (2, 'left')),
         ([1], [1], (2, 'right')),
         ([2], [1], (1, 'left')),
         ([2], [1], (1, 'right')),
         ([2], [], (2, 'left')),
         ([2], [], (2, 'right')),
         ([], [1], (1, 'left')),
         ([], [1], (1, 'right')),
         ([], [2], (2, 'left')),
         ([], [2], (2, 'right'))]
        sage: list(map(sorted, M.j_classes()))
        [[[1], [1, 1]], [[], [2]]]
        sage: M.j_classes_of_idempotents()
        [[[1, 1]], [[]]]
        sage: M.j_transversal_of_idempotents()
        [[1, 1], []]

        sage: list(map(attrcall('pseudo_order'), M.list()))                             # needs sage.graphs
        [[1, 0], [3, 1], [2, 0], [2, 1]]

    We can also use it to get submonoids from groups. We check that in the
    symmetric group, a transposition and a long cycle generate the whole group::

        sage: # needs sage.groups
        sage: G5 = SymmetricGroup(5)
        sage: N = AutomaticSemigroup(Family({1: G5([2,1,3,4,5]), 2: G5([2,3,4,5,1])}),
        ....:                        one=G5.one())
        sage: N.repr_element_method("reduced_word")
        sage: N.cardinality() == G5.cardinality()
        True
        sage: N.retract(G5((1,4,3,5,2)))
        [1, 2, 1, 2, 2, 1, 2, 1, 2, 2]
        sage: N.from_reduced_word([1, 2, 1, 2, 2, 1, 2, 1, 2, 2]).lift()
        (1,4,3,5,2)

    We can also create a semigroup of matrices, where we define the
    multiplication as matrix multiplication::

        sage: # needs sage.modules
        sage: M1 = matrix([[0,0,1],[1,0,0],[0,1,0]])
        sage: M2 = matrix([[0,0,0],[1,1,0],[0,0,1]])
        sage: M1.set_immutable()
        sage: M2.set_immutable()
        sage: def prod_m(x,y):
        ....:     z=x*y
        ....:     z.set_immutable()
        ....:     return z
        sage: Mon = AutomaticSemigroup([M1,M2], mul=prod_m,
        ....:                          category=Monoids().Finite().Subobjects())
        sage: Mon.cardinality()
        24
        sage: C = Mon.cayley_graph()                                                    # needs sage.graphs
        sage: C.is_directed_acyclic()                                                   # needs sage.graphs
        False

    Let us construct and play with the 0-Hecke Monoid::

        sage: # needs sage.graphs sage.modules
        sage: W = WeylGroup(['A',4]); W.rename('W')
        sage: ambient_monoid = FiniteSetMaps(W, action='right')
        sage: pi = W.simple_projections(length_increasing=True).map(ambient_monoid)
        sage: M = AutomaticSemigroup(pi, one=ambient_monoid.one()); M
        A submonoid of (Maps from W to itself) with 4 generators
        sage: M.repr_element_method("reduced_word")
        sage: sorted(M._elements_set, key=str)
        [[1], [2], [3], [4], []]
        sage: M.construct(n=10)
        sage: sorted(M._elements_set, key=str)
        [[1, 2], [1, 3], [1, 4], [1], [2, 1], [2, 3], [2], [3], [4], []]
        sage: elt = M.from_reduced_word([3,1,2,4,2])
        sage: M.construct(up_to=elt)
        sage: len(M._elements_set)
        36
        sage: M.cardinality()
        120

    We check that the 0-Hecke monoid is `J`-trivial and contains `2^4`
    idempotents::

        sage: len(M.idempotents())                                                      # needs sage.graphs sage.modules
        16
        sage: all(len(j) == 1 for j in M.j_classes())                                   # needs sage.graphs sage.modules
        True

    TESTS::

        sage: (g[1]).__hash__() == (g[1]*g[1]*g[1]).__hash__()
        True
        sage: g[1] == g[1]*g[1]*g[1]
        True
        sage: M.__class__                                                               # needs sage.graphs sage.modules
        <class 'sage.monoids.automatic_semigroup.AutomaticMonoid_with_category'>
        sage: TestSuite(M).run()                                                        # needs sage.graphs sage.modules

        sage: from sage.monoids.automatic_semigroup import AutomaticSemigroup
        sage: R = IntegerModRing(34)
        sage: M = AutomaticSemigroup(Family({1: R(3), 2: R(7)}), one=R.one())
        sage: M[3] in M
        True

    We need to pass in the ambient monoid to ``__init__`` to guarantee
    :class:`UniqueRepresentation` works properly::

        sage: R1 = IntegerModRing(12)
        sage: R2 = IntegerModRing(16)
        sage: M1 = AutomaticSemigroup(Family({1: R1(3), 2: R1(5)}), one=R1.one())
        sage: M2 = AutomaticSemigroup(Family({1: R2(3), 2: R2(5)}), one=R2.one())
        sage: M1 is M2
        False

    .. NOTE::

        Unlike what the name of the class may suggest, this currently
        implements only a subclass of automatic semigroups;
        essentially the finite ones. See :wikipedia:`Automatic_semigroup`.

    .. WARNING::

        :class:`AutomaticSemigroup` is designed primarily for finite
        semigroups. This property is not checked automatically (this
        would be too costly, if not undecidable). Use with care for an
        infinite semigroup, as certain features may require
        constructing all of it::

            sage: M = AutomaticSemigroup([2], category = Monoids().Subobjects()); M
            A submonoid of (Integer Ring) with 1 generators
            sage: M.retract(2)
            2
            sage: M.retract(3)   # not tested: runs forever trying to find 3
    """
    @staticmethod
    def __classcall_private__(cls, generators, ambient=None, one=None, mul=operator.mul, category=None):
        """
        Parse and straighten the arguments; figure out the category.

        TESTS::

            sage: from sage.monoids.automatic_semigroup import AutomaticSemigroup
            sage: R = IntegerModRing(9)
            sage: M = AutomaticSemigroup((), one=R.one())
            sage: M.ambient() == R
            True
            sage: AutomaticSemigroup((0,)).category()
            Join of Category of finitely generated semigroups
                and Category of subquotients of semigroups
                and Category of commutative magmas
                and Category of subobjects of sets
            sage: AutomaticSemigroup((0,), one=1).category()
            Join of Category of subquotients of monoids and
            Category of commutative monoids and
            Category of finitely generated semigroups and
            Category of subobjects of sets
            sage: AutomaticSemigroup((0,), one=0).category()
            Join of Category of commutative monoids and
            Category of finitely generated semigroups and
            Category of subquotients of semigroups and
            Category of subobjects of sets
            sage: AutomaticSemigroup((0,), mul=operator.add).category()
            Join of Category of semigroups and Category of subobjects of sets
            sage: AutomaticSemigroup((0,), one=0, mul=operator.add).category()
            Join of Category of monoids and Category of subobjects of sets

            sage: S5 = SymmetricGroup(5)                                                # needs sage.groups
            sage: AutomaticSemigroup([S5((1,2))]).category()                            # needs sage.groups
            Join of Category of finite groups and
            Category of subquotients of monoids and
            Category of finite finitely generated semigroups and
            Category of subquotients of finite sets and
            Category of subobjects of sets

        .. TODO::

            One would want a subsemigroup of a group to be
            automatically a subgroup (in ``Groups().Subobjects()``).
        """
        generators = Family(generators)
        if ambient is None:
            # Try to guess the ambient semigroup from the generators or the unit
            if generators.cardinality() > 0:
                ambient = generators.first().parent()
            elif one is not None:
                ambient = one.parent()
            else:
                raise ValueError("AutomaticSemigroup requires at least one generator or `one` to determine the ambient space")
        elif ambient not in Sets():
            raise ValueError("ambient (=%s) should be a set" % ambient)

        # if mul is not operator.mul  and category.is_subcategory(Monoids().Subobjects())  error

        if one is None and category is not None:
            if category.is_subcategory(Monoids().Subobjects()):
                one = ambient.one()
            elif category.is_subcategory(Monoids()):
                raise ValueError("For a monoid which is just a subsemigroup, the unit should be specified")

        # Try to determine the most specific category
        # This logic should be in the categories
        if mul is operator.mul:
            default_category = Semigroups().FinitelyGenerated()
            if one is not None and one == ambient.one():
                default_category = default_category.Unital()
            if ambient in Semigroups().Commutative():
                default_category = default_category.Commutative()
            if ambient in Groups().Finite():
                default_category = default_category & Groups()
        else:
            default_category = Sets()

        if ambient in Sets().Finite():
            default_category = default_category.Finite()

        default_category = default_category.Subobjects() & Semigroups()
        if one is not None:
            default_category = default_category.Unital()
            cls = AutomaticMonoid

        if category is None:
            category = default_category
        else:
            category = default_category & category
        return super().__classcall__(cls, generators, ambient=ambient,
                                     one=one, mul=mul, category=category)

    def __init__(self, generators, ambient, one, mul, category):
        """
        Initialize this semigroup.

        TESTS::

            sage: from sage.monoids.automatic_semigroup import AutomaticSemigroup
            sage: R = IntegerModRing(21)
            sage: M = AutomaticSemigroup(Family(()), one=R.one())
            sage: M.ambient() == R
            True
            sage: M = AutomaticSemigroup(Family(()))
            Traceback (most recent call last):
            ...
            ValueError: AutomaticSemigroup requires at least one generator or `one` to determine the ambient space
        """
        Parent.__init__(self, category=category)

        # Attributes for the multiplicative structure
        self._ambient = ambient
        self._mul = mul
        if one is not None:
            self._one = self._retract(one)
            self._one._reduced_word = []
        self._generators_in_ambient = generators
        self._generators = generators.map(self._retract)
        for e in self._generators:
            e._reduced_word = [self._generators.inverse_family()[e]]

        # Attributes for the lazy construction of the elements
        self._constructed = False
        self._done = 0
        if one is not None and self._one not in self._generators:
            self._elements = [self._one] + list(self._generators)
        else:
            self._elements = list(self._generators)
        self._elements_set = set(self._elements)
        self._iter = self.__init__iter()

        # Customization
        self._repr_element_method = "ambient"

    def _repr_(self):
        """
        Return the string representation for ``self``.

        EXAMPLES::

            sage: from sage.monoids.automatic_semigroup import AutomaticSemigroup
            sage: R = IntegerModRing(12)
            sage: AutomaticSemigroup(Family({1: R(3), 2: R(5)}), one=R.one())
            A submonoid of (Ring of integers modulo 12) with 2 generators
            sage: AutomaticSemigroup(Family({1: R(3), 2: R(5)}))
            A subsemigroup of (Ring of integers modulo 12) with 2 generators

            sage: AutomaticSemigroup(Family({1: R(3), 2: R(5)}), mul=operator.add)
            A semigroup with 2 generators
            sage: AutomaticSemigroup(Family({1: R(3), 2: R(5)}), mul=operator.add, one=R.zero())
            A semigroup with 2 generators

            sage: S5 = SymmetricGroup(5); S5.rename('S5')                               # needs sage.groups
            sage: AutomaticSemigroup(Family({1: S5((1,2))}),                            # needs sage.groups
            ....:                    category=Groups().Finite().Subobjects())
            A subgroup of (S5) with 1 generators
        """
        categories = [Groups(), Monoids(), Semigroups()]
        for category in categories:
            if self in category:
                typ = "A " + category._repr_object_names()[:-1]
        for category in [Groups(), Monoids(), Semigroups()]:
            if self.ambient() in category and self in category.Subobjects():
                typ = "A sub" + category._repr_object_names()[:-1]
                break
        if self._mul is operator.mul:
            of = " of (%s)" % self.ambient()
        else:
            of = ""

        return f"{typ}{of} with {len(self._generators)} generators"

    def repr_element_method(self, style='ambient'):
        """
        Set the representation of the elements of the monoid.

        INPUT:

        - ``style`` -- "ambient" or "reduced_word"

        EXAMPLES::

            sage: from sage.monoids.automatic_semigroup import AutomaticSemigroup
            sage: R = IntegerModRing(17)
            sage: M = AutomaticSemigroup(Family({1: R(3), 2: R(5)}), one=R.one())
            sage: M.list()
            [1, 3, 5, 9, 15, 8, 10, 11, 7, 6, 13, 16, 4, 14, 12, 2]
            sage: M.repr_element_method("reduced_word")
            sage: M.list()
            [[], [1], [2], [1, 1], [1, 2], [2, 2], [1, 1, 1], [1, 1, 2], [1, 2, 2],
             [2, 2, 2], [1, 1, 1, 1], [1, 1, 1, 2], [1, 1, 2, 2], [1, 1, 1, 1, 2],
             [1, 1, 1, 2, 2], [1, 1, 1, 1, 2, 2]]
        """
        self._repr_element_method = style

    def an_element(self):
        """
        Return the first given generator of ``self``.

        EXAMPLES::

            sage: from sage.monoids.automatic_semigroup import AutomaticSemigroup
            sage: R = IntegerModRing(16)
            sage: M = AutomaticSemigroup(Family({1: R(3), 2: R(5)}), one=R.one())
            sage: M.an_element()
            3
        """
        return self._generators.first()

    def ambient(self):
        """
        Return the ambient semigroup of ``self``.

        EXAMPLES::

            sage: from sage.monoids.automatic_semigroup import AutomaticSemigroup
            sage: R = IntegerModRing(12)
            sage: M = AutomaticSemigroup(Family({1: R(3), 2: R(5)}), one=R.one())
            sage: M.ambient()
            Ring of integers modulo 12

            sage: # needs sage.modules
            sage: M1 = matrix([[0,0,1],[1,0,0],[0,1,0]])
            sage: M2 = matrix([[0,0,0],[1,1,0],[0,0,1]])
            sage: M1.set_immutable()
            sage: M2.set_immutable()
            sage: def prod_m(x,y):
            ....:     z=x*y
            ....:     z.set_immutable()
            ....:     return z
            sage: Mon = AutomaticSemigroup([M1,M2], mul=prod_m)
            sage: Mon.ambient()
            Full MatrixSpace of 3 by 3 dense matrices over Integer Ring
        """
        return self._ambient

    def retract(self, ambient_element, check=True):
        """
        Retract an element of the ambient semigroup into ``self``.

        EXAMPLES::

            sage: # needs sage.groups
            sage: from sage.monoids.automatic_semigroup import AutomaticSemigroup
            sage: S5 = SymmetricGroup(5); S5.rename('S5')
            sage: M = AutomaticSemigroup(Family({1:S5((1,2)), 2:S5((1,2,3,4))}),
            ....:                        one=S5.one())
            sage: m = M.retract(S5((3,1))); m
            (1,3)
            sage: m.parent() is M
            True
            sage: M.retract(S5((4,5)), check=False)
            (4,5)
            sage: M.retract(S5((4,5)))
            Traceback (most recent call last):
            ...
            ValueError: (4,5) not in A subgroup of (S5) with 2 generators

        TESTS::

            sage: len(M._retract.cache.keys())                                          # needs sage.groups
            24
        """
        element = self._retract(ambient_element)
        if check:
            self.construct(up_to=ambient_element)
            if element not in self._elements_set:
                cache = self._retract.cache
                del cache[((ambient_element,), ())]
                raise ValueError(f"{ambient_element} not in {self}")
        return element

    @cached_method
    def _retract(self, ambient_element):
        r"""
        Retract an element of the ambient semigroup into ``self``.

        This is an internal method which does not check that
        ``ambient_element`` is indeed in this semigroup.

        EXAMPLES::

            sage: from sage.monoids.automatic_semigroup import AutomaticSemigroup
            sage: S5 = SymmetricGroup(5)                                                # needs sage.groups
            sage: S4 = AutomaticSemigroup(Family({1:S5((1,2)), 2:S5((1,2,3,4))}),       # needs sage.groups
            ....:                         one=S5.one())
            sage: S4._retract(S5((3,1)))                                                # needs sage.groups
            (1,3)

        No check is done::

            sage: S4._retract(S5((4,5)))                                                # needs sage.groups
            (4,5)
        """
        return self.element_class(self, ambient_element)

    def lift(self, x):
        """
        Lift an element of ``self`` into its ambient space.

        EXAMPLES::

            sage: from sage.monoids.automatic_semigroup import AutomaticSemigroup
            sage: R = IntegerModRing(15)
            sage: M = AutomaticSemigroup(Family({1: R(3), 2: R(5)}), one=R.one())
            sage: a = M.an_element()
            sage: a.lift() in R
            True
            sage: a.lift()
            3
            sage: [m.lift() for m in M]
            [1, 3, 5, 9, 0, 10, 12, 6]
        """
        assert x in self
        return x.lift()

    def semigroup_generators(self):
        """
        Return the family of generators of ``self``.

        EXAMPLES::

            sage: from sage.monoids.automatic_semigroup import AutomaticSemigroup
            sage: R = IntegerModRing(28)
            sage: M = AutomaticSemigroup(Family({1: R(3), 2: R(5)}))
            sage: M.semigroup_generators()
            Finite family {1: 3, 2: 5}
        """
        return self._generators
    gens = semigroup_generators

    def __init__iter(self):
        """
        Iterator on the elements of ``self``.

        EXAMPLES::

            sage: from sage.monoids.automatic_semigroup import AutomaticSemigroup
            sage: R = IntegerModRing(18)
            sage: M = AutomaticSemigroup([R(3), R(5)], one=R.one())
            sage: M.repr_element_method("reduced_word")
            sage: next(M.__iter__())
            []
            sage: list(M)
            [[], [0], [1], [0, 0], [0, 1], [1, 1], [1, 1, 1], [1, 1, 1, 1], [1,
            1, 1, 1, 1]]

        ALGORITHM:

            Breadth first search on the elements generated by the generators.
            The algorithm stops when all branches have been fully explored.
        """
        while self._done < len(self._elements):
            x = self._elements[self._done]
            for i in self._generators.keys():
                y = x.transition(i)
                if y in self._elements_set:
                    continue
                self._elements.append(y)
                self._elements_set.add(y)
                y._reduced_word = x.reduced_word() + [i]
                yield y
            self._done += 1
        self._constructed = True

    def __iter__(self):
        """
        Return iterator over elements of the semigroup.

        EXAMPLES::

            sage: from sage.monoids.automatic_semigroup import AutomaticSemigroup
            sage: R = IntegerModRing(5)
            sage: M = AutomaticSemigroup([R(3), R(4)], one=R.one())
            sage: I = M.__iter__()
            sage: next(I)
            1
            sage: M.list()
            [1, 3, 4, 2]
            sage: next(I)
            3
        """
        if self._constructed:
            return iter(self._elements)
        else:
            return self._iter_concurrent()

    def _iter_concurrent(self):
        """
        We need to take special care since several iterators may run
        concurrently.

        TESTS::

            sage: from sage.monoids.automatic_semigroup import AutomaticSemigroup
            sage: R = IntegerModRing(11)
            sage: M = AutomaticSemigroup(Family({1: R(3), 2: R(5)}), one=R.one())
            sage: f = iter(M)            # indirect doctest
            sage: g = iter(M)
            sage: next(f), next(g)
            (1, 1)
            sage: next(g), next(f)
            (3, 3)
            sage: next(f), next(g)
            (5, 5)
            sage: next(f), next(g)
            (9, 9)
            sage: h = iter(M)
            sage: next(h), next(h), next(h), next(h), next(h)
            (1, 3, 5, 9, 4)
            sage: next(f), next(g)
            (4, 4)
            sage: M._constructed
            False
            sage: next(f)
            Traceback (most recent call last):
            ...
            StopIteration
            sage: next(g)
            Traceback (most recent call last):
            ...
            StopIteration
            sage: next(h)
            Traceback (most recent call last):
            ...
            StopIteration
            sage: M._constructed
            True
        """
        i = 0
        # self._elements is never empty; so we are sure
        for x in self._elements:
            yield x
            # some other iterator/ method of the semigroup may have
            # been called before we move on to the next line
            i += 1
            if i == len(self._elements) and not self._constructed:
                try:
                    next(self._iter)
                except StopIteration:
                    # Don't allow StopIteration to bubble up from generator
                    # see PEP-479
                    break

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: from sage.monoids.automatic_semigroup import AutomaticSemigroup
            sage: R = IntegerModRing(12)
            sage: M = AutomaticSemigroup(Family({1: R(3), 2: R(5)}), one=R.one())
            sage: M.cardinality()
            4

        TESTS::

            sage: assert isinstance(M.cardinality(), Integer)  # This did fail at some point
        """
        if not self._constructed:
            self.construct()
        return Integer(len(self._elements))

    def list(self):
        """
        Return the list of elements of ``self``.

        EXAMPLES::

            sage: from sage.monoids.automatic_semigroup import AutomaticSemigroup
            sage: R = IntegerModRing(12)
            sage: M = AutomaticSemigroup(Family({1: R(3), 2: R(5)}), one=R.one())
            sage: M.repr_element_method("reduced_word")
            sage: M.list()
            [[], [1], [2], [1, 1]]

        TESTS::

            sage: assert isinstance(M.cardinality(), Integer)  # This did fail at some point
        """
        if not self._constructed:
            self.construct()
        return list(self._elements)

    def product(self, x, y):
        """
        Return the product of two elements in ``self``. It is done by
        retracting the multiplication in the ambient semigroup.

        EXAMPLES::

            sage: from sage.monoids.automatic_semigroup import AutomaticSemigroup
            sage: R = IntegerModRing(12)
            sage: M = AutomaticSemigroup(Family({1: R(3), 2: R(5)}), one=R.one())
            sage: a = M[1]
            sage: b = M[2]
            sage: a*b
            [1]
        """
        assert x in self
        assert y in self
        red = y._reduced_word
        if red is None:
            return self._retract(self._mul(x.lift(), y.lift()))

        for i in red:
            x = x.transition(i)
        return x

    def from_reduced_word(self, l):
        """
        Return the element of ``self`` obtained from the reduced word ``l``.

        INPUT:

        - ``l`` -- list of indices of the generators

        .. NOTE::

            We do not save the given reduced word ``l`` as an attribute of the
            element, as some elements above in the branches may have not been
            explored by the iterator yet.

        EXAMPLES::

            sage: # needs sage.groups
            sage: from sage.monoids.automatic_semigroup import AutomaticSemigroup
            sage: G4 = SymmetricGroup(4)
            sage: M = AutomaticSemigroup(Family({1:G4((1,2)), 2:G4((1,2,3,4))}),
            ....:                        one=G4.one())
            sage: M.from_reduced_word([2, 1, 2, 2, 1]).lift()
            (1,3)
            sage: M.from_reduced_word([2, 1, 2, 2, 1]) == M.retract(G4((3,1)))
            True
        """
        result = self.one()
        for i in l:
            result = result.transition(i)
        return result

    def construct(self, up_to=None, n=None):
        """
        Construct the elements of ``self``.

        INPUT:

        - ``up_to`` -- an element of ``self`` or of the ambient semigroup

        - ``n`` -- integer or ``None`` (default: ``None``)

        This construct all the elements of this semigroup, their
        reduced words, and the right Cayley graph. If `n` is
        specified, only the `n` first elements of the semigroup are
        constructed. If ``element`` is specified, only the elements up
        to ``ambient_element`` are constructed.

        EXAMPLES::

            sage: # needs sage.groups sage.modules
            sage: from sage.monoids.automatic_semigroup import AutomaticSemigroup
            sage: W = WeylGroup(['A',3]); W.rename('W')
            sage: ambient_monoid = FiniteSetMaps(W, action='right')
            sage: pi = W.simple_projections(length_increasing=True).map(ambient_monoid)
            sage: M = AutomaticSemigroup(pi, one=ambient_monoid.one()); M
            A submonoid of (Maps from W to itself) with 3 generators
            sage: M.repr_element_method("reduced_word")
            sage: sorted(M._elements_set, key=str)
            [[1], [2], [3], []]
            sage: elt = M.from_reduced_word([2,3,1,2])
            sage: M.construct(up_to=elt)
            sage: len(M._elements_set)
            19
            sage: M.cardinality()
            24
        """
        if self._constructed:
            return
        if n is not None:
            if up_to is not None:
                raise ValueError("Only one of the options `up_to` or `n` should be specified")
            i = len(self._elements)
            while i < n and not self._constructed:
                next(self._iter)
                i += 1
        elif up_to is not None:
            if up_to.parent() is self._ambient:
                up_to = self._retract(up_to)
                # TODO: remove up_to from the cache if not found at the end
            if up_to in self._elements_set:
                return
            for x in self._iter:
                if up_to is x:
                    return
        else:
            for x in self._iter:
                pass

    class Element(ElementWrapper):

        def __init__(self, ambient_element, parent):
            """
            TESTS::

                sage: from sage.monoids.automatic_semigroup import AutomaticSemigroup
                sage: R = IntegerModRing(21)
                sage: M = AutomaticSemigroup(Family([2]))
                sage: m = M(2); m
                2
                sage: type(m)
                <class 'sage.monoids.automatic_semigroup.AutomaticSemigroup_with_category.element_class'>
            """
            ElementWrapper.__init__(self, ambient_element, parent)
            self._reduced_word = None

        def reduced_word(self):
            r"""
            Return the length-lexicographic shortest word of ``self``.

            OUTPUT: list of indexes of the generators

            Obtaining the reduced word requires having constructed the
            Cayley graph of the semigroup up to ``self``. If this is
            not the case, an error is raised.

            EXAMPLES::

                sage: from sage.monoids.automatic_semigroup import AutomaticSemigroup
                sage: R = IntegerModRing(15)
                sage: M = AutomaticSemigroup(Family({1: R(3), 2: R(5)}), one=R.one())
                sage: M.construct()
                sage: for m in M: print((m, m.reduced_word()))
                (1, [])
                (3, [1])
                (5, [2])
                (9, [1, 1])
                (0, [1, 2])
                (10, [2, 2])
                (12, [1, 1, 1])
                (6, [1, 1, 1, 1])

            TESTS:

            We check that :issue:`19631` is fixed::

                sage: R = IntegerModRing(101)
                sage: M = AutomaticSemigroup(Family({1: R(3), 2: R(5)}), one=R.one())
                sage: e = M.from_reduced_word([1, 1, 1, 2, 2, 2])
                sage: e.reduced_word()
                [1, 1, 1, 2, 2, 2]
            """
            if self._reduced_word is None:
                self.parent().construct(up_to=self)
            return self._reduced_word

        def lift(self):
            """
            Lift the element ``self`` into its ambient semigroup.

            EXAMPLES::

                sage: from sage.monoids.automatic_semigroup import AutomaticSemigroup
                sage: R = IntegerModRing(18)
                sage: M = AutomaticSemigroup(Family({1: R(3), 2: R(5)}))
                sage: M.repr_element_method("reduced_word")
                sage: m = M.an_element(); m
                [1]
                sage: type(m)
                <class 'sage.monoids.automatic_semigroup.AutomaticSemigroup_with_category.element_class'>
                sage: m.lift()
                3
                sage: type(m.lift())
                <class 'sage.rings.finite_rings.integer_mod.IntegerMod_int'>
            """
            return self.value

        @cached_method
        def transition(self, i):
            """
            The multiplication on the right by a generator.

            INPUT:

            - ``i`` -- an element from the indexing set of the generators

            This method computes ``self * self._generators[i]``.

            EXAMPLES::

                sage: from sage.monoids.automatic_semigroup import AutomaticSemigroup
                sage: R = IntegerModRing(17)
                sage: M = AutomaticSemigroup(Family({1: R(3), 2: R(5)}), one=R.one())
                sage: M.repr_element_method("reduced_word")
                sage: M.construct()
                sage: a = M.an_element()
                sage: a.transition(1)
                [1, 1]
                sage: a.transition(2)
                [1, 2]
            """
            parent = self.parent()
            assert i in parent._generators.keys()
            return parent._retract(parent._mul(self.lift(), parent._generators_in_ambient[i]))

        def _repr_(self):
            """
            EXAMPLES::

                sage: from sage.monoids.automatic_semigroup import AutomaticSemigroup
                sage: R = IntegerModRing(19)
                sage: M = AutomaticSemigroup(Family({1: R(3), 2: R(5)}), one=R.one())
                sage: a = M.an_element(); a
                3
                sage: b = M.from_reduced_word([1,2,1]); b
                7
                sage: M.repr_element_method("reduced_word")
                sage: a
                [1]
                sage: b
                7
                sage: M.construct(up_to=b)
                sage: b
                [1, 1, 2]
            """
            if self.parent()._repr_element_method == "ambient" or self._reduced_word is None:
                return ElementWrapper._repr_(self)
            return str(self._reduced_word)

        def __copy__(self, memo=None):
            r"""
            Return ``self`` since this has unique representation.

            INPUT:

            - ``memo`` -- ignored, but required by the deepcopy API

            EXAMPLES::

                sage: from sage.monoids.automatic_semigroup import AutomaticSemigroup
                sage: R = IntegerModRing(12)
                sage: M = AutomaticSemigroup(Family({1: R(3), 2: R(5)}), one=R.one())
                sage: m = M.an_element()
                sage: copy(m) is m
                True
                sage: from copy import deepcopy
                sage: deepcopy(m) is m
                True
            """
            return self

        __deepcopy__ = __copy__


class AutomaticMonoid(AutomaticSemigroup):

    def one(self):
        """
        Return the unit of ``self``.

        EXAMPLES::

            sage: from sage.monoids.automatic_semigroup import AutomaticSemigroup
            sage: R = IntegerModRing(21)
            sage: M = R.submonoid(())
            sage: M.one()
            1
            sage: M.one().parent() is M
            True
        """
        return self._one

    # This method takes the monoid generators and adds the unit
    semigroup_generators = raw_getattr(Monoids.ParentMethods, "semigroup_generators")

    def monoid_generators(self):
        """
        Return the family of monoid generators of ``self``.

        EXAMPLES::

            sage: from sage.monoids.automatic_semigroup import AutomaticSemigroup
            sage: R = IntegerModRing(28)
            sage: M = R.submonoid(Family({1: R(3), 2: R(5)}))
            sage: M.monoid_generators()
            Finite family {1: 3, 2: 5}

        Note that the monoid generators do not include the unit,
        unlike the semigroup generators::

            sage: M.semigroup_generators()
            Family (1, 3, 5)
        """
        return self._generators
    gens = monoid_generators
