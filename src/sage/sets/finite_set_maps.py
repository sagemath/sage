r"""
Maps between finite sets

This module implements parents modeling the set of all maps between
two finite sets. At the user level, any such parent should be
constructed using the factory class :class:`FiniteSetMaps` which
properly selects which of its subclasses to use.

AUTHORS:

- Florent Hivert
"""
#*****************************************************************************
#       Copyright (C) 2010 Florent Hivert <Florent.Hivert@univ-rouen.fr>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import itertools

from sage.structure.parent import Parent
from sage.rings.integer import Integer
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.sets_cat import Sets, EmptySetError
from sage.categories.monoids import Monoids
from sage.categories.enumerated_sets import EnumeratedSets
from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
from sage.sets.integer_range import IntegerRange
from sage.sets.finite_set_map_cy import (
    FiniteSetMap_MN, FiniteSetMap_Set,
    FiniteSetEndoMap_N, FiniteSetEndoMap_Set )
from sage.misc.cachefunc import cached_method

# TODO: finite set maps should be morphisms in the category of finite sets


class FiniteSetMaps(UniqueRepresentation, Parent):
    r"""
    Maps between finite sets

    Constructs the set of all maps between two sets. The sets can be
    given using any of the three following ways:

    1. an object in the category ``Sets()``.

    2. a finite iterable. In this case, an object of the class
       :class:`~sage.sets.finite_enumerated_set.FiniteEnumeratedSet`
       is constructed from the iterable.

    3. an integer ``n`` designing the set `\{0, 1, \dots, n-1\}`. In this case
       an object of the class :class:`~sage.sets.integer_range.IntegerRange` is
       constructed.

    INPUT:

    - ``domain`` -- a set, finite iterable, or integer.

    - ``codomain`` -- a set, finite iterable, integer, or ``None``
      (default). In this last case, the maps are endo-maps of the domain.

    - ``action`` -- ``"left"`` (default) or ``"right"``. The side
      where the maps act on the domain. This is used in particular to
      define the meaning of the product (composition) of two maps.

    - ``category`` -- the category in which the sets of maps is
      constructed. By default, this is ``FiniteMonoids()`` if the domain and
      codomain coincide, and ``FiniteEnumeratedSets()`` otherwise.

    OUTPUT:

        an instance of a subclass of :class:`FiniteSetMaps` modeling
        the set of all maps between ``domain`` and ``codomain``.

    EXAMPLES:

    We construct the set ``M`` of all maps from `\{a,b\}` to `\{3,4,5\}`::

        sage: M = FiniteSetMaps(["a", "b"], [3, 4, 5]); M
        Maps from {'a', 'b'} to {3, 4, 5}
        sage: M.cardinality()
        9
        sage: M.domain()
        {'a', 'b'}
        sage: M.codomain()
        {3, 4, 5}
        sage: for f in M: print(f)
        map: a -> 3, b -> 3
        map: a -> 3, b -> 4
        map: a -> 3, b -> 5
        map: a -> 4, b -> 3
        map: a -> 4, b -> 4
        map: a -> 4, b -> 5
        map: a -> 5, b -> 3
        map: a -> 5, b -> 4
        map: a -> 5, b -> 5

    Elements can be constructed from functions and dictionaries::

        sage: M(lambda c: ord(c)-94)
        map: a -> 3, b -> 4

        sage: M.from_dict({'a':3, 'b':5})
        map: a -> 3, b -> 5

    If the domain is equal to the codomain, then maps can be
    composed::

        sage: M = FiniteSetMaps([1, 2, 3])
        sage: f = M.from_dict({1:2, 2:1, 3:3}); f
        map: 1 -> 2, 2 -> 1, 3 -> 3
        sage: g = M.from_dict({1:2, 2:3, 3:1}); g
        map: 1 -> 2, 2 -> 3, 3 -> 1

        sage: f * g
        map: 1 -> 1, 2 -> 3, 3 -> 2

    This makes `M` into a monoid::

        sage: M.category()
        Category of finite enumerated monoids
        sage: M.one()
        map: 1 -> 1, 2 -> 2, 3 -> 3

    By default, composition is from right to left, which corresponds
    to an action on the left. If one specifies ``action`` to right,
    then the composition is from left to right::

        sage: M = FiniteSetMaps([1, 2, 3], action = 'right')
        sage: f = M.from_dict({1:2, 2:1, 3:3})
        sage: g = M.from_dict({1:2, 2:3, 3:1})
        sage: f * g
        map: 1 -> 3, 2 -> 2, 3 -> 1

    If the domains and codomains are both of the form `\{0,\dots\}`,
    then one can use the shortcut::

        sage: M = FiniteSetMaps(2,3); M
        Maps from {0, 1} to {0, 1, 2}
        sage: M.cardinality()
        9

    For a compact notation, the elements are then printed as lists
    `[f(i), i=0,\dots]`::

        sage: list(M)
        [[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2], [2, 0], [2, 1], [2, 2]]

    TESTS::

        sage: TestSuite(FiniteSetMaps(0)).run()
        sage: TestSuite(FiniteSetMaps(0, 2)).run()
        sage: TestSuite(FiniteSetMaps(2, 0)).run()
        sage: TestSuite(FiniteSetMaps([], [])).run()
        sage: TestSuite(FiniteSetMaps([1, 2], [])).run()
        sage: TestSuite(FiniteSetMaps([], [1, 2])).run()
    """
    @staticmethod
    def __classcall_private__(cls, domain, codomain = None, action = "left", category = None):
        """
        TESTS::

            sage: FiniteSetMaps(3)
            Maps from {0, 1, 2} to itself
            sage: FiniteSetMaps(4, 2)
            Maps from {0, 1, 2, 3} to {0, 1}
            sage: FiniteSetMaps(4, ["a","b","c"])
            Maps from {0, 1, 2, 3} to {'a', 'b', 'c'}
            sage: FiniteSetMaps([1,2], ["a","b","c"])
            Maps from {1, 2} to {'a', 'b', 'c'}
            sage: FiniteSetMaps([1,2,4], 3)
            Maps from {1, 2, 4} to {0, 1, 2}
        """
        if codomain is None:
            if isinstance(domain, (int, Integer)):
                return FiniteSetEndoMaps_N(domain, action, category)
            else:
                if domain not in Sets():
                    domain = FiniteEnumeratedSet(domain)
                return FiniteSetEndoMaps_Set(domain, action, category)

        if isinstance(domain, (int, Integer)):
            if isinstance(codomain, (int, Integer)):
                return FiniteSetMaps_MN(domain, codomain, category)
            else:
                domain = IntegerRange(domain)
        if isinstance(codomain, (int, Integer)):
            codomain = IntegerRange(codomain)

        if domain not in Sets():
            domain = FiniteEnumeratedSet(domain)
        if codomain not in Sets():
            codomain = FiniteEnumeratedSet(codomain)
        return FiniteSetMaps_Set(domain, codomain, category)

    def cardinality(self):
        """
        The cardinality of ``self``

        EXAMPLES::

            sage: FiniteSetMaps(4, 3).cardinality()
            81
        """
        return self.codomain().cardinality()**self.domain().cardinality()


class FiniteSetMaps_MN(FiniteSetMaps):
    r"""
    The set of all maps from `\{1, 2, \dots, m\}` to `\{1, 2, \dots, n\}`.

    Users should use the factory class :class:`FiniteSetMaps` to
    create instances of this class.

    INPUT:

    - ``m``, ``n`` -- integers

    - ``category`` -- the category in which the sets of maps is
      constructed. It must be a sub-category of
      ``EnumeratedSets().Finite()`` which is the default value.
    """

    def __init__(self, m, n, category=None):
        """
        TESTS::

            sage: M = FiniteSetMaps(2,3)
            sage: M.category()
            Category of finite enumerated sets
            sage: M.__class__
            <class 'sage.sets.finite_set_maps.FiniteSetMaps_MN_with_category'>
            sage: TestSuite(M).run()
        """
        Parent.__init__(self,
                        category=EnumeratedSets().Finite().or_subcategory(category))
        self._m = Integer(m)
        self._n = Integer(n)

    def domain(self):
        """
        The domain of ``self``

        EXAMPLES::

            sage: FiniteSetMaps(3,2).domain()
            {0, 1, 2}
        """
        return IntegerRange(self._m)

    def codomain(self):
        """
        The codomain of ``self``

        EXAMPLES::

            sage: FiniteSetMaps(3,2).codomain()
            {0, 1}
        """
        return IntegerRange(self._n)

    def _repr_(self):
        """
        TESTS::

            sage: FiniteSetMaps(2,3)
            Maps from {0, 1} to {0, 1, 2}
        """
        return "Maps from %s to %s"%(self.domain(), self.codomain())

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: M = FiniteSetMaps(3,2)
            sage: [0,1,1] in M
            True
            sage: [1,2,4] in M
            False
        """
        if isinstance(x, self.element_class):
            return x.parent() is self and len(x) == self._m
        else:
            x = list(x)
            if len(x) != self._m:
                return False
            for i in x:
                if not (0 <= i < self._n):
                    return False
            return True

    def an_element(self):
        """
        Returns a map in ``self``

        EXAMPLES::

            sage: M = FiniteSetMaps(4, 2)
            sage: M.an_element()
            [0, 0, 0, 0]

            sage: M = FiniteSetMaps(0, 0)
            sage: M.an_element()
            []

        An exception :class:`~sage.categories.sets_cat.EmptySetError`
        is raised if this set is empty, that is if the codomain is
        empty and the domain is not.

            sage: M = FiniteSetMaps(4, 0)
            sage: M.cardinality()
            0
            sage: M.an_element()
            Traceback (most recent call last):
            ...
            EmptySetError
        """
        if self._m > 0 and self._n == 0:
            raise EmptySetError
        return self._from_list_([0]*self._m)

    def __iter__(self):
        """
        EXAMPLES::

            sage: M = FiniteSetMaps(2,2)
            sage: M.list()
            [[0, 0], [0, 1], [1, 0], [1, 1]]

        TESTS::

            sage: FiniteSetMaps(0,0).list()
            [[]]
            sage: FiniteSetMaps(0,1).list()
            [[]]
            sage: FiniteSetMaps(0,10).list()
            [[]]
            sage: FiniteSetMaps(1,0).list()
            []
            sage: FiniteSetMaps(1,1).list()
            [[0]]
        """
        for v in itertools.product(range(self._n), repeat=self._m):
            yield self._from_list_(v)

    def _from_list_(self, v):
        """
        EXAMPLES::

            sage: M = FiniteSetMaps(4,3)
            sage: M._from_list_([2,1,1,0])
            [2, 1, 1, 0]
        """
        return self.element_class(self, v, check=False)

    def _element_constructor_(self, *args, **keywords):
        """
        EXAMPLES::

            sage: M = FiniteSetMaps(4,3)
            sage: M([2,1,1,0])
            [2, 1, 1, 0]
        """
        return self.element_class(self, *args, **keywords)

    Element = FiniteSetMap_MN


class FiniteSetMaps_Set(FiniteSetMaps_MN):
    """
    The sets of all maps between two sets

    Users should use the factory class :class:`FiniteSetMaps` to
    create instances of this class.

    INPUT:

    - ``domain`` -- an object in the category ``FiniteSets()``.

    - ``codomain`` -- an object in the category ``FiniteSets()``.

    - ``category`` -- the category in which the sets of maps is
      constructed. It must be a sub-category of
      ``EnumeratedSets().Finite()`` which is the default value.
    """
    def __init__(self, domain, codomain, category=None):
        """
        EXAMPLES::

            sage: M = FiniteSetMaps(["a", "b"], [3, 4, 5])
            sage: M
            Maps from {'a', 'b'} to {3, 4, 5}
            sage: M.cardinality()
            9
            sage: for f in M: print(f)
            map: a -> 3, b -> 3
            map: a -> 3, b -> 4
            map: a -> 3, b -> 5
            map: a -> 4, b -> 3
            map: a -> 4, b -> 4
            map: a -> 4, b -> 5
            map: a -> 5, b -> 3
            map: a -> 5, b -> 4
            map: a -> 5, b -> 5

        TESTS::

            sage: M.__class__
            <class 'sage.sets.finite_set_maps.FiniteSetMaps_Set_with_category'>
            sage: M.category()
            Category of finite enumerated sets
            sage: TestSuite(M).run()
        """
        FiniteSetMaps_MN.__init__(self, domain.cardinality(), codomain.cardinality(),
                                 category=category)

        self._domain = domain
        self._codomain = codomain

        import sage.combinat.ranker as ranker
        ldomain = domain.list()
        lcodomain = codomain.list()
        self._unrank_domain = ranker.unrank_from_list(ldomain)
        self._rank_domain = ranker.rank_from_list(ldomain)
        self._unrank_codomain = ranker.unrank_from_list(lcodomain)
        self._rank_codomain = ranker.rank_from_list(lcodomain)

    def domain(self):
        """
        The domain of ``self``

        EXAMPLES::

            sage: FiniteSetMaps(["a", "b"], [3, 4, 5]).domain()
            {'a', 'b'}
        """
        return self._domain

    def codomain(self):
        """
        The codomain of ``self``

        EXAMPLES::

            sage: FiniteSetMaps(["a", "b"], [3, 4, 5]).codomain()
            {3, 4, 5}
        """
        return self._codomain

    # TODO: consistency from_dict / from_list
    def _from_list_(self, v):
        """
        Create a function from a list

        The list gives in the order of the element of the domain the
        rank (index) of its image in the codomain.

        EXAMPLES::

            sage: M = FiniteSetMaps(["a", "b"], [3, 4, 5])
            sage: M._from_list_([2,1])
            map: a -> 5, b -> 4
        """
        return self.element_class.from_list(self, v)

    def from_dict(self, d):
        """
        Create a map from a dictionary

        EXAMPLES::

            sage: M = FiniteSetMaps(["a", "b"], [3, 4, 5])
            sage: M.from_dict({"a": 4, "b": 3})
            map: a -> 4, b -> 3
        """
        return self.element_class.from_dict(self, d)

    Element = FiniteSetMap_Set


class FiniteSetEndoMaps_N(FiniteSetMaps_MN):
    r"""
    The sets of all maps from `\{1, 2, \dots, n\}` to itself

    Users should use the factory class :class:`FiniteSetMaps` to
    create instances of this class.

    INPUT:

    - ``n`` -- an integer.

    - ``category`` -- the category in which the sets of maps is
      constructed. It must be a sub-category of ``Monoids().Finite()``
      and ``EnumeratedSets().Finite()`` which is the default value.
    """

    def __init__(self, n, action, category=None):
        """
        EXAMPLES::

            sage: M = FiniteSetMaps(3)
            sage: M.category()
            Category of finite enumerated monoids
            sage: M.__class__
            <class 'sage.sets.finite_set_maps.FiniteSetEndoMaps_N_with_category'>
            sage: TestSuite(M).run()
        """
        category = (EnumeratedSets() & Monoids().Finite()).or_subcategory(category)
        FiniteSetMaps_MN.__init__(self, n, n, category=category)
        self._action = action

    @cached_method
    def one(self):
        """
        EXAMPLES::

            sage: M = FiniteSetMaps(4)
            sage: M.one()
            [0, 1, 2, 3]
        """
        return self._from_list_(range(self._n))

    def an_element(self):
        """
        Returns a map in ``self``

        EXAMPLES::

            sage: M = FiniteSetMaps(4)
            sage: M.an_element()
            [3, 2, 1, 0]
        """
        return self._from_list_(range(self._n-1, -1, -1))

    def _repr_(self):
        """
        TESTS::

            sage: FiniteSetMaps(2)
            Maps from {0, 1} to itself
        """
        return "Maps from %s to itself"%(self.domain())

    Element = FiniteSetEndoMap_N

class FiniteSetEndoMaps_Set(FiniteSetMaps_Set, FiniteSetEndoMaps_N):
    """
    The sets of all maps from a set to itself

    Users should use the factory class :class:`FiniteSetMaps` to
    create instances of this class.

    INPUT:

    - ``domain`` -- an object in the category ``FiniteSets()``.

    - ``category`` -- the category in which the sets of maps is
      constructed. It must be a sub-category of ``Monoids().Finite()``
      and ``EnumeratedSets().Finite()`` which is the default value.
     """
    def __init__(self, domain, action, category=None):
        """
        TESTS::

            sage: M = FiniteSetMaps(["a", "b", "c"])
            sage: M.category()
            Category of finite enumerated monoids
            sage: M.__class__
            <class 'sage.sets.finite_set_maps.FiniteSetEndoMaps_Set_with_category'>
            sage: TestSuite(M).run()
        """
        category = (EnumeratedSets() & Monoids().Finite()).or_subcategory(category)
        FiniteSetMaps_MN.__init__(self, domain.cardinality(), domain.cardinality(),
                                 category=category)

        self._domain = domain
        self._codomain = domain

        import sage.combinat.ranker as ranker
        ldomain = domain.list()
        self._unrank_domain = ranker.unrank_from_list(ldomain)
        self._rank_domain = ranker.rank_from_list(ldomain)
        self._unrank_codomain = self._unrank_domain
        self._rank_codomain = self._rank_domain
        self._action = action

    Element = FiniteSetEndoMap_Set

