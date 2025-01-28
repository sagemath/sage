r"""
Families

A Family is an associative container which models a family
`(f_i)_{i \in I}`. Then, ``f[i]`` returns the element of the family indexed by
``i``. Whenever available, set and combinatorial class operations (counting,
iteration, listing) on the family are induced from those of the index
set. Families should be created through the :func:`Family` function.

AUTHORS:

- Nicolas Thiery (2008-02): initial release

- Florent Hivert (2008-04): various fixes, cleanups and improvements.

TESTS:

Check :issue:`12482` (shall be run in a fresh session)::

    sage: P = Partitions(3)                                                             # needs sage.combinat
    sage: Family(P, lambda x: x).category()                                             # needs sage.combinat
    Category of finite enumerated sets
"""

# *****************************************************************************
#       Copyright (C) 2008-2017 Nicolas Thiery <nthiery at users.sf.net>
#                     2008-2009 Mike Hansen <mhansen@gmail.com>
#                     2008-2010 Florent Hivert <Florent.Hivert@univ-rouen.fr>
#                     2013-2021 Travis Scrimshaw
#                     2014      Nathann Cohen
#                     2017      Erik M. Bray
#                     2018      Frédéric Chapoton
#                     2019      Markus Wageringel
#                     2022-2023 Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import types
from copy import copy
from pprint import pformat, saferepr
from collections.abc import Iterable

from sage.categories.enumerated_sets import EnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.misc.cachefunc import cached_method
from sage.misc.call import AttrCallObject
from sage.rings.infinity import Infinity
from sage.rings.integer import Integer
from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
from sage.sets.non_negative_integers import NonNegativeIntegers


def Family(indices, function=None, hidden_keys=[], hidden_function=None, lazy=False, name=None):
    r"""
    A Family is an associative container which models a family
    `(f_i)_{i \in I}`. Then, ``f[i]`` returns the element of the family
    indexed by `i`. Whenever available, set and combinatorial class
    operations (counting, iteration, listing) on the family are induced
    from those of the index set.

    There are several available implementations (classes) for different
    usages; Family serves as a factory, and will create instances of
    the appropriate classes depending on its arguments.

    INPUT:

    - ``indices`` -- the indices for the family
    - ``function`` -- (optional) the function `f` applied to all visible
      indices; the default is the identity function
    - ``hidden_keys`` -- (optional) a list of hidden indices that can be
      accessed through ``my_family[i]``
    - ``hidden_function`` -- (optional) a function for the hidden indices
    - ``lazy`` -- boolean (default: ``False``); whether the family is lazily
      created or not; if the indices are infinite, then this is automatically
      made ``True``
    - ``name`` -- (optional) the name of the function; only used when the
      family is lazily created via a function

    EXAMPLES:

    In its simplest form, a list `l = [l_0, l_1, \ldots, l_{\ell}]` or a
    tuple by itself is considered as the family `(l_i)_{i \in I}` where
    `I` is the set `\{0, \ldots, \ell\}` where `\ell` is ``len(l) - 1``.
    So ``Family(l)`` returns the corresponding family::

        sage: f = Family([1,2,3])
        sage: f
        Family (1, 2, 3)
        sage: f = Family((1,2,3))
        sage: f
        Family (1, 2, 3)

    Instead of a list you can as well pass any iterable object::

        sage: f = Family(2*i+1 for i in [1,2,3])
        sage: f
        Family (3, 5, 7)

    A family can also be constructed from a dictionary ``t``. The resulting
    family is very close to ``t``, except that the elements of the family
    are the values of ``t``. Here, we define the family
    `(f_i)_{i \in \{3,4,7\}}` with `f_3 = a`, `f_4 = b`, and `f_7 = d`::

        sage: f = Family({3: 'a', 4: 'b', 7: 'd'})
        sage: f
        Finite family {3: 'a', 4: 'b', 7: 'd'}
        sage: f[7]
        'd'
        sage: len(f)
        3
        sage: list(f)
        ['a', 'b', 'd']
        sage: [ x for x in f ]
        ['a', 'b', 'd']
        sage: f.keys()
        [3, 4, 7]
        sage: 'b' in f
        True
        sage: 'e' in f
        False

    A family can also be constructed by its index set `I` and
    a function `f`, as in `(f(i))_{i \in I}`::

        sage: f = Family([3,4,7], lambda i: 2*i)
        sage: f
        Finite family {3: 6, 4: 8, 7: 14}
        sage: f.keys()
        [3, 4, 7]
        sage: f[7]
        14
        sage: list(f)
        [6, 8, 14]
        sage: [x for x in f]
        [6, 8, 14]
        sage: len(f)
        3

    By default, all images are computed right away, and stored in an internal
    dictionary::

        sage: f = Family((3,4,7), lambda i: 2*i)
        sage: f
        Finite family {3: 6, 4: 8, 7: 14}

    Note that this requires all the elements of the list to be
    hashable. One can ask instead for the images `f(i)` to be computed
    lazily, when needed::

        sage: f = Family([3,4,7], lambda i: 2*i, lazy=True)
        sage: f
        Lazy family (<lambda>(i))_{i in [3, 4, 7]}
        sage: f[7]
        14
        sage: list(f)
        [6, 8, 14]
        sage: [x for x in f]
        [6, 8, 14]

    This allows in particular for modeling infinite families::

        sage: f = Family(ZZ, lambda i: 2*i, lazy=True)
        sage: f
        Lazy family (<lambda>(i))_{i in Integer Ring}
        sage: f.keys()
        Integer Ring
        sage: f[1]
        2
        sage: f[-5]
        -10
        sage: i = iter(f)
        sage: next(i), next(i), next(i), next(i), next(i)
        (0, 2, -2, 4, -4)

    Note that the ``lazy`` keyword parameter is only needed to force
    laziness. Usually it is automatically set to a correct default value (ie:
    ``False`` for finite data structures and ``True`` for enumerated sets::

        sage: f == Family(ZZ, lambda i: 2*i)
        True

    Beware that for those kind of families len(f) is not supposed to
    work. As a replacement, use the .cardinality() method::

       sage: f = Family(Permutations(3), attrcall("to_lehmer_code"))
       sage: list(f)
       [[0, 0, 0], [0, 1, 0], [1, 0, 0], [1, 1, 0], [2, 0, 0], [2, 1, 0]]
       sage: f.cardinality()
       6

    Caveat: Only certain families with lazy behavior can be pickled. In
    particular, only functions that work with Sage's pickle_function
    and unpickle_function (in sage.misc.fpickle) will correctly
    unpickle. The following two work::

       sage: f = Family(Permutations(3), lambda p: p.to_lehmer_code()); f
       Lazy family (<lambda>(i))_{i in Standard permutations of 3}
       sage: f == loads(dumps(f))
       True

       sage: f = Family(Permutations(3), attrcall("to_lehmer_code")); f
       Lazy family (i.to_lehmer_code())_{i in Standard permutations of 3}
       sage: f == loads(dumps(f))
       True

    But this one does not::

       sage: def plus_n(n): return lambda x: x+n
       sage: f = Family([1,2,3], plus_n(3), lazy=True); f
       Lazy family (<lambda>(i))_{i in [1, 2, 3]}
       sage: f == loads(dumps(f))
       Traceback (most recent call last):
       ...
       ValueError: Cannot pickle code objects from closures

    Finally, it can occasionally be useful to add some hidden elements
    in a family, which are accessible as ``f[i]``, but do not appear in the
    keys or the container operations::

        sage: f = Family([3,4,7], lambda i: 2*i, hidden_keys=[2])
        sage: f
        Finite family {3: 6, 4: 8, 7: 14}
        sage: f.keys()
        [3, 4, 7]
        sage: f.hidden_keys()
        [2]
        sage: f[7]
        14
        sage: f[2]
        4
        sage: list(f)
        [6, 8, 14]
        sage: [x for x in f]
        [6, 8, 14]
        sage: len(f)
        3

    The following example illustrates when the function is actually
    called::

        sage: def compute_value(i):
        ....:     print('computing 2*'+str(i))
        ....:     return 2*i
        sage: f = Family([3,4,7], compute_value, hidden_keys=[2])
        computing 2*3
        computing 2*4
        computing 2*7
        sage: f
        Finite family {3: 6, 4: 8, 7: 14}
        sage: f.keys()
        [3, 4, 7]
        sage: f.hidden_keys()
        [2]
        sage: f[7]
        14
        sage: f[2]
        computing 2*2
        4
        sage: f[2]
        4
        sage: list(f)
        [6, 8, 14]
        sage: [x for x in f]
        [6, 8, 14]
        sage: len(f)
        3

    Here is a close variant where the function for the hidden keys is
    different from that for the other keys::

        sage: f = Family([3,4,7], lambda i: 2*i, hidden_keys=[2], hidden_function = lambda i: 3*i)
        sage: f
        Finite family {3: 6, 4: 8, 7: 14}
        sage: f.keys()
        [3, 4, 7]
        sage: f.hidden_keys()
        [2]
        sage: f[7]
        14
        sage: f[2]
        6
        sage: list(f)
        [6, 8, 14]
        sage: [x for x in f]
        [6, 8, 14]
        sage: len(f)
        3

    Family accept finite and infinite EnumeratedSets as input::

        sage: f = Family(FiniteEnumeratedSet([1,2,3]))
        sage: f
        Family (1, 2, 3)
        sage: f = Family(NonNegativeIntegers())
        sage: f
        Family (Non negative integers)

    ::

        sage: f = Family(FiniteEnumeratedSet([3,4,7]), lambda i: 2*i)
        sage: f
        Finite family {3: 6, 4: 8, 7: 14}
        sage: f.keys()
        {3, 4, 7}
        sage: f[7]
        14
        sage: list(f)
        [6, 8, 14]
        sage: [x for x in f]
        [6, 8, 14]
        sage: len(f)
        3

    TESTS::

        sage: f = Family({1:'a', 2:'b', 3:'c'})
        sage: f
        Finite family {1: 'a', 2: 'b', 3: 'c'}
        sage: f[2]
        'b'
        sage: loads(dumps(f)) == f
        True

    ::

        sage: f = Family({1:'a', 2:'b', 3:'c'}, lazy=True)
        Traceback (most recent call last):
        ...
        ValueError: lazy keyword only makes sense together with function keyword

    ::

        sage: f = Family(list(range(1,27)), lambda i: chr(i+96))
        sage: f
        Finite family {1: 'a', 2: 'b', 3: 'c', 4: 'd', 5: 'e', 6: 'f', 7: 'g',
        8: 'h', 9: 'i', 10: 'j', 11: 'k', 12: 'l', 13: 'm', 14: 'n', 15: 'o',
        16: 'p', 17: 'q', 18: 'r', 19: 's', 20: 't', 21: 'u', 22: 'v', 23: 'w',
        24: 'x', 25: 'y', 26: 'z'}
        sage: f[2]
        'b'

    The factory ``Family`` is supposed to be idempotent. We test this feature here::

        sage: from sage.sets.family import FiniteFamily, LazyFamily, TrivialFamily
        sage: f = FiniteFamily({3: 'a', 4: 'b', 7: 'd'})
        sage: g = Family(f)
        sage: f == g
        True

        sage: f = Family([3,4,7], lambda i: 2*i, hidden_keys=[2])
        sage: g = Family(f)
        sage: f == g
        True

        sage: f = LazyFamily([3,4,7], lambda i: 2*i)
        sage: g = Family(f)
        sage: f == g
        True

        sage: f = TrivialFamily([3,4,7])
        sage: g = Family(f)
        sage: f == g
        True

    A family should keep the order of the keys::

        sage: f = Family(["c", "a", "b"], lambda i: 2*i)
        sage: list(f)
        ['cc', 'aa', 'bb']

    Even with hidden keys (see :issue:`22955`)::

        sage: f = Family(["c", "a", "b"], lambda i: 2*i,
        ....:           hidden_keys=[5], hidden_function=lambda i: i%2)
        sage: list(f)
        ['cc', 'aa', 'bb']

    Only the hidden function is applied to the hidden keys::

        sage: f[5]
        1
    """
    assert isinstance(hidden_keys, list)
    assert isinstance(lazy, bool)

    if not hidden_keys:
        if hidden_function is not None:
            raise ValueError("hidden_function keyword only makes sense "
                             "together with hidden_keys keyword !")
        if function is None:
            if lazy:
                raise ValueError("lazy keyword only makes sense together with function keyword")
            if isinstance(indices, dict):
                return FiniteFamily(indices)
            if isinstance(indices, (list, tuple)):
                return TrivialFamily(indices)
            if isinstance(indices, (FiniteFamily, LazyFamily, TrivialFamily)):
                return indices
            if indices in EnumeratedSets():
                return EnumeratedFamily(indices)
            if isinstance(indices, Iterable):
                return TrivialFamily(indices)

            raise NotImplementedError
        if (isinstance(indices, (list, tuple, FiniteEnumeratedSet))
                and not lazy):
            return FiniteFamily({i: function(i) for i in indices},
                                keys=indices)

        return LazyFamily(indices, function, name)
    if lazy:
        raise ValueError("lazy keyword is incompatible with hidden keys")
    if hidden_function is None:
        hidden_function = function
    return FiniteFamilyWithHiddenKeys({i: function(i) for i in indices},
                                      hidden_keys, hidden_function,
                                      keys=indices)


cdef class AbstractFamily(Parent):
    """
    The abstract class for family.

    Any family belongs to a class which inherits from :class:`AbstractFamily`.
    """
    def hidden_keys(self):
        """
        Return the hidden keys of the family, if any.

        EXAMPLES::

            sage: f = Family({3: 'a', 4: 'b', 7: 'd'})
            sage: f.hidden_keys()
            []
        """
        return []

    def keys(self):
        """
        Return the keys of the family.

        EXAMPLES::

            sage: f = Family({3: 'a', 4: 'b', 7: 'd'})
            sage: sorted(f.keys())
            [3, 4, 7]
        """
        raise NotImplementedError

    def values(self):
        """
        Return the elements (values) of this family.

        EXAMPLES::

            sage: f = Family(["c", "a", "b"], lambda x: x + x)
            sage: sorted(f.values())
            ['aa', 'bb', 'cc']
        """
        raise NotImplementedError

    def items(self):
        """
        Return an iterator for key-value pairs.

        A key can only appear once, but if the function is not injective, values may
        appear multiple times.

        EXAMPLES::

            sage: f = Family([-2, -1, 0, 1, 2], abs)
            sage: list(f.items())
            [(-2, 2), (-1, 1), (0, 0), (1, 1), (2, 2)]
        """
        return zip(self.keys(), self.values())

    def zip(self, f, other, name=None):
        r"""
        Given two families with same index set `I` (and same hidden
        keys if relevant), returns the family
        `( f(self[i], other[i]) )_{i \in I}`

        .. TODO:: generalize to any number of families and merge with map?

        EXAMPLES::

            sage: f = Family({3: 'a', 4: 'b', 7: 'd'})
            sage: g = Family({3: '1', 4: '2', 7: '3'})
            sage: h = f.zip(lambda x,y: x+y, g)
            sage: list(h)
            ['a1', 'b2', 'd3']
        """
        assert self.keys() == other.keys()
        assert self.hidden_keys() == other.hidden_keys()
        return Family(self.keys(), lambda i: f(self[i], other[i]),
                      hidden_keys=self.hidden_keys(), name=name)

    def map(self, f, name=None):
        r"""
        Return the family `( f(\mathtt{self}[i]) )_{i \in I}`, where
        `I` is the index set of ``self``.

        EXAMPLES::

            sage: f = Family({3: 'a', 4: 'b', 7: 'd'})
            sage: g = f.map(lambda x: x+'1')
            sage: list(g)
            ['a1', 'b1', 'd1']
        """
        return Family(self.keys(), lambda i: f(self[i]), hidden_keys=self.hidden_keys(), name=name)

    # temporary; tested by TestSuite.
    _an_element_ = EnumeratedSets.ParentMethods._an_element_

    @cached_method
    def inverse_family(self):
        """
        Return the inverse family, with keys and values exchanged. This
        presumes that there are no duplicate values in ``self``.

        This default implementation is not lazy and therefore will
        only work with not too big finite families. It is also cached
        for the same reason::

            sage: Family({3: 'a', 4: 'b', 7: 'd'}).inverse_family()
            Finite family {'a': 3, 'b': 4, 'd': 7}

            sage: Family((3,4,7)).inverse_family()
            Finite family {3: 0, 4: 1, 7: 2}
        """
        return Family({self[k]: k for k in self.keys()})


cdef class FiniteFamily(AbstractFamily):
    r"""
    A :class:`FiniteFamily` is an associative container which models a finite
    family `(f_i)_{i \in I}`. Its elements `f_i` are therefore its
    values. Instances should be created via the :func:`Family` factory. See its
    documentation for examples and tests.

    EXAMPLES:

    We define the family `(f_i)_{i \in \{3,4,7\}}` with `f_3=a`,
    `f_4=b`, and `f_7=d`::

        sage: from sage.sets.family import FiniteFamily
        sage: f = FiniteFamily({3: 'a', 4: 'b', 7: 'd'})

    Individual elements are accessible as in a usual dictionary::

        sage: f[7]
        'd'

    And the other usual dictionary operations are also available::

        sage: len(f)
        3
        sage: f.keys()
        [3, 4, 7]

    However f behaves as a container for the `f_i`'s::

        sage: list(f)
        ['a', 'b', 'd']
        sage: [ x for x in f ]
        ['a', 'b', 'd']

    The order of the elements can be specified using the ``keys`` optional argument::

        sage: f = FiniteFamily({"a": "aa", "b": "bb", "c" : "cc" }, keys = ["c", "a", "b"])
        sage: list(f)
        ['cc', 'aa', 'bb']
    """

    def __init__(self, dictionary, keys=None):
        """
        TESTS::

            sage: from sage.sets.family import FiniteFamily
            sage: f = FiniteFamily({3: 'a', 4: 'b', 7: 'd'})
            sage: TestSuite(f).run()

        Check for bug :issue:`5538`::

            sage: d = {1:"a", 3:"b", 4:"c"}
            sage: f = Family(d)
            sage: d[2] = 'DD'
            sage: f
            Finite family {1: 'a', 3: 'b', 4: 'c'}
            """
        # TODO: use keys to specify the order of the elements
        Parent.__init__(self, category=FiniteEnumeratedSets())
        self._dictionary = dict(dictionary)
        self._keys = keys

    @cached_method
    def __hash__(self):
        """
        Return a hash value for ``self``.

        EXAMPLES::

            sage: f = Family(["c", "a", "b"], lambda x: x+x)
            sage: hash(f) == hash(f)
            True
            sage: f2 = Family(["a", "c", "b"], lambda x: x+x)
            sage: hash(f) == hash(f2)
            True
            sage: g = Family(["b", "c", "a"], lambda x: x+x+x)
            sage: hash(f) == hash(g)
            False

        ::

            sage: f = Family({1:[1,2]})
            sage: hash(f) == hash(f)
            True
        """
        try:
            return hash(frozenset(self._dictionary.items()))
        except (TypeError, ValueError):
            return hash(frozenset(self.keys() +
                                  [repr(v) for v in self.values()]))

    def __bool__(self):
        r"""
        Return if ``self`` is empty or not.

        EXAMPLES::

            sage: from sage.sets.family import TrivialFamily
            sage: f = Family(["c", "a", "b"], lambda x: x+x)
            sage: bool(f)
            True
            sage: g = Family({})
            sage: bool(g)
            False
            sage: h = Family([], lambda x: x+x)
            sage: bool(h)
            False
        """
        return bool(self._dictionary)

    def keys(self):
        """
        Return the index set of this family.

        EXAMPLES::

            sage: f = Family(["c", "a", "b"], lambda x: x+x)
            sage: f.keys()
            ['c', 'a', 'b']
        """
        return (self._keys if self._keys is not None
                else list(self._dictionary))

    def values(self):
        """
        Return the elements of this family.

        EXAMPLES::

            sage: f = Family(["c", "a", "b"], lambda x: x+x)
            sage: f.values()
            ['cc', 'aa', 'bb']
        """
        if self._keys is not None:
            return [self._dictionary[key] for key in self._keys]
        else:
            return list(self._dictionary.values())

    def has_key(self, k):
        """
        Return whether ``k`` is a key of ``self``.

        EXAMPLES::

            sage: Family({"a":1, "b":2, "c":3}).has_key("a")
            True
            sage: Family({"a":1, "b":2, "c":3}).has_key("d")
            False
        """
        return k in self._dictionary

    def __eq__(self, other):
        """
        EXAMPLES::

            sage: f = Family({1:'a', 2:'b', 3:'c'})
            sage: g = Family({1:'a', 2:'b', 3:'c'})
            sage: f == g
            True

        TESTS::

            sage: from sage.sets.family import FiniteFamily

            sage: f1 = FiniteFamily({1:'a', 2:'b', 3:'c'}, keys = [1,2,3])
            sage: g1 = FiniteFamily({1:'a', 2:'b', 3:'c'}, keys = [1,2,3])
            sage: h1 = FiniteFamily({1:'a', 2:'b', 3:'c'}, keys = [2,1,3])

            sage: f1 == g1
            True
            sage: f1 == h1
            False
            sage: f1 == f
            False
        """
        return (isinstance(other, self.__class__) and
                self._keys == (<FiniteFamily> other)._keys and
                self._dictionary == (<FiniteFamily> other)._dictionary)

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.sets.family import FiniteFamily
            sage: FiniteFamily({3: 'a'}) # indirect doctest
            Finite family {3: 'a'}

            sage: FiniteFamily({3: 'a', 4: 'b'}) # indirect doctest
            Finite family {3: 'a', 4: 'b'}

            sage: FiniteFamily({3: 'a', 4: 'b'}, keys=[4,3]) # indirect doctest
            Finite family {4: 'b', 3: 'a'}
        """
        if self._keys is None:
            d = ' '.join(pformat(self._dictionary)[1:-1].splitlines())
        else:
            d = ', '.join('{}: {}'.format(saferepr(key),
                                          saferepr(self._dictionary[key]))
                          for key in self._keys)
        return 'Finite family {{{}}}'.format(d)

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: from sage.sets.family import FiniteFamily
            sage: f = FiniteFamily({3: 'a'})
            sage: 'a' in f
            True
            sage: 'b' in f
            False
        """
        return x in self.values()

    def __len__(self):
        """
        Return the number of elements in ``self``.

        EXAMPLES::

            sage: from sage.sets.family import FiniteFamily
            sage: f = FiniteFamily({3: 'a', 4: 'b', 7: 'd'})
            sage: len(f)
            3
        """
        return len(self._dictionary)

    def cardinality(self):
        """
        Return the number of elements in ``self``.

        EXAMPLES::

            sage: from sage.sets.family import FiniteFamily
            sage: f = FiniteFamily({3: 'a', 4: 'b', 7: 'd'})
            sage: f.cardinality()
            3
        """
        return Integer(len(self._dictionary))

    def __iter__(self):
        """
        EXAMPLES::

            sage: from sage.sets.family import FiniteFamily
            sage: f = FiniteFamily({3: 'a'})
            sage: i = iter(f)
            sage: next(i)
            'a'
        """
        return iter(self.values())

    def __getitem__(self, i):
        """
        Note that we can't just do self.__getitem__ =
        dictionary.__getitem__ in the __init__ method since Python
        queries the object's type/class for the special methods rather than
        querying the object itself.

        EXAMPLES::

            sage: from sage.sets.family import FiniteFamily
            sage: f = FiniteFamily({3: 'a', 4: 'b', 7: 'd'})
            sage: f[3]
            'a'
        """
        return self._dictionary[i]

    # For the pickle and copy modules
    def __getstate__(self):
        """
        TESTS::

            sage: from sage.sets.family import FiniteFamily
            sage: f = FiniteFamily({3: 'a'})
            sage: f.__getstate__()
            {'dictionary': {3: 'a'}, 'keys': None}
        """
        return {'dictionary': self._dictionary, 'keys': self._keys}

    def __setstate__(self, state):
        """
        TESTS::

            sage: from sage.sets.family import FiniteFamily
            sage: f = FiniteFamily({3: 'a'})
            sage: f.__setstate__({'dictionary': {4:'b'}})
            sage: f
            Finite family {4: 'b'}
        """
        self.__init__(state['dictionary'], keys=state.get("keys"))


class FiniteFamilyWithHiddenKeys(FiniteFamily):
    r"""
    A close variant of :class:`FiniteFamily` where the family contains some
    hidden keys whose corresponding values are computed lazily (and
    remembered). Instances should be created via the :func:`Family` factory.
    See its documentation for examples and tests.

    Caveat: Only instances of this class whose functions are compatible
    with :mod:`sage.misc.fpickle` can be pickled.
    """
    def __init__(self, dictionary, hidden_keys, hidden_function, keys=None):
        """
        EXAMPLES::

            sage: f = Family([3,4,7], lambda i: 2*i, hidden_keys=[2])
            sage: TestSuite(f).run()
        """
        FiniteFamily.__init__(self, dictionary, keys=keys)
        self._hidden_keys = hidden_keys
        self.hidden_function = hidden_function
        self.hidden_dictionary = {}

        # would be better to define as usual method
        # any better to unset the def of __getitem__ by FiniteFamily?
        # self.__getitem__ = lambda i: dictionary[i] if dictionary.has_key(i) else hidden_dictionary[i]

    def __getitem__(self, i):
        """
        EXAMPLES::

            sage: f = Family([3,4,7], lambda i: 2*i, hidden_keys=[2])
            sage: f[3]
            6
            sage: f[2]
            4
            sage: f[5]
            Traceback (most recent call last):
            ...
            KeyError
        """
        try:
            return FiniteFamily.__getitem__(self, i)
        except KeyError:
            try:
                return self.hidden_dictionary[i]
            except KeyError:
                if i not in self._hidden_keys:
                    raise KeyError
                v = self.hidden_function(i)
                self.hidden_dictionary[i] = v
                return v

    def hidden_keys(self):
        """
        Return ``self``'s hidden keys.

        EXAMPLES::

            sage: f = Family([3,4,7], lambda i: 2*i, hidden_keys=[2])
            sage: f.hidden_keys()
            [2]
        """
        return self._hidden_keys

    def __getstate__(self):
        """
        TESTS::

            sage: f = Family([3,4,7], lambda i: 2*i, hidden_keys=[2])
            sage: d = f.__getstate__()
            sage: d['hidden_keys']
            [2]
        """
        from sage.misc.fpickle import pickle_function
        f = pickle_function(self.hidden_function)
        state = super().__getstate__()
        state.update({'hidden_keys': self._hidden_keys,
                      'hidden_dictionary': self.hidden_dictionary,
                      'hidden_function': f})
        return state

    def __setstate__(self, d):
        """
        TESTS::

            sage: f = Family([3,4,7], lambda i: 2*i, hidden_keys=[2])
            sage: d = f.__getstate__()
            sage: f = Family([4,5,6], lambda i: 2*i, hidden_keys=[2])
            sage: f.__setstate__(d)
            sage: f.keys()
            [3, 4, 7]
            sage: f[3]
            6
        """
        hidden_function = d['hidden_function']
        if isinstance(hidden_function, (str, bytes)):
            # Let's assume that hidden_function is an unpickled function.
            from sage.misc.fpickle import unpickle_function
            hidden_function = unpickle_function(hidden_function)
        self.__init__(d['dictionary'], d['hidden_keys'], hidden_function)
        self.hidden_dictionary = d['hidden_dictionary']
        # Old pickles from before Issue #22955 may not have a 'keys'
        if 'keys' in d:
            self._keys = d['keys']
        else:
            self._keys = None


class LazyFamily(AbstractFamily):
    r"""
    A LazyFamily(I, f) is an associative container which models the
    (possibly infinite) family `(f(i))_{i \in I}`.

    Instances should be created via the :func:`Family` factory. See its
    documentation for examples and tests.
    """
    def __init__(self, set, function, name=None):
        """
        TESTS::

            sage: from sage.sets.family import LazyFamily
            sage: f = LazyFamily([3,4,7], lambda i: 2*i); f
            Lazy family (<lambda>(i))_{i in [3, 4, 7]}
            sage: TestSuite(f).run()

        Check for :issue:`5538`::

            sage: l = [3,4,7]
            sage: f = LazyFamily(l, lambda i: 2*i)
            sage: l[1] = 18
            sage: f
            Lazy family (<lambda>(i))_{i in [3, 4, 7]}
        """
        if set in FiniteEnumeratedSets():
            category = FiniteEnumeratedSets()
        elif set in InfiniteEnumeratedSets():
            category = InfiniteEnumeratedSets()
        elif isinstance(set, (list, tuple, range)):
            category = FiniteEnumeratedSets()
        else:
            category = EnumeratedSets()

        Parent.__init__(self, category=category)

        self.set = copy(set)
        self.function = function
        self.function_name = name

    def __bool__(self):
        r"""
        Return if ``self`` is empty or not.

        EXAMPLES::

            sage: from sage.sets.family import LazyFamily
            sage: f = LazyFamily([3,4,7], lambda i: 2*i)
            sage: bool(f)
            True
            sage: g = LazyFamily([], lambda i: 2*i)
            sage: bool(g)
            False
            sage: h = Family(ZZ, lambda x: x+x)
            sage: bool(h)
            True
        """
        return bool(self.set)

    @cached_method
    def __hash__(self):
        """
        Return a hash value for ``self``.

        EXAMPLES::

            sage: from sage.sets.family import LazyFamily
            sage: f = LazyFamily([3,4,7], lambda i: 2*i)
            sage: hash(f) == hash(f)
            True
            sage: g = LazyFamily(ZZ, lambda i: 2*i)
            sage: hash(g) == hash(g)
            True
            sage: h = LazyFamily(ZZ, lambda i: 2*i, name='foo')
            sage: hash(h) == hash(h)
            True

        ::

            sage: class X():
            ....:     def __call__(self, x):
            ....:         return x
            ....:     __hash__ = None
            sage: f = Family([1,2,3], X())
            sage: hash(f) == hash(f)
            True
        """
        try:
            return hash(self.keys()) + hash(self.function)
        except (TypeError, ValueError):
            return super().__hash__()

    def __eq__(self, other):
        """
        WARNING: Since there is no way to compare function, we only compare
        their name.

        TESTS::

            sage: from sage.sets.family import LazyFamily
            sage: fun = lambda i: 2*i
            sage: f = LazyFamily([3,4,7], fun)
            sage: g = LazyFamily([3,4,7], fun)
            sage: f == g
            True
        """
        if not isinstance(other, self.__class__):
            return False
        if not self.set == other.set:
            return False
        return repr(self) == repr(other)

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.sets.family import LazyFamily
            sage: def fun(i): 2*i
            sage: f = LazyFamily([3,4,7], fun); f
            Lazy family (fun(i))_{i in [3, 4, 7]}

            sage: f = Family(Permutations(3), attrcall("to_lehmer_code"), lazy=True); f
            Lazy family (i.to_lehmer_code())_{i in Standard permutations of 3}

            sage: f = LazyFamily([3,4,7], lambda i: 2*i); f
            Lazy family (<lambda>(i))_{i in [3, 4, 7]}

            sage: f = LazyFamily([3,4,7], lambda i: 2*i, name='foo'); f
            Lazy family (foo(i))_{i in [3, 4, 7]}

        TESTS:

            Check that using a class as the function is correctly handled::

                sage: Family(NonNegativeIntegers(), PerfectMatchings)                   # needs sage.combinat
                Lazy family (<class 'sage.combinat.perfect_matching.PerfectMatchings'>(i))_{i in Non negative integers}
        """
        if self.function_name is not None:
            name = self.function_name + "(i)"
        elif isinstance(self.function, types.LambdaType):
            name = self.function.__name__
            name = name + "(i)"
        else:
            name = repr(self.function)
            if isinstance(self.function, AttrCallObject):
                name = "i" + name[1:]
            else:
                name = name + "(i)"
        return "Lazy family ({})_{{i in {}}}".format(name, self.set)

    def keys(self):
        """
        Return ``self``'s keys.

        EXAMPLES::

            sage: from sage.sets.family import LazyFamily
            sage: f = LazyFamily([3,4,7], lambda i: 2*i)
            sage: f.keys()
            [3, 4, 7]
        """
        return self.set

    def cardinality(self):
        """
        Return the number of elements in ``self``.

        EXAMPLES::

            sage: from sage.sets.family import LazyFamily
            sage: f = LazyFamily([3,4,7], lambda i: 2*i)
            sage: f.cardinality()
            3
            sage: l = LazyFamily(NonNegativeIntegers(), lambda i: 2*i)
            sage: l.cardinality()
            +Infinity

        TESTS:

        Check that :issue:`15195` is fixed::

            sage: C = cartesian_product([PositiveIntegers(), [1,2,3]])
            sage: C.cardinality()
            +Infinity
            sage: F = Family(C, lambda x: x)
            sage: F.cardinality()
            +Infinity
            """
        try:
            return Integer(len(self.set))
        except (AttributeError, NotImplementedError, TypeError):
            return self.set.cardinality()

    def __iter__(self):
        """
        EXAMPLES::

            sage: from sage.sets.family import LazyFamily
            sage: f = LazyFamily([3,4,7], lambda i: 2*i)
            sage: [i for i in f]
            [6, 8, 14]
        """
        for i in self.set:
            yield self[i]

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: from sage.sets.family import LazyFamily
            sage: f = LazyFamily([3,4,7], lambda i: 2*i)
            sage: 3 in f, 14 in f
            (False, True)

        By default this expands the lazy family, which is only done for
        families known to be finite::

            sage: 5 in LazyFamily(NonNegativeIntegers(), lambda i: 2*i)
            Traceback (most recent call last):
            ...
            ValueError: family must be finite to check containment
        """
        if self not in FiniteEnumeratedSets():
            raise ValueError('family must be finite to check containment')
        return x in iter(self)

    def __getitem__(self, i):
        """
        EXAMPLES::

            sage: from sage.sets.family import LazyFamily
            sage: f = LazyFamily([3,4,7], lambda i: 2*i)
            sage: f[3]
            6

        TESTS::

            sage: f[5]
            10
        """
        return self.function(i)

    def __getstate__(self):
        """
        EXAMPLES::

            sage: from sage.sets.family import LazyFamily
            sage: f = LazyFamily([3,4,7], lambda i: 2*i)
            sage: d = f.__getstate__()
            sage: d['set']
            [3, 4, 7]

            sage: f = LazyFamily(Permutations(3), lambda p: p.to_lehmer_code())
            sage: f == loads(dumps(f))
            True

            sage: f = LazyFamily(Permutations(3), attrcall("to_lehmer_code"))
            sage: f == loads(dumps(f))
            True
        """
        f = self.function
        # This should be done once for all by registering
        # sage.misc.fpickle.pickle_function to copyreg
        if isinstance(f, types.FunctionType):
            from sage.misc.fpickle import pickle_function
            f = pickle_function(f)

        return {'set': self.set,
                'function': f}

    def __setstate__(self, d):
        """
        EXAMPLES::

            sage: from sage.sets.family import LazyFamily
            sage: f = LazyFamily([3,4,7], lambda i: 2*i)
            sage: d = f.__getstate__()
            sage: f = LazyFamily([4,5,6], lambda i: 2*i)
            sage: f.__setstate__(d)
            sage: f.keys()
            [3, 4, 7]
            sage: f[3]
            6
        """
        function = d['function']
        if isinstance(function, bytes):
            # Let's assume that function is an unpickled function.
            from sage.misc.fpickle import unpickle_function
            function = unpickle_function(function)

        self.__init__(d['set'], function)


class TrivialFamily(AbstractFamily):
    r"""
    :class:`TrivialFamily` turns a list/tuple `c` into a family indexed by the
    set `\{0, \dots, |c|-1\}`.

    Instances should be created via the :func:`Family` factory. See its
    documentation for examples and tests.
    """
    def __init__(self, enumeration):
        """
        EXAMPLES::

            sage: from sage.sets.family import TrivialFamily
            sage: f = TrivialFamily((3,4,7)); f
            Family (3, 4, 7)
            sage: f = TrivialFamily([3,4,7]); f
            Family (3, 4, 7)
            sage: TestSuite(f).run()
        """
        Parent.__init__(self, category=FiniteEnumeratedSets())
        self._enumeration = tuple(enumeration)

    def __bool__(self):
        r"""
        Return if ``self`` is empty or not.

        EXAMPLES::

            sage: from sage.sets.family import TrivialFamily
            sage: f = TrivialFamily((3,4,7))
            sage: bool(f)
            True
            sage: g = Family([])
            sage: bool(g)
            False
        """
        return bool(self._enumeration)

    def __eq__(self, other):
        """
        TESTS::

            sage: f = Family((3,4,7))
            sage: g = Family([3,4,7])
            sage: f == g
            True
        """
        return (isinstance(other, self.__class__) and
                self._enumeration == other._enumeration)

    def __hash__(self):
        """
        Return a hash value for ``self``.

        EXAMPLES::

            sage: from sage.sets.family import TrivialFamily
            sage: f = TrivialFamily((3,4,7))
            sage: hash(f) == hash(f)
            True
        """
        return hash(self._enumeration)

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.sets.family import TrivialFamily
            sage: f = TrivialFamily([3,4,7]); f # indirect doctest
            Family (3, 4, 7)
        """
        return "Family %s" % ((self._enumeration),)

    def keys(self):
        """
        Return ``self``'s keys.

        EXAMPLES::

            sage: from sage.sets.family import TrivialFamily
            sage: f = TrivialFamily([3,4,7])
            sage: f.keys()
            [0, 1, 2]
        """
        return list(range(len(self._enumeration)))

    def cardinality(self):
        """
        Return the number of elements in ``self``.

        EXAMPLES::

            sage: from sage.sets.family import TrivialFamily
            sage: f = TrivialFamily([3,4,7])
            sage: f.cardinality()
            3
        """
        return Integer(len(self._enumeration))

    def __iter__(self):
        """
        EXAMPLES::

            sage: from sage.sets.family import TrivialFamily
            sage: f = TrivialFamily([3,4,7])
            sage: [i for i in f]
            [3, 4, 7]
        """
        return iter(self._enumeration)

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: from sage.sets.family import TrivialFamily
            sage: f = TrivialFamily([3,4,7])
            sage: 3 in f
            True
            sage: 5 in f
            False
        """
        return x in self._enumeration

    def __getitem__(self, i):
        """
        EXAMPLES::

            sage: from sage.sets.family import TrivialFamily
            sage: f = TrivialFamily([3,4,7])
            sage: f[1]
            4
        """
        return self._enumeration[i]

    def __getstate__(self):
        """
        TESTS::

            sage: from sage.sets.family import TrivialFamily
            sage: f = TrivialFamily([3,4,7])
            sage: f.__getstate__()
            {'_enumeration': (3, 4, 7)}
        """
        return {'_enumeration': self._enumeration}

    def __setstate__(self, state):
        """
        TESTS::

            sage: from sage.sets.family import TrivialFamily
            sage: f = TrivialFamily([3,4,7])
            sage: f.__setstate__({'_enumeration': (2, 4, 8)})
            sage: f
            Family (2, 4, 8)
        """
        self.__init__(state['_enumeration'])

    def map(self, f, name=None):
        r"""
        Return the family `( f(\mathtt{self}[i]) )_{i \in I}`,
        where `I` is the index set of ``self``.

        The result is again a :class:`TrivialFamily`.

        EXAMPLES::

            sage: from sage.sets.family import TrivialFamily
            sage: f = TrivialFamily(['a', 'b', 'd'])
            sage: g = f.map(lambda x: x + '1'); g
            Family ('a1', 'b1', 'd1')
        """
        # tuple([... for ...]) is faster than tuple(... for ...)
        return Family(tuple([f(x) for x in self._enumeration]), name=name)


class EnumeratedFamily(LazyFamily):
    r"""
    :class:`EnumeratedFamily` turns an enumerated set ``c`` into a family
    indexed by the set `\{0,\dots, |c|-1\}` (or ``NN`` if `|c|` is
    countably infinite).

    Instances should be created via the :func:`Family` factory. See its
    documentation for examples and tests.
    """
    def __init__(self, enumset):
        """
        EXAMPLES::

            sage: from sage.sets.family import EnumeratedFamily
            sage: f = EnumeratedFamily(Permutations(3))
            sage: TestSuite(f).run()

            sage: f = Family(NonNegativeIntegers())
            sage: TestSuite(f).run()

        TESTS:

        Check that category and keys are set correctly (:issue:`28274`)::

            sage: from sage.sets.family import EnumeratedFamily
            sage: f = EnumeratedFamily(Permutations(4))
            sage: f.category()
            Category of finite enumerated sets
            sage: list(f.keys()) == list(range(f.cardinality()))
            True
            sage: Family(Permutations()).keys()
            Non negative integers
            sage: type(Family(NN))
            <class 'sage.sets.family.EnumeratedFamily_with_category'>
        """
        if enumset.cardinality() == Infinity:
            baseset = NonNegativeIntegers()
        else:
            baseset = range(enumset.cardinality())
        LazyFamily.__init__(self, baseset, enumset.unrank)
        self.enumset = enumset

    def __eq__(self, other):
        """
        EXAMPLES::

            sage: f = Family(Permutations(3))
            sage: g = Family(Permutations(3))
            sage: f == g
            True
        """
        return (isinstance(other, self.__class__) and
                self.enumset == other.enumset)

    def __repr__(self):
        """
        EXAMPLES::

            sage: f = Family(Permutations(3)); f # indirect doctest
            Family (Standard permutations of 3)

            sage: f = Family(NonNegativeIntegers()); f
            Family (Non negative integers)
        """
#        return "Family ((%s)[i])_(i=1...%s)"%(self.enumset, self.enumset.cardinality())
        if isinstance(self.enumset, FiniteEnumeratedSet):
            return "Family %s" % (self.enumset._elements,)
        return "Family (%s)" % (self.enumset)

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: f = Family(Permutations(3))
            sage: [2,1,3] in f
            True
        """
        return x in self.enumset

    def cardinality(self):
        """
        Return the number of elements in ``self``.

        EXAMPLES::

            sage: from sage.sets.family import EnumeratedFamily
            sage: f = EnumeratedFamily(Permutations(3))
            sage: f.cardinality()
            6

            sage: f = Family(NonNegativeIntegers())
            sage: f.cardinality()
            +Infinity
        """
        return self.enumset.cardinality()

    def __iter__(self):
        """
        EXAMPLES::

            sage: from sage.sets.family import EnumeratedFamily
            sage: f = EnumeratedFamily(Permutations(3))
            sage: [i for i in f]
            [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]]
        """
        yield from self.enumset

    def __getitem__(self, i):
        """
        EXAMPLES::

            sage: from sage.sets.family import EnumeratedFamily
            sage: f = EnumeratedFamily(Permutations(3))
            sage: f[1]
            [1, 3, 2]
        """
        return self.enumset.unrank(i)

    def __getstate__(self):
        """
        EXAMPLES::

            sage: from sage.sets.family import EnumeratedFamily
            sage: f = EnumeratedFamily(Permutations(3))
            sage: f.__getstate__()
            {'enumset': Standard permutations of 3}
            sage: loads(dumps(f)) == f
            True
        """
        return {'enumset': self.enumset}

    def __setstate__(self, state):
        """
        EXAMPLES::

            sage: from sage.sets.family import EnumeratedFamily
            sage: f = EnumeratedFamily(Permutations(0))
            sage: f.__setstate__({'enumset': Permutations(3)})
            sage: f
            Family (Standard permutations of 3)
        """
        self.__init__(state['enumset'])
