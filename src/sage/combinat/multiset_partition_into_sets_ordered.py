r"""
Ordered Multiset Partitions into Sets and the Minimaj Crystal

This module provides element and parent classes for ordered multiset
partitions. It also implements the minimaj crystal of Benkart et al.
[BCHOPSY2017]_. (See :class:`MinimajCrystal`.)

AUTHORS:

- Aaron Lauve (2018): initial implementation. First draft of minimaj crystal
  code provided by Anne Schilling.

REFERENCES:

- [BCHOPSY2017]_
- [HRW2015]_
- [HRS2016]_
- [LM2018]_

EXAMPLES:

An ordered multiset partition into sets of the multiset `\{\{1, 3, 3, 5\}\}`::

    sage: OrderedMultisetPartitionIntoSets([[5, 3], [1, 3]])
    [{3,5}, {1,3}]

Ordered multiset partitions into sets of the multiset `\{\{1, 3, 3\}\}`::

    sage: OrderedMultisetPartitionsIntoSets([1,1,3]).list()
    [[{1}, {1}, {3}], [{1}, {1,3}], [{1}, {3}, {1}], [{1,3}, {1}], [{3}, {1}, {1}]]

Ordered multiset partitions into sets of the integer 4::

    sage: OrderedMultisetPartitionsIntoSets(4).list()
    [[{4}], [{1,3}], [{3}, {1}], [{1,2}, {1}], [{2}, {2}], [{2}, {1}, {1}],
     [{1}, {3}], [{1}, {1,2}], [{1}, {2}, {1}], [{1}, {1}, {2}], [{1}, {1}, {1}, {1}]]

Ordered multiset partitions into sets on the alphabet `\{1, 4\}` of order 3::

    sage: OrderedMultisetPartitionsIntoSets([1,4], 3).list()
    [[{1,4}, {1}], [{1,4}, {4}], [{1}, {1,4}], [{4}, {1,4}], [{1}, {1}, {1}],
     [{1}, {1}, {4}], [{1}, {4}, {1}], [{1}, {4}, {4}], [{4}, {1}, {1}],
     [{4}, {1}, {4}], [{4}, {4}, {1}], [{4}, {4}, {4}]]

Crystal of ordered multiset partitions into sets on the alphabet `\{1,2,3\}`
with 4 letters divided into 2 blocks::

    sage: crystals.Minimaj(3, 4, 2).list()                                              # needs sage.modules
    [((2, 3, 1), (1,)), ((2, 3), (1, 2)), ((2, 3), (1, 3)), ((2, 1), (1, 2)),
     ((3, 1), (1, 2)), ((3, 1, 2), (2,)), ((3, 1), (1, 3)), ((3, 1), (2, 3)),
     ((3, 2), (2, 3)), ((2, 1), (1, 3)), ((2,), (1, 2, 3)), ((3,), (1, 2, 3)),
     ((1,), (1, 2, 3)), ((1, 2), (2, 3)), ((1, 2, 3), (3,))]
"""

# ****************************************************************************
#       Copyright (C) 2018 Aaron Lauve       <lauve at math.luc.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from functools import reduce
from itertools import chain, product

from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.classical_crystals import ClassicalCrystals
from sage.categories.tensor import tensor
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.list_clone import ClonableArray
from sage.structure.parent import Parent
from sage.structure.element_wrapper import ElementWrapper
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.misc_c import prod, running_total
from sage.misc.latex import latex
from sage.misc.lazy_import import lazy_import
from sage.sets.set import Set_object
from sage.rings.infinity import infinity
from sage.rings.integer_ring import ZZ
from sage.rings.power_series_ring import PowerSeriesRing
from sage.arith.misc import binomial

from sage.combinat.subset import Subsets_sk
from sage.combinat.composition import Composition, Compositions, composition_iterator_fast
from sage.combinat.permutation import Permutations_mset
from sage.combinat.integer_lists.invlex import IntegerListsLex
from sage.combinat.combinatorial_map import combinatorial_map
from sage.combinat.shuffle import ShuffleProduct, ShuffleProduct_overlapping

lazy_import('sage.combinat.crystals.letters', 'CrystalOfLetters', as_='Letters')
lazy_import('sage.combinat.root_system.cartan_type', 'CartanType')
lazy_import('sage.combinat.sf.sf', 'SymmetricFunctions')


class OrderedMultisetPartitionIntoSets(ClonableArray,
        metaclass=InheritComparisonClasscallMetaclass):
    r"""
    Ordered Multiset Partition into sets.

    An *ordered multiset partition into sets* `c` of a multiset `X` is a list
    `[c_1, \ldots, c_r]` of nonempty subsets of `X` (note: not
    sub-multisets), called the *blocks* of `c`, whose multi-union is `X`.

    EXAMPLES:

    The simplest way to create an ordered multiset partition into sets is by
    specifying its blocks as a list or tuple::

        sage: OrderedMultisetPartitionIntoSets([[3],[2,1]])
        [{3}, {1,2}]
        sage: OrderedMultisetPartitionIntoSets(((3,), (1,2)))
        [{3}, {1,2}]
        sage: OrderedMultisetPartitionIntoSets([set([i]) for i in range(2,5)])
        [{2}, {3}, {4}]

    REFERENCES:

    - [HRW2015]_
    - [HRS2016]_
    - [LM2018]_
    """
    @staticmethod
    def __classcall_private__(cls, co):
        """
        Create an ordered multiset partition into sets (i.e., a list of sets)
        from the passed arguments with the appropriate parent.

        EXAMPLES::

            sage: OrderedMultisetPartitionIntoSets([[3], [2,1]])
            [{3}, {1,2}]
            sage: c = OrderedMultisetPartitionsIntoSets()([{2}, {3}, {4}, {5}]); c
            [{2}, {3}, {4}, {5}]
            sage: d = OrderedMultisetPartitionsIntoSets((1,1,1,2,3,5))([{1}, {5, 1, 3}, {2, 1}]); d
            [{1}, {1,3,5}, {1,2}]

        TESTS::

            sage: c.parent() == OrderedMultisetPartitionsIntoSets([2,3,4,5])
            False
            sage: d.parent() == OrderedMultisetPartitionsIntoSets([1,1,1,2,3,5])
            True
            sage: repr(OrderedMultisetPartitionIntoSets([]).parent())
            'Ordered Multiset Partitions into Sets of multiset {{}}'
        """
        if not co:
            P = OrderedMultisetPartitionsIntoSets([])
            return P.element_class(P, [])
        else:
            X = _concatenate(co)
            P = OrderedMultisetPartitionsIntoSets(_get_weight(X))
            return P.element_class(P, co)

    def __init__(self, parent, data):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: c = OrderedMultisetPartitionsIntoSets(7)([[1,3], [1,2]])
            sage: OrderedMultisetPartitionIntoSets([[1,3], [1,2]]) == c
            True
            sage: c.weight()
            {1: 2, 2: 1, 3: 1}

        TESTS::

            sage: OMP = OrderedMultisetPartitionIntoSets
            sage: c0 = OMP([])
            sage: OMP([[]]) == c0
            True
            sage: TestSuite(c0).run()

            sage: d = OMP([[1, 3], [1, 'a', 'b']])
            sage: TestSuite(d).run()

            sage: OMPs = OrderedMultisetPartitionsIntoSets()
            sage: d = OMPs([['a','b','c'],['a','b'],['a']])
            sage: TestSuite(d).run()

            sage: c.size() == 7
            True
            sage: d.size() == None
            True
        """
        # Delete empty blocks
        co = [block for block in data if block]
        if not _has_nonempty_sets(co):
            raise ValueError("cannot view %s as an ordered partition of %s" % (co, parent._Xtup))

        ClonableArray.__init__(self, parent, [frozenset(k) for k in co])
        self._multiset = _get_multiset(co)
        self._weight = _get_weight(self._multiset)
        self._order = sum(len(block) for block in self)
        if all((a in ZZ and a > 0) for a in self._multiset):
            self._n = ZZ(sum(self._multiset))
        else:
            self._n = None

    def check(self):
        """
        Check that we are a valid ordered multiset partition into sets.

        EXAMPLES::

            sage: c = OrderedMultisetPartitionsIntoSets(4)([[1], [1,2]])
            sage: c.check()

            sage: OMPs = OrderedMultisetPartitionsIntoSets()
            sage: c = OMPs([[1], [1], ['a']])
            sage: c.check()

        TESTS::

            sage: c = OMPs([[1, 1], [1, 4]])
            Traceback (most recent call last):
            ...
            ValueError: cannot convert [[1, 1], [1, 4]] into an element
             of Ordered Multiset Partitions into Sets
        """
        if self not in self.parent():
            raise ValueError("{} not an element of {}".format(self, self.parent()))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: A = OrderedMultisetPartitionIntoSets([[4], [1,2,4], [2,3], [1]])
            sage: A
            [{4}, {1,2,4}, {2,3}, {1}]
        """
        return self._repr_tight()

    def _repr_normal(self):
        r"""
        Viewing ``self`` as a list `[A_1, \ldots, A_r]` of sets,
        return the standard Sage string representation of `[A_1, \ldots, A_r]`.

        EXAMPLES::

            sage: OrderedMultisetPartitionIntoSets([[4,1,3], [3,2,5]])._repr_normal()
            '[{1, 3, 4}, {2, 3, 5}]'
            sage: OrderedMultisetPartitionIntoSets([[4,1,3,11], [3,'a',5]])._repr_normal()
            "[{1, 11, 3, 4}, {3, 5, 'a'}]"
        """
        # TODO: simplify if/once ``_repr_`` method for ``Set`` sorts its elements.
        if self._n:
            string_parts = (str(sorted(k)) for k in self)
        else:
            string_parts = (str(sorted(k, key=str)) for k in self)
        string_parts = ", ".join(string_parts).replace("[","{").replace("]","}")
        return "[" + string_parts + "]"

    def _repr_tight(self):
        r"""
        Starting from the standard Sage string representation of ``self``
        as a list `[A_1, \ldots, A_r]` of sets, return the shorter string
        gotten by deleting spaces within ``repr(A_i)``.

        EXAMPLES::

            sage: A = OrderedMultisetPartitionIntoSets([[4], [1,2,4], [2,3], [1]])
            sage: A._repr_normal()
            '[{4}, {1, 2, 4}, {2, 3}, {1}]'
            sage: A._repr_tight()
            '[{4}, {1,2,4}, {2,3}, {1}]'
        """
        repr = self._repr_normal()
        # eliminate spacing within blocks
        return repr.replace(", ", ",").replace("},{", "}, {")

    def __hash__(self):
        """
        Return the hash of ``self``.

        The parent is not included as part of the hash.

        EXAMPLES::

            sage: OMP = OrderedMultisetPartitionsIntoSets(4)
            sage: A = OMP([[1], [1, 2]])
            sage: B = OMP([{1}, {1, 2}])
            sage: hash(A) == hash(B)
            True
        """
        return sum(hash(x) for x in self)

    def __eq__(self, y):
        """
        Check equality of ``self`` and ``y``.

        The parent is not included as part of the equality check.

        TESTS::

            sage: OMP_n = OrderedMultisetPartitionsIntoSets(4)
            sage: OMP_X = OrderedMultisetPartitionsIntoSets([1,1,2])
            sage: OMP_Ad = OrderedMultisetPartitionsIntoSets(2, 3)
            sage: mu = [[1], [1, 2]]
            sage: OMP_n(mu) == OMP_X(mu) == OMP_Ad(mu)
            True
            sage: OMP_n(mu) == mu
            False
            sage: OMP_n(mu) == OMP_n([{1}, {3}])
            False
            sage: OMP_n(mu) == OMP_X([[1], [1,2]])
            True
        """
        if not isinstance(y, OrderedMultisetPartitionIntoSets):
            return False
        return list(self) == list(y)

    def __ne__(self, y):
        """
        Check lack of equality of ``self`` and ``y``.

        The parent is not included as part of the equality check.

        TESTS::

            sage: OMP = OrderedMultisetPartitionsIntoSets(4)
            sage: mu = [[1], [1, 2]]
            sage: OMP(mu).__ne__(mu)
            True
            sage: nu = [[1], [2], [1]]
            sage: OMP(mu).__ne__(OMP(nu))
            True
        """
        return not (self == y)

    def __add__(self, other):
        """
        Return the concatenation of two ordered multiset partitions into sets.

        This operation represents the product in Hopf algebra of ordered multiset
        partitions into sets in its natural basis [LM2018]_.

        EXAMPLES::

            sage: OMP = OrderedMultisetPartitionIntoSets
            sage: OMP([[1],[1],[1,3]]) + OMP([[4,1],[2]])
            [{1}, {1}, {1,3}, {1,4}, {2}]

        TESTS::

            sage: OMP([]) + OMP([]) == OMP([])
            True
            sage: OMP([[1],[1],[1,3]]) + OMP([]) == OMP([[1],[1],[1,3]])
            True
        """
        co = list(self) + list(other)
        X = _concatenate(co)
        return OrderedMultisetPartitionsIntoSets(_get_weight(X))(co)

    @combinatorial_map(order=2, name='reversal')
    def reversal(self):
        r"""
        Return the reverse ordered multiset partition into sets of ``self``.

        Given an ordered multiset partition into sets `(B_1, B_2, \ldots, B_k)`,
        its reversal is defined to be the ordered multiset partition into sets
        `(B_k, \ldots, B_2, B_1)`.

        EXAMPLES::

            sage: C = OrderedMultisetPartitionIntoSets([[1], [1, 3], [2, 3, 4]]); C
            [{1}, {1,3}, {2,3,4}]
            sage: C.reversal()
            [{2,3,4}, {1,3}, {1}]
        """
        return self.parent()(list(reversed(self)))

    def shape_from_cardinality(self):
        """
        Return a composition that records the cardinality of each block of ``self``.

        EXAMPLES::

            sage: C = OrderedMultisetPartitionIntoSets([[3, 4, 1], [2], [1, 2, 3, 7]]); C
            [{1,3,4}, {2}, {1,2,3,7}]
            sage: C.shape_from_cardinality()
            [3, 1, 4]
            sage: OrderedMultisetPartitionIntoSets([]).shape_from_cardinality() == Composition([])
            True
        """
        return Composition([len(k) for k in self])

    def shape_from_size(self):
        """
        Return a composition that records the sum of entries of each
        block of ``self``.

        EXAMPLES::

            sage: C = OrderedMultisetPartitionIntoSets([[3, 4, 1], [2], [1, 2, 3, 7]]); C
            [{1,3,4}, {2}, {1,2,3,7}]
            sage: C.shape_from_size()
            [8, 2, 13]

        TESTS::

            sage: OrderedMultisetPartitionIntoSets([]).shape_from_size() == Composition([])
            True
            sage: D = OrderedMultisetPartitionIntoSets([['a', 'b'], ['a']]); D
            [{'a','b'}, {'a'}]
            sage: D.shape_from_size() == None
            True
        """
        if self._n is not None:
            return Composition([sum(k) for k in self])

    def letters(self):
        """
        Return the set of distinct elements occurring within the blocks
        of ``self``.

        EXAMPLES::

            sage: C = OrderedMultisetPartitionIntoSets([[3, 4, 1], [2], [1, 2, 3, 7]]); C
            [{1,3,4}, {2}, {1,2,3,7}]
            sage: C.letters()
            frozenset({1, 2, 3, 4, 7})
        """
        return _union_of_sets(list(self))

    def multiset(self, as_dict=False):
        """
        Return the multiset corresponding to ``self``.

        INPUT:

        - ``as_dict`` -- boolean (default: ``False``); whether to return the multiset
          as a tuple of a dict of multiplicities

        EXAMPLES::

            sage: C = OrderedMultisetPartitionIntoSets([[3, 4, 1], [2], [1, 2, 3, 7]]); C
            [{1,3,4}, {2}, {1,2,3,7}]
            sage: C.multiset()
            (1, 1, 2, 2, 3, 3, 4, 7)
            sage: C.multiset(as_dict=True)
            {1: 2, 2: 2, 3: 2, 4: 1, 7: 1}
            sage: OrderedMultisetPartitionIntoSets([]).multiset() == ()
            True
        """
        if as_dict:
            return self._weight
        else:
            return self._multiset

    def max_letter(self):
        """
        Return the maximum letter appearing in ``self.letters()`` of ``self``.

        EXAMPLES::

            sage: C = OrderedMultisetPartitionIntoSets([[3, 4, 1], [2], [1, 2, 3, 7]])
            sage: C.max_letter()
            7
            sage: D = OrderedMultisetPartitionIntoSets([['a','b','c'],['a','b'],['a'],['b','c','f'],['c','d']])
            sage: D.max_letter()
            'f'
            sage: C = OrderedMultisetPartitionIntoSets([])
            sage: C.max_letter()
        """
        if not self.letters():
            return None
        else:
            return max(self.letters())

    def size(self):
        """
        Return the size of ``self`` (that is, the sum of all integers in
        all blocks) if ``self`` is a list of subsets of positive integers.

        Else, return ``None``.

        EXAMPLES::

            sage: C = OrderedMultisetPartitionIntoSets([[3, 4, 1], [2], [1, 2, 3, 7]]); C
            [{1,3,4}, {2}, {1,2,3,7}]
            sage: C.size()
            23
            sage: C.size() == sum(k for k in C.shape_from_size())
            True
            sage: OrderedMultisetPartitionIntoSets([[7,1],[3]]).size()
            11

        TESTS::

            sage: OrderedMultisetPartitionIntoSets([]).size() == 0
            True
            sage: OrderedMultisetPartitionIntoSets([['a','b'],['a','b','c']]).size() is None
            True
        """
        return self._n

    def order(self):
        """
        Return the total number of elements in all blocks of ``self``.

        EXAMPLES::

            sage: C = OrderedMultisetPartitionIntoSets([[3, 4, 1], [2], [1, 2, 3, 7]]); C
            [{1,3,4}, {2}, {1,2,3,7}]
            sage: C.order()
            8
            sage: C.order() == sum(C.weight().values())
            True
            sage: C.order() == sum(k for k in C.shape_from_cardinality())
            True
            sage: OrderedMultisetPartitionIntoSets([[7,1],[3]]).order()
            3
        """
        return self._order

    def length(self):
        """
        Return the number of blocks of ``self``.

        EXAMPLES::

            sage: OrderedMultisetPartitionIntoSets([[7,1],[3]]).length()
            2
        """
        return len(self)

    def weight(self, as_weak_comp=False):
        r"""
        Return a dictionary, with keys being the letters in ``self.letters()``
        and values being their (positive) frequency.

        Alternatively, if ``as_weak_comp`` is ``True``, count the number of instances
        `n_i` for each distinct positive integer `i` across all blocks of ``self``.
        Return as a list `[n_1, n_2, n_3, ..., n_k]`, where `k` is the max letter
        appearing in ``self.letters()``.

        EXAMPLES::

            sage: c = OrderedMultisetPartitionIntoSets([[6,1],[1,3],[1,3,6]])
            sage: c.weight()
            {1: 3, 3: 2, 6: 2}
            sage: c.weight(as_weak_comp=True)
            [3, 0, 2, 0, 0, 2]

        TESTS::

            sage: OrderedMultisetPartitionIntoSets([]).weight() == {}
            True

            sage: c = OrderedMultisetPartitionIntoSets([['a','b'],['a','b','c'],['b'],['b'],['c']])
            sage: c.weight()
            {'a': 2, 'b': 4, 'c': 2}
            sage: c.weight(as_weak_comp=True)
            Traceback (most recent call last):
            ...
            ValueError: {'a': 2, 'b': 4, 'c': 2} is not a numeric multiset
        """
        from pprint import pformat
        w = self._weight
        if as_weak_comp:
            if all(v in ZZ for v in w):
                w = [w.get(i, 0) for i in range(1, self.max_letter() + 1)]
            else:
                raise ValueError("%s is not a numeric multiset" % pformat(w))
        return w

    def deconcatenate(self, k=2):
        r"""
        Return the list of `k`-deconcatenations of ``self``.

        A `k`-tuple `(C_1, \ldots, C_k)` of ordered multiset partitions into sets
        represents a `k`-deconcatenation of an ordered multiset partition into sets
        `C` if `C_1 + \cdots + C_k = C`.

        .. NOTE::

            This is not to be confused with ``self.split_blocks()``,
            which splits each block of ``self`` before making `k`-tuples
            of ordered multiset partitions into sets.

        EXAMPLES::

            sage: OrderedMultisetPartitionIntoSets([[7,1],[3,4,5]]).deconcatenate()
            [([{1,7}, {3,4,5}], []), ([{1,7}], [{3,4,5}]), ([], [{1,7}, {3,4,5}])]
            sage: OrderedMultisetPartitionIntoSets([['b','c'],['a']]).deconcatenate()
            [([{'b','c'}, {'a'}], []), ([{'b','c'}], [{'a'}]), ([], [{'b','c'}, {'a'}])]
            sage: OrderedMultisetPartitionIntoSets([['a','b','c']]).deconcatenate(3)
            [([{'a','b','c'}], [], []),
             ([], [{'a','b','c'}], []),
             ([], [], [{'a','b','c'}])]

        TESTS::

            sage: C = OrderedMultisetPartitionIntoSets([['a'],['b'],['c'],['d'],['e']]); C
            [{'a'}, {'b'}, {'c'}, {'d'}, {'e'}]
            sage: all( len(C.deconcatenate(k))
            ....:      == binomial(C.length() + k-1, k-1)
            ....:      for k in range(1, 5) )
            True
        """
        P = OrderedMultisetPartitionsIntoSets(alphabet=self.letters(),
                                              max_length=self.length())
        out = []
        for c in IntegerListsLex(self.length(), length=k):
            ps = [sum(c[:i]) for i in range(k+1)]
            out.append(tuple([P(self[ps[i]:ps[i+1]]) for i in range(len(ps)-1)]))
        return out

    def split_blocks(self, k=2):
        r"""
        Return a dictionary representing the `k`-splittings of ``self``.

        A `k`-tuple `(A^1, \ldots, A^k)` of ordered multiset partitions into sets
        represents a `k`-splitting of an ordered multiset partition into sets
        `A = [b_1, \ldots, b_r]` if one can express each block `b_i` as
        an (ordered) disjoint union of sets `b_i = b^1_i \sqcup \cdots
        \sqcup b^k_i` (some possibly empty) so that each `A^j` is the
        ordered multiset partition into sets corresponding to the list `[b^j_1,
        b^j_2, \ldots, b^j_r]`, excising empty sets appearing therein.

        This operation represents the coproduct in Hopf algebra of ordered
        multiset partitions into sets in its natural basis [LM2018]_.

        EXAMPLES::

            sage: sorted(OrderedMultisetPartitionIntoSets([[1,2],[3,4]]).split_blocks(), key=str)
            [([], [{1,2}, {3,4}]),
             ([{1,2}, {3,4}], []),
             ([{1,2}, {3}], [{4}]),
             ([{1,2}, {4}], [{3}]),
             ([{1,2}], [{3,4}]),
             ([{1}, {3,4}], [{2}]),
             ([{1}, {3}], [{2}, {4}]),
             ([{1}, {4}], [{2}, {3}]),
             ([{1}], [{2}, {3,4}]),
             ([{2}, {3,4}], [{1}]),
             ([{2}, {3}], [{1}, {4}]),
             ([{2}, {4}], [{1}, {3}]),
             ([{2}], [{1}, {3,4}]),
             ([{3,4}], [{1,2}]),
             ([{3}], [{1,2}, {4}]),
             ([{4}], [{1,2}, {3}])]
            sage: sorted(OrderedMultisetPartitionIntoSets([[1,2]]).split_blocks(3), key=str)
            [([], [], [{1,2}]), ([], [{1,2}], []), ([], [{1}], [{2}]),
             ([], [{2}], [{1}]), ([{1,2}], [], []), ([{1}], [], [{2}]),
             ([{1}], [{2}], []), ([{2}], [], [{1}]), ([{2}], [{1}], [])]
            sage: OrderedMultisetPartitionIntoSets([[4],[4]]).split_blocks()
            {([], [{4}, {4}]): 1, ([{4}], [{4}]): 2, ([{4}, {4}], []): 1}

        TESTS::

            sage: C = OrderedMultisetPartitionIntoSets([[1,2],[4,5,6]]); C
            [{1,2}, {4,5,6}]
            sage: sum(C.split_blocks().values()) == 2**len(C[0]) * 2**len(C[1])
            True
            sage: sum(C.split_blocks(3).values()) == (1+2)**len(C[0]) * (1+2)**len(C[1])
            True
            sage: C = OrderedMultisetPartitionIntoSets([])
            sage: C.split_blocks(3) == {(C, C, C): 1}
            True
        """
        P = OrderedMultisetPartitionsIntoSets(alphabet=self.letters(),
                                              max_length=self.length())

        # corner case
        if not self:
            return {tuple([self]*k): 1}

        out = {}
        for t in product(*[_split_block(block, k) for block in self]):
            tt = tuple([P([l for l in c if l]) for c in zip(*t)])
            out[tt] = out.get(tt, 0) + 1
        return out

    def finer(self, strong=False):
        r"""
        Return the set of ordered multiset partitions into sets that are finer
        than ``self``.

        An ordered multiset partition into sets `A` is finer than another `B`
        if, reading left-to-right, every block of `B` is the union of some
        consecutive blocks of `A`.

        If optional argument ``strong`` is set to ``True``, then return
        only those `A` whose blocks are deconcatenations of blocks of `B`.
        (Here, we view blocks of `B` as sorted lists instead of sets.)

        EXAMPLES::

            sage: C = OrderedMultisetPartitionIntoSets([[3,2]]).finer()
            sage: len(C)
            3
            sage: sorted(C, key=str)
            [[{2,3}], [{2}, {3}], [{3}, {2}]]
            sage: OrderedMultisetPartitionIntoSets([]).finer()
            {[]}
            sage: O = OrderedMultisetPartitionsIntoSets([1, 1, 'a', 'b'])
            sage: o = O([{1}, {'a', 'b'}, {1}])
            sage: sorted(o.finer(), key=str)
            [[{1}, {'a','b'}, {1}], [{1}, {'a'}, {'b'}, {1}], [{1}, {'b'}, {'a'}, {1}]]
            sage: o.finer() & o.fatter() == set([o])
            True
        """
        P = OrderedMultisetPartitionsIntoSets(self._multiset)

        if not self:
            return set([self])

        CP = product(*[_refine_block(block, strong) for block in self])
        return set(P(_concatenate(map(list, c))) for c in CP)

    def is_finer(self, co):
        """
        Return ``True`` if the ordered multiset partition into sets ``self``
        is finer than the composition ``co``; otherwise, return ``False``.

        EXAMPLES::

            sage: OrderedMultisetPartitionIntoSets([[4],[1],[2]]).is_finer([[1,4],[2]])
            True
            sage: OrderedMultisetPartitionIntoSets([[1],[4],[2]]).is_finer([[1,4],[2]])
            True
            sage: OrderedMultisetPartitionIntoSets([[1,4],[1],[1]]).is_finer([[1,4],[2]])
            False
        """
        X = _concatenate(co)
        if self.weight() != OrderedMultisetPartitionsIntoSets(_get_weight(X))(co).weight():
            return False

        # trim common prefix and suffix to make the search-space smaller
        co1 = list(map(set, self))
        co2 = list(map(set, co))
        while co1[0] == co2[0]:
            co1 = co1[1:]
            co2 = co2[1:]
        while co1[-1] == co2[-1]:
            co1 = co1[:-1]
            co2 = co2[:-1]

        co1 = OrderedMultisetPartitionIntoSets(co1)
        co2 = OrderedMultisetPartitionIntoSets(co2)
        return co1 in co2.finer()

    def fatten(self, grouping):
        r"""
        Return the ordered multiset partition into sets fatter than ``self``,
        obtained by grouping together consecutive parts according to ``grouping``
        (whenever this does not violate the strictness condition).

        INPUT:

        - ``grouping`` -- a composition (or list) whose sum is the length
          of ``self``

        EXAMPLES:

        Let us start with the composition::

            sage: C = OrderedMultisetPartitionIntoSets([[4,1,5], [2], [7,1]]); C
            [{1,4,5}, {2}, {1,7}]

        With ``grouping`` equal to `(1, 1, 1)`, `C` is left unchanged::

            sage: C.fatten([1,1,1])
            [{1,4,5}, {2}, {1,7}]

        With ``grouping`` equal to `(2,1)` or `(1,2)`, a union of consecutive
        parts is achieved::

            sage: C.fatten([2,1])
            [{1,2,4,5}, {1,7}]
            sage: C.fatten([1,2])
            [{1,4,5}, {1,2,7}]

        However, the ``grouping`` `(3)` will throw an error, as `1` cannot
        appear twice in any block of ``C``::

            sage: C.fatten(Composition([3]))
            Traceback (most recent call last):
            ...
            ValueError: [{1,4,5,2,1,7}] is not a valid ordered multiset partition into sets
        """
        if sum(list(grouping)) != self.length():
            raise ValueError("%s is not a composition of ``self.length()`` (=%s)"
                             % (grouping, self.length()))

        valid = True
        result = []
        for i in range(len(grouping)):
            result_i = self[sum(grouping[:i]) : sum(grouping[:i+1])]
            # check that grouping[i] is allowed, i.e., `|A\cup B| = |A| + |B|`
            strict_size = sum(map(len, result_i))
            size = len(_union_of_sets(result_i))
            if size < strict_size:
                valid = False
            result.append(_concatenate(result_i))
        if not valid:
            str_rep = '['
            for i in range(len(grouping)):
                st = ",".join(str(k) for k in result[i])
                str_rep += "{" + st + "}"
            str_rep = str_rep.replace("}{", "}, {") + "]"
            raise ValueError("%s is not a valid ordered multiset partition into sets" % (str_rep))
        else:
            return OrderedMultisetPartitionsIntoSets(self._multiset)(result)

    def fatter(self):
        """
        Return the set of ordered multiset partitions into sets which are fatter
        than ``self``.

        An ordered multiset partition into sets `A` is fatter than another `B`
        if, reading left-to-right, every block of `A` is the union of some
        consecutive blocks of `B`.

        EXAMPLES::

            sage: C = OrderedMultisetPartitionIntoSets([{1,4,5}, {2}, {1,7}]).fatter()
            sage: len(C)
            3
            sage: sorted(C)
            [[{1,4,5}, {2}, {1,7}], [{1,4,5}, {1,2,7}], [{1,2,4,5}, {1,7}]]
            sage: sorted(OrderedMultisetPartitionIntoSets([['a','b'],['c'],['a']]).fatter())
            [[{'a','b'}, {'c'}, {'a'}], [{'a','b'}, {'a','c'}], [{'a','b','c'}, {'a'}]]

        Some extreme cases::

            sage: list(OrderedMultisetPartitionIntoSets([['a','b','c']]).fatter())
            [[{'a','b','c'}]]
            sage: list(OrderedMultisetPartitionIntoSets([]).fatter())
            [[]]
            sage: A = OrderedMultisetPartitionIntoSets([[1], [2], [3], [4]])
            sage: B = OrderedMultisetPartitionIntoSets([[1,2,3,4]])
            sage: A.fatter().issubset(B.finer())
            True
        """
        out = set()
        for c in composition_iterator_fast(self.length()):
            try:
                out.add(self.fatten(c))
            except ValueError:
                pass
        return out

    def minimaj(self):
        r"""
        Return the minimaj statistic on ordered multiset partitions into sets.

        We define `minimaj` via an example:

        1. Sort the block in ``self`` as prescribed by ``self.minimaj_word()``,
           keeping track of the original separation into blocks::

             in:   [{1,5,7}, {2,4}, {5,6}, {4,6,8}, {1,3}, {1,2,3}]
             out:  ( 5,7,1 /  2,4 /  5,6 /  4,6,8 /  3,1 /  1,2,3 )

        2. Record the indices where descents in this word occur::

             word:      (5, 7, 1 / 2, 4 / 5, 6 / 4, 6, 8 / 3, 1 / 1, 2, 3)
             indices:    1  2  3   4  5   6  7   8  9 10  11 12  13 14 15
             descents:  {   2,               7,       10, 11             }

        3. Compute the sum of the descents::

             minimaj = 2 + 7 + 10 + 11 = 30

        REFERENCES:

        - [HRW2015]_

        EXAMPLES::

            sage: C = OrderedMultisetPartitionIntoSets([{1,5,7}, {2,4}, {5,6}, {4,6,8}, {1,3}, {1,2,3}])
            sage: C, C.minimaj_word()
            ([{1,5,7}, {2,4}, {5,6}, {4,6,8}, {1,3}, {1,2,3}],
             (5, 7, 1, 2, 4, 5, 6, 4, 6, 8, 3, 1, 1, 2, 3))
            sage: C.minimaj()
            30
            sage: C = OrderedMultisetPartitionIntoSets([{2,4}, {1,2,3}, {1,6,8}, {2,3}])
            sage: C, C.minimaj_word()
            ([{2,4}, {1,2,3}, {1,6,8}, {2,3}], (2, 4, 1, 2, 3, 6, 8, 1, 2, 3))
            sage: C.minimaj()
            9
            sage: OrderedMultisetPartitionIntoSets([]).minimaj()
            0
            sage: C = OrderedMultisetPartitionIntoSets([['b','d'],['a','b','c'],['b']])
            sage: C, C.minimaj_word()
            ([{'b','d'}, {'a','b','c'}, {'b'}], ('d', 'b', 'c', 'a', 'b', 'b'))
            sage: C.minimaj()
            4
        """
        D = _descents(self.minimaj_word())
        return sum(D) + len(D)

    def minimaj_word(self):
        """
        Return an ordering of ``self._multiset`` derived from the minimaj
        ordering on blocks of ``self``.

        .. SEEALSO::

            :meth:`OrderedMultisetPartitionIntoSets.minimaj_blocks()`.

        EXAMPLES::

            sage: C = OrderedMultisetPartitionIntoSets([[2,1], [1,2,3], [1,2], [3], [1]]); C
            [{1,2}, {1,2,3}, {1,2}, {3}, {1}]
            sage: C.minimaj_blocks()
            ((1, 2), (2, 3, 1), (1, 2), (3,), (1,))
            sage: C.minimaj_word()
            (1, 2, 2, 3, 1, 1, 2, 3, 1)
        """
        return _concatenate(self.minimaj_blocks())

    def minimaj_blocks(self):
        r"""
        Return the minimaj ordering on blocks of ``self``.

        We define the ordering via the example below.

        Sort the blocks `[B_1,...,B_k]` of ``self`` from right to left via:

        1. Sort the last block `B_k` in increasing order, call it the word `W_k`

        2. If blocks `B_{i+1}, \ldots, B_k` have been converted to words
           `W_{i+1}, \ldots, W_k`, use the letters in `B_i` to make the unique
           word `W_i` that has a factorization `W_i = (u, v)` satisfying:

           - letters of `u` and `v` appear in increasing order, with `v`
             possibly empty;
           - letters in `vu` appear in increasing order;
           - ``v[-1]`` is the largest letter `a \in B_i` satisfying
             ``a <= W_{i+1}[0]``.

        EXAMPLES::

            sage: OrderedMultisetPartitionIntoSets([[1,5,7], [2,4], [5,6], [4,6,8], [1,3], [1,2,3]])
            [{1,5,7}, {2,4}, {5,6}, {4,6,8}, {1,3}, {1,2,3}]
            sage: _.minimaj_blocks()
            ((5, 7, 1), (2, 4), (5, 6), (4, 6, 8), (3, 1), (1, 2, 3))
            sage: OrderedMultisetPartitionIntoSets([]).minimaj_blocks()
            ()
        """
        if not self:
            return ()

        C = [sorted(self[-1])]
        for i in range(1, len(self)):
            lower = []
            upper = []
            for j in self[-1 - i]:
                if j <= C[0][0]:
                    lower.append(j)
                else:
                    upper.append(j)
            C = [sorted(upper) + sorted(lower)] + C
        return tuple(map(tuple, C))

    def to_tableaux_words(self):
        r"""
        Return a sequence of lists corresponding to row words
        of (skew-)tableaux.

        OUTPUT:

        The minimaj bijection `\phi` of [BCHOPSY2017]_
        applied to ``self``.

        .. TODO::

            Implement option for mapping to sequence of (skew-)tableaux?

        EXAMPLES::

            sage: co = ((1,2,4),(4,5),(3,),(4,6,1),(2,3,1),(1,),(2,5))
            sage: OrderedMultisetPartitionIntoSets(co).to_tableaux_words()
            [[5, 1], [3, 1], [6], [5, 4, 2], [1, 4, 3, 4, 2, 1, 2]]
        """
        if not self:
            return []
        bb = self.minimaj_blocks()
        b = [block[0] for block in bb]
        beginning = [0]+running_total(self.shape_from_cardinality())
        w = _concatenate(bb)
        D = [0] + _descents(w) + [len(w)]
        pieces = [b]
        for i in range(len(D)-1):
            p = [w[j] for j in range(D[i]+1,D[i+1]+1) if j not in beginning]
            pieces = [p[::-1]] + pieces
        return pieces

    def major_index(self):
        r"""
        Return the major index of ``self``.

        The major index is a statistic on ordered multiset partitions into sets,
        which we define here via an example.

        1. Sort each block in the list ``self`` in descending order to create
           a word `w`, keeping track of the original separation into blocks::

             in:  [{3,4,5}, {2,3,4}, {1}, {4,5}]
             out: [ 5,4,3 /  4,3,2 /  1 /  5,4 ]

        2. Create a sequence `v = (v_0, v_1, v_2, \ldots)`  of length
           ``self.order()+1``, built recursively by:

           1. `v_0 = 0`
           2. `v_j = v_{j-1} + \delta(j)`, where `\delta(j) = 1` if `j` is
              the index of an end of a block, and zero otherwise.

           ::

             in:    [ 5,4,3 /  4,3,2 /  1 /  5,4]
             out: (0, 0,0,1,   1,1,2,   3,   3,4)

        3. Compute `\sum_j v_j`, restricted to descent positions in `w`, i.e.,
           sum over those `j` with `w_j > w_{j+1}`::

             in:  w:   [5, 4, 3, 4, 3, 2, 1, 5, 4]
                  v: (0 0, 0, 1, 1, 1, 2, 3, 3, 4)
             maj :=     0 +0    +1 +1 +2    +3     = 7

        REFERENCES:

        - [HRW2015]_

        EXAMPLES::

            sage: C = OrderedMultisetPartitionIntoSets([{1,5,7}, {2,4}, {5,6}, {4,6,8}, {1,3}, {1,2,3}])
            sage: C.major_index()
            27
            sage: C = OrderedMultisetPartitionIntoSets([{3,4,5}, {2,3,4}, {1}, {4,5}])
            sage: C.major_index()
            7
        """
        ew = [enumerate(sorted(k)) for k in self]
        w = []
        v = [0]
        for eblock in ew:
            for (i,wj) in sorted(eblock, reverse=True):
                vj = v[-1]
                if i == 0:
                    vj += 1
                v.append(vj)
                w.append(wj)
        maj = [v[j+1] for j in range(len(w)-1) if w[j] > w[j+1]]
        return sum(maj)

    def shuffle_product(self, other, overlap=False):
        r"""
        Return the shuffles (with multiplicity) of blocks of ``self``
        with blocks of ``other``.

        In case optional argument ``overlap`` is ``True``, instead return
        the allowable overlapping shuffles. An overlapping shuffle `C` is
        allowable if, whenever one of its blocks `c` comes from the union
        `c = a \cup b` of a block of ``self`` and a block of ``other``,
        then this union is disjoint.

        .. SEEALSO::

            :meth:`Composition.shuffle_product()`

        EXAMPLES::

            sage: A = OrderedMultisetPartitionIntoSets([[2,1,3], [1,2]]); A
            [{1,2,3}, {1,2}]
            sage: B = OrderedMultisetPartitionIntoSets([[3,4]]); B
            [{3,4}]
            sage: C = OrderedMultisetPartitionIntoSets([[4,5]]); C
            [{4,5}]
            sage: list(A.shuffle_product(B))
            [[{1,2,3}, {1,2}, {3,4}], [{3,4}, {1,2,3}, {1,2}], [{1,2,3}, {3,4}, {1,2}]]
            sage: list(A.shuffle_product(B, overlap=True))
            [[{1,2,3}, {1,2}, {3,4}], [{1,2,3}, {3,4}, {1,2}],
             [{3,4}, {1,2,3}, {1,2}], [{1,2,3}, {1,2,3,4}]]
            sage: list(A.shuffle_product(C, overlap=True))
            [[{1,2,3}, {1,2}, {4,5}], [{1,2,3}, {4,5}, {1,2}], [{4,5}, {1,2,3}, {1,2}],
             [{1,2,3,4,5}, {1,2}], [{1,2,3}, {1,2,4,5}]]
        """
        other = OrderedMultisetPartitionIntoSets(other)
        P = OrderedMultisetPartitionsIntoSets(self._multiset + other._multiset)
        if not overlap:
            for term in ShuffleProduct(self, other, element_constructor=P):
                yield term
        else:
            A = list(map(tuple, self))
            B = list(map(tuple, other))
            for term in ShuffleProduct_overlapping(A, B):
                if len(_concatenate(map(frozenset, term))) == len(P._Xtup):
                    yield P(term)

##############################################################


class OrderedMultisetPartitionsIntoSets(UniqueRepresentation, Parent):
    r"""
    Ordered Multiset Partitions into Sets.

    An *ordered multiset partition into sets* `c` of a multiset `X` is
    a list of nonempty subsets (not multisets), called the *blocks* of `c`,
    whose multi-union is `X`.

    The number of blocks of `c` is called its *length*. The *order* of `c`
    is the cardinality of the multiset `X`. If, additionally, `X` is a
    multiset of positive integers, then the *size* of `c` is the sum of
    all elements of `X`.

    The user may wish to focus on ordered multiset partitions into sets
    of a given size, or over a given alphabet. Hence, this class allows
    a variety of arguments as input.

    INPUT:

    Expects one or two arguments, with different behaviors resulting:

    - One Argument:

      + `X` -- a dictionary or list or tuple
        (representing a multiset for `c`),
        or an integer (representing the size of `c`)

    - Two Arguments:

      + `A` -- list (representing allowable letters within blocks of `c`),
        or a positive integer (representing the maximal allowable letter)
      + `n` -- a nonnegative integer (the total number of letters within `c`)

    Optional keyword arguments are as follows:
    (See corresponding methods in see :class:`OrderedMultisetPartitionIntoSets` for more details.)

    - ``weight=X``     (list or dictionary `X`) specifies the multiset for `c`
    - ``size=n``       (integer `n`) specifies the size of `c`
    - ``alphabet=A``   (iterable `A`) specifies allowable elements for the blocks of `c`
    - ``length=k``     (integer `k`) specifies the number of blocks in the partition
    - ``min_length=k`` (integer `k`) specifies minimum number of blocks in the partition
    - ``max_length=k`` (integer `k`) specifies maximum number of blocks in the partition
    - ``order=n``      (integer `n`) specifies the cardinality of the multiset that `c` partitions
    - ``min_order=n``  (integer `n`) specifies minimum number of elements in the partition
    - ``max_order=n``  (integer `n`) specifies maximum number of elements in the partition

    EXAMPLES:

    Passing one argument to :class:`OrderedMultisetPartitionsIntoSets`:

    There are 5 ordered multiset partitions into sets of the multiset
    `\{\{1, 1, 4\}\}`::

        sage: OrderedMultisetPartitionsIntoSets([1,1,4]).cardinality()
        5

    Here is the list of them::

        sage: OrderedMultisetPartitionsIntoSets([1,1,4]).list()
        [[{1}, {1}, {4}], [{1}, {1,4}], [{1}, {4}, {1}], [{1,4}, {1}], [{4}, {1}, {1}]]

    By chance, there are also 5 ordered multiset partitions into sets of
    the integer 3::

        sage: OrderedMultisetPartitionsIntoSets(3).cardinality()
        5

    Here is the list of them::

        sage: OrderedMultisetPartitionsIntoSets(3).list()
        [[{3}], [{1,2}], [{2}, {1}], [{1}, {2}], [{1}, {1}, {1}]]

    Passing two arguments to :class:`OrderedMultisetPartitionsIntoSets`:

    There are also 5 ordered multiset partitions into sets of order 2
    over the alphabet `\{1, 4\}`::

        sage: OrderedMultisetPartitionsIntoSets([1, 4], 2)
        Ordered Multiset Partitions into Sets of order 2 over alphabet {1, 4}
        sage: OrderedMultisetPartitionsIntoSets([1, 4], 2).cardinality()
        5

    Here is the list of them::

        sage: OrderedMultisetPartitionsIntoSets([1, 4], 2).list()
        [[{1,4}], [{1}, {1}], [{1}, {4}], [{4}, {1}], [{4}, {4}]]

    If no arguments are passed to :class:`OrderedMultisetPartitionsIntoSets`,
    then the code returns all ordered multiset partitions into sets::

        sage: OrderedMultisetPartitionsIntoSets()
        Ordered Multiset Partitions into Sets
        sage: [] in OrderedMultisetPartitionsIntoSets()
        True
        sage: [[2,3], [1]] in OrderedMultisetPartitionsIntoSets()
        True
        sage: [['a','b'], ['a']] in OrderedMultisetPartitionsIntoSets()
        True
        sage: [[-2,3], [3]] in OrderedMultisetPartitionsIntoSets()
        True
        sage: [[2], [3,3]] in OrderedMultisetPartitionsIntoSets()
        False

    The following examples show how to test whether or not an object
    is an ordered multiset partition into sets::

        sage: [[3,2],[2]] in OrderedMultisetPartitionsIntoSets()
        True
        sage: [[3,2],[2]] in OrderedMultisetPartitionsIntoSets(7)
        True
        sage: [[3,2],[2]] in OrderedMultisetPartitionsIntoSets([2,2,3])
        True
        sage: [[3,2],[2]] in OrderedMultisetPartitionsIntoSets(5)
        False

    .. RUBRIC:: Optional keyword arguments

    Passing keyword arguments that are incompatible with required requirements
    results in an error; otherwise, the collection of ordered multiset partitions
    into sets is restricted accordingly:

    *The* ``weight`` *keyword:*

    This is used to specify which multiset `X` is to be considered,
    if this multiset was not passed as one of the required arguments for
    :class:`OrderedMultisetPartitionsIntoSets`. In principle, it is a dictionary,
    but weak compositions are also allowed. For example, the ordered multiset
    partitions into sets of integer 4 are listed by weight below::

        sage: OrderedMultisetPartitionsIntoSets(4, weight=[0,0,0,1])
        Ordered Multiset Partitions into Sets of integer 4 with constraint: weight={4: 1}
        sage: OrderedMultisetPartitionsIntoSets(4, weight=[0,0,0,1]).list()
        [[{4}]]
        sage: OrderedMultisetPartitionsIntoSets(4, weight=[1,0,1]).list()
        [[{1}, {3}], [{1,3}], [{3}, {1}]]
        sage: OrderedMultisetPartitionsIntoSets(4, weight=[0,2]).list()
        [[{2}, {2}]]
        sage: OrderedMultisetPartitionsIntoSets(4, weight=[0,1,1]).list()
        []
        sage: OrderedMultisetPartitionsIntoSets(4, weight=[2,1]).list()
        [[{1}, {1}, {2}], [{1}, {1,2}], [{1}, {2}, {1}], [{1,2}, {1}], [{2}, {1}, {1}]]
        sage: O1 = OrderedMultisetPartitionsIntoSets(weight=[2,0,1])
        sage: O2 = OrderedMultisetPartitionsIntoSets(weight={1:2, 3:1})
        sage: O1 == O2
        True
        sage: OrderedMultisetPartitionsIntoSets(4, weight=[4]).list()
        [[{1}, {1}, {1}, {1}]]

    *The* ``size`` *keyword:*

    This is used to constrain the sum of entries across all blocks of the ordered
    multiset partition into sets. (This size is not pre-determined when alphabet
    `A` and order `d` are passed as required arguments.) For example, the ordered
    multiset partitions into sets of order 3 over the alphabet `[1,2,4]` that have
    size equal to 5 are as follows::

        sage: OMPs = OrderedMultisetPartitionsIntoSets
        sage: OMPs([1,2,4], 3, size=5).list()
        [[{1,2}, {2}], [{2}, {1,2}], [{2}, {2}, {1}],
         [{2}, {1}, {2}], [{1}, {2}, {2}]]

    *The* ``alphabet`` *option:*

    This is used to constrain which integers appear across all blocks of the
    ordered multiset partition into sets. For example, the ordered multiset
    partitions into sets of integer 4 are listed for different choices of alphabet
    below. Note that ``alphabet`` is allowed to be an integer or an iterable::

        sage: OMPs = OrderedMultisetPartitionsIntoSets
        sage: OMPs(4, alphabet=3).list()
        [[{1,3}], [{3}, {1}],
         [{1,2}, {1}], [{2}, {2}],
         [{2}, {1}, {1}], [{1}, {3}],
         [{1}, {1,2}], [{1}, {2}, {1}],
         [{1}, {1}, {2}], [{1}, {1}, {1}, {1}]]
        sage: OMPs(4, alphabet=3) == OMPs(4, alphabet=[1,2,3])
        True
        sage: OMPs(4, alphabet=[3]).list()
        []
        sage: OMPs(4, alphabet=[1,3]).list()
        [[{1,3}], [{3}, {1}], [{1}, {3}], [{1}, {1}, {1}, {1}]]
        sage: OMPs(4, alphabet=[2]).list()
        [[{2}, {2}]]
        sage: OMPs(4, alphabet=[1,2]).list()
        [[{1,2}, {1}], [{2}, {2}], [{2}, {1}, {1}], [{1}, {1,2}],
         [{1}, {2}, {1}], [{1}, {1}, {2}], [{1}, {1}, {1}, {1}]]
        sage: OMPs(4, alphabet=4).list() == OMPs(4).list()
        True

    *The* ``length``, ``min_length``, *and* ``max_length`` *options:*

    These are used to constrain the number of blocks within the ordered multiset
    partitions into sets. For example, the ordered multiset partitions into sets
    of integer 4 of length exactly 2, at least 2, and at most 2 are given by::

        sage: OrderedMultisetPartitionsIntoSets(4, length=2).list()
        [[{3}, {1}], [{1,2}, {1}], [{2}, {2}], [{1}, {3}], [{1}, {1,2}]]
        sage: OrderedMultisetPartitionsIntoSets(4, min_length=3).list()
        [[{2}, {1}, {1}], [{1}, {2}, {1}], [{1}, {1}, {2}], [{1}, {1}, {1}, {1}]]
        sage: OrderedMultisetPartitionsIntoSets(4, max_length=2).list()
        [[{4}], [{1,3}], [{3}, {1}], [{1,2}, {1}], [{2}, {2}], [{1}, {3}],
         [{1}, {1,2}]]

    *The* ``order``, ``min_order``, *and* ``max_order`` *options:*

    These are used to constrain the number of elements across all blocks of the
    ordered multiset partitions into sets. For example, the ordered multiset
    partitions into sets of integer 4 are listed by order below::

        sage: OrderedMultisetPartitionsIntoSets(4, order=1).list()
        [[{4}]]
        sage: OrderedMultisetPartitionsIntoSets(4, order=2).list()
        [[{1,3}], [{3}, {1}], [{2}, {2}], [{1}, {3}]]
        sage: OrderedMultisetPartitionsIntoSets(4, order=3).list()
        [[{1,2}, {1}], [{2}, {1}, {1}], [{1}, {1,2}], [{1}, {2}, {1}], [{1}, {1}, {2}]]
        sage: OrderedMultisetPartitionsIntoSets(4, order=4).list()
        [[{1}, {1}, {1}, {1}]]

    Also, here is a use of ``max_order``, giving the ordered multiset
    partitions into sets of integer 4 with order 1 or 2::

        sage: OrderedMultisetPartitionsIntoSets(4, max_order=2).list()
        [[{4}], [{1,3}], [{3}, {1}], [{2}, {2}], [{1}, {3}]]

    TESTS::

        sage: C = OrderedMultisetPartitionsIntoSets(8, length=3); C.cardinality()
        72
        sage: TestSuite(C).run()
    """
    @staticmethod
    def __classcall_private__(self, *args, **constraints):
        """
        Return the correct parent based upon the input:

        EXAMPLES::

            sage: OrderedMultisetPartitionsIntoSets()
            Ordered Multiset Partitions into Sets
            sage: OrderedMultisetPartitionsIntoSets(4)
            Ordered Multiset Partitions into Sets of integer 4
            sage: OrderedMultisetPartitionsIntoSets(4, max_order=2)
            Ordered Multiset Partitions into Sets of integer 4 with constraint: max_order=2

            sage: OrderedMultisetPartitionsIntoSets({1:2, 3:1})
            Ordered Multiset Partitions into Sets of multiset {{1, 1, 3}}
            sage: OrderedMultisetPartitionsIntoSets({1:2, 3:1}) == OrderedMultisetPartitionsIntoSets([1,1,3])
            True
            sage: OrderedMultisetPartitionsIntoSets({'a':2, 'c':1}, length=2)
            Ordered Multiset Partitions into Sets of multiset {{a, a, c}} with constraint: length=2
            sage: OrderedMultisetPartitionsIntoSets({'a':2, 'c':1}, length=4).list()
            []

            sage: OrderedMultisetPartitionsIntoSets(4, 3)
            Ordered Multiset Partitions into Sets of order 3 over alphabet {1, 2, 3, 4}
            sage: OrderedMultisetPartitionsIntoSets(['a', 'd'], 3)
            Ordered Multiset Partitions into Sets of order 3 over alphabet {a, d}
            sage: OrderedMultisetPartitionsIntoSets([2,4], 3, min_length=2)
            Ordered Multiset Partitions into Sets of order 3 over alphabet {2, 4}
             with constraint: min_length=2

        TESTS:

        The alphabet and order keywords cannot be used if they are also passed
        as required arguments, even if the values are compatible::

            sage: OrderedMultisetPartitionsIntoSets([1,2,4], 4, alphabet=[2,4], order=3)
            Traceback (most recent call last):
            ...
            ValueError: cannot pass alphabet as first argument and keyword argument
            sage: OrderedMultisetPartitionsIntoSets([1,2,4], 4, order=4)
            Traceback (most recent call last):
            ...
            ValueError: cannot pass order as second argument and keyword argument

        The weight, size, and order keywords cannot be used if a multiset is
        passed as a required argument, even if the values are compatible::

            sage: OrderedMultisetPartitionsIntoSets([1,1,4], weight={1:3, 2:1}).list()
            Traceback (most recent call last):
            ...
            ValueError: cannot pass multiset as first argument and weight as keyword argument
            sage: OrderedMultisetPartitionsIntoSets([1,1,4], size=6).list()
            Traceback (most recent call last):
            ...
            ValueError: cannot pass multiset as first argument and size as keyword argument
            sage: OrderedMultisetPartitionsIntoSets([1,1,4], weight={1:3, 2:1}, order=2).list()
            Traceback (most recent call last):
            ...
            ValueError: cannot pass multiset as first argument and ['order', 'weight'] as keyword arguments

        The size keyword cannot be used if it is also passed as a required argument,
        even if the value is compatible::

            sage: OrderedMultisetPartitionsIntoSets(5, size=5)
            Traceback (most recent call last):
            ...
            ValueError: cannot pass size as first argument and keyword argument
        """
        constraints = dict(constraints)
        if "weight" in constraints:
            # Should be a 'dictionary' of letter-frequencies, but accept a weak composition
            w = constraints["weight"]
            if not isinstance(w, dict):
                # make sure we didn't receive ``some_dict.items()``
                if len(w) > 0 and isinstance(w[0], (list, tuple)):
                    w = dict(w)
                else:
                    w = {i+1: w[i] for i in range(len(w)) if w[i] > 0}
            if not all((a in ZZ and a > 0) for a in w.values()):
                raise ValueError("%s must be a dictionary of letter-frequencies or a weak composition" % w)
            else:
                constraints["weight"] = tuple(w.items())

        if "alphabet" in constraints:
            A = constraints["alphabet"]
            if A in ZZ:
                A = range(1, A + 1)
            constraints["alphabet"] = frozenset(A)

        if len(args) == 2:  # treat as `alphabet` & `order`
            alph = args[0]
            order = args[1]
            if alph in ZZ:
                alph = range(1, alph + 1)
            if (alph and len(set(alph)) == len(alph)) and (order in ZZ and order >= 0):
                if "alphabet" in constraints:
                    raise ValueError("cannot pass alphabet as first argument and keyword argument")
                elif "order" in constraints:
                    raise ValueError("cannot pass order as second argument and keyword argument")
                if constraints == {}:
                    return OrderedMultisetPartitionsIntoSets_alph_d(frozenset(alph), order)
                else:
                    return OrderedMultisetPartitionsIntoSets_alph_d_constraints(frozenset(alph), order, **constraints)
            elif frozenset(alph) == frozenset() and order == 0:
                return OrderedMultisetPartitionsIntoSets_alph_d_constraints(frozenset(alph), order, **constraints)
            else:
                raise ValueError("alphabet=%s must be a nonempty set and order=%s must be a nonnegative integer" % (alph, order))

        elif len(args) == 1: # treat as `size` or `multiset`
            X = args[0]
            if isinstance(X, (list, tuple)):
                tmp = {}
                for i in X:
                    tmp[i] = tmp.get(i, 0) + 1
                X = tmp
            if isinstance(X, dict):
                over_determined = set(["size", "weight", "alphabet", "order", "min_order", "max_order"]).intersection(set(constraints))
                if over_determined:
                    if len(over_determined) > 1:
                        suff = "s"
                        offenses = str(sorted(over_determined))
                    else:
                        suff = ""
                        offenses = str(over_determined.pop())
                    raise ValueError("cannot pass multiset as first argument and %s as keyword argument%s" % (offenses, suff))
                X_items = tuple(X.items())
                if constraints == {}:
                    return OrderedMultisetPartitionsIntoSets_X(X_items)
                else:
                    return OrderedMultisetPartitionsIntoSets_X_constraints(X_items, **constraints)

            elif X in ZZ and X >= 0:
                if "size" in constraints:
                    raise ValueError("cannot pass size as first argument and keyword argument")
                if constraints == {}:
                    return OrderedMultisetPartitionsIntoSets_n(X)
                else:
                    return OrderedMultisetPartitionsIntoSets_n_constraints(X, **constraints)

            else:
                # zero arguments are passed?
                raise ValueError("%s must be a nonnegative integer or a list or dictionary representing a multiset" % X)

        elif len(args) > 2:
            raise ValueError("OrderedMultisetPartitonsIntoSets takes 1, 2, or 3 arguments")
        else:
            # try to do better than a generic parent
            if "weight" in constraints:
                X = constraints.pop("weight")
                return OrderedMultisetPartitionsIntoSets(dict(X), **constraints)
            elif "size" in constraints:
                n = constraints.pop("size")
                return OrderedMultisetPartitionsIntoSets(n, **constraints)
            elif "alphabet" in constraints and "order" in constraints:
                A = constraints.pop("alphabet")
                d = constraints.pop("order")
                return OrderedMultisetPartitionsIntoSets(A, d, **constraints)

            # generic parent
            return OrderedMultisetPartitionsIntoSets_all_constraints(**constraints)

    def __init__(self, is_finite=None, **constraints):
        """
        Initialize ``self``.

        TESTS::

            sage: c = {"length":4, "max_order":6, "alphabet":[2,4,5,6]}
            sage: OrderedMultisetPartitionsIntoSets(**c).constraints
            {'alphabet': frozenset({2, 4, 5, 6}), 'length': 4, 'max_order': 6}
            sage: OrderedMultisetPartitionsIntoSets(17, **c).constraints
            {'alphabet': frozenset({2, 4, 5, 6}), 'length': 4, 'max_order': 6}
            sage: OrderedMultisetPartitionsIntoSets(17, **c).full_constraints
            {'alphabet': frozenset({2, 4, 5, 6}), 'length': 4, 'max_order': 6, 'size': 17}

            sage: c = {"length":4, "min_length":5, "max_order":6, "order":5, "alphabet":4}
            sage: OrderedMultisetPartitionsIntoSets(**c).full_constraints
            {'alphabet': frozenset({1, 2, 3, 4}), 'length': 4, 'order': 5}
            sage: OrderedMultisetPartitionsIntoSets(**c).constraints
            {'length': 4}
            sage: OrderedMultisetPartitionsIntoSets(4, 5, **c).constraints
            Traceback (most recent call last):
            ...
            ValueError: cannot pass alphabet as first argument and keyword argument

            sage: c = {"weight":[2,2,0,3], "min_length":5, "max_order":6, "order":5, "alphabet":4}
            sage: OrderedMultisetPartitionsIntoSets(**c).constraints
            Traceback (most recent call last):
            ...
            ValueError: cannot pass multiset as first argument and ['alphabet', 'max_order', 'order'] as keyword arguments
        """
        constraints = dict(constraints)

        # standardize values for certain keywords
        if "alphabet" in constraints:
            if constraints["alphabet"] in ZZ:
                constraints["alphabet"] = frozenset(range(1, constraints["alphabet"]+1))
            else:
                constraints["alphabet"] = frozenset(constraints["alphabet"])

        if "weight" in constraints:
            X = dict(constraints["weight"])
            constraints["weight"] = X
            constraints.pop("alphabet", None)
            constraints.pop("min_order", None)
            constraints.pop("order", None)
            constraints.pop("max_order", None)
            constraints.pop("size", None)

        if "length" in constraints:
            constraints.pop("min_length", None)
            constraints.pop("max_length", None)
        min_k = constraints.get("min_length", 0)
        max_k = constraints.get("max_length", infinity)
        assert min_k <= max_k, "min_length=%s <= max_length=%s" % (min_k, max_k)
        if min_k == max_k:
            constraints["length"] = constraints.pop("min_length",
                                                    constraints.pop("max_length"))

        if "order" in constraints:
            constraints.pop("min_order", None)
            constraints.pop("max_order", None)
        min_ord = constraints.get("min_order", 0)
        max_ord = constraints.get("max_order", infinity)
        assert min_ord <= max_ord, "min_order=%s <= max_order=%s" % (min_ord, max_ord)
        if min_ord == max_ord:
            constraints["order"] = constraints.pop("min_order",
                                                   constraints.pop("max_order"))

        # pop keys with empty values, with the exception of 'size' or 'order'
        self.constraints = {}
        for (key,val) in constraints.items():
            if val:
                self.constraints[key] = val
            elif key in ("size", "order", "length") and val is not None:
                self.constraints[key] = val

        self.full_constraints = dict(self.constraints)
        if hasattr(self, "_X"):
            self.full_constraints["weight"] = dict(self._X)
            self.constraints.pop("weight", None)
        if hasattr(self, "_n"):
            self.full_constraints["size"] = self._n
            self.constraints.pop("size", None)
        if hasattr(self, "_alphabet"):
            self.full_constraints["alphabet"] = self._alphabet
            self.constraints.pop("alphabet", None)
            self.full_constraints["order"] = self._order
            self.constraints.pop("order", None)

        if is_finite or _is_finite(constraints):
            Parent.__init__(self, category=FiniteEnumeratedSets())
        else:
            Parent.__init__(self, category=InfiniteEnumeratedSets())

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: OrderedMultisetPartitionsIntoSets()
            Ordered Multiset Partitions into Sets
        """
        return "Ordered Multiset Partitions into Sets"

    def _constraint_repr_(self, cdict=None):
        """
        Return a string representation of all constraints
        appearing within ``self.constraints``.

        A helper method for ``self._repr_()``.

        EXAMPLES::

            sage: OMPs = OrderedMultisetPartitionsIntoSets()
            sage: c = {"length":4, "max_order":6, "alphabet":frozenset([2,4,5,6])}
            sage: OMPs._constraint_repr_(c)
            ' with constraints: alphabet={2, 4, 5, 6}, length=4, max_order=6'
            sage: c = {"size":14}
            sage: OMPs._constraint_repr_(c)
            ' with constraint: size=14'
        """
        if not cdict:
            cdict = dict(self.constraints)
        if "alphabet" in cdict:
            # make, e.g., `set([2,3,4])` print as `{2, 3, 4}`
            if not all(l in ZZ for l in cdict["alphabet"]):
                A = sorted(cdict["alphabet"], key=str)
            else:
                A = sorted(cdict["alphabet"])
            cdict["alphabet"] = "{" + repr(A)[1:-1] + "}"
        constr = ""
        ss = ['%s=%s' % item for item in cdict.items()]
        ss = sorted(ss)
        if len(ss) > 1:
            constr = " with constraints: " + ", ".join(ss)
        elif len(ss) == 1:
            constr = " with constraint: " + ", ".join(ss)
        return constr

    def _element_constructor_(self, lst):
        """
        Construct an element of ``self`` from ``lst``.

        EXAMPLES::

            sage: P = OrderedMultisetPartitionsIntoSets()
            sage: A = P([[3],[3,1]]) ; A # indirect doctest
            [{3}, {1,3}]
            sage: P1 = OrderedMultisetPartitionsIntoSets(7, alphabet=3)
            sage: A1 = P1([[3],[3,1]]); A1
            [{3}, {1,3}]
            sage: P2 = OrderedMultisetPartitionsIntoSets(alphabet=3)
            sage: A2 = P2([[3],[3,1]]); A2
            [{3}, {1,3}]
            sage: A == A1 == A2
            True
            sage: P = OrderedMultisetPartitionsIntoSets(3)
            sage: P([[3],[3,1]])
            Traceback (most recent call last):
            ...
            ValueError: cannot convert [[3], [3, 1]] into an element of
             Ordered Multiset Partitions into Sets of integer 3
        """
        if not lst:
            omp = []
        else:
            omp = [list(z) for z in lst]

        if omp in self:
            return self.element_class(self, list(map(frozenset, omp)))
        else:
            raise ValueError("cannot convert %s into an element of %s" % (lst, self))

    Element = OrderedMultisetPartitionIntoSets

    def __contains__(self, x):
        """
        Return if ``x`` is contained in ``self``.

        TESTS::

            sage: [[2,1], [1,3]] in OrderedMultisetPartitionsIntoSets()
            True
            sage: [[2,1], [1,3]] in OrderedMultisetPartitionsIntoSets(7)
            True
            sage: [[2,2], [1,3]] in OrderedMultisetPartitionsIntoSets()
            False
            sage: [] in OrderedMultisetPartitionsIntoSets()
            True
            sage: [] in OrderedMultisetPartitionsIntoSets(0)
            True
            sage: [] in OrderedMultisetPartitionsIntoSets(2)
            False
            sage: [[2, 1]] in OrderedMultisetPartitionsIntoSets(3, length=2)
            False
            sage: [[2, -1]] in OrderedMultisetPartitionsIntoSets()
            True
        """
        if not isinstance(x, (OrderedMultisetPartitionIntoSets, list, tuple)):
            return False
        return _has_nonempty_sets(x) and self._satisfies_constraints(x)

    def _satisfies_constraints(self, x):
        """
        Check whether or not ``x`` satisfies all of the constraints
        appearing within ``self.full_constraints`` (Boolean output).

        .. NOTE::

            This test will cause an infinite recursion with
            ``self._element_constructor_()`` if the ``__contains__``
            method in ``OrderedMultisetPartitionsIntoSets_X`` is removed.

        TESTS::

            sage: c = {"length":3, "max_order":5, "alphabet":[1,2,4], "size":12}
            sage: OMPs = OrderedMultisetPartitionsIntoSets(**c)
            sage: OMPs._satisfies_constraints([{2,4}, {1}, {1,4}])
            True
            sage: failures = {((2,4), (2,4)), ((1,2,4), (1,), (1,4)),
            ....:             ((2,4), (3,), (3,)), ((2,4), (1,), (2,4))}
            sage: any(OMPs._satisfies_constraints(x) for x in failures)
            False
            sage: c = {"max_length":4, "weight":{1:2, 2:1, 4:2}}
            sage: OMPs = OrderedMultisetPartitionsIntoSets(**c)
            sage: OMPs._satisfies_constraints([{2,4}, {1}, {1,4}])
            True
            sage: failures = {((2,), (4,), (1,), (1,), (4,)), ((1,), (1,), (2,4), (2,4))}
            sage: any(OMPs._satisfies_constraints(x) for x in failures)
            False
        """
        X = _concatenate(x)
        P = OrderedMultisetPartitionsIntoSets_X(tuple(_get_weight(X).items()))
        x = P.element_class(P, [frozenset(block) for block in x])
        constr = self.full_constraints
        tsts = []
        if 'size' in constr:
            tsts.append( x.size() == constr['size'] )
        if 'weight' in constr:
            tsts.append(  x.weight() == constr['weight'] )
        if 'alphabet' in constr:
            tsts.append(  frozenset(x.letters()).issubset(constr['alphabet']) )
        if 'length' in constr:
            tsts.append( x.length() == constr['length'] )
        if 'min_length' in constr:
            tsts.append(  x.length() >= constr['min_length'] )
        if 'max_length' in constr:
            tsts.append(  x.length() <= constr['max_length'] )
        if 'order' in constr:
            tsts.append( x.order() == constr['order'] )
        if 'min_order' in constr:
            tsts.append(  x.order() >= constr['min_order'] )
        if 'max_order' in constr:
            tsts.append(  x.order() <= constr['max_order'] )

        return all(tsts)

    def _from_list(self, lst):
        """
        Return an ordered multiset partition into sets of singleton blocks, whose
        singletons are the elements ``lst``.

        If any of the elements of ``lst`` are zero (or '0'), then use
        these as breaks points for the blocks.

        .. SEEALSO::

            :meth:`OrderedMultisetPartitionsIntoSets._from_list_with_zeros()`.

        INPUT:

        - ``lst`` -- an iterable

        EXAMPLES::

            sage: OMPs = OrderedMultisetPartitionsIntoSets()
            sage: OMPs._from_list([1,4,0,8])
            [{1,4}, {8}]
            sage: OMPs._from_list([1,4,8])
            [{1}, {4}, {8}]
            sage: OMPs._from_list([1,4,8,0]) == OrderedMultisetPartitionIntoSets([[1,4,8]])
            True
            sage: OMPs._from_list('abaa')
            [{'a'}, {'b'}, {'a'}, {'a'}]
            sage: OMPs._from_list('ab0a0a')
            [{'a','b'}, {'a'}, {'a'}]

        TESTS::

            sage: OMPs._from_list([1,0,2,3,1]) == OrderedMultisetPartitionIntoSets([[1], [2,3,1]])
            True
            sage: OMPs._from_list([1,2,'3',0,1]) == OrderedMultisetPartitionIntoSets([{1,2,'3'}, [1]])
            True
        """
        if all(a in ZZ for a in lst) and any(a < 0 for a in lst):
            raise ValueError("`_from_list` does not expect to see negative integers; received {}".format(str(lst)))
        if 0 in list(lst) or '0' in list(lst):
            return self._from_list_with_zeros(lst)

        d = [frozenset([x]) for x in lst]
        c = self.element_class(self, d)
        # give a better parent, if self is generic
        if isinstance(self, OrderedMultisetPartitionsIntoSets_all_constraints):
            P = OrderedMultisetPartitionsIntoSets(_get_weight(lst))
            return P.element_class(P, c)
        else:
            return self.element_class(self, c)

    def _from_list_with_zeros(self, lst_with_zeros):
        r"""
        Return an ordered multiset partition into sets from a list of nonnegative
        integers (or their string equivalents).

        Blocks are separated by zeros. Consecutive zeros are ignored.

        EXAMPLES::

            sage: OrderedMultisetPartitionsIntoSets()._from_list([1,2,4])
            [{1}, {2}, {4}]
            sage: OrderedMultisetPartitionsIntoSets()._from_list_with_zeros([1,2,4])
            [{1,2,4}]
            sage: OrderedMultisetPartitionsIntoSets()._from_list_with_zeros([1,0,2,0,0,4])
            [{1}, {2}, {4}]
            sage: OrderedMultisetPartitionsIntoSets()._from_list_with_zeros('abc00a0b')
            [{'a','b','c'}, {'a'}, {'b'}]
        """
        from_zero_lst = list(lst_with_zeros)
        if from_zero_lst[-1] not in {0, '0'}:
            from_zero_lst += [0]
        co = []
        block = []
        for a in from_zero_lst:
            if a in {0, '0'}:
                if block:
                    co.append(block)
                    block = []
            else:
                block.append(a)
        if co in self:
            c = self.element_class(self, map(frozenset, co))
            # give a better parent, if `self` is generic
            if isinstance(self, OrderedMultisetPartitionsIntoSets_all_constraints):
                P = OrderedMultisetPartitionsIntoSets(c.weight())
                return P.element_class(P, c)
            else:
                return c
        else:
            raise ValueError("ordered multiset partitions into sets do not have repeated entries within blocks (%s received)" % str(co))

    def __iter__(self):
        """
        Iterate over ordered multiset partitions into sets.

        EXAMPLES::

            sage: OrderedMultisetPartitionsIntoSets(3).list()
            [[{3}], [{1,2}], [{2}, {1}], [{1}, {2}], [{1}, {1}, {1}]]
            sage: OrderedMultisetPartitionsIntoSets(0).list()
            [[]]
            sage: C = OrderedMultisetPartitionsIntoSets()
            sage: it = C.__iter__()
            sage: [next(it) for i in range(16)]
            [[], [{1}], [{2}], [{1}, {1}], [{3}], [{1,2}], [{2}, {1}],
             [{1}, {2}], [{1}, {1}, {1}], [{4}], [{1,3}], [{3}, {1}],
             [{1,2}, {1}], [{2}, {2}], [{2}, {1}, {1}], [{1}, {3}]]

        TESTS::

            sage: OrderedMultisetPartitionsIntoSets(alphabet=[1,3], max_length=2).list()
            [[], [{1}], [{3}], [{1,3}], [{1}, {1}], [{1}, {3}],
             [{3}, {1}], [{3}, {3}], [{1,3}, {1}], [{1,3}, {3}],
             [{1}, {1,3}], [{3}, {1,3}], [{1,3}, {1,3}]]
            sage: C = OrderedMultisetPartitionsIntoSets(min_length=2, max_order=2)
            sage: it = C.__iter__()
            sage: [next(it) for i in range(15)]
            [[{1}, {1}], [{2}, {1}], [{1}, {2}], [{3}, {1}], [{2}, {2}],
             [{1}, {3}], [{4}, {1}], [{3}, {2}], [{2}, {3}], [{1}, {4}],
             [{5}, {1}], [{4}, {2}], [{3}, {3}], [{2}, {4}], [{1}, {5}]]
            sage: OrderedMultisetPartitionsIntoSets(alphabet=[1,3], min_length=2).list()
            Traceback (most recent call last):
            ...
            NotImplementedError: cannot list an infinite set
        """
        # Look for evidence of ``FiniteEnumeratedSets()`` among constraints.
        # ``_base_iterator`` ignores most constraints in ``self.full_constraints``.
        iterator = _base_iterator(self.full_constraints)
        if iterator:
            for co in iterator:
                if self._satisfies_constraints(co):
                    yield self.element_class(self, co)
        else:
            # iterate over blocks of letters over an alphabet
            if "alphabet" in self.constraints:
                A = self.constraints["alphabet"]
                # establish a cutoff order `max_ell`
                max = self.constraints.get("max_length", infinity)
                max = self.constraints.get("length", max)
                max = max * len(A)
                max = self.constraints.get("max_order", max)
                max_ell = self.constraints.get("order", max)
                ell = 0
                while True and ell <= max_ell:
                    for co in _iterator_order(A, ell):
                        if self._satisfies_constraints(co):
                            yield self.element_class(self, co)
                    ell += 1
            # or iterate over partitions of multisets of positive integers
            else:
                n = 0
                while True:
                    for co in _iterator_size(n):
                        if self._satisfies_constraints(co):
                            yield self.element_class(self, co)
                    n += 1

    def subset(self, size):
        """
        Return a subset of all ordered multiset partitions into sets.

        INPUT:

        - ``size`` -- integer representing a slice of all ordered
          multiset partitions into sets

        The slice alluded to above is taken with respect to length, or
        to order, or to size, depending on the constraints of  ``self``.

        EXAMPLES::

            sage: C = OrderedMultisetPartitionsIntoSets(weight={2:2, 3:1, 5:1})
            sage: C.subset(3)
            Ordered Multiset Partitions into Sets of multiset {{2, 2, 3, 5}} with constraint: length=3
            sage: C = OrderedMultisetPartitionsIntoSets(weight={2:2, 3:1, 5:1}, min_length=2)
            sage: C.subset(3)
            Ordered Multiset Partitions into Sets of multiset {{2, 2, 3, 5}} with constraint: length=3
            sage: C = OrderedMultisetPartitionsIntoSets(alphabet=[2,3,5])
            sage: C.subset(3)
            Ordered Multiset Partitions into Sets of order 3 over alphabet {2, 3, 5}
            sage: C = OrderedMultisetPartitionsIntoSets(order=5)
            sage: C.subset(3)
            Ordered Multiset Partitions into Sets of integer 3 with constraint: order=5
            sage: C = OrderedMultisetPartitionsIntoSets(alphabet=[2,3,5], order=5, length=3)
            sage: C.subset(3)
            Ordered Multiset Partitions into Sets of order 3 over alphabet {2, 3, 5} with constraint: length=3
            sage: C = OrderedMultisetPartitionsIntoSets()
            sage: C.subset(3)
            Ordered Multiset Partitions into Sets of integer 3
            sage: C.subset(3) == OrderedMultisetPartitionsIntoSets(3)
            True
        """
        fc = self.full_constraints

        # slice by 'length'
        if "weight" in fc:
            return OrderedMultisetPartitionsIntoSets(fc["weight"], length=size, **self.constraints)
        elif "alphabet" in fc and "size" in fc:
            add_length = dict(self.constraints)
            add_length["length"] = size
            return OrderedMultisetPartitionsIntoSets(fc["alphabet"], fc["order"], **add_length)

        # slice by 'order'
        if "alphabet" in fc:
            no_alpha = {k: v for (k, v) in self.constraints.items() if k != "alphabet"}
            return OrderedMultisetPartitionsIntoSets(fc["alphabet"], size, **no_alpha)

        # slice by 'size'
        return OrderedMultisetPartitionsIntoSets(size, **self.constraints)

###############


class OrderedMultisetPartitionsIntoSets_all_constraints(OrderedMultisetPartitionsIntoSets):
    r"""
    All ordered multiset partitions into sets (with or without constraints).

    EXAMPLES::

        sage: C = OrderedMultisetPartitionsIntoSets(); C
        Ordered Multiset Partitions into Sets
        sage: [[1],[1,'a']] in C
        True

        sage: OrderedMultisetPartitionsIntoSets(weight=[2,0,1], length=2)
        Ordered Multiset Partitions into Sets of multiset {{1, 1, 3}} with constraint: length=2

    TESTS::

        sage: OMP = OrderedMultisetPartitionsIntoSets()
        sage: TestSuite(OMP).run()  # long time

        sage: C = OrderedMultisetPartitionsIntoSets(weight=[2,0,1], length=2)
        sage: TestSuite(C).run()

        sage: D1 = OrderedMultisetPartitionsIntoSets(weight={1:2, 3:1}, min_length=2, max_length=2)
        sage: D2 = OrderedMultisetPartitionsIntoSets({1:2, 3:1}, min_length=2, max_length=2)
        sage: D3 = OrderedMultisetPartitionsIntoSets(5, weight={1:2, 3:1}, length=2)
        sage: D4 = OrderedMultisetPartitionsIntoSets([1,3], 3, weight={1:2, 3:1}, length=2)
        sage: D5 = OrderedMultisetPartitionsIntoSets([1,3], 3, size=5, length=2)
        sage: all(C != D for D in [D1, D2, D3, D4, D5])
        True
        sage: all(Set(C) == Set(D) for D in [D1, D2, D3, D4, D5])
        True
        sage: E = OrderedMultisetPartitionsIntoSets({1:2, 3:1}, min_length=2)
        sage: Set(C) == Set(E)
        False
    """

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: OrderedMultisetPartitionsIntoSets(min_length=3, max_order=5)
            Ordered Multiset Partitions into Sets with constraints: max_order=5, min_length=3
            sage: OrderedMultisetPartitionsIntoSets(min_length=3, max_order=5, alphabet=[1,'a'])
            Ordered Multiset Partitions into Sets with constraints:
             alphabet={1, 'a'}, max_order=5, min_length=3
        """
        return "Ordered Multiset Partitions into Sets" + self._constraint_repr_()

###############


class OrderedMultisetPartitionsIntoSets_n(OrderedMultisetPartitionsIntoSets):
    """
    Ordered multiset partitions into sets of a fixed integer `n`.
    """

    def __init__(self, n):
        """
        Initialize ``self``.

        TESTS::

            sage: C = OrderedMultisetPartitionsIntoSets(Integer(4))
            sage: TestSuite(C).run()
            sage: C2 = OrderedMultisetPartitionsIntoSets(int(4))
            sage: C is C2
            True
            sage: C3 = OrderedMultisetPartitionsIntoSets(7/2)
            Traceback (most recent call last):
            ...
            ValueError:  7/2 must be a nonnegative integer or a list or
             dictionary representing a multiset
        """
        self._n = n
        OrderedMultisetPartitionsIntoSets.__init__(self, True)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: OrderedMultisetPartitionsIntoSets(3)
            Ordered Multiset Partitions into Sets of integer 3
        """
        return "Ordered Multiset Partitions into Sets of integer %s" % self._n

    def cardinality(self):
        """
        Return the number of elements in ``self``.

        TESTS::

            sage: len(OrderedMultisetPartitionsIntoSets(10).list())
            1500
            sage: OrderedMultisetPartitionsIntoSets(10).cardinality()
            1500
        """
        # Dispense with the complex computation for small orders.
        if self._n <= 5:
            orders = {0: 1, 1: 1, 2: 2, 3: 5, 4: 11, 5: 25}
            return ZZ(orders[self._n])

        # We view an ordered multiset partition into sets as a list of 2-regular integer partitions.
        #
        # The 2-regular partitions have a nice generating function (see OEIS:A000009).
        # Below, we take (products of) coefficients of polynomials to compute cardinality.
        t = PowerSeriesRing(ZZ, 't').gen().O(self._n + 1)
        partspoly = prod(1 + t**k for k in range(1, self._n + 1)).dict()
        deg = 0
        for alpha in composition_iterator_fast(self._n):
            deg += prod(partspoly[d] for d in alpha)
        return ZZ(deg)

    def _an_element_(self):
        """
        Return a typical element of ``self``.

        EXAMPLES::

            sage: OrderedMultisetPartitionsIntoSets(13).an_element()
            [{2,3}, {2,3}, {1,2}]
            sage: OrderedMultisetPartitionsIntoSets(14).an_element()
            [{2,3}, {2,3}, {4}]
        """
        #output will have at most three blocks, each of size 1, 2, or 3.
        alpha = Compositions(self._n, max_part=self._n//3+1).an_element()
        out = []
        for a in alpha:
            if a in {1, 2, 4}:
                out.append([a])
            else:
                if a % 2:
                    out.append([a//2+1, a//2])
                else:
                    out.append([a//2, a//2-1, 1])
        return self.element_class(self, map(frozenset, out))

    def random_element(self):
        """
        Return a random element of ``self``.

        This method does not return elements of ``self`` with uniform probability,
        but it does cover all elements. The scheme is as follows:

        - produce a random composition `C`;
        - choose a random partition of `c` into distinct parts for each `c` in `C`.

        EXAMPLES::

            sage: OrderedMultisetPartitionsIntoSets(5).random_element()  # random
            [{1,2}, {1}, {1}]
            sage: OrderedMultisetPartitionsIntoSets(5).random_element()  # random
            [{2}, {1,2}]

            sage: OMP = OrderedMultisetPartitionsIntoSets(5)
            sage: d = {}
            sage: for _ in range(1100):
            ....:     x = OMP.random_element()
            ....:     d[x] = d.get(x, 0) + 1
            sage: d.values()  # random
            [72, 73, 162, 78, 135, 75, 109, 65, 135, 134, 62]
        """
        C = Compositions(self._n).random_element()
        co = [IntegerListsLex(c, min_part=1, max_part=c,
                              min_slope=1).random_element() for c in C]
        return self.element_class(self, map(frozenset, co))

    def __iter__(self):
        """
        Iterate over ``self``.

        TESTS::

            sage: O = OrderedMultisetPartitionsIntoSets(6)
            sage: it = O.__iter__()
            sage: [next(it) for _ in range(10)]
            [[{6}], [{2,4}], [{1,5}], [{1,2,3}],
             [{5}, {1}], [{2,3}, {1}], [{1,4}, {1}],
             [{4}, {2}], [{1,3}, {2}], [{4}, {1}, {1}]]
        """
        for co in _iterator_size(self._n):
            yield self.element_class(self, co)


class OrderedMultisetPartitionsIntoSets_n_constraints(OrderedMultisetPartitionsIntoSets):
    """
    Class of ordered multiset partitions into sets of a fixed integer `n`
    satisfying constraints.
    """

    def __init__(self, n, **constraints):
        """
        Mimic class ``OrderedMultisetPartitionsIntoSets_n`` to initialize.

        TESTS::

            sage: C = OrderedMultisetPartitionsIntoSets(6, length=3)
            sage: TestSuite(C).run()

            sage: C = OrderedMultisetPartitionsIntoSets(6, weight=[3,0,1], length=3)
            sage: TestSuite(C).run()
        """
        self._n = n
        OrderedMultisetPartitionsIntoSets.__init__(self, True, size=n, **constraints)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: O = OrderedMultisetPartitionsIntoSets(14, length=4, max_order=6, alphabet={2,4,5,6})
            sage: O
            Ordered Multiset Partitions into Sets of integer 14 with constraints:
             alphabet={2, 4, 5, 6}, length=4, max_order=6
        """
        cdict = dict(self.constraints)
        cdict.pop("size", None)
        base_repr = "Ordered Multiset Partitions into Sets of integer %s" % self._n
        return base_repr + self._constraint_repr_(cdict)

###############


class OrderedMultisetPartitionsIntoSets_X(OrderedMultisetPartitionsIntoSets):
    """
    Class of ordered multiset partitions into sets of a fixed multiset `X`.
    """

    def __init__(self, X):
        """
        Initialize ``self``.

        TESTS::

            sage: C = OrderedMultisetPartitionsIntoSets([1,1,4])
            sage: TestSuite(C).run()

            sage: C2 = OrderedMultisetPartitionsIntoSets({1:2, 4:1})
            sage: C is C2
            True
        """
        self._X = X
        # sort the multiset
        if all((k in ZZ and k > 0) for (k,v) in X):
            self._Xtup = tuple([k for (k,v) in sorted(X) for _ in range(v)])
        else:
            self._Xtup = tuple([k for (k,v) in sorted(X, key=str) for _ in range(v)])
        OrderedMultisetPartitionsIntoSets.__init__(self, True)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: repr(OrderedMultisetPartitionsIntoSets([1,1,4]))
            'Ordered Multiset Partitions into Sets of multiset {{1, 1, 4}}'
        """
        ms_rep = "{{" + ", ".join(map(str, self._Xtup)) + "}}"
        return "Ordered Multiset Partitions into Sets" + " of multiset %s" % ms_rep

    def __contains__(self, x):
        """
        Return if ``x`` is contained in ``self``.

        TESTS::

            sage: from sage.combinat.multiset_partition_into_sets_ordered import OrderedMultisetPartitionsIntoSets_X as OMPX
            sage: [[2,1], [1,3]] in OMPX(((1,2), (2,1), (3,1)))
            True
            sage: co = OrderedMultisetPartitionIntoSets([[2,1], [1,3]])
            sage: co in OMPX(((1,2), (2,1), (3,1)))
            True
            sage: [[2,1], [2,3]] in OMPX(((1,2), (2,1), (3,1)))
            False
            sage: [] in OMPX(())
            True
            sage: [[2, -1], [2,'a']] in OMPX(((2,2), (-1,1), ('a',1)))
            True
        """
        if not isinstance(x, (OrderedMultisetPartitionIntoSets, list, tuple)):
            return False

        x_Xtup = sorted(_concatenate(x), key=str)
        self_Xtup = sorted(self._Xtup, key=str)
        return _has_nonempty_sets(x) and x_Xtup == self_Xtup

    def cardinality(self):
        """
        Return the number of ordered partitions of multiset ``X``.

        TESTS::

            sage: len(OrderedMultisetPartitionsIntoSets([2,2,2,3,4,5]).list())
            535
            sage: OrderedMultisetPartitionsIntoSets([2,2,2,3,4,5]).cardinality()
            535
        """
        if self._Xtup == ():
            return ZZ(0)

        # We build ordered multiset partitions into sets of `X` by permutation + deconcatenation
        deg = 0
        for alpha in Permutations_mset(self._Xtup):
            fattest = _break_at_descents(alpha)
            deg += prod(2**(len(k)-1) for k in fattest)
        return ZZ(deg)

    def _an_element_(self):
        """
        Return a typical element of ``self``.

        EXAMPLES::

            sage: OrderedMultisetPartitionsIntoSets([2,2,2,3,4,5]).an_element()
            [{2}, {2}, {2,3,4}, {5}]
            sage: OrderedMultisetPartitionsIntoSets([2,2,2,3,4,4,5]).an_element()
            [{2}, {2}, {2,3}, {4}, {4,5}]
        """
        if not self._Xtup:
            return self.element_class(self, [])
        alpha = Permutations_mset(self._Xtup).an_element()
        co = _break_at_descents(alpha)

        # construct "an element" by breaking the first fat block of `co` in two
        elt = []
        for i in range(len(co)):
            if len(co[i]) == 1:
                elt.append(co[i])
            else:
                break
        elt.append(co[i][:len(co[i])//2 + 1])
        elt.append(co[i][len(co[i])//2 + 1:])
        elt.extend(co[i+1:])
        return self.element_class(self, map(frozenset, elt))

    def random_element(self):
        """
        Return a random element of ``self``.

        This method does not return elements of ``self`` with uniform probability,
        but it does cover all elements. The scheme is as follows:

        - produce a random permutation ``p`` of the multiset;
        - create blocks of an OMP ``fat`` by breaking ``p`` after non-ascents;
        - take a random element of ``fat.finer()``.

        EXAMPLES::

            sage: OrderedMultisetPartitionsIntoSets([1,1,3]).random_element()  # random
            [{1}, {1,3}]
            sage: OrderedMultisetPartitionsIntoSets([1,1,3]).random_element()  # random
            [{3}, {1}, {1}]

            sage: OMP = OrderedMultisetPartitionsIntoSets([1,1,3,3])
            sage: d = {}
            sage: for _ in range(1000):
            ....:     x = OMP.random_element()
            ....:     d[x] = d.get(x, 0) + 1
            sage: d.values()  # random
            [102, 25, 76, 24, 66, 88, 327, 27, 83, 83, 239, 72, 88]
        """
        if not self._Xtup:
            return self.element_class(self, [])

        alpha = Permutations_mset(self._Xtup).random_element()
        co = _break_at_descents(alpha)
        finer = self.element_class(self, map(frozenset,co)).finer()
        return FiniteEnumeratedSets()(finer).random_element()

    def __iter__(self):
        """
        Iterate over ``self``.

        TESTS::

            sage: O = OrderedMultisetPartitionsIntoSets(['a', 'b', 'a'])
            sage: sorted(O, key=str)
            [[{'a','b'}, {'a'}],
             [{'a'}, {'a','b'}],
             [{'a'}, {'a'}, {'b'}],
             [{'a'}, {'b'}, {'a'}],
             [{'b'}, {'a'}, {'a'}]]

            sage: O = OrderedMultisetPartitionsIntoSets([1, 1, 2])
            sage: list(O)
            [[{1}, {1}, {2}], [{1}, {1,2}], [{1}, {2}, {1}],
            [{1,2}, {1}], [{2}, {1}, {1}]]
        """
        for co in _iterator_weight(weight=dict(self._X)):
            yield self.element_class(self, co)


class OrderedMultisetPartitionsIntoSets_X_constraints(OrderedMultisetPartitionsIntoSets):
    """
    Class of ordered multiset partitions into sets of a fixed multiset `X`
    satisfying constraints.
    """

    def __init__(self, X, **constraints):
        """
        Mimic class ``OrderedMultisetPartitionsIntoSets_X`` to initialize.

        TESTS::

            sage: C = OrderedMultisetPartitionsIntoSets([1,1,2,4], length=3)
            sage: TestSuite(C).run()

            sage: C = OrderedMultisetPartitionsIntoSets([1,1,2,4], max_length=3)
            sage: TestSuite(C).run()
        """
        self._X = X
        self._Xtup = tuple(k for (k,v) in sorted(X) for _ in range(v))
        OrderedMultisetPartitionsIntoSets.__init__(self, True, weight=X, **constraints)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: O = OrderedMultisetPartitionsIntoSets([2,2,2,3,4,4,5], min_length=4, max_length=5)
            sage: O
            Ordered Multiset Partitions into Sets of multiset {{2, 2, 2, 3, 4, 4, 5}}
             with constraints: max_length=5, min_length=4
        """
        cdict = dict(self.constraints)
        cdict.pop("weight", None)
        ms_rep = "{{" + ", ".join(map(str, self._Xtup)) + "}}"
        base_repr = "Ordered Multiset Partitions into Sets" + " of multiset %s" % ms_rep
        return base_repr + self._constraint_repr_(cdict)

###############


class OrderedMultisetPartitionsIntoSets_alph_d(OrderedMultisetPartitionsIntoSets):
    """
    Class of ordered multiset partitions into sets of specified order `d`
    over a fixed alphabet `A`.
    """

    def __init__(self, A, d):
        """
        Initialize ``self``.

        TESTS::

            sage: C = OrderedMultisetPartitionsIntoSets(3, 2)
            sage: TestSuite(C).run()

            sage: C2 = OrderedMultisetPartitionsIntoSets([1,2,3], 2)
            sage: C is C2
            True

            sage: list(OrderedMultisetPartitionsIntoSets([1,2,3], 2))
            [[{1,2}], [{1,3}], [{2,3}], [{1}, {1}], [{1}, {2}], [{1}, {3}], [{2}, {1}],
             [{2}, {2}], [{2}, {3}], [{3}, {1}], [{3}, {2}], [{3}, {3}]]
        """
        self._alphabet = A
        self._order = d
        OrderedMultisetPartitionsIntoSets.__init__(self, True)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: repr(OrderedMultisetPartitionsIntoSets(3, 2))
            'Ordered Multiset Partitions into Sets of order 2 over alphabet {1, 2, 3}'
            sage: repr(OrderedMultisetPartitionsIntoSets([1,3], 2))
            'Ordered Multiset Partitions into Sets of order 2 over alphabet {1, 3}'
        """
        A_rep = "Ordered Multiset Partitions into Sets of order " + str(self._order)
        A_rep += " over alphabet {%s}" % (", ".join(map(str, sorted(self._alphabet))))
        return A_rep

    def _an_element_(self):
        """
        Return a typical element of ``OrderedMultisetPartitionIntoSets_alph_d``.

        EXAMPLES::

            sage: OrderedMultisetPartitionsIntoSets([2,3,4,5], 3).an_element()
            [{2,4,5}]
        """
        alpha = Compositions(self._order, max_part=len(self._alphabet)).an_element()
        co = [Subsets_sk(self._alphabet, a).an_element() for a in alpha]
        return self.element_class(self, map(frozenset, co))

    def random_element(self):
        r"""
        Return a random element of ``self``.

        This method does not return elements of ``self`` with uniform probability,
        but it does cover all elements. The scheme is as follows:

        - produce a random composition `C`;
        - choose random subsets of ``self._alphabet`` of size `c` for each `c` in `C`.

        EXAMPLES::

            sage: OrderedMultisetPartitionsIntoSets([1,4], 3).random_element()  # random
            [{4}, {1,4}]
            sage: OrderedMultisetPartitionsIntoSets([1,3], 4).random_element()  # random
            [{1,3}, {1}, {3}]

            sage: OMP = OrderedMultisetPartitionsIntoSets([2,3,4], 2)
            sage: d = {}
            sage: for _ in range(1200):
            ....:     x = OMP.random_element()
            ....:     d[x] = d.get(x, 0) + 1
            sage: d.values()  # random
            [192, 68, 73, 61, 69, 60, 77, 204, 210, 66, 53, 67]
        """
        if not self._alphabet:
            return self.element_class(self, [])

        alpha = Compositions(self._order, max_part=len(self._alphabet)).random_element()
        co = [Subsets_sk(self._alphabet, a).random_element() for a in alpha]
        return self.element_class(self, map(frozenset, co))

    def __iter__(self):
        """
        Iterate over ``self``.

        TESTS::

            sage: O = OrderedMultisetPartitionsIntoSets(['a', 'b'], 3)
            sage: it = O.__iter__()
            sage: sorted([next(it) for _ in range(O.cardinality())], key=str)
            [[{'a','b'}, {'a'}], [{'a','b'}, {'b'}], [{'a'}, {'a','b'}],
             [{'a'}, {'a'}, {'a'}], [{'a'}, {'a'}, {'b'}], [{'a'}, {'b'}, {'a'}],
             [{'a'}, {'b'}, {'b'}], [{'b'}, {'a','b'}], [{'b'}, {'a'}, {'a'}],
             [{'b'}, {'a'}, {'b'}], [{'b'}, {'b'}, {'a'}], [{'b'}, {'b'}, {'b'}]]
        """
        for co in _iterator_order(self._alphabet, self._order):
            yield self.element_class(self, co)

    def cardinality(self):
        """
        Return the number of ordered partitions of order ``self._order`` on
        alphabet ``self._alphabet``.

        TESTS::

            sage: len(OrderedMultisetPartitionsIntoSets([1, 'a'], 3).list())
            12
            sage: OrderedMultisetPartitionsIntoSets([1, 'a'], 3).cardinality()
            12
        """
        if self._order == 0:
            return ZZ(0)

        # iteration scheme:
        # - start from an integer composition ``alpha`` of ``self._order``.
        # - for each ``a`` in ``alpha``, pick ``a`` letters from ``alphabet``
        min_length = self._order // len(self._alphabet)
        max_length = self._order

        deg = 0
        for k in range(min_length, max_length+1):
            for alpha in IntegerListsLex(self._order, length=k, min_part=1, max_part=len(self._alphabet)):
                deg += prod(binomial(len(self._alphabet), a) for a in alpha)
        return ZZ(deg)


class OrderedMultisetPartitionsIntoSets_alph_d_constraints(OrderedMultisetPartitionsIntoSets):
    """
    Class of ordered multiset partitions into sets of specified order `d`
    over a fixed alphabet `A` satisfying constraints.
    """

    def __init__(self, A, d, **constraints):
        """
        Mimic class ``OrderedMultisetPartitionsIntoSets_alph_d`` to initialize.

        EXAMPLES::

            sage: list(OrderedMultisetPartitionsIntoSets(3, 2, length=3))
            []
            sage: list(OrderedMultisetPartitionsIntoSets([1,2,4], 2, length=1))
            [[{1,2}], [{1,4}], [{2,4}]]

        TESTS::

            sage: C = OrderedMultisetPartitionsIntoSets(3, 2, length=3)
            sage: TestSuite(C).run()

            sage: C = OrderedMultisetPartitionsIntoSets([1,2,4], 4, min_length=3)
            sage: TestSuite(C).run()
        """
        self._alphabet = A
        self._order = d
        OrderedMultisetPartitionsIntoSets.__init__(self, True, alphabet=A,
                                                   order=d, **constraints)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: O = OrderedMultisetPartitionsIntoSets([2,3,4,5], 4, length=3)
            sage: O
            Ordered Multiset Partitions into Sets of order 4 over alphabet {2, 3, 4, 5}
             with constraint: length=3
            sage: O = OrderedMultisetPartitionsIntoSets([2,3,4,5], 4, max_length=4, size=10)
            sage: O
            Ordered Multiset Partitions into Sets of order 4 over alphabet {2, 3, 4, 5}
             with constraints: max_length=4, size=10
        """
        cdict = dict(self.constraints)
        cdict.pop("alphabet", None)
        cdict.pop("order", None)
        base_repr = "Ordered Multiset Partitions into Sets of order " + str(self._order)
        base_repr += " over alphabet {%s}" % (", ".join(map(str, sorted(self._alphabet))))
        return base_repr + self._constraint_repr_(cdict)

###############


def _get_multiset(co):
    """
    Construct the multiset (as a sorted tuple) suggested by the lists
    of lists ``co``.

    EXAMPLES::

        sage: from sage.combinat.multiset_partition_into_sets_ordered import _get_multiset
        sage: L = ((1,), (1, 6), (6, 7), (1,), (1, 3))
        sage: _get_multiset(L)
        (1, 1, 1, 1, 3, 6, 6, 7)
    """
    return tuple(sorted(_concatenate(co), key=str))


def _get_weight(lst):
    """
    Construct the multiset (as a dictionary) suggested by the
    multiset-as-list ``lst``.

    EXAMPLES::

        sage: from sage.combinat.multiset_partition_into_sets_ordered import _get_weight
        sage: L = (1, 1, 1, 3, 1, 6, 6, 7)
        sage: _get_weight(L)
        {1: 4, 3: 1, 6: 2, 7: 1}
    """
    out = {}
    for k in lst:
        out[k] = out.get(k,0) + 1
    return out


def _has_nonempty_sets(x):
    """
    Blocks should be nonempty sets/lists/tuples of distinct elements.

    TESTS::

        sage: from sage.combinat.multiset_partition_into_sets_ordered import _has_nonempty_sets
        sage: _has_nonempty_sets([[2,4], {1}, (1,4)])
        True
        sage: _has_nonempty_sets([[2,4], {}, (1,4)])
        False
        sage: _has_nonempty_sets([(2,4), (1,1), (1,4)])
        False
    """
    return all((isinstance(block, (list, tuple, set, frozenset, Set_object))
                and block and len(set(block)) == len(block))
               for block in x)


def _union_of_sets(list_of_sets):
    """
    Return the union of a list of iterables as a frozenset.

    EXAMPLES::

        sage: from sage.combinat.multiset_partition_into_sets_ordered import _union_of_sets
        sage: L = ([1,2,3], Set([1,5,6]), [], range(5,8))
        sage: _union_of_sets(L)
        frozenset({1, 2, 3, 5, 6, 7})
    """
    return reduce(lambda a, b: frozenset(a) | frozenset(b),
                  list_of_sets, frozenset())


def _concatenate(list_of_iters):
    """
    Return the concatenation of a list of iterables as a tuple.

    EXAMPLES::

        sage: from sage.combinat.multiset_partition_into_sets_ordered import _concatenate
        sage: L = ([1,2,3], Set([4,5,6]), [], range(7,11))
        sage: _concatenate(L)
        (1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
    """
    return tuple([val for block in list_of_iters for val in block])


def _is_finite(constraints):
    """
    Return ``True`` if the dictionary ``constraints`` corresponds to
    a finite collection of ordered multiset partitions into sets.

    If either ``weight`` or ``size`` is among the constraints, then
    the constraints represent a finite collection of ordered multiset
    partitions into sets. If both are absent, one needs ``alphabet`` to be
    present (plus a bound on length or order) in order to have a
    finite collection of ordered multiset partitions into sets.

    EXAMPLES::

        sage: from sage.combinat.multiset_partition_into_sets_ordered import _is_finite
        sage: W = {"weight": {1:3, 2:3, 4:1}, "length": 5}
        sage: S = {"size": 44, "min_length": 5}
        sage: AO = {"alphabet": range(44), "max_order": 5}
        sage: all(_is_finite(constr) for constr in (W, S, AO))
        True
        sage: AL = {"alphabet": range(44), "min_order": 5}
        sage: _is_finite(AL)
        False
    """
    if "weight" in constraints or "size" in constraints:
        return True
    elif "alphabet" in constraints:
        # Assume the alphabet is finite
        Bounds = set(["length", "max_length", "order", "max_order"])
        return Bounds.intersection(set(constraints)) != set()


def _base_iterator(constraints):
    """
    Return a base iterator for ordered multiset partitions into sets or ``None``.

    If the keys within ``constraints`` dictionary correspond to a finite set
    of ordered multiset partitions into sets, return an iterator. Else,
    return ``None``.

    OUTPUT:

    Tuples of ``frozenset`` objects representing ordered multiset partitions
    into sets.

    EXAMPLES:

    If key ``weight`` is present, ignore all other constraints
    (passes to ``_iterator_weight``)::

        sage: from sage.combinat.multiset_partition_into_sets_ordered import _base_iterator
        sage: OMP = OrderedMultisetPartitionIntoSets
        sage: constraints = {"weight": {1:3, 2:3, 4:1}, "length": 5}
        sage: it = _base_iterator(constraints)
        sage: sorted(OMP(next(it)) for _ in range(4)) # note the partitions of length 6 and 7
        [[{1}, {1}, {1}, {2}, {2}, {2}, {4}],
         [{1}, {1}, {1}, {2}, {2}, {2,4}],
         [{1}, {1}, {1,2}, {2}, {2}, {4}],
         [{1}, {1}, {1,2}, {2}, {2,4}]]

    If key ``size`` is present, pass to ``_iterator_size``, which then
    takes into account whether or not keys ``length`` and ``alphabet``
    are among the constraints::

        sage: constraints = {"size": 5}
        sage: it = _base_iterator(constraints)
        sage: [OMP(next(it)) for _ in range(8)]
        [[{5}], [{2,3}], [{1,4}], [{4}, {1}], [{1,3}, {1}],
         [{3}, {2}], [{1,2}, {2}], [{3}, {1}, {1}]]

        sage: constraintsL = {"size": 6, "length":2}
        sage: it = _base_iterator(constraintsL)
        sage: [OMP(next(it)) for _ in range(8)]
        [[{5}, {1}], [{2,3}, {1}], [{1,4}, {1}], [{4}, {2}],
         [{1,3}, {2}], [{3}, {3}], [{3}, {1,2}], [{1,2}, {3}]]

        sage: constraintsA = {"size": 6, "alphabet":frozenset([2, 3])}
        sage: it = _base_iterator(constraintsA)
        sage: list(it)
        [(frozenset({3}), frozenset({3})),
         (frozenset({2}), frozenset({2}), frozenset({2}))]

    If key ``alphabet`` is present, the slice may still be infinite, in
    which case ``None`` is returned. Else, use to ``_iterator_order``::

        sage: constraints = {"alphabet": frozenset([3, 4]), "min_length":2}
        sage: _base_iterator(constraints)  is None
        True
        sage: constraints = {"alphabet": frozenset([3, 4]), "max_length":2}
        sage: it = _base_iterator(constraints)
        sage: list(map(OMP, it))
        [[], [{3}], [{4}], [{3,4}], [{3}, {3}], [{3}, {4}],
         [{4}, {3}], [{4}, {4}], [{3,4}, {3}], [{3,4}, {4}],
         [{3}, {3,4}], [{4}, {3,4}], [{3,4}, {3,4}]]
    """
    if "weight" in constraints:
        return _iterator_weight(constraints["weight"])
    elif "size" in constraints:
        return _iterator_size(constraints["size"],
            constraints.get("length",None), constraints.get("alphabet",None))
    elif "alphabet" in constraints:
        A = constraints["alphabet"]
        # assumes `alphabet` is finite
        min_k = constraints.get("min_length", 0)
        max_k = constraints.get("max_length", infinity)
        min_ord = constraints.get("min_order", 0)
        max_ord = constraints.get("max_order", max_k * len(A))
        max_k = min(max_k, max_ord)
        if "length" in constraints:
            min_k = max_k = constraints["length"]
            min_ord = max(min_ord, min_k)
            max_ord = min(max_ord, len(A) * max_k)
        if "order" in constraints:
            min_ord = max_ord = constraints["order"]
            max_k = min(max_k, max_ord)
            if min_ord:
                min_k = max(1, min_k, min_ord // len(A))
        if infinity not in (max_k, max_ord):
            return chain(*(_iterator_order(A, ord, range(min_k, max_k + 1))
                           for ord in range(min_ord, max_ord + 1)))
    # else
    return None


def _iterator_weight(weight):
    """
    An iterator for the ordered multiset partitions into sets with weight given by
    the dictionary (or weak composition) ``weight``.

    The dictionary ``weight`` may contain values equal to `0`;
    the corresponding keys are ignored.

    OUTPUT:

    Tuples of ``frozenset`` objects representing ordered multiset partitions
    into sets.

    EXAMPLES::

        sage: from sage.combinat.multiset_partition_into_sets_ordered import _iterator_weight
        sage: weight = {1:2, 'b':1}
        sage: OMP = OrderedMultisetPartitionsIntoSets(weight)
        sage: l = list(_iterator_weight(weight))

        sage: sorted(map(OMP, l), key=str) == sorted(map(OMP,
        ....:     [[{1}, {1}, {'b'}], [{1}, {1,'b'}], [{1}, {'b'}, {1}],
        ....:     [{1,'b'}, {1}], [{'b'}, {1}, {1}]]), key=str)
        True
        sage: OMP = OrderedMultisetPartitionsIntoSets({1:3, 3:1})
        sage: list(map(OMP, _iterator_weight([3,0,1])))
        [[{1}, {1}, {1}, {3}], [{1}, {1}, {1,3}], [{1}, {1}, {3}, {1}],
         [{1}, {1,3}, {1}], [{1}, {3}, {1}, {1}],
         [{1,3}, {1}, {1}], [{3}, {1}, {1}, {1}]]

    TESTS::

        sage: list(_iterator_weight([2,0,1]))
        [(frozenset({1}), frozenset({1}), frozenset({3})),
         (frozenset({1}), frozenset({1, 3})),
         (frozenset({1}), frozenset({3}), frozenset({1})),
         (frozenset({1, 3}), frozenset({1})),
         (frozenset({3}), frozenset({1}), frozenset({1}))]
        sage: list(_iterator_weight([]))
        [()]
    """
    # "weight" should be a dict mapping keys to weights
    if isinstance(weight, (list, tuple)):
        weight = {k+1: val for k, val in enumerate(weight) if val}

    # We first map the arbitrary keys to integers to combat unreliable
    # sorting behavior.
    keys = tuple(set(weight))
    multiset = []
    for i, key in enumerate(keys):
        multiset += [i] * weight[key]

    # We build ordered multiset partitions into sets of `X` by
    # permutation + deconcatenation
    for alpha in Permutations_mset(multiset):
        co = _break_at_descents(alpha, weak=True)
        for A in OrderedMultisetPartitionIntoSets(co).finer(strong=True):
            B = tuple([frozenset([keys[i] for i in block]) for block in A])
            yield B


def _iterator_size(size, length=None, alphabet=None):
    r"""
    An iterator for the ordered multiset partitions into sets of integer `n`.

    The degree `n` part of ordered multiset partitions into sets contains all
    sequences of subsets of `\NN_+` whose total sum adds up to `n`.

    If optional argument ``alphabet`` is given, it should be a ``Set`` object.
    Then only yield those `c` with all letters taken from ``alphabet``.

    OUTPUT:

    Tuples of ``frozenset`` objects representing ordered multiset partitions
    into sets.

    TESTS::

        sage: from sage.combinat.multiset_partition_into_sets_ordered import _iterator_size
        sage: OMP = OrderedMultisetPartitionsIntoSets(3)
        sage: list(map(OMP, _iterator_size(3)))
        [[{3}], [{1,2}], [{2}, {1}], [{1}, {2}], [{1}, {1}, {1}]]

        sage: OMP = OrderedMultisetPartitionsIntoSets(5, alphabet=(1,3))
        sage: list(map(OMP, _iterator_size(5, alphabet={1,3})))
        [[{1,3}, {1}], [{3}, {1}, {1}], [{1}, {1,3}], [{1}, {3}, {1}],
         [{1}, {1}, {3}], [{1}, {1}, {1}, {1}, {1}]]

        sage: list(_iterator_size(2))
        [(frozenset({2}),), (frozenset({1}), frozenset({1}))]
    """
    # iteration scheme:
    # - start from an integer composition ``alpha``.
    # - for each ``a`` in ``alpha``, pick distinct integers that sum to ``a``
    if alphabet:
        min_p = min(alphabet)
        max_p = max(alphabet)
        for alpha in IntegerListsLex(size, length=length, min_part=1,
                                     max_part=min(size, sum(alphabet))):
            for p in product(*[IntegerListsLex(a, min_slope=1,
                                               min_part=min_p,
                                               max_part=min(a, max_p))
                               for a in alpha]):
                if frozenset(_concatenate(p)).issubset(frozenset(alphabet)):
                    yield tuple(frozenset(k) for k in p)
    else:
        for alpha in IntegerListsLex(size, length=length, min_part=1, max_part=size):
            for p in product(*[IntegerListsLex(a, min_slope=1,
                                               min_part=1) for a in alpha]):
                yield tuple(frozenset(k) for k in p)


def _iterator_order(A, d, lengths=None):
    """
    An iterator for the ordered multiset partitions into sets of order `d`
    over alphabet `A`.

    If optional argument ``lengths`` is given, it should be a list of integers.
    Then only yield those ordered multiset partitions into sets with length
    in ``lengths``.

    OUTPUT:

    Tuples of ``frozenset`` objects representing ordered multiset partitions
    into sets.

    TESTS::

        sage: from sage.combinat.multiset_partition_into_sets_ordered import _iterator_order
        sage: OMP = OrderedMultisetPartitionsIntoSets([1,4], 3)
        sage: list(map(OMP, _iterator_order({1,4}, 3)))
        [[{1,4}, {1}], [{1,4}, {4}], [{1}, {1,4}], [{4}, {1,4}], [{1}, {1}, {1}],
         [{1}, {1}, {4}], [{1}, {4}, {1}], [{1}, {4}, {4}], [{4}, {1}, {1}],
         [{4}, {1}, {4}], [{4}, {4}, {1}], [{4}, {4}, {4}]]
        sage: list(map(OMP, _iterator_order([1,4], 3, [3])))
        [[{1}, {1}, {1}], [{1}, {1}, {4}], [{1}, {4}, {1}], [{1}, {4}, {4}],
         [{4}, {1}, {1}], [{4}, {1}, {4}], [{4}, {4}, {1}], [{4}, {4}, {4}]]

        sage: OMP = OrderedMultisetPartitionsIntoSets([1,2,4], 3)
        sage: list(map(OMP, _iterator_order([1,2,4], 3, [1,2])))[:10]
        [[{1,2,4}],  [{1,2}, {1}], [{1,2}, {2}], [{1,2}, {4}], [{1,4}, {1}],
         [{1,4}, {2}], [{1,4}, {4}], [{2,4}, {1}], [{2,4}, {2}], [{2,4}, {4}]]

        sage: list(_iterator_order([1,4], 3, [1]))
        []
        sage: list(_iterator_order([1,4], 3, [2]))
        [(frozenset({1, 4}), frozenset({1})), (frozenset({1, 4}), frozenset({4})),
         (frozenset({1}), frozenset({1, 4})), (frozenset({4}), frozenset({1, 4}))]
        sage: list(_iterator_order([1,4], 3, [4]))
        []
        sage: list(_iterator_order([1,4], 0, [3]))
        []
        sage: list(_iterator_order([1,4], 0, [0,3]))
        [()]
        sage: list(_iterator_order([1,4], 0))
        [()]
    """
    A = frozenset(A)

    # iteration scheme:
    # start from an integer composition ``alpha`` of ``d``.
    # for each ``a`` in ``alpha``, pick ``a`` letters from ``A``
    n = len(A)
    if not lengths:
        if d:
            lengths = range(max(1, d // n), d+1)
        else:
            lengths = (0,)

    for k in lengths:
        if not k and not d:
            yield ()
        else:
            for alpha in IntegerListsLex(d, length=k, min_part=1, max_part=n):
                for co in product(*[Subsets_sk(A, a) for a in alpha]):
                    yield tuple(frozenset(X) for X in co)


def _descents(w) -> list:
    r"""
    Return descent positions in the word ``w``.

    EXAMPLES::

        sage: from sage.combinat.multiset_partition_into_sets_ordered import _descents
        sage: _descents([1, 2, 3, 2, 2, 1, 4, 3]) == [2, 4, 6]
        True
        sage: _descents([])
        []
    """
    return [j for j in range(len(w) - 1) if w[j] > w[j + 1]]


def _break_at_descents(alpha, weak=True):
    r"""
    Return the deconcatenation of the composition ``alpha`` at its
    set of descents.

    OUTPUT:

    A list `[a_1, \ldots, a_r]` of nonempty lists whose concatenation
    is ``list(alpha)`` with the property that ``alpha[i] >= alpha[i+1]``
    if and only if positions `i` and `i+1` correspond to different
    lists. (Specifically, ``alpha[i]`` is the last letter of some
    `a_j` and ``alpha[i+1]`` is the first letter of `a_{j+1}`.)

    If the optional argument ``weak`` is ``False``, then only make
    breaks when ``alpha[i] > alpha[i+1]``.

    EXAMPLES::

        sage: from sage.combinat.multiset_partition_into_sets_ordered import _break_at_descents
        sage: _break_at_descents([1, 2, 3, 2, 2, 1, 4, 3])
        [[1, 2, 3], [2], [2], [1, 4], [3]]
        sage: _break_at_descents([1, 2, 3, 2, 2, 1, 4, 3], weak=False)
        [[1, 2, 3], [2, 2], [1, 4], [3]]
        sage: _break_at_descents([])
        []
    """
    if not alpha:
        return []

    Blocks = []
    block = [alpha[0]]
    for i in range(1,len(alpha)):
        if (alpha[i-1] > alpha[i]) or (alpha[i-1] == alpha[i] and weak):
            Blocks.append(block)
            block = [alpha[i]]
        else:
            block.append(alpha[i])
    if block:
        Blocks.append(block)
    return Blocks


def _refine_block(S, strong=False):
    r"""
    Return the list of all possible refinements of a set `S`.

    A refinement of `S` is a tuple of nonempty subsets whose union is `S`.

    If optional argument ``strong`` is set to ``True``, then only those
    refinements that are deconcatenations of the list ``sorted(S)`` are returned.

    (The subsets involved are stored as ``frozenset`` objects.)

    EXAMPLES::

        sage: from sage.combinat.multiset_partition_into_sets_ordered import _refine_block
        sage: _refine_block([1, 2], strong=True)
        [(frozenset({1}), frozenset({2})), (frozenset({1, 2}),)]

        sage: [tuple(Set(x) for x in tupl) for tupl in _refine_block([1, 2, 3], strong=True)]
        [({1}, {2}, {3}), ({1}, {2, 3}), ({1, 2}, {3}), ({1, 2, 3},)]
        sage: [tuple(Set(x) for x in tupl) for tupl in _refine_block([1, 2, 3])]
        [({3}, {2}, {1}), ({2}, {3}, {1}), ({3}, {1}, {2}),
         ({3}, {1, 2}), ({2}, {1}, {3}), ({2}, {1, 3}),
         ({2, 3}, {1}), ({1}, {3}, {2}), ({1}, {2}, {3}),
         ({1}, {2, 3}), ({1, 3}, {2}), ({1, 2}, {3}),
         ({1, 2, 3},)]

    TESTS::

        sage: len(_refine_block([1, 2, 3, 4])) == 1 + binomial(4,1)*2 + binomial(4,2) + binomial(4,2)*factorial(3) + factorial(4)
        True
        sage: len(_refine_block([1, 2, 3, 4], strong=True)) == 2**3
        True
        sage: _refine_block([])
        Traceback (most recent call last):
        ...
        ValueError: S (=[]) must be nonempty
    """
    if not S:
        raise ValueError("S (=%s) must be nonempty" % S)

    if all(s in ZZ for s in S):
        X = sorted(S)
    else:
        X = sorted(S, key=str)
    n = len(X)
    out = []
    if not strong:
        WordSet = IntegerListsLex(min_part=0, max_part=n-1, length=n)
    else:
        WordSet = IntegerListsLex(min_part=0, max_part=n-1, length=n, min_slope=0)

    for w in WordSet:
        if _is_initial_segment(sorted(set(w))):
            a = [frozenset() for _ in range(max(w)+1)]
            for pos in range(n):
                a[w[pos]] = a[w[pos]].union({X[pos]})
            out.append(tuple(a))
    return out


def _is_initial_segment(lst):
    r"""
    Return ``True`` if ``lst`` is an interval in `\ZZ` of the form `[0, 1, \ldots, n]`.

    EXAMPLES::

        sage: from sage.combinat.multiset_partition_into_sets_ordered import _is_initial_segment
        sage: _is_initial_segment([1, 2, 3])
        False
        sage: _is_initial_segment([0, 1, 2, 3])
        True
        sage: _is_initial_segment([0])
        True
    """
    return list(range(max(lst)+1)) == lst


def _split_block(S, k=2):
    """
    Return the list of all possible splittings of a set `S` into `k` parts.

    A splitting of `S` is a tuple of (possibly empty) subsets whose union is `S`.

    (The subsets involved are stored as ``frozenset`` objects.)

    EXAMPLES::

        sage: from sage.combinat.multiset_partition_into_sets_ordered import _split_block
        sage: S = [1, 2, 3]
        sage: _split_block(S, 1)
        [(frozenset({1, 2, 3}),)]

        sage: [tuple(Set(x) for x in tupl) for tupl in _split_block(S, 2)]
        [({}, {1, 2, 3}), ({3}, {1, 2}), ({2}, {1, 3}), ({2, 3}, {1}),
         ({1}, {2, 3}), ({1, 3}, {2}), ({1, 2}, {3}), ({1, 2, 3}, {})]
        sage: [tuple(Set(x) for x in tupl) for tupl in _split_block({2, 4}, 3)]
        [({}, {}, {2, 4}), ({}, {4}, {2}), ({4}, {}, {2}),
         ({}, {2}, {4}), ({}, {2, 4}, {}), ({4}, {2}, {}),
         ({2}, {}, {4}), ({2}, {4}, {}), ({2, 4}, {}, {})]
    """
    if all(s in ZZ for s in S):
        X = sorted(S)
    else:
        X = sorted(S, key=str)
    n = len(X)
    out = []
    for w in IntegerListsLex(min_part=0, max_part=k-1, length=n):
        a = [frozenset() for _ in range(k)]
        for pos in range(n):
            a[w[pos]] = a[w[pos]].union({X[pos]})
        out.append(tuple(a))
    return out


def _to_minimaj_blocks(T):
    r"""
    Return a tuple of tuples, representing an ordered multiset partition into sets
    in the minimaj ordering on blocks

    INPUT:

    - ``T`` -- a sequence of row words corresponding to (skew-)tableaux

    OUTPUT:

    The minimaj bijection `\phi^{-1}` of [BCHOPSY2017]_ applied to ``T``.

    EXAMPLES::

        sage: from sage.combinat.multiset_partition_into_sets_ordered import _to_minimaj_blocks
        sage: co = OrderedMultisetPartitionsIntoSets(14).an_element(); co
        [{2,3}, {2,3}, {4}]
        sage: co.to_tableaux_words()
        [[3, 2], [], [3, 2, 4]]
        sage: _to_minimaj_blocks(co.to_tableaux_words()) == co.minimaj_blocks()
        True
    """
    mu = [(i,) for i in T[-1]]
    breaks = [0] + _descents(T[-1]) + [len(mu)-1]
    T = [T[i][::-1] for i in range(len(T)-1)][::-1]
    for f in range(len(breaks)-1):
        for j in range(breaks[f],breaks[f+1]+1):
            mu[j] += tuple(i for i in T[f] if (mu[j][0] < i or j == breaks[f])
                                               and (j == breaks[f+1] or i <= mu[j+1][0]))
    return tuple(mu)


###############

class MinimajCrystal(UniqueRepresentation, Parent):
    r"""
    Crystal of ordered multiset partitions into sets with `ell` letters from
    alphabet `\{1, 2, \ldots, n\}` divided into `k` blocks.

    Elements are represented in the minimaj ordering of blocks as in
    Benkart et al. [BCHOPSY2017]_.

    .. NOTE::

        Elements are not stored internally as ordered multiset partitions
        into sets, but as certain (pairs of) words stemming from the minimaj
        bijection `\phi` of [BCHOPSY2017]_. See
        :class:`sage.combinat.multiset_partition_into_sets_ordered.MinimajCrystal.Element`
        for further details.

    AUTHORS:

    - Anne Schilling (2018): initial draft
    - Aaron Lauve (2018): changed to use ``Letters`` crystal for elements

    EXAMPLES::

        sage: list(crystals.Minimaj(2,3,2))                                             # needs sage.modules
        [((2, 1), (1,)), ((2,), (1, 2)), ((1,), (1, 2)), ((1, 2), (2,))]

        sage: b = crystals.Minimaj(3, 5, 2).an_element(); b                             # needs sage.modules
        ((2, 3, 1), (1, 2))
        sage: b.f(2)                                                                    # needs sage.modules
        ((2, 3, 1), (1, 3))
        sage: b.e(2)                                                                    # needs sage.modules
    """

    def __init__(self, n, ell, k):
        """
        Initialize ``self``.

        TESTS::

            sage: # needs sage.modules
            sage: B = crystals.Minimaj(2,3,2)
            sage: TestSuite(B).run()
            sage: B = crystals.Minimaj(3, 5, 2)
            sage: TestSuite(B).run()
            sage: list(crystals.Minimaj(2,6,3))
            [((1, 2), (2, 1), (1, 2))]
            sage: list(crystals.Minimaj(2,5,2))  # blocks too fat for alphabet
            []
            sage: list(crystals.Minimaj(4,2,3))  # more blocks than letters
            Traceback (most recent call last):
            ...
            ValueError: n (=4), ell (=2), and k (=3) must all be positive integers
        """
        Parent.__init__(self, category=ClassicalCrystals())
        self.n = n
        self.ell = ell
        self.k = k
        if not all([n in ZZ, ell in ZZ, k in ZZ]):
            raise TypeError("n (=%s), ell (=%s), and k (=%s) must all be positive integers" % (n, ell, k))
        if not all([n > 0, ell >= k, k > 0]):
            raise ValueError("n (=%s), ell (=%s), and k (=%s) must all be positive integers" % (n, ell, k))
        self._cartan_type = CartanType(['A',n-1])
        B = Letters(['A', n-1])
        T = tensor([B]*ell)
        self._BT = (B, T)
        self._OMPs = OrderedMultisetPartitionsIntoSets(n, ell, length=k)
        self.module_generators = []
        for co in self._OMPs:
            t = co.to_tableaux_words()
            word = T(*[B(a) for a in _concatenate(t)])
            blocks = [len(h) for h in t]
            breaks = tuple([0]+running_total(blocks))
            mu = self.element_class(self, (word, breaks))
            self.module_generators.append(mu)

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: B = crystals.Minimaj(3,4,2); B                                        # needs sage.modules
            Minimaj Crystal of type A_2 of words of length 4 into 2 blocks
        """
        return ("Minimaj Crystal of type A_%s of words of length %s into %s blocks"
                % (self.n-1, self.ell, self.k))

    def _an_element_(self):
        """
        Return a typical element of ``self``.

        EXAMPLES::

            sage: # needs sage.modules
            sage: B = crystals.Minimaj(4,5,3)
            sage: B.an_element()
            ((2, 3, 1), (1,), (1,))
            sage: B = crystals.Minimaj(2,2,1)
            sage: B.an_element()
            ((1, 2),)
            sage: B = crystals.Minimaj(1,2,1)
            sage: B.an_element()
            Traceback (most recent call last):
            ...
            EmptySetError
        """
        t = self._OMPs.an_element().to_tableaux_words()
        breaks = tuple([0]+running_total([len(h) for h in t]))
        B,T = self._BT
        return self.element_class(self, (T(*[B(a) for a in _concatenate(t)]), breaks))

    def _element_constructor_(self, x):
        """
        Build an element of Minimaj from an ordered multiset partition into sets.

        EXAMPLES::

            sage: # needs sage.modules
            sage: B1 = crystals.Minimaj(4,5,3); b = B1.an_element(); b
            ((2, 3, 1), (1,), (1,))
            sage: B1._element_constructor_(list(b))
            ((2, 3, 1), (1,), (1,))
            sage: B1._element_constructor_([[1,2,3], [2], [2]])
            ((3, 1, 2), (2,), (2,))
            sage: B2 = crystals.Minimaj(5,5,3)
            sage: B2._element_constructor_(b)
            ((2, 3, 1), (1,), (1,))
        """
        # Allow ``x`` to be either of:
        # - an ordered multiset partition into sets in ``self._OMPs``;
        # - an element of another Minimaj crystal with
        #   + same `ell` and `k`, and
        #   + all letters smaller or equal to ``self._n``.
        x = list(x)
        if x in self:
            t = self._OMPs(x).to_tableaux_words()
            breaks = tuple([0]+running_total([len(h) for h in t]))
            B,T = self._BT
            return self.element_class(self, (T(*[B(a) for a in _concatenate(t)]), breaks))
        else:
            raise ValueError("cannot convert %s into an element of %s" % (x, self))

    def __contains__(self, x):
        """
        Return ``True`` if ``x`` is an element of ``self`` or an ordered
        multiset partition into sets.

        EXAMPLES::

            sage: # needs sage.modules
            sage: B1 = crystals.Minimaj(2,5,3); b1 = B1.an_element(); b1
            ((1, 2), (2, 1), (1,))
            sage: B2 = crystals.Minimaj(5,5,3); b2 = B2.an_element(); b2
            ((2, 3, 1), (1,), (1,))
            sage: b2a = B2(((1,2), (1,), (1,2))); b2a
            ((2, 1), (1,), (1, 2))
            sage: b1 in B2
            True
            sage: b2 in B1
            False
            sage: b2a in B1
            True
        """
        if isinstance(x, MinimajCrystal.Element):
            if x.parent() == self:
                return True
            else:
                return list(x) in self._OMPs
        else:
            return x in self._OMPs

    def from_tableau(self, t):
        r"""
        Return the bijection `\phi^{-1}` of [BCHOPSY2017]_ applied to ``t``.

        INPUT:

        - ``t`` -- a sequence of column tableaux and a ribbon tableau

        EXAMPLES::

            sage: # needs sage.modules
            sage: B = crystals.Minimaj(3,6,3)
            sage: b = B.an_element(); b
            ((3, 1, 2), (2, 1), (1,))
            sage: t = b.to_tableaux_words(); t
            [[1], [2, 1], [], [3, 2, 1]]
            sage: B.from_tableau(t)
            ((3, 1, 2), (2, 1), (1,))
            sage: B.from_tableau(t) == b
            True

        TESTS::

            sage: # needs sage.modules
            sage: B = crystals.Minimaj(3,6,3)
            sage: all(mu == B.from_tableau(mu.to_tableaux_words()) for mu in B)
            True
            sage: t = B.an_element().to_tableaux_words()
            sage: B1 = crystals.Minimaj(3,6,2)
            sage: B1.from_tableau(t)
            Traceback (most recent call last):
            ...
            ValueError: ((3, 1, 2), (2, 1), (1,)) is not an element of
             Minimaj Crystal of type A_2 of words of length 6 into 2 blocks
        """
        mu = _to_minimaj_blocks(t)
        if mu in self:
            return self(mu)
        else:
            raise ValueError("%s is not an element of %s" % (mu, self))

    def val(self, q='q'):
        r"""
        Return the `Val` polynomial corresponding to ``self``.

        EXAMPLES:

        Verifying Example 4.5 from [BCHOPSY2017]_::

            sage: B = crystals.Minimaj(3, 4, 2)  # for `Val_{4,1}^{(3)}`                # needs sage.modules
            sage: B.val()                                                               # needs sage.modules
            (q^2+q+1)*s[2, 1, 1] + q*s[2, 2]
        """
        H = [self._OMPs(list(b)) for b in self.highest_weight_vectors()]
        Sym = SymmetricFunctions(ZZ[q])
        q = Sym.base_ring().gens()[0]
        s = Sym.schur()
        return sum((q**(t.minimaj()) * s[sorted(t.weight().values(), reverse=True)]
                   for t in H), Sym.zero())

    class Element(ElementWrapper):
        r"""
        An element of a Minimaj crystal.

        .. NOTE::

            Minimaj elements `b` are stored internally as pairs
            ``(w, breaks)``, where:

            - ``w`` -- a word of length ``self.parent().ell`` over the
              letters `1` up to ``self.parent().n``;
            - ``breaks`` is a list of de-concatenation points to turn ``w``
              into a list of row words of (skew-)tableaux that represent
              `b` under the minimaj bijection `\phi` of [BCHOPSY2017]_.

            The pair ``(w, breaks)`` may be recovered via ``b.value``.
        """

        def _repr_(self):
            """
            Return the string representation of ``self``.

            EXAMPLES::

                sage: crystals.Minimaj(4,5,3).an_element()                              # needs sage.modules
                ((2, 3, 1), (1,), (1,))
            """
            return repr(self._minimaj_blocks_from_word_pair())

        def __iter__(self):
            """
            Iterate over ``self._minimaj_blocks_from_word_pair()``.

            EXAMPLES::

                sage: b = crystals.Minimaj(4,5,3).an_element(); b                       # needs sage.modules
                ((2, 3, 1), (1,), (1,))
                sage: b.value                                                           # needs sage.modules
                ([1, 3, 2, 1, 1], (0, 1, 2, 5))
                sage: list(b)                                                           # needs sage.modules
                [(2, 3, 1), (1,), (1,)]
            """
            return self._minimaj_blocks_from_word_pair().__iter__()

        def _latex_(self):
            r"""
            Return the latex representation of ``self``.

            EXAMPLES::

                sage: b = crystals.Minimaj(4,5,3).an_element(); b                       # needs sage.modules
                ((2, 3, 1), (1,), (1,))
                sage: latex(b)                                                          # needs sage.modules
                \left(\left(2, 3, 1\right), \left(1\right), \left(1\right)\right)
            """
            return latex(self._minimaj_blocks_from_word_pair())

        def _minimaj_blocks_from_word_pair(self):
            """
            Return the tuple of tuples (in the minimaj ordering on blocks of
            ordered multiset partitions into sets) corresponding to ``self``.

            EXAMPLES::

                sage: B = crystals.Minimaj(4,5,3)                                       # needs sage.modules
                sage: b = B.an_element(); b.value                                       # needs sage.modules
                ([1, 3, 2, 1, 1], (0, 1, 2, 5))
                sage: b._minimaj_blocks_from_word_pair()                                # needs sage.modules
                ((2, 3, 1), (1,), (1,))
            """
            return _to_minimaj_blocks(self.to_tableaux_words())

        def to_tableaux_words(self):
            r"""
            Return the image of the ordered multiset partition into sets ``self``
            under the minimaj bijection `\phi` of [BCHOPSY2017]_.

            EXAMPLES::

                sage: # needs sage.modules
                sage: B = crystals.Minimaj(4,5,3)
                sage: b = B.an_element(); b
                ((2, 3, 1), (1,), (1,))
                sage: b.to_tableaux_words()
                [[1], [3], [2, 1, 1]]
                sage: b = B([[1,3,4], [3], [3]]); b
                ((4, 1, 3), (3,), (3,))
                sage: b.to_tableaux_words()
                [[3, 1], [], [4, 3, 3]]
            """
            w, breaks = self.value
            return [[ZZ(w[a].value) for a in range(breaks[j], breaks[j+1])]
                        for j in range(len(breaks)-1)]

        def e(self, i):
            r"""
            Return `e_i` on ``self``.

            EXAMPLES::

                sage: B = crystals.Minimaj(4,3,2)                                       # needs sage.modules
                sage: b = B([[2,3], [3]]); b                                            # needs sage.modules
                ((2, 3), (3,))
                sage: [b.e(i) for i in range(1,4)]                                      # needs sage.modules
                [((1, 3), (3,)), ((2,), (2, 3)), None]
            """
            P = self.parent()
            w, breaks = self.value
            if w.e(i) is None:
                return None
            w = w.e(i)
            return P.element_class(P, (w, breaks))

        def f(self,i):
            r"""
            Return `f_i` on ``self``.

            EXAMPLES::

                sage: B = crystals.Minimaj(4,3,2)                                       # needs sage.modules
                sage: b = B([[2,3], [3]]); b                                            # needs sage.modules
                ((2, 3), (3,))
                sage: [b.f(i) for i in range(1,4)]                                      # needs sage.modules
                [None, None, ((2, 3), (4,))]
            """
            P = self.parent()
            w, breaks = self.value
            if w.f(i) is None:
                return None
            w = w.f(i)
            return P.element_class(P, (w, breaks))
