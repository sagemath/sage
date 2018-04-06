r"""
Ordered Set Partitions

AUTHORS:

- Mike Hansen

- MuPAD-Combinat developers (for algorithms and design inspiration)

- Travis Scrimshaw (2013-02-28): Removed ``CombinatorialClass`` and added
  entry point through :class:`OrderedSetPartition`.
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from six import add_metaclass

from sage.arith.all import factorial
import sage.rings.integer
from sage.sets.set import Set, is_Set
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.all import prod
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.list_clone import ClonableArray
from sage.combinat.combinatorial_map import combinatorial_map
from sage.combinat.combinat import stirling_number2
from sage.combinat.composition import Composition, Compositions
import sage.combinat.permutation as permutation
from functools import reduce


@add_metaclass(InheritComparisonClasscallMetaclass)
class OrderedSetPartition(ClonableArray):
    """
    An ordered partition of a set.

    An ordered set partition `p` of a set `s` is a list of pairwise
    disjoint nonempty subsets of `s` such that the union of these
    subsets is `s`. These subsets are called the parts of the partition.
    We represent an ordered set partition as a list of sets. By
    extension, an ordered set partition of a nonnegative integer `n` is
    the set partition of the integers from `1` to `n`. The number of
    ordered set partitions of `n` is called the `n`-th ordered Bell
    number.

    There is a natural integer composition associated with an ordered
    set partition, that is the sequence of sizes of all its parts in
    order.

    The number `T_n` of ordered set partitions of
    `\{ 1, 2, ..., n \}` is the so-called `n`-th *Fubini number* (also
    known as the `n`-th ordered Bell number; see
    :wikipedia:`Ordered Bell number`). Its exponential generating
    function is
    
    .. MATH::

        \sum_n {T_n \over n!} x^n = {1 \over 2-e^x}.

    (See sequence A000670 in OEIS.)

    EXAMPLES:

    There are 13 ordered set partitions of `\{1,2,3\}`::

        sage: OrderedSetPartitions(3).cardinality()
        13

    Here is the list of them::

        sage: OrderedSetPartitions(3).list()
        [[{1}, {2}, {3}],
         [{1}, {3}, {2}],
         [{2}, {1}, {3}],
         [{3}, {1}, {2}],
         [{2}, {3}, {1}],
         [{3}, {2}, {1}],
         [{1}, {2, 3}],
         [{2}, {1, 3}],
         [{3}, {1, 2}],
         [{1, 2}, {3}],
         [{1, 3}, {2}],
         [{2, 3}, {1}],
         [{1, 2, 3}]]

    There are 12 ordered set partitions of `\{1,2,3,4\}` whose underlying
    composition is `[1,2,1]`::

        sage: OrderedSetPartitions(4,[1,2,1]).list()
        [[{1}, {2, 3}, {4}],
         [{1}, {2, 4}, {3}],
         [{1}, {3, 4}, {2}],
         [{2}, {1, 3}, {4}],
         [{2}, {1, 4}, {3}],
         [{3}, {1, 2}, {4}],
         [{4}, {1, 2}, {3}],
         [{3}, {1, 4}, {2}],
         [{4}, {1, 3}, {2}],
         [{2}, {3, 4}, {1}],
         [{3}, {2, 4}, {1}],
         [{4}, {2, 3}, {1}]]

    Since :trac:`14140`, we can create an ordered set partition directly by
    :class:`OrderedSetPartition` which creates the parent object by taking the
    union of the partitions passed in. However it is recommended and
    (marginally) faster to create the parent first and then create the ordered
    set partition from that. ::

        sage: s = OrderedSetPartition([[1,3],[2,4]]); s
        [{1, 3}, {2, 4}]
        sage: s.parent()
        Ordered set partitions of {1, 2, 3, 4}

    REFERENCES:

    :wikipedia:`Ordered_partition_of_a_set`
    """
    @staticmethod
    def __classcall_private__(cls, parts):
        """
        Create a set partition from ``parts`` with the appropriate parent.

        EXAMPLES::

            sage: s = OrderedSetPartition([[1,3],[2,4]]); s
            [{1, 3}, {2, 4}]
            sage: s.parent()
            Ordered set partitions of {1, 2, 3, 4}
            sage: t = OrderedSetPartition([[2,4],[1,3]]); t
            [{2, 4}, {1, 3}]
            sage: s != t
            True
            sage: OrderedSetPartition([])
            []
        """
        P = OrderedSetPartitions( reduce(lambda x,y: x.union(y), map(Set, parts), Set([])) )
        return P.element_class(P, parts)

    def __init__(self, parent, s):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: OS = OrderedSetPartitions(4)
            sage: s = OS([[1, 3], [2, 4]])
            sage: TestSuite(s).run()
        """
        ClonableArray.__init__(self, parent, [Set(_) for _ in s])

    def check(self):
        """
        Check that we are a valid ordered set partition.

        EXAMPLES::

            sage: OS = OrderedSetPartitions(4)
            sage: s = OS([[1, 3], [2, 4]])
            sage: s.check()
        """
        assert self in self.parent()

    @combinatorial_map(name='to composition')
    def to_composition(self):
        r"""
        Return the integer composition whose parts are the sizes of the sets
        in ``self``.

        EXAMPLES::

            sage: S = OrderedSetPartitions(5)
            sage: x = S([[3,5,4], [1, 2]])
            sage: x.to_composition()
            [3, 2]
            sage: y = S([[3,1], [2], [5,4]])
            sage: y.to_composition()
            [2, 1, 2]
        """
        return Composition([len(_) for _ in self])

class OrderedSetPartitions(UniqueRepresentation, Parent):
    """
    Return the combinatorial class of ordered set partitions of ``s``.

    EXAMPLES::

        sage: OS = OrderedSetPartitions([1,2,3,4]); OS
        Ordered set partitions of {1, 2, 3, 4}
        sage: OS.cardinality()
        75
        sage: OS.first()
        [{1}, {2}, {3}, {4}]
        sage: OS.last()
        [{1, 2, 3, 4}]
        sage: OS.random_element()
        [{3}, {1}, {2}, {4}]

    ::

        sage: OS = OrderedSetPartitions([1,2,3,4], [2,2]); OS
        Ordered set partitions of {1, 2, 3, 4} into parts of size [2, 2]
        sage: OS.cardinality()
        6
        sage: OS.first()
        [{1, 2}, {3, 4}]
        sage: OS.last()
        [{3, 4}, {1, 2}]
        sage: OS.list()
        [[{1, 2}, {3, 4}],
         [{1, 3}, {2, 4}],
         [{1, 4}, {2, 3}],
         [{2, 3}, {1, 4}],
         [{2, 4}, {1, 3}],
         [{3, 4}, {1, 2}]]

    ::

        sage: OS = OrderedSetPartitions("cat"); OS
        Ordered set partitions of {'a', 'c', 't'}
        sage: OS.list()
        [[{'a'}, {'c'}, {'t'}],
         [{'a'}, {'t'}, {'c'}],
         [{'c'}, {'a'}, {'t'}],
         [{'t'}, {'a'}, {'c'}],
         [{'c'}, {'t'}, {'a'}],
         [{'t'}, {'c'}, {'a'}],
         [{'a'}, {'c', 't'}],
         [{'c'}, {'a', 't'}],
         [{'t'}, {'a', 'c'}],
         [{'a', 'c'}, {'t'}],
         [{'a', 't'}, {'c'}],
         [{'c', 't'}, {'a'}],
         [{'a', 'c', 't'}]]
    """
    @staticmethod
    def __classcall_private__(cls, s, c=None):
        """
        Choose the correct parent based upon input.

        EXAMPLES::

            sage: OrderedSetPartitions(4)
            Ordered set partitions of {1, 2, 3, 4}
            sage: OrderedSetPartitions(4, [1, 2, 1])
            Ordered set partitions of {1, 2, 3, 4} into parts of size [1, 2, 1]
        """
        if isinstance(s, (int, sage.rings.integer.Integer)):
            if s < 0:
                raise ValueError("s must be non-negative")
            s = frozenset(range(1, s+1))
        else:
            s = frozenset(s)

        if c is None:
            return OrderedSetPartitions_s(s)

        if isinstance(c, (int, sage.rings.integer.Integer)):
            return OrderedSetPartitions_sn(s, c)
        if c not in Compositions(len(s)):
            raise ValueError("c must be a composition of %s"%len(s))
        return OrderedSetPartitions_scomp(s, Composition(c))

    def __init__(self, s):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: OS = OrderedSetPartitions(4)
            sage: TestSuite(OS).run()
        """
        self._set = s
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def _element_constructor_(self, s):
        """
        Construct an element of ``self`` from ``s``.

        EXAMPLES::

            sage: OS = OrderedSetPartitions(4)
            sage: OS([[1,3],[2,4]])
            [{1, 3}, {2, 4}]
        """
        if isinstance(s, OrderedSetPartition):
            if s.parent() == self:
                return s
            raise ValueError("Cannot convert %s into an element of %s"%(s, self))
        return self.element_class(self, list(s))

    Element = OrderedSetPartition

    def __contains__(self, x):
        """
        TESTS::

            sage: OS = OrderedSetPartitions([1,2,3,4])
            sage: all(sp in OS for sp in OS)
            True
        """
        #x must be a list
        if not isinstance(x, (OrderedSetPartition, list, tuple)):
            return False

        #The total number of elements in the list
        #should be the same as the number is self._set
        if sum(map(len, x)) != len(self._set):
            return False

        #Check to make sure each element of the list
        #is a set
        u = Set([])
        for s in x:
            if not isinstance(s, (set, frozenset)) and not is_Set(s):
                return False
            u = u.union(s)

        #Make sure that the union of all the
        #sets is the original set
        if u != Set(self._set):
            return False

        return True

class OrderedSetPartitions_s(OrderedSetPartitions):
    """
    Class of ordered partitions of a set `S`.
    """
    def _repr_(self):
        """
        TESTS::

            sage: OrderedSetPartitions([1,2,3,4])
            Ordered set partitions of {1, 2, 3, 4}
        """
        return "Ordered set partitions of %s"%Set(self._set)

    def cardinality(self):
        """
        EXAMPLES::

            sage: OrderedSetPartitions(0).cardinality()
            1
            sage: OrderedSetPartitions(1).cardinality()
            1
            sage: OrderedSetPartitions(2).cardinality()
            3
            sage: OrderedSetPartitions(3).cardinality()
            13
            sage: OrderedSetPartitions([1,2,3]).cardinality()
            13
            sage: OrderedSetPartitions(4).cardinality()
            75
            sage: OrderedSetPartitions(5).cardinality()
            541
        """
        return sum([factorial(k)*stirling_number2(len(self._set),k) for k in range(len(self._set)+1)])

    def __iter__(self):
        """
        EXAMPLES::

            sage: [ p for p in OrderedSetPartitions([1,2,3]) ]
            [[{1}, {2}, {3}],
             [{1}, {3}, {2}],
             [{2}, {1}, {3}],
             [{3}, {1}, {2}],
             [{2}, {3}, {1}],
             [{3}, {2}, {1}],
             [{1}, {2, 3}],
             [{2}, {1, 3}],
             [{3}, {1, 2}],
             [{1, 2}, {3}],
             [{1, 3}, {2}],
             [{2, 3}, {1}],
             [{1, 2, 3}]]
        """
        for x in Compositions(len(self._set)):
            for z in OrderedSetPartitions(self._set, x):
                yield self.element_class(self, z)

class OrderedSetPartitions_sn(OrderedSetPartitions):
    def __init__(self, s, n):
        """
        TESTS::

            sage: OS = OrderedSetPartitions([1,2,3,4], 2)
            sage: OS == loads(dumps(OS))
            True
        """
        OrderedSetPartitions.__init__(self, s)
        self.n = n

    def __contains__(self, x):
        """
        TESTS::

            sage: OS = OrderedSetPartitions([1,2,3,4], 2)
            sage: all(sp in OS for sp in OS)
            True
            sage: OS.cardinality()
            14
            sage: len(filter(lambda x: x in OS, OrderedSetPartitions([1,2,3,4])))
            14
        """
        return OrderedSetPartitions.__contains__(self, x) and len(x) == self.n

    def __repr__(self):
        """
        TESTS::

            sage: OrderedSetPartitions([1,2,3,4], 2)
            Ordered set partitions of {1, 2, 3, 4} into 2 parts
        """
        return "Ordered set partitions of %s into %s parts"%(Set(self._set),self.n)

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        The number of ordered partitions of a set of size `n` into `k`
        parts is equal to `k! S(n,k)` where `S(n,k)` denotes the Stirling
        number of the second kind.

        EXAMPLES::

            sage: OrderedSetPartitions(4,2).cardinality()
            14
            sage: OrderedSetPartitions(4,1).cardinality()
            1
        """
        return factorial(self.n)*stirling_number2(len(self._set), self.n)

    def __iter__(self):
        """
        EXAMPLES::

            sage: [ p for p in OrderedSetPartitions([1,2,3,4], 2) ]
            [[{1, 2, 3}, {4}],
             [{1, 2, 4}, {3}],
             [{1, 3, 4}, {2}],
             [{2, 3, 4}, {1}],
             [{1, 2}, {3, 4}],
             [{1, 3}, {2, 4}],
             [{1, 4}, {2, 3}],
             [{2, 3}, {1, 4}],
             [{2, 4}, {1, 3}],
             [{3, 4}, {1, 2}],
             [{1}, {2, 3, 4}],
             [{2}, {1, 3, 4}],
             [{3}, {1, 2, 4}],
             [{4}, {1, 2, 3}]]
        """
        for x in Compositions(len(self._set),length=self.n):
            for z in OrderedSetPartitions_scomp(self._set,x):
                yield self.element_class(self, z)

class OrderedSetPartitions_scomp(OrderedSetPartitions):
    def __init__(self, s, comp):
        """
        TESTS::

            sage: OS = OrderedSetPartitions([1,2,3,4], [2,1,1])
            sage: OS == loads(dumps(OS))
            True
        """
        OrderedSetPartitions.__init__(self, s)
        self.c = Composition(comp)

    def __repr__(self):
        """
        TESTS::

            sage: OrderedSetPartitions([1,2,3,4], [2,1,1])
            Ordered set partitions of {1, 2, 3, 4} into parts of size [2, 1, 1]
        """
        return "Ordered set partitions of %s into parts of size %s"%(Set(self._set), self.c)

    def __contains__(self, x):
        """
        TESTS::

            sage: OS = OrderedSetPartitions([1,2,3,4], [2,1,1])
            sage: all(sp in OS for sp in OS)
            True
            sage: OS.cardinality()
            12
            sage: len(filter(lambda x: x in OS, OrderedSetPartitions([1,2,3,4])))
            12
        """
        return OrderedSetPartitions.__contains__(self, x) and [len(_) for _ in x] == self.c

    def cardinality(self):
        r"""
        Return the cardinality of ``self``.

        The number of ordered set partitions of a set of length `k` with
        composition shape `\mu` is equal to

        .. MATH::

            \frac{k!}{\prod_{\mu_i \neq 0} \mu_i!}.

        EXAMPLES::

            sage: OrderedSetPartitions(5,[2,3]).cardinality()
            10
            sage: OrderedSetPartitions(0, []).cardinality()
            1
            sage: OrderedSetPartitions(0, [0]).cardinality()
            1
            sage: OrderedSetPartitions(0, [0,0]).cardinality()
            1
            sage: OrderedSetPartitions(5, [2,0,3]).cardinality()
            10
        """
        return factorial(len(self._set))/prod([factorial(i) for i in self.c])

    def __iter__(self):
        """
        TESTS::

            sage: [ p for p in OrderedSetPartitions([1,2,3,4], [2,1,1]) ]
            [[{1, 2}, {3}, {4}],
             [{1, 2}, {4}, {3}],
             [{1, 3}, {2}, {4}],
             [{1, 4}, {2}, {3}],
             [{1, 3}, {4}, {2}],
             [{1, 4}, {3}, {2}],
             [{2, 3}, {1}, {4}],
             [{2, 4}, {1}, {3}],
             [{3, 4}, {1}, {2}],
             [{2, 3}, {4}, {1}],
             [{2, 4}, {3}, {1}],
             [{3, 4}, {2}, {1}]]

            sage: len(OrderedSetPartitions([1,2,3,4], [1,1,1,1]))
            24

            sage: [ x for x in OrderedSetPartitions([1,4,7], [3]) ]
            [[{1, 4, 7}]]

            sage: [ x for x in OrderedSetPartitions([1,4,7], [1,2]) ]
            [[{1}, {4, 7}], [{4}, {1, 7}], [{7}, {1, 4}]]

            sage: [ p for p in OrderedSetPartitions([], []) ]
            [[]]

            sage: [ p for p in OrderedSetPartitions([1], [1]) ]
            [[{1}]]

        Let us check that it works for large size (:trac:`16646`)::

            sage: OrderedSetPartitions(42).first()
            [{1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}, {9}, {10}, {11}, {12},
            {13}, {14}, {15}, {16}, {17}, {18}, {19}, {20}, {21}, {22}, {23},
            {24}, {25}, {26}, {27}, {28}, {29}, {30}, {31}, {32}, {33}, {34},
            {35}, {36}, {37}, {38}, {39}, {40}, {41}, {42}]
        """
        comp = self.c
        lset = [x for x in self._set]
        l = len(self.c)
        dcomp = [-1] + comp.descents(final_descent=True)

        p = []
        for j in range(l):
            p += [j + 1] * comp[j]

        for x in permutation.Permutations(p):
            res = permutation.to_standard(x).inverse()
            res = [lset[x - 1] for x in res]
            yield self.element_class(self, [Set(res[dcomp[i]+1:dcomp[i+1]+1])
                                            for i in range(l)])

##########################################################
# Deprecations


class SplitNK(OrderedSetPartitions_scomp):
    def __setstate__(self, state):
        r"""
        For unpickling old ``SplitNK`` objects.

        TESTS::

            sage: loads("x\x9ck`J.NLO\xd5K\xce\xcfM\xca\xccK,\xd1+.\xc8\xc9,"
            ....:   "\x89\xcf\xcb\xe6\n\x061\xfc\xbcA\xccBF\xcd\xc6B\xa6\xda"
            ....:   "Bf\x8dP\xa6\xf8\xbcB\x16\x88\x96\xa2\xcc\xbc\xf4b\xbd\xcc"
            ....:   "\xbc\x92\xd4\xf4\xd4\"\xae\xdc\xc4\xec\xd4x\x18\xa7\x905"
            ....:   "\x94\xd1\xb45\xa8\x90\r\xa8>\xbb\x90=\x03\xc85\x02r9J\x93"
            ....:   "\xf4\x00\xb4\xc6%f")
            Ordered set partitions of {0, 1, 2, 3, 4} into parts of size [2, 3]
        """
        self.__class__ = OrderedSetPartitions_scomp
        n = state['_n']
        k = state['_k']
        OrderedSetPartitions_scomp.__init__(self, range(state['_n']), (k,n-k))

from sage.structure.sage_object import register_unpickle_override
register_unpickle_override("sage.combinat.split_nk", "SplitNK_nk", SplitNK)

