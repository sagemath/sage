# sage.doctest: needs sage.combinat sage.modules
r"""
Diagram and Partition Algebras

AUTHORS:

- Mike Hansen (2007): Initial version
- Stephen Doty, Aaron Lauve, George H. Seelinger (2012): Implementation of
  partition, Brauer, Temperley--Lieb, and ideal partition algebras
- Stephen Doty, Aaron Lauve, George H. Seelinger (2015): Implementation of
  ``*Diagram`` classes and other methods to improve diagram algebras.
- Mike Zabrocki (2018): Implementation of individual element diagram classes
- Aaron Lauve, Mike Zabrocki (2018): Implementation of orbit basis for Partition algebra.
- Travis Scrimshaw (2024): Implemented the Potts representation of the
  partition algebra.
"""

# ****************************************************************************
#  Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#                2012 Stephen Doty <doty@math.luc.edu>,
#                     Aaron Lauve <lauve@math.luc.edu>,
#                     George H. Seelinger <ghseeli@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ***************************************************************************

import itertools

from sage.categories.associative_algebras import AssociativeAlgebras
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.arith.power import generic_power
from sage.combinat.free_module import CombinatorialFreeModule
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.combinat.combinat import bell_number, catalan_number
from sage.structure.global_options import GlobalOptions
from sage.combinat.combinat_cython import (perfect_matchings_iterator,
                                           set_partition_composition)
from sage.combinat.set_partition import SetPartitions, AbstractSetPartition
from sage.combinat.set_partition_iterator import set_partition_iterator
from sage.combinat.permutation import Permutations
from sage.misc.cachefunc import cached_method
from sage.misc.flatten import flatten
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.lazy_import import lazy_import
from sage.misc.misc_c import prod
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.arith.misc import integer_floor as floor
from sage.arith.misc import integer_ceil as ceil

lazy_import('sage.graphs.graph', 'Graph')
lazy_import('sage.combinat.symmetric_group_algebra', 'SymmetricGroupAlgebra_n')


def partition_diagrams(k):
    r"""
    Return a generator of all partition diagrams of order ``k``.

    A partition diagram of order `k \in \ZZ` to is a set partition of
    `\{1, \ldots, k, -1, \ldots, -k\}`. If we have `k - 1/2 \in ZZ`, then
    a partition diagram of order `k \in 1/2 \ZZ` is a set partition of
    `\{1, \ldots, k+1/2, -1, \ldots, -(k+1/2)\}` with `k+1/2` and `-(k+1/2)`
    in the same block. See [HR2005]_.

    INPUT:

    - ``k`` -- the order of the partition diagrams

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: [SetPartition(p) for p in da.partition_diagrams(2)]
        [{{-2, -1, 1, 2}},
         {{-2, 1, 2}, {-1}},
         {{-2}, {-1, 1, 2}},
         {{-2, -1}, {1, 2}},
         {{-2}, {-1}, {1, 2}},
         {{-2, -1, 1}, {2}},
         {{-2, 1}, {-1, 2}},
         {{-2, 1}, {-1}, {2}},
         {{-2, 2}, {-1, 1}},
         {{-2, -1, 2}, {1}},
         {{-2, 2}, {-1}, {1}},
         {{-2}, {-1, 1}, {2}},
         {{-2}, {-1, 2}, {1}},
         {{-2, -1}, {1}, {2}},
         {{-2}, {-1}, {1}, {2}}]
        sage: [SetPartition(p) for p in da.partition_diagrams(3/2)]
        [{{-2, -1, 1, 2}},
         {{-2, 1, 2}, {-1}},
         {{-2, 2}, {-1, 1}},
         {{-2, -1, 2}, {1}},
         {{-2, 2}, {-1}, {1}}]
    """
    if k in ZZ:
        S = set_partition_iterator(list(range(1, k+1)) + list(range(-k,0)))
        for p in S:
            yield p
    elif k + ZZ(1)/ZZ(2) in ZZ:  # Else k in 1/2 ZZ
        k = ZZ(k + ZZ(1) / ZZ(2))
        S = set_partition_iterator(list(range(1, k+1)) + list(range(-k+1,0)))
        for p in S:
            yield [b + [-k] if k in b else b for b in p]
    else:
        raise ValueError("argument %s must be a half-integer" % k)


def brauer_diagrams(k):
    r"""
    Return a generator of all Brauer diagrams of order ``k``.

    A Brauer diagram of order `k` is a partition diagram of order `k`
    with block size 2.

    INPUT:

    - ``k`` -- the order of the Brauer diagrams

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: [SetPartition(p) for p in da.brauer_diagrams(2)]
        [{{-2, -1}, {1, 2}}, {{-2, 1}, {-1, 2}}, {{-2, 2}, {-1, 1}}]
        sage: [SetPartition(p) for p in da.brauer_diagrams(5/2)]
        [{{-3, 3}, {-2, -1}, {1, 2}},
         {{-3, 3}, {-2, 1}, {-1, 2}},
         {{-3, 3}, {-2, 2}, {-1, 1}}]
    """
    if k in ZZ:
        s = list(range(1, k+1)) + list(range(-k,0))
        for p in perfect_matchings_iterator(k):
            yield [(s[a],s[b]) for a,b in p]
    elif k + ZZ.one() / 2 in ZZ:  # Else k in 1/2 ZZ
        k = ZZ(k + ZZ.one() / 2)
        s = list(range(1, k)) + list(range(-k+1,0))
        for p in perfect_matchings_iterator(k-1):
            yield [(s[a], s[b]) for a, b in p] + [[k, -k]]


def temperley_lieb_diagrams(k):
    r"""
    Return a generator of all Temperley--Lieb diagrams of order ``k``.

    A Temperley--Lieb diagram of order `k` is a partition diagram of order `k`
    with block size  2 and is planar.

    INPUT:

    - ``k`` -- the order of the Temperley--Lieb diagrams

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: [SetPartition(p) for p in da.temperley_lieb_diagrams(2)]
        [{{-2, -1}, {1, 2}}, {{-2, 2}, {-1, 1}}]
        sage: [SetPartition(p) for p in da.temperley_lieb_diagrams(5/2)]
        [{{-3, 3}, {-2, -1}, {1, 2}}, {{-3, 3}, {-2, 2}, {-1, 1}}]
    """
    for i in brauer_diagrams(k):
        if is_planar(i):
            yield i


def planar_diagrams(k):
    r"""
    Return a generator of all planar diagrams of order ``k``.

    A planar diagram of order `k` is a partition diagram of order `k`
    that has no crossings.

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: all_diagrams = [SetPartition(p) for p in da.partition_diagrams(2)]
        sage: da2 = [SetPartition(p) for p in da.planar_diagrams(2)]
        sage: [p for p in all_diagrams if p not in da2]
        [{{-2, 1}, {-1, 2}}]
        sage: all_diagrams = [SetPartition(p) for p in da.partition_diagrams(5/2)]
        sage: da5o2 = [SetPartition(p) for p in da.planar_diagrams(5/2)]
        sage: [p for p in all_diagrams if p not in da5o2]
        [{{-3, -1, 3}, {-2, 1, 2}},
         {{-3, -2, 1, 3}, {-1, 2}},
         {{-3, -1, 1, 3}, {-2, 2}},
         {{-3, 1, 3}, {-2, -1, 2}},
         {{-3, 1, 3}, {-2, 2}, {-1}},
         {{-3, 1, 3}, {-2}, {-1, 2}},
         {{-3, -1, 2, 3}, {-2, 1}},
         {{-3, 3}, {-2, 1}, {-1, 2}},
         {{-3, -1, 3}, {-2, 1}, {2}},
         {{-3, -1, 3}, {-2, 2}, {1}}]
    """
    if k in ZZ:
        X = list(range(1,k+1)) + list(range(-k,0))
        yield from planar_partitions_rec(X)
    elif k + ZZ(1)/ZZ(2) in ZZ:  # Else k in 1/2 ZZ
        k = ZZ(k + ZZ(1) / ZZ(2))
        X = list(range(1,k+1)) + list(range(-k+1,0))
        for Y in planar_partitions_rec(X):
            Y = list(Y)
            for part in Y:
                if k in part:
                    part.append(-k)
                    break
            yield Y
    else:
        raise ValueError("argument %s must be a half-integer" % k)


def planar_partitions_rec(X):
    r"""
    Iterate over all planar set partitions of ``X`` by using a
    recursive algorithm.

    ALGORITHM:

    To construct the set partition `\rho = \{\rho_1, \ldots, \rho_k\}` of
    `[n]`, we remove the part of the set partition containing the last
    element of ``X``, which, we consider to be `\rho_k = \{i_1, \ldots, i_m\}`
    without loss of generality. The remaining parts come from the planar set
    partitions of `\{1, \ldots, i_1-1\}, \{i_1+1, \ldots, i_2-1\}, \ldots,
    \{i_m+1, \ldots, n\}`.

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: list(da.planar_partitions_rec([1,2,3]))
        [([1, 2], [3]), ([1], [2], [3]), ([2], [1, 3]), ([1], [2, 3]), ([1, 2, 3],)]
    """
    if not X:
        return
    if len(X) <= 2:
        # Direct implementation of small cases
        yield (X,)
        if len(X) > 1:
            yield ([X[0]], [X[1]])
        return
    from sage.combinat.subset import powerset
    from itertools import product
    for S in powerset(range(len(X)-1)):
        if not S:
            for Y in planar_partitions_rec(X[:-1]):
                yield Y + ([X[-1]],)
            continue
        last = [X[i] for i in S]
        last.append(X[-1])
        pt = []
        if S[0] != 0:
            pt += [X[:S[0]]]
        pt = [X[S[i]+1:S[i+1]] for i in range(len(S)-1) if S[i]+1 != S[i+1]]
        if S[-1] + 1 != len(X) - 1:
            pt += [X[S[-1]+1:-1]]
        parts = [planar_partitions_rec(X[S[i]+1:S[i+1]]) for i in range(len(S)-1)
                 if S[i] + 1 != S[i+1]]
        if S[0] != 0:
            parts.append(planar_partitions_rec(X[:S[0]]))
        if S[-1] + 1 != len(X) - 1:
            parts.append(planar_partitions_rec(X[S[-1]+1:-1]))
        for Y in product(*parts):
            yield sum(Y, ()) + (last,)


def ideal_diagrams(k):
    r"""
    Return a generator of all "ideal" diagrams of order ``k``.

    An ideal diagram of order `k` is a partition diagram of order `k` with
    propagating number less than `k`.

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: all_diagrams = da.partition_diagrams(2)
        sage: [SetPartition(p) for p in all_diagrams if p not in da.ideal_diagrams(2)]
        [{{-2, 1}, {-1, 2}}, {{-2, 2}, {-1, 1}}]

        sage: all_diagrams = da.partition_diagrams(3/2)
        sage: [SetPartition(p) for p in all_diagrams if p not in da.ideal_diagrams(3/2)]
        [{{-2, 2}, {-1, 1}}]
    """
    for i in partition_diagrams(k):
        if propagating_number(i) < k:
            yield i


class AbstractPartitionDiagram(AbstractSetPartition):
    r"""
    Abstract base class for partition diagrams.

    This class represents a single partition diagram, that is used as a
    basis key for a diagram algebra element. A partition diagram should
    be a partition of the set  `\{1, \ldots, k, -1, \ldots, -k\}`. Each
    such set partition is regarded as a graph on nodes
    `\{1, \ldots, k, -1, \ldots, -k\}` arranged in two rows, with nodes
    `1, \ldots, k` in the top row from left to right and with nodes
    `-1, \ldots, -k` in the bottom row from left to right, and an edge
    connecting two nodes if and only if the nodes lie in the same
    subset of the set partition.

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: pd = da.AbstractPartitionDiagrams(2)
        sage: pd1 = da.AbstractPartitionDiagram(pd, [[1,2],[-1,-2]])
        sage: pd2 = da.AbstractPartitionDiagram(pd, [[1,2],[-1,-2]])
        sage: pd1
        {{-2, -1}, {1, 2}}
        sage: pd1 == pd2
        True
        sage: pd1 == [[1,2],[-1,-2]]
        True
        sage: pd1 == ((-2,-1),(2,1))
        True
        sage: pd1 == SetPartition([[1,2],[-1,-2]])
        True
        sage: pd3 = da.AbstractPartitionDiagram(pd, [[1,-2],[-1,2]])
        sage: pd1 == pd3
        False
        sage: pd4 = da.AbstractPartitionDiagram(pd, [[1,2],[3,4]])
        Traceback (most recent call last):
        ...
        ValueError: {{1, 2}, {3, 4}} does not represent two rows of vertices of order 2
    """

    def __init__(self, parent, d, check=True):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: import sage.combinat.diagram_algebras as da
            sage: pd = da.AbstractPartitionDiagrams(2)
            sage: pd1 = da.AbstractPartitionDiagram(pd, ((-2,-1),(1,2)) )
        """
        self._base_diagram = tuple(sorted(tuple(sorted(i)) for i in d))
        super().__init__(parent, self._base_diagram, check=check)

    def check(self):
        r"""
        Check the validity of the input for the diagram.

        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: pd = da.AbstractPartitionDiagrams(2)
            sage: pd1 = da.AbstractPartitionDiagram(pd, [[1,2],[-1,-2]]) # indirect doctest
            sage: pd2 = da.AbstractPartitionDiagram(pd, [[1,2],[3,4]]) # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: {{1, 2}, {3, 4}} does not represent two rows of vertices of order 2
            sage: pd2 = da.AbstractPartitionDiagram(pd, [[1],[-1]]) # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: {{-1}, {1}} does not represent two rows of vertices of order 2

            sage: pd2 = da.AbstractPartitionDiagram(pd, [[[1,2],[-1,-2]]]) # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: unhashable type: 'list'
        """
        if self._base_diagram:
            tst = frozenset(e for B in self._base_diagram for e in B)
            if tst != self.parent()._set:
                raise ValueError("{} does not represent two rows of vertices of order {}".format(
                                     self, self.parent().order))

    def __hash__(self):
        """
        Return the hash of ``self``.

        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: pd = da.AbstractPartitionDiagrams(2)
            sage: pd1 = da.AbstractPartitionDiagram(pd, [[1,2],[-1,-2]])
            sage: pd2 = da.AbstractPartitionDiagram(pd, [[1,2],[-1,-2]])
            sage: hash(pd1) == hash(pd2)
            True
            sage: hash(pd1) == hash( ((-2,-1), (1,2)) )
            True
        """
        return hash(self._base_diagram)

    def __eq__(self, other):
        r"""
        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: pd = da.AbstractPartitionDiagrams(2)
            sage: pd1 = da.AbstractPartitionDiagram(pd, [[1,2],[-1,-2]])
            sage: pd2 = da.AbstractPartitionDiagram(pd, [[1,2],[-1,-2]])
            sage: pd1 == pd2
            True
            sage: pd1 == [[1,2],[-1,-2]]
            True
            sage: pd1 == ((-2,-1),(2,1))
            True
            sage: pd1 == SetPartition([[1,2],[-1,-2]])
            True
            sage: pd3 = da.AbstractPartitionDiagram(pd, [[1,-2],[-1,2]])
            sage: pd1 == pd3
            False

        Check the inherited inequality::

            sage: pd1 = da.AbstractPartitionDiagram(pd, [[1,2],[-1,-2]])
            sage: pd2 = da.AbstractPartitionDiagram(pd, [[1,-2],[-1,2]])
            sage: pd1 != pd2
            True
            sage: pd1 != ((-2,-1),(2,1))
            False
        """
        try:
            return self._base_diagram == other._base_diagram
        except AttributeError:
            pass

        try:
            other2 = self.parent(other)
            return self._base_diagram == other2._base_diagram
        except (TypeError, ValueError, AttributeError):
            return False

    def __lt__(self, other):
        """
        Compare less than.

        EXAMPLES::

            sage: import sage.combinat.diagram_algebras as da
            sage: pd = da.AbstractPartitionDiagrams(2)
            sage: pd1 = da.AbstractPartitionDiagram(pd, [[1,2],[-1,-2]])
            sage: pd2 = da.AbstractPartitionDiagram(pd, [[1,-2],[-1,2]])
            sage: pd1 < pd2
            True
            sage: pd2 < pd1
            False
            sage: pd2 > pd1
            True
            sage: pd1 > pd2
            False
        """
        if not isinstance(other, AbstractPartitionDiagram):
            return False
        return self._base_diagram < other._base_diagram

    def base_diagram(self):
        r"""
        Return the underlying implementation of the diagram.

        OUTPUT: tuple of tuples of integers

        EXAMPLES::

            sage: import sage.combinat.diagram_algebras as da
            sage: pd = da.AbstractPartitionDiagrams(2)
            sage: pd([[1,2],[-1,-2]]).base_diagram() == ((-2,-1),(1,2))
            True
        """
        return self._base_diagram # note, this works because self._base_diagram is immutable

    diagram = base_diagram

    def set_partition(self):
        r"""
        Return the underlying implementation of the diagram as a set of sets.

        EXAMPLES::

            sage: import sage.combinat.diagram_algebras as da
            sage: pd = da.AbstractPartitionDiagrams(2)
            sage: X = pd([[1,2],[-1,-2]]).set_partition(); X
            {{-2, -1}, {1, 2}}
            sage: X.parent()
            Set partitions
        """
        return SetPartitions()(self)

    def compose(self, other, check=True):
        r"""
        Compose ``self`` with ``other``.

        The composition of two diagrams `X` and `Y` is given by placing
        `X` on top of `Y` and removing all loops.

        OUTPUT:

        A tuple where the first entry is the composite diagram and the
        second entry is how many loop were removed.

        .. NOTE::

            This is not really meant to be called directly, but it works
            to call it this way if desired.

        EXAMPLES::

            sage: import sage.combinat.diagram_algebras as da
            sage: pd = da.AbstractPartitionDiagrams(2)
            sage: pd([[1,2],[-1,-2]]).compose(pd([[1,2],[-1,-2]]))
            ({{-2, -1}, {1, 2}}, 1)
        """
        (composite_diagram, loops_removed) = set_partition_composition(self._base_diagram, other._base_diagram)
        return (self.__class__(self.parent(), composite_diagram, check=check), loops_removed)

    def propagating_number(self):
        r"""
        Return the propagating number of the diagram.

        The propagating number is the number of blocks with both a
        positive and negative number.

        EXAMPLES::

            sage: import sage.combinat.diagram_algebras as da
            sage: pd = da.AbstractPartitionDiagrams(2)
            sage: d1 = pd([[1,-2],[2,-1]])
            sage: d1.propagating_number()
            2
            sage: d2 = pd([[1,2],[-2,-1]])
            sage: d2.propagating_number()
            0
        """
        return ZZ(sum(1 for part in self._base_diagram if min(part) < 0 and max(part) > 0))

    def count_blocks_of_size(self, n):
        r"""
        Count the number of blocks of a given size.

        INPUT:

        - ``n`` -- positive integer

        EXAMPLES::

            sage: from sage.combinat.diagram_algebras import PartitionDiagram
            sage: pd = PartitionDiagram([[1,-3,-5],[2,4],[3,-1,-2],[5],[-4]])
            sage: pd.count_blocks_of_size(1)
            2
            sage: pd.count_blocks_of_size(2)
            1
            sage: pd.count_blocks_of_size(3)
            2
        """
        return sum(ZZ(len(block) == n) for block in self)

    def order(self):
        r"""
        Return the maximum entry in the diagram element.

        A diagram element will be a partition of the set
        `\{-1, -2, \ldots, -k, 1, 2, \ldots, k\}`.  The order of
        the diagram element is the value `k`.

        EXAMPLES::

            sage: from sage.combinat.diagram_algebras import PartitionDiagram
            sage: PartitionDiagram([[1,-1],[2,-2,-3],[3]]).order()
            3
            sage: PartitionDiagram([[1,-1]]).order()
            1
            sage: PartitionDiagram([[1,-3,-5],[2,4],[3,-1,-2],[5],[-4]]).order()
            5
        """
        return self.parent().order

    def is_planar(self):
        r"""
        Test if the diagram ``self`` is planar.

        A diagram element is planar if the graph of the nodes is planar.

        EXAMPLES::

            sage: from sage.combinat.diagram_algebras import BrauerDiagram
            sage: BrauerDiagram([[1,-2],[2,-1]]).is_planar()
            False
            sage: BrauerDiagram([[1,-1],[2,-2]]).is_planar()
            True
        """
        return is_planar(self)

    def dual(self):
        """
        Return the dual diagram of ``self`` by flipping it top-to-bottom.

        EXAMPLES::

            sage: from sage.combinat.diagram_algebras import PartitionDiagram
            sage: D = PartitionDiagram([[1,-1],[2,-2,-3],[3]])
            sage: D.dual()
            {{-3}, {-2, 2, 3}, {-1, 1}}
        """
        return self.parent([[-i for i in part] for part in self])


class IdealDiagram(AbstractPartitionDiagram):
    r"""
    The element class for a ideal diagram.

    An ideal diagram for an integer `k` is a partition of the set
    `\{1, \ldots, k, -1, \ldots, -k\}` where the propagating number is
    strictly smaller than the order.

    EXAMPLES::

        sage: from sage.combinat.diagram_algebras import IdealDiagrams as IDs
        sage: IDs(2)
        Ideal diagrams of order 2
        sage: IDs(2).list()
        [{{-2, -1, 1, 2}},
         {{-2, 1, 2}, {-1}},
         {{-2}, {-1, 1, 2}},
         {{-2, -1}, {1, 2}},
         {{-2}, {-1}, {1, 2}},
         {{-2, -1, 1}, {2}},
         {{-2, 1}, {-1}, {2}},
         {{-2, -1, 2}, {1}},
         {{-2, 2}, {-1}, {1}},
         {{-2}, {-1, 1}, {2}},
         {{-2}, {-1, 2}, {1}},
         {{-2, -1}, {1}, {2}},
         {{-2}, {-1}, {1}, {2}}]

        sage: from sage.combinat.diagram_algebras import PartitionDiagrams as PDs
        sage: PDs(4).cardinality() == factorial(4) + IDs(4).cardinality()
        True
    """
    @staticmethod
    def __classcall_private__(cls, diag):
        """
        Normalize input to initialize diagram.

        The order of the diagram element is the maximum value found in
        the list of lists.

        EXAMPLES::

            sage: from sage.combinat.diagram_algebras import IdealDiagram
            sage: IdealDiagram([[1],[-1]])
            {{-1}, {1}}
            sage: IdealDiagram([[1], [-1]]).parent()
            Ideal diagrams of order 1
        """
        order = max(v for p in diag for v in p)
        return IdealDiagrams(order)(diag)

    def check(self):
        r"""
        Check the validity of the input for ``self``.

        TESTS::

            sage: from sage.combinat.diagram_algebras import IdealDiagram
            sage: pd1 = IdealDiagram([[1,2],[-1,-2]])  # indirect doctest
            sage: pd2 = IdealDiagram([[1,-2],[2,-1]])  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: the diagram {{-2, 1}, {-1, 2}} must have a propagating number smaller than the order
            sage: pd3 = IdealDiagram([[1,2,-1,-3]])    # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: {{-3, -1, 1, 2}} does not represent two rows of vertices of order 2
            sage: pd4 = IdealDiagram([[1,-2,-1],[2]])  # indirect doctest
        """
        super().check()
        if self.propagating_number() >= self.order():
            raise ValueError("the diagram %s must have a propagating number smaller than the order" % self)


class PlanarDiagram(AbstractPartitionDiagram):
    r"""
    The element class for a planar diagram.

    A planar diagram for an integer `k` is a partition of the set
    `\{1, \ldots, k, -1, \ldots, -k\}` so that the diagram is non-crossing.

    EXAMPLES::

        sage: from sage.combinat.diagram_algebras import PlanarDiagrams
        sage: PlanarDiagrams(2)
        Planar diagrams of order 2
        sage: PlanarDiagrams(2).list()
        [{{-2}, {-1}, {1, 2}},
         {{-2}, {-1}, {1}, {2}},
         {{-2, 1}, {-1}, {2}},
         {{-2, 2}, {-1}, {1}},
         {{-2, 1, 2}, {-1}},
         {{-2, 2}, {-1, 1}},
         {{-2}, {-1, 1}, {2}},
         {{-2}, {-1, 2}, {1}},
         {{-2}, {-1, 1, 2}},
         {{-2, -1}, {1, 2}},
         {{-2, -1}, {1}, {2}},
         {{-2, -1, 1}, {2}},
         {{-2, -1, 2}, {1}},
         {{-2, -1, 1, 2}}]
    """
    @staticmethod
    def __classcall_private__(cls, diag):
        """
        Normalize input to initialize diagram.

        The order of the diagram element is the maximum value found in
        the list of lists.

        EXAMPLES::

            sage: from sage.combinat.diagram_algebras import PlanarDiagram
            sage: PlanarDiagram([[1,-1]])
            {{-1, 1}}
            sage: PlanarDiagram([[1, -1]]).parent()
            Planar diagrams of order 1
        """
        order = max(v for p in diag for v in p)
        PD = PlanarDiagrams(order)
        return PD(diag)

    def check(self):
        r"""
        Check the validity of the input for ``self``.

        TESTS::

            sage: from sage.combinat.diagram_algebras import PlanarDiagram
            sage: pd1 = PlanarDiagram([[1,2],[-1,-2]])  # indirect doctest
            sage: pd2 = PlanarDiagram([[1,-2],[2,-1]])  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: the diagram {{-2, 1}, {-1, 2}} must be planar
            sage: pd3 = PlanarDiagram([[1,2,-1,-3]])    # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: {{-3, -1, 1, 2}} does not represent two rows of vertices of order 2
            sage: pd4 = PlanarDiagram([[1,-2,-1],[2]])  # indirect doctest
        """
        super().check()
        if not self.is_planar():
            raise ValueError("the diagram %s must be planar" % self)


class TemperleyLiebDiagram(AbstractPartitionDiagram):
    r"""
    The element class for a Temperley-Lieb diagram.

    A Temperley-Lieb diagram for an integer `k` is a partition of the set
    `\{1, \ldots, k, -1, \ldots, -k\}` so that the blocks are all of size
    2 and the diagram is planar.

    EXAMPLES::

        sage: from sage.combinat.diagram_algebras import TemperleyLiebDiagrams
        sage: TemperleyLiebDiagrams(2)
        Temperley Lieb diagrams of order 2
        sage: TemperleyLiebDiagrams(2).list()
        [{{-2, -1}, {1, 2}}, {{-2, 2}, {-1, 1}}]
    """
    @staticmethod
    def __classcall_private__(cls, diag):
        """
        Normalize input to initialize diagram.

        The order of the diagram element is the maximum value found in
        the list of lists.

        EXAMPLES::

            sage: from sage.combinat.diagram_algebras import TemperleyLiebDiagram
            sage: TemperleyLiebDiagram([[1,-1]])
            {{-1, 1}}
            sage: TemperleyLiebDiagram([[1, -1]]).parent()
            Temperley Lieb diagrams of order 1
        """
        order = max(v for p in diag for v in p)
        TLD = TemperleyLiebDiagrams(order)
        return TLD(diag)

    def check(self):
        r"""
        Check the validity of the input for ``self``.

        TESTS::

            sage: from sage.combinat.diagram_algebras import TemperleyLiebDiagram
            sage: pd1 = TemperleyLiebDiagram([[1,2],[-1,-2]])  # indirect doctest
            sage: pd2 = TemperleyLiebDiagram([[1,-2],[2,-1]])  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: the diagram {{-2, 1}, {-1, 2}} must be planar
            sage: pd3 = TemperleyLiebDiagram([[1,2,-1,-3]])    # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: {{-3, -1, 1, 2}} does not represent two rows of vertices of order 2
            sage: pd4 = TemperleyLiebDiagram([[1,-2,-1],[2]])  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: all blocks of {{-2, -1, 1}, {2}} must be of size 2
        """
        super().check()
        if any(len(block) != 2 for block in self):
            raise ValueError("all blocks of %s must be of size 2" % self )
        if not self.is_planar():
            raise ValueError("the diagram %s must be planar" % self)


class PartitionDiagram(AbstractPartitionDiagram):
    r"""
    The element class for a partition diagram.

    A partition diagram for an integer `k` is a partition of the set
    `\{1, \ldots, k, -1, \ldots, -k\}`

    EXAMPLES::

        sage: from sage.combinat.diagram_algebras import PartitionDiagram, PartitionDiagrams
        sage: PartitionDiagrams(1)
        Partition diagrams of order 1
        sage: PartitionDiagrams(1).list()
        [{{-1, 1}}, {{-1}, {1}}]
        sage: PartitionDiagram([[1,-1]])
        {{-1, 1}}
        sage: PartitionDiagram(((1,-2),(2,-1))).parent()
        Partition diagrams of order 2
    """
    @staticmethod
    def __classcall_private__(cls, diag):
        """
        Normalize input to initialize diagram.

        The order of the diagram element is the maximum value found in
        the list of lists.

        EXAMPLES::

            sage: from sage.combinat.diagram_algebras import PartitionDiagram
            sage: PartitionDiagram([[1],[-1]])
            {{-1}, {1}}
            sage: PartitionDiagram([[1],[-1]]).parent()
            Partition diagrams of order 1
        """
        order = max(v for p in diag for v in p)
        PD = PartitionDiagrams(order)
        return PD(diag)


class BrauerDiagram(AbstractPartitionDiagram):
    r"""
    A Brauer diagram.

    A Brauer diagram for an integer `k` is a partition of the set
    `\{1, \ldots, k, -1, \ldots, -k\}` with block size 2.

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: bd = da.BrauerDiagrams(2)
        sage: bd1 = bd([[1,2],[-1,-2]])
        sage: bd2 = bd([[1,2,-1,-2]])
        Traceback (most recent call last):
        ...
        ValueError: all blocks of {{-2, -1, 1, 2}} must be of size 2

    TESTS::

        sage: import sage.combinat.diagram_algebras as da
        sage: bd = da.BrauerDiagrams(2)( ((-2,-1),(1,2)) )
        sage: TestSuite(bd).run()
    """
    @staticmethod
    def __classcall_private__(cls, diag):
        """
        Normalize input to initialize diagram.

        The order of the diagram element is the maximum value found in
        the list of lists.

        EXAMPLES::

            sage: from sage.combinat.diagram_algebras import BrauerDiagram
            sage: bd = BrauerDiagram([[1,-1]]); bd
            {{-1, 1}}
            sage: bd.parent()
            Brauer diagrams of order 1
        """
        order = max(v for p in diag for v in p)
        BD = BrauerDiagrams(order)
        return BD(diag)

    def check(self):
        r"""
        Check the validity of the input for ``self``.

        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: bd = da.BrauerDiagrams(2)
            sage: bd1 = bd([[1,2],[-1,-2]])  # indirect doctest
            sage: bd2 = bd([[1,2,-1,-2]])    # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: all blocks of {{-2, -1, 1, 2}} must be of size 2
        """
        super().check()
        if any(len(i) != 2 for i in self):
            raise ValueError("all blocks of %s must be of size 2" % self)

    def _repr_(self):
        r"""
        Return a string representation of a Brauer diagram.

        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: bd = da.BrauerDiagrams(2)
            sage: bd1 = bd([[1,2],[-1,-2]]); bd1
            {{-2, -1}, {1, 2}}
        """
        return self.parent().options._dispatch(self, '_repr_', 'display')

    # add options to class
    class options(GlobalOptions):
        r"""
        Set and display the global options for Brauer diagram (algebras). If no
        parameters are set, then the function returns a copy of the options
        dictionary.

        The ``options`` to diagram algebras can be accessed as the method
        :obj:`BrauerAlgebra.options` of :class:`BrauerAlgebra` and
        related classes.

        @OPTIONS@

        The compact representation ``[A/B;pi]`` of the Brauer algebra diagram
        (see [GL1996]_) has the following components:

        - ``A`` -- list of pairs of positive elements (upper row) that
          are connected

        - ``B`` -- list of pairs of negative elements (lower row) that
          are connected

        - ``pi`` -- a permutation that is to be interpreted as the relative
          order of the remaining elements in the top row and the bottom row

        EXAMPLES::

            sage: R.<q> = QQ[]
            sage: BA = BrauerAlgebra(2, q)
            sage: E = BA([[1,2],[-1,-2]])
            sage: E
            B{{-2, -1}, {1, 2}}
            sage: BA8 = BrauerAlgebra(8, q)
            sage: BA8([[1,-4],[2,4],[3,8],[-7,-2],[5,7],[6,-1],[-3,-5],[-6,-8]])
            B{{-8, -6}, {-7, -2}, {-5, -3}, {-4, 1}, {-1, 6}, {2, 4}, {3, 8}, {5, 7}}
            sage: BrauerAlgebra.options.display = "compact"
            sage: E
            B[12/12;]
            sage: BA8([[1,-4],[2,4],[3,8],[-7,-2],[5,7],[6,-1],[-3,-5],[-6,-8]])
            B[24.38.57/35.27.68;21]
            sage: BrauerAlgebra.options._reset()
        """
        NAME = 'Brauer diagram'
        module = 'sage.combinat.diagram_algebras'
        option_class = 'BrauerDiagram'
        display = dict(default='normal',
                       description='Specifies how the Brauer diagrams should be printed',
                       values=dict(normal="Using the normal representation",
                                   compact="Using the compact representation"),
                       case_sensitive=False)

    def _repr_normal(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: import sage.combinat.diagram_algebras as da
            sage: bd = da.BrauerDiagrams(2)
            sage: bd([[1,2],[-1,-2]])._repr_normal()
            '{{-2, -1}, {1, 2}}'
        """
        return super()._repr_()

    def _repr_compact(self):
        """
        Return a compact string representation of ``self``.

        EXAMPLES::

            sage: import sage.combinat.diagram_algebras as da
            sage: bd = da.BrauerDiagrams(2)
            sage: bd([[1,2],[-1,-2]])._repr_compact()
            '[12/12;]'
            sage: bd([[1,-2],[2,-1]])._repr_compact()
            '[/;21]'
            sage: bd = da.BrauerDiagrams(7)
            sage: bd([[1,4],[6,7], [-2,-6],[-5,-7], [2,-4],[3,-1],[5,-3]])._repr_compact()
            '[14.67/26.57;312]'
        """
        (top, bot, thru) = self.involution_permutation_triple()
        bot.reverse()
        s1 = ".".join("".join(str(b) for b in block) for block in top)
        s2 = ".".join("".join(str(abs(k)) for k in sorted(block, reverse=True))
                      for block in bot)
        s3 = "".join(str(x) for x in thru)
        return "[{}/{};{}]".format(s1, s2, s3)

    def involution_permutation_triple(self, curt=True):
        r"""
        Return the involution permutation triple of ``self``.

        From Graham-Lehrer (see :class:`BrauerDiagrams`), a Brauer diagram
        is a triple `(D_1, D_2, \pi)`, where:

        - `D_1` is a partition of the top nodes;
        - `D_2` is a partition of the bottom nodes;
        - `\pi` is the induced permutation on the free nodes.

        INPUT:

        - ``curt`` -- boolean (default: ``True``); if ``True``, then return
          bijection on free nodes as a one-line notation (standardized to look
          like a permutation), else, return the honest mapping, a list of pairs
          `(i, -j)` describing the bijection on free nodes

        EXAMPLES::

            sage: import sage.combinat.diagram_algebras as da
            sage: bd = da.BrauerDiagrams(3)
            sage: elm = bd([[1,2],[-2,-3],[3,-1]])
            sage: elm.involution_permutation_triple()
            ([(1, 2)], [(-3, -2)], [1])
            sage: elm.involution_permutation_triple(curt=False)
            ([(1, 2)], [(-3, -2)], [[3, -1]])
        """
        top = []
        bottom = []
        for v in self.diagram():
            if min(v) > 0:
                top += [v]
            if max(v) < 0:
                bottom += [v]
        if curt:
            perm = self.perm()
        else:
            perm = self.bijection_on_free_nodes()
        return (top, bottom, perm)

    def bijection_on_free_nodes(self, two_line=False):
        r"""
        Return the induced bijection - as a list of `(x,f(x))` values -
        from the free nodes on the top at the Brauer diagram to the free
        nodes at the bottom of ``self``.

        OUTPUT:

        If ``two_line`` is ``True``, then the output is the induced
        bijection as a two-row list ``(inputs, outputs)``.

        EXAMPLES::

            sage: import sage.combinat.diagram_algebras as da
            sage: bd = da.BrauerDiagrams(3)
            sage: elm = bd([[1,2],[-2,-3],[3,-1]])
            sage: elm.bijection_on_free_nodes()
            [[3, -1]]
            sage: elm2 = bd([[1,-2],[2,-3],[3,-1]])
            sage: elm2.bijection_on_free_nodes(two_line=True)
            [[1, 2, 3], [-2, -3, -1]]
        """
        terms = sorted(sorted(v, reverse=True) for v in self.diagram()
                       if max(v) > 0 and min(v) < 0)
        if two_line:
            terms = [[t[i] for t in terms] for i in range(2)]
        return terms

    def perm(self):
        r"""
        Return the induced bijection on the free nodes of ``self`` in
        one-line notation, re-indexed and treated as a permutation.

        .. SEEALSO::

            :meth:`bijection_on_free_nodes`

        EXAMPLES::

            sage: import sage.combinat.diagram_algebras as da
            sage: bd = da.BrauerDiagrams(3)
            sage: elm = bd([[1,2],[-2,-3],[3,-1]])
            sage: elm.perm()
            [1]
        """
        long_form = self.bijection_on_free_nodes()
        if not long_form:
            return long_form

        short_form = [abs(v[1]) for v in long_form]
        # given any list [i1,i2,...,ir] with distinct positive integer entries,
        # return naturally associated permutation of [r].
        # probably already defined somewhere in Permutations/Compositions/list/etc.
        std = list(range(1, len(short_form) + 1))
        j = 0
        for i in range(max(short_form)+1):
            if i in short_form:
                j += 1
                std[short_form.index(i)] = j
        return std

    def is_elementary_symmetric(self):
        r"""
        Check if is elementary symmetric.

        Let `(D_1, D_2, \pi)` be the Graham-Lehrer representation
        of the Brauer diagram `d`. We say `d` is *elementary symmetric*
        if `D_1 = D_2` and `\pi` is the identity.

        EXAMPLES::

            sage: import sage.combinat.diagram_algebras as da
            sage: bd = da.BrauerDiagrams(3)
            sage: elm = bd([[1,2],[-1,-2],[3,-3]])
            sage: elm.is_elementary_symmetric()
            True
            sage: elm2 = bd([[1,2],[-1,-3],[3,-2]])
            sage: elm2.is_elementary_symmetric()
            False
        """
        (D1,D2,pi) = self.involution_permutation_triple()
        D1 = sorted(sorted(abs(y) for y in x) for x in D1)
        D2 = sorted(sorted(abs(y) for y in x) for x in D2)
        return D1 == D2 and pi == list(range(1,len(pi)+1))


class AbstractPartitionDiagrams(Parent, UniqueRepresentation):
    r"""
    This is an abstract base class for partition diagrams.

    The primary use of this class is to serve as basis keys for
    diagram algebras, but diagrams also have properties in their
    own right. Furthermore, this class is meant to be extended to
    create more efficient contains methods.

    INPUT:

    - ``order`` -- integer or integer `+ 1/2`; the order of the diagrams
    - ``category`` -- (default: ``FiniteEnumeratedSets()``) the category

    All concrete classes should implement attributes

    - ``_name`` -- the name of the class
    - ``_diagram_func`` -- an iterator function that takes the order
      as its only input

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: pd = da.PartitionDiagrams(2)
        sage: pd
        Partition diagrams of order 2
        sage: pd.an_element() in pd
        True
        sage: elm = pd([[1,2],[-1,-2]])
        sage: elm in pd
        True
    """
    Element = AbstractPartitionDiagram

    def __init__(self, order, category=None):
        r"""
        Initialize ``self``.

        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: pd = da.AbstractPartitionDiagrams(2)
            sage: pd.category()
            Category of finite enumerated sets
            sage: pd = da.AbstractPartitionDiagrams(2, Sets().Finite())
            sage: pd.category()
            Category of finite sets

            sage: pd = da.PartitionDiagrams(2)
            sage: TestSuite(pd).run()

            sage: bd = da.BrauerDiagrams(2)
            sage: TestSuite(bd).run()

            sage: td = da.TemperleyLiebDiagrams(2)
            sage: TestSuite(td).run()

            sage: pld = da.PlanarDiagrams(2)
            sage: TestSuite(pld).run()

            sage: id = da.IdealDiagrams(2)
            sage: TestSuite(id).run()
        """
        if category is None:
            category = FiniteEnumeratedSets()
        Parent.__init__(self, category=category)
        if order in ZZ:
            self.order = ZZ(order)
            base_set = frozenset(list(range(1,order+1)) + list(range(-order,0)))
        else:
            # order is a half-integer.
            self.order = QQ(order)
            base_set = frozenset(list(range(1,ZZ(ZZ(1)/ZZ(2) + order)+1))
                                 + list(range(ZZ(-ZZ(1)/ZZ(2) - order),0)))
        self._set = base_set

    def _repr_(self):
        r"""
        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: da.PartitionDiagrams(2)
            Partition diagrams of order 2
        """
        return "{} diagrams of order {}".format(self._name, self.order)

    def __iter__(self):
        r"""
        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: list(da.PartitionDiagrams(2))
            [{{-2, -1, 1, 2}},
             {{-2, 1, 2}, {-1}},
             {{-2}, {-1, 1, 2}},
             {{-2, -1}, {1, 2}},
             {{-2}, {-1}, {1, 2}},
             {{-2, -1, 1}, {2}},
             {{-2, 1}, {-1, 2}},
             {{-2, 1}, {-1}, {2}},
             {{-2, 2}, {-1, 1}},
             {{-2, -1, 2}, {1}},
             {{-2, 2}, {-1}, {1}},
             {{-2}, {-1, 1}, {2}},
             {{-2}, {-1, 2}, {1}},
             {{-2, -1}, {1}, {2}},
             {{-2}, {-1}, {1}, {2}}]

            sage: list(da.PartitionDiagrams(3/2))
            [{{-2, -1, 1, 2}},
             {{-2, 1, 2}, {-1}},
             {{-2, 2}, {-1, 1}},
             {{-2, -1, 2}, {1}},
             {{-2, 2}, {-1}, {1}}]

            sage: list(da.BrauerDiagrams(5/2))
            [{{-3, 3}, {-2, -1}, {1, 2}},
             {{-3, 3}, {-2, 1}, {-1, 2}},
             {{-3, 3}, {-2, 2}, {-1, 1}}]
            sage: list(da.BrauerDiagrams(2))
            [{{-2, -1}, {1, 2}}, {{-2, 1}, {-1, 2}}, {{-2, 2}, {-1, 1}}]

            sage: list(da.TemperleyLiebDiagrams(5/2))
            [{{-3, 3}, {-2, -1}, {1, 2}}, {{-3, 3}, {-2, 2}, {-1, 1}}]
            sage: list(da.TemperleyLiebDiagrams(2))
            [{{-2, -1}, {1, 2}}, {{-2, 2}, {-1, 1}}]

            sage: list(da.PlanarDiagrams(3/2))
            [{{-2, 1, 2}, {-1}},
             {{-2, 2}, {-1}, {1}},
             {{-2, 2}, {-1, 1}},
             {{-2, -1, 2}, {1}},
             {{-2, -1, 1, 2}}]

            sage: list(da.PlanarDiagrams(2))
            [{{-2}, {-1}, {1, 2}},
             {{-2}, {-1}, {1}, {2}},
             {{-2, 1}, {-1}, {2}},
             {{-2, 2}, {-1}, {1}},
             {{-2, 1, 2}, {-1}},
             {{-2, 2}, {-1, 1}},
             {{-2}, {-1, 1}, {2}},
             {{-2}, {-1, 2}, {1}},
             {{-2}, {-1, 1, 2}},
             {{-2, -1}, {1, 2}},
             {{-2, -1}, {1}, {2}},
             {{-2, -1, 1}, {2}},
             {{-2, -1, 2}, {1}},
             {{-2, -1, 1, 2}}]

            sage: list(da.IdealDiagrams(3/2))
            [{{-2, -1, 1, 2}},
             {{-2, 1, 2}, {-1}},
             {{-2, -1, 2}, {1}},
             {{-2, 2}, {-1}, {1}}]
            sage: list(da.IdealDiagrams(2))
            [{{-2, -1, 1, 2}},
             {{-2, 1, 2}, {-1}},
             {{-2}, {-1, 1, 2}},
             {{-2, -1}, {1, 2}},
             {{-2}, {-1}, {1, 2}},
             {{-2, -1, 1}, {2}},
             {{-2, 1}, {-1}, {2}},
             {{-2, -1, 2}, {1}},
             {{-2, 2}, {-1}, {1}},
             {{-2}, {-1, 1}, {2}},
             {{-2}, {-1, 2}, {1}},
             {{-2, -1}, {1}, {2}},
             {{-2}, {-1}, {1}, {2}}]
        """
        # The _diagram_func gets set as a method, but we want to
        #   treat it like an attribute, so we call the underlying
        #   __func__.
        for i in self._diagram_func.__func__(self.order):
            yield self.element_class(self, i, check=False)

    def __contains__(self, obj):
        r"""
        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: pd = da.PartitionDiagrams(2)
            sage: pd.an_element() in pd
            True
            sage: elm = pd([[1,2],[-1,-2]])
            sage: elm in pd # indirect doctest
            True
        """
        if not hasattr(obj, '_base_diagram'):
            try:
                obj = self._element_constructor_(obj)
            except (ValueError, TypeError):
                return False
        if obj.base_diagram():
            tst = sorted(flatten(obj.base_diagram()))
            if len(tst) % 2 or tst != list(range(-len(tst)//2,0)) + list(range(1,len(tst)//2+1)):
                return False
            return True
        return self.order == 0

    def _element_constructor_(self, d):
        r"""
        Construct an element of ``self``.

        EXAMPLES::

            sage: import sage.combinat.diagram_algebras as da
            sage: pd = da.AbstractPartitionDiagrams(2)
            sage: elm = pd( [[1,2], [-1,-2]] ); elm
            {{-2, -1}, {1, 2}}
            sage: pd( [{1,2}, {-1,-2}] ) == elm
            True
            sage: pd( ((1,2), (-1,-2)) ) == elm
            True
            sage: pd( SetPartition([[1,2], [-1,-2]]) ) == elm
            True

            sage: bd = da.BrauerDiagrams(2)
            sage: bd( [[1,2],[-1,-2]] )
            {{-2, -1}, {1, 2}}
        """
        return self.element_class(self, d)


class PartitionDiagrams(AbstractPartitionDiagrams):
    r"""
    This class represents all partition diagrams of integer or integer
    `+ 1/2` order.

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: pd = da.PartitionDiagrams(1); pd
        Partition diagrams of order 1
        sage: pd.list()
        [{{-1, 1}}, {{-1}, {1}}]

        sage: pd = da.PartitionDiagrams(3/2); pd
        Partition diagrams of order 3/2
        sage: pd.list()
        [{{-2, -1, 1, 2}},
         {{-2, 1, 2}, {-1}},
         {{-2, 2}, {-1, 1}},
         {{-2, -1, 2}, {1}},
         {{-2, 2}, {-1}, {1}}]

    TESTS::

        sage: import sage.combinat.diagram_algebras as da
        sage: pd = da.PartitionDiagrams(3)
        sage: pd.an_element() in pd
        True
        sage: pd.cardinality() == len(pd.list())
        True

        sage: pd = da.PartitionDiagrams(5/2)
        sage: pd.an_element() in pd
        True
        sage: pd.cardinality() == len(pd.list())
        True
    """
    Element = PartitionDiagram
    _name = "Partition"
    _diagram_func = partition_diagrams

    def cardinality(self):
        r"""
        The cardinality of partition diagrams of half-integer order `n` is
        the `2n`-th Bell number.

        EXAMPLES::

            sage: import sage.combinat.diagram_algebras as da
            sage: pd = da.PartitionDiagrams(3)
            sage: pd.cardinality()
            203

            sage: pd = da.PartitionDiagrams(7/2)
            sage: pd.cardinality()
            877
        """
        return bell_number(ZZ(2 * self.order))


class BrauerDiagrams(AbstractPartitionDiagrams):
    r"""
    This class represents all Brauer diagrams of integer or integer
    `+1/2` order. For more information on Brauer diagrams,
    see :class:`BrauerAlgebra`.

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: bd = da.BrauerDiagrams(2); bd
        Brauer diagrams of order 2
        sage: bd.list()
        [{{-2, -1}, {1, 2}}, {{-2, 1}, {-1, 2}}, {{-2, 2}, {-1, 1}}]

        sage: bd = da.BrauerDiagrams(5/2); bd
        Brauer diagrams of order 5/2
        sage: bd.list()
        [{{-3, 3}, {-2, -1}, {1, 2}},
         {{-3, 3}, {-2, 1}, {-1, 2}},
         {{-3, 3}, {-2, 2}, {-1, 1}}]

    TESTS::

        sage: import sage.combinat.diagram_algebras as da
        sage: bd = da.BrauerDiagrams(3)
        sage: bd.an_element() in bd
        True
        sage: bd.cardinality() == len(bd.list())
        True

        sage: bd = da.BrauerDiagrams(5/2)
        sage: bd.an_element() in bd
        True
        sage: bd.cardinality() == len(bd.list())
        True

    These diagrams also come equipped with a compact representation based
    on their bipartition triple representation. See the
    :meth:`from_involution_permutation_triple` method for more information.

    ::

        sage: bd = da.BrauerDiagrams(3)
        sage: bd.options.display="compact"
        sage: bd.list()
        [[12/12;1],
         [13/12;1],
         [23/12;1],
         [23/13;1],
         [23/23;1],
         [/;132],
         [/;231],
         [/;321],
         [13/13;1],
         [12/13;1],
         [12/23;1],
         [13/23;1],
         [/;312],
         [/;213],
         [/;123]]
        sage: bd.options._reset()
    """
    Element = BrauerDiagram
    options = BrauerDiagram.options
    _name = "Brauer"
    _diagram_func = brauer_diagrams

    def __contains__(self, obj):
        r"""
        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: bd = da.BrauerDiagrams(2)
            sage: bd.an_element() in bd
            True
            sage: bd([[1,2],[-1,-2]]) in bd
            True
            sage: [[1,2,-1,-2]] in bd
            False
            sage: bd = da.BrauerDiagrams(3/2)
            sage: bd.an_element() in bd
            True
        """
        if self.order in ZZ:
            r = ZZ(self.order)
        else:
            r = ZZ(self.order + ZZ(1)/ZZ(2))
        return super().__contains__(obj) and [len(i) for i in obj] == [2]*r

    def cardinality(self):
        r"""
        Return the cardinality of ``self``.

        The number of Brauer diagrams of integer order `k` is `(2k-1)!!`.

        EXAMPLES::

            sage: import sage.combinat.diagram_algebras as da
            sage: bd = da.BrauerDiagrams(3)
            sage: bd.cardinality()
            15

            sage: bd = da.BrauerDiagrams(7/2)
            sage: bd.cardinality()
            15
        """
        if self.order in ZZ:
            return (2 * ZZ(self.order) - 1).multifactorial(2)
        else:
            return (2 * ZZ(self.order - 1 / 2) - 1).multifactorial(2)

    def symmetric_diagrams(self, l=None, perm=None):
        r"""
        Return the list of Brauer diagrams with symmetric placement of `l` arcs,
        and with free nodes permuted according to `perm`.

        EXAMPLES::

            sage: import sage.combinat.diagram_algebras as da
            sage: bd = da.BrauerDiagrams(4)
            sage: bd.symmetric_diagrams(l=1, perm=[2,1])
            [{{-4, -2}, {-3, 1}, {-1, 3}, {2, 4}},
             {{-4, -3}, {-2, 1}, {-1, 2}, {3, 4}},
             {{-4, -1}, {-3, 2}, {-2, 3}, {1, 4}},
             {{-4, 2}, {-3, -1}, {-2, 4}, {1, 3}},
             {{-4, 3}, {-3, 4}, {-2, -1}, {1, 2}},
             {{-4, 1}, {-3, -2}, {-1, 4}, {2, 3}}]

        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: bd = da.BrauerDiagrams(3/2)
            sage: bd.symmetric_diagrams(l=1, perm=[2,1])
            Traceback (most recent call last):
            ...
            NotImplementedError: only implemented for integer order, not for order 3/2
        """
        # perm = permutation on free nodes
        # l = number of arcs
        if self.order not in ZZ:
            raise NotImplementedError("only implemented for integer order,"
                                      " not for order %s" % (self.order))
        n = ZZ(self.order)
        if l is None:
            l = 0
        if perm is None:
            perm = list(range(1, n+1-2*l))
        out = []
        partition_shape = [2]*l + [1]*(n-2*l)
        for sp in SetPartitions(n, partition_shape):
            sp0 = [block for block in sp if len(block) == 2]
            diag = self.from_involution_permutation_triple((sp0,sp0,perm))
            out.append(diag)
        return out

    def from_involution_permutation_triple(self, D1_D2_pi):
        r"""
        Construct a Brauer diagram of ``self`` from an involution
        permutation triple.

        A Brauer diagram can be represented as a triple where the first
        entry is a list of arcs on the top row of the diagram, the second
        entry is a list of arcs on the bottom row of the diagram, and the
        third entry is a permutation on the remaining nodes. This triple
        is called the *involution permutation triple*. For more
        information, see [GL1996]_.

        INPUT:

        - ``D1_D2_pi`` -- list or tuple where the first entry is a list of
          arcs on the top of the diagram, the second entry is a list of arcs
          on the bottom of the diagram, and the third entry is a permutation
          on the free nodes.

        REFERENCES:

        .. [GL1996] \J.J. Graham and G.I. Lehrer, Cellular algebras.
           Inventiones mathematicae 123 (1996), 1--34.

        EXAMPLES::

            sage: import sage.combinat.diagram_algebras as da
            sage: bd = da.BrauerDiagrams(4)
            sage: bd.from_involution_permutation_triple([[[1,2]],[[3,4]],[2,1]])
            {{-4, -3}, {-2, 3}, {-1, 4}, {1, 2}}

        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: bd = da.BrauerDiagrams(5/2)
            sage: bd.from_involution_permutation_triple([[[1,2]],[[3,4]],[2,1]])
            Traceback (most recent call last):
            ...
            NotImplementedError: only implemented for integer order, not for order 5/2
        """
        if self.order not in ZZ:
            raise NotImplementedError("only implemented for integer order,"
                                      " not for order %s" % (self.order))
        try:
            (D1,D2,pi) = tuple(D1_D2_pi)
        except ValueError:
            raise ValueError("argument %s not in correct form; must be a tuple (D1, D2, pi)" % D1_D2_pi)
        D1 = [[abs(x) for x in b] for b in D1 if len(b) == 2] # not needed if argument correctly passed at outset.
        D2 = [[abs(x) for x in b] for b in D2 if len(b) == 2] # ditto.
        nD2 = [[-i for i in b] for b in D2]
        pi = list(pi)
        nn = set(range(1, self.order+1))
        dom = sorted(nn.difference(flatten([list(x) for x in D1])))
        rng = sorted(nn.difference(flatten([list(x) for x in D2])))
        SP0 = D1 + nD2
        if len(pi) != len(dom) or pi not in Permutations():
            raise ValueError("in the tuple (D1, D2, pi)={}, pi must be a permutation of {} (indicating a permutation on the free nodes of the diagram)".format(
                                    (D1,D2,pi), self.order-2*len(D1)))
        Perm = [[dom[i], -rng[val-1]] for i,val in enumerate(pi)]
        SP = SP0 + Perm
        return self(SP) # could pass 'SetPartition' ?


class TemperleyLiebDiagrams(AbstractPartitionDiagrams):
    r"""
    All Temperley-Lieb diagrams of integer or integer `+1/2` order.

    For more information on Temperley-Lieb diagrams, see
    :class:`TemperleyLiebAlgebra`.

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: td = da.TemperleyLiebDiagrams(3); td
        Temperley Lieb diagrams of order 3
        sage: td.list()
        [{{-3, 3}, {-2, -1}, {1, 2}},
         {{-3, 1}, {-2, -1}, {2, 3}},
         {{-3, -2}, {-1, 1}, {2, 3}},
         {{-3, -2}, {-1, 3}, {1, 2}},
         {{-3, 3}, {-2, 2}, {-1, 1}}]

        sage: td = da.TemperleyLiebDiagrams(5/2); td
        Temperley Lieb diagrams of order 5/2
        sage: td.list()
        [{{-3, 3}, {-2, -1}, {1, 2}}, {{-3, 3}, {-2, 2}, {-1, 1}}]

    TESTS::

        sage: import sage.combinat.diagram_algebras as da
        sage: td = da.TemperleyLiebDiagrams(3)
        sage: td.an_element() in td
        True
        sage: td.cardinality() == len(td.list())
        True

        sage: td = da.TemperleyLiebDiagrams(7/2)
        sage: td.an_element() in td
        True
        sage: td.cardinality() == len(td.list())
        True
    """
    Element = TemperleyLiebDiagram
    _name = "Temperley Lieb"
    _diagram_func = temperley_lieb_diagrams

    def cardinality(self):
        r"""
        Return the cardinality of ``self``.

        The number of Temperley--Lieb diagrams of integer order `k` is the
        `k`-th Catalan number.

        EXAMPLES::

            sage: import sage.combinat.diagram_algebras as da
            sage: td = da.TemperleyLiebDiagrams(3)
            sage: td.cardinality()
            5
        """
        if self.order in ZZ:
            return catalan_number(ZZ(self.order))
        else:
            return catalan_number(ZZ(self.order - 1/2))

    def __contains__(self, obj):
        r"""
        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: td = da.TemperleyLiebDiagrams(2)
            sage: td.an_element() in td
            True
            sage: td([[1,2],[-1,-2]]) in td
            True
            sage: [[1,2],[-1,-2]] in td
            True
            sage: [[1,-2],[-1,2]] in td
            False
        """
        if not hasattr(obj, '_base_diagram'):
            try:
                obj = self._element_constructor_(obj)
            except (ValueError, TypeError):
                return False
        return obj in BrauerDiagrams(self.order) and obj.is_planar()


class PlanarDiagrams(AbstractPartitionDiagrams):
    r"""
    All planar diagrams of integer or integer `+1/2` order.

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: pld = da.PlanarDiagrams(1); pld
        Planar diagrams of order 1
        sage: pld.list()
        [{{-1, 1}}, {{-1}, {1}}]

        sage: pld = da.PlanarDiagrams(3/2); pld
        Planar diagrams of order 3/2
        sage: pld.list()
        [{{-2, 1, 2}, {-1}},
         {{-2, 2}, {-1}, {1}},
         {{-2, 2}, {-1, 1}},
         {{-2, -1, 2}, {1}},
         {{-2, -1, 1, 2}}]

    TESTS::

        sage: import sage.combinat.diagram_algebras as da
        sage: pld = da.PlanarDiagrams(3)
        sage: pld.an_element() in pld
        True
        sage: pld.cardinality() == len(pld.list())
        True
        sage: pld = da.PlanarDiagrams(5/2)
        sage: pld.an_element() in pld
        True
        sage: pld.cardinality() == len(pld.list())
        True
    """
    Element = PlanarDiagram
    _name = "Planar"
    _diagram_func = planar_diagrams

    def cardinality(self):
        r"""
        Return the cardinality of ``self``.

        The number of all planar diagrams of order `k` is the
        `2k`-th Catalan number.

        EXAMPLES::

            sage: import sage.combinat.diagram_algebras as da
            sage: pld = da.PlanarDiagrams(3)
            sage: pld.cardinality()
            132
        """
        return catalan_number(2*self.order)

    def __contains__(self, obj):
        r"""
        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: pld = da.PlanarDiagrams(2)
            sage: pld.an_element() in pld
            True
            sage: pld([[1,2],[-1,-2]]) in pld
            True
            sage: [[1,2],[-1,-2]] in pld
            True
            sage: [[1,-2],[-1,2]] in pld
            False
        """
        if not hasattr(obj, '_base_diagram'):
            try:
                obj = self._element_constructor_(obj)
            except (ValueError, TypeError):
                return False
        return super().__contains__(obj)


class IdealDiagrams(AbstractPartitionDiagrams):
    r"""
    All "ideal" diagrams of integer or integer `+1/2` order.

    If `k` is an integer then an ideal diagram of order `k` is a partition
    diagram of order `k` with propagating number less than `k`.

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: id = da.IdealDiagrams(3)
        sage: id.an_element() in id
        True
        sage: id.cardinality() == len(id.list())
        True
        sage: da.IdealDiagrams(3/2).list()
        [{{-2, -1, 1, 2}},
         {{-2, 1, 2}, {-1}},
         {{-2, -1, 2}, {1}},
         {{-2, 2}, {-1}, {1}}]
    """
    Element = IdealDiagram
    _name = "Ideal"
    _diagram_func = ideal_diagrams

    def __contains__(self, obj):
        r"""
        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: id = da.IdealDiagrams(2)
            sage: id.an_element() in id
            True
            sage: id([[1,2],[-1,-2]]) in id
            True
            sage: [[1,2],[-1,-2]] in id
            True
            sage: [[1,-2],[-1,2]] in id
            False
        """
        if not hasattr(obj, '_base_diagram'):
            try:
                obj = self._element_constructor_(obj)
            except (ValueError, TypeError):
                return False
        return super().__contains__(obj) and obj.propagating_number() < self.order


class DiagramAlgebra(CombinatorialFreeModule):
    r"""
    Abstract class for diagram algebras and is not designed to be used
    directly.

    TESTS::

        sage: import sage.combinat.diagram_algebras as da
        sage: R.<x> = QQ[]
        sage: D = da.DiagramAlgebra(2, x, R, 'P', da.PartitionDiagrams(2))
        sage: list(D.basis())
        [P{{-2, -1, 1, 2}},
         P{{-2, 1, 2}, {-1}},
         P{{-2}, {-1, 1, 2}},
         P{{-2, -1}, {1, 2}},
         P{{-2}, {-1}, {1, 2}},
         P{{-2, -1, 1}, {2}},
         P{{-2, 1}, {-1, 2}},
         P{{-2, 1}, {-1}, {2}},
         P{{-2, 2}, {-1, 1}},
         P{{-2, -1, 2}, {1}},
         P{{-2, 2}, {-1}, {1}},
         P{{-2}, {-1, 1}, {2}},
         P{{-2}, {-1, 2}, {1}},
         P{{-2, -1}, {1}, {2}},
         P{{-2}, {-1}, {1}, {2}}]
    """

    def __init__(self, k, q, base_ring, prefix, diagrams, category=None):
        r"""
        Initialize ``self``.

        INPUT:

        - ``k`` -- the rank
        - ``q`` -- the deformation parameter
        - ``base_ring`` -- the base ring
        - ``prefix`` -- the prefix of our monomials
        - ``diagrams`` -- the object representing all the diagrams
          (i.e. indices for the basis elements)

        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: R.<x> = QQ[]
            sage: D = da.DiagramBasis(2, x, R, 'P', da.PartitionDiagrams(2))
            sage: TestSuite(D).run()
        """
        self._prefix = prefix
        self._q = base_ring(q)
        self._k = k
        self._base_diagrams = diagrams
        cat = AssociativeAlgebras(base_ring.category()).FiniteDimensional().WithBasis()
        if isinstance(self, UnitDiagramMixin):
            cat = cat.Unital()
        category = cat.or_subcategory(category)
        CombinatorialFreeModule.__init__(self, base_ring, diagrams,
                                         category=category, prefix=prefix,
                                         bracket=False)

    def _element_constructor_(self, set_partition):
        r"""
        Construct an element of ``self``.

        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: R.<x> = QQ[]
            sage: D = da.DiagramAlgebra(2, x, R, 'P', da.PartitionDiagrams(2))
            sage: sp = da.to_set_partition( [[1,2], [-1,-2]] )
            sage: b_elt = D(sp); b_elt
            P{{-2, -1}, {1, 2}}
            sage: b_elt in D
            True
            sage: D([[1,2],[-1,-2]]) == b_elt
            True
            sage: D([{1,2},{-1,-2}]) == b_elt
            True
            sage: S = SymmetricGroupAlgebra(R,2)
            sage: D(S([2,1]))
            P{{-2, 1}, {-1, 2}}
            sage: D2 = da.DiagramAlgebra(2, x, R, 'P', da.PlanarDiagrams(2))
            sage: D2(S([1,2]))
            P{{-2, 2}, {-1, 1}}
            sage: D2(S([2,1]))
            Traceback (most recent call last):
            ...
            ValueError: the diagram {{-2, 1}, {-1, 2}} must be planar
        """
        if self.basis().keys().is_parent_of(set_partition):
            return self.basis()[set_partition]
        if isinstance(set_partition, SymmetricGroupAlgebra_n.Element):
            return self._apply_module_morphism(set_partition, self._perm_to_Blst, self)
        sp = self._base_diagrams(set_partition) # attempt conversion
        if sp in self.basis().keys():
            return self.basis()[sp]

        raise ValueError("invalid input of {0}".format(set_partition))

    def __getitem__(self, d):
        """
        Get the basis item of ``self`` indexed by ``d``.

        EXAMPLES::

            sage: import sage.combinat.diagram_algebras as da
            sage: R.<x> = QQ[]
            sage: D = da.DiagramAlgebra(2, x, R, 'P', da.PartitionDiagrams(2))
            sage: sp = da.PartitionDiagrams(2)( [[1,2], [-1,-2]] )
            sage: D[sp]
            P{{-2, -1}, {1, 2}}
            sage: D[[1,-1,2,-2]]
            P{{-2, -1, 1, 2}}
            sage: D3 = da.DiagramAlgebra(3, x, R, 'P', da.PartitionDiagrams(3))
            sage: da.PartitionDiagrams(3)( [[1,2], [-1,-2]] )
            Traceback (most recent call last):
            ...
            ValueError: {{-2, -1}, {1, 2}} does not represent two rows of vertices of order 3
            sage: D3[sp]
            P{{-3, 3}, {-2, -1}, {1, 2}}
            sage: D3[[1,-1,2,-2]]
            P{{-3, 3}, {-2, -1, 1, 2}}
            sage: D3[[1,2,-2]]
            P{{-3, 3}, {-2, 1, 2}, {-1}}
            sage: P = PartitionAlgebra(3,x)
            sage: P[[1]]
            P{{-3, 3}, {-2, 2}, {-1}, {1}}
        """
        if isinstance(d, (list, tuple)) and all(a in ZZ for a in d):
            d = [d]
        d = self._base_diagrams(to_set_partition(d, self._k))
        if d in self.basis().keys():
            return self.basis()[d]
        raise ValueError("{0} is not an index of a basis element".format(d))

    def _perm_to_Blst(self, w):
        """
        Convert the permutation ``w`` to an element of ``self``.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: S = SymmetricGroupAlgebra(R,2)
            sage: import sage.combinat.diagram_algebras as da
            sage: D2 = da.DiagramAlgebra(2, x, R, 'P', da.PlanarDiagrams(2))
            sage: D2._perm_to_Blst([2,1])
            Traceback (most recent call last):
            ...
            ValueError: the diagram {{-2, 1}, {-1, 2}} must be planar
        """
        # 'perm' is a permutation in one-line notation
        # turns w into an expression suitable for the element constructor.
        u = sorted(w)
        p = [[u[i], -x] for i, x in enumerate(w)]
        if len(u) < self.order():
            p1 = [[j, -j] for j in range(len(u) + 1, self.order() + 1)]
            p.extend(p1)
        return self[p]

    def _diag_to_Blst(self, d):
        r"""
        Return an element of ``self`` from the input ``d``.

        If ``d`` is a partial diagram of `\{1,2,\ldots,k,-1,-2,\ldots,-k\}`
        then the set partition is filled in by adding the parts `\{i,-i\}`
        if possible, and singletons sets for the remaining parts.

        INPUT:

        - ``d`` -- an iterable that behaves like
          :class:`AbstractPartitionDiagram` or :class:`Permutation`

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: PartitionAlgebra(3, x, R)._diag_to_Blst([[1,2], [-3,-1]])
            P{{-3, -1}, {-2}, {1, 2}, {3}}
            sage: BrauerAlgebra(4, x, R)._diag_to_Blst([3,1,2])
            B{{-4, 4}, {-3, 1}, {-2, 3}, {-1, 2}}
            sage: import sage.combinat.diagram_algebras as da
            sage: D3 = da.DiagramAlgebra(3, x, R, 'P', da.PlanarDiagrams(3))
            sage: D3._diag_to_Blst([[1, 2], [-2,-1]])
            P{{-3, 3}, {-2, -1}, {1, 2}}
            sage: D3._diag_to_Blst([[-1,2], [-2,1]])
            Traceback (most recent call last):
            ...
            ValueError: the diagram {{-3, 3}, {-2, 1}, {-1, 2}} must be planar
            sage: D3._diag_to_Blst([[-1,2], [-3,1]])
            Traceback (most recent call last):
            ...
            ValueError: the diagram {{-3, 1}, {-2}, {-1, 2}, {3}} must be planar
        """
        d = list(d)
        if not d:
            return self.one()
        if d[0] in ZZ:
            return self._perm_to_Blst(d)
        d = to_set_partition(d, self._k)
        return self[self._base_diagrams(d)]

    def order(self):
        r"""
        Return the order of ``self``.

        The order of a partition algebra is defined as half of the number
        of nodes in the diagrams.

        EXAMPLES::

            sage: q = var('q')                                                          # needs sage.symbolic
            sage: PA = PartitionAlgebra(2, q)                                           # needs sage.symbolic
            sage: PA.order()                                                            # needs sage.symbolic
            2
        """
        return self._k

    def set_partitions(self):
        r"""
        Return the collection of underlying set partitions indexing the
        basis elements of a given diagram algebra.

        .. TODO:: Is this really necessary? deprecate?

        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: R.<x> = QQ[]
            sage: D = da.DiagramAlgebra(2, x, R, 'P', da.PartitionDiagrams(2))
            sage: list(D.set_partitions()) == list(da.PartitionDiagrams(2))
            True
        """
        return self.basis().keys()

    def _latex_term(self, diagram):
        r"""
        Return `\LaTeX` representation of ``diagram`` to draw
        diagram algebra element in latex using tikz.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: P = PartitionAlgebra(2, x, R)
            sage: latex(P([[1,2],[-2,-1]])) # indirect doctest
            \begin{tikzpicture}[scale = 0.5,thick, baseline={(0,-1ex/2)}]
            \tikzstyle{vertex} = [shape = circle, minimum size = 7pt, inner sep = 1pt]
            \node[vertex] (G--2) at (1.5, -1) [shape = circle, draw] {};
            \node[vertex] (G--1) at (0.0, -1) [shape = circle, draw] {};
            \node[vertex] (G-1) at (0.0, 1) [shape = circle, draw] {};
            \node[vertex] (G-2) at (1.5, 1) [shape = circle, draw] {};
            \draw[] (G--2) .. controls +(-0.5, 0.5) and +(0.5, 0.5) .. (G--1);
            \draw[] (G-1) .. controls +(0.5, -0.5) and +(-0.5, -0.5) .. (G-2);
            \end{tikzpicture}

            sage: latex(P.orbit_basis()([[1,2],[-2,-1]])) # indirect doctest
            \begin{tikzpicture}[scale = 0.5,thick, baseline={(0,-1ex/2)}]
            \tikzstyle{vertex} = [shape = circle, minimum size = 7pt, inner sep = 1pt]
            \node[vertex] (G--2) at (1.5, -1) [shape = circle, draw, fill] {};
            \node[vertex] (G--1) at (0.0, -1) [shape = circle, draw, fill] {};
            \node[vertex] (G-1) at (0.0, 1) [shape = circle, draw, fill] {};
            \node[vertex] (G-2) at (1.5, 1) [shape = circle, draw, fill] {};
            \draw[] (G--2) .. controls +(-0.5, 0.5) and +(0.5, 0.5) .. (G--1);
            \draw[] (G-1) .. controls +(0.5, -0.5) and +(-0.5, -0.5) .. (G-2);
            \end{tikzpicture}
        """
        return diagram_latex(diagram, fill=hasattr(self, '_fill'))

    # The following subclass provides a few additional methods for
    # (sub)partition algebra elements.
    class Element(CombinatorialFreeModule.Element):
        r"""
        An element of a diagram algebra.

        This subclass provides a few additional methods for
        partition algebra elements. Most element methods are
        already implemented elsewhere.
        """

        def diagram(self):
            r"""
            Return the underlying diagram of ``self`` if ``self`` is a basis
            element. Raises an error if ``self`` is not a basis element.

            EXAMPLES::

                sage: R.<x> = ZZ[]
                sage: P = PartitionAlgebra(2, x, R)
                sage: elt = 3*P([[1,2],[-2,-1]])
                sage: elt.diagram()
                {{-2, -1}, {1, 2}}
            """
            if len(self) != 1:
                raise ValueError("this is only defined for basis elements")
            PA = self.parent()
            ans = self.support_of_term()
            if ans not in PA.basis().keys():
                raise ValueError("element should be keyed by a diagram")
            return ans

        def diagrams(self):
            r"""
            Return the diagrams in the support of ``self``.

            EXAMPLES::

                sage: R.<x> = ZZ[]
                sage: P = PartitionAlgebra(2, x, R)
                sage: elt = 3*P([[1,2],[-2,-1]]) + P([[1,2],[-2], [-1]])
                sage: sorted(elt.diagrams(), key=str)
                [{{-2, -1}, {1, 2}}, {{-2}, {-1}, {1, 2}}]
            """
            return self.support()


class UnitDiagramMixin:
    """
    Mixin class for diagram algebras that have the unit indexed by
    the :func:`identity_set_partition`.
    """
    @cached_method
    def one_basis(self):
        r"""
        The following constructs the identity element of ``self``.

        It is not called directly; instead one should use ``DA.one()`` if
        ``DA`` is a defined diagram algebra.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: P = PartitionAlgebra(2, x, R)
            sage: P.one_basis()
            {{-2, 2}, {-1, 1}}
        """
        return self._base_diagrams(identity_set_partition(self._k))


class DiagramBasis(DiagramAlgebra):
    """
    Abstract base class for diagram algebras in the diagram basis.
    """

    def product_on_basis(self, d1, d2):
        r"""
        Return the product `D_{d_1} D_{d_2}` by two basis diagrams.

        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: R.<x> = QQ[]
            sage: D = da.DiagramBasis(2, x, R, 'P', da.PartitionDiagrams(2))
            sage: sp = da.PartitionDiagrams(2)([[1,2],[-1,-2]])
            sage: D.product_on_basis(sp, sp)
            x*P{{-2, -1}, {1, 2}}
        """
        if not self._indices.is_parent_of(d1):
            d1 = self._indices(d1)
        if not self._indices.is_parent_of(d2):
            d2 = self._indices(d2)
        (composite_diagram, loops_removed) = d1.compose(d2, check=False)
        return self.term(composite_diagram, self._q**loops_removed)


class PartitionAlgebra(DiagramBasis, UnitDiagramMixin):
    r"""
    A partition algebra.

    A partition algebra of rank `k` over a given ground ring `R` is an
    algebra with (`R`-module) basis indexed by the collection of set
    partitions of `\{1, \ldots, k, -1, \ldots, -k\}`. Each such set
    partition can be represented by a graph on nodes `\{1, \ldots, k, -1,
    \ldots, -k\}` arranged in two rows, with nodes `1, \ldots, k` in the
    top row from left to right and with nodes `-1, \ldots, -k` in the
    bottom row from left to right, and edges drawn such that the connected
    components of the graph are precisely the parts of the set partition.
    (This choice of edges is often not unique, and so there are often many
    graphs representing one and the same set partition; the representation
    nevertheless is useful and vivid. We often speak of "diagrams" to mean
    graphs up to such equivalence of choices of edges; of course, we could
    just as well speak of set partitions.)

    There is not just one partition algebra of given rank over a given
    ground ring, but rather a whole family of them, indexed by the
    elements of `R`. More precisely, for every `q \in R`, the partition
    algebra of rank `k` over `R` with parameter `q` is defined to be the
    `R`-algebra with basis the collection of all set partitions of
    `\{1, \ldots, k, -1, \ldots, -k\}`, where the product of two basis
    elements is given by the rule

    .. MATH::

        a \cdot b = q^N (a \circ b),

    where `a \circ b` is the composite set partition obtained by placing
    the diagram (i.e., graph) of `a` above the diagram of `b`, identifying
    the bottom row nodes of `a` with the top row nodes of `b`, and
    omitting any closed "loops" in the middle. The number `N` is the
    number of connected components formed by the omitted loops.

    The parameter `q` is a deformation parameter. Taking `q = 1` produces
    the semigroup algebra (over the base ring) of the partition monoid,
    in which the product of two set partitions is simply given by their
    composition.

    The partition algebra is regarded as an example of a "diagram algebra"
    due to the fact that its natural basis is given by certain graphs
    often called diagrams.

    There are a number of predefined elements for the partition algebra.
    We define the cup/cap pair by :meth:`a()`. The simple transpositions
    are denoted :meth:`s()`. Finally, we define elements :meth:`e()`,
    where if `i = (2r+1)/2`, then ``e(i)`` contains the blocks `\{r+1\}`
    and `\{-r-1\}` and if `i \in \ZZ`, then `e_i` contains the block
    `\{-i, -i-1, i, i+1\}`, with all other blocks being `\{-j, j\}`.
    So we have::

        sage: P = PartitionAlgebra(4, 0)
        sage: P.a(2)
        P{{-4, 4}, {-3, -2}, {-1, 1}, {2, 3}}
        sage: P.e(3/2)
        P{{-4, 4}, {-3, 3}, {-2}, {-1, 1}, {2}}
        sage: P.e(2)
        P{{-4, 4}, {-3, -2, 2, 3}, {-1, 1}}
        sage: P.e(5/2)
        P{{-4, 4}, {-3}, {-2, 2}, {-1, 1}, {3}}
        sage: P.s(2)
        P{{-4, 4}, {-3, 2}, {-2, 3}, {-1, 1}}

    An excellent reference for partition algebras and their various
    subalgebras (Brauer algebra, Temperley--Lieb algebra, etc) is the
    paper [HR2005]_.

    INPUT:

    - ``k`` -- rank of the algebra

    - ``q`` -- the deformation parameter `q`

    - ``base_ring`` -- (default: ``None``) a ring containing ``q``; if
      ``None``, then Sage automatically chooses the parent of ``q``

    - ``prefix`` -- (default: ``'P'``) a label for the basis elements

    EXAMPLES:

    The following shorthand simultaneously defines the univariate polynomial
    ring over the rationals as well as the variable ``x``::

        sage: R.<x> = PolynomialRing(QQ)
        sage: R
        Univariate Polynomial Ring in x over Rational Field
        sage: x
        x
        sage: x.parent() is R
        True

    We now define the partition algebra of rank `2` with parameter ``x``
    over `\ZZ` in the usual (diagram) basis::

        sage: R.<x> = ZZ[]
        sage: A2 = PartitionAlgebra(2, x, R)
        sage: A2
        Partition Algebra of rank 2 with parameter x
         over Univariate Polynomial Ring in x over Integer Ring
        sage: A2.basis().keys()
        Partition diagrams of order 2
        sage: A2.basis().keys()([[-2, 1, 2], [-1]])
        {{-2, 1, 2}, {-1}}
        sage: A2.basis().list()
        [P{{-2, -1, 1, 2}}, P{{-2, 1, 2}, {-1}},
         P{{-2}, {-1, 1, 2}}, P{{-2, -1}, {1, 2}},
         P{{-2}, {-1}, {1, 2}}, P{{-2, -1, 1}, {2}},
         P{{-2, 1}, {-1, 2}}, P{{-2, 1}, {-1}, {2}},
         P{{-2, 2}, {-1, 1}}, P{{-2, -1, 2}, {1}},
         P{{-2, 2}, {-1}, {1}}, P{{-2}, {-1, 1}, {2}},
         P{{-2}, {-1, 2}, {1}}, P{{-2, -1}, {1}, {2}},
         P{{-2}, {-1}, {1}, {2}}]
        sage: E = A2([[1,2],[-2,-1]]); E
        P{{-2, -1}, {1, 2}}
        sage: E in A2.basis().list()
        True
        sage: E^2
        x*P{{-2, -1}, {1, 2}}
        sage: E^5
        x^4*P{{-2, -1}, {1, 2}}
        sage: (A2([[2,-2],[-1,1]]) - 2*A2([[1,2],[-1,-2]]))^2
        (4*x-4)*P{{-2, -1}, {1, 2}} + P{{-2, 2}, {-1, 1}}

    Next, we construct an element::

        sage: a2 = A2.an_element(); a2
        3*P{{-2}, {-1, 1, 2}} + 2*P{{-2, -1, 1, 2}} + 2*P{{-2, 1, 2}, {-1}}

    There is a natural embedding into partition algebras on more
    elements, by adding identity strands::

        sage: A4 = PartitionAlgebra(4, x, R)
        sage: A4(a2)
        3*P{{-4, 4}, {-3, 3}, {-2}, {-1, 1, 2}}
         + 2*P{{-4, 4}, {-3, 3}, {-2, -1, 1, 2}}
         + 2*P{{-4, 4}, {-3, 3}, {-2, 1, 2}, {-1}}

    Thus, the empty partition corresponds to the identity::

        sage: A4([])
        P{{-4, 4}, {-3, 3}, {-2, 2}, {-1, 1}}
        sage: A4(5)
        5*P{{-4, 4}, {-3, 3}, {-2, 2}, {-1, 1}}

    The group algebra of the symmetric group is a subalgebra::

        sage: S3 = SymmetricGroupAlgebra(ZZ, 3)
        sage: s3 = S3.an_element(); s3
        [1, 2, 3] + 2*[1, 3, 2] + 3*[2, 1, 3] + [3, 1, 2]
        sage: A4(s3)
        P{{-4, 4}, {-3, 1}, {-2, 3}, {-1, 2}}
         + 2*P{{-4, 4}, {-3, 2}, {-2, 3}, {-1, 1}}
         + 3*P{{-4, 4}, {-3, 3}, {-2, 1}, {-1, 2}}
         + P{{-4, 4}, {-3, 3}, {-2, 2}, {-1, 1}}
        sage: A4([2,1])
        P{{-4, 4}, {-3, 3}, {-2, 1}, {-1, 2}}

    Be careful not to confuse the embedding of the group algebra of
    the symmetric group with the embedding of partial set partitions.
    The latter are embedded by adding the parts `\{i,-i\}` if
    possible, and singletons sets for the remaining parts::

        sage: A4([[2,1]])
        P{{-4, 4}, {-3, 3}, {-2}, {-1}, {1, 2}}
        sage: A4([[-1,3],[-2,-3,1]])
        P{{-4, 4}, {-3, -2, 1}, {-1, 3}, {2}}

    Another subalgebra is the Brauer algebra, which has perfect
    matchings as basis elements.  The group algebra of the
    symmetric group is in fact a subalgebra of the Brauer algebra::

        sage: B3 = BrauerAlgebra(3, x, R)
        sage: b3 = B3(s3); b3
        B{{-3, 1}, {-2, 3}, {-1, 2}} + 2*B{{-3, 2}, {-2, 3}, {-1, 1}}
         + 3*B{{-3, 3}, {-2, 1}, {-1, 2}} + B{{-3, 3}, {-2, 2}, {-1, 1}}

    An important basis of the partition algebra is the
    :meth:`orbit basis <orbit_basis>`::

        sage: O2 = A2.orbit_basis()
        sage: o2 = O2([[1,2],[-1,-2]]) + O2([[1,2,-1,-2]]); o2
        O{{-2, -1}, {1, 2}} + O{{-2, -1, 1, 2}}

    The diagram basis element corresponds to the sum of all orbit
    basis elements indexed by coarser set partitions::

        sage: A2(o2)
        P{{-2, -1}, {1, 2}}

    We can convert back from the orbit basis to the diagram basis::

        sage: o2 = O2.an_element(); o2
        3*O{{-2}, {-1, 1, 2}} + 2*O{{-2, -1, 1, 2}} + 2*O{{-2, 1, 2}, {-1}}
        sage: A2(o2)
        3*P{{-2}, {-1, 1, 2}} - 3*P{{-2, -1, 1, 2}} + 2*P{{-2, 1, 2}, {-1}}

    One can work with partition algebras using a symbol for the parameter,
    leaving the base ring unspecified. This implies that the underlying
    base ring is Sage's symbolic ring.

    ::

        sage: # needs sage.symbolic
        sage: q = var('q')
        sage: PA = PartitionAlgebra(2, q); PA
        Partition Algebra of rank 2 with parameter q over Symbolic Ring
        sage: PA([[1,2],[-2,-1]])^2 == q*PA([[1,2],[-2,-1]])
        True
        sage: ((PA([[2, -2], [1, -1]]) - 2*PA([[-2, -1], [1, 2]]))^2
        ....:   == (4*q-4)*PA([[1, 2], [-2, -1]]) + PA([[2, -2], [1, -1]]))
        True

    The identity element of the partition algebra is the set
    partition `\{\{1,-1\}, \{2,-2\}, \ldots, \{k,-k\}\}`::

        sage: # needs sage.symbolic
        sage: P = PA.basis().list()
        sage: PA.one()
        P{{-2, 2}, {-1, 1}}
        sage: PA.one() * P[7] == P[7]
        True
        sage: P[7] * PA.one() == P[7]
        True

    We now give some further examples of the use of the other arguments.
    One may wish to "specialize" the parameter to a chosen element of
    the base ring::

        sage: R.<q> = RR[]
        sage: PA = PartitionAlgebra(2, q, R, prefix='B')
        sage: PA
        Partition Algebra of rank 2 with parameter q over
         Univariate Polynomial Ring in q over Real Field with 53 bits of precision
        sage: PA([[1,2],[-1,-2]])
        1.00000000000000*B{{-2, -1}, {1, 2}}
        sage: PA = PartitionAlgebra(2, 5, base_ring=ZZ, prefix='B')
        sage: PA
        Partition Algebra of rank 2 with parameter 5 over Integer Ring
        sage: ((PA([[2, -2], [1, -1]]) - 2*PA([[-2, -1], [1, 2]]))^2
        ....:   == 16*PA([[-2, -1], [1, 2]]) + PA([[2, -2], [1, -1]]))
        True

    Symmetric group algebra elements and elements from other subalgebras
    of the partition algebra (e.g., ``BrauerAlgebra`` and
    ``TemperleyLiebAlgebra``) can also be coerced into the partition algebra::

        sage: # needs sage.symbolic
        sage: S = SymmetricGroupAlgebra(SR, 2)
        sage: B = BrauerAlgebra(2, x, SR)
        sage: A = PartitionAlgebra(2, x, SR)
        sage: S([2,1]) * A([[1,-1],[2,-2]])
        P{{-2, 1}, {-1, 2}}
        sage: B([[-1,-2],[2,1]]) * A([[1],[-1],[2,-2]])
        P{{-2}, {-1}, {1, 2}}
        sage: A([[1],[-1],[2,-2]]) * B([[-1,-2],[2,1]])
        P{{-2, -1}, {1}, {2}}

    The same is true if the elements come from a subalgebra of a partition
    algebra of smaller order, or if they are defined over a different
    base ring::

        sage: R = FractionField(ZZ['q']); q = R.gen()
        sage: S = SymmetricGroupAlgebra(ZZ, 2)
        sage: B = BrauerAlgebra(2, q, ZZ[q])
        sage: A = PartitionAlgebra(3, q, R)
        sage: S([2,1]) * A([[1,-1],[2,-3],[3,-2]])
        P{{-3, 1}, {-2, 3}, {-1, 2}}
        sage: A(B([[-1,-2],[2,1]]))
        P{{-3, 3}, {-2, -1}, {1, 2}}

    TESTS:

    A computation that returned an incorrect result until :issue:`15958`::

        sage: A = PartitionAlgebra(1,17)
        sage: g = SetPartitionsAk(1).list()
        sage: a = A[g[1]]
        sage: a
        P{{-1}, {1}}
        sage: a*a
        17*P{{-1}, {1}}

    Shorthands for working with basis elements are as follows::

        sage: # needs sage.symbolic
        sage: S = SymmetricGroupAlgebra(ZZ, 3)
        sage: A = PartitionAlgebra(3, x, SR)
        sage: A([[1,3],[-1],[-3]])  # pair up the omitted nodes as `{-i, i}`, if possible
        P{{-3}, {-2, 2}, {-1}, {1, 3}}
        sage: A([[1,3],[-1],[-3]]) == A[[1,3],[-1],[-3]]
        True
        sage: A([[1,2]])
        P{{-3, 3}, {-2}, {-1}, {1, 2}}
        sage: A([[1,2]]) == A[[1,2]]
        True
        sage: A([2,3,1])  # permutations in one-line notation are imported as well
        P{{-3, 2}, {-2, 1}, {-1, 3}}
        sage: A([2,3,1]) == A(S([2,3,1]))
        True
    """
    @staticmethod
    def __classcall_private__(cls, k, q, base_ring=None, prefix='P'):
        r"""
        Standardize the input by getting the base ring from the parent of
        the parameter ``q`` if no ``base_ring`` is given.

        TESTS::

            sage: R.<q> = QQ[]
            sage: PA1 = PartitionAlgebra(2, q)
            sage: PA2 = PartitionAlgebra(2, q, R, 'P')
            sage: PA1 is PA2
            True
        """
        if base_ring is None:
            base_ring = q.parent()
        return super().__classcall__(cls, k, q, base_ring, prefix)

    # The following is the basic constructor method for the class.
    # The purpose of the "prefix" is to label the basis elements
    def __init__(self, k, q, base_ring, prefix):
        r"""
        Initialize ``self``.

        TESTS::

            sage: R.<q> = QQ[]
            sage: PA = PartitionAlgebra(2, q, R)
            sage: TestSuite(PA).run()
        """
        self._k = k
        self._prefix = prefix
        self._q = base_ring(q)
        DiagramAlgebra.__init__(self, k, q, base_ring, prefix, PartitionDiagrams(k))

    def _element_constructor_(self, x):
        r"""
        Construct an element of ``self``.

        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: R.<x> = QQ[]
            sage: PA = PartitionAlgebra(2, x, R, 'P')
            sage: PA([]) == PA.one()
            True
            sage: D = da.DiagramAlgebra(2, x, R, 'P', da.PartitionDiagrams(2))
            sage: D([]) == D.one()
            Traceback (most recent call last):
            ...
            ValueError: invalid input of []
            sage: sp = da.to_set_partition( [[1,2], [-1,-2]] )
            sage: b_elt = D(sp); b_elt
            P{{-2, -1}, {1, 2}}
            sage: b_elt in D
            True
            sage: D([[1,2],[-1,-2]]) == b_elt
            True
            sage: D([{1,2},{-1,-2}]) == b_elt
            True
            sage: S = SymmetricGroupAlgebra(R,2)
            sage: D(S([2,1]))
            P{{-2, 1}, {-1, 2}}
            sage: D2 = da.DiagramAlgebra(2, x, R, 'P', da.PlanarDiagrams(2))
            sage: D2(S([1,2]))
            P{{-2, 2}, {-1, 1}}
            sage: D2(S([2,1]))
            Traceback (most recent call last):
            ...
            ValueError: the diagram {{-2, 1}, {-1, 2}} must be planar
        """
        # coercion from basis keys
        if self.basis().keys().is_parent_of(x):
            return self.basis()[x]

        # conversion from (smaller) diagram or permutation
        if isinstance(x, (AbstractPartitionDiagram, list, tuple, Permutations.Element)):
            return self._diag_to_Blst(x)

        # conversion from orbit basis
        if (isinstance(x, OrbitBasis.Element)
                and self.base_ring().has_coerce_map_from(x.parent().base_ring())):
            return self(x.parent().to_diagram_basis(x))

        # conversion from SubPartitionAlgebra
        if (isinstance(x, (PartitionAlgebra.Element, SubPartitionAlgebra.Element))
                and self.has_coerce_map_from(x.parent().base_ring())):
            return sum(a * self._diag_to_Blst(d) for (d,a) in x)

        return super()._element_constructor_(x)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: R.<q> = QQ[]
            sage: PartitionAlgebra(2, q, R)
            Partition Algebra of rank 2 with parameter q
             over Univariate Polynomial Ring in q over Rational Field
        """
        return "Partition Algebra of rank {} with parameter {} over {}".format(
                self._k, self._q, self.base_ring())

    def _coerce_map_from_(self, R):
        """
        Return a coerce map from ``R`` if one exists and ``None`` otherwise.

        .. TODO::

            - Refactor some of these generic morphisms as compositions
              of morphisms.
            - Allow for coercion if base_rings and parameters are the same,
              up to relabeling/isomorphism?

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: S = SymmetricGroupAlgebra(R, 4)
            sage: A = PartitionAlgebra(4, x, R)
            sage: O = A.orbit_basis()
            sage: A._coerce_map_from_(S)
            Generic morphism:
              From: Symmetric group algebra of order 4 over Univariate Polynomial Ring in x over Rational Field
              To:   Partition Algebra of rank 4 with parameter x over Univariate Polynomial Ring in x over Rational Field
            sage: A._coerce_map_from_(O)
            Generic morphism:
              From: Orbit basis of Partition Algebra of rank 4 with parameter x over Univariate Polynomial Ring in x over Rational Field
              To:   Partition Algebra of rank 4 with parameter x over Univariate Polynomial Ring in x over Rational Field
            sage: Sp3 = SymmetricGroupAlgebra(ZZ, 3)
            sage: A._coerce_map_from_(Sp3)
            Generic morphism:
              From: Symmetric group algebra of order 3 over Integer Ring
              To:   Partition Algebra of rank 4 with parameter x over Univariate Polynomial Ring in x over Rational Field
            sage: B3 = BrauerAlgebra(3, x, R)
            sage: A._coerce_map_from_(B3)
            Generic morphism:
              From: Brauer Algebra of rank 3 with parameter x over Univariate Polynomial Ring in x over Rational Field
              To:   Partition Algebra of rank 4 with parameter x over Univariate Polynomial Ring in x over Rational Field
            sage: A3 = PartitionAlgebra(3, x, R)
            sage: A._coerce_map_from_(A3)
            Generic morphism:
              From: Partition Algebra of rank 3 with parameter x over Univariate Polynomial Ring in x over Rational Field
              To:   Partition Algebra of rank 4 with parameter x over Univariate Polynomial Ring in x over Rational Field
            sage: O3 = A3.orbit_basis()
            sage: A._coerce_map_from_(O3)
            Generic morphism:
              From: Orbit basis of Partition Algebra of rank 3 with parameter x over Univariate Polynomial Ring in x over Rational Field
              To:   Partition Algebra of rank 4 with parameter x over Univariate Polynomial Ring in x over Rational Field

        TESTS::

            sage: elt = O3.an_element(); elt
            2*O{{-3, -2, -1, 1, 2, 3}} + 2*O{{-3, -2, 1, 2, 3}, {-1}}
             + 3*O{{-3, -1, 1, 2, 3}, {-2}}
            sage: A._coerce_map_from_(O3)(elt)
            -3*P{{-4, 4}, {-3, -2, -1, 1, 2, 3}}
             + 2*P{{-4, 4}, {-3, -2, 1, 2, 3}, {-1}}
             + 3*P{{-4, 4}, {-3, -1, 1, 2, 3}, {-2}}
        """
        # coerce from Orbit basis.
        if isinstance(R, OrbitBasis):
            if R._k <= self._k and self.base_ring().has_coerce_map_from(R.base_ring()):
                return R.module_morphism(self._orbit_to_diagram_on_basis, codomain=self)
            return None

        # coerce from sub-partition algebras.
        if isinstance(R, (PartitionAlgebra, SubPartitionAlgebra)):
            if R._k <= self._k and self.base_ring().has_coerce_map_from(R.base_ring()):
                return R.module_morphism(self._diag_to_Blst, codomain=self)
            return None

        # coerce from Symmetric group algebras.
        if isinstance(R, SymmetricGroupAlgebra_n):
            if R.n <= self._k and self.base_ring().has_coerce_map_from(R.base_ring()):
                return R.module_morphism(self._perm_to_Blst, codomain=self)
            return None
        return super()._coerce_map_from_(R)

    def orbit_basis(self):
        r"""
        Return the orbit basis of ``self``.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: P2 = PartitionAlgebra(2, x, R)
            sage: O2 = P2.orbit_basis(); O2
            Orbit basis of Partition Algebra of rank 2 with parameter x over
             Univariate Polynomial Ring in x over Rational Field
            sage: pp = 7 * P2[{-1}, {-2, 1, 2}] - 2 * P2[{-2}, {-1, 1}, {2}]; pp
            -2*P{{-2}, {-1, 1}, {2}} + 7*P{{-2, 1, 2}, {-1}}
            sage: op = pp.to_orbit_basis(); op
            -2*O{{-2}, {-1, 1}, {2}} - 2*O{{-2}, {-1, 1, 2}}
             - 2*O{{-2, -1, 1}, {2}} + 5*O{{-2, -1, 1, 2}}
             + 7*O{{-2, 1, 2}, {-1}} - 2*O{{-2, 2}, {-1, 1}}
            sage: op == O2(op)
            True
            sage: pp * op.leading_term()
            4*P{{-2}, {-1, 1}, {2}} - 4*P{{-2, -1, 1}, {2}}
             + 14*P{{-2, -1, 1, 2}} - 14*P{{-2, 1, 2}, {-1}}
        """
        return OrbitBasis(self)

    def _orbit_to_diagram_on_basis(self, d):
        r"""
        Return the orbit basis element, indexed by the partition
        diagram ``d``, in the diagram basis of the partition algebra.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: P2 = PartitionAlgebra(2, x, R)
            sage: from sage.combinat.diagram_algebras import PartitionDiagrams
            sage: PD = PartitionDiagrams(2)
            sage: P2._orbit_to_diagram_on_basis(PD([[1,2,-2],[-1]]))
            -P{{-2, -1, 1, 2}} + P{{-2, 1, 2}, {-1}}
        """
        # Moebius inversion in the poset of coarsenings of ``d``
        SPd = SetPartitions(len(d))
        return self.sum((-1)**(len(d)-len(sp)) * prod(ZZ(len(p)-1).factorial() for p in sp)
                        * self([sum((list(d[i-1]) for i in p),[]) for p in sp])
                        for sp in SPd)

    @cached_method
    def a(self, i):
        r"""
        Return the element `a_i` in ``self``.

        The element `a_i` is the cap and cup at `(i, i+1)`, so it contains
        the blocks `\{i, i+1\}`, `\{-i, -i-1\}`.  Other blocks are of the
        form `\{-j, j\}`.

        INPUT:

        - ``i`` -- integer between 1 and `k-1`

        EXAMPLES::

            sage: R.<n> = QQ[]
            sage: P3 = PartitionAlgebra(3, n)
            sage: P3.a(1)
            P{{-3, 3}, {-2, -1}, {1, 2}}
            sage: P3.a(2)
            P{{-3, -2}, {-1, 1}, {2, 3}}

            sage: P3 = PartitionAlgebra(5/2, n)
            sage: P3.a(1)
            P{{-3, 3}, {-2, -1}, {1, 2}}
            sage: P3.a(2)
            Traceback (most recent call last):
            ...
            ValueError: i must be an integer between 1 and 1
        """
        if i <= 0 or i >= floor(self._k):
            raise ValueError("i must be an integer between 1 and {}".format(floor(self._k)-1))
        B = self.basis()
        SP = B.keys()
        D = [[-j, j] for j in range(1, ceil(self._k)+1)]
        D[i-1] = [i,i+1]
        D[i] = [-i,-(i+1)]
        return B[SP(D)]

    generator_a = a

    @cached_method
    def e(self, i):
        r"""
        Return the element `e_i` in ``self``.

        If `i = (2r+1)/2`, then `e_i` contains the blocks `\{r+1\}` and
        `\{-r-1\}`.  If `i \in \ZZ`, then `e_i` contains the block
        `\{-i, -i-1, i, i+1\}`.  Other blocks are of the form `\{-j, j\}`.

        INPUT:

        - ``i`` -- half integer between `1/2` and `k-1/2`

        EXAMPLES::

            sage: R.<n> = QQ[]
            sage: P3 = PartitionAlgebra(3, n)
            sage: P3.e(1)
            P{{-3, 3}, {-2, -1, 1, 2}}
            sage: P3.e(2)
            P{{-3, -2, 2, 3}, {-1, 1}}
            sage: P3.e(1/2)
            P{{-3, 3}, {-2, 2}, {-1}, {1}}
            sage: P3.e(5/2)
            P{{-3}, {-2, 2}, {-1, 1}, {3}}
            sage: P3.e(0)
            Traceback (most recent call last):
            ...
            ValueError: i must be an (half) integer between 1/2 and 5/2
            sage: P3.e(3)
            Traceback (most recent call last):
            ...
            ValueError: i must be an (half) integer between 1/2 and 5/2

            sage: P2h = PartitionAlgebra(5/2,n)
            sage: [P2h.e(k/2) for k in range(1,5)]
            [P{{-3, 3}, {-2, 2}, {-1}, {1}},
             P{{-3, 3}, {-2, -1, 1, 2}},
             P{{-3, 3}, {-2}, {-1, 1}, {2}},
             P{{-3, -2, 2, 3}, {-1, 1}}]
        """
        if i <= 0 or i >= self._k:
            raise ValueError("i must be an (half) integer between 1/2 and {}".format((2*self._k-1)/2))
        B = self.basis()
        SP = B.keys()
        if i in ZZ:
            i -= 1
            D = [[-j, j] for j in range(1, ceil(self._k)+1)]
            D[i] += D.pop(i+1)
            return B[SP(D)]
        else:
            i = ceil(i)
            D = [[-j, j] for j in range(1, ceil(self._k)+1)]
            D[i-1] = [-i]
            D.append([i])
            return B[SP(D)]

    generator_e = e

    @cached_method
    def s(self, i):
        r"""
        Return the ``i``-th simple transposition `s_i` in ``self``.

        Borrowing the notation from the symmetric group, the `i`-th
        simple transposition `s_i` has blocks of the form `\{-i, i+1\}`,
        `\{-i-1, i\}`.  Other blocks are of the form `\{-j, j\}`.

        INPUT:

        - ``i`` -- integer between 1 and `k-1`

        EXAMPLES::

            sage: R.<n> = QQ[]
            sage: P3 = PartitionAlgebra(3, n)
            sage: P3.s(1)
            P{{-3, 3}, {-2, 1}, {-1, 2}}
            sage: P3.s(2)
            P{{-3, 2}, {-2, 3}, {-1, 1}}

            sage: R.<n> = ZZ[]
            sage: P2h = PartitionAlgebra(5/2,n)
            sage: P2h.s(1)
            P{{-3, 3}, {-2, 1}, {-1, 2}}
        """
        if i not in ZZ or i <= 0 or i >= self._k:
            raise ValueError("i must be an integer between 1 and {}".format(self._k-1))
        B = self.basis()
        SP = B.keys()
        D = [[-j, j] for j in range(1, ceil(self._k)+1)]
        D[i-1] = [-(i+1), i]
        D[i] = [-i, i+1]
        return B[SP(D)]

    generator_s = s

    @cached_method
    def sigma(self, i):
        r"""
        Return the element `\sigma_i` from [Eny2012]_ of ``self``.

        INPUT:

        - ``i`` -- half integer between `1/2` and `k-1/2`

        .. NOTE::

            In [Cre2020]_ and [Eny2013]_, these are the elements `\sigma_{2i}`.

        EXAMPLES::

            sage: R.<n> = QQ[]
            sage: P3 = PartitionAlgebra(3, n)
            sage: P3.sigma(1)
            P{{-3, 3}, {-2, 2}, {-1, 1}}
            sage: P3.sigma(3/2)
            P{{-3, 3}, {-2, 1}, {-1, 2}}
            sage: P3.sigma(2)
            -P{{-3, -1, 1, 3}, {-2, 2}} + P{{-3, -1, 3}, {-2, 1, 2}}
             + P{{-3, 1, 3}, {-2, -1, 2}} - P{{-3, 3}, {-2, -1, 1, 2}}
             + P{{-3, 3}, {-2, 2}, {-1, 1}}
            sage: P3.sigma(5/2)
            -P{{-3, -1, 1, 2}, {-2, 3}} + P{{-3, -1, 2}, {-2, 1, 3}}
             + P{{-3, 1, 2}, {-2, -1, 3}} - P{{-3, 2}, {-2, -1, 1, 3}}
             + P{{-3, 2}, {-2, 3}, {-1, 1}}

        We test the relations in Lemma 2.2.3(1) in [Cre2020]_ (v1)::

            sage: k = 4
            sage: R.<x> = QQ[]
            sage: P = PartitionAlgebra(k, x)
            sage: all(P.sigma(i/2).dual() == P.sigma(i/2)
            ....:     for i in range(1,2*k))
            True
            sage: all(P.sigma(i)*P.sigma(i+1/2) == P.sigma(i+1/2)*P.sigma(i) == P.s(i)
            ....:     for i in range(1,floor(k)))
            True
            sage: all(P.sigma(i)*P.e(i) == P.e(i)*P.sigma(i) == P.e(i)
            ....:     for i in range(1,floor(k)))
            True
            sage: all(P.sigma(i+1/2)*P.e(i) == P.e(i)*P.sigma(i+1/2) == P.e(i)
            ....:     for i in range(1,floor(k)))
            True

            sage: k = 9/2
            sage: R.<x> = QQ[]
            sage: P = PartitionAlgebra(k, x)
            sage: all(P.sigma(i/2).dual() == P.sigma(i/2)
            ....:     for i in range(1,2*k-1))
            True
            sage: all(P.sigma(i)*P.sigma(i+1/2) == P.sigma(i+1/2)*P.sigma(i) == P.s(i)
            ....:     for i in range(1,k-1/2))
            True
            sage: all(P.sigma(i)*P.e(i) == P.e(i)*P.sigma(i) == P.e(i)
            ....:     for i in range(1,floor(k)))
            True
            sage: all(P.sigma(i+1/2)*P.e(i) == P.e(i)*P.sigma(i+1/2) == P.e(i)
            ....:     for i in range(1,floor(k)))
            True
        """
        if i <= 0 or i >= self._k:
            raise ValueError("i must be an (half) integer between 1 and {}".format((2*self._k-1)/2))

        half = QQ.one() / 2
        if i in ZZ:
            if i == 1:
                return self.one()
            si = self.s(i)
            sim = self.s(i-1)
            x = self.e(i-1) * self.jucys_murphy_element(i-1) * si * self.e(i-1)
            return (sim * si * self.sigma(i-1) * si * sim
                    + x * si + si * x
                    - self.e(i-1) * self.jucys_murphy_element(i-1) * sim
                      * self.e(i) * self.e(i-half) * self.e(i-1)
                    - si * self.e(i-1) * self.e(i-half) * self.e(i) * sim
                      * self.jucys_murphy_element(i-1) * self.e(i-1) * si)
        else:
            j = ceil(i) - 1
            if j == 0:
                return self.zero()
            if j == 1:
                return self.s(1)
            si = self.s(j)
            sim = self.s(j-1)
            x = self.e(j-1) * self.jucys_murphy_element(j-1) * si * self.e(j-1)
            return (sim * si * self.sigma(i-1) * si * sim
                    + si * x * si + x
                    - si * self.e(j-1) * self.jucys_murphy_element(j-1) * sim
                      * self.e(j) * self.e(i-1) * self.e(j-1)
                    - self.e(j-1) * self.e(i-1) * self.e(j) * sim
                      * self.jucys_murphy_element(j-1) * self.e(j-1) * si)

    @cached_method
    def jucys_murphy_element(self, i):
        r"""
        Return the ``i``-th Jucys-Murphy element `L_i` from [Eny2012]_.

        INPUT:

        - ``i`` -- half integer between `1/2` and `k`

        ALGORITHM:

        We use the recursive definition for `L_{2i}` given in [Cre2020]_.
        See also [Eny2012]_ and [Eny2013]_.

        .. NOTE::

            `L_{1/2}` and `L_1` differs from [HR2005]_.

        EXAMPLES::

            sage: R.<n> = QQ[]
            sage: P3 = PartitionAlgebra(3, n)
            sage: P3.jucys_murphy_element(1/2)
            0
            sage: P3.jucys_murphy_element(1)
            P{{-3, 3}, {-2, 2}, {-1}, {1}}
            sage: P3.jucys_murphy_element(2)
            P{{-3, 3}, {-2}, {-1, 1}, {2}} - P{{-3, 3}, {-2}, {-1, 1, 2}}
             + P{{-3, 3}, {-2, -1}, {1, 2}} - P{{-3, 3}, {-2, -1, 1}, {2}}
             + P{{-3, 3}, {-2, 1}, {-1, 2}}
            sage: P3.jucys_murphy_element(3/2)
            n*P{{-3, 3}, {-2, -1, 1, 2}} - P{{-3, 3}, {-2, -1, 2}, {1}}
             - P{{-3, 3}, {-2, 1, 2}, {-1}} + P{{-3, 3}, {-2, 2}, {-1, 1}}
            sage: P3.L(3/2) * P3.L(2) == P3.L(2) * P3.L(3/2)
            True

        We test the relations in Lemma 2.2.3(2) in [Cre2020]_ (v1)::

            sage: k = 4
            sage: R.<n> = QQ[]
            sage: P = PartitionAlgebra(k, n)
            sage: L = [P.L(i/2) for i in range(1,2*k+1)]
            sage: all(x.dual() == x for x in L)
            True
            sage: all(x * y == y * x for x, y in Subsets(L, 2))  # long time
            True
            sage: Lsum = sum(L)
            sage: gens = [P.s(i) for i in range(1,k)]
            sage: gens += [P.e(i/2) for i in range(1,2*k)]
            sage: all(x * Lsum == Lsum * x for x in gens)
            True

        Also the relations in Lemma 2.2.3(3) in [Cre2020]_ (v1)::

            sage: all(P.e((2*i+1)/2) * P.sigma(2*i/2) * P.e((2*i+1)/2)
            ....:     == (n - P.L((2*i-1)/2)) * P.e((2*i+1)/2) for i in range(1,k))
            True
            sage: all(P.e(i/2) * (P.L(i/2) + P.L((i+1)/2))
            ....:     == (P.L(i/2) + P.L((i+1)/2)) * P.e(i/2)
            ....:     == n * P.e(i/2) for i in range(1,2*k))
            True
            sage: all(P.sigma(2*i/2) * P.e((2*i-1)/2) * P.e(2*i/2)
            ....:     == P.L(2*i/2) * P.e(2*i/2) for i in range(1,k))
            True
            sage: all(P.e(2*i/2) * P.e((2*i-1)/2) * P.sigma(2*i/2)
            ....:     == P.e(2*i/2) * P.L(2*i/2) for i in range(1,k))
            True
            sage: all(P.sigma((2*i+1)/2) * P.e((2*i+1)/2) * P.e(2*i/2)
            ....:     == P.L(2*i/2) * P.e(2*i/2) for i in range(1,k))
            True
            sage: all(P.e(2*i/2) * P.e((2*i+1)/2) * P.sigma((2*i+1)/2)
            ....:     == P.e(2*i/2) * P.L(2*i/2) for i in range(1,k))
            True

        The same tests for a half integer partition algebra::

            sage: k = 7/2
            sage: R.<n> = QQ[]
            sage: P = PartitionAlgebra(k, n)
            sage: L = [P.L(i/2) for i in range(1,2*k+1)]
            sage: all(x.dual() == x for x in L)
            True
            sage: all(x * y == y * x for x, y in Subsets(L, 2))  # long time
            True
            sage: Lsum = sum(L)
            sage: gens = [P.s(i) for i in range(1,k-1/2)]
            sage: gens += [P.e(i/2) for i in range(1,2*k)]
            sage: all(x * Lsum == Lsum * x for x in gens)
            True
            sage: all(P.e((2*i+1)/2) * P.sigma(2*i/2) * P.e((2*i+1)/2)
            ....:     == (n - P.L((2*i-1)/2)) * P.e((2*i+1)/2) for i in range(1,floor(k)))
            True
            sage: all(P.e(i/2) * (P.L(i/2) + P.L((i+1)/2))
            ....:     == (P.L(i/2) + P.L((i+1)/2)) * P.e(i/2)
            ....:     == n * P.e(i/2) for i in range(1,2*k))
            True
            sage: all(P.sigma(2*i/2) * P.e((2*i-1)/2) * P.e(2*i/2)
            ....:     == P.L(2*i/2) * P.e(2*i/2) for i in range(1,ceil(k)))
            True
            sage: all(P.e(2*i/2) * P.e((2*i-1)/2) * P.sigma(2*i/2)
            ....:     == P.e(2*i/2) * P.L(2*i/2) for i in range(1,ceil(k)))
            True
            sage: all(P.sigma((2*i+1)/2) * P.e((2*i+1)/2) * P.e(2*i/2)
            ....:     == P.L(2*i/2) * P.e(2*i/2) for i in range(1,floor(k)))
            True
            sage: all(P.e(2*i/2) * P.e((2*i+1)/2) * P.sigma((2*i+1)/2)
            ....:     == P.e(2*i/2) * P.L(2*i/2) for i in range(1,floor(k)))
            True
        """
        if i <= 0 or i > self._k:
            raise ValueError("i must be an (half) integer between 1/2 and {}".format(self._k))

        half = QQ.one() / 2
        if i in ZZ:
            if i == 1:
                return self.e(half)
            i -= 1
            L = self.jucys_murphy_element
            return ((self.s(i) * L(i)) * (self.s(i) - self.e(i))
                    - (self.e(i) * L(i)) * (self.s(i) - self.e(i+half)*self.e(i))
                    + self.sigma(i+half))
        else:
            j = ceil(i) - 1
            if j == 0:
                return self.zero()
            L = self.jucys_murphy_element
            return (self.s(j) * L(i-1) * self.s(j)
                    - self.e(j)*L(j)
                    + (self._q*self.one() - L(i-1) - L(j))*self.e(j)
                    + self.sigma(j))

    L = jucys_murphy_element

    def potts_representation(self, y=None):
        r"""
        Return the :class:`PottsRepresentation` with magnetic field
        direction ``y`` of ``self``.

        .. NOTE::

            The deformation parameter `d` of ``self`` must be a
            positive integer.

        INPUT:

        - ``y`` -- (optional) an integer between 1 and `d`; ignored
          if the order of ``self`` is an integer, otherwise the
          default is `1`

        EXAMPLES::

            sage: PA = algebras.Partition(5/2, QQ(4))
            sage: PR = PA.potts_representation()

            sage: PA = algebras.Partition(5/2, 3/2)
            sage: PA.potts_representation()
            Traceback (most recent call last):
            ...
            ValueError: the partition algebra deformation parameter must
             be a positive integer
        """
        if self.order() not in ZZ and y is None:
            y = ZZ.one()
        return PottsRepresentation(self, y)

    class Element(DiagramBasis.Element):
        def to_orbit_basis(self):
            """
            Return ``self`` in the orbit basis of the associated
            partition algebra.

            EXAMPLES::

                sage: R.<x> = QQ[]
                sage: P = PartitionAlgebra(2, x, R)
                sage: pp = P.an_element();
                sage: pp.to_orbit_basis()
                3*O{{-2}, {-1, 1, 2}} + 7*O{{-2, -1, 1, 2}} + 2*O{{-2, 1, 2}, {-1}}
                sage: pp = (3*P([[-2], [-1, 1, 2]]) + 2*P([[-2, -1, 1, 2]])
                ....:       + 2*P([[-2, 1, 2], [-1]])); pp
                3*P{{-2}, {-1, 1, 2}} + 2*P{{-2, -1, 1, 2}} + 2*P{{-2, 1, 2}, {-1}}
                sage: pp.to_orbit_basis()
                3*O{{-2}, {-1, 1, 2}} + 7*O{{-2, -1, 1, 2}} + 2*O{{-2, 1, 2}, {-1}}
            """
            OP = self.parent().orbit_basis()
            return OP(self)

        def dual(self):
            r"""
            Return the dual of ``self``.

            The dual of an element in the partition algebra is formed
            by taking the dual of each diagram in the support.

            EXAMPLES::

                sage: R.<x> = QQ[]
                sage: P = PartitionAlgebra(2, x, R)
                sage: elt = P.an_element(); elt
                3*P{{-2}, {-1, 1, 2}} + 2*P{{-2, -1, 1, 2}} + 2*P{{-2, 1, 2}, {-1}}
                sage: elt.dual()
                3*P{{-2, -1, 1}, {2}} + 2*P{{-2, -1, 1, 2}} + 2*P{{-2, -1, 2}, {1}}
            """
            P = self.parent()
            return P._from_dict({D.dual(): c for D, c in self._monomial_coefficients.items()},
                                remove_zeros=False)


class OrbitBasis(DiagramAlgebra):
    r"""
    The orbit basis of the partition algebra.

    Let `D_\pi` represent the diagram basis element indexed by the
    partition `\pi`, then (see equations (2.14), (2.17) and (2.18) of [BH2017]_)

    .. MATH::

        D_\pi = \sum_{\tau \geq \pi} O_\tau,

    where the sum is over all partitions `\tau` which are coarser than `\pi`
    and `O_\tau` is the orbit basis element indexed by the partition `\tau`.

    If `\mu_{2k}(\pi,\tau)` represents the Moebius function of the partition
    lattice, then

    .. MATH::

        O_\pi = \sum_{\tau \geq \pi} \mu_{2k}(\pi, \tau) D_\tau.

    If `\tau` is a partition of `\ell` blocks and the `i`-th block of
    `\tau` is a union of `b_i` blocks of `\pi`, then

    .. MATH::

        \mu_{2k}(\pi, \tau) = \prod_{i=1}^\ell (-1)^{b_i-1} (b_i-1)! .

    EXAMPLES::

        sage: R.<x> = QQ[]
        sage: P2 = PartitionAlgebra(2, x, R)
        sage: O2 = P2.orbit_basis(); O2
        Orbit basis of Partition Algebra of rank 2 with parameter x over
         Univariate Polynomial Ring in x over Rational Field
        sage: oa = O2([[1],[-1],[2,-2]]); ob = O2([[-1,-2,2],[1]]); oa, ob
        (O{{-2, 2}, {-1}, {1}}, O{{-2, -1, 2}, {1}})
        sage: oa * ob
        (x-2)*O{{-2, -1, 2}, {1}}

    We can convert between the two bases::

        sage: pa = P2(oa); pa
        2*P{{-2, -1, 1, 2}} - P{{-2, -1, 2}, {1}} - P{{-2, 1, 2}, {-1}}
         + P{{-2, 2}, {-1}, {1}} - P{{-2, 2}, {-1, 1}}
        sage: pa * ob
        (-x+2)*P{{-2, -1, 1, 2}} + (x-2)*P{{-2, -1, 2}, {1}}
        sage: _ == pa * P2(ob)
        True
        sage: O2(pa * ob)
        (x-2)*O{{-2, -1, 2}, {1}}

    Note that the unit in the orbit basis is not a single diagram,
    in contrast to the natural diagram basis::

        sage: P2.one()
        P{{-2, 2}, {-1, 1}}
        sage: O2.one()
        O{{-2, -1, 1, 2}} + O{{-2, 2}, {-1, 1}}
        sage: O2.one() == P2.one()
        True

    TESTS:

    Check that going between the two bases is the identity::

        sage: R.<x> = QQ[]
        sage: P2 = PartitionAlgebra(2, x, R)
        sage: O2 = P2.orbit_basis(); O2
        Orbit basis of Partition Algebra of rank 2 with parameter x over
         Univariate Polynomial Ring in x over Rational Field
        sage: PD = P2.basis().keys()
        sage: all(O2(P2(O2(m))) == O2(m) for m in PD)
        True
        sage: all(P2(O2(P2(m))) == P2(m) for m in PD)
        True
    """
    @staticmethod
    def __classcall_private__(cls, *args):
        """
        Normalize input to ensure a unique representation.

        INPUT:

        Either:

        - ``A`` -- an abstract diagram algebra

        or the arguments to construct a diagram algebra:

        - ``k`` -- the rank
        - ``q`` -- the parameter
        - ``R`` -- the base ring

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: P2 = PartitionAlgebra(2, x, R)
            sage: from sage.combinat.diagram_algebras import OrbitBasis
            sage: O2a = P2.orbit_basis()
            sage: O2b = OrbitBasis(P2)
            sage: O2c = OrbitBasis(2, x, R)
            sage: O2a is O2b and O2a is O2c
            True
            sage: O2d = OrbitBasis(2, x, QQ[x])
            sage: O2a is O2d
            True
        """
        if len(args) == 1:
            PA = args[0]
            if not isinstance(PA, DiagramAlgebra):
                raise ValueError("{} is not a partition algebra".format(PA))
            alg = PA
        elif len(args) != 3:
            raise ValueError("expected 1 or 3 arguments, received %s: %s" % (len(args), args))
        else:
            k, q, R = args
            q = R(q)
            alg = PartitionAlgebra(k, q, R)
        return super().__classcall__(cls, alg)

    def __init__(self, alg):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: O2 = PartitionAlgebra(2, -1, QQ).orbit_basis()
            sage: TestSuite(O2).run()
        """
        base_ring = alg.base_ring()
        k = alg._k
        q = alg._q
        diagrams = alg._base_diagrams
        # TODO: should we add additional categories?
        category = alg.category()
        DiagramAlgebra.__init__(self, k, q, base_ring, "O", diagrams, category)
        self._fill = True
        self._alg = alg

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: PartitionAlgebra(2, -1, QQ).orbit_basis()
            Orbit basis of Partition Algebra of rank 2 with parameter -1 over Rational Field
        """
        return "Orbit basis of {}".format(self._alg)

    def _element_constructor_(self, x):
        """
        Convert ``x`` into ``self``.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: P2 = PartitionAlgebra(2, x, R)
            sage: O2 = P2.orbit_basis()
            sage: O2(P2([]))
            O{{-2, -1, 1, 2}} + O{{-2, 2}, {-1, 1}}
            sage: O2(3).to_diagram_basis() == 3 * P2.one()
            True
            sage: O2(P2([[1,2,-2],[-1]]))
            O{{-2, -1, 1, 2}} + O{{-2, 1, 2}, {-1}}
        """
        if isinstance(x, (PartitionAlgebra.Element, SubPartitionAlgebra.Element)):
            return self._alg(x).to_orbit_basis()
        d = self._alg._diag_to_Blst(x).diagram()
        return CombinatorialFreeModule._element_constructor_(self, d)

    def _coerce_map_from_(self, R):
        r"""
        Return a coerce map from ``R`` if one exists and ``None`` otherwise.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: P2 = PartitionAlgebra(2, x, R)
            sage: O2 = P2.orbit_basis()
            sage: O2(P2([]))
            O{{-2, -1, 1, 2}} + O{{-2, 2}, {-1, 1}}
            sage: O2(3)
            3*O{{-2, -1, 1, 2}} + 3*O{{-2, 2}, {-1, 1}}
            sage: O2([[1,2,-2],[-1]])
            O{{-2, 1, 2}, {-1}}
        """
        if R is self._alg:
            return self._alg.module_morphism(self._diagram_to_orbit_on_basis, codomain=self)
        if self._alg.coerce_map_from(R):
            return self._coerce_map_via([self._alg], R)
        return super()._coerce_map_from_(R)

    @cached_method
    def one(self):
        """
        Return the element `1` of the partition algebra in the orbit basis.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: P2 = PartitionAlgebra(2, x, R)
            sage: O2 = P2.orbit_basis()
            sage: O2.one()
            O{{-2, -1, 1, 2}} + O{{-2, 2}, {-1, 1}}
        """
        PDs = self._base_diagrams
        base = SetPartitions()(identity_set_partition(self._k))
        brone = self.base_ring().one()
        return self._from_dict({PDs(d): brone for d in base.coarsenings()},
                               coerce=False, remove_zeros=False)

    def diagram_basis(self):
        """
        Return the associated partition algebra of ``self``
        in the diagram basis.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: O2 = PartitionAlgebra(2, x, R).orbit_basis()
            sage: P2 = O2.diagram_basis(); P2
            Partition Algebra of rank 2 with parameter x over Univariate
            Polynomial Ring in x over Rational Field
            sage: o2 = O2.an_element(); o2
            3*O{{-2}, {-1, 1, 2}} + 2*O{{-2, -1, 1, 2}} + 2*O{{-2, 1, 2}, {-1}}
            sage: P2(o2)
            3*P{{-2}, {-1, 1, 2}} - 3*P{{-2, -1, 1, 2}} + 2*P{{-2, 1, 2}, {-1}}

        TESTS::

            sage: R.<x> = QQ[]
            sage: P2 = PartitionAlgebra(2, x, R)
            sage: O2 = P2.orbit_basis()
            sage: op = O2([]); op
            O{{-2, 2}, {-1, 1}}
            sage: PA = O2.diagram_basis()
            sage: P2 == PA
            True
            sage: PA([]) == P2.one()
            True
            sage: PA(op)
            -P{{-2, -1, 1, 2}} + P{{-2, 2}, {-1, 1}}
            sage: op == PA(op).to_orbit_basis()
            True
        """
        return self._alg

    def _diagram_to_orbit_on_basis(self, diag):
        """
        Return the element ``diag`` in the orbit basis.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: P2 = PartitionAlgebra(2, x, R)
            sage: O2 = P2.orbit_basis()
            sage: from sage.combinat.diagram_algebras import PartitionDiagrams
            sage: PD = PartitionDiagrams(2)
            sage: O2._diagram_to_orbit_on_basis(PD([[1,2,-2],[-1]]))
            O{{-2, -1, 1, 2}} + O{{-2, 1, 2}, {-1}}
            sage: P2.one().to_orbit_basis()
            O{{-2, -1, 1, 2}} + O{{-2, 2}, {-1, 1}}
            sage: pp = P2[{-2}, {-1, 1}, {2}]
            sage: O2(pp)
            O{{-2}, {-1, 1}, {2}} + O{{-2}, {-1, 1, 2}} + O{{-2, -1, 1}, {2}}
             + O{{-2, -1, 1, 2}} + O{{-2, 2}, {-1, 1}}

        TESTS::

            sage: P2([]).to_orbit_basis() == O2.one()
            True
            sage: O2([]) == O2.one()
            False
            sage: op = O2.an_element()
            sage: op == op.to_diagram_basis().to_orbit_basis()
            True
        """
        PDs = PartitionDiagrams(self._alg._k)
        one = self.base_ring().one()
        return self._from_dict({PDs(d): one for d in diag.set_partition().coarsenings()},
                               coerce=False, remove_zeros=False)

    def product_on_basis(self, d1, d2):
        r"""
        Return the product `O_{d_1} O_{d_2}` of two elements
        in the orbit basis ``self``.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: OP = PartitionAlgebra(2, x, R).orbit_basis()
            sage: SP = OP.basis().keys(); sp = SP([[-2, -1, 1, 2]])
            sage: OP.product_on_basis(sp, sp)
            O{{-2, -1, 1, 2}}
            sage: o1 = OP.one(); o2 = OP([]); o3 = OP.an_element()
            sage: o2 == o1
            False
            sage: o1 * o1 == o1
            True
            sage: o3 * o1 == o1 * o3 and o3 * o1 == o3
            True
            sage: o4 = (3*OP([[-2, -1, 1], [2]]) + 2*OP([[-2, -1, 1, 2]])
            ....:       + 2*OP([[-2, -1, 2], [1]]))
            sage: o4 * o4
            6*O{{-2, -1, 1}, {2}} + 4*O{{-2, -1, 1, 2}} + 4*O{{-2, -1, 2}, {1}}

        We compute Examples 4.5 in [BH2017]_::

            sage: R.<x> = QQ[]
            sage: P = PartitionAlgebra(3,x); O = P.orbit_basis()
            sage: O[[1,2,3],[-1,-2,-3]] * O[[1,2,3],[-1,-2,-3]]
            (x-2)*O{{-3, -2, -1}, {1, 2, 3}} + (x-1)*O{{-3, -2, -1, 1, 2, 3}}

            sage: P = PartitionAlgebra(4,x); O = P.orbit_basis()
            sage: O[[1],[-1],[2,3],[4,-2],[-3,-4]] * O[[1],[2,-2],[3,4],[-1,-3],[-4]]
            (x^2-11*x+30)*O{{-4}, {-3, -1}, {-2, 4}, {1}, {2, 3}}
             + (x^2-9*x+20)*O{{-4}, {-3, -1, 1}, {-2, 4}, {2, 3}}
             + (x^2-9*x+20)*O{{-4}, {-3, -1, 2, 3}, {-2, 4}, {1}}
             + (x^2-9*x+20)*O{{-4, 1}, {-3, -1}, {-2, 4}, {2, 3}}
             + (x^2-7*x+12)*O{{-4, 1}, {-3, -1, 2, 3}, {-2, 4}}
             + (x^2-9*x+20)*O{{-4, 2, 3}, {-3, -1}, {-2, 4}, {1}}
             + (x^2-7*x+12)*O{{-4, 2, 3}, {-3, -1, 1}, {-2, 4}}

            sage: O[[1,-1],[2,-2],[3],[4,-3],[-4]] * O[[1,-2],[2],[3,-1],[4],[-3],[-4]]
            (x-6)*O{{-4}, {-3}, {-2, 1}, {-1, 4}, {2}, {3}}
             + (x-5)*O{{-4}, {-3, 3}, {-2, 1}, {-1, 4}, {2}}
             + (x-5)*O{{-4, 3}, {-3}, {-2, 1}, {-1, 4}, {2}}

            sage: P = PartitionAlgebra(6,x); O = P.orbit_basis()
            sage: (O[[1,-2,-3],[2,4],[3,5,-6],[6],[-1],[-4,-5]]
            ....:  * O[[1,-2],[2,3],[4],[5],[6,-4,-5,-6],[-1,-3]])
            0

            sage: (O[[1,-2],[2,-3],[3,5],[4,-5],[6,-4],[-1],[-6]]
            ....:  * O[[1,-2],[2,-1],[3,-4],[4,-6],[5,-3],[6,-5]])
            O{{-6, 6}, {-5}, {-4, 2}, {-3, 4}, {-2}, {-1, 1}, {3, 5}}

        TESTS:

        Check that multiplication agrees with the multiplication in the
        partition algebra::

            sage: R.<x> = QQ[]
            sage: OP = PartitionAlgebra(2, x).orbit_basis()
            sage: P = OP.diagram_basis()
            sage: o1 = OP.one(); o2 = OP([]); o3 = OP.an_element()
            sage: p1 = P(o1); p2 = P(o2); p3 = P(o3)
            sage: (p2 * p3).to_orbit_basis() == o2 * o3
            True
            sage: (3*p3 * (p1 - 2*p2)).to_orbit_basis() == 3*o3 * (o1 - 2*o2)
            True

            sage: R.<x> = QQ[]
            sage: P = PartitionAlgebra(2,x); O = P.orbit_basis()
            sage: all(b * bp == OP(P(b) * P(bp)) for b in OP.basis() # long time
            ....:     for bp in OP.basis())
            True

        REFERENCES:

        - [BH2017]_
        """
        # According to Corollary 4.12 in [BH2017]_, product is zero unless the
        # stacked diagrams "exactly match" in the middle.
        pi_1 = [frozenset([-i for i in part if i < 0]) for part in d1]
        pi_2 = [frozenset([i for i in part if i > 0]) for part in d2]
        if set([part for part in pi_1 if part]) != set([part for part in pi_2 if part]):
            return self.zero()

        q = self._q
        R = q.parent()
        PDs = self._base_diagrams

        def matchings(A, B):
            for i in range(min(len(A), len(B)) + 1):
                for X in itertools.combinations(A, i):
                    restA = list(A.difference(X))
                    for Y in itertools.combinations(B, i):
                        restB = list(B.difference(Y))
                        for sigma in Permutations(Y):
                            yield [x.union(y) for x, y in zip(X, sigma)] + restA + restB

        D, removed = d1.compose(d2, check=False)
        only_top = set([frozenset(part) for part in d1
                        if all(i > 0 for i in part)])
        only_bottom = set([frozenset(part) for part in d2
                           if all(i < 0 for i in part)])
        only_both = only_top.union(only_bottom)
        restD = [P for P in D if frozenset(P) not in only_both]
        term_dict = {PDs(restD + X):
                     R.prod(q - t for t in range(len(X) + len(restD),
                                                 len(X) + len(restD) + removed))
                     for X in matchings(only_top, only_bottom)}
        return self._from_dict(term_dict)

    class Element(PartitionAlgebra.Element):
        def to_diagram_basis(self):
            """
            Expand ``self`` in the natural diagram basis of the
            partition algebra.

            EXAMPLES::

                sage: R.<x> = QQ[]
                sage: P = PartitionAlgebra(2, x, R)
                sage: O = P.orbit_basis()
                sage: elt = O.an_element(); elt
                3*O{{-2}, {-1, 1, 2}} + 2*O{{-2, -1, 1, 2}} + 2*O{{-2, 1, 2}, {-1}}
                sage: elt.to_diagram_basis()
                3*P{{-2}, {-1, 1, 2}} - 3*P{{-2, -1, 1, 2}} + 2*P{{-2, 1, 2}, {-1}}
                sage: pp = P.an_element(); pp
                3*P{{-2}, {-1, 1, 2}} + 2*P{{-2, -1, 1, 2}} + 2*P{{-2, 1, 2}, {-1}}
                sage: op = pp.to_orbit_basis(); op
                3*O{{-2}, {-1, 1, 2}} + 7*O{{-2, -1, 1, 2}} + 2*O{{-2, 1, 2}, {-1}}
                sage: pp == op.to_diagram_basis()
                True
            """
            # use _coerce_map_from_
            return self.parent()._alg(self)
            #return self._alg.coerce_map_from(self)


class SubPartitionAlgebra(DiagramBasis):
    """
    A subalgebra of the partition algebra in the diagram basis indexed
    by a subset of the diagrams.
    """

    def __init__(self, k, q, base_ring, prefix, diagrams, category=None):
        """
        Initialize ``self`` by adding a coercion to the ambient space.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: BA = BrauerAlgebra(2, x, R)
            sage: BA.ambient().has_coerce_map_from(BA)
            True
        """
        DiagramBasis.__init__(self, k, q, base_ring, prefix, diagrams, category)

    #These methods allow for a subalgebra to be correctly identified in a partition algebra
    def ambient(self):
        r"""
        Return the partition algebra ``self`` is a sub-algebra of.

        EXAMPLES::

            sage: x = var('x')                                                          # needs sage.symbolic
            sage: BA = BrauerAlgebra(2, x)                                              # needs sage.symbolic
            sage: BA.ambient()                                                          # needs sage.symbolic
            Partition Algebra of rank 2 with parameter x over Symbolic Ring
        """
        return self.lift.codomain()

    @lazy_attribute
    def lift(self):
        r"""
        Return the lift map from diagram subalgebra to the ambient space.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: BA = BrauerAlgebra(2, x, R)
            sage: E = BA([[1,2],[-1,-2]])
            sage: lifted = BA.lift(E); lifted
            B{{-2, -1}, {1, 2}}
            sage: lifted.parent() is BA.ambient()
            True
        """
        amb = PartitionAlgebra(self._k, self._q, self.base_ring(), prefix=self._prefix)
        phi = self.module_morphism(amb.monomial, codomain=amb, category=self.category())
        phi.register_as_coercion()
        return phi

    def retract(self, x):
        r"""
        Retract an appropriate partition algebra element to the
        corresponding element in the partition subalgebra.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: BA = BrauerAlgebra(2, x, R)
            sage: PA = BA.ambient()
            sage: E = PA([[1,2], [-1,-2]])
            sage: BA.retract(E) in BA
            True
        """
        if ( x not in self.ambient()
                or any(i not in self._indices for i in x.support()) ):
            raise ValueError("{0} cannot retract to {1}".format(x, self))
        return self._from_dict(x._monomial_coefficients, remove_zeros=False)

    class Element(DiagramBasis.Element):
        def to_orbit_basis(self):
            """
            Return ``self`` in the orbit basis of the associated
            ambient partition algebra.

            EXAMPLES::

                sage: R.<x> = QQ[]
                sage: B = BrauerAlgebra(2, x, R)
                sage: bb = B([[-2, -1], [1, 2]]); bb
                B{{-2, -1}, {1, 2}}
                sage: bb.to_orbit_basis()
                O{{-2, -1}, {1, 2}} + O{{-2, -1, 1, 2}}
            """
            P = self.parent().lift.codomain()
            OP = P.orbit_basis()
            return OP(P(self))


class BrauerAlgebra(SubPartitionAlgebra, UnitDiagramMixin):
    r"""
    A Brauer algebra.

    The Brauer algebra of rank `k` is an algebra with basis indexed by the
    collection of set partitions of `\{1, \ldots, k, -1, \ldots, -k\}`
    with block size 2.

    This algebra is a subalgebra of the partition algebra.
    For more information, see :class:`PartitionAlgebra`.

    INPUT:

    - ``k`` -- rank of the algebra

    - ``q`` -- the deformation parameter `q`

    - ``base_ring`` -- (default: ``None``) a ring containing ``q``; if ``None``
      then just takes the parent of ``q``

    - ``prefix`` -- (default: ``'B'``) a label for the basis elements

    EXAMPLES:

    We now define the Brauer algebra of rank `2` with parameter ``x``
    over `\ZZ`::

        sage: R.<x> = ZZ[]
        sage: B = BrauerAlgebra(2, x, R)
        sage: B
        Brauer Algebra of rank 2 with parameter x
         over Univariate Polynomial Ring in x over Integer Ring
        sage: B.basis()
        Lazy family (Term map from Brauer diagrams of order 2 to Brauer Algebra
         of rank 2 with parameter x over Univariate Polynomial Ring in x
         over Integer Ring(i))_{i in Brauer diagrams of order 2}
        sage: B.basis().keys()
        Brauer diagrams of order 2
        sage: B.basis().keys()([[-2, 1], [2, -1]])
        {{-2, 1}, {-1, 2}}
        sage: b = B.basis().list(); b
        [B{{-2, -1}, {1, 2}}, B{{-2, 1}, {-1, 2}}, B{{-2, 2}, {-1, 1}}]
        sage: b[0]
        B{{-2, -1}, {1, 2}}
        sage: b[0]^2
        x*B{{-2, -1}, {1, 2}}
        sage: b[0]^5
        x^4*B{{-2, -1}, {1, 2}}

    Note, also that since the symmetric group algebra is contained in
    the Brauer algebra, there is also a conversion between the two. ::

        sage: R.<x> = ZZ[]
        sage: B = BrauerAlgebra(2, x, R)
        sage: S = SymmetricGroupAlgebra(R, 2)
        sage: S([2,1])*B([[1,-1],[2,-2]])
        B{{-2, 1}, {-1, 2}}
    """

    @staticmethod
    def __classcall_private__(cls, k, q, base_ring=None, prefix='B'):
        r"""
        Standardize the input by getting the base ring from the parent of
        the parameter ``q`` if no ``base_ring`` is given.

        TESTS::

            sage: R.<q> = QQ[]
            sage: BA1 = BrauerAlgebra(2, q)
            sage: BA2 = BrauerAlgebra(2, q, R, 'B')
            sage: BA1 is BA2
            True
        """
        if base_ring is None:
            base_ring = q.parent()
        return super().__classcall__(cls, k, q, base_ring, prefix)

    def __init__(self, k, q, base_ring, prefix):
        r"""
        Initialize ``self``.

        TESTS::

            sage: R.<q> = QQ[]
            sage: BA = BrauerAlgebra(2, q, R)
            sage: TestSuite(BA).run()
        """
        SubPartitionAlgebra.__init__(self, k, q, base_ring, prefix, BrauerDiagrams(k))

    options = BrauerDiagram.options

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: R.<q> = QQ[]
            sage: BrauerAlgebra(2, q, R)
            Brauer Algebra of rank 2 with parameter q
             over Univariate Polynomial Ring in q over Rational Field
        """
        return "Brauer Algebra of rank {} with parameter {} over {}".format(
                self._k, self._q, self.base_ring())

    # TODO: Make a mixin class for diagram algebras that have coercions from SGA?
    def _coerce_map_from_(self, R):
        """
        Return a coerce map from ``R`` if one exists and ``None`` otherwise.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: S = SymmetricGroupAlgebra(R, 4)
            sage: A = BrauerAlgebra(4, x, R)
            sage: A._coerce_map_from_(S)
            Generic morphism:
              From: Symmetric group algebra of order 4 over Univariate Polynomial Ring in x over Rational Field
              To:   Brauer Algebra of rank 4 with parameter x over Univariate Polynomial Ring in x over Rational Field
            sage: Sp = SymmetricGroupAlgebra(QQ, 4)
            sage: A._coerce_map_from_(Sp)
            Generic morphism:
              From: Symmetric group algebra of order 4 over Rational Field
              To:   Brauer Algebra of rank 4 with parameter x over Univariate Polynomial Ring in x over Rational Field
            sage: Sp3 = SymmetricGroupAlgebra(QQ, 3)
            sage: A._coerce_map_from_(Sp3)
            Generic morphism:
              From: Symmetric group algebra of order 3 over Rational Field
              To:   Brauer Algebra of rank 4 with parameter x over Univariate Polynomial Ring in x over Rational Field
        """
        if isinstance(R, SymmetricGroupAlgebra_n):
            if R.n <= self._k and self.base_ring().has_coerce_map_from(R.base_ring()):
                return R.module_morphism(self._perm_to_Blst, codomain=self)
            return None
        return super()._coerce_map_from_(R)

    def _element_constructor_(self, set_partition):
        r"""
        Construct an element of ``self``.

        EXAMPLES::

            sage: R.<q> = QQ[]
            sage: BA = BrauerAlgebra(2, q, R)
            sage: sp = SetPartition([[1,2], [-1,-2]])
            sage: b_elt = BA(sp); b_elt
            B{{-2, -1}, {1, 2}}
            sage: b_elt in BA
            True
            sage: BA([[1,2],[-1,-2]]) == b_elt
            True
            sage: BA([{1,2},{-1,-2}]) == b_elt
            True
        """
        set_partition = to_Brauer_partition(set_partition, k=self.order())
        return DiagramAlgebra._element_constructor_(self, set_partition)

    def jucys_murphy(self, j):
        r"""
        Return the ``j``-th generalized Jucys-Murphy element of ``self``.

        The `j`-th Jucys-Murphy element of a Brauer algebra is simply
        the `j`-th Jucys-Murphy element of the symmetric group algebra
        with an extra `(z-1)/2` term, where ``z`` is the parameter
        of the Brauer algebra.

        REFERENCES:

        .. [Naz96] Maxim Nazarov, Young's Orthogonal Form for Brauer's
           Centralizer Algebra. Journal of Algebra 182 (1996), 664--693.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: z = var('z')
            sage: B = BrauerAlgebra(3,z)
            sage: B.jucys_murphy(1)
            (1/2*z-1/2)*B{{-3, 3}, {-2, 2}, {-1, 1}}
            sage: B.jucys_murphy(3)
            -B{{-3, -2}, {-1, 1}, {2, 3}} - B{{-3, -1}, {-2, 2}, {1, 3}}
             + B{{-3, 1}, {-2, 2}, {-1, 3}} + B{{-3, 2}, {-2, 3}, {-1, 1}}
             + (1/2*z-1/2)*B{{-3, 3}, {-2, 2}, {-1, 1}}
        """
        if j < 1:
            raise ValueError("Jucys-Murphy index must be positive")
        k = self.order()
        if j > k:
            raise ValueError("Jucys-Murphy index cannot be greater than the order of the algebra")

        def convertI(x):
            return self._indices(to_Brauer_partition(x, k=k))
        R = self.base_ring()
        one = R.one()
        d = {self.one_basis(): R((self._q - 1)/2)}
        for i in range(1, j):
            d[convertI([[i, -j], [j, -i]])] = one
            d[convertI([[i, j], [-i, -j]])] = -one
        return self._from_dict(d, remove_zeros=True)


class HalfTemperleyLiebDiagrams(UniqueRepresentation, Parent):
    r"""
    Half diagrams for the Temperley-Lieb algebra cell modules.
    """
    def __init__(self, order, defects):
        r"""
        Initialize ``self``.

        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: htld = da.HalfTemperleyLiebDiagrams(7, 3)
            sage: TestSuite(htld).run()
        """
        Parent.__init__(self, category=FiniteEnumeratedSets())
        self._order = ZZ(order)
        self._defects = ZZ(defects)
        if (self._order - self._defects) % 2:
            raise ValueError("the number of non-defects must be even")

    def _repr_(self):
        r"""
        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: da.HalfTemperleyLiebDiagrams(7, 5)
            Half Temperley-Lieb diagrams of order 7 with 5 defects
        """
        return "Half Temperley-Lieb diagrams of order {} with {} defects".format(self._order, self._defects)

    def __iter__(self):
        r"""
        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: list(da.HalfTemperleyLiebDiagrams(5, 3))
            [{{1, 2}}, {{2, 3}}, {{3, 4}}, {{4, 5}}]
        """
        n = self._order
        k = self._defects
        b = (n - k) // 2
        from sage.combinat.dyck_word import DyckWords
        for dw in DyckWords(b+k, b):
            ret = []
            offset = 0
            for D in dw.catalan_factorization():
                ret.extend((offset+a+1, offset+b) for (a, b) in D.tunnels())
                offset += len(D) + 1
            yield self.element_class(self, ret)

    def __contains__(self, obj):
        r"""
        Check containment.

        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: htld = da.HalfTemperleyLiebDiagrams(7, 3)
            sage: htld.an_element() in htld
            True
        """
        return isinstance(obj, self.Element) and obj.parent() is self

    def cardinality(self):
        r"""
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: import sage.combinat.diagram_algebras as da
            sage: htld = da.HalfTemperleyLiebDiagrams(7, 3)
            sage: htld.cardinality()
            14
        """
        from sage.functions.other import binomial
        n = self._order
        k = self._defects
        b = (n - k) // 2
        return (k + 1) * binomial(n, b) // (b + k + 1)

    def _element_constructor_(self, d):
        r"""
        Construct an element of ``self``.

        EXAMPLES::

            sage: import sage.combinat.diagram_algebras as da
            sage: htld = da.HalfTemperleyLiebDiagrams(7, 3)
            sage: htld([[1, 4], [2, 3]])
            {{1, 4}, {2, 3}}
        """
        return self.element_class(self, d)

    class Element(AbstractPartitionDiagram):
        def _latex_(self):
            r"""
            Return a latex representation of ``self``.

            EXAMPLES::

                sage: import sage.combinat.diagram_algebras as da
                sage: htld = da.HalfTemperleyLiebDiagrams(7, 3)
                sage: d = htld([[1, 2], [4, 5]])
                sage: latex(d)
                \begin{tikzpicture}[scale = 0.5,thick, baseline={(0,-1ex/2)}]
                \tikzstyle{vertex} = [shape = circle, minimum size = 7pt, inner sep = 1pt]
                \node[vertex] (G--7) at (9.0, -1) [shape = circle, draw] {};
                \node[vertex] (G--7) at (9.0, -1) [shape = circle, draw] {};
                \node[vertex] (G--6) at (7.5, -1) [shape = circle, draw] {};
                \node[vertex] (G--6) at (7.5, -1) [shape = circle, draw] {};
                \node[vertex] (G--5) at (6.0, -1) [shape = circle, draw] {};
                \node[vertex] (G--4) at (4.5, -1) [shape = circle, draw] {};
                \node[vertex] (G--3) at (3.0, -1) [shape = circle, draw] {};
                \node[vertex] (G--3) at (3.0, -1) [shape = circle, draw] {};
                \node[vertex] (G--2) at (1.5, -1) [shape = circle, draw] {};
                \node[vertex] (G--1) at (0.0, -1) [shape = circle, draw] {};
                \draw[] (G--7) .. controls +(0.0, 0.4) and +(-0.0, 0.4) .. (G--7);
                \draw[] (G--6) .. controls +(0.0, 0.4) and +(-0.0, 0.4) .. (G--6);
                \draw[] (G--5) .. controls +(-0.5, 0.5) and +(0.5, 0.5) .. (G--4);
                \draw[] (G--3) .. controls +(0.0, 0.4) and +(-0.0, 0.4) .. (G--3);
                \draw[] (G--2) .. controls +(-0.5, 0.5) and +(0.5, 0.5) .. (G--1);
                \end{tikzpicture}
            """
            defects = self.defects()
            temp = [[-b, -a] for (a, b) in self] + [[-d, -d] for d in defects]
            temp.sort()
            return diagram_latex(temp)

        def _ascii_art_(self):
            r"""
            Return an ascii art representation of ``self``.

            EXAMPLES::

                sage: import sage.combinat.diagram_algebras as da
                sage: htld = da.HalfTemperleyLiebDiagrams(7, 3)
                sage: d = htld([[1, 2], [4, 5]])
                sage: ascii_art(d)
                 .-. | .-. | |
                 o o o o o o o
                sage: htld = da.HalfTemperleyLiebDiagrams(8, 0)
                sage: d = htld([[1, 6], [2, 3], [4, 5], [7, 8]])
                sage: ascii_art(d)
                 .---------.
                 | .-. .-. | .-.
                 o o o o o o o o
            """
            defects = self.defects()
            temp = [[-b, -a] for (a, b) in self] + [[-d, d] for d in defects]
            rank = self.parent()._order
            temp.append([rank, rank])
            temp.sort()
            ret = TL_diagram_ascii_art(temp)
            from sage.typeset.ascii_art import AsciiArt
            return AsciiArt(ret[2:])

        def _unicode_art_(self):
            r"""
            Return an ascii art representation of ``self``.

            EXAMPLES::

                sage: import sage.combinat.diagram_algebras as da
                sage: htld = da.HalfTemperleyLiebDiagrams(7, 3)
                sage: d = htld([[1, 2], [4, 5]])
                sage: unicode_art(d)
                 ╭─╮ │ ╭─╮ │ │
                 ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬
                sage: htld = da.HalfTemperleyLiebDiagrams(8, 0)
                sage: d = htld([[1, 6], [2, 3], [4, 5], [7, 8]])
                sage: unicode_art(d)
                 ╭─────────╮
                 │ ╭─╮ ╭─╮ │ ╭─╮
                 ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬
            """
            defects = self.defects()
            temp = [[-b, -a] for (a, b) in self] + [[-d, d] for d in defects]
            rank = self.parent()._order
            temp.append([rank, rank])
            temp.sort()
            ret = TL_diagram_ascii_art(temp, use_unicode=True)
            from sage.typeset.unicode_art import UnicodeArt
            return UnicodeArt(ret[2:])

        def defects(self):
            r"""
            Return the defects of ``self``.

            EXAMPLES::

                sage: import sage.combinat.diagram_algebras as da
                sage: htld = da.HalfTemperleyLiebDiagrams(7, 3)
                sage: d = htld([[1, 2], [4, 5]])
                sage: d.defects()
                frozenset({3, 6, 7})
            """
            order = self.parent()._order
            return frozenset(range(1, order+1)) - frozenset(e for B in self for e in B)

        def check(self):
            r"""
            Check the validity of the input of ``self``.

            EXAMPLES::

                sage: import sage.combinat.diagram_algebras as da
                sage: htld = da.HalfTemperleyLiebDiagrams(7, 3)
                sage: htld([[1,2], [3,4]])  # indirect doctest
                {{1, 2}, {3, 4}}
                sage: htld([[1,2], [-1, -2]])  # indirect doctest
                Traceback (most recent call last):
                ...
                ValueError: {{-2, -1}, {1, 2}} does not represent a half TL diagram of order 7
                sage: htld([[1,2,3], [4,5]])  # indirect doctest
                Traceback (most recent call last):
                ...
                ValueError: all blocks of {{1, 2, 3}, {4, 5}} must be of size 2
                sage: htld([[1,2], [3,4], [5,6]])  # indirect doctest
                Traceback (most recent call last):
                ...
                ValueError: {{1, 2}, {3, 4}, {5, 6}} does not have 3 defects
                sage: htld([[1,3], [2,4]])  # indirect doctest
                Traceback (most recent call last):
                ...
                ValueError: {{1, 3}, {2, 4}} is not planar
            """
            tst = frozenset(e for B in self._base_diagram for e in B)
            P = self.parent()
            if not (tst <= frozenset(range(1, P._order+1))):
                raise ValueError("{} does not represent a half TL diagram of order {}".format(
                                     self, self.parent()._order))
            if any(len(block) != 2 for block in self):
                raise ValueError("all blocks of {} must be of size 2".format(self))
            if len(tst) != P._order - P._defects:
                raise ValueError("{} does not have {} defects".format(self, P._defects))
            if not self.is_planar():
                raise ValueError("{} is not planar".format(self))


class TemperleyLiebAlgebra(SubPartitionAlgebra, UnitDiagramMixin):
    r"""
    A Temperley--Lieb algebra.

    The Temperley--Lieb algebra of rank `k` is an algebra with basis
    indexed by the collection of planar set partitions of
    `\{1, \ldots, k, -1, \ldots, -k\}` with block size 2.

    This algebra is thus a subalgebra of the partition algebra.
    For more information, see :class:`PartitionAlgebra`.

    INPUT:

    - ``k`` -- rank of the algebra

    - ``q`` -- the deformation parameter `q`

    - ``base_ring`` -- (default: ``None``) a ring containing ``q``; if ``None``
      then just takes the parent of ``q``

    - ``prefix`` -- (default: ``'T'``) a label for the basis elements

    EXAMPLES:

    We define the Temperley--Lieb algebra of rank `2` with parameter
    `x` over `\ZZ`::

        sage: R.<x> = ZZ[]
        sage: T = TemperleyLiebAlgebra(2, x, R); T
        Temperley-Lieb Algebra of rank 2 with parameter x
         over Univariate Polynomial Ring in x over Integer Ring
        sage: T.basis()
        Lazy family (Term map from Temperley Lieb diagrams of order 2
         to Temperley-Lieb Algebra of rank 2 with parameter x over
         Univariate Polynomial Ring in x over Integer
         Ring(i))_{i in Temperley Lieb diagrams of order 2}
        sage: T.basis().keys()
        Temperley Lieb diagrams of order 2
        sage: T.basis().keys()([[-1, 1], [2, -2]])
        {{-2, 2}, {-1, 1}}
        sage: b = T.basis().list(); b
        [T{{-2, -1}, {1, 2}}, T{{-2, 2}, {-1, 1}}]
        sage: b[0]
        T{{-2, -1}, {1, 2}}
        sage: b[0]^2 == x*b[0]
        True
        sage: b[0]^5 == x^4*b[0]
        True

    The Temperley-Lieb algebra is a cellular algebra, and we verify that
    the dimensions of the simple modules at `q = 0` is given by
    :oeis:`A050166`::

        sage: for k in range(1,5):
        ....:     TL = TemperleyLiebAlgebra(2*k, 0, QQ)
        ....:     print("".join("{:3}".format(TL.cell_module(la).simple_module().dimension())
        ....:                   for la in reversed(TL.cell_poset()) if la != 0))
          1
          1  2
          1  4  5
          1  6 14 14
        sage: for k in range(1,4):
        ....:     TL = TemperleyLiebAlgebra(2*k+1, 0, QQ)
        ....:     print("".join("{:3}".format(TL.cell_module(la).simple_module().dimension())
        ....:                   for la in reversed(TL.cell_poset()) if la != 0))
          1  2
          1  4  5
          1  6 14 14

    Additional examples when the Temperley-Lieb algebra is not semisimple::

        sage: TL = TemperleyLiebAlgebra(8, -1, QQ)
        sage: for la in TL.cell_poset():
        ....:     CM = TL.cell_module(la)
        ....:     if not CM.nonzero_bilinear_form():
        ....:         continue
        ....:     print(la, CM.dimension(), CM.simple_module().dimension())
        ....:
        0 14 1
        2 28 28
        4 20 13
        6 7 7
        8 1 1
        sage: for k in range(1,5):
        ....:     TL = TemperleyLiebAlgebra(2*k, -1, QQ)
        ....:     print("".join("{:3}".format(TL.cell_module(la).simple_module().dimension())
        ....:                   for la in reversed(TL.cell_poset())
        ....:                    if TL.cell_module(la).nonzero_bilinear_form()))
          1  1
          1  3  1
          1  4  9  1
          1  7 13 28  1
        sage: C5.<z5> = CyclotomicField(5)
        sage: for k in range(1,5):
        ....:     TL = TemperleyLiebAlgebra(2*k, z5+~z5, C5)
        ....:     print("".join("{:3}".format(TL.cell_module(la).simple_module().dimension())
        ....:                   for la in reversed(TL.cell_poset())
        ....:                    if TL.cell_module(la).nonzero_bilinear_form()))
          1  1
          1  3  2
          1  5  8  5
          1  7 20 21 13
    """
    @staticmethod
    def __classcall_private__(cls, k, q, base_ring=None, prefix='T'):
        r"""
        Standardize the input by getting the base ring from the parent of
        the parameter ``q`` if no ``base_ring`` is given.

        TESTS::

            sage: R.<q> = QQ[]
            sage: T1 = TemperleyLiebAlgebra(2, q)
            sage: T2 = TemperleyLiebAlgebra(2, q, R, 'T')
            sage: T1 is T2
            True
        """
        if base_ring is None:
            base_ring = q.parent()
        return super().__classcall__(cls, k, q, base_ring, prefix)

    def __init__(self, k, q, base_ring, prefix):
        r"""
        Initialize ``self``.

        TESTS::

            sage: R.<q> = QQ[]
            sage: TL = TemperleyLiebAlgebra(2, q, R)
            sage: TestSuite(TL).run()

            sage: TL = TemperleyLiebAlgebra(3, 0, QQ)
            sage: TestSuite(TL).run()
        """
        cat = AssociativeAlgebras(base_ring.category()).Unital().FiniteDimensional().WithBasis().Cellular()
        SubPartitionAlgebra.__init__(self, k, q, base_ring, prefix, TemperleyLiebDiagrams(k), category=cat)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: R.<q> = QQ[]
            sage: TemperleyLiebAlgebra(2, q, R)
            Temperley-Lieb Algebra of rank 2 with parameter q
             over Univariate Polynomial Ring in q over Rational Field
        """
        return "Temperley-Lieb Algebra of rank {} with parameter {} over {}".format(
                self._k, self._q, self.base_ring())

    def _element_constructor_(self, set_partition):
        r"""
        Construct an element of ``self``.

        EXAMPLES::

            sage: R.<q> = QQ[]
            sage: TL = TemperleyLiebAlgebra(2, q, R)
            sage: sp = SetPartition([[1,2], [-1,-2]])
            sage: b_elt = TL(sp); b_elt
            T{{-2, -1}, {1, 2}}
            sage: b_elt in TL
            True
            sage: TL([[1,2],[-1,-2]]) == b_elt
            True
            sage: TL([{1,2},{-1,-2}]) == b_elt
            True
            sage: S = SymmetricGroupAlgebra(R, 2)
            sage: TL(S([1,2]))
            T{{-2, 2}, {-1, 1}}
            sage: TL(S([2,1]))
            Traceback (most recent call last):
            ...
            ValueError: the diagram {{-2, 1}, {-1, 2}} must be planar
        """
        if isinstance(set_partition, SymmetricGroupAlgebra_n.Element):
            return SubPartitionAlgebra._element_constructor_(self, set_partition)
        set_partition = to_Brauer_partition(set_partition, k=self.order())
        return SubPartitionAlgebra._element_constructor_(self, set_partition)

    def _ascii_art_term(self, diagram):
        r"""
        Return an ascii art representation of ``diagram``.

        EXAMPLES::

            sage: R.<q> = QQ[]
            sage: TL = TemperleyLiebAlgebra(4, q, R)
            sage: x = TL.an_element()
            sage: ascii_art(x)  # indirect doctest
                            o o o o       o o o o
               o o o o      | `-` |       | `-` |
            2* `-` `-` + 2* `-----`  + 3* `---. |
               .-. .-.      .-. .-.       .-. | |
               o o o o      o o o o       o o o o
        """
        return TL_diagram_ascii_art(diagram, use_unicode=False)

    def _unicode_art_term(self, diagram):
        r"""
        Return a unicode art representation of ``diagram``.

        EXAMPLES::

            sage: R.<q> = QQ[]
            sage: TL = TemperleyLiebAlgebra(4, q, R)
            sage: x = TL.an_element()
            sage: unicode_art(x)  # indirect doctest
                            ⚬ ⚬ ⚬ ⚬       ⚬ ⚬ ⚬ ⚬
               ⚬ ⚬ ⚬ ⚬      │ ╰─╯ │       │ ╰─╯ │
            2* ╰─╯ ╰─╯ + 2* ╰─────╯  + 3* ╰───╮ │
               ╭─╮ ╭─╮      ╭─╮ ╭─╮       ╭─╮ │ │
               ⚬ ⚬ ⚬ ⚬      ⚬ ⚬ ⚬ ⚬       ⚬ ⚬ ⚬ ⚬
        """
        return TL_diagram_ascii_art(diagram, use_unicode=True)

    @cached_method
    def cell_poset(self):
        """
        Return the cell poset of ``self``.

        EXAMPLES::

            sage: R.<q> = QQ[]
            sage: TL = TemperleyLiebAlgebra(7, q, R)
            sage: TL.cell_poset().cover_relations()
            [[1, 3], [3, 5], [5, 7]]

            sage: TL = TemperleyLiebAlgebra(8, q, R)
            sage: TL.cell_poset().cover_relations()
            [[0, 2], [2, 4], [4, 6], [6, 8]]
        """
        from sage.combinat.posets.posets import Poset
        return Poset({k-2: [k] for k in range(self._k, 1, -2)})

    def cell_module_indices(self, la):
        r"""
        Return the indices of the cell module of ``self``
        indexed by ``la`` .

        This is the finite set `M(\lambda)`.

        EXAMPLES::

            sage: R.<q> = QQ[]
            sage: TL = TemperleyLiebAlgebra(8, q, R)
            sage: TL.cell_module_indices(4)
            Half Temperley-Lieb diagrams of order 8 with 4 defects
        """
        return HalfTemperleyLiebDiagrams(self._k, la)

    def _to_cellular_element(self, d):
        r"""
        Return the image in the cellular basis of the basis element
        of ``self`` indexed by ``d``.

        EXAMPLES::

            sage: R.<q> = QQ[]
            sage: TL = TemperleyLiebAlgebra(3, q, R)
            sage: for d in TL.basis().keys():
            ....:     print(d)
            ....:     print(TL._to_cellular_element(d))
            {{-3, 3}, {-2, -1}, {1, 2}}
            C(1, {{1, 2}}, {{1, 2}})
            {{-3, 1}, {-2, -1}, {2, 3}}
            C(1, {{2, 3}}, {{1, 2}})
            {{-3, -2}, {-1, 1}, {2, 3}}
            C(1, {{2, 3}}, {{2, 3}})
            {{-3, -2}, {-1, 3}, {1, 2}}
            C(1, {{1, 2}}, {{2, 3}})
            {{-3, 3}, {-2, 2}, {-1, 1}}
            C(3, {}, {})
        """
        C = self.cellular_basis()
        top = []
        bottom = []
        defects = ZZ.zero()
        for (a, b) in d:
            if b < 0:
                bottom.append((-b, -a))
            elif a > 0:
                top.append((a, b))
            else:
                defects += ZZ.one()
        CMI = self.cell_module_indices(defects)
        tup = (defects, CMI(top), CMI(bottom))
        return C.monomial(C._indices(tup))

    def _from_cellular_index(self, x):
        r"""
        Return the image in ``self`` from the index of the
        cellular basis ``x``.

        EXAMPLES::

            sage: R.<q> = QQ[]
            sage: TL = TemperleyLiebAlgebra(3, q, R)
            sage: C = TL.cellular_basis()
            sage: for i in C.basis().keys():
            ....:     print(i)
            ....:     print(TL._from_cellular_index(i))
            (1, {{1, 2}}, {{1, 2}})
            T{{-3, 3}, {-2, -1}, {1, 2}}
            (1, {{1, 2}}, {{2, 3}})
            T{{-3, -2}, {-1, 3}, {1, 2}}
            (1, {{2, 3}}, {{1, 2}})
            T{{-3, 1}, {-2, -1}, {2, 3}}
            (1, {{2, 3}}, {{2, 3}})
            T{{-3, -2}, {-1, 1}, {2, 3}}
            (3, {}, {})
            T{{-3, 3}, {-2, 2}, {-1, 1}}

            sage: TL = TemperleyLiebAlgebra(4, QQ.zero(), QQ)
            sage: C = TL.cellular_basis()
            sage: [TL._from_cellular_index(i) for i in C.basis().keys()]
            [T{{-4, -3}, {-2, -1}, {1, 2}, {3, 4}},
             T{{-4, -1}, {-3, -2}, {1, 2}, {3, 4}},
             T{{-4, -3}, {-2, -1}, {1, 4}, {2, 3}},
             T{{-4, -1}, {-3, -2}, {1, 4}, {2, 3}},
             T{{-4, 4}, {-3, 3}, {-2, -1}, {1, 2}},
             T{{-4, 4}, {-3, -2}, {-1, 3}, {1, 2}},
             T{{-4, -3}, {-2, 4}, {-1, 3}, {1, 2}},
             T{{-4, 4}, {-3, 1}, {-2, -1}, {2, 3}},
             T{{-4, 4}, {-3, -2}, {-1, 1}, {2, 3}},
             T{{-4, -3}, {-2, 4}, {-1, 1}, {2, 3}},
             T{{-4, 2}, {-3, 1}, {-2, -1}, {3, 4}},
             T{{-4, 2}, {-3, -2}, {-1, 1}, {3, 4}},
             T{{-4, -3}, {-2, 2}, {-1, 1}, {3, 4}},
             T{{-4, 4}, {-3, 3}, {-2, 2}, {-1, 1}}]
        """
        _, top, bottom = x
        bottom = [[-b, -a] for (a, b) in bottom]
        tmiss = frozenset(range(1, self._k+1)) - frozenset(e for B in top for e in B)
        bmiss = frozenset(range(-1, -self._k-1, -1)) - frozenset(e for B in bottom for e in B)
        prop = list(zip(sorted(tmiss, reverse=True), sorted(bmiss)))
        return self.monomial(self._indices(bottom + prop + list(top)))

    def cellular_involution(self, x):
        r"""
        Return the cellular involution of ``x`` in ``self``.

        EXAMPLES::

            sage: TL = TemperleyLiebAlgebra(4, QQ.zero(), QQ)
            sage: ascii_art(TL.an_element())
                            o o o o       o o o o
               o o o o      | `-` |       | `-` |
            2* `-` `-` + 2* `-----`  + 3* `---. |
               .-. .-.      .-. .-.       .-. | |
               o o o o      o o o o       o o o o
            sage: ascii_art(TL.cellular_involution(TL.an_element()))
                            o o o o       o o o o
               o o o o      `-` `-`       `-` | |
            2* `-` `-` + 2* .-----.  + 3* .---` |
               .-. .-.      | .-. |       | .-. |
               o o o o      o o o o       o o o o
        """
        M = x.monomial_coefficients(copy=False)
        I = self._indices
        return self._from_dict({d.dual(): c for d, c in M.items()},
                               remove_zeros=False)


class PlanarAlgebra(SubPartitionAlgebra, UnitDiagramMixin):
    r"""
    A planar algebra.

    The planar algebra of rank `k` is an algebra with basis indexed by the
    collection of all planar set partitions of
    `\{1, \ldots, k, -1, \ldots, -k\}`.

    This algebra is thus a subalgebra of the partition algebra. For more
    information, see :class:`PartitionAlgebra`.

    INPUT:

    - ``k`` -- rank of the algebra

    - ``q`` -- the deformation parameter `q`

    - ``base_ring`` -- (default: ``None``) a ring containing ``q``; if ``None``
      then just takes the parent of ``q``

    - ``prefix`` -- (default: ``'Pl'``) a label for the basis elements

    EXAMPLES:

    We define the planar algebra of rank `2` with parameter
    `x` over `\ZZ`::

        sage: R.<x> = ZZ[]
        sage: Pl = PlanarAlgebra(2, x, R); Pl
        Planar Algebra of rank 2 with parameter x over Univariate Polynomial Ring in x over Integer Ring
        sage: Pl.basis().keys()
        Planar diagrams of order 2
        sage: Pl.basis().keys()([[-1, 1], [2, -2]])
        {{-2, 2}, {-1, 1}}
        sage: Pl.basis().list()
        [Pl{{-2}, {-1}, {1, 2}},
         Pl{{-2}, {-1}, {1}, {2}},
         Pl{{-2, 1}, {-1}, {2}},
         Pl{{-2, 2}, {-1}, {1}},
         Pl{{-2, 1, 2}, {-1}},
         Pl{{-2, 2}, {-1, 1}},
         Pl{{-2}, {-1, 1}, {2}},
         Pl{{-2}, {-1, 2}, {1}},
         Pl{{-2}, {-1, 1, 2}},
         Pl{{-2, -1}, {1, 2}},
         Pl{{-2, -1}, {1}, {2}},
         Pl{{-2, -1, 1}, {2}},
         Pl{{-2, -1, 2}, {1}},
         Pl{{-2, -1, 1, 2}}]
        sage: E = Pl([[1,2],[-1,-2]])
        sage: E^2 == x*E
        True
        sage: E^5 == x^4*E
        True
    """
    @staticmethod
    def __classcall_private__(cls, k, q, base_ring=None, prefix='Pl'):
        r"""
        Standardize the input by getting the base ring from the parent of
        the parameter ``q`` if no ``base_ring`` is given.

        TESTS::

            sage: R.<q> = QQ[]
            sage: Pl1 = PlanarAlgebra(2, q)
            sage: Pl2 = PlanarAlgebra(2, q, R, 'Pl')
            sage: Pl1 is Pl2
            True
        """
        if base_ring is None:
            base_ring = q.parent()
        return super().__classcall__(cls, k, q, base_ring, prefix)

    def __init__(self, k, q, base_ring, prefix):
        r"""
        Initialize ``self``.

        TESTS::

            sage: R.<q> = QQ[]
            sage: PlA = PlanarAlgebra(2, q, R)
            sage: TestSuite(PlA).run()
        """
        SubPartitionAlgebra.__init__(self, k, q, base_ring, prefix, PlanarDiagrams(k))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: Pl = PlanarAlgebra(2, x, R); Pl
            Planar Algebra of rank 2 with parameter x
             over Univariate Polynomial Ring in x over Integer Ring
        """
        txt = "Planar Algebra of rank {} with parameter {} over {}"
        return txt.format(self._k, self._q, self.base_ring())


class PropagatingIdeal(SubPartitionAlgebra):
    r"""
    A propagating ideal.

    The propagating ideal of rank `k` is a non-unital algebra with basis
    indexed by the collection of ideal set partitions of `\{1, \ldots, k, -1,
    \ldots, -k\}`. We say a set partition is *ideal* if its propagating
    number is less than `k`.

    This algebra is a non-unital subalgebra and an ideal of the partition
    algebra.
    For more information, see :class:`PartitionAlgebra`.

    EXAMPLES:

    We now define the propagating ideal of rank `2` with parameter
    `x` over `\ZZ`::

        sage: R.<x> = QQ[]
        sage: I = PropagatingIdeal(2, x, R); I
        Propagating Ideal of rank 2 with parameter x
         over Univariate Polynomial Ring in x over Rational Field
        sage: I.basis().keys()
        Ideal diagrams of order 2
        sage: I.basis().list()
        [I{{-2, -1, 1, 2}},
         I{{-2, 1, 2}, {-1}},
         I{{-2}, {-1, 1, 2}},
         I{{-2, -1}, {1, 2}},
         I{{-2}, {-1}, {1, 2}},
         I{{-2, -1, 1}, {2}},
         I{{-2, 1}, {-1}, {2}},
         I{{-2, -1, 2}, {1}},
         I{{-2, 2}, {-1}, {1}},
         I{{-2}, {-1, 1}, {2}},
         I{{-2}, {-1, 2}, {1}},
         I{{-2, -1}, {1}, {2}},
         I{{-2}, {-1}, {1}, {2}}]
        sage: E = I([[1,2],[-1,-2]])
        sage: E^2 == x*E
        True
        sage: E^5 == x^4*E
        True
    """
    @staticmethod
    def __classcall_private__(cls, k, q, base_ring=None, prefix='I'):
        r"""
        Standardize the input by getting the base ring from the parent of
        the parameter ``q`` if no ``base_ring`` is given.

        TESTS::

            sage: R.<q> = QQ[]
            sage: IA1 = PropagatingIdeal(2, q)
            sage: IA2 = PropagatingIdeal(2, q, R, 'I')
            sage: IA1 is IA2
            True
        """
        if base_ring is None:
            base_ring = q.parent()
        return super().__classcall__(cls, k, q, base_ring, prefix)

    def __init__(self, k, q, base_ring, prefix):
        r"""
        Initialize ``self``.

        TESTS::

            sage: R.<q> = QQ[]
            sage: I = PropagatingIdeal(2, q, R)
            sage: TestSuite(I).run()
        """
        category = AssociativeAlgebras(base_ring.category()).FiniteDimensional().WithBasis()
        SubPartitionAlgebra.__init__(self, k, q, base_ring, prefix,
                                     IdealDiagrams(k), category)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: PropagatingIdeal(2, x, R)
            Propagating Ideal of rank 2 with parameter x over Univariate
             Polynomial Ring in x over Rational Field
        """
        return "Propagating Ideal of rank {} with parameter {} over {}".format(
                self._k, self._q, self.base_ring())

    class Element(SubPartitionAlgebra.Element):
        """
        An element of a propagating ideal.

        We need to take care of exponents since we are not unital.
        """

        def __pow__(self, n):
            """
            Return ``self`` to the `n`-th power.

            INPUT:

            - ``n`` -- positive integer

            EXAMPLES::

                sage: R.<x> = QQ[]
                sage: I = PropagatingIdeal(2, x, R)
                sage: E = I([[1,2],[-1,-2]])
                sage: E^2
                x*I{{-2, -1}, {1, 2}}
                sage: E^0
                Traceback (most recent call last):
                ...
                ValueError: can only take positive integer powers
            """
            if n <= 0:
                raise ValueError("can only take positive integer powers")
            return generic_power(self, n)


def TL_diagram_ascii_art(diagram, use_unicode=False, blobs=[]):
    r"""
    Return ascii art for a Temperley-Lieb diagram ``diagram``.

    INPUT:

    - ``diagram`` -- list of pairs of matchings of the set
      `\{-1, \ldots, -n, 1, \ldots, n\}`
    - ``use_unicode`` -- boolean (default: ``False``); whether or not
      to use unicode art instead of ascii art
    - ``blobs`` -- (optional) a list of matchings with blobs on them

    EXAMPLES::

        sage: from sage.combinat.diagram_algebras import TL_diagram_ascii_art
        sage: TL = [(-15,-12), (-14,-13), (-11,15), (-10,14), (-9,-6),
        ....:       (-8,-7), (-5,-4), (-3,1), (-2,-1), (2,3), (4,5),
        ....:       (6,11), (7, 8), (9,10), (12,13)]
        sage: TL_diagram_ascii_art(TL, use_unicode=False)
         o o o o o o o o o o o o o o o
         | `-` `-` | `-` `-` | `-` | |
         |         `---------`     | |
         |                 .-------` |
         `---.             | .-------`
             |     .-----. | | .-----.
         .-. | .-. | .-. | | | | .-. |
         o o o o o o o o o o o o o o o
        sage: TL_diagram_ascii_art(TL, use_unicode=True)
         ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬
         │ ╰─╯ ╰─╯ │ ╰─╯ ╰─╯ │ ╰─╯ │ │
         │         ╰─────────╯     │ │
         │                 ╭───────╯ │
         ╰───╮             │ ╭───────╯
             │     ╭─────╮ │ │ ╭─────╮
         ╭─╮ │ ╭─╮ │ ╭─╮ │ │ │ │ ╭─╮ │
         ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬

        sage: TL = [(-20,-9), (-19,-10), (-18,-11), (-17,-16), (-15,-12), (2,3),
        ....:       (-14,-13), (-8,16), (-7,7), (-6,6), (-5,1), (-4,-3), (-2,-1),
        ....:       (4,5), (8,15), (9,10), (11,14), (12,13), (17,20), (18,19)]
        sage: TL_diagram_ascii_art(TL, use_unicode=False, blobs=[(-2,-1), (-5,1)])
         o o o o o o o o o o o o o o o o o o o o
         | `-` `-` | | | `-` | `-` | | | | `-` |
         |         | | |     `-----` | | `-----`
         |         | | `-------------` |
         `---0---. | | .---------------`
                 | | | | .---------------------.
                 | | | | | .-----------------. |
                 | | | | | | .-------------. | |
                 | | | | | | | .-----.     | | |
         .0. .-. | | | | | | | | .-. | .-. | | |
         o o o o o o o o o o o o o o o o o o o o
        sage: TL_diagram_ascii_art(TL, use_unicode=True, blobs=[(-2,-1), (-5,1)])
         ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬
         │ ╰─╯ ╰─╯ │ │ │ ╰─╯ │ ╰─╯ │ │ │ │ ╰─╯ │
         │         │ │ │     ╰─────╯ │ │ ╰─────╯
         │         │ │ ╰─────────────╯ │
         ╰───●───╮ │ │ ╭───────────────╯
                 │ │ │ │ ╭─────────────────────╮
                 │ │ │ │ │ ╭─────────────────╮ │
                 │ │ │ │ │ │ ╭─────────────╮ │ │
                 │ │ │ │ │ │ │ ╭─────╮     │ │ │
         ╭●╮ ╭─╮ │ │ │ │ │ │ │ │ ╭─╮ │ ╭─╮ │ │ │
         ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬ ⚬
    """
    def insert_pairing(cur, intervals):
        """
        Helper function to insert a possibly nested interval
        and push the others up, assuming inserting points from
        right-to-left.
        """
        for level in intervals:
            for j, I in enumerate(level):
                # Singleton intervals are vertical lines,
                #   so we don't need to worry about them
                if len(I) > 1 and I[0] < cur[0]:
                    cur, level[j] = level[j], cur
                    level.append([cur[0]])
                    level.append([cur[1]])
                    break
            else:
                level.append(cur)
                return  # We have stopped
        else:
            intervals.append([cur])
    # Build a set of intervals that defines where to draw the diagram
    intervals = [[]]
    propogating = []
    vertical = []
    top_intervals = [[]]
    num_left = 0
    num_right = 0

    def key_func(P):
        if P[1] < 0:  # cap
            return (0, P[0], P[1])
        elif P[0] > 0:  # cup
            return (3, -P[1], -P[0])
        else:
            bot, top = -P[0], P[1]
            if top < bot:  # left moving
                return (1, top, bot)
            elif top > bot:  # right moving
                return (2, -bot, -top)
            else:  # vertical line
                return (1, top, bot)
    diagram = sorted(diagram, key=key_func)
    # Since diagram is sorted in lex order, we will first do the matchings
    #   from right-to-left on the bottom, then the propogating lines, and
    #   then the matchings on the top from right-to-left.
    # Note that we need the top to go from right-to-left so the
    #   insert_pairing() function's assumptions are satisfied.
    for P in diagram:
        if P[1] < 0:  # Bottom matching
            insert_pairing([-P[1], -P[0], False, False], intervals)
        elif P[0] > 0:  # Top matching
            insert_pairing([P[0], P[1], True, True], top_intervals)
        else:  # Propogating line
            if -P[0] == P[1]:
                vertical.append(P[1])
            else:
                if -P[0] < P[1]:
                    num_right += 1
                else:
                    num_left += 1
                propogating.append(P)

    # Now piece together the intervals together
    total_prop = max(num_left, num_right)
    prop_intervals = [[] for _ in range(total_prop)]
    count_left = 0
    # Recall that the left-moving propogating lines come before the right-moving
    for i, P in enumerate(propogating):
        bot, top = P
        bot = -bot  # This makes it equal to its x-coordinate
        for level in intervals:
            level.append([bot])
        for level in top_intervals:
            level.append([top])
        left_moving = count_left < num_left
        if not left_moving:
            i -= num_left
        else:
            count_left += 1
        for j in range(i):
            prop_intervals[j].append([bot])
        for j in range(i+1, total_prop):
            prop_intervals[j].append([top])
        if not left_moving:
            top, bot = bot, top
        prop_intervals[i].append([top, bot, left_moving, not left_moving])
    intervals += prop_intervals
    intervals += reversed(top_intervals)
    for level in intervals:
        level.extend([i] for i in vertical)

    n = max(max(P) for P in diagram)

    # Finally, convert to a picture
    if use_unicode:
        from sage.typeset.unicode_art import UnicodeArt
        d = ["╭", "╮", "╰", "╯", "─", "│"]
        #db = ["┏", "┓", "┗", "┛", "━", "┃"]
        blob = '●'
        ret = [" ⚬" * n]
        char_art = UnicodeArt
    else:
        from sage.typeset.ascii_art import AsciiArt
        d = [".", ".", "`", "`", "-", "|"]
        # db = [".", ".", "`", "`", "=", "|"]
        blob = '0'
        ret = [" o" * n]
        char_art = AsciiArt

    def signed(val, pos):
        return val if pos else -val
    for level in reversed(intervals):
        cur = ""
        for I in sorted(level):
            cur += ' '*(2*I[0]-1 - len(cur))
            if len(I) == 1:
                cur += d[5] + ' '
            else:
                cur += d[2] if I[2] else d[0]
                if tuple(sorted([signed(I[0], I[2]), signed(I[1], I[3])])) in blobs:
                    cur += d[4] * (I[1]-I[0]-1)
                    cur += blob
                    cur += d[4] * (I[1]-I[0]-1)
                else:
                    cur += d[4] * (2*(I[1]-I[0])-1)
                cur += d[3] if I[3] else d[1]
        ret.append(cur)
    # Note that the top row and bottom row will be the same
    ret.append(ret[0])
    return char_art(ret, baseline=len(ret)//2)


def diagram_latex(diagram, fill=False, edge_options=None, edge_additions=None):
    r"""
    Return latex code for the diagram ``diagram`` using tikz.

    EXAMPLES::

        sage: from sage.combinat.diagram_algebras import PartitionDiagrams, diagram_latex
        sage: P = PartitionDiagrams(2)
        sage: D = P([[1,2],[-2,-1]])
        sage: print(diagram_latex(D)) # indirect doctest
        \begin{tikzpicture}[scale = 0.5,thick, baseline={(0,-1ex/2)}]
        \tikzstyle{vertex} = [shape = circle, minimum size = 7pt, inner sep = 1pt]
        \node[vertex] (G--2) at (1.5, -1) [shape = circle, draw] {};
        \node[vertex] (G--1) at (0.0, -1) [shape = circle, draw] {};
        \node[vertex] (G-1) at (0.0, 1) [shape = circle, draw] {};
        \node[vertex] (G-2) at (1.5, 1) [shape = circle, draw] {};
        \draw[] (G--2) .. controls +(-0.5, 0.5) and +(0.5, 0.5) .. (G--1);
        \draw[] (G-1) .. controls +(0.5, -0.5) and +(-0.5, -0.5) .. (G-2);
        \end{tikzpicture}
    """
    # these allow the view command to work (maybe move them
    # somewhere more appropriate?)
    from sage.misc.latex import latex
    latex.add_package_to_preamble_if_available('tikz')

    if fill:
        filled_str = ", fill"
    else:
        filled_str = ""

    if edge_options is None:
        edge_options = lambda P: ''
    if edge_additions is None:
        edge_additions = lambda P: ''

    def sgn(x):
        # Define the sign function
        if x > 0:
            return 1
        if x < 0:
            return -1
        return 0

    l1 = []  # list of blocks
    l2 = []  # list of nodes
    for i in list(diagram):
        l1.append(list(i))
        l2.extend(list(i))
    output = "\\begin{tikzpicture}[scale = 0.5,thick, baseline={(0,-1ex/2)}] \n\\tikzstyle{vertex} = [shape = circle, minimum size = 7pt, inner sep = 1pt] \n"  # setup beginning of picture
    for i in l2:  # add nodes
        output = output + "\\node[vertex] (G-{}) at ({}, {}) [shape = circle, draw{}] {{}}; \n".format(i, (abs(i)-1)*1.5, sgn(i), filled_str)
    for i in l1:  # add edges
        if len(i) > 1:
            l4 = list(i)
            posList = []
            negList = []
            for j in l4:  # sort list so rows are grouped together
                if j > 0:
                    posList.append(j)
                elif j < 0:
                    negList.append(j)
            posList.sort()
            negList.sort()
            l4 = posList + negList
            l5 = l4[:]  # deep copy
            for j in range(len(l5)):
                l5[j-1] = l4[j]  # create a permuted list
            if len(l4) == 2:
                l4.pop()
                l5.pop()  # pops to prevent duplicating edges
            for j in zip(l4, l5):
                xdiff = abs(j[1])-abs(j[0])
                y1 = sgn(j[0])
                y2 = sgn(j[1])
                if y2-y1 == 0 and abs(xdiff) < 5:  # if nodes are close to each other on same row
                    diffCo = (0.5+0.1*(abs(xdiff)-1))  # gets bigger as nodes are farther apart; max value of 1; min value of 0.5.
                    outVec = (sgn(xdiff)*diffCo, -1*diffCo*y1)
                    inVec = (-1*diffCo*sgn(xdiff), -1*diffCo*y2)
                elif y2-y1 != 0 and abs(xdiff) == 1:  # if nodes are close enough curviness looks bad.
                    outVec = (sgn(xdiff)*0.75, -1*y1)
                    inVec = (-1*sgn(xdiff)*0.75, -1*y2)
                else:
                    outVec = (sgn(xdiff)*1, -1*y1)
                    inVec = (-1*sgn(xdiff), -1*y2)
                output = output + "\\draw[{}] (G-{}) .. controls +{} and +{} .. {}(G-{}); \n".format(
                            edge_options(j), j[0], outVec, inVec, edge_additions(j), j[1])
    output = output + "\\end{tikzpicture}"  # end picture
    return output


class PottsRepresentation(CombinatorialFreeModule):
    r"""
    The Potts representation of the partition algebra.

    Let `P_n(d)` be the :class:`PartitionAlgebra` over `R` with the
    deformation parameter `d \in \ZZ_{>0}` being a positive integer.
    Recall the multiplication convention of diagrams in `P_n(d)`
    computing `D D'` by placing `D` above `D'`.

    The *Potts representation* is the  right `P_n(d)`-module on
    `M = V^{\otimes n}`, with `V = R^d`, with the action given as follows.
    We identify the natural basis vectors in `M` with words of length `n`
    in the alphabet `\{1, \dotsc, d\}` (which we call colors). For a basis
    vector `w` and diagram `D`, define `w \cdot D` as the sum over all `v`
    such that every part in `w D v` (consider this as coloring the nodes
    of `D`) is given by the same color.

    If `n` is a half integer, then there is an extra fixed color for the
    node `\lceil n \rceil`, which is called the *magnetic field direction*
    from the physics interpretation of this representation.

    EXAMPLES:

    In this example, we consider `R = \QQ` and use the Potts representation
    to construct the centralizer algebra of the left `S_{d-1}`-action on
    `V^{\otimes n}` with `V = \QQ^d` being the permutation action. ::

        sage: PA = algebras.Partition(5/2, QQ(2))
        sage: PR = PA.potts_representation(2)
        sage: mats = [PR.representation_matrix(x) for x in PA.basis()]
        sage: MS = mats[0].parent()
        sage: CM = MS.submodule(mats)
        sage: CM.dimension()
        16

    We check that this commutes with the `S_{d-1}`-action::

        sage: all((g * v) * x == g * (v * x) for g in PR.symmetric_group()
        ....:     for v in PR.basis() for x in PA.basis())
        True

    Next, we see that the centralizer of the `S_d`-action is smaller
    than the semisimple quotient of the partition algebra::

        sage: PA.dimension()
        52
        sage: len(PA.radical_basis())
        9
        sage: SQ = PA.semisimple_quotient()
        sage: SQ.dimension()
        43

    Next, we get orthogonal idempotents that project onto the central
    orthogonal idempotents in the semisimple quotient and construct
    the corresponding Peirce summands `e_i P_n(d) e_i`::

        sage: # long time
        sage: orth_idems = PA.orthogonal_idempotents_central_mod_radical()
        sage: algs = [PA.peirce_summand(idm, idm) for idm in orth_idems]
        sage: [A.dimension() for A in algs]
        [16, 2, 1, 25]

    We saw that we obtain the entire endomorphism algebra since `d = 2`
    and `S_{d-1}` is the trivial group. Hence, the 16 dimensional Peirce
    summand computed above is isomorphic to this endomorphism algebra
    (both are `4 \times 4` matrix algebras over `\QQ`). Hence, we have a
    natural quotient construction of the centralizer algebra from the
    partition algebra.

    Next, we consider a case with a nontrivial `S_d`-action (now it is `S_d`
    since the partition algebra has integer rank). We perform the same
    computations as before::

        sage: PA = algebras.Partition(2, QQ(2))
        sage: PA.dimension()
        15
        sage: PA.semisimple_quotient().dimension()
        10
        sage: orth_idems = PA.orthogonal_idempotents_central_mod_radical()
        sage: algs = [PA.peirce_summand(idm, idm) for idm in orth_idems]
        sage: [A.dimension() for A in algs]
        [4, 2, 4, 1]

        sage: PR = PA.potts_representation()
        sage: mats = [PR.representation_matrix(x) for x in PA.basis()]
        sage: MS = mats[0].parent()
        sage: cat = Algebras(QQ).WithBasis().Subobjects()
        sage: CM = MS.submodule(mats, category=cat)
        sage: CM.dimension()
        8

    To do the remainder of the computation, we need to monkey patch a
    ``product_on_basis`` method::

        sage: CM.product_on_basis
        NotImplemented
        sage: CM.product_on_basis = lambda x,y: CM.retract(CM.basis()[x].lift() * CM.basis()[y].lift())
        sage: CM.orthogonal_idempotents_central_mod_radical()
        (1/2*B[0] + 1/2*B[3] + 1/2*B[5] + 1/2*B[6],
         1/2*B[0] - 1/2*B[3] + 1/2*B[5] - 1/2*B[6])
        sage: CM.peirce_decomposition()
        [[Free module generated by {0, 1, 2, 3} over Rational Field,
          Free module generated by {} over Rational Field],
         [Free module generated by {} over Rational Field,
          Free module generated by {0, 1, 2, 3} over Rational Field]]

    Hence, we see that the centralizer algebra is isomorphic to a product
    of two `2 \times 2` matrix algebras (over `\QQ`), which are naturally
    a part of the partition algebra decomposition.

    Lastly, we verify the commuting actions::

        sage: all((g * v) * x == g * (v * x) for g in PR.symmetric_group()
        ....:     for v in PR.basis() for x in PA.basis())
        True

    REFERENCES:

    - [MR1998]_
    """
    def __init__(self, PA, y):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: PA = algebras.Partition(5/2, QQ(4))
            sage: PR = PA.potts_representation()
            sage: TestSuite(PR).run()

            sage: PA = algebras.Partition(5/2, QQ(2))  # not semisimple
            sage: PR = PA.potts_representation()
            sage: TestSuite(PR).run()
            sage: PR = PA.potts_representation(2)
            sage: TestSuite(PR).run()

        ::

            sage: PA = algebras.Partition(2, QQ(4))
            sage: PR = PA.potts_representation()
            sage: TestSuite(PR).run()

            sage: PA = algebras.Partition(3, QQ(2))  # not semisimple
            sage: PR = PA.potts_representation()
            sage: TestSuite(PR).run()

        TESTS::

            sage: PA = algebras.Partition(3, QQ(0))
            sage: PA.potts_representation()
            Traceback (most recent call last):
            ...
            ValueError: the partition algebra deformation parameter must be a positive integer

            sage: PA = algebras.Partition(3, QQ(2))
            sage: PA.potts_representation(1)
            Traceback (most recent call last):
            ...
            ValueError: the magnetic field direction should not be given for integer rank

            sage: PA = algebras.Partition(5/2, QQ(4))
            sage: PA.potts_representation(6)
            Traceback (most recent call last):
            ...
            ValueError: the magnetic field direction must be an integer in [1, 4]
            sage: PA.potts_representation(0)
            Traceback (most recent call last):
            ...
            ValueError: the magnetic field direction must be an integer in [1, 4]
            sage: PA.potts_representation(3/2)
            Traceback (most recent call last):
            ...
            ValueError: the magnetic field direction must be an integer in [1, 4]
        """
        if PA._q not in ZZ or PA._q <= 0:
            raise ValueError("the partition algebra deformation parameter must be a positive integer")
        self._d = ZZ(PA._q)
        self._PA = PA
        order = PA.order()
        if order not in ZZ:
            order = order.floor()
            if y not in ZZ or y < 1 or y > self._d:
                raise ValueError(f"the magnetic field direction must be an integer in [1, {self._d}]")
            y = ZZ(y)
        else:
            if y is not None:
                raise ValueError("the magnetic field direction should not be given for integer rank")
        self._y = y
        # _order is used to define the ordering of the basis elements
        self._num_factors = ZZ(order)
        from sage.combinat.words.words import Words
        indices = Words(self._d, self._num_factors)
        R = PA.base_ring()
        from sage.categories.modules_with_basis import ModulesWithBasis
        cat = ModulesWithBasis(R).FiniteDimensional()
        CombinatorialFreeModule.__init__(self, R, indices, prefix='P', category=cat)

        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        if self._y is None:
            self._symgp = SymmetricGroup(self._d)
        else:
            self._symgp = SymmetricGroup([i for i in range(1, self._d+1) if i != self._y])
        self._sga = self._symgp.algebra(R)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: PA = algebras.Partition(5/2, QQ(4))
            sage: PA.potts_representation(2)
            Potts representation with 4 colors and magnetic field direction 2 of
             Partition Algebra of rank 5/2 with parameter 4 over Rational Field

            sage: PA = algebras.Partition(3, QQ(4))
            sage: PA.potts_representation()
            Potts representation with 4 colors of
             Partition Algebra of rank 3 with parameter 4 over Rational Field
        """
        if self._y is not None:
            return "Potts representation with {} colors and magnetic field direction {} of {}".format(self._d, self._y, self._PA)
        return "Potts representation with {} colors of {}".format(self._d, self._PA)

    def _test_representation(self, **options):
        r"""
        Test that ``self`` is a representation.

        EXAMPLES::

            sage: PA = algebras.Partition(4, QQ(3))
            sage: PR = PA.potts_representation()
            sage: PR._test_representation()
        """
        tester = self._tester(**options)
        S = tester.some_elements()
        B = self._PA.basis()
        from sage.misc.misc import some_tuples
        for x, y in some_tuples(B, 2, tester._max_runs // len(S)):
            for v in S:
                tester.assertEqual((v * x) * y, v * (x * y))

    def _monomial(self, index):
        """
        Consturct a monomial of ``self`` from ``index``.

        EXAMPLES::

            sage: PA = algebras.Partition(4, QQ(3))
            sage: PR = PA.potts_representation()
            sage: PR._monomial([1, 2, 1, 3])
            P[word: 1213]
        """
        return super()._monomial(self._indices(index))

    def __getitem__(self, ind):
        r"""
        Return the basis element indexed by ``ind``.

        EXAMPLES::

            sage: PA = algebras.Partition(4, QQ(3))
            sage: PR = PA.potts_representation()
            sage: PR[1,2,1,3]
            P[word: 1213]

            sage: PA = algebras.Partition(1, QQ(4))
            sage: PR = PA.potts_representation()
            sage: PR[2]
            P[word: 2]
            sage: PR[10]
            Traceback (most recent call last):
            ...
            AttributeError: 'PottsRepresentation_with_category' object has no attribute 'list'
        """
        try:
            if self._num_factors == 1 and ind in ZZ:
                ind = self._indices([ind])
            else:
                ind = self._indices(ind)
        except (TypeError, ValueError):
            return super().__getitem__(ind)
        return self.monomial(ind)

    def partition_algebra(self):
        r"""
        Return the partition algebra that ``self`` is a representation of.

        EXAMPLES::

            sage: PA = algebras.Partition(3, QQ(4))
            sage: PR = PA.potts_representation()
            sage: PR.partition_algebra() is PA
            True
        """
        return self._PA

    def number_of_factors(self):
        r"""
        Return the number of factors defining ``self``.

        EXAMPLES::

            sage: PA = algebras.Partition(7/2, QQ(4))
            sage: PR = PA.potts_representation()
            sage: PR.number_of_factors()
            3
        """
        return self._num_factors

    def number_of_colors(self):
        r"""
        Return the number of colors defining ``self``.

        EXAMPLES::

            sage: PA = algebras.Partition(3, QQ(4))
            sage: PR = PA.potts_representation()
            sage: PR.number_of_colors()
            4
        """
        return self._d

    def magnetic_field_direction(self):
        r"""
        Return the magnetic field direction defining ``self``.

        EXAMPLES::

            sage: PA = algebras.Partition(7/2, QQ(4))
            sage: PR = PA.potts_representation(2)
            sage: PR.magnetic_field_direction()
            2
        """
        return self._y

    def symmetric_group(self):
        r"""
        Return the symmetric group that naturally acts on ``self``.

        EXAMPLES::

            sage: PA = algebras.Partition(3, QQ(4))
            sage: PR = PA.potts_representation()
            sage: PR.symmetric_group()
            Symmetric group of order 4! as a permutation group

            sage: PA = algebras.Partition(7/2, QQ(4))
            sage: PR = PA.potts_representation()
            sage: PR.symmetric_group().domain()
            {2, 3, 4}
            sage: PR = PA.potts_representation(2)
            sage: PR.symmetric_group().domain()
            {1, 3, 4}
            sage: PR = PA.potts_representation(4)
            sage: PR.symmetric_group().domain()
            {1, 2, 3}
        """
        return self._symgp

    def _basis_action(self, word, diagram):
        r"""
        Return the action of the partition algebra basis element indexed
        by ``diagram`` on the basis element of ``self`` indexed by ``word``.

        EXAMPLES::

            sage: PA = algebras.Partition(9/2, QQ(4))
            sage: PR = PA.potts_representation(2)
            sage: D = PA._indices([[1,2],[3,-3],[-1,-2,-4],[4,5,-5]])
            sage: PR._basis_action(PR._indices([1,1,3,4]), D)
            0
            sage: PR._basis_action(PR._indices([1,2,3,2]), D)
            0
            sage: PR._basis_action(PR._indices([1,1,3,2]), D)
            P[word: 1131] + P[word: 2232] + P[word: 3333] + P[word: 4434]
        """
        order = self._num_factors
        fixed = [None] * self._num_factors
        word = [None] + list(word) + [self._y]  # make this 1-based indexing
        neg_parts = []
        for part in diagram:
            if part[-1] < 0:
                # This should never happen since it should be a propagating block
                assert -part[0] != order + 1
                neg_parts.append(part)
                continue
            color = word[part[-1]]
            for i in reversed(part):
                if i > 0:
                    if word[i] != color:
                        return self.zero()
                else:  # i < 0
                    if -i == order + 1:
                        assert color == self._y
                    else:
                        fixed[-i-1] = color  # convert 1-based to 0-based
        if not neg_parts:
            return self._monomial(fixed)

        import itertools
        ret = []
        for c in itertools.product(range(1, self._d+1), repeat=len(neg_parts)):
            temp = list(fixed)  # make a copy
            for color, part in zip(c, neg_parts):
                for j in part:
                    temp[-j-1] = color  # convert 1-based to 0-based
            ret.append(self._indices(temp))
        return self.sum_of_monomials(ret)

    def _sym_group_action(self, g, word):
        r"""
        Return the action of the symmetric group element ``g``
        on the basis element of ``self`` indexed by ``word``.

        EXAMPLES::

            sage: PA = algebras.Partition(9/2, QQ(4))
            sage: PR = PA.potts_representation(2)
            sage: G = PR.symmetric_group()
            sage: g = G.an_element(); g
            (3,4)
            sage: PR._sym_group_action(g, PR._indices([1,2,3,4]))
            P[word: 1243]
        """
        return self._monomial([g(i) for i in word])

    def representation_matrix(self, elt):
        r"""
        Return the representation matrix of ``self`` in ``self``.

        EXAMPLES::

            sage: PA = algebras.Partition(7/2, QQ(2))
            sage: PR = PA.potts_representation()
            sage: PR.representation_matrix(PA.an_element())
            [7 0 3 0 2 0 0 0]
            [0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0]

            sage: all(b.to_vector() * PR.representation_matrix(x)  # long time
            ....:     == (b * x).to_vector()
            ....:     for b in PR.basis() for x in PA.basis())
            True

            sage: PA = algebras.Partition(2, QQ(2))
            sage: PR = PA.potts_representation()
            sage: [PR.representation_matrix(x) for x in PA.basis()]
            [
            [1 0 0 0]  [1 0 1 0]  [1 1 0 0]  [1 0 0 1]  [1 1 1 1]  [1 0 0 0]
            [0 0 0 0]  [0 0 0 0]  [0 0 0 0]  [0 0 0 0]  [0 0 0 0]  [1 0 0 0]
            [0 0 0 0]  [0 0 0 0]  [0 0 0 0]  [0 0 0 0]  [0 0 0 0]  [0 0 0 1]
            [0 0 0 1], [0 1 0 1], [0 0 1 1], [1 0 0 1], [1 1 1 1], [0 0 0 1],
            <BLANKLINE>
            [1 0 0 0]  [1 0 1 0]  [1 0 0 0]  [1 0 0 0]  [1 0 1 0]  [1 1 0 0]
            [0 0 1 0]  [1 0 1 0]  [0 1 0 0]  [0 0 0 1]  [0 1 0 1]  [1 1 0 0]
            [0 1 0 0]  [0 1 0 1]  [0 0 1 0]  [1 0 0 0]  [1 0 1 0]  [0 0 1 1]
            [0 0 0 1], [0 1 0 1], [0 0 0 1], [0 0 0 1], [0 1 0 1], [0 0 1 1],
            <BLANKLINE>
            [1 1 0 0]  [1 0 0 1]  [1 1 1 1]
            [0 0 1 1]  [1 0 0 1]  [1 1 1 1]
            [1 1 0 0]  [1 0 0 1]  [1 1 1 1]
            [0 0 1 1], [1 0 0 1], [1 1 1 1]
            ]

            sage: PA = algebras.Partition(5/2, QQ(2))
            sage: PR = PA.potts_representation()
            sage: all(PR.representation_matrix(x) * PR.representation_matrix(y)  # long time
            ....:     == PR.representation_matrix(x * y)
            ....:     for x in PA.basis() for y in PA.basis())
            True

            sage: PA = algebras.Partition(2, QQ(4))
            sage: PR = PA.potts_representation()
            sage: all(PR.representation_matrix(x) * PR.representation_matrix(y)
            ....:     == PR.representation_matrix(x * y)
            ....:     for x in PA.basis() for y in PA.basis())
            True
        """
        from sage.matrix.constructor import matrix
        return matrix([(b * elt).to_vector() for b in self.basis()])

    class Element(CombinatorialFreeModule.Element):
        def _acted_upon_(self, scalar, self_on_left=True):
            """
            Return the action of ``scalar`` on ``self``.

            EXAMPLES::

                sage: PA = algebras.Partition(5/2, QQ(4))
                sage: PR = PA.potts_representation()
                sage: v = PR.an_element(); v
                P[word: 11] + 3*P[word: 12] + 3*P[word: 13]
                sage: 3 * v
                3*P[word: 11] + 9*P[word: 12] + 9*P[word: 13]
                sage: v * 3
                3*P[word: 11] + 9*P[word: 12] + 9*P[word: 13]
                sage: x = PA.an_element(); x
                2*P{{-3, -2, -1, 1, 2, 3}} + 2*P{{-3, -2, 1, 2, 3}, {-1}}
                 + 3*P{{-3, -1, 1, 2, 3}, {-2}}
                sage: v * x
                7*P[word: 11] + 3*P[word: 12] + 3*P[word: 13] + 3*P[word: 14]
                 + 2*P[word: 21] + 2*P[word: 31] + 2*P[word: 41]
                sage: x * v
                Traceback (most recent call last):
                ...
                TypeError: unsupported operand parent(s) for *:
                 'Partition Algebra ...' and 'Potts representation ...'

                sage: PA = algebras.Partition(2, QQ(4))
                sage: PR = PA.potts_representation()
                sage: v = PR.an_element(); v
                P[word: 11] + 3*P[word: 12] + 3*P[word: 13]
                sage: G = PR.symmetric_group()
                sage: g = G.an_element(); g
                (2,3,4)
                sage: g * v
                P[word: 11] + 3*P[word: 13] + 3*P[word: 14]
                sage: SGA = G.algebra(ZZ)
                sage: a = SGA.an_element(); a
                () + (2,3,4) + 2*(1,3)(2,4) + 3*(1,4)(2,3)
                sage: a * v
                2*P[word: 11] + 3*P[word: 12] + 6*P[word: 13] + 3*P[word: 14]
                 + 6*P[word: 31] + 2*P[word: 33] + 6*P[word: 34]
                 + 9*P[word: 42] + 9*P[word: 43] + 3*P[word: 44]
            """
            # Check for scalars first
            ret = super()._acted_upon_(scalar, self_on_left)
            if ret is not None:
                return ret

            par = self.parent()

            if not self_on_left:
                # left action, so convert to the symmetric group algebra
                sP = scalar.parent()
                if sP is par._symgp:
                    scalar = par._sga.monomial(scalar)
                    sP = par._sga
                if sP is not par._sga:
                    try:
                        scalar = par._sga(scalar)
                    except (ValueError, TypeError):
                        return None
                return par.linear_combination((par._sym_group_action(g, wd), cg * cr)
                                              for (wd, cr) in self._monomial_coefficients.items()
                                              for (g, cg) in scalar.monomial_coefficients(copy=False).items())

            # right action, so convert to the partition algebra
            try:
                scalar = par._PA(scalar)
            except (ValueError, TypeError):
                return None
            return par.linear_combination((par._basis_action(wd, diag), cr * ca)
                                          for (wd, cr) in self._monomial_coefficients.items()
                                          for (diag, ca) in scalar.monomial_coefficients(copy=False).items())


#########################################################################
# START BORROWED CODE
#########################################################################
# Borrowed from Mike Hansen's original code -- global methods for dealing
# with partition diagrams, compositions of partition diagrams, and so on.
# --> CHANGED 'identity' to 'identity_set_partition' for enhanced clarity.
#########################################################################


def is_planar(sp):
    r"""
    Return ``True`` if the diagram corresponding to the set partition ``sp``
    is planar; otherwise, return ``False``.

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: da.is_planar( da.to_set_partition([[1,-2],[2,-1]]))
        False
        sage: da.is_planar( da.to_set_partition([[1,-1],[2,-2]]))
        True
    """
    # Singletons don't affect planarity
    to_consider = [x for x in map(list, sp) if len(x) > 1]
    n = len(to_consider)

    for i in range(n):
        # Get the positive and negative entries of this part
        ap = [x for x in to_consider[i] if x > 0]
        an = [abs(x) for x in to_consider[i] if x < 0]

        # Check if a includes numbers in both the top and bottom rows
        if ap and an:
            for j in range(n):
                if i == j:
                    continue
                # Get the positive and negative entries of this part
                bp = [x for x in to_consider[j] if x > 0]
                bn = [abs(x) for x in to_consider[j] if x < 0]

                # Skip the ones that don't involve numbers in both
                # the bottom and top rows
                if not bn or not bp:
                    continue

                # Make sure that if min(bp) > max(ap)
                # then min(bn) >  max(an)
                if max(bp) > max(ap):
                    if min(bn) < min(an):
                        return False

        # Go through the bottom and top rows
        for row in [ap, an]:
            if len(row) > 1:
                row.sort()
                for s in range(len(row)-1):
                    if row[s] + 1 == row[s+1]:
                        # No gap, continue on
                        continue

                    rng = list(range(row[s] + 1, row[s+1]))

                    # Go through and make sure any parts that
                    # contain numbers in this range are completely
                    # contained in this range
                    for j in range(n):
                        if i == j:
                            continue

                        # Make sure we make the numbers negative again
                        # if we are in the bottom row
                        if row is ap:
                            sr = set(rng)
                        else:
                            sr = set(-x for x in rng)

                        sj = set(to_consider[j])
                        intersection = sr.intersection(sj)
                        if intersection:
                            if sj != intersection:
                                return False

    return True


def to_graph(sp):
    r"""
    Return a graph representing the set partition ``sp``.

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: g = da.to_graph( da.to_set_partition([[1,-2],[2,-1]])); g
        Graph on 4 vertices

        sage: g.vertices(sort=True)
        [-2, -1, 1, 2]
        sage: g.edges(sort=True)
        [(-2, 1, None), (-1, 2, None)]
    """
    g = Graph()
    for part in sp:
        part_list = list(part)
        if len(part_list) > 0:
            g.add_vertex(part_list[0])
        for i in range(1, len(part_list)):
            g.add_vertex(part_list[i])
            g.add_edge(part_list[i-1], part_list[i])
    return g


def pair_to_graph(sp1, sp2):
    r"""
    Return a graph consisting of the disjoint union of the graphs of set
    partitions ``sp1`` and ``sp2`` along with edges joining the bottom
    row (negative numbers) of ``sp1`` to the top row (positive numbers)
    of ``sp2``.

    The vertices of the graph ``sp1`` appear in the result as pairs
    ``(k, 1)``, whereas the vertices of the graph ``sp2`` appear as
    pairs ``(k, 2)``.

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: sp1 = da.to_set_partition([[1,-2],[2,-1]])
        sage: sp2 = da.to_set_partition([[1,-2],[2,-1]])
        sage: g = da.pair_to_graph( sp1, sp2 ); g
        Graph on 8 vertices

        sage: g.vertices(sort=True)
        [(-2, 1), (-2, 2), (-1, 1), (-1, 2), (1, 1), (1, 2), (2, 1), (2, 2)]
        sage: g.edges(sort=True)
        [((-2, 1), (1, 1), None), ((-2, 1), (2, 2), None),
         ((-2, 2), (1, 2), None), ((-1, 1), (1, 2), None),
         ((-1, 1), (2, 1), None), ((-1, 2), (2, 2), None)]

    Another example which used to be wrong until :issue:`15958`::

        sage: sp3 = da.to_set_partition([[1, -1], [2], [-2]])
        sage: sp4 = da.to_set_partition([[1], [-1], [2], [-2]])
        sage: g = da.pair_to_graph( sp3, sp4 ); g
        Graph on 8 vertices

        sage: g.vertices(sort=True)
        [(-2, 1), (-2, 2), (-1, 1), (-1, 2), (1, 1), (1, 2), (2, 1), (2, 2)]
        sage: g.edges(sort=True)
        [((-2, 1), (2, 2), None), ((-1, 1), (1, 1), None),
         ((-1, 1), (1, 2), None)]
    """
    g = Graph()

    #Add the first set partition to the graph
    for part in sp1:
        part_list = list(part)
        if part_list:
            g.add_vertex( (part_list[0], 1) )

            #Add the edge to the second part of the graph
            if part_list[0] < 0:
                g.add_edge( (part_list[0], 1), (abs(part_list[0]), 2)  )

        for i in range(1, len(part_list)):
            g.add_vertex( (part_list[i], 1) )

            #Add the edge to the second part of the graph
            if part_list[i] < 0:
                g.add_edge( (part_list[i], 1), (abs(part_list[i]), 2) )

            #Add the edge between adjacent elements of a part
            g.add_edge( (part_list[i-1], 1), (part_list[i], 1) )

    #Add the second set partition to the graph
    for part in sp2:
        part_list = list(part)
        if part_list:
            g.add_vertex( (part_list[0], 2) )
        for i in range(1, len(part_list)):
            g.add_vertex( (part_list[i], 2) )
            g.add_edge( (part_list[i-1], 2), (part_list[i], 2) )

    return g


def propagating_number(sp):
    r"""
    Return the propagating number of the set partition ``sp``.

    The propagating number is the number of blocks with both a positive and
    negative number.

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: sp1 = da.to_set_partition([[1,-2],[2,-1]])
        sage: sp2 = da.to_set_partition([[1,2],[-2,-1]])
        sage: da.propagating_number(sp1)
        2
        sage: da.propagating_number(sp2)
        0
    """
    return sum(1 for part in sp if min(part) < 0 < max(part))


def to_set_partition(l, k=None):
    r"""
    Convert input to a set partition of `\{1, \ldots, k, -1, \ldots, -k\}`.

    Convert a list of a list of numbers to a set partitions. Each list
    of numbers in the outer list specifies the numbers contained in one
    of the blocks in the set partition.

    If `k` is specified, then the set partition will be a set partition
    of `\{1, \ldots, k, -1, \ldots, -k\}`. Otherwise, `k` will default to
    the minimum number needed to contain all of the specified numbers.

    INPUT:

    - ``l`` -- list of lists of integers
    - ``k`` -- integer (default: ``None``)

    OUTPUT: list of sets

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: f = lambda sp: SetPartition(da.to_set_partition(sp))
        sage: f([[1,-1],[2,-2]]) == SetPartition(da.identity_set_partition(2))
        True
        sage: da.to_set_partition([[1]])
        [{1}, {-1}]
        sage: da.to_set_partition([[1,-1],[-2,3]],9/2)
        [{-1, 1}, {-2, 3}, {2}, {-4, 4}, {-5, 5}, {-3}]
    """
    if k is None:
        if not l:
            return []
        k = max(max(map(abs, x)) for x in l)

    to_be_added = set(list(range(1, ceil(k+1))) + [-x for x in range(1, ceil(k+1))])

    sp = []
    for part in l:
        spart = set(part)
        to_be_added -= spart
        sp.append(spart)

    while to_be_added:
        i = to_be_added.pop()
        if -i in to_be_added:
            to_be_added.remove(-i)
            sp.append(set([i, -i]))
        else:
            sp.append(set([i]))

    return sp


def to_Brauer_partition(l, k=None):
    r"""
    Same as :func:`to_set_partition` but assumes omitted elements are
    connected straight through.

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: f = lambda sp: SetPartition(da.to_Brauer_partition(sp))
        sage: f([[1,2],[-1,-2]]) == SetPartition([[1,2],[-1,-2]])
        True
        sage: f([[1,3],[-1,-3]]) == SetPartition([[1,3],[-3,-1],[2,-2]])
        True
        sage: f([[1,-4],[-3,-1],[3,4]]) == SetPartition([[-3,-1],[2,-2],[1,-4],[3,4]])
        True
        sage: p = SetPartition([[1,2],[-1,-2],[3,-3],[4,-4]])
        sage: SetPartition(da.to_Brauer_partition([[1,2],[-1,-2]], k=4)) == p
        True
    """
    L = to_set_partition(l, k=k)
    paired = []
    not_paired = []
    L2 = (list(i) for i in L)
    for i in L2:
        if len(i) > 2:
            raise ValueError("blocks must have size at most 2, but {} has {}".format(i, len(i)))
        if len(i) == 2:
            paired.append(i)
        if len(i) == 1:
            not_paired.append(i)
    if any(i[0] in j or -1*i[0] in j for i in not_paired for j in paired):
        raise ValueError("unable to convert {} to a Brauer partition due to the invalid block {}".format(l, i))
    for i in not_paired:
        if [-i[0]] in not_paired:
            not_paired.remove([-i[0]])
        paired.append([i[0], -i[0]])
    return to_set_partition(paired)


def identity_set_partition(k):
    r"""
    Return the identity set partition `\{\{1, -1\}, \ldots, \{k, -k\}\}`.

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: SetPartition(da.identity_set_partition(2))
        {{-2, 2}, {-1, 1}}
    """
    if k in ZZ:
        return [[i,-i] for i in range(1, k + 1)]
    # Else k in 1/2 ZZ
    return [[i, -i] for i in range(1, k + ZZ(3)/ZZ(2))]

##########################################################################
# END BORROWED CODE
##########################################################################
