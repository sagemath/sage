r"""
Diagram and Partition Algebras

AUTHORS:

- Mike Hansen (2007): Initial version
- Stephen Doty, Aaron Lauve, George H. Seelinger (2012): Implementation of
  partition, Brauer, Temperley--Lieb, and ideal partition algebras
- Stephen Doty, Aaron Lauve, George H. Seelinger (2015): Implementation of
  ``*Diagram`` classes and other methods to improve diagram algebras.
"""

#*****************************************************************************
#  Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#                2012 Stephen Doty <doty@math.luc.edu>,
#                     Aaron Lauve <lauve@math.luc.edu>,
#                     George H. Seelinger <ghseeli@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#****************************************************************************
# python3
from __future__ import division
from six.moves import range

from sage.categories.algebras import Algebras
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.structure.element import generic_power
from sage.combinat.free_module import (CombinatorialFreeModule,
    CombinatorialFreeModuleElement)
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.combinat.combinat import bell_number, catalan_number
from sage.structure.global_options import GlobalOptions
from sage.combinat.set_partition import SetPartitions, SetPartition
from sage.combinat.partition import Partitions
from sage.combinat.symmetric_group_algebra import SymmetricGroupAlgebra_n
from sage.combinat.permutation import Permutations
from sage.sets.set import Set
from sage.graphs.graph import Graph
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.flatten import flatten
from sage.rings.all import ZZ


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
        [{{-2, -1, 1, 2}}, {{-2, -1, 2}, {1}}, {{-2, -1, 1}, {2}},
         {{-2}, {-1, 1, 2}}, {{-2, 1, 2}, {-1}}, {{-2, 1}, {-1, 2}},
         {{-2, 2}, {-1, 1}}, {{-2, -1}, {1, 2}}, {{-2, -1}, {1}, {2}},
         {{-2}, {-1, 2}, {1}}, {{-2, 2}, {-1}, {1}}, {{-2}, {-1, 1}, {2}},
         {{-2, 1}, {-1}, {2}}, {{-2}, {-1}, {1, 2}}, {{-2}, {-1}, {1}, {2}}]
        sage: [SetPartition(p) for p in da.partition_diagrams(3/2)]
        [{{-2, -1, 1, 2}}, {{-2, -1, 2}, {1}}, {{-2, 2}, {-1, 1}},
         {{-2, 1, 2}, {-1}}, {{-2, 2}, {-1}, {1}}]
    """
    if k in ZZ:
        S = SetPartitions(list(range(1, k+1)) + [-j for j in range(1, k+1)] )
        for p in Partitions(2*k):
            for i in S._iterator_part(p):
                yield i
    elif k + ZZ(1)/ZZ(2) in ZZ: # Else k in 1/2 ZZ
        k = ZZ(k + ZZ(1) / ZZ(2))
        S = SetPartitions(list(range(1, k+1)) + [-j for j in range(1, k)] )
        for p in Partitions(2*k-1):
            for sp in S._iterator_part(p):
                sp = list(sp)
                for i in range(len(sp)):
                    if k in sp[i]:
                        sp[i] += Set([-k])
                        break
                yield sp

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
        [{{-2, 1}, {-1, 2}}, {{-2, 2}, {-1, 1}}, {{-2, -1}, {1, 2}}]
        sage: [SetPartition(p) for p in da.brauer_diagrams(5/2)]
        [{{-3, 3}, {-2, 1}, {-1, 2}}, {{-3, 3}, {-2, 2}, {-1, 1}}, {{-3, 3}, {-2, -1}, {1, 2}}]
    """
    if k in ZZ:
        S = SetPartitions(list(range(1,k+1)) + [-j for j in range(1,k+1)],
                           [2 for j in range(1,k+1)] )
        for i in S._iterator_part(S.parts):
            yield list(i)
    elif k + ZZ(1) / ZZ(2) in ZZ: # Else k in 1/2 ZZ
        k = ZZ(k + ZZ(1) / ZZ(2))
        S = SetPartitions(list(range(1, k)) + [-j for j in range(1, k)],
                           [2 for j in range(1, k)] )
        for i in S._iterator_part(S.parts):
            yield list(i) + [[k, -k]]

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
        [{{-2, 2}, {-1, 1}}, {{-2, -1}, {1, 2}}]
        sage: [SetPartition(p) for p in da.temperley_lieb_diagrams(5/2)]
        [{{-3, 3}, {-2, 2}, {-1, 1}}, {{-3, 3}, {-2, -1}, {1, 2}}]
    """
    B = brauer_diagrams(k)
    for i in B:
        if is_planar(i):
            yield i

def planar_diagrams(k):
    r"""
    Return a generator of all planar diagrams of order ``k``.

    A planar diagram of order `k` is a partition diagram of order `k`
    that has no crossings.

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: [SetPartition(p) for p in da.planar_diagrams(2)]
        [{{-2, -1, 1, 2}}, {{-2, -1, 2}, {1}}, {{-2, -1, 1}, {2}},
         {{-2}, {-1, 1, 2}}, {{-2, 1, 2}, {-1}}, {{-2, 2}, {-1, 1}},
         {{-2, -1}, {1, 2}}, {{-2, -1}, {1}, {2}}, {{-2}, {-1, 2}, {1}},
         {{-2, 2}, {-1}, {1}}, {{-2}, {-1, 1}, {2}}, {{-2, 1}, {-1}, {2}},
         {{-2}, {-1}, {1, 2}}, {{-2}, {-1}, {1}, {2}}]
        sage: [SetPartition(p) for p in da.planar_diagrams(3/2)]
        [{{-2, -1, 1, 2}}, {{-2, -1, 2}, {1}}, {{-2, 2}, {-1, 1}},
         {{-2, 1, 2}, {-1}}, {{-2, 2}, {-1}, {1}}]
    """
    A = partition_diagrams(k)
    for i in A:
        if is_planar(i):
            yield i

def ideal_diagrams(k):
    r"""
    Return a generator of all "ideal" diagrams of order ``k``.

    An ideal diagram of order `k` is a partition diagram of order `k` with
    propagating number less than `k`.

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: [SetPartition(p) for p in da.ideal_diagrams(2)]
        [{{-2, -1, 1, 2}}, {{-2, -1, 2}, {1}}, {{-2, -1, 1}, {2}}, {{-2}, {-1, 1, 2}},
         {{-2, 1, 2}, {-1}}, {{-2, -1}, {1, 2}}, {{-2, -1}, {1}, {2}},
         {{-2}, {-1, 2}, {1}}, {{-2, 2}, {-1}, {1}}, {{-2}, {-1, 1}, {2}}, {{-2, 1},
         {-1}, {2}}, {{-2}, {-1}, {1, 2}}, {{-2}, {-1}, {1}, {2}}]
        sage: [SetPartition(p) for p in da.ideal_diagrams(3/2)]
        [{{-2, -1, 1, 2}}, {{-2, -1, 2}, {1}}, {{-2, 1, 2}, {-1}}, {{-2, 2}, {-1}, {1}}]
    """
    A = partition_diagrams(k)
    for i in A:
        if propagating_number(i) < k:
            yield i

class AbstractPartitionDiagram(SetPartition):
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
        sage: pd = da.AbstractPartitionDiagrams(da.partition_diagrams, 2)
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
        ValueError: this does not represent two rows of vertices
    """
    def __init__(self, parent, d):
        r"""
        Initialize ``self``.

            sage: import sage.combinat.diagram_algebras as da
            sage: pd = da.AbstractPartitionDiagrams(da.partition_diagrams, 2)
            sage: pd1 = da.AbstractPartitionDiagram(pd, ((-2,-1),(1,2)) )
        """
        self._base_diagram = tuple(sorted(tuple(sorted(i)) for i in d))
        super(AbstractPartitionDiagram, self).__init__(parent, self._base_diagram)

    # add options to class
    options=GlobalOptions('Brauer diagram', option_class='AbstractPartitionDiagram',
        module='sage.combinat.diagram_algebras',
        doc=r"""
        Set and display the global options for Brauer diagram (algebras). If no
        parameters are set, then the function returns a copy of the options
        dictionary.

        The ``options`` to diagram algebras can be accessed as the method
        :obj:`BrauerAlgebra.options` of :class:`BrauerAlgebra` and
        related classes.
        """,
        end_doc=r"""
        EXAMPLES::

            sage: R.<q> = QQ[]
            sage: BA = BrauerAlgebra(2, q)
            sage: E = BA([[1,2],[-1,-2]])
            sage: E
            B{{-2, -1}, {1, 2}}
            sage: BrauerAlgebra.options.display="compact"
            sage: E
            B[12/12;]
            sage: BrauerAlgebra.options._reset()
        """,
        display=dict(default="normal",
                       description='Specifies how the Brauer diagrams should be printed',
                       values=dict(normal="Using the normal representation",
                                   compact="Using the compact representation"),
                                   case_sensitive=False),
    )

    def check(self):
        r"""
        Check the validity of the input for the diagram.

        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: pd = da.AbstractPartitionDiagrams(da.partition_diagrams, 2)
            sage: pd1 = da.AbstractPartitionDiagram(pd, [[1,2],[-1,-2]]) # indirect doctest
            sage: pd2 = da.AbstractPartitionDiagram(pd, [[1,2],[3,4]]) # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: this does not represent two rows of vertices
        """
        if self._base_diagram:
            tst = sorted(flatten(self._base_diagram))
            if len(tst) % 2 or tst != list(range(-len(tst)//2,0)) + list(range(1,len(tst)//2+1)):
                raise ValueError("this does not represent two rows of vertices")

    def __eq__(self, other):
        r"""
        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: pd = da.AbstractPartitionDiagrams(da.partition_diagrams, 2)
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

    def __ne__(self, other):
        """
        Check not equals.

        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: pd = da.AbstractPartitionDiagrams(da.partition_diagrams, 2)
            sage: pd1 = da.AbstractPartitionDiagram(pd, [[1,2],[-1,-2]])
            sage: pd2 = da.AbstractPartitionDiagram(pd, [[1,-2],[-1,2]])
            sage: pd1 != pd2
            True
            sage: pd1 != ((-2,-1),(2,1))
            False
        """
        return not self == other

    def base_diagram(self):
        r"""
        Return the underlying implementation of the diagram.

        EXAMPLES::

            sage: import sage.combinat.diagram_algebras as da
            sage: pd = da.AbstractPartitionDiagrams(da.partition_diagrams, 2)
            sage: pd([[1,2],[-1,-2]]).base_diagram() == ((-2,-1),(1,2))
            True
        """
        return self._base_diagram # note, this works because self._base_diagram is immutable

    def diagram(self):
        r"""
        Return the underlying implementation of the diagram.

        EXAMPLES::

            sage: import sage.combinat.diagram_algebras as da
            sage: pd = da.AbstractPartitionDiagrams(da.partition_diagrams, 2)
            sage: pd([[1,2],[-1,-2]]).base_diagram() == pd([[1,2],[-1,-2]]).diagram()
            True
        """
        return self.base_diagram()

    def compose(self, other):
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
            sage: pd = da.AbstractPartitionDiagrams(da.partition_diagrams, 2)
            sage: pd([[1,2],[-1,-2]]).compose(pd([[1,2],[-1,-2]]))
            ({{-2, -1}, {1, 2}}, 1)
        """
        (composite_diagram, loops_removed) = set_partition_composition(self._base_diagram, other._base_diagram)
        return (self.__class__(self.parent(), composite_diagram), loops_removed)

    def propagating_number(self):
        r"""
        Return the propagating number of the diagram.

        The propagating number is the number of blocks with both a
        positive and negative number.

        EXAMPLES::

            sage: import sage.combinat.diagram_algebras as da
            sage: pd = da.AbstractPartitionDiagrams(da.partition_diagrams, 2)
            sage: d1 = pd([[1,-2],[2,-1]])
            sage: d1.propagating_number()
            2
            sage: d2 = pd([[1,2],[-2,-1]])
            sage: d2.propagating_number()
            0
        """
        return ZZ(sum(1 for part in self._base_diagram if min(part) < 0 and max(part) > 0))

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
        ValueError: all blocks must be of size 2
    """
    def __init__(self, parent, d):
        r"""
        Initialize ``self``.

            sage: import sage.combinat.diagram_algebras as da
            sage: bd = da.BrauerDiagrams(2)
            sage: bd1 = da.BrauerDiagram(bd, ((-2,-1),(1,2)) )
        """
        super(BrauerDiagram, self).__init__(parent,d)

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
            ValueError: all blocks must be of size 2
        """
        super(BrauerDiagram, self).check()
        if any(len(i) != 2 for i in self):
            raise ValueError("all blocks must be of size 2")

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

    def _repr_normal(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: import sage.combinat.diagram_algebras as da
            sage: bd = da.BrauerDiagrams(2)
            sage: bd([[1,2],[-1,-2]])._repr_normal()
            '{{-2, -1}, {1, 2}}'
        """
        return super(BrauerDiagram, self)._repr_()

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
        """
        (top, bot, thru) = self.involution_permutation_triple()
        bot.reverse()
        s1 = ".".join("".join(str(b) for b in block) for block in top)
        s2 = ".".join("".join(str(abs(k)) for k in sorted(block,reverse=True))
                              for block in bot)
        s3 = "".join(str(x) for x in thru)
        return "[{}/{};{}]".format(s1,s2,s3)

    def involution_permutation_triple(self, curt=True):
        r"""
        Return the involution permutation triple of ``self``.

        From Graham-Lehrer (see :class:`BrauerDiagrams`), a Brauer diagram
        is a triple `(D_1, D_2, \pi)`, where:

        - `D_1` is a partition of the top nodes;
        - `D_2` is a partition of the bottom nodes;
        - `\pi` is the induced permutation on the free nodes.

        INPUT:

        - ``curt`` -- (default: ``True``) if ``True``, then return bijection
          on free nodes as a one-line notation (standardized to look like a
          permutation), else, return the honest mapping, a list of pairs
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
        diagram = self.diagram()
        top = []
        bottom = []
        for v in diagram:
            if min(v)>0:
                top+=[v]
            if max(v)<0:
                bottom+=[v]
        if curt:
            perm = self.perm()
        else:
            perm = self.bijection_on_free_nodes()
        return (top,bottom,perm)

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
        terms = sorted(sorted(list(v), reverse=True) for v in self.diagram()
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

        .. TODO:: Come up with a better name?

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
    This is a class that generates partition diagrams.

    Thee primary use of this class is to serve as basis keys for
    diagram algebras, but diagrams also have properties in their
    own right. Furthermore, this class is meant to be extended to
    create more efficient contains methods.

    INPUT:

    - ``diagram_func`` -- generator; a function that can create the type
      of diagram desired
    - ``order`` -- integer or integer `+ 1/2`; the order of the diagrams

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: pd = da.AbstractPartitionDiagrams(da.partition_diagrams, 2)
        sage: pd
        Partition diagrams of order 2
        sage: [i for i in pd]
        [{{-2, -1, 1, 2}},
         {{-2, -1, 2}, {1}},
         {{-2, -1, 1}, {2}},
         {{-2}, {-1, 1, 2}},
         {{-2, 1, 2}, {-1}},
         {{-2, 1}, {-1, 2}},
         {{-2, 2}, {-1, 1}},
         {{-2, -1}, {1, 2}},
         {{-2, -1}, {1}, {2}},
         {{-2}, {-1, 2}, {1}},
         {{-2, 2}, {-1}, {1}},
         {{-2}, {-1, 1}, {2}},
         {{-2, 1}, {-1}, {2}},
         {{-2}, {-1}, {1, 2}},
         {{-2}, {-1}, {1}, {2}}]
        sage: pd.an_element() in pd
        True
        sage: elm = pd([[1,2],[-1,-2]])
        sage: elm in pd
        True
    """
    Element = AbstractPartitionDiagram

    def __init__(self, diagram_func, order, category=None):
        r"""
        See :class:`AbstractPartitionDiagrams` for full documentation.

        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: pd = da.AbstractPartitionDiagrams(da.partition_diagrams, 2)
            sage: TestSuite(pd).run() # long time
        """
        if category is None:
            category = FiniteEnumeratedSets()
        Parent.__init__(self, category=category)
        self.diagram_func = diagram_func
        self.order = order

    def __iter__(self):
        r"""
        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: pd = da.AbstractPartitionDiagrams(da.partition_diagrams, 2)
            sage: for i in pd: print(i) # indirect doctest
            {{-2, -1, 1, 2}}
            {{-2, -1, 2}, {1}}
            {{-2, -1, 1}, {2}}
            {{-2}, {-1, 1, 2}}
            {{-2, 1, 2}, {-1}}
            {{-2, 1}, {-1, 2}}
            {{-2, 2}, {-1, 1}}
            {{-2, -1}, {1, 2}}
            {{-2, -1}, {1}, {2}}
            {{-2}, {-1, 2}, {1}}
            {{-2, 2}, {-1}, {1}}
            {{-2}, {-1, 1}, {2}}
            {{-2, 1}, {-1}, {2}}
            {{-2}, {-1}, {1, 2}}
            {{-2}, {-1}, {1}, {2}}
        """
        for i in self.diagram_func(self.order):
            yield self.element_class(self, i)

    def _repr_(self):
        r"""
        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: da.AbstractPartitionDiagrams(da.partition_diagrams, 2)
            Partition diagrams of order 2
        """
        name = self.diagram_func.__name__.replace("_diagrams","").replace("_","").title()
        return "{} diagrams of order {}".format(name, self.order)

    def __contains__(self, obj):
        r"""
        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: pd = da.AbstractPartitionDiagrams(da.partition_diagrams, 2)
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
        if len(obj.base_diagram()) > 0:
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
            sage: pd = da.AbstractPartitionDiagrams(da.partition_diagrams, 2)
            sage: elm = pd([[1,2],[-1,-2]]); elm
            {{-2, -1}, {1, 2}}
            sage: pd([{1,2},{-1,-2}]) == elm
            True
            sage: pd( ((1,2),(-1,-2)) ) == elm
            True
            sage: pd( SetPartition([[1,2],[-1,-2]]) ) == elm
            True
        """
        return self.element_class(self, d)

class PartitionDiagrams(AbstractPartitionDiagrams):
    r"""
    This class represents all partition diagrams of integer or integer
    `+ 1/2` order.

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: pd = da.PartitionDiagrams(3)
        sage: pd.an_element() in pd
        True
        sage: pd.cardinality() == len(pd.list())
        True
    """
    def __init__(self, order, category=None):
        r"""
        Initialize ``self``.

        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: pd = da.PartitionDiagrams(2)
            sage: TestSuite(pd).run() # long time
        """
        super(PartitionDiagrams, self).__init__(partition_diagrams, order, category=category)

    def cardinality(self):
        r"""
        The cardinality of partition diagrams of integer order `n` is
        the `2n`-th Bell number.

        EXAMPLES::

            sage: import sage.combinat.diagram_algebras as da
            sage: pd = da.PartitionDiagrams(3)
            sage: pd.cardinality()
            203
        """
        if self.order in ZZ:
            return bell_number(2*self.order)
        return bell_number(2*(self.order-1/2))

class BrauerDiagrams(AbstractPartitionDiagrams):
    r"""
    This class represents all Brauer diagrams of integer or integer
    `+1/2` order. For more information on Brauer diagrams,
    see :class:`BrauerAlgebra`.

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: bd = da.BrauerDiagrams(3)
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
        [[/;321],
         [/;312],
         [23/12;1],
         [/;231],
         [/;132],
         [13/12;1],
         [/;213],
         [/;123],
         [12/12;1],
         [23/23;1],
         [13/23;1],
         [12/23;1],
         [23/13;1],
         [13/13;1],
         [12/13;1]]
        sage: bd.options._reset()
    """
    Element = BrauerDiagram
    options = AbstractPartitionDiagram.options

    def __init__(self, order, category=None):
        r"""
        Initialize ``self``.

        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: bd = da.BrauerDiagrams(2)
            sage: TestSuite(bd).run() # long time
        """
        super(BrauerDiagrams, self).__init__(brauer_diagrams, order, category=category)

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
        """
        return super(BrauerDiagrams, self).__contains__(obj) and [len(i) for i in obj] == [2]*self.order

    def _element_constructor_(self, d):
        r"""
        Construct an element of ``self``.

        EXAMPLES::

            sage: import sage.combinat.diagram_algebras as da
            sage: bd = da.BrauerDiagrams(2)
            sage: bd([[1,2],[-1,-2]])
            {{-2, -1}, {1, 2}}
        """
        return self.element_class(self, d)

    def cardinality(self):
        r"""
        Return the cardinality of ``self``.

        The number of Brauer diagrams of integer order `k` is `(2k-1)!!`.

        EXAMPLES::

            sage: import sage.combinat.diagram_algebras as da
            sage: bd = da.BrauerDiagrams(3)
            sage: bd.cardinality()
            15
        """
        if self.order in ZZ:
            return (2*self.order-1).multifactorial(2)
        else:
            return (2*(self.order-1/2)-1).multifactorial(2)

    def symmetric_diagrams(self,l=None,perm=None):
        r"""
        Return the list of brauer diagrams with symmetric placement of `l` arcs,
        and with free nodes permuted according to `perm`.

        EXAMPLES::

            sage: import sage.combinat.diagram_algebras as da
            sage: bd = da.BrauerDiagrams(4)
            sage: bd.symmetric_diagrams(l=1,perm=[2,1])
            [{{-4, -3}, {-2, 1}, {-1, 2}, {3, 4}},
             {{-4, -2}, {-3, 1}, {-1, 3}, {2, 4}},
             {{-4, 1}, {-3, -2}, {-1, 4}, {2, 3}},
             {{-4, -1}, {-3, 2}, {-2, 3}, {1, 4}},
             {{-4, 2}, {-3, -1}, {-2, 4}, {1, 3}},
             {{-4, 3}, {-3, 4}, {-2, -1}, {1, 2}}]
        """
        # perm = permutation on free nodes
        # l = number of arcs
        n = self.order
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
        Construct a Bruaer diagram of ``self`` from an involution
        permutation triple.

        A Brauer diagram can be represented as a triple where the first
        entry is a list of arcs on the top row of the diagram, the second
        entry is a list of arcs on the bottom row of the diagram, and the
        third entry is a permutation on the remaining nodes. This triple
        is called the *involution permutation triple*. For more
        information, see [GL1996]_.

        INPUT:

        - ``D1_D2_pi``-- a list or tuple where the first entry is a list of
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
        """
        try:
            (D1,D2,pi) = tuple(D1_D2_pi)
        except ValueError:
            raise ValueError("argument %s not in correct form; must be a tuple (D1, D2, pi)" % D1_D2_pi)
        D1 = [[abs(x) for x in b] for b in D1 if len(b) == 2] # not needed if argument correctly passed at outset.
        D2 = [[abs(x) for x in b] for b in D2 if len(b) == 2] # ditto.
        nD2 = [map(lambda i: -i,b) for b in D2]
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
        sage: td = da.TemperleyLiebDiagrams(3)
        sage: td.an_element() in td
        True
        sage: td.cardinality() == len(td.list())
        True
    """
    def __init__(self, order):
        r"""
        See :class:`TemperleyLiebDiagrams` for full documentation.

        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: td = da.TemperleyLiebDiagrams(2)
            sage: TestSuite(td).run() # long time
        """
        super(TemperleyLiebDiagrams, self).__init__(temperley_lieb_diagrams, order)

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
            return catalan_number(self.order)
        else:
            return catalan_number(self.order-1/2)

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
            obj = self._element_constructor_(obj)
        if obj not in BrauerDiagrams(self.order):
            return False
        if not is_planar(obj):
            return False
        return True

class PlanarDiagrams(AbstractPartitionDiagrams):
    r"""
    All planar diagrams of integer or integer `+1/2` order.

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: pld = da.PlanarDiagrams(3)
        sage: pld.an_element() in pld
        True
        sage: pld.cardinality() == len(pld.list())
        True
    """
    def __init__(self, order):
        r"""
        See :class:`PlanarDiagrams` for full documentation.

        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: pld = da.PlanarDiagrams(2)
            sage: TestSuite(pld).run() # long time
        """
        super(PlanarDiagrams, self).__init__(planar_diagrams, order)

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
        if self.order in ZZ:
            return catalan_number(2*self.order)
        else:
            return catalan_number(2*self.order-1)

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
            obj = self._element_constructor_(obj)
        return super(PlanarDiagrams, self).__contains__(obj) and is_planar(obj)

class IdealDiagrams(AbstractPartitionDiagrams):
    r"""
    All "ideal" diagrams of integer or integer `+1/2` order.

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: id = da.IdealDiagrams(3)
        sage: id.an_element() in id
        True
        sage: id.cardinality() == len(id.list())
        True
    """
    def __init__(self, order):
        r"""
        See :class:`TemperleyLiebDiagrams` for full documentation.

        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: id = da.IdealDiagrams(2)
            sage: TestSuite(id).run() # long time
        """
        super(IdealDiagrams, self).__init__(ideal_diagrams, order)

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
            obj = self._element_constructor_(obj)
        return super(IdealDiagrams, self).__contains__(obj) and obj.propagating_number() < self.order

class DiagramAlgebra(CombinatorialFreeModule):
    r"""
    Abstract class for diagram algebras and is not designed to be used
    directly. If used directly, the class could create an "algebra"
    that is not actually an algebra.

    TESTS::

        sage: import sage.combinat.diagram_algebras as da
        sage: R.<x> = QQ[]
        sage: D = da.DiagramAlgebra(2, x, R, 'P', da.PartitionDiagrams(2))
        sage: sorted(D.basis())
        [P{{-2}, {-1}, {1}, {2}},
         P{{-2}, {-1}, {1, 2}},
         P{{-2}, {-1, 1}, {2}},
         P{{-2}, {-1, 1, 2}},
         P{{-2}, {-1, 2}, {1}},
         P{{-2, -1}, {1}, {2}},
         P{{-2, -1}, {1, 2}},
         P{{-2, -1, 1}, {2}},
         P{{-2, -1, 1, 2}},
         P{{-2, -1, 2}, {1}},
         P{{-2, 1}, {-1}, {2}},
         P{{-2, 1}, {-1, 2}},
         P{{-2, 1, 2}, {-1}},
         P{{-2, 2}, {-1}, {1}},
         P{{-2, 2}, {-1, 1}}]
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
            sage: D = da.DiagramAlgebra(2, x, R, 'P', da.PartitionDiagrams(2))
            sage: TestSuite(D).run()
        """
        self._prefix = prefix
        self._q = base_ring(q)
        self._k = k
        self._base_diagrams = diagrams
        category = Algebras(base_ring.category()).FiniteDimensional().WithBasis()
        category = category.or_subcategory(category)
        CombinatorialFreeModule.__init__(self, base_ring, diagrams,
                    category=category, prefix=prefix, bracket=False)

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
            ValueError: {{-2, 1}, {-1, 2}} is not an index of a basis element
        """
        if self.basis().keys().is_parent_of(set_partition):
            return self.basis()[set_partition]
        if isinstance(set_partition, SymmetricGroupAlgebra_n.Element):
            return self._apply_module_morphism(set_partition, self._perm_to_Blst, self)
        sp = self._base_diagrams(set_partition) # attempt conversion
        if sp in self.basis().keys():
            return self.basis()[sp]

        raise ValueError("invalid input of {0}".format(set_partition))

    def __getitem__(self, i):
        """
        Get the basis item of ``self`` indexed by ``i``.

        EXAMPLES::

            sage: import sage.combinat.diagram_algebras as da
            sage: R.<x> = QQ[]
            sage: D = da.DiagramAlgebra(2, x, R, 'P', da.PartitionDiagrams(2))
            sage: sp = da.PartitionDiagrams(2)( [[1,2], [-1,-2]] )
            sage: D[sp]
            P{{-2, -1}, {1, 2}}
        """
        i = self._base_diagrams(i)
        if i in self.basis().keys():
            return self.basis()[i]
        raise ValueError("{0} is not an index of a basis element".format(i))

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
            ValueError: {{-2, 1}, {-1, 2}} is not an index of a basis element
        """
        ## 'perm' is a permutation in one-line notation
        ## turns w into an expression suitable for the element constructor.
        u = sorted(w)
        p = [[u[i],-x] for i,x in enumerate(w)]
        return self[p]

    def order(self):
        r"""
        Return the order of ``self``.

        The order of a partition algebra is defined as half of the number
        of nodes in the diagrams.

        EXAMPLES::

            sage: q = var('q')
            sage: PA = PartitionAlgebra(2, q)
            sage: PA.order()
            2
        """
        return self._k

    def set_partitions(self):
        r"""
        Return the collection of underlying set partitions indexing the
        basis elements of a given diagram algebra.

        .. TODO:: Is this really necessary?

        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: R.<x> = QQ[]
            sage: D = da.DiagramAlgebra(2, x, R, 'P', da.PartitionDiagrams(2))
            sage: list(D.set_partitions()) == list(da.PartitionDiagrams(2))
            True
        """
        return self.basis().keys()

    def product_on_basis(self, d1, d2):
        r"""
        Return the product `D_{d_1} D_{d_2}` by two basis diagrams.

        TESTS::

            sage: import sage.combinat.diagram_algebras as da
            sage: R.<x> = QQ[]
            sage: D = da.DiagramAlgebra(2, x, R, 'P', da.PartitionDiagrams(2))
            sage: sp = da.PartitionDiagrams(2)([[1,2],[-1,-2]])
            sage: D.product_on_basis(sp, sp)
            x*P{{-2, -1}, {1, 2}}
        """
        if not self._indices.is_parent_of(d1):
            d1 = self._indices(d1)
        if not self._indices.is_parent_of(d2):
            d2 = self._indices(d2)
        (composite_diagram, loops_removed) = d1.compose(d2)
        return self.term(composite_diagram, self._q**loops_removed)

    @cached_method
    def one_basis(self):
        r"""
        The following constructs the identity element of ``self``.

        It is not called directly; instead one should use ``DA.one()`` if
        ``DA`` is a defined diagram algebra.

        EXAMPLES::

            sage: import sage.combinat.diagram_algebras as da
            sage: R.<x> = QQ[]
            sage: D = da.DiagramAlgebra(2, x, R, 'P', da.PartitionDiagrams(2))
            sage: D.one_basis()
            {{-2, 2}, {-1, 1}}
        """
        return self._base_diagrams(identity_set_partition(self._k))

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
            \draw (G--2) .. controls +(-0.5, 0.5) and +(0.5, 0.5) .. (G--1);
            \draw (G-1) .. controls +(0.5, -0.5) and +(-0.5, -0.5) .. (G-2);
            \end{tikzpicture}
        """
        # these allow the view command to work (maybe move them somewhere more appropriate?)
        from sage.misc.latex import latex
        latex.add_to_mathjax_avoid_list('tikzpicture')
        latex.add_package_to_preamble_if_available('tikz')
        # Define the sign function
        def sgn(x):
            if x > 0:
                return 1
            if x < 0:
                return -1
            return 0
        l1 = [] #list of blocks
        l2 = [] #lsit of nodes
        for i in list(diagram):
            l1.append(list(i))
            for j in list(i):
                l2.append(j)
        output = "\\begin{tikzpicture}[scale = 0.5,thick, baseline={(0,-1ex/2)}] \n\\tikzstyle{vertex} = [shape = circle, minimum size = 7pt, inner sep = 1pt] \n" #setup beginning of picture
        for i in l2: #add nodes
            output = output + "\\node[vertex] (G-{}) at ({}, {}) [shape = circle, draw] {{}}; \n".format(i, (abs(i)-1)*1.5, sgn(i))
        for i in l1: #add edges
            if len(i) > 1:
                l4 = list(i)
                posList = []
                negList = []
                for i in l4: #sort list so rows are grouped together
                    if i > 0:
                        posList.append(i)
                    elif i < 0:
                        negList.append(i)
                posList.sort()
                negList.sort()
                l4 = posList + negList
                l5 = l4[:] #deep copy
                for j in range(len(l5)):
                    l5[j-1] = l4[j] #create a permuted list
                if len(l4) == 2:
                    l4.pop()
                    l5.pop() #pops to prevent duplicating edges
                for j in zip(l4, l5):
                    xdiff = abs(j[1])-abs(j[0])
                    y1 = sgn(j[0])
                    y2 = sgn(j[1])
                    if y2-y1 == 0 and abs(xdiff) < 5: #if nodes are close to each other on same row
                        diffCo = (0.5+0.1*(abs(xdiff)-1)) #gets bigger as nodes are farther apart; max value of 1; min value of 0.5.
                        outVec = (sgn(xdiff)*diffCo, -1*diffCo*y1)
                        inVec = (-1*diffCo*sgn(xdiff), -1*diffCo*y2)
                    elif y2-y1 != 0 and abs(xdiff) == 1: #if nodes are close enough curviness looks bad.
                        outVec = (sgn(xdiff)*0.75, -1*y1)
                        inVec = (-1*sgn(xdiff)*0.75, -1*y2)
                    else:
                        outVec = (sgn(xdiff)*1, -1*y1)
                        inVec = (-1*sgn(xdiff), -1*y2)
                    output = output + "\\draw (G-{}) .. controls +{} and +{} .. (G-{}); \n".format(j[0], outVec, inVec, j[1])
        output = output + "\\end{tikzpicture} \n" #end picture
        return output

    # The following subclass provides a few additional methods for
    # partition algebra elements.
    class Element(CombinatorialFreeModuleElement):
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
                sage: elt.diagrams()
                [{{-2}, {-1}, {1, 2}}, {{-2, -1}, {1, 2}}]
            """
            return self.support()

class PartitionAlgebra(DiagramAlgebra):
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

    The Iwahori--Hecke algebra of type `A` (with a single parameter) is
    naturally a subalgebra of the partition algebra.

    The partition algebra is regarded as an example of a "diagram algebra"
    due to the fact that its natural basis is given by certain graphs
    often called diagrams.

    An excellent reference for partition algebras and their various
    subalgebras (Brauer algebra, Temperley--Lieb algebra, etc) is the
    paper [HR2005]_.

    INPUT:

    - ``k`` -- rank of the algebra

    - ``q`` -- the deformation parameter `q`

    OPTIONAL ARGUMENTS:

    - ``base_ring`` -- (default ``None``) a ring containing ``q``; if
      ``None``, then Sage automatically chooses the parent of ``q``

    - ``prefix`` -- (default ``"P"``) a label for the basis elements

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
    over `\ZZ`::

        sage: R.<x> = ZZ[]
        sage: P = PartitionAlgebra(2, x, R)
        sage: P
        Partition Algebra of rank 2 with parameter x
         over Univariate Polynomial Ring in x over Integer Ring
        sage: P.basis().list()
        [P{{-2, -1, 1, 2}}, P{{-2, -1, 2}, {1}},
         P{{-2, -1, 1}, {2}}, P{{-2}, {-1, 1, 2}},
         P{{-2, 1, 2}, {-1}}, P{{-2, 1}, {-1, 2}},
         P{{-2, 2}, {-1, 1}}, P{{-2, -1}, {1, 2}},
         P{{-2, -1}, {1}, {2}}, P{{-2}, {-1, 2}, {1}},
         P{{-2, 2}, {-1}, {1}}, P{{-2}, {-1, 1}, {2}},
         P{{-2, 1}, {-1}, {2}}, P{{-2}, {-1}, {1, 2}},
         P{{-2}, {-1}, {1}, {2}}]
        sage: E = P([[1,2],[-2,-1]]); E
        P{{-2, -1}, {1, 2}}
        sage: E in P.basis().list()
        True
        sage: E^2
        x*P{{-2, -1}, {1, 2}}
        sage: E^5
        x^4*P{{-2, -1}, {1, 2}}
        sage: (P([[2,-2],[-1,1]]) - 2*P([[1,2],[-1,-2]]))^2
        (4*x-4)*P{{-2, -1}, {1, 2}} + P{{-2, 2}, {-1, 1}}

    One can work with partition algebras using a symbol for the parameter,
    leaving the base ring unspecified. This implies that the underlying
    base ring is Sage's symbolic ring.

    ::

        sage: q = var('q')
        sage: PA = PartitionAlgebra(2, q); PA
        Partition Algebra of rank 2 with parameter q over Symbolic Ring
        sage: PA([[1,2],[-2,-1]])^2 == q*PA([[1,2],[-2,-1]])
        True
        sage: (PA([[2, -2], [1, -1]]) - 2*PA([[-2, -1], [1, 2]]))^2 == (4*q-4)*PA([[1, 2], [-2, -1]]) + PA([[2, -2], [1, -1]])
        True

    The identity element of the partition algebra is the set
    partition `\{\{1,-1\}, \{2,-2\}, \ldots, \{k,-k\}\}`::

        sage: P = PA.basis().list()
        sage: PA.one()
        P{{-2, 2}, {-1, 1}}
        sage: PA.one()*P[7] == P[7]
        True
        sage: P[7]*PA.one() == P[7]
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
        sage: (PA([[2, -2], [1, -1]]) - 2*PA([[-2, -1], [1, 2]]))^2 == 16*PA([[-2, -1], [1, 2]]) + PA([[2, -2], [1, -1]])
        True

    TESTS:

    A computation that returned an incorrect result until :trac:`15958`::

        sage: A = PartitionAlgebra(1,17)
        sage: g = SetPartitionsAk(1).list()
        sage: a = A[g[1]]
        sage: a
        P{{-1}, {1}}
        sage: a*a
        17*P{{-1}, {1}}

    Symmetric group algebra elements can also be coerced into the
    partition algebra::

        sage: S = SymmetricGroupAlgebra(SR, 2)
        sage: A = PartitionAlgebra(2, x, SR)
        sage: S([2,1])*A([[1,-1],[2,-2]])
        P{{-2, 1}, {-1, 2}}

    REFERENCES:

    .. [HR2005] Tom Halverson and Arun Ram, *Partition algebras*. European
       Journal of Combinatorics **26** (2005), 869--921.
    """
    @staticmethod
    def __classcall_private__(cls, k, q, base_ring=None, prefix="P"):
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
        return super(PartitionAlgebra, cls).__classcall__(cls, k, q, base_ring, prefix)

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

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: S = SymmetricGroupAlgebra(R, 4)
            sage: A = PartitionAlgebra(4, x, R)
            sage: A._coerce_map_from_(S)
            Generic morphism:
              From: Symmetric group algebra of order 4 over Univariate Polynomial Ring in x over Rational Field
              To:   Partition Algebra of rank 4 with parameter x over Univariate Polynomial Ring in x over Rational Field
            sage: Sp = SymmetricGroupAlgebra(QQ, 4)
            sage: A._coerce_map_from_(Sp)
            Generic morphism:
              From: Symmetric group algebra of order 4 over Rational Field
              To:   Partition Algebra of rank 4 with parameter x over Univariate Polynomial Ring in x over Rational Field
        """
        if isinstance(R, SymmetricGroupAlgebra_n):
            if R.n == self._k and self.base_ring().has_coerce_map_from(R.base_ring()):
                return R.module_morphism(self._perm_to_Blst, codomain=self)
            return None
        return super(PartitionAlgebra, self)._coerce_map_from_(R)

class SubPartitionAlgebra(DiagramAlgebra):
    """
    A subalgebra of the partition algebra indexed by a subset of the diagrams.
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
        DiagramAlgebra.__init__(self, k, q, base_ring, prefix, diagrams, category)

    #These methods allow for a subalgebra to be correctly identified in a partition algebra
    def ambient(self):
        r"""
        Return the partition algebra ``self`` is a sub-algebra of.

        EXAMPLES::

            sage: x = var('x')
            sage: BA = BrauerAlgebra(2, x)
            sage: BA.ambient()
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
        phi = self.module_morphism(lambda d: amb.monomial(d),
                                   codomain=amb, category=self.category())
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

class BrauerAlgebra(SubPartitionAlgebra):
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

    OPTIONAL ARGUMENTS:

    - ``base_ring`` -- (default ``None``) a ring containing ``q``; if ``None``
      then just takes the parent of ``q``

    - ``prefix`` -- (default ``"B"``) a label for the basis elements

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
        sage: b = B.basis().list()
        sage: b
        [B{{-2, 1}, {-1, 2}}, B{{-2, 2}, {-1, 1}}, B{{-2, -1}, {1, 2}}]
        sage: b[2]
        B{{-2, -1}, {1, 2}}
        sage: b[2]^2
        x*B{{-2, -1}, {1, 2}}
        sage: b[2]^5
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
    def __classcall_private__(cls, k, q, base_ring=None, prefix="B"):
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
        return super(BrauerAlgebra, cls).__classcall__(cls, k, q, base_ring, prefix)

    def __init__(self, k, q, base_ring, prefix):
        r"""
        Initialize ``self``.

        TESTS::

            sage: R.<q> = QQ[]
            sage: BA = BrauerAlgebra(2, q, R)
            sage: TestSuite(BA).run()
        """
        SubPartitionAlgebra.__init__(self, k, q, base_ring, prefix, BrauerDiagrams(k))

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
        """
        if isinstance(R, SymmetricGroupAlgebra_n):
            if R.n == self._k and self.base_ring().has_coerce_map_from(R.base_ring()):
                return R.module_morphism(self._perm_to_Blst, codomain=self)
            return None
        return super(BrauerAlgebra, self)._coerce_map_from_(R)

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
        set_partition = to_Brauer_partition(set_partition, k = self.order())
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
        I = lambda x: self._indices(to_Brauer_partition(x, k=k))
        R = self.base_ring()
        one = R.one()
        d = {self.one_basis(): R( (self._q-1) / 2 )}
        for i in range(1,j):
            d[I([[i,-j],[j,-i]])] = one
            d[I([[i,j],[-i,-j]])] = -one
        return self._from_dict(d, remove_zeros=True)

class TemperleyLiebAlgebra(SubPartitionAlgebra):
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

    OPTIONAL ARGUMENTS:

    - ``base_ring`` -- (default ``None``) a ring containing ``q``; if ``None``
      then just takes the parent of ``q``

    - ``prefix`` -- (default ``"T"``) a label for the basis elements

    EXAMPLES:

    We define the Temperley--Lieb algebra of rank `2` with parameter
    `x` over `\ZZ`::

        sage: R.<x> = ZZ[]
        sage: T = TemperleyLiebAlgebra(2, x, R); T
        Temperley-Lieb Algebra of rank 2 with parameter x
         over Univariate Polynomial Ring in x over Integer Ring
        sage: T.basis()
        Lazy family (Term map from Temperleylieb diagrams of order 2
         to Temperley-Lieb Algebra of rank 2 with parameter x
         over Univariate Polynomial Ring in x over
         Integer Ring(i))_{i in Temperleylieb diagrams of order 2}
        sage: b = T.basis().list()
        sage: b
        [T{{-2, 2}, {-1, 1}}, T{{-2, -1}, {1, 2}}]
        sage: b[1]
        T{{-2, -1}, {1, 2}}
        sage: b[1]^2 == x*b[1]
        True
        sage: b[1]^5 == x^4*b[1]
        True
    """
    @staticmethod
    def __classcall_private__(cls, k, q, base_ring=None, prefix="T"):
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
        return super(TemperleyLiebAlgebra, cls).__classcall__(cls, k, q, base_ring, prefix)

    def __init__(self, k, q, base_ring, prefix):
        r"""
        Initialize ``self``

        TESTS::

            sage: R.<q> = QQ[]
            sage: TL = TemperleyLiebAlgebra(2, q, R)
            sage: TestSuite(TL).run()
        """
        SubPartitionAlgebra.__init__(self, k, q, base_ring, prefix, TemperleyLiebDiagrams(k))

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
            ValueError: {{-2, 1}, {-1, 2}} is not an index of a basis element
        """
        if isinstance(set_partition, SymmetricGroupAlgebra_n.Element):
            return SubPartitionAlgebra._element_constructor_(self, set_partition)
        set_partition = to_Brauer_partition(set_partition, k = self.order())
        return SubPartitionAlgebra._element_constructor_(self, set_partition)

class PlanarAlgebra(SubPartitionAlgebra):
    """
    A planar algebra.

    The planar algebra of rank `k` is an algebra with basis indexed by the
    collection of all planar set partitions of
    `\{1, \ldots, k, -1, \ldots, -k\}`.

    This algebra is thus a subalgebra of the partition algebra. For more
    information, see :class:`PartitionAlgebra`.

    INPUT:

    - ``k`` -- rank of the algebra

    - ``q`` -- the deformation parameter `q`

    OPTIONAL ARGUMENTS:

    - ``base_ring`` -- (default ``None``) a ring containing ``q``; if ``None``
      then just takes the parent of ``q``

    - ``prefix`` -- (default ``"Pl"``) a label for the basis elements

    EXAMPLES:

    We define the planar algebra of rank `2` with parameter
    `x` over `\ZZ`::

        sage: R.<x> = ZZ[]
        sage: Pl = PlanarAlgebra(2, x, R); Pl
        Planar Algebra of rank 2 with parameter x over Univariate Polynomial Ring in x over Integer Ring
        sage: Pl.basis().list()
        [Pl{{-2, -1, 1, 2}}, Pl{{-2, -1, 2}, {1}},
         Pl{{-2, -1, 1}, {2}}, Pl{{-2}, {-1, 1, 2}},
         Pl{{-2, 1, 2}, {-1}}, Pl{{-2, 2}, {-1, 1}},
         Pl{{-2, -1}, {1, 2}}, Pl{{-2, -1}, {1}, {2}},
         Pl{{-2}, {-1, 2}, {1}}, Pl{{-2, 2}, {-1}, {1}},
         Pl{{-2}, {-1, 1}, {2}}, Pl{{-2, 1}, {-1}, {2}},
         Pl{{-2}, {-1}, {1, 2}}, Pl{{-2}, {-1}, {1}, {2}}]
        sage: E = Pl([[1,2],[-1,-2]])
        sage: E^2 == x*E
        True
        sage: E^5 == x^4*E
        True
    """
    @staticmethod
    def __classcall_private__(cls, k, q, base_ring=None, prefix="Pl"):
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
        return super(PlanarAlgebra, cls).__classcall__(cls, k, q, base_ring, prefix)

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
        return "Planar Algebra of rank {} with parameter {} over {}".format(self._k,
                self._q, self.base_ring())

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
        sage: I.basis().list()
        [I{{-2, -1, 1, 2}}, I{{-2, -1, 2}, {1}},
         I{{-2, -1, 1}, {2}}, I{{-2}, {-1, 1, 2}},
         I{{-2, 1, 2}, {-1}}, I{{-2, -1}, {1, 2}},
         I{{-2, -1}, {1}, {2}}, I{{-2}, {-1, 2}, {1}},
         I{{-2, 2}, {-1}, {1}}, I{{-2}, {-1, 1}, {2}},
         I{{-2, 1}, {-1}, {2}}, I{{-2}, {-1}, {1, 2}},
         I{{-2}, {-1}, {1}, {2}}]
        sage: E = I([[1,2],[-1,-2]])
        sage: E^2 == x*E
        True
        sage: E^5 == x^4*E
        True
    """
    @staticmethod
    def __classcall_private__(cls, k, q, base_ring=None, prefix="I"):
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
        return super(PropagatingIdeal, cls).__classcall__(cls, k, q, base_ring, prefix)

    def __init__(self, k, q, base_ring, prefix):
        r"""
        Initialize ``self``.

        TESTS::

            sage: R.<q> = QQ[]
            sage: I = PropagatingIdeal(2, q, R)
            sage: TestSuite(I).run() # Not tested -- needs non-unital algebras category
        """
        # This should be the category of non-unital fin-dim algebras with basis
        category = Algebras(base_ring.category()).FiniteDimensional().WithBasis()
        SubPartitionAlgebra.__init__(self, k, q, base_ring, prefix,
                                     IdealDiagrams(k), category)

    @cached_method
    def one_basis(self):
        r"""
        The propagating ideal is a non-unital algebra, i.e. it does not have a
        multiplicative identity.

        EXAMPLES::

            sage: R.<q> = QQ[]
            sage: I = PropagatingIdeal(2, q, R)
            sage: I.one_basis()
            Traceback (most recent call last):
            ...
            ValueError: The ideal partition algebra is not unital
            sage: I.one()
            Traceback (most recent call last):
            ...
            ValueError: The ideal partition algebra is not unital
        """
        raise ValueError("The ideal partition algebra is not unital")
        #return identity_set_partition(self._k)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: PropagatingIdeal(2, x, R)
            Propagating Ideal of rank 2 with parameter x over Univariate Polynomial Ring in x over Rational Field
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

            - ``n`` -- a positive integer

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
    #Singletons don't affect planarity
    to_consider = [x for x in (list(_) for _ in sp) if len(x) > 1]
    n = len(to_consider)

    for i in range(n):
        #Get the positive and negative entries of this part
        ap = [x for x in to_consider[i] if x>0]
        an = [abs(x) for x in to_consider[i] if x<0]

        #Check if a includes numbers in both the top and bottom rows
        if len(ap) > 0 and len(an) > 0:
            for j in range(n):
                if i == j:
                    continue
                #Get the positive and negative entries of this part
                bp = [x for x in to_consider[j] if x>0]
                bn = [abs(x) for x in to_consider[j] if x<0]

                #Skip the ones that don't involve numbers in both
                #the bottom and top rows
                if not bn or not bp:
                    continue

                #Make sure that if min(bp) > max(ap)
                #then min(bn) >  max(an)
                if max(bp) > max(ap):
                    if min(bn) < min(an):
                        return False

        #Go through the bottom and top rows
        for row in [ap, an]:
            if len(row) > 1:
                row.sort()
                for s in range(len(row)-1):
                    if row[s] + 1 == row[s+1]:
                        #No gap, continue on
                        continue

                    rng = list(range(row[s] + 1, row[s+1]))

                    #Go through and make sure any parts that
                    #contain numbers in this range are completely
                    #contained in this range
                    for j in range(n):
                        if i == j:
                            continue

                        #Make sure we make the numbers negative again
                        #if we are in the bottom row
                        if row is ap:
                            sr = set(rng)
                        else:
                            sr = set((-1*x for x in rng))

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

        sage: g.vertices()
        [-2, -1, 1, 2]
        sage: g.edges()
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

        sage: g.vertices()
        [(-2, 1), (-2, 2), (-1, 1), (-1, 2), (1, 1), (1, 2), (2, 1), (2, 2)]
        sage: g.edges()
        [((-2, 1), (1, 1), None), ((-2, 1), (2, 2), None),
         ((-2, 2), (1, 2), None), ((-1, 1), (1, 2), None),
         ((-1, 1), (2, 1), None), ((-1, 2), (2, 2), None)]

    Another example which used to be wrong until :trac:`15958`::

        sage: sp3 = da.to_set_partition([[1, -1], [2], [-2]])
        sage: sp4 = da.to_set_partition([[1], [-1], [2], [-2]])
        sage: g = da.pair_to_graph( sp3, sp4 ); g
        Graph on 8 vertices

        sage: g.vertices()
        [(-2, 1), (-2, 2), (-1, 1), (-1, 2), (1, 1), (1, 2), (2, 1), (2, 2)]
        sage: g.edges()
        [((-2, 1), (2, 2), None), ((-1, 1), (1, 1), None),
         ((-1, 1), (1, 2), None)]
    """
    g = Graph()

    #Add the first set partition to the graph
    for part in sp1:
        part_list = list(part)
        if len(part_list) > 0:
            g.add_vertex( (part_list[0],1) )

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
        if len(part_list) > 0:
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
    pn = 0
    for part in sp:
        if min(part) < 0  and max(part) > 0:
            pn += 1
    return pn

def to_set_partition(l, k=None):
    r"""
    Convert a list of a list of numbers to a set partitions. Each list
    of numbers in the outer list specifies the numbers contained in one
    of the blocks in the set partition.

    If `k` is specified, then the set partition will be a set partition
    of `\{1, \ldots, k, -1, \ldots, -k\}`. Otherwise, `k` will default to
    the minimum number needed to contain all of the specified numbers.

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: f = lambda sp: SetPartition(da.to_set_partition(sp))
        sage: f([[1,-1],[2,-2]]) == SetPartition(da.identity_set_partition(2))
        True
    """
    if k is None:
        if l == []:
            return []
        else:
            k = max( (max( map(abs, x) ) for x in l) )

    to_be_added = set( list(range(1, k+1)) + [-1*x for x in range(1, k+1)] )

    sp = []
    for part in l:
        spart = set(part)
        to_be_added -= spart
        sp.append(spart)

    for singleton in to_be_added:
        sp.append(set([singleton]))

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
    L2 = []
    paired = []
    not_paired = []
    for i in L:
        L2.append(list(i))
    for i in L2:
        if len(i) >= 3:
            raise ValueError("blocks must have size at most 2, but {0} has {1}".format(i, len(i)))
        if len(i) == 2:
            paired.append(i)
        if len(i) == 1:
            not_paired.append(i)
    if any(i[0] in j or -1*i[0] in j for i in not_paired for j in paired):
        raise ValueError("unable to convert {0} to a Brauer partition due to the invalid block {1}".format(l, i))
    for i in not_paired:
        if [-1*i[0]] in not_paired:
            not_paired.remove([-1*i[0]])
        paired.append([i[0], -1*i[0]])
    return to_set_partition(paired)

def identity_set_partition(k):
    """
    Return the identity set partition `\{\{1, -1\}, \ldots, \{k, -k\}\}`

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: SetPartition(da.identity_set_partition(2))
        {{-2, 2}, {-1, 1}}
    """
    if k in ZZ:
        return [[i,-i] for i in range(1, k + 1)]
    # Else k in 1/2 ZZ
    return [[i, -i] for i in range(1, k + ZZ(3)/ZZ(2))]

def set_partition_composition(sp1, sp2):
    r"""
    Return a tuple consisting of the composition of the set partitions
    ``sp1`` and ``sp2`` and the number of components removed from the middle
    rows of the graph.

    EXAMPLES::

        sage: import sage.combinat.diagram_algebras as da
        sage: sp1 = da.to_set_partition([[1,-2],[2,-1]])
        sage: sp2 = da.to_set_partition([[1,-2],[2,-1]])
        sage: p, c = da.set_partition_composition(sp1, sp2)
        sage: (SetPartition(p), c) == (SetPartition(da.identity_set_partition(2)), 0)
        True
    """
    g = pair_to_graph(sp1, sp2)
    connected_components = g.connected_components()

    res = []
    total_removed = 0
    for cc in connected_components:
        #Remove the vertices that live in the middle two rows
        new_cc = [x for x in cc if not ( (x[0]<0 and x[1] == 1) or (x[0]>0 and x[1]==2))]

        if new_cc == []:
            if len(cc) > 1:
                total_removed += 1
        else:
            res.append( set((x[0] for x in new_cc)) )

    return (res, total_removed)

##########################################################################
# END BORROWED CODE
##########################################################################

# Deprecations from trac:18555. July 2016
from sage.misc.superseded import deprecated_function_alias
AbstractPartitionDiagram.global_options=deprecated_function_alias(18555, AbstractPartitionDiagram.options)
BrauerDiagramOptions = deprecated_function_alias(18555, AbstractPartitionDiagram.options)
BrauerDiagrams.global_options = deprecated_function_alias(18555, BrauerDiagrams.options)
