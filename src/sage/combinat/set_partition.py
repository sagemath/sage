r"""
Set partitions

This module defines a class for immutable partitioning of a set. For
mutable version see :func:`DisjointSet`.

AUTHORS:

- Mike Hansen
- MuPAD-Combinat developers (for algorithms and design inspiration).
- Travis Scrimshaw (2013-02-28): Removed ``CombinatorialClass`` and added
  entry point through :class:`SetPartition`.
- Martin Rubey (2017-10-10): Cleanup, add crossings and nestings, add
  random generation.
"""
# ****************************************************************************
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
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import itertools
from itertools import repeat
from sage.sets.set import Set, Set_generic

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.list_clone import ClonableArray
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.lazy_import import lazy_import
from sage.rings.infinity import infinity
from sage.rings.integer import Integer
from sage.combinat.combinatorial_map import combinatorial_map
from sage.combinat.set_partition_iterator import (set_partition_iterator,
                                                  set_partition_iterator_blocks)
from sage.combinat.partition import Partition, Partitions
from sage.combinat.combinat import bell_number, stirling_number2 as stirling2
from sage.combinat.permutation import Permutation
from sage.arith.misc import factorial
from sage.misc.prandom import random, randint, sample
from sage.sets.disjoint_set import DisjointSet

lazy_import('sage.combinat.posets.hasse_diagram', 'HasseDiagram')
lazy_import('sage.probability.probability_distribution', 'GeneralDiscreteDistribution')


class AbstractSetPartition(ClonableArray,
        metaclass=InheritComparisonClasscallMetaclass):
    r"""
    Methods of set partitions which are independent of the base set
    """

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: S = SetPartitions(4)
            sage: S([[1,3],[2,4]])
            {{1, 3}, {2, 4}}
        """
        try:
            s = [sorted(x) for x in self]
        except TypeError:
            s = [sorted(x, key=str) for x in self]
        return '{' + ', '.join('{' + repr(x)[1:-1] + '}' for x in s) + '}'

    def __hash__(self):
        """
        Return the hash of ``self``.

        The parent is not included as part of the hash.

        EXAMPLES::

            sage: P = SetPartitions(4)
            sage: A = SetPartition([[1], [2,3], [4]])
            sage: B = P([[1], [2,3], [4]])
            sage: hash(A) == hash(B)
            True
        """
        return sum(hash(x) for x in self)

    def __eq__(self, y):
        """
        Check equality of ``self`` and ``y``.

        The parent is not included as part of the equality check.

        EXAMPLES::

            sage: P = SetPartitions(4)
            sage: A = SetPartition([[1], [2,3], [4]])
            sage: B = P([[1], [2,3], [4]])
            sage: A == B
            True
            sage: C = P([[2, 3], [1], [4]])
            sage: A == C
            True
            sage: D = P([[1], [2, 4], [3]])
            sage: A == D
            False

        Note that this may give incorrect answers if the base set is not totally ordered::

            sage: a,b = frozenset([0,1]), frozenset([2,3])
            sage: p1 = SetPartition([[a], [b]])
            sage: p2 = SetPartition([[b], [a]])
            sage: p1 == p2
            False
        """
        if not isinstance(y, AbstractSetPartition):
            return False
        return list(self) == list(y)

    def __ne__(self, y):
        """
        Check lack of equality of ``self`` and ``y``.

        The parent is not included as part of the equality check.

        EXAMPLES::

            sage: P = SetPartitions(4)
            sage: A = SetPartition([[1], [2,3], [4]])
            sage: B = P([[1], [2,3], [4]])
            sage: A != B
            False
            sage: C = P([[2, 3], [1], [4]])
            sage: A != C
            False
            sage: D = P([[1], [2, 4], [3]])
            sage: A != D
            True

        Note that this may give incorrect answers if the base set is not totally ordered::

            sage: a,b = frozenset([0,1]), frozenset([2,3])
            sage: p1 = SetPartition([[a], [b]])
            sage: p2 = SetPartition([[b], [a]])
            sage: p1 != p2
            True
        """
        return not (self == y)

    def __lt__(self, y):
        """
        Check that ``self`` is less than ``y``.

        The ordering used is lexicographic, where:

        - a set partition is considered as the list of its parts
          sorted by increasing smallest element;

        - each part is regarded as a list of its elements, sorted
          in increasing order;

        - the parts themselves are compared lexicographically.

        EXAMPLES::

            sage: P = SetPartitions(4)
            sage: A = P([[1], [2,3], [4]])
            sage: B = SetPartition([[1,2,3], [4]])
            sage: A < B
            True
            sage: C = P([[1,2,4], [3]])
            sage: B < C
            True
            sage: B < B
            False
            sage: D = P([[1,4], [2], [3]])
            sage: E = P([[1,4], [2,3]])
            sage: D < E
            True
            sage: F = P([[1,2,4], [3]])
            sage: E < C
            False
            sage: A < E
            True
            sage: A < C
            True
        """
        if not isinstance(y, AbstractSetPartition):
            return False
        return [sorted(i) for i in self] < [sorted(i) for i in y]

    def __gt__(self, y):
        """
        Check that ``self`` is greater than ``y``.

        The ordering used is lexicographic, where:

        - a set partition is considered as the list of its parts
          sorted by increasing smallest element;

        - each part is regarded as a list of its elements, sorted
          in increasing order;

        - the parts themselves are compared lexicographically.

        EXAMPLES::

            sage: P = SetPartitions(4)
            sage: A = P([[1], [2,3], [4]])
            sage: B = SetPartition([[1,2,3], [4]])
            sage: B > A
            True
            sage: A > B
            False
        """
        if not isinstance(y, AbstractSetPartition):
            return False
        return [sorted(i) for i in self] > [sorted(i) for i in y]

    def __le__(self, y):
        """
        Check that ``self`` is less than or equals ``y``.

        The ordering used is lexicographic, where:

        - a set partition is considered as the list of its parts
          sorted by increasing smallest element;

        - each part is regarded as a list of its elements, sorted
          in increasing order;

        - the parts themselves are compared lexicographically.

        EXAMPLES::

            sage: P = SetPartitions(4)
            sage: A = P([[1], [2,3], [4]])
            sage: B = SetPartition([[1,2,3], [4]])
            sage: A <= B
            True
            sage: A <= A
            True
        """
        return self == y or self < y

    def __ge__(self, y):
        """
        Check that ``self`` is greater than or equals ``y``.

        The ordering used is lexicographic, where:

        - a set partition is considered as the list of its parts
          sorted by increasing smallest element;

        - each part is regarded as a list of its elements, sorted
          in increasing order;

        - the parts themselves are compared lexicographically.

        EXAMPLES::

            sage: P = SetPartitions(4)
            sage: A = P([[1], [2,3], [4]])
            sage: B = SetPartition([[1,2,3], [4]])
            sage: B >= A
            True
            sage: B >= B
            True
        """
        return self == y or self > y

    def __mul__(self, other):
        r"""
        The product of the set partitions ``self`` and ``other``.

        The product of two set partitions `B` and `C` is defined as the
        set partition whose parts are the nonempty intersections between
        each part of `B` and each part of `C`. This product is also
        the infimum of `B` and `C` in the classical set partition
        lattice (that is, the coarsest set partition which is finer than
        each of `B` and `C`). Consequently, ``inf`` acts as an alias for
        this method.

        .. SEEALSO::

            :meth:`sup`

        EXAMPLES::

            sage: x = SetPartition([ [1,2], [3,5,4] ])
            sage: y = SetPartition(( (3,1,2), (5,4) ))
            sage: x * y
            {{1, 2}, {3}, {4, 5}}

            sage: S = SetPartitions(4)
            sage: sp1 = S([[2,3,4], [1]])
            sage: sp2 = S([[1,3], [2,4]])
            sage: s = S([[2,4], [3], [1]])
            sage: sp1.inf(sp2) == s
            True

        TESTS:

        Here is a different implementation of the ``__mul__`` method
        (one that was formerly used for the ``inf`` method, before it
        was realized that the methods do the same thing)::

            sage: def mul2(s, t):
            ....:     temp = [ss.intersection(ts) for ss in s for ts in t]
            ....:     temp = filter(bool, temp)
            ....:     return s.__class__(s.parent(), temp)

        Let us check that this gives the same as ``__mul__`` on set
        partitions of `\{1, 2, 3, 4\}`::

            sage: all( all( mul2(s, t) == s * t for s in SetPartitions(4) )
            ....:      for t in SetPartitions(4) )
            True
        """
        new_composition = []
        for B in self:
            for C in other:
                BintC = B.intersection(C)
                if BintC:
                    new_composition.append(BintC)
        return SetPartition(new_composition)

    inf = __mul__

    def sup(self, t):
        """
        Return the supremum of ``self`` and ``t`` in the classical set
        partition lattice.

        The supremum of two set partitions `B` and `C` is obtained as the
        transitive closure of the relation which relates `i` to `j` if
        and only if `i` and `j` are in the same part in at least
        one of the set partitions `B` and `C`.

        .. SEEALSO::

            :meth:`__mul__`

        EXAMPLES::

            sage: S = SetPartitions(4)
            sage: sp1 = S([[2,3,4], [1]])
            sage: sp2 = S([[1,3], [2,4]])
            sage: s = S([[1,2,3,4]])
            sage: sp1.sup(sp2) == s
            True
        """
        res = list(self)
        for p in t:
            # find blocks in res which intersect p
            inters = [(i, q) for i, q in enumerate(res)
                      if any(a in q for a in p)]
            # remove these blocks from res
            for i, _ in reversed(inters):
                del res[i]
            # add the union
            res.append([e for _, q in inters for e in q])
        return self.parent()(res)

    def standard_form(self):
        r"""
        Return ``self`` as a list of lists.

        When the ground set is totally ordered, the elements of each
        block are listed in increasing order.

        This is not related to standard set partitions (which simply
        means set partitions of `[n] = \{ 1, 2, \ldots , n \}` for some
        integer `n`) or standardization (:meth:`standardization`).

        EXAMPLES::

            sage: [x.standard_form() for x in SetPartitions(4, [2,2])]                  # needs sage.graphs sage.rings.finite_rings
            [[[1, 2], [3, 4]], [[1, 4], [2, 3]], [[1, 3], [2, 4]]]

        TESTS::

            sage: SetPartition([(1, 9, 8), (2, 3, 4, 5, 6, 7)]).standard_form()
            [[1, 8, 9], [2, 3, 4, 5, 6, 7]]
        """
        return [sorted(i) for i in self]

    def base_set(self):
        """
        Return the base set of ``self``, which is the union of all parts
        of ``self``.

        EXAMPLES::

            sage: SetPartition([[1], [2,3], [4]]).base_set()
            {1, 2, 3, 4}
            sage: SetPartition([[1,2,3,4]]).base_set()
            {1, 2, 3, 4}
            sage: SetPartition([]).base_set()
            {}
        """
        return Set(e for p in self for e in p)

    def base_set_cardinality(self):
        """
        Return the cardinality of the base set of ``self``, which is the sum
        of the sizes of the parts of ``self``.

        This is also known as the *size* (sometimes the *weight*) of
        a set partition.

        EXAMPLES::

            sage: SetPartition([[1], [2,3], [4]]).base_set_cardinality()
            4
            sage: SetPartition([[1,2,3,4]]).base_set_cardinality()
            4
        """
        return sum(len(x) for x in self)

    def coarsenings(self):
        """
        Return a list of coarsenings of ``self``.

        .. SEEALSO::

            :meth:`refinements`

        EXAMPLES::

            sage: SetPartition([[1,3],[2,4]]).coarsenings()
            [{{1, 2, 3, 4}}, {{1, 3}, {2, 4}}]
            sage: SetPartition([[1],[2,4],[3]]).coarsenings()
            [{{1, 2, 3, 4}},
             {{1, 2, 4}, {3}},
             {{1, 3}, {2, 4}},
             {{1}, {2, 3, 4}},
             {{1}, {2, 4}, {3}}]
            sage: SetPartition([]).coarsenings()
            [{}]
        """
        SP = SetPartitions(len(self))

        def union(s):
            # Return the partition obtained by combining, for every
            # part of s, those parts of self which are indexed by
            # the elements of this part of s into a single part.
            ret = []
            for part in s:
                cur = []
                for i in part:
                    cur.extend(self[i - 1])  # -1 for indexing
                ret.append(cur)
            return ret
        return [self.parent()(union(s)) for s in SP]

    def max_block_size(self):
        r"""
        The maximum block size of the diagram.

        EXAMPLES::

            sage: # needs sage.modules
            sage: from sage.combinat.diagram_algebras import PartitionDiagram, PartitionDiagrams
            sage: pd = PartitionDiagram([[1,-3,-5],[2,4],[3,-1,-2],[5],[-4]])
            sage: pd.max_block_size()
            3
            sage: sorted(d.max_block_size() for d in PartitionDiagrams(2))
            [1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4]
            sage: sorted(sp.max_block_size() for sp in SetPartitions(3))
            [1, 2, 2, 2, 3]
        """
        return max(len(block) for block in self)

    def conjugate(self):
        r"""
        An involution exchanging singletons and circular adjacencies.

        This method implements the definition of the conjugate of
        a set partition defined in [Cal2005]_.

        INPUT:

        - ``self`` -- set partition of an ordered set

        OUTPUT: a set partition

        EXAMPLES::

            sage: SetPartition([[1,6,7],[2,8],[3,4,5]]).conjugate()
            {{1, 4, 7}, {2, 8}, {3}, {5}, {6}}
            sage: all(sp.conjugate().conjugate()==sp for sp in SetPartitions([1,3,5,7]))
            True
            sage: SetPartition([]).conjugate()
            {}
        """
        def next_one(a, support):
            return support[(support.index(a) + 1) % len(support)]

        def addback(S, terminals, rsupport):
            out = list(S)
            for a in terminals * 2:
                if a not in out and next_one(a, rsupport) in out:
                    out.append(a)
            return out

        def pre_conjugate(sp):
            if len(sp) <= 1:
                return SetPartition([[a] for S in sp for a in S])
            if sp.max_block_size() == 1:
                return SetPartition([sp.base_set()])
            support = sorted(a for S in sp for a in S)
            initials = [a for S in sp for a in S if next_one(a, support) in S]
            singletons = [a for S in sp for a in S if len(S) == 1]
            if not initials and not singletons:
                return sp
            rho = pre_conjugate(
                SetPartition([[a for a in S if a not in initials]
                for S in sp if len(S) > 1 and any(a not in initials for a in S)]))
            # add back initials as singletons and singletons as terminals
            return SetPartition([addback(S, singletons, support[::-1])
                for S in rho] + [[a] for a in initials])
        support = sorted(a for S in self for a in S)
        return SetPartition([[support[-support.index(a) - 1] for a in S]
            for S in pre_conjugate(self)])


class SetPartition(AbstractSetPartition,
                   metaclass=InheritComparisonClasscallMetaclass):
    r"""
    A partition of a set.

    A set partition `p` of a set `S` is a partition of `S` into subsets
    called parts and represented as a set of sets. By extension, a set
    partition of a nonnegative integer `n` is the set partition of the
    integers from 1 to `n`. The number of set partitions of `n` is called
    the `n`-th Bell number.

    There is a natural integer partition associated with a set partition,
    namely the nonincreasing sequence of sizes of all its parts.

    There is a classical lattice associated with all set partitions of
    `n`. The infimum of two set partitions is the set partition obtained
    by intersecting all the parts of both set partitions. The supremum
    is obtained by transitive closure of the relation `i` related to `j`
    if and only if they are in the same part in at least one of the set
    partitions.

    We will use terminology from partitions, in particular the *length* of
    a set partition `A = \{A_1, \ldots, A_k\}` is the number of parts of `A`
    and is denoted by `|A| := k`. The *size* of `A` is the cardinality of `S`.
    We will also sometimes use the notation `[n] := \{1, 2, \ldots, n\}`.

    EXAMPLES:

    There are 5 set partitions of the set `\{1,2,3\}`::

        sage: SetPartitions(3).cardinality()                                            # needs sage.libs.flint
        5

    Here is the list of them::

        sage: SetPartitions(3).list()                                                   # needs sage.graphs
        [{{1, 2, 3}}, {{1, 2}, {3}}, {{1, 3}, {2}}, {{1}, {2, 3}}, {{1}, {2}, {3}}]

    There are 6 set partitions of `\{1,2,3,4\}` whose underlying partition is
    `[2, 1, 1]`::

        sage: SetPartitions(4, [2,1,1]).list()                                          # needs sage.graphs sage.rings.finite_rings
        [{{1}, {2, 4}, {3}},
         {{1}, {2}, {3, 4}},
         {{1, 4}, {2}, {3}},
         {{1, 3}, {2}, {4}},
         {{1, 2}, {3}, {4}},
         {{1}, {2, 3}, {4}}]

    Since :issue:`14140`, we can create a set partition directly by
    :class:`SetPartition`, which creates the base set by taking the
    union of the parts passed in::

        sage: s = SetPartition([[1,3],[2,4]]); s
        {{1, 3}, {2, 4}}
        sage: s.parent()
        Set partitions
    """
    @staticmethod
    def __classcall_private__(cls, parts, check=True):
        """
        Create a set partition from ``parts`` with the appropriate parent.

        EXAMPLES::

            sage: s = SetPartition([[1,3],[2,4]]); s
            {{1, 3}, {2, 4}}
            sage: s.parent()
            Set partitions
        """
        P = SetPartitions()
        return P.element_class(P, parts, check=check)

    def __init__(self, parent, s, check=True):
        """
        Initialize ``self``.

        Internally, a set partition is stored as iterable of blocks,
        sorted by minimal element.

        EXAMPLES::

            sage: S = SetPartitions(4)
            sage: s = S([[1,3],[2,4]])
            sage: TestSuite(s).run()
            sage: SetPartition([])
            {}
        """
        self._latex_options = {}
        try:
            s = sorted(map(frozenset, s), key=min)
        except TypeError:
            s = sorted(map(frozenset, s), key=lambda b: min(str(b)))
        ClonableArray.__init__(self, parent, s, check=check)

    def check(self):
        """
        Check that we are a valid set partition.

        EXAMPLES::

            sage: S = SetPartitions(4)
            sage: s = S([[1, 3], [2, 4]])
            sage: s.check()

        TESTS::

            sage: s = S([[1, 2, 3]], check=False)
            sage: s.check()
            Traceback (most recent call last):
            ...
            ValueError: {{1, 2, 3}} is not an element of Set partitions of {1, 2, 3, 4}

            sage: s = S([1, 2, 3])
            Traceback (most recent call last):
            ...
            TypeError: 'sage.rings.integer.Integer' object is not iterable
        """
        if self not in self.parent():
            raise ValueError(f"{self} is not an element of {self.parent()}")

    def set_latex_options(self, **kwargs):
        r"""
        Set the latex options for use in the ``_latex_`` function.

        - ``tikz_scale`` -- (default: 1) scale for use with tikz package

        - ``plot`` -- (default: ``None``) ``None`` returns the set notation,
          ``linear`` returns a linear plot, ``cyclic`` returns a cyclic
          plot

        - ``color`` -- (default: ``'black'``) the arc colors

        - ``fill`` -- boolean (default: ``False``); if ``True`` then fills
          ``color``, else you can pass in a color to alter the fill color -
          *only works with cyclic plot*

        - ``show_labels`` -- boolean (default: ``True``); if ``True`` shows
          labels (*only works with plots*)

        - ``radius`` -- (default: ``'1cm'``) radius of circle for cyclic
          plot - *only works with cyclic plot*

        - ``angle`` -- (default: 0) angle for linear plot

        EXAMPLES::

            sage: SP = SetPartition([[1,6], [3,5,4]])
            sage: SP.set_latex_options(tikz_scale=2,plot='linear',fill=True,color='blue',angle=45)
            sage: SP.set_latex_options(plot='cyclic')
            sage: SP.latex_options()
            {'angle': 45,
             'color': 'blue',
             'fill': True,
             'plot': 'cyclic',
             'radius': '1cm',
             'show_labels': True,
             'tikz_scale': 2}
        """
        valid_args = ['tikz_scale', 'plot', 'color', 'fill', 'show_labels',
                      'radius', 'angle']

        for key in kwargs:
            if key not in valid_args:
                raise ValueError(f"unknown keyword argument: {key}")
            if key == 'plot':
                if not (kwargs['plot'] == 'cyclic'
                        or kwargs['plot'] == 'linear'
                        or kwargs['plot'] is None):
                    raise ValueError("plot must be None, 'cyclic', or 'linear'")

        self._latex_options.update(kwargs)

    def latex_options(self):
        r"""
        Return the latex options for use in the ``_latex_`` function as a
        dictionary. The default values are set using the global options.

        Options can be found in :meth:`set_latex_options`

        EXAMPLES::

            sage: SP = SetPartition([[1,6], [3,5,4]]); SP.latex_options()
            {'angle': 0,
             'color': 'black',
             'fill': False,
             'plot': None,
             'radius': '1cm',
             'show_labels': True,
             'tikz_scale': 1}
        """
        opts = self._latex_options.copy()
        if "tikz_scale" not in opts:
            opts["tikz_scale"] = 1
        if "plot" not in opts:
            opts["plot"] = None
        if "color" not in opts:
            opts['color'] = 'black'
        if "fill" not in opts:
            opts["fill"] = False
        if "show_labels" not in opts:
            opts['show_labels'] = True
        if "radius" not in opts:
            opts['radius'] = "1cm"
        if "angle" not in opts:
            opts['angle'] = 0
        return opts

    def _latex_(self):
        r"""
        Return a `\LaTeX` string representation of ``self``.

        EXAMPLES::

            sage: x = SetPartition([[1,2], [3,5,4]])
            sage: latex(x)
            \{\{1, 2\}, \{3, 4, 5\}\}

            sage: x.set_latex_options(plot='linear', angle=25, color='red')
            sage: latex(x)
            \begin{tikzpicture}[scale=1]
            \node[below=.05cm] at (0,0) {$1$};
            \node[draw,circle, inner sep=0pt, minimum width=4pt, fill=black] (0) at (0,0) {};
            \node[below=.05cm] at (1,0) {$2$};
            \node[draw,circle, inner sep=0pt, minimum width=4pt, fill=black] (1) at (1,0) {};
            \node[below=.05cm] at (2,0) {$3$};
            \node[draw,circle, inner sep=0pt, minimum width=4pt, fill=black] (2) at (2,0) {};
            \node[below=.05cm] at (3,0) {$4$};
            \node[draw,circle, inner sep=0pt, minimum width=4pt, fill=black] (3) at (3,0) {};
            \node[below=.05cm] at (4,0) {$5$};
            \node[draw,circle, inner sep=0pt, minimum width=4pt, fill=black] (4) at (4,0) {};
            \draw[color=red] (1) to [out=115,in=65] (0);
            \draw[color=red] (3) to [out=115,in=65] (2);
            \draw[color=red] (4) to [out=115,in=65] (3);
            \end{tikzpicture}

            sage: p = SetPartition([['a','c'],['b','d'],['e']])
            sage: p.set_latex_options(plot='cyclic', color='blue', fill=True, tikz_scale=2)
            sage: latex(p)
            \begin{tikzpicture}[scale=2]
            \draw (0,0) circle [radius=1cm];
            \node[label=90:a] (0) at (90:1cm) {};
            \node[label=18:b] (1) at (18:1cm) {};
            \node[label=-54:c] (2) at (-54:1cm) {};
            \node[label=-126:d] (3) at (-126:1cm) {};
            \node[label=-198:e] (4) at (-198:1cm) {};
            \draw[-,thick,color=blue,fill=blue,fill opacity=0.1] ...
            \draw[-,thick,color=blue,fill=blue,fill opacity=0.1] ...
            \draw[-,thick,color=blue,fill=blue,fill opacity=0.1] ...
            \fill[color=black] (0) circle (1.5pt);
            \fill[color=black] (1) circle (1.5pt);
            \fill[color=black] (2) circle (1.5pt);
            \fill[color=black] (3) circle (1.5pt);
            \fill[color=black] (4) circle (1.5pt);
            \end{tikzpicture}
        """
        latex_options = self.latex_options()
        if latex_options["plot"] is None:
            return repr(self).replace("{", r"\{").replace("}", r"\}")

        from sage.misc.latex import latex
        latex.add_package_to_preamble_if_available("tikz")
        res = "\\begin{{tikzpicture}}[scale={}]\n".format(latex_options['tikz_scale'])

        cardinality = self.base_set_cardinality()
        from sage.rings.integer_ring import ZZ
        if all(x in ZZ for x in self.base_set()):
            sort_key = ZZ
        else:
            sort_key = str
        base_set = sorted(self.base_set(), key=sort_key)
        color = latex_options['color']

        # If we want cyclic plots
        if latex_options['plot'] == 'cyclic':
            degrees = 360 // cardinality
            radius = latex_options['radius']

            res += "\\draw (0,0) circle [radius={}];\n".format(radius)

            # Add nodes
            for k, i in enumerate(base_set):
                location = (cardinality - k) * degrees - 270
                if latex_options['show_labels']:
                    res += "\\node[label={}:{}]".format(location, i)
                else:
                    res += "\\node"
                res += " ({}) at ({}:{}) {{}};\n".format(k, location, radius)

            # Setup partitions
            for partition in sorted(self, key=str):
                res += "\\draw[-,thick,color=" + color
                if latex_options['fill'] is not False:
                    if isinstance(latex_options['fill'], str):
                        res += ",fill=" + latex_options['fill']
                    else:
                        res += ",fill={},fill opacity=0.1".format(color)
                res += "] "
                res += " -- ".join("({}.center)".format(base_set.index(j))
                                   for j in sorted(partition, key=sort_key))
                res += " -- cycle;\n"

            # Draw the circles on top
            for k in range(len(base_set)):
                res += "\\fill[color=black] ({}) circle (1.5pt);\n".format(k)

        # If we want line plots
        elif latex_options['plot'] == 'linear':
            angle = latex_options['angle']
            # setup line
            for k, i in enumerate(base_set):
                if latex_options['show_labels']:
                    res += "\\node[below=.05cm] at ({},0) {{${}$}};\n".format(k, i)
                res += "\\node[draw,circle, inner sep=0pt, minimum width=4pt, fill=black] "
                res += "({k}) at ({k},0) {{}};\n".format(k=k)

            # setup arcs
            for partition in sorted(self, key=str):
                p = sorted(partition, key=sort_key)
                if len(p) <= 1:
                    continue
                for k in range(1, len(p)):
                    res += "\\draw[color={}] ({})".format(color, base_set.index(p[k]))
                    res += " to [out={},in={}] ".format(90 + angle, 90 - angle)
                    res += "({});\n".format(base_set.index(p[k - 1]))
        else:
            raise ValueError("plot must be None, 'cyclic', or 'linear'")

        res += "\\end{tikzpicture}"
        return res

    cardinality = ClonableArray.__len__

    size = AbstractSetPartition.base_set_cardinality

    def pipe(self, other):
        r"""
        Return the pipe of the set partitions ``self`` and ``other``.

        The pipe of two set partitions is defined as follows:

        For any integer `k` and any subset `I` of `\ZZ`, let `I + k`
        denote the subset of `\ZZ` obtained by adding `k` to every
        element of `k`.

        If `B` and `C` are set partitions of `[n]` and `[m]`,
        respectively, then the pipe of `B` and `C` is defined as the
        set partition

        .. MATH::

            \{ B_1, B_2, \ldots, B_b,
            C_1 + n, C_2 + n, \ldots, C_c + n \}

        of `[n+m]`, where `B = \{ B_1, B_2, \ldots, B_b \}` and
        `C = \{ C_1, C_2, \ldots, C_c \}`. This pipe is denoted by
        `B | C`.

        EXAMPLES::

            sage: SetPartition([[1,3],[2,4]]).pipe(SetPartition([[1,3],[2]]))
            {{1, 3}, {2, 4}, {5, 7}, {6}}
            sage: SetPartition([]).pipe(SetPartition([[1,2],[3,5],[4]]))
            {{1, 2}, {3, 5}, {4}}
            sage: SetPartition([[1,2],[3,5],[4]]).pipe(SetPartition([]))
            {{1, 2}, {3, 5}, {4}}
            sage: SetPartition([[1,2],[3]]).pipe(SetPartition([[1]]))
            {{1, 2}, {3}, {4}}
        """
        # Note: GIGO if self and other are not standard.
        parts = list(self)
        n = self.base_set_cardinality()
        for newpart in other:
            raised_newpart = Set(i + n for i in newpart)
            parts.append(raised_newpart)
        return SetPartition(parts)

    @combinatorial_map(name='shape')
    def shape(self):
        r"""
        Return the integer partition whose parts are the sizes of the sets
        in ``self``.

        EXAMPLES::

            sage: S = SetPartitions(5)
            sage: x = S([[1,2], [3,5,4]])
            sage: x.shape()
            [3, 2]
            sage: y = S([[2], [3,1], [5,4]])
            sage: y.shape()
            [2, 2, 1]
        """
        return Partition(sorted(map(len, self), reverse=True))

    # we define aliases for shape()
    shape_partition = shape
    to_partition = shape

    @combinatorial_map(name='to permutation')
    def to_permutation(self):
        r"""
        Convert a set partition of `\{1,...,n\}` to a permutation by considering
        the blocks of the partition as cycles.

        The cycles are such that the number of excedences is maximised, that is,
        each cycle is of the form `(a_1,a_2, ...,a_k)` with `a_1<a_2<...<a_k`.

        EXAMPLES::

            sage: s = SetPartition([[1,3],[2,4]])
            sage: s.to_permutation()
            [3, 4, 1, 2]
        """
        return Permutation(tuple(map(tuple, self.standard_form())))

    def to_restricted_growth_word(self, bijection='blocks'):
        r"""
        Convert a set partition of `\{1,...,n\}` to a word of length `n`
        with letters in the nonnegative integers such that each
        letter is at most 1 larger than all the letters before.

        INPUT:

        - ``bijection`` -- (default: ``blocks``) defines the map from
          set partitions to restricted growth functions.  These are
          currently:

          - ``blocks``: :meth:`to_restricted_growth_word_blocks`.

          - ``intertwining``: :meth:`to_restricted_growth_word_intertwining`.

        OUTPUT: a restricted growth word

        .. SEEALSO::

            :meth:`SetPartitions.from_restricted_growth_word`

        EXAMPLES::

            sage: P = SetPartition([[1,4],[2,8],[3,5,6,9],[7]])
            sage: P.to_restricted_growth_word()
            [0, 1, 2, 0, 2, 2, 3, 1, 2]

            sage: P.to_restricted_growth_word("intertwining")
            [0, 1, 2, 2, 1, 0, 3, 3, 2]

            sage: P = SetPartition([[1,2,4,7],[3,9],[5,6,10,11,13],[8],[12]])
            sage: P.to_restricted_growth_word()
            [0, 0, 1, 0, 2, 2, 0, 3, 1, 2, 2, 4, 2]

            sage: P.to_restricted_growth_word("intertwining")
            [0, 0, 1, 1, 2, 0, 1, 3, 3, 3, 0, 4, 1]

        TESTS::

            sage: P = SetPartition([])
            sage: P.to_restricted_growth_word()
            []
            sage: P.to_restricted_growth_word("intertwining")
            []
            sage: S = SetPartitions(5, 2)
            sage: all(S.from_restricted_growth_word(P.to_restricted_growth_word()) == P for P in S)
            True

            sage: S = SetPartitions(5, 2)
            sage: all(S.from_restricted_growth_word(P.to_restricted_growth_word("intertwining"), "intertwining") == P for P in S)
            True
        """
        if bijection == "blocks":
            return self.to_restricted_growth_word_blocks()
        if bijection == "intertwining":
            return self.to_restricted_growth_word_intertwining()
        raise ValueError("the given bijection is not valid")

    def to_restricted_growth_word_blocks(self):
        r"""
        Convert a set partition of `\{1,...,n\}` to a word of length `n`
        with letters in the nonnegative integers such that each
        letter is at most 1 larger than all the letters before.

        The word is obtained by sorting the blocks by their minimal
        element and setting the letters at the positions of the
        elements in the `i`-th block to `i`.

        OUTPUT: a restricted growth word

        .. SEEALSO::

            :meth:`to_restricted_growth_word`
            :meth:`SetPartitions.from_restricted_growth_word`

        EXAMPLES::

            sage: P = SetPartition([[1,4],[2,8],[3,5,6,9],[7]])
            sage: P.to_restricted_growth_word_blocks()
            [0, 1, 2, 0, 2, 2, 3, 1, 2]
        """
        w = [0] * self.size()
        # we can assume that the blocks are sorted by minimal element
        for i, B in enumerate(self):
            for j in B:
                w[j - 1] = i
        return w

    def to_restricted_growth_word_intertwining(self):
        r"""
        Convert a set partition of `\{1,...,n\}` to a word of length `n`
        with letters in the nonnegative integers such that each
        letter is at most 1 larger than all the letters before.

        The `i`-th letter of the word is the numbers of crossings of
        the arc (or half-arc) in the extended arc diagram ending at
        `i`, with arcs (or half-arcs) beginning at a smaller element
        and ending at a larger element.

        OUTPUT: a restricted growth word

        .. SEEALSO::

            :meth:`to_restricted_growth_word`
            :meth:`SetPartitions.from_restricted_growth_word`

        EXAMPLES::

            sage: P = SetPartition([[1,4],[2,8],[3,5,6,9],[7]])
            sage: P.to_restricted_growth_word_intertwining()
            [0, 1, 2, 2, 1, 0, 3, 3, 2]
        """
        A = sorted(self.arcs())
        O = (min(B) for B in self)  # openers
        C = [max(B) for B in self]  # closers
        I = [0] * self.size()
        for i in O:
            I[i - 1] = sum(1 for k, l in A if k < i < l) + sum(1 for k in C if k < i)
        for i, j in A:
            I[j - 1] = sum(1 for k, l in A if i < k < j < l) + sum(1 for k in C if i < k < j)
        return I

    def openers(self):
        """
        Return the minimal elements of the blocks.

        EXAMPLES::

            sage: P = SetPartition([[1,2,4,7],[3,9],[5,6,10,11,13],[8],[12]])
            sage: P.openers()
            [1, 3, 5, 8, 12]
        """
        return sorted([min(B) for B in self])

    def closers(self):
        """
        Return the maximal elements of the blocks.

        EXAMPLES::

            sage: P = SetPartition([[1,2,4,7],[3,9],[5,6,10,11,13],[8],[12]])
            sage: P.closers()
            [7, 8, 9, 12, 13]
        """
        return sorted([max(B) for B in self])

    def to_rook_placement(self, bijection='arcs'):
        r"""
        Return a set of pairs defining a placement of non-attacking rooks
        on a triangular board.

        The cells of the board corresponding to a set partition of
        `\{1,...,n\}` are the pairs `(i,j)` with `0 < i < j < n+1`.

        INPUT:

        - ``bijection`` -- (default: ``arcs``) defines the bijection
          from set partitions to rook placements.  These are
          currently:

          - ``arcs``: :meth:`arcs`
          - ``gamma``: :meth:`to_rook_placement_gamma`
          - ``rho``: :meth:`to_rook_placement_rho`
          - ``psi``: :meth:`to_rook_placement_psi`

        .. SEEALSO::

            :meth:`SetPartitions.from_rook_placement`

        EXAMPLES::

            sage: P = SetPartition([[1,2,4,7],[3,9],[5,6,10,11,13],[8],[12]])
            sage: P.to_rook_placement()
            [(1, 2), (2, 4), (4, 7), (3, 9), (5, 6), (6, 10), (10, 11), (11, 13)]
            sage: P.to_rook_placement("gamma")
            [(1, 4), (3, 5), (4, 6), (5, 8), (7, 11), (8, 9), (10, 12), (12, 13)]
            sage: P.to_rook_placement("rho")
            [(1, 2), (2, 6), (3, 4), (4, 10), (5, 9), (6, 7), (10, 11), (11, 13)]
            sage: P.to_rook_placement("psi")
            [(1, 2), (2, 6), (3, 4), (5, 9), (6, 7), (7, 10), (9, 11), (11, 13)]
        """
        if bijection == "arcs":
            return self.arcs()
        if bijection == "gamma":
            return self.to_rook_placement_gamma()
        if bijection == "rho":
            return self.to_rook_placement_rho()
        if bijection == "psi":
            return self.to_rook_placement_psi()
        raise ValueError("the given map is not valid")

    def to_rook_placement_gamma(self):
        """
        Return the rook diagram obtained by placing rooks according to
        Wachs and White's bijection gamma.

        Note that our index convention differs from the convention in
        [WW1991]_: regarding the rook board as a lower-right
        triangular grid, we refer with `(i,j)` to the cell in the
        `i`-th column from the right and the `j`-th row from the top.

        The algorithm proceeds as follows: non-attacking rooks are
        placed beginning at the left column.  If `n+1-i` is an
        opener, column `i` remains empty.  Otherwise, we place a rook
        into column `i`, such that the number of cells below the
        rook, which are not yet attacked by another rook, equals the
        index of the block to which `n+1-i` belongs.

        OUTPUT: list of coordinates

        .. SEEALSO::

            - :meth:`to_rook_placement`
            - :meth:`SetPartitions.from_rook_placement`
            - :meth:`SetPartitions.from_rook_placement_gamma`

        EXAMPLES::

            sage: P = SetPartition([[1,4],[2,8],[3,5,6,9],[7]])
            sage: P.to_rook_placement_gamma()
            [(1, 3), (2, 7), (4, 5), (5, 6), (6, 9)]

        Figure 5 in [WW1991]_::

            sage: P = SetPartition([[1,2,4,7],[3,9],[5,6,10,11,13],[8],[12]])
            sage: r = P.to_rook_placement_gamma(); r
            [(1, 4), (3, 5), (4, 6), (5, 8), (7, 11), (8, 9), (10, 12), (12, 13)]

        TESTS::

            sage: P = SetPartition([])
            sage: P.to_rook_placement_gamma()
            []
            sage: S = SetPartitions(5, 2)
            sage: all(S.from_rook_placement(P.to_rook_placement("gamma"), "gamma") == P for P in S)
            True
        """
        n = self.size()
        if n == 0:
            return []
        w = self.to_restricted_growth_word_blocks()
        # the set of openers - leftmost occurrences of a letter in w
        EC = sorted([w.index(i) for i in range(max(w) + 1)])
        rooks = []  # pairs (row i, column j)
        R = []  # attacked rows
        for c in range(n):  # columns from left to right
            if c not in EC:
                r = 0
                w_c = w[c]
                while w_c > 0 or r in R:
                    if r not in R:
                        w_c -= 1
                    r += 1
                rooks.append((n - c, n - r))
                R.append(r)
        return sorted(rooks)

    def to_rook_placement_rho(self):
        """
        Return the rook diagram obtained by placing rooks according to
        Wachs and White's bijection rho.

        Note that our index convention differs from the convention in
        [WW1991]_: regarding the rook board as a lower-right
        triangular grid, we refer with `(i,j)` to the cell in the
        `i`-th column from the right and the `j`-th row from the top.

        The algorithm proceeds as follows: non-attacking rooks are
        placed beginning at the top row.  The columns corresponding
        to the closers of the set partition remain empty.  Let `rs_j`
        be the number of closers which are larger than `j` and
        whose block is before the block of `j`.

        We then place a rook into row `j`, such that the number of
        cells to the left of the rook, which are not yet attacked by
        another rook and are not in a column corresponding to a
        closer, equals `rs_j`, unless there are not enough cells in
        this row available, in which case the row remains empty.

        One can show that the precisely those rows which correspond
        to openers of the set partition remain empty.

        OUTPUT: list of coordinates

        .. SEEALSO::

            - :meth:`to_rook_placement`
            - :meth:`SetPartitions.from_rook_placement`
            - :meth:`SetPartitions.from_rook_placement_rho`

        EXAMPLES::

            sage: P = SetPartition([[1,4],[2,8],[3,5,6,9],[7]])
            sage: P.to_rook_placement_rho()
            [(1, 5), (2, 6), (3, 4), (5, 9), (6, 8)]

        Figure 6 in [WW1991]_::

            sage: P = SetPartition([[1,2,4,7],[3,9],[5,6,10,11,13],[8],[12]])
            sage: r = P.to_rook_placement_rho(); r
            [(1, 2), (2, 6), (3, 4), (4, 10), (5, 9), (6, 7), (10, 11), (11, 13)]

            sage: sorted(P.closers() + [i for i, _ in r]) == list(range(1,14))
            True
            sage: sorted(P.openers() + [j for _, j in r]) == list(range(1,14))
            True

        TESTS::

            sage: P = SetPartition([])
            sage: P.to_rook_placement_rho()
            []
            sage: S = SetPartitions(5, 2)
            sage: all(S.from_rook_placement(P.to_rook_placement("rho"), "rho") == P for P in S)
            True
        """
        n = self.size()
        if n == 0:
            return []
        w = self.to_restricted_growth_word_blocks()
        w_rev = w[::-1]
        # the set of closers - rightmost occurrences of a letter in w
        R = sorted([n - w_rev.index(i) - 1 for i in range(max(w) + 1)])
        # the number of closers which are larger than i and whose
        # block is before the block of i
        rs = [sum(1 for j in R if j > i and w[j] < w[i]) for i in range(n)]
        EC = [n - j for j in R]  # empty columns
        rooks = []  # pairs (row i, column j)
        for i in range(1, n):
            U = [j for j in range(n + 1 - i, n + 1) if j not in EC]
            if rs[i] < len(U):
                j = U[rs[i]]
                rooks.append((n + 1 - j, i + 1))
                EC.append(j)
        return sorted(rooks)

    def to_rook_placement_psi(self):
        r"""
        Return the rook diagram obtained by placing rooks according to
        Yip's bijection psi.

        OUTPUT: list of coordinates

        .. SEEALSO::

            - :meth:`to_rook_placement`
            - :meth:`SetPartitions.from_rook_placement`
            - :meth:`SetPartitions.from_rook_placement_psi`

        EXAMPLES:

        Example 36 (arXiv version: Example 4.5) in [Yip2018]_::

            sage: P = SetPartition([[1, 5], [2], [3, 8, 9], [4], [6, 7]])
            sage: P.to_rook_placement_psi()
            [(1, 7), (3, 8), (4, 5), (7, 9)]

        Note that the columns corresponding to the minimal elements
        of the blocks remain empty.

        TESTS::

            sage: P = SetPartition([])
            sage: P.to_rook_placement_psi()
            []
            sage: S = SetPartitions(5,2)
            sage: all(S.from_rook_placement(P.to_rook_placement("psi"), "psi") == P for P in S)
            True
        """
        # Yip draws the diagram as an upper triangular matrix, thus
        # we refer to the cell in row i and column j with (i, j)
        n = self.size()
        degrees = []
        P = [sorted(e) for e in self]
        for j in range(n, 0, -1):
            # find the block number into which c was placed by first
            # removing j and then sorting the blocks
            B = next(B for B in P if B[-1] == j)
            if len(B) == 1:
                P.remove(B)
            else:
                del B[-1]
                P = sorted(P, key=lambda B: (-len(B), min(B)))
                b = P.index(B)
                i = j - b - 1
                degrees.append((j, i))
        # reconstruct rooks from degree sequence
        rooks = []
        attacked_rows = []
        for j, d in reversed(degrees):
            i = 1
            while d > i + sum(1 for r in attacked_rows if r > i):
                i += 1
            attacked_rows.append(i)
            rooks.append((i, j))
        return sorted(rooks)

    def apply_permutation(self, p):
        r"""
        Apply ``p`` to the underlying set of ``self``.

        INPUT:

        - ``p`` -- a permutation

        EXAMPLES::

            sage: x = SetPartition([[1,2], [3,5,4]])
            sage: p = Permutation([2,1,4,5,3])
            sage: x.apply_permutation(p)
            {{1, 2}, {3, 4, 5}}
            sage: q = Permutation([3,2,1,5,4])
            sage: x.apply_permutation(q)
            {{1, 4, 5}, {2, 3}}

            sage: m = PerfectMatching([(1,4),(2,6),(3,5)])
            sage: m.apply_permutation(Permutation([4,1,5,6,3,2]))
            [(1, 2), (3, 5), (4, 6)]
        """
        return self.__class__(self.parent(), [Set(map(p, B)) for B in self])

    def crossings_iterator(self):
        r"""
        Return the crossing arcs of a set partition on a totally ordered set.

        OUTPUT:

        We place the elements of the ground set in order on a
        line and draw the set partition by linking consecutive
        elements of each block in the upper half-plane. This
        function returns an iterator over the pairs of crossing
        lines (as a line correspond to a pair, the iterator
        produces pairs of pairs).

        EXAMPLES::

            sage: p = SetPartition([[1,4],[2,5,7],[3,6]])
            sage: next(p.crossings_iterator())
            ((1, 4), (2, 5))

        TESTS::

            sage: p = SetPartition([]);  p.crossings()
            []
        """
        # each arc is sorted, but the set of arcs might not be
        arcs = sorted(self.arcs(), key=min)
        while arcs:
            i1, j1 = arcs.pop(0)
            for i2, j2 in arcs:
                # we know that i1 < i2 and i1 < j1 and i2 < j2
                if i2 < j1 < j2:
                    yield ((i1, j1), (i2, j2))

    def crossings(self):
        r"""
        Return the crossing arcs of a set partition on a totally ordered set.

        OUTPUT:

        We place the elements of the ground set in order on a
        line and draw the set partition by linking consecutive
        elements of each block in the upper half-plane. This
        function returns a list of the pairs of crossing lines
        (as a line correspond to a pair, it returns a list of
        pairs of pairs).

        EXAMPLES::

            sage: p = SetPartition([[1,4],[2,5,7],[3,6]])
            sage: p.crossings()
            [((1, 4), (2, 5)), ((1, 4), (3, 6)), ((2, 5), (3, 6)), ((3, 6), (5, 7))]

        TESTS::

            sage: p = SetPartition([]);  p.crossings()
            []
        """
        return list(self.crossings_iterator())

    def number_of_crossings(self):
        r"""
        Return the number of crossings.

        OUTPUT:

        We place the elements of the ground set in order on a
        line and draw the set partition by linking consecutive
        elements of each block in the upper half-plane. This
        function returns the number the pairs of crossing lines.

        EXAMPLES::

            sage: p = SetPartition([[1,4],[2,5,7],[3,6]])
            sage: p.number_of_crossings()
            4

            sage: n = PerfectMatching([3,8,1,7,6,5,4,2]); n
            [(1, 3), (2, 8), (4, 7), (5, 6)]
            sage: n.number_of_crossings()
            1
        """
        return Integer(len(list(self.crossings_iterator())))

    def is_noncrossing(self) -> bool:
        r"""
        Check if ``self`` is noncrossing.

        OUTPUT:

        We place the elements of the ground set in order on a
        line and draw the set partition by linking consecutive
        elements of each block in the upper half-plane.  This
        function returns ``True`` if the picture obtained this
        way has no crossings.

        EXAMPLES::

            sage: p = SetPartition([[1,4],[2,5,7],[3,6]])
            sage: p.is_noncrossing()
            False

            sage: n = PerfectMatching([3,8,1,7,6,5,4,2]); n
            [(1, 3), (2, 8), (4, 7), (5, 6)]
            sage: n.is_noncrossing()
            False
            sage: PerfectMatching([(1, 4), (2, 3), (5, 6)]).is_noncrossing()
            True
        """
        it = self.crossings_iterator()
        try:
            next(it)
        except StopIteration:
            return True
        return False

    def nestings_iterator(self):
        r"""
        Iterate over the nestings of ``self``.

        OUTPUT:

        We place the elements of the ground set in order on a
        line and draw the set partition by linking consecutive
        elements of each block in the upper half-plane. This
        function returns an iterator over the pairs of nesting
        lines (as a line correspond to a pair, the iterator
        produces pairs of pairs).

        EXAMPLES::

            sage: n = PerfectMatching([(1, 6), (2, 7), (3, 5), (4, 8)])
            sage: it = n.nestings_iterator()
            sage: next(it)
            ((1, 6), (3, 5))
            sage: next(it)
            ((2, 7), (3, 5))
            sage: next(it)
            Traceback (most recent call last):
            ...
            StopIteration
        """
        # each arc is sorted, but the set of arcs might not be
        arcs = sorted(self.arcs(), key=min)
        while arcs:
            i1, j1 = arcs.pop(0)
            for i2, j2 in arcs:
                # we know that i1 < i2 and i1 < j1 and i2 < j2
                if i2 < j2 < j1:
                    yield ((i1, j1), (i2, j2))

    def nestings(self):
        r"""
        Return the nestings of ``self``.

        OUTPUT:

        We place the elements of the ground set in order on a
        line and draw the set partition by linking consecutive
        elements of each block in the upper half-plane. This
        function returns the list of the pairs of nesting lines
        (as a line correspond to a pair, it returns a list of
        pairs of pairs).

        EXAMPLES::

            sage: m = PerfectMatching([(1, 6), (2, 7), (3, 5), (4, 8)])
            sage: m.nestings()
            [((1, 6), (3, 5)), ((2, 7), (3, 5))]

            sage: n = PerfectMatching([3,8,1,7,6,5,4,2]); n
            [(1, 3), (2, 8), (4, 7), (5, 6)]
            sage: n.nestings()
            [((2, 8), (4, 7)), ((2, 8), (5, 6)), ((4, 7), (5, 6))]

        TESTS::

            sage: m = PerfectMatching([]); m.nestings()
            []
        """
        return list(self.nestings_iterator())

    def number_of_nestings(self):
        r"""
        Return the number of nestings of ``self``.

        OUTPUT:

        We place the elements of the ground set in order on a
        line and draw the set partition by linking consecutive
        elements of each block in the upper half-plane. This
        function returns the number the pairs of nesting lines.

        EXAMPLES::

            sage: n = PerfectMatching([3,8,1,7,6,5,4,2]); n
            [(1, 3), (2, 8), (4, 7), (5, 6)]
            sage: n.number_of_nestings()
            3
        """
        c = Integer(0)
        one = Integer(1)
        for _ in self.nestings_iterator():
            c += one
        return c

    def is_nonnesting(self) -> bool:
        r"""
        Return if ``self`` is nonnesting or not.

        OUTPUT:

        We place the elements of the ground set in order on a
        line and draw the set partition by linking consecutive
        elements of each block in the upper half-plane. This
        function returns ``True`` if the picture obtained this
        way has no nestings.

        EXAMPLES::

            sage: n = PerfectMatching([3,8,1,7,6,5,4,2]); n
            [(1, 3), (2, 8), (4, 7), (5, 6)]
            sage: n.is_nonnesting()
            False
            sage: PerfectMatching([(1, 3), (2, 5), (4, 6)]).is_nonnesting()
            True
        """
        it = self.nestings_iterator()
        try:
            next(it)
        except StopIteration:
            return True
        return False

    def is_atomic(self) -> bool:
        r"""
        Return if ``self`` is an atomic set partition.

        A (standard) set partition `A` can be split if there exist `j < i`
        such that `\max(A_j) < \min(A_i)` where `A` is ordered by minimal
        elements. This means we can write `A = B | C` for some nonempty set
        partitions `B` and `C`. We call a set partition *atomic* if it
        cannot be split and is nonempty. Here, the pipe symbol
        `|` is as defined in method :meth:`pipe`.

        EXAMPLES::

            sage: SetPartition([[1,3], [2]]).is_atomic()
            True
            sage: SetPartition([[1,3], [2], [4]]).is_atomic()
            False
            sage: SetPartition([[1], [2,4], [3]]).is_atomic()
            False
            sage: SetPartition([[1,2,3,4]]).is_atomic()
            True
            sage: SetPartition([[1, 4], [2], [3]]).is_atomic()
            True
            sage: SetPartition([]).is_atomic()
            False
        """
        if len(self) == 0:
            return False
        maximum_so_far = max(self[0])
        for S in self[1:]:
            if maximum_so_far < min(S):
                return False
            maximum_so_far = max(maximum_so_far, max(S))
        return True

    def standardization(self):
        r"""
        Return the standardization of ``self``.

        Given a set partition `A = \{A_1, \ldots, A_n\}` of an ordered
        set `S`, the standardization of `A` is the set partition of
        `\{1, 2, \ldots, |S|\}` obtained by replacing the elements of
        the parts of `A` by the integers `1, 2, \ldots, |S|` in such
        a way that their relative order is preserved (i. e., the
        smallest element in the whole set partition is replaced by
        `1`, the next-smallest by `2`, and so on).

        EXAMPLES::

            sage: SetPartition([[4], [1, 3]]).standardization()
            {{1, 2}, {3}}
            sage: SetPartition([[4], [6, 3]]).standardization()
            {{1, 3}, {2}}
            sage: SetPartition([]).standardization()
            {}
            sage: SetPartition([('c','b'),('d','f'),('e','a')]).standardization()
            {{1, 5}, {2, 3}, {4, 6}}
        """
        r = {e: i for i, e in enumerate(sorted(self.base_set()), 1)}
        return SetPartitions(len(r))([[r[e] for e in b] for b in self])

    def restriction(self, I):
        """
        Return the restriction of ``self`` to a subset ``I``
        (which is given as a set or list or any other iterable).

        EXAMPLES::

            sage: A = SetPartition([[1], [2,3]])
            sage: A.restriction([1,2])
            {{1}, {2}}
            sage: A.restriction([2,3])
            {{2, 3}}
            sage: A.restriction([])
            {}
            sage: A.restriction([4])
            {}
        """
        ret = []
        for part in self:
            newpart = [i for i in part if i in I]
            if len(newpart) != 0:
                ret.append(newpart)
        return SetPartition(ret)

    def ordered_set_partition_action(self, s):
        r"""
        Return the action of an ordered set partition ``s`` on ``self``.

        Let `A = \{A_1, A_2, \ldots, A_k\}` be a set partition of some
        set `S` and `s` be an ordered set partition (i.e., set composition)
        of a subset of `[k]`. Let `A^{\downarrow}` denote the standardization
        of `A`, and `A_{\{ i_1, i_2, \ldots, i_m \}}` denote the sub-partition
        `\{A_{i_1}, A_{i_2}, \ldots, A_{i_m}\}` for any subset
        `\{i_1, \ldots, i_m\}` of `\{1, \ldots, k\}`. We define the set
        partition `s(A)` by

        .. MATH::

            s(A) = A_{s_1}^{\downarrow} | A_{s_2}^{\downarrow} | \cdots
            | A_{s_q}^{\downarrow}.

        where `s = (s_1, s_2, \ldots, s_q)`. Here, the pipe symbol
        `|` is as defined in method :meth:`pipe`.

        This is `s[A]` in section 2.3 in [LM2011]_.

        INPUT:

        - ``s`` -- an ordered set partition with base set a subset
          of `\{1, \ldots, k\}`

        EXAMPLES::

            sage: A = SetPartition([[1], [2,4], [3]])
            sage: s = OrderedSetPartition([[1,3], [2]])
            sage: A.ordered_set_partition_action(s)
            {{1}, {2}, {3, 4}}
            sage: s = OrderedSetPartition([[2,3], [1]])
            sage: A.ordered_set_partition_action(s)
            {{1, 3}, {2}, {4}}

        We create Figure 1 in [LM2011]_ (we note that there is a typo in the
        lower-left corner of the table in the published version of the
        paper, whereas the arXiv version gives the correct partition)::

            sage: A = SetPartition([[1,3], [2,9], [4,5,8], [7]])
            sage: B = SetPartition([[1,3], [2,8], [4,5,6], [7]])
            sage: C = SetPartition([[1,5], [2,8], [3,4,6], [7]])
            sage: s = OrderedSetPartition([[1,3], [2]])
            sage: t = OrderedSetPartition([[2], [3,4]])
            sage: u = OrderedSetPartition([[1], [2,3,4]])
            sage: A.ordered_set_partition_action(s)
            {{1, 2}, {3, 4, 5}, {6, 7}}
            sage: A.ordered_set_partition_action(t)
            {{1, 2}, {3, 4, 6}, {5}}
            sage: A.ordered_set_partition_action(u)
            {{1, 2}, {3, 8}, {4, 5, 7}, {6}}
            sage: B.ordered_set_partition_action(s)
            {{1, 2}, {3, 4, 5}, {6, 7}}
            sage: B.ordered_set_partition_action(t)
            {{1, 2}, {3, 4, 5}, {6}}
            sage: B.ordered_set_partition_action(u)
            {{1, 2}, {3, 8}, {4, 5, 6}, {7}}
            sage: C.ordered_set_partition_action(s)
            {{1, 4}, {2, 3, 5}, {6, 7}}
            sage: C.ordered_set_partition_action(t)
            {{1, 2}, {3, 4, 5}, {6}}
            sage: C.ordered_set_partition_action(u)
            {{1, 2}, {3, 8}, {4, 5, 6}, {7}}

        REFERENCES:

        - [LM2011]_
        """
        cur = 1
        ret = []
        for part in s:
            sub_parts = [list(self[i - 1]) for i in part]  # -1 for indexing
            # Standardizing sub_parts (the cur variable not being reset
            # to 1 gives us the offset we want):
            mins = [min(i) for i in sub_parts]
            over_max = max(map(max, sub_parts)) + 1
            temp = [[] for _ in repeat(None, len(part))]
            while min(mins) != over_max:
                m = min(mins)
                i = mins.index(m)
                temp[i].append(cur)
                cur += 1
                sub_parts[i].pop(sub_parts[i].index(m))
                if len(sub_parts[i]) != 0:
                    mins[i] = min(sub_parts[i])
                else:
                    mins[i] = over_max
            ret += temp
        return SetPartition(ret)

    def refinements(self):
        """
        Return a list of refinements of ``self``.

        .. SEEALSO::

            :meth:`coarsenings`

        EXAMPLES::

            sage: SetPartition([[1,3],[2,4]]).refinements()                             # needs sage.graphs sage.libs.flint
            [{{1, 3}, {2, 4}},
             {{1, 3}, {2}, {4}},
             {{1}, {2, 4}, {3}},
             {{1}, {2}, {3}, {4}}]
            sage: SetPartition([[1],[2,4],[3]]).refinements()                           # needs sage.graphs sage.libs.flint
            [{{1}, {2, 4}, {3}}, {{1}, {2}, {3}, {4}}]
            sage: SetPartition([]).refinements()                                        # needs sage.graphs sage.libs.flint
            [{}]
        """
        L = [SetPartitions(part) for part in self]
        return [SetPartition(sum(map(list, x), [])) for x in itertools.product(*L)]

    def strict_coarsenings(self):
        r"""
        Return all strict coarsenings of ``self``.

        Strict coarsening is the binary relation on set partitions
        defined as the transitive-and-reflexive closure of the
        relation `\prec` defined as follows: For two set partitions
        `A` and `B`, we have `A \prec B` if there exist parts
        `A_i, A_j` of `A` such that `\max(A_i) < \min(A_j)` and
        `B = A \setminus \{A_i, A_j\} \cup \{ A_i \cup A_j \}`.

        EXAMPLES::

            sage: A = SetPartition([[1],[2,3],[4]])
            sage: A.strict_coarsenings()
            [{{1}, {2, 3}, {4}}, {{1, 2, 3}, {4}}, {{1, 4}, {2, 3}},
             {{1}, {2, 3, 4}}, {{1, 2, 3, 4}}]
            sage: SetPartition([[1],[2,4],[3]]).strict_coarsenings()
            [{{1}, {2, 4}, {3}}, {{1, 2, 4}, {3}}, {{1, 3}, {2, 4}}]
            sage: SetPartition([]).strict_coarsenings()
            [{}]
        """
        # This is more or less generic code for computing a
        # transitive-and-reflexive closure by depth-first search.
        todo = [self]
        visited = set([self])
        ret = [self]
        while todo:
            A = todo.pop()
            for i, part in enumerate(A):
                for j, other in enumerate(A[i + 1:]):
                    if max(part) < min(other):
                        next_pi = A[:i]
                        next_pi.append(part.union(other))
                        next_pi += A[i + 1:i + 1 + j] + A[i + j + 2:]
                        next_pi = SetPartition(next_pi)
                        if next_pi not in visited:
                            todo.append(next_pi)
                            visited.add(next_pi)
                            ret.append(next_pi)
        return ret

    def arcs(self):
        r"""
        Return ``self`` as a list of arcs.

        Assuming that the blocks are sorted, the arcs are the pairs
        of consecutive elements in the blocks.

        EXAMPLES::

            sage: A = SetPartition([[1],[2,3],[4]])
            sage: A.arcs()
            [(2, 3)]
            sage: B = SetPartition([[1,3,6,7],[2,5],[4]])
            sage: B.arcs()
            [(1, 3), (3, 6), (6, 7), (2, 5)]
        """
        arcs = []
        for p in self:
            p = sorted(p)
            arcs.extend((p[i], p[i + 1]) for i in range(len(p) - 1))
        return arcs

    def plot(self, angle=None, color='black', base_set_dict=None):
        r"""
        Return a plot of ``self``.

        INPUT:

        - ``angle`` -- (default: `\pi/4`) the angle at which the arcs take off
          (if angle is negative, the arcs are drawn below the horizontal line)

        - ``color`` -- (default: ``'black'``) color of the arcs

        - ``base_set_dict`` -- (optional) dictionary with keys elements
          of :meth:`base_set()` and values as integer or float

        EXAMPLES::

            sage: p = SetPartition([[1,10,11],[2,3,7],[4,5,6],[8,9]])
            sage: p.plot()                                                              # needs sage.plot sage.symbolic
            Graphics object consisting of 29 graphics primitives

        .. PLOT::

            p = SetPartition([[1,10,11],[2,3,7],[4,5,6],[8,9]])
            sphinx_plot(p.plot())

        ::

            sage: p = SetPartition([[1,3,4],[2,5]])
            sage: print(p.plot().description())                                         # needs sage.plot sage.symbolic
            Point set defined by 1 point(s):    [(0.0, 0.0)]
            Point set defined by 1 point(s):    [(1.0, 0.0)]
            Point set defined by 1 point(s):    [(2.0, 0.0)]
            Point set defined by 1 point(s):    [(3.0, 0.0)]
            Point set defined by 1 point(s):    [(4.0, 0.0)]
            Text '1' at the point (0.0,-0.1)
            Text '2' at the point (1.0,-0.1)
            Text '3' at the point (2.0,-0.1)
            Text '4' at the point (3.0,-0.1)
            Text '5' at the point (4.0,-0.1)
            Arc with center (1.0,-1.0) radii (1.41421356237...,1.41421356237...)
             angle 0.0 inside the sector (0.785398163397...,2.35619449019...)
            Arc with center (2.5,-0.5) radii (0.70710678118...,0.70710678118...)
             angle 0.0 inside the sector (0.785398163397...,2.35619449019...)
            Arc with center (2.5,-1.5) radii (2.1213203435...,2.1213203435...)
             angle 0.0 inside the sector (0.785398163397...,2.35619449019...)
            sage: p = SetPartition([['a','c'],['b','d'],['e']])
            sage: print(p.plot().description())                                         # needs sage.plot sage.symbolic
            Point set defined by 1 point(s):  [(0.0, 0.0)]
            Point set defined by 1 point(s):    [(1.0, 0.0)]
            Point set defined by 1 point(s):    [(2.0, 0.0)]
            Point set defined by 1 point(s):    [(3.0, 0.0)]
            Point set defined by 1 point(s):    [(4.0, 0.0)]
            Text 'a' at the point (0.0,-0.1)
            Text 'b' at the point (1.0,-0.1)
            Text 'c' at the point (2.0,-0.1)
            Text 'd' at the point (3.0,-0.1)
            Text 'e' at the point (4.0,-0.1)
            Arc with center (1.0,-1.0) radii (1.41421356237...,1.41421356237...)
             angle 0.0 inside the sector (0.785398163397...,2.35619449019...)
            Arc with center (2.0,-1.0) radii (1.41421356237...,1.41421356237...)
             angle 0.0 inside the sector (0.785398163397...,2.35619449019...)
            sage: p = SetPartition([['a','c'],['b','d'],['e']])
            sage: print(p.plot(base_set_dict={'a':0,'b':1,'c':2,                        # needs sage.plot sage.symbolic
            ....:                             'd':-2.3,'e':5.4}).description())
            Point set defined by 1 point(s):    [(-2.3, 0.0)]
            Point set defined by 1 point(s):    [(0.0, 0.0)]
            Point set defined by 1 point(s):    [(1.0, 0.0)]
            Point set defined by 1 point(s):    [(2.0, 0.0)]
            Point set defined by 1 point(s):    [(5.4, 0.0)]
            Text 'a' at the point (0.0,-0.1)
            Text 'b' at the point (1.0,-0.1)
            Text 'c' at the point (2.0,-0.1)
            Text 'd' at the point (-2.3,-0.1)
            Text 'e' at the point (5.4,-0.1)
            Arc with center (-0.6...,-1.65) radii (2.3334523779...,2.3334523779...)
             angle 0.0 inside the sector (0.785398163397...,2.35619449019...)
            Arc with center (1.0,-1.0) radii (1.4142135623...,1.4142135623...)
             angle 0.0 inside the sector (0.785398163397...,2.35619449019...)
        """
        from sage.plot.graphics import Graphics
        from sage.plot.point import point
        from sage.plot.text import text
        from sage.plot.arc import arc
        from sage.symbolic.constants import pi
        from sage.functions.trig import tan, sin
        from sage.functions.generalized import sgn

        diag = Graphics()
        sorted_vertices_list = list(self.base_set())
        sorted_vertices_list.sort()

        if angle is None:
            angle = pi / 4

        if base_set_dict is not None:
            vertices_dict = base_set_dict
        else:
            vertices_dict = {val: pos for pos, val in enumerate(sorted_vertices_list)}

        for elt in vertices_dict:
            pos = vertices_dict[elt]
            diag += point((pos, 0), size=30, color=color)
            diag += text(elt, (pos, -sgn(angle) * 0.1), color=color)
            # TODO: change 0.1 to something proportional to the height of the picture

        for k, j in self.arcs():
            pos_k, pos_j = float(vertices_dict[k]), float(vertices_dict[j])
            center = ((pos_k + pos_j) / 2,
                      -abs(pos_j - pos_k) / (2 * tan(angle)))
            r1 = abs((pos_j - pos_k) / (2 * sin(angle)))
            sector = (sgn(angle) * (pi / 2 - angle),
                      sgn(angle) * (pi / 2 + angle))
            diag += arc(center=center, r1=r1, sector=sector, color=color)

        diag.axes(False)
        return diag


class SetPartitions(UniqueRepresentation, Parent):
    r"""
    An (unordered) partition of a set `S` is a set of pairwise
    disjoint nonempty subsets with union `S`, and is represented
    by a sorted list of such subsets.

    ``SetPartitions(s)`` returns the class of all set partitions of the set
    ``s``, which can be given as a set or a string; if a string, each
    character is considered an element.

    ``SetPartitions(n)``, where ``n`` is an integer, returns the class of
    all set partitions of the set `\{1, 2, \ldots, n\}`.

    You may specify a second argument `k`. If `k` is an integer,
    :class:`SetPartitions` returns the class of set partitions into `k` parts;
    if it is an integer partition, :class:`SetPartitions` returns the class of
    set partitions whose block sizes correspond to that integer partition.

    The Bell number `B_n`, named in honor of Eric Temple Bell,
    is the number of different partitions of a set with `n` elements.

    EXAMPLES::

        sage: S = [1,2,3,4]
        sage: SetPartitions(S, 2)
        Set partitions of {1, 2, 3, 4} with 2 parts
        sage: SetPartitions([1,2,3,4], [3,1]).list()                                    # needs sage.graphs sage.rings.finite_rings
        [{{1}, {2, 3, 4}}, {{1, 2, 3}, {4}}, {{1, 2, 4}, {3}}, {{1, 3, 4}, {2}}]
        sage: SetPartitions(7, [3,3,1]).cardinality()                                   # needs sage.libs.flint
        70

    In strings, repeated letters are not considered distinct as of
    :issue:`14140`::

        sage: SetPartitions('abcde').cardinality()                                      # needs sage.libs.flint
        52
        sage: SetPartitions('aabcd').cardinality()                                      # needs sage.libs.flint
        15

    If the number of parts exceeds the length of the set,
    an empty iterator is returned (:issue:`37643`)::

        sage: SetPartitions(range(3), 4).list()
        []
        sage: SetPartitions('abcd', 6).list()
        []

    REFERENCES:

    - :wikipedia:`Partition_of_a_set`
    """
    @staticmethod
    def __classcall_private__(cls, s=None, part=None):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: S = SetPartitions(4)
            sage: T = SetPartitions([1,2,3,4])
            sage: S is T
            True
        """
        if s is None:
            return SetPartitions_all()
        if isinstance(s, (int, Integer)):
            s = frozenset(range(1, s + 1))
        else:
            try:
                if s.cardinality() == infinity:
                    raise ValueError("the set must be finite")
            except AttributeError:
                pass
            s = frozenset(s)

        if part is None:
            return SetPartitions_set(s)
        else:
            if isinstance(part, (int, Integer)):
                return SetPartitions_setn(s, part)
            else:
                part = sorted(part, reverse=True)
                if part not in Partitions(len(s)):
                    raise ValueError("part must be an integer partition of %s" % len(s))
                return SetPartitions_setparts(s, Partition(part))

    def __contains__(self, x):
        """
        TESTS::

            sage: S = SetPartitions(4, [2,2])
            sage: SA = SetPartitions()
            sage: all(sp in SA for sp in S)                                             # needs sage.graphs sage.modules sage.rings.finite_rings
            True
            sage: Set([Set([1,2]),Set([3,7])]) in SA                                    # needs sage.graphs
            True
            sage: Set([Set([1,2]),Set([2,3])]) in SA                                    # needs sage.graphs
            False
            sage: Set([]) in SA                                                         # needs sage.graphs
            True
        """
        # x must be a set
        if not isinstance(x, (SetPartition, set, frozenset, Set_generic)):
            return False

        # Check that all parts are disjoint
        base_set = set(e for p in x for e in p)
        if len(base_set) != sum(map(len, x)):
            return False

        # Check to make sure each element of x is a set
        for s in x:
            if not isinstance(s, (set, frozenset, Set_generic)):
                return False

        return True

    def _element_constructor_(self, s, check=True):
        """
        Construct an element of ``self`` from ``s``.

        INPUT:

        - ``s`` -- set of sets

        EXAMPLES::

            sage: S = SetPartitions(4)
            sage: elt = S([[1,3],[2,4]]); elt
            {{1, 3}, {2, 4}}
            sage: P = SetPartitions()
            sage: P(elt).parent() is P
            True
            sage: S = SetPartitions([])
            sage: S([])
            {}
        """
        if isinstance(s, SetPartition):
            if isinstance(s.parent(), SetPartitions):
                return self.element_class(self, s, check=check)
            raise ValueError("cannot convert %s into an element of %s" % (s, self))
        return self.element_class(self, s, check=check)

    Element = SetPartition

    def from_restricted_growth_word(self, w, bijection='blocks'):
        r"""
        Convert a word of length `n` with letters in the nonnegative
        integers such that each letter is at most 1 larger than all
        the letters before to a set partition of `\{1,...,n\}`.

        INPUT:

        - ``w`` -- a restricted growth word

        - ``bijection`` -- (default: ``blocks``) defines the map from
          restricted growth functions to set partitions.  These are
          currently:

          - ``blocks``: .

          - ``intertwining``: :meth:`from_restricted_growth_word_intertwining`.

        OUTPUT: a set partition

        .. SEEALSO::

            :meth:`SetPartition.to_restricted_growth_word`

        EXAMPLES::

            sage: SetPartitions().from_restricted_growth_word([0, 1, 2, 0, 2, 2, 3, 1, 2])
            {{1, 4}, {2, 8}, {3, 5, 6, 9}, {7}}

            sage: SetPartitions().from_restricted_growth_word([0, 0, 1, 0, 2, 2, 0, 3, 1, 2, 2, 4, 2])
            {{1, 2, 4, 7}, {3, 9}, {5, 6, 10, 11, 13}, {8}, {12}}

            sage: SetPartitions().from_restricted_growth_word([0, 0, 1, 0, 2, 2, 0, 3, 1, 2, 2, 4, 2], "intertwining")
            {{1, 2, 6, 7, 9}, {3, 4}, {5, 10, 13}, {8, 11}, {12}}
        """
        if bijection == "blocks":
            return self.from_restricted_growth_word_blocks(w)
        if bijection == "intertwining":
            return self.from_restricted_growth_word_intertwining(w)
        raise ValueError("the given bijection is not valid")

    def from_restricted_growth_word_blocks(self, w):
        r"""
        Convert a word of length `n` with letters in the nonnegative
        integers such that each letter is at most 1 larger than all
        the letters before to a set partition of `\{1,...,n\}`.

        ``w[i]`` is the index of the block containing ``i+1`` when
        sorting the blocks by their minimal element.

        INPUT:

        - ``w`` -- a restricted growth word

        OUTPUT: a set partition

        .. SEEALSO::

            :meth:`from_restricted_growth_word`
            :meth:`SetPartition.to_restricted_growth_word`

        EXAMPLES::

            sage: SetPartitions().from_restricted_growth_word_blocks([0, 0, 1, 0, 2, 2, 0, 3, 1, 2, 2, 4, 2])
            {{1, 2, 4, 7}, {3, 9}, {5, 6, 10, 11, 13}, {8}, {12}}
        """
        R = []
        for i, B in enumerate(w, 1):
            if len(R) <= B:
                R.append([i])
            else:
                R[B].append(i)
        return self.element_class(self, R)

    def from_restricted_growth_word_intertwining(self, w):
        r"""
        Convert a word of length `n` with letters in the nonnegative
        integers such that each letter is at most 1 larger than all
        the letters before to a set partition of `\{1,...,n\}`.

        The `i`-th letter of the word is the numbers of crossings of
        the arc (or half-arc) in the extended arc diagram ending at
        `i`, with arcs (or half-arcs) beginning at a smaller element
        and ending at a larger element.

        INPUT:

        - ``w`` -- a restricted growth word

        OUTPUT: a set partition

        .. SEEALSO::

            :meth:`from_restricted_growth_word`
            :meth:`SetPartition.to_restricted_growth_word`

        EXAMPLES::

            sage: SetPartitions().from_restricted_growth_word_intertwining([0, 0, 1, 0, 2, 2, 0, 3, 1, 2, 2, 4, 2])
            {{1, 2, 6, 7, 9}, {3, 4}, {5, 10, 13}, {8, 11}, {12}}
        """
        if len(w) == 0:
            return self.element_class(self, [])
        R = [[1]]
        C = [1]  # closers, always reverse sorted
        m = 0  # max
        for i in range(1, len(w)):
            if w[i] == 1 + m:  # i+1 is an opener
                m += 1
                R.append([i + 1])
            else:
                # add i+1 to the block, such that there are I[i] closers thereafter
                l = C[w[i]]
                B = next(B for B in R if l in B)
                B.append(i + 1)
                C.remove(l)
            C = [i + 1] + C
        return self.element_class(self, R)

    def from_rook_placement(self, rooks, bijection='arcs', n=None):
        r"""
        Convert a rook placement of the triangular grid to a set
        partition of `\{1,...,n\}`.

        If ``n`` is not given, it is first checked whether it can be
        determined from the parent, otherwise it is the maximal
        occurring integer in the set of rooks.

        INPUT:

        - ``rooks`` -- list of pairs `(i,j)` satisfying
          `0 < i < j < n+1`

        - ``bijection`` -- (default: ``arcs``) defines the map from
          rook placements to set partitions.  These are currently:

          - ``arcs``: :meth:`from_arcs`.
          - ``gamma``: :meth:`from_rook_placement_gamma`.
          - ``rho``: :meth:`from_rook_placement_rho`.
          - ``psi``: :meth:`from_rook_placement_psi`.

        - ``n`` -- (optional) the size of the ground set

        .. SEEALSO::

            :meth:`SetPartition.to_rook_placement`

        EXAMPLES::

            sage: SetPartitions(9).from_rook_placement([[1,4],[2,8],[3,5],[5,6],[6,9]])
            {{1, 4}, {2, 8}, {3, 5, 6, 9}, {7}}

            sage: SetPartitions(13).from_rook_placement([[12,13],[10,12],[8,9],[7,11],[5,8],[4,6],[3,5],[1,4]], "gamma")
            {{1, 2, 4, 7}, {3, 9}, {5, 6, 10, 11, 13}, {8}, {12}}

        TESTS::

            sage: SetPartitions().from_rook_placement([])
            {}
            sage: SetPartitions().from_rook_placement([], "gamma")
            {}
            sage: SetPartitions().from_rook_placement([], "rho")
            {}
            sage: SetPartitions().from_rook_placement([], "psi")
            {}
            sage: SetPartitions().from_rook_placement([], n=2)
            {{1}, {2}}
            sage: SetPartitions().from_rook_placement([], "gamma", 2)
            {{1}, {2}}
            sage: SetPartitions().from_rook_placement([], "rho", 2)
            {{1}, {2}}
            sage: SetPartitions().from_rook_placement([], "psi", 2)
            {{1}, {2}}
        """
        if n is None:
            try:
                n = self.base_set_cardinality()
            except AttributeError:
                n = max((max(r) for r in rooks), default=0)

        if bijection == "arcs":
            return self.from_arcs(rooks, n)
        if bijection == "rho":
            return self.from_rook_placement_rho(rooks, n)
        if bijection == "gamma":
            return self.from_rook_placement_gamma(rooks, n)
        if bijection == "psi":
            return self.from_rook_placement_psi(rooks, n)
        raise ValueError("the given bijection is not valid")

    def from_arcs(self, arcs, n):
        r"""
        Return the coarsest set partition of `\{1,...,n\}` such that any
        two elements connected by an arc are in the same block.

        INPUT:

        - ``n`` -- integer specifying the size of the set
          partition to be produced

        - ``arcs`` -- list of pairs specifying which elements are
          in the same block

        .. SEEALSO::

            - :meth:`from_rook_placement`
            - :meth:`SetPartition.to_rook_placement`
            - :meth:`SetPartition.arcs`

        EXAMPLES::

            sage: SetPartitions().from_arcs([(2,3)], 5)
            {{1}, {2, 3}, {4}, {5}}
        """
        P = DisjointSet(range(1, n + 1))
        for i, j in arcs:
            P.union(i, j)
        return self.element_class(self, P)

    def from_rook_placement_gamma(self, rooks, n):
        r"""
        Return the set partition of `\{1,...,n\}` corresponding to the
        given rook placement by applying Wachs and White's bijection
        gamma.

        Note that our index convention differs from the convention in
        [WW1991]_: regarding the rook board as a lower-right
        triangular grid, we refer with `(i,j)` to the cell in the
        `i`-th column from the right and the `j`-th row from the top.

        INPUT:

        - ``n`` -- integer specifying the size of the set
          partition to be produced

        - ``rooks`` -- list of pairs `(i,j)` such that `0 < i < j < n+1`

        OUTPUT: a set partition

        .. SEEALSO::

            - :meth:`from_rook_placement`
            - :meth:`SetPartition.to_rook_placement`
            - :meth:`SetPartition.to_rook_placement_gamma`

        EXAMPLES:

        Figure 5 in [WW1991]_ concerns the following rook placement::

            sage: r = [(1, 4), (3, 5), (4, 6), (5, 8), (7, 11), (8, 9), (10, 12), (12, 13)]

        Note that the rook `(1, 4)`, translated into Wachs and
        White's convention, is a rook in row 4 from the top and
        column 13 from the left.  The corresponding set partition
        is::

            sage: SetPartitions().from_rook_placement_gamma(r, 13)
            {{1, 2, 4, 7}, {3, 9}, {5, 6, 10, 11, 13}, {8}, {12}}
        """
        if n == 0:
            return self.element_class(self, [])
        # the columns of the board, beginning with column n-1
        C = [set(range(n + 1 - j, n + 1)) for j in range(1, n)]
        # delete cells north and east of each rook
        for j, i in rooks:
            # north
            C[n - j - 1].difference_update(range(j + 1, i + 1))
            # east
            for l in range(n + 1 - j, n + 1):
                C[l - 2].discard(i)
        w = [0] + [len(c) for c in C]
        return self.from_restricted_growth_word_blocks(w)

    def from_rook_placement_rho(self, rooks, n):
        r"""
        Return the set partition of `\{1,...,n\}` corresponding to the
        given rook placement by applying Wachs and White's bijection
        rho.

        Note that our index convention differs from the convention in
        [WW1991]_: regarding the rook board as a lower-right
        triangular grid, we refer with `(i,j)` to the cell in the
        `i`-th column from the right and the `j`-th row from the top.

        INPUT:

        - ``n`` -- integer specifying the size of the set
          partition to be produced

        - ``rooks`` -- list of pairs `(i,j)` such that `0 < i < j < n+1`

        OUTPUT: a set partition

        .. SEEALSO::

            - :meth:`from_rook_placement`
            - :meth:`SetPartition.to_rook_placement`
            - :meth:`SetPartition.to_rook_placement_rho`

        EXAMPLES:

        Figure 5 in [WW1991]_ concerns the following rook placement::

            sage: r = [(1, 2), (2, 6), (3, 4), (4, 10), (5, 9), (6, 7), (10, 11), (11, 13)]

        Note that the rook `(1, 2)`, translated into Wachs and
        White's convention, is a rook in row 2 from the top and
        column 13 from the left.  The corresponding set partition
        is::

            sage: SetPartitions().from_rook_placement_rho(r, 13)
            {{1, 2, 4, 7}, {3, 9}, {5, 6, 10, 11, 13}, {8}, {12}}
        """
        # the closers correspond to the empty columns
        cols = [j for j, _ in rooks]
        R = [j for j in range(1, n + 1) if j not in cols]
        # the columns of the board, beginning with column n-1
        C = [set(range(n + 1 - j, n + 1)) if n - j not in R else set()
             for j in range(1, n)]
        for j, i in rooks:  # column j from right, row i from top
            # south
            C[n - j - 1].difference_update(range(i, n + 1))
            # east
            for l in range(n + 1 - j, n + 1):
                C[l - 2].discard(i)

        C_flat = [i for c in C for i in c]
        # the number of closers which are larger than i and whose
        # block is before the block of i
        rs = [C_flat.count(i) for i in range(1, n + 1)]
        # create the blocks
        P = [[] for _ in R]
        for i in range(1, n + 1):
            k = rs[i - 1]
            # find k-th block which does not yet have a closer
            b = 0
            while k > 0 or (P[b] and P[b][-1] in R):
                if P[b][-1] not in R:
                    k -= 1
                b += 1
            P[b].append(i)
        return self.element_class(self, P)

    def from_rook_placement_psi(self, rooks, n):
        r"""
        Return the set partition of `\{1,...,n\}` corresponding to the
        given rook placement by applying Yip's bijection psi.

        INPUT:

        - ``n`` -- integer specifying the size of the set
          partition to be produced

        - ``rooks`` -- list of pairs `(i,j)` such that `0 < i < j < n+1`

        OUTPUT: a set partition

        .. SEEALSO::

            - :meth:`from_rook_placement`
            - :meth:`SetPartition.to_rook_placement`
            - :meth:`SetPartition.to_rook_placement_psi`

        EXAMPLES:

        Example 36 (arXiv version: Example 4.5) in [Yip2018]_
        concerns the following rook placement::

            sage: r = [(4,5), (1,7), (3, 8), (7,9)]
            sage: SetPartitions().from_rook_placement_psi(r, 9)
            {{1, 5}, {2}, {3, 8, 9}, {4}, {6, 7}}
        """
        # Yip draws the diagram as an upper triangular matrix, thus
        # we refer to the cell in row i and column j with (i, j)
        P = []
        rooks_by_column = {j: i for (i, j) in rooks}
        for c in range(1, n + 1):
            # determine the weight of column c
            try:
                r = rooks_by_column[c]
                n_rooks = 1
                ne = r - 1 + sum(1 for i, j in rooks if i > r and j < c)
            except KeyError:
                n_rooks = 0
                ne = sum(1 for i, j in rooks if j < c)

            b = c - n_rooks - ne
            if len(P) == b - 1:
                P.append([c])
            else:
                P[b - 1].append(c)
            P = sorted(P, key=lambda B: (-len(B), min(B)))
        return self.element_class(self, P)

    def is_less_than(self, s, t) -> bool:
        r"""
        Check if `s < t` in the refinement ordering on set partitions.

        This means that `s` is a refinement of `t` and satisfies
        `s \neq t`.

        A set partition `s` is said to be a refinement of a set
        partition `t` of the same set if and only if each part of
        `s` is a subset of a part of `t`.

        EXAMPLES::

            sage: S = SetPartitions(4)
            sage: s = S([[1,3],[2,4]])
            sage: t = S([[1],[2],[3],[4]])
            sage: S.is_less_than(t, s)
            True
            sage: S.is_less_than(s, t)
            False
            sage: S.is_less_than(s, s)
            False
        """
        if hasattr(s.parent(), "_set"):
            S = s.parent()._set
        else:
            S = s.base_set()
        if hasattr(t.parent(), "_set"):
            T = t.parent()._set
        else:
            T = t.base_set()
        if S != T:
            raise ValueError("cannot compare partitions of different sets")

        if s == t:
            return False

        for p in s:
            x = next(iter(p))
            for t_ in t:
                if x in t_:
                    break
            for p_ in p:
                if p_ not in t_:
                    return False
        return True

    lt = is_less_than

    def is_strict_refinement(self, s, t) -> bool:
        r"""
        Return ``True`` if ``s`` is a strict refinement of ``t`` and
        satisfies `s \neq t`.

        A set partition `s` is said to be a strict refinement of a set
        partition `t` of the same set if and only if one can obtain
        `t` from `s` by repeatedly combining pairs of parts whose
        convex hulls don't intersect (i. e., whenever we are combining
        two parts, the maximum of each of them should be smaller than
        the minimum of the other).

        EXAMPLES::

            sage: S = SetPartitions(4)
            sage: s = S([[1],[2],[3],[4]])
            sage: t = S([[1,3],[2,4]])
            sage: u = S([[1,2,3,4]])
            sage: S.is_strict_refinement(s, t)
            True
            sage: S.is_strict_refinement(t, u)
            False
            sage: A = SetPartition([[1,3],[2,4]])
            sage: B = SetPartition([[1,2,3,4]])
            sage: S.is_strict_refinement(s, A)
            True
            sage: S.is_strict_refinement(t, B)
            False
        """
        if hasattr(s.parent(), "_set"):
            S = s.parent()._set
        else:
            S = frozenset(s.base_set())
        if hasattr(t.parent(), "_set"):
            T = t.parent()._set
        else:
            T = frozenset(t.base_set())
        if S != T:
            raise ValueError("cannot compare partitions of different sets")

        if s == t:
            return False

        for p in t:
            L = [x for x in s if x.issubset(p)]
            if sum(len(x) for x in L) != len(p) \
                    or any(max(L[i]) > min(L[i + 1])
                           for i in range(len(L) - 1)):
                return False
        return True


class SetPartitions_all(SetPartitions):
    r"""
    All set partitions.
    """

    def __init__(self):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: S = SetPartitions()
            sage: TestSuite(S).run()
        """
        SetPartitions.__init__(self, category=InfiniteEnumeratedSets())

    def subset(self, size=None, **kwargs):
        """
        Return the subset of set partitions of a given size and
        additional keyword arguments.

        EXAMPLES::

            sage: P = SetPartitions()
            sage: P.subset(4)
            Set partitions of {1, 2, 3, 4}
        """
        if size is None:
            return self
        return SetPartitions(size, **kwargs)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: SetPartitions()
            Set partitions
        """
        return "Set partitions"

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: it = SetPartitions().__iter__()
            sage: [next(it) for x in range(10)]
            [{}, {{1}}, {{1, 2}}, {{1}, {2}}, {{1, 2, 3}}, {{1, 2}, {3}},
             {{1, 3}, {2}}, {{1}, {2, 3}}, {{1}, {2}, {3}}, {{1, 2, 3, 4}}]
        """
        n = 0
        while True:
            for x in SetPartitions_set(frozenset(range(1, n + 1))):
                yield self.element_class(self, list(x))
            n += 1


class SetPartitions_set(SetPartitions):
    """
    Set partitions of a fixed set `S`.
    """
    @staticmethod
    def __classcall_private__(cls, s):
        """
        Normalize ``s`` to ensure a unique representation.

        EXAMPLES::

            sage: S1 = SetPartitions(set([2,1,4]))
            sage: S2 = SetPartitions([4,1,2])
            sage: S3 = SetPartitions((1,2,4))
            sage: S1 is S2, S1 is S3
            (True, True)
        """
        return super().__classcall__(cls, frozenset(s))

    def __init__(self, s):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: S = SetPartitions(3)
            sage: TestSuite(S).run()
            sage: SetPartitions(0).list()
            [{}]
            sage: SetPartitions([]).list()
            [{}]
        """
        self._set = s
        SetPartitions.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self):
        """
        TESTS::

            sage: SetPartitions([1,2,3])
            Set partitions of {1, 2, 3}
        """
        return "Set partitions of %s" % (Set(self._set))

    def __contains__(self, x):
        """
        TESTS::

            sage: S = SetPartitions(4, [2,2])
            sage: all(sp in S for sp in S)                                              # needs sage.graphs sage.rings.finite_rings
            True
            sage: SetPartition([[1,3],[2,4]]) in SetPartitions(3)                       # needs sage.graphs
            False
            sage: SetPartition([[1,3],[2,4]]) in SetPartitions(4, [3,1])                # needs sage.graphs
            False
            sage: SetPartition([[2],[1,3,4]]) in SetPartitions(4, [3,1])                # needs sage.graphs
            True
        """
        # Must pass the general check
        if not SetPartitions.__contains__(self, x):
            return False

        # Make sure that the number of elements match up
        if sum(map(len, x)) != len(self._set):
            return False

        # Make sure that the union of all the sets is the original set
        return Set(e for p in x for e in p) == Set(self._set)

    def random_element(self):
        r"""
        Return a random set partition.

        This is a very naive implementation of Knuths outline in F3B,
        7.2.1.5.

        EXAMPLES::

            sage: S = SetPartitions(10)
            sage: s = S.random_element()                                                # needs sage.symbolic
            sage: s.parent() is S                                                       # needs sage.symbolic
            True
            sage: assert s in S, s                                                      # needs sage.symbolic

            sage: S = SetPartitions(["a", "b", "c"])
            sage: s = S.random_element()                                                # needs sage.symbolic
            sage: s.parent() is S                                                       # needs sage.symbolic
            True
            sage: assert s in S, s                                                      # needs sage.symbolic
        """
        base_set = list(self.base_set())
        N = len(base_set)
        from sage.symbolic.constants import e
        c = float(e) * bell_number(N)
        # it would be much better to generate M in the way Knuth
        # recommends, the following is a waste
        G = GeneralDiscreteDistribution([float(m)**N / (c * factorial(m)) for m in range(4 * N)])
        M = G.get_random_element() - 1
        l = (randint(0, M) for i in range(N))
        p = {}
        for i, b in enumerate(l):
            if b in p:
                p[b].append(base_set[i])
            else:
                p[b] = [base_set[i]]

        return self.element_class(self, p.values(), check=False)

    def cardinality(self):
        """
        Return the number of set partitions of the set `S`.

        The cardinality is given by the `n`-th Bell number where `n` is the
        number of elements in the set `S`.

        EXAMPLES::

            sage: # needs sage.libs.flint
            sage: SetPartitions([1,2,3,4]).cardinality()
            15
            sage: SetPartitions(3).cardinality()
            5
            sage: SetPartitions(3,2).cardinality()
            3
            sage: SetPartitions([]).cardinality()
            1
        """
        return bell_number(len(self._set))

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: SetPartitions(3).list()
            [{{1, 2, 3}}, {{1, 2}, {3}}, {{1, 3}, {2}}, {{1}, {2, 3}}, {{1}, {2}, {3}}]

            sage: SetPartitions(["a", "b"]).list()
            [{{'a', 'b'}}, {{'a'}, {'b'}}]
        """
        try:
            s = sorted(self._set)
        except TypeError:
            s = sorted(self._set, key=str)
        for sp in set_partition_iterator(s):
            yield self.element_class(self, sp, check=False)

    def base_set(self):
        """
        Return the base set of ``self``.

        EXAMPLES::

            sage: SetPartitions(3).base_set()
            {1, 2, 3}

            sage: sorted(SetPartitions(["a", "b", "c"]).base_set())
            ['a', 'b', 'c']
        """
        return Set(self._set)

    def base_set_cardinality(self):
        """
        Return the cardinality of the base set of ``self``.

        EXAMPLES::

            sage: SetPartitions(3).base_set_cardinality()
            3
        """
        return len(self._set)


class SetPartitions_setparts(SetPartitions_set):
    r"""
    Set partitions with fixed partition sizes corresponding to an
    integer partition `\lambda`.
    """
    @staticmethod
    def __classcall_private__(cls, s, parts):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: S = SetPartitions(4, [2,2])
            sage: T = SetPartitions([1,2,3,4], Partition([2,2]))
            sage: S is T
            True

            sage: S = SetPartitions(4, [3,1])
            sage: T = SetPartitions(4, (1,3))
            sage: S is T
            True
        """
        if isinstance(s, (int, Integer)):
            s = list(range(1, s + 1))
        return super().__classcall__(cls, frozenset(s), Partition(parts))

    def __init__(self, s, parts):
        """
        Initialize the data structure.

        We can assume here that ``parts`` is a :class:`Partition`.

        TESTS::

            sage: S = SetPartitions(4, [2,2])
            sage: TestSuite(S).run()                                                    # needs sage.graphs sage.libs.flint
        """
        SetPartitions_set.__init__(self, s)
        self._parts = parts

    def _repr_(self):
        """
        TESTS::

            sage: SetPartitions(4, [2,2])
            Set partitions of {1, 2, 3, 4} with sizes in [2, 2]
        """
        return "Set partitions of %s with sizes in %s" % (Set(self._set), self._parts)

    def shape(self):
        r"""
        Return the partition of block sizes of the set partitions in ``self``.

        EXAMPLES::

            sage: SetPartitions(5, [2,2,1]).shape()
            [2, 2, 1]
        """
        return self._parts

    def cardinality(self):
        r"""
        Return the cardinality of ``self``.

        This algorithm counts for each block of the partition the
        number of ways to fill it using values from the set.  Then,
        for each distinct value `v` of block size, we divide the result by
        the number of ways to arrange the blocks of size `v` in the
        set partition.

        For example, if we want to count the number of set partitions
        of size 13 having [3,3,3,2,2] as underlying partition we
        compute the number of ways to fill each block of the
        partition, which is `\binom{13}{3} \binom{10}{3} \binom{7}{3}
        \binom{4}{2}\binom{2}{2}` and as we have three blocks of size
        `3` and two blocks of size `2`, we divide the result by
        `3!2!` which gives us `600600`.

        EXAMPLES::

            sage: SetPartitions(3, [2,1]).cardinality()
            3
            sage: SetPartitions(13, Partition([3,3,3,2,2])).cardinality()
            600600

        TESTS::

            sage: all((len(SetPartitions(size, part)) == SetPartitions(size, part).cardinality()
            ....:     for size in range(8) for part in Partitions(size)))
            True
            sage: sum((SetPartitions(13, p).cardinality()                               # needs sage.libs.flint
            ....:     for p in Partitions(13))) == SetPartitions(13).cardinality()
            True
        """
        from sage.misc.misc_c import prod

        remaining_subset_size = Integer(len(self._set))
        cardinal = Integer(1)
        for subset_size in self._parts:
            cardinal *= remaining_subset_size.binomial(subset_size)
            remaining_subset_size -= subset_size

        repetitions = (Integer(rep).factorial()
                       for rep in self._parts.to_exp_dict().values()
                       if rep != 1)
        cardinal /= prod(repetitions)
        return Integer(cardinal)

    def _set_partition_poset(self):
        r"""
        Return the Hasse diagram of a poset whose linear extensions correspond
        to the set partitions with specified block sizes.

        TESTS::

            sage: P = SetPartitions(["a", "b", "c", "d", "e"],                          # needs sage.graphs
            ....:                   [2,2,1])._set_partition_poset()
            sage: P.cover_relations()                                                   # needs sage.graphs
            [(1, 2), (1, 3), (3, 4)]

            sage: n = 9
            sage: all(SetPartitions(n, mu).cardinality() ==                             # needs sage.graphs sage.modules
            ....:     len(list(SetPartitions(n, mu)._set_partition_poset().linear_extensions()))
            ....:     for mu in Partitions(n))
            True
        """
        c = self._parts.to_exp_dict()
        covers = {}
        i = 0
        for s in sorted(c):
            # s is the block size
            # each block is one tree in the poset
            for m in range(c[s]):
                # m is the multiplicity of blocks with size s
                #
                # the first element in each non-final block has an
                # additional cover
                first = i
                if s == 1:
                    covers[i] = []
                else:
                    for _ in range(s - 1):
                        covers[i] = [i + 1]
                        i += 1
                i += 1
                if m < c[s] - 1:
                    covers[first].append(i)
        return HasseDiagram(covers)

    def __iter__(self):
        """
        An iterator for all the set partitions of the given set with
        the given sizes.

        EXAMPLES::

            sage: SetPartitions(3, [2,1]).list()                                        # needs sage.graphs sage.rings.finite_rings
            [{{1}, {2, 3}}, {{1, 2}, {3}}, {{1, 3}, {2}}]

            sage: SetPartitions(["a", "b", "c"], [2,1]).list()                          # needs sage.graphs sage.rings.finite_rings
            [{{'a'}, {'b', 'c'}}, {{'a', 'b'}, {'c'}}, {{'a', 'c'}, {'b'}}]

        TESTS::

            sage: n = 8
            sage: all(SetPartitions(n, mu).cardinality()                                # needs sage.graphs sage.rings.finite_rings
            ....:      == len(list(SetPartitions(n, mu))) for mu in Partitions(n))
            True
        """
        # Ruskey, Combinatorial Generation, sec. 5.10.1 and Knuth TAOCP 4A 7.2.1.5, Exercise 6
        k = len(self._parts)
        n = len(self._set)
        P = self._set_partition_poset()
        try:
            s = sorted(self._set)
        except TypeError:
            s = sorted(self._set, key=str)

        sums = [0]
        for b in sorted(self._parts):
            sums.append(sums[-1] + b)

        for ext in P.linear_extensions():
            pi = [None] * n
            for i in range(n):
                pi[ext[i]] = s[i]
            sp = [[pi[j] for j in range(sums[i], sums[i + 1])]
                  for i in range(k)]
            yield self.element_class(self, sp, check=False)

    def __contains__(self, x):
        """
        Check containment.

        TESTS::

            sage: S = SetPartitions(4, [3,1])
            sage: Set([Set([1,2,3]), Set([4])]) in S
            True
            sage: Set([Set([1,3]), Set([2,4])]) in S
            False
            sage: Set([Set([1,2,3,4])]) in S
            False
        """
        if not SetPartitions_set.__contains__(self, x):
            return False
        return sorted(map(len, x), reverse=True) == self._parts

    def random_element(self):
        r"""
        Return a random set partition of ``self``.

        ALGORITHM:

        Based on the cardinality method. For each block size `k_i`,
        we choose a uniformly random subset `X_i \subseteq S_i` of
        size `k_i` of the elements `S_i` that have not yet been selected.
        Thus, we define `S_{i+1} = S_i \setminus X_i` with `S_i = S`
        being the defining set. This is not yet proven to be uniformly
        distributed, but numerical tests show this is likely uniform.

        EXAMPLES::

            sage: S = SetPartitions(10, [4,3,2,1])
            sage: s = S.random_element()
            sage: s.parent() is S
            True
            sage: assert s in S, s

            sage: S = SetPartitions(["a", "b", "c", "d"], [2,2])
            sage: s = S.random_element()
            sage: s.parent() is S
            True
            sage: assert s in S, s
        """
        base_set = list(self.base_set())
        N = len(base_set)
        ret = []
        for p in self._parts:
            X = sample(range(N), p)
            ret.append([base_set[i] for i in X])
            for i in sorted(X, reverse=True):
                del base_set[i]
            N -= p

        return self.element_class(self, ret, check=False)


class SetPartitions_setn(SetPartitions_set):
    """
    Set partitions with a given number of blocks.
    """
    @staticmethod
    def __classcall_private__(cls, s, k):
        """
        Normalize ``s`` to ensure a unique representation.

        EXAMPLES::

            sage: S1 = SetPartitions(set([2,1,4]), 2)
            sage: S2 = SetPartitions([4,1,2], 2)
            sage: S3 = SetPartitions((1,2,4), 2)
            sage: S1 is S2, S1 is S3
            (True, True)
        """
        return super().__classcall__(cls, frozenset(s), k)

    def __init__(self, s, k):
        """
        TESTS::

            sage: S = SetPartitions(5, 3)
            sage: TestSuite(S).run()
        """
        self._k = k
        SetPartitions_set.__init__(self, s)

    def _repr_(self):
        """
        TESTS::

            sage: SetPartitions(5, 3)
            Set partitions of {1, 2, 3, 4, 5} with 3 parts
        """
        return "Set partitions of %s with %s parts" % (Set(self._set), self._k)

    def number_of_blocks(self):
        r"""
        Return the number of blocks of the set partitions in ``self``.

        EXAMPLES::

            sage: SetPartitions(5, 3).number_of_blocks()
            3
        """
        return self._k

    def cardinality(self):
        """
        The Stirling number of the second kind is the number of partitions
        of a set of size `n` into `k` blocks.

        EXAMPLES::

            sage: SetPartitions(5, 3).cardinality()
            25
            sage: stirling_number2(5,3)
            25
        """
        return stirling2(len(self._set), self._k)

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: SetPartitions(4, 2).list()
            [{{1, 3, 4}, {2}},
             {{1, 4}, {2, 3}},
             {{1, 2, 4}, {3}},
             {{1, 3}, {2, 4}},
             {{1}, {2, 3, 4}},
             {{1, 2}, {3, 4}},
             {{1, 2, 3}, {4}}]

            sage: SetPartitions(["a", "b", "c"], 2).list()
            [{{'a', 'c'}, {'b'}}, {{'a'}, {'b', 'c'}}, {{'a', 'b'}, {'c'}}]
        """
        try:
            s = sorted(self._set)
        except TypeError:
            s = sorted(self._set, key=str)
        for sp in set_partition_iterator_blocks(s, self._k):
            yield self.element_class(self, sp, check=False)

    def __contains__(self, x):
        """
        Check containment.

        TESTS::

            sage: S = SetPartitions(4, 2)
            sage: Set([Set([1,2,3]), Set([4])]) in S
            True
            sage: Set([Set([1,3]), Set([2,4])]) in S
            True
            sage: Set([Set([1,2,3,4])]) in S
            False
        """
        if not SetPartitions_set.__contains__(self, x):
            return False
        return len(x) == self._k

    def random_element(self):
        r"""
        Return a random set partition of ``self``.

        See https://mathoverflow.net/questions/141999.

        EXAMPLES::

            sage: S = SetPartitions(10, 4)
            sage: s = S.random_element()
            sage: s.parent() is S
            True
            sage: assert s in S, s

            sage: S = SetPartitions(["a", "b", "c"], 2)
            sage: s = S.random_element()
            sage: s.parent() is S
            True
            sage: assert s in S, s
        """
        def re(N, k):
            if N == 0:
                return [[]]
            if N == 1:
                return [[0]]
            if stirling2(N - 1, k - 1) > random() * stirling2(N, k):
                return [[N - 1]] + re(N - 1, k - 1)

            p = re(N - 1, k)
            p[randint(0, len(p) - 1)].append(N - 1)
            return p

        base_set = list(self.base_set())
        N = len(base_set)
        k = self._k
        p = re(N, k)
        return self.element_class(self, [[base_set[e] for e in b] for b in p], check=False)


def cyclic_permutations_of_set_partition(set_part):
    """
    Return all combinations of cyclic permutations of each cell of the
    set partition.

    AUTHORS:

    - Robert L. Miller

    EXAMPLES::

        sage: from sage.combinat.set_partition import cyclic_permutations_of_set_partition
        sage: cyclic_permutations_of_set_partition([[1,2,3,4],[5,6,7]])
        [[[1, 2, 3, 4], [5, 6, 7]],
         [[1, 2, 4, 3], [5, 6, 7]],
         [[1, 3, 2, 4], [5, 6, 7]],
         [[1, 3, 4, 2], [5, 6, 7]],
         [[1, 4, 2, 3], [5, 6, 7]],
         [[1, 4, 3, 2], [5, 6, 7]],
         [[1, 2, 3, 4], [5, 7, 6]],
         [[1, 2, 4, 3], [5, 7, 6]],
         [[1, 3, 2, 4], [5, 7, 6]],
         [[1, 3, 4, 2], [5, 7, 6]],
         [[1, 4, 2, 3], [5, 7, 6]],
         [[1, 4, 3, 2], [5, 7, 6]]]
    """
    return list(cyclic_permutations_of_set_partition_iterator(set_part))


def cyclic_permutations_of_set_partition_iterator(set_part):
    """
    Iterate over all combinations of cyclic permutations of each cell
    of the set partition.

    AUTHORS:

    - Robert L. Miller

    EXAMPLES::

        sage: from sage.combinat.set_partition import cyclic_permutations_of_set_partition_iterator
        sage: list(cyclic_permutations_of_set_partition_iterator([[1,2,3,4],[5,6,7]]))
        [[[1, 2, 3, 4], [5, 6, 7]],
         [[1, 2, 4, 3], [5, 6, 7]],
         [[1, 3, 2, 4], [5, 6, 7]],
         [[1, 3, 4, 2], [5, 6, 7]],
         [[1, 4, 2, 3], [5, 6, 7]],
         [[1, 4, 3, 2], [5, 6, 7]],
         [[1, 2, 3, 4], [5, 7, 6]],
         [[1, 2, 4, 3], [5, 7, 6]],
         [[1, 3, 2, 4], [5, 7, 6]],
         [[1, 3, 4, 2], [5, 7, 6]],
         [[1, 4, 2, 3], [5, 7, 6]],
         [[1, 4, 3, 2], [5, 7, 6]]]
    """
    from sage.combinat.permutation import CyclicPermutations
    if len(set_part) == 1:
        for i in CyclicPermutations(set_part[0]):
            yield [i]
    else:
        for right in cyclic_permutations_of_set_partition_iterator(set_part[1:]):
            for perm in CyclicPermutations(set_part[0]):
                yield [perm] + right
