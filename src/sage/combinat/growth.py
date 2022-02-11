# -*- coding: utf-8 -*-
r"""
Growth diagrams and dual graded graphs

AUTHORS:

- Martin Rubey (2016-09): Initial version
- Martin Rubey (2017-09): generalize, more rules, improve documentation
- Travis Scrimshaw (2017-09): switch to rule-based framework

.. TODO::

    - provide examples for the P and Q-symbol in the skew case
    - implement a method providing a visualization of the growth
      diagram with all labels, perhaps as LaTeX code
    - when shape is given, check that it is compatible with filling
      or labels
    - optimize rules, mainly for :class:`RuleRSK` and
      :class:`RuleBurge`
    - implement backward rules for :class:`GrowthDiagram.rules.Domino`
    - implement backward rule from [LLMSSZ2013]_, [LS2007]_
    - make semistandard extension generic
    - accommodate dual filtered graphs

A guided tour
=============

Growth diagrams, invented by Sergey Fomin [Fom1994]_, [Fom1995]_,
provide a vast generalization of the Robinson-Schensted-Knuth (RSK)
correspondence between matrices with non-negative integer entries and
pairs of semistandard Young tableaux of the same shape.

The main fact is that many correspondences similar to RSK can be
defined by providing a pair of so-called local rules: a 'forward'
rule, whose input are three vertices `y`, `t` and `x` of a certain
directed graph (in the case of Robinson-Schensted: the directed graph
corresponding to Young's lattice) and an integer (in the case of
Robinson-Schensted: `0` or `1`), and whose output is a fourth vertex
`z`.  This rule should be invertible in the following sense: there is
a so-called 'backward' rule that recovers the integer and `t` given
`y`, `z` and `x`.

As an example, the growth rules for the classical RSK correspondence
are provided by :class:`RuleRSK`.  To produce a growth diagram, pass
the desired rule and a permutation to :class:`GrowthDiagram`::

    sage: RuleRSK = GrowthDiagram.rules.RSK()
    sage: w = [2,3,6,1,4,5]; G = GrowthDiagram(RuleRSK, w); G
    0  0  0  1  0  0
    1  0  0  0  0  0
    0  1  0  0  0  0
    0  0  0  0  1  0
    0  0  0  0  0  1
    0  0  1  0  0  0

The forward rule just mentioned assigns 49 partitions to the corners
of each of the 36 cells of this matrix (i.e., 49 the vertices of a
`(6+1) \times (6+1)` grid graph), with the exception of the corners
on the left and top boundary, which are initialized with the empty
partition. More precisely, for each cell, the
:meth:`~sage.combinat.growth.RuleRSK.forward_rule` computes the
partition `z` labelling the lower right corner, given the content `c`
of a cell and the other three partitions::

    t --- x
    |  c  |
    y --- z

.. WARNING::

    Note that a growth diagram is printed with matrix coordinates,
    the origin being in the top-left corner.  Therefore, the growth
    is from the top left to the bottom right!

The partitions along the boundary opposite of the origin, reading
from the bottom left to the top right, are obtained by using the
method :meth:`~sage.combinat.growth.GrowthDiagram.out_labels`::

    sage: G.out_labels()
    [[],
     [1],
     [2],
     [3],
     [3, 1],
     [3, 2],
     [4, 2],
     [4, 1],
     [3, 1],
     [2, 1],
     [1, 1],
     [1],
     []]

However, in the case of a rectangular filling, it is more practical
to split this sequence of labels in two.  Interpreting the sequence
of partitions along the right boundary as a standard Young tableau,
we then obtain the so-called
:meth:`~sage.combinat.growth.RulePartitions.P_symbol`, the partitions
along the bottom boundary yield the so-called
:meth:`~sage.combinat.growth.RulePartitions.Q_symbol`.  These
coincide with the output of the classical
:func:`~sage.combinat.rsk.RSK` insertion algorithm::

    sage: ascii_art([G.P_symbol(), G.Q_symbol()])
    [   1  3  4  5    1  2  3  6 ]
    [   2  6      ,   4  5       ]
    sage: ascii_art(RSK(w))
    [   1  3  4  5    1  2  3  6 ]
    [   2  6      ,   4  5       ]

The filling can be recovered knowing the partitions labelling the
corners of the bottom and the right boundary alone, by repeatedly
applying the :meth:`~sage.combinat.growth.RuleRSK.backward_rule`.
Therefore, to initialize a :class:`GrowthDiagram`, we can provide
these labels instead of the filling::

    sage: GrowthDiagram(RuleRSK, labels=G.out_labels())
    0  0  0  1  0  0
    1  0  0  0  0  0
    0  1  0  0  0  0
    0  0  0  0  1  0
    0  0  0  0  0  1
    0  0  1  0  0  0

Invocation
----------

In general, growth diagrams are defined for `0-1`-fillings of
arbitrary skew shapes.  In the case of the Robinson-Schensted-Knuth
correspondence, even arbitrary non-negative integers are allowed.  In
other cases, entries may be either zero or an `r`-th root of unity -
for example, :class:`~sage.combinat.growth.RuleDomino` insertion is
defined for signed permutations, that is, `r=2`.  Traditionally, words
and permutations are also used to specify a filling in special cases.

To accommodate all this, the filling may be passed in various ways.
The most general possibility is to pass a dictionary of coordinates
to (signed) entries, where zeros can be omitted.  In this case, when
the parameter ``shape`` is not explicitly specified, it is assumed
to be the minimal rectangle containing the origin and all coordinates
with non-zero entries.

For example, consider the following generalized permutation::

    1 2 2 2 4 4
    4 2 3 3 2 3

that we encode as the dictionary::

    sage: P = {(1-1,4-1): 1, (2-1,2-1): 1, (2-1,3-1): 2, (4-1,2-1): 1, (4-1,3-1): 1}

Note that we are subtracting `1` from all entries because of
zero-based indexing, we obtain::

    sage: GrowthDiagram(RuleRSK, P)
    0  0  0  0
    0  1  0  1
    0  2  0  1
    1  0  0  0

Alternatively, we could create the same growth diagram using a
matrix.

Let us also mention that one can pass the arguments specifying
a growth diagram directly to the rule::

    sage: RuleRSK(P)
    0  0  0  0
    0  1  0  1
    0  2  0  1
    1  0  0  0

In contrast to the classical insertion algorithms, growth diagrams
immediately generalize to fillings whose shape is an arbitrary skew
partition::

    sage: GrowthDiagram(RuleRSK, [3,1,2], shape=SkewPartition([[3,3,2],[1,1]]))
    .  1  0
    .  0  1
    1  0

As an important example, consider the Stanley-Sundaram correspondence
between oscillating tableaux and (partial) perfect matchings.
Perfect matchings of `\{1, \ldots, 2r\}` are in bijection with
`0-1`-fillings of a triangular shape with `2r-1` rows, such that for
each `k` there is either exactly one non-zero entry in row `k` or
exactly one non-zero entry in column `2r-k`.  Explicitly, if `(i,j)`
is a pair in the perfect matching, the entry in column `i-1` and row
`2r-j` equals `1`.  For example::

    sage: m = [[1,5],[3,4],[2,7],[6,8]]
    sage: G = RuleRSK({(i-1, 8-j): 1 for i,j in m}, shape=[7,6,5,4,3,2,1]); G
    0  0  0  0  0  1  0
    0  1  0  0  0  0
    0  0  0  0  0
    1  0  0  0
    0  0  1
    0  0
    0

The partitions labelling the bottom-right corners along the boundary
opposite of the origin then form a so-called oscillating tableau -
the remaining partitions along the bottom-right boundary are
redundant::

    sage: G.out_labels()[1::2]
    [[1], [1, 1], [2, 1], [1, 1], [1], [1, 1], [1]]

Another great advantage of growth diagrams is that we immediately
have access to a skew version of the correspondence, by providing
different initialization for the labels on the side of the origin.
We reproduce the original example of Bruce Sagan and Richard Stanley,
see also Tom Roby's thesis [Rob1991]_::

    sage: w = {(1-1,4-1): 1, (2-1,2-1): 1, (4-1,3-1): 1}
    sage: T = SkewTableau([[None, None], [None, 5], [1]])
    sage: U = SkewTableau([[None, None], [None, 3], [5]])
    sage: labels = T.to_chain()[::-1] + U.to_chain()[1:]
    sage: G = GrowthDiagram(RuleRSK, filling=w, shape=[5,5,5,5,5], labels=labels); G
    0  0  0  0  0
    0  1  0  0  0
    0  0  0  1  0
    1  0  0  0  0
    0  0  0  0  0
    sage: ascii_art([G.P_symbol(), G.Q_symbol()])
    [   .  .  2  3    .  .  1  4 ]
    [   .  .          .  .       ]
    [   .  4          .  2       ]
    [   1             3          ]
    [   5         ,   5          ]

Similarly, there is a correspondence for skew oscillating tableau.
Let us conclude by reproducing Example 4.2.6 from [Rob1991]_.  The
oscillating tableau, as given, is::

    sage: o = [[2,1],[2,2],[3,2],[4,2],[4,1],[4,1,1],[3,1,1],[3,1],[3,2],[3,1],[2,1]]

From this, we have to construct the list of labels of the corners
along the bottom-right boundary.  The labels with odd indices are
given by the oscillating tableau, the other labels are obtained by
taking the smaller of the two neighbouring partitions::

    sage: l = [o[i//2] if is_even(i) else min(o[(i-1)//2],o[(i+1)//2])
    ....:      for i in range(2*len(o)-1)]
    sage: la = list(range(len(o)-2, 0, -1))
    sage: G = RuleRSK(labels=l[1:-1], shape=la); G
    0  0  0  0  0  0  0  1  0
    0  1  0  0  0  0  0  0
    0  0  0  0  0  0  0
    0  0  0  0  0  0
    0  0  1  0  0
    0  0  0  0
    0  0  0
    0  0
    0

The skew tableaux can now be read off the partitions labelling the
left and the top boundary.  These can be accessed using the method
:meth:`~sage.combinat.growth.GrowthDiagram.in_labels`::

    sage: ascii_art(SkewTableau(chain=G.in_labels()[len(o)-2:]),
    ....:           SkewTableau(chain=G.in_labels()[len(o)-2::-1]))
    .  1  .  7
    5     4

Rules currently available
-------------------------

As mentioned at the beginning, the Robinson-Schensted-Knuth
correspondence is just a special case of growth diagrams.  In
particular, we have implemented the following local rules:

- RSK (:class:`~sage.combinat.growth.RuleRSK`).
- A variation of RSK originally due to Burge
  (:class:`~sage.combinat.growth.RuleBurge`).
- A correspondence producing binary words originally due to Viennot
  (:class:`~sage.combinat.growth.RuleBinaryWord`).
- A correspondence producing domino tableaux
  (:class:`~sage.combinat.growth.RuleDomino`) originally due
  to Barbasch and Vogan.
- A correspondence for shifted shapes
  (:class:`~sage.combinat.growth.RuleShiftedShapes`), where the
  original insertion algorithm is due to Sagan and Worley, and Haiman.
- The Sylvester correspondence, producing binary trees
  (:class:`~sage.combinat.growth.RuleSylvester`).
- The Young-Fibonacci correspondence
  (:class:`~sage.combinat.growth.RuleYoungFibonacci`).
- LLMS insertion (:class:`~sage.combinat.growth.RuleLLMS`).

Background
----------

At the heart of Fomin's framework is the notion of dual graded
graphs.  This is a pair of digraphs `P, Q` (multiple edges being
allowed) on the same set of vertices `V`, that satisfy the following
conditions:

* the graphs are graded, that is, there is a function `\rho : V \to
  \NN`, such that for any edge `(v, w)` of `P` and also of `Q` we
  have `\rho(w) = \rho(v) + 1`,

* there is a vertex `0` with rank zero, and

* there is a positive integer `r` such that `DU = UD + rI` on the
  free `\ZZ`-module `\ZZ[V]`, where `D` is the down operator of `Q`,
  assigning to each vertex the formal sum of its predecessors, `U` is
  the up operator of `P`, assigning to each vertex the formal sum of
  its successors, and `I` is the identity operator.

Note that the condition `DU = UD + rI` is symmetric with respect to
the interchange of the graphs `P` and `Q`, because the up operator of
a graph is the transpose of its down operator.

For example, taking for both `P` and `Q` to be Young's lattice and
`r=1`, we obtain the dual graded graphs for classical Schensted
insertion.

Given such a pair of graphs, there is a bijection between the
`r`-colored permutations on `k` letters and pairs `(p, q)`, where `p`
is a path in `P` from zero to a vertex of rank `k` and `q` is a path
in `Q` from zero to the same vertex.

It turns out that - in principle - this bijection can always be
described by so-called local forward and backward rules, see
[Fom1995]_ for a detailed description.  Knowing at least the forward
rules, or the backward rules, you can implement your own growth
diagram class.

Implementing your own growth diagrams
-------------------------------------

The class :class:`GrowthDiagram` is written so that it is easy to
implement growth diagrams you come across in your research.
Moreover, the class tolerates some deviations from Fomin's
definitions.  For example, although the general
Robinson-Schensted-Knuth correspondence between integer matrices and
semistandard tableaux is, strictly speaking, not a growth on dual
graded graphs, it is supported by our framework.

For illustration, let us implement a growth diagram class with the
backward rule only.  Suppose that the vertices of the graph are the
non-negative integers, the rank is given by the integer itself, and
the backward rule is `(y, z, x) \mapsto (\min(x,y), 0)` if `y = z`
or `x = z` and `(y, z, x) \mapsto (\min(x,y), 1)` otherwise.

We first need to import the base class for a rule::

    sage: from sage.combinat.growth import Rule

Next, we implement the backward rule and the rank function and
provide the bottom element ``zero`` of the graph.  For more
information, see :class:`~sage.combinat.growth.Rule`. ::

    sage: class RulePascal(Rule):
    ....:     zero = 0
    ....:     def rank(self, v): return v
    ....:     def backward_rule(self, y, z, x):
    ....:         return (min(x,y), 0 if y==z or x==z else 1)

We can now compute the filling corresponding to a sequence of labels
as follows::

    sage: GrowthDiagram(RulePascal(), labels=[0,1,2,1,2,1,0])
    1  0  0
    0  0  1
    0  1

Of course, since we have not provided the forward rule, we cannot
compute the labels belonging to a filling::

    sage: GrowthDiagram(RulePascal(), [3,1,2])
    Traceback (most recent call last):
    ...
    AttributeError: 'RulePascal' object has no attribute 'forward_rule'

We now re-implement the rule where we provide the dual graded graphs::

    sage: class RulePascal(Rule):
    ....:     zero = 0
    ....:     def rank(self, v): return v
    ....:     def backward_rule(self, y, z, x):
    ....:         return (min(x,y), 0 if y==z or x==z else 1)
    ....:     def vertices(self, n): return [n]
    ....:     def is_P_edge(self, v, w): return w == v + 1
    ....:     def is_Q_edge(self, v, w): return w == v + 1

Are they really dual? ::

    sage: RulePascal()._check_duality(3)
    Traceback (most recent call last):
    ...
    ValueError: D U - U D differs from 1 I for vertex 3:
    D U = [3]
    U D + 1 I = [3, 3]

With our current definition, duality fails - in fact, there are no
dual graded graphs on the integers without multiple edges.
Consequently, also the backward rule cannot work as ``backward_rule``
requires additional information (the edge labels as arguments).

Let us thus continue with the example from Section 4.7 of [Fom1995]_
instead, which defines dual graded graphs with multiple edges on the
integers.  The color ``self.zero_edge``, which defaults to ``0`` is
reserved for degenerate edges, but may be abused for the unique edge
if one of the graphs has no multiple edges.  For greater clarity in
this example we set it to ``None``::

    sage: class RulePascal(Rule):
    ....:     zero = 0
    ....:     has_multiple_edges = True
    ....:     zero_edge = None
    ....:     def rank(self, v): return v
    ....:     def vertices(self, n): return [n]
    ....:     def is_P_edge(self, v, w): return [0] if w == v + 1 else []
    ....:     def is_Q_edge(self, v, w): return list(range(w)) if w == v+1 else []

We verify these are `1` dual at level `5`::

    sage: RulePascal()._check_duality(5)

Finally, let us provide the backward rule.  The arguments of the rule
are vertices together with the edge labels now, specifying the path
from the lower left to the upper right of the cell.  The horizontal
edges come from `Q`, whereas the vertical edges come from `P`.

Thus, the definition in Section 4.7 of [Fom1995]_ translates as
follows::

    sage: class RulePascal(Rule):
    ....:     zero = 0
    ....:     has_multiple_edges = True
    ....:     zero_edge = None
    ....:     def rank(self, v): return v
    ....:     def vertices(self, n): return [n]
    ....:     def is_P_edge(self, v, w): return [0] if w == v + 1 else []
    ....:     def is_Q_edge(self, v, w): return list(range(w)) if w == v+1 else []
    ....:     def backward_rule(self, y, g, z, h, x):
    ....:         if g is None:
    ....:             return (0, x, None, 0)
    ....:         if h is None:
    ....:             return (None, y, g, 0)
    ....:         if g == 0:
    ....:             return (None, y, None, 1)
    ....:         else:
    ....:             return (0, x-1, g-1, 0)

The labels are now alternating between vertices and edge-colors::

    sage: GrowthDiagram(RulePascal(), labels=[0,0,1,0,2,0,1,0,0])
    1  0
    0  1

    sage: GrowthDiagram(RulePascal(), labels=[0,0,1,1,2,0,1,0,0])
    0  1
    1  0
"""

# ****************************************************************************
#       Copyright (C) 2017 Martin Rubey <martin.rubey at tuwien.ac.at>
#                     2017 Travis Scrimshaw <tcscrims at gmail.com>
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
# ***************************************************************************
from __future__ import annotations
from copy import copy
from itertools import zip_longest

from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation
from sage.combinat.posets.posets import Poset
from sage.combinat.words.word import Word
from sage.combinat.words.words import Words
from sage.combinat.binary_tree import BinaryTree, BinaryTrees, LabelledBinaryTree
from sage.combinat.composition import Compositions
from sage.combinat.partition import _Partitions, Partitions
from sage.combinat.skew_partition import SkewPartition
from sage.combinat.skew_tableau import SkewTableau
from sage.combinat.core import Core, Cores
from sage.combinat.k_tableau import WeakTableau, StrongTableau
from sage.combinat.shifted_primed_tableau import ShiftedPrimedTableau
from sage.graphs.digraph import DiGraph


def _make_partition(l):
    """
    Return the list as a partition.

    This is intended to be fast, so checks are bypassed.

    TESTS::

        sage: from sage.combinat.growth import _make_partition
        sage: p = _make_partition([3,2,1,0]); p
        [3, 2, 1]

        sage: p.parent()
        Partitions
    """
    return _Partitions.element_class(_Partitions, l)


class GrowthDiagram(SageObject):
    r"""
    A generalized Schensted growth diagram in the sense of Fomin.

    Growth diagrams were introduced by Sergey Fomin [Fom1994]_,
    [Fom1995]_ and provide a vast generalization of the
    Robinson-Schensted-Knuth (RSK) correspondence between matrices
    with non-negative integer entries and pairs of semistandard Young
    tableaux of the same shape.

    A growth diagram is based on the notion of *dual graded graphs*,
    a pair of digraphs `P, Q` (multiple edges being allowed) on the
    same set of vertices `V`, that satisfy the following conditions:

    * the graphs are graded, that is, there is a function `\rho:
      V \to \NN`, such that for any edge `(v, w)` of `P` and also
      of `Q` we have `\rho(w) = \rho(v) + 1`,

    * there is a vertex `0` with rank zero, and

    * there is a positive integer `r` such that `DU = UD + rI` on the
      free `\ZZ`-module `\ZZ[V]`, where `D` is the down operator of
      `Q`, assigning to each vertex the formal sum of its
      predecessors, `U` is the up operator of `P`, assigning to each
      vertex the formal sum of its successors, and `I` is the
      identity operator.

    Growth diagrams are defined by providing a pair of local rules: a
    'forward' rule, whose input are three vertices `y`, `t` and `x`
    of the dual graded graphs and an integer, and whose output is a
    fourth vertex `z`.  This rule should be invertible in the
    following sense: there is a so-called 'backward' rule that
    recovers the integer and `t` given `y`, `z` and `x`.

    All implemented growth diagram rules are available by
    ``GrowthDiagram.rules.<tab>``. The current list is:

    - :class:`~sage.combinat.growth.RuleRSK` -- RSK
    - :class:`~sage.combinat.growth.RuleBurge` -- a variation of RSK
      originally due to Burge
    - :class:`~sage.combinat.growth.RuleBinaryWord` -- a correspondence
      producing binary words originally due to Viennot
    - :class:`~sage.combinat.growth.RuleDomino` -- a correspondence
      producing domino tableaux originally due to Barbasch and Vogan
    - :class:`~sage.combinat.growth.RuleShiftedShapes` -- a correspondence
      for shifted shapes, where the original insertion algorithm is due
      to Sagan and Worley, and Haiman.
    - :class:`~sage.combinat.growth.RuleSylvester` -- the Sylvester
      correspondence, producing binary trees
    - :class:`~sage.combinat.growth.RuleYoungFibonacci` -- the
      Young-Fibonacci correspondence
    - :class:`~sage.combinat.growth.RuleLLMS` -- LLMS insertion

    INPUT:

    - ``rule`` -- :class:`~sage.combinat.growth.Rule`;
      the growth diagram rule

    - ``filling`` -- (optional) a dictionary whose keys are coordinates
      and values are integers, a list of lists of integers, or a word
      with integer values; if a word, then negative letters but without
      repetitions are allowed and interpreted as coloured permutations

    - ``shape`` -- (optional) a (possibly skew) partition

    - ``labels`` -- (optional) a list that specifies a path whose length
      in the half-perimeter of the shape; more details given below

    If ``filling`` is not given, then the growth diagram is determined
    by applying the backward rule to ``labels`` decorating the
    boundary opposite of the origin of the ``shape``. In this case,
    ``labels`` are interpreted as labelling the boundary opposite of
    the origin.

    Otherwise, ``shape`` is inferred from ``filling`` or ``labels`` if
    possible and ``labels`` is set to ``rule.zero`` if not specified.
    Here, ``labels`` are labelling the boundary on the side of the origin.

    For ``labels``, if ``rule.has_multiple_edges`` is ``True``, then the
    elements should be of the form `(v_1, e_1, \ldots, e_{n-1}, v_n)`,
    where `n` is the half-perimeter of ``shape``, and `(v_{i-1}, e_i, v_i)`
    is an edge in the dual graded graph for all `i`.  Otherwise, it is a
    list of `n` vertices.

    .. NOTE::

        Coordinates are of the form ``(col, row)`` where the origin is
        in the upper left, to be consistent with permutation matrices
        and skew tableaux (in English convention).  This is different
        from Fomin's convention, who uses a Cartesian coordinate system.

        Conventions are chosen such that for permutations, the same
        growth diagram is constructed when passing the permutation
        matrix instead.

    EXAMPLES:

    We create a growth diagram using the forward RSK rule and a permutation::

        sage: RuleRSK = GrowthDiagram.rules.RSK()
        sage: pi = Permutation([4, 1, 2, 3])
        sage: G = GrowthDiagram(RuleRSK, pi); G
        0  1  0  0
        0  0  1  0
        0  0  0  1
        1  0  0  0
        sage: G.out_labels()
        [[], [1], [1, 1], [2, 1], [3, 1], [3], [2], [1], []]

    Passing the permutation matrix instead gives the same result::

        sage: G = GrowthDiagram(RuleRSK, pi.to_matrix())
        sage: ascii_art([G.P_symbol(), G.Q_symbol()])
        [   1  2  3    1  3  4 ]
        [   4      ,   2       ]

    We give the same example but using a skew shape::

        sage: shape = SkewPartition([[4,4,4,2],[1,1]])
        sage: G = GrowthDiagram(RuleRSK, pi, shape=shape); G
        .  1  0  0
        .  0  1  0
        0  0  0  1
        1  0
        sage: G.out_labels()
        [[], [1], [1, 1], [1], [2], [3], [2], [1], []]

    We construct a growth diagram using the backwards RSK rule by
    specifying the labels::

        sage: GrowthDiagram(RuleRSK, labels=G.out_labels())
        0  1  0  0
        0  0  1  0
        0  0  0  1
        1  0
    """
    def __init__(self, rule, filling=None, shape=None, labels=None):
        r"""
        Initialize ``self``.

        TESTS::

            sage: RuleRSK = GrowthDiagram.rules.RSK()
            sage: w = [3,3,2,4,1]; G = GrowthDiagram(RuleRSK, w)
            sage: [G.P_symbol(), G.Q_symbol()]
            [[[1, 3, 4], [2], [3]], [[1, 2, 4], [3], [5]]]
            sage: RSK(w)
            [[[1, 3, 4], [2], [3]], [[1, 2, 4], [3], [5]]]

            sage: TestSuite(G).run()

            sage: GrowthDiagram(RuleRSK)
            Traceback (most recent call last):
            ...
            ValueError: please provide a filling or a sequence of labels
        """
        if not isinstance(rule, Rule):
            raise TypeError("the rule must be an instance of Rule")
        self.rule = rule

        if filling is None:
            if labels is None:
                raise ValueError("please provide a filling or a sequence of labels")

            labels = self._process_labels(labels)

            if shape is None:
                shape = self._shape_from_labels(labels)

            self._lambda, self._mu = self._process_shape(shape)
            self._out_labels = labels
            self._check_labels(self._out_labels)
            self._shrink()
        else:
            self._filling, (self._lambda, self._mu) = self._process_filling_and_shape(filling, shape)

            if labels is None:
                rule = self.rule
                if rule.has_multiple_edges:
                    self._in_labels = [rule.zero, rule.zero_edge]*(self.half_perimeter()-1) + [rule.zero]
                else:
                    self._in_labels = [rule.zero] * self.half_perimeter()
            else:
                labels = self._process_labels(labels)
                self._in_labels = labels

            self._check_labels(self._in_labels)
            self._grow()

    def filling(self):
        r"""
        Return the filling of the diagram as a dictionary.

        EXAMPLES::

            sage: RuleRSK = GrowthDiagram.rules.RSK()
            sage: G = GrowthDiagram(RuleRSK, [[0,1,0], [1,0,2]])
            sage: G.filling()
            {(0, 1): 1, (1, 0): 1, (2, 1): 2}
        """
        return self._filling

    def conjugate(self):
        r"""
        Return the conjugate growth diagram of ``self``.

        This is the growth diagram with the filling reflected over the
        main diagonal.

        The sequence of labels along the boundary on the side of the
        origin is the reversal of the corresponding sequence of the
        original growth diagram.

        When the filling is a permutation, the conjugate filling
        corresponds to its inverse.

        EXAMPLES::

            sage: RuleRSK = GrowthDiagram.rules.RSK()
            sage: G = GrowthDiagram(RuleRSK, [[0,1,0], [1,0,2]])
            sage: Gc = G.conjugate()
            sage: (Gc.P_symbol(), Gc.Q_symbol()) == (G.Q_symbol(), G.P_symbol())
            True

        TESTS:

        Check that labels and shape are handled correctly::

            sage: o = [[2,1],[2,2],[3,2],[4,2],[4,1],[4,1,1],[3,1,1],[3,1],[3,2],[3,1],[2,1]]
            sage: l = [o[i//2] if is_even(i) else min(o[(i-1)//2],o[(i+1)//2])
            ....:      for i in range(2*len(o)-1)]
            sage: la = list(range(len(o)-2, 0, -1))
            sage: G = RuleRSK(labels=l[1:-1], shape=la)
            sage: G.out_labels() == G.conjugate().out_labels()[::-1]
            True
        """
        F = {(j,i): v for (i,j),v in self._filling.items()}
        return GrowthDiagram(self.rule,
                             filling=F,
                             shape=self.shape().conjugate(),
                             labels=self.in_labels()[::-1])

    def rotate(self):
        r"""
        Return the growth diagram with the filling rotated by 180 degrees.

        The rotated growth diagram is initialized with
        ``labels=None``, that is, all labels along the boundary on
        the side of the origin are set to ``rule.zero``.

        For RSK-growth diagrams and rectangular fillings, this
        corresponds to evacuation of the `P`- and the `Q`-symbol.

        EXAMPLES::

            sage: RuleRSK = GrowthDiagram.rules.RSK()
            sage: G = GrowthDiagram(RuleRSK, [[0,1,0], [1,0,2]])
            sage: Gc = G.rotate()
            sage: ascii_art([Gc.P_symbol(), Gc.Q_symbol()])
            [   1  1  1    1  1  2 ]
            [   2      ,   3       ]

            sage: ascii_art([Tableau(t).evacuation()
            ....:            for t in [G.P_symbol(), G.Q_symbol()]])
            [   1  1  1    1  1  2 ]
            [   2      ,   3       ]

        TESTS:

        Check that shape is handled correctly::

            sage: RuleRSK = GrowthDiagram.rules.RSK()
            sage: G = GrowthDiagram(RuleRSK,
            ....:                   filling={(0,2):1, (3,1):2, (2,1):3},
            ....:                   shape=SkewPartition([[5,5,5,3],[3,1]]))
            sage: G
            .  .  .  0  0
            .  0  3  2  0
            1  0  0  0  0
            0  0  0
            sage: G.rotate()
            .  .  0  0  0
            0  0  0  0  1
            0  2  3  0
            0  0
        """
        l = self._lambda[0]
        h = len(self._lambda)
        shape_lambda = [l-p for p in self._mu] + [l]*(h-len(self._mu))
        shape_mu     = [l-p for p in self._lambda]
        shape = SkewPartition([shape_lambda[::-1], shape_mu[::-1]])
        F = {(l-i-1, h-j-1): v for (i,j),v in self._filling.items()}
        return GrowthDiagram(self.rule,
                             filling=F,
                             shape=shape)

    def half_perimeter(self):
        r"""
        Return half the perimeter of the shape of the growth diagram.

        TESTS::

            sage: RuleRSK = GrowthDiagram.rules.RSK()
            sage: G = GrowthDiagram(RuleRSK, {(0,1):1, (2,0):1}, SkewPartition([[3,1],[1]])); G
            .  0  1
            1
            sage: G.half_perimeter()
            6
        """
        # Assume that ``self._lambda`` is already set.
        if not self._lambda:
            return 1
        return self._lambda[0] + len(self._lambda) + 1

    def shape(self):
        r"""
        Return the shape of the growth diagram as a skew partition.

        .. WARNING::

            In the literature the label on the corner opposite of the
            origin of a rectangular filling is often called the shape
            of the filling.  This method returns the shape of the
            region instead.

        EXAMPLES::

            sage: RuleRSK = GrowthDiagram.rules.RSK()
            sage: GrowthDiagram(RuleRSK, [1]).shape()
            [1] / []
        """
        return SkewPartition([self._lambda, self._mu])

    def out_labels(self):
        r"""
        Return the labels along the boundary opposite of the origin.

        EXAMPLES::

            sage: RuleRSK = GrowthDiagram.rules.RSK()
            sage: G = GrowthDiagram(RuleRSK, [[0,1,0], [1,0,2]])
            sage: G.out_labels()
            [[], [1], [1, 1], [3, 1], [1], []]
        """
        return self._out_labels

    def in_labels(self):
        r"""
        Return the labels along the boundary on the side of the origin.

        EXAMPLES::

            sage: RuleRSK = GrowthDiagram.rules.RSK()
            sage: G = GrowthDiagram(RuleRSK, labels=[[2,2],[3,2],[3,3],[3,2]]); G
            1 0
            sage: G.in_labels()
            [[2, 2], [2, 2], [2, 2], [3, 2]]
        """
        return self._in_labels

    def P_symbol(self):
        r"""
        Return the labels along the vertical boundary of a rectangular
        growth diagram as a generalized standard tableau.

        EXAMPLES::

            sage: RuleRSK = GrowthDiagram.rules.RSK()
            sage: G = GrowthDiagram(RuleRSK, [[0,1,0], [1,0,2]])
            sage: ascii_art([G.P_symbol(), G.Q_symbol()])
            [   1  2  2    1  3  3 ]
            [   2      ,   2       ]
        """
        return self.rule.P_symbol(self.P_chain())

    def Q_symbol(self):
        r"""
        Return the labels along the horizontal boundary of a rectangular
        growth diagram as a generalized standard tableau.

        EXAMPLES::

            sage: RuleRSK = GrowthDiagram.rules.RSK()
            sage: G = GrowthDiagram(RuleRSK, [[0,1,0], [1,0,2]])
            sage: ascii_art([G.P_symbol(), G.Q_symbol()])
            [   1  2  2    1  3  3 ]
            [   2      ,   2       ]
        """
        return self.rule.Q_symbol(self.Q_chain())

    def P_chain(self):
        r"""
        Return the labels along the vertical boundary of a rectangular
        growth diagram.

        EXAMPLES::

            sage: BinaryWord = GrowthDiagram.rules.BinaryWord()
            sage: G = GrowthDiagram(BinaryWord, [4, 1, 2, 3])
            sage: G.P_chain()
            [word: , word: 1, word: 11, word: 111, word: 1011]

        Check that :trac:`25631` is fixed::

            sage: BinaryWord = GrowthDiagram.rules.BinaryWord()
            sage: BinaryWord(filling = {}).P_chain()
            [word: ]

        """
        if not self.is_rectangular():
            raise ValueError("the P symbol is only defined for rectangular shapes")
        if self._lambda:
            if self.rule.has_multiple_edges:
                r = 2*self._lambda[0]
            else:
                r = self._lambda[0]
        else:
            r = 0
        return self._out_labels[r:][::-1]

    def Q_chain(self):
        r"""
        Return the labels along the horizontal boundary of a rectangular
        growth diagram.

        EXAMPLES::

            sage: BinaryWord = GrowthDiagram.rules.BinaryWord()
            sage: G = GrowthDiagram(BinaryWord, [[0,1,0,0], [0,0,1,0], [0,0,0,1], [1,0,0,0]])
            sage: G.Q_chain()
            [word: , word: 1, word: 10, word: 101, word: 1011]

        Check that :trac:`25631` is fixed::

            sage: BinaryWord = GrowthDiagram.rules.BinaryWord()
            sage: BinaryWord(filling = {}).Q_chain()
            [word: ]

        """
        if not self.is_rectangular():
            raise ValueError("the Q symbol is only defined for rectangular shapes")
        if self._lambda:
            if self.rule.has_multiple_edges:
                r = 2*self._lambda[0]+1
            else:
                r = self._lambda[0]+1
        else:
            r = 1
        return self._out_labels[:r]

    def is_rectangular(self):
        r"""
        Return ``True`` if the shape of the growth diagram is rectangular.

        EXAMPLES::

            sage: RuleRSK = GrowthDiagram.rules.RSK()
            sage: GrowthDiagram(RuleRSK, [2,3,1]).is_rectangular()
            True
            sage: GrowthDiagram(RuleRSK, [[1,0,1],[0,1]]).is_rectangular()
            False
        """
        return (all(x == 0 for x in self._mu)
                and all(x == self._lambda[0] for x in self._lambda))

    def to_word(self):
        r"""
        Return the filling as a word, if the shape is rectangular and
        there is at most one nonzero entry in each column, which must
        be 1.

        EXAMPLES::

            sage: RuleRSK = GrowthDiagram.rules.RSK()
            sage: w = [3,3,2,4,1]; G = GrowthDiagram(RuleRSK, w)
            sage: G
            0  0  0  0  1
            0  0  1  0  0
            1  1  0  0  0
            0  0  0  1  0
            sage: G.to_word()
            [3, 3, 2, 4, 1]
        """
        if not self.is_rectangular():
            raise ValueError("can only convert fillings of rectangular shapes to words")
        w = [0] * self._lambda[0]
        for ((i,j), v) in self._filling.items():
            if v != 0:
                if v == 1:
                    if w[i] == 0:
                        w[i] = j+1
                    else:
                        raise ValueError("can only convert fillings with at"
                                         " most one entry per column to words")
                elif v == -1:
                    if w[i] == 0:
                        w[i] = -(j+1)
                    else:
                        raise ValueError("can only convert fillings with at"
                                         " most one entry per column to words")
                else:
                    raise ValueError("can only convert 0-1 fillings to words;"
                                     " try 'to_biword'")
        return w

    def to_biword(self):
        r"""
        Return the filling as a biword, if the shape is rectangular.

        EXAMPLES::

            sage: RuleRSK = GrowthDiagram.rules.RSK()
            sage: P = Tableau([[1,2,2],[2]])
            sage: Q = Tableau([[1,3,3],[2]])
            sage: bw = RSK_inverse(P, Q); bw
            [[1, 2, 3, 3], [2, 1, 2, 2]]
            sage: G = GrowthDiagram(RuleRSK, labels=Q.to_chain()[:-1]+P.to_chain()[::-1]); G
            0  1  0
            1  0  2

            sage: P = SemistandardTableau([[1, 1, 2], [2]])
            sage: Q = SemistandardTableau([[1, 2, 2], [2]])
            sage: G = GrowthDiagram(RuleRSK, labels=Q.to_chain()[:-1]+P.to_chain()[::-1]); G
            0  2
            1  1
            sage: G.to_biword()
            ([1, 2, 2, 2], [2, 1, 1, 2])
            sage: RSK([1, 2, 2, 2], [2, 1, 1, 2])
            [[[1, 1, 2], [2]], [[1, 2, 2], [2]]]
        """
        if not self.is_rectangular():
            raise ValueError("can only convert fillings of rectangular shapes to words")
        w1 = []
        w2 = []
        for ((i,j), v) in sorted(self._filling.items()):
            if v >= 0:
                w1.extend([i+1]*v)
                w2.extend([j+1]*v)
            else:
                raise ValueError("can only convert fillings with"
                                 " non-negative entries to words")
        return (w1, w2)

    def __iter__(self):
        r"""
        Return the rows of the filling.

        TESTS::

            sage: RuleRSK = GrowthDiagram.rules.RSK()
            sage: G = GrowthDiagram(RuleRSK, {(0,1):1, (1,0):1}, SkewPartition([[2,1],[1]]))
            sage: list(G)
            [[None, 1], [1]]

            sage: pi = Permutation([2,3,1,6,4,5])
            sage: G = GrowthDiagram(RuleRSK, pi)
            sage: list(G)
            [[0, 0, 1, 0, 0, 0],
             [1, 0, 0, 0, 0, 0],
             [0, 1, 0, 0, 0, 0],
             [0, 0, 0, 0, 1, 0],
             [0, 0, 0, 0, 0, 1],
             [0, 0, 0, 1, 0, 0]]
        """
        return ([None]*self._mu[r] + [self._filling.get((self._mu[r]+j,r), 0)
                                      for j in range(self._lambda[r]-self._mu[r])]
                for r in range(len(self._lambda)))

    def _repr_(self):
        r"""
        Return a string with the filling of the growth diagram
        as a skew tableau.

        TESTS::

            sage: RuleRSK = GrowthDiagram.rules.RSK()
            sage: GrowthDiagram(RuleRSK, {(0,1):1, (1,0):1}, SkewPartition([[2,1],[1]]))
            .  1
            1

            sage: GrowthDiagram(RuleRSK, {(0,1):1, (2,0):1}, SkewPartition([[3,1],[1]]))
            .  0  1
            1
        """
        return SkewTableau(expr=[self._mu,
                                 [[self._filling.get((self._mu[r]+j,r), 0)
                                   for j in range(self._lambda[r]-self._mu[r])]
                                  for r in range(len(self._lambda))][::-1]])._repr_diagram()

    def __eq__(self, other):
        r"""
        Return ``True`` if the growth diagram ``other`` has the same
        shape and the same filling as ``self``.

        EXAMPLES:

        Equality ignores zeros in fillings::

            sage: RuleRSK = GrowthDiagram.rules.RSK()
            sage: G1 = GrowthDiagram(RuleRSK, {(0, 1): 1, (1, 0): 1})
            sage: G2 = GrowthDiagram(RuleRSK, {(0, 0): 0, (0, 1): 1, (1, 0): 1})
            sage: G1 == G2
            True

        Growth diagrams with different shapes are different::

            sage: G1 = GrowthDiagram(RuleRSK, [[0,1,0],[1,0]])
            sage: G2 = GrowthDiagram(RuleRSK, [[0,1,0],[1]])
            sage: G1 == G2
            False

        Growth diagrams with different rules are different::

            sage: G1 = GrowthDiagram(RuleRSK, {(0, 1): 1, (1, 0): 1})
            sage: BinaryWord = GrowthDiagram.rules.BinaryWord()
            sage: G2 = GrowthDiagram(BinaryWord, {(0, 1): 1, (1, 0): 1})
            sage: G1 == G2
            False
        """
        return (type(self) == type(other) and
                self.rule == other.rule and
                self._lambda == other._lambda and
                self._mu == other._mu and
                self._filling == other._filling)

    def __ne__(self, other):
        r"""
        Return ``True`` if the growth diagram ``other`` does not have the
        same shape and the same filling as ``self``.

        TESTS:

        Equality ignores zeros in fillings::

            sage: RuleRSK = GrowthDiagram.rules.RSK()
            sage: G1 = GrowthDiagram(RuleRSK, {(0, 1): 1, (1, 0): 1})
            sage: G2 = GrowthDiagram(RuleRSK, {(0, 0): 0, (0, 1): 1, (1, 0): 1})
            sage: G1 != G2
            False

        Growth diagrams with different shapes are different::

            sage: G1 = GrowthDiagram(RuleRSK, [[0,1,0],[1,0]])
            sage: G2 = GrowthDiagram(RuleRSK, [[0,1,0],[1]])
            sage: G1 != G2
            True

        Growth diagrams with different rules are different::

            sage: RuleRSK = GrowthDiagram.rules.RSK()
            sage: BinaryWord = GrowthDiagram.rules.BinaryWord()
            sage: G1 = GrowthDiagram(RuleRSK, {(0, 1): 1, (1, 0): 1})
            sage: G2 = GrowthDiagram(BinaryWord, {(0, 1): 1, (1, 0): 1})
            sage: G1 != G2
            True
        """
        return not (self == other)

    def _process_labels(self, labels):
        r"""
        Return the list of labels such that each element has the
        correct type from the rule.

        .. WARNING::

            Assumes that ``self.rule`` is set.

        EXAMPLES::

            sage: RuleRSK = GrowthDiagram.rules.RSK()
            sage: labels = [[], [1], [1,1], [1], []]
            sage: G = GrowthDiagram(RuleRSK, labels=labels)  # indirect doctest
            sage: G.out_labels()[2].parent()
            Partitions
        """
        rule = self.rule
        if rule.has_multiple_edges:
            return [rule.normalize_vertex(val) if i % 2 == 0 else val
                    for i, val in enumerate(labels)]
        else:
            return [rule.normalize_vertex(la) for la in labels]

    def _shape_from_labels(self, labels):
        r"""
        Determine the shape of the growth diagram given a list of labels
        during initialization.

        The shape can be determined from the labels if the size of
        each label differs from the size of its successor.

        Otherwise raise an error.

        .. WARNING::

            Assumes that ``self.rule`` and ``self.rank`` is set.

        TESTS::

            sage: RuleRSK = GrowthDiagram.rules.RSK()
            sage: labels = [[],[2],[1],[],[1],[]]
            sage: G = GrowthDiagram(RuleRSK, labels=labels); G
            0 1
            1
            1
            sage: G._shape_from_labels(G.out_labels())
            [2, 1, 1]

            sage: Shifted = GrowthDiagram.rules.ShiftedShapes()
            sage: Shifted({(0, 0): 1}).out_labels()
            [[], 1, [1], 0, []]
            sage: Shifted(labels=[[], 1, [2], 0, []])  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: [] has smaller rank than [2] but there is no edge of color 1 in Q
        """
        # we can determine the shape even if is_P_edge is not implemented
        rule = self.rule
        is_P_edge = getattr(rule, "is_P_edge", None)
        is_Q_edge = getattr(rule, "is_Q_edge", None)
        if rule.has_multiple_edges:
            def right_left_multi(la, mu, e) -> int:
                if rule.rank(la) < rule.rank(mu):
                    if is_Q_edge is not None and e not in is_Q_edge(la, mu):
                        raise ValueError("%s has smaller rank than %s but there is no edge of color %s in Q" % (la, mu, e))
                    return 1
                elif rule.rank(la) > rule.rank(mu):
                    if is_P_edge is not None and e not in is_P_edge(mu, la):
                        raise ValueError("%s has smaller rank than %s but there is no edge of color %s in P" % (mu, la, e))
                    return 0
                raise ValueError("can only determine the shape of the growth"
                                 " diagram if ranks of successive labels differ")
            return _Partitions.from_zero_one([right_left_multi(labels[i], labels[i+2], labels[i+1])
                                              for i in range(0, len(labels)-2, 2)])
        else:
            def right_left(la, mu) -> int:
                if rule.rank(la) < rule.rank(mu):
                    if is_Q_edge is not None and not is_Q_edge(la, mu):
                        raise ValueError("%s has smaller rank than %s but is not covered by it in Q" % (la, mu))
                    return 1
                elif rule.rank(la) > rule.rank(mu):
                    if is_P_edge is not None and not is_P_edge(mu, la):
                        raise ValueError("%s has smaller rank than %s but is not covered by it in P" % (mu, la))
                    return 0
                raise ValueError("can only determine the shape of the growth"
                                 " diagram if ranks of successive labels differ")
            return _Partitions.from_zero_one([right_left(labels[i], labels[i+1])
                                              for i in range(len(labels)-1)])

    def _check_labels(self, labels):
        r"""
        Check sanity of the parameter ``labels``.

        .. WARNING::

            Assumes that ``self.rule`` and ``self._lambda`` is set.

        TESTS::

            sage: RuleRSK = GrowthDiagram.rules.RSK()
            sage: GrowthDiagram(RuleRSK, shape=[1], labels=[[], [1]])  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: the number of labels is 2, but for this shape we need 3

            sage: GrowthDiagram(RuleRSK, labels=[[], [1], [2], [2,1]])  # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: the number of labels is 4, but for this shape we need 1

        .. TODO::

            Can we do something more sensible when the chain
            of labels is strictly increasing?
        """
        half_perimeter = self.half_perimeter()
        if self.rule.has_multiple_edges:
            if not (len(labels) % 2):
                raise ValueError("only a list of odd length can specify a path, but %s has even length" % len(labels))
            path_length = (len(labels) + 1) / 2
        else:
            path_length = len(labels)

        if path_length != half_perimeter:
            raise ValueError("the number of labels is %s, but for this shape we need %s"
                             % (path_length, half_perimeter))

    def _process_shape(self, shape):
        r"""
        Return a pair of partitions as lists describing the region
        of the growth diagram.

        TESTS:

        ``shape`` is a skew partition::

            sage: RuleRSK = GrowthDiagram.rules.RSK()
            sage: filling = []
            sage: shape = SkewPartition([[4,2,1,1],[2,1,1]])
            sage: G = GrowthDiagram(RuleRSK, filling, shape)  # indirect doctest
            sage: G._lambda, G._mu
            ([4, 2, 1, 1], [2, 1, 1, 0])

        ``shape`` is a partition::

            sage: filling = []
            sage: shape = Partition([3,2,1,1])
            sage: G = GrowthDiagram(RuleRSK, filling, shape)  # indirect doctest
            sage: G._lambda, G._mu
            ([3, 2, 1, 1], [0, 0, 0, 0])
        """
        try:
            shape = _Partitions(shape)
        except ValueError:
            try:
                shape = SkewPartition(shape)
            except ValueError:
                raise ValueError("cannot make sense of shape %s" % shape)
            return ( list(shape[0]),
                     list(shape[1]) + [0]*(len(shape[0])-len(shape[1])) )
        return (list(shape), [0]*len(shape))

    def _process_filling_and_shape(self, filling, shape):
        r"""
        Return a dict ``F`` such that ``F[(i,j)]`` is the element in row
        ``i`` and column ``j`` and a pair of partitions describing the
        region of the growth diagram.

        TESTS:

        ``filling`` is a dict of coordinates::

            sage: RuleRSK = GrowthDiagram.rules.RSK()
            sage: pi = Permutation([2,3,1,6,4,5])
            sage: G = GrowthDiagram(RuleRSK, {(i,pi[i]-1):1 for i in range(len(pi))})  # indirect doctest
            sage: G._filling
            {(0, 1): 1, (1, 2): 1, (2, 0): 1, (3, 5): 1, (4, 3): 1, (5, 4): 1}
            sage: G.shape()
            [6, 6, 6, 6, 6, 6] / []

        ``filling`` is a dict of dicts::

            sage: G = GrowthDiagram(RuleRSK, {i:{pi[i]-1:1} for i in range(len(pi))})  # indirect doctest
            sage: G._filling
            {(0, 1): 1, (1, 2): 1, (2, 0): 1, (3, 5): 1, (4, 3): 1, (5, 4): 1}
            sage: G.shape()
            [6, 6, 6, 6, 6, 6] / []

        ``filling`` is a matrix::

            sage: G = GrowthDiagram(RuleRSK, pi.to_matrix())  # indirect doctest
            sage: G._filling
            {(0, 1): 1, (1, 2): 1, (2, 0): 1, (3, 5): 1, (4, 3): 1, (5, 4): 1}
            sage: G.shape()
            [6, 6, 6, 6, 6, 6] / []

        ``filling`` is a permutation::

            sage: G = GrowthDiagram(RuleRSK, pi)  # indirect doctest
            sage: G._filling
            {(0, 1): 1, (1, 2): 1, (2, 0): 1, (3, 5): 1, (4, 3): 1, (5, 4): 1}
            sage: G.shape()
            [6, 6, 6, 6, 6, 6] / []

        ``filling`` is a list::

            sage: G = GrowthDiagram(RuleRSK, [3,1,4,1,5])  # indirect doctest
            sage: G._filling
            {(0, 2): 1, (1, 0): 1, (2, 3): 1, (3, 0): 1, (4, 4): 1}
            sage: G.shape()
            [5, 5, 5, 5, 5] / []

        ``filling`` is a list of lists::

            sage: G = GrowthDiagram(RuleRSK, [[1,0,1],[0,1]])  # indirect doctest
            sage: G._filling
            {(0, 0): 1, (1, 1): 1, (2, 0): 1}
            sage: G.shape()
            [3, 2] / []

        ``filling`` is a list of lists and shape is given::

            sage: G = GrowthDiagram(RuleRSK, [[1,0,1],[0,1]], shape=SkewPartition([[3,2],[1]]))  # indirect doctest
            sage: G._filling
            {(0, 0): 1, (1, 1): 1, (2, 0): 1}
            sage: G.shape()
            [3, 2] / [1]

        ``filling`` is empty and shape is ``None``::

            sage: G = GrowthDiagram(RuleRSK, {})
            sage: (G.filling(), G.shape())
            ({}, [] / [])
        """
        if isinstance(filling, dict):
            try:
                v = next(iter(filling.values()))
                if isinstance(v, dict):
                    # it is a dict of dicts
                    F = dict()
                    for (i, row) in filling.items():
                        for (j, v) in row.items():
                            if v != 0:
                                F[(i,j)] = int(v)
                else:
                    # it is dict of coordinates
                    F = {(i,j): v for ((i,j), v) in filling.items()
                         if v != 0}
            except StopIteration:
                # it is an empty dict of coordinates
                F = filling

        else:
            # it is a sequence
            F = dict()
            try:
                # it is a sequence of sequences
                for i, row in enumerate(filling):
                    for j, v in enumerate(row):
                        if v != 0:
                            F[j,i] = int(v)
                if shape is None:
                    shape = [len(row) for row in filling]

            except TypeError:
                # it is a word - for convenience we allow signed words
                for i, l in enumerate(filling):
                    if l > 0:
                        F[i, l-1] = 1
                    else:
                        F[i, -l-1] = -1

        if shape is None:
            if F == {}:
                shape = []
            else:
                # find bounding rectangle of ``filling``
                max_row = max(i for i, _ in F)+1
                max_col = max(j for _, j in F)+1
                shape = [max_row] * max_col

        return (F, self._process_shape(shape))

    def _grow(self):
        r"""
        Compute the labels on the boundary opposite of the origin, given
        the filling.

        TESTS::

            sage: RuleRSK = GrowthDiagram.rules.RSK()
            sage: pi = Permutation([1])
            sage: G = GrowthDiagram(RuleRSK, pi)  # indirect doctest
            sage: G._out_labels
            [[], [1], []]

            sage: pi = Permutation([1,2])
            sage: G = GrowthDiagram(RuleRSK, pi)  # indirect doctest
            sage: G._out_labels
            [[], [1], [2], [1], []]

            sage: pi = Permutation([2,1])
            sage: G = GrowthDiagram(RuleRSK, pi)  # indirect doctest
            sage: G._out_labels
            [[], [1], [1, 1], [1], []]

            sage: G = GrowthDiagram(RuleRSK, {(0,1):1, (1,0):1}, SkewPartition([[2,1],[1]]))  # indirect doctest
            sage: G._out_labels
            [[], [1], [], [1], []]

            sage: G = GrowthDiagram(RuleRSK, {(1,1):1}, SkewPartition([[2,2],[1]]), labels=[[],[],[1],[],[]])  # indirect doctest
            sage: G._out_labels
            [[], [1], [2], [1], []]

            sage: BinaryWord = GrowthDiagram.rules.BinaryWord()
            sage: G = GrowthDiagram(BinaryWord, {(1,1):1}, SkewPartition([[2,2],[1]]), labels=[[],[],[1],[],[]])  # indirect doctest
            sage: G._out_labels
            [word: , word: 1, word: 11, word: 1, word: ]
        """
        labels = list(self._in_labels)
        l = len(self._lambda)
        rule = self.rule
        if rule.has_multiple_edges:
            for r in range(l):
                for c in range(self._mu[r]+l-r, self._lambda[r]+l-r):
                    j = r
                    i = c-l+r
                    (labels[2*c-1],
                     labels[2*c],
                     labels[2*c+1]) = rule.forward_rule(labels[2*c-2],
                                                        labels[2*c-1],
                                                        labels[2*c],
                                                        labels[2*c+1],
                                                        labels[2*c+2],
                                                        self._filling.get((i,j), 0))
        else:
            for r in range(l):
                for c in range(self._mu[r]+l-r, self._lambda[r]+l-r):
                    j = r
                    i = c-l+r
                    labels[c] = rule.forward_rule(labels[c-1],
                                                  labels[c],
                                                  labels[c+1],
                                                  self._filling.get((i,j), 0))

        self._out_labels = labels

    def _shrink(self):
        r"""
        Compute the labels on the boundary near the origin, and the filling.

        TESTS::

            sage: RuleRSK = GrowthDiagram.rules.RSK()
            sage: filling = [[0,0,1,0,0,0,0], [0,1,0,0,0,0,0], [1,0,0,0,0,0,0],
            ....:            [0,0,0,1,0,0,0], [0,0,0,0,0,0,1],
            ....:            [0,0,0,0,0,1,0], [0,0,0,0,1,0,0]]
            sage: G = GrowthDiagram(RuleRSK, filling)
            sage: list(GrowthDiagram(RuleRSK, labels=G._out_labels)) == filling
            True

            sage: labels = [[], [1], []]
            sage: G = GrowthDiagram(RuleRSK, labels=labels)  # indirect doctest
            sage: G._filling
            {(0, 0): 1}
            sage: G._in_labels
            [[], [], []]

            sage: labels = [[], [1], [2], [2,1], [1,1], [1], []]
            sage: G = GrowthDiagram(RuleRSK, labels=labels)  # indirect doctest
            sage: G._filling
            {(0, 1): 1, (1, 2): 1, (2, 0): 1}
            sage: G._in_labels
            [[], [], [], [], [], [], []]

            sage: labels = [[], [1], [2], [3], [3, 1], [3, 2], [4, 2], [4, 1], [3, 1], [2, 1], [1, 1], [1], []]
            sage: G = GrowthDiagram(RuleRSK, labels=labels)  # indirect doctest
            sage: G._filling
            {(0, 1): 1, (1, 2): 1, (2, 5): 1, (3, 0): 1, (4, 3): 1, (5, 4): 1}

            sage: labels = [[],[1],[1],[2],[2],[2,1],[2]]
            sage: G = GrowthDiagram(RuleRSK, labels=labels)
            Traceback (most recent call last):
            ...
            ValueError: can only determine the shape of the growth diagram
             if ranks of successive labels differ
            sage: G = GrowthDiagram(RuleRSK, shape=[3,2,1], labels=labels)  # indirect doctest
            sage: G._filling
            {(1, 0): 1}
            sage: G._in_labels
            [[], [], [], [], [1], [1], [2]]

            sage: labels = [[], [1],[1],[2],[2],[2,1],[2],[2,1],[1,1],[2,1],[1,1]]
            sage: G = GrowthDiagram(RuleRSK, shape=[5,4,3,2,1], labels=labels)  # indirect doctest
            sage: G._filling
            {(1, 2): 1, (2, 1): 1, (4, 0): 1}
            sage: G._in_labels
            [[], [], [], [], [], [], [1], [1], [1], [1, 1], [1, 1]]

            sage: labels = [[], [1],[1],[2],[2],[2,1],[2],[2,1],[1,1],[2,1],[1,1]]
            sage: G = GrowthDiagram(RuleRSK, shape=SkewPartition([[5,4,3,2,1],[3,2,1]]), labels=labels)  # indirect doctest
            sage: G._filling
            {(1, 2): 1, (2, 1): 1, (4, 0): 1}
            sage: G._in_labels
            [[], [], [], [1], [1], [1], [1], [1], [1], [1, 1], [1, 1]]
        """
        F = dict()
        labels = list(self._out_labels)
        l = len(self._lambda)
        rule = self.rule
        if rule.has_multiple_edges:
            for r in range(l):
                for c in range(self._lambda[l-r-1]+r, self._mu[l-r-1]+r, -1):
                    j = l-r-1
                    i = c-r-1
                    (labels[2*c-1],
                     labels[2*c],
                     labels[2*c+1], v) = rule.backward_rule(labels[2*c-2],
                                                            labels[2*c-1],
                                                            labels[2*c],
                                                            labels[2*c+1],
                                                            labels[2*c+2])
                    if v != 0:
                        F[(i,j)] = v

        else:
            for r in range(l):
                for c in range(self._lambda[l-r-1]+r, self._mu[l-r-1]+r, -1):
                    j = l-r-1
                    i = c-r-1
                    labels[c], v = rule.backward_rule(labels[c-1],
                                                      labels[c],
                                                      labels[c+1])
                    if v != 0:
                        F[(i,j)] = v

        self._in_labels = labels
        self._filling = F

######################################################################
# ABC for rules of growth diagrams
######################################################################

class Rule(UniqueRepresentation):
    r"""
    Generic base class for a rule for a growth diagram.

    Subclasses may provide the following attributes:

    - ``zero`` -- the zero element of the vertices of the graphs

    - ``r`` -- (default: 1) the parameter in the equation `DU - UD = rI`

    - ``has_multiple_edges`` -- (default: ``False``) if the dual
      graded graph has multiple edges and therefore edges are
      triples consisting of two vertices and a label.

    - ``zero_edge`` -- (default: 0) the zero label of the
      edges of the graphs used for degenerate edges.  It is
      allowed to use this label also for other edges.

    Subclasses may provide the following methods:

    - ``normalize_vertex`` -- a function that converts its input to a
      vertex.

    - ``vertices`` -- a function that takes a non-negative integer
      as input and returns the list of vertices on this rank.

    - ``rank`` -- the rank function of the dual graded graphs.

    - ``forward_rule`` -- a function with input ``(y, t, x,
      content)`` or ``(y, e, t, f, x, content)`` if
      ``has_multiple_edges`` is ``True``.  ``(y, e, t)`` is an
      edge in the graph `P`, ``(t, f, x)`` an edge in the graph
      ``Q``.  It should return the fourth vertex ``z``, or, if
      ``has_multiple_edges`` is ``True``, the path ``(g, z, h)``
      from ``y`` to ``x``.

    - ``backward_rule`` -- a function with input ``(y, z, x)`` or
      ``(y, g, z, h, x)`` if ``has_multiple_edges`` is ``True``.
      ``(y, g, z)`` is an edge in the graph `Q`, ``(z, h, x)`` an
      edge in the graph ``P``.  It should return the fourth
      vertex and the content ``(t, content)``, or, if
      ``has_multiple_edges`` is ``True``, the path from ``y`` to
      ``x`` and the content as ``(e, t, f, content)``.

    - ``is_P_edge``, ``is_Q_edge`` -- functions that take two
      vertices as arguments and return ``True`` or ``False``, or,
      if multiple edges are allowed, the list of edge labels of
      the edges from the first vertex to the second in the
      respective graded graph.  These are only used for checking
      user input and providing the dual graded graph, and are
      therefore not mandatory.

    Note that the class :class:`GrowthDiagram` is able to use
    partially implemented subclasses just fine.  Suppose that
    ``MyRule`` is such a subclass.  Then:

    - ``GrowthDiagram(MyRule, my_filling)`` requires only an
      implementation of ``forward_rule``, ``zero`` and possibly
      ``has_multiple_edges``.

    - ``GrowthDiagram(MyRule, labels=my_labels, shape=my_shape)``
      requires only an implementation of ``backward_rule`` and
      possibly ``has_multiple_edges``, provided that the labels
      ``my_labels`` are given as needed by ``backward_rule``.

    - ``GrowthDiagram(MyRule, labels=my_labels)`` additionally needs
      an implementation of ``rank`` to deduce the shape.

    In particular, this allows to implement rules which do not quite
    fit Fomin's notion of dual graded graphs.  An example would be
    Bloom and Saracino's variant of the RSK correspondence [BS2012]_,
    where a backward rule is not available.

    Similarly:

    - ``MyRule.P_graph`` only requires an implementation of
      ``vertices``, ``is_P_edge`` and possibly ``has_multiple_edges``
      is required, mutatis mutandis for ``MyRule.Q_graph``.

    - ``MyRule._check_duality`` requires ``P_graph`` and ``Q_graph``.

    In particular, this allows to work with dual graded graphs
    without local rules.
    """
    has_multiple_edges = False          # override when necessary
    zero_edge = 0                       # override when necessary
    r = 1                               # override when necessary

    def normalize_vertex(self, v):      # override when necessary
        r"""
        Return ``v`` as a vertex of the dual graded graph.

        This is a default implementation, returning its argument.

        EXAMPLES::

            sage: from sage.combinat.growth import Rule
            sage: Rule().normalize_vertex("hello") == "hello"
            True
        """
        return v

    def __call__(self, *args, **kwds):
        r"""
        Return the growth diagram corresponding to the parameters.

        This provides a shorthand for calling :class:`GrowthDiagram`
        directly.

        EXAMPLES::

            sage: RuleRSK = GrowthDiagram.rules.RSK()
            sage: RuleRSK([2,3,1], shape=[3,2,2])
            0  0  1
            1  0
            0  1

            sage: RuleRSK(labels=[[], [1], [2], [1], [], [1], []])
            0  0  1
            1  0
            0  1
        """
        return GrowthDiagram(self, *args, **kwds)

    def _check_duality(self, n):
        r"""
        Raise an error if the graphs are not `r`-dual at level ``n``.

        `P` and `Q` are `r`-dual if `DU = UD + rI` on the free
        `\ZZ`-module `\ZZ[V]`, where `D` is the down operator of `Q`,
        assigning to each vertex the formal sum of its predecessors,
        `U` is the up operator of `P`, assigning to each vertex the
        formal sum of its successors, and `I` is the identity
        operator.

        INPUT:

        - ``n`` -- a positive integer specifying which rank of
          the graph to test

        EXAMPLES:

        For binary words, we have indeed provided dual graded graphs::

            sage: BinaryWord = GrowthDiagram.rules.BinaryWord()
            sage: BinaryWord._check_duality(3)

        The following two graphs are not `1`-dual::

            sage: from sage.combinat.growth import Rule
            sage: class RuleWrong(Rule):
            ....:     def vertices(self, n): return Partitions(n)
            ....:     def is_Q_edge(self, v, w):
            ....:         return (v, w) in [([1],[2]), ([2],[3])]
            ....:     def is_P_edge(self, v, w):
            ....:         return (v, w) in [([1],[2]), ([1],[1,1]), ([2],[3])]

            sage: RuleWrong()._check_duality(2)
            Traceback (most recent call last):
            ...
            ValueError: D U - U D differs from 1 I for vertex [2]:
            D U = [[2]]
            U D + 1 I = [[1, 1], [2], [2]]
        """
        if self.has_multiple_edges:
            def check_vertex(w, P, Q):
                DUw = [v[0] for uw in P.outgoing_edges(w) for v in Q.incoming_edges(uw[1])]
                UDw = [v[1] for lw in Q.incoming_edges(w) for v in P.outgoing_edges(lw[0])]
                UDw.extend([w]*self.r)
                if sorted(DUw) != sorted(UDw):
                    raise ValueError("D U - U D differs from %s I for vertex %s:\n"
                                     "D U = %s\n"
                                     "U D + %s I = %s"
                                     % (self.r, w, DUw, self.r, UDw))
        else:
            def check_vertex(w, P, Q):
                DUw = [v for uw in P.upper_covers(w) for v in Q.lower_covers(uw)]
                UDw = [v for lw in Q.lower_covers(w) for v in P.upper_covers(lw)]
                UDw.extend([w]*self.r)
                if sorted(DUw) != sorted(UDw):
                    raise ValueError("D U - U D differs from %s I for vertex %s:\n"
                                     "D U = %s\n"
                                     "U D + %s I = %s"
                                     % (self.r, w, DUw, self.r, UDw))

        P = self.P_graph(n + 2)
        Q = self.Q_graph(n + 2)
        for w in self.vertices(n):
            check_vertex(w, P, Q)

    def P_graph(self, n):
        r"""
        Return the first ``n`` levels of the first dual graded graph.

        The non-degenerate edges in the vertical direction come from
        this graph.

        EXAMPLES::

            sage: Domino = GrowthDiagram.rules.Domino()
            sage: Domino.P_graph(3)
            Finite poset containing 8 elements
        """
        if self.has_multiple_edges:
            D = DiGraph([(x,y,e) for k in range(n-1)
                         for x in self.vertices(k)
                         for y in self.vertices(k+1)
                         for e in self.is_P_edge(x, y)], multiedges=True)
            # unfortunately, layout_acyclic will not show multiple edges
            # D.layout_default = D.layout_acyclic
            return D
        else:
            return Poset(([w for k in range(n) for w in self.vertices(k)],
                          self.is_P_edge), cover_relations=True)

    def Q_graph(self, n):
        r"""
        Return the first ``n`` levels of the second dual graded graph.

        The non-degenerate edges in the horizontal direction come
        from this graph.

        EXAMPLES::

            sage: Domino = GrowthDiagram.rules.Domino()
            sage: Q = Domino.Q_graph(3); Q
            Finite poset containing 8 elements

            sage: Q.upper_covers(Partition([1,1]))
            [[1, 1, 1, 1], [3, 1], [2, 2]]
        """
        if self.has_multiple_edges:
            D = DiGraph([(x, y, e) for k in range(n - 1)
                         for x in self.vertices(k)
                         for y in self.vertices(k + 1)
                         for e in self.is_Q_edge(x, y)], multiedges=True)
            # unfortunately, layout_acyclic will not show multiple edges
            # D.layout_default = D.layout_acyclic
            return D
        else:
            return Poset(([w for k in range(n) for w in self.vertices(k)],
                          self.is_Q_edge), cover_relations=True)

######################################################################
# Specific rules of growth diagrams
######################################################################


class RuleShiftedShapes(Rule):
    r"""
    A class modelling the Schensted correspondence for shifted
    shapes.

    This agrees with Sagan [Sag1987]_ and Worley's [Wor1984]_, and
    Haiman's [Hai1989]_ insertion algorithms, see Proposition 4.5.2
    of [Fom1995]_.

    EXAMPLES::

        sage: Shifted = GrowthDiagram.rules.ShiftedShapes()
        sage: GrowthDiagram(Shifted, [3,1,2])
        0  1  0
        0  0  1
        1  0  0

    The vertices of the dual graded graph are shifted shapes::

        sage: Shifted.vertices(3)
        Partitions of the integer 3 satisfying constraints max_slope=-1

    Let us check the example just before Corollary 3.2 in [Sag1987]_.
    Note that, instead of passing the rule to :class:`GrowthDiagram`,
    we can also call the rule to create growth diagrams::

        sage: G = Shifted([2,6,5,1,7,4,3])
        sage: G.P_chain()
        [[], 0, [1], 0, [2], 0, [3], 0, [3, 1], 0, [3, 2], 0, [4, 2], 0, [5, 2]]
        sage: G.Q_chain()
        [[], 1, [1], 2, [2], 1, [2, 1], 3, [3, 1], 2, [4, 1], 3, [4, 2], 3, [5, 2]]

    TESTS::

        sage: Shifted = GrowthDiagram.rules.ShiftedShapes()
        sage: Shifted.zero
        []

        sage: Shifted._check_duality(4)

    Check that the rules are bijective::

        sage: all(Shifted(labels=Shifted(pi).out_labels()).to_word() == pi
        ....:     for pi in Permutations(5))
        True
        sage: pi = Permutations(10).random_element()
        sage: G = Shifted(pi)
        sage: list(Shifted(labels=G.out_labels())) == list(G)
        True
    """
    zero = _make_partition([])
    has_multiple_edges = True

    def normalize_vertex(self, v):
        r"""
        Return ``v`` as a partition.

        EXAMPLES::

            sage: Shifted = GrowthDiagram.rules.ShiftedShapes()
            sage: Shifted.normalize_vertex([3,1]).parent()
            Partitions
        """
        return _make_partition(v)

    def vertices(self, n):
        r"""
        Return the vertices of the dual graded graph on level ``n``.

        EXAMPLES::

            sage: Shifted = GrowthDiagram.rules.ShiftedShapes()
            sage: Shifted.vertices(3)
            Partitions of the integer 3 satisfying constraints max_slope=-1
        """
        if n == 0:
            return [self.zero]
        else:
            return Partitions(n, max_slope=-1)

    def rank(self, v):
        r"""
        Return the rank of ``v``: the size of the shifted partition.

        EXAMPLES::

            sage: Shifted = GrowthDiagram.rules.ShiftedShapes()
            sage: Shifted.rank(Shifted.vertices(3)[0])
            3
        """
        return v.size()

    def is_Q_edge(self, v, w):
        r"""
        Return whether ``(v, w)`` is a `Q`-edge of ``self``.

        ``(v, w)`` is an edge if ``w`` is obtained from ``v`` by adding a
        cell.  It is a black (color 1) edge, if the cell is on the
        diagonal, otherwise it can be blue or red (color 2 or 3).

        EXAMPLES::

            sage: Shifted = GrowthDiagram.rules.ShiftedShapes()
            sage: v = Shifted.vertices(2)[0]; v
            [2]
            sage: [(w, Shifted.is_Q_edge(v, w)) for w in Shifted.vertices(3)]
            [([3], [2, 3]), ([2, 1], [1])]
            sage: all(Shifted.is_Q_edge(v, w) == [] for w in Shifted.vertices(4))
            True
        """
        if self.rank(v) + 1 != self.rank(w):
            return []
        try:
            l = SkewPartition([w, v]).cells()
        except ValueError:
            return []
        else:
            if l[0][1] == 0:
                return [1]   # black
            else:
                return [2,3] # blue, red

    def is_P_edge(self, v, w):
        r"""
        Return whether ``(v, w)`` is a `P`-edge of ``self``.

        ``(v, w)`` is an edge if ``w`` contains ``v``.

        EXAMPLES::

            sage: Shifted = GrowthDiagram.rules.ShiftedShapes()
            sage: v = Shifted.vertices(2)[0]; v
            [2]
            sage: [w for w in Shifted.vertices(3) if Shifted.is_P_edge(v, w)]
            [[3], [2, 1]]
        """
        if self.rank(v) + 1 != self.rank(w):
            return []
        return [0] if w.contains(v) else []

    def P_symbol(self, P_chain):
        r"""
        Return the labels along the vertical boundary of a rectangular
        growth diagram as a shifted tableau.

        EXAMPLES:

        Check the example just before Corollary 3.2 in [Sag1987]_::

            sage: Shifted = GrowthDiagram.rules.ShiftedShapes()
            sage: G = Shifted([2,6,5,1,7,4,3])
            sage: G.P_symbol().pp()
            1  2  3  6  7
               4  5

        Check the example just before Corollary 8.2 in [SS1990]_::

            sage: T = ShiftedPrimedTableau([[4],[1],[5]], skew=[3,1])
            sage: T.pp()
             .  .  .  4
                .  1
                   5
            sage: U = ShiftedPrimedTableau([[1],[3.5],[5]], skew=[3,1])
            sage: U.pp()
             .  .  .  1
                .  4'
                   5
            sage: Shifted = GrowthDiagram.rules.ShiftedShapes()
            sage: labels = [mu if is_even(i) else 0 for i, mu in enumerate(T.to_chain()[::-1])] + U.to_chain()[1:]
            sage: G = Shifted({(1,2):1, (2,1):1}, shape=[5,5,5,5,5], labels=labels)
            sage: G.P_symbol().pp()
             .  .  .  .  2
                .  .  1  3
                   .  4  5

        """
        chain = P_chain[::2]
        shape = chain[-1]
        T = [[None for _ in range(r)] for r in shape]
        for i in range(1,len(chain)):
            la = chain[i]
            mu = chain[i-1]
            mu += [0]*(len(la) - len(mu))

            for r in range(len(la)):
                for c in range(mu[r], la[r]):
                    T[r][c] = i

        skew = _make_partition([row.count(None) for row in T])
        T = [[e for e in row if e is not None] for row in T]
        return ShiftedPrimedTableau(T, skew=skew)

    def Q_symbol(self, Q_chain):
        r"""
        Return the labels along the horizontal boundary of a rectangular
        growth diagram as a skew tableau.

        EXAMPLES:

        Check the example just before Corollary 3.2 in [Sag1987]_::

            sage: Shifted = GrowthDiagram.rules.ShiftedShapes()
            sage: G = Shifted([2,6,5,1,7,4,3])
            sage: G.Q_symbol().pp()
            1  2  4' 5  7'
               3  6'

        Check the example just before Corollary 8.2 in [SS1990]_::

            sage: T = ShiftedPrimedTableau([[4],[1],[5]], skew=[3,1])
            sage: T.pp()
             .  .  .  4
                .  1
                   5
            sage: U = ShiftedPrimedTableau([[1],[3.5],[5]], skew=[3,1])
            sage: U.pp()
             .  .  .  1
                .  4'
                   5
            sage: Shifted = GrowthDiagram.rules.ShiftedShapes()
            sage: labels = [mu if is_even(i) else 0 for i, mu in enumerate(T.to_chain()[::-1])] + U.to_chain()[1:]
            sage: G = Shifted({(1,2):1, (2,1):1}, shape=[5,5,5,5,5], labels=labels)
            sage: G.Q_symbol().pp()
             .  .  .  .  2
                .  .  1  4'
                   .  3' 5'

        """
        chain = Q_chain
        shape = chain[-1]
        T = [[None for _ in range(r)] for r in shape]
        for i in range(1,(len(chain)+1)//2):
            la = chain[2*i]
            if chain[2*i-1] == 3:
                prime = 0.5
            else:
                prime = 0
            mu = chain[2*(i-1)]
            mu += [0]*(len(la) - len(mu))

            for r in range(len(la)):
                for c in range(mu[r], la[r]):
                    T[r][c] = i - prime

        skew = _make_partition([row.count(None) for row in T])
        T = [[e for e in row if e is not None] for row in T]
        return ShiftedPrimedTableau(T, skew=skew)

    def forward_rule(self, y, e, t, f, x, content):
        r"""
        Return the output path given two incident edges and the content.

        See [Fom1995]_ Lemma 4.5.1, page 38.

        INPUT:

        - ``y, e, t, f, x`` -- a path of three partitions and two
          colors from a cell in a growth diagram, labelled as::

              t f x
              e
              y

        - ``content`` -- `0` or `1`; the content of the cell

        OUTPUT:

        The two colors and the fourth partition ``g``, ``z``, ``h``
        according to Sagan-Worley insertion.

        EXAMPLES::

            sage: Shifted = GrowthDiagram.rules.ShiftedShapes()
            sage: Shifted.forward_rule([], 0, [], 0, [], 1)
            (1, [1], 0)

            sage: Shifted.forward_rule([1], 0, [1], 0, [1], 1)
            (2, [2], 0)

        if ``x != y``::

            sage: Shifted.forward_rule([3], 0, [2], 1, [2,1], 0)
            (1, [3, 1], 0)

            sage: Shifted.forward_rule([2,1], 0, [2], 2, [3], 0)
            (2, [3, 1], 0)

        if ``x == y != t``::

            sage: Shifted.forward_rule([3], 0, [2], 2, [3], 0)
            (1, [3, 1], 0)

            sage: Shifted.forward_rule([3,1], 0, [2,1], 2, [3,1], 0)
            (2, [3, 2], 0)

            sage: Shifted.forward_rule([2,1], 0, [2], 1, [2,1], 0)
            (3, [3, 1], 0)

            sage: Shifted.forward_rule([3], 0, [2], 3, [3], 0)
            (3, [4], 0)

        """
        if e != 0:
            raise ValueError("the P-graph should not be colored")
        h = 0
        if x == t == y:
            if f != 0:
                raise ValueError("degenerate edge f should have color 0")
            if content == 0:
                g, z = 0, x
            elif content == 1:
                if not x:
                    g, z = 1, _Partitions(x).add_cell(0) # black
                else:
                    g, z = 2, _make_partition(x).add_cell(0) # blue
            else:
                raise NotImplementedError
        elif content != 0:
            raise ValueError("for y=%s, t=%s, x=%s, the content should be 0 but is %s"
                             % (y, t, x, content))
        elif x != t == y:
            g, z = f, x
        elif x == t != y:
            if f != 0:
                raise ValueError("degenerate edge f should have color 0")
            g, z = f, y
        else:
            if x != y:
                row = SkewPartition([x, t]).cells()[0][0]
                g, z = f, _make_partition(y).add_cell(row)
            elif x == y != t and f == 2: # blue
                row = 1+SkewPartition([x, t]).cells()[0][0]
                if row == len(y):
                    g, z = 1, _make_partition(y).add_cell(row) # black
                else:
                    g, z = 2, _make_partition(y).add_cell(row) # blue
            elif x == y != t and f in [1, 3]: # black or red
                c = SkewPartition([x, t]).cells()[0]
                col = c[0] + c[1] + 1
                # print y, t, x, c, col
                for i in range(len(y)):
                    if i + y[i] == col:
                        z = y[:i] + [y[i]+1] + y[i+1:]
                        break
                g = 3
            else:
                raise NotImplementedError
        return g, _make_partition(z), h

    def backward_rule(self, y, g, z, h, x):
        r"""
        Return the input path and the content given two incident edges.

        See [Fom1995]_ Lemma 4.5.1, page 38.

        INPUT:

        - ``y, g, z, h, x`` -- a path of three partitions and two
          colors from a cell in a growth diagram, labelled as::

                  x
                  h
              y g z

        OUTPUT:

        A tuple ``(e, t, f, content)`` consisting of the shape ``t``
        of the fourth word, the colours of the incident edges and the
        content of the cell according to Sagan - Worley insertion.

        EXAMPLES::

            sage: Shifted = GrowthDiagram.rules.ShiftedShapes()
            sage: Shifted.backward_rule([], 1, [1], 0, [])
            (0, [], 0, 1)

            sage: Shifted.backward_rule([1], 2, [2], 0, [1])
            (0, [1], 0, 1)

        if ``x != y``::

            sage: Shifted.backward_rule([3], 1, [3, 1], 0, [2,1])
            (0, [2], 1, 0)

            sage: Shifted.backward_rule([2,1], 2, [3, 1], 0, [3])
            (0, [2], 2, 0)

        if ``x == y != t``::

            sage: Shifted.backward_rule([3], 1, [3, 1], 0, [3])
            (0, [2], 2, 0)

            sage: Shifted.backward_rule([3,1], 2, [3, 2], 0, [3,1])
            (0, [2, 1], 2, 0)

            sage: Shifted.backward_rule([2,1], 3, [3, 1], 0, [2,1])
            (0, [2], 1, 0)

            sage: Shifted.backward_rule([3], 3, [4], 0, [3])
            (0, [2], 3, 0)
        """
        if h != 0:
            raise ValueError("the P-graph should not be colored")

        if x == y == z:
            if g != 0:
                raise ValueError("degenerate edge g should have color 0")
            return (0, x, 0, 0)
        elif x == z != y:
            return (0, y, g, 0)
        elif x != z == y:
            if g != 0:
                raise ValueError("degenerate edge g should have color 0")
            return (0, x, 0, 0)
        else:
            if x != y:
                row = SkewPartition([z, x]).cells()[0][0]
                return (0, _make_partition(y).remove_cell(row), g, 0)
            else:
                row, col = SkewPartition([z, x]).cells()[0]
                if row > 0 and g in [1, 2]: # black or blue
                    return (0, _make_partition(y).remove_cell(row-1), 2, 0)
                elif row == 0 and g in [1, 2]: # black or blue
                    return (0, y, 0, 1)
                else:
                    # find last cell in column col-1
                    for i in range(len(y)-1,-1,-1):
                        if i + y[i] == col + row:
                            if y[i] == 1:
                                t = y[:i]
                                return (0, t, 1, 0)
                            else:
                                t = y[:i] + [y[i]-1] + y[i+1:]
                                return (0, t, 3, 0)
                    raise ValueError("this should not happen")

class RuleLLMS(Rule):
    r"""
    A rule modelling the Schensted correspondence for affine
    permutations.

    EXAMPLES::

        sage: LLMS3 = GrowthDiagram.rules.LLMS(3)
        sage: GrowthDiagram(LLMS3, [3,1,2])
        0  1  0
        0  0  1
        1  0  0

    The vertices of the dual graded graph are
    :class:`~sage.combinat.core.Cores`::

        sage: LLMS3.vertices(4)
        3-Cores of length 4

    Let us check example of Figure 1 in [LS2007]_.  Note that,
    instead of passing the rule to :class:`GrowthDiagram`, we can
    also call the rule to create growth diagrams::

        sage: G = LLMS3([4,1,2,6,3,5]); G
        0  1  0  0  0  0
        0  0  1  0  0  0
        0  0  0  0  1  0
        1  0  0  0  0  0
        0  0  0  0  0  1
        0  0  0  1  0  0

    The :meth:`P_symbol` is a
    :class:`~sage.combinat.k_tableau.StrongTableau`::

        sage: G.P_symbol().pp()
        -1 -2 -3 -5
         3  5
        -4 -6
         5
         6

    The :meth:`Q_symbol` is a
    :class:`~sage.combinat.k_tableau.WeakTableau`::

        sage: G.Q_symbol().pp()
        1  3  4  5
        2  5
        3  6
        5
        6

    Let us also check Example 6.2 in [LLMSSZ2013]_::

        sage: G = LLMS3([4,1,3,2])
        sage: G.P_symbol().pp()
        -1 -2  3
        -3
        -4

        sage: G.Q_symbol().pp()
        1  3  4
        2
        3

    TESTS::

        sage: LLMS3 = GrowthDiagram.rules.LLMS(3)
        sage: LLMS3.zero
        []
    """
    zero_edge = None  # to prevent confusion with the edge labelled with content 0
    has_multiple_edges = True

    def __init__(self, k):
        r"""
        Initialize ``self``.

        TESTS::

            sage: LLMS3 = GrowthDiagram.rules.LLMS(3)
            sage: TestSuite(LLMS3).run()
        """
        self.k = k
        self.zero = Core([], k)

    def normalize_vertex(self, v):
        r"""
        Convert ``v`` to a `k`-core.

        EXAMPLES::

            sage: LLMS3 = GrowthDiagram.rules.LLMS(3)
            sage: LLMS3.normalize_vertex([3,1]).parent()
            3-Cores of length 3
        """
        return Core(v, self.k)

    def rank(self, v):
        r"""
        Return the rank of ``v``: the length of the core.

        EXAMPLES::

            sage: LLMS3 = GrowthDiagram.rules.LLMS(3)
            sage: LLMS3.rank(LLMS3.vertices(3)[0])
            3
        """
        return v.length()

    def vertices(self, n):
        r"""
        Return the vertices of the dual graded graph on level ``n``.

        EXAMPLES::

            sage: LLMS3 = GrowthDiagram.rules.LLMS(3)
            sage: LLMS3.vertices(2)
            3-Cores of length 2
        """
        return Cores(self.k, length=n)

    def is_Q_edge(self, v, w):
        r"""
        Return whether ``(v, w)`` is a `Q`-edge of ``self``.

        ``(v, w)`` is an edge if ``w`` is a weak cover of ``v``, see
        :meth:`~sage.combinat.core.Core.weak_covers()`.

        EXAMPLES::

            sage: LLMS4 = GrowthDiagram.rules.LLMS(4)
            sage: v = LLMS4.vertices(3)[1]; v
            [2, 1]
            sage: [w for w in LLMS4.vertices(4) if len(LLMS4.is_Q_edge(v, w)) > 0]
            [[2, 2], [3, 1, 1]]
            sage: all(LLMS4.is_Q_edge(v, w) == [] for w in LLMS4.vertices(5))
            True
        """
        return [None] if w in v.weak_covers() else []

    def is_P_edge(self, v, w):
        r"""
        Return whether ``(v, w)`` is a `P`-edge of ``self``.

        For two k-cores v and w containing v, there are as many edges as
        there are components in the skew partition w/v.  These
        components are ribbons, and therefore contain a unique cell
        with maximal content.  We index the edge with this content.

        EXAMPLES::

            sage: LLMS4 = GrowthDiagram.rules.LLMS(4)
            sage: v = LLMS4.vertices(2)[0]; v
            [2]
            sage: [(w, LLMS4.is_P_edge(v, w)) for w in LLMS4.vertices(3)]
            [([3], [2]), ([2, 1], [-1]), ([1, 1, 1], [])]
            sage: all(LLMS4.is_P_edge(v, w) == [] for w in LLMS4.vertices(4))
            True
        """
        if w in v.strong_covers():
            T = SkewPartition([w.to_partition(), v.to_partition()])
            return [max([j-i for i,j in c]) for c in T.cell_poset().connected_components()]
        else:
            return []

    def P_symbol(self, P_chain):
        r"""
        Return the labels along the vertical boundary of a
        rectangular growth diagram as a skew
        :class:`~sage.combinat.k_tableau.StrongTableau`.

        EXAMPLES::

            sage: LLMS4 = GrowthDiagram.rules.LLMS(4)
            sage: G = LLMS4([3,4,1,2])
            sage: G.P_symbol().pp()
            -1 -2
            -3 -4
        """
        C = P_chain
        T = SkewTableau(chain=C[::2])
        S = T.to_list()
        for entry, content in enumerate(C[1::2], 1):
            for i,j in T.cells_containing(entry):
                if j-i == content:
                    S[i][j] = -S[i][j]
                    break
        return StrongTableau(S, self.k-1)

    def Q_symbol(self, Q_chain):
        r"""
        Return the labels along the horizontal boundary of a
        rectangular growth diagram as a skew
        :class:`~sage.combinat.k_tableau.WeakTableau`.

        EXAMPLES::

            sage: LLMS4 = GrowthDiagram.rules.LLMS(4)
            sage: G = LLMS4([3,4,1,2])
            sage: G.Q_symbol().pp()
            1 2
            3 4
        """
        return WeakTableau(SkewTableau(chain=Q_chain[::2]), self.k-1)

    def forward_rule(self, y, e, t, f, x, content):
        r"""
        Return the output path given two incident edges and the content.

        See [LS2007]_ Section 3.4 and [LLMSSZ2013]_ Section 6.3.

        INPUT:

        - ``y, e, t, f, x`` -- a path of three partitions and two
          colors from a cell in a growth diagram, labelled as::

              t f x
              e
              y

        - ``content`` -- `0` or `1`; the content of the cell

        OUTPUT:

        The two colors and the fourth partition g, z, h according to
        LLMS insertion.

        EXAMPLES::

            sage: LLMS3 = GrowthDiagram.rules.LLMS(3)
            sage: LLMS4 = GrowthDiagram.rules.LLMS(4)

            sage: Z = LLMS3.zero
            sage: LLMS3.forward_rule(Z, None, Z, None, Z, 0)
            (None, [], None)

            sage: LLMS3.forward_rule(Z, None, Z, None, Z, 1)
            (None, [1], 0)

            sage: Y = Core([3,1,1], 3)
            sage: LLMS3.forward_rule(Y, None, Y, None, Y, 1)
            (None, [4, 2, 1, 1], 3)

        if ``x != y``::

            sage: Y = Core([1,1], 3); T = Core([1], 3); X = Core([2], 3)
            sage: LLMS3.forward_rule(Y, -1, T, None, X, 0)
            (None, [2, 1, 1], -1)

            sage: Y = Core([2], 4); T = Core([1], 4); X = Core([1,1], 4)
            sage: LLMS4.forward_rule(Y, 1, T, None, X, 0)
            (None, [2, 1], 1)

            sage: Y = Core([2,1,1], 3); T = Core([2], 3); X = Core([3,1], 3)
            sage: LLMS3.forward_rule(Y, -1, T, None, X, 0)
            (None, [3, 1, 1], -2)

        if ``x == y != t``::

            sage: Y = Core([1], 3); T = Core([], 3); X = Core([1], 3)
            sage: LLMS3.forward_rule(Y, 0, T, None, X, 0)
            (None, [1, 1], -1)

            sage: Y = Core([1], 4); T = Core([], 4); X = Core([1], 4)
            sage: LLMS4.forward_rule(Y, 0, T, None, X, 0)
            (None, [1, 1], -1)

            sage: Y = Core([2,1], 4); T = Core([1,1], 4); X = Core([2,1], 4)
            sage: LLMS4.forward_rule(Y, 1, T, None, X, 0)
            (None, [2, 2], 0)
        """
        if f is not None:
            raise ValueError("the Q-graph should not be colored")
        g = None
        if x == t == y:
            if e is not None:
                raise ValueError("degenerate edge e should have color None")
            if content == 0:
                z, h = x, None
            elif content == 1:
                if t.size() == 0:
                    z = t.affine_symmetric_group_simple_action(0)
                else:
                    z = t.affine_symmetric_group_simple_action(t[0] % self.k)
                h = z[0] - 1
            else:
                assert False, "BUG in RuleLLMS"
        elif content != 0:
            raise ValueError("for y=%s, t=%s, x=%s, the content should be 0 but is %s"
                             % (y, t, x, content))
        elif x != t == y:
            if e is not None:
                raise ValueError("degenerate edge e should have color None")
            z, h = x, e
        elif x == t != y:
            z, h = y, e
        else: #  x != t and y != t
            qx = SkewPartition([x.to_partition(), t.to_partition()])
            qy = SkewPartition([y.to_partition(), t.to_partition()])
            if not all(c in qx.cells() for c in qy.cells()):
                res = [(j-i) % self.k for i,j in qx.cells()]
                assert len(set(res)) == 1
                r = res[0]
                z = y.affine_symmetric_group_simple_action(r)
                if e % self.k == r:
                    h = e-1
                else:
                    h = e
            elif x == y != t:
                # the addable cell with largest content at most e
                cprime = sorted([c for c in y.to_partition().addable_cells()
                                 if c[1]-c[0] <= e],
                                key = lambda c: -(c[1]-c[0]))[0]
                h = cprime[1] - cprime[0]
                z = y.affine_symmetric_group_simple_action(h % self.k)

        return g, z, h

class RuleBinaryWord(Rule):
    r"""
    A rule modelling a Schensted-like correspondence for binary words.

    EXAMPLES::

        sage: BinaryWord = GrowthDiagram.rules.BinaryWord()
        sage: GrowthDiagram(BinaryWord, [3,1,2])
        0  1  0
        0  0  1
        1  0  0

    The vertices of the dual graded graph are binary words::

        sage: BinaryWord.vertices(3)
        [word: 100, word: 101, word: 110, word: 111]

    Note that, instead of passing the rule to :class:`GrowthDiagram`,
    we can also use call the rule to create growth diagrams.  For
    example::

        sage: BinaryWord([2,4,1,3]).P_chain()
        [word: , word: 1, word: 10, word: 101, word: 1101]
        sage: BinaryWord([2,4,1,3]).Q_chain()
        [word: , word: 1, word: 11, word: 110, word: 1101]

    The Kleitman Greene invariant is the descent word, encoded by the
    positions of the zeros::

        sage: pi = Permutation([4,1,8,3,6,5,2,7,9])
        sage: G = BinaryWord(pi); G
        0  1  0  0  0  0  0  0  0
        0  0  0  0  0  0  1  0  0
        0  0  0  1  0  0  0  0  0
        1  0  0  0  0  0  0  0  0
        0  0  0  0  0  1  0  0  0
        0  0  0  0  1  0  0  0  0
        0  0  0  0  0  0  0  1  0
        0  0  1  0  0  0  0  0  0
        0  0  0  0  0  0  0  0  1
        sage: pi.descents()
        [1, 3, 5, 6]

    TESTS::

        sage: BinaryWord = GrowthDiagram.rules.BinaryWord()
        sage: BinaryWord.zero
        word:
        sage: G = BinaryWord(labels=[[1,1],[1,1,0],[0,1]])
        Traceback (most recent call last):
        ...
        ValueError: 01 has smaller rank than 110 but is not covered by it in P

        sage: G = BinaryWord(labels=[[1,1],[1,0,1],[0,1]])
        Traceback (most recent call last):
        ...
        ValueError: 11 has smaller rank than 101 but is not covered by it in Q

    Check duality::

        sage: BinaryWord._check_duality(4)

    Check that the rules are bijective::

        sage: all(BinaryWord(labels=BinaryWord(pi).out_labels()).to_word()
        ....:      == pi for pi in Permutations(4))
        True
        sage: pi = Permutations(10).random_element()
        sage: G = BinaryWord(pi)
        sage: list(BinaryWord(labels=G.out_labels())) == list(G)
        True

    Test that the Kleitman Greene invariant is indeed the descent word::

        sage: r = 4
        sage: all(Word([0 if i in w.descents() else 1 for i in range(r)])
        ....:      == BinaryWord(w).out_labels()[r]
        ....:     for w in Permutations(r))
        True
    """
    zero = Word([], alphabet=[0,1])

    def normalize_vertex(self, v):
        r"""
        Return ``v`` as a binary word.

        EXAMPLES::

            sage: BinaryWord = GrowthDiagram.rules.BinaryWord()
            sage: BinaryWord.normalize_vertex([0,1]).parent()
            Finite words over {0, 1}
        """
        return Word(v, alphabet=[0,1])

    def vertices(self, n):
        r"""
        Return the vertices of the dual graded graph on level ``n``.

        EXAMPLES::

            sage: BinaryWord = GrowthDiagram.rules.BinaryWord()
            sage: BinaryWord.vertices(3)
            [word: 100, word: 101, word: 110, word: 111]
        """
        if n == 0:
            return [self.zero]
        else:
            w1 = Word([1], [0,1])
            return [w1 + w for w in Words([0,1], n-1)]

    def rank(self, v):
        r"""
        Return the rank of ``v``: number of letters of the word.

        EXAMPLES::

            sage: BinaryWord = GrowthDiagram.rules.BinaryWord()
            sage: BinaryWord.rank(BinaryWord.vertices(3)[0])
            3
        """
        return len(v)

    def is_Q_edge(self, v, w):
        r"""
        Return whether ``(v, w)`` is a `Q`-edge of ``self``.

        ``(w, v)`` is an edge if ``w`` is obtained from ``v`` by
        appending a letter.

        EXAMPLES::

            sage: BinaryWord = GrowthDiagram.rules.BinaryWord()
            sage: v = BinaryWord.vertices(2)[0]; v
            word: 10
            sage: [w for w in BinaryWord.vertices(3) if BinaryWord.is_Q_edge(v, w)]
            [word: 100, word: 101]
            sage: [w for w in BinaryWord.vertices(4) if BinaryWord.is_Q_edge(v, w)]
            []
        """
        return w[:-1] == v

    def is_P_edge(self, v, w):
        r"""
        Return whether ``(v, w)`` is a `P`-edge of ``self``.

        ``(v, w)`` is an edge if ``v`` is obtained from ``w`` by
        deleting a letter.

        EXAMPLES::

            sage: BinaryWord = GrowthDiagram.rules.BinaryWord()
            sage: v = BinaryWord.vertices(2)[1]; v
            word: 11
            sage: [w for w in BinaryWord.vertices(3) if BinaryWord.is_P_edge(v, w)]
            [word: 101, word: 110, word: 111]
            sage: [w for w in BinaryWord.vertices(4) if BinaryWord.is_P_edge(v, w)]
            []
        """
        return len(w) == len(v) + 1 and v.is_subword_of(w)

    def forward_rule(self, y, t, x, content):
        r"""
        Return the output shape given three shapes and the content.

        See [Fom1995]_ Lemma 4.6.1, page 40.

        INPUT:

        - ``y, t, x`` -- three binary words from a cell in a growth
          diagram, labelled as::

              t x
              y

        - ``content`` -- `0` or `1`; the content of the cell

        OUTPUT:

        The fourth binary word ``z`` according to Viennot's
        bijection [Vie1983]_.

        EXAMPLES::

            sage: BinaryWord = GrowthDiagram.rules.BinaryWord()

            sage: BinaryWord.forward_rule([], [], [], 1)
            word: 1

            sage: BinaryWord.forward_rule([1], [1], [1], 1)
            word: 11

        if ``x != y`` append last letter of ``x`` to ``y``::

            sage: BinaryWord.forward_rule([1,0], [1], [1,1], 0)
            word: 101

        if ``x == y != t`` append ``0`` to ``y``::

            sage: BinaryWord.forward_rule([1,1], [1], [1,1], 0)
            word: 110
        """
        if x == t == y:
            if content == 0:
                z = x
            elif content == 1:
                z = Word(list(y) + [1], alphabet=[0,1])
            else:
                raise NotImplementedError
        elif content != 0:
            raise ValueError("for y=%s, t=%s, x=%s, the content should be 0 but is %s"
                             % (y, t, x, content))
        elif x != t == y:
            z = x
        elif x == t != y:
            z = y
        else:
            if x != y:
                z = Word(list(y) + [x[-1]], alphabet=[0,1])
            elif x == y != t:
                z = Word(list(y) + [0], alphabet=[0,1])
            else:
                raise NotImplementedError
        return z

    def backward_rule(self, y, z, x):
        r"""
        Return the content and the input shape.

        See [Fom1995]_ Lemma 4.6.1, page 40.

        - ``y, z, x`` -- three binary words from a cell in a growth diagram,
          labelled as::

                x
              y z

        OUTPUT:

        A pair ``(t, content)`` consisting of the shape of the fourth
        word and the content of the cell according to Viennot's
        bijection [Vie1983]_.

        TESTS::

            sage: BinaryWord = GrowthDiagram.rules.BinaryWord()
            sage: w = [4,1,8,3,6,5,2,7,9]; G = GrowthDiagram(BinaryWord, w)
            sage: BinaryWord(labels=G.out_labels()).to_word() == w  # indirect doctest
            True
        """
        if x == y == z:
            return (x, 0)
        elif x == z != y:
            return (y, 0)
        elif x != z == y:
            return (x, 0)
        else:
            if x == y and len(z) > 0 and z[-1] == 1:
                return (x, 1)
            else:
                return (x[:-1], 0)

class RuleSylvester(Rule):
    r"""
    A rule modelling a Schensted-like correspondence for binary trees.

    EXAMPLES::

        sage: Sylvester = GrowthDiagram.rules.Sylvester()
        sage: GrowthDiagram(Sylvester, [3,1,2])
        0  1  0
        0  0  1
        1  0  0

    The vertices of the dual graded graph are
    :class:`~sage.combinat.binary_tree.BinaryTrees`::

        sage: Sylvester.vertices(3)
        Binary trees of size 3

    The :meth:`~sage.combinat.growth.Rule.P_graph` is also known as
    the bracket tree, the :meth:`~sage.combinat.growth.Rule.Q_graph`
    is the lattice of finite order ideals of the infinite binary
    tree, see Example 2.4.6 in [Fom1994]_.

    For a permutation, the :meth:`P_symbol` is the binary search
    tree, the :meth:`Q_symbol` is the increasing tree corresponding
    to the inverse permutation.  Note that, instead of passing the
    rule to :class:`GrowthDiagram`, we can also call the rule to
    create growth diagrams.  From [Nze2007]_::

        sage: pi = Permutation([3,5,1,4,2,6]); G = Sylvester(pi); G
        0  0  1  0  0  0
        0  0  0  0  1  0
        1  0  0  0  0  0
        0  0  0  1  0  0
        0  1  0  0  0  0
        0  0  0  0  0  1
        sage: ascii_art(G.P_symbol())
          __3__
         /     \
        1       5
         \     / \
          2   4   6
        sage: ascii_art(G.Q_symbol())
          __1__
         /     \
        3       2
         \     / \
          5   4   6

        sage: all(Sylvester(pi).P_symbol() == pi.binary_search_tree()
        ....:     for pi in Permutations(5))
        True

        sage: all(Sylvester(pi).Q_symbol() == pi.inverse().increasing_tree()
        ....:     for pi in Permutations(5))
        True

    TESTS::

        sage: Sylvester.zero
        .

        sage: B = BinaryTree; R = B([None,[]]); L = B([[],None])
        sage: T = B([[],[]]); S = B([L,None])
        sage: G = Sylvester(labels=[R, T, R])
        Traceback (most recent call last):
        ...
        ValueError: [., [., .]] has smaller rank than [[., .], [., .]]
         but is not covered by it in P

        sage: G = Sylvester(labels=[R, S, R])
        Traceback (most recent call last):
        ...
        ValueError: [., [., .]] has smaller rank than [[[., .], .], .]
         but is not covered by it in Q

    Check duality::

        sage: Sylvester._check_duality(4)

    Check that the rules are bijective::

        sage: all(Sylvester(labels=GrowthDiagram(Sylvester, pi).out_labels()).to_word()
        ....:      == pi for pi in Permutations(4))
        True
        sage: pi = Permutations(10).random_element()
        sage: G = GrowthDiagram(Sylvester, pi)
        sage: list(Sylvester(labels=G.out_labels())) == list(G)
        True
    """
    zero = BinaryTree()  # type: ignore

    def normalize_vertex(self, v):
        r"""
        Return ``v`` as a binary tree.

        EXAMPLES::

            sage: Sylvester = GrowthDiagram.rules.Sylvester()
            sage: Sylvester.normalize_vertex([[],[]]).parent()
            Binary trees
        """
        return BinaryTree(v)

    def vertices(self, n):
        r"""
        Return the vertices of the dual graded graph on level ``n``.

        EXAMPLES::

            sage: Sylvester = GrowthDiagram.rules.Sylvester()
            sage: Sylvester.vertices(3)
            Binary trees of size 3
        """
        return BinaryTrees(n)

    def rank(self, v):
        r"""
        Return the rank of ``v``: the number of nodes of the tree.

        EXAMPLES::

            sage: Sylvester = GrowthDiagram.rules.Sylvester()
            sage: Sylvester.rank(Sylvester.vertices(3)[0])
            3
        """
        return v.node_number()

    def is_Q_edge(self, v, w):
        r"""
        Return whether ``(v, w)`` is a `Q`-edge of ``self``.

        ``(v, w)`` is an edge if ``v`` is a sub-tree of ``w`` with one
        node less.

        EXAMPLES::

            sage: Sylvester = GrowthDiagram.rules.Sylvester()
            sage: v = Sylvester.vertices(2)[1]; ascii_art(v)
              o
             /
            o
            sage: ascii_art([w for w in Sylvester.vertices(3) if Sylvester.is_Q_edge(v, w)])
            [   o  ,   o,     o ]
            [  / \    /      /  ]
            [ o   o  o      o   ]
            [         \    /    ]
            [          o  o     ]
            sage: [w for w in Sylvester.vertices(4) if Sylvester.is_Q_edge(v, w)]
            []
        """
        def is_subtree(T1, T2):
            if T2.is_empty():
                return False
            elif T2[0].is_empty() and T2[1].is_empty():
                return T1.is_empty()
            elif T1.is_empty():
                return False
            else:
                return ((T1[0] == T2[0] and is_subtree(T1[1], T2[1])) or
                        (T1[1] == T2[1] and is_subtree(T1[0], T2[0])))
        return is_subtree(v, w)

    def is_P_edge(self, v, w):
        r"""
        Return whether ``(v, w)`` is a `P`-edge of ``self``.

        ``(v, w)`` is an edge if ``v`` is obtained from ``w`` by deleting
        its right-most node.

        EXAMPLES::

            sage: Sylvester = GrowthDiagram.rules.Sylvester()
            sage: v = Sylvester.vertices(2)[1]; ascii_art(v)
              o
             /
            o

            sage: ascii_art([w for w in Sylvester.vertices(3) if Sylvester.is_P_edge(v, w)])
            [   o  ,     o ]
            [  / \      /  ]
            [ o   o    o   ]
            [         /    ]
            [        o     ]

            sage: [w for w in Sylvester.vertices(4) if Sylvester.is_P_edge(v, w)]
            []
        """
        if w.is_empty():
            return False
        else:
            return v == RuleSylvester._delete_right_most_node(w)

    def P_symbol(self, P_chain):
        r"""
        Return the labels along the vertical boundary of a rectangular
        growth diagram as a labelled binary tree.

        For permutations, this coincides with the binary search tree.

        EXAMPLES::

            sage: Sylvester = GrowthDiagram.rules.Sylvester()
            sage: pi = Permutation([2,4,3,1])
            sage: ascii_art(Sylvester(pi).P_symbol())
              _2_
             /   \
            1     4
                 /
                3
            sage: Sylvester(pi).P_symbol() == pi.binary_search_tree()
            True

        We can also do the skew version::

            sage: B = BinaryTree; E = B(); N = B([])
            sage: ascii_art(Sylvester([3,2], shape=[3,3,3], labels=[N,N,N,E,E,E,N]).P_symbol())
              __1___
             /      \
            None     3
                    /
                   2
        """
        def add_label(L, S, T, m):
            if T[0] == S:
                L = LabelledBinaryTree([L, None], m)
            else:
                assert T[0] == S[0]
                l = L.label()
                L = LabelledBinaryTree([L[0], add_label(L[1], S[1], T[1], m)], l)
            return L

        L = LabelledBinaryTree(P_chain[0])
        for i in range(1, len(P_chain)):
            S, T = P_chain[i-1], P_chain[i]
            L = add_label(L, S, T, i)
        return L

    def Q_symbol(self, Q_chain):
        r"""
        Return the labels along the vertical boundary of a rectangular
        growth diagram as a labelled binary tree.

        For permutations, this coincides with the increasing tree.

        EXAMPLES::

            sage: Sylvester = GrowthDiagram.rules.Sylvester()
            sage: pi = Permutation([2,4,3,1])
            sage: ascii_art(Sylvester(pi).Q_symbol())
              _1_
             /   \
            4     2
                 /
                3
            sage: Sylvester(pi).Q_symbol() == pi.inverse().increasing_tree()
            True

        We can also do the skew version::

            sage: B = BinaryTree; E = B(); N = B([])
            sage: ascii_art(Sylvester([3,2], shape=[3,3,3], labels=[N,N,N,E,E,E,N]).Q_symbol())
              _None_
             /      \
            3        1
                    /
                   2
        """
        def add_label(L, S, T, m):
            if L.is_empty():
                assert T.node_number() == 1
                return LabelledBinaryTree([], m)
            l = L.label()
            if T[0] == S[0]:
                return LabelledBinaryTree([L[0], add_label(L[1], S[1], T[1], m)], l)
            else:
                return LabelledBinaryTree([add_label(L[0], S[0], T[0], m), L[1]], l)

        L = LabelledBinaryTree(Q_chain[0])
        for i in range(1, len(Q_chain)):
            S, T = Q_chain[i-1], Q_chain[i]
            L = add_label(L, S, T, i)
        return L


    @staticmethod
    def _delete_right_most_node(b):
        r"""
        Return the tree obtained by deleting the right most node from ``b``.

        TESTS::

            sage: Sylvester = GrowthDiagram.rules.Sylvester()
            sage: b = BinaryTree([]); b
            [., .]
            sage: Sylvester._delete_right_most_node(b)
            .
            sage: b = BinaryTree([[[], []], None]); ascii_art(b)
                o
               /
              o
             / \
            o   o
            sage: ascii_art(Sylvester._delete_right_most_node(b))
              o
             / \
            o   o
        """
        if b.is_empty():
            raise ValueError("cannot delete right most node from empty tree")
        elif b[1].is_empty():
            return b[0]
        else:
            return BinaryTree([b[0], RuleSylvester._delete_right_most_node(b[1])])

    def forward_rule(self, y, t, x, content):
        r"""
        Return the output shape given three shapes and the content.

        See [Nze2007]_, page 9.

        INPUT:

        - ``y, t, x`` -- three binary trees from a cell in a growth
          diagram, labelled as::

              t x
              y

        - ``content`` -- `0` or `1`; the content of the cell

        OUTPUT:

        The fourth binary tree ``z``.

        EXAMPLES::

            sage: Sylvester = GrowthDiagram.rules.Sylvester()
            sage: B = BinaryTree; E = B(); N = B([]); L = B([[],None])
            sage: R = B([None,[]]); T = B([[],[]])

            sage: ascii_art(Sylvester.forward_rule(E, E, E, 1))
            o
            sage: ascii_art(Sylvester.forward_rule(N, N, N, 1))
            o
             \
              o
            sage: ascii_art(Sylvester.forward_rule(L, L, L, 1))
              o
             / \
            o   o
            sage: ascii_art(Sylvester.forward_rule(R, R, R, 1))
            o
             \
              o
               \
                o

        If ``y != x``, obtain ``z`` from ``y`` adding a node such
        that deleting the right most gives ``x``::

            sage: ascii_art(Sylvester.forward_rule(R, N, L, 0))
              o
             / \
            o   o

            sage: ascii_art(Sylvester.forward_rule(L, N, R, 0))
              o
             /
            o
             \
              o

        If ``y == x != t``, obtain ``z`` from ``x`` by adding a node
        as left child to the right most node::

            sage: ascii_art(Sylvester.forward_rule(N, E, N, 0))
              o
             /
            o
            sage: ascii_art(Sylvester.forward_rule(T, L, T, 0))
              _o_
             /   \
            o     o
                 /
                o
            sage: ascii_art(Sylvester.forward_rule(L, N, L, 0))
                o
               /
              o
             /
            o
            sage: ascii_art(Sylvester.forward_rule(R, N, R, 0))
            o
             \
              o
             /
            o
        """
        def successors(b):
            r"""
            Return all trees obtained from ``b`` by adding a node.
            """
            if b.is_empty():
                yield BinaryTree([])
            else:
                for t in successors(b[0]):
                    yield BinaryTree([t, b[1]])
                for t in successors(b[1]):
                    yield BinaryTree([b[0], t])

        def union(y, x):
            r"""
            Return the unique tree obtained by adding a node to ``y`` such
            that deleting the right most node gives ``x``.
            """
            for t in successors(y):
                if RuleSylvester._delete_right_most_node(t) == x:
                    return t
            raise ValueError("could not find union of %s and %s" % (y,x))

        if y == t == x:
            if content == 0:
                z = y
            elif content == 1:
                z = t.over(BinaryTree([]))
            else:
                raise NotImplementedError
        elif content != 0:
            raise ValueError("for y=%s, t=%s, x=%s, the content should be 0 but is %s" % (y, t, x, content))
        elif y != t == x:
            z = y
        elif y == t != x:
            z = x
        else:
            z = union(y, x)
        return z

    def backward_rule(self, y, z, x):
        r"""
        Return the output shape given three shapes and the content.

        See [Nze2007]_, page 9.

        INPUT:

        - ``y, z, x`` -- three binary trees from a cell in a growth
          diagram, labelled as::

                x
              y z

        OUTPUT:

        A pair ``(t, content)`` consisting of the shape of the fourth
        binary tree t and the content of the cell.

        EXAMPLES::

            sage: Sylvester = GrowthDiagram.rules.Sylvester()
            sage: B = BinaryTree; E = B(); N = B([]); L = B([[],None])
            sage: R = B([None,[]]); T = B([[],[]])

            sage: ascii_art(Sylvester.backward_rule(E, E, E))
            ( , 0 )
            sage: ascii_art(Sylvester.backward_rule(N, N, N))
            ( o, 0 )
        """
        if x == y == z:
            return (x, 0)
        elif x == z != y:
            return (y, 0)
        elif x != z == y:
            return (x, 0)
        else:
            if x == y and z == x.over(BinaryTree([])):
                return (x, 1)
            else:
                t = RuleSylvester._delete_right_most_node(y)
                return (t, 0)

class RuleYoungFibonacci(Rule):
    r"""
    A rule modelling a Schensted-like correspondence for
    Young-Fibonacci-tableaux.

    EXAMPLES::

        sage: YF = GrowthDiagram.rules.YoungFibonacci()
        sage: GrowthDiagram(YF, [3,1,2])
        0  1  0
        0  0  1
        1  0  0

    The vertices of the dual graded graph are Fibonacci words -
    compositions into parts of size at most two::

        sage: YF.vertices(4)
        [word: 22, word: 211, word: 121, word: 112, word: 1111]

    Note that, instead of passing the rule to :class:`GrowthDiagram`,
    we can also use call the rule to create growth diagrams.  For
    example::

        sage: G = YF([2, 3, 7, 4, 1, 6, 5]); G
        0  0  0  0  1  0  0
        1  0  0  0  0  0  0
        0  1  0  0  0  0  0
        0  0  0  1  0  0  0
        0  0  0  0  0  0  1
        0  0  0  0  0  1  0
        0  0  1  0  0  0  0

    The Kleitman Greene invariant is: take the last letter and the
    largest letter of the permutation and remove them.  If they
    coincide write 1, otherwise write 2::

        sage: G.P_chain()[-1]
        word: 21211

    TESTS::

        sage: YF = GrowthDiagram.rules.YoungFibonacci()
        sage: YF.zero
        word:

    Check duality::

        sage: YF._check_duality(4)

        sage: G = YF(labels=[[1],[1,0],[1]])
        Traceback (most recent call last):
        ...
        ValueError: 0 not in alphabet!

        sage: G = YF(labels=[[1,1],[1,2]])
        Traceback (most recent call last):
        ...
        ValueError: 11 has smaller rank than 12 but is not covered by it in Q

        sage: G = YF(labels=[[1,2],[1,1]])
        Traceback (most recent call last):
        ...
        ValueError: 11 has smaller rank than 12 but is not covered by it in P

        sage: all(YF(labels=YF(pi).out_labels()).to_word()
        ....:      == pi for pi in Permutations(4))
        True
        sage: pi = Permutations(10).random_element()
        sage: G = YF(pi)
        sage: list(YF(labels=G.out_labels())) == list(G)
        True
    """
    zero = Word([], alphabet=[1,2])

    def normalize_vertex(self, v):
        r"""
        Return ``v`` as a word with letters 1 and 2.

        EXAMPLES::

            sage: YF = GrowthDiagram.rules.YoungFibonacci()
            sage: YF.normalize_vertex([1,2,1]).parent()
            Finite words over {1, 2}
        """
        return Word(v, alphabet=[1,2])

    def vertices(self, n):
        r"""
        Return the vertices of the dual graded graph on level ``n``.

        EXAMPLES::

            sage: YF = GrowthDiagram.rules.YoungFibonacci()
            sage: YF.vertices(3)
            [word: 21, word: 12, word: 111]
        """
        if n == 0:
            return [self.zero]
        else:
            return [Word(list(w), [1,2]) for w in Compositions(n, max_part=2)]

    def rank(self, v):
        r"""
        Return the rank of ``v``: the size of the corresponding composition.

        EXAMPLES::

            sage: YF = GrowthDiagram.rules.YoungFibonacci()
            sage: YF.rank(YF.vertices(3)[0])
            3
        """
        return sum(v)

    def is_P_edge(self, v, w):
        r"""
        Return whether ``(v, w)`` is a `P`-edge of ``self``.

        ``(v, w)`` is an edge if ``v`` is obtained from ``w`` by deleting
        a ``1`` or replacing the left-most ``2`` by a ``1``.

        EXAMPLES::

            sage: YF = GrowthDiagram.rules.YoungFibonacci()
            sage: v = YF.vertices(5)[5]; v
            word: 1121
            sage: [w for w in YF.vertices(6) if YF.is_P_edge(v, w)]
            [word: 2121, word: 11121]
            sage: [w for w in YF.vertices(7) if YF.is_P_edge(v, w)]
            []
        """
        if sum(w) != sum(v) + 1:
            return False
        ell = len(v)
        w = list(w)
        for i in range(ell+1):
            d = list(v)
            d.insert(i, 1)
            if w == d:
                return True
            if i < ell and v[i] == 1:
                d = list(v)
                d[i] = 2
                if w == d:
                    return True
                break
        return False

    is_Q_edge = is_P_edge

    def forward_rule(self, y, t, x, content):
        r"""
        Return the output shape given three shapes and the content.

        See [Fom1995]_ Lemma 4.4.1, page 35.

        INPUT:

        - ``y, t, x`` -- three Fibonacci words from a
          cell in a growth diagram, labelled as::

              t x
              y

        - ``content`` -- `0` or `1`; the content of the cell

        OUTPUT:

        The fourth Fibonacci word.

        EXAMPLES::

            sage: YF = GrowthDiagram.rules.YoungFibonacci()

            sage: YF.forward_rule([], [], [], 1)
            word: 1

            sage: YF.forward_rule([1], [1], [1], 1)
            word: 11

            sage: YF.forward_rule([1,2], [1], [1,1], 0)
            word: 21

            sage: YF.forward_rule([1,1], [1], [1,1], 0)
            word: 21
        """
        if x == t == y:
            if content == 0:
                r = x
            elif content == 1:
                r = Word([1] + list(y), alphabet=[1,2])
            else:
                raise NotImplementedError
        elif content != 0:
            raise ValueError("for y=%s, t=%s, x=%s, the content should be 0 but is %s"
                             % (y, t, x, content))
        elif x == t:
            r = y
        elif y == t:
            r = x
        else:
            if x != t != y:
                r = Word([2] + list(t), alphabet=[1,2])
            else:
                raise NotImplementedError("for y=%s, t=%s, x=%s, content %s we have no rule"
                                          % (y, t, x, content))
        return r

    def backward_rule(self, y, z, x):
        r"""
        Return the content and the input shape.

        See [Fom1995]_ Lemma 4.4.1, page 35.

        - ``y, z, x`` -- three Fibonacci words from a cell in a
          growth diagram, labelled as::

                x
              y z

        OUTPUT:

        A pair ``(t, content)`` consisting of the shape of the fourth
        word and the content of the cell.

        TESTS::

            sage: YF = GrowthDiagram.rules.YoungFibonacci()
            sage: w = [4,1,8,3,6,5,2,7,9]; G = YF(w)
            sage: GrowthDiagram(YF, labels=G.out_labels()).to_word() == w  # indirect doctest
            True
        """
        if x == y == z:
            return (x, 0)
        elif x == z != y:
            return (y, 0)
        elif x != z == y:
            return (x, 0)
        else:
            if z[0] == 1:
                return (z[1:], 1)
            elif z[0] == 2:
                return (z[1:], 0)

class RulePartitions(Rule):
    r"""
    A rule for growth diagrams on Young's lattice on integer
    partitions graded by size.

    TESTS::

        sage: Burge = GrowthDiagram.rules.Burge()
        sage: Burge.zero
        []
        sage: G = GrowthDiagram(Burge, labels=[[1],[1]])
        Traceback (most recent call last):
        ...
        ValueError: can only determine the shape of the growth diagram
         if ranks of successive labels differ
    """
    zero = _make_partition([])

    def vertices(self, n):
        r"""
        Return the vertices of the dual graded graph on level ``n``.

        EXAMPLES::

            sage: RSK = GrowthDiagram.rules.RSK()
            sage: RSK.vertices(3)
            Partitions of the integer 3
        """
        return Partitions(n)

    def normalize_vertex(self, v):
        r"""
        Return ``v`` as a partition.

        EXAMPLES::

            sage: RSK = GrowthDiagram.rules.RSK()
            sage: RSK.normalize_vertex([3,1]).parent()
            Partitions
        """
        return _make_partition(v)

    def rank(self, v):
        r"""
        Return the rank of ``v``: the size of the partition.

        EXAMPLES::

            sage: RSK = GrowthDiagram.rules.RSK()
            sage: RSK.rank(RSK.vertices(3)[0])
            3
        """
        return v.size()

    def P_symbol(self, P_chain):
        r"""
        Return the labels along the vertical boundary of a rectangular
        growth diagram as a (skew) tableau.

        EXAMPLES::

            sage: RuleRSK = GrowthDiagram.rules.RSK()
            sage: G = RuleRSK([[0,1,0], [1,0,2]])
            sage: G.P_symbol().pp()
            1  2  2
            2
        """
        return SkewTableau(chain=P_chain)

    def Q_symbol(self, Q_chain):
        r"""
        Return the labels along the horizontal boundary of a rectangular
        growth diagram as a skew tableau.

        EXAMPLES::

            sage: RuleRSK = GrowthDiagram.rules.RSK()
            sage: G = RuleRSK([[0,1,0], [1,0,2]])
            sage: G.Q_symbol().pp()
            1  3  3
            2
        """
        return SkewTableau(chain=Q_chain)

class RuleRSK(RulePartitions):
    r"""
    A rule modelling Robinson-Schensted-Knuth insertion.

    EXAMPLES::

        sage: RuleRSK = GrowthDiagram.rules.RSK()
        sage: GrowthDiagram(RuleRSK, [3,2,1,2,3])
        0  0  1  0  0
        0  1  0  1  0
        1  0  0  0  1

    The vertices of the dual graded graph are integer partitions::

        sage: RuleRSK.vertices(3)
        Partitions of the integer 3

    The local rules implemented provide the RSK correspondence
    between matrices with non-negative integer entries and pairs of
    semistandard tableaux, the
    :meth:`~sage.combinat.growth.RulePartitions.P_symbol` and the
    :meth:`~sage.combinat.growth.RulePartitions.Q_symbol`.  For
    permutations, it reduces to classical Schensted insertion.

    Instead of passing the rule to :class:`GrowthDiagram`, we can
    also call the rule to create growth diagrams.  For example::

        sage: m = matrix([[0,0,0,0,1],[1,1,0,2,0], [0,3,0,0,0]])
        sage: G = RuleRSK(m); G
        0  0  0  0  1
        1  1  0  2  0
        0  3  0  0  0

        sage: ascii_art([G.P_symbol(), G.Q_symbol()])
        [   1  2  2  2  3    1  2  2  2  2 ]
        [   2  3             4  4          ]
        [   3            ,   5             ]

    For rectangular fillings, the Kleitman-Greene invariant is the
    shape of the :meth:`P_symbol` (or the :meth:`Q_symbol`).  Put
    differently, it is the partition labelling the lower right corner
    of the filling (recall that we are using matrix coordinates).  It
    can be computed alternatively as the partition
    `(\mu_1,\dots,\mu_n)`, where `\mu_1 + \dots + \mu_i` is the
    maximal sum of entries in a collection of `i` pairwise disjoint
    sequences of cells with weakly increasing coordinates.

    For rectangular fillings, we could also use the (faster)
    implementation provided via :func:`~sage.combinat.rsk.RSK`.
    Because the of the coordinate conventions in
    :func:`~sage.combinat.rsk.RSK`, we have to transpose matrices::

        sage: [G.P_symbol(), G.Q_symbol()] == RSK(m.transpose())
        True

        sage: n=5; l=[(pi, RuleRSK(pi)) for pi in Permutations(n)]
        sage: all([G.P_symbol(), G.Q_symbol()] == RSK(pi) for pi, G in l)
        True

        sage: n=5; l=[(w, RuleRSK(w)) for w in Words([1,2,3], 5)]
        sage: all([G.P_symbol(), G.Q_symbol()] == RSK(pi) for pi, G in l)
        True
    """
    def forward_rule(self, y, t, x, content):
        r"""
        Return the output shape given three shapes and the content.

        See [Kra2006]_ `(F^1 0)-(F^1 2)`.

        INPUT:

        - ``y, t, x`` -- three partitions from a cell in a
          growth diagram, labelled as::

              t x
              y

        - ``content`` -- a non-negative integer; the content of the cell

        OUTPUT:

        The fourth partition according to the Robinson-Schensted-Knuth
        correspondence.

        EXAMPLES::

            sage: RuleRSK = GrowthDiagram.rules.RSK()
            sage: RuleRSK.forward_rule([2,1],[2,1],[2,1],1)
            [3, 1]

            sage: RuleRSK.forward_rule([1],[],[2],2)
            [4, 1]
        """
        carry = content
        z = []
        while True:
            if x == []:
                row1 = 0
            else:
                row1 = x[0]
            if y == []:
                row3 = 0
            else:
                row3 = y[0]
            newPart = max(row1, row3) + carry
            if newPart == 0:
                # returning this as a Partition costs a lot of time
                return z[::-1]
            else:
                z = [newPart] + z
                if t == []:
                    carry = min(row1, row3)
                else:
                    carry = min(row1, row3) - t[0]
                x = x[1:]
                t = t[1:]
                y = y[1:]

    def backward_rule(self, y, z, x):
        r"""
        Return the content and the input shape.

        See [Kra2006]_ `(B^1 0)-(B^1 2)`.

        INPUT:

        - ``y, z, x`` -- three partitions from a cell in a
          growth diagram, labelled as::

              x
            y z

        OUTPUT:

        A pair ``(t, content)`` consisting of the shape of the fourth
        word according to the Robinson-Schensted-Knuth correspondence
        and the content of the cell.

        TESTS::

            sage: RuleRSK = GrowthDiagram.rules.RSK()
            sage: w = [4,1,8,3,6,5,2,7,9]; G = RuleRSK(w)
            sage: GrowthDiagram(RuleRSK, labels=G._out_labels).to_word() == w  # indirect doctest
            True
        """
        carry = 0
        i = len(z)
        t = []
        while i > 0:
            if len(x) < i:
                row1 = 0
            else:
                row1 = x[i-1]
            if len(y) < i:
                row3 = 0
            else:
                row3 = y[i-1]
            t = [min(row1, row3) - carry] + t
            carry = z[i-1] - max(row1, row3)
            i = i-1
        return (_make_partition(t), carry)


class RuleBurge(RulePartitions):
    r"""
    A rule modelling Burge insertion.

    EXAMPLES::

        sage: Burge = GrowthDiagram.rules.Burge()
        sage: GrowthDiagram(Burge, labels=[[],[1,1,1],[2,1,1,1],[2,1,1],[2,1],[1,1],[]])
        1  1
        0  1
        1  0
        1  0

    The vertices of the dual graded graph are integer partitions::

        sage: Burge.vertices(3)
        Partitions of the integer 3

    The local rules implemented provide Burge's correspondence
    between matrices with non-negative integer entries and pairs of
    semistandard tableaux, the
    :meth:`~sage.combinat.growth.RulePartitions.P_symbol` and the
    :meth:`~sage.combinat.growth.RulePartitions.Q_symbol`.  For
    permutations, it reduces to classical Schensted insertion.

    Instead of passing the rule to :class:`GrowthDiagram`, we can
    also call the rule to create growth diagrams.  For example::

        sage: m = matrix([[2,0,0,1,0],[1,1,0,0,0], [0,0,0,0,3]])
        sage: G = Burge(m); G
        2  0  0  1  0
        1  1  0  0  0
        0  0  0  0  3

        sage: ascii_art([G.P_symbol(), G.Q_symbol()])
        [   1  2  3    1  2  5 ]
        [   1  3       1  5    ]
        [   1  3       1  5    ]
        [   2      ,   4       ]

    For rectangular fillings, the Kleitman-Greene invariant is the
    shape of the
    :meth:`~sage.combinat.growth.RulePartitions.P_symbol`.  Put
    differently, it is the partition labelling the lower right corner
    of the filling (recall that we are using matrix coordinates).  It
    can be computed alternatively as the transpose of the partition
    `(\mu_1, \ldots, \mu_n)`, where `\mu_1 + \cdots + \mu_i` is the
    maximal sum of entries in a collection of `i` pairwise disjoint
    sequences of cells with weakly decreasing row indices and weakly
    increasing column indices.
    """
    def forward_rule(self, y, t, x, content):
        r"""
        Return the output shape given three shapes and the content.

        See [Kra2006]_ `(F^4 0)-(F^4 2)`.

        INPUT:

        - ``y, t, x`` -- three  from a cell in a growth diagram,
          labelled as::

              t x
              y

        - ``content`` -- a non-negative integer; the content of the cell

        OUTPUT:

        The fourth partition according to the Burge correspondence.

        EXAMPLES::

            sage: Burge = GrowthDiagram.rules.Burge()
            sage: Burge.forward_rule([2,1],[2,1],[2,1],1)
            [3, 1]

            sage: Burge.forward_rule([1],[],[2],2)
            [2, 1, 1, 1]
        """
        # n is the maximal length of longest decreasing chain by
        # Kleitman-Greene's theorem
        n = content + len(x) + len(y)
        x += [0]*(n-len(x))
        y += [0]*(n-len(y))
        t += [0]*(n-len(t))
        z = [0]*n
        carry = content
        for i, (row1, row2, row3) in enumerate(zip(x, t, y)):
            s = min(int(row1 == row2 == row3), carry)
            new_part = max(row1, row3) + s
            if new_part:
                z[i] = new_part
                carry += -s + min(row1, row3) - row2
            else:
                break
        return _make_partition(z)

    def backward_rule(self, y, z, x):
        r"""
        Return the content and the input shape.

        See [Kra2006]_ `(B^4 0)-(B^4 2)`.  (In the arXiv version of
        the article there is a typo: in the computation of carry in
        `(B^4 2)` , `\rho` must be replaced by `\lambda`).

        INPUT:

        - ``y, z, x`` -- three partitions from a cell in a
          growth diagram, labelled as::

              x
            y z

        OUTPUT:

        A pair ``(t, content)`` consisting of the shape of the fourth
        partition according to the Burge correspondence and the content of
        the cell.

        EXAMPLES::

            sage: Burge = GrowthDiagram.rules.Burge()
            sage: Burge.backward_rule([1,1,1],[2,1,1,1],[2,1,1])
            ([1, 1], 0)

        TESTS::

            sage: w = [4,1,8,3,6,5,2,7,9]; G = Burge(w)
            sage: GrowthDiagram(Burge, labels=G._out_labels).to_word() == w  # indirect doctest
            True
        """
        t = [0]*len(z) # z must be the longest partition
        mu = [0]*(len(z)-len(x)) + x[::-1]
        nu = [0]*(len(z)-len(y)) + y[::-1]
        la = z[::-1]
        carry = 0
        for i, (mu_i, la_i, nu_i) in enumerate(zip(mu, la, nu)):
            s = min(int(mu_i == nu_i == la_i), carry)
            t[i] = min(mu_i, nu_i) - s
            carry += -s + la_i - max(mu_i, nu_i)
        t.reverse()
        return (_make_partition(t), carry)

class RuleDomino(Rule):
    r"""
    A rule modelling domino insertion.

    EXAMPLES::

        sage: Domino = GrowthDiagram.rules.Domino()
        sage: GrowthDiagram(Domino, [[1,0,0],[0,0,1],[0,-1,0]])
        1  0  0
        0  0  1
        0 -1  0

    The vertices of the dual graded graph are integer partitions
    whose Ferrers diagram can be tiled with dominoes::

        sage: Domino.vertices(2)
        [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]

    Instead of passing the rule to :class:`GrowthDiagram`, we can
    also call the rule to create growth diagrams.  For example, let
    us check Figure 3 in [Lam2004]_::

        sage: G = Domino([[0,0,0,-1],[0,0,1,0],[-1,0,0,0],[0,1,0,0]]); G
         0  0  0 -1
         0  0  1  0
        -1  0  0  0
         0  1  0  0

        sage: ascii_art([G.P_symbol(), G.Q_symbol()])
        [   1  2  4    1  2  2 ]
        [   1  2  4    1  3  3 ]
        [   3  3   ,   4  4    ]

    The spin of a domino tableau is half the number of vertical dominoes::

        sage: def spin(T):
        ....:     return sum(2*len(set(row)) - len(row) for row in T)/4

    According to [Lam2004]_, the number of negative entries in the
    signed permutation equals the sum of the spins of the two
    associated tableaux::

        sage: pi = [3,-1,2,4,-5]
        sage: G = Domino(pi)
        sage: list(G.filling().values()).count(-1) == spin(G.P_symbol()) + spin(G.Q_symbol())
        True

    Negating all signs transposes all the partitions::

        sage: G.P_symbol() == Domino([-e for e in pi]).P_symbol().conjugate()
        True

    TESTS:

    Check duality::

        sage: Domino = GrowthDiagram.rules.Domino()
        sage: Domino._check_duality(3)

        sage: G = Domino([[0,1,0],[0,0,-1],[1,0,0]]); G
        0  1  0
        0  0 -1
        1  0  0

        sage: ascii_art([G.P_symbol(), G.Q_symbol()])
        [   1  1    1  1 ]
        [   2  3    2  2 ]
        [   2  3,   3  3 ]

        sage: l = {pi: Domino(pi) for pi in SignedPermutations(4)}
        sage: S = Set([(G.P_symbol(), G.Q_symbol()) for G in l.values()])
        sage: S.cardinality()
        384

    Check the color-to-spin property for all permutations of size 4::

        sage: all(list(G.filling().values()).count(-1) == spin(G.P_symbol()) + spin(G.Q_symbol())
        ....:     for G in l.values())
        True

    Negating all signs transposes all the partitions::

        sage: W = SignedPermutations(4)
        sage: all(l[pi].P_symbol() == l[W([-e for e in pi])].P_symbol().conjugate()
        ....:     for pi in l)
        True

    Check part of Theorem 4.2.3 in [Lee1996]_::

        sage: def to_permutation(pi):
        ....:     pi1 = list(pi)
        ....:     n = len(pi1)
        ....:     pi2 = [-e for e in pi][::-1] + pi1
        ....:     return Permutation([e+n+1 if e<0 else e+n for e in pi2])
        sage: RuleRSK = GrowthDiagram.rules.RSK()
        sage: def good(pi):
        ....:     return Domino(pi).P_chain()[-1] == RuleRSK(to_permutation(pi)).P_chain()[-1]
        sage: all(good(pi) for pi in SignedPermutations(4))
        True

        sage: G = Domino(labels=[[1],[2,1]])
        Traceback (most recent call last):
        ...
        ValueError: [1] has smaller rank than [2, 1] but is not covered by it in Q

        sage: G = Domino(labels=[[2,1],[1]])
        Traceback (most recent call last):
        ...
        ValueError: [1] has smaller rank than [2, 1] but is not covered by it in P
    """
    r = 2
    zero = _make_partition([])

    def normalize_vertex(self, v):
        """
        Return ``v`` as a partition.

        EXAMPLES::

            sage: Domino = GrowthDiagram.rules.Domino()
            sage: Domino.normalize_vertex([3,1]).parent()
            Partitions
        """
        return _make_partition(v)

    def vertices(self, n):
        r"""
        Return the vertices of the dual graded graph on level ``n``.

        EXAMPLES::

            sage: Domino = GrowthDiagram.rules.Domino()
            sage: Domino.vertices(2)
            [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]
        """
        return [la for la in Partitions(2*n) if len(la.core(2)) == 0]

    def rank(self, v):
        r"""
        Return the rank of ``v``.

        The rank of a vertex is half the size of the partition,
        which equals the number of dominoes in any filling.

        EXAMPLES::

            sage: Domino = GrowthDiagram.rules.Domino()
            sage: Domino.rank(Domino.vertices(3)[0])
            3
        """
        return v.size() // 2

    def is_P_edge(self, v, w):
        r"""
        Return whether ``(v, w)`` is a `P`-edge of ``self``.

        ``(v, w)`` is an edge if ``v`` is obtained from ``w`` by deleting
        a domino.

        EXAMPLES::

            sage: Domino = GrowthDiagram.rules.Domino()
            sage: v = Domino.vertices(2)[1]; ascii_art(v)
            ***
            *
            sage: ascii_art([w for w in Domino.vertices(3) if Domino.is_P_edge(v, w)])
            [             *** ]
            [             *   ]
            [ *****  ***  *   ]
            [ *    , ***, *   ]
            sage: [w for w in Domino.vertices(4) if Domino.is_P_edge(v, w)]
            []
        """
        try:
            (row_1, col_1), (row_2, col_2) = SkewPartition([w, v]).cells()
        except ValueError:
            return False
        return row_1 == row_2 or col_1 == col_2

    is_Q_edge = is_P_edge

    def P_symbol(self, P_chain):
        r"""
        Return the labels along the vertical boundary of a rectangular
        growth diagram as a (skew) domino tableau.

        EXAMPLES::

            sage: Domino = GrowthDiagram.rules.Domino()
            sage: G = Domino([[0,1,0],[0,0,-1],[1,0,0]])
            sage: G.P_symbol().pp()
            1  1
            2  3
            2  3
        """
        return SkewTableau(chain=P_chain)

    Q_symbol = P_symbol

    def forward_rule(self, y, t, x, content):
        r"""
        Return the output shape given three shapes and the content.

        See [Lam2004]_ Section 3.1.

        INPUT:

        - ``y, t, x`` -- three partitions from a cell in a
          growth diagram, labelled as::

              t x
              y

        - ``content`` -- `-1`, `0` or `1`; the content of the cell

        OUTPUT:

        The fourth partition according to domino insertion.

        EXAMPLES::

            sage: Domino = GrowthDiagram.rules.Domino()

        Rule 1::

            sage: Domino.forward_rule([], [], [], 1)
            [2]

            sage: Domino.forward_rule([1,1], [1,1], [1,1], 1)
            [3, 1]

        Rule 2::

            sage: Domino.forward_rule([1,1], [1,1], [1,1], -1)
            [1, 1, 1, 1]

        Rule 3::

            sage: Domino.forward_rule([1,1], [1,1], [2,2], 0)
            [2, 2]

        Rule 4::

            sage: Domino.forward_rule([2,2,2], [2,2], [3,3], 0)
            [3, 3, 2]

            sage: Domino.forward_rule([2], [], [1,1], 0)
            [2, 2]

            sage: Domino.forward_rule([1,1], [], [2], 0)
            [2, 2]

            sage: Domino.forward_rule([2], [], [2], 0)
            [2, 2]

            sage: Domino.forward_rule([4], [2], [4], 0)
            [4, 2]

            sage: Domino.forward_rule([1,1,1,1], [1,1], [1,1,1,1], 0)
            [2, 2, 1, 1]

            sage: Domino.forward_rule([2,1,1], [2], [4], 0)
            [4, 1, 1]
        """
        def union(la, mu):
            r"""
            Return the union of the two partitions.
            """
            return [max(p,q) for (p,q) in zip_longest(la, mu, fillvalue=0)]

        if content not in [0,1,-1]:
            raise ValueError("domino: the content of the filling must be in {-1,0,1}")

        if content == 1:
            if not (x == t == y):
                raise ValueError("all shapes must be equal")
            if t == []:
                z = [2]
            else:
                z = [t[0] + 2] + t[1:]

        elif content == -1:
            if not (x == t == y):
                raise ValueError("all shapes must be equal")
            z = t + [1,1]

        elif content == 0 and (t == x or t == y):
            z = union(x, y)

        else:
            # content == 0 and t differs from x and y by
            # domino's gamma1 and gamma3

            # the following is certainly very slow
            gamma3 = set(SkewPartition([y, t]).cells())
            gamma1 = set(SkewPartition([x, t]).cells())
            diff = gamma1.intersection(gamma3)
            cell1, cell2 = gamma3
            if len(diff) == 0:
                z = union(x, y)

            elif len(diff) == 1:
                z = copy(x)
                # diff is a single cell
                (k,l) = diff.pop()
                # add (k+1, l+1) to x
                # either (k, l+1) or (k+1, l) must also be added
                if z[k] <= l + 1:
                    z[k] += 1
                    z[k+1] += 1
                else:
                    if len(z) <= k + 1:
                        z += [2]
                    else:
                        z[k+1] += 2

            # diff has size 2, that is x == y
            elif cell1[0] == cell2[0]:
                z = copy(x)
                # a horizontal domino - add 2 to row below of gamma
                if len(z) <= cell1[0] + 1:
                    z += [2]
                else:
                    z[cell1[0]+1] += 2

            else:
                z = copy(x)
                # a vertical domino - add 2 to column right of gamma
                # find first row shorter than cell1[1]+1
                for r, p in enumerate(z):
                    if p <= cell1[1] + 1:
                        z[r] += 1
                        z[r+1] += 1
                        break
                else:
                    raise NotImplementedError("domino: cannot call forward rule with shapes %s and content %s"
                                              % ((y, t, x), content))

        return z

#####################################################################
## Set the rules available from GrowthDiagram.rules.<tab>
#####################################################################

class Rules(object):
    """
    Catalog of rules for growth diagrams.
    """
    ShiftedShapes = RuleShiftedShapes
    LLMS = RuleLLMS
    BinaryWord = RuleBinaryWord
    Sylvester = RuleSylvester
    YoungFibonacci = RuleYoungFibonacci
    RSK = RuleRSK
    Burge = RuleBurge
    Domino = RuleDomino

GrowthDiagram.rules = Rules
