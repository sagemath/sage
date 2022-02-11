r"""
An introduction to crystals
===========================

Informally, a crystal `\mathcal{B}` is an oriented graph with edges
colored in some set `I` such that, for each `i\in I`, each node `x`
has:

- at most one `i`-successor, denoted `f_i x`;

- at most one `i`-predecessor, denoted `e_i x`.

By convention, one writes `f_i x=\emptyset` and `e_i x=\emptyset` when
`x` has no successor resp. predecessor.

One may think of `\mathcal{B}` as essentially a deterministic
automaton whose dual is also deterministic; in this context, the
`f_i`'s and `e_i`'s are respectively the transition functions of the
automaton and of its dual, and `\emptyset` is the sink.

A crystal comes further endowed with a weight function
`\operatorname{wt} : \mathcal{B} \to L` which satisfies
appropriate conditions.

In combinatorial representation theory, crystals are used as
combinatorial data to model representations of Lie algebra.

Axiomatic definition
--------------------

Let `C` be a Cartan type (:class:`CartanType`) with index set `I`,
and `L` be a realization of the weight lattice of the type `C`.
Let `\alpha_i` and `\alpha^{\vee}_i` denote the simple roots and
coroots respectively.

A type `C` crystal is a non-empty set `\mathcal{B}` endowed with maps
`\operatorname{wt} : \mathcal{B} \to L`,
`e_i, f_i : \mathcal{B} \to \mathcal{B} \cup \{\emptyset\}`, and
`\varepsilon_i, \varphi_i : \mathcal{B} \to \ZZ \cup \{-\infty\}`
for `i \in I` satisfying the following properties for all `i \in I`:

- for `b, b^{\prime} \in \mathcal{B}`, we have
  `f_i b^{\prime} = b` if and only if `e_i b = b^{\prime}`;

- if `e_i b \in \mathcal{B}`, then:

  * `\operatorname{wt}(e_i b) = \operatorname{wt}(b) + \alpha_i`,
  * `\varepsilon_i(e_i b) = \varepsilon_i(b) - 1`,
  * `\varphi_i(e_i b) = \varphi_i(b) + 1`;

- if `f_i b \in \mathcal{B}`, then:

  * `\operatorname{wt}(f_i b) = \operatorname{wt}(b) - \alpha_i`,
  * `\varepsilon_i(f_i b) = \varepsilon_i(b) + 1`,
  * `\varphi_i(f_i b) = \varphi_i(b) - 1`;

- `\varphi_i(b) = \varepsilon_i(b) + \langle \alpha^{\vee}_i,
  \operatorname{wt}(b) \rangle`,

- if `\varphi_i(b) = -\infty` for `b \in \mathcal{B}`,
  then `e_i b = f_i b = \emptyset`.

Some further conditions are required to guarantee that this data
indeed models a representation of a Lie algebra. For finite simply
laced types a complete characterization is given by Stembridge's local
axioms [Ste2003]_.

EXAMPLES:

We construct the type `A_5` crystal on letters (or in representation
theoretic terms, the highest weight crystal of type `A_5`
corresponding to the highest weight `\Lambda_1`)::

    sage: C = crystals.Letters(['A',5]); C
    The crystal of letters for type ['A', 5]

It has a single highest weight element::

    sage: C.highest_weight_vectors()
    (1,)

A crystal is an enumerated set (see :class:`EnumeratedSets`); and we
can count and list its elements in the usual way::

    sage: C.cardinality()
    6
    sage: C.list()
    [1, 2, 3, 4, 5, 6]

as well as use it in for loops::

    sage: [x for x in C]
    [1, 2, 3, 4, 5, 6]

Here are some more elaborate crystals (see their respective
documentations)::

    sage: Tens = crystals.TensorProduct(C, C)
    sage: Spin = crystals.Spins(['B', 3])
    sage: Tab  = crystals.Tableaux(['A', 3], shape = [2,1,1])
    sage: Fast = crystals.FastRankTwo(['B', 2], shape = [3/2, 1/2])
    sage: KR = crystals.KirillovReshetikhin(['A',2,1],1,1)

One can get (currently) crude plotting via::

    sage: Tab.plot()
    Graphics object consisting of 52 graphics primitives

If dot2tex is installed, one can obtain nice latex pictures via::

    sage: K = crystals.KirillovReshetikhin(['A',3,1], 1,1)
    sage: view(K, pdflatex=True) # optional - dot2tex graphviz, not tested (opens external window)

or with colored edges::

    sage: K = crystals.KirillovReshetikhin(['A',3,1], 1,1)
    sage: G = K.digraph()
    sage: G.set_latex_options(color_by_label={0:"black", 1:"red", 2:"blue", 3:"green"})
    sage: view(G, pdflatex=True) # optional - dot2tex graphviz, not tested (opens external window)

For rank two crystals, there is an alternative method of getting
metapost pictures. For more information see ``C.metapost?``.

.. SEEALSO:: :ref:`The overview of crystal features in Sage<sage.combinat.crystals.all>`

.. TODO::

    -  Vocabulary and conventions:

       -  For a classical crystal: connected / highest weight /
          irreducible

       -  ...

    -  Layout instructions for plot() for rank 2 types

    -  RestrictionOfCrystal


The crystals library in Sage grew up from an initial implementation in
MuPAD-Combinat (see <MuPAD-Combinat>/lib/COMBINAT/crystals.mu).
"""

#*****************************************************************************
#       Copyright (C) 2007 Anne Schilling <anne at math.ucdavis.edu>
#                          Nicolas Thiery <nthiery at users.sf.net>
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
#****************************************************************************
# Acknowledgment: most of the design and implementation of this
# library is heavily inspired from MuPAD-Combinat.
#****************************************************************************

#from sage.structure.unique_representation import UniqueRepresentation
#from sage.structure.parent import Parent
#from sage.structure.element import Element
#from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
#from sage.graphs.all import DiGraph
#from sage.combinat import ranker
#from sage.combinat.root_system.weyl_characters import WeylCharacter
from sage.combinat.backtrack import GenericBacktracker

class CrystalBacktracker(GenericBacktracker):
    def __init__(self, crystal, index_set=None):
        r"""
        Time complexity: `O(nF)` amortized for each produced
        element, where `n` is the size of the index set, and `F` is
        the cost of computing `e` and `f` operators.

        Memory complexity: `O(D)` where `D` is the depth of the crystal.

        Principle of the algorithm:

        Let `C` be a classical crystal. It's an acyclic graph where each
        connected component has a unique element without predecessors (the
        highest weight element for this component). Let's assume for
        simplicity that `C` is irreducible (i.e. connected) with highest
        weight element `u`.

        One can define a natural spanning tree of `C` by taking
        `u` as the root of the tree, and for any other element
        `y` taking as ancestor the element `x` such that
        there is an `i`-arrow from `x` to `y` with
        `i` minimal. Then, a path from `u` to `y`
        describes the lexicographically smallest sequence
        `i_1,\dots,i_k` such that
        `(f_{i_k} \circ f_{i_1})(u)=y`.

        Morally, the iterator implemented below just does a depth first
        search walk through this spanning tree. In practice, this can be
        achieved recursively as follows: take an element `x`, and
        consider in turn each successor `y = f_i(x)`, ignoring
        those such that `y = f_j(x^{\prime})` for some `x^{\prime}` and
        `j<i` (this can be tested by computing `e_j(y)`
        for `j<i`).

        EXAMPLES::

            sage: from sage.combinat.crystals.crystals import CrystalBacktracker
            sage: C = crystals.Tableaux(['B',3],shape=[3,2,1])
            sage: CB = CrystalBacktracker(C)
            sage: len(list(CB))
            1617
            sage: CB = CrystalBacktracker(C, [1,2])
            sage: len(list(CB))
            8
        """
        GenericBacktracker.__init__(self, None, None)
        self._crystal = crystal
        if index_set is None:
            self._index_set = crystal.index_set()
        else:
            self._index_set = index_set

    def _rec(self, x, state):
        """
        Return an iterator for the (immediate) children of ``x`` in the search
        tree.

        EXAMPLES::

            sage: from sage.combinat.crystals.crystals import CrystalBacktracker
            sage: C = crystals.Letters(['A', 5])
            sage: CB = CrystalBacktracker(C)
            sage: list(CB._rec(C(1), 'n/a'))
            [(2, 'n/a', True)]
        """
        #We will signal the initial case by having a object and state
        #of None and consider it separately.
        if x is None and state is None:
            for gen in self._crystal.highest_weight_vectors():
                yield gen, "n/a", True
            return

        # Run through the children y of x
        for i in self._index_set:
            y = x.f(i)
            if y is None:
                continue
            # Ignore those which can be reached by an arrow with smaller label
            hasParent = False
            for j in self._index_set:
                if j == i:
                    break
                if not y.e(j) is None:
                    hasParent = True
                    break
            if hasParent:
                continue

            # yield y and all elements further below
            yield y, "n/a", True

