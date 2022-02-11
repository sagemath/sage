# -*- coding: utf-8 -*-
r"""
Tamari Interval-posets

This module implements Tamari interval-posets: combinatorial objects which
represent intervals of the Tamari order. They have been introduced in
[CP2015]_ and allow for many combinatorial operations on Tamari intervals.
In particular, they are linked to :class:`DyckWords` and :class:`BinaryTrees`.
An introduction into Tamari interval-posets is given in Chapter 7
of [Pons2013]_.

The Tamari lattice can be defined as a lattice structure on either of several
classes of Catalan objects, especially binary trees and Dyck paths
[Tam1962]_ [HT1972]_ [Sta-EC2]_. An interval can be seen as
a pair of comparable elements. The number of intervals has been given in
[Cha2008]_.

AUTHORS:

- Viviane Pons 2014: initial implementation
- Frédéric Chapoton 2014: review
- Darij Grinberg 2014: review
- Travis Scrimshaw 2014: review
"""
# ****************************************************************************
#       Copyright (C) 2013 Viviane Pons <viviane.pons@univie.ac.at>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from __future__ import annotations
from collections.abc import Iterator

from sage.categories.enumerated_sets import EnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.posets import Posets
from sage.categories.all import Monoids
from sage.combinat.posets.posets import Poset, FinitePoset
from sage.categories.finite_posets import FinitePosets
from sage.combinat.binary_tree import BinaryTrees
from sage.combinat.binary_tree import LabelledBinaryTrees, LabelledBinaryTree
from sage.combinat.dyck_word import DyckWords
from sage.combinat.permutation import Permutation
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.cachefunc import cached_method
from sage.misc.latex import latex
from sage.misc.lazy_attribute import lazy_attribute
from sage.rings.integer import Integer
from sage.rings.all import NN
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.family import Family
from sage.structure.element import Element
from sage.structure.global_options import GlobalOptions
from sage.structure.parent import Parent
from sage.structure.richcmp import richcmp, op_NE, op_EQ
from sage.structure.unique_representation import UniqueRepresentation
from sage.graphs.digraph import DiGraph


class TamariIntervalPoset(Element,
        metaclass=InheritComparisonClasscallMetaclass):
    r"""
    The class of Tamari interval-posets.

    An interval-poset is a labelled poset of size `n`, with labels
    `1, 2, \ldots, n`, satisfying the following conditions:

    - if `a < c` (as integers) and `a` precedes `c` in the poset, then,
      for all `b` such that `a < b < c`, `b` precedes `c`,

    - if `a < c` (as integers) and `c` precedes `a` in the poset, then,
      for all `b` such that `a < b < c`, `b` precedes `a`.

    We use the word "precedes" here to distinguish the poset order and
    the natural order on numbers. "Precedes" means "is smaller than
    with respect to the poset structure"; this does not imply a
    covering relation.

    Interval-posets of size `n` are in bijection with intervals of
    the Tamari lattice of binary trees of size `n`. Specifically, if
    `P` is an interval-poset of size `n`, then the set of linear
    extensions of `P` (as permutations in `S_n`) is an interval in the
    right weak order (see
    :meth:`~sage.combinat.permutation.Permutation.permutohedron_lequal`),
    and is in fact the preimage of an interval in the Tamari lattice (of
    binary trees of size `n`) under the operation which sends a
    permutation to its right-to-left binary search tree
    (:meth:`~sage.combinat.permutation.Permutation.binary_search_tree`
    with the ``left_to_right`` variable set to ``False``)
    without its labelling.

    INPUT:

    - ``size`` -- an integer, the size of the interval-posets (number of
      vertices)

    - ``relations`` -- a list (or tuple) of pairs ``(a,b)`` (themselves
      lists or tuples), each representing a relation of the form
      '`a` precedes `b`' in the poset.

    - ``check`` -- (default: ``True``) whether to check the interval-poset
      condition or not.

    .. WARNING::

        The ``relations`` input can be a list or tuple, but not an
        iterator (nor should its entries be iterators).

    NOTATION:

    Here and in the following, the signs `<` and `>` always refer to
    the natural ordering on integers, whereas the word "precedes" refers
    to the order of the interval-poset. "Minimal" and "maximal" refer
    to the natural ordering on integers.

    The *increasing relations* of an interval-poset `P` mean the pairs
    `(a, b)` of elements of `P` such that `a < b` as integers and `a`
    precedes `b` in `P`. The *initial forest* of `P` is the poset
    obtained by imposing (only) the increasing relations on the ground
    set of `P`. It is a sub-interval poset of `P`, and is a forest with
    its roots on top. This forest is usually given the structure of a
    planar forest by ordering brother nodes by their labels; it then has
    the property that if its nodes are traversed in post-order
    (see :meth:`~sage.combinat.abstract_tree.AbstractTree.post_order_traversal`,
    and traverse the trees of the forest from left to right as well),
    then the labels encountered are `1, 2, \ldots, n` in this order.

    The *decreasing relations* of an interval-poset `P` mean the pairs
    `(a, b)` of elements of `P` such that `b < a` as integers and `a`
    precedes `b` in `P`. The *final forest* of `P` is the poset
    obtained by imposing (only) the decreasing relations on the ground
    set of `P`. It is a sub-interval poset of `P`, and is a forest with
    its roots on top. This forest is usually given the structure of a
    planar forest by ordering brother nodes by their labels; it then has
    the property that if its nodes are traversed in pre-order
    (see :meth:`~sage.combinat.abstract_tree.AbstractTree.pre_order_traversal`,
    and traverse the trees of the forest from left to right as well),
    then the labels encountered are `1, 2, \ldots, n` in this order.

    EXAMPLES::

        sage: TamariIntervalPoset(0,[])
        The Tamari interval of size 0 induced by relations []
        sage: TamariIntervalPoset(3,[])
        The Tamari interval of size 3 induced by relations []
        sage: TamariIntervalPoset(3,[(1,2)])
        The Tamari interval of size 3 induced by relations [(1, 2)]
        sage: TamariIntervalPoset(3,[(1,2),(2,3)])
        The Tamari interval of size 3 induced by relations [(1, 2), (2, 3)]
        sage: TamariIntervalPoset(3,[(1,2),(2,3),(1,3)])
        The Tamari interval of size 3 induced by relations [(1, 2), (2, 3)]
        sage: TamariIntervalPoset(3,[(1,2),(3,2)])
        The Tamari interval of size 3 induced by relations [(1, 2), (3, 2)]
        sage: TamariIntervalPoset(3,[[1,2],[2,3]])
        The Tamari interval of size 3 induced by relations [(1, 2), (2, 3)]
        sage: TamariIntervalPoset(3,[[1,2],[2,3],[1,2],[1,3]])
        The Tamari interval of size 3 induced by relations [(1, 2), (2, 3)]

        sage: TamariIntervalPoset(3,[(3,4)])
        Traceback (most recent call last):
        ...
        ValueError: the relations do not correspond to the size of the poset

        sage: TamariIntervalPoset(2,[(2,1),(1,2)])
        Traceback (most recent call last):
        ...
        ValueError: The graph is not directed acyclic

        sage: TamariIntervalPoset(3,[(1,3)])
        Traceback (most recent call last):
        ...
        ValueError: this does not satisfy the Tamari interval-poset condition

    It is also possible to transform a poset directly into an interval-poset::

        sage: TIP = TamariIntervalPosets()
        sage: p = Poset(([1,2,3], [(1,2)]))
        sage: TIP(p)
        The Tamari interval of size 3 induced by relations [(1, 2)]
        sage: TIP(Poset({1: []}))
        The Tamari interval of size 1 induced by relations []
        sage: TIP(Poset({}))
        The Tamari interval of size 0 induced by relations []
    """
    @staticmethod
    def __classcall_private__(cls, *args, **opts):
        r"""
        Ensure that interval-posets created by the enumerated sets and
        directly are the same and that they are instances of
        :class:`TamariIntervalPoset`.

        TESTS::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: ip.parent()
            Interval-posets
            sage: type(ip)
            <class 'sage.combinat.interval_posets.TamariIntervalPosets_all_with_category.element_class'>

            sage: ip2 = TamariIntervalPosets()(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: ip2.parent() is ip.parent()
            True
            sage: type(ip) is type(ip2)
            True

            sage: ip3 = TamariIntervalPosets(4)([(2,4),(3,4),(2,1),(3,1)])
            sage: ip3.parent() is ip.parent()
            False
            sage: type(ip3) is type(ip)
            True
        """
        P = TamariIntervalPosets_all()
        return P.element_class(P, *args, **opts)

    def __init__(self, parent, size, relations=None, check=True):
        r"""
        TESTS::

            sage: TamariIntervalPoset(3,[(1,2),(3,2)]).parent()
            Interval-posets
            sage: P = Poset(DiGraph([(4,1),(3,1),(2,1)]))
            sage: TamariIntervalPoset(P).parent()
            Interval-posets
        """
        if relations is None:
            relations = []
        if isinstance(size, FinitePoset):
            # first argument is a poset
            self._poset = size
            self._size = size.cardinality()
        else:
            # arguments are size and relations
            self._size = size
            self._poset = Poset((list(range(1, size + 1)), relations))
            if self._poset.cardinality() != size:
                # This can happen as the Poset constructor automatically adds
                # in elements from the relations.
                raise ValueError("the relations do not correspond to the size of the poset")

        if check and not TamariIntervalPosets.check_poset(self._poset):
            raise ValueError("this does not satisfy the Tamari interval-poset condition")

        Element.__init__(self, parent)

        self._cover_relations = tuple(self._poset.cover_relations())
        self._latex_options = dict()

    def set_latex_options(self, D):
        r"""
        Set the latex options for use in the ``_latex_`` function.

        The default values are set in the ``__init__`` function.

        - ``tikz_scale`` -- (default: 1) scale for use with the tikz package

        - ``line_width`` -- (default: 1 * ``tikz_scale``) value representing the
          line width

        - ``color_decreasing`` -- (default: red) the color for decreasing
          relations

        - ``color_increasing`` -- (default: blue) the color for increasing
          relations

        - ``hspace`` -- (default: 1) the difference between horizontal
          coordinates of adjacent vertices

        - ``vspace`` -- (default: 1) the difference between vertical
          coordinates of adjacent vertices

        INPUT:

        - ``D`` -- a dictionary with a list of latex parameters to change

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: ip.latex_options()["color_decreasing"]
            'red'
            sage: ip.set_latex_options({"color_decreasing":'green'})
            sage: ip.latex_options()["color_decreasing"]
            'green'
            sage: ip.set_latex_options({"color_increasing":'black'})
            sage: ip.latex_options()["color_increasing"]
            'black'

        To change the default options for all interval-posets, use the
        parent's latex options::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: ip2 = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: ip.latex_options()["color_decreasing"]
            'red'
            sage: ip2.latex_options()["color_decreasing"]
            'red'
            sage: TamariIntervalPosets.options(latex_color_decreasing='green')
            sage: ip.latex_options()["color_decreasing"]
            'green'
            sage: ip2.latex_options()["color_decreasing"]
            'green'

        Next we set a local latex option and show the global option does not
        override it::

            sage: ip.set_latex_options({"color_decreasing": 'black'})
            sage: ip.latex_options()["color_decreasing"]
            'black'
            sage: TamariIntervalPosets.options(latex_color_decreasing='blue')
            sage: ip.latex_options()["color_decreasing"]
            'black'
            sage: ip2.latex_options()["color_decreasing"]
            'blue'
            sage: TamariIntervalPosets.options._reset()
        """
        for opt in D:
            self._latex_options[opt] = D[opt]

    def latex_options(self) -> dict:
        r"""
        Return the latex options for use in the ``_latex_`` function as a
        dictionary.

        The default values are set using the options.

        - ``tikz_scale`` -- (default: 1) scale for use with the tikz package

        - ``line_width`` -- (default: 1) value representing the line width
          (additionally scaled by ``tikz_scale``)

        - ``color_decreasing`` -- (default: ``'red'``) the color for
          decreasing relations

        - ``color_increasing`` -- (default: ``'blue'``) the color for
          increasing relations

        - ``hspace`` -- (default: 1) the difference between horizontal
          coordinates of adjacent vertices

        - ``vspace`` -- (default: 1) the difference between vertical
          coordinates of adjacent vertices

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: ip.latex_options()['color_decreasing']
            'red'
            sage: ip.latex_options()['hspace']
            1
        """
        d = self._latex_options.copy()
        if "tikz_scale" not in d:
            d["tikz_scale"] = self.parent().options["latex_tikz_scale"]
        if "line_width" not in d:
            d["line_width"] = self.parent().options["latex_line_width_scalar"] * d["tikz_scale"]
        if "color_decreasing" not in d:
            d["color_decreasing"] = self.parent().options["latex_color_decreasing"]
        if "color_increasing" not in d:
            d["color_increasing"] = self.parent().options["latex_color_increasing"]
        if "hspace" not in d:
            d["hspace"] = self.parent().options["latex_hspace"]
        if "vspace" not in d:
            d["vspace"] = self.parent().options["latex_vspace"]
        return d

    def _find_node_positions(self, hspace=1, vspace=1) -> dict[int, list]:
        """
        Compute a nice embedding.

        If `x` precedes `y`, then `y` will always be placed on top of `x`
        and/or to the right of `x`.
        Decreasing relations tend to be drawn vertically and increasing
        relations horizontally.
        The algorithm tries to avoid superposition but on big
        interval-posets, it might happen.

        OUTPUT:

        a dictionary {vertex: (x,y)}

        EXAMPLES::

            sage: ti = TamariIntervalPosets(4)[2]
            sage: list(ti._find_node_positions().values())
            [[0, 0], [0, -1], [0, -2], [1, -2]]
        """
        node_positions = {}

        to_draw = [(1, 0)]
        current_parent = [self.increasing_parent(1)]
        parenty = [0]
        x = 0
        y = 0
        for i in range(2, self.size() + 1):
            decreasing_parent = self.decreasing_parent(i)
            increasing_parent = self.increasing_parent(i)
            while to_draw and (decreasing_parent is None or
                               decreasing_parent < to_draw[-1][0]):
                n = to_draw.pop()
                node_positions[n[0]] = [x, n[1]]
            if i != current_parent[-1]:
                if (not self.le(i, i - 1) and decreasing_parent is not None):
                    x += hspace
                    if current_parent[-1] is not None:
                        y -= vspace
                else:
                    y -= vspace
                if increasing_parent != current_parent[-1]:
                    current_parent.append(increasing_parent)
                    parenty.append(y)
                nodey = y
            else:
                current_parent.pop()
                x += hspace
                nodey = parenty.pop()
                if not current_parent or increasing_parent != current_parent[-1]:
                    current_parent.append(increasing_parent)
                    parenty.append(nodey)
            to_draw.append((i, nodey))

        for n in to_draw:
            node_positions[n[0]] = [x, n[1]]
        return node_positions

    def plot(self, **kwds):
        """
        Return a picture.

        The picture represents the Hasse diagram, where the covers are
        colored in blue if they are increasing and in red if they are
        decreasing.

        This uses the same coordinates as the latex view.

        EXAMPLES::

            sage: ti = TamariIntervalPosets(4)[2]
            sage: ti.plot()
            Graphics object consisting of 6 graphics primitives

        TESTS::

            sage: ti = TamariIntervalPoset(3, [[2,1], [2,3]])
            sage: ti.plot()
            Graphics object consisting of 6 graphics primitives
        """
        c0 = 'blue'   # self.latex_options()["color_increasing"]
        c1 = 'red'    # self.latex_options()["color_decreasing"]
        G = self.poset().hasse_diagram()
        G.set_pos(self._find_node_positions())
        for a, b, c in G.edges():
            if a < b:
                G.set_edge_label(a, b, 0)
            else:
                G.set_edge_label(a, b, 1)
        return G.plot(color_by_label={0: c0, 1: c1}, **kwds)

    def _latex_(self) -> str:
        r"""
        A latex representation of ``self`` using the tikzpicture package.

        This picture shows the union of the Hasse diagrams of the
        initial and final forests.

        If `x` precedes `y`, then `y` will always be placed on top of `x`
        and/or to the right of `x`.
        Decreasing relations tend to be drawn vertically and increasing
        relations horizontally.
        The algorithm tries to avoid superposition but on big
        interval-posets, it might happen.

        You can use ``self.set_latex_options()`` to change default latex
        options. Or you can use the parent's options.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: latex(ip)
            \begin{tikzpicture}[scale=1]
            \node(T1) at (1,0) {1};
            \node(T2) at (0,-1) {2};
            \node(T3) at (1,-2) {3};
            \node(T4) at (2,-1) {4};
            \draw[line width = 0.5, color=red] (T3) -- (T1);
            \draw[line width = 0.5, color=red] (T2) -- (T1);
            \draw[line width = 0.5, color=blue] (T2) -- (T4);
            \draw[line width = 0.5, color=blue] (T3) -- (T4);
            \end{tikzpicture}

        TESTS::

            sage: ip = TamariIntervalPoset(0,[])
            sage: latex(ip)
            \begin{tikzpicture}[scale=1]
            \node(T0) at (0,0){$\emptyset$};\end{tikzpicture}
        """
        latex.add_package_to_preamble_if_available("tikz")
        latex_options = self.latex_options()
        start = "\\begin{tikzpicture}[scale=" + str(latex_options['tikz_scale']) + "]\n"
        end = "\\end{tikzpicture}"
        vspace = latex_options["vspace"]
        hspace = latex_options["hspace"]

        def draw_node(j, x, y) -> str:
            r"""
            Internal method to draw vertices
            """
            return "\\node(T" + str(j) + ") at (" + str(x) + "," + str(y) + ") {" + str(j) + "};\n"

        def draw_increasing(i, j) -> str:
            r"""
            Internal method to draw increasing relations
            """
            return "\\draw[line width = " + str(latex_options["line_width"]) + ", color=" + latex_options["color_increasing"] + "] (T" + str(i) + ") -- (T" + str(j) + ");\n"

        def draw_decreasing(i, j) -> str:
            r"""
            Internal method to draw decreasing relations
            """
            return "\\draw[line width = " + str(latex_options["line_width"]) + ", color=" + latex_options["color_decreasing"] + "] (T" + str(i) + ") -- (T" + str(j) + ");\n"

        if self.size() == 0:
            nodes = "\\node(T0) at (0,0){$\\emptyset$};"
            relations = ""
        else:
            positions = self._find_node_positions(hspace, vspace)
            nodes = ""  # latex for node declarations
            relations = ""  # latex for drawing relations
            for i in range(1, self.size() + 1):
                nodes += draw_node(i, *positions[i])
            for i, j in self.decreasing_cover_relations():
                relations += draw_decreasing(i, j)
            for i, j in self.increasing_cover_relations():
                relations += draw_increasing(i, j)

        return start + nodes + relations + end

    def poset(self):
        r"""
        Return ``self`` as a labelled poset.

        An interval-poset is indeed constructed from a labelled poset which
        is stored internally. This method allows to access the poset and
        all the associated methods.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(1,2),(3,2),(2,4),(3,4)])
            sage: pos = ip.poset(); pos
            Finite poset containing 4 elements
            sage: pos.maximal_chains()
            [[3, 2, 4], [1, 2, 4]]
            sage: pos.maximal_elements()
            [4]
            sage: pos.is_lattice()
            False
        """
        return self._poset

    def _mul_(self, other: TIP) -> TIP:
        """
        Return the associative product of ``self`` and ``other``.

        This is defined by taking the disjoint union of the relations
        of ``self`` with the relations of ``other`` shifted by `n`,
        where `n` is the size of ``self``.

        EXAMPLES::

            sage: T1 = TamariIntervalPoset(1,[])
            sage: T2 = TamariIntervalPoset(2,[[1,2]])
            sage: T1*T1
            The Tamari interval of size 2 induced by relations []
            sage: T2*T1
            The Tamari interval of size 3 induced by relations [(1, 2)]
            sage: T1*T2
            The Tamari interval of size 3 induced by relations [(2, 3)]
            sage: T2*T2
            The Tamari interval of size 4 induced by relations [(1, 2), (3, 4)]

        TESTS::

            sage: U = TamariIntervalPoset(0,[])
            sage: U*T1 == T1
            True
            sage: T2*U == T2
            True
        """
        n = self._size
        m = other.size()
        relations = self._poset.cover_relations()
        relations.extend([(i + n, j + n)
                          for i, j in other._poset.cover_relations_iterator()])
        P = FinitePoset(DiGraph([list(range(1, n + m + 1)), relations],
                                format='vertices_and_edges'))  # type:ignore
        return TamariIntervalPoset(P, check=False)  # type:ignore

    def factor(self) -> list[TamariIntervalPoset]:
        """
        Return the unique decomposition as a list of connected components.

        EXAMPLES::

            sage: factor(TamariIntervalPoset(2,[]))  # indirect doctest
            [The Tamari interval of size 1 induced by relations [],
             The Tamari interval of size 1 induced by relations []]

        .. SEEALSO:: :meth:`is_connected`

        TESTS::

            sage: T = TamariIntervalPosets(20).random_element()
            sage: facs = factor(T)
            sage: all(U.is_connected() for U in facs)
            True
            sage: T == prod(facs)
            True
        """
        hasse = self.poset().hasse_diagram()
        cc = hasse.connected_components_subgraphs()
        resu = []
        for comp in sorted(cc, key=min):
            shift = 1 - min(comp)
            comp.relabel(lambda i: i + shift)
            resu.append(TamariIntervalPoset(len(comp),
                                            comp.edges(labels=False)))
        return resu

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: len(set(hash(u) for u in TamariIntervalPosets(4)))
            68
        """
        pair = (self.size(), tuple(tuple(e) for e in self._cover_relations))
        return hash(pair)

    @cached_method
    def increasing_cover_relations(self) -> list[tuple[int, int]]:
        r"""
        Return the cover relations of the initial forest of ``self``.

        This is the poset formed by keeping only the relations of the form
        `a` precedes `b` with `a < b`.

        The initial forest of ``self`` is a forest with its roots
        being on top. It is also called the increasing poset of ``self``.

        .. WARNING::

            This method computes the cover relations of the initial
            forest. This is not identical with the cover relations of
            ``self`` which happen to be increasing!

        .. SEEALSO::

            :meth:`initial_forest`

        EXAMPLES::

            sage: TamariIntervalPoset(4,[(1,2),(3,2),(2,4),(3,4)]).increasing_cover_relations()
            [(1, 2), (2, 4), (3, 4)]
            sage: TamariIntervalPoset(3,[(1,2),(1,3),(2,3)]).increasing_cover_relations()
            [(1, 2), (2, 3)]
        """
        relations = []
        size = self.size()
        for i in range(1, size):
            for j in range(i + 1, size + 1):
                if self.le(i, j):
                    relations.append((i, j))
                    break
        return relations

    def increasing_roots(self) -> list[int]:
        r"""
        Return the root vertices of the initial forest of ``self``.

        These are the vertices `a` of ``self`` such that there is no
        `b > a` with `a` precedes `b`.

        OUTPUT:

        The list of all roots of the initial forest of ``self``, in
        decreasing order.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(3,5),(4,5)]); ip
            The Tamari interval of size 6 induced by relations [(1, 2), (3, 5), (4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]
            sage: ip.increasing_roots()
            [6, 5, 2]
            sage: ip.initial_forest().increasing_roots()
            [6, 5, 2]

        TESTS::

            sage: TamariIntervalPoset(0,[]).increasing_roots()
            []
        """
        size = self.size()
        if size == 0:
            return []
        roots = [size]
        root = size
        for i in range(size - 1, 0, -1):
            if not self.le(i, root):
                roots.append(i)
                root = i
        return roots

    def increasing_children(self, v) -> list[int]:
        r"""
        Return the children of ``v`` in the initial forest of ``self``.

        INPUT:

        - ``v`` -- an integer representing a vertex of ``self``
          (between 1 and ``size``)

        OUTPUT:

        The list of all children of ``v`` in the initial forest of
        ``self``, in decreasing order.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(3,5),(4,5)]); ip
            The Tamari interval of size 6 induced by relations [(1, 2), (3, 5), (4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]
            sage: ip.increasing_children(2)
            [1]
            sage: ip.increasing_children(5)
            [4, 3]
            sage: ip.increasing_children(1)
            []
        """
        children = []
        root = None
        for i in range(v - 1, 0, -1):
            if not self.le(i, v):
                break
            if root is None or not self.le(i, root):
                children.append(i)
                root = i
        return children

    def increasing_parent(self, v) -> None | int:
        r"""
        Return the vertex parent of ``v`` in the initial forest of ``self``.

        This is the lowest (as integer!) vertex `b > v` such that `v`
        precedes `b`. If there is no such vertex (that is, `v` is an
        increasing root), then ``None`` is returned.

        INPUT:

        - ``v`` -- an integer representing a vertex of ``self``
          (between 1 and ``size``)

        EXAMPLES::

            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(3,5),(4,5)]); ip
            The Tamari interval of size 6 induced by relations [(1, 2), (3, 5), (4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]
            sage: ip.increasing_parent(1)
            2
            sage: ip.increasing_parent(3)
            5
            sage: ip.increasing_parent(4)
            5
            sage: ip.increasing_parent(5) is None
            True
        """
        parent = None
        for i in range(self.size(), v, -1):
            if self.le(v, i):
                parent = i
        return parent

    @cached_method
    def decreasing_cover_relations(self) -> list[tuple[int, int]]:
        r"""
        Return the cover relations of the final forest of ``self``.

        This is the poset formed by keeping only the relations of the form
        `a` precedes `b` with `a > b`.

        The final forest of ``self`` is a forest with its roots
        being on top. It is also called the decreasing poset of ``self``.

        .. WARNING::

            This method computes the cover relations of the final
            forest. This is not identical with the cover relations of
            ``self`` which happen to be decreasing!

        .. SEEALSO::

            :meth:`final_forest`

        EXAMPLES::

            sage: TamariIntervalPoset(4,[(2,1),(3,2),(3,4),(4,2)]).decreasing_cover_relations()
            [(4, 2), (3, 2), (2, 1)]
            sage: TamariIntervalPoset(4,[(2,1),(4,3),(2,3)]).decreasing_cover_relations()
            [(4, 3), (2, 1)]
            sage: TamariIntervalPoset(3,[(2,1),(3,1),(3,2)]).decreasing_cover_relations()
            [(3, 2), (2, 1)]
        """
        relations = []
        for i in range(self.size(), 1, -1):
            for j in range(i - 1, 0, -1):
                if self.le(i, j):
                    relations.append((i, j))
                    break
        return relations

    def decreasing_roots(self) -> list[int]:
        r"""
        Return the root vertices of the final forest of ``self``.

        These are the vertices `b` such that there is no `a < b` with `b`
        preceding `a`.

        OUTPUT:

        The list of all roots of the final forest of ``self``, in
        increasing order.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(3,5),(4,5)]); ip
            The Tamari interval of size 6 induced by relations [(1, 2), (3, 5), (4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]
            sage: ip.decreasing_roots()
            [1, 2]
            sage: ip.final_forest().decreasing_roots()
            [1, 2]
        """
        if self.size() == 0:
            return []
        roots = [1]
        root = 1
        for i in range(2, self.size() + 1):
            if not self.le(i, root):
                roots.append(i)
                root = i
        return roots

    def decreasing_children(self, v) -> list[int]:
        r"""
        Return the children of ``v`` in the final forest of ``self``.

        INPUT:

        - ``v`` -- an integer representing a vertex of ``self``
          (between 1 and ``size``)

        OUTPUT:

        The list of all children of ``v`` in the final forest of ``self``,
        in increasing order.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(3,5),(4,5)]); ip
            The Tamari interval of size 6 induced by relations [(1, 2), (3, 5), (4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]
            sage: ip.decreasing_children(2)
            [3, 5]
            sage: ip.decreasing_children(3)
            [4]
            sage: ip.decreasing_children(1)
            []
        """
        children = []
        root = None
        for i in range(v + 1, self.size() + 1):
            if not self.le(i, v):
                break
            if root is None or not self.le(i, root):
                children.append(i)
                root = i
        return children

    def decreasing_parent(self, v) -> None | int:
        r"""
        Return the vertex parent of ``v`` in the final forest of ``self``.

        This is the highest (as integer!) vertex `a < v` such that ``v``
        precedes ``a``. If there is no such vertex (that is, `v` is a
        decreasing root), then ``None`` is returned.

        INPUT:

        - ``v`` -- an integer representing a vertex of ``self`` (between
          1 and ``size``)

        EXAMPLES::

            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(3,5),(4,5)]); ip
            The Tamari interval of size 6 induced by relations [(1, 2), (3, 5), (4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]
            sage: ip.decreasing_parent(4)
            3
            sage: ip.decreasing_parent(3)
            2
            sage: ip.decreasing_parent(5)
            2
            sage: ip.decreasing_parent(2) is None
            True
        """
        parent = None
        for i in range(1, v):
            if self.le(v, i):
                parent = i
        return parent

    def le(self, e1, e2) -> bool:
        r"""
        Return whether ``e1`` precedes or equals ``e2`` in ``self``.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: ip.le(1,2)
            True
            sage: ip.le(1,3)
            True
            sage: ip.le(2,3)
            True
            sage: ip.le(3,4)
            False
            sage: ip.le(1,1)
            True
        """
        return self._poset.le(e1, e2)

    def lt(self, e1, e2) -> bool:
        r"""
        Return whether ``e1`` strictly precedes ``e2`` in ``self``.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: ip.lt(1,2)
            True
            sage: ip.lt(1,3)
            True
            sage: ip.lt(2,3)
            True
            sage: ip.lt(3,4)
            False
            sage: ip.lt(1,1)
            False
        """
        return self._poset.lt(e1, e2)

    def ge(self, e1, e2) -> bool:
        r"""
        Return whether ``e2`` precedes or equals ``e1`` in ``self``.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: ip.ge(2,1)
            True
            sage: ip.ge(3,1)
            True
            sage: ip.ge(3,2)
            True
            sage: ip.ge(4,3)
            False
            sage: ip.ge(1,1)
            True
        """
        return self._poset.ge(e1, e2)

    def gt(self, e1, e2) -> bool:
        r"""
        Return whether ``e2`` strictly precedes ``e1`` in ``self``.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: ip.gt(2,1)
            True
            sage: ip.gt(3,1)
            True
            sage: ip.gt(3,2)
            True
            sage: ip.gt(4,3)
            False
            sage: ip.gt(1,1)
            False
        """
        return self._poset.gt(e1, e2)

    def size(self) -> Integer:
        r"""
        Return the size (number of vertices) of the interval-poset.

        EXAMPLES::

            sage: TamariIntervalPoset(3,[(2,1),(3,1)]).size()
            3
        """
        return self._size

    @cached_method
    def cubical_coordinates(self) -> tuple[int, ...]:
        """
        Return the cubical coordinates of ``self``.

        This provides a fast and natural way to order
        the set of interval-posets of a given size.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: ip.cubical_coordinates()
            (-1, -2, 0)

        TESTS::

            sage: ip = TamariIntervalPoset(1,[])
            sage: ip.cubical_coordinates()
            ()
            sage: ip = TamariIntervalPoset(3,[])
            sage: ip.cubical_coordinates()
            (0, 0)
            sage: ip = TamariIntervalPosets(10).random_element()
            sage: len(ip.cubical_coordinates())
            9
            sage: sorted(ip.cubical_coordinates() for ip in TamariIntervalPosets(2))
            [(-1,), (0,), (1,)]

        REFERENCES:

        - [Com2019]_
        """
        tup = [0] * (self.size() - 1)
        for i, j in self._poset.relations_iterator(strict=True):
            if i < j:
                tup[j - 2] -= 1
            else:
                tup[j - 1] += 1
        return tuple(tup)

    def complement(self) -> TIP:
        r"""
        Return the complement of the interval-poset ``self``.

        If `P` is a Tamari interval-poset of size `n`, then the
        *complement* of `P` is defined as the interval-poset `Q` whose
        base set is `[n] = \{1, 2, \ldots, n\}` (just as for `P`), but
        whose order relation has `a` precede `b` if and only if
        `n + 1 - a` precedes `n + 1 - b` in `P`.

        In terms of the Tamari lattice, the *complement* is the symmetric
        of ``self``. It is formed from the left-right symmetrized of
        the binary trees of the interval (switching left and right
        subtrees, see
        :meth:`~sage.combinat.binary_tree.BinaryTree.left_right_symmetry`).
        In particular, initial intervals are sent to final intervals and
        vice-versa.

        EXAMPLES::

            sage: TamariIntervalPoset(3, [(2, 1), (3, 1)]).complement()
            The Tamari interval of size 3 induced by relations [(1, 3), (2, 3)]
            sage: TamariIntervalPoset(0, []).complement()
            The Tamari interval of size 0 induced by relations []
            sage: ip = TamariIntervalPoset(4, [(1, 2), (2, 4), (3, 4)])
            sage: ip.complement() == TamariIntervalPoset(4, [(2, 1), (3, 1), (4, 3)])
            True
            sage: ip.lower_binary_tree() == ip.complement().upper_binary_tree().left_right_symmetry()
            True
            sage: ip.upper_binary_tree() == ip.complement().lower_binary_tree().left_right_symmetry()
            True
            sage: ip.is_initial_interval()
            True
            sage: ip.complement().is_final_interval()
            True
        """
        N = self._size + 1
        new_covers = [[N - i, N - j]
                      for i, j in self._poset.cover_relations_iterator()]
        P = FinitePoset(DiGraph([list(range(1, N)), new_covers],
                                format='vertices_and_edges'))  # type:ignore
        return TamariIntervalPoset(P, check=False)  # type:ignore

    def left_branch_involution(self) -> TIP:
        """
        Return the image of ``self`` by the left-branch involution.

        OUTPUT: an interval-poset

        .. SEEALSO:: :meth:`rise_contact_involution`

        EXAMPLES::

            sage: tip = TamariIntervalPoset(8, [(1,2), (2,4), (3,4), (6,7), (3,2), (5,4), (6,4), (8,7)])
            sage: t = tip.left_branch_involution(); t
            The Tamari interval of size 8 induced by relations [(1, 6), (2, 6),
            (3, 5), (4, 5), (5, 6), (6, 8), (7, 8), (7, 6), (4, 3), (3, 1),
            (2, 1)]
            sage: t.left_branch_involution() == tip
            True

        REFERENCES:

        - [Pons2018]_
        """
        gt = self.grafting_tree().left_border_symmetry()
        return TamariIntervalPosets.from_grafting_tree(gt)

    def rise_contact_involution(self) -> TIP:
        """
        Return the image of ``self`` by the rise-contact involution.

        OUTPUT: an interval-poset

        This is defined by conjugating the complement involution
        by the left-branch involution.

        .. SEEALSO:: :meth:`left_branch_involution`, :meth:`complement`

        EXAMPLES::

            sage: tip = TamariIntervalPoset(8, [(1,2), (2,4), (3,4), (6,7), (3,2), (5,4), (6,4), (8,7)])
            sage: t = tip.rise_contact_involution(); t
            The Tamari interval of size 8 induced by relations [(2, 8), (3, 8),
            (4, 5), (5, 7), (6, 7), (7, 8), (8, 1), (7, 2), (6, 2), (5, 3),
            (4, 3), (3, 2), (2, 1)]
            sage: t.rise_contact_involution() == tip
            True
            sage: tip.lower_dyck_word().number_of_touch_points() == t.upper_dyck_word().number_of_initial_rises()
            True
            sage: tip.number_of_tamari_inversions() == t.number_of_tamari_inversions()
            True

        REFERENCES:

        - [Pons2018]_
        """
        t = self.left_branch_involution().complement()
        return t.left_branch_involution()

    def insertion(self, i) -> TIP:
        r"""
        Return the Tamari insertion of an integer `i` into the
        interval-poset ``self``.

        If `P` is a Tamari interval-poset of size `n` and `i` is an
        integer with `1 \leq i \leq n+1`, then the Tamari insertion of
        `i` into `P` is defined as the Tamari interval-poset of size
        `n+1` which corresponds to the interval `[C_1, C_2]` on the
        Tamari lattice, where the binary trees `C_1` and `C_2` are
        defined as follows: We write the interval-poset `P` as
        `[B_1, B_2]` for two binary trees `B_1` and `B_2`. We label
        the vertices of each of these two trees with the integers
        `1, 2, \ldots, i-1, i+1, i+2, \ldots, n+1` in such a way that
        the trees are binary search trees (this labelling is unique).
        Then, we insert `i` into each of these trees (in the way as
        explained in
        :meth:`~sage.combinat.binary_tree.LabelledBinaryTree.binary_search_insert`).
        The shapes of the resulting two trees are denoted `C_1` and
        `C_2`.

        An alternative way to construct the insertion of `i` into
        `P` is by relabeling each vertex `u` of `P` satisfying
        `u \geq i` (as integers) as `u+1`, and then adding a vertex
        `i` which should precede `i-1` and `i+1`.

        .. TODO::

            To study this, it would be more natural to define
            interval-posets on arbitrary ordered sets rather than just
            on `\{1, 2, \ldots, n\}`.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4, [(2, 3), (4, 3)]); ip
            The Tamari interval of size 4 induced by relations [(2, 3), (4, 3)]
            sage: ip.insertion(1)
            The Tamari interval of size 5 induced by relations [(1, 2), (3, 4), (5, 4)]
            sage: ip.insertion(2)
            The Tamari interval of size 5 induced by relations [(2, 3), (3, 4), (5, 4), (2, 1)]
            sage: ip.insertion(3)
            The Tamari interval of size 5 induced by relations [(2, 4), (3, 4), (5, 4), (3, 2)]
            sage: ip.insertion(4)
            The Tamari interval of size 5 induced by relations [(2, 3), (4, 5), (5, 3), (4, 3)]
            sage: ip.insertion(5)
            The Tamari interval of size 5 induced by relations [(2, 3), (5, 4), (4, 3)]

            sage: ip = TamariIntervalPoset(0, [])
            sage: ip.insertion(1)
            The Tamari interval of size 1 induced by relations []

            sage: ip = TamariIntervalPoset(1, [])
            sage: ip.insertion(1)
            The Tamari interval of size 2 induced by relations [(1, 2)]
            sage: ip.insertion(2)
            The Tamari interval of size 2 induced by relations [(2, 1)]

        TESTS:

        Verifying that the two ways of computing insertion are
        equivalent::

            sage: def insert_alternative(T, i):
            ....:     # Just another way to compute the insertion of i into T.
            ....:     from sage.combinat.binary_tree import LabelledBinaryTree
            ....:     B1 = T.lower_binary_tree().canonical_labelling()
            ....:     B2 = T.upper_binary_tree().canonical_labelling()
            ....:     C1 = B1.binary_search_insert(i)
            ....:     C2 = B2.binary_search_insert(i)
            ....:     return TamariIntervalPosets.from_binary_trees(C1, C2)

        We should have relabelled the trees to "make space" for a label i,
        but we did not, because it does not make a difference: The
        binary search insertion will go precisely the same, because
        an integer equal to the label of the root gets sent onto
        the left branch.

            sage: def test_equivalence(n):
            ....:     for T in TamariIntervalPosets(n):
            ....:         for i in range(1, n + 2):
            ....:             if not (insert_alternative(T, i) == T.insertion(i)):
            ....:                 print(T, i)
            ....:                 return False
            ....:     return True
            sage: test_equivalence(3)
            True

            sage: ti = TamariIntervalPosets(3).an_element()
            sage: ti.insertion(6)
            Traceback (most recent call last):
            ...
            ValueError: integer to be inserted not in the appropriate interval
        """
        n = self._size
        if not 0 < i <= n + 1:
            raise ValueError("integer to be inserted not "
                             "in the appropriate interval")

        def add1(u):
            if u >= i:
                return u + 1
            return u
        rels = [(add1(a), add1(b))
                for (a, b) in self.decreasing_cover_relations()]
        rels += [(add1(a), add1(b))
                 for (a, b) in self.increasing_cover_relations()]
        rels += [(k, k - 1) for k in [i] if i > 1]
        rels += [(k, k + 1) for k in [i] if i <= n]
        return TamariIntervalPoset(n + 1, rels)

    def _repr_(self) -> str:
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: TamariIntervalPoset(3,[(2,1),(3,1)])
            The Tamari interval of size 3 induced by relations [(3, 1), (2, 1)]
            sage: TamariIntervalPoset(3,[(3,1),(2,1)])
            The Tamari interval of size 3 induced by relations [(3, 1), (2, 1)]
            sage: TamariIntervalPoset(3,[(2,3),(2,1)])
            The Tamari interval of size 3 induced by relations [(2, 3), (2, 1)]
        """
        msg = "The Tamari interval of size {} induced by relations {}"
        return msg.format(self.size(),
                          self.increasing_cover_relations() +
                          self.decreasing_cover_relations())

    def _ascii_art_(self):
        """
        Return an ascii art picture of ``self``.

        This is a picture of the Hasse diagram. Vertices from `1` to `n`
        are placed on the diagonal from top-left to bottom-right.
        Then increasing covers are drawn above the diagonal
        and decreasing covers are drawn below the diagonal.

        EXAMPLES::

            sage: T = TamariIntervalPosets(5)[56]
            sage: ascii_art(T)
             O-----------+
                O--------+
                +--O--+  |
                      O--+
                         O
            sage: T.poset().cover_relations()
            [[3, 4], [3, 2], [4, 5], [2, 5], [1, 5]]
        """
        n = self.size()
        M = [[' O ' if i == j else '   ' for i in range(n)] for j in range(n)]

        def superpose(x, y, b):
            # put symbol b at position x, y
            # on top of existing symbols there
            i = x - 1
            j = y - 1
            a = M[i][j]
            if a == '   ':
                M[i][j] = b
            elif a == '-+ ':
                if b == a:
                    pass
                elif b == '---':
                    M[i][j] = '-+-'
                elif b == ' | ':
                    M[i][j] = '-+ '
            elif a == ' +-':
                if b == a:
                    pass
                elif b == '---':
                    M[i][j] = '-+-'
                elif b == ' | ':
                    M[i][j] = ' +-'
            elif a == '---':
                if b == a:
                    pass
                elif b == '-+ ':
                    M[i][j] = '-+-'
                elif b == ' +-':
                    M[i][j] = '-+-'
            elif a == ' | ':
                if b == a:
                    pass
                elif b == '-+ ':
                    M[i][j] = '-+ '
                elif b == ' +-':
                    M[i][j] = ' +-'

        def superpose_node(i, right=True):
            i -= 1  # for indexing
            if M[i][i] == ' O ':
                if right:
                    M[i][i] = ' O-'
                else:
                    M[i][i] = '-O '
            elif M[i][i] == ' O-' and not right:
                M[i][i] = '-O-'
            elif M[i][i] == '-O ' and right:
                M[i][i] = '-O-'

        for i, j in self.poset().hasse_diagram().edges(labels=False):
            if i > j:
                superpose_node(i, False)
                superpose(i, j, ' +-')
                for k in range(j + 1, i):
                    superpose(k, j, ' | ')
                    superpose(i, k, '---')
            else:
                superpose_node(i, True)
                superpose(i, j, '-+ ')
                for k in range(i + 1, j):
                    superpose(i, k, '---')
                    superpose(k, j, ' | ')

        from sage.typeset.ascii_art import AsciiArt
        return AsciiArt([''.join(ligne) for ligne in M])

    def _unicode_art_(self):
        """
        Return an unicode picture of ``self``.

        This is a picture of the Hasse diagram. Vertices from `1` to `n` are
        placed on the diagonal from top-left to bottom-right.
        Then increasing covers are drawn above the diagonal
        and decreasing covers are drawn below the diagonal.

        EXAMPLES::

            sage: T = TamariIntervalPosets(5)[56]
            sage: unicode_art(T)
            o───╮
             o──┤
             ╰o╮│
               o┤
                o
            sage: T.poset().cover_relations()
            [[3, 4], [3, 2], [4, 5], [2, 5], [1, 5]]
        """
        n = self.size()
        M = [['o' if i == j else ' ' for i in range(n)] for j in range(n)]

        def superpose(x, y, b):
            # put symbol b at position x, y
            # on top of existing symbols there
            i = x - 1
            j = y - 1
            a = M[i][j]
            if a == ' ':
                M[i][j] = b
            elif a == '╮':
                if b == a:
                    pass
                elif b == '─':
                    M[i][j] = '┬'
                elif b == '│':
                    M[i][j] = '┤'
            elif a == '╰':
                if b == a:
                    pass
                elif b == '─':
                    M[i][j] = '┴'
                elif b == '│':
                    M[i][j] = '├'
            elif a == '─':
                if b == a:
                    pass
                elif b == '╮':
                    M[i][j] = '┬'
                elif b == '╰':
                    M[i][j] = '┴'
            elif a == '│':
                if b == a:
                    pass
                elif b == '╮':
                    M[i][j] = '┤'
                elif b == '╰':
                    M[i][j] = '├'

        for i, j in self.poset().hasse_diagram().edges(labels=False):
            if i > j:
                superpose(i, j, '╰')
                for k in range(j + 1, i):
                    superpose(k, j, '│')
                    superpose(i, k, '─')
            else:
                superpose(i, j, '╮')
                for k in range(i + 1, j):
                    superpose(i, k, '─')
                    superpose(k, j, '│')

        from sage.typeset.unicode_art import UnicodeArt
        return UnicodeArt([''.join(ligne) for ligne in M])

    def _richcmp_(self, other, op) -> bool:
        r"""
        Comparison.

        The comparison is first by size, then
        using cubical coordinates.

        .. SEEALSO:: :meth:`cubical_coordinates`

        TESTS::

            sage: TamariIntervalPoset(0,[]) == TamariIntervalPoset(0,[])
            True
            sage: TamariIntervalPoset(1,[]) == TamariIntervalPoset(0,[])
            False
            sage: TamariIntervalPoset(3,[(1,2),(3,2)]) == TamariIntervalPoset(3,[(3,2),(1,2)])
            True
            sage: TamariIntervalPoset(3,[(1,2),(3,2)]) == TamariIntervalPoset(3,[(1,2)])
            False
            sage: TamariIntervalPoset(3,[(1,2),(3,2)]) != TamariIntervalPoset(3,[(3,2),(1,2)])
            False

            sage: ip1 = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)])
            sage: ip2 = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: ip1 <= ip2
            False
            sage: ip1 <= ip1
            True
            sage: ip2 <= ip1
            True

            sage: ip1 != 33
            True
        """
        if not isinstance(other, TamariIntervalPoset):
            return NotImplemented
        if op == op_EQ:
            return (self.size() == other.size() and
                    self._cover_relations == other._cover_relations)
        if op == op_NE:
            return not(self.size() == other.size() and
                       self._cover_relations == other._cover_relations)
        return richcmp((self.size(), self.cubical_coordinates()),
                       (other.size(), other.cubical_coordinates()), op)

    def __iter__(self) -> Iterator[int]:
        r"""
        Iterate through the vertices of ``self``.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(1,2),(3,2)])
            sage: [i for i in ip]
            [1, 2, 3, 4]
        """
        return iter(range(1, self.size() + 1))

    def contains_interval(self, other) -> bool:
        r"""
        Return whether the interval represented by ``other`` is contained
        in ``self`` as an interval of the Tamari lattice.

        In terms of interval-posets, it means that all relations of ``self``
        are relations of ``other``.

        INPUT:

        - ``other`` -- an interval-poset

        EXAMPLES::

            sage: ip1 = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)])
            sage: ip2 = TamariIntervalPoset(4,[(2,3)])
            sage: ip2.contains_interval(ip1)
            True
            sage: ip3 = TamariIntervalPoset(4,[(2,1)])
            sage: ip2.contains_interval(ip3)
            False
            sage: ip4 = TamariIntervalPoset(3,[(2,3)])
            sage: ip2.contains_interval(ip4)
            False
        """
        if other.size() != self.size():
            return False
        return all(other.le(i, j) for (i, j) in self._cover_relations)

    def lower_contains_interval(self, other) -> bool:
        r"""
        Return whether the interval represented by ``other`` is contained
        in ``self`` as an interval of the Tamari lattice and if they share
        the same lower bound.

        As interval-posets, it means that ``other`` contains the relations
        of ``self`` plus some extra increasing relations.

        INPUT:

        - ``other`` -- an interval-poset

        EXAMPLES::

            sage: ip1 = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)])
            sage: ip2 = TamariIntervalPoset(4,[(4,3)])
            sage: ip2.lower_contains_interval(ip1)
            True
            sage: ip2.contains_interval(ip1) and ip2.lower_binary_tree() == ip1.lower_binary_tree()
            True
            sage: ip3 = TamariIntervalPoset(4,[(4,3),(2,1)])
            sage: ip2.contains_interval(ip3)
            True
            sage: ip2.lower_binary_tree() == ip3.lower_binary_tree()
            False
            sage: ip2.lower_contains_interval(ip3)
            False

        TESTS::

            sage: ip1 = TamariIntervalPoset(3,[(1,2),(2,3)])
            sage: ip2 = TamariIntervalPoset(3,[(2,1),(3,2)])
            sage: ip2.lower_contains_interval(ip1)
            False
        """
        if not self.contains_interval(other):
            return False
        return all(self.le(i, j)
                   for (i, j) in other.decreasing_cover_relations())

    def upper_contains_interval(self, other) -> bool:
        r"""
        Return whether the interval represented by ``other`` is contained
        in ``self`` as an interval of the Tamari lattice and if they share
        the same upper bound.

        As interval-posets, it means that ``other`` contains the relations
        of ``self`` plus some extra decreasing relations.

        INPUT:

        - ``other`` -- an interval-poset

        EXAMPLES::

            sage: ip1 = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)])
            sage: ip2 = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: ip2.upper_contains_interval(ip1)
            True
            sage: ip2.contains_interval(ip1) and ip2.upper_binary_tree() == ip1.upper_binary_tree()
            True
            sage: ip3 = TamariIntervalPoset(4,[(1,2),(2,3),(3,4)])
            sage: ip2.upper_contains_interval(ip3)
            False
            sage: ip2.contains_interval(ip3)
            True
            sage: ip2.upper_binary_tree() == ip3.upper_binary_tree()
            False

        TESTS::

            sage: ip1 = TamariIntervalPoset(3,[(1,2),(2,3)])
            sage: ip2 = TamariIntervalPoset(3,[(2,1),(3,2)])
            sage: ip2.lower_contains_interval(ip1)
            False
        """
        if not self.contains_interval(other):
            return False
        return all(self.le(i, j)
                   for (i, j) in other.increasing_cover_relations())

    def is_linear_extension(self, perm) -> bool:
        r"""
        Return whether the permutation ``perm`` is a linear extension
        of ``self``.

        INPUT:

        - ``perm`` -- a permutation of the size of ``self``

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)])
            sage: ip.is_linear_extension([1,4,2,3])
            True
            sage: ip.is_linear_extension(Permutation([1,4,2,3]))
            True
            sage: ip.is_linear_extension(Permutation([1,4,3,2]))
            False
        """
        return self._poset.is_linear_extension(perm)

    def contains_binary_tree(self, binary_tree) -> bool:
        r"""
        Return whether the interval represented by ``self`` contains
        the binary tree ``binary_tree``.

        INPUT:

        - ``binary_tree`` -- a binary tree

        .. SEEALSO:: :meth:`contains_dyck_word`

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: ip.contains_binary_tree(BinaryTree([[None,[None,[]]],None]))
            True
            sage: ip.contains_binary_tree(BinaryTree([None,[[[],None],None]]))
            True
            sage: ip.contains_binary_tree(BinaryTree([[],[[],None]]))
            False
            sage: ip.contains_binary_tree(ip.lower_binary_tree())
            True
            sage: ip.contains_binary_tree(ip.upper_binary_tree())
            True
            sage: all(ip.contains_binary_tree(bt) for bt in ip.binary_trees())
            True
        """
        return self.is_linear_extension(binary_tree.to_132_avoiding_permutation())

    def contains_dyck_word(self, dyck_word) -> bool:
        r"""
        Return whether the interval represented by ``self`` contains
        the Dyck word ``dyck_word``.

        INPUT:

        - ``dyck_word`` -- a Dyck word

        .. SEEALSO:: :meth:`contains_binary_tree`

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: ip.contains_dyck_word(DyckWord([1,1,1,0,0,0,1,0]))
            True
            sage: ip.contains_dyck_word(DyckWord([1,1,0,1,0,1,0,0]))
            True
            sage: ip.contains_dyck_word(DyckWord([1,0,1,1,0,1,0,0]))
            False
            sage: ip.contains_dyck_word(ip.lower_dyck_word())
            True
            sage: ip.contains_dyck_word(ip.upper_dyck_word())
            True
            sage: all(ip.contains_dyck_word(bt) for bt in ip.dyck_words())
            True
        """
        return self.contains_binary_tree(dyck_word.to_binary_tree_tamari())

    def intersection(self, other: TIP) -> TIP:
        r"""
        Return the interval-poset formed by combining the relations from
        both ``self`` and ``other``. It corresponds to the intersection
        of the two corresponding intervals of the Tamari lattice.

        INPUT:

        - ``other`` -- an interval-poset of the same size as ``self``

        EXAMPLES::

            sage: ip1 = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: ip2 = TamariIntervalPoset(4,[(4,3)])
            sage: ip1.intersection(ip2)
            The Tamari interval of size 4 induced by relations [(1, 2), (2, 3), (4, 3)]
            sage: ip3 = TamariIntervalPoset(4,[(2,1)])
            sage: ip1.intersection(ip3)
            Traceback (most recent call last):
            ...
            ValueError: this intersection is empty, it does not correspond to an interval-poset
            sage: ip4 = TamariIntervalPoset(3,[(2,3)])
            sage: ip2.intersection(ip4)
            Traceback (most recent call last):
            ...
            ValueError: intersections are only possible on interval-posets of the same size
        """
        if other.size() != self.size():
            raise ValueError("intersections are only possible on interval-posets of the same size")
        try:
            return TamariIntervalPoset(self.size(), self._cover_relations + other._cover_relations)
        except ValueError:
            raise ValueError("this intersection is empty, it does not correspond to an interval-poset")

    def initial_forest(self) -> TIP:
        r"""
        Return the initial forest of ``self``, i.e., the interval-poset
        formed from only the increasing relations of ``self``.

        .. SEEALSO:: :meth:`final_forest`

        EXAMPLES::

            sage: TamariIntervalPoset(4,[(1,2),(3,2),(2,4),(3,4)]).initial_forest()
            The Tamari interval of size 4 induced by relations [(1, 2), (2, 4), (3, 4)]
            sage: ip = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: ip.initial_forest() == ip
            True
        """
        relations = self.increasing_cover_relations()
        P = FinitePoset(DiGraph([list(range(1, self._size + 1)), relations],
                                format='vertices_and_edges'))  # type:ignore
        return TamariIntervalPoset(P, check=False)  # type:ignore

    def final_forest(self) -> TIP:
        r"""
        Return the final forest of ``self``, i.e., the interval-poset
        formed with only the decreasing relations of ``self``.

        .. SEEALSO:: :meth:`initial_forest`

        EXAMPLES::

            sage: TamariIntervalPoset(4,[(2,1),(3,2),(3,4),(4,2)]).final_forest()
            The Tamari interval of size 4 induced by relations [(4, 2), (3, 2), (2, 1)]
            sage: ip = TamariIntervalPoset(3,[(2,1),(3,1)])
            sage: ip.final_forest() == ip
            True
        """
        relations = self.decreasing_cover_relations()
        P = FinitePoset(DiGraph([list(range(1, self._size + 1)), relations],
                                format='vertices_and_edges'))  # type:ignore
        return TamariIntervalPoset(P, check=False)  # type:ignore

    def is_initial_interval(self) -> bool:
        r"""
        Return if ``self`` corresponds to an initial interval of the Tamari
        lattice.

        This means that its lower end is the smallest element of the lattice.
        It consists of checking that ``self`` does not contain any decreasing
        relations.

        .. SEEALSO:: :meth:`is_final_interval`

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4, [(1, 2), (2, 4), (3, 4)])
            sage: ip.is_initial_interval()
            True
            sage: ip.lower_dyck_word()
            [1, 0, 1, 0, 1, 0, 1, 0]
            sage: ip = TamariIntervalPoset(4, [(1, 2), (2, 4), (3, 4), (3, 2)])
            sage: ip.is_initial_interval()
            False
            sage: ip.lower_dyck_word()
            [1, 0, 1, 1, 0, 0, 1, 0]
            sage: all(DyckWord([1,0,1,0,1,0]).tamari_interval(dw).is_initial_interval() for dw in DyckWords(3))
            True
        """
        return not self.decreasing_cover_relations()

    def is_final_interval(self) -> bool:
        r"""
        Return if ``self`` corresponds to a final interval of the Tamari
        lattice.

        This means that its upper end is the largest element of the lattice.
        It consists of checking that ``self`` does not contain any increasing
        relations.

        .. SEEALSO:: :meth:`is_initial_interval`

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4, [(4, 3), (3, 1), (2, 1)])
            sage: ip.is_final_interval()
            True
            sage: ip.upper_dyck_word()
            [1, 1, 1, 1, 0, 0, 0, 0]
            sage: ip = TamariIntervalPoset(4, [(4, 3), (3, 1), (2, 1), (2, 3)])
            sage: ip.is_final_interval()
            False
            sage: ip.upper_dyck_word()
            [1, 1, 0, 1, 1, 0, 0, 0]
            sage: all(dw.tamari_interval(DyckWord([1, 1, 1, 0, 0, 0])).is_final_interval() for dw in DyckWords(3))
            True
        """
        return not self.increasing_cover_relations()

    def lower_binary_tree(self):
        r"""
        Return the lowest binary tree in the interval of the Tamari
        lattice represented by ``self``.

        This is a binary tree. It is the shape of the unique binary
        search tree whose left-branch ordered forest (i.e., the result
        of applying
        :meth:`~sage.combinat.binary_tree.BinaryTree.to_ordered_tree_left_branch`
        and cutting off the root) is the final forest of ``self``.

        .. SEEALSO:: :meth:`lower_dyck_word`

        EXAMPLES::

            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(4,5)]); ip
            The Tamari interval of size 6 induced by relations [(1, 2), (4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]
            sage: ip.lower_binary_tree()
            [[., .], [[., [., .]], [., .]]]
            sage: TamariIntervalPosets.final_forest(ip.lower_binary_tree()) == ip.final_forest()
            True
            sage: ip == TamariIntervalPosets.from_binary_trees(ip.lower_binary_tree(),ip.upper_binary_tree())
            True
        """
        return self.min_linear_extension().binary_search_tree_shape(left_to_right=False)

    def lower_dyck_word(self):
        r"""
        Return the lowest Dyck word in the interval of the Tamari lattice
        represented by ``self``.

        .. SEEALSO:: :meth:`lower_binary_tree`

        EXAMPLES::

            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(4,5)]); ip
            The Tamari interval of size 6 induced by relations [(1, 2), (4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]
            sage: ip.lower_dyck_word()
            [1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0]
            sage: TamariIntervalPosets.final_forest(ip.lower_dyck_word()) == ip.final_forest()
            True
            sage: ip == TamariIntervalPosets.from_dyck_words(ip.lower_dyck_word(),ip.upper_dyck_word())
            True
        """
        return self.lower_binary_tree().to_dyck_word_tamari()

    def upper_binary_tree(self):
        r"""
        Return the highest binary tree in the interval of the Tamari
        lattice represented by ``self``.

        This is a binary tree. It is the shape of the unique binary
        search tree whose right-branch ordered forest (i.e., the result
        of applying
        :meth:`~sage.combinat.binary_tree.BinaryTree.to_ordered_tree_right_branch`
        and cutting off the root) is the initial forest of ``self``.

        .. SEEALSO:: :meth:`upper_dyck_word`

        EXAMPLES::

            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(4,5)]); ip
            The Tamari interval of size 6 induced by relations [(1, 2), (4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]
            sage: ip.upper_binary_tree()
            [[., .], [., [[., .], [., .]]]]
            sage: TamariIntervalPosets.initial_forest(ip.upper_binary_tree()) == ip.initial_forest()
            True
            sage: ip == TamariIntervalPosets.from_binary_trees(ip.lower_binary_tree(),ip.upper_binary_tree())
            True
        """
        return self.max_linear_extension().binary_search_tree_shape(left_to_right=False)

    def upper_dyck_word(self):
        r"""
        Return the highest Dyck word in the interval of the Tamari lattice
        represented by ``self``.

        .. SEEALSO:: :meth:`upper_binary_tree`

        EXAMPLES::

            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(4,5)]); ip
            The Tamari interval of size 6 induced by relations [(1, 2), (4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]
            sage: ip.upper_dyck_word()
            [1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0]
            sage: TamariIntervalPosets.initial_forest(ip.upper_dyck_word()) == ip.initial_forest()
            True
            sage: ip == TamariIntervalPosets.from_dyck_words(ip.lower_dyck_word(),ip.upper_dyck_word())
            True
        """
        return self.upper_binary_tree().to_dyck_word_tamari()

    def subposet(self, start, end) -> TIP:
        r"""
        Return the renormalized subposet of ``self`` consisting solely
        of integers from ``start`` (inclusive) to ``end`` (not inclusive).

        "Renormalized" means that these integers are relabelled
        `1,2,\ldots,k` in the obvious way (i.e., by subtracting
        ``start - 1``).

        INPUT:

        - ``start`` -- an integer, the starting vertex (inclusive)
        - ``end`` -- an integer, the ending vertex (not inclusive)

        EXAMPLES::

            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(3,5),(4,5)]); ip
            The Tamari interval of size 6 induced by relations [(1, 2), (3, 5), (4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]
            sage: ip.subposet(1,3)
            The Tamari interval of size 2 induced by relations [(1, 2)]
            sage: ip.subposet(1,4)
            The Tamari interval of size 3 induced by relations [(1, 2), (3, 2)]
            sage: ip.subposet(1,5)
            The Tamari interval of size 4 induced by relations [(1, 2), (4, 3), (3, 2)]
            sage: ip.subposet(1,7) == ip
            True
            sage: ip.subposet(1,1)
            The Tamari interval of size 0 induced by relations []

        TESTS::

            sage: ip.sub_poset(1,1)
            The Tamari interval of size 0 induced by relations []
            sage: ip = TamariIntervalPosets(4).an_element()
            sage: ip.subposet(2,9)
            Traceback (most recent call last):
            ...
            ValueError: invalid starting or ending value
        """
        if start < 1 or start > end or end > self.size() + 1:
            raise ValueError("invalid starting or ending value")
        if start == end:
            return TamariIntervalPoset(0, [])
        relations = [(i - start + 1, j - start + 1)
                     for (i, j) in self.increasing_cover_relations()
                     if i >= start and j < end]
        relations.extend([(j - start + 1, i - start + 1)
                          for (j, i) in self.decreasing_cover_relations()
                          if i >= start and j < end])
        return TamariIntervalPoset(end - start, relations, check=False)

    sub_poset = subposet

    def min_linear_extension(self) -> Permutation:
        r"""
        Return the minimal permutation for the right weak order which is
        a linear extension of ``self``.

        This is also the minimal permutation in the sylvester
        class of ``self.lower_binary_tree()`` and is a 312-avoiding
        permutation.

        The right weak order is also known as the right permutohedron
        order. See
        :meth:`~sage.combinat.permutation.Permutation.permutohedron_lequal`
        for its definition.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)])
            sage: ip.min_linear_extension()
            [1, 2, 4, 3]
            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(4,5)])
            sage: ip.min_linear_extension()
            [1, 4, 3, 6, 5, 2]
            sage: ip = TamariIntervalPoset(0,[])
            sage: ip.min_linear_extension()
            []
            sage: ip = TamariIntervalPoset(5, [(1, 4), (2, 4), (3, 4), (5, 4)]); ip
            The Tamari interval of size 5 induced by relations [(1, 4), (2, 4), (3, 4), (5, 4)]
            sage: ip.min_linear_extension()
            [1, 2, 3, 5, 4]
        """
        # The min linear extension is built by postfix-reading the
        # final forest of ``self``.
        final_forest = DiGraph([list(self),
                                self.decreasing_cover_relations()],
                               format='vertices_and_edges')

        def add(perm: list, i):
            r"""
            Internal recursive method to compute the min linear extension.
            """
            for j in sorted(final_forest.neighbors_in(i)):
                add(perm, j)
            perm.append(i)
        perm: list[int] = []
        for i in sorted(final_forest.sinks()):
            add(perm, i)
        return Permutation(perm)  # type:ignore

    def max_linear_extension(self) -> Permutation:
        r"""
        Return the maximal permutation for the right weak order which is
        a linear extension of ``self``.

        This is also the maximal permutation in the sylvester
        class of ``self.upper_binary_tree()`` and is a 132-avoiding
        permutation.

        The right weak order is also known as the right permutohedron
        order. See
        :meth:`~sage.combinat.permutation.Permutation.permutohedron_lequal`
        for its definition.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)])
            sage: ip.max_linear_extension()
            [4, 1, 2, 3]
            sage: ip = TamariIntervalPoset(6,[(3,2),(4,3),(5,2),(6,5),(1,2),(4,5)]); ip
            The Tamari interval of size 6 induced by relations [(1, 2), (4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]
            sage: ip.max_linear_extension()
            [6, 4, 5, 3, 1, 2]
            sage: ip = TamariIntervalPoset(0,[]); ip
            The Tamari interval of size 0 induced by relations []
            sage: ip.max_linear_extension()
            []
            sage: ip = TamariIntervalPoset(5, [(1, 4), (2, 4), (3, 4), (5, 4)]); ip
            The Tamari interval of size 5 induced by relations [(1, 4), (2, 4), (3, 4), (5, 4)]
            sage: ip.max_linear_extension()
            [5, 3, 2, 1, 4]
        """
        # The max linear extension is built by right-to-left
        # postfix-reading the initial forest of ``self``.
        initial_forest = DiGraph([list(self),
                                  self.increasing_cover_relations()],
                                 format='vertices_and_edges')

        def add(perm: list, i):
            r"""
            Internal recursive method to compute the max linear extension.
            """
            for j in sorted(initial_forest.neighbors_in(i), reverse=True):
                add(perm, j)
            perm.append(i)
        perm: list[int] = []
        for i in sorted(initial_forest.sinks(), reverse=True):
            add(perm, i)
        return Permutation(perm)  # type:ignore

    def linear_extensions(self) -> Iterator[Permutation]:
        r"""
        Return an iterator on the permutations which are linear
        extensions of ``self``.

        They form an interval of the right weak order (also called the
        right permutohedron order -- see
        :meth:`~sage.combinat.permutation.Permutation.permutohedron_lequal`
        for a definition).

        EXAMPLES::

            sage: ip = TamariIntervalPoset(3,[(1,2),(3,2)])
            sage: list(ip.linear_extensions())
            [[3, 1, 2], [1, 3, 2]]
            sage: ip = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)])
            sage: list(ip.linear_extensions())
            [[4, 1, 2, 3], [1, 2, 4, 3], [1, 4, 2, 3]]
        """
        for ext in self._poset.linear_extensions():
            yield Permutation(ext)  # type:ignore

    def lower_contained_intervals(self) -> Iterator[TIP]:
        r"""
        If ``self`` represents the interval `[t_1, t_2]` of the Tamari
        lattice, return an iterator on all intervals `[t_1,t]` with
        `t \leq t_2` for the Tamari lattice.

        In terms of interval-posets, it corresponds to adding all possible
        relations of the form `n` precedes `m` with `n<m`.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: list(ip.lower_contained_intervals())
            [The Tamari interval of size 4 induced by relations [(2, 4), (3, 4), (3, 1), (2, 1)],
             The Tamari interval of size 4 induced by relations [(1, 4), (2, 4), (3, 4), (3, 1), (2, 1)],
             The Tamari interval of size 4 induced by relations [(2, 3), (3, 4), (3, 1), (2, 1)],
             The Tamari interval of size 4 induced by relations [(1, 4), (2, 3), (3, 4), (3, 1), (2, 1)]]
            sage: ip = TamariIntervalPoset(4,[])
            sage: len(list(ip.lower_contained_intervals()))
            14
        """
        size = self._size
        yield self
        r"""
        we try to add links recursively in this order:
        1 -> 2
        2 -> 3
        1 -> 3
        3 -> 4
        2 -> 4
        1 -> 4
        ...
        ("Link" means "relation of the poset".)

        One useful feature of interval-posets is that if you add a single
        new relation -- say, `x` precedes `y` -- to an existing
        interval-poset and take the transitive closure, and if the axioms
        of an interval-poset are still satisfied for `(a,c) = (x,y)` and
        for `(a,c) = (y,x)`, then the transitive closure is an
        interval-poset (i.e., roughly speaking, the other new relations
        forced by `x` preceding `y` under transitive closure cannot
        invalidate the axioms). This is helpful when extending
        interval-posets, and is the reason why this and other iterators
        don't yield invalid interval-posets.
        """
        def add_relations(poset, n, m):
            r"""
            Internal recursive method to generate all possible intervals.
            At every step during the iteration, we have n < m and every
            i satisfying n < i < m satisfies that i precedes m in the
            poset ``poset`` (except when m > size).
            """
            if n <= 0:
                # if n<=0, then we go to the next m
                n = m
                m += 1
            if m > size:
                # if m>size, it's finished
                return

            if poset.le(n, m):
                # there is already a link n->m, so we go to the next n
                for pos in add_relations(poset, n - 1, m):
                    yield pos
            elif poset.le(m, n):
                # there is an inverse link m->n, we know we won't be able
                # to create a link i->m with i<=n, so we go to the next m
                for pos in add_relations(poset, m, m + 1):
                    yield pos
            else:
                # there is no link n->m
                # first option : we don't create the link and go to the next m
                # (since the lack of a link n->m forbids any links i->m
                # with i<n)
                for pos in add_relations(poset, m, m + 1):
                    yield pos
                # second option : we create the link
                # (this is allowed because links i->m already exist for all
                # n<i<m, or else we wouldn't be here)
                poset = TamariIntervalPoset(poset.size(), poset._cover_relations + ((n, m),))
                yield poset
                # and then, we go to the next n
                for pos in add_relations(poset, n - 1, m):
                    yield pos

        for inter in add_relations(self, 1, 2):
            yield inter

    def interval_cardinality(self) -> Integer:
        r"""
        Return the cardinality of the interval, i.e., the number of elements
        (binary trees or Dyck words) in the interval represented by ``self``.

        Not to be confused with :meth:`size` which is the number of
        vertices.

        .. SEEALSO:: :meth:`binary_trees`

        EXAMPLES::

            sage: TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)]).interval_cardinality()
            4
            sage: TamariIntervalPoset(4,[]).interval_cardinality()
            14
            sage: TamariIntervalPoset(4,[(1,2),(2,3),(3,4)]).interval_cardinality()
            1
        """
        return len(list(self.lower_contained_intervals()))

    def binary_trees(self) -> Iterator:
        r"""
        Return an iterator on all the binary trees in the interval
        represented by ``self``.

        .. SEEALSO:: :meth:`interval_cardinality`

        EXAMPLES::

            sage: list(TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)]).binary_trees())
            [[., [[., [., .]], .]],
             [[., [., [., .]]], .],
             [., [[[., .], .], .]],
             [[., [[., .], .]], .]]
            sage: set(TamariIntervalPoset(4,[]).binary_trees()) == set(BinaryTrees(4))
            True
        """
        for ip in self.lower_contained_intervals():
            yield ip.upper_binary_tree()

    def dyck_words(self) -> Iterator:
        r"""
        Return an iterator on all the Dyck words in the interval
        represented by ``self``.

        EXAMPLES::

            sage: list(TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)]).dyck_words())
            [[1, 1, 1, 0, 0, 1, 0, 0],
             [1, 1, 1, 0, 0, 0, 1, 0],
             [1, 1, 0, 1, 0, 1, 0, 0],
             [1, 1, 0, 1, 0, 0, 1, 0]]
            sage: set(TamariIntervalPoset(4,[]).dyck_words()) == set(DyckWords(4))
            True
        """
        for ip in self.lower_contained_intervals():
            yield ip.upper_dyck_word()

    def maximal_chain_tamari_intervals(self) -> Iterator[TIP]:
        r"""
        Return an iterator on the upper contained intervals of one
        longest chain of the Tamari interval represented by ``self``.

        If ``self`` represents the interval `[T_1,T_2]` of the Tamari
        lattice, this returns intervals `[T',T_2]` with `T'` following
        one longest chain between `T_1` and `T_2`.

        To obtain a longest chain, we use the Tamari inversions of ``self``.
        The elements of the chain are obtained by adding one by one the
        relations `(b,a)` from each Tamari inversion `(a,b)` to ``self``,
        where the Tamari inversions are taken in lexicographic order.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: list(ip.maximal_chain_tamari_intervals())
            [The Tamari interval of size 4 induced by relations [(2, 4), (3, 4), (3, 1), (2, 1)],
             The Tamari interval of size 4 induced by relations [(2, 4), (3, 4), (4, 1), (3, 1), (2, 1)],
             The Tamari interval of size 4 induced by relations [(2, 4), (3, 4), (4, 1), (3, 2), (2, 1)]]
            sage: ip = TamariIntervalPoset(4,[])
            sage: list(ip.maximal_chain_tamari_intervals())
            [The Tamari interval of size 4 induced by relations [],
             The Tamari interval of size 4 induced by relations [(2, 1)],
             The Tamari interval of size 4 induced by relations [(3, 1), (2, 1)],
             The Tamari interval of size 4 induced by relations [(4, 1), (3, 1), (2, 1)],
             The Tamari interval of size 4 induced by relations [(4, 1), (3, 2), (2, 1)],
             The Tamari interval of size 4 induced by relations [(4, 2), (3, 2), (2, 1)],
             The Tamari interval of size 4 induced by relations [(4, 3), (3, 2), (2, 1)]]
        """
        yield self
        n = self.size()
        cover_relations = list(self._cover_relations)
        for inv in self.tamari_inversions_iter():
            cover_relations.append((inv[1], inv[0]))
            yield TamariIntervalPoset(n, cover_relations, check=False)

    def maximal_chain_binary_trees(self) -> Iterator:
        r"""
        Return an iterator on the binary trees forming a longest chain of
        ``self`` (regarding ``self`` as an interval of the Tamari
        lattice).

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: list(ip.maximal_chain_binary_trees())
            [[[., [[., .], .]], .], [., [[[., .], .], .]], [., [[., [., .]], .]]]
            sage: ip = TamariIntervalPoset(4,[])
            sage: list(ip.maximal_chain_binary_trees())
            [[[[[., .], .], .], .],
             [[[., [., .]], .], .],
             [[., [[., .], .]], .],
             [., [[[., .], .], .]],
             [., [[., [., .]], .]],
             [., [., [[., .], .]]],
             [., [., [., [., .]]]]]
        """
        for it in self.maximal_chain_tamari_intervals():
            yield it.lower_binary_tree()

    def maximal_chain_dyck_words(self) -> Iterator:
        r"""
        Return an iterator on the Dyck words forming a longest chain of
        ``self`` (regarding ``self`` as an interval of the Tamari
        lattice).

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: list(ip.maximal_chain_dyck_words())
            [[1, 1, 0, 1, 0, 0, 1, 0], [1, 1, 0, 1, 0, 1, 0, 0], [1, 1, 1, 0, 0, 1, 0, 0]]
            sage: ip = TamariIntervalPoset(4,[])
            sage: list(ip.maximal_chain_dyck_words())
            [[1, 0, 1, 0, 1, 0, 1, 0],
             [1, 1, 0, 0, 1, 0, 1, 0],
             [1, 1, 0, 1, 0, 0, 1, 0],
             [1, 1, 0, 1, 0, 1, 0, 0],
             [1, 1, 1, 0, 0, 1, 0, 0],
             [1, 1, 1, 0, 1, 0, 0, 0],
             [1, 1, 1, 1, 0, 0, 0, 0]]
        """
        for it in self.maximal_chain_tamari_intervals():
            yield it.lower_dyck_word()

    def tamari_inversions(self) -> list[tuple[int, int]]:
        r"""
        Return the Tamari inversions of ``self``.

        A Tamari inversion is
        a pair of vertices `(a,b)` with `a < b` such that:

        - the decreasing parent of `b` is strictly smaller than `a` (or
          does not exist), and
        - the increasing parent of `a` is strictly bigger than `b` (or
          does not exist).

        "Smaller" and "bigger" refer to the numerical values of the
        elements, not to the poset order.

        This method returns the list of all Tamari inversions in
        lexicographic order.

        The number of Tamari inversions is the length of the
        longest chain of the Tamari interval represented by ``self``.

        Indeed, when an interval consists of just one binary tree, it has
        no inversion. One can also prove that if a Tamari interval
        `I' = [T_1', T_2']` is a proper subset of a Tamari interval
        `I = [T_1, T_2]`, then the inversion number of `I'` is strictly
        lower than the inversion number of `I`. And finally, by adding
        the relation `(b,a)` to the interval-poset where `(a,b)` is the
        first inversion of `I` in lexicographic order, one reduces the
        inversion number by exactly `1`.

        .. SEEALSO::

            :meth:`tamari_inversions_iter`, :meth:`number_of_tamari_inversions`

        EXAMPLES::

            sage: ip = TamariIntervalPoset(3,[])
            sage: ip.tamari_inversions()
            [(1, 2), (1, 3), (2, 3)]
            sage: ip = TamariIntervalPoset(3,[(2,1)])
            sage: ip.tamari_inversions()
            [(1, 3), (2, 3)]
            sage: ip = TamariIntervalPoset(3,[(1,2)])
            sage: ip.tamari_inversions()
            [(2, 3)]
            sage: ip = TamariIntervalPoset(3,[(1,2),(3,2)])
            sage: ip.tamari_inversions()
            []
            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: ip.tamari_inversions()
            [(1, 4), (2, 3)]
            sage: ip = TamariIntervalPoset(4,[])
            sage: ip.tamari_inversions()
            [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
            sage: all(len(TamariIntervalPosets.from_binary_trees(bt,bt).tamari_inversions())==0 for bt in BinaryTrees(3))
            True
            sage: all(len(TamariIntervalPosets.from_binary_trees(bt,bt).tamari_inversions())==0 for bt in BinaryTrees(4))
            True
        """
        return list(self.tamari_inversions_iter())

    def tamari_inversions_iter(self) -> Iterator[tuple[int, int]]:
        r"""
        Iterate over the Tamari inversions of ``self``, in
        lexicographic order.

        See :meth:`tamari_inversions` for the definition of the terms
        involved.

        EXAMPLES::

            sage: T = TamariIntervalPoset(5, [[1,2],[3,4],[3,2],[5,2],[4,2]])
            sage: list(T.tamari_inversions_iter())
            [(4, 5)]

            sage: T = TamariIntervalPoset(8, [(2, 7), (3, 7), (4, 7), (5, 7), (6, 7), (8, 7), (6, 4), (5, 4), (4, 3), (3, 2)])
            sage: list(T.tamari_inversions_iter())
            [(1, 2), (1, 7), (5, 6)]

            sage: T = TamariIntervalPoset(1, [])
            sage: list(T.tamari_inversions_iter())
            []

            sage: T = TamariIntervalPoset(0, [])
            sage: list(T.tamari_inversions_iter())
            []
        """
        final_forest = DiGraph([list(self),
                                self.decreasing_cover_relations()],
                               format='vertices_and_edges')
        initial_forest = DiGraph([list(self),
                                  self.increasing_cover_relations()],
                                 format='vertices_and_edges')
        n1 = self.size() + 1
        for a in range(1, self.size()):   # a == n will never work
            try:
                ipa = next(initial_forest.neighbor_out_iterator(a))
                max_b_1 = ipa
            except StopIteration:
                max_b_1 = n1
            for b in range(a + 1, max_b_1):
                try:
                    dpb = next(final_forest.neighbor_out_iterator(b))
                    if dpb < a:
                        yield (a, b)
                except StopIteration:
                    yield (a, b)

    def number_of_tamari_inversions(self) -> Integer:
        r"""
        Return the number of Tamari inversions of ``self``.

        This is also the length the longest chain of the Tamari
        interval represented by ``self``.

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(2,4),(3,4),(2,1),(3,1)])
            sage: ip.number_of_tamari_inversions()
            2
            sage: ip = TamariIntervalPoset(4,[])
            sage: ip.number_of_tamari_inversions()
            6
            sage: ip = TamariIntervalPoset(3,[])
            sage: ip.number_of_tamari_inversions()
            3
        """
        return len(self.tamari_inversions())

    def number_of_new_components(self) -> Integer:
        """
        Return the number of terms in the decomposition in new interval-posets.

        Every interval-poset has a unique decomposition as a planar tree
        of new interval-posets, as explained in [Cha2008]_. This function
        just computes the number of terms, not the planar tree nor
        the terms themselves.

        .. SEEALSO:: :meth:`is_new`, :meth:`new_decomposition`

        EXAMPLES::

            sage: TIP4 = TamariIntervalPosets(4)
            sage: nb = [u.number_of_new_components() for u in TIP4]
            sage: [nb.count(i) for i in range(1, 5)]
            [12, 21, 21, 14]
        """
        t_low = self.lower_binary_tree().to_tilting()
        t_up = self.upper_binary_tree().to_tilting()
        return sum(1 for p in t_low if p in t_up)

    def new_decomposition(self) -> list[TIP]:
        """
        Return the decomposition of the interval-poset into
        new interval-posets.

        Every interval-poset has a unique decomposition as a planar
        tree of new interval-posets, as explained in
        [Cha2008]_. This function computes the terms of this
        decomposition, but not the planar tree.

        For the number of terms, you can use instead the method
        :meth:`number_of_new_components`.

        OUTPUT:

        a list of new interval-posets.

        .. SEEALSO::

            :meth:`number_of_new_components`, :meth:`is_new`

        EXAMPLES::

            sage: ex = TamariIntervalPosets(4)[11]
            sage: ex.number_of_new_components()
            3
            sage: ex.new_decomposition()
            [The Tamari interval of size 1 induced by relations [],
             The Tamari interval of size 2 induced by relations [],
             The Tamari interval of size 1 induced by relations []]

        TESTS::

            sage: ex = TamariIntervalPosets(4).random_element()
            sage: dec = ex.new_decomposition()
            sage: len(dec) == ex.number_of_new_components()
            True
            sage: all(u.is_new() for u in dec)
            True
        """
        from sage.combinat.binary_tree import BinaryTree
        t_low = self.lower_binary_tree().to_tilting()
        t_up = self.upper_binary_tree().to_tilting()
        common = [p for p in t_low if p in t_up]

        def extract_tree(x, y, tilt, common):
            """
            Extract a tree with root at position xy (recursive).
            """
            left_tree = None
            for k in range(y - 1, x, -1):
                if (x, k) in tilt:
                    if (x, k) not in common:
                        left_tree = extract_tree(x, k, tilt, common)
                    break
            right_tree = None
            for k in range(x + 1, y):
                if (k, y) in tilt:
                    if (k, y) not in common:
                        right_tree = extract_tree(k, y, tilt, common)
                    break
            return BinaryTree([left_tree, right_tree], check=False)

        TIP = self.parent()
        return [TIP.from_binary_trees(extract_tree(cx, cy, t_low, common),
                                      extract_tree(cx, cy, t_up, common))
                for cx, cy in common]

    def decomposition_to_triple(self) -> None | tuple[TIP, TIP, int]:
        """
        Decompose an interval-poset into a triple (``left``, ``right``, ``r``).

        For the inverse method, see
        :meth:`TamariIntervalPosets.recomposition_from_triple`.

        OUTPUT:

        a triple (``left``, ``right``, ``r``) where ``left`` and
        ``right`` are interval-posets and ``r`` (an integer) is the
        parameter of the decomposition.

        EXAMPLES::

            sage: tip = TamariIntervalPoset(8, [(1,2), (2,4), (3,4), (6,7), (3,2), (5,4), (6,4), (8,7)])
            sage: tip.decomposition_to_triple()
            (The Tamari interval of size 3 induced by relations [(1, 2), (3, 2)],
            The Tamari interval of size 4 induced by relations [(2, 3), (4, 3)],
            2)
            sage: tip == TamariIntervalPosets.recomposition_from_triple(*tip.decomposition_to_triple())
            True

        TESTS::

            sage: tip = TamariIntervalPoset(0, [])
            sage: tip.decomposition_to_triple()

        REFERENCES:

        - [CP2015]_
        """
        n = self.size()
        if n == 0:
            return None
        root = self.increasing_roots()[-1]
        r = len(self.decreasing_children(root))
        left = self.subposet(1, root)
        right = self.subposet(root + 1, n + 1)
        return left, right, r

    def grafting_tree(self) -> LabelledBinaryTree:
        """
        Return the grafting tree of the interval-poset.

        For the inverse method,
        see :meth:`TamariIntervalPosets.from_grafting_tree`.

        EXAMPLES::

            sage: tip = TamariIntervalPoset(8, [(1,2), (2,4), (3,4), (6,7), (3,2), (5,4), (6,4), (8,7)])
            sage: tip.grafting_tree()
            2[1[0[., .], 0[., .]], 0[., 1[0[., .], 0[., .]]]]
            sage: tip == TamariIntervalPosets.from_grafting_tree(tip.grafting_tree())
            True

        REFERENCES:

        - [Pons2018]_
        """
        n = self.size()
        if n == 0:
            return LabelledBinaryTree(None)
        triplet = self.decomposition_to_triple()
        assert triplet is not None
        left, right, r = triplet
        return LabelledBinaryTree([left.grafting_tree(),
                                   right.grafting_tree()], label=r)

    def is_new(self) -> bool:
        """
        Return whether ``self`` is a new Tamari interval.

        Here 'new' means that the interval is not contained in any
        facet of the associahedron.
        This condition is invariant under complementation.

        They have been considered in section 9 of [Cha2008]_.

        .. SEEALSO:: :meth:`is_modern`

        EXAMPLES::

            sage: TIP4 = TamariIntervalPosets(4)
            sage: len([u for u in TIP4 if u.is_new()])
            12

            sage: TIP3 = TamariIntervalPosets(3)
            sage: len([u for u in TIP3 if u.is_new()])
            3
        """
        c_up = self.upper_binary_tree().single_edge_cut_shapes()
        c_down = self.lower_binary_tree().single_edge_cut_shapes()
        return not any(x in c_up for x in c_down)

    def is_simple(self) -> bool:
        """
        Return whether ``self`` is a simple Tamari interval.

        Here 'simple' means that the interval contains a unique binary tree.

        These intervals define the simple modules over the
        incidence algebras of the Tamari lattices.

        .. SEEALSO:: :meth:`is_final_interval`, :meth:`is_initial_interval`

        EXAMPLES::

            sage: TIP4 = TamariIntervalPosets(4)
            sage: len([u for u in TIP4 if u.is_simple()])
            14

            sage: TIP3 = TamariIntervalPosets(3)
            sage: len([u for u in TIP3 if u.is_simple()])
            5
        """
        return self.upper_binary_tree() == self.lower_binary_tree()

    def is_synchronized(self) -> bool:
        """
        Return whether ``self`` is a synchronized Tamari interval.

        This means that the upper and lower binary trees have the same canopee.
        This condition is invariant under complementation.

        This has been considered in [FPR2015]_. The numbers of
        synchronized intervals are given by the sequence :oeis:`A000139`.

        EXAMPLES::

            sage: len([T for T in TamariIntervalPosets(3)
            ....:     if T.is_synchronized()])
            6
        """
        up = self.upper_binary_tree()
        down = self.lower_binary_tree()
        return down.canopee() == up.canopee()

    def is_modern(self) -> bool:
        r"""
        Return whether ``self`` is a modern Tamari interval.

        This is defined by exclusion of a simple pattern in the Hasse diagram,
        namely there is no configuration `y \rightarrow x \leftarrow z`
        with `1 \leq y < x < z \leq n`.

        This condition is invariant under complementation.

        .. SEEALSO:: :meth:`is_new`, :meth:`is_infinitely_modern`

        EXAMPLES::

            sage: len([T for T in TamariIntervalPosets(3) if T.is_modern()])
            12

        REFERENCES:

        - [Rog2018]_
        """
        G = self.poset().hasse_diagram()
        for x in G:
            nx = G.neighbors_in(x)
            nx.append(x)
            if min(nx) < x < max(nx):
                return False
        return True

    def is_infinitely_modern(self) -> bool:
        r"""
        Return whether ``self`` is an infinitely-modern Tamari interval.

        This is defined by the exclusion of the configuration
        `i \rightarrow i + 1` and `j + 1 \rightarrow j` with `i < j`.

        This condition is invariant under complementation.

        .. SEEALSO:: :meth:`is_new`, :meth:`is_modern`

        EXAMPLES::

            sage: len([T for T in TamariIntervalPosets(3)
            ....:     if T.is_infinitely_modern()])
            12

        REFERENCES:

        - [Rog2018]_
        """
        n = self.size()
        found = False
        for i in range(1, n):
            if self.le(i, i + 1):
                found = True
                continue
            if self.le(i + 1, i) and found:
                return False
        return True

    def is_exceptional(self) -> bool:
        r"""
        Return whether ``self`` is an exceptional Tamari interval.

        This is defined by exclusion of a simple pattern in the Hasse diagram,
        namely there is no configuration ``y <-- x --> z``
        with `1 \leq y < x < z \leq n`.

        This condition is invariant under complementation.

        EXAMPLES::

            sage: len([T for T in TamariIntervalPosets(3)
            ....:     if T.is_exceptional()])
            12
        """
        G = self.poset().hasse_diagram()
        for x in G:
            nx = G.neighbors_out(x)
            nx.append(x)
            if min(nx) < x < max(nx):
                return False
        return True

    def is_dexter(self) -> bool:
        r"""
        Return whether ``self`` is a dexter Tamari interval.

        This is defined by an exclusion pattern in the Hasse diagram.
        See the code for the exact description.

        This condition is not invariant under complementation.

        EXAMPLES::

            sage: len([T for T in TamariIntervalPosets(3) if T.is_dexter()])
            12
        """
        G = self.poset().hasse_diagram()
        n = self.size()
        for x in range(2, n):
            nx = G.neighbors_out(x)
            nx.append(x)
            y = min(nx)
            if y < x and any(z > x for z in G.neighbor_out_iterator(y)):
                return False
        return True

    def is_indecomposable(self) -> bool:
        """
        Return whether ``self`` is an indecomposable Tamari interval.

        This is the terminology of [Cha2008]_.

        This means that the upper binary tree has an empty left subtree.

        This condition is not invariant under complementation.

        .. SEEALSO:: :meth:`is_connected`

        EXAMPLES::

            sage: len([T for T in TamariIntervalPosets(3)
            ....:      if T.is_indecomposable()])
            8
        """
        return not self.upper_binary_tree()[0]

    def is_connected(self) -> bool:
        """
        Return whether ``self`` is a connected Tamari interval.

        This means that the Hasse diagram is connected.

        This condition is invariant under complementation.

        .. SEEALSO:: :meth:`is_indecomposable`, :meth:`factor`

        EXAMPLES::

            sage: len([T for T in TamariIntervalPosets(3) if T.is_connected()])
            8
        """
        return self.poset().is_connected()


# Abstract class to serve as a Factory ; no instances are created.
class TamariIntervalPosets(UniqueRepresentation, Parent):
    r"""
    Factory for interval-posets.

    INPUT:

    - ``size`` -- (optional) an integer

    OUTPUT:

    - the set of all interval-posets (of the given ``size`` if specified)

    EXAMPLES::

        sage: TamariIntervalPosets()
        Interval-posets

        sage: TamariIntervalPosets(2)
        Interval-posets of size 2

    .. NOTE::

        This is a factory class whose constructor returns instances of
        subclasses.
    """
    @staticmethod
    def __classcall_private__(cls, n=None):
        r"""
        TESTS::

            sage: from sage.combinat.interval_posets import TamariIntervalPosets_all, TamariIntervalPosets_size
            sage: isinstance(TamariIntervalPosets(2), TamariIntervalPosets_size)
            True
            sage: isinstance(TamariIntervalPosets(), TamariIntervalPosets_all)
            True
            sage: TamariIntervalPosets(2) is TamariIntervalPosets_size(2)
            True
            sage: TamariIntervalPosets() is TamariIntervalPosets_all()
            True

            sage: TamariIntervalPosets(-2)
            Traceback (most recent call last):
            ...
            ValueError: n must be a nonnegative integer
        """
        if n is None:
            return TamariIntervalPosets_all()

        if n not in NN:
            raise ValueError("n must be a nonnegative integer")
        return TamariIntervalPosets_size(Integer(n))

    # add options to class
    class options(GlobalOptions):
        r"""
        Set and display the options for Tamari interval-posets.

        If no parameters are set, then the function returns a copy of
        the options dictionary.

        The ``options`` to Tamari interval-posets can be accessed as the method
        :meth:`TamariIntervalPosets.options` of :class:`TamariIntervalPosets`
        and related parent classes.

        @OPTIONS@

        EXAMPLES::

            sage: TIP = TamariIntervalPosets
            sage: TIP.options.latex_color_decreasing
            red
            sage: TIP.options.latex_color_decreasing='green'
            sage: TIP.options.latex_color_decreasing
            green
            sage: TIP.options._reset()
            sage: TIP.options.latex_color_decreasing
            red
        """
        NAME = 'TamariIntervalPosets'
        module = 'sage.combinat.interval_posets'
        latex_tikz_scale = dict(default=1,
                              description='the default value for the tikz scale when latexed',
                              checker=lambda x: True)  # More trouble than it's worth to check
        latex_line_width_scalar = dict(default=0.5,
                                     description='the default value for the line width as a'
                                                 'multiple of the tikz scale when latexed',
                                     checker=lambda x: True)  # More trouble than it's worth to check
        latex_color_decreasing = dict(default="red",
                                    description='the default color of decreasing relations when latexed',
                                    checker=lambda x: True)  # More trouble than it's worth to check
        latex_color_increasing = dict(default="blue",
                                    description='the default color of increasing relations when latexed',
                                    checker=lambda x: True)  # More trouble than it's worth to check
        latex_hspace = dict(default=1,
                          description='the default difference between horizontal'
                                      ' coordinates of vertices when latexed',
                          checker=lambda x: True)  # More trouble than it's worth to check
        latex_vspace = dict(default=1,
                          description='the default difference between vertical'
                                      ' coordinates of vertices when latexed',
                          checker=lambda x: True)   # More trouble than it's worth to check

    @staticmethod
    def check_poset(poset) -> bool:
        r"""
        Check if the given poset ``poset`` is a interval-poset, that is,
        if it satisfies the following properties:

        - Its labels are exactly `1, \ldots, n` where `n` is its size.
        - If `a < c` (as numbers) and `a` precedes `c`, then `b` precedes
          `c` for all `b` such that `a < b < c`.
        - If `a < c` (as numbers) and `c` precedes `a`, then `b` precedes
          `a` for all `b` such that `a < b < c`.

        INPUT:

        - ``poset`` -- a finite labeled poset

        EXAMPLES::

            sage: p = Poset(([1,2,3],[(1,2),(3,2)]))
            sage: TamariIntervalPosets.check_poset(p)
            True
            sage: p = Poset(([2,3],[(3,2)]))
            sage: TamariIntervalPosets.check_poset(p)
            False
            sage: p = Poset(([1,2,3],[(3,1)]))
            sage: TamariIntervalPosets.check_poset(p)
            False
            sage: p = Poset(([1,2,3],[(1,3)]))
            sage: TamariIntervalPosets.check_poset(p)
            False
        """
        if not set(poset._elements) == set(range(1, poset.cardinality() + 1)):
            return False

        for i in range(1, poset.cardinality() + 1):
            stop = False
            for j in range(i - 1, 0, -1):
                if not poset.le(j, i):
                    stop = True  # j does not precede i so no j'<j should
                elif stop:
                    return False
            stop = False
            for j in range(i + 1, poset.cardinality() + 1):
                if not poset.le(j, i):
                    stop = True  # j does not precede i so no j'>j should
                elif stop:
                    return False
        return True

    @staticmethod
    def final_forest(element) -> TIP:
        r"""
        Return the final forest of a binary tree, an interval-poset or a
        Dyck word.

        A final forest is an interval-poset corresponding to a final
        interval of the Tamari lattice, i.e., containing only decreasing
        relations.

        It can be constructed from a binary tree by its binary
        search tree labeling with the rule: `b` precedes
        `a` in the final forest iff `b` is in the right subtree of `a`
        in the binary search tree.

        INPUT:

        - ``element`` -- a binary tree, a Dyck word or an interval-poset

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)])
            sage: TamariIntervalPosets.final_forest(ip)
            The Tamari interval of size 4 induced by relations [(4, 3)]

        From binary trees::

            sage: bt = BinaryTree(); bt
            .
            sage: TamariIntervalPosets.final_forest(bt)
            The Tamari interval of size 0 induced by relations []
            sage: bt = BinaryTree([]); bt
            [., .]
            sage: TamariIntervalPosets.final_forest(bt)
            The Tamari interval of size 1 induced by relations []
            sage: bt = BinaryTree([[],None]); bt
            [[., .], .]
            sage: TamariIntervalPosets.final_forest(bt)
            The Tamari interval of size 2 induced by relations []
            sage: bt = BinaryTree([None,[]]); bt
            [., [., .]]
            sage: TamariIntervalPosets.final_forest(bt)
            The Tamari interval of size 2 induced by relations [(2, 1)]
            sage: bt = BinaryTree([[],[]]); bt
            [[., .], [., .]]
            sage: TamariIntervalPosets.final_forest(bt)
            The Tamari interval of size 3 induced by relations [(3, 2)]
            sage: bt = BinaryTree([[None,[[],None]],[]]); bt
            [[., [[., .], .]], [., .]]
            sage: TamariIntervalPosets.final_forest(bt)
            The Tamari interval of size 5 induced by relations [(5, 4), (3, 1), (2, 1)]

        From Dyck words::

            sage: dw = DyckWord([1,0])
            sage: TamariIntervalPosets.final_forest(dw)
            The Tamari interval of size 1 induced by relations []
            sage: dw = DyckWord([1,1,0,1,0,0,1,1,0,0])
            sage: TamariIntervalPosets.final_forest(dw)
            The Tamari interval of size 5 induced by relations [(5, 4), (3, 1), (2, 1)]

        TESTS::

            sage: TamariIntervalPosets.final_forest('mont')
            Traceback (most recent call last):
            ...
            ValueError: do not know how to construct the final forest of mont
        """
        if isinstance(element, TamariIntervalPoset):
            return element.final_forest()
        elif element in DyckWords():
            binary_tree = element.to_binary_tree_tamari()
        elif element in BinaryTrees() or element in LabelledBinaryTrees():
            binary_tree = element
        else:
            raise ValueError("do not know how to construct the final forest of {}".format(element))

        def get_relations(bt, start=1):
            r"""
            Recursive method to get the binary tree final forest relations
            with only one recursive reading of the tree.

            The vertices are being labelled with integers starting with
            ``start``.

            OUTPUT:

            - the indices of the nodes on the left border of the tree
              (these become the roots of the forest)
            - the relations of the final forest (as a list of tuples)
            - the next available index for a node (size of tree +
              ``start``)
            """
            if not bt:
                return [], [], start  # leaf
            roots, relations, index = get_relations(bt[0], start=start)
            rroots, rrelations, rindex = get_relations(bt[1], start=index + 1)
            roots.append(index)
            relations.extend(rrelations)
            relations.extend([(j, index) for j in rroots])
            return roots, relations, rindex

        roots, relations, index = get_relations(binary_tree)
        P = FinitePoset(DiGraph([list(range(1, index)), relations],
                                format='vertices_and_edges'))  # type:ignore
        return TamariIntervalPoset(P, check=False)  # type:ignore

    @staticmethod
    def initial_forest(element) -> TIP:
        r"""
        Return the initial forest of a binary tree, an interval-poset or
        a Dyck word.

        An initial forest is an interval-poset corresponding to an initial
        interval of the Tamari lattice, i.e., containing only increasing
        relations.

        It can be constructed from a binary tree by its binary
        search tree labeling with the rule: `a` precedes `b` in the
        initial forest iff `a` is in the left subtree of `b` in the
        binary search tree.

        INPUT:

        - ``element`` -- a binary tree, a Dyck word or an interval-poset

        EXAMPLES::

            sage: ip = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)])
            sage: TamariIntervalPosets.initial_forest(ip)
            The Tamari interval of size 4 induced by relations [(1, 2), (2, 3)]

        with binary trees::

            sage: bt = BinaryTree(); bt
            .
            sage: TamariIntervalPosets.initial_forest(bt)
            The Tamari interval of size 0 induced by relations []
            sage: bt = BinaryTree([]); bt
            [., .]
            sage: TamariIntervalPosets.initial_forest(bt)
            The Tamari interval of size 1 induced by relations []
            sage: bt = BinaryTree([[],None]); bt
            [[., .], .]
            sage: TamariIntervalPosets.initial_forest(bt)
            The Tamari interval of size 2 induced by relations [(1, 2)]
            sage: bt = BinaryTree([None,[]]); bt
            [., [., .]]
            sage: TamariIntervalPosets.initial_forest(bt)
            The Tamari interval of size 2 induced by relations []
            sage: bt = BinaryTree([[],[]]); bt
            [[., .], [., .]]
            sage: TamariIntervalPosets.initial_forest(bt)
            The Tamari interval of size 3 induced by relations [(1, 2)]
            sage: bt = BinaryTree([[None,[[],None]],[]]); bt
            [[., [[., .], .]], [., .]]
            sage: TamariIntervalPosets.initial_forest(bt)
            The Tamari interval of size 5 induced by relations [(1, 4), (2, 3), (3, 4)]

        from Dyck words::

            sage: dw = DyckWord([1,0])
            sage: TamariIntervalPosets.initial_forest(dw)
            The Tamari interval of size 1 induced by relations []
            sage: dw = DyckWord([1,1,0,1,0,0,1,1,0,0])
            sage: TamariIntervalPosets.initial_forest(dw)
            The Tamari interval of size 5 induced by relations [(1, 4), (2, 3), (3, 4)]

        TESTS::

            sage: TamariIntervalPosets.initial_forest('mont')
            Traceback (most recent call last):
            ...
            ValueError: do not know how to construct the initial forest of mont

        """
        if isinstance(element, TamariIntervalPoset):
            return element.initial_forest()
        elif element in DyckWords():
            binary_tree = element.to_binary_tree_tamari()
        elif element in BinaryTrees() or element in LabelledBinaryTrees():
            binary_tree = element
        else:
            raise ValueError("do not know how to construct the initial forest of {}".format(element))

        def get_relations(bt, start=1):
            r"""
            Recursive method to get the binary tree initial forest
            relations with only one recursive reading of the tree.

            The vertices are being labelled with integers starting with
            ``start``.

            OUTPUT:

            - the indices of the nodes on the right border of the tree
              (these become the roots of the forest)
            - the relations of the initial forest (as a list of tuples)
            - the next available index for a node (size of tree +
              ``start``)
            """
            if not bt:
                return [], [], start  # leaf
            lroots, lrelations, index = get_relations(bt[0], start=start)
            roots, relations, rindex = get_relations(bt[1], start=index + 1)
            roots.append(index)
            relations.extend(lrelations)
            relations.extend([(j, index) for j in lroots])
            return roots, relations, rindex

        roots, relations, index = get_relations(binary_tree)
        P = FinitePoset(DiGraph([list(range(1, index)), relations],
                                format='vertices_and_edges'))  # type:ignore
        return TamariIntervalPoset(P, check=False)  # type:ignore

    @staticmethod
    def from_binary_trees(tree1, tree2) -> TIP:
        r"""
        Return the interval-poset corresponding to the interval
        [``tree1``, ``tree2``] of the Tamari lattice.

        Raise an exception if
        ``tree1`` is not `\leq` ``tree2`` in the Tamari lattice.

        INPUT:

        - ``tree1`` -- a binary tree
        - ``tree2`` -- a binary tree greater or equal than ``tree1`` for
          the Tamari lattice

        EXAMPLES::

            sage: tree1 = BinaryTree([[],None])
            sage: tree2 = BinaryTree([None,[]])
            sage: TamariIntervalPosets.from_binary_trees(tree1,tree2)
            The Tamari interval of size 2 induced by relations []
            sage: TamariIntervalPosets.from_binary_trees(tree1,tree1)
            The Tamari interval of size 2 induced by relations [(1, 2)]
            sage: TamariIntervalPosets.from_binary_trees(tree2,tree2)
            The Tamari interval of size 2 induced by relations [(2, 1)]

            sage: tree1 = BinaryTree([[],[[None,[]],[]]])
            sage: tree2 = BinaryTree([None,[None,[None,[[],[]]]]])
            sage: TamariIntervalPosets.from_binary_trees(tree1,tree2)
            The Tamari interval of size 6 induced by relations [(4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]

            sage: tree3 = BinaryTree([None,[None,[[],[None,[]]]]])
            sage: TamariIntervalPosets.from_binary_trees(tree1,tree3)
            Traceback (most recent call last):
            ...
            ValueError: the two binary trees are not comparable on the Tamari lattice
            sage: TamariIntervalPosets.from_binary_trees(tree1,BinaryTree())
            Traceback (most recent call last):
            ...
            ValueError: the two binary trees are not comparable on the Tamari lattice
        """
        initial_forest = TamariIntervalPosets.initial_forest(tree2)
        final_forest = TamariIntervalPosets.final_forest(tree1)
        try:
            return initial_forest.intersection(final_forest)
        except Exception:
            raise ValueError("the two binary trees are not comparable on the Tamari lattice")

    @staticmethod
    def from_dyck_words(dw1, dw2) -> TIP:
        r"""
        Return the interval-poset corresponding to the interval
        [``dw1``, ``dw2``] of the Tamari lattice.

        Raise an exception if the
        two Dyck words ``dw1`` and ``dw2`` do not satisfy
        ``dw1`` `\leq` ``dw2`` in the Tamari lattice.

        INPUT:

        - ``dw1`` -- a Dyck word
        - ``dw2`` -- a Dyck word greater or equal than ``dw1`` for
          the Tamari lattice

        EXAMPLES::

            sage: dw1 = DyckWord([1,0,1,0])
            sage: dw2 = DyckWord([1,1,0,0])
            sage: TamariIntervalPosets.from_dyck_words(dw1,dw2)
            The Tamari interval of size 2 induced by relations []
            sage: TamariIntervalPosets.from_dyck_words(dw1,dw1)
            The Tamari interval of size 2 induced by relations [(1, 2)]
            sage: TamariIntervalPosets.from_dyck_words(dw2,dw2)
            The Tamari interval of size 2 induced by relations [(2, 1)]

            sage: dw1 = DyckWord([1,0,1,1,1,0,0,1,1,0,0,0])
            sage: dw2 = DyckWord([1,1,1,1,0,1,1,0,0,0,0,0])
            sage: TamariIntervalPosets.from_dyck_words(dw1,dw2)
            The Tamari interval of size 6 induced by relations [(4, 5), (6, 5), (5, 2), (4, 3), (3, 2)]

            sage: dw3 = DyckWord([1,1,1,0,1,1,1,0,0,0,0,0])
            sage: TamariIntervalPosets.from_dyck_words(dw1,dw3)
            Traceback (most recent call last):
            ...
            ValueError: the two Dyck words are not comparable on the Tamari lattice
            sage: TamariIntervalPosets.from_dyck_words(dw1,DyckWord([1,0]))
            Traceback (most recent call last):
            ...
            ValueError: the two Dyck words are not comparable on the Tamari lattice
        """
        tree1 = dw1.to_binary_tree_tamari()
        tree2 = dw2.to_binary_tree_tamari()
        try:
            return TamariIntervalPosets.from_binary_trees(tree1, tree2)
        except Exception:
            raise ValueError("the two Dyck words are not comparable on the Tamari lattice")

    @staticmethod
    def recomposition_from_triple(left: TIP, right: TIP, r) -> TIP:
        """
        Recompose an interval-poset from a triple (``left``, ``right``, ``r``).

        For the inverse method,
        see :meth:`TamariIntervalPoset.decomposition_to_triple`.

        INPUT:

        - ``left`` -- an interval-poset
        - ``right`` -- an interval-poset
        - ``r`` -- the parameter of the decomposition, an integer

        OUTPUT: an interval-poset

        EXAMPLES::

            sage: T1 = TamariIntervalPoset(3, [(1, 2), (3, 2)])
            sage: T2 = TamariIntervalPoset(4, [(2, 3), (4, 3)])
            sage: TamariIntervalPosets.recomposition_from_triple(T1, T2, 2)
            The Tamari interval of size 8 induced by relations [(1, 2), (2, 4),
            (3, 4), (6, 7), (8, 7), (6, 4), (5, 4), (3, 2)]

        REFERENCES:

        - [Pons2018]_
        """
        root = left.size() + 1
        rel = left.poset().cover_relations()
        rel.extend((i, root) for i in left)
        rel.extend((root + a, root + b)
                   for a, b in right.poset().cover_relations())
        decroot = right.decreasing_roots()[:r]
        rel.extend((root + i, root) for i in decroot)
        # does this describe only cover relations ?
        # if yes, create the poset here and pass it as argument
        return TamariIntervalPoset(left.size() + right.size() + 1, rel)

    @staticmethod
    def from_grafting_tree(tree) -> TIP:
        """
        Return an interval-poset from a grafting tree.

        For the inverse method,
        see :meth:`TamariIntervalPoset.grafting_tree`.

        EXAMPLES::

            sage: tip = TamariIntervalPoset(8, [(1,2), (2,4), (3,4), (6,7), (3,2), (5,4), (6,4), (8,7)])
            sage: t = tip.grafting_tree()
            sage: TamariIntervalPosets.from_grafting_tree(t) == tip
            True

        REFERENCES:

        - [Pons2018]_
        """
        if tree.is_empty():
            return TamariIntervalPoset(0, [])
        r = tree.label()
        left = TamariIntervalPosets.from_grafting_tree(tree[0])
        right = TamariIntervalPosets.from_grafting_tree(tree[1])
        return TamariIntervalPosets.recomposition_from_triple(left, right, r)

    @staticmethod
    def from_minimal_schnyder_wood(graph) -> TIP:
        """
        Return a Tamari interval built from a minimal Schnyder wood.

        This is an implementation of Bernardi and Bonichon's bijection
        [BeBo2009]_.

        INPUT:

        a minimal Schnyder wood, given as a graph with colored and
        oriented edges, without the three exterior unoriented edges

        The three boundary vertices must be -1, -2 and -3.

        One assumes moreover that the embedding around -1 is the
        list of neighbors of -1 and not just a cyclic permutation of that.

        Beware that the embedding convention used here is the opposite of
        the one used by the plot method.

        OUTPUT:

        a Tamari interval-poset

        EXAMPLES:

        A small example::

            sage: TIP = TamariIntervalPosets
            sage: G = DiGraph([(0,-1,0),(0,-2,1),(0,-3,2)], format='list_of_edges')
            sage: G.set_embedding({-1:[0],-2:[0],-3:[0],0:[-1,-2,-3]})
            sage: TIP.from_minimal_schnyder_wood(G)
            The Tamari interval of size 1 induced by relations []

        An example from page 14 of [BeBo2009]_::

            sage: c0 = [(0,-1),(1,0),(2,0),(4,3),(3,-1),(5,3)]
            sage: c1 = [(5,-2),(3,-2),(4,5),(1,3),(2,3),(0,3)]
            sage: c2 = [(0,-3),(1,-3),(3,-3),(4,-3),(5,-3),(2,1)]
            sage: ed = [(u,v,0) for u,v in c0]
            sage: ed += [(u,v,1) for u,v in c1]
            sage: ed += [(u,v,2) for u,v in c2]
            sage: G = DiGraph(ed, format='list_of_edges')
            sage: embed = {-1:[3,0],-2:[5,3],-3:[0,1,3,4,5]}
            sage: data_emb = [[3,2,1,-3,-1],[2,3,-3,0],[3,1,0]]
            sage: data_emb += [[-2,5,4,-3,1,2,0,-1],[5,-3,3],[-2,-3,4,3]]
            sage: for k in range(6):
            ....:     embed[k] = data_emb[k]
            sage: G.set_embedding(embed)
            sage: TIP.from_minimal_schnyder_wood(G)
            The Tamari interval of size 6 induced by relations [(1, 4), (2, 4), (3, 4), (5, 6), (6, 4), (5, 4), (3, 1), (2, 1)]

        An example from page 18 of [BeBo2009]_::

            sage: c0 = [(0,-1),(1,0),(2,-1),(3,2),(4,2),(5,-1)]
            sage: c1 = [(5,-2),(2,-2),(4,-2),(3,4),(1,2),(0,2)]
            sage: c2 = [(0,-3),(1,-3),(3,-3),(4,-3),(2,-3),(5,2)]
            sage: ed = [(u,v,0) for u,v in c0]
            sage: ed += [(u,v,1) for u,v in c1]
            sage: ed += [(u,v,2) for u,v in c2]
            sage: G = DiGraph(ed, format='list_of_edges')
            sage: embed = {-1:[5,2,0],-2:[4,2,5],-3:[0,1,2,3,4]}
            sage: data_emb = [[2,1,-3,-1],[2,-3,0],[3,-3,1,0,-1,5,-2,4]]
            sage: data_emb += [[4,-3,2],[-2,-3,3,2],[-2,2,-1]]
            sage: for k in range(6):
            ....:     embed[k] = data_emb[k]
            sage: G.set_embedding(embed)
            sage: TIP.from_minimal_schnyder_wood(G)
            The Tamari interval of size 6 induced by relations [(1, 3), (2, 3), (4, 5), (5, 3), (4, 3), (2, 1)]

        Another small example::

            sage: c0 = [(0,-1),(2,-1),(1,0)]
            sage: c1 = [(2,-2),(1,-2),(0,2)]
            sage: c2 = [(0,-3),(1,-3),(2,1)]
            sage: ed = [(u,v,0) for u,v in c0]
            sage: ed += [(u,v,1) for u,v in c1]
            sage: ed += [(u,v,2) for u,v in c2]
            sage: G = DiGraph(ed, format='list_of_edges')
            sage: embed = {-1:[2,0],-2:[1,2],-3:[0,1]}
            sage: data_emb = [[2,1,-3,-1],[-3,0,2,-2],[-2,1,0,-1]]
            sage: for k in range(3):
            ....:     embed[k] = data_emb[k]
            sage: G.set_embedding(embed)
            sage: TIP.from_minimal_schnyder_wood(G)
            The Tamari interval of size 3 induced by relations [(2, 3), (2, 1)]
        """
        from sage.graphs.digraph import DiGraph
        from sage.combinat.dyck_word import DyckWord
        color_a = graph.incoming_edges(-1)[0][2]
        color_b = graph.incoming_edges(-2)[0][2]

        embedding = graph.get_embedding()
        graph0 = DiGraph([e for e in graph.edges(sort=False)
                          if e[2] == color_a],
                         format='list_of_edges')
        restricted_embedding = {u: [v for v in embedding[u]
                                    if v in graph0.neighbors_in(u) or
                                    v in graph0.neighbors_out(u)]
                                for u in graph0}

        voisins_in = {}
        for u in graph0:
            if u != -1:
                bad_emb = restricted_embedding[u]
                sortie = graph0.neighbors_out(u)[0]
                idx = bad_emb.index(sortie)
                restricted_embedding[u] = bad_emb[idx:] + bad_emb[:idx]
                voisins_in[u] = restricted_embedding[u][1:]
            else:
                voisins_in[u] = list(restricted_embedding[u])
            voisins_in[u].reverse()  # To have them in the right order

        graph0.set_embedding(restricted_embedding)

        def clockwise_labelling(gr, vertex):
            if len(gr) == 1:
                return [vertex]
            else:
                lbl = [vertex]
                for w in voisins_in[vertex]:
                    lbl += clockwise_labelling(gr, w)
                return lbl

        def profil(gr, vertex):
            if len(gr) == 1:
                return []
            else:
                lbl = []
                for w in voisins_in[vertex]:
                    lbl += [1] + profil(gr, w) + [0]
                return lbl

        dyckword_bottom = profil(graph0, -1)
        # this is the profile of the planar graph graph0

        liste = clockwise_labelling(graph0, -1)[1:]
        relabelling = {l: i for i, l in enumerate(liste)}
        for l in [-1, -2, -3]:
            relabelling[l] = l
        new_graph = graph.relabel(relabelling, inplace=False)

        dyckword_top = []
        for i in range(1, len(graph) - 3):
            indegree1 = len([u for u in new_graph.incoming_edges(i)
                             if u[2] == color_b])
            dyckword_top += [1] + [0] * indegree1
        indegree1 = len([u for u in new_graph.incoming_edges(-2)
                         if u[2] == color_b])
        dyckword_top += [1] + [0] * indegree1

        dyckword_bottom = DyckWord(dyckword_bottom)  # type:ignore
        dyckword_top = DyckWord(dyckword_top)  # type:ignore
        TIP = TamariIntervalPosets(len(dyckword_bottom) // 2)
        return TIP.from_dyck_words(dyckword_bottom, dyckword_top)

    def __call__(self, *args, **keywords):
        r"""
        Allows for a poset to be directly transformed into an interval-poset.

        It is some kind of coercion but cannot be made through the coercion
        system because posets do not have parents.

        EXAMPLES::

            sage: TIP = TamariIntervalPosets()
            sage: p = Poset( ([1,2,3], [(1,2)]))
            sage: TIP(p)
            The Tamari interval of size 3 induced by relations [(1, 2)]
            sage: TIP(TIP(p))
            The Tamari interval of size 3 induced by relations [(1, 2)]
            sage: TIP(3,[(1,2)])
            The Tamari interval of size 3 induced by relations [(1, 2)]
            sage: p = Poset(([1,2,3],[(1,3)]))
            sage: TIP(p)
            Traceback (most recent call last):
            ...
            ValueError: this does not satisfy the Tamari interval-poset condition
        """
        if isinstance(args[0], TamariIntervalPoset):
            return args[0]
        if len(args) == 1 and isinstance(args[0], FinitePoset):
            return self.element_class(self, args[0])

        return super(TamariIntervalPosets, self).__call__(*args, **keywords)

    def le(self, el1, el2) -> bool:
        r"""
        Poset structure on the set of interval-posets.

        The comparison is first by size, then using the
        cubical coordinates.

        .. SEEALSO:: :meth:`cubical_coordinates`

        INPUT:

        - ``el1`` -- an interval-poset
        - ``el2`` -- an interval-poset

        EXAMPLES::

            sage: ip1 = TamariIntervalPoset(4,[(1,2),(2,3),(4,3)])
            sage: ip2 = TamariIntervalPoset(4,[(1,2),(2,3)])
            sage: TamariIntervalPosets().le(ip1,ip2)
            False
            sage: TamariIntervalPosets().le(ip2,ip1)
            True
        """
        cc1 = el1.cubical_coordinates()
        cc2 = el2.cubical_coordinates()
        return all(x1 <= x2 for x1, x2 in zip(cc1, cc2))

#################################################################
# Enumerated set of all Tamari Interval-posets
#################################################################


class TamariIntervalPosets_all(DisjointUnionEnumeratedSets, TamariIntervalPosets):
    r"""
    The enumerated set of all Tamari interval-posets.
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.combinat.interval_posets import TamariIntervalPosets_all
            sage: S = TamariIntervalPosets_all()
            sage: S.cardinality()
            +Infinity

            sage: it = iter(S)
            sage: [next(it) for i in range(5)]
            [The Tamari interval of size 0 induced by relations [],
             The Tamari interval of size 1 induced by relations [],
             The Tamari interval of size 2 induced by relations [],
             The Tamari interval of size 2 induced by relations [(2, 1)],
             The Tamari interval of size 2 induced by relations [(1, 2)]]
            sage: next(it).parent()
            Interval-posets
            sage: S(0,[])
            The Tamari interval of size 0 induced by relations []

            sage: S is TamariIntervalPosets_all()
            True
            sage: TestSuite(S).run()  # long time (7s)
            """
        DisjointUnionEnumeratedSets.__init__(
            self, Family(NonNegativeIntegers(), TamariIntervalPosets_size),
            facade=True, keepkey=False,
            category=(Posets(), EnumeratedSets(), Monoids()))

    def _repr_(self) -> str:
        r"""
        TESTS::

            sage: TamariIntervalPosets()
            Interval-posets
        """
        return "Interval-posets"

    def one(self) -> TIP:
        """
        Return the unit of the monoid.

        This is the empty interval poset, of size 0.

        EXAMPLES::

            sage: TamariIntervalPosets().one()
            The Tamari interval of size 0 induced by relations []
        """
        return TamariIntervalPoset(0, [])

    def _element_constructor_(self, size, relations) -> TIP:
        r"""
        EXAMPLES::

            sage: TIP = TamariIntervalPosets()
            sage: TIP(3,[(1,2)])
            The Tamari interval of size 3 induced by relations [(1, 2)]
        """
        return self.element_class(self, size, relations)

    def __contains__(self, x) -> bool:
        r"""
        TESTS::

            sage: S = TamariIntervalPosets()
            sage: 1 in S
            False
            sage: S(0,[]) in S
            True
        """
        return isinstance(x, self.element_class)

    Element = TamariIntervalPoset


TIP = TamariIntervalPoset


#################################################################
# Enumerated set of Tamari interval-posets of a given size
#################################################################
class TamariIntervalPosets_size(TamariIntervalPosets):
    r"""
    The enumerated set of interval-posets of a given size.
    """
    def __init__(self, size):
        r"""
        TESTS::

            sage: S = TamariIntervalPosets(3)
            sage: assert S is TamariIntervalPosets(3)
            sage: for i in range(5): TestSuite(TamariIntervalPosets(i)).run()
        """
        # there is a natural order on interval-posets through inclusions
        # that is why we use the FinitePosets category
        super(TamariIntervalPosets_size, self).__init__(category=(FinitePosets(), FiniteEnumeratedSets()))

        self._size = size

    def _repr_(self) -> str:
        r"""
        TESTS::

            sage: TamariIntervalPosets(3)
            Interval-posets of size 3
        """
        return "Interval-posets of size {}".format(self._size)

    def __contains__(self, x) -> bool:
        r"""
        TESTS::

            sage: S = TamariIntervalPosets(3)
            sage: 1 in S
            False
            sage: S([]) in S
            True
        """
        return isinstance(x, self.element_class) and x.size() == self._size

    def cardinality(self) -> Integer:
        r"""
        The cardinality of ``self``. That is, the number of
        interval-posets of size `n`.

        The formula was given in [Cha2008]_:

        .. MATH::

            \frac{2(4n+1)!}{(n+1)!(3n+2)!}
            = \frac{2}{n(n+1)} \binom{4n+1}{n-1}.

        EXAMPLES::

            sage: [TamariIntervalPosets(i).cardinality() for i in range(6)]
            [1, 1, 3, 13, 68, 399]
        """
        from sage.arith.all import binomial
        n = self._size
        if n == 0:
            return Integer(1)
        return (2 * binomial(4 * n + 1, n - 1)) // (n * (n + 1))
        # return Integer(2 * factorial(4*n+1)/(factorial(n+1)*factorial(3*n+2)))

    def __iter__(self) -> Iterator[TIP]:
        r"""
        Recursive generation: we iterate through all interval-posets of
        size ``size - 1`` and add all possible relations to the last
        vertex.

        This works thanks to the fact that the restriction of an
        interval-poset of size `n` to the subset `\{1, 2, \ldots, k\}` for
        a fixed `k \leq n` is an interval-poset.

        TESTS::

            sage: TIP1 = TamariIntervalPosets(1)
            sage: list(TIP1)
            [The Tamari interval of size 1 induced by relations []]
            sage: TIP2 = TamariIntervalPosets(2)
            sage: list(TIP2)
            [The Tamari interval of size 2 induced by relations [],
             The Tamari interval of size 2 induced by relations [(2, 1)],
             The Tamari interval of size 2 induced by relations [(1, 2)]]
            sage: TIP3 = TamariIntervalPosets(3)
            sage: list(TIP3)
            [The Tamari interval of size 3 induced by relations [],
             The Tamari interval of size 3 induced by relations [(3, 2)],
             The Tamari interval of size 3 induced by relations [(2, 3)],
             The Tamari interval of size 3 induced by relations [(1, 3), (2, 3)],
             The Tamari interval of size 3 induced by relations [(2, 1)],
             The Tamari interval of size 3 induced by relations [(3, 2), (2, 1)],
             The Tamari interval of size 3 induced by relations [(3, 1), (2, 1)],
             The Tamari interval of size 3 induced by relations [(2, 3), (2, 1)],
             The Tamari interval of size 3 induced by relations [(2, 3), (3, 1), (2, 1)],
             The Tamari interval of size 3 induced by relations [(1, 3), (2, 3), (2, 1)],
             The Tamari interval of size 3 induced by relations [(1, 2)],
             The Tamari interval of size 3 induced by relations [(1, 2), (3, 2)],
             The Tamari interval of size 3 induced by relations [(1, 2), (2, 3)]]
            sage: all(len(list(TamariIntervalPosets(i)))==TamariIntervalPosets(i).cardinality() for i in range(6))
            True
        """
        n = self._size
        if n <= 1:
            yield TamariIntervalPoset(n, [], check=False)
            return

        for tip in TamariIntervalPosets(n - 1):
            new_tip = TamariIntervalPoset(n, tip._cover_relations, check=False)
            yield new_tip  # we have added an extra vertex but no relations

            # adding a decreasing relation n>>m2 with m2<n and no
            # increasing relations
            for m2 in range(n - 1, 0, -1):
                if new_tip.le(n - 1, m2):
                    yield TamariIntervalPoset(n, new_tip._cover_relations + ((n, m2),), check=False)

            for m in range(n - 1, 0, -1):
                # adding an increasing relation m>>n
                if not new_tip.le(m, n):
                    new_tip = TamariIntervalPoset(n, new_tip._cover_relations + ((m, n),), check=False)
                    yield new_tip
                else:
                    continue

                # further adding a decreasing relation n>>m2 with m2<m
                for m2 in range(m - 1, 0, -1):
                    if new_tip.le(n - 1, m2):
                        yield TamariIntervalPoset(n, new_tip._cover_relations + ((n, m2),), check=False)

    def random_element(self) -> TIP:
        """
        Return a random Tamari interval of fixed size.

        This is obtained by first creating a random rooted
        planar triangulation, then computing its unique
        minimal Schnyder wood, then applying a bijection
        of Bernardi and Bonichon [BeBo2009]_.

        Because the random rooted planar triangulation is
        chosen uniformly at random, the Tamari interval is
        also chosen according to the uniform distribution.

        EXAMPLES::

            sage: T = TamariIntervalPosets(4).random_element()
            sage: T.parent()
            Interval-posets
            sage: u = T.lower_dyck_word(); u   # random
            [1, 1, 0, 1, 0, 0, 1, 0]
            sage: v = T.lower_dyck_word(); v   # random
            [1, 1, 0, 1, 0, 0, 1, 0]
            sage: len(u)
            8
        """
        from sage.graphs.schnyder import minimal_schnyder_wood
        from sage.graphs.generators.random import RandomTriangulation
        n = self._size
        tri = RandomTriangulation(n + 3)
        TIP = TamariIntervalPosets
        schnyder = minimal_schnyder_wood(tri, root_edge=(-1, -2),
                                         check=False)
        return TIP.from_minimal_schnyder_wood(schnyder)

    @lazy_attribute
    def _parent_for(self):
        r"""
        The parent of the element generated by ``self``.

        TESTS::

            sage: TIP3 = TamariIntervalPosets(3)
            sage: TIP3._parent_for
            Interval-posets
        """
        return TamariIntervalPosets_all()

    # This is needed because this is a facade parent via DisjointUnionEnumeratedSets
    @lazy_attribute
    def element_class(self):
        r"""
        TESTS::

            sage: S = TamariIntervalPosets(3)
            sage: S.element_class
            <class 'sage.combinat.interval_posets.TamariIntervalPosets_all_with_category.element_class'>
            sage: S.first().__class__ == TamariIntervalPosets().first().__class__
            True
        """
        return self._parent_for.element_class

    def _element_constructor_(self, relations) -> TIP:
        r"""
        EXAMPLES::

            sage: TIP3 = TamariIntervalPosets(3)
            sage: TIP3([(1,2)])
            The Tamari interval of size 3 induced by relations [(1, 2)]
            sage: TIP3([(3,4)])
            Traceback (most recent call last):
            ...
            ValueError: the relations do not correspond to the size of the poset
        """
        return self.element_class(self, self._size, relations)
