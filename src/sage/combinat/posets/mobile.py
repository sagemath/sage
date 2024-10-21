# sage_setup: distribution = sagemath-graphs
"""
Mobile posets
"""
# ****************************************************************************
#       Copyright (C) 2020 Stefan Grosser <stefan.grosser1@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.combinat.posets.posets import Poset, FinitePoset
from sage.misc.lazy_attribute import lazy_attribute
from .linear_extensions import LinearExtensionsOfMobile


class MobilePoset(FinitePoset):
    r"""
    A mobile poset.

    Mobile posets are an extension of d-complete posets which permit a determinant
    formula for counting linear extensions. They are formed by having a ribbon
    poset with d-complete posets 'hanging' below it and at most one
    d-complete poset above it, known as the anchor. See [GGMM2020]_
    for the definition.

    EXAMPLES::

        sage: P = posets.MobilePoset(posets.RibbonPoset(7, [1,3]),                      # needs sage.combinat sage.modules
        ....:                        {1: [posets.YoungDiagramPoset([3, 2], dual=True)],
        ....:                         3: [posets.DoubleTailedDiamond(6)]},
        ....:                        anchor=(4, 2, posets.ChainPoset(6)))
        sage: len(P._ribbon)                                                            # needs sage.combinat sage.modules
        8
        sage: P._anchor                                                                 # needs sage.combinat sage.modules
        (4, 5)

    This example is Example 5.9 in [GGMM2020]_::

        sage: P1 = posets.MobilePoset(posets.RibbonPoset(8, [2,3,4]),
        ....:                         {4: [posets.ChainPoset(1)]},
        ....:                         anchor=(3, 0, posets.ChainPoset(1)))
        sage: sorted([P1._element_to_vertex(i) for i in P1._ribbon])
        [0, 1, 2, 6, 7, 9]
        sage: P1._anchor
        (3, 2)

        sage: P2 = posets.MobilePoset(posets.RibbonPoset(15, [1,3,5,7,9,11,13]),
        ....:                         {}, anchor=(8, 0, posets.ChainPoset(1)))
        sage: sorted(P2._ribbon)
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]
        sage: P2._anchor
        (8, (8, 0))
        sage: P2.linear_extensions().cardinality()                                      # needs sage.modules
        21399440939

        sage: EP = posets.MobilePoset(posets.ChainPoset(0), {})
        Traceback (most recent call last):
        ...
        ValueError: the empty poset is not a mobile poset
    """
    _lin_ext_type = LinearExtensionsOfMobile
    _desc = 'Finite mobile poset'

    def __init__(self, hasse_diagram, elements, category, facade, key, ribbon=None, check=True) -> None:
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: P = posets.MobilePoset(posets.RibbonPoset(15, [1,3,5,7,9,11,13]),
            ....:                        {}, anchor=(8, 0, posets.ChainPoset(1)))
            sage: TestSuite(P).run()
        """
        FinitePoset.__init__(self, hasse_diagram=hasse_diagram, elements=elements,
                             category=category, facade=facade, key=key)
        if not self._hasse_diagram:
            raise ValueError("the empty poset is not a mobile poset")

        if ribbon:
            if check and not self._is_valid_ribbon(ribbon):
                raise ValueError("invalid ribbon")
            self._ribbon = ribbon

    def _is_valid_ribbon(self, ribbon):
        r"""
        Return ``True`` if a ribbon has at most one anchor, no vertex has two
        or more anchors, and every hanging poset is d-complete.

        INPUT:

        - ``ribbon`` -- list of elements that form a ribbon in your poset

        TESTS::

            sage: P = posets.RibbonPoset(5, [2])
            sage: P._is_valid_ribbon([0,1,2,3,4])
            True
            sage: P._is_valid_ribbon([2])
            False
            sage: P._is_valid_ribbon([2,3,4])
            True
            sage: P._is_valid_ribbon([2,3])
            True
        """
        ribbon = [self._element_to_vertex(x) for x in ribbon]
        G = self._hasse_diagram
        G_un = G.to_undirected().copy(immutable=False)
        R = G.subgraph(ribbon)
        num_anchors = 0

        for r in ribbon:
            anchor_neighbors = set(G.neighbor_out_iterator(r)).difference(set(R.neighbor_out_iterator(r)))
            if len(anchor_neighbors) == 1:
                num_anchors += 1
            elif len(anchor_neighbors) > 1:
                return False

            for lc in G.neighbor_in_iterator(r):
                if lc in ribbon:
                    continue

                G_un.delete_edge(lc, r)
                P = Poset(G.subgraph(G_un.connected_component_containing_vertex(lc, sort=False)))
                if P.top() != lc or not P.is_d_complete():
                    return False
                G_un.add_edge(lc, r)

        return True

    @lazy_attribute
    def _anchor(self):
        r"""
        The anchor of the mobile poset.

        TESTS::

            sage: from sage.combinat.posets.mobile import MobilePoset
            sage: M = MobilePoset(DiGraph([[0,1,2,3,4,5,6,7,8],
            ....:   [(1,0),(3,0),(2,1),(2,3),(4,3), (5,4),(5,6),(7,4),(7,8)]]))
            sage: M._anchor
            (4, 3)
        """
        ribbon = [self._element_to_vertex(x) for x in self._ribbon]
        H = self._hasse_diagram
        R = H.subgraph(ribbon)

        anchor = None

        # Find the anchor vertex, if it exists, and return the edge
        for r in ribbon:
            anchor_neighbors = set(H.neighbor_out_iterator(r)).difference(set(R.neighbor_out_iterator(r)))
            if len(anchor_neighbors) == 1:
                anchor = (r, anchor_neighbors.pop())
                break
        return (self._vertex_to_element(anchor[0]), self._vertex_to_element(anchor[1])) if anchor is not None else None

    @lazy_attribute
    def _ribbon(self):
        r"""
        The ribbon of the mobile poset.

        TESTS::

            sage: from sage.combinat.posets.mobile import MobilePoset
            sage: M = MobilePoset(DiGraph([[0,1,2,3,4,5,6,7,8],
            ....:   [(1,0),(3,0),(2,1),(2,3),(4,3), (5,4),(5,6),(7,4),(7,8)]]))
            sage: sorted(M._ribbon)
            [4, 5, 6, 7, 8]
            sage: M._is_valid_ribbon(M._ribbon)
            True
            sage: M2 = MobilePoset(Poset([[0,1,2,3,4,5,6,7,8],
            ....:          [(1,0),(3,0),(2,1),(2,3),(4,3),(5,4),(7,4),(7,8)]]))
            sage: sorted(M2._ribbon)
            [4, 7, 8]
            sage: M2._is_valid_ribbon(M2._ribbon)
            True
        """
        H = self._hasse_diagram
        H_un = H.to_undirected()
        max_elmts = H.sinks()
        # Compute anchor, ribbon
        ribbon = []  # In order list of elements on zigzag

        if len(max_elmts) == 1:
            return [self._vertex_to_element(max_elmts[0])]
        # Compute max element tree by merging shortest paths
        start = max_elmts[0]

        zigzag_elmts = set()
        for m in max_elmts[1:]:
            sp = H_un.shortest_path(start, m)
            zigzag_elmts.update(sp)
        max_elmt_graph = H.subgraph(zigzag_elmts)
        G = max_elmt_graph.to_undirected()

        if G.is_path():
            # Check if there is a anchor by seeing if there is more than one acyclic path to the next max
            ends = max_elmt_graph.vertices(sort=True, degree=1)
            # Form ribbon
            ribbon = G.shortest_path(ends[0], ends[1])
            for end_count, end in enumerate(ends):
                if not (H_un.is_cut_vertex(end) or H_un.degree(end) == 1):
                    traverse_ribbon = ribbon if end_count == 0 else ribbon[::-1]
                    for ind, p in enumerate(traverse_ribbon):
                        if H_un.is_cut_edge(p, traverse_ribbon[ind + 1]):
                            return [self._vertex_to_element(r)
                                    for r in G.shortest_path(ends[(end_count + 1) % 2], traverse_ribbon[ind + 1])]
            return [self._vertex_to_element(r) for r in ribbon]

        # First check path counts between ends and deg3 vertex
        # Then check if more than one max elmt on way to degree 3 vertex.
        # Then check if the edge going to a max element is down from the degree 3 vertex
        # Arbitrarily choose between ones with just 1

        ends = max_elmt_graph.vertices(sort=True, degree=1)
        deg3 = max_elmt_graph.vertices(sort=True, degree=3)[0]

        anchoredEnd = None
        for end in ends:
            if not (H_un.is_cut_vertex(end) or H_un.degree(end) == 1):
                anchoredEnd = end
                break

        if anchoredEnd is not None:
            ends.remove(anchoredEnd)
            return [self._vertex_to_element(r) for r in G.shortest_path(ends[0], ends[1])]

        possible_anchors = ends[:]
        for end in ends:
            path = G.shortest_path(end, deg3)
            if sum(bool(z in max_elmts) for z in path) != 1:
                possible_anchors.remove(end)

        for p in possible_anchors:
            path = G.shortest_path(p, deg3)
            if max_elmt_graph.has_edge(path[-2], path[-1]):
                possible_anchors.remove(p)

        anchoredEnd = possible_anchors[0]
        ends.remove(anchoredEnd)
        return [self._vertex_to_element(r) for r in G.shortest_path(ends[0], ends[1])]

    def ribbon(self):
        r"""
        Return the ribbon of the mobile poset.

        EXAMPLES::

            sage: from sage.combinat.posets.mobile import MobilePoset
            sage: M3 = MobilePoset(Posets.RibbonPoset(5, [1,2]))
            sage: sorted(M3.ribbon())
            [1, 2, 3, 4]
        """
        return self._ribbon

    def anchor(self):
        r"""
        Return the anchor of the mobile poset.

        EXAMPLES::

            sage: from sage.combinat.posets.mobile import MobilePoset
            sage: M2 = MobilePoset(Poset([[0,1,2,3,4,5,6,7,8],
            ....:          [(1,0),(3,0),(2,1),(2,3),(4,3),(5,4),(7,4),(7,8)]]))
            sage: M2.anchor()
            (4, 3)
            sage: M3 = MobilePoset(Posets.RibbonPoset(5, [1,2]))
            sage: M3.anchor() is None
            True
        """
        return self._anchor
