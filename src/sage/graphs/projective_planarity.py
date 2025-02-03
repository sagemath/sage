r"""
Minimal forbidden minors for projective plane and function for finding them in G

Contain the 35 minimal forbidden minors for the projective plane,
a function for finding if a graph, G, has one of them as a minor.

EXAMPLES::

The Peterson graph is a known projective planar graph so it doesn't have a
forbidden minor::

    sage: P = graphs.PetersenGraph()
    sage: _ = get_p2_forbidden_minor(P); type(_)        # long
    <class 'NoneType'>

K_{4,4} has a projective plane crossing number of 2. One of the minimal
forbidden minors is K_{4,4} - e, so we get a 1-to-1 dictionary from
:meth:`~Graph.minor`::

    sage: K44 = graphs.CompleteBipartiteGraph(4, 4)
    sage: minor_map = get_p2_forbidden_minor(K44); minor_map
    {0: [0], 1: [1], 2: [2], 3: [4], 4: [8], 5: [6], 6: [5], 7: [7]}

AUTHORS:

- Juan M. Lazaro Ruiz (2025-01-27): initial version

"""

# ****************************************************************************
#       Copyright (C) 2024 Juan M. Lazaro Ruiz <juan.m.lazaro.ruiz@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from copy import deepcopy

from sage.graphs.graph import Graph
from sage.graphs.generators.basic import (
    CompleteGraph,
    CompleteBipartiteGraph,
    CompleteMultipartiteGraph,
    CycleGraph
)


# Kuratowski graphs - don't use any functions/methods that operate in place
K33 = CompleteBipartiteGraph(3, 3)
K5 = CompleteGraph(5)

# Common verticies to merge
VERTEX_PAIR_1 = [(0, 0), (1, 0)]
VERTEX_PAIR_2 = [(0, 1), (1, 1)]
FIRST_EDGE = [(0, 0), (0, 1)]

# The 35 minimal forbidden minors for the projective plane as listed in
# Graphs on Surfaces by Mohar and Thomassen, IBSN: 0-8018-6689-8
# Some are also named according to:
# https://mathworld.wolfram.com/ProjectivePlanarGraph.html#:~:text=A%20graph%20with%20projective%20plane%20crossing%20number%20equal%20to%200%20may%20be%20said%20to%20be%20projective%20planar

G1 = K33.disjoint_union(K33)
G2 = K5.disjoint_union(K33)
G3 = K5.disjoint_union(K5)

# K_{3,3} \one_vertex_identification K_{3,3}
G4 = K33.disjoint_union(K33)
G4.merge_vertices(VERTEX_PAIR_1)

# K_5 \one_vertex_identification K_{3,3}
G5 = K5.disjoint_union(K33)
G5.merge_vertices(VERTEX_PAIR_1)

# K_5 \one_vertex_identification K_5
G6 = K5.disjoint_union(K5)
G6.merge_vertices(VERTEX_PAIR_1)

# \mathcal{B}_3 = K_5 \two_vertex_identification_minus_shared_edge K_5
G7 = K5.disjoint_union(K5)
G7.merge_vertices(VERTEX_PAIR_1)
G7.merge_vertices(VERTEX_PAIR_2)
G7.delete_edge(FIRST_EDGE)

# \mathcal{D}_4 = K_5 \two_vertex_identification_minus_shared_edge K_{3,3}
G8 = K5.disjoint_union(K33)
G8.merge_vertices(VERTEX_PAIR_1)
G8.merge_vertices([(0, 1), (1, 3)])
G8.delete_edge(FIRST_EDGE)

G9 = Graph({
    0: [1, 4, 5],
    1: [0, 2, 3],
    2: [1, 4, 5, 6, 7, 8],
    3: [1, 4, 5, 6, 7, 8],
    4: [0, 2, 3],
    5: [0, 2, 3],
    6: [2, 3, 7, 8],
    7: [2, 3, 6, 8],
    8: [2, 3, 6, 7],
})
G10 = Graph({
    0: [1, 4, 5],
    1: [0, 2, 3],
    2: [1, 4, 5, 6, 8, 9],
    3: [1, 4, 5, 6, 8, 9],
    4: [0, 2, 3],
    5: [0, 2, 3],
    6: [2, 3, 7],
    7: [6, 8, 9],
    8: [2, 3, 7],
    9: [2, 3, 7],
})
G11 = Graph({
    0: [3, 4, 5],
    1: [3, 4, 5],
    2: [3, 4, 6, 8, 9],
    3: [0, 1, 2],
    4: [0, 1, 2],
    5: [0, 1, 6, 8, 9],
    6: [2, 5, 7],
    7: [6, 8, 9],
    8: [2, 5, 7],
    9: [2, 5, 7],
})

# \mathcal{F}_6 = K_{3,3} \two_vertex_identification_minus_shared_edge K_{3,3}
G12 = K33.disjoint_union(K33)
G12.merge_vertices(VERTEX_PAIR_1)
G12.merge_vertices([(0, 3), (1, 3)])
G12.delete_edge([(0, 0), (0, 3)])

# K_7 - C_4
G13 = CompleteGraph(7)
G13.delete_edges(CycleGraph(4).edge_iterator())

G14 = Graph({
    0: [1, 3, 4, 5],
    1: [0, 2, 4, 5, 6],
    2: [1, 3, 6, 7],
    3: [0, 2, 4, 6, 7],
    4: [0, 1, 3, 5],
    5: [0, 1, 4, 7],
    6: [1, 2, 3, 7],
    7: [2, 3, 5, 6],
})
G15 = Graph({
    0: [1, 4, 6],
    1: [0, 2, 3, 5, 8],
    2: [1, 4, 7],
    3: [1, 4, 7],
    4: [0, 2, 3, 5, 8],
    5: [1, 4, 6],
    6: [0, 5, 8],
    7: [2, 3, 8],
    8: [1, 4, 6, 7],
})
G16 = Graph({
    0: [1, 4, 6],
    1: [0, 2, 3, 5, 7],
    2: [1, 3, 4, 7],
    3: [1, 2, 4, 7],
    4: [0, 2, 3, 5, 7],
    5: [1, 4, 6],
    6: [0, 5, 7],
    7: [1, 2, 3, 4, 6],
})
G17 = Graph({
    0: [1, 5, 6, 7],
    1: [0, 2, 6, 7],
    2: [1, 3, 4],
    3: [2, 5, 8],
    4: [2, 5, 8],
    5: [0, 3, 4],
    6: [0, 1, 7, 8],
    7: [0, 1, 6, 8],
    8: [3, 4, 6, 7],
})
G18 = Graph({
    0: [1, 5, 7],
    1: [0, 2, 3, 6],
    2: [1, 4, 8],
    3: [1, 4, 8],
    4: [2, 3, 5, 7],
    5: [0, 4, 6],
    6: [1, 5, 7],
    7: [0, 4, 6, 8],
    8: [2, 3, 7],
})
G19 = Graph({
    0: [3, 4, 7],
    1: [3, 4, 7],
    2: [3, 4, 8],
    3: [0, 1, 2, 5],
    4: [0, 1, 2, 6],
    5: [3, 6, 7, 8],
    6: [4, 5, 7, 8],
    7: [0, 1, 5, 6, 8],
    8: [2, 5, 6, 7],
})
G20 = Graph({
    0: [2, 3, 4, 5],
    1: [2, 3, 4, 7],
    2: [0, 1, 8],
    3: [0, 1, 6],
    4: [0, 1, 8],
    5: [0, 8, 9],
    6: [3, 8, 9],
    7: [1, 8, 9],
    8: [2, 4, 5, 6, 7],
    9: [5, 6, 7],
})
G21 = Graph({
    0: [5, 7, 8, 9],
    1: [2, 7, 8, 9],
    2: [1, 3, 4],
    3: [2, 5, 6],
    4: [2, 5, 6, 8],
    5: [0, 3, 4],
    6: [3, 4, 7, 9],
    7: [0, 1, 6],
    8: [0, 1, 4],
    9: [0, 1, 6],
})
G22 = Graph({
    0: [2, 5, 7],
    1: [2, 5, 7],
    2: [0, 1, 3, 4, 8],
    3: [2, 5, 9],
    4: [2, 5, 9],
    5: [0, 1, 3, 4, 6],
    6: [5, 7, 9],
    7: [0, 1, 6, 8],
    8: [2, 7, 9],
    9: [3, 4, 6, 8],
})
G23 = CompleteBipartiteGraph(3, 5)
G24 = Graph({
    0: [2, 3, 4],
    1: [2, 3, 4],
    2: [0, 1, 5],
    3: [0, 1, 6],
    4: [0, 1, 7],
    5: [2, 8, 9],
    6: [3, 8, 9],
    7: [4, 8, 9],
    8: [5, 6, 7],
    9: [5, 6, 7],
})

# K_{4,4} - e
G25 = CompleteBipartiteGraph(4, 4)
G25.delete_edge((0, 4))

# K_{4,5} - 4 K_2
G26 = CompleteBipartiteGraph(4, 5)
G26.delete_edges([(u, u + 4) for u in range(4)])

G27 = Graph({
    0: [1, 3, 4, 5],
    1: [0, 2, 4, 5],
    2: [1, 3, 6, 7],
    3: [0, 2, 6, 7],
    4: [0, 1, 5, 7],
    5: [0, 1, 4, 6],
    6: [2, 3, 5, 7],
    7: [2, 3, 4, 6],
})
G28 = Graph({
    0: [1, 5, 7, 8],
    1: [0, 2, 7, 8],
    2: [1, 3, 4],
    3: [2, 5, 6],
    4: [2, 5, 6, 7],
    5: [0, 3, 4],
    6: [3, 4, 8],
    7: [0, 1, 4, 8],
    8: [0, 1, 6, 7],
})
G29 = Graph({
    0: [2, 7, 8, 9],
    1: [2, 7, 8],
    2: [0, 1, 3],
    3: [2, 4, 5],
    4: [3, 6, 8, 9],
    5: [3, 6, 9],
    6: [4, 5, 7],
    7: [0, 1, 6],
    8: [0, 1, 4],
    9: [0, 4, 5],
})
G30 = CompleteMultipartiteGraph([1,2,2,2])
G31 = Graph({
    0: [1, 3, 4, 5, 7],
    1: [0, 2, 5, 7],
    2: [1, 3, 6],
    3: [0, 2, 4, 7],
    4: [0, 3, 5, 6, 7],
    5: [0, 1, 4, 6, 7],
    6: [2, 4, 5, 7],
    7: [0, 1, 3, 4, 5, 6],
})
G32 = Graph({
    0: [1, 3, 4, 8],
    1: [0, 2, 5, 8],
    2: [1, 3, 6],
    3: [0, 2, 7, 8],
    4: [0, 5, 7],
    5: [1, 4, 6, 8],
    6: [2, 5, 7, 8],
    7: [3, 4, 6, 8],
    8: [0, 1, 3, 5, 6, 7],
})
G33 = Graph({
    0: [1, 4, 8],
    1: [0, 2, 6],
    2: [1, 3, 7, 8],
    3: [2, 4, 7, 8],
    4: [0, 3, 5],
    5: [4, 6, 7, 8],
    6: [1, 5, 7, 8],
    7: [2, 3, 5, 6, 8],
    8: [0, 2, 3, 5, 6, 7],
})
G34 = Graph({
    0: [1, 5, 9],
    1: [0, 2, 7],
    2: [1, 3, 9],
    3: [2, 4, 8],
    4: [3, 5, 9],
    5: [0, 4, 6],
    6: [5, 7, 8, 9],
    7: [1, 6, 8, 9],
    8: [3, 6, 7, 9],
    9: [0, 2, 4, 6, 7, 8],
})
G35 = Graph({
    0: [1, 5, 7],
    1: [0, 2, 6],
    2: [1, 3, 7],
    3: [2, 4, 10],
    4: [3, 5, 7],
    5: [0, 4, 8],
    6: [1, 7, 9],
    7: [0, 2, 4, 6, 8, 10],
    8: [5, 7, 9],
    9: [6, 8, 10],
    10: [3, 7, 9],
})


P2_FORBIDDEN_MINORS = [locals()[f'G{i}'] for i in range(1, 36)]

def get_p2_forbidden_minor(G, **minor_kwargs):
    """
    Check if one of the minimal forbidden minors of the projective plane is an
    induced minor of G.

    INPUT:

    - ``G`` -- Graph

    - ``**minor_kwargs`` -- Optional keyword arguments to be passed to
      :meth:`~Graph.minor`

    OUTPUT:

    Return :meth:`~Graph.minor` output if an element of P2_FORBIDDEN_MINORS is
    a minor of G else None

    EXAMPLES:

    #. The Peterson graph is a known projective planar graph so it doesn't have
       a forbidden minor::

        sage: P = graphs.PetersenGraph()
        sage: _ = get_p2_forbidden_minor(P); type(_)        # long
        <class 'NoneType'>

    #. `K_{4,4}` has a projective plane crossing number of 2. One of the minimal
       forbidden minors is `K_{4,4} - e`, so we get a one-to-one dictionary from
       :meth:`~Graph.minor`::

        sage: K44 = graphs.CompleteBipartiteGraph(4, 4)
        sage: minor_map = get_p2_forbidden_minor(K44); minor_map
        {0: [0], 1: [1], 2: [2], 3: [4], 4: [8], 5: [6], 6: [5], 7: [7]}

    .. SEEALSO::

        - :meth:`~Graph.minor`

        - :meth:`~Graph.is_projective_planar` -- if you only care about whether
          G is projective planar and don't need the output of
          :meth:`~Graph.minor`
    """
    num_verts_G = G.num_verts()
    num_edges_G = G.num_edges()

    for forbidden_minor in P2_FORBIDDEN_MINORS:
        # Can't be a minor if it has more vertices or edges than G
        if (
            forbidden_minor.num_verts() > num_verts_G
            or forbidden_minor.num_edges() > num_edges_G
        ):
            continue

        try:
            minor_map = G.minor(forbidden_minor, **minor_kwargs)
            if minor_map is not None:
                return minor_map
        except ValueError:
            continue

    return None
