# sage_setup: distribution = sagemath-graphs
r"""
Various small graphs

The methods defined here appear in :mod:`sage.graphs.graph_generators`.
"""
# ****************************************************************************
#           Copyright (C) 2006 Robert L. Miller <rlmillster@gmail.com>
#                              and Emily A. Kirkman
#           Copyright (C) 2009 Michael C. Yurko <myurko@gmail.com>
#           Copyright (C) 2024 Janmenjaya Panda <janmenjaya.panda.22@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

# import from Sage library
from sage.graphs.graph import Graph
from sage.rings.rational_field import QQ
from sage.misc.functional import sqrt
from math import sin, cos, pi


# ****************************************************************************
#   Named Graphs
# ****************************************************************************

def HarborthGraph():
    r"""
    Return the Harborth Graph.

    The Harborth graph has 104 edges and 52 vertices, and is the smallest known
    example of a 4-regular matchstick graph. For more information, see the
    :wikipedia:`Harborth_graph`.

    EXAMPLES::

        sage: g = graphs.HarborthGraph(); g
        Harborth Graph: Graph on 52 vertices
        sage: g.is_regular(4)
        True
    """
    g = Graph(':s_OGKI?@_?g[QABAo__YEFCp@?iIEbqHWuWLbbh?}[OfcXpGhNHdYPY_SgdYX]'
              'pZkfJPuo[lfZHys^mFcDs}`pG{UNNgoHC}DIgrI[qjMhTyDQrQlVydrBYmWkn',
              loops=False, multiedges=False)

    g.set_pos(
        {
            0: (51.5, 400.0),
            1: (90.6, 308.0),
            2: (90.6, 492.0),
            3: (129.8, 216.0),
            4: (129.8, 584.0),
            5: (150.7, 387.9),
            6: (150.7, 412.1),
            7: (169.0, 124.0),
            8: (169.0, 676.0),
            9: (189.9, 295.9),
            10: (189.9, 504.1),
            11: (229.1, 203.9),
            12: (229.1, 596.1),
            13: (250.0, 400.0),
            14: (251.4, 180.6),
            15: (251.4, 619.4),
            16: (256.1, 300.2),
            17: (256.1, 499.8),
            18: (259.3, 80.9),
            19: (259.3, 719.1),
            20: (333.8, 237.2),
            21: (333.8, 562.8),
            22: (341.7, 137.5),
            23: (341.7, 662.5),
            24: (350.0, 37.9),
            25: (350.0, 336.0),
            26: (350.0, 464.0),
            27: (350.0, 762.1),
            28: (358.3, 137.5),
            29: (358.3, 662.5),
            30: (366.2, 237.2),
            31: (366.2, 562.8),
            32: (440.7, 80.9),
            33: (440.7, 719.1),
            34: (443.9, 300.2),
            35: (443.9, 499.8),
            36: (448.6, 180.6),
            37: (448.6, 619.4),
            38: (450.0, 400.0),
            39: (470.9, 203.9),
            40: (470.9, 596.1),
            41: (510.1, 295.9),
            42: (510.1, 504.1),
            43: (531.0, 124.0),
            44: (531.0, 676.0),
            45: (549.3, 387.9),
            46: (549.3, 412.1),
            47: (570.2, 216.0),
            48: (570.2, 584.0),
            49: (609.4, 308.0),
            50: (609.4, 492.0),
            51: (648.5, 400.0),
        }
    )
    g.name("Harborth Graph")
    return g


def HarriesGraph(embedding=1):
    r"""
    Return the Harries Graph.

    The Harries graph is a Hamiltonian 3-regular graph on 70 vertices.
    See the :wikipedia:`Harries_graph`.

    The default embedding here is to emphasize the graph's 4 orbits. This graph
    actually has a funny construction. The following procedure gives an idea of
    it, though not all the adjacencies are being properly defined.

    #. Take two disjoint copies of a :meth:`Petersen graph
       <PetersenGraph>`. Their vertices will form an orbit of the final graph.

    #. Subdivide all the edges once, to create 15+15=30 new vertices, which
       together form another orbit.

    #. Create 15 vertices, each of them linked to 2 corresponding vertices of
       the previous orbit, one in each of the two subdivided Petersen graphs. At
       the end of this step all vertices from the previous orbit have degree 3,
       and the only vertices of degree 2 in the graph are those that were just
       created.

    #. Create 5 vertices connected only to the ones from the previous orbit so
       that the graph becomes 3-regular.

    INPUT:

    - ``embedding`` -- integer (default: `1`); two embeddings are available,
      and can be selected by setting ``embedding`` to 1 or 2

    EXAMPLES::

        sage: # needs networkx
        sage: g = graphs.HarriesGraph()
        sage: g.order()
        70
        sage: g.size()
        105
        sage: g.girth()
        10
        sage: g.diameter()
        6
        sage: g.show(figsize=[10, 10])          # long time                             # needs sage.plot
        sage: graphs.HarriesGraph(embedding=2).show(figsize=[10, 10])   # long time, needs sage.plot

    TESTS::

        sage: graphs.HarriesGraph(embedding=3)                                          # needs networkx
        Traceback (most recent call last):
        ...
        ValueError: the value of embedding must be 1 or 2
    """
    from sage.graphs.generators.families import LCFGraph
    g = LCFGraph(70, [-29, -19, -13, 13, 21, -27, 27, 33, -13, 13,
                      19, -21, -33, 29], 5)
    g.name("Harries Graph")

    if embedding == 1:
        gpos = g.get_pos()
        ppos = PetersenGraph().get_pos()

        # The graph's four orbits
        o = [None] * 4
        o[0] = [0, 2, 6, 8, 14, 16, 20, 22, 28, 30, 34, 36, 42, 44, 48, 50,
                56, 58, 62, 64]
        o[1] = [1, 3, 5, 7, 9, 13, 15, 17, 19, 21, 23, 27, 29, 31, 33, 35,
                37, 41, 43, 45, 47, 49, 51, 55, 57, 59, 61, 63, 65, 69]
        o[2] = [60, 10, 12, 4, 24, 26, 18, 38, 40, 32, 52, 54, 46, 66, 68]
        o[3] = [11, 25, 39, 53, 67]

        # Correspondence between the vertices of one of the two Petersen graphs
        # on o[0] and the vertices of a standard Petersen graph object
        g_to_p = {0: 0, 2: 1, 42: 5, 44: 8, 14: 7, 16: 2, 56: 9, 58: 6,
                  28: 4, 30: 3}

        # Correspondence between the vertices of the other Petersen graph on
        # o[0] and the vertices of the first one
        g_to_g = {64: 44, 34: 0, 36: 28, 6: 2, 8: 58, 48: 16, 50: 30,
                  20: 14, 22: 56, 62: 42}

        # Position for the vertices from the first copy
        for v, i in g_to_p.items():
            gpos[v] = ppos[i]

        # Position for the vertices in the second copy. Moves the first, too.
        offset = 3.5
        for v, i in g_to_g.items():
            x, y = gpos[i]
            gpos[v] = (x + offset*0.5, y)
            gpos[i] = (x - offset*0.5, y)

        # Vertices from o[1]. These are actually the "edges" of the copies of
        # Petersen.
        for v in o[1]:
            p1, p2 = (gpos[x] for x in g.neighbors(v) if x in o[0])
            gpos[v] = ((p1[0] + p2[0])/2, (p1[1] + p2[1])/2)

        # 15 vertices from o[2]
        for i, v in enumerate(o[2]):
            gpos[v] = (-1.75 + i*.25, 2)

        # 5 vertices from o[3]
        for i, v in enumerate(o[3]):
            gpos[v] = (-1 + i*.5, 2.5)

        return g

    elif embedding == 2:
        return g
    else:
        raise ValueError("the value of embedding must be 1 or 2")


def HarriesWongGraph(embedding=1):
    r"""
    Return the Harries-Wong Graph.

    See the :wikipedia:`Harries-Wong_graph`.

    *About the default embedding:*

    The default embedding is an attempt to emphasize the graph's 8 (!!!)
    different orbits. In order to understand this better, one can picture the
    graph as being built in the following way.

    #. One first creates a 3-dimensional cube (8 vertices, 12 edges),
       whose vertices define the first orbit of the final graph.

    #. The edges of this graph are subdivided once, to create 12 new vertices
       which define a second orbit.

    #. The edges of the graph are subdivided once more, to create 24 new
       vertices giving a third orbit.

    #. 4 vertices are created and made adjacent to the vertices of the second
       orbit so that they have degree 3. These 4 vertices also define a new
       orbit.

    #. In order to make the vertices from the third orbit 3-regular (they all
       miss one edge), one creates a binary tree on 1 + 3 + 6 + 12 vertices. The
       leaves of this new tree are made adjacent to the 12 vertices of the third
       orbit, and the graph is now 3-regular. This binary tree contributes 4 new
       orbits to the Harries-Wong graph.

    INPUT:

    - ``embedding`` -- integer (default: `1`); two embeddings are available,
      and can be selected by setting ``embedding`` to 1 or 2

    EXAMPLES::

        sage: # needs networkx
        sage: g = graphs.HarriesWongGraph()
        sage: g.order()
        70
        sage: g.size()
        105
        sage: g.girth()
        10
        sage: g.diameter()
        6
        sage: orbits = g.automorphism_group(orbits=True)[-1]    # long time             # needs sage.groups
        sage: g.show(figsize=[15, 15], partition=orbits)        # long time             # needs sage.groups sage.plot

    Alternative embedding::

        sage: graphs.HarriesWongGraph(embedding=2).show()       # long time             # needs networkx sage.plot

    TESTS::

        sage: graphs.HarriesWongGraph(embedding=3)                                      # needs networkx
        Traceback (most recent call last):
        ...
        ValueError: the value of embedding must be 1 or 2
    """
    L = [9, 25, 31, -17, 17, 33, 9, -29, -15, -9, 9, 25, -25, 29, 17, -9,
         9, -27, 35, -9, 9, -17, 21, 27, -29, -9, -25, 13, 19, -9, -33,
         -17, 19, -31, 27, 11, -25, 29, -33, 13, -13, 21, -29, -21, 25,
         9, -11, -19, 29, 9, -27, -19, -13, -35, -9, 9, 17, 25, -9, 9, 27,
         -27, -21, 15, -9, 29, -29, 33, -9, -25]

    from sage.graphs.generators.families import LCFGraph
    g = LCFGraph(70, L, 1)
    g.name("Harries-Wong graph")

    if embedding == 1:
        d = g.get_pos()

        # Binary tree (left side)
        d[66] = (-9.5, 0)
        g._line_embedding([37, 65, 67], first=(-8, 2.25), last=(-8, -2.25))
        g._line_embedding([36, 38, 64, 24, 68, 30], first=(-7, 3),
                          last=(-7, -3))
        g._line_embedding([35, 39, 63, 25, 59, 29, 11, 5, 55, 23, 69, 31],
                          first=(-6, 3.5), last=(-6, -3.5))

        # Cube, corners: [9, 15, 21, 27, 45, 51, 57, 61]
        g._circle_embedding([61, 9], center=(0, -1.5), shift=.2, radius=4)
        g._circle_embedding([27, 15], center=(0, -1.5), shift=.7, radius=4*.707)
        g._circle_embedding([51, 21], center=(0, 2.5), shift=.2, radius=4)
        g._circle_embedding([45, 57], center=(0, 2.5), shift=.7, radius=4*.707)

        # Cube, subdivision
        g._line_embedding([21, 22, 43, 44, 45], first=d[21], last=d[45])
        g._line_embedding([21, 4, 3, 56, 57], first=d[21], last=d[57])
        g._line_embedding([57, 12, 13, 14, 15], first=d[57], last=d[15])
        g._line_embedding([15, 6, 7, 8, 9], first=d[15], last=d[9])
        g._line_embedding([9, 10, 19, 20, 21], first=d[9], last=d[21])
        g._line_embedding([45, 54, 53, 52, 51], first=d[45], last=d[51])
        g._line_embedding([51, 50, 49, 58, 57], first=d[51], last=d[57])
        g._line_embedding([51, 32, 33, 34, 61], first=d[51], last=d[61])
        g._line_embedding([61, 62, 41, 40, 27], first=d[61], last=d[27])
        g._line_embedding([9, 0, 1, 26, 27], first=d[9], last=d[27])
        g._line_embedding([27, 28, 47, 46, 45], first=d[27], last=d[45])
        g._line_embedding([15, 16, 17, 60, 61], first=d[15], last=d[61])

        # Top vertices
        g._line_embedding([2, 18, 42, 48], first=(-1, 7), last=(3, 7))

        return g

    elif embedding == 2:
        return g
    else:
        raise ValueError("the value of embedding must be 1 or 2")


def WellsGraph():
    r"""
    Return the Wells graph.

    For more information on the Wells graph (also called Armanios-Wells graph),
    see `this page <https://www.win.tue.nl/~aeb/graphs/Wells.html>`_.

    The implementation follows the construction given on page 266 of [BCN1989]_.
    This requires to create intermediate graphs and run a small isomorphism
    test, while everything could be replaced by a pre-computed list of edges.
    I believe that it is better to keep "the recipe" in the code, however, as it
    is quite unlikely that this could become the most time-consuming operation
    in any sensible algorithm, and .... "preserves knowledge", which is what
    open-source software is meant to do.

    EXAMPLES::

        sage: g = graphs.WellsGraph(); g
        Wells graph: Graph on 32 vertices
        sage: g.order()
        32
        sage: g.size()
        80
        sage: g.girth()
        5
        sage: g.diameter()
        4
        sage: g.chromatic_number()
        4
        sage: g.is_regular(k=5)
        True
    """
    from .platonic_solids import DodecahedralGraph
    from .basic import CompleteBipartiteGraph

    # Following the construction from the book "Distance-regular graphs"
    dodecahedron = DodecahedralGraph()

    # Vertices at distance 3 in the Dodecahedron
    distance3 = dodecahedron.distance_graph([3])

    # Building the graph whose line graph is the dodecahedron.
    b = CompleteBipartiteGraph(5, 5)
    b.delete_edges([(0, 5), (1, 6), (2, 7), (3, 8), (4, 9)])

    # Computing the isomorphism between the two
    b = b.line_graph(labels=False)
    _, labels = distance3.is_isomorphic(b, certificate=True)

    # The relabeling that the books claims to exist.
    for v, new_name in labels.items():
        x, y = new_name
        labels[v] = (x % 5, y % 5)

    dodecahedron.relabel(labels)

    # Checking that the above computations indeed produces a good labeling.
    for u in dodecahedron:
        for v in dodecahedron:
            if u == v:
                continue

            if (u[0] != v[0]) and (u[1] != v[1]):
                continue

            if dodecahedron.distance(u, v) != 3:
                raise ValueError("there is something wrong going on !")

    # The graph we will return, starting from the dodecahedron
    g = dodecahedron

    # Good ! Now adding 12 new vertices
    for i in range(5):
        g.add_edge((i, '+'), ('inf', '+'))
        g.add_edge((i, '-'), ('inf', '-'))
        for k in range(5):
            if k == i:
                continue
            g.add_edge((i, '+'), (i, k))
            g.add_edge((i, '-'), (k, i))

    g.name("Wells graph")

    # Giving our graph a "not-so-bad" layout
    g.relabel({
            (1, 3): 8, (3, 0): 18, (3, '+'): 22, (2, 1): 13,
            (1, '+'): 10, (0, 3): 2, (2, '+'): 16, ('inf', '-'): 31,
            (4, 0): 24, (1, 2): 7, (4, '+'): 28, (0, '-'): 5,
            (0, 4): 3, (4, 1): 25, (2, '-'): 17, (3, 2): 20,
            (3, '-'): 23, (1, '-'): 11, (1, 4): 9, (2, 3): 14,
            ('inf', '+'): 30, (4, 2): 26, (1, 0): 6, (0, 1): 0,
            (3, 1): 19, (0, 2): 1, (2, 0): 12, (4, '-'): 29,
            (0, '+'): 4, (4, 3): 27, (3, 4): 21, (2, 4): 15})

    p = [(1, 29, 20, 13, 12, 28, 14, 7),
         (2, 5, 30, 23, 18, 4, 31, 22),
         (3, 17, 21, 9, 24, 16, 27, 25),
         (6, 10, 8, 15, 0, 11, 19, 26)]

    g._circle_embedding(p[0], radius=1)
    g._circle_embedding(p[1], radius=.9)
    g._circle_embedding(p[2], radius=.8)
    g._circle_embedding(p[3], radius=.7)

    return g


def Cell600(embedding=1):
    r"""
    Return the 600-Cell graph.

    This is the adjacency graph of the 600-cell. It has 120 vertices and 720
    edges. For more information, see the :wikipedia:`600-cell`.

    INPUT:

    - ``embedding`` -- integer (default: `1`); two different embeddings for a
      plot

    EXAMPLES::

        sage: # long time, needs sage.rings.number_field
        sage: g = graphs.Cell600()
        sage: g.size()
        720
        sage: g.is_regular(12)
        True
        sage: g.is_vertex_transitive()
        True
    """
    from sage.rings.polynomial.polynomial_ring import polygen
    from sage.rings.number_field.number_field import NumberField
    from sage.modules.free_module import VectorSpace
    from sage.groups.perm_gps.permgroup_named import AlternatingGroup

    x = polygen(QQ, 'x')
    K = NumberField(x ** 2 - x - 1, 'f')
    f = K.gen()
    K4 = VectorSpace(K, 4)

    # first 96 vertices
    step = [[a * f / 2, b * K(1) / 2, c * (f - 1) / 2, 0]
            for a in [-1, 1] for b in [-1, 1] for c in [-1, 1]]
    vert96 = [K4([v[s(1) - 1], v[s(2) - 1], v[s(3) - 1], v[s(4) - 1]])
              for v in step for s in AlternatingGroup(4)]

    # 16 more vertices
    vert16 = [K4([K(a) / 2, K(b) / 2, K(c) / 2, K(d) / 2])
              for a in [-1, 1] for b in [-1, 1]
              for c in [-1, 1] for d in [-1, 1]]

    # 8 last vertices
    vert8 = [K4([1, 0, 0, 0]), K4([-1, 0, 0, 0]),
             K4([0, 1, 0, 0]), K4([0, -1, 0, 0]),
             K4([0, 0, 1, 0]), K4([0, 0, -1, 0]),
             K4([0, 0, 0, 1]), K4([0, 0, 0, -1])]

    # all vertices together
    U = vert96 + vert16 + vert8

    g = Graph([list(range(120)),
               lambda i, j: U[i].inner_product(U[j]) == f / 2])

    # Embedding
    if embedding == 1:
        pos = [0, 1, 3, 13, 78, 90, 93, 110, 29, 104, 11, 48, 107, 83, 92, 55,
               32, 16, 117, 24, 26, 56, 52, 47, 75, 72, 66, 112, 27, 115, 21,
               33, 118, 79, 91, 37, 2, 5, 96, 31, 82, 88, 94, 74, 50, 28, 20,
               105, 45, 99, 70, 25, 101, 54, 46, 51, 17, 35, 98, 41, 84, 85,
               87, 73, 18, 6, 9, 97, 65, 103, 95, 36, 100, 23, 8, 43, 68, 76,
               116, 60, 62, 44, 40, 59, 15, 12, 30, 113, 63, 114, 81, 69, 119,
               19, 7, 49, 86, 89, 111, 67, 22, 4, 10, 14, 38, 64, 80, 102, 57,
               108, 34, 61, 106, 42, 58, 39, 77, 71, 109, 53]
    else:
        pos = [0, 1, 2, 3, 4, 6, 7, 8, 10, 13, 14, 21, 37, 103, 36, 65, 113,
               25, 80, 26, 12, 78, 24, 83, 54, 66, 114, 46, 63, 101, 109, 93,
               79, 75, 51, 44, 31, 119, 43, 5, 57, 100, 11, 108, 34, 41, 69,
               96, 82, 116, 68, 64, 47, 102, 52, 35, 17, 76, 110, 38, 84, 85,
               86, 87, 88, 90, 91, 92, 94, 73, 74, 81, 49, 104, 48, 29, 112,
               61, 20, 62, 72, 18, 60, 23, 42, 30, 115, 58, 27, 106, 98, 9, 19,
               15, 39, 56, 67, 118, 55, 89, 45, 107, 95, 99, 70, 53, 33, 111,
               22, 117, 32, 28, 59, 105, 40, 71, 77, 16, 97, 50]

    g._circle_embedding(pos)

    return g


def Cell120():
    r"""
    Return the 120-Cell graph.

    This is the adjacency graph of the 120-cell. It has 600 vertices and 1200
    edges. For more information, see the :wikipedia:`120-cell`.

    EXAMPLES::

        sage: # long time, needs sage.rings.number_field
        sage: g = graphs.Cell120()
        sage: g.size()
        1200
        sage: g.is_regular(4)
        True
        sage: g.is_vertex_transitive()
        True
    """
    from sage.rings.polynomial.polynomial_ring import polygen
    from sage.rings.number_field.number_field import NumberField
    from sage.modules.free_module import VectorSpace
    from sage.groups.perm_gps.permgroup_named import AlternatingGroup
    from sage.combinat.permutation import Permutations

    x = polygen(QQ, 'x')
    K = NumberField(x ** 2 - x - 1, 'f')
    f = K.gen()
    K4 = VectorSpace(K, 4)

    # first 216 vertices
    step = [(0, 0, K(a) * 2, K(b) * 2)
            for a in [-1, 1] for b in [-1, 1]]
    step += [(a * K(1), b * K(1), c * K(1), d * (2 * f - 1))
             for a in [-1, 1] for b in [-1, 1]
             for c in [-1, 1] for d in [-1, 1]]
    step += [(a * (2 - f), b * f, c * f, d * f)
             for a in [-1, 1] for b in [-1, 1]
             for c in [-1, 1] for d in [-1, 1]]
    step += [(a * (f - 1), b * (f - 1), c * (f - 1), d * (f + 1))
             for a in [-1, 1] for b in [-1, 1]
             for c in [-1, 1] for d in [-1, 1]]
    ens1 = frozenset([(v[s(1) - 1], v[s(2) - 1], v[s(3) - 1], v[s(4) - 1])
                      for v in step for s in Permutations(4)])
    vert1 = [K4(w) for w in ens1]

    # 384 more vertices
    step = [(0, a * (2 - f), b * K(1), c * (f + 1))
            for a in [-1, 1] for b in [-1, 1] for c in [-1, 1]]
    step += [(0, a * (f - 1), b * f, c * (2 * f - 1))
             for a in [-1, 1] for b in [-1, 1] for c in [-1, 1]]
    step += [(a * (f - 1), b * K(1), c * f, d * K(2))
             for a in [-1, 1] for b in [-1, 1]
             for c in [-1, 1] for d in [-1, 1]]
    vert2 = [K4([v[s(1) - 1], v[s(2) - 1], v[s(3) - 1], v[s(4) - 1]])
             for v in step for s in AlternatingGroup(4)]

    # all vertices together
    U = vert1 + vert2

    g = Graph([list(range(600)),
               lambda i, j: U[i].inner_product(U[j]) == 6*f-2])

    pos = [0, 1, 3, 5, 6, 7, 8, 9, 11, 12, 14, 15, 16, 17, 20, 21, 23, 24, 25,
           27, 33, 40, 47, 49, 76, 77, 216, 217, 218, 219, 220, 222, 224, 225,
           226, 230, 231, 232, 233, 235, 238, 241, 242, 245, 247, 249, 251, 253,
           260, 261, 211, 66, 26, 307, 598, 305, 187, 374, 311, 205, 296, 108,
           366, 172, 255, 89, 229, 81, 529, 548, 439, 382, 166, 496, 313, 484,
           402, 234, 530, 256, 358, 406, 553, 577, 583, 401, 334, 417, 257, 438,
           373, 544, 509, 365, 378, 487, 377, 390, 349, 325, 65, 78, 184, 13,
           185, 18, 210, 84, 145, 83, 180, 158, 118, 109, 103, 130, 105, 51,
           178, 155, 110, 85, 206, 95, 204, 190, 514, 513, 515, 466, 467, 441,
           442, 587, 585, 576, 565, 564, 566, 540, 506, 436, 435, 424, 507, 543,
           545, 547, 582, 440, 169, 63, 29, 575, 237, 549, 37, 375, 430, 159,
           457, 61, 331, 208, 498, 39, 578, 48, 244, 486, 411, 364, 73, 455,
           321, 240, 381, 542, 243, 500, 343, 333, 271, 518, 552, 357, 314, 299,
           499, 412, 376, 596, 561, 319, 400, 264, 388, 362, 355, 386, 87, 186,
           52, 99, 125, 113, 36, 121, 41, 127, 149, 100, 31, 137, 177, 43, 32,
           45, 62, 191, 188, 106, 195, 141, 142, 96, 489, 491, 490, 475, 474,
           447, 448, 589, 588, 517, 472, 473, 471, 450, 419, 519, 521, 468, 562,
           594, 595, 488, 554, 413, 167, 116, 4, 557, 504, 536, 170, 389, 410,
           128, 559, 203, 348, 147, 477, 22, 516, 162, 423, 266, 274, 320, 144,
           246, 395, 437, 363, 452, 425, 478, 315, 312, 428, 288, 270, 344, 323,
           493, 479, 275, 387, 286, 284, 347, 359, 462, 336, 368, 392, 324, 44,
           75, 69, 46, 57, 138, 35, 80, 88, 199, 70, 152, 161, 181, 34, 207,
           164, 71, 115, 55, 163, 72, 171, 93, 165, 124, 300, 301, 302, 303,
           304, 306, 308, 309, 310, 290, 291, 292, 293, 295, 298, 277, 278, 281,
           283, 285, 287, 265, 272, 273, 19, 10, 107, 223, 418, 221, 67, 338,
           227, 196, 236, 91, 354, 154, 267, 30, 289, 215, 469, 464, 571, 346,
           151, 508, 397, 520, 318, 294, 470, 268, 370, 322, 445, 421, 427, 317,
           394, 597, 269, 570, 337, 460, 497, 353, 342, 523, 341, 330, 361, 385,
           126, 92, 94, 176, 135, 117, 114, 197, 214, 179, 60, 42, 198, 202,
           102, 101, 174, 104, 146, 90, 38, 111, 122, 157, 153, 133, 502, 501,
           503, 550, 551, 573, 574, 431, 429, 420, 433, 432, 434, 456, 494, 568,
           567, 580, 495, 459, 461, 463, 426, 572, 182, 58, 82, 443, 297, 465,
           86, 339, 586, 209, 541, 140, 391, 143, 510, 28, 422, 213, 280, 522,
           591, 352, 120, 563, 405, 276, 345, 458, 279, 512, 379, 393, 259, 482,
           444, 369, 398, 239, 511, 592, 340, 416, 453, 403, 316, 252, 328, 350,
           367, 326, 2, 175, 97, 139, 74, 131, 173, 134, 193, 192, 132, 79, 50,
           200, 64, 150, 201, 194, 212, 183, 54, 56, 98, 123, 112, 156, 525,
           527, 526, 535, 534, 555, 556, 409, 408, 481, 532, 533, 531, 558, 599,
           483, 485, 528, 454, 414, 415, 524, 446, 593, 160, 59, 68, 449, 492,
           476, 148, 329, 590, 119, 451, 189, 360, 53, 537, 129, 480, 136, 579,
           254, 262, 404, 168, 282, 335, 569, 351, 560, 581, 538, 399, 396, 584,
           228, 258, 380, 407, 505, 539, 263, 327, 250, 248, 383, 371, 546, 372,
           356, 332, 384]
    g._circle_embedding(pos)

    return g


def SuzukiGraph():
    r"""
    Return the Suzuki Graph.

    The Suzuki graph has 1782 vertices, and is strongly regular with parameters
    `(1782,416,100,96)`. Known as S.15 in [Hub1975]_.

    .. NOTE::

        It takes approximately 50 seconds to build this graph. Do not be too
        impatient.

    EXAMPLES::

        sage: g = graphs.SuzukiGraph(); g  # optional internet # not tested
        Suzuki graph: Graph on 1782 vertices
        sage: g.is_strongly_regular(parameters=True)  # optional internet # not tested
        (1782, 416, 100, 96)
    """
    from sage.groups.perm_gps.permgroup_named import SuzukiSporadicGroup
    g = Graph()
    g.add_edges(SuzukiSporadicGroup().orbit((1, 2), "OnSets"))
    g.relabel()
    g.name("Suzuki graph")
    return g


def HallJankoGraph(from_string=True):
    r"""
    Return the Hall-Janko graph.

    For more information on the Hall-Janko graph, see the
    :wikipedia:`Hall-Janko_graph`.

    The construction used to generate this graph in Sage is by a 100-point
    permutation representation of the Janko group `J_2`, as described in version
    3 of the ATLAS of Finite Group representations, in particular on the page
    `ATLAS: J2 -- Permutation representation on 100 points
    <http://brauer.maths.qmul.ac.uk/Atlas/v3/permrep/J2G1-p100B0>`_.

    INPUT:

    - ``from_string`` -- boolean (default: ``True``); whether to build the graph
      from its sparse6 string or through GAP. The two methods return the same
      graph though doing it through GAP takes more time.

    EXAMPLES::

        sage: g = graphs.HallJankoGraph()
        sage: g.is_regular(36)
        True
        sage: g.is_vertex_transitive()                                                  # needs sage.groups
        True

    Is it really strongly regular with parameters 14, 12? ::

        sage: nu = set(g.neighbors(0))
        sage: for v in range(1, 100):
        ....:     if v in nu:
        ....:         expected = 14
        ....:     else:
        ....:         expected = 12
        ....:     nv = set(g.neighbors(v))
        ....:     nv.discard(0)
        ....:     if len(nu & nv) != expected:
        ....:         print("Something is wrong here!!!")
        ....:         break

    Some other properties that we know how to check::

        sage: g.diameter()
        2
        sage: g.girth()
        3
        sage: factor(g.characteristic_polynomial())                                     # needs sage.libs.pari sage.modules
        (x - 36) * (x - 6)^36 * (x + 4)^63

    TESTS::

        sage: gg = graphs.HallJankoGraph(from_string=False)  # long time # optional - internet
        sage: g.is_isomorphic(gg)  # long time # optional - internet
        True
    """
    if from_string:
        string = (":~?@c__E@?g?A?w?A@GCA_?CA`OWF`W?EAW?@?_OD@_[GAgcIaGGB@OcIA"
                  "wCE@o_K_?GB@?WGAouC@OsN_?GB@O[GB`A@@_e?@OgLB_{Q_?GC@O[GAOs"
                  "OCWGBA?kKBPA@?_[KB_{OCPKT`o_RD`]A?o[HBOwODW?DA?cIB?wRDP[X`"
                  "ogKB_{QD@]B@o_KBPWXE`mC@o_JB?{PDPq@?oWGA_{OCPKTDp_YEwCA@_c"
                  "IBOwOC`OX_OGB@?WPDPcYFg?C@_gKBp?SE@cYF`{_`?SGAOoOC`_\\FwCE"
                  "A?gKBO{QD@k[FqI??_OFA_oQE@k\\Fq?`GgCB@pGRD@_XFP{a_?SE@ocIA"
                  "ooNCPOUEqU@?oODA?cJB_{UEqYC@_kLC@CREPk]GAGbHgCA@?SMBpCSD`["
                  "YFq?`Ga]BA?gPC`KSD`_\\Fa?cHWGB@?[IAooPD`[WF@s^HASeIg?@@OcP"
                  "C`KYF@w^GQ[h`O[HAooMC@CQCpSVEPk\\GaSeIG?FA?kLB_{OC`OVE@cYG"
                  "QUA@?WLBp?PC`KVEqKgJg?DA?sMBpCSDP[WEQKfIay@?_KD@_[GC`SUE@k"
                  "[FaKdHa[k_?OLC@CRD@WVEpo^HAWfIAciIqoo_?CB@?kMCpOUE`o\\GAKg"
                  "IQgq_?GD@_[GB?{OCpWVE@cYFACaHAWhJR?q_?CC@_kKBpC\\GACdHa[kJ"
                  "a{o_?CA?oOFBpGRD@o\\GaKdIQonKrOt_?WHA`?PC`KTD`k]FqSeIaolJr"
                  "CqLWCA@OkKCPGRDpcYGAKdIAgjJAsmJr?t__OE@ogJB_{XEps`HA[gIQwn"
                  "KWKGAOoMBpGUE`k[Fa?aHqckJbSuLw?@?_SHA_kLC@OTFPw^GaOkLg?B@?"
                  "[HA_{PDP_XFaCbHa[gIqooKRWx_?CFBpOTE@cZFPw^GACcHQgoKrSvMwWG"
                  "BOwQCp_YFP{`HASfJAwnKRSx_OSSDP[WEq?aGqSfIQsoKR_zNWCE@o_HA_"
                  "sREPg^GAGcHQWfIAciKbOxNg?A@__IAooMC`KTD`g\\GAKcIasoKrOtLb["
                  "wMbyCA?cKBp?TD`[WE`s^GQGbHqcjJrK{NRw~_oODA?sNC@CQCpOZF@s]G"
                  "QOfIaolJrGsLbk}_?OFA_sRD@SVE`k[HQcjJa{qLb[xMb|?_OOFA?cIAos"
                  "RDP_ZFa?aGqOfIAsuMbk{Ns@@OsQAA_sPDPWXE`o\\FqKdIQkkJrCuLr_x"
                  "Mro}NsDAPG?@@OWFApKUE@o`IQolKRKsLrc|NsQC@OWGAOgJCpOWE`o_GQ"
                  "KiIqwnKr_~OcLCPS]A?oWHA_oMBpKSDP[\\FagjKBWxMbk{OSQ@@O_IAoo"
                  "LBpCSD`g\\FaGbHQWgIQgmKRKwMRl?PgGC@OWHB@KSE@c[FqCaGqSeIAkk"
                  "KBCqLBSuMBpGQWCA@?cKBOwRDPWVE@k^GqOfJr?pKbKtLrs}OSHDQwKIBO"
                  "wPD@WWEQ?`HQWfIQglKBOtLbo}Ns@@OsTE_?kLCpWWHA[gIqomKBGwMRgz"
                  "NBw~OSPDPc\\H_?CFAOoLCPSVE`o\\GAOeJAwpKbKtMrx?Qcq??OKFA?gJ"
                  "B`?QDpcYEpo]FqKfIAgjJB?qKr_{NS@A__SE@o_HBO{PC`OTD`{_HaciIq"
                  "{vMbt?OcPFQCeB@?SKBOwRD@SXE`k[FPw`HQ_lKRKxNRxBPC\\HQclK_?K"
                  "EB?sOC`OTDa?`GqWgJRCrNBw~OSHFQStMRtDQ_?KC@OoQE`k_GaOdHa[gI"
                  "q{tMBg|Nb|?OcPMSDDQSwCB@_cJB_{OCpOVFP{dHa[jJQwqKrk}NsHBQCd"
                  "MRtMA?oSEA_wPDp_YEpo]GAOeIq{pLBk}NsLEQCtNTDU??OKEA_oLC@[[G"
                  "aKnKBOtLbk~OCPFQStNSDLSTgGKC@GSD`[WEpw_GQGcIAciJAwpKb_xMbk"
                  "~QShJRc|R`_wNCPcZF@s^GAGbHA_hJR?qKrOvMRg|NsDEPsxTTgCB@?gJB"
                  "?sMC@CUDp_]FqCaHQcjJQwtLrhCPS\\IRCtQTw?B@?SHA_wPC`_aGqOiJa"
                  "{oKRKvMRpFQChKRtXVUTi??ocNC@KUE@cYFaGdHa_mJrKsLb[yMro|OcXI"
                  "RdPTTddZaOgJB@?UEPk[FQCfIaolJrSvMBczNR|AOsXFQCtOTtaB@?WGAP"
                  "?TEPo\\GAGdHqgmKBCqLR[xMb|?PC`HQs|TTt`XUtu@?o[HB?sNCPGXF@{"
                  "_GQKcIqolJb_yNCLDPs`MRtDRTTdYUwSEA?kLB`CWF@s]FqGgIqooLRgzN"
                  "RxFQSlMSDDQTDXVUTi@?_KDAOoLBpKUEQOfIa{oLB_xMrt?Os\\HQcpMST"
                  "HSTtl[VT}A@ocJBOwSD`_XEpo_Ha_mJrKtLbgzNSTGQspLRtDUUDp\\WG["
                  "HB`CQCp[WFQGgIQgkJQ{rLbc{Nc@APsdLRt@PSt\\WUtt_Wn")
        g = Graph(string, loops=False, multiedges=False)
    else:
        # The following construction is due to version 3 of the ATLAS of
        # Finite Group Representations, specifically the page at
        # http://brauer.maths.qmul.ac.uk/Atlas/v5/permrep/J2G1-p100B0 .

        from sage.libs.gap.libgap import libgap
        libgap.load_package("AtlasRep")  # representation of HJ on 100 points
        G = libgap.AtlasGroup("HJ", libgap.NrMovedPoints, 100)
        edges = G.Orbit([1, 5], libgap.OnSets)
        g = Graph()
        g.add_edges(edges)
        g.relabel()

    g._circle_embedding(list(range(100)))
    g.name("Hall-Janko graph")
    return g


def Balaban10Cage(embedding=1):
    r"""
    Return the Balaban 10-cage.

    The Balaban 10-cage is a 3-regular graph with 70 vertices and 105 edges. See
    the :wikipedia:`Balaban_10-cage`.

    The default embedding gives a deeper understanding of the graph's
    automorphism group. It is divided into 4 layers (each layer being a set of
    points at equal distance from the drawing's center). From outside to inside:

    - L1: The outer layer (vertices which are the furthest from the origin) is
      actually the disjoint union of two cycles of length 10.

    - L2: The second layer is an independent set of 20 vertices.

    - L3: The third layer is a matching on 10 vertices.

    - L4: The inner layer (vertices which are the closest from the origin) is
      also the disjoint union of two cycles of length 10.

    This graph is not vertex-transitive, and its vertices are partitioned into 3
    orbits: L2, L3, and the union of L1 of L4 whose elements are equivalent.

    INPUT:

    - ``embedding`` -- integer (default: `1`); two embeddings are available,
      and can be selected by setting ``embedding`` to be either 1 or 2

    EXAMPLES::

        sage: # needs networkx
        sage: g = graphs.Balaban10Cage()
        sage: g.girth()
        10
        sage: g.chromatic_number()
        2
        sage: g.diameter()
        6
        sage: g.is_hamiltonian()                                                        # needs sage.numerical.mip
        True
        sage: g.show(figsize=[10,10])           # long time                             # needs sage.plot

    TESTS::

        sage: graphs.Balaban10Cage(embedding='foo')                                     # needs networkx
        Traceback (most recent call last):
        ...
        ValueError: the value of embedding must be 1 or 2
    """
    L = [-9, -25, -19, 29, 13, 35, -13, -29, 19, 25, 9, -29, 29, 17, 33,
         21, 9, -13, -31, -9, 25, 17, 9, -31, 27, -9, 17, -19, -29, 27,
         -17, -9, -29, 33, -25, 25, -21, 17, -17, 29, 35, -29, 17, -17,
         21, -25, 25, -33, 29, 9, 17, -27, 29, 19, -17, 9, -27, 31, -9,
         -17, -25, 9, 31, 13, -9, -21, -33, -17, -29, 29]

    from sage.graphs.generators.families import LCFGraph
    g = LCFGraph(70, L, 1)
    g.name("Balaban 10-cage")

    if embedding == 2:
        return g
    elif embedding != 1:
        raise ValueError("the value of embedding must be 1 or 2")

    L3 = [5, 24, 35, 46, 29, 40, 51, 34, 45, 56]
    g._circle_embedding(L3, center=(0, 0), radius=4.3)

    L2 = [6, 4, 23, 25, 60, 36, 1, 47, 28, 30, 39, 41, 50, 52, 33, 9, 44,
          20, 55, 57]
    g._circle_embedding(L2, center=(0, 0), radius=5, shift=-.5)

    L1a = [69, 68, 67, 66, 65, 64, 63, 62, 61, 0]
    L1b = [19, 18, 17, 16, 15, 14, 13, 12, 11, 10]
    g._circle_embedding(L1a, center=(0, 0), radius=6, shift=3.25)
    g._circle_embedding(L1b, center=(0, 0), radius=6, shift=-1.25)

    L4a = [37, 2, 31, 38, 53, 32, 21, 54, 3, 22]
    g._circle_embedding(L4a, center=(0, 0), radius=3, shift=1.9)

    L4b = [26, 59, 48, 27, 42, 49, 8, 43, 58, 7]
    g._circle_embedding(L4b, center=(0, 0), radius=3, shift=1.1)

    return g


def Balaban11Cage(embedding=1):
    r"""
    Return the Balaban 11-cage.

    For more information, see the :wikipedia:`Balaban_11-cage`.

    INPUT:

    - ``embedding`` -- integer (default: `1`); three embeddings are available,
      and can be selected by setting ``embedding`` to be 1, 2, or 3

      - The first embedding is the one appearing on page 9 of the Fifth Annual
        Graph Drawing Contest report [EMMN1998]_. It separates vertices based on
        their eccentricity (see :meth:`eccentricity()
        <sage.graphs.generic_graph.GenericGraph.eccentricity>`).

      - The second embedding has been produced just for Sage and is meant to
        emphasize the automorphism group's 6 orbits.

      - The last embedding is the default one produced by the :meth:`LCFGraph`
        constructor.

    .. NOTE::

        The vertex labeling changes according to the value of ``embedding=1``.

    EXAMPLES:

    Basic properties::

        sage: g = graphs.Balaban11Cage()
        sage: g.order()
        112
        sage: g.size()
        168
        sage: g.girth()
        11
        sage: g.diameter()
        8
        sage: g.automorphism_group().cardinality()                                      # needs sage.groups
        64

    Our many embeddings::

        sage: g1 = graphs.Balaban11Cage(embedding=1)
        sage: g2 = graphs.Balaban11Cage(embedding=2)                                    # needs networkx
        sage: g3 = graphs.Balaban11Cage(embedding=3)                                    # needs networkx
        sage: g1.show(figsize=[10,10])          # long time                             # needs sage.plot
        sage: g2.show(figsize=[10,10])          # long time                             # needs networkx sage.plot
        sage: g3.show(figsize=[10,10])          # long time                             # needs networkx sage.plot

    Proof that the embeddings are the same graph::

        sage: g1.is_isomorphic(g2)  # g2 and g3 are obviously isomorphic                # needs networkx
        True

    TESTS::

        sage: graphs.Balaban11Cage(embedding='xyzzy')
        Traceback (most recent call last):
        ...
        ValueError: the value of embedding must be 1, 2, or 3
    """
    if embedding == 1:
        pos_dict = {}
        for j in range(8):
            for i in range(8):
                pos_dict[str(j) + str(i)] = [
                        0.8 * float(cos(2*((8*j + i)*pi/64 + pi/128))),
                        0.8 * float(sin(2*((8*j + i)*pi/64 + pi/128)))
                ]
            for i in range(4):
                pos_dict['1' + str(j) + str(i)] = [
                        1.1 * float(cos(2*((4*j + i)*pi/32 + pi/64))),
                        1.1 * float(sin(2*((4*j + i)*pi/32 + pi/64)))
                ]
            for i in range(2):
                pos_dict['1' + str(j) + str(i + 4)] = [
                        1.4 * float(cos(2*((2*j + i)*pi/16 + pi/32))),
                        1.4 * float(sin(2*((2*j + i)*pi/16 + pi/32)))
                ]

        edge_dict = {
            "00": ["11"], "01": ["10"],   "02": ["53"], "03": ["52"],
            "11": ["20"], "10": ["21"],   "53": ["22"], "52": ["23"],
            "20": ["31"], "21": ["30"],   "22": ["33"], "23": ["32"],
            "31": ["40"], "30": ["41"],   "33": ["43"], "32": ["42"],
            "40": ["50"], "41": ["51"],   "43": ["12"], "42": ["13"],
            "50": ["61"], "51": ["60"],   "12": ["63"], "13": ["62"],
            "61": ["70"], "60": ["71"],   "63": ["72"], "62": ["73"],
            "70": ["01"], "71": ["00"],   "72": ["03"], "73": ["02"],

            "04": ["35"], "05": ["34"],   "06": ["37"], "07": ["36"],
            "35": ["64"], "34": ["65"],   "37": ["66"], "36": ["67"],
            "64": ["55"], "65": ["54"],   "66": ["17"], "67": ["16"],
            "55": ["45"], "54": ["44"],   "17": ["46"], "16": ["47"],
            "45": ["74"], "44": ["75"],   "46": ["76"], "47": ["77"],
            "74": ["25"], "75": ["24"],   "76": ["27"], "77": ["26"],
            "25": ["14"], "24": ["15"],   "27": ["56"], "26": ["57"],
            "14": ["05"], "15": ["04"],   "56": ["07"], "57": ["06"],

            "100": ["03", "04"],   "110": ["10", "12"],
            "101": ["01", "06"],   "111": ["11", "13"],
            "102": ["00", "07"],   "112": ["14", "16"],
            "103": ["02", "05"],   "113": ["15", "17"],

            "120": ["22", "24"],   "130": ["33", "36"],
            "121": ["20", "26"],   "131": ["32", "37"],
            "122": ["21", "27"],   "132": ["31", "34"],
            "123": ["23", "25"],   "133": ["30", "35"],

            "140": ["43", "45"],   "150": ["50", "52"],
            "141": ["40", "46"],   "151": ["51", "53"],
            "142": ["41", "47"],   "152": ["54", "56"],
            "143": ["42", "44"],   "153": ["55", "57"],

            "160": ["60", "66"],   "170": ["73", "76"],
            "161": ["63", "65"],   "171": ["72", "77"],
            "162": ["62", "64"],   "172": ["71", "74"],
            "163": ["61", "67"],   "173": ["70", "75"],

            "104": ["100", "102", "105"],   "114": ["110", "111", "115"],
            "105": ["101", "103", "104"],   "115": ["112", "113", "114"],

            "124": ["120", "121", "125"],   "134": ["130", "131", "135"],
            "125": ["122", "123", "124"],   "135": ["132", "133", "134"],

            "144": ["140", "141", "145"],   "154": ["150", "151", "155"],
            "145": ["142", "143", "144"],   "155": ["152", "153", "154"],

            "164": ["160", "161", "165"],   "174": ["170", "171", "175"],
            "165": ["162", "163", "164"],   "175": ["172", "173", "174"]
        }

        return Graph(edge_dict, pos=pos_dict, name="Balaban 11-cage")

    elif embedding == 2 or embedding == 3:
        L = [44, 26, -47, -15, 35, -39, 11, -27, 38, -37, 43, 14, 28, 51,
             -29, -16, 41, -11, -26, 15, 22, -51, -35, 36, 52, -14, -33,
             -26, -46, 52, 26, 16, 43, 33, -15, 17, -53, 23, -42, -35, -28,
             30, -22, 45, -44, 16, -38, -16, 50, -55, 20, 28, -17, -43,
             47, 34, -26, -41, 11, -36, -23, -16, 41, 17, -51, 26, -33,
             47, 17, -11, -20, -30, 21, 29, 36, -43, -52, 10, 39, -28, -17,
             -52, 51, 26, 37, -17, 10, -10, -45, -34, 17, -26, 27, -21,
             46, 53, -10, 29, -50, 35, 15, -47, -29, -41, 26, 33, 55, -17,
             42, -26, -36, 16]

        from sage.graphs.generators.families import LCFGraph
        g = LCFGraph(112, L, 1)
        g.name("Balaban 11-cage")

        if embedding == 3:
            return g

        v1 = [34, 2, 54, 43, 66, 20, 89, 100, 72, 76, 6, 58, 16, 78, 74,
              70, 36, 94, 27, 25, 10, 8, 45, 60, 14, 64, 80, 82, 109, 107,
              49, 98]
        v2 = [88, 3, 19, 55, 67, 42, 101, 33, 77, 5, 17, 57, 69, 71, 73,
              75, 11, 61, 28, 9, 37, 26, 46, 95, 13, 63, 81, 83, 108, 106,
              48, 97]
        l1 = [35, 93, 1, 24, 53, 7, 44, 59, 15, 65, 79, 21, 110, 90, 50,
              99]
        l2 = [87, 4, 18, 56, 68, 41, 102, 32, 12, 62, 29, 84, 38, 105, 47,
              96]

        d = g.get_pos()
        for i, v in enumerate(v1):
            d[v] = (-2, 16.5 - i)

        for i, v in enumerate(l1):
            d[v] = (-10, 8 - i)

        for i, v in enumerate(l2):
            d[v] = (10, 8.5 - i)

        for i, v in enumerate(v2):
            d[v] = (2, 16.5 - i)

        for i, v in enumerate([0, 111, 92, 91, 52, 51, 23, 22]):
            d[v] = (-20, 14.5 - 4*i)

        for i, v in enumerate([104, 103, 86, 85, 40, 39, 31, 30]):
            d[v] = (20, 14.5 - 4*i)

        return g

    else:
        raise ValueError("the value of embedding must be 1, 2, or 3")


def BidiakisCube():
    r"""
    Return the Bidiakis cube.

    For more information, see the :wikipedia:`Bidiakis_cube`.

    EXAMPLES:

    The Bidiakis cube is a 3-regular graph having 12 vertices and 18 edges. This
    means that each vertex has a degree of 3::

        sage: g = graphs.BidiakisCube(); g
        Bidiakis cube: Graph on 12 vertices
        sage: g.show()                          # long time                             # needs sage.plot
        sage: g.order()
        12
        sage: g.size()
        18
        sage: g.is_regular(3)
        True

    It is a Hamiltonian graph with diameter 3 and girth 4::

        sage: g.is_hamiltonian()                                                        # needs sage.numerical.mip
        True
        sage: g.diameter()
        3
        sage: g.girth()
        4

    It is a planar graph with characteristic polynomial
    `(x - 3) (x - 2) (x^4) (x + 1) (x + 2) (x^2 + x - 4)^2` and
    chromatic number 3::

        sage: g.is_planar()
        True
        sage: char_poly = g.characteristic_polynomial()                                 # needs sage.modules
        sage: x = char_poly.parent()('x')                                               # needs sage.modules
        sage: char_poly == (x - 3) * (x - 2) * (x^4) * (x + 1) * (x + 2) * (x^2 + x - 4)^2          # needs sage.modules
        True
        sage: g.chromatic_number()                                                      # needs sage.modules
        3
    """
    edge_dict = {
        0: [1, 6, 11], 1: [2, 5], 2: [3, 10], 3: [4, 9], 4: [5, 8],
        5: [6], 6: [7], 7: [8, 11], 8: [9], 9: [10], 10: [11]}
    g = Graph(edge_dict, format='dict_of_lists', name="Bidiakis cube")
    g._circle_embedding(range(12), angle=pi/2)
    return g


def BiggsSmithGraph(embedding=1):
    r"""
    Return the Biggs-Smith graph.

    For more information, see the :wikipedia:`Biggs-Smith_graph`.

    INPUT:

    - ``embedding`` -- integer (default: `1`); two embeddings are available,
      and can be selected by setting ``embedding`` to be 1 or 2

    EXAMPLES:

    Basic properties::

        sage: # needs networkx
        sage: g = graphs.BiggsSmithGraph()
        sage: g.order()
        102
        sage: g.size()
        153
        sage: g.girth()
        9
        sage: g.diameter()
        7
        sage: g.automorphism_group().cardinality()      # long time
        2448
        sage: g.show(figsize=[10, 10])          # long time                             # needs sage.plot

    The other embedding::

        sage: graphs.BiggsSmithGraph(embedding=2).show()        # long time             # needs networkx

    TESTS::

        sage: graphs.BiggsSmithGraph(embedding='xyzzy')                                 # needs networkx
        Traceback (most recent call last):
        ...
        ValueError: the value of embedding must be 1 or 2
    """
    L = [16, 24, -38, 17, 34, 48, -19, 41, -35, 47, -20, 34, -36,
         21, 14, 48, -16, -36, -43, 28, -17, 21, 29, -43, 46, -24,
         28, -38, -14, -50, -45, 21, 8, 27, -21, 20, -37, 39, -34,
         -44, -8, 38, -21, 25, 15, -34, 18, -28, -41, 36, 8, -29,
         -21, -48, -28, -20, -47, 14, -8, -15, -27, 38, 24, -48, -18,
         25, 38, 31, -25, 24, -46, -14, 28, 11, 21, 35, -39, 43, 36,
         -38, 14, 50, 43, 36, -11, -36, -24, 45, 8, 19, -25, 38, 20,
         -24, -14, -21, -8, 44, -31, -38, -28, 37]

    from sage.graphs.generators.families import LCFGraph
    g = LCFGraph(102, L, 1)
    g.name("Biggs-Smith graph")

    if embedding == 1:

        orbs = [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 0],
                [17, 101, 25, 66, 20, 38, 53, 89, 48, 75, 56, 92, 45, 78,
                 34, 28, 63],
                [18, 36, 26, 65, 19, 37, 54, 90, 47, 76, 55, 91, 46, 77,
                 35, 27, 64],
                [21, 39, 52, 88, 49, 74, 57, 93, 44, 79, 33, 29, 62, 83,
                 100, 24, 67],
                [22, 97, 51, 96, 50, 95, 58, 94, 59, 80, 60, 81, 61, 82,
                 99, 23, 98],
                [30, 86, 84, 72, 70, 68, 42, 40, 31, 87, 85, 73, 71, 69,
                 43, 41, 32]]

        # central orbits
        g._circle_embedding(orbs[1], center=(-.4, 0), radius=.2)
        g._circle_embedding(orbs[3], center=(.4, 0), radius=.2, shift=4)

        # lower orbits
        g._circle_embedding(orbs[0], center=(-.9, -.5), radius=.3, shift=2)
        g._circle_embedding(orbs[2], center=(-.9, .5), radius=.3)

        # upper orbits
        g._circle_embedding(orbs[4], center=(.9, -.5), radius=.3, shift=4)
        g._circle_embedding(orbs[5], center=(.9, .5), radius=.3, shift=-2)

    elif embedding == 2:
        pass
    else:
        raise ValueError("the value of embedding must be 1 or 2")

    return g


def BlanusaFirstSnarkGraph():
    r"""
    Return the first Blanusa Snark Graph.

    The Blanusa graphs are two snarks on 18 vertices and 27 edges. For more
    information on them, see the :wikipedia:`Blanusa_snarks`.

    .. SEEALSO::

        * :meth:`~sage.graphs.graph_generators.GraphGenerators.BlanusaSecondSnarkGraph`.

    EXAMPLES::

        sage: g = graphs.BlanusaFirstSnarkGraph()
        sage: g.order()
        18
        sage: g.size()
        27
        sage: g.diameter()
        4
        sage: g.girth()
        5
        sage: g.automorphism_group().cardinality()                                      # needs sage.groups
        8
    """
    g = Graph({17: [4, 7, 1], 0: [5], 3: [8], 13: [9], 12: [16],
               10: [15], 11: [6], 14: [2]},
              name="Blanusa First Snark Graph")

    g.add_cycle(list(range(17)))
    g._circle_embedding(list(range(17)), shift=0.25)
    g._pos[17] = (0, 0)
    return g


def BlanusaSecondSnarkGraph():
    r"""
    Return the second Blanusa Snark Graph.

    The Blanusa graphs are two snarks on 18 vertices and 27 edges. For more
    information on them, see the :wikipedia:`Blanusa_snarks`.

    .. SEEALSO::

        * :meth:`~sage.graphs.graph_generators.GraphGenerators.BlanusaFirstSnarkGraph`.

    EXAMPLES::

        sage: g = graphs.BlanusaSecondSnarkGraph()
        sage: g.order()
        18
        sage: g.size()
        27
        sage: g.diameter()
        4
        sage: g.girth()
        5
        sage: g.automorphism_group().cardinality()                                      # needs sage.groups
        4
    """
    c0 = (-1, 0)
    c1 = (-1, 1)
    g = Graph({c0: [(0, 0), (1, 4), c1], c1: [(0, 3), (1, 1)],
               (0, 2): [(0, 5)], (0, 6): [(0, 4)],
               (0, 7): [(0, 1)], (1, 7): [(1, 2)],
               (1, 0): [(1, 6)], (1, 3): [(1, 5)]},
              name="Blanusa Second Snark Graph")

    g.add_cycle([(0, i) for i in range(5)])
    g.add_cycle([(1, i) for i in range(5)])
    g.add_cycle([(0, 5), (0, 6), (0, 7), (1, 5), (1, 6), (1, 7)])

    g._circle_embedding([(0, (2 * i) % 5) for i in range(5)],
                        center=(-1.5, 0),
                        shift=.5)
    g._circle_embedding([(1, (2 * i) % 5) for i in range(5)],
                        center=(1.5, 0))

    g._circle_embedding([(0, i) for i in range(5, 8)] + [c0] * 4,
                        center=(-1.2, 0),
                        shift=2.5,
                        radius=2.2)
    g._circle_embedding([(1, i) for i in range(5, 8)] + [c0] * 4,
                        center=(1.2, 0),
                        shift=-1,
                        radius=2.2)

    g._circle_embedding([c0, c1], shift=.5)
    g.relabel()
    return g


def BrinkmannGraph():
    r"""
    Return the Brinkmann graph.

    For more information, see the :wikipedia:`Brinkmann_graph`.

    EXAMPLES:

    The Brinkmann graph is a 4-regular graph having 21 vertices and 42
    edges. This means that each vertex has degree 4::

        sage: G = graphs.BrinkmannGraph(); G
        Brinkmann graph: Graph on 21 vertices
        sage: G.show()                          # long time                             # needs sage.plot
        sage: G.order()
        21
        sage: G.size()
        42
        sage: G.is_regular(4)
        True

    It is an Eulerian graph with radius 3, diameter 3, and girth 5::

        sage: G.is_eulerian()
        True
        sage: G.radius()
        3
        sage: G.diameter()
        3
        sage: G.girth()
        5

    The Brinkmann graph is also Hamiltonian with chromatic number 4::

        sage: G.is_hamiltonian()                                                        # needs sage.numerical.mip
        True
        sage: G.chromatic_number()
        4

    Its automorphism group is isomorphic to `D_7`::

        sage: ag = G.automorphism_group()                                               # needs sage.groups
        sage: ag.is_isomorphic(DihedralGroup(7))                                        # needs sage.groups
        True
    """
    edge_dict = {
        0: [2, 5, 7, 13],
        1: [3, 6, 7, 8],
        2: [4, 8, 9],
        3: [5, 9, 10],
        4: [6, 10, 11],
        5: [11, 12],
        6: [12, 13],
        7: [15, 20],
        8: [14, 16],
        9: [15, 17],
        10: [16, 18],
        11: [17, 19],
        12: [18, 20],
        13: [14, 19],
        14: [17, 18],
        15: [18, 19],
        16: [19, 20],
        17: [20]}
    g = Graph(edge_dict, format='dict_of_lists', name="Brinkmann graph")
    g._circle_embedding(range(7), radius=4, angle=pi/2)
    g._circle_embedding(range(7, 14), radius=2, angle=pi/2 + pi/7)
    g._circle_embedding(range(14, 21), radius=1, angle=pi/2 + pi/7)
    return g


def BrouwerHaemersGraph():
    r"""
    Return the Brouwer-Haemers Graph.

    The Brouwer-Haemers is the only strongly regular graph of parameters
    `(81,20,1,6)`. It is build in Sage as the Affine Orthogonal graph
    `VO^-(6,3)`. For more information on this graph, see its `corresponding page
    on Andries Brouwer's website
    <https://www.win.tue.nl/~aeb/graphs/Brouwer-Haemers.html>`_.

    EXAMPLES::

        sage: g = graphs.BrouwerHaemersGraph(); g                                       # needs sage.modules
        Brouwer-Haemers: Graph on 81 vertices

    It is indeed strongly regular with parameters `(81,20,1,6)`::

        sage: g.is_strongly_regular(parameters=True)    # long time                     # needs sage.modules sage.rings.finite_rings
        (81, 20, 1, 6)

    Its has as eigenvalues `20,2` and `-7`::

        sage: set(g.spectrum()) == {20,2,-7}                                            # needs sage.modules sage.rings.finite_rings
        True
    """
    from sage.rings.finite_rings.finite_field_constructor import FiniteField
    from sage.modules.free_module import VectorSpace
    from sage.matrix.constructor import Matrix
    from sage.matrix.constructor import identity_matrix

    d = 4
    q = 3
    F = FiniteField(q, "x")
    V = VectorSpace(F, d)
    M = Matrix(F, identity_matrix(d))
    M[1, 1] = -1
    G = Graph([[tuple(_) for _ in V], lambda x, y: (V(x) - V(y))*(M*(V(x) - V(y))) == 0], loops=False)
    G.relabel()
    ordering = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
                18, 19, 20, 21, 22, 23, 24, 25, 26, 48, 49, 50, 51, 52, 53,
                45, 46, 47, 30, 31, 32, 33, 34, 35, 27, 28, 29, 39, 40, 41,
                42, 43, 44, 36, 37, 38, 69, 70, 71, 63, 64, 65, 66, 67, 68,
                78, 79, 80, 72, 73, 74, 75, 76, 77, 60, 61, 62, 54, 55, 56,
                57, 58, 59]
    G._circle_embedding(ordering)
    G.name("Brouwer-Haemers")
    return G


def BuckyBall():
    r"""
    Return the Bucky Ball graph.

    This graph is a 3-regular 60-vertex planar graph. Its vertices and edges
    correspond precisely to the carbon atoms and bonds in buckminsterfullerene.
    When embedded on a sphere, its 12 pentagon and 20 hexagon faces are arranged
    exactly as the sections of a soccer ball.

    EXAMPLES:

    The Bucky Ball is planar::

        sage: g = graphs.BuckyBall()
        sage: g.is_planar()
        True

    The Bucky Ball can also be created by extracting the 1-skeleton of the Bucky
    Ball polyhedron, but this is much slower::

        sage: # needs sage.geometry.polyhedron sage.groups sage.rings.number_field
        sage: g = polytopes.buckyball().vertex_graph()
        sage: g.remove_loops()
        sage: h = graphs.BuckyBall()
        sage: g.is_isomorphic(h)
        True

    The graph is returned along with an attractive embedding::

        sage: g = graphs.BuckyBall()  # long time
        sage: g.plot(vertex_labels=False, vertex_size=10).show()        # long time, needs sage.plot
    """
    edges = [(0, 2), (0, 48), (0, 59), (1, 3), (1, 9), (1, 58),
             (2, 3), (2, 36), (3, 17), (4, 6), (4, 8), (4, 12),
             (5, 7), (5, 9), (5, 16), (6, 7), (6, 20), (7, 21),
             (8, 9), (8, 56), (10, 11), (10, 12), (10, 20), (11, 27),
             (11, 47), (12, 13), (13, 46), (13, 54), (14, 15), (14, 16),
             (14, 21), (15, 25), (15, 41), (16, 17), (17, 40), (18, 19),
             (18, 20), (18, 26), (19, 21), (19, 24), (22, 23), (22, 31),
             (22, 34), (23, 25), (23, 38), (24, 25), (24, 30), (26, 27),
             (26, 30), (27, 29), (28, 29), (28, 31), (28, 35), (29, 44),
             (30, 31), (32, 34), (32, 39), (32, 50), (33, 35), (33, 45),
             (33, 51), (34, 35), (36, 37), (36, 40), (37, 39), (37, 52),
             (38, 39), (38, 41), (40, 41), (42, 43), (42, 46), (42, 55),
             (43, 45), (43, 53), (44, 45), (44, 47), (46, 47), (48, 49),
             (48, 52), (49, 53), (49, 57), (50, 51), (50, 52), (51, 53),
             (54, 55), (54, 56), (55, 57), (56, 58), (57, 59), (58, 59)
             ]
    g = Graph()
    g.add_edges(edges)
    g.name("Bucky Ball")

    pos = {
        0: (1.00000000000000, 0.000000000000000),
        1: (-1.00000000000000, 0.000000000000000),
        2: (0.500000000000000, 0.866025403784439),
        3: (-0.500000000000000, 0.866025403784439),
        4: (-0.252886764483159, -0.146004241548845),
        5: (-0.368953972399043, 0.0928336233191176),
        6: (-0.217853192651371, -0.0480798425451855),
        7: (-0.255589950938772, 0.0495517623332213),
        8: (-0.390242139418333, -0.225306404242310),
        9: (-0.586398703939125, -0.0441575936410641),
        10: (-0.113926229169631, -0.101751920396670),
        11: (-0.0461308635969359, -0.0928422349110366),
        12: (-0.150564961379772, -0.164626477859040),
        13: (-0.0848818904865275, -0.246123271631605),
        14: (-0.170708060452244, 0.196571509298384),
        15: (-0.0672882312715990, 0.212706320404226),
        16: (-0.264873262319233, 0.273106701265196),
        17: (-0.254957754106411, 0.529914971178085),
        18: (-0.103469165775548, 0.00647061768205703),
        19: (-0.113590051906687, 0.0655812470455896),
        20: (-0.145082862532183, -0.0477870484199328),
        21: (-0.179962687765901, 0.103901506225732),
        22: (0.0573383021786124, 0.0863716172289798),
        23: (0.0311566333625530, 0.149538968816603),
        24: (-0.0573383021786121, 0.0863716172289799),
        25: (-0.0311566333625527, 0.149538968816603),
        26: (-0.0517345828877740, 0.00161765442051429),
        27: (-0.0244663616211774, -0.0456122902452611),
        28: (0.0517345828877743, 0.00161765442051431),
        29: (0.0244663616211777, -0.0456122902452611),
        30: (-0.0272682212665964, 0.0439946358247470),
        31: (0.0272682212665968, 0.0439946358247470),
        32: (0.179962687765901, 0.103901506225732),
        33: (0.145082862532184, -0.0477870484199329),
        34: (0.113590051906687, 0.0655812470455895),
        35: (0.103469165775548, 0.00647061768205698),
        36: (0.254957754106411, 0.529914971178085),
        37: (0.264873262319233, 0.273106701265196),
        38: (0.0672882312715993, 0.212706320404226),
        39: (0.170708060452245, 0.196571509298384),
        40: (1.59594559789866e-16, 0.450612808484620),
        41: (2.01227923213310e-16, 0.292008483097691),
        42: (0.0848818904865278, -0.246123271631605),
        43: (0.150564961379773, -0.164626477859040),
        44: (0.0461308635969362, -0.0928422349110366),
        45: (0.113926229169631, -0.101751920396670),
        46: (1.66533453693773e-16, -0.207803012451463),
        47: (1.80411241501588e-16, -0.131162494091179),
        48: (0.586398703939126, -0.0441575936410641),
        49: (0.390242139418333, -0.225306404242310),
        50: (0.255589950938772, 0.0495517623332212),
        51: (0.217853192651372, -0.0480798425451855),
        52: (0.368953972399044, 0.0928336233191175),
        53: (0.252886764483159, -0.146004241548845),
        54: (-0.104080710079810, -0.365940324584313),
        55: (0.104080710079811, -0.365940324584313),
        56: (-0.331440949832714, -0.485757377537020),
        57: (0.331440949832715, -0.485757377537021),
        58: (-0.500000000000000, -0.866025403784438),
        59: (0.500000000000000, -0.866025403784439)
    }

    g.set_pos(pos)

    return g


def GossetGraph():
    r"""
    Return the Gosset graph.

    The Gosset graph is the skeleton of the
    :meth:`~sage.geometry.polyhedron.library.Polytopes.Gosset_3_21` polytope. It
    has with 56 vertices and degree 27. For more information, see the
    :wikipedia:`Gosset_graph`.

    EXAMPLES::

        sage: g = graphs.GossetGraph(); g
        Gosset Graph: Graph on 56 vertices
        sage: g.order(), g.size()
        (56, 756)

    TESTS::

        sage: g.is_isomorphic(polytopes.Gosset_3_21().graph())  # not tested (~16s)
        True
    """
    string = ('w~~~~rt{~Z\\ZxnvYZYmlfrb}|hDuhLlcmmMNf_^zzQGNYcP\\kcRZbaJjoNBx{'
              '?N~o^}?A`}F_Kbbm_[QZ\\_]Cj\\oN_dm{BzB{?]WIMM@tPQRYBYRPIuAyJgQv?'
              '|Bxb_M[kWIR@jTQcciDjShXCkFMgpwqBKxeKoS`TYqdTCcKtkdKwWQXrbEZ@OdU'
              'mITZ@_e[{KXn?YPABzvY?IcO`zvYg@caC\\zlf?BaGR]zb{?@wOjv`~w??N_n_~'
              '~w???^_^~~{')

    G = Graph(string, name="Gosset Graph")

    ordering = [0, 2, 4, 6, 43, 23, 50, 18, 28, 9, 8, 7, 44, 3, 26, 35, 16, 14,
                33, 15, 54, 30, 17, 21, 10, 13, 36, 31, 55, 53, 51, 49, 12, 32,
                5, 37, 27, 46, 47, 48, 11, 52, 29, 20, 39, 41, 22, 40, 1, 25, 38,
                34, 45, 42, 19, 24]

    G._circle_embedding(ordering)

    return G


def DoubleStarSnark():
    r"""
    Return the double star snark.

    The double star snark is a 3-regular graph on 30 vertices. See the
    :wikipedia:`Double-star_snark`.

    EXAMPLES::

        sage: g = graphs.DoubleStarSnark()
        sage: g.order()
        30
        sage: g.size()
        45
        sage: g.chromatic_number()
        3
        sage: g.is_hamiltonian()                                                        # needs sage.numerical.mip
        False
        sage: g.automorphism_group().cardinality()                                      # needs sage.groups
        80
        sage: g.show()                                                                  # needs sage.plot
    """
    d = {0: [1, 14, 15],
         1: [0, 2, 11],
         2: [1, 3, 7],
         3: [2, 4, 18],
         4: [3, 5, 14],
         5: [10, 4, 6],
         6: [5, 21, 7],
         7: [8, 2, 6],
         8: [9, 13, 7],
         9: [24, 8, 10],
         10: [9, 11, 5],
         11: [1, 10, 12],
         12: [11, 27, 13],
         13: [8, 12, 14],
         14: [0, 4, 13],
         15: [0, 16, 29],
         16: [15, 20, 23],
         17: [25, 18, 28],
         18: [3, 17, 19],
         19: [18, 26, 23],
         20: [16, 28, 21],
         21: [20, 6, 22],
         22: [26, 21, 29],
         23: [16, 24, 19],
         24: [25, 9, 23],
         25: [24, 17, 29],
         26: [27, 19, 22],
         27: [12, 26, 28],
         28: [17, 27, 20],
         29: [25, 22, 15]}

    g = Graph(d, format='dict_of_lists', name="Double star snark")
    g._circle_embedding(list(range(15)), radius=2)
    g._circle_embedding(list(range(15, 30)), radius=1.4)

    return g


def MeredithGraph():
    r"""
    Return the Meredith Graph.

    The Meredith Graph is a 4-regular 4-connected non-hamiltonian graph. For
    more information on the Meredith Graph, see the :wikipedia:`Meredith_graph`.

    EXAMPLES::

        sage: g = graphs.MeredithGraph()
        sage: g.is_regular(4)
        True
        sage: g.order()
        70
        sage: g.size()
        140
        sage: g.radius()
        7
        sage: g.diameter()
        8
        sage: g.girth()
        4
        sage: g.chromatic_number()
        3
        sage: g.is_hamiltonian()                # long time                             # needs sage.numerical.mip
        False
    """
    g = Graph(name="Meredith Graph")
    g.add_vertex(0)

    # Edges between copies of K_{4,3}
    for i in range(5):
        g.add_edge(('outer', i, 3), ('outer', (i + 1) % 5, 0))
        g.add_edge(('inner', i, 3), ('inner', (i + 2) % 5, 0))
        g.add_edge(('outer', i, 1), ('inner', i, 1))
        g.add_edge(('outer', i, 2), ('inner', i, 2))

    # Edges inside of the K_{4,3}s.
    for i in range(5):
        for j in range(4):
            for k in range(3):
                g.add_edge(('inner', i, j), ('inner', i, k + 4))
                g.add_edge(('outer', i, j), ('outer', i, k + 4))

    g._circle_embedding(sum([[('outer', i, j) for j in range(4)] + 10 * [0] for i in range(5)], []),
                        radius=1, shift=2)
    g._circle_embedding(sum([[('outer', i, j) for j in range(4, 7)] + 10 * [0] for i in range(5)], []),
                        radius=1.2, shift=2.2)
    g._circle_embedding(sum([[('inner', i, j) for j in range(4)] + 7 * [0] for i in range(5)], []),
                        radius=.6, shift=1.24)
    g._circle_embedding(sum([[('inner', i, j) for j in range(4, 7)] + 5 * [0] for i in range(5)], []),
                        radius=.4, shift=1.05)

    g.delete_vertex(0)
    g.relabel()
    return g


def KittellGraph():
    r"""
    Return the Kittell Graph.

    For more information, see the `Wolfram page about the Kittel Graph
    <http://mathworld.wolfram.com/KittellGraph.html>`_.

    EXAMPLES::

        sage: g = graphs.KittellGraph()
        sage: g.order()
        23
        sage: g.size()
        63
        sage: g.radius()
        3
        sage: g.diameter()
        4
        sage: g.girth()
        3
        sage: g.chromatic_number()
        4
    """
    g = Graph({0: [1, 2, 4, 5, 6, 7], 1: [0, 2, 7, 10, 11, 13],
               2: [0, 1, 11, 4, 14], 3: [16, 12, 4, 5, 14], 4: [0, 2, 3, 5, 14],
               5: [0, 16, 3, 4, 6], 6: [0, 5, 7, 15, 16, 17, 18],
               7: [0, 1, 6, 8, 13, 18], 8: [9, 18, 19, 13, 7],
               9: [8, 10, 19, 20, 13], 10: [1, 9, 11, 13, 20, 21],
               11: [1, 2, 10, 12, 14, 15, 21], 12: [11, 16, 3, 14, 15],
               13: [8, 1, 10, 9, 7], 14: [11, 12, 2, 3, 4],
               15: [6, 11, 12, 16, 17, 21, 22],
               16: [3, 12, 5, 6, 15], 17: [18, 19, 22, 6, 15],
               18: [8, 17, 19, 6, 7], 19: [8, 9, 17, 18, 20, 22],
               20: [9, 10, 19, 21, 22], 21: [10, 11, 20, 22, 15],
               22: [17, 19, 20, 21, 15]},
              name="Kittell Graph")

    g._circle_embedding(list(range(3)), shift=.75)
    g._circle_embedding(list(range(3, 13)), radius=.4)
    g._circle_embedding(list(range(15, 22)), radius=.2, shift=-.15)
    pos = g._pos
    pos[13] = (-.65, -.35)
    pos[14] = (.65, -.35)
    pos[22] = (0, 0)

    return g


def CameronGraph():
    r"""
    Return the Cameron graph.

    The Cameron graph is strongly regular with parameters `v = 231, k = 30,
    \lambda = 9, \mu = 3`.

    For more information on the Cameron graph, see
    `<https://www.win.tue.nl/~aeb/graphs/Cameron.html>`_.

    EXAMPLES::

        sage: # needs sage.groups
        sage: g = graphs.CameronGraph()
        sage: g.order()
        231
        sage: g.size()
        3465
        sage: g.is_strongly_regular(parameters=True)    # long time
        (231, 30, 9, 3)
    """
    from sage.groups.perm_gps.permgroup_named import MathieuGroup
    from itertools import combinations
    g = Graph(name="Cameron Graph")
    sets = MathieuGroup(22).orbit((1, 2, 3, 7, 10, 20), action='OnSets')
    for s in sets:
        for a, b, c, d in combinations(set(s), 4):
            g.add_edges([((a, b), (c, d)), ((a, c), (b, d)), ((a, d), (b, c))])

    g.relabel()
    ordering = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 18, 19, 20,
                21, 24, 25, 26, 27, 29, 31, 34, 35, 38, 39, 96, 97, 101, 105,
                51, 117, 198, 32, 196, 201, 131, 167, 199, 197, 86, 102, 195,
                200, 186, 144, 202, 177, 44, 53, 58, 45, 48, 54, 43, 57, 50,
                46, 59, 133, 169, 104, 188, 118, 208, 157, 52, 207, 209, 132,
                204, 13, 187, 33, 203, 70, 145, 103, 168, 178, 87, 124, 123,
                125, 111, 120, 116, 119, 112, 95, 114, 115, 137, 218, 213, 108,
                76, 77, 74, 62, 64, 67, 63, 68, 69, 61, 41, 75, 73, 66, 71, 72,
                60, 22, 230, 151, 184, 138, 193, 109, 228, 174, 214, 219, 93,
                126, 143, 150, 146, 224, 181, 16, 223, 171, 90, 135, 106, 205,
                211, 121, 148, 160, 216, 222, 190, 36, 55, 185, 175, 94, 139,
                110, 215, 152, 220, 229, 194, 40, 128, 99, 141, 173, 154, 82,
                156, 164, 159, 28, 127, 158, 65, 162, 163, 153, 161, 155, 140,
                98, 47, 113, 84, 180, 30, 129, 179, 183, 165, 176, 142, 100,
                49, 134, 210, 170, 147, 91, 37, 206, 182, 191, 56, 136, 225,
                221, 149, 227, 217, 17, 107, 172, 212, 122, 226, 23, 85, 42,
                80, 92, 81, 89, 78, 83, 88, 79, 130, 192, 189, 166]

    g._circle_embedding(ordering)
    return g


def ChvatalGraph():
    r"""
    Return the Chvatal graph.

    Chvatal graph is one of the few known graphs to satisfy Grunbaum's
    conjecture that for every `m`, `n`, there is an `m`-regular, `m`-chromatic
    graph of girth at least `n`. For more information, see the
    :wikipedia:`Chv%C3%A1tal_graph`.

    EXAMPLES:

    The Chvatal graph has 12 vertices and 24 edges. It is a 4-regular,
    4-chromatic graph with radius 2, diameter 2, and girth 4::

        sage: G = graphs.ChvatalGraph(); G
        Chvatal graph: Graph on 12 vertices
        sage: G.order(); G.size()
        12
        24
        sage: G.degree()
        [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]
        sage: G.chromatic_number()
        4
        sage: G.radius(); G.diameter(); G.girth()
        2
        2
        4

    TESTS::

        sage: import networkx                                                           # needs networkx
        sage: G = graphs.ChvatalGraph()
        sage: G.is_isomorphic(Graph(networkx.chvatal_graph()))                          # needs networkx
        True
    """
    edges = {0: [1, 4, 6, 9], 1: [2, 5, 7], 2: [3, 6, 8], 3: [4, 7, 9],
             4: [5, 8], 5: [10, 11], 6: [10, 11], 7: [8, 11], 8: [10],
             9: [10, 11]}
    g = Graph(edges, format='dict_of_lists', name="Chvatal graph")
    g._circle_embedding(range(5), radius=4, angle=pi/2)
    g._circle_embedding(range(5, 10), radius=2, angle=pi/2)
    g._circle_embedding(range(10, 12), radius=1)
    return g


def ClebschGraph():
    r"""
    Return the Clebsch graph.

    See the :wikipedia:`Clebsch_graph` for more information.

    EXAMPLES::

        sage: g = graphs.ClebschGraph()
        sage: g.automorphism_group().cardinality()                                      # needs sage.groups
        1920
        sage: g.girth()
        4
        sage: g.chromatic_number()
        4
        sage: g.diameter()
        2
        sage: g.show(figsize=[10, 10])          # long time                             # needs sage.plot
    """
    g = Graph(pos={})
    x = 0
    for i in range(8):
        g.add_edge(x % 16, (x + 1) % 16)
        g.add_edge(x % 16, (x + 6) % 16)
        g.add_edge(x % 16, (x + 8) % 16)
        x += 1
        g.add_edge(x % 16, (x + 3) % 16)
        g.add_edge(x % 16, (x + 2) % 16)
        g.add_edge(x % 16, (x + 8) % 16)
        x += 1

    g._circle_embedding(list(range(16)), shift=.5)
    g.name("Clebsch graph")

    return g


def CoxeterGraph():
    r"""
    Return the Coxeter graph.

    See the :wikipedia:`Coxeter_graph`.

    EXAMPLES::

        sage: g = graphs.CoxeterGraph()
        sage: g.automorphism_group().cardinality()                                      # needs sage.groups
        336
        sage: g.girth()
        7
        sage: g.chromatic_number()
        3
        sage: g.diameter()
        4
        sage: g.show(figsize=[10, 10])          # long time                             # needs sage.plot
    """
    g = Graph({
            27: [6, 22, 14],
            24: [0, 7, 18],
            25: [8, 15, 2],
            26: [10, 16, 23],
            }, pos={})

    g.add_cycle(list(range(24)))
    g.add_edges([(5, 11), (9, 20), (12, 1), (13, 19), (17, 4), (3, 21)])

    g._circle_embedding(list(range(24)))
    g._circle_embedding([24, 25, 26], radius=.5)
    g._pos[27] = (0, 0)

    g.name("Coxeter Graph")

    return g


def CubeplexGraph(embedding='LM'):
    r"""
    Return the Cubeplex graph.

    The Cubeplex graph is the cubic hamiltonian graph of order 12 that
    corresponds to the graph labeled as `\Gamma_1` in Fischer and Little
    [FiLi2001]_. It has LCF notation `[-6, -5, -3, -6, 3, 5, -6, -3, 5, -6, -5,
    3]`.

    The Fischer-Little Theorem [FiLi2001]_ may be stated as follows [LM2024]_:

    A near-bipartite graph is non-Pfaffian if and only if it contains one of
    the graphs `K_{3, 3}`, `\Gamma_1` and `\Gamma_2` as an `S`-minor.

    Norine and Thomas [NT2007]_ use the term ``Cubeplex`` to describe one of
    the 12-vertex cubic graphs, `\Gamma_1` and `\Gamma_2`, as defined by
    Fischer and Little [FiLi2001]_. However, the figure in their paper that
    supposedly provides embeddings for the graphs labeled Cubeplex and Twinplex
    actually shows both embeddings corresponding to Fischer and Little's
    `\Gamma_1`, which is the Cubeplex graph. Followingly, for
    ``embedding='NT'``, we present only the embedding that is shown by the
    labeling ``Cubeplex`` in the paper of Norine and Thomas [NT2007]_.

    PLOTTING:

    Upon construction, the position dictionary is filled to override
    the spring-layout algorithm. For different values of the parameter
    ``embedding``, the Cubeplex graph is displayed as it is mentioned in the
    respective paper/ book.

    INPUT:

    - ``embedding`` -- string (default: ``'LM'``)

      - ``'LM'`` displays the embedding as shown for `\Gamma_1` by Lucchesi and
        Murty [LM2024]_

      - ``'FL'`` displays the embedding as shown for `\Gamma_1` by Fischer and
        Little [FiLi2001]_

      - ``'NT'`` displays the embedding as shown for the ``Cubeplex`` by Norine
        and Thomas [NT2007]_

    OUTPUT:

    - ``G`` -- the Cubeplex graph; note that a :class:`ValueError` is returned
      if ``embedding`` is none of ``'FT'``, ``'NT'`` or ``'LM'``

    EXAMPLES:

    Construct and show the Cubeplex graph::

        sage: g = graphs.CubeplexGraph()
        sage: g.name()
        'Cubeplex Graph'
        sage: g.order()
        12
        sage: g.size()
        18
        sage: g.girth()
        4
        sage: g.diameter()
        3
        sage: g.is_hamiltonian()
        True
        sage: g.crossing_number()
        1
        sage: g.show()                          # long time                             # needs sage.plot

    TESTS:

    Note that all three embeddings refer to the same graph, the Cubeplex graph,
    aka `\Gamma_1`::

        sage: fl = graphs.CubeplexGraph(embedding='FL')
        sage: nt = graphs.CubeplexGraph(embedding='NT')
        sage: lm = graphs.CubeplexGraph(embedding='LM')
        sage: fl.is_isomorphic(nt) and fl.is_isomorphic(lm)
        True

    The input parameter must be one of 'FL', 'NT' or 'LM'::

        sage: g = graphs.CubeplexGraph(embedding='embedding')
        Traceback (most recent call last):
        ...
        ValueError: parameter 'embedding' must be 'FL', 'NT' or 'LM'

    .. SEEALSO::

        :meth:`~sage.graphs.graph_generators.GraphGenerators.TwinplexGraph`

    AUTHORS:

    - Janmenjaya Panda (2024-08-03)
    """
    if embedding == 'FL':
        from math import pi

        G = Graph(12, name='Cubeplex Graph')
        G.add_cycle(list(range(12)))

        G.add_edges([
            (0, 3), (1, 6), (2, 8),
            (4, 9), (5, 11), (7, 10)
        ])

        G._circle_embedding(list(range(12)), angle=2*pi/3)

    elif embedding == 'NT':
        pos_dict = {
            0: (1, 2),
            1: (3, 2),
            2: (0, 1),
            3: (1, 1),
            4: (2, 1),
            5: (3, 1),
            6: (4, 1),
            7: (0, -1),
            8: (1, 0),
            9: (2, 0),
            10: (3, 0),
            11: (4, -1),
        }

        G = Graph(12, pos=pos_dict, name='Cubeplex Graph')
        G.add_edges([
            (0, 2), (0, 4), (0, 6),
            (1, 3), (1, 5), (1, 6),
            (2, 7), (2, 8), (3, 7),
            (3, 8), (4, 9), (4, 10),
            (5, 9), (5, 10), (6, 11),
            (7, 11), (8, 9), (10, 11)
        ])

    elif embedding == 'LM':
        from math import pi

        pos_dict = {
            8: (0, 1),
            9: (1, 0),
            10: (-3*cos(pi/16), -3*sin(pi/16)),
            11: (3*cos(pi/16), -3*sin(pi/16))
        }

        for v in range(8):
            t = pi * (v+2)/4
            pos_dict[v] = (-2*cos(t), 2*sin(t))

        G = Graph(12, pos=pos_dict, name='Cubeplex Graph')

        G.add_cycle(list(range(8)))
        G.add_edges([
            (0, 8), (1, 11), (2, 9), (3, 11), (4, 8),
            (5, 10), (6, 9), (7, 10), (8, 9), (10, 11)
        ])

    else:
        raise ValueError("parameter 'embedding' must be 'FL', 'NT' or 'LM'")

    return G


def DejterGraph():
    r"""
    Return the Dejter graph.

    The Dejter graph is obtained from the binary 7-cube by deleting a copy of
    the Hamming code of length 7. It is 6-regular, with 112 vertices and 336
    edges. For more information, see the :wikipedia:`Dejter_graph`.

    EXAMPLES::

        sage: g = graphs.DejterGraph(); g                                               # needs sage.rings.finite_rings
        Dejter Graph: Graph on 112 vertices
        sage: g.is_regular(k=6)                                                         # needs sage.rings.finite_rings
        True
        sage: g.girth()                                                                 # needs sage.rings.finite_rings
        4
    """
    from sage.graphs.generators.families import CubeGraph
    from sage.coding.hamming_code import HammingCode
    from sage.rings.finite_rings.finite_field_constructor import FiniteField

    g = CubeGraph(7)
    g.delete_vertices(["".join(map(str, x))
                       for x in HammingCode(FiniteField(2), 3)])
    g.name("Dejter Graph")
    return g


def DesarguesGraph():
    """
    Return the Desargues graph.

    PLOTTING: The layout chosen is the same as on the cover of [Har1994]_.

    EXAMPLES::

        sage: D = graphs.DesarguesGraph()
        sage: L = graphs.LCFGraph(20,[5,-5,9,-9],5)                                     # needs networkx
        sage: D.is_isomorphic(L)                                                        # needs networkx
        True
        sage: D.show()                          # long time                             # needs sage.plot
    """
    from sage.graphs.generators.families import GeneralizedPetersenGraph
    G = GeneralizedPetersenGraph(10, 3)
    G.name("Desargues Graph")
    return G


def DurerGraph():
    r"""
    Return the Dürer graph.

    For more information, see the :wikipedia:`D%C3%BCrer_graph`.

    EXAMPLES:

    The Dürer graph is named after Albrecht Dürer. It is a planar graph
    with 12 vertices and 18 edges::

        sage: G = graphs.DurerGraph(); G
        Durer graph: Graph on 12 vertices
        sage: G.is_planar()
        True
        sage: G.order()
        12
        sage: G.size()
        18

    The Dürer graph has chromatic number 3, diameter 4, and girth 3::

        sage: G.chromatic_number()
        3
        sage: G.diameter()
        4
        sage: G.girth()
        3

    Its automorphism group is isomorphic to `D_6`::

        sage: ag = G.automorphism_group()                                               # needs sage.groups
        sage: ag.is_isomorphic(DihedralGroup(6))                                        # needs sage.groups
        True
    """
    from sage.graphs.generators.families import GeneralizedPetersenGraph
    G = GeneralizedPetersenGraph(6, 2)
    G.name("Durer graph")
    return G


def DyckGraph():
    """
    Return the Dyck graph.

    For more information, see the `MathWorld article on the Dyck graph
    <http://mathworld.wolfram.com/DyckGraph.html>`_ or the
    :wikipedia:`Dyck_graph`.

    EXAMPLES:

    The Dyck graph was defined by Walther von Dyck in 1881. It has `32` vertices
    and `48` edges, and is a cubic graph (regular of degree `3`)::

        sage: G = graphs.DyckGraph(); G
        Dyck graph: Graph on 32 vertices
        sage: G.order()
        32
        sage: G.size()
        48
        sage: G.is_regular()
        True
        sage: G.is_regular(3)
        True

    It is non-planar and Hamiltonian, as well as bipartite (making it a bicubic
    graph)::

        sage: G.is_planar()
        False
        sage: G.is_hamiltonian()                                                        # needs sage.numerical.mip
        True
        sage: G.is_bipartite()
        True

    It has radius `5`, diameter `5`, and girth `6`::

        sage: G.radius()
        5
        sage: G.diameter()
        5
        sage: G.girth()
        6

    Its chromatic number is `2` and its automorphism group is of order `192`::

        sage: G.chromatic_number()
        2
        sage: G.automorphism_group().cardinality()                                      # needs sage.groups
        192

    It is a non-integral graph as it has irrational eigenvalues::

        sage: G.characteristic_polynomial().factor()                                    # needs sage.libs.pari sage.modules
        (x - 3) * (x + 3) * (x - 1)^9 * (x + 1)^9 * (x^2 - 5)^6

    It is a toroidal graph, and its embedding on a torus is dual to an embedding
    of the Shrikhande graph (:meth:`ShrikhandeGraph
    <GraphGenerators.ShrikhandeGraph>`).
    """
    pos_dict = {}
    for i in range(8):
        pos_dict[i] = [float(cos((2*i) * pi/8)),
                       float(sin((2*i) * pi/8))]
        pos_dict[8 + i] = [0.75 * pos_dict[i][0],
                           0.75 * pos_dict[i][1]]
        pos_dict[16 + i] = [0.50 * pos_dict[i][0],
                            0.50 * pos_dict[i][1]]
        pos_dict[24 + i] = [0.25 * pos_dict[i][0],
                            0.25 * pos_dict[i][1]]

    edge_dict = {
        0O00: [0O07, 0O01,   0O10], 0O10: [0O00,   0O27, 0O21],
        0O01: [0O00, 0O02,   0O11], 0O11: [0O01,   0O20, 0O22],
        0O02: [0O01, 0O03,   0O12], 0O12: [0O02,   0O21, 0O23],
        0O03: [0O02, 0O04,   0O13], 0O13: [0O03,   0O22, 0O24],
        0O04: [0O03, 0O05,   0O14], 0O14: [0O04,   0O23, 0O25],
        0O05: [0O04, 0O06,   0O15], 0O15: [0O05,   0O24, 0O26],
        0O06: [0O05, 0O07,   0O16], 0O16: [0O06,   0O25, 0O27],
        0O07: [0O06, 0O00,   0O17], 0O17: [0O07,   0O26, 0O20],

        0O20: [0O17, 0O11,   0O30], 0O30: [0O20,   0O35, 0O33],
        0O21: [0O10, 0O12,   0O31], 0O31: [0O21,   0O36, 0O34],
        0O22: [0O11, 0O13,   0O32], 0O32: [0O22,   0O37, 0O35],
        0O23: [0O12, 0O14,   0O33], 0O33: [0O23,   0O30, 0O36],
        0O24: [0O13, 0O15,   0O34], 0O34: [0O24,   0O31, 0O37],
        0O25: [0O14, 0O16,   0O35], 0O35: [0O25,   0O32, 0O30],
        0O26: [0O15, 0O17,   0O36], 0O36: [0O26,   0O33, 0O31],
        0O27: [0O16, 0O10,   0O37], 0O37: [0O27,   0O34, 0O32],
    }

    return Graph(edge_dict, pos=pos_dict, name="Dyck graph")


def HortonGraph():
    r"""
    Return the Horton Graph.

    The Horton graph is a cubic 3-connected non-hamiltonian graph. For more
    information, see the :wikipedia:`Horton_graph`.

    EXAMPLES::

        sage: # needs networkx
        sage: g = graphs.HortonGraph()
        sage: g.order()
        96
        sage: g.size()
        144
        sage: g.radius()
        10
        sage: g.diameter()
        10
        sage: g.girth()
        6
        sage: g.automorphism_group().cardinality()
        96
        sage: g.chromatic_number()
        2
        sage: g.is_hamiltonian()                # not tested (veeeery long)             # needs sage.numerical.mip
        False
    """
    g = Graph(name="Horton Graph")

    # Each group of the 6 groups of vertices is based on the same 3-regular
    # graph.
    from sage.graphs.generators.families import LCFGraph
    lcf = LCFGraph(16, [5, -5], 8)
    lcf.delete_edge(15, 0)
    lcf.delete_edge(7, 8)

    for i in range(6):
        for u, v in lcf.edge_iterator(labels=False):
            g.add_edge((i, u), (i, v))

    # Modifying the groups and linking them together
    c0 = (-1, 0)
    c1 = (-1, 1)
    c2 = (-1, 2)
    for i in range(3):
        g.add_edge((2 * i, 0), (2 * i + 1, 7))
        g.add_edge((2 * i + 1, 8), (2 * i, 7))
        g.add_edge((2 * i, 15), (2 * i + 1, 0))
        g.add_edge((2 * i, 8), c1)
        g.add_edge((2 * i + 1, 14), c2)
        g.add_edge((2 * i + 1, 10), c0)

    # Embedding
    for i in range(6):
        g._circle_embedding([(i, j) for j in range(16)],
                            center=(cos(2 * i * pi / 6), sin(2 * i * pi / 6)),
                            radius=.3)

    for i in range(3):
        g.delete_vertex((2 * i + 1, 15))

    g._circle_embedding([c0, c1, c2], radius=.2, shift=-0.75)

    g.relabel()

    return g


def EllinghamHorton54Graph():
    r"""
    Return the Ellingham-Horton 54-graph.

    For more information, see the :wikipedia:`Ellingham-Horton_graph`.

    EXAMPLES:

    This graph is 3-regular::

        sage: g = graphs.EllinghamHorton54Graph()
        sage: g.is_regular(k=3)
        True

    It is 3-connected and bipartite::

        sage: g.vertex_connectivity()   # not tested - too long
        3
        sage: g.is_bipartite()
        True

    It is not Hamiltonian::

        sage: g.is_hamiltonian()                # not tested                            # needs sage.numerical.mip
        False

    ... and it has a nice drawing ::

        sage: g.show(figsize=[10, 10])  # not tested - too long

    TESTS::

        sage: g.show()                          # long time                             # needs sage.plot
    """
    edge_dict = {
        0: [1, 11, 15], 1: [2, 47], 2: [3, 13], 3: [4, 8], 4: [5, 15],
        5: [6, 10], 6: [7, 30], 7: [8, 12], 8: [9], 9: [10, 29], 10: [11],
        11: [12], 12: [13], 13: [14], 14: [48, 15], 16: [17, 21, 28],
        17: [24, 29], 18: [19, 23, 30], 19: [20, 31], 20: [32, 21], 21: [33],
        22: [23, 27, 28], 23: [29], 24: [25, 30], 25: [26, 31], 26: [32, 27],
        27: [33], 28: [31], 32: [52], 33: [53], 34: [35, 39, 46], 35: [42, 47],
        36: [48, 37, 41], 37: [49, 38], 38: [50, 39], 39: [51],
        40: [41, 45, 46], 41: [47], 42: [48, 43], 43: [49, 44], 44: [50, 45],
        45: [51], 46: [49], 50: [52], 51: [53], 52: [53]}

    g = Graph(data=edge_dict, format='dict_of_lists',
              name="Ellingham-Horton 54-graph")

    # The set of vertices on top is 0..15
    # Bottom left is 16..33
    # Bottom right is 34..51
    # The two other vertices are 52, 53

    # Top
    g._circle_embedding(list(range(16)), center=(0, .5), shift=.5, radius=.5)

    # Bottom-left
    g._circle_embedding(list(range(16, 22)), center=(-1.5, -1))
    g._circle_embedding(list(range(22, 28)), center=(-1.5, -1), radius=.5)
    g._circle_embedding(list(range(28, 34)), center=(-1.5, -1), radius=.7)

    # Bottom right
    g._circle_embedding(list(range(34, 40)), center=(1.5, -1))
    g._circle_embedding(list(range(40, 46)), center=(1.5, -1), radius=.5)
    g._circle_embedding(list(range(46, 52)), center=(1.5, -1), radius=.7)

    d = g._pos
    d[52] = (-.3, -2.5)
    d[53] = (.3, -2.5)
    d[31] = (-2.2, -.9)
    d[28] = (-.8, -.9)
    d[46] = (2.2, -.9)
    d[49] = (.8, -.9)

    return g


def EllinghamHorton78Graph():
    r"""
    Return the Ellingham-Horton 78-graph.

    For more information, see the :wikipedia:`Ellingham%E2%80%93Horton_graph`

    EXAMPLES:

    This graph is 3-regular::

        sage: g = graphs.EllinghamHorton78Graph()
        sage: g.is_regular(k=3)
        True

    It is 3-connected and bipartite::

        sage: g.vertex_connectivity()   # not tested (too long)
        3
        sage: g.is_bipartite()
        True

    It is not Hamiltonian::

        sage: g.is_hamiltonian()                # not tested                            # needs sage.numerical.mip
        False

    ... and it has a nice drawing ::

        sage: g.show(figsize=[10,10])   # not tested (too long)

    TESTS::

        sage: g.show(figsize=[10, 10])  # not tested (too long)
    """
    g = Graph({
            0: [1, 5, 60], 1: [2, 12], 2: [3, 7], 3: [4, 14], 4: [5, 9],
            5: [6], 6: [7, 11], 7: [15], 8: [9, 13, 22], 9: [10],
            10: [11, 72], 11: [12], 12: [13], 13: [14], 14: [72],
            15: [16, 20], 16: [17, 27], 17: [18, 22], 18: [19, 29],
            19: [20, 24], 20: [21], 21: [22, 26], 23: [24, 28, 72],
            24: [25], 25: [26, 71], 26: [27], 27: [28], 28: [29],
            29: [69], 30: [31, 35, 52], 31: [32, 42], 32: [33, 37],
            33: [34, 43], 34: [35, 39], 35: [36], 36: [41, 63],
            37: [65, 66], 38: [39, 59, 74], 39: [40], 40: [41, 44],
            41: [42], 42: [74], 43: [44, 74], 44: [45], 45: [46, 50],
            46: [47, 57], 47: [48, 52], 48: [49, 75], 49: [50, 54],
            50: [51], 51: [52, 56], 53: [54, 58, 73], 54: [55],
            55: [56, 59], 56: [57], 57: [58], 58: [75], 59: [75],
            60: [61, 64], 61: [62, 71], 62: [63, 77], 63: [67],
            64: [65, 69], 65: [77], 66: [70, 73], 67: [68, 73],
            68: [69, 76], 70: [71, 76], 76: [77]}, pos={})

    g._circle_embedding(list(range(15)), center=(-2.5, 1.5))
    g._circle_embedding(list(range(15, 30)), center=(-2.5, -1.5))
    g._circle_embedding([30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41,
                         42, 74, 43, 44], center=(2.5, 1.5))
    g._circle_embedding([45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56,
                         57, 58, 75, 59], center=(2.5, -1.5))

    d = g._pos

    d[76] = (-.2, -.1)
    d[77] = (.2, .1)
    d[38] = (2.2, .1)
    d[52] = (2.3, -.1)
    d[15] = (-2.1, -.1)
    d[72] = (-2.1, .1)

    g._line_embedding([60, 61, 62, 63], first=(-1, 2), last=(1, 2))
    g._line_embedding([64, 65, 37], first=(-.5, 1.5), last=(1.2, 1.5))
    g._line_embedding([66, 73, 67, 68, 69], first=(1.2, -2), last=(-.8, -2))
    g._line_embedding([66, 70, 71], first=(.7, -1.5), last=(-1, -1.5))

    g.name("Ellingham-Horton 78-graph")
    return g


def ErreraGraph():
    r"""
    Return the Errera graph.

    For more information, see the :wikipedia:`Errera_graph`.

    EXAMPLES:

    The Errera graph is named after Alfred Errera. It is a planar graph on 17
    vertices and having 45 edges::

        sage: G = graphs.ErreraGraph(); G
        Errera graph: Graph on 17 vertices
        sage: G.is_planar()
        True
        sage: G.order()
        17
        sage: G.size()
        45

    The Errera graph is Hamiltonian with radius 3, diameter 4, girth 3, and
    chromatic number 4::

        sage: G.is_hamiltonian()                                                        # needs sage.numerical.mip
        True
        sage: G.radius()
        3
        sage: G.diameter()
        4
        sage: G.girth()
        3
        sage: G.chromatic_number()
        4

    Each vertex degree is either 5 or 6. That is, if `f` counts the number of
    vertices of degree 5 and `s` counts the number of vertices of degree 6, then
    `f + s` is equal to the order of the Errera graph::

        sage: D = G.degree_sequence()
        sage: D.count(5) + D.count(6) == G.order()
        True

    The automorphism group of the Errera graph is isomorphic to the dihedral
    group of order 20::

        sage: ag = G.automorphism_group()                                               # needs sage.groups
        sage: ag.is_isomorphic(DihedralGroup(10))                                       # needs sage.groups
        True
    """
    edge_dict = {
        0: [1, 7, 14, 15, 16],
        1: [2, 9, 14, 15],
        2: [3, 8, 9, 10, 14],
        3: [4, 9, 10, 11],
        4: [5, 10, 11, 12],
        5: [6, 11, 12, 13],
        6: [7, 8, 12, 13, 16],
        7: [13, 15, 16],
        8: [10, 12, 14, 16],
        9: [11, 13, 15],
        10: [12],
        11: [13],
        13: [15],
        14: [16]}
    return Graph(edge_dict, name="Errera graph")


def F26AGraph():
    r"""
    Return the F26A graph.

    The F26A graph is a symmetric bipartite cubic graph with 26 vertices and 39
    edges. For more information, see the :wikipedia:`F26A_graph`.

    EXAMPLES::

        sage: # needs networkx
        sage: g = graphs.F26AGraph(); g
        F26A Graph: Graph on 26 vertices
        sage: g.order(), g.size()
        (26, 39)
        sage: g.automorphism_group().cardinality()
        78
        sage: g.girth()
        6
        sage: g.is_bipartite()
        True
        sage: g.characteristic_polynomial().factor()
        (x - 3) * (x + 3) * (x^4 - 5*x^2 + 3)^6
    """
    from sage.graphs.generators.families import LCFGraph
    g = LCFGraph(26, [7, -7], 13)
    g.name("F26A Graph")
    return g


def FlowerSnark():
    """
    Return a Flower Snark.

    A flower snark has 20 vertices. It is part of the class of biconnected cubic
    graphs with edge chromatic number = 4, known as snarks. (i.e.: the Petersen
    graph). All snarks are not Hamiltonian, non-planar and have Petersen graph
    graph minors. See the :wikipedia:`Flower_snark`.

    PLOTTING: Upon construction, the position dictionary is filled to override
    the spring-layout algorithm. By convention, the nodes are drawn 0-14 on the
    outer circle, and 15-19 in an inner pentagon.

    EXAMPLES: Inspect a flower snark::

        sage: F = graphs.FlowerSnark()
        sage: F
        Flower Snark: Graph on 20 vertices
        sage: F.graph6_string()
        'ShCGHC@?GGg@?@?Gp?K??C?CA?G?_G?Cc'

    Now show it::

        sage: F.show()                          # long time                             # needs sage.plot
    """
    d = {0: [1, 14, 15], 1: [2, 11], 2: [3, 7], 3: [2, 4, 16], 4: [5, 14],
         5: [6, 10], 6: [5, 7, 17], 8: [7, 9, 13], 9: [10, 18], 11: [10, 12],
         12: [13, 19], 13: [14], 15: [19], 16: [15, 17], 18: [17, 19]}
    g = Graph(d, format='dict_of_lists', name="Flower Snark")
    g._circle_embedding(range(15), radius=2.5, angle=pi/2)
    g._circle_embedding(range(15, 20), radius=1, angle=pi/2)
    return g


def FolkmanGraph():
    """
    Return the Folkman graph.

    See the :wikipedia:`Folkman_graph`.

    EXAMPLES::

        sage: # needs networkx
        sage: g = graphs.FolkmanGraph()
        sage: g.order()
        20
        sage: g.size()
        40
        sage: g.diameter()
        4
        sage: g.girth()
        4
        sage: g.charpoly().factor()
        (x - 4) * (x + 4) * x^10 * (x^2 - 6)^4
        sage: g.chromatic_number()
        2
        sage: g.is_eulerian()
        True
        sage: g.is_hamiltonian()                                                        # needs sage.numerical_mip
        True
        sage: g.is_vertex_transitive()
        False
        sage: g.is_bipartite()
        True
    """
    from sage.graphs.generators.families import LCFGraph
    g = LCFGraph(20, [5, -7, -7, 5], 5)
    g.name("Folkman Graph")
    return g


def FosterGraph():
    """
    Return the Foster graph.

    See the :wikipedia:`Foster_graph`.

    EXAMPLES::

        sage: # needs networkx
        sage: g = graphs.FosterGraph()
        sage: g.order()
        90
        sage: g.size()
        135
        sage: g.diameter()
        8
        sage: g.girth()
        10
        sage: g.automorphism_group().cardinality()
        4320
        sage: g.is_hamiltonian()                                                        # needs sage.numerical_mip
        True
    """
    from sage.graphs.generators.families import LCFGraph
    g = LCFGraph(90, [17, -9, 37, -37, 9, -17], 15)
    g.name("Foster Graph")
    return g


def FranklinGraph():
    r"""
    Return the Franklin graph.

    For more information, see the :wikipedia:`Franklin_graph`.

    EXAMPLES:

    The Franklin graph is named after Philip Franklin. It is a 3-regular graph
    on 12 vertices and having 18 edges::

        sage: G = graphs.FranklinGraph(); G
        Franklin graph: Graph on 12 vertices
        sage: G.is_regular(3)
        True
        sage: G.order()
        12
        sage: G.size()
        18

    The Franklin graph is a Hamiltonian, bipartite graph with radius 3, diameter
    3, and girth 4::

        sage: G.is_hamiltonian()                                                        # needs sage.numerical_mip
        True
        sage: G.is_bipartite()
        True
        sage: G.radius()
        3
        sage: G.diameter()
        3
        sage: G.girth()
        4

    It is a perfect, triangle-free graph having chromatic number 2::

        sage: G.is_perfect()
        True
        sage: G.is_triangle_free()
        True
        sage: G.chromatic_number()
        2
    """
    edge_dict = {
        0: [1, 5, 6],
        1: [2, 7],
        2: [3, 8],
        3: [4, 9],
        4: [5, 10],
        5: [11],
        6: [7, 9],
        7: [10],
        8: [9, 11],
        10: [11]}
    g = Graph(edge_dict, format='dict_of_lists', name="Franklin graph")
    g._circle_embedding(range(6), radius=2)
    g._circle_embedding(range(6, 12), radius=1)
    return g


def FruchtGraph():
    """
    Return a Frucht Graph.

    A Frucht graph has 12 nodes and 18 edges. It is the smallest cubic identity
    graph. It is planar and Hamiltonian. See the :wikipedia:`Frucht_graph`.

    PLOTTING: Upon construction, the position dictionary is filled to override
    the spring-layout algorithm. By convention, the first seven nodes are on the
    outer circle, with the next four on an inner circle and the last in the
    center.

    EXAMPLES::

        sage: FRUCHT = graphs.FruchtGraph()
        sage: FRUCHT
        Frucht graph: Graph on 12 vertices
        sage: FRUCHT.graph6_string()
        'KhCKM?_EGK?L'
        sage: (graphs.FruchtGraph()).show()     # long time                             # needs networkx

    TESTS::

        sage: import networkx                                                           # needs networkx
        sage: G = graphs.FruchtGraph()
        sage: G.is_isomorphic(Graph(networkx.frucht_graph()))                           # needs networkx
        True
    """
    edges = {0: [1, 6, 7], 1: [2, 7], 2: [3, 8], 3: [4, 9], 4: [5, 9],
             5: [6, 10], 6: [10], 7: [11], 8: [9, 11], 10: [11]}
    g = Graph(edges, format='dict_of_lists', name="Frucht graph")
    g._circle_embedding(range(7), radius=2, angle=pi/2)
    g._circle_embedding(range(7, 11), radius=1, angle=pi/2)
    g._pos[11] = (0, 0)
    return g


def GoldnerHararyGraph():
    r"""
    Return the Goldner-Harary graph.

    For more information, see the :wikipedia:`Goldner%E2%80%93Harary_graph`.

    EXAMPLES:

    The Goldner-Harary graph is named after A. Goldner and Frank Harary. It is
    a planar graph having 11 vertices and 27 edges::

        sage: G = graphs.GoldnerHararyGraph(); G
        Goldner-Harary graph: Graph on 11 vertices
        sage: G.is_planar()
        True
        sage: G.order()
        11
        sage: G.size()
        27

    The Goldner-Harary graph is chordal with radius 2, diameter 2, and girth 3::

        sage: G.is_chordal()
        True
        sage: G.radius()
        2
        sage: G.diameter()
        2
        sage: G.girth()
        3

    Its chromatic number is 4 and its automorphism group is isomorphic to the
    dihedral group `D_6`::

        sage: G.chromatic_number()
        4
        sage: ag = G.automorphism_group()                                               # needs sage.groups
        sage: ag.is_isomorphic(DihedralGroup(6))                                        # needs sage.groups
        True
    """
    edge_dict = {
        0: [1, 3, 4],
        1: [2, 3, 4, 5, 6, 7, 10],
        2: [3, 7],
        3: [7, 8, 9, 10],
        4: [3, 5, 9, 10],
        5: [10],
        6: [7, 10],
        7: [8, 10],
        8: [10],
        9: [10]}

    pos = {
        0: (-2, 0),
        1: (0, 1.5),
        2: (2, 0),
        3: (0, -1.5),
        4: (-1.5, 0),
        5: (-0.5, 0.5),
        6: (0.5, 0.5),
        7: (1.5, 0),
        8: (0.5, -0.5),
        9: (-0.5, -0.5),
        10: (0, 0)}

    return Graph(edge_dict, pos=pos, name="Goldner-Harary graph")


def GolombGraph():
    r"""
    Return the Golomb graph.

    See the :wikipedia:`Golomb_graph` for more information.

    EXAMPLES:

    The Golomb graph is a planar and Hamiltonian graph with 10 vertices
    and 18 edges. It has chromatic number 4, diameter 3, radius 2 and
    girth 3. It can be drawn in the plane as a unit distance graph::

        sage: G = graphs.GolombGraph(); G                                               # needs sage.symbolic
        Golomb graph: Graph on 10 vertices
        sage: pos = G.get_pos()                                                         # needs sage.symbolic
        sage: def dist2(u, v):
        ....:     return (u[0]-v[0])**2 + (u[1]-v[1])**2
        sage: all(dist2(pos[u], pos[v]) == 1 for u, v in G.edge_iterator(labels=None))  # needs sage.symbolic
        True
    """
    edge_dict = {
        0: [1, 2, 3],
        1: [2, 5],
        2: [7],
        3: [4, 8, 9],
        4: [5, 9],
        5: [6, 9],
        6: [7, 9],
        7: [8, 9],
        8: [9]}
    pos_dict = {
        0: [QQ('1/6'), QQ('1/6') * sqrt(11)],
        1: [QQ('1/12') * sqrt(33) - QQ('1/12'), - sqrt(QQ('1/72') * sqrt(33) + QQ('7/72'))],
        2: [- QQ('1/12') * sqrt(33) - QQ('1/12'), - sqrt(- QQ('1/72') * sqrt(33) + QQ('7/72'))],
        3: [1, 0],
        4: [QQ('1/2'), - QQ('1/2') * sqrt(3)],
        5: [- QQ('1/2'), - QQ('1/2') * sqrt(3)],
        6: [-1, 0],
        7: [- QQ('1/2'), QQ('1/2') * sqrt(3)],
        8: [QQ('1/2'), QQ('1/2') * sqrt(3)],
        9: [0, 0]}
    return Graph(edge_dict, pos=pos_dict, name="Golomb graph")


def GrayGraph(embedding=1):
    r"""
    Return the Gray graph.

    See the :wikipedia:`Gray_graph`.

    INPUT:

    - ``embedding`` -- integer (default: `1`); two embeddings are available,
      and can be selected by setting ``embedding`` to 1 or 2

    EXAMPLES::

        sage: # needs networkx
        sage: g = graphs.GrayGraph()
        sage: g.order()
        54
        sage: g.size()
        81
        sage: g.girth()
        8
        sage: g.diameter()
        6
        sage: g.show(figsize=[10, 10])          # long time                             # needs sage.plot
        sage: graphs.GrayGraph(embedding=2).show(figsize=[10, 10])      # long time, needs sage.plot

    TESTS::

        sage: graphs.GrayGraph(embedding=3)                                             # needs networkx
        Traceback (most recent call last):
        ...
        ValueError: the value of embedding must be 1, 2, or 3
    """
    from sage.graphs.generators.families import LCFGraph
    g = LCFGraph(54, [-25, 7, -7, 13, -13, 25], 9)
    g.name("Gray graph")

    if embedding == 1:
        o = g.automorphism_group(orbits=True)[-1]
        g._circle_embedding(o[0], center=(0, 0), radius=1)
        g._circle_embedding(o[1], center=(0, 0), radius=.6, shift=-.5)

    elif embedding != 2:
        raise ValueError("the value of embedding must be 1, 2, or 3")

    return g


def GrotzschGraph():
    r"""
    Return the Grötzsch graph.

    The Grötzsch graph is an example of a triangle-free graph with chromatic
    number equal to 4. For more information, see the
    :wikipedia:`Gr%C3%B6tzsch_graph`.

    EXAMPLES:

    The Grötzsch graph is named after Herbert Grötzsch. It is a Hamiltonian
    graph with 11 vertices and 20 edges::

        sage: G = graphs.GrotzschGraph(); G
        Grotzsch graph: Graph on 11 vertices
        sage: G.is_hamiltonian()                                                        # needs sage.numerical.mip
        True
        sage: G.order()
        11
        sage: G.size()
        20

    The Grötzsch graph is triangle-free and having radius 2, diameter 2, and
    girth 4::

        sage: G.is_triangle_free()
        True
        sage: G.radius()
        2
        sage: G.diameter()
        2
        sage: G.girth()
        4

    Its chromatic number is 4 and its automorphism group is isomorphic to the
    dihedral group `D_5`::

        sage: G.chromatic_number()
        4
        sage: ag = G.automorphism_group()                                               # needs sage.groups
        sage: ag.is_isomorphic(DihedralGroup(5))                                        # needs sage.groups
        True
    """
    edges = [(0, u) for u in range(1, 6)]
    edges.append((10, 6))
    edges.append((10, 1))
    edges.append((6, 5))

    for u in range(6, 10):
        edges.append((u, u + 1))
        edges.append((u, u - 4))

    for u in range(7, 11):
        edges.append((u, u - 6))

    g = Graph(edges, format='list_of_edges', name="Grotzsch graph")
    g._circle_embedding(range(1, 6), radius=1, angle=pi/2)
    g._circle_embedding(range(6, 11), radius=2, angle=pi/2)
    g._pos[0] = (0, 0)

    return g


def HeawoodGraph():
    """
    Return a Heawood graph.

    The Heawood graph is a cage graph that has 14 nodes. It is a cubic symmetric
    graph. (See also the Möbius-Kantor graph, :meth:`~MobiusKantorGraph`). It is
    nonplanar and Hamiltonian. It has diameter 3, radius 3, girth 6, and
    chromatic number 2. It is 4-transitive but not 5-transitive.
    See the :wikipedia:`Heawood_graph`.

    PLOTTING: Upon construction, the position dictionary is filled to override
    the spring-layout algorithm. By convention, the nodes are positioned in a
    circular layout with the first node appearing at the top, and then
    continuing counterclockwise.

    EXAMPLES::

        sage: H = graphs.HeawoodGraph()
        sage: H
        Heawood graph: Graph on 14 vertices
        sage: H.graph6_string()
        'MhEGHC@AI?_PC@_G_'
        sage: (graphs.HeawoodGraph()).show()    # long time                             # needs sage.plot

    TESTS::

        sage: import networkx                                                           # needs networkx
        sage: G = graphs.HeawoodGraph()
        sage: G.is_isomorphic(Graph(networkx.heawood_graph()))                          # needs networkx
        True
    """
    edges = {0: [1, 5, 13], 1: [2, 10], 2: [3, 7], 3: [4, 12], 4: [5, 9],
             5: [6], 6: [7, 11], 7: [8], 8: [9, 13], 9: [10], 10: [11],
             11: [12], 12: [13]}
    g = Graph(edges, format='dict_of_lists', name="Heawood graph")
    g._circle_embedding(range(14), radius=1, angle=pi/2)
    return g


def HerschelGraph():
    r"""
    Return the Herschel graph.

    For more information, see the :wikipedia:`Herschel_graph`.

    EXAMPLES:

    The Herschel graph is named after Alexander Stewart Herschel. It is a
    planar, bipartite graph with 11 vertices and 18 edges::

        sage: G = graphs.HerschelGraph(); G
        Herschel graph: Graph on 11 vertices
        sage: G.is_planar()
        True
        sage: G.is_bipartite()
        True
        sage: G.order()
        11
        sage: G.size()
        18

    The Herschel graph is a perfect graph with radius 3, diameter 4, and girth
    4::

        sage: G.is_perfect()
        True
        sage: G.radius()
        3
        sage: G.diameter()
        4
        sage: G.girth()
        4

    Its chromatic number is 2 and its automorphism group is isomorphic to the
    dihedral group `D_6`::

        sage: G.chromatic_number()
        2
        sage: ag = G.automorphism_group()                                               # needs sage.groups
        sage: ag.is_isomorphic(DihedralGroup(6))                                        # needs sage.groups
        True
    """
    edge_dict = {
        0: [1, 3, 4],
        1: [2, 5, 6],
        2: [3, 7],
        3: [8, 9],
        4: [5, 9],
        5: [10],
        6: [7, 10],
        7: [8],
        8: [10],
        9: [10]}
    g = Graph(edge_dict, format='dict_of_lists', name="Herschel graph")
    g._circle_embedding(range(4), radius=2)
    g._circle_embedding(range(4, 10), radius=1)
    g._pos[10] = (0, 0)
    return g


def GritsenkoGraph():
    r"""
    Return SRG(65, 32, 15, 16) constructed by Gritsenko.

    We took the adjacency matrix from O. Gritsenko's [Gri2021]_ and extracted
    orbits of the automorphism group on the edges.

    EXAMPLES::

        sage: H = graphs.GritsenkoGraph(); H                                            # needs sage.groups
        Gritsenko strongly regular graph: Graph on 65 vertices
        sage: H.is_strongly_regular(parameters=True)                                    # needs sage.groups
        (65, 32, 15, 16)
    """
    from sage.groups.perm_gps.permgroup import PermutationGroup
    from functools import reduce
    p = '(0)(1,17,2,18)(3,6,4,5)(7,31,8,32)(9,25,10,26)(11,14,12,13)'
    p += '(15,24,16,23)(19,22,20,21)(27,29,28,30)(33,35,34,36)(37,61,38,62)'
    p += '(39,55,40,56)(41,43,42,44)(45,53,46,54)(47,63,48,64)(49,52,50,51)'
    p += '(57,59,58,60)'
    q = '(0)(1,20,32,6,9,27,23,13,2,19,31,5,10,28,24,14)'
    q += '(3,8,22,18,11,15,29,25,4,7,21,17,12,16,30,26)'
    q += '(33,37,52,63,41,46,60,55,34,38,51,64,42,45,59,56)'
    q += '(35,39,58,53,44,47,49,62,36,40,57,54,43,48,50,61)'
    a = PermutationGroup([p, q])
    edges = [(0, 1), (1, 2), (1, 6), (1, 7), (1, 9), (1, 11), (1, 14), (1, 21),
             (1, 24), (1, 36), (1, 38), (1, 40), (1, 42), (1, 44), (1, 47),
             (1, 48), (1, 50), (1, 52), (1, 54), (1, 55), (1, 56), (1, 58),
             (1, 62), (1, 63), (1, 64), (33, 35), (33, 38), (33, 46), (33, 47),
             (33, 49), (33, 51), (33, 57), (33, 61)]
    # use the union of the orbits of a on the edges
    return Graph(reduce(lambda x, y: x + y,
                        (a.orbit(o, action='OnSets') for o in edges)),
                 format='list_of_edges',
                 name="Gritsenko strongly regular graph")


def HigmanSimsGraph(relabel=True):
    r"""
    Return the Higman-Sims graph.

    The Higman-Sims graph is a remarkable strongly regular graph of degree 22 on
    100 vertices.  For example, it can be split into two sets of 50 vertices
    each, so that each half induces a subgraph isomorphic to the
    Hoffman-Singleton graph (:meth:`~HoffmanSingletonGraph`). This can be done
    in 352 ways (see `Higman-Sims graph
    <https://www.win.tue.nl/~aeb/graphs/Higman-Sims.html>`_ by Andries
    E. Brouwer, accessed 24 October 2009.)

    Its most famous property is that the automorphism group has an index 2
    subgroup which is one of the 26 sporadic groups [HS1968]_.

    The construction used here follows [Haf2004]_.

    See also the :wikipedia:`Higman–Sims_graph`.

    INPUT:

    - ``relabel`` -- boolean (default: ``True``); whether to relabel the
      vertices with consecutive integers. If ``False`` the labels are strings
      that are three digits long. "xyz" means the vertex is in group `x` (zero
      through three), pentagon or pentagram `y` (zero through four), and is
      vertex `z` (zero through four) of that pentagon or pentagram. See
      [Haf2004]_ for more.

    OUTPUT: the Higman-Sims graph

    EXAMPLES:

    A split into the first 50 and last 50 vertices will induce two copies of the
    Hoffman-Singleton graph, and we illustrate another such split, which is
    obvious based on the construction used::

        sage: H = graphs.HigmanSimsGraph()
        sage: A = H.subgraph(range(0,50))
        sage: B = H.subgraph(range(50,100))
        sage: K = graphs.HoffmanSingletonGraph()
        sage: K.is_isomorphic(A) and K.is_isomorphic(B)
        True
        sage: C = H.subgraph(range(25,75))
        sage: D = H.subgraph(list(range(0,25))+list(range(75,100)))
        sage: K.is_isomorphic(C) and K.is_isomorphic(D)
        True

    The automorphism group contains only one nontrivial proper normal subgroup,
    which is of index 2 and is simple.  It is known as the Higman-Sims group::

        sage: H = graphs.HigmanSimsGraph()
        sage: G = H.automorphism_group()                                                # needs sage.groups
        sage: g = G.order(); g                                                          # needs sage.groups
        88704000
        sage: K = G.normal_subgroups()[1]                                               # needs sage.groups
        sage: K.is_simple()                                                             # needs sage.groups
        True
        sage: g//K.order()                                                              # needs sage.groups
        2

    AUTHOR:

        - Rob Beezer (2009-10-24)
    """
    HS = Graph()
    HS.name('Higman-Sims graph')

    # Four groups of either five pentagons, or five pentagrams 4 x 5 x 5 = 100
    # vertices
    # First digit is "group", second is "penta{gon|gram}", third is "vertex"
    vlist = ['%d%d%d' % (g, p, v)
             for g in range(4) for p in range(5) for v in range(5)]
    HS.add_vertices(vlist)

    # Edges: Within groups 0 and 2, joined as pentagons
    # Edges: Within groups 1 and 3, joined as pentagrams
    for g in range(4):
        shift = 1
        if g in [1, 3]:
            shift += 1
        for p in range(5):
            for v in range(5):
                HS.add_edge(('%d%d%d' % (g, p, v), '%d%d%d' % (g, p, (v + shift) % 5)))

    # Edges: group 0 to group 1
    for x in range(5):
        for m in range(5):
            for c in range(5):
                y = (m * x + c) % 5
                HS.add_edge(('0%d%d' % (x, y), '1%d%d' % (m, c)))

    # Edges: group 1 to group 2
    for m in range(5):
        for A in range(5):
            for B in range(5):
                c = (2 * (m - A) * (m - A) + B) % 5
                HS.add_edge(('1%d%d' % (m, c), '2%d%d' % (A, B)))

    # Edges: group 2 to group 3
    for A in range(5):
        for a in range(5):
            for b in range(5):
                B = (2*A*A + 3*a*A - a*a + b) % 5
                HS.add_edge(('2%d%d' % (A, B), '3%d%d' % (a, b)))

    # Edges: group 3 to group 0
    for a in range(5):
        for b in range(5):
            for x in range(5):
                y = ((x - a) * (x - a) + b) % 5
                HS.add_edge(('3%d%d' % (a, b), '0%d%d' % (x, y)))

    # Edges: group 0 to group 2
    for x in range(5):
        for A in range(5):
            for B in range(5):
                y = (3*x*x + A*x + B + 1) % 5
                HS.add_edge(('0%d%d' % (x, y), '2%d%d' % (A, B)))
                y = (3*x*x + A*x + B - 1) % 5
                HS.add_edge(('0%d%d' % (x, y), '2%d%d' % (A, B)))

    # Edges: group 1 to group 3
    for m in range(5):
        for a in range(5):
            for b in range(5):
                c = (m*(m - a) + b + 2) % 5
                HS.add_edge(('1%d%d' % (m, c), '3%d%d' % (a, b)))
                c = (m*(m - a) + b - 2) % 5
                HS.add_edge(('1%d%d' % (m, c), '3%d%d' % (a, b)))

    # Layout vertices in a circle, in the order given in vlist
    HS._circle_embedding(vlist, radius=10, angle=pi/2)
    if relabel:
        HS.relabel(range(100))
    return HS


def HoffmanSingletonGraph():
    r"""
    Return the Hoffman-Singleton graph.

    The Hoffman-Singleton graph is the Moore graph of degree 7, diameter 2 and
    girth 5. The Hoffman-Singleton theorem states that any Moore graph with
    girth 5 must have degree 2, 3, 7 or 57. The first three respectively are the
    pentagon, the Petersen graph, and the Hoffman-Singleton graph. The existence
    of a Moore graph with girth 5 and degree 57 is still open.

    A Moore graph is a graph with diameter `d` and girth `2d + 1`. This implies
    that the graph is regular, and distance regular.

    For more details, see [GR2001]_ and the
    :wikipedia:`Hoffman–Singleton_graph`.

    PLOTTING: Upon construction, the position dictionary is filled to override
    the spring-layout algorithm. A novel algorithm written by Tom Boothby gives
    a random layout which is pleasing to the eye.

    EXAMPLES::

        sage: HS = graphs.HoffmanSingletonGraph()
        sage: Set(HS.degree())
        {7}
        sage: HS.girth()
        5
        sage: HS.diameter()
        2
        sage: HS.num_verts()
        50

    Note that you get a different layout each time you create the graph.  ::

        sage: HS.layout()[1]  # random
        (-0.844..., 0.535...)
        sage: HS = graphs.HoffmanSingletonGraph()
        sage: HS.layout()[1]  # random
        (-0.904..., 0.425...)
    """
    H = Graph({
        'q00': ['q01'], 'q01': ['q02'], 'q02': ['q03'], 'q03': ['q04'], 'q04': ['q00'],
        'q10': ['q11'], 'q11': ['q12'], 'q12': ['q13'], 'q13': ['q14'], 'q14': ['q10'],
        'q20': ['q21'], 'q21': ['q22'], 'q22': ['q23'], 'q23': ['q24'], 'q24': ['q20'],
        'q30': ['q31'], 'q31': ['q32'], 'q32': ['q33'], 'q33': ['q34'], 'q34': ['q30'],
        'q40': ['q41'], 'q41': ['q42'], 'q42': ['q43'], 'q43': ['q44'], 'q44': ['q40'],
        'p00': ['p02'], 'p02': ['p04'], 'p04': ['p01'], 'p01': ['p03'], 'p03': ['p00'],
        'p10': ['p12'], 'p12': ['p14'], 'p14': ['p11'], 'p11': ['p13'], 'p13': ['p10'],
        'p20': ['p22'], 'p22': ['p24'], 'p24': ['p21'], 'p21': ['p23'], 'p23': ['p20'],
        'p30': ['p32'], 'p32': ['p34'], 'p34': ['p31'], 'p31': ['p33'], 'p33': ['p30'],
        'p40': ['p42'], 'p42': ['p44'], 'p44': ['p41'], 'p41': ['p43'], 'p43': ['p40']})
    for j in range(5):
        for i in range(5):
            for k in range(5):
                con = (i + j*k) % 5
                H.add_edge(('q%d%d' % (k, con), 'p%d%d' % (j, i)))
    H.name('Hoffman-Singleton graph')
    from sage.combinat.permutation import Permutations
    from sage.misc.prandom import randint
    P = Permutations([1, 2, 3, 4])
    qpp = [0] + list(P[randint(0, 23)])
    ppp = [0] + list(P[randint(0, 23)])

    def qcycle(i, s):
        return ['q%s%s' % (i, (j + s) % 5) for j in qpp]

    def pcycle(i, s):
        return ['p%s%s' % (i, (j + s) % 5) for j in ppp]

    l = 0
    s = 0
    D = []
    while l < 5:
        for q in qcycle(l, s):
            D.append(q)
        vv = 'p%s' % q[1]
        s = int([v[-1] for v in H.neighbors(q) if v[:2] == vv][0])
        for p in pcycle(l, s):
            D.append(p)
        vv = 'q%s' % (int(p[1]) + 1)
        v = [v[-1] for v in H.neighbors(p) if v[:2] == vv]
        if len(v):
            s = int(v[0])
        l += 1
    map = H.relabel(range(50), return_map=True)
    H._circle_embedding([map[d] for d in D], angle=pi/2)
    return H


def HoffmanGraph():
    r"""
    Return the Hoffman Graph.

    See the :wikipedia:`Hoffman_graph`.

    EXAMPLES::

        sage: g = graphs.HoffmanGraph()
        sage: g.is_bipartite()
        True
        sage: g.is_hamiltonian()                # long time                             # needs sage.numerical.mip
        True
        sage: g.radius()
        3
        sage: g.diameter()
        4
        sage: g.automorphism_group().cardinality()                                      # needs sage.groups
        48
    """
    g = Graph({
            0: [1, 7, 8, 13],
            1: [2, 9, 14],
            2: [3, 8, 10],
            3: [4, 9, 15],
            4: [5, 10, 11],
            5: [6, 12, 14],
            6: [7, 11, 13],
            7: [12, 15],
            8: [12, 14],
            9: [11, 13],
            10: [12, 15],
            11: [14],
            13: [15]})

    g._circle_embedding(list(range(8)))
    g._circle_embedding(list(range(8, 14)), radius=.7, shift=.5)
    g._circle_embedding([14, 15], radius=.1)

    g.name("Hoffman Graph")

    return g


def HoltGraph():
    r"""
    Return the Holt graph (also called the Doyle graph).

    See the :wikipedia:`Holt_graph`.

    EXAMPLES::

        sage: g = graphs.HoltGraph();g
        Holt graph: Graph on 27 vertices
        sage: g.is_regular()
        True
        sage: g.is_vertex_transitive()                                                  # needs sage.groups
        True
        sage: g.chromatic_number()
        3
        sage: g.is_hamiltonian()                # long time                             # needs sage.numerical.mip
        True
        sage: g.radius()
        3
        sage: g.diameter()
        3
        sage: g.girth()
        5
        sage: g.automorphism_group().cardinality()                                      # needs sage.groups
        54
    """
    g = Graph(loops=False, name="Holt graph", pos={})
    for x in range(9):
        for y in range(3):
            g.add_edge((x, y), ((4 * x + 1) % 9, (y - 1) % 3))
            g.add_edge((x, y), ((4 * x - 1) % 9, (y - 1) % 3))
            g.add_edge((x, y), ((7 * x + 7) % 9, (y + 1) % 3))
            g.add_edge((x, y), ((7 * x - 7) % 9, (y + 1) % 3))

    for j in range(0, 6, 2):
        g._line_embedding([(x, j / 2) for x in range(9)],
                          first=(cos(2 * j * pi / 6), sin(2 * j * pi / 6)),
                          last=(cos(2 * (j + 1) * pi / 6), sin(2 * (j + 1) * pi / 6)))

    return g


def KrackhardtKiteGraph():
    """
    Return a Krackhardt kite graph with 10 nodes.

    The Krackhardt kite graph was originally developed by David Krackhardt for
    the purpose of studying social networks (see [Kre2002]_ and
    the :wikipedia:`Krackhardt_kite_graph`). It is used to show the distinction
    between degree centrality, betweenness centrality, and closeness
    centrality. For more information read the plotting section below in
    conjunction with the example.

    PLOTTING: Upon construction, the position dictionary is filled to override
    the spring-layout algorithm. By convention, the graph is drawn left to
    right, in top to bottom row sequence of [2, 3, 2, 1, 1, 1] nodes on each
    row. This places the fourth node (3) in the center of the kite, with the
    highest degree. But the fourth node only connects nodes that are otherwise
    connected, or those in its clique (i.e.: Degree Centrality). The eighth (7)
    node is where the kite meets the tail. It has degree = 3, less than the
    average, but is the only connection between the kite and tail (i.e.:
    Betweenness Centrality). The sixth and seventh nodes (5 and 6) are drawn in
    the third row and have degree = 5. These nodes have the shortest path to all
    other nodes in the graph (i.e.: Closeness Centrality).  Please execute the
    example for visualization.

    EXAMPLES:

    Construct and show a Krackhardt kite graph ::

        sage: g = graphs.KrackhardtKiteGraph()
        sage: g.show()                          # long time                             # needs sage.plot

    TESTS::

        sage: import networkx                                                           # needs networkx
        sage: G = graphs.KrackhardtKiteGraph()
        sage: G.is_isomorphic(Graph(networkx.krackhardt_kite_graph()))                  # needs networkx
        True
    """
    edges = {0: [1, 2, 3, 5], 1: [3, 4, 6], 2: [3, 5], 3: [4, 5, 6],
             4: [6], 5: [6, 7], 6: [7], 7: [8], 8: [9]}
    pos_dict = {0: (-1, 4), 1: (1, 4), 2: (-2, 3), 3: (0, 3), 4: (2, 3),
                5: (-1, 2), 6: (1, 2), 7: (0, 1), 8: (0, 0), 9: (0, -1)}
    return Graph(edges, pos=pos_dict, name="Krackhardt Kite Graph")


def Klein3RegularGraph():
    r"""
    Return the Klein 3-regular graph.

    The cubic Klein graph has 56 vertices and can be embedded on a
    surface of genus 3. It is the dual of
    :meth:`~sage.graphs.graph_generators.GraphGenerators.Klein7RegularGraph`.
    For more information, see the :wikipedia:`Klein_graphs`.

    EXAMPLES::

        sage: g = graphs.Klein3RegularGraph(); g
        Klein 3-regular Graph: Graph on 56 vertices
        sage: g.order(), g.size()
        (56, 84)
        sage: g.girth()
        7
        sage: g.automorphism_group().cardinality()                                      # needs sage.groups
        336
        sage: g.chromatic_number()
        3
    """
    g3 = Graph(':w`_GKWDBap`CMWFCpWsQUNdBwwuXPHrg`U`RIqypehVLqgHupYcFJyAv^Prk]'
               'EcarHwIVHAKh|\\tLVUxT]`ZDTJ{Af[o_AuKs{r_?ef',
               loops=False, multiedges=False)
    g3._circle_embedding([0, 2, 3, 4, 6, 8, 14, 1, 37, 30, 34, 48, 55, 43, 40,
                          45, 18, 20, 47, 42, 23, 17, 16, 10, 41, 11, 49, 25,
                          51, 26, 54, 9, 22, 15, 21, 12, 24, 7, 52, 31, 32, 36,
                          46, 35, 29, 50, 27, 19, 28, 5, 33, 13, 53, 39, 38, 44])
    g3.name("Klein 3-regular Graph")
    return g3


def Klein7RegularGraph():
    r"""
    Return the Klein 7-regular graph.

    The 7-valent Klein graph has 24 vertices and can be embedded on a surface of
    genus 3. It is the dual of
    :meth:`~sage.graphs.graph_generators.GraphGenerators.Klein3RegularGraph`.
    For more information, see the :wikipedia:`Klein_graphs`.

    EXAMPLES::

        sage: g = graphs.Klein7RegularGraph(); g
        Klein 7-regular Graph: Graph on 24 vertices
        sage: g.order(), g.size()
        (24, 84)
        sage: g.girth()
        3
        sage: g.automorphism_group().cardinality()                                      # needs sage.groups
        336
        sage: g.chromatic_number()
        4
    """
    g7 = Graph(':W__@`AaBbC_CDbDcE`F_AG_@DEH_IgHIJbFGIKaFHILeFGHMdFKN_EKOPaCNP'
               'Q`HOQRcGLRS`BKMSTdJKLPTU', loops=False, multiedges=False)
    g7._circle_embedding([0, 2, 3, 1, 9, 16, 20, 21, 4, 19, 17, 7, 15,
                          10, 8, 13, 11, 5, 23, 22, 14, 12, 18, 6])
    g7.name("Klein 7-regular Graph")
    return g7


def LocalMcLaughlinGraph():
    r"""
    Return the local McLaughlin graph.

    The local McLaughlin graph is a strongly regular graph with parameters
    `(162,56,10,24)`. It can be obtained from
    :meth:`~sage.graphs.graph_generators.GraphGenerators.McLaughlinGraph` by
    considering the stabilizer of a point: one of its orbits has cardinality
    162.

    EXAMPLES::

        sage: g = graphs.LocalMcLaughlinGraph(); g      # long time, optional - gap_package_design
        Local McLaughlin Graph: Graph on 162 vertices
        sage: g.is_strongly_regular(parameters=True)    # long time, optional - gap_package_design
        (162, 56, 10, 24)
    """
    g = McLaughlinGraph()
    orbits = g.automorphism_group().stabilizer(1).orbits()
    orbit = [x for x in orbits if len(x) == 162][0]
    g = g.subgraph(vertices=orbit)
    g.relabel()
    g.name("Local McLaughlin Graph")
    return g


def LjubljanaGraph(embedding=1):
    r"""
    Return the Ljubljana Graph.

    The Ljubljana graph is a bipartite 3-regular graph on 112 vertices and 168
    edges. It is not vertex-transitive as it has two orbits which are also
    independent sets of size 56. See the :wikipedia:`Ljubljana_graph`.

    The default embedding is obtained from the Heawood graph.

    INPUT:

    - ``embedding`` -- integer (default: `1`); two embeddings are available,
      and can be selected by setting ``embedding`` to 1 or 2

    EXAMPLES::

        sage: # needs networkx
        sage: g = graphs.LjubljanaGraph()
        sage: g.order()
        112
        sage: g.size()
        168
        sage: g.girth()
        10
        sage: g.diameter()
        8
        sage: g.show(figsize=[10, 10])          # long time                             # needs sage.plot
        sage: graphs.LjubljanaGraph(embedding=2).show(figsize=[10, 10])         # long time, needs sage.plot

    TESTS::

        sage: graphs.LjubljanaGraph(embedding=3)                                        # needs networkx
        Traceback (most recent call last):
        ...
        ValueError: the value of embedding must be 1 or 2
    """
    L = [47, -23, -31, 39, 25, -21, -31, -41, 25, 15, 29, -41, -19, 15,
         -49, 33, 39, -35, -21, 17, -33, 49, 41, 31, -15, -29, 41, 31,
         -15, -25, 21, 31, -51, -25, 23, 9, -17, 51, 35, -29, 21, -51,
         -39, 33, -9, -51, 51, -47, -33, 19, 51, -21, 29, 21, -31, -39]

    from sage.graphs.generators.families import LCFGraph
    g = LCFGraph(112, L, 2)
    g.name("Ljubljana graph")

    if embedding == 1:
        dh = HeawoodGraph().get_pos()

        # Correspondence between the vertices of the Heawood Graph and 8-sets of
        # the Ljubljana Graph.
        d = {
            0: [1, 21, 39, 57, 51, 77, 95, 107],
            1: [2, 22, 38, 58, 50, 78, 94, 106],
            2: [3, 23, 37, 59, 49, 79, 93, 105],
            3: [4, 24, 36, 60, 48, 80, 92, 104],
            4: [5, 25, 35, 61, 15, 81, 91, 71],
            9: [6, 26, 44, 62, 16, 82, 100, 72],
            10: [7, 27, 45, 63, 17, 83, 101, 73],
            11: [8, 28, 46, 64, 18, 84, 102, 74],
            12: [9, 29, 47, 65, 19, 85, 103, 75],
            13: [10, 30, 0, 66, 20, 86, 56, 76],
            8: [11, 31, 111, 67, 99, 87, 55, 43],
            7: [12, 32, 110, 68, 98, 88, 54, 42],
            6: [13, 33, 109, 69, 97, 89, 53, 41],
            5: [14, 34, 108, 70, 96, 90, 52, 40]
            }

        # The vertices of each 8-set are plotted on a circle, and the
        # circles are slowly shifted to obtain a symmetric drawing.
        for i, (u, vertices) in enumerate(d.items()):
            g._circle_embedding(vertices, center=dh[u], radius=.1,
                                shift=8.*i/14)

    elif embedding != 2:
        raise ValueError("the value of embedding must be 1 or 2")

    return g


def LivingstoneGraph():
    r"""
    Return the Livingstone Graph.

    The Livingstone graph is a distance-transitive graph on 266 vertices whose
    automorphism group is the :class:`J1 group
    <sage.groups.perm_gps.permgroup_named.JankoGroup>`. For more information,
    see the :wikipedia:`Livingstone_graph`.

    EXAMPLES::

        sage: # optional - internet
        sage: g = graphs.LivingstoneGraph()
        sage: g.order()
        266
        sage: g.size()
        1463
        sage: g.girth()
        5
        sage: g.is_vertex_transitive()
        True
        sage: g.is_distance_regular()
        True
    """
    from sage.groups.perm_gps.permgroup_named import JankoGroup
    from sage.graphs.graph import Graph
    G = JankoGroup(1)
    g = Graph(name="Livingstone Graph")
    g.add_edges(map(tuple, G.orbit((1, 24), action='OnSets')))
    return g


def M22Graph():
    r"""
    Return the M22 graph.

    The `M_{22}` graph is the unique strongly regular graph with parameters
    `v = 77, k = 16, \lambda = 0, \mu = 4`.

    For more information on the `M_{22}` graph, see
    `<https://www.win.tue.nl/~aeb/graphs/M22.html>`_.

    EXAMPLES::

        sage: # needs sage.groups
        sage: g = graphs.M22Graph()
        sage: g.order()
        77
        sage: g.size()
        616
        sage: g.is_strongly_regular(parameters=True)
        (77, 16, 0, 4)
    """
    from sage.groups.perm_gps.permgroup_named import MathieuGroup
    sets = [tuple(_) for _ in MathieuGroup(22).orbit((1, 2, 3, 7, 10, 20), action='OnSets')]
    g = Graph([sets, lambda x, y: not any(xx in y for xx in x)], name="M22 Graph")
    g.relabel()
    ordering = [0, 1, 3, 4, 5, 6, 7, 10, 12, 19, 20, 31, 2, 24, 35, 34, 22, 32,
                36, 23, 27, 25, 40, 26, 16, 71, 61, 63, 50, 68, 39, 52, 48, 44,
                69, 28, 9, 64, 60, 17, 38, 49, 45, 65, 14, 70, 72, 21, 43, 56,
                33, 73, 58, 55, 41, 29, 66, 54, 76, 46, 67, 11, 51, 47, 62, 53,
                15, 8, 18, 13, 59, 37, 30, 57, 75, 74, 42]

    g._circle_embedding(ordering)

    return g


def MarkstroemGraph():
    r"""
    Return the Markström Graph.

    The Markström Graph is a cubic planar graph with no cycles of length 4 nor
    8, but containing cycles of length 16. For more information, see the
    `Wolfram page about the Markström Graph
    <http://mathworld.wolfram.com/MarkstroemGraph.html>`_.

    EXAMPLES::

        sage: g = graphs.MarkstroemGraph()
        sage: g.order()
        24
        sage: g.size()
        36
        sage: g.is_planar()
        True
        sage: g.is_regular(3)
        True
        sage: g.subgraph_search(graphs.CycleGraph(4)) is None                           # needs sage.modules
        True
        sage: g.subgraph_search(graphs.CycleGraph(8)) is None                           # needs sage.modules
        True
        sage: g.subgraph_search(graphs.CycleGraph(16))                                  # needs sage.modules
        Subgraph of (Markstroem Graph): Graph on 16 vertices
    """
    g = Graph(name="Markstroem Graph")

    g.add_cycle(list(range(9)))
    g.add_path([0, 9, 10, 11, 2, 1, 11])
    g.add_path([3, 12, 13, 14, 5, 4, 14])
    g.add_path([6, 15, 16, 17, 8, 7, 17])
    g.add_cycle([10, 9, 18])
    g.add_cycle([12, 13, 19])
    g.add_cycle([15, 16, 20])
    g.add_cycle([21, 22, 23])
    g.add_edges([(19, 22), (18, 21), (20, 23)])

    g._circle_embedding(sum([[9 + 3*i + j for j in range(3)] + [0]*2 for i in range(3)], []),
                        radius=.6, shift=.7)
    g._circle_embedding([18, 19, 20], radius=.35, shift=.25)
    g._circle_embedding([21, 22, 23], radius=.15, shift=.25)
    g._circle_embedding(list(range(9)))

    return g


def McGeeGraph(embedding=2):
    r"""
    Return the McGee Graph.

    See the :wikipedia:`McGee_graph`.

    INPUT:

    - ``embedding`` -- integer (default: `2`); two embeddings are available,
      and can be selected by setting ``embedding`` to 1 or 2

    EXAMPLES::

        sage: # needs networkx
        sage: g = graphs.McGeeGraph()
        sage: g.order()
        24
        sage: g.size()
        36
        sage: g.girth()
        7
        sage: g.diameter()
        4
        sage: g.show()                                                                  # needs sage.plot
        sage: graphs.McGeeGraph(embedding=1).show()     # long time                     # needs sage.plot

    TESTS::

        sage: graphs.McGeeGraph(embedding=3)                                            # needs networkx
        Traceback (most recent call last):
        ...
        ValueError: the value of embedding must be 1 or 2
    """
    from sage.graphs.generators.families import LCFGraph
    g = LCFGraph(24, [12, 7, -7], 8)
    g.name('McGee graph')

    if embedding == 1:
        return g

    elif embedding == 2:

        o = [[7, 2, 13, 8, 19, 14, 1, 20],
             [5, 4, 11, 10, 17, 16, 23, 22],
             [3, 12, 9, 18, 15, 0, 21, 6]]

        g._circle_embedding(o[0], radius=1.5)
        g._circle_embedding(o[1], radius=3, shift=-.5)
        g._circle_embedding(o[2], radius=2.25, shift=.5)

        return g

    else:
        raise ValueError("the value of embedding must be 1 or 2")


def McLaughlinGraph():
    r"""
    Return the McLaughlin Graph.

    The McLaughlin Graph is the unique strongly regular graph of parameters
    `(275, 112, 30, 56)`.

    For more information on the McLaughlin Graph, see its web page on `Andries
    Brouwer's website <https://www.win.tue.nl/~aeb/graphs/McL.html>`_ which
    gives the definition that this method implements.

    .. NOTE::

        To create this graph you must have the gap_packages spkg installed.

    EXAMPLES::

        sage: g = graphs.McLaughlinGraph()              # optional - gap_package_design
        sage: g.is_strongly_regular(parameters=True)    # optional - gap_package_design
        (275, 112, 30, 56)
        sage: set(g.spectrum()) == {112, 2, -28}        # optional - gap_package_design
        True
    """
    from sage.combinat.designs.block_design import WittDesign
    from itertools import combinations
    from sage.sets.set import Set

    blocks = [Set(_) for _ in WittDesign(23).blocks()]

    B = [b for b in blocks if 0 in b]
    C = [b for b in blocks if 0 not in b]
    g = Graph()
    for b in B:
        for x in range(1, 23):
            if x not in b:
                g.add_edge(b, x)

    for b in C:
        for x in b:
            g.add_edge(b, x)

    for b, bb in combinations(B, 2):
        if len(b & bb) == 1:
            g.add_edge(b, bb)

    for c, cc in combinations(C, 2):
        if len(c & cc) == 1:
            g.add_edge(c, cc)

    for b in B:
        for c in C:
            if len(b & c) == 3:
                g.add_edge(b, c)

    # Here we relabel the elements of g in an architecture-independent way
    g.relabel({v: i for i, v in enumerate(list(range(1, 23)) +
                                          sorted(blocks, key=sorted))})
    g.name("McLaughlin")
    return g


def MoebiusKantorGraph():
    """
    Return a Möbius-Kantor Graph.

    A Möbius-Kantor graph is a cubic symmetric graph. (See also the Heawood
    graph). It has 16 nodes and 24 edges. It is nonplanar and Hamiltonian. It
    has diameter 4, girth 6, and chromatic number 2. It is identical to the
    Generalized Petersen graph, P[8, 3].

    For more details, see `Möbius-Kantor Graph - from Wolfram MathWorld
    <http://mathworld.wolfram.com/Moebius-KantorGraph.html>`_.

    PLOTTING: See the plotting section for the generalized Petersen graphs.

    EXAMPLES::

        sage: MK = graphs.MoebiusKantorGraph()
        sage: MK
        Moebius-Kantor Graph: Graph on 16 vertices
        sage: MK.graph6_string()
        'OhCGKE?O@?ACAC@I?Q_AS'
        sage: (graphs.MoebiusKantorGraph()).show()      # long time                     # needs sage.plot
    """
    from sage.graphs.generators.families import GeneralizedPetersenGraph
    G = GeneralizedPetersenGraph(8, 3)
    G.name("Moebius-Kantor Graph")
    return G


def MoserSpindle():
    r"""
    Return the Moser spindle.

    For more information, see the :wikipedia:`Moser_spindle`.

    EXAMPLES:

    The Moser spindle is a planar graph having 7 vertices and 11 edges::

        sage: # needs sage.symbolic
        sage: G = graphs.MoserSpindle(); G
        Moser spindle: Graph on 7 vertices
        sage: G.is_planar()
        True
        sage: G.order()
        7
        sage: G.size()
        11

    It is a Hamiltonian graph with radius 2, diameter 2, and girth 3::

        sage: # needs sage.symbolic
        sage: G.is_hamiltonian()                                                        # needs sage.numerical.mip
        True
        sage: G.radius()
        2
        sage: G.diameter()
        2
        sage: G.girth()
        3

    The Moser spindle can be drawn in the plane as a unit distance graph,
    has chromatic number 4, and its automorphism group is isomorphic to
    the dihedral group `D_4`::

        sage: # needs sage.symbolic
        sage: pos = G.get_pos()
        sage: all(sum((ui-vi)**2 for ui, vi in zip(pos[u], pos[v])) == 1
        ....:         for u, v in G.edge_iterator(labels=None))
        True
        sage: G.chromatic_number()
        4
        sage: ag = G.automorphism_group()
        sage: ag.is_isomorphic(DihedralGroup(4))
        True
    """
    edge_dict = {
        0: [1, 4, 6],
        1: [2, 5],
        2: [3, 5],
        3: [4, 5, 6],
        4: [6]}
    pos_dict = {
        0: [QQ('1/2'), 0],
        1: [- QQ('1/2'), 0],
        2: [- QQ('1/12') * sqrt(33) - QQ('1/4'),
            QQ('1/2') * sqrt(QQ('1/6') * sqrt(33) + QQ('17/6'))],
        3: [0, QQ('1/2') * sqrt(11)],
        4: [QQ('1/12') * sqrt(33) + QQ('1/4'),
            QQ('1/2') * sqrt(QQ('1/6') * sqrt(33) + QQ('17/6'))],
        5: [QQ('1/12') * sqrt(33) - QQ('1/4'),
            QQ('1/2') * sqrt(- QQ('1/6') * sqrt(33) + QQ('17/6'))],
        6: [- QQ('1/12') * sqrt(33) + QQ('1/4'),
            QQ('1/2') * sqrt(- QQ('1/6') * sqrt(33) + QQ('17/6'))]}
    return Graph(edge_dict, pos=pos_dict, name="Moser spindle")


def MurtyGraph():
    r"""
    Return the Murty graph.

    Consider the complete bipartite graph `K_{3, 3}`. There is a set of three
    black vertices and a set of three white vertices. Now, consider splicing
    the complete graph `K_4` with one of the black vertices, this generates the
    graph `K_4 \odot K_{3, 3}`. The Murty graph is obtained from
    `K_4 \odot K_{3, 3}` with the addition of an edge joining the remaining two
    black vertices. The Murty graph is free of conformal bicycles; in
    other words, the Murty graph is an example of a graph that is Birkhoff-von
    Neumann as well as PM-compact.

    This is the smallest brick that is Birkhoff-von Neumann, aka a solid
    brick, but is not odd-intercyclic. It is in this context that
    Prof. U.S.R. Murty first stumbled upon this graph, and it also appears in
    the work of Carvalho, Lucchesi, and Murty [CLM2006]_.

    PLOTTING:

    Upon construction, the position dictionary is filled to override
    the spring-layout algorithm. By convention, the Murty graph is
    displayed as mentioned in the paper [CKWL2019]_, with the first two
    (noncubic) vertices on the top row, the second three vertices (that form a
    stable set) in the middle row, and the remaining three vertices (that form
    a triangle) at the bottom.

    OUTPUT:

    - ``G`` -- the Murty graph

    EXAMPLES:

    Construct and show the Murty graph::

        sage: g = graphs.MurtyGraph()
        sage: g.name()
        'Murty Graph'
        sage: g.order()
        8
        sage: g.size()
        13
        sage: g.girth()
        3
        sage: g.diameter()
        2
        sage: g.is_hamiltonian()
        True
        sage: g.show()                          # long time                             # needs sage.plot

    REFERENCES:

    - [CKWL2019]_
    - [CLM2006]_
    - [LM2024]_

    AUTHORS:

    - Janmenjaya Panda (2024-08-03)
    """
    pos_dict = {
        0: (-0.5, sqrt(3)/2),
        1: (0.5, sqrt(3)/2),
        2: (-1, 0),
        3: (0, 0),
        4: (1, 0),
        5: (-0.5, -1 - sqrt(3)/2),
        6: (0, -1),
        7: (0.5, -1 - sqrt(3)/2)
    }

    G = Graph(8, pos=pos_dict, name="Murty Graph")

    G.add_edge(0, 1)
    for v in range(2, 5):
        G.add_edges([(0, v), (1, v), (v, v+3)])

    G.add_edges([(5, 6), (5, 7), (6, 7)])

    return G


def NauruGraph(embedding=2):
    """
    Return the Nauru Graph.

    See the :wikipedia:`Nauru_graph`.

    INPUT:

    - ``embedding`` -- integer (default: `2`); two embeddings are available,
      and can be selected by setting ``embedding`` to 1 or 2

    EXAMPLES::

        sage: g = graphs.NauruGraph()
        sage: g.order()
        24
        sage: g.size()
        36
        sage: g.girth()
        6
        sage: g.diameter()
        4
        sage: g.show()                                                                  # needs sage.plot
        sage: graphs.NauruGraph(embedding=1).show()     # long time                     # needs sage.plot

    TESTS::

        sage: graphs.NauruGraph(embedding=3)
        Traceback (most recent call last):
        ...
        ValueError: the value of embedding must be 1 or 2
        sage: graphs.NauruGraph(embedding=1).is_isomorphic(g)                           # needs networkx
        True
    """

    if embedding == 1:
        from sage.graphs.generators.families import LCFGraph
        g = LCFGraph(24, [5, -9, 7, -7, 9, -5], 4)
        g.name('Nauru Graph')
        return g
    elif embedding == 2:
        from sage.graphs.generators.families import GeneralizedPetersenGraph
        g = GeneralizedPetersenGraph(12, 5)
        g.name("Nauru Graph")
        return g
    else:
        raise ValueError("the value of embedding must be 1 or 2")


def PappusGraph():
    """
    Return the Pappus graph, a graph on 18 vertices.

    The Pappus graph is cubic, symmetric, and distance-regular.

    EXAMPLES::

        sage: G = graphs.PappusGraph()
        sage: G.show()                          # long time                             # needs sage.plot
        sage: L = graphs.LCFGraph(18, [5,7,-7,7,-7,-5], 3)                              # needs networkx
        sage: L.show()                          # long time                             # needs networkx sage.plot
        sage: G.is_isomorphic(L)                                                        # needs networkx
        True
    """
    edges = {0: [1, 5, 6], 1: [2, 7], 2: [3, 8], 3: [4, 9], 4: [5, 10], 5: [11],
             6: [13, 17], 7: [12, 14], 8: [13, 15], 9: [14, 16], 10: [15, 17],
             11: [12, 16], 12: [15], 13: [16], 14: [17]}
    g = Graph(edges, format='dict_of_lists', name="Pappus Graph")
    g._circle_embedding(range(6), radius=3, angle=pi/2)
    g._circle_embedding(range(6, 12), radius=2, angle=pi/2)
    g._circle_embedding(range(12, 18), radius=1, angle=pi/2)
    return g


def PoussinGraph():
    r"""
    Return the Poussin Graph.

    For more information on the Poussin Graph, see its corresponding `Wolfram
    page <http://mathworld.wolfram.com/PoussinGraph.html>`_.

    EXAMPLES::

        sage: g = graphs.PoussinGraph()
        sage: g.order()
        15
        sage: g.is_planar()
        True
    """
    g = Graph({2: [7, 8, 3, 4], 1: [7, 6], 0: [6, 5, 4], 3: [5]},
              name="Poussin Graph")

    g.add_cycle(list(range(3)))
    g.add_cycle(list(range(3, 9)))
    g.add_cycle(list(range(9, 14)))
    g.add_path([8, 12, 7, 11, 6, 10, 5, 9, 3, 13, 8, 12])
    g.add_edges([(14, i) for i in range(9, 14)])
    g._circle_embedding(list(range(3)), shift=.75)
    g._circle_embedding(list(range(3, 9)), radius=.4, shift=0)
    g._circle_embedding(list(range(9, 14)), radius=.2, shift=.4)
    g._pos[14] = (0,0)

    return g


def PetersenGraph():
    r"""
    Return the Petersen Graph.

    The Petersen Graph is a named graph that consists of 10 vertices and 15
    edges, usually drawn as a five-point star embedded in a pentagon.

    The Petersen Graph is a common counterexample. For example, it is not
    Hamiltonian.

    PLOTTING: See the plotting section for the generalized Petersen graphs.

    EXAMPLES: We compare below the Petersen graph with the default spring-layout
    versus a planned position dictionary of `(x, y)` tuples::

        sage: petersen_spring = Graph({0:[1,4,5], 1:[0,2,6], 2:[1,3,7],
        ....:                          3:[2,4,8], 4:[0,3,9], 5:[0,7,8],
        ....:                          6:[1,8,9], 7:[2,5,9], 8:[3,5,6],
        ....:                          9:[4,6,7]})
        sage: petersen_spring.show()            # long time                             # needs sage.plot
        sage: petersen_database = graphs.PetersenGraph()
        sage: petersen_database.show()          # long time                             # needs sage.plot
    """
    from sage.graphs.generators.families import GeneralizedPetersenGraph
    P = GeneralizedPetersenGraph(5, 2)
    P.name("Petersen graph")
    return P


def PerkelGraph():
    r"""
    Return the Perkel Graph.

    The Perkel Graph is a 6-regular graph with `57` vertices and `171` edges. It
    is the unique distance-regular graph with intersection array
    `(6,5,2;1,1,3)`. For more information, see the :wikipedia:`Perkel_graph` or
    https://www.win.tue.nl/~aeb/graphs/Perkel.html.

    EXAMPLES::

        sage: g = graphs.PerkelGraph(); g
        Perkel Graph: Graph on 57 vertices
        sage: g.is_distance_regular(parameters=True)
        ([6, 5, 2, None], [None, 1, 1, 3])
    """
    g = Graph(name="Perkel Graph")
    for i in range(19):
        g.add_edges(((0, i), (1, (i + j) % 19)) for j in [2, 5, 7])
        g.add_edges(((0, i), (2, (i + j) % 19)) for j in [5, -4, -8])
        g.add_edges(((1, i), (2, (i + j) % 19)) for j in [7, -4, -5])
    g.relabel()
    g._circle_embedding([0, 2, 3, 35, 8, 33, 45, 5, 53, 51, 18, 50, 29, 46, 30,
                         48, 40, 17, 20, 27, 43, 16, 7, 14, 6, 4, 15, 41, 24, 37,
                         28, 9, 55, 38, 19, 34, 39, 36, 54, 52, 44, 23, 12, 22,
                         32, 10, 13, 26, 1, 21, 42, 56, 49, 31, 47, 11, 25])
    return g


def RobertsonGraph():
    """
    Return the Robertson graph.

    See the :wikipedia:`Robertson_graph`.

    EXAMPLES::

        sage: # needs networkx
        sage: g = graphs.RobertsonGraph()
        sage: g.order()
        19
        sage: g.size()
        38
        sage: g.diameter()
        3
        sage: g.girth()
        5
        sage: g.charpoly().factor()
        (x - 4) * (x - 1)^2 * (x^2 + x - 5) * (x^2 + x - 1)
         * (x^2 - 3)^2 * (x^2 + x - 4)^2 * (x^2 + x - 3)^2
        sage: g.chromatic_number()
        3
        sage: g.is_hamiltonian()                                                        # needs sage.numerical.mip
        True
        sage: g.is_vertex_transitive()
        False
    """
    from sage.graphs.generators.families import LCFGraph
    lcf = [8, 4, 7, 4, 8, 5, 7, 4, 7, 8, 4, 5, 7, 8, 4, 8, 4, 8, 4]
    g = LCFGraph(19, lcf, 1)
    g.name("Robertson Graph")
    return g


def SchlaefliGraph():
    r"""
    Return the Schläfli graph.

    The Schläfli graph is the only strongly regular graphs of parameters
    `(27,16,10,8)` (see [GR2001]_).

    For more information, see the :wikipedia:`Schläfli_graph`.

    .. SEEALSO::

        :meth:`Graph.is_strongly_regular` -- tests whether a graph is strongly
        regular and/or returns its parameters.

    .. TODO::

        Find a beautiful layout for this beautiful graph.

    EXAMPLES:

    Checking that the method actually returns the Schläfli graph::

        sage: S = graphs.SchlaefliGraph()
        sage: S.is_strongly_regular(parameters=True)
        (27, 16, 10, 8)

    The graph is vertex-transitive::

        sage: S.is_vertex_transitive()                                                  # needs sage.groups
        True

    The neighborhood of each vertex is isomorphic to the complement of the
    Clebsch graph::

        sage: neighborhood = S.subgraph(vertices=S.neighbors(0))
        sage: graphs.ClebschGraph().complement().is_isomorphic(neighborhood)
        True
    """
    from sage.graphs.graph import Graph
    G = Graph('ZBXzr|}^z~TTitjLth|dmkrmsl|if}TmbJMhrJX]YfFyTbmsseztKTvyhDvw')
    order = [1, 8, 5, 10, 2, 6, 11, 15, 17, 13, 18, 12, 9, 24, 25, 3, 26, 7,
             16, 20, 23, 0, 21, 14, 22, 4, 19]
    G._circle_embedding(order)
    G.name("Schläfli graph")
    return G


def ShrikhandeGraph():
    """
    Return the Shrikhande graph.

    For more information, see the `MathWorld article on the Shrikhande graph
    <http://mathworld.wolfram.com/ShrikhandeGraph.html>`_ or the
    :wikipedia:`Shrikhande_graph`.

    .. SEEALSO::

        :meth:`Graph.is_strongly_regular` -- tests whether a graph is strongly
        regular and/or returns its parameters.

    EXAMPLES:

    The Shrikhande graph was defined by S. S. Shrikhande in 1959. It has `16`
    vertices and `48` edges, and is strongly regular of degree `6` with
    parameters `(2,2)`::

        sage: G = graphs.ShrikhandeGraph(); G
        Shrikhande graph: Graph on 16 vertices
        sage: G.order()
        16
        sage: G.size()
        48
        sage: G.is_regular(6)
        True
        sage: set([ len([x for x in G.neighbors(i) if x in G.neighbors(j)])
        ....:     for i in range(G.order())
        ....:     for j in range(i) ])
        {2}

    It is non-planar, and both Hamiltonian and Eulerian::

        sage: G.is_planar()
        False
        sage: G.is_hamiltonian()                                                        # needs sage.numerical.mip
        True
        sage: G.is_eulerian()
        True

    It has radius `2`, diameter `2`, and girth `3`::

        sage: G.radius()
        2
        sage: G.diameter()
        2
        sage: G.girth()
        3

    Its chromatic number is `4` and its automorphism group is of order `192`::

        sage: G.chromatic_number()
        4
        sage: G.automorphism_group().cardinality()                                      # needs sage.groups
        192

    It is an integral graph since it has only integral eigenvalues::

        sage: G.characteristic_polynomial().factor()                                    # needs sage.libs.pari sage.modules
        (x - 6) * (x - 2)^6 * (x + 2)^9

    It is a toroidal graph, and its embedding on a torus is dual to an
    embedding of the Dyck graph (:meth:`DyckGraph <GraphGenerators.DyckGraph>`).
    """
    pos_dict = {}
    for i in range(8):
        pos_dict[i] = [float(cos((2*i) * pi/8)),
                       float(sin((2*i) * pi/8))]
        pos_dict[8 + i] = [0.5 * pos_dict[i][0],
                           0.5 * pos_dict[i][1]]
    edge_dict = {
        0O00: [0O06, 0O07, 0O01, 0O02,   0O11, 0O17],
        0O01: [0O07, 0O00, 0O02, 0O03,   0O12, 0O10],
        0O02: [0O00, 0O01, 0O03, 0O04,   0O13, 0O11],
        0O03: [0O01, 0O02, 0O04, 0O05,   0O14, 0O12],
        0O04: [0O02, 0O03, 0O05, 0O06,   0O15, 0O13],
        0O05: [0O03, 0O04, 0O06, 0O07,   0O16, 0O14],
        0O06: [0O04, 0O05, 0O07, 0O00,   0O17, 0O15],
        0O07: [0O05, 0O06, 0O00, 0O01,   0O10, 0O16],

        0O10: [0O12, 0O13, 0O15, 0O16,   0O07, 0O01],
        0O11: [0O13, 0O14, 0O16, 0O17,   0O00, 0O02],
        0O12: [0O14, 0O15, 0O17, 0O10,   0O01, 0O03],
        0O13: [0O15, 0O16, 0O10, 0O11,   0O02, 0O04],
        0O14: [0O16, 0O17, 0O11, 0O12,   0O03, 0O05],
        0O15: [0O17, 0O10, 0O12, 0O13,   0O04, 0O06],
        0O16: [0O10, 0O11, 0O13, 0O14,   0O05, 0O07],
        0O17: [0O11, 0O12, 0O14, 0O15,   0O06, 0O00]
    }

    return Graph(edge_dict, pos=pos_dict, name="Shrikhande graph")


def SylvesterGraph():
    """
    Return the Sylvester Graph.

    This graph is obtained from the Hoffman Singleton graph by considering the
    graph induced by the vertices at distance two from the vertices of an (any)
    edge.

    For more information on the Sylvester graph, see
    `<https://www.win.tue.nl/~aeb/graphs/Sylvester.html>`_.

    .. SEEALSO::

        * :meth:`~sage.graphs.graph_generators.GraphGenerators.HoffmanSingletonGraph`.

    EXAMPLES::

        sage: g = graphs.SylvesterGraph(); g
        Sylvester Graph: Graph on 36 vertices
        sage: g.order()
        36
        sage: g.size()
        90
        sage: g.is_regular(k=5)
        True
    """
    g = HoffmanSingletonGraph()
    e = next(g.edge_iterator(labels=False))
    g.delete_vertices(g.neighbors(e[0]) + g.neighbors(e[1]))
    g.relabel()
    ordering = [0, 1, 2, 4, 5, 9, 16, 35, 15, 18, 20, 30, 22, 6, 33, 32, 14,
                10, 28, 29, 7, 24, 23, 26, 19, 12, 13, 21, 11, 31, 3, 27, 25,
                17, 8, 34]
    g._circle_embedding(ordering, shift=.5)
    g.name("Sylvester Graph")
    return g


def SimsGewirtzGraph():
    r"""
    Return the Sims-Gewirtz Graph.

    This graph is obtained from the Higman Sims graph by considering the graph
    induced by the vertices at distance two from the vertices of an (any)
    edge. It is the only strongly regular graph with parameters `v = 56`,
    `k = 10`, `\lambda = 0`, `\mu = 2`

    For more information on the Sylvester graph, see
    `<https://www.win.tue.nl/~aeb/graphs/Sims-Gewirtz.html>`_ or its
    :wikipedia:`Gewirtz_graph`.

    .. SEEALSO::

        * :meth:`~sage.graphs.graph_generators.GraphGenerators.HigmanSimsGraph`.

    EXAMPLES::

        sage: g = graphs.SimsGewirtzGraph(); g
        Sims-Gewirtz Graph: Graph on 56 vertices
        sage: g.order()
        56
        sage: g.size()
        280
        sage: g.is_strongly_regular(parameters = True)
        (56, 10, 0, 2)
    """
    g = HigmanSimsGraph()
    e = next(g.edge_iterator(labels=False))
    g.delete_vertices(g.neighbors(e[0]) + g.neighbors(e[1]))
    g.relabel()
    ordering = [0, 2, 3, 4, 6, 7, 8, 17, 1, 41, 49, 5, 22, 26, 11, 27, 15, 47,
                53, 52, 38, 43, 44, 18, 20, 32, 19, 42, 54, 36, 51, 30, 33, 35,
                37, 28, 34, 12, 29, 23, 55, 25, 40, 24, 9, 14, 48, 39, 45, 16,
                13, 21, 31, 50, 10, 46]
    g._circle_embedding(ordering)
    g.name("Sims-Gewirtz Graph")
    return g


def SousselierGraph():
    r"""
    Return the Sousselier Graph.

    The Sousselier graph is a hypohamiltonian graph on 16 vertices and 27
    edges. For more information, see :wikipedia:`Sousselier_graph` or
    the corresponding French
    `Wikipedia page <https://fr.wikipedia.org/wiki/Graphe_de_Sousselier>`_.

    EXAMPLES::

        sage: g = graphs.SousselierGraph()
        sage: g.order()
        16
        sage: g.size()
        27
        sage: g.radius()
        2
        sage: g.diameter()
        3
        sage: g.automorphism_group().cardinality()                                      # needs sage.groups
        2
        sage: g.is_hamiltonian()                                                        # needs sage.numerical.mip
        False
        sage: g.delete_vertex(g.random_vertex())
        sage: g.is_hamiltonian()                                                        # needs sage.numerical.mip
        True
    """
    g = Graph(name="Sousselier Graph")

    g.add_cycle(list(range(15)))
    g.add_path([12, 8, 3, 14])
    g.add_path([9, 5, 0, 11])
    g.add_edge(6, 2)
    g.add_edges([(15, i) for i in range(15) if i % 3 == 1])

    g._circle_embedding(list(range(15)), shift=-.25)
    g._pos[15] = (0, 0)

    return g


def SzekeresSnarkGraph():
    r"""
    Return the Szekeres Snark Graph.

    The Szekeres graph is a snark with 50 vertices and 75 edges. For more
    information on this graph, see the :wikipedia:`Szekeres_snark`.

    EXAMPLES::

        sage: g = graphs.SzekeresSnarkGraph()
        sage: g.order()
        50
        sage: g.size()
        75
        sage: g.chromatic_number()
        3
    """
    g = Graph(name="Szekeres Snark Graph")

    c = [(-1, i) for i in range(5)]
    for i in range(5):
        g.add_cycle([(i, j) for j in range(9)])
        g.delete_edge((i, 0), (i, 8))
        g.add_edge((i, 1), c[i])
        g.add_edge((i, 4), c[i])
        g.add_edge((i, 7), c[i])
        g.add_edge((i, 0), (i, 5))
        g.add_edge((i, 8), (i, 3))

        g.add_edge((i, 0), ((i + 1) % 5, 8))
        g.add_edge((i, 6), ((i + 2) % 5, 2))
        g._circle_embedding([(i, j) for j in range(9)],
                            radius=.3,
                            center=(cos(2 * (i + .25) * pi / 5),
                                    sin(2 * (i + .25) * pi / 5)),
                            shift=5.45 + 1.8 * i)

    g._circle_embedding(c, radius=1, shift=.25)

    g.relabel()
    return g


def ThomsenGraph():
    """
    Return the Thomsen Graph.

    The Thomsen Graph is actually a complete bipartite graph with `(n1, n2) =
    (3, 3)`. It is also called the Utility graph.

    PLOTTING: See CompleteBipartiteGraph.

    EXAMPLES::

        sage: T = graphs.ThomsenGraph()
        sage: T
        Thomsen graph: Graph on 6 vertices
        sage: T.graph6_string()
        'EFz_'
        sage: (graphs.ThomsenGraph()).show()    # long time                             # needs sage.plot
    """
    from sage.graphs.generators.basic import CompleteBipartiteGraph
    G = CompleteBipartiteGraph(3, 3)
    G.name("Thomsen graph")
    return G


def TietzeGraph():
    r"""
    Return the Tietze Graph.

    For more information on the Tietze Graph, see the
    :wikipedia:`Tietze%27s_graph`.

    EXAMPLES::

        sage: g = graphs.TietzeGraph()
        sage: g.order()
        12
        sage: g.size()
        18
        sage: g.diameter()
        3
        sage: g.girth()
        3
        sage: g.automorphism_group().cardinality()                                      # needs sage.groups
        12
        sage: g.automorphism_group().is_isomorphic(groups.permutation.Dihedral(6))      # needs sage.groups
        True
    """
    g = Graph([(0, 9), (3, 10), (6, 11), (1, 5), (2, 7), (4, 8)],
              name="Tietze Graph")
    g.add_cycle(list(range(9)))
    g.add_cycle([9, 10, 11])
    g._circle_embedding(list(range(9)))
    g._circle_embedding([9, 10, 11], radius=.5)
    return g


def TricornGraph():
    r"""
    Return the Tricorn graph.

    The Tricorn graph is obtained by splicing a complete graph `K_4` with the
    the triangular circular ladder graph `\overline{C_6}`. (Note that this
    generates a unqiue graph as both of the graphs `K_4` and `\overline{C_6}`
    are vertex-transitive). It is a nonsolid brick. This matching covered graph
    is one of the ten extremal cubic bricks. (A matching covered graph `G` is
    *extremal* if `\Phi(G) = dim(\mathcal{Lin}(G))`, where `\Phi(G)` denotes
    the number of perfect matchings of `G`, and `dim(\mathcal{Lin}(G))` stands
    for the dimension of the linear space of `G`).

    The Tricorn graph has no removable doubletons and has precisely three
    removable edges. The wheel graph `W_5` and the complete graph `K_4` are
    matching minors of the Tricorn graph.

    As per a theorem of Lovász [Lov1983]_, each non bipartite matching covered
    graph has a conformal subgraph which is either a bi-subdivision of `K_4` or
    of `\overline{C_6}` or both. In their paper, Kothari and Murty [KM2015]_
    characterized those planar bricks that are free of `\overline{C_6}` (that
    is, the planar bricks that do not contain a bi-subdivision of
    `\overline{C_6}` as a conformal subgraph). Besides two infinite families of
    matching covered graphs (odd wheel graphs and staircase graphs of order
    *4k*), the Tricorn graph is the only exception brick that is simple, planar
    and free of `\overline{C_6}`.

    PLOTTING:

    Upon construction, the position dictionary is filled to override
    the spring-layout algorithm. By convention, the Tricorn graph is
    displayed as mentioned in the book [LM2024]_, with the central vertex being
    the `0`-th one. Rest of the nine vertices are shown in groups of three,
    one on the top, rest two on the bottom left and on the bottom right
    corners respectively.

    OUTPUT:

    - ``G`` -- the Tricorn graph

    EXAMPLES:

    Construct and show the Tricorn graph; note that the edges `(2, 3)`,
    `(5, 6)` and `(8, 9)` are the only removable edges of the Tricorn
    graph::

        sage: g = graphs.TricornGraph()
        sage: g.name()
        'Tricorn Graph'
        sage: g.order()
        10
        sage: g.size()
        15
        sage: g.girth()
        3
        sage: g.diameter()
        3
        sage: g.is_hamiltonian()
        True
        sage: g.show()                          # long time                             # needs sage.plot

    REFERENCES:

    - [KM2015]_
    - [LM2024]_
    - [Lov1983]_

    AUTHORS:

    - Janmenjaya Panda (2024-08-02)
    """
    pos_dict = {
        0: (0, 0),
        1: (0, 1),
        2: (1/2, 1 + sqrt(3)/2),
        3: (-1/2, 1 + sqrt(3)/2),
        4: (-sqrt(3)/2, -1/2),
        5: (-sqrt(3)/2 - 1, -1/2),
        6: (-sqrt(3)/2 - 1/2, -1/2 - sqrt(3)/2),
        7: (sqrt(3)/2, -1/2),
        8: (sqrt(3)/2 + 1/2, -1/2 - sqrt(3)/2),
        9: (sqrt(3)/2 + 1, -1/2)
    }

    G = Graph(10, pos=pos_dict, name="Tricorn Graph")

    for v in range(1, 8, 3):
        G.add_edges([
            (0, v), (v, v+1),
            (v, v+2), (v+1, v+2),
            (v+2, int((-v**2 + 7*v + 4)/2))
        ])

    return G


def TruncatedIcosidodecahedralGraph():
    r"""
    Return the truncated icosidodecahedron.

    The truncated icosidodecahedron is an Archimedean solid with 30 square
    faces, 20 regular hexagonal faces, 12 regular decagonal faces, 120 vertices
    and 180 edges. For more information, see the
    :wikipedia:`Truncated_icosidodecahedron`.

    EXAMPLES:

    Unfortunately, this graph can not be constructed currently, due to numerical issues::

        sage: g = graphs.TruncatedIcosidodecahedralGraph(); g                           # needs sage.geometry.polyhedron sage.groups sage.rings.number_field
        Traceback (most recent call last):
        ...
        ValueError: *Error: Numerical inconsistency is found.  Use the GMP exact arithmetic.
        sage: g.order(), g.size()               # not tested                            # needs sage.geometry.polyhedron sage.groups sage.rings.number_field
        (120, 180)
    """
    from sage.geometry.polyhedron.library import polytopes
    # note that dropping exact=False here makes the construction take forever
    G = polytopes.icosidodecahedron(exact=False).truncation().graph()
    G.name("Truncated Icosidodecahedron")
    return G


def TruncatedTetrahedralGraph():
    r"""
    Return the truncated tetrahedron.

    The truncated tetrahedron is an Archimedean solid with 12 vertices and 18
    edges. For more information, see the :wikipedia:`Truncated_tetrahedron`.

    EXAMPLES::

        sage: g = graphs.TruncatedTetrahedralGraph(); g
        Truncated Tetrahedron: Graph on 12 vertices
        sage: g.order(), g.size()
        (12, 18)
        sage: g.is_isomorphic(polytopes.simplex(3).truncation().graph())                # needs sage.geometry.polyhedron
        True
    """
    g = Graph(':K`ESwC_EOyDl\\MCi', loops=False, multiedges=False)
    g._circle_embedding(list(range(6)), radius=1)
    g._circle_embedding(list(range(6, 9)), radius=.6, shift=.25)
    g._circle_embedding(list(range(9, 12)), radius=.2, shift=.25)
    g.name("Truncated Tetrahedron")
    return g


def Tutte12Cage():
    r"""
    Return the Tutte 12-Cage.

    See the :wikipedia:`Tutte_12-cage`.

    EXAMPLES::

        sage: # needs networkx
        sage: g = graphs.Tutte12Cage()
        sage: g.order()
        126
        sage: g.size()
        189
        sage: g.girth()
        12
        sage: g.diameter()
        6
        sage: g.show()                                                                  # needs sage.plot
    """
    L = [17, 27, -13, -59, -35, 35, -11, 13, -53, 53, -27, 21, 57, 11,
         -21, -57, 59, -17]

    from sage.graphs.generators.families import LCFGraph
    g = LCFGraph(126, L, 7)
    g.name("Tutte 12-Cage")
    return g


def TutteCoxeterGraph(embedding=2):
    r"""
    Return the Tutte-Coxeter graph.

    See the :wikipedia:`Tutte-Coxeter_graph`.

    INPUT:

    - ``embedding`` -- integer (default: `2`); two embeddings are available,
      and can be selected by setting ``embedding`` to 1 or 2

    EXAMPLES::

        sage: # needs networkx
        sage: g = graphs.TutteCoxeterGraph()
        sage: g.order()
        30
        sage: g.size()
        45
        sage: g.girth()
        8
        sage: g.diameter()
        4
        sage: g.show()                                                                  # needs sage.plot
        sage: graphs.TutteCoxeterGraph(embedding=1).show()      # long time             # needs sage.plot

    TESTS::

        sage: graphs.TutteCoxeterGraph(embedding=3)                                     # needs networkx
        Traceback (most recent call last):
        ...
        ValueError: the value of embedding must be 1 or 2
    """
    from sage.graphs.generators.families import LCFGraph
    g = LCFGraph(30, [-13, -9, 7, -7, 9, 13], 5)
    g.name("Tutte-Coxeter graph")

    if embedding == 1:
        d = {
            0: [1, 3, 5, 7, 29],
            1: [2, 4, 6, 28, 0],
            2: [8, 18, 26, 22, 12],
            3: [9, 13, 23, 27, 17],
            4: [11, 15, 21, 25, 19],
            5: [10, 14, 24, 20, 16]
            }

        g._circle_embedding(d[0], center=(-1, 1), radius=.25)
        g._circle_embedding(d[1], center=(1, 1), radius=.25)
        g._circle_embedding(d[2], center=(-.8, 0), radius=.25, shift=2.5)
        g._circle_embedding(d[3], center=(1.2, 0), radius=.25)
        g._circle_embedding(d[4], center=(-1, -1), radius=.25, shift=2)
        g._circle_embedding(d[5], center=(1, -1), radius=.25)

    elif embedding != 2:
        raise ValueError("the value of embedding must be 1 or 2")

    return g


def TutteGraph():
    r"""
    Return the Tutte Graph.

    The Tutte graph is a 3-regular, 3-connected, and planar non-hamiltonian
    graph. For more information on the Tutte Graph, see the
    :wikipedia:`Tutte_graph`.

    EXAMPLES::

        sage: g = graphs.TutteGraph()
        sage: g.order()
        46
        sage: g.size()
        69
        sage: g.is_planar()
        True
        sage: g.vertex_connectivity()           # long time                             # needs sage.numerical.mip
        3
        sage: g.girth()
        4
        sage: g.automorphism_group().cardinality()                                      # needs sage.groups
        3
        sage: g.is_hamiltonian()                                                        # needs sage.numerical.mip
        False
    """
    g = Graph(name="Tutte Graph")

    g.add_cycle([(i, j) for i in range(3) for j in range(3)])
    for i in range(3):
        g.add_cycle([(i, j) for j in range(9)])
        g.add_cycle([(i, j) for j in range(9, 14)])
        g.add_edge((i, 5), 0)
        g.add_edge((i, 13), (i, 3))
        g.add_edge((i, 12), (i, 1))
        g.add_edge((i, 11), (i, 8))
        g.add_edge((i, 10), (i, 7))
        g.add_edge((i, 6), (i, 14))
        g.add_edge((i, 4), (i, 14))
        g.add_edge((i, 9), (i, 14))

    g._circle_embedding([(i, j) for i in range(3) for j in range(6)], shift=.5)
    g._circle_embedding([(i, 14) for i in range(3)], radius=.3, shift=.25)

    for i in range(3):
        g._circle_embedding([(i, j) for j in range(3, 9)] + [0]*5,
                            shift=3.7*(i - 2) + .75,
                            radius=.4,
                            center=(.6*cos(2*(i + .25)*pi/3), .6*sin(2*(i + .25)*pi/3)))
        g._circle_embedding([(i, j) for j in range(9, 14)],
                            shift=1.7*(i - 2) + 1,
                            radius=.2,
                            center=(.6*cos(2*(i + .25)*pi/3), .6*sin(2*(i + .25)*pi/3)))

    g._pos[0] = (0,0)

    return g


def TwinplexGraph(embedding='LM'):
    r"""
    Return the Twinplex graph.

    The Twinplex graph is a cubic hamiltonian graph of order 12 with the graph
    crossing number 2 and has a girth 5 (that is the maximal girth among all
    cubic graphs on 12 vertices [CHNP2020]_). It corresponds to the graph
    labeled as `\Gamma_2` by Fischer and Little [FiLi2001]_. The Twinplex graph
    has LCF notation `[-5, -4, 4, -4, 4, 5, -4, 5, -4, 4, -5, 4]`.

    The Fischer-Little Theorem [FiLi2001]_ may be stated as follows [LM2024]_:

    A near-bipartite graph is non-Pfaffian if and only if it contains one of
    the graphs `K_{3, 3}`, `\Gamma_1` and `\Gamma_2` as an `S`-minor.

    Norine and Thomas [NT2007]_ use the term ``Twinplex`` to describe one of
    the 12-vertex cubic graphs, `\Gamma_1` and `\Gamma_2`, as defined by
    Fischer and Little [FiLi2001]_. However, the figure in their paper that
    supposedly provides embeddings for the graphs labeled Cubeplex and Twinplex
    actually shows both embeddings corresponding to Fischer and Little's
    `\Gamma_1`, which is the Cubeplex graph. Followingly, for
    ``embedding='NT'``, we present a correct version of the Twinplex graph
    with a slight modification of the embedding that is labeled as ``Twinplex``
    in the paper of Norine and Thomas [NT2007]_.

    PLOTTING:

    Upon construction, the position dictionary is filled to override
    the spring-layout algorithm. For different values of the parameter
    ``embedding``, the Twinplex graph is displayed as it is mentioned in the
    respective paper/ book. Note that for ``embedding='NT'``, a correct
    embedding of the Twinplex graph is displayed with a minor modification to
    the (incorrect) embedding shown in the paper [NT2007]_.

    INPUT:

    - ``embedding`` -- string (default: ``'LM'``)

      - ``'LM'`` displays the embedding as shown for `\Gamma_2` by Lucchesi and
        Murty [LM2024]_

      - ``'FL'`` displays the embedding as shown for `\Gamma_2` by Fischer and
        Little [FiLi2001]_

      - ``'NT'`` displays the correct embedding with a minor modification to
        the one shown as the (incorrect) ``Twinplex`` by Norine and Thomas
        [NT2007]_

      - ``'RST'`` displays the embedding as shown for the ``Twinplex`` by
        Robertson, Seymour and Thomas [RST2019]_

    OUTPUT:

    - ``G`` -- the Twinplex graph; note that a :class:`ValueError` is returned
      if ``embedding`` is none of ``'FT'``, ``'NT'``, ``'RST'`` or ``'LM'``

    EXAMPLES:

    Construct and show the Twinplex graph::

        sage: g = graphs.TwinplexGraph()
        sage: g.name()
        'Twinplex Graph'
        sage: g.order()
        12
        sage: g.size()
        18
        sage: g.girth()
        5
        sage: g.diameter()
        3
        sage: g.is_hamiltonian()
        True
        sage: g.crossing_number()
        2
        sage: g.show()                          # long time                             # needs sage.plot

    TESTS:

    Note that all four embeddings refer to the same graph, the Twinplex graph,
    aka `\Gamma_2`::

        sage: fl = graphs.TwinplexGraph(embedding='FL')
        sage: nt = graphs.TwinplexGraph(embedding='NT')
        sage: rst = graphs.TwinplexGraph(embedding='RST')
        sage: lm = graphs.TwinplexGraph(embedding='LM')
        sage: all(fl.is_isomorphic(g) for g in (nt, rst, lm))
        True

    The input parameter must be one of 'FL', 'NT', 'RST' or 'LM'::

        sage: g = graphs.TwinplexGraph(embedding='embedding')
        Traceback (most recent call last):
        ...
        ValueError: parameter 'embedding' must be 'FL', 'NT', 'LM' or 'RST'

    .. SEEALSO::

        :meth:`~sage.graphs.graph_generators.GraphGenerators.CubeplexGraph`

    AUTHORS:

    - Janmenjaya Panda (2024-08-03)
    """
    if embedding == 'FL':
        from math import pi

        G = Graph(12, name='Twinplex Graph')
        G.add_cycle(list(range(12)))

        G.add_edges([
            (0, 8), (1, 5), (2, 9),
            (3, 7), (4, 11), (6, 10)
        ])

        G._circle_embedding(list(range(12)), angle=5*pi/12)

    elif embedding == 'NT':
        pos_dict = {
            0: (1, 2),
            1: (3, 2),
            2: (0, 1),
            3: (1, 1),
            4: (2, 1),
            5: (3, 1),
            6: (4, 1),
            7: (0, -1),
            8: (1, 0),
            9: (2, 0),
            10: (3, 0),
            11: (4, -1),
        }

        G = Graph(12, pos=pos_dict, name='Twinplex Graph')
        G.add_edges([
            (0, 2), (0, 4), (0, 6),
            (1, 3), (1, 5), (1, 6),
            (2, 7), (2, 9), (3, 7),
            (3, 8), (4, 8), (4, 10),
            (5, 9), (5, 10), (6, 11),
            (7, 11), (8, 9), (10, 11)
        ])

    elif embedding == 'RST':
        pos_dict = {
            0: (-1, 3),
            1: (1, 3),
            2: (3, 1),
            3: (3, -1),
            4: (1, -3),
            5: (-1, -3),
            6: (-3, -1),
            7: (-3, 1),
            8: (-1, 1),
            9: (1, 1),
            10: (1, -1),
            11: (-1, -1)
        }

        G = Graph(12, pos=pos_dict, name='Twinplex Graph')

        G.add_cycle(list(range(8)))
        G.add_edges([
            (0, 4), (1, 8), (2, 10),
            (3, 9), (5, 10), (6, 8),
            (7, 11), (8, 9), (9, 11),
            (10, 11)
        ])

    elif embedding == 'LM':
        from math import pi

        pos_dict = {
            8: (0, 1),
            9: (1, 0),
            10: (-3*cos(pi/16), -3*sin(pi/16)),
            11: (3*cos(pi/16), -3*sin(pi/16))
        }

        for v in range(8):
            t = pi * (v+2)/4
            pos_dict[v] = (-2*cos(t), 2*sin(t))

        G = Graph(12, pos=pos_dict, name='Twinplex Graph')

        G.add_cycle(list(range(8)))
        G.add_edges([
            (0, 8), (1, 11), (2, 9), (3, 10), (4, 8),
            (5, 11), (6, 9), (7, 10), (8, 9), (10, 11)
        ])

    else:
        raise ValueError("parameter 'embedding' must be 'FL', 'NT',"
                         " 'LM' or 'RST'")

    return G


def WagnerGraph():
    """
    Return the Wagner Graph.

    See the :wikipedia:`Wagner_graph`.

    EXAMPLES::

        sage: # needs networkx
        sage: g = graphs.WagnerGraph()
        sage: g.order()
        8
        sage: g.size()
        12
        sage: g.girth()
        4
        sage: g.diameter()
        2
        sage: g.show()                                                                  # needs sage.plot
    """
    from sage.graphs.generators.families import LCFGraph
    g = LCFGraph(8, [4], 8)
    g.name("Wagner Graph")
    return g


def WatkinsSnarkGraph():
    r"""
    Return the Watkins Snark Graph.

    The Watkins Graph is a snark with 50 vertices and 75 edges. For more
    information, see the :wikipedia:`Watkins_snark`.

    EXAMPLES::

        sage: g = graphs.WatkinsSnarkGraph()
        sage: g.order()
        50
        sage: g.size()
        75
        sage: g.chromatic_number()
        3
    """
    g = Graph(name="Watkins Snark Graph")

    for i in range(5):
        g.add_cycle([(i, j) for j in range(9)])
        L = [(i, j) for j in range(4)] + [0, 0, (i, 4), 0, 0]
        L.extend((i, j) for j in range(5, 9))
        g._circle_embedding(L,
                            radius=.3,
                            center=(cos(2*(i + .25)*pi/5), sin(2*(i + .25)*pi/5)),
                            shift=2.7*i + 7.55)
        g.add_edge((i, 5), ((i + 1) % 5, 0))
        g.add_edge((i, 8), ((i + 2) % 5, 3))
        g.add_edge((i, 1), i)
        g.add_edge((i, 7), i)
        g.add_edge((i, 4), i)
        g.add_edge((i, 6), (i, 2))

    g._circle_embedding(list(range(5)), shift=.25, radius=1.1)
    return g


def WienerArayaGraph():
    r"""
    Return the Wiener-Araya Graph.

    The Wiener-Araya Graph is a planar hypohamiltonian graph on 42 vertices and
    67 edges. For more information, see the `Wolfram Page on the Wiener-Araya
    Graph <http://mathworld.wolfram.com/Wiener-ArayaGraph.html>`_ or
    :wikipedia:`Wiener-Araya_graph`.

    EXAMPLES::

        sage: g = graphs.WienerArayaGraph()
        sage: g.order()
        42
        sage: g.size()
        67
        sage: g.girth()
        4
        sage: g.is_planar()
        True
        sage: g.is_hamiltonian()                # not tested (30s)                      # needs sage.numerical.mip
        False
        sage: g.delete_vertex(g.random_vertex())
        sage: g.is_hamiltonian()                                                        # needs sage.numerical.mip
        True
    """
    g = Graph(name="Wiener-Araya Graph")

    g.add_cycle([(0, i) for i in range(4)])
    g.add_cycle([(1, i) for i in range(12)])
    g.add_cycle([(2, i) for i in range(20)])
    g.add_cycle([(3, i) for i in range(6)])
    g._circle_embedding([(0, i) for i in range(4)], shift=.5)
    g._circle_embedding(sum([[(1, 3 * i), (1, 3 * i + 1), 0, 0, 0, (1, 3 * i + 2), 0, 0, 0]
                             for i in range(4)], []),
                        shift=4,
                        radius=.65)
    g._circle_embedding([(2, i) for i in range(20)], radius=.5)
    g._circle_embedding([(3, i) for i in range(6)], radius=.3, shift=.5)

    for i in range(4):
        g.delete_edge((1, 3 * i), (1, 3 * i + 1))
        g.add_edge((1, 3 * i), (0, i))
        g.add_edge((1, 3 * i + 1), (0, i))
        g.add_edge((2, 5 * i + 2), (1, 3 * i))
        g.add_edge((2, 5 * i + 3), (1, 3 * i + 1))
        g.add_edge((2, (5 * i + 5) % 20), (1, 3 * i + 2))
        g.add_edge((2, (5 * i + 1) % 20), (3, i + (i >= 1) + (i >= 3)))
        g.add_edge((2, (5 * i + 4) % 20), (3, i + (i >= 1) + (i >= 3)))

    g.delete_edge((3, 1), (3, 0))
    g.add_edge((3, 1), (2, 4))
    g.delete_edge((3, 4), (3, 3))
    g.add_edge((3, 4), (2, 14))
    g.add_edge((3, 1), (3, 4))

    g.relabel()
    return g


def _EllipticLinesProjectivePlaneScheme(k):
    r"""
    Pseudo-cyclic association scheme for action of `O(3,2^k)` on elliptic lines.

    The group `O(3,2^k)` acts naturally on the `q(q-1)/2` lines of `PG(2,2^k)`
    skew to the conic preserved by it, see Sect. 12.7.B of [BCN1989]_ and
    Sect. 6.D in [BL1984]_. Compute the orbitals of this action and return them.

    This is a helper for
    :func:`sage.graphs.generators.smallgraphs.MathonStronglyRegularGraph`.

    INPUT:

    - ``k`` -- integer; the exponent of 2 to get the field size

    TESTS::

        sage: from sage.graphs.generators.smallgraphs import _EllipticLinesProjectivePlaneScheme
        sage: _EllipticLinesProjectivePlaneScheme(2)                                    # needs sage.libs.gap
        [
        [1 0 0 0 0 0]  [0 1 1 1 1 0]  [0 0 0 0 0 1]
        [0 1 0 0 0 0]  [1 0 1 1 0 1]  [0 0 0 0 1 0]
        [0 0 1 0 0 0]  [1 1 0 0 1 1]  [0 0 0 1 0 0]
        [0 0 0 1 0 0]  [1 1 0 0 1 1]  [0 0 1 0 0 0]
        [0 0 0 0 1 0]  [1 0 1 1 0 1]  [0 1 0 0 0 0]
        [0 0 0 0 0 1], [0 1 1 1 1 0], [1 0 0 0 0 0]
        ]
    """
    from sage.libs.gap.libgap import libgap
    from sage.matrix.constructor import matrix
    from itertools import product
    q = 2**k
    g0 = libgap.GeneralOrthogonalGroup(3, q)  # invariant form x0^2 + x1*x2
    g = libgap.Group(libgap.List(g0.GeneratorsOfGroup(), libgap.TransposedMat))
    W = libgap.FullRowSpace(libgap.GF(q), 3)
    l = sum(libgap.Elements(libgap.Basis(W)))
    gp = libgap.Action(g, libgap.Orbit(g, l, libgap.OnLines), libgap.OnLines)
    orbitals = gp.Orbits(list(product(gp.Orbit(1), gp.Orbit(1))),
                         libgap.OnTuples)
    mats = ([(int(x[0]) - 1, int(x[1]) - 1) for x in o] for o in orbitals)
    return [matrix((q * (q - 1)) // 2, lambda i, j: 1 if (i, j) in x else 0)
            for x in mats]


def MathonStronglyRegularGraph(t):
    r"""
    Return one of Mathon's graphs on 784 vertices.

    INPUT:

    - ``t`` -- integer; the number of the graph, from 0 to 2

    EXAMPLES::

        sage: # long time, needs sage.libs.gap
        sage: from sage.graphs.generators.smallgraphs import MathonStronglyRegularGraph
        sage: G = MathonStronglyRegularGraph(0)
        sage: G.is_strongly_regular(parameters=True)
        (784, 243, 82, 72)

    TESTS::

        sage: # long time, needs sage.libs.gap
        sage: G = graphs.MathonStronglyRegularGraph(1)
        sage: G.is_strongly_regular(parameters=True)
        (784, 270, 98, 90)
        sage: G = graphs.MathonStronglyRegularGraph(2)
        sage: G.is_strongly_regular(parameters=True)
        (784, 297, 116, 110)
    """
    from sage.graphs.generators.families import MathonPseudocyclicMergingGraph
    ES = _EllipticLinesProjectivePlaneScheme(3)
    return MathonPseudocyclicMergingGraph(ES, t)


def JankoKharaghaniGraph(v):
    r"""
    Return a `(936, 375, 150, 150)`-srg or a `(1800, 1029, 588, 588)`-srg.

    This functions returns a strongly regular graph for the two sets of
    parameters shown to be realizable in [JK2002]_. The paper also uses a
    construction from [GM1987]_.

    INPUT:

    - ``v`` -- integer; one of 936 or 1800

    EXAMPLES::

        sage: g = graphs.JankoKharaghaniGraph(936)      # long time                     # needs sage.libs.pari
        sage: g.is_strongly_regular(parameters=True)    # long time                     # needs sage.libs.pari
        (936, 375, 150, 150)

        sage: g = graphs.JankoKharaghaniGraph(1800)     # not tested (30s)
        sage: g.is_strongly_regular(parameters=True)    # not tested (30s)
        (1800, 1029, 588, 588)
    """
    from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
    from sage.matrix.constructor import matrix

    # The notations of [JK02] are rather tricky, and so this code attempts to
    # stick as much as possible to the paper's variable names.

    assert v in [1800, 936]

    J = matrix.ones
    I = matrix.identity

    # Definition of the 36x36 matrix H ([JK02], section 2)
    A = J(6)
    B = ("111---", "1---11", "1--1-1", "--111-", "-1-11-", "-11--1")
    C = ("-1-1-1", "1---11", "--11-1", "1-1-1-", "-1-11-", "111---")
    D = ("--1-11", "-11-1-", "11-1--", "--11-1", "11---1", "1--11-")
    E = ("-1--11", "1-1--1", "-11-1-", "---111", "1-11--", "11-1--")
    F = ("-1-1-1", "11--1-", "--111-", "1-11--", "-11--1", "1---11")
    B, C, D, E, F = (matrix([map({'1': 1, '-': -1}.get, r) for r in m])
                     for m in [B, C, D, E, F])

    H = [A, B, C, D, E, F]
    H = [[-x for x in H[6-i:]] + H[:6-i] for i in range(6)]
    H = matrix.block(H)

    # Definition of the BGW matrix W with the cyclotomic method
    # ([JK02] Lemma 1, and [GM87] Construction 1)
    m = 12
    t = (2 if v == 936 else 4)
    k = m
    q = m * t + 1
    K = GF(q, 'alpha')
    a = K.primitive_element()
    Ci = [[K(0)]] + [set(a**(k*j + i) for j in range(t)) for i in range(m)]
    Kelem_to_Ci = {v: i for i, s in enumerate(Ci) for v in s}  # maps v to [0,...,12]

    W = ([[0] + [1]*(len(K))] +
         [[1] + [Kelem_to_Ci[aj-ai] for aj in K] for ai in K])

    # The nonzero elements of W are considered as elements of C_12, generated by
    # a matrix Omega of order 12
    n = 18
    U = matrix.circulant([int(i == 1) for i in range(2 * n)])
    N = matrix.diagonal([1 if i else -1 for i in range(2 * n)])
    Omega = (U * N)**6
    assert Omega**12 == I(36)

    # The value w_{ij} is understood in the paper as matrix generated by Omega
    # acting on the left of a matrix L, which we now define.
    M = H - I(6).tensor_product(J(6))
    L = matrix(list(reversed(I(6).rows()))).tensor_product(I(6))

    # w_ij represents in the paper the matrix w_{ij}*L. We perform this action while
    # computing what is noted '[ M w_{ij} ]' in the paper.
    D = [[M*0 if w == 0 else M*(Omega**w)*L for w in R]
         for R in W]
    D = matrix.block(D)

    # for v=1800 the construction is slightly different, and we must add to D a
    # matrix which we now compute.
    if v == 1800:
        def my_abs(M):
            return matrix([[1 if x else 0 for x in R] for R in M.rows()])

        M = (J(6) + I(6)).tensor_product(J(6))  # we define M = (J(6)+I(6)) x J(6)
        D2 = [[M*0 if w == 0 else M*my_abs((Omega**w)*L) for w in R]  # '[ (J(6)+I(6)) x J(6) |w_{ij}| ]'
              for R in W]
        D = (D + matrix.block(D2))/2

    return Graph([e for e, v in D.dict().items() if v == 1],
                 multiedges=False,
                 name='Janko-Kharaghani')


def JankoKharaghaniTonchevGraph():
    r"""
    Return a `(324,153,72,72)`-strongly regular graph from [JKT2001]_.

    Build the graph using the description given in [JKT2001]_, taking sets B1
    and B163 in the text as adjacencies of vertices 1 and 163, respectively, and
    taking the edge orbits of the group `G` provided.

    EXAMPLES::

        sage: Gamma = graphs.JankoKharaghaniTonchevGraph()      # long time             # needs sage.libs.gap
        sage: Gamma.is_strongly_regular(parameters=True)        # long time             # needs sage.libs.gap
        (324, 153, 72, 72)
    """
    from sage.misc.misc_c import prod
    from sage.combinat.permutation import Permutation as P
    from sage.libs.gap.libgap import libgap

    m1 = prod(P((9 * x + k, 9 * x + k + 3, 9 * x + k + 6))
              for k in range(1, 4) for x in range(36))
    m2 = prod(P((3 * x + 1, 3 * x + 2, 3 * x + 3)) for x in range(108))
    t = prod(prod(map(P, [(9 * x + 2, 9 * x + 3), (9 * x + 4, 9 * x + 7),
                          (9 * x + 5, 9 * x + 9), (9 * x + 6, 9 * x + 8)]))
             for x in range(36))
    n1 = prod(prod(map(P, [(1 + x, 19 + x, 37 + x), (55 + x, 73 + x, 91 + x),
                           (109 + x, 127 + x, 145 + x), (163 + x, 181 + x, 199 + x),
                           (217 + x, 235 + x, 253 + x), (271 + x, 289 + x, 307 + x)]))
              for x in range(18))
    n2 = prod(prod(map(P, [(1 + x, 55 + x, 109 + x), (19 + x, 73 + x, 127 + x),
                           (37 + x, 91 + x, 145 + x), (163 + x, 217 + x, 271 + x),
                           (181 + x, 235 + x, 289 + x), (199 + x, 253 + x, 307 + x)]))
              for x in range(18))
    s = prod(prod(map(P, [(19 + x, 37 + x), (55 + x, 109 + x), (73 + x, 145 + x),
                          (91 + x, 127 + x), (181 + x, 199 + x), (217 + x, 271 + x),
                          (235 + x, 307 + x), (253 + x, 289 + x)]))
             for x in range(18))
    k = prod(prod(map(P, [(18 * x + 1, 18 * x + 10), (18 * x + 2, 18 * x + 11),
                          (18 * x + 3, 18 * x + 12), (18 * x + 4, 18 * x + 13),
                          (18 * x + 5, 18 * x + 14), (18 * x + 6, 18 * x + 15),
                          (18 * x + 7, 18 * x + 16), (18 * x + 8, 18 * x + 17),
                          (18 * x + 9, 18 * x + 18)]))
             for x in range(18))
    G = libgap.Group([libgap.PermList(p) for p in [m1, m2, t, n1, n2, s, k]])
    st = libgap.Group([libgap.PermList(p) for p in [t, s]])
    B1 = (19, 22, 25, 29, 30, 31, 33, 34, 35, 37, 40, 43, 47, 48, 49, 51, 52,
          53, 55, 56, 57, 65, 66, 67, 68, 70, 72, 76, 77, 78, 79, 80, 81, 82,
          86, 90, 92, 93, 95, 96, 98, 99, 100, 105, 107, 109, 110, 111, 119,
          120, 121, 122, 124, 126, 128, 129, 131, 132, 134, 135, 136, 141, 143,
          148, 149, 150, 151, 152, 153, 154, 158, 162, 167, 168, 170, 171, 172,
          176, 177, 179, 180, 184, 186, 187, 188, 190, 191, 192, 193, 196, 202,
          204, 205, 206, 208, 209, 210, 211, 214, 218, 219, 221, 225, 226, 227,
          228, 229, 232, 236, 237, 238, 241, 244, 245, 246, 249, 251, 254, 255,
          256, 259, 262, 265, 266, 268, 270, 272, 273, 275, 279, 280, 281, 282,
          283, 286, 290, 291, 292, 295, 298, 301, 302, 304, 306, 308, 309, 310,
          313, 316, 317, 318, 321, 323)
    B163 = (5, 6, 8, 9, 10, 14, 15, 17, 18, 22, 24, 25, 26, 28, 29, 30, 31, 34,
            40, 42, 43, 44, 46, 47, 48, 49, 52, 56, 57, 59, 63, 64, 65, 66, 67,
            70, 74, 75, 76, 79, 82, 83, 84, 87, 89, 92, 93, 94, 97, 100, 103,
            104, 106, 108, 110, 111, 113, 117, 118, 119, 120, 121, 124, 128,
            129, 130, 133, 136, 139, 140, 142, 144, 146, 147, 148, 151, 154,
            155, 156, 159, 161, 181, 185, 189, 191, 192, 194, 195, 197, 198,
            199, 203, 207, 209, 210, 212, 213, 215, 216, 217, 222, 224, 229,
            230, 231, 232, 233, 234, 236, 237, 238, 240, 241, 242, 244, 245,
            246, 254, 255, 256, 257, 259, 261, 262, 265, 268, 271, 276, 278,
            283, 284, 285, 286, 287, 288, 290, 291, 292, 293, 295, 297, 298,
            301, 304, 308, 309, 310, 312, 313, 314, 316, 317, 318)
    Gamma = Graph(multiedges=False, name='Janko-Kharaghani-Tonchev')
    for i, b in ((1, B1), (163, B163)):
        for j in (x[0] for x in st.OrbitsDomain(b)):
            Gamma.add_edges(map(tuple, G.Orbit(libgap.Set([i, j]), libgap.OnSets)))
    Gamma.relabel(range(Gamma.order()))
    return Gamma


def IoninKharaghani765Graph():
    r"""
    Return a `(765, 192, 48, 48)`-strongly regular graph.

    Existence of a strongly regular graph with these parameters was claimed in
    [IK2003]_.  Implementing the construction in the latter did not work,
    however. This function implements the following instructions, shared by Yury
    Ionin and Hadi Kharaghani.

        Let `A` be the affine plane over the field `GF(3)=\{-1,0,1\}`. Let

        .. MATH::

            \phi_1(x,y) &= x\\
            \phi_2(x,y) &= y\\
            \phi_3(x,y) &= x+y\\
            \phi_4(x,y) &= x-y\\

        For `i=1,2,3,4` and `j\in GF(3)`, let `L_{i,j}` be the line in `A`
        defined by `\phi_i(x,y)=j`. Let `\mathcal M` be the set of all 12 lines
        `L_{i,j}`, plus the empty set. Let `\pi` be the permutation defined on
        `\mathcal M` by `\pi(L_{i,j}) = L_{i,j+1}` and `\pi(\emptyset) =
        \emptyset`, so that `\pi` has three orbits of cardinality 3 and one of
        cardinality 1.

        Let `A=(p_1,...,p_9)` with `p_1=(-1,1)`, `p_2=(-1,0)`, `p_3=(-1,1)`,
        `p_4=(0,-1)`, `p_5=(0,0)`, `p_6=(0,1)`, `p_7=(1,-1)`, `p_8=(1,0)`,
        `p_9=(1,1)`. Note that `p_i+p_{10-i}=(0,0)`. For any subset `X` of `A`,
        let `M(X)` be the `(0,1)`-matrix of order 9 whose `(i,j)`-entry equals 1
        if and only if `p_{10-i}-p_j\in X`. Note that `M` is a symmetric matrix.

        An `MF`-tuple is an ordered quintuple `(X_1, X_2, X_3, X_4, X_5)` of
        subsets of `A`, of which one is the empty set and the other four are
        pairwise non-parallel lines. Such a quintuple generates the following
        block matrix:

        .. MATH::

            N(X_1, X_2, X_3, X_4, X_5) = \left( \begin{array}{ccccc}
                M(X_1) & M(X_2) & M(X_3) & M(X_4) & M(X_5)\\
                M(X_2) & M(X_3) & M(X_4) & M(X_5) & M(X_1)\\
                M(X_3) & M(X_4) & M(X_5) & M(X_1) & M(X_2)\\
                M(X_4) & M(X_5) & M(X_1) & M(X_2) & M(X_3)\\
                M(X_5) & M(X_1) & M(X_2) & M(X_3) & M(X_4)
                \end{array}\right)

        Observe that if `(X_1, X_2, X_3, X_4, X_5)` is an `MF`-tuple, then
        `N(X_1, X_2, X_3, X_4, X_5)` is the symmetric incidence matrix of a
        symmetric `(45, 12, 3)`-design.

        Let `\mathcal F` be the set of all `MF`-tuples and let `\sigma` be the
        following permutation of `\mathcal F`:

        .. MATH::

            \sigma(X_1, X_2, X_3, X_4, X_5) & = (X_2, X_3, X_4, X_5, X_1)\\
            \pi(X_1, X_2, X_3, X_4, X_5) & = (\pi(X_1), \pi(X_2), \pi(X_3), \pi(X_4), \pi(X_5))\\

        Observe that `\sigma` and `\pi` commute, and generate a (cyclic) group
        `G` of order 15. We will from now on identify `G` with the (cyclic)
        multiplicative group of the field `GF(16)` equal to
        `\{\omega^0,...,\omega^{14}\}`. Let `W=[w_{ij}]` be the following matrix
        of order 17 over `GF(16)=\{a_1,...,a_16\}`:

        .. MATH::

            w_{ij}=\left\{\begin{array}{ll}
                a_i+a_j & \text{if }1\leq i\leq 16, 1\leq j\leq 16,\\
                1       & \text{if }i=17, j\neq 17,\\
                1       & \text{if }i\neq 17, j= 17,\\
                0       & \text{if }i=j=17
                \end{array}\right.

        The diagonal entries of `W` are equal to 0, each off-diagonal entry can
        be represented as `\omega^k` with `0\leq k\leq 14`. Matrix `W` is a
        symmetric `BGW(17,16,15; G)`.

        Fix an `MF`-tuple `(X_1, X_2, X_3, X_4, X_5)` and let `S` be the block
        matrix obtained from `W` by replacing every diagonal entry of `W` by the
        zero matrix of order 45, and every off-diagonal entry `\omega^k` by the
        matrix `N(\sigma^k(X_1, X_2, X_3, X_4, X_5))` (through the association
        of `\omega^k` with an element of `G`). Then `S` is a symmetric incidence
        matrix of a symmetric `(765, 192, 48)`-design with zero diagonal, and
        therefore `S` is an adjacency matrix of a strongly regular graph with
        parameters `(765, 192, 48, 48)`.

    EXAMPLES::

        sage: g = graphs.IoninKharaghani765Graph(); g                                   # needs sage.modules sage.rings.finite_rings
        Ionin-Kharaghani: Graph on 765 vertices

    TESTS::

        sage: graphs.strongly_regular_graph(765, 192, 48, 48)                           # needs sage.modules sage.rings.finite_rings
        Ionin-Kharaghani: Graph on 765 vertices

    .. TODO::

        An update to [IK2003]_ meant to fix the problem encountered became available
        2016/02/24, see http://www.cs.uleth.ca/~hadi/research/IoninKharaghani.pdf
    """
    from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
    K = GF(3)

    # the four φ functions
    phi = [lambda xy: 1*xy[0] + 0*xy[1],
           lambda xy: 0*xy[0] + 1*xy[1],
           lambda xy: 1*xy[0] + 1*xy[1],
           lambda xy: 1*xy[0] - 1*xy[1]]

    # Defining L_{i,j}
    L = {(i, j): set() for i in range(4) for j in K}
    from itertools import product
    for p in product(K, K):
        for i in range(4):
            L[i, phi[i](p)].add(p)

    L = {k: frozenset(v) for k, v in L.items()}

    # Defining pi
    pi = {L[i, j]: L[i, (j + 1) % 3] for i, j in L}
    pi[frozenset()] = frozenset()

    # Defining A
    A = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 0), (0, 1), (1, -1), (1, 0), (1, 1)]

    def M(S):
        S = set((K(x), K(y)) for x, y in S)

        def difference(xy, xxyy):
            return (K(xy[0] - xxyy[0]), K(xy[1] - xxyy[1]))
        return matrix([[1 if difference(A[8-i], A[j]) in S else 0
                        for i in range(9)]
                       for j in range(9)])

    def N(Xi):
        Xi = [M(x) for x in Xi]
        return matrix.block([Xi[i:] + Xi[:i]
                             for i in range(len(Xi))])

    # sigma = lambda Xi: Xi[1:] + [pi[Xi[0]]]
    def f_pow(f, i, X):
        return f_pow(f, i-1, f(X)) if i else X

    def sigma2(Xi):
        return Xi[1:] + [Xi[0]]

    def pi_vec(x):
        return [pi.get(_) for _ in x]

    # The matrix W, with off-diagonal entries equal to integers 1,...,15
    # (instead of x^1,...,x^15)
    from sage.matrix.constructor import matrix
    GF16 = GF(16, 'x')
    W = matrix([[x + y for x in GF16] + [1] for y in GF16] +
               [[1]*16 + [0]])
    x = GF16.primitive_element()
    log_x = {x**i: i for i in range(15)}
    W = W.apply_map(lambda x: log_x[x] + 1 if x else 0)

    # Associate a matrix to every entry of W
    int_to_matrix = {0: matrix.zero(45)}
    for i in range(15):
        vec = [frozenset([]), L[0, 0], L[1, 0], L[2, 0], L[3, 0]]
        vec = f_pow(pi_vec, i % 3, vec)
        vec = f_pow(sigma2, i % 5, vec)
        int_to_matrix[i + 1] = N(vec)

    M2 = matrix.block([[int_to_matrix[x] for x in R] for R in W.rows()])
    g = Graph(M2, name='Ionin-Kharaghani')
    return g


def U42Graph216():
    r"""
    Return a (216,40,4,8)-strongly regular graph from [CRS2016]_.

    Build the graph, interpreting the `U_4(2)`-action considered in [CRS2016]_
    as the one on the hyperbolic lines of the corresponding unitary polar space,
    and then doing the unique merging of the orbitals leading to a graph with
    the parameters in question.

    EXAMPLES::

        sage: G=graphs.U42Graph216()                    # optional - gap_package_grape
        sage: G.is_strongly_regular(parameters=True)    # optional - gap_package_grape
        (216, 40, 4, 8)
    """
    from sage.libs.gap.libgap import libgap
    from sage.features.gap import GapPackage

    GapPackage("grape", spkg='gap_packages').require()

    adj_list = libgap.function_factory("""function()
                local gg, hl, o216, a216, x, h, re, G;
                LoadPackage("grape");
                gg:=SpecialUnitaryGroup(4,2);
                hl:=Z(2)*[
                [0,0,1,0],
                [1,1,0,0],
                [0,1,0,1],
                [0,1,1,0],
                [1,1,0,1]];
                o216:=Orbit(gg,Set(hl),OnSets);
                a216:=Action(gg,o216,OnSets);
                h:=Stabilizer(a216,1);
                re:=Filtered(Orbits(h,[1..216]),x->Length(x)=20);
                G:=EdgeOrbitsGraph(a216, [[1,re[1][1]], [1,re[2][1]]]);
                return List([1..216],x->Adjacency(G,x));
                end;""")

    adj = adj_list()  # for each vertex, we get the list of vertices it is adjacent to
    G = Graph(((i, int(j - 1)) for i, ni in enumerate(adj) for j in ni),
              format='list_of_edges', multiedges=False)
    G.name('U42Graph216')
    return G


def U42Graph540():
    r"""
    Return a (540,187,58,68)-strongly regular graph from [CRS2016]_.

    Build the graph, interpreting the `U_4(2)`-action considered in [CRS2016]_
    as the action of `U_4(2)=Sp_4(3)<U_4(3)` on the nonsingular, w.r.t.  to the
    Hermitean form stabilised by `U_4(3)`, points of the 3-dimensional
    projective space over `GF(9)`. There are several possible mergings of
    orbitals, some leading to non-isomorphic graphs with the same parameters. We
    found the merging here using [FK1991]_.

    EXAMPLES::

        sage: G = graphs.U42Graph540()                  # optional - gap_package_grape
        sage: G.is_strongly_regular(parameters=True)    # optional - gap_package_grape
        (540, 187, 58, 68)
    """
    from sage.libs.gap.libgap import libgap
    from sage.features.gap import GapPackage

    GapPackage("grape", spkg='gap_packages').require()

    adj_list = libgap.function_factory("""function()
                local f, o540, a540, x, oh, h, lo, G;
                LoadPackage("grape");
                f:=Sp(4,3);
                o540:=Orbit(f,Z(3)^0*[1,0,0,Z(9)],OnLines);
                a540:=Action(f,o540,OnLines);
                h:=Stabilizer(a540,1);
                oh:=Orbits(h,[1..540]);
                lo:=List([8,9,10,11,12,16,19,22,23,24],x->[1,oh[x+1][1]]);
                G:=EdgeOrbitsGraph(a540,lo);
                return List([1..540],x->Adjacency(G,x));
                end;""")

    adj = adj_list()  # for each vertex, we get the list of vertices it is adjacent to
    G = Graph(((i, int(j - 1)) for i, ni in enumerate(adj) for j in ni),
              format='list_of_edges', multiedges=False)
    G.name('U42Graph540')
    return G
