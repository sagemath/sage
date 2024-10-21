# sage_setup: distribution = sagemath-graphs
r"""
1-skeletons of Platonic solids

The methods defined here appear in :mod:`sage.graphs.graph_generators`.
"""
# ****************************************************************************
#
#           Copyright (C) 2006 Robert L. Miller <rlmillster@gmail.com>
#                              and Emily A. Kirkman
#           Copyright (C) 2009 Michael C. Yurko <myurko@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
# ****************************************************************************

# import from Sage library
from sage.graphs.graph import Graph
from math import sin, cos, pi


def TetrahedralGraph():
    """
    Return a tetrahedral graph (with 4 nodes).

    A tetrahedron is a 4-sided triangular pyramid. The tetrahedral graph
    corresponds to the connectivity of the vertices of the tetrahedron. This
    graph is equivalent to a wheel graph with 4 nodes and also a complete graph
    on four nodes. (See examples below).

    PLOTTING: The Tetrahedral graph should be viewed in 3 dimensions. We choose
    to use a planar embedding of the graph. We hope to add rotatable,
    3-dimensional viewing in the future. In such a case, a argument will be
    added to select the desired layout.

    EXAMPLES:

    Construct and show a Tetrahedral graph::

        sage: g = graphs.TetrahedralGraph()
        sage: g.show()                          # long time                             # needs sage.plot

    The following example requires networkx::

        sage: import networkx as NX                                                     # needs networkx

    Compare this Tetrahedral, Wheel(4), Complete(4), and the Tetrahedral plotted
    with the spring-layout algorithm below in a Sage graphics array::

        sage: # needs networkx sage.plot
        sage: tetra_pos = graphs.TetrahedralGraph()
        sage: tetra_spring = Graph(NX.tetrahedral_graph())
        sage: wheel = graphs.WheelGraph(4)
        sage: complete = graphs.CompleteGraph(4)
        sage: g = [tetra_pos, tetra_spring, wheel, complete]
        sage: j = []
        sage: for i in range(2):
        ....:     n = []
        ....:     for m in range(2):
        ....:         n.append(g[i + m].plot(vertex_size=50, vertex_labels=False))
        ....:     j.append(n)
        sage: G = graphics_array(j)
        sage: G.show()                          # long time
    """
    edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    pos = {0: (0, 0),
           1: (0, 1),
           2: (cos(3.5*pi/3), sin(3.5*pi/3)),
           3: (cos(5.5*pi/3), sin(5.5*pi/3))}
    return Graph(edges, name='Tetrahedron', pos=pos)


def HexahedralGraph():
    """
    Return a hexahedral graph (with 8 nodes).

    A regular hexahedron is a 6-sided cube. The hexahedral graph corresponds to
    the connectivity of the vertices of the hexahedron. This graph is
    equivalent to a 3-cube.

    PLOTTING: The Hexahedral graph should be viewed in 3 dimensions. We choose
    to use a planar embedding of the graph. We hope to add rotatable,
    3-dimensional viewing in the future. In such a case, a argument will be
    added to select the desired layout.

    EXAMPLES:

    Construct and show a Hexahedral graph::

        sage: g = graphs.HexahedralGraph()
        sage: g.show()                          # long time                             # needs sage.plot

    Create several hexahedral graphs in a Sage graphics array. They will be
    drawn differently due to the use of the spring-layout algorithm::

        sage: # needs sage.plot
        sage: g = []
        sage: j = []
        sage: for i in range(9):
        ....:     k = graphs.HexahedralGraph()
        ....:     g.append(k)
        sage: for i in range(3):
        ....:     n = []
        ....:     for m in range(3):
        ....:         n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ....:     j.append(n)
        sage: G = graphics_array(j)
        sage: G.show()                          # long time
    """
    adj = {0: [1, 3, 4], 1: [2, 5], 2: [3, 6], 3: [7], 4: [5, 7], 5: [6], 6: [7]}
    pos = {
        0: (0, 0),
        1: (1, 0),
        3: (0, 1),
        2: (1, 1),
        4: (.5, .5),
        5: (1.5, .5),
        7: (.5, 1.5),
        6: (1.5, 1.5)
        }
    return Graph(adj, name='Hexahedron', pos=pos)


def OctahedralGraph():
    """
    Return an Octahedral graph (with 6 nodes).

    The regular octahedron is an 8-sided polyhedron with triangular faces. The
    octahedral graph corresponds to the connectivity of the vertices of the
    octahedron. It is the line graph of the tetrahedral graph. The octahedral is
    symmetric, so the spring-layout algorithm will be very effective for
    display.

    PLOTTING: The Octahedral graph should be viewed in 3 dimensions. We choose
    to use a planar embedding of the graph. We hope to add rotatable,
    3-dimensional viewing in the future. In such a case, a argument will be
    added to select the desired layout.

    EXAMPLES:

    Construct and show an Octahedral graph::

        sage: g = graphs.OctahedralGraph()
        sage: g.show()                          # long time                             # needs sage.plot

    Create several octahedral graphs in a Sage graphics array They will be drawn
    differently due to the use of the spring-layout algorithm::

        sage: # needs sage.plot
        sage: g = []
        sage: j = []
        sage: for i in range(9):
        ....:     k = graphs.OctahedralGraph()
        ....:     g.append(k)
        sage: for i in range(3):
        ....:     n = []
        ....:     for m in range(3):
        ....:         n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ....:     j.append(n)
        sage: G = graphics_array(j)
        sage: G.show()                          # long time
    """
    adj = {0: [1, 2, 3, 4], 1: [2, 3, 5], 2: [4, 5], 3: [4, 5], 4: [5]}
    G = Graph(adj, format='dict_of_lists', name='Octahedron')
    G._circle_embedding([0, 1, 2], radius=5, angle=pi/2)
    G._circle_embedding([4, 3, 5], radius=1, angle=pi/6)
    return G


def IcosahedralGraph():
    """
    Return an Icosahedral graph (with 12 nodes).

    The regular icosahedron is a 20-sided triangular polyhedron. The icosahedral
    graph corresponds to the connectivity of the vertices of the icosahedron. It
    is dual to the dodecahedral graph. The icosahedron is symmetric, so the
    spring-layout algorithm will be very effective for display.

    PLOTTING: The Icosahedral graph should be viewed in 3 dimensions. We choose
    to use a planar embedding of the graph. We hope to add rotatable,
    3-dimensional viewing in the future. In such a case, a argument will be
    added to select the desired layout.

    EXAMPLES:

    Construct and show an Octahedral graph::

        sage: g = graphs.IcosahedralGraph()
        sage: g.show()                          # long time                             # needs sage.plot

    Create several icosahedral graphs in a Sage graphics array. They will be
    drawn differently due to the use of the spring-layout algorithm::

        sage: # needs sage.plot
        sage: g = []
        sage: j = []
        sage: for i in range(9):
        ....:     k = graphs.IcosahedralGraph()
        ....:     g.append(k)
        sage: for i in range(3):
        ....:     n = []
        ....:     for m in range(3):
        ....:         n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ....:     j.append(n)
        sage: G = graphics_array(j)
        sage: G.show()                          # long time
    """
    adj = {0: [1, 5, 7, 8, 11], 1: [2, 5, 6, 8], 2: [3, 6, 8, 9],
           3: [4, 6, 9, 10], 4: [5, 6, 10, 11], 5: [6, 11],
           7: [8, 9, 10, 11], 8: [9], 9: [10], 10: [11]}
    G = Graph(adj, format='dict_of_lists', name='Icosahedron')
    G._circle_embedding([2, 8, 7, 11, 4, 6], radius=5, angle=pi/6)
    G._circle_embedding([1, 9, 0, 10, 5, 3], radius=2, angle=pi/6)
    return G


def DodecahedralGraph():
    """
    Return a Dodecahedral graph (with 20 nodes).

    The dodecahedral graph is cubic symmetric, so the spring-layout algorithm
    will be very effective for display. It is dual to the icosahedral graph.

    PLOTTING: The Dodecahedral graph should be viewed in 3 dimensions. We
    choose to use a planar embedding of the graph. We hope to add rotatable,
    3-dimensional viewing in the future. In such a case, a argument will be
    added to select the desired layout.

    EXAMPLES:

    Construct and show a Dodecahedral graph::

        sage: g = graphs.DodecahedralGraph()
        sage: g.show()                          # long time                             # needs sage.plot

    Create several dodecahedral graphs in a Sage graphics array They will be
    drawn differently due to the use of the spring-layout algorithm::

        sage: # needs sage.plot
        sage: g = []
        sage: j = []
        sage: for i in range(9):
        ....:     k = graphs.DodecahedralGraph()
        ....:     g.append(k)
        sage: for i in range(3):
        ....:     n = []
        ....:     for m in range(3):
        ....:         n.append(g[3*i + m].plot(vertex_size=50, vertex_labels=False))
        ....:     j.append(n)
        sage: G = graphics_array(j)
        sage: G.show()                          # long time
    """
    adj = {0: [1, 10, 19], 1: [2, 8], 2: [3, 6], 3: [4, 19], 4: [5, 17],
           5: [6, 15], 6: [7], 7: [8, 14], 8: [9], 9: [10, 13], 10: [11],
           11: [12, 18], 12: [13, 16], 13: [14], 14: [15], 15: [16], 16: [17],
           17: [18], 18: [19]}
    G = Graph(adj, format='dict_of_lists', name='Dodecahedron')
    G._circle_embedding([19, 0, 1, 2, 3], radius=7, angle=pi/10)
    G._circle_embedding([18, 10, 8, 6, 4], radius=4.7, angle=pi/10)
    G._circle_embedding([11, 9, 7, 5, 17], radius=3.8, angle=3*pi/10)
    G._circle_embedding([12, 13, 14, 15, 16], radius=1.5, angle=3*pi/10)
    return G
