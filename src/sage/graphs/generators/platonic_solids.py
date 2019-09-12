# -*- coding: utf-8 -*-
r"""
Platonic solids

The methods defined here appear in :mod:`sage.graphs.graph_generators`.

"""
###########################################################################
#
#           Copyright (C) 2006 Robert L. Miller <rlmillster@gmail.com>
#                              and Emily A. Kirkman
#           Copyright (C) 2009 Michael C. Yurko <myurko@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
###########################################################################

# import from Sage library
from sage.graphs.graph import Graph
from math import sin, cos, pi

def TetrahedralGraph():
    """
    Returns a tetrahedral graph (with 4 nodes).

    A tetrahedron is a 4-sided triangular pyramid. The tetrahedral
    graph corresponds to the connectivity of the vertices of the
    tetrahedron. This graph is equivalent to a wheel graph with 4 nodes
    and also a complete graph on four nodes. (See examples below).

    PLOTTING: The tetrahedral graph should be viewed in 3 dimensions.
    We chose to use the default spring-layout algorithm here, so that
    multiple iterations might yield a different point of reference for
    the user. We hope to add rotatable, 3-dimensional viewing in the
    future. In such a case, a string argument will be added to select
    the flat spring-layout over a future implementation.

    EXAMPLES: Construct and show a Tetrahedral graph

    ::

        sage: g = graphs.TetrahedralGraph()
        sage: g.show() # long time

    The following example requires networkx::

        sage: import networkx as NX

    Compare this Tetrahedral, Wheel(4), Complete(4), and the
    Tetrahedral plotted with the spring-layout algorithm below in a
    Sage graphics array::

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
        sage: G.show() # long time
    """
    import networkx
    G = networkx.tetrahedral_graph()
    return Graph(G, name="Tetrahedron", pos =
                       { 0 : (0, 0),
                         1 : (0, 1),
                         2 : (cos(3.5*pi/3), sin(3.5*pi/3)),
                         3 : (cos(5.5*pi/3), sin(5.5*pi/3))}
                       )

def HexahedralGraph():
    """
    Returns a hexahedral graph (with 8 nodes).

    A regular hexahedron is a 6-sided cube. The hexahedral graph
    corresponds to the connectivity of the vertices of the hexahedron.
    This graph is equivalent to a 3-cube.

    PLOTTING: The hexahedral graph should be viewed in 3 dimensions. We
    chose to use the default spring-layout algorithm here, so that
    multiple iterations might yield a different point of reference for
    the user. We hope to add rotatable, 3-dimensional viewing in the
    future. In such a case, a string argument will be added to select
    the flat spring-layout over a future implementation.

    EXAMPLES: Construct and show a Hexahedral graph

    ::

        sage: g = graphs.HexahedralGraph()
        sage: g.show() # long time

    Create several hexahedral graphs in a Sage graphics array. They
    will be drawn differently due to the use of the spring-layout
    algorithm.

    ::

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
        sage: G.show() # long time
    """
    return Graph({0:[1,3,4], 1:[2,5], 2:[3,6], 3:[7], 4:[5,7], 5:[6], 6:[7]},
                       name="Hexahedron",
                       pos = {
                          0 : (0,0),
                          1 : (1,0),
                          3 : (0,1),
                          2 : (1,1),
                          4 : (.5,.5),
                          5 : (1.5,.5),
                          7 : (.5,1.5),
                          6 : (1.5,1.5)
                          })

def OctahedralGraph():
    """
    Returns an Octahedral graph (with 6 nodes).

    The regular octahedron is an 8-sided polyhedron with triangular
    faces. The octahedral graph corresponds to the connectivity of the
    vertices of the octahedron. It is the line graph of the tetrahedral
    graph. The octahedral is symmetric, so the spring-layout algorithm
    will be very effective for display.

    PLOTTING: The Octahedral graph should be viewed in 3 dimensions. We
    chose to use the default spring-layout algorithm here, so that
    multiple iterations might yield a different point of reference for
    the user. We hope to add rotatable, 3-dimensional viewing in the
    future. In such a case, a string argument will be added to select
    the flat spring-layout over a future implementation.

    EXAMPLES: Construct and show an Octahedral graph

    ::

        sage: g = graphs.OctahedralGraph()
        sage: g.show() # long time

    Create several octahedral graphs in a Sage graphics array They will
    be drawn differently due to the use of the spring-layout algorithm

    ::

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
        sage: G.show() # long time
    """
    import networkx
    G = networkx.octahedral_graph()

    pos = {}
    r1 = 5
    r2 = 1
    for i,v in enumerate([0,1,2]):
        i = i + 0.75
        pos[v] = (r1*cos(i*2*pi/3),r1*sin(i*2*pi/3))

    for i,v in enumerate([4,3,5]):
        i = i + .25
        pos[v] = (r2*cos(i*2*pi/3),r2*sin(i*2*pi/3))

    return Graph(G, name="Octahedron", pos=pos)

def IcosahedralGraph():
    """
    Returns an Icosahedral graph (with 12 nodes).

    The regular icosahedron is a 20-sided triangular polyhedron. The
    icosahedral graph corresponds to the connectivity of the vertices
    of the icosahedron. It is dual to the dodecahedral graph. The
    icosahedron is symmetric, so the spring-layout algorithm will be
    very effective for display.

    PLOTTING: The Icosahedral graph should be viewed in 3 dimensions.
    We chose to use the default spring-layout algorithm here, so that
    multiple iterations might yield a different point of reference for
    the user. We hope to add rotatable, 3-dimensional viewing in the
    future. In such a case, a string argument will be added to select
    the flat spring-layout over a future implementation.

    EXAMPLES: Construct and show an Octahedral graph

    ::

        sage: g = graphs.IcosahedralGraph()
        sage: g.show() # long time

    Create several icosahedral graphs in a Sage graphics array. They
    will be drawn differently due to the use of the spring-layout
    algorithm.

    ::

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
        sage: G.show() # long time
    """
    import networkx
    G = networkx.icosahedral_graph()

    pos = {}
    r1 = 5
    r2 = 2
    for i,v in enumerate([2,8,7,11,4,6]):
        i = i + .5
        pos[v] = (r1*cos(i*pi/3),r1*sin(i*pi/3))

    for i,v in enumerate([1,9,0,10,5,3]):
        i = i + .5
        pos[v] = (r2*cos(i*pi/3),r2*sin(i*pi/3))

    return Graph(G, name="Icosahedron", pos = pos)

def DodecahedralGraph():
    """
    Returns a Dodecahedral graph (with 20 nodes)

    The dodecahedral graph is cubic symmetric, so the spring-layout
    algorithm will be very effective for display. It is dual to the
    icosahedral graph.

    PLOTTING: The Dodecahedral graph should be viewed in 3 dimensions.
    We chose to use the default spring-layout algorithm here, so that
    multiple iterations might yield a different point of reference for
    the user. We hope to add rotatable, 3-dimensional viewing in the
    future. In such a case, a string argument will be added to select
    the flat spring-layout over a future implementation.

    EXAMPLES: Construct and show a Dodecahedral graph

    ::

        sage: g = graphs.DodecahedralGraph()
        sage: g.show() # long time

    Create several dodecahedral graphs in a Sage graphics array They
    will be drawn differently due to the use of the spring-layout
    algorithm

    ::

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
        sage: G.show() # long time
    """
    import networkx
    G = networkx.dodecahedral_graph()

    pos = {}
    r1 = 7
    r2 = 4.7
    r3 = 3.8
    r4 = 1.5

    for i,v in enumerate([19,0,1,2,3]):
        i = i + .25
        pos[v] = (r1*cos(i*2*pi/5),r1*sin(i*2*pi/5))

    for i,v in enumerate([18,10,8,6,4]):
        i = i + .25
        pos[v] = (r2*cos(i*2*pi/5),r2*sin(i*2*pi/5))

    for i,v in enumerate([17,11,9,7,5]):
        i = i - .25
        pos[v] = (r3*cos(i*2*pi/5),r3*sin(i*2*pi/5))

    for i,v in enumerate([12,13,14,15,16]):
        i = i + .75
        pos[v] = (r4*cos(i*2*pi/5),r4*sin(i*2*pi/5))

    return Graph(G, name="Dodecahedron", pos=pos)

