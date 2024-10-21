# sage_setup: distribution = sagemath-graphs
"""
Print graph

Copying old scheme print-graph methods over to Sage

AUTHORS:

- Craig Citro (2007-06-12) --  initial version
"""

# ****************************************************************************
#                        Copyright (C) 2007 Craig Citro
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
# ****************************************************************************

import time
from math import floor
from copy import copy


def print_header_ps(s):
    """
    Give the header for a postscript file.

    EXAMPLES::

        sage: from sage.graphs.print_graphs import print_header_ps
        sage: print(print_header_ps(''))
        %% --- Auto-generated PostScript ---
        %% Generated on:
        %%...
    """
    s += "%% --- Auto-generated PostScript ---\n\n\n"
    s += "%% Generated on: \n"
    s += "%%" + time.asctime() + "\n"
    return s


def print_header_eps(s, xmin, ymin, xmax, ymax):
    """
    Give the header for an encapsulated postscript file.

    EXAMPLES::

        sage: from sage.graphs.print_graphs import print_header_eps
        sage: print(print_header_eps('',0,0,1,1))
        %!PS-Adobe-3.0 EPSF-3.0
        %%BoundingBox: 0 0 1 1
    """
    s += "%!PS-Adobe-3.0 EPSF-3.0\n"
    s += "%" + "%" + "BoundingBox: %s %s %s %s \n" % (xmin, ymin, xmax, ymax)

    return s


def print_functions(s):
    """
    Define edge and point drawing functions.

    EXAMPLES::

        sage: from sage.graphs.print_graphs import print_functions
        sage: print(print_functions(''))
        /point %% input: x y
        { moveto
          gsave
          currentpoint translate
          0 0 2 0 360 arc
          fill
          grestore
          } def
        /edge %% input: x1 y1 x2 y2
        { moveto
          lineto
          stroke
          } def
    """
    s += "/point %% input: x y\n"
    s += "{ moveto\n"
    s += "  gsave\n"
    s += "  currentpoint translate\n"
    s += "  0 0 2 0 360 arc\n"
    s += "  fill\n"
    s += "  grestore\n"
    s += "  } def\n\n\n"
    s += "/edge %% input: x1 y1 x2 y2\n"
    s += "{ moveto\n"
    s += "  lineto\n"
    s += "  stroke\n"
    s += "  } def\n\n"

    return s


def print_graph_ps(vert_ls, edge_iter, pos_dict):
    """
    Give postscript text for drawing a graph.

    EXAMPLES::

        sage: from sage.graphs.print_graphs import print_graph_ps
        sage: P = graphs.PetersenGraph()
        sage: print(print_graph_ps(P.vertices(sort=True), P.edges(sort=True), sage.graphs.generic_graph_pyx.spring_layout_fast(P)))
        %% --- Auto-generated PostScript ---
        %% Generated on:
        %%...
        /point %% input: x y
        { moveto
          gsave
          currentpoint translate
          0 0 2 0 360 arc
          fill
          grestore
          } def
        /edge %% input: x1 y1 x2 y2
        { moveto
          lineto
          stroke
          } def
        ... point
        ...
        ... point
        ... edge
        ...
        ... edge
    """
    pos_dict = copy(pos_dict)  # assumption: all pos's are -1 <= ... <= 1

    s = ""

    s = print_header_ps(s)
    s = print_functions(s)

    for v in vert_ls:
        x, y = pos_dict[v]
        pos_dict[v] = int(floor(50*x))+50, int(floor(50*y))+50
        x, y = pos_dict[v]
        s += "%s %s point\n" % (x, y)

    for (u, v, l) in edge_iter:
        ux, uy = pos_dict[u]
        vx, vy = pos_dict[v]
        s += "%s %s %s %s edge\n" % (ux, uy, vx, vy)

    return s


def print_graph_eps(vert_ls, edge_iter, pos_dict):
    """
    Give postscript text for drawing a graph.

    EXAMPLES::

        sage: from sage.graphs.print_graphs import print_graph_eps
        sage: P = graphs.PetersenGraph()
        sage: print(print_graph_eps(P.vertices(sort=True), P.edges(sort=True), sage.graphs.generic_graph_pyx.spring_layout_fast(P)))
        %!PS-Adobe-3.0 EPSF-3.0
        %%BoundingBox: 0 0 100 100
        /point %% input: x y
        { moveto
          gsave
          currentpoint translate
          0 0 2 0 360 arc
          fill
          grestore
          } def
        /edge %% input: x1 y1 x2 y2
        { moveto
          lineto
          stroke
          } def
        ... point
        ...
        ... point
        ... edge
        ...
        ... edge
    """
    pos_dict = copy(pos_dict)  # assumption: all pos's are -1 <= ... <= 1

    t = ""
    s = ""

    for v in vert_ls:
        x, y = pos_dict[v]
        pos_dict[v] = int(floor(50*x)) + 50, int(floor(50*y)) + 50
        x, y = pos_dict[v]
        s += "%s %s point\n" % (x, y)

    for (u, v, l) in edge_iter:
        ux, uy = pos_dict[u]
        vx, vy = pos_dict[v]
        s += "%s %s %s %s edge\n" % (ux, uy, vx, vy)

    t = print_header_eps(t, 0, 0, 100, 100)
    t = print_functions(t)

    return t + s
