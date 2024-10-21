# sage_setup: distribution = sagemath-graphs
r"""
Graph editor widget

This module adds an interface to ``phitigra``, a graph editor widget for
Jupyter and JupyterLab. The ``phitigra`` optional package should be installed
on your Sage installation.

AUTHORS:

- Radoslav Kirov (2009): initial editor for use with the old sage notebook

- Jean-Florent Raymond (2022-04-12): replacement with the ``phitigra`` package
"""

# ****************************************************************************
#      Copyright (C) 2022 Jean-Florent Raymond
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.lazy_import import lazy_import
from sage.features.phitigra import Phitigra
lazy_import('phitigra', 'GraphEditor', feature=Phitigra())


def graph_editor(graph=None, **display_options):
    """
    Return a graph editor widget.

    The graph editor widget can be displayed with Jupyter or JupyterLab.
    It is provided by the ``phitigra`` optional package, see
    https://github.com/jfraymond/phitigra for details about the
    possible options (changing the width/height of the canvas, the
    default size and color of vertices, etc.).

    INPUT:

    - ``graph`` -- a graph to edit (default: ``None``)

    - ``display_options`` -- options for the widget

    EXAMPLES::

        sage: e = graph_editor()            # optional - phitigra
        sage: e.show()                      # not tested

    Opening an existing graph::

        sage: G = graphs.RandomGNP(10, 0.5)
        sage: e = graph_editor(G)           # optional - phitigra
        sage: e.show()                      # not tested

    Retrieving a copy of the drawn graph::

        sage: G = graphs.RandomGNP(10, 0.5)
        sage: e = graph_editor(G)           # optional - phitigra
        sage: H = e.get_graph()             # optional - phitigra
        sage: H == G and not H is G         # optional - phitigra
        True

    Using different display options::

        sage: e = graph_editor(graphs.PetersenGraph(), width=300, height=300,       # optional - phitigra
        ....:                  default_radius=12, default_vertex_color='orange',
        ....:                  default_edge_color='#666', show_vertex_labels=False)
        sage: e.show()                                                              # not tested

    .. NOTE::

        The editor does not support multigraphs.
    """
    from .graph import Graph

    if graph is None:
        graph = Graph(0)

    return GraphEditor(graph, **display_options)
