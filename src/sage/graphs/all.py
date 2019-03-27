from __future__ import absolute_import

from sage.misc.lazy_import import lazy_import

lazy_import("sage.graphs.graph_generators", "graphs")
lazy_import("sage.graphs.digraph_generators", "digraphs")
lazy_import("sage.graphs.hypergraph_generators", "hypergraphs")
from .graph_database import GraphDatabase, GenericGraphQuery, GraphQuery
from .graph import Graph
from .digraph import DiGraph
from .bipartite_graph import BipartiteGraph
import sage.graphs.weakly_chordal
import sage.graphs.lovasz_theta
import sage.graphs.partial_cube
from . import graph_list as graphs_list
lazy_import("sage.graphs", "graph_coloring")
lazy_import("sage.graphs.cliquer", ['all_max_clique', 'max_clique',
                                    'clique_number'],
            deprecation=26200)
from .graph_database import graph_db_info
lazy_import("sage.graphs.graph_editor", "graph_editor")

from sage.graphs.isgci import graph_classes

"""
TESTS:

Test that sagenb.misc.support is not imported (see :trac:`22941`)::

    sage: import sage.graphs.graph_editor
    sage: 'sagenb.misc.support' in sys.modules
    False

Test that methods all_max_clique, max_clique and clique_number from
sage.graphs.cliquer are deprecated from the global namespace (:trac:`26200`)::

    sage: all_max_clique(Graph())
    doctest:...: DeprecationWarning: Importing all_max_clique from here is deprecated. If you need to use it, please import it directly from sage.graphs.cliquer
    See https://trac.sagemath.org/26200 for details.
    [[]]
    sage: max_clique(Graph())
    doctest:...: DeprecationWarning: Importing max_clique from here is deprecated. If you need to use it, please import it directly from sage.graphs.cliquer
    See https://trac.sagemath.org/26200 for details.
    []
    sage: clique_number(Graph())
    doctest:...: DeprecationWarning: Importing clique_number from here is deprecated. If you need to use it, please import it directly from sage.graphs.cliquer
    See https://trac.sagemath.org/26200 for details.
    0
"""
