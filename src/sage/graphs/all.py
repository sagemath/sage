
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
from .graph_database import graph_db_info
lazy_import("sage.graphs.graph_editor", "graph_editor")

from sage.graphs.isgci import graph_classes

"""
TESTS:

Test that methods all_max_clique, max_clique and clique_number from
sage.graphs.cliquer are deprecated from the global namespace (:trac:`26200`)::

    sage: all_max_clique(Graph())
    Traceback (most recent call last):
    ...
    NameError: name 'all_max_clique' is not defined
    sage: max_clique(Graph())
    Traceback (most recent call last):
    ...
    NameError: name 'max_clique' is not defined
    sage: clique_number(Graph())
    Traceback (most recent call last):
    ...
    NameError: name 'clique_number' is not defined
"""
