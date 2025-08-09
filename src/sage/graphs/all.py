
from sage.graphs.isgci import graph_classes
from sage.graphs import graph_list as graphs_list
import sage.graphs.partial_cube
import sage.graphs.lovasz_theta
import sage.graphs.weakly_chordal
from sage.graphs.maps.all import *
from sage.graphs.matching_covered_graph import MatchingCoveredGraph
from sage.graphs.bipartite_graph import BipartiteGraph
from sage.graphs.digraph import DiGraph
from sage.graphs.graph import Graph
from sage.misc.lazy_import import lazy_import

lazy_import("sage.graphs.graph_generators", "graphs")
lazy_import("sage.graphs.digraph_generators", "digraphs")
lazy_import("sage.graphs.hypergraph_generators", "hypergraphs")
lazy_import("sage.graphs.graph_database", [
            "GraphDatabase", "GenericGraphQuery", "GraphQuery"])
lazy_import("sage.graphs", "graph_coloring")
lazy_import("sage.graphs.graph_database", "graph_db_info")
lazy_import("sage.graphs.graph_editor", "graph_editor")

"""
TESTS:

Test that methods all_max_clique, max_clique and clique_number from
sage.graphs.cliquer are deprecated from the global namespace (:issue:`26200`)::

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
del lazy_import
