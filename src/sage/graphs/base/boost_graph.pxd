# distutils: language = c++
# distutils: extra_compile_args = -std=c++11

#*****************************************************************************
#       Copyright (C) 2015 Michele Borassi michele.borassi@imtlucca.it
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from libcpp.vector cimport vector
from libcpp.pair cimport pair

cdef extern from "boost/graph/adjacency_list.hpp" namespace "boost":
    cdef cppclass vecS:
        pass
    cdef cppclass listS:
        pass
    cdef cppclass slistS:
        pass
    cdef cppclass setS:
        pass
    cdef cppclass multisetS:
        pass
    cdef cppclass hash_setS:
        pass
    cdef cppclass undirectedS:
        pass
    cdef cppclass directedS:
        pass
    cdef cppclass bidirectionalS:
        pass
    cdef cppclass no_property:
        pass
    cdef cppclass edge_weight_t:
        pass
    cdef cppclass property[T, U]:
        pass


cdef extern from "boost_interface.cpp":
    ctypedef unsigned long v_index
    ctypedef unsigned long e_index

    cdef cppclass result_ec:
        v_index ec
        vector[v_index] edges

    cdef cppclass result_cc:
        float average_clustering_coefficient
        vector[float] clust_of_v

    cdef cppclass result_distances:
        vector[double] distances
        vector[v_index] predecessors

    cdef cppclass BoostGraph[OutEdgeListS, VertexListS, DirectedS, EdgeListS, EdgeProperty]:
        BoostGraph()
        void add_vertex()
        v_index num_verts()
        void add_edge(v_index u, v_index v)
        void add_edge(v_index u, v_index v, double w)
        vector[pair[v_index, pair[v_index, double]]] edge_list()
        e_index num_edges()
        result_ec edge_connectivity()
        double clustering_coeff(v_index v)
        result_cc clustering_coeff_all()
        vector[v_index] dominator_tree(v_index v)
        vector[v_index] bandwidth_ordering(bool)
        vector[v_index] kruskal_min_spanning_tree()
        vector[v_index] prim_min_spanning_tree()
        vector[vector[v_index]] blocks_and_cut_vertices()
        result_distances dijkstra_shortest_paths(v_index s)
        result_distances bellman_ford_shortest_paths(v_index s)
        vector[vector[double]] johnson_shortest_paths()
        vector[vector[double]] floyd_warshall_shortest_paths()

ctypedef property[edge_weight_t, double] EdgeWeight

ctypedef BoostGraph[vecS, vecS, undirectedS,    vecS, no_property] BoostVecGraph
ctypedef BoostGraph[vecS, vecS, bidirectionalS, vecS, no_property] BoostVecDiGraph

ctypedef BoostGraph[vecS, vecS, undirectedS,    vecS, EdgeWeight] BoostVecWeightedGraph
ctypedef BoostGraph[vecS, vecS, directedS,      vecS, EdgeWeight] BoostVecWeightedDiGraphU
ctypedef BoostGraph[vecS, vecS, bidirectionalS, vecS, EdgeWeight] BoostVecWeightedDiGraph

ctypedef BoostGraph[setS, vecS, undirectedS,    vecS, no_property] BoostSetGraph

ctypedef fused BoostVecGenGraph:
    BoostVecGraph
    BoostVecDiGraph

ctypedef fused BoostWeightedGraph:
    BoostVecWeightedGraph
    BoostVecWeightedDiGraph
    BoostVecWeightedDiGraphU

ctypedef fused BoostGenGraph:
    BoostVecGraph
    BoostVecDiGraph
    BoostGraph[vecS, vecS, directedS,      vecS, no_property]
    BoostSetGraph
    BoostGraph[setS, vecS, directedS,      vecS, no_property]
    BoostGraph[setS, vecS, bidirectionalS, vecS, no_property]
