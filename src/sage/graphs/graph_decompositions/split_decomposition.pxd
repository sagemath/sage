from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "split_decomposition.hpp":

    cdef cppclass GltExport:
        GltExport() except +
        int root_vertex
        vector[int] order
        vector[int] node_kind
        vector[int] node_parent
        vector[int] node_star_centre
        vector[int] node_marker_start
        vector[int] node_marker_neighbour
        vector[int] node_edge_start
        vector[int] node_edge_pairs

    bool glt_build(const vector[vector[int]]& graph,
                   const vector[int]& given_order,
                   int start_vertex,
                   bool run_verify,
                   GltExport& out,
                   string& error) except +
