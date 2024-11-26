from libcpp cimport bool
from libcpp.list cimport list as cpplist
from libcpp.vector cimport vector

cdef extern from "modular_decomposition.hpp":
    cdef cppclass SDData:
        void set_from_data(size_t lex_label_offset, const int* sigma,
                           const size_t *xslice_len,
                           const vector[int] *lex_label)

    cdef cppclass md_tree_node:
        bool is_leaf() const
        bool is_prime() const
        bool is_parallel() const
        bool is_series() const
        cpplist[md_tree_node *] children
        # For a leaf, the corresponding vertex, for a internal node, any vertex
        # corresponding to a any leaf below the node
        int vertex

    void dealloc_md_tree_nodes_recursively(md_tree_node *)

    cdef md_tree_node * corneil_habib_paul_tedder_inner(const SDData &SD)
