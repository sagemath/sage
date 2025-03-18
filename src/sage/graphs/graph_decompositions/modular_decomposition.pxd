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
        # If is_leaf() is true, the attribute 'vertex' contains the id of the
        # corresponding vertex. If is_leaf() is false, the attribute 'vertex'
        # contains the id of a vertex corresponding to any leaf below
        # the node (i.e., 'vertex' contains the id of a vertex belonging to the
        # module corresponding to the node).
        int vertex

    void dealloc_md_tree_nodes_recursively(md_tree_node *)

    cdef md_tree_node * corneil_habib_paul_tedder_inner(const SDData &SD)
