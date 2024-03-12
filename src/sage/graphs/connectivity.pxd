from memory_allocator cimport MemoryAllocator
from sage.graphs.generic_graph_pyx cimport GenericGraph_pyx

ctypedef struct _LinkedListNode:
    _LinkedListNode * prev
    _LinkedListNode * next
    Py_ssize_t data

ctypedef struct _LinkedList:
    _LinkedListNode * head
    _LinkedListNode * tail
    Py_ssize_t length

cdef class _Component:
    cdef MemoryAllocator mem
    cdef _LinkedList * edge_list
    cdef int component_type

    cdef add_edge(self, Py_ssize_t e_index)
    cdef finish_tric_or_poly(self, Py_ssize_t e_index)
    cdef list get_edge_list(self)

cdef class TriconnectivitySPQR:
    cdef MemoryAllocator mem
    cdef Py_ssize_t n
    cdef Py_ssize_t m
    cdef Py_ssize_t max_number_of_edges
    cdef str graph_name

    # Relabel vertices in range [0..n-1]
    cdef list int_to_vertex
    cdef dict vertex_to_int

    # We associate a unique identifier (int) to each edge and store extremities
    # and original labels in distinct arrays.
    cdef int* edge_extremity_first
    cdef int* edge_extremity_second
    cdef list int_to_original_edge_label
    cdef int virtual_edge_num  # number of created virtual edges

    cdef int* edge_status
    cdef bint* reverse_edges

    cdef int* dfs_number

    cdef _LinkedList ** highpt
    cdef _LinkedListNode ** in_high

    # Translates DFS number of a vertex to its new number
    cdef int* old_to_new
    cdef int* newnum   # new number of vertex i
    cdef int* node_at  # node at dfs number of i
    cdef int* lowpt1   # lowpt1 number of vertex i
    cdef int* lowpt2   # lowpt2 number of vertex i

    cdef _LinkedList ** adj
    cdef _LinkedListNode ** in_adj
    cdef int* nd        # Number of descendants of vertex i
    cdef int* parent    # Parent vertex of vertex i in the palm tree
    cdef int* degree     # Degree of vertex i
    cdef int* tree_arc   # Tree arc entering the vertex i
    cdef int* vertex_at  # vertex with DFS number of i

    cdef int dfs_counter
    cdef list components_list  # list of components of `graph_copy`
    cdef list graph_copy_adjacency  # Stores adjacency list

    cdef bint* starts_path   # Does edge e_index start a path
    cdef int start_vertex    # First vertex of exploration

    # Stacks used in `path_search` function
    cdef list e_stack
    cdef int* t_stack_h
    cdef int* t_stack_a
    cdef int* t_stack_b
    cdef int t_stack_top

    cdef list comp_final_edge_list  # entry i is list of edges in component i
    cdef list comp_type             # entry i is type of component i
    cdef dict final_edge_to_edge_index  # associate final edge e to its index in int_to_edge
    cdef GenericGraph_pyx spqr_tree  # The final SPQR tree is stored

    # Arrays used in different methods. Allocated only once
    cdef int* tmp_array_n_int_1
    cdef int* tmp_array_n_int_2
    cdef int* tmp_array_n_int_3
    cdef bint* tmp_array_n_bint_1

    # === Methods ===

    cdef inline __tstack_push(self, int h, int a, int b):
        """
        Push ``(h, a, b)`` triple on ``Tstack``.
        """
        self.t_stack_top += 1
        self.t_stack_h[self.t_stack_top] = h
        self.t_stack_a[self.t_stack_top] = a
        self.t_stack_b[self.t_stack_top] = b

    cdef inline __tstack_push_eos(self):
        """
        Push end-of-stack marker on ``Tstack``.
        """
        self.t_stack_top += 1
        self.t_stack_a[self.t_stack_top] = -1

    cdef inline bint __tstack_not_eos(self) noexcept:
        """
        Return ``True`` iff end-of-stack marker is not on top of ``Tstack``.
        """
        return self.t_stack_a[self.t_stack_top] != -1

    cdef inline int __estack_pop(self) noexcept:
        """
        Pop from estack and return the popped element
        """
        return <int> self.e_stack.pop()

    cdef inline __new_component(self, list edges, int type_c):
        """
        Create a new component and add ``edges`` to it.

        ``type_c = 0`` for bond, ``1`` for polygon, ``2`` for triconnected.
        """
        self.components_list.append(_Component(edges, type_c))

    cdef inline bint __is_virtual_edge(self, int e_index) noexcept:
        """
        Return ``True`` if edge number ``e_index`` is a virtual edge.

        By construction, the first ``m`` edge indices are the original edges of
        the graph. The edges created during the execution of the algorithm,
        i.e., virtual edges, have indices ``>= m``.
        """
        return e_index >= self.m

    cdef inline int __edge_other_extremity(self, int e_index, int u) noexcept:
        """
        Return the other extremity of the edge
        """
        if self.edge_extremity_first[e_index] == u:
            return self.edge_extremity_second[e_index]
        return self.edge_extremity_first[e_index]

    cdef int __new_virtual_edge(self, int u, int v) noexcept
    cdef _LinkedListNode * __new_LinkedListNode(self, Py_ssize_t e_index) noexcept
    cdef Py_ssize_t __high(self, Py_ssize_t v) noexcept
    cdef __del_high(self, int e_index)
    cdef __split_multiple_edges(self)
    cdef int __dfs1(self, int start, bint check=*) noexcept
    cdef __build_acceptable_adj_struct(self)
    cdef __path_finder(self, int start)
    cdef __dfs2(self)
    cdef int __path_search(self, int start) except -1
    cdef __assemble_triconnected_components(self)
    cdef __build_spqr_tree(self)
