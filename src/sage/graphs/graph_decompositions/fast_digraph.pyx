r"""
Compact structure for fast operations on less than 32 vertices

This module implements a digraph structure meant to be used in Cython in
**highly enumerative** algorithms. It can store graphs on less than
``sizeof(int)`` vertices and perform several basic operations **quickly**
(add/remove arcs, count the out-neighborhood of a set of vertices or return its
cardinality).

**Sets and integers :**

In the following code, sets are represented as integers, where the `i`-th bit is
set if element `i` belongs to the set.
"""
from cysignals.memory cimport check_allocarray, check_calloc, sig_free


cdef class FastDigraph:

    def __cinit__(self, D, vertex_list=None):
        r"""
        Constructor for ``FastDigraph``.

        If the input parameter ``D`` is a Graph, it is handled as a symmetric
        DiGraph.

        INPUT:

        - ``D`` -- a (Di)Graph

        - ``vertex_list`` -- list (default: ``None``); specifies a mapping
          between `[0..n-1]` and the set of vertices of the input (Di)Graph,
          ``list(D)`` by default

        EXAMPLES::

            sage: cython_code = [
            ....: 'from sage.graphs.graph import Graph',
            ....: 'from sage.graphs.graph_decompositions.fast_digraph cimport FastDigraph',
            ....: 'G = Graph([(0, 1), (1, 2)])',
            ....: 'cdef FastDigraph F = FastDigraph(G)',
            ....: 'cdef int i',
            ....: 'print([F.degree[i] for i in range(F.n)])']
            sage: cython(os.linesep.join(cython_code))                                  # needs sage.misc.cython
            [1, 2, 1]
        """
        if D.order() > 8*sizeof(int):
            raise OverflowError("Too many vertices. This structure can only "
                                "encode digraphs on at most "
                                "%i vertices" % (8 * sizeof(int)))

        self.n = D.order()
        self.graph = <int *>check_calloc(self.n, sizeof(int))

        cdef int i
        cdef int tmp

        # When the vertices are not consecutive integers
        if vertex_list is None:
            self.int_to_vertices = list(D)
        elif len(vertex_list) == self.n and not set(vertex_list).symmetric_difference(D):
            self.int_to_vertices = list(vertex_list)
        else:
            raise ValueError("the input vertex_list is incorrect")
        cdef dict vertices_to_int = {v: i for i, v in enumerate(self.int_to_vertices)}

        if D.is_directed():
            for u in D:
                tmp = 0
                for v in D.neighbors_out(u):
                    tmp |= 1 << vertices_to_int[v]
                self.graph[vertices_to_int[u]] = tmp
        else:
            for u in D:
                tmp = 0
                for v in D.neighbors(u):
                    tmp |= 1 << vertices_to_int[v]
                self.graph[vertices_to_int[u]] = tmp

        self.degree = <int *>check_allocarray(self.n, sizeof(int))
        for i in range(self.n):
            self.degree[i] = popcount32(self.graph[i])

    def __dealloc__(self):
        r"""
        Destructor.
        """
        sig_free(self.graph)
        sig_free(self.degree)

    def print_adjacency_matrix(self):
        r"""
        Displays the adjacency matrix of ``self``.

        EXAMPLES::

            sage: cython_code = [
            ....: 'from sage.graphs.graph import Graph',
            ....: 'from sage.graphs.graph_decompositions.fast_digraph cimport FastDigraph',
            ....: 'FastDigraph(Graph([(0, 1), (1, 2)])).print_adjacency_matrix()']
            sage: cython(os.linesep.join(cython_code))                                  # needs sage.misc.cython
            010
            101
            010
        """
        cdef int i, j
        for i in range(self.n):
            for j in range(self.n):
                print(((self.graph[i] >> j) & 1), end="")
            print("")

cdef inline int compute_out_neighborhood_cardinality(FastDigraph g, int S) noexcept:
    r"""
    Return the cardinality of `N^+(S)\S`.

    INPUT:

    - ``g`` -- a FastDigraph

    - ``S`` -- integer describing the set

    EXAMPLES::

        sage: cython_code = [
        ....: 'from sage.graphs.graph import Graph',
        ....: 'from sage.graphs.graph_decompositions.fast_digraph cimport FastDigraph',
        ....: 'from sage.graphs.graph_decompositions.fast_digraph cimport compute_out_neighborhood_cardinality',
        ....: 'cdef FastDigraph F = FastDigraph(Graph([(0, 1), (1, 2)]))',
        ....: 'cdef int i',
        ....: 'print([compute_out_neighborhood_cardinality(F, 1<<i) for i in range(F.n)])']
        sage: cython(os.linesep.join(cython_code))                                      # needs sage.misc.cython
        [1, 2, 1]
    """
    cdef int i
    cdef int tmp = 0
    for i in range(g.n):
        tmp |= g.graph[i] & (-((S >> i) & 1))

    tmp &= (~S)
    return popcount32(tmp)

cdef inline int popcount32(int i) noexcept:
    r"""
    Return the number of '1' bits in a 32-bits integer.

    If ``sizeof(int) > 4``, this function only returns the number of '1'
    bits in ``(i & ((1<<32) - 1))``.

    EXAMPLES::

        sage: cython_code = [
        ....: 'from sage.graphs.graph_decompositions.fast_digraph cimport popcount32',
        ....: 'cdef int i',
        ....: 'print([popcount32(i) for i in range(16)])']
        sage: cython(os.linesep.join(cython_code))                                      # needs sage.misc.cython
        [0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4]
    """
    i = i - ((i >> 1) & 0x55555555)
    i = (i & 0x33333333) + ((i >> 2) & 0x33333333)
    return ((i + (i >> 4) & 0x0F0F0F0F) * 0x01010101) >> 24


# If you happened to be doubting the consistency of the popcount32 function
# above, you can give the following doctest a try. It is not tested
# automatically by Sage as it takes a *LONG* time to run (around 5 minutes), but
# it would report any problem if it finds one.

def test_popcount():
    """
    Correction test for popcount32.

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.fast_digraph import test_popcount
        sage: test_popcount() # not tested
    """
    cdef int i = 1
    # While the last 32 bits of i are not equal to 0
    while (i & ((1 << 32) - 1)):
        if popcount32(i) != slow_popcount32(i):
            print("Error for i = ", str(i))
            print("Result with popcount32 : " + str(popcount32(i)))
            print("Result with slow_popcount32 : " + str(slow_popcount32(i)))
        i += 1


cdef inline int slow_popcount32(int i) noexcept:
    """
    Return the number of '1' bits in a 32-bits integer.

    If ``sizeof(int) > 4``, this function only returns the number of '1'
    bits in ``(i & ((1<<32) - 1))``.

    EXAMPLES::

        sage: cython_code = [
        ....: 'from sage.graphs.graph_decompositions.fast_digraph cimport popcount32',
        ....: 'from sage.graphs.graph_decompositions.fast_digraph cimport slow_popcount32',
        ....: 'cdef int i',
        ....: 'print(all(popcount32(i) == slow_popcount32(i) for i in range(16)))']
        sage: cython(os.linesep.join(cython_code))                                      # needs sage.misc.cython
        True
    """
    # Slow popcount for 32bits integers
    cdef int j = 0
    cdef int k

    for k in range(32):
        j += (i >> k) & 1

    return j
