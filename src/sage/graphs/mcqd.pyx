# distutils: language = c++
# sage_setup: distribution = sagemath-mcqd

from cysignals.signals cimport sig_on, sig_off
from memory_allocator cimport MemoryAllocator


def mcqd(G):
    """
    Compute the max clique using MCQD.

    INPUT:

    - ``G`` -- a graph

    TESTS::

        sage: from sage.graphs.mcqd import mcqd         # optional - mcqd
        sage: for i in range(10):                       # optional - mcqd
        ....:     g = graphs.RandomGNP(15,.5)
        ....:     if g.clique_number() != len(mcqd(g)):
        ....:         print("This is dead wrong !")
    """
    cdef int n = G.order()

    # - c0 is the adjacency matrix
    # - c points toward each row of the matrix
    # - qmax stores the max clique
    cdef MemoryAllocator mem = MemoryAllocator()
    cdef bool ** c = <bool **> mem.allocarray(n, sizeof(bool *))
    cdef bool * c0 = <bool *> mem.calloc(n*n, sizeof(bool))
    cdef int * qmax = <int *> mem.allocarray(n, sizeof(int))

    c[0] = c0

    # Defines c
    cdef int i, ui, vi
    for i in range(1, n):
        c[i] = c[i - 1] + n

    # Defines the adjacency matrix
    cdef list vertices = G.vertices(sort=False)
    cdef dict vertex_to_id = {v: i for i, v in enumerate(vertices)}

    for u in G:
        ui = vertex_to_id[u]
        for v in G.neighbors(u):
            vi = vertex_to_id[v]
            c[ui][vi] = 1
            c[vi][ui] = 1

    # Calls the solver
    cdef int clique_number
    cdef Maxclique * C = new Maxclique(c, n)
    sig_on()
    C.mcqdyn(qmax, clique_number)
    sig_off()

    # Returns the answer
    cdef list answer = [vertices[qmax[i]] for i in range(clique_number)]
    del C

    return answer
