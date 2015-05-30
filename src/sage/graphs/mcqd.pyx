include "sage/ext/interrupt.pxi"
include 'sage/ext/stdsage.pxi'
include 'sage/ext/cdefs.pxi'

def mcqd(G):
    """
    Computes the max clique using MCQD

    INPUT:

    - ``G`` -- A graph

    TESTS::

        sage: from sage.graphs.mcqd import mcqd         # optional - mcqd
        sage: for i in range(10):                       # optional - mcqd
        ...       g = graphs.RandomGNP(15,.5)           # optional - mcqd
        ...       if g.clique_number() != len(mcqd(g)): # optional - mcqd
        ...           print "This is dead wrong !"      # optional - mcqd
    """
    cdef int n = G.order()

    # - c0 is the adjacency matrix
    # - c points toward each row of the matrix
    # - qmax stores the max clique
    cdef bool ** c = <bool **> sage_malloc(n*sizeof(bool *))
    cdef bool * c0 = <bool *> sage_calloc(n*n,sizeof(bool))
    cdef int * qmax = <int *> sage_malloc(n*sizeof(int))
    sage_free(NULL)
    if c == NULL or c0 == NULL or qmax == NULL:
        sage_free(c)
        sage_free(c0)
        sage_free(qmax)
        raise MemoryError("Allocation Failed")

    c[0] = c0

    # Defines c
    cdef int i,ui,vi
    for i in range(1,n):
        c[i] = c[i-1] + n

    # Defines the adjacency matrix
    cdef list vertices = G.vertices()
    cdef dict vertex_to_id = {v:i for i,v in enumerate(vertices)}

    for u in G:
        ui = vertex_to_id[u]
        for v in G.neighbors(u):
            vi = vertex_to_id[v]
            c[ui][vi] = 1
            c[vi][ui] = 1

    # Calls the solver
    cdef int clique_number
    cdef Maxclique * C = new Maxclique(c,n)
    C.mcqdyn(qmax, clique_number)

    # Returns the answer
    cdef list answer = [vertices[qmax[i]] for i in range(clique_number)]
    sage_free(c[0])
    sage_free(c)
    sage_free(qmax)
    del C

    return answer
