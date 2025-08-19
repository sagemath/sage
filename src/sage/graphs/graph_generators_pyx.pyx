r"""
Common graphs and digraphs generators (Cython)

AUTHORS:

- David Coudert (2012)
"""

# #############################################################################
#           Copyright (C) 2012 David Coudert <david.coudert@inria.fr>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         https://www.gnu.org/licenses/
# #############################################################################

from sage.misc.randstate cimport random
from sage.misc.randstate import set_random_seed


def RandomGNP(n, p, bint directed=False, bint loops=False, seed=None,
              immutable=False):
    r"""
    Return a random graph or a digraph on `n` nodes.

    Each edge is inserted independently with probability `p`.

    INPUT:

    - ``n`` -- number of nodes of the digraph

    - ``p`` -- probability of an edge

    - ``directed`` -- boolean (default: ``False``); whether the random graph is
      directed or undirected (default)

    - ``loops`` -- boolean (default: ``False``); whether the random digraph may
      have loops or not. This value is used only when ``directed == True``

    - ``seed`` -- a ``random.Random`` seed or a Python ``int`` for the random
      number generator (default: ``None``)

    - ``immutable`` -- boolean (default: ``False``); whether to return an
      immutable or mutable (di)graph.

    REFERENCES:

    - [ER1959]_

    - [Gil1959]_

    EXAMPLES::

        sage: from sage.graphs.graph_generators_pyx import RandomGNP
        sage: D = RandomGNP(10, .2, directed=True, seed=0)
        sage: D.num_verts()
        10
        sage: D.edges(sort=True, labels=False)
        [(0, 3), (0, 6), (1, 7), (1, 9), (4, 6), (4, 7), (5, 4), (5, 6),
         (5, 8), (5, 9), (6, 3), (7, 2), (7, 9), (8, 5), (9, 1), (9, 5)]

    TESTS::

        sage: from numpy import mean                                                    # needs numpy
        sage: abs(mean([RandomGNP(200, .2).density() for i in range(30)]) - .2) < .001  # needs numpy
        ...True...
        sage: RandomGNP(150, .2, loops=True)
        Traceback (most recent call last):
        ...
        ValueError: parameter 'loops' can be set to True only when 'directed' is True
    """
    if seed is not None:
        set_random_seed(seed)

    # according the sage.misc.randstate.pyx documentation, random
    # integers are on 31 bits. We thus set the pivot value to p*2^31
    cdef float RAND_MAX_f = float(1<<31)
    cdef int pp = int(round(float(p * RAND_MAX_f)))

    if directed:
        from sage.graphs.digraph import DiGraph as GT
    else:
        if loops:
            raise ValueError("parameter 'loops' can be set to True only when 'directed' is True")
        from sage.graphs.graph import Graph as GT

    name = 'Random' + ('Directed' if directed else '') + 'GNP(%s,%s)' % (n, p)

    cdef int i, j
    edges = ((i, j) for i in range(n)
                 for j in range((0 if directed else i + 1), n)
                 if (i != j or loops) and random() < pp)

    return GT([range(n), edges], format='vertices_and_edges',
              loops=directed and loops, name=name, immutable=immutable)
