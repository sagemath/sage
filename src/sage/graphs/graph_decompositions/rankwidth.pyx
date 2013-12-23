r"""
Rank Decompositions of graphs

This modules wraps a C code from Philipp Klaus Krause computing a an optimal
rank-decomposition [RWKlause]_.

**Definitions :**

Given a graph `G` and a subset `S\subseteq V(G)` of vertices, the *rank-width*
of `S` in `G`, denoted `rw_G(S)`, is equal to the rank in `GF(2)` of the `|S|
\times (|V|-|S|)` matrix of the adjacencies between the vertices of `S` and
`V\backslash S`. By definition, `rw_G(S)` is qual to `rw_G(\overline S)` where
`\overline S` is the complement of `S` in `V(G)`.

A *rank-decomposition* of `G` is a tree whose `n` leaves are the elements of
`V(G)`, and whose internal noes have degree 3. In a tree, ay edge naturally
corresponds to a bipartition of the vertex set : indeed, the reoal of any edge
splits the tree into two connected components, thus splitting the set of leaves
(i.e. vertices of `G`) into two sets. Hence we can define for any edge `e\in
E(G)` a width equal to the value `rw_G(S)` or `rw_G(\overline S)`, where
`S,\overline S` is the bipartition obtained from `e`. The *rank-width*
associated to the whole decomposition is then set to the maximum of the width of
all the edges it contains.

A *rank-decomposition* is said to be optimal for `G` if it is the decomposition
achieving the minimal *rank-width*.

**RW -- The original source code :**

RW [RWKlause]_ is a program that calculates rank-width and
rank-decompositions. It is based on ideas from :

    * "Computing rank-width exactly" by Sang-il Oum [Oum]_
    * "Sopra una formula numerica" by Ernesto Pascal
    * "Generation of a Vector from the Lexicographical Index" by B.P. Buckles and M. Lybanon [BL]_
    * "Fast additions on masked integers" by Michael D. Adams and David S. Wise [AW]_

**OUTPUT:**

The rank decomposition is returned as a tree whose vertices are subsets of
`V(G)`. Its leaves, corresponding to the vertices of `G` are sets of 1 elements,
i.e. singletons.

::

    sage: g = graphs.PetersenGraph()
    sage: rw, tree = g.rank_decomposition()
    sage: all(len(v)==1 for v in tree if tree.degree(v) == 1)
    True

The internal nodes are sets of the decomposition. This way, it is easy to deduce
the bipartition associated to an edge from the tree. Indeed, two adjacent
vertices of the tree are comarable sets : they yield the bipartition obtained
from the smaller of the two and its complement.

::

    sage: g = graphs.PetersenGraph()
    sage: rw, tree = g.rank_decomposition()
    sage: u = Set([8, 9, 3, 7])
    sage: v = Set([8, 9])
    sage: tree.has_edge(u,v)
    True
    sage: m = min(u,v)
    sage: bipartition = (m, Set(g.vertices()) - m)
    sage: bipartition
    ({8, 9}, {0, 1, 2, 3, 4, 5, 6, 7})

.. WARNING::

    * The current implementation cannot handle graphs of `\geq 32` vertices.

    * A bug that has been reported upstream make the code crash immediately on
      instances of size `30`. If you experience this kind of bug please report
      it to us, what we need is some information on the hardware you run to know
      where it comes from !

EXAMPLE::

        sage: g = graphs.PetersenGraph()
        sage: g.rank_decomposition()
        (3, Graph on 19 vertices)

AUTHORS:

- Philipp Klaus Krause : Implementation of the C algorithm [RWKlause]_.
- Nathann Cohen : Interface with Sage and documentation.

REFERENCES:

  .. [RWKlause] Philipp Klaus Krause -- rw v0.2
    http://pholia.tdi.informatik.uni-frankfurt.de/~philipp/software/rw.shtml

  .. [Oum] Sang-il Oum
    Computing rank-width exactly
    Information Processing Letters, 2008
    vol. 109, n. 13, p. 745--748
    Elsevier
    http://mathsci.kaist.ac.kr/~sangil/pdf/2008exp.pdf

  .. [BL] Buckles, B.P. and Lybanon, M.
    Algorithm 515: generation of a vector from the lexicographical index
    ACM Transactions on Mathematical Software (TOMS), 1977
    vol. 3, n. 2, pages 180--182
    ACM

  .. [AW] Adams, M.D. and Wise, D.S.
    Fast additions on masked integers
    ACM SIGPLAN Notices, 2006
    vol. 41, n.5, pages 39--45
    ACM
    http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.86.1801&rep=rep1&type=pdf

Methods
-------
"""

#*****************************************************************************
#      Copyright (C) 2011 Nathann Cohen <nathann.cohen@gail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************


include 'sage/ext/stdsage.pxi'
include 'sage/ext/cdefs.pxi'
include "sage/ext/interrupt.pxi"

cdef list id_to_vertices
cdef dict vertices_to_id

def rank_decomposition(G, verbose = False):
    r"""
    Computes an optimal rank-decomposition of the given graph.

    This function is available as a method of the :class:`Graph
    <sage.graphs.graph>` class. See :meth:`rank_decomposition
    <sage.graphs.graph.Graph.rank_decomposition>`.

    INPUT:

    - ``verbose`` (boolean) -- whether to display progress information while
      computing the decomposition.

    OUTPUT:

    A pair ``(rankwidth, decomposition_tree)``, where ``rankwidth`` is a
    numerical value and ``decomposition_tree`` is a ternary tree describing the
    decomposition (cf. the module's documentation).

    EXAMPLE::

        sage: from sage.graphs.graph_decompositions.rankwidth import rank_decomposition
        sage: g = graphs.PetersenGraph()
        sage: rank_decomposition(g)
        (3, Graph on 19 vertices)

    On more than 32 vertices::

        sage: g = graphs.RandomGNP(40, .5)
        sage: rank_decomposition(g)
        Traceback (most recent call last):
        ...
        RuntimeError: the rank decomposition cannot be computed on graphs of >= 32 vertices

    The empty graph::

        sage: g = Graph()
        sage: rank_decomposition(g)
        (0, Graph on 0 vertices)
    """
    cdef int n = G.order()

    if n >= 32:
        raise RuntimeError("the rank decomposition cannot be computed "+
                        "on graphs of >= 32 vertices")

    elif n == 0:
        from sage.graphs.graph import Graph
        return (0, Graph())

    global num_vertices

    cdef int i

    if sage_graph_to_matrix(G):
        raise RuntimeError("There has been a mistake while converting the Sage "+
                        "graph to a C structure. The memory is probably "+
                        "insufficient (2^(n+1) is a *LOT*).")

    for 0 <= i < n+1:

        if verbose:
            print "Calculating for subsets of size ", i, "/",(n+1)

        # We want to properly deal with exceptions, in particular
        # KeyboardInterrupt. Whatever happens, when this code fails the memory
        # *HAS* to be freed, as it is particularly greedy (a graph of 29
        # vertices already requires more than 1GB of memory).

        if not sig_on_no_except():
            destroy_rw()
        cython_check_exception()

        # Actual computation
        calculate_level(i)
        sig_off()

    cdef int rank_width = <int> get_rw()

    #Original way of displaying the decomposition
    #print_rank_dec(0x7ffffffful >> (31 - num_vertices), 0)
    g = mkgraph()

    # Free the memory
    destroy_rw()

    return (rank_width, g)

cdef int sage_graph_to_matrix(G):
    r"""
    Converts the given Sage graph as an adjacency matrix.
    """
    global id_to_vertices
    global vertices_to_id
    global adjacency_matrix
    global slots
    global cslots

    id_to_vertices = []
    vertices_to_id = {}

    global num_vertices
    num_vertices = G.order()

    cdef int i,j

    # Prepares the C structure for the computation
    if init_rw_dec(num_vertices):
        # The original code does not completely frees the memory in case of
        # error
        if cslots != NULL:
            sage_free(cslots)
        if slots != NULL:
            sage_free(slots)
        if adjacency_matrix != NULL:
            sage_free(adjacency_matrix)
        return 1

    memset(adjacency_matrix, 0, sizeof(subset_t) * num_vertices)

    # Initializing the lists of vertices
    for i,v in enumerate(G.vertices()):
        id_to_vertices.append(v)
        vertices_to_id[v] = i

    # Filling the matrix
    for u,v in G.edges(labels = False):
        if u==v:
            continue
        set_am(vertices_to_id[u], vertices_to_id[v], 1)

    # All is fine.
    return 0

cdef uint_fast32_t bitmask(int i):
    return(1ul << i)

cdef void set_am(int i, int j, int val):
    r"""
    Set/Unset an arc between vertices i and j

    (this function is a copy of what can be found in rankwidth/rw.c)
    """
    global adjacency_matrix

    adjacency_matrix[i] &= ~bitmask(j)
    adjacency_matrix[j] &= ~bitmask(i)

    if val:
        adjacency_matrix[i] |= bitmask(j)
        adjacency_matrix[j] |= bitmask(i)

cdef void print_rank_dec(subset_t s, int l):
    r"""
    Prints the current rank decomposition as a text

    This function is a copy of the C routine printing the rank-decomposition is
    the original source code. It s not used at the moment, but can still prove
    useful when debugging the code.
    """
    global cslots

    print ('\t'*l),

    print "cslot: ", <unsigned int> s
    if cslots[s] == 0:
        return
    print_rank_dec(cslots[s], l + 1)
    print_rank_dec(s & ~cslots[s], l + 1)

def mkgraph():
    r"""
    Returns the graph corresponding the the current rank-decomposition.

    (This function is for internal use)

    EXAMPLE::

        sage: from sage.graphs.graph_decompositions.rankwidth import rank_decomposition
        sage: g = graphs.PetersenGraph()
        sage: rank_decomposition(g)
        (3, Graph on 19 vertices)
    """
    global cslots
    global num_vertices
    global id_to_vertices

    cdef subset_t s

    from sage.graphs.graph import Graph
    g = Graph()

    cdef subset_t * tab = <subset_t *> sage_malloc(sizeof(subset_t) * (2*num_vertices -1))
    tab[0] = 0x7ffffffful >> (31 - num_vertices)

    cdef int beg = 0
    cdef int end = 1

    while beg != end:

        s = tab[beg]
        beg += 1

        if cslots[s] == 0:
            continue

        g.add_edge(bitset_to_vertex_set(s), bitset_to_vertex_set(s&~cslots[s]))
        g.add_edge(bitset_to_vertex_set(s), bitset_to_vertex_set(cslots[s]))

        tab[end] = s&~cslots[s]
        end += 1
        tab[end] = cslots[s]
        end += 1

    sage_free(tab)
    return g

cdef bitset_to_vertex_set(subset_t s):
    """
    Returns as a Set object the set corresponding to the given subset_t
    variable.
    """
    from sage.rings.integer import Integer
    from sage.sets.set import Set
    from sage.misc.bitset import FrozenBitset
    return Set(list(FrozenBitset((Integer(<unsigned int> s).binary())[::-1])))
