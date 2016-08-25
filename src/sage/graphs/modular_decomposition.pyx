r"""
Modular decomposition
"""

include "cysignals/memory.pxi"

from libc.string cimport memset

#####################################################
# The following code is mainly a Cythonized
# copy of code found in src/random.c and src/dm.c
#####################################################

############################
# CONSTANTS
############################

cdef int FEUILLE = 0  # the node is a leaf
cdef int UNKN = 1
cdef int MODULE = 2
cdef int ARTEFACT = 3
cdef int SERIE = 4    # series composition
cdef int PARALLELE = 5  # parallel composition
cdef int PREMIER = 6  # prime composition

cdef dict codes = {
    UNKN: "Node",
    MODULE: "Module",
    ARTEFACT: "Artefact",
    SERIE: "Serie",
    PARALLELE: "Parallel",
    PREMIER: "Prime"
    }

cpdef modular_decomposition(g):
    r"""
    Returns a modular decomposition of the given graph.

    INPUT:

    - ``g`` -- a graph

    OUTPUT:

    A pair of two values (recursively encoding the decomposition) :

        * The type of the current module :

            * ``"Parallel"``
            * ``"Prime"``
            * ``"Serie"``

        * The list of submodules (as list of pairs ``(type, list)``,
          recursively...)  or the vertex's name if the module is a
          singleton.

    .. NOTE::

        As this function could be used by efficient C routines, the
        vertices returned are not labels but identifiants from ``[0,
        ..., g.order()-1]``

    ALGORITHM:

    This function uses a C implementation of a 2-step algorithm
    implemented by Fabien de Montgolfier [FMDecb]_ :

        * Computation of a factorizing permutation [HabibViennot1999b]_.

        * Computation of the tree itself [CapHabMont02b]_.

    EXAMPLES:

    The Bull Graph is prime::

        sage: from sage.graphs.modular_decomposition import modular_decomposition # optional -- modular_decomposition
        sage: modular_decomposition(graphs.BullGraph()) # optional -- modular_decomposition
        ('Prime', [3, 4, 0, 1, 2])

    The Petersen Graph too::

        sage: modular_decomposition(graphs.PetersenGraph()) # optional -- modular_decomposition
        ('Prime', [2, 6, 3, 9, 7, 8, 0, 1, 5, 4])

    This a clique on 5 vertices with 2 pendant edges, though, has a more
    interesting decomposition ::

        sage: g = graphs.CompleteGraph(5)
        sage: g.add_edge(0,5)
        sage: g.add_edge(0,6)
        sage: modular_decomposition(g) # optional -- modular_decomposition
        ('Serie', [0, ('Parallel', [5, ('Serie', [1, 4, 3, 2]), 6])])


    REFERENCES:

    .. [FMDecb] Fabien de Montgolfier
      http://www.liafa.jussieu.fr/~fm/algos/index.html

    .. [HabibViennot1999b] Michel Habib, Christiphe Paul, Laurent Viennot
      Partition refinement techniques: An interesting algorithmic tool kit
      International Journal of Foundations of Computer Science
      vol. 10 n2 pp.147--170, 1999

    .. [CapHabMont02b] \C. Capelle, M. Habib et F. de Montgolfier
      Graph decomposition and Factorising Permutations
      Discrete Mathematics and Theoretical Computer Sciences, vol 5 no. 1 , 2002.
    """
    cdef c_graphe G
    cdef c_adj * a
    cdef c_noeud * R

    cdef dict label_id = {}

    cdef int i
    cdef int j

    # Building the dictionary associating id to labels
    for id,label in enumerate(g.vertices()):
        label_id[label] = id

    G.n = g.order()
    G.G = <c_adj **> sig_malloc(G.n*sizeof(c_adj *))

    memset( G.G, 0, G.n*sizeof(c_adj *))

    # Creating the graph structure
    for u,v in g.edges(labels = False):

        i = label_id[u]
        j = label_id[v]

        a=<c_adj *> sig_malloc(sizeof(c_adj))
        a.s = j
        a.suiv = G.G[i]
        G.G[i] = a

        a=<c_adj *> sig_malloc(sizeof(c_adj))
        a.s = i
        a.suiv = G.G[j]
        G.G[j] = a

    R = decomposition_modulaire(G)

    return build_dict_from_decomposition(R)

cdef build_dict_from_decomposition(c_noeud * N):
    r"""
    Returns the dictionary corresponding to a decomposition
    """

    # recursively building the decomposition

    cdef c_fils *ffils
    cdef c_noeud *nfils
    cdef int i

    ffils = N.fils

    cdef list value = []

    while ffils != NULL:

        nfils = <c_noeud *> ffils.pointe
        if nfils.type == FEUILLE:
            value.append(nfils.nom)
        else:
            value.append(build_dict_from_decomposition(nfils))

        ffils = ffils.suiv

    return (codes[N.type], value)
