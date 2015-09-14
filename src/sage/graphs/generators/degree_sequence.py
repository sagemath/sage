# -*- coding: utf-8 -*-
r"""
Graphs with a given degree sequence

The methods defined here appear in :mod:`sage.graphs.graph_generators`.
"""

###########################################################################
#
#           Copyright (C) 2006 Robert L. Miller <rlmillster@gmail.com>
#                              and Emily A. Kirkman
#           Copyright (C) 2009 Michael C. Yurko <myurko@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
###########################################################################

# import from Sage library
from sage.graphs.graph import Graph
from sage.misc.randstate import current_randstate


def DegreeSequence(deg_sequence):
    """
    Returns a graph with the given degree sequence. Raises a NetworkX
    error if the proposed degree sequence cannot be that of a graph.

    Graph returned is the one returned by the Havel-Hakimi algorithm,
    which constructs a simple graph by connecting vertices of highest
    degree to other vertices of highest degree, resorting the remaining
    vertices by degree and repeating the process. See Theorem 1.4 in
    [CharLes1996]_.

    INPUT:

    -  ``deg_sequence`` - a list of integers with each
       entry corresponding to the degree of a different vertex.


    EXAMPLES::

        sage: G = graphs.DegreeSequence([3,3,3,3])
        sage: G.edges(labels=False)
        [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
        sage: G.show()  # long time

    ::

        sage: G = graphs.DegreeSequence([3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3])
        sage: G.show()  # long time

    ::

        sage: G = graphs.DegreeSequence([4,4,4,4,4,4,4,4])
        sage: G.show()  # long time

    ::

        sage: G = graphs.DegreeSequence([1,2,3,4,3,4,3,2,3,2,1])
        sage: G.show()  # long time

    REFERENCE:

    .. [CharLes1996] Chartrand, G. and Lesniak, L.: Graphs and Digraphs.
      Chapman and Hall/CRC, 1996.
    """
    import networkx
    return Graph(networkx.havel_hakimi_graph([int(i) for i in deg_sequence]))

def DegreeSequenceBipartite(s1 ,s2 ):
    r"""
    Returns a bipartite graph whose two sets have the given
    degree sequences.

    Given two different sequences of degrees `s_1` and `s_2`,
    this functions returns ( if possible ) a bipartite graph
    on sets `A` and `B` such that the vertices in `A` have
    `s_1` as their degree sequence, while `s_2` is the degree
    sequence of the vertices in `B`.

    INPUT:

    - ``s_1`` -- list of integers corresponding to the degree
      sequence of the first set.
    - ``s_2`` -- list of integers corresponding to the degree
      sequence of the second set.

    ALGORITHM:

    This function works through the computation of the matrix
    given by the Gale-Ryser theorem, which is in this case
    the adjacency matrix of the bipartite graph.

    EXAMPLES:

    If we are given as sequences ``[2,2,2,2,2]`` and ``[5,5]``
    we are given as expected the complete bipartite
    graph `K_{2,5}` ::

        sage: g = graphs.DegreeSequenceBipartite([2,2,2,2,2],[5,5])
        sage: g.is_isomorphic(graphs.CompleteBipartiteGraph(5,2))
        True

    Some sequences being incompatible if, for example, their sums
    are different, the functions raises a ``ValueError`` when no
    graph corresponding to the degree sequences exists. ::

        sage: g = graphs.DegreeSequenceBipartite([2,2,2,2,1],[5,5])
        Traceback (most recent call last):
        ...
        ValueError: There exists no bipartite graph corresponding to the given degree sequences

    TESTS:

    Trac ticket #12155::

        sage: graphs.DegreeSequenceBipartite([2,2,2,2,2],[5,5]).complement()
        Graph on 7 vertices
    """

    from sage.combinat.integer_vector import gale_ryser_theorem
    from sage.graphs.bipartite_graph import BipartiteGraph

    s1 = sorted(s1, reverse = True)
    s2 = sorted(s2, reverse = True)

    m = gale_ryser_theorem(s1,s2)

    if m is False:
        raise ValueError("There exists no bipartite graph corresponding to the given degree sequences")
    else:
        return Graph(BipartiteGraph(m))

def DegreeSequenceConfigurationModel(deg_sequence, seed=None):
    """
    Returns a random pseudograph with the given degree sequence. Raises
    a NetworkX error if the proposed degree sequence cannot be that of
    a graph with multiple edges and loops.

    One requirement is that the sum of the degrees must be even, since
    every edge must be incident with two vertices.

    INPUT:

    -  ``deg_sequence`` - a list of integers with each
       entry corresponding to the expected degree of a different vertex.

    -  ``seed`` - for the random number generator.


    EXAMPLES::

        sage: G = graphs.DegreeSequenceConfigurationModel([1,1])
        sage: G.adjacency_matrix()
        [0 1]
        [1 0]

    Note: as of this writing, plotting of loops and multiple edges is
    not supported, and the output is allowed to contain both types of
    edges.

    ::

        sage: G = graphs.DegreeSequenceConfigurationModel([3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3])
        sage: G.edges(labels=False)
        [(0, 2), (0, 10), (0, 15), (1, 6), (1, 16), (1, 17), (2, 5), (2, 19), (3, 7), (3, 14), (3, 14), (4, 9), (4, 13), (4, 19), (5, 6), (5, 15), (6, 11), (7, 11), (7, 17), (8, 11), (8, 18), (8, 19), (9, 12), (9, 13), (10, 15), (10, 18), (12, 13), (12, 16), (14, 17), (16, 18)]
        sage: G.show()  # long time

    REFERENCE:

    .. [Newman2003] Newman, M.E.J. The Structure and function of complex
      networks, SIAM Review vol. 45, no. 2 (2003), pp. 167-256.
    """
    if seed is None:
        seed = current_randstate().long_seed()
    import networkx
    return Graph(networkx.configuration_model([int(i) for i in deg_sequence], seed=seed), loops=True, multiedges=True, sparse=True)

def DegreeSequenceTree(deg_sequence):
    """
    Returns a tree with the given degree sequence. Raises a NetworkX
    error if the proposed degree sequence cannot be that of a tree.

    Since every tree has one more vertex than edge, the degree sequence
    must satisfy len(deg_sequence) - sum(deg_sequence)/2 == 1.

    INPUT:

    -  ``deg_sequence`` - a list of integers with each
       entry corresponding to the expected degree of a different vertex.


    EXAMPLE::

        sage: G = graphs.DegreeSequenceTree([3,1,3,3,1,1,1,2,1])
        sage: G.show()  # long time
    """
    import networkx
    return Graph(networkx.degree_sequence_tree([int(i) for i in deg_sequence]))

def DegreeSequenceExpected(deg_sequence, seed=None):
    """
    Returns a random graph with expected given degree sequence. Raises
    a NetworkX error if the proposed degree sequence cannot be that of
    a graph.

    One requirement is that the sum of the degrees must be even, since
    every edge must be incident with two vertices.

    INPUT:

    -  ``deg_sequence`` - a list of integers with each
       entry corresponding to the expected degree of a different vertex.

    -  ``seed`` - for the random number generator.


    EXAMPLE::

        sage: G = graphs.DegreeSequenceExpected([1,2,3,2,3])
        sage: G.edges(labels=False)
        [(0, 2), (0, 3), (1, 1), (1, 4), (2, 3), (2, 4), (3, 4), (4, 4)]
        sage: G.show()  # long time

    REFERENCE:

    .. [ChungLu2002] Chung, Fan and Lu, L. Connected components in random
      graphs with given expected degree sequences.
      Ann. Combinatorics (6), 2002 pp. 125-145.
    """
    if seed is None:
        seed = current_randstate().long_seed()
    import networkx
    return Graph(networkx.expected_degree_graph([int(i) for i in deg_sequence], seed=seed), loops=True)
