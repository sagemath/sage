# sage_setup: distribution = sagemath-graphs
r"""
Graphs with a given degree sequence

The methods defined here appear in :mod:`sage.graphs.graph_generators`.
"""

# ****************************************************************************
#       Copyright (C) 2006 Robert L. Miller <rlmillster@gmail.com>
#                          Emily A. Kirkman
#                     2009 Michael C. Yurko <myurko@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import sys
from sage.graphs.graph import Graph
from sage.misc.randstate import current_randstate


def DegreeSequence(deg_sequence):
    """
    Return a graph with the given degree sequence.

    This method raises a NetworkX error if the proposed degree sequence cannot
    be that of a graph.

    Graph returned is the one returned by the Havel-Hakimi algorithm, which
    constructs a simple graph by connecting vertices of highest degree to other
    vertices of highest degree, resorting the remaining vertices by degree and
    repeating the process. See Theorem 1.4 in [CL1996]_.

    INPUT:

    - ``deg_sequence`` -- list of integers with each entry corresponding to the
      degree of a different vertex

    EXAMPLES::

        sage: G = graphs.DegreeSequence([3,3,3,3])                                      # needs networkx
        sage: G.edges(sort=True, labels=False)                                          # needs networkx
        [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
        sage: G.show()                          # long time                             # needs networkx sage.plot

    ::

        sage: G = graphs.DegreeSequence([3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3])  # needs networkx
        sage: G.show()                          # long time                             # needs networkx sage.plot

    ::

        sage: G = graphs.DegreeSequence([4,4,4,4,4,4,4,4])                              # needs networkx
        sage: G.show()                          # long time                             # needs networkx sage.plot

    ::

        sage: G = graphs.DegreeSequence([1,2,3,4,3,4,3,2,3,2,1])                        # needs networkx
        sage: G.show()                          # long time                             # needs networkx sage.plot
    """
    import networkx
    return Graph(networkx.havel_hakimi_graph([int(i) for i in deg_sequence]))


def DegreeSequenceBipartite(s1, s2):
    r"""
    Return a bipartite graph whose two sets have the given degree sequences.

    Given two different sequences of degrees `s_1` and `s_2`, this functions
    returns ( if possible ) a bipartite graph on sets `A` and `B` such that the
    vertices in `A` have `s_1` as their degree sequence, while `s_2` is the
    degree sequence of the vertices in `B`.

    INPUT:

    - ``s_1`` -- list of integers corresponding to the degree sequence of the
      first set of vertices

    - ``s_2`` -- list of integers corresponding to the degree sequence of the
      second set of vertices

    ALGORITHM:

    This function works through the computation of the matrix given by the
    Gale-Ryser theorem, which is in this case the adjacency matrix of the
    bipartite graph.

    EXAMPLES:

    If we are given as sequences ``[2,2,2,2,2]`` and ``[5,5]`` we are given as
    expected the complete bipartite graph `K_{2,5}`::

        sage: g = graphs.DegreeSequenceBipartite([2,2,2,2,2],[5,5])                     # needs sage.combinat sage.modules
        sage: g.is_isomorphic(graphs.CompleteBipartiteGraph(5,2))                       # needs sage.combinat sage.modules
        True

    Some sequences being incompatible if, for example, their sums are different,
    the functions raises a :exc:`ValueError` when no graph corresponding
    to the degree sequences exists::

        sage: g = graphs.DegreeSequenceBipartite([2,2,2,2,1],[5,5])                     # needs sage.combinat sage.modules
        Traceback (most recent call last):
        ...
        ValueError: there exists no bipartite graph corresponding to the given degree sequences

    TESTS:

    :issue:`12155`::

        sage: graphs.DegreeSequenceBipartite([2,2,2,2,2],[5,5]).complement()            # needs sage.combinat sage.modules
        Graph on 7 vertices
    """
    from sage.combinat.integer_vector import gale_ryser_theorem
    from sage.graphs.bipartite_graph import BipartiteGraph

    s1 = sorted(s1, reverse=True)
    s2 = sorted(s2, reverse=True)

    m = gale_ryser_theorem(s1, s2)

    if m is False:
        raise ValueError("there exists no bipartite graph corresponding to "
                         "the given degree sequences")
    return Graph(BipartiteGraph(m))


def DegreeSequenceConfigurationModel(deg_sequence, seed=None):
    """
    Return a random pseudograph with the given degree sequence.

    This method raises a NetworkX error if the proposed degree sequence cannot
    be that of a graph with multiple edges and loops.

    One requirement is that the sum of the degrees must be even, since every
    edge must be incident with two vertices.

    INPUT:

    - ``deg_sequence`` -- list of integers with each entry corresponding to the
      expected degree of a different vertex

    - ``seed`` -- (optional) a ``random.Random`` seed or a Python ``int`` for
      the random number generator

    EXAMPLES::

        sage: G = graphs.DegreeSequenceConfigurationModel([1,1])                        # needs networkx
        sage: G.adjacency_matrix()                                                      # needs networkx sage.modules
        [0 1]
        [1 0]

    The output is allowed to contain both loops and multiple edges::

        sage: # needs networkx
        sage: deg_sequence = [3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3]
        sage: G = graphs.DegreeSequenceConfigurationModel(deg_sequence)
        sage: G.order(), G.size()
        (20, 30)
        sage: G.has_loops() or G.has_multiple_edges()  # random
        True
        sage: G.show()                          # long time                             # needs sage.plot

    REFERENCE:

    [New2003]_
    """
    if seed is None:
        seed = int(current_randstate().long_seed() % sys.maxsize)
    import networkx
    deg_sequence = [int(i) for i in deg_sequence]
    return Graph(networkx.configuration_model(deg_sequence, seed=seed),
                 loops=True, multiedges=True, sparse=True)


def DegreeSequenceTree(deg_sequence):
    """
    Return a tree with the given degree sequence.

    This method raises a NetworkX error if the proposed degree sequence cannot
    be that of a tree.

    Since every tree has one more vertex than edge, the degree sequence
    must satisfy ``len(deg_sequence) - sum(deg_sequence)/2 == 1``.

    INPUT:

    - ``deg_sequence`` -- list of integers with each entry corresponding to the
      expected degree of a different vertex

    EXAMPLES::

        sage: G = graphs.DegreeSequenceTree([3,1,3,3,1,1,1,2,1]); G                     # needs networkx
        Graph on 9 vertices
        sage: G.show()                          # long time                             # needs networkx sage.plot
    """
    import networkx
    return Graph(networkx.degree_sequence_tree([int(i) for i in deg_sequence]))


def DegreeSequenceExpected(deg_sequence, seed=None):
    """
    Return a random graph with expected given degree sequence.

    This method raises a NetworkX error if the proposed degree sequence cannot
    be that of a graph.

    One requirement is that the sum of the degrees must be even, since every
    edge must be incident with two vertices.

    INPUT:

    - ``deg_sequence`` -- list of integers with each entry corresponding to the
      expected degree of a different vertex

    - ``seed`` -- (optional) a ``random.Random`` seed or a Python ``int`` for
      the random number generator

    EXAMPLES::

        sage: G = graphs.DegreeSequenceExpected([1,2,3,2,3]); G                         # needs networkx
        Looped graph on 5 vertices
        sage: G.show()                          # long time                             # needs networkx sage.plot

    REFERENCE:

    [CL2002]_
    """
    if seed is None:
        seed = int(current_randstate().long_seed() % sys.maxsize)
    import networkx
    deg_sequence = [int(i) for i in deg_sequence]
    return Graph(networkx.expected_degree_graph(deg_sequence, seed=seed), loops=True)
