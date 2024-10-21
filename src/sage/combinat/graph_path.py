# sage_setup: distribution = sagemath-graphs
r"""
Paths in Directed Acyclic Graphs
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.structure.parent import Parent
import sage.graphs.digraph as digraph


def GraphPaths(g, source=None, target=None):
    """
    Return the combinatorial class of paths in the directed acyclic graph g.

    EXAMPLES::

        sage: G = DiGraph({1:[2,2,3], 2:[3,4], 3:[4], 4:[5,5]}, multiedges=True)

    If source and target are not given, then the returned class
    contains all paths (including trivial paths containing only one
    vertex).

    ::

        sage: p = GraphPaths(G); p
        Paths in Multi-digraph on 5 vertices
        sage: p.cardinality()
        37
        sage: path = p.random_element()
        sage: all(G.has_edge(*path[i:i+2]) for i in range(len(path) -1))
        True

    If the source is specified, then the returned class contains all of
    the paths starting at the vertex source (including the trivial
    path).

    ::

        sage: p = GraphPaths(G, source=3); p
        Paths in Multi-digraph on 5 vertices starting at 3
        sage: p.list()
        [[3], [3, 4], [3, 4, 5], [3, 4, 5]]

    If the target is specified, then the returned class contains all of
    the paths ending at the vertex target (including the trivial
    path).

    ::

        sage: p = GraphPaths(G, target=3); p
        Paths in Multi-digraph on 5 vertices ending at 3
        sage: p.cardinality()
        5
        sage: p.list()
        [[3], [1, 3], [2, 3], [1, 2, 3], [1, 2, 3]]

    If both the target and source are specified, then the returned
    class contains all of the paths from source to target.

    ::

        sage: p = GraphPaths(G, source=1, target=3); p
        Paths in Multi-digraph on 5 vertices starting at 1 and ending at 3
        sage: p.cardinality()
        3
        sage: p.list()
        [[1, 2, 3], [1, 2, 3], [1, 3]]

    Note that G must be a directed acyclic graph.

    ::

        sage: G = DiGraph({1:[2,2,3,5], 2:[3,4], 3:[4], 4:[2,5,7], 5:[6]}, multiedges=True)
        sage: GraphPaths(G)
        Traceback (most recent call last):
        ...
        TypeError: g must be a directed acyclic graph
    """
    if not isinstance(g, digraph.DiGraph):
        raise TypeError("g must be a DiGraph")
    elif not g.is_directed_acyclic():
        raise TypeError("g must be a directed acyclic graph")

    if source is None and target is None:
        return GraphPaths_all(g)
    elif source is not None and target is None:
        if source not in g:
            raise ValueError("source must be in g")
        return GraphPaths_s(g, source)
    elif source is None and target is not None:
        if target not in g:
            raise ValueError("target must be in g")
        return GraphPaths_t(g, target)
    else:
        if source not in g:
            raise ValueError("source must be in g")
        if target not in g:
            raise ValueError("target must be in g")
        return GraphPaths_st(g, source, target)


class GraphPaths_common:
    def __eq__(self, other):
        """
        Test for equality.

        EXAMPLES::

            sage: G1 = DiGraph({1:[2,3], 2:[3]})
            sage: p1 = GraphPaths(G1)
            sage: G2 = DiGraph({2:[3], 3:[4]})
            sage: p2 = GraphPaths(G2)
            sage: p1 == p1
            True
            sage: p1 == p2
            False
        """
        if not isinstance(other, GraphPaths_common):
            return False
        return self.graph == other.graph

    def __ne__(self, other):
        """
        Test for unequality.

        EXAMPLES::

            sage: G1 = DiGraph({1:[2,3], 2:[3]})
            sage: p1 = GraphPaths(G1)
            sage: G2 = DiGraph({2:[3], 3:[4]})
            sage: p2 = GraphPaths(G2)
            sage: p1 != p2
            True
            sage: p1 != p1
            False
        """
        return not (self == other)

    def outgoing_edges(self, v):
        """
        Return a list of v's outgoing edges.

        EXAMPLES::

            sage: G = DiGraph({1:[2,2,3], 2:[3,4], 3:[4], 4:[5,5]}, multiedges=True)
            sage: p = GraphPaths(G)
            sage: p.outgoing_edges(2)
            [(2, 3, None), (2, 4, None)]
        """
        return list(self.graph.outgoing_edge_iterator(v))

    def incoming_edges(self, v):
        """
        Return a list of v's incoming edges.

        EXAMPLES::

            sage: G = DiGraph({1:[2,2,3], 2:[3,4], 3:[4], 4:[5,5]}, multiedges=True)
            sage: p = GraphPaths(G)
            sage: p.incoming_edges(2)
            [(1, 2, None), (1, 2, None)]
        """
        return list(self.graph.incoming_edge_iterator(v))

    def outgoing_paths(self, v):
        """
        Return a list of the paths that start at v.

        EXAMPLES::

            sage: G = DiGraph({1:[2,2,3], 2:[3,4], 3:[4], 4:[5,5]}, multiedges=True)
            sage: gp = GraphPaths(G)
            sage: gp.outgoing_paths(3)
            [[3], [3, 4], [3, 4, 5], [3, 4, 5]]
            sage: gp.outgoing_paths(2)
            [[2],
             [2, 3],
             [2, 3, 4],
             [2, 3, 4, 5],
             [2, 3, 4, 5],
             [2, 4],
             [2, 4, 5],
             [2, 4, 5]]
        """
        source_paths = [[v]]
        for e in self.outgoing_edges(v):
            target = e[1]
            target_paths = self.outgoing_paths(target)
            target_paths = [[v] + path for path in target_paths]

            source_paths += target_paths

        return source_paths

    def incoming_paths(self, v):
        """
        Return a list of paths that end at v.

        EXAMPLES::

            sage: G = DiGraph({1:[2,2,3], 2:[3,4], 3:[4], 4:[5,5]}, multiedges=True)
            sage: gp = GraphPaths(G)
            sage: gp.incoming_paths(2)
            [[2], [1, 2], [1, 2]]
        """
        target_paths = [[v]]
        for e in self.incoming_edges(v):
            source = e[0]
            source_paths = self.incoming_paths(source)
            source_paths = [path + [v] for path in source_paths]
            target_paths += source_paths
        return target_paths

    def paths_from_source_to_target(self, source, target):
        """
        Return a list of paths from source to target.

        EXAMPLES::

            sage: G = DiGraph({1:[2,2,3], 2:[3,4], 3:[4], 4:[5,5]}, multiedges=True)
            sage: gp = GraphPaths(G)
            sage: gp.paths_from_source_to_target(2,4)
            [[2, 3, 4], [2, 4]]
        """
        source_paths = self.outgoing_paths(source)
        return [path for path in source_paths if path[-1] == target]

    def paths(self):
        """
        Return a list of all the paths of ``self``.

        EXAMPLES::

            sage: G = DiGraph({1:[2,2,3], 2:[3,4], 3:[4], 4:[5,5]}, multiedges=True)
            sage: gp = GraphPaths(G)
            sage: len(gp.paths())
            37
        """
        paths = []
        for source in self.graph.vertex_iterator():
            paths += self.outgoing_paths(source)
        return paths


class GraphPaths_all(Parent, GraphPaths_common):
    """
    EXAMPLES::

        sage: G = DiGraph({1:[2,2,3], 2:[3,4], 3:[4], 4:[5,5]}, multiedges=True)
        sage: p = GraphPaths(G)
        sage: p.cardinality()
        37
    """

    def __init__(self, g):
        """
        TESTS::

            sage: G = DiGraph({1:[2,2,3], 2:[3,4], 3:[4], 4:[5,5]}, multiedges=True)
            sage: p = GraphPaths(G)
            sage: p == loads(dumps(p))
            True
        """
        self.graph = g
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def __repr__(self):
        """
        TESTS::

            sage: G = DiGraph({1:[2,2,3], 2:[3,4], 3:[4], 4:[5,5]}, multiedges=True)
            sage: p = GraphPaths(G)
            sage: repr(p)
            'Paths in Multi-digraph on 5 vertices'
        """
        return "Paths in %s" % repr(self.graph)

    def list(self):
        """
        Return a list of the paths of ``self``.

        EXAMPLES::

            sage: G = DiGraph({1:[2,2,3], 2:[3,4], 3:[4], 4:[5,5]}, multiedges=True)
            sage: len(GraphPaths(G).list())
            37
        """
        return self.paths()


class GraphPaths_t(Parent, GraphPaths_common):
    def __init__(self, g, target):
        """
        TESTS::

            sage: G = DiGraph({1:[2,2,3], 2:[3,4], 3:[4], 4:[5,5]}, multiedges=True)
            sage: p = GraphPaths(G, target=4)
            sage: p == loads(dumps(p))
            True
        """
        self.graph = g
        self.target = target
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def __repr__(self):
        """
        TESTS::

            sage: G = DiGraph({1:[2,2,3], 2:[3,4], 3:[4], 4:[5,5]}, multiedges=True)
            sage: p = GraphPaths(G, target=4)
            sage: repr(p)
            'Paths in Multi-digraph on 5 vertices ending at 4'
        """
        return "Paths in %s ending at %s" % (repr(self.graph), self.target)

    def list(self):
        """
        EXAMPLES::

            sage: G = DiGraph({1:[2,2,3], 2:[3,4], 3:[4], 4:[5,5]}, multiedges=True)
            sage: p = GraphPaths(G, target=4)
            sage: p.list()
            [[4],
             [2, 4],
             [1, 2, 4],
             [1, 2, 4],
             [3, 4],
             [1, 3, 4],
             [2, 3, 4],
             [1, 2, 3, 4],
             [1, 2, 3, 4]]
        """
        return self.incoming_paths(self.target)


class GraphPaths_s(Parent, GraphPaths_common):
    def __init__(self, g, source):
        """
        TESTS::

            sage: G = DiGraph({1:[2,2,3], 2:[3,4], 3:[4], 4:[5,5]}, multiedges=True)
            sage: p = GraphPaths(G, 4)
            sage: p == loads(dumps(p))
            True
        """
        self.graph = g
        self.source = source
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def __repr__(self):
        """
        TESTS::

            sage: G = DiGraph({1:[2,2,3], 2:[3,4], 3:[4], 4:[5,5]}, multiedges=True)
            sage: p = GraphPaths(G, 4)
            sage: repr(p)
            'Paths in Multi-digraph on 5 vertices starting at 4'
        """
        return "Paths in %s starting at %s" % (repr(self.graph), self.source)

    def list(self):
        """
        EXAMPLES::

            sage: G = DiGraph({1:[2,2,3], 2:[3,4], 3:[4], 4:[5,5]}, multiedges=True)
            sage: p = GraphPaths(G, 4)
            sage: p.list()
            [[4], [4, 5], [4, 5]]
        """
        return self.outgoing_paths(self.source)


class GraphPaths_st(Parent, GraphPaths_common):
    """
    EXAMPLES::

        sage: G = DiGraph({1:[2,2,3], 2:[3,4], 3:[4], 4:[5,5]}, multiedges=True)
        sage: GraphPaths(G,1,2).cardinality()
        2
        sage: GraphPaths(G,1,3).cardinality()
        3
        sage: GraphPaths(G,1,4).cardinality()
        5
        sage: GraphPaths(G,1,5).cardinality()
        10
        sage: GraphPaths(G,2,3).cardinality()
        1
        sage: GraphPaths(G,2,4).cardinality()
        2
        sage: GraphPaths(G,2,5).cardinality()
        4
        sage: GraphPaths(G,3,4).cardinality()
        1
        sage: GraphPaths(G,3,5).cardinality()
        2
        sage: GraphPaths(G,4,5).cardinality()
        2
    """

    def __init__(self, g, source, target):
        """
        TESTS::

            sage: G = DiGraph({1:[2,2,3], 2:[3,4], 3:[4], 4:[5,5]}, multiedges=True)
            sage: p = GraphPaths(G,1,2)
            sage: p == loads(dumps(p))
            True
        """
        self.graph = g
        self.source = source
        self.target = target
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def __repr__(self):
        """
        TESTS::

            sage: G = DiGraph({1:[2,2,3], 2:[3,4], 3:[4], 4:[5,5]}, multiedges=True)
            sage: p = GraphPaths(G,1,2)
            sage: repr(p)
            'Paths in Multi-digraph on 5 vertices starting at 1 and ending at 2'
        """
        return "Paths in %s starting at %s and ending at %s" % (repr(self.graph), self.source, self.target)

    def list(self):
        """
        EXAMPLES::

            sage: G = DiGraph({1:[2,2,3], 2:[3,4], 3:[4], 4:[5,5]}, multiedges=True)
            sage: p = GraphPaths(G,1,2)
            sage: p.list()
            [[1, 2], [1, 2]]
        """
        return self.paths_from_source_to_target(self.source, self.target)
