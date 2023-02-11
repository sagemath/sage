r"""
Acyclic Orientations

This module implements the :class:`AcyclicOrientations` class that can be used to enumerate acyclic orientations of a given graph.
Return an iterable over acyclic orientations: Each acyclic orientation is given as a list of edges `A` such that `DiGraph(A)` raises the DiGraph endowed with the wished acyclic orientation of ``self``.

Classes and methods
-------------------

"""

from sage.graphs.digraph import DiGraph
from sage.graphs.graph import Graph

class AcyclicOrientations():
    r"""
    Class of acyclic orientations of a (di)graph.
    
    Made for enumerating acyclic orientations of a given graph. Return an iterable over acyclic orientations: Each acyclic orientation is given as a list of edges `A` such that `DiGraph(A)` raises the DiGraph endowed with the wished acyclic orientation of ``self``.
    If `as_reorientations` is set at `True`, then acyclic orientations are given as a dictionary that maps edges to a boolean indicating if the edge is oriented in the same direction as in the given (di)graph.
    
    EXAMPLES::
            
        sage: from sage.graphs.acyclic_orientations import AcyclicOrientations
        sage: G = graphs.CompleteGraph(3)
        sage: list(AcyclicOrientations(G))
        [[(0, 1), (0, 2), (1, 2)], [(0, 1), (0, 2), (2, 1)], [(0, 1), (2, 0), (2, 1)], [(1, 0), (0, 2), (1, 2)], [(1, 0), (2, 0), (1, 2)], [(1, 0), (2, 0), (2, 1)]]

        sage: G = graphs.PetersenGraph()
        sage: A = AcyclicOrientations(G, as_reorientations=True)  #Put these 2 lines in a Jupyter cell and the 3 next in an other one, then run the second one several times.
        sage: O = next(A)  #Get the next acyclic orientation of G.
        sage: D = DiGraph([[(e[1], e[0]), e][O[e]] for e in O], pos = G.get_pos())  #Give to D the same embedding as G.
        sage: D.graphplot(edge_colors = {'blue':[e for e in O if O[e]], 'red':[(e[1], e[0]) for e in O if not O[e]]}).plot()
        Graphics object consisting of 26 graphics primitives
        
        sage: G = graphs.CompleteGraph(6)
        sage: sum(1 for _ in AcyclicOrientations(G))
        720

        sage: [G.is_clique() for G in graphs(5) if sum(1 for _ in AcyclicOrientations(G)) == factorial(5)]
        [True]
    
    The algorithm is a backtracking: it adds edges one by one checking if there is a cycle at each step. After finding a cycle or yielding an acyclic orientation, it backtracks.
    This could be improved either by implementing a block-and-cut graph decomposition (i.e. 2-connected component decomposition), or by maintaining a topological order.

    .. WARNING::
        
        Each acyclic orientation is a list of edges, if your graph has multiedges or isolated vertices, it will disapear.
        Loops are taken in account, thus a graph with a loop will have no acyclic orientation.
    
    .. SEEALSO::
    
        - :meth:`~sage.graphs.graph.orientations`
        - :meth:`~sage.graphs.graph.Graph.strong_orientation`
        - :meth:`~sage.graphs.orientations.strong_orientations_iterator`
        - :meth:`~sage.graphs.digraph_generators.DiGraphGenerators.nauty_directg`
        - :meth:`~sage.graphs.orientations.random_orientation`
    
    TESTS::
        
        sage: list(AcyclicOrientations(Graph()))
        []
        
        sage: list(AcyclicOrientations(Graph(5)))
        []
    """
    def __init__(self, G, as_reorientations=False):
        self.V, self.E = G.vertices(sort=False), (Graph(G, multiedges=False)).edges(sort=False, labels=False)
        self.l_V, self.l_E = len(self.V), len(self.E)
        self.has_loops = G.has_loops()
        self.D = DiGraph()
        self.D.add_vertices(self.V)
        self.mark, self.mem = 0, 0
        self.direction = 'f'
        self.A, self.Q, self.o = [], [], 1
        self.stop = False
        self.as_reorientations = as_reorientations
    
    def __iter__(self):
        return self
    
    def __next__(self):
        if not self.E or self.has_loops:
            raise StopIteration
        while not self.stop:
            if self.direction == 'f':
                if self.mark == self.l_E:
                    if self.as_reorientations:
                        output = {self.E[i] : bool(self.Q[i]+1) for i in range(self.l_E)}
                    else:
                        output = list(self.A)
                    if self.Q[-1] == 1:
                        self.Q.pop(-1)
                        self.D.delete_edge(self.A.pop(-1))
                        self.o = -1
                        self.mark -= 1
                    else:
                        self.direction = 'b'
                    return output
                else:
                    u, v = self.E[self.mark]
                    if self.o == -1:
                        u, v = v, u
                    
                    if self.is_flipable((u, v)):
                        self.D.add_edge((u, v))
                        self.A.append((u, v))
                        self.Q.append(self.o)
                        self.mark += 1
                        self.o = 1
                    
                    elif self.o == 1:
                        self.o = -1
                    else:
                        self.direction = 'b'
            elif self.direction == 'b':
                to_del = []
                while self.mark > 0 and self.Q[-1] == -1:
                    self.Q.pop(-1)
                    to_del.append(self.A.pop(-1))
                    self.mark -= 1
                if self.mark > 0:
                    self.Q.pop(-1)
                    self.D.delete_edges(to_del + [self.A.pop(-1)])
                    self.o = -1
                    self.mark -= 1
                    self.direction = 'f'
                else:
                    self.stop = True
        raise StopIteration
    
    def is_flipable(self, e):
        # Need a method 'has_directed_path(v, u)' to improve (i.e. look if there is a path v -> u): doesn't exists yet (maybe use dfs implemented methods?).
        return not self.D.shortest_path(e[1], e[0])