import warnings
from collections import deque
from sage.graphs.planar_maps.CustomSwap import CustomSwap
from sage.graphs.planar_maps.TopologicalDemiEdge import TopologicalDemiEdge
from sage.graphs.planar_maps.PermutationUtilsAbstractor import PermutationUtilsAbstractor
from sage.all import Permutation, Graph  # Import sage library
from sage.graphs.planar_maps.MapPermutation import *

try:
    import networkx as nx
except ImportError:
    nx = None

try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None


def transitiveCouplePermutation(sigma, alpha):
    """
    Check that sigma and alpha act transitively
    """
    assert alpha.size() == sigma.size()
    size = sigma.size()
    seen = [False] * (size + 1)
    seen[0] = seen[1] = True
    # Half-edges are numbered from 1 to size, included

    todo = [1]
    while todo:
        i = todo.pop()
        if not seen[alpha(i)]:
            todo.append(alpha(i))
            seen[alpha(i)] = True
        if not seen[sigma(i)]:
            todo.append(sigma(i))
            seen[sigma(i)] = True
    return False not in seen


class LabelledMap:
    """
    A class to represent a labelled map.

    Attributes
    ----------
    sigma : Permutation or MapPermutation
        Fixed-point free involution whose cycles are given by the edges.
    alpha : Permutation or MapPermutation
        Permutation that maps a half-edge to the half-edge incident to
        it in a clockwise direction around the vertex it belongs to.

    Methods
    -------
    """

    def __init__(
        self,
        sigma: Permutation | MapPermutation = None,
        alpha: Permutation | MapPermutation = None,
        adj=None,
        trust=False,
    ):
        r"""
        Initializes the labelled map from either the permutations alpha
        and sigma, or an adjacency list (giving each vertex a list of
        its neighbors in order; vertices must be numbered from 1 to n).

        INPUT:
        - ``sigma`` -- Permutation ; Fixed-point free involution whose
          cycles are given by the edges.
        - ``alpha`` -- Permutation ; Permutation that maps a half-edge
          to the half-edge incident to it in clockwise direction around
          the vertex it belongs to.
        - ``ajd``-- and adjacency list be careful the order of the
            node in your adjaceny will be used to choose the embedding
        - ``trust`` -- A parameter that indicates whether the validity check (i.e., whether the map is connex, etc.)
          should be skipped when initializing the map. It makes initialization faster but can be dangerous because
          if the map isn't well-formed, all the other methods become unsafe. You should be absolutely sure of your
          map's validity if you set this to true.

          - The advantage of setting `trust` to true is that it makes the initialization faster, which is useful when
          you are initializing a lot of big maps (like in long bijections).

          - Therefore, the best workflow is to leave it at the default during testing, and when you are 100% sure that
          your code works, set `trust = true` to gain a constant factor boost. By default, it is set to false.

        EXAMPLES::

            sage: sigma = Permutation([1, 3, 2, 5, 4, 6])
            sage: alpha = Permutation([(1, 2), (3, 4), (5, 6)])
            sage: LabelledMap(sigma, alpha)
            Labelled map | Sigma : [1, 3, 2, 5, 4, 6],
                           Alpha : [2, 1, 4, 3, 6, 5]

        TESTS::

            sage: sigma = Permutation([3, 4, 1, 2, 6, 5])
            sage: alpha = Permutation([(1, 2), (3, 4)])
            sage: map = LabelledMap(sigma, alpha)
            Traceback (most recent call last):
            ...
            ValueError: The two permutations do not have the same size

            sage: sigma = Permutation([3, 4, 1, 2, 5])
            sage: alpha = Permutation([(1, 2), (3, 4, 5)])
            sage: LabelledMap(sigma, alpha)
            Traceback (most recent call last):
            ...
            ValueError: The permutation alpha is not an involution

            sage: sigma = Permutation([3, 4, 1, 2, 5])
            sage: alpha = Permutation([2, 1, 3, 5, 4])
            sage: LabelledMap(sigma, alpha)
            Traceback (most recent call last):
            ...
            ValueError: The permutation alpha should not have any fixed points

            sage: sigma = Permutation([1, 2, 3, 4, 5, 6])
            sage: alpha = Permutation([(1, 2), (3, 4), (5, 6)])
            sage: LabelledMap(sigma, alpha)
            Traceback (most recent call last):
            ...
            ValueError: The graph is not connected

            sage: adj = [(3,), (1, 3), (2,)]
            sage: LabelledMap(adj=adj)
            Traceback (most recent call last):
            ...
            ValueError: Invalid adjacency list

            sage: sigma = Permutation([1, 3, 2, 5, 4, 6])
            sage: adj = [(3,), (1, 3), (2,)]
            sage: LabelledMap(sigma, adj=adj)
            Traceback (most recent call last):
            ...
            ValueError: Cannot build the map from both an adjacency list
            and permutations
        """

        # Set it to false during
        # Debugging
        self._production = True
        if not self._production:
            trust = False

        if sigma is None and alpha is None and adj is None:
            raise ValueError("sigma and alpha and adj cannot all be none")

        if adj is not None and (sigma is not None or alpha is not None):
            raise ValueError(
                """Cannot build the map from both an adjacency list
                 and permutations"""
            )

        self.topologicalMap = {}
        if adj is None:
            self._build_from_permutations(sigma, alpha, trust)
        else:
            self._build_from_adj(adj, trust)
        self.extend()

    def extend(self):
        self.sigmaUtilsAbstractor = PermutationUtilsAbstractor(self.sigma)
        self.phiUtilsAbstractor = PermutationUtilsAbstractor(self.phi)
        # Initialising the topologicalMap
        for i in range(1, self.q + 1):
            self.topologicalMap[i] = TopologicalDemiEdge(self, i)

    def _build_from_permutations(self, sigma, alpha, trust):
        r"""
        Initializes the labelled map from the underlying permutations.
        """
        self.alpha = alpha
        self.sigma = sigma

        if isinstance(self.sigma, Permutation):
            self.sigma = MapPermutation(self.sigma, trust=trust)

        if isinstance(self.alpha, Permutation):
            self.alpha = MapPermutation(self.alpha, trust=trust)

        self.phi = self.alpha.right_action_product(self.sigma)

        if trust:
            return
        size = self.sigma.size()

        if self.sigma.size() != self.alpha.size():
            raise ValueError("The two permutations do not have the same size")

        if (
            self.alpha.right_action_product(self.alpha)
            != MapPermutation(size)
        ):
            raise ValueError("The permutation alpha is not an involution")

        if self.alpha.number_of_fixed_points() != 0:
            raise ValueError(
                "The permutation alpha should not have any fixed points"
            )

        if not transitiveCouplePermutation(self.sigma, self.alpha):
            raise ValueError("The graph is not connected")

    def _build_from_adj(self, adj, trust):
        """
        Initializes the labelled map from an adjacency list.

        Raises:
            ValueError if the adjacency list is invalid.
        """
        n = len(adj)

        if sum(map(len, adj)) % 2 != 0:
            raise ValueError("Invalid adjacency list")

        pairs = []  # Pairs of half-edges representing the edges
        cycles = []  # Cycles of outgoing half-edges for each vertex

        edges = {}
        iEdge = 1

        for u in range(1, n + 1):
            c = []
            for v in adj[u - 1]:
                other = None
                if (v, u) in edges:
                    # Check before setting (u, v) for loops
                    other = edges[(v, u)]
                    edges.pop((v, u))

                if other:
                    pairs.append((iEdge, other))
                else:
                    edges[(u, v)] = iEdge

                c.append(iEdge)
                iEdge += 1

            cycles.append(tuple(c))

        if edges:
            raise ValueError("Invalid adjacency list")

        self._build_from_permutations(
            MapPermutation(cycles), MapPermutation(pairs), trust)

    def buildGraph(self):
        """
        Builds the multigraph corresponding to this labelled map.
        Vertices are numbered from 1 to n.

        Returns:
            A multigraph corresponding to the labelled map.
        """
        vertices = self.sigma.to_cycles()
        # Associates a vertex to each half-edge
        corres = [0] * (2 * self.m + 1)
        for i in range(1, len(vertices) + 1):
            for k in vertices[i - 1]:
                corres[k] = i

        edges = []
        for i in range(1, 2 * self.m + 1):
            if i < self.alpha(i):  # Avoid adding edges twice
                edges.append((corres[i], corres[self.alpha(i)]))

        return Graph(edges, loops=True, multiedges=True)

    def show(
        self,
        show_halfedges=True,
        show_vertices=False,
        ax=None,
        use_sage_viewer=False
    ):
        """
        Show the planar map using networkx (unless unavailable or
        use_sage_viewer is set to True, in which case the default
        sage viewer is used). Half-edges numbers are displayed
        closest to the node they depart from.

        Note: The order of the edges may not be displayed correctly
        if the genus is not 0. For half-edges to be shown, the
        networkx viewer must be used.

        INPUT:

        - ``show_halfedges`` -- bool; whether to show half-edges
        numbers on the plot (default: True).
        - ``show_vertices`` -- bool; whether to show vertex labels
        (default: False).
        - ``ax`` -- matplotlib.axes._axes.Axes; if specified, draw
        the graph into the given matplotlib axes (default: None).
        - ``use_sage_viewer`` -- bool; use the default sage viewer
        if True (default: False).
        """
        vertices = self.sigma.to_cycles()

        real_n_vertices = len(vertices)  # Remove multiedges and loops
        real_n_halfedges = self.m        # These half-edges should not be drawn

        alpha = self.alpha
        sigma = self.sigma
        m = self.m

        def minmax(i, j):
            """Ensure edges always go from lowest to highest vertex id."""
            return min(i, j), max(i, j)

        # Three dictionaries to store half-edge IDs
        edge_labels_head = {}  # (i, j): half-edge from i to j
        edge_labels_tail = {}  # (i, j): half-edge from j to i
        edge_labels_middle = {}  # Used for loops

        # Map half-edge i to its corresponding vertex
        corres = [0] * (2 * self.m + 1)
        for i in range(1, len(vertices) + 1):
            for k in vertices[i - 1]:
                corres[k] = i

        def break_down(i, write_labels):
            nonlocal alpha, sigma, corres, vertices, m
            if write_labels:
                # Avoid writing half-edge numbers during the first step
                edge_labels_middle[(corres[i], len(vertices) + 1)] = i
                edge_labels_middle[(
                    corres[alpha(i)], len(vertices) + 1)] = alpha(i)

            # Add a new vertex v, and break down the edge whose half-edges
            # are i & alpha(i) into 2 edges (i, 2*m+1) and (2*m+2, alpha(i)).
            alpha *= MapPermutation([(alpha(i), 2 *
                                    m + 1, i, 2 * m + 2)], trust=self._production)
            sigma *= CustomSwap([(2 * m + 1, 2 * m + 2)])

            corres.append(len(vertices) + 1)
            corres.append(len(vertices) + 1)
            vertices.append((2 * m + 1, 2 * m + 2))
            m += 1

        def break_loop(i):
            j = alpha(i)
            vertex = len(vertices)

            break_down(i, False)
            break_down(2 * m, False)

            edge_labels_middle[(corres[i], vertex + 1)] = i
            edge_labels_middle[(corres[j], vertex + 2)] = j

        # For each loop a-a, add a new vertex v and replace the edge a-a
        # with two edges a-v, v-a.
        for i in range(1, 2 * m + 1):
            if corres[i] == corres[alpha(i)]:
                break_loop(i)

        # Handle each vertex and break down edges if needed.
        for v in range(1, len(vertices) + 1):
            seen_vertices = set()
            for i in vertices[v - 1]:
                if corres[alpha(i)] in seen_vertices:
                    break_down(i, True)
                else:
                    seen_vertices.add(corres[alpha(i)])
                    if (
                        corres[i] <= real_n_vertices
                        and corres[alpha(i)] <= real_n_vertices
                    ):
                        if corres[i] < corres[alpha(i)]:
                            edge_labels_head[minmax(
                                corres[i], corres[alpha(i)])] = i
                        else:
                            edge_labels_tail[minmax(
                                corres[i], corres[alpha(i)])] = i

        # Build the graph embedding
        embedding = {
            i: list(corres[alpha(k)] for k in vertices[i - 1])[::-1]
            for i in range(1, len(vertices) + 1)
        }

        edges = [
            (corres[i], corres[alpha(i)]) for i in range(1, 2 * m + 1) if i < alpha(i)
        ]

        g = Graph(edges, loops=False, multiedges=False)
        g.set_embedding(embedding)

        if not use_sage_viewer and nx is None:
            warnings.warn(
                "Package networkx not found; falling back to default sage viewer. "
                "Consider installing networkx using sage --pip install networkx"
            )
            use_sage_viewer = True
        if not use_sage_viewer and plt is None:
            warnings.warn(
                "Package matplotlib not found; falling back to default sage viewer. "
                "Consider installing matplotlib using sage --pip install matplotlib"
            )
            use_sage_viewer = True

        vertex_size = min(300, 1000 / len(vertices) **.5)

        if use_sage_viewer:
            layout = "planar" if self.genus() == 0 else "spring"
            g.show(
                layout=layout,
                vertex_size=vertex_size,
                vertex_labels={
                    i: str(i)
                    if i <= real_n_vertices and show_vertices
                    else ""
                    for i in range(1, len(vertices) + 1)
                },
                vertex_colors={
                    "red": list(range(1, real_n_vertices + 1)),
                    "white": list(
                        range(real_n_vertices + 1, len(vertices) + 1)
                    ),
                },
                figsize=(8, 8),
            )
        else:
            layout = g.layout_planar() if self.genus() == 0 else g.layout()
            G = nx.DiGraph()
            G.add_edges_from(
                minmax(i, j) for (i, j, _) in g.edges()
            )

            nx.draw(
                G,
                layout,
                ax=ax,
                labels={
                    i: str(i)
                    if i <= real_n_vertices
                    else ""
                    for i in range(1, len(vertices) + 1)
                },
                node_size=[vertex_size] * real_n_vertices
                + [0] * (len(vertices) - real_n_vertices),
                nodelist=list(range(1, len(vertices) + 1)),
                arrows=False,
                with_labels=show_vertices,
                node_color="red",
            )

            if show_halfedges:
                nx.draw_networkx_edge_labels(
                    G,
                    layout,
                    ax=ax,
                    rotate=False,
                    edge_labels=edge_labels_head,
                    label_pos=0.7,
                )
                nx.draw_networkx_edge_labels(
                    G,
                    layout,
                    ax=ax,
                    rotate=False,
                    edge_labels=edge_labels_tail,
                    label_pos=0.3,
                )
                nx.draw_networkx_edge_labels(
                    G,
                    layout,
                    ax=ax,
                    rotate=False,
                    edge_labels=edge_labels_middle,
                    label_pos=0.5,
                )

            if ax is None:
                plt.show()

    def __repr__(self):
        r"""
        Return string representation of this labelled map

        OUTPUT: The string representation of this labelled map

        EXAMPLES::

            sage: sigma = Permutation([1,3,2,5,4,6])
            sage: alpha = Permutation([(1,2),(3,4),(5,6)])
            sage: str(LabelledMap(sigma, alpha))
            'Labelled map | Sigma : [1, 3, 2, 5, 4, 6], Alpha : [2, 1, 4, 3, 6, 5]'

        """
        return "Labelled map | Sigma : " + \
            str(self.sigma) + ", Alpha : " + str(self.alpha)

    def _numberOfFaces(self):
        """
        A method that return the number of faces of the  map
        -------
        Returns:
             The number of faces of self
        -------
        O(1)
        """
        return self.phiUtilsAbstractor.numberOfCycles()

    def numberOfFaces(self):
        """
        A method that return the number of faces of the labelled map

        OUTPUT: The number of faces of this labelled map

        EXAMPLES::

            sage: sigma = Permutation([1,3,2,5,4,6])
            sage: alpha = Permutation([(1,2),(3,4),(5,6)])
            sage: LabelledMap(sigma, alpha).numberOfFaces()
            1

            sage: sigma = Permutation([(1,6),(2,3),(4,5)])
            sage: alpha = Permutation([(1,2),(3,4),(5,6)])
            sage: LabelledMap(sigma, alpha).numberOfFaces()
            2

        TESTS::

            sage: sigma = Permutation([1,2])
            sage: alpha = Permutation([(1,2)])
            sage: LabelledMap(sigma, alpha).numberOfFaces()
            1

            sage: sigma = Permutation([1,3,2,4])
            sage: alpha = Permutation([(1,2),(3,4)])
            sage: LabelledMap(sigma, alpha).numberOfFaces()
            1

            sage: sigma = Permutation([2,1])
            sage: alpha = Permutation([(1,2)])
            sage: LabelledMap(sigma, alpha).numberOfFaces()
            2

            sage: sigma = Permutation([(1,3,5),(2,6,4)])
            sage: alpha = Permutation([(1,2),(3,4),(5,6)])
            sage: LabelledMap(sigma, alpha).numberOfFaces()
            3

            sage: sigma = Permutation([(1,7,3,5),(2,6,4,8)])
            sage: alpha = Permutation([(1,2),(3,4),(5,6),(7,8)])
            sage: LabelledMap(sigma, alpha).numberOfFaces()
            4

        .. NOTE::

            Complexity is O(1) where m is the number of edges
        """
        return self._numberOfFaces()

    def _numberOfNodes(self):
        """
        A method that returns the number of vertices of the map
        -------
        Returns:
             The number of nodes of self
        -------
        O(1)
        where m is the number of edges
        """

        return self.sigmaUtilsAbstractor.numberOfCycles()

    def numberOfNodes(self):
        """
        A method that returns the number of nodes
        or vertices of this labelled map

        OUTPUT: The number of nodes of this labelled map

        EXAMPLES::

            sage: sigma = Permutation([1,3,2,5,4,6])
            sage: alpha = Permutation([(1,2),(3,4),(5,6)])
            sage: LabelledMap(sigma, alpha).numberOfNodes()
            4

            sage: sigma = Permutation([(1,6),(2,3),(4,5)])
            sage: alpha = Permutation([(1,2),(3,4),(5,6)])
            sage: LabelledMap(sigma, alpha).numberOfNodes()
            3

        TESTS::

            sage: sigma = Permutation([1,2])
            sage: alpha = Permutation([(1,2)])
            sage: LabelledMap(sigma, alpha).numberOfNodes()
            2

            sage: sigma = Permutation([1,3,2,4])
            sage: alpha = Permutation([(1,2),(3,4)])
            sage: LabelledMap(sigma, alpha).numberOfNodes()
            3

            sage: sigma = Permutation([2,1])
            sage: alpha = Permutation([(1,2)])
            sage: LabelledMap(sigma, alpha).numberOfNodes()
            1

            sage: sigma = Permutation([(1,3,5),(2,6,4)])
            sage: alpha = Permutation([(1,2),(3,4),(5,6)])
            sage: LabelledMap(sigma, alpha).numberOfNodes()
            2

            sage: sigma = Permutation([(1,7,3,5),(2,6,4,8)])
            sage: alpha = Permutation([(1,2),(3,4),(5,6),(7,8)])
            sage: LabelledMap(sigma, alpha).numberOfNodes()
            2

        .. NOTE::

            Complexity is O(m) where m is the number of edges
        """
        return self._numberOfNodes()

    def numberOfEdges(self):
        """
        A method that returns the number of edges of this labelled map

        OUTPUT: The number of edges of this labelled map

        EXAMPLES::

            sage: sigma = Permutation([1,3,2,5,4,6])
            sage: alpha = Permutation([(1,2),(3,4),(5,6)])
            sage: LabelledMap(sigma, alpha).numberOfEdges()
            3

            sage: sigma = Permutation([(1,6),(2,3),(4,5)])
            sage: alpha = Permutation([(1,2),(3,4),(5,6)])
            sage: LabelledMap(sigma, alpha).numberOfEdges()
            3

        TESTS::

            sage: sigma = Permutation([1,2])
            sage: alpha = Permutation([(1,2)])
            sage: LabelledMap(sigma, alpha).numberOfEdges()
            1

            sage: sigma = Permutation([1,3,2,4])
            sage: alpha = Permutation([(1,2),(3,4)])
            sage: LabelledMap(sigma, alpha).numberOfEdges()
            2

            sage: sigma = Permutation([2,1])
            sage: alpha = Permutation([(1,2)])
            sage: LabelledMap(sigma, alpha).numberOfEdges()
            1

            sage: sigma = Permutation([(1,3,5),(2,6,4)])
            sage: alpha = Permutation([(1,2),(3,4),(5,6)])
            sage: LabelledMap(sigma, alpha).numberOfEdges()
            3

            sage: sigma = Permutation([(1,7,3,5),(2,6,4,8)])
            sage: alpha = Permutation([(1,2),(3,4),(5,6),(7,8)])
            sage: LabelledMap(sigma, alpha).numberOfEdges()
            4

        .. NOTE::

            Complexity is O(1)
        """
        return self.sigma.size() // 2

    def genus(self):
        """
        Returns the genus of this labelled map.

        The genus is the minimum number of handles that must be added
        to a sphere to embed the map without edge crossings.

        OUTPUT:

        int: The genus of this labelled map.

        EXAMPLES::

            sage: sigma = Permutation([(2, 3, 1, 5), (4,), (6,)])
            sage: alpha = Permutation([(1, 2), (3, 4), (5, 6)])
            sage: LabelledMap(sigma, alpha).genus()
            0
            sage: adj = [(5, 4, 2), (1, 3, 6), (4, 7, 2), (8, 3, 1),\
            (8, 1, 6), (5, 2, 7), (3, 8, 6), (7, 4, 5)]
            sage: LabelledMap(adj=adj).genus()
            0

            sage: adj = [(4, 5, 6), (4, 5, 6), (4, 5, 6), (1, 2, 3),\
            (1, 2, 3), (1, 2, 3)]
            sage: LabelledMap(adj=adj).genus()
            1

        .. NOTE::

            Complexity is O(1)
        """
        return (
            self.numberOfEdges()
            + 2
            - self.numberOfFaces()
            - self.numberOfNodes()
        ) // 2

    def force_planar(self):
        """
        Returns a map of genus 0 with the same underlying graph
        if it is planar.

        For example, the adjacency list [(2, 3, 4), (3, 4, 1),
        (4, 1, 2), (1, 2, 3)] is a valid adjacency list for the
        complete graph with 4 nodes, but its genus is 1.
        This method could, for instance, build a map from the list
        [(2, 4, 3), (3, 4, 1), (4, 2, 1), (1, 2, 3)],
        which is the same graph  but with genus 0.

        Raises an error if the underlying graph is not planar.

        Returns
        -------
        LabelledMap
            Another map corresponding to the above description.
        """
        g = self.buildGraph()

        if not g.is_planar(set_embedding=True):
            raise ValueError(
                """The force_planar method can only be used on maps
                 whose underlying graph is planar."""
            )

        e = g.get_embedding()
        adj = [tuple(reversed(e[i])) for i in range(1, len(e) + 1)]

        return LabelledMap(adj=adj, trust=self._production)

    def getSpanningTree(self):
        """
        A method that returns any spanning tree of the  map
        -------
        Returns:
             A spanning tree of self
        -------
        O(m)
        """

        g = self.buildGraph()
        n = g.order()

        tree = Graph()

        seen = [False] * (n + 1)
        seen[0] = True

        def dfs(u):
            seen[u] = True
            for v in g.neighbor_iterator(u):
                if not seen[v]:
                    tree.add_edge(u, v)
                    dfs(v)

        dfs(1)

        assert False not in seen

        return tree

    def dual(self):
        """
        A method that return the dual of this map

        OUTPUT:

        The dual map of this map.

        EXAMPLES::

            sage: sigma = Permutation([1,3,2,5,4,6])
            sage: alpha = Permutation([(1,2),(3,4),(5,6)])
            sage: Map = LabelledMap(sigma, alpha)
            sage: dualMap = Map.dual()
            sage: dualMap.buildGraph().edges(labels=False)
            [(1, 1), (1, 1), (1, 1)]

            sage: sigma = Permutation([(1,6),(2,3),(4,5)])
            sage: alpha = Permutation([(1,2),(3,4),(5,6)])
            sage: Map = LabelledMap(sigma, alpha)
            sage: dualMap = Map.dual()
            sage: dualMap.buildGraph().edges(labels=False)
            [(1, 2), (1, 2), (1, 2)]

        .. NOTE::

            Complexity is O(m) where m is the number of edges
        """
        return LabelledMap(self.phi.inverse(), self.alpha, trust=self._production)

    def diameter(self):
        """
        A method that return the diameter of this map,
        i.e. the maximum length of a simple path in the map.

        OUTPUT:

        The diameter of this map

        EXAMPLES::

            sage: sigma = Permutation([1,3,2,5,4,6])
            sage: alpha = Permutation([(1,2),(3,4),(5,6)])
            sage: LabelledMap(sigma, alpha).diameter()
            3

            sage: sigma = Permutation([(1,6),(2,3),(4,5)])
            sage: alpha = Permutation([(1,2),(3,4),(5,6)])
            sage: LabelledMap(sigma, alpha).diameter()
            1

        .. NOTE::

            Complexity is O(m*n) where m is the number of edges
            and n is the number of nodes.
        """
        graph = self.buildGraph()
        return Graph.diameter(graph)

    def derivedMap(self):
        """
        A method that return the derived map of this map

        OUTPUT:

        The canonical representant of the derived map of this map

        EXAMPLES::

            sage: sigma = Permutation([1,3,2,5,4,6])
            sage: alpha = Permutation([(1,2),(3,4),(5,6)])
            sage: derivedMap = LabelledMap(sigma, alpha).derivedMap()
            sage: derivedMap.numberOfFaces()
            6
            sage: derivedMap.numberOfEdges()
            12
            sage: derivedMap.numberOfNodes()
            8
            sage: derivedMap.buildGraph().edges(labels=False)
            [(1, 2), (2, 3), (2, 3), (2, 4), (3, 5), (3, 5),
            (3, 6), (3, 6), (4, 5), (5, 7), (6, 7), (6, 8)]

            sage: sigma = Permutation([(1,6),(2,3),(4,5)])
            sage: alpha = Permutation([(1,2),(3,4),(5,6)])
            sage: derivedMap = LabelledMap(sigma, alpha).derivedMap()
            sage: derivedMap.numberOfFaces()
            6
            sage: derivedMap.numberOfEdges()
            12
            sage: derivedMap.numberOfNodes()
            8
            sage: derivedMap.buildGraph().edges(labels=False)
            [(1, 2), (1, 3), (2, 4), (2, 5), (2, 6), (3, 4),
            (3, 5), (3, 7), (4, 8), (5, 8), (6, 8), (7, 8)]

        TESTS::

            sage: sigma = Permutation([1,2])
            sage: alpha = Permutation([(1,2)])
            sage: derivedMap = LabelledMap(sigma, alpha).derivedMap()
            sage: derivedMap.numberOfFaces()
            2
            sage: derivedMap.numberOfEdges()
            4
            sage: derivedMap.numberOfNodes()
            4
            sage: derivedMap.buildGraph().edges(labels=False)
            [(1, 2), (2, 3), (2, 3), (2, 4)]

            sage: alpha = Permutation([(1,2),(3,4)])
            sage: sigma = Permutation([1,2])
            sage: alpha = Permutation([(1,2)])
            sage: derivedMap = LabelledMap(sigma, alpha).derivedMap()
            sage: derivedMap.numberOfFaces()
            2
            sage: derivedMap.numberOfEdges()
            4
            sage: derivedMap.numberOfNodes()
            4
            sage: derivedMap.buildGraph().edges(labels=False)
            [(1, 2), (2, 3), (2, 3), (2, 4)]

            sage: sigma = Permutation([2,1])
            sage: alpha = Permutation([(1,2)])
            sage: derivedMap = LabelledMap(sigma, alpha).derivedMap()
            sage: derivedMap.numberOfFaces()
            2
            sage: derivedMap.numberOfEdges()
            4
            sage: derivedMap.numberOfNodes()
            4
            sage: derivedMap.buildGraph().edges(labels=False)
            [(1, 2), (1, 2), (2, 3), (2, 4)]

            sage: sigma = Permutation([(1,3,5),(2,6,4)])
            sage: alpha = Permutation([(1,2),(3,4),(5,6)])
            sage: derivedMap = LabelledMap(sigma, alpha).derivedMap()
            sage: derivedMap.numberOfFaces()
            6
            sage: derivedMap.numberOfEdges()
            12
            sage: derivedMap.numberOfNodes()
            8
            sage: derivedMap.buildGraph().edges(labels=False)
            [(1, 2), (1, 3), (1, 4), (2, 5), (2, 6), (2, 7), (3, 6),
             (3, 7), (3, 8), (4, 5), (4, 7), (4, 8)]

            sage: sigma = Permutation( [(1,7,3,5),(2,6,4,8)])
            sage: alpha = Permutation( [(1,2),(3,4),(5,6),(7,8)])
            sage: derivedMap = LabelledMap(sigma, alpha).derivedMap()
            sage: derivedMap.numberOfFaces()
            8
            sage: derivedMap.numberOfEdges()
            16
            sage: derivedMap.numberOfNodes()
            10
            sage: derivedMap.buildGraph().edges(labels=False)
            [(4, 10), (6, 10), (1, 2), (1, 3), (1, 4), (1, 6), (2, 5),
            (2, 7), (2, 8), (3, 7), (3, 8), (3, 9), (4, 8), (4, 9),
            (5, 6), (6, 8)]

        .. NOTE::

            Complexity is O(m) where m is the number of edges
        """
        K = 8 * self.m + 1

        derivedAlphaList = list(range(1, K))
        derivedSigmaList = list(range(1, K))

        invPhi = self.phi.inverse()

        m = int(self.m)

        for i in list(range(1, K)):
            if i <= 2 * m:
                derivedAlphaList[i - 1] = i + 2 * m
                derivedSigmaList[i - 1] = self.sigma(i)
            elif i > 2 * m and i <= 4 * m:
                derivedAlphaList[i - 1] = i - 2 * m
                derivedSigmaList[i - 1] = i + 4 * m
            elif i > 4 * m and i <= 6 * m:
                derivedAlphaList[i - 1] = i + 2 * m
                derivedSigmaList[i - 1] = invPhi(i - 4 * m) + 4 * m
            else:
                derivedAlphaList[i - 1] = i - 2 * m
                derivedSigmaList[i - 1] = self.alpha(i - 6 * m) + 2 * m

        derivedSigma = MapPermutation(derivedSigmaList, trust=self._production)
        derivedAlpha = MapPermutation(derivedAlphaList, trust=self._production)
        return LabelledMap(derivedSigma, derivedAlpha, trust=self._production).canonicalRepresentant()

    def quadrangulation(self):
        """
        There is a bijection between rooted maps with m edges of genus
        g and bipartite quadrangulations with m faces and genus g.
        This function returns the canonical representant of the rooted
        bipartite quadrangulation associated to self.

        Returns:
            The canonical representant of the bipartite rooted
            quadrangulation associated to rooted(self).

        Complexity: O(m), where m is the number of edges.
        """
        return self.incidenceMap()

    def inverseQuadrangulation(self):
        """
        This function is the inverse of quadrangulation, assuming self
        is a bipartite quadrangulation. It returns a map M such that
        M.quadrangulation() = self.canonicalRepresentant() and
        M = M.canonicalRepresentant(). If self isn't a bipartite
        quadrangulation, it raises an error.

        Returns:
            The canonical representant of the inverse of rooted(self)
            by quadrangulation, if self is a bipartite quadrangulation.
            Otherwise, it raises an error.

        Complexity: O(m), where m is the number of edges.
        """
        bipartition = self.getBipartition()

        if bipartition is None or not self.isQuandrangulation():
            raise ValueError("Self isn't a bipartite quadrangulation")

        alpha, sigma, phi = self.alpha, self.sigma, self.phi
        alphaInvList = [-1 for _ in range(self.m)]
        sigmaInvList = [-1 for _ in range(self.m)]
        colorFace = bipartition[1]

        corres = [-1 for _ in range(self.m + 1)]
        invCorres = [-1 for _ in range(2 * self.m + 1)]
        corres[1], invCorres[1] = 1, 1

        cnt = 2
        for i in range(2, 2 * self.m + 1):
            if bipartition[i] == colorFace:
                corres[cnt], invCorres[i] = i, cnt
                cnt += 1

        for invDemiEdge in range(1, self.m + 1):
            alphaInvList[invDemiEdge - 1] = invCorres[
                sigma(alpha(sigma(alpha(corres[invDemiEdge]))))
            ]
            sigmaInvList[invDemiEdge - 1] = invCorres[
                alpha(sigma(alpha(corres[invDemiEdge])))
            ]

        alphaInv = MapPermutation(alphaInvList, trust=self._production)
        sigmaInv = MapPermutation(sigmaInvList, trust=self._production)

        return LabelledMap(
            sigma=sigmaInv,
            alpha=alphaInv, trust=self._production).canonicalRepresentant()

    def incidenceMap(self):
        """
        A method that returns the incidence map of this map
        as its canonical representant.

        OUTPUT:

        Incidence map of self.

        EXAMPLES::

            sage: sigma = Permutation([1, 3, 2, 5, 4, 6])
            sage: alpha = Permutation([(1, 2), (3, 4), (5, 6)])
            sage: incidenceMap = LabelledMap(sigma, alpha).incidenceMap()
            sage: incidenceMap.numberOfNodes()
            5
            sage: incidenceMap.numberOfFaces()
            3
            sage: incidenceMap.numberOfEdges()
            6
            sage: incidenceMap.buildGraph().edges(labels=False)
            [(1, 2), (1, 3), (1, 3), (1, 4), (1, 4), (1, 5)]

        .. NOTE::

            Complexity is O(m), where m is the number of edges.
        """

        invPhi = self.phi.inverse()
        invPhiCycles = invPhi.to_cycles()

        quadDemiEdge = 1
        corres = [-1]
        invCorres = list(range(2 * self.m + 1))
        sigmaQuadList = []

        for cycle in invPhiCycles:
            startQuadDemiEdge = quadDemiEdge
            for demiEdge in cycle:
                if quadDemiEdge != startQuadDemiEdge:
                    sigmaQuadList.append(quadDemiEdge)
                corres.append(demiEdge)
                invCorres[demiEdge] = quadDemiEdge
                quadDemiEdge += 1
            sigmaQuadList.append(startQuadDemiEdge)

        numberOfQuadEdge = quadDemiEdge - 1
        alphaQuadList = list(range(2 * numberOfQuadEdge))

        for quadDemiEdge in range(1, numberOfQuadEdge + 1):
            demiEdge = corres[quadDemiEdge]
            turnedDemiEdge = self.sigma(demiEdge)
            quadDemiEdgePrime = invCorres[turnedDemiEdge]

            sigmaQuadList.append(quadDemiEdgePrime + numberOfQuadEdge)
            alphaQuadList[quadDemiEdge - 1] = quadDemiEdge + numberOfQuadEdge
            alphaQuadList[quadDemiEdge + numberOfQuadEdge - 1] = quadDemiEdge

        alphaQuad = MapPermutation(alphaQuadList, trust=self._production)
        sigmaQuad = MapPermutation(sigmaQuadList, trust=self._production)
        relabelList = [i + 1 for i in range(2 * numberOfQuadEdge)]

        for quadDemiEdge in range(1, numberOfQuadEdge + 1):
            relabelList[quadDemiEdge - 1] = corres[quadDemiEdge]
            relabelList[quadDemiEdge + numberOfQuadEdge - 1] = (
                relabelList[quadDemiEdge - 1] + numberOfQuadEdge
            )

        relabelPerm = MapPermutation(relabelList, trust=self._production)
        return LabelledMap(sigmaQuad, alphaQuad, trust=self._production).relabel(
            relabelPerm
        ).canonicalRepresentant()

    def getRootedMapCorrespondance(
            self, otherMap, rootDemiEdge, return_map_perm=False, trust=False):
        """
        A method that returns a labelling of the demi-edges of this map
        giving `otherMap` while keeping `rootDemiEdge` invariant.
        It checks if this map and `otherMap` represent the same
        rooted map at `rootDemiEdge`.

        INPUT:

        - ``otherMap`` -- LabelledMap; the other map
        - ``rootDemiEdge`` -- int; the edge on which to root
        - ``return_map_perm`` -- ; whether or not to return a MapPermutation default to False
        - ``trust`` -- ; whether or not to trust that there is a correspondence default to
            False
        OUTPUT:

        Returns `t`, a permutation mapping the demi-edges of `self`
        to those of `otherMap`, or None if they don't represent
        the same rooted map at `rootDemiEdge`.

        EXAMPLES::

            sage: sigma = Permutation([1, 3, 2, 5, 4, 6])
            sage: alpha = Permutation([(1, 2), (3, 4), (5, 6)])
            sage: tau = Permutation([(1, 3)])
            sage: Map = LabelledMap(sigma, alpha)
            sage: relabelMap = Map.relabel(tau)
            sage: Map.getRootedMapcorrespondence(relabelMap, 2)
            [3, 2, 1, 4, 5, 6]

        .. NOTE::

            Complexity is O(m), where m is the number of edges.
        """
        if not trust and otherMap.numberOfEdges() != self.numberOfEdges():
            return None

        m = self.numberOfEdges()
        tList = [-1 for _ in range(2 * m)]
        seen = [False for _ in range(2 * m)]

        alpha = self.alpha
        sigma = self.sigma
        sigmaOther = otherMap.sigma
        alphaOther = otherMap.alpha

        tList[rootDemiEdge - 1] = rootDemiEdge
        p = [rootDemiEdge]
        seen[rootDemiEdge - 1] = True

        while p:
            u = p.pop()
            if not seen[alpha(u) - 1]:
                seen[alpha(u) - 1] = True
                tList[alpha(u) - 1] = alphaOther(tList[u - 1])
                p.append(alpha(u))
            if not seen[sigma(u) - 1]:
                seen[sigma(u) - 1] = True
                tList[sigma(u) - 1] = sigmaOther(tList[u - 1])
                p.append(sigma(u))

        try:
            if not return_map_perm:
                t = Permutation(tList)
            else:
                t = MapPermutation(tList, trust=trust)
        except ValueError:
            return None

        if not trust and self.relabel(t) != otherMap:
            return None

        return t

    def relabel(self, tau):
        """
        A method that returns this map with demi-edge `i`
        relabeled by `tau(i)`.

        INPUT:

        - ``tau`` -- Permutation or Permutation; a permutation on the demi-edges
          representing the relabelling.

        OUTPUT:

        The relabeled map.

        EXAMPLES::

            sage: sigma = Permutation([1, 3, 2, 5, 4, 6])
            sage: alpha = Permutation([(1, 2), (3, 4), (5, 6)])
            sage: tau = Permutation([(1, 3)])
            sage: LabelledMap(sigma, alpha).relabel(tau)
            Labelled map | Sigma : [2, 1, 3, 5, 4, 6],
            Alpha : [4, 3, 2, 1, 6, 5]

        .. NOTE::

            Complexity is O(m), where m is the number of edges.
        """
        if isinstance(tau, Permutation):
            tau = MapPermutation(tau, trust=self._production)
        invTau = tau.inverse()
        relabeledSigma = tau.left_action_product(
            invTau.right_action_product(self.sigma)
        )
        relabeledAlpha = tau.left_action_product(
            invTau.right_action_product(self.alpha)
        )
        return LabelledMap(relabeledSigma, relabeledAlpha, trust=self._production)

    def __eq__(self, other):
        """Checks equality between two LabelledMap instances."""
        if isinstance(other, self.__class__):
            return self.sigma == other.sigma and self.alpha == other.alpha
        return False

    def tetravalance(self):
        """
        Returns the canonical representative of the rooted,
        face-bicolorable tetravalent map associated with the current map.

        There is a bijection between rooted maps with m edges of
        genus g and face-bicolorable tetravalent rooted maps of
        genus g with m vertices.

        OUTPUT:

        The canonical representative of a tetravalent
        bicolorable rooted map.

        .. NOTE::

            Complexity is O(m), where m is the number of edges.
        """
        return self.edgeMap()

    def edgeMap(self):
        """
        A method that return the edge map of this map
        as its canonical representant

        OUTPUT:

        A canonical representant of the edge map of self

        EXAMPLES::

            sage: sigma = Permutation([1,3,2,5,4,6])
            sage: alpha = Permutation([(1,2),(3,4),(5,6)])
            sage: edgeMap = LabelledMap(sigma, alpha).edgeMap()
            sage: edgeMap.numberOfNodes()
            3
            sage: edgeMap.numberOfFaces()
            5
            sage: edgeMap.numberOfEdges()
            6
            sage: edgeMap.buildGraph().edges(labels=False)
            [(1, 1), (1, 2), (1, 2), (2, 3), (2, 3), (3, 3)]

            sage: sigma = Permutation( [(1,6),(2,3),(4,5)])
            sage: alpha = Permutation( [(1,2),(3,4),(5,6)])
            sage: edgeMap = LabelledMap(sigma, alpha).edgeMap()
            sage: edgeMap.numberOfNodes()
            3
            sage: edgeMap.numberOfFaces()
            5
            sage: edgeMap.numberOfEdges()
            6
            sage: edgeMap.buildGraph().edges(labels=False)
            [(1, 2), (1, 2), (1, 3), (1, 3), (2, 3), (2, 3)]

        TESTS::

            sage: sigma = Permutation([1,2])
            sage: alpha = Permutation([(1,2)])
            sage: edgeMap = LabelledMap(sigma, alpha).edgeMap()
            sage: edgeMap.numberOfNodes()
            1
            sage: edgeMap.numberOfFaces()
            3
            sage: edgeMap.numberOfEdges()
            2
            sage: edgeMap.buildGraph().edges(labels=False)
            [(1, 1), (1, 1)]

            sage: sigma = Permutation([1,3,2,4])
            sage: alpha = Permutation([(1,2),(3,4)])
            sage: edgeMap = LabelledMap(sigma, alpha).edgeMap()
            sage: edgeMap.numberOfNodes()
            2
            sage: edgeMap.numberOfFaces()
            4
            sage: edgeMap.numberOfEdges()
            4
            sage: edgeMap.buildGraph().edges(labels=False)
            [(1, 1), (1, 2), (1, 2), (2, 2)]

            sage: sigma = Permutation([2,1])
            sage: alpha = Permutation([(1,2)])
            sage: edgeMap = LabelledMap(sigma, alpha).edgeMap()
            sage: edgeMap.numberOfNodes()
            1
            sage: edgeMap.numberOfFaces()
            3
            sage: edgeMap.numberOfEdges()
            2
            sage: edgeMap.buildGraph().edges(labels=False)
            [(1, 1), (1, 1)]

            sage: sigma = Permutation( [(1,3,5),(2,6,4)])
            sage: alpha = Permutation( [(1,2),(3,4),(5,6)])
            sage: edgeMap = LabelledMap(sigma, alpha).edgeMap()
            sage: edgeMap.numberOfNodes()
            3
            sage: edgeMap.numberOfFaces()
            5
            sage: edgeMap.numberOfEdges()
            6
            sage: edgeMap.buildGraph().edges(labels=False)
            [(1, 2), (1, 2), (1, 3), (1, 3), (2, 3), (2, 3)]

            sage: sigma = Permutation( [(1,7,3,5),(2,6,4,8)])
            sage: alpha = Permutation( [(1,2),(3,4),(5,6),(7,8)])
            sage: edgeMap = LabelledMap(sigma, alpha).edgeMap()
            sage: edgeMap.numberOfNodes()
            4
            sage: edgeMap.numberOfFaces()
            6
            sage: edgeMap.numberOfEdges()
            8
            sage: edgeMap.buildGraph().edges(labels=False)
            [(1, 2), (1, 2), (1, 3), (1, 3), (2, 4), (2, 4),
            (3, 4), (3, 4)]

        .. NOTE::

            Complexity is O(m) where m is the number of edges
        """
        invSigma = self.sigma.inverse()
        alpha = self.alpha
        sigma = self.sigma
        m = self.m

        # The number of edge in the edge map
        L = int(2 * m)

        alphaListEdgeMap = [-1 for k in range(2 * L)]
        sigmaListEdgeMap = [-1 for k in range(2 * L)]

        # Construction of alpha and sigma for the edge map
        for k in range(1, L + 1):
            alphaListEdgeMap[k - 1] = k + L
            alphaListEdgeMap[k + L - 1] = k

            t = invSigma(k)
            sigmaListEdgeMap[k - 1] = L + t

            j = sigma(k)

            sigmaListEdgeMap[k + L - 1] = alpha(j)

        alphaEdgeMap = MapPermutation(alphaListEdgeMap, trust=self._production)
        sigmaEdgeMap = MapPermutation(sigmaListEdgeMap, trust=self._production)

        return LabelledMap(sigmaEdgeMap, alphaEdgeMap, trust=self._production, ).canonicalRepresentant()

    def isQuandrangulation(self):
        """
        A function to test whether or not self is a quandrangulation
        ---
        Returns:
            A boolean indicating if self is a quandrangulation or not
        """
        phi_cycles = self.phi.to_cycles()

        for i in range(len(phi_cycles)):
            if len(phi_cycles[i]) != 4:
                return False

        return True

    def isBipartite(self):
        """
        Returns : A boolean indicating whether or not self is bipartite
        -------
        O(m)
        where m is the number of edges
        """
        return self.getBipartition() is not None

    def getBipartition(self):
        """
        If self isn't bipartite this method will return None.
        Otherwise, it will return a tab clr such that clr[i](=0,1)
        for a demi-edge i gives the color of the node on which it is
        attached. A node (i.e., a cycle of sigma) will be white
        (resp. black) if all its elements are of color 0 (resp. 1).
        clr[0] = -1 because 0 isn't a valid demi-edge.
        -------
        Returns:
            clr: None if self isn't bipartite; otherwise, it gives
            the color description mentioned above.
        -------

        O(m)
        where m is the number of edges
        """
        clr = [-1 for i in range(2 * self.m + 1)]
        clr[1] = 0
        alpha = self.alpha
        sigma = self.sigma
        phi = self.phi
        p = []
        p.append(1)

        seen = [False for i in range(2 * self.m)]

        seen[0] = True
        cnt = 2
        while len(p) > 0:
            u = p.pop()
            if not seen[sigma(u) - 1]:
                seen[sigma(u) - 1] = True
                p.append(sigma(u))
                clr[sigma(u)] = clr[u]

            if not seen[alpha(u) - 1]:
                seen[alpha(u) - 1] = True
                p.append(alpha(u))
                clr[alpha(u)] = (1 + clr[u]) % 2
        sigma_cycles = sigma.to_cycles()
        for i in range(len(sigma_cycles)):
            r = clr[sigma_cycles[i][0]]
            for j in range(len(sigma_cycles[i])):
                if r != clr[sigma_cycles[i][j]]:
                    return None
        for i in range(1, 2 * self.m + 1):
            if clr[i] == clr[alpha(i)]:
                return None

        return clr

    def canonicalRepresentant(self):
        """
        This function returns the canonical representant of rooted(self),
        i.e a labelled  map such that M and self are representant of
        the same rooted map and M is the canonical representant.
        -------
        Returns:
             The canonical representant of rooted(self)
        -------
        O(m)
        where m is the number of edges
        """
        relabelList = [-1 for i in range(2 * self.m)]
        alphaCycles = self.alpha.to_cycles()

        sigma = self.sigma
        alpha = self.alpha
        rootDemiEdge = 1
        relabelList[rootDemiEdge - 1] = rootDemiEdge

        p = deque()

        p.append(rootDemiEdge)

        seen = [False for i in range(2 * self.m)]

        seen[rootDemiEdge - 1] = True
        cnt = 2
        while len(p) > 0:
            u = p.popleft()
            if not seen[sigma(u) - 1]:
                seen[sigma(u) - 1] = True
                p.append(sigma(u))
                relabelList[sigma(u) - 1] = cnt
                cnt += 1

            if not seen[alpha(u) - 1]:
                seen[alpha(u) - 1] = True
                p.append(alpha(u))
                relabelList[alpha(u) - 1] = cnt
                cnt += 1

        relabel = MapPermutation(relabelList, trust=self._production)

        return self.relabel(relabel)

    def isPlaneTree(self):
        """
        A method that returns a boolean indicating
        if self is a plane Tree or not
        -------
        Returns:
             A boolean indicating if self is plane tree or not
        -------
        O(m)
        where m is the number of edges
        """
        return self.numberOfFaces() == 1 \
            and self.numberOfEdges() == self.numberOfNodes() - 1

    def schaefferTree(self, markedDemiEdge):
        """
        The Schaeffer surjection from rooted bipartite quadrangulation of genus g with k face and a marked node to
        rooted one face map (tree in the case g=0) of genus g with k edges and a labelling of its nodes (i.e a function on the nodes of the tree considered up to translation
        such that if u and v are adjacent f(u) and f(v) differs by at most one) such that for every rooted one face map T only two rooted marked bipartite quadrangulation give T.
        Given a markDemiEdge which is the corresponding marked node(a node is just a cycle of self.sigma) , this method will return the canonical representant of the
        rooted one face map associated to rooted(self) and a labelling on its demi edge such that f(node) is the common value of all its demi edge(note that labelling[0] is present but it deosn't have any meaning).
        If self isn't a bipartite quandrangulation this function will raise an error.
        -------
        Args:
            -markedDemiEdge a demi edge on the node which is marked
        Returns:
            - tree: The canonical representant of the rooted one face map corresponding to the above description
            - labelling: A list of labelling on the demi edge of tree corresponding to the above description
        -------
        O(m)
        where m is the number of edges
        """
        if not self.isBipartite() or not self.isQuandrangulation():
            raise ValueError("Self isn't a bipartite quadrangulation")

        sigma = self.sigma
        alpha = self.alpha
        phi = self.phi
        labellingQuad = [-1 for i in range(2 * self.m + 1)]

        p = deque()

        nodes = sigma.to_cycles()

        demiEdgeToNodeId = [-1 for i in range(2 * self.m + 1)]

        for i in range(len(nodes)):
            for j in range(len(nodes[i])):
                demiEdgeToNodeId[nodes[i][j]] = i

        startNodeId = demiEdgeToNodeId[markedDemiEdge]

        p.append(startNodeId)

        distNodes = [-1 for i in range(len(nodes))]

        distNodes[startNodeId] = 0

        seen = [False for i in range(len(nodes))]

        seen[startNodeId] = True

        while len(p) > 0:
            nodeId = p.popleft()

            for demiEdge in nodes[nodeId]:
                labellingQuad[demiEdge] = distNodes[nodeId]

                alphaDemiEdge = alpha(demiEdge)
                alphaNodeId = demiEdgeToNodeId[alphaDemiEdge]
                if not seen[alphaNodeId]:
                    distNodes[alphaNodeId] = 1 + distNodes[nodeId]
                    seen[alphaNodeId] = True
                    p.append(alphaNodeId)

        phi_cycles = phi.to_cycles()

        corres = [-1 for i in range(2 * len(phi_cycles) + 1)]
        invCorres = [-1 for i in range(2 * self.m + 1)]

        cnt = 1
        for i in range(len(phi_cycles)):
            A, D, C, B = phi_cycles[i][0], phi_cycles[i][1], phi_cycles[i][2], phi_cycles[i][3]
            link = None
            if labellingQuad[A] == labellingQuad[C]:
                if labellingQuad[B] == labellingQuad[D]:
                    if labellingQuad[A] > labellingQuad[B]:
                        link = (A, C)
                    else:
                        link = (D, B)
                else:
                    if labellingQuad[B] > labellingQuad[D]:
                        link = (C, B)
                    else:
                        link = (D, A)
            else:
                if labellingQuad[A] > labellingQuad[C]:
                    link = (A, B)
                else:
                    link = (C, D)

            corres[cnt] = link[0]
            invCorres[link[0]] = cnt
            cnt += 1
            corres[cnt] = link[1]
            invCorres[link[1]] = cnt
            cnt += 1

        numberOfTreeDemiEdge = cnt - 1

        alphaTreeList = [-1 for i in range(numberOfTreeDemiEdge)]
        sigmaTreeList = [-1 for i in range(numberOfTreeDemiEdge)]

        prelabelling = [-1 for i in range(numberOfTreeDemiEdge + 1)]

        for treeDemiEdge in range(1, numberOfTreeDemiEdge + 1):
            if treeDemiEdge % 2 == 0:
                alphaTreeList[treeDemiEdge - 1] = treeDemiEdge - 1
            else:
                alphaTreeList[treeDemiEdge - 1] = treeDemiEdge + 1

            U = corres[treeDemiEdge]
            prelabelling[treeDemiEdge] = labellingQuad[U]

            turnU = sigma(U)
            while invCorres[turnU] == -1:
                turnU = sigma(turnU)

            sigmaTreeList[treeDemiEdge - 1] = invCorres[turnU]

        A, D, C, B = 1, phi(1), phi(phi(1)), phi(phi(phi(1)))

        treeRoot = None

        if invCorres[A] != -1 and invCorres[B] != -1:
            if labellingQuad[A] > labellingQuad[B]:
                treeRoot = invCorres[A]
            else:
                treeRoot = invCorres[B]

        elif invCorres[D] != -1 and invCorres[A] != -1:
            if labellingQuad[D] > labellingQuad[A]:
                treeRoot = invCorres[D]
            else:
                treeRoot = invCorres[A]

        elif invCorres[B] != -1 and invCorres[C] != -1:
            if labellingQuad[B] > labellingQuad[C]:
                treeRoot = invCorres[C]
            else:
                treeRoot = invCorres[B]

        elif invCorres[C] != -1 and invCorres[D] != -1:
            if labellingQuad[C] > labellingQuad[D]:
                treeRoot = invCorres[D]
            else:
                treeRoot = invCorres[C]

        elif invCorres[A] != -1 and invCorres[C] != -1:
            treeRoot = invCorres[A]

        else:
            treeRoot = invCorres[B]

        tau = CustomSwap([(1, treeRoot)])
        alphaTree = MapPermutation(alphaTreeList, trust=self._production)
        sigmaTree = MapPermutation(sigmaTreeList, trust=self._production)

        tree = LabelledMap(alpha=alphaTree, sigma=sigmaTree,
                           trust=self._production, ).relabel(tau)

        canonicalTree = tree.canonicalRepresentant()

        tauCanonical = tree.getRootedMapCorrespondance(
            canonicalTree, rootDemiEdge=1, return_map_perm=True, trust=self._production)

        labelling = [-1 for i in range(numberOfTreeDemiEdge + 1)]

        for i in range(1, numberOfTreeDemiEdge + 1):
            labelling[tauCanonical(tau(i))] = prelabelling[i]

        return canonicalTree, labelling

    def inverseShaefferTree(self, labelled, returnMarkedDemiEdge=True):
        """
        This method is the inverse of the schaefferTree method given that self is a one face map it will return a quadruple
        (quadA,quadB,markedDemiEdgeA,markedDemiEdgeB) where quadA and quaB are in canonical form and we have
        the following( we will note nodeA and nodeB the node on which markedDemiEdgeA(resp markedDemiEdgeB) is attached in A(resp B))
        (rooted(quadA),nodeA) are (rooted(quadB),nodeB) are the only marked rooted quandrangulation such that calling schaefferTree with quadA (resp quadB) with
        any demi edge attached to nodeA (resp nodeB) give rooted(self) in particular quadA.schaefferTree(markedDemiEdgeA) = self.canonicalForm() same for B.
        Note that if returnMarkedDemiEdge = False it will only return (quadA,quadB)
        -------
        Args:
            -labelled a list of size 2*m+1 such that for the demiEdge i labelled[i] is the labelled of its attached node,
             0 isn't a valid demiEdge so labelled[0] can take any value it will be ignored.
            -returnMarkedDemiEdge : a parameter indicating whether or not to return the markedDemiEdge default to true
        Returns:
            -(quadA,quadB,markedDemiEdgeA,markedDemiEdgeB) as in the above description if returnMarkedDemiEdge = True otherwise (quadA,quadB) corresponding to the above description
            ,if self isn't a one face map it will raise an error
        -------
        O(m)
        where m is the number of edges
        -------
        TESTS::

            sage: sigma = Permutation( [(1,6),(2,3),(4,5)])
            sage: alpha = Permutation( [(1,2),(3,4),(5,6)])
            sage: tri = LabelledMap(sigma,alpha)
            sage: bigQuad = tri.derivedMap().derivedMap(
            ).derivedMap().derivedMap().quadrangulation()
            sage: bigQuad.numberOfEdges()
            1536
            sage: markedDemiEdge = 750
            sage: sct,labelled = bigQuad.schaefferTree(markedDemiEdge = markedDemiEdge)
            sage: quadA,quadB,markedDemiEdgeA,markedDemiEdgeB = sct.inverseShaefferTree(labelled)
            sage: quadA.schaefferTree(markedDemiEdge = markedDemiEdgeA)[0] == sct.canonicalRepresentant(
            ) and quadB.schaefferTree(markedDemiEdge = markedDemiEdgeB)[0] == sct.canonicalRepresentant()
            True
        """
        alpha = self.alpha
        sigma = self.sigma

        phi = self.phi

        nextAction = alpha.right_action_product(sigma.inverse())
        backAction = nextAction.inverse()
        # In case the labelling isn't such that the minimum labelling is 1
        # We update the labelling
        # maxLabel = labelled[1]
        # minLabel = labelled[1]
        # for i in range(1, len(labelled)):
        #     # if labelled[i] < minLabel:
        #     #     minLabel = labelled[i]
        #     # if labelled[i] > maxLabel:
        #     #     maxLabel = labelled[i]

        #     minLabel = min(minLabel, labelled[i])
        #     maxLabel = max(maxLabel, labelled[i])

        maxLabel = max(labelled[1:])
        minLabel = min(labelled[1:])

        for i in range(1, len(labelled)):
            labelled[i] -= minLabel - 1
        maxLabel -= minLabel - 1
        minLabel = 1
        nodes = sigma.to_cycles()

        # We get a correspondence between nodes and demiEdge
        nodesId = [-1 for i in range(2 * self.m + 1)]
        for i in range(len(nodes)):
            for j in range(len(nodes[i])):
                nodesId[nodes[i][j]] = i + 1

        alphaQuadList = [-1 for i in range(4 * self.m)]
        p = [[] for i in range(maxLabel + 1)]
        sigmaQuadCycleDemiEdge = [deque() for i in range(2 * self.m + 1)]
        sigmaQuadCycle = [[]]

        check = [-1 for i in range(2 * self.m + 1)]

        oneOrdered = []
        backPartner = [-1 for i in range(2 * self.m + 1)]
        corres = [-1 for i in range(2 * self.m + 1)]

        # We're turning around the tree in a  clockwise manner
        # When we connect each demi edge A the first demiEdge before it
        # Which label is -1 less than the label of A
        # If we see that the demi edge is of label one we just add it
        # to a list of demi edge of label 1 and continue
        # And similarly during this loop we compute the back partner of
        # each demi edge A of not maximal label which is define as the first
        # demi edge B such that B has label on more than the label of A
        # And such that B is first one to connect to A
        # We also construct the cycle around each demiEdge
        root = 1
        cnt = 1
        curDemiEdge = root
        N = 2 * self.m
        while N != 0:
            curLabel = labelled[curDemiEdge]
            curNode = -1
            otherNodeId = -1
            if curLabel > 1 and len(
                    p[curLabel - 1]) > 0 and check[curDemiEdge] == -1:
                check[curDemiEdge] = -2

                otherDemiEdge = p[curLabel - 1][-1]

                N -= 1

                if backPartner[otherDemiEdge] == -1:
                    backPartner[otherDemiEdge] = curDemiEdge

                corres[curDemiEdge] = cnt

                sigmaQuadCycleDemiEdge[curDemiEdge].appendleft(cnt)
                sigmaQuadCycleDemiEdge[otherDemiEdge].append(cnt + 1)

                alphaQuadList[cnt - 1] = cnt + 1
                alphaQuadList[cnt] = cnt

                cnt += 2

            if labelled[curDemiEdge] == 1 and check[curDemiEdge] == -1:
                check[curDemiEdge] == -2

                corres[curDemiEdge] = cnt

                otherNodeId = 0

                sigmaQuadCycleDemiEdge[curDemiEdge].appendleft(cnt)

                sigmaQuadCycle[otherNodeId].append(cnt + 1)

                alphaQuadList[cnt - 1] = cnt + 1
                alphaQuadList[cnt] = cnt
                cnt += 2
                N -= 1
            p[curLabel].append(curDemiEdge)
            curDemiEdge = nextAction(curDemiEdge)

        # For demi edge  A such that sigmaQuadCycleDemiEdge[A]>1 we need
        # We need to make a shift on sigmaQuadCycleDemiEdge[A] such that the first demi edge in the cycle is
        # The one corresponding to his back partner
        # This is needed to have a correct merge later on when we merge them by
        # node
        for curDemiEdge in range(1, 2 * self.m + 1):

            sigmaQuadCycleDemiEdge[curDemiEdge] = list(
                sigmaQuadCycleDemiEdge[curDemiEdge])
            if len(sigmaQuadCycleDemiEdge[curDemiEdge]) == 1:
                continue
            curLabel = labelled[curDemiEdge]

            firstDemiEdgeQuad = sigmaQuadCycleDemiEdge[curDemiEdge][0]
            restDemiEdgeQuad = sigmaQuadCycleDemiEdge[curDemiEdge][1:]

            startId = -1

            otherDemiEdge = backPartner[curDemiEdge]

            P = len(restDemiEdgeQuad)

            for i in range(P):
                if alphaQuadList[restDemiEdgeQuad[i] -
                                 1] == corres[otherDemiEdge]:
                    startId = i

            newDemiEdgeQuadCycle = []
            newDemiEdgeQuadCycle.append(firstDemiEdgeQuad)

            for i in range(P):
                newDemiEdgeQuadCycle.append(
                    restDemiEdgeQuad[(i + startId) % P])

            sigmaQuadCycleDemiEdge[curDemiEdge] = newDemiEdgeQuadCycle

        sigmaQuadCycle[0] = tuple(sigmaQuadCycle[0])
        # Merging cycle on demi edge per  node
        # And also making transforming them into tuple
        for node in nodes:
            accum = []
            for demiEdge in node:
                for e in sigmaQuadCycleDemiEdge[demiEdge][1:]:
                    accum.append(e)
                accum.append(sigmaQuadCycleDemiEdge[demiEdge][0])

            sigmaQuadCycle.append(tuple(accum))

        sigmaQuad = MapPermutation(sigmaQuadCycle, trust=self._production)
        alphaQuad = MapPermutation(alphaQuadList, trust=self._production)
        quad = LabelledMap(sigma=sigmaQuad, alpha=alphaQuad,
                           trust=self._production, )
        numberOfQuadDemiEdge = len(alphaQuadList)

        phiQuad = quad.phi
        alphaQuad = quad.alpha
        U = corres[root]
        X, W, V = phiQuad(U), phiQuad(phiQuad(U)), phiQuad(phiQuad(phiQuad(U)))

        # There is two quadragulation possible depending on which root we choose
        # Here we calculate the two possible permutation corresponding to the
        # Two possible root
        tauA = None
        tauB = None
        if labelled[root] == labelled[alpha(root)]:
            tauA = CustomSwap([(U, root)])
            tauB = CustomSwap([(X, root)])

        elif labelled[root] > labelled[alpha(root)]:
            tauA = CustomSwap([(U, root)])
            tauB = CustomSwap([(V, root)])

        else:
            U = corres[alpha(root)]
            X, W, V = phiQuad(U), phiQuad(
                phiQuad(U)), phiQuad(phiQuad(phiQuad(U)))
            tauA = CustomSwap([(W, root)])
            tauB = CustomSwap([(X, root)])

        quadA = quad.relabel(tauA)
        quadB = quad.relabel(tauB)

        if returnMarkedDemiEdge is True:
            quadACanonical = quadA.canonicalRepresentant()
            quadBCanonical = quadB.canonicalRepresentant()

            canonicalTauA = quadA.getRootedMapCorrespondance(
                otherMap=quadACanonical, rootDemiEdge=root, return_map_perm=True, trust=self._production)
            canonicalTauB = quadB.getRootedMapCorrespondance(
                otherMap=quadBCanonical, rootDemiEdge=root, return_map_perm=True, trust=self._production)

            markedDemiEdge = sigmaQuadCycle[0][0]

            markedDemiEdgeA = tauA(markedDemiEdge)
            markedDemiEdgeB = tauB(markedDemiEdge)

            markedDemiEdgeA = canonicalTauA(markedDemiEdgeA)
            markedDemiEdgeB = canonicalTauB(markedDemiEdgeB)

            return quadACanonical, quadBCanonical, markedDemiEdgeA, markedDemiEdgeB

        return quadA.canonicalRepresentant(), quadB.canonicalRepresentant()

    def nodes(self):
        """
        This function return the nodes of self as cycle of self.sigma
        ----
        Returns :
            A list of cycle of sigma representing nodes of self
        ----
        O(m)
        """
        return self.sigma.to_cycles()

    def faces(self):
        """
        This function return the faces of self as cycle of self.phi
        ----
        Returns :
            A list of cycle of phi representing faces of self
        ----
        O(m)
        """
        return self.phi.to_cycles()

    def getDyckPath(self, isCanonical=False):
        """
        There is a canonical bijection between rooted planar trees
        with m edges and dyck path of size m, this method return the
        associated dyck path to rooted(self) if self is a tree
        if self isn't a tree it will raise an error.
        -----
        Args:
            isCanonical: A boolean indicating if self
            is already in canonical form
        Returns :
            A list of size 2*m representing the dyck path associated to
            rooted(self) +1 for step up and -1 for step down if self is
            a plane tree otherwise it will raise an error
        -----
        O(m)
        """
        if not self.isPlaneTree():
            raise ValueError("Self isn't a plane tree.")

        canonicalTree = self.canonicalRepresentant()
        canonicalAlpha = canonicalTree.alpha
        canonicalPhi = canonicalTree.phi

        seen = [False for i in range(2 * self.m + 1)]

        curDemiEdge = 1

        dyckPath = []

        while not seen[curDemiEdge]:
            dyckPath.append(-1 if seen[canonicalAlpha(curDemiEdge)] else 1)
            seen[curDemiEdge] = True
            curDemiEdge = canonicalPhi(curDemiEdge)
        return dyckPath

    @property
    def q(self):
        """
        Returns: The number of demi edge of self
        ----
        O(1)
        """
        return 2 * self.m

    def getTopologicalDemiEdge(self, demiEdge):
        """
        Args:
            demiEdge an index
        Returns:
            The TopologicalDemiEdge associated to demiEdge
        ----
        O(1)
        """
        return self.topologicalMap[demiEdge]

    def getListTopologicalDemiEdge(self):
        """
        Returns:
            The list of TopologicalDemiEdge in self such that the ith element
            is the TopologicalDemiEdge associated to the i+1 index
        ----
        O(m)
        """
        lst = []
        for i in range(1, self.q + 1):
            lst.append(self.getTopologicalDemiEdge(i))
        return lst

    def X(self, demiEdge):
        """
        Args:
            demiEdge an index
        Returns:
            The TopologicalDemiEdge associated to demiEdge
        ----
        O(1)
        """
        return self.getTopologicalDemiEdge(demiEdge)

    def XList(self):
        """
        Returns:
            The list of TopologicalDemiEdge in self such that the ith element
            is the TopologicalDemiEdge associated to the i+1 index
        ----
        O(m)
        """
        return self.getListTopologicalDemiEdge()

    @property
    def m(self):
        """
        Returns: The number of edge of self
        ----
        O(1)
        """
        return self.numberOfEdges()

    @property
    def f(self):
        """
        Returns: The number of face of self
        ----
        O(1)
        """
        return self.numberOfFaces()

    @property
    def n(self):
        """
        Returns: The number of node of self
        ----
        O(1)
        """
        return self.numberOfNodes()

    def pretty_print(self):
        """
        A prettier print
        """
        print(f"""
            Alpha: {self.alpha.to_cycles()}
            Sigma (Node): {self.sigma.to_cycles()}
            Phi (Face): {self.phi.to_cycles()}
        """)

    def copy(self):
        """
        Returns : A copy of self
        -------
        O(1)
        """
        return LabelledMap(self.sigma, self.alpha, trust=self._production, )

    def areOnTheSameNode(self, demiEdgeA, demiEdgeB):
        """
        Args:
            demiEdgeA,demiEdgeB: Index representing demi edges of self
        Returns:
            A boolean indicating whether or note demiEdgeA and demiEdgeB are on the node
        -------
        O(1)
        """
        return self.sigmaUtilsAbstractor.sameCycle(demiEdgeA, demiEdgeB)

    def areOnTheSameFace(self, demiEdgeA, demiEdgeB):
        """
        Args:
            demiEdgeA,demiEdgeB: Index representing demi edges of self
        Returns:
            A boolean indicating whether or note demiEdgeA and demiEdgeB are on the same face
        -------
        O(1)
        """
        return self.phiUtilsAbstractor.sameCycle(demiEdgeA, demiEdgeB)

    def demiEdgesOnTheSameNode(self, demiEdge):
        """
        Args:
            demiEdge: An index representing a demi edge of self
        Returns:
            A list of demiEdge on the same node as demiEdge
        -------
        O(d) where d is the number of demi edge on  the node
        """
        lst = []
        lst.append(demiEdge)
        curDemiEdge = self.sigma(demiEdge)
        while curDemiEdge != demiEdge:
            lst.append(curDemiEdge)
            curDemiEdge = self.sigma(curDemiEdge)
        return lst

    def demiEdgesOnTheSameFace(self, demiEdge):
        """
        Args:
            demiEdge: An index representing a demi edge of self
        Returns:
            A list of demiEdge on the same face as demiEdge
        -------
        O(f) where f is the number of demi edge on the face
        """
        lst = []
        lst.append(demiEdge)
        curDemiEdge = self.phi(demiEdge)
        while curDemiEdge != demiEdge:
            lst.append(curDemiEdge)
            curDemiEdge = self.phi(curDemiEdge)
        return lst

    def numberInTheSameFace(self, demiEdge):
        """
        Args:
            demiEdge
        Returns:
            The number of  demi edge on the same face as demi edge
        ----
        O(1)
        """
        return self.phiUtilsAbstractor.numberInCycle(demiEdge)

    def numberInTheSameNode(self, demiEdge):
        """
        Args:
            demiEdge
        Returns:
            The number of  demi edge on the same node as demi edge
        ----
        O(1)
        """

        return self.sigmaUtilsAbstractor.numberInCycle(demiEdge)

    def checkTwoInTheSameFace(self, listDemiEdges):
        """
        A method that will return a boolean indicating whether or not
        two demiEdge are on the same face
        -----
        Args:
            listDemiEdges: A list of demi edges
        OUTPUT:
            a boolean indicating whether or not there is two demi edge on the
            same face
        -----
        O(len(listDemiEdges))
        """
        return self.phiUtilsAbstractor.checkTwoInTheSameCycle(listDemiEdges)

    def checkTwoInTheSameNode(self, listDemiEdges):
        """
        A method that will return a boolean indicating whether or not
        two demiEdge are on the same node
        -----
        Args:
            listDemiEdges: A list of demi edges
        OUTPUT:
            a boolean indicating whether or not there is two demi edge on the
            same node
        -----
        O(len(listDemiEdges))
        """

        return self.sigmaUtilsAbstractor.checkTwoInTheSameCycle(listDemiEdges)

    @property
    def g(self):
        """
        The genus of self
        ----
        O(1)
        """
        return self.genus()

    def isTriangulation(self):
        """
        Returns: A boolean indicating if self is a triangulation or not
        ----
        O(m)
        """

        faces = self.faces()
        for face in faces:
            if len(face) != 3:
                return False
        return True
