from sage.graphs.graph import Graph

from sage.graphs.homomorphisms.helper_functions import *

from sage.graphs.graph_decompositions.tree_decomposition import make_nice_tree_decomposition, label_nice_tree_decomposition

# In integer rep, the DP table is of the following form:
# { node_index: [1, 2, 3, 4, 5],
#   second_node_index: [10, 20, 30, 40, 50], ...}

class GraphHomomorphismCounter:
    def __init__(self, graph, target_graph, density_threshold=0.25, graph_clr=None, target_clr=None, colourful=False):
        r"""
        INPUT:

        - ``graph`` -- a Sage graph

        - ``target_graph`` -- the graph to which ``graph`` is sent

        - ``density_threshold`` (default: 0.25) -- the desnity threshold for `target_graph` representation

        - ``graph_clr`` (default: None) -- a list of integers representing the colours of the vertices of `graph`

        - ``target_clr`` (default: None) -- a list of integers representing the colours of the vertices of `target_graph`

        - ``colourful`` (default: False) -- whether the graph homomorphism is colour-preserving
        """
        self.graph = graph
        self.target_graph = target_graph
        self.density_threshold = density_threshold
        self.graph_clr = graph_clr
        self.target_clr = target_clr
        self.colourful = colourful

        # Bookkeeping for colourful mappings
        self.actual_target_graph = target_graph
        self.actual_target_size = len(self.actual_target_graph)

        if not isinstance(graph, Graph):
            raise ValueError("first argument must be a sage Graph")
        if not isinstance(target_graph, Graph):
            raise ValueError("second argument must be a sage Graph")

        if colourful and (graph_clr is None or target_clr is None):
            raise ValueError("Both graph_clr and target_clr must be provided when colourful is True")

        self.graph._scream_if_not_simple()
        self.target_graph._scream_if_not_simple()

        self.tree_decomp = graph.treewidth(certificate=True)
        self.nice_tree_decomp = make_nice_tree_decomposition(graph, self.tree_decomp)
        self.root = sorted(self.nice_tree_decomp)[0]

        # Make it into directed graph for better access
        # to children and parent, if needed
        #
        # Each node in a labelled nice tree decomposition
        # has the following form:
        #
        # (node_index, bag_vertices) node_type
        #
        # Example: (5, {0, 4}) intro
        self.dir_labelled_TD = label_nice_tree_decomposition(self.nice_tree_decomp, self.root, directed=True)

        # `node_changes_dict` is responsible for recording introduced and
        # forgotten vertices in a nice tree decomposition
        self.node_changes_dict = node_changes(self.dir_labelled_TD)

        # `DP_table` is a vector/list of dictionaries
        # Each element (dict) corresponds to the (induced) hom's of a tree node.
        # For each pair inside a dict, the key is the hom, and the value is the number of hom's:
        #
        # An example of K2 to K3, the values are arbitrary for demo purpose:
        #
        # [{((4, 0), (5, 1)): 10,
        #   ((4, 0), (5, 2)): 20,
        #   ((4, 1), (5, 0)): 30,
        #   ((4, 1), (5, 2)): 40,
        #   ((4, 2), (5, 0)): 50,
        #   ((4, 2), (5, 1)): 60}, {}, ...]
        self.DP_table = [{} for _ in range(len(self.dir_labelled_TD))]


    def count_homomorphisms(self):
        r"""
        Return the number of homomorphisms from the graph `G` to the graph `H`.

        A homomorphism from a graph `G` to a graph `H` is a function
        `\varphi : V(G) \mapsto V(H)`, such that for any edge `uv \in E(G)` the
        pair `\varphi(u) \varphi(v)` is an edge of `H`.

        For more information, see the :wikipedia:`Graph_homomorphism`.

        ALGORITHM:

        This is an implementation based on the proof of Prop. 1.6 in [CDM2017]_.

        OUTPUT:

        - an integer, the number of homomorphisms from `graph` to `target_graph`

        EXAMPLES::

            sage: from sage.graphs.homomorphisms.count_homomorphisms import GraphHomomorphismCounter
            sage: square = graphs.CycleGraph(4)
            sage: bip = graphs.CompleteBipartiteGraph(2, 4)
            sage: counter = GraphHomomorphismCounter(square, bip)
            sage: counter.count_homomorphisms()
            128
        """
        # Whether it's BFS or DFS, every node below join node(s) would be
        # computed already, so we can go bottom-up safely.
        for node in reversed(self.dir_labelled_TD.vertices()):
            node_type = self.dir_labelled_TD.get_vertex(node)

            match node_type:
                case 'intro':
                    self._add_intro_node(node)
                case 'forget':
                    self._add_forget_node(node)
                case 'join':
                    self._add_join_node(node)

                case _: 
                    self._add_leaf_node(node)

        return self.DP_table[0][0]

    ### Main adding functions

    def _add_leaf_node(self, node):
        r"""
        Add the computation result of processing leaf node(s) to the DP table.
        """
        node_index = get_node_index(node)
        self.DP_table[node_index] = [1]

    def _add_intro_node(self, node):
        r"""
        Add the computation result of processing introduce node(s) to the DP table.
        """
        # Basic setup
        node_index, node_vertices = node
        node_vtx_tuple = tuple(node_vertices)

        child_node_index, child_node_vtx = self.dir_labelled_TD.neighbors_out(node)[0]
        child_node_vtx_tuple = tuple(child_node_vtx)

        mappings_length = self.actual_target_size ** len(node_vtx_tuple)

        mappings_count = [0] * mappings_length

        # Use the adjacency matrix when dense, otherwise the graph itself
        target_density = self.actual_target_graph.density()
        target = self.actual_target_graph.adjacency_matrix() if target_density >= self.density_threshold else self.actual_target_graph

        # Intro node specifically
        intro_vertex = self.node_changes_dict[node_index]
        intro_vtx_index = node_vtx_tuple.index(intro_vertex) # Index of the intro vertex in the node/bag

        if self.colourful:
            intro_vtx_clr = self.graph_clr[intro_vertex]

        # Neighborhood of intro vertex in the bag
        intro_vtx_nbhs = [child_node_vtx_tuple.index(vtx) for vtx in child_node_vtx_tuple if self.graph.has_edge(intro_vertex, vtx)]

        child_DP_entry = self.DP_table[child_node_index]

        for mapped in range(len(child_DP_entry)):
            # Neighborhood of the mapped vertices of intro vertex in the target graph
            mapped_intro_nbhs = [extract_bag_vertex(mapped, vtx, self.actual_target_size) for vtx in intro_vtx_nbhs]

            mapping = add_vertex_into_mapping(0, mapped, intro_vtx_index, self.actual_target_size)

            for target_vtx in self.actual_target_graph:
                # If the colours do not match, skip current iteration and
                # move on to the next vertex
                if self.colourful:
                    target_vtx_clr = self.target_clr[target_vtx]
                    if intro_vtx_clr != target_vtx_clr:
                        mapping += self.actual_target_size ** intro_vtx_index
                        continue

                if is_valid_mapping(target_vtx, mapped_intro_nbhs, target):
                    mappings_count[mapping] = child_DP_entry[mapped]

                mapping += self.actual_target_size ** intro_vtx_index

        self.DP_table[node_index] = mappings_count

    def _add_forget_node(self, node):
        r"""
        Add the computation result of processing forget node(s) to the DP table.
        """
        # Basic setup
        node_index, node_vertices = node
        node_vtx_tuple = tuple(node_vertices)

        child_node_index, child_node_vtx = self.dir_labelled_TD.neighbors_out(node)[0]
        child_node_vtx_tuple = tuple(child_node_vtx)

        target_graph_size = len(self.target_graph)
        mappings_length = self.actual_target_size ** len(node_vtx_tuple)
        mappings_count = [0] * mappings_length

        # Forget node specifically
        forgotten_vtx = self.node_changes_dict[node_index]
        forgotten_vtx_index = child_node_vtx_tuple.index(forgotten_vtx)

        child_DP_entry = self.DP_table[child_node_index]

        for mapping in range(mappings_length):
            sum = 0
            extended_mapping = add_vertex_into_mapping(0, mapping, forgotten_vtx_index, self.actual_target_size)

            for _ in range(target_graph_size):
                sum += child_DP_entry[extended_mapping]
                extended_mapping += self.actual_target_size ** forgotten_vtx_index

            mappings_count[mapping] = sum

        self.DP_table[node_index] = mappings_count

    def _add_join_node(self, node):
        r"""
        Add the computation result of processing join node(s) to the DP table.
        """
        node_index, node_vertices = node
        left_child, right_child  = [vtx for vtx in self.dir_labelled_TD.neighbors_out(node)
                                        if get_node_content(vtx) == node_vertices]
        left_child_index = get_node_index(left_child)
        right_child_index = get_node_index(right_child)

        mappings_count = [left_count * right_count for left_count, right_count
                            in zip(self.DP_table[left_child_index], self.DP_table[right_child_index])]

        self.DP_table[node_index] = mappings_count
