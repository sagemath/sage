from sage.graphs.graph import Graph


class GraphHomomorphismCounter:
    def __init__(self, graph, target_graph, density_threshold=0.25, graph_clr=None, target_clr=None, colourful=False):
        r"""
        INPUT:

        - ``graph`` -- a Sage graph

        - ``target_graph`` -- the graph to which ``graph`` is sent

        - ``density_threshold`` -- float (default: ``0.25``); the density threshold for ``target_graph`` representation

        - ``graph_clr`` -- list (default: ``None``); a list of integers representing the colours of the vertices of ``graph``

        - ``target_clr`` -- list (default: ``None``); a list of integers representing the colors of the vertices of ``target_graph``

        - ``colourful`` -- boolean (default: ``False``); whether the graph homomorphism is colour-preserving
        """
        # Check for valid graphs and settings
        if not isinstance(graph, Graph):
            raise ValueError("first argument must be a sage Graph")
        if not isinstance(target_graph, Graph):
            raise ValueError("second argument must be a sage Graph")

        if colourful and (graph_clr is None or target_clr is None):
            raise ValueError("Both graph_clr and target_clr must be provided when colourful is True")

        graph._scream_if_not_simple()
        target_graph._scream_if_not_simple()

        # Only used in `__init__`
        from sage.graphs.graph_decompositions.tree_decomposition import make_nice_tree_decomposition, label_nice_tree_decomposition

        self.graph = graph
        self.target_graph = target_graph
        self.density_threshold = density_threshold
        self.graph_clr = graph_clr
        self.target_clr = target_clr
        self.colourful = colourful

        # Bookkeeping for colourful mappings
        self.actual_target_graph = target_graph
        self.actual_target_size = self.actual_target_graph.order()

        self.tree_decomp = graph.treewidth(certificate=True)
        self.nice_tree_decomp = make_nice_tree_decomposition(graph, self.tree_decomp)
        self.root = min(self.nice_tree_decomp)
        # self.root = sorted(self.nice_tree_decomp)[0]

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

        # `DP_table` is a list of lists of natural numbers (0 included)
        # Each element (list) corresponds to the (induced) homomorphisms of
        # a tree node.
        self.DP_table = [{} for _ in range(self.dir_labelled_TD.order())]

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
        for node in reversed(list(self.dir_labelled_TD.breadth_first_search(min(self.dir_labelled_TD)))):
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

        child_node_index, child_node_vtx = next(self.dir_labelled_TD.neighbor_out_iterator(node))
        child_node_vtx_tuple = tuple(child_node_vtx)

        mappings_length = self.actual_target_size ** len(node_vtx_tuple)

        mappings_count = [0] * mappings_length

        # Use the adjacency matrix when dense, otherwise the graph itself
        target_density = self.actual_target_graph.density()

        if target_density >= self.density_threshold:
            target = self.actual_target_graph.adjacency_matrix(
                vertices=self.actual_target_graph.vertices()
            )
        else:
            target = self.actual_target_graph

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

        child_node_index, child_node_vtx = next(self.dir_labelled_TD.neighbor_out_iterator(node))
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


### General helper functions

def get_node_index(node):
    r"""
    Return the index in the host graph of `node` in a (nice) tree decomposition
    """
    return node[0]

def get_node_content(node):
    r"""
    Return the content of `node` in a (nice) tree decomposition
    """
    return node[1]

def node_changes(labelled_TD):
    r"""
    Record introduced and forgotten nodes in a directed labelled nice tree
    decomposition.

    INPUT:

    - ``labelled_TD`` -- a directed labelled nice tree decomposition,
      with the root as source of the dirgraph

    OUTPUT:

    - A dictionary of recorded nodes, where the `key` is node index, and
      the `value` is the introduced/forgotten node

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.tree_decomposition import label_nice_tree_decomposition
        sage: bip_one_four = graphs.CompleteBipartiteGraph(1, 4)
        sage: nice_tree_decomp = bip_one_four.treewidth(certificate=True, nice=True)
        sage: root = sorted(nice_tree_decomp)[0]
        sage: dir_labelled_TD = label_nice_tree_decomposition(nice_tree_decomp, root, directed=True)
        sage: node_changes(dir_labelled_TD)
        {0: 0, 1: 1, 2: 1, 3: 4, 5: 4, 6: 4, 7: 3, 8: 2, 9: 0, 10: 0, 11: 3, 12: 2}
    """
    node_changes_dict = {}

    for node in sorted(labelled_TD):
        node_index, node_vertex_set = node

        node_type = labelled_TD.get_vertex(node)
        match node_type:
            case 'intro' | 'forget':
                child_vertex_set = labelled_TD.neighbors_out(node)[0][1]
                # Get one element from the one-element set
                (extra_vertex,) = node_vertex_set.symmetric_difference(child_vertex_set)
                node_changes_dict[node_index] = extra_vertex

    return node_changes_dict

def is_valid_mapping(mapped_vtx, mapped_nbhrs, target_graph):
    r"""
    Check if the mapping is valid.
    """
    from sage.graphs.graph import Graph

    if isinstance(target_graph, Graph):
        return all(target_graph.has_edge(mapped_vtx, vtx) for vtx in mapped_nbhrs)
    else:
        # Assume that `target_graph` is an adjacency matrix
        return all(target_graph[mapped_vtx, vtx] for vtx in mapped_nbhrs)


### For integer representation of mappings

def extract_bag_vertex(mapping, index, graph_size):
    r"""
    Extract the bag vertex at `index` from `mapping`
    """
    # Equivalent to taking the floor
    return mapping // (graph_size ** index) % graph_size

def add_vertex_into_mapping(new_vertex, mapping, index, graph_size):
    r"""
    Insert `new_vertex` at `index` into `mapping`
    """
    temp = graph_size ** index
    right_digits = mapping % temp
    left_digits = mapping - right_digits

    return graph_size * left_digits + temp * new_vertex + right_digits

def remove_vertex_from_mapping(mapping, index, graph_size):
    r"""
    Return a new mapping from removing vertex of `index` from `mapping`
    """
    left_digits = mapping - (mapping % (graph_size ** (index + 1)))
    right_digits = mapping % (graph_size ** index)

    return left_digits // graph_size + right_digits
