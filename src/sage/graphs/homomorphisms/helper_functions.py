from sage.graphs.graph import Graph

from collections import Counter


### General helper functions

def get_node_index(node):
    return node[0]

def get_node_content(node):
    return node[1]

def node_changes(labelled_TD):
    r"""
    Record introduced and forgotten nodes in a directed labelled nice tree decomposition.

    INPUT:

    - ``labelled_TD`` -- a directed labelled nice tree decomposition,
      with the root as source of the dirgraph

    OUTPUT:

    - A dictionary of recorded nodes, where the `key` is node index, and
      the `value` is the introduced/forgotten node

    EXAMPLES::

        sage: from sage.graphs.graph_decompositions.tree_decomposition import label_nice_tree_decomposition
        sage: from sage.graphs.hom_count import node_changes
        sage: bip_one_four = graphs.CompleteBipartiteGraph(1, 4)
        sage: nice_tree_decomp = bip_one_four.treewidth(certificate=True, nice=True)
        sage: root = sorted(nice_tree_decomp)[0]
        sage: dir_labelled_TD = label_nice_tree_decomposition(nice_tree_decomp, root, directed=True)
        sage: node_changes(dir_labelled_TD)
        {1: 1, 2: 1, 3: 4, 5: 4, 6: 4, 7: 3, 8: 2, 9: 0, 10: 0, 11: 3, 12: 2}
    """
    node_changes_dict = {}

    for node in sorted(labelled_TD):
        node_index, node_vertex_set = node

        node_type = labelled_TD.get_vertex(node)
        match node_type:
            case 'intro':
                child_vertex_set = labelled_TD.neighbors_out(node)[0][1]
                # Get one element from the one-element set
                (extra_vertex,) = node_vertex_set.symmetric_difference(child_vertex_set)
                node_changes_dict[node_index] = extra_vertex
            case 'forget':
                child_vertex_set = labelled_TD.neighbors_out(node)[0][1]
                # Get one element from the one-element set
                (extra_vertex,) = node_vertex_set.symmetric_difference(child_vertex_set)
                node_changes_dict[node_index] = extra_vertex

    return node_changes_dict

def is_valid_mapping(mapped_vtx, mapped_nbhrs, target_graph):
    r"""
    Check if the mapping is valid.
    """
    if isinstance(target_graph, Graph):
        return all(target_graph.has_edge(mapped_vtx, vtx) for vtx in mapped_nbhrs)
    else:
        # Assume that `target_graph` is the adjacency matrix
        return all(target_graph[mapped_vtx, vtx] for vtx in mapped_nbhrs)

def count_occurrences(lst):
    r"""
    Count the occurrences of each element in `lst`
    """
    return Counter(lst)

def list_intersection(lst1, lst2):
    r"""
    Return the intersection of `lst1` and `lst2`
    """
    return list(set(lst1) & set(lst2))


### For integer representation

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


### For integer representation in colourful case

def encode_clr_list(clr_list, base):
    r"""
    Converts a list of integers to an integer in base-k representation.
    """
    return sum(val * base**idx for idx, val in enumerate(clr_list))

def decode_clr_int(num, base, nth):
    r"""
    Retrieve the nth element from the base-k representation.
    """
    num //= base ** nth
    return num % base
