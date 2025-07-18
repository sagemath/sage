from dataclasses import dataclass
from sage.graphs.graph import Graph
from sage.graphs.views import EdgesView
from typing import cast, Any, List, Tuple, Hashable, NamedTuple

Edge = Tuple[Hashable, Hashable, Any]


def get_maximum_cardinality_matching(G: Graph) -> EdgesView:
    r"""
    Compute a maximum cardinality matching in a simple undirected graph using the Micali-Vazirani algorithm.

    INPUT:
    - ``G`` -- a simple undirected graph

    OUTPUT:
    - An `EdgesView` of the maximum cardinality matching in `G`
    """
    if G.has_loops() or G.has_multiple_edges():
        raise ValueError("Micali-Vazirani algorithm is only applicable to simple undirected graphs")

    # Return Empty EdgesView if G has no edges
    if not G.size():
        return EdgesView(Graph())

    # Global constants
    INFINITY = float('inf')

    @dataclass
    class Petal(NamedTuple):
        base: Hashable
        peaks: Tuple[Hashable, Hashable]

    # ******************************
    # Initialization of data structures
    # ******************************
    def initialize() -> None:
        global search_level_vertices, H, M
        for vertex in H:
            if not H.degree(vertex):
                H.delete_vertex(vertex)
                continue

            deletion_phase[vertex] = -1
            vertex_petal_map[vertex] = None
            vertex_bud_map[vertex] = vertex
            level[vertex] = [0, INFINITY] # (even_level, odd_level)
            min_level[vertex] = 0
            max_level[vertex] = INFINITY
            predecessor[vertex] = []
            successor[vertex] = []
            color[vertex] = None
            visit_mark[vertex] = None
            search_level_vertices.append(vertex)

        for (u, v, l) in H.edge_iterator():
            edge_scanned[(u, v, l)], edge_scanned[(v, u, l)] = -1, -1
            is_prop[(u, v, l)], is_prop[(v, u, l)] = None, None

        for index in range(1, H.order() + 2):
            tenacity_bridges_map[index] = []

    # *************************************
    # Greedy initial maximal matching (so as to reduce the total number of phases)
    # *************************************
    def compute_initial_maximal_matching() -> None:
        global matching_cardinality, M

        # 1) Make a mutable copy J of H
        J = H.copy()

        # 2) Initialize the degree map and the <degree, vertex_set> bucket structre
        degree_map = {vertex: J.degree(vertex) for vertex in J}
        maximum_degree = max(cast(list[int], degree_map.values())) if degree_map else 0

        buckets = [set() for _ in range(maximum_degree + 1)]

        # 3) Start at the smallest possible positive degree
        minimum_degree = 1
        while minimum_degree <= maximum_degree and not buckets[minimum_degree]:
            minimum_degree += 1

        # 4) Main loop
        while minimum_degree <= maximum_degree:
            # 4.0) If there are no vertices with minimum degree, increment the minimum degree
            if not buckets[minimum_degree]:
                minimum_degree += 1
                continue

            # 4.1) Get a vertex u of minimum positive degree
            u = buckets[minimum_degree].pop()

            # 4.2) Find a neighbor v of u with minimum degree
            v = min(J.neighbors(u), key=lambda x: J.degree(x))

            # 4.3) Add the edge (u, v) to the matching M
            M.add_edge(u, v, J.edge_label(u, v))

            # 4.4) Update the degree map and buckets
            vertices_to_remove = {u, v}
            vertices_to_remove.update(J.neighbors(u))
            vertices_to_remove.update(J.neighbors(v))

            for vertex in vertices_to_remove:
                if vertex in degree_map:
                    del degree_map[vertex]

                if J.degree(vertex) in buckets:
                    if vertex in buckets[J.degree(vertex)]:
                        buckets[J.degree(vertex)].remove(vertex)

            # 4.5) Gather all other neighbors whose degree will drop by 1
            vertices_to_update = set()
            for vertex in vertices_to_remove:
                vertices_to_update.update(J.neighbors(vertex))

            for vertex in vertices_to_update:
                if vertex in degree_map:
                    degree_map[vertex] -= 1
                    if degree_map[vertex] > 0:
                        buckets[J.degree(vertex)].remove(vertex)
                        buckets[J.degree(vertex) - 1].add(vertex)

            # 4.6) Update J by removing u, v, and their neighbors
            J.delete_vertices(vertices_to_remove)

            # 4.7) Note that we might have a vertex of minimum degree in J now, so continue the loop

    # ******************************
    # Start a new phase
    # ******************************
    def start_new_phase() -> None:
        global search_level_vertices, H, M
        search_level_vertices = []

        for u in H:
            if u in M:
                min_level[u] = INFINITY
                max_level[u] = INFINITY
                level[u] = [INFINITY, INFINITY]

            else:
                search_level_vertices.append(u)
                min_level[u] = 0
                max_level[u] = INFINITY
                level[u] = [0, INFINITY]

            predecessor[u] = []
            successor[u] = []
            vertex_petal_map[u] = None
            vertex_bud_map[u] = u
            color[u] = None
            visit_mark[u] = None

        for (u, v, l) in H.edge_iterator():
            is_prop[(u, v, l)], is_prop[(v, u, l)] = None, None
            edge_scanned[(u, v, l)], edge_scanned[(v, u, l)] = -1, -1

        for index in range(1, int(2*H.order()+1)):
            tenacity_bridges_map[index] = []

    # ******************************
    # Primary Subroutine: Find min_level of vertices
    # ******************************
    def MIN(search_level: int) -> bool:
        global search_level_vertices, phase_index, H, M
        next_search_level_vertices = []
        parity = search_level % 2

        if not search_level_vertices or search_level > H.order():
            return True

        for u in search_level_vertices:
            if level[u][parity] != search_level and level[u][parity] < INFINITY:
                next_search_level_vertices.append(u)
                continue

            if deletion_phase[u] == phase_index:
                continue

            for v in H.neighbor_iterator(u):
                l = H.edge_label(u, v)

                if edge_scanned[(u, v, l)] != phase_index and M.has_edge(u, v, l) == parity and deletion_phase[v] != phase_index:
                    edge_scanned[(u, v, l)], edge_scanned[(v, u, l)] = phase_index, phase_index

                    if min_level[v] > search_level:
                        min_level[v] = search_level + 1
                        level[v][1 - parity] = search_level + 1
                        next_search_level_vertices.append(v)
                        predecessor[v].append(u)
                        successor[u].append(v)
                        is_prop[(u, v, l)], is_prop[(v, u, l)] = True, True

                    else:
                        tenacity = level[u][parity]+level[v][parity]+1

                        # In the case where tenacity is defined and thus we know which level the bridge will be processed
                        if tenacity < INFINITY:
                            tenacity_bridges_map[tenacity].append((u, v, l))

                        # The case where tenacity is not yet known (possibly due to the even/ odd level of the blossom not yet labeled
                        else:
                            is_prop[(u, v, l)], is_prop[(v, u, l)] = False, False

        search_level_vertices = next_search_level_vertices
        return False

    # ******************************
    # Primary Subroutine: Find max_level of vertices
    # ******************************
    def MAX(search_level: int) -> bool:
        global H, M, phase_index, num_augmentations, previous_search_level
        is_augmented = False

        for (u, v, l) in tenacity_bridges_map[2*search_level + 1]:
            if deletion_phase[u] == phase_index or deletion_phase[v] == phase_index:
                continue

            left_support, right_support, bottleneck, encountered_deleted_vertex = DDFS(u, v)

            # if the bridge has been augmented
            if bottleneck is None:
                if not encountered_deleted_vertex:
                    augmentation_success = augment(left_support, right_support, (u, v, l), search_level)
                    if augmentation_success:
                        is_augmented = True
                        if M.size() == H.order() / 2:
                            return is_augmented

            else:
                if not encountered_deleted_vertex:
                    form_blossom(left_support, right_support, bottleneck, (u, v, l))
                    label_max(left_support, search_level)
                    label_max(right_support, search_level)

        if is_augmented:
            previous_search_level = search_level
            num_augmentations += 1
            is_augmented = False
        return is_augmented

    # After identifying support, we assign the vertex its max level label.
    # Note: This step is skipped during augmentation, as max levels are reset regardless.
    # ******************************
    # Label vertices after forming a blossom
    # ******************************
    def label_max(support: List[Hashable], search_level: int) -> None:
        global search_level_vertices, H

        next_search_level_vertices = []
        for vertex in support:
            max_level[vertex] = 2*search_level + 1 - min_level[vertex]
            level_parity = max_level[vertex] % 2
            level[vertex][level_parity] = max_level[vertex] % 2
            next_search_level_vertices.append(vertex)

            if not level_parity:
                for neighbor in H.neighbor_iterator(vertex):
                    edge = (vertex, neighbor, H.edge_label(vertex, neighbor))

                    # In the case were the tenacity of a tenacity_bridges_map was not yet found
                    if not is_prop[edge]:
                        tenacity_bridges_map[max_level[vertex] + level[neighbor][0] + 1].append(edge)
        search_level_vertices += next_search_level_vertices

    # ******************************
    # Double DFS to locate augmenting paths
    # ******************************
    def DDFS(source_red_vertex: Hashable, source_green_vertex: Hashable) -> Tuple[List[Hashable], List[Hashable], Hashable, bool]:
        global phase_index
        encountered_deleted_vertex = False

        # Set the starting point for each of red and green DFS's
        red_stack, green_stack = [], []  # Stack saves previously traversed vertices
        red_vertex, green_vertex = get_bud(source_red_vertex), get_bud(source_green_vertex)  # Set the initial point for both DFS's

        red_predecessors, green_predecessors = \
            predecessor[red_vertex][:], predecessor[green_vertex][:]  # Copy predecessor list over for the current vertex
        red_support, green_support = [red_vertex], [green_vertex]  # the lists holding the support of the current bridge

        # Following is used to save the data for DFS's for when they backtrack in the case a bottleneck is reached
        previous_red_support, previous_green_support = [red_vertex], [green_vertex]

        # Boolean variables are initiated
        no_augmentation_found = False if not min_level[red_vertex] and not min_level[green_vertex] else True
        collision = True if red_vertex == green_vertex else False

        # Returns nothing if there is no support for the petal
        if collision and not no_augmentation_found:
            return [], [], red_vertex, encountered_deleted_vertex

        # Label is used to track if vertices have been visit_mark in the current DDFS
        label = (source_red_vertex, source_green_vertex)
        visit_mark[red_vertex], visit_mark[green_vertex] = label, label

        # DDFS continues to run while an augmenting path still isn't found
        while no_augmentation_found:

            # Checks for when the two DFS's land on the same vertex
            if collision:

                # The the levels of the vertices are the same, we reverse the green DFS
                if min_level[red_vertex] == min_level[green_vertex]:
                    previous_green_support = green_support[:]
                    green_vertex, green_predecessors, reverse_check = reverse_DFS(green_vertex, green_predecessors, green_stack, green_support)

                    if reverse_check:
                        previous_red_support, red_bottleneck = red_support[:], red_vertex
                        red_vertex, red_predecessors, reverse_check = reverse_DFS(red_vertex, red_predecessors, red_stack, red_support)

                elif min_level[red_vertex] > min_level[green_vertex]:
                    red_vertex, red_predecessors, reverse_check = reverse_DFS(red_vertex, red_predecessors, red_stack, red_support)

                elif min_level[red_vertex] < min_level[green_vertex]:
                    green_vertex, green_predecessors, reverse_check = reverse_DFS(green_vertex, green_predecessors, green_stack, green_support)

                if red_vertex == green_vertex:
                    previous_red_support.pop()
                    green_support.pop()
                    return previous_red_support, green_support, red_vertex, encountered_deleted_vertex

                collision = False

            # Case where red DFS advances in search
            elif min_level[red_vertex] >= min_level[green_vertex]:

                # Advance the red DFS, will reverse if no vertices to travel to
                red_vertex, red_predecessors, collision = advance_DFS(red_vertex, red_predecessors, red_stack, red_support, label)

                # If stack is cleared and no vertices left to explore, bottleneck is found
                if not red_stack and not red_predecessors:
                    previous_red_support.pop()
                    green_support.pop()
                    return previous_red_support, green_support, green_vertex, encountered_deleted_vertex

            # Case where green DFS advances in search
            else:

                # Advance the green DFS, will reverse if no vertices to travel to
                green_vertex, green_predecessors, collision = advance_DFS(green_vertex, green_predecessors, green_stack, green_support, label)

                # If stack is clearned and no vertices left to explore, reverse red DFS
                if not green_stack and not green_predecessors:
                    green_support = previous_green_support
                    previous_green_support, green_vertex, green_predecessors = [green_vertex], red_vertex, red_predecessors[:]
                    previous_red_support, red_bottleneck = red_support[:], red_vertex
                    red_vertex, red_predecessors, reverse_check = reverse_DFS(red_vertex, red_predecessors, red_stack, red_support)

                    if reverse_check:
                        previous_red_support.pop()
                        green_support.pop()
                        return previous_red_support, green_support, red_bottleneck, encountered_deleted_vertex

            # Checks if vertex was removed in previous augmentation during current search_level
            if deletion_phase[red_vertex] == phase_index or deletion_phase[green_vertex] == phase_index:
                encountered_deleted_vertex = True

            # Checks if augmenting path has been found
            if not min_level[red_vertex] and not min_level[green_vertex] and red_vertex != green_vertex:
                no_augmentation_found = False

        return red_support, green_support, None, encountered_deleted_vertex

    # ******************************
    # Auxiliary Subroutine: advance DFS along predecessors
    # ******************************
    def advance_DFS(vertex: Hashable, predecessor_list: List[Hashable], stack: List[Tuple[Hashable, List[Hashable]]], support: List[Hashable], label: Tuple[Hashable, Hashable]) -> Tuple[Hashable, List[Hashable], bool]:
        reverse_check = False
        if predecessor_list:

            next_vertex = get_bud(predecessor_list.pop())
            # Save the previous vertex with it's predecessor list to the stack
            stack.append((vertex, predecessor_list))
            # Add next vertex to support
            support.append(next_vertex)
            predecessor_list = predecessor[next_vertex][:]

            if visit_mark[next_vertex] == label:
                return next_vertex, predecessor_list, True
            visit_mark[next_vertex] = label
        # If next vertex not found reverse path
        else:
            next_vertex, predecessor_list, reverse_check = reverse_DFS(vertex, predecessor_list, stack, support)
        return next_vertex, predecessor_list, reverse_check

    # ******************************
    # Auxiliary Subroutine: backtrack in DFS stack
    # ******************************
    def reverse_DFS(vertex: Hashable, predecessor_list: List[Hashable], stack: List[Tuple[Hashable, List[Hashable]]], support: List[Hashable]) -> Tuple[Hashable, List[Hashable], bool]:
        failure = False
        if stack:
            previous_vertex = stack.pop()
            vertex = previous_vertex[0]
            predecessor_list = previous_vertex[1]
            support.pop()
        else:
            failure = True
        return vertex, predecessor_list, failure

    # Each vertex can only belong to one petal
    # The bud cannot be part of the petal
    # Each vertex in the petal points to the bud

    # ******************************
    # Contract a blossom (petal)
    # ******************************
    def form_blossom(left_support: List[Hashable], right_support: List[Hashable], bud: Hashable, bridge: Edge):
        """
        Create a new blossom centered at 'bud' with supports from both sides.
        """
        global phase_index
        petal_ = Petal(base=bud, peaks=(bridge[0], bridge[1]))
        form_petal(left_support, bud, petal_, 0)
        form_petal(right_support, bud, petal_, 1)

    def form_petal(support: List[Hashable], bud: Hashable, petal: Petal, direction: int):
        """
        Assign each vertex in support to the given petal and direction.
        """
        for vertex in support:
            vertex_bud_map[vertex] = get_bud(bud)
            vertex_petal_map[vertex] = petal
            color[vertex] = direction

    # ******************************
    # Path compression: find the bud of a vertex
    # ******************************
    def get_bud(vertex: Hashable) -> Hashable:
        if vertex != vertex_bud_map[vertex]:
            vertex_bud_map[vertex] = get_bud(vertex_bud_map[vertex])
        return vertex_bud_map[vertex]

    # ******************************
    # Unfold a blossom (petal)
    # ******************************
    # unfolds a petal to get the path from the vertex that is part of the petal to the bud
    def unfold_petal(vertex: Hashable, target: Hashable) -> List[Hashable]:
        path = list()
        petal = vertex_petal_map[vertex]
        bud = petal.base
        if max_level[vertex] % 2:
            path = unfold_path_in_petal(vertex, bud, petal)
        else:
            red_vertex = petal.peaks[0]
            green_vertex = petal.peaks[1]
            if not color[vertex]:
                left_path = unfold_path_in_petal(red_vertex, vertex, petal)
                right_path = unfold_path_in_petal(green_vertex, bud, petal)
                if left_path and right_path:
                    left_path.reverse()
                    path = left_path + right_path
                else:
                    return []
            elif color[vertex]:
                left_path = unfold_path_in_petal(red_vertex, bud, petal)
                right_path = unfold_path_in_petal(green_vertex, vertex, petal)
                if left_path and right_path:
                    right_path.reverse()
                    path = right_path + left_path
                else:
                    return []
        if bud == target:
            return path
        else:
            path.pop()
            petal_path = unfold_petal(bud, target)
            return path + petal_path

    def unfold_path_in_petal(start_vertex: Hashable, end_vertex: Hashable, petal: Petal) -> List[Hashable]:
        global M
        if start_vertex == end_vertex:
            return [start_vertex]
        if vertex_petal_map[start_vertex] != petal:
            new_target = vertex_petal_map[start_vertex].base
            path = unfold_petal(start_vertex, new_target)
            current_vertex = new_target
        else:
            path = [start_vertex]
            current_vertex = start_vertex
        while current_vertex != end_vertex:
            previous_vertex = current_vertex
            predecessor_list = predecessor[current_vertex][:]
            new_petal = None
            next_petal_vertex = None
            wrong_petal_vertex = None
            for vertex in predecessor_list:
                if vertex_petal_map[vertex] is not None:
                    if vertex == end_vertex:
                        current_vertex = vertex
                        path.append(end_vertex)
                        break
                    elif vertex_petal_map[vertex] == petal and color[current_vertex] == color[vertex]:
                        next_petal_vertex = vertex
                    elif vertex_petal_map[vertex] == petal and color[current_vertex] != color[vertex]:
                        wrong_petal_vertex = vertex
                    else:
                        new_petal = vertex
                else:
                    if vertex == petal.base:
                        current_vertex = vertex
                        path.append(end_vertex)
                        break
            if previous_vertex == current_vertex:
                if next_petal_vertex is not None:
                    current_vertex = next_petal_vertex
                    path.append(current_vertex)
                elif wrong_petal_vertex is not None:
                    current_vertex = vertex
                    if vertex_petal_map[vertex] != petal and vertex_petal_map[vertex] is not None:
                        petal_path = unfold_petal(current_vertex, vertex_petal_map[vertex].base)
                        path += petal_path
                        current_vertex = path[-1]
                    else:
                        path.append(current_vertex)
                elif new_petal is None:
                    return []
                else:
                    if not M.has_edge(current_vertex, new_petal, H.edge_label(current_vertex, new_petal)):
                        path_addition = unfold_petal(new_petal, vertex_petal_map[new_petal].base)
                        if not path_addition:
                            return []
                        path += path_addition
                        current_vertex = vertex_petal_map[new_petal].base
                        while vertex_petal_map[current_vertex] != petal and current_vertex != end_vertex:
                            # digging deeper into the petal
                            path.pop()
                            if vertex_petal_map[current_vertex]:
                                path += unfold_petal(current_vertex, vertex_petal_map[current_vertex].base)
                                current_vertex = vertex_petal_map[current_vertex].base
                            else:
                                # bud failure
                                return []
                    else:
                        path += unfold_path_in_petal(new_petal, vertex_petal_map[new_petal].base, vertex_petal_map[new_petal])
                        current_vertex = vertex_petal_map[new_petal].base

        return path

    # ******************************
    # Augment along a found path
    # ******************************
    def augment(left_support: List[Hashable], right_support: List[Hashable], bridge: Edge, search_level: int) -> bool:
        """
        Augment the matching M with the agumenting path found.
        """
        global M, phase_index
        left_path = get_path(left_support, bridge[0])
        right_path = get_path(right_support, bridge[1])
        if not left_path or not right_path:
            return False
        left_path.reverse()
        path = left_path + right_path
        if len(path) != search_level * 2 + 2:
            return False # Error in augmentation path length

        for index in range(len(path) - 1):
            edge = (path[index], path[index+1], H.edge_label(path[index], path[index+1]))
            if M.has_edge(edge):
                M.delete_edge(edge)
            else:
                M.add_edge(edge)

        # Erase vertex based on search level
        erase_vertex_list = []
        for vertex in path:
            deletion_phase[vertex] = phase_index
            erase_vertex_list.append(vertex)

        while erase_vertex_list:
            current_vertex = erase_vertex_list.pop()
            successors = successor[current_vertex][:]
            for vertex in successors:
                if deletion_phase[vertex] != phase_index:
                    predecessor[vertex].remove(current_vertex)
                    successor[current_vertex].remove(vertex)
                    if not predecessor[vertex]:
                        erase_vertex_list.append(vertex)
                        deletion_phase[vertex] = phase_index
        return True

    # Procedure to find path after discovering an augmenting path via DDFS
    def get_path(support: List[Hashable], peak: Hashable) -> List[Hashable]:
        """
        Build the vertex sequence of an augmenting path.
        """
        path = []
        current_vertex = peak
        predecessor_list = predecessor[current_vertex][:]

        # Follow the support given to get the path
        for vertex in support:

            # Do a search following the trail given by the support
            while get_bud(current_vertex) != vertex:
                # If it is not the correct vertex pop the next vertex in the predecessor list
                current_vertex = predecessor_list.pop()

            predecessor_list = predecessor[current_vertex][:]
            # If the vertex is not part of a petal, add it to the path
            # Otherwise unfold the petal
            if vertex_petal_map[current_vertex] is None:
                path.append(current_vertex)
            else:
                petal_path = unfold_petal(current_vertex, get_bud(current_vertex))
                if not petal_path:
                    return []
                path += petal_path
                current_vertex = get_bud(current_vertex)
                predecessor_list = predecessor[current_vertex][:]
        return path

    # ******************************
    # Main: search phases
    # This is the main loop that finds and aguments phases (a maximal set of minimum length disjoin augmenting paths) and erase those vertices judiciously.
    # ******************************
    def search() -> bool:
        """
        Perform one complete search phase to find and augment any disjoint augmenting paths.
        Returns True if any augmentation was found.
        """
        global search_level, previous_search_level, num_augmentations
        search_level = 0
        previous_search_level = 0
        augmentation_found = False
        search_complete = False
        while not augmentation_found and not search_complete:
            search_complete = MIN(search_level)
            augmentation_found = MAX(search_level)
            if search_complete and not num_augmentations:
                return False

            search_level += 1
        return True

    # ******************************
    # Set up global state containers
    # ******************************
    tenacity_bridges_map = {}
    deletion_phase = {}
    visit_mark = {}
    vertex_petal_map = {}
    vertex_bud_map = {}
    level = {}
    min_level = {}
    max_level = {}
    predecessor = {}
    successor = {}
    color = {}
    edge_scanned = {}
    is_prop = {}

    global H, M, search_level_vertices, phase_index, num_augmentations
    H = G.copy()
    M = Graph()
    search_level_vertices = []
    phase_index = 0
    num_augmentations = 0

    # ******************************
    # Algorithm execution flow
    # ******************************

    there_exists_a_phase = True
    initialize()
    compute_initial_maximal_matching()

    start_new_phase()
    while there_exists_a_phase:
        num_augmentations = 0
        there_exists_a_phase = search()
        phase_index += 1
        start_new_phase()

        # Stop early if perfect matching found
        if M.size() == H.order() / 2:
            break

    return EdgesView(M)