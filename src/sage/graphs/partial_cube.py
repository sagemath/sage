r"""
Partial cubes

The code in this module that recognizes partial cubes is originally
from the PADS library by David Eppstein, which is available at
http://www.ics.uci.edu/~eppstein/PADS/ under the MIT license. It has a
quadratic runtime and has been described in [Epp2008]_.

For more information on partial cubes, see the :wikipedia:`Partial cube`.

Recognition algorithm
---------------------

Definitions
^^^^^^^^^^^

A **partial cube** is an isometric subgraph `G` of a
:meth:`~sage.graphs.graph_generators.GraphGenerators.CubeGraph` (of
possibly high dimension). Consequently, the vertices of `G` can be
labelled with binary sequences in such a way that the distance between
two vertices `u,v\in G` is the Hamming distance between their labels.

**Tokens** and their **action**: in the terminology of
[Epp2008]_, a token represents a transition of the form:

    *switch the k-th bit of the binary string from 0 to 1*

Each token can be matched with a 'reversed' token that performs the
same switch in the opposite direction. Alternatively, a token can be
seen as a set of disjoint (directed) edges of `G`, corresponding to
the transitions. When a vertex `v\in G` is the source of such an edge,
it is said that the token *acts* on `v`.

Observations
^^^^^^^^^^^^

**Shortest paths**: in a hypercube, a shortest path between two
vertices uses each token at most once. Furthermore, it cannot use both
a token and it reverse.

**Cycles**: a cycle in a partial cube is necessarily even, as
hypercubes are bipartite. If an edge `e` of a cycle `C` belongs to a
token `T`, then the edge opposite to `e` in `C` belongs to the reverse
of `T`.

**Incident edges**: all `2d_G(v)` arcs incident to a given vertex
belong to as many different tokens.

Algorithm
^^^^^^^^^

**Labeling**: Iteratively, the algorithm selects a vertex `v\in G`,
which is naturally associated to `2d(v)` tokens. It then performs a
breadth-first search from `v`, applying the previous observation on
cycles to attribute a token to some of the edges it meets. None of the
edges whose token remains undecided after this step can belong to one
of those `2d(v)` tokens, by virtue of the observation on shortest
paths.

The labeled edges can then be simplified (contracted) if the previous
step did not lead to a contradiction, and the procedure is applied
again until the graph is contracted to a single vertex and all edges
are labeled.

A partial cube is correctly labeled at this step, but some other
graphs can also satisfy the procedure.

**Checking the labeling**: once all tokens are defined and the
vertices are labeled with a binary string, we check that they define
an isometric subgraph of the hypercube. To ensure that the distance
`d(v_0,u)` is what we expect for any vertex `u`, it is sufficient to
find, for any vertex `u`, a neighbor `n_u` of `u` whose Hamming
distance with `v_0` is strictly less than the Hamming distance between
`u` and `v_0`. Here is the algorithm used to check the labeling:

* For an initial vertex `v`, run a BFS starting from `v`, and
  associate to every other vertex `u` a token that brings `u` closer
  to `v`. This yields shortest paths from every vertex to `v`.

* Assuming that the information is computed (and correct) for `v`, it
  is easy to update it for a neighbor `v'` of `v`. Indeed, if we write
  `T` the token that turns `v` into `v'`, only the vertices which were
  associated with the reverse of `T` need to select a new neighbour. All
  others can remain as they were previously.

  With this second observation, one can efficiently check that the
  distance between all pairs of vertices are what they should be. In
  the implementation, the sequence of the sources `(v, v', ...)` is
  given by a depth-first search.

Functions
---------
"""

def breadth_first_level_search(G, start):
    r"""
    Generate a sequence of dictionaries, each mapping the vertices at
    distance ``i`` from ``start`` to the set of their neighbours at
    distance ``i+1``.

    Originally written by D. Eppstein for the PADS library
    (http://www.ics.uci.edu/~eppstein/PADS/).

    INPUT:

    - ``G`` -- a graph to perform the search on.

    - ``start`` -- vertex or list of vertices from which to start the traversal.

    EXAMPLES::

        sage: H = digraphs.DeBruijn(3,2)
        sage: list(sage.graphs.partial_cube.breadth_first_level_search(H, '00'))
        [{'00': {'01', '02'}},
         {'01': {'10', '11', '12'}, '02': {'20', '21', '22'}},
         {'10': set(),
          '11': set(),
          '12': set(),
          '20': set(),
          '21': set(),
          '22': set()}]

    """
    neighbors = G.neighbor_out_iterator
    visited = set()
    if isinstance(start, list):
        currentLevel = start
    else:
        currentLevel = [start]
    while currentLevel:
        visited.update(currentLevel)
        nextLevel = set()
        levelGraph = {v: set() for v in currentLevel}
        for v in currentLevel:
            for w in neighbors(v):
                if w not in visited:
                    levelGraph[v].add(w)
                    nextLevel.add(w)
        yield levelGraph
        currentLevel = nextLevel

def depth_first_traversal(G, start):
    r"""
    Generate a sequence of triples (v,w,edgetype) for DFS of graph G.

    Originally written by D. Eppstein for the PADS library
    (http://www.ics.uci.edu/~eppstein/PADS/).

    INPUT:

    - ``G`` -- a graph to perform the search on.

    - ``start`` -- vertex or list of vertices from which to start the traversal.

    OUTPUT:

    - a generator of triples ``(v,w,edgetype)``, where ``edgetype`` is ``True``
      if the algorithm is progressing via the edge ``vw``, or ``False`` if the
      algorithm is backtracking via the edge ``wv``.

    EXAMPLES::

        sage: H = digraphs.DeBruijn(3,2)
        sage: t = list(sage.graphs.partial_cube.depth_first_traversal(H, '00'))
        sage: len(t)
        16

    """
    neighbors = G.neighbor_out_iterator
    seen = set()
    if not isinstance(start, list):
        start = [start]

    for v in start:
        if v in seen:
            continue
        seen.add(v)
        stack = [(v, neighbors(v))]
        while stack:
            parent, children = stack[-1]
            try:
                child = next(children)
                if child not in seen:
                    yield (parent, child, True)
                    seen.add(child)
                    stack.append((child, neighbors(child)))
            except StopIteration:
                stack.pop()
                if stack:
                    yield (stack[-1][0], parent, False)

def is_partial_cube(G, certificate=False):
    r"""
    Test whether the given graph is a partial cube.

    A partial cube is a graph that can be isometrically embedded into a
    hypercube, i.e., its vertices can be labelled with (0,1)-vectors of some
    fixed length such that the distance between any two vertices in the graph
    equals the Hamming distance of their labels.

    Originally written by D. Eppstein for the PADS library
    (http://www.ics.uci.edu/~eppstein/PADS/), see also
    [Epp2008]_.  The algorithm runs in `O(n^2)` time, where `n`
    is the number of vertices. See the documentation of
    :mod:`~sage.graphs.partial_cube` for an overview of the algorithm.

    INPUT:

    - ``certificate`` -- boolean (default: ``False``); this function returns
      ``True`` or ``False`` according to the graph, when ``certificate =
      False``. When ``certificate = True`` and the graph is a partial cube, the
      function returns ``(True, mapping)``, where ``mapping`` is an isometric
      mapping of the vertices of the graph to the vertices of a hypercube
      ((0, 1)-strings of a fixed length). When ``certificate = True`` and the
      graph is not a partial cube, ``(False, None)`` is returned.

    EXAMPLES:

    The Petersen graph is not a partial cube::

        sage: g = graphs.PetersenGraph()
        sage: g.is_partial_cube()
        False

    All prisms are partial cubes::

        sage: g = graphs.CycleGraph(10).cartesian_product(graphs.CompleteGraph(2))
        sage: g.is_partial_cube()
        True

    TESTS:

    The returned mapping is an isometric embedding into a hypercube::

        sage: g = graphs.DesarguesGraph()
        sage: _, m = g.is_partial_cube(certificate=True)
        sage: m # random
        {0: '00000',
         1: '00001',
         2: '00011',
         3: '01011',
         4: '11011',
         5: '11111',
         6: '11110',
         7: '11100',
         8: '10100',
         9: '00100',
         10: '01000',
         11: '10001',
         12: '00111',
         13: '01010',
         14: '11001',
         15: '10111',
         16: '01110',
         17: '11000',
         18: '10101',
         19: '00110'}
        sage: all(all(g.distance(u, v) == len([i for i in range(len(m[u])) if m[u][i] != m[v][i]]) for v in m) for u in m)
        True

    A graph without vertices is trivially a partial cube::

        sage: Graph().is_partial_cube(certificate=True)
        (True, {})

    """
    G._scream_if_not_simple()

    if not G.order():
        if certificate:
            return (True, {})
        else:
            return True

    if certificate:
        fail = (False, None)
    else:
        fail = False

    if not G.is_connected():
        return fail
    n = G.order()

    # Initial sanity check: are there few enough edges?
    # Needed so that we don't try to use union-find on a dense
    # graph and incur superquadratic runtimes.
    if 1 << (2*G.size()//n) > n:
        return fail

    # Check for bipartiteness.
    # This ensures also that each contraction will be bipartite.
    if not G.is_bipartite():
        return fail

    # Set up data structures for algorithm:
    # - contracted: contracted graph at current stage of algorithm
    # - unionfind: union find data structure representing known edge equivalences
    # - available: limit on number of remaining available labels
    from sage.graphs.digraph import DiGraph
    from sage.graphs.graph import Graph
    from sage.sets.disjoint_set import DisjointSet
    contracted = DiGraph({v: {w: (v, w) for w in G[v]} for v in G})
    unionfind = DisjointSet(contracted.edges(labels=False))
    available = n - 1

    # Main contraction loop in place of the original algorithm's recursion
    while contracted.order() > 1:
        # Find max degree vertex in contracted, and update label limit
        deg, root = max([(contracted.out_degree(v), v) for v in contracted], key=lambda x: x[0])
        if deg > available:
            return fail
        available -= deg

        # Set up bitvectors on vertices
        bitvec = {v: 0 for v in contracted}
        neighbors = {}
        for i, neighbor in enumerate(contracted[root]):
            bitvec[neighbor] = 1 << i
            neighbors[1 << i] = neighbor

        # Breadth first search to propagate bitvectors to the rest of the graph
        for level in breadth_first_level_search(contracted, root):
            for v in level:
                for w in level[v]:
                    bitvec[w] |= bitvec[v]

        # Make graph of labeled edges and union them together
        labeled = Graph([contracted.vertices(), []])
        for v, w in contracted.edge_iterator(labels=False):
            diff = bitvec[v]^bitvec[w]
            if not diff or not bitvec[w] &~ bitvec[v]:
                continue    # zero edge or wrong direction
            if diff not in neighbors:
                return fail
            neighbor = neighbors[diff]
            unionfind.union(contracted.edge_label(v, w),
                            contracted.edge_label(root, neighbor))
            unionfind.union(contracted.edge_label(w, v),
                            contracted.edge_label(neighbor, root))
            labeled.add_edge(v, w)

        # Map vertices to components of labeled-edge graph
        component = {}
        for i, SCC in enumerate(labeled.connected_components()):
            for v in SCC:
                component[v] = i

        # generate new compressed subgraph
        newgraph = DiGraph()
        for v, w, t in contracted.edge_iterator():
            if bitvec[v] == bitvec[w]:
                vi = component[v]
                wi = component[w]
                if vi == wi:
                    return fail
                if newgraph.has_edge(vi, wi):
                    unionfind.union(newgraph.edge_label(vi, wi), t)
                else:
                    newgraph.add_edge(vi, wi, t)
        contracted = newgraph

    # Make a digraph with edges labeled by the equivalence classes in unionfind
    g = DiGraph({v: {w: unionfind.find((v, w)) for w in G[v]} for v in G})

    # Associates to a vertex the token that acts on it, an check that
    # no two edges on a single vertex have the same label
    action = {}
    for v in g:
        action[v] = set(t for _, _, t in g.edge_iterator(v))
        if len(action[v]) != g.out_degree(v):
            return fail

    # Associate every token to its reverse
    reverse = {}
    for v, w, t in g.edge_iterator():
        rt = g.edge_label(w, v)
        reverse[t] = rt
        reverse[rt] = t

    current = initialState = next(g.vertex_iterator())

    # A token T is said to be 'active' for a vertex u if it takes u
    # one step closer to the source in terms of distance. The 'source'
    # is initially 'initialState'. See the module's documentation for
    # more explanations.

    # Find list of tokens that lead to the initial state
    activeTokens = set()
    for level in breadth_first_level_search(g, initialState):
        for v in level:
            for w in level[v]:
                activeTokens.add(g.edge_label(w, v))
    for t in activeTokens:
        if reverse[t] in activeTokens:
            return fail
    activeTokens = list(activeTokens)

    # Rest of data structure: point from states to list and list to states
    state_to_active_token = {v: -1 for v in g}
    token_to_states = [[] for i in activeTokens] # (i.e. vertices on which each token acts)

    def scan(v):
        """Find the next token that is effective for v."""
        a = next(i for i in range(state_to_active_token[v]+1, len(activeTokens))
                 if activeTokens[i] is not None
                    and activeTokens[i] in action[v])
        state_to_active_token[v] = a
        token_to_states[a].append(v)

    # Initialize isometric embedding into a hypercube
    if certificate:
        dim = 0
        tokmap = {}
        for t in reverse:
            if t not in tokmap:
                tokmap[t] = tokmap[reverse[t]] = 1 << dim
                dim += 1
        embed = {initialState: 0}

    # Set initial active states
    for v in g:
        if v != current:
            try:
                scan(v)
            except StopIteration:
                return fail

    # Traverse the graph, maintaining active tokens
    for prev, current, fwd in depth_first_traversal(g, initialState):
        if not fwd:
            prev, current = current, prev
        elif certificate:
            embed[current] = embed[prev] ^ tokmap[g.edge_label(prev, current)]

        # Add token to end of list, point to it from old state
        activeTokens.append(g.edge_label(prev, current))
        state_to_active_token[prev] = len(activeTokens) - 1
        token_to_states.append([prev])

        # Deactivate reverse token, find new token for its states
        #
        # (the 'active' token of 'current' is necessarily the label of
        #  (current, previous))
        activeTokens[state_to_active_token[current]] = None
        for v in token_to_states[state_to_active_token[current]]:
            if v != current:
                try:
                    scan(v)
                except StopIteration:
                    return fail

    # All checks passed, return the result
    if certificate:
        format = "{0:0%db}" % dim
        return (True, {v: format.format(l) for v, l in embed.items()})
    else:
        return True
