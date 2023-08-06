# distutils: libraries = planarity
"""
Wrapper for Boyer's (C) planarity algorithm
"""

cdef extern from "planarity/graph.h":
    ctypedef struct vertexRec:
        int link[2]
        int index
    ctypedef vertexRec * vertexRecP

    ctypedef struct edgeRec:
        int link[2]
        int neighbor
    ctypedef edgeRec * edgeRecP

    ctypedef struct BM_graph:
        vertexRecP V
        edgeRecP E
        int N
    ctypedef BM_graph * graphP

    cdef int OK, EMBEDFLAGS_PLANAR, NONEMBEDDABLE, NOTOK

    cdef graphP gp_New()
    cdef void gp_Free(graphP *pGraph)
    cdef int gp_InitGraph(graphP theGraph, int N)
    cdef int gp_AddEdge(graphP theGraph, int u, int ulink, int v, int vlink)
    cdef int gp_Embed(graphP theGraph, int embedFlags)
    cdef int gp_SortVertices(graphP theGraph)


def is_planar(g, kuratowski=False, set_pos=False, set_embedding=False, circular=None):
    r"""
    Check whether ``g`` is planar using Boyer's planarity algorithm.

    If ``kuratowski`` is ``False``, returns ``True`` if ``g`` is planar,
    ``False`` otherwise.  If ``kuratowski`` is ``True``, returns a tuple, first
    entry is a boolean (whether or not the graph is planar) and second entry is
    a Kuratowski subgraph, i.e. an edge subdivision of `K_5` or `K_{3,3}` (if
    not planar) or ``None`` (if planar).  Also, will set an ``_embedding``
    attribute for the graph ``g`` if ``set_embedding`` is set to ``True``.

    INPUT:

    - ``kuratowski`` -- boolean (default: ``False``); when set to ``True``,
      return a tuple of a boolean and either ``None`` or a Kuratowski subgraph
      (i.e. an edge subdivision of `K_5` or `K_{3,3}`). When set to ``False``,
      returns ``True`` if ``g`` is planar, ``False`` otherwise.

    - ``set_pos`` -- boolean (default: ``False``); whether to use Schnyder's
      algorithm to determine and set positions

    - ``set_embedding`` -- boolean (default: ``False``); whether to record the
      combinatorial embedding returned (see
      :meth:`~sage.graphs.generic_graph.GenericGraph.get_embedding`)

    - ``circular`` -- deprecated argument

    EXAMPLES::

        sage: G = graphs.DodecahedralGraph()
        sage: from sage.graphs.planarity import is_planar
        sage: is_planar(G)
        True
        sage: Graph('@').is_planar()
        True

    TESTS:

    We try checking the planarity of all graphs on 7 or fewer
    vertices.  In fact, to try to track down a segfault, we do it
    twice. ::

        sage: # long time, needs networkx
        sage: import networkx.generators.atlas
        sage: atlas_graphs = [Graph(i) for i in networkx.generators.atlas.graph_atlas_g()]
        sage: a = [i for i in [1..1252] if atlas_graphs[i].is_planar()]
        sage: b = [i for i in [1..1252] if atlas_graphs[i].is_planar()]
        sage: a == b
        True

    There were some problems with ``set_pos`` stability in the past,
    so let's check if this runs without exception::

        sage: for i, g in enumerate(atlas_graphs):      # long time                     # needs networkx
        ....:     if (not g.is_connected() or i == 0):
        ....:         continue
        ....:     _ = g.is_planar(set_embedding=True, set_pos=True)

        Argument saving::

            sage: G = Graph([(1, 2)])
            sage: for set_embedding, set_pos in ((True,True), (True,False), (False, True), (False, False)):
            ....:     G = Graph([(1, 2)])
            ....:     assert is_planar(G, set_embedding=set_embedding, set_pos=set_pos)
            ....:     assert (hasattr(G, '_embedding') and G._embedding is not None) == set_embedding, (set_embedding, set_pos)
            ....:     assert (hasattr(G, '_pos') and G._pos is not None) == set_pos, (set_embedding, set_pos)

    """
    if circular is not None:
        from sage.misc.superseded import deprecation
        deprecation(33759, 'the circular argument of is_planar is deprecated and has no effect')

    if set_pos and not g.is_connected():
        raise ValueError("is_planar() cannot set vertex positions for a disconnected graph")

    # First take care of a trivial cases
    if not g.size():
        # There are no edges
        if set_embedding:
            g._embedding = {v: [] for v in g}
        return (True, None) if kuratowski else True
    if g.order() == 2 and g.is_connected():
        # P_2 is too small to be triangulated
        u, v = list(g)
        if set_embedding:
            g._embedding = {u: [v], v: [u]}
        if set_pos:
            g._pos = {u: [0, 0], v: [0, 1]}
        return (True, None) if kuratowski else True

    # Create to and from mappings to relabel vertices to the set {1,...,n}
    # (planarity 3 uses 1-based array indexing, with 0 representing NIL)
    cdef int i
    cdef list listto = list(g)
    cdef dict ffrom = {vvv: i + 1 for i, vvv in enumerate(listto)}
    cdef dict to = {i + 1: vvv for i, vvv in enumerate(listto)}
    g.relabel(ffrom)

    cdef graphP theGraph
    theGraph = gp_New()
    cdef int status
    status = gp_InitGraph(theGraph, g.order())
    if status != OK:
        raise RuntimeError("gp_InitGraph status is not ok")
    for u, v in g.edge_iterator(labels=False):
        status = gp_AddEdge(theGraph, u, 0, v, 0)
        if status == NOTOK:
            raise RuntimeError("gp_AddEdge status is not ok")
        elif status == NONEMBEDDABLE:
            # We now know that the graph is nonplanar.
            if not kuratowski:
                return False
            # With just the current edges, we have a nonplanar graph,
            # so to isolate a kuratowski subgraph, just keep going.
            break

    status = gp_Embed(theGraph, EMBEDFLAGS_PLANAR)
    gp_SortVertices(theGraph)

    # Use to and from mappings to relabel vertices back from the set {1,...,n}
    g.relabel(to)

    if status == NOTOK:
        raise RuntimeError("status is not ok")
    elif status == NONEMBEDDABLE:
        # Kuratowski subgraph isolator
        g_dict = {}
        from sage.graphs.graph import Graph
        for i in range(1, theGraph.N + 1):
            linked_list = []
            j = theGraph.V[i].link[1]
            while j:
                linked_list.append(to[theGraph.E[j].neighbor])
                j = theGraph.E[j].link[1]
            if linked_list:
                g_dict[to[i]] = linked_list
        G = Graph(g_dict)
        gp_Free(&theGraph)
        if kuratowski:
            return (False, G)
        else:
            return False
    else:
        if set_pos or set_embedding:
            emb_dict = {}
            for i in range(1, theGraph.N + 1):
                linked_list = []
                j = theGraph.V[i].link[1]
                while j:
                    linked_list.append(to[theGraph.E[j].neighbor])
                    j = theGraph.E[j].link[1]
                emb_dict[to[i]] = linked_list
            if set_embedding:
                g._embedding = emb_dict
            if set_pos:
                g.layout(layout='planar', save_pos=True, on_embedding=emb_dict)

        gp_Free(&theGraph)
        if kuratowski:
            return (True, None)
        else:
            return True
