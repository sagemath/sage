r"""
Backends for Sage (di)graphs.

This module implements :class:`GenericGraphBackend` (the base class for
backends) and :class:`NetworkXGraphBackend` (a wrapper for `NetworkX
<http://networkx.lanl.gov/>`__ graphs)

Any graph backend must redefine the following methods (for which
:class:`GenericGraphBackend` raises a ``NotImplementedError``)

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~GenericGraphBackend.add_edge` | Add an edge `(u,v)` to ``self``, with label `l`.
    :meth:`~GenericGraphBackend.add_edges` | Add a sequence of edges to ``self``.
    :meth:`~GenericGraphBackend.add_vertex` | Add a labelled vertex to ``self``.
    :meth:`~GenericGraphBackend.add_vertices` | Add labelled vertices to ``self``.
    :meth:`~GenericGraphBackend.degree` | Return the total number of vertices incident to `v`.
    :meth:`~GenericGraphBackend.in_degree` | Return the in-degree of `v`
    :meth:`~GenericGraphBackend.out_degree` | Return the out-degree of `v`
    :meth:`~GenericGraphBackend.del_edge` | Delete the edge `(u,v)` with label `l`.
    :meth:`~GenericGraphBackend.del_vertex` | Delete a labelled vertex in ``self``.
    :meth:`~GenericGraphBackend.del_vertices` | Delete labelled vertices in ``self``.
    :meth:`~GenericGraphBackend.get_edge_label` | Return the edge label of `(u,v)`.
    :meth:`~GenericGraphBackend.has_edge` | True if ``self`` has an edge `(u,v)` with label `l`.
    :meth:`~GenericGraphBackend.has_vertex` | True if ``self`` has a vertex with label `v`.
    :meth:`~GenericGraphBackend.iterator_edges` | Iterate over the edges incident to a sequence of vertices.
    :meth:`~GenericGraphBackend.iterator_in_edges` | Iterate over the incoming edges incident to a sequence of vertices.
    :meth:`~GenericGraphBackend.iterator_out_edges` | Iterate over the outbound edges incident to a sequence of vertices.
    :meth:`~GenericGraphBackend.iterator_nbrs` | Iterate over the vertices adjacent to `v`.
    :meth:`~GenericGraphBackend.iterator_in_nbrs` | Iterate over the vertices u such that the edge `(u,v)` is in ``self`` (that is, predecessors of `v`).
    :meth:`~GenericGraphBackend.iterator_out_nbrs` | Iterate over the vertices u such that the edge `(v,u)` is in ``self`` (that is, successors of `v`).
    :meth:`~GenericGraphBackend.iterator_verts` | Iterate over the vertices `v` with labels in verts.
    :meth:`~GenericGraphBackend.loops` | Get/set whether or not ``self`` allows loops.
    :meth:`~GenericGraphBackend.multiple_edges` | Get/set whether or not ``self`` allows multiple edges.
    :meth:`~GenericGraphBackend.name` | Get/set name of ``self``.
    :meth:`~GenericGraphBackend.num_edges` | The number of edges in ``self``
    :meth:`~GenericGraphBackend.num_verts` | The number of vertices in ``self``
    :meth:`~GenericGraphBackend.relabel` | Relabel the vertices of ``self`` by a permutation.
    :meth:`~GenericGraphBackend.set_edge_label` | Label the edge `(u,v)` by `l`.

For an overview of graph data structures in sage, see
:mod:`~sage.graphs.base.overview`.

Classes and methods
-------------------
"""

#*******************************************************************************
#        Copyright (C) 2008 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*******************************************************************************
from __future__ import absolute_import

from .c_graph cimport CGraphBackend
from .c_graph cimport CGraph

from six import itervalues


cdef class GenericGraphBackend(SageObject):
    """
    A generic wrapper for the backend of a graph.  Various graph classes use
    extensions of this class.  Note, this graph has a number of placeholder
    functions, so the doctests are rather silly.

    TESTS::

        sage: import sage.graphs.base.graph_backends

    """
    _loops = False
    _multiple_edges = False
    _name = ''
    def add_edge(self, u, v, l, directed):
        """
        Add an edge (u,v) to self, with label l.  If directed is True, this is
        interpreted as an arc from u to v.

        INPUT:

        - ``u,v`` -- vertices
        - ``l`` -- edge label
        - ``directed`` -- boolean

        TESTS::

            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.add_edge(1,2,'a',True)
            Traceback (most recent call last):
            ...
            NotImplementedError
         """
        raise NotImplementedError()
    def add_edges(self, edges, directed):
        """
        Add a sequence of edges to self.  If directed is True, these are
        interpreted as arcs.

        INPUT:

        - ``edges`` -- list/iterator of edges to be added.

        - ``directed`` -- boolean

        TESTS::

            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.add_edges([],True)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()

    def add_vertex(self, name):
        """
        Add a labelled vertex to self.

        INPUT:

        - ``name`` -- vertex label

        OUTPUT:

        If ``name=None``, the new vertex name is returned, ``None`` otherwise.

        TESTS::

            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.add_vertex(0)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def add_vertices(self, vertices):
        """
        Add labelled vertices to self.

        INPUT:

        - ``vertices`` -- iterator of vertex labels. A new label is created,
          used and returned in the output list for all ``None`` values in
          ``vertices``.

        OUTPUT:

        Generated names of new vertices if there is at least one ``None`` value
        present in ``vertices``. ``None`` otherwise.

        EXAMPLES::

            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.add_vertices([1,2,3])
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def degree(self, v, directed):
        """
        Return the total number of vertices incident to `v`.

        INPUT:

        - ``v`` -- a vertex label
        - ``directed`` -- boolean

        OUTPUT:

            degree of v

        TESTS::

            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.degree(1, False)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def in_degree(self, v):
        r"""
        Return the in-degree of `v`

        INPUT:

        - ``v`` -- a vertex label

        TESTS::

            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.in_degree(1)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def out_degree(self, v):
        r"""
        Return the out-degree of `v`

        INPUT:

        - ``v`` -- a vertex label

        TESTS::

            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.out_degree(1)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def del_edge(self, u, v, l, directed):
        """
        Delete the edge `(u,v)` with label `l`.

        INPUT:

        - ``u,v`` -- vertices
        - ``l`` -- edge label
        - ``directed`` -- boolean

        TESTS::

            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.del_edge(1,2,'a',True)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def del_vertex(self, v):
        """
        Delete a labelled vertex in self.

        INPUT:

        - ``v`` -- vertex label

        TESTS::

            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.del_vertex(0)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def del_vertices(self, vertices):
        """
        Delete labelled vertices in self.

        INPUT:

        - ``vertices`` -- iterator of vertex labels

        TESTS::

            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.del_vertices([1,2,3])
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def get_edge_label(self, u, v):
        """
        Return the edge label of `(u,v)`.

        INPUT:

        - ``u,v`` -- vertex labels

        OUTPUT:

            label of `(u,v)`

        TESTS::

            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.get_edge_label(1,2)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def has_edge(self, u, v, l):
        """
        True if self has an edge (u,v) with label l.

        INPUT:

        - ``u,v`` -- vertex labels
        - ``l`` -- label

        OUTPUT:

            boolean

        TESTS::

            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.has_edge(1,2,'a')
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def has_vertex(self, v):
        """
        True if self has a vertex with label v.

        INPUT:

        - ``v`` -- vertex label

        OUTPUT:
            boolean

        TESTS::

            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.has_vertex(0)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def iterator_edges(self, vertices, labels):
        """
        Iterate over the edges incident to a sequence of vertices. Edges are
        assumed to be undirected.

        INPUT:

        - ``vertices`` -- a list of vertex labels
        - ``labels`` -- boolean

        OUTPUT:

            a generator which yields edges, with or without labels
            depending on the labels parameter.

        TESTS::

            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.iterator_edges([],True)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def iterator_in_edges(self, vertices, labels):
        """
        Iterate over the incoming edges incident to a sequence of vertices.

        INPUT:

        - ``vertices`` -- a list of vertex labels
        - ``labels`` -- boolean

        OUTPUT:
            a generator which yields edges, with or without labels
            depending on the labels parameter.

        TESTS::

            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.iterator_in_edges([],True)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def iterator_out_edges(self, vertices, labels):
        """
        Iterate over the outbound edges incident to a sequence of vertices.

        INPUT:

        - ``vertices`` -- a list of vertex labels
        - ``labels`` -- boolean

        OUTPUT:

            a generator which yields edges, with or without labels
            depending on the labels parameter.

        TESTS::

            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.iterator_out_edges([],True)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def iterator_nbrs(self, v):
        """
        Iterate over the vertices adjacent to v.

        INPUT:

        - ``v`` -- vertex label

        OUTPUT:

            a generator which yields vertex labels

        TESTS::

            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.iterator_nbrs(0)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def iterator_in_nbrs(self, v):
        """
        Iterate over the vertices u such that the edge (u,v) is in self
        (that is, predecessors of v).

        INPUT:

        - ``v`` -- vertex label

        OUTPUT:

            a generator which yields vertex labels

        TESTS::

            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.iterator_in_nbrs(0)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def iterator_out_nbrs(self, v):
        """
        Iterate over the vertices u such that the edge (v,u) is in self
        (that is, successors of v).

        INPUT:

        - ``v`` -- vertex label

        OUTPUT:

            a generator which yields vertex labels

        TESTS::

            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.iterator_out_nbrs(0)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def iterator_verts(self, verts):
        """
        Iterate over the vertices v with labels in verts.

        INPUT:

        - ``vertex`` -- vertex labels

        OUTPUT:

            a generator which yields vertices

        TESTS::

            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.iterator_verts(0)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()

    def loops(self, new=None):
        """
        Get/set whether or not self allows loops.

        INPUT:

        - ``new`` -- can be a boolean (in which case it sets the value) or
          ``None``, in which case the current value is returned. It is set to
          ``None`` by default.

        TESTS::

            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.loops(True)
            Traceback (most recent call last):
            ...
            NotImplementedError
            sage: G.loops(None)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def multiple_edges(self, new=None):
        """
        Get/set whether or not self allows multiple edges.

        INPUT:

        - ``new`` -- can be a boolean (in which case it sets the value) or
          ``None``, in which case the current value is returned. It is set to
          ``None`` by default.

        TESTS::

            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.multiple_edges(True)
            Traceback (most recent call last):
            ...
            NotImplementedError
            sage: G.multiple_edges(None)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def name(self, new=None):
        """
        Get/set name of self.

        INPUT:

        - ``new`` -- can be a string (in which case it sets the value) or
          ``None``, in which case the current value is returned. It is set to
          ``None`` by default.

        TESTS::

            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.name("A Generic Graph")
            Traceback (most recent call last):
            ...
            NotImplementedError
            sage: G.name(None)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def num_edges(self, directed):
        """
        The number of edges in self

        INPUT:

        - ``directed`` -- boolean

        TESTS::

            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.num_edges(True)
            Traceback (most recent call last):
            ...
            NotImplementedError
            sage: G.num_edges(False)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def num_verts(self):
        """
        The number of vertices in self

        TESTS::

            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.num_verts()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def relabel(self, perm, directed):
        """
        Relabel the vertices of self by a permutation.

        INPUT:

        - ``perm`` -- permutation
        - ``directed`` -- boolean

        TESTS::

            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.relabel([],False)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()
    def set_edge_label(self, u, v, l, directed):
        """
        Label the edge (u,v) by l.

        INPUT:

        - ``u,v`` -- vertices
        - ``l`` -- edge label
        - ``directed`` -- boolean

        TESTS::

            sage: G = sage.graphs.base.graph_backends.GenericGraphBackend()
            sage: G.set_edge_label(1,2,'a',True)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError()

    def __reduce__(self):
        r""""
        Return a tuple used for pickling this graph.

        This function returns a pair ``(f, args)`` such that ``f(*args)``
        produces a copy of ``self``. The function returned is always
        :func:`unpickle_graph_backend`.

        Pickling of the static graph backend makes pickling of immutable
        graphs and digraphs work::

            sage: G = Graph(graphs.PetersenGraph(), immutable=True)
            sage: G == loads(dumps(G))
            True
            sage: uc = [[2,3], [], [1], [1], [1], [3,4]]
            sage: D = DiGraph(dict([[i,uc[i]] for i in range(len(uc))]), immutable=True)
            sage: loads(dumps(D)) == D
            True

        No problems with loops and multiple edges, with Labels::

            sage: g = Graph(multiedges = True, loops = True)
            sage: g.add_edges(2*graphs.PetersenGraph().edges())
            sage: g.add_edge(0,0)
            sage: g.add_edge(1,1, "a label")
            sage: g.add_edge([(0,1,"labellll"), (0,1,"labellll"), (0,1,"LABELLLL")])
            sage: g.add_vertex("isolated vertex")
            sage: gi = g.copy(immutable=True)
            sage: loads(dumps(gi)) == gi
            True

        Similar, with a directed graph::

            sage: g = DiGraph(multiedges = True, loops = True)
            sage: H = 2*(digraphs.Circuit(15)+DiGraph(graphs.PetersenGraph()))
            sage: g.add_edges(H.edges())
            sage: g.add_edge(0,0)
            sage: g.add_edge(1,1, "a label")
            sage: g.add_edge([(0,1,"labellll"), (0,1,"labellll"), (0,1,"LABELLLL")])
            sage: g.add_vertex("isolated vertex")
            sage: gi = g.copy(immutable=True)
            sage: loads(dumps(gi)) == gi
            True
        """
        from .static_sparse_backend import StaticSparseBackend
        from .sparse_graph import SparseGraphBackend
        from .dense_graph import DenseGraphBackend

        # implementation, data_structure, multiedges, directed, loops
        if isinstance(self, CGraphBackend):
            implementation = "c_graph"
            if isinstance(self,SparseGraphBackend):
                data_structure = "sparse"
            elif isinstance(self,DenseGraphBackend):
                data_structure = "dense"
            elif isinstance(self,StaticSparseBackend):
                implementaton = "static_sparse"
            else:
                raise Exception
            multiedges = (<CGraphBackend> self)._multiple_edges
            directed   = (<CGraphBackend> self)._directed
            loops      = (<CGraphBackend> self)._loops
        elif isinstance(self, NetworkXGraphBackend):
            data_structure = "implementation"
            implementation = "networkx"
            multiedges =  self._nxg.is_multigraph()
            directed   =  self._nxg.is_directed()
            loops      =  bool(self._nxg.number_of_selfloops)
        else:
            raise Exception

        # Vertices and edges
        vertices = list(self.iterator_verts(None))
        if directed:
            edges    = list(self.iterator_out_edges(vertices,True))
        else:
            edges    = list(self.iterator_edges(vertices,True))

        return (unpickle_graph_backend,
                (directed, vertices, edges,
                 {'loops': loops,
                  'multiedges': multiedges}))

def unpickle_graph_backend(directed,vertices,edges,kwds):
    r"""
    Return a backend from its pickled data

    This methods is defined because Python's pickling mechanism can only build
    objects from a pair ``(f,args)`` by running ``f(*args)``. In particular,
    there is apparently no way to define a ``**kwargs`` (i.e. define the value
    of keyword arguments of ``f``), which means that one must know the order of
    all arguments of ``f`` (here, ``f`` is :class:`Graph` or :class:`DiGraph`).

    As a consequence, this means that the order cannot change in the future,
    which is something we cannot swear.

    INPUT:

    - ``directed`` (boolean)

    - ``vertices`` -- list of vertices.

    - ``edges`` -- list of edges

    - ``kwds`` -- any dictionary whose keywords will be forwarded to the graph
      constructor.

    This function builds a :class:`Graph` or :class:`DiGraph` from its data, and
    returns the ``_backend`` attribute of this object.

    EXAMPLE::

        sage: from sage.graphs.base.graph_backends import unpickle_graph_backend
        sage: b = unpickle_graph_backend(0,[0,1,2,3],[(0,3,'label'),(0,0,1)],{'loops':True})
        sage: b
        <type 'sage.graphs.base.sparse_graph.SparseGraphBackend'>
        sage: list(b.iterator_edges(range(4),1))
        [(0, 0, 1), (0, 3, 'label')]
    """
    if directed:
        from sage.graphs.digraph import DiGraph as constructor
    else:
        from sage.graphs.graph import Graph as constructor

    G = constructor(data=edges,**kwds)
    G.add_vertices(vertices)
    return G._backend

class NetworkXGraphDeprecated(SageObject):
    """
    Class for unpickling old networkx.XGraph formats

    TESTS::

        sage: from sage.graphs.base.graph_backends import NetworkXGraphDeprecated as NXGD
        sage: X = NXGD()
        doctest:...

    """

    def __init__(self):
        """
        Issue deprecation warnings for the old networkx XGraph formats

        EXAMPLE::

            sage: from sage.graphs.base.graph_backends import NetworkXGraphDeprecated
            sage: NetworkXGraphDeprecated()
            <class 'sage.graphs.base.graph_backends.NetworkXGraphDeprecated'>
        """
        from sage.misc.superseded import deprecation
        deprecation(10900, "Your graph object is saved in an old format since networkx "+
                    "was updated to 1.0.1. You can re-save your graph by typing "+
                    "graph.save(filename) to make this warning go away.")

    def mutate(self):
        """
        Change the old networkx XGraph format into the new one.

        OUTPUT:

        - The networkx.Graph or networkx.MultiGraph corresponding to the
          unpickled data.

        EXAMPLES::

            sage: from sage.graphs.base.graph_backends import NetworkXGraphDeprecated as NXGD
            sage: X = NXGD()
            sage: X.adj = {1:{2:7}, 2:{1:7}, 3:{2:[4,5,6,7]}, 2:{3:[4,5,6,7]}}
            sage: X.multiedges = True
            sage: G = X.mutate()
            sage: G.edges()
            [(1, 2), (2, 3)]
            sage: G.edges(data=True)
            [(1, 2, {'weight': 7}), (2, 3, {4: {}, 5: {}, 6: {}, 7: {}})]

        """
        import networkx
        new_adj = {}

        for node1, edges in self.adj.iteritems():
            new_adj[node1] = {}
            for node2, weights in edges.iteritems():
                new_adj[node1][node2] = {}
                if weights is not None:
                    try:
                        for weight in weights:
                            new_adj[node1][node2][weight] = {}
                    except TypeError:
                        new_adj[node1][node2]['weight'] = weights

        if self.multiedges:
            G = networkx.MultiGraph(new_adj)
        else:
            G = networkx.Graph(new_adj)

        return G

class NetworkXDiGraphDeprecated(SageObject):
    """
    Class for unpickling old networkx.XDiGraph formats

    TESTS::

        sage: import sage.graphs.base.graph_backends
    """

    def __init__(self):
        """
        Issue deprecation warnings for the old networkx XDiGraph formats

        EXAMPLE::

            sage: from sage.graphs.base.graph_backends import NetworkXDiGraphDeprecated
            sage: NetworkXDiGraphDeprecated()
            doctest:...
            <class 'sage.graphs.base.graph_backends.NetworkXDiGraphDeprecated'>
        """
        from sage.misc.superseded import deprecation
        deprecation(10900, "Your digraph object is saved in an old format since networkx "+
                    "was updated to 1.0.1. You can re-save your digraph by typing "+
                    "digraph.save(filename) to make this warning go away.")

    def mutate(self):
        """
        Change the old networkx XDiGraph format into the new one.

        OUTPUT:

        - The networkx.DiGraph or networkx.MultiDiGraph corresponding to the
          unpickled data.

        EXAMPLES::

            sage: from sage.graphs.base.graph_backends import NetworkXDiGraphDeprecated as NXDGD
            sage: X = NXDGD()
            sage: X.adj = {1:{2:7}, 2:{1:[7,8], 3:[4,5,6,7]}}
            sage: X.multiedges = True
            sage: G = X.mutate()
            sage: G.edges()
            [(1, 2), (2, 1), (2, 3)]
            sage: G.edges(data=True)
            [(1, 2, {'weight': 7}),
             (2, 1, {7: {}, 8: {}}),
             (2, 3, {4: {}, 5: {}, 6: {}, 7: {}})]

        """
        import networkx
        new_adj = {}

        for node1, edges in self.adj.iteritems():
            new_adj[node1] = {}
            for node2, weights in edges.iteritems():
                new_adj[node1][node2] = {}
                if weights is not None:
                    try:
                        for weight in weights:
                            new_adj[node1][node2][weight] = {}
                    except TypeError:
                        new_adj[node1][node2]['weight'] = weights

        if self.multiedges:
            G = networkx.MultiDiGraph(new_adj)
        else:
            G = networkx.DiGraph(new_adj)

        return G

from sage.structure.sage_object import register_unpickle_override
register_unpickle_override('networkx.xgraph','XGraph', NetworkXGraphDeprecated)
register_unpickle_override('networkx.xdigraph','XDiGraph', NetworkXDiGraphDeprecated)

class NetworkXGraphBackend(GenericGraphBackend):
    """
    A wrapper for NetworkX as the backend of a graph.

    TESTS::

        sage: import sage.graphs.base.graph_backends

    """

    _nxg = None

    def __init__(self, N=None):
        """
        Initialize the backend with NetworkX graph N.

        TESTS::

            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            doctest:...: DeprecationWarning: This class is not supported anymore and will soon be removed
            See http://trac.sagemath.org/18375 for details.
            sage: G.iterator_edges([],True)
            <generator object at ...>
        """
        if N is None:
            import networkx
            N = networkx.MultiGraph()
        self._nxg = N
        from sage.misc.superseded import deprecation
        deprecation(18375,"This class is not supported anymore and will "
                    "soon be removed")

    def __setstate__(self,state):
        r"""
        Fix the deprecated class if necessary.

        EXAMPLE::

            sage: sage.structure.sage_object.unpickle_all() # indirect random
        """
        for k,v in state.iteritems():
            self.__dict__[k] = v
        if isinstance(self._nxg, (NetworkXGraphDeprecated, NetworkXDiGraphDeprecated)):
            from sage.misc.superseded import deprecation
            deprecation(18375,"You unpickled an object which relies on an old "
                        "data structure. Save it again to update it, for it "
                        "may break in the future.")
            self._nxg = self._nxg.mutate()

    def add_edge(self, u, v, l, directed):
        """
        Add an edge (u,v) to self, with label l.  If directed is True, this is
        interpreted as an arc from u to v.

        INPUT:

        - ``u,v`` -- vertices
        - ``l`` -- edge label
        - ``directed`` -- boolean

        TESTS::

            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.add_edge(1,2,'a',True)
        """

        if u is None: u = self.add_vertex(None)
        if v is None: v = self.add_vertex(None)

        if l:
            self._nxg.add_edge(u, v, weight=l)
        else:
            self._nxg.add_edge(u, v)

    def add_edges(self, edges, directed):
        """
        Add a sequence of edges to self.  If directed is True, these are
        interpreted as arcs.

        INPUT:

        - ``edges`` -- list/iterator of edges to be added.
        - ``directed`` -- boolean

        TESTS::

            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.add_edges([],True)
        """
        for e in edges:
            try:
                u,v,l = e
            except ValueError:
                u,v = e
                l = None
            self.add_edge(u,v,l,directed)

    def add_vertex(self, name):
        """
        Add a labelled vertex to self.

        INPUT:

        - ``name``: vertex label

        OUTPUT:

        If ``name=None``, the new vertex name is returned. ``None`` otherwise.

        TESTS::

            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.add_vertex(0)
        """
        if isinstance(self._nxg, (NetworkXGraphDeprecated, NetworkXDiGraphDeprecated)):
            self._nxg = self._nxg.mutate()

        retval = None
        if name is None: # then find an integer to use as a key
            i = 0
            while self.has_vertex(i):
                i=i+1
            name = i
            retval = name

        self._nxg.add_node(name)

        return retval

    def add_vertices(self, vertices):
        """
        Add labelled vertices to self.

        INPUT:

        - ``vertices``: iterator of vertex labels. A new label is created, used and returned in
          the output list for all ``None`` values in ``vertices``.

        OUTPUT:

        Generated names of new vertices if there is at least one ``None`` value
        present in ``vertices``. ``None`` otherwise.

        EXAMPLES::

            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.add_vertices([1,2,3])
            sage: G.add_vertices([4,None,None,5])
            [0, 6]
        """
        vertices = list(vertices)
        nones = vertices.count(None)
        vertices = [v for v in vertices if v is not None]
        self._nxg.add_nodes_from(vertices)

        new_names = []
        i = 0
        while nones > 0:
            while self.has_vertex(i):
                i += 1
            self._nxg.add_node(i)
            new_names.append(i)

            nones -= 1
            i += 1

        return new_names if new_names != [] else None

    def degree(self, v, directed):
        """
        Return the total number of vertices incident to `v`.

        INPUT:

        - ``v`` -- a vertex label
        - ``directed`` -- boolean

        OUTPUT:

            degree of v

        TESTS::

            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.add_vertices(range(3))
            sage: G.degree(1, False)
            0
        """
        return self._nxg.degree(v)

    def in_degree(self, v):
        """
        Return the in-degree of `v`.

        INPUT:

        - ``v`` -- a vertex label

        OUTPUT:

            degree of v

        TESTS::

            sage: G = DiGraph(digraphs.Path(5),implementation="networkx")
            doctest:...: DeprecationWarning: The 'implementation' keyword is deprecated, and the graphs has been stored as a 'c_graph'
            See http://trac.sagemath.org/18375 for details.
            sage: G = G._backend
            sage: G.in_degree(0)
            0
            sage: G.in_degree(4)
            1
        """
        return self._nxg.in_degree(v)

    def out_degree(self, v):
        """
        Return the out-degree of `v`.

        INPUT:

        - ``v`` -- a vertex label

        OUTPUT:

            degree of v

        TESTS::

            sage: G = DiGraph(digraphs.Path(5),implementation="networkx")
            doctest:...: DeprecationWarning: The 'implementation' keyword is deprecated, and the graphs has been stored as a 'c_graph'
            See http://trac.sagemath.org/18375 for details.
            sage: G = G._backend
            sage: G.out_degree(0)
            1
            sage: G.out_degree(4)
            0
        """
        return self._nxg.out_degree(v)

    def del_edge(self, u, v, l, directed):
        """
        Delete the edge `(u,v)` with label `l`.

        INPUT:

        - ``u,v`` -- vertices
        - ``l`` -- edge label
        - ``directed`` -- boolean

        TESTS::

            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.del_edge(1,2,'a',True)
        """
        import networkx
        try:
            if self._nxg.is_multigraph():
                for k,d in self._nxg.edge[u][v].iteritems():
                    if d.get('weight',None) == l:
                        self._nxg.remove_edge(u,v,k)
                        break
            else:
                if l is None or self._nxg.edge[u][v].get('weight',None) == l:
                    self._nxg.remove_edge(u,v)
        except (KeyError, networkx.NetworkXError):
            pass


    def del_vertex(self, v):
        """
        Delete a labelled vertex in self.

        INPUT:

        - ``v`` -- vertex label

        TESTS::

            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.del_vertex(0)
            Traceback (most recent call last):
            ...
            NetworkXError: The node 0 is not in the graph.
        """
        self._nxg.remove_node(v)

    def del_vertices(self, vertices):
        """
        Delete labelled vertices in self.

        INPUT:

        - ``vertices`` -- iterator of vertex labels

        TESTS::

            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.del_vertices([1,2,3])
            Traceback (most recent call last):
            ...
            NetworkXError: The node 1 is not in the graph.
        """
        for v in vertices:
            self._nxg.remove_node(v)

    def get_edge_label(self, u, v):
        """
        Return the edge label of `(u,v)`.

        INPUT:

        - ``u,v`` -- vertex labels

        OUTPUT:
            label of `(u,v)`

        TESTS::

            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.get_edge_label(1,2)
            Traceback (most recent call last):
            ...
            NetworkXError: Edge (1,2) requested via get_edge_label does not exist.
        """
        try:
            E = self._nxg.edge[u][v]
        except KeyError:
            from networkx import NetworkXError
            raise NetworkXError("Edge (%s,%s) requested via get_edge_label does not exist."%(u,v))

        if self._nxg.is_multigraph():
            return [ e.get('weight',None) for e in itervalues(E) ]
        else:
            return E.get('weight',None)

    def has_edge(self, u, v, l):
        """
        True if self has an edge (u,v) with label l.

        INPUT:

        - ``u,v`` -- vertex labels
        - ``l`` -- label

        OUTPUT:
            boolean

        TESTS::

            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.has_edge(1,2,'a')
            False
        """
        if not self._nxg.has_edge(u, v):
            return False
        if l is None:
            return True
        if self._nxg.is_multigraph():
            return any(e.get('weight', None) == l
                       for e in itervalues(self._nxg.adj[u][v]))
        else:
            return any(e == l for e in itervalues(self._nxg.adj[u][v]))

    def has_vertex(self, v):
        """
        True if self has a vertex with label v.

        INPUT:

        - ``v`` -- vertex label

        OUTPUT:
            boolean

        TESTS::

            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.has_vertex(0)
            False
        """
        return self._nxg.has_node(v)

    def iterator_edges(self, vertices, labels):
        """
        Iterate over the edges incident to a sequence of vertices. Edges are
        assumed to be undirected.

        INPUT:

        - ``vertices`` -- a list of vertex labels
        - ``labels`` -- boolean

        OUTPUT:
            a generator which yields edges, with or without labels
            depending on the labels parameter.

        TESTS::

            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.iterator_edges([],True)
            <generator object at ...>
        """
        if labels:
            for u,v,d in self._nxg.edges_iter(data=True):
                if u in vertices or v in vertices:
                    yield (u,v,d.get('weight',None))
        else:
            for u,v in self._nxg.edges_iter():
                if u in vertices or v in vertices:
                    yield (u,v)

    def _iterator_in_edges_with_labels(self, vertices):
        """
        Iterate over the incoming edges incident to a sequence of vertices.
        Special case, only for internal use.

        EXAMPLE::

            sage: g = DiGraph(graphs.PetersenGraph(), implementation="networkx")._backend
            doctest:...: DeprecationWarning: The 'implementation' keyword is deprecated, and the graphs has been stored as a 'c_graph'
            See http://trac.sagemath.org/18375 for details.
            sage: sorted(list(g.iterator_in_edges([0,1], True)))
            [(0, 1, None), (1, 0, None), (2, 1, None), (4, 0, None), (5, 0, None), (6, 1, None)]
        """
        for u,v,d in self._nxg.in_edges_iter(vertices,data=True):
            yield (u,v,d.get('weight',None))

    def iterator_in_edges(self, vertices, labels):
        """
        Iterate over the incoming edges incident to a sequence of vertices.

        INPUT:

        - ``vertices`` -- a list of vertex labels
        - ``labels`` -- boolean

        OUTPUT:
            a generator which yields edges, with or without labels
            depending on the labels parameter.

        TESTS::

            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: i = G.iterator_in_edges([],True)
        """
        if self._nxg.is_directed():
            if labels:
                return self._iterator_in_edges_with_labels(vertices)
            else:
                return self._nxg.in_edges_iter(vertices)
        else:
            return self.iterator_edges(vertices,labels)

    def _iterator_out_edges_with_labels(self, vertices):
        """
        Iterate over the outbound edges incident to a sequence of vertices.
        Special case, only for internal use.

        EXAMPLE::

            sage: g = DiGraph(graphs.PetersenGraph(), implementation="networkx")._backend
            doctest:...: DeprecationWarning: The 'implementation' keyword is deprecated, and the graphs has been stored as a 'c_graph'
            See http://trac.sagemath.org/18375 for details.
            sage: sorted(list(g.iterator_out_edges([0,1], True)))
            [(0, 1, None), (0, 4, None), (0, 5, None), (1, 0, None), (1, 2, None), (1, 6, None)]
        """
        for u,v,d in self._nxg.out_edges_iter(vertices,data=True):
            yield (u,v,d.get('weight',None))

    def iterator_out_edges(self, vertices, labels):
        """
        Iterate over the outbound edges incident to a sequence of vertices.

        INPUT:

        - ``vertices`` -- a list of vertex labels
        - ``labels`` -- boolean

        OUTPUT:
            a generator which yields edges, with or without labels
            depending on the labels parameter.

        TESTS::

            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: i = G.iterator_out_edges([],True)
        """
        if self._nxg.is_directed():
            if labels:
                return self._iterator_out_edges_with_labels(vertices)
            else:
                return self._nxg.out_edges_iter(vertices)
        else:
            return self.iterator_edges(vertices,labels)

    def iterator_nbrs(self, v):
        """
        Iterate over the vertices adjacent to v.

        INPUT:

        - ``v`` -- vertex label

        OUTPUT:
            a generator which yields vertex labels

        TESTS::

            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.add_vertex(0)
            sage: G.iterator_nbrs(0)
            <dictionary-keyiterator object at ...>
        """
        return self._nxg.neighbors_iter(v)

    def iterator_in_nbrs(self, v):
        """
        Iterate over the vertices u such that the edge (u,v) is in self
        (that is, predecessors of v).

        INPUT:

        - ``v`` -- vertex label

        OUTPUT:
            a generator which yields vertex labels

        TESTS::

            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.iterator_in_nbrs(0)
            Traceback (most recent call last):
            ...
            AttributeError: 'MultiGraph' object has no attribute 'predecessors_iter'
        """
        return self._nxg.predecessors_iter(v)

    def iterator_out_nbrs(self, v):
        """
        Iterate over the vertices u such that the edge (v,u) is in self
        (that is, successors of v).

        INPUT:

        - ``v`` -- vertex label

        OUTPUT:
            a generator which yields vertex labels

        TESTS::

            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.iterator_out_nbrs(0)
            Traceback (most recent call last):
            ...
            AttributeError: 'MultiGraph' object has no attribute 'successors_iter'
        """
        return self._nxg.successors_iter(v)

    def iterator_verts(self, verts):
        """
        Iterate over the vertices v with labels in verts.

        INPUT:

        - ``vertex`` -- vertex labels

        OUTPUT:
            a generator which yields vertices

        TESTS::

            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.iterator_verts(0)
            <generator object bunch_iter at ...>
        """
        return self._nxg.nbunch_iter(verts)

    def loops(self, new=None):
        """
        Get/set whether or not self allows loops.

        INPUT:

        - ``new`` -- can be a boolean (in which case it sets the value) or
          ``None``, in which case the current value is returned. It is set to
          ``None`` by default.

        TESTS::

            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.loops(True)
            sage: G.loops(None)
            True
        """
        if new is None:
            return self._loops
        if new:
            self._loops = True
        else:
            self._loops = False

    def multiple_edges(self, new=None):
        """
        Get/set whether or not self allows multiple edges.

        INPUT:

        - ``new`` -- can be a boolean (in which case it sets the value) or
          ``None``, in which case the current value is returned. It is set to
          ``None`` by default.

        TESTS::

            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.multiple_edges(True)
            sage: G.multiple_edges(None)
            True
        """
        from networkx import Graph,MultiGraph,DiGraph,MultiDiGraph
        if new is None:
            return self._nxg.is_multigraph()
        if new == self._nxg.is_multigraph():
            return
        if new:
            if self._nxg.is_directed():
                self._nxg = MultiDiGraph(self._nxg)
            else:
                self._nxg = MultiGraph(self._nxg)
        else:
            if self._nxg.is_directed():
                self._nxg = DiGraph(self._nxg)
            else:
                self._nxg = Graph(self._nxg)

    def name(self, new=None):
        """
        Get/set name of self.

        INPUT:

        - ``new`` -- can be a string (in which case it sets the value) or
          ``None``, in which case the current value is returned. It is set to
          ``None`` by default.

        TESTS::

            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.name("A NetworkX Graph")
            sage: G.name(None)
            'A NetworkX Graph'
        """
        if new is None:
            return self._nxg.name
        self._nxg.name = new

    def num_edges(self, directed):
        """
        The number of edges in self

        INPUT:

        - ``directed`` -- boolean

        TESTS::

            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.num_edges(True)
            0
            sage: G.num_edges(False)
            0
        """
        return self._nxg.size()

    def num_verts(self):
        """
        The number of vertices in self

        TESTS::

            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.num_verts()
            0
        """
        return self._nxg.order()

    def relabel(self, perm, directed):
        """
        Relabel the vertices of self by a permutation.

        INPUT:

        - ``perm`` -- permutation
        - ``directed`` -- boolean

        TESTS::

            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.relabel([],False)
        """
        from networkx import relabel_nodes
        name = self._nxg.name
        self._nxg = relabel_nodes(self._nxg,perm)
        self._nxg.name = name

    def set_edge_label(self, u, v, l, directed):
        """
        Label the edge (u,v) by l.

        INPUT:

        - ``u,v`` -- vertices
        - ``l`` -- edge label
        - ``directed`` -- boolean

        TESTS::

            sage: G = sage.graphs.base.graph_backends.NetworkXGraphBackend()
            sage: G.set_edge_label(1,2,'a',True)
        """
        if not self.has_edge(u, v, None):
            return
        if self.multiple_edges(None):
            self._nxg[u][v].clear()
            self._nxg[u][v][0] = dict(weight=l)
            if directed is False:
                self._nxg[v][u].clear()
                self._nxg[v][u][0] = dict(weight=l)
        else:
            self._nxg[u][v]['weight'] = l
            if directed is False:
                self._nxg[v][u]['weight'] = l
