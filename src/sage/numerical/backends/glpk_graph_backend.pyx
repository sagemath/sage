# sage.doctest: needs sage.graphs
"""
GLPK Backend for access to GLPK graph functions

AUTHORS:

- Christian Kuper (2012-11): Initial implementation

Methods index
-------------

**Graph creation and modification operations:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~GLPKGraphBackend.add_vertex`          | Add an isolated vertex to the graph.
    :meth:`~GLPKGraphBackend.add_vertices`        | Add vertices from an iterable container of vertices.
    :meth:`~GLPKGraphBackend.set_vertex_demand`   | Set the vertex parameters.
    :meth:`~GLPKGraphBackend.set_vertices_demand` | Set the parameters of selected vertices.
    :meth:`~GLPKGraphBackend.get_vertex`          | Return a specific vertex as a dictionary Object.
    :meth:`~GLPKGraphBackend.get_vertices`        | Return a dictionary of the dictionaries associated to each vertex.
    :meth:`~GLPKGraphBackend.vertices`            | Return a ``list`` of all vertices.
    :meth:`~GLPKGraphBackend.delete_vertex`       | Remove a vertex from the graph.
    :meth:`~GLPKGraphBackend.delete_vertices`     | Remove vertices from the graph.
    :meth:`~GLPKGraphBackend.add_edge`            | Add an edge between vertices ``u`` and ``v``.
    :meth:`~GLPKGraphBackend.add_edges`           | Add edges to the graph.
    :meth:`~GLPKGraphBackend.get_edge`            | Return an edge connecting two vertices.
    :meth:`~GLPKGraphBackend.edges`               | Return a ``list`` of all edges in the graph.
    :meth:`~GLPKGraphBackend.delete_edge`         | Delete an edge from the graph.
    :meth:`~GLPKGraphBackend.delete_edges`        | Delete edges from the graph.

**Graph writing operations:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~GLPKGraphBackend.write_graph`    | Write the graph to a plain text file.
    :meth:`~GLPKGraphBackend.write_ccdata`   | Write the graph to a text file in DIMACS format.
    :meth:`~GLPKGraphBackend.write_mincost`  | Write the mincost flow problem data to a text file in DIMACS format.
    :meth:`~GLPKGraphBackend.write_maxflow`  | Write the maximum flow problem data to a text file in DIMACS format.

**Network optimization operations:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~GLPKGraphBackend.mincost_okalg`  | Find solution to the mincost problem with the out-of-kilter algorithm.
    :meth:`~GLPKGraphBackend.maxflow_ffalg`  | Find solution to the maxflow problem with Ford-Fulkerson algorithm.
    :meth:`~GLPKGraphBackend.cpp`            | Solve the critical path problem of a project network.

Classes and methods
-------------------
"""

#*****************************************************************************
#       Copyright (C) 2012 Christian Kuper <christian.kuper@t-online.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cysignals.memory cimport check_allocarray, sig_free

from sage.cpython.string cimport str_to_bytes, char_to_str
from sage.cpython.string import FS_ENCODING
from sage.libs.glpk.constants cimport *
from sage.libs.glpk.graph cimport *
from sage.numerical.mip import MIPSolverException

cdef class GLPKGraphBackend():
    """
    GLPK Backend for access to GLPK graph functions.

    The constructor can either be called without arguments (which results in an
    empty graph) or with arguments to read graph data from a file.

    INPUT:

    - ``data`` -- a filename or a :class:`Graph` object

    - ``format`` -- when ``data`` is a filename, specifies the format of the
      data read from a file.  The ``format`` parameter is a string and can take
      values as described in the table below.

    **Format parameters:**

    .. list-table::
     :widths: 10 70

     * - ``plain``

       - Read data from a plain text file containing the following information:

           | nv na
           | i[1] j[1]
           | i[2] j[2]
           | . . .
           | i[na] j[na]

         where:

         * nv is the number of vertices (nodes);

         * na is the number of arcs;

         * i[k], k = 1, . . . , na, is the index of tail vertex of arc k;

         * j[k], k = 1, . . . , na, is the index of head vertex of arc k.


     * - ``dimacs``

       - Read data from a plain ASCII text file in DIMACS format.
         A description of the DIMACS format can be found at
         http://dimacs.rutgers.edu/Challenges/.

     * - ``mincost``

       - Reads the mincost flow problem data from a text file in DIMACS format

     * - ``maxflow``

       - Reads the maximum flow problem data from a text file in DIMACS format

    .. NOTE::

        When ``data`` is a :class:`Graph`, the following restrictions are
        applied.

            * vertices -- the value of the demand of each vertex (see
              :meth:`set_vertex_demand`) is obtained from the numerical
              value associated with the key "rhs" if it is a dictionary.

            * edges -- The edge values used in the algorithms are read from the
              edges labels (and left undefined if the edge labels are equal to
              ``None``). To be defined, the labels must be dictionary objects with
              keys "low", "cap" and "cost". See :meth:`get_edge` for details.

    EXAMPLES:

    The following example creates an empty graph::

        sage: from sage.numerical.backends.glpk_graph_backend import GLPKGraphBackend
        sage: gbe = GLPKGraphBackend()

    The following example creates an empty graph, adds some data, saves the data
    to a file and loads it::

        sage: from sage.numerical.backends.glpk_graph_backend import GLPKGraphBackend
        sage: gbe = GLPKGraphBackend()
        sage: gbe.add_vertices([None, None])
        ['0', '1']
        sage: a = gbe.add_edge('0', '1')
        sage: import tempfile
        sage: with tempfile.NamedTemporaryFile() as f:
        ....:     _ = gbe.write_graph(f.name)
        ....:     gbe1 = GLPKGraphBackend(f.name, "plain")
        Writing graph to ...
        4 lines were written
        Reading graph from ...
        Graph has 2 vertices and 1 edge
        3 lines were read

    The following example imports a Sage ``Graph`` and then uses it to solve a
    maxflow problem::

        sage: from sage.numerical.backends.glpk_graph_backend import GLPKGraphBackend
        sage: g = graphs.PappusGraph()
        sage: for ed in g.edges(sort=False):
        ....:     g.set_edge_label(ed[0], ed[1], {"cap":1})
        sage: gbe = GLPKGraphBackend(g)
        sage: gbe.maxflow_ffalg('1', '2')
        3.0
    """

    def __cinit__(self, data = None, format = "plain"):
        """
        Constructor.

        The constructor can either be called without arguments creating an empty
        graph or with arguments to read graph data from a file or a Sage
        :class:`Graph`. See documentation of :class:`GLPKGraphBackend` for
        details.

        EXAMPLES::

            sage: from sage.numerical.backends.glpk_graph_backend import GLPKGraphBackend
            sage: gbe = GLPKGraphBackend()
        """

        from sage.graphs.graph import Graph

        self.graph = <glp_graph*> glp_create_graph(sizeof(c_v_data),
                       sizeof(c_a_data))

        if self.graph is NULL:
            raise MemoryError("Error allocating memory.")

        self.s = 1
        self.t = 1

        if isinstance(data, str):
            fname = str_to_bytes(data, FS_ENCODING, 'surrogateescape')
            res = 0
            if format == "plain":
                res = glp_read_graph(self.graph, fname)
            elif format == "dimacs":
                res = glp_read_ccdata(self.graph, 0, fname)
            elif format == "mincost":
                res = glp_read_mincost(self.graph, 0, 0, sizeof(double),
                   sizeof(double) + sizeof(double), fname)
            elif format == "maxflow":
                res = glp_read_maxflow(self.graph, &self.s, &self.t,
                   sizeof(double), fname)
            if res != 0:
                raise IOError("Could not read graph from file %s" % (fname))

        elif isinstance(data, Graph):
            self.__add_vertices_sage(data)
            self.__add_edges_sage(data)
        else:
            ValueError("Input data is not supported")

    cpdef add_vertex(self, name=None):
        """
        Add an isolated vertex to the graph.

        If the vertex already exists, nothing is done.

        INPUT:

        - ``name`` -- string of max 255 chars length. If no name is
          specified, then the vertex will be represented by the string
          representation of the ID of the vertex or - if this already exists -
          a string representation of the least integer not already representing
          a vertex.

        OUTPUT:

        If no ``name`` is passed as an argument, the new vertex name is
        returned. ``None`` otherwise.

        EXAMPLES::

            sage: from sage.numerical.backends.glpk_graph_backend import GLPKGraphBackend
            sage: gbe = GLPKGraphBackend()
            sage: gbe.add_vertex()
            '0'
            sage: gbe.add_vertex("2")
            sage: gbe.add_vertex()
            '1'
        """
        cdef int n
        cdef vn_t = 0

        if name is not None and self._find_vertex(name) >= 0:
            return None

        cdef int vn = glp_add_vertices(self.graph, 1)

        if name is not None:
            glp_set_vertex_name(self.graph, vn, str_to_bytes(name))
            return None

        else:
            s = str(vn - 1)
            n = self._find_vertex(s)

            # This is costly, but hopefully will not happen often.
            while n >= 0:
                vn_t += 1
                s = str(vn_t - 1)
                n = self._find_vertex(s)

            glp_set_vertex_name(self.graph, vn, str_to_bytes(s))
            return s

    cpdef __add_vertices_sage(self, g):
        """
        Add vertices to the GLPK Graph.

        This function is only used when importing a
        :class:`~sage.graphs.generic_graph.GenericGraph` object.

        EXAMPLES::

            sage: from sage.numerical.backends.glpk_graph_backend import GLPKGraphBackend
            sage: g = graphs.PappusGraph()
            sage: for ed in g.edges(sort=False):
            ....:     g.set_edge_label(ed[0], ed[1], {"cap":1})
            sage: gbe = GLPKGraphBackend(g)
            sage: gbe.maxflow_ffalg('1', '2')
            3.0
        """
        cdef int n
        cdef int i
        cdef double rhs
        cdef glp_vertex* vert

        verts = g.vertices(sort=True)
        n = len(verts)
        if n < 1:
            raise ValueError("Graph must contain vertices")

        glp_add_vertices(self.graph, n)

        for i in range(n):
            vert = self.graph.v[i+1]
            s = str(verts[i])
            glp_set_vertex_name(self.graph, i + 1, str_to_bytes(s))

            if g.get_vertex(verts[i]) is not None:
                try:
                    (<c_v_data *>vert.data).rhs = g.get_vertex(verts[i])["rhs"]
                except AttributeError:
                    pass

        glp_create_v_index(self.graph)

    cpdef list add_vertices(self, vertices):
        """
        Add vertices from an iterable container of vertices.

        Vertices that already exist in the graph will not be added again.

        INPUT:

        - ``vertices`` -- iterator of vertex labels (string); a label can be
          ``None``

        OUTPUT:

        Generated names of new vertices if there is at least one ``None`` value
        present in ``vertices``. ``None`` otherwise.

        EXAMPLES::

            sage: from sage.numerical.backends.glpk_graph_backend import GLPKGraphBackend
            sage: gbe = GLPKGraphBackend()
            sage: vertices = [None for i in range(3)]
            sage: gbe.add_vertices(vertices)
            ['0', '1', '2']
            sage: gbe.add_vertices(['A', 'B', None])
            ['5']
            sage: gbe.add_vertices(['A', 'B', 'C'])
            sage: gbe.vertices()
            ['0', '1', '2', 'A', 'B', '5', 'C']

        TESTS::

            sage: from sage.numerical.backends.glpk_graph_backend import GLPKGraphBackend
            sage: gbe = GLPKGraphBackend()
            sage: gbe.add_vertices([None, None, None, '1'])
            ['0', '2', '3']
        """

        # We do not want to have [None,None,None,1] as input as a vertex named
        # "1" would be created twice (a first time when adding a 'None' vertex,
        # and a second time when reading the last item of the list).
        nonecount = 0
        for v in vertices:
            if v is None:
                nonecount += 1
            else:
                self.add_vertex(v)

        if nonecount:
            return [self.add_vertex() for i in range(nonecount)]
        else:
            return None

    cpdef set_vertex_demand(self, vertex, demand):
        """
        Set the demand of the vertex in a mincost flow algorithm.

        INPUT:

        - ``vertex`` -- name of the vertex

        - ``demand`` -- the numerical value representing demand of the vertex in
          a mincost flow algorithm (it could be for instance `-1` to represent a
          sink, or `1` to represent a source and `0` for a neutral vertex). This
          can either be an ``int`` or ``float`` value.

        EXAMPLES::

            sage: from sage.numerical.backends.glpk_graph_backend import GLPKGraphBackend
            sage: gbe = GLPKGraphBackend()
            sage: vertices = [None for i in range(3)]
            sage: gbe.add_vertices(vertices)
            ['0', '1', '2']
            sage: gbe.set_vertex_demand('0', 2)
            sage: gbe.get_vertex('0')['rhs']
            2.0
            sage: gbe.set_vertex_demand('3', 2)
            Traceback (most recent call last):
            ...
            KeyError: 'Vertex 3 does not exist.'
        """
        cdef int n = self._find_vertex(vertex)

        if n < 0:
            raise KeyError("Vertex " + vertex + " does not exist.")

        cdef glp_vertex* vert = self.graph.v[n+1]
        cdef double val = demand
        (<c_v_data *>vert.data).rhs = val

    cpdef set_vertices_demand(self, list pairs):
        """
        Set the parameters of selected vertices.

        INPUT:

        - ``pairs`` -- list of pairs ``(vertex, demand)`` associating a demand
          to each vertex. For more information, see the documentation of
          :meth:`set_vertex_demand`.

        EXAMPLES::

            sage: from sage.numerical.backends.glpk_graph_backend import GLPKGraphBackend
            sage: gbe = GLPKGraphBackend()
            sage: vertices = [None for i in range(3)]
            sage: gbe.add_vertices(vertices)
            ['0', '1', '2']
            sage: gbe.set_vertices_demand([('0', 2), ('1', 3), ('3', 4)])
            sage: sorted(gbe.get_vertex('1').items())
            [('cut', 0), ('es', 0.0), ('ls', 0.0), ('pi', 0.0), ('rhs', 3.0)]
        """

        for v, param in pairs:
            try:
                self.set_vertex_demand(v, param)
            except KeyError:
                pass

    cpdef dict get_vertex(self, vertex):
        """
        Return a specific vertex as a dictionary Object.

        INPUT:

        - ``vertex`` -- the vertex label as string

        OUTPUT:

        The vertex as a dictionary object or ``None`` if the vertex does not
        exist. The dictionary contains the values used or created by the different
        algorithms. The values associated with the keys following keys contain:

            * "rhs"  -- The supply / demand value the vertex (mincost alg)
            * "pi"   -- The node potential (mincost alg)
            * "cut"  -- The cut flag of the vertex (maxflow alg)
            * "es"   -- The earliest start of task (cpp alg)
            * "ls"   -- The latest start of task (cpp alg)

        EXAMPLES::

            sage: from sage.numerical.backends.glpk_graph_backend import GLPKGraphBackend
            sage: gbe = GLPKGraphBackend()
            sage: verts = ["A", "B", "C", "D"]
            sage: gbe.add_vertices(verts)
            sage: sorted(gbe.get_vertex("A").items())
            [('cut', 0), ('es', 0.0), ('ls', 0.0), ('pi', 0.0), ('rhs', 0.0)]
            sage: gbe.get_vertex("F") is None
            True
        """

        cdef int i = self._find_vertex(vertex)
        if i < 0:
            return None

        cdef glp_vertex* vert = self.graph.v[i+1]
        cdef c_v_data * vdata = <c_v_data *> vert.data

        return {
            "rhs" : vdata.rhs,
            "pi"  : vdata.pi,
            "cut" : vdata.cut,
            "es"  : vdata.es,
            "ls"  : vdata.ls
            }

    cpdef dict get_vertices(self, verts):
        """
        Return a dictionary of the dictionaries associated to each vertex.

        INPUT:

        - ``verts`` -- iterable container of vertices

        OUTPUT:

        A list of pairs ``(vertex, properties)`` where ``properties`` is a
        dictionary containing the numerical values associated with a vertex. For
        more information, see the documentation of
        :meth:`GLPKGraphBackend.get_vertex`.

        EXAMPLES::

            sage: from sage.numerical.backends.glpk_graph_backend import GLPKGraphBackend
            sage: gbe = GLPKGraphBackend()
            sage: verts = ['A', 'B']
            sage: gbe.add_vertices(verts)
            sage: sorted(gbe.get_vertices(verts)['B'].items())
            [('cut', 0), ('es', 0.0), ('ls', 0.0), ('pi', 0.0), ('rhs', 0.0)]
            sage: gbe.get_vertices(["C", "D"])
            {}
        """
        vl = [(v, self.get_vertex(v)) for v in verts]
        return dict([(v, p) for v, p in vl if p is not None])

    cpdef list vertices(self):
        """
        Return the list of all vertices.

        .. NOTE::

           Changing elements of the ``list`` will not change anything in the
           the graph.

        .. NOTE::

           If a vertex in the graph does not have a name / label it will appear
           as ``None`` in the resulting ``list``.

        EXAMPLES::

            sage: from sage.numerical.backends.glpk_graph_backend import GLPKGraphBackend
            sage: gbe = GLPKGraphBackend()
            sage: verts = ["A", "B", "C"]
            sage: gbe.add_vertices(verts)
            sage: a = gbe.vertices(); a
            ['A', 'B', 'C']
            sage: a.pop(0)
            'A'
            sage: gbe.vertices()
            ['A', 'B', 'C']
        """

        return [char_to_str(self.graph.v[i+1].name)
                if self.graph.v[i+1].name is not NULL else None
                for i in range(self.graph.nv)]

    cpdef add_edge(self, u, v, dict params=None):
        """
        Add an edge between vertices ``u`` and ``v``.

        Allows adding an edge and optionally providing parameters used by the
        algorithms. If a vertex does not exist it is created.

        INPUT:

        - ``u`` -- the name (as string) of the tail vertex

        - ``v`` -- the name (as string) of the head vertex

        - ``params`` -- an optional dictionary containing the edge parameters used
          for the algorithms. The following keys are used:

            * ``low`` -- the minimum flow through the edge

            * ``cap`` -- the maximum capacity of the edge

            * ``cost`` -- the cost of transporting one unit through the edge

        EXAMPLES::

            sage: from sage.numerical.backends.glpk_graph_backend import GLPKGraphBackend
            sage: gbe = GLPKGraphBackend()
            sage: gbe.add_edge("A", "B", {"low":0.0, "cap":10.0, "cost":5})
            sage: gbe.vertices()
            ['A', 'B']
            sage: for ed in gbe.edges():
            ....:     print((ed[0], ed[1], ed[2]['cap'], ed[2]['cost'], ed[2]['low']))
            ('A', 'B', 10.0, 5.0, 0.0)
            sage: gbe.add_edge("B", "C", {"low":0.0, "cap":10.0, "cost":'5'})
            Traceback (most recent call last):
            ...
            TypeError: Invalid edge parameter.
        """
        cdef int i = self._find_vertex(u)
        cdef int j = self._find_vertex(v)

        if i < 0:
            self.add_vertex(u)
            i = self._find_vertex(u)

        if j < 0:
            self.add_vertex(v)
            j = self._find_vertex(v)

        cdef glp_arc *a

        a = glp_add_arc(self.graph, i+1, j+1)

        if params is not None:
            try:
                if "low" in params:
                    (<c_a_data *>a.data).low = params["low"]
                if "cap" in params:
                    (<c_a_data *>a.data).cap = params["cap"]
                if "cost" in params:
                    (<c_a_data *>a.data).cost = params["cost"]
            except TypeError:
                glp_del_arc(self.graph, a)
                raise TypeError("Invalid edge parameter.")

    cpdef list add_edges(self, edges):
        """
        Add edges to the graph.

        INPUT:

        - ``edges`` -- an iterable container of pairs of the form ``(u, v)``,
          where ``u`` is name (as string) of the tail vertex and ``v`` is the
          name (as string) of the head vertex or an iterable container of
          triples of the form ``(u, v, params)`` where params is a dictionary as
          described in ``add_edge``.

        EXAMPLES::

            sage: from sage.numerical.backends.glpk_graph_backend import GLPKGraphBackend
            sage: gbe = GLPKGraphBackend()
            sage: edges = [("A", "B", {"low":0.0, "cap":10.0, "cost":5})]
            sage: edges.append(("B", "C"))
            sage: gbe.add_edges(edges)
            sage: for ed in gbe.edges():
            ....:     print((ed[0], ed[1], ed[2]['cap'], ed[2]['cost'], ed[2]['low']))
            ('A', 'B', 10.0, 5.0, 0.0)
            ('B', 'C', 0.0, 0.0, 0.0)
            sage: edges = [("C", "D", {"low":0.0, "cap":10.0, "cost":5})]
            sage: edges.append(("C", "E", 5))
            sage: gbe.add_edges(edges)
            Traceback (most recent call last):
            ...
            TypeError: Argument 'params' has incorrect type ...
            sage: for ed in gbe.edges():
            ....:     print((ed[0], ed[1], ed[2]['cap'], ed[2]['cost'], ed[2]['low']))
            ('A', 'B', 10.0, 5.0, 0.0)
            ('B', 'C', 0.0, 0.0, 0.0)
            ('C', 'D', 10.0, 5.0, 0.0)
        """
        for ed in edges:
            self.add_edge(*ed)

    cpdef __add_edges_sage(self, g):
        """
        Add edges to the Graph.

        This function is only used when importing a ``GenericGraph``.

        EXAMPLES::

            sage: from sage.numerical.backends.glpk_graph_backend import GLPKGraphBackend
            sage: g = graphs.PappusGraph()
            sage: for ed in g.edges(sort=False):
            ....:     g.set_edge_label(ed[0], ed[1], {"cap":1})
            sage: gbe = GLPKGraphBackend(g)
            sage: gbe.maxflow_ffalg('1', '2')
            3.0
        """
        cdef glp_arc* a
        cdef int u
        cdef int v
        cdef double cost
        cdef double cap
        cdef double low
        cdef int isdirected = g.is_directed()

        for eu, ev, label in g.edges(sort=False):
            u_name = str(eu)
            v_name = str(ev)
            u = glp_find_vertex(self.graph, str_to_bytes(u_name))
            v = glp_find_vertex(self.graph, str_to_bytes(v_name))
            if u < 1 or v < 1:
                raise IndexError(u_name + " or " + v_name + " not found")

            a = glp_add_arc(self.graph, u, v)

            if isinstance(label, dict):
                if "cost" in label:
                    cost = label["cost"]
                    (<c_a_data *>a.data).cost = cost
                if "cap" in label:
                    cap = label["cap"]
                    (<c_a_data *>a.data).cap = cap
                if "low" in label:
                    low = label["low"]
                    (<c_a_data *>a.data).low = low

            if not isdirected:
                a = glp_add_arc(self.graph, v, u)
                if isinstance(label, dict):
                    if "cost" in label:
                        (<c_a_data *>a.data).cost = cost
                    if "cap" in label:
                        (<c_a_data *>a.data).cap = cap
                    if "low" in label:
                        (<c_a_data *>a.data).low = low

    cpdef tuple get_edge(self, u, v):
        """
        Return an edge connecting two vertices.

        .. NOTE::

           If multiple edges connect the two vertices only the first edge
           found is returned.

        INPUT:

        - ``u`` -- name (as string) of the tail vertex
        - ``v`` -- name (as string) of the head vertex

        OUTPUT:

        A ``triple`` describing if edge was found or ``None`` if not. The third
        value of the triple is a dictionary containing the following edge
        parameters:

            * ``low`` -- the minimum flow through the edge
            * ``cap`` -- the maximum capacity of the edge
            * ``cost`` -- the cost of transporting one unit through the edge
            * ``x`` -- the actual flow through the edge after solving

        EXAMPLES::

            sage: from sage.numerical.backends.glpk_graph_backend import GLPKGraphBackend
            sage: gbe = GLPKGraphBackend()
            sage: edges = [("A", "B"), ("A", "C"), ("B", "C")]
            sage: gbe.add_edges(edges)
            sage: ed = gbe.get_edge("A", "B")
            sage: ed[0], ed[1], ed[2]['x']
            ('A', 'B', 0.0)
            sage: gbe.get_edge("A", "F") is None
            True
        """
        cdef int i = self._find_vertex(u)
        cdef int j = self._find_vertex(v)

        if i < 0 or j < 0:
            return None

        cdef glp_vertex* vert_u = self.graph.v[i+1]
        cdef glp_vertex* vert_v = self.graph.v[j+1]
        cdef glp_arc* a = vert_u.out
        while a is not NULL:
            if a.head == vert_v:
                return (u, v, {"low":(<c_a_data *>a.data).low,
                               "cap":(<c_a_data *>a.data).cap,
                               "cost":(<c_a_data *>a.data).cost,
                               "x":(<c_a_data *>a.data).x})
            a = a.t_next

        return None

    cpdef list edges(self):
        """
        Return a ``list`` of all edges in the graph.

        OUTPUT: a ``list`` of ``triples`` representing the edges of the graph

        EXAMPLES::

            sage: from sage.numerical.backends.glpk_graph_backend import GLPKGraphBackend
            sage: gbe = GLPKGraphBackend()
            sage: edges = [("A", "B", {"low":0.0, "cap":10.0, "cost":5})]
            sage: edges.append(("B", "C"))
            sage: gbe.add_edges(edges)
            sage: for ed in gbe.edges():
            ....:     print((ed[0], ed[1], ed[2]['cost']))
            ('A', 'B', 5.0)
            ('B', 'C', 0.0)
        """

        cdef int i = 1
        cdef glp_vertex* vert_u
        cdef glp_vertex* vert_v
        cdef glp_arc* a
        edge_list = []

        while i <= self.graph.nv:
            vert_u = self.graph.v[i]
            a = vert_u.out
            while a is not NULL:
                vert_v = a.head
                if vert_u.name is NULL:
                    u_name = None
                else:
                    u_name = char_to_str(vert_u.name)
                if vert_v.name is NULL:
                    v_name = None
                else:
                    v_name = char_to_str(vert_v.name)
                edge_list.append((u_name, v_name,
                              {"low":(<c_a_data *>a.data).low,
                               "cap":(<c_a_data *>a.data).cap,
                               "cost":(<c_a_data *>a.data).cost,
                               "x":(<c_a_data *>a.data).x}))
                a = a.t_next
            i += 1
        return edge_list

    cpdef delete_vertex(self, vert):
        r"""
        Removes a vertex from the graph.

        Trying to delete a non existing vertex will raise an exception.

        INPUT:

        - ``vert`` -- the name (as string) of the vertex to delete

        EXAMPLES::

            sage: from sage.numerical.backends.glpk_graph_backend import GLPKGraphBackend
            sage: gbe = GLPKGraphBackend()
            sage: verts = ["A", "D"]
            sage: gbe.add_vertices(verts)
            sage: gbe.delete_vertex("A")
            sage: gbe.vertices()
            ['D']
            sage: gbe.delete_vertex("A")
            Traceback (most recent call last):
            ...
            RuntimeError: Vertex A does not exist.
        """

        cdef int i = self._find_vertex(vert)

        if i < 0:
            raise RuntimeError("Vertex %s does not exist."%(vert))

        cdef int num[2]
        num[1] = i + 1
        cdef int ndel = 1

        glp_del_vertices(self.graph, ndel, num)

    cpdef delete_vertices(self, list verts):
        r"""
        Removes vertices from the graph.

        Trying to delete a non existing vertex will raise an exception.

        INPUT:

        - ``verts`` -- iterable container containing names (as string) of the
          vertices to delete

        EXAMPLES::

            sage: from sage.numerical.backends.glpk_graph_backend import GLPKGraphBackend
            sage: gbe = GLPKGraphBackend()
            sage: verts = ["A", "B", "C", "D"]
            sage: gbe.add_vertices(verts)
            sage: v_d = ["A", "B"]
            sage: gbe.delete_vertices(v_d)
            sage: gbe.vertices()
            ['C', 'D']
            sage: gbe.delete_vertices(["C", "A"])
            Traceback (most recent call last):
            ...
            RuntimeError: Vertex A does not exist.
            sage: gbe.vertices()
            ['C', 'D']
        """

        verts_val = [self._find_vertex(v) for v in verts]
        if -1 in verts_val:
            i = verts_val.index(-1)
            raise RuntimeError("Vertex %s does not exist."%(verts[i]))

        cdef int * num = <int *>check_allocarray(len(verts_val) + 1, sizeof(int))
        cdef int ndel = len(verts_val)

        for i,(v) in enumerate(verts_val):
            num[i+1] = v+1

        glp_del_vertices(self.graph, ndel, num)

        sig_free(num)

    cpdef delete_edge(self, u, v, dict params=None):
        """
        Delete an edge from the graph.

        If an edge does not exist it is ignored.

        INPUT:

        - ``u`` -- the name (as string) of the tail vertex of the edge
        - ``v`` -- the name (as string) of the tail vertex of the edge
        - ``params`` -- an optional dictionary containing the edge
          parameters (see :meth:`add_edge`). If this parameter
          is not provided, all edges connecting ``u`` and ``v`` are deleted.
          Otherwise only edges with matching parameters are deleted.

        .. SEEALSO::

            :meth:`delete_edges`

        EXAMPLES::

            sage: from sage.numerical.backends.glpk_graph_backend import GLPKGraphBackend
            sage: gbe = GLPKGraphBackend()
            sage: edges = [("A", "B", {"low":0.0, "cap":10.0, "cost":5})]
            sage: edges.append(("A", "B", {"low":0.0, "cap":15.0, "cost":10}))
            sage: edges.append(("B", "C", {"low":0.0, "cap":20.0, "cost":1}))
            sage: edges.append(("B", "C", {"low":0.0, "cap":10.0, "cost":20}))
            sage: gbe.add_edges(edges)
            sage: gbe.delete_edge("A", "B")
            sage: gbe.delete_edge("B", "C", {"low":0.0, "cap":10.0, "cost":20})
            sage: gbe.edges()[0][0], gbe.edges()[0][1], gbe.edges()[0][2]['cost']
            ('B', 'C', 1.0)
        """

        cdef int i = self._find_vertex(u)
        cdef int j = self._find_vertex(v)
        if i < 0 or j < 0:
            return

        cdef glp_vertex* vert_u = self.graph.v[i+1]
        cdef glp_vertex* vert_v = self.graph.v[j+1]
        cdef glp_arc* a = vert_u.out
        cdef glp_arc* a2 = a

        cdef double low, cap, cost, x

        if params is not None:
            if "low" in params:
                low = params["low"]
            if "cap" in params:
                cap = params["cap"]
            if "cost" in params:
                cost = params["cost"]
            if "x" in params:
                x = params["x"]

        while a is not NULL:
            a2 = a.t_next
            if a.head == vert_v and params is None:
                glp_del_arc(self.graph, a)
            elif a.head == vert_v:
                del_it = True
                if "low" in params:
                    if (<c_a_data *>a.data).low != low:
                        del_it = False
                if "cap" in params:
                    if (<c_a_data *>a.data).cap != cap:
                        del_it = False
                if "cost" in params:
                    if (<c_a_data *>a.data).cost != cost:
                        del_it = False
                if "x" in params:
                    if (<c_a_data *>a.data).x != x:
                        del_it = False
                if del_it:
                    glp_del_arc(self.graph, a)

            a = a2

    def delete_edges(self, edges):
        """
        Delete edges from the graph.

        Non existing edges are ignored.

        INPUT:

        - ``edges`` -- an iterable container of edges

        .. SEEALSO::

            :meth:`delete_edge`

        EXAMPLES::

            sage: from sage.numerical.backends.glpk_graph_backend import GLPKGraphBackend
            sage: gbe = GLPKGraphBackend()
            sage: edges = [("A", "B", {"low":0.0, "cap":10.0, "cost":5})]
            sage: edges.append(("A", "B", {"low":0.0, "cap":15.0, "cost":10}))
            sage: edges.append(("B", "C", {"low":0.0, "cap":20.0, "cost":1}))
            sage: edges.append(("B", "C", {"low":0.0, "cap":10.0, "cost":20}))
            sage: gbe.add_edges(edges)
            sage: gbe.delete_edges(edges[1:])
            sage: len(gbe.edges())
            1
            sage: gbe.edges()[0][0], gbe.edges()[0][1], gbe.edges()[0][2]['cap']
            ('A', 'B', 10.0)
        """

        for edge in edges:
            self.delete_edge(*edge)

    cpdef int _find_vertex(self, name) noexcept:
        """
        Return the index of a vertex specified by a name

        INPUT:

        - ``name`` -- name of the vertex

        OUTPUT: the index of the vertex or ``-1`` if the vertex is not found

        EXAMPLES::

            sage: from sage.numerical.backends.glpk_graph_backend import GLPKGraphBackend
            sage: gbe = GLPKGraphBackend()
            sage: verts = ["A", "B", "C", "D"]
            sage: gbe.add_vertices(verts)
            sage: gbe._find_vertex("A")
            0
            sage: gbe._find_vertex("F")
            -1
        """

        glp_create_v_index(self.graph)
        return glp_find_vertex(self.graph, str_to_bytes(name)) - 1

    cpdef int write_graph(self, fname) noexcept:
        r"""
        Writes the graph to a plain text file

        INPUT:

        - ``fname`` -- full name of the file

        OUTPUT: zero if the operations was successful otherwise nonzero

        EXAMPLES::

            sage: from sage.numerical.backends.glpk_graph_backend import GLPKGraphBackend
            sage: gbe = GLPKGraphBackend()
            sage: a = gbe.add_edge("0", "1")
            sage: import tempfile
            sage: with tempfile.NamedTemporaryFile() as f:
            ....:     gbe.write_graph(f.name)
            Writing graph to ...
            4 lines were written
            0
        """

        fname = str_to_bytes(fname, FS_ENCODING, 'surrogateescape')
        return glp_write_graph(self.graph, fname)

    cpdef int write_ccdata(self, fname) noexcept:
        r"""
        Writes the graph to a text file in DIMACS format.

        Writes the data to plain ASCII text file in DIMACS format.
        A description of the DIMACS format can be found at
        http://dimacs.rutgers.edu/Challenges/.

        INPUT:

        - ``fname`` -- full name of the file

        OUTPUT: zero if the operations was successful otherwise nonzero

        EXAMPLES::

            sage: from sage.numerical.backends.glpk_graph_backend import GLPKGraphBackend
            sage: gbe = GLPKGraphBackend()
            sage: a = gbe.add_edge("0", "1")
            sage: import tempfile
            sage: with tempfile.NamedTemporaryFile() as f:
            ....:     gbe.write_ccdata(f.name)
            Writing graph to ...
            6 lines were written
            0
        """

        fname = str_to_bytes(fname, FS_ENCODING, 'surrogateescape')
        return glp_write_ccdata(self.graph, 0, fname)

    cpdef int write_mincost(self, fname) noexcept:
        """
        Writes the mincost flow problem data to a text file in DIMACS format

        INPUT:

        - ``fname`` -- full name of file

        OUTPUT: zero if successful, otherwise nonzero

        EXAMPLES::

            sage: from sage.numerical.backends.glpk_graph_backend import GLPKGraphBackend
            sage: gbe = GLPKGraphBackend()
            sage: a = gbe.add_edge("0", "1")
            sage: import tempfile
            sage: with tempfile.NamedTemporaryFile() as f:
            ....:     gbe.write_mincost(f.name)
            Writing min-cost flow problem data to ...
            4 lines were written
            0
        """

        fname = str_to_bytes(fname, FS_ENCODING, 'surrogateescape')
        return glp_write_mincost(self.graph, 0, 0, sizeof(double),
                   sizeof(double) + sizeof(double), fname)

    cpdef double mincost_okalg(self) except -1:
        r"""
        Finds solution to the mincost problem with the out-of-kilter algorithm.

        The out-of-kilter algorithm requires all problem data to be integer
        valued.

        OUTPUT:

        The solution to the mincost problem, i.e. the total cost, if operation
        was successful.

        .. NOTE::

           This method raises ``MIPSolverException`` exceptions when
           the solution cannot be computed for any reason (none
           exists, or the LP solver was not able to find it, etc...)

        EXAMPLES::

            sage: from sage.numerical.backends.glpk_graph_backend import GLPKGraphBackend
            sage: gbe = GLPKGraphBackend()
            sage: vertices = (35, 50, 40, -45, -20, -30, -30)
            sage: vs = gbe.add_vertices([None for i in range(len(vertices))])
            sage: v_dict = {}
            sage: for i, v in enumerate(vs):
            ....:     v_dict[v] = vertices[i]
            sage: gbe.set_vertices_demand(list(v_dict.items()))
            sage: cost = ((8, 6, 10, 9), (9, 12, 13, 7), (14, 9, 16, 5))

            sage: for i in range(len(cost)):
            ....:     for j in range(len(cost[0])):
            ....:         gbe.add_edge(str(i), str(j + len(cost)), {"cost":cost[i][j], "cap":100})
            sage: gbe.mincost_okalg()
            1020.0
            sage: for ed in gbe.edges():
            ....:     print("{} -> {} {}".format(ed[0], ed[1], ed[2]["x"]))
            0 -> 6 0.0
            0 -> 5 25.0
            0 -> 4 10.0
            0 -> 3 0.0
            1 -> 6 0.0
            1 -> 5 5.0
            1 -> 4 0.0
            1 -> 3 45.0
            2 -> 6 30.0
            2 -> 5 0.0
            2 -> 4 10.0
            2 -> 3 0.0
        """
        cdef double graph_sol
        cdef int status = glp_mincost_okalg(self.graph, 0, 0, sizeof(double),
              2 * sizeof(double), &graph_sol,
              3 * sizeof(double), sizeof(double))
        if status == 0:
            pass
        elif status == GLP_ENOPFS:
            raise MIPSolverException("No (primal) feasible solution exists")
        elif status == GLP_EDATA:
            raise MIPSolverException("Unable to start search due to " +
                                     "problem data")
        elif status == GLP_ERANGE:
            raise MIPSolverException("The search was prematurely terminated " +
                                     "because of integer overflow")
        elif status == GLP_EFAIL:
            raise MIPSolverException("An error has been detected" +
                                     "in the program logic")

        return graph_sol

    cpdef int write_maxflow(self, fname) except -1:
        """
        Writes the maximum flow problem data to a text file in DIMACS format.

        INPUT:

        - ``fname`` -- full name of file

        OUTPUT:

        ``Zero`` if successful, otherwise ``nonzero``

        EXAMPLES::

            sage: from sage.numerical.backends.glpk_graph_backend import GLPKGraphBackend
            sage: gbe = GLPKGraphBackend()
            sage: import tempfile
            sage: with tempfile.NamedTemporaryFile() as f:
            ....:     gbe.write_maxflow(f.name)
            Traceback (most recent call last):
            ...
            OSError: Cannot write empty graph
            sage: gbe.add_vertices([None for i in range(2)])
            ['0', '1']
            sage: a = gbe.add_edge('0', '1')
            sage: gbe.maxflow_ffalg('0', '1')
            0.0
            sage: with tempfile.NamedTemporaryFile() as f:
            ....:     gbe.write_maxflow(f.name)
            Writing maximum flow problem data to ...
            6 lines were written
            0
        """

        if self.graph.nv <= 0:
            raise IOError("Cannot write empty graph")

        fname = str_to_bytes(fname, FS_ENCODING, 'surrogateescape')
        return glp_write_maxflow(self.graph, self.s+1, self.t+1,
                   sizeof(double), fname)

    cpdef double maxflow_ffalg(self, u=None, v=None) except -1:
        r"""
        Finds solution to the maxflow problem with Ford-Fulkerson algorithm.

        INPUT:

        - ``u`` -- name (as string) of the tail vertex; default is ``None``
        - ``v`` -- name (as string) of the head vertex; default is ``None``

        If ``u`` or ``v`` are ``None``, the currently stored values for the
        head or tail vertex are used. This behavior is useful when reading
        maxflow data from a file. When calling this function with values for
        ``u`` and ``v``, the head and tail vertex are stored for
        later use.

        OUTPUT: the solution to the maxflow problem, i.e. the maximum flow

        .. NOTE::

            * If the source or sink vertex does not exist, an :exc:`IndexError` is
              raised.

            * If the source and sink are identical, a :exc:`ValueError` is raised.

            * This method raises ``MIPSolverException`` exceptions when the
              solution cannot be computed for any reason (none exists, or the
              LP solver was not able to find it, etc...)

        EXAMPLES::

            sage: from sage.numerical.backends.glpk_graph_backend import GLPKGraphBackend
            sage: gbe = GLPKGraphBackend()
            sage: v = gbe.add_vertices([None for i in range(5)])
            sage: edges = ((0, 1, 2), (0, 2, 3), (1, 2, 3), (1, 3, 4),
            ....:         (3, 4, 1), (2, 4, 2))
            sage: for a in edges:
            ....:     edge = gbe.add_edge(str(a[0]), str(a[1]), {"cap":a[2]})
            sage: gbe.maxflow_ffalg('0', '4')
            3.0
            sage: gbe.maxflow_ffalg()
            3.0
            sage: gbe.maxflow_ffalg('0', '8')
            Traceback (most recent call last):
            ...
            IndexError: Source or sink vertex does not exist
        """
        cdef int s, t

        if u is not None and v is not None:
            s = self._find_vertex(u)
            t = self._find_vertex(v)
        else:
            s = self.s
            t = self.t

        if s < 0 or t < 0:
            raise IndexError("Source or sink vertex does not exist")
        if s == t:
            raise ValueError("Source and sink are identical")

        self.s = s
        self.t = t

        s += 1
        t += 1

        cdef double graph_sol
        cdef int status = glp_maxflow_ffalg(self.graph, s, t, sizeof(double),
                          &graph_sol, 3 * sizeof(double),
                          4 * sizeof(double))
        if status == 0:
            pass
        elif status == GLP_ENOPFS:
            raise MIPSolverException("No (primal) feasible solution exists")
        elif status == GLP_EDATA:
            raise MIPSolverException("Unable to start search due " +
                                     "to problem data")
        elif status == GLP_ERANGE:
            raise MIPSolverException("The search was prematurely terminated " +
                                     "because of integer overflow")
        elif status == GLP_EFAIL:
            raise MIPSolverException("An error has been detected in the " +
                                     "program logic")

        return graph_sol

    cpdef double cpp(self) noexcept:
        r"""
        Solve the critical path problem of a project network.

        OUTPUT: the length of the critical path of the network

        EXAMPLES::

            sage: from sage.numerical.backends.glpk_graph_backend import GLPKGraphBackend
            sage: gbe = GLPKGraphBackend()
            sage: gbe.add_vertices([None for i in range(3)])
            ['0', '1', '2']
            sage: gbe.set_vertex_demand('0', 3)
            sage: gbe.set_vertex_demand('1', 1)
            sage: gbe.set_vertex_demand('2', 4)
            sage: a = gbe.add_edge('0', '2')
            sage: a = gbe.add_edge('1', '2')
            sage: gbe.cpp()
            7.0
            sage: v = gbe.get_vertex('1')
            sage: 1, v["rhs"], v["es"], v["ls"] # abs tol 1e-6
            (1, 1.0, 0.0, 2.0)
        """

        return glp_cpp(self.graph, 0, 2 * sizeof(double),
                   3 * sizeof(double))

    def __dealloc__(self):
        """
        Destructor
        """
        if self.graph is not NULL:
            glp_delete_graph(self.graph)
        self.graph = NULL
