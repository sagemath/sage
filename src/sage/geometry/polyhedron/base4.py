# sage.doctest: optional - sage.graphs
r"""
Base class for polyhedra, part 4

Define methods relying on :mod:`sage.graphs`.
"""

# ****************************************************************************
#       Copyright (C) 2008-2012 Marshall Hampton <hamptonio@gmail.com>
#       Copyright (C) 2011-2015 Volker Braun <vbraun.name@gmail.com>
#       Copyright (C) 2012-2018 Frederic Chapoton
#       Copyright (C) 2013      Andrey Novoseltsev
#       Copyright (C) 2014-2017 Moritz Firsching
#       Copyright (C) 2014-2019 Thierry Monteil
#       Copyright (C) 2015      Nathann Cohen
#       Copyright (C) 2015-2017 Jeroen Demeyer
#       Copyright (C) 2015-2017 Vincent Delecroix
#       Copyright (C) 2015-2018 Dima Pasechnik
#       Copyright (C) 2015-2020 Jean-Philippe Labbe <labbe at math.huji.ac.il>
#       Copyright (C) 2015-2021 Matthias Koeppe
#       Copyright (C) 2016-2019 Daniel Krenn
#       Copyright (C) 2017      Marcelo Forets
#       Copyright (C) 2017-2018 Mark Bell
#       Copyright (C) 2019      Julian Ritter
#       Copyright (C) 2019-2020 Laith Rastanawi
#       Copyright (C) 2019-2020 Sophia Elia
#       Copyright (C) 2019-2021 Jonathan Kliem <jonathan.kliem@fu-berlin.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from .base3 import Polyhedron_base3

class Polyhedron_base4(Polyhedron_base3):
    """
    Methods relying on :mod:`sage.graphs`.

    See :class:`sage.geometry.polyhedron.base.Polyhedron_base`.

    TESTS::

        sage: from sage.geometry.polyhedron.base4 import Polyhedron_base4
        sage: P = polytopes.cube()
        sage: Polyhedron_base4.vertex_facet_graph.f(P)
        Digraph on 14 vertices
        sage: Polyhedron_base4.vertex_graph(P)
        Graph on 8 vertices
        sage: Polyhedron_base4.face_lattice(P)
        Finite lattice containing 28 elements
        sage: Polyhedron_base4.flag_f_vector(P, 0, 2)
        24
        sage: Polyhedron_base4.is_self_dual(P)
        False
        sage: Q = polytopes.cube(intervals='zero_one')
        sage: P == Q
        False
        sage: Polyhedron_base4.is_combinatorially_isomorphic(P, Q)
        True
    """

    @cached_method
    def vertex_facet_graph(self, labels=True):
        r"""
        Return the vertex-facet graph.

        This function constructs a directed bipartite graph.
        The nodes of the graph correspond to the vertices of the polyhedron
        and the facets of the polyhedron. There is an directed edge
        from a vertex to a face if and only if the vertex is incident to the face.

        INPUT:

        - ``labels`` -- boolean (default: ``True``); decide how the nodes
          of the graph are labelled. Either with the original vertices/facets
          of the Polyhedron or with integers.

        OUTPUT:

        - a bipartite DiGraph. If ``labels`` is ``True``, then the nodes
          of the graph will actually be the vertices and facets of ``self``,
          otherwise they will be integers.

        .. SEEALSO::

            :meth:`combinatorial_automorphism_group`,
            :meth:`is_combinatorially_isomorphic`.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: G = P.vertex_facet_graph(); G
            Digraph on 14 vertices
            sage: G.vertices(key = lambda v: str(v))
            [A vertex at (-1, -1, -1),
             A vertex at (-1, -1, 1),
             A vertex at (-1, 1, -1),
             A vertex at (-1, 1, 1),
             A vertex at (1, -1, -1),
             A vertex at (1, -1, 1),
             A vertex at (1, 1, -1),
             A vertex at (1, 1, 1),
             An inequality (-1, 0, 0) x + 1 >= 0,
             An inequality (0, -1, 0) x + 1 >= 0,
             An inequality (0, 0, -1) x + 1 >= 0,
             An inequality (0, 0, 1) x + 1 >= 0,
             An inequality (0, 1, 0) x + 1 >= 0,
             An inequality (1, 0, 0) x + 1 >= 0]
            sage: G.automorphism_group().is_isomorphic(P.hasse_diagram().automorphism_group())
            True
            sage: O = polytopes.octahedron(); O
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 6 vertices
            sage: O.vertex_facet_graph()
            Digraph on 14 vertices
            sage: H = O.vertex_facet_graph()
            sage: G.is_isomorphic(H)
            False
            sage: G2 = copy(G)
            sage: G2.reverse_edges(G2.edges())
            sage: G2.is_isomorphic(H)
            True

        TESTS:

        Check that :trac:`28828` is fixed::

            sage: G._immutable
            True

        Check that :trac:`29188` is fixed::

            sage: P = polytopes.cube()
            sage: P.vertex_facet_graph().is_isomorphic(P.vertex_facet_graph(False))
            True
        """
        return self.combinatorial_polyhedron().vertex_facet_graph(names=labels)

    def vertex_graph(self):
        """
        Return a graph in which the vertices correspond to vertices
        of the polyhedron, and edges to edges.

        ..NOTE::

            The graph of a polyhedron with lines has no vertices,
            as the polyhedron has no vertices (`0`-faces).

            The method :meth:`Polyhedron_base:vertices` returns
            the defining points in this case.

        EXAMPLES::

            sage: g3 = polytopes.hypercube(3).vertex_graph(); g3
            Graph on 8 vertices
            sage: g3.automorphism_group().cardinality()
            48
            sage: s4 = polytopes.simplex(4).vertex_graph(); s4
            Graph on 5 vertices
            sage: s4.is_eulerian()
            True

        The graph of an unbounded polyhedron
        is the graph of the bounded complex::

            sage: open_triangle = Polyhedron(vertices=[[1,0], [0,1]],
            ....:                            rays    =[[1,1]])
            sage: open_triangle.vertex_graph()
            Graph on 2 vertices

        The graph of a polyhedron with lines has no vertices::

            sage: line = Polyhedron(lines=[[0,1]])
            sage: line.vertex_graph()
            Graph on 0 vertices

        TESTS:

        Check for a line segment (:trac:`30545`)::

            sage: polytopes.simplex(1).graph().edges()
            [(A vertex at (0, 1), A vertex at (1, 0), None)]
        """
        return self.combinatorial_polyhedron().vertex_graph()

    graph = vertex_graph

    def vertex_digraph(self, f, increasing=True):
        r"""
        Return the directed graph of the polyhedron according to a linear form.

        The underlying undirected graph is the graph of vertices and edges.

        INPUT:

        - ``f`` -- a linear form. The linear form can be provided as:

            - a vector space morphism with one-dimensional codomain, (see
              :meth:`sage.modules.vector_space_morphism.linear_transformation`
              and
              :class:`sage.modules.vector_space_morphism.VectorSpaceMorphism`)
            - a vector ; in this case the linear form is obtained by duality
              using the dot product: ``f(v) = v.dot_product(f)``.

        - ``increasing`` -- boolean (default ``True``) whether to orient
          edges in the increasing or decreasing direction.

        By default, an edge is oriented from `v` to `w` if
        `f(v) \leq f(w)`.

        If `f(v)=f(w)`, then two opposite edges are created.

        EXAMPLES::

            sage: penta = Polyhedron([[0,0],[1,0],[0,1],[1,2],[3,2]])
            sage: G = penta.vertex_digraph(vector([1,1])); G
            Digraph on 5 vertices
            sage: G.sinks()
            [A vertex at (3, 2)]

            sage: A = matrix(ZZ, [[1], [-1]])
            sage: f = linear_transformation(A)
            sage: G = penta.vertex_digraph(f) ; G
            Digraph on 5 vertices
            sage: G.is_directed_acyclic()
            False

        .. SEEALSO::

            :meth:`vertex_graph`
        """
        from sage.modules.vector_space_morphism import VectorSpaceMorphism
        if isinstance(f, VectorSpaceMorphism):
            if f.codomain().dimension() == 1:
                orientation_check = lambda v: f(v) >= 0
            else:
                raise TypeError('the linear map f must have '
                                'one-dimensional codomain')
        else:
            try:
                if f.is_vector():
                    orientation_check = lambda v: v.dot_product(f) >= 0
                else:
                    raise TypeError('f must be a linear map or a vector')
            except AttributeError:
                raise TypeError('f must be a linear map or a vector')
        if not increasing:
            f = -f
        from sage.graphs.digraph import DiGraph
        dg = DiGraph()
        for j in range(self.n_vertices()):
            vj = self.Vrepresentation(j)
            for vi in vj.neighbors():
                if orientation_check(vj.vector() - vi.vector()):
                    dg.add_edge(vi, vj)
        return dg

    def face_lattice(self):
        """
        Return the face-lattice poset.

        OUTPUT:

        A :class:`~sage.combinat.posets.posets.FinitePoset`. Elements
        are given as
        :class:`~sage.geometry.polyhedron.face.PolyhedronFace`.

        In the case of a full-dimensional polytope, the faces are
        pairs (vertices, inequalities) of the spanning vertices and
        corresponding saturated inequalities. In general, a face is
        defined by a pair (V-rep. objects, H-rep. objects). The
        V-representation objects span the face, and the corresponding
        H-representation objects are those inequalities and equations
        that are saturated on the face.

        The bottom-most element of the face lattice is the "empty
        face". It contains no V-representation object. All
        H-representation objects are incident.

        The top-most element is the "full face". It is spanned by all
        V-representation objects. The incident H-representation
        objects are all equations and no inequalities.

        In the case of a full-dimensional polytope, the "empty face"
        and the "full face" are the empty set (no vertices, all
        inequalities) and the full polytope (all vertices, no
        inequalities), respectively.

        ALGORITHM:

        See :mod:`sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator`.

        .. NOTE::

            The face lattice is not cached, as long as this creates a memory leak, see :trac:`28982`.

        EXAMPLES::

            sage: square = polytopes.hypercube(2)
            sage: fl = square.face_lattice();fl
            Finite lattice containing 10 elements
            sage: list(f.ambient_V_indices() for f in fl)
            [(), (0,), (1,), (0, 1), (2,), (1, 2), (3,), (0, 3), (2, 3), (0, 1, 2, 3)]
            sage: poset_element = fl[5]
            sage: a_face = poset_element
            sage: a_face
            A 1-dimensional face of a Polyhedron in ZZ^2 defined as the convex hull of 2 vertices
            sage: a_face.ambient_V_indices()
            (1, 2)
            sage: set(a_face.ambient_Vrepresentation()) == \
            ....: set([square.Vrepresentation(1), square.Vrepresentation(2)])
            True
            sage: a_face.ambient_Vrepresentation()
            (A vertex at (1, 1), A vertex at (-1, 1))
            sage: a_face.ambient_Hrepresentation()
            (An inequality (0, -1) x + 1 >= 0,)

        A more complicated example::

            sage: c5_10 = Polyhedron(vertices = [[i,i^2,i^3,i^4,i^5] for i in range(1,11)])
            sage: c5_10_fl = c5_10.face_lattice()
            sage: [len(x) for x in c5_10_fl.level_sets()]
            [1, 10, 45, 100, 105, 42, 1]

        Note that if the polyhedron contains lines then there is a
        dimension gap between the empty face and the first non-empty
        face in the face lattice::

            sage: line = Polyhedron(vertices=[(0,)], lines=[(1,)])
            sage: [ fl.dim() for fl in line.face_lattice() ]
            [-1, 1]

        TESTS::

            sage: c5_20 = Polyhedron(vertices = [[i,i^2,i^3,i^4,i^5]
            ....:     for i in range(1,21)])
            sage: c5_20_fl = c5_20.face_lattice() # long time
            sage: [len(x) for x in c5_20_fl.level_sets()] # long time
            [1, 20, 190, 580, 680, 272, 1]
            sage: polytopes.hypercube(2).face_lattice().plot()  # optional - sage.plot
            Graphics object consisting of 27 graphics primitives
            sage: level_sets = polytopes.cross_polytope(2).face_lattice().level_sets()
            sage: level_sets[0][0].ambient_V_indices(), level_sets[-1][0].ambient_V_indices()
            ((), (0, 1, 2, 3))

        Various degenerate polyhedra::

            sage: [[ls.ambient_V_indices() for ls in lss] for lss in Polyhedron(vertices=[[0,0,0],[1,0,0],[0,1,0]]).face_lattice().level_sets()]
            [[()], [(0,), (1,), (2,)], [(0, 1), (0, 2), (1, 2)], [(0, 1, 2)]]
            sage: [[ls.ambient_V_indices() for ls in lss] for lss in Polyhedron(vertices=[(1,0,0),(0,1,0)], rays=[(0,0,1)]).face_lattice().level_sets()]
            [[()], [(1,), (2,)], [(0, 1), (0, 2), (1, 2)], [(0, 1, 2)]]
            sage: [[ls.ambient_V_indices() for ls in lss] for lss in Polyhedron(rays=[(1,0,0),(0,1,0)], vertices=[(0,0,1)]).face_lattice().level_sets()]
            [[()], [(0,)], [(0, 1), (0, 2)], [(0, 1, 2)]]
            sage: [[ls.ambient_V_indices() for ls in lss] for lss in Polyhedron(rays=[(1,0),(0,1)], vertices=[(0,0)]).face_lattice().level_sets()]
            [[()], [(0,)], [(0, 1), (0, 2)], [(0, 1, 2)]]
            sage: [[ls.ambient_V_indices() for ls in lss] for lss in Polyhedron(vertices=[(1,),(0,)]).face_lattice().level_sets()]
            [[()], [(0,), (1,)], [(0, 1)]]
            sage: [[ls.ambient_V_indices() for ls in lss] for lss in Polyhedron(vertices=[(1,0,0),(0,1,0)], lines=[(0,0,1)]).face_lattice().level_sets()]
            [[()], [(0, 1), (0, 2)], [(0, 1, 2)]]
            sage: [[ls.ambient_V_indices() for ls in lss] for lss in Polyhedron(lines=[(1,0,0)], vertices=[(0,0,1)]).face_lattice().level_sets()]
            [[()], [(0, 1)]]
            sage: [[ls.ambient_V_indices() for ls in lss] for lss in Polyhedron(lines=[(1,0),(0,1)], vertices=[(0,0)]).face_lattice().level_sets()]
            [[()], [(0, 1, 2)]]
            sage: [[ls.ambient_V_indices() for ls in lss] for lss in Polyhedron(lines=[(1,0)], rays=[(0,1)], vertices=[(0,0)]).face_lattice().level_sets()]
            [[()], [(0, 1)], [(0, 1, 2)]]
            sage: [[ls.ambient_V_indices() for ls in lss] for lss in Polyhedron(vertices=[(0,)], lines=[(1,)]).face_lattice().level_sets()]
            [[()], [(0, 1)]]
            sage: [[ls.ambient_V_indices() for ls in lss] for lss in Polyhedron(lines=[(1,0)], vertices=[(0,0)]).face_lattice().level_sets()]
            [[()], [(0, 1)]]

        """
        from sage.combinat.posets.lattices import FiniteLatticePoset
        return FiniteLatticePoset(self.hasse_diagram())

    @cached_method
    def hasse_diagram(self):
        r"""
        Return the Hasse diagram of the face lattice of ``self``.

        This is the Hasse diagram of the poset of the faces of ``self``.

        OUTPUT: a directed graph

        EXAMPLES::

            sage: P = polytopes.regular_polygon(4).pyramid()                    # optional - sage.rings.number_field
            sage: D = P.hasse_diagram(); D                                      # optional - sage.rings.number_field
            Digraph on 20 vertices
            sage: D.degree_polynomial()                                         # optional - sage.rings.number_field
            x^5 + x^4*y + x*y^4 + y^5 + 4*x^3*y + 8*x^2*y^2 + 4*x*y^3

        Faces of an mutable polyhedron are not hashable. Hence those are not suitable as
        vertices of the hasse diagram. Use the combinatorial polyhedron instead::

            sage: P = polytopes.regular_polygon(4).pyramid()                    # optional - sage.rings.number_field
            sage: parent = P.parent()                                           # optional - sage.rings.number_field
            sage: parent = parent.change_ring(QQ, backend='ppl')                # optional - sage.rings.number_field
            sage: Q = parent._element_constructor_(P, mutable=True)             # optional - sage.rings.number_field
            sage: Q.hasse_diagram()                                             # optional - sage.rings.number_field
            Traceback (most recent call last):
            ...
            TypeError: mutable polyhedra are unhashable
            sage: C = Q.combinatorial_polyhedron()                              # optional - sage.rings.number_field
            sage: D = C.hasse_diagram()                                         # optional - sage.rings.number_field
            sage: set(D.vertices()) == set(range(20))                           # optional - sage.rings.number_field
            True
            sage: def index_to_combinatorial_face(n):
            ....:     return C.face_by_face_lattice_index(n)
            sage: D.relabel(index_to_combinatorial_face, inplace=True)          # optional - sage.rings.number_field
            sage: D.vertices()                                                  # optional - sage.rings.number_field
            [A -1-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 0-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 0-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 0-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 0-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 0-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 1-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 1-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 1-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 1-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 1-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 1-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 1-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 1-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 2-dimensional face of a 3-dimensional combinatorial polyhedron,
             A 3-dimensional face of a 3-dimensional combinatorial polyhedron]
            sage: D.degree_polynomial()                                         # optional - sage.rings.number_field
            x^5 + x^4*y + x*y^4 + y^5 + 4*x^3*y + 8*x^2*y^2 + 4*x*y^3
        """

        from sage.geometry.polyhedron.face import combinatorial_face_to_polyhedral_face
        C = self.combinatorial_polyhedron()
        D = C.hasse_diagram()

        def index_to_polyhedron_face(n):
            return combinatorial_face_to_polyhedral_face(
                    self, C.face_by_face_lattice_index(n))

        return D.relabel(index_to_polyhedron_face, inplace=False, immutable=True)

    def flag_f_vector(self, *args):
        r"""
        Return the flag f-vector.

        For each `-1 < i_0 < \dots < i_n < d` the flag f-vector
        counts the number of flags `F_0 \subset \dots \subset F_n`
        with `F_j` of dimension `i_j` for each `0 \leq j \leq n`,
        where `d` is the dimension of the polyhedron.

        INPUT:

        - ``args`` -- integers (optional); specify an entry of the
          flag-f-vector; must be an increasing sequence of integers

        OUTPUT:

        - a dictionary, if no arguments were given

        - an Integer, if arguments were given

        EXAMPLES:

        Obtain the entire flag-f-vector::

            sage: P = polytopes.twenty_four_cell()
            sage: P.flag_f_vector()
                {(-1,): 1,
                 (0,): 24,
                 (0, 1): 192,
                 (0, 1, 2): 576,
                 (0, 1, 2, 3): 1152,
                 (0, 1, 3): 576,
                 (0, 2): 288,
                 (0, 2, 3): 576,
                 (0, 3): 144,
                 (1,): 96,
                 (1, 2): 288,
                 (1, 2, 3): 576,
                 (1, 3): 288,
                 (2,): 96,
                 (2, 3): 192,
                 (3,): 24,
                 (4,): 1}

        Specify an entry::

            sage: P.flag_f_vector(0,3)
            144
            sage: P.flag_f_vector(2)
            96

        Leading ``-1`` and trailing entry of dimension are allowed::

            sage: P.flag_f_vector(-1,0,3)
            144
            sage: P.flag_f_vector(-1,0,3,4)
            144

        One can get the number of trivial faces::

            sage: P.flag_f_vector(-1)
            1
            sage: P.flag_f_vector(4)
            1

        Polyhedra with lines, have ``0`` entries accordingly::

            sage: P = (Polyhedron(lines=[[1]]) * polytopes.cross_polytope(3))
            sage: P.flag_f_vector()
            {(-1,): 1,
             (0, 1): 0,
             (0, 1, 2): 0,
             (0, 1, 3): 0,
             (0, 2): 0,
             (0, 2, 3): 0,
             (0, 3): 0,
             (0,): 0,
             (1, 2): 24,
             (1, 2, 3): 48,
             (1, 3): 24,
             (1,): 6,
             (2, 3): 24,
             (2,): 12,
             (3,): 8,
             4: 1}

        If the arguments are not stricly increasing or out of range, a key error is raised::

            sage: P.flag_f_vector(-1,0,3,6)
            Traceback (most recent call last):
            ...
            KeyError: (0, 3, 6)
            sage: P.flag_f_vector(-1,3,0)
            Traceback (most recent call last):
            ...
            KeyError: (3, 0)
        """
        flag = self._flag_f_vector()
        if len(args) == 0:
            return flag
        elif len(args) == 1:
            return flag[(args[0],)]
        else:
            dim = self.dimension()
            if args[0] == -1:
                args = args[1:]
            if args[-1] == dim:
                args = args[:-1]
            return flag[tuple(args)]

    @cached_method(do_pickle=True)
    def _flag_f_vector(self):
        r"""
        Return the flag-f-vector.

        See :meth:`flag_f_vector`.

        TESTS::

            sage: polytopes.hypercube(4)._flag_f_vector()
            {(-1,): 1,
            (0,): 16,
            (0, 1): 64,
            (0, 1, 2): 192,
            (0, 1, 2, 3): 384,
            (0, 1, 3): 192,
            (0, 2): 96,
            (0, 2, 3): 192,
            (0, 3): 64,
            (1,): 32,
            (1, 2): 96,
            (1, 2, 3): 192,
            (1, 3): 96,
            (2,): 24,
            (2, 3): 48,
            (3,): 8,
            (4,): 1}
        """
        return self.combinatorial_polyhedron()._flag_f_vector()

    @cached_method
    def combinatorial_automorphism_group(self, vertex_graph_only=False):
        """
        Computes the combinatorial automorphism group.

        If ``vertex_graph_only`` is ``True``,  the automorphism group
        of the vertex-edge graph of the polyhedron is returned. Otherwise
        the automorphism group of the vertex-facet graph, which is
        isomorphic to the automorphism group of the face lattice is returned.

        INPUT:

        - ``vertex_graph_only`` -- boolean (default: ``False``); whether
          to return the automorphism group of the vertex edges graph or
          of the lattice

        OUTPUT:

        A
        :class:`PermutationGroup<sage.groups.perm_gps.permgroup.PermutationGroup_generic_with_category'>`
        that is isomorphic to the combinatorial automorphism group is
        returned.

        - if ``vertex_graph_only`` is ``True``:
          The automorphism group of the vertex-edge graph of the polyhedron

        - if ``vertex_graph_only`` is ``False`` (default):
          The automorphism group of the vertex-facet graph of the polyhedron,
          see :meth:`vertex_facet_graph`. This group is isomorphic to the
          automorphism group of the face lattice of the polyhedron.

        NOTE:

            Depending on ``vertex_graph_only``, this method returns groups
            that are not necessarily isomorphic, see the examples below.

        .. SEEALSO::

            :meth:`is_combinatorially_isomorphic`,
            :meth:`graph`,
            :meth:`vertex_facet_graph`.

        EXAMPLES::

            sage: quadrangle = Polyhedron(vertices=[(0,0),(1,0),(0,1),(2,3)])
            sage: quadrangle.combinatorial_automorphism_group().is_isomorphic(groups.permutation.Dihedral(4))
            True
            sage: quadrangle.restricted_automorphism_group()
            Permutation Group with generators [()]

        Permutations of the vertex graph only exchange vertices with vertices::

            sage: P = Polyhedron(vertices=[(1,0), (1,1)], rays=[(1,0)])
            sage: P.combinatorial_automorphism_group(vertex_graph_only=True)
            Permutation Group with generators [(A vertex at (1,0),A vertex at (1,1))]

        This shows an example of two polytopes whose vertex-edge graphs are isomorphic,
        but their face_lattices are not isomorphic::

            sage: Q=Polyhedron([[-123984206864/2768850730773, -101701330976/922950243591, -64154618668/2768850730773, -2748446474675/2768850730773],
            ....: [-11083969050/98314591817, -4717557075/98314591817, -32618537490/98314591817, -91960210208/98314591817],
            ....: [-9690950/554883199, -73651220/554883199, 1823050/554883199, -549885101/554883199], [-5174928/72012097, 5436288/72012097, -37977984/72012097, 60721345/72012097],
            ....: [-19184/902877, 26136/300959, -21472/902877, 899005/902877], [53511524/1167061933, 88410344/1167061933, 621795064/1167061933, 982203941/1167061933],
            ....: [4674489456/83665171433, -4026061312/83665171433, 28596876672/83665171433, -78383796375/83665171433], [857794884940/98972360190089, -10910202223200/98972360190089, 2974263671400/98972360190089, -98320463346111/98972360190089]])
            sage: C = polytopes.cyclic_polytope(4,8)
            sage: C.is_combinatorially_isomorphic(Q)
            False
            sage: C.combinatorial_automorphism_group(vertex_graph_only=True).is_isomorphic(Q.combinatorial_automorphism_group(vertex_graph_only=True))
            True
            sage: C.combinatorial_automorphism_group(vertex_graph_only=False).is_isomorphic(Q.combinatorial_automorphism_group(vertex_graph_only=False))
            False

        The automorphism group of the face lattice is isomorphic to the combinatorial automorphism group::

            sage: CG = C.hasse_diagram().automorphism_group()
            sage: C.combinatorial_automorphism_group().is_isomorphic(CG)
            True
            sage: QG = Q.hasse_diagram().automorphism_group()
            sage: Q.combinatorial_automorphism_group().is_isomorphic(QG)
            True

        """
        if vertex_graph_only:
            G = self.graph()
        else:
            G = self.vertex_facet_graph()
        return G.automorphism_group(edge_labels=True)

    @cached_method
    def restricted_automorphism_group(self, output="abstract"):
        r"""
        Return the restricted automorphism group.

        First, let the linear automorphism group be the subgroup of
        the affine group `AGL(d,\RR) = GL(d,\RR) \ltimes \RR^d`
        preserving the `d`-dimensional polyhedron. The affine group
        acts in the usual way `\vec{x}\mapsto A\vec{x}+b` on the
        ambient space.

        The restricted automorphism group is the subgroup of the linear
        automorphism group generated by permutations of the generators
        of the same type. That is, vertices can only be permuted with
        vertices, ray generators with ray generators, and line
        generators with line generators.

        For example, take the first quadrant

        .. MATH::

            Q = \Big\{ (x,y) \Big| x\geq 0,\; y\geq0 \Big\}
            \subset \QQ^2

        Then the linear automorphism group is

        .. MATH::

            \mathrm{Aut}(Q) =
            \left\{
            \begin{pmatrix}
            a & 0 \\ 0 & b
            \end{pmatrix}
            ,~
            \begin{pmatrix}
            0 & c \\ d & 0
            \end{pmatrix}
            :~
            a, b, c, d \in \QQ_{>0}
            \right\}
            \subset
            GL(2,\QQ)
            \subset
            E(d)

        Note that there are no translations that map the quadrant `Q`
        to itself, so the linear automorphism group is contained in
        the general linear group (the subgroup of transformations
        preserving the origin). The restricted automorphism group is

        .. MATH::

            \mathrm{Aut}(Q) =
            \left\{
            \begin{pmatrix}
            1 & 0 \\ 0 & 1
            \end{pmatrix}
            ,~
            \begin{pmatrix}
            0 & 1 \\ 1 & 0
            \end{pmatrix}
            \right\}
            \simeq \ZZ_2

        INPUT:

        - ``output`` -- how the group should be represented:

          - ``"abstract"`` (default) -- return an abstract permutation
            group without further meaning.

          - ``"permutation"`` -- return a permutation group on the
            indices of the polyhedron generators. For example, the
            permutation ``(0,1)`` would correspond to swapping
            ``self.Vrepresentation(0)`` and ``self.Vrepresentation(1)``.

          - ``"matrix"`` -- return a matrix group representing affine
            transformations. When acting on affine vectors, you should
            append a `1` to every vector. If the polyhedron is not full
            dimensional, the returned matrices act as the identity on
            the orthogonal complement of the affine space spanned by
            the polyhedron.

          - ``"matrixlist"`` -- like ``matrix``, but return the list of
            elements of the matrix group. Useful for fields without a
            good implementation of matrix groups or to avoid the
            overhead of creating the group.

        OUTPUT:

        - For ``output="abstract"`` and ``output="permutation"``:
          a :class:`PermutationGroup<sage.groups.perm_gps.permgroup.PermutationGroup_generic>`.

        - For ``output="matrix"``: a :class:`MatrixGroup`.

        - For ``output="matrixlist"``: a list of matrices.

        REFERENCES:

        - [BSS2009]_

        EXAMPLES:

        A cross-polytope example::

            sage: P = polytopes.cross_polytope(3)
            sage: P.restricted_automorphism_group() == PermutationGroup([[(3,4)], [(2,3),(4,5)],[(2,5)],[(1,2),(5,6)],[(1,6)]])
            True
            sage: P.restricted_automorphism_group(output="permutation") == PermutationGroup([[(2,3)],[(1,2),(3,4)],[(1,4)],[(0,1),(4,5)],[(0,5)]])
            True
            sage: mgens = [[[1,0,0,0],[0,1,0,0],[0,0,-1,0],[0,0,0,1]], [[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]], [[0,1,0,0],[1,0,0,0],[0,0,1,0],[0,0,0,1]]]

        We test groups for equality in a fool-proof way; they can have different generators, etc::

            sage: poly_g = P.restricted_automorphism_group(output="matrix")
            sage: matrix_g = MatrixGroup([matrix(QQ,t) for t in mgens])
            sage: all(t.matrix() in poly_g for t in matrix_g.gens())
            True
            sage: all(t.matrix() in matrix_g for t in poly_g.gens())
            True

        24-cell example::

            sage: P24 = polytopes.twenty_four_cell()
            sage: AutP24 = P24.restricted_automorphism_group()
            sage: PermutationGroup([
            ....:     '(1,20,2,24,5,23)(3,18,10,19,4,14)(6,21,11,22,7,15)(8,12,16,17,13,9)',
            ....:     '(1,21,8,24,4,17)(2,11,6,15,9,13)(3,20)(5,22)(10,16,12,23,14,19)'
            ....: ]).is_isomorphic(AutP24)
            True
            sage: AutP24.order()
            1152

        Here is the quadrant example mentioned in the beginning::

            sage: P = Polyhedron(rays=[(1,0),(0,1)])
            sage: P.Vrepresentation()
            (A vertex at (0, 0), A ray in the direction (0, 1), A ray in the direction (1, 0))
            sage: P.restricted_automorphism_group(output="permutation")
            Permutation Group with generators [(1,2)]

        Also, the polyhedron need not be full-dimensional::

            sage: P = Polyhedron(vertices=[(1,2,3,4,5),(7,8,9,10,11)])
            sage: P.restricted_automorphism_group()
            Permutation Group with generators [(1,2)]
            sage: G = P.restricted_automorphism_group(output="matrixlist")
            sage: G
            (
            [1 0 0 0 0 0]  [ -87/55  -82/55    -2/5   38/55   98/55   12/11]
            [0 1 0 0 0 0]  [-142/55  -27/55    -2/5   38/55   98/55   12/11]
            [0 0 1 0 0 0]  [-142/55  -82/55     3/5   38/55   98/55   12/11]
            [0 0 0 1 0 0]  [-142/55  -82/55    -2/5   93/55   98/55   12/11]
            [0 0 0 0 1 0]  [-142/55  -82/55    -2/5   38/55  153/55   12/11]
            [0 0 0 0 0 1], [      0       0       0       0       0       1]
            )
            sage: g = AffineGroup(5, QQ)(G[1])
            sage: g
                  [ -87/55  -82/55    -2/5   38/55   98/55]     [12/11]
                  [-142/55  -27/55    -2/5   38/55   98/55]     [12/11]
            x |-> [-142/55  -82/55     3/5   38/55   98/55] x + [12/11]
                  [-142/55  -82/55    -2/5   93/55   98/55]     [12/11]
                  [-142/55  -82/55    -2/5   38/55  153/55]     [12/11]
            sage: g^2
                  [1 0 0 0 0]     [0]
                  [0 1 0 0 0]     [0]
            x |-> [0 0 1 0 0] x + [0]
                  [0 0 0 1 0]     [0]
                  [0 0 0 0 1]     [0]
            sage: g(list(P.vertices()[0]))
            (7, 8, 9, 10, 11)
            sage: g(list(P.vertices()[1]))
            (1, 2, 3, 4, 5)

        Affine transformations do not change the restricted automorphism
        group. For example, any non-degenerate triangle has the
        dihedral group with 6 elements, `D_6`, as its automorphism
        group::

            sage: initial_points = [vector([1,0]), vector([0,1]), vector([-2,-1])]
            sage: points = initial_points
            sage: Polyhedron(vertices=points).restricted_automorphism_group()
            Permutation Group with generators [(2,3), (1,2)]
            sage: points = [pt - initial_points[0] for pt in initial_points]
            sage: Polyhedron(vertices=points).restricted_automorphism_group()
            Permutation Group with generators [(2,3), (1,2)]
            sage: points = [pt - initial_points[1] for pt in initial_points]
            sage: Polyhedron(vertices=points).restricted_automorphism_group()
            Permutation Group with generators [(2,3), (1,2)]
            sage: points = [pt - 2*initial_points[1] for pt in initial_points]
            sage: Polyhedron(vertices=points).restricted_automorphism_group()
            Permutation Group with generators [(2,3), (1,2)]

        The ``output="matrixlist"`` can be used over fields without a
        complete implementation of matrix groups::

            sage: P = polytopes.dodecahedron(); P
            A 3-dimensional polyhedron in (Number Field in sqrt5 with defining polynomial x^2 - 5 with sqrt5 = 2.236067977499790?)^3 defined as the convex hull of 20 vertices
            sage: G = P.restricted_automorphism_group(output="matrixlist")
            sage: len(G)
            120

        Floating-point computations are supported with a simple fuzzy
        zero implementation::

            sage: P = Polyhedron(vertices=[(1/3,0,0,1),(0,1/4,0,1),(0,0,1/5,1)], base_ring=RDF)
            sage: P.restricted_automorphism_group()
            Permutation Group with generators [(2,3), (1,2)]
            sage: len(P.restricted_automorphism_group(output="matrixlist"))
            6

        TESTS::

            sage: P = Polyhedron(vertices=[(1,0), (1,1)], rays=[(1,0)])
            sage: P.restricted_automorphism_group(output="permutation")
            Permutation Group with generators [(1,2)]
            sage: P.restricted_automorphism_group(output="matrix")
            Matrix group over Rational Field with 1 generators (
            [ 1  0  0]
            [ 0 -1  1]
            [ 0  0  1]
            )
            sage: P.restricted_automorphism_group(output="foobar")
            Traceback (most recent call last):
            ...
            ValueError: unknown output 'foobar', valid values are ('abstract', 'permutation', 'matrix', 'matrixlist')

        Check that :trac:`28828` is fixed::

            sage: P.restricted_automorphism_group(output="matrixlist")[0].is_immutable()
            True
        """
        # The algorithm works as follows:
        #
        # Let V be the matrix where every column is a homogeneous
        # coordinate of a V-representation object (vertex, ray, line).
        # Let us assume that V has full rank, that the polyhedron is
        # full dimensional.
        #
        # Let Q = V Vt and C = Vt Q^-1 V. The rows and columns of C
        # can be thought of as being indexed by the V-rep objects of the
        # polytope.
        #
        # It turns out that we can identify the restricted automorphism
        # group with the automorphism group of the edge-colored graph
        # on the V-rep objects with colors determined by the symmetric
        # matrix C.
        #
        # An automorphism of this graph is equivalent to a permutation
        # matrix P such that C = Pt C P. If we now define
        # A = V P Vt Q^-1, then one can check that V P = A V.
        # In other words: permuting the generators is the same as
        # applying the affine transformation A on the generators.
        #
        # If the given polyhedron is not fully-dimensional,
        # then Q will be not invertible. In this case, we use a
        # pseudoinverse Q+ instead of Q^-1. The formula for A acting on
        # the space spanned by V then simplifies to A = V P V+ where V+
        # denotes the pseudoinverse of V, which also equals V+ = Vt Q+.
        #
        # If we are asked to return the (group of) transformation
        # matrices to the user, we also require that those
        # transformations act as the identity on the orthogonal
        # complement of the space spanned by V. This complement is the
        # space spanned by the columns of W = 1 - V V+. One can check
        # that B = (V P V+) + W is the correct matrix: it acts the same
        # as A on V and it satisfies B W = W.

        outputs = ("abstract", "permutation", "matrix", "matrixlist")
        if output not in outputs:
            raise ValueError("unknown output {!r}, valid values are {}".format(output, outputs))

        # For backwards compatibility, we treat "abstract" as
        # "permutation", but where we add 1 to the indices of the
        # permutations.
        index0 = 0
        if output == "abstract":
            index0 = 1
            output = "permutation"

        if self.base_ring().is_exact():
            def rational_approximation(c):
                return c
        else:
            c_list = []

            def rational_approximation(c):
                # Implementation detail: Return unique integer if two
                # c-values are the same up to machine precision. But
                # you can think of it as a uniquely-chosen rational
                # approximation.
                for i, x in enumerate(c_list):
                    if self._is_zero(x - c):
                        return i
                c_list.append(c)
                return len(c_list) - 1

        if self.is_compact():
            def edge_label(i, j, c_ij):
                return c_ij
        else:
            # In the non-compact case, we also label the edges by the
            # type of the V-representation object. This ensures that
            # vertices, rays, and lines are only permuted amongst
            # themselves.
            def edge_label(i, j, c_ij):
                return (self.Vrepresentation(i).type(), c_ij, self.Vrepresentation(j).type())

        # Homogeneous coordinates for the V-representation objects.
        # Mathematically, V is a matrix. For efficiency however, we
        # represent it as a list of column vectors.
        V = [v.homogeneous_vector() for v in self.Vrepresentation()]

        # Pseudoinverse of V Vt
        Qplus = sum(v.column() * v.row() for v in V).pseudoinverse()

        # Construct the graph.
        from sage.graphs.graph import Graph
        G = Graph()
        for i in range(len(V)):
            for j in range(i+1, len(V)):
                c_ij = rational_approximation(V[i] * Qplus * V[j])
                G.add_edge(index0+i, index0+j, edge_label(i, j, c_ij))

        permgroup = G.automorphism_group(edge_labels=True)
        if output == "permutation":
            return permgroup
        elif output == "matrix":
            permgroup = permgroup.gens()

        # Compute V+ = Vt Q+ as list of row vectors
        from sage.matrix.constructor import matrix
        Vplus = list(matrix(V) * Qplus)  # matrix(V) is Vt

        # Compute W = 1 - V V+
        W = 1 - sum(V[i].column() * Vplus[i].row() for i in range(len(V)))

        # Convert the permutation group to a matrix group.
        # If P is a permutation, then we return the matrix
        # B = (V P V+) + W.
        #
        # If output == "matrix", we loop over the generators of the group.
        # Otherwise, we loop over all elements.
        matrices = []
        for perm in permgroup:
            A = sum(V[perm(i)].column() * Vplus[i].row() for i in range(len(V)))
            matrices.append(A + W)

        for mat in matrices:
            mat.set_immutable()

        if output == "matrixlist":
            return tuple(matrices)
        else:
            from sage.groups.matrix_gps.finitely_generated import MatrixGroup
            return MatrixGroup(matrices)

    def is_combinatorially_isomorphic(self, other, algorithm='bipartite_graph'):
        r"""
        Return whether the polyhedron is combinatorially isomorphic to another polyhedron.

        We only consider bounded polyhedra. By definition, they are
        combinatorially isomorphic if their face lattices are isomorphic.

        INPUT:

        - ``other`` -- a polyhedron object
        - ``algorithm`` (default = ``bipartite_graph``) -- the algorithm to use.
          The other possible value is ``face_lattice``.

        OUTPUT:

        - ``True`` if the two polyhedra are combinatorially isomorphic
        - ``False`` otherwise

        .. SEEALSO::

            :meth:`combinatorial_automorphism_group`,
            :meth:`vertex_facet_graph`.

        REFERENCES:

        For the equivalence of the two algorithms see [KK1995]_, p. 877-878

        EXAMPLES:

        The square is combinatorially isomorphic to the 2-dimensional cube::

            sage: polytopes.hypercube(2).is_combinatorially_isomorphic(polytopes.regular_polygon(4))
            True

        All the faces of the 3-dimensional permutahedron are either
        combinatorially isomorphic to a square or a hexagon::

            sage: H = polytopes.regular_polygon(6)                              # optional - sage.rings.number_field
            sage: S = polytopes.hypercube(2)
            sage: P = polytopes.permutahedron(4)
            sage: all(F.as_polyhedron().is_combinatorially_isomorphic(S)        # optional - sage.rings.number_field
            ....:       or F.as_polyhedron().is_combinatorially_isomorphic(H)
            ....:     for F in P.faces(2))
            True

        Checking that a regular simplex intersected with its reflection
        through the origin is combinatorially isomorphic to the intersection
        of a cube with a hyperplane perpendicular to its long diagonal::

            sage: def simplex_intersection(k):
            ....:   S1 = Polyhedron([vector(v)-vector(polytopes.simplex(k).center()) for v in polytopes.simplex(k).vertices_list()])
            ....:   S2 = Polyhedron([-vector(v) for v in S1.vertices_list()])
            ....:   return S1.intersection(S2)
            sage: def cube_intersection(k):
            ....:    C = polytopes.hypercube(k+1)
            ....:    H = Polyhedron(eqns=[[0]+[1 for i in range(k+1)]])
            ....:    return C.intersection(H)
            sage: [simplex_intersection(k).is_combinatorially_isomorphic(cube_intersection(k)) for k in range(2,5)]
            [True, True, True]
            sage: simplex_intersection(2).is_combinatorially_isomorphic(polytopes.regular_polygon(6))   # optional - sage.rings.number_field
            True
            sage: simplex_intersection(3).is_combinatorially_isomorphic(polytopes.octahedron())
            True

        Two polytopes with the same `f`-vector, but different combinatorial types::

            sage: P = Polyhedron([[-605520/1525633, -605520/1525633, -1261500/1525633, -52200/1525633, 11833/1525633],\
             [-720/1769, -600/1769, 1500/1769, 0, -31/1769], [-216/749, 240/749, -240/749, -432/749, 461/749], \
             [-50/181, 50/181, 60/181, -100/181, -119/181], [-32/51, -16/51, -4/51, 12/17, 1/17],\
             [1, 0, 0, 0, 0], [16/129, 128/129, 0, 0, 1/129], [64/267, -128/267, 24/89, -128/267, 57/89],\
             [1200/3953, -1200/3953, -1440/3953, -360/3953, -3247/3953], [1512/5597, 1512/5597, 588/5597, 4704/5597, 2069/5597]])
            sage: C = polytopes.cyclic_polytope(5,10)
            sage: C.f_vector() == P.f_vector(); C.f_vector()
            True
            (1, 10, 45, 100, 105, 42, 1)
            sage: C.is_combinatorially_isomorphic(P)
            False

            sage: S = polytopes.simplex(3)
            sage: S = S.face_truncation(S.faces(0)[3])
            sage: S = S.face_truncation(S.faces(0)[4])
            sage: S = S.face_truncation(S.faces(0)[5])
            sage: T = polytopes.simplex(3)
            sage: T = T.face_truncation(T.faces(0)[3])
            sage: T = T.face_truncation(T.faces(0)[4])
            sage: T = T.face_truncation(T.faces(0)[4])
            sage: T.is_combinatorially_isomorphic(S)
            False
            sage: T.f_vector(), S.f_vector()
            ((1, 10, 15, 7, 1), (1, 10, 15, 7, 1))

            sage: C = polytopes.hypercube(5)
            sage: C.is_combinatorially_isomorphic(C)
            True
            sage: C.is_combinatorially_isomorphic(C, algorithm='magic')
            Traceback (most recent call last):
            ...
            AssertionError: `algorithm` must be 'bipartite graph' or 'face_lattice'

            sage: G = Graph()
            sage: C.is_combinatorially_isomorphic(G)
            Traceback (most recent call last):
            ...
            AssertionError: input `other` must be a polyhedron

            sage: H = Polyhedron(eqns=[[0,1,1,1,1]]); H
            A 3-dimensional polyhedron in QQ^4 defined as the convex hull of 1 vertex and 3 lines
            sage: C.is_combinatorially_isomorphic(H)
            Traceback (most recent call last):
            ...
            AssertionError: polyhedron `other` must be bounded

        """
        assert isinstance(other, Polyhedron_base4), "input `other` must be a polyhedron"
        assert self.is_compact(), "polyhedron `self` must be bounded"
        assert other.is_compact(), "polyhedron `other` must be bounded"
        assert algorithm in ['bipartite_graph', 'face_lattice'], "`algorithm` must be 'bipartite graph' or 'face_lattice'"

        # For speed, we check if the polyhedra have the same number of facets and vertices.
        # This is faster than building the bipartite graphs first and
        # then check that they won't be isomorphic.
        if self.n_vertices() != other.n_vertices() or self.n_facets() != other.n_facets():
            return False

        if algorithm == 'bipartite_graph':
            G_self = self.vertex_facet_graph(False)
            G_other = other.vertex_facet_graph(False)

            return G_self.is_isomorphic(G_other)
        else:
            return self.face_lattice().is_isomorphic(other.face_lattice())

    def _test_is_combinatorially_isomorphic(self, tester=None, **options):
        """
        Run tests on the method :meth:`.is_combinatorially_isomorphic`.

        TESTS::

            sage: polytopes.cross_polytope(3)._test_is_combinatorially_isomorphic()
        """
        if tester is None:
            tester = self._tester(**options)

        if not self.is_compact():
            with tester.assertRaises(AssertionError):
                self.is_combinatorially_isomorphic(self)
            return

        if self.n_vertices() > 200 or self.n_facets() > 200:
            # Avoid very long doctests.
            return

        try:
            import sage.graphs.graph
        except ImportError:
            return

        from sage.rings.integer_ring import ZZ
        tester.assertTrue(self.is_combinatorially_isomorphic(ZZ(4)*self))
        if self.n_vertices():
            tester.assertTrue(self.is_combinatorially_isomorphic(self + self.center()))

        if self.n_vertices() < 20 and self.n_facets() < 20 and self.is_immutable():
            tester.assertTrue(self.is_combinatorially_isomorphic(ZZ(4)*self, algorithm='face_lattice'))
            if self.n_vertices():
                tester.assertTrue(self.is_combinatorially_isomorphic(self + self.center(), algorithm='face_lattice'))

    def is_self_dual(self):
        r"""
        Return whether the polytope is self-dual.

        A polytope is self-dual if its face lattice is isomorphic to the face
        lattice of its dual polytope.

        EXAMPLES::

            sage: polytopes.simplex().is_self_dual()
            True
            sage: polytopes.twenty_four_cell().is_self_dual()
            True
            sage: polytopes.cube().is_self_dual()
            False
            sage: polytopes.hypersimplex(5,2).is_self_dual()
            False
            sage: P = Polyhedron(vertices=[[1/2, 1/3]], rays=[[1, 1]]).is_self_dual()
            Traceback (most recent call last):
            ...
            ValueError: polyhedron has to be compact

        """
        if not self.is_compact():
            raise ValueError("polyhedron has to be compact")

        n = self.n_vertices()
        m = self.n_facets()
        if n != m:
            return False

        G1 = self.vertex_facet_graph()
        G2 = G1.reverse()
        return G1.is_isomorphic(G2)
