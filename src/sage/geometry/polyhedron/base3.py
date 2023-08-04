r"""
Base class for polyhedra: Methods regarding the combinatorics of a polyhedron

Excluding methods relying on :mod:`sage.graphs`.
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
#       Copyright (C) 2019-2021 Jonathan Kliem <jonathan.kliem@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from .base2 import Polyhedron_base2

class Polyhedron_base3(Polyhedron_base2):
    """
    Methods related to the combinatorics of a polyhedron.

    See :class:`sage.geometry.polyhedron.base.Polyhedron_base`.

    TESTS::

        sage: from sage.geometry.polyhedron.base3 import Polyhedron_base3
        sage: P = polytopes.cube()
        sage: Polyhedron_base3.is_simple(P)
        True
        sage: Polyhedron_base3.is_simplicial(P)
        False
        sage: Polyhedron_base3.is_prism(P)
        True
        sage: Polyhedron_base3.is_pyramid(P)
        False
        sage: Polyhedron_base3.combinatorial_polyhedron.f(P)
        A 3-dimensional combinatorial polyhedron with 6 facets
        sage: Polyhedron_base3.facets(P)
        (A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices,
        A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices,
        A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices,
        A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices,
        A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices,
        A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices)
        sage: Polyhedron_base3.f_vector.f(P)
        (1, 8, 12, 6, 1)
        sage: next(Polyhedron_base3.face_generator(P))
        A 3-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 8 vertices
    """

    def _init_empty_polyhedron(self):
        """
        Initializes an empty polyhedron.

        TESTS::

            sage: Polyhedron().vertex_adjacency_matrix()  # indirect doctest
            []
            sage: Polyhedron().facet_adjacency_matrix()
            [0]
        """
        Polyhedron_base2._init_empty_polyhedron(self)

        V_matrix = matrix(ZZ, 0, 0, 0)
        V_matrix.set_immutable()
        self.vertex_adjacency_matrix.set_cache(V_matrix)

        H_matrix = matrix(ZZ, 1, 1, 0)
        H_matrix.set_immutable()
        self.facet_adjacency_matrix.set_cache(H_matrix)

    @cached_method
    def slack_matrix(self):
        r"""
        Return the slack matrix.

        The entries correspond to the evaluation of the Hrepresentation
        elements on the  Vrepresentation elements.

        .. NOTE::

            The columns correspond to inequalities/equations in the
            order :meth:`~sage.geometry.polyhedron.base0.Polyhedron_base0.Hrepresentation`, the rows correspond to
            vertices/rays/lines in the order
            :meth:`~sage.geometry.polyhedron.base0.Polyhedron_base0.Vrepresentation`.

        .. SEEALSO::

            :meth:`incidence_matrix`.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: P.slack_matrix()
            [0 2 2 2 0 0]
            [0 0 2 2 0 2]
            [0 0 0 2 2 2]
            [0 2 0 2 2 0]
            [2 2 0 0 2 0]
            [2 2 2 0 0 0]
            [2 0 2 0 0 2]
            [2 0 0 0 2 2]

            sage: P = polytopes.cube(intervals='zero_one')
            sage: P.slack_matrix()
            [0 1 1 1 0 0]
            [0 0 1 1 0 1]
            [0 0 0 1 1 1]
            [0 1 0 1 1 0]
            [1 1 0 0 1 0]
            [1 1 1 0 0 0]
            [1 0 1 0 0 1]
            [1 0 0 0 1 1]

            sage: # needs sage.rings.number_field
            sage: P = polytopes.dodecahedron().faces(2)[0].as_polyhedron()
            sage: P.slack_matrix()
            [1/2*sqrt5 - 1/2               0               0               1 1/2*sqrt5 - 1/2               0]
            [              0               0 1/2*sqrt5 - 1/2 1/2*sqrt5 - 1/2               1               0]
            [              0 1/2*sqrt5 - 1/2               1               0 1/2*sqrt5 - 1/2               0]
            [              1 1/2*sqrt5 - 1/2               0 1/2*sqrt5 - 1/2               0               0]
            [1/2*sqrt5 - 1/2               1 1/2*sqrt5 - 1/2               0               0               0]

            sage: P = Polyhedron(rays=[[1, 0], [0, 1]])
            sage: P.slack_matrix()
            [0 0]
            [0 1]
            [1 0]

        TESTS::

            sage: Polyhedron().slack_matrix()
            []
            sage: Polyhedron(base_ring=QuadraticField(2)).slack_matrix().base_ring()    # needs sage.rings.number_field
            Number Field in a with defining polynomial x^2 - 2 with a = 1.41...
        """
        if not self.n_Vrepresentation() or not self.n_Hrepresentation():
            slack_matrix = matrix(self.base_ring(), self.n_Vrepresentation(),
                                  self.n_Hrepresentation(), 0)
        else:
            Vrep_matrix = matrix(self.base_ring(), self.Vrepresentation())
            Hrep_matrix = matrix(self.base_ring(), self.Hrepresentation())

            # Getting homogeneous coordinates of the Vrepresentation.
            hom_helper = matrix(self.base_ring(), [1 if v.is_vertex() else 0 for v in self.Vrepresentation()])
            hom_Vrep = hom_helper.stack(Vrep_matrix.transpose())

            slack_matrix = (Hrep_matrix * hom_Vrep).transpose()

        slack_matrix.set_immutable()
        return slack_matrix

    @cached_method
    def incidence_matrix(self):
        """
        Return the incidence matrix.

        .. NOTE::

            The columns correspond to inequalities/equations in the
            order :meth:`~sage.geometry.polyhedron.base0.Polyhedron_base0.Hrepresentation`, the rows correspond to
            vertices/rays/lines in the order
            :meth:`~sage.geometry.polyhedron.base0.Polyhedron_base0.Vrepresentation`.

        .. SEEALSO::

            :meth:`slack_matrix`.

        EXAMPLES::

            sage: p = polytopes.cuboctahedron()
            sage: p.incidence_matrix()
            [0 0 1 1 0 1 0 0 0 0 1 0 0 0]
            [0 0 0 1 0 0 1 0 1 0 1 0 0 0]
            [0 0 1 1 1 0 0 1 0 0 0 0 0 0]
            [1 0 0 1 1 0 1 0 0 0 0 0 0 0]
            [0 0 0 0 0 1 0 0 1 1 1 0 0 0]
            [0 0 1 0 0 1 0 1 0 0 0 1 0 0]
            [1 0 0 0 0 0 1 0 1 0 0 0 1 0]
            [1 0 0 0 1 0 0 1 0 0 0 0 0 1]
            [0 1 0 0 0 1 0 0 0 1 0 1 0 0]
            [0 1 0 0 0 0 0 0 1 1 0 0 1 0]
            [0 1 0 0 0 0 0 1 0 0 0 1 0 1]
            [1 1 0 0 0 0 0 0 0 0 0 0 1 1]
            sage: v = p.Vrepresentation(0)
            sage: v
            A vertex at (-1, -1, 0)
            sage: h = p.Hrepresentation(2)
            sage: h
            An inequality (1, 1, -1) x + 2 >= 0
            sage: h.eval(v)        # evaluation (1, 1, -1) * (-1/2, -1/2, 0) + 1
            0
            sage: h*v              # same as h.eval(v)
            0
            sage: p.incidence_matrix() [0,2]   # this entry is (v,h)
            1
            sage: h.contains(v)
            True
            sage: p.incidence_matrix() [2,0]   # note: not symmetric
            0

        The incidence matrix depends on the ambient dimension::

            sage: simplex = polytopes.simplex(); simplex
            A 3-dimensional polyhedron in ZZ^4 defined as the convex hull of 4 vertices
            sage: simplex.incidence_matrix()
            [1 1 1 1 0]
            [1 1 1 0 1]
            [1 1 0 1 1]
            [1 0 1 1 1]
            sage: simplex = simplex.affine_hull_projection(); simplex
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 4 vertices
            sage: simplex.incidence_matrix()
            [1 1 1 0]
            [1 1 0 1]
            [1 0 1 1]
            [0 1 1 1]

        An incidence matrix does not determine a unique
        polyhedron::

            sage: P = Polyhedron(vertices=[[0,1],[1,1],[1,0]])
            sage: P.incidence_matrix()
            [1 1 0]
            [1 0 1]
            [0 1 1]

            sage: Q = Polyhedron(vertices=[[0,1], [1,0]], rays=[[1,1]])
            sage: Q.incidence_matrix()
            [1 1 0]
            [1 0 1]
            [0 1 1]


        An example of two polyhedra with isomorphic face lattices
        but different incidence matrices::

            sage: Q.incidence_matrix()
            [1 1 0]
            [1 0 1]
            [0 1 1]

            sage: R = Polyhedron(vertices=[[0,1], [1,0]], rays=[[1,3/2], [3/2,1]])
            sage: R.incidence_matrix()
            [1 1 0]
            [1 0 1]
            [0 1 0]
            [0 0 1]

        The incidence matrix has base ring integers. This way one can express various
        counting questions::

            sage: P = polytopes.twenty_four_cell()
            sage: M = P.incidence_matrix()
            sage: sum(sum(x) for x in M) == P.flag_f_vector(0, 3)                       # needs sage.combinat
            True

        TESTS:

        Check that :trac:`28828` is fixed::

            sage: R.incidence_matrix().is_immutable()
            True

        Test that this method works for inexact base ring
        (``cdd`` sets the cache already)::

            sage: # needs sage.groups
            sage: P = polytopes.dodecahedron(exact=False)
            sage: M = P.incidence_matrix.cache
            sage: P.incidence_matrix.clear_cache()
            sage: M == P.incidence_matrix()
            True
        """
        if self.base_ring() in (ZZ, QQ):
            # Much faster for integers or rationals.
            incidence_matrix = self.slack_matrix().zero_pattern_matrix(ZZ)
            incidence_matrix.set_immutable()
            return incidence_matrix

        incidence_matrix = matrix(ZZ, self.n_Vrepresentation(),
                                  self.n_Hrepresentation(), 0)

        Vvectors_vertices = tuple((v.vector(), v.index())
                                  for v in self.Vrep_generator()
                                  if v.is_vertex())
        Vvectors_rays_lines = tuple((v.vector(), v.index())
                                    for v in self.Vrep_generator()
                                    if not v.is_vertex())

        # Determine ``is_zero`` to save lots of time.
        if self.base_ring().is_exact():
            def is_zero(x):
                return not x
        else:
            is_zero = self._is_zero

        for H in self.Hrep_generator():
            Hconst = H.b()
            Hvec = H.A()
            Hindex = H.index()
            for Vvec, Vindex in Vvectors_vertices:
                if is_zero(Hvec*Vvec + Hconst):
                    incidence_matrix[Vindex, Hindex] = 1

            # A ray or line is considered incident with a hyperplane,
            # if it is orthogonal to the normal vector of the hyperplane.
            for Vvec, Vindex in Vvectors_rays_lines:
                if is_zero(Hvec*Vvec):
                    incidence_matrix[Vindex, Hindex] = 1

        incidence_matrix.set_immutable()
        return incidence_matrix

    @cached_method
    def combinatorial_polyhedron(self):
        r"""
        Return the combinatorial type of ``self``.

        See :class:`sage.geometry.polyhedron.combinatorial_polyhedron.base.CombinatorialPolyhedron`.

        EXAMPLES::

            sage: polytopes.cube().combinatorial_polyhedron()
            A 3-dimensional combinatorial polyhedron with 6 facets

            sage: polytopes.cyclic_polytope(4,10).combinatorial_polyhedron()
            A 4-dimensional combinatorial polyhedron with 35 facets

            sage: Polyhedron(rays=[[0,1], [1,0]]).combinatorial_polyhedron()
            A 2-dimensional combinatorial polyhedron with 2 facets
        """
        from sage.geometry.polyhedron.combinatorial_polyhedron.base import CombinatorialPolyhedron
        return CombinatorialPolyhedron(self)

    def _test_combinatorial_polyhedron(self, tester=None, **options):
        """
        Run test suite of combinatorial polyhedron.

        TESTS::

            sage: polytopes.cross_polytope(3)._test_combinatorial_polyhedron()
        """
        from sage.misc.sage_unittest import TestSuite

        tester = self._tester(tester=tester, **options)
        tester.info("\n  Running the test suite of self.combinatorial_polyhedron()")
        TestSuite(self.combinatorial_polyhedron()).run(verbose=tester._verbose,
                                                       prefix=tester._prefix+"  ")
        tester.info(tester._prefix + " ", newline=False)

    def face_generator(self, face_dimension=None, algorithm=None, **kwds):
        r"""
        Return an iterator over the faces of given dimension.

        If dimension is not specified return an iterator over all faces.

        INPUT:

        - ``face_dimension`` -- integer (default ``None``),
          yield only faces of this dimension if specified

        - ``algorithm`` -- string (optional);
          specify whether to start with facets or vertices:

          * ``'primal'`` -- start with the facets
          * ``'dual'`` -- start with the vertices
          * ``None`` -- choose automatically

        OUTPUT:

        A :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator.FaceIterator_geom`.
        This class iterates over faces as
        :class:`~sage.geometry.polyhedron.face.PolyhedronFace`. See
        :mod:`~sage.geometry.polyhedron.face` for details. The order
        is random but fixed.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: it = P.face_generator()
            sage: it
            Iterator over the faces of a 3-dimensional polyhedron in ZZ^3
            sage: list(it)
            [A 3-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 8 vertices,
             A -1-dimensional face of a Polyhedron in ZZ^3,
             A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices,
             A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices,
             A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices,
             A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices,
             A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices,
             A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex,
             A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex,
             A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex,
             A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex,
             A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices]

             sage: P = polytopes.hypercube(4)
             sage: list(P.face_generator(2))[:4]
             [A 2-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 4 vertices,
              A 2-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 4 vertices,
              A 2-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 4 vertices,
              A 2-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 4 vertices]

        If a polytope has more facets than vertices, the dual mode is chosen::

            sage: P = polytopes.cross_polytope(3)
            sage: list(P.face_generator())
            [A 3-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 6 vertices,
             A -1-dimensional face of a Polyhedron in ZZ^3,
             A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex,
             A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex,
             A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex,
             A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex,
             A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex,
             A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 3 vertices,
             A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 3 vertices,
             A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 3 vertices,
             A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 3 vertices,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 3 vertices,
             A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 3 vertices,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 3 vertices,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 3 vertices,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices]

        The face iterator can also be slightly modified.
        In non-dual mode we can skip subfaces of the current (proper) face::

            sage: P = polytopes.cube()
            sage: it = P.face_generator(algorithm='primal')
            sage: _ = next(it), next(it)
            sage: face = next(it)
            sage: face.ambient_H_indices()
            (5,)
            sage: it.ignore_subfaces()
            sage: face = next(it)
            sage: face.ambient_H_indices()
            (4,)
            sage: it.ignore_subfaces()
            sage: [face.ambient_H_indices() for face in it]
            [(3,),
             (2,),
             (1,),
             (0,),
             (2, 3),
             (1, 3),
             (1, 2, 3),
             (1, 2),
             (0, 2),
             (0, 1, 2),
             (0, 1)]

        In dual mode we can skip supfaces of the current (proper) face::

            sage: P = polytopes.cube()
            sage: it = P.face_generator(algorithm='dual')
            sage: _ = next(it), next(it)
            sage: face = next(it)
            sage: face.ambient_V_indices()
            (7,)
            sage: it.ignore_supfaces()
            sage: next(it)
            A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex
            sage: face = next(it)
            sage: face.ambient_V_indices()
            (5,)
            sage: it.ignore_supfaces()
            sage: [face.ambient_V_indices() for face in it]
            [(4,),
             (3,),
             (2,),
             (1,),
             (0,),
             (1, 6),
             (3, 4),
             (2, 3),
             (0, 3),
             (0, 1, 2, 3),
             (1, 2),
             (0, 1)]

        In non-dual mode, we cannot skip supfaces::

            sage: it = P.face_generator(algorithm='primal')
            sage: _ = next(it), next(it)
            sage: next(it)
            A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices
            sage: it.ignore_supfaces()
            Traceback (most recent call last):
            ...
            ValueError: only possible when in dual mode

        In dual mode, we cannot skip subfaces::

            sage: it = P.face_generator(algorithm='dual')
            sage: _ = next(it), next(it)
            sage: next(it)
            A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex
            sage: it.ignore_subfaces()
            Traceback (most recent call last):
            ...
            ValueError: only possible when not in dual mode

        We can only skip sub-/supfaces of proper faces::

            sage: it = P.face_generator(algorithm='primal')
            sage: next(it)
            A 3-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 8 vertices
            sage: it.ignore_subfaces()
            Traceback (most recent call last):
            ...
            ValueError: iterator not set to a face yet

        .. SEEALSO::

            :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator.FaceIterator_geom`.

        ALGORITHM:

        See :class:`~sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator.FaceIterator`.

        TESTS::

            sage: P = polytopes.simplex()
            sage: list(P.face_generator(-2))
            []
            sage: list(P.face_generator(-1))
            [A -1-dimensional face of a Polyhedron in ZZ^4]
            sage: list(P.face_generator(3))
            [A 3-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 4 vertices]

            sage: list(Polyhedron().face_generator())
            [A -1-dimensional face of a Polyhedron in ZZ^0]

        Check that :trac:`29155` is fixed::

            sage: P = polytopes.permutahedron(3)
            sage: [f] = P.face_generator(2)
            sage: f.ambient_Hrepresentation()
            (An equation (1, 1, 1) x - 6 == 0,)

        The ``dual`` keyword is deprecated::

             sage: P = polytopes.hypercube(4)
             sage: list(P.face_generator(dual=False))[:4]
             doctest:...: DeprecationWarning: the keyword dual is deprecated; use algorithm instead
             See https://github.com/sagemath/sage/issues/33646 for details.
             [A 4-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 16 vertices,
              A -1-dimensional face of a Polyhedron in ZZ^4,
              A 3-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 8 vertices,
              A 3-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 8 vertices]
             sage: list(P.face_generator(True))[:4]
             [A 1-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 2 vertices,
              A 1-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 2 vertices,
              A 1-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 2 vertices,
              A 1-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 2 vertices]

        Check that we catch incorrect algorithms::

             sage: list(P.face_generator(2, algorithm='integrate'))[:4]
             Traceback (most recent call last):
             ...
             ValueError: algorithm must be 'primal', 'dual' or None
        """
        dual = None
        if algorithm == 'primal':
            dual = False
        elif algorithm == 'dual':
            dual = True
        elif algorithm in (False, True):
            from sage.misc.superseded import deprecation
            deprecation(33646, "the keyword dual is deprecated; use algorithm instead")
            dual = algorithm
        elif algorithm is not None:
            raise ValueError("algorithm must be 'primal', 'dual' or None")

        if kwds:
            from sage.misc.superseded import deprecation
            deprecation(33646, "the keyword dual is deprecated; use algorithm instead")
            if 'dual' in kwds and dual is None:
                dual = kwds['dual']

        from sage.geometry.polyhedron.combinatorial_polyhedron.face_iterator import FaceIterator_geom
        return FaceIterator_geom(self, output_dimension=face_dimension, dual=dual)

    def faces(self, face_dimension):
        """
        Return the faces of given dimension

        INPUT:

        - ``face_dimension`` -- integer. The dimension of the faces
          whose representation will be returned.

        OUTPUT:

        A tuple of
        :class:`~sage.geometry.polyhedron.face.PolyhedronFace`. See
        module :mod:`sage.geometry.polyhedron.face` for details. The order
        is random but fixed.

        .. SEEALSO::

            :meth:`face_generator`,
            :meth:`~sage.geometry.polyhedron.base0.Polyhedron_base0.facet`.

        EXAMPLES:

        Here we find the vertex and face indices of the eight three-dimensional
        facets of the four-dimensional hypercube::

            sage: p = polytopes.hypercube(4)
            sage: list(f.ambient_V_indices() for f in p.faces(3))
            [(0, 5, 6, 7, 8, 9, 14, 15),
             (1, 4, 5, 6, 10, 13, 14, 15),
             (1, 2, 6, 7, 8, 10, 11, 15),
             (8, 9, 10, 11, 12, 13, 14, 15),
             (0, 3, 4, 5, 9, 12, 13, 14),
             (0, 2, 3, 7, 8, 9, 11, 12),
             (1, 2, 3, 4, 10, 11, 12, 13),
             (0, 1, 2, 3, 4, 5, 6, 7)]

            sage: face = p.faces(3)[3]
            sage: face.ambient_Hrepresentation()
            (An inequality (1, 0, 0, 0) x + 1 >= 0,)
            sage: face.vertices()
            (A vertex at (-1, -1, 1, -1),
             A vertex at (-1, -1, 1, 1),
             A vertex at (-1, 1, -1, -1),
             A vertex at (-1, 1, 1, -1),
             A vertex at (-1, 1, 1, 1),
             A vertex at (-1, 1, -1, 1),
             A vertex at (-1, -1, -1, 1),
             A vertex at (-1, -1, -1, -1))

        You can use the
        :meth:`~sage.geometry.polyhedron.representation.PolyhedronRepresentation.index`
        method to enumerate vertices and inequalities::

            sage: def get_idx(rep): return rep.index()
            sage: [get_idx(_) for _ in face.ambient_Hrepresentation()]
            [4]
            sage: [get_idx(_) for _ in face.ambient_Vrepresentation()]
            [8, 9, 10, 11, 12, 13, 14, 15]

            sage: [ ([get_idx(_) for _ in face.ambient_Vrepresentation()],
            ....:    [get_idx(_) for _ in face.ambient_Hrepresentation()])
            ....:   for face in p.faces(3) ]
            [([0, 5, 6, 7, 8, 9, 14, 15], [7]),
             ([1, 4, 5, 6, 10, 13, 14, 15], [6]),
             ([1, 2, 6, 7, 8, 10, 11, 15], [5]),
             ([8, 9, 10, 11, 12, 13, 14, 15], [4]),
             ([0, 3, 4, 5, 9, 12, 13, 14], [3]),
             ([0, 2, 3, 7, 8, 9, 11, 12], [2]),
             ([1, 2, 3, 4, 10, 11, 12, 13], [1]),
             ([0, 1, 2, 3, 4, 5, 6, 7], [0])]

        TESTS::

            sage: pr = Polyhedron(rays = [[1,0,0],[-1,0,0],[0,1,0]], vertices = [[-1,-1,-1]], lines=[(0,0,1)])
            sage: pr.faces(4)
            ()
            sage: pr.faces(3)[0].ambient_V_indices()
            (0, 1, 2, 3)
            sage: pr.facets()[0].ambient_V_indices()
            (0, 1, 2)
            sage: pr.faces(1)
            ()
            sage: pr.faces(0)
            ()
            sage: pr.faces(-1)
            (A -1-dimensional face of a Polyhedron in QQ^3,)
        """
        return tuple(self.face_generator(face_dimension))

    def facets(self):
        r"""
        Return the facets of the polyhedron.

        Facets are the maximal nontrivial faces of polyhedra.
        The empty face and the polyhedron itself are trivial.

        A facet of a `d`-dimensional polyhedron is a face of dimension
        `d-1`. For `d \neq 0` the converse is true as well.

        OUTPUT:

        A tuple of
        :class:`~sage.geometry.polyhedron.face.PolyhedronFace`. See
        :mod:`~sage.geometry.polyhedron.face` for details. The order
        is random but fixed.

        .. SEEALSO:: :meth:`facets`

        EXAMPLES:

        Here we find the eight three-dimensional facets of the
        four-dimensional hypercube::

            sage: p = polytopes.hypercube(4)
            sage: p.facets()
            (A 3-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 8 vertices,
             A 3-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 8 vertices,
             A 3-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 8 vertices,
             A 3-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 8 vertices,
             A 3-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 8 vertices,
             A 3-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 8 vertices,
             A 3-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 8 vertices,
             A 3-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 8 vertices)

        This is the same result as explicitly finding the
        three-dimensional faces::

            sage: dim = p.dimension()
            sage: p.faces(dim-1)
            (A 3-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 8 vertices,
             A 3-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 8 vertices,
             A 3-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 8 vertices,
             A 3-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 8 vertices,
             A 3-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 8 vertices,
             A 3-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 8 vertices,
             A 3-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 8 vertices,
             A 3-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 8 vertices)

        The ``0``-dimensional polyhedron does not have facets::

            sage: P = Polyhedron([[0]])
            sage: P.facets()
            ()
        """
        if self.dimension() == 0:
            return ()
        return self.faces(self.dimension()-1)

    @cached_method(do_pickle=True, key=lambda self, x, y, z: None)
    def f_vector(self, num_threads=None, parallelization_depth=None, algorithm=None):
        r"""
        Return the f-vector.

        INPUT:

        - ``num_threads`` -- integer (optional); specify the number of threads;
          otherwise determined by :func:`~sage.parallel.ncpus.ncpus`

        - ``parallelization_depth`` -- integer (optional); specify
          how deep in the lattice the parallelization is done

        - ``algorithm`` -- string (optional);
          specify whether the face generator starts with facets or vertices:

          * ``'primal'`` -- start with the facets
          * ``'dual'`` -- start with the vertices
          * ``None`` -- choose automatically

        OUTPUT:

        Return a vector whose `i`-th entry is the number of
        `i-2`-dimensional faces of the polytope.

        .. NOTE::

            The ``vertices`` as given by :meth:`~sage.geometry.polyhedron.base0.Polyhedron_base0.vertices`
            do not need to correspond to `0`-dimensional faces. If a polyhedron
            contains `k` lines they correspond to `k`-dimensional faces.
            See example below.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[[1, 2, 3], [1, 3, 2],
            ....:     [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1], [0, 0, 0]])
            sage: p.f_vector()
            (1, 7, 12, 7, 1)

            sage: polytopes.cyclic_polytope(4,10).f_vector()
            (1, 10, 45, 70, 35, 1)

            sage: polytopes.hypercube(5).f_vector()
            (1, 32, 80, 80, 40, 10, 1)

        Polyhedra with lines do not have `0`-faces::

            sage: Polyhedron(ieqs=[[1,-1,0,0],[1,1,0,0]]).f_vector()
            (1, 0, 0, 2, 1)

        However, the method :meth:`Polyhedron_base.vertices` returns
        two points that belong to the ``Vrepresentation``::

            sage: P = Polyhedron(ieqs=[[1,-1,0],[1,1,0]])
            sage: P.vertices()
            (A vertex at (1, 0), A vertex at (-1, 0))
            sage: P.f_vector()
            (1, 0, 2, 1)

        TESTS:

        Check that :trac:`28828` is fixed::

            sage: P.f_vector().is_immutable()
            True

        The cache of the f-vector is being pickled::

            sage: P = polytopes.cube()
            sage: P.f_vector()
            (1, 8, 12, 6, 1)
            sage: Q = loads(dumps(P))
            sage: Q.f_vector.is_in_cache()
            True
        """
        return self.combinatorial_polyhedron().f_vector(num_threads, parallelization_depth, algorithm=algorithm)

    def bounded_edges(self):
        """
        Return the bounded edges (excluding rays and lines).

        OUTPUT:

        A generator for pairs of vertices, one pair per edge.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[[1,0],[0,1]], rays=[[1,0],[0,1]])
            sage: [ e for e in p.bounded_edges() ]
            [(A vertex at (0, 1), A vertex at (1, 0))]
            sage: for e in p.bounded_edges(): print(e)
            (A vertex at (0, 1), A vertex at (1, 0))
        """
        obj = self.Vrepresentation()
        for i in range(len(obj)):
            if not obj[i].is_vertex():
                continue
            for j in range(i+1, len(obj)):
                if not obj[j].is_vertex():
                    continue
                if self.vertex_adjacency_matrix()[i, j] == 0:
                    continue
                yield (obj[i], obj[j])

    @cached_method
    def vertex_adjacency_matrix(self, algorithm=None):
        """
        Return the binary matrix of vertex adjacencies.

        INPUT:

        - ``algorithm`` -- string (optional);
          specify whether the face generator starts with facets or vertices:

          * ``'primal'`` -- start with the facets
          * ``'dual'`` -- start with the vertices
          * ``None`` -- choose automatically

        EXAMPLES::

            sage: polytopes.simplex(4).vertex_adjacency_matrix()
            [0 1 1 1 1]
            [1 0 1 1 1]
            [1 1 0 1 1]
            [1 1 1 0 1]
            [1 1 1 1 0]

        The rows and columns of the vertex adjacency matrix correspond
        to the :meth:`Vrepresentation` objects: vertices, rays, and
        lines. The `(i,j)` matrix entry equals `1` if the `i`-th and
        `j`-th V-representation object are adjacent.

        Two vertices are adjacent if they are the endpoints of an
        edge, that is, a one-dimensional face. For unbounded polyhedra
        this clearly needs to be generalized and we define two
        V-representation objects (see
        :mod:`sage.geometry.polyhedron.constructor`) to be adjacent if
        they together generate a one-face. There are three possible
        combinations:

        * Two vertices can bound a finite-length edge.

        * A vertex and a ray can generate a half-infinite edge
          starting at the vertex and with the direction given by the
          ray.

        * A vertex and a line can generate an infinite edge. The
          position of the vertex on the line is arbitrary in this
          case, only its transverse position matters. The direction of
          the edge is given by the line generator.

        For example, take the half-plane::

            sage: half_plane = Polyhedron(ieqs=[(0,1,0)])
            sage: half_plane.Hrepresentation()
            (An inequality (1, 0) x + 0 >= 0,)

        Its (non-unique) V-representation consists of a vertex, a ray,
        and a line. The only edge is spanned by the vertex and the
        line generator, so they are adjacent::

            sage: half_plane.Vrepresentation()
            (A line in the direction (0, 1), A ray in the direction (1, 0), A vertex at (0, 0))
            sage: half_plane.vertex_adjacency_matrix()
            [0 0 1]
            [0 0 0]
            [1 0 0]

        In one dimension higher, that is for a half-space in 3
        dimensions, there is no one-dimensional face. Hence nothing is
        adjacent::

            sage: Polyhedron(ieqs=[(0,1,0,0)]).vertex_adjacency_matrix()
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]
            [0 0 0 0]

        EXAMPLES:

        In a bounded polygon, every vertex has precisely two adjacent ones::

            sage: P = Polyhedron(vertices=[(0, 1), (1, 0), (3, 0), (4, 1)])
            sage: for v in P.Vrep_generator():
            ....:     print("{} {}".format(P.adjacency_matrix().row(v.index()), v))
            (0, 1, 0, 1) A vertex at (0, 1)
            (1, 0, 1, 0) A vertex at (1, 0)
            (0, 1, 0, 1) A vertex at (3, 0)
            (1, 0, 1, 0) A vertex at (4, 1)

        If the V-representation of the polygon contains vertices and
        one ray, then each V-representation object is adjacent to two
        V-representation objects::

            sage: P = Polyhedron(vertices=[(0, 1), (1, 0), (3, 0), (4, 1)],
            ....:                rays=[(0,1)])
            sage: for v in P.Vrep_generator():
            ....:       print("{} {}".format(P.adjacency_matrix().row(v.index()), v))
            (0, 1, 0, 0, 1) A ray in the direction (0, 1)
            (1, 0, 1, 0, 0) A vertex at (0, 1)
            (0, 1, 0, 1, 0) A vertex at (1, 0)
            (0, 0, 1, 0, 1) A vertex at (3, 0)
            (1, 0, 0, 1, 0) A vertex at (4, 1)

        If the V-representation of the polygon contains vertices and
        two distinct rays, then each vertex is adjacent to two
        V-representation objects (which can now be vertices or
        rays). The two rays are not adjacent to each other::

            sage: P = Polyhedron(vertices=[(0, 1), (1, 0), (3, 0), (4, 1)],
            ....:                rays=[(0,1), (1,1)])
            sage: for v in P.Vrep_generator():
            ....:     print("{} {}".format(P.adjacency_matrix().row(v.index()), v))
            (0, 1, 0, 0, 0) A ray in the direction (0, 1)
            (1, 0, 1, 0, 0) A vertex at (0, 1)
            (0, 1, 0, 0, 1) A vertex at (1, 0)
            (0, 0, 0, 0, 1) A ray in the direction (1, 1)
            (0, 0, 1, 1, 0) A vertex at (3, 0)

        The vertex adjacency matrix has base ring integers. This way one can express various
        counting questions::

            sage: P = polytopes.cube()
            sage: Q = P.stack(P.faces(2)[0])
            sage: M = Q.vertex_adjacency_matrix()
            sage: sum(M)
            (4, 4, 3, 3, 4, 4, 4, 3, 3)
            sage: G = Q.vertex_graph()                                                  # needs sage.graphs
            sage: G.degree()                                                            # needs sage.graphs
            [4, 4, 3, 3, 4, 4, 4, 3, 3]

        TESTS:

        Check that :trac:`28828` is fixed::

                sage: P.adjacency_matrix().is_immutable()
                True
        """
        return self.combinatorial_polyhedron().vertex_adjacency_matrix(algorithm=algorithm)

    adjacency_matrix = vertex_adjacency_matrix

    @cached_method
    def facet_adjacency_matrix(self, algorithm=None):
        """
        Return the adjacency matrix for the facets.

        INPUT:

        - ``algorithm`` -- string (optional);
          specify whether the face generator starts with facets or vertices:

          * ``'primal'`` -- start with the facets
          * ``'dual'`` -- start with the vertices
          * ``None`` -- choose automatically

        EXAMPLES::

            sage: s4 = polytopes.simplex(4, project=True)
            sage: s4.facet_adjacency_matrix()
            [0 1 1 1 1]
            [1 0 1 1 1]
            [1 1 0 1 1]
            [1 1 1 0 1]
            [1 1 1 1 0]

            sage: p = Polyhedron(vertices=[(0,0),(1,0),(0,1)])
            sage: p.facet_adjacency_matrix()
            [0 1 1]
            [1 0 1]
            [1 1 0]


        The facet adjacency matrix has base ring integers. This way one can express various
        counting questions::

            sage: P = polytopes.cube()
            sage: Q = P.stack(P.faces(2)[0])
            sage: M = Q.facet_adjacency_matrix()
            sage: sum(M)
            (4, 4, 4, 4, 3, 3, 3, 3, 4)

        TESTS:

        Check that :trac:`28828` is fixed::

            sage: s4.facet_adjacency_matrix().is_immutable()
            True

        Checks that :trac:`22455` is fixed::

            sage: s = polytopes.simplex(2)
            sage: s.facet_adjacency_matrix()
            [0 1 1]
            [1 0 1]
            [1 1 0]
        """
        return self.combinatorial_polyhedron().facet_adjacency_matrix(algorithm=algorithm)

    def a_maximal_chain(self):
        r"""
        Return a maximal chain of the face lattice in increasing order.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: P.a_maximal_chain()
            [A -1-dimensional face of a Polyhedron in ZZ^3,
             A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices,
             A 3-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 8 vertices]
            sage: P = polytopes.cube()
            sage: chain = P.a_maximal_chain(); chain
            [A -1-dimensional face of a Polyhedron in ZZ^3,
             A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex,
             A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices,
             A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices,
             A 3-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 8 vertices]
            sage: [face.ambient_V_indices() for face in chain]
            [(), (5,), (0, 5), (0, 3, 4, 5), (0, 1, 2, 3, 4, 5, 6, 7)]

        TESTS:

        Check output for the empty polyhedron::

            sage: P = Polyhedron()
            sage: P.a_maximal_chain()
            [A -1-dimensional face of a Polyhedron in ZZ^0]
        """
        comb_chain = self.combinatorial_polyhedron().a_maximal_chain()

        from sage.geometry.polyhedron.face import combinatorial_face_to_polyhedral_face
        empty_face = self.faces(-1)[0]
        universe = self.faces(self.dim())[0]

        if self.dim() == -1:
            return [empty_face]

        return [empty_face] + \
               [combinatorial_face_to_polyhedral_face(self, face)
                for face in comb_chain] + \
               [universe]

    def is_simplex(self):
        r"""
        Return whether the polyhedron is a simplex.

        A simplex is a bounded polyhedron with `d+1` vertices, where
        `d` is the dimension.

        EXAMPLES::

            sage: Polyhedron([(0,0,0), (1,0,0), (0,1,0)]).is_simplex()
            True
            sage: polytopes.simplex(3).is_simplex()
            True
            sage: polytopes.hypercube(3).is_simplex()
            False
        """
        return self.is_compact() and (self.dim()+1 == self.n_vertices())

    def simplicity(self):
        r"""
        Return the largest integer `k` such that the polytope is `k`-simple.

        A polytope `P` is `k`-simple, if every `(d-1-k)`-face
        is contained in exactly `k+1` facets of `P` for `1 \leq k \leq d-1`.
        Equivalently it is `k`-simple if the polar/dual polytope is `k`-simplicial.
        If ``self`` is a simplex, it returns its dimension.

        EXAMPLES::

            sage: polytopes.hypersimplex(4,2).simplicity()
            1
            sage: polytopes.hypersimplex(5,2).simplicity()
            2
            sage: polytopes.hypersimplex(6,2).simplicity()
            3
            sage: polytopes.simplex(3).simplicity()
            3
            sage: polytopes.simplex(1).simplicity()
            1

        The method is not implemented for unbounded polyhedra::

            sage: p = Polyhedron(vertices=[(0,0)],rays=[(1,0),(0,1)])
            sage: p.simplicity()
            Traceback (most recent call last):
            ...
            NotImplementedError: this function is implemented for polytopes only
        """
        if not self.is_compact():
            raise NotImplementedError("this function is implemented for polytopes only")
        return self.combinatorial_polyhedron().simplicity()

    def is_simple(self):
        """
        Test for simplicity of a polytope.

        See :wikipedia:`Simple_polytope`

        EXAMPLES::

            sage: p = Polyhedron([[0,0,0],[1,0,0],[0,1,0],[0,0,1]])
            sage: p.is_simple()
            True
            sage: p = Polyhedron([[0,0,0],[4,4,0],[4,0,0],[0,4,0],[2,2,2]])
            sage: p.is_simple()
            False
        """
        if not self.is_compact():
            return False
        return self.combinatorial_polyhedron().is_simple()

    def simpliciality(self):
        r"""
        Return the largest integer `k` such that the polytope is `k`-simplicial.

        A polytope is `k`-simplicial, if every `k`-face is a simplex.
        If ``self`` is a simplex, returns its dimension.

        EXAMPLES::

            sage: polytopes.cyclic_polytope(10,4).simpliciality()
            3
            sage: polytopes.hypersimplex(5,2).simpliciality()
            2
            sage: polytopes.cross_polytope(4).simpliciality()
            3
            sage: polytopes.simplex(3).simpliciality()
            3
            sage: polytopes.simplex(1).simpliciality()
            1

        The method is not implemented for unbounded polyhedra::

            sage: p = Polyhedron(vertices=[(0,0)],rays=[(1,0),(0,1)])
            sage: p.simpliciality()
            Traceback (most recent call last):
            ...
            NotImplementedError: this function is implemented for polytopes only
        """
        if not self.is_compact():
            raise NotImplementedError("this function is implemented for polytopes only")
        return self.combinatorial_polyhedron().simpliciality()

    def is_simplicial(self):
        """
        Tests if the polytope is simplicial

        A polytope is simplicial if every facet is a simplex.

        See :wikipedia:`Simplicial_polytope`

        EXAMPLES::

            sage: p = polytopes.hypercube(3)
            sage: p.is_simplicial()
            False
            sage: q = polytopes.simplex(5, project=True)
            sage: q.is_simplicial()
            True
            sage: p = Polyhedron([[0,0,0],[1,0,0],[0,1,0],[0,0,1]])
            sage: p.is_simplicial()
            True
            sage: q = Polyhedron([[1,1,1],[-1,1,1],[1,-1,1],[-1,-1,1],[1,1,-1]])
            sage: q.is_simplicial()
            False
            sage: P = polytopes.simplex(); P
            A 3-dimensional polyhedron in ZZ^4 defined as the convex hull of 4 vertices
            sage: P.is_simplicial()
            True

        The method is not implemented for unbounded polyhedra::

            sage: p = Polyhedron(vertices=[(0,0)],rays=[(1,0),(0,1)])
            sage: p.is_simplicial()
            Traceback (most recent call last):
            ...
            NotImplementedError: this function is implemented for polytopes only
        """
        if not self.is_compact():
            raise NotImplementedError("this function is implemented for polytopes only")
        return self.combinatorial_polyhedron().is_simplicial()

    def is_pyramid(self, certificate=False):
        """
        Test whether the polytope is a pyramid over one of its facets.

        INPUT:

        - ``certificate`` -- boolean (default: ``False``); specifies whether
          to return a vertex of the polytope which is the apex of a pyramid,
          if found

        OUTPUT:

        If ``certificate`` is ``True``, returns a tuple containing:

        1. Boolean.
        2. The apex of the pyramid or ``None``.

        If ``certificate`` is ``False`` returns a boolean.

        EXAMPLES::

            sage: P = polytopes.simplex(3)
            sage: P.is_pyramid()
            True
            sage: P.is_pyramid(certificate=True)
            (True, A vertex at (1, 0, 0, 0))
            sage: egyptian_pyramid = polytopes.regular_polygon(4).pyramid()             # needs sage.rings.number_field
            sage: egyptian_pyramid.is_pyramid()                                         # needs sage.rings.number_field
            True
            sage: Q = polytopes.octahedron()
            sage: Q.is_pyramid()
            False

        For the `0`-dimensional polyhedron, the output is ``True``,
        but it cannot be constructed as a pyramid over the empty polyhedron::

            sage: P = Polyhedron([[0]])
            sage: P.is_pyramid()
            True
            sage: Polyhedron().pyramid()
            Traceback (most recent call last):
            ...
            ZeroDivisionError: rational division by zero
        """
        if not self.is_compact():
            raise ValueError("polyhedron has to be compact")

        return self.combinatorial_polyhedron().is_pyramid(certificate)

    def is_bipyramid(self, certificate=False):
        r"""
        Test whether the polytope is combinatorially equivalent to a
        bipyramid over some polytope.

        INPUT:

        - ``certificate`` -- boolean (default: ``False``); specifies whether
          to return two vertices of the polytope which are the apices of a
          bipyramid, if found

        OUTPUT:

        If ``certificate`` is ``True``, returns a tuple containing:

        1. Boolean.
        2. ``None`` or a tuple containing:
            a. The first apex.
            b. The second apex.

        If ``certificate`` is ``False`` returns a boolean.

        EXAMPLES::

            sage: P = polytopes.octahedron()
            sage: P.is_bipyramid()
            True
            sage: P.is_bipyramid(certificate=True)
            (True, [A vertex at (1, 0, 0), A vertex at (-1, 0, 0)])
            sage: Q = polytopes.cyclic_polytope(3,7)
            sage: Q.is_bipyramid()
            False
            sage: R = Q.bipyramid()
            sage: R.is_bipyramid(certificate=True)
            (True, [A vertex at (1, 3, 13, 63), A vertex at (-1, 3, 13, 63)])

        TESTS::

            sage: P = polytopes.permutahedron(4).bipyramid()
            sage: P.is_bipyramid()
            True

            sage: P = polytopes.cube()
            sage: P.is_bipyramid()
            False

            sage: P = Polyhedron(vertices=[[0,1], [1,0]], rays=[[1,1]])
            sage: P.is_bipyramid()
            Traceback (most recent call last):
            ...
            ValueError: polyhedron has to be compact
        """
        if not self.is_compact():
            raise ValueError("polyhedron has to be compact")

        return self.combinatorial_polyhedron().is_bipyramid(certificate)

    def is_prism(self, certificate=False):
        """
        Test whether the polytope is combinatorially equivalent to a prism of
        some polytope.

        INPUT:

        - ``certificate`` -- boolean (default: ``False``); specifies whether
          to return two facets of the polytope which are the bases of a prism,
          if found

        OUTPUT:

        If ``certificate`` is ``True``, returns a tuple containing:

        1. Boolean.
        2. ``None`` or a tuple containing:
            a. List of the vertices of the first base facet.
            b. List of the vertices of the second base facet.

        If ``certificate`` is ``False`` returns a boolean.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: P.is_prism()
            True
            sage: P.is_prism(certificate=True)
            (True,
             [(A vertex at (1, -1, -1),
               A vertex at (1, -1, 1),
               A vertex at (-1, -1, 1),
               A vertex at (-1, -1, -1)),
              (A vertex at (1, 1, -1),
               A vertex at (1, 1, 1),
               A vertex at (-1, 1, -1),
               A vertex at (-1, 1, 1))])
            sage: Q = polytopes.cyclic_polytope(3,8)
            sage: Q.is_prism()
            False
            sage: R = Q.prism()
            sage: R.is_prism(certificate=True)
            (True,
             [(A vertex at (0, 3, 9, 27),
               A vertex at (0, 6, 36, 216),
               A vertex at (0, 0, 0, 0),
               A vertex at (0, 7, 49, 343),
               A vertex at (0, 5, 25, 125),
               A vertex at (0, 1, 1, 1),
               A vertex at (0, 2, 4, 8),
               A vertex at (0, 4, 16, 64)),
              (A vertex at (1, 6, 36, 216),
               A vertex at (1, 0, 0, 0),
               A vertex at (1, 7, 49, 343),
               A vertex at (1, 5, 25, 125),
               A vertex at (1, 1, 1, 1),
               A vertex at (1, 2, 4, 8),
               A vertex at (1, 4, 16, 64),
               A vertex at (1, 3, 9, 27))])

        TESTS::

            sage: P = polytopes.cross_polytope(5)
            sage: P.is_prism()
            False

            sage: P = polytopes.permutahedron(4).prism()
            sage: P.is_prism()
            True

            sage: P = Polyhedron(vertices=[[0,1], [1,0]], rays=[[1,1]])
            sage: P.is_prism()
            Traceback (most recent call last):
            ...
            NotImplementedError: polyhedron has to be compact
        """
        if not self.is_compact():
            raise NotImplementedError("polyhedron has to be compact")

        return self.combinatorial_polyhedron().is_prism(certificate)

    def is_lawrence_polytope(self):
        """
        Return ``True`` if ``self`` is a Lawrence polytope.

        A polytope is called a Lawrence polytope if it has a centrally
        symmetric (normalized) Gale diagram.

        EXAMPLES::

            sage: P = polytopes.hypersimplex(5,2)
            sage: L = P.lawrence_polytope()
            sage: L.is_lattice_polytope()
            True

            sage: egyptian_pyramid = polytopes.regular_polygon(4).pyramid()             # needs sage.number_field
            sage: egyptian_pyramid.is_lawrence_polytope()                               # needs sage.number_field
            True

            sage: polytopes.octahedron().is_lawrence_polytope()
            False

        REFERENCES:

            For more information, see [BaSt1990]_.
        """
        if not self.is_compact():
            raise NotImplementedError("self must be a polytope")

        return self.combinatorial_polyhedron().is_lawrence_polytope()

    def neighborliness(self):
        r"""
        Return the largest ``k``, such that the polyhedron is ``k``-neighborly.

        A polyhedron is `k`-neighborly if every set of `n` vertices forms a face
        for `n` up to `k`.

        In case of the `d`-dimensional simplex, it returns `d + 1`.

        .. SEEALSO::

            :meth:`is_neighborly`

        EXAMPLES::

            sage: cube = polytopes.cube()
            sage: cube.neighborliness()
            1
            sage: P = Polyhedron(); P
            The empty polyhedron in ZZ^0
            sage: P.neighborliness()
            0
            sage: P = Polyhedron([[0]]); P
            A 0-dimensional polyhedron in ZZ^1 defined as the convex hull of 1 vertex
            sage: P.neighborliness()
            1
            sage: S = polytopes.simplex(5); S
            A 5-dimensional polyhedron in ZZ^6 defined as the convex hull of 6 vertices
            sage: S.neighborliness()
            6
            sage: C = polytopes.cyclic_polytope(7,10); C
            A 7-dimensional polyhedron in QQ^7 defined as the convex hull of 10 vertices
            sage: C.neighborliness()
            3
            sage: C = polytopes.cyclic_polytope(6,11); C
            A 6-dimensional polyhedron in QQ^6 defined as the convex hull of 11 vertices
            sage: C.neighborliness()
            3
            sage: [polytopes.cyclic_polytope(5,n).neighborliness() for n in range(6,10)]
            [6, 2, 2, 2]

        """
        return self.combinatorial_polyhedron().neighborliness()

    def is_neighborly(self, k=None):
        r"""
        Return whether the polyhedron is neighborly.

        If the input ``k`` is provided, then return whether the polyhedron is ``k``-neighborly

        A polyhedron is neighborly if every set of `n` vertices forms a face
        for `n` up to floor of half the dimension of the polyhedron.
        It is `k`-neighborly if this is true for `n` up to `k`.

        INPUT:

        - ``k`` -- the dimension up to which to check if every set of ``k``
          vertices forms a face. If no ``k`` is provided, check up to floor
          of half the dimension of the polyhedron.

        OUTPUT:

        - ``True`` if every set of up to ``k`` vertices forms a face,
        - ``False`` otherwise

        .. SEEALSO::

            :meth:`neighborliness`

        EXAMPLES::

            sage: cube = polytopes.hypercube(3)
            sage: cube.is_neighborly()
            True
            sage: cube = polytopes.hypercube(4)
            sage: cube.is_neighborly()
            False

        Cyclic polytopes are neighborly::

            sage: all(polytopes.cyclic_polytope(i, i + 1 + j).is_neighborly() for i in range(5) for j in range(3))
            True

        The neighborliness of a polyhedron equals floor of dimension half
        (or larger in case of a simplex) if and only if the polyhedron
        is neighborly::

            sage: testpolys = [polytopes.cube(), polytopes.cyclic_polytope(6, 9), polytopes.simplex(6)]
            sage: [(P.neighborliness() >= P.dim() // 2) == P.is_neighborly() for P in testpolys]
            [True, True, True]

        """
        return self.combinatorial_polyhedron().is_neighborly()

    def join_of_Vrep(self, *Vrepresentatives):
        r"""
        Return the smallest face that contains ``Vrepresentatives``.

        INPUT:

        - ``Vrepresentatives`` -- vertices/rays/lines of ``self`` or indices of such

        OUTPUT: a :class:`~sage.geometry.polyhedron.face.PolyhedronFace`

        .. NOTE::

            In the case of unbounded polyhedra, the join of rays etc. may not be well-defined.

        EXAMPLES::

            sage: P = polytopes.permutahedron(5)
            sage: P.join_of_Vrep(1)
            A 0-dimensional face of a Polyhedron in ZZ^5 defined as the convex hull of 1 vertex
            sage: P.join_of_Vrep()
            A -1-dimensional face of a Polyhedron in ZZ^5
            sage: P.join_of_Vrep(0,12,13).ambient_V_indices()
            (0, 12, 13, 68)

        The input is flexible::

            sage: P.join_of_Vrep(2, P.vertices()[3], P.Vrepresentation(4))
            A 2-dimensional face of a Polyhedron in ZZ^5 defined as the convex hull of 6 vertices

        ::

            sage: P = polytopes.cube()
            sage: a, b = P.faces(0)[:2]
            sage: P.join_of_Vrep(a, b)
            A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices

        In the case of an unbounded polyhedron, the join may not be well-defined::

            sage: P = Polyhedron(vertices=[[1,0], [0,1]], rays=[[1,1]])
            sage: P.join_of_Vrep(0)
            A 0-dimensional face of a Polyhedron in QQ^2 defined as the convex hull of 1 vertex
            sage: P.join_of_Vrep(0,1)
            A 1-dimensional face of a Polyhedron in QQ^2 defined as the convex hull of 2 vertices
            sage: P.join_of_Vrep(0,2)
            A 1-dimensional face of a Polyhedron in QQ^2 defined as the convex hull of 1 vertex and 1 ray
            sage: P.join_of_Vrep(1,2)
            A 1-dimensional face of a Polyhedron in QQ^2 defined as the convex hull of 1 vertex and 1 ray
            sage: P.join_of_Vrep(2)
            Traceback (most recent call last):
            ...
            ValueError: the join is not well-defined

        The ``Vrepresentatives`` must be of ``self``::

            sage: P = polytopes.cube(backend='ppl')
            sage: Q = polytopes.cube(backend='field')
            sage: v = P.vertices()[0]
            sage: P.join_of_Vrep(v)
            A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex
            sage: Q.join_of_Vrep(v)
            Traceback (most recent call last):
            ...
            ValueError: not a Vrepresentative of ``self``
            sage: f = P.faces(0)[0]
            sage: P.join_of_Vrep(v)
            A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex
            sage: Q.join_of_Vrep(v)
            Traceback (most recent call last):
            ...
            ValueError: not a Vrepresentative of ``self``

        TESTS:

        ``least_common_superface_of_Vrep`` is an alias::

            sage: P.least_common_superface_of_Vrep(v)
            A 0-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 1 vertex
            sage: P.least_common_superface_of_Vrep == P.join_of_Vrep
            True

        Error message for invalid input::

            sage: P.join_of_Vrep('foo')
            Traceback (most recent call last):
            ...
            ValueError: foo is not a Vrepresentative
        """
        from sage.geometry.polyhedron.representation import Vrepresentation
        from sage.geometry.polyhedron.face import PolyhedronFace

        new_indices = [0]*len(Vrepresentatives)
        for i, v in enumerate(Vrepresentatives):
            if isinstance(v, PolyhedronFace) and v.dim() == 0:
                if v.polyhedron() is not self:
                    raise ValueError("not a Vrepresentative of ``self``")
                new_indices[i] = v.ambient_V_indices()[0]
            elif v in ZZ:
                new_indices[i] = v
            elif isinstance(v, Vrepresentation):
                if v.polyhedron() is not self:
                    raise ValueError("not a Vrepresentative of ``self``")
                new_indices[i] = v.index()
            else:
                raise ValueError("{} is not a Vrepresentative".format(v))

        return self.face_generator().join_of_Vrep(*new_indices)

    least_common_superface_of_Vrep = join_of_Vrep

    def meet_of_Hrep(self, *Hrepresentatives):
        r"""
        Return the largest face that is contained in ``Hrepresentatives``.

        INPUT:

        - ``Hrepresentatives`` -- facets or indices of Hrepresentatives;
          the indices are assumed to be the indices of the
          :meth:`~sage.geometry.polyhedron.base0.Polyhedron_base0.Hrepresentation`

        OUTPUT: a :class:`~sage.geometry.polyhedron.face.PolyhedronFace`

        EXAMPLES::

            sage: P = polytopes.permutahedron(5)
            sage: P.meet_of_Hrep()
            A 4-dimensional face of a Polyhedron in ZZ^5 defined as the convex hull of 120 vertices
            sage: P.meet_of_Hrep(1)
            A 3-dimensional face of a Polyhedron in ZZ^5 defined as the convex hull of 24 vertices
            sage: P.meet_of_Hrep(4)
            A 3-dimensional face of a Polyhedron in ZZ^5 defined as the convex hull of 12 vertices
            sage: P.meet_of_Hrep(1,3,7)
            A 1-dimensional face of a Polyhedron in ZZ^5 defined as the convex hull of 2 vertices
            sage: P.meet_of_Hrep(1,3,7).ambient_H_indices()
            (0, 1, 3, 7)

        The indices are the indices of the
        :meth:`~sage.geometry.polyhedron.base0.Polyhedron_base0.Hrepresentation`.
        ``0`` corresponds to an equation and is ignored::

            sage: P.meet_of_Hrep(0)
            A 4-dimensional face of a Polyhedron in ZZ^5 defined as the convex hull of 120 vertices

        The input is flexible::

            sage: P.meet_of_Hrep(P.facets()[-1], P.inequalities()[2], 7)
            A 1-dimensional face of a Polyhedron in ZZ^5 defined as the convex hull of 2 vertices

        The ``Hrepresentatives`` must belong to ``self``::

            sage: P = polytopes.cube(backend='ppl')
            sage: Q = polytopes.cube(backend='field')
            sage: f = P.facets()[0]
            sage: P.meet_of_Hrep(f)
            A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices
            sage: Q.meet_of_Hrep(f)
            Traceback (most recent call last):
            ...
            ValueError: not a facet of ``self``
            sage: f = P.inequalities()[0]
            sage: P.meet_of_Hrep(f)
            A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 4 vertices
            sage: Q.meet_of_Hrep(f)
            Traceback (most recent call last):
            ...
            ValueError: not a facet of ``self``

        TESTS:

        Equations are not considered by the combinatorial polyhedron.
        We check that the index corresponds to the
        :meth:`~sage.geometry.polyhedron.base0.Polyhedron_base0.Hrepresentation` index::

            sage: P = polytopes.permutahedron(3, backend='field')
            sage: P.Hrepresentation()
            (An inequality (0, 0, 1) x - 1 >= 0,
             An inequality (0, 1, 0) x - 1 >= 0,
             An inequality (0, 1, 1) x - 3 >= 0,
             An inequality (1, 0, 0) x - 1 >= 0,
             An inequality (1, 0, 1) x - 3 >= 0,
             An inequality (1, 1, 0) x - 3 >= 0,
             An equation (1, 1, 1) x - 6 == 0)
            sage: P.meet_of_Hrep(0).ambient_Hrepresentation()
            (An inequality (0, 0, 1) x - 1 >= 0, An equation (1, 1, 1) x - 6 == 0)

            sage: P = polytopes.permutahedron(3, backend='ppl')
            sage: P.Hrepresentation()
            (An equation (1, 1, 1) x - 6 == 0,
             An inequality (1, 1, 0) x - 3 >= 0,
             An inequality (-1, -1, 0) x + 5 >= 0,
             An inequality (0, 1, 0) x - 1 >= 0,
             An inequality (-1, 0, 0) x + 3 >= 0,
             An inequality (1, 0, 0) x - 1 >= 0,
             An inequality (0, -1, 0) x + 3 >= 0)
            sage: P.meet_of_Hrep(1).ambient_Hrepresentation()
            (An equation (1, 1, 1) x - 6 == 0, An inequality (1, 1, 0) x - 3 >= 0)

        ``greatest_common_subface_of_Hrep`` is an alias::

            sage: P.greatest_common_subface_of_Hrep(1).ambient_Hrepresentation()
            (An equation (1, 1, 1) x - 6 == 0, An inequality (1, 1, 0) x - 3 >= 0)
            sage: P.greatest_common_subface_of_Hrep == P.meet_of_Hrep
            True

        Error message for invalid input::

            sage: P.meet_of_Hrep('foo')
            Traceback (most recent call last):
            ...
            ValueError: foo is not a Hrepresentative
        """
        from sage.geometry.polyhedron.representation import Inequality, Equation
        from sage.geometry.polyhedron.face import PolyhedronFace

        # Equations are ignored by combinatorial polyhedron for indexing.
        offset = 0
        if self.n_equations() and self.Hrepresentation(0).is_equation():
            offset = self.n_equations()

        new_indices = []
        for i, facet in enumerate(Hrepresentatives):
            if isinstance(facet, PolyhedronFace) and facet.dim() + 1 == self.dim():
                if facet.polyhedron() is not self:
                    raise ValueError("not a facet of ``self``")
                H_indices = facet.ambient_H_indices()
                facet = H_indices[0] if H_indices[0] >= offset else H_indices[-1]

            if facet in ZZ and facet >= offset:
                # Note that ``CombinatorialPolyhedron`` ignores indices of equations
                # and has equations last.
                new_indices.append(facet - offset)
            elif isinstance(facet, Inequality):
                if facet.polyhedron() is not self:
                    raise ValueError("not a facet of ``self``")
                new_indices.append(facet.index() - offset)
            elif isinstance(facet, Equation):
                # Ignore equations.
                continue
            elif facet in ZZ and 0 <= facet < offset:
                # Ignore equations.
                continue
            else:
                raise ValueError("{} is not a Hrepresentative".format(facet))

        return self.face_generator().meet_of_Hrep(*new_indices)

    greatest_common_subface_of_Hrep = meet_of_Hrep

    def _test_combinatorial_face_as_combinatorial_polyhedron(self, tester=None, **options):
        """
        Run tests on obtaining the combinatorial face as combinatorial polyhedron.

        TESTS::

            sage: polytopes.cross_polytope(3)._test_combinatorial_face_as_combinatorial_polyhedron()
        """
        if not self.is_compact():
            return
        if self.dim() > 7 or self.n_vertices() > 100 or self.n_facets() > 100:
            # Avoid very long tests.
            return
        if self.dim() < 1:
            # Avoid trivial cases.
            return

        from sage.misc.prandom import random

        if tester is None:
            tester = self._tester(**options)

        it = self.face_generator()
        _ = next(it), next(it)  # get rid of non_proper faces
        C1 = self.combinatorial_polyhedron()
        it1 = C1.face_iter()
        C2 = C1.dual()
        algorithm = 'primal' if it1.dual else 'dual'
        it2 = C2.face_iter(algorithm=algorithm)

        for f in it:
            f1 = next(it1)
            f2 = next(it2)
            if random() < 0.95:
                # Only test a random 5 percent of the faces.
                continue

            P = f.as_polyhedron()
            D1 = f1.as_combinatorial_polyhedron()
            D2 = f2.as_combinatorial_polyhedron(quotient=True).dual()
            D1._test_bitsets(tester, **options)
            D2._test_bitsets(tester, **options)
            try:
                import sage.graphs.graph
                assert sage.graphs.graph  # to muffle pyflakes
            except ImportError:
                pass
            else:
                tester.assertTrue(P.combinatorial_polyhedron().vertex_facet_graph().is_isomorphic(D1.vertex_facet_graph()))
                tester.assertTrue(P.combinatorial_polyhedron().vertex_facet_graph().is_isomorphic(D2.vertex_facet_graph()))
