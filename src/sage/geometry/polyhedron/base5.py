r"""
Base class for polyhedra: Methods for constructing new polyhedra

Except for affine hull and affine hull projection.
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

from sage.structure.element import coerce_binop, Vector, Matrix

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.matrix.constructor import matrix

from sage.modules.free_module_element import vector

from .base4 import Polyhedron_base4


class Polyhedron_base5(Polyhedron_base4):
    """
    Methods constructing new polyhedra
    except for affine hull and affine hull projection.

    See :class:`sage.geometry.polyhedron.base.Polyhedron_base`.

    TESTS::

        sage: from sage.geometry.polyhedron.base5 import Polyhedron_base5
        sage: P = polytopes.cube()
        sage: Q = polytopes.cross_polytope(3)

    Unary operations::

        sage: Polyhedron_base5.__neg__(P)
        A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 8 vertices
        sage: Polyhedron_base5.polar(P)
        A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 6 vertices
        sage: Polyhedron_base5.pyramid(P)
        A 4-dimensional polyhedron in QQ^4 defined as the convex hull of 9 vertices
        sage: Polyhedron_base5.bipyramid(P)
        A 4-dimensional polyhedron in QQ^4 defined as the convex hull of 10 vertices
        sage: Polyhedron_base5.prism(P)
        A 4-dimensional polyhedron in ZZ^4 defined as the convex hull of 16 vertices
        sage: Polyhedron_base5.truncation(P)
        A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 24 vertices
        sage: Polyhedron_base5.lawrence_polytope(P)
        A 11-dimensional polyhedron in ZZ^11 defined as the convex hull of 16 vertices

    Binary operations::

        sage: Polyhedron_base5.minkowski_sum(P, Q)
        A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 24 vertices
        sage: Polyhedron_base5.minkowski_difference(P, Q)
        A 0-dimensional polyhedron in QQ^3 defined as the convex hull of 1 vertex
        sage: Polyhedron_base5.product(P, Q)
        A 6-dimensional polyhedron in ZZ^6 defined as the convex hull of 48 vertices
        sage: Polyhedron_base5.join(P, Q)
        A 7-dimensional polyhedron in ZZ^7 defined as the convex hull of 14 vertices
        sage: Polyhedron_base5.subdirect_sum(P, Q)
        A 6-dimensional polyhedron in ZZ^6 defined as the convex hull of 14 vertices
        sage: Polyhedron_base5.direct_sum(P, Q)
        A 6-dimensional polyhedron in QQ^6 defined as the convex hull of 14 vertices
        sage: Polyhedron_base5.convex_hull(P, Q)
        A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 8 vertices
        sage: Polyhedron_base5.intersection(P, Q)
        A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 6 vertices

    Actions::

        sage: Polyhedron_base5.translation(P, vector([1, 1, 1]))
        A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 8 vertices
        sage: Polyhedron_base5.dilation(P, 2)
        A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 8 vertices
        sage: Polyhedron_base5.__truediv__(P, 2)
        A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 8 vertices
        sage: Polyhedron_base5.linear_transformation(P, matrix([[2, 1, 2], [1, 2, 4], [3, 1, 2]]))
        A 2-dimensional polyhedron in ZZ^3 defined as the convex hull of 4 vertices

    Methods using a face::

        sage: Polyhedron_base5.face_truncation(P, P.facets()[0])
        A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 8 vertices
        sage: Polyhedron_base5.stack(P, P.facets()[0])
        A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 9 vertices
        sage: Polyhedron_base5.wedge(P, P.facets()[0])
        A 4-dimensional polyhedron in QQ^4 defined as the convex hull of 12 vertices
        sage: Polyhedron_base5.face_split(P, P.faces(2)[0])
        A 4-dimensional polyhedron in QQ^4 defined as the convex hull of 10 vertices

    Methods using a vertex or vector::

        sage: Polyhedron_base5.lawrence_extension(P, P.vertices()[0])
        A 4-dimensional polyhedron in ZZ^4 defined as the convex hull of 9 vertices
        sage: Polyhedron_base5.one_point_suspension(P, P.vertices()[0])
        A 4-dimensional polyhedron in QQ^4 defined as the convex hull of 9 vertices
    """

    ###########################################################
    # Unary operations.
    ###########################################################

    def __neg__(self):
        """
        Negation of a polytope is defined as inverting the coordinates.

        EXAMPLES::

            sage: t = polytopes.simplex(3,project=False);  t.vertices()
            (A vertex at (0, 0, 0, 1), A vertex at (0, 0, 1, 0),
             A vertex at (0, 1, 0, 0), A vertex at (1, 0, 0, 0))
            sage: neg_ = -t
            sage: neg_.vertices()
            (A vertex at (-1, 0, 0, 0), A vertex at (0, -1, 0, 0),
             A vertex at (0, 0, -1, 0), A vertex at (0, 0, 0, -1))

        TESTS::

            sage: p = Polyhedron(ieqs=[[1,1,0]])
            sage: p.rays()
            (A ray in the direction (1, 0),)
            sage: pneg = p.__neg__()
            sage: pneg.rays()
            (A ray in the direction (-1, 0),)
        """
        return self.dilation(-1)

    def polar(self, in_affine_span=False):
        """
        Return the polar (dual) polytope.

        The original vertices are translated so that their barycenter
        is at the origin, and then the vertices are used as the
        coefficients in the polar inequalities.

        The polytope must be full-dimensional, unless ``in_affine_span`` is ``True``.
        If ``in_affine_span`` is ``True``, then the operation will be performed in the
        linear/affine span of the polyhedron (after translation).

        EXAMPLES::

            sage: p = Polyhedron(vertices=[[0,0,1], [0,1,0], [1,0,0], [0,0,0], [1,1,1]],
            ....:                base_ring=QQ); p
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 5 vertices
            sage: p.polar()
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 6 vertices

            sage: cube = polytopes.hypercube(3)
            sage: octahedron = polytopes.cross_polytope(3)
            sage: cube_dual = cube.polar()
            sage: octahedron == cube_dual
            True

        ``in_affine_span`` somewhat ignores equations, performing the polar in the
        spanned subspace (after translating barycenter to origin)::

            sage: P = polytopes.simplex(3, base_ring=QQ)
            sage: P.polar(in_affine_span=True)
            A 3-dimensional polyhedron in QQ^4 defined as the convex hull of 4 vertices

        Embedding the polytope in a higher dimension, commutes with polar in this case::

            sage: point = Polyhedron([[0]])
            sage: P = polytopes.cube().change_ring(QQ)
            sage: (P*point).polar(in_affine_span=True) == P.polar()*point
            True

        TESTS:

        Check that :issue:`25081` is fixed::

            sage: C = polytopes.hypercube(4,backend='cdd')
            sage: C.polar().backend()
            'cdd'

        Check that :issue:`28850` is fixed::

            sage: P = polytopes.simplex(3, base_ring=QQ)
            sage: P.polar()
            Traceback (most recent call last):
            ...
            ValueError: not full-dimensional; try with 'in_affine_span=True'

        Check that the double description is set up correctly::

            sage: P = Polyhedron([[1,0],[0,1],[-1,-1]], backend='field')
            sage: Q = P.change_ring(QQ, backend='ppl')
            sage: P.polar() == Q.polar()
            True

            sage: P = polytopes.simplex(4, backend='field')
            sage: Q = P.change_ring(QQ, backend='ppl')
            sage: P.polar(in_affine_span=True) == Q.polar(in_affine_span=True)
            True

        Check that it works, even when the equations are not orthogonal to each other::

            sage: P = polytopes.cube()*Polyhedron([[0,0,0]])
            sage: P = P.change_ring(QQ)

            sage: from sage.geometry.polyhedron.backend_field import Polyhedron_field
            sage: from sage.geometry.polyhedron.parent import Polyhedra_field
            sage: parent = Polyhedra_field(QQ, 6, 'field')
            sage: equations = [[0, 0, 0, 0, 1, 1, 1], [0, 0, 0, 0, -1, 1, -1], [0, 0, 0, 0, 1, -1, -1]]
            sage: Q = Polyhedron_field(parent, [P.vertices(), [], []], [P.inequalities(), equations],
            ....:                      Vrep_minimal=True, Hrep_minimal=True)
            sage: Q == P
            True
            sage: Q.polar(in_affine_span=True) == P.polar(in_affine_span=True)
            True
        """
        if not self.is_compact():
            raise ValueError("not a polytope")
        if not in_affine_span and not self.dim() == self.ambient_dim():
            raise ValueError("not full-dimensional; try with 'in_affine_span=True'")

        t_Vrep, t_Hrep, parent = self._translation_double_description(-self.center())
        t_verts = t_Vrep[0]
        t_ieqs = t_Hrep[0]
        t_eqns = t_Hrep[1]

        new_ieqs = ((1,) + tuple(-v) for v in t_verts)
        if self.n_vertices() == 1:
            new_verts = self.vertices()
        elif not self.n_equations():
            new_verts = ((-h/h[0])[1:] for h in t_ieqs)
        else:
            # Transform the equations such that the normals are pairwise orthogonal.
            t_eqns = list(t_eqns)
            for i, h in enumerate(t_eqns):
                for h1 in t_eqns[:i]:
                    a = h[1:]*h1[1:]
                    if a:
                        b = h1[1:]*h1[1:]
                        t_eqns[i] = b*h - a*h1

            def move_vertex_to_subspace(vertex):
                for h in t_eqns:
                    offset = vertex*h[1:]+h[0]
                    vertex = vertex-h[1:]*offset/(h[1:]*h[1:])
                return vertex

            new_verts = (move_vertex_to_subspace((-h/h[0])[1:]) for h in t_ieqs)

        pref_rep = 'Hrep' if self.n_vertices() <= self.n_inequalities() else 'Vrep'

        return parent.element_class(parent, [new_verts, [], []],
                                    [new_ieqs, t_eqns],
                                    Vrep_minimal=True, Hrep_minimal=True, pref_rep=pref_rep)

    def pyramid(self):
        """
        Return a polyhedron that is a pyramid over the original.

        EXAMPLES::

            sage: square = polytopes.hypercube(2);  square
            A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 4 vertices
            sage: egyptian_pyramid = square.pyramid();  egyptian_pyramid
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 5 vertices
            sage: egyptian_pyramid.n_vertices()
            5
            sage: for v in egyptian_pyramid.vertex_generator(): print(v)
            A vertex at (0, -1, -1)
            A vertex at (0, -1, 1)
            A vertex at (0, 1, -1)
            A vertex at (0, 1, 1)
            A vertex at (1, 0, 0)

        TESTS::

            sage: polytopes.simplex(backend='cdd').pyramid().backend()
            'cdd'
        """
        assert self.is_compact(), "Not a polytope."
        c = self.center()

        from itertools import chain
        new_verts = chain(([0] + x for x in self.Vrep_generator()),
                          [[1] + list(c)])
        new_ieqs = chain(([i.b()] + [-c*i.A() - i.b()] + list(i.A()) for i in self.inequalities()),
                         [[0, 1] + [0]*self.ambient_dim()])
        new_eqns = ([e.b()] + [0] + list(e.A()) for e in self.equations())

        pref_rep = 'Hrep' if self.n_vertices() > self.n_inequalities() else 'Vrep'
        parent = self.parent().base_extend(self.center().parent(), ambient_dim=self.ambient_dim()+1)

        if self.n_vertices() == 1:
            # Fix the polyhedron with one vertex.
            return parent.element_class(parent, [new_verts, [], []], None)

        return parent.element_class(parent, [new_verts, [], []],
                                    [new_ieqs, new_eqns],
                                    Vrep_minimal=True, Hrep_minimal=True, pref_rep=pref_rep)

    def _test_pyramid(self, tester=None, **options):
        """
        Run tests on the methods related to pyramids.

        TESTS:

            sage: polytopes.regular_polygon(4)._test_pyramid()                          # needs sage.rings.number_field
        """
        if tester is None:
            tester = self._tester(**options)

        def check_pyramid_certificate(P, cert):
            others = set(v for v in P.vertices() if not v == cert)
            if others:
                tester.assertTrue(any(set(f.ambient_Vrepresentation()) == others for f in P.facets()))

        if self.is_compact():
            b, cert = self.is_pyramid(certificate=True)
            if b:
                check_pyramid_certificate(self, cert)

            if 1 < self.n_vertices() < 50 and self.n_facets() < 50:
                pyr = self.pyramid()
                b, cert = pyr.is_pyramid(certificate=True)
                tester.assertTrue(b)
                check_pyramid_certificate(pyr, cert)
        else:
            with tester.assertRaises(AssertionError):
                pyr = self.pyramid()

        if self.is_compact() and 1 < self.n_vertices() < 50 and self.n_facets() < 50:
            # Check the pyramid of the polar.
            self_fraction_field = self.base_extend(QQ)

            polar = self_fraction_field.polar(in_affine_span=True)
            pyr_polar = polar.pyramid()
            b, cert = pyr_polar.is_pyramid(certificate=True)
            tester.assertTrue(b)
            check_pyramid_certificate(pyr_polar, cert)

            pyr = self_fraction_field.pyramid()
            polar_pyr = pyr.polar(in_affine_span=True)
            b, cert = polar_pyr.is_pyramid(certificate=True)
            tester.assertTrue(b)
            check_pyramid_certificate(polar_pyr, cert)

            try:
                import sage.graphs.graph
            except ImportError:
                pass
            else:
                tester.assertTrue(pyr_polar.is_combinatorially_isomorphic(pyr_polar))

            # Basic properties of the pyramid.

            # Check that the prism preserves the backend.
            tester.assertEqual(pyr.backend(), self.backend())

            tester.assertEqual(1 + self.n_vertices(), pyr.n_vertices())
            tester.assertEqual(self.n_equations(), pyr.n_equations())
            tester.assertEqual(1 + self.n_inequalities(), pyr.n_inequalities())

            if self.n_vertices() < 15 and self.n_facets() < 15:
                pyr._test_basic_properties()

    def bipyramid(self):
        """
        Return a polyhedron that is a bipyramid over the original.

        EXAMPLES::

            sage: octahedron = polytopes.cross_polytope(3)
            sage: cross_poly_4d = octahedron.bipyramid()
            sage: cross_poly_4d.n_vertices()
            8
            sage: q = [list(v) for v in cross_poly_4d.vertex_generator()]; q
            [[-1, 0, 0, 0],
             [0, -1, 0, 0],
             [0, 0, -1, 0],
             [0, 0, 0, -1],
             [0, 0, 0, 1],
             [0, 0, 1, 0],
             [0, 1, 0, 0],
             [1, 0, 0, 0]]

        Now check that bipyramids of cross-polytopes are cross-polytopes::

            sage: q2 = [list(v) for v in polytopes.cross_polytope(4).vertex_generator()]
            sage: [v in q2 for v in q]
            [True, True, True, True, True, True, True, True]

        TESTS::

            sage: polytopes.simplex(backend='cdd').bipyramid().backend()
            'cdd'
        """
        c = self.center()
        from itertools import chain
        new_verts = chain(([0] + list(x) for x in self.vertex_generator()),
                          [[1] + list(c), [-1] + list(c)])
        new_rays = ([0] + r for r in self.rays())
        new_lines = ([0] + l for l in self.lines())
        new_ieqs = chain(([i.b()] + [ c*i.A() + i.b()] + list(i.A()) for i in self.inequalities()),
                         ([i.b()] + [-c*i.A() - i.b()] + list(i.A()) for i in self.inequalities()))
        new_eqns = ([e.b()] + [0] + list(e.A()) for e in self.equations())

        pref_rep = 'Hrep' if 2 + (self.n_vertices() + self.n_rays()) >= 2*self.n_inequalities() else 'Vrep'
        parent = self.parent().base_extend(self.center().parent(), ambient_dim=self.ambient_dim()+1)

        if c not in self.relative_interior():
            # Fix polyhedra with non-proper center.
            return parent.element_class(parent, [new_verts, new_rays, new_lines], None)

        return parent.element_class(parent, [new_verts, new_rays, new_lines],
                                    [new_ieqs, new_eqns],
                                    Vrep_minimal=True, Hrep_minimal=True, pref_rep=pref_rep)

    def _test_bipyramid(self, tester=None, **options):
        """
        Run tests on the method :meth:`.bipyramid`.

        TESTS::

            sage: polytopes.cross_polytope(3)._test_bipyramid()
        """
        if tester is None:
            tester = self._tester(**options)

        if (self.n_vertices() + self.n_rays() >= 40
                or self.n_facets() >= 40
                or self.n_vertices() <= 1):
            return

        bipyramid = self.bipyramid()

        # Check that the double description is set up correctly.
        if self.base_ring().is_exact() and self.n_vertices() + self.n_rays() < 15 and self.n_facets() < 15:
            bipyramid._test_basic_properties(tester)

        # Check that the bipyramid preserves the backend.
        tester.assertEqual(bipyramid.backend(), self.backend())

        if self.center() not in self.relative_interior():
            # In this case (unbounded) the bipyramid behaves a bit different.
            return

        tester.assertEqual(2 + self.n_vertices(), bipyramid.n_vertices())
        tester.assertEqual(self.n_rays(), bipyramid.n_rays())
        tester.assertEqual(self.n_lines(), bipyramid.n_lines())
        tester.assertEqual(self.n_equations(), bipyramid.n_equations())
        tester.assertEqual(2*self.n_inequalities(), bipyramid.n_inequalities())

        if not self.is_compact():
            # ``is_bipyramid`` is only implemented for compact polyhedra.
            return

        b, cert = bipyramid.is_bipyramid(certificate=True)
        tester.assertTrue(b)

        if not self.is_bipyramid() and self.base_ring().is_exact():
            # In this case the certificate is unique.

            R = self.base_ring()
            a = (R(1),) + tuple(self.center())
            b = (R(-1),) + tuple(self.center())
            c, d = (tuple(v) for v in cert)
            tester.assertEqual(sorted([a, b]), sorted([c, d]))

    def prism(self):
        """
        Return a prism of the original polyhedron.

        EXAMPLES::

            sage: square = polytopes.hypercube(2)
            sage: cube = square.prism()
            sage: cube
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 8 vertices
            sage: hypercube = cube.prism()
            sage: hypercube.n_vertices()
            16

        TESTS::

            sage: polytopes.simplex(backend='cdd').prism().backend()
            'cdd'
        """
        from itertools import chain
        new_verts = chain(([0] + v for v in self.vertices()),
                          ([1] + v for v in self.vertices()))
        new_rays = ([0] + r for r in self.rays())
        new_lines = ([0] + l for l in self.lines())
        new_eqns = ([e.b()] + [0] + list(e[1:]) for e in self.equations())
        new_ieqs = chain(([i.b()] + [0] + list(i[1:]) for i in self.inequalities()),
                         [[0, 1] + [0]*self.ambient_dim(), [1, -1] + [0]*self.ambient_dim()])

        pref_rep = 'Hrep' if 2*(self.n_vertices() + self.n_rays()) >= self.n_inequalities() + 2 else 'Vrep'
        parent = self.parent().change_ring(self.base_ring(), ambient_dim=self.ambient_dim()+1)

        if not self.vertices():
            # Fix the empty polyhedron.
            return parent.element_class(parent, [[], [], []], None)

        return parent.element_class(parent, [new_verts, new_rays, new_lines],
                                    [new_ieqs, new_eqns],
                                    Vrep_minimal=True, Hrep_minimal=True, pref_rep=pref_rep)

    def _test_prism(self, tester=None, **options):
        """
        Run tests on the method :meth:`.prism`.

        TESTS::

            sage: polytopes.cross_polytope(3)._test_prism()
        """
        if tester is None:
            tester = self._tester(**options)

        if self.n_vertices() + self.n_rays() < 40 and self.n_facets() < 40:
            prism = self.prism()

            # Check that the double description is set up correctly.
            if self.base_ring().is_exact() and self.n_vertices() + self.n_rays() < 15 and self.n_facets() < 15:
                prism._test_basic_properties(tester)

            # Check that the prism preserves the backend.
            tester.assertEqual(prism.backend(), self.backend())

            tester.assertEqual(2*self.n_vertices(), prism.n_vertices())
            tester.assertEqual(self.n_rays(), prism.n_rays())
            tester.assertEqual(self.n_lines(), prism.n_lines())
            tester.assertEqual(self.n_equations(), prism.n_equations())
            if self.is_empty():
                return

            tester.assertEqual(2 + self.n_inequalities(), prism.n_inequalities())

            if not self.is_compact():
                # ``is_prism`` only implemented for compact polyhedra.
                return

            b, cert = prism.is_prism(certificate=True)
            tester.assertTrue(b)

            if not self.is_prism() and self.base_ring().is_exact():
                # In this case the certificate is unique.

                R = self.base_ring()
                cert_set = set(frozenset(tuple(v) for v in f) for f in cert)
                expected_cert = set(frozenset((i,) + tuple(v)
                                              for v in self.vertices())
                                    for i in (R(0), R(1)))
                tester.assertEqual(cert_set, expected_cert)

    def truncation(self, cut_frac=None):
        r"""
        Return a new polyhedron formed from two points on each edge
        between two vertices.

        INPUT:

        - ``cut_frac`` -- integer; how deeply to cut into the edge
          Default is `\frac{1}{3}`

        OUTPUT: a Polyhedron object, truncated as described above

        EXAMPLES::

            sage: cube = polytopes.hypercube(3)
            sage: trunc_cube = cube.truncation()
            sage: trunc_cube.n_vertices()
            24
            sage: trunc_cube.n_inequalities()
            14

        TESTS::

            sage: polytopes.simplex(backend='field').truncation().backend()
            'field'
        """
        if cut_frac is None:
            cut_frac = ZZ.one() / 3

        new_vertices = []
        for e in self.bounded_edges():
            new_vertices.append((1 - cut_frac) * e[0]() + cut_frac * e[1]())
            new_vertices.append(cut_frac * e[0]() + (1 - cut_frac) * e[1]())

        new_vertices = [list(v) for v in new_vertices]
        new_rays = self.rays()
        new_lines = self.lines()

        parent = self.parent().base_extend(cut_frac)
        return parent.element_class(parent, [new_vertices, new_rays, new_lines], None)

    def lawrence_polytope(self):
        r"""
        Return the Lawrence polytope of ``self``.

        Let `P` be a `d`-polytope in `\RR^r` with `n` vertices. The Lawrence
        polytope of `P` is the polytope whose vertices are the columns of the
        following `(r+n)`-by-`2n` matrix.

        .. MATH::

            \begin{pmatrix}
             V      &   V    \\
             I_n    &   2I_n
            \end{pmatrix},

        where `V` is the `r`-by-`n` vertices matrix of `P`.

        EXAMPLES::

            sage: P = polytopes.octahedron()
            sage: L = P.lawrence_polytope(); L
            A 9-dimensional polyhedron in ZZ^9 defined as the convex hull of 12 vertices
            sage: V = P.vertices_list()
            sage: for i, v in enumerate(V):
            ....:     v = v + i*[0]
            ....:     P = P.lawrence_extension(v)
            sage: P == L
            True

        REFERENCES:

            For more information, see Section 6.6 of [Zie2007]_.
        """
        from sage.matrix.constructor import block_matrix

        if not self.is_compact():
            raise NotImplementedError("self must be a polytope")

        V = self.vertices_matrix().transpose()
        n = self.n_vertices()
        I_n = matrix.identity(n)
        lambda_V = block_matrix([[V, I_n], [V, 2*I_n]])
        parent = self.parent().change_ring(self.base_ring(), ambient_dim=self.ambient_dim() + n)
        return parent.element_class(parent, [lambda_V, [], []], None)

    def deformation_cone(self):
        r"""
        Return the deformation cone of ``self``.

        Let `P` be a `d`-polytope in `\RR^r` with `n` facets. The deformation
        cone is a polyhedron in `\RR^n` who points are the right-hand side `b`
        in `Ax\leq b` where `A` is the matrix of facet normals of ``self``, so
        that the resulting polytope has a normal fan which is a coarsening of
        the normal fan of ``self``.

        EXAMPLES:

        Let's examine the deformation cone of the square with one truncated
        vertex::

            sage: tc = Polyhedron([(1, -1), (1/3, 1), (1, 1/3), (-1, 1), (-1, -1)])
            sage: dc = tc.deformation_cone()
            sage: dc.an_element()
            (2, 1, 1, 0, 0)
            sage: [_.A() for _ in tc.Hrepresentation()]
            [(1, 0), (0, 1), (0, -1), (-3, -3), (-1, 0)]
            sage: P = Polyhedron(rays=[(1, 0, 2), (0, 1, 1), (0, -1, 1), (-3, -3, 0), (-1, 0, 0)])
            sage: P.rays()
            (A ray in the direction (-1, -1, 0),
             A ray in the direction (-1, 0, 0),
             A ray in the direction (0, -1, 1),
             A ray in the direction (0, 1, 1),
             A ray in the direction (1, 0, 2))

        Now, let's compute the deformation cone of the pyramid over a square
        and verify that it is not full dimensional::

            sage: py = Polyhedron([(0, -1, -1), (0, -1, 1), (0, 1, -1), (0, 1, 1), (1, 0, 0)])
            sage: dc_py = py.deformation_cone(); dc_py
            A 4-dimensional polyhedron in QQ^5 defined as the convex hull of 1 vertex, 1 ray, 3 lines
            sage: [_.b() for _ in py.Hrepresentation()]
            [0, 1, 1, 1, 1]
            sage: r = dc_py.rays()[0]
            sage: l1,l2,l3 = dc_py.lines()
            sage: r.vector()-l1.vector()/2-l2.vector()-l3.vector()/2
            (0, 1, 1, 1, 1)

        .. SEEALSO::

            :meth:`~sage.schemes.toric.variety.Kaehler_cone`

        REFERENCES:

            For more information, see Section 5.4 of [DLRS2010]_ and Section
            2.2 of [ACEP2020].
        """
        from .constructor import Polyhedron
        A = matrix([_.A() for _ in self.Hrepresentation()])
        A = A.transpose()
        A_ker = A.right_kernel_matrix(basis='computed')
        gale = tuple(A_ker.columns())
        collection = [f.ambient_H_indices() for f in self.faces(0)]
        n = len(gale)
        K = None
        for cone_indices in collection:
            dual_cone = Polyhedron(rays=[gale[i] for i in range(n) if i not in
                                         cone_indices])
            K = K.intersection(dual_cone) if K is not None else dual_cone
        preimages = [A_ker.solve_right(r.vector()) for r in K.rays()]
        return Polyhedron(lines=A.rows(), rays=preimages)

    ###########################################################
    # Binary operations.
    ###########################################################

    @coerce_binop
    def minkowski_sum(self, other):
        r"""
        Return the Minkowski sum.

        Minkowski addition of two subsets of a vector space is defined
        as

        .. MATH::

            X \oplus Y =
            \bigcup_{y\in Y} (X+y) =
            \bigcup_{x\in X, y\in Y} (x+y)

        See :meth:`minkowski_difference` for a partial inverse operation.

        INPUT:

        - ``other`` -- a :class:`~sage.geometry.polyhedron.base.Polyhedron_base`

        OUTPUT: the Minkowski sum of ``self`` and ``other``

        EXAMPLES::

            sage: X = polytopes.hypercube(3)
            sage: Y = Polyhedron(vertices=[(0,0,0), (0,0,1/2), (0,1/2,0), (1/2,0,0)])
            sage: X+Y
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 13 vertices

            sage: four_cube = polytopes.hypercube(4)
            sage: four_simplex = Polyhedron(vertices=[[0, 0, 0, 1], [0, 0, 1, 0],
            ....:                                     [0, 1, 0, 0], [1, 0, 0, 0]])
            sage: four_cube + four_simplex
            A 4-dimensional polyhedron in ZZ^4 defined as the convex hull of 36 vertices
            sage: four_cube.minkowski_sum(four_simplex) == four_cube + four_simplex
            True

            sage: poly_spam = Polyhedron([[3,4,5,2], [1,0,0,1], [0,0,0,0],
            ....:                         [0,4,3,2], [-3,-3,-3,-3]], base_ring=ZZ)
            sage: poly_eggs = Polyhedron([[5,4,5,4], [-4,5,-4,5],
            ....:                         [4,-5,4,-5], [0,0,0,0]], base_ring=QQ)
            sage: poly_spam + poly_spam + poly_eggs
            A 4-dimensional polyhedron in QQ^4 defined as the convex hull of 12 vertices
        """
        new_vertices = []
        for v1 in self.vertex_generator():
            for v2 in other.vertex_generator():
                new_vertices.append(list(v1() + v2()))
        if new_vertices != []:
            new_rays = self.rays() + other.rays()
            new_lines = self.lines() + other.lines()
            return self.parent().element_class(self.parent(), [new_vertices, new_rays, new_lines], None)
        else:
            return self.parent().element_class(self.parent(), None, None)

    _add_ = minkowski_sum

    @coerce_binop
    def minkowski_difference(self, other):
        r"""
        Return the Minkowski difference.

        Minkowski subtraction can equivalently be defined via
        Minkowski addition (see :meth:`minkowski_sum`) or as
        set-theoretic intersection via

        .. MATH::

            X \ominus Y =
            (X^c \oplus Y)^c =
            \bigcap_{y\in Y} (X-y)

        where superscript-"c" means the complement in the ambient
        vector space. The Minkowski difference of convex sets is
        convex, and the difference of polyhedra is again a
        polyhedron. We only consider the case of polyhedra in the
        following. Note that it is not quite the inverse of
        addition. In fact:

        * `(X+Y)-Y = X` for any polyhedra `X`, `Y`.

        * `(X-Y)+Y \subseteq X`

        * `(X-Y)+Y = X` if and only if Y is a Minkowski summand of X.

        INPUT:

        - ``other`` -- a :class:`~sage.geometry.polyhedron.base.Polyhedron_base`

        OUTPUT:

        The Minkowski difference of ``self`` and ``other``. Also known
        as Minkowski subtraction of ``other`` from ``self``.

        EXAMPLES::

            sage: X = polytopes.hypercube(3)
            sage: Y = Polyhedron(vertices=[(0,0,0), (0,0,1), (0,1,0), (1,0,0)]) / 2
            sage: (X+Y)-Y == X
            True
            sage: (X-Y)+Y < X
            True

        The polyhedra need not be full-dimensional::

            sage: X2 = Polyhedron(vertices=[(-1,-1,0), (1,-1,0), (-1,1,0), (1,1,0)])
            sage: Y2 = Polyhedron(vertices=[(0,0,0), (0,1,0), (1,0,0)]) / 2
            sage: (X2+Y2)-Y2 == X2
            True
            sage: (X2-Y2)+Y2 < X2
            True

        Minus sign is really an alias for :meth:`minkowski_difference`
        ::

            sage: four_cube = polytopes.hypercube(4)
            sage: four_simplex = Polyhedron(vertices=[[0, 0, 0, 1], [0, 0, 1, 0],
            ....:                                     [0, 1, 0, 0], [1, 0, 0, 0]])
            sage: four_cube - four_simplex
            A 4-dimensional polyhedron in QQ^4 defined as the convex hull of 16 vertices
            sage: four_cube.minkowski_difference(four_simplex) == four_cube - four_simplex
            True

        Coercion of the base ring works::

            sage: poly_spam = Polyhedron([[3,4,5,2], [1,0,0,1], [0,0,0,0],
            ....:                         [0,4,3,2], [-3,-3,-3,-3]], base_ring=ZZ)
            sage: poly_eggs = Polyhedron([[5,4,5,4], [-4,5,-4,5],
            ....:                         [4,-5,4,-5], [0,0,0,0]], base_ring=QQ) / 100
            sage: poly_spam - poly_eggs
            A 4-dimensional polyhedron in QQ^4 defined as the convex hull of 5 vertices

        TESTS::

            sage: X = polytopes.hypercube(2)
            sage: Y = Polyhedron(vertices=[(1,1)])
            sage: (X-Y).Vrepresentation()
            (A vertex at (0, -2), A vertex at (0, 0), A vertex at (-2, 0), A vertex at (-2, -2))

            sage: Y = Polyhedron(vertices=[(1,1), (0,0)])
            sage: (X-Y).Vrepresentation()
            (A vertex at (0, -1), A vertex at (0, 0), A vertex at (-1, 0), A vertex at (-1, -1))

            sage: X = X + Y   # now Y is a Minkowski summand of X
            sage: (X+Y)-Y == X
            True
            sage: (X-Y)+Y == X
            True

        Testing that :issue:`28506` is fixed::

            sage: Q = Polyhedron([[1,0],[0,1]])
            sage: S = Polyhedron([[0,0],[1,2]])
            sage: S.minkowski_difference(Q)
            A 1-dimensional polyhedron in QQ^2 defined as the convex hull of 2 vertices
        """
        if other.is_empty():
            return self.parent().universe()   # empty intersection = everything
        if not other.is_compact():
            raise NotImplementedError('only subtracting compact polyhedra is implemented')
        new_eqns = []
        for eq in self.equations():
            values = [ eq.A() * v.vector() for v in other.vertices() ]
            eq = list(eq)
            eq[0] += min(values)   # shift constant term
            new_eqns.append(eq)
        P = self.parent()
        new_ieqs = []
        for ieq in self.inequalities():
            values = [ ieq.A() * v.vector() for v in other.vertices() ]
            ieq = list(ieq)
            ieq[0] += min(values)   # shift constant term
            new_ieqs.append(ieq)

        # Some vertices might need fractions.
        P = self.parent().change_ring(self.base_ring().fraction_field())
        return P.element_class(P, None, [new_ieqs, new_eqns])

    def __sub__(self, other):
        r"""
        Implement minus binary operation.

        Polyhedra are not a ring with respect to dilatation and
        Minkowski sum, for example `X\oplus(-1)*Y \not= X\ominus Y`.

        INPUT:

        - ``other`` -- a translation vector or a polyhedron

        OUTPUT:

        Either translation by the negative of the given vector or
        Minkowski subtraction by the given polyhedron.

        EXAMPLES::

            sage: X = polytopes.hypercube(2)
            sage: v = vector([1,1])
            sage: (X - v/2).Vrepresentation()
            (A vertex at (-3/2, -3/2), A vertex at (-3/2, 1/2),
             A vertex at (1/2, -3/2), A vertex at (1/2, 1/2))
            sage: (X-v)+v == X
            True

            sage: Y = Polyhedron(vertices=[(1/2,0), (0,1/2)])
            sage: (X-Y).Vrepresentation()
            (A vertex at (1/2, -1), A vertex at (1/2, 1/2),
             A vertex at (-1, 1/2), A vertex at (-1, -1))
            sage: (X+Y)-Y == X
            True
        """
        if isinstance(other, Polyhedron_base5):
            return self.minkowski_difference(other)
        return self + (-other)

    def product(self, other):
        """
        Return the Cartesian product.

        INPUT:

        - ``other`` -- a :class:`~sage.geometry.polyhedron.base.Polyhedron_base`

        OUTPUT:

        The Cartesian product of ``self`` and ``other`` with a
        suitable base ring to encompass the two.

        EXAMPLES::

            sage: P1 = Polyhedron([[0], [1]], base_ring=ZZ)
            sage: P2 = Polyhedron([[0], [1]], base_ring=QQ)
            sage: P1.product(P2)
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 4 vertices

        The Cartesian product is the product in the semiring of polyhedra::

            sage: P1 * P1
            A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 4 vertices
            sage: P1 * P2
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 4 vertices
            sage: P2 * P2
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 4 vertices
            sage: 2 * P1
            A 1-dimensional polyhedron in ZZ^1 defined as the convex hull of 2 vertices
            sage: P1 * 2.0
            A 1-dimensional polyhedron in RDF^1 defined as the convex hull of 2 vertices

        An alias is :meth:`cartesian_product`::

            sage: P1.cartesian_product(P2) == P1.product(P2)
            True

        TESTS:

        Check that :issue:`15253` is fixed::

            sage: polytopes.hypercube(1) * polytopes.hypercube(2)
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 8 vertices
        """
        try:
            new_ring = self.parent()._coerce_base_ring(other)
        except TypeError:
            raise TypeError("no common canonical parent for objects with parents: " + str(self.parent())
                             + " and " + str(other.parent()))

        from itertools import chain

        new_vertices = (tuple(x) + tuple(y)
                        for x in self.vertex_generator() for y in other.vertex_generator())

        self_zero = tuple(0 for _ in range( self.ambient_dim()))
        other_zero = tuple(0 for _ in range(other.ambient_dim()))

        rays = chain((tuple(r) + other_zero for r in self.ray_generator()),
                     (self_zero + tuple(r) for r in other.ray_generator()))

        lines = chain((tuple(l) + other_zero for l in self.line_generator()),
                      (self_zero + tuple(l) for l in other.line_generator()))

        if self.n_vertices() == 0 or other.n_vertices() == 0:
            # In this case we obtain the empty polyhedron.
            # There is not vertex to attach the rays or lines to.
            # By our convention, in this case the polyhedron shall also not have rays or lines.
            rays = ()
            lines = ()

        ieqs = chain((tuple(i) + other_zero
                      for i in self.inequality_generator()),
                     ((i.b(),) + self_zero + tuple(i.A())
                      for i in other.inequality_generator()))

        eqns = chain((tuple(e) + other_zero
                      for e in self.equation_generator()),
                     ((e.b(),) + self_zero + tuple(e.A())
                      for e in other.equation_generator()))

        pref_rep = 'Vrep' if self.n_vertices() + self.n_rays() + other.n_vertices() + other.n_rays() \
                             <= self.n_inequalities() + other.n_inequalities() else 'Hrep'

        parent = self.parent().change_ring(new_ring, ambient_dim=self.ambient_dim() + other.ambient_dim())
        return parent.element_class(parent, [new_vertices, rays, lines],
                                    [ieqs, eqns],
                                    Vrep_minimal=True, Hrep_minimal=True, pref_rep=pref_rep)

    _mul_ = product

    cartesian_product = product

    def _test_product(self, tester=None, **options):
        """
        Run tests on the method :meth:`.product`.

        TESTS::

            sage: polytopes.cross_polytope(3)._test_product()
        """
        from sage.rings.real_double import RDF
        from .library import polytopes

        if tester is None:
            tester = self._tester(**options)

        if self.n_vertices() + self.n_rays() < 40 and self.n_facets() < 40:
            # Check that the product preserves the backend, where possible.
            P = polytopes.simplex(backend='cdd')
            tester.assertEqual((self*P).backend(), self.backend())
            Q = polytopes.simplex(backend='ppl')
            tester.assertEqual((self*Q).backend(), self.backend())

            # And that it changes the backend correctly where necessary.
            try:
                from sage.rings.qqbar import AA
            except ImportError:
                pass
            else:
                if self.base_ring() is not AA and AA.has_coerce_map_from(self.base_ring()):
                    R = self*polytopes.regular_polygon(5, exact=True)
                    assert R
                if RDF.has_coerce_map_from(self.base_ring()):
                    R = self*polytopes.regular_polygon(5, exact=False)
                    assert R

        if self.base_ring() in (ZZ, QQ):
            # Check that the double description is set up correctly.
            self_field = self.base_extend(self.base_ring(), backend='field')
            try:
                P = polytopes.permutahedron(4, backend='field').base_extend(QQ)
            except ImportError:
                pass
            else:
                (self_field * P)._test_basic_properties(tester)
            from .constructor import Polyhedron
            Q = Polyhedron(rays=[[1,0,0,0],[0,1,1,0]], lines=[[0,1,0,1]], backend='field')
            (self_field * Q)._test_basic_properties(tester)

    def join(self, other):
        """
        Return the join of ``self`` and ``other``.

        The join of two polyhedra is obtained by first placing the two objects in
        two non-intersecting affine subspaces `V`, and `W` whose affine hull is
        the whole ambient space, and finally by taking the convex hull of their
        union. The dimension of the join is the sum of the dimensions of the
        two polyhedron plus 1.

        INPUT:

        - ``other`` -- a polyhedron

        EXAMPLES::

            sage: P1 = Polyhedron([[0],[1]], base_ring=ZZ)
            sage: P2 = Polyhedron([[0],[1]], base_ring=QQ)
            sage: P1.join(P2)
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 4 vertices
            sage: P1.join(P1)
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 4 vertices
            sage: P2.join(P2)
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 4 vertices

        An unbounded example::

            sage: R1 = Polyhedron(rays=[[1]])
            sage: R1.join(R1)
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 2 vertices and 2 rays

        TESTS::

            sage: C = polytopes.hypercube(5)
            sage: S = Polyhedron([[1]])
            sage: C.join(S).is_combinatorially_isomorphic(C.pyramid())                  # needs sage.graphs
            True

            sage: P = polytopes.simplex(backend='cdd')
            sage: Q = polytopes.simplex(backend='ppl')
            sage: P.join(Q).backend()
            'cdd'
            sage: Q.join(P).backend()
            'ppl'

        Check that the double description is set up correctly::

            sage: P = polytopes.cross_polytope(4)
            sage: P1 = polytopes.cross_polytope(4, backend='field')
            sage: P.join(P) == P1.join(P1)
            True

            sage: P = 4*polytopes.hypercube(4)
            sage: P1 = 4*polytopes.hypercube(4, backend='field')
            sage: P.join(P) == P1.join(P1)
            True

            sage: P = polytopes.permutahedron(4)
            sage: P1 = polytopes.permutahedron(4, backend='field')
            sage: P.join(P) == P1.join(P1)
            True
        """
        try:
            new_ring = self.parent()._coerce_base_ring(other)
        except TypeError:
            raise TypeError("no common canonical parent for objects with parents: " + str(self.parent())
                     + " and " + str(other.parent()))

        from itertools import chain

        dim_self = self.ambient_dim()
        dim_other = other.ambient_dim()
        parent = self.parent().change_ring(new_ring, ambient_dim=self.ambient_dim() + other.ambient_dim() + 1)

        new_vertices1 = (list(x) + [0]*dim_other + [0] for x in self.vertex_generator())
        new_vertices2 = ([0]*dim_self + list(x) + [1] for x in other.vertex_generator())
        new_vertices = chain(new_vertices1, new_vertices2)

        new_rays1 = (list(r) + [0]*dim_other + [0] for r in self.ray_generator())
        new_rays2 = ([0]*dim_self + list(r) + [1] for r in other.ray_generator())
        new_rays = chain(new_rays1, new_rays2)

        new_lines1 = (list(l) + [0]*dim_other + [0] for l in self.line_generator())
        new_lines2 = ([0]*dim_self + list(l) + [1] for l in other.line_generator())
        new_lines = chain(new_lines1, new_lines2)

        if not self.is_compact() or not other.is_compact() or self.n_vertices() <= 1 or other.n_vertices() <= 1:
            # Cases for which the below double description does not work.
            return parent.element_class(parent, [new_vertices, new_rays, new_lines], None)

        # Facet defining inequalities that contain the corresponding vertices from ``new_vertices1``
        # and all vertices from ``new_vertices2``.
        new_inequalities1 = ([i[0]] + list(i[1:]) + [0]*dim_other + [-i[0]] for i in self.inequality_generator())

        # Facet defining inequalities that contain the corresponding vertices from ``new_vertices2``
        # and all vertices from ``new_vertices1``.
        new_inequalities2 = ([0] + [0]*dim_self + list(i[1:]) + [i[0]] for i in other.inequality_generator())

        new_inequalities = chain(new_inequalities1, new_inequalities2)

        # Equations that all vertices corresponding to ``new_vertices1`` satisfy.
        # For any vertex from ``new_vertices2`` the condition is trivial.
        new_equations1 = ([e[0]] + list(e[1:]) + [0]*dim_other + [-e[0]] for e in self.equation_generator())

        # Equations that all vertices corresponding to ``new_vertices2`` satisfy.
        # For any vertex from ``new_vertices1`` the condition is trivial.
        new_equations2 = ([0] + [0]*dim_self + list(e[1:]) + [e[0]] for e in other.equation_generator())

        new_equations = chain(new_equations1, new_equations2)

        new_n_inequalities = self.n_inequalities() + other.n_inequalities()
        new_n_vertices = self.n_vertices() + other.n_vertices()
        new_n_rays = self.n_rays() + other.n_rays()

        pref_rep = 'Vrep' if new_n_vertices + new_n_rays <= new_n_inequalities else 'Hrep'

        return parent.element_class(parent,
                                    [new_vertices, new_rays, new_lines],
                                    [new_inequalities, new_equations],
                                    Vrep_minimal=True, Hrep_minimal=True,
                                    pref_rep=pref_rep)

    def subdirect_sum(self, other):
        """
        Return the subdirect sum of ``self`` and ``other``.

        The subdirect sum of two polyhedron is a projection of the join of the
        two polytopes. It is obtained by placing the two objects in orthogonal subspaces
        intersecting at the origin.

        INPUT:

        - ``other`` -- a :class:`~sage.geometry.polyhedron.base.Polyhedron_base`

        EXAMPLES::

            sage: P1 = Polyhedron([[1], [2]], base_ring=ZZ)
            sage: P2 = Polyhedron([[3], [4]], base_ring=QQ)
            sage: sds = P1.subdirect_sum(P2); sds
            A 2-dimensional polyhedron in QQ^2
             defined as the convex hull of 4 vertices
            sage: sds.vertices()
            (A vertex at (0, 3),
             A vertex at (0, 4),
             A vertex at (1, 0),
             A vertex at (2, 0))

        .. SEEALSO::

            :meth:`join`
            :meth:`direct_sum`

        TESTS::

            sage: P = polytopes.simplex(backend='cdd')
            sage: Q = polytopes.simplex(backend='ppl')
            sage: P.subdirect_sum(Q).backend()
            'cdd'
            sage: Q.subdirect_sum(P).backend()
            'ppl'
        """
        try:
            new_ring = self.parent()._coerce_base_ring(other)
        except TypeError:
            raise TypeError("no common canonical parent for objects with parents: " + str(self.parent())
                     + " and " + str(other.parent()))

        dim_self = self.ambient_dim()
        dim_other = other.ambient_dim()

        new_vertices = [list(x)+[0]*dim_other for x in self.vertex_generator()] + \
                       [[0]*dim_self+list(x) for x in other.vertex_generator()]
        new_rays = []
        new_rays.extend( [ r+[0]*dim_other
                           for r in self.ray_generator() ] )
        new_rays.extend( [ [0]*dim_self+r
                           for r in other.ray_generator() ] )
        new_lines = []
        new_lines.extend( [ l+[0]*dim_other
                            for l in self.line_generator() ] )
        new_lines.extend( [ [0]*dim_self+l
                            for l in other.line_generator() ] )

        parent = self.parent().change_ring(new_ring, ambient_dim=self.ambient_dim() + other.ambient_dim())
        return parent.element_class(parent, [new_vertices, new_rays, new_lines], None)

    def direct_sum(self, other):
        """
        Return the direct sum of ``self`` and ``other``.

        The direct sum of two polyhedron is the subdirect sum of the two, when
        they have the origin in their interior. To avoid checking if the origin
        is contained in both, we place the affine subspace containing ``other``
        at the center of ``self``.

        INPUT:

        - ``other`` -- a :class:`~sage.geometry.polyhedron.base.Polyhedron_base`

        EXAMPLES::

            sage: P1 = Polyhedron([[1], [2]], base_ring=ZZ)
            sage: P2 = Polyhedron([[3], [4]], base_ring=QQ)
            sage: ds = P1.direct_sum(P2);ds
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 4 vertices
            sage: ds.vertices()
            (A vertex at (1, 0),
             A vertex at (2, 0),
             A vertex at (3/2, -1/2),
             A vertex at (3/2, 1/2))

        .. SEEALSO::

            :meth:`join`
            :meth:`subdirect_sum`

        TESTS:

        Check that the backend is preserved::

            sage: P = polytopes.simplex(backend='cdd')
            sage: Q = polytopes.simplex(backend='ppl')
            sage: P.direct_sum(Q).backend()
            'cdd'
            sage: Q.direct_sum(P).backend()
            'ppl'

        Check that :issue:`28506` is fixed::

            sage: s2 = polytopes.simplex(2)
            sage: s3 = polytopes.simplex(3)
            sage: s2.direct_sum(s3)
            A 5-dimensional polyhedron in QQ^7 defined as the convex hull of 7 vertices
        """
        try:
            # Some vertices might need fractions.
            new_ring = self.parent()._coerce_base_ring(other).fraction_field()
        except TypeError:
            raise TypeError("no common canonical parent for objects with parents: " + str(self.parent())
                     + " and " + str(other.parent()))

        dim_self = self.ambient_dim()
        dim_other = other.ambient_dim()

        new_vertices = [list(x) + [0]*dim_other for x in self.vertex_generator()] + \
                       [list(self.center()) + list(x.vector() - other.center()) for x in other.vertex_generator()]
        new_rays = []
        new_rays.extend( [ r + [0]*dim_other
                           for r in self.ray_generator() ] )
        new_rays.extend( [ [0]*dim_self + r
                           for r in other.ray_generator() ] )
        new_lines = []
        new_lines.extend( [ l + [0]*dim_other
                            for l in self.line_generator() ] )
        new_lines.extend( [ [0]*dim_self + l
                            for l in other.line_generator() ] )

        parent = self.parent().change_ring(new_ring, ambient_dim=self.ambient_dim() + other.ambient_dim())
        return parent.element_class(parent, [new_vertices, new_rays, new_lines], None)

    @coerce_binop
    def convex_hull(self, other):
        """
        Return the convex hull of the set-theoretic union of the two
        polyhedra.

        INPUT:

        - ``other`` -- a :class:`Polyhedron`

        OUTPUT: the convex hull

        EXAMPLES::

            sage: a_simplex = polytopes.simplex(3, project=True)
            sage: verts = a_simplex.vertices()
            sage: verts = [[x[0]*3/5 + x[1]*4/5, -x[0]*4/5 + x[1]*3/5, x[2]] for x in verts]
            sage: another_simplex = Polyhedron(vertices=verts)
            sage: simplex_union = a_simplex.convex_hull(another_simplex)
            sage: simplex_union.n_vertices()
            7
        """
        hull_vertices = self.vertices() + other.vertices()
        hull_rays = self.rays() + other.rays()
        hull_lines = self.lines() + other.lines()
        return self.parent().element_class(self.parent(), [hull_vertices, hull_rays, hull_lines], None)

    @coerce_binop
    def intersection(self, other):
        r"""
        Return the intersection of one polyhedron with another.

        INPUT:

        - ``other`` -- a :class:`Polyhedron`

        OUTPUT: the intersection

        Note that the intersection of two `\ZZ`-polyhedra might not be
        a `\ZZ`-polyhedron. In this case, a `\QQ`-polyhedron is
        returned.

        EXAMPLES::

            sage: cube = polytopes.hypercube(3)
            sage: oct = polytopes.cross_polytope(3)
            sage: cube.intersection(oct*2)
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 12 vertices

        As a shorthand, one may use::

            sage: cube & oct*2
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 12 vertices

        The intersection of two `\ZZ`-polyhedra is not necessarily a `\ZZ`-polyhedron::

            sage: P = Polyhedron([(0,0),(1,1)], base_ring=ZZ)
            sage: P.intersection(P)
            A 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices
            sage: Q = Polyhedron([(0,1),(1,0)], base_ring=ZZ)
            sage: P.intersection(Q)
            A 0-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex
            sage: _.Vrepresentation()
            (A vertex at (1/2, 1/2),)

        TESTS:

        Check that :issue:`19012` is fixed::

            sage: # needs sage.rings.number_field
            sage: K.<a> = QuadraticField(5)
            sage: P = Polyhedron([[0, 0], [0, a], [1, 1]])
            sage: Q = Polyhedron(ieqs=[[-1, a, 1]])
            sage: P.intersection(Q)
            A 2-dimensional polyhedron in
             (Number Field in a with defining polynomial x^2 - 5 with a = 2.236067977499790?)^2
             defined as the convex hull of 4 vertices
        """
        new_ieqs = self.inequalities() + other.inequalities()
        new_eqns = self.equations() + other.equations()
        parent = self.parent()
        try:
            intersection = parent.element_class(parent, None, [new_ieqs, new_eqns])

            # Force calculation of the vertices.
            _ = intersection.n_vertices()
            return intersection
        except TypeError as msg:
            if self.base_ring() is ZZ:
                parent = parent.base_extend(QQ)
                return parent.element_class(parent, None, [new_ieqs, new_eqns])
            else:
                raise TypeError(msg)

    __and__ = intersection

    ###########################################################
    # Actions.
    ###########################################################

    def _acted_upon_(self, actor, self_on_left):
        """
        Implement the action by scalars, vectors, matrices or other polyhedra.

        INPUT:

        - ``actor`` -- one of the following:
          - a scalar, not necessarily in :meth:`base_ring`,
          - a :class:`Polyhedron`,
          - a :class:`sage.modules.free_module_element.vector`,
          - a :class:`sage.matrix.constructor.matrix`,
        - ``self_on_right`` -- must be ``False`` for actor a matrix;
          ignored otherwise

        OUTPUT:

        - Dilation for a scalar
        - Product for a polyhedron
        - Translation for a vector
        - Linear transformation for a matrix

        EXAMPLES:

        ``actor`` is a scalar::

             sage: p = Polyhedron(vertices = [[t,t^2,t^3] for t in srange(2,6)])
             sage: p._acted_upon_(2, True) == p.dilation(2)
             True
             sage: p*2 == p.dilation(2)
             True

        ``actor`` is a polyhedron::

             sage: p*p == p.product(p)
             True

        ``actor`` is a vector::

             sage: p + vector(ZZ,[1,2,3]) == p.translation([1,2,3])
             True

        ``actor`` is a matrix::

             sage: matrix(ZZ,[[1,2,3]]) * p
             A 1-dimensional polyhedron in ZZ^1 defined as the convex hull of 2 vertices

        A matrix must act from the left::

             sage: p * matrix(ZZ, [[1,2,3]]*3)
             Traceback (most recent call last):
             ...
             ValueError: matrices should act on the left
        """
        if isinstance(actor, Polyhedron_base5):
            return self.product(actor)
        elif isinstance(actor, Vector):
            return self.translation(actor)
        elif isinstance(actor, Matrix):
            if self_on_left:
                raise ValueError("matrices should act on the left")
            else:
                return self.linear_transformation(actor)
        else:
            return self.dilation(actor)

    def translation(self, displacement):
        """
        Return the translated polyhedron.

        INPUT:

        - ``displacement`` -- a displacement vector or a list/tuple of
          coordinates that determines a displacement vector

        OUTPUT: the translated polyhedron

        EXAMPLES::

            sage: P = Polyhedron([[0,0], [1,0], [0,1]], base_ring=ZZ)
            sage: P.translation([2,1])
            A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 3 vertices
            sage: P.translation(vector(QQ, [2,1]))
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 3 vertices

        TESTS::

            sage: P = Polyhedron([[0,0], [1,0], [0,1]], base_ring=ZZ, backend='field')
            sage: P.translation([2,1]).backend()
            'field'

        Check that precomputed data is set up correctly::

            sage: P = polytopes.permutahedron(4)*Polyhedron(lines=[[1]])
            sage: Q = P.change_ring(P.base_ring(), backend='field')
            sage: P + vector([1,2,3,4,5]) == Q + vector([1,2,3,4,5])
            True
            sage: P + vector([1,2,3,4,5/2]) == Q + vector([1,2,3,4,5/2])
            True
        """
        Vrep, Hrep, parent = self._translation_double_description(displacement)

        pref_rep = 'Vrep' if self.n_vertices() + self.n_rays() <= self.n_inequalities() else 'Hrep'

        return parent.element_class(parent, Vrep, Hrep,
                                    Vrep_minimal=True, Hrep_minimal=True, pref_rep=pref_rep)

    def _translation_double_description(self, displacement):
        r"""
        Return the input parameters for the translation.

        INPUT:

        - ``displacement`` -- a displacement vector or a list/tuple of
          coordinates that determines a displacement vector

        OUTPUT: tuple of consisting of new Vrepresentation, Hrepresentation and parent

        .. SEEALSO::

            :meth:`translation`

        EXAMPLES::

            sage: P = Polyhedron([[0,0], [1,0], [0,1]], base_ring=ZZ)
            sage: Vrep, Hrep, parent = P._translation_double_description([2,1])
            sage: [tuple(x) for x in Vrep], [tuple(x) for x in Hrep], parent
            ([((2, 1), (2, 2), (3, 1)), (), ()],
             [((-2, 1, 0), (-1, 0, 1), (4, -1, -1)), ()],
             Polyhedra in ZZ^2)
        """
        displacement = vector(displacement)
        new_vertices = (x.vector()+displacement for x in self.vertex_generator())
        new_rays = self.rays()
        new_lines = self.lines()
        parent = self.parent().base_extend(displacement)

        # Replace a hyperplane of the form A*x + b >= 0 by
        # A(x-displacement) + b >= 0 <=> Ax + b - A*displacement >= 0.
        # Likewise for equations.
        def get_new(x):
            y = x.vector().change_ring(parent.base_ring())
            y[0] -= x.A()*displacement
            return y

        new_ieqs = (get_new(x) for x in self.inequality_generator())
        new_eqns = (get_new(x) for x in self.equation_generator())
        return [new_vertices, new_rays, new_lines], [new_ieqs, new_eqns], parent

    def dilation(self, scalar):
        """
        Return the dilated (uniformly stretched) polyhedron.

        INPUT:

        - ``scalar`` -- a scalar, not necessarily in :meth:`base_ring`

        OUTPUT:

        The polyhedron dilated by that scalar, possibly coerced to a
        bigger base ring.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[[t,t^2,t^3] for t in srange(2,6)])
            sage: next(p.vertex_generator())
            A vertex at (2, 4, 8)
            sage: p2 = p.dilation(2)
            sage: next(p2.vertex_generator())
            A vertex at (4, 8, 16)
            sage: p.dilation(2) == p * 2
            True

        TESTS:

        Dilation of empty polyhedra works, see :issue:`14987`::

            sage: p = Polyhedron(ambient_dim=2); p
            The empty polyhedron in ZZ^2
            sage: p.dilation(3)
            The empty polyhedron in ZZ^2

            sage: p = Polyhedron(vertices=[(1,1)], rays=[(1,0)], lines=[(0,1)])
            sage: (-p).rays()
            (A ray in the direction (-1, 0),)
            sage: (-p).lines()
            (A line in the direction (0, 1),)

            sage: (0*p).rays()
            ()
            sage: (0*p).lines()
            ()
        """
        parent = self.parent().base_extend(scalar)

        if scalar == 0:
            new_vertices = tuple(self.ambient_space().zero() for v in self.vertex_generator())
            new_rays = []
            new_lines = []
            return parent.element_class(parent, [new_vertices, new_rays, new_lines], None)

        one = parent.base_ring().one()
        sign = one if scalar > 0 else -one

        make_new_Hrep = lambda h: tuple(scalar*sign*x if i == 0 else sign*x
                                        for i, x in enumerate(h._vector))

        new_vertices = (tuple(scalar*x for x in v._vector) for v in self.vertex_generator())
        new_rays = (tuple(sign*x for x in r._vector) for r in self.ray_generator())
        new_lines = self.line_generator()
        new_inequalities = map(make_new_Hrep, self.inequality_generator())
        new_equations = map(make_new_Hrep, self.equation_generator())

        pref_rep = 'Vrep' if self.n_vertices() + self.n_rays() <= self.n_inequalities() else 'Hrep'

        return parent.element_class(parent, [new_vertices, new_rays, new_lines],
                                    [new_inequalities, new_equations],
                                    Vrep_minimal=True, Hrep_minimal=True, pref_rep=pref_rep)

    def __truediv__(self, scalar):
        """
        Divide by a scalar factor.

        See :meth:`dilation` for details.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[t,t^2,t^3] for t in srange(2,4)])
            sage: (p/5).Vrepresentation()
            (A vertex at (2/5, 4/5, 8/5), A vertex at (3/5, 9/5, 27/5))
            sage: (p/int(5)).Vrepresentation()
            (A vertex at (0.4, 0.8, 1.6), A vertex at (0.6, 1.8, 5.4))
        """
        return self.dilation(1/scalar)

    def _test_dilation(self, tester=None, **options):
        """
        Run tests on the method :meth:`.dilation`.

        TESTS::

            sage: polytopes.cross_polytope(3)._test_dilation()
        """
        from sage.rings.real_double import RDF
        from .base import Polyhedron_base

        if tester is None:
            tester = self._tester(**options)

        # Testing that the backend is preserved.
        tester.assertEqual(self.dilation(2*self.base_ring().gen()).backend(), self.backend())
        tester.assertEqual(self.dilation(ZZ(3)).backend(), self.backend())

        if self.n_vertices() + self.n_rays() > 40:
            # Avoid long time computations.
            return

        # Testing that the double description is set up correctly.
        if self.base_ring().is_exact():
            if self.base_ring() in (QQ, ZZ):
                p = self.base_extend(self.base_ring(), backend='field')
                (ZZ(2) * p)._test_basic_properties(tester)
                (ZZ(2)/2 * p)._test_basic_properties(tester)
                (ZZ(-3) * p)._test_basic_properties(tester)
                (ZZ(-1)/2 * p)._test_basic_properties(tester)
        else:
            tester.assertIsInstance(ZZ(1)/3*self, Polyhedron_base)

        try:
            from sage.rings.qqbar import AA
        except ImportError:
            return

        if self.n_vertices() > 20 or self.base_ring() is AA:
            # Avoid long time computations.
            return

        # Some sanity check on the volume (only run for relatively small instances).
        if self.dim() > -1 and self.is_compact() and self.base_ring().is_exact():
            tester.assertEqual(self.dilation(3).volume(measure='induced'), self.volume(measure='induced')*3**self.dim())

        # Testing coercion with algebraic numbers.
        from sage.rings.number_field.number_field import QuadraticField
        K1 = QuadraticField(2, embedding=AA(2).sqrt())
        sqrt2 = K1.gen()
        K2 = QuadraticField(3, embedding=AA(3).sqrt())
        sqrt3 = K2.gen()

        if self.base_ring() in (QQ, ZZ, AA, RDF):
            tester.assertIsInstance(sqrt2*self, Polyhedron_base)
            tester.assertIsInstance(sqrt3*self, Polyhedron_base)
        elif hasattr(self.base_ring(), "composite_fields"):
            for scalar, K in ((sqrt2, K1), (sqrt3, K2)):
                new_ring = None
                try:
                    new_ring = self.base_ring().composite_fields()[0]
                except (KeyError, AttributeError, TypeError):
                    # This isn't about testing composite fields.
                    pass
                if new_ring:
                    p = self.change_ring(new_ring)
                    tester.assertIsInstance(scalar*p, Polyhedron_base)

    def linear_transformation(self, linear_transf, new_base_ring=None):
        """
        Return the linear transformation of ``self``.

        INPUT:

        - ``linear_transf`` -- a matrix, not necessarily in :meth:`base_ring`
        - ``new_base_ring`` -- ring (optional); specify the new base ring;
          may avoid coercion failure

        OUTPUT:

        The polyhedron transformed by that matrix, possibly coerced to a
        bigger base ring.

        EXAMPLES::

            sage: b3 = polytopes.Birkhoff_polytope(3)
            sage: proj_mat = matrix([[0,1,0,0,0,0,0,0,0], [0,0,0,1,0,0,0,0,0],
            ....:                    [0,0,0,0,0,1,0,0,0], [0,0,0,0,0,0,0,1,0]])
            sage: b3_proj = proj_mat * b3; b3_proj
            A 3-dimensional polyhedron in ZZ^4 defined as the convex hull of 5 vertices

            sage: # needs sage.rings.number_field
            sage: square = polytopes.regular_polygon(4)
            sage: square.vertices_list()
            [[0, -1], [1, 0], [-1, 0], [0, 1]]
            sage: transf = matrix([[1,1], [0,1]])
            sage: sheared = transf * square
            sage: sheared.vertices_list()
            [[-1, -1], [1, 0], [-1, 0], [1, 1]]
            sage: sheared == square.linear_transformation(transf)
            True

        Specifying the new base ring may avoid coercion failure::

            sage: # needs sage.rings.number_field
            sage: K.<sqrt2> = QuadraticField(2)
            sage: L.<sqrt3> = QuadraticField(3)
            sage: P = polytopes.cube()*sqrt2
            sage: M = matrix([[sqrt3, 0, 0], [0, sqrt3, 0], [0, 0, 1]])
            sage: P.linear_transformation(M, new_base_ring=K.composite_fields(L)[0])
            A 3-dimensional polyhedron in
             (Number Field in sqrt2sqrt3 with defining polynomial x^4 - 10*x^2 + 1
              with sqrt2sqrt3 = 0.3178372451957823?)^3
             defined as the convex hull of 8 vertices

        Linear transformation without specified new base ring fails in this case::

            sage: M*P                                                                   # needs sage.rings.number_field
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for *:
            'Full MatrixSpace of 3 by 3 dense matrices over Number Field in sqrt3
            with defining polynomial x^2 - 3 with sqrt3 = 1.732050807568878?' and
            'Full MatrixSpace of 3 by 8 dense matrices over Number Field in sqrt2
            with defining polynomial x^2 - 2 with sqrt2 = 1.414213562373095?'

        TESTS:

        Linear transformation respects backend::

            sage: P = polytopes.simplex(backend='field')
            sage: t = matrix([[1,1,1,1],[0,1,1,1],[0,0,1,1],[0,0,0,1]])
            sage: P.linear_transformation(t).backend()
            'field'

        Check that coercion works::

            sage: (1.0 * proj_mat) * b3
            A 3-dimensional polyhedron in RDF^4 defined as the convex hull of 5 vertices
            sage: (1/1 * proj_mat) * b3
            A 3-dimensional polyhedron in QQ^4 defined as the convex hull of 5 vertices
            sage: (AA(2).sqrt() * proj_mat) * b3                                        # needs sage.rings.number_field
            A 3-dimensional polyhedron in AA^4 defined as the convex hull of 5 vertices

        Check that zero-matrices act correctly::

            sage: Matrix([]) * b3
            A 0-dimensional polyhedron in ZZ^0 defined as the convex hull of 1 vertex
            sage: Matrix([[0 for _ in range(9)]]) * b3
            A 0-dimensional polyhedron in ZZ^1 defined as the convex hull of 1 vertex
            sage: Matrix([[0 for _ in range(9)] for _ in range(4)]) * b3
            A 0-dimensional polyhedron in ZZ^4 defined as the convex hull of 1 vertex
            sage: Matrix([[0 for _ in range(8)]]) * b3
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for *:
            'Full MatrixSpace of 1 by 8 dense matrices over Integer Ring' and
            'Full MatrixSpace of 9 by 6 dense matrices over Integer Ring'
            sage: Matrix(ZZ, []) * b3
            A 0-dimensional polyhedron in ZZ^0 defined as the convex hull of 1 vertex
            sage: Matrix(ZZ, [[],[]]) * b3
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for *:
            'Full MatrixSpace of 2 by 0 dense matrices over Integer Ring' and
            'Full MatrixSpace of 9 by 6 dense matrices over Integer Ring'

        Check that the precomputed double description is correct::

            sage: P = polytopes.permutahedron(4)
            sage: Q = P.change_ring(QQ, backend='field')
            sage: P.affine_hull_projection() == Q.affine_hull_projection()
            True

            sage: M = matrix([[1, 2, 3, 4], [2, 3, 4, 5], [0, 0, 5, 1], [0, 2, 0, 3]])
            sage: M*P == M*Q
            True

            sage: M = matrix([[1, 2, 3, 4], [2, 3, 4, 5], [0, 0, 5, 1], [0, 2, 0, 3], [0, 1, 0, -3]])
            sage: M*P == M*Q
            True
        """
        is_injective = False
        if linear_transf.nrows() != 0:
            if new_base_ring:
                R = new_base_ring
            else:
                R = self.base_ring()

            # Multiplying a matrix with a vector is slow.
            # So we multiply the entire vertex matrix etc.
            # Still we create generators, as possibly the Vrepresentation will be discarded later on.
            if self.n_vertices():
                new_vertices = ( v for v in ((linear_transf*self.vertices_matrix(R)).transpose()) )
            else:
                new_vertices = ()
            if self.n_rays():
                new_rays = ( r for r in matrix(R, self.rays())*linear_transf.transpose() )
            else:
                new_rays = ()
            if self.n_lines():
                new_lines = ( l for l in matrix(R, self.lines())*linear_transf.transpose() )
            else:
                new_lines = ()

            if self.is_compact() and self.n_vertices() and self.n_inequalities():
                homogeneous_basis = matrix(R, ( [1] + list(v) for v in self.an_affine_basis() )).transpose()

                # To convert first to a list and then to a matrix seems to be necessary to obtain a meaningful error,
                # in case the number of columns doesn't match the dimension.
                new_homogeneous_basis = matrix(list( [1] + list(linear_transf*vector(R, v)) for v in self.an_affine_basis()) ).transpose()

                if self.dim() + 1 == new_homogeneous_basis.rank():
                    # The transformation is injective on the polytope.
                    is_injective = True

                    # Let V be the homogeneous vertex matrix (each vertex a column)
                    # and M the linear transformation.
                    # Then M*V is the new homogeneous vertex matrix.

                    # Let H be the inequalities matrix (each inequality a row).
                    # If we find N such that N*M*V = V than the new inequalities are
                    # given by H*N.

                    # Note that such N must exist, as our map is injective on the polytope.
                    # It is uniquely defined by considering a basis of the homogeneous vertices.
                    N = new_homogeneous_basis.solve_left(homogeneous_basis)
                    new_inequalities = ( h for h in matrix(R, self.inequalities())*N )

                    # The equations are the left kernel matrix of the homogeneous vertices
                    # or equivalently a basis thereof.
                    new_equations = (new_homogeneous_basis.transpose()).right_kernel_matrix()

        else:
            new_vertices = [[] for v in self.vertex_generator() ]
            new_rays = []
            new_lines = []

        new_dim = linear_transf.nrows()
        par = self.parent()

        if new_base_ring:
            new_parent = par.change_ring(new_base_ring, ambient_dim=new_dim)
        else:
            new_parent = par.base_extend(linear_transf.base_ring(), ambient_dim=new_dim)

        if is_injective:
            # Set up with both Vrepresentation and Hrepresentation.
            pref_rep = 'Vrep' if self.n_vertices() <= self.n_inequalities() else 'Hrep'

            return new_parent.element_class(new_parent, [new_vertices, new_rays, new_lines],
                                            [new_inequalities, new_equations],
                                            Vrep_minimal=True, Hrep_minimal=True, pref_rep=pref_rep)

        return new_parent.element_class(new_parent, [tuple(new_vertices), tuple(new_rays), tuple(new_lines)], None)

    def _test_linear_transformation(self, tester=None, **options):
        """
        Run some tests on linear transformation.

        TESTS::

            sage: Polyhedron(rays=[(0,1)])._test_linear_transformation()
        """
        if tester is None:
            tester = self._tester(**options)

        if self.n_vertices() > 200 or self.n_facets() > 200:
            # Avoid very long doctests.
            return

        # Check that :issue:`30146` is fixed.
        from sage.matrix.special import identity_matrix
        tester.assertEqual(self, self.linear_transformation(identity_matrix(self.ambient_dim())))

    ###########################################################
    # Methods using a face.
    ###########################################################

    def face_truncation(self, face, linear_coefficients=None, cut_frac=None):
        r"""
        Return a new polyhedron formed by truncating a face by an hyperplane.

        By default, the normal vector of the hyperplane used to truncate the
        polyhedron is obtained by taking the barycenter vector of the cone
        corresponding to the truncated face in the normal fan of the
        polyhedron. It is possible to change the direction using the option
        ``linear_coefficients``.

        To determine how deep the truncation is done, the method uses the
        parameter ``cut_frac``. By default it is equal to `\frac{1}{3}`. Once
        the normal vector of the cutting hyperplane is chosen, the vertices of
        polyhedron are evaluated according to the corresponding linear
        function. The parameter `\frac{1}{3}` means that the cutting
        hyperplane is placed `\frac{1}{3}` of the way from the vertices of the
        truncated face to the next evaluated vertex.

        INPUT:

        - ``face`` -- a :class:`~sage.geometry.polyhedron.face.PolyhedronFace`
        - ``linear_coefficients`` -- tuple of integer. Specifies the coefficient
          of the normal vector of the cutting hyperplane used to truncate the
          face.
          The default direction is determined using the normal fan of the
          polyhedron.
        - ``cut_frac`` -- number between 0 and 1. Determines where the
          hyperplane cuts the polyhedron. A value close to 0 cuts very close
          to the face, whereas a value close to 1 cuts very close to the next
          vertex (according to the normal vector of the cutting hyperplane).
          Default is `\frac{1}{3}`.

        OUTPUT: a Polyhedron object, truncated as described above

        EXAMPLES::

            sage: Cube = polytopes.hypercube(3)
            sage: vertex_trunc1 = Cube.face_truncation(Cube.faces(0)[0])
            sage: vertex_trunc1.f_vector()
            (1, 10, 15, 7, 1)
            sage: tuple(f.ambient_V_indices() for f in vertex_trunc1.faces(2))
            ((4, 5, 6, 7, 9),
             (0, 3, 4, 8, 9),
             (0, 1, 6, 7, 8),
             (7, 8, 9),
             (2, 3, 4, 5),
             (1, 2, 5, 6),
             (0, 1, 2, 3))
            sage: vertex_trunc1.vertices()
            (A vertex at (1, -1, -1),
             A vertex at (1, 1, -1),
             A vertex at (1, 1, 1),
             A vertex at (1, -1, 1),
             A vertex at (-1, -1, 1),
             A vertex at (-1, 1, 1),
             A vertex at (-1, 1, -1),
             A vertex at (-1, -1/3, -1),
             A vertex at (-1/3, -1, -1),
             A vertex at (-1, -1, -1/3))
            sage: vertex_trunc2 = Cube.face_truncation(Cube.faces(0)[0], cut_frac=1/2)
            sage: vertex_trunc2.f_vector()
            (1, 10, 15, 7, 1)
            sage: tuple(f.ambient_V_indices() for f in vertex_trunc2.faces(2))
            ((4, 5, 6, 7, 9),
             (0, 3, 4, 8, 9),
             (0, 1, 6, 7, 8),
             (7, 8, 9),
             (2, 3, 4, 5),
             (1, 2, 5, 6),
             (0, 1, 2, 3))
            sage: vertex_trunc2.vertices()
            (A vertex at (1, -1, -1),
             A vertex at (1, 1, -1),
             A vertex at (1, 1, 1),
             A vertex at (1, -1, 1),
             A vertex at (-1, -1, 1),
             A vertex at (-1, 1, 1),
             A vertex at (-1, 1, -1),
             A vertex at (-1, 0, -1),
             A vertex at (0, -1, -1),
             A vertex at (-1, -1, 0))
            sage: vertex_trunc3 = Cube.face_truncation(Cube.faces(0)[0], cut_frac=0.3)
            sage: vertex_trunc3.vertices()
            (A vertex at (-1.0, -1.0, 1.0),
             A vertex at (-1.0, 1.0, -1.0),
             A vertex at (-1.0, 1.0, 1.0),
             A vertex at (1.0, 1.0, -1.0),
             A vertex at (1.0, 1.0, 1.0),
             A vertex at (1.0, -1.0, 1.0),
             A vertex at (1.0, -1.0, -1.0),
             A vertex at (-0.4, -1.0, -1.0),
             A vertex at (-1.0, -0.4, -1.0),
             A vertex at (-1.0, -1.0, -0.4))
            sage: edge_trunc = Cube.face_truncation(Cube.faces(1)[11])
            sage: edge_trunc.f_vector()
            (1, 10, 15, 7, 1)
            sage: tuple(f.ambient_V_indices() for f in edge_trunc.faces(2))
            ((0, 5, 6, 7),
             (1, 4, 5, 6, 8),
             (6, 7, 8, 9),
             (0, 2, 3, 7, 9),
             (1, 2, 8, 9),
             (0, 3, 4, 5),
             (1, 2, 3, 4))
             sage: face_trunc = Cube.face_truncation(Cube.faces(2)[2])
             sage: face_trunc.vertices()
             (A vertex at (1, -1, -1),
              A vertex at (1, 1, -1),
              A vertex at (1, 1, 1),
              A vertex at (1, -1, 1),
              A vertex at (-1/3, -1, 1),
              A vertex at (-1/3, 1, 1),
              A vertex at (-1/3, 1, -1),
              A vertex at (-1/3, -1, -1))
             sage: face_trunc.face_lattice().is_isomorphic(Cube.face_lattice())         # needs sage.combinat sage.graphs
             True

        TESTS:

        Testing that the backend is preserved::

            sage: Cube = polytopes.cube(backend='field')
            sage: face_trunc = Cube.face_truncation(Cube.faces(2)[0])
            sage: face_trunc.backend()
            'field'

        Testing that :issue:`28506` is fixed::

            sage: P = polytopes.twenty_four_cell()
            sage: P = P.dilation(6)
            sage: P = P.change_ring(ZZ)
            sage: P.face_truncation(P.faces(2)[0], cut_frac=1)
            A 4-dimensional polyhedron in QQ^4 defined as the convex hull of 27 vertices
        """
        if cut_frac is None:
            cut_frac = ZZ.one() / 3

        face_vertices = face.vertices()

        normal_vectors = []

        for facet in self.Hrepresentation():
            if all(facet.contains(x) and not facet.interior_contains(x)
                   for x in face_vertices):
                # The facet contains the face
                normal_vectors.append(facet.A())

        if linear_coefficients is not None:
            normal_vector = sum(linear_coefficients[i]*normal_vectors[i]
                                for i in range(len(normal_vectors)))
        else:
            normal_vector = sum(normal_vectors)

        B = - normal_vector * (face_vertices[0].vector())

        linear_evaluation = set(-normal_vector * (v.vector()) for v in self.vertices())

        if B == max(linear_evaluation):
            C = max(linear_evaluation.difference(set([B])))
        else:
            C = min(linear_evaluation.difference(set([B])))

        cut_height = (1 - cut_frac) * B + cut_frac * C
        ineq_vector = tuple([cut_height]) + tuple(normal_vector)

        new_ieqs = self.inequalities_list() + [ineq_vector]
        new_eqns = self.equations_list()

        # Some vertices might need fractions.
        parent = self.parent().base_extend(cut_frac/1)
        return parent.element_class(parent, None, [new_ieqs, new_eqns])

    def stack(self, face, position=None):
        r"""
        Return a new polyhedron formed by stacking onto a ``face``. Stacking a
        face adds a new vertex located slightly outside of the designated face.

        INPUT:

        - ``face`` -- a PolyhedronFace

        - ``position`` -- a positive number. Determines a relative distance
          from the barycenter of ``face``. A value close to 0 will place the
          new vertex close to the face and a large value further away. Default
          is `1`. If the given value is too large, an error is returned.

        OUTPUT: a Polyhedron object

        EXAMPLES::

            sage: cube = polytopes.cube()
            sage: square_face = cube.facets()[2]
            sage: stacked_square = cube.stack(square_face)
            sage: stacked_square.f_vector()
            (1, 9, 16, 9, 1)

            sage: edge_face = cube.faces(1)[3]
            sage: stacked_edge = cube.stack(edge_face)
            sage: stacked_edge.f_vector()
            (1, 9, 17, 10, 1)

            sage: cube.stack(cube.faces(0)[0])
            Traceback (most recent call last):
            ...
            ValueError: cannot stack onto a vertex

            sage: stacked_square_half = cube.stack(square_face, position=1/2)
            sage: stacked_square_half.f_vector()
            (1, 9, 16, 9, 1)
            sage: stacked_square_large = cube.stack(square_face, position=10)

            sage: # needs sage.rings.number_field
            sage: hexaprism = polytopes.regular_polygon(6).prism()
            sage: hexaprism.f_vector()
            (1, 12, 18, 8, 1)
            sage: square_face = hexaprism.faces(2)[2]
            sage: stacked_hexaprism = hexaprism.stack(square_face)
            sage: stacked_hexaprism.f_vector()
            (1, 13, 22, 11, 1)

            sage: hexaprism.stack(square_face, position=4)                              # needs sage.rings.number_field
            Traceback (most recent call last):
            ...
            ValueError: the chosen position is too large

            sage: s = polytopes.simplex(7)
            sage: f = s.faces(3)[69]
            sage: sf = s.stack(f); sf
            A 7-dimensional polyhedron in QQ^8 defined as the convex hull of 9 vertices
            sage: sf.vertices()
            (A vertex at (-4, -4, -4, -4, 17/4, 17/4, 17/4, 17/4),
             A vertex at (0, 0, 0, 0, 0, 0, 0, 1),
             A vertex at (0, 0, 0, 0, 0, 0, 1, 0),
             A vertex at (0, 0, 0, 0, 0, 1, 0, 0),
             A vertex at (0, 0, 0, 0, 1, 0, 0, 0),
             A vertex at (0, 0, 0, 1, 0, 0, 0, 0),
             A vertex at (0, 0, 1, 0, 0, 0, 0, 0),
             A vertex at (0, 1, 0, 0, 0, 0, 0, 0),
             A vertex at (1, 0, 0, 0, 0, 0, 0, 0))

        It is possible to stack on unbounded faces::

            sage: Q = Polyhedron(vertices=[[0,1], [1,0]], rays=[[1,1]])
            sage: E = Q.faces(1)
            sage: Q.stack(E[0],1/2).Vrepresentation()
            (A vertex at (0, 1),
             A vertex at (1, 0),
             A ray in the direction (1, 1),
             A vertex at (2, 0))
            sage: Q.stack(E[1],1/2).Vrepresentation()
            (A vertex at (0, 1),
             A vertex at (0, 2),
             A vertex at (1, 0),
             A ray in the direction (1, 1))
            sage: Q.stack(E[2],1/2).Vrepresentation()
            (A vertex at (0, 0),
             A vertex at (0, 1),
             A vertex at (1, 0),
             A ray in the direction (1, 1))

        Stacking requires a proper face::

            sage: Q.stack(Q.faces(2)[0])
            Traceback (most recent call last):
            ...
            ValueError: can only stack on proper face

        TESTS:

        Checking that the backend is preserved::

            sage: Cube = polytopes.cube(backend='field')
            sage: stack = Cube.stack(Cube.faces(2)[0])
            sage: stack.backend()
            'field'

        Taking the stacking vertex too far with the parameter ``position``
        may result in a failure to produce the desired
        (combinatorial type of) polytope.
        The interval of permitted values is always open.
        This is the smallest unpermitted value::

            sage: P = polytopes.octahedron()
            sage: P.stack(P.faces(2)[0], position=4)
            Traceback (most recent call last):
            ...
            ValueError: the chosen position is too large

        Testing that :issue:`29057` is fixed::

            sage: P = polytopes.cross_polytope(4)
            sage: P.stack(P.faces(3)[0])
            A 4-dimensional polyhedron in QQ^4 defined as the convex hull of 9 vertices
        """
        from sage.geometry.polyhedron.face import PolyhedronFace
        if not isinstance(face, PolyhedronFace):
            raise TypeError("{} should be a PolyhedronFace of {}".format(face, self))
        elif face.dim() == 0:
            raise ValueError("cannot stack onto a vertex")
        elif face.dim() == -1 or face.dim() == self.dim():
            raise ValueError("can only stack on proper face")
        if position is None:
            position = 1

        barycenter = ZZ.one()*sum([v.vector() for v in face.vertices()]) / len(face.vertices())
        locus_polyhedron = face.stacking_locus()
        repr_point = locus_polyhedron.representative_point()
        new_vertex = (1-position)*barycenter + position*repr_point
        if not locus_polyhedron.relative_interior_contains(new_vertex):
            raise ValueError("the chosen position is too large")

        parent = self.parent().base_extend(new_vertex)
        return parent.element_class(parent, [self.vertices() + (new_vertex,), self.rays(), self.lines()], None)

    def wedge(self, face, width=1):
        r"""
        Return the wedge over a ``face`` of the polytope ``self``.

        The wedge over a face `F` of a polytope `P` with width `w \not= 0`
        is defined as:

        .. MATH::

            (P \times \mathbb{R}) \cap \{a^\top x + |w x_{d+1}| \leq b\}

        where `\{x | a^\top x = b\}` is a supporting hyperplane defining `F`.

        INPUT:

        - ``face`` -- a PolyhedronFace of ``self``, the face which we take
          the wedge over
        - ``width`` -- a nonzero number (default: ``1``);
          specifies how wide the wedge will be

        OUTPUT:

        A (bounded) polyhedron

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: P_4 = polytopes.regular_polygon(4)
            sage: W1 = P_4.wedge(P_4.faces(1)[0]); W1
            A 3-dimensional polyhedron in AA^3 defined as the convex hull of 6 vertices
            sage: triangular_prism = polytopes.regular_polygon(3).prism()
            sage: W1.is_combinatorially_isomorphic(triangular_prism)                    # needs sage.graphs
            True

            sage: Q = polytopes.hypersimplex(4,2)
            sage: W2 = Q.wedge(Q.faces(2)[7]); W2
            A 4-dimensional polyhedron in QQ^5 defined as the convex hull of 9 vertices
            sage: W2.vertices()
            (A vertex at (1, 1, 0, 0, 1),
             A vertex at (1, 1, 0, 0, -1),
             A vertex at (1, 0, 1, 0, 1),
             A vertex at (1, 0, 1, 0, -1),
             A vertex at (1, 0, 0, 1, 1),
             A vertex at (1, 0, 0, 1, -1),
             A vertex at (0, 0, 1, 1, 0),
             A vertex at (0, 1, 1, 0, 0),
             A vertex at (0, 1, 0, 1, 0))

            sage: W3 = Q.wedge(Q.faces(1)[11]); W3
            A 4-dimensional polyhedron in QQ^5 defined as the convex hull of 10 vertices
            sage: W3.vertices()
            (A vertex at (1, 1, 0, 0, -2),
             A vertex at (1, 1, 0, 0, 2),
             A vertex at (1, 0, 1, 0, -2),
             A vertex at (1, 0, 1, 0, 2),
             A vertex at (1, 0, 0, 1, 1),
             A vertex at (1, 0, 0, 1, -1),
             A vertex at (0, 1, 0, 1, 0),
             A vertex at (0, 1, 1, 0, 1),
             A vertex at (0, 0, 1, 1, 0),
             A vertex at (0, 1, 1, 0, -1))

            sage: C_3_7 = polytopes.cyclic_polytope(3,7)
            sage: P_6 = polytopes.regular_polygon(6)                                    # needs sage.rings.number_field
            sage: W4 = P_6.wedge(P_6.faces(1)[0])                                       # needs sage.rings.number_field
            sage: W4.is_combinatorially_isomorphic(C_3_7.polar())                       # needs sage.graphs sage.rings.number_field
            True

        REFERENCES:

        For more information, see Chapter 15 of [HoDaCG17]_.

        TESTS:

        The backend should be preserved as long as the value of width permits.
        The base_ring will change to the field of fractions of the current
        base_ring, unless ``width`` forces a different ring. ::

            sage: P = polytopes.cyclic_polytope(3,7, base_ring=ZZ, backend='field')
            sage: W1 = P.wedge(P.faces(2)[0]); W1.base_ring(); W1.backend()
            Rational Field
            'field'
            sage: W2 = P.wedge(P.faces(2)[0], width=5/2); W2.base_ring(); W2.backend()
            Rational Field
            'field'
            sage: W2 = P.wedge(P.faces(2)[9], width=4/2); W2.base_ring(); W2.backend()
            Rational Field
            'field'
            sage: W2.vertices()
            (A vertex at (3, 9, 27, -1/2),
             A vertex at (4, 16, 64, -2),
             A vertex at (6, 36, 216, -10),
             A vertex at (5, 25, 125, -5),
             A vertex at (2, 4, 8, 0),
             A vertex at (1, 1, 1, 0),
             A vertex at (0, 0, 0, 0),
             A vertex at (3, 9, 27, 1/2),
             A vertex at (4, 16, 64, 2),
             A vertex at (6, 36, 216, 10),
             A vertex at (5, 25, 125, 5))
            sage: W2 = P.wedge(P.faces(2)[2], width=1.0); W2.base_ring(); W2.backend()
            Real Double Field
            'cdd'
        """
        width = width*ZZ.one()

        if not self.is_compact():
            raise ValueError("polyhedron 'self' must be a polytope")

        if width == 0:
            raise ValueError("the width should be nonzero")

        from sage.geometry.polyhedron.face import PolyhedronFace
        if not isinstance(face, PolyhedronFace):
            raise TypeError("{} should be a PolyhedronFace of {}".format(face, self))

        F_Hrep = vector([0]*(self.ambient_dim()+1))
        for facet in face.ambient_Hrepresentation():
            if facet.is_inequality():
                F_Hrep = F_Hrep + facet.vector()
        F_Hrep = list(F_Hrep)

        parent = self.parent()
        parent1 = parent.base_extend(self.base_ring(), ambient_dim=1)
        parent2 = parent.base_extend(width.base_ring().fraction_field(), ambient_dim=1 + self.ambient_dim())

        L = parent1.element_class(parent1, [[[0]], [], [[1]]], None)
        Q = self.product(L)
        ieqs = [F_Hrep + [width], F_Hrep + [-width]]

        H = parent2.element_class(parent2, None, [ieqs, []])
        return Q.intersection(H)

    def face_split(self, face):
        """
        Return the face splitting of the face ``face``.

        Splitting a face correspond to the bipyramid (see :meth:`bipyramid`)
        of ``self`` where the two new vertices are placed above and below
        the center of ``face`` instead of the center of the whole polyhedron.
        The two new vertices are placed in the new dimension at height `-1` and
        `1`.

        INPUT:

        - ``face`` -- a PolyhedronFace or a Vertex

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: pentagon  = polytopes.regular_polygon(5)
            sage: f = pentagon.faces(1)[0]
            sage: fsplit_pentagon = pentagon.face_split(f)
            sage: fsplit_pentagon.f_vector()
            (1, 7, 14, 9, 1)

        TESTS:

        Check that :issue:`28668` is fixed::

            sage: P = polytopes.octahedron()
            sage: P.face_split(P.faces(2)[0])
            A 4-dimensional polyhedron in QQ^4 defined as the convex hull of 8 vertices

        .. SEEALSO::

            :meth:`one_point_suspension`
        """
        from sage.geometry.polyhedron.representation import Vertex
        from sage.geometry.polyhedron.face import PolyhedronFace
        if isinstance(face, Vertex):
            new_vertices = [list(x) + [0] for x in self.vertex_generator()] + \
                           [list(face) + [x] for x in [-1, 1]]  # Splitting the vertex
        elif isinstance(face, PolyhedronFace):
            new_vertices = [list(x) + [0] for x in self.vertex_generator()] + \
                           [list(face.as_polyhedron().center()) + [x] for x in [-1, 1]]  # Splitting the face
        else:
            raise TypeError("the face {} should be a Vertex or PolyhedronFace".format(face))

        new_rays = []
        new_rays.extend( [ r + [0] for r in self.ray_generator() ] )

        new_lines = []
        new_lines.extend( [ l + [0] for l in self.line_generator() ] )

        parent = self.parent().change_ring(self.base_ring().fraction_field(), ambient_dim=self.ambient_dim()+1)
        return parent.element_class(parent, [new_vertices, new_rays, new_lines], None)

    ###########################################################
    # Methods using a vertex or vector.
    ###########################################################

    def lawrence_extension(self, v):
        """
        Return the Lawrence extension of ``self`` on the point ``v``.

        Let `P` be a polytope and `v` be a vertex of `P` or a point outside
        `P`. The Lawrence extension of `P` on `v` is the convex hull of
        `(v,1),(v,2)` and `(u,0)` for all vertices `u` in `P` other than `v`
        if `v` is a vertex.

        INPUT:

        - ``v`` -- a vertex of ``self`` or a point outside it

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: P.lawrence_extension(P.vertices()[0])
            A 4-dimensional polyhedron in ZZ^4 defined as the convex hull of 9 vertices
            sage: P.lawrence_extension([-1,-1,-1])
            A 4-dimensional polyhedron in ZZ^4 defined as the convex hull of 9 vertices

        REFERENCES:

            For more information, see Section 6.6 of [Zie2007]_.
        """
        if not self.is_compact():
            raise NotImplementedError("self must be a polytope")

        V = self.vertices_list()
        v = list(v)

        if self.contains(v) and (v not in V):
            raise ValueError("{} must not be a vertex or outside self".format(v))

        lambda_V = [u + [0] for u in V if u != v] + [v+[1]] + [v+[2]]
        parent = self.parent().base_extend(vector(v), ambient_dim=self.ambient_dim() + 1)
        return parent.element_class(parent, [lambda_V, [], []], None)

    def _test_lawrence(self, tester=None, **options):
        """
        Run tests on the methods related to lawrence extensions.

        TESTS:

        Check that :issue:`28725` is fixed::

            sage: polytopes.regular_polygon(3)._test_lawrence()                         # needs sage.rings.number_field

        Check that :issue:`30293` is fixed::

            sage: polytopes.cube()._test_lawrence()
        """
        if tester is None:
            tester = self._tester(**options)

        if self.backend() == 'normaliz' and self.base_ring() not in (ZZ, QQ):
            # Speeds up the doctest for significantly.
            self = self.change_ring(self._internal_base_ring)

        if not self.is_compact():
            with tester.assertRaises(NotImplementedError):
                self.lawrence_polytope()
            with tester.assertRaises(NotImplementedError):
                self.lawrence_extension(self.vertices()[0])
            return

        if self.n_vertices() > 1:
            # ``v`` must be a vertex or outside ``self``.
            with tester.assertRaises(ValueError):
                self.lawrence_extension(self.center())

        if self.n_vertices() >= 40 or self.n_facets() > 40:
            # Avoid very long tests.
            return

        if self.n_vertices():
            from sage.misc.prandom import randint
            v = self.vertices()[randint(0, self.n_vertices()-1)].vector()

            # A lawrence extension with a vertex.
            P = self.lawrence_extension(v)
            tester.assertEqual(self.dim() + 1, P.dim())
            tester.assertEqual(self.n_vertices() + 1, P.n_vertices())
            tester.assertEqual(self.backend(), P.backend())

            if self.n_vertices() > 1:
                # A lawrence extension with a point outside of the polyhedron.
                Q = self.lawrence_extension(2*v - self.center())
                tester.assertEqual(self.dim() + 1, Q.dim())
                tester.assertEqual(self.n_vertices() + 2, Q.n_vertices())
                tester.assertEqual(self.backend(), Q.backend())  # Any backend should handle the fraction field.

                import warnings

                with warnings.catch_warnings():
                    warnings.simplefilter("error")
                    try:
                        from sage.rings.real_double_field import RDF
                        two = RDF(2.0)
                        # Implicitly checks :issue:`30328`.
                        R = self.lawrence_extension(two * v - self.center())
                        tester.assertEqual(self.dim() + 1, R.dim())
                        tester.assertEqual(self.n_vertices() + 2, R.n_vertices())

                        tester.assertTrue(Q.is_combinatorially_isomorphic(R))
                    except ImportError:
                        # RDF not available
                        pass
                    except UserWarning:
                        # Data is numerically complicated.
                        pass
                    except ValueError as err:
                        if "Numerical inconsistency" not in err.args[0]:
                            raise err

        if self.n_vertices() >= 12 or (self.base_ring() not in (ZZ, QQ) and self.backend() == 'field'):
            # Avoid very long tests.
            return

        P = self.lawrence_polytope()
        tester.assertEqual(self.dim() + self.n_vertices(), P.dim())
        tester.assertEqual(self.n_vertices()*2, P.n_vertices())
        tester.assertEqual(self.backend(), P.backend())
        tester.assertTrue(P.is_lawrence_polytope())

        # Construct the lawrence polytope iteratively by lawrence extensions.
        V = self.vertices_list()
        Q = self
        i = 0
        for v in V:
            v = v + i*[0]
            Q = Q.lawrence_extension(v)
            i = i + 1
        tester.assertEqual(P, Q)

    def one_point_suspension(self, vertex):
        """
        Return the one-point suspension of ``self`` by splitting the vertex
        ``vertex``.

        The resulting polyhedron has one more vertex and its dimension
        increases by one.

        INPUT:

        - ``vertex`` -- a Vertex of ``self``

        EXAMPLES::

            sage: cube = polytopes.cube()
            sage: v = cube.vertices()[0]
            sage: ops_cube = cube.one_point_suspension(v)
            sage: ops_cube.f_vector()
            (1, 9, 24, 24, 9, 1)

            sage: # needs sage.rings.number_field
            sage: pentagon  = polytopes.regular_polygon(5)
            sage: v = pentagon.vertices()[0]
            sage: ops_pentagon = pentagon.one_point_suspension(v)
            sage: ops_pentagon.f_vector()
            (1, 6, 12, 8, 1)

        It works with a polyhedral face as well::

            sage: vv = cube.faces(0)[1]
            sage: ops_cube2 = cube.one_point_suspension(vv)
            sage: ops_cube == ops_cube2
            True

        .. SEEALSO::

            :meth:`face_split`

        TESTS::

            sage: e = cube.faces(1)[0]
            sage: cube.one_point_suspension(e)
            Traceback (most recent call last):
            ...
            TypeError: the vertex
            A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices
            should be a Vertex or PolyhedronFace of dimension 0
        """
        from sage.geometry.polyhedron.representation import Vertex
        from sage.geometry.polyhedron.face import PolyhedronFace
        if isinstance(vertex, Vertex):
            return self.face_split(vertex)
        elif isinstance(vertex, PolyhedronFace) and vertex.dim() == 0:
            return self.face_split(vertex)
        else:
            raise TypeError("the vertex {} should be a Vertex or PolyhedronFace of dimension 0".format(vertex))
