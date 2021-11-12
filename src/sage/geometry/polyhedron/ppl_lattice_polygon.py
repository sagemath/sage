"""
Fast Lattice Polygons using PPL

See :mod:`ppl_lattice_polytope` for the implementation of
arbitrary-dimensional lattice polytopes. This module is about the
specialization to 2 dimensions. To be more precise, the
:class:`LatticePolygon_PPL_class` is used if the ambient space is of
dimension 2 or less. These all allow you to cyclically order (see
:meth:`LatticePolygon_PPL_class.ordered_vertices`) the vertices, which
is in general not possible in higher dimensions.
"""

########################################################################
#       Copyright (C) 2012 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################

from sage.rings.integer_ring import ZZ
from sage.misc.cachefunc import cached_method, cached_function
from sage.modules.free_module_element import vector, zero_vector
from sage.matrix.constructor import (matrix, zero_matrix, block_matrix)
from ppl import C_Polyhedron, Poly_Con_Relation
from sage.geometry.polyhedron.lattice_euclidean_group_element import (
    LatticeEuclideanGroupElement)
from sage.geometry.polyhedron.ppl_lattice_polytope import (
    LatticePolytope_PPL, LatticePolytope_PPL_class)


########################################################################
class LatticePolygon_PPL_class(LatticePolytope_PPL_class):
    """
    A lattice polygon

    This includes 2-dimensional polytopes as well as degenerate (0 and
    1-dimensional) lattice polygons. Any polytope in 2d is a polygon.
    """

    @cached_method
    def ordered_vertices(self):
        """
        Return the vertices of a lattice polygon in cyclic order.

        OUTPUT:

        A tuple of vertices ordered along the perimeter of the
        polygon. The first point is arbitrary.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: square = LatticePolytope_PPL((0,0), (1,1), (0,1), (1,0))
            sage: square.vertices()
            ((0, 0), (0, 1), (1, 0), (1, 1))
            sage: square.ordered_vertices()
            ((0, 0), (1, 0), (1, 1), (0, 1))
        """
        neighbors = dict()
        if self.affine_dimension() < 2:
            return self.vertices()
        for c in self.minimized_constraints():
            v1, v2 = self.vertices_saturating(c)
            neighbors[v1] = [v2] + neighbors.get(v1, [])
            neighbors[v2] = [v1] + neighbors.get(v2, [])
        v_prev = self.vertices()[0]
        v_curr = neighbors[v_prev][0]
        result = [v_prev, v_curr]
        while len(result) < self.n_vertices():
            v1, v2 = neighbors[v_curr]
            if v1 == v_prev:
                v_next = v2
            else:
                v_next = v1
            result.append(v_next)
            v_prev = v_curr
            v_curr = v_next
        return tuple(result)

    def _find_isomorphism_degenerate(self, polytope):
        """
        Helper to pick an isomorphism of degenerate polygons

        INPUT:

        - ``polytope`` -- a :class:`LatticePolytope_PPL_class`. The
          polytope to compare with.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL, C_Polyhedron
            sage: L1 = LatticePolytope_PPL(C_Polyhedron(2, 'empty'))
            sage: L2 = LatticePolytope_PPL(C_Polyhedron(3, 'empty'))
            sage: iso = L1.find_isomorphism(L2)   # indirect doctest
            sage: iso(L1) == L2
            True
            sage: iso = L1._find_isomorphism_degenerate(L2)
            sage: iso(L1) == L2
            True

            sage: L1 = LatticePolytope_PPL((-1,4))
            sage: L2 = LatticePolytope_PPL((2,1,5))
            sage: iso = L1.find_isomorphism(L2)
            sage: iso(L1) == L2
            True

            sage: L1 = LatticePolytope_PPL((-1,), (3,))
            sage: L2 = LatticePolytope_PPL((2,1,5), (2,-3,5))
            sage: iso = L1.find_isomorphism(L2)
            sage: iso(L1) == L2
            True

            sage: L1 = LatticePolytope_PPL((-1,-1), (3,-1))
            sage: L2 = LatticePolytope_PPL((2,1,5), (2,-3,5))
            sage: iso = L1.find_isomorphism(L2)
            sage: iso(L1) == L2
            True

            sage: L1 = LatticePolytope_PPL((-1,2), (3,1))
            sage: L2 = LatticePolytope_PPL((1,2,3),(1,2,4))
            sage: iso = L1.find_isomorphism(L2)
            sage: iso(L1) == L2
            True

            sage: L1 = LatticePolytope_PPL((-1,2), (3,2))
            sage: L2 = LatticePolytope_PPL((1,2,3),(1,2,4))
            sage: L1.find_isomorphism(L2)
            Traceback (most recent call last):
            ...
            LatticePolytopesNotIsomorphicError: different number of integral points

            sage: L1 = LatticePolytope_PPL((-1,2), (3,1))
            sage: L2 = LatticePolytope_PPL((1,2,3),(1,2,5))
            sage: L1.find_isomorphism(L2)
            Traceback (most recent call last):
            ...
            LatticePolytopesNotIsomorphicError: different number of integral points
        """
        from sage.geometry.polyhedron.lattice_euclidean_group_element import \
            LatticePolytopesNotIsomorphicError
        polytope_vertices = polytope.vertices()
        self_vertices = self.ordered_vertices()
        # handle degenerate cases
        if self.n_vertices() == 0:
            A = zero_matrix(ZZ, polytope.space_dimension(), self.space_dimension())
            b = zero_vector(ZZ, polytope.space_dimension())
            return LatticeEuclideanGroupElement(A, b)
        if self.n_vertices() == 1:
            A = zero_matrix(ZZ, polytope.space_dimension(), self.space_dimension())
            b = polytope_vertices[0]
            return LatticeEuclideanGroupElement(A, b)
        if self.n_vertices() == 2:
            self_origin = self_vertices[0]
            self_ray = self_vertices[1] - self_origin
            polytope_origin = polytope_vertices[0]
            polytope_ray = polytope_vertices[1] - polytope_origin
            Ds, Us, Vs = self_ray.column().smith_form()
            Dp, Up, Vp = polytope_ray.column().smith_form()
            assert Vs.nrows() == Vs.ncols() == Vp.nrows() == Vp.ncols() == 1
            assert abs(Vs[0, 0]) == abs(Vp[0, 0]) == 1
            A = zero_matrix(ZZ, Dp.nrows(), Ds.nrows())
            A[0, 0] = 1
            A = Up.inverse() * A * Us * (Vs[0, 0] * Vp[0, 0])
            b = polytope_origin - A*self_origin
            try:
                A = matrix(ZZ, A)
                b = vector(ZZ, b)
            except TypeError:
                raise LatticePolytopesNotIsomorphicError('different lattice')
            hom = LatticeEuclideanGroupElement(A, b)
            if hom(self) == polytope:
                return hom
            raise LatticePolytopesNotIsomorphicError('different polygons')

    def _find_cyclic_isomorphism_matching_edge(self, polytope,
                                               polytope_origin, p_ray_left,
                                               p_ray_right):
        r"""
        Helper to find an isomorphism of polygons

        INPUT:

        - ``polytope`` -- the lattice polytope to compare to.

        - ``polytope_origin`` -- `\ZZ`-vector. a vertex of ``polytope``

        - ``p_ray_left`` - vector. the vector from ``polytope_origin``
          to one of its neighboring vertices.

        - ``p_ray_right`` - vector. the vector from
          ``polytope_origin`` to the other neighboring vertices.

        OUTPUT:

        The element of the lattice Euclidean group that maps ``self``
        to ``polytope`` with given origin and left/right neighboring
        vertex. A
        :class:`~sage.geometry.polyhedron.lattice_euclidean_group_element.LatticePolytopesNotIsomorphicError`
        is raised if no such isomorphism exists.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: L1 = LatticePolytope_PPL((1,0),(0,1),(0,0))
            sage: L2 = LatticePolytope_PPL((1,0,3),(0,1,0),(0,0,1))
            sage: v0, v1, v2 = L2.vertices()
            sage: L1._find_cyclic_isomorphism_matching_edge(L2, v0, v1-v0, v2-v0)
            The map A*x+b with A=
            [ 0  1]
            [-1 -1]
            [ 1  3]
            b =
            (0, 1, 0)
        """
        from sage.geometry.polyhedron.lattice_euclidean_group_element import \
            LatticePolytopesNotIsomorphicError
        polytope_matrix = block_matrix(1, 2, [p_ray_left.column(),
                                              p_ray_right.column()])
        self_vertices = self.ordered_vertices()
        for i in range(len(self_vertices)):
            # three consecutive vertices
            v_left = self_vertices[(i+0) % len(self_vertices)]
            v_origin = self_vertices[(i+1) % len(self_vertices)]
            v_right = self_vertices[(i+2) % len(self_vertices)]
            r_left = v_left-v_origin
            r_right = v_right-v_origin
            self_matrix = block_matrix(1, 2, [r_left.column(),
                                              r_right.column()])
            A = self_matrix.solve_left(polytope_matrix)
            b = polytope_origin - A*v_origin
            try:
                A = matrix(ZZ, A)
                b = vector(ZZ, b)
            except TypeError:
                continue
            if A.elementary_divisors()[0:2] != [1, 1]:
                continue
            hom = LatticeEuclideanGroupElement(A, b)
            if hom(self) == polytope:
                return hom
        raise LatticePolytopesNotIsomorphicError('different polygons')

    def find_isomorphism(self, polytope):
        r"""
        Return a lattice isomorphism with ``polytope``.

        INPUT:

        - ``polytope`` -- a polytope, potentially higher-dimensional.

        OUTPUT:

        A
        :class:`~sage.geometry.polyhedron.lattice_euclidean_group_element.LatticeEuclideanGroupElement`. It
        is not necessarily invertible if the affine dimension of
        ``self`` or ``polytope`` is not two. A
        :class:`~sage.geometry.polyhedron.lattice_euclidean_group_element.LatticePolytopesNotIsomorphicError`
        is raised if no such isomorphism exists.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: L1 = LatticePolytope_PPL((1,0),(0,1),(0,0))
            sage: L2 = LatticePolytope_PPL((1,0,3),(0,1,0),(0,0,1))
            sage: iso = L1.find_isomorphism(L2)
            sage: iso(L1) == L2
            True

            sage: L1 = LatticePolytope_PPL((0, 1), (3, 0), (0, 3), (1, 0))
            sage: L2 = LatticePolytope_PPL((0,0,2,1),(0,1,2,0),(2,0,0,3),(2,3,0,0))
            sage: iso = L1.find_isomorphism(L2)
            sage: iso(L1) == L2
            True

        The following polygons are isomorphic over `\QQ`, but not as
        lattice polytopes::

            sage: L1 = LatticePolytope_PPL((1,0),(0,1),(-1,-1))
            sage: L2 = LatticePolytope_PPL((0, 0), (0, 1), (1, 0))
            sage: L1.find_isomorphism(L2)
            Traceback (most recent call last):
            ...
            LatticePolytopesNotIsomorphicError: different number of integral points
            sage: L2.find_isomorphism(L1)
            Traceback (most recent call last):
            ...
            LatticePolytopesNotIsomorphicError: different number of integral points
        """
        from sage.geometry.polyhedron.lattice_euclidean_group_element import \
            LatticePolytopesNotIsomorphicError
        if polytope.affine_dimension() != self.affine_dimension():
            raise LatticePolytopesNotIsomorphicError('different dimension')
        polytope_vertices = polytope.vertices()
        if len(polytope_vertices) != self.n_vertices():
            raise LatticePolytopesNotIsomorphicError('different number of vertices')
        self_vertices = self.ordered_vertices()
        if len(polytope.integral_points()) != len(self.integral_points()):
            raise LatticePolytopesNotIsomorphicError('different number of integral points')

        if len(self_vertices) < 3:
            return self._find_isomorphism_degenerate(polytope)

        polytope_origin = polytope_vertices[0]
        origin_P = C_Polyhedron(next(iter(polytope.minimized_generators())))

        neighbors = []
        for c in polytope.minimized_constraints():
            if not c.is_inequality():
                continue
            if origin_P.relation_with(c).implies(Poly_Con_Relation.saturates()):
                for i, g in enumerate(polytope.minimized_generators()):
                    if i == 0:
                        continue
                    g = C_Polyhedron(g)
                    if g.relation_with(c).implies(Poly_Con_Relation.saturates()):
                        neighbors.append(polytope_vertices[i])
                        break

        p_ray_left = neighbors[0] - polytope_origin
        p_ray_right = neighbors[1] - polytope_origin
        try:
            return self._find_cyclic_isomorphism_matching_edge(polytope, polytope_origin,
                                                               p_ray_left, p_ray_right)
        except LatticePolytopesNotIsomorphicError:
            pass
        try:
            return self._find_cyclic_isomorphism_matching_edge(polytope, polytope_origin,
                                                               p_ray_right, p_ray_left)
        except LatticePolytopesNotIsomorphicError:
            pass
        raise LatticePolytopesNotIsomorphicError('different polygons')

    def is_isomorphic(self, polytope):
        """
        Test if ``self`` and ``polytope`` are isomorphic.

        INPUT:

        - ``polytope`` -- a lattice polytope.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: L1 = LatticePolytope_PPL((1,0),(0,1),(0,0))
            sage: L2 = LatticePolytope_PPL((1,0,3),(0,1,0),(0,0,1))
            sage: L1.is_isomorphic(L2)
            True
        """
        from sage.geometry.polyhedron.lattice_euclidean_group_element import \
            LatticePolytopesNotIsomorphicError
        try:
            self.find_isomorphism(polytope)
            return True
        except LatticePolytopesNotIsomorphicError:
            return False

    def sub_polytopes(self):
        """
        Return a list of all lattice sub-polygons up to isomorphism.

        OUTPUT:

        All non-empty sub-lattice polytopes up to isomorphism. This
        includes ``self`` as improper sub-polytope, but excludes the
        empty polytope. Isomorphic sub-polytopes that can be embedded
        in different places are only returned once.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: P1xP1 = LatticePolytope_PPL((1,0), (0,1), (-1,0), (0,-1))
            sage: P1xP1.sub_polytopes()
            (A 2-dimensional lattice polytope in ZZ^2 with 4 vertices,
             A 2-dimensional lattice polytope in ZZ^2 with 3 vertices,
             A 2-dimensional lattice polytope in ZZ^2 with 3 vertices,
             A 1-dimensional lattice polytope in ZZ^2 with 2 vertices,
             A 1-dimensional lattice polytope in ZZ^2 with 2 vertices,
             A 0-dimensional lattice polytope in ZZ^2 with 1 vertex)
        """
        subpolytopes = [self]
        todo = list(subpolytopes)
        while todo:
            polytope = todo.pop()
            for p in polytope.sub_polytope_generator():
                if p.is_empty():
                    continue
                if any(p.is_isomorphic(q) for q in subpolytopes):
                    continue
                subpolytopes.append(p)
                todo.append(p)
        return tuple(subpolytopes)

    def plot(self):
        """
        Plot the lattice polygon.

        OUTPUT:

        A graphics object.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: P = LatticePolytope_PPL((1,0), (0,1), (0,0), (2,2))
            sage: P.plot()  # optional - sage.plot
            Graphics object consisting of 6 graphics primitives
            sage: LatticePolytope_PPL([0], [1]).plot()  # optional - sage.plot
            Graphics object consisting of 3 graphics primitives
            sage: LatticePolytope_PPL([0]).plot()  # optional - sage.plot
            Graphics object consisting of 2 graphics primitives
        """
        from sage.plot.point import point2d
        from sage.plot.polygon import polygon2d
        vertices = self.ordered_vertices()
        points = self.integral_points()
        if self.space_dimension() == 1:
            vertices = [vector(ZZ, (v[0], 0)) for v in vertices]
            points = [vector(ZZ, (p[0], 0)) for p in points]
        point_plot = sum(point2d(p, pointsize=100, color='red')
                         for p in points)
        polygon_plot = polygon2d(vertices, alpha=0.2, color='green',
                                 zorder=-1, thickness=2)
        return polygon_plot + point_plot


########################################################################
#
#  Reflexive lattice polygons and their subpolygons
#
########################################################################

@cached_function
def polar_P2_polytope():
    """
    The polar of the `P^2` polytope

    EXAMPLES::

        sage: from sage.geometry.polyhedron.ppl_lattice_polygon import polar_P2_polytope
        sage: polar_P2_polytope()
        A 2-dimensional lattice polytope in ZZ^2 with 3 vertices
        sage: _.vertices()
        ((0, 0), (0, 3), (3, 0))
    """
    return LatticePolytope_PPL((0, 0), (3, 0), (0, 3))


@cached_function
def polar_P1xP1_polytope():
    r"""
    The polar of the `P^1 \times P^1` polytope

    EXAMPLES::

        sage: from sage.geometry.polyhedron.ppl_lattice_polygon import polar_P1xP1_polytope
        sage: polar_P1xP1_polytope()
        A 2-dimensional lattice polytope in ZZ^2 with 4 vertices
        sage: _.vertices()
        ((0, 0), (0, 2), (2, 0), (2, 2))
    """
    return LatticePolytope_PPL((0, 0), (2, 0), (0, 2), (2, 2))


@cached_function
def polar_P2_112_polytope():
    """
    The polar of the `P^2[1,1,2]` polytope

    EXAMPLES::

        sage: from sage.geometry.polyhedron.ppl_lattice_polygon import polar_P2_112_polytope
        sage: polar_P2_112_polytope()
        A 2-dimensional lattice polytope in ZZ^2 with 3 vertices
        sage: _.vertices()
        ((0, 0), (0, 2), (4, 0))
    """
    return LatticePolytope_PPL((0, 0), (4, 0), (0, 2))


@cached_function
def subpolygons_of_polar_P2():
    """
    The lattice sub-polygons of the polar `P^2` polytope

    OUTPUT:

    A tuple of lattice polytopes.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.ppl_lattice_polygon import subpolygons_of_polar_P2
        sage: len(subpolygons_of_polar_P2())
        27
    """
    return polar_P2_polytope().sub_polytopes()


@cached_function
def subpolygons_of_polar_P2_112():
    """
    The lattice sub-polygons of the polar `P^2[1,1,2]` polytope

    OUTPUT:

    A tuple of lattice polytopes.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.ppl_lattice_polygon import subpolygons_of_polar_P2_112
        sage: len(subpolygons_of_polar_P2_112())
        28
    """
    return polar_P2_112_polytope().sub_polytopes()


@cached_function
def subpolygons_of_polar_P1xP1():
    r"""
    The lattice sub-polygons of the polar `P^1 \times P^1` polytope

    OUTPUT:

    A tuple of lattice polytopes.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.ppl_lattice_polygon import subpolygons_of_polar_P1xP1
        sage: len(subpolygons_of_polar_P1xP1())
        20
    """
    return polar_P1xP1_polytope().sub_polytopes()


@cached_function
def sub_reflexive_polygons():
    """
    Return all lattice sub-polygons of reflexive polygons.

    OUTPUT:

    A tuple of all lattice sub-polygons. Each sub-polygon is returned
    as a pair sub-polygon, containing reflexive polygon.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.ppl_lattice_polygon import sub_reflexive_polygons
        sage: l = sub_reflexive_polygons(); l[5]
        (A 2-dimensional lattice polytope in ZZ^2 with 6 vertices,
         A 2-dimensional lattice polytope in ZZ^2 with 3 vertices)
        sage: len(l)
        33
    """
    result = []

    def add_result(subpolygon, ambient):
        if not any(subpolygon.is_isomorphic(p[0]) for p in result):
            result.append((subpolygon, ambient))
    for p in subpolygons_of_polar_P2():
        add_result(p, polar_P2_polytope())
    for p in subpolygons_of_polar_P2_112():
        add_result(p, polar_P2_112_polytope())
    for p in subpolygons_of_polar_P1xP1():
        add_result(p, polar_P1xP1_polytope())
    return tuple(result)
