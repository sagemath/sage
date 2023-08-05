r"""
Fast Lattice Polytopes using PPL.

The :func:`LatticePolytope_PPL` class is a thin wrapper around PPL
polyhedra. Its main purpose is to be fast to construct, at the cost of
being much less full-featured than the usual polyhedra. This makes it
possible to iterate with it over the list of all 473800776 reflexive
polytopes in 4 dimensions.

.. NOTE::

    For general lattice polyhedra you should use
    :func:`~sage.geometry.polyhedron.constructor.Polyhedron` with
    ``base_ring=ZZ``.

The class derives from the PPL :class:`ppl.polyhedron.C_Polyhedron`
class, so you can work with the underlying generator and constraint
objects. However, integral points are generally represented by
`\ZZ`-vectors. In the following, we always use *generator* to refer
the PPL generator objects and *vertex* (or integral point) for the
corresponding `\ZZ`-vector.

EXAMPLES::

    sage: vertices = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1), (-9, -6, -1, -1)]
    sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
    sage: P = LatticePolytope_PPL(vertices);  P
    A 4-dimensional lattice polytope in ZZ^4 with 5 vertices
    sage: P.integral_points()
    ((-9, -6, -1, -1), (-3, -2, 0, 0), (-2, -1, 0, 0), (-1, -1, 0, 0),
     (-1, 0, 0, 0), (0, 0, 0, 0), (1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 0, 1), (0, 0, 1, 0))
    sage: P.integral_points_not_interior_to_facets()
    ((-9, -6, -1, -1), (-3, -2, 0, 0), (0, 0, 0, 0), (1, 0, 0, 0),
     (0, 1, 0, 0), (0, 0, 0, 1), (0, 0, 1, 0))

Fibrations of the lattice polytopes are defined as lattice
sub-polytopes and give rise to fibrations of toric varieties for
suitable fan refinements. We can compute them using
:meth:`~LatticePolytope_PPL.fibration_generator` ::

    sage: F = next(P.fibration_generator(2))
    sage: F.vertices()
    ((1, 0, 0, 0), (0, 1, 0, 0), (-3, -2, 0, 0))

Finally, we can compute automorphisms and identify fibrations that
only differ by a lattice automorphism::

    sage: square = LatticePolytope_PPL((-1,-1), (-1,1), (1,-1), (1,1))
    sage: fibers = [ f.vertices() for f in square.fibration_generator(1) ];  fibers
    [((1, 0), (-1, 0)), ((0, 1), (0, -1)), ((-1, -1), (1, 1)), ((-1, 1), (1, -1))]
    sage: square.pointsets_mod_automorphism(fibers)                                     # needs sage.groups
    (frozenset({(-1, -1), (1, 1)}), frozenset({(-1, 0), (1, 0)}))

AUTHORS:

    - Volker Braun: initial version, 2012
"""

# ****************************************************************************
#       Copyright (C) 2012 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import copy
from sage.rings.integer import GCD_list, Integer
from sage.rings.integer_ring import ZZ
from sage.misc.cachefunc import cached_method
from sage.modules.free_module_element import vector
from sage.matrix.constructor import matrix
from ppl import (
    C_Polyhedron, Linear_Expression, Variable,
    point, line, Generator, Generator_System,
    Poly_Con_Relation )


########################################################################
def _class_for_LatticePolytope(dim):
    """
    Return the appropriate class in the given dimension.

    Helper function for :func:`LatticePolytope_PPL`. You should not
    have to use this function manually.

    INPUT:

    - ``dim`` -- integer. The ambient space dimension.

    OUTPUT:

    The appropriate class for the lattice polytope.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.ppl_lattice_polytope import _class_for_LatticePolytope
        sage: _class_for_LatticePolytope(2)
        <class 'sage.geometry.polyhedron.ppl_lattice_polygon.LatticePolygon_PPL_class'>
        sage: _class_for_LatticePolytope(3)
        <class 'sage.geometry.polyhedron.ppl_lattice_polytope.LatticePolytope_PPL_class'>
    """
    if dim <= 2:
        from sage.geometry.polyhedron.ppl_lattice_polygon import LatticePolygon_PPL_class
        return LatticePolygon_PPL_class
    return LatticePolytope_PPL_class


########################################################################
def LatticePolytope_PPL(*args):
    """
    Construct a new instance of the PPL-based lattice polytope class.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
        sage: LatticePolytope_PPL((0,0), (1,0), (0,1))
        A 2-dimensional lattice polytope in ZZ^2 with 3 vertices

        sage: from ppl import point, Generator_System, C_Polyhedron, Linear_Expression  # needs pplpy
        sage: p = point(Linear_Expression([2,3],0));  p                                 # needs pplpy
        point(2/1, 3/1)
        sage: LatticePolytope_PPL(p)                                                    # needs pplpy
        A 0-dimensional lattice polytope in ZZ^2 with 1 vertex

        sage: P = C_Polyhedron(Generator_System(p));  P                                 # needs pplpy
        A 0-dimensional polyhedron in QQ^2 defined as the convex hull of 1 point
        sage: LatticePolytope_PPL(P)                                                    # needs pplpy
        A 0-dimensional lattice polytope in ZZ^2 with 1 vertex

    A ``TypeError`` is raised if the arguments do not specify a lattice polytope::

        sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
        sage: LatticePolytope_PPL((0,0), (1/2,1))                                       # needs pplpy
        Traceback (most recent call last):
        ...
        TypeError: unable to convert rational 1/2 to an integer

        sage: from ppl import point, Generator_System, C_Polyhedron, Linear_Expression  # needs pplpy
        sage: p = point(Linear_Expression([2,3],0), 5);  p                              # needs pplpy
        point(2/5, 3/5)
        sage: LatticePolytope_PPL(p)                                                    # needs pplpy
        Traceback (most recent call last):
         ...
        TypeError: generator is not a lattice polytope generator

        sage: P = C_Polyhedron(Generator_System(p));  P                                 # needs pplpy
        A 0-dimensional polyhedron in QQ^2 defined as the convex hull of 1 point
        sage: LatticePolytope_PPL(P)                                                    # needs pplpy
        Traceback (most recent call last):
        ...
        TypeError: polyhedron has non-integral generators
    """
    polytope_class = LatticePolytope_PPL_class
    if len(args)==1 and isinstance(args[0], C_Polyhedron):
        polyhedron = args[0]
        polytope_class = _class_for_LatticePolytope(polyhedron.space_dimension())
        if not all(p.is_point() and p.divisor() == 1 for p in polyhedron.generators()):
            raise TypeError('polyhedron has non-integral generators')
        return polytope_class(polyhedron)
    if len(args)==1 \
            and isinstance(args[0], (list, tuple)) \
            and isinstance(args[0][0], (list,tuple)):
        vertices = args[0]
    else:
        vertices = args
    gs = Generator_System()
    for v in vertices:
        if isinstance(v, Generator):
            if (not v.is_point()) or v.divisor() != 1:
                raise TypeError('generator is not a lattice polytope generator')
            gs.insert(v)
        else:
            gs.insert(point(Linear_Expression(v, 0)))
    if not gs.empty():
        dim = next(iter(gs)).space_dimension()
        polytope_class = _class_for_LatticePolytope(dim)
    return polytope_class(gs)


########################################################################
class LatticePolytope_PPL_class(C_Polyhedron):
    """
    The lattice polytope class.

    You should use :func:`LatticePolytope_PPL` to construct instances.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
        sage: LatticePolytope_PPL((0,0), (1,0), (0,1))
        A 2-dimensional lattice polytope in ZZ^2 with 3 vertices
    """

    def __repr__(self):
        """
        Return the string representation

        OUTPUT:

        A string.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: P = LatticePolytope_PPL((0,0), (1,0), (0,1))
            sage: P
            A 2-dimensional lattice polytope in ZZ^2 with 3 vertices
            sage: P.__repr__()
            'A 2-dimensional lattice polytope in ZZ^2 with 3 vertices'

            sage: LatticePolytope_PPL()
            The empty lattice polytope in ZZ^0
        """
        desc = ''
        if self.n_vertices() == 0:
            desc += 'The empty lattice polytope'
        else:
            desc += 'A ' + repr(self.affine_dimension()) + '-dimensional lattice polytope'
        desc += ' in ZZ^' + repr(self.space_dimension())

        if self.n_vertices() > 0:
            desc += ' with '
            desc += repr(self.n_vertices())
            if self.n_vertices() == 1:
                desc += ' vertex'
            else:
                desc += ' vertices'
        return desc

    def is_bounded(self):
        """
        Return whether the lattice polytope is compact.

        OUTPUT:

        Always ``True``, since polytopes are by definition compact.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: LatticePolytope_PPL((0,0), (1,0), (0,1)).is_bounded()
            True
        """
        return True

    @cached_method
    def n_vertices(self):
        """
        Return the number of vertices.

        OUTPUT:

        An integer, the number of vertices.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: LatticePolytope_PPL((0,0,0), (1,0,0), (0,1,0)).n_vertices()
            3
        """
        return len(self.minimized_generators())

    @cached_method
    def is_simplex(self):
        r"""
        Return whether the polyhedron is a simplex.

        OUTPUT:

        Boolean, whether the polyhedron is a simplex (possibly of
        strictly smaller dimension than the ambient space).

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: LatticePolytope_PPL((0,0,0), (1,0,0), (0,1,0)).is_simplex()
            True
        """
        return self.affine_dimension()+1==self.n_vertices()

    @cached_method
    def bounding_box(self):
        r"""
        Return the coordinates of a rectangular box containing the non-empty polytope.

        OUTPUT:

        A pair of tuples ``(box_min, box_max)`` where ``box_min`` are
        the coordinates of a point bounding the coordinates of the
        polytope from below and ``box_max`` bounds the coordinates
        from above.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: LatticePolytope_PPL((0,0), (1,0), (0,1)).bounding_box()
            ((0, 0), (1, 1))
        """
        box_min = []
        box_max = []
        if self.is_empty():
            raise ValueError('empty polytope is not allowed')
        for i in range(self.space_dimension()):
            x = Variable(i)
            coords = [Integer(v.coefficient(x)) for v in self.generators()]
            max_coord = max(coords)
            min_coord = min(coords)
            box_max.append(max_coord)
            box_min.append(min_coord)
        return (tuple(box_min), tuple(box_max))

    @cached_method
    def n_integral_points(self):
        """
        Return the number of integral points.

        OUTPUT:

        Integer. The number of integral points contained in the
        lattice polytope.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: LatticePolytope_PPL((0,0), (1,0), (0,1)).n_integral_points()
            3
        """
        if self.is_empty():
            return tuple()
        box_min, box_max = self.bounding_box()
        from sage.geometry.integral_points import rectangular_box_points
        return rectangular_box_points(list(box_min), list(box_max), self, count_only=True)

    @cached_method
    def integral_points(self):
        r"""
        Return the integral points in the polyhedron.

        Uses the naive algorithm (iterate over a rectangular bounding
        box).

        OUTPUT:

        The list of integral points in the polyhedron. If the
        polyhedron is not compact, a ``ValueError`` is raised.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: LatticePolytope_PPL((-1,-1), (1,0), (1,1), (0,1)).integral_points()
            ((-1, -1), (0, 0), (0, 1), (1, 0), (1, 1))

            sage: simplex = LatticePolytope_PPL((1,2,3), (2,3,7), (-2,-3,-11))
            sage: simplex.integral_points()
            ((-2, -3, -11), (0, 0, -2), (1, 2, 3), (2, 3, 7))

        The polyhedron need not be full-dimensional::

            sage: simplex = LatticePolytope_PPL((1,2,3,5), (2,3,7,5), (-2,-3,-11,5))
            sage: simplex.integral_points()
            ((-2, -3, -11, 5), (0, 0, -2, 5), (1, 2, 3, 5), (2, 3, 7, 5))

            sage: point = LatticePolytope_PPL((2,3,7))
            sage: point.integral_points()
            ((2, 3, 7),)

            sage: empty = LatticePolytope_PPL()
            sage: empty.integral_points()
            ()

        Here is a simplex where the naive algorithm of running over
        all points in a rectangular bounding box no longer works fast
        enough::

            sage: v = [(1,0,7,-1), (-2,-2,4,-3), (-1,-1,-1,4), (2,9,0,-5), (-2,-1,5,1)]
            sage: simplex = LatticePolytope_PPL(v); simplex
            A 4-dimensional lattice polytope in ZZ^4 with 5 vertices
            sage: len(simplex.integral_points())
            49

        Finally, the 3-d reflexive polytope number 4078::

            sage: v = [(1,0,0), (0,1,0), (0,0,1), (0,0,-1), (0,-2,1),
            ....:      (-1,2,-1), (-1,2,-2), (-1,1,-2), (-1,-1,2), (-1,-3,2)]
            sage: P = LatticePolytope_PPL(*v)
            sage: pts1 = P.integral_points()            # Sage's own code
            sage: pts2 = LatticePolytope(v).points()                                    # needs palp
            sage: for p in pts1: p.set_immutable()
            sage: set(pts1) == set(pts2)                                                # needs palp
            True

            sage: len(Polyhedron(v).integral_points())  # takes about 1 ms
            23
            sage: len(LatticePolytope(v).points())      # takes about 13 ms             # needs palp
            23
            sage: len(LatticePolytope_PPL(*v).integral_points())  # takes about 0.5 ms
            23
        """
        if self.is_empty():
            return tuple()
        box_min, box_max = self.bounding_box()
        from sage.geometry.integral_points import rectangular_box_points
        points = rectangular_box_points(list(box_min), list(box_max), self)
        if not self.n_integral_points.is_in_cache():
            self.n_integral_points.set_cache(len(points))
        return points

    @cached_method
    def _integral_points_saturating(self):
        """
        Return the integral points together with information about
        which inequalities are saturated.

        See :func:`~sage.geometry.integral_points.rectangular_box_points`.

        OUTPUT:

        A tuple of pairs (one for each integral point) consisting of a
        pair ``(point, Hrep)``, where ``point`` is the coordinate
        vector of the integral point and ``Hrep`` is the set of
        indices of the :meth:`minimized_constraints` that are
        saturated at the point.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: quad = LatticePolytope_PPL((-1,-1), (0,1), (1,0), (1,1))
            sage: quad._integral_points_saturating()
            (((-1, -1), frozenset({0, 1})),
             ((0, 0), frozenset()),
             ((0, 1), frozenset({0, 3})),
             ((1, 0), frozenset({1, 2})),
             ((1, 1), frozenset({2, 3})))
        """
        if self.is_empty():
            return tuple()
        box_min, box_max = self.bounding_box()
        from sage.geometry.integral_points import rectangular_box_points
        points = rectangular_box_points(list(box_min), list(box_max), self,
                                        return_saturated=True)
        if not self.n_integral_points.is_in_cache():
            self.n_integral_points.set_cache(len(points))
        if not self.integral_points.is_in_cache():
            self.integral_points.set_cache(tuple(p[0] for p in points))
        return points

    @cached_method
    def integral_points_not_interior_to_facets(self):
        """
        Return the integral points not interior to facets.

        OUTPUT:

        A tuple whose entries are the coordinate vectors of integral
        points not interior to facets (codimension one faces) of the
        lattice polytope.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: square = LatticePolytope_PPL((-1,-1), (-1,1), (1,-1), (1,1))
            sage: square.n_integral_points()
            9
            sage: square.integral_points_not_interior_to_facets()
            ((-1, -1), (-1, 1), (0, 0), (1, -1), (1, 1))
        """
        n = 1 + self.space_dimension() - self.affine_dimension()
        return tuple(p[0] for p in self._integral_points_saturating() if len(p[1])!=n)

    @cached_method
    def vertices(self):
        r"""
        Return the vertices as a tuple of `\ZZ`-vectors.

        OUTPUT:

        A tuple of `\ZZ`-vectors. Each entry is the coordinate vector
        of an integral points of the lattice polytope.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: p = LatticePolytope_PPL((-9,-6,-1,-1),
            ....:                         (0,0,0,1), (0,0,1,0), (0,1,0,0), (1,0,0,0))
            sage: p.vertices()
            ((-9, -6, -1, -1), (0, 0, 0, 1), (0, 0, 1, 0), (0, 1, 0, 0), (1, 0, 0, 0))
            sage: p.minimized_generators()
            Generator_System {point(-9/1, -6/1, -1/1, -1/1), point(0/1, 0/1, 0/1, 1/1),
            point(0/1, 0/1, 1/1, 0/1), point(0/1, 1/1, 0/1, 0/1), point(1/1, 0/1, 0/1, 0/1)}
        """
        d = self.space_dimension()
        v = vector(ZZ, d)
        points = []
        for g in self.minimized_generators():
            for i in range(d):
                v[i] = g.coefficient(Variable(i))
            v_copy = copy.copy(v)
            v_copy.set_immutable()
            points.append(v_copy)
        return tuple(points)

    def vertices_saturating(self, constraint):
        r"""
        Return the vertices saturating the constraint.

        INPUT:

        - ``constraint`` -- a constraint (inequality or equation) of
          the polytope.

        OUTPUT:

        The tuple of vertices saturating the constraint. The vertices
        are returned as `\ZZ`-vectors, as in :meth:`vertices`.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: p = LatticePolytope_PPL((0,0), (0,1), (1,0))
            sage: ieq = next(iter(p.constraints()));  ieq
            x0>=0
            sage: p.vertices_saturating(ieq)
            ((0, 0), (0, 1))
        """
        from ppl import C_Polyhedron, Poly_Con_Relation
        result = []
        for i,v in enumerate(self.minimized_generators()):
            v = C_Polyhedron(v)
            if v.relation_with(constraint).implies(Poly_Con_Relation.saturates()):
                result.append(self.vertices()[i])
        return tuple(result)

    @cached_method
    def is_full_dimensional(self):
        """
        Return whether the lattice polytope is full dimensional.

        OUTPUT:

        Boolean. Whether the :meth:`affine_dimension` equals the
        ambient space dimension.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: p = LatticePolytope_PPL((0,0), (0,1))
            sage: p.is_full_dimensional()
            False
            sage: q = LatticePolytope_PPL((0,0), (0,1), (1,0))
            sage: q.is_full_dimensional()
            True
        """

        return self.affine_dimension() == self.space_dimension()

    def fibration_generator(self, dim):
        """
        Generate the lattice polytope fibrations.

        For the purposes of this function, a lattice polytope fiber is
        a sub-lattice polytope. Projecting the plane spanned by the
        subpolytope to a point yields another lattice polytope, the
        base of the fibration.

        INPUT:

        - ``dim`` -- integer. The dimension of the lattice polytope
          fiber.

        OUTPUT:

        A generator yielding the distinct lattice polytope fibers of
        given dimension.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: p = LatticePolytope_PPL((-9,-6,-1,-1),
            ....:                         (0,0,0,1), (0,0,1,0), (0,1,0,0), (1,0,0,0))
            sage: list(p.fibration_generator(2))
            [A 2-dimensional lattice polytope in ZZ^4 with 3 vertices]
        """
        assert self.is_full_dimensional()
        # "points" are the potential vertices of the fiber. They are
        # in the $codim$-skeleton of the polytope, which is contained
        # in the points that saturate at least $dim$ equations.
        points = [ p for p in self._integral_points_saturating() if len(p[1])>=dim ]
        points = sorted(points, key=lambda x:len(x[1]))

        # iterate over point combinations subject to all points being on one facet.
        def point_combinations_iterator(n, i0=0, saturated=None):
            for i in range(i0, len(points)):
                p, ieqs = points[i]
                if saturated is None:
                    saturated_ieqs = ieqs
                else:
                    saturated_ieqs = saturated.intersection(ieqs)
                if len(saturated_ieqs)==0:
                    continue
                if n == 1:
                    yield [i]
                else:
                    for c in point_combinations_iterator(n-1, i+1, saturated_ieqs):
                        yield [i] + c

        point_lines = [ line(Linear_Expression(p[0].list(),0)) for p in points ]
        origin = point()
        fibers = set()
        gs = Generator_System()
        for indices in point_combinations_iterator(dim):
            gs.clear()
            gs.insert(origin)
            for i in indices:
                gs.insert(point_lines[i])
            plane = C_Polyhedron(gs)
            if plane.affine_dimension() != dim:
                continue
            plane.intersection_assign(self)
            if (not self.is_full_dimensional()) and (plane.affine_dimension() != dim):
                continue
            try:
                fiber = LatticePolytope_PPL(plane)
            except TypeError:   # not a lattice polytope
                continue
            fiber_vertices = tuple(sorted(fiber.vertices()))
            if fiber_vertices not in fibers:
                yield fiber
                fibers.update([fiber_vertices])

    def pointsets_mod_automorphism(self, pointsets):
        """
        Return ``pointsets`` modulo the automorphisms of ``self``.

        INPUT:

        - ``polytopes`` -- a tuple/list/iterable of subsets of the
          integral points of ``self``.

        OUTPUT:

        Representatives of the point sets modulo the
        :meth:`lattice_automorphism_group` of ``self``.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: square = LatticePolytope_PPL((-1,-1), (-1,1), (1,-1), (1,1))
            sage: fibers = [f.vertices() for f in square.fibration_generator(1)]
            sage: square.pointsets_mod_automorphism(fibers)                             # needs sage.graphs sage.groups
            (frozenset({(-1, -1), (1, 1)}), frozenset({(-1, 0), (1, 0)}))

            sage: cell24 = LatticePolytope_PPL(
            ....:     (1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,1), (1,-1,-1,1), (0,0,-1,1),
            ....:     (0,-1,0,1), (-1,0,0,1), (1,0,0,-1), (0,1,0,-1), (0,0,1,-1), (-1,1,1,-1),
            ....:     (1,-1,-1,0), (0,0,-1,0), (0,-1,0,0), (-1,0,0,0), (1,-1,0,0), (1,0,-1,0),
            ....:     (0,1,1,-1), (-1,1,1,0), (-1,1,0,0), (-1,0,1,0), (0,-1,-1,1), (0,0,0,-1))
            sage: fibers = [f.vertices() for f in cell24.fibration_generator(2)]
            sage: cell24.pointsets_mod_automorphism(fibers)     # long time             # needs sage.graphs sage.groups
            (frozenset({(-1, 0, 0, 0),
                        (-1, 0, 0, 1),
                        (0, 0, 0, -1),
                        (0, 0, 0, 1),
                        (1, 0, 0, -1),
                        (1, 0, 0, 0)}),
             frozenset({(-1, 0, 0, 0), (-1, 1, 1, 0), (1, -1, -1, 0), (1, 0, 0, 0)}))
        """
        points = set()
        for ps in pointsets:
            points.update(ps)
        points = tuple(sorted(points))
        Aut = self.lattice_automorphism_group(points,
                                              point_labels=tuple(range(len(points))))
        point_to_index = {p: i for i, p in enumerate(points)}
        indexsets = set(frozenset(point_to_index[p] for p in ps)
                        for ps in pointsets)
        orbits = []
        while indexsets:
            idx = indexsets.pop()
            min_orb = idx
            for g in Aut:
                g_idx = frozenset(g(i) for i in idx)
                if sorted(g_idx) < sorted(min_orb):
                    min_orb = g_idx
                indexsets.difference_update([g_idx])
            orbits.append(frozenset(points[i] for i in min_orb))
        return tuple(sorted(orbits, key=sorted))

    @cached_method
    def ambient_space(self):
        r"""
        Return the ambient space.

        OUTPUT:

        The free module `\ZZ^d`, where `d` is the ambient space
        dimension.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: point = LatticePolytope_PPL((1,2,3))
            sage: point.ambient_space()
            Ambient free module of rank 3 over the principal ideal domain Integer Ring
        """
        from sage.modules.free_module import FreeModule
        return FreeModule(ZZ, self.space_dimension())

    def contains(self, point_coordinates):
        r"""
        Test whether point is contained in the polytope.

        INPUT:

        - ``point_coordinates`` -- a list/tuple/iterable of rational
          numbers. The coordinates of the point.

        OUTPUT: Boolean.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: line = LatticePolytope_PPL((1,2,3), (-1,-2,-3))
            sage: line.contains([0,0,0])
            True
            sage: line.contains([1,0,0])
            False
        """
        p = C_Polyhedron(point(Linear_Expression(list(point_coordinates), 1)))
        is_included = Poly_Con_Relation.is_included()
        for c in self.constraints():
            if not p.relation_with(c).implies(is_included):
                return False
        return True

    @cached_method
    def contains_origin(self):
        """
        Test whether the polytope contains the origin

        OUTPUT: Boolean.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: LatticePolytope_PPL((1,2,3), (-1,-2,-3)).contains_origin()
            True
            sage: LatticePolytope_PPL((1,2,5), (-1,-2,-3)).contains_origin()
            False
        """
        return self.contains(self.ambient_space().zero())

    @cached_method
    def affine_space(self):
        r"""
        Return the affine space spanned by the polytope.

        OUTPUT:

        The free module `\ZZ^n`, where `n` is the dimension of the
        affine space spanned by the points of the polytope.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: point = LatticePolytope_PPL((1,2,3))
            sage: point.affine_space()
            Free module of degree 3 and rank 0 over Integer Ring
            Echelon basis matrix:
            []
            sage: line = LatticePolytope_PPL((1,1,1), (1,2,3))
            sage: line.affine_space()
            Free module of degree 3 and rank 1 over Integer Ring
            Echelon basis matrix:
            [0 1 2]
        """
        vertices = self.vertices()
        if not self.contains_origin():
            v0 = vertices[0]
            vertices = [v-v0 for v in vertices]
        return self.ambient_space().span(vertices).saturation()

    def affine_lattice_polytope(self):
        """
        Return the lattice polytope restricted to
        :meth:`affine_space`.

        OUTPUT:

        A new, full-dimensional lattice polytope.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: poly_4d = LatticePolytope_PPL((-9,-6,0,0), (0,1,0,0), (1,0,0,0));  poly_4d
            A 2-dimensional lattice polytope in ZZ^4 with 3 vertices
            sage: poly_4d.space_dimension()
            4
            sage: poly_2d = poly_4d.affine_lattice_polytope();  poly_2d
            A 2-dimensional lattice polytope in ZZ^2 with 3 vertices
            sage: poly_2d.space_dimension()
            2
        """
        V = self.affine_space()
        if self.contains_origin():
            vertices = [ V.coordinates(v) for v in self.vertices() ]
        else:
            v0 = vertices[0]
            vertices = [ V.coordinates(v-v0) for v in self.vertices() ]
        return LatticePolytope_PPL(*vertices)

    def base_projection(self, fiber):
        """
        The projection that maps the sub-polytope ``fiber`` to a
        single point.

        OUTPUT:

        The quotient module of the ambient space modulo the
        :meth:`affine_space` spanned by the fiber.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: poly = LatticePolytope_PPL((-9,-6,-1,-1),
            ....:                            (0,0,0,1), (0,0,1,0), (0,1,0,0), (1,0,0,0))
            sage: fiber = next(poly.fibration_generator(2))
            sage: poly.base_projection(fiber)
            Finitely generated module V/W over Integer Ring with invariants (0, 0)
        """
        return self.ambient_space().quotient(fiber.affine_space())

    def base_projection_matrix(self, fiber):
        """
        The projection that maps the sub-polytope ``fiber`` to a
        single point.

        OUTPUT:

        An integer matrix that represents the projection to the
        base.

        .. SEEALSO::

            The :meth:`base_projection` yields equivalent information,
            and is easier to use. However, just returning the matrix
            has lower overhead.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: poly = LatticePolytope_PPL((-9,-6,-1,-1),
            ....:                            (0,0,0,1), (0,0,1,0), (0,1,0,0), (1,0,0,0))
            sage: fiber = next(poly.fibration_generator(2))
            sage: poly.base_projection_matrix(fiber)
            [ 0  0 -1  0]
            [ 0  0  0 -1]

        Note that the basis choice in :meth:`base_projection` for the
        quotient is usually different::

            sage: proj = poly.base_projection(fiber)
            sage: proj_matrix = poly.base_projection_matrix(fiber)
            sage: [proj(p) for p in poly.integral_points()]
            [(-1, -1), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 1), (1, 0)]
            sage: [proj_matrix*p for p in poly.integral_points()]
            [(1, 1), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, -1), (-1, 0)]
        """
        return matrix(ZZ, fiber.vertices()).right_kernel_matrix()

    def base_rays(self, fiber, points):
        r"""
        Return the primitive lattice vectors that generate the
        direction given by the base projection of points.

        INPUT:

        - ``fiber`` -- a sub-lattice polytope defining the
          :meth:`base_projection`.

        - ``points`` -- the points to project to the base.

        OUTPUT:

        A tuple of primitive `\ZZ`-vectors.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: poly = LatticePolytope_PPL((-9,-6,-1,-1),
            ....:                            (0,0,0,1), (0,0,1,0), (0,1,0,0), (1,0,0,0))
            sage: fiber = next(poly.fibration_generator(2))
            sage: poly.base_rays(fiber, poly.integral_points_not_interior_to_facets())
            ((-1, -1), (0, 1), (1, 0))

            sage: p = LatticePolytope_PPL((1,0), (1,2), (-1,0))
            sage: f = LatticePolytope_PPL((1,0), (-1,0))
            sage: p.base_rays(f, p.integral_points())
            ((1),)
        """
        quo = self.base_projection(fiber)
        vertices = set()
        for p in points:
            v = quo(p).vector()
            if v.is_zero():
                continue
            d = GCD_list(v.list())
            if d > 1:
                v = v.__copy__()
                for i in range(v.degree()):
                    v[i] /= d
                v.set_immutable()
            vertices.add(v)
        return tuple(sorted(vertices))

    @cached_method
    def has_IP_property(self):
        """
        Whether the lattice polytope has the IP property.

        That is, the polytope is full-dimensional and the origin is a
        interior point not on the boundary.

        OUTPUT: Boolean.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: LatticePolytope_PPL((-1,-1), (0,1), (1,0)).has_IP_property()
            True
            sage: LatticePolytope_PPL((-1,-1), (1,1)).has_IP_property()
            False
        """
        origin = C_Polyhedron(point(0*Variable(self.space_dimension())))
        is_included = Poly_Con_Relation.is_included()
        saturates = Poly_Con_Relation.saturates()
        for c in self.constraints():
            rel = origin.relation_with(c)
            if (not rel.implies(is_included)) or rel.implies(saturates):
                return False
        return True

    @cached_method
    def restricted_automorphism_group(self, vertex_labels=None):
        r"""
        Return the restricted automorphism group.

        First, let the linear automorphism group be the subgroup of
        the Euclidean group `E(d) = GL(d,\RR) \ltimes \RR^d`
        preserving the `d`-dimensional polyhedron. The Euclidean group
        acts in the usual way `\vec{x}\mapsto A\vec{x}+b` on the
        ambient space. The restricted automorphism group is the
        subgroup of the linear automorphism group generated by
        permutations of vertices. If the polytope is full-dimensional,
        it is equal to the full (unrestricted) automorphism group.

        INPUT:

        - ``vertex_labels`` -- a tuple or ``None`` (default). The
          labels of the vertices that will be used in the output
          permutation group. By default, the vertices are used
          themselves.

        OUTPUT:

        A
        :class:`PermutationGroup<sage.groups.perm_gps.permgroup.PermutationGroup_generic>`
        acting on the vertices (or the ``vertex_labels``, if
        specified).

        REFERENCES:

        [BSS2009]_

        EXAMPLES::

            sage: # needs sage.graphs sage.groups
            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: Z3square = LatticePolytope_PPL((0,0), (1,2), (2,1), (3,3))
            sage: G1234 = Z3square.restricted_automorphism_group(
            ....:     vertex_labels=(1,2,3,4))
            sage: G1234 == PermutationGroup([[(2,3)], [(1,2),(3,4)]])
            True
            sage: G = Z3square.restricted_automorphism_group()
            sage: G == PermutationGroup([[((1,2),(2,1))], [((0,0),(1,2)),
            ....:                         ((2,1),(3,3))], [((0,0),(3,3))]])
            True
            sage: set(G.domain()) == set(Z3square.vertices())
            True
            sage: (set(tuple(x) for x in G.orbit(Z3square.vertices()[0]))
            ....:   == set([(0, 0), (1, 2), (3, 3), (2, 1)]))
            True
            sage: cell24 = LatticePolytope_PPL(
            ....:     (1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,1), (1,-1,-1,1), (0,0,-1,1),
            ....:     (0,-1,0,1), (-1,0,0,1), (1,0,0,-1), (0,1,0,-1), (0,0,1,-1), (-1,1,1,-1),
            ....:     (1,-1,-1,0), (0,0,-1,0), (0,-1,0,0), (-1,0,0,0), (1,-1,0,0), (1,0,-1,0),
            ....:     (0,1,1,-1), (-1,1,1,0), (-1,1,0,0), (-1,0,1,0), (0,-1,-1,1), (0,0,0,-1))
            sage: cell24.restricted_automorphism_group().cardinality()
            1152
        """
        if not self.is_full_dimensional():
            return self.affine_lattice_polytope().\
                restricted_automorphism_group(vertex_labels=vertex_labels)
        if vertex_labels is None:
            vertex_labels = self.vertices()
        from sage.graphs.graph import Graph
        # good coordinates for the vertices
        v_list = []
        for v in self.minimized_generators():
            assert v.divisor() == 1
            v_coords = (1,) + tuple(Integer(mpz) for mpz in v.coefficients())
            v_list.append(vector(v_coords))

        # Finally, construct the graph
        Qinv = sum( v.column() * v.row() for v in v_list ).inverse()
        G = Graph()
        for i in range(len(v_list)):
            for j in range(i+1,len(v_list)):
                v_i = v_list[i]
                v_j = v_list[j]
                G.add_edge(vertex_labels[i], vertex_labels[j], v_i * Qinv * v_j)
        return G.automorphism_group(edge_labels=True)

    @cached_method
    def lattice_automorphism_group(self, points=None, point_labels=None):
        """
        The integral subgroup of the restricted automorphism group.

        INPUT:

        - ``points`` -- A tuple of coordinate vectors or ``None``
          (default). If specified, the points must form complete
          orbits under the lattice automorphism group. If ``None`` all
          vertices are used.

        - ``point_labels`` -- A tuple of labels for the ``points`` or
          ``None`` (default). These will be used as labels for the do
          permutation group. If ``None``, the ``points`` will be used
          themselves.

        OUTPUT:

        The integral subgroup of the restricted automorphism group
        acting on the given ``points``, or all vertices if not
        specified.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: Z3square = LatticePolytope_PPL((0,0), (1,2), (2,1), (3,3))
            sage: Z3square.lattice_automorphism_group()                                 # needs sage.graphs sage.groups
            Permutation Group with generators [(), ((1,2),(2,1)),
            ((0,0),(3,3)), ((0,0),(3,3))((1,2),(2,1))]

            sage: G1 = Z3square.lattice_automorphism_group(point_labels=(1,2,3,4))      # needs sage.graphs sage.groups
            sage: G1                                                                    # needs sage.graphs sage.groups
            Permutation Group with generators [(), (2,3), (1,4), (1,4)(2,3)]
            sage: G1.cardinality()                                                      # needs sage.graphs sage.groups
            4

            sage: G2 = Z3square.restricted_automorphism_group(vertex_labels=(1,2,3,4))  # needs sage.graphs sage.groups
            sage: G2 == PermutationGroup([[(2,3)], [(1,2),(3,4)], [(1,4)]])             # needs sage.graphs sage.groups
            True
            sage: G2.cardinality()                                                      # needs sage.graphs sage.groups
            8

            sage: points = Z3square.integral_points(); points
            ((0, 0), (1, 1), (1, 2), (2, 1), (2, 2), (3, 3))
            sage: Z3square.lattice_automorphism_group(points,                           # needs sage.graphs sage.groups
            ....:                                     point_labels=(1,2,3,4,5,6))
            Permutation Group with generators [(), (3,4), (1,6)(2,5), (1,6)(2,5)(3,4)]

        Point labels also work for lattice polytopes that are not
        full-dimensional, see :trac:`16669`::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: lp = LatticePolytope_PPL((1,0,0), (0,1,0), (-1,-1,0))
            sage: lp.lattice_automorphism_group(point_labels=(0,1,2))                   # needs sage.graphs sage.groups
            Permutation Group with generators [(), (1,2), (0,1), (0,1,2), (0,2,1), (0,2)]
        """
        if not self.is_full_dimensional():
            return self.affine_lattice_polytope().lattice_automorphism_group(
                point_labels=point_labels)

        if points is None:
            points = self.vertices()
        if point_labels is None:
            point_labels = tuple(points)
        points = [ vector(ZZ, [1]+v.list()) for v in points ]
        for p in points:
            p.set_immutable()

        vertices = [ vector(ZZ, [1]+v.list()) for v in self.vertices() ]
        pivots = matrix(ZZ, vertices).pivot_rows()
        basis = matrix(ZZ, [ vertices[i] for i in pivots ])
        Mat_ZZ = basis.parent()
        basis_inverse = basis.inverse()

        from sage.groups.perm_gps.permgroup import PermutationGroup
        lattice_gens = []
        G = self.restricted_automorphism_group(
            vertex_labels=tuple(range(len(vertices))))
        for g in G:
            image = matrix(ZZ, [ vertices[g(i)] for i in pivots ])
            m = basis_inverse*image
            if m not in Mat_ZZ:
                continue
            perm_list = [ point_labels[points.index(p*m)]
                          for p in points ]
            lattice_gens.append(perm_list)
        return PermutationGroup(lattice_gens, domain=point_labels)

    def sub_polytope_generator(self):
        """
        Generate the maximal lattice sub-polytopes.

        OUTPUT:

        A generator yielding the maximal (with respect to inclusion)
        lattice sub polytopes. That is, each can be gotten as the
        convex hull of the integral points of ``self`` with one vertex
        removed.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: P = LatticePolytope_PPL((1,0,0), (0,1,0), (0,0,1), (-1,-1,-1))
            sage: for p in P.sub_polytope_generator():
            ....:     print(p.vertices())
            ((0, 0, 0), (0, 0, 1), (0, 1, 0), (1, 0, 0))
            ((-1, -1, -1), (0, 0, 0), (0, 1, 0), (1, 0, 0))
            ((-1, -1, -1), (0, 0, 0), (0, 0, 1), (1, 0, 0))
            ((-1, -1, -1), (0, 0, 0), (0, 0, 1), (0, 1, 0))
        """
        pointset = set(self.integral_points())
        for v in self.vertices():
            sub = list(pointset.difference([v]))
            yield LatticePolytope_PPL(*sub)

    @cached_method
    def _find_isomorphism_to_subreflexive_polytope(self):
        """
        Find an isomorphism to a sub-polytope of a maximal reflexive
        polytope.

        OUTPUT:

        A tuple consisting of the ambient reflexive polytope, the
        subpolytope, and the embedding of ``self`` into the ambient
        polytope.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: polygon = LatticePolytope_PPL((0,0,2,1), (0,1,2,0), (2,3,0,0), (2,0,0,3))
            sage: polygon._find_isomorphism_to_subreflexive_polytope()
            (A 2-dimensional lattice polytope in ZZ^2 with 3 vertices,
             A 2-dimensional lattice polytope in ZZ^2 with 4 vertices,
             The map A*x+b with
               A=
                 [ 1  1]
                 [ 0  1]
                 [-1 -1]
                 [ 1  0]
               b = (-1, 0, 3, 0))
            sage: ambient, sub, embedding = _
            sage: ambient.vertices()
            ((0, 0), (0, 3), (3, 0))
            sage: sub.vertices()
            ((0, 1), (3, 0), (0, 3), (1, 0))
        """
        from .ppl_lattice_polygon import sub_reflexive_polygons
        from sage.geometry.polyhedron.lattice_euclidean_group_element import \
            LatticePolytopesNotIsomorphicError, LatticePolytopeNoEmbeddingError
        for p, ambient in sub_reflexive_polygons():
            try:
                return (ambient, p, p.find_isomorphism(self))
            except LatticePolytopesNotIsomorphicError:
                pass
        raise LatticePolytopeNoEmbeddingError('not a sub-polytope of a reflexive polygon')

    def embed_in_reflexive_polytope(self, output='hom'):
        """
        Find an embedding as a sub-polytope of a maximal reflexive
        polytope.

        INPUT:

        - ``hom`` -- string. One of ``'hom'`` (default),
          ``'polytope'``, or ``points``. How the embedding is
          returned. See the output section for details.

        OUTPUT:

        An embedding into a reflexive polytope. Depending on the
        ``output`` option slightly different data is returned.

        - If ``output='hom'``, a map from a reflexive polytope onto
          ``self`` is returned.

        - If ``output='polytope'``, a reflexive polytope that contains
          ``self`` (up to a lattice linear transformation) is
          returned. That is, the domain of the ``output='hom'`` map is
          returned. If the affine span of ``self`` is less or equal
          2-dimensional, the output is one of the following three
          possibilities:

          :func:`~sage.geometry.polyhedron.ppl_lattice_polygon.polar_P2_polytope`,
          :func:`~sage.geometry.polyhedron.ppl_lattice_polygon.polar_P1xP1_polytope`,
          or
          :func:`~sage.geometry.polyhedron.ppl_lattice_polygon.polar_P2_112_polytope`.

        - If ``output='points'``, a dictionary containing the integral
          points of ``self`` as keys and the corresponding integral
          point of the reflexive polytope as value.

        If there is no such embedding, a
        :class:`~sage.geometry.polyhedron.lattice_euclidean_group_element.LatticePolytopeNoEmbeddingError`
        is raised. Even if it exists, the ambient reflexive polytope
        is usually not uniquely determined and a random but fixed
        choice will be returned.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
            sage: polygon = LatticePolytope_PPL((0,0,2,1), (0,1,2,0), (2,3,0,0), (2,0,0,3))
            sage: polygon.embed_in_reflexive_polytope()
            The map A*x+b with
              A=
                [ 1  1]
                [ 0  1]
                [-1 -1]
                [ 1  0]
              b = (-1, 0, 3, 0)
            sage: polygon.embed_in_reflexive_polytope('polytope')
            A 2-dimensional lattice polytope in ZZ^2 with 3 vertices
            sage: polygon.embed_in_reflexive_polytope('points')
            {(0, 0, 2, 1): (1, 0),
             (0, 1, 2, 0): (0, 1),
             (1, 0, 1, 2): (2, 0),
             (1, 1, 1, 1): (1, 1),
             (1, 2, 1, 0): (0, 2),
             (2, 0, 0, 3): (3, 0),
             (2, 1, 0, 2): (2, 1),
             (2, 2, 0, 1): (1, 2),
             (2, 3, 0, 0): (0, 3)}

            sage: LatticePolytope_PPL((0,0), (4,0), (0,4)).embed_in_reflexive_polytope()
            Traceback (most recent call last):
            ...
            LatticePolytopeNoEmbeddingError: not a sub-polytope of a reflexive polygon
        """
        if self.affine_dimension() > 2:
            raise NotImplementedError('can only embed in reflexive polygons')
        ambient, subreflexive, hom = self._find_isomorphism_to_subreflexive_polytope()
        if output == 'hom':
            return hom
        elif output == 'polytope':
            return ambient
        elif output == 'points':
            points = dict()
            for p in subreflexive.integral_points():
                points[ tuple(hom(p)) ] = p
            return points
        else:
            raise ValueError('output='+str(output)+' is not valid.')
