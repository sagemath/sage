r"""
Base class for polyhedra: Implementation of the :class:`ConvexSet_base` API

Define methods that exist for convex sets,
but not constructions such as dilation or product.
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

from sage.structure.element import coerce_binop
from sage.structure.richcmp import rich_to_bool, op_NE
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.modules.free_module_element import vector
from .base0 import Polyhedron_base0
from sage.geometry.convex_set import ConvexSet_closed
from sage.geometry.relative_interior import RelativeInterior


class Polyhedron_base1(Polyhedron_base0, ConvexSet_closed):
    """
    Convex set methods for polyhedra,
    but not constructions such as dilation or product.

    See :class:`sage.geometry.polyhedron.base.Polyhedron_base`.

    TESTS::

        sage: from sage.geometry.polyhedron.base1 import Polyhedron_base1
        sage: P = polytopes.cube()
        sage: Q = polytopes.cube()
        sage: Polyhedron_base1.__hash__(P) == Polyhedron_base1.__hash__(Q)
        True
        sage: Polyhedron_base1.__repr__(P)
        'A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 8 vertices'
        sage: Polyhedron_base1.is_empty(P)
        False
        sage: Polyhedron_base1.is_universe(P)
        False
        sage: Polyhedron_base1.dim(P)
        3
        sage: Polyhedron_base1.ambient_vector_space(P)
        Vector space of dimension 3 over Rational Field
        sage: Polyhedron_base1.ambient_dim(P)
        3
        sage: Polyhedron_base1.an_affine_basis(P)
        [A vertex at (-1, -1, -1),
        A vertex at (1, -1, -1),
        A vertex at (1, -1, 1),
        A vertex at (1, 1, -1)]
        sage: list(Polyhedron_base1._some_elements_(P))
        [(0, 0, 0),
        (1, -1, -1),
        (1, 0, -1),
        (1, 1/2, 0),
        (1, -1/4, 1/2),
        (0, -5/8, 3/4)]
        sage: Polyhedron_base1.contains(P, vector([1, 1, 1]))
        True
        sage: Polyhedron_base1.interior_contains(P, vector([1, 1, 1]))
        False
        sage: Polyhedron_base1.is_relatively_open(P)
        False
        sage: Polyhedron_base1.relative_interior.f(P) == Polyhedron_base1.interior.f(P)
        True
    """

    def __hash__(self):
        r"""
        TESTS::

            sage: # needs sage.rings.number_field
            sage: K.<a> = QuadraticField(2)
            sage: p = Polyhedron(vertices=[(0, 1, a), (3, a, 5)],
            ....:                rays=[(a, 2, 3), (0, 0, 1)],
            ....:                base_ring=K)
            sage: q = Polyhedron(vertices=[(3, a, 5), (0, 1, a)],
            ....:                rays=[(0, 0, 1), (a, 2, 3)],
            ....:                base_ring=K)
            sage: hash(p) == hash(q)
            True
        """
        # TODO: find something better *but* fast
        return hash((self.dim(),
                     self.ambient_dim(),
                     self.n_Hrepresentation(),
                     self.n_Vrepresentation(),
                     self.n_equations(),
                     self.n_facets(),
                     self.n_inequalities(),
                     self.n_lines(),
                     self.n_rays(),
                     self.n_vertices()))

    def _repr_(self):
        """
        Return a description of the polyhedron.

        EXAMPLES::

            sage: poly_test = Polyhedron(vertices = [[1,2,3,4],[2,1,3,4],[4,3,2,1]])
            sage: poly_test._repr_()
            'A 2-dimensional polyhedron in ZZ^4 defined as the convex hull of 3 vertices'
            sage: grammar_test = Polyhedron(vertices = [[1,1,1,1,1,1]])
            sage: grammar_test._repr_()
            'A 0-dimensional polyhedron in ZZ^6 defined as the convex hull of 1 vertex'
        """
        desc = ''
        if self.n_vertices() == 0:
            desc += 'The empty polyhedron'
        else:
            desc += 'A ' + repr(self.dim()) + '-dimensional polyhedron'
        desc += ' in '
        desc += self.parent()._repr_ambient_module()

        if self.n_vertices() > 0:
            desc += ' defined as the convex hull of '
            desc += repr(self.n_vertices())
            if self.n_vertices() == 1:
                desc += ' vertex'
            else:
                desc += ' vertices'

            if self.n_rays() > 0:
                if self.n_lines() > 0:
                    desc += ", "
                else:
                    desc += " and "
                desc += repr(self.n_rays())
                if self.n_rays() == 1:
                    desc += ' ray'
                else:
                    desc += ' rays'

            if self.n_lines() > 0:
                if self.n_rays() > 0:
                    desc += ", "
                else:
                    desc += " and "
                desc += repr(self.n_lines())
                if self.n_lines() == 1:
                    desc += ' line'
                else:
                    desc += ' lines'

        return desc

    def _richcmp_(self, other, op):
        """
        Compare ``self`` and ``other``.

        INPUT:

        - ``other`` -- a polyhedron

        OUTPUT:

        If ``other`` is a polyhedron, then the comparison
        operator "less or equal than" means "is contained in", and
        "less than" means "is strictly contained in".

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(1,0), (0,1)], rays=[(1,1)])
            sage: Q = Polyhedron(vertices=[(1,0), (0,1)])
            sage: P >= Q
            True
            sage: Q <= P
            True
            sage: P == P
            True

        The polytope ``Q`` is strictly contained in ``P``::

            sage: P > Q
            True
            sage: P < Q
            False
            sage: P == Q
            False

        Test that we have fixed a problem revealed in :issue:`31701`,
        where neither of the two polyhedra contains the other::

            sage: P = Polyhedron(vertices=[(1, 1), (0, 0), (1, 2)])
            sage: Q = Polyhedron(vertices=[(1, 2), (0, 0), (0, 2)])
            sage: Q < P
            False
            sage: P > Q
            False
        """
        if self.Vrepresentation() is None or other.Vrepresentation() is None:
            raise RuntimeError('some V representation is missing')
            # make sure deleted polyhedra are not used in cache

        if self.ambient_dim() != other.ambient_dim():
            return op == op_NE

        c0 = self._is_subpolyhedron(other)
        c1 = other._is_subpolyhedron(self)
        if c0 and c1:
            return rich_to_bool(op, 0)
        elif c0:
            return rich_to_bool(op, -1)
        elif c1:
            return rich_to_bool(op, 1)
        else:
            return op == op_NE

    @coerce_binop
    def _is_subpolyhedron(self, other):
        """
        Test whether ``self`` is a (not necessarily strict)
        sub-polyhedron of ``other``.

        INPUT:

        - ``other`` -- a :class:`Polyhedron`

        OUTPUT: boolean

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(1,0), (0,1)], rays=[(1,1)])
            sage: Q = Polyhedron(vertices=[(1,0), (0,1)])
            sage: P._is_subpolyhedron(Q)
            False
            sage: Q._is_subpolyhedron(P)
            True
        """
        return all(other_H.contains(self_V)
                   for other_H in other.Hrepresentation()
                   for self_V in self.Vrepresentation())

    def is_empty(self):
        """
        Test whether the polyhedron is the empty polyhedron.

        OUTPUT: boolean

        EXAMPLES::

            sage: P = Polyhedron(vertices=[[1,0,0],[0,1,0],[0,0,1]]);  P
            A 2-dimensional polyhedron in ZZ^3 defined as the convex hull of 3 vertices
            sage: P.is_empty(), P.is_universe()
            (False, False)

            sage: Q = Polyhedron(vertices=());  Q
            The empty polyhedron in ZZ^0
            sage: Q.is_empty(), Q.is_universe()
            (True, False)

            sage: R = Polyhedron(lines=[(1,0),(0,1)]);  R
            A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 1 vertex and 2 lines
            sage: R.is_empty(), R.is_universe()
            (False, True)
        """
        return self.n_Vrepresentation() == 0

    def is_universe(self):
        """
        Test whether the polyhedron is the whole ambient space.

        OUTPUT: boolean

        EXAMPLES::

            sage: P = Polyhedron(vertices=[[1,0,0],[0,1,0],[0,0,1]]);  P
            A 2-dimensional polyhedron in ZZ^3 defined as the convex hull of 3 vertices
            sage: P.is_empty(), P.is_universe()
            (False, False)

            sage: Q = Polyhedron(vertices=());  Q
            The empty polyhedron in ZZ^0
            sage: Q.is_empty(), Q.is_universe()
            (True, False)

            sage: R = Polyhedron(lines=[(1,0),(0,1)]);  R
            A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 1 vertex and 2 lines
            sage: R.is_empty(), R.is_universe()
            (False, True)
        """
        return self.n_Hrepresentation() == 0

    def dim(self):
        """
        Return the dimension of the polyhedron.

        OUTPUT: -1 if the polyhedron is empty, otherwise a nonnegative integer

        EXAMPLES::

            sage: simplex = Polyhedron(vertices = [[1,0,0,0],[0,0,0,1],[0,1,0,0],[0,0,1,0]])
            sage: simplex.dim()
            3
            sage: simplex.ambient_dim()
            4

        The empty set is a special case (:issue:`12193`)::

            sage: P1=Polyhedron(vertices=[[1,0,0],[0,1,0],[0,0,1]])
            sage: P2=Polyhedron(vertices=[[2,0,0],[0,2,0],[0,0,2]])
            sage: P12 = P1.intersection(P2)
            sage: P12
            The empty polyhedron in ZZ^3
            sage: P12.dim()
            -1
        """
        if self.n_Vrepresentation() == 0:
            return -1   # the empty set
        else:
            return self.ambient_dim() - self.n_equations()

    dimension = dim

    def Vrepresentation_space(self):
        r"""
        Return the ambient free module.

        OUTPUT: a free module over the base ring of dimension :meth:`ambient_dim`

        EXAMPLES::

            sage: poly_test = Polyhedron(vertices = [[1,0,0,0],[0,1,0,0]])
            sage: poly_test.Vrepresentation_space()
            Ambient free module of rank 4 over the principal ideal domain Integer Ring
            sage: poly_test.ambient_space() is poly_test.Vrepresentation_space()
            True
        """
        return self.parent().Vrepresentation_space()

    def Hrepresentation_space(self):
        r"""
        Return the linear space containing the H-representation vectors.

        OUTPUT: a free module over the base ring of dimension :meth:`ambient_dim` + 1

        EXAMPLES::

            sage: poly_test = Polyhedron(vertices = [[1,0,0,0],[0,1,0,0]])
            sage: poly_test.Hrepresentation_space()
            Ambient free module of rank 5 over the principal ideal domain Integer Ring
        """
        return self.parent().Hrepresentation_space()

    ambient_space = Vrepresentation_space

    def ambient_vector_space(self, base_field=None):
        r"""
        Return the ambient vector space.

        It is the ambient free module (:meth:`Vrepresentation_space`) tensored
        with a field.

        INPUT:

        - ``base_field`` -- a field (default: the fraction field of the base ring)

        EXAMPLES::

            sage: poly_test = Polyhedron(vertices = [[1,0,0,0],[0,1,0,0]])
            sage: poly_test.ambient_vector_space()
            Vector space of dimension 4 over Rational Field
            sage: poly_test.ambient_vector_space() is poly_test.ambient()
            True

            sage: poly_test.ambient_vector_space(AA)                                    # needs sage.rings.number_field
            Vector space of dimension 4 over Algebraic Real Field
            sage: poly_test.ambient_vector_space(RDF)
            Vector space of dimension 4 over Real Double Field
            sage: poly_test.ambient_vector_space(SR)                                    # needs sage.symbolic
            Vector space of dimension 4 over Symbolic Ring
        """
        return self.Vrepresentation_space().vector_space(base_field=base_field)

    ambient = ambient_vector_space

    def ambient_dim(self):
        r"""
        Return the dimension of the ambient space.

        EXAMPLES::

            sage: poly_test = Polyhedron(vertices = [[1,0,0,0],[0,1,0,0]])
            sage: poly_test.ambient_dim()
            4
        """
        return self.parent().ambient_dim()

    def an_affine_basis(self):
        """
        Return points in ``self`` that form a basis for the affine span of ``self``.

        This implementation of the method :meth:`~sage.geometry.convex_set.ConvexSet_base.an_affine_basis`
        for polytopes guarantees the following:

        - All points are vertices.

        - The basis is obtained by considering a maximal chain of faces
          in the face lattice and picking for each cover relation
          one vertex that is in the difference. Thus this method
          is independent of the concrete realization of the polytope.

        For unbounded polyhedra, the result may contain arbitrary points of ``self``,
        not just vertices.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: P.an_affine_basis()
            [A vertex at (-1, -1, -1),
             A vertex at (1, -1, -1),
             A vertex at (1, -1, 1),
             A vertex at (1, 1, -1)]

            sage: P = polytopes.permutahedron(5)
            sage: P.an_affine_basis()
            [A vertex at (1, 2, 3, 5, 4),
             A vertex at (2, 1, 3, 5, 4),
             A vertex at (1, 3, 2, 5, 4),
             A vertex at (4, 1, 3, 5, 2),
             A vertex at (4, 2, 5, 3, 1)]

        Unbounded polyhedra::

            sage: p = Polyhedron(vertices=[(0, 0)], rays=[(1,0), (0,1)])
            sage: p.an_affine_basis()
            [A vertex at (0, 0), (1, 0), (0, 1)]
            sage: p = Polyhedron(vertices=[(2, 1)], rays=[(1,0), (0,1)])
            sage: p.an_affine_basis()
            [A vertex at (2, 1), (3, 1), (2, 2)]
            sage: p = Polyhedron(vertices=[(2, 1)], rays=[(1,0)], lines=[(0,1)])
            sage: p.an_affine_basis()
            [(2, 1), A vertex at (2, 0), (3, 0)]
        """
        chain = self.a_maximal_chain()[1:]  # we exclude the empty face
        chain_indices = [face.ambient_V_indices() for face in chain]
        basis_indices = []

        # We use in the following that elements in ``chain_indices`` are sorted lists
        # of V-indices.
        # Thus for each two faces we can easily find the first Vrep that differs.
        for dim, face in enumerate(chain_indices):
            if dim == 0:
                # Append the V-indices of the minimal face.
                basis_indices.extend(face[:])
                continue

            prev_face = chain_indices[dim-1]
            for i in range(len(prev_face)):
                if prev_face[i] != face[i]:
                    # We found a Vrep that ``face`` has, but its facet does not.
                    basis_indices.append(face[i])
                    break
            else:  # no break
                # ``prev_face`` contains all the same vertices as ``face`` until now.
                # But ``face`` is guaranteed to contain one more Vrep (at least).
                basis_indices.append(face[len(prev_face)])

        Vreps = [self.Vrepresentation()[i] for i in basis_indices]
        if all(vrep.is_vertex() for vrep in Vreps):
            return Vreps
        for vrep in Vreps:
            if vrep.is_vertex():
                vertex = vrep
        return [vrep if vrep.is_vertex() else vertex.vector() + vrep.vector()
                for vrep in Vreps]

    @abstract_method
    def a_maximal_chain(self):
        r"""
        Return a maximal chain of the face lattice in increasing order.

        Subclasses must provide an implementation of this method.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.base1 import Polyhedron_base1
            sage: P = polytopes.cube()
            sage: Polyhedron_base1.a_maximal_chain
            <abstract method a_maximal_chain at ...>
        """

    @cached_method
    def representative_point(self):
        """
        Return a "generic" point.

        .. SEEALSO::

            :meth:`sage.geometry.polyhedron.base.Polyhedron_base.center`.

        OUTPUT:

        A point as a coordinate vector. The point is chosen to be
        interior if possible. If the polyhedron is not
        full-dimensional, the point is in the relative interior. If
        the polyhedron is zero-dimensional, its single point is
        returned.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[(3,2)], rays=[(1,-1)])
            sage: p.representative_point()
            (4, 1)
            sage: p.center()
            (3, 2)

            sage: Polyhedron(vertices=[(3,2)]).representative_point()
            (3, 2)
        """
        accumulator = vector(self.base_ring(), [0]*self.ambient_dim())
        for v in self.vertex_generator():
            accumulator += v.vector()
        accumulator /= self.n_vertices()
        for r in self.ray_generator():
            accumulator += r.vector()
        accumulator.set_immutable()
        return accumulator

    def _some_elements_(self):
        r"""
        Generate some points of ``self``.

        If ``self`` is empty, no points are generated; no exception will be raised.

        EXAMPLES::

            sage: P = polytopes.simplex()
            sage: P.an_element()            # indirect doctest
            (1/4, 1/4, 1/4, 1/4)
            sage: P.some_elements()         # indirect doctest
            [(1/4, 1/4, 1/4, 1/4),
             (0, 0, 0, 1),
             (0, 0, 1/2, 1/2),
             (0, 1/2, 1/4, 1/4),
             (1/2, 1/4, 1/8, 1/8)]
        """
        if self.is_empty():
            return
        yield self.representative_point()
        vertex_iter = iter(self.vertex_generator())
        try:
            p = next(vertex_iter).vector()
            yield vector(p, immutable=True)
            for i in range(4):
                p = (p + next(vertex_iter).vector()) / 2
                yield vector(p, immutable=True)
        except StopIteration:
            pass

    def contains(self, point):
        """
        Test whether the polyhedron contains the given ``point``.

        .. SEEALSO::

            :meth:`interior_contains`, :meth:`relative_interior_contains`.

        INPUT:

        - ``point`` -- coordinates of a point (an iterable)

        OUTPUT: boolean

        EXAMPLES::

            sage: P = Polyhedron(vertices=[[1,1],[1,-1],[0,0]])
            sage: P.contains( [1,0] )
            True
            sage: P.contains( P.center() )  # true for any convex set
            True

        As a shorthand, one may use the usual ``in`` operator::

            sage: P.center() in P
            True
            sage: [-1,-1] in P
            False

        The point need not have coordinates in the same field as the
        polyhedron::

            sage: # needs sage.symbolic
            sage: ray = Polyhedron(vertices=[(0,0)], rays=[(1,0)], base_ring=QQ)
            sage: ray.contains([sqrt(2)/3,0])        # irrational coordinates are ok
            True
            sage: a = var('a')
            sage: ray.contains([a,0])                # a might be negative!
            False
            sage: assume(a>0)
            sage: ray.contains([a,0])
            True
            sage: ray.contains(['hello', 'kitty'])   # no common ring for coordinates
            False

        The empty polyhedron needs extra care, see :issue:`10238`::

            sage: empty = Polyhedron(); empty
            The empty polyhedron in ZZ^0
            sage: empty.contains([])
            False
            sage: empty.contains([0])               # not a point in QQ^0
            False
            sage: full = Polyhedron(vertices=[()]); full
            A 0-dimensional polyhedron in ZZ^0 defined as the convex hull of 1 vertex
            sage: full.contains([])
            True
            sage: full.contains([0])
            False

        TESTS:

        Passing non-iterable objects does not cause an exception, see :issue:`32013`::

            sage: None in Polyhedron(vertices=[(0,0)], rays=[(1,0)], base_ring=QQ)
            False
        """
        try:
            p = vector(point)
        except TypeError:  # point not iterable or no common ring for elements
            try:
                l = len(point)
            except TypeError:
                return False
            if l > 0:
                return False
            else:
                p = vector(self.base_ring(), [])

        if len(p) != self.ambient_dim():
            return False

        for H in self.Hrep_generator():
            if not H.contains(p):
                return False
        return True

    __contains__ = contains

    @cached_method
    def interior(self):
        """
        The interior of ``self``.

        OUTPUT:

        - either an empty polyhedron or an instance of
          :class:`~sage.geometry.relative_interior.RelativeInterior`

        EXAMPLES:

        If the polyhedron is full-dimensional, the result is the
        same as that of :meth:`relative_interior`::

            sage: P_full = Polyhedron(vertices=[[0,0],[1,1],[1,-1]])
            sage: P_full.interior()
            Relative interior of
             a 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 3 vertices

        If the polyhedron is of strictly smaller dimension than the
        ambient space, its interior is empty::

            sage: P_lower = Polyhedron(vertices=[[0,1], [0,-1]])
            sage: P_lower.interior()
            The empty polyhedron in ZZ^2

        TESTS::

            sage: Empty = Polyhedron(ambient_dim=2); Empty
            The empty polyhedron in ZZ^2
            sage: Empty.interior() is Empty
            True
        """
        if self.is_open():
            return self
        if not self.is_full_dimensional():
            return self.parent().element_class(self.parent(), None, None)
        return self.relative_interior()

    def interior_contains(self, point):
        """
        Test whether the interior of the polyhedron contains the
        given ``point``.

        .. SEEALSO::

            :meth:`contains`, :meth:`relative_interior_contains`.

        INPUT:

        - ``point`` -- coordinates of a point

        OUTPUT: boolean

        EXAMPLES::

            sage: P = Polyhedron(vertices=[[0,0],[1,1],[1,-1]])
            sage: P.contains( [1,0] )
            True
            sage: P.interior_contains( [1,0] )
            False

        If the polyhedron is of strictly smaller dimension than the
        ambient space, its interior is empty::

            sage: P = Polyhedron(vertices=[[0,1],[0,-1]])
            sage: P.contains( [0,0] )
            True
            sage: P.interior_contains( [0,0] )
            False

        The empty polyhedron needs extra care, see :issue:`10238`::

            sage: empty = Polyhedron(); empty
            The empty polyhedron in ZZ^0
            sage: empty.interior_contains([])
            False
        """
        try:
            p = vector(point)
        except TypeError:  # point not iterable or no common ring for elements
            try:
                l = len(point)
            except TypeError:
                return False
            if l > 0:
                return False
            else:
                p = vector(self.base_ring(), [])

        if len(p) != self.ambient_dim():
            return False

        for H in self.Hrep_generator():
            if not H.interior_contains(p):
                return False
        return True

    def is_relatively_open(self):
        r"""
        Return whether ``self`` is relatively open.

        OUTPUT: boolean

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(1,0), (-1,0)]); P
            A 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices
            sage: P.is_relatively_open()
            False

            sage: P0 = Polyhedron(vertices=[[1, 2]]); P0
            A 0-dimensional polyhedron in ZZ^2 defined as the convex hull of 1 vertex
            sage: P0.is_relatively_open()
            True

            sage: Empty = Polyhedron(ambient_dim=2); Empty
            The empty polyhedron in ZZ^2
            sage: Empty.is_relatively_open()
            True

            sage: Line = Polyhedron(vertices=[(1, 1)], lines=[(1, 0)]); Line
            A 1-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 1 line
            sage: Line.is_relatively_open()
            True
        """
        return not self.inequalities()

    @cached_method
    def relative_interior(self):
        """
        Return the relative interior of ``self``.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(1,0), (-1,0)])
            sage: ri_P = P.relative_interior(); ri_P
            Relative interior of
             a 1-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices
            sage: (0, 0) in ri_P
            True
            sage: (1, 0) in ri_P
            False

            sage: P0 = Polyhedron(vertices=[[1, 2]])
            sage: P0.relative_interior() is P0
            True

            sage: Empty = Polyhedron(ambient_dim=2)
            sage: Empty.relative_interior() is Empty
            True

            sage: Line = Polyhedron(vertices=[(1, 1)], lines=[(1, 0)])
            sage: Line.relative_interior() is Line
            True
        """
        if self.is_relatively_open():
            return self
        return RelativeInterior(self)

    def relative_interior_contains(self, point):
        """
        Test whether the relative interior of the polyhedron
        contains the given ``point``.

        .. SEEALSO::

            :meth:`contains`, :meth:`interior_contains`.

        INPUT:

        - ``point`` -- coordinates of a point

        OUTPUT: boolean

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(1,0), (-1,0)])
            sage: P.contains( (0,0) )
            True
            sage: P.interior_contains( (0,0) )
            False
            sage: P.relative_interior_contains( (0,0) )
            True
            sage: P.relative_interior_contains( (1,0) )
            False

        The empty polyhedron needs extra care, see :issue:`10238`::

            sage: empty = Polyhedron(); empty
            The empty polyhedron in ZZ^0
            sage: empty.relative_interior_contains([])
            False
        """
        try:
            p = vector(point)
        except TypeError:  # point not iterable or no common ring for elements
            try:
                l = len(point)
            except TypeError:
                return False
            if l > 0:
                return False
            else:
                p = vector(self.base_ring(), [])

        if len(p) != self.ambient_dim():
            return False

        for eq in self.equation_generator():
            if not eq.contains(p):
                return False

        for ine in self.inequality_generator():
            if not ine.interior_contains(p):
                return False

        return True
