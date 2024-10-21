# sage_setup: distribution = sagemath-polyhedra
"""
A class to keep information about faces of a polyhedron

This module gives you a tool to work with the faces of a polyhedron
and their relative position. First, you need to find the faces. To get
the faces in a particular dimension, use the
:meth:`~sage.geometry.polyhedron.base.face` method::

    sage: P = polytopes.cross_polytope(3)
    sage: P.faces(3)
    (A 3-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 6 vertices,)
    sage: [f.ambient_V_indices() for f in P.facets()]
    [(3, 4, 5),
     (2, 4, 5),
     (1, 3, 5),
     (1, 2, 5),
     (0, 3, 4),
     (0, 2, 4),
     (0, 1, 3),
     (0, 1, 2)]
    sage: [f.ambient_V_indices() for f in P.faces(1)]
    [(4, 5),
     (3, 5),
     (2, 5),
     (1, 5),
     (3, 4),
     (2, 4),
     (0, 4),
     (1, 3),
     (0, 3),
     (1, 2),
     (0, 2),
     (0, 1)]

or :meth:`~sage.geometry.polyhedron.base.face_lattice` to get the
whole face lattice as a poset::

    sage: P.face_lattice()                                                              # needs sage.combinat
    Finite lattice containing 28 elements

The faces are printed in shorthand notation where each integer is the
index of a vertex/ray/line in the same order as the containing
Polyhedron's :meth:`~sage.geometry.polyhedron.base.Vrepresentation` ::

    sage: face = P.faces(1)[8];  face
    A 1-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 2 vertices
    sage: face.ambient_V_indices()
    (0, 3)
    sage: P.Vrepresentation(0)
    A vertex at (-1, 0, 0)
    sage: P.Vrepresentation(3)
    A vertex at (0, 0, 1)
    sage: face.vertices()
    (A vertex at (-1, 0, 0), A vertex at (0, 0, 1))

The face itself is not represented by Sage's
:func:`sage.geometry.polyhedron.constructor.Polyhedron` class, but by
an auxiliary class to keep the information. You can get the face as a
polyhedron with the :meth:`PolyhedronFace.as_polyhedron` method::

    sage: face.as_polyhedron()
    A 1-dimensional polyhedron in ZZ^3 defined as the convex hull of 2 vertices
    sage: _.equations()
    (An equation (0, 1, 0) x + 0 == 0,
     An equation (1, 0, -1) x + 1 == 0)
"""

########################################################################
#       Copyright (C) 2012 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################

from sage.structure.richcmp import richcmp_method, richcmp
from sage.misc.cachefunc import cached_method
from sage.modules.free_module_element import vector
from sage.matrix.constructor import matrix
from sage.geometry.convex_set import ConvexSet_closed


#########################################################################
@richcmp_method
class PolyhedronFace(ConvexSet_closed):
    r"""
    A face of a polyhedron.

    This class is for use in
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.face_lattice`.

    INPUT:

    No checking is performed whether the H/V-representation indices
    actually determine a face of the polyhedron. You should not
    manually create :class:`PolyhedronFace` objects unless you know
    what you are doing.

    OUTPUT: a :class:`PolyhedronFace`

    EXAMPLES::

        sage: octahedron = polytopes.cross_polytope(3)
        sage: inequality = octahedron.Hrepresentation(2)
        sage: face_h = tuple([ inequality ])
        sage: face_v = tuple( inequality.incident() )
        sage: face_h_indices = [ h.index() for h in face_h ]
        sage: face_v_indices = [ v.index() for v in face_v ]
        sage: from sage.geometry.polyhedron.face import PolyhedronFace
        sage: face = PolyhedronFace(octahedron, face_v_indices, face_h_indices)
        sage: face
        A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 3 vertices
        sage: face.dim()
        2
        sage: face.ambient_V_indices()
        (0, 1, 2)
        sage: face.ambient_Hrepresentation()
        (An inequality (1, 1, 1) x + 1 >= 0,)
        sage: face.ambient_Vrepresentation()
        (A vertex at (-1, 0, 0), A vertex at (0, -1, 0), A vertex at (0, 0, -1))

    TESTS::

        sage: TestSuite(face).run()
    """

    def __init__(self, polyhedron, V_indices, H_indices):
        r"""
        The constructor.

        See :class:`PolyhedronFace` for more information.

        INPUT:

        - ``polyhedron`` -- a :class:`Polyhedron`; the ambient
          polyhedron

        - ``V_indices`` -- list of sorted integers; the indices of the
          face-spanning V-representation objects in the ambient
          polyhedron

        - ``H_indices`` -- list of sorted integers; the indices of the
          H-representation objects of the ambient polyhedron that are
          saturated on the face

        TESTS::

            sage: from sage.geometry.polyhedron.face import PolyhedronFace
            sage: PolyhedronFace(Polyhedron(), [], [])   # indirect doctest
            A -1-dimensional face of a Polyhedron in ZZ^0
            sage: TestSuite(_).run(skip='_test_pickling')
        """
        self._polyhedron = polyhedron
        self._ambient_Vrepresentation_indices = tuple(V_indices)
        self._ambient_Hrepresentation_indices = tuple(H_indices)
        self._ambient_Vrepresentation = tuple( polyhedron.Vrepresentation(i) for i in V_indices )
        self._ambient_Hrepresentation = tuple( polyhedron.Hrepresentation(i) for i in H_indices )
        if polyhedron.is_mutable():
            polyhedron._add_dependent_object(self)

    def __hash__(self):
        r"""
        TESTS::

            sage: P = Polyhedron([[0,0],[0,1],[23,3],[9,12]])
            sage: list(map(hash, P.faces(1)))  # random
            [2377119663630407734,
             2377136578164722109,
             5966674064902575359,
             4795242501625591634]
        """
        return hash((self._polyhedron, self._ambient_Vrepresentation_indices))

    def vertex_generator(self):
        """
        Return a generator for the vertices of the face.

        EXAMPLES::

            sage: triangle = Polyhedron(vertices=[[1,0],[0,1],[1,1]])
            sage: face = triangle.facets()[0]
            sage: for v in face.vertex_generator(): print(v)
            A vertex at (1, 0)
            A vertex at (1, 1)
            sage: type(face.vertex_generator())
            <... 'generator'>

        TESTS::

            sage: TestSuite(face).run(skip='_test_pickling')
        """
        for V in self.ambient_Vrepresentation():
            if V.is_vertex():
                yield V

    @cached_method
    def vertices(self):
        """
        Return all vertices of the face.

        OUTPUT: a tuple of vertices

        EXAMPLES::

            sage: triangle = Polyhedron(vertices=[[1,0],[0,1],[1,1]])
            sage: face = triangle.faces(1)[2]
            sage: face.vertices()
            (A vertex at (0, 1), A vertex at (1, 0))

        TESTS::

            sage: TestSuite(face).run(skip='_test_pickling')
        """
        return tuple(self.vertex_generator())

    @cached_method
    def n_vertices(self):
        """
        Return the number of vertices of the face.

        OUTPUT: integer

        EXAMPLES::

            sage: Q = polytopes.cross_polytope(3)
            sage: face = Q.faces(2)[0]
            sage: face.n_vertices()
            3

        TESTS::

            sage: TestSuite(face).run(skip='_test_pickling')
        """
        return len(self.vertices())

    def ray_generator(self):
        """
        Return a generator for the rays of the face.

        EXAMPLES::

            sage: pi = Polyhedron(ieqs = [[1,1,0],[1,0,1]])
            sage: face = pi.faces(1)[1]
            sage: next(face.ray_generator())
            A ray in the direction (1, 0)

        TESTS::

            sage: TestSuite(face).run(skip='_test_pickling')
        """
        for V in self.ambient_Vrepresentation():
            if V.is_ray():
                yield V

    @cached_method
    def rays(self):
        """
        Return the rays of the face.

        OUTPUT: a tuple of rays

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,0,0,1],[0,0,1,0],[1,1,0,0]])
            sage: face = p.faces(2)[2]
            sage: face.rays()
            (A ray in the direction (1, 0, 0), A ray in the direction (0, 1, 0))
        """
        return tuple(self.ray_generator())

    @cached_method
    def n_rays(self):
        """
        Return the number of rays of the face.

        OUTPUT: integer

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,0,0,1],[0,0,1,0],[1,1,0,0]])
            sage: face = p.faces(2)[0]
            sage: face.n_rays()
            2
        """
        return len(self.rays())

    def line_generator(self):
        """
        Return a generator for the lines of the face.

        EXAMPLES::

            sage: pr = Polyhedron(rays = [[1,0],[-1,0],[0,1]], vertices = [[-1,-1]])
            sage: face = pr.faces(1)[0]
            sage: next(face.line_generator())
            A line in the direction (1, 0)

        TESTS::

            sage: TestSuite(face).run(skip='_test_pickling')
        """
        for V in self.ambient_Vrepresentation():
            if V.is_line():
                yield V

    @cached_method
    def lines(self):
        """
        Return all lines of the face.

        OUTPUT: a tuple of lines

        EXAMPLES::

            sage: p = Polyhedron(rays = [[1,0],[-1,0],[0,1],[1,1]], vertices = [[-2,-2],[2,3]])
            sage: p.lines()
            (A line in the direction (1, 0),)
        """
        return tuple(self.line_generator())

    @cached_method
    def n_lines(self):
        """
        Return the number of lines of the face.

        OUTPUT: integer

        EXAMPLES::

            sage: p = Polyhedron(rays = [[1,0],[-1,0],[0,1],[1,1]], vertices = [[-2,-2],[2,3]])
            sage: p.n_lines()
            1
        """
        return len(self.lines())

    def __richcmp__(self, other, op):
        """
        Compare ``self`` and ``other``.

        INPUT:

        - ``other`` -- anything

        OUTPUT:

        Two faces test equal if and only if they are faces of the same
        (not just isomorphic) polyhedron and their generators have the
        same indices.

        EXAMPLES::

            sage: square = polytopes.hypercube(2)
            sage: f = square.faces(1)
            sage: matrix(4,4, lambda i,j: ZZ(f[i] <= f[j]))
            [1 1 1 0]
            [0 1 0 0]
            [0 1 1 0]
            [1 1 1 1]
            sage: matrix(4,4, lambda i,j: ZZ(f[i] == f[j])) == 1
            True
        """
        if not isinstance(other, PolyhedronFace):
            return NotImplemented
        if self._polyhedron is not other._polyhedron:
            if (self._polyhedron.Vrepresentation() != other._polyhedron.Vrepresentation()
                    or self._polyhedron.Hrepresentation() != other._polyhedron.Hrepresentation()):
                return NotImplemented
        return richcmp(self._ambient_Vrepresentation_indices,
                       other._ambient_Vrepresentation_indices, op)

    def ambient_Hrepresentation(self, index=None):
        r"""
        Return the H-representation objects of the ambient polytope
        defining the face.

        INPUT:

        - ``index`` -- integer or ``None`` (default)

        OUTPUT:

        If the optional argument is not present, a tuple of
        H-representation objects. Each entry is either an inequality
        or an equation.

        If the optional integer ``index`` is specified, the
        ``index``-th element of the tuple is returned.

        EXAMPLES::

            sage: square = polytopes.hypercube(2)
            sage: for face in square.face_lattice():                                    # needs sage.combinat
            ....:     print(face.ambient_Hrepresentation())
            (An inequality (-1, 0) x + 1 >= 0, An inequality (0, -1) x + 1 >= 0,
             An inequality (1, 0) x + 1 >= 0, An inequality (0, 1) x + 1 >= 0)
            (An inequality (-1, 0) x + 1 >= 0, An inequality (0, 1) x + 1 >= 0)
            (An inequality (-1, 0) x + 1 >= 0, An inequality (0, -1) x + 1 >= 0)
            (An inequality (-1, 0) x + 1 >= 0,)
            (An inequality (0, -1) x + 1 >= 0, An inequality (1, 0) x + 1 >= 0)
            (An inequality (0, -1) x + 1 >= 0,)
            (An inequality (1, 0) x + 1 >= 0, An inequality (0, 1) x + 1 >= 0)
            (An inequality (0, 1) x + 1 >= 0,)
            (An inequality (1, 0) x + 1 >= 0,)
            ()
        """
        if index is None:
            return self._ambient_Hrepresentation
        else:
            return self._ambient_Hrepresentation[index]

    def ambient_Vrepresentation(self, index=None):
        r"""
        Return the V-representation objects of the ambient polytope
        defining the face.

        INPUT:

        - ``index`` -- integer or ``None`` (default)

        OUTPUT:

        If the optional argument is not present, a tuple of
        V-representation objects. Each entry is either a vertex, a
        ray, or a line.

        If the optional integer ``index`` is specified, the
        ``index``-th element of the tuple is returned.

        EXAMPLES::

            sage: square = polytopes.hypercube(2)
            sage: for fl in square.face_lattice():                                      # needs sage.combinat
            ....:     print(fl.ambient_Vrepresentation())
            ()
            (A vertex at (1, -1),)
            (A vertex at (1, 1),)
            (A vertex at (1, -1), A vertex at (1, 1))
            (A vertex at (-1, 1),)
            (A vertex at (1, 1), A vertex at (-1, 1))
            (A vertex at (-1, -1),)
            (A vertex at (1, -1), A vertex at (-1, -1))
            (A vertex at (-1, 1), A vertex at (-1, -1))
            (A vertex at (1, -1), A vertex at (1, 1),
             A vertex at (-1, 1), A vertex at (-1, -1))
        """
        if index is None:
            return self._ambient_Vrepresentation
        else:
            return self._ambient_Vrepresentation[index]

    def n_ambient_Hrepresentation(self):
        """
        Return the number of objects that make up the ambient
        H-representation of the polyhedron.

        See also :meth:`ambient_Hrepresentation`.

        OUTPUT: integer

        EXAMPLES::

            sage: p = polytopes.cross_polytope(4)
            sage: face = p.face_lattice()[5]; face                                      # needs sage.combinat
            A 1-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 2 vertices
            sage: face.ambient_Hrepresentation()                                        # needs sage.combinat
            (An inequality (1, -1, 1, -1) x + 1 >= 0,
             An inequality (1, 1, 1, 1) x + 1 >= 0,
             An inequality (1, 1, 1, -1) x + 1 >= 0,
             An inequality (1, -1, 1, 1) x + 1 >= 0)
            sage: face.n_ambient_Hrepresentation()                                      # needs sage.combinat
            4
        """
        return len(self.ambient_Hrepresentation())

    def n_ambient_Vrepresentation(self):
        """
        Return the number of objects that make up the ambient
        V-representation of the polyhedron.

        See also :meth:`ambient_Vrepresentation`.

        OUTPUT: integer

        EXAMPLES::

            sage: p = polytopes.cross_polytope(4)
            sage: face = p.face_lattice()[5]; face                                      # needs sage.combinat
            A 1-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 2 vertices
            sage: face.ambient_Vrepresentation()                                        # needs sage.combinat
            (A vertex at (-1, 0, 0, 0), A vertex at (0, 0, -1, 0))
            sage: face.n_ambient_Vrepresentation()                                      # needs sage.combinat
            2
        """
        return len(self.ambient_Vrepresentation())

    def ambient_H_indices(self):
        """
        Return the indices of the H-representation objects of the
        ambient polyhedron that make up the H-representation of ``self``.

        See also :meth:`ambient_Hrepresentation`.

        OUTPUT: tuple of indices

        EXAMPLES::

            sage: Q = polytopes.cross_polytope(3)
            sage: F = Q.faces(1)
            sage: [f.ambient_H_indices() for f in F]
            [(4, 5),
             (5, 6),
             (4, 7),
             (6, 7),
             (0, 5),
             (3, 4),
             (0, 3),
             (1, 6),
             (0, 1),
             (2, 7),
             (2, 3),
             (1, 2)]
        """
        return self._ambient_Hrepresentation_indices

    def ambient_V_indices(self):
        """
        Return the indices of the V-representation objects of the
        ambient polyhedron that make up the V-representation of ``self``.

        See also :meth:`ambient_Vrepresentation`.

        OUTPUT: tuple of indices

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: F = P.faces(2)
            sage: [f.ambient_V_indices() for f in F]
            [(0, 3, 4, 5),
             (0, 1, 5, 6),
             (4, 5, 6, 7),
             (2, 3, 4, 7),
             (1, 2, 6, 7),
             (0, 1, 2, 3)]
        """
        return self._ambient_Vrepresentation_indices

    def ambient_dim(self):
        r"""
        Return the dimension of the containing polyhedron.

        EXAMPLES::

            sage: P = Polyhedron(vertices = [[1,0,0,0],[0,1,0,0]])
            sage: face = P.faces(1)[0]
            sage: face.ambient_dim()
            4
        """
        return self._polyhedron.ambient_dim()

    @cached_method
    def dim(self):
        """
        Return the dimension of the face.

        OUTPUT: integer

        EXAMPLES::

            sage: fl = polytopes.dodecahedron().face_lattice()                          # needs sage.combinat sage.rings.number_field
            sage: sorted(x.dim() for x in fl)                                           # needs sage.combinat sage.rings.number_field
            [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
              1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3]

        TESTS:

        Check that :issue:`28650` is fixed::

            sage: P = Polyhedron(vertices=[[1,0]], rays=[[1,0], [0,1]])
            sage: P.faces(2)
            (A 2-dimensional face of a Polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 rays,)
        """
        if self.n_ambient_Vrepresentation() == 0:
            return -1
        else:
            origin = self.vertices()[0].vector()
            v_list = [vector(v) - origin for v in
                     self.ambient_Vrepresentation() if v.is_vertex()]
            v_list += [vector(v) for v in self.ambient_Vrepresentation()
                      if v.is_ray() or v.is_line()]
            return matrix(v_list).rank()

    def _repr_(self):
        r"""
        Return a string representation.

        OUTPUT: a string listing the V-representation indices of the face

        EXAMPLES::

            sage: square = polytopes.hypercube(2)
            sage: a_face = list( square.face_lattice() )[8]                             # needs sage.combinat
            sage: a_face.__repr__()                                                     # needs sage.combinat
            'A 1-dimensional face of a Polyhedron in ZZ^2 defined as the convex hull of 2 vertices'
        """
        desc = ''
        desc += 'A ' + repr(self.dim()) + '-dimensional face'
        desc += ' of a Polyhedron in '
        desc += self.polyhedron().parent()._repr_ambient_module()

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

    def polyhedron(self):
        """
        Return the containing polyhedron.

        EXAMPLES::

            sage: P = polytopes.cross_polytope(3); P
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 6 vertices
            sage: face = P.facets()[3]; face
            A 2-dimensional face of a Polyhedron in ZZ^3 defined as the convex hull of 3 vertices
            sage: face.polyhedron()
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 6 vertices
        """
        return self._polyhedron

    ambient = polyhedron

    def ambient_vector_space(self, base_field=None):
        r"""
        Return the ambient vector space.

        It is the ambient free module of the containing polyhedron tensored
        with a field.

        INPUT:

        - ``base_field`` -- a field (default: the fraction field of the base ring)

        EXAMPLES::

            sage: half_plane = Polyhedron(ieqs=[(0,1,0)])
            sage: line = half_plane.faces(1)[0]; line
            A 1-dimensional face of a
             Polyhedron in QQ^2 defined as the convex hull of 1 vertex and 1 line
            sage: line.ambient_vector_space()
            Vector space of dimension 2 over Rational Field
            sage: line.ambient_vector_space(AA)                                         # needs sage.rings.number_field
            Vector space of dimension 2 over Algebraic Real Field
        """
        return self.polyhedron().ambient_vector_space(base_field=base_field)

    def is_relatively_open(self):
        r"""
        Return whether ``self`` is relatively open.

        OUTPUT: boolean

        EXAMPLES::

            sage: half_plane = Polyhedron(ieqs=[(0,1,0)])
            sage: line = half_plane.faces(1)[0]; line
            A 1-dimensional face of a
             Polyhedron in QQ^2 defined as the convex hull of 1 vertex and 1 line
            sage: line.is_relatively_open()
            True
        """
        return self.as_polyhedron().is_relatively_open()

    def is_compact(self):
        r"""
        Return whether ``self`` is compact.

        OUTPUT: boolean

        EXAMPLES::

            sage: half_plane = Polyhedron(ieqs=[(0,1,0)])
            sage: line = half_plane.faces(1)[0]; line
            A 1-dimensional face of a
             Polyhedron in QQ^2 defined as the convex hull of 1 vertex and 1 line
            sage: line.is_compact()
            False
        """
        return not any(V.is_ray() or V.is_line()
                       for V in self.ambient_Vrepresentation())

    @cached_method
    def as_polyhedron(self, **kwds):
        """
        Return the face as an independent polyhedron.

        OUTPUT: a polyhedron

        EXAMPLES::

            sage: P = polytopes.cross_polytope(3);  P
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 6 vertices
            sage: face = P.faces(2)[3]; face
            A 2-dimensional face of a
             Polyhedron in ZZ^3 defined as the convex hull of 3 vertices
            sage: face.as_polyhedron()
            A 2-dimensional polyhedron in ZZ^3 defined as the convex hull of 3 vertices

            sage: P.intersection(face.as_polyhedron()) == face.as_polyhedron()
            True
        """
        P = self._polyhedron
        parent = P.parent()
        Vrep = (self.vertices(), self.rays(), self.lines())
        result = P.__class__(parent, Vrep, None)
        if any(kwds.get(kwd) is not None
               for kwd in ('base_ring', 'backend')):
            from .constructor import Polyhedron
            return Polyhedron(result, **kwds)
        return result

    def _some_elements_(self):
        r"""
        Generate some points of ``self``.

        If ``self`` is empty, no points are generated; no exception will be raised.

        EXAMPLES::

            sage: P = polytopes.cross_polytope(3);  P
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 6 vertices
            sage: face = P.faces(2)[3]; face
            A 2-dimensional face of a
             Polyhedron in ZZ^3 defined as the convex hull of 3 vertices
            sage: face.as_polyhedron().vertices()
            (A vertex at (0, -1, 0), A vertex at (0, 0, -1), A vertex at (1, 0, 0))
            sage: face.an_element()              # indirect doctest
            (1/3, -1/3, -1/3)
            sage: face.some_elements()           # indirect doctest
            [(1/3, -1/3, -1/3), (0, -1, 0), (0, -1/2, -1/2), (1/2, -1/4, -1/4)]
        """
        yield from self.as_polyhedron().some_elements()

    def contains(self, point):
        """
        Test whether the polyhedron contains the given ``point``.

        INPUT:

        - ``point`` -- a point or its coordinates

        EXAMPLES::

            sage: half_plane = Polyhedron(ieqs=[(0,1,0)])
            sage: line = half_plane.faces(1)[0]; line
            A 1-dimensional face of a
             Polyhedron in QQ^2 defined as the convex hull of 1 vertex and 1 line
            sage: line.contains([0, 1])
            True

        As a shorthand, one may use the usual ``in`` operator::

            sage: [5, 7] in line
            False
        """
        # preprocess in the same way as Polyhedron_base.contains
        try:
            p = vector(point)
        except TypeError:  # point not iterable or no common ring for elements
            if len(point) > 0:
                return False
            else:
                p = vector(self.polyhedron().base_ring(), [])

        if len(p) != self.ambient_dim():
            return False

        if not self.polyhedron().contains(p):
            return False

        for H in self.ambient_Hrepresentation():
            if H.eval(p) != 0:
                return False
        return True

    __contains__ = contains

    @cached_method
    def normal_cone(self, direction='outer'):
        """
        Return the polyhedral cone consisting of normal vectors to
        hyperplanes supporting ``self``.

        INPUT:

        - ``direction`` -- string (default: ``'outer'``); the direction in
          which to consider the normals. The other allowed option is
          ``'inner'``.

        OUTPUT: a polyhedron

        EXAMPLES::

            sage: p = Polyhedron(vertices=[[1,2], [2,1], [-2,2], [-2,-2], [2,-2]])
            sage: for v in p.face_generator(0):
            ....:     vect = v.vertices()[0].vector()
            ....:     nc = v.normal_cone().rays_list()
            ....:     print("{} has outer normal cone spanned by {}".format(vect,nc))
            ....:
            (2, 1) has outer normal cone spanned by [[1, 0], [1, 1]]
            (1, 2) has outer normal cone spanned by [[0, 1], [1, 1]]
            (2, -2) has outer normal cone spanned by [[0, -1], [1, 0]]
            (-2, -2) has outer normal cone spanned by [[-1, 0], [0, -1]]
            (-2, 2) has outer normal cone spanned by [[-1, 0], [0, 1]]

            sage: for v in p.face_generator(0):
            ....:     vect = v.vertices()[0].vector()
            ....:     nc = v.normal_cone(direction='inner').rays_list()
            ....:     print("{} has inner normal cone spanned by {}".format(vect,nc))
            ....:
            (2, 1) has inner normal cone spanned by [[-1, -1], [-1, 0]]
            (1, 2) has inner normal cone spanned by [[-1, -1], [0, -1]]
            (2, -2) has inner normal cone spanned by [[-1, 0], [0, 1]]
            (-2, -2) has inner normal cone spanned by [[0, 1], [1, 0]]
            (-2, 2) has inner normal cone spanned by [[0, -1], [1, 0]]

        The function works for polytopes that are not full-dimensional::

            sage: p = polytopes.permutahedron(3)
            sage: f1 = p.faces(0)[0]
            sage: f2 = p.faces(1)[0]
            sage: f3 = p.faces(2)[0]
            sage: f1.normal_cone()
            A 3-dimensional polyhedron in ZZ^3 defined as
             the convex hull of 1 vertex, 2 rays, 1 line
            sage: f2.normal_cone()
            A 2-dimensional polyhedron in ZZ^3 defined as
             the convex hull of 1 vertex, 1 ray, 1 line
            sage: f3.normal_cone()
            A 1-dimensional polyhedron in ZZ^3 defined as
             the convex hull of 1 vertex and 1 line

        Normal cones are only defined for non-empty faces::

            sage: f0 = p.faces(-1)[0]
            sage: f0.normal_cone()
            Traceback (most recent call last):
            ...
            ValueError: the empty face does not have a normal cone
        """
        if self.dim() == -1:
            raise ValueError("the empty face does not have a normal cone")
        elif direction not in ['outer','inner']:
            raise ValueError("the direction should be either 'outer' or 'inner'")
        rays = []
        lines = []
        for facet in self.ambient_Hrepresentation():
            if facet.is_equation():
                lines += [facet.A()]
            elif direction == 'outer':
                rays += [-facet.A()]
            else:  # 'inner'
                rays += [facet.A()]
        parent = self.polyhedron().parent()
        origin = self.polyhedron().ambient_space().zero()
        return parent.element_class(parent, [[origin], rays, lines], None)

    @cached_method
    def affine_tangent_cone(self):
        """
        Return the affine tangent cone of ``self`` as a polyhedron.

        It is equal to the sum of ``self`` and the cone of feasible directions
        at any point of the relative interior of ``self``.

        OUTPUT: a polyhedron

        EXAMPLES::

            sage: half_plane_in_space = Polyhedron(ieqs=[(0,1,0,0)], eqns=[(0,0,0,1)])
            sage: line = half_plane_in_space.faces(1)[0]; line
            A 1-dimensional face of a
             Polyhedron in QQ^3 defined as the convex hull of 1 vertex and 1 line
            sage: T_line = line.affine_tangent_cone()
            sage: T_line == half_plane_in_space
            True

            sage: c = polytopes.cube()
            sage: edge = min(c.faces(1))
            sage: edge.vertices()
            (A vertex at (1, -1, -1), A vertex at (1, 1, -1))
            sage: T_edge = edge.affine_tangent_cone()
            sage: T_edge.Vrepresentation()
            (A line in the direction (0, 1, 0),
             A ray in the direction (0, 0, 1),
             A vertex at (1, 0, -1),
             A ray in the direction (-1, 0, 0))

        TESTS:

        Check that :issue:`32658` is fixed::

            sage: P = polytopes.hypercube(2)
            sage: P.faces(-1)[0].affine_tangent_cone()
            Traceback (most recent call last):
            ...
            ValueError: affine tangent cone of the empty face not defined
        """
        if self.dim() == -1:
            raise ValueError("affine tangent cone of the empty face not defined")
        parent = self.polyhedron().parent()
        new_ieqs = [H for H in self.ambient_Hrepresentation()
                    if H.is_inequality()]
        new_eqns = [H for H in self.ambient_Hrepresentation()
                    if H.is_equation()]
        return parent.element_class(parent, None, [new_ieqs, new_eqns])

    @cached_method
    def stacking_locus(self):
        """
        Return the polyhedron containing the points that sees every facet
        containing ``self``.

        OUTPUT: a polyhedron

        EXAMPLES::

            sage: cp = polytopes.cross_polytope(4)
            sage: facet = cp.facets()[0]
            sage: facet.stacking_locus().vertices()
            (A vertex at (1/2, 1/2, 1/2, 1/2),
             A vertex at (1, 0, 0, 0),
             A vertex at (0, 0, 0, 1),
             A vertex at (0, 0, 1, 0),
             A vertex at (0, 1, 0, 0))
            sage: face = cp.faces(2)[0]
            sage: face.stacking_locus().vertices()
            (A vertex at (0, 1, 0, 0),
             A vertex at (0, 0, 1, 0),
             A vertex at (1, 0, 0, 0),
             A vertex at (1, 1, 1, 0),
             A vertex at (1/2, 1/2, 1/2, 1/2),
             A vertex at (1/2, 1/2, 1/2, -1/2))
        """
        # Taking all facets that contain the face
        if self.dim() == self.polyhedron().dim() - 1:
            face_star = set([self.ambient_Hrepresentation()[-1]])
        else:
            face_star = set(facet for facet in self.ambient_Hrepresentation() if facet.is_inequality()
                            if all(not facet.interior_contains(x) for x in self.vertices()))

        neighboring_facets = set()
        for facet in face_star:
            for neighbor_facet in facet.neighbors():
                if neighbor_facet not in face_star:
                    neighboring_facets.add(neighbor_facet)

        # Create the polyhedron where we can put the new vertex
        locus_ieqs = [facet.vector() for facet in neighboring_facets]
        locus_ieqs += [-facet.vector() for facet in face_star]
        locus_eqns = self.polyhedron().equations_list()
        parent = self.polyhedron().parent().change_ring(self.polyhedron().base_ring().fraction_field())

        return parent.element_class(parent, None, [locus_ieqs, locus_eqns])

def combinatorial_face_to_polyhedral_face(polyhedron, combinatorial_face):
    r"""
    Convert a combinatorial face to a face of a polyhedron.

    INPUT:

    - ``polyhedron`` -- a polyhedron containing ``combinatorial_face``
    - ``combinatorial_face`` -- a :class:`CombinatorialFace`

    OUTPUT: a :class:`PolyhedronFace`

    EXAMPLES::

        sage: from sage.geometry.polyhedron.face import combinatorial_face_to_polyhedral_face
        sage: P = polytopes.simplex()
        sage: C = P.combinatorial_polyhedron()
        sage: it = C.face_iter()
        sage: comb_face = next(it)
        sage: combinatorial_face_to_polyhedral_face(P, comb_face)
        A 2-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 3 vertices

    TESTS:

    Making sure that backends do not change their order of inequalites/equations
    without applying the changes to this method::

        sage: polytopes.simplex(backend='field').equations()[0].index()
        4
        sage: polytopes.simplex(backend='ppl').equations()[0].index()
        0
        sage: polytopes.simplex(backend='cdd').equations()[0].index()
        4
        sage: polytopes.simplex(backend='normaliz').equations()[0].index()   # optional - pynormaliz
        4
        sage: polytopes.simplex(backend='polymake').equations()[0].index()   # optional - jupymake
        4
    """
    V_indices = combinatorial_face.ambient_V_indices()
    n_equations = polyhedron.n_equations()

    if polyhedron.backend() in ('ppl',):
        # Equations before inequalities in Hrep.
        H_indices = tuple(range(n_equations))
        H_indices += tuple(x+n_equations for x in combinatorial_face.ambient_H_indices(add_equations=False))
    elif polyhedron.backend() in ('normaliz', 'cdd', 'field', 'number_field', 'polymake'):
        # Equations after the inequalities in Hrep.
        n_ieqs = polyhedron.n_inequalities()
        H_indices = tuple(combinatorial_face.ambient_H_indices(add_equations=False))
        H_indices += tuple(range(n_ieqs, n_ieqs + n_equations))
    else:
        raise NotImplementedError("unknown backend")

    if polyhedron.dimension() == 0:
        # Taking care of a special case:
        # In this case the face lattice has a coatom,
        # but the polyhedron does not have a facet
        # (a facet is defined to be non-empty).

        # More important, there is no inequality for that coatom.
        # So the above would produce an index error.
        # Instead, any case of the 0-dimensional polyhedron
        # satisfies all of the equations.
        H_indices = tuple(range(n_equations))

    return PolyhedronFace(polyhedron, V_indices, H_indices)
