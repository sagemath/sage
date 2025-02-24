"""
Functions for plotting polyhedra
"""

########################################################################
#       Copyright (C) 2008 Marshall Hampton <hamptonio@gmail.com>
#       Copyright (C) 2011 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
########################################################################

from math import pi

from sage.misc.lazy_import import lazy_import
from sage.structure.sage_object import SageObject
from sage.modules.free_module_element import vector
from sage.matrix.constructor import matrix, identity_matrix
from sage.matrix.special import diagonal_matrix
from sage.misc.functional import norm
from sage.misc.latex import LatexExpr
from sage.structure.sequence import Sequence

lazy_import("sage.plot.all", ["Graphics", "point2d", "line2d", "arrow", "polygon2d", "rainbow"])
lazy_import("sage.plot.plot3d.all", ["point3d", "line3d", "arrow3d", "polygons3d"])
lazy_import("sage.plot.plot3d.transform", "rotate_arbitrary")
lazy_import("sage.plot.plot3d.texture", "Texture")


#############################################################


def cyclic_sort_vertices_2d(Vlist):
    """
    Return the vertices/rays in cyclic order if possible.

    .. NOTE::

        This works if and only if each vertex/ray is adjacent to exactly
        two others. For example, any 2-dimensional polyhedron satisfies
        this.

    See
    :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.vertex_adjacency_matrix`
    for a discussion of "adjacent".

    EXAMPLES::

        sage: from sage.geometry.polyhedron.plot import cyclic_sort_vertices_2d
        sage: square = Polyhedron([[1,0],[-1,0],[0,1],[0,-1]])
        sage: vertices = [v for v in square.vertex_generator()]
        sage: vertices
        [A vertex at (-1, 0),
         A vertex at (0, -1),
         A vertex at (0, 1),
         A vertex at (1, 0)]
        sage: cyclic_sort_vertices_2d(vertices)
        [A vertex at (1, 0),
         A vertex at (0, -1),
         A vertex at (-1, 0),
         A vertex at (0, 1)]

    Rays are allowed, too::

        sage: P = Polyhedron(vertices=[(0, 1), (1, 0), (2, 0), (3, 0), (4, 1)], rays=[(0,1)])
        sage: P.adjacency_matrix()
        [0 1 0 1 0]
        [1 0 1 0 0]
        [0 1 0 0 1]
        [1 0 0 0 1]
        [0 0 1 1 0]
        sage: cyclic_sort_vertices_2d(P.Vrepresentation())
        [A vertex at (3, 0),
         A vertex at (1, 0),
         A vertex at (0, 1),
         A ray in the direction (0, 1),
         A vertex at (4, 1)]

        sage: P = Polyhedron(vertices=[(0, 1), (1, 0), (2, 0), (3, 0), (4, 1)], rays=[(0,1), (1,1)])
        sage: P.adjacency_matrix()
        [0 1 0 0 0]
        [1 0 1 0 0]
        [0 1 0 0 1]
        [0 0 0 0 1]
        [0 0 1 1 0]
        sage: cyclic_sort_vertices_2d(P.Vrepresentation())
        [A ray in the direction (1, 1),
         A vertex at (3, 0),
         A vertex at (1, 0),
         A vertex at (0, 1),
         A ray in the direction (0, 1)]

        sage: P = Polyhedron(vertices=[(1,2)], rays=[(0,1)], lines=[(1,0)])
        sage: P.adjacency_matrix()
        [0 0 1]
        [0 0 0]
        [1 0 0]
        sage: cyclic_sort_vertices_2d(P.Vrepresentation())
        [A vertex at (0, 2),
         A line in the direction (1, 0),
         A ray in the direction (0, 1)]
    """
    if not Vlist:
        return Vlist
    Vlist = list(Vlist)
    result = []
    adjacency_matrix = Vlist[0].polyhedron().vertex_adjacency_matrix()

    # Any object in Vlist has 0,1, or 2 adjacencies. Break into connected chains:
    chain = [Vlist.pop()]
    while Vlist:
        first_index = chain[0].index()
        last_index = chain[-1].index()
        for v in Vlist:
            v_index = v.index()
            if adjacency_matrix[last_index, v_index] == 1:
                chain = chain + [v]
                Vlist.remove(v)
                break
            if adjacency_matrix[first_index, v.index()] == 1:
                chain = [v] + chain
                Vlist.remove(v)
                break
        else:
            result += chain
            chain = [Vlist.pop()]
    result += chain
    return result

#########################################################################


def projection_func_identity(x):
    """
    The identity projection.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.plot import projection_func_identity
        sage: projection_func_identity((1,2,3))
        [1, 2, 3]
    """
    return list(x)


class ProjectionFuncStereographic:
    """
    The stereographic (or perspective) projection onto a codimension-1 linear
    subspace with respect to a sphere centered at the origin.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.plot import ProjectionFuncStereographic
        sage: cube = polytopes.hypercube(3).vertices()
        sage: proj = ProjectionFuncStereographic([1.2, 3.4, 5.6])
        sage: ppoints = [proj(vector(x)) for x in cube]
        sage: ppoints[5]
        (-0.0918273..., -0.036375...)
    """
    def __init__(self, projection_point):
        """
        Create a stereographic projection function.

        INPUT:

        - ``projection_point`` -- list of coordinates in the
          appropriate dimension, which is the point projected from

        EXAMPLES::

            sage: from sage.geometry.polyhedron.plot import ProjectionFuncStereographic
            sage: proj = ProjectionFuncStereographic([1.0,1.0])
            sage: proj.__init__([1.0,1.0])
            sage: proj.house
            [ 0.7071067811... -0.7071067811...]
            [ 0.7071067811...  0.7071067811...]
            sage: TestSuite(proj).run(skip='_test_pickling')
        """
        from sage.rings.real_double import RDF

        self.projection_point = vector(projection_point)
        self.dim = self.projection_point.degree()

        pproj = vector(RDF, self.projection_point)
        self.psize = norm(pproj)
        if (self.psize).is_zero():
            raise ValueError("projection direction must be a nonzero vector.")
        v = vector(RDF, [0.0] * (self.dim - 1) + [-self.psize]) - pproj
        polediff = matrix(RDF, v).transpose()
        denom = RDF((polediff.transpose() * polediff)[0][0])
        if denom.is_zero():
            self.house = identity_matrix(RDF, self.dim)
        else:
            house = identity_matrix(RDF, self.dim) \
            - 2*polediff*polediff.transpose()/denom   # Householder reflector
            # Make it preserve orientation (chirality):
            self.house = diagonal_matrix(RDF, [1] * (self.dim - 1) + [-1]) * house

    def __call__(self, x):
        """
        Action of the stereographic projection.

        INPUT:

        - ``x`` -- a vector or anything convertible to a vector

        OUTPUT:

        First reflects ``x`` with a Householder reflection which takes
        the projection point to ``(0,...,0,self.psize)`` where
        ``psize`` is the length of the projection point, and then
        dilates by ``1/(zdiff)`` where ``zdiff`` is the difference
        between the last coordinate of ``x`` and ``psize``.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.plot import ProjectionFuncStereographic
            sage: proj = ProjectionFuncStereographic([1.0,1.0])
            sage: proj.__call__(vector([1,2]))
            (1.0000000000000002)
            sage: proj = ProjectionFuncStereographic([2.0,1.0])
            sage: proj.__call__(vector([1,2]))  # abs tol 1e-14
            (-2.999999999999997)
            sage: proj = ProjectionFuncStereographic([0,0,2])
            sage: proj.__call__(vector([0,0,1]))
            (0.0, 0.0)
            sage: proj.__call__(vector([1,0,0]))
            (0.5, 0.0)
        """
        from sage.rings.real_double import RDF

        img = self.house * x
        denom = self.psize - img[self.dim - 1]
        if denom.is_zero():
            raise ValueError('Point cannot coincide with '
                             'coordinate singularity at ' + repr(x))
        return vector(RDF, [img[i] / denom for i in range(self.dim - 1)])


class ProjectionFuncSchlegel:
    """
    The Schlegel projection from the given input point.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.plot import ProjectionFuncSchlegel
        sage: fcube = polytopes.hypercube(4)
        sage: facet = fcube.facets()[0]
        sage: proj = ProjectionFuncSchlegel(facet,[0,-1.5,0,0])
        sage: proj([0,0,0,0])[0]
        1.0
    """
    def __init__(self, facet, projection_point):
        """
        Initialize the projection.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.plot import ProjectionFuncSchlegel
            sage: fcube = polytopes.hypercube(4)
            sage: facet = fcube.facets()[0]
            sage: proj = ProjectionFuncSchlegel(facet,[0,-1.5,0,0])
            sage: proj.facet
            A 3-dimensional face of a Polyhedron in ZZ^4 defined as the convex hull of 8 vertices
            sage: proj.A
            [1.0 0.0 0.0]
            [0.0 0.0 0.0]
            [0.0 0.0 1.0]
            [0.0 1.0 0.0]
            sage: proj.b
            (1.0, 1.0, 1.0)
            sage: proj.projection_point
            (0.0, -1.5, 0.0, 0.0)
            sage: proj([-1,1,1,1])
            (0.8, 1.2, 1.2)
            sage: proj([1,1,1,1])
            (1.2, 1.2, 1.2)
            sage: proj([1,-1,1,1])
            (2.0, 2.0, 2.0)
            sage: proj([1,-1,-1,1])
            (2.0, 2.0, 0.0)
            sage: TestSuite(proj).run(skip='_test_pickling')
        """
        from sage.rings.real_double import RDF

        self.facet = facet
        ineq = [h for h in facet.ambient_Hrepresentation() if h.is_inequality()][0]
        self.full_A = ineq.A()
        self.full_b = ineq.b()
        A, b = self.facet.as_polyhedron().affine_hull_projection(as_affine_map=True, orthonormal=True, extend=True)
        self.A = A.change_ring(RDF).matrix()
        self.b = b.change_ring(RDF)
        self.projection_point = vector(RDF, projection_point)

    def __call__(self, x):
        """
        Apply the projection to a vector.

        - ``x`` -- a vector or anything convertible to a vector

        EXAMPLES::

            sage: from sage.geometry.polyhedron.plot import ProjectionFuncSchlegel
            sage: fcube = polytopes.hypercube(4)
            sage: facet = fcube.facets()[0]
            sage: proj = ProjectionFuncSchlegel(facet,[0,-1.5,0,0])
            sage: proj.__call__([0,0,0,0])
            (1.0, 1.0, 1.0)
        """
        # The intersection of the segment with the facet
        # See Ziegler's "Lectures on Polytopes" p.133
        vx = vector(x)
        z = (self.full_b)
        a = -(self.full_A)
        y = self.projection_point
        preimage = y + ((z-a*y)/(a*vx-a*y))*(vx - y)
        # The transformation matrix acts on the right:
        return preimage*self.A + self.b

#########################################################################


class Projection(SageObject):
    """
    The projection of a :class:`Polyhedron`.

    This class keeps track of the necessary data to plot the input
    polyhedron.
    """

    def __init__(self, polyhedron, proj=projection_func_identity):
        """
        Initialize the projection of a Polyhedron() object.

        INPUT:

        - ``polyhedron`` -- a ``Polyhedron()`` object

        - ``proj`` -- a projection function for the points

        .. NOTE::

            Once initialized, the polyhedral data is fixed. However, the
            projection can be changed later on.

        EXAMPLES::

            sage: # needs sage.groups
            sage: p = polytopes.icosahedron(exact=False)
            sage: from sage.geometry.polyhedron.plot import Projection
            sage: Projection(p)
            The projection of a polyhedron into 3 dimensions
            sage: def pr_12(x): return [x[1],x[2]]
            sage: Projection(p, pr_12)
            The projection of a polyhedron into 2 dimensions
            sage: Projection(p,  lambda x: [x[1],x[2]] )   # another way of doing the same projection
            The projection of a polyhedron into 2 dimensions
            sage: _.plot()   # plot of the projected icosahedron in 2d                  # needs sage.plot
            Graphics object consisting of 51 graphics primitives
            sage: proj = Projection(p)
            sage: proj.stereographic([1,2,3])
            The projection of a polyhedron into 2 dimensions
            sage: proj.plot()                                                           # needs sage.plot
            Graphics object consisting of 51 graphics primitives
            sage: TestSuite(proj).run(skip='_test_pickling')
        """
        self.parent_polyhedron = polyhedron
        self.coords = Sequence([])
        self.points = Sequence([])
        self.lines = Sequence([])
        self.arrows = Sequence([])
        self.polygons = Sequence([])
        self.polyhedron_ambient_dim = polyhedron.ambient_dim()
        self.polyhedron_dim = polyhedron.dim()

        if polyhedron.ambient_dim() == 2:
            self._init_from_2d(polyhedron)
        elif polyhedron.ambient_dim() == 3:
            self._init_from_3d(polyhedron)
        else:
            self._init_points(polyhedron)
            self._init_lines_arrows(polyhedron)

        self.coords.set_immutable()
        self.points.set_immutable()
        self.lines.set_immutable()
        self.arrows.set_immutable()
        self.polygons.set_immutable()

        self(proj)

    def _repr_(self):
        """
        Return a string describing the projection.

        EXAMPLES::

            sage: p = polytopes.hypercube(3)
            sage: from sage.geometry.polyhedron.plot import Projection
            sage: proj = Projection(p)
            sage: print(proj._repr_())
            The projection of a polyhedron into 3 dimensions
        """
        s = 'The projection of a polyhedron into '
        s += repr(self.dimension) + ' dimensions'
        return s

    def __call__(self, proj=projection_func_identity):
        """
        Apply a projection.

        EXAMPLES::

            sage: # needs sage.groups
            sage: p = polytopes.icosahedron(exact=False)
            sage: from sage.geometry.polyhedron.plot import Projection
            sage: pproj = Projection(p)
            sage: from sage.geometry.polyhedron.plot import ProjectionFuncStereographic
            sage: pproj_stereo = pproj.__call__(proj=ProjectionFuncStereographic([1, 2, 3]))
            sage: sorted(pproj_stereo.polygons)
            [[2, 0, 9],
             [3, 1, 10],
             [4, 0, 8],
             ...
             [11, 1, 3],
             [11, 3, 7],
             [11, 7, 9]]
        """
        self.transformed_coords = Sequence([proj(p) for p in self.coords])
        self._init_dimension()
        return self

    def identity(self):
        """
        Return the identity projection of the polyhedron.

        EXAMPLES::

            sage: # needs sage.groups
            sage: p = polytopes.icosahedron(exact=False)
            sage: from sage.geometry.polyhedron.plot import Projection
            sage: pproj = Projection(p)
            sage: ppid = pproj.identity()
            sage: ppid.dimension
            3
        """
        return self(projection_func_identity)

    def stereographic(self, projection_point=None):
        r"""
        Return the stereographic projection.

        INPUT:

        - ``projection_point`` -- the projection point. This must be
          distinct from the polyhedron's vertices. Default is `(1,0,\dots,0)`.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.plot import Projection
            sage: proj = Projection(polytopes.buckyball()); proj  # long time
            The projection of a polyhedron into 3 dimensions
            sage: proj.stereographic([5,2,3]).plot()    # long time                     # needs sage.plot
            Graphics object consisting of 123 graphics primitives
            sage: Projection(polytopes.twenty_four_cell()).stereographic([2,0,0,0])
            The projection of a polyhedron into 3 dimensions
        """
        if projection_point is None:
            projection_point = [1] + [0]*(self.polyhedron_ambient_dim-1)
        return self(ProjectionFuncStereographic(projection_point))

    def schlegel(self, facet=None, position=None):
        r"""
        Return the Schlegel projection.

        * The facet is orthonormally transformed into its affine hull.

        * The position specifies a point coming out of the barycenter of the
          facet from which the other vertices will be projected into the facet.

        INPUT:

        - ``facet`` -- a PolyhedronFace; the facet into which the Schlegel
          diagram is created. The default is the first facet.

        - ``position`` -- a positive number. Determines a relative distance
          from the barycenter of ``facet``. A value close to 0 will place the
          projection point close to the facet and a large value further away.
          If the given value is too large, an error is returned.
          If no position is given, it takes the midpoint of the possible
          point of views along a line spanned by the barycenter of the facet
          and a valid point outside the facet.

        EXAMPLES::

            sage: cube4 = polytopes.hypercube(4)
            sage: from sage.geometry.polyhedron.plot import Projection
            sage: Projection(cube4).schlegel()
            The projection of a polyhedron into 3 dimensions
            sage: _.plot()                                                              # needs sage.plot
            Graphics3d Object

        The 4-cube with a truncated vertex seen into the resulting tetrahedron
        facet::

            sage: tcube4 = cube4.face_truncation(cube4.faces(0)[0])
            sage: tcube4.facets()[4]
            A 3-dimensional face of a Polyhedron in QQ^4 defined as the convex hull of 4 vertices
            sage: into_tetra = Projection(tcube4).schlegel(tcube4.facets()[4])          # needs sage.symbolic
            sage: into_tetra.plot()                                                     # needs sage.plot sage.symbolic
            Graphics3d Object

        Taking a larger value for the position changes the image::

            sage: into_tetra_far = Projection(tcube4).schlegel(tcube4.facets()[4], 4)   # needs sage.symbolic
            sage: into_tetra_far.plot()                                                 # needs sage.plot sage.symbolic
            Graphics3d Object

        A value which is too large or negative give a projection point that
        sees more than one facet resulting in a error::

            sage: Projection(tcube4).schlegel(tcube4.facets()[4], 5)
            Traceback (most recent call last):
            ...
            ValueError: the chosen position is too large
            sage: Projection(tcube4).schlegel(tcube4.facets()[4], -1)
            Traceback (most recent call last):
            ...
            ValueError: 'position' should be a positive number
        """
        from sage.geometry.polyhedron.face import PolyhedronFace
        from sage.rings.integer_ring import ZZ
        if facet is None:
            facet = self.parent_polyhedron.facets()[0]
        elif not isinstance(facet, PolyhedronFace):
            raise TypeError("{} should be a PolyhedronFace of {}".format(facet, self))
        elif facet.dim() != self.parent_polyhedron.dim() - 1:
            raise ValueError("The face should be a facet of the polyhedron")
        if position is not None and position <= 0:
            raise ValueError("'position' should be a positive number")

        barycenter = ZZ.one() * sum([v.vector() for v in facet.vertices()]) / len(facet.vertices())
        locus_polyhedron = facet.stacking_locus()
        repr_point = locus_polyhedron.representative_point()
        if position is None:
            # Figure out a somehow canonical point of view inside the locus
            # polyhedron
            from sage.geometry.polyhedron.constructor import Polyhedron
            the_ray = Polyhedron(vertices=[barycenter],
                                 rays=[repr_point - barycenter],
                                 backend=locus_polyhedron.backend()) & locus_polyhedron
            projection_point = the_ray.representative_point()
        else:
            projection_point = (1-position)*barycenter + position*repr_point
            if not locus_polyhedron.relative_interior_contains(projection_point):
                raise ValueError("the chosen position is too large")

        return self(ProjectionFuncSchlegel(facet, projection_point))

    def coord_index_of(self, v):
        """
        Convert a coordinate vector to its internal index.

        EXAMPLES::

            sage: p = polytopes.hypercube(3)
            sage: proj = p.projection()
            sage: proj.coord_index_of(vector((1,1,1)))
            2
        """
        try:
            return self.coords.index(v)
        except ValueError:
            self.coords.append(v)
            return len(self.coords)-1

    def coord_indices_of(self, v_list):
        """
        Convert list of coordinate vectors to the corresponding list
        of internal indices.

        EXAMPLES::

            sage: p = polytopes.hypercube(3)
            sage: proj = p.projection()
            sage: proj.coord_indices_of([vector((1,1,1)), vector((1,-1,1))])
            [2, 3]
        """
        return [self.coord_index_of(v) for v in v_list]

    def coordinates_of(self, coord_index_list):
        """
        Given a list of indices, return the projected coordinates.

        EXAMPLES::

            sage: p = polytopes.simplex(4, project=True).projection()
            sage: p.coordinates_of([1])
            [[-0.7071067812, 0.4082482905, 0.2886751346, 0.2236067977]]
        """
        return [self.transformed_coords[i] for i in coord_index_list]

    def _init_dimension(self):
        """
        Internal function: Initialize from polyhedron with
        projected coordinates. Must always be called after
        a coordinate projection.

        TESTS::

            sage: p = polytopes.simplex(2, project=True).projection()
            sage: test = p._init_dimension()
            sage: p.plot.__doc__ == p.render_2d.__doc__
            True
        """
        if self.transformed_coords:
            self.dimension = len(self.transformed_coords[0])
        else:
            self.dimension = 0
        if self.dimension == 0:
            self.plot = self.render_0d
        elif self.dimension == 1:
            self.plot = self.render_1d
        elif self.dimension == 2:
            self.plot = self.render_2d
        elif self.dimension == 3:
            self.plot = self.render_3d
        else:
            try:
                del self.plot
            except AttributeError:
                pass

    def _init_from_2d(self, polyhedron):
        """
        Internal function: Initialize from polyhedron in
        2-dimensional space. The polyhedron could be lower
        dimensional.

        TESTS::

            sage: p = Polyhedron(vertices=[[0,0],[0,1],[1,0],[1,1]])
            sage: proj = p.projection()
            sage: [proj.coordinates_of([i]) for i in proj.points]
            [[[0, 0]], [[0, 1]], [[1, 0]], [[1, 1]]]
            sage: proj._init_from_2d
            <bound method Projection._init_from_2d of The projection
            of a polyhedron into 2 dimensions>
        """
        assert polyhedron.ambient_dim() == 2, "Requires polyhedron in 2d"
        self.dimension = 2
        self._init_points(polyhedron)
        self._init_lines_arrows(polyhedron)
        self._init_area_2d(polyhedron)

    def _init_from_3d(self, polyhedron):
        """
        Internal function: Initialize from polyhedron in
        3-dimensional space. The polyhedron could be
        lower dimensional.

        TESTS::

            sage: p = Polyhedron(vertices=[[0,0,1],[0,1,2],[1,0,3],[1,1,5]])
            sage: proj = p.projection()
            sage: [proj.coordinates_of([i]) for i in proj.points]
            [[[0, 0, 1]], [[0, 1, 2]], [[1, 0, 3]], [[1, 1, 5]]]
            sage: proj._init_from_3d
            <bound method Projection._init_from_3d of The projection
            of a polyhedron into 3 dimensions>
        """
        assert polyhedron.ambient_dim() == 3, "Requires polyhedron in 3d"
        self.dimension = 3
        self._init_points(polyhedron)
        self._init_lines_arrows(polyhedron)
        self._init_solid_3d(polyhedron)

    def _init_points(self, polyhedron):
        """
        Internal function: Initialize points (works in arbitrary
        dimensions).

        TESTS::

            sage: p = polytopes.hypercube(2)
            sage: pp = p.projection()
            sage: del pp.points
            sage: pp.points = Sequence([])
            sage: pp._init_points(p)
            sage: pp.points
            [0, 1, 2, 3]
        """
        for v in polyhedron.vertex_generator():
            self.points.append(self.coord_index_of(v.vector()))

    def _init_lines_arrows(self, polyhedron):
        """
        Internal function: Initialize compact and non-compact edges
        (works in arbitrary dimensions).

        TESTS::

            sage: p = Polyhedron(ieqs=[[1, 0, 0, 1], [1, 1, 0, 0]])
            sage: pp = p.projection()
            sage: pp.arrows
            [[0, 1], [0, 2], [0, 3], [0, 4]]
            sage: del pp.arrows
            sage: pp.arrows = Sequence([])
            sage: pp._init_lines_arrows(p)
            sage: pp.arrows
            [[0, 1], [0, 2], [0, 3], [0, 4]]

        We check that :issue:`31802` is fixed::

            sage: x = Polyhedron(lines=[(1, 0, 0), (0, 1, 0)], rays=[(0, 0, 1)])
            sage: y = x.projection()
            sage: del y.arrows
            sage: y.arrows = Sequence([])
            sage: y._init_lines_arrows(x)
            sage: y.arrows
            [[0, 1], [0, 2], [0, 3], [0, 4], [0, 5]]
        """
        obj = polyhedron.Vrepresentation()
        adj_matrix = polyhedron.vertex_adjacency_matrix()
        for i in range(len(obj)):
            if not obj[i].is_vertex():
                if any(adj_matrix[i, j] != 0
                       for j in range(len(obj))):
                    continue
                # obj[i] is ray or line
                v = polyhedron.vertices()[0].vector()
                r = obj[i].vector()
                self.arrows.append( [ self.coord_index_of(v),
                                      self.coord_index_of(v + r) ] )
                if obj[i].is_line():
                    self.arrows.append( [ self.coord_index_of(v),
                                          self.coord_index_of(v - r) ] )
            for j in range(len(obj)):
                if adj_matrix[i, j] == 0:
                    continue
                if i < j and obj[j].is_vertex():
                    l = [obj[i].vector(), obj[j].vector()]
                    self.lines.append( [ self.coord_index_of(l[0]),
                                         self.coord_index_of(l[1]) ] )
                if obj[j].is_ray():
                    l = [obj[i].vector(), obj[i].vector() + obj[j].vector()]
                    self.arrows.append( [ self.coord_index_of(l[0]),
                                          self.coord_index_of(l[1]) ] )
                if obj[j].is_line():
                    l1 = [obj[i].vector(), obj[i].vector() + obj[j].vector()]
                    l2 = [obj[i].vector(), obj[i].vector() - obj[j].vector()]
                    self.arrows.append( [ self.coord_index_of(l1[0]),
                                          self.coord_index_of(l1[1]) ] )
                    self.arrows.append( [ self.coord_index_of(l2[0]),
                                          self.coord_index_of(l2[1]) ] )

    def _init_area_2d(self, polyhedron):
        """
        Internal function: Initialize polygon area for 2d polyhedron.

        TESTS::

            sage: p = polytopes.cyclic_polytope(2,4)
            sage: proj = p.projection()
            sage: proj.polygons = Sequence([])
            sage: proj._init_area_2d(p)
            sage: proj.polygons
            [[3, 0, 1, 2]]
        """
        assert polyhedron.ambient_dim() == 2, "Requires polyhedron in 2d"
        vertices = list(polyhedron.Vrep_generator())
        vertices = cyclic_sort_vertices_2d(vertices)
        coords = []

        def adjacent_vertices(i):
            n = len(vertices)
            if vertices[(i-1) % n].is_vertex():
                yield vertices[(i-1) % n]
            if vertices[(i+1) % n].is_vertex():
                yield vertices[(i+1) % n]

        for i in range(len(vertices)):
            v = vertices[i]
            if v.is_vertex():
                coords.append(v())
            if v.is_ray():
                for a in adjacent_vertices(i):
                    coords.append(a() + v())

        if polyhedron.n_lines() == 0:
            self.polygons.append(self.coord_indices_of(coords))
            return

        polygons = []

        if polyhedron.n_lines() == 1:
            a_line = next(polyhedron.line_generator())
            for shift in [a_line(), -a_line()]:
                for i in range(len(coords)):
                    polygons.append([coords[i - 1], coords[i],
                                     coords[i] + shift, coords[i - 1] + shift])

        if polyhedron.n_lines() == 2:
            line1, line2 = polyhedron.lines()
            assert len(coords) == 1, "Can have only a single vertex!"
            v = coords[0]
            l1 = line1()
            l2 = line2()
            polygons = [[v-l1-l2, v+l1-l2, v+l1+l2, v-l1+l2]]

        polygons = [self.coord_indices_of(p) for p in polygons]
        self.polygons.extend(polygons)

    def _init_solid_3d(self, polyhedron):
        """
        Internal function: Initialize facet polygons for 3d polyhedron.

        TESTS::

            sage: p = polytopes.cyclic_polytope(3,4)
            sage: proj = p.projection()
            sage: proj.polygons = Sequence([])
            sage: proj._init_solid_3d(p)
            sage: proj.polygons
            [[1, 0, 2], [3, 0, 1], [2, 0, 3], [3, 1, 2]]

            sage: x = Polyhedron(rays=[(-1, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)])
            sage: y = x.projection()
            sage: y.polygons
            [[5, 2, 1, 6], [2, 7, 8, 1]]

            sage: cylinder = Polyhedron(vertices=[(0, 0, 0), (1, 0, 0), (0, 1, 0)], lines=[(0, 0, 1)])
            sage: len(cylinder.projection().polygons)
            3
        """
        assert polyhedron.ambient_dim() == 3, "Requires polyhedron in 3d"

        if polyhedron.dim() <= 1:  # empty or 0d or 1d polyhedron => no polygon
            return None

        def defining_equation():  # corresponding to a polygon
            if polyhedron.dim() < 3:
                yield next(polyhedron.equation_generator())
            else:
                yield from polyhedron.inequality_generator()

        faces = []
        face_inequalities = []
        for facet_equation in defining_equation():
            vertices = list(facet_equation.incident())
            face_inequalities.append(facet_equation)
            vertices = cyclic_sort_vertices_2d(vertices)
            if len(vertices) >= 3:
                v0, v1, v2 = (vector(v) for v in vertices[:3])
                normal = (v2 - v0).cross_product(v1 - v0)
                if normal.dot_product(facet_equation.A()) < 0:
                    vertices.reverse()
            coords = []

            def adjacent_vertices(i):
                n = len(vertices)
                if vertices[(i-1) % n].is_vertex():
                    yield vertices[(i-1) % n]
                if vertices[(i+1) % n].is_vertex():
                    yield vertices[(i+1) % n]

            for i in range(len(vertices)):
                v = vertices[i]
                if v.is_vertex():
                    coords.append(v())
                if v.is_ray():
                    for a in adjacent_vertices(i):
                        coords.append(a() + v())

            faces.append(coords)
        self.face_inequalities = face_inequalities

        if polyhedron.n_lines() == 0:
            assert faces, "no vertices?"
            self.polygons.extend( [self.coord_indices_of(f) for f in faces] )
            return

        # now some special cases if there are lines (dim < ambient_dim)
        polygons = []

        if polyhedron.n_lines() == 1:
            assert faces, "no vertices?"
            a_line = next(polyhedron.line_generator())
            shift = a_line()
            for coords in faces:
                assert len(coords) == 2, "There must be two points."
                polygons.append([coords[0] - shift, coords[1] - shift,
                                 coords[1] + shift, coords[0] + shift])

        if polyhedron.n_lines() == 2:
            line1, line2 = polyhedron.line_generator()
            l1 = line1()
            l2 = line2()
            for v in polyhedron.vertex_generator():
                polygons.append([v()-l1-l2, v()+l1-l2, v()+l1+l2, v()-l1+l2])

        self.polygons.extend([self.coord_indices_of(p) for p in polygons])

    def render_points_1d(self, **kwds):
        """
        Return the points of a polyhedron in 1d.

        INPUT:

        - ``**kwds`` -- options passed through to
          :func:`~sage.plot.point.point2d`

        OUTPUT: a 2-d graphics object

        EXAMPLES::

            sage: cube1 = polytopes.hypercube(1)
            sage: proj = cube1.projection()
            sage: points = proj.render_points_1d()                                      # needs sage.plot
            sage: points._objects                                                       # needs sage.plot
            [Point set defined by 2 point(s)]
        """
        return point2d([c + [0] for c in self.coordinates_of(self.points)], **kwds)

    def render_line_1d(self, **kwds):
        """
        Return the line of a polyhedron in 1d.

        INPUT:

        - ``**kwds`` -- options passed through to
          :func:`~sage.plot.line.line2d`

        OUTPUT: a 2-d graphics object

        EXAMPLES::

            sage: outline = polytopes.hypercube(1).projection().render_line_1d()        # needs sage.plot
            sage: outline._objects[0]                                                   # needs sage.plot
            Line defined by 2 points
        """
        if len(self.lines) == 0:
            return Graphics()
        elif len(self.lines) == 1:
            line = self.coordinates_of(self.lines[0])
            return line2d([line[0] + [0], line[1] + [0]], **kwds)
        else:
            assert False   # unreachable

    def render_points_2d(self, **kwds):
        """
        Return the points of a polyhedron in 2d.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: hex = polytopes.regular_polygon(6)
            sage: proj = hex.projection()
            sage: hex_points = proj.render_points_2d()                                  # needs sage.plot
            sage: hex_points._objects                                                   # needs sage.plot
            [Point set defined by 6 point(s)]
        """
        return point2d(self.coordinates_of(self.points), **kwds)

    def render_outline_2d(self, **kwds):
        """
        Return the outline (edges) of a polyhedron in 2d.

        EXAMPLES::

            sage: penta = polytopes.regular_polygon(5)                                  # needs sage.rings.number_field
            sage: outline = penta.projection().render_outline_2d()                      # needs sage.plot sage.rings.number_field
            sage: outline._objects[0]                                                   # needs sage.plot sage.rings.number_field
            Line defined by 2 points
        """
        wireframe = []
        for l in self.lines:
            l_coords = self.coordinates_of(l)
            wireframe.append(line2d(l_coords, **kwds))
        for a in self.arrows:
            a_coords = self.coordinates_of(a)
            wireframe.append(arrow(a_coords[0], a_coords[1], **kwds))
        return sum(wireframe)

    def render_fill_2d(self, **kwds):
        """
        Return the filled interior (a polygon) of a polyhedron in 2d.

        EXAMPLES::

            sage: cps = [i^3 for i in srange(-2, 2, 1/5)]
            sage: p = Polyhedron(vertices=[[(t^2-1)/(t^2+1), 2*t/(t^2+1)] for t in cps])
            sage: proj = p.projection()
            sage: filled_poly = proj.render_fill_2d()                                   # needs sage.plot
            sage: filled_poly.axes_width()                                              # needs sage.plot
            0.8
        """
        poly = [polygon2d(self.coordinates_of(p), **kwds)
                for p in self.polygons]
        return sum(poly)

    def render_vertices_3d(self, **kwds):
        """
        Return the 3d rendering of the vertices.

        EXAMPLES::

            sage: p = polytopes.cross_polytope(3)
            sage: proj = p.projection()
            sage: verts = proj.render_vertices_3d()                                     # needs sage.plot
            sage: verts.bounding_box()                                                  # needs sage.plot
            ((-1.0, -1.0, -1.0), (1.0, 1.0, 1.0))
        """
        return point3d(self.coordinates_of(self.points), **kwds)

    def render_wireframe_3d(self, **kwds):
        r"""
        Return the 3d wireframe rendering.

        EXAMPLES::

            sage: cube = polytopes.hypercube(3)
            sage: cube_proj = cube.projection()
            sage: wire = cube_proj.render_wireframe_3d()                                # needs sage.plot
            sage: print(wire.tachyon().split('\n')[77])  # for testing                  # needs sage.plot
            FCylinder base 1.0 1.0 -1.0 apex -1.0 1.0 -1.0 rad 0.005 texture...
        """
        wireframe = []
        for l in self.lines:
            l_coords = self.coordinates_of(l)
            wireframe.append(line3d(l_coords, **kwds))
        for a in self.arrows:
            a_coords = self.coordinates_of(a)
            wireframe.append(arrow3d(a_coords[0], a_coords[1], **kwds))
        return sum(wireframe)

    def render_solid_3d(self, **kwds):
        """
        Return solid 3d rendering of a 3d polytope.

        EXAMPLES::

            sage: p = polytopes.hypercube(3).projection()
            sage: p_solid = p.render_solid_3d(opacity=.7)                               # needs sage.plot
            sage: type(p_solid)                                                         # needs sage.plot
            <class 'sage.plot.plot3d.index_face_set.IndexFaceSet'>
        """
        polys = self.polygons
        n = len(polys)
        N = max([-1] + [i for p in polys for i in p]) + 1
        coords = self.coordinates_of(range(N))
        col = kwds.pop('color', (0, 0, 1))
        if col == 'rainbow':
            t_list = [Texture(rainbow(n, 'rgbtuple')[i]) for i in range(n)]
            return polygons3d(polys, coords, texture_list=t_list, **kwds)
        else:
            return polygons3d(polys, coords, color=col, **kwds)

    def render_0d(self, point_opts=None, line_opts=None, polygon_opts=None):
        """
        Return 0d rendering of the projection of a polyhedron into
        2-dimensional ambient space.

        INPUT:

        See
        :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.plot`.

        OUTPUT: a 2-d graphics object

        EXAMPLES::

            sage: print(Polyhedron([]).projection().render_0d().description())          # needs sage.plot
            <BLANKLINE>
            sage: P = Polyhedron(ieqs=[(1,)])
            sage: print(P.projection().render_0d().description())                       # needs sage.plot
            Point set defined by 1 point(s):    [(0.0, 0.0)]
        """
        if point_opts is None:
            point_opts = {}
        if isinstance(point_opts, dict):
            point_opts.setdefault('zorder', 2)
            point_opts.setdefault('pointsize', 10)
        if self.points:
            return point2d([0, 0], **point_opts)
        else:
            return Graphics()

    def render_1d(self, point_opts=None, line_opts=None, polygon_opts=None):
        """
        Return 1d rendering of the projection of a polyhedron into
        2-dimensional ambient space.

        INPUT:

        See
        :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.plot`.

        OUTPUT: a 2-d graphics object

        EXAMPLES::

            sage: Polyhedron([(0,), (1,)]).projection().render_1d()                     # needs sage.plot
            Graphics object consisting of 2 graphics primitives
        """
        plt = Graphics()
        if point_opts is None:
            point_opts = {}
        if line_opts is None:
            line_opts = {}
        if isinstance(point_opts, dict):
            point_opts.setdefault('zorder', 2)
            point_opts.setdefault('pointsize', 10)
            plt += self.render_points_1d(**point_opts)
        if isinstance(line_opts, dict):
            line_opts.setdefault('zorder', 1)
            plt += self.render_line_1d(**line_opts)
        return plt

    def render_2d(self, point_opts=None, line_opts=None, polygon_opts=None):
        """
        Return 2d rendering of the projection of a polyhedron into
        2-dimensional ambient space.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[[1,1]], rays=[[1,1]])
            sage: q1 = p1.projection()
            sage: p2 = Polyhedron(vertices=[[1,0], [0,1], [0,0]])
            sage: q2 = p2.projection()
            sage: p3 = Polyhedron(vertices=[[1,2]])
            sage: q3 = p3.projection()
            sage: p4 = Polyhedron(vertices=[[2,0]], rays=[[1,-1]], lines=[[1,1]])
            sage: q4 = p4.projection()
            sage: q1.plot() + q2.plot() + q3.plot() + q4.plot()                         # needs sage.plot
            Graphics object consisting of 18 graphics primitives
        """
        plt = Graphics()
        if point_opts is None:
            point_opts = {}
        if line_opts is None:
            line_opts = {}
        if polygon_opts is None:
            polygon_opts = {}
        if isinstance(point_opts, dict):
            point_opts.setdefault('zorder', 2)
            point_opts.setdefault('pointsize', 10)
            plt += self.render_points_2d(**point_opts)
        if isinstance(line_opts, dict):
            line_opts.setdefault('zorder', 1)
            plt += self.render_outline_2d(**line_opts)
        if isinstance(polygon_opts, dict):
            polygon_opts.setdefault('zorder', 0)
            plt += self.render_fill_2d(**polygon_opts)
        return plt

    def render_3d(self, point_opts=None, line_opts=None, polygon_opts=None):
        """
        Return 3d rendering of a polyhedron projected into
        3-dimensional ambient space.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[[1,1,1]], rays=[[1,1,1]])
            sage: p2 = Polyhedron(vertices=[[2,0,0], [0,2,0], [0,0,2]])
            sage: p3 = Polyhedron(vertices=[[1,0,0], [0,1,0], [0,0,1]],
            ....:                 rays=[[-1,-1,-1]])
            sage: (p1.projection().plot() + p2.projection().plot()                      # needs sage.plot
            ....:   + p3.projection().plot())
            Graphics3d Object

        It correctly handles various degenerate cases::

            sage: # needs sage.plot
            sage: Polyhedron(lines=[[1,0,0], [0,1,0], [0,0,1]]).plot()  # whole space
            Graphics3d Object
            sage: Polyhedron(vertices=[[1,1,1]], rays=[[1,0,0]],
            ....:            lines=[[0,1,0], [0,0,1]]).plot()           # half space
            Graphics3d Object
            sage: Polyhedron(lines=[[0,1,0], [0,0,1]],
            ....:            vertices=[[1,1,1]]).plot()      # R^2 in R^3
            Graphics3d Object
            sage: Polyhedron(rays=[[0,1,0], [0,0,1]],        # quadrant wedge in R^2
            ....:            lines=[[1,0,0]]).plot()
            Graphics3d Object
            sage: Polyhedron(rays=[[0,1,0]],                 # upper half plane in R^3
            ....:            lines=[[1,0,0]]).plot()
            Graphics3d Object
            sage: Polyhedron(lines=[[1,0,0]]).plot()         # R^1 in R^2
            Graphics3d Object
            sage: Polyhedron(rays=[[0,1,0]]).plot()          # Half-line in R^3
            Graphics3d Object
            sage: Polyhedron(vertices=[[1,1,1]]).plot()      # point in R^3
            Graphics3d Object

        The origin is not included, if it is not in the polyhedron (:issue:`23555`)::

            sage: Q = Polyhedron([[100],[101]])
            sage: P = Q*Q*Q; P
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 8 vertices
            sage: p = P.plot()                                                          # needs sage.plot
            sage: p.bounding_box()                                                      # needs sage.plot
            ((100.0, 100.0, 100.0), (101.0, 101.0, 101.0))

        Plot 3d polytope with rainbow colors::

            sage: polytopes.hypercube(3).plot(polygon='rainbow', alpha=0.4)             # needs sage.plot
            Graphics3d Object
        """
        pplt = None
        lplt = None
        pgplt = None
        if point_opts is None:
            point_opts = {}
        if line_opts is None:
            line_opts = {}
        if polygon_opts is None:
            polygon_opts = {}
        if isinstance(point_opts, dict):
            point_opts.setdefault('size', 10)
            pplt = self.render_vertices_3d(**point_opts)
        if isinstance(line_opts, dict):
            line_opts.setdefault('width', 1) # controls the width of arrow3d for a ray
            line_opts.setdefault('thickness', 1) # controls the thickness of line3d
            lplt = self.render_wireframe_3d(**line_opts)
        if isinstance(polygon_opts, dict):
            if 'threejs_flat_shading' not in polygon_opts:
                polygon_opts['threejs_flat_shading'] = True
            pgplt = self.render_solid_3d(**polygon_opts)
        # zorder is not available
        return sum(_ for _ in [pplt, lplt, pgplt] if _ is not None)

    def tikz(self, view=[0, 0, 1], angle=0, scale=1,
             edge_color='blue!95!black', facet_color='blue!95!black',
             opacity=0.8, vertex_color='green', axis=False,
             output_type=None):
        r"""
        Return a tikz picture of ``self`` as a string or as a
        :class:`~sage.misc.latex_standalone.TikzPicture`
        according to a projection ``view`` and an angle ``angle``
        obtained via the threejs viewer.

        INPUT:

        - ``view`` -- list (default: [0,0,1]) representing the rotation axis (see note below)
        - ``angle`` -- integer (default: 0); angle of rotation in degree from 0 to 360 (see note
          below)
        - ``scale`` -- integer (default: 1); the scaling of the tikz picture
        - ``edge_color`` -- string (default: ``'blue!95!black'``); representing colors which tikz
          recognizes
        - ``facet_color`` -- string (default: ``'blue!95!black'``); representing colors which tikz
          recognizes
        - ``vertex_color`` -- string (default: ``'green'``); representing colors which tikz
          recognizes
        - ``opacity`` -- real number (default: 0.8) between 0 and 1 giving the opacity of
          the front facets
        - ``axis`` -- boolean (default: ``False``); draw the axes at the origin or not
        - ``output_type`` -- string (default: ``None``); valid values
          are ``None`` (deprecated), ``'LatexExpr'`` and ``'TikzPicture'``,
          whether to return a :class:`LatexExpr` object (which inherits from Python
          :class:`str`) or a :class:`TikzPicture` object from module
          :mod:`sage.misc.latex_standalone`

        OUTPUT: :class:`LatexExpr` object or :class:`TikzPicture` object

        .. NOTE::

            The inputs ``view`` and ``angle`` can be obtained by visualizing it
            using ``.show(aspect_ratio=1)``. This will open an interactive view
            in your default browser, where you can rotate the polytope. Once
            the desired view angle is found, click on the information icon in
            the lower right-hand corner and select *Get Viewpoint*. This will
            copy a string of the form '[x,y,z],angle' to your local clipboard.
            Go back to Sage and type ``Img = P.projection().tikz([x,y,z],angle)``.

            The inputs ``view`` and ``angle`` can also be obtained from the
            viewer Jmol::

                1) Right click on the image
                2) Select ``Console``
                3) Select the tab ``State``
                4) Scroll to the line ``moveto``

            It reads something like::

                moveto 0.0 {x y z angle} Scale

            The ``view`` is then [x,y,z] and ``angle`` is angle.
            The following number is the scale.

            Jmol performs a rotation of ``angle`` degrees along the
            vector [x,y,z] and show the result from the z-axis.

        EXAMPLES::

            sage: # needs sage.plot sage.rings.number_field
            sage: P1 = polytopes.small_rhombicuboctahedron()
            sage: Image1 = P1.projection().tikz([1,3,5], 175, scale=4,
            ....:                               output_type='TikzPicture')
            sage: type(Image1)
            <class 'sage.misc.latex_standalone.TikzPicture'>
            sage: Image1
            \documentclass[tikz]{standalone}
            \begin{document}
            \begin{tikzpicture}%
                    [x={(-0.939161cm, 0.244762cm)},
                    y={(0.097442cm, -0.482887cm)},
                    z={(0.329367cm, 0.840780cm)},
                    scale=4.000000,
            ...
            Use print to see the full content.
            ...
            \node[vertex] at (-2.41421, 1.00000, -1.00000)     {};
            \node[vertex] at (-2.41421, -1.00000, 1.00000)     {};
            %%
            %%
            \end{tikzpicture}
            \end{document}
            sage: _ = Image1.tex('polytope-tikz1.tex')          # not tested
            sage: _ = Image1.png('polytope-tikz1.png')          # not tested
            sage: _ = Image1.pdf('polytope-tikz1.pdf')          # not tested
            sage: _ = Image1.svg('polytope-tikz1.svg')          # not tested

        A second example::

            sage: P2 = Polyhedron(vertices=[[1, 1], [1, 2], [2, 1]])
            sage: Image2 = P2.projection().tikz(scale=3, edge_color='blue!95!black',
            ....:                               facet_color='orange!95!black', opacity=0.4,
            ....:                               vertex_color='yellow', axis=True,
            ....:                               output_type='TikzPicture')
            sage: Image2
            \documentclass[tikz]{standalone}
            \begin{document}
            \begin{tikzpicture}%
                    [scale=3.000000,
                    back/.style={loosely dotted, thin},
                    edge/.style={color=blue!95!black, thick},
                    facet/.style={fill=orange!95!black,fill opacity=0.400000},
            ...
            Use print to see the full content.
            ...
            \node[vertex] at (1.00000, 2.00000)     {};
            \node[vertex] at (2.00000, 1.00000)     {};
            %%
            %%
            \end{tikzpicture}
            \end{document}

        The second example using a LatexExpr as output type::

            sage: # needs sage.plot
            sage: Image2 = P2.projection().tikz(scale=3, edge_color='blue!95!black',
            ....:                               facet_color='orange!95!black', opacity=0.4,
            ....:                               vertex_color='yellow', axis=True,
            ....:                               output_type='LatexExpr')
            sage: type(Image2)
            <class 'sage.misc.latex.LatexExpr'>
            sage: print('\n'.join(Image2.splitlines()[:4]))
            \begin{tikzpicture}%
                [scale=3.000000,
                back/.style={loosely dotted, thin},
                edge/.style={color=blue!95!black, thick},
            sage: with open('polytope-tikz2.tex', 'w') as f:    # not tested
            ....:     _ = f.write(Image2)

        A third example::

            sage: # needs sage.plot
            sage: P3 = Polyhedron(vertices=[[-1, -1, 2], [-1, 2, -1], [2, -1, -1]]); P3
            A 2-dimensional polyhedron in ZZ^3 defined as the convex hull of 3 vertices
            sage: Image3 = P3.projection().tikz([0.5, -1, -0.1], 55, scale=3,
            ....:                               edge_color='blue!95!black',
            ....:                               facet_color='orange!95!black', opacity=0.7,
            ....:                               vertex_color='yellow', axis=True,
            ....:                               output_type='TikzPicture')
            sage: Image3
            \documentclass[tikz]{standalone}
            \begin{document}
            \begin{tikzpicture}%
                    [x={(0.658184cm, -0.242192cm)},
                    y={(-0.096240cm, 0.912008cm)},
                    z={(-0.746680cm, -0.331036cm)},
                    scale=3.000000,
            ...
            Use print to see the full content.
            ...
            \node[vertex] at (-1.00000, 2.00000, -1.00000)     {};
            \node[vertex] at (2.00000, -1.00000, -1.00000)     {};
            %%
            %%
            \end{tikzpicture}
            \end{document}
            sage: _ = Image3.tex('polytope-tikz3.tex')          # not tested
            sage: _ = Image3.png('polytope-tikz3.png')          # not tested
            sage: _ = Image3.pdf('polytope-tikz3.pdf')          # not tested
            sage: _ = Image3.svg('polytope-tikz3.svg')          # not tested

        A fourth example::

            sage: P = Polyhedron(vertices=[[1,1,0,0], [1,2,0,0],
            ....:                          [2,1,0,0], [0,0,1,0], [0,0,0,1]]); P
            A 4-dimensional polyhedron in ZZ^4 defined as the convex hull of 5 vertices
            sage: P.projection().tikz(output_type='TikzPicture')
            Traceback (most recent call last):
            ...
            NotImplementedError: The polytope has to live in 2 or 3 dimensions.

        .. TODO::

            Make it possible to draw Schlegel diagram for 4-polytopes. ::

                sage: P = Polyhedron(vertices=[[1,1,0,0], [1,2,0,0],
                ....:                          [2,1,0,0], [0,0,1,0], [0,0,0,1]]); P
                A 4-dimensional polyhedron in ZZ^4 defined as the convex hull of 5 vertices
                sage: P.projection().tikz(output_type='TikzPicture')
                Traceback (most recent call last):
                ...
                NotImplementedError: The polytope has to live in 2 or 3 dimensions.

            Make it possible to draw 3-polytopes living in higher dimension.
        """
        if self.polyhedron_ambient_dim > 3 or self.polyhedron_ambient_dim < 2:
            raise NotImplementedError("The polytope has to live in 2 or 3 dimensions.")
        elif self.polyhedron_dim < 2 or self.polyhedron_dim > 3:
            raise NotImplementedError("The polytope has to be 2 or 3-dimensional.")
        elif self.polyhedron_ambient_dim == 2:  # self is a polygon in 2-space
            tikz_string = self._tikz_2d(scale, edge_color, facet_color, opacity,
                                 vertex_color, axis)
        elif self.polyhedron_dim == 2:  # self is a polygon in 3-space
            tikz_string = self._tikz_2d_in_3d(view, angle, scale, edge_color,
                                       facet_color, opacity, vertex_color, axis)
        else:  # self is a 3-polytope in 3-space
            tikz_string = self._tikz_3d_in_3d(view, angle, scale, edge_color,
                                       facet_color, opacity, vertex_color, axis)

        # set default value
        if output_type is None:
            from sage.misc.superseded import deprecation
            msg = ("The default type of the returned object will soon be "
                   "changed from `sage.misc.latex.LatexExpr` to "
                   "`sage.misc.latex_standalone.TikzPicture`.  Please "
                   "update your code to specify the desired output type as "
                   "`.tikz(output_type='LatexExpr')` to keep the old "
                   "behavior or `.tikz(output_type='TikzPicture')` to use "
                   "the future default behavior.")
            deprecation(33002, msg)
            output_type = 'LatexExpr'

        # return
        if output_type == 'LatexExpr':
            return tikz_string
        elif output_type == 'TikzPicture':
            from sage.misc.latex_standalone import TikzPicture
            return TikzPicture(tikz_string, standalone_config=None,
                    usepackage=None, usetikzlibrary=None, macros=None,
                    use_sage_preamble=False)
        else:
            raise ValueError("output_type (='{}') must be 'LatexExpr' or"
                    " 'TikzPicture'".format(output_type))

    def _tikz_2d(self, scale, edge_color, facet_color, opacity, vertex_color, axis):
        r"""
        Return a string ``tikz_pic`` consisting of a tikz picture of
        ``self``, which is assumed to be a polygon on the plane.

        INPUT:

        - ``scale`` -- integer specifying the scaling of the tikz picture
        - ``edge_color`` -- string representing colors which tikz
          recognizes
        - ``facet_color`` -- string representing colors which tikz
          recognizes
        - ``vertex_color`` -- string representing colors which tikz
          recognizes
        - ``opacity`` -- real number between 0 and 1 giving the opacity of
          the front facets
        - ``axis`` -- boolean (default: ``False``); draw the axes at the origin or not

        OUTPUT:

        :class:`LatexExpr` -- containing the TikZ picture.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[[1, 1], [1, 2], [2, 1]])
            sage: Image = P.projection()._tikz_2d(scale=3, edge_color='black',
            ....:                                 facet_color='orange', opacity=0.75,
            ....:                                 vertex_color='yellow', axis=True)
            sage: type(Image)
            <class 'sage.misc.latex.LatexExpr'>
            sage: print('\n'.join(Image.splitlines()[:4]))
            \begin{tikzpicture}%
                [scale=3.000000,
                back/.style={loosely dotted, thin},
                edge/.style={color=black, thick},
            sage: with open('polytope-tikz2.tex', 'w') as f:  # not tested
            ....:     _ = f.write(Image)

        Scientific notation is not used in the output (:issue:`16519`)::

            sage: P = Polyhedron([[2*10^-10,0], [0,1], [1,0]], base_ring=QQ)
            sage: tikz = P.projection().tikz(output_type='TikzPicture')
            sage: 'e-10' in tikz.content()
            False

        .. NOTE::

            The ``facet_color`` is the filing color of the polytope (polygon).
        """

        # Creates the nodes, coordinate and tag for every vertex of the polytope.
        # The tag is used to draw the front facets later on.

        dict_drawing = {}
        edges = ''
        for vert in self.points:
            v = self.coords[vert]
            v_vect = str(['%.5f' % i for i in v]).replace('\'', '')
            v_vect = v_vect.replace('[', '(')
            v_vect = v_vect.replace(']', ')')
            tag = '%s' % v_vect
            node = "\\node[%s] at %s     {};\n" % ('vertex', tag)
            coord = '\\coordinate %s at %s;\n' % (tag, tag)
            dict_drawing[vert] = node, coord, tag

        for index1, index2 in self.lines:
            # v1 = self.coords[index1]
            # v2 = self.coords[index2]
            edges += "\\draw[%s] %s -- %s;\n" % ('edge',
                                                 dict_drawing[index1][2],
                                                 dict_drawing[index2][2])

        # Start to write the output
        tikz_pic = ''
        tikz_pic += '\\begin{tikzpicture}%\n'
        tikz_pic += '\t[scale=%f,\n' % scale
        tikz_pic += '\tback/.style={loosely dotted, thin},\n'
        tikz_pic += '\tedge/.style={color=%s, thick},\n' % edge_color
        tikz_pic += '\tfacet/.style={fill=%s,fill opacity=%f},\n' % (facet_color, opacity)
        tikz_pic += '\tvertex/.style={inner sep=1pt,circle,draw=%s!25!black,' % vertex_color
        tikz_pic += 'fill=%s!75!black,thick}]\n%%\n%%\n' % vertex_color

        # Gives the reproduction information
        from sage.env import SAGE_VERSION
        tikz_pic += "%% This TikZ-picture was produced with Sagemath version {}\n".format(SAGE_VERSION)
        tikz_pic += "%% with the command: ._tikz_2d and parameters:\n"
        tikz_pic += "%% scale = {}\n".format(scale)
        tikz_pic += "%% edge_color = {}\n".format(edge_color)
        tikz_pic += "%% facet_color = {}\n".format(facet_color)
        tikz_pic += "%% opacity = {}\n".format(opacity)
        tikz_pic += "%% vertex_color = {}\n".format(vertex_color)
        tikz_pic += "%% axis = {}\n%%\n".format(axis)

        # Draws the axes if True
        if axis:
            tikz_pic += '%% Drawing the axes\n'
            tikz_pic += '\\draw[color=black,thick,->] (0,0,0) -- (1,0,0) node[anchor=north east]{$x$};\n'
            tikz_pic += '\\draw[color=black,thick,->] (0,0,0) -- (0,1,0) node[anchor=north west]{$y$};\n\n'

        # Create the coordinate of the vertices:
        tikz_pic += '%% Coordinate of the vertices:\n%%\n'
        for v in dict_drawing:
            tikz_pic += dict_drawing[v][1]

        # Draw the interior by going in a cycle
        vertices = list(self.parent_polyhedron.Vrep_generator())
        tikz_pic += '%%\n%%\n%% Drawing the interior\n%%\n'
        cyclic_vert = cyclic_sort_vertices_2d(list(self.parent_polyhedron.Vrep_generator()))
        cyclic_indices = [vertices.index(v) for v in cyclic_vert]
        tikz_pic += '\\fill[facet] '
        for v in cyclic_indices:
            if v in dict_drawing:
                tikz_pic += '%s -- ' % dict_drawing[v][2]
        tikz_pic += 'cycle {};\n'

        # Draw the edges
        tikz_pic += '%%\n%%\n%% Drawing edges\n%%\n'
        tikz_pic += edges

        # Finally, the vertices in front are drawn on top of everything.
        tikz_pic += '%%\n%%\n%% Drawing the vertices\n%%\n'
        for v in dict_drawing:
            tikz_pic += dict_drawing[v][0]
        tikz_pic += '%%\n%%\n\\end{tikzpicture}'

        return LatexExpr(tikz_pic)

    def _tikz_2d_in_3d(self, view, angle, scale, edge_color, facet_color,
                       opacity, vertex_color, axis):
        r"""
        Return a string ``tikz_pic`` consisting of a tikz picture of ``self``
        according to a projection ``view`` and an angle ``angle``
        obtained via Jmol through the current state property. ``self`` is
        assumed to be a polygon in 3-space.

        INPUT:

        - ``view`` -- list (default: [0,0,1]) representing the rotation axis
        - ``angle`` -- integer angle of rotation in degree from 0 to 360
        - ``scale`` -- integer specifying the scaling of the tikz picture
        - ``edge_color`` -- string representing colors which tikz
          recognizes
        - ``facet_color`` -- string representing colors which tikz
          recognizes
        - ``vertex_color`` -- string representing colors which tikz
          recognizes
        - ``opacity`` -- real number between 0 and 1 giving the opacity of
          the front facets
        - ``axis`` -- boolean draw the axes at the origin or not

        OUTPUT: :class:`LatexExpr` -- containing the TikZ picture

        EXAMPLES::

            sage: # needs sage.plot
            sage: P = Polyhedron(vertices=[[-1, -1, 2], [-1, 2, -1], [2, -1, -1]]); P
            A 2-dimensional polyhedron in ZZ^3 defined as the convex hull of 3 vertices
            sage: Image = P.projection()._tikz_2d_in_3d(view=[0.5, -1, -0.5], angle=55, scale=3,
            ....:                                       edge_color='blue!95!black', facet_color='orange',
            ....:                                       opacity=0.5, vertex_color='yellow', axis=True)
            sage: print('\n'.join(Image.splitlines()[:4]))
            \begin{tikzpicture}%
                [x={(0.644647cm, -0.476559cm)},
                y={(0.192276cm, 0.857859cm)},
                z={(-0.739905cm, -0.192276cm)},
            sage: with open('polytope-tikz3.tex', 'w') as f:    # not tested
            ....:     _ = f.write(Image)

        ::

            sage: p = Polyhedron(vertices=[[1, 0, 0], [0, 1, 0], [0, 0, 1]])
            sage: proj = p.projection()
            sage: Img = proj.tikz([1, 1, 1], 130, axis=True, output_type='LatexExpr')   # needs sage.plot
            sage: print('\n'.join(Img.splitlines()[12:21]))                             # needs sage.plot
            %% with the command: ._tikz_2d_in_3d and parameters:
            %% view = [1, 1, 1]
            %% angle = 130
            %% scale = 1
            %% edge_color = blue!95!black
            %% facet_color = blue!95!black
            %% opacity = 0.8
            %% vertex_color = green
            %% axis = True

        .. NOTE::

            The ``facet_color`` is the filing color of the polytope (polygon).
        """
        from sage.rings.real_double import RDF

        view_vector = vector(RDF, view)
        rot = rotate_arbitrary(view_vector, -(angle/360)*2*pi)
        rotation_matrix = rot[:2].transpose()

        # Creates the nodes, coordinate and tag for every vertex of the polytope.
        # The tag is used to draw the front facets later on.
        dict_drawing = {}
        edges = ''
        for vert in self.points:
            v = self.coords[vert]
            v_vect = str(['%.5f' % i for i in v]).replace('\'', '')
            v_vect = v_vect.replace('[', '(')
            v_vect = v_vect.replace(']', ')')
            tag = '%s' % v_vect
            node = "\\node[%s] at %s     {};\n" % ('vertex', tag)
            coord = '\\coordinate %s at %s;\n' % (tag, tag)
            dict_drawing[vert] = node, coord, tag

        for index1, index2 in self.lines:
            # v1 = self.coords[index1]
            # v2 = self.coords[index2]
            edges += "\\draw[%s] %s -- %s;\n" % ('edge',
                                                 dict_drawing[index1][2],
                                                 dict_drawing[index2][2])

        # Start to write the output
        tikz_pic = ''
        tikz_pic += '\\begin{tikzpicture}%\n'
        tikz_pic += '\t[x={(%fcm, %fcm)},\n' % (RDF(rotation_matrix[0][0]),
                                                RDF(rotation_matrix[0][1]))
        tikz_pic += '\ty={(%fcm, %fcm)},\n' % (RDF(rotation_matrix[1][0]),
                                               RDF(rotation_matrix[1][1]))
        tikz_pic += '\tz={(%fcm, %fcm)},\n' % (RDF(rotation_matrix[2][0]),
                                               RDF(rotation_matrix[2][1]))
        tikz_pic += '\tscale=%f,\n' % scale
        tikz_pic += '\tback/.style={loosely dotted, thin},\n'
        tikz_pic += '\tedge/.style={color=%s, thick},\n' % edge_color
        tikz_pic += '\tfacet/.style={fill=%s,fill opacity=%f},\n' % (facet_color, opacity)
        tikz_pic += '\tvertex/.style={inner sep=1pt,circle,draw=%s!25!black,' % vertex_color
        tikz_pic += 'fill=%s!75!black,thick}]\n%%\n%%\n' % vertex_color

        # Gives the reproduction information
        from sage.env import SAGE_VERSION
        tikz_pic += "%% This TikZ-picture was produced with Sagemath version {}\n".format(SAGE_VERSION)
        tikz_pic += "%% with the command: ._tikz_2d_in_3d and parameters:\n"
        tikz_pic += "%% view = {}\n".format(view)
        tikz_pic += "%% angle = {}\n".format(angle)
        tikz_pic += "%% scale = {}\n".format(scale)
        tikz_pic += "%% edge_color = {}\n".format(edge_color)
        tikz_pic += "%% facet_color = {}\n".format(facet_color)
        tikz_pic += "%% opacity = {}\n".format(opacity)
        tikz_pic += "%% vertex_color = {}\n".format(vertex_color)
        tikz_pic += "%% axis = {}\n%%\n".format(axis)

        # Draws the axes if True
        if axis:
            tikz_pic += '%% Drawing the axes\n'
            tikz_pic += '\\draw[color=black,thick,->] (0,0,0) -- (1,0,0) node[anchor=north east]{$x$};\n'
            tikz_pic += '\\draw[color=black,thick,->] (0,0,0) -- (0,1,0) node[anchor=north west]{$y$};\n'
            tikz_pic += '\\draw[color=black,thick,->] (0,0,0) -- (0,0,1) node[anchor=south]{$z$};\n'

        # Create the coordinate of the vertices:
        tikz_pic += '%% Coordinate of the vertices:\n%%\n'
        for v in dict_drawing:
            tikz_pic += dict_drawing[v][1]

        # Draw the interior by going in a cycle
        vertices = list(self.parent_polyhedron.Vrep_generator())
        tikz_pic += '%%\n%%\n%% Drawing the interior\n%%\n'
        cyclic_vert = cyclic_sort_vertices_2d(list(self.parent_polyhedron.Vrep_generator()))
        cyclic_indices = [vertices.index(v) for v in cyclic_vert]
        tikz_pic += '\\fill[facet] '
        for v in cyclic_indices:
            if v in dict_drawing:
                tikz_pic += '%s -- ' % dict_drawing[v][2]
        tikz_pic += 'cycle {};\n'

        # Draw the edges in the front
        tikz_pic += '%%\n%%\n%% Drawing edges\n%%\n'
        tikz_pic += edges

        # Finally, the vertices in front are drawn on top of everything.
        tikz_pic += '%%\n%%\n%% Drawing the vertices\n%%\n'
        for v in dict_drawing:
            tikz_pic += dict_drawing[v][0]
        tikz_pic += '%%\n%%\n\\end{tikzpicture}'

        return LatexExpr(tikz_pic)

    def _tikz_3d_in_3d(self, view, angle, scale, edge_color,
                       facet_color, opacity, vertex_color, axis):
        r"""
        Return a string ``tikz_pic`` consisting of a tikz picture of ``self``
        according to a projection ``view`` and an angle ``angle``
        obtained via Jmol through the current state property. ``self`` is
        assumed to be a 3-polytope in 3-space.

        INPUT:

        - ``view`` -- list (default: [0,0,1]) representing the rotation axis
        - ``angle`` -- integer angle of rotation in degree from 0 to 360
        - ``scale`` -- integer specifying the scaling of the tikz picture
        - ``edge_color`` -- string representing colors which tikz
          recognizes
        - ``facet_color`` -- string representing colors which tikz
          recognizes
        - ``vertex_color`` -- string representing colors which tikz
          recognizes
        - ``opacity`` -- real number between 0 and 1 giving the opacity of
          the front facets
        - ``axis`` -- boolean draw the axes at the origin or not

        OUTPUT: :class:`LatexExpr` -- containing the TikZ picture.

        EXAMPLES::

            sage: # needs sage.plot sage.rings.number_field
            sage: P = polytopes.small_rhombicuboctahedron()
            sage: Image = P.projection()._tikz_3d_in_3d([3, 7, 5], 100, scale=3,
            ....:                                       edge_color='blue', facet_color='orange',
            ....:                                       opacity=0.5, vertex_color='green', axis=True)
            sage: type(Image)
            <class 'sage.misc.latex.LatexExpr'>
            sage: print('\n'.join(Image.splitlines()[:4]))
            \begin{tikzpicture}%
                [x={(-0.046385cm, 0.837431cm)},
                y={(-0.243536cm, 0.519228cm)},
                z={(0.968782cm, 0.170622cm)},
            sage: with open('polytope-tikz1.tex', 'w') as f:    # not tested
            ....:     _ = f.write(Image)

        ::

            sage: # needs sage.plot
            sage: Associahedron = Polyhedron(vertices=[[1, 0, 1], [1, 0, 0], [1, 1, 0],
            ....:                                      [0, 0, -1], [0, 1, 0], [-1, 0, 0],
            ....:                                      [0, 1, 1], [0, 0, 1], [0, -1, 0]]).polar()
            sage: ImageAsso = Associahedron.projection().tikz([-15, -755, -655], 116, scale=1,
            ....:                                             output_type='LatexExpr')
            sage: print('\n'.join(ImageAsso.splitlines()[12:30]))
            %% with the command: ._tikz_3d_in_3d and parameters:
            %% view = [-15, -755, -655]
            %% angle = 116
            %% scale = 1
            %% edge_color = blue!95!black
            %% facet_color = blue!95!black
            %% opacity = 0.8
            %% vertex_color = green
            %% axis = False
            %%
            %% Coordinate of the vertices:
            %%
            \coordinate (0.00000, 1.00000, -1.00000) at (0.00000, 1.00000, -1.00000);
            \coordinate (1.00000, 1.00000, -1.00000) at (1.00000, 1.00000, -1.00000);
            \coordinate (1.00000, 1.00000, 1.00000) at (1.00000, 1.00000, 1.00000);
            \coordinate (1.00000, -1.00000, 1.00000) at (1.00000, -1.00000, 1.00000);
            \coordinate (1.00000, -1.00000, 0.00000) at (1.00000, -1.00000, 0.00000);
            \coordinate (1.00000, 0.00000, -1.00000) at (1.00000, 0.00000, -1.00000);
        """
        from sage.rings.real_double import RDF

        view_vector = vector(RDF, view)
        rot = rotate_arbitrary(view_vector, -(angle/360)*2*pi)
        rotation_matrix = rot[:2].transpose()
        proj_vector = (rot**(-1))*vector(RDF, [0, 0, 1])

        # First compute the back and front vertices and facets
        facets = self.face_inequalities
        front_facets = []
        back_facets = []
        for index_facet in range(len(facets)):
            A = facets[index_facet].vector()[1:]
            B = facets[index_facet].vector()[0]
            if A*(2000*proj_vector)+B < 0:
                front_facets += [index_facet]
            else:
                back_facets += [index_facet]

        vertices = list(self.parent_polyhedron.Vrep_generator())
        front_vertices = []
        for index_facet in front_facets:
            A = facets[index_facet].vector()[1:]
            B = facets[index_facet].vector()[0]
            for v in self.points:
                if A*self.coords[v]+B < 0.0005 and v not in front_vertices:
                    front_vertices += [v]

        back_vertices = []
        for index_facet in back_facets:
            A = facets[index_facet].vector()[1:]
            B = facets[index_facet].vector()[0]
            for v in self.points:
                if A * self.coords[v] + B < 0.0005 and v not in back_vertices:
                    back_vertices += [v]

        # Creates the nodes, coordinate and tag for every vertex of the polytope.
        # The tag is used to draw the front facets later on.
        dict_drawing = {}
        back_part = ''
        front_part = ''

        for vert in self.points:
            v = self.coords[vert]
            v_vect = str(['%.5f' % i for i in v]).replace('\'', '')
            v_vect = v_vect.replace('[', '(')
            v_vect = v_vect.replace(']', ')')
            tag = '%s' % v_vect
            node = "\\node[%s] at %s     {};\n" % ('vertex', tag)
            coord = '\\coordinate %s at %s;\n' % (tag, tag)
            dict_drawing[vert] = node, coord, tag

        # Separate the edges between back and front

        for index1, index2 in self.lines:
            # v1 = self.coords[index1]
            # v2 = self.coords[index2]

            H_v1 = set(self.parent_polyhedron.Vrepresentation(index1).incident())
            H_v2 = set(self.parent_polyhedron.Vrepresentation(index2).incident())
            H_v12 = [h for h in H_v1.intersection(H_v2) if facets.index(h) in back_facets]

            # The back edge has to be between two vertices in the Back
            # AND such that the 2 facets touching them are in the Back
            if index1 in back_vertices and index2 in back_vertices and len(H_v12) == 2:
                back_part += "\\draw[%s,back] %s -- %s;\n" % ('edge',
                                                              dict_drawing[index1][2],
                                                              dict_drawing[index2][2])
            else:
                front_part += "\\draw[%s] %s -- %s;\n" % ('edge',
                                                          dict_drawing[index1][2],
                                                          dict_drawing[index2][2])

        # Start to write the output
        tikz_pic = ''
        tikz_pic += '\\begin{tikzpicture}%\n'
        tikz_pic += '\t[x={(%fcm, %fcm)},\n' % (RDF(rotation_matrix[0][0]),
                                                RDF(rotation_matrix[0][1]))
        tikz_pic += '\ty={(%fcm, %fcm)},\n' % (RDF(rotation_matrix[1][0]),
                                               RDF(rotation_matrix[1][1]))
        tikz_pic += '\tz={(%fcm, %fcm)},\n' % (RDF(rotation_matrix[2][0]),
                                               RDF(rotation_matrix[2][1]))
        tikz_pic += '\tscale=%f,\n' % scale
        tikz_pic += '\tback/.style={loosely dotted, thin},\n'
        tikz_pic += '\tedge/.style={color=%s, thick},\n' % edge_color
        tikz_pic += '\tfacet/.style={fill=%s,fill opacity=%f},\n' % (facet_color, opacity)
        tikz_pic += '\tvertex/.style={inner sep=1pt,circle,draw=%s!25!black,' % vertex_color
        tikz_pic += 'fill=%s!75!black,thick}]\n%%\n%%\n' % vertex_color

        # Gives the reproduction information
        from sage.env import SAGE_VERSION
        tikz_pic += "%% This TikZ-picture was produced with Sagemath version {}\n".format(SAGE_VERSION)
        tikz_pic += "%% with the command: ._tikz_3d_in_3d and parameters:\n"
        tikz_pic += "%% view = {}\n".format(view)
        tikz_pic += "%% angle = {}\n".format(angle)
        tikz_pic += "%% scale = {}\n".format(scale)
        tikz_pic += "%% edge_color = {}\n".format(edge_color)
        tikz_pic += "%% facet_color = {}\n".format(facet_color)
        tikz_pic += "%% opacity = {}\n".format(opacity)
        tikz_pic += "%% vertex_color = {}\n".format(vertex_color)
        tikz_pic += "%% axis = {}\n%%\n".format(axis)

        # Draws the axes if True
        if axis:
            tikz_pic += '%% Drawing the axes\n'
            tikz_pic += '\\draw[color=black,thick,->] (0,0,0) -- (1,0,0) node[anchor=north east]{$x$};\n'
            tikz_pic += '\\draw[color=black,thick,->] (0,0,0) -- (0,1,0) node[anchor=north west]{$y$};\n'
            tikz_pic += '\\draw[color=black,thick,->] (0,0,0) -- (0,0,1) node[anchor=south]{$z$};\n'

        # Create the coordinate of the vertices
        tikz_pic += '%% Coordinate of the vertices:\n%%\n'
        for v in dict_drawing:
            tikz_pic += dict_drawing[v][1]

        # Draw the edges in the back
        tikz_pic += '%%\n%%\n%% Drawing edges in the back\n%%\n'
        tikz_pic += back_part

        # Draw the vertices on top of the back-edges
        tikz_pic += '%%\n%%\n%% Drawing vertices in the back\n%%\n'
        for v in back_vertices:
            if v not in front_vertices and v in dict_drawing:
                tikz_pic += dict_drawing[v][0]

        # Draw the facets in the front by going in cycles for every facet.
        tikz_pic += '%%\n%%\n%% Drawing the facets\n%%\n'
        vertex_to_index = {v: i for i, v in enumerate(vertices)}
        for index_facet in front_facets:
            cyclic_vert = cyclic_sort_vertices_2d(list(facets[index_facet].incident()))
            cyclic_indices = [vertex_to_index[v] for v in cyclic_vert]
            tikz_pic += '\\fill[facet] '
            for v in cyclic_indices:
                if v in dict_drawing:
                    tikz_pic += '%s -- ' % dict_drawing[v][2]
            tikz_pic += 'cycle {};\n'

        # Draw the edges in the front
        tikz_pic += '%%\n%%\n%% Drawing edges in the front\n%%\n'
        tikz_pic += front_part

        # Finally, the vertices in front are drawn on top of everything.
        tikz_pic += '%%\n%%\n%% Drawing the vertices in the front\n%%\n'
        for v in self.points:
            if v in front_vertices:
                if v in dict_drawing:
                    tikz_pic += dict_drawing[v][0]
        tikz_pic += '%%\n%%\n\\end{tikzpicture}'
        return LatexExpr(tikz_pic)
