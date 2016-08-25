"""
Functions for plotting polyhedra
"""

########################################################################
#       Copyright (C) 2008 Marshall Hampton <hamptonio@gmail.com>
#       Copyright (C) 2011 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################
from __future__ import print_function
from __future__ import absolute_import

from sage.rings.all import RDF
from sage.structure.sage_object import SageObject
from sage.modules.free_module_element import vector
from sage.matrix.constructor import matrix, identity_matrix
from sage.misc.functional import norm
from sage.misc.latex import LatexExpr
from sage.symbolic.constants import pi
from sage.structure.sequence import Sequence

from sage.plot.all import Graphics, point2d, line2d, arrow, polygon2d
from sage.plot.plot3d.all import point3d, line3d, arrow3d, polygon3d
from sage.plot.plot3d.transform import rotate_arbitrary

from .base import is_Polyhedron



#############################################################
def render_2d(projection, *args, **kwds):
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
        sage: q1.plot() + q2.plot() + q3.plot() + q4.plot()
        Graphics object consisting of 17 graphics primitives
        sage: from sage.geometry.polyhedron.plot import render_2d
        sage: q = render_2d(p1.projection())
        doctest:...: DeprecationWarning: use Projection.render_2d instead
        See http://trac.sagemath.org/16625 for details.
        sage: q._objects
        [Point set defined by 1 point(s),
         Arrow from (1.0,1.0) to (2.0,2.0),
         Polygon defined by 3 points]
    """
    from sage.misc.superseded import deprecation
    deprecation(16625, 'use Projection.render_2d instead')
    if is_Polyhedron(projection):
        projection = Projection(projection)
    return projection.render_2d(*args, **kwds)

def render_3d(projection, *args, **kwds):
    """
    Return 3d rendering of a polyhedron projected into
    3-dimensional ambient space.

    .. NOTE::

        This method, ``render_3d``, is used in the ``show()``
        method of a polyhedron if it is in 3 dimensions.

    EXAMPLES::

        sage: p1 = Polyhedron(vertices=[[1,1,1]], rays=[[1,1,1]])
        sage: p2 = Polyhedron(vertices=[[2,0,0], [0,2,0], [0,0,2]])
        sage: p3 = Polyhedron(vertices=[[1,0,0], [0,1,0], [0,0,1]], rays=[[-1,-1,-1]])
        sage: p1.projection().plot() + p2.projection().plot() + p3.projection().plot()
        Graphics3d Object

    It correctly handles various degenerate cases::

        sage: Polyhedron(lines=[[1,0,0],[0,1,0],[0,0,1]]).plot()                              # whole space
        Graphics3d Object
        sage: Polyhedron(vertices=[[1,1,1]], rays=[[1,0,0]], lines=[[0,1,0],[0,0,1]]).plot()  # half space
        Graphics3d Object
        sage: Polyhedron(vertices=[[1,1,1]], lines=[[0,1,0],[0,0,1]]).plot()                  # R^2 in R^3
        Graphics3d Object
        sage: Polyhedron(rays=[[0,1,0],[0,0,1]], lines=[[1,0,0]]).plot()                      # quadrant wedge in R^2
        Graphics3d Object
        sage: Polyhedron(rays=[[0,1,0]], lines=[[1,0,0]]).plot()                              # upper half plane in R^3
        Graphics3d Object
        sage: Polyhedron(lines=[[1,0,0]]).plot()                                              # R^1 in R^2
        Graphics3d Object
        sage: Polyhedron(rays=[[0,1,0]]).plot()                                               # Half-line in R^3
        Graphics3d Object
        sage: Polyhedron(vertices=[[1,1,1]]).plot()                                           # point in R^3
        Graphics3d Object
    """
    from sage.misc.superseded import deprecation
    deprecation(16625, 'use Projection.render_3d instead')
    if is_Polyhedron(projection):
        projection = Projection(projection)
    return projection.render_3d(*args, **kwds)

def render_4d(polyhedron, point_opts={}, line_opts={}, polygon_opts={}, projection_direction=None):
    """
    Return a 3d rendering of the Schlegel projection of a 4d
    polyhedron projected into 3-dimensional space.

    .. NOTE::

        The ``show()`` method of ``Polyhedron()`` uses this to draw itself
        if the ambient dimension is 4.

    INPUT:

    - ``polyhedron`` -- A
      :mod:`~sage.geometry.polyhedron.constructor.Polyhedron` object.

    - ``point_opts``, ``line_opts``, ``polygon_opts`` -- dictionaries
      of plot keywords or ``False`` to disable.

    - ``projection_direction`` -- list/tuple/iterable of coordinates
      or ``None`` (default). Sets the projetion direction of the
      Schlegel projection. If it is not given, the center of a facet
      is used.

    EXAMPLES::

        sage: poly = polytopes.twenty_four_cell()
        sage: poly
        A 4-dimensional polyhedron in QQ^4 defined as the convex hull of 24 vertices
        sage: poly.plot()
        Graphics3d Object
        sage: poly.plot(projection_direction=[2,5,11,17])
        Graphics3d Object
        sage: type( poly.plot() )
        <class 'sage.plot.plot3d.base.Graphics3dGroup'>

    TESTS::

        sage: from sage.geometry.polyhedron.plot import render_4d
        sage: p = polytopes.hypercube(4)
        sage: q = render_4d(p)
        doctest:...: DeprecationWarning: use Polyhedron.schlegel_projection instead
        See http://trac.sagemath.org/16625 for details.
        doctest:...: DeprecationWarning: use Projection.render_3d instead
        See http://trac.sagemath.org/16625 for details.
        sage: tach_str = q.tachyon()
        sage: tach_str.count('FCylinder')
        32
    """
    from sage.misc.superseded import deprecation
    deprecation(16625, 'use Polyhedron.schlegel_projection instead')
    if projection_direction is None:
        for ineq in polyhedron.inequality_generator():
            center = [v() for v in ineq.incident() if v.is_vertex()]
            center = sum(center) / len(center)
            if not center.is_zero():
                projection_direction = center
                break
    projection_3d = Projection(polyhedron).schlegel(projection_direction)
    return render_3d(projection_3d, point_opts, line_opts, polygon_opts)



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
    chain = [ Vlist.pop() ]
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
            chain = [ Vlist.pop() ]
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



class ProjectionFuncStereographic():
    """
    The stereographic (or perspective) projection.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.plot import ProjectionFuncStereographic
        sage: cube = polytopes.hypercube(3).vertices()
        sage: proj = ProjectionFuncStereographic([1.2, 3.4, 5.6])
        sage: ppoints = [proj(vector(x)) for x in cube]
        sage: ppoints[0]
        (-0.0486511..., 0.0859565...)
    """
    def __init__(self, projection_point):
        """
        Create a stereographic projection function.

        INPUT:

        - ``projection_point`` -- a list of coordinates in the
          appropriate dimension, which is the point projected from.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.plot import ProjectionFuncStereographic
            sage: proj = ProjectionFuncStereographic([1.0,1.0])
            sage: proj.__init__([1.0,1.0])
            sage: proj.house
            [-0.7071067811...  0.7071067811...]
            [ 0.7071067811...  0.7071067811...]
            sage: TestSuite(proj).run(skip='_test_pickling')
        """
        self.projection_point = vector(projection_point)
        self.dim = self.projection_point.degree()

        pproj = vector(RDF, self.projection_point)
        self.psize = norm(pproj)
        if (self.psize).is_zero():
            raise ValueError("projection direction must be a non-zero vector.")
        v = vector(RDF, [0.0]*(self.dim-1) + [self.psize]) - pproj
        polediff = matrix(RDF, v).transpose()
        denom = RDF((polediff.transpose()*polediff)[0][0])
        if denom.is_zero():
            self.house = identity_matrix(RDF, self.dim)
        else:
            self.house = identity_matrix(RDF, self.dim) \
            - 2*polediff*polediff.transpose()/denom   # Householder reflector

    def __call__(self, x):
        """
        Action of the stereographic projection.

        INPUT:

        - ``x`` -- a vector or anything convertible to a vector.

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
            (-1.0000000000000002)
            sage: proj = ProjectionFuncStereographic([2.0,1.0])
            sage: proj.__call__(vector([1,2]))  # abs tol 1e-14
            (2.9999999999999996)
            sage: proj = ProjectionFuncStereographic([0,0,2])
            sage: proj.__call__(vector([0,0,1]))
            (0.0, 0.0)
            sage: proj.__call__(vector([1,0,0]))
            (0.5, 0.0)
        """
        img = self.house * x
        denom = self.psize-img[self.dim-1]
        if denom.is_zero():
            raise ValueError('Point cannot coincide with ' \
                'coordinate singularity at ' + repr(x))
        return vector(RDF, [img[i]/denom for i in range(self.dim-1)])


class ProjectionFuncSchlegel():
    """
    The Schlegel projection from the given input point.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.plot import ProjectionFuncSchlegel
        sage: proj = ProjectionFuncSchlegel([2,2,2])
        sage: proj(vector([1.1,1.1,1.11]))[0]
        0.0302...
    """
    def __init__(self, projection_direction, height=1.1, center=0):
        """
        Initializes the projection.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.plot import ProjectionFuncSchlegel
            sage: proj = ProjectionFuncSchlegel([2,2,2])
            sage: proj.__init__([2,2,2])
            sage: proj(vector([1.1,1.1,1.11]))[0]
            0.0302...
            sage: TestSuite(proj).run(skip='_test_pickling')
        """
        self.center = center
        self.projection_dir = vector(RDF, projection_direction)
        if norm(self.projection_dir).is_zero():
            raise ValueError("projection direction must be a non-zero vector.")
        self.dim = self.projection_dir.degree()
        spcenter = height * self.projection_dir/norm(self.projection_dir)
        self.height = height
        v = vector(RDF, [0.0]*(self.dim-1) + [self.height]) - spcenter
        polediff = matrix(RDF, v).transpose()
        denom = (polediff.transpose()*polediff)[0][0]
        if denom.is_zero():
            self.house = identity_matrix(RDF, self.dim)
        else:
            self.house = identity_matrix(RDF, self.dim) \
            - 2*polediff*polediff.transpose()/denom  # Householder reflector

    def __call__(self, x):
        """
        Apply the projection to a vector.

        - ``x`` -- a vector or anything convertible to a vector.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.plot import ProjectionFuncSchlegel
            sage: proj = ProjectionFuncSchlegel([2,2,2])
            sage: proj.__call__([1,2,3])
            (0.56162854..., 2.09602626...)
        """
        v = vector(RDF, x) - self.center
        if v.is_zero():
            raise ValueError("The origin must not be a vertex.")
        v = v/norm(v)         # normalize vertices to unit sphere
        v = self.house*v      # reflect so self.projection_dir is at "north pole"
        denom = self.height-v[self.dim-1]
        if denom.is_zero():
            raise ValueError('Point cannot coincide with ' \
                'coordinate singularity at ' + repr(x))
        return vector(RDF, [ v[i]/denom for i in range(self.dim-1) ])



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

            sage: p = polytopes.icosahedron(exact=False)
            sage: from sage.geometry.polyhedron.plot import Projection
            sage: Projection(p)
            The projection of a polyhedron into 3 dimensions
            sage: def pr_12(x): return [x[1],x[2]]
            sage: Projection(p, pr_12)
            The projection of a polyhedron into 2 dimensions
            sage: Projection(p,  lambda x: [x[1],x[2]] )   # another way of doing the same projection
            The projection of a polyhedron into 2 dimensions
            sage: _.plot()   # plot of the projected icosahedron in 2d
            Graphics object consisting of 51 graphics primitives
            sage: proj = Projection(p)
            sage: proj.stereographic([1,2,3])
            The projection of a polyhedron into 2 dimensions
            sage: proj.plot()
            Graphics object consisting of 51 graphics primitives
            sage: TestSuite(proj).run(skip='_test_pickling')
        """
        self.parent_polyhedron = polyhedron
        self.coords = Sequence([])
        self.points = Sequence([])
        self.lines  = Sequence([])
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
        s = 'The projection of a polyhedron into ' + \
            repr(self.dimension) + ' dimensions'
        return s


    def __call__(self, proj=projection_func_identity):
        """
        Apply a projection.

        EXAMPLES::

            sage: p = polytopes.icosahedron(exact=False)
            sage: from sage.geometry.polyhedron.plot import Projection
            sage: pproj = Projection(p)
            sage: from sage.geometry.polyhedron.plot import ProjectionFuncStereographic
            sage: pproj_stereo = pproj.__call__(proj = ProjectionFuncStereographic([1,2,3]))
            sage: pproj_stereo.polygons[0]
            [10, 1, 4]
        """
        self.transformed_coords = \
            Sequence([proj(p) for p in self.coords])
        self._init_dimension()
        return self


    def identity(self):
        """
        Return the identity projection of the polyhedron.

        EXAMPLES::

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

        - ``projection_point`` - The projection point. This must be
          distinct from the polyhedron's vertices. Default is `(1,0,\dots,0)`

        EXAMPLES::

            sage: from sage.geometry.polyhedron.plot import Projection
            sage: proj = Projection(polytopes.buckyball())  #long time
            sage: proj                                      #long time
            The projection of a polyhedron into 3 dimensions
            sage: proj.stereographic([5,2,3]).plot()        #long time
            Graphics object consisting of 123 graphics primitives
            sage: Projection( polytopes.twenty_four_cell() ).stereographic([2,0,0,0])
            The projection of a polyhedron into 3 dimensions
        """
        if projection_point is None:
            projection_point = [1] + [0]*(self.polyhedron_ambient_dim-1)
        return self(ProjectionFuncStereographic(projection_point))


    def schlegel(self, projection_direction=None, height=1.1):
        """
        Return the Schlegel projection.

        * The polyhedron is translated such that its
          :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.center`
          is at the origin.

        * The vertices are then normalized to the unit sphere

        * The normalized points are stereographically projected from a
          point slightly outside of the sphere.

        INPUT:

        - ``projection_direction`` -- coordinate list/tuple/iterable
          or ``None`` (default). The direction of the Schlegel
          projection. For a full-dimensional polyhedron, the default
          is the first facet normal; Otherwise, the vector consisting
          of the first n primes is chosen.

        - ``height`` -- float (default: `1.1`). How far outside of the
          unit sphere the focal point is.

        EXAMPLES::

            sage: cube4 = polytopes.hypercube(4)
            sage: from sage.geometry.polyhedron.plot import Projection
            sage: Projection(cube4).schlegel([1,0,0,0])
            The projection of a polyhedron into 3 dimensions
            sage: _.plot()
            Graphics3d Object

        TESTS::

            sage: Projection(cube4).schlegel()
            The projection of a polyhedron into 3 dimensions
        """
        center = self.parent_polyhedron.center()
        if projection_direction is None:
            if self.parent_polyhedron.is_full_dimensional():
                projection_direction = next(self.parent_polyhedron.inequality_generator()).A()
            else:
                from sage.arith.all import primes_first_n
                projection_direction = primes_first_n(self.polyhedron_ambient_dim)
        return self(ProjectionFuncSchlegel(
            projection_direction, height=height, center=center))

    def coord_index_of(self, v):
        """
        Convert a coordinate vector to its internal index.

        EXAMPLES::

            sage: p = polytopes.hypercube(3)
            sage: proj = p.projection()
            sage: proj.coord_index_of(vector((1,1,1)))
            7
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
            sage: proj.coord_indices_of([vector((1,1,1)),vector((1,-1,1))])
            [7, 5]
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

            sage: from sage.geometry.polyhedron.plot import Projection, render_2d
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

    def show(self, *args, **kwds):
        """
        Deprecated method to show the projection as a graphics
        object.  Use ``Projection.plot()`` instead.

        EXAMPLE::

            sage: P8 = polytopes.hypercube(4)
            sage: P8.schlegel_projection([2,5,11,17]).show()
            doctest:...: DeprecationWarning: use Projection.plot instead
            See http://trac.sagemath.org/16625 for details.
            Graphics3d Object
        """
        from sage.misc.superseded import deprecation
        deprecation(16625, 'use Projection.plot instead')
        return self.plot(*args, **kwds)

    def _init_from_2d(self, polyhedron):
        """
        Internal function: Initialize from polyhedron in
        2-dimensional space. The polyhedron could be lower
        dimensional.

        TESTS::

            sage: p = Polyhedron(vertices = [[0,0],[0,1],[1,0],[1,1]])
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

            sage: p = Polyhedron(vertices = [[0,0,1],[0,1,2],[1,0,3],[1,1,5]])
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
            self.points.append( self.coord_index_of(v.vector()) )


    def _init_lines_arrows(self, polyhedron):
        """
        Internal function: Initialize compact and non-compact edges
        (works in arbitrary dimensions).

        TESTS::

            sage: p = Polyhedron(ieqs = [[1, 0, 0, 1],[1,1,0,0]])
            sage: pp = p.projection()
            sage: pp.arrows
            [[0, 1], [0, 2]]
            sage: del pp.arrows
            sage: pp.arrows = Sequence([])
            sage: pp._init_lines_arrows(p)
            sage: pp.arrows
            [[0, 1], [0, 2]]
        """
        obj = polyhedron.Vrepresentation()
        for i in range(len(obj)):
            if not obj[i].is_vertex(): continue
            for j in range(len(obj)):
                if polyhedron.vertex_adjacency_matrix()[i,j] == 0: continue
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
        vertices = [v for v in polyhedron.Vrep_generator()]
        vertices = cyclic_sort_vertices_2d(vertices)
        coords = []

        def adjacent_vertices(i):
            n = len(vertices)
            if vertices[(i-1) % n].is_vertex(): yield vertices[(i-1) % n]
            if vertices[(i+1) % n].is_vertex(): yield vertices[(i+1) % n]

        for i in range(len(vertices)):
            v = vertices[i]
            if v.is_vertex():
                coords.append(v())
            if v.is_ray():
                for a in adjacent_vertices(i):
                    coords.append(a() + v())

        if polyhedron.n_lines() == 0:
            self.polygons.append( self.coord_indices_of(coords) )
            return

        polygons = []

        if polyhedron.n_lines() == 1:
            aline = next(polyhedron.line_generator())
            for shift in [aline(), -aline()]:
                for i in range(len(coords)):
                    polygons.append( [ coords[i-1],coords[i],
                                       coords[i]+shift, coords[i-1]+shift ] )

        if polyhedron.n_lines() == 2:
            [line1, line2] = [l for l in polyhedron.lines()]
            assert len(coords) == 1, "Can have only a single vertex!"
            v = coords[0]
            l1 = line1()
            l2 = line2()
            polygons = [ [v-l1-l2, v+l1-l2, v+l1+l2, v-l1+l2] ]

        polygons = [ self.coord_indices_of(p) for p in polygons ]
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
            [[2, 0, 1], [3, 0, 1], [3, 0, 2], [3, 1, 2]]
        """
        assert polyhedron.ambient_dim() == 3, "Requires polyhedron in 3d"

        if polyhedron.dim() <= 1:  # empty or 0d or 1d polyhedron => no polygon
            return None

        def defining_equation():  # corresponding to a polygon
            if polyhedron.dim() < 3:
                yield next(polyhedron.equation_generator())
            else:
                for ineq in polyhedron.inequality_generator():
                    yield ineq

        faces = []
        face_inequalities = []
        for facet_equation in defining_equation():
            vertices = [v for v in facet_equation.incident()]
            face_inequalities.append(facet_equation)
            vertices = cyclic_sort_vertices_2d(vertices)
            coords = []

            def adjacent_vertices(i):
                n = len(vertices)
                if vertices[(i-1) % n].is_vertex(): yield vertices[(i-1) % n]
                if vertices[(i+1) % n].is_vertex(): yield vertices[(i+1) % n]

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
            assert len(faces) > 0, "no vertices?"
            self.polygons.extend( [self.coord_indices_of(f) for f in faces] )
            return

        # now some special cases if there are lines (dim < ambient_dim)
        polygons = []

        if polyhedron.n_lines() == 1:
            assert len(faces) > 0, "no vertices?"
            aline = next(polyhedron.line_generator())
            for shift in [aline(), -aline()]:
                for coords in faces:
                    assert len(coords) == 2, "There must be two points."
                    polygons.append( [ coords[0],coords[1],
                                       coords[1]+shift, coords[0]+shift ] )

        if polyhedron.n_lines() == 2:
            [line1, line2] = [l for l in polyhedron.line_generator()]
            l1 = line1()
            l2 = line2()
            for v in polyhedron.vertex_generator():
                polygons.append( [v()-l1-l2, v()+l1-l2, v()+l1+l2, v()-l1+l2] )

        self.polygons.extend( [self.coord_indices_of(p) for p in polygons] )


    def render_points_1d(self, **kwds):
        """
        Return the points of a polyhedron in 1d.

        INPUT:

        - ``**kwds`` -- options passed through to
          :func:`~sage.plot.point.point2d`.

        OUTPUT:

        A 2-d graphics object.

        EXAMPLES::

            sage: cube1 = polytopes.hypercube(1)
            sage: proj = cube1.projection()
            sage: points = proj.render_points_1d()
            sage: points._objects
            [Point set defined by 2 point(s)]
        """
        return point2d([c + [0] for c in self.coordinates_of(self.points)], **kwds)

    def render_line_1d(self, **kwds):
        """
        Return the line of a polyhedron in 1d.

        INPUT:

        - ``**kwds`` -- options passed through to
          :func:`~sage.plot.line.line2d`.

        OUTPUT:

        A 2-d graphics object.

        EXAMPLES::

            sage: outline = polytopes.hypercube(1).projection().render_line_1d()
            sage: outline._objects[0]
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

            sage: hex = polytopes.regular_polygon(6)
            sage: proj = hex.projection()
            sage: hex_points = proj.render_points_2d()
            sage: hex_points._objects
            [Point set defined by 6 point(s)]
        """
        return point2d(self.coordinates_of(self.points), **kwds)


    def render_outline_2d(self, **kwds):
        """
        Return the outline (edges) of a polyhedron in 2d.

        EXAMPLES::

            sage: penta = polytopes.regular_polygon(5)
            sage: outline = penta.projection().render_outline_2d()
            sage: outline._objects[0]
            Line defined by 2 points
        """
        wireframe = [];
        for l in self.lines:
            l_coords = self.coordinates_of(l)
            wireframe.append( line2d(l_coords, **kwds) )
        for a in self.arrows:
            a_coords = self.coordinates_of(a)
            wireframe.append( arrow(a_coords[0], a_coords[1], **kwds) )
        return sum(wireframe)


    def render_fill_2d(self, **kwds):
        """
        Return the filled interior (a polygon) of a polyhedron in 2d.

        EXAMPLES::

            sage: cps = [i^3 for i in srange(-2,2,1/5)]
            sage: p = Polyhedron(vertices = [[(t^2-1)/(t^2+1),2*t/(t^2+1)] for t in cps])
            sage: proj = p.projection()
            sage: filled_poly = proj.render_fill_2d()
            sage: filled_poly.axes_width()
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
            sage: verts = proj.render_vertices_3d()
            sage: verts.bounding_box()
            ((-1.0, -1.0, -1.0), (1.0, 1.0, 1.0))
        """
        return point3d(self.coordinates_of(self.points), **kwds)


    def render_wireframe_3d(self, **kwds):
        r"""
        Return the 3d wireframe rendering.

        EXAMPLES::

            sage: cube = polytopes.hypercube(3)
            sage: cube_proj = cube.projection()
            sage: wire = cube_proj.render_wireframe_3d()
            sage: print(wire.tachyon().split('\n')[77])  # for testing
            FCylinder base -1.0 1.0 -1.0 apex -1.0 -1.0 -1.0 rad 0.005 texture...
        """
        wireframe = []
        for l in self.lines:
            l_coords = self.coordinates_of(l)
            wireframe.append( line3d(l_coords, **kwds))
        for a in self.arrows:
            a_coords = self.coordinates_of(a)
            wireframe.append(arrow3d(a_coords[0], a_coords[1], **kwds))
        return sum(wireframe)

    def render_solid_3d(self, **kwds):
        """
        Return solid 3d rendering of a 3d polytope.

        EXAMPLES::

            sage: p = polytopes.hypercube(3).projection()
            sage: p_solid = p.render_solid_3d(opacity = .7)
            sage: type(p_solid)
            <class 'sage.plot.plot3d.base.Graphics3dGroup'>
        """
        return sum([ polygon3d(self.coordinates_of(f), **kwds)
                     for f in self.polygons ])

    def render_0d(self, point_opts={}, line_opts={}, polygon_opts={}):
        """
        Return 0d rendering of the projection of a polyhedron into
        2-dimensional ambient space.

        INPUT:

        See
        :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.plot`.

        OUTPUT:

        A 2-d graphics object.

        EXAMPLES::

            sage: print(Polyhedron([]).projection().render_0d().description())
            <BLANKLINE>
            sage: print(Polyhedron(ieqs=[(1,)]).projection().render_0d().description())
            Point set defined by 1 point(s):    [(0.0, 0.0)]
        """
        if isinstance(point_opts, dict):
            point_opts.setdefault('zorder', 2)
            point_opts.setdefault('pointsize', 10)
        if self.points:
            return point2d([0,0], **point_opts)
        else:
            return Graphics()


    def render_1d(self, point_opts={}, line_opts={}, polygon_opts={}):
        """
        Return 1d rendering of the projection of a polyhedron into
        2-dimensional ambient space.

        INPUT:

        See
        :meth:`~sage.geometry.polyhedron.base.Polyhedron_base.plot`.

        OUTPUT:

        A 2-d graphics object.

        EXAMPLES::

            sage: Polyhedron([(0,), (1,)]).projection().render_1d()
            Graphics object consisting of 2 graphics primitives
        """
        plt = Graphics()
        if isinstance(point_opts, dict):
            point_opts.setdefault('zorder', 2)
            point_opts.setdefault('pointsize', 10)
            plt += self.render_points_1d(**point_opts)
        if isinstance(line_opts, dict):
            line_opts.setdefault('zorder', 1)
            plt += self.render_line_1d(**line_opts)
        return plt

    def render_2d(self, point_opts={}, line_opts={}, polygon_opts={}):
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
            sage: q1.plot() + q2.plot() + q3.plot() + q4.plot()
            Graphics object consisting of 17 graphics primitives
         """
        plt = Graphics()
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

    def render_3d(self, point_opts={}, line_opts={}, polygon_opts={}):
        """
        Return 3d rendering of a polyhedron projected into
        3-dimensional ambient space.

        EXAMPLES::

            sage: p1 = Polyhedron(vertices=[[1,1,1]], rays=[[1,1,1]])
            sage: p2 = Polyhedron(vertices=[[2,0,0], [0,2,0], [0,0,2]])
            sage: p3 = Polyhedron(vertices=[[1,0,0], [0,1,0], [0,0,1]], rays=[[-1,-1,-1]])
            sage: p1.projection().plot() + p2.projection().plot() + p3.projection().plot()
            Graphics3d Object

        It correctly handles various degenerate cases::

            sage: Polyhedron(lines=[[1,0,0],[0,1,0],[0,0,1]]).plot()           # whole space
            Graphics3d Object
            sage: Polyhedron(vertices=[[1,1,1]], rays=[[1,0,0]],
            ....:            lines=[[0,1,0],[0,0,1]]).plot()                   # half space
            Graphics3d Object
            sage: Polyhedron(vertices=[[1,1,1]],
            ....:            lines=[[0,1,0],[0,0,1]]).plot()                   # R^2 in R^3
            Graphics3d Object
            sage: Polyhedron(rays=[[0,1,0],[0,0,1]], lines=[[1,0,0]]).plot()   # quadrant wedge in R^2
            Graphics3d Object
            sage: Polyhedron(rays=[[0,1,0]], lines=[[1,0,0]]).plot()           # upper half plane in R^3
            Graphics3d Object
            sage: Polyhedron(lines=[[1,0,0]]).plot()                           # R^1 in R^2
            Graphics3d Object
            sage: Polyhedron(rays=[[0,1,0]]).plot()                            # Half-line in R^3
            Graphics3d Object
            sage: Polyhedron(vertices=[[1,1,1]]).plot()                        # point in R^3
            Graphics3d Object
        """
        from sage.plot.plot3d.base import Graphics3d
        plt = Graphics3d()
        if isinstance(point_opts, dict):
            point_opts.setdefault('width', 3)
            plt += self.render_vertices_3d(**point_opts)
        if isinstance(line_opts, dict):
            line_opts.setdefault('width', 3)
            plt += self.render_wireframe_3d(**line_opts)
        if isinstance(polygon_opts, dict):
            plt += self.render_solid_3d(**polygon_opts)
        return plt

    def tikz(self, view=[0, 0, 1], angle=0, scale=2,
             edge_color='blue!95!black', facet_color='blue!95!black',
             opacity=0.8, vertex_color='green', axis=False):
        r"""
        Return a string ``tikz_pic`` consisting of a tikz picture of ``self``
        according to a projection ``view`` and an angle ``angle``
        obtained via Jmol through the current state property.

        INPUT:

        - ``view`` - list (default: [0,0,1]) representing the rotation axis (see note below).
        - ``angle`` - integer (default: 0) angle of rotation in degree from 0 to 360 (see note
          below).
        - ``scale`` - integer (default: 2) specifying the scaling of the tikz picture.
        - ``edge_color`` - string (default: 'blue!95!black') representing colors which tikz
          recognize.
        - ``facet_color`` - string (default: 'blue!95!black') representing colors which tikz
          recognize.
        - ``vertex_color`` - string (default: 'green') representing colors which tikz
          recognize.
        - ``opacity`` - real number (default: 0.8) between 0 and 1 giving the opacity of
          the front facets.
        - ``axis`` - Boolean (default: False) draw the axes at the origin or not.

        OUTPUT:

        - LatexExpr -- containing the TikZ picture.

        .. NOTE::

            The inputs ``view`` and ``angle`` can be obtained from the
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

            sage: P1 = polytopes.small_rhombicuboctahedron()
            sage: Image1 = P1.projection().tikz([1,3,5], 175, scale=4)
            sage: type(Image1)
            <class 'sage.misc.latex.LatexExpr'>
            sage: print('\n'.join(Image1.splitlines()[:4]))
            \begin{tikzpicture}%
                [x={(-0.939161cm, 0.244762cm)},
                y={(0.097442cm, -0.482887cm)},
                z={(0.329367cm, 0.840780cm)},
            sage: open('polytope-tikz1.tex', 'w').write(Image1)    # not tested

            sage: P2 = Polyhedron(vertices=[[1, 1],[1, 2],[2, 1]])
            sage: Image2 = P2.projection().tikz(scale=3, edge_color='blue!95!black', facet_color='orange!95!black', opacity=0.4, vertex_color='yellow', axis=True)
            sage: type(Image2)
            <class 'sage.misc.latex.LatexExpr'>
            sage: print('\n'.join(Image2.splitlines()[:4]))
            \begin{tikzpicture}%
                [scale=3.000000,
                back/.style={loosely dotted, thin},
                edge/.style={color=blue!95!black, thick},
            sage: open('polytope-tikz2.tex', 'w').write(Image2)    # not tested

            sage: P3 = Polyhedron(vertices=[[-1, -1, 2],[-1, 2, -1],[2, -1, -1]])
            sage: P3
            A 2-dimensional polyhedron in ZZ^3 defined as the convex hull of 3 vertices
            sage: Image3 = P3.projection().tikz([0.5,-1,-0.1], 55, scale=3, edge_color='blue!95!black',facet_color='orange!95!black', opacity=0.7, vertex_color='yellow', axis=True)
            sage: print('\n'.join(Image3.splitlines()[:4]))
            \begin{tikzpicture}%
                [x={(0.658184cm, -0.242192cm)},
                y={(-0.096240cm, 0.912008cm)},
                z={(-0.746680cm, -0.331036cm)},
            sage: open('polytope-tikz3.tex', 'w').write(Image3)    # not tested

            sage: P=Polyhedron(vertices=[[1,1,0,0],[1,2,0,0],[2,1,0,0],[0,0,1,0],[0,0,0,1]])
            sage: P
            A 4-dimensional polyhedron in ZZ^4 defined as the convex hull of 5 vertices
            sage: P.projection().tikz()
            Traceback (most recent call last):
            ...
            NotImplementedError: The polytope has to live in 2 or 3 dimensions.

        .. TODO::

            Make it possible to draw Schlegel diagram for 4-polytopes. ::

                sage: P=Polyhedron(vertices=[[1,1,0,0],[1,2,0,0],[2,1,0,0],[0,0,1,0],[0,0,0,1]])
                sage: P
                A 4-dimensional polyhedron in ZZ^4 defined as the convex hull of 5 vertices
                sage: P.projection().tikz()
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
            return self._tikz_2d(scale, edge_color, facet_color, opacity,
                                 vertex_color, axis)
        elif self.polyhedron_dim == 2:  # self is a polygon in 3-space
            return self._tikz_2d_in_3d(view, angle, scale, edge_color,
                                       facet_color, opacity, vertex_color, axis)
        else:  # self is a 3-polytope in 3-space
            return self._tikz_3d_in_3d(view, angle, scale, edge_color,
                                       facet_color, opacity, vertex_color, axis)

    def _tikz_2d(self, scale, edge_color, facet_color, opacity, vertex_color, axis):
        r"""
        Return a string ``tikz_pic`` consisting of a tikz picture of
        ``self``, which is assumed to be a polygon on the plane.

        INPUT:

        - ``scale`` - integer specifying the scaling of the tikz picture.
        - ``edge_color`` - string representing colors which tikz
          recognize.
        - ``facet_color`` - string representing colors which tikz
          recognize.
        - ``vertex_color`` - string representing colors which tikz
          recognize.
        - ``opacity`` - real number between 0 and 1 giving the opacity of
          the front facets.
        - ``axis`` - Boolean (default: False) draw the axes at the origin or not.

        OUTPUT:

        - LatexExpr -- containing the TikZ picture.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[[1, 1],[1, 2],[2, 1]])
            sage: Image = P.projection()._tikz_2d(scale=3, edge_color='black', facet_color='orange', opacity=0.75, vertex_color='yellow', axis=True)
            sage: type(Image)
            <class 'sage.misc.latex.LatexExpr'>
            sage: print('\n'.join(Image.splitlines()[:4]))
            \begin{tikzpicture}%
                [scale=3.000000,
                back/.style={loosely dotted, thin},
                edge/.style={color=black, thick},
            sage: open('polytope-tikz2.tex', 'w').write(Image)    # not tested


        Scientific notation is not used in the output (:trac:`16519`)::

            sage: P=Polyhedron([[2*10^-10,0],[0,1],[1,0]],base_ring=QQ)
            sage: tikzstr=P.projection().tikz()
            sage: 'e-10' in tikzstr
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
            v_vect = str(['%.5f' % i for i in v]).replace('\'','')
            v_vect = v_vect.replace('[', '(')
            v_vect = v_vect.replace(']', ')')
            tag = '%s' %v_vect
            node = "\\node[%s] at %s     {};\n" % ('vertex', tag)
            coord = '\coordinate %s at %s;\n' % (tag, tag)
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
        tikz_pic += '\tfacet/.style={fill=%s,fill opacity=%f},\n' % (facet_color,opacity)
        tikz_pic += '\tvertex/.style={inner sep=1pt,circle,draw=%s!25!black,' % vertex_color
        tikz_pic += 'fill=%s!75!black,thick,anchor=base}]\n%%\n%%\n' % vertex_color

        # Draws the axes if True
        if axis:
            tikz_pic += '\\draw[color=black,thick,->] (0,0,0) -- (1,0,0) node[anchor=north east]{$x$};\n'
            tikz_pic += '\\draw[color=black,thick,->] (0,0,0) -- (0,1,0) node[anchor=north west]{$y$};\n'

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

        - ``view`` - list (default: [0,0,1]) representing the rotation axis.
        - ``angle`` - integer angle of rotation in degree from 0 to 360.
        - ``scale`` - integer specifying the scaling of the tikz picture.
        - ``edge_color`` - string representing colors which tikz
          recognize.
        - ``facet_color`` - string representing colors which tikz
          recognize.
        - ``vertex_color`` - string representing colors which tikz
          recognize.
        - ``opacity`` - real number between 0 and 1 giving the opacity of
          the front facets.
        - ``axis`` - Boolean draw the axes at the origin or not.

        OUTPUT:

        - LatexExpr -- containing the TikZ picture.

        EXAMPLE::

            sage: P = Polyhedron(vertices=[[-1, -1, 2],[-1, 2, -1],[2, -1, -1]])
            sage: P
            A 2-dimensional polyhedron in ZZ^3 defined as the convex hull of 3 vertices
            sage: Image = P.projection()._tikz_2d_in_3d(view=[0.5,-1,-0.5], angle=55, scale=3, edge_color='blue!95!black', facet_color='orange', opacity=0.5, vertex_color='yellow', axis=True)
            sage: print('\n'.join(Image.splitlines()[:4]))
            \begin{tikzpicture}%
                [x={(0.644647cm, -0.476559cm)},
                y={(0.192276cm, 0.857859cm)},
                z={(-0.739905cm, -0.192276cm)},
            sage: open('polytope-tikz3.tex', 'w').write(Image)    # not tested

            sage: p = Polyhedron(vertices=[[1,0,0],[0,1,0],[0,0,1]])
            sage: proj = p.projection()
            sage: Img = proj.tikz([1,1,1],130,axis=True)
            sage: print('\n'.join(Img.splitlines()[21:25]))
            %% Drawing the interior
            %%
            \fill[facet] (1.00000, 0.00000, 0.00000) -- (0.00000, 0.00000, 1.00000) -- (0.00000, 1.00000, 0.00000) -- cycle {};
            %%

        .. NOTE::

            The ``facet_color`` is the filing color of the polytope (polygon).
        """
        view_vector = vector(RDF, view)
        rot = rotate_arbitrary(view_vector, -(angle/360)*2*pi)
        rotation_matrix = rot[:2].transpose()

        # Creates the nodes, coordinate and tag for every vertex of the polytope.
        # The tag is used to draw the front facets later on.
        dict_drawing = {}
        edges = ''
        for vert in self.points:
            v = self.coords[vert]
            v_vect = str(['%.5f' % i for i in v]).replace('\'','')
            v_vect = v_vect.replace('[','(')
            v_vect = v_vect.replace(']',')')
            tag = '%s' %v_vect
            node = "\\node[%s] at %s     {};\n" % ('vertex', tag)
            coord = '\coordinate %s at %s;\n' % (tag, tag)
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
        tikz_pic += '\tfacet/.style={fill=%s,fill opacity=%f},\n' % (facet_color,opacity)
        tikz_pic += '\tvertex/.style={inner sep=1pt,circle,draw=%s!25!black,' % vertex_color
        tikz_pic += 'fill=%s!75!black,thick,anchor=base}]\n%%\n%%\n' % vertex_color

        # Draws the axes if True
        if axis:
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

        - ``view`` - list (default: [0,0,1]) representing the rotation axis.
        - ``angle`` - integer angle of rotation in degree from 0 to 360.
        - ``scale`` - integer specifying the scaling of the tikz picture.
        - ``edge_color`` - string representing colors which tikz
          recognize.
        - ``facet_color`` - string representing colors which tikz
          recognize.
        - ``vertex_color`` - string representing colors which tikz
          recognize.
        - ``opacity`` - real number between 0 and 1 giving the opacity of
          the front facets.
        - ``axis`` - Boolean draw the axes at the origin or not.

        OUTPUT:

        - LatexExpr -- containing the TikZ picture.

        EXAMPLES::

            sage: P = polytopes.small_rhombicuboctahedron()
            sage: Image = P.projection()._tikz_3d_in_3d([3,7,5], 100, scale=3, edge_color='blue', facet_color='orange', opacity=0.5, vertex_color='green', axis=True)
            sage: type(Image)
            <class 'sage.misc.latex.LatexExpr'>
            sage: print('\n'.join(Image.splitlines()[:4]))
            \begin{tikzpicture}%
                [x={(-0.046385cm, 0.837431cm)},
                y={(-0.243536cm, 0.519228cm)},
                z={(0.968782cm, 0.170622cm)},
            sage: open('polytope-tikz1.tex', 'w').write(Image)    # not tested

            sage: Associahedron = Polyhedron(vertices=[[1,0,1],[1,0,0],[1,1,0],[0,0,-1],[0,1,0],[-1,0,0],[0,1,1],[0,0,1],[0,-1,0]]).polar()
            sage: ImageAsso = Associahedron.projection().tikz([-15,-755,-655], 116, scale=1)
            sage: print('\n'.join(ImageAsso.splitlines()[29:41]))
            %% Drawing edges in the back
            %%
            \draw[edge,back] (-0.50000, -0.50000, -0.50000) -- (-1.00000, 0.00000, 0.00000);
            \draw[edge,back] (-0.50000, -0.50000, -0.50000) -- (0.00000, -1.00000, 0.00000);
            \draw[edge,back] (-0.50000, -0.50000, -0.50000) -- (0.00000, 0.00000, -1.00000);
            \draw[edge,back] (-1.00000, 0.00000, 0.00000) -- (-1.00000, 0.00000, 1.00000);
            \draw[edge,back] (-1.00000, 0.00000, 0.00000) -- (-1.00000, 1.00000, 0.00000);
            \draw[edge,back] (0.00000, -1.00000, 0.00000) -- (0.00000, -1.00000, 1.00000);
            \draw[edge,back] (0.00000, -1.00000, 0.00000) -- (1.00000, -1.00000, 0.00000);
            \draw[edge,back] (0.00000, 0.00000, -1.00000) -- (0.00000, 1.00000, -1.00000);
            \draw[edge,back] (0.00000, 0.00000, -1.00000) -- (1.00000, 0.00000, -1.00000);
            %%
        """
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
                if A*self.coords[v]+B < 0.0005 and v not in back_vertices:
                    back_vertices += [v]

        # Creates the nodes, coordinate and tag for every vertex of the polytope.
        # The tag is used to draw the front facets later on.
        dict_drawing = {}
        back_part = ''
        front_part = ''

        for vert in self.points:
            v = self.coords[vert]
            v_vect = str(['%.5f' % i for i in v]).replace('\'','')
            v_vect = v_vect.replace('[','(')
            v_vect = v_vect.replace(']',')')
            tag = '%s' %v_vect
            node = "\\node[%s] at %s     {};\n" % ('vertex', tag)
            coord = '\coordinate %s at %s;\n' %(tag, tag)
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
            if index1 in back_vertices and index2 in back_vertices and len(H_v12)==2:
                back_part += "\\draw[%s,back] %s -- %s;\n" % ('edge', dict_drawing[index1][2], dict_drawing[index2][2])
            else:
                front_part += "\\draw[%s] %s -- %s;\n" % ('edge',dict_drawing[index1][2],dict_drawing[index2][2])

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
        tikz_pic += '\tfacet/.style={fill=%s,fill opacity=%f},\n' % (facet_color,opacity)
        tikz_pic += '\tvertex/.style={inner sep=1pt,circle,draw=%s!25!black,' % vertex_color
        tikz_pic += 'fill=%s!75!black,thick,anchor=base}]\n%%\n%%\n' % vertex_color

        # Draws the axes if True
        if axis:
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
            if not v in front_vertices and v in dict_drawing:
                tikz_pic += dict_drawing[v][0]

        # Draw the facets in the front by going in cycles for every facet.
        tikz_pic += '%%\n%%\n%% Drawing the facets\n%%\n'
        for index_facet in front_facets:
            cyclic_vert = cyclic_sort_vertices_2d(list(facets[index_facet].incident()))
            cyclic_indices = [vertices.index(v) for v in cyclic_vert]
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
