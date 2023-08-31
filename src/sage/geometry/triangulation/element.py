"""
A triangulation

In Sage, the
:class:`~sage.geometry.triangulation.point_configuration.PointConfiguration`
and :class:`Triangulation` satisfy a parent/element relationship. In
particular, each triangulation refers back to its point
configuration. If you want to triangulate a point configuration, you
should construct a point configuration first and then use one of its
methods to triangulate it according to your requirements. You should
never have to construct a :class:`Triangulation` object directly.

EXAMPLES:

First, we select the internal implementation for enumerating
triangulations::

    sage: PointConfiguration.set_engine('internal')   # to make doctests independent of TOPCOM

Here is a simple example of how to triangulate a point configuration::

    sage: p = [[0,-1,-1],[0,0,1],[0,1,0], [1,-1,-1],[1,0,1],[1,1,0]]
    sage: points = PointConfiguration(p)
    sage: triang = points.triangulate();  triang
    (<0,1,2,5>, <0,1,3,5>, <1,3,4,5>)
    sage: triang.plot(axes=False)                                                       # needs sage.plot
    Graphics3d Object

See :mod:`sage.geometry.triangulation.point_configuration` for more details.
"""

#*****************************************************************************
#       Copyright (C) 2010 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.richcmp import richcmp
from sage.structure.element import Element
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.modules.free_module_element import vector
from sage.misc.cachefunc import cached_method
from sage.sets.set import Set
from sage.graphs.graph import Graph


########################################################################
def triangulation_render_2d(triangulation, **kwds):
    r"""
    Return a graphical representation of a 2-d triangulation.

    INPUT:

    - ``triangulation`` -- a :class:`Triangulation`.

    - ``**kwds`` -- keywords that are passed on to the graphics primitives.

    OUTPUT:

    A 2-d graphics object.

    EXAMPLES::

        sage: points = PointConfiguration([[0,0],[0,1],[1,0],[1,1],[-1,-1]])
        sage: triang = points.triangulate()
        sage: triang.plot(axes=False, aspect_ratio=1)   # indirect doctest              # needs sage.plot
        Graphics object consisting of 12 graphics primitives
    """
    from sage.plot.all import point2d, line2d, polygon2d
    points = [point.reduced_affine() for point in triangulation.point_configuration()]
    coord = [ [p[0], p[1]] for p in points ]
    plot_points = sum([ point2d(p,
                                zorder=2, pointsize=10, **kwds)
                        for p in coord ])

    tmp_lines = []
    for t in triangulation:
        if len(t)>=2:
            tmp_lines.append([t[0], t[1]])
        if len(t)>=3:
            tmp_lines.append([t[0], t[2]])
            tmp_lines.append([t[1], t[2]])
    all_lines = []
    interior_lines = []
    for l in tmp_lines:
        if l not in all_lines:
            all_lines.append(l)
        else:
            interior_lines.append(l)
    exterior_lines = [l for l in all_lines if l not in interior_lines]

    plot_interior_lines = sum([ line2d([ coord[l[0]], coord[l[1]] ],
                                       zorder=1, rgbcolor=(0,1,0), **kwds)
                                for l in interior_lines ])
    plot_exterior_lines = sum([ line2d([ coord[l[0]], coord[l[1]] ],
                                       zorder=1, rgbcolor=(0,0,1), **kwds)
                                for l in exterior_lines ])

    plot_triangs = sum([ polygon2d([coord[t[0]], coord[t[1]], coord[t[2]]],
                                   zorder=0, rgbcolor=(0.8, 1, 0.8), **kwds)
                         for t in triangulation if len(t)>=3 ])

    return \
        plot_points + \
        plot_interior_lines + plot_exterior_lines + \
        plot_triangs


def triangulation_render_3d(triangulation, **kwds):
    r"""
    Return a graphical representation of a 3-d triangulation.

    INPUT:

    - ``triangulation`` -- a :class:`Triangulation`.

    - ``**kwds`` -- keywords that are  passed on to the graphics primitives.

    OUTPUT:

    A 3-d graphics object.

    EXAMPLES::

        sage: p = [[0,-1,-1],[0,0,1],[0,1,0], [1,-1,-1],[1,0,1],[1,1,0]]
        sage: points = PointConfiguration(p)
        sage: triang = points.triangulate()
        sage: triang.plot(axes=False)     # indirect doctest                            # needs sage.plot
        Graphics3d Object
    """
    from sage.plot.plot3d.all import point3d, line3d, polygon3d
    points = [ point.reduced_affine() for point in triangulation.point_configuration() ]
    coord = [ [p[0], p[1], p[2] ] for p in points ]
    plot_points = sum([ point3d(p, size=15,
                                **kwds)
                        for p in coord ])

    tmp_lines = []
    for t in triangulation:
        if len(t)>=2:
            tmp_lines.append([t[0], t[1]])
        if len(t)>=3:
            tmp_lines.append([t[0], t[2]])
            tmp_lines.append([t[1], t[2]])
        if len(t)>=4:
            tmp_lines.append([t[0], t[3]])
            tmp_lines.append([t[1], t[3]])
            tmp_lines.append([t[2], t[3]])
    all_lines = []
    interior_lines = []
    for l in tmp_lines:
        if l not in all_lines:
            all_lines.append(l)
        else:
            interior_lines.append(l)
    exterior_lines = [l for l in all_lines if l not in interior_lines]

    from sage.plot.plot3d.texture import Texture
    line_int = Texture(color='darkblue', ambient=1, diffuse=0)
    line_ext = Texture(color='green', ambient=1, diffuse=0)
    triang_int = Texture(opacity=0.3, specular=0, shininess=0, diffuse=0, ambient=1, color='yellow')
    triang_ext = Texture(opacity=0.6, specular=0, shininess=0, diffuse=0, ambient=1, color='green')

    plot_interior_lines = sum([ line3d([ coord[l[0]], coord[l[1]] ],
                                       thickness=2, texture=line_int, **kwds)
                                for l in interior_lines ])
    plot_exterior_lines = sum([ line3d([ coord[l[0]], coord[l[1]] ],
                                       thickness=3, texture=line_ext, **kwds)
                                for l in exterior_lines ])

    tmp_triangs = []
    for t in triangulation:
        if len(t)>=3:
            tmp_triangs.append([t[0], t[1], t[2]])
        if len(t)>=4:
            tmp_triangs.append([t[0], t[1], t[3]])
            tmp_triangs.append([t[0], t[2], t[3]])
            tmp_triangs.append([t[1], t[2], t[3]])
    all_triangs = []
    interior_triangs = []
    for l in tmp_triangs:
        if l not in all_triangs:
            all_triangs.append(l)
        else:
            interior_triangs.append(l)
    exterior_triangs = [l for l in all_triangs if l not in interior_triangs]

    plot_interior_triangs = \
        sum([polygon3d([coord[t[0]], coord[t[1]], coord[t[2]]],
                       texture=triang_int, **kwds)
             for t in interior_triangs])
    plot_exterior_triangs = \
        sum([polygon3d([coord[t[0]], coord[t[1]], coord[t[2]]],
                       texture=triang_ext, **kwds)
              for t in exterior_triangs])

    return plot_points + \
        plot_interior_lines + plot_exterior_lines + \
        plot_interior_triangs + plot_exterior_triangs


########################################################################
class Triangulation(Element):
    """
    A triangulation of a
    :class:`~sage.geometry.triangulation.point_configuration.PointConfiguration`.

    .. WARNING::

        You should never create :class:`Triangulation` objects
        manually. See
        :meth:`~sage.geometry.triangulation.point_configuration.PointConfiguration.triangulate`
        and
        :meth:`~sage.geometry.triangulation.point_configuration.PointConfiguration.triangulations`
        to triangulate point configurations.
    """
    def __init__(self, triangulation, parent, check=True):
        """
        The constructor of a ``Triangulation`` object.

        Note that an internal reference to the underlying ``PointConfiguration`` is
        kept.

        INPUT:

        - ``parent`` -- a
          :class:`~sage.geometry.triangulation.point_configuration.PointConfiguration`

        - ``triangulation`` -- an iterable of integers or an iterable of
          iterables (e.g. a list of lists), specifying the maximal simplices
          of the triangulation. In the first case, each integer specifies a simplex
          by the correspondence :meth:`PointConfiguration.simplex_to_int`. In the second
          case, a simplex is specified by listing the indices of the included points.

        - ``check`` -- boolean. Whether to perform checks that the
          triangulation is, indeed, a triangulation of the point
          configuration.

        NOTE:

        Passing ``check=False`` allows you to create triangulations of
        subsets of the points of the configuration, see
        :meth:`~sage.geometry.triangulation.point_configuration.PointConfiguration.bistellar_flips`.

        EXAMPLES::

            sage: p = [[0,1],[0,0],[1,0]]
            sage: points = PointConfiguration(p)
            sage: from sage.geometry.triangulation.point_configuration import Triangulation
            sage: Triangulation([(0,1,2)], points)
            (<0,1,2>)
            sage: Triangulation([1], points)
            (<0,1,2>)
        """
        Element.__init__(self, parent=parent)
        self._point_configuration = parent

        try:
            triangulation = tuple(sorted( tuple(sorted(t)) for t in triangulation))
        except TypeError:
            triangulation = tuple( self.point_configuration().int_to_simplex(i)
                                   for i in triangulation )
        assert not check or all( len(t)==self.point_configuration().dim()+1
                                 for t in triangulation)
        self._triangulation = triangulation

    def point_configuration(self):
        """
        Returns the point configuration underlying the triangulation.

        EXAMPLES::

            sage: pconfig = PointConfiguration([[0,0],[0,1],[1,0]])
            sage: pconfig
            A point configuration in affine 2-space over Integer Ring
            consisting of 3 points. The triangulations of this point
            configuration are assumed to be connected, not necessarily
            fine, not necessarily regular.
            sage: triangulation = pconfig.triangulate()
            sage: triangulation
            (<0,1,2>)
            sage: triangulation.point_configuration()
            A point configuration in affine 2-space over Integer Ring
            consisting of 3 points. The triangulations of this point
            configuration are assumed to be connected, not necessarily
            fine, not necessarily regular.
            sage: pconfig == triangulation.point_configuration()
            True
        """
        return self._point_configuration

    def _richcmp_(self, right, op):
        r"""
        Compare ``self`` and ``right``.

        INPUT:

        - ``right`` -- a triangulation

        TESTS::

            sage: pc = PointConfiguration([[0,0],[0,1],[1,0]])
            sage: t1 = pc.triangulate()
            sage: from sage.geometry.triangulation.point_configuration import Triangulation
            sage: t2 = Triangulation([[2,1,0]], pc)
            sage: t1 is t2
            False
            sage: t1 == t2    # indirect doctest
            True
            sage: t1 != Triangulation(((0,1),(1,2)), pc, check=False)
            True
        """
        return richcmp(self._triangulation, right._triangulation, op)

    def __iter__(self):
        """
        Iterate through the simplices of the triangulation.

        EXAMPLES::

            sage: PointConfiguration.set_engine('internal')   # to make doctests independent of TOPCOM
            sage: pc = PointConfiguration([[0,0],[0,1],[1,0],[1,1],[-1,-1]])
            sage: triangulation = pc.triangulate()
            sage: iter = triangulation.__iter__()
            sage: next(iter)
            (1, 3, 4)
            sage: next(iter)
            (2, 3, 4)
            sage: next(iter)
            Traceback (most recent call last):
            ...
            StopIteration
        """
        for p in self._triangulation:
            yield p

    def __getitem__(self, i):
        """
        Access the point indices of the i-th simplex of the triangulation.

        INPUT:

        - ``i`` -- integer. The index of a simplex.

        OUTPUT:

        A tuple of integers. The vertex indices of the i-th simplex.

        EXAMPLES::

            sage: PointConfiguration.set_engine('internal')   # to make doctests independent of TOPCOM
            sage: pc = PointConfiguration([[0,0],[0,1],[1,0],[1,1],[-1,-1]])
            sage: triangulation = pc.triangulate()
            sage: triangulation[1]
            (2, 3, 4)
        """
        return self._triangulation[i]

    def __len__(self):
        """
        Return the length of the triangulation.

        TESTS::

            sage: PointConfiguration.set_engine('internal')   # to make doctests independent of TOPCOM
            sage: pc = PointConfiguration([[0,0],[0,1],[1,0],[1,1],[-1,-1]])
            sage: triangulation = next(pc.triangulations())
            sage: triangulation.__len__()
            2
            sage: len(triangulation)    # equivalent
            2
        """
        return len(self._triangulation)

    def _repr_(self):
        r"""
        Return a string representation.

        TESTS::

            sage: PointConfiguration.set_engine('internal')   # to make doctests independent of TOPCOM
            sage: pc = PointConfiguration([[0,0],[0,1],[1,0],[1,1],[-1,-1],[2,2]])
            sage: t = pc.triangulations()
            sage: next(t)._repr_()
            '(<1,4,5>, <2,4,5>)'
        """
        #s = 'A triangulation'
        #s += ' in QQ^'+str(self.point_configuration().ambient_dim())
        #s += ' consisting of '+str(len(self))+' simplices.'
        s = '('
        s += ', '.join([ '<'+','.join(map(str,t))+'>' for t in self._triangulation])
        s += ')'
        return s

    def plot(self, **kwds):
        r"""
        Produce a graphical representation of the triangulation.

        EXAMPLES::

            sage: p = PointConfiguration([[0,0],[0,1],[1,0],[1,1],[-1,-1]])
            sage: triangulation = p.triangulate()
            sage: triangulation
            (<1,3,4>, <2,3,4>)
            sage: triangulation.plot(axes=False)                                        # needs sage.plot
            Graphics object consisting of 12 graphics primitives
        """
        dim = self.point_configuration().dim()

        if dim == 2:
            return triangulation_render_2d(self, **kwds)

        if dim == 3:
            return triangulation_render_3d(self, **kwds)

        raise NotImplementedError('Plotting '+str(dim)+'-dimensional triangulations not implemented!')

    def gkz_phi(self):
        r"""
        Calculate the GKZ phi vector of the triangulation.

        The phi vector is a vector of length equals to the number of
        points in the point configuration. For a fixed triangulation
        `T`, the entry corresponding to the `i`-th point `p_i` is

        .. MATH::

            \phi_T(p_i) = \sum_{t\in T, t\owns p_i} Vol(t)

        that is, the total volume of all simplices containing `p_i`.
        See also [GKZ1994]_ page 220 equation 1.4.

        OUTPUT:

        The phi vector of self.

        EXAMPLES::

            sage: p = PointConfiguration([[0,0],[1,0],[2,1],[1,2],[0,1]])
            sage: p.triangulate().gkz_phi()
            (3, 1, 5, 2, 4)
            sage: p.lexicographic_triangulation().gkz_phi()
            (1, 3, 4, 2, 5)
        """
        vec = vector(ZZ, self.point_configuration().n_points())
        for simplex in self:
            vol = self.point_configuration().volume(simplex)
            for i in simplex:
                vec[i] = vec[i] + vol
        return vec

    def enumerate_simplices(self):
        r"""
        Return the enumerated simplices.

        OUTPUT:

        A tuple of integers that uniquely specifies the triangulation.

        EXAMPLES::

            sage: pc = PointConfiguration(matrix([
            ....:    [ 0, 0, 0, 0, 0, 2, 4,-1, 1, 1, 0, 0, 1, 0],
            ....:    [ 0, 0, 0, 1, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0],
            ....:    [ 0, 2, 0, 0, 0, 0,-1, 0, 1, 0, 1, 0, 0, 1],
            ....:    [ 0, 1, 1, 0, 0, 1, 0,-2, 1, 0, 0,-1, 1, 1],
            ....:    [ 0, 0, 0, 0, 1, 0,-1, 0, 0, 0, 0, 0, 0, 0]
            ....: ]).columns())
            sage: triangulation = pc.lexicographic_triangulation()
            sage: triangulation.enumerate_simplices()
            (1678, 1688, 1769, 1779, 1895, 1905, 2112, 2143, 2234, 2360, 2555, 2580,
             2610, 2626, 2650, 2652, 2654, 2661, 2663, 2667, 2685, 2755, 2757, 2759,
             2766, 2768, 2772, 2811, 2881, 2883, 2885, 2892, 2894, 2898)

        You can recreate the triangulation from this list by passing
        it to the constructor::

            sage: from sage.geometry.triangulation.point_configuration import Triangulation
            sage: Triangulation([1678, 1688, 1769, 1779, 1895, 1905, 2112, 2143,
            ....:  2234, 2360, 2555, 2580, 2610, 2626, 2650, 2652, 2654, 2661, 2663,
            ....:  2667, 2685, 2755, 2757, 2759, 2766, 2768, 2772, 2811, 2881, 2883,
            ....:  2885, 2892, 2894, 2898], pc)
            (<1,3,4,7,10,13>, <1,3,4,8,10,13>, <1,3,6,7,10,13>, <1,3,6,8,10,13>,
             <1,4,6,7,10,13>, <1,4,6,8,10,13>, <2,3,4,6,7,12>, <2,3,4,7,12,13>,
             <2,3,6,7,12,13>, <2,4,6,7,12,13>, <3,4,5,6,9,12>, <3,4,5,8,9,12>,
             <3,4,6,7,11,12>, <3,4,6,9,11,12>, <3,4,7,10,11,13>, <3,4,7,11,12,13>,
             <3,4,8,9,10,12>, <3,4,8,10,12,13>, <3,4,9,10,11,12>, <3,4,10,11,12,13>,
             <3,5,6,8,9,12>, <3,6,7,10,11,13>, <3,6,7,11,12,13>, <3,6,8,9,10,12>,
             <3,6,8,10,12,13>, <3,6,9,10,11,12>, <3,6,10,11,12,13>, <4,5,6,8,9,12>,
             <4,6,7,10,11,13>, <4,6,7,11,12,13>, <4,6,8,9,10,12>, <4,6,8,10,12,13>,
             <4,6,9,10,11,12>, <4,6,10,11,12,13>)
        """
        pc = self._point_configuration
        return tuple( pc.simplex_to_int(t) for t in self )

    def fan(self, origin=None):
        r"""
        Construct the fan of cones over the simplices of the triangulation.

        INPUT:

        - ``origin`` -- ``None`` (default) or coordinates of a
          point. The common apex of all cones of the fan. If ``None``,
          the triangulation must be a star triangulation and the
          distinguished central point is used as the origin.

        OUTPUT:

        A :class:`~sage.geometry.fan.RationalPolyhedralFan`. The
        coordinates of the points are shifted so that the apex of the
        fan is the origin of the coordinate system.

        .. note:: If the set of cones over the simplices is not a fan, a
            suitable exception is raised.

        EXAMPLES::

            sage: pc = PointConfiguration([(0,0), (1,0), (0,1), (-1,-1)], star=0, fine=True)
            sage: triangulation = pc.triangulate()
            sage: fan = triangulation.fan(); fan
            Rational polyhedral fan in 2-d lattice N
            sage: fan.is_equivalent(toric_varieties.P2().fan())                         # needs palp sage.graphs
            True

        Toric diagrams (the `\ZZ_5` hyperconifold)::

            sage: vertices=[(0, 1, 0), (0, 3, 1), (0, 2, 3), (0, 0, 2)]
            sage: interior=[(0, 1, 1), (0, 1, 2), (0, 2, 1), (0, 2, 2)]
            sage: points = vertices + interior
            sage: pc = PointConfiguration(points, fine=True)
            sage: triangulation = pc.triangulate()
            sage: fan = triangulation.fan((-1,0,0)); fan
            Rational polyhedral fan in 3-d lattice N
            sage: fan.rays()
            N(1, 1, 0),
            N(1, 3, 1),
            N(1, 2, 3),
            N(1, 0, 2),
            N(1, 1, 1),
            N(1, 1, 2),
            N(1, 2, 1),
            N(1, 2, 2)
            in 3-d lattice N
        """
        from sage.geometry.fan import Fan
        if origin is None:
            origin = self.point_configuration().star_center()
        R = self.base_ring()
        origin = vector(R, origin)
        points = self.point_configuration().points()
        return Fan(self, (vector(R, p) - origin for p in points))

    @cached_method
    def simplicial_complex(self):
        r"""
        Return ``self`` as an (abstract) simplicial complex.

        OUTPUT:

        A :class:`~sage.topology.simplicial_complex.SimplicialComplex`.

        EXAMPLES::

            sage: p = polytopes.cuboctahedron()
            sage: sc = p.triangulate(engine='internal').simplicial_complex(); sc        # needs sage.graphs
            Simplicial complex with 12 vertices and 16 facets

        Any convex set is contractable, so its reduced homology groups vanish::

            sage: sc.homology()                                                         # needs sage.graphs
            {0: 0, 1: 0, 2: 0, 3: 0}
        """
        from sage.topology.simplicial_complex import SimplicialComplex
        return SimplicialComplex(self)

    @cached_method
    def _boundary_simplex_dictionary(self):
        """
        Return facets and the simplices they bound.

        TESTS::

            sage: triangulation = polytopes.hypercube(2).triangulate(engine='internal')
            sage: triangulation._boundary_simplex_dictionary()
            {(0, 1): ((0, 1, 3),),
             (0, 3): ((0, 1, 3),),
             (1, 2): ((1, 2, 3),),
             (1, 3): ((0, 1, 3), (1, 2, 3)),
             (2, 3): ((1, 2, 3),)}

            sage: triangulation = polytopes.cube().triangulate(engine='internal')
            sage: triangulation._boundary_simplex_dictionary()
            {(0, 1, 2): ((0, 1, 2, 7),),
             (0, 1, 5): ((0, 1, 5, 7),),
             (0, 1, 7): ((0, 1, 2, 7), (0, 1, 5, 7)),
             (0, 2, 3): ((0, 2, 3, 7),),
             (0, 2, 7): ((0, 1, 2, 7), (0, 2, 3, 7)),
             (0, 3, 4): ((0, 3, 4, 7),),
             (0, 3, 7): ((0, 2, 3, 7), (0, 3, 4, 7)),
             (0, 4, 5): ((0, 4, 5, 7),),
             (0, 4, 7): ((0, 3, 4, 7), (0, 4, 5, 7)),
             (0, 5, 7): ((0, 1, 5, 7), (0, 4, 5, 7)),
             (1, 2, 7): ((0, 1, 2, 7),),
             (1, 5, 6): ((1, 5, 6, 7),),
             (1, 5, 7): ((0, 1, 5, 7), (1, 5, 6, 7)),
             (1, 6, 7): ((1, 5, 6, 7),),
             (2, 3, 7): ((0, 2, 3, 7),),
             (3, 4, 7): ((0, 3, 4, 7),),
             (4, 5, 7): ((0, 4, 5, 7),),
             (5, 6, 7): ((1, 5, 6, 7),)}
        """
        result = dict()
        for simplex in self:
            for i in range(len(simplex)):
                facet = simplex[:i] + simplex[i+1:]
                result[facet] = result.get(facet, tuple()) + (simplex,)
        return result

    @cached_method
    def boundary(self):
        """
        Return the boundary of the triangulation.

        OUTPUT:

        The outward-facing boundary simplices (of dimension `d-1`) of
        the `d`-dimensional triangulation as a set. Each boundary is
        returned by a tuple of point indices.

        EXAMPLES::

            sage: triangulation = polytopes.cube().triangulate(engine='internal')
            sage: triangulation
            (<0,1,2,7>, <0,1,5,7>, <0,2,3,7>, <0,3,4,7>, <0,4,5,7>, <1,5,6,7>)
            sage: triangulation.boundary()
            frozenset({(0, 1, 2),
                       (0, 1, 5),
                       (0, 2, 3),
                       (0, 3, 4),
                       (0, 4, 5),
                       (1, 2, 7),
                       (1, 5, 6),
                       (1, 6, 7),
                       (2, 3, 7),
                       (3, 4, 7),
                       (4, 5, 7),
                       (5, 6, 7)})
            sage: triangulation.interior_facets()
            frozenset({(0, 1, 7), (0, 2, 7), (0, 3, 7), (0, 4, 7), (0, 5, 7), (1, 5, 7)})
        """
        return frozenset(facet for facet, bounded_simplices
                         in self._boundary_simplex_dictionary().items()
                         if len(bounded_simplices) == 1)

    @cached_method
    def boundary_simplicial_complex(self):
        r"""
        Return the boundary of ``self`` as an (abstract) simplicial complex.

        OUTPUT:

        A :class:`~sage.topology.simplicial_complex.SimplicialComplex`.

        EXAMPLES::

            sage: p = polytopes.cuboctahedron()
            sage: triangulation = p.triangulate(engine='internal')
            sage: bd_sc = triangulation.boundary_simplicial_complex(); bd_sc            # needs sage.graphs
            Simplicial complex with 12 vertices and 20 facets

        The boundary of every convex set is a topological sphere, so it has
        spherical homology::

            sage: bd_sc.homology()                                                      # needs sage.graphs
            {0: 0, 1: 0, 2: Z}

        It is a subcomplex of ``self`` as a :meth:`simplicial_complex`::

            sage: sc = triangulation.simplicial_complex()                               # needs sage.graphs
            sage: all(f in sc for f in bd_sc.maximal_faces())                           # needs sage.graphs
            True
        """
        from sage.topology.simplicial_complex import SimplicialComplex
        return SimplicialComplex(self.boundary(), maximality_check=False)

    @cached_method
    def interior_facets(self):
        """
        Return the interior facets of the triangulation.

        OUTPUT:

        The inward-facing boundary simplices (of dimension `d-1`) of
        the `d`-dimensional triangulation as a set. Each boundary is
        returned by a tuple of point indices.

        EXAMPLES::

            sage: triangulation = polytopes.cube().triangulate(engine='internal')
            sage: triangulation
            (<0,1,2,7>, <0,1,5,7>, <0,2,3,7>, <0,3,4,7>, <0,4,5,7>, <1,5,6,7>)
            sage: triangulation.boundary()
            frozenset({(0, 1, 2),
                       (0, 1, 5),
                       (0, 2, 3),
                       (0, 3, 4),
                       (0, 4, 5),
                       (1, 2, 7),
                       (1, 5, 6),
                       (1, 6, 7),
                       (2, 3, 7),
                       (3, 4, 7),
                       (4, 5, 7),
                       (5, 6, 7)})
            sage: triangulation.interior_facets()
            frozenset({(0, 1, 7), (0, 2, 7), (0, 3, 7), (0, 4, 7), (0, 5, 7), (1, 5, 7)})
        """
        return frozenset(facet for facet, bounded_simplices
                         in self._boundary_simplex_dictionary().items()
                         if len(bounded_simplices) == 2)

    def polyhedral_complex(self, **kwds):
        """
        Return ``self`` as a :class:`~sage.geometry.polyhedral_complex.PolyhedralComplex`.

        OUTPUT:

        A :class:`~sage.geometry.polyhedral_complex.PolyhedralComplex` whose maximal cells
        are the simplices of the triangulation.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: pc = PointConfiguration(P.vertices())
            sage: T = pc.placing_triangulation(); T
            (<0,1,2,7>, <0,1,5,7>, <0,2,3,7>, <0,3,4,7>, <0,4,5,7>, <1,5,6,7>)
            sage: C = T.polyhedral_complex(); C                                         # needs sage.graphs
            Polyhedral complex with 6 maximal cells
            sage: [P.vertices_list() for P in C.maximal_cells_sorted()]                 # needs sage.graphs
            [[[-1, -1, -1], [-1, -1, 1], [-1, 1, 1], [1, -1, -1]],
             [[-1, -1, -1], [-1, 1, -1], [-1, 1, 1], [1, 1, -1]],
             [[-1, -1, -1], [-1, 1, 1], [1, -1, -1], [1, 1, -1]],
             [[-1, -1, 1], [-1, 1, 1], [1, -1, -1], [1, -1, 1]],
             [[-1, 1, 1], [1, -1, -1], [1, -1, 1], [1, 1, 1]],
             [[-1, 1, 1], [1, -1, -1], [1, 1, -1], [1, 1, 1]]]
        """
        from sage.geometry.polyhedral_complex import PolyhedralComplex
        from sage.geometry.polyhedron.constructor import Polyhedron
        ambient_dim = self.point_configuration().ambient_dim()
        points = self.point_configuration().points()
        return PolyhedralComplex([Polyhedron(vertices=[points[i] for i in simplex])
                                  for simplex in self],
                                 ambient_dim=ambient_dim,
                                 maximality_check=False,
                                 face_to_face_check=False,
                                 **kwds)

    def boundary_polyhedral_complex(self, **kwds):
        r"""
        Return the boundary of ``self`` as a :class:`~sage.geometry.polyhedral_complex.PolyhedralComplex`.

        OUTPUT:

        A :class:`~sage.geometry.polyhedral_complex.PolyhedralComplex` whose maximal cells
        are the simplices of the boundary of ``self``.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: pc = PointConfiguration(P.vertices())
            sage: T = pc.placing_triangulation(); T
            (<0,1,2,7>, <0,1,5,7>, <0,2,3,7>, <0,3,4,7>, <0,4,5,7>, <1,5,6,7>)
            sage: bd_C = T.boundary_polyhedral_complex(); bd_C                          # needs sage.graphs
            Polyhedral complex with 12 maximal cells
            sage: [P.vertices_list() for P in bd_C.maximal_cells_sorted()]              # needs sage.graphs
            [[[-1, -1, -1], [-1, -1, 1], [-1, 1, 1]],
             [[-1, -1, -1], [-1, -1, 1], [1, -1, -1]],
             [[-1, -1, -1], [-1, 1, -1], [-1, 1, 1]],
             [[-1, -1, -1], [-1, 1, -1], [1, 1, -1]],
             [[-1, -1, -1], [1, -1, -1], [1, 1, -1]],
             [[-1, -1, 1], [-1, 1, 1], [1, -1, 1]],
             [[-1, -1, 1], [1, -1, -1], [1, -1, 1]],
             [[-1, 1, -1], [-1, 1, 1], [1, 1, -1]],
             [[-1, 1, 1], [1, -1, 1], [1, 1, 1]],
             [[-1, 1, 1], [1, 1, -1], [1, 1, 1]],
             [[1, -1, -1], [1, -1, 1], [1, 1, 1]],
             [[1, -1, -1], [1, 1, -1], [1, 1, 1]]]

        It is a subcomplex of ``self`` as a :meth:`polyhedral_complex`::

            sage: C = T.polyhedral_complex()                                            # needs sage.graphs
            sage: bd_C.is_subcomplex(C)                                                 # needs sage.graphs
            True
        """
        from sage.geometry.polyhedral_complex import PolyhedralComplex
        from sage.geometry.polyhedron.constructor import Polyhedron
        ambient_dim = self.point_configuration().ambient_dim()
        points = self.point_configuration().points()
        return PolyhedralComplex([Polyhedron(vertices=[points[i] for i in simplex])
                                  for simplex in self.boundary()],
                                 ambient_dim=ambient_dim,
                                 maximality_check=False,
                                 face_to_face_check=False,
                                 **kwds)

    @cached_method
    def normal_cone(self):
        r"""
        Return the (closure of the) normal cone of the triangulation.

        Recall that a regular triangulation is one that equals the
        "crease lines" of a convex piecewise-linear function. This
        support function is not unique, for example, you can scale it
        by a positive constant. The set of all piecewise-linear
        functions with fixed creases forms an open cone. This cone can
        be interpreted as the cone of normal vectors at a point of the
        secondary polytope, which is why we call it normal cone. See
        [GKZ1994]_ Section 7.1 for details.

        OUTPUT:

        The closure of the normal cone. The `i`-th entry equals the
        value of the piecewise-linear function at the `i`-th point of
        the configuration.

        For an irregular triangulation, the normal cone is empty. In
        this case, a single point (the origin) is returned.

        EXAMPLES::

            sage: triangulation = polytopes.hypercube(2).triangulate(engine='internal')
            sage: triangulation
            (<0,1,3>, <1,2,3>)
            sage: N = triangulation.normal_cone();  N
            4-d cone in 4-d lattice
            sage: N.rays()
            ( 0,  0,  0, -1),
            ( 0,  0,  1,  1),
            ( 0,  0, -1, -1),
            ( 1,  0,  0,  1),
            (-1,  0,  0, -1),
            ( 0,  1,  0, -1),
            ( 0, -1,  0,  1)
            in Ambient free module of rank 4
            over the principal ideal domain Integer Ring
            sage: N.dual().rays()
            (1, -1, 1, -1)
            in Ambient free module of rank 4
            over the principal ideal domain Integer Ring

        TESTS::

            sage: polytopes.simplex(2).triangulate().normal_cone()
            3-d cone in 3-d lattice
            sage: _.dual().is_trivial()
            True
        """
        if not self.point_configuration().base_ring().is_subring(QQ):
            raise NotImplementedError('Only base rings ZZ and QQ are supported')
        from ppl import Constraint_System, Linear_Expression, C_Polyhedron
        from sage.matrix.constructor import matrix
        from sage.arith.functions import lcm
        pc = self.point_configuration()
        cs = Constraint_System()
        for facet in self.interior_facets():
            s0, s1 = self._boundary_simplex_dictionary()[facet]
            p = set(s0).difference(facet).pop()
            q = set(s1).difference(facet).pop()
            origin = pc.point(p).reduced_affine_vector()
            base_indices = [i for i in s0 if i != p]
            base = matrix([ pc.point(i).reduced_affine_vector()-origin for i in base_indices ])
            sol = base.solve_left( pc.point(q).reduced_affine_vector()-origin )
            relation = [0]*pc.n_points()
            relation[p] = sum(sol)-1
            relation[q] = 1
            for i, base_i in enumerate(base_indices):
                relation[base_i] = -sol[i]
            rel_denom = lcm([QQ(r).denominator() for r in relation])
            relation = [ ZZ(r*rel_denom) for r in relation ]
            ex = Linear_Expression(relation,0)
            cs.insert(ex >= 0)
        from sage.modules.free_module import FreeModule
        ambient = FreeModule(ZZ, self.point_configuration().n_points())
        if cs.empty():
            cone = C_Polyhedron(ambient.dimension(), 'universe')
        else:
            cone = C_Polyhedron(cs)
        from sage.geometry.cone import _Cone_from_PPL
        return _Cone_from_PPL(cone, lattice=ambient)

    def adjacency_graph(self):
        """
        Return a graph showing which simplices are adjacent in the triangulation.

        OUTPUT:

        A graph consisting of vertices referring to the simplices in the
        triangulation, and edges showing which simplices are adjacent to each
        other.

        .. SEEALSO::

            * To obtain the triangulation's 1-skeleton, use
              :meth:`SimplicialComplex.graph` through
              ``MyTriangulation.simplicial_complex().graph()``.

        AUTHORS:

        * Stephen Farley (2013-08-10): initial version

        EXAMPLES::

            sage: p = PointConfiguration([[1,0,0], [0,1,0], [0,0,1], [-1,0,1],
            ....:                         [1,0,-1], [-1,0,0], [0,-1,0], [0,0,-1]])
            sage: t = p.triangulate()
            sage: t.adjacency_graph()                                                   # needs sage.graphs
            Graph on 8 vertices

        """
        vertices = [Set(_) for _ in list(self)]
        return Graph([vertices,
                  lambda x,y: len(x-y)==1])
