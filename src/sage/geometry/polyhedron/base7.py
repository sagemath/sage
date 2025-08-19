r"""
Base class for polyhedra: Methods for triangulation and volume computation
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

from sage.cpython.string import bytes_to_str
from sage.misc.cachefunc import cached_method
from sage.modules.free_module_element import vector
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from .base6 import Polyhedron_base6


class Polyhedron_base7(Polyhedron_base6):
    r"""
    Methods related to triangulation and volume.

    TESTS::

        sage: # needs sage.combinat
        sage: from sage.geometry.polyhedron.base7 import Polyhedron_base7
        sage: P = polytopes.associahedron(['A', 3])
        sage: Polyhedron_base7.centroid(P)
        (81/632, 36/79, 81/632)
        sage: Polyhedron_base7.triangulate(P)
        (<0,1,2,13>, <0,1,7,13>, <0,2,5,13>, <0,6,7,12>, <0,6,8,13>,
         <0,6,12,13>, <0,7,12,13>, <1,2,7,12>, <1,2,12,13>, <1,7,12,13>,
         <2,3,7,12>, <2,3,12,13>, <3,4,7,12>, <3,11,12,13>, <6,8,9,12>,
         <6,8,12,13>, <6,9,10,12>, <8,9,12,13>)
        sage: Polyhedron_base7.volume(P, measure='induced')
        79/3
    """
    @cached_method(do_pickle=True)
    def centroid(self, engine='auto', **kwds):
        r"""
        Return the center of the mass of the polytope.

        The mass is taken with respect to the induced Lebesgue measure,
        see :meth:`volume`.

        If the polyhedron is not compact, a :exc:`NotImplementedError` is
        raised.

        INPUT:

        - ``engine`` -- either 'auto' (default), 'internal',
          'TOPCOM', or 'normaliz'.  The 'internal' and 'TOPCOM' instruct
          this package to always use its own triangulation algorithms
          or TOPCOM's algorithms, respectively. By default ('auto'),
          TOPCOM is used if it is available and internal routines otherwise.

        - ``**kwds`` -- keyword arguments that are passed to the
          triangulation engine (see :meth:`triangulate`)

        OUTPUT: the centroid as vector

        ALGORITHM:

        We triangulate the polytope and find the barycenter of the simplices.
        We add the individual barycenters weighted by the fraction of the total
        mass.

        EXAMPLES::

            sage: P = polytopes.hypercube(2).pyramid()
            sage: P.centroid()
            (1/4, 0, 0)

            sage: P = polytopes.associahedron(['A', 2])                                 # needs sage.combinat
            sage: P.centroid()                                                          # needs sage.combinat
            (2/21, 2/21)

            sage: P = polytopes.permutahedron(4, backend='normaliz')    # optional - pynormaliz
            sage: P.centroid()                                          # optional - pynormaliz
            (5/2, 5/2, 5/2, 5/2)

        The method is not implemented for unbounded polyhedra::

            sage: P = Polyhedron(vertices=[(0, 0)], rays=[(1, 0), (0, 1)])
            sage: P.centroid()
            Traceback (most recent call last):
            ...
            NotImplementedError: the polyhedron is not compact

        The centroid of an empty polyhedron is not defined::

            sage: Polyhedron().centroid()
            Traceback (most recent call last):
            ...
            ZeroDivisionError: rational division by zero

        TESTS::

            sage: Polyhedron(vertices=[[0,1]]).centroid()
            (0, 1)
        """
        if not self.is_compact():
            raise NotImplementedError("the polyhedron is not compact")
        if self.n_vertices() == self.dim() + 1:
            # The centroid of a simplex is its center.
            return self.center()

        triangulation = self.triangulate(engine=engine, **kwds)

        if self.ambient_dim() == self.dim():
            pc = triangulation.point_configuration()
        else:
            from sage.geometry.triangulation.point_configuration import PointConfiguration
            A, b = self.affine_hull_projection(as_affine_map=True, orthogonal=True, orthonormal=True, extend=True)
            pc = PointConfiguration(A(v.vector()) for v in self.Vrep_generator())

        barycenters = [sum(self.Vrepresentation(i).vector() for i in simplex)/(self.dim() + 1) for simplex in triangulation]
        volumes = [pc.volume(simplex) for simplex in triangulation]

        centroid = sum(volumes[i]*barycenters[i] for i in range(len(volumes)))/sum(volumes)
        if self.ambient_dim() != self.dim():
            # By the affine hull projection, the centroid has base ring ``AA``,
            # we try return the centroid in a reasonable ring.
            try:
                return centroid.change_ring(self.base_ring().fraction_field())
            except ValueError:
                pass
        return centroid

    def _triangulate_normaliz(self):
        r"""
        Give a triangulation of the polyhedron using normaliz.

        OUTPUT:

        A tuple of pairs ``(simplex,simplex_volume)`` used in the
        triangulation.

        .. NOTE::

            This function depends on Normaliz (i.e. the ``pynormaliz`` optional
            package). See the Normaliz documentation for further details.

        TESTS::

            sage: K = Polyhedron(vertices=[[1,1]], rays=[[1,0],[1,2]])
            sage: K._triangulate_normaliz()
            Traceback (most recent call last):
            ...
            TypeError: the polyhedron's backend should be 'normaliz'
        """
        raise TypeError("the polyhedron's backend should be 'normaliz'")

    def triangulate(self, engine='auto', connected=True, fine=False, regular=None, star=None):
        r"""
        Return a triangulation of the polytope.

        INPUT:

        - ``engine`` -- either 'auto' (default), 'internal',
          'TOPCOM', or 'normaliz'.  The 'internal' and 'TOPCOM' instruct
          this package to always use its own triangulation algorithms
          or TOPCOM's algorithms, respectively. By default ('auto'),
          TOPCOM is used if it is available and internal routines otherwise.

        The remaining keyword parameters are passed through to the
        :class:`~sage.geometry.triangulation.point_configuration.PointConfiguration`
        constructor:

        - ``connected`` -- boolean (default: ``True``); whether the
          triangulations should be connected to the regular
          triangulations via bistellar flips. These are much easier to
          compute than all triangulations.

        - ``fine`` -- boolean (default: ``False``); whether the
          triangulations must be fine, that is, make use of all points
          of the configuration

        - ``regular`` -- boolean or ``None`` (default:
          ``None``); whether the triangulations must be regular. A
          regular triangulation is one that is induced by a
          piecewise-linear convex support function. In other words,
          the shadows of the faces of a polyhedron in one higher
          dimension.

          * ``True``: Only regular triangulations.

          * ``False``: Only non-regular triangulations.

          * ``None`` (default): Both kinds of triangulation.

        - ``star`` -- either ``None`` (default) or a point. Whether
          the triangulations must be star. A triangulation is star if
          all maximal simplices contain a common point. The central
          point can be specified by its index (an integer) in the
          given points or by its coordinates (anything iterable.)

        OUTPUT:

        A triangulation of the convex hull of the vertices as a
        :class:`~sage.geometry.triangulation.element.Triangulation`. The
        indices in the triangulation correspond to the
        :meth:`~sage.geometry.polyhedron.base0.Polyhedron_base0.Vrepresentation` objects.

        EXAMPLES::

            sage: cube = polytopes.hypercube(3)
            sage: triangulation = cube.triangulate(
            ....:    engine='internal')  # to make doctest independent of TOPCOM
            sage: triangulation
            (<0,1,2,7>, <0,1,5,7>, <0,2,3,7>, <0,3,4,7>, <0,4,5,7>, <1,5,6,7>)
            sage: simplex_indices = triangulation[0]; simplex_indices
            (0, 1, 2, 7)
            sage: simplex_vertices = [cube.Vrepresentation(i) for i in simplex_indices]
            sage: simplex_vertices
            [A vertex at (1, -1, -1),
             A vertex at (1, 1, -1),
             A vertex at (1, 1, 1),
             A vertex at (-1, 1, 1)]
            sage: Polyhedron(simplex_vertices)
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 4 vertices

        It is possible to use ``'normaliz'`` as an engine. For this, the
        polyhedron should have the backend set to normaliz::

            sage: P = Polyhedron(vertices=[[0,0,1], [1,0,1],            # optional - pynormaliz
            ....:                          [0,1,1], [1,1,1]],
            ....:                backend='normaliz')
            sage: P.triangulate(engine='normaliz')                      # optional - pynormaliz
            (<0,1,2>, <1,2,3>)

            sage: P = Polyhedron(vertices=[[0,0,1], [1,0,1],
            ....:                          [0,1,1], [1,1,1]])
            sage: P.triangulate(engine='normaliz')
            Traceback (most recent call last):
            ...
            TypeError: the polyhedron's backend should be 'normaliz'

        The normaliz engine can triangulate pointed cones::

            sage: # optional - pynormaliz
            sage: C1 = Polyhedron(rays=[[0,0,1], [1,0,1],
            ....:                       [0,1,1], [1,1,1]],
            ....:                 backend='normaliz')
            sage: C1.triangulate(engine='normaliz')
            (<0,1,2>, <1,2,3>)
            sage: C2 = Polyhedron(rays=[[1,0,1], [0,0,1],
            ....:                       [0,1,1], [1,1,10/9]],
            ....:                 backend='normaliz')
            sage: C2.triangulate(engine='normaliz')
            (<0,1,2>, <1,2,3>)

        They can also be affine cones::

            sage: K = Polyhedron(vertices=[[1,1,1]],                    # optional - pynormaliz
            ....:                rays=[[1,0,0], [0,1,0], [1,1,-1], [1,1,1]],
            ....:                backend='normaliz')
            sage: K.triangulate(engine='normaliz')                      # optional - pynormaliz
            (<0,1,2>, <0,1,3>)
        """
        if self.lines():
            raise NotImplementedError('triangulation of polyhedra with lines is not supported')
        if len(self.vertices_list()) >= 2 and self.rays_list():
            raise NotImplementedError('triangulation of non-compact polyhedra that are not cones is not supported')
        if not self.is_compact() and engine != 'normaliz':
            raise NotImplementedError("triangulation of pointed polyhedra requires 'normaliz'")
        from sage.geometry.triangulation.point_configuration import PointConfiguration
        if self.is_compact():
            pc = PointConfiguration((v.vector() for v in self.vertex_generator()),
                                    connected=connected, fine=fine, regular=regular, star=star)
            # If the engine is not normaliz, we pass directly to the
            # PointConfiguration module.
            if engine != 'normaliz':
                pc.set_engine(engine)
                return pc.triangulate()
            else:
                return pc(self._triangulate_normaliz())
        else:  # From above, we have a pointed cone and the engine is normaliz
            try:
                pc = PointConfiguration((v.vector() for v in self.ray_generator()),
                                        connected=connected, fine=fine, regular=regular, star=star)
                return pc(self._triangulate_normaliz())
            except AssertionError:
                # PointConfiguration is not adapted to inhomogeneous cones
                # This is a hack. TODO: Implement the necessary things in
                # PointConfiguration to accept such cases.
                c = self.representative_point()
                normed_v = ((1/(r.vector()*c))*r.vector() for r in self.ray_generator())
                pc = PointConfiguration(normed_v, connected=connected, fine=fine, regular=regular, star=star)
                return pc(self._triangulate_normaliz())

    def _volume_lrs(self, verbose=False):
        """
        Compute the volume of a polytope using lrs.

        OUTPUT: the exact volume as a rational number

        EXAMPLES::

            sage: polytopes.hypercube(3)._volume_lrs()                  # optional - lrslib
            8
            sage: (polytopes.hypercube(3)*2)._volume_lrs()              # optional - lrslib
            64
            sage: polytopes.twenty_four_cell()._volume_lrs()            # optional - lrslib
            2

        REFERENCES:

        - David Avis's lrs program.
        """
        from sage.features.lrs import Lrs
        Lrs().require()

        from sage.misc.temporary_file import tmp_filename
        from subprocess import Popen, PIPE

        in_str = self.cdd_Vrepresentation()
        in_str += 'volume'
        in_filename = tmp_filename()
        in_file = open(in_filename, 'w')
        in_file.write(in_str)
        in_file.close()
        if verbose:
            print(in_str)

        lrs_procs = Popen([Lrs().absolute_filename(), in_filename],
                          stdin=PIPE, stdout=PIPE, stderr=PIPE)
        ans, err = lrs_procs.communicate()
        ans = bytes_to_str(ans)
        err = bytes_to_str(err)
        if verbose:
            print(ans)
        # FIXME: check err

        for a_line in ans.splitlines():
            if 'Volume=' in a_line:
                volume = a_line.split('Volume=')[1]
                volume = QQ(volume)
                return volume

        raise ValueError("lrs did not return a volume")

    def _volume_latte(self, verbose=False, algorithm='triangulate', **kwargs):
        """
        Compute the volume of a polytope using LattE integrale.

        INPUT:

        - ``arg`` -- a cdd or LattE description string

        - ``algorithm`` -- (default: ``'triangulate'``) the integration method;
          use 'triangulate' for polytope triangulation or 'cone-decompose' for
          tangent cone decomposition method

        - ``raw_output`` -- if ``True`` then return directly the output string
          from LattE

        - ``verbose`` -- if ``True`` then return directly verbose output from
          LattE

        - For all other options, consult the LattE manual.

        OUTPUT: a rational value, or a string if ``raw_output`` if set to ``True``

        .. NOTE::

            This function depends on LattE (i.e., the ``latte_int`` optional
            package). See the LattE documentation for further details.

        EXAMPLES::

            sage: # optional - latte_int
            sage: polytopes.hypercube(3)._volume_latte()
            8
            sage: (polytopes.hypercube(3)*2)._volume_latte()
            64
            sage: polytopes.twenty_four_cell()._volume_latte()
            2
            sage: polytopes.cuboctahedron()._volume_latte()
            20/3

        TESTS:

        Testing triangulate algorithm::

            sage: polytopes.cuboctahedron()._volume_latte(algorithm='triangulate')      # optional - latte_int
            20/3

        Testing cone decomposition algorithm::

            sage: polytopes.cuboctahedron()._volume_latte(algorithm='cone-decompose')   # optional - latte_int
            20/3

        Testing raw output::

            sage: polytopes.cuboctahedron()._volume_latte(raw_output=True)              # optional - latte_int
            '20/3'

        Testing inexact rings::

            sage: P = Polyhedron(vertices=[[0,0], [1,0], [0,1]],base_ring=RDF)
            sage: P.volume(engine='latte')
            Traceback (most recent call last):
            ...
            ValueError: LattE integrale cannot be applied over inexact rings
        """
        from sage.interfaces.latte import integrate
        from sage.rings.real_double import RDF

        if self.base_ring() == RDF:
            raise ValueError("LattE integrale cannot be applied over inexact rings")
        else:
            return integrate(self.cdd_Hrepresentation(), algorithm=algorithm, cdd=True, verbose=verbose, **kwargs)

    def _volume_normaliz(self, measure='induced'):
        r"""
        Compute the volume of a polytope using normaliz.

        INPUT:

        - ``measure`` -- (default: ``'induced'``) the measure to take;
          'induced' correspond to ``EuclideanVolume`` in normaliz and
          'induced_lattice' correspond to ``Volume`` in normaliz

        OUTPUT:

        A float value (when ``measure`` is 'induced') or a rational number
        (when ``measure`` is 'induced_lattice')

        .. NOTE::

            This function depends on Normaliz (i.e., the ``pynormaliz`` optional
            package). See the Normaliz documentation for further details.

        TESTS::

            sage: P = Polyhedron(vertices=[[0,0], [1,0], [0,1], [1,1]])
            sage: P._volume_normaliz()
            Traceback (most recent call last):
            ...
            TypeError: the backend should be normaliz
        """
        raise TypeError("the backend should be normaliz")

    @cached_method(do_pickle=True)
    def volume(self, measure='ambient', engine='auto', **kwds):
        """
        Return the volume of the polytope.

        INPUT:

        - ``measure`` -- string. The measure to use. Allowed values are:

          * ``ambient`` (default): Lebesgue measure of ambient space (volume)
          * ``induced``: Lebesgue measure of the affine hull (relative volume)
          * ``induced_rational``: Scaling of the Lebesgue measure for rational
            polytopes, such that the unit hypercube has volume 1
          * ``induced_lattice``: Scaling of the Lebesgue measure, such that the
            volume of the hypercube is factorial(n)

        - ``engine`` -- string. The backend to use. Allowed values are:

          * ``'auto'`` (default): choose engine according to measure
          * ``'internal'``: see :meth:`triangulate`
          * ``'TOPCOM'``: see :meth:`triangulate`
          * ``'lrs'``: use David Avis's lrs program (optional)
          * ``'latte'``: use LattE integrale program (optional)
          * ``'normaliz'``: use Normaliz program (optional)

        - ``**kwds`` -- keyword arguments that are passed to the
          triangulation engine

        OUTPUT: the volume of the polytope

        EXAMPLES::

            sage: polytopes.hypercube(3).volume()
            8
            sage: (polytopes.hypercube(3)*2).volume()
            64
            sage: polytopes.twenty_four_cell().volume()
            2

        Volume of the same polytopes, using the optional package lrslib
        (which requires a rational polytope)::

            sage: I3 = polytopes.hypercube(3)
            sage: I3.volume(engine='lrs')                               # optional - lrslib
            8
            sage: C24 = polytopes.twenty_four_cell()
            sage: C24.volume(engine='lrs')                              # optional - lrslib
            2

        If the base ring is exact, the answer is exact::

            sage: P5 = polytopes.regular_polygon(5)                                     # needs sage.rings.number_field
            sage: P5.volume()                                                           # needs sage.rings.number_field
            2.377641290737884?

            sage: polytopes.icosahedron().volume()                                      # needs sage.groups sage.rings.number_field
            5/12*sqrt5 + 5/4
            sage: numerical_approx(_)  # abs tol 1e9                                    # needs sage.groups sage.rings.number_field
            2.18169499062491

        When considering lower-dimensional polytopes, we can ask for the
        ambient (full-dimensional), the induced measure (of the affine
        hull) or, in the case of lattice polytopes, for the induced rational measure.
        This is controlled by the parameter ``measure``. Different engines
        may have different ideas on the definition of volume of a
        lower-dimensional object::

            sage: P = Polyhedron([[0, 0], [1, 1]])
            sage: P.volume()
            0
            sage: P.volume(measure='induced')                                           # needs sage.rings.number_field
            1.414213562373095?
            sage: P.volume(measure='induced_rational')                  # optional - latte_int
            1

            sage: # needs sage.rings.number_field
            sage: S = polytopes.regular_polygon(6); S
            A 2-dimensional polyhedron in AA^2 defined as the convex hull of 6 vertices
            sage: edge = S.faces(1)[4].as_polyhedron()
            sage: edge.vertices()
            (A vertex at (0.866025403784439?, 1/2), A vertex at (0, 1))
            sage: edge.volume()
            0
            sage: edge.volume(measure='induced')
            1

            sage: # optional - pynormaliz
            sage: P = Polyhedron(backend='normaliz',
            ....:                vertices=[[1,0,0], [0,0,1],
            ....:                          [-1,1,1], [-1,2,0]])
            sage: P.volume()
            0
            sage: P.volume(measure='induced')                                           # needs sage.rings.number_field
            2.598076211353316?
            sage: P.volume(measure='induced', engine='normaliz')
            2.598076211353316
            sage: P.volume(measure='induced_rational')                  # optional - latte_int
            3/2
            sage: P.volume(measure='induced_rational',
            ....:          engine='normaliz')
            3/2
            sage: P.volume(measure='induced_lattice')
            3

        The same polytope without normaliz backend::

            sage: P = Polyhedron(vertices=[[1,0,0], [0,0,1], [-1,1,1], [-1,2,0]])
            sage: P.volume(measure='induced_lattice', engine='latte')   # optional - latte_int
            3

            sage: # needs sage.groups sage.rings.number_field
            sage: Dexact = polytopes.dodecahedron()
            sage: F0 = Dexact.faces(2)[0].as_polyhedron()
            sage: v = F0.volume(measure='induced', engine='internal'); v
            1.53406271079097?
            sage: F4 = Dexact.faces(2)[4].as_polyhedron()
            sage: v = F4.volume(measure='induced', engine='internal'); v
            1.53406271079097?
            sage: RDF(v)    # abs tol 1e-9
            1.53406271079044

            sage: # needs sage.groups
            sage: Dinexact = polytopes.dodecahedron(exact=False)
            sage: F2 = Dinexact.faces(2)[2].as_polyhedron()
            sage: w = F2.volume(measure='induced', engine='internal')
            sage: RDF(w)    # abs tol 1e-9
            1.5340627082974878

            sage: all(polytopes.simplex(d).volume(measure='induced')                    # needs sage.rings.number_field sage.symbolic
            ....:        == sqrt(d+1)/factorial(d)
            ....:     for d in range(1,5))
            True

            sage: I = Polyhedron([[-3, 0], [0, 9]])
            sage: I.volume(measure='induced')                                           # needs sage.rings.number_field
            9.48683298050514?
            sage: I.volume(measure='induced_rational')                  # optional - latte_int
            3

            sage: T = Polyhedron([[3, 0, 0], [0, 4, 0], [0, 0, 5]])
            sage: T.volume(measure='induced')                                           # needs sage.rings.number_field
            13.86542462386205?
            sage: T.volume(measure='induced_rational')                  # optional - latte_int
            1/2

            sage: Q = Polyhedron(vertices=[(0, 0, 1, 1), (0, 1, 1, 0), (1, 1, 0, 0)])
            sage: Q.volume(measure='induced')
            1
            sage: Q.volume(measure='induced_rational')                  # optional - latte_int
            1/2

        The volume of a full-dimensional unbounded polyhedron is infinity::

            sage: P = Polyhedron(vertices=[[1, 0], [0, 1]], rays=[[1, 1]])
            sage: P.volume()
            +Infinity

        The volume of a non full-dimensional unbounded polyhedron depends on the measure used::

            sage: P = Polyhedron(ieqs = [[1,1,1], [-1,-1,-1], [3,1,0]]); P
            A 1-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 1 ray
            sage: P.volume()
            0
            sage: P.volume(measure='induced')
            +Infinity
            sage: P.volume(measure='ambient')
            0
            sage: P.volume(measure='induced_rational')                  # optional - pynormaliz
            +Infinity
            sage: P.volume(measure='induced_rational',engine='latte')
            +Infinity

        The volume in `0`-dimensional space is taken by counting measure::

            sage: P = Polyhedron(vertices=[[]]); P
            A 0-dimensional polyhedron in ZZ^0 defined as the convex hull of 1 vertex
            sage: P.volume()
            1
            sage: P = Polyhedron(vertices=[]); P
            The empty polyhedron in ZZ^0
            sage: P.volume()
            0

        TESTS:

        The cache of the volume is being pickled::

            sage: P = polytopes.cube()
            sage: P.volume()
            8
            sage: Q = loads(dumps(P))
            sage: Q.volume.is_in_cache()
            True

        Induced volumes work with lrs (:issue:`33410`)::

            sage: P = Polyhedron([[0, 0], [1, 1]])
            sage: P.volume(measure='induced', engine='lrs')             # optional - lrslib
            1.414213562373095?
        """
        from sage.features import FeatureNotPresentError
        if measure == 'induced_rational' and engine not in ['auto', 'latte', 'normaliz']:
            raise RuntimeError("the induced rational measure can only be computed with the engine set to `auto`, `latte`, or `normaliz`")
        if measure == 'induced_lattice' and engine not in ['auto', 'latte', 'normaliz']:
            raise RuntimeError("the induced lattice measure can only be computed with the engine set to `auto`, `latte`, or `normaliz`")
        if engine == 'auto' and measure == 'induced_rational':
            # Enforce a default choice, change if a better engine is found.
            from sage.features.latte import Latte
            try:
                Latte().require()
                engine = 'latte'
            except FeatureNotPresentError:
                from sage.features.normaliz import PyNormaliz
                try:
                    PyNormaliz().require()
                    engine = 'normaliz'
                except FeatureNotPresentError:
                    raise RuntimeError("the induced rational measure can only be computed with the optional packages `latte_int`, or `pynormaliz`")

        if engine == 'auto' and measure == 'induced_lattice':
            # Enforce a default choice, change if a better engine is found.
            from sage.features.normaliz import PyNormaliz
            try:
                PyNormaliz().require()
                engine = 'normaliz'
            except FeatureNotPresentError:
                try:
                    from sage.features.latte import Latte
                    Latte().require()
                    engine = 'latte'
                except FeatureNotPresentError:
                    raise RuntimeError("the induced rational measure can only be computed with the optional packages `latte_int`, or `pynormaliz`")

        if engine == 'auto' and measure == 'ambient' and self.backend() == 'normaliz':
            engine = 'normaliz'

        if measure == 'ambient':
            if self.dim() < self.ambient_dim():
                return self.base_ring().zero()
            elif self.dim() == 0:
                return 1
            # if the polyhedron is unbounded, return infinity
            if not self.is_compact():
                from sage.rings.infinity import infinity
                return infinity
            if engine == 'lrs':
                return self._volume_lrs(**kwds)
            elif engine == 'latte':
                return self._volume_latte(**kwds)
            elif engine == 'normaliz':
                return self._volume_normaliz(measure='ambient')

            triangulation = self.triangulate(engine=engine, **kwds)
            pc = triangulation.point_configuration()
            return sum([pc.volume(simplex) for simplex in triangulation]) / ZZ(self.dim()).factorial()
        elif measure == 'induced':
            # if polyhedron is actually full-dimensional, return volume with ambient measure
            if self.dim() == self.ambient_dim():
                return self.volume(measure='ambient', engine=engine, **kwds)
            # if the polyhedron is unbounded, return infinity
            if not self.is_compact():
                from sage.rings.infinity import infinity
                return infinity
            if engine == 'normaliz':
                return self._volume_normaliz(measure='euclidean')
            # use an orthogonal transformation, which preserves volume up to a factor provided by the transformation matrix
            affine_hull_data = self.affine_hull_projection(orthogonal=True, as_polyhedron=True, as_affine_map=True)
            A = affine_hull_data.projection_linear_map.matrix()
            Adet = (A.transpose() * A).det()
            scaled_volume = affine_hull_data.image.volume(measure='ambient', engine=engine, **kwds)
            if Adet.is_square():
                sqrt_Adet = Adet.sqrt()
            else:
                from sage.rings.qqbar import AA
                sqrt_Adet = AA(Adet).sqrt()
                scaled_volume = AA(scaled_volume)
            return scaled_volume / sqrt_Adet
        elif measure == 'induced_rational':
            # if the polyhedron is unbounded, return infinity
            if not self.is_compact():
                from sage.rings.infinity import infinity
                return infinity
            if engine == 'latte':
                return self._volume_latte(**kwds)
            else:  # engine is 'normaliz'
                return self._volume_normaliz(measure='induced_lattice') / ZZ(self.dim()).factorial()
        elif measure == 'induced_lattice':
            # if the polyhedron is unbounded, return infinity
            if not self.is_compact():
                from sage.rings.infinity import infinity
                return infinity
            if engine == 'latte':
                return self._volume_latte(**kwds) * ZZ(self.dim()).factorial()
            else:  # engine is 'normaliz'
                return self._volume_normaliz(measure='induced_lattice')
        else:
            raise TypeError("the measure should be `ambient`, `induced`, `induced_rational`, or `induced_lattice`")

    def integrate(self, function, measure='ambient', **kwds):
        r"""
        Return the integral of ``function`` over this polytope.

        INPUT:

        - ``self`` -- Polyhedron

        - ``function`` -- a multivariate polynomial or
          a valid LattE description string for polynomials

        - ``measure`` -- string, the measure to use

          Allowed values are:

          * ``ambient`` (default): Lebesgue measure of ambient space,
          * ``induced``: Lebesgue measure of the affine hull,
          * ``induced_nonnormalized``: Lebesgue measure of the affine hull
            without the normalization by `\sqrt{\det(A^\top A)}` (with
            `A` being the affine transformation matrix; see :meth:`affine_hull`).

        - ``**kwds`` -- additional keyword arguments that
          are passed to the engine

        OUTPUT: the integral of the polynomial over the polytope

        .. NOTE::

            The polytope triangulation algorithm is used. This function depends
            on LattE (i.e., the ``latte_int`` optional package).

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: x, y, z = polygens(QQ, 'x, y, z')
            sage: P.integrate(x^2*y^2*z^2)                              # optional - latte_int
            8/27

        If the polyhedron has floating point coordinates, an inexact result can
        be obtained if we transform to rational coordinates::

            sage: P = 1.4142*polytopes.cube()
            sage: P_QQ = Polyhedron(vertices=[[QQ(vi) for vi in v] for v in P.vertex_generator()])
            sage: RDF(P_QQ.integrate(x^2*y^2*z^2))                      # optional - latte_int
            6.703841212195228

        Integral over a non full-dimensional polytope::

            sage: x, y = polygens(QQ, 'x, y')
            sage: P = Polyhedron(vertices=[[0,0], [1,1]])
            sage: P.integrate(x*y)
            0
            sage: ixy = P.integrate(x*y, measure='induced'); ixy        # optional - latte_int
            0.4714045207910317?
            sage: ixy.parent()                                          # optional - latte_int
            Algebraic Real Field

        Convert to a symbolic expression::

            sage: ixy.radical_expression()                              # optional - latte_int
            1/3*sqrt(2)

        Another non full-dimensional polytope integration::

            sage: R.<x, y, z> = QQ[]
            sage: P = polytopes.simplex(2)
            sage: V = AA(P.volume(measure='induced'))                                   # needs sage.rings.number_field
            sage: V.radical_expression()                                                # needs sage.rings.number_field sage.symbolic
            1/2*sqrt(3)
            sage: P.integrate(R(1), measure='induced') == V             # optional - latte_int, needs sage.rings.number_field sage.symbolic
            True

        Computing the mass center::

            sage: (P.integrate(x, measure='induced')                    # optional - latte_int, needs sage.rings.number_field sage.symbolic
            ....:     / V).radical_expression()
            1/3
            sage: (P.integrate(y, measure='induced')                    # optional - latte_int, needs sage.rings.number_field sage.symbolic
            ....:     / V).radical_expression()
            1/3
            sage: (P.integrate(z, measure='induced')                    # optional - latte_int, needs sage.rings.number_field sage.symbolic
            ....:     / V).radical_expression()
            1/3

        TESTS:

        Testing a three-dimensional integral::

            sage: P = polytopes.octahedron()
            sage: x, y, z = polygens(QQ, 'x, y, z')
            sage: P.integrate(2*x^2*y^4*z^6 + z^2)                      # optional - latte_int
            630632/4729725

        Testing a polytope with non-rational vertices::

            sage: P = polytopes.icosahedron()                                           # needs sage.groups sage.rings.number_field
            sage: P.integrate(x^2*y^2*z^2)                              # optional - latte_int, needs sage.groups sage.rings.number_field
            Traceback (most recent call last):
            ...
            TypeError: the base ring must be ZZ, QQ, or RDF

        Testing a univariate polynomial::

            sage: P = Polyhedron(vertices=[[0], [1]])
            sage: x = polygen(QQ, 'x')
            sage: P.integrate(x)                                        # optional - latte_int
            1/2

        Testing a polytope with floating point coordinates::

            sage: P = Polyhedron(vertices=[[0, 0], [1, 0], [1.1, 1.1], [0, 1]])
            sage: P.integrate('[[1,[2,2]]]')
            Traceback (most recent call last):
            ...
            TypeError: LattE integrale cannot be applied over inexact rings

        Integration of zero-polynomial::

            sage: R.<x, y, z> = QQ[]
            sage: P = polytopes.simplex(2)
            sage: P.integrate(R(0))
            0
            sage: P.integrate('[]')  # with LattE description string
            0

        ::

            sage: R.<x, y, z> = QQ[]
            sage: P = Polyhedron(vertices=[(0, 0, 1), (0, 1, 0)])
            sage: P.integrate(x^2)
            0
        """
        if function == 0 or function == '[]':
            return self.base_ring().zero()

        if not self.is_compact():
            raise NotImplementedError(
                'integration over non-compact polyhedra not allowed')

        if measure == 'ambient':
            if not self.is_full_dimensional():
                return self.base_ring().zero()

            return self._integrate_latte_(function, **kwds)

        elif measure == 'induced' or measure == 'induced_nonnormalized':
            # if polyhedron is actually full-dimensional,
            # return with ambient measure
            if self.is_full_dimensional():
                return self.integrate(function, measure='ambient', **kwds)

            if isinstance(function, str):
                raise NotImplementedError(
                    'LattE description strings for polynomials not allowed '
                    'when using measure="induced"')

            # use an orthogonal transformation
            affine_hull_data = self.affine_hull_projection(orthogonal=True, return_all_data=True)
            polyhedron = affine_hull_data.image
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            R = PolynomialRing(affine_hull_data.section_linear_map.base_ring(), 'x', self.dim())
            coordinate_images = affine_hull_data.section_linear_map.matrix().transpose() * vector(R.gens()) + affine_hull_data.section_translation

            hom = function.parent().hom(coordinate_images)
            function_in_affine_hull = hom(function)

            I = polyhedron.integrate(function_in_affine_hull,
                                     measure='ambient', **kwds)
            if measure == 'induced_nonnormalized':
                return I
            else:
                A = affine_hull_data.projection_linear_map.matrix()
                Adet = (A.transpose() * A).det()
                try:
                    from sage.rings.qqbar import AA
                    Adet = AA.coerce(Adet)
                except TypeError:
                    pass
                return I / Adet.sqrt()

        else:
            raise ValueError('unknown measure "{}"'.format(measure))

    def _integrate_latte_(self, polynomial, **kwds):
        r"""
        Return the integral of a polynomial over this polytope by calling LattE.

        INPUT:

        - ``polynomial`` -- a multivariate polynomial or
          a valid LattE description string for polynomials

        - ``**kwds`` -- additional keyword arguments that are passed
          to the engine

        OUTPUT: the integral of the polynomial over the polytope

        .. NOTE::

            The polytope triangulation algorithm is used. This function depends
            on LattE (i.e., the ``latte_int`` optional package).

        TESTS::

            sage: P = polytopes.cube()
            sage: x, y, z = polygens(QQ, 'x, y, z')
            sage: P._integrate_latte_(x^2 + y^2*z^2)                    # optional - latte_int
            32/9

        ::

            sage: R = PolynomialRing(QQ, '', 0)
            sage: Polyhedron(vertices=[()]).integrate(R(42))
            42
        """
        from sage.rings.real_double import RDF

        if self.base_ring() == RDF:
            raise TypeError("LattE integrale cannot be applied over inexact rings")
        if self.dimension() == 0:
            vertices = self.vertices()
            assert len(self.vertices()) == 1
            vertex = tuple(vertices[0])
            return polynomial(vertex)

        from sage.interfaces.latte import integrate
        return integrate(self.cdd_Hrepresentation(),
                         polynomial,
                         cdd=True, **kwds)
