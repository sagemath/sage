# -*- coding: utf-8 -*- 
"""
The Normaliz backend for polyhedral computations

.. NOTE::

    This backend requires `PyNormaliz <https://pypi.python.org/pypi/PyNormaliz/1.5>`_.
    To install PyNormaliz, type :code:`sage -i pynormaliz` in the terminal.

AUTHORS:

- Matthias Köppe (2016-12): initial version
"""

#*****************************************************************************
#  Copyright (C) 2016 Matthias Köppe <mkoeppe at math.ucdavis.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import absolute_import, print_function

from sage.structure.element import Element
from sage.misc.all import prod

from sage.rings.all import ZZ, QQ
from sage.arith.functions import LCM_list
from sage.misc.functional import denominator
from sage.matrix.constructor import matrix, vector

from .base import Polyhedron_base
from .base_QQ import Polyhedron_QQ
from .base_ZZ import Polyhedron_ZZ


#########################################################################
class Polyhedron_normaliz(Polyhedron_base):
    """
    Polyhedra with normaliz

    INPUT:

    - ``parent`` -- :class:`~sage.geometry.polyhedron.parent.Polyhedra`
      the parent

    - ``Vrep`` -- a list ``[vertices, rays, lines]`` or ``None``; the
      V-representation of the polyhedron; if ``None``, the polyhedron
      is determined by the H-representation

    - ``Hrep`` -- a list ``[ieqs, eqns]`` or ``None``; the
      H-representation of the polyhedron; if ``None``, the polyhedron
      is determined by the V-representation

    - ``normaliz_cone`` -- a PyNormaliz wrapper of a normaliz cone

    Only one of ``Vrep``, ``Hrep``, or ``normaliz_cone`` can be different
    from ``None``.

    EXAMPLES::

        sage: p = Polyhedron(vertices=[(0,0),(1,0),(0,1)], rays=[(1,1)],   # optional - pynormaliz
        ....:                lines=[], backend='normaliz')
        sage: TestSuite(p).run(skip='_test_pickling')                      # optional - pynormaliz

    Two ways to get the full space::

        sage: Polyhedron(eqns=[[0, 0, 0]], backend='normaliz')             # optional - pynormaliz
        A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 lines
        sage: Polyhedron(ieqs=[[0, 0, 0]], backend='normaliz')             # optional - pynormaliz
        A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 lines

    A lower-dimensional affine cone; we test that there are no mysterious
    inequalities coming in from the homogenization::

        sage: P = Polyhedron(vertices=[(1, 1)], rays=[(0, 1)],             # optional - pynormaliz
        ....:                backend='normaliz')
        sage: P.n_inequalities()                                           # optional - pynormaliz
        1
        sage: P.equations()                                                # optional - pynormaliz
        (An equation (-1, 0) x + 1 == 0,)

    The empty polyhedron::

        sage: P=Polyhedron(ieqs=[[-2, 1, 1], [-3, -1, -1], [-4, 1, -2]],   # optional - pynormaliz
        ....:              backend='normaliz')
        sage: P                                                            # optional - pynormaliz
        The empty polyhedron in QQ^2
        sage: P.Vrepresentation()                                          # optional - pynormaliz
        ()
        sage: P.Hrepresentation()                                          # optional - pynormaliz
        (An equation -1 == 0,)

    TESTS:

    Tests copied from various methods in :mod:`sage.geometry.polyhedron.base`::

        sage: p = Polyhedron(vertices = [[1,0,0], [0,1,0], [0,0,1]],       # optional - pynormaliz
        ....:                backend='normaliz')
        sage: p.n_equations()                                              # optional - pynormaliz
        1
        sage: p.n_inequalities()                                           # optional - pynormaliz
        3

        sage: p = Polyhedron(vertices = [[t,t^2,t^3] for t in range(6)],   # optional - pynormaliz
        ....:                backend='normaliz')
        sage: p.n_facets()                                                 # optional - pynormaliz
        8

        sage: p = Polyhedron(vertices = [[1,0],[0,1],[1,1]], rays=[[1,1]], # optional - pynormaliz
        ....:                backend='normaliz')
        sage: p.n_vertices()                                               # optional - pynormaliz
        2

        sage: p = Polyhedron(vertices = [[1,0],[0,1]], rays=[[1,1]],       # optional - pynormaliz
        ....:                backend='normaliz')
        sage: p.n_rays()                                                   # optional - pynormaliz
        1

        sage: p = Polyhedron(vertices = [[0,0]], rays=[[0,1],[0,-1]],      # optional - pynormaliz
        ....:                backend='normaliz')
        sage: p.n_lines()                                                  # optional - pynormaliz
        1

    """
    def __init__(self, parent, Vrep, Hrep, normaliz_cone=None, **kwds):
        """
        Initializes the polyhedron.

        See :class:`Polyhedron_normaliz` for a description of the input
        data.

        TESTS:

        We skip the pickling test because pickling is currently
        not implemented::

            sage: p = Polyhedron(backend='normaliz')                 # optional - pynormaliz
            sage: TestSuite(p).run(skip="_test_pickling")            # optional - pynormaliz
            sage: p = Polyhedron(vertices=[(1, 1)], rays=[(0, 1)],   # optional - pynormaliz
            ....:                backend='normaliz')
            sage: TestSuite(p).run(skip="_test_pickling")            # optional - pynormaliz
            sage: p = Polyhedron(vertices=[(-1,-1), (1,0), (1,1), (0,1)],  # optional - pynormaliz
            ....:                backend='normaliz')
            sage: TestSuite(p).run(skip="_test_pickling")            # optional - pynormaliz
        """
        if normaliz_cone:
            if Hrep is not None or Vrep is not None:
                raise ValueError("only one of Vrep, Hrep, or normaliz_cone can be different from None")
            Element.__init__(self, parent=parent)
            self._init_from_normaliz_cone(normaliz_cone)
        else:
            Polyhedron_base.__init__(self, parent, Vrep, Hrep, **kwds)

    def _init_from_normaliz_cone(self, normaliz_cone):
        """
        Construct polyhedron from a PyNormaliz wrapper of a normaliz cone.

        TESTS::

            sage: p = Polyhedron(backend='normaliz')                       # optional - pynormaliz
            sage: from sage.geometry.polyhedron.backend_normaliz import Polyhedron_normaliz   # optional - pynormaliz
            sage: Polyhedron_normaliz._init_from_Hrepresentation(p, [], [])  # indirect doctest  # optional - pynormaliz
        """
        import PyNormaliz
        if normaliz_cone and PyNormaliz.NmzResult(normaliz_cone, "AffineDim") < 0:
            # Empty polyhedron. Special case because Normaliz defines the
            # recession cone of an empty polyhedron given by an
            # H-representation as the cone defined by the homogenized system.
            self._init_empty_polyhedron()
        else:
            self._normaliz_cone = normaliz_cone
            self._init_Vrepresentation_from_normaliz()
            self._init_Hrepresentation_from_normaliz()

    def _init_from_normaliz_data(self, data, verbose=False):
        """
        Construct polyhedron from normaliz ``data`` (a dictionary).

        TESTS::

            sage: p = Polyhedron(backend='normaliz', ambient_dim=2)                           # optional - pynormaliz
            sage: from sage.geometry.polyhedron.backend_normaliz import Polyhedron_normaliz   # optional - pynormaliz
            sage: data = {'inhom_inequalities': [[-1L, 2L, 0L], [0L, 0L, 1L], [2L, -1L, 0L]]} # optional - pynormaliz
            sage: Polyhedron_normaliz._init_from_normaliz_data(p, data)                       # optional - pynormaliz
            sage: p.inequalities_list()                                                       # optional - pynormaliz
            [[0, -1, 2], [0, 2, -1]]
        """
        import PyNormaliz
        if verbose:
            print("# Calling PyNormaliz.NmzCone(**{})".format(data))
            import six
            if isinstance(verbose, six.string_types):
                print("# Wrote equivalent Normaliz input file to {}".format(verbose))
                self._normaliz_format(data, file_output=verbose)
            else:
                print("# ----8<---- Equivalent Normaliz input file ----8<----")
                print(self._normaliz_format(data), end='')
                print("# ----8<-------------------8<-------------------8<----")
        cone = PyNormaliz.NmzCone(**data)
        assert cone, "NmzCone(**{}) did not return a cone".format(data)
        self._init_from_normaliz_cone(cone)

    def _init_from_Vrepresentation(self, vertices, rays, lines, minimize=True, verbose=False):
        r"""
        Construct polyhedron from V-representation data.

        INPUT:

        - ``vertices`` -- list of point; each point can be specified
           as any iterable container of
           :meth:`~sage.geometry.polyhedron.base.base_ring` elements

        - ``rays`` -- list of rays; each ray can be specified as any
          iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements

        - ``lines`` -- list of lines; each line can be specified as
          any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements

        - ``verbose`` -- boolean (default: ``False``); whether to print
          verbose output for debugging purposes

        EXAMPLES::

            sage: p = Polyhedron(backend='normaliz')                       # optional - pynormaliz
            sage: from sage.geometry.polyhedron.backend_normaliz import Polyhedron_normaliz   # optional - pynormaliz
            sage: Polyhedron_normaliz._init_from_Vrepresentation(p, [], [], [])   # optional - pynormaliz
        """
        if vertices is None:
            vertices = []
        nmz_vertices = []
        for v in vertices:
            d = LCM_list([denominator(v_i) for v_i in v])
            dv = [ d*v_i for v_i in v ]
            nmz_vertices.append(dv + [d])
        if rays is None:
            rays = []
        nmz_rays = []
        for r in rays:
            d = LCM_list([denominator(r_i) for r_i in r])
            dr = [ d*r_i for r_i in r ]
            nmz_rays.append(dr)
        if lines is None: lines = []
        nmz_lines = []
        for l in lines:
            d = LCM_list([denominator(l_i) for l_i in l])
            dl = [ d*l_i for l_i in l ]
            nmz_lines.append(dl)
        if not nmz_vertices and not nmz_rays and not nmz_lines:
            # Special case to avoid:
            #   error: Some error in the normaliz input data detected:
            #   All input matrices empty!
            self._init_empty_polyhedron()
        else:
            data = {"vertices": nmz_vertices,
                    "cone": nmz_rays,
                    "subspace": nmz_lines}
            self._init_from_normaliz_data(data, verbose=verbose)

    def _init_from_Hrepresentation(self, ieqs, eqns, minimize=True, verbose=False):
        r"""
        Construct polyhedron from H-representation data.

        INPUT:

        - ``ieqs`` -- list of inequalities; each line can be specified
          as any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements

        - ``eqns`` -- list of equalities; each line can be specified
          as any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements

        - ``minimize`` -- boolean (default: ``True``); ignored

        - ``verbose`` -- boolean (default: ``False``); whether to print
          verbose output for debugging purposes

        EXAMPLES::

            sage: p = Polyhedron(backend='normaliz')                       # optional - pynormaliz
            sage: from sage.geometry.polyhedron.backend_normaliz import Polyhedron_normaliz   # optional - pynormaliz
            sage: Polyhedron_normaliz._init_from_Hrepresentation(p, [], [])   # optional - pynormaliz
        """
        if ieqs is None: ieqs = []
        nmz_ieqs = []
        for ieq in ieqs:
            d = LCM_list([denominator(ieq_i) for ieq_i in ieq])
            dieq = [ ZZ(d*ieq_i) for ieq_i in ieq ]
            b = dieq[0]
            A = dieq[1:]
            nmz_ieqs.append(A + [b])
        if not nmz_ieqs:
            # If normaliz gets an empty list of inequalities, it adds
            # nonnegativities. So let's add a tautological inequality to work
            # around this.
            nmz_ieqs.append([0]*self.ambient_dim() + [0])
        if eqns is None: eqns = []
        nmz_eqns = []
        for eqn in eqns:
            d = LCM_list([denominator(eqn_i) for eqn_i in eqn])
            deqn = [ ZZ(d*eqn_i) for eqn_i in eqn ]
            b = deqn[0]
            A = deqn[1:]
            nmz_eqns.append(A + [b])
        data = {"inhom_equations": nmz_eqns,
                "inhom_inequalities": nmz_ieqs}
        self._init_from_normaliz_data(data, verbose=verbose)

    def _init_Vrepresentation_from_normaliz(self):
        r"""
        Create the Vrepresentation objects from the normaliz polyhedron.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[(0,1/2),(2,0),(4,5/6)],  # indirect doctest # optional - pynormaliz
            ....:                backend='normaliz')
            sage: p.Hrepresentation()                               # optional - pynormaliz
            (An inequality (-5, 12) x + 10 >= 0,
             An inequality (1, -12) x + 6 >= 0,
             An inequality (1, 4) x - 2 >= 0)
            sage: p.Vrepresentation()                               # optional - pynormaliz
            (A vertex at (0, 1/2), A vertex at (2, 0), A vertex at (4, 5/6))
        """
        import PyNormaliz
        self._Vrepresentation = []
        parent = self.parent()
        base_ring = self.base_ring()
        cone = self._normaliz_cone
        for g in PyNormaliz.NmzResult(cone, "VerticesOfPolyhedron"):
            d = g[-1]
            if d == 1:
                parent._make_Vertex(self, g[:-1])
            else:
                parent._make_Vertex(self, [base_ring(x)/d for x in g[:-1]])
        for g in PyNormaliz.NmzResult(cone, "ExtremeRays"):
            parent._make_Ray(self, g[:-1])
        for g in PyNormaliz.NmzResult(cone, "MaximalSubspace"):
            parent._make_Line(self, g[:-1])
        self._Vrepresentation = tuple(self._Vrepresentation)

    def _init_Hrepresentation_from_normaliz(self):
        r"""
        Create the Hrepresentation objects from the normaliz polyhedron.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[(0,1/2), (2,0), (4,5/6)],  # indirect doctest # optional - pynormaliz
            ....:                backend='normaliz')
            sage: p.Hrepresentation()                                 # optional - pynormaliz
            (An inequality (-5, 12) x + 10 >= 0,
             An inequality (1, -12) x + 6 >= 0,
             An inequality (1, 4) x - 2 >= 0)
            sage: p.Vrepresentation()                                 # optional - pynormaliz
            (A vertex at (0, 1/2), A vertex at (2, 0), A vertex at (4, 5/6))
        """
        import PyNormaliz
        self._Hrepresentation = []
        cone = self._normaliz_cone
        parent = self.parent()
        for g in PyNormaliz.NmzResult(cone, "SupportHyperplanes"):
            if all(x == 0 for x in g[:-1]):
                # Ignore vertical inequality
                pass
            else:
                parent._make_Inequality(self, (g[-1],) + tuple(g[:-1]))
        for g in PyNormaliz.NmzResult(cone, "Equations"):
            parent._make_Equation(self, (g[-1],) + tuple(g[:-1]))
        self._Hrepresentation = tuple(self._Hrepresentation)

    def _init_empty_polyhedron(self):
        r"""
        Initializes an empty polyhedron.

        TESTS::

            sage: empty = Polyhedron(backend='normaliz'); empty            # optional - pynormaliz
            The empty polyhedron in ZZ^0
            sage: empty.Vrepresentation()                                  # optional - pynormaliz
            ()
            sage: empty.Hrepresentation()                                  # optional - pynormaliz
            (An equation -1 == 0,)
            sage: Polyhedron(vertices = [], backend='normaliz')            # optional - pynormaliz
            The empty polyhedron in ZZ^0
            sage: Polyhedron(backend='normaliz')._init_empty_polyhedron()  # optional - pynormaliz
        """
        super(Polyhedron_normaliz, self)._init_empty_polyhedron()
        # Can't seem to set up an empty _normaliz_cone.
        # For example, PyNormaliz.NmzCone(vertices=[]) gives
        # error: Some error in the normaliz input data detected: All input matrices empty!
        self._normaliz_cone = None

    @classmethod
    def _from_normaliz_cone(cls, parent, normaliz_cone):
        r"""
        Initializes a polyhedron from a PyNormaliz wrapper of a normaliz cone.

        TESTS::

            sage: P=Polyhedron(ieqs=[[1, 0, 2], [3, 0, -2], [3, 2, -2]],   # optional - pynormaliz
            ....:              backend='normaliz')
            sage: PI = P.integral_hull()                 # indirect doctest; optional - pynormaliz
        """
        return cls(parent, None, None, normaliz_cone=normaliz_cone)

    @staticmethod
    def _make_normaliz_cone(data, verbose=False):
        r"""
        Returns a normaliz cone from ``data``.

        INPUT:

        - ``data`` -- a dictionary

        - ``verbose`` -- a boolean (default: ``False``)

        TESTS::

            sage: from sage.geometry.polyhedron.backend_normaliz import Polyhedron_normaliz    # optional - pynormaliz
            sage: data = {'inhom_inequalities': [[-1L, 2L, 0L], [0L, 0L, 1L], [2L, -1L, 0L]]}  # optional - pynormaliz
            sage: nmz_cone = Polyhedron_normaliz._make_normaliz_cone(data,verbose=False)       # optional - pynormaliz
            sage: from PyNormaliz import NmzResult                                             # optional - pynormaliz
            sage: NmzResult(nmz_cone, "ExtremeRays")                                           # optional - pynormaliz
            [[1L, 2L, 0L], [2L, 1L, 0L]]
        """
        import PyNormaliz
        if verbose:
            print("# Calling PyNormaliz.NmzCone(**{})".format(data))
        cone = PyNormaliz.NmzCone(**data)
        assert cone, "NmzCone(**{}) did not return a cone".format(data)
        return cone

    def _get_nmzcone_data(self):
        r"""
        Get the data necessary to reproduce the normaliz cone.

        OUTPUT:

        - ``data`` -- a dictionary.

        TESTS:

        The empty polyhedron::

            sage: P = Polyhedron(backend='normaliz')                               # optional - pynormaliz
            sage: P._get_nmzcone_data()                                            # optional - pynormaliz
            {}

        Another simple example::

            sage: C = Polyhedron(backend='normaliz',rays=[[1,2],[2,1]])            # optional - pynormaliz
            sage: C._get_nmzcone_data()                                            # optional - pynormaliz
            {'cone': [[1L, 2L], [2L, 1L]],
             'inhom_equations': [],
             'inhom_inequalities': [[-1L, 2L, 0L], [0L, 0L, 1L], [2L, -1L, 0L]],
             'subspace': [],
             'vertices': [[0L, 0L, 1L]]}
        """
        import PyNormaliz
        if self.is_empty():
            return {}

        vertices = PyNormaliz.NmzResult(self._normaliz_cone, "VerticesOfPolyhedron")
        # get rid of the last 0 in rays:
        rays = [r[:-1] for r in PyNormaliz.NmzResult(self._normaliz_cone, "ExtremeRays")]
        lines = PyNormaliz.NmzResult(self._normaliz_cone, "MaximalSubspace")
        ineqs = PyNormaliz.NmzResult(self._normaliz_cone, "SupportHyperplanes")
        eqs = PyNormaliz.NmzResult(self._normaliz_cone, "Equations")

        data = {'vertices': vertices,
                'cone': rays,
                'subspace': lines,
                'inhom_equations': eqs,
                'inhom_inequalities': ineqs}

        return data

    def _normaliz_format(self, data, file_output=None):
        r"""
        Return a string containing normaliz format.

        INPUT:

        - ``data`` -- a dictionary of PyNormaliz cone input properties

        - ``file_output`` (string; optional) -- a filename to which the
          representation should be written. If set to ``None`` (default),
          representation is returned as a string.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[[0, 0], [0, 1], [1, 0]], # indirect doctest; optional - pynormaliz
            ....:                backend='normaliz', verbose=True)
            # Calling ...
            # ----8<---- Equivalent Normaliz input file ----8<----
            amb_space 2
            subspace 0
            vertices 3
             0 0 1
             0 1 1
             1 0 1
            cone 0
            # ----8<-------------------8<-------------------8<----
        """
        s = 'amb_space {}\n'.format(self.ambient_dim())
        for key, value in data.iteritems():
            s += '{} {}\n'.format(key, len(value))
            for e in value:
                for x in e:
                    s += ' ' + repr(x)
                s += '\n'

        if file_output is not None:
            in_file = open(file_output, 'w')
            in_file.write(s)
            in_file.close()
        else:
            return s

    def integral_hull(self):
        r"""
        Return the integral hull in the polyhedron.

        This is a new polyhedron that is the convex hull of all integral
        points.

        EXAMPLES:

        Unbounded example from Normaliz manual, "a dull polyhedron"::

            sage: P = Polyhedron(ieqs=[[1, 0, 2], [3, 0, -2], [3, 2, -2]], # optional - pynormaliz
            ....:              backend='normaliz')
            sage: PI = P.integral_hull()                                   # optional - pynormaliz
            sage: P.plot(color='yellow') + PI.plot(color='green')          # optional - pynormaliz
            Graphics object consisting of 10 graphics primitives
            sage: PI.Vrepresentation()                                     # optional - pynormaliz
            (A vertex at (-1, 0), A vertex at (0, 1), A ray in the direction (1, 0))

        Nonpointed case::

            sage: P = Polyhedron(vertices=[[1/2, 1/3]], rays=[[1, 1]],     # optional - pynormaliz
            ....:              lines=[[-1, 1]], backend='normaliz')
            sage: PI = P.integral_hull()                                   # optional - pynormaliz
            sage: PI.Vrepresentation()                                     # optional - pynormaliz
            (A vertex at (1, 0),
             A ray in the direction (1, 0),
             A line in the direction (-1, 1))

        Empty polyhedron::

            sage: P = Polyhedron(backend='normaliz')                       # optional - pynormaliz
            sage: PI = P.integral_hull()                                   # optional - pynormaliz
            sage: PI.Vrepresentation()                                     # optional - pynormaliz
            ()
        """
        import PyNormaliz
        if self.is_empty():
            return self
        cone = PyNormaliz.NmzResult(self._normaliz_cone, "IntegerHull")
        return self.parent().element_class._from_normaliz_cone(parent=self.parent(),
                                                               normaliz_cone=cone)

    def integral_points(self, threshold=10000):
        r"""
        Return the integral points in the polyhedron.

        Uses either the naive algorithm (iterate over a rectangular
        bounding box) or triangulation + Smith form.

        INPUT:

        - ``threshold`` -- integer (default: 10000); use the naïve
          algorithm as long as the bounding box is smaller than this

        OUTPUT:

        The list of integral points in the polyhedron. If the
        polyhedron is not compact, a ``ValueError`` is raised.

        EXAMPLES::

            sage: Polyhedron(vertices=[(-1,-1), (1,0), (1,1), (0,1)],      # optional - pynormaliz
            ....:            backend='normaliz').integral_points()
            ((-1, -1), (0, 0), (0, 1), (1, 0), (1, 1))

            sage: simplex = Polyhedron([(1,2,3), (2,3,7), (-2,-3,-11)],    # optional - pynormaliz
            ....:                      backend='normaliz')
            sage: simplex.integral_points()                                # optional - pynormaliz
            ((-2, -3, -11), (0, 0, -2), (1, 2, 3), (2, 3, 7))

        The polyhedron need not be full-dimensional::

            sage: simplex = Polyhedron([(1,2,3,5), (2,3,7,5), (-2,-3,-11,5)],   # optional - pynormaliz
            ....:                      backend='normaliz')
            sage: simplex.integral_points()                                # optional - pynormaliz
            ((-2, -3, -11, 5), (0, 0, -2, 5), (1, 2, 3, 5), (2, 3, 7, 5))

            sage: point = Polyhedron([(2,3,7)],                            # optional - pynormaliz
            ....:                    backend='normaliz')
            sage: point.integral_points()                                  # optional - pynormaliz
            ((2, 3, 7),)

            sage: empty = Polyhedron(backend='normaliz')                   # optional - pynormaliz
            sage: empty.integral_points()                                  # optional - pynormaliz
            ()

        Here is a simplex where the naive algorithm of running over
        all points in a rectangular bounding box no longer works fast
        enough::

            sage: v = [(1,0,7,-1), (-2,-2,4,-3), (-1,-1,-1,4), (2,9,0,-5), (-2,-1,5,1)]
            sage: simplex = Polyhedron(v, backend='normaliz'); simplex     # optional - pynormaliz
            A 4-dimensional polyhedron in ZZ^4 defined as the convex hull of 5 vertices
            sage: len(simplex.integral_points())                           # optional - pynormaliz
            49

        A rather thin polytope for which the bounding box method would
        be a very bad idea (note this is a rational (non-lattice)
        polytope, so the other backends use the bounding box method)::

            sage: P = Polyhedron(vertices=((0, 0), (178933,37121))) + 1/1000*polytopes.hypercube(2)
            sage: P = Polyhedron(vertices=P.vertices_list(),               # optional - pynormaliz
            ....:                backend='normaliz')
            sage: len(P.integral_points())                                 # optional - pynormaliz
            434

        Finally, the 3-d reflexive polytope number 4078::

            sage: v = [(1,0,0), (0,1,0), (0,0,1), (0,0,-1), (0,-2,1),
            ....:      (-1,2,-1), (-1,2,-2), (-1,1,-2), (-1,-1,2), (-1,-3,2)]
            sage: P = Polyhedron(v, backend='normaliz')                    # optional - pynormaliz
            sage: pts1 = P.integral_points()                               # optional - pynormaliz
            sage: all(P.contains(p) for p in pts1)                         # optional - pynormaliz
            True
            sage: pts2 = LatticePolytope(v).points()          # PALP
            sage: for p in pts1: p.set_immutable()                         # optional - pynormaliz
            sage: set(pts1) == set(pts2)                                   # optional - pynormaliz
            True

            sage: timeit('Polyhedron(v, backend='normaliz').integral_points()')   # not tested - random
            625 loops, best of 3: 1.41 ms per loop
            sage: timeit('LatticePolytope(v).points()')       # not tested - random
            25 loops, best of 3: 17.2 ms per loop

        TESTS:

        Test some trivial cases (see :trac:`17937`):

        Empty polyhedron in 1 dimension::

            sage: P = Polyhedron(ambient_dim=1, backend='normaliz')        # optional - pynormaliz
            sage: P.integral_points()                                      # optional - pynormaliz
            ()

        Empty polyhedron in 0 dimensions::

            sage: P = Polyhedron(ambient_dim=0, backend='normaliz')        # optional - pynormaliz
            sage: P.integral_points()                                      # optional - pynormaliz
            ()

        Single point in 1 dimension::

            sage: P = Polyhedron([[3]], backend='normaliz')                # optional - pynormaliz
            sage: P.integral_points()                                      # optional - pynormaliz
            ((3),)

        Single non-integral point in 1 dimension::

            sage: P = Polyhedron([[1/2]], backend='normaliz')              # optional - pynormaliz
            sage: P.integral_points()                                      # optional - pynormaliz
            ()

        Single point in 0 dimensions::

            sage: P = Polyhedron([[]], backend='normaliz')                 # optional - pynormaliz
            sage: P.integral_points()                                      # optional - pynormaliz
            ((),)

        A polytope with no integral points (:trac:`22938`)::

            sage: ieqs = [[1, 2, -1, 0], [0, -1, 2, -1], [0, 0, -1, 2],
            ....:         [0, -1, 0, 0], [0, 0, -1, 0],  [0, 0, 0, -1],
            ....:         [-1, -1, -1, -1], [1, 1, 0, 0], [1, 0, 1, 0],
            ....:         [1, 0, 0, 1]]
            sage: P = Polyhedron(ieqs=ieqs, backend='normaliz')            # optional - pynormaliz
            sage: P.bounding_box()                                         # optional - pynormaliz
            ((-3/4, -1/2, -1/4), (-1/2, -1/4, 0))
            sage: P.bounding_box(integral_hull=True)                       # optional - pynormaliz
            (None, None)
            sage: P.integral_points()                                      # optional - pynormaliz
            ()

        Check the polytopes from :trac:`22984`::

            sage: base = [[0, 2, 0, -1, 0, 0, 0, 0, 0],
            ....:         [0, 0, 2, 0, -1, 0, 0, 0, 0],
            ....:         [1, -1, 0, 2, -1, 0, 0, 0, 0],
            ....:         [0, 0, -1, -1, 2, -1, 0, 0, 0],
            ....:         [0, 0, 0, 0, -1, 2, -1, 0, 0],
            ....:         [0, 0, 0, 0, 0, -1, 2, -1, 0],
            ....:         [1, 0, 0, 0, 0, 0, -1, 2, -1],
            ....:         [0, 0, 0, 0, 0, 0, 0, -1, 2],
            ....:         [0, -1, 0, 0, 0, 0, 0, 0, 0],
            ....:         [0, 0, -1, 0, 0, 0, 0, 0, 0],
            ....:         [0, 0, 0, -1, 0, 0, 0, 0, 0],
            ....:         [0, 0, 0, 0, -1, 0, 0, 0, 0],
            ....:         [0, 0, 0, 0, 0, -1, 0, 0, 0],
            ....:         [0, 0, 0, 0, 0, 0, -1, 0, 0],
            ....:         [0, 0, 0, 0, 0, 0, 0, -1, 0],
            ....:         [0, 0, 0, 0, 0, 0, 0, 0, -1],
            ....:         [-1, -1, -1, -1, -1, -1, -1, -1, -1]]

            sage: ieqs = base + [
            ....:         [2, 1, 0, 0, 0, 0, 0, 0, 0],
            ....:         [4, 0, 1, 0, 0, 0, 0, 0, 0],
            ....:         [4, 0, 0, 1, 0, 0, 0, 0, 0],
            ....:         [7, 0, 0, 0, 1, 0, 0, 0, 0],
            ....:         [6, 0, 0, 0, 0, 1, 0, 0, 0],
            ....:         [4, 0, 0, 0, 0, 0, 1, 0, 0],
            ....:         [2, 0, 0, 0, 0, 0, 0, 1, 0],
            ....:         [1, 0, 0, 0, 0, 0, 0, 0, 1]]
            sage: P = Polyhedron(ieqs=ieqs, backend='normaliz')            # optional - pynormaliz
            sage: P.integral_points()                                      # optional - pynormaliz
            ((-2, -2, -4, -5, -4, -3, -2, -1),
             (-2, -2, -4, -5, -4, -3, -2, 0),
             (-1, -2, -3, -4, -3, -2, -2, -1),
             (-1, -2, -3, -4, -3, -2, -1, 0),
             (-1, -1, -2, -2, -2, -2, -2, -1),
             (-1, -1, -2, -2, -1, -1, -1, 0),
             (-1, -1, -2, -2, -1, 0, 0, 0),
             (-1, 0, -2, -2, -2, -2, -2, -1),
             (0, -1, -1, -2, -2, -2, -2, -1),
             (0, 0, -1, -1, -1, -1, -1, 0))

            sage: ieqs = base + [
            ....:         [3, 1, 0, 0, 0, 0, 0, 0, 0],
            ....:         [4, 0, 1, 0, 0, 0, 0, 0, 0],
            ....:         [6, 0, 0, 1, 0, 0, 0, 0, 0],
            ....:         [8, 0, 0, 0, 1, 0, 0, 0, 0],
            ....:         [6, 0, 0, 0, 0, 1, 0, 0, 0],
            ....:         [4, 0, 0, 0, 0, 0, 1, 0, 0],
            ....:         [2, 0, 0, 0, 0, 0, 0, 1, 0],
            ....:         [1, 0, 0, 0, 0, 0, 0, 0, 1]]
            sage: P = Polyhedron(ieqs=ieqs, backend='normaliz')            # optional - pynormaliz
            sage: P.integral_points()                                      # optional - pynormaliz
            ((-3, -4, -6, -8, -6, -4, -2, -1),
             (-3, -4, -6, -8, -6, -4, -2, 0),
             (-2, -2, -4, -5, -4, -3, -2, -1),
             (-2, -2, -4, -5, -4, -3, -2, 0),
             (-1, -2, -3, -4, -3, -2, -2, -1),
             (-1, -2, -3, -4, -3, -2, -1, 0),
             (-1, -1, -2, -2, -2, -2, -2, -1),
             (-1, -1, -2, -2, -1, -1, -1, 0),
             (-1, -1, -2, -2, -1, 0, 0, 0),
             (-1, 0, -2, -2, -2, -2, -2, -1),
             (0, -1, -1, -2, -2, -2, -2, -1),
             (0, 0, -1, -1, -1, -1, -1, 0))
        """
        import PyNormaliz
        if not self.is_compact():
            raise ValueError('can only enumerate points in a compact polyhedron')
        # Trivial cases: polyhedron with 0 or 1 vertices
        if self.n_vertices() == 0:
            return ()
        if self.n_vertices() == 1:
            v = self.vertices_list()[0]
            try:
                return (vector(ZZ, v),)
            except TypeError:  # vertex not integral
                return ()
        # for small bounding boxes, it is faster to naively iterate over the points of the box
        if threshold > 1:
            box_min, box_max = self.bounding_box(integral_hull=True)
            if box_min is None:
                return ()
            box_points = prod(max_coord-min_coord+1 for min_coord, max_coord in zip(box_min, box_max))
            if  box_points<threshold:
                from sage.geometry.integral_points import rectangular_box_points
                return rectangular_box_points(list(box_min), list(box_max), self)
        # Compute with normaliz
        points = []
        cone = self._normaliz_cone
        assert cone
        for g in PyNormaliz.NmzResult(cone, "ModuleGenerators"):
            assert g[-1] == 1
            points.append(vector(ZZ, g[:-1]))
        return tuple(points)


#########################################################################
class Polyhedron_QQ_normaliz(Polyhedron_normaliz, Polyhedron_QQ):
    r"""
    Polyhedra over `\QQ` with normaliz.

    INPUT:

    - ``Vrep`` -- a list ``[vertices, rays, lines]`` or ``None``
    - ``Hrep`` -- a list ``[ieqs, eqns]`` or ``None``

    EXAMPLES::

        sage: p = Polyhedron(vertices=[(0,0),(1,0),(0,1)],                 # optional - pynormaliz
        ....:                rays=[(1,1)], lines=[],
        ....:                backend='normaliz', base_ring=QQ)
        sage: TestSuite(p).run(skip='_test_pickling')                      # optional - pynormaliz
    """
    pass


#########################################################################
class Polyhedron_ZZ_normaliz(Polyhedron_normaliz, Polyhedron_ZZ):
    r"""
    Polyhedra over `\ZZ` with normaliz.

    INPUT:

    - ``Vrep`` -- a list ``[vertices, rays, lines]`` or ``None``
    - ``Hrep`` -- a list ``[ieqs, eqns]`` or ``None``

    EXAMPLES::

        sage: p = Polyhedron(vertices=[(0,0),(1,0),(0,1)],                 # optional - pynormaliz
        ....:                rays=[(1,1)], lines=[],
        ....:                backend='normaliz', base_ring=ZZ)
        sage: TestSuite(p).run(skip='_test_pickling')                      # optional - pynormaliz
    """
    pass

