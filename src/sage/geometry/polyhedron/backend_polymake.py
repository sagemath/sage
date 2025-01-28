"""
The polymake backend for polyhedral computations

.. NOTE::

    This backend requires polymake.
    To install it, type :code:`sage -i polymake` in the terminal.

AUTHORS:

- Matthias Köppe (2017-03): initial version
"""

#*****************************************************************************
#  Copyright (C) 2017 Matthias Köppe <mkoeppe at math.ucdavis.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import itertools

from sage.structure.element import Element

from .base import Polyhedron_base
from .base_QQ import Polyhedron_QQ
from .base_ZZ import Polyhedron_ZZ


#########################################################################
class Polyhedron_polymake(Polyhedron_base):
    """
    Polyhedra with polymake.

    INPUT:

    - ``parent`` -- :class:`~sage.geometry.polyhedron.parent.Polyhedra`
      the parent

    - ``Vrep`` -- list ``[vertices, rays, lines]`` or ``None``; the
      V-representation of the polyhedron. If ``None``, the polyhedron
      is determined by the H-representation.

    - ``Hrep`` -- list ``[ieqs, eqns]`` or ``None``; the
      H-representation of the polyhedron. If ``None``, the polyhedron
      is determined by the V-representation.

    - ``polymake_polytope`` -- a polymake polytope object

    Only one of ``Vrep``, ``Hrep``, or ``polymake_polytope`` can be different
    from ``None``.

    EXAMPLES::

        sage: p = Polyhedron(vertices=[(0,0),(1,0),(0,1)], rays=[(1,1)],   # optional - jupymake
        ....:                lines=[], backend='polymake')
        sage: TestSuite(p).run()                                           # optional - jupymake

    A lower-dimensional affine cone; we test that there are no mysterious
    inequalities coming in from the homogenization::

        sage: P = Polyhedron(vertices=[(1, 1)], rays=[(0, 1)],             # optional - jupymake
        ....:                backend='polymake')
        sage: P.n_inequalities()                                           # optional - jupymake
        1
        sage: P.equations()                                                # optional - jupymake
        (An equation (1, 0) x - 1 == 0,)

    The empty polyhedron::

        sage: Polyhedron(eqns=[[1, 0, 0]], backend='polymake')             # optional - jupymake
        The empty polyhedron in QQ^2

    It can also be obtained differently::

        sage: # optional - jupymake
        sage: P=Polyhedron(ieqs=[[-2, 1, 1], [-3, -1, -1], [-4, 1, -2]],
        ....:              backend='polymake')
        sage: P
        The empty polyhedron in QQ^2
        sage: P.Vrepresentation()
        ()
        sage: P.Hrepresentation()
        (An equation -1 == 0,)

    The full polyhedron::

        sage: Polyhedron(eqns=[[0, 0, 0]], backend='polymake')             # optional - jupymake
        A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 lines
        sage: Polyhedron(ieqs=[[0, 0, 0]], backend='polymake')             # optional - jupymake
        A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 lines


    Quadratic fields work::

        sage: V = polytopes.dodecahedron().vertices_list()                              # needs sage.groups sage.rings.number_field
        sage: Polyhedron(vertices=V, backend='polymake')        # optional - jupymake, needs sage.groups sage.rings.number_field
        A 3-dimensional polyhedron
         in (Number Field in sqrt5 with defining polynomial x^2 - 5
             with sqrt5 = 2.236067977499790?)^3
         defined as the convex hull of 20 vertices

    TESTS:

    Tests copied from various methods in :mod:`sage.geometry.polyhedron.base`::

        sage: p = Polyhedron(vertices = [[1,0,0], [0,1,0], [0,0,1]],       # optional - jupymake
        ....:                backend='polymake')
        sage: p.n_equations()                                              # optional - jupymake
        1
        sage: p.n_inequalities()                                           # optional - jupymake
        3

        sage: p = Polyhedron(vertices = [[t,t^2,t^3] for t in range(6)],   # optional - jupymake
        ....:                backend='polymake')
        sage: p.n_facets()                                                 # optional - jupymake
        8

        sage: p = Polyhedron(vertices = [[1,0],[0,1],[1,1]], rays=[[1,1]], # optional - jupymake
        ....:                backend='polymake')
        sage: p.n_vertices()                                               # optional - jupymake
        2

        sage: p = Polyhedron(vertices = [[1,0],[0,1]], rays=[[1,1]],       # optional - jupymake
        ....:                backend='polymake')
        sage: p.n_rays()                                                   # optional - jupymake
        1

        sage: p = Polyhedron(vertices = [[0,0]], rays=[[0,1],[0,-1]],      # optional - jupymake
        ....:                backend='polymake')
        sage: p.n_lines()                                                  # optional - jupymake
        1
    """

    def _is_zero(self, x):
        """
        Test whether ``x`` is zero.

        INPUT:

        - ``x`` -- a number in the base ring

        OUTPUT: boolean

        EXAMPLES::

            sage: p = Polyhedron([(0,0)], backend='polymake')   # optional - jupymake
            sage: p._is_zero(0)           # optional - jupymake
            True
            sage: p._is_zero(1/100000)    # optional - jupymake
            False
        """
        return x == 0

    def _is_nonneg(self, x):
        """
        Test whether ``x`` is nonnegative.

        INPUT:

        - ``x`` -- a number in the base ring

        OUTPUT: boolean

        EXAMPLES::

            sage: p = Polyhedron([(0,0)], backend='polymake')   # optional - jupymake
            sage: p._is_nonneg(1)         # optional - jupymake
            True
            sage: p._is_nonneg(-1/100000)   # optional - jupymake
            False
        """
        return x >= 0

    def _is_positive(self, x):
        """
        Test whether ``x`` is positive.

        INPUT:

        - ``x`` -- a number in the base ring

        OUTPUT: boolean

        EXAMPLES::

            sage: p = Polyhedron([(0,0)], backend='polymake')   # optional - jupymake
            sage: p._is_positive(1)       # optional - jupymake
            True
            sage: p._is_positive(0)       # optional - jupymake
            False
        """
        return x > 0

    def __init__(self, parent, Vrep, Hrep, polymake_polytope=None, **kwds):
        """
        Initialize the polyhedron.

        See :class:`Polyhedron_polymake` for a description of the input
        data.

        TESTS:

            sage: # optional - jupymake
            sage: p = Polyhedron(backend='polymake')
            sage: TestSuite(p).run()
            sage: p = Polyhedron(vertices=[(1, 1)], rays=[(0, 1)],
            ....:                backend='polymake')
            sage: TestSuite(p).run()

        We skip the Lawrence test because it involves numerically unstable
        floating point arithmetic::

            sage: p = Polyhedron(vertices=[(-1,-1), (1,0), (1,1), (0,1)],  # optional - jupymake
            ....:                backend='polymake')
            sage: TestSuite(p).run(skip='_test_lawrence')            # optional - jupymake

        ::

            sage: # optional - jupymake
            sage: p = Polyhedron(rays=[[1,1]], backend='polymake')
            sage: TestSuite(p).run()
            sage: p = Polyhedron(rays=[[1]], backend='polymake')
            sage: TestSuite(p).run()
            sage: p = Polyhedron(rays=[[1,1,1]], lines=[[1,0,0]], backend='polymake')
            sage: TestSuite(p).run()
            sage: p = Polyhedron(vertices=[[]], backend='polymake')
            sage: TestSuite(p).run()
        """
        if polymake_polytope is not None:
            if Hrep is not None or Vrep is not None:
                raise ValueError("only one of Vrep, Hrep, or polymake_polytope can be different from None")
            Element.__init__(self, parent=parent)
            self._init_from_polymake_polytope(polymake_polytope)
        else:
            Polyhedron_base.__init__(self, parent, Vrep, Hrep, **kwds)

    def _init_from_polymake_polytope(self, polymake_polytope):
        """
        Construct polyhedron from a Polymake polytope object.

        TESTS::

            sage: p = Polyhedron(backend='polymake')                       # optional - jupymake
            sage: from sage.geometry.polyhedron.backend_polymake import Polyhedron_polymake
            sage: Polyhedron_polymake._init_from_Hrepresentation(p, [], [])  # indirect doctest  # optional - jupymake
        """
        self._polymake_polytope = polymake_polytope
        self._init_Vrepresentation_from_polymake()
        self._init_Hrepresentation_from_polymake()

    def _init_from_Vrepresentation(self, vertices, rays, lines, minimize=True, verbose=False):
        r"""
        Construct polyhedron from V-representation data.

        INPUT:

        - ``vertices`` -- list of points; each point can be specified
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

            sage: p = Polyhedron(backend='polymake')                       # optional - jupymake
            sage: from sage.geometry.polyhedron.backend_polymake import Polyhedron_polymake
            sage: Polyhedron_polymake._init_from_Vrepresentation(p, [], [], [])   # optional - jupymake
        """
        from sage.interfaces.polymake import polymake
        data = self._polymake_Vrepresentation_data(vertices, rays, lines)
        polymake_field = polymake(self.base_ring().fraction_field())
        p = polymake.new_object("Polytope<{}>".format(polymake_field), **data)
        self._init_from_polymake_polytope(p)

    def _polymake_Vrepresentation_data(self, vertices, rays, lines, minimal=False):
        r"""
        Compute a dictionary with V-representation input for a polymake Polytope object.

        INPUT:

        - ``vertices`` -- list of points; each point can be specified
          as any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements

        - ``rays`` -- list of rays; each ray can be specified as any
          iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements

        - ``lines`` -- list of lines; each line can be specified as
          any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements

        - ``minimal`` -- boolean (default: ``False``); whether the input is already minimal

        .. WARNING::

            If ``minimal`` the representation is assumed to be correct.
            It is not checked.
        """
        if not minimal:
            return dict(CONE_AMBIENT_DIM=1+self.parent().ambient_dim(),
                        POINTS=(  [ [1] + list(v) for v in vertices ]
                                + [ [0] + list(r) for r in rays ]),
                        INPUT_LINEALITY=[ [0] + list(l) for l in lines ])
        else:
            return dict(CONE_AMBIENT_DIM=1+self.parent().ambient_dim(),
                        VERTICES=(  [ [1] + list(v) for v in vertices ]
                                  + [ [0] + list(r) for r in rays ]),
                        LINEALITY_SPACE=[ [0] + list(l) for l in lines ])

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

            sage: p = Polyhedron(backend='polymake')                       # optional - jupymake
            sage: from sage.geometry.polyhedron.backend_polymake import Polyhedron_polymake
            sage: Polyhedron_polymake._init_from_Hrepresentation(p, [], [])   # optional - jupymake
        """
        from sage.interfaces.polymake import polymake
        data = self._polymake_Hrepresentation_data(ieqs, eqns)
        polymake_field = polymake(self.base_ring().fraction_field())
        p = polymake.new_object("Polytope<{}>".format(polymake_field), **data)
        self._init_from_polymake_polytope(p)

    def _polymake_Hrepresentation_data(self, ieqs, eqns, minimal=False):
        r"""
        Compute a dictionary with H-representation input for a polymake Polytope object.

        INPUT:

        - ``ieqs`` -- list of inequalities; each line can be specified
          as any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements

        - ``eqns`` -- list of equalities; each line can be specified
          as any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements

        - ``minimal`` -- boolean (default: ``False``); whether the input is already minimal

        .. WARNING::

            If ``minimal`` the representation is assumed to be correct.
            It is not checked.
        """
        if ieqs is None:
            ieqs = []
        if eqns is None:
            eqns = []
        # Polymake 3.0r2 and 3.1 crash with a segfault for a test case
        # using QuadraticExtension, when some all-zero inequalities are input.
        # https://forum.polymake.org/viewtopic.php?f=8&t=547
        # Filter them out.
        ieqs = [ list(v) for v in ieqs if not all(self._is_zero(x) for x in v) ]
        # We do a similar filtering for equations.
        # Since Polymake 3.2, we can not give all zero vectors in equations
        eqns = [ list(v) for v in eqns if not all(self._is_zero(x) for x in v) ]
        if not ieqs:
            # Put in one trivial (all-zero) inequality.  This is so that
            # the ambient dimension is set correctly.
            # Since Polymake 3.2, the constant should not be zero.
            ieqs.append([1] + [0]*self.ambient_dim())
        if not minimal:
            return dict(EQUATIONS=eqns,
                        INEQUALITIES=ieqs)
        else:
            return dict(AFFINE_HULL=eqns,
                        FACETS=ieqs)

    def _init_from_Vrepresentation_and_Hrepresentation(self, Vrep, Hrep):
        """
        Construct polyhedron from V-representation and H-representation data.

        See :class:`Polyhedron_base` for a description of ``Vrep`` and ``Hrep``.

        .. WARNING::

            The representation is assumed to be correct.
            It is not checked.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.parent import Polyhedra_polymake
            sage: from sage.geometry.polyhedron.backend_polymake import Polyhedron_polymake
            sage: parent = Polyhedra_polymake(ZZ, 1, 'polymake')
            sage: Vrep = [[[0], [1]], [], []]
            sage: Hrep = [[[0, 1], [1, -1]], []]
            sage: p = Polyhedron_polymake(parent, Vrep, Hrep,  # indirect doctest  # optional - jupymake
            ....:                         Vrep_minimal=True, Hrep_minimal=True)
            sage: p  # optional - jupymake
            A 1-dimensional polyhedron in ZZ^1 defined as the convex hull of 2 vertices
        """
        Vrep = [list(x) for x in Vrep]
        Hrep = [list(x) for x in Hrep]
        p = self._polymake_polytope_from_Vrepresentation_and_Hrepresentation(Vrep, Hrep)
        if p is None:
            self._init_empty_polyhedron()
            return

        self._polymake_polytope = p

        # As the conversion from polymake to sage is slow,
        # we skip it.
        parent = self.parent()
        vertices, rays, lines = Vrep
        inequalities, equations = Hrep
        self._Vrepresentation = []
        self._Hrepresentation = []
        for x in vertices:
            parent._make_Vertex(self, x)
        for x in rays:
            parent._make_Ray(self, x)
        for x in lines:
            parent._make_Line(self, x)
        for x in inequalities:
            parent._make_Inequality(self, x)
        for x in equations:
            parent._make_Equation(self, x)
        self._Vrepresentation = tuple(self._Vrepresentation)
        self._Hrepresentation = tuple(self._Hrepresentation)

    def _polymake_polytope_from_Vrepresentation_and_Hrepresentation(self, Vrep, Hrep):
        if not any(Vrep):
            # The empty polyhedron.
            return

        from sage.interfaces.polymake import polymake
        data = self._polymake_Vrepresentation_data(*Vrep, minimal=True)

        if any(Vrep[1:]):
            from sage.matrix.constructor import Matrix
            polymake_rays = [r for r in data['VERTICES'] if r[0] == 0]
            if Matrix(data['VERTICES']).rank() == Matrix(polymake_rays).rank() + 1:
                # The recession cone is full-dimensional.
                # In this case the homogenized inequalities
                # do not ensure nonnegativy in the last coordinate.
                # In the homogeneous cone the far face is a facet.
                Hrep[0] += [[1] + [0]*self.ambient_dim()]
        data.update(self._polymake_Hrepresentation_data(*Hrep, minimal=True))

        polymake_field = polymake(self.base_ring().fraction_field())
        return polymake.new_object("Polytope<{}>".format(polymake_field), **data)

    def _init_Vrepresentation_from_polymake(self):
        r"""
        Create the Vrepresentation objects from the polymake polyhedron.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[(0,1/2),(2,0),(4,5/6)],  # indirect doctest # optional - jupymake
            ....:                backend='polymake')
            sage: set(p.Hrepresentation())                                 # optional - jupymake
            {An inequality (1, 4) x - 2 >= 0,
             An inequality (1, -12) x + 6 >= 0,
             An inequality (-5, 12) x + 10 >= 0}
            sage: set(p.Vrepresentation())                                 # optional - jupymake
            {A vertex at (0, 1/2), A vertex at (2, 0), A vertex at (4, 5/6)}
        """
        self._Vrepresentation = []
        parent = self.parent()
        p = self._polymake_polytope
        for g in p.VERTICES.sage():
            d = g[0]
            if d == 0:
                parent._make_Ray(self, g[1:])
            elif d == 1:
                parent._make_Vertex(self, g[1:])
            else:
                raise NotImplementedError("Non-normalized vertex encountered: {}".format(g))
        for g in p.LINEALITY_SPACE.sage():
            d = g[0]
            if d == 0:
                parent._make_Line(self, g[1:])
            else:
                raise NotImplementedError("Non-homogeneous line encountered: {}".format(g))
        self._Vrepresentation = tuple(self._Vrepresentation)

    def _init_Hrepresentation_from_polymake(self):
        r"""
        Create the Hrepresentation objects from the polymake polyhedron.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[(0,1/2), (2,0), (4,5/6)],  # indirect doctest # optional - jupymake
            ....:                backend='polymake')
            sage: set(p.Hrepresentation())                                 # optional - jupymake
            {An inequality (1, 4) x - 2 >= 0,
             An inequality (1, -12) x + 6 >= 0,
             An inequality (-5, 12) x + 10 >= 0}
            sage: set(p.Vrepresentation())                                 # optional - jupymake
            {A vertex at (0, 1/2), A vertex at (2, 0), A vertex at (4, 5/6)}
        """
        p = self._polymake_polytope
        if not p.FEASIBLE:
            self._init_empty_polyhedron()
        else:
            self._Hrepresentation = []
            parent = self.parent()
            for g in p.FACETS.sage():
                if all(x == 0 for x in g[1:]):
                    # Ignore vertical inequality
                    pass
                else:
                    parent._make_Inequality(self, g)
            for g in p.AFFINE_HULL.sage():
                parent._make_Equation(self, g)
            self._Hrepresentation = tuple(self._Hrepresentation)

    @classmethod
    def _from_polymake_polytope(cls, parent, polymake_polytope):
        r"""
        Initialize a polyhedron from a polymake Polytope object.

        TESTS::

            sage: from sage.geometry.polyhedron.backend_polymake import Polyhedron_polymake
            sage: from sage.geometry.polyhedron.parent import Polyhedra
            sage: P=Polyhedron(ieqs=[[1, 0, 2], [3, 0, -2], [3, 2, -2]])
            sage: PP = polymake(P)        # optional - jupymake
            sage: parent = Polyhedra(QQ, 2, backend='polymake')   # optional - jupymake
            sage: Q=Polyhedron_polymake._from_polymake_polytope(parent, PP)   # optional - jupymake
        """
        if parent is None:
            from .parent import Polyhedra
            from sage.rings.rational_field import QQ
            from sage.rings.qqbar import AA
            if polymake_polytope.typeof()[0] == 'Polymake::polytope::Polytope__Rational':
                base_ring = QQ
            else:
                from sage.structure.element import coercion_model
                data = [g.sage()
                        for g in itertools.chain(polymake_polytope.VERTICES,
                                                 polymake_polytope.LINEALITY_SPACE,
                                                 polymake_polytope.FACETS,
                                                 polymake_polytope.AFFINE_HULL)]
                if data:
                    base_ring = coercion_model.common_parent(*data).base_ring()
                else:
                    base_ring = QQ
            ambient_dim = polymake_polytope.AMBIENT_DIM().sage()
            parent = Polyhedra(base_ring=base_ring, ambient_dim=ambient_dim, backend='polymake')
        return cls(parent, None, None, polymake_polytope=polymake_polytope)

    def _polymake_(self, polymake):
        """
        Return a polymake "Polytope" object corresponding to ``self``.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[[1, 0], [0, 1], [0, 0]], backend='polymake')   # optional - jupymake
            sage: PP = polymake(P)         # optional - jupymake
            sage: PP.N_VERTICES            # optional - jupymake
            3
        """
        if self._polymake_polytope.parent() is polymake:
            # Same polymake interface, can just return our object
            return self._polymake_polytope
        return super()._polymake_(polymake)

    def __getstate__(self):
        r"""
        Remove the polymake polytope object for pickling.

        TESTS::

            sage: P = polytopes.simplex(backend='polymake')   # optional - jupymake
            sage: P.__getstate__()                            # optional - jupymake
            (Polyhedra in QQ^4,
             {'_Hrepresentation': (An inequality (0, -1, -1, -1) x + 1 >= 0,
               An inequality (0, 1, 0, 0) x + 0 >= 0,
               An inequality (0, 0, 1, 0) x + 0 >= 0,
               An inequality (0, 0, 0, 1) x + 0 >= 0,
               An equation (1, 1, 1, 1) x - 1 == 0),
              '_Vrepresentation': (A vertex at (1, 0, 0, 0),
               A vertex at (0, 1, 0, 0),
               A vertex at (0, 0, 1, 0),
               A vertex at (0, 0, 0, 1)),
              '_pickle_equations': [(-1, 1, 1, 1, 1)],
              '_pickle_inequalities': [(1, 0, -1, -1, -1),
               (0, 0, 1, 0, 0),
               (0, 0, 0, 1, 0),
               (0, 0, 0, 0, 1)],
              '_pickle_lines': [],
              '_pickle_rays': [],
              '_pickle_vertices': [(1, 0, 0, 0),
               (0, 1, 0, 0),
               (0, 0, 1, 0),
               (0, 0, 0, 1)]})
        """
        state = super().__getstate__()
        state = (state[0], state[1].copy())
        # Remove the unpicklable entries.
        if '_polymake_polytope' in state[1]:
            del state[1]['_polymake_polytope']
        state[1]["_pickle_vertices"] = [v._vector for v in self.vertices()]
        state[1]["_pickle_rays"] = [v._vector for v in self.rays()]
        state[1]["_pickle_lines"] = [v._vector for v in self.lines()]
        state[1]["_pickle_inequalities"] = [v._vector for v in self.inequalities()]
        state[1]["_pickle_equations"] = [v._vector for v in self.equations()]
        return state

    def __setstate__(self, state):
        r"""
        Initialize the polymake polytope object after pickling.

        TESTS:

        Test that the obtained cone is valid::

            sage: # optional - jupymake
            sage: from sage.geometry.polyhedron.backend_polymake import Polyhedron_polymake
            sage: P = polytopes.permutahedron(4, backend='polymake')
            sage: P1 = loads(dumps(P))
            sage: P2 = Polyhedron_polymake(P1.parent(), None, None, P1._polymake_polytope)
            sage: P._test_polymake_pickling(other=P2)

            sage: # optional - jupymake
            sage: P = Polyhedron(lines=[[1,0], [0,1]], backend='polymake')
            sage: P1 = loads(dumps(P))
            sage: P2 = Polyhedron_polymake(P1.parent(), None, None, P1._polymake_polytope)
            sage: P._test_polymake_pickling(other=P2)

            sage: P = Polyhedron(backend='polymake')                                        # optional - jupymake
            sage: P1 = loads(dumps(P))                                                      # optional - jupymake
            sage: P._test_polymake_pickling(other=P1)                                       # optional - jupymake

            sage: # optional - jupymake
            sage: P = polytopes.permutahedron(4, backend='polymake') * Polyhedron(lines=[[1]], backend='polymake')
            sage: P1 = loads(dumps(P))
            sage: P2 = Polyhedron_polymake(P1.parent(), None, None, P1._polymake_polytope)
            sage: P._test_polymake_pickling(other=P2)

            sage: # optional - jupymake, needs sage.rings.number_field
            sage: print("Possible output"); P = polytopes.dodecahedron(backend='polymake')
            Possible output...
            sage: P1 = loads(dumps(P))
            sage: P2 = Polyhedron_polymake(P1.parent(), None, None, P1._polymake_polytope)
            sage: P._test_polymake_pickling(other=P2)
        """
        if "_pickle_vertices" in state[1]:
            vertices = state[1].pop("_pickle_vertices")
            rays = state[1].pop("_pickle_rays")
            lines = state[1].pop("_pickle_lines")
            inequalities = state[1].pop("_pickle_inequalities")
            equations = state[1].pop("_pickle_equations")
        else:
            vertices = None

        super().__setstate__(state)

        if vertices is None:
            vertices = self.vertices()
            rays = self.rays()
            lines = self.lines()
            inequalities = self.inequalities()
            equations = self.equations()

        p = self._polymake_polytope_from_Vrepresentation_and_Hrepresentation([vertices, rays, lines], [inequalities, equations])
        if p is not None:
            self._polymake_polytope = p

    def _test_polymake_pickling(self, tester=None, other=None, **options):
        """
        Run tests to see that our polymake pickling/unpickling works.

        INPUT:

        - ``other`` -- a pickling polytope of ``self`` to be tested against

        TESTS::

            sage: polytopes.cross_polytope(3, backend='polymake')._test_polymake_pickling()  # optional - jupymake
        """
        if tester is None:
            tester = self._tester(**options)

        if other is None:
            from sage.misc.persist import loads, dumps
            other = loads(dumps(self))

        tester.assertEqual(self, other)

        if not hasattr(self, '_polymake_polytope'):
            tester.assertFalse(hasattr(other, '_polymake_polytope'))
            return

        P = self._polymake_polytope
        P1 = other._polymake_polytope

        tester.assertEqual(P.F_VECTOR, P1.F_VECTOR)
        tester.assertEqual(P.VERTICES, P1.VERTICES)
        tester.assertEqual(P.LINEALITY_SPACE, P1.LINEALITY_SPACE)
        tester.assertEqual(P.FACETS, P1.FACETS)
        tester.assertEqual(P.AFFINE_HULL, P1.AFFINE_HULL)

#########################################################################


class Polyhedron_QQ_polymake(Polyhedron_polymake, Polyhedron_QQ):
    r"""
    Polyhedra over `\QQ` with polymake.

    INPUT:

    - ``Vrep`` -- list ``[vertices, rays, lines]`` or ``None``
    - ``Hrep`` -- list ``[ieqs, eqns]`` or ``None``

    EXAMPLES::

        sage: p = Polyhedron(vertices=[(0,0),(1,0),(0,1)],                 # optional - jupymake
        ....:                rays=[(1,1)], lines=[],
        ....:                backend='polymake', base_ring=QQ)
        sage: TestSuite(p).run()                                           # optional - jupymake
    """
    pass


#########################################################################
class Polyhedron_ZZ_polymake(Polyhedron_polymake, Polyhedron_ZZ):
    r"""
    Polyhedra over `\ZZ` with polymake.

    INPUT:

    - ``Vrep`` -- list ``[vertices, rays, lines]`` or ``None``
    - ``Hrep`` -- list ``[ieqs, eqns]`` or ``None``

    EXAMPLES::

        sage: p = Polyhedron(vertices=[(0,0),(1,0),(0,1)],                 # optional - jupymake
        ....:                rays=[(1,1)], lines=[],
        ....:                backend='polymake', base_ring=ZZ)
        sage: TestSuite(p).run()                                           # optional - jupymake
    """
    pass
