# -*- coding: utf-8 -*- 
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


from sage.structure.element import Element

from .base import Polyhedron_base
from .base_QQ import Polyhedron_QQ
from .base_ZZ import Polyhedron_ZZ


#########################################################################
class Polyhedron_polymake(Polyhedron_base):
    """
    Polyhedra with polymake

    INPUT:

    - ``parent`` -- :class:`~sage.geometry.polyhedron.parent.Polyhedra`
      the parent

    - ``Vrep`` -- a list ``[vertices, rays, lines]`` or ``None``; the
      V-representation of the polyhedron; if ``None``, the polyhedron
      is determined by the H-representation

    - ``Hrep`` -- a list ``[ieqs, eqns]`` or ``None``; the
      H-representation of the polyhedron; if ``None``, the polyhedron
      is determined by the V-representation

    - ``polymake_polytope`` -- a polymake polytope object

    Only one of ``Vrep``, ``Hrep``, or ``polymake_polytope`` can be different
    from ``None``.

    EXAMPLES::

        sage: p = Polyhedron(vertices=[(0,0),(1,0),(0,1)], rays=[(1,1)],   # optional - polymake
        ....:                lines=[], backend='polymake')
        sage: TestSuite(p).run(skip='_test_pickling')                      # optional - polymake

    A lower-dimensional affine cone; we test that there are no mysterious
    inequalities coming in from the homogenization::

        sage: P = Polyhedron(vertices=[(1, 1)], rays=[(0, 1)],             # optional - polymake
        ....:                backend='polymake')
        sage: P.n_inequalities()                                           # optional - polymake
        1
        sage: P.equations()                                                # optional - polymake
        (An equation (1, 0) x - 1 == 0,)

    The empty polyhedron::

        sage: Polyhedron(eqns=[[1, 0, 0]], backend='polymake')             # optional - polymake
        The empty polyhedron in QQ^2

    It can also be obtained differently::

        sage: P=Polyhedron(ieqs=[[-2, 1, 1], [-3, -1, -1], [-4, 1, -2]],   # optional - polymake
        ....:              backend='polymake')
        sage: P                                                            # optional - polymake
        The empty polyhedron in QQ^2
        sage: P.Vrepresentation()                                          # optional - polymake
        ()
        sage: P.Hrepresentation()                                          # optional - polymake
        (An equation -1 == 0,)

    The full polyhedron::

        sage: Polyhedron(eqns=[[0, 0, 0]], backend='polymake')             # optional - polymake
        A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 lines
        sage: Polyhedron(ieqs=[[0, 0, 0]], backend='polymake')             # optional - polymake
        A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex and 2 lines


    Quadratic fields work::

        sage: V = polytopes.dodecahedron().vertices_list()
        sage: Polyhedron(vertices=V, backend='polymake')                   # optional - polymake
        A 3-dimensional polyhedron in (Number Field in sqrt5 with defining polynomial x^2 - 5)^3 defined as the convex hull of 20 vertices

    TESTS:

    Tests copied from various methods in :mod:`sage.geometry.polyhedron.base`::

        sage: p = Polyhedron(vertices = [[1,0,0], [0,1,0], [0,0,1]],       # optional - polymake
        ....:                backend='polymake')
        sage: p.n_equations()                                              # optional - polymake
        1
        sage: p.n_inequalities()                                           # optional - polymake
        3

        sage: p = Polyhedron(vertices = [[t,t^2,t^3] for t in range(6)],   # optional - polymake
        ....:                backend='polymake')
        sage: p.n_facets()                                                 # optional - polymake
        8

        sage: p = Polyhedron(vertices = [[1,0],[0,1],[1,1]], rays=[[1,1]], # optional - polymake
        ....:                backend='polymake')
        sage: p.n_vertices()                                               # optional - polymake
        2

        sage: p = Polyhedron(vertices = [[1,0],[0,1]], rays=[[1,1]],       # optional - polymake
        ....:                backend='polymake')
        sage: p.n_rays()                                                   # optional - polymake
        1

        sage: p = Polyhedron(vertices = [[0,0]], rays=[[0,1],[0,-1]],      # optional - polymake
        ....:                backend='polymake')
        sage: p.n_lines()                                                  # optional - polymake
        1

    """

    def _is_zero(self, x):
        """
        Test whether ``x`` is zero.

        INPUT:

        - ``x`` -- a number in the base ring.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: p = Polyhedron([(0,0)], backend='polymake')   # optional - polymake
            sage: p._is_zero(0)           # optional - polymake
            True
            sage: p._is_zero(1/100000)    # optional - polymake
            False
        """
        return x == 0

    def _is_nonneg(self, x):
        """
        Test whether ``x`` is nonnegative.

        INPUT:

        - ``x`` -- a number in the base ring.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: p = Polyhedron([(0,0)], backend='polymake')   # optional - polymake
            sage: p._is_nonneg(1)         # optional - polymake
            True
            sage: p._is_nonneg(-1/100000)   # optional - polymake
            False
        """
        return x >= 0

    def _is_positive(self, x):
        """
        Test whether ``x`` is positive.

        INPUT:

        - ``x`` -- a number in the base ring.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: p = Polyhedron([(0,0)], backend='polymake')   # optional - polymake
            sage: p._is_positive(1)       # optional - polymake
            True
            sage: p._is_positive(0)       # optional - polymake
            False
        """
        return x > 0

    def __init__(self, parent, Vrep, Hrep, polymake_polytope=None, **kwds):
        """
        Initializes the polyhedron.

        See :class:`Polyhedron_polymake` for a description of the input
        data.

        TESTS:

        We skip the pickling test because pickling is currently
        not implemented::

            sage: p = Polyhedron(backend='polymake')                 # optional - polymake
            sage: TestSuite(p).run(skip="_test_pickling")            # optional - polymake
            sage: p = Polyhedron(vertices=[(1, 1)], rays=[(0, 1)],   # optional - polymake
            ....:                backend='polymake')
            sage: TestSuite(p).run(skip="_test_pickling")            # optional - polymake
            sage: p = Polyhedron(vertices=[(-1,-1), (1,0), (1,1), (0,1)],  # optional - polymake
            ....:                backend='polymake')
            sage: TestSuite(p).run(skip="_test_pickling")            # optional - polymake
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

            sage: p = Polyhedron(backend='polymake')                       # optional - polymake
            sage: from sage.geometry.polyhedron.backend_polymake import Polyhedron_polymake   # optional - polymake
            sage: Polyhedron_polymake._init_from_Hrepresentation(p, [], [])  # indirect doctest  # optional - polymake
        """
        self._polymake_polytope = polymake_polytope
        self._init_Vrepresentation_from_polymake()
        self._init_Hrepresentation_from_polymake()

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

            sage: p = Polyhedron(backend='polymake')                       # optional - polymake
            sage: from sage.geometry.polyhedron.backend_polymake import Polyhedron_polymake   # optional - polymake
            sage: Polyhedron_polymake._init_from_Vrepresentation(p, [], [], [])   # optional - polymake
        """
        from sage.interfaces.polymake import polymake
        polymake_field = polymake(self.base_ring().fraction_field())
        p = polymake.new_object("Polytope<{}>".format(polymake_field),
                                CONE_AMBIENT_DIM=1+self.parent().ambient_dim(),
                                POINTS=  [ [1] + v for v in vertices ] \
                                       + [ [0] + r for r in rays ],
                                INPUT_LINEALITY=[ [0] + l for l in lines ])
        self._init_from_polymake_polytope(p)

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

            sage: p = Polyhedron(backend='polymake')                       # optional - polymake
            sage: from sage.geometry.polyhedron.backend_polymake import Polyhedron_polymake   # optional - polymake
            sage: Polyhedron_polymake._init_from_Hrepresentation(p, [], [])   # optional - polymake
        """
        from sage.interfaces.polymake import polymake
        if ieqs is None: ieqs = []
        if eqns is None: eqns = []
        # Polymake 3.0r2 and 3.1 crash with a segfault for a test case
        # using QuadraticExtension, when some all-zero inequalities are input.
        # https://forum.polymake.org/viewtopic.php?f=8&t=547
        # Filter them out.
        ieqs = [ v for v in ieqs if not all(self._is_zero(x) for x in v) ]
        # We do a similar filtering for equations.
        # Since Polymake 3.2, we can not give all zero vectors in equations
        eqns = [ v for v in eqns if not all(self._is_zero(x) for x in v) ]
        if not ieqs:
            # Put in one trivial (all-zero) inequality.  This is so that
            # the ambient dimension is set correctly.
            # Since Polymake 3.2, the constant should not be zero.
            ieqs.append([1] + [0]*self.ambient_dim())
        polymake_field = polymake(self.base_ring().fraction_field())
        p = polymake.new_object("Polytope<{}>".format(polymake_field),
                                EQUATIONS=eqns,
                                INEQUALITIES=ieqs)
        self._init_from_polymake_polytope(p)

    def _init_Vrepresentation_from_polymake(self):
        r"""
        Create the Vrepresentation objects from the polymake polyhedron.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[(0,1/2),(2,0),(4,5/6)],  # indirect doctest # optional - polymake
            ....:                backend='polymake')
            sage: set(p.Hrepresentation())                                 # optional - polymake
            {An inequality (1, 4) x - 2 >= 0,
             An inequality (1, -12) x + 6 >= 0,
             An inequality (-5, 12) x + 10 >= 0}
            sage: set(p.Vrepresentation())                                 # optional - polymake
            {A vertex at (0, 1/2), A vertex at (2, 0), A vertex at (4, 5/6)}

        """
        self._Vrepresentation = []
        parent = self.parent()
        p = self._polymake_polytope
        for g in p.VERTICES:
            g = g.sage()
            d = g[0]
            if d == 0:
                parent._make_Ray(self, g[1:])
            elif d == 1:
                parent._make_Vertex(self, g[1:])
            else:
                raise NotImplementedError("Non-normalized vertex encountered: {}".format(g))
        for g in p.LINEALITY_SPACE:
            g = g.sage()
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

            sage: p = Polyhedron(vertices=[(0,1/2), (2,0), (4,5/6)],  # indirect doctest # optional - polymake
            ....:                backend='polymake')
            sage: set(p.Hrepresentation())                                 # optional - polymake
            {An inequality (1, 4) x - 2 >= 0,
             An inequality (1, -12) x + 6 >= 0,
             An inequality (-5, 12) x + 10 >= 0}
            sage: set(p.Vrepresentation())                                 # optional - polymake
            {A vertex at (0, 1/2), A vertex at (2, 0), A vertex at (4, 5/6)}

        """
        p = self._polymake_polytope
        if not p.FEASIBLE:
            self._init_empty_polyhedron()
        else:
            self._Hrepresentation = []
            parent = self.parent()
            for g in p.FACETS:
                if all(x==0 for x in g[1:]):
                    # Ignore vertical inequality
                    pass
                else:
                    parent._make_Inequality(self, g.sage())
            for g in p.AFFINE_HULL:
                parent._make_Equation(self, g.sage())
            self._Hrepresentation = tuple(self._Hrepresentation)

    @classmethod
    def _from_polymake_polytope(cls, parent, polymake_polytope):
        r"""
        Initializes a polyhedron from a polymake Polytope object.

        TESTS::

            sage: from sage.geometry.polyhedron.backend_polymake import Polyhedron_polymake
            sage: from sage.geometry.polyhedron.parent import Polyhedra
            sage: P=Polyhedron(ieqs=[[1, 0, 2], [3, 0, -2], [3, 2, -2]])
            sage: PP = polymake(P)        # optional - polymake
            sage: parent = Polyhedra(QQ, 2, backend='polymake')   # optional - polymake
            sage: Q=Polyhedron_polymake._from_polymake_polytope(parent, PP)   # optional - polymake
        """
        return cls(parent, None, None, polymake_polytope=polymake_polytope)

    def _polymake_(self, polymake):
        """
        Return a polymake "Polytope" object corresponding to ``self``.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[[1, 0], [0, 1], [0, 0]], backend='polymake')   # optional - polymake
            sage: PP = polymake(P)         # optional - polymake
            sage: PP.N_VERTICES            # optional - polymake
            3
        """
        if self._polymake_polytope.parent() is polymake:
            # Same polymake interface, can just return our object
            return self._polymake_polytope
        else:
            return super(Polyhedron_polymake, self)._polymake_(polymake)


#########################################################################
class Polyhedron_QQ_polymake(Polyhedron_polymake, Polyhedron_QQ):
    r"""
    Polyhedra over `\QQ` with polymake.

    INPUT:

    - ``Vrep`` -- a list ``[vertices, rays, lines]`` or ``None``
    - ``Hrep`` -- a list ``[ieqs, eqns]`` or ``None``

    EXAMPLES::

        sage: p = Polyhedron(vertices=[(0,0),(1,0),(0,1)],                 # optional - polymake
        ....:                rays=[(1,1)], lines=[],
        ....:                backend='polymake', base_ring=QQ)
        sage: TestSuite(p).run(skip='_test_pickling')                      # optional - polymake
    """
    pass


#########################################################################
class Polyhedron_ZZ_polymake(Polyhedron_polymake, Polyhedron_ZZ):
    r"""
    Polyhedra over `\ZZ` with polymake.

    INPUT:

    - ``Vrep`` -- a list ``[vertices, rays, lines]`` or ``None``
    - ``Hrep`` -- a list ``[ieqs, eqns]`` or ``None``

    EXAMPLES::

        sage: p = Polyhedron(vertices=[(0,0),(1,0),(0,1)],                 # optional - polymake
        ....:                rays=[(1,1)], lines=[],
        ....:                backend='polymake', base_ring=ZZ)
        sage: TestSuite(p).run(skip='_test_pickling')                      # optional - polymake
    """
    pass

