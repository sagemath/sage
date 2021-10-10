"""
The Python backend

While slower than specialized C/C++ implementations, the
implementation is general and works with any exact field in Sage that
allows you to define polyhedra.

EXAMPLES::

    sage: p0 = (0, 0)
    sage: p1 = (1, 0)
    sage: p2 = (1/2, AA(3).sqrt()/2)                                            # optional - sage.rings.number_field
    sage: equilateral_triangle = Polyhedron([p0, p1, p2])                       # optional - sage.rings.number_field
    sage: equilateral_triangle.vertices()                                       # optional - sage.rings.number_field
    (A vertex at (0, 0),
     A vertex at (1, 0),
     A vertex at (0.500000000000000?, 0.866025403784439?))
    sage: equilateral_triangle.inequalities()                                   # optional - sage.rings.number_field
    (An inequality (-1, -0.5773502691896258?) x + 1 >= 0,
     An inequality (1, -0.5773502691896258?) x + 0 >= 0,
     An inequality (0, 1.154700538379252?) x + 0 >= 0)
"""
#*****************************************************************************
#       Copyright (C) 2014 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from .base import Polyhedron_base


class Polyhedron_field(Polyhedron_base):
    """
    Polyhedra over all fields supported by Sage

    INPUT:

    - ``Vrep`` -- a list ``[vertices, rays, lines]`` or ``None``.

    - ``Hrep`` -- a list ``[ieqs, eqns]`` or ``None``.

    EXAMPLES::

        sage: p = Polyhedron(vertices=[(0,0),(AA(2).sqrt(),0),(0,AA(3).sqrt())],    # optional - sage.rings.number_field
        ....:                rays=[(1,1)], lines=[], backend='field', base_ring=AA)
        sage: TestSuite(p).run()                                                    # optional - sage.rings.number_field

    TESTS::

        sage: K.<sqrt3> = QuadraticField(3)                                     # optional - sage.rings.number_field
        sage: p = Polyhedron([(0,0), (1,0), (1/2, sqrt3/2)])                    # optional - sage.rings.number_field
        sage: TestSuite(p).run()                                                # optional - sage.rings.number_field

    Check that :trac:`19013` is fixed::

        sage: K.<phi> = NumberField(x^2-x-1, embedding=1.618)                   # optional - sage.rings.number_field
        sage: P1 = Polyhedron([[0,1],[1,1],[1,-phi+1]])                         # optional - sage.rings.number_field
        sage: P2 = Polyhedron(ieqs=[[-1,-phi,0]])                               # optional - sage.rings.number_field
        sage: P1.intersection(P2)                                               # optional - sage.rings.number_field
        The empty polyhedron in (Number Field in phi with defining polynomial x^2 - x - 1 with phi = 1.618033988749895?)^2

    Check that :trac:`28654` is fixed::

        sage: Polyhedron(lines=[[1]], backend='field')
        A 1-dimensional polyhedron in QQ^1 defined as the convex hull of 1 vertex and 1 line
    """
    def _is_zero(self, x):
        """
        Test whether ``x`` is zero.

        INPUT:

        - ``x`` -- a number in the base ring.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: p = Polyhedron([(sqrt(3),sqrt(2))], base_ring=AA)             # optional - sage.rings.number_field
            sage: p._is_zero(0)                                                 # optional - sage.rings.number_field
            True
            sage: p._is_zero(1/100000)                                          # optional - sage.rings.number_field
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

            sage: p = Polyhedron([(sqrt(3),sqrt(2))], base_ring=AA)             # optional - sage.rings.number_field
            sage: p._is_nonneg(1)                                               # optional - sage.rings.number_field
            True
            sage: p._is_nonneg(-1/100000)                                       # optional - sage.rings.number_field
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

            sage: p = Polyhedron([(sqrt(3),sqrt(2))], base_ring=AA)             # optional - sage.rings.number_field
            sage: p._is_positive(1)                                             # optional - sage.rings.number_field
            True
            sage: p._is_positive(0)                                             # optional - sage.rings.number_field
            False
        """
        return x > 0

    def _init_from_Vrepresentation_and_Hrepresentation(self, Vrep, Hrep):
        """
        Construct polyhedron from V-representation and H-representation data.

        See :class:`Polyhedron_base` for a description of ``Vrep`` and ``Hrep``.

        .. WARNING::

            The representation is assumed to be correct.
            It is not checked.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.parent import Polyhedra_field
            sage: from sage.geometry.polyhedron.backend_field import Polyhedron_field
            sage: parent = Polyhedra_field(AA, 1, 'field')                      # optional - sage.rings.number_field
            sage: Vrep = [[[0], [1]], [], []]
            sage: Hrep = [[[0, 1], [1, -1]], []]
            sage: p = Polyhedron_field(parent, Vrep, Hrep,  # indirect doctest  # optional - sage.rings.number_field
            ....:                      Vrep_minimal=True, Hrep_minimal=True)
            sage: p                                                             # optional - sage.rings.number_field
            A 1-dimensional polyhedron in AA^1 defined as the convex hull of 2 vertices
        """
        self._init_Vrepresentation(*Vrep)
        self._init_Hrepresentation(*Hrep)

    def _init_from_Vrepresentation(self, vertices, rays, lines,
                                   minimize=True, verbose=False):
        """
        Construct polyhedron from V-representation data.

        INPUT:

        - ``vertices`` -- list of points. Each point can be specified
           as any iterable container of
           :meth:`~sage.geometry.polyhedron.base.base_ring` elements.

        - ``rays`` -- list of rays. Each ray can be specified as any
          iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements.

        - ``lines`` -- list of lines. Each line can be specified as
          any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements.

        - ``verbose`` -- boolean (default: ``False``). Whether to print
          verbose output for debugging purposes.

        EXAMPLES::

            sage: p = Polyhedron(ambient_dim=2, backend='field')
            sage: from sage.geometry.polyhedron.backend_field import Polyhedron_field
            sage: Polyhedron_field._init_from_Vrepresentation(p, [(0,0)], [], [])
        """
        from sage.geometry.polyhedron.double_description_inhomogeneous import Hrep2Vrep, Vrep2Hrep
        H = Vrep2Hrep(self.base_ring(), self.ambient_dim(), vertices, rays, lines)
        V = Hrep2Vrep(self.base_ring(), self.ambient_dim(),
                      H.inequalities, H.equations)
        self._init_Vrepresentation_backend(V)
        self._init_Hrepresentation_backend(H)

    def _init_from_Hrepresentation(self, ieqs, eqns, minimize=True, verbose=False):
        """
        Construct polyhedron from H-representation data.

        INPUT:

        - ``ieqs`` -- list of inequalities. Each line can be specified
          as any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements.

        - ``eqns`` -- list of equalities. Each line can be specified
          as any iterable container of
          :meth:`~sage.geometry.polyhedron.base.base_ring` elements.

        - ``verbose`` -- boolean (default: ``False``). Whether to print
          verbose output for debugging purposes.

        TESTS::

            sage: p = Polyhedron(ambient_dim=2, backend='field')
            sage: from sage.geometry.polyhedron.backend_field import Polyhedron_field
            sage: Polyhedron_field._init_from_Hrepresentation(p, [(1, 2, 3)], [])
        """
        from sage.geometry.polyhedron.double_description_inhomogeneous import Hrep2Vrep, Vrep2Hrep
        V = Hrep2Vrep(self.base_ring(), self.ambient_dim(), ieqs, eqns)
        H = Vrep2Hrep(self.base_ring(), self.ambient_dim(),
                      V.vertices, V.rays, V.lines)
        self._init_Vrepresentation_backend(V)
        self._init_Hrepresentation_backend(H)

    def _init_Vrepresentation(self, vertices, rays, lines):
        """
        Create the Vrepresentation objects from the given minimal data.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.parent import Polyhedra_field
            sage: from sage.geometry.polyhedron.backend_field import Polyhedron_field
            sage: parent = Polyhedra_field(AA, 1, 'field')                      # optional - sage.rings.number_field
            sage: Vrep = [[[0], [1]], [], []]
            sage: Hrep = [[[0, 1], [1, -1]], []]
            sage: p = Polyhedron_field(parent, Vrep, Hrep,  # indirect doctest  # optional - sage.rings.number_field
            ....:                      Vrep_minimal=True, Hrep_minimal=True)
            sage: p.vertices_list()                                             # optional - sage.rings.number_field
            [[0], [1]]
        """
        self._Vrepresentation = []
        parent = self.parent()
        for v in vertices:
            parent._make_Vertex(self, v)
        for r in rays:
            parent._make_Ray(self, r)
        for l in lines:
            parent._make_Line(self, l)
        self._Vrepresentation = tuple(self._Vrepresentation)

    def _init_Vrepresentation_backend(self, Vrep):
        """
        Create the V-representation objects from the double description.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[(0,1/sqrt(2)),(sqrt(2),0),(4,sqrt(5)/6)],   # optional - sage.rings.number_field
            ....:                base_ring=AA, backend='field')  # indirect doctest
            sage: p.Hrepresentation()                                           # optional - sage.rings.number_field
            (An inequality (-0.1582178750233332?, 1.097777812326429?) x + 0.2237538646678492? >= 0,
             An inequality (-0.1419794359520263?, -1.698172434277148?) x + 1.200789243901438? >= 0,
             An inequality (0.3001973109753594?, 0.600394621950719?) x - 0.4245431085692869? >= 0)
            sage: p.Vrepresentation()                                           # optional - sage.rings.number_field
            (A vertex at (0.?e-15, 0.707106781186548?),
             A vertex at (1.414213562373095?, 0),
             A vertex at (4.000000000000000?, 0.372677996249965?))
        """
        self._init_Vrepresentation(Vrep.vertices, Vrep.rays, Vrep.lines)

    def _init_Hrepresentation(self, inequalities, equations):
        """
        Create the Vrepresentation objects from the given minimal data.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.parent import Polyhedra_field
            sage: from sage.geometry.polyhedron.backend_field import Polyhedron_field
            sage: parent = Polyhedra_field(AA, 1, 'field')                      # optional - sage.rings.number_field
            sage: Vrep = [[[0], [1]], [], []]
            sage: Hrep = [[[0, 1], [1, -1]], []]
            sage: p = Polyhedron_field(parent, Vrep, Hrep,  # indirect doctest  # optional - sage.rings.number_field
            ....:                      Vrep_minimal=True, Hrep_minimal=True)
            sage: p.inequalities_list()                                         # optional - sage.rings.number_field
            [[0, 1], [1, -1]]
        """
        self._Hrepresentation = []
        parent = self.parent()
        for ieq in inequalities:
            parent._make_Inequality(self, ieq)
        for eqn in equations:
            parent._make_Equation(self, eqn)
        self._Hrepresentation = tuple(self._Hrepresentation)

    def _init_Hrepresentation_backend(self, Hrep):
        """
        Create the H-representation objects from the double description.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[(0,1/sqrt(2)),(sqrt(2),0),(4,sqrt(5)/6)],   # optional - sage.rings.number_field
            ....:                base_ring=AA, backend='field')  # indirect doctest
            sage: p.Hrepresentation()                                           # optional - sage.rings.number_field
            (An inequality (-0.1582178750233332?, 1.097777812326429?) x + 0.2237538646678492? >= 0,
             An inequality (-0.1419794359520263?, -1.698172434277148?) x + 1.200789243901438? >= 0,
             An inequality (0.3001973109753594?, 0.600394621950719?) x - 0.4245431085692869? >= 0)
            sage: p.Vrepresentation()                                           # optional - sage.rings.number_field
            (A vertex at (0.?e-15, 0.707106781186548?),
             A vertex at (1.414213562373095?, 0),
             A vertex at (4.000000000000000?, 0.372677996249965?))
        """
        self._init_Hrepresentation(Hrep.inequalities, Hrep.equations)

    def _init_empty_polyhedron(self):
        """
        Initializes an empty polyhedron.

        TESTS::

            sage: empty = Polyhedron(backend='field', base_ring=AA); empty      # optional - sage.rings.number_field
            The empty polyhedron in AA^0
            sage: empty.Vrepresentation()                                       # optional - sage.rings.number_field
            ()
            sage: empty.Hrepresentation()                                       # optional - sage.rings.number_field
            (An equation -1 == 0,)
            sage: Polyhedron(vertices=[], backend='field')
            The empty polyhedron in QQ^0
            sage: Polyhedron(backend='field')._init_empty_polyhedron()
        """
        super(Polyhedron_field, self)._init_empty_polyhedron()
