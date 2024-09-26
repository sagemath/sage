r"""
Base class for polyhedra: Miscellaneous methods
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

from sage.misc.cachefunc import cached_method

import sage.rings.abc
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector

from .base7 import Polyhedron_base7

#########################################################################
# Notes if you want to implement your own backend:
#
#  * derive from Polyhedron_base
#
#  * you must implement _init_from_Vrepresentation and
#    _init_from_Hrepresentation
#
#  * You might want to override _init_empty_polyhedron
#
#  * You may implement _init_from_Vrepresentation_and_Hrepresentation
#
#  * You can of course also override any other method for which you
#    have a faster implementation.
#########################################################################


#########################################################################
def is_Polyhedron(X):
    """
    Test whether ``X`` is a Polyhedron.

    INPUT:

    - ``X`` -- anything

    OUTPUT: boolean

    EXAMPLES::

        sage: p = polytopes.hypercube(2)
        sage: from sage.geometry.polyhedron.base import is_Polyhedron
        sage: is_Polyhedron(p)
        doctest:warning...
        DeprecationWarning: is_Polyhedron is deprecated, use isinstance instead
        See https://github.com/sagemath/sage/issues/34307 for details.
        True
        sage: is_Polyhedron(123456)
        False
    """
    from sage.misc.superseded import deprecation
    deprecation(34307, "is_Polyhedron is deprecated, use isinstance instead")
    return isinstance(X, Polyhedron_base)


#########################################################################
class Polyhedron_base(Polyhedron_base7):
    """
    Base class for Polyhedron objects.

    INPUT:

    - ``parent`` -- the parent, an instance of
      :class:`~sage.geometry.polyhedron.parent.Polyhedra`

    - ``Vrep`` -- list ``[vertices, rays, lines]`` or ``None``. The
      V-representation of the polyhedron; if ``None``, the polyhedron
      is determined by the H-representation

    - ``Hrep`` -- list ``[ieqs, eqns]`` or ``None``. The
      H-representation of the polyhedron; if ``None``, the polyhedron
      is determined by the V-representation

    - ``Vrep_minimal`` -- (optional) see below

    - ``Hrep_minimal`` -- (optional) see below

    - ``pref_rep`` -- string (default: ``None``);
      one of ``Vrep`` or ``Hrep`` to pick this in case the backend
      cannot initialize from complete double description

    - ``mutable`` -- ignored

    If both ``Vrep`` and ``Hrep`` are provided, then
    ``Vrep_minimal`` and ``Hrep_minimal`` must be set to ``True``.

    TESTS::

        sage: p = Polyhedron()
        sage: TestSuite(p).run()

    ::

        sage: p = Polyhedron(vertices=[(1,0), (0,1)], rays=[(1,1)], base_ring=ZZ)
        sage: TestSuite(p).run()

    ::

        sage: p = polytopes.flow_polytope(digraphs.DeBruijn(3,2))                       # needs sage.combinat sage.graphs
        sage: TestSuite(p).run()

    ::

        sage: TestSuite(Polyhedron([[]])).run()
        sage: TestSuite(Polyhedron([[0]])).run()
        sage: TestSuite(Polyhedron([[1]])).run()

    ::

        sage: P3 = polytopes.permutahedron(3)
        sage: P = P3 * Polyhedron(rays=[[0,0,1], [0,1,1], [1,2,3]])
        sage: TestSuite(P).run()

    ::

        sage: P = P3 * Polyhedron(rays=[[0,0,1], [0,1,1]], lines=[[1,0,0]])
        sage: TestSuite(P).run()

    ::

        sage: M = random_matrix(ZZ, 5, 5, distribution='uniform')
        sage: while True:
        ....:     M = random_matrix(ZZ, 5, 5, distribution='uniform')
        ....:     if M.rank() != 5:
        ....:         break
        ....:
        sage: P = Polyhedron(M)
        sage: TestSuite(P).run()
    """

    def _test_basic_properties(self, tester=None, **options):
        """
        Run some basic tests to see, that some general assertion on polyhedra hold.

        TESTS::

            sage: polytopes.cross_polytope(3)._test_basic_properties()
        """
        from .constructor import Polyhedron

        if tester is None:
            tester = self._tester(**options)

        tester.assertEqual(self.n_vertices() + self.n_rays() + self.n_lines(), self.n_Vrepresentation())
        tester.assertEqual(self.n_inequalities() + self.n_equations(), self.n_Hrepresentation())
        if self.n_vertices():
            # Depending on the backend, this does not hold for the empty polyhedron.
            tester.assertEqual(self.dim() + self.n_equations(), self.ambient_dim())

        tester.assertTrue(all(len(v[::]) == self.ambient_dim() for v in self.Vrep_generator()))
        tester.assertTrue(all(len(h[::]) == self.ambient_dim() + 1 for h in self.Hrep_generator()))

        if self.n_vertices() + self.n_rays() < 40:
            tester.assertEqual(self, Polyhedron(vertices=self.vertices(), rays=self.rays(), lines=self.lines(), ambient_dim=self.ambient_dim()))
        if self.n_inequalities() < 40:
            tester.assertEqual(self, Polyhedron(ieqs=self.inequalities(), eqns=self.equations(), ambient_dim=self.ambient_dim()))

    def to_linear_program(self, solver=None, return_variable=False, base_ring=None):
        r"""
        Return a linear optimization problem over the polyhedron in the form of
        a :class:`MixedIntegerLinearProgram`.

        INPUT:

        - ``solver`` -- select a solver (MIP backend). See the documentation
          of for :class:`MixedIntegerLinearProgram`. Set to ``None`` by default.

        - ``return_variable`` -- boolean (default: ``False``); if ``True``, return a tuple
          ``(p, x)``, where ``p`` is the :class:`MixedIntegerLinearProgram` object
          and ``x`` is the vector-valued MIP variable in this problem, indexed
          from 0.  If ``False``, only return ``p``.

        - ``base_ring`` -- select a field over which the linear program should be
          set up.  Use ``RDF`` to request a fast inexact (floating point) solver
          even if ``self`` is exact.

        Note that the :class:`MixedIntegerLinearProgram` object will have the
        null function as an objective to be maximized.

        .. SEEALSO::

            :meth:`~MixedIntegerLinearProgram.polyhedron` -- return the
            polyhedron associated with a :class:`MixedIntegerLinearProgram`
            object.

        EXAMPLES:

        Exact rational linear program::

            sage: p = polytopes.cube()
            sage: p.to_linear_program()
            Linear Program (no objective, 3 variables, 6 constraints)
            sage: lp, x = p.to_linear_program(return_variable=True)
            sage: lp.set_objective(2*x[0] + 1*x[1] + 39*x[2])
            sage: lp.solve()
            42
            sage: lp.get_values(x[0], x[1], x[2])
            [1, 1, 1]

        Floating-point linear program::

            sage: lp, x = p.to_linear_program(return_variable=True, base_ring=RDF)
            sage: lp.set_objective(2*x[0] + 1*x[1] + 39*x[2])
            sage: lp.solve()
            42.0

        Irrational algebraic linear program over an embedded number field::

            sage: # needs sage.groups sage.rings.number_field
            sage: p = polytopes.icosahedron()
            sage: lp, x = p.to_linear_program(return_variable=True)
            sage: lp.set_objective(x[0] + x[1] + x[2])
            sage: lp.solve()
            1/4*sqrt5 + 3/4

        Same example with floating point::

            sage: # needs sage.groups sage.rings.number_field
            sage: lp, x = p.to_linear_program(return_variable=True, base_ring=RDF)
            sage: lp.set_objective(x[0] + x[1] + x[2])
            sage: lp.solve()                                               # tol 1e-5
            1.3090169943749475

        Same example with a specific floating point solver::

            sage: # needs sage.groups sage.rings.number_field
            sage: lp, x = p.to_linear_program(return_variable=True, solver='GLPK')
            sage: lp.set_objective(x[0] + x[1] + x[2])
            sage: lp.solve()                                               # tol 1e-8
            1.3090169943749475

        Irrational algebraic linear program over `AA`::

            sage: # needs sage.groups sage.rings.number_field
            sage: p = polytopes.icosahedron(base_ring=AA)
            sage: lp, x = p.to_linear_program(return_variable=True)
            sage: lp.set_objective(x[0] + x[1] + x[2])
            sage: lp.solve()                    # long time
            1.309016994374948?

        TESTS::

            sage: p = polytopes.flow_polytope(digraphs.DeBruijn(3,2)); p                # needs sage.combinat sage.graphs
            A 19-dimensional polyhedron in QQ^27
             defined as the convex hull of 1 vertex and 148 rays
            sage: p.to_linear_program().polyhedron() == p
            True

            sage: p = polytopes.icosahedron()                                           # needs sage.groups sage.rings.number_field
            sage: p.to_linear_program(solver='PPL')                                     # needs sage.groups sage.rings.number_field
            Traceback (most recent call last):
            ...
            TypeError: The PPL backend only supports rational data.

        Test that equations are handled correctly (:issue:`24154`)::

            sage: p = Polyhedron(vertices=[[19]])
            sage: lp, x = p.to_linear_program(return_variable=True)
            sage: lp.set_objective(x[0])
            sage: lp.solve()
            19
        """
        if base_ring is None:
            base_ring = self.base_ring()
        base_ring = base_ring.fraction_field()
        from sage.numerical.mip import MixedIntegerLinearProgram
        p = MixedIntegerLinearProgram(solver=solver, base_ring=base_ring)
        x = p.new_variable(real=True, nonnegative=False)

        for ineqn in self.inequalities_list():
            b = -ineqn.pop(0)
            p.add_constraint(p.sum([x[i] * ineqn[i] for i in range(len(ineqn))]) >= b)

        for eqn in self.equations_list():
            b = -eqn.pop(0)
            p.add_constraint(p.sum([x[i] * eqn[i] for i in range(len(eqn))]) == b)

        if return_variable:
            return p, x
        else:
            return p

    def boundary_complex(self):
        """
        Return the simplicial complex given by the boundary faces of ``self``,
        if it is simplicial.

        OUTPUT:

        A (spherical) simplicial complex

        EXAMPLES:

        The boundary complex of the octahedron::

            sage: # needs sage.graphs
            sage: oc = polytopes.octahedron()
            sage: sc_oc = oc.boundary_complex()
            sage: fl_oc = oc.face_lattice()                                             # needs sage.combinat
            sage: fl_sc = sc_oc.face_poset()                                            # needs sage.combinat
            sage: [len(x) for x in fl_oc.level_sets()]                                  # needs sage.combinat
            [1, 6, 12, 8, 1]
            sage: [len(x) for x in fl_sc.level_sets()]                                  # needs sage.combinat
            [6, 12, 8]
            sage: sc_oc.euler_characteristic()
            2
            sage: sc_oc.homology()
            {0: 0, 1: 0, 2: Z}

        The polyhedron should be simplicial::

            sage: c = polytopes.cube()
            sage: c.boundary_complex()
            Traceback (most recent call last):
            ...
            NotImplementedError: this function is only implemented for simplicial polytopes

        TESTS::

            sage: p = Polyhedron(rays=[[1,1]])
            sage: p.boundary_complex()
            Traceback (most recent call last):
            ...
            ValueError: self should be compact
        """
        if not self.is_compact():
            raise ValueError("self should be compact")

        if self.is_simplicial():
            from sage.topology.simplicial_complex import SimplicialComplex
            inc_mat_cols = self.incidence_matrix().columns()
            ineq_indices = [inc_mat_cols[i].nonzero_positions()
                            for i in range(self.n_Hrepresentation())
                            if self.Hrepresentation()[i].is_inequality()]
            return SimplicialComplex(ineq_indices, maximality_check=False)
        else:
            raise NotImplementedError("this function is only implemented for simplicial polytopes")

    @cached_method
    def center(self):
        """
        Return the average of the vertices.

        .. SEEALSO::

            :meth:`sage.geometry.polyhedron.base1.Polyhedron_base1.representative_point`.

        OUTPUT:

        The center of the polyhedron. All rays and lines are
        ignored. Raises a :exc:`ZeroDivisionError` for the empty
        polytope.

        EXAMPLES::

            sage: p = polytopes.hypercube(3)
            sage: p = p + vector([1,0,0])
            sage: p.center()
            (1, 0, 0)
        """
        if self.dim() == 0:
            return self.vertices()[0].vector()
        else:
            vertex_sum = vector(self.base_ring(), [0] * self.ambient_dim())
            for v in self.vertex_generator():
                vertex_sum += v.vector()
            vertex_sum.set_immutable()
            return vertex_sum / self.n_vertices()

    @cached_method
    def radius_square(self):
        """
        Return the square of the maximal distance from the
        :meth:`center` to a vertex. All rays and lines are ignored.

        OUTPUT:

        The square of the radius, which is in
        :meth:`~sage.geometry.polyhedron.base0.Polyhedron_base0.base_ring`.

        EXAMPLES::

            sage: p = polytopes.permutahedron(4, project = False)
            sage: p.radius_square()
            5
        """
        vertices = [v.vector() - self.center() for v in self.vertex_generator()]
        return max(v.dot_product(v) for v in vertices)

    def radius(self):
        """
        Return the maximal distance from the center to a vertex. All
        rays and lines are ignored.

        OUTPUT:

        The radius for a rational polyhedron is, in general, not
        rational.  use :meth:`radius_square` if you need a rational
        distance measure.

        EXAMPLES::

            sage: p = polytopes.hypercube(4)
            sage: p.radius()
            2
        """
        return self.radius_square().sqrt()

    def is_inscribed(self, certificate=False):
        """
        This function tests whether the vertices of the polyhedron are
        inscribed on a sphere.

        The polyhedron is expected to be compact and full-dimensional.
        A full-dimensional compact polytope is inscribed if there exists
        a point in space which is equidistant to all its vertices.

        ALGORITHM:

        The function first computes the circumsphere of a full-dimensional
        simplex with vertices of ``self``. It is found by lifting the points on a
        paraboloid to find the hyperplane on which the circumsphere is lifted.
        Then, it checks if all other vertices are equidistant to the
        circumcenter of that simplex.

        INPUT:

        - ``certificate`` -- boolean (default: ``False``); specifies whether to
          return the circumcenter, if found

        OUTPUT: if ``certificate`` is true, returns a tuple containing:

        1. Boolean.
        2. The circumcenter of the polytope or None.

        If ``certificate`` is false:

        - a Boolean.

        EXAMPLES::

            sage: q = Polyhedron(vertices=[[1,1,1,1],[-1,-1,1,1],[1,-1,-1,1],
            ....:                            [-1,1,-1,1],[1,1,1,-1],[-1,-1,1,-1],
            ....:                            [1,-1,-1,-1],[-1,1,-1,-1],[0,0,10/13,-24/13],
            ....:                            [0,0,-10/13,-24/13]])
            sage: q.is_inscribed(certificate=True)
            (True, (0, 0, 0, 0))

            sage: cube = polytopes.cube()
            sage: cube.is_inscribed()
            True

            sage: translated_cube = Polyhedron(vertices=[v.vector() + vector([1,2,3])
            ....:                                        for v in cube.vertices()])
            sage: translated_cube.is_inscribed(certificate=True)
            (True, (1, 2, 3))

            sage: truncated_cube = cube.face_truncation(cube.faces(0)[0])
            sage: truncated_cube.is_inscribed()
            False

        The method is not implemented for non-full-dimensional polytope or
        unbounded polyhedra::

            sage: square = Polyhedron(vertices=[[1,0,0],[0,1,0],[1,1,0],[0,0,0]])
            sage: square.is_inscribed()
            Traceback (most recent call last):
            ...
            NotImplementedError: this function is implemented for full-dimensional polyhedra only

            sage: p = Polyhedron(vertices=[(0,0)],rays=[(1,0),(0,1)])
            sage: p.is_inscribed()
            Traceback (most recent call last):
            ...
            NotImplementedError: this function is not implemented for unbounded polyhedra

        TESTS:

        We check that :issue:`28464` is fixed::

            sage: P = Polyhedron(vertices=[(-130658298093891402635075/416049251842505144482473,
            ....: 177469511761879509172000/1248147755527515433447419,
            ....: 485550543257132133136169/2496295511055030866894838,
            ....: 2010744967797898733758669/2496295511055030866894838),
            ....: (-146945725603929909850/706333405676769433081,
            ....: -84939725782618445000/706333405676769433081,
            ....: 560600045283000988081/1412666811353538866162,
            ....: 969778382942371268081/1412666811353538866162),
            ....: (-46275018824497300/140422338198040641,
            ....: -5747688262110000/46807446066013547, 1939357556329/7033601552658,
            ....: 1939357556329/7033601552658), (-17300/59929, -10000/59929, 39929/119858,
            ....: 39929/119858), (-4700/32209, -10000/32209, 12209/64418, 12209/64418),
            ....: (QQ(0), QQ(0), QQ(0), QQ(1)), (QQ(0), QQ(0), 1/2, 1/2), (300/10027,
            ....: -10000/30081, 10081/60162, 10081/60162), (112393975400/1900567733649,
            ....: 117311600000/633522577883, 43678681/95197362, 43678681/95197362),
            ....: (6109749955400/133380598418321, 37106807920000/133380598418321,
            ....: 2677964249/6680888498, 2677964249/6680888498),
            ....: (29197890764005600/402876806828660641,
            ....: -2150510776960000/402876806828660641,
            ....: 398575785274740641/805753613657321282,
            ....: 398575785274740641/805753613657321282),
            ....: (5576946899441759759983005325/110078073300232813237456943251,
            ....: -29071211718677797926570478000/110078073300232813237456943251,
            ....: 59439312069347378584317232001/220156146600465626474913886502,
            ....: 181346577228466312205473034501/220156146600465626474913886502),
            ....: (150040732779124914266530235300/6774574358246204311268446913881,
            ....: -2813827375989039189507000218000/6774574358246204311268446913881,
            ....: 1260217414021285074925933133881/13549148716492408622536893827762,
            ....: 3232518047094242684574253773881/13549148716492408622536893827762),
            ....: (3816349407976279597850158016285000/88842127448735433741180809504357161,
            ....: 27965821247423216557301387453968000/88842127448735433741180809504357161,
            ....: 68546256000224819256028677086357161/177684254897470867482361619008714322,
            ....: 86062257922545755787315412690197161/177684254897470867482361619008714322)])
            sage: P.is_inscribed()
            True

            sage: P = Polyhedron(vertices=[[0, -1, 0, 0],
            ....:                          [0, 0, -1, 0],
            ....:                          [0, 0, 0, -1],
            ....:                          [0, 0, +1, 0],
            ....:                          [0, 0, 0, +1],
            ....:                          [+1, 0, 0, 0]])
            sage: P.is_inscribed()
            True

        We check that :issue:`29125` is fixed::

            sage: P = Polyhedron(vertices=[[-2,-1], [-2,1], [0,-1], [0,1]], backend='field')
            sage: P.is_inscribed()
            True
            sage: V = P.Vrepresentation()
            sage: H = P.Hrepresentation()
            sage: parent = P.parent()
            sage: for V1 in Permutations(V):
            ....:     P1 = parent._element_constructor_(
            ....:         [V1, [], []], [H, []], Vrep_minimal=True, Hrep_minimal=True)
            ....:     assert P1.is_inscribed()
        """

        if not self.is_compact():
            raise NotImplementedError("this function is not implemented for unbounded polyhedra")

        if not self.is_full_dimensional():
            raise NotImplementedError("this function is implemented for full-dimensional polyhedra only")

        dimension = self.dimension()
        vertices = self.vertices()

        # We obtain vertices that are an affine basis of the affine hull.
        affine_basis = self.an_affine_basis()
        raw_data = []
        for vertex in affine_basis:
            vertex_vector = vertex.vector()
            raw_data += [[sum(i**2 for i in vertex_vector)] +
                         [i for i in vertex_vector] + [1]]
        matrix_data = matrix(raw_data)

        # The determinant "a" should not be zero because
        # the vertices in ``affine_basis`` are an affine basis.
        a = matrix_data.matrix_from_columns(range(1, dimension+2)).determinant()

        minors = [(-1)**(i)*matrix_data.matrix_from_columns([j for j in range(dimension+2) if j != i]).determinant()
                  for i in range(1, dimension+1)]
        c = (-1)**(dimension+1)*matrix_data.matrix_from_columns(range(dimension+1)).determinant()

        circumcenter = vector([minors[i]/(2*a) for i in range(dimension)])
        squared_circumradius = (sum(m**2 for m in minors) - 4 * a * c) / (4*a**2)

        # Checking if the circumcenter has the correct sign
        if not all(sum(i**2 for i in v.vector() - circumcenter) == squared_circumradius
                   for v in vertices if v in affine_basis):
            circumcenter = - circumcenter

        is_inscribed = all(sum(i**2 for i in v.vector() - circumcenter) == squared_circumradius
                           for v in vertices if v not in affine_basis)

        if certificate:
            if is_inscribed:
                return (True, circumcenter)
            else:
                return (False, None)
        else:
            return is_inscribed

    def hyperplane_arrangement(self):
        """
        Return the hyperplane arrangement defined by the equations and
        inequalities.

        OUTPUT:

        A :class:`hyperplane arrangement
        <sage.geometry.hyperplane_arrangement.arrangement.HyperplaneArrangementElement>`
        consisting of the hyperplanes defined by the
        :meth:`~sage.geometry.polyhedron.base0.Polyhedron_base0.Hrepresentation`.
        If the polytope is full-dimensional, this is the hyperplane
        arrangement spanned by the facets of the polyhedron.

        EXAMPLES::

            sage: p = polytopes.hypercube(2)
            sage: p.hyperplane_arrangement()
            Arrangement <-t0 + 1 | -t1 + 1 | t1 + 1 | t0 + 1>
        """
        names = tuple('t' + str(i) for i in range(self.ambient_dim()))
        from sage.geometry.hyperplane_arrangement.arrangement import HyperplaneArrangements
        field = self.base_ring().fraction_field()
        H = HyperplaneArrangements(field, names)
        return H(self)

    @cached_method
    def normal_fan(self, direction='inner'):
        r"""
        Return the normal fan of a compact full-dimensional rational polyhedron.

        This returns the inner normal fan of ``self``. For the outer normal fan,
        use ``direction='outer'``.

        INPUT:

        - ``direction`` -- either ``'inner'`` (default) or ``'outer'``; if
          set to ``'inner'``, use the inner normal vectors to span the cones of
          the fan, if set to ``'outer'``, use the outer normal vectors.

        OUTPUT:

        A complete fan of the ambient space as a
        :class:`~sage.geometry.fan.RationalPolyhedralFan`.

        .. SEEALSO::

            :meth:`face_fan`.

        EXAMPLES::

            sage: S = Polyhedron(vertices=[[0, 0], [1, 0], [0, 1]])
            sage: S.normal_fan()
            Rational polyhedral fan in 2-d lattice N

            sage: C = polytopes.hypercube(4)
            sage: NF = C.normal_fan(); NF
            Rational polyhedral fan in 4-d lattice N

        Currently, it is only possible to get the normal fan of a bounded rational polytope::

            sage: P = Polyhedron(rays=[[1, 0], [0, 1]])
            sage: P.normal_fan()
            Traceback (most recent call last):
            ...
            NotImplementedError: the normal fan is only supported for polytopes (compact polyhedra).

            sage: Q = Polyhedron(vertices=[[1, 0, 0], [0, 1, 0], [0, 0, 1]])
            sage: Q.normal_fan()
            Traceback (most recent call last):
            ...
            ValueError: the normal fan is only defined for full-dimensional polytopes

            sage: R = Polyhedron(vertices=[[0, 0],                                      # needs sage.rings.number_field sage.symbolic
            ....:                          [AA(sqrt(2)), 0],
            ....:                          [0, AA(sqrt(2))]])
            sage: R.normal_fan()                                                        # needs sage.rings.number_field sage.symbolic
            Traceback (most recent call last):
            ...
            NotImplementedError: normal fan handles only polytopes over the rationals

            sage: P = Polyhedron(vertices=[[0,0], [2,0], [0,2], [2,1], [1,2]])
            sage: P.normal_fan(direction=None)
            Traceback (most recent call last):
            ...
            TypeError: the direction should be 'inner' or 'outer'

            sage: inner_nf = P.normal_fan()
            sage: inner_nf.rays()
            N( 1,  0),
            N( 0, -1),
            N( 0,  1),
            N(-1,  0),
            N(-1, -1)
            in 2-d lattice N

            sage: outer_nf = P.normal_fan(direction='outer')
            sage: outer_nf.rays()
            N( 1,  0),
            N( 1,  1),
            N( 0,  1),
            N(-1,  0),
            N( 0, -1)
            in 2-d lattice N

        REFERENCES:

        For more information, see Chapter 7 of [Zie2007]_.
        """
        from sage.geometry.fan import NormalFan

        if not QQ.has_coerce_map_from(self.base_ring()):
            raise NotImplementedError('normal fan handles only polytopes over the rationals')
        if direction == 'inner':
            return NormalFan(self)
        elif direction == 'outer':
            return NormalFan(-self)
        else:
            raise TypeError("the direction should be 'inner' or 'outer'")

    @cached_method
    def face_fan(self):
        r"""
        Return the face fan of a compact rational polyhedron.

        OUTPUT:

        A fan of the ambient space as a
        :class:`~sage.geometry.fan.RationalPolyhedralFan`.

        .. SEEALSO::

            :meth:`normal_fan`.

        EXAMPLES::

            sage: T = polytopes.cuboctahedron()
            sage: T.face_fan()
            Rational polyhedral fan in 3-d lattice M

        The polytope should contain the origin in the interior::

            sage: P = Polyhedron(vertices=[[1/2, 1], [1, 1/2]])
            sage: P.face_fan()
            Traceback (most recent call last):
            ...
            ValueError: face fans are defined only for polytopes
            containing the origin as an interior point!

            sage: Q = Polyhedron(vertices=[[-1, 1/2], [1, -1/2]])
            sage: Q.contains([0,0])
            True
            sage: FF = Q.face_fan(); FF
            Rational polyhedral fan in 2-d lattice M

        The polytope has to have rational coordinates::

            sage: S = polytopes.dodecahedron()                                          # needs sage.groups sage.rings.number_field
            sage: S.face_fan()                                                          # needs sage.groups sage.rings.number_field
            Traceback (most recent call last):
            ...
            NotImplementedError: face fan handles only polytopes over the rationals

        REFERENCES:

        For more information, see Chapter 7 of [Zie2007]_.
        """
        from sage.geometry.fan import FaceFan

        if not QQ.has_coerce_map_from(self.base_ring()):
            raise NotImplementedError('face fan handles only polytopes over the rationals')

        return FaceFan(self)

    def is_minkowski_summand(self, Y):
        r"""
        Test whether ``Y`` is a Minkowski summand.

        See :meth:`~sage.geometry.polyhedron.base5.Polyhedron_base5.minkowski_sum`.

        OUTPUT: boolean; whether there exists another polyhedron `Z` such that
        ``self`` can be written as `Y\oplus Z`

        EXAMPLES::

            sage: A = polytopes.hypercube(2)
            sage: B = Polyhedron(vertices=[(0,1), (1/2,1)])
            sage: C = Polyhedron(vertices=[(1,1)])
            sage: A.is_minkowski_summand(B)
            True
            sage: A.is_minkowski_summand(C)
            True
            sage: B.is_minkowski_summand(C)
            True
            sage: B.is_minkowski_summand(A)
            False
            sage: C.is_minkowski_summand(A)
            False
            sage: C.is_minkowski_summand(B)
            False
        """
        return self.minkowski_difference(Y).minkowski_sum(Y) == self

    def barycentric_subdivision(self, subdivision_frac=None):
        r"""
        Return the barycentric subdivision of a compact polyhedron.

        DEFINITION:

        The barycentric subdivision of a compact polyhedron is a standard way
        to triangulate its faces in such a way that maximal faces correspond to
        flags of faces of the starting polyhedron (i.e. a maximal chain in the
        face lattice of the polyhedron). As a simplicial complex, this is known
        as the order complex of the face lattice of the polyhedron.

        REFERENCE:

        See :wikipedia:`Barycentric_subdivision`

        Section 6.6, Handbook of Convex Geometry, Volume A, edited by P.M. Gruber and J.M.
        Wills. 1993, North-Holland Publishing Co..

        INPUT:

        - ``subdivision_frac`` -- number. Gives the proportion how far the new
          vertices are pulled out of the polytope. Default is `\frac{1}{3}` and
          the value should be smaller than `\frac{1}{2}`. The subdivision is
          computed on the polar polyhedron.

        OUTPUT: a Polyhedron object, subdivided as described above

        EXAMPLES::

            sage: P = polytopes.hypercube(3)
            sage: P.barycentric_subdivision()
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull
            of 26 vertices
            sage: P = Polyhedron(vertices=[[0,0,0],[0,1,0],[1,0,0],[0,0,1]])
            sage: P.barycentric_subdivision()
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull
            of 14 vertices
            sage: P = Polyhedron(vertices=[[0,1,0],[0,0,1],[1,0,0]])
            sage: P.barycentric_subdivision()
            A 2-dimensional polyhedron in QQ^3 defined as the convex hull
            of 6 vertices
            sage: P = polytopes.regular_polygon(4, base_ring=QQ)                        # needs sage.rings.number_field
            sage: P.barycentric_subdivision()                                           # needs sage.rings.number_field
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 8
            vertices

        TESTS::

            sage: P.barycentric_subdivision(1/2)
            Traceback (most recent call last):
            ...
            ValueError: the subdivision fraction should be between 0 and 1/2
            sage: P = Polyhedron(ieqs=[[1,0,1],[0,1,0],[1,0,0],[0,0,1]])
            sage: P.barycentric_subdivision()
            Traceback (most recent call last):
            ...
            ValueError: the polytope has to be compact
            sage: P = Polyhedron(vertices=[[0,0,0],[0,1,0],[1,0,0],[0,0,1]], backend='field')
            sage: P.barycentric_subdivision()
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 14 vertices

            sage: polytopes.simplex(backend='field').barycentric_subdivision().backend()
            'field'
            sage: polytopes.cube(backend='cdd').barycentric_subdivision().backend()
            'cdd'
        """
        if subdivision_frac is None:
            subdivision_frac = ZZ.one() / 3

        if not self.is_compact():
            raise ValueError("the polytope has to be compact")
        if not (0 < subdivision_frac < ZZ.one() / 2):
            raise ValueError("the subdivision fraction should be "
                             "between 0 and 1/2")

        barycenter = self.center()
        parent = self.parent().base_extend(subdivision_frac)

        start_polar = (self - barycenter).polar(in_affine_span=True)
        polar = (self - barycenter).polar(in_affine_span=True)

        for i in range(self.dimension() - 1):

            new_ineq = []
            subdivided_faces = list(start_polar.faces(i))
            Hrep = polar.Hrepresentation()

            for face in subdivided_faces:

                face_vertices = face.vertices()
                normal_vectors = []

                for facet in Hrep:
                    if all(facet.contains(v) and not facet.interior_contains(v)
                           for v in face_vertices):
                        # The facet contains the face
                        normal_vectors.append(facet.A())

                normal_vector = sum(normal_vectors)
                B = - normal_vector * (face_vertices[0].vector())
                linear_evaluation = set([-normal_vector * (v.vector())
                                         for v in polar.vertices()])

                if B == max(linear_evaluation):
                    C = max(linear_evaluation.difference(set([B])))
                else:
                    C = min(linear_evaluation.difference(set([B])))

                ineq_vector = [(1 - subdivision_frac) * B + subdivision_frac * C] + list(normal_vector)
                new_ineq += [ineq_vector]

            new_ieqs = polar.inequalities_list() + new_ineq
            new_eqns = polar.equations_list()

            polar = parent.element_class(parent, None, [new_ieqs, new_eqns])

        return (polar.polar(in_affine_span=True)) + barycenter

    def permutations_to_matrices(self, conj_class_reps, acting_group=None, additional_elts=None):
        r"""
        Return a dictionary between different representations of elements in
        the ``acting_group``, with group elements represented as permutations
        of the vertices of this polytope (keys) or matrices (values).

        The dictionary has entries for the generators of the ``acting_group``
        and the representatives of conjugacy classes in ``conj_class_reps``. By
        default, the ``acting_group`` is the
        :meth:`~sage.geometry.polyhedron.base4.Polyhedron_base4.restricted_automorphism_group`
        of the polytope. Each element in ``additional_elts`` also becomes a key.

        INPUT:

        - ``conj_class_reps`` -- list; a list of representatives of the
          conjugacy classes of the ``acting_group``

        - ``acting_group`` -- a subgroup of polytope's
          :meth:`~sage.geometry.polyhedron.base4.Polyhedron_base4.restricted_automorphism_group`

        - ``additional_elts`` -- list (default: ``None``); a subset of the
          :meth:`~sage.geometry.polyhedron.base4.Polyhedron_base4.restricted_automorphism_group`
          of the polytope expressed as permutations.

        OUTPUT:

        A dictionary between elements of the ``acting_group`` expressed as permutations
        (keys) and matrices (values).

        EXAMPLES:

        This example shows the dictionary between permutations and matrices
        for the generators of the ``restricted_automorphism_group`` of the
        `\pm 1` 2-dimensional square. The permutations are written in terms
        of the vertices of the square::

            sage: # optional - pynormaliz, needs sage.groups
            sage: square = Polyhedron(vertices=[[1,1], [-1,1],
            ....:                               [-1,-1], [1,-1]],
            ....:                     backend='normaliz')
            sage: square.vertices()
            (A vertex at (-1, -1),
             A vertex at (-1, 1),
             A vertex at (1, -1),
             A vertex at (1, 1))
            sage: aut_square = square.restricted_automorphism_group(output='permutation')
            sage: conj_reps = aut_square.conjugacy_classes_representatives()
            sage: gens_dict = square.permutations_to_matrices(conj_reps)
            sage: rotation_180 = aut_square([(0,3),(1,2)])
            sage: rotation_180, gens_dict[rotation_180]
            (
                        [-1  0  0]
                        [ 0 -1  0]
            (0,3)(1,2), [ 0  0  1]
            )

        This example tests the functionality for additional elements::

            sage: # needs sage.groups sage.rings.real_mpfr
            sage: C = polytopes.cross_polytope(2)
            sage: G = C.restricted_automorphism_group(output='permutation')
            sage: conj_reps = G.conjugacy_classes_representatives()
            sage: add_elt = G([(0, 2, 3, 1)])
            sage: dict = C.permutations_to_matrices(conj_reps,
            ....:                                   additional_elts=[add_elt])
            sage: dict[add_elt]
             [ 0  1  0]
             [-1  0  0]
             [ 0  0  1]
        """
        if self.is_empty():
            raise NotImplementedError('empty polyhedra are not supported')
        if not self.is_compact():
            raise NotImplementedError('unbounded polyhedra are not supported')
        V = [v.homogeneous_vector() for v in self.Vrepresentation()]
        Qplus = sum(v.column() * v.row() for v in V).pseudoinverse()
        Vplus = list(matrix(V) * Qplus)
        W = 1 - sum(V[i].column() * Vplus[i].row() for i in range(len(V)))

        G = self.restricted_automorphism_group(output='permutation')
        if acting_group is not None:
            G = acting_group

        group_dict = {}

        def permutation_to_matrix(permutation, V, Vplus, W):
            A = sum(V[permutation(i)].column() * Vplus[i].row()
                    for i in range(len(V)))
            return A + W

        for perm in G.gens():
            group_dict[perm] = permutation_to_matrix(perm, V, Vplus, W)

        for perm in conj_class_reps:
            group_dict[perm] = permutation_to_matrix(perm, V, Vplus, W)

        if additional_elts is not None:
            for perm in additional_elts:
                group_dict[perm] = permutation_to_matrix(perm, V, Vplus, W)
        return group_dict

    @cached_method
    def bounding_box(self, integral=False, integral_hull=False):
        r"""
        Return the coordinates of a rectangular box containing the non-empty polytope.

        INPUT:

        - ``integral`` -- boolean (default: ``False``); whether to
          only allow integral coordinates in the bounding box

        - ``integral_hull`` -- boolean (default: ``False``); if ``True``, return a
          box containing the integral points of the polytope, or ``None, None`` if it
          is known that the polytope has no integral points

        OUTPUT:

        A pair of tuples ``(box_min, box_max)`` where ``box_min`` are
        the coordinates of a point bounding the coordinates of the
        polytope from below and ``box_max`` bounds the coordinates
        from above.

        EXAMPLES::

            sage: Polyhedron([(1/3,2/3), (2/3, 1/3)]).bounding_box()
            ((1/3, 1/3), (2/3, 2/3))
            sage: Polyhedron([(1/3,2/3), (2/3, 1/3)]).bounding_box(integral=True)
            ((0, 0), (1, 1))
            sage: Polyhedron([(1/3,2/3), (2/3, 1/3)]).bounding_box(integral_hull=True)
            (None, None)
            sage: Polyhedron([(1/3,2/3), (3/3, 4/3)]).bounding_box(integral_hull=True)
            ((1, 1), (1, 1))
            sage: polytopes.buckyball(exact=False).bounding_box()                       # needs sage.groups
            ((-0.8090169944, -0.8090169944, -0.8090169944),
             (0.8090169944, 0.8090169944, 0.8090169944))

        TESTS::

            sage: Polyhedron().bounding_box()
            Traceback (most recent call last):
            ...
            ValueError: empty polytope is not allowed
        """
        from sage.arith.misc import integer_ceil as ceil
        from sage.arith.misc import integer_floor as floor
        box_min = []
        box_max = []
        if not self.is_compact():
            raise ValueError("only polytopes (compact polyhedra) are allowed")
        if self.n_vertices() == 0:
            raise ValueError("empty polytope is not allowed")
        for i in range(self.ambient_dim()):
            coords = [v[i] for v in self.vertex_generator()]
            max_coord = max(coords)
            min_coord = min(coords)
            if integral_hull:
                a = ceil(min_coord)
                b = floor(max_coord)
                if a > b:
                    return None, None
                box_max.append(b)
                box_min.append(a)
            elif integral:
                box_max.append(ceil(max_coord))
                box_min.append(floor(min_coord))
            else:
                box_max.append(max_coord)
                box_min.append(min_coord)
        return (tuple(box_min), tuple(box_max))

    def _polymake_init_(self):
        """
        Return a polymake "Polytope" object corresponding to ``self``.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: PP = polymake(P)         # optional - jupymake
            sage: PP.N_VERTICES            # optional - jupymake
            8

        Lower-dimensional polyhedron::

            sage: P = Polyhedron(vertices=[[1, 0], [0, 1]])
            sage: PP = polymake(P)         # optional - jupymake
            sage: PP.COMBINATORIAL_DIM     # optional - jupymake
            1
            sage: PP.AFFINE_HULL           # optional - jupymake
            -1 1 1

        Empty polyhedron::

            sage: P = Polyhedron(ambient_dim=2, vertices=[])
            sage: PP = polymake(P)         # optional - jupymake
            sage: PP.COMBINATORIAL_DIM     # optional - jupymake
            -1

        Pointed unbounded polyhedron::

            sage: P = Polyhedron(vertices=[[1, 0], [0, 1]], rays=[[1, 0]])
            sage: PP = polymake(P)         # optional - jupymake
            sage: PP.VERTICES              # optional - jupymake
            1 0 1
            1 1 0
            0 1 0
            sage: PP.FACETS                # optional - jupymake
            1 0 -1
            -1 1 1
            0 0 1

        Non-pointed polyhedron::

            sage: # optional - jupymake
            sage: P = Polyhedron(vertices=[[1, 0], [0, 1]], lines=[[1, 0]])
            sage: PP = polymake(P)
            sage: PP.VERTICES
            1 0 1
            1 0 0
            sage: PP.FACETS
            1 0 -1
            0 0 1
            sage: PP.LINEALITY_SPACE
            0 1 0

        Algebraic polyhedron::

            sage: # needs sage.groups sage.rings.number_field
            sage: P = polytopes.dodecahedron(); P
            A 3-dimensional polyhedron
             in (Number Field in sqrt5 with defining polynomial x^2 - 5
                 with sqrt5 = 2.236067977499790?)^3
             defined as the convex hull of 20 vertices
            sage: print("Maybe recompile warning"); PP = polymake(P); PP        # optional - jupymake
            Maybe recompile warning...
            Polytope<QuadraticExtension<Rational>>[...]
            sage: sorted(PP.VERTICES[:], key=repr)[0]                           # optional - jupymake
            1 -1+1r5 -4+2r5 0

        Floating-point polyhedron::

            sage: # optional - jupymake, needs sage.groups
            sage: P = polytopes.dodecahedron(exact=False); P
            A 3-dimensional polyhedron in RDF^3 defined as the convex hull of 20 vertices
            sage: print("There may be a recompilation warning"); PP = polymake(P); PP
            There may be a recompilation warning...
            Polytope<Float>[...]
            sage: sorted(PP.VERTICES[:], key=repr)[0]
            1 -0.472135955 0 -1.236067978
        """
        from sage.interfaces.polymake import polymake
        polymake_field = polymake(self.base_ring().fraction_field())
        polymake_class = "Polytope<{}>".format(polymake_field)
        if self.is_empty():
            # Polymake 3.1 cannot enter an empty polyhedron using
            # FACETS and AFFINE_HULL.
            # Use corresponding input properties instead.
            # https://forum.polymake.org/viewtopic.php?f=8&t=545
            return polymake.new_object(polymake_class,
                                       INEQUALITIES=self.inequalities_list(),
                                       EQUATIONS=self.equations_list())

        verts_and_rays = [[1] + v for v in self.vertices_list()]
        verts_and_rays += [[0] + r for r in self.rays_list()]
        return polymake.new_object(polymake_class,
                                   FACETS=self.inequalities_list(),
                                   AFFINE_HULL=self.equations_list(),
                                   VERTICES=verts_and_rays,
                                   LINEALITY_SPACE=[[0] + l for l in self.lines_list()])
