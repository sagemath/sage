r"""
Manifold Subsets Defined as Pullbacks of Subsets under Continuous Maps
"""


# ****************************************************************************
#       Copyright (C) 2021 Matthias Koeppe <mkoeppe@math.ucdavis.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.sets_cat import Sets, EmptySetError
from sage.categories.metric_spaces import MetricSpaces
from sage.misc.lazy_import import lazy_import
from sage.modules.free_module import FreeModule_generic
from sage.rings.infinity import infinity, minus_infinity
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.complex_double import CDF
from sage.rings.real_double import RDF
from sage.rings.real_lazy import CLF, RLF
from sage.symbolic.ring import SR
from sage.modules.free_module_element import vector
from sage.manifolds.subset import ManifoldSubset
from sage.manifolds.chart import Chart
from sage.manifolds.scalarfield import ScalarField
from sage.sets.real_set import RealSet
import sage.geometry.abc

lazy_import('sage.geometry.relative_interior', 'RelativeInterior')


class ManifoldSubsetPullback(ManifoldSubset):
    """
    Manifold subset defined as a pullback of a subset under a continuous map.

    INPUT:

    - ``map`` -- an instance of :class:`~sage.manifolds.continuous_map.ContinuousMap`,
      :class:`ScalarField`, or :class:`Chart`

    - ``codomain_subset`` -- an instance of :class:`~sage.manifolds.subset.ManifoldSubset`,
      :class:`RealSet`, or :class:`~sage.geometry.convex_set.ConvexSet_base`

    EXAMPLES::

        sage: from sage.manifolds.subsets.pullback import ManifoldSubsetPullback
        sage: M = Manifold(2, 'R^2', structure='topological')
        sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2

    Pulling back a real interval under a scalar field::

        sage: r_squared = M.scalar_field(x^2+y^2)
        sage: r_squared.set_immutable()
        sage: cl_I = RealSet([1, 4]); cl_I
        [1, 4]
        sage: cl_O = ManifoldSubsetPullback(r_squared, cl_I); cl_O
        Subset f_inv_[1, 4] of the 2-dimensional topological manifold R^2
        sage: M.point((0, 0)) in cl_O
        False
        sage: M.point((0, 1)) in cl_O
        True

    Pulling back an open real interval gives an open subset::

        sage: I = RealSet((1, 4)); I
        (1, 4)
        sage: O = ManifoldSubsetPullback(r_squared, I); O
        Open subset f_inv_(1, 4) of the 2-dimensional topological manifold R^2
        sage: M.point((1, 0)) in O
        False
        sage: M.point((1, 1)) in O
        True

    Pulling back a polytope under a chart::

        sage: # needs sage.geometry.polyhedron
        sage: P = Polyhedron(vertices=[[0, 0], [1, 2], [2, 1]]); P
        A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 3 vertices
        sage: S = ManifoldSubsetPullback(c_cart, P); S
        Subset x_y_inv_P of the 2-dimensional topological manifold R^2
        sage: M((1, 2)) in S
        True
        sage: M((2, 0)) in S
        False

    Pulling back the interior of a polytope under a chart::

        sage: # needs sage.geometry.polyhedron
        sage: int_P = P.interior(); int_P
        Relative interior of a
         2-dimensional polyhedron in ZZ^2 defined as the convex hull of 3 vertices
        sage: int_S = ManifoldSubsetPullback(c_cart, int_P, name='int_S'); int_S
        Open subset int_S of the 2-dimensional topological manifold R^2
        sage: M((0, 0)) in int_S
        False
        sage: M((1, 1)) in int_S
        True

    Using the embedding map of a submanifold::

        sage: M = Manifold(3, 'M', structure='topological')
        sage: N = Manifold(2, 'N', ambient=M, structure='topological'); N
        2-dimensional topological submanifold N
         immersed in the 3-dimensional topological manifold M
        sage: CM.<x,y,z> = M.chart()
        sage: CN.<u,v> = N.chart()
        sage: t = var('t')
        sage: phi = N.continuous_map(M, {(CN,CM): [u,v,t+u^2+v^2]})
        sage: phi_inv = M.continuous_map(N, {(CM,CN): [x,y]})
        sage: phi_inv_t = M.scalar_field({CM: z-x^2-y^2})
        sage: N.set_immersion(phi, inverse=phi_inv, var=t,
        ....:                 t_inverse={t: phi_inv_t})
        sage: N.declare_embedding()

        sage: from sage.manifolds.subsets.pullback import ManifoldSubsetPullback
        sage: S = M.open_subset('S', coord_def={CM: z<1})
        sage: phi_without_t = N.continuous_map(M, {(CN, CM): [expr.subs(t=0)
        ....:                                                 for expr in phi.expr()]})
        sage: phi_without_t
        Continuous map
         from the 2-dimensional topological submanifold N
          embedded in the 3-dimensional topological manifold M
         to the 3-dimensional topological manifold M
        sage: phi_without_t.expr()
        (u, v, u^2 + v^2)
        sage: D = ManifoldSubsetPullback(phi_without_t, S); D
        Subset f_inv_S of the 2-dimensional topological submanifold N
         embedded in the 3-dimensional topological manifold M
        sage: N.point((2,0)) in D
        False
    """
    @staticmethod
    def __classcall_private__(cls, map, codomain_subset, inverse=None,
                              name=None, latex_name=None):
        """
        Normalize arguments and delegate to other constructors.

        TESTS::

            sage: # needs sage.geometry.polyhedron
            sage: from sage.manifolds.subsets.pullback import ManifoldSubsetPullback
            sage: M = Manifold(2, 'R^2', structure='topological')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
            sage: P = Polyhedron(vertices=[[0, 0], [1, 2], [3, 4]]); P
            A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 3 vertices
            sage: S = ManifoldSubsetPullback(c_cart, P); S
            Subset x_y_inv_P of the 2-dimensional topological manifold R^2
            sage: S is ManifoldSubsetPullback(c_cart, P)
            True
        """

        try:
            is_mutable = map.is_mutable()
        except AttributeError:
            pass
        else:
            if is_mutable:
                map = map.copy()
                map.set_immutable()
        try:
            is_mutable = inverse.is_mutable()
        except AttributeError:
            pass
        else:
            if is_mutable:
                inverse = inverse.copy()
                inverse.set_immutable()

        if inverse is None:
            if isinstance(map, Chart):
                from sage.misc.latex import latex
                inverse_latex_name = '(' + ','.join(str(latex(x)) + '^{-1}' for x in map) + ')'
                inverse_name = '_'.join(repr(x) for x in map) + '_inv'
            else:
                map_name = map._name or 'f'
                map_latex_name = map._latex_name or map_name
                inverse_name = map_name + '_inv'
                inverse_latex_name = map_latex_name + r'^{-1}'
        else:
            inverse_name = inverse._name
            inverse_latex_name = inverse._latex_name
        try:
            codomain_subset_latex_name = codomain_subset._latex_name
            codomain_subset_name = codomain_subset._name
        except AttributeError:
            from sage.misc.latex import latex
            codomain_subset_latex_name = str(latex(codomain_subset))
            s = repr(codomain_subset)
            if len(s) > 10:
                codomain_subset_name = 'P'
            else:
                codomain_subset_name = s
        if latex_name is None:
            if name is None:
                latex_name = inverse_latex_name + '(' + codomain_subset_latex_name + ')'
            else:
                latex_name = name
        if name is None:
            name = inverse_name + '_' + codomain_subset_name

        if cls._is_open(codomain_subset):

            try:
                coord_def = cls._coord_def(map, codomain_subset)
            except NotImplementedError:
                pass
            else:
                return map.domain().open_subset(name=name, latex_name=latex_name,
                                                coord_def=coord_def)

        self = super().__classcall__(cls, map, codomain_subset, inverse, name, latex_name)

        return self

    @staticmethod
    def _is_open(codomain_subset):
        """
        Return whether ``codomain_subset`` is (known to be) an open subset of its ambient space.

        EXAMPLES:

        Manifolds and subsets::

            sage: from sage.manifolds.subsets.pullback import ManifoldSubsetPullback
            sage: R2 = Manifold(2, 'R^2', structure='topological'); R2
            2-dimensional topological manifold R^2
            sage: ManifoldSubsetPullback._is_open(R2)
            True
            sage: A = R2.subset('A'); A
            Subset A of the 2-dimensional topological manifold R^2
            sage: ManifoldSubsetPullback._is_open(A)
            False

        :class:`RealSet` instances::

            sage: I = RealSet.open(1, 2); I
            (1, 2)
            sage: ManifoldSubsetPullback._is_open(I)
            True
            sage: cl_I = RealSet.closed(1, 2); cl_I
            [1, 2]
            sage: ManifoldSubsetPullback._is_open(cl_I)
            False

        Polyhedra::

            sage: # needs sage.geometry.polyhedron
            sage: Empty = Polyhedron(ambient_dim=2); Empty
            The empty polyhedron in ZZ^2
            sage: ManifoldSubsetPullback._is_open(Empty)
            True
            sage: C = polytopes.cube(); C
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 8 vertices
            sage: ManifoldSubsetPullback._is_open(C)
            False

        Interiors of polyhedra::

            sage: int_C = C.interior(); int_C                                           # needs sage.geometry.polyhedron
            Relative interior of a
             3-dimensional polyhedron in ZZ^3 defined as the convex hull of 8 vertices
            sage: ManifoldSubsetPullback._is_open(int_C)                                # needs sage.geometry.polyhedron
            True

        PPL polyhedra and not-necessarily-closed polyhedra::

            sage: # needs pplpy
            sage: from ppl import Variable, C_Polyhedron, NNC_Polyhedron, Constraint_System
            sage: u = Variable(0)
            sage: v = Variable(1)
            sage: CS = Constraint_System()
            sage: CS.insert(0 < u)
            sage: CS.insert(u < 1)
            sage: CS.insert(0 < v)
            sage: CS.insert(v < 1)
            sage: CS.insert(u + v <= 3)       # redundant inequality
            sage: P = NNC_Polyhedron(CS); P
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 point, 4 closure_points
            sage: ManifoldSubsetPullback._is_open(P)
            True
            sage: CS.insert(u + v <= 1)
            sage: T = NNC_Polyhedron(CS); T
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 point, 3 closure_points
            sage: ManifoldSubsetPullback._is_open(T)
            False
        """

        if isinstance(codomain_subset, ManifoldSubset):
            return codomain_subset.is_open()

        if isinstance(codomain_subset, RealSet):
            return codomain_subset.is_open()

        if isinstance(codomain_subset, sage.geometry.abc.Polyhedron):
            return codomain_subset.is_empty() or codomain_subset.is_universe()

        if isinstance(codomain_subset, RelativeInterior):
            return codomain_subset.closure().is_full_dimensional()

        if codomain_subset in Sets().Finite():
            return codomain_subset.cardinality() == 0

        if hasattr(codomain_subset, 'minimized_constraints'):
            try:
                from ppl import NNC_Polyhedron, C_Polyhedron
            except ImportError:
                pass
            else:
                if isinstance(codomain_subset, (NNC_Polyhedron, C_Polyhedron)):
                    cs = codomain_subset.minimized_constraints()
                    if cs.has_equalities():
                        return False
                    if any(constraint.is_nonstrict_inequality()
                           for constraint in cs):
                        return False
                    return True

        return False

    @staticmethod
    def _interval_restriction(expr, interval):
        """
        Return a restriction expressing that ``expr`` lies in ``interval``.

        INPUT:

        - ``expr`` -- a symbolic expression
        - ``interval`` -- an instance of :class:`~sage.sets.real_set.InternalRealInterval`

        OUTPUT:

        - A restriction suitable as input to :meth:`~sage.manifolds.chart.restrict`:
          lists are conjunctions, tuples are disjunctions

        EXAMPLES::

            sage: from sage.manifolds.subsets.pullback import ManifoldSubsetPullback
            sage: _interval_restriction = ManifoldSubsetPullback._interval_restriction
            sage: var('t')
            t
            sage: assume(t >= -2)
            sage: assume(t <= 5)
            sage: _interval_restriction(t, RealSet(3, 4)[0])
            [t > 3, t < 4]
            sage: _interval_restriction(t, RealSet.unbounded_below_closed(2)[0])
            t <= 2
            sage: _interval_restriction(t, RealSet.closed(-5, 5)[0])
            []
            sage: _interval_restriction(t, RealSet.unbounded_below_closed(-5)[0])
            ()
            sage: _interval_restriction(t, RealSet.unbounded_above_closed(6)[0])
            ()
            sage: _interval_restriction(t^2, RealSet.unbounded_above_closed(0)[0])
            []
        """

        conjunction = []
        if interval.lower() != minus_infinity:
            if interval.lower_closed():
                condition = (expr >= interval.lower())
                negation = (expr < interval.lower())
            else:
                condition = (expr > interval.lower())
                negation = (expr <= interval.lower())
            if negation:
                # known to be false
                return ()
            if not condition:
                # not known to be true
                conjunction.append(condition)

        if interval.upper() != infinity:
            if interval.upper_closed():
                condition = (expr <= interval.upper())
                negation = (expr > interval.upper())
            else:
                condition = (expr < interval.upper())
                negation = (expr >= interval.upper())
            if negation:
                # known to be false
                return ()
            if not condition:
                # not known to be true
                conjunction.append(condition)

        if len(conjunction) == 1:
            return conjunction[0]
        else:
            # lists express 'and'
            return conjunction

    @staticmethod
    def _realset_restriction(expr, realset):
        """
        Return a restriction expressing that ``expr`` lies in ``realset``.

        INPUT:

        - ``expr`` -- a symbolic expression
        - ``interval`` -- an instance of :class:`~sage.sets.real_set.RealSet`

        OUTPUT:

        - A restriction suitable as input to :meth:`~sage.manifolds.chart.restrict`:
          lists are conjunctions, tuples are disjunctions

        EXAMPLES::

            sage: from sage.manifolds.subsets.pullback import ManifoldSubsetPullback
            sage: _realset_restriction = ManifoldSubsetPullback._realset_restriction
            sage: var('t')
            t
            sage: assume(t >= -2)
            sage: assume(t <= 5)
            sage: _realset_restriction(t, RealSet(-oo, oo))
            []
            sage: _realset_restriction(t, RealSet())
            ()
            sage: _realset_restriction(t, RealSet([-5, -4], (-1, 1), [3, 4], [6, 7]))
            ([t > -1, t < 1], [t >= 3, t <= 4])
        """
        disjunction = []
        for interval in realset:
            condition = ManifoldSubsetPullback._interval_restriction(expr, interval)
            if condition == []:
                return []
            if condition != ():
                disjunction.append(condition)

        if len(disjunction) == 1:
            return disjunction[0]
        else:
            # tuples express 'or'
            return tuple(disjunction)

    @staticmethod
    def _polyhedron_restriction(expr, polyhedron, relint=False):
        """
        Return a restriction expressing that ``expr`` lies in ``polyhedron`` or its relative interior.

        INPUT:

        - ``expr`` -- a symbolic expression
        - ``polyhedron`` -- an instance of :class:`~sage.geometry.polyhedron.base.Polyhedron_base`
        - ``relint`` -- whether the restriction should use the relative interior

        OUTPUT:

        - A restriction suitable as input to :meth:`~sage.manifolds.chart.restrict`:
          lists are conjunctions, tuples are disjunctions

        EXAMPLES::

            sage: from sage.manifolds.subsets.pullback import ManifoldSubsetPullback
            sage: _polyhedron_restriction = ManifoldSubsetPullback._polyhedron_restriction
            sage: var('x y z')
            (x, y, z)
            sage: c = polytopes.cube()                                                  # needs sage.geometry.polyhedron
            sage: _polyhedron_restriction((x, y, z), c)                                 # needs sage.geometry.polyhedron
            [-x + 1 >= 0, -y + 1 >= 0, -z + 1 >= 0, x + 1 >= 0, z + 1 >= 0, y + 1 >= 0]
            sage: _polyhedron_restriction((x, y, z), c, relint=True)                    # needs sage.geometry.polyhedron
            [-x + 1 > 0, -y + 1 > 0, -z + 1 > 0, x + 1 > 0, z + 1 > 0, y + 1 > 0]
        """
        conjunction = []

        expr = vector(SR, expr)
        for constraint in polyhedron.Hrepresentation():

            if constraint.is_inequality():
                if relint:
                    condition = (constraint.eval(expr) > 0)
                else:
                    condition = (constraint.eval(expr) >= 0)
            else:
                condition = (constraint.eval(expr) == 0)
            if not condition:
                # not known to be true
                conjunction.append(condition)

        if len(conjunction) == 1:
            return conjunction[0]
        else:
            # lists express 'and'
            return conjunction

    @staticmethod
    def _coord_def(map, codomain_subset):
        r"""
        Return a coordinate definition of the open subset that is the pullback of ``codomain_subset``.

        INPUT:

        - ``map`` -- an instance of :class:`ScalarField` or :class:`Chart`

        - ``codomain_subset`` -- if ``map`` is a :class:`ScalarField`, an
          instance of :class:`RealSet`; if ``map`` is a :class:`Chart`, the
          relative interior of a polyhedron

        For other inputs, a :exc:`NotImplementedError` will be raised.

        OUTPUT:

        - an object suitable for the parameter ``coord_def`` of
          :meth:`sage.manifolds.manifold.TopologicalManifold.open_subset`.

        EXAMPLES::

            sage: from sage.manifolds.subsets.pullback import ManifoldSubsetPullback
            sage: _coord_def = ManifoldSubsetPullback._coord_def
            sage: M = Manifold(2, 'R^2', structure='topological')

        Coordinate definition of an open chart polyhedron::

            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
            sage: P = Polyhedron(vertices=[[0, 0], [1, 2], [3, 4]]); P                  # needs sage.geometry.polyhedron
            A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 3 vertices
            sage: ri_P = P.relative_interior(); ri_P                                    # needs sage.geometry.polyhedron
            Relative interior of a 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 3 vertices
            sage: _coord_def(c_cart, ri_P)                                              # needs sage.geometry.polyhedron
            {Chart (R^2, (x, y)): [2*x - y > 0, -4*x + 3*y > 0, x - y + 1 > 0]}

        Coordinate definition of the pullback of an open interval under a scalar field::

            sage: r_squared = M.scalar_field(x^2 + y^2)
            sage: I = RealSet((1, 4)); I
            (1, 4)
            sage: _coord_def(r_squared, I)
            {Chart (R^2, (x, y)): [x^2 + y^2 > 1, x^2 + y^2 < 4]}
        """
        if isinstance(map, ScalarField) and isinstance(codomain_subset, RealSet):

            return {chart: ManifoldSubsetPullback._realset_restriction(func.expr(),
                                                                       codomain_subset)
                    for chart, func in map._express.items()}

        if isinstance(map, Chart):

            chart = map

            if isinstance(codomain_subset, RealSet):
                return {chart: ManifoldSubsetPullback._realset_restriction(chart[0],
                                                                           codomain_subset)}

            if isinstance(codomain_subset, RelativeInterior) and isinstance(codomain_subset.closure(), sage.geometry.abc.Polyhedron):
                return {chart: ManifoldSubsetPullback._polyhedron_restriction(
                                   chart, codomain_subset.closure(), relint=True)}

        raise NotImplementedError

    def __init__(self, map, codomain_subset, inverse, name, latex_name):
        r"""
        Construct a manifold subset that is a pullback.

        TESTS::

            sage: from sage.manifolds.subsets.pullback import ManifoldSubsetPullback
            sage: M = Manifold(2, 'R^2', structure='topological')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
            sage: r_squared = M.scalar_field(x^2+y^2)
            sage: r_squared.set_immutable()
            sage: cl_I = RealSet([1, 4]); cl_I
            [1, 4]
            sage: cl_O = ManifoldSubsetPullback(r_squared, cl_I); cl_O
            Subset f_inv_[1, 4] of the 2-dimensional topological manifold R^2
            sage: TestSuite(cl_O).run(skip='_test_elements')
        """
        if inverse is None and isinstance(map, Chart):
            chart = map
            scalar_codomain = (isinstance(codomain_subset, RealSet)
                               or any(field.has_coerce_map_from(codomain_subset)
                                      for field in (CDF, RDF, CLF, RLF)))
            if scalar_codomain:
                if chart.domain().dimension() != 1:
                    raise ValueError('to pull back a set of scalars by a chart, the manifold must be 1-dimensional')
                map = chart.domain().scalar_field({chart: chart[0]})

                def _inverse(coord):
                    return self.point((coord,), chart=chart)
            else:
                def _inverse(coords):
                    return self.point(coords, chart=map)
                inverse = _inverse

        self._map = map
        self._inverse = inverse

        self._codomain_subset = codomain_subset
        base_manifold = map.domain()
        ManifoldSubset.__init__(self, base_manifold, name, latex_name=latex_name)

    def _an_element_(self):
        r"""
        Construct some point in ``self``.

        EXAMPLES::

            sage: # needs sage.geometry.polyhedron
            sage: from sage.manifolds.subsets.pullback import ManifoldSubsetPullback
            sage: M = Manifold(3, 'R^3', structure='topological')
            sage: c_cart.<x,y,z> = M.chart() # Cartesian coordinates on R^3
            sage: Cube = polytopes.cube(); Cube
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 8 vertices
            sage: McCube = ManifoldSubsetPullback(c_cart, Cube, name='McCube'); McCube
            Subset McCube of the 3-dimensional topological manifold R^3
            sage: p = McCube.an_element(); p
            Point on the 3-dimensional topological manifold R^3
            sage: p.coordinates(c_cart)
            (0, 0, 0)

            sage: # needs sage.geometry.polyhedron
            sage: Empty = Polyhedron(ambient_dim=3)
            sage: McEmpty = ManifoldSubsetPullback(c_cart, Empty, name='McEmpty')
            sage: McEmpty
            Subset McEmpty of the 3-dimensional topological manifold R^3
            sage: McEmpty.an_element()
            Traceback (most recent call last):
            ...
            sage.categories.sets_cat.EmptySetError
        """
        try:
            return next(iter(self.some_elements()))
        except StopIteration:
            raise EmptySetError

    def some_elements(self):
        r"""
        Generate some elements of ``self``.

        EXAMPLES::

            sage: # needs sage.geometry.polyhedron
            sage: from sage.manifolds.subsets.pullback import ManifoldSubsetPullback
            sage: M = Manifold(3, 'R^3', structure='topological')
            sage: c_cart.<x,y,z> = M.chart() # Cartesian coordinates on R^3
            sage: Cube = polytopes.cube(); Cube
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 8 vertices
            sage: McCube = ManifoldSubsetPullback(c_cart, Cube, name='McCube'); McCube
            Subset McCube of the 3-dimensional topological manifold R^3
            sage: L = list(McCube.some_elements()); L
            [Point on the 3-dimensional topological manifold R^3,
             Point on the 3-dimensional topological manifold R^3,
             Point on the 3-dimensional topological manifold R^3,
             Point on the 3-dimensional topological manifold R^3,
             Point on the 3-dimensional topological manifold R^3,
             Point on the 3-dimensional topological manifold R^3]
            sage: list(p.coordinates(c_cart) for p in L)
            [(0, 0, 0),
             (1, -1, -1),
             (1, 0, -1),
             (1, 1/2, 0),
             (1, -1/4, 1/2),
             (0, -5/8, 3/4)]

            sage: # needs sage.geometry.polyhedron
            sage: Empty = Polyhedron(ambient_dim=3)
            sage: McEmpty = ManifoldSubsetPullback(c_cart, Empty, name='McEmpty')
            sage: McEmpty
            Subset McEmpty of the 3-dimensional topological manifold R^3
            sage: list(McEmpty.some_elements())
            []
        """
        if self._inverse is not None:
            for y in self._codomain_subset.some_elements():
                yield self._inverse(y)
        elif self.is_empty():
            return
        else:
            # Fallback
            p = super()._an_element_()
            if p in self:
                yield p

    def __contains__(self, point):
        r"""
        Check whether ``point`` is contained in ``self``.

        EXAMPLES::

            sage: # needs sage.geometry.polyhedron
            sage: from sage.manifolds.subsets.pullback import ManifoldSubsetPullback
            sage: M = Manifold(3, 'R^3', structure='topological')
            sage: c_cart.<x,y,z> = M.chart() # Cartesian coordinates on R^3
            sage: Cube = polytopes.cube(); Cube
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 8 vertices
            sage: Cube.vertices_list()
            [[1, -1, -1],
            [1, 1, -1],
            [1, 1, 1],
            [1, -1, 1],
            [-1, -1, 1],
            [-1, -1, -1],
            [-1, 1, -1],
            [-1, 1, 1]]
            sage: McCube = ManifoldSubsetPullback(c_cart, Cube, name='McCube'); McCube
            Subset McCube of the 3-dimensional topological manifold R^3
            sage: p = M.point((0, 0, 0)); p
            Point on the 3-dimensional topological manifold R^3
            sage: p in McCube
            True
            sage: q = M.point((2, 3, 4)); q
            Point on the 3-dimensional topological manifold R^3
            sage: q in McCube
            False
        """
        if super().__contains__(point):
            return True
        coords = self._map(point)
        if isinstance(coords, (tuple, list)):
            coords = vector(coords)
        return coords in self._codomain_subset

    def is_open(self):
        """
        Return if ``self`` is (known to be) an open set.

        This version of the method always returns ``False``.

        Because the map is continuous, the pullback is open if the
        ``codomain_subset`` is open.

        However, the design of :class:`~sage.manifolds.subset.ManifoldSubset` requires that open subsets
        are instances of the subclass :class:`sage.manifolds.manifold.TopologicalManifold`.
        The constructor of :class:`ManifoldSubsetPullback` delegates to a subclass
        of :class:`sage.manifolds.manifold.TopologicalManifold` for some open subsets.

        EXAMPLES::

            sage: from sage.manifolds.subsets.pullback import ManifoldSubsetPullback
            sage: M = Manifold(2, 'R^2', structure='topological')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2

            sage: # needs sage.geometry.polyhedron
            sage: P = Polyhedron(vertices=[[0, 0], [1, 2], [3, 4]]); P
            A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 3 vertices
            sage: P.is_open()
            False
            sage: McP = ManifoldSubsetPullback(c_cart, P, name='McP'); McP
            Subset McP of the 2-dimensional topological manifold R^2
            sage: McP.is_open()
            False
        """
        return super().is_open()

    def is_closed(self):
        """
        Return if ``self`` is (known to be) a closed subset of the manifold.

        EXAMPLES::

            sage: from sage.manifolds.subsets.pullback import ManifoldSubsetPullback
            sage: M = Manifold(2, 'R^2', structure='topological')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2

        The pullback of a closed real interval under a scalar field is closed::

            sage: r_squared = M.scalar_field(x^2+y^2)
            sage: r_squared.set_immutable()
            sage: cl_I = RealSet([1, 2]); cl_I
            [1, 2]
            sage: cl_O = ManifoldSubsetPullback(r_squared, cl_I); cl_O
            Subset f_inv_[1, 2] of the 2-dimensional topological manifold R^2
            sage: cl_O.is_closed()
            True

        The pullback of a (closed convex) polyhedron under a chart is closed::

            sage: # needs sage.geometry.polyhedron
            sage: P = Polyhedron(vertices=[[0, 0], [1, 2], [3, 4]]); P
            A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 3 vertices
            sage: McP = ManifoldSubsetPullback(c_cart, P, name='McP'); McP
            Subset McP of the 2-dimensional topological manifold R^2
            sage: McP.is_closed()
            True

        The pullback of real vector subspaces under a chart is closed::

            sage: V = span([[1, 2]], RR); V
            Vector space of degree 2 and dimension 1 over Real Field with 53 bits of precision
            Basis matrix:
            [1.00000000000000 2.00000000000000]
            sage: McV = ManifoldSubsetPullback(c_cart, V, name='McV'); McV
            Subset McV of the 2-dimensional topological manifold R^2
            sage: McV.is_closed()
            True

        The pullback of point lattices under a chart is closed::

            sage: W = span([[1, 0], [3, 5]], ZZ); W
            Free module of degree 2 and rank 2 over Integer Ring
            Echelon basis matrix:
            [1 0]
            [0 5]
            sage: McW = ManifoldSubsetPullback(c_cart, W, name='McW'); McW
            Subset McW of the 2-dimensional topological manifold R^2
            sage: McW.is_closed()
            True

        The pullback of finite sets is closed::

            sage: F = Family([vector(QQ, [1, 2], immutable=True), vector(QQ, [2, 3], immutable=True)])
            sage: McF = ManifoldSubsetPullback(c_cart, F, name='McF'); McF
            Subset McF of the 2-dimensional topological manifold R^2
            sage: McF.is_closed()
            True
        """
        if self.manifold().dimension() == 0:
            return True
        if isinstance(self._codomain_subset, ManifoldSubset):
            if self._codomain_subset.is_closed():
                # known closed
                return True
        elif isinstance(self._codomain_subset, RealSet):
            # RealSet can decide closedness authoritatively
            return self._codomain_subset.is_closed()
        elif isinstance(self._codomain_subset, sage.geometry.abc.Polyhedron):
            # Regardless of their base_ring, we treat polyhedra as closed
            # convex subsets of R^n
            return True
        elif isinstance(self._codomain_subset, FreeModule_generic) and self._codomain_subset.rank() != infinity:
            if self._codomain_subset.base_ring() in MetricSpaces().Complete():
                # Closed topological vector subspace
                return True
            if self._codomain_subset.base_ring() == ZZ:
                if self._codomain_subset.coordinate_ring().is_subring(QQ):
                    # Discrete subgroup of R^n
                    return True
                if self._codomain_subset.rank() == self._codomain_subset.base_extend(RR).dimension():
                    # Discrete subgroup of R^n
                    return True
        elif self._codomain_subset in Sets().Finite():
            return True
        else:
            if hasattr(self._codomain_subset, 'is_topologically_closed'):
                try:
                    from ppl import NNC_Polyhedron, C_Polyhedron
                except ImportError:
                    pass
                else:
                    if isinstance(self._codomain_subset, (NNC_Polyhedron, C_Polyhedron)):
                        # ppl polyhedra can decide closedness authoritatively
                        return self._codomain_subset.is_topologically_closed()
        return super().is_closed()

    def closure(self, name=None, latex_name=None):
        """
        Return the topological closure of ``self`` in the manifold.

        Because ``self`` is a pullback of some subset under a continuous map,
        the closure of ``self`` is the pullback of the closure.

        EXAMPLES::

            sage: from sage.manifolds.subsets.pullback import ManifoldSubsetPullback
            sage: M = Manifold(2, 'R^2', structure='topological')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
            sage: r_squared = M.scalar_field(x^2+y^2)
            sage: r_squared.set_immutable()
            sage: I = RealSet.open_closed(1, 2); I
            (1, 2]
            sage: O = ManifoldSubsetPullback(r_squared, I); O
            Subset f_inv_(1, 2] of the 2-dimensional topological manifold R^2
            sage: latex(O)
            f^{-1}((1, 2])
            sage: cl_O = O.closure(); cl_O
            Subset f_inv_[1, 2] of the 2-dimensional topological manifold R^2
            sage: cl_O.is_closed()
            True
        """
        if self.is_closed():
            return self
        try:
            codomain_subset_closure = self._codomain_subset.closure()
        except AttributeError:
            return super().closure()
        closure = ManifoldSubsetPullback(self._map, codomain_subset_closure,
                                         inverse=self._inverse,
                                         name=name, latex_name=latex_name)
        closure.declare_superset(self)
        return closure
