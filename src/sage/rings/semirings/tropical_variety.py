r"""
Tropical Varieties

A tropical variety is a piecewise-linear geometric object derived from
a classical algebraic variety by using tropical mathematics, where the
tropical semiring replaces the usual arithmetic operations.

AUTHORS:

- Verrel Rievaldo Wijaya (2024-06): initial version

REFERENCES:

- [Bru2014]_
- [Mac2015]_
- [Fil2017]_
"""

# ****************************************************************************
#       Copyright (C) 2024 Verrel Rievaldo Wijaya <verrelrievaldo@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation

from sage.rings.rational_field import QQ
from sage.rings.infinity import infinity

class TropicalVariety(UniqueRepresentation, SageObject):
    r"""
    A tropical variety in `\RR^n`.

    A tropical variety is defined as a corner locus of tropical polynomial
    function. This means it consist of all points in `\RR^n` for which
    the minimum (maximum) of the function is attained at least twice.

    We represent the tropical variety as a list of lists, where the
    inner list consist of three parts. The first one is a parametric
    equations for tropical roots. The second one is the condition
    for parameters. The third one is the order of the corresponding
    component.

    INPUT:

    - ``poly`` -- a :class:`TropicalMPolynomial`

    ALGORITHM:

    We need to determine a corner locus of this tropical polynomial
    function, which is all points `(x_1, x_2, \ldots, x_n)` for which
    the maximum (minimum) is obtained at least twice. First, we convert
    each monomial to its corresponding linear function. Then for each two
    monomials of polynomial, we find the points where their values are
    equal. Since we attempt to solve the equality of two equations in `n`
    variables, the solution set will be described by `n-1` parameters.

    Next, we need to check if the value of previous two monomials at the
    points in solution set is really the maximum (minimum) of function.
    We do this by solving the inequality of the previous monomial with all
    other monomials in the polynomial after substituting the parameter.
    This will give us the condition of parameters. Each of this condition
    is then combined by union operator. If this final condition is not an
    empty set, then it represent one component of tropical root. Then we
    calculate the weight of this particular component by the maximum of
    gcd of the numbers `|i-k|` and `|j-l|` for all pairs `(i,j)` and
    `(k,l)` such that the value of on this component is given by the
    corresponding monomials.

    EXAMPLES:

    We construct a tropical variety in `\RR^2`, where it is called a
    tropical curve::

        sage: T = TropicalSemiring(QQ, use_min=False)
        sage: R.<x,y> = PolynomialRing(T)
        sage: p1 = R(1)*x + x*y + R(0); p1
        0*x*y + 1*x + 0
        sage: tv = p1.tropical_variety(); tv
        Tropical curve of 0*x*y + 1*x + 0
        sage: tv.components()
        [[(t1, 1), [t1 >= -1], 1], [(-1, t1), [t1 <= 1], 1], [(-t1, t1), [t1 >= 1], 1]]
        sage: tv.vertices()
        {(-1, 1)}
        sage: tv.plot()
        Graphics object consisting of 3 graphics primitives

    .. PLOT::
        :width: 300 px

        T = TropicalSemiring(QQ, use_min=False)
        R = PolynomialRing(T, ('x,y'))
        x, y = R.gen(), R.gen(1)
        p1 = R(1)*x + x*y + R(0)
        sphinx_plot(p1.tropical_variety().plot())

    A slightly different result will be obtained if we use min-plus algebra
    for the base tropical semiring::

        sage: T = TropicalSemiring(QQ, use_min=True)
        sage: R.<x,y> = PolynomialRing(T)
        sage: p1 = R(1)*x + x*y + R(0)
        sage: tv = p1.tropical_variety(); tv
        Tropical curve of 0*x*y + 1*x + 0
        sage: tv.components()
        [[(t1, 1), [t1 <= -1], 1], [(-1, t1), [t1 >= 1], 1], [(-t1, t1), [t1 <= 1], 1]]
        sage: tv.plot()
        Graphics object consisting of 3 graphics primitives

    .. PLOT::
        :width: 300 px

        T = TropicalSemiring(QQ, use_min=True)
        R = PolynomialRing(T, ('x,y'))
        x, y = R.gen(), R.gen(1)
        p1 = R(1)*x + x*y + R(0)
        sphinx_plot(p1.tropical_variety().plot())

    Tropical variety can consist of multiple components with varying orders::

        sage: T = TropicalSemiring(QQ, use_min=False)
        sage: R.<x,y> = PolynomialRing(T)
        sage: p1 = R(7) + T(4)*x + y + R(4)*x*y + R(3)*y^2 + R(-3)*x^2
        sage: tv = p1.tropical_variety(); tv
        Tropical curve of (-3)*x^2 + 4*x*y + 3*y^2 + 4*x + 0*y + 7
        sage: tv.components()
        [[(3, t1), [t1 <= 0], 1],
        [(-t1 + 3, t1), [0 <= t1, t1 <= 2], 1],
        [(t1, 2), [t1 <= 1], 2],
        [(t1, 0), [3 <= t1, t1 <= 7], 1],
        [(7, t1), [t1 <= 0], 1],
        [(t1 - 1, t1), [2 <= t1], 1],
        [(t1 + 7, t1), [0 <= t1], 1]]
        sage: tv.plot()
        Graphics object consisting of 8 graphics primitives

    .. PLOT::
        :width: 300 px

        T = TropicalSemiring(QQ, use_min=False)
        R = PolynomialRing(T, ('x,y'))
        x, y = R.gen(), R.gen(1)
        p1 = R(7) + T(4)*x + y + R(4)*x*y + R(3)*y**2 + R(-3)*x**2
        sphinx_plot(p1.tropical_variety().plot())

    If the tropical polynomial have `n>2` variables, then the result will be
    a tropical hypersurface embedded in a real space `\RR^n`::

        sage: T = TropicalSemiring(QQ)
        sage: R.<w,x,y,z> = PolynomialRing(T)
        sage: p1 = x*y + R(-1/2)*x*z + R(4)*z^2 + w*x
        sage: tv = p1.tropical_variety(); tv
        Tropical hypersurface of 0*w*x + 0*x*y + (-1/2)*x*z + 4*z^2
        sage: tv.components()
        [[(t1, t2, t3 - 1/2, t3), [t2 - 9/2 <= t3, t3 <= t1 + 1/2, t2 - 5 <= t1], 1],
        [(t1, 2*t2 - t3 + 4, t3, t2), [t3 + 1/2 <= t2, t3 <= t1], 1],
        [(t1, t2, t1, t3), [max(t1 + 1/2, 1/2*t1 + 1/2*t2 - 2) <= t3], 1],
        [(t1, t2 + 9/2, t3, t2), [t2 <= min(t3 + 1/2, t1 + 1/2)], 1],
        [(t1 - 1/2, t2, t3, t1), [t2 - 9/2 <= t1, t1 <= t3 + 1/2, t2 - 5 <= t3], 1],
        [(2*t1 - t2 + 4, t2, t3, t1), [t1 <= min(1/2*t2 + 1/2*t3 - 2, t2 - 9/2)], 1]]
    """
    def __init__(self, poly):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R.<x,y> = PolynomialRing(T)
            sage: tv = (x+y).tropical_variety()
            sage: TestSuite(tv).run()

        TESTS::

            sage: from sage.rings.semirings.tropical_variety import TropicalVariety
            sage: R.<x,y> = QQ[]
            sage: p1 = x + y
            sage: TropicalVariety(p1)
            Traceback (most recent call last):
            ...
            ValueError: x + y is not a multivariate tropical polynomial
        """
        import operator
        from itertools import combinations
        from sage.symbolic.ring import SR
        from sage.symbolic.relation import solve
        from sage.arith.misc import gcd
        from sage.rings.semirings.tropical_mpolynomial import TropicalMPolynomial

        if not isinstance(poly, TropicalMPolynomial):
            raise ValueError(f"{poly} is not a multivariate tropical polynomial")

        self._poly = poly
        self._hypersurface = []
        tropical_roots = []
        variables = []
        for name in poly.parent().variable_names():
            variables.append(SR.var(name))

        # Convert each term to its linear function
        linear_eq = {}
        pd = poly.dict()
        for key in pd:
            eq = sum(variables[i] * e for i, e in enumerate(key))
            eq += pd[key].lift()
            linear_eq[key] = eq
        temp_keys = []
        temp_order = []

        # Checking for all possible combinations of two terms
        for keys in combinations(pd, 2):
            sol = solve(linear_eq[keys[0]] == linear_eq[keys[1]], variables)

            # Parametric solution of the chosen two terms
            final_sol = []
            for s in sol[0]:
                final_sol.append(s.right())
            xy_interval = []
            xy_interval.append(tuple(final_sol))

            # Comparing with other terms
            min_max = linear_eq[keys[0]]
            for i, v in enumerate(variables):
                min_max = min_max.subs(**{str(v): final_sol[i]})
            all_sol_compare = []
            no_solution = False
            for compare in pd:
                if compare not in keys:
                    temp_compare = linear_eq[compare]
                    for i, v in enumerate(variables):
                        temp_compare = temp_compare.subs(**{str(v): final_sol[i]})
                    if min_max == temp_compare:
                        sol_compare = [[]]
                    elif poly.parent().base()._use_min:
                        sol_compare = solve(min_max < temp_compare, variables)
                    else:
                        sol_compare = solve(min_max > temp_compare, variables)
                    if sol_compare:
                        if isinstance(sol_compare[0], list):
                            if sol_compare[0]:
                                all_sol_compare.append(sol_compare[0][0])
                        else:  # solution is unbounded on one side
                            all_sol_compare.append(sol_compare[0])
                    else:
                        no_solution = True
                        break

            # Solve the condition for parameter
            if not no_solution:
                parameter = set()
                for sol in all_sol_compare:
                    parameter = parameter.union(set(sol.variables()))
                parameter_solution = solve(all_sol_compare, list(parameter))
                if parameter_solution:
                    xy_interval.append(parameter_solution[0])
                    tropical_roots.append(xy_interval)
                    # Calculate the order
                    index_diff = []
                    for i in range(len(keys[0])):
                        index_diff.append(abs(keys[0][i] - keys[1][i]))
                    order = gcd(index_diff)
                    temp_order.append(order)
                    temp_keys.append(keys)

        # Changing all the operator's symbol to <= or >=
        self._keys = []
        components = []
        dim_param = 0
        if tropical_roots:
            dim_param = len(tropical_roots[0][0]) - 1
        vars = [SR.var('t{}'.format(i)) for i in range(1, dim_param+1)]
        for arg in tropical_roots:
            subs_dict = {}
            index_vars = 0
            new_eq = []
            for eq in arg[0]:
                var_eq = eq.variables()
                for var in var_eq:
                    if var not in subs_dict:
                        subs_dict[var] = vars[index_vars]
                        index_vars += 1
                new_eq.append(eq.subs(subs_dict))
            new_eq = tuple(new_eq)
            arg.remove(arg[0])
            arg.insert(0, new_eq)
            if not arg[1] or not isinstance(arg[1], list):
                arg[1] = []
                for var in vars:
                    expr1 = -infinity < var
                    expr2 = var < infinity
                    arg[1].append(expr1)
                    arg[1].append(expr2)
            else:
                params = arg[1]
                arg.remove(params)
                new_param = []
                for param in params:
                    lhs = param.lhs().subs(subs_dict)
                    rhs = param.rhs().subs(subs_dict)
                    if param.operator() == operator.gt:
                        expr = lhs >= rhs
                    else:
                        expr = lhs <= rhs
                    new_param.append(expr)
                arg.insert(1, new_param)
            components.append(arg)

        # Determine the order of each component
        self._vars = vars
        final_order = []
        for i, component in enumerate(components):
            if component not in self._hypersurface:
                self._hypersurface.append(component)
                final_order.append(temp_order[i])
                self._keys.append(temp_keys[i])
            else:
                index = self._hypersurface.index(component)
                if temp_order[i] > final_order[index]:
                    final_order[index] = temp_order[i]
                    self._keys[index] = temp_keys[i]
        for i in range(len(self._hypersurface)):
            self._hypersurface[i].append(final_order[i])

    def dimension(self):
        """
        Return the dimension of ``self``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R.<a,x,y,z> = PolynomialRing(T)
            sage: p1 = x*y + R(-1)*x*z
            sage: p1.tropical_variety().dimension()
            4
        """
        return self._poly.parent().ngens()

    def number_of_components(self):
        """
        Return the number of components that make up ``self``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R.<a,x,y,z> = PolynomialRing(T)
            sage: p1 = x*y*a + x*z + y^2 + a*x + y + z
            sage: p1.tropical_variety().number_of_components()
            13
        """
        from sage.rings.integer_ring import ZZ
        return ZZ(len(self._hypersurface))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R.<w,x,y,z> = PolynomialRing(T)
            sage: (w).tropical_variety()
            Tropical hypersurface of 0*w
        """
        return f"Tropical hypersurface of {self._poly}"

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R.<w,x,y,z> = PolynomialRing(T)
            sage: tv = (R(1)*w^2 + x*y*z + R(-1)).tropical_variety()
            sage: latex(tv)
            TV\left(0 x y z + 1 w^{2} + \left(-1\right)\right)
        """
        return f"TV\\left({self._poly._latex_()}\\right)"

    def components(self):
        """
        Return all components of ``self``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R.<a,x,y,z> = PolynomialRing(T)
            sage: tv = (a+x+y+z).tropical_variety()
            sage: tv.components()
            [[(t1, t1, t2, t3), [t1 <= min(t3, t2)], 1],
             [(t1, t2, t1, t3), [t1 <= t3, t1 <= t2], 1],
             [(t1, t2, t3, t1), [t1 <= min(t3, t2)], 1],
             [(t1, t2, t2, t3), [t2 <= t3, t2 <= t1], 1],
             [(t1, t2, t3, t2), [t2 <= min(t3, t1)], 1],
             [(t1, t2, t3, t3), [t3 <= min(t1, t2)], 1]]
        """
        return self._hypersurface

    def _components_intersection(self):
        r"""
        Return the intersection of three or more components of ``self``.

        For a tropical variety in `\RR^n`, the intersection is characterized
        by a linear equation in `\RR^{n-1}`. Specifically, this becomes a
        vertex for tropical curve and an edges for tropical surface.

        OUTPUT:

        A dictionary where the keys represent component indices and the
        values are lists of tuples. Each tuple contains a parametric
        equation of points and the corresponding parameter's condition.

        EXAMPLES:

        In two dimension, it will provide vertices that are incident with
        each component::

            sage: T = TropicalSemiring(QQ)
            sage: R.<x,y> = PolynomialRing(T)
            sage: p1 = R(2)*x^2 + x*y + R(2)*y^2 + x + R(-1)*y + R(3)
            sage: tv = p1.tropical_variety()
            sage: tv._components_intersection()
            {0: [((-2, 0), {})],
            1: [((-2, 0), {})],
            2: [((-1, -3), {})],
            3: [((-2, 0), {}), ((-1, 0), {})],
            4: [((-1, -3), {}), ((-1, 0), {})],
            5: [((-1, -3), {})],
            6: [((-1, 0), {}), ((3, 4), {})],
            7: [((3, 4), {})],
            8: [((3, 4), {})]}

        In three dimensions, it will provide all parametric equations of
        lines that lie within each component::

            sage: T = TropicalSemiring(QQ)
            sage: R.<x,y,z> = PolynomialRing(T)
            sage: p1 = x + y + z + x^2
            sage: tv = p1.tropical_variety()
            sage: tv._components_intersection()
            {0: [((t2, t2, t2), {0 <= t2}), ((0, 0, t2), {0 <= t2})],
            1: [((0, t2, 0), {0 <= t2}), ((t2, t2, t2), {0 <= t2})],
            2: [((0, t1, 0), {0 <= t1}), ((0, 0, t2), {0 <= t2})],
            3: [((t1, t1, t1), {0 <= t1}), ((t1, 2*t1, 2*t1), {t1 <= 0})],
            4: [((1/2*t2, t2, t2), {t2 <= 0}), ((0, 0, t2), {0 <= t2})],
            5: [((0, t2, 0), {0 <= t2}), ((1/2*t2, t2, t2), {t2 <= 0})]}
        """
        import operator
        from sage.functions.min_max import max_symbolic, min_symbolic
        from sage.symbolic.relation import solve
        from sage.symbolic.expression import Expression
        from sage.sets.set import Set

        def update_result(result):
            sol_param = solve(new_expr, vars)
            sol_param_sim = set()
            for sol in sol_param:
                if sol == []:
                    for v in vars:
                        if v != var:
                            sol_param_sim.add(v < infinity)
                elif isinstance(sol, list):
                    for eqn in sol:
                        if eqn.operator() == operator.eq:
                            if not eqn.rhs().is_numeric():
                                eqn_var = eqn.rhs().variables()
                                param_var = [v for v in eqn_var if v in vars]
                                if not param_var:
                                    v = eqn.lhs()
                                    if v != var:
                                        sol_param_sim.add(v < infinity)
                        elif eqn.operator() == operator.lt:
                            sol_param_sim.add(eqn.lhs() <= eqn.rhs())
                        elif eqn.operator() == operator.gt:
                            sol_param_sim.add(eqn.lhs() >= eqn.rhs())
                else:
                    sol_param_sim.add(sol)

            # Checking there are no conditions with the same variables
            # that use the <= and >= operators simultaneously
            unique_sol_param = set()
            temp = [s for s in sol_param_sim]
            op_temp = {i: set(temp[i].operands()) for i in range(len(temp))}
            for s_value in op_temp.values():
                match_keys = [k for k, v in op_temp.items() if v == s_value]
                if len(match_keys) == 1:
                    for i in match_keys:
                        unique_sol_param.add(temp[i])

            if (unique_sol_param) or (self.dimension() == 2):
                if not unique_sol_param:
                    unique_sol_param = Set()
                if index not in result:
                    result[index] = [(tuple(points), unique_sol_param)]
                else:
                    result[index].append((tuple(points), unique_sol_param))

        result = {}
        vars = self._vars
        for index, comp in enumerate(self._hypersurface):
            for expr in comp[1]:
                left = expr.lhs()
                right = expr.rhs()
                # If the lhs contains a min or max operator
                if (left.operator() == max_symbolic) or (left.operator() == min_symbolic):
                    for operand in expr.lhs().operands():
                        points = list(comp[0])
                        new_expr = [e.subs(**{str(right): operand}) for e in comp[1]]
                        for i, p in enumerate(points):
                            new_eq = p.subs(**{str(right): operand})
                            points[i] = new_eq
                        update_result(result)
                # If the rhs contains a min or max operator
                elif (right.operator() == max_symbolic) or (right.operator() == min_symbolic):
                    for operand in expr.rhs().operands():
                        points = list(comp[0])
                        new_expr = [e.subs(**{str(left): operand}) for e in comp[1]]
                        for i, p in enumerate(points):
                            new_eq = p.subs(**{str(left): operand})
                            points[i] = new_eq
                        update_result(result)
                else:
                    var = expr.variables()[0]
                    points = list(comp[0])
                    subs_expr = solve(left == right, var)[0].rhs()
                    new_expr = [e.subs(**{str(var): subs_expr}) for e in comp[1]]
                    for i, p in enumerate(points):
                        new_eq = p.subs(**{str(var): subs_expr})
                        points[i] = new_eq
                    update_result(result)
        return result

    def dual_subdivision(self):
        """
        Return the dual subdivision of ``self``.

        Dual subdivision refers to a specific decomposition of the
        Newton polygon of a tropical polynomial. This Newton polygon
        is the convex hull of all the points corresponding to the
        exponents of the terms of the tropical polynomial. The term
        "dual" is used in the sense that the combinatorial structure
        of the tropical variety is reflected in the dual subdivision.
        Vertices of the dual subdivision correspond to the intersection
        of multiple components. Edges of the dual subdivision correspond
        to the individual components.

        OUTPUT: :class:`sage.geometry.polyhedral_complex.PolyhedralComplex`

        EXAMPLES:

        Dual subdivision of a tropical curve::

            sage: T = TropicalSemiring(QQ, use_min=False)
            sage: R.<x,y> = PolynomialRing(T)
            sage: p1 = R(3) + R(2)*x + R(2)*y + R(3)*x*y + x^2 + y^2
            sage: tv = p1.tropical_variety()
            sage: pc = tv.dual_subdivision()
            sage: pc.plot()
            Graphics object consisting of 20 graphics primitives

        .. PLOT::
            :width: 300 px

            T = TropicalSemiring(QQ,  use_min=False)
            R = PolynomialRing(T, ('x,y'))
            x, y = R.gen(), R.gen(1)
            p1 = R(3) + R(2)*x + R(2)*y + R(3)*x*y + x**2 + y**2
            tv = p1.tropical_variety()
            pc = tv.dual_subdivision()
            sphinx_plot(pc.plot())

        Dual subdivision of a tropical surface::

            sage: T = TropicalSemiring(QQ)
            sage: R.<x,y,z> = PolynomialRing(T)
            sage: p1 = x + y + z + x^2 + R(1)
            sage: tv = p1.tropical_variety()
            sage: pc = tv.dual_subdivision()
            sage: pc.plot()
            Graphics3d Object

        .. PLOT::
            :width: 300 px

            T = TropicalSemiring(QQ,  use_min=False)
            R = PolynomialRing(T, ('x,y,z'))
            x, y, z = R.gen(), R.gen(1), R.gen(2)
            p1 = x + y + z + x**2 + R(1)
            tv = p1.tropical_variety()
            pc = tv.dual_subdivision()
            sphinx_plot(pc.plot())

        Dual subdivision of a tropical hypersurface::

            sage: T = TropicalSemiring(QQ)
            sage: R.<a,b,c,d> = PolynomialRing(T)
            sage: p1 = a^2 + b^2 + c^2 + d^2 + a*b*c*d
            sage: tv = p1.tropical_variety()
            sage: pc = tv.dual_subdivision(); pc
            Polyhedral complex with 6 maximal cells
        """
        from sage.graphs.graph import Graph
        from sage.geometry.polyhedron.constructor import Polyhedron
        from sage.geometry.polyhedral_complex import PolyhedralComplex

        G = Graph()
        edges = [e for e in self._keys]
        G.add_edges(edges)

        polyhedron_lst = []
        for cycle in G.cycle_basis():
            polyhedron = Polyhedron(vertices=cycle)
            polyhedron_lst.append(polyhedron)
        pc = PolyhedralComplex(polyhedron_lst)
        return pc

    def weight_vectors(self):
        r"""
        Return the weight vectors for each unique intesection of
        components of ``self``.

        Weight vectors are a list of vectors associated with each
        unique intersection of the components of tropical variety.
        Each vector is a normal vector to a component with respect
        to the unique intersection lying within that component.

        Assume ``self`` is a `n`-dimensional tropical variety.
        Suppose `L` is an intersection lying within the components
        `S_1, ldots, S_k` with respective weights `w_1, ldots, w_k`.
        This `L` is a linear structure in `\RR^{n-1}` and has `n-1`
        direction vectors `d_1,d_2,\dots, d_{n-1}`. Each component
        `S_1, ldots, S_k` has a normal vector `n_1, \ldots, n_k`.
        Then, we scale each normal vector to an integer vector such
        that the greatest common divisor of its elements is 1.

        The weight vector of a component `S_i` with respect to `L`
        can be found by calculating the cross product between direction
        vectors of `L` and normal vector `n_i`.These vectors will
        satisfy the balancing condition `\sum_{i=1}^k w_k v_k = 0`.

        OUTPUT:

        A tuple of two dictionaries. The first dictionary contains
        equations representing the intersections. The second dictionary
        contains lists of vectors.

        EXAMPLES:

        Weight vectors of tropical surface::

            sage: T = TropicalSemiring(QQ)
            sage: R.<x,y,z> = PolynomialRing(T)
            sage: p = x^2 + R(-1)*y + z + R(1)
            sage: tv = p.tropical_variety()
            sage: tv.weight_vectors()
            ({0: ((1/2*u2, u2 + 1, u2), {u2 <= 1}),
             1: ((1/2, 2, u2), {1 <= u2}),
             2: ((1/2, u2, 1), {2 <= u2}),
             3: ((u1, 2, 1), {(1/2) <= u1})},
            {0: [(1, 2, -5/2), (1, -5/2, 2), (-2, 1/2, 1/2)],
             1: [(-1, -2, 0), (0, 2, 0), (1, 0, 0)],
             2: [(1, 0, 2), (0, 0, -2), (-1, 0, 0)],
             3: [(0, 1, 1), (0, 0, -1), (0, -1, 0)]})

        Weight vectors of tropical hypersurface::

            sage: T = TropicalSemiring(QQ)
            sage: R.<a,b,c,d> = PolynomialRing(T)
            sage: p1 = R(2)*a*b + R(3)*a*c + R(-1)*c^2 + R(-1/3)*a*d
            sage: tv = p1.tropical_variety()
            sage: tv.weight_vectors()
            ({0: ((u1, u3 - 7/3, u3 - 10/3, u3), {u1 <= u3 - 22/3}),
             1: ((u2 - 4, u2 + 1, u2, u3), {u2 <= u3 - 10/3}),
             2: ((2*u1 - u3 - 2/3, u3 - 7/3, u1, u3), {u3 - 10/3 <= u1}),
             3: ((u3 - 22/3, u2, u3 - 10/3, u3), {u3 - 7/3 <= u2})},
            {0: [(0, 1, 1, -2), (0, 1, -2, 1), (0, -2, 1, 1)],
             1: [(-2, 1, 1, 0), (3, -3, 0, 0), (-1, 2, -1, 0)],
             2: [(-1, 5, 2, -6), (2, 1, -4, 1), (-1, -6, 2, 5)],
             3: [(-1, 0, -1, 2), (-2, 0, 1, 1), (3, 0, 0, -3)]})
        """
        from sage.symbolic.ring import SR
        from sage.symbolic.relation import solve
        from sage.calculus.functional import diff
        from sage.arith.misc import gcd
        from sage.matrix.constructor import matrix
        from sage.modules.free_module_element import vector, zero_vector

        dim = self.dimension()
        t = SR.var('t')
        t_vars = [SR.var('t{}'.format(i)) for i in range(dim)]
        u_vars = [SR.var('u{}'.format(i)) for i in range(dim)]
        convert_tu = {ti: ui for ti, ui in zip(t_vars, u_vars)}
        CI = self._components_intersection()
        unique_line = set()
        index_line = {}
        line_comps = {}
        index = 0

        # Find the unique intersection between multiple components and
        # the indices of the components containing this intersection.
        for i, lines in CI.items():
            for line in lines:
                eqn = line[0]
                is_unique = True
                for uniq in unique_line:
                    subs_index = -1
                    for j in range(dim):
                        if eqn[j] != uniq[j]:
                            subs_index = j
                            break
                    if subs_index == -1:
                        new_line = eqn
                        is_unique = False
                        break
                    subs_dict = {}
                    while len(subs_dict) != dim-2 and subs_index < dim:
                        eq1 = eqn[subs_index].subs(subs_dict)
                        vib = None
                        for unk in eq1.variables():
                            if unk not in subs_dict:
                                if unk in t_vars:
                                    vib = unk
                                    break
                        if vib:
                            eq1 = eq1.subs(**{str(vib): t})
                            eq2 = uniq[subs_index]
                            temp_sol = solve(eq1 == eq2, t)
                            if temp_sol:
                                temp_sol = temp_sol[0].rhs()
                                if not temp_sol.is_numeric():
                                    subs_dict[vib] = temp_sol
                        subs_index += 1
                    if subs_dict:
                        new_line = []
                        for l in eqn:
                            for key, value in subs_dict.items():
                                l = l.subs(key == value)
                            new_line.append(l)
                        if tuple(new_line) in unique_line:
                            is_unique = False
                            break
                if is_unique:
                    new_eqn = [eq.subs(convert_tu) for eq in eqn]
                    new_eqn = tuple(new_eqn)
                    cdns = line[1]
                    new_cdn = [cdn.subs(convert_tu) for cdn in cdns]
                    new_cdn = set(new_cdn)
                    unique_line.add(new_eqn)
                    index_line[index] = tuple([new_eqn, new_cdn])
                    line_comps[index] = [i]
                    index += 1
                else:
                    match_key = [k for k, v in index_line.items() if v[0] == tuple(new_line)][0]
                    line_comps[match_key].append(i)

        WV = {i: [] for i in range(len(line_comps))}
        for k, index in line_comps.items():

            # Calculate direction vector of the line
            dir_vecs = []
            line = index_line[k][0]
            all_var = set()
            for l in line:
                for v in l.variables():
                    all_var.add(v)
            for vpar in all_var:
                par_drv = []
                for l in line:
                    par_drv.append(QQ(diff(l, vpar)))
                par_drv = vector(par_drv)
                dir_vecs.append(par_drv)

            # Calculate the outgoing normal vector of each surface in the
            # direction of the line
            for i in index:
                surface = self._hypersurface[i][0]
                drv_vectors = []
                for vpar in self._vars:
                    temp_vec = []
                    for s in surface:
                        temp_vec.append(QQ(diff(s, vpar)))
                    temp_vec = vector(temp_vec)
                    drv_vectors.append(temp_vec)
                temp = [t_vars]
                for vec in drv_vectors:
                    temp.append(vec)
                vec_matrix = matrix(SR, temp)
                normal_vec = vec_matrix.det()
                temp_nor = []
                for tvar in t_vars:
                    temp_nor.append(QQ(diff(normal_vec, tvar)))
                normal_vec = vector(temp_nor)
                normal_vec *= 1/gcd(normal_vec)

                # Calculate the weight vector
                temp_final = [t_vars]
                for v in dir_vecs:
                    temp_final.append(v)
                temp_final.append(normal_vec)
                vec_matrix = matrix(SR, temp_final)
                weight_vec = vec_matrix.det()
                temp_weight = []
                for tvar in t_vars:
                    temp_weight.append(QQ(diff(weight_vec, tvar)))
                weight_vec = vector(temp_weight)
                order = self._hypersurface[i][2]
                weight_vec *= order
                WV[k].append(weight_vec)

            for i in range(len(WV[k])):
                test_vectors = [v for v in WV[k]]
                test_vectors[i] = -test_vectors[i]
                if sum(test_vectors) == zero_vector(QQ, dim):
                    WV[k] = test_vectors
                    break

        return index_line, WV


class TropicalSurface(TropicalVariety):
    r"""
    A tropical surface in `\RR^3`.

    The tropical surface consists of planar regions and facets, which we
    can call cells. These cells are connected in such a way that they form
    a piecewise linear structure embedded in three-dimensional space. These
    cells meet along edges, where the balancing condition is satisfied.
    This balancing condition ensures that the sum of the outgoing normal
    vectors at each edge is zero, reflecting the equilibrium.

    EXAMPLES::

        sage: T = TropicalSemiring(QQ, use_min=False)
        sage: R.<x,y,z> = PolynomialRing(T)
        sage: p1 = x + y + z + R(0)
        sage: tv = p1.tropical_variety(); tv
        Tropical surface of 0*x + 0*y + 0*z + 0
        sage: tv.components()
        [[(t1, t1, t2), [t2 <= t1, 0 <= t1], 1],
        [(t1, t2, t1), [max(0, t2) <= t1], 1],
        [(0, t1, t2), [t2 <= 0, t1 <= 0], 1],
        [(t1, t2, t2), [max(0, t1) <= t2], 1],
        [(t1, 0, t2), [t2 <= 0, t1 <= 0], 1],
        [(t1, t2, 0), [t1 <= 0, t2 <= 0], 1]]
    """
    def _axes(self):
        r"""
        Set the default axes for ``self``.

        This default axes is used for the 3d plot. The axes is centered
        around where the intersection of the components occured so it
        gives a nice visual representation for the interactions between
        different components of the surface. Additionally, it enhances
        the visibility and interpretation of how the components align
        and interact in three-dimensional space.

        OUTPUT:

        A list of three lists, where the first inner list represent value
        of x-axis, the second inner list represent value of y-axis, and
        the third inner list represent value of z-axis. If there are
        either no components or only one component, the axis will be set
        to `[[-1, 1], [-1, 1], [-1, 1]]`.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R.<x,y,z> = PolynomialRing(T)
            sage: p1 = x
            sage: p1.tropical_variety()._axes()
            [[-1, 1], [-1, 1], [-1, 1]]
            sage: p2 = x + y + z + x^2 + R(1)
            sage: p2.tropical_variety()._axes()
            [[-1, 2], [-1, 2], [-1, 2]]
        """
        from sage.symbolic.relation import solve
        from sage.arith.srange import srange

        if not self._hypersurface:
            return [[-1, 1], [-1, 1], [-1, 1]]
        elif len(self._hypersurface) == 1:
            bound = 1
            for eqn in self._hypersurface[0][0]:
                for op in eqn.operands():
                    if op.is_numeric():
                        if op > bound:
                            bound = op
            return [[-bound, bound]] * 3

        u_set = set()
        v_set = set()
        for comp in self._hypersurface:
            list_expr = []
            temp_u = set()
            temp_v = set()
            for expr in comp[1]:
                if expr.lhs().is_numeric():
                    if bool(expr.rhs() == self._vars[0]):
                        temp_u.add(expr.lhs())
                    else:
                        temp_v.add(expr.lhs())
                elif expr.rhs().is_numeric():
                    if bool(expr.lhs() == self._vars[0]):
                        temp_u.add(expr.rhs())
                    else:
                        temp_v.add(expr.rhs())
                else:
                    list_expr.append(expr)
            if not temp_u:
                temp_u.add(0)
            if not temp_v:
                temp_v.add(0)
            for expr in list_expr:
                for u in temp_u:
                    sol = solve(expr.subs(**{str(self._vars[0]): u}), self._vars[1])
                    if not sol:
                        temp_v.add(0)
                    elif not sol[0]:
                        temp_v.add(0)
                    else:
                        temp_v.add(sol[0][0].rhs())
                for v in temp_v:
                    sol = solve(expr.subs(**{str(self._vars[1]): v}), self._vars[0])
                    if not sol:
                        temp_u.add(0)
                    elif not sol[0]:
                        temp_u.add(0)
                    else:
                        temp_u.add(sol[0][0].rhs())
            u_set = u_set.union(temp_u)
            v_set = v_set.union(temp_v)
        axes = [[min(u_set)-1, max(u_set)+1], [min(v_set)-1, max(v_set)+1]]

        # Calculate the z-axis
        step = 10
        du = (axes[0][1]-axes[0][0]) / step
        dv = (axes[1][1]-axes[1][0]) / step
        u_range = srange(axes[0][0], axes[0][1]+du, du)
        v_range = srange(axes[1][0], axes[1][1]+dv, dv)
        zmin, zmax = None, None
        for comp in self._hypersurface:
            for u in u_range:
                for v in v_range:
                    checkpoint = True
                    for exp in comp[1]:
                        final_exp = exp.subs(**{str(self._vars[0]): u, str(self._vars[1]): v})
                        if not final_exp:
                            checkpoint = False
                            break
                    if checkpoint:
                        z = comp[0][2].subs(**{str(self._vars[0]): u, str(self._vars[1]): v})
                        if (zmin is None) and (zmax is None):
                            zmin = z
                            zmax = z
                        else:
                            if z < zmin:
                                zmin = z
                            if z > zmax:
                                zmax = z
        axes.append([zmin, zmax])
        return axes

    def _polygon_vertices(self):
        r"""
        Return the vertices of the polygon for each components of ``self``
        to be used for plotting.

        OUTPUT:

        A dictionary where the keys represent component indices and the
        values are a set of points in three dimensional space.

        EXAMPLES:

        A tropical surface with only one component::

            sage: T = TropicalSemiring(QQ)
            sage: R.<x,y,z> = PolynomialRing(T)
            sage: p1 = x + z
            sage: tv1 = p1.tropical_variety()
            sage: tv1._polygon_vertices()
            {0: {(-1, -1, -1), (-1, 1, -1), (1, -1, 1), (1, 1, 1)}}

        A tropical surface with multiple components::

            sage: p2 = x^2 + x + y + z + R(1)
            sage: tv2 = p2.tropical_variety()
            sage: tv2._polygon_vertices()
            {0: {(0, 0, 0), (0, 0, 2), (1, 1, 1), (2, 2, 2)},
            1: {(0, 0, 0), (0, 2, 0), (1, 1, 1), (2, 2, 2)},
            2: {(0, 0, 0), (0, 0, 2), (0, 2, 0), (0, 2, 2)},
            3: {(1, 1, 1), (1, 1, 2), (1, 2, 1), (1, 2, 2)},
            4: {(-1/2, -1, -1), (0, 0, 0), (1, 1, 1), (2, -1, -1), (2, 1, 1)},
            5: {(-1/2, -1, -1), (-1/2, -1, 2), (0, 0, 0), (0, 0, 2)},
            6: {(1, 1, 1), (1, 1, 2), (2, 1, 1), (2, 1, 2)},
            7: {(-1/2, -1, -1), (-1/2, 2, -1), (0, 0, 0), (0, 2, 0)},
            8: {(1, 1, 1), (1, 2, 1), (2, 1, 1), (2, 2, 1)}}
        """
        from sage.symbolic.relation import solve
        from sage.sets.real_set import RealSet

        poly_verts = {i: set() for i in range(self.number_of_components())}
        axes = self._axes()
        comps = self.components()
        vars = self._vars
        comps_int = self._components_intersection()

        # Find the inside vertices (intersection of components)
        for index, lines in comps_int.items():
            for line in lines:
                v = list(line[1])[0].variables()[0]
                for param in line[1]:
                    left = param.lhs()
                    right = param.rhs()
                    if left.is_numeric():
                        vertex = [QQ(e.subs(**{str(v): left})) for e in line[0]]
                        poly_verts[index].add(tuple(vertex))
                    elif right.is_numeric():
                        vertex = [QQ(e.subs(**{str(v): right})) for e in line[0]]
                        poly_verts[index].add(tuple(vertex))

        def find_edge_vertices(i):
            j = (i+1) % 2
            if i == 0:  # interval for t1
                interval = interval1
            else:  # interval for t2
                interval = interval2
            for p in [interval.inf(), interval.sup()]:
                new_param = [e.subs(**{str(vars[i]): p}) for e in comps[index][1]]
                sol = solve(new_param, vars[j])
                if sol:
                    interval_param = RealSet()
                    for s in sol:
                        if s != []:
                            # Handle cases where 's' is not a list (s = r1 < +Infinity),
                            # or if 's' is a list, ensure that its elements define a real
                            # interval (catch invalid cases like s = [t1 == r100]).
                            try:
                                interval_param += RealSet(s[0])
                            except (IndexError, ValueError):
                                interval_param += RealSet(-infinity, infinity)
                        else:
                            interval_param += RealSet(-infinity, infinity)
                    interval_param = interval_param.intersection(interval2)
                    if is_doublevar:
                        int1 = RealSet()
                        for s1 in sol1:
                            subs1 = solve(s1[0].subs(**{str(vars[i]): p}), vars[j])
                            try:
                                int1 += RealSet(subs1[0])
                            except TypeError:
                                int1 += RealSet(subs1[0][0])
                        int2 = RealSet()
                        for s2 in sol2:
                            subs2 = solve(s2[0].subs(**{str(vars[i]): p}), vars[j])
                            try:
                                int2 += RealSet(subs2[0])
                            except TypeError:
                                int2 += RealSet(subs2[0][0])
                        final_int = int1.intersection(int2)
                        interval_param = interval_param.intersection(final_int)
                    if interval_param:
                        vertex1 = [QQ(e.subs(**{str(vars[i]): p, str(vars[j]): interval_param.inf()})) for e in comps[index][0]]
                        vertex2 = [QQ(e.subs(**{str(vars[i]): p, str(vars[j]): interval_param.sup()})) for e in comps[index][0]]
                        poly_verts[index].add(tuple(vertex1))
                        poly_verts[index].add(tuple(vertex2))

        # Find the interval of parameter for outer vertex
        for index in range(len(comps)):
            interval1 = RealSet(-infinity,infinity)  # represent t1
            interval2 = RealSet(-infinity,infinity)  # represent t2
            is_doublevar = False
            for i, point in enumerate(comps[index][0]):
                pv = point.variables()
                if len(pv) == 1:
                    temp1 = RealSet(solve(point >= axes[i][0], pv[0])[0][0])
                    temp2 = RealSet(solve(point <= axes[i][1], pv[0])[0][0])
                    temp = temp1.intersection(temp2)
                    if pv[0] == vars[0]:
                        interval1 = interval1.intersection(temp)
                    else:
                        interval2 = interval2.intersection(temp)
                elif len(pv) == 2:
                    sol1 = solve(point >= axes[i][0], pv)
                    sol2 = solve(point <= axes[i][1], pv)
                    is_doublevar = True

            # Find the edge vertices (those that touch the axes)
            find_edge_vertices(0)  # t1 fixed
            find_edge_vertices(1)  # t2 fixed
        return poly_verts

    def plot(self, color='random'):
        """
        Return the plot of ``self`` by constructing a polyhedron from
        vertices in ``self.polygon_vertices()``.

        INPUT:

        - ``color`` -- string or tuple that represent a color (default:
          ``random``); ``random`` means each polygon will be assigned
          a different color. If instead a specific ``color`` is provided,
          then all polygon will be given the same color.

        OUTPUT: Graphics3d Object

        EXAMPLES:

        A tropical surface that consist of only one cell::

            sage: T = TropicalSemiring(QQ)
            sage: R.<x,y,z> = PolynomialRing(T)
            sage: p1 = x + z
            sage: tv = p1.tropical_variety()
            sage: tv.plot()
            Graphics3d Object

        .. PLOT::
            :width: 300 px

            T = TropicalSemiring(QQ)
            R = PolynomialRing(T, ('x,y,z'))
            x, y, z = R.gen(), R.gen(1), R.gen(2)
            p1 = x + z
            sphinx_plot(p1.tropical_variety().plot())

        A tropical surface with multiple cells that exhibit complex and
        intriguing geometric structures::

            sage: p2 = x^2 + x + y + z + R(1)
            sage: tv = p2.tropical_variety()
            sage: tv.plot()
            Graphics3d Object

        .. PLOT::
            :width: 300 px

            T = TropicalSemiring(QQ)
            R = PolynomialRing(T, ('x,y,z'))
            x, y, z = R.gen(), R.gen(1), R.gen(2)
            p2 = x**2 + x + y + z + R(1)
            sphinx_plot(p2.tropical_variety().plot())
        """
        from random import random
        from sage.plot.graphics import Graphics
        from sage.geometry.polyhedron.constructor import Polyhedron

        if color == 'random':
            colors = []
            for _ in range(self.number_of_components()):
                color = (random(), random(), random())
                colors.append(color)
        elif isinstance(color, str):
            colors = [color] * self.number_of_components()
        else:
            colors = color

        combined_plot = Graphics()
        for i, vertex in self._polygon_vertices().items():
            points = list(vertex)
            plot = Polyhedron(vertices=points).plot(color=colors[i])
            combined_plot += plot
        return combined_plot

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R.<x,y,z> = PolynomialRing(T)
            sage: (x^4+z^2).tropical_variety()
            Tropical surface of 0*x^4 + 0*z^2
        """
        return f"Tropical surface of {self._poly}"


class TropicalCurve(TropicalVariety):
    r"""
    A tropical curve in `\RR^2`.

    The tropical curve consists of line segments and half-lines, which we
    call edges. These edges are connected in such a way that they form a
    piecewise linear graph embedded in the plane. These edges meet at
    a vertices, where the balancing condition is satisfied. This balancing
    condition ensures that the sum of the outgoing slopes at each vertex
    is zero, reflecting the equilibrium.

    EXAMPLES::

        sage: T = TropicalSemiring(QQ, use_min=False)
        sage: R.<x,y> = PolynomialRing(T)
        sage: p1 = x + y + R(0)
        sage: tv = p1.tropical_variety(); tv
        Tropical curve of 0*x + 0*y + 0
        sage: tv.components()
        [[(t1, t1), [t1 >= 0], 1], [(0, t1), [t1 <= 0], 1], [(t1, 0), [t1 <= 0], 1]]
    """
    def _axes(self):
        """
        Set the default axes for ``self``.

        This default axes is used for plot of tropical curve and also the
        3d plot of tropical polynomial function. The axes is chosen by first
        find all vertices of this tropical curve. Then we choose the minimum
        and maximum of all x-component in this vertices to be the x-axis.
        The same apply to the y-axis.

        OUTPUT:

        A list of two lists, where the first inner list represent value of
        x-axis and the second inner list represent value of y-axis.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R.<x,y> = PolynomialRing(T)
            sage: p1 = x^2
            sage: p1.tropical_variety()._axes()
            [[-1, 1], [-1, 1]]
            sage: p2 = R(12)*x*y + R(-2)*y^2 + R(16)*y + R(25)
            sage: p2.tropical_variety()._axes()
            [[-3/2, 1/2], [25/2, 29/2]]
        """
        if self.number_of_components() == 0:
            return [[-1, 1], [-1, 1]]
        if self.number_of_components() <= 2:
            bound = 1
            for comps in self._hypersurface:
                eqns = comps[0]
                temp_operands = []
                for eq in eqns:
                    if not eq.operator():
                        temp_operands.append(eq)
                    else:
                        temp_operands += eq.operands()
                for op in temp_operands:
                    if op.is_numeric():
                        if abs(op) > bound:
                            bound = abs(op)
            return [[-bound, bound]] * 2

        verts = self.vertices()
        xmin = xmax = list(verts)[0][0]
        for vertex in verts:
            if vertex[0] < xmin:
                xmin = vertex[0]
            elif vertex[0] > xmax:
                xmax = vertex[0]
        ymin = ymax = list(verts)[0][1]
        for vertex in verts:
            if vertex[1] < ymin:
                ymin = vertex[1]
            elif vertex[1] > ymax:
                ymax = vertex[1]
        return [[xmin-1, xmax+1], [ymin-1, ymax+1]]

    def vertices(self):
        r"""
        Return all vertices of ``self``, which is the point where three or
        more edges intersect.

        OUTPUT: A set of `(x,y)` points

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R.<x,y> = PolynomialRing(T)
            sage: p1 = x + y
            sage: p1.tropical_variety().vertices()
            set()
            sage: p2 = R(-2)*x^2 + R(-1)*x + R(1/2)*y + R(1/6)
            sage: p2.tropical_variety().vertices()
            {(1, -1/2), (7/6, -1/3)}
        """
        if len(self._hypersurface) < 3:
            return set()
        return set(self._vertices_components().keys())

    def _vertices_components(self):
        """
        Return the index of components adjacent to each vertex of ``self``.

        OUTPUT:

        A dictionary where the keys represent the vertices, and the
        values are lists of tuples. Each tuple consists of the index
        of an adjacent edge (component) `e_i` and a string indicating
        the directionality of `e_i` relative to the vertex. The string
        is either "pos" or "neg", specifying whether it is positive or
        negative.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R.<x,y> = PolynomialRing(T)
            sage: p1 = R(0) + x + y + x*y + x^2*y + x*y^2
            sage: p1.tropical_variety()._vertices_components()
            {(0, 0): [(0, 'pos'), (1, 'pos'), (2, 'pos'), (3, 'neg'), (4, 'neg')]}
            sage: p2 = R(2)*x^2 + x*y + R(2)*y^2 + x + R(-1)*y + R(3)
            sage: p2.tropical_variety()._vertices_components()
            {(-2, 0): [(0, 'neg'), (1, 'pos'), (3, 'pos')],
             (-1, -3): [(2, 'neg'), (4, 'pos'), (5, 'pos')],
             (-1, 0): [(3, 'neg'), (4, 'neg'), (6, 'pos')],
             (3, 4): [(6, 'neg'), (7, 'pos'), (8, 'pos')]}
        """
        comp_vert = {}
        if len(self._hypersurface) >= 3:
            for i, component in enumerate(self._hypersurface):
                parametric_function = component[0]
                v = component[1][0].variables()[0]
                interval = self._parameter_intervals()[i]
                lower = interval[0].lower()
                upper = interval[0].upper()
                if lower != -infinity:
                    x = parametric_function[0].subs(**{str(v): lower})
                    y = parametric_function[1].subs(**{str(v): lower})
                    if (x,y) not in comp_vert:
                        comp_vert[(x,y)] = [(i, 'pos')]
                    else:
                        comp_vert[(x,y)].append((i, 'pos'))
                if upper != infinity:
                    x = parametric_function[0].subs(**{str(v): upper})
                    y = parametric_function[1].subs(**{str(v): upper})
                    if (x,y) not in comp_vert:
                        comp_vert[(x,y)] = [(i, 'neg')]
                    else:
                        comp_vert[(x,y)].append((i, 'neg'))
        return comp_vert

    def weight_vectors(self):
        r"""
        Return the weight vectors for all vertices of ``self``.

        Weight vectors are a list of vectors associated with each vertex
        of the curve. Each vector corresponds to an edge emanating from
        that vertex and points in the direction of the edge.

        Suppose `v` is a vertex adjacent to the edges `e_1, ldots, e_k`
        with respective weights `w_1, ldots, w_k`. Every edge `e_i` is
        contained in a line (component) defined by an equation. Therefore,
        there exists a unique integer vector `v_i=(\alpha, \beta)` in
        the direction of `e_i` such that `\gcd(\alpha, \beta)=1`. Then,
        each vertex `v` yield the vectors `w_1v_1,ldots,w_kv_k`.
        These vectors will satisfy the following balancing condition:
        `\sum_{i=1}^k w_i v_i = 0`.

        OUTPUT:

        A dictionary where the keys represent the vertices, and the values
        are lists of vectors.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R.<x,y> = PolynomialRing(T)
            sage: p1 = R(-2)*x^2 + R(-1)*x + R(1/2)*y + R(1/6)
            sage: p1.tropical_variety().weight_vectors()
            {(1, -1/2): [(0, 1), (-1, -2), (1, 1)],
             (7/6, -1/3): [(-1, -1), (0, 1), (1, 0)]}

            sage: p2 = R(2)*x^2 + x*y + R(2)*y^2 + x + R(-1)*y + R(3)
            sage: p2.tropical_variety().weight_vectors()
            {(-2, 0): [(-1, -1), (0, 1), (1, 0)],
             (-1, -3): [(-1, -1), (0, 1), (1, 0)],
             (-1, 0): [(-1, 0), (0, -1), (1, 1)],
             (3, 4): [(-1, -1), (0, 1), (1, 0)]}
        """
        from sage.calculus.functional import diff
        from sage.modules.free_module_element import vector
        from sage.arith.misc import gcd

        if not self._vertices_components():
            return {}

        # Calculate the base vector in the direction of each edge
        temp_vectors = []
        par = self._hypersurface[0][1][0].variables()[0]
        for comp in self._hypersurface:
            dx = diff(comp[0][0], par)
            dy = diff(comp[0][1], par)
            multiplier = gcd(QQ(dx), QQ(dy))
            temp_vectors.append(vector([dx/multiplier, dy/multiplier]))

        # Calculate the weight vectors of each vertex
        cov = self._vertices_components()
        result = {}
        for vertex in cov:
            vectors = []
            for comp in cov[vertex]:
                weight = self._hypersurface[comp[0]][2]
                if comp[1] == 'pos':
                    vectors.append(weight*temp_vectors[comp[0]])
                else:
                    vectors.append(weight*(-temp_vectors[comp[0]]))
            result[vertex] = vectors
        return result

    def is_smooth(self):
        r"""
        Return ``True`` if ``self`` is smooth and ``False`` otherwise.

        Suppose `C` is a tropical curve of degree `d`. A tropical curve
        `C` is smooth if the dual subdivision of `C` consists of `d^2`
        triangles each having unit area `1/2`. This is equivalent with
        `C` having `d^2` vertices. These vertices are necessarily
        trivalent (has three adjacent edges).

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R.<x,y> = PolynomialRing(T)
            sage: p1 = x^2 + x + R(1)
            sage: p1.tropical_variety().is_smooth()
            False
            sage: p2 = R(2)*x^2 + x*y + R(2)*y^2 + x + R(-1)*y + R(3)
            sage: p2.tropical_variety().is_smooth()
            True
        """
        if len(self.vertices()) == self._poly.degree()**2:
            return True
        return False

    def is_simple(self):
        r"""
        Return ``True`` if ``self`` is simple and ``False`` otherwise.

        A tropical curve `C` is called simple if each vertex is either
        trivalent or is locally the intersection of two line segments.
        Equivalently, `C` is simple if the corresponding subdivision
        consists only of triangles and parallelograms.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R.<x,y> = PolynomialRing(T)
            sage: p1 = R(0) + x + y + x*y + x^2*y + x*y^2
            sage: p1.tropical_variety().is_simple()
            False
            sage: p2 = R(2)*x^2 + x*y + R(2)*y^2 + x + R(-1)*y + R(3)
            sage: p2.tropical_variety().is_simple()
            True
        """
        vov = self.weight_vectors()
        for vertex in self.vertices():
            if len(vov[vertex]) > 4:
                return False
            elif len(vov[vertex]) == 4:
                for v in vov[vertex]:
                    if -v not in vov[vertex]:
                        return False
        return True

    def genus(self):
        r"""
        Return the genus of ``self``.

        Let `t(C)` be the number of trivalent vertices, and let `r(C)` be
        the number of unbounded edges of `C`. The genus of simple tropical
        curve `C` is defined by the formula:
        `g(C) = \frac{1}{2}t(C) - \frac{1}{2}r(C) + 1`.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R.<x,y> = PolynomialRing(T)
            sage: p1 = x^2 + y^2 + x*y
            sage: p1.tropical_variety().genus()
            1
            sage: p2 = R(2)*x^2 + x*y + R(2)*y^2 + x + R(-1)*y + R(3)
            sage: p2.tropical_variety().genus()
            0

        TESTS::

            sage: T = TropicalSemiring(QQ)
            sage: R.<x,y> = PolynomialRing(T)
            sage: p1 = R(0) + y + x^2*y + x*y^2
            sage: p1.tropical_variety().genus()
            Traceback (most recent call last):
            ...
            ValueError: Tropical curve of 0*x^2*y + 0*x*y^2 + 0*y + 0 is not simple
        """
        if not self.is_simple():
            raise ValueError(f"{self} is not simple")
        trivalent = 0  # number of trivalent vertices
        for vectors in self.weight_vectors().values():
            if len(vectors) == 3:
                trivalent += 1
        unbounded = 0  # number of unbounded edges
        for component in self._hypersurface:
            if len(component[1]) == 1:
                unbounded += 1
        return trivalent//2 - unbounded//2 + 1

    def contribution(self):
        r"""
        Return the contribution of ``self``.

        The contribution of a simple curve `C` is defined as the product
        of the normalized areas of all triangles in the corresponding
        dual subdivision. We just multiply positive integers attached to
        the trivalent vertices. The contribution of a trivalent vertex
        equals `w_1w_2|\det(v_1,v_2)|`, with `w_i` are the weights of
        the adjacent edges and `v_i` are their weight vectors. That
        formula is independent of the choice made because of the
        balancing condition `w_1v_1+w_2v_2+w_3v_3=0`.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R.<x,y> = PolynomialRing(T)
            sage: p1 = R(2)*x^2 + x*y + R(2)*y^2 + x + R(-1)*y + R(3)
            sage: p1.tropical_variety().contribution()
            1
            sage: p2 = R(-1/3)*x^2 + R(1)*x*y + R(1)*y^2 + R(-1/3)*x + R(1/3)
            sage: p2.tropical_variety().contribution()
            16

        TESTS::

            sage: T = TropicalSemiring(QQ)
            sage: R.<x,y> = PolynomialRing(T)
            sage: p1 = R(0) + x + x^2*y + x*y^2
            sage: p1.tropical_variety().contribution()
            Traceback (most recent call last):
            ...
            ValueError: Tropical curve of 0*x^2*y + 0*x*y^2 + 0*x + 0 is not simple
        """
        if not self.is_simple():
            raise ValueError(f"{self} is not simple")
        result = 1
        voc = self._vertices_components()
        vov = self.weight_vectors()
        for vertex in vov:
            if len(vov[vertex]) == 3:
                u1 = vov[vertex][0]
                u2 = vov[vertex][1]
                index1 = voc[vertex][0][0]
                index2 = voc[vertex][1][0]
                w1 = self._hypersurface[index1][2]
                w2 = self._hypersurface[index2][2]
                det = u1[0]*u2[1] - u1[1]*u2[0]
                result *= w1 * w2 * abs(det)
        return result

    def _parameter_intervals(self):
        r"""
        Return the intervals of each component's parameter of ``self``.

        OUTPUT: A list of ``RealSet``

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R.<x,y> = PolynomialRing(T)
            sage: p1 = y + y^2
            sage: p1.tropical_variety()._parameter_intervals()
            [(-oo, +oo)]
            sage: p2 = x^2 + R(-1)*x*y + R(-1)*x + R(1/3)
            sage: p2.tropical_variety()._parameter_intervals()
            [(-oo, 0], [0, +oo), [-1, 4/3], (-oo, 0], [0, +oo)]
        """
        from sage.sets.real_set import RealSet

        intervals = []
        R = self._poly.parent().base().base_ring()
        for component in self._hypersurface:
            if len(component[1]) == 1:
                interval = RealSet(component[1][0])
            else:
                lower = component[1][0].left()
                upper = component[1][1].right()
                if lower == -infinity:
                    interval = RealSet(-infinity, infinity)
                else:
                    interval = RealSet([R(lower), R(upper)])
            intervals.append(interval)
        return intervals

    def plot(self):
        """
        Return the plot of ``self``.

        Generates a visual representation of the tropical curve in cartesian
        coordinates. The plot shows piecewise-linear segments representing
        each components. The axes are centered around the vertices.

        OUTPUT:

        A Graphics object. The weight of the component will be written if it
        is greater or equal than 2. The weight is written near the vertex.

        EXAMPLES:

        A polynomial with only two terms will give one straight line::

            sage: T = TropicalSemiring(QQ)
            sage: R.<x,y> = PolynomialRing(T)
            sage: (y+R(1)).tropical_variety().components()
            [[(t1, 1), [-Infinity < t1, t1 < +Infinity], 1]]
            sage: (y+R(1)).tropical_variety().plot()
            Graphics object consisting of 1 graphics primitive

        .. PLOT::
            :width: 300 px

            T = TropicalSemiring(QQ)
            R = PolynomialRing(T, ('x,y'))
            x, y = R.gen(), R.gen(1)
            sphinx_plot((y+R(1)).tropical_variety().plot())

        An intriguing and fascinating tropical curve can be obtained with
        a more complex tropical polynomial::

            sage: p1 = R(1) + R(2)*x + R(3)*y + R(6)*x*y + R(10)*x*y^2
            sage: p1.tropical_variety().components()
            [[(-1, t1), [-2 <= t1], 1],
            [(t1, -2), [-1 <= t1], 1],
            [(t1 + 1, t1), [-4 <= t1, t1 <= -2], 1],
            [(t1, -4), [t1 <= -3], 2],
            [(-t1 - 7, t1), [t1 <= -4], 1]]
            sage: p1.tropical_variety().plot()
            Graphics object consisting of 6 graphics primitives

        .. PLOT::
            :width: 300 px

            T = TropicalSemiring(QQ)
            R = PolynomialRing(T, ('x,y'))
            x, y = R.gen(), R.gen(1)
            p1 = R(1) + R(2)*x + R(3)*y + R(6)*x*y + R(10)*x*y**2
            sphinx_plot(p1.tropical_variety().plot())

        Another tropical polynomial with numerous components, resulting
        in a more intricate structure::

            sage: p2 = (R(8) + R(4)*x + R(2)*y + R(1)*x^2 + x*y + R(1)*y^2
            ....:       + R(2)*x^3 + x^2*y + x*y^2 + R(4)*y^3 + R(8)*x^4
            ....:       + R(4)*x^3*y + x^2*y^2 + R(2)*x*y^3 + y^4)
            sage: p2.tropical_variety().plot()
            Graphics object consisting of 23 graphics primitives

        .. PLOT::
            :width: 300 px

            T = TropicalSemiring(QQ)
            R = PolynomialRing(T, ('x,y'))
            x, y = R.gen(), R.gen(1)
            p2 = R(8) + R(4)*x + R(2)*y + R(1)*x**2 + x*y + R(1)*y**2 \
            + R(2)*x**3 + x**2*y + x*y**2 + R(4)*y**3 + R(8)*x**4 \
            + R(4)*x**3*y + x**2*y**2 + R(2)*x*y**3 + y**4
            sphinx_plot(p2.tropical_variety().plot())
        """
        from sage.plot.plot import plot
        from sage.plot.text import text
        from sage.plot.graphics import Graphics
        from sage.plot.plot import parametric_plot

        if not self._hypersurface:
            return plot(lambda x: float('nan'), {-1, 1})

        combined_plot = Graphics()
        large_int = 100
        intervals = self._parameter_intervals()
        for i, component in enumerate(self._hypersurface):
            var = component[1][0].variables()[0]
            parametric_function = component[0]
            order = component[2]
            interval = intervals[i]
            if interval[0].lower() == -infinity:
                lower = interval[0].upper() - large_int
                upper = interval[0].upper()
                midpoint = upper - 0.5
            elif interval[0].upper() == infinity:
                lower = interval[0].lower()
                upper = interval[0].lower() + large_int
                midpoint = lower + 0.5
            else:
                lower = interval[0].lower()
                upper = interval[0].upper()
                midpoint = (lower+upper) / 2

            if (lower == infinity) and (upper == infinity):
                midpoint = 0
                plot = parametric_plot(parametric_function, (var, -large_int,
                                       large_int), color='red')
            else:
                plot = parametric_plot(parametric_function, (var, lower, upper),
                                       color='red')

            # Add the order if it is greater than or equal to 2
            if component[2] > 1:
                point = []
                for eq in component[0]:
                    value = eq.subs(**{str(var): midpoint})
                    point.append(value)
                text_order = text(str(order), (point[0], point[1]),
                                  fontsize=16, color='black')
                combined_plot += plot + text_order
            else:
                combined_plot += plot

        # Set default axes
        axes = self._axes()
        xmin, xmax = axes[0][0], axes[0][1]
        ymin, ymax = axes[1][0], axes[1][1]
        combined_plot.set_axes_range(xmin=xmin, xmax=xmax,
                                     ymin=ymin, ymax=ymax)
        return combined_plot

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: T = TropicalSemiring(QQ)
            sage: R.<x,y> = PolynomialRing(T)
            sage: (x^2+R(0)).tropical_variety()
            Tropical curve of 0*x^2 + 0
        """
        return f"Tropical curve of {self._poly}"
