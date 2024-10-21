# sage_setup: distribution = sagemath-symbolics
r"""
Integrated Curves and Geodesics in Manifolds

Given a differentiable manifold `M`, an *integrated curve* in `M`
is a differentiable curve constructed as a solution to a system of
second order differential equations.

Integrated curves are implemented by the class :class:`IntegratedCurve`, from
which the classes :class:`IntegratedAutoparallelCurve` and
:class:`IntegratedGeodesic` inherit.

.. RUBRIC:: Example: a geodesic in the hyperbolic plane

First declare the hyperbolic plane as a 2-dimensional Riemannian manifold ``M``
and introduce the chart ``X`` corresponding to the Poincaré half-plane model::

    sage: M = Manifold(2, 'M', structure='Riemannian')
    sage: X.<x,y> = M.chart('x y:(0,+oo)')

Then set the metric to be the hyperbolic one::

    sage: g = M.metric()
    sage: g[0,0], g[1,1] = 1/y^2, 1/y^2
    sage: g.display()
    g = y^(-2) dx⊗dx + y^(-2) dy⊗dy

Pick an initial point and an initial tangent vector::

    sage: p = M((0,1), name='p')
    sage: v = M.tangent_space(p)((1,3/2), name='v')
    sage: v.display()
    v = ∂/∂x + 3/2 ∂/∂y

Declare a geodesic with such initial conditions, denoting by `t` the
corresponding affine parameter::

    sage: t = var('t')
    sage: c = M.integrated_geodesic(g, (t, 0, 10), v, name='c')

Numerically integrate the geodesic (see :meth:`~IntegratedCurve.solve` for
all possible options, including the choice of the numerical algorithm)::

    sage: sol = c.solve()                                                               # needs scipy

Plot the geodesic after interpolating the solution ``sol``::

    sage: interp = c.interpolate()

    sage: # needs sage.plot
    sage: graph = c.plot_integrated()
    sage: p_plot = p.plot(size=30, label_offset=-0.07, fontsize=20)
    sage: v_plot = v.plot(label_offset=0.05, fontsize=20)
    sage: graph + p_plot + v_plot
    Graphics object consisting of 5 graphics primitives

.. PLOT::

    M = Manifold(2, 'M', structure='Riemannian')
    X = M.chart('x y'); x, y = X[:]
    g = M.metric()
    g[0,0], g[1,1] = 1/y**2, 1/y**2
    p = M((0,1), name='p')
    v = M.tangent_space(p)((1,3/2), name='v')
    t = var('t')
    c = M.integrated_geodesic(g, (t, 0, 10), v, name='c')
    sol = c.solve()
    interp = c.interpolate()
    graph = c.plot_integrated()
    p_plot = p.plot(size=30, label_offset=-0.07, fontsize=20)
    v_plot = v.plot(label_offset=0.05, fontsize=20)
    sphinx_plot(graph + p_plot + v_plot)

`c` is a differentiable curve in `M` and inherits from the properties of
:class:`~sage.manifolds.differentiable.curve.DifferentiableCurve`::

    sage: c.domain()
    Real interval (0, 10)
    sage: c.codomain()
    2-dimensional Riemannian manifold M
    sage: c.display()
    c: (0, 10) → M

In particular, its value at `t=1` is::

    sage: c(1)
    Point on the 2-dimensional Riemannian manifold M

which corresponds to the following `(x, y)` coordinates::

    sage: X(c(1))  # abs tol 1e-12
    (2.4784140715580136, 1.5141683866138937)

AUTHORS:

- Karim Van Aelst (2017): initial version
- Florentin Jaffredo (2018): integration over multiple charts, use of
  ``fast_callable`` to improve the computation speed
"""

# **********************************************************************
#       Copyright (C) 2017 Karim Van Aelst <karim.van-aelst@obspm.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# **********************************************************************

from sage.symbolic.expression import Expression
from sage.rings.infinity import Infinity
from sage.calculus.desolvers import desolve_system_rk4
from sage.calculus.desolvers import desolve_odeint
from sage.manifolds.chart import Chart
from sage.manifolds.differentiable.curve import DifferentiableCurve
from sage.manifolds.differentiable.tangent_vector import TangentVector
from sage.calculus.interpolation import Spline
from sage.misc.decorators import options
from sage.misc.functional import numerical_approx
from sage.misc.lazy_import import lazy_import
from sage.arith.srange import srange
from sage.ext.fast_callable import fast_callable
from sage.symbolic.ring import SR
from random import shuffle

lazy_import('scipy.integrate', 'ode')


class IntegratedCurve(DifferentiableCurve):
    r"""
    Given a chart with coordinates denoted `(x_{1}, \ldots, x_{n})`,
    an instance of :class:`IntegratedCurve` is a curve
    `t \mapsto (x_{1}(t), \ldots, x_{n}(t))` constructed as a
    solution to a system of second order differential equations
    satisfied by the coordinate curves `t \mapsto x_{i}(t)`.

    INPUT:

    - ``parent`` --
      :class:`~sage.manifolds.differentiable.manifold_homset.IntegratedCurveSet`
      the set of curves `\mathrm{Hom_{integrated}}(I, M)` to which the
      curve belongs
    - ``equations_rhs`` -- list of the right-hand sides of the equations
      on the velocities only (the term *velocity* referring to the
      derivatives `d x_{i} / dt` of the coordinate curves)
    - ``velocities`` -- list of the symbolic expressions used in
      ``equations_rhs`` to denote the velocities
    - ``curve_parameter`` -- symbolic expression used in
      ``equations_rhs`` to denote the parameter of the curve (denoted
      `t` in the descriptions above)
    - ``initial_tangent_vector`` --
      :class:`~sage.manifolds.differentiable.tangent_vector.TangentVector`
      initial tangent vector of the curve
    - ``chart`` -- (default: ``None``) chart on the manifold in
      which the equations are given; if ``None`` the default chart
      of the manifold is assumed
    - ``name`` -- (default: ``None``) string; symbol given to the curve
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to
      denote the curve; if none is provided, ``name`` will be used

    EXAMPLES:

    Motion of a charged particle in an axial magnetic field linearly
    increasing in time and exponentially decreasing in space:

    .. MATH::

        \mathbf{B}(t,\mathbf{x}) = \frac{B_{0}t}{T} \exp \left(
        -\frac{ x_{1}^{2} + x_{2}^{2} }{ L^{2} } \right) \mathbf{e_{3}}.

    Equations of motion are:

    .. MATH::

        \begin{aligned}
        \ddot{x}_{1}(t) &= \frac{qB(t,\mathbf{x}(t))}{m} \dot{x}_{2}(t), \\
        \ddot{x}_{2}(t) &= -\frac{qB(t, \mathbf{x}(t))}{m} \dot{x}_{1}(t), \\
        \ddot{x}_{3}(t) &= 0.
        \end{aligned}

    Start with declaring a chart on a 3-dimensional manifold and the
    symbolic expressions denoting the velocities and the various
    parameters::

        sage: M = Manifold(3, 'M', start_index=1)
        sage: X.<x1,x2,x3> = M.chart()
        sage: var('t B_0 m q L T')
        (t, B_0, m, q, L, T)
        sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
        sage: D = X.symbolic_velocities(); D
        [Dx1, Dx2, Dx3]
        sage: eqns = [q*B/m*D[1], -q*B/m*D[0], 0]

    Set the initial conditions::

        sage: p = M.point((0,0,0), name='p')
        sage: Tp = M.tangent_space(p)
        sage: v = Tp((1,0,1))

    Declare an integrated curve and display information relative to it::

        sage: c = M.integrated_curve(eqns, D, (t, 0, 5), v, name='c',
        ....:                                              verbose=True)
        The curve was correctly set.
        Parameters appearing in the differential system defining the
         curve are [B_0, L, T, m, q].
        sage: c
        Integrated curve c in the 3-dimensional differentiable
         manifold M
        sage: sys = c.system(verbose=True)
        Curve c in the 3-dimensional differentiable manifold M
         integrated over the Real interval (0, 5) as a solution to the
         following system, written with respect to
         Chart (M, (x1, x2, x3)):
        <BLANKLINE>
        Initial point: Point p on the 3-dimensional differentiable
         manifold M with coordinates [0, 0, 0] with respect to
         Chart (M, (x1, x2, x3))
        Initial tangent vector: Tangent vector at Point p on
         the 3-dimensional differentiable manifold M with
         components [1, 0, 1] with respect to Chart (M, (x1, x2, x3))
        <BLANKLINE>
        d(x1)/dt = Dx1
        d(x2)/dt = Dx2
        d(x3)/dt = Dx3
        d(Dx1)/dt = B_0*Dx2*q*t*e^(-(x1^2 + x2^2)/L^2)/(T*m)
        d(Dx2)/dt = -B_0*Dx1*q*t*e^(-(x1^2 + x2^2)/L^2)/(T*m)
        d(Dx3)/dt = 0
        <BLANKLINE>

    Generate a solution of the system and an interpolation of this
    solution::

        sage: sol = c.solve(step=0.2,                                                   # needs scipy
        ....:         parameters_values={B_0:1, m:1, q:1, L:10, T:1},
        ....:         solution_key='carac time 1', verbose=True)
        Performing numerical integration with method 'odeint'...
        Numerical integration completed.
        <BLANKLINE>
        Checking all points are in the chart domain...
        All points are in the chart domain.
        <BLANKLINE>
        The resulting list of points was associated with the key
         'carac time 1' (if this key already referred to a former
         numerical solution, such a solution was erased).
        sage: interp = c.interpolate(solution_key='carac time 1',                       # needs scipy
        ....:                interpolation_key='interp 1', verbose=True)
        Performing cubic spline interpolation by default...
        Interpolation completed and associated with the key 'interp 1'
         (if this key already referred to a former interpolation,
         such an interpolation was erased).

    Such an interpolation is required to evaluate the curve and the
    vector tangent to the curve for any value of the curve parameter::

        sage: # needs scipy
        sage: p = c(1.9, verbose=True)
        Evaluating point coordinates from the interpolation associated
         with the key 'interp 1' by default...
        sage: p
        Point on the 3-dimensional differentiable manifold M
        sage: p.coordinates()     # abs tol 1e-12
        (1.377689074756845, -0.900114533011232, 1.9)
        sage: v2 = c.tangent_vector_eval_at(4.3, verbose=True)
        Evaluating tangent vector components from the interpolation
         associated with the key 'interp 1' by default...
        sage: v2
        Tangent vector at Point on the 3-dimensional differentiable
         manifold M
        sage: v2[:]     # abs tol 1e-12
        [-0.9425156073651124, -0.33724314284285434, 1.0]

    Plotting a numerical solution (with or without its tangent vector
    field) also requires the solution to be interpolated at least once::

        sage: c_plot_2d_1 = c.plot_integrated(ambient_coords=[x1, x2],                  # needs scipy
        ....:               interpolation_key='interp 1', thickness=2.5,
        ....:               display_tangent=True, plot_points=200,
        ....:               plot_points_tangent=10, scale=0.5,
        ....:               color='blue', color_tangent='red',
        ....:               verbose=True)
        A tiny final offset equal to 0.000251256281407035 was introduced
         for the last point in order to safely compute it from the
         interpolation.
        sage: c_plot_2d_1                                                               # needs scipy sage.plot
        Graphics object consisting of 11 graphics primitives

    .. PLOT::

        M = Manifold(3, 'M')
        X = M.chart('x1 x2 x3'); x1, x2, x3 = X[:]
        t, B_0, m, q, L, T = var('t B_0 m q L T')
        B = B_0*t/T*exp(-(x1**2 + x2**2)/L**2)
        D = X.symbolic_velocities()
        eqns = [q*B/m*D[1], -q*B/m*D[0], 0]
        p = M.point((0,0,0), name='p')
        Tp = M.tangent_space(p)
        v = Tp((1,0,1))
        c = M.integrated_curve(eqns, D, (t, 0, 5), v, name='c')
        sol = c.solve(step=0.2,
                      parameters_values={B_0:1, m:1, q:1, L:10, T:1},
                      solution_key='carac time 1')
        interp = c.interpolate(solution_key='carac time 1',
                               interpolation_key='interp 1')
        c_plot_2d_1 = c.plot_integrated(ambient_coords=[x1, x2],
                        interpolation_key='interp 1', thickness=2.5,
                        display_tangent=True, plot_points=200,
                        plot_points_tangent=10, scale=0.5, color='blue',
                        color_tangent='red')
        sphinx_plot(c_plot_2d_1)

    An instance of :class:`IntegratedCurve` may store several numerical
    solutions and interpolations::

        sage: # needs scipy
        sage: sol = c.solve(step=0.2,
        ....:         parameters_values={B_0:1, m:1, q:1, L:10, T:100},
        ....:         solution_key='carac time 100')
        sage: interp = c.interpolate(solution_key='carac time 100',
        ....:                            interpolation_key='interp 100')
        sage: c_plot_3d_100 = c.plot_integrated(interpolation_key='interp 100',         # needs sage.plot
        ....:                   thickness=2.5, display_tangent=True,
        ....:                   plot_points=200, plot_points_tangent=10,
        ....:                   scale=0.5, color='green',
        ....:                   color_tangent='orange')
        sage: c_plot_3d_1 = c.plot_integrated(interpolation_key='interp 1',             # needs sage.plot
        ....:                   thickness=2.5, display_tangent=True,
        ....:                   plot_points=200, plot_points_tangent=10,
        ....:                   scale=0.5, color='blue',
        ....:                   color_tangent='red')
        sage: c_plot_3d_1 + c_plot_3d_100                                               # needs sage.plot
        Graphics3d Object

    .. PLOT::

        M = Manifold(3, 'M')
        X = M.chart('x1 x2 x3'); x1, x2, x3 = X[:]
        t, B_0, m, q, L, T = var('t B_0 m q L T')
        B = B_0*t/T*exp(-(x1**2 + x2**2)/L**2)
        D = X.symbolic_velocities()
        eqns = [q*B/m*D[1], -q*B/m*D[0], 0]
        p = M.point((0,0,0), name='p')
        Tp = M.tangent_space(p)
        v = Tp((1,0,1))
        c = M.integrated_curve(eqns, D, (t, 0, 5), v, name='c')
        sol = c.solve(step=0.2, parameters_values={B_0:1, m:1, q:1, L:10, T:1},
                      solution_key='carac time 1')
        interp = c.interpolate(solution_key='carac time 1',
                               interpolation_key='interp 1')
        sol = c.solve(step=0.2, parameters_values={B_0:1, m:1, q:1, L:10, T:100},
                      solution_key='carac time 100')
        interp = c.interpolate(solution_key='carac time 100',
                               interpolation_key='interp 100')
        c_plot_3d_1 = c.plot_integrated(interpolation_key='interp 1',
                        thickness=2.5, display_tangent=True,
                        plot_points=200, plot_points_tangent=10,
                        scale=0.5, color='blue', color_tangent='red')
        c_plot_3d_100 = c.plot_integrated(interpolation_key='interp 100',
                            thickness=2.5, display_tangent=True,
                            plot_points=200, plot_points_tangent=10,
                            scale=0.5, color='green',
                            color_tangent='orange')
        graph = c_plot_3d_1 + c_plot_3d_100
        sphinx_plot(graph)
    """

    def __init__(self, parent, equations_rhs, velocities,
                 curve_parameter, initial_tangent_vector, chart=None,
                 name=None, latex_name=None, verbose=False,
                 across_charts=False):
        r"""
        Construct a curve defined by a system of second order
        differential equations in the coordinate functions.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: [t, B_0, m, q, L, T] = var('t B_0 m q L T')
            sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B/m*D[1], -q*B/m*D[0], 0]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns + [x1], D, (t, 0, 5), v,
            ....:                                              name='c')
            Traceback (most recent call last):
            ...
            ValueError: number of equations should equal codomain
             dimension
            sage: c = M.integrated_curve(eqns, D + [x1], (t, 0, 5), v,
            ....:                                              name='c')
            Traceback (most recent call last):
            ...
            ValueError: number of velocities should equal codomain
             dimension
            sage: c = M.integrated_curve(eqns, D,(t,-oo,5), v, name='c')
            Traceback (most recent call last):
            ...
            ValueError: both boundaries of the interval defining the
             domain of a Homset of integrated curves need to be finite
            sage: c = M.integrated_curve(eqns, D, (t,0,5), x1, name='c')
            Traceback (most recent call last):
            ...
            TypeError: x1 should be a tangent vector

            sage: c = M.integrated_curve(eqns, D, (x1,0,5), v, name='c')
            Traceback (most recent call last):
            ...
            ValueError: x1 should not be used as the curve parameter
             since it also denotes a coordinate or a velocity
            sage: c = M.integrated_curve(eqns, D, (t,0,5), v, name='c'); c
            Integrated curve c in the 3-dimensional differentiable
             manifold M
            sage: TestSuite(c).run()

        Check that :issue:`28669` is fixed::

            sage: E.<r,phi> = EuclideanSpace(coordinates='polar')
            sage: p = E((1, 0))  # the initial point
            sage: v = E.tangent_space(p)((2, 1))  # the initial vector
            sage: t = var('t')
            sage: c = E.integrated_geodesic(E.metric(), (t, 0, 10), v); c
            Integrated geodesic in the Euclidean plane E^2
        """
        from sage.symbolic.ring import SR

        # start with parent class method to initialize the four last
        # arguments:
        DifferentiableCurve.__init__(self, parent, name=name,
                                     latex_name=latex_name)

        # check argument 'parent': 't_min' and 't_max' below are only
        # allowed to be either expressions of finite real values:
        domain = self.domain()
        t_min = domain.lower_bound()
        t_max = domain.upper_bound()
        if t_min == -Infinity or t_max == +Infinity:
            raise ValueError("both boundaries of the interval " +
                             "need to be finite")

        codomain = self.codomain()

        # check argument 'equations_rhs':
        dim = codomain.dim()

        if not isinstance(equations_rhs, dict):
            if len(equations_rhs) != dim:
                raise ValueError("number of equations should equal " +
                                 "codomain dimension")
        else:
            for eq in equations_rhs.values():
                if len(eq) != dim:
                    raise ValueError("number of equations should equal " +
                                     "codomain dimension")

        # check the chart:
        if chart is not None:
            if chart not in codomain.atlas():
                raise ValueError("{} should be a chart ".format(chart) +
                                 "on the {}".format(codomain))
        else:
            chart = codomain.default_chart()

        # check argument 'velocities':
        if len(velocities) != dim:
            raise ValueError("number of velocities should equal " +
                             "codomain dimension")
        # in particular, check that no velocity coincides with a
        # coordinate:
        for vel in velocities:
            if vel in chart[:]:
                str_error = "{} should not be used as a ".format(vel)
                str_error += "velocity since it also denotes "
                str_error += "a coordinate"
                raise ValueError(str_error)

        # check argument 'curve_parameter':
        if not isinstance(curve_parameter, Expression):
            raise TypeError("{} should be ".format(curve_parameter) +
                             "a symbolic expression")
        # in particular, check that it does not coincide with a
        # coordinate or a velocity:
        coords_vels = list(chart[:]) + list(velocities)
        if curve_parameter in coords_vels:
            str_error = "{} should not be used ".format(curve_parameter)
            str_error += "as the curve parameter since it also denotes "
            str_error += "a coordinate or a velocity"
            raise ValueError(str_error)
        # the various algorithms called in 'solve' method are in charge
        # of raising errors about possibly remaining problems regarding
        # 'curve_parameter'

        # check argument 'initial_tangent_vector':
        if not isinstance(initial_tangent_vector, TangentVector):
            raise TypeError("{} ".format(initial_tangent_vector) +
                            "should be a tangent vector")
        initial_pt = initial_tangent_vector.parent().base_point()
        # line above retrieves the initial point as the base point of
        # the tangent space to which the initial tangent vector belongs
        initial_pt_coords = initial_pt.coordinates(chart)
        # prepare attribute '_parameters':
        announced_variables = set(coords_vels + [curve_parameter])
        parameters = set()
        # extract all the variables appearing in the equations:
        for eqn in equations_rhs:
            if isinstance(eqn, Expression): # some right hand sides
            # might merely be real numbers and not expressions, so that
            # they do not contain any variable, and method 'variables'
            # could not be called on them
                parameters = parameters.union(eqn.variables())
        # remove the Expressions that should not be treated as
        # parameters (i.e. the coordinate functions, the velocities and
        # the curve parameter):
        parameters = parameters.difference(announced_variables)
        # extract all the variables appearing in the boundaries:
        if isinstance(t_min, Expression):
            parameters = parameters.union(t_min.variables())
        if isinstance(t_max, Expression):
            parameters = parameters.union(t_max.variables())
        # extract all the variables appearing in the initial point
        # coordinates:
        for coord in initial_pt_coords:
            if isinstance(coord,Expression):
                parameters = parameters.union(coord.variables())
        # extract all the variables appearing in the initial tangent
        # vector components:
        initial_coord_basis = chart.frame().at(initial_pt)
        initial_tgt_vec_comps = initial_tangent_vector[initial_coord_basis,:]
        for comp in initial_tgt_vec_comps:
            if isinstance(comp, Expression):
                parameters = parameters.union(comp.variables())

        # check at this stage that no parameter coincides with a
        # coordinate, a velocity, or the curve parameter; this would
        # mean that an Expression used to denote either a bound, a
        # coordinate of the initial point or a component of the initial
        # tangent vector coincides with a coordinate, a velocity or the
        # curve parameter (which would make no sense):
        if len(parameters) != 0:
            for param in parameters:
                if param in announced_variables:
                    str_error = "{} should not be used ".format(param)
                    str_error += "as a parameter since it also denotes "
                    str_error += "a coordinate, a velocity or the "
                    str_error += "curve parameter"
                    raise ValueError(str_error)

        # define all attributes
        if not isinstance(equations_rhs, dict):
            self._equations_rhs = list(equations_rhs) # converts to list
            # since might not already be a list (which is later required)
        else: # case multi charts
            self._equations_rhs = equations_rhs

        self._across_charts = across_charts
        if across_charts:
            # pre-compute the changes of chart for faster switching
            # approx gain : 200 ms per switch
            self._fast_changes_of_frame = {}
            self._fast_changes_of_chart = {}
            for CoF in self._codomain.changes_of_frame():
                M = self._codomain.changes_of_frame()[CoF][CoF[1], :, CoF[1]._chart]
                M = M.apply_map(lambda e: e.expr())
                M = M.numpy()
                for i in range(dim):
                    for j in range(dim):
                        M[i,j] = fast_callable(SR(M[i, j]), vars=list(CoF[1]._chart[:]), domain=float)

                import numpy as np
                def fast_CoF(pos, vel, M=M):
                # using default arguments for binding (ugly python)
                    #print(det(*pos))
                    return list(np.dot( [[M[j, i](*pos) for i in range(dim)]
                                    for j in range(dim)], vel))

                self._fast_changes_of_frame[CoF] = fast_CoF

            for CoC in self._codomain._coord_changes:
                transf = self._codomain._coord_changes[CoC]._transf
                fast_transf = [fast_callable(f.expr(), vars=list(CoC[0][:]), domain=float)
                               for f in transf]
                self._fast_changes_of_chart[CoC] = fast_transf

        self._velocities = list(velocities) # converts to list
        # since might not already be a list (which is later required)
        self._curve_parameter = curve_parameter
        self._initial_tangent_vector = initial_tangent_vector
        self._chart = chart
        self._parameters = parameters
        self._ode_solver = None # if needed, becomes an instance of
        # 'ode_solver', which performs most of the numerical integrations
        # offered by method 'solve'
        self._solutions = {} # dictionary containing all numerically
        # computed lists of points of the curve, the keys being chosen
        # by the user when calling method 'solve'
        self._interpolations = {} # dictionary containing lists of
        # interpolation objects, each interpolation object implementing
        # the interpolation of one of the numerical coordinate curves,
        # and the keys being chosen by the user when calling
        # method 'interpolate'

        if verbose:
            print("The curve was correctly set.")
            if self._parameters:
                print("Parameters appearing in the differential " +
                      "system defining the curve are " +
                      "{}.".format(sorted(self._parameters, key=str)))
            else:
                print("No parameter appears in the differential " +
                      "system defining the curve.")

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: [t, B_0, m, q, L, T] = var('t B_0 m q L T')
            sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B/m*D[1], -q*B/m*D[0], 0]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t, 0, 5), v) ; c
            Integrated curve in the 3-dimensional differentiable
             manifold M
            sage: c = M.integrated_curve(eqns, D, (t,0,5), v, name='c'); c
            Integrated curve c in the 3-dimensional differentiable
             manifold M
        """

        description = "Integrated curve "
        if self._name is not None:
            description += self._name + " "
        description += "in the {}".format(self._codomain)
        return description

    def __reduce__(self):
        r"""
        Reduction function for the pickle protocole.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: [t, B_0, m, q, L, T] = var('t B_0 m q L T')
            sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B/m*D[1], -q*B/m*D[0], 0]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t,0,5), v, name='c')
            sage: c.__reduce__()
            (<class 'sage.manifolds.differentiable.manifold_homset.IntegratedCurveSet_with_category.element_class'>,
             (Set of Morphisms from Real interval (0, 5) to
              3-dimensional differentiable manifold M in Category of homsets of
              topological spaces which actually are integrated curves,
              [B_0*Dx2*q*t*e^(-(x1^2 + x2^2)/L^2)/(T*m),
               -B_0*Dx1*q*t*e^(-(x1^2 + x2^2)/L^2)/(T*m),
               0],
              [Dx1, Dx2, Dx3],
              t,
              Tangent vector at Point p on the 3-dimensional
               differentiable manifold M,
              Chart (M, (x1, x2, x3)),
              'c',
              'c',
              False,
              False))

        Test of pickling::

            sage: loads(dumps(c))
            Integrated curve c in the 3-dimensional differentiable
             manifold M
        """
        return (type(self), (self.parent(), self._equations_rhs,
                self._velocities, self._curve_parameter,
                self._initial_tangent_vector, self._chart,
                self._name, self._latex_name, False, self._across_charts))

    def system(self, verbose=False):
        r"""
        Provide a detailed description of the system defining the curve
        and return the system defining it: chart, equations and initial
        conditions.

        INPUT:

        - ``verbose`` -- boolean (default: ``False``); prints a detailed
          description of the curve

        OUTPUT: list containing

          * the equations
          * the initial conditions
          * the chart

        EXAMPLES:

        System defining an integrated curve::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: [t, B_0, m, q, L, T] = var('t B_0 m q L T')
            sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B/m*D[1], -q*B/m*D[0], 0]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t,0,5), v, name='c')
            sage: sys = c.system(verbose=True)
            Curve c in the 3-dimensional differentiable manifold M
             integrated over the Real interval (0, 5) as a solution to
             the following system, written with respect to
             Chart (M, (x1, x2, x3)):
            <BLANKLINE>
            Initial point: Point p on the 3-dimensional differentiable
             manifold M with coordinates [0, 0, 0] with respect to
             Chart (M, (x1, x2, x3))
            Initial tangent vector: Tangent vector at Point p on the
             3-dimensional differentiable manifold M with
             components [1, 0, 1] with respect to Chart (M, (x1, x2, x3))
            <BLANKLINE>
            d(x1)/dt = Dx1
            d(x2)/dt = Dx2
            d(x3)/dt = Dx3
            d(Dx1)/dt = B_0*Dx2*q*t*e^(-(x1^2 + x2^2)/L^2)/(T*m)
            d(Dx2)/dt = -B_0*Dx1*q*t*e^(-(x1^2 + x2^2)/L^2)/(T*m)
            d(Dx3)/dt = 0
            <BLANKLINE>
            sage: sys_mute = c.system()
            sage: sys_mute == sys
            True
        """

        v0 = self._initial_tangent_vector
        chart = self._chart

        if verbose:
            initial_tgt_space = v0.parent()
            initial_pt = initial_tgt_space.base_point() # retrieves
            # the initial point as the base point of the tangent space
            # to which initial tangent vector belongs
            initial_pt_coords = list(initial_pt.coordinates(chart))
            # previous line converts to list since would otherwise be a
            # tuple ; will raise error if coordinates in chart are not
            # known

            initial_coord_basis = chart.frame().at(initial_pt)
            initial_tgt_vec_comps = v0[initial_coord_basis,:] # will
            # raise error if components in coordinate basis are not
            # known

            description = "Curve "
            if self._name is not None:
                description += self._name + " "
            description += "in the {} ".format(self.codomain())
            description += "integrated over the "
            description += "{} ".format(self.domain())
            description += "as a solution to the following system, "
            description += "written with respect to "
            description += "{}:\n\n".format(chart)

            description += "Initial point: {} ".format(initial_pt)
            description += "with coordinates "
            description += "{} ".format(initial_pt_coords)
            description += "with respect to {}\n".format(chart)

            description += "Initial tangent vector: {} ".format(v0)
            description += "with components "
            description += "{}".format(initial_tgt_vec_comps)
            description += " with respect to {}\n\n".format(chart)

            for coord_func,velocity in zip(chart[:],self._velocities):
                description += "d({})/d{} = {}\n".format(coord_func,
                                                  self._curve_parameter,
                                                  velocity)

            for velocity,eqn in zip(self._velocities,self._equations_rhs):
                description += "d({})/d{} = {}\n".format(velocity,
                                                  self._curve_parameter,
                                                  eqn)

            print(description)

        return [self._equations_rhs, v0, chart]

    def solve_analytical(self, verbose=False):
        r"""
        Solve the differential system defining ``self`` analytically.

        Solve analytically the differential system defining a curve
        using Maxima via Sage solver ``desolve_system``.
        In case of success, the analytical expressions are added to the
        dictionary of expressions representing the curve.
        Pay attention to the fact that ``desolve_system`` only considers
        initial conditions given at an initial parameter value equal to
        zero, although the parameter range may not contain zero.
        Yet, assuming that it does, values of the coordinates functions
        at such zero initial parameter value are denoted by the name of
        the coordinate function followed by the string ``'_0'``.

        OUTPUT:

        - list of the analytical expressions of the coordinate functions
          (when the differential system could be solved analytically),
          or boolean ``False`` (in case the differential system could
          not be solved analytically)

        EXAMPLES:

        Analytical expression of the trajectory of a charged particle in
        a uniform, stationary magnetic field::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: [t, B_0, m, q] = var('t B_0 m q')
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B_0/m*D[1], -q*B_0/m*D[0], 0]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t,0,5), v, name='c')
            sage: sys = c.system(verbose=True)
            Curve c in the 3-dimensional differentiable manifold M
             integrated over the Real interval (0, 5) as a solution to
             the following system, written with respect to
             Chart (M, (x1, x2, x3)):
            <BLANKLINE>
            Initial point: Point p on the 3-dimensional differentiable
             manifold M with coordinates [0, 0, 0] with respect to
             Chart (M, (x1, x2, x3))
            Initial tangent vector: Tangent vector at Point p on the
             3-dimensional differentiable manifold M with components
             [1, 0, 1] with respect to Chart (M, (x1, x2, x3))
            <BLANKLINE>
            d(x1)/dt = Dx1
            d(x2)/dt = Dx2
            d(x3)/dt = Dx3
            d(Dx1)/dt = B_0*Dx2*q/m
            d(Dx2)/dt = -B_0*Dx1*q/m
            d(Dx3)/dt = 0
            <BLANKLINE>
            sage: sol = c.solve_analytical()
            sage: c.expr()
            ((B_0*q*x1_0 - Dx2_0*m*cos(B_0*q*t/m) +
               Dx1_0*m*sin(B_0*q*t/m) + Dx2_0*m)/(B_0*q),
             (B_0*q*x2_0 + Dx1_0*m*cos(B_0*q*t/m) +
              Dx2_0*m*sin(B_0*q*t/m) - Dx1_0*m)/(B_0*q),
             Dx3_0*t + x3_0)
        """

        from sage.calculus.var import function
        from sage.calculus.functional import diff
        from sage.calculus.desolvers import desolve_system
        from sage.symbolic.assumptions import assume, forget
        from sage.symbolic.ring import var

        dim = self.codomain().dim()
        i0 = self.codomain().start_index()
        des = self._velocities + self._equations_rhs
        par = self._curve_parameter

        for param in self._parameters:
            assume(param != 0)

        y = []
        for i in range(2*dim):
            name = "y{}".format(i+i0)
            y += [function(name)(par)]

        for i in range(dim):
            vel = self._velocities[i]
            des[i] = des[i].substitute({vel: y[dim+i]})
            des[i] = diff(y[i],par) == des[i]
            for j in range(dim):
                coord = self._chart[:][j] # important to use '[:]' on
                # 'chart' to avoid problems due to nonzero starting
                # index (i0)
                veloc = self._velocities[j]
                des[dim+i] = des[dim+i].substitute({coord: y[j]})
                des[dim+i] = des[dim+i].substitute({veloc: y[dim+j]})
            des[dim+i] = (diff(y[dim+i], par) == des[dim+i])

        dvars = y
        ics = [0]
        y_ics_first_half = []
        y_ics_second_half = []
        for i in range(dim):
            coord = self._chart[:][i] # important to use '[:]'
            # on 'chart' to avoid problems due to nonzero
            # starting index (i0)
            veloc = self._velocities[i]
            str_var_coord = "{}_0".format(coord)
            str_var_veloc = "{}_0".format(veloc)
            y_coord_0 = var(str_var_coord)
            y_veloc_0 = var(str_var_veloc)
            y_ics_first_half += [y_coord_0]
            y_ics_second_half += [y_veloc_0]
        ics += y_ics_first_half + y_ics_second_half

        try:
            sol = desolve_system(des, dvars, ivar=self._curve_parameter, ics=ics)
        except NotImplementedError:
            coords_sol_expr = False
            if verbose:
                print("The system could not be solved analytically.")
        else:
            coords_sol_expr = []
            for relation in sol[:dim]:
                expr = relation.rhs().simplify_full()
                coords_sol_expr += [expr]
            self.add_expr(self.domain().default_chart(), self._chart,
                                                        coords_sol_expr)

        for param in self._parameters:
            forget(param != 0)

        return tuple(coords_sol_expr)

    def solve(self, step=None, method='odeint', solution_key=None,
              parameters_values=None, verbose=False, **control_param):
        r"""
        Integrate the curve numerically over the domain of definition.

        INPUT:

        - ``step`` -- (default: ``None``) step of integration; default
          value is a hundredth of the domain of integration if none is
          provided
        - ``method`` -- (default: ``'odeint'``) numerical scheme to
          use for the integration of the curve; available algorithms are:

          * ``'odeint'`` -- makes use of
            :func:`scipy:scipy.integrate.odeint`
            via Sage solver
            :func:`~sage.calculus.desolvers.desolve_odeint`; ``odeint`` invokes
            the LSODA algorithm of the
            `ODEPACK suite <https://www.netlib.org/odepack/>`_, which
            automatically selects between implicit Adams method (for non-stiff
            problems) and a method based on backward differentiation formulas
            (BDF) (for stiff problems).
          * ``'rk4_maxima'`` -- 4th order classical Runge-Kutta, which
            makes use of Maxima's dynamics package via Sage solver
            :func:`~sage.calculus.desolvers.desolve_system_rk4` (quite slow)
          * ``'dopri5'`` -- Dormand-Prince Runge-Kutta of order (4)5 provided by
            :obj:`scipy:scipy.integrate.ode`
          * ``'dop853'`` -- Dormand-Prince Runge-Kutta of order 8(5,3) provided by
            :obj:`scipy:scipy.integrate.ode`

          and those provided by ``GSL`` via Sage class
          :class:`~sage.calculus.ode.ode_solver`:

          * ``'rk2'`` -- embedded Runge-Kutta (2,3)
          * ``'rk4'`` -- 4th order classical Runge-Kutta
          * ``'rkf45'`` -- Runge-Kutta-Felhberg (4,5)
          * ``'rkck'`` -- embedded Runge-Kutta-Cash-Karp (4,5)
          * ``'rk8pd'`` -- Runge-Kutta Prince-Dormand (8,9)
          * ``'rk2imp'`` -- implicit 2nd order Runge-Kutta at Gaussian points
          * ``'rk4imp'`` -- implicit 4th order Runge-Kutta at Gaussian points
          * ``'gear1'`` -- `M=1` implicit Gear
          * ``'gear2'`` -- `M=2` implicit Gear
          * ``'bsimp'`` -- implicit Bulirsch-Stoer (requires Jacobian)

        - ``solution_key`` -- (default: ``None``) key which the
          resulting numerical solution will be associated to; a default
          value is given if none is provided
        - ``parameters_values`` -- (default: ``None``) list of numerical
          values of the parameters present in the system defining the
          curve, to be substituted in the equations before integration
        - ``verbose`` -- boolean (default: ``False``); prints information about
          the computation in progress
        - ``**control_param`` -- extra control parameters to be passed to the
          chosen solver; see the example with ``rtol`` and ``atol`` below

        OUTPUT: list of the numerical points of the computed solution

        EXAMPLES:

        Computing a numerical solution::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: [t, B_0, m, q, L, T] = var('t B_0 m q L T')
            sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B/m*D[1], -q*B/m*D[0], 0]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t,0,5), v, name='c')
            sage: sol = c.solve(parameters_values={B_0:1, m:1, q:1, L:10, T:1},         # needs scipy
            ....:               verbose=True)
            Performing numerical integration with method 'odeint'...
            Resulting list of points will be associated with the key
             'odeint' by default.
            Numerical integration completed.
            <BLANKLINE>
            Checking all points are in the chart domain...
            All points are in the chart domain.
            <BLANKLINE>
            The resulting list of points was associated with the key
             'odeint' (if this key already referred to a former
             numerical solution, such a solution was erased).

        The first 3 points of the solution, in the form ``[t, x1, x2, x3]``::

            sage: sol[:3]  # abs tol 1e-12                                              # needs scipy
            [[0.0, 0.0, 0.0, 0.0],
             [0.05, 0.04999999218759271, -2.083327338392213e-05, 0.05],
             [0.1, 0.09999975001847655, -0.00016666146190783666, 0.1]]

        The default is ``verbose=False``::

            sage: sol_mute = c.solve(parameters_values={B_0:1, m:1, q:1,                # needs scipy
            ....:                                       L:10, T:1})
            sage: sol_mute == sol                                                       # needs scipy
            True

        Specifying the relative and absolute error tolerance parameters to
        be used in :func:`~sage.calculus.desolvers.desolve_odeint`::

            sage: sol = c.solve(parameters_values={B_0:1, m:1, q:1, L:10, T:1},         # needs scipy
            ....:               rtol=1e-12, atol=1e-12)

        Using a numerical method different from the default one::

            sage: sol = c.solve(parameters_values={B_0:1, m:1, q:1, L:10, T:1},         # needs scipy
            ....:               method='rk8pd')


        TESTS::

            sage: sol = c.solve(parameters_values={m:1, q:1, L:10, T:1})                # needs scipy
            Traceback (most recent call last):
            ...
            ValueError: numerical values should be provided for each of
             the parameters [B_0, L, T, m, q]
            sage: sol = c.solve(method='my method',                                     # needs scipy
            ....:        parameters_values={B_0:1, m:1, q:1, L:10, T:1})
            Traceback (most recent call last):
            ...
            ValueError: no available method of integration referred to
             as 'my method'
        """
        from sage.symbolic.ring import SR

        if verbose:
            print("Performing numerical integration with method '" +
                  method + "'...")

        if solution_key is None:
            solution_key = method
            if verbose:
                print("Resulting list of points will be associated " +
                      "with the key '{}' ".format(solution_key) +
                      "by default.")

        t_min = self.domain().lower_bound()
        t_max = self.domain().upper_bound()

        eqns_num = list(self._equations_rhs)
        # 'self._equations_rhs' needs not to be modified ever, because we
        # want to keep track of the most general form of the equations
        # defining self, since those may contain parameters (which, for
        # instance, we want to display as their original expressions
        # when calling 'system' method with option 'verbose', and not
        # substituted with some numerical values).
        # This is why 'eqns_num' is declared: it will contain copies of
        # the equations of 'self._equations_rhs' in which the parameters
        # will be substituted with numerical values.
        # It was then important to declare it as above, in order to make
        # independent copies of each equations of 'self._equations_rhs',
        # rather than declaring 'eqns_num = self._equations_rhs', in which
        # case making substitutions in 'eqns_num' would have meant making
        # the same substitutions in the original equations of
        # 'self._equations_rhs'

        v0 = self._initial_tangent_vector
        chart = self._chart

        initial_tgt_space = v0.parent()
        initial_pt = initial_tgt_space.base_point() # retrieves
        # the initial point as the base point of the tangent space
        # to which the initial tangent vector belongs
        initial_pt_coords = list(initial_pt.coordinates(chart))
        # previous line converts to list since would otherwise be a
        # tuple (yet might need to be added to [t_min] later); will
        # raise error if coordinates in chart cannot be obtained

        initial_coord_basis = chart.frame().at(initial_pt)
        initial_tgt_vec_comps = list(v0[initial_coord_basis,:]) #idem

        dim = self.codomain().dim()

        if self._parameters:
            if parameters_values is None or len(parameters_values) != len(self._parameters):
                raise ValueError("numerical values should be " +
                                 "provided for each of the " +
                                 "parameters "
                                 "{}".format(sorted(self._parameters, key=str)))
            for key in parameters_values:
                # Get numerical values in case some parameters values
                # contain expressions such as pi; will raise error if
                # any element of parameters_values is not numerical
                parameters_values[key] = numerical_approx(parameters_values[key])

            if isinstance(t_min, Expression):
                t_min = parameters_values[t_min]
                if t_min == -Infinity or t_min == +Infinity:
                    raise ValueError("both boundaries of the " +
                                      "interval need to be finite")

            if isinstance(t_max, Expression):
                t_max = parameters_values[t_max]
                if t_max == -Infinity or t_max == +Infinity:
                    raise ValueError("both boundaries of the " +
                                     "interval need to be finite")

            for i in range(dim):
                if isinstance(eqns_num[i], Expression): # some right
                # hand sides might merely be real numbers and not
                # expressions, so that they do not contain any variable,
                # and hence no substitution is required
                    eqns_num[i] = eqns_num[i].substitute(parameters_values)

            for i in range(dim):
                if isinstance(initial_pt_coords[i], Expression):
                    AUX = initial_pt_coords[i]
                    AUX = AUX.substitute(parameters_values)
                    initial_pt_coords[i] = AUX
                if isinstance(initial_tgt_vec_comps[i], Expression):
                    AUX2 = initial_tgt_vec_comps[i]
                    AUX2 = AUX2.substitute(parameters_values)
                    initial_tgt_vec_comps[i] = AUX2
                # 'AUX' and 'AUX2' only used for the lines of
                # source code to be shorter

        t_min = numerical_approx(t_min)
        t_max = numerical_approx(t_max)

        for i in range(dim):
            if not isinstance(eqns_num[i], Expression): # in case of a
            # right hand side that is not an Expression (and then is a
            # number), it is needed to be converted to an Expression
            # since some solvers called below require only expressions
                eqns_num[i] = SR(eqns_num[i])

        if step is None:
            step = (t_max - t_min) / 100

        step = numerical_approx(step)

        initial_pt_coords = [numerical_approx(coord) for coord
                             in initial_pt_coords]
        initial_tgt_vec_comps = [numerical_approx(comp) for comp
                                 in initial_tgt_vec_comps]
        # the last two instructions retrieve numerical values even
        # if no parameters had to be substituted, in case some
        # coordinates or components contain expressions such as pi,
        # or are not RealNumber, since variable 'ics' of
        # 'desolve_system_rk4' used below needs to be a list of
        # RealNumber

        if not chart.valid_coordinates(*initial_pt_coords):
            raise ValueError("initial point should be in the " +
                             "domain of the chart")

        ode_solver_methods = ["rk2", "rk4", "rkf45", "rkck", "rk8pd",
                              "rk2imp", "rk4imp", "gear1", "gear2", "bsimp"]

        if method == 'rk4_maxima':
            des = self._velocities + eqns_num
            dvars = list(chart[:]) + self._velocities
            ics = [t_min] + initial_pt_coords + initial_tgt_vec_comps

            sol = desolve_system_rk4(des, dvars,
                                     ivar=self._curve_parameter,
                                     ics=ics,
                                     end_points=[t_min, t_max],
                                     step=step)

            # The value of 'step' being set by the user when calling
            # method 'solve', the value of (t_max - tmin)/step is not
            # necessarily an integer.
            # As a result, when the solver 'desolve_system_rk4' reaches
            # a curve parameter that is distant to 't_max' by less than
            # 'step', it computes one last point evaluated for a curve
            # parameter exactly equal to 't_max'.
            # Therefore, the difference between the curve parameter
            # corresponding to this last point and that corresponding
            # to the previous one is strictly less than 'step'. If this
            # difference is too small (that is, if the solver considered
            # that it did not reach 't_max', and hence computed one more
            # point, although it was already very close to 't_max'),
            # problems arise when using an interpolation of this
            # solution (such as getting points with coordinates 'nan').
            # As a result, we choose to remove the last point of a
            # solution when it is a point that was added by the solver
            # and threatens to be too close to the previous one
            # (arbitrarily, we consider two points to be too close if
            # their curve parameters are separated by less than 90% of a
            # step).
            if len(sol) > 1 and abs(sol[-1][0] - sol[-2][0]) < 0.9 * step:
                del sol[-1]

        elif method in ["odeint", "ode_int"]:
            # "ode_int" is here only for backward compatibility
            des = [fast_callable(eq, vars=tuple(list(self._chart[:])
                                                + self._velocities
                                                + [self._curve_parameter]),
                                 domain=float)
                   for eq in (self._velocities + eqns_num)]
            ics = initial_pt_coords + initial_tgt_vec_comps
            times = srange(t_min, t_max, step, include_endpoint=True)
            dvars = list(chart[:]) + self._velocities
            # Setting 1.e-10 as default value for the error control
            # parameters rtol and atol:
            if 'rtol' not in control_param:
                control_param['rtol'] = 1.e-10
            if 'atol' not in control_param:
                control_param['atol'] = 1.e-10
            sol0 = desolve_odeint(des, ics, times, dvars,
                                  ivar=self._curve_parameter, **control_param)

            # rewrite the solution to prepare for the extraction (which
            # removes information about the velocities), and convert
            # elements of type 'numpy.float64' to standard type 'float'

            import numpy as np
            sol = np.column_stack((times, sol0)) # tolist() done later

        elif method in ["dopri5", "dop853"]:
            import numpy as np
            des = [fast_callable(eq, vars=tuple(list(self._chart[:])
                                                + self._velocities), domain=float)
                   for eq in (self._velocities + eqns_num)]
            ics = initial_pt_coords + initial_tgt_vec_comps
            times = np.linspace(t_min, t_max, int((t_max-t_min)/step) + 1,
                                endpoint=True)
            # ode accepts a function returning a list, and not a list of functions
            r = ode(lambda t, y: [de(*y) for de in des]).set_integrator(method,
                                                               **control_param)
            r.set_initial_value(ics, t_min)
            r.set_solout(lambda t, y: 0 if chart.valid_coordinates_numerical(*y[0:dim])
                                      else -1)

            nt = len(times)
            sol0 = np.zeros((nt, 2*dim))
            sol0[0,:] = np.array(ics)
            for i in range(1, nt):
                sol0[i,:] = r.integrate(times[i])
                if not r.successful():
                    break
            sol = np.column_stack((times, sol0)) # tolist() done later

        elif method in ode_solver_methods:
            T = self._ode_solver

            if T is None:
                def system(t, y):
                    syst = self._velocities + eqns_num
                    par = self._curve_parameter
                    for i in range(dim):
                        vel = self._velocities[i]
                        syst[i] = syst[i].substitute({vel:y[dim+i]})
                        syst[dim+i] = syst[dim+i].substitute({par:t})
                        for j in range(dim):
                            coord = chart[:][j] # important to use '[:]'
                            # on 'chart' to avoid problems due to non
                            # zero starting index (i0)
                            veloc = self._velocities[j]
                            syst[dim+i] = syst[dim+i].substitute({coord:y[j]})
                            syst[dim+i] = syst[dim+i].substitute({veloc:y[dim+j]})
                    return syst
                from sage.calculus.ode import ode_solver
                T = ode_solver(function=system, **control_param)

            T.algorithm = method
            y_0 = initial_pt_coords + initial_tgt_vec_comps
            t_span = srange(t_min, t_max, step, include_endpoint=True)

            if method == "bsimp":
                # this method requires the expression of the Jacobian
                # matrix of the application defining the right-hand side
                # of the system to be provided

                if T.jacobian is None:
                    def jacobian(t,y):
                        jac = []
                        par = self._curve_parameter
                        for i in range(dim):
                            new_row = [0] * (2*dim)
                            new_row[dim + i] = 1
                            jac += [new_row]

                        for i in range(dim):
                            semi_row_coords = []
                            semi_row_vels = []
                            for j in range(dim):
                                coord = chart[:][j] # important to use
                                # '[:]' on 'chart' to avoid problems due
                                # to nonzero starting index (i0)
                                vel = self._velocities[j]
                                AUX = eqns_num[i].derivative(coord)
                                AUX2 = eqns_num[i].derivative(vel)
                                AUX = AUX.substitute({par: t})
                                AUX2 = AUX2.substitute({par: t})
                                for k in range(dim):
                                    coordin = chart[:][k] # important to
                                    # use '[:]' on 'chart' to avoid
                                    # problems due to nonzero starting
                                    # index (i0)
                                    veloc = self._velocities[k]
                                    AUX = AUX.substitute({coordin: y[k]})
                                    AUX = AUX.substitute({veloc: y[dim+k]})
                                    AUX2 = AUX2.substitute({coordin: y[k]})
                                    AUX2 = AUX2.substitute({veloc: y[dim+k]})
                                semi_row_coords += [AUX]
                                semi_row_vels += [AUX2]
                            jac += [semi_row_coords + semi_row_vels]

                        last_semi_row_coords = [0] * dim
                        last_semi_row_vels = []
                        for j in range(dim):
                            AUX3 = eqns_num[j].derivative(par)
                            AUX3 = AUX3.substitute({par: t})
                            for m in range(dim):
                                coordin = chart[:][m] # important to use
                                # '[:]' on 'chart' to avoid problems due
                                # to nonzero starting index (i0)
                                veloc = self._velocities[m]
                                AUX3 = AUX3.substitute({coordin: y[m]})
                                AUX3 = AUX3.substitute({veloc: y[dim+m]})
                            last_semi_row_vels += [AUX3]
                        jac += [last_semi_row_coords + last_semi_row_vels]
                        # 'AUX', 'AUX2' and 'AUX3' only used for the lines
                        # of source code to be shorter
                        return jac
                    T.jacobian = jacobian

            T.ode_solve(y_0=y_0, t_span=t_span)

            sol0 = T.solution
            sol = []
            for point in sol0:
                sol += [[point[0]] + point[1]]
            # above loop rewrites the solution in the same form than
            # that provided by other methods ('rk4_maxima' and
            # 'odeint'), in order to extract the time and corresponding
            # coordinate values a few lines below, in the same way for
            # all methods

        else:
            raise ValueError("no available method of integration " +
                             "referred to as '{}'".format(method))

        # eventually, extract the time and corresponding coordinate
        # values from each point of the solution computed (thus removing
        # information about the values of the velocities ; should the
        # latter be conserved ? They could turn useful in method
        # 'tangent_vector_eval_at', and in 'plot' when plotting the
        # tangent vectors.)

        if isinstance(sol, list):
            coords_sol = [point[0:dim + 1] for point in sol]
        else:
            coords_sol = sol[:, 0:dim + 1].tolist() # far faster in numpy

        if verbose:
            print("Numerical integration completed.\n\n" +
                  "Checking all points are in the chart domain...")

        N = len(coords_sol)
        n = 0
        while n < N and chart.valid_coordinates_numerical(*coords_sol[n][1:dim+1]):
            n += 1

        if n < N:
            raise ValueError("the {}th point ".format(n) +
                             "(initial point being the '0th' point) " +
                             "of the numerical solution (obtained " +
                             "for a curve parameter equal " +
                             "to {}) is out ".format(sol[n][0]) +
                             "of the chart domain; a curve with a " +
                             "smaller maximal value of the curve " +
                             "parameter, or a smaller initial tangent " +
                             "vector, might be considered. You can also try " +
                             "'solve_across_charts' in order not to be " +
                             "confined to a single chart")
        else:
            self._solutions[solution_key] = coords_sol
            if verbose:
                print("All points are in the chart domain.\n\n" +
                      "The resulting list of points was associated " +
                      "with the key '{}' ".format(solution_key) +
                      "(if this key already referred to a former " +
                      "numerical solution, such a solution was erased).")
            return self._solutions[solution_key]

    def solve_across_charts(self, charts=None, step=None, solution_key=None,
                            parameters_values=None, verbose=False,
                            **control_param):
        r"""
        Integrate the curve numerically over the domain of integration, with
        the ability to switch chart mid-integration.

        The only supported solver is :obj:`scipy:scipy.integrate.ode`,
        because it supports basic event handling, needed to detect when the
        curve is reaching the frontier of the chart. This is an adaptive step
        solver. So the ``step`` is not the step of integration but instead the
        step used to peak at the current chart, and switch if needed.

        INPUT:

        - ``step`` -- (default: ``None``) step of chart checking; default
          value is a hundredth of the domain of integration if none is
          provided. If your curve can't find a new frame on exiting the
          current frame, consider reducing this parameter.
        - ``charts`` -- (default: ``None``) list of chart allowed. The
          integration stops once it leaves those charts. By default the whole
          atlas is taken (only the top-charts).
        - ``solution_key`` -- (default: ``None``) key which the
          resulting numerical solution will be associated to; a default
          value is given if none is provided
        - ``parameters_values`` -- (default: ``None``) list of numerical
          values of the parameters present in the system defining the
          curve, to be substituted in the equations before integration
        - ``verbose`` -- boolean (default: ``False``); prints information about
          the computation in progress
        - ``**control_param`` -- extra control parameters to be passed to the
          solver

        OUTPUT: list of the numerical points of the computed solution

        EXAMPLES:

        Let us use :meth:`solve_across_charts` to integrate a geodesic of the
        Euclidean plane (a straight line) in polar coordinates.

        In pure polar coordinates `(r, \theta)`, artefacts can appear near
        the origin because of the fast variation of `\theta`, resulting in
        the direction of the geodesic being different before and after
        getting close to the origin.

        The solution to this problem is to switch to Cartesian coordinates
        near `(0,0)` to avoid any singularity.

        First let's declare the plane as a 2-dimensional manifold, with two
        charts `P` en `C` (for "Polar" and "Cartesian") and their transition
        maps::

            sage: M = Manifold(2, 'M', structure='Riemannian')
            sage: C.<x,y> = M.chart(coord_restrictions=lambda x,y: x**2+y**2 < 3**2)
            sage: P.<r,th> = M.chart(coord_restrictions=lambda r, th: r > 2)
            sage: P_to_C = P.transition_map(C,(r*cos(th), r*sin(th)))
            sage: C_to_P = C.transition_map(P,(sqrt(x**2+y**2), atan2(y,x)))

        Here we added restrictions on those charts, to avoid any
        singularity.  The intersection is the donut region `2 < r < 3`.

        We still have to define the metric. This is done in the Cartesian
        frame. The metric in the polar frame is computed automatically::

            sage: g = M.metric()
            sage: g[0,0,C]=1
            sage: g[1,1,C]=1
            sage: g[P.frame(), : ,P]
            [  1   0]
            [  0 r^2]

        To visualize our manifold, let's declare a mapping between every chart
        and the Cartesian chart, and then plot each chart in term of this
        mapping::

            sage: phi = M.diff_map(M, {(C,C): [x, y], (P,C): [r*cos(th), r*sin(th)]})
            sage: fig = P.plot(number_values=9, chart=C, mapping=phi,                   # needs sage.plot
            ....:              color='grey', ranges= {r:(2, 6), th:(0,2*pi)})
            sage: fig += C.plot(number_values=13, chart=C, mapping=phi,                 # needs sage.plot
            ....:               color='grey', ranges= {x:(-3, 3), y:(-3, 3)})

        There is a clear non-empty intersection between the two
        charts. This is the key point to successfully switch chart during the
        integration. Indeed, at least 2 points must fall in the intersection.

        .. RUBRIC:: Geodesic integration

        Let's define the time as `t`, the initial point as `p`, and the
        initial velocity vector as `v` (define as a member of the tangent
        space `T_p`). The chosen geodesic should enter the central region
        from the left and leave it to the right::

            sage: t = var('t')
            sage: p = M((5,pi+0.3), P)
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((-1,-0.03), P.frame().at(p))

        While creating the integrated geodesic, we need to specify the
        optional argument ``across_chart=True``, to prepare the compiled
        version of the changes of charts::

            sage: c = M.integrated_geodesic(g, (t, 0, 10), v, across_charts=True)

        The integration is done as usual, but using the method
        :meth:`solve_across_charts` instead of :meth:`solve`. This forces the
        use of :obj:`scipy:scipy.integrate.ode` as the solver, because of event
        handling support.

        The argument ``verbose=True`` will cause the solver to write a small
        message each time it is switching chart::

            sage: sol = c.solve_across_charts(step=0.1, verbose=True)
            Performing numerical integration with method 'ode'.
            Integration will take place on the whole manifold domain.
            Resulting list of points will be associated with the key 'ode_multichart' by default.
               ...
            Exiting chart, trying to switch to another chart.
            New chart found. Resuming integration.
            Exiting chart, trying to switch to another chart.
            New chart found. Resuming integration.
            Integration successful.

        As expected, two changes of chart occur.

        The returned solution is a list of pairs ``(chart, solution)``,
        where each solution is given on a unique chart, and the last
        point of a solution is the first of the next.

        The following code prints the corresponding charts::

            sage: for chart, solution in sol:
            ....:     print(chart)
            Chart (M, (r, th))
            Chart (M, (x, y))
            Chart (M, (r, th))

        The interpolation is done as usual::

            sage: interp = c.interpolate()

        To plot the result, you must first be sure that the mapping
        encompasses all the chart, which is the case here.
        You must also specify ``across_charts=True`` in order to call
        :meth:`plot_integrated` again on each part.
        Finally, ``color`` can be a list, which will be cycled through::

            sage: fig += c.plot_integrated(mapping=phi, color=["green","red"],          # needs sage.plot
            ....: thickness=3, plot_points=100, across_charts=True)
            sage: fig                                                                   # needs sage.plot
            Graphics object consisting of 43 graphics primitives

        .. PLOT::

            M = Manifold(2, 'M', structure='Riemannian')
            C= M.chart(names = ("x", "y"), coord_restrictions=lambda x,y: x**2+y**2 < 3**2)
            x, y = C[:]
            P = M.chart(names = ("r", "th"), coord_restrictions=lambda r,th: r > 2)
            r, th = P[:]
            P_to_C = P.transition_map(C,(r*cos(th), r*sin(th)))
            C_to_P = C.transition_map(P,(sqrt(x**2+y**2), atan2(y,x)))
            g = M.metric()
            g[0,0,C] = 1
            g[1,1,C] = 1
            g[P.frame(), : , P]
            phi = M.diff_map(M, {(C,C): [x, y], (P,C): [r*cos(th), r*sin(th)]})
            fig = P.plot(number_values=9, chart=C, mapping=phi, color='grey',
                         ranges= {r:(2, 6), th:(0,2*pi)})
            fig += C.plot(number_values=13, chart=C, mapping=phi, color='grey',
                          ranges= {x:(-3, 3), y:(-3, 3)})
            t = var('t')
            p = M((5,pi+0.3), P)
            Tp = M.tangent_space(p)
            v = Tp((-1,-0.03), P.frame().at(p))
            c = M.integrated_geodesic(g, (t, 0, 10), v, across_charts=True)
            sol = c.solve_across_charts(step=0.1, verbose=True)
            interp = c.interpolate()
            fig += c.plot_integrated(mapping=phi, color=["green","red"],
                        thickness=3, plot_points=100, across_charts=True)
            sphinx_plot(fig)
        """
        import numpy as np

        if verbose:
            print("Performing numerical integration with method 'ode'.")

        if charts is None:
            charts = self._codomain.top_charts()
            if verbose:
                print("Integration will take place on the whole manifold domain.")
        else:
            for c in charts:
                if not isinstance(c, Chart) or c.domain() is not self._codomain:
                    raise ValueError("'charts' needs to be a list of "
                                     "charts of the manifold")
            print("Integration will take place on {} charts.".format(len(charts)))

        if solution_key is None:
            solution_key = "ode_multichart"
            if verbose:
                print("Resulting list of points will be associated " +
                      "with the key '{}' ".format(solution_key) +
                      "by default.")
                print("   ...")

        t_min = self.domain().lower_bound()
        t_max = self.domain().upper_bound()

        eqns_num = self._equations_rhs.copy()

        v0 = self._initial_tangent_vector

        initial_tgt_space = v0.parent()
        initial_pt = initial_tgt_space.base_point()

        # Find a suitable initial chart, ie top chart in which the coordinates
        # of the initial point are known.

        for ichart in set(initial_pt._coordinates.keys()).intersection(charts):

            initial_chart = ichart

            initial_pt_coords = list(initial_pt.coordinates(initial_chart))
            initial_coord_basis = initial_chart.frame().at(initial_pt)
            initial_tgt_vec_comps = list(v0[initial_coord_basis, :])

            if step is None:
                step = (t_max - t_min) / 100

            dim = self.codomain().dim()

            if self._parameters:
                if parameters_values is None or len(parameters_values) != len(self._parameters):
                    raise ValueError("numerical values should be " +
                                     "provided for each of the " +
                                     "parameters "
                                     "{}".format(sorted(self._parameters, key=str)))
                for key in parameters_values:
                    parameters_values[key] = numerical_approx(parameters_values[key])

                if isinstance(t_min, Expression):
                    t_min = parameters_values[t_min]
                    if t_min == -Infinity or t_min == +Infinity:
                        raise ValueError("both boundaries of the " +
                                          "interval need to be finite")

                if isinstance(t_max, Expression):
                    t_max = parameters_values[t_max]
                    if t_max == -Infinity or t_max == +Infinity:
                        raise ValueError("both boundaries of the " +
                                         "interval need to be finite")

                for i in range(dim):
                    for chart in eqns_num:
                        if isinstance(eqns_num[chart][i], Expression):
                            eqns_num[chart][i] = eqns_num[chart][i].substitute(parameters_values)

                for i in range(dim):
                    if isinstance(initial_pt_coords[i], Expression):
                        AUX = initial_pt_coords[i]
                        AUX = AUX.substitute(parameters_values)
                        initial_pt_coords[i] = AUX
                    if isinstance(initial_tgt_vec_comps[i], Expression):
                        AUX2 = initial_tgt_vec_comps[i]
                        AUX2 = AUX2.substitute(parameters_values)
                        initial_tgt_vec_comps[i] = AUX2

            step = numerical_approx(step)

            initial_pt_coords = [numerical_approx(coord) for coord
                                 in initial_pt_coords]
            initial_tgt_vec_comps = [numerical_approx(comp) for comp
                                     in initial_tgt_vec_comps]

            t_min = numerical_approx(t_min)
            t_max = numerical_approx(t_max)

            if initial_chart.valid_coordinates(*initial_pt_coords):
                # found acceptable initial chart
                break

        else:
            # No initial chart found
            raise ValueError("initial point should be in the " +
                             "domain of its chart")

        # Transformation to fast_callable happens here
        des = {chart: [fast_callable(SR(eq), vars=tuple(
            list(chart[:]) + chart.symbolic_velocities()), domain=float)
               for eq in (chart.symbolic_velocities() + eqns_num[chart])]
               for chart in charts}

        ics = initial_pt_coords + initial_tgt_vec_comps
        times = np.linspace(t_min, t_max, int((t_max - t_min) / step) + 1,
                            endpoint=True)
        nt = len(times)

        sol = []

        chart = initial_chart

        start_index = 0  # current index while entering each new chart
        sol_chart = np.zeros((nt, 2 * dim))  # current chart solution
        sol_chart[0, :] = np.array(ics)  # starting with initial condition

        # Current equation to integrate, with initial and stop conditions
        r = ode(lambda t, y: [de(*y) for de in des[chart]]).set_integrator('dopri5',
                                                               **control_param)
        r.set_initial_value(ics, t_min)
        r.set_solout(lambda t, y: 0 if chart.valid_coordinates_numerical(*y[0:dim]) else -1)

        i = 1
        tried_charts = set()  # set of charts already searched at this step

        # Integration loop
        while i < nt:

            current_sol = r.integrate(times[i])
            if not r.successful():
                raise RuntimeError("unsuccessful integration")

            # step leads outside of the chart domain
            if abs(r.t-times[i]) > 1e-8:
                if verbose:
                    print("Exiting chart, trying to switch to another chart.")

                # Last known point
                last_pts = sol_chart[i-2-start_index, :dim]
                last_vel = sol_chart[i-2-start_index, dim:]

                random_order = list(set(charts).difference(tried_charts))
                shuffle(random_order)
                for new_chart in random_order:
                    tried_charts.add(new_chart)
                    if new_chart not in chart._subcharts:  # includes new != old

                        inter = chart.domain().intersection(new_chart.domain())

                        # The change of chart is performed here
                        new_pts = [f(*last_pts) for f in
                            self._fast_changes_of_chart[(chart.restrict(inter),
                                        new_chart.restrict(inter))]]
                        # If this line throws an error, check your changes
                        # of chart

                        if new_chart.valid_coordinates_numerical(*new_pts):
                            if verbose:
                                print("New chart found. Resuming integration.")
                            if start_index != i - 1:  # len(1) solution are ditched
                                # col-stack the times
                                sol_stacked = np.column_stack((times[start_index:i-1],
                                                sol_chart[:i-start_index-1, :]))
                                # add it to the global solution
                                sol.append((chart, sol_stacked))

                            # unfortunately building the tangent space is too
                            # slow, so we have to cheat a little and apply the
                            # change of frame manually (with a precompiled
                            # function)

                            new_vel = self._fast_changes_of_frame[(new_chart.frame().restrict(inter),
                               chart.frame().restrict(inter))](last_pts, last_vel)

                            ics = new_pts + new_vel
                            chart = new_chart

                            start_index = i - 1
                            sol_chart = np.zeros((nt, 2 * dim))
                            sol_chart[0, :] = np.array(ics)

                            r = ode(lambda t, y: [de(*y) for de in des[chart]])\
                                .set_integrator('dopri5')
                            r.set_initial_value(ics, times[i - 1])
                            r.set_solout(lambda t, y: 0 if chart.
                                valid_coordinates_numerical(*y[0:dim]) else -1)
                            i -= 1  # go back in the past to redo failed step
                            break
                # every chart was tried
                else:
                    if verbose:
                        print("No chart found, stopping integration.")
                        # col-stack the times
                    sol_chart = np.column_stack((times[start_index:i-1],
                                        sol_chart[:i-start_index-1, :]))
                    # add it to the global solution
                    sol.append((chart, sol_chart))
                    break

            # the integration step was successful
            else:
                sol_chart[i-start_index, :] = current_sol  # register the result
                tried_charts.clear()  # the set is reset.

            i += 1

        else:           # integration finishes successfully
            if verbose:
                print("Integration successful.")
            # col-stack the times
            sol_chart = np.column_stack((times[start_index:i-1],
                                         sol_chart[:i-start_index-1, :]))
            # add it to the global solution
            sol.append((chart, sol_chart))

        coords_sol = []
        for chart, chart_sol in sol:
            coords_sol.append((chart, chart_sol[:, 0:dim + 1])) # remove velocities

        self._solutions[solution_key] = coords_sol

        return self._solutions[solution_key]

    def solution(self, solution_key=None, verbose=False):
        r"""
        Return the solution (list of points) associated with the given
        key.

        INPUT:

        - ``solution_key`` -- (default: ``None``) key which the
          requested numerical solution is associated to; a default
          value is chosen if none is provided
        - ``verbose`` -- boolean (default: ``False``); prints information about
          the solution returned

        OUTPUT: list of the numerical points of the solution requested

        EXAMPLES:

        Requesting a numerical solution previously computed::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: [t, B_0, m, q, L, T] = var('t B_0 m q L T')
            sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B/m*D[1], -q*B/m*D[0], 0]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t,0,5), v, name='c')

            sage: # needs scipy
            sage: sol = c.solve(solution_key='sol_T1',
            ....:        parameters_values={B_0:1, m:1, q:1, L:10, T:1})
            sage: sol_bis = c.solution(verbose=True)
            Returning the numerical solution associated with the key
             'sol_T1' by default...
            sage: sol_bis == sol
            True
            sage: sol_ter = c.solution(solution_key='sol_T1')
            sage: sol_ter == sol
            True
            sage: sol_mute = c.solution()
            sage: sol_mute == sol
            True
        """
        if solution_key is None:
            if 'odeint' in self._solutions:
                solution_key = 'odeint'
            else:
                solution_key = next(iter(self._solutions))
                # will raise an error if self._solutions is empty
            if verbose:
                print("Returning the numerical solution associated " +
                      "with the key '{}' ".format(solution_key) +
                      "by default...")
        elif solution_key not in self._solutions:
            raise ValueError("no existing key " +
                             "'{}' ".format(solution_key) +
                             "referring to any numerical solution")

        return self._solutions[solution_key]

    def interpolate(self, solution_key=None, method=None,
                    interpolation_key=None, verbose=False):
        r"""
        Interpolate the chosen numerical solution using the given
        interpolation method.

        INPUT:

        - ``solution_key`` -- (default: ``None``) key which the
          numerical solution to interpolate is associated to ; a default
          value is chosen if none is provided
        - ``method`` -- (default: ``None``) interpolation scheme to use;
          algorithms available are

          * ``'cubic spline'``, which makes use of ``GSL`` via
            :class:`~sage.calculus.interpolation.Spline`

        - ``interpolation_key`` -- (default: ``None``) key which the
          resulting interpolation will be associated to ; a default
          value is given if none is provided
        - ``verbose`` -- boolean (default: ``False``); prints information about
          the interpolation in progress

        OUTPUT: built interpolation object

        EXAMPLES:

        Interpolating a numerical solution previously computed::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: [t, B_0, m, q, L, T] = var('t B_0 m q L T')
            sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B/m*D[1], -q*B/m*D[0], 0]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t,0,5), v, name='c')

            sage: # needs scipy
            sage: sol = c.solve(method='odeint',
            ....:        solution_key='sol_T1',
            ....:        parameters_values={B_0:1, m:1, q:1, L:10, T:1})
            sage: interp = c.interpolate(method='cubic spline',
            ....:                        solution_key='sol_T1',
            ....:                        interpolation_key='interp_T1',
            ....:                        verbose=True)
            Interpolation completed and associated with the key
             'interp_T1' (if this key already referred to a former
             interpolation, such an interpolation was erased).
            sage: interp = c.interpolate(verbose=True)
            Interpolating the numerical solution associated with the
             key 'sol_T1' by default...
            Performing cubic spline interpolation by default...
            Resulting interpolation will be associated with the key
             'cubic spline-interp-sol_T1' by default.
            Interpolation completed and associated with the key
             'cubic spline-interp-sol_T1' (if this key already referred
             to a former interpolation, such an interpolation was
             erased).

        TESTS::

            sage: interp = c.interpolate(solution_key='my solution')
            Traceback (most recent call last):
            ...
            ValueError: no existing key 'my solution' referring to any
             numerical solution
            sage: interp = c.interpolate(solution_key='sol_T1',
            ....:                        method='my method')
            Traceback (most recent call last):
            ...
            ValueError: no available method of interpolation referred to
             as 'my method'
        """
        if solution_key is None:
            if 'odeint' in self._solutions:
                solution_key = 'odeint'
            else:
                solution_key = next(iter(self._solutions)) # will raise
                # error if self._solutions empty
            if verbose:
                print("Interpolating the numerical solution " +
                      "associated with the key " +
                      "'{}' ".format(solution_key) +
                      "by default...")
        elif solution_key not in self._solutions:
            raise ValueError("no existing key " +
                             "'{}' ".format(solution_key) +
                             "referring to any numerical solution")

        if method is None:
            method = 'cubic spline'
            if verbose:
                print("Performing cubic spline interpolation by "
                      "default...")

        if interpolation_key is None:
            interpolation_key = "{}-interp-".format(method)
            interpolation_key += "{}".format(solution_key)
            if verbose:
                print("Resulting interpolation will be associated " +
                      "with the key '{}' ".format(interpolation_key) +
                      "by default.")

        if method == 'cubic spline':
            self._interpolations[interpolation_key] = []
            dim = self.codomain().dim()
            if not isinstance(self._solutions[solution_key][0], tuple):
                for i in range(dim):
                    coordinate_curve = []
                    for point in self._solutions[solution_key]:
                        coordinate_curve += [[point[0], point[i+1]]]
                    self._interpolations[interpolation_key] += [Spline(coordinate_curve)]
            else:   # case multi charts
                j = 0
                for chart, sol in self._solutions[solution_key]:
                    interp_chart = []
                    for i in range(dim):
                        coordinate_curve = []
                        for point in sol:
                            coordinate_curve += [[point[0], point[i + 1]]]
                        interp_chart += [Spline(coordinate_curve)]
                    self._interpolations[interpolation_key] += [(chart, interp_chart)]
                    self._interpolations[interpolation_key+"_chart_"+str(j)] = interp_chart
                    j += 1
        else:
            raise ValueError("no available method of interpolation " +
                             "referred to as '{}'".format(method))

        if verbose:
            print("Interpolation completed and associated with the " +
                  "key '{}' ".format(interpolation_key) +
                  "(if this key already referred to a former " +
                  "interpolation, such an interpolation was erased).")

        return self._interpolations[interpolation_key]

    def interpolation(self, interpolation_key=None, verbose=False):
        r"""
        Return the interpolation object associated with the given key.

        INPUT:

        - ``interpolation_key`` -- (default: ``None``) key which the
          requested interpolation is associated to; a default
          value is chosen if none is provided
        - ``verbose`` -- boolean (default: ``False``); prints information about
          the interpolation object returned

        OUTPUT: requested interpolation object

        EXAMPLES:

        Requesting an interpolation object previously computed::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: [t, B_0, m, q, L, T] = var('t B_0 m q L T')
            sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B/m*D[1], -q*B/m*D[0], 0]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t,0,5), v, name='c')

            sage: # needs scipy
            sage: sol = c.solve(method='odeint',
            ....:        solution_key='sol_T1',
            ....:        parameters_values={B_0:1, m:1, q:1, L:10, T:1})
            sage: interp = c.interpolate(method='cubic spline',
            ....:                         solution_key='sol_T1',
            ....:                         interpolation_key='interp_T1')
            sage: default_interp = c.interpolation(verbose=True)
            Returning the interpolation associated with the key
             'interp_T1' by default...
            sage: default_interp == interp
            True
            sage: interp_mute = c.interpolation()
            sage: interp_mute == interp
            True

        TESTS::

            sage: c.interpolation(interpolation_key='my interp')
            Traceback (most recent call last):
            ...
            ValueError: no existing key 'my interp' referring to any
             interpolation
        """

        if interpolation_key is None:
            if 'cubic spline' in self._interpolations:
                interpolation_key = 'cubic spline'
            else:
                interpolation_key = next(iter(self._interpolations))  # will
                # raise error if self._interpolations empty
            if verbose:
                print("Returning the interpolation associated with " +
                      "the key '{}' ".format(interpolation_key) +
                      "by default...")
        elif interpolation_key not in self._interpolations:
            raise ValueError("no existing key " +
                             "'{}' ".format(interpolation_key) +
                             "referring to any interpolation")

        return self._interpolations[interpolation_key]

    def __call__(self, t, interpolation_key=None,
                 verbose=False):
        r"""
        Return the image of the curve for the given value of the curve
        parameter, using the chosen interpolation.

        INPUT:

        - ``t`` -- curve parameter value at which the curve is evaluated
        - ``interpolation_key`` -- (default: ``None``) key which the
          interpolation requested to compute the point is associated to;
          a default value is chosen if none is provided
        - ``verbose`` -- boolean (default: ``False``); prints information about
          the interpolation used

        OUTPUT:

        - :class:`~sage.manifolds.point.ManifoldPoint` on a
          manifold (codomain) with numerical coordinates

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: [t, B_0, m, q, L, T] = var('t B_0 m q L T')
            sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B/m*D[1], -q*B/m*D[0], 0]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t,0,5), v, name='c')

            sage: # needs scipy
            sage: sol = c.solve(method='odeint',
            ....:        solution_key='sol_T1',
            ....:        parameters_values={B_0:1, m:1, q:1, L:10, T:1})
            sage: interp = c.interpolate(method='cubic spline',
            ....:                         solution_key='sol_T1',
            ....:                         interpolation_key='interp_T1')
            sage: c(1.1, interpolation_key='my interp')
            Traceback (most recent call last):
            ...
            ValueError: no existing key 'my interp' referring to any
             interpolation
            sage: p = c(1.1, verbose=True)
            Evaluating point coordinates from the interpolation
             associated with the key 'interp_T1' by default...
            sage: p.coordinates()     # abs tol 1e-12
            (1.060743337877276, -0.21538352256822146, 1.1)
        """

        if interpolation_key is None:
            if 'cubic spline' in self._interpolations:
                interpolation_key = 'cubic spline'
            else:
                # will raise error if self._interpolations empty
                interpolation_key = next(iter(self._interpolations))
            if verbose:
                print("Evaluating point coordinates from the " +
                  "interpolation associated with the key " +
                  "'{}' by default...".format(interpolation_key))
        elif interpolation_key not in self._interpolations:
            raise ValueError("no existing key " +
                             "'{}' ".format(interpolation_key) +
                             "referring to any interpolation")

        interpolation = self._interpolations[interpolation_key]

        if not isinstance(interpolation[0], Spline):
            # partial test, in case future interpolation objects do not
            # contain lists of instances of the Spline class
            raise TypeError("unexpected type of interpolation object")

        interpolated_coordinates = [coord_curve_spline(t)
                                    for coord_curve_spline in interpolation]
        return self.codomain().point(coords=interpolated_coordinates,
                                                  chart=self._chart)

    def tangent_vector_eval_at(self, t,
                               interpolation_key=None, verbose=False):
        r"""
        Return the vector tangent to ``self`` at the given curve
        parameter with components evaluated from the given
        interpolation.

        INPUT:

        - ``t`` -- curve parameter value at which the tangent vector is
          evaluated
        - ``interpolation_key`` -- (default: ``None``) key which the
          interpolation requested to compute the tangent vector is
          associated to; a default value is chosen if none is provided
        - ``verbose`` -- boolean (default: ``False``); prints information about
          the interpolation used

        OUTPUT:

        - :class:`~sage.manifolds.differentiable.tangent_vector.TangentVector`
          tangent vector with numerical components

        EXAMPLES:

        Evaluating a vector tangent to the curve::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: [t, B_0, m, q, L, T] = var('t B_0 m q L T')
            sage: B = B_0*t/T*exp(-(x1^2 + x2^2)/L^2)
            sage: D = X.symbolic_velocities()
            sage: eqns = [q*B/m*D[1], -q*B/m*D[0], 0]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t,0,5), v, name='c')

            sage: # needs scipy
            sage: sol = c.solve(method='odeint',
            ....:        solution_key='sol_T1',
            ....:        parameters_values={B_0:1, m:1, q:1, L:10, T:1})
            sage: interp = c.interpolate(method='cubic spline',
            ....:                         solution_key='sol_T1',
            ....:                         interpolation_key='interp_T1')
            sage: tg_vec = c.tangent_vector_eval_at(1.22, verbose=True)
            Evaluating tangent vector components from the interpolation
             associated with the key 'interp_T1' by default...
            sage: tg_vec
            Tangent vector at Point on the 3-dimensional differentiable
             manifold M
            sage: tg_vec[:]     # abs tol 1e-12
            [0.7392640422917979, -0.6734182509826023, 1.0]
            sage: tg_vec_mute = c.tangent_vector_eval_at(1.22,
            ....:                         interpolation_key='interp_T1')
            sage: tg_vec_mute == tg_vec
            True

        TESTS::

            sage: tg_vec = c.tangent_vector_eval_at(1.22,
            ....:                         interpolation_key='my interp')
            Traceback (most recent call last):
            ...
            ValueError: no existing key 'my interp' referring to any
             interpolation
        """
        if interpolation_key is None:
            if 'cubic spline' in self._interpolations:
                interpolation_key = 'cubic spline'
            else:
                # will raise error if self._interpolations empty
                interpolation_key = next(iter(self._interpolations))
            if verbose:
                print("Evaluating tangent vector components from the " +
                      "interpolation associated with the key " +
                      "'{}' by default...".format(interpolation_key))
        elif interpolation_key not in self._interpolations:
            raise ValueError("no existing key " +
                             "'{}' ".format(interpolation_key) +
                             "referring to any interpolation")

        interpolation = self._interpolations[interpolation_key]

        if not isinstance(interpolation[0], Spline):
            # partial test, in case future interpolation objects do not
            # contain lists of instances of the Spline class
            raise TypeError("unexpected type of interpolation object")

        interpolated_coordinates = [coordinate_curve_spline(t)
                       for coordinate_curve_spline in interpolation]
        M = self.codomain()
        p = M.point(interpolated_coordinates, chart=self._chart, name=None)
        Tp = M.tangent_space(p)

        # by default, order=1 in method 'derivative' of a class Spline
        evaluated_tgt_vec_comp = [coord_curve_spline.derivative(t)
                                  for coord_curve_spline in interpolation]
        basis = self._chart.frame().at(p)
        return Tp(evaluated_tgt_vec_comp, basis=basis)

    @options(thickness=1, plot_points=75, aspect_ratio='automatic',
             plot_points_tangent=10, width_tangent=1, scale=1)
    def plot_integrated(self, chart=None, ambient_coords=None,
             mapping=None, prange=None, interpolation_key=None,
             include_end_point=(True, True),
             end_point_offset=(0.001, 0.001), verbose=False, color='red',
             style='-', label_axes=True, display_tangent=False,
             color_tangent='blue', across_charts=False, **kwds):
        r"""
        Plot the 2D or 3D projection of ``self`` onto the space of the
        chosen two or three ambient coordinates, based on the
        interpolation of a numerical solution previously computed.

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.curve.DifferentiableCurve.plot`
            for complete information about the input.

        ADDITIONAL INPUT:

        - ``interpolation_key`` -- (default: ``None``) key associated to
          the interpolation object used for the plot; a default value
          is chosen if none is provided
        - ``verbose`` -- boolean (default: ``False``); prints information about
          the interpolation object used and the plotting in progress
        - ``display_tangent`` -- boolean (default: ``False``); determines whether
          some tangent vectors should also be plotted
        - ``color_tangent`` -- (default: ``blue``) color of the tangent
          vectors when these are plotted
        - ``plot_points_tangent`` -- (default: 10) number of tangent
          vectors to display when these are plotted
        - ``width_tangent`` -- (default: 1) sets the width of the arrows
          representing the tangent vectors
        - ``scale`` -- (default: 1) scale applied to the tangent vectors
          before displaying them

        EXAMPLES:

        Trajectory of a particle of unit mass and unit charge in an
        unit, axial, uniform, stationary magnetic field::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: var('t')
            t
            sage: D = X.symbolic_velocities()
            sage: eqns = [D[1], -D[0], 0]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t,0,6), v, name='c')

            sage: # needs scipy
            sage: sol = c.solve()
            sage: interp = c.interpolate()
            sage: c_plot_2d = c.plot_integrated(ambient_coords=[x1, x2],
            ....:                 thickness=2.5,
            ....:                 display_tangent=True, plot_points=200,
            ....:                 plot_points_tangent=10, scale=0.5,
            ....:                 color='blue', color_tangent='red',
            ....:                 verbose=True)
            Plotting from the interpolation associated with the key
             'cubic spline-interp-odeint' by default...
            A tiny final offset equal to 0.000301507537688442 was
             introduced for the last point in order to safely compute it
             from the interpolation.
            sage: c_plot_2d
            Graphics object consisting of 11 graphics primitives

        .. PLOT::

            M = Manifold(3, 'M')
            X = M.chart('x1 x2 x3'); x1, x2, x3 = X[:]
            D = X.symbolic_velocities()
            eqns = [D[1], -D[0], 0]
            p = M.point((0,0,0), name='p')
            Tp = M.tangent_space(p)
            v = Tp((1,0,1))
            t = var('t')
            c = M.integrated_curve(eqns, D, (t, 0, 6), v, name='c')
            sol = c.solve()
            interp = c.interpolate()
            c_plot_2d_1 = c.plot_integrated(ambient_coords=[x1, x2],
                            thickness=2.5,
                            display_tangent=True, plot_points=200,
                            plot_points_tangent=10, scale=0.5,
                            color='blue', color_tangent='red')
            sphinx_plot(c_plot_2d_1)
        """
        from sage.manifolds.chart import RealChart

        #
        # Get the @options plot_points from kwds
        #
        plot_points = kwds.pop('plot_points')

        #
        # Interpolation to use
        #
        if interpolation_key is None:
            if 'cubic spline' in self._interpolations:
                interpolation_key = 'cubic spline'
            else:
                if across_charts:
                    for key in self._interpolations:
                        if key[-8:-1] != '_chart_':       # check if not a subplot
                            interpolation_key = key
                            break
                    else:
                        raise ValueError("Did you forget to "
                                         "integrate or interpolate the result?")
                else:
                    interpolation_key = next(iter(self._interpolations)) #will
                # raise error if self._interpolations empty

            if verbose:
                print("Plotting from the interpolation associated " +
                      "with the key '{}' ".format(interpolation_key) +
                      "by default...")
        elif interpolation_key not in self._interpolations:
            raise ValueError("no existing key '{}' ".format(interpolation_key)
                             + "referring to any interpolation")

        interpolation = self._interpolations[interpolation_key]

        if across_charts:
            len_tot = sum(len(interp[1][0]) for interp in interpolation)
            if isinstance(color, list):
                color = color * (len(interpolation) // 3 + 1)
            else:
                color = color * len(interpolation)
            res = 0
            for i in range(len(interpolation)):
                nb_pts = int(float(plot_points)*len(interpolation[i][1][0])/len_tot)
                self._chart = interpolation[i][0]
                res += self.plot_integrated(chart=chart, ambient_coords=ambient_coords,
                                            mapping=mapping, prange=prange,
                                            interpolation_key=interpolation_key+"_chart_"+str(i),
                                            include_end_point=include_end_point,
                                            end_point_offset=end_point_offset,
                                            verbose=verbose, color=color[i],
                                            style=style, label_axes=False,
                                            display_tangent=display_tangent,
                                            color_tangent=color_tangent,
                                            across_charts=False,
                                            plot_points=nb_pts, **kwds)

            return res

        #
        # Get the remaining @options from kwds
        #
        thickness = kwds.pop('thickness')
        aspect_ratio = kwds.pop('aspect_ratio')

        #
        # The mapping, if present, and the chart with respect to which the curve
        # is plotted
        #
        if mapping is None:
            i0 = self.codomain().start_index()
            if chart is None:
                chart = self._chart
            else:
                if not isinstance(chart, RealChart):
                    raise TypeError("{} is not a real chart".format(chart))
                mapping = self.codomain().identity_map()
        else:
            i0 = mapping.codomain().start_index()
            if chart is None:
                chart = mapping.codomain().default_chart()
            elif not isinstance(chart, RealChart):
                raise TypeError("{} is not a real chart".format(chart))

        #
        # Coordinates of the above chart with respect to which the curve is
        # plotted
        #
        if ambient_coords is None:
            ambient_coords = chart[:]  # all chart coordinates are used
        n_pc = len(ambient_coords)
        if n_pc != 2 and n_pc != 3:
            raise ValueError("the number of coordinates involved in " +
                             "the plot must be either 2 or 3, " +
                             "not {}".format(n_pc))

        # From now on, 'pc' will denote coordinates in terms of which
        # the curve is plotted (i.e. the "ambient coordinates"), while
        # 'coord' will denote coordinates on self.domain(); of course,
        # when, for instance, the mapping is the identity map, these may
        # be the same.

        # indices of plot coordinates
        # will raise an error if ambient_coords are not associated with chart
        ind_pc = [chart[:].index(pc) + i0 for pc in ambient_coords]

        #
        # Maximal parameter range for the plot of the chosen
        # interpolation
        #
        # these two lines are the only general way to get the maximal
        # parameter range since, at this point, there is no clue about
        # the solution from which 'interpolation' was build, and it would
        # be an obvious error to declare param_min=self.domain().lower_bound()
        # for instance, since this might be an expression
        param_min = interpolation[0][0][0]
        param_max = interpolation[0][-1][0]

        if prange is None:
            prange = (param_min, param_max)
        elif not isinstance(prange, (tuple, list)):
            raise TypeError("{} is neither ".format(prange) +
                            "a tuple nor a list")
        elif len(prange) != 2:
            raise ValueError("the argument prange must be a " +
                             "tuple/list of 2 elements")
        else:
            p = prange #'p' declared only for the line below to be shorter
            if p[0] < param_min or p[0] > param_max or p[1] < param_min or p[1] > param_max:
                raise ValueError("parameter range should be a " +
                                 "subinterval of the curve domain " +
                                 "({})".format(self.domain()))

        tmin = numerical_approx(prange[0])
        tmax = numerical_approx(prange[1])

        if not include_end_point[0]:
            tmin += numerical_approx(end_point_offset[0])

        if not include_end_point[1]:
            tmax -= numerical_approx(end_point_offset[1])

        if mapping is None:
            if not isinstance(interpolation[0], Spline):
                # partial test in case future interpolation objects do not
                # contain lists of instances of the Spline class
                raise TypeError("unexpected type of interpolation object")

            #
            # List of points for the plot curve
            #
            plot_curve = []
            dt = (tmax - tmin) / (plot_points - 1)
            t = tmin

            for k in range(plot_points):
                if k == 0 and t < param_min:
                    # This might happen for the first point (i.e. k = 0)
                    # when prange[0], and hence tmin, should equal param_min;
                    # but mere numerical rounding coming from having taken
                    # tmin = numerical_approx(prange[0]) might
                    # raise errors from trying to evaluate the
                    # interpolation at a time smaller than
                    # self.domain.lower_bound(). Hence the line below
                    # that adds 1% of the step to compute even more
                    # safely the first point
                    t = param_min + 0.01*dt
                    if verbose:
                        print("A tiny initial offset equal to " +
                              "{} ".format(0.01*dt) +
                              "was introduced for the first point " +
                              "only, in order to safely compute " +
                              "it from the interpolation.")

                if k == plot_points-1 and t > param_max:
                    # This might happen for the last point
                    # (i.e. k = plot_points-1) when prange[1], and hence
                    # tmax, should equal param_max; but mere numerical
                    # rounding coming from having taken
                    # tmax = numerical_approx(prange[1) might raise errors
                    # from trying to evaluate the interpolation at a time
                    # greater than self.domain.upper_bound().
                    # Hence the line below that subtract 1% of the
                    # step to compute even more safely the last point
                    t = param_max - 0.01*dt
                    if verbose:
                        print("A tiny final offset equal to " +
                              "{} ".format(0.01*dt) +
                              "was introduced for the last point " +
                              "in order to safely compute " +
                              "it from the interpolation.")

                plot_curve.append([interpolation[j-i0](t) for j in ind_pc])

                if k == 0 and t > tmin:
                    # in case an initial offset was earlier added to
                    # 'tmin' in order to avoid errors, it is now needed
                    # to cancel this offset for the next steps
                    t = tmin
                t += dt

            if display_tangent:
                from sage.plot.graphics import Graphics
                from sage.plot.arrow import arrow2d
                from sage.plot.plot3d.shapes import arrow3d

                scale = kwds.pop('scale')
                plot_points_tangent = kwds.pop('plot_points_tangent')
                width_tangent = kwds.pop('width_tangent')

                plot_vectors = Graphics()
                dt = (tmax - tmin) / (plot_points_tangent - 1)
                t = tmin

                for k in range(plot_points_tangent):
                    if k == 0 and t < param_min:
                        # This might happen for the first point
                        # (i.e. k = 0) when prange[0], and hence tmin
                        # should equal param_min; but mere numerical
                        # rounding coming from having taken
                        # tmin = numerical_approx(prange[0]) might
                        # raise errors from trying to evaluate the
                        # interpolation at a time smaller than
                        # self.domain.lower_bound().
                        # Hence the line below that add 1% of the step
                        # to compute even more safely the first point.
                        t = param_min + 0.01*dt
                        if verbose:
                            print("A tiny initial offset equal to " +
                                  "{} ".format(0.01*dt) +
                                  "was introduced for the first point " +
                                  "only, in order to safely compute " +
                                  "it from the interpolation.")

                    if k == plot_points_tangent - 1 and t > param_max:
                        # This might happen for the last point
                        # (i.e. k = plot_points_tangent-1) when
                        # prange[1], and hence tmax, should equal
                        # param_max; but mere numerical rounding coming from
                        # having taken tmax = numerical_approx(prange[1)
                        # might raise errors from trying to evaluate the
                        # interpolation at a time greater than
                        # self.domain.upper_bound(). Hence the line below
                        # that subtracts 1% of the step to compute even
                        # more safely the last point.
                        t = param_max - 0.01*dt
                        if verbose:
                            print("A tiny final offset equal to " +
                                  "{} ".format(0.01*dt) +
                                  "was introduced for the last point " +
                                  "in order to safely compute " +
                                  "it from the interpolation.")

                    # interpolated ambient coordinates:
                    xp = [interpolation[j-i0](t) for j in ind_pc]

                    # tangent vector ambiant components evaluated
                    # from the interpolation:
                    vec = [coordinate_curve_spline.derivative(t)
                       for coordinate_curve_spline in interpolation]

                    coord_tail = xp
                    coord_head = [xp[j] + scale*vec[j]
                                  for j in range(len(xp))]

                    if coord_head != coord_tail:
                        if n_pc == 2:
                            plot_vectors += arrow2d(tailpoint=coord_tail,
                                                    headpoint=coord_head,
                                                    color=color_tangent,
                                                    width=width_tangent)
                        else:
                            plot_vectors += arrow3d(coord_tail,
                                                    coord_head,
                                                    color=color_tangent,
                                                    width=width_tangent)

                    if k == 0 and t > tmin:
                        # in case an initial offset was earlier added
                        # to 'tmin' in order to avoid errors, it is now
                        # needed to cancel this offset for the next steps
                        t = tmin
                    t += dt

                return plot_vectors + DifferentiableCurve._graphics(self,
                                         plot_curve, ambient_coords,
                                         thickness=thickness,
                                         aspect_ratio=aspect_ratio,
                                         color=color,
                                         style=style,
                                         label_axes=label_axes)

            return DifferentiableCurve._graphics(self, plot_curve,
                             ambient_coords, thickness=thickness,
                             aspect_ratio=aspect_ratio, color=color,
                             style=style, label_axes=label_axes)
        else:
            #
            # The coordinate expressions of the mapping and the
            # coordinates involved
            #
            for chart_pair in mapping._coord_expression:
                subs = (chart_pair[0]._subcharts, chart_pair[1]._subcharts)
                # 'subs' declared only for the line below to be shorter
                if self._chart in subs[0] and chart in subs[1]:
                    transf = {}
                    required_coords = set()
                    for pc in ambient_coords:
                        jpc = chart[:].index(pc)
                        AUX = mapping._coord_expression[chart_pair]
                        # 'AUX' used only for the lines of source code
                        # to be shorter
                        transf[pc] = AUX.expr()[jpc]
                        AUX2 = transf[pc].variables() # idem
                        required_coords = required_coords.union(AUX2)
                    break
            else:
                raise ValueError("no expression has been found for " +
                                 "{} in terms of {}".format(self,chart))

            # fastf is the fast version of a substitution + numerical evaluation
            # using fast_callable.
            fastf = [fast_callable(transf[chart[i]], vars=tuple(self._chart[:]))
                     for i in ind_pc]

            if not isinstance(interpolation[0], Spline):
                # partial test, in case future interpolation objects do not
                # contain lists of instances of the Spline class
                raise TypeError("unexpected type of interpolation object")

            #
            # List of points for the plot curve
            #
            plot_curve = []
            dt = (tmax - tmin) / (plot_points - 1)
            t = tmin
            required_coords_values = {}

            for k in range(plot_points):
                if k == 0 and t < param_min:
                    # This might happen for the first point (i.e. k = 0)
                    # when prange[0], and hence tmin, should equal param_min;
                    # but mere numerical rounding coming from having taken
                    # tmin = numerical_approx(prange[0]) might
                    # raise errors from trying to evaluate the
                    # interpolation at a time smaller than
                    # self.domain.lower_bound(). Hence the line below that adds
                    # 1% of the step to compute even more safely the first point
                    t = param_min + 0.01*dt
                    if verbose:
                        print("A tiny initial offset equal to " +
                              "{} ".format(0.01*dt) +
                              "was introduced for the first point " +
                              "only, in order to safely compute " +
                              "it from the interpolation.")

                if k == plot_points - 1 and t > param_max:
                    # This might happen for the last point (i.e. k = plot_points-1)
                    # when prange[1], and hence tmax, should equal
                    # param_max; but mere numerical rounding coming from
                    # having taken tmax = numerical_approx(prange[1)
                    # might raise errors from trying to evaluate the
                    # interpolation at a time greater than
                    # self.domain.upper_bound(). Hence the line below that
                    # subtracts 1% of the step to compute even more safely
                    # the last point.
                    t = param_max - 0.01*dt
                    if verbose:
                        print("A tiny final offset equal to " +
                              "{} ".format(0.01*dt) +
                              "was introduced for the last point " +
                              "in order to safely compute " +
                              "it from the interpolation.")

                # list of coordinates, argument of fastf, the fast diff_map
                arg = [inter(t) for inter in interpolation]
                # evaluation of fastf
                xp = [fastf[j](*arg) for j in range(len(ambient_coords))]
                plot_curve.append(xp)

                if k == 0 and t > tmin:
                    # in case an initial offset was earlier added to
                    # 'tmin' in order to avoid errors, it is now needed
                    # to cancel this offset for the next steps
                    t = tmin

                t += dt

            if display_tangent:
                from sage.plot.graphics import Graphics
                from sage.plot.arrow import arrow2d
                from sage.plot.plot3d.shapes import arrow3d

                scale = kwds.pop('scale')
                plot_points_tangent = kwds.pop('plot_points_tangent')
                width_tangent = kwds.pop('width_tangent')

                plot_vectors = Graphics()
                dt = (tmax - tmin) / (plot_points_tangent - 1)
                t = tmin
                Dcoord_Dt = {}

                Dpc_Dcoord = {}
                for pc in ambient_coords:
                    Dpc_Dcoord[pc] = {}
                    for coord in transf[pc].variables():
                        Dpc_Dcoord[pc][coord] = transf[pc].derivative(coord)

                for k in range(plot_points_tangent):
                    if k == 0 and t < param_min:
                        # This might happen for the first point (i.e. k = 0)
                        # when prange[0], and hence tmin, should equal param_min;
                        # but mere numerical rounding coming from having taken
                        # tmin = numerical_approx(prange[0]) might
                        # raise errors from trying to evaluate the
                        # interpolation at a time smaller than
                        # self.domain.lower_bound(). Hence the line below
                        # that adds 1% of the step to compute even more
                        # safely the first point
                        t = param_min + 0.01*dt
                        if verbose:
                            print("A tiny initial offset equal to " +
                                  "{} ".format(0.01*dt) +
                                  "was introduced for the first point " +
                                  "only, in order to safely compute " +
                                  "it from the interpolation.")

                    if k == plot_points_tangent - 1 and t > param_max:
                        # This might happen for the last point
                        # (i.e. k = plot_points_tangent-1) when
                        # when prange[1], and hence tmax, should equal
                        # param_max; but mere numerical rounding coming from
                        # having taken tmax = numerical_approx(prange[1)
                        # might raise errors from trying to evaluate the
                        # interpolation at a time greater than
                        # self.domain.upper_bound(). Hence the line below
                        # that subtracts 1% of the step to compute even
                        # more safely the last point
                        t = param_max - 0.01*dt
                        if verbose:
                            print("A tiny final offset equal to " +
                                  "{} ".format(0.01*dt) +
                                  "was introduced for the last point " +
                                  "in order to safely compute " +
                                  "it from the interpolation.")

                    for coord in required_coords:
                        i = self._chart[:].index(coord)
                        AUX = interpolation[i] # 'AUX' only used
                        # for the lines below to be shorter
                        required_coords_values[coord] = AUX(t)
                        Dcoord_Dt[coord] = AUX.derivative(t)

                    xp = []
                    pushed_vec = []
                    for j in ind_pc:
                        pc = chart[j]
                        AUX = transf[pc]
                        AUX = AUX.substitute(required_coords_values)
                        # 'AUX' only used for the lines of code to
                        # be shorter
                        xp += [numerical_approx(AUX)]

                        pushed_comp = 0
                        for coord in transf[pc].variables():
                            D = Dpc_Dcoord[pc][coord]
                            D = D.substitute(required_coords_values)
                            D = numerical_approx(D)
                            pushed_comp += Dcoord_Dt[coord] * D

                        pushed_vec += [pushed_comp]

                    coord_tail = xp
                    coord_head = [val + scale*pushed_vec[j]
                                  for j, val in enumerate(xp)]

                    if coord_head != coord_tail:
                        if n_pc == 2:
                            plot_vectors += arrow2d(tailpoint=coord_tail,
                                                    headpoint=coord_head,
                                                    color=color_tangent,
                                                    width=width_tangent)
                        else:
                            plot_vectors += arrow3d(coord_tail,
                                                    coord_head,
                                                    color=color_tangent,
                                                    width=width_tangent)

                    if k == 0 and t > tmin:
                        # in case an initial offset was earlier added to
                        # 'tmin' in order to avoid errors, it is now needed
                        # to cancel this offset for the next steps
                        t = tmin

                    t += dt
                return plot_vectors + DifferentiableCurve._graphics(self,
                                         plot_curve, ambient_coords,
                                         thickness=thickness,
                                         aspect_ratio=aspect_ratio,
                                         color=color,
                                         style=style,
                                         label_axes=label_axes)
            return DifferentiableCurve._graphics(self, plot_curve,
                             ambient_coords, thickness=thickness,
                             aspect_ratio=aspect_ratio, color=color,
                             style=style, label_axes=label_axes)

class IntegratedAutoparallelCurve(IntegratedCurve):
    r"""
    Autoparallel curve on the manifold with respect to a given
    affine connection.

    INPUT:

    - ``parent`` --
      :class:`~sage.manifolds.differentiable.manifold_homset.IntegratedAutoparallelCurveSet`
      the set of curves `\mathrm{Hom_{autoparallel}}(I, M)` to which the
      curve belongs
    - ``affine_connection`` --
      :class:`~sage.manifolds.differentiable.affine_connection.AffineConnection`
      affine connection with respect to which the curve is autoparallel
    - ``curve_parameter`` -- symbolic expression to be used as the
      parameter of the curve (the equations defining an instance of
      IntegratedAutoparallelCurve are such that ``t`` will actually be
      an affine parameter of the curve)
    - ``initial_tangent_vector`` --
      :class:`~sage.manifolds.differentiable.tangent_vector.TangentVector`
      initial tangent vector of the curve
    - ``chart`` -- (default: ``None``) chart on the manifold in terms of
      which the equations are expressed; if ``None`` the default chart
      of the manifold is assumed
    - ``name`` -- (default: ``None``) string; symbol given to the curve
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to
      denote the curve; if none is provided, ``name`` will be used

    EXAMPLES:

    Autoparallel curves associated with the Mercator projection of the
    unit 2-sphere `\mathbb{S}^{2}`.

    .. SEEALSO::

        https://idontgetoutmuch.wordpress.com/2016/11/24/mercator-a-connection-with-torsion/
        for more details about Mercator projection.

    On the Mercator projection, the lines of longitude all appear
    vertical and then all parallel with respect to each other.
    Likewise, all the lines of latitude appear horizontal and parallel
    with respect to each other.
    These curves may be recovered as autoparallel curves of a certain
    connection `\nabla` to be made explicit.

    Start with declaring the standard polar coordinates
    `(\theta, \phi)` on `\mathbb{S}^{2}` and the
    corresponding coordinate frame `(e_{\theta}, e_{\phi})`::

        sage: S2 = Manifold(2, 'S^2', start_index=1)
        sage: polar.<th,ph>=S2.chart()
        sage: epolar = polar.frame()

    Normalizing `e_{\phi}` provides an orthonormal basis::

        sage: ch_basis = S2.automorphism_field()
        sage: ch_basis[1,1], ch_basis[2,2] = 1, 1/sin(th)
        sage: epolar_ON = epolar.new_frame(ch_basis,'epolar_ON')

    Denote `(\hat{e}_{\theta}, \hat{e}_{\phi})` such an orthonormal frame
    field. In any point, the vector field `\hat{e}_{\theta}` is
    normalized and tangent to the line of longitude through the point.
    Likewise, `\hat{e}_{\phi}` is normalized and tangent to the
    line of latitude.

    Now, set an affine connection with respect to such fields that are
    parallelly transported in all directions, that is:
    `\nabla \hat{e}_{\theta} = \nabla \hat{e}_{\phi} = 0`.
    This is equivalent to setting all the connection coefficients to
    zero with respect to this frame::

        sage: nab = S2.affine_connection('nab')
        sage: nab.set_coef(frame=epolar_ON)[:]
        [[[0, 0], [0, 0]], [[0, 0], [0, 0]]]

    This connection is such that two vectors are parallel if their
    angles to a given meridian are the same.
    Check that this connection is compatible with the Euclidean
    metric tensor `g` induced on `\mathbb{S}^{2}`::

        sage: g = S2.metric('g')
        sage: g[1,1], g[2,2] = 1, (sin(th))^2
        sage: nab(g)[:]
        [[[0, 0], [0, 0]], [[0, 0], [0, 0]]]

    Yet, this connection is not the Levi-Civita connection, which
    implies that it has non-vanishing torsion::

        sage: nab.torsion()[:]
        [[[0, 0], [0, 0]], [[0, cos(th)/sin(th)], [-cos(th)/sin(th), 0]]]

    Set generic initial conditions for the autoparallel curves to
    compute::

        sage: [th0, ph0, v_th0, v_ph0] = var('th0 ph0 v_th0 v_ph0')
        sage: p = S2.point((th0, ph0), name='p')
        sage: Tp = S2.tangent_space(p)
        sage: v = Tp((v_th0, v_ph0), basis=epolar_ON.at(p))

    Note here that the components ``(v_th0, v_ph0)`` of the initial
    tangent vector ``v`` refer to the basis
    ``epolar_ON`` `= (\hat{e}_{\theta}, \hat{e}_{\phi})`
    and not the coordinate basis ``epolar`` `= (e_{\theta}, e_{\phi})`.
    This is merely to help picture the aspect of the tangent vector in
    the usual embedding of `\mathbb{S}^{2}` in
    `\mathbb{R}^{3}` thanks to using an orthonormal frame,
    since providing the components with respect to the coordinate basis
    would require multiplying the second component (i.e. the `\phi`
    component) in order to picture the vector in the same way.
    This subtlety will need to be taken into account later when the
    numerical curve will be compared to the analytical solution.

    Now, declare the corresponding integrated autoparallel curve and display
    the differential system it satisfies::

        sage: [t, tmin, tmax] = var('t tmin tmax')
        sage: c = S2.integrated_autoparallel_curve(nab, (t, tmin, tmax),
        ....:                                  v, chart=polar, name='c')
        sage: sys = c.system(verbose=True)
        Autoparallel curve c in the 2-dimensional differentiable
         manifold S^2 equipped with Affine connection nab on the
         2-dimensional differentiable manifold S^2, and integrated over
         the Real interval (tmin, tmax) as a solution to the following
         equations, written with respect to Chart (S^2, (th, ph)):
        <BLANKLINE>
        Initial point: Point p on the 2-dimensional differentiable
         manifold S^2 with coordinates [th0, ph0] with respect to
         Chart (S^2, (th, ph))
        Initial tangent vector: Tangent vector at Point p on the
         2-dimensional differentiable manifold S^2 with
         components [v_th0, v_ph0/sin(th0)] with respect to Chart (S^2, (th, ph))
        <BLANKLINE>
        d(th)/dt = Dth
        d(ph)/dt = Dph
        d(Dth)/dt = 0
        d(Dph)/dt = -Dph*Dth*cos(th)/sin(th)
        <BLANKLINE>

    Set a dictionary providing the parameter range and the initial
    conditions for a line of latitude and a line of longitude::

        sage: dict_params={'latit':{tmin:0,tmax:3,th0:pi/4,ph0:0.1,v_th0:0,v_ph0:1},
        ....:   'longi':{tmin:0,tmax:3,th0:0.1,ph0:0.1,v_th0:1,v_ph0:0}}

    Declare the Mercator coordinates `(\xi, \zeta)` and the
    corresponding coordinate change from the polar coordinates::

        sage: mercator.<xi,ze> = S2.chart(r'xi:(-oo,oo):\xi ze:(0,2*pi):\zeta')
        sage: polar.transition_map(mercator, (log(tan(th/2)), ph))
        Change of coordinates from Chart (S^2, (th, ph)) to Chart
         (S^2, (xi, ze))

    Ask for the identity map in terms of these charts in order to add
    this coordinate change to its dictionary of expressions. This is
    required to plot the curve with respect to the Mercator chart::

        sage: identity = S2.identity_map()
        sage: identity.coord_functions(polar, mercator)
        Coordinate functions (log(sin(1/2*th)/cos(1/2*th)), ph) on the
         Chart (S^2, (th, ph))

    Solve, interpolate and prepare the plot for the solutions
    corresponding to the two initial conditions previously set::

        sage: graph2D_mercator = Graphics()
        sage: for key in dict_params:
        ....:     sol = c.solve(solution_key='sol-'+key,
        ....:                        parameters_values=dict_params[key])
        ....:     interp = c.interpolate(solution_key='sol-'+key,
        ....:                           interpolation_key='interp-'+key)
        ....:     graph2D_mercator+=c.plot_integrated(interpolation_key='interp-'+key,
        ....:                               chart=mercator, thickness=2)

    Prepare a grid of Mercator coordinates lines, and plot the curves
    over it::

        sage: graph2D_mercator_coords=mercator.plot(chart=mercator,
        ....:                            number_values=8,color='yellow')
        sage: graph2D_mercator + graph2D_mercator_coords
        Graphics object consisting of 18 graphics primitives

    .. PLOT::

        S2 = Manifold(2, 'S^2', start_index=1)
        polar = S2.chart('th ph'); th, ph = polar[:]
        epolar = polar.frame()
        ch_basis = S2.automorphism_field()
        ch_basis[1,1], ch_basis[2,2] = 1, 1/sin(th)
        epolar_ON = epolar.new_frame(ch_basis, 'epolar_ON')
        nab = S2.affine_connection('nab')
        _ = nab.set_coef(frame=epolar_ON)
        t,tmin,tmax,th0,ph0,v_th0,v_ph0 = var('t tmin tmax th0 ph0 v_th0 v_ph0')
        p = S2.point((th0, ph0), name='p')
        Tp = S2.tangent_space(p)
        v = Tp((v_th0, v_ph0), basis=epolar_ON.at(p))
        c = S2.integrated_autoparallel_curve(nab, (t, tmin, tmax), v,
                                                  chart=polar, name='c')
        dict_params={'latit':{tmin:0,tmax:3,th0:pi/4,ph0:0.1,v_th0:0,v_ph0:1},
                'longi':{tmin:0,tmax:3,th0:0.1,ph0:0.1,v_th0:1,v_ph0:0}}
        mercator = S2.chart(r'xi:(-oo,oo):\xi ze:(0,2*pi):\zeta')
        xi,ze = var('xi ze')
        _ = polar.transition_map(mercator, (log(tan(th/2)), ph))
        identity = S2.identity_map()
        identity.coord_functions(polar, mercator)
        graph2D_mercator = Graphics()
        for key in dict_params:
            sol = c.solve(solution_key='sol-'+key,
                          parameters_values=dict_params[key])
            interp = c.interpolate(solution_key='sol-'+key,
                                   interpolation_key='interp-'+key)
            graph2D_mercator += c.plot_integrated(interpolation_key='interp-'+key,
                                            chart=mercator, thickness=2)
        graph2D_mercator_coords = mercator.plot(chart=mercator,
                                                number_values=8, color='yellow')
        sphinx_plot(graph2D_mercator + graph2D_mercator_coords)

    The resulting curves are horizontal and vertical as expected.
    It is easier to check that these are latitude and longitude lines
    respectively when plotting them on `\mathbb{S}^{2}`.
    To do so, use `\mathbb{R}^{3}` as the codomain of the standard
    map embedding `(\mathbb{S}^{2}, (\theta, \phi))` in the
    3-dimensional Euclidean space::

        sage: R3 = Manifold(3, 'R3', start_index=1)
        sage: cart.<X,Y,Z> = R3.chart()
        sage: euclid_embedding = S2.diff_map(R3,
        ....:  {(polar, cart):[sin(th)*cos(ph),sin(th)*sin(ph),cos(th)]})

    Plot the resulting curves on the grid of polar coordinates lines on
    `\mathbb{S}^{2}`::

        sage: graph3D_embedded_curves = Graphics()
        sage: for key in dict_params:
        ....:     graph3D_embedded_curves += c.plot_integrated(interpolation_key='interp-'+key,
        ....:            mapping=euclid_embedding, thickness=5,
        ....:            display_tangent=True, scale=0.4, width_tangent=0.5)
        sage: graph3D_embedded_polar_coords = polar.plot(chart=cart,
        ....:                          mapping=euclid_embedding,
        ....:                          number_values=15, color='yellow')
        sage: graph3D_embedded_curves + graph3D_embedded_polar_coords
        Graphics3d Object

    .. PLOT::

        S2 = Manifold(2, 'S^2', start_index=1)
        polar = S2.chart('th ph'); th, ph = polar[:]
        epolar = polar.frame()
        ch_basis = S2.automorphism_field()
        ch_basis[1,1], ch_basis[2,2] = 1, 1/sin(th)
        epolar_ON = epolar.new_frame(ch_basis, 'epolar_ON')
        nab = S2.affine_connection('nab')
        _ = nab.set_coef(frame=epolar_ON)
        t,tmin,tmax,th0,ph0,v_th0,v_ph0 = var('t tmin tmax th0 ph0 v_th0 v_ph0')
        p = S2.point((th0, ph0), name='p')
        Tp = S2.tangent_space(p)
        v = Tp((v_th0, v_ph0), basis=epolar_ON.at(p))
        c = S2.integrated_autoparallel_curve(nab, (t, tmin, tmax), v,
                                             chart=polar, name='c')
        dict_params = {'latit':{tmin:0,tmax:3,th0:pi/4,ph0:0.1,v_th0:0,v_ph0:1},
                       'longi':{tmin:0,tmax:3,th0:0.1,ph0:0.1,v_th0:1,v_ph0:0}}
        R3 = Manifold(3, 'R3', start_index=1)
        cart = R3.chart('X Y Z'); X, Y, Z = cart[:]
        euclid_embedding = S2.diff_map(R3,
              {(polar, cart): [sin(th)*cos(ph),sin(th)*sin(ph),cos(th)]})
        graph3D_embedded_curves = Graphics()
        for key in dict_params:
            sol = c.solve(solution_key='sol-'+key,
                          parameters_values=dict_params[key])
            interp = c.interpolate(solution_key='sol-'+key,
                                   interpolation_key='interp-'+key)
            graph3D_embedded_curves += c.plot_integrated(interpolation_key='interp-'+key,
                     mapping=euclid_embedding, thickness=5,
                     display_tangent=True, scale=0.4, width_tangent=0.5)
        graph3D_embedded_polar_coords = polar.plot(chart=cart,
                                       mapping=euclid_embedding,
                                       number_values=15, color='yellow')
        graph = graph3D_embedded_curves+graph3D_embedded_polar_coords
        sphinx_plot(graph)

    Finally, one may plot a general autoparallel curve with respect to
    `\nabla` that is neither a line of latitude or longitude.
    The vectors tangent to such a curve make an angle different from 0
    or `\pi/2` with the lines of latitude and longitude.
    Then, compute a curve such that both components of its initial
    tangent vectors are nonzero::

        sage: sol = c.solve(solution_key='sol-angle',
        ....:  parameters_values={tmin:0,tmax:2,th0:pi/4,ph0:0.1,v_th0:1,v_ph0:8})
        sage: interp = c.interpolate(solution_key='sol-angle',
        ....:                          interpolation_key='interp-angle')

    Plot the resulting curve in the Mercator plane.
    This generates a straight line, as expected::

        sage: c.plot_integrated(interpolation_key='interp-angle',
        ....:         chart=mercator, thickness=1, display_tangent=True,
        ....:         scale=0.2, width_tangent=0.2)
        Graphics object consisting of 11 graphics primitives

    .. PLOT::

        S2 = Manifold(2, 'S^2', start_index=1)
        polar = S2.chart('th ph'); th, ph = polar[:]
        epolar = polar.frame()
        ch_basis = S2.automorphism_field()
        ch_basis[1,1], ch_basis[2,2] = 1, 1/sin(th)
        epolar_ON = epolar.new_frame(ch_basis, 'epolar_ON')
        nab = S2.affine_connection('nab')
        _ = nab.set_coef(frame=epolar_ON)
        t,tmin,tmax,th0,ph0,v_th0,v_ph0 = var('t tmin tmax th0 ph0 v_th0 v_ph0')
        p = S2.point((th0, ph0), name='p')
        Tp = S2.tangent_space(p)
        v = Tp((v_th0, v_ph0), basis=epolar_ON.at(p))
        c = S2.integrated_autoparallel_curve(nab, (t, tmin, tmax), v,
                                                  chart=polar, name='c')
        mercator = S2.chart(r'xi:(-oo,oo):\xi ze:(0,2*pi):\zeta')
        xi, ze = mercator[:]
        trans_map = polar.transition_map(mercator, (log(tan(th/2)), ph))
        identity = S2.identity_map()
        _ = identity.coord_functions(polar, mercator)
        sol = c.solve(solution_key='sol-angle',
            parameters_values={tmin:0,tmax:2,th0:pi/4,ph0:0.1,v_th0:1,v_ph0:8})
        interp = c.interpolate(solution_key='sol-angle',
                               interpolation_key='interp-angle')
        graph2D_mercator_angle_curve=c.plot_integrated(
                      interpolation_key='interp-angle',
                      chart=mercator, thickness=1, display_tangent=True,
                      scale=0.2, width_tangent=0.2)
        sphinx_plot(graph2D_mercator_angle_curve)

    One may eventually plot such a curve on `\mathbb{S}^{2}`::

        sage: graph3D_embedded_angle_curve=c.plot_integrated(interpolation_key='interp-angle',
        ....:        mapping=euclid_embedding, thickness=5,
        ....:        display_tangent=True, scale=0.1, width_tangent=0.5)
        sage: graph3D_embedded_angle_curve + graph3D_embedded_polar_coords
        Graphics3d Object

    .. PLOT::

        S2 = Manifold(2, 'S^2', start_index=1)
        polar = S2.chart('th ph'); th, ph = polar[:]
        epolar = polar.frame()
        ch_basis = S2.automorphism_field()
        ch_basis[1,1], ch_basis[2,2] = 1, 1/sin(th)
        epolar_ON = epolar.new_frame(ch_basis, 'epolar_ON')
        nab = S2.affine_connection('nab')
        _ = nab.set_coef(frame=epolar_ON)
        t,tmin,tmax,th0,ph0,v_th0,v_ph0 = var('t tmin tmax th0 ph0 v_th0 v_ph0')
        p = S2.point((th0, ph0), name='p')
        Tp = S2.tangent_space(p)
        v = Tp((v_th0, v_ph0), basis=epolar_ON.at(p))
        c = S2.integrated_autoparallel_curve(nab, (t, tmin, tmax), v,
                                                  chart=polar, name='c')
        R3 = Manifold(3, 'R3', start_index=1)
        cart = R3.chart('X Y Z')
        euclid_embedding = S2.diff_map(R3,
              {(polar, cart):[sin(th)*cos(ph),sin(th)*sin(ph),cos(th)]})
        sol = c.solve(solution_key='sol-angle',
            parameters_values={tmin:0,tmax:2,th0:pi/4,ph0:0.1,v_th0:1,v_ph0:8})
        interp = c.interpolate(solution_key='sol-angle',
                               interpolation_key='interp-angle')
        graph3D_embedded_angle_curve = c.plot_integrated(interpolation_key='interp-angle',
            mapping=euclid_embedding, thickness=5, display_tangent=True,
            scale=0.1, width_tangent=0.5)
        graph3D_embedded_polar_coords = polar.plot(chart=cart,
             mapping=euclid_embedding, number_values=15, color='yellow')
        graph = graph3D_embedded_angle_curve + graph3D_embedded_polar_coords
        sphinx_plot(graph)

    All the curves presented are loxodromes, and the differential system
    defining them (displayed above) may be solved analytically,
    providing the following expressions:

    .. MATH::

        \begin{aligned}
        \theta(t) &= \theta_{0} + \dot{\theta}_{0} (t - t_{0}),      \\
        \phi(t) &= \phi_{0} - \frac{1}{\tan \alpha} \left(
        \ln \tan \frac{\theta_{0} + \dot{\theta}_{0} (t - t_{0})}{2} -
        \ln \tan \frac{\theta_{0}}{2} \right),
        \end{aligned}

    where `\alpha` is the angle between the curve and any latitude
    line it crosses; then, one finds
    `\tan \alpha = - \dot{\theta}_{0} / (\dot{\phi}_{0} \sin \theta_{0})`
    (then `\tan \alpha \leq 0` when the initial tangent vector
    points towards the southeast).

    In order to use these expressions to compare with the result
    provided by the numerical integration, remember that the components
    ``(v_th0, v_ph0)`` of the initial
    tangent vector ``v`` refer to the basis
    ``epolar_ON`` `= (\hat{e}_{\theta}, \hat{e}_{\phi})` and not the
    coordinate basis
    ``epolar`` `= (e_{\theta}, e_{\phi})`.
    Therefore, the following relations hold:
    ``v_ph0`` `= \dot{\phi}_{0} \sin \theta_{0}` (and not merely
    `\dot{\phi}_{0}`), while ``v_th0`` clearly is `\dot{\theta}_{0}`.

    With this in mind, plot an analytical curve to compare with a
    numerical solution::

        sage: graph2D_mercator_angle_curve=c.plot_integrated(interpolation_key='interp-angle',
        ....:                               chart=mercator, thickness=1)
        sage: expr_ph = ph0+v_ph0/v_th0*(ln(tan((v_th0*t+th0)/2))-ln(tan(th0/2)))
        sage: c_loxo = S2.curve({polar:[th0+v_th0*t, expr_ph]}, (t,0,2),
        ....:                                             name='c_loxo')

    Ask for the expression of the loxodrome in terms of the Mercator
    chart in order to add it to its dictionary of expressions.
    It is a particularly long expression, and there is no particular
    need to display it, which is why it may simply be affected to an
    arbitrary variable ``expr_mercator``, which will never be used
    again.
    But adding the expression to the dictionary is required to plot the
    curve with respect to the Mercator chart::

        sage: expr_mercator = c_loxo.expression(chart2=mercator)

    Plot the curves (for clarity, set a 2 degrees shift in the initial
    value of `\theta_{0}` so that the curves do not overlap)::

        sage: graph2D_mercator_loxo = c_loxo.plot(chart=mercator,
        ....:  parameters={th0:pi/4+2*pi/180, ph0:0.1, v_th0:1, v_ph0:8},
        ....:  thickness=1, color='blue')
        sage: graph2D_mercator_angle_curve + graph2D_mercator_loxo
        Graphics object consisting of 2 graphics primitives

    .. PLOT::

        S2 = Manifold(2, 'S^2', start_index=1)
        polar = S2.chart('th ph'); th, ph = polar[:]
        epolar = polar.frame()
        ch_basis = S2.automorphism_field()
        ch_basis[1,1], ch_basis[2,2] = 1, 1/sin(th)
        epolar_ON = epolar.new_frame(ch_basis, 'epolar_ON')
        nab = S2.affine_connection('nab')
        _ = nab.set_coef(frame=epolar_ON)
        t, tmin, tmax, th0, ph0 = var('t tmin tmax th0 ph0')
        v_th0, v_ph0, alpha = var('v_th0 v_ph0 alpha')
        p = S2.point((th0, ph0), name='p')
        Tp = S2.tangent_space(p)
        v = Tp((v_th0, v_ph0), basis=epolar_ON.at(p))
        c = S2.integrated_autoparallel_curve(nab, (t, tmin, tmax), v,
                                                  chart=polar, name='c')
        mercator = S2.chart(r'xi:(-oo,oo):\xi ze:(0,2*pi):\zeta')
        xi, ze = mercator[:]
        trans_map = polar.transition_map(mercator, (log(tan(th/2)), ph))
        identity = S2.identity_map()
        _ = identity.coord_functions(polar, mercator)
        sol = c.solve(solution_key='sol-angle',
            parameters_values={tmin:0,tmax:2,th0:pi/4,ph0:0.1,v_th0:1,v_ph0:8})
        interp = c.interpolate(solution_key='sol-angle',
                               interpolation_key='interp-angle')
        graph2D_mercator_angle_curve = c.plot_integrated(interpolation_key='interp-angle',
                                                         chart=mercator, thickness=1)
        expr_ph = ph0+v_ph0/v_th0*(ln(tan((v_th0*t+th0)/2))-ln(tan(th0/2)))
        c_loxo = S2.curve({polar: [th0+v_th0*t, expr_ph]}, (t,0,2), name='c_loxo')
        expr = c_loxo.expression(chart2=mercator)
        graph2D_mercator_loxo = c_loxo.plot(chart=mercator,
              parameters={th0:pi/4+2*pi/180, ph0:0.1, v_th0:1, v_ph0:8},
              thickness=1, color='blue')
        sphinx_plot(graph2D_mercator_angle_curve+graph2D_mercator_loxo)

    Both curves do have the same aspect.
    One may eventually compare these curves on `\mathbb{S}^{2}`::

        sage: graph3D_embedded_angle_curve=c.plot_integrated(interpolation_key='interp-angle',
        ....:                     mapping=euclid_embedding, thickness=3)
        sage: graph3D_embedded_loxo = c_loxo.plot(mapping=euclid_embedding,
        ....:  parameters={th0:pi/4+2*pi/180, ph0:0.1, v_th0:1, v_ph0:8},
        ....:  thickness=3, color = 'blue')
        sage: (graph3D_embedded_angle_curve + graph3D_embedded_loxo
        ....:  + graph3D_embedded_polar_coords)
        Graphics3d Object

    .. PLOT::

        S2 = Manifold(2, 'S^2', start_index=1)
        polar = S2.chart('th ph'); th, ph = polar[:]
        epolar = polar.frame()
        ch_basis = S2.automorphism_field()
        ch_basis[1,1], ch_basis[2,2] = 1, 1/sin(th)
        epolar_ON = epolar.new_frame(ch_basis, 'epolar_ON')
        nab = S2.affine_connection('nab')
        _ = nab.set_coef(frame=epolar_ON)
        t, tmin, tmax, th0, ph0 = var('t tmin tmax th0 ph0')
        v_th0, v_ph0, alpha = var('v_th0 v_ph0 alpha')
        p = S2.point((th0, ph0), name='p')
        Tp = S2.tangent_space(p)
        v = Tp((v_th0, v_ph0), basis=epolar_ON.at(p))
        c = S2.integrated_autoparallel_curve(nab, (t, tmin, tmax), v,
                                             chart=polar, name='c')
        R3 = Manifold(3, 'R3', start_index=1)
        cart = R3.chart('X Y Z')
        euclid_embedding = S2.diff_map(R3,
              {(polar, cart):[sin(th)*cos(ph),sin(th)*sin(ph),cos(th)]})
        sol = c.solve(solution_key='sol-angle',
            parameters_values={tmin:0,tmax:2,th0:pi/4,ph0:0.1,v_th0:1,v_ph0:8})
        interp = c.interpolate(solution_key='sol-angle',
                               interpolation_key='interp-angle')
        graph3D_embedded_angle_curve = c.plot_integrated(interpolation_key='interp-angle',
                                  mapping=euclid_embedding, thickness=3)
        expr_ph = ph0+v_ph0/v_th0*(ln(tan((v_th0*t+th0)/2))-ln(tan(th0/2)))
        c_loxo = S2.curve({polar: [th0+v_th0*t, expr_ph]}, (t,0,2), name='c_loxo')
        graph3D_embedded_loxo = c_loxo.plot(mapping=euclid_embedding,
              parameters={th0:pi/4+2*pi/180, ph0:0.1, v_th0:1, v_ph0:8},
              thickness=3, color='blue')
        graph3D_embedded_polar_coords = polar.plot(chart=cart,
             mapping=euclid_embedding, number_values=15, color='yellow')
        graph = graph3D_embedded_angle_curve + graph3D_embedded_loxo
        graph += graph3D_embedded_polar_coords
        sphinx_plot(graph)
    """

    def __init__(self, parent, affine_connection, curve_parameter,
                 initial_tangent_vector, chart=None, name=None,
                 latex_name=None, verbose=False, across_charts=False):
        r"""
        Construct an autoparallel curve with respect to the given affine
        connection with the given initial tangent vector.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: [t, A, B] = var('t A B')
            sage: nab = M.affine_connection('nabla', r'\nabla')
            sage: nab[X.frame(),0,0,1],nab[X.frame(),2,1,2]=A*x1^2,B*x2*x3
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_autoparallel_curve(nab, (t, 0, 5), v,
            ....:                                          name='c') ; c
            Integrated autoparallel curve c in the 3-dimensional
             differentiable manifold M
            sage: TestSuite(c).run()
        """

        # setting the chart to gain access to the coordinate functions
        if chart is None:
            chart = parent.codomain().default_chart()

        velocities = chart.symbolic_velocities()

        dim = parent.codomain().dim()
        i0 = parent.codomain().start_index()

        self._across_charts = across_charts
        if not across_charts:

            equations_rhs = []

            gamma = affine_connection.coef(frame=chart.frame())

            for rho in range(dim):
                rhs = 0
                for mu in range(dim):
                    for nu in range(dim):
                        vMUvNU = velocities[mu] * velocities[nu]
                        gammaRHO_mu_nu = gamma[[rho+i0, mu+i0, nu+i0]].expr(chart=chart)
                        # line above is the expression of the scalar
                        # field 'gamma[[rho+i0, mu+i0, nu+i0]]' in terms
                        # of 'chart' (here, in any point of the manifold,
                        # the scalar field 'gamma[[rho+i0, mu+i0, nu+i0]]'
                        # provides the coefficient [rho+i0, mu+i0, nu+i0]
                        # of the affine connection with respect to frame
                        # 'chart.frame()')
                        rhs -= gammaRHO_mu_nu * vMUvNU
                        # 'vMUvNU' and 'gammaRHO_mu_nu' only used for the
                        # line above to be shorter
                equations_rhs += [rhs.simplify_full()]
        else:
            equations_rhs = {}          # Dict of all equation in all top_charts
            for chart in parent.codomain().top_charts():
                velocities = chart.symbolic_velocities()
                equations_rhs_chart = []  # Equation in one chart
                gamma = affine_connection.coef(frame=chart.frame())
                for rho in range(dim):
                    rhs = 0
                    for mu in range(dim):
                        for nu in range(dim):
                            vMUvNU = velocities[mu] * velocities[nu]
                            gammaRHO_mu_nu = gamma[
                                [rho + i0, mu + i0, nu + i0]].expr(chart=chart)
                            rhs -= gammaRHO_mu_nu * vMUvNU
                    equations_rhs_chart += [rhs.simplify_full()]
                equations_rhs[chart] = equations_rhs_chart

        IntegratedCurve.__init__(self, parent, equations_rhs,
                                 velocities, curve_parameter,
                                 initial_tangent_vector, chart=chart,
                                 name=name, latex_name=latex_name,
                                 verbose=verbose, across_charts=across_charts)

        self._affine_connection = affine_connection

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: [t, A, B] = var('t A B')
            sage: nab = M.affine_connection('nabla', r'\nabla')
            sage: nab[X.frame(),0,0,1],nab[X.frame(),2,1,2]=A*x1^2,B*x2*x3
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_autoparallel_curve(nab, (t,0,5), v); c
            Integrated autoparallel curve in the 3-dimensional
             differentiable manifold M
            sage: c = M.integrated_autoparallel_curve(nab, (t, 0, 5), v,
            ....:                                          name='c') ; c
            Integrated autoparallel curve c in the 3-dimensional
             differentiable manifold M
        """

        description = "Integrated autoparallel curve "
        if self._name is not None:
            description += self._name + " "
        description += "in the {}".format(self._codomain)
        return description

    def __reduce__(self):
        r"""
        Reduction function for the pickle protocole.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: [t, A, B] = var('t A B')
            sage: nab = M.affine_connection('nabla', r'\nabla')
            sage: nab[X.frame(),0,0,1],nab[X.frame(),2,1,2]=A*x1^2,B*x2*x3
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_autoparallel_curve(nab, (t, 0, 5), v,
            ....:                                              name='c')
            sage: c.__reduce__()
            (<class 'sage.manifolds.differentiable.manifold_homset.IntegratedAutoparallelCurveSet_with_category.element_class'>,
             (Set of Morphisms from Real interval (0, 5) to
              3-dimensional differentiable manifold M in Category of homsets of
              topological spaces which actually are integrated autoparallel
              curves with respect to a certain affine connection,
              Affine connection nabla on the 3-dimensional
              differentiable manifold M,
              t,
              Tangent vector at Point p on the 3-dimensional
               differentiable manifold M,
              Chart (M, (x1, x2, x3)),
              'c',
              'c',
              False,
              False))

        Test of pickling::

            sage: loads(dumps(c))
            Integrated autoparallel curve c in the 3-dimensional differentiable manifold M
        """

        return (type(self), (self.parent(), self._affine_connection,
                self._curve_parameter, self._initial_tangent_vector,
                self._chart, self._name, self._latex_name, False,
                self._across_charts))

    def system(self, verbose=False):
        r"""
        Provide a detailed description of the system defining the
        autoparallel curve and returns the system defining it: chart,
        equations and initial conditions.

        INPUT:

        - ``verbose`` -- boolean (default: ``False``); prints a detailed
          description of the curve

        OUTPUT:

        - list containing the

          * the equations
          * the initial conditions
          * the chart

        EXAMPLES:

        System defining an autoparallel curve::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: [t, A, B] = var('t A B')
            sage: nab = M.affine_connection('nabla', r'\nabla')
            sage: nab[X.frame(),0,0,1],nab[X.frame(),2,1,2]=A*x1^2,B*x2*x3
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_autoparallel_curve(nab, (t, 0, 5), v)
            sage: sys = c.system(verbose=True)
            Autoparallel curve in the 3-dimensional differentiable
             manifold M equipped with Affine connection nabla on the
             3-dimensional differentiable manifold M, and integrated
             over the Real interval (0, 5) as a solution to the
             following equations, written with respect to
             Chart (M, (x1, x2, x3)):
            <BLANKLINE>
            Initial point: Point p on the 3-dimensional differentiable
             manifold M with coordinates [0, 0, 0] with respect to
             Chart (M, (x1, x2, x3))
            Initial tangent vector: Tangent vector at Point p on the
             3-dimensional differentiable manifold M with
             components [1, 0, 1] with respect to Chart (M, (x1, x2, x3))
            <BLANKLINE>
            d(x1)/dt = Dx1
            d(x2)/dt = Dx2
            d(x3)/dt = Dx3
            d(Dx1)/dt = -A*Dx1*Dx2*x1^2
            d(Dx2)/dt = 0
            d(Dx3)/dt = -B*Dx2*Dx3*x2*x3
            <BLANKLINE>
            sage: sys_bis = c.system()
            sage: sys_bis == sys
            True
        """

        v0 = self._initial_tangent_vector
        chart = self._chart

        if verbose:
            initial_tgt_space = v0.parent()
            initial_pt = initial_tgt_space.base_point() # retrieves
            # the initial point as the base point of the tangent space
            # to which initial tangent vector belongs
            initial_pt_coords = list(initial_pt.coordinates(chart))
            # previous line converts to list since would otherwise be a
            # tuple ; will raise error if coordinates in chart are not
            # known

            initial_coord_basis = chart.frame().at(initial_pt)
            initial_tgt_vec_comps = v0[initial_coord_basis,:] # will
            # raise error if components in coordinate basis are not
            # known

            description = "Autoparallel curve "
            if self._name is not None:
                description += self._name + " "
            description += "in the {} ".format(self.codomain())
            description += "equipped with "
            description += "{}, ".format(self._affine_connection)
            description += "and integrated over the "
            description += "{} ".format(self.domain())
            description += "as a solution to the following equations, "
            description += "written with respect to "
            description += "{}:\n\n".format(chart)

            description += "Initial point: {} ".format(initial_pt)
            description += "with coordinates "
            description += "{} ".format(initial_pt_coords)
            description += "with respect to {}\n".format(chart)

            description += "Initial tangent vector: {} ".format(v0)
            description += "with components "
            description += "{}".format(initial_tgt_vec_comps)
            description += " with respect to {}\n\n".format(chart)

            for coord_func,velocity in zip(chart[:],self._velocities):
                description += "d({})/d{} = {}\n".format(coord_func,
                                                  self._curve_parameter,
                                                  velocity)

            for velocity,eqn in zip(self._velocities,self._equations_rhs):
                description += "d({})/d{} = {}\n".format(velocity,
                                                  self._curve_parameter,
                                                  eqn)

            print(description)

        return [self._equations_rhs, v0, chart]

class IntegratedGeodesic(IntegratedAutoparallelCurve):
    r"""
    Geodesic on the manifold with respect to a given metric.

    INPUT:

    - ``parent`` --
      :class:`~sage.manifolds.differentiable.manifold_homset.IntegratedGeodesicSet`
      the set of curves `\mathrm{Hom_{geodesic}}(I, M)` to which the
      curve belongs
    - ``metric`` --
      :class:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric`
      metric with respect to which the curve is a geodesic
    - ``curve_parameter`` -- symbolic expression to be used as the
      parameter of the curve (the equations defining an instance of
      IntegratedGeodesic are such that ``t`` will actually be an affine
      parameter of the curve);
    - ``initial_tangent_vector`` --
      :class:`~sage.manifolds.differentiable.tangent_vector.TangentVector`
      initial tangent vector of the curve
    - ``chart`` -- (default: ``None``) chart on the manifold in terms of
      which the equations are expressed; if ``None`` the default chart
      of the manifold is assumed
    - ``name`` -- (default: ``None``) string; symbol given to the curve
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote
      the curve; if none is provided, ``name`` will be used

    EXAMPLES:

    Geodesics of the unit 2-sphere `\mathbb{S}^{2}`.
    Start with declaring the standard polar coordinates
    `(\theta, \phi)` on `\mathbb{S}^{2}` and the
    corresponding coordinate frame `(e_{\theta}, e_{\phi})`::

        sage: S2 = Manifold(2, 'S^2', structure='Riemannian', start_index=1)
        sage: polar.<th,ph>=S2.chart('th ph')
        sage: epolar = polar.frame()

    Set the standard round metric::

        sage: g = S2.metric()
        sage: g[1,1], g[2,2] = 1, (sin(th))^2

    Set generic initial conditions for the geodesics to compute::

        sage: [th0, ph0, v_th0, v_ph0] = var('th0 ph0 v_th0 v_ph0')
        sage: p = S2.point((th0, ph0), name='p')
        sage: Tp = S2.tangent_space(p)
        sage: v = Tp((v_th0, v_ph0), basis=epolar.at(p))

    Declare the corresponding integrated geodesic and display the
    differential system it satisfies::

        sage: [t, tmin, tmax] = var('t tmin tmax')
        sage: c = S2.integrated_geodesic(g, (t, tmin, tmax), v,
        ....:                            chart=polar, name='c')
        sage: sys = c.system(verbose=True)
        Geodesic c in the 2-dimensional Riemannian manifold S^2
         equipped with Riemannian metric g on the 2-dimensional
         Riemannian manifold S^2, and integrated over the Real
         interval (tmin, tmax) as a solution to the following geodesic
         equations, written with respect to Chart (S^2, (th, ph)):
        <BLANKLINE>
        Initial point: Point p on the 2-dimensional Riemannian
        manifold S^2 with coordinates [th0, ph0] with respect to
        Chart (S^2, (th, ph))
        Initial tangent vector: Tangent vector at Point p on the
        2-dimensional Riemannian manifold S^2 with
        components [v_th0, v_ph0] with respect to Chart (S^2, (th, ph))
        <BLANKLINE>
        d(th)/dt = Dth
        d(ph)/dt = Dph
        d(Dth)/dt = Dph^2*cos(th)*sin(th)
        d(Dph)/dt = -2*Dph*Dth*cos(th)/sin(th)
        <BLANKLINE>

    Set a dictionary providing the parameter range and the initial
    conditions for various geodesics::

        sage: dict_params={'equat':{tmin:0,tmax:3,th0:pi/2,ph0:0.1,v_th0:0,v_ph0:1},
        ....:   'longi':{tmin:0,tmax:3,th0:0.1,ph0:0.1,v_th0:1,v_ph0:0},
        ....:   'angle':{tmin:0,tmax:3,th0:pi/4,ph0:0.1,v_th0:1,v_ph0:1}}

    Use `\mathbb{R}^{3}` as the codomain of the standard map
    embedding `(\mathbb{S}^{2}, (\theta, \phi))` in the
    3-dimensional Euclidean space::

        sage: R3 = Manifold(3, 'R3', start_index=1)
        sage: cart.<X,Y,Z> = R3.chart()
        sage: euclid_embedding = S2.diff_map(R3,
        ....:  {(polar, cart):[sin(th)*cos(ph),sin(th)*sin(ph),cos(th)]})

    Solve, interpolate and prepare the plot for the solutions
    corresponding to the three initial conditions previously set::

        sage: # needs scipy sage.plot
        sage: graph3D_embedded_geods = Graphics()
        sage: for key in dict_params:
        ....:     sol = c.solve(solution_key='sol-'+key,
        ....:                        parameters_values=dict_params[key])
        ....:     interp = c.interpolate(solution_key='sol-'+key,
        ....:                           interpolation_key='interp-'+key)
        ....:     graph3D_embedded_geods += c.plot_integrated(interpolation_key='interp-'+key,
        ....:                      mapping=euclid_embedding, thickness=5,
        ....:                      display_tangent=True, scale=0.3,
        ....:                      width_tangent=0.5)

    Plot the resulting geodesics on the grid of polar coordinates lines
    on `\mathbb{S}^{2}` and check that these are great circles::

        sage: # needs scipy sage.plot
        sage: graph3D_embedded_polar_coords = polar.plot(chart=cart,
        ....:                          mapping=euclid_embedding,
        ....:                          number_values=15, color='yellow')
        sage: graph3D_embedded_geods + graph3D_embedded_polar_coords
        Graphics3d Object

    .. PLOT::

        S2 = Manifold(2, 'S^2', structure='Riemannian', start_index=1)
        polar = S2.chart('th ph'); th, ph = polar[:]
        epolar = polar.frame()
        g = S2.metric()
        g[1,1], g[2,2] = 1, (sin(th))**2
        t,tmin,tmax,th0,ph0,v_th0,v_ph0 = var('t tmin tmax th0 ph0 v_th0 v_ph0')
        p = S2.point((th0, ph0), name='p')
        Tp = S2.tangent_space(p)
        v = Tp((v_th0, v_ph0), basis=epolar.at(p))
        c = S2.integrated_geodesic(g, (t, tmin, tmax), v, chart=polar,
                                                               name='c')
        dict_params={'equat':{tmin:0,tmax:3,th0:pi/2,ph0:0.1,v_th0:0,v_ph0:1},
                     'longi':{tmin:0,tmax:3,th0:0.1,ph0:0.1,v_th0:1,v_ph0:0},
                     'angle':{tmin:0,tmax:3,th0:pi/4,ph0:0.1,v_th0:1,v_ph0:1}}
        R3 = Manifold(3, 'R3', start_index=1)
        cart = R3.chart('X Y Z')
        euclid_embedding = S2.diff_map(R3,
              {(polar, cart): [sin(th)*cos(ph),sin(th)*sin(ph),cos(th)]})
        graph3D_embedded_geods = Graphics()
        for key in dict_params:
            sol = c.solve(solution_key='sol-'+key,
                          parameters_values=dict_params[key])
            interp = c.interpolate(solution_key='sol-'+key,
                                   interpolation_key='interp-'+key)
            graph3D_embedded_geods += c.plot_integrated(interpolation_key='interp-'+key,
                                  mapping=euclid_embedding, thickness=5,
                                  display_tangent=True, scale=0.3,
                                  width_tangent=0.5)
        graph3D_embedded_polar_coords = polar.plot(chart=cart,
                                       mapping=euclid_embedding,
                                       number_values=15, color='yellow')
        graph = graph3D_embedded_geods + graph3D_embedded_polar_coords
        sphinx_plot(graph)
    """

    def __init__(self, parent, metric, curve_parameter,
                 initial_tangent_vector, chart=None, name=None,
                 latex_name=None, verbose=False, across_charts=False):
        r"""
        Construct a geodesic curve with respect to the given metric with the
        given initial tangent vector.

        TESTS::

            sage: S2 = Manifold(2, 'S^2', structure='Riemannian')
            sage: X.<theta,phi> = S2.chart()
            sage: t, A = var('t A')
            sage: g = S2.metric()
            sage: g[0,0] = A
            sage: g[1,1] = A*sin(theta)^2
            sage: p = S2.point((pi/2,0), name='p')
            sage: Tp = S2.tangent_space(p)
            sage: v = Tp((1/sqrt(2),1/sqrt(2)))
            sage: c = S2.integrated_geodesic(g, (t,0,pi), v, name='c'); c
            Integrated geodesic c in the 2-dimensional Riemannian
             manifold S^2
            sage: TestSuite(c).run()
        """

        affine_connection = metric.connection()

        IntegratedAutoparallelCurve.__init__(self, parent,
                                             affine_connection, curve_parameter,
                                             initial_tangent_vector, chart=chart,
                                             name=name, latex_name=latex_name,
                                             verbose=verbose, across_charts=across_charts)

        self._metric = metric
        self._across_charts = across_charts

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: S2 = Manifold(2, 'S^2', structure='Riemannian')
            sage: X.<theta,phi> = S2.chart()
            sage: t, A = var('t A')
            sage: g = S2.metric()
            sage: g[0,0] = A
            sage: g[1,1] = A*sin(theta)^2
            sage: p = S2.point((pi/2,0), name='p')
            sage: Tp = S2.tangent_space(p)
            sage: v = Tp((1/sqrt(2),1/sqrt(2)))
            sage: c = S2.integrated_geodesic(g, (t, 0, pi), v) ; c
            Integrated geodesic in the 2-dimensional Riemannian
             manifold S^2
            sage: c = S2.integrated_geodesic(g, (t,0,pi), v, name='c'); c
            Integrated geodesic c in the 2-dimensional Riemannian
             manifold S^2
        """

        description = "Integrated geodesic "
        if self._name is not None:
            description += self._name + " "
        description += "in the {}".format(self._codomain)
        return description

    def __reduce__(self):
        r"""
        Reduction function for the pickle protocole.

        TESTS::

            sage: S2 = Manifold(2, 'S^2', structure='Riemannian')
            sage: X.<theta,phi> = S2.chart()
            sage: t, A = var('t A')
            sage: g = S2.metric()
            sage: g[0,0] = A
            sage: g[1,1] = A*sin(theta)^2
            sage: p = S2.point((pi/2,0), name='p')
            sage: Tp = S2.tangent_space(p)
            sage: v = Tp((1/sqrt(2),1/sqrt(2)))
            sage: c = S2.integrated_geodesic(g, (t, 0, pi), v, name='c')
            sage: c.__reduce__()
            (<...IntegratedGeodesicSet_with_category.element_class'>,
             (Set of Morphisms from Real interval (0, pi) to
              2-dimensional Riemannian manifold S^2 in Category of homsets of
              topological spaces which actually are integrated geodesics with
              respect to a certain metric,
              Riemannian metric g on the 2-dimensional Riemannian
              manifold S^2,
              t,
              Tangent vector at Point p on the 2-dimensional
               Riemannian manifold S^2,
              Chart (S^2, (theta, phi)),
              'c',
              'c',
              False,
              False))

        Test of pickling::

            sage: loads(dumps(c))
            Integrated geodesic c in the 2-dimensional Riemannian manifold S^2
        """

        return (type(self), (self.parent(), self._metric,
                self._curve_parameter, self._initial_tangent_vector,
                self._chart, self._name, self._latex_name, False,
                self._across_charts))

    def system(self, verbose=False):
        r"""
        Return the system defining the geodesic: chart, equations and
        initial conditions.

        INPUT:

        - ``verbose`` -- boolean (default: ``False``); prints a detailed
          description of the curve

        OUTPUT: list containing

          * the equations
          * the initial equations
          * the chart

        EXAMPLES:

        System defining a geodesic::

            sage: S2 = Manifold(2, 'S^2',structure='Riemannian')
            sage: X.<theta,phi> = S2.chart()
            sage: t, A = var('t A')
            sage: g = S2.metric()
            sage: g[0,0] = A
            sage: g[1,1] = A*sin(theta)^2
            sage: p = S2.point((pi/2,0), name='p')
            sage: Tp = S2.tangent_space(p)
            sage: v = Tp((1/sqrt(2),1/sqrt(2)))
            sage: c = S2.integrated_geodesic(g, (t, 0, pi), v, name='c')
            sage: sys = c.system(verbose=True)
            Geodesic c in the 2-dimensional Riemannian manifold S^2
             equipped with Riemannian metric g on the 2-dimensional
             Riemannian manifold S^2, and integrated over the Real
             interval (0, pi) as a solution to the following geodesic
             equations, written with respect to Chart (S^2, (theta, phi)):
            <BLANKLINE>
            Initial point: Point p on the 2-dimensional Riemannian
             manifold S^2 with coordinates [1/2*pi, 0] with respect to
             Chart (S^2, (theta, phi))
            Initial tangent vector: Tangent vector at Point p on the
             2-dimensional Riemannian manifold S^2 with
             components [1/2*sqrt(2), 1/2*sqrt(2)] with respect to
             Chart (S^2, (theta, phi))
            <BLANKLINE>
            d(theta)/dt = Dtheta
            d(phi)/dt = Dphi
            d(Dtheta)/dt = Dphi^2*cos(theta)*sin(theta)
            d(Dphi)/dt = -2*Dphi*Dtheta*cos(theta)/sin(theta)
            <BLANKLINE>
            sage: sys_bis = c.system()
            sage: sys_bis == sys
            True
        """

        v0 = self._initial_tangent_vector
        chart = self._chart

        if verbose:
            initial_tgt_space = v0.parent()
            initial_pt = initial_tgt_space.base_point()#retrieves
            # the initial point as the base point of the tangent space
            # to which initial tangent vector belongs
            initial_pt_coords = list(initial_pt.coordinates(chart))
            # previous line converts to list since would otherwise be a
            # tuple ; will raise error if coordinates in chart are not
            # known

            initial_coord_basis = chart.frame().at(initial_pt)
            initial_tgt_vec_comps = v0[initial_coord_basis,:] # will
            # raise error if components in coordinate basis are not
            # known

            description = "Geodesic "
            if self._name is not None:
                description += self._name + " "
            description += "in the {} ".format(self.codomain())
            description += "equipped with "
            description += "{}, ".format(self._metric)
            description += "and integrated over the "
            description += "{} ".format(self.domain())
            description += "as a solution to the following "
            description += "geodesic equations, written with respect to "
            description += "{}:\n\n".format(chart)

            description += "Initial point: {} ".format(initial_pt)
            description += "with coordinates "
            description += "{} ".format(initial_pt_coords)
            description += "with respect to {}\n".format(chart)

            description += "Initial tangent vector: {} ".format(v0)
            description += "with components "
            description += "{}".format(initial_tgt_vec_comps)
            description += " with respect to {}\n\n".format(chart)

            for coord_func,velocity in zip(chart[:],self._velocities):
                description += "d({})/d{} = {}\n".format(coord_func,
                                                  self._curve_parameter,
                                                  velocity)

            for velocity,eqn in zip(self._velocities,self._equations_rhs):
                description += "d({})/d{} = {}\n".format(velocity,
                                                  self._curve_parameter,
                                                  eqn)

            print(description)

        return [self._equations_rhs, v0, chart]
