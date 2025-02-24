r"""
Curves in Manifolds

Given a differentiable manifold `M`, a *differentiable curve* in
`M` is a differentiable mapping

.. MATH::

    \gamma: I \longrightarrow M,

where `I` is an interval of `\RR`.

Differentiable curves are implemented by :class:`DifferentiableCurve`.

AUTHORS:

- Eric Gourgoulhon (2015): initial version
- Travis Scrimshaw (2016): review tweaks

REFERENCES:

- Chap. 1 of [KN1963]_
- Chap. 3 of [Lee2013]_
"""

#*****************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.latex import latex
from sage.misc.decorators import options
from sage.manifolds.point import ManifoldPoint
from sage.manifolds.differentiable.diff_map import DiffMap


class DifferentiableCurve(DiffMap):
    r"""
    Curve in a differentiable manifold.

    Given a differentiable manifold `M`, a *differentiable curve* in
    `M` is a differentiable map

    .. MATH::

        \gamma: I \longrightarrow M,

    where `I` is an interval of `\RR`.

    INPUT:

    - ``parent`` --
      :class:`~sage.manifolds.differentiable.manifold_homset.DifferentiableCurveSet`
      the set of curves `\mathrm{Hom}(I, M)` to which the curve belongs
    - ``coord_expression`` -- (default: ``None``) dictionary (possibly empty)
      of the functions of the curve parameter `t` expressing the curve in
      various charts of `M`, the keys of the dictionary being the charts and
      the values being lists or tuples of `n` symbolic expressions of `t`,
      where `n` is the dimension of `M`
    - ``name`` -- (default: ``None``) string; symbol given to the curve
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote
      the curve; if none is provided, ``name`` will be used
    - ``is_isomorphism`` -- boolean (default: ``False``); determines whether the
      constructed object is a diffeomorphism; if set to ``True``,
      then `M` must have dimension one
    - ``is_identity`` -- boolean (default: ``False``); determines whether the
      constructed object is the identity map; if set to ``True``,
      then `M` must be the interval `I`

    EXAMPLES:

    The lemniscate of Gerono in the 2-dimensional Euclidean plane::

        sage: M = Manifold(2, 'M')
        sage: X.<x,y> = M.chart()
        sage: t = var('t')
        sage: c = M.curve({X: [sin(t), sin(2*t)/2]}, (t, 0, 2*pi), name='c') ; c
        Curve c in the 2-dimensional differentiable manifold M
        sage: type(c)
        <class 'sage.manifolds.differentiable.manifold_homset.DifferentiableCurveSet_with_category.element_class'>

    Instead of declaring the parameter `t`  as a symbolic variable by means
    of ``var('t')``, it is equivalent to get it as the canonical coordinate
    of the real number line (see
    :class:`~sage.manifolds.differentiable.examples.real_line.RealLine`)::

        sage: R.<t> = manifolds.RealLine()
        sage: c = M.curve({X: [sin(t), sin(2*t)/2]}, (t, 0, 2*pi), name='c') ; c
        Curve c in the 2-dimensional differentiable manifold M

    A graphical view of the curve is provided by the method :meth:`plot`::

        sage: c.plot(aspect_ratio=1)                                                    # needs sage.plot
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        M = Manifold(2, 'M')
        X = M.chart('x y')
        t = manifolds.RealLine().canonical_coordinate()
        c = M.curve([sin(t), sin(2*t)/2], (t, 0, 2*pi), name='c')
        g = c.plot(aspect_ratio=1)
        sphinx_plot(g)

    Curves are considered as (manifold) morphisms from real intervals to
    differentiable manifolds::

        sage: c.parent()
        Set of Morphisms from Real interval (0, 2*pi) to 2-dimensional
         differentiable manifold M in Category of smooth manifolds over Real
         Field with 53 bits of precision
        sage: I = R.open_interval(0, 2*pi)
        sage: c.parent() is Hom(I, M)
        True
        sage: c.domain()
        Real interval (0, 2*pi)
        sage: c.domain() is I
        True
        sage: c.codomain()
        2-dimensional differentiable manifold M

    Accordingly, all methods of
    :class:`~sage.manifolds.differentiable.diff_map.DiffMap` are available
    for them. In particular, the method
    :meth:`~sage.manifolds.continuous_map.ContinuousMap.display`
    shows the coordinate representations in various charts of manifold ``M``::

        sage: c.display()
        c: (0, 2*pi) → M
           t ↦ (x, y) = (sin(t), 1/2*sin(2*t))

    Another map method is using the usual call syntax, which returns
    the image of a point in the curve's domain::

        sage: t0 = pi/2
        sage: I(t0)
        Point on the Real number line ℝ
        sage: c(I(t0))
        Point on the 2-dimensional differentiable manifold M
        sage: c(I(t0)).coord(X)
        (1, 0)

    For curves, the value of the parameter, instead of the corresponding
    point in the real line manifold, can be passed directly::

        sage: c(t0)
        Point c(1/2*pi) on the 2-dimensional differentiable manifold M
        sage: c(t0).coord(X)
        (1, 0)
        sage: c(t0) == c(I(t0))
        True

    Instead of a dictionary of coordinate expressions, the curve can be
    defined by a single coordinate expression in a given chart::

        sage: c = M.curve([sin(t), sin(2*t)/2], (t, 0, 2*pi), chart=X, name='c') ; c
        Curve c in the 2-dimensional differentiable manifold M
        sage: c.display()
        c: (0, 2*pi) → M
           t ↦ (x, y) = (sin(t), 1/2*sin(2*t))

    Since ``X`` is the default chart on ``M``, it can be omitted::

        sage: c = M.curve([sin(t), sin(2*t)/2], (t, 0, 2*pi), name='c') ; c
        Curve c in the 2-dimensional differentiable manifold M
        sage: c.display()
        c: (0, 2*pi) → M
           t ↦ (x, y) = (sin(t), 1/2*sin(2*t))

    Note that a curve in `M` can also be created as a differentiable
    map `I \to M`::

        sage: c1 = I.diff_map(M, coord_functions={X: [sin(t), sin(2*t)/2]},
        ....:                 name='c') ; c1
        Curve c in the 2-dimensional differentiable manifold M
        sage: c1.parent() is c.parent()
        True
        sage: c1 == c
        True

    LaTeX symbols representing a curve::

        sage: c = M.curve([sin(t), sin(2*t)/2], (t, 0, 2*pi))
        sage: latex(c)
        \text{Curve in the 2-dimensional differentiable manifold M}
        sage: c = M.curve([sin(t), sin(2*t)/2], (t, 0, 2*pi), name='c')
        sage: latex(c)
        c
        sage: c = M.curve([sin(t), sin(2*t)/2], (t, 0, 2*pi), name='c',
        ....:             latex_name=r'\gamma')
        sage: latex(c)
        \gamma

    The curve's tangent vector field (velocity vector)::

        sage: v = c.tangent_vector_field() ; v
        Vector field c' along the Real interval (0, 2*pi) with values on the
         2-dimensional differentiable manifold M
        sage: v.display()
        c' = cos(t) ∂/∂x + (2*cos(t)^2 - 1) ∂/∂y

    Plot of the curve and its tangent vector field::

        sage: show(c.plot(thickness=2, aspect_ratio=1) +
        ....:      v.plot(chart=X, number_values=17, scale=0.5))

    .. PLOT::

        M = Manifold(2, 'M')
        X = M.chart('x y')
        t = manifolds.RealLine().canonical_coordinate()
        c = M.curve([sin(t), sin(2*t)/2], (t, 0, 2*pi), name='c')
        v = c.tangent_vector_field()
        g = c.plot(thickness=2, aspect_ratio=1) + v.plot(chart=X, number_values=17, scale=0.5)
        sphinx_plot(g)

    Value of the tangent vector field at `t = \pi`::

        sage: v.at(R(pi))
        Tangent vector c' at Point on the 2-dimensional differentiable
         manifold M
        sage: v.at(R(pi)) in M.tangent_space(c(R(pi)))
        True
        sage: v.at(R(pi)).display()
        c' = -∂/∂x + ∂/∂y

    Curves `\RR \to \RR` can be composed: the operator `\circ` is
    given by ``*``::

        sage: f = R.curve(t^2, (t,-oo,+oo))
        sage: g = R.curve(cos(t), (t,-oo,+oo))
        sage: s = g*f ; s
        Differentiable map from the Real number line ℝ to itself
        sage: s.display()
        ℝ → ℝ
           t ↦ cos(t^2)
        sage: s = f*g ; s
        Differentiable map from the Real number line ℝ to itself
        sage: s.display()
        ℝ → ℝ
           t ↦ cos(t)^2

    .. RUBRIC:: Curvature and torsion of a curve in a Riemannian manifold

    Let us consider a helix `C` in the Euclidean space `\mathbb{E}^3`
    parametrized by its arc length `s`::

        sage: E.<x,y,z> = EuclideanSpace()
        sage: R.<s> = manifolds.RealLine()
        sage: C = E.curve((2*cos(s/3), 2*sin(s/3), sqrt(5)*s/3), (s, 0, 6*pi),
        ....:             name='C')

    Its tangent vector field is::

        sage: T = C.tangent_vector_field()
        sage: T.display()
        C' = -2/3*sin(1/3*s) e_x + 2/3*cos(1/3*s) e_y + 1/3*sqrt(5) e_z

    Since `C` is parametrized by its arc length `s`, `T` is a unit vector (with
    respect to the Euclidean metric of `\mathbb{E}^3`)::

        sage: norm(T)
        Scalar field |C'| on the Real interval (0, 6*pi)
        sage: norm(T).display()
        |C'|: (0, 6*pi) → ℝ
           s ↦ 1

    Vector fields along `C` are defined by the method
    :meth:`~sage.manifolds.differentiable.manifold.DifferentiableManifold.vector_field`
    of the domain of `C` with the keyword argument ``dest_map`` set to `C`. For
    instance the derivative vector `T'=\mathrm{d}T/\mathrm{d}s` is::

        sage: I = C.domain(); I
        Real interval (0, 6*pi)
        sage: Tp = I.vector_field([diff(T[i], s) for i in E.irange()], dest_map=C,
        ....:                     name="T'")
        sage: Tp.display()
        T' = -2/9*cos(1/3*s) e_x - 2/9*sin(1/3*s) e_y

    The norm of `T'` is the curvature of the helix::

        sage: kappa = norm(Tp)
        sage: kappa
        Scalar field |T'| on the Real interval (0, 6*pi)
        sage: kappa.expr()
        2/9

    The unit normal vector along `C` is::

        sage: N = Tp / kappa
        sage: N.display()
        -cos(1/3*s) e_x - sin(1/3*s) e_y

    while the binormal vector along `C` is `B = T \times N`::

        sage: B = T.cross_product(N)
        sage: B
        Vector field along the Real interval (0, 6*pi) with values on the
         Euclidean space E^3
        sage: B.display()
        1/3*sqrt(5)*sin(1/3*s) e_x - 1/3*sqrt(5)*cos(1/3*s) e_y + 2/3 e_z

    The three vector fields `(T, N, B)` form the **Frenet-Serret frame** along
    `C`::

        sage: FS = I.vector_frame(('T', 'N', 'B'), (T, N, B),
        ....:                     symbol_dual=('t', 'n', 'b'))
        sage: FS
        Vector frame ((0, 6*pi), (T,N,B)) with values on the Euclidean space E^3

    The Frenet-Serret frame is orthonormal::

        sage: matrix([[u.dot(v).expr() for v in FS] for u in FS])
        [1 0 0]
        [0 1 0]
        [0 0 1]

    The derivative vectors `N'` and `B'`::

        sage: Np = I.vector_field([diff(N[i], s) for i in E.irange()],
        ....:                     dest_map=C, name="N'")
        sage: Np.display()
        N' = 1/3*sin(1/3*s) e_x - 1/3*cos(1/3*s) e_y
        sage: Bp = I.vector_field([diff(B[i], s) for i in E.irange()],
        ....:                     dest_map=C, name="B'")
        sage: Bp.display()
        B' = 1/9*sqrt(5)*cos(1/3*s) e_x + 1/9*sqrt(5)*sin(1/3*s) e_y

    The Frenet-Serret formulas::

        sage: for v in (Tp, Np, Bp):
        ....:     v.display(FS)
        ....:
        T' = 2/9 N
        N' = -2/9 T + 1/9*sqrt(5) B
        B' = -1/9*sqrt(5) N

    The torsion of `C` is obtained as the third component of `N'` in the
    Frenet-Serret frame::

        sage: tau = Np[FS, 3]
        sage: tau
        1/9*sqrt(5)
    """
    def __init__(self, parent, coord_expression=None, name=None,
                 latex_name=None, is_isomorphism=False, is_identity=False):
        r"""
        Construct a curve.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: R.<t> = manifolds.RealLine()
            sage: I = R.open_interval(0, 2*pi)
            sage: c = Hom(I,M)({X: (cos(t), sin(2*t))},  name='c') ; c
            Curve c in the 2-dimensional differentiable manifold M
            sage: TestSuite(c).run()

        The identity of interval ``I``::

            sage: c = Hom(I,I)({}, is_identity=True) ; c
            Identity map Id_(0, 2*pi) of the Real interval (0, 2*pi)
            sage: TestSuite(c).run()
        """
        if coord_expression is None:
            coord_functions = None
        else:
            if not isinstance(coord_expression, dict):
                raise TypeError("{} is not a dictionary".format(
                                                             coord_expression))
            param_chart = parent.domain().canonical_chart()
            coord_functions = {}
            for chart, expr in coord_expression.items():
                if isinstance(chart, tuple):
                    # a pair of charts is passed:
                    coord_functions[chart] = expr
                else:
                    coord_functions[(param_chart, chart)] = expr
        DiffMap.__init__(self, parent, coord_functions=coord_functions,
                         name=name, latex_name=latex_name,
                         is_isomorphism=is_isomorphism,
                         is_identity=is_identity)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: R.<t> = manifolds.RealLine()
            sage: M.curve([cos(t), sin(2*t)], (t, 0, 2*pi))
            Curve in the 2-dimensional differentiable manifold M
            sage: M.curve([cos(t), sin(2*t)], (t, 0, 2*pi), name='c')
            Curve c in the 2-dimensional differentiable manifold M
        """
        if self._codomain._dim == 1:
            return DiffMap._repr_(self)
        description = "Curve "
        if self._name is not None:
            description += self._name + " "
        description += "in the {}".format(self._codomain)
        return description

    def __reduce__(self):
        r"""
        Reduction function for the pickle protocole.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: R.<t> = manifolds.RealLine()
            sage: c = M.curve([cos(t), sin(2*t)], (t, 0, 2*pi))
            sage: c.__reduce__()
            (<class 'sage.manifolds.differentiable.manifold_homset.DifferentiableCurveSet_with_category.element_class'>,
             (Set of Morphisms from Real interval (0, 2*pi) to 2-dimensional
              differentiable manifold M in Category of smooth manifolds over
              Real Field with 53 bits of precision,
              None,
              None,
              None,
              False,
              False))

        Test of pickling::

            sage: loads(dumps(c))
            Curve in the 2-dimensional differentiable manifold M
        """
        return (type(self), (self.parent(), None, self._name, self._latex_name,
                             self._is_isomorphism, self._is_identity))

    def coord_expr(self, chart=None):
        r"""
        Return the coordinate functions expressing the curve in a given chart.

        INPUT:

        - ``chart`` -- (default: ``None``) chart on the curve's codomain; if
          ``None``, the codomain's default chart is assumed

        OUTPUT: symbolic expression representing the curve in the above chart

        EXAMPLES:

        Cartesian and polar expression of a curve in the Euclidean plane::

            sage: M = Manifold(2, 'R^2', r'\RR^2')  # the Euclidean plane R^2
            sage: c_xy.<x,y> = M.chart() # Cartesian coordinate on R^2
            sage: U = M.open_subset('U', coord_def={c_xy: (y!=0, x<0)}) # the complement of the segment y=0 and x>0
            sage: c_cart = c_xy.restrict(U) # Cartesian coordinates on U
            sage: c_spher.<r,ph> = U.chart(r'r:(0,+oo) ph:(0,2*pi):\phi') # spherical coordinates on U

        Links between spherical coordinates and Cartesian ones::

            sage: ch_cart_spher = c_cart.transition_map(c_spher, [sqrt(x*x+y*y), atan2(y,x)])
            sage: ch_cart_spher.set_inverse(r*cos(ph), r*sin(ph))
            Check of the inverse coordinate transformation:
              x == x  *passed*
              y == y  *passed*
              r == r  *passed*
              ph == arctan2(r*sin(ph), r*cos(ph))  **failed**
            NB: a failed report can reflect a mere lack of simplification.
            sage: R.<t> = manifolds.RealLine()
            sage: c = U.curve({c_spher: (1,t)}, (t, 0, 2*pi), name='c')
            sage: c.coord_expr(c_spher)
            (1, t)
            sage: c.coord_expr(c_cart)
            (cos(t), sin(t))

        Since ``c_cart`` is the default chart on ``U``, it can be omitted::

            sage: c.coord_expr()
            (cos(t), sin(t))

        Cartesian expression of a cardioid::

            sage: c = U.curve({c_spher: (2*(1+cos(t)), t)}, (t, 0, 2*pi), name='c')
            sage: c.coord_expr(c_cart)
            (2*cos(t)^2 + 2*cos(t), 2*(cos(t) + 1)*sin(t))
        """
        return self.expr(chart1=self._domain.canonical_chart(), chart2=chart)

    def __call__(self, t, simplify=True):
        r"""
        Image for a given value of the curve parameter.

        This is a redefinition of :meth:`sage.categories.map.Map.__call__`
        to allow for the direct call with some value of the parameter
        (numerical value or symbolic expression) instead of the element
        (ManifoldPoint) of the domain corresponding to that value.

        EXAMPLES:

        Points on circle in the Euclidean plane::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: R.<t> = manifolds.RealLine()
            sage: c = M.curve([cos(t), sin(t)], (t, 0, 2*pi), name='c')
            sage: c(0)
            Point c(0) on the 2-dimensional differentiable manifold M
            sage: c(0) in M
            True
            sage: c(0).coord(X)
            (1, 0)
            sage: c(pi/4).coord(X)
            (1/2*sqrt(2), 1/2*sqrt(2))
            sage: c(t)
            Point c(t) on the 2-dimensional differentiable manifold M
            sage: c(t).coord(X)
            (cos(t), sin(t))
        """
        # Case of a point in the domain:
        if isinstance(t, ManifoldPoint):
            return DiffMap.__call__(self, t)

        # Case of a value of the canonical coordinate in the domain:
        codom = self._codomain
        canon_chart = self._domain._canon_chart
        canon_coord = canon_chart._xx[0]
        if (canon_chart, codom._def_chart) in self._coord_expression:
            chart_pair = (canon_chart, codom._def_chart)
        else:
            chart_pair = next(iter(self._coord_expression.keys()))
            # a chart is picked at random
        coord_functions = self._coord_expression[chart_pair]._functions
        n = codom._dim
        dict_subs = {canon_coord: t}
        coords = [coord_functions[i].expr().substitute(dict_subs)
                  for i in range(n)]
        if simplify:
            coords = [chart_pair[0].simplify(coords[i]) for i in range(n)]
        if self._name is not None:
            name = "{}({})".format(self._name, t)
        else:
            name = None
        if self._latex_name is not None:
            latex_name = r"{}\left({}\right)".format(self._latex_name, latex(t))
        else:
            latex_name = None
        return codom.element_class(codom, coords=coords, chart=chart_pair[1],
                                   name=name, latex_name=latex_name,
                                   check_coords=False)

    def tangent_vector_field(self, name=None, latex_name=None):
        r"""
        Return the tangent vector field to the curve (velocity vector).

        INPUT:

        - ``name`` -- (default: ``None``) string; symbol given to the tangent
          vector field; if none is provided, the primed curve symbol (if any)
          will be used
        - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote
          the tangent vector field; if ``None`` then (i) if ``name`` is
          ``None`` as well, the primed curve LaTeX symbol (if any) will be
          used or (ii) if ``name`` is not ``None``, ``name`` will be used

        OUTPUT:

        - the tangent vector field, as an instance of
          :class:`~sage.manifolds.differentiable.vectorfield.VectorField`

        EXAMPLES:

        Tangent vector field to a circle curve in `\RR^2`::

            sage: M = Manifold(2, 'R^2')
            sage: X.<x,y> = M.chart()
            sage: R.<t> = manifolds.RealLine()
            sage: c = M.curve([cos(t), sin(t)], (t, 0, 2*pi), name='c')
            sage: v = c.tangent_vector_field() ; v
            Vector field c' along the Real interval (0, 2*pi) with values on
             the 2-dimensional differentiable manifold R^2
            sage: v.display()
            c' = -sin(t) ∂/∂x + cos(t) ∂/∂y
            sage: latex(v)
            {c'}
            sage: v.parent()
            Free module X((0, 2*pi),c) of vector fields along the Real interval
             (0, 2*pi) mapped into the 2-dimensional differentiable manifold R^2

        Value of the tangent vector field for some specific value of the
        curve parameter (`t = \pi`)::

            sage: R(pi) in c.domain()  # pi in (0, 2*pi)
            True
            sage: vp = v.at(R(pi)) ; vp
            Tangent vector c' at Point on the 2-dimensional differentiable
             manifold R^2
            sage: vp.parent() is M.tangent_space(c(R(pi)))
            True
            sage: vp.display()
            c' = -∂/∂y

        Tangent vector field to a curve in a non-parallelizable manifold (the
        2-sphere `S^2`): first, we introduce the 2-sphere::

            sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                   intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                   restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: W = U.intersection(V)
            sage: A = W.open_subset('A', coord_def={c_xy.restrict(W): (y!=0, x<0)})
            sage: c_spher.<th,ph> = A.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi') # spherical coordinates
            sage: spher_to_xy = c_spher.transition_map(c_xy.restrict(A),
            ....:           (sin(th)*cos(ph)/(1-cos(th)), sin(th)*sin(ph)/(1-cos(th))) )
            sage: spher_to_xy.set_inverse(2*atan(1/sqrt(x^2+y^2)), atan2(y, x), check=False)

        Then we define a curve (a loxodrome) by its expression in terms of
        spherical coordinates and evaluate the tangent vector field::

            sage: R.<t> = manifolds.RealLine()
            sage: c = M.curve({c_spher: [2*atan(exp(-t/10)), t]}, (t, -oo, +oo),
            ....:              name='c') ; c
            Curve c in the 2-dimensional differentiable manifold M
            sage: vc = c.tangent_vector_field() ; vc
            Vector field c' along the Real number line ℝ with values on
             the 2-dimensional differentiable manifold M
            sage: vc.parent()
            Module X(ℝ,c) of vector fields along the Real number line ℝ
             mapped into the 2-dimensional differentiable manifold M
            sage: vc.display(c_spher.frame().along(c.restrict(R,A)))
            c' = -1/5*e^(1/10*t)/(e^(1/5*t) + 1) ∂/∂th + ∂/∂ph
        """
        vmodule = self._domain.vector_field_module(dest_map=self)
        if latex_name is None:
            if name is None:
                if self._latex_name is not None:
                    latex_name = r"{%s'}" % (self._latex_name)
            else:
                latex_name = name
        if name is None and self._name is not None:
            name = self._name + "'"
        resu = vmodule.element_class(vmodule, name=name, latex_name=latex_name)
        canon_chart = self._domain.canonical_chart()
        codom = self._codomain
        dim = codom._dim
        codom_top_charts = codom._top_charts
        for chart in codom_top_charts:
            try:
                jacob = self.differential_functions(canon_chart, chart)
                restrict = self.restrict(canon_chart.domain(),
                                     subcodomain=chart.domain())
                fmodule = restrict._domain.vector_field_module(dest_map=restrict)
                frame = fmodule.basis(from_frame=chart.frame())
                resu_rest = resu.restrict(canon_chart.domain(), dest_map=restrict)
                resu_rest.add_comp(frame)[:, canon_chart] = [jacob[i][0]
                                                             for i in range(dim)]
            except ValueError:
                pass
        return resu

    @options(thickness=1, plot_points=75, max_range=8, aspect_ratio='automatic')
    def plot(self, chart=None, ambient_coords=None, mapping=None, prange=None,
             include_end_point=(True, True), end_point_offset=(0.001, 0.001),
             parameters=None, color='red', style='-', label_axes=True, **kwds):
        r"""
        Plot the current curve in a Cartesian graph based on the
        coordinates of some ambient chart.

        The curve is drawn in terms of two (2D graphics) or three (3D graphics)
        coordinates of a given chart, called hereafter the *ambient chart*.
        The ambient chart's domain must overlap with the curve's codomain or
        with the codomain of the composite curve `\Phi\circ c`, where `c` is
        the current curve and `\Phi` some manifold differential map (argument
        ``mapping`` below).

        INPUT:

        - ``chart`` -- (default: ``None``) the ambient chart (see above);
          if ``None``, the default chart of the codomain of the curve (or of
          the curve composed with `\Phi`) is used

        - ``ambient_coords`` -- (default: ``None``) tuple containing the 2
          or 3 coordinates of the ambient chart in terms of which the plot
          is performed; if ``None``, all the coordinates of the ambient chart
          are considered

        - ``mapping`` -- (default: ``None``) differentiable mapping `\Phi`
          (instance of
          :class:`~sage.manifolds.differentiable.diff_map.DiffMap`)
          providing the link between the curve and the ambient chart ``chart``
          (cf. above); if ``None``, the ambient chart is supposed to be defined
          on the codomain of the curve.

        - ``prange`` -- (default: ``None``) range of the curve parameter for
          the plot; if ``None``, the entire parameter range declared during the
          curve construction is considered (with -Infinity
          replaced by ``-max_range`` and +Infinity by ``max_range``)

        - ``include_end_point`` -- (default: ``(True, True)``) determines
          whether the end points of ``prange`` are included in the plot

        - ``end_point_offset`` -- (default: ``(0.001, 0.001)``) offsets from
          the end points when they are not included in the plot: if
          ``include_end_point[0] == False``, the minimal value of the curve
          parameter used for the plot is ``prange[0] + end_point_offset[0]``,
          while if ``include_end_point[1] == False``, the maximal value is
          ``prange[1] - end_point_offset[1]``.

        - ``max_range`` -- (default: 8) numerical value substituted to
          +Infinity if the latter is the upper bound of the parameter range;
          similarly ``-max_range`` is the numerical valued substituted for
          -Infinity

        - ``parameters`` -- (default: ``None``) dictionary giving the numerical
          values of the parameters that may appear in the coordinate expression
          of the curve

        - ``color`` -- (default: ``'red'``) color of the drawn curve

        - ``style`` -- (default: ``'-'``) color of the drawn curve; NB: ``style``
          is effective only for 2D plots

        - ``thickness`` -- (default: 1) thickness of the drawn curve

        - ``plot_points`` -- (default: 75) number of points to plot the curve

        - ``label_axes`` -- boolean (default: ``True``); determining whether the
          labels of the coordinate axes of ``chart`` shall be added to the
          graph; can be set to ``False`` if the graph is 3D and must be
          superposed with another graph.

        - ``aspect_ratio`` -- (default: ``'automatic'``) aspect ratio of the
          plot; the default value (``'automatic'``) applies only for 2D plots;
          for 3D plots, the default value is ``1`` instead

        OUTPUT:

        - a graphic object, either an instance of
          :class:`~sage.plot.graphics.Graphics` for a 2D plot (i.e. based on
          2 coordinates of ``chart``) or an instance of
          :class:`~sage.plot.plot3d.base.Graphics3d` for a 3D plot (i.e.
          based on 3 coordinates of ``chart``)

        EXAMPLES:

        Plot of the lemniscate of Gerono::

            sage: R2 = Manifold(2, 'R^2')
            sage: X.<x,y> = R2.chart()
            sage: R.<t> = manifolds.RealLine()
            sage: c = R2.curve([sin(t), sin(2*t)/2], (t, 0, 2*pi), name='c')
            sage: c.plot()  # 2D plot
            Graphics object consisting of 1 graphics primitive

        .. PLOT::

            R2 = Manifold(2, 'R^2')
            X = R2.chart('x y')
            t = manifolds.RealLine().canonical_coordinate()
            c = R2.curve([sin(t), sin(2*t)/2], (t, 0, 2*pi), name='c')
            g = c.plot()
            sphinx_plot(g)

        Plot for a subinterval of the curve's domain::

            sage: c.plot(prange=(0,pi))
            Graphics object consisting of 1 graphics primitive

        .. PLOT::

            R2 = Manifold(2, 'R^2')
            X = R2.chart('x y')
            t = manifolds.RealLine().canonical_coordinate()
            c = R2.curve([sin(t), sin(2*t)/2], (t, 0, 2*pi), name='c')
            g = c.plot(prange=(0,pi))
            sphinx_plot(g)

        Plot with various options::

            sage: c.plot(color='green', style=':', thickness=3, aspect_ratio=1)
            Graphics object consisting of 1 graphics primitive

        .. PLOT::

            R2 = Manifold(2, 'R^2')
            X = R2.chart('x y')
            t = manifolds.RealLine().canonical_coordinate()
            c = R2.curve([sin(t), sin(2*t)/2], (t, 0, 2*pi), name='c')
            g = c.plot(color='green', style=':', thickness=3, aspect_ratio=1)
            sphinx_plot(g)

        Cardioid defined in terms of polar coordinates and plotted with respect
        to Cartesian coordinates, as an example of use of the optional argument
        ``chart``::

            sage: E.<r,ph> = EuclideanSpace(coordinates='polar')
            sage: c = E.curve((1 + cos(ph), ph), (ph, 0, 2*pi))
            sage: c.plot(chart=E.cartesian_coordinates(), aspect_ratio=1)
            Graphics object consisting of 1 graphics primitive

        .. PLOT::

            E = EuclideanSpace(2, coordinates='polar')
            r, ph = E.polar_coordinates()[:]
            c = E.curve((1 + cos(ph), ph), (ph, 0, 2*pi))
            g = c.plot(chart=E.cartesian_coordinates(), aspect_ratio=1)
            sphinx_plot(g)

        Plot via a mapping to another manifold: loxodrome of a sphere viewed
        in `\RR^3`::

            sage: S2 = Manifold(2, 'S^2')
            sage: U = S2.open_subset('U')
            sage: XS.<th,ph> = U.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi')
            sage: R3 = Manifold(3, 'R^3')
            sage: X3.<x,y,z> = R3.chart()
            sage: F = S2.diff_map(R3, {(XS, X3): [sin(th)*cos(ph),
            ....:                      sin(th)*sin(ph), cos(th)]}, name='F')
            sage: F.display()
            F: S^2 → R^3
            on U: (th, ph) ↦ (x, y, z) = (cos(ph)*sin(th), sin(ph)*sin(th), cos(th))
            sage: c = S2.curve([2*atan(exp(-t/10)), t], (t, -oo, +oo), name='c')
            sage: graph_c = c.plot(mapping=F, max_range=40,
            ....:                  plot_points=200, thickness=2, label_axes=False)  # 3D plot
            sage: graph_S2 = XS.plot(X3, mapping=F, number_values=11, color='black') # plot of the sphere
            sage: show(graph_c + graph_S2) # the loxodrome + the sphere

        .. PLOT::

            S2 = Manifold(2, 'S^2')
            U = S2.open_subset('U')
            XS = U.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi')
            th, ph = XS[:]
            R3 = Manifold(3, 'R^3')
            X3 = R3.chart('x y z')
            F = S2.diff_map(R3, {(XS, X3): [sin(th)*cos(ph), sin(th)*sin(ph),
                                            cos(th)]}, name='F')
            t = manifolds.RealLine().canonical_coordinate()
            c = S2.curve([2*atan(exp(-t/10)), t], (t, -oo, +oo), name='c')
            graph_c = c.plot(mapping=F, max_range=40, plot_points=200,
                             thickness=2, label_axes=False)
            graph_S2 = XS.plot(X3, mapping=F, number_values=11, color='black')
            sphinx_plot(graph_c + graph_S2)

        Example of use of the argument ``parameters``: we define a curve with
        some symbolic parameters ``a`` and ``b``::

            sage: a, b = var('a b')
            sage: c = R2.curve([a*cos(t) + b, a*sin(t)], (t, 0, 2*pi), name='c')

        To make a plot, we set specific values for ``a`` and ``b`` by means
        of the Python dictionary ``parameters``::

            sage: c.plot(parameters={a: 2, b: -3}, aspect_ratio=1)
            Graphics object consisting of 1 graphics primitive

        .. PLOT::

            R2 = Manifold(2, 'R^2')
            X = R2.chart('x y')
            t = manifolds.RealLine().canonical_coordinate()
            a, b = var('a b')
            c = R2.curve([a*cos(t) + b, a*sin(t)], (t, 0, 2*pi), name='c')
            g = c.plot(parameters={a: 2, b: -3}, aspect_ratio=1)
            sphinx_plot(g)
        """
        from sage.rings.infinity import Infinity
        from sage.misc.functional import numerical_approx
        from sage.manifolds.chart import RealChart

        #
        # Get the @options from kwds
        #
        thickness = kwds.pop('thickness')
        plot_points = kwds.pop('plot_points')
        max_range = kwds.pop('max_range')
        aspect_ratio = kwds.pop('aspect_ratio')
        #
        # The "effective" curve to be plotted
        #
        if mapping is None:
            eff_curve = self
        else:
            eff_curve = mapping.restrict(self.codomain()) * self
        #
        # The chart w.r.t. which the curve is plotted
        #
        if chart is None:
            chart = eff_curve._codomain.default_chart()
        elif not isinstance(chart, RealChart):
            raise TypeError("{} is not a real chart".format(chart))
        #
        # Coordinates of the above chart w.r.t. which the curve is plotted
        #
        if ambient_coords is None:
            ambient_coords = chart[:]  # all chart coordinates are used
        n_pc = len(ambient_coords)
        if n_pc != 2 and n_pc != 3:
            raise ValueError("the number of coordinates involved in the " +
                             "plot must be either 2 or 3, not {}".format(n_pc))
        # indices of plot coordinates
        ind_pc = [chart[:].index(pc) for pc in ambient_coords]
        #
        # Parameter range for the plot
        #
        if prange is None:
            prange = (self._domain.lower_bound(), self._domain.upper_bound())
        elif not isinstance(prange, (tuple, list)):
            raise TypeError("{} is neither a tuple nor a list".format(prange))
        elif len(prange) != 2:
            raise ValueError("the argument prange must be a tuple/list " +
                             "of 2 elements")
        tmin = prange[0]
        tmax = prange[1]
        if tmin == -Infinity:
            tmin = -max_range
        elif not include_end_point[0]:
            tmin = tmin + end_point_offset[0]
        if tmax == Infinity:
            tmax = max_range
        elif not include_end_point[1]:
            tmax = tmax - end_point_offset[1]
        tmin = numerical_approx(tmin)
        tmax = numerical_approx(tmax)
        #
        # The coordinate expression of the effective curve
        #
        transf = eff_curve.coord_functions(chart1=self._domain.canonical_chart(),
                                           chart2=chart)
        #
        # List of points for the plot curve
        #
        plot_curve = []
        dt = (tmax - tmin) / (plot_points - 1)
        t = tmin
        if parameters is None:
            for i in range(plot_points):
                x = transf(t, simplify=False)
                plot_curve.append( [numerical_approx(x[j]) for j in ind_pc] )
                t += dt
        else:
            for i in range(plot_points):
                x = transf(t, simplify=False)
                plot_curve.append(
                               [numerical_approx( x[j].substitute(parameters) )
                                for j in ind_pc] )
                t += dt

        return self._graphics(plot_curve, ambient_coords,
                              thickness=thickness,
                              aspect_ratio=aspect_ratio, color=color,
                              style=style, label_axes=label_axes)

    def _graphics(self, plot_curve, ambient_coords, thickness=1,
                  aspect_ratio='automatic', color='red', style='-',
                  label_axes=True):
        r"""
        Plot a 2D or 3D curve in a Cartesian graph with axes labeled by
        the ambient coordinates; it is invoked by the methods
        :meth:`plot` of
        :class:`~sage.manifolds.differentiable.curve.DifferentiableCurve`,
        and its subclasses
        (:class:`~sage.manifolds.differentiable.integrated_curve.IntegratedCurve`,
        :class:`~sage.manifolds.differentiable.integrated_curve.IntegratedAutoparallelCurve`,
        and
        :class:`~sage.manifolds.differentiable.integrated_curve.IntegratedGeodesic`).

        TESTS::

            sage: M = Manifold(2, 'R^2')
            sage: X.<x,y> = M.chart()
            sage: R.<t> = manifolds.RealLine()
            sage: c = M.curve([cos(t), sin(t)], (t, 0, 2*pi), name='c')
            sage: graph = c._graphics([[1,2], [3,4]], [x,y])
            sage: graph._objects[0].xdata == [1,3]
            True
            sage: graph._objects[0].ydata == [2,4]
            True
            sage: graph._objects[0]._options['thickness'] == 1
            True
            sage: graph._extra_kwds['aspect_ratio'] == 'automatic'
            True
            sage: graph._objects[0]._options['rgbcolor'] == 'red'
            True
            sage: graph._objects[0]._options['linestyle'] == '-'
            True
            sage: l = [r'$'+latex(x)+r'$', r'$'+latex(y)+r'$']
            sage: graph._extra_kwds['axes_labels'] == l
            True
        """
        from sage.plot.graphics import Graphics
        from sage.plot.line import line
        from sage.manifolds.utilities import set_axes_labels

        #
        # The plot
        #
        n_pc = len(ambient_coords)
        resu = Graphics()
        resu += line(plot_curve, color=color, linestyle=style,
                     thickness=thickness)
        if n_pc == 2:  # 2D graphic
            resu.set_aspect_ratio(aspect_ratio)
            if label_axes:
                # We update the dictionary _extra_kwds (options to be passed
                # to show()), instead of using the method
                # Graphics.axes_labels() since the latter is not robust w.r.t.
                # graph addition
                resu._extra_kwds['axes_labels'] = [r'$'+latex(pc)+r'$'
                                                   for pc in ambient_coords]
        else: # 3D graphic
            if aspect_ratio == 'automatic':
                aspect_ratio = 1
            resu.aspect_ratio(aspect_ratio)
            if label_axes:
                labels = [str(pc) for pc in ambient_coords]
                resu = set_axes_labels(resu, *labels)
        return resu
