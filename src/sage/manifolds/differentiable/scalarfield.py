# sage_setup: distribution = sagemath-symbolics
r"""
Differentiable Scalar Fields

Given a differentiable manifold `M` of class `C^k` over a topological field `K`
(in most applications, `K = \RR` or `K = \CC`), a *differentiable scalar field*
on `M` is a map

.. MATH::

    f: M \longrightarrow K

of class `C^k`.

Differentiable scalar fields are implemented by the class
:class:`DiffScalarField`.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013-2015): initial version
- Eric Gourgoulhon (2018): operators gradient, Laplacian and d'Alembertian

REFERENCES:

- [KN1963]_
- [Lee2013]_
- [ONe1983]_
"""

#******************************************************************************
#       Copyright (C) 2015, 2018 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from __future__ import annotations

from typing import TYPE_CHECKING, Union

from sage.manifolds.scalarfield import ScalarField

if TYPE_CHECKING:
    from sage.manifolds.differentiable.diff_form import DiffForm
    from sage.manifolds.differentiable.metric import PseudoRiemannianMetric
    from sage.manifolds.differentiable.symplectic_form import SymplecticForm


class DiffScalarField(ScalarField):
    r"""
    Differentiable scalar field on a differentiable manifold.

    Given a differentiable manifold `M` of class `C^k` over a topological field
    `K` (in most applications, `K = \RR` or `K = \CC`), a *differentiable
    scalar field* defined on `M` is a map

    .. MATH::

        f: M \longrightarrow K

    that is `k`-times continuously differentiable.

    The class :class:`DiffScalarField` is a Sage *element* class, whose
    *parent* class is
    :class:`~sage.manifolds.differentiable.scalarfield_algebra.DiffScalarFieldAlgebra`.
    It inherits from the class :class:`~sage.manifolds.scalarfield.ScalarField`
    devoted to generic continuous scalar fields on topological manifolds.

    INPUT:

    - ``parent`` -- the algebra of scalar fields containing the scalar field
      (must be an instance of class
      :class:`~sage.manifolds.differentiable.scalarfield_algebra.DiffScalarFieldAlgebra`)
    - ``coord_expression`` -- (default: ``None``) coordinate expression(s) of
      the scalar field; this can be either

      - a dictionary of coordinate expressions in various charts on the domain,
        with the charts as keys;
      - a single coordinate expression; if the argument ``chart`` is
        ``'all'``, this expression is set to all the charts defined
        on the open set; otherwise, the expression is set in the
        specific chart provided by the argument ``chart``

      NB: If ``coord_expression`` is ``None`` or incomplete, coordinate
      expressions can be added after the creation of the object, by means of
      the methods
      :meth:`~sage.manifolds.scalarfield.ScalarField.add_expr`,
      :meth:`~sage.manifolds.scalarfield.ScalarField.add_expr_by_continuation`
      and :meth:`~sage.manifolds.scalarfield.ScalarField.set_expr`
    - ``chart`` -- (default: ``None``) chart defining the coordinates used
      in ``coord_expression`` when the latter is a single coordinate
      expression; if none is provided (default), the default chart of the
      open set is assumed. If ``chart=='all'``, ``coord_expression`` is
      assumed to be independent of the chart (constant scalar field).
    - ``name`` -- (default: ``None``) string; name (symbol) given to the
      scalar field
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote the
      scalar field; if none is provided, the LaTeX symbol is set to ``name``

    EXAMPLES:

    A scalar field on the 2-sphere::

        sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
        sage: U = M.open_subset('U') # complement of the North pole
        sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
        sage: V = M.open_subset('V') # complement of the South pole
        sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
        sage: M.declare_union(U,V)   # S^2 is the union of U and V
        sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
        ....:                                intersection_name='W',
        ....:                                restrictions1= x^2+y^2!=0,
        ....:                                restrictions2= u^2+v^2!=0)
        sage: uv_to_xy = xy_to_uv.inverse()
        sage: f = M.scalar_field({c_xy: 1/(1+x^2+y^2), c_uv: (u^2+v^2)/(1+u^2+v^2)},
        ....:                    name='f') ; f
        Scalar field f on the 2-dimensional differentiable manifold M
        sage: f.display()
        f: M → ℝ
        on U: (x, y) ↦ 1/(x^2 + y^2 + 1)
        on V: (u, v) ↦ (u^2 + v^2)/(u^2 + v^2 + 1)

    For scalar fields defined by a single coordinate expression, the latter
    can be passed instead of the dictionary over the charts::

        sage: g = U.scalar_field(x*y, chart=c_xy, name='g') ; g
        Scalar field g on the Open subset U of the 2-dimensional differentiable
         manifold M

    The above is indeed equivalent to::

        sage: g = U.scalar_field({c_xy: x*y}, name='g') ; g
        Scalar field g on the Open subset U of the 2-dimensional differentiable
         manifold M

    Since ``c_xy`` is the default chart of ``U``, the argument ``chart`` can
    be skipped::

        sage: g = U.scalar_field(x*y, name='g') ; g
        Scalar field g on the Open subset U of the 2-dimensional differentiable
         manifold M

    The scalar field `g` is defined on `U` and has an expression in terms of
    the coordinates `(u,v)` on `W=U\cap V`::

        sage: g.display()
        g: U → ℝ
           (x, y) ↦ x*y
        on W: (u, v) ↦ u*v/(u^4 + 2*u^2*v^2 + v^4)

    Scalar fields on `M` can also be declared with a single chart::

        sage: f = M.scalar_field(1/(1+x^2+y^2), chart=c_xy, name='f') ; f
        Scalar field f on the 2-dimensional differentiable manifold M

    Their definition must then be completed by providing the expressions on
    other charts, via the method
    :meth:`~sage.manifolds.scalarfield.ScalarField.add_expr`, to get a global
    cover of the manifold::

        sage: f.add_expr((u^2+v^2)/(1+u^2+v^2), chart=c_uv)
        sage: f.display()
        f: M → ℝ
        on U: (x, y) ↦ 1/(x^2 + y^2 + 1)
        on V: (u, v) ↦ (u^2 + v^2)/(u^2 + v^2 + 1)

    We can even first declare the scalar field without any coordinate
    expression and provide them subsequently::

        sage: f = M.scalar_field(name='f')
        sage: f.add_expr(1/(1+x^2+y^2), chart=c_xy)
        sage: f.add_expr((u^2+v^2)/(1+u^2+v^2), chart=c_uv)
        sage: f.display()
        f: M → ℝ
        on U: (x, y) ↦ 1/(x^2 + y^2 + 1)
        on V: (u, v) ↦ (u^2 + v^2)/(u^2 + v^2 + 1)

    We may also use the method
    :meth:`~sage.manifolds.scalarfield.ScalarField.add_expr_by_continuation`
    to complete the coordinate definition using the analytic continuation from
    domains in which charts overlap::

        sage: f = M.scalar_field(1/(1+x^2+y^2), chart=c_xy, name='f') ; f
        Scalar field f on the 2-dimensional differentiable manifold M
        sage: f.add_expr_by_continuation(c_uv, U.intersection(V))
        sage: f.display()
        f: M → ℝ
        on U: (x, y) ↦ 1/(x^2 + y^2 + 1)
        on V: (u, v) ↦ (u^2 + v^2)/(u^2 + v^2 + 1)

    A scalar field can also be defined by some unspecified function of the
    coordinates::

        sage: h = U.scalar_field(function('H')(x, y), name='h') ; h
        Scalar field h on the Open subset U of the 2-dimensional differentiable
         manifold M
        sage: h.display()
        h: U → ℝ
           (x, y) ↦ H(x, y)
        on W: (u, v) ↦ H(u/(u^2 + v^2), v/(u^2 + v^2))

    We may use the argument ``latex_name`` to specify the LaTeX symbol denoting
    the scalar field if the latter is different from ``name``::

        sage: latex(f)
        f
        sage: f = M.scalar_field({c_xy: 1/(1+x^2+y^2), c_uv: (u^2+v^2)/(1+u^2+v^2)},
        ....:                    name='f', latex_name=r'\mathcal{F}')
        sage: latex(f)
        \mathcal{F}

    The coordinate expression in a given chart is obtained via the method
    :meth:`~sage.manifolds.scalarfield.ScalarField.expr`, which returns a
    symbolic expression::

        sage: f.expr(c_uv)
        (u^2 + v^2)/(u^2 + v^2 + 1)
        sage: type(f.expr(c_uv))
        <class 'sage.symbolic.expression.Expression'>

    The method :meth:`~sage.manifolds.scalarfield.ScalarField.coord_function`
    returns instead a function of the chart coordinates, i.e. an instance of
    :class:`~sage.manifolds.chart_func.ChartFunction`::

        sage: f.coord_function(c_uv)
        (u^2 + v^2)/(u^2 + v^2 + 1)
        sage: type(f.coord_function(c_uv))
        <class 'sage.manifolds.chart_func.ChartFunctionRing_with_category.element_class'>
        sage: f.coord_function(c_uv).display()
        (u, v) ↦ (u^2 + v^2)/(u^2 + v^2 + 1)

    The value returned by the method
    :meth:`~sage.manifolds.scalarfield.ScalarField.expr`
    is actually the coordinate expression of the chart function::

        sage: f.expr(c_uv) is f.coord_function(c_uv).expr()
        True

    A constant scalar field is declared by setting the argument ``chart`` to
    ``'all'``::

        sage: c = M.scalar_field(2, chart='all', name='c') ; c
        Scalar field c on the 2-dimensional differentiable manifold M
        sage: c.display()
        c: M → ℝ
        on U: (x, y) ↦ 2
        on V: (u, v) ↦ 2

    A shortcut is to use the method
    :meth:`~sage.manifolds.manifold.TopologicalManifold.constant_scalar_field`::

        sage: c == M.constant_scalar_field(2)
        True

    The constant value can be some unspecified parameter::

        sage: var('a')
        a
        sage: c = M.constant_scalar_field(a, name='c') ; c
        Scalar field c on the 2-dimensional differentiable manifold M
        sage: c.display()
        c: M → ℝ
        on U: (x, y) ↦ a
        on V: (u, v) ↦ a

    A special case of constant field is the zero scalar field::

        sage: zer = M.constant_scalar_field(0) ; zer
        Scalar field zero on the 2-dimensional differentiable manifold M
        sage: zer.display()
        zero: M → ℝ
        on U: (x, y) ↦ 0
        on V: (u, v) ↦ 0

    It can be obtained directly by means of the function
    :meth:`~sage.manifolds.manifold.TopologicalManifold.zero_scalar_field`::

        sage: zer is M.zero_scalar_field()
        True

    A third way is to get it as the zero element of the algebra `C^k(M)`
    of scalar fields on `M` (see below)::

        sage: zer is M.scalar_field_algebra().zero()
        True

    By definition, a scalar field acts on the manifold's points, sending
    them to elements of the manifold's base field (real numbers in the
    present case)::

        sage: N = M.point((0,0), chart=c_uv) # the North pole
        sage: S = M.point((0,0), chart=c_xy) # the South pole
        sage: E = M.point((1,0), chart=c_xy) # a point at the equator
        sage: f(N)
        0
        sage: f(S)
        1
        sage: f(E)
        1/2
        sage: h(E)
        H(1, 0)
        sage: c(E)
        a
        sage: zer(E)
        0

    A scalar field can be compared to another scalar field::

        sage: f == g
        False

    ...to a symbolic expression::

        sage: f == x*y
        False
        sage: g == x*y
        True
        sage: c == a
        True

    ...to a number::

        sage: f == 2
        False
        sage: zer == 0
        True

    ...to anything else::

        sage: f == M
        False

    Standard mathematical functions are implemented::

        sage: sqrt(f)
        Scalar field sqrt(f) on the 2-dimensional differentiable manifold M
        sage: sqrt(f).display()
        sqrt(f): M → ℝ
        on U: (x, y) ↦ 1/sqrt(x^2 + y^2 + 1)
        on V: (u, v) ↦ sqrt(u^2 + v^2)/sqrt(u^2 + v^2 + 1)

    ::

        sage: tan(f)
        Scalar field tan(f) on the 2-dimensional differentiable manifold M
        sage: tan(f).display()
        tan(f): M → ℝ
        on U: (x, y) ↦ sin(1/(x^2 + y^2 + 1))/cos(1/(x^2 + y^2 + 1))
        on V: (u, v) ↦ sin((u^2 + v^2)/(u^2 + v^2 + 1))/cos((u^2 + v^2)/(u^2 + v^2 + 1))

    .. RUBRIC:: Arithmetics of scalar fields

    Scalar fields on `M` (resp. `U`) belong to the algebra `C^k(M)`
    (resp. `C^k(U)`)::

        sage: f.parent()
        Algebra of differentiable scalar fields on the 2-dimensional
         differentiable manifold M
        sage: f.parent() is M.scalar_field_algebra()
        True
        sage: g.parent()
        Algebra of differentiable scalar fields on the Open subset U of the
         2-dimensional differentiable manifold M
        sage: g.parent() is U.scalar_field_algebra()
        True

    Consequently, scalar fields can be added::

        sage: s = f + c ; s
        Scalar field f+c on the 2-dimensional differentiable manifold M
        sage: s.display()
        f+c: M → ℝ
        on U: (x, y) ↦ (a*x^2 + a*y^2 + a + 1)/(x^2 + y^2 + 1)
        on V: (u, v) ↦ ((a + 1)*u^2 + (a + 1)*v^2 + a)/(u^2 + v^2 + 1)

    and subtracted::

        sage: s = f - c ; s
        Scalar field f-c on the 2-dimensional differentiable manifold M
        sage: s.display()
        f-c: M → ℝ
        on U: (x, y) ↦ -(a*x^2 + a*y^2 + a - 1)/(x^2 + y^2 + 1)
        on V: (u, v) ↦ -((a - 1)*u^2 + (a - 1)*v^2 + a)/(u^2 + v^2 + 1)

    Some tests::

        sage: f + zer == f
        True
        sage: f - f == zer
        True
        sage: f + (-f) == zer
        True
        sage: (f+c)-f == c
        True
        sage: (f-c)+c == f
        True

    We may add a number (interpreted as a constant scalar field) to a scalar
    field::

        sage: s = f + 1 ; s
        Scalar field f+1 on the 2-dimensional differentiable manifold M
        sage: s.display()
        f+1: M → ℝ
        on U: (x, y) ↦ (x^2 + y^2 + 2)/(x^2 + y^2 + 1)
        on V: (u, v) ↦ (2*u^2 + 2*v^2 + 1)/(u^2 + v^2 + 1)
        sage: (f+1)-1 == f
        True

    The number can represented by a symbolic variable::

        sage: s = a + f ; s
        Scalar field on the 2-dimensional differentiable manifold M
        sage: s == c + f
        True

    However if the symbolic variable is a chart coordinate, the addition
    is performed only on the chart domain::

        sage: s = f + x; s
        Scalar field on the 2-dimensional differentiable manifold M
        sage: s.display()
        M → ℝ
        on U: (x, y) ↦ (x^3 + x*y^2 + x + 1)/(x^2 + y^2 + 1)
        on W: (u, v) ↦ (u^4 + v^4 + u^3 + (2*u^2 + u)*v^2 + u)/(u^4 + v^4 + (2*u^2 + 1)*v^2 + u^2)
        sage: s = f + u; s
        Scalar field on the 2-dimensional differentiable manifold M
        sage: s.display()
        M → ℝ
        on W: (x, y) ↦ (x^3 + (x + 1)*y^2 + x^2 + x)/(x^4 + y^4 + (2*x^2 + 1)*y^2 + x^2)
        on V: (u, v) ↦ (u^3 + (u + 1)*v^2 + u^2 + u)/(u^2 + v^2 + 1)

    The addition of two scalar fields with different domains is possible if
    the domain of one of them is a subset of the domain of the other; the
    domain of the result is then this subset::

        sage: f.domain()
        2-dimensional differentiable manifold M
        sage: g.domain()
        Open subset U of the 2-dimensional differentiable manifold M
        sage: s = f + g ; s
        Scalar field f+g on the Open subset U of the 2-dimensional
         differentiable manifold M
        sage: s.domain()
        Open subset U of the 2-dimensional differentiable manifold M
        sage: s.display()
        f+g: U → ℝ
           (x, y) ↦ (x*y^3 + (x^3 + x)*y + 1)/(x^2 + y^2 + 1)
        on W: (u, v) ↦ (u^6 + 3*u^4*v^2 + 3*u^2*v^4 + v^6 + u*v^3
         + (u^3 + u)*v)/(u^6 + v^6 + (3*u^2 + 1)*v^4 + u^4 + (3*u^4 + 2*u^2)*v^2)

    The operation actually performed is `f|_U + g`::

        sage: s == f.restrict(U) + g
        True

    In Sage framework, the addition of `f` and `g` is permitted because
    there is a *coercion* of the parent of `f`, namely `C^k(M)`, to
    the parent of `g`, namely `C^k(U)` (see
    :class:`~sage.manifolds.differentiable.scalarfield_algebra.DiffScalarFieldAlgebra`)::

        sage: CM = M.scalar_field_algebra()
        sage: CU = U.scalar_field_algebra()
        sage: CU.has_coerce_map_from(CM)
        True

    The coercion map is nothing but the restriction to domain `U`::

        sage: CU.coerce(f) == f.restrict(U)
        True

    Since the algebra `C^k(M)` is a vector space over `\RR`, scalar fields
    can be multiplied by a number, either an explicit one::

        sage: s = 2*f ; s
        Scalar field on the 2-dimensional differentiable manifold M
        sage: s.display()
        M → ℝ
        on U: (x, y) ↦ 2/(x^2 + y^2 + 1)
        on V: (u, v) ↦ 2*(u^2 + v^2)/(u^2 + v^2 + 1)

    or a symbolic one::

        sage: s = a*f ; s
        Scalar field on the 2-dimensional differentiable manifold M
        sage: s.display()
        M → ℝ
        on U: (x, y) ↦ a/(x^2 + y^2 + 1)
        on V: (u, v) ↦ (u^2 + v^2)*a/(u^2 + v^2 + 1)

    However, if the symbolic variable is a chart coordinate, the multiplication
    is performed only in the corresponding chart::

        sage: s = x*f; s
        Scalar field on the 2-dimensional differentiable manifold M
        sage: s.display()
        M → ℝ
        on U: (x, y) ↦ x/(x^2 + y^2 + 1)
        on W: (u, v) ↦ u/(u^2 + v^2 + 1)
        sage: s = u*f; s
        Scalar field on the 2-dimensional differentiable manifold M
        sage: s.display()
        M → ℝ
        on W: (x, y) ↦ x/(x^4 + y^4 + (2*x^2 + 1)*y^2 + x^2)
        on V: (u, v) ↦ (u^2 + v^2)*u/(u^2 + v^2 + 1)

    Some tests::

        sage: 0*f == 0
        True
        sage: 0*f == zer
        True
        sage: 1*f == f
        True
        sage: (-2)*f == - f - f
        True

    The ring multiplication of the algebras `C^k(M)` and `C^k(U)`
    is the pointwise multiplication of functions::

        sage: s = f*f ; s
        Scalar field f*f on the 2-dimensional differentiable manifold M
        sage: s.display()
        f*f: M → ℝ
        on U: (x, y) ↦ 1/(x^4 + y^4 + 2*(x^2 + 1)*y^2 + 2*x^2 + 1)
        on V: (u, v) ↦ (u^4 + 2*u^2*v^2 + v^4)/(u^4 + v^4 + 2*(u^2 + 1)*v^2 + 2*u^2 + 1)
        sage: s = g*h ; s
        Scalar field g*h on the Open subset U of the 2-dimensional
         differentiable manifold M
        sage: s.display()
        g*h: U → ℝ
           (x, y) ↦ x*y*H(x, y)
        on W: (u, v) ↦ u*v*H(u/(u^2 + v^2), v/(u^2 + v^2))/(u^4 + 2*u^2*v^2 + v^4)

    Thanks to the coercion `C^k(M)\rightarrow C^k(U)` mentioned
    above, it is possible to multiply a scalar field defined on `M` by a
    scalar field defined on `U`, the result being a scalar field defined on
    `U`::

        sage: f.domain(), g.domain()
        (2-dimensional differentiable manifold M,
         Open subset U of the 2-dimensional differentiable manifold M)
        sage: s = f*g ; s
        Scalar field f*g on the Open subset U of the 2-dimensional
         differentiable manifold M
        sage: s.display()
        f*g: U → ℝ
           (x, y) ↦ x*y/(x^2 + y^2 + 1)
        on W: (u, v) ↦ u*v/(u^4 + v^4 + (2*u^2 + 1)*v^2 + u^2)
        sage: s == f.restrict(U)*g
        True

    Scalar fields can be divided (pointwise division)::

        sage: s = f/c ; s
        Scalar field f/c on the 2-dimensional differentiable manifold M
        sage: s.display()
        f/c: M → ℝ
        on U: (x, y) ↦ 1/(a*x^2 + a*y^2 + a)
        on V: (u, v) ↦ (u^2 + v^2)/(a*u^2 + a*v^2 + a)
        sage: s = g/h ; s
        Scalar field g/h on the Open subset U of the 2-dimensional
         differentiable manifold M
        sage: s.display()
        g/h: U → ℝ
           (x, y) ↦ x*y/H(x, y)
        on W: (u, v) ↦ u*v/((u^4 + 2*u^2*v^2 + v^4)*H(u/(u^2 + v^2), v/(u^2 + v^2)))
        sage: s = f/g ; s
        Scalar field f/g on the Open subset U of the 2-dimensional
         differentiable manifold M
        sage: s.display()
        f/g: U → ℝ
           (x, y) ↦ 1/(x*y^3 + (x^3 + x)*y)
        on W: (u, v) ↦ (u^6 + 3*u^4*v^2 + 3*u^2*v^4 + v^6)/(u*v^3 + (u^3 + u)*v)
        sage: s == f.restrict(U)/g
        True

    For scalar fields defined on a single chart domain, we may perform some
    arithmetics with symbolic expressions involving the chart coordinates::

        sage: s = g + x^2 - y ; s
        Scalar field on the Open subset U of the 2-dimensional differentiable
         manifold M
        sage: s.display()
        U → ℝ
        (x, y) ↦ x^2 + (x - 1)*y
        on W: (u, v) ↦ -(v^3 - u^2 + (u^2 - u)*v)/(u^4 + 2*u^2*v^2 + v^4)

    ::

        sage: s = g*x ; s
        Scalar field on the Open subset U of the 2-dimensional differentiable
         manifold M
        sage: s.display()
        U → ℝ
        (x, y) ↦ x^2*y
        on W: (u, v) ↦ u^2*v/(u^6 + 3*u^4*v^2 + 3*u^2*v^4 + v^6)

    ::

        sage: s = g/x ; s
        Scalar field on the Open subset U of the 2-dimensional differentiable
         manifold M
        sage: s.display()
        U → ℝ
        (x, y) ↦ y
        on W: (u, v) ↦ v/(u^2 + v^2)
        sage: s = x/g ; s
        Scalar field on the Open subset U of the 2-dimensional differentiable
         manifold M
        sage: s.display()
        U → ℝ
        (x, y) ↦ 1/y
        on W: (u, v) ↦ (u^2 + v^2)/v

    The test suite is passed::

        sage: TestSuite(f).run()
        sage: TestSuite(zer).run()
    """
    def __init__(self, parent, coord_expression=None, chart=None, name=None,
                 latex_name=None):
        r"""
        Construct a scalar field.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x+y}, name='f') ; f
            Scalar field f on the 2-dimensional differentiable manifold M
            sage: from sage.manifolds.scalarfield import ScalarField
            sage: isinstance(f, ScalarField)
            True
            sage: f.parent()
            Algebra of differentiable scalar fields on the 2-dimensional
             differentiable manifold M
            sage: TestSuite(f).run()
        """
        ScalarField.__init__(self, parent, coord_expression=coord_expression,
                             chart=chart, name=name, latex_name=latex_name)
        self._tensor_type = (0,0)
        self._tensor_rank = 0

    ####### Required methods for an algebra element (beside arithmetic) #######

    def _init_derived(self):
        r"""
        Initialize the derived quantities.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x+y})
            sage: f._init_derived()
        """
        ScalarField._init_derived(self) # derived quantities of the parent class
        self._differential = None  # differential 1-form of the scalar field
        self._lie_derivatives = {} # dict. of Lie derivatives of self, (keys: id(vector))

    def _del_derived(self):
        r"""
        Delete the derived quantities.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x+y})
            sage: U = M.open_subset('U', coord_def={X: x>0})
            sage: f.restrict(U)
            Scalar field on the Open subset U of the 2-dimensional
             differentiable manifold M
            sage: f._restrictions
            {Open subset U of the 2-dimensional differentiable manifold M:
             Scalar field on the Open subset U of the 2-dimensional
             differentiable manifold M}
            sage: f._del_derived()
            sage: f._restrictions  # restrictions are derived quantities
            {}
        """
        ScalarField._del_derived(self) # derived quantities of the mother class
        self._differential = None  # reset of the differential
        # First deletes any reference to self in the vectors' dictionaries:
        for val in self._lie_derivatives.values():
            del val[0]._lie_der_along_self[id(self)]
        # Then clears the dictionary of Lie derivatives
        self._lie_derivatives.clear()

    def tensor_type(self):
        r"""
        Return the tensor type of ``self``, when the latter is considered
        as a tensor field on the manifold. This is always `(0, 0)`.

        OUTPUT:

        - always `(0, 0)`

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: f = M.scalar_field(x+2*y)
            sage: f.tensor_type()
            (0, 0)
        """
        return self._tensor_type

    def differential(self) -> DiffForm:
        r"""
        Return the differential of ``self``.

        OUTPUT:

        - a :class:`~sage.manifolds.differentiable.diff_form.DiffForm` (or of
          :class:`~sage.manifolds.differentiable.diff_form.DiffFormParal` if
          the scalar field's domain is parallelizable) representing the 1-form
          that is the differential of the scalar field

        EXAMPLES:

        Differential of a scalar field on a 3-dimensional differentiable
        manifold::

            sage: M = Manifold(3, 'M')
            sage: c_xyz.<x,y,z> = M.chart()
            sage: f = M.scalar_field(cos(x)*z^3 + exp(y)*z^2, name='f')
            sage: df = f.differential() ; df
            1-form df on the 3-dimensional differentiable manifold M
            sage: df.display()
            df = -z^3*sin(x) dx + z^2*e^y dy + (3*z^2*cos(x) + 2*z*e^y) dz
            sage: latex(df)
            \mathrm{d}f
            sage: df.parent()
            Free module Omega^1(M) of 1-forms on the 3-dimensional
             differentiable manifold M

        The result is cached, i.e. is not recomputed unless ``f`` is changed::

            sage: f.differential() is df
            True

        Instead of invoking the method :meth:`differential`, one may apply the
        function ``diff`` to the scalar field::

            sage: diff(f) is f.differential()
            True

        Since the exterior derivative of a scalar field (considered a 0-form)
        is nothing but its differential, ``exterior_derivative()`` is an
        alias of ``differential()``::

            sage: df = f.exterior_derivative() ; df
            1-form df on the 3-dimensional differentiable manifold M
            sage: df.display()
            df = -z^3*sin(x) dx + z^2*e^y dy + (3*z^2*cos(x) + 2*z*e^y) dz
            sage: latex(df)
            \mathrm{d}f

        Differential computed on a chart that is not the default one::

            sage: c_uvw.<u,v,w> = M.chart()
            sage: g = M.scalar_field(u*v^2*w^3, c_uvw, name='g')
            sage: dg = g.differential() ; dg
            1-form dg on the 3-dimensional differentiable manifold M
            sage: dg._components
            {Coordinate frame (M, (∂/∂u,∂/∂v,∂/∂w)): 1-index components w.r.t.
             Coordinate frame (M, (∂/∂u,∂/∂v,∂/∂w))}
            sage: dg.comp(c_uvw.frame())[:, c_uvw]
            [v^2*w^3, 2*u*v*w^3, 3*u*v^2*w^2]
            sage: dg.display(c_uvw)
            dg = v^2*w^3 du + 2*u*v*w^3 dv + 3*u*v^2*w^2 dw

        The exterior derivative is nilpotent::

            sage: ddf = df.exterior_derivative() ; ddf
            2-form ddf on the 3-dimensional differentiable manifold M
            sage: ddf == 0
            True
            sage: ddf[:] # for the incredule
            [0 0 0]
            [0 0 0]
            [0 0 0]
            sage: ddg = dg.exterior_derivative() ; ddg
            2-form ddg on the 3-dimensional differentiable manifold M
            sage: ddg == 0
            True
        """
        from sage.tensor.modules.format_utilities import (format_unop_txt,
                                                          format_unop_latex)
        if self._differential is None:
            # A new computation is necessary:
            rname = format_unop_txt('d', self._name)
            rlname = format_unop_latex(r'\mathrm{d}', self._latex_name)
            self._differential = self._domain.one_form(name=rname,
                                                       latex_name=rlname)
            if self._is_zero:
                for chart in self._domain._atlas:
                    self._differential.add_comp(chart._frame) # since a newly
                                            # created set of components is zero
            else:
                for chart, func in self._express.items():
                    diff_func = self._differential.add_comp(chart._frame)
                    for i in self._manifold.irange():
                        diff_func[i, chart] = func.diff(i)
        return self._differential

    exterior_derivative = differential  # a scalar field being a 0-form
    derivative = differential  # allows one to use functional notation,
                               # e.g. diff(f) for f.differential()

    def lie_derivative(self, vector):
        r"""
        Compute the Lie derivative with respect to a vector field.

        In the present case (scalar field), the Lie derivative is equal to
        the scalar field resulting from the action of the vector field on
        the scalar field.

        INPUT:

        - ``vector`` -- vector field with respect to which the Lie derivative
          is to be taken

        OUTPUT:

        - the scalar field that is the Lie derivative of the scalar field with
          respect to ``vector``

        EXAMPLES:

        Lie derivative on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: f = M.scalar_field(x^2*cos(y))
            sage: v = M.vector_field(name='v')
            sage: v[:] = (-y, x)
            sage: f.lie_derivative(v)
            Scalar field on the 2-dimensional differentiable manifold M
            sage: f.lie_derivative(v).expr()
            -x^3*sin(y) - 2*x*y*cos(y)

        The result is cached::

            sage: f.lie_derivative(v) is f.lie_derivative(v)
            True

        An alias is ``lie_der``::

            sage: f.lie_der(v) is f.lie_derivative(v)
            True

        Alternative expressions of the Lie derivative of a scalar field::

            sage: f.lie_der(v) == v(f)  # the vector acting on f
            True
            sage: f.lie_der(v) == f.differential()(v)  # the differential of f acting on the vector
            True

        A vanishing Lie derivative::

            sage: f.set_expr(x^2 + y^2)
            sage: f.lie_der(v).display()
            M → ℝ
            (x, y) ↦ 0
        """
        # The Lie derivative is cached in _lie_derivatives if neither
        #   the scalar field nor ``vector`` have been modified.
        if id(vector) not in self._lie_derivatives:
            # A new computation must be performed
            res = vector(self)
            self._lie_derivatives[id(vector)] = (vector, res)
            vector._lie_der_along_self[id(self)] = self
        return self._lie_derivatives[id(vector)][1]

    lie_der = lie_derivative

    def hodge_dual(
        self, nondegenerate_tensor: Union[PseudoRiemannianMetric, SymplecticForm]
    ) -> DiffForm:
        r"""
        Compute the Hodge dual of the scalar field with respect to some non-degenerate
        bilinear form (Riemannian metric or symplectic form).

        If `M` is the domain of the scalar field (denoted by `f`), `n` is the
        dimension of `M` and `g` is a non-degenerate bilinear form on `M`, the
        *Hodge dual* of `f` w.r.t. `g` is the `n`-form `*f` defined by

        .. MATH::

            *f = f \epsilon,

        where `\epsilon` is the volume `n`-form associated with `g` (see
        :meth:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric.volume_form`).

        INPUT:

        - ``nondegenerate_tensor`` -- a non-degenerate bilinear form defined on the same manifold
          as the current differential form; must be an instance of
          :class:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric` or
          :class:`~sage.manifolds.differentiable.symplectic_form.SymplecticForm`.

        OUTPUT: the `n`-form `*f`

        EXAMPLES:

        Hodge dual of a scalar field in the Euclidean space `R^3`::

            sage: M = Manifold(3, 'M', start_index=1)
            sage: X.<x,y,z> = M.chart()
            sage: g = M.metric('g')
            sage: g[1,1], g[2,2], g[3,3] = 1, 1, 1
            sage: f = M.scalar_field(function('F')(x,y,z), name='f')
            sage: sf = f.hodge_dual(g) ; sf
            3-form *f on the 3-dimensional differentiable manifold M
            sage: sf.display()
            *f = F(x, y, z) dx∧dy∧dz
            sage: ssf = sf.hodge_dual(g) ; ssf
            Scalar field **f on the 3-dimensional differentiable manifold M
            sage: ssf.display()
            **f: M → ℝ
               (x, y, z) ↦ F(x, y, z)
            sage: ssf == f # must hold for a Riemannian metric
            True

        Instead of calling the method :meth:`hodge_dual` on the scalar field,
        one can invoke the method
        :meth:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric.hodge_star`
        of the metric::

            sage: f.hodge_dual(g) == g.hodge_star(f)
            True
        """
        from sage.tensor.modules.format_utilities import (
            format_unop_txt,
            format_unop_latex,
        )

        result = self * nondegenerate_tensor.volume_form()
        result.set_name(
            name=format_unop_txt("*", self._name),
            latex_name=format_unop_latex(r"\star ", self._latex_name),
        )
        return result

    def bracket(self, other):
        r"""
        Return the Schouten-Nijenhuis bracket of ``self``, considered as a
        multivector field of degree 0, with a multivector field.

        See
        :meth:`~sage.manifolds.differentiable.multivectorfield.MultivectorFieldParal.bracket`
        for details.

        INPUT:

        - ``other`` -- a multivector field of degree `p`

        OUTPUT:

        - if `p=0`, a zero scalar field
        - if `p=1`, an instance of :class:`DiffScalarField` representing
          the Schouten-Nijenhuis bracket ``[self,other]``
        - if `p\geq 2`, an instance of
          :class:`~sage.manifolds.differentiable.multivectorfield.MultivectorField`
          representing the Schouten-Nijenhuis bracket ``[self,other]``

        EXAMPLES:

        The Schouten-Nijenhuis bracket of two scalar fields is identically
        zero::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x+y^2}, name='f')
            sage: g = M.scalar_field({X: y-x}, name='g')
            sage: s = f.bracket(g); s
            Scalar field zero on the 2-dimensional differentiable manifold M
            sage: s.display()
            zero: M → ℝ
               (x, y) ↦ 0

        while the Schouten-Nijenhuis bracket of a scalar field `f` with a
        multivector field `a` is equal to minus the interior product of the
        differential of `f` with `a`::

            sage: a = M.multivector_field(2, name='a')
            sage: a[0,1] = x*y ; a.display()
            a = x*y ∂/∂x∧∂/∂y
            sage: s = f.bracket(a); s
            Vector field -i_df a on the 2-dimensional differentiable manifold M
            sage: s.display()
            -i_df a = 2*x*y^2 ∂/∂x - x*y ∂/∂y

        See
        :meth:`~sage.manifolds.differentiable.multivectorfield.MultivectorFieldParal.bracket`
        for other examples.
        """
        if isinstance(other, DiffScalarField):
            return self._domain.intersection(other._domain).zero_scalar_field()
        return - self.differential().interior_product(other)

    def wedge(self, other):
        r"""
        Return the exterior product of ``self``, considered as a differential
        form of degree 0 or a multivector field of degree 0, with ``other``.

        See
        :meth:`~sage.manifolds.differentiable.diff_form.DiffFormParal.wedge`
        (exterior product of differential forms) or
        :meth:`~sage.manifolds.differentiable.multivectorfield.MultivectorFieldParal.wedge`
        (exterior product of multivector fields) for details.

        For a scalar field `f` and a `p`-form (or `p`-vector field) `a`, the
        exterior product reduces to the standard product on the left by an
        element of the base ring of the module of `p`-forms (or `p`-vector
        fields): `f\wedge a = f a`.

        INPUT:

        - ``other`` -- a differential form or a multivector field `a`

        OUTPUT: the product `f a`, where `f` is ``self``

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x+y^2}, name='f')
            sage: a = M.diff_form(2, name='a')
            sage: a[0,1] = x*y
            sage: s = f.wedge(a); s
            2-form f*a on the 2-dimensional differentiable manifold M
            sage: s.display()
            f*a = (x*y^3 + x^2*y) dx∧dy
        """
        return self * other

    def degree(self):
        r"""
        Return the degree of ``self``, considered as a differential
        form or a multivector field, i.e. zero.

        This trivial method is provided for consistency with the exterior
        calculus scheme, cf. the methods
        :meth:`~sage.manifolds.differentiable.diff_form.DiffForm.degree`
        (differential forms) and
        :meth:`~sage.manifolds.differentiable.multivectorfield.MultivectorField.degree`
        (multivector fields).

        OUTPUT: 0

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x+y^2})
            sage: f.degree()
            0
        """
        return self._tensor_rank

    def gradient(self, metric=None):
        r"""
        Return the gradient of ``self`` (with respect to a given metric).

        The *gradient* of a scalar field `f` with respect to a metric `g`
        is the vector field `\mathrm{grad}\, f` whose components in any
        coordinate frame are

        .. MATH::

            (\mathrm{grad}\, f)^i = g^{ij} \frac{\partial F}{\partial x^j}

        where the `x^j`'s are the coordinates with respect to which the
        frame is defined and `F` is the chart function representing `f` in
        these coordinates: `f(p) = F(x^1(p),\ldots,x^n(p))` for any point `p`
        in the chart domain.
        In other words, the gradient of `f` is the vector field that is the
        `g`-dual of the differential of `f`.

        INPUT:

        - ``metric`` -- (default: ``None``) the pseudo-Riemannian metric `g`
          involved in the definition of the gradient; if none is provided, the
          domain of ``self`` is supposed to be endowed with a default metric
          (i.e. is supposed to be pseudo-Riemannian manifold, see
          :class:`~sage.manifolds.differentiable.pseudo_riemannian.PseudoRiemannianManifold`)
          and the latter is used to define the gradient

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.differentiable.vectorfield.VectorField`
          representing the gradient of ``self``

        EXAMPLES:

        Gradient of a scalar field in the Euclidean plane::

            sage: M.<x,y> = EuclideanSpace()
            sage: f = M.scalar_field(cos(x*y), name='f')
            sage: v = f.gradient(); v
            Vector field grad(f) on the Euclidean plane E^2
            sage: v.display()
            grad(f) = -y*sin(x*y) e_x - x*sin(x*y) e_y
            sage: v[:]
            [-y*sin(x*y), -x*sin(x*y)]

        Gradient in polar coordinates::

            sage: M.<r,phi> = EuclideanSpace(coordinates='polar')
            sage: f = M.scalar_field(r*cos(phi), name='f')
            sage: f.gradient().display()
            grad(f) = cos(phi) e_r - sin(phi) e_phi
            sage: f.gradient()[:]
            [cos(phi), -sin(phi)]

        Note that ``(e_r, e_phi)`` is the orthonormal vector frame associated
        with polar coordinates (see
        :meth:`~sage.manifolds.differentiable.examples.euclidean.EuclideanPlane.polar_frame`);
        the gradient expressed in the coordinate frame is::

            sage: f.gradient().display(M.polar_coordinates().frame())
            grad(f) = cos(phi) ∂/∂r - sin(phi)/r ∂/∂phi

        The function :func:`~sage.manifolds.operators.grad` from the
        :mod:`~sage.manifolds.operators` module can be used instead of the
        method :meth:`gradient`::

            sage: from sage.manifolds.operators import grad
            sage: grad(f) == f.gradient()
            True

        The gradient can be taken with respect to a metric tensor that is
        not the default one::

            sage: h = M.lorentzian_metric('h')
            sage: h[1,1], h[2,2] = -1, 1/(1+r^2)
            sage: h.display(M.polar_coordinates().frame())
            h = -dr⊗dr + r^2/(r^2 + 1) dphi⊗dphi
            sage: v = f.gradient(h); v
            Vector field grad_h(f) on the Euclidean plane E^2
            sage: v.display()
            grad_h(f) = -cos(phi) e_r + (-r^2*sin(phi) - sin(phi)) e_phi
        """
        default_metric = metric is None
        if default_metric:
            metric = self._domain.metric()
        resu = self.differential().up(metric)
        if self._name is not None:
            if default_metric:
                resu._name = "grad({})".format(self._name)
                resu._latex_name = r"\mathrm{grad}\left(" + \
                                   self._latex_name + r"\right)"
            else:
                resu._name = "grad_{}({})".format(metric._name, self._name)
                resu._latex_name = r"\mathrm{grad}_{" + metric._latex_name + \
                                   r"}\left(" + self._latex_name + r"\right)"
            # The name is propagated to possible restrictions of self:
            for restrict in resu._restrictions.values():
                restrict.set_name(resu._name, latex_name=resu._latex_name)
        return resu

    def laplacian(self, metric=None):
        r"""
        Return the Laplacian of ``self`` with respect to a given
        metric (Laplace-Beltrami operator).

        The *Laplacian* of a scalar field `f` with respect to a metric `g`
        is the scalar field

        .. MATH::

            \Delta f  = g^{ij} \nabla_i \nabla_j f = \nabla_i \nabla^i f

        where `\nabla` is the Levi-Civita connection of `g`.
        `\Delta` is also called the *Laplace-Beltrami operator*.

        INPUT:

        - ``metric`` -- (default: ``None``) the pseudo-Riemannian metric `g`
          involved in the definition of the Laplacian; if none is provided, the
          domain of ``self`` is supposed to be endowed with a default metric
          (i.e. is supposed to be pseudo-Riemannian manifold, see
          :class:`~sage.manifolds.differentiable.pseudo_riemannian.PseudoRiemannianManifold`)
          and the latter is used to define the Laplacian

        OUTPUT:

        - instance of :class:`DiffScalarField` representing the Laplacian of
          ``self``

        EXAMPLES:

        Laplacian of a scalar field on the Euclidean plane::

            sage: M.<x,y> = EuclideanSpace()
            sage: f = M.scalar_field(function('F')(x,y), name='f')
            sage: s = f.laplacian(); s
            Scalar field Delta(f) on the Euclidean plane E^2
            sage: s.display()
            Delta(f): E^2 → ℝ
               (x, y) ↦ d^2(F)/dx^2 + d^2(F)/dy^2

        The function :func:`~sage.manifolds.operators.laplacian` from the
        :mod:`~sage.manifolds.operators` module can be used instead of the
        method :meth:`laplacian`::

            sage: from sage.manifolds.operators import laplacian
            sage: laplacian(f) == s
            True

        The Laplacian can be taken with respect to a metric tensor that is
        not the default one::

            sage: h = M.lorentzian_metric('h')
            sage: h[1,1], h[2,2] = -1, 1/(1+x^2+y^2)
            sage: s = f.laplacian(h); s
            Scalar field Delta_h(f) on the Euclidean plane E^2
            sage: s.display()
            Delta_h(f): E^2 → ℝ
               (x, y) ↦ (y^4*d^2(F)/dy^2 + y^3*d(F)/dy
               + (2*(x^2 + 1)*d^2(F)/dy^2 - d^2(F)/dx^2)*y^2
               + (x^2 + 1)*y*d(F)/dy + x*d(F)/dx - (x^2 + 1)*d^2(F)/dx^2
               + (x^4 + 2*x^2 + 1)*d^2(F)/dy^2)/(x^2 + y^2 + 1)

        The Laplacian of `f` is equal to the divergence of the gradient of `f`:

        .. MATH::

            \Delta f = \mathrm{div}( \mathrm{grad}\,  f )

        Let us check this formula::

            sage: s == f.gradient(h).div(h)
            True
        """
        default_metric = metric is None
        if default_metric:
            metric = self._domain.metric()
        nabla = metric.connection()
        resu = nabla(self.differential().up(metric)).trace()
        if self._name is not None:
            if default_metric:
                resu._name = "Delta({})".format(self._name)
                resu._latex_name = r"\Delta\left(" + self._latex_name + \
                                   r"\right)"
            else:
                resu._name = "Delta_{}({})".format(metric._name, self._name)
                resu._latex_name = r"\Delta_{" + metric._latex_name + \
                                   r"}\left(" + self._latex_name + r"\right)"
            # The name is propagated to possible restrictions of self:
            for restrict in resu._restrictions.values():
                restrict.set_name(resu._name, latex_name=resu._latex_name)
        return resu

    def dalembertian(self, metric=None):
        r"""
        Return the d'Alembertian of ``self`` with respect to a given
        Lorentzian metric.

        The *d'Alembertian* of a scalar field `f` with respect to a Lorentzian
        metric `g` is nothing but the Laplacian (see :meth:`laplacian`) of `f`
        with respect to that metric:

        .. MATH::

            \Box f  = g^{ij} \nabla_i \nabla_j f = \nabla_i \nabla^i f

        where `\nabla` is the Levi-Civita connection of `g`.

        .. NOTE::

            If the metric `g` is not Lorentzian, the name *d'Alembertian* is
            not appropriate and one should use :meth:`laplacian` instead.

        INPUT:

        - ``metric`` -- (default: ``None``) the Lorentzian metric `g`
          involved in the definition of the d'Alembertian; if none is provided,
          the domain of ``self`` is supposed to be endowed with a default
          Lorentzian metric (i.e. is supposed to be Lorentzian manifold, see
          :class:`~sage.manifolds.differentiable.pseudo_riemannian.PseudoRiemannianManifold`)
          and the latter is used to define the d'Alembertian

        OUTPUT:

        - instance of :class:`DiffScalarField` representing the d'Alembertian
          of ``self``

        EXAMPLES:

        d'Alembertian of a scalar field in Minkowski spacetime::

            sage: M = Manifold(4, 'M', structure='Lorentzian')
            sage: X.<t,x,y,z> = M.chart()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1, 1, 1, 1
            sage: f = M.scalar_field(t + x^2 + t^2*y^3 - x*z^4, name='f')
            sage: s = f.dalembertian(); s
            Scalar field Box(f) on the 4-dimensional Lorentzian manifold M
            sage: s.display()
            Box(f): M → ℝ
               (t, x, y, z) ↦ 6*t^2*y - 2*y^3 - 12*x*z^2 + 2

        The function :func:`~sage.manifolds.operators.dalembertian` from the
        :mod:`~sage.manifolds.operators` module can be used instead of the
        method :meth:`dalembertian`::

            sage: from sage.manifolds.operators import dalembertian
            sage: dalembertian(f) == s
            True
        """
        default_metric = metric is None
        if default_metric:
            metric = self._domain.metric()
        nm2 = self._manifold.dim() - 2
        if metric.signature() not in [nm2, -nm2]:
            raise TypeError("the {} is not a Lorentzian ".format(metric) +
                            "metric; use laplacian() instead")
        nabla = metric.connection()
        resu = nabla(self.differential().up(metric)).trace()
        if self._name is not None:
            if default_metric:
                resu._name = "Box({})".format(self._name)
                resu._latex_name = r"\Box\left(" + self._latex_name + \
                                   r"\right)"
            else:
                resu._name = "Box_{}({})".format(metric._name, self._name)
                resu._latex_name = r"\Box_{" + metric._latex_name + \
                                   r"}\left(" + self._latex_name + r"\right)"
            # The name is propagated to possible restrictions of self:
            for restrict in resu._restrictions.values():
                restrict.set_name(resu._name, latex_name=resu._latex_name)
        return resu
