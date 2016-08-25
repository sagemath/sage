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

REFERENCES:

.. [1] \S. Kobayashi & K. Nomizu : *Foundations of Differential Geometry*,
   vol. 1, Interscience Publishers (New York) (1963)
.. [2] \J.M. Lee : *Introduction to Smooth Manifolds*, 2nd ed., Springer
   (New York) (2013)
.. [3] \B. O'Neill : *Semi-Riemannian Geometry*, Academic Press (San Diego)
   (1983)

"""

#******************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.manifolds.scalarfield import ScalarField

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
        f: M --> R
        on U: (x, y) |--> 1/(x^2 + y^2 + 1)
        on V: (u, v) |--> (u^2 + v^2)/(u^2 + v^2 + 1)

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
        g: U --> R
           (x, y) |--> x*y
        on W: (u, v) |--> u*v/(u^4 + 2*u^2*v^2 + v^4)

    Scalar fields on `M` can also be declared with a single chart::

        sage: f = M.scalar_field(1/(1+x^2+y^2), chart=c_xy, name='f') ; f
        Scalar field f on the 2-dimensional differentiable manifold M

    Their definition must then be completed by providing the expressions on
    other charts, via the method
    :meth:`~sage.manifolds.scalarfield.ScalarField.add_expr`, to get a global
    cover of the manifold::

        sage: f.add_expr((u^2+v^2)/(1+u^2+v^2), chart=c_uv)
        sage: f.display()
        f: M --> R
        on U: (x, y) |--> 1/(x^2 + y^2 + 1)
        on V: (u, v) |--> (u^2 + v^2)/(u^2 + v^2 + 1)

    We can even first declare the scalar field without any coordinate
    expression and provide them subsequently::

        sage: f = M.scalar_field(name='f')
        sage: f.add_expr(1/(1+x^2+y^2), chart=c_xy)
        sage: f.add_expr((u^2+v^2)/(1+u^2+v^2), chart=c_uv)
        sage: f.display()
        f: M --> R
        on U: (x, y) |--> 1/(x^2 + y^2 + 1)
        on V: (u, v) |--> (u^2 + v^2)/(u^2 + v^2 + 1)

    We may also use the method
    :meth:`~sage.manifolds.scalarfield.ScalarField.add_expr_by_continuation`
    to complete the coordinate definition using the analytic continuation from
    domains in which charts overlap::

        sage: f = M.scalar_field(1/(1+x^2+y^2), chart=c_xy, name='f') ; f
        Scalar field f on the 2-dimensional differentiable manifold M
        sage: f.add_expr_by_continuation(c_uv, U.intersection(V))
        sage: f.display()
        f: M --> R
        on U: (x, y) |--> 1/(x^2 + y^2 + 1)
        on V: (u, v) |--> (u^2 + v^2)/(u^2 + v^2 + 1)

    A scalar field can also be defined by some unspecified function of the
    coordinates::

        sage: h = U.scalar_field(function('H')(x, y), name='h') ; h
        Scalar field h on the Open subset U of the 2-dimensional differentiable
         manifold M
        sage: h.display()
        h: U --> R
           (x, y) |--> H(x, y)
        on W: (u, v) |--> H(u/(u^2 + v^2), v/(u^2 + v^2))

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
        <type 'sage.symbolic.expression.Expression'>

    The method :meth:`~sage.manifolds.scalarfield.ScalarField.coord_function`
    returns instead a function of the chart coordinates, i.e. an instance of
    :class:`~sage.manifolds.coord_func.CoordFunction`::

        sage: f.coord_function(c_uv)
        (u^2 + v^2)/(u^2 + v^2 + 1)
        sage: type(f.coord_function(c_uv))
        <class 'sage.manifolds.coord_func_symb.CoordFunctionSymbRing_with_category.element_class'>
        sage: f.coord_function(c_uv).display()
        (u, v) |--> (u^2 + v^2)/(u^2 + v^2 + 1)

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
        c: M --> R
        on U: (x, y) |--> 2
        on V: (u, v) |--> 2

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
        c: M --> R
        on U: (x, y) |--> a
        on V: (u, v) |--> a

    A special case of constant field is the zero scalar field::

        sage: zer = M.constant_scalar_field(0) ; zer
        Scalar field zero on the 2-dimensional differentiable manifold M
        sage: zer.display()
        zero: M --> R
        on U: (x, y) |--> 0
        on V: (u, v) |--> 0

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
        sqrt(f): M --> R
        on U: (x, y) |--> 1/sqrt(x^2 + y^2 + 1)
        on V: (u, v) |--> sqrt(u^2 + v^2)/sqrt(u^2 + v^2 + 1)

    ::

        sage: tan(f)
        Scalar field tan(f) on the 2-dimensional differentiable manifold M
        sage: tan(f).display()
        tan(f): M --> R
        on U: (x, y) |--> sin(1/(x^2 + y^2 + 1))/cos(1/(x^2 + y^2 + 1))
        on V: (u, v) |--> sin((u^2 + v^2)/(u^2 + v^2 + 1))/cos((u^2 + v^2)/(u^2 + v^2 + 1))

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
        f+c: M --> R
        on U: (x, y) |--> (a*x^2 + a*y^2 + a + 1)/(x^2 + y^2 + 1)
        on V: (u, v) |--> ((a + 1)*u^2 + (a + 1)*v^2 + a)/(u^2 + v^2 + 1)

    and subtracted::

        sage: s = f - c ; s
        Scalar field f-c on the 2-dimensional differentiable manifold M
        sage: s.display()
        f-c: M --> R
        on U: (x, y) |--> -(a*x^2 + a*y^2 + a - 1)/(x^2 + y^2 + 1)
        on V: (u, v) |--> -((a - 1)*u^2 + (a - 1)*v^2 + a)/(u^2 + v^2 + 1)

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
        Scalar field on the 2-dimensional differentiable manifold M
        sage: s.display()
        M --> R
        on U: (x, y) |--> (x^2 + y^2 + 2)/(x^2 + y^2 + 1)
        on V: (u, v) |--> (2*u^2 + 2*v^2 + 1)/(u^2 + v^2 + 1)
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
        M --> R
        on U: (x, y) |--> (x^3 + x*y^2 + x + 1)/(x^2 + y^2 + 1)
        sage: s = f + u; s
        Scalar field on the 2-dimensional differentiable manifold M
        sage: s.display()
        M --> R
        on V: (u, v) |--> (u^3 + (u + 1)*v^2 + u^2 + u)/(u^2 + v^2 + 1)

    The addition of two scalar fields with different domains is possible if
    the domain of one of them is a subset of the domain of the other; the
    domain of the result is then this subset::

        sage: f.domain()
        2-dimensional differentiable manifold M
        sage: g.domain()
        Open subset U of the 2-dimensional differentiable manifold M
        sage: s = f + g ; s
        Scalar field on the Open subset U of the 2-dimensional differentiable
         manifold M
        sage: s.domain()
        Open subset U of the 2-dimensional differentiable manifold M
        sage: s.display()
        U --> R
        (x, y) |--> (x*y^3 + (x^3 + x)*y + 1)/(x^2 + y^2 + 1)
        on W: (u, v) |--> (u^6 + 3*u^4*v^2 + 3*u^2*v^4 + v^6 + u*v^3
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
        M --> R
        on U: (x, y) |--> 2/(x^2 + y^2 + 1)
        on V: (u, v) |--> 2*(u^2 + v^2)/(u^2 + v^2 + 1)

    or a symbolic one::

        sage: s = a*f ; s
        Scalar field on the 2-dimensional differentiable manifold M
        sage: s.display()
        M --> R
        on U: (x, y) |--> a/(x^2 + y^2 + 1)
        on V: (u, v) |--> (u^2 + v^2)*a/(u^2 + v^2 + 1)

    However, if the symbolic variable is a chart coordinate, the multiplication
    is performed only in the corresponding chart::

        sage: s = x*f; s
        Scalar field on the 2-dimensional differentiable manifold M
        sage: s.display()
        M --> R
        on U: (x, y) |--> x/(x^2 + y^2 + 1)
        sage: s = u*f; s
        Scalar field on the 2-dimensional differentiable manifold M
        sage: s.display()
        M --> R
        on V: (u, v) |--> (u^2 + v^2)*u/(u^2 + v^2 + 1)

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
        f*f: M --> R
        on U: (x, y) |--> 1/(x^4 + y^4 + 2*(x^2 + 1)*y^2 + 2*x^2 + 1)
        on V: (u, v) |--> (u^4 + 2*u^2*v^2 + v^4)/(u^4 + v^4 + 2*(u^2 + 1)*v^2 + 2*u^2 + 1)
        sage: s = g*h ; s
        Scalar field g*h on the Open subset U of the 2-dimensional
         differentiable manifold M
        sage: s.display()
        g*h: U --> R
           (x, y) |--> x*y*H(x, y)
        on W: (u, v) |--> u*v*H(u/(u^2 + v^2), v/(u^2 + v^2))/(u^4 + 2*u^2*v^2 + v^4)

    Thanks to the coercion `C^k(M)\rightarrow C^k(U)` mentionned
    above, it is possible to multiply a scalar field defined on `M` by a
    scalar field defined on `U`, the result being a scalar field defined on
    `U`::

        sage: f.domain(), g.domain()
        (2-dimensional differentiable manifold M,
         Open subset U of the 2-dimensional differentiable manifold M)
        sage: s = f*g ; s
        Scalar field on the Open subset U of the 2-dimensional differentiable
         manifold M
        sage: s.display()
        U --> R
        (x, y) |--> x*y/(x^2 + y^2 + 1)
        on W: (u, v) |--> u*v/(u^4 + v^4 + (2*u^2 + 1)*v^2 + u^2)
        sage: s == f.restrict(U)*g
        True

    Scalar fields can be divided (pointwise division)::

        sage: s = f/c ; s
        Scalar field f/c on the 2-dimensional differentiable manifold M
        sage: s.display()
        f/c: M --> R
        on U: (x, y) |--> 1/(a*x^2 + a*y^2 + a)
        on V: (u, v) |--> (u^2 + v^2)/(a*u^2 + a*v^2 + a)
        sage: s = g/h ; s
        Scalar field g/h on the Open subset U of the 2-dimensional
         differentiable manifold M
        sage: s.display()
        g/h: U --> R
           (x, y) |--> x*y/H(x, y)
        on W: (u, v) |--> u*v/((u^4 + 2*u^2*v^2 + v^4)*H(u/(u^2 + v^2), v/(u^2 + v^2)))
        sage: s = f/g ; s
        Scalar field on the Open subset U of the 2-dimensional differentiable
         manifold M
        sage: s.display()
        U --> R
        (x, y) |--> 1/(x*y^3 + (x^3 + x)*y)
        on W: (u, v) |--> (u^6 + 3*u^4*v^2 + 3*u^2*v^4 + v^6)/(u*v^3 + (u^3 + u)*v)
        sage: s == f.restrict(U)/g
        True

    For scalar fields defined on a single chart domain, we may perform some
    arithmetics with symbolic expressions involving the chart coordinates::

        sage: s = g + x^2 - y ; s
        Scalar field on the Open subset U of the 2-dimensional differentiable
         manifold M
        sage: s.display()
        U --> R
        (x, y) |--> x^2 + (x - 1)*y
        on W: (u, v) |--> -(v^3 - u^2 + (u^2 - u)*v)/(u^4 + 2*u^2*v^2 + v^4)

    ::

        sage: s = g*x ; s
        Scalar field on the Open subset U of the 2-dimensional differentiable
         manifold M
        sage: s.display()
        U --> R
        (x, y) |--> x^2*y
        on W: (u, v) |--> u^2*v/(u^6 + 3*u^4*v^2 + 3*u^2*v^4 + v^6)

    ::

        sage: s = g/x ; s
        Scalar field on the Open subset U of the 2-dimensional differentiable
         manifold M
        sage: s.display()
        U --> R
        (x, y) |--> y
        on W: (u, v) |--> v/(u^2 + v^2)
        sage: s = x/g ; s
        Scalar field on the Open subset U of the 2-dimensional differentiable
         manifold M
        sage: s.display()
        U --> R
        (x, y) |--> 1/y
        on W: (u, v) |--> (u^2 + v^2)/v

    The test suite is passed::

        sage: TestSuite(f).run()
        sage: TestSuite(zer).run()

    """
    def __init__(self, parent, coord_expression=None, chart=None, name=None,
                 latex_name=None):
        r"""
        Construct a scalar field.

        TEST::

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

    ####### Required methods for an algebra element (beside arithmetic) #######

    def _init_derived(self):
        r"""
        Initialize the derived quantities.

        TEST::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x+y})
            sage: f._init_derived()

        """
        ScalarField._init_derived(self) # derived quantities of the mother class
        # More to come in ticket #18843

    def _del_derived(self):
        r"""
        Delete the derived quantities.

        TEST::

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
        # More to come in ticket #18843

