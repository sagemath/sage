r"""
Scalar Fields

Given a topological manifold `M` over a topological field `K` (in most
applications, `K = \RR` or `K = \CC`), a *scalar field* on `M` is a
continuous map

.. MATH::

    f:  M \longrightarrow K

Scalar fields are implemented by the class :class:`ScalarField`.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013-2015): initial version
- Travis Scrimshaw (2016): review tweaks

REFERENCES:

- [Lee11]_ \J.M. Lee : *Introduction to Topological Manifolds*, 2nd ed.,
  Springer (New York) (2011)
- [KN63]_ \S. Kobayashi & K. Nomizu : *Foundations of Differential Geometry*,
  vol. 1, Interscience Publishers (New York) (1963)

"""

#******************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#       Copyright (C) 2016 Travis Scrimshaw <tscrimsh@umn.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from six import itervalues

from sage.structure.element import CommutativeAlgebraElement
from sage.symbolic.expression import Expression
from sage.manifolds.coord_func import CoordFunction

class ScalarField(CommutativeAlgebraElement):
    r"""
    Scalar field on a topological manifold.

    Given a topological manifold `M` over a topological field `K` (in most
    applications, `K = \RR` or `K = \CC`), a *scalar field on* `M` is a
    continuous map

    .. MATH::

        f: M \longrightarrow K.

    A scalar field on `M` is an element of the commutative algebra
    `C^0(M)` (see
    :class:`~sage.manifolds.scalarfield_algebra.ScalarFieldAlgebra`).

    INPUT:

    - ``parent`` -- the algebra of scalar fields containing the scalar field
      (must be an instance of class
      :class:`~sage.manifolds.scalarfield_algebra.ScalarFieldAlgebra`)

    - ``coord_expression`` -- (default: ``None``) coordinate expression(s) of
      the scalar field; this can be either

      * a dictionary of coordinate expressions in various charts on
        the domain, with the charts as keys;
      * a single coordinate expression; if the argument ``chart`` is
        ``'all'``, this expression is set to all the charts defined
        on the open set; otherwise, the expression is set in the
        specific chart provided by the argument ``chart``

    - ``chart`` -- (default: ``None``) chart defining the coordinates used
      in ``coord_expression`` when the latter is a single coordinate
      expression; if none is provided (default), the default chart of the
      open set is assumed. If ``chart=='all'``, ``coord_expression`` is
      assumed to be independent of the chart (constant scalar field).

    - ``name`` -- (default: ``None``) string; name (symbol) given to the
      scalar field

    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote the
      scalar field; if none is provided, the LaTeX symbol is set to ``name``

    If ``coord_expression`` is ``None`` or incomplete, coordinate
    expressions can be added after the creation of the object, by means of
    the methods :meth:`add_expr`, :meth:`add_expr_by_continuation` and
    :meth:`set_expr`.

    EXAMPLES:

    A scalar field on the 2-sphere::

        sage: M = Manifold(2, 'M', structure='topological') # the 2-dimensional sphere S^2
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
        Scalar field f on the 2-dimensional topological manifold M
        sage: f.display()
        f: M --> R
        on U: (x, y) |--> 1/(x^2 + y^2 + 1)
        on V: (u, v) |--> (u^2 + v^2)/(u^2 + v^2 + 1)

    For scalar fields defined by a single coordinate expression, the latter
    can be passed instead of the dictionary over the charts::

        sage: g = U.scalar_field(x*y, chart=c_xy, name='g') ; g
        Scalar field g on the Open subset U of the 2-dimensional topological
         manifold M

    The above is indeed equivalent to::

        sage: g = U.scalar_field({c_xy: x*y}, name='g') ; g
        Scalar field g on the Open subset U of the 2-dimensional topological
         manifold M

    Since ``c_xy`` is the default chart of ``U``, the argument ``chart`` can
    be skipped::

        sage: g = U.scalar_field(x*y, name='g') ; g
        Scalar field g on the Open subset U of the 2-dimensional topological
         manifold M

    The scalar field `g` is defined on `U` and has an expression in terms of
    the coordinates `(u,v)` on `W=U\cap V`::

        sage: g.display()
        g: U --> R
           (x, y) |--> x*y
        on W: (u, v) |--> u*v/(u^4 + 2*u^2*v^2 + v^4)

    Scalar fields on `M` can also be declared with a single chart::

        sage: f = M.scalar_field(1/(1+x^2+y^2), chart=c_xy, name='f') ; f
        Scalar field f on the 2-dimensional topological manifold M

    Their definition must then be completed by providing the expressions on
    other charts, via the method :meth:`add_expr`, to get a global cover of
    the manifold::

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

    We may also use the method :meth:`add_expr_by_continuation` to complete
    the coordinate definition using the analytic continuation from domains in
    which charts overlap::

        sage: f = M.scalar_field(1/(1+x^2+y^2), chart=c_xy, name='f') ; f
        Scalar field f on the 2-dimensional topological manifold M
        sage: f.add_expr_by_continuation(c_uv, U.intersection(V))
        sage: f.display()
        f: M --> R
        on U: (x, y) |--> 1/(x^2 + y^2 + 1)
        on V: (u, v) |--> (u^2 + v^2)/(u^2 + v^2 + 1)

    A scalar field can also be defined by some unspecified function of the
    coordinates::

        sage: h = U.scalar_field(function('H')(x, y), name='h') ; h
        Scalar field h on the Open subset U of the 2-dimensional topological
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
    :meth:`expr`, which returns a symbolic expression::

        sage: f.expr(c_uv)
        (u^2 + v^2)/(u^2 + v^2 + 1)
        sage: type(f.expr(c_uv))
        <type 'sage.symbolic.expression.Expression'>

    The method :meth:`coord_function` returns instead a function of the
    chart coordinates, i.e. an instance of
    :class:`~sage.manifolds.coord_func.CoordFunction`::

        sage: f.coord_function(c_uv)
        (u^2 + v^2)/(u^2 + v^2 + 1)
        sage: type(f.coord_function(c_uv))
        <class 'sage.manifolds.coord_func_symb.CoordFunctionSymbRing_with_category.element_class'>
        sage: f.coord_function(c_uv).display()
        (u, v) |--> (u^2 + v^2)/(u^2 + v^2 + 1)

    The value returned by the method :meth:`expr` is actually the coordinate
    expression of the chart function::

        sage: f.expr(c_uv) is f.coord_function(c_uv).expr()
        True

    A constant scalar field is declared by setting the argument ``chart`` to
    ``'all'``::

        sage: c = M.scalar_field(2, chart='all', name='c') ; c
        Scalar field c on the 2-dimensional topological manifold M
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
        Scalar field c on the 2-dimensional topological manifold M
        sage: c.display()
        c: M --> R
        on U: (x, y) |--> a
        on V: (u, v) |--> a

    A special case of constant field is the zero scalar field::

        sage: zer = M.constant_scalar_field(0) ; zer
        Scalar field zero on the 2-dimensional topological manifold M
        sage: zer.display()
        zero: M --> R
        on U: (x, y) |--> 0
        on V: (u, v) |--> 0

    It can be obtained directly by means of the function
    :meth:`~sage.manifolds.manifold.TopologicalManifold.zero_scalar_field`::

        sage: zer is M.zero_scalar_field()
        True

    A third way is to get it as the zero element of the algebra `C^0(M)`
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
        Scalar field sqrt(f) on the 2-dimensional topological manifold M
        sage: sqrt(f).display()
        sqrt(f): M --> R
        on U: (x, y) |--> 1/sqrt(x^2 + y^2 + 1)
        on V: (u, v) |--> sqrt(u^2 + v^2)/sqrt(u^2 + v^2 + 1)

    ::

        sage: tan(f)
        Scalar field tan(f) on the 2-dimensional topological manifold M
        sage: tan(f).display()
        tan(f): M --> R
        on U: (x, y) |--> sin(1/(x^2 + y^2 + 1))/cos(1/(x^2 + y^2 + 1))
        on V: (u, v) |--> sin((u^2 + v^2)/(u^2 + v^2 + 1))/cos((u^2 + v^2)/(u^2 + v^2 + 1))

    .. RUBRIC:: Arithmetics of scalar fields

    Scalar fields on `M` (resp. `U`) belong to the algebra `C^0(M)`
    (resp. `C^0(U)`)::

        sage: f.parent()
        Algebra of scalar fields on the 2-dimensional topological manifold M
        sage: f.parent() is M.scalar_field_algebra()
        True
        sage: g.parent()
        Algebra of scalar fields on the Open subset U of the 2-dimensional
         topological manifold M
        sage: g.parent() is U.scalar_field_algebra()
        True

    Consequently, scalar fields can be added::

        sage: s = f + c ; s
        Scalar field f+c on the 2-dimensional topological manifold M
        sage: s.display()
        f+c: M --> R
        on U: (x, y) |--> (a*x^2 + a*y^2 + a + 1)/(x^2 + y^2 + 1)
        on V: (u, v) |--> ((a + 1)*u^2 + (a + 1)*v^2 + a)/(u^2 + v^2 + 1)

    and subtracted::

        sage: s = f - c ; s
        Scalar field f-c on the 2-dimensional topological manifold M
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
        Scalar field on the 2-dimensional topological manifold M
        sage: s.display()
        M --> R
        on U: (x, y) |--> (x^2 + y^2 + 2)/(x^2 + y^2 + 1)
        on V: (u, v) |--> (2*u^2 + 2*v^2 + 1)/(u^2 + v^2 + 1)
        sage: (f+1)-1 == f
        True

    The number can represented by a symbolic variable::

        sage: s = a + f ; s
        Scalar field on the 2-dimensional topological manifold M
        sage: s == c + f
        True

    However if the symbolic variable is a chart coordinate, the addition
    is performed only on the chart domain::

        sage: s = f + x; s
        Scalar field on the 2-dimensional topological manifold M
        sage: s.display()
        M --> R
        on U: (x, y) |--> (x^3 + x*y^2 + x + 1)/(x^2 + y^2 + 1)
        sage: s = f + u; s
        Scalar field on the 2-dimensional topological manifold M
        sage: s.display()
        M --> R
        on V: (u, v) |--> (u^3 + (u + 1)*v^2 + u^2 + u)/(u^2 + v^2 + 1)

    The addition of two scalar fields with different domains is possible if
    the domain of one of them is a subset of the domain of the other; the
    domain of the result is then this subset::

        sage: f.domain()
        2-dimensional topological manifold M
        sage: g.domain()
        Open subset U of the 2-dimensional topological manifold M
        sage: s = f + g ; s
        Scalar field on the Open subset U of the 2-dimensional topological
         manifold M
        sage: s.domain()
        Open subset U of the 2-dimensional topological manifold M
        sage: s.display()
        U --> R
        (x, y) |--> (x*y^3 + (x^3 + x)*y + 1)/(x^2 + y^2 + 1)
        on W: (u, v) |--> (u^6 + 3*u^4*v^2 + 3*u^2*v^4 + v^6 + u*v^3
         + (u^3 + u)*v)/(u^6 + v^6 + (3*u^2 + 1)*v^4 + u^4 + (3*u^4 + 2*u^2)*v^2)

    The operation actually performed is `f|_U + g`::

        sage: s == f.restrict(U) + g
        True

    In Sage framework, the addition of `f` and `g` is permitted because
    there is a *coercion* of the parent of `f`, namely `C^0(M)`, to
    the parent of `g`, namely `C^0(U)` (see
    :class:`~sage.manifolds.scalarfield_algebra.ScalarFieldAlgebra`)::

        sage: CM = M.scalar_field_algebra()
        sage: CU = U.scalar_field_algebra()
        sage: CU.has_coerce_map_from(CM)
        True

    The coercion map is nothing but the restriction to domain `U`::

        sage: CU.coerce(f) == f.restrict(U)
        True

    Since the algebra `C^0(M)` is a vector space over `\RR`, scalar fields
    can be multiplied by a number, either an explicit one::

        sage: s = 2*f ; s
        Scalar field on the 2-dimensional topological manifold M
        sage: s.display()
        M --> R
        on U: (x, y) |--> 2/(x^2 + y^2 + 1)
        on V: (u, v) |--> 2*(u^2 + v^2)/(u^2 + v^2 + 1)

    or a symbolic one::

        sage: s = a*f ; s
        Scalar field on the 2-dimensional topological manifold M
        sage: s.display()
        M --> R
        on U: (x, y) |--> a/(x^2 + y^2 + 1)
        on V: (u, v) |--> (u^2 + v^2)*a/(u^2 + v^2 + 1)

    However, if the symbolic variable is a chart coordinate, the
    multiplication is performed only in the corresponding chart::

        sage: s = x*f; s
        Scalar field on the 2-dimensional topological manifold M
        sage: s.display()
        M --> R
        on U: (x, y) |--> x/(x^2 + y^2 + 1)
        sage: s = u*f; s
        Scalar field on the 2-dimensional topological manifold M
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

    The ring multiplication of the algebras `C^0(M)` and `C^0(U)`
    is the pointwise multiplication of functions::

        sage: s = f*f ; s
        Scalar field f*f on the 2-dimensional topological manifold M
        sage: s.display()
        f*f: M --> R
        on U: (x, y) |--> 1/(x^4 + y^4 + 2*(x^2 + 1)*y^2 + 2*x^2 + 1)
        on V: (u, v) |--> (u^4 + 2*u^2*v^2 + v^4)/(u^4 + v^4 + 2*(u^2 + 1)*v^2
         + 2*u^2 + 1)
        sage: s = g*h ; s
        Scalar field g*h on the Open subset U of the 2-dimensional topological
         manifold M
        sage: s.display()
        g*h: U --> R
           (x, y) |--> x*y*H(x, y)
        on W: (u, v) |--> u*v*H(u/(u^2 + v^2), v/(u^2 + v^2))/(u^4 + 2*u^2*v^2 + v^4)

    Thanks to the coercion `C^0(M) \to C^0(U)` mentioned above,
    it is possible to multiply a scalar field defined on `M` by a
    scalar field defined on `U`, the result being a scalar field
    defined on `U`::

        sage: f.domain(), g.domain()
        (2-dimensional topological manifold M,
         Open subset U of the 2-dimensional topological manifold M)
        sage: s = f*g ; s
        Scalar field on the Open subset U of the 2-dimensional topological
         manifold M
        sage: s.display()
        U --> R
        (x, y) |--> x*y/(x^2 + y^2 + 1)
        on W: (u, v) |--> u*v/(u^4 + v^4 + (2*u^2 + 1)*v^2 + u^2)
        sage: s == f.restrict(U)*g
        True

    Scalar fields can be divided (pointwise division)::

        sage: s = f/c ; s
        Scalar field f/c on the 2-dimensional topological manifold M
        sage: s.display()
        f/c: M --> R
        on U: (x, y) |--> 1/(a*x^2 + a*y^2 + a)
        on V: (u, v) |--> (u^2 + v^2)/(a*u^2 + a*v^2 + a)
        sage: s = g/h ; s
        Scalar field g/h on the Open subset U of the 2-dimensional topological
         manifold M
        sage: s.display()
        g/h: U --> R
           (x, y) |--> x*y/H(x, y)
        on W: (u, v) |--> u*v/((u^4 + 2*u^2*v^2 + v^4)*H(u/(u^2 + v^2), v/(u^2 + v^2)))
        sage: s = f/g ; s
        Scalar field on the Open subset U of the 2-dimensional topological
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
        Scalar field on the Open subset U of the 2-dimensional topological
         manifold M
        sage: s.display()
        U --> R
        (x, y) |--> x^2 + (x - 1)*y
        on W: (u, v) |--> -(v^3 - u^2 + (u^2 - u)*v)/(u^4 + 2*u^2*v^2 + v^4)

    ::

        sage: s = g*x ; s
        Scalar field on the Open subset U of the 2-dimensional topological
         manifold M
        sage: s.display()
        U --> R
        (x, y) |--> x^2*y
        on W: (u, v) |--> u^2*v/(u^6 + 3*u^4*v^2 + 3*u^2*v^4 + v^6)

    ::

        sage: s = g/x ; s
        Scalar field on the Open subset U of the 2-dimensional topological
         manifold M
        sage: s.display()
        U --> R
        (x, y) |--> y
        on W: (u, v) |--> v/(u^2 + v^2)
        sage: s = x/g ; s
        Scalar field on the Open subset U of the 2-dimensional topological
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

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x+y}, name='f') ; f
            Scalar field f on the 2-dimensional topological manifold M
            sage: from sage.manifolds.scalarfield import ScalarField
            sage: isinstance(f, ScalarField)
            True
            sage: f.parent()
            Algebra of scalar fields on the 2-dimensional topological
             manifold M
            sage: TestSuite(f).run()

        """
        CommutativeAlgebraElement.__init__(self, parent)
        domain = parent._domain
        self._domain = domain
        self._manifold = domain.manifold()
        self._is_zero = False # a priori, may be changed below or via
                              # method __nonzero__()
        self._name = name
        if latex_name is None:
            self._latex_name = self._name
        else:
            self._latex_name = latex_name
        self._express = {} # dict of coordinate expressions (CoordFunction
                           # instances) with charts as keys
        if coord_expression is not None:
            if isinstance(coord_expression, dict):
                for chart, expression in coord_expression.iteritems():
                    if isinstance(expression, CoordFunction):
                        self._express[chart] = expression
                    else:
                        self._express[chart] = chart.function(expression)
            elif isinstance(coord_expression, CoordFunction):
                self._express[coord_expression.chart()] = coord_expression
            else:
                if chart is None:
                    chart = self._domain.default_chart()
                if chart == 'all':
                    # coord_expression is the same in all charts (constant
                    # scalar field)
                    for ch in self._domain.atlas():
                        self._express[ch] = ch.function(coord_expression)
                else:
                    self._express[chart] = chart.function(coord_expression)
        self._init_derived()   # initialization of derived quantities

    ####### Required methods for an algebra element (beside arithmetic) #######

    def __nonzero__(self):
        r"""
        Return ``True`` if ``self`` is nonzero and ``False`` otherwise.

        This method is called by :meth:`is_zero()`.

        EXAMPLES:

        Tests on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: c_xy.<x,y> = M.chart()
            sage: f = M.scalar_field(x*y)
            sage: f.is_zero()
            False
            sage: f.set_expr(0)
            sage: f.is_zero()
            True
            sage: g = M.scalar_field(0)
            sage: g.is_zero()
            True
            sage: M.zero_scalar_field().is_zero()
            True

        """
        if self._is_zero:
            return False
        if not self._express:
            # undefined scalar field
            return True
        iszero = True
        for funct in itervalues(self._express):
            iszero = iszero and funct.is_zero()
        self._is_zero = iszero
        return not iszero

    def __eq__(self, other):
        r"""
        Comparison (equality) operator.

        INPUT:

        - ``other`` -- a scalar field (or something else)

        OUTPUT:

        - ``True`` if ``self`` is equal to ``other``, ``False`` otherwise

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x+y})
            sage: f == 1
            False
            sage: f == M.zero_scalar_field()
            False
            sage: g = M.scalar_field({X: x+y})
            sage: f == g
            True
            sage: h = M.scalar_field({X: 1})
            sage: h == M.one_scalar_field()
            True
            sage: h == 1
            True

        """
        if other is self:
            return True
        if not isinstance(other, ScalarField):
            # We try a conversion of other to a scalar field, except if
            # other is None (since this would generate an undefined scalar
            # field)
            if other is None:
                return False
            try:
                other = self.parent()(other)  # conversion to a scalar field
            except TypeError:
                return False
        if other._domain != self._domain:
            return False
        if other.is_zero():
            return self.is_zero()
        com_charts = self.common_charts(other)
        if com_charts is None:
            raise ValueError("no common chart for the comparison")
        resu = True
        for chart in com_charts:
            resu = resu and (self._express[chart] == other._express[chart])
        return resu

    def __ne__(self, other):
        r"""
        Non-equality operator.

        INPUT:

        - ``other`` -- a scalar field

        OUTPUT:

        - ``True`` if ``self`` differs from ``other``, ``False`` otherwise

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x+y})
            sage: f != 1
            True
            sage: f != M.zero_scalar_field()
            True
            sage: g = M.scalar_field({X: x+y})
            sage: f != g
            False

        """
        return not (self == other)

    ####### End of required methods for an algebra element (beside arithmetic) #######

    def _init_derived(self):
        r"""
        Initialize the derived quantities.

        TEST::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x+y})
            sage: f._init_derived()

        """
        self._restrictions = {} # dict. of restrictions of self on subsets
                                # of self._domain, with the subsets as keys

    def _del_derived(self):
        r"""
        Delete the derived quantities.

        TEST::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x+y})
            sage: U = M.open_subset('U', coord_def={X: x>0})
            sage: f.restrict(U)
            Scalar field on the Open subset U of the 2-dimensional topological
             manifold M
            sage: f._restrictions
            {Open subset U of the 2-dimensional topological manifold M:
             Scalar field on the Open subset U of the 2-dimensional topological
             manifold M}
            sage: f._del_derived()
            sage: f._restrictions  # restrictions are derived quantities
            {}

        """
        self._restrictions.clear()

    def _repr_(self):
        r"""
        String representation of the object.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x+y})
            sage: f._repr_()
            'Scalar field on the 2-dimensional topological manifold M'
            sage: f = M.scalar_field({X: x+y}, name='f')
            sage: f._repr_()
            'Scalar field f on the 2-dimensional topological manifold M'
            sage: f
            Scalar field f on the 2-dimensional topological manifold M

        """
        description = "Scalar field"
        if self._name is not None:
            description += " " + self._name
        description += " on the {}".format(self._domain)
        return description

    def _latex_(self):
        r"""
        LaTeX representation of the object.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x+y})
            sage: f._latex_()
            '\\mbox{Scalar field on the 2-dimensional topological manifold M}'
            sage: f = M.scalar_field({X: x+y}, name='f')
            sage: f._latex_()
            'f'
            sage: f = M.scalar_field({X: x+y}, name='f', latex_name=r'\Phi')
            sage: f._latex_()
            '\\Phi'
            sage: latex(f)
            \Phi

        """
        if self._latex_name is None:
            return r'\mbox{' + str(self) + r'}'
        else:
           return self._latex_name

    def set_name(self, name=None, latex_name=None):
        r"""
        Set (or change) the text name and LaTeX name of the scalar field.

        INPUT:

        - ``name`` -- (string; default: ``None``) name given to the scalar
          field
        - ``latex_name`` -- (string; default: ``None``) LaTeX symbol to denote
          the scalar field; if ``None`` while ``name`` is provided, the LaTeX
          symbol is set to ``name``

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x+y})
            sage: f = M.scalar_field({X: x+y}); f
            Scalar field on the 2-dimensional topological manifold M
            sage: f.set_name('f'); f
            Scalar field f on the 2-dimensional topological manifold M
            sage: latex(f)
            f
            sage: f.set_name('f', latex_name=r'\Phi'); f
            Scalar field f on the 2-dimensional topological manifold M
            sage: latex(f)
            \Phi

        """
        if name is not None:
            self._name = name
            if latex_name is None:
                self._latex_name = self._name
        if latex_name is not None:
            self._latex_name = latex_name

    def domain(self):
        r"""
        Return the open subset on which the scalar field is defined.

        OUTPUT:

        - instance of class
          :class:`~sage.manifolds.manifold.TopologicalManifold`
          representing the manifold's open subset on which the
          scalar field is defined

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: c_xy.<x,y> = M.chart()
            sage: f = M.scalar_field(x+2*y)
            sage: f.domain()
            2-dimensional topological manifold M
            sage: U = M.open_subset('U', coord_def={c_xy: x<0})
            sage: g = f.restrict(U)
            sage: g.domain()
            Open subset U of the 2-dimensional topological manifold M

        """
        return self._domain

    def copy(self):
        r"""
        Return an exact copy of the scalar field.

        EXAMPLES:

        Copy on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: c_xy.<x,y> = M.chart()
            sage: f = M.scalar_field(x*y^2)
            sage: g = f.copy()
            sage: type(g)
            <class 'sage.manifolds.scalarfield.ScalarFieldAlgebra_with_category.element_class'>
            sage: g.expr()
            x*y^2
            sage: g == f
            True
            sage: g is f
            False

        """
        result = type(self)(self.parent(), name=self._name,
                            latex_name=self._latex_name)
        for chart, funct in self._express.iteritems():
            result._express[chart] = funct.copy()
        return result

    def coord_function(self, chart=None, from_chart=None):
        r"""
        Return the function of the coordinates representing the scalar field
        in a given chart.

        INPUT:

        - ``chart`` -- (default: ``None``) chart with respect to which the
          coordinate expression is to be returned; if ``None``, the
          default chart of the scalar field's domain will be used
        - ``from_chart`` -- (default: ``None``) chart from which the
          required expression is computed if it is not known already in the
          chart ``chart``; if ``None``, a chart is picked in the known
          expressions

        OUTPUT:

        - instance of :class:`~sage.manifolds.coord_func.CoordFunction`
          representing the coordinate function of the scalar field in the
          given chart

        EXAMPLES:

        Coordinate function on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: c_xy.<x,y> = M.chart()
            sage: f = M.scalar_field(x*y^2)
            sage: f.coord_function()
            x*y^2
            sage: f.coord_function(c_xy)  # equivalent form (since c_xy is the default chart)
            x*y^2
            sage: type(f.coord_function())
            <class 'sage.manifolds.coord_func_symb.CoordFunctionSymbRing_with_category.element_class'>

        Expression via a change of coordinates::

            sage: c_uv.<u,v> = M.chart()
            sage: c_uv.transition_map(c_xy, [u+v, u-v])
            Change of coordinates from Chart (M, (u, v)) to Chart (M, (x, y))
            sage: f._express # at this stage, f is expressed only in terms of (x,y) coordinates
            {Chart (M, (x, y)): x*y^2}
            sage: f.coord_function(c_uv) # forces the computation of the expression of f in terms of (u,v) coordinates
            u^3 - u^2*v - u*v^2 + v^3
            sage: f.coord_function(c_uv) == (u+v)*(u-v)^2  # check
            True
            sage: f._express  # random (dict. output); f has now 2 coordinate expressions:
            {Chart (M, (x, y)): x*y^2, Chart (M, (u, v)): u^3 - u^2*v - u*v^2 + v^3}

        Usage in a physical context (simple Lorentz transformation - boost in
        ``x`` direction, with relative velocity ``v`` between ``o1``
        and ``o2`` frames)::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: o1.<t,x> = M.chart()
            sage: o2.<T,X> = M.chart()
            sage: f = M.scalar_field(x^2 - t^2)
            sage: f.coord_function(o1)
            -t^2 + x^2
            sage: v = var('v'); gam = 1/sqrt(1-v^2)
            sage: o2.transition_map(o1, [gam*(T - v*X), gam*(X - v*T)])
            Change of coordinates from Chart (M, (T, X)) to Chart (M, (t, x))
            sage: f.coord_function(o2)
            -T^2 + X^2

        """
        if chart is None:
            chart = self._domain._def_chart
        else:
            if chart not in self._domain._atlas:
                raise ValueError("the {} is not a chart ".format(chart) +
                                 "defined on the {}".format(self._domain))
        if chart not in self._express:
            # Check whether chart corresponds to a subchart of a chart
            # where the expression of self is known:
            for known_chart in self._express:
                if chart in known_chart._subcharts:
                    new_expr = self._express[known_chart].expr()
                    self._express[chart] = chart.function(new_expr)
                    return self._express[chart]
            # If this point is reached, the expression must be computed
            # from that in the chart from_chart, by means of a
            # change-of-coordinates formula:
            if from_chart is None:
                # from_chart in searched among the charts of known expressions
                # and subcharts of them
                known_express = self._express.copy()
                found = False
                for kchart in known_express:
                    for skchart in kchart._subcharts:
                        if (chart, skchart) in self._domain._coord_changes:
                            from_chart = skchart
                            found = True
                            if skchart not in self._express:
                                self._express[skchart] = skchart.function(
                                                  self._express[kchart].expr())
                            break
                    if found:
                        break
                if not found:
                    raise ValueError("no starting chart could be found to " +
                                     "compute the expression in the {}".format(chart))
            change = self._domain._coord_changes[(chart, from_chart)]
            # old coordinates expressed in terms of the new ones:
            coords = [ change._transf._functions[i]._express
                       for i in range(self._manifold.dim()) ]
            new_expr = self._express[from_chart](*coords)
            self._express[chart] = chart.function(new_expr)
            self._del_derived()
        return self._express[chart]

    def function_chart(self, chart=None, from_chart=None):
        r"""
        Deprecated.

        Use :meth:`coord_function` instead.

        EXAMPLE::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: c_xy.<x,y> = M.chart()
            sage: f = M.scalar_field(x*y^2)
            sage: fc = f.function_chart()
            doctest:...: DeprecationWarning: Use coord_function() instead.
            See http://trac.sagemath.org/18640 for details.
            sage: fc
            x*y^2

        """
        from sage.misc.superseded import deprecation
        deprecation(18640, 'Use coord_function() instead.')
        return self.coord_function(chart=chart, from_chart=from_chart)

    def expr(self, chart=None, from_chart=None):
        r"""
        Return the coordinate expression of the scalar field in a given
        chart.

        INPUT:

        - ``chart`` -- (default: ``None``) chart with respect to which the
          coordinate expression is required; if ``None``, the default
          chart of the scalar field's domain will be used
        - ``from_chart`` -- (default: ``None``) chart from which the
          required expression is computed if it is not known already in the
          chart ``chart``; if ``None``, a chart is picked in ``self._express``

        OUTPUT:

        - symbolic expression representing the coordinate
          expression of the scalar field in the given chart.

        EXAMPLES:

        Expression of a scalar field on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: c_xy.<x,y> = M.chart()
            sage: f = M.scalar_field(x*y^2)
            sage: f.expr()
            x*y^2
            sage: f.expr(c_xy)  # equivalent form (since c_xy is the default chart)
            x*y^2
            sage: type(f.expr())
            <type 'sage.symbolic.expression.Expression'>

        Expression via a change of coordinates::

            sage: c_uv.<u,v> = M.chart()
            sage: c_uv.transition_map(c_xy, [u+v, u-v])
            Change of coordinates from Chart (M, (u, v)) to Chart (M, (x, y))
            sage: f._express # at this stage, f is expressed only in terms of (x,y) coordinates
            {Chart (M, (x, y)): x*y^2}
            sage: f.expr(c_uv) # forces the computation of the expression of f in terms of (u,v) coordinates
            u^3 - u^2*v - u*v^2 + v^3
            sage: bool( f.expr(c_uv) == (u+v)*(u-v)^2 ) # check
            True
            sage: f._express  # random (dict. output); f has now 2 coordinate expressions:
            {Chart (M, (x, y)): x*y^2, Chart (M, (u, v)): u^3 - u^2*v - u*v^2 + v^3}

        """
        return self.coord_function(chart, from_chart)._express

    def set_expr(self, coord_expression, chart=None):
        r"""
        Set the coordinate expression of the scalar field.

        The expressions with respect to other charts are deleted, in order to
        avoid any inconsistency. To keep them, use :meth:`add_expr` instead.

        INPUT:

        - ``coord_expression`` -- coordinate expression of the scalar field
        - ``chart`` -- (default: ``None``) chart in which ``coord_expression``
          is defined; if ``None``, the default chart of the scalar field's
          domain is assumed

        EXAMPLES:

        Setting scalar field expressions on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: c_xy.<x,y> = M.chart()
            sage: f = M.scalar_field(x^2 + 2*x*y +1)
            sage: f._express
            {Chart (M, (x, y)): x^2 + 2*x*y + 1}
            sage: f.set_expr(3*y)
            sage: f._express  # the (x,y) expression has been changed:
            {Chart (M, (x, y)): 3*y}
            sage: c_uv.<u,v> = M.chart()
            sage: f.set_expr(cos(u)-sin(v), c_uv)
            sage: f._express # the (x,y) expression has been lost:
            {Chart (M, (u, v)): cos(u) - sin(v)}
            sage: f.set_expr(3*y)
            sage: f._express # the (u,v) expression has been lost:
            {Chart (M, (x, y)): 3*y}

        """
        if chart is None:
            chart = self._domain._def_chart
        self._is_zero = False # a priori
        self._express.clear()
        self._express[chart] = chart.function(coord_expression)
        self._del_derived()

    def add_expr(self, coord_expression, chart=None):
        r"""
        Add some coordinate expression to the scalar field.

        The previous expressions with respect to other charts are kept. To
        clear them, use :meth:`set_expr` instead.

        INPUT:

        - ``coord_expression`` -- coordinate expression of the scalar field
        - ``chart`` -- (default: ``None``) chart in which ``coord_expression``
          is defined; if ``None``, the default chart of the scalar field's
          domain is assumed

        .. WARNING::

            If the scalar field has already expressions in other charts, it
            is the user's responsibility to make sure that the expression
            to be added is consistent with them.

        EXAMPLES:

        Adding scalar field expressions on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: c_xy.<x,y> = M.chart()
            sage: f = M.scalar_field(x^2 + 2*x*y +1)
            sage: f._express
            {Chart (M, (x, y)): x^2 + 2*x*y + 1}
            sage: f.add_expr(3*y)
            sage: f._express  # the (x,y) expression has been changed:
            {Chart (M, (x, y)): 3*y}
            sage: c_uv.<u,v> = M.chart()
            sage: f.add_expr(cos(u)-sin(v), c_uv)
            sage: f._express # random (dict. output); f has now 2 expressions:
            {Chart (M, (x, y)): 3*y, Chart (M, (u, v)): cos(u) - sin(v)}

        """
        if chart is None:
            chart = self._domain._def_chart
        self._express[chart] = chart.function(coord_expression)
        self._is_zero = False # a priori
        self._del_derived()

    def add_expr_by_continuation(self, chart, subdomain):
        r"""
        Set coordinate expression in a chart by continuation of the
        coordinate expression in a subchart.

        The continuation is performed by demanding that the coordinate
        expression is identical to that in the restriction of the chart to
        a given subdomain.

        INPUT:

        - ``chart`` -- coordinate chart `(U,(x^i))` in which the expression of
          the scalar field is to set
        - ``subdomain`` -- open subset `V\subset U` in which the expression
          in terms of the restriction of the coordinate chart `(U,(x^i))` to
          `V` is already known or can be evaluated by a change of coordinates.

        EXAMPLES:

        Scalar field on the sphere `S^2`::

            sage: M = Manifold(2, 'S^2', structure='topological')
            sage: U = M.open_subset('U') ; V = M.open_subset('V') # the complement of resp. N pole and S pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: c_xy.<x,y> = U.chart() ; c_uv.<u,v> = V.chart() # stereographic coordinates
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                intersection_name='W', restrictions1= x^2+y^2!=0,
            ....:                restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: W =  U.intersection(V)  # S^2 minus the two poles
            sage: f = M.scalar_field(atan(x^2+y^2), chart=c_xy, name='f')

        The scalar field has been defined only on the domain covered by the
        chart ``c_xy``, i.e. `U`::

            sage: f.display()
            f: S^2 --> R
            on U: (x, y) |--> arctan(x^2 + y^2)

        We note that on `W = U \cap V`, the expression of `f` in terms of
        coordinates `(u,v)` can be deduced from that in the coordinates
        `(x,y)` thanks to the transition map between the two charts::

            sage: f.display(c_uv.restrict(W))
            f: S^2 --> R
            on W: (u, v) |--> arctan(1/(u^2 + v^2))

        We use this fact to extend the definition of `f` to the open
        subset `V`, covered by the chart ``c_uv``::

            sage: f.add_expr_by_continuation(c_uv, W)

        Then, `f` is known on the whole sphere::

            sage: f.display()
            f: S^2 --> R
            on U: (x, y) |--> arctan(x^2 + y^2)
            on V: (u, v) |--> arctan(1/(u^2 + v^2))

        """
        if not chart._domain.is_subset(self._domain):
            raise ValueError("the chart is not defined on a subset of " +
                             "the scalar field domain")
        schart = chart.restrict(subdomain)
        self._express[chart] = chart.function(self.expr(schart))
        self._is_zero = False # a priori
        self._del_derived()

    def display(self, chart=None):
        r"""
        Display the expression of the scalar field in a given chart.

        Without any argument, this function displays the expressions of the
        scalar field in all the charts defined on the scalar field's domain
        that are not restrictions of another chart to some subdomain
        (the "top charts").

        INPUT:

        - ``chart`` -- (default: ``None``) chart with respect to which
          the coordinate expression is to be displayed; if ``None``, the
          display is performed in all the top charts in which the
          coordinate expression is known

        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode).

        EXAMPLES:

        Various displays::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: c_xy.<x,y> = M.chart()
            sage: f = M.scalar_field(sqrt(x+1), name='f')
            sage: f.display()
            f: M --> R
               (x, y) |--> sqrt(x + 1)
            sage: latex(f.display())
            \begin{array}{llcl} f:& M & \longrightarrow & \mathbb{R} \\ & \left(x, y\right) & \longmapsto & \sqrt{x + 1} \end{array}
            sage: g = M.scalar_field(function('G')(x, y), name='g')
            sage: g.display()
            g: M --> R
               (x, y) |--> G(x, y)
            sage: latex(g.display())
            \begin{array}{llcl} g:& M & \longrightarrow & \mathbb{R} \\ & \left(x, y\right) & \longmapsto & G\left(x, y\right) \end{array}

        A shortcut of ``display()`` is ``disp()``::

            sage: f.disp()
            f: M --> R
               (x, y) |--> sqrt(x + 1)

        """
        from sage.misc.latex import latex
        from sage.tensor.modules.format_utilities import FormattedExpansion

        def _display_expression(self, chart, result):
            r"""
            Helper function for :meth:`display`.
            """
            try:
                expression = self.coord_function(chart)
                coords = chart[:]
                if len(coords) == 1:
                    coords = coords[0]
                if chart._domain == self._domain:
                    if self._name is not None:
                        result._txt += "   "
                    result._latex += " & "
                else:
                    result._txt += "on " + chart._domain._name + ": "
                    result._latex += r"\mbox{on}\ " + latex(chart._domain) + \
                                     r": & "
                result._txt += repr(coords) + " |--> " + repr(expression) + "\n"
                result._latex += latex(coords) + r"& \longmapsto & " + \
                                 latex(expression) + r"\\"
            except (TypeError, ValueError):
                pass

        # Name of the base field:
        field = self._domain.base_field()
        field_type = self._domain.base_field_type()
        if field_type == 'real':
            field_name = 'R'
            field_latex_name = r'\mathbb{R}'
        elif field_type == 'complex':
            field_name = 'C'
            field_latex_name = r'\mathbb{C}'
        else:
            field_name = str(field)
            field_latex_name = latex(field)
        #
        result = FormattedExpansion()
        if self._name is None:
            symbol = ""
        else:
            symbol = self._name + ": "
        result._txt = symbol + self._domain._name + " --> " + field_name + "\n"
        if self._latex_name is None:
            symbol = ""
        else:
            symbol = self._latex_name + ":"
        result._latex = r"\begin{array}{llcl} " + symbol + r"&" + \
                        latex(self._domain) + r"& \longrightarrow & " + \
                        field_latex_name + r" \\"
        if chart is None:
            for ch in self._domain._top_charts:
                _display_expression(self, ch, result)
        else:
            _display_expression(self, chart, result)
        result._txt = result._txt[:-1]
        result._latex = result._latex[:-2] + r"\end{array}"
        return result

    disp = display

    def restrict(self, subdomain):
        r"""
        Restriction of the scalar field to an open subset of its domain of
        definition.

        INPUT:

        - ``subdomain`` -- an open subset of the scalar field's domain

        OUTPUT:

        - instance of :class:`ScalarField` representing the restriction of
          the scalar field to ``subdomain``

        EXAMPLES:

        Restriction of a scalar field defined on `\RR^2` to the
        unit open disc::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()  # Cartesian coordinates
            sage: U = M.open_subset('U', coord_def={X: x^2+y^2 < 1}) # U unit open disc
            sage: f = M.scalar_field(cos(x*y), name='f')
            sage: f_U = f.restrict(U) ; f_U
            Scalar field f on the Open subset U of the 2-dimensional
             topological manifold M
            sage: f_U.display()
            f: U --> R
               (x, y) |--> cos(x*y)
            sage: f.parent()
            Algebra of scalar fields on the 2-dimensional topological
             manifold M
            sage: f_U.parent()
            Algebra of scalar fields on the Open subset U of the 2-dimensional
             topological manifold M

        The restriction to the whole domain is the identity::

            sage: f.restrict(M) is f
            True
            sage: f_U.restrict(U) is f_U
            True

        Restriction of the zero scalar field::

            sage: M.zero_scalar_field().restrict(U)
            Scalar field zero on the Open subset U of the 2-dimensional
             topological manifold M
            sage: M.zero_scalar_field().restrict(U) is U.zero_scalar_field()
            True

        """
        if subdomain == self._domain:
            return self
        if subdomain not in self._restrictions:
            if not subdomain.is_subset(self._domain):
                raise ValueError("the specified domain is not a subset of " +
                                 "the domain of definition of the scalar field")
            # Special case of the zero scalar field:
            if self._is_zero:
                return subdomain._zero_scalar_field
            # First one tries to get the restriction from a tighter domain:
            for dom, rst in self._restrictions.iteritems():
                if subdomain.is_subset(dom):
                    self._restrictions[subdomain] = rst.restrict(subdomain)
                    break
            else:
            # If this fails, the restriction is obtained via coercion
                resu = subdomain.scalar_field_algebra()(self)
                resu._name = self._name
                resu._latex_name = self._latex_name
                self._restrictions[subdomain] = resu
        return self._restrictions[subdomain]

    def common_charts(self, other):
        r"""
        Find common charts for the expressions of the scalar field and
        ``other``.

        INPUT:

        - ``other`` -- a scalar field

        OUTPUT:

        - list of common charts; if no common chart is found, ``None`` is
          returned (instead of an empty list)

        EXAMPLES:

        Search for common charts on a 2-dimensional manifold with 2
        overlapping domains::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: U = M.open_subset('U')
            sage: c_xy.<x,y> = U.chart()
            sage: V = M.open_subset('V')
            sage: c_uv.<u,v> = V.chart()
            sage: M.declare_union(U,V)   # M is the union of U and V
            sage: f = U.scalar_field(x^2)
            sage: g = M.scalar_field(x+y)
            sage: f.common_charts(g)
            [Chart (U, (x, y))]
            sage: g.add_expr(u, c_uv)
            sage: f._express
            {Chart (U, (x, y)): x^2}
            sage: g._express  # random (dictionary output)
            {Chart (U, (x, y)): x + y, Chart (V, (u, v)): u}
            sage: f.common_charts(g)
            [Chart (U, (x, y))]

        Common charts found as subcharts: the subcharts are introduced via
        a transition map between charts c_xy and c_uv on the intersecting
        subdomain `W = U\cap V`::

            sage: trans = c_xy.transition_map(c_uv, (x+y, x-y), 'W', x<0, u+v<0)
            sage: M.atlas()
            [Chart (U, (x, y)), Chart (V, (u, v)), Chart (W, (x, y)),
             Chart (W, (u, v))]
            sage: c_xy_W = M.atlas()[2]
            sage: c_uv_W = M.atlas()[3]
            sage: trans.inverse()
            Change of coordinates from Chart (W, (u, v)) to Chart (W, (x, y))
            sage: f.common_charts(g)
            [Chart (U, (x, y))]
            sage: f.expr(c_xy_W)
            x^2
            sage: f._express  # random (dictionary output)
            {Chart (U, (x, y)): x^2, Chart (W, (x, y)): x^2}
            sage: g._express  # random (dictionary output)
            {Chart (U, (x, y)): x + y, Chart (V, (u, v)): u}
            sage: g.common_charts(f)  # c_xy_W is not returned because it is subchart of 'xy'
            [Chart (U, (x, y))]
            sage: f.expr(c_uv_W)
            1/4*u^2 + 1/2*u*v + 1/4*v^2
            sage: f._express  # random (dictionary output)
            {Chart (U, (x, y)): x^2, Chart (W, (x, y)): x^2,
             Chart (W, (u, v)): 1/4*u^2 + 1/2*u*v + 1/4*v^2}
            sage: g._express  # random (dictionary output)
            {Chart (U, (x, y)): x + y, Chart (V, (u, v)): u}
            sage: f.common_charts(g)
            [Chart (U, (x, y)), Chart (W, (u, v))]
            sage: # the expressions have been updated on the subcharts
            sage: g._express #  random (dictionary output)
            {Chart (U, (x, y)): x + y, Chart (V, (u, v)): u,
             Chart (W, (u, v)): u}

        Common charts found by computing some coordinate changes::

            sage: W = U.intersection(V)
            sage: f = W.scalar_field(x^2, c_xy_W)
            sage: g = W.scalar_field(u+1, c_uv_W)
            sage: f._express
            {Chart (W, (x, y)): x^2}
            sage: g._express
            {Chart (W, (u, v)): u + 1}
            sage: f.common_charts(g)
            [Chart (W, (u, v)), Chart (W, (x, y))]
            sage: f._express # random (dictionary output)
            {Chart (W, (u, v)): 1/4*u^2 + 1/2*u*v + 1/4*v^2,
             Chart (W, (x, y)): x^2}
            sage: g._express # random (dictionary output)
            {Chart (W, (u, v)): u + 1, Chart (W, (x, y)): x + y + 1}

        """
        if not isinstance(other, ScalarField):
            raise TypeError("the second argument must be a scalar field")
        coord_changes = self._manifold._coord_changes
        resu = []
        #
        # 1/ Search for common charts among the existing expressions, i.e.
        #    without performing any expression transformation.
        #    -------------------------------------------------------------
        for chart1 in self._express:
            if chart1 in other._express:
                resu.append(chart1)
        # Search for a subchart:
        known_expr1 = self._express.copy()
        known_expr2 = other._express.copy()
        for chart1 in known_expr1:
            if chart1 not in resu:
                for chart2 in known_expr2:
                    if chart2 not in resu:
                        if chart2 in chart1._subcharts:
                            self.expr(chart2)
                            resu.append(chart2)
                        if chart1 in chart2._subcharts:
                            other.expr(chart1)
                            resu.append(chart1)
        #
        # 2/ Search for common charts via one expression transformation
        #    ----------------------------------------------------------
        for chart1 in known_expr1:
            if chart1 not in resu:
                for chart2 in known_expr2:
                    if chart2 not in resu:
                        if (chart1, chart2) in coord_changes:
                            self.coord_function(chart2, from_chart=chart1)
                            resu.append(chart2)
                        if (chart2, chart1) in coord_changes:
                            other.coord_function(chart1, from_chart=chart2)
                            resu.append(chart1)
        if resu == []:
            return None
        else:
            return resu

    def __call__(self, p, chart=None):
        r"""
        Compute the value of the scalar field at a given point.

        INPUT:

        - ``p`` -- point in the scalar field's domain
        - ``chart`` -- (default: ``None``) chart in which the coordinates
          of ``p`` are to be considered; if ``None``, a chart in which
          both ``p``'s coordinates and the expression of the scalar field
          are known is searched, starting from the default chart
          of ``self._domain``

        OUTPUT:

        - value at ``p``

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x+y}, name='f')
            sage: p = M((2,-5), name='p'); p
            Point p on the 2-dimensional topological manifold M
            sage: f.__call__(p)
            -3
            sage: f(p)
            -3
            sage: M.zero_scalar_field()(p)
            0
            sage: M.one_scalar_field()(p)
            1

        Example with a change of chart::

            sage: Y.<u,v> = M.chart()
            sage: X_to_Y = X.transition_map(Y, [x+y, x-y])
            sage: Y_to_X = X_to_Y.inverse()
            sage: g = M.scalar_field({Y: u*v}, name='g')
            sage: g(p)
            -21
            sage: p.coord(Y)
            (-3, 7)

        """
        #!# it should be "if p not in self_domain:" instead, but this test is
        # skipped for efficiency
        if p not in self._manifold:
            raise ValueError("the {} ".format(p) + "does not belong " +
                             "to the {}".format(self._manifold))
        if self._is_zero:
            return 0
        if chart is None:
            # A common chart is searched:
            def_chart = self._domain._def_chart
            if def_chart in p._coordinates and def_chart in self._express:
                chart = def_chart
            else:
                for chart_p in p._coordinates:
                    if chart_p in self._express:
                        chart = chart_p
                        break
        if chart is None:
            # A change of coordinates is attempted for p:
            for chart_s in self._express:
                try:
                    p.coord(chart_s)
                    chart = chart_s
                    break
                except ValueError:
                    pass
            else:
                # A change of coordinates is attempted on the scalar field
                # expressions:
                for chart_p in p._coordinates:
                    try:
                        self.coord_function(chart_p)
                        chart = chart_p
                        break
                    except (TypeError, ValueError):
                        pass
        if chart is None:
            raise ValueError("no common chart has been found to evaluate " +
                             "the action of {} on the {}".format(self, p))
        return self._express[chart](*(p._coordinates[chart]))

    def __pos__(self):
        r"""
        Unary plus operator.

        OUTPUT:

        - an exact copy of the scalar field

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x+y}, name='f')
            sage: g = f.__pos__(); g
            Scalar field +f on the 2-dimensional topological manifold M
            sage: g == f
            True

        """
        result = type(self)(self.parent())
        for chart in self._express:
            result._express[chart] = + self._express[chart]
        if self._name is not None:
            result._name = '+' + self._name
        if self._latex_name is not None:
            result._latex_name = '+' + self._latex_name
        return result

    def __neg__(self):
        r"""
        Unary minus operator.

        OUTPUT:

        - the negative of the scalar field

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x+y}, name='f')
            sage: g = f.__neg__(); g
            Scalar field -f on the 2-dimensional topological manifold M
            sage: g.display()
            -f: M --> R
               (x, y) |--> -x - y
            sage: g.__neg__() == f
            True

        """
        result = type(self)(self.parent())
        for chart in self._express:
            result._express[chart] = - self._express[chart]
        if self._name is not None:
            result._name = '-' + self._name
        if self._latex_name is not None:
            result._latex_name = '-' + self._latex_name
        return result


    #########  CommutativeAlgebraElement arithmetic operators ########

    def _add_(self, other):
        r"""
        Scalar field addition.

        INPUT:

        - ``other`` -- a scalar field (in the same algebra as ``self``)

        OUTPUT:

        - the scalar field resulting from the addition of ``self`` and
          ``other``

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x+y}, name='f')
            sage: g = M.scalar_field({X: x*y}, name='g')
            sage: s = f._add_(g); s
            Scalar field f+g on the 2-dimensional topological manifold M
            sage: s.display()
            f+g: M --> R
               (x, y) |--> (x + 1)*y + x
            sage: s == f+g
            True
            sage: f._add_(M.zero_scalar_field()) == f
            True

        """
        # Special cases:
        if self._is_zero:
            return other.copy()
        if other._is_zero:
            return self.copy()
        # Generic case:
        com_charts = self.common_charts(other)
        if com_charts is None:
            raise ValueError("no common chart for the addition")
        result = type(self)(self.parent())
        for chart in com_charts:
            # CoordFunction addition:
            result._express[chart] = self._express[chart] + other._express[chart]
        if result.is_zero():
            return self._domain.zero_scalar_field()
        if self._name is not None and other._name is not None:
            result._name = self._name + '+' + other._name
        if self._latex_name is not None and other._latex_name is not None:
            result._latex_name = self._latex_name + '+' + other._latex_name
        return result

    def _sub_(self, other):
        r"""
        Scalar field subtraction.

        INPUT:

        - ``other`` -- a scalar field (in the same algebra as ``self``)

        OUTPUT:

        - the scalar field resulting from the subtraction of ``other`` from
          ``self``

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x+y}, name='f')
            sage: g = M.scalar_field({X: x*y}, name='g')
            sage: s = f._sub_(g); s
            Scalar field f-g on the 2-dimensional topological manifold M
            sage: s.display()
            f-g: M --> R
               (x, y) |--> -(x - 1)*y + x
            sage: s == f-g
            True
            sage: f._sub_(M.zero_scalar_field()) == f
            True

        """
        # Special cases:
        if self._is_zero:
            return -other
        if other._is_zero:
            return self.copy()
        # Generic case:
        com_charts = self.common_charts(other)
        if com_charts is None:
            raise ValueError("no common chart for the subtraction")
        result = type(self)(self.parent())
        for chart in com_charts:
            # CoordFunction subtraction:
            result._express[chart] = self._express[chart] - other._express[chart]
        if result.is_zero():
            return self._domain.zero_scalar_field()
        if self._name is not None and other._name is not None:
            result._name = self._name + '-' + other._name
        if self._latex_name is not None and other._latex_name is not None:
            result._latex_name = self._latex_name + '-' + other._latex_name
        return result


    def _mul_(self, other):
        r"""
        Scalar field multiplication.

        INPUT:

        - ``other`` -- a scalar field (in the same algebra as ``self``)

        OUTPUT:

        - the scalar field resulting from the multiplication of ``self`` by
          ``other``

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x+y}, name='f')
            sage: g = M.scalar_field({X: x*y}, name='g')
            sage: s = f._mul_(g); s
            Scalar field f*g on the 2-dimensional topological manifold M
            sage: s.display()
            f*g: M --> R
               (x, y) |--> x^2*y + x*y^2
            sage: s == f*g
            True
            sage: f._mul_(M.zero_scalar_field()) == M.zero_scalar_field()
            True
            sage: f._mul_(M.one_scalar_field()) == f
            True

        """
        from sage.tensor.modules.format_utilities import format_mul_txt, \
                                                         format_mul_latex
        # Special cases:
        if self._is_zero or other._is_zero:
            return self._domain.zero_scalar_field()
        # Generic case:
        com_charts = self.common_charts(other)
        if com_charts is None:
            raise ValueError("no common chart for the multiplication")
        result = type(self)(self.parent())
        for chart in com_charts:
            # CoordFunction multiplication:
            result._express[chart] = self._express[chart] * other._express[chart]
        result._name = format_mul_txt(self._name, '*', other._name)
        result._latex_name = format_mul_latex(self._latex_name, ' ',
                                             other._latex_name)
        return result

    def _div_(self, other):
        r"""
        Scalar field division.

        INPUT:

        - ``other`` -- a scalar field (in the same algebra as self)

        OUTPUT:

        - the scalar field resulting from the division of ``self`` by
          ``other``

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x+y}, name='f')
            sage: g = M.scalar_field({X: x*y}, name='g')
            sage: s = f._div_(g); s
            Scalar field f/g on the 2-dimensional topological manifold M
            sage: s.display()
            f/g: M --> R
               (x, y) |--> (x + y)/(x*y)
            sage: s == f/g
            True
            sage: f._div_(M.zero_scalar_field())
            Traceback (most recent call last):
            ...
            ZeroDivisionError: division of a scalar field by zero

        """
        from sage.tensor.modules.format_utilities import format_mul_txt, \
                                                         format_mul_latex
        # Special cases:
        if other._is_zero:
            raise ZeroDivisionError("division of a scalar field by zero")
        if self._is_zero:
            return self._domain.zero_scalar_field()
        # Generic case:
        com_charts = self.common_charts(other)
        if com_charts is None:
            raise ValueError("no common chart for the division")
        result = type(self)(self.parent())
        for chart in com_charts:
            # CoordFunction division:
            result._express[chart] = self._express[chart] / other._express[chart]
        result._name = format_mul_txt(self._name, '/', other._name)
        result._latex_name = format_mul_latex(self._latex_name, '/',
                                             other._latex_name)
        return result

    def _lmul_(self, number):
        r"""
        Scalar multiplication operator: return ``number * self`` or
        ``self * number``.

        This differs from ``_mul_(self, other)`` by the fact that ``number``
        is not assumed to be a scalar field defined on the same domain as
        ``self``, contrary to ``other`` in ``_mul_(self, other)``. In
        practice, ``number`` is a an element of the field on which the
        scalar field algebra is defined.

        INPUT:

        - ``number`` -- an element of the ring on which the scalar field
          algebra is defined; this should be an element of the topological
          field on which the manifold is constructed (possibly represented
          by a symbolic expression)

        OUTPUT:

        - the scalar field ``number * self``

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x+y}, name='f')
            sage: s = f._lmul_(2); s
            Scalar field on the 2-dimensional topological manifold M
            sage: s.display()
            M --> R
            (x, y) |--> 2*x + 2*y
            sage: s == 2 * f
            True
            sage: s == f * 2
            True
            sage: f._lmul_(pi).display()
            M --> R
            (x, y) |--> pi*(x + y)
            sage: f._lmul_(pi) == pi*f
            True
            sage: f._lmul_(0) == M.zero_scalar_field()
            True
            sage: f._lmul_(1) == f
            True

        """
        if number == 0:
            return self.parent().zero()
        result = type(self)(self.parent())
        if isinstance(number, Expression):
            var = number.variables()  # possible symbolic variables in number
            if var:
                # There are symbolic variables in number
                # Are any of them a chart coordinate ?
                chart_var = False
                for chart in self._express:
                    if any(s in chart[:] for s in var):
                        chart_var = True
                        break
                if chart_var:
                    # Some symbolic variables in number are chart coordinates
                    for chart, expr in self._express.iteritems():
                        # The multiplication is performed only if
                        # either
                        # (i) all the symbolic variables in number are
                        # coordinates of this chart
                        # or (ii) no symbolic variable in number belongs to a
                        # different chart
                        chart_coords = chart[:]
                        var_not_in_chart = [s for s in var
                                            if not s in chart_coords]
                        any_in_other_chart = False
                        if var_not_in_chart != []:
                            for other_chart in self._domain.atlas():
                                other_chart_coords = other_chart[:]
                                for s in var_not_in_chart:
                                    if s in other_chart_coords:
                                        any_in_other_chart = True
                                        break
                                if any_in_other_chart:
                                    break
                        if not any_in_other_chart:
                            result._express[chart] = number * expr
                    return result
        # General case: the multiplication is performed on all charts:
        for chart, expr in self._express.iteritems():
            result._express[chart] = number * expr
        return result

    #########  End of CommutativeAlgebraElement arithmetic operators ########

    def _function_name(self, func, func_latex, parentheses=True):
        r"""
        Helper function to set the symbol of a function applied to the
        scalar field.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x+y}, name='f', latex_name=r"\Phi")
            sage: f._function_name("cos", r"\cos")
            ('cos(f)', '\\cos\\left(\\Phi\\right)')
            sage: f._function_name("sqrt", r"\sqrt", parentheses=False)
            ('sqrt(f)', '\\sqrt{\\Phi}')
            sage: f = M.scalar_field({X: x+y})  # no name given to f
            sage: f._function_name("cos", r"\cos")
            (None, None)

        """
        if self._name is None:
            name = None
        else:
            name = func + "(" + self._name + ")"
        if self._latex_name is None:
            latex_name = None
        else:
            if parentheses:
                latex_name = func_latex + r"\left(" + self._latex_name + \
                             r"\right)"
            else:
                latex_name = func_latex + r"{" + self._latex_name + r"}"
        return name, latex_name

    def exp(self):
        r"""
        Exponential of the scalar field.

        OUTPUT:

        - the scalar field `\exp f`, where `f` is the current scalar field

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x+y}, name='f', latex_name=r"\Phi")
            sage: g = exp(f) ; g
            Scalar field exp(f) on the 2-dimensional topological manifold M
            sage: g.display()
            exp(f): M --> R
               (x, y) |--> e^(x + y)
            sage: latex(g)
            \exp\left(\Phi\right)

        Automatic simplifications occur::

            sage: f = M.scalar_field({X: 2*ln(1+x^2)}, name='f')
            sage: exp(f).display()
            exp(f): M --> R
               (x, y) |--> x^4 + 2*x^2 + 1

        The inverse function is :meth:`log`::

            sage: log(exp(f)) == f
            True

        Some tests::

            sage: exp(M.zero_scalar_field()) == M.constant_scalar_field(1)
            True
            sage: exp(M.constant_scalar_field(1)) == M.constant_scalar_field(e)
            True

        """
        name, latex_name = self._function_name("exp", r"\exp")
        resu = type(self)(self.parent(), name=name, latex_name=latex_name)
        for chart, func in self._express.iteritems():
            resu._express[chart] = func.exp()
        return resu

    def log(self):
        r"""
        Natural logarithm of the scalar field.

        OUTPUT:

        - the scalar field `\ln f`, where `f` is the current scalar field

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x+y}, name='f', latex_name=r"\Phi")
            sage: g = log(f) ; g
            Scalar field ln(f) on the 2-dimensional topological manifold M
            sage: g.display()
            ln(f): M --> R
               (x, y) |--> log(x + y)
            sage: latex(g)
            \ln\left(\Phi\right)

        The inverse function is :meth:`exp`::

            sage: exp(log(f)) == f
            True

        """
        name, latex_name = self._function_name("ln", r"\ln")
        resu = type(self)(self.parent(), name=name, latex_name=latex_name)
        for chart, func in self._express.iteritems():
            resu._express[chart] = func.log()
        return resu

    def __pow__(self, exponent):
        r"""
        The scalar field to a given power.

        INPUT:

        - ``exponent`` -- the exponent

        OUTPUT:

        - the scalar field `f^a`, where `f` is the current scalar field and
          `a` the exponent

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x+y}, name='f', latex_name=r'\Phi')
            sage: g = f.__pow__(pi) ; g
            Scalar field f^pi on the 2-dimensional topological manifold M
            sage: latex(g)
            {\Phi}^{ \pi }
            sage: g.display()
            f^pi: M --> R
               (x, y) |--> (x + y)^pi

        The global function ``pow`` can be used::

            sage: pow(f, pi) == f.__pow__(pi)
            True

        as well as the exponent notation::

            sage: f^pi == f.__pow__(pi)
            True

        Some checks::

            sage: pow(f, 2) == f*f
            True
            sage: pow(pow(f, 1/2), 2) == f
            True

        """
        from sage.misc.latex import latex
        if self._name is None:
            name = None
        else:
            name = self._name + "^{}".format(exponent)
        if self._latex_name is None:
            latex_name = None
        else:
            latex_name = r"{" + self._latex_name + r"}^{" + \
                         latex(exponent) + r"}"
        resu = type(self)(self.parent(), name=name, latex_name=latex_name)
        for chart, func in self._express.iteritems():
            resu._express[chart] = func.__pow__(exponent)
        return resu

    def sqrt(self):
        r"""
        Square root of the scalar field.

        OUTPUT:

        - the scalar field `\sqrt f`, where `f` is the current scalar field

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: 1+x^2+y^2}, name='f',
            ....:                    latex_name=r"\Phi")
            sage: g = sqrt(f) ; g
            Scalar field sqrt(f) on the 2-dimensional topological manifold M
            sage: latex(g)
            \sqrt{\Phi}
            sage: g.display()
            sqrt(f): M --> R
               (x, y) |--> sqrt(x^2 + y^2 + 1)

        Some tests::

            sage: g^2 == f
            True
            sage: sqrt(M.zero_scalar_field()) == M.zero_scalar_field()
            True

        """
        name, latex_name = self._function_name("sqrt", r"\sqrt",
                                               parentheses=False)
        resu = type(self)(self.parent(), name=name, latex_name=latex_name)
        for chart, func in self._express.iteritems():
            resu._express[chart] = func.sqrt()
        return resu

    def cos(self):
        r"""
        Cosine of the scalar field.

        OUTPUT:

        - the scalar field `\cos f`, where `f` is the current scalar field

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x*y}, name='f', latex_name=r"\Phi")
            sage: g = cos(f) ; g
            Scalar field cos(f) on the 2-dimensional topological manifold M
            sage: latex(g)
            \cos\left(\Phi\right)
            sage: g.display()
            cos(f): M --> R
               (x, y) |--> cos(x*y)

        Some tests::

            sage: cos(M.zero_scalar_field()) == M.constant_scalar_field(1)
            True
            sage: cos(M.constant_scalar_field(pi/2)) == M.zero_scalar_field()
            True

        """
        name, latex_name = self._function_name("cos", r"\cos")
        resu = type(self)(self.parent(), name=name, latex_name=latex_name)
        for chart, func in self._express.iteritems():
            resu._express[chart] = func.cos()
        return resu

    def sin(self):
        r"""
        Sine of the scalar field.

        OUTPUT:

        - the scalar field `\sin f`, where `f` is the current scalar field

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x*y}, name='f', latex_name=r"\Phi")
            sage: g = sin(f) ; g
            Scalar field sin(f) on the 2-dimensional topological manifold M
            sage: latex(g)
            \sin\left(\Phi\right)
            sage: g.display()
            sin(f): M --> R
               (x, y) |--> sin(x*y)

        Some tests::

            sage: sin(M.zero_scalar_field()) == M.zero_scalar_field()
            True
            sage: sin(M.constant_scalar_field(pi/2)) == M.constant_scalar_field(1)
            True

        """
        name, latex_name = self._function_name("sin", r"\sin")
        resu = type(self)(self.parent(), name=name, latex_name=latex_name)
        for chart, func in self._express.iteritems():
            resu._express[chart] = func.sin()
        return resu

    def tan(self):
        r"""
        Tangent of the scalar field.

        OUTPUT:

        - the scalar field `\tan f`, where `f` is the current scalar field

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x*y}, name='f', latex_name=r"\Phi")
            sage: g = tan(f) ; g
            Scalar field tan(f) on the 2-dimensional topological manifold M
            sage: latex(g)
            \tan\left(\Phi\right)
            sage: g.display()
            tan(f): M --> R
               (x, y) |--> sin(x*y)/cos(x*y)

        Some tests::

            sage: tan(f) == sin(f) / cos(f)
            True
            sage: tan(M.zero_scalar_field()) == M.zero_scalar_field()
            True
            sage: tan(M.constant_scalar_field(pi/4)) == M.constant_scalar_field(1)
            True

        """
        name, latex_name = self._function_name("tan", r"\tan")
        resu = type(self)(self.parent(), name=name, latex_name=latex_name)
        for chart, func in self._express.iteritems():
            resu._express[chart] = func.tan()
        return resu

    def arccos(self):
        r"""
        Arc cosine of the scalar field.

        OUTPUT:

        - the scalar field `\arccos f`, where `f` is the current scalar field

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x*y}, name='f', latex_name=r"\Phi")
            sage: g = arccos(f) ; g
            Scalar field arccos(f) on the 2-dimensional topological manifold M
            sage: latex(g)
            \arccos\left(\Phi\right)
            sage: g.display()
            arccos(f): M --> R
               (x, y) |--> arccos(x*y)

        The notation ``acos`` can be used as well::

            sage: acos(f)
            Scalar field arccos(f) on the 2-dimensional topological manifold M
            sage: acos(f) == g
            True

        Some tests::

            sage: cos(g) == f
            True
            sage: arccos(M.constant_scalar_field(1)) == M.zero_scalar_field()
            True
            sage: arccos(M.zero_scalar_field()) == M.constant_scalar_field(pi/2)
            True

        """
        name, latex_name = self._function_name("arccos", r"\arccos")
        resu = type(self)(self.parent(), name=name, latex_name=latex_name)
        for chart, func in self._express.iteritems():
            resu._express[chart] = func.arccos()
        return resu

    def arcsin(self):
        r"""
        Arc sine of the scalar field.

        OUTPUT:

        - the scalar field `\arcsin f`, where `f` is the current scalar field

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x*y}, name='f', latex_name=r"\Phi")
            sage: g = arcsin(f) ; g
            Scalar field arcsin(f) on the 2-dimensional topological manifold M
            sage: latex(g)
            \arcsin\left(\Phi\right)
            sage: g.display()
            arcsin(f): M --> R
               (x, y) |--> arcsin(x*y)

        The notation ``asin`` can be used as well::

            sage: asin(f)
            Scalar field arcsin(f) on the 2-dimensional topological manifold M
            sage: asin(f) == g
            True

        Some tests::

            sage: sin(g) == f
            True
            sage: arcsin(M.zero_scalar_field()) == M.zero_scalar_field()
            True
            sage: arcsin(M.constant_scalar_field(1)) == M.constant_scalar_field(pi/2)
            True

        """
        name, latex_name = self._function_name("arcsin", r"\arcsin")
        resu = type(self)(self.parent(), name=name, latex_name=latex_name)
        for chart, func in self._express.iteritems():
            resu._express[chart] = func.arcsin()
        return resu

    def arctan(self):
        r"""
        Arc tangent of the scalar field.

        OUTPUT:

        - the scalar field `\arctan f`, where `f` is the current scalar field

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x*y}, name='f', latex_name=r"\Phi")
            sage: g = arctan(f) ; g
            Scalar field arctan(f) on the 2-dimensional topological manifold M
            sage: latex(g)
            \arctan\left(\Phi\right)
            sage: g.display()
            arctan(f): M --> R
               (x, y) |--> arctan(x*y)

        The notation ``atan`` can be used as well::

            sage: atan(f)
            Scalar field arctan(f) on the 2-dimensional topological manifold M
            sage: atan(f) == g
            True

        Some tests::

            sage: tan(g) == f
            True
            sage: arctan(M.zero_scalar_field()) == M.zero_scalar_field()
            True
            sage: arctan(M.constant_scalar_field(1)) == M.constant_scalar_field(pi/4)
            True

        """
        name, latex_name = self._function_name("arctan", r"\arctan")
        resu = type(self)(self.parent(), name=name, latex_name=latex_name)
        for chart, func in self._express.iteritems():
            resu._express[chart] = func.arctan()
        return resu

    def cosh(self):
        r"""
        Hyperbolic cosine of the scalar field.

        OUTPUT:

        - the scalar field `\cosh f`, where `f` is the current scalar field

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x*y}, name='f', latex_name=r"\Phi")
            sage: g = cosh(f) ; g
            Scalar field cosh(f) on the 2-dimensional topological manifold M
            sage: latex(g)
            \cosh\left(\Phi\right)
            sage: g.display()
            cosh(f): M --> R
               (x, y) |--> cosh(x*y)

        Some test::

            sage: cosh(M.zero_scalar_field()) == M.constant_scalar_field(1)
            True

        """
        name, latex_name = self._function_name("cosh", r"\cosh")
        resu = type(self)(self.parent(), name=name, latex_name=latex_name)
        for chart, func in self._express.iteritems():
            resu._express[chart] = func.cosh()
        return resu

    def sinh(self):
        r"""
        Hyperbolic sine of the scalar field.

        OUTPUT:

        - the scalar field `\sinh f`, where `f` is the current scalar field

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x*y}, name='f', latex_name=r"\Phi")
            sage: g = sinh(f) ; g
            Scalar field sinh(f) on the 2-dimensional topological manifold M
            sage: latex(g)
            \sinh\left(\Phi\right)
            sage: g.display()
            sinh(f): M --> R
               (x, y) |--> sinh(x*y)

        Some test::

            sage: sinh(M.zero_scalar_field()) == M.zero_scalar_field()
            True

        """
        name, latex_name = self._function_name("sinh", r"\sinh")
        resu = type(self)(self.parent(), name=name, latex_name=latex_name)
        for chart, func in self._express.iteritems():
            resu._express[chart] = func.sinh()
        return resu

    def tanh(self):
        r"""
        Hyperbolic tangent of the scalar field.

        OUTPUT:

        - the scalar field `\tanh f`, where `f` is the current scalar field

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x*y}, name='f', latex_name=r"\Phi")
            sage: g = tanh(f) ; g
            Scalar field tanh(f) on the 2-dimensional topological manifold M
            sage: latex(g)
            \tanh\left(\Phi\right)
            sage: g.display()
            tanh(f): M --> R
               (x, y) |--> sinh(x*y)/cosh(x*y)

        Some tests::

            sage: tanh(f) == sinh(f) / cosh(f)
            True
            sage: tanh(M.zero_scalar_field()) == M.zero_scalar_field()
            True

        """
        name, latex_name = self._function_name("tanh", r"\tanh")
        resu = type(self)(self.parent(), name=name, latex_name=latex_name)
        for chart, func in self._express.iteritems():
            resu._express[chart] = func.tanh()
        return resu

    def arccosh(self):
        r"""
        Inverse hyperbolic cosine of the scalar field.

        OUTPUT:

        - the scalar field `\mathrm{arcosh}\, f`, where `f` is the current
          scalar field

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x*y}, name='f', latex_name=r"\Phi")
            sage: g = arccosh(f) ; g
            Scalar field arccosh(f) on the 2-dimensional topological manifold M
            sage: latex(g)
            \,\mathrm{arcosh}\left(\Phi\right)
            sage: g.display()
            arccosh(f): M --> R
               (x, y) |--> arccosh(x*y)

        The notation ``acosh`` can be used as well::

            sage: acosh(f)
            Scalar field arccosh(f) on the 2-dimensional topological manifold M
            sage: acosh(f) == g
            True

        Some tests::

            sage: cosh(g) == f
            True
            sage: arccosh(M.constant_scalar_field(1)) == M.zero_scalar_field()
            True

        """
        name, latex_name = self._function_name("arccosh", r"\,\mathrm{arcosh}")
        resu = type(self)(self.parent(), name=name, latex_name=latex_name)
        for chart, func in self._express.iteritems():
            resu._express[chart] = func.arccosh()
        return resu

    def arcsinh(self):
        r"""
        Inverse hyperbolic sine of the scalar field.

        OUTPUT:

        - the scalar field `\mathrm{arsinh}\, f`, where `f` is the current
          scalar field

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x*y}, name='f', latex_name=r"\Phi")
            sage: g = arcsinh(f) ; g
            Scalar field arcsinh(f) on the 2-dimensional topological manifold M
            sage: latex(g)
            \,\mathrm{arsinh}\left(\Phi\right)
            sage: g.display()
            arcsinh(f): M --> R
               (x, y) |--> arcsinh(x*y)

        The notation ``asinh`` can be used as well::

            sage: asinh(f)
            Scalar field arcsinh(f) on the 2-dimensional topological manifold M
            sage: asinh(f) == g
            True

        Some tests::

            sage: sinh(g) == f
            True
            sage: arcsinh(M.zero_scalar_field()) == M.zero_scalar_field()
            True

        """
        name, latex_name = self._function_name("arcsinh", r"\,\mathrm{arsinh}")
        resu = type(self)(self.parent(), name=name, latex_name=latex_name)
        for chart, func in self._express.iteritems():
            resu._express[chart] = func.arcsinh()
        return resu

    def arctanh(self):
        r"""
        Inverse hyperbolic tangent of the scalar field.

        OUTPUT:

        - the scalar field `\mathrm{artanh}\, f`, where `f` is the current
          scalar field

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field({X: x*y}, name='f', latex_name=r"\Phi")
            sage: g = arctanh(f) ; g
            Scalar field arctanh(f) on the 2-dimensional topological manifold M
            sage: latex(g)
            \,\mathrm{artanh}\left(\Phi\right)
            sage: g.display()
            arctanh(f): M --> R
               (x, y) |--> arctanh(x*y)

        The notation ``atanh`` can be used as well::

            sage: atanh(f)
            Scalar field arctanh(f) on the 2-dimensional topological manifold M
            sage: atanh(f) == g
            True

        Some tests::

            sage: tanh(g) == f
            True
            sage: arctanh(M.zero_scalar_field()) == M.zero_scalar_field()
            True
            sage: arctanh(M.constant_scalar_field(1/2)) == M.constant_scalar_field(log(3)/2)
            True

        """
        name, latex_name = self._function_name("arctanh", r"\,\mathrm{artanh}")
        resu = type(self)(self.parent(), name=name, latex_name=latex_name)
        for chart, func in self._express.iteritems():
            resu._express[chart] = func.arctanh()
        return resu

