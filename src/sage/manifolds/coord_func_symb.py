r"""
Symbolic Coordinate Functions

In the context of a topological manifold `M` over a topological field `K`,
a *coordinate function*  is a function from a chart codomain
to `K`. In other words, a coordinate function is a `K`-valued function of
the coordinates associated to some chart.

More precisely, let `(U, \varphi)` be a chart on `M`, i.e. `U` is an open
subset of `M` and `\varphi: U \to V \subset K^n` is a homeomorphism
from `U` to an open subset `V` of `K^n`. A *coordinate function associated
to the chart* `(U, \varphi)` is a function

.. MATH::

    \begin{array}{cccc}
        f:&  V \subset K^n & \longrightarrow & K \\
          &  (x^1, \ldots, x^n) & \longmapsto & f(x^1, \ldots, x^n)
    \end{array}

This module implements symbolic coordinate functions via the class
:class:`CoordFunctionSymb`.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013-2015) : initial version
- Travis Scrimshaw (2016) : make coordinate functions elements of
  :class:`CoordFunctionSymbRing`.

"""
#*****************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#       Copyright (C) 2016 Travis Scrimshaw <tscrimsh@umn.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.symbolic.ring import SR
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.commutative_algebras import CommutativeAlgebras
from sage.manifolds.coord_func import CoordFunction
from sage.manifolds.utilities import (ExpressionNice, simplify_chain_real,
                                      simplify_chain_generic)

class CoordFunctionSymb(CoordFunction):
    r"""
    Coordinate function with symbolic representation.

    If `(U, \varphi)` is a chart on a topological manifold `M` of
    dimension `n` over a topological field `K`,  a *coordinate function*
    associated to `(U, \varphi)` is a map

    .. MATH::

        \begin{array}{llcl}
        f:& V \subset K^n & \longrightarrow & K \\
          & (x^1, \ldots, x^n) & \longmapsto & f(x^1, \ldots, x^n),
        \end{array}

    where `V` is the codomain of `\varphi`. In other words, `f` is a
    `K`-valued function of the
    coordinates associated to the chart `(U, \varphi)`.

    INPUT:

    - ``parent`` -- the algebra of coordinate functions on the chart
      `(U, \varphi)`
    - ``expression`` -- a symbolic expression representing
      `f(x^1, \ldots, x^n)`, where `(x^1, \ldots, x^n)` are the
      coordinates of the chart `(U, \varphi)`

    EXAMPLES:

    A symbolic coordinate function associated with a 2-dimensional chart::

        sage: M = Manifold(2, 'M', structure='topological')
        sage: X.<x,y> = M.chart()
        sage: f = X.function(x^2+3*y+1)
        sage: type(f)
        <class 'sage.manifolds.coord_func_symb.CoordFunctionSymbRing_with_category.element_class'>
        sage: f.display()
        (x, y) |--> x^2 + 3*y + 1
        sage: f(x,y)
        x^2 + 3*y + 1

    The symbolic expression is returned when asking the direct display of
    the function::

        sage: f
        x^2 + 3*y + 1
        sage: latex(f)
        x^{2} + 3 \, y + 1

    A similar output is obtained by means of the method :meth:`expr`::

        sage: f.expr()
        x^2 + 3*y + 1

    The value of the function at specified coordinates is obtained by means
    of the standard parentheses notation::

        sage: f(2,-1)
        2
        sage: var('a b')
        (a, b)
        sage: f(a,b)
        a^2 + 3*b + 1

    An unspecified coordinate function::

        sage: g = X.function(function('G')(x, y))
        sage: g
        G(x, y)
        sage: g.display()
        (x, y) |--> G(x, y)
        sage: g.expr()
        G(x, y)
        sage: g(2,3)
        G(2, 3)

    Coordinate functions can be compared to other values::

        sage: f = X.function(x^2+3*y+1)
        sage: f == 2
        False
        sage: f == x^2 + 3*y + 1
        True
        sage: g = X.function(x*y)
        sage: f == g
        False
        sage: h = X.function(x^2+3*y+1)
        sage: f == h
        True

    .. RUBRIC:: Differences between ``CoordFunctionSymb`` and callable
      symbolic expressions

    Callable symbolic expressions are defined directly from symbolic
    expressions of the coordinates::

        sage: f0(x,y) = x^2 + 3*y + 1
        sage: type(f0)
        <type 'sage.symbolic.expression.Expression'>
        sage: f0
        (x, y) |--> x^2 + 3*y + 1
        sage: f0(x,y)
        x^2 + 3*y + 1

    To get an output similar to that of ``f0`` for the coordinate function
    ``f``, we must use the method :meth:`display`::

        sage: f
        x^2 + 3*y + 1
        sage: f.display()
        (x, y) |--> x^2 + 3*y + 1
        sage: f(x,y)
        x^2 + 3*y + 1

    More importantly, instances of :class:`CoordFunctionSymb` differ from
    callable symbolic expression by the automatic simplifications in all
    operations. For instance, adding the two callable symbolic expressions::

        sage: f0(x,y,z) = cos(x)^2 ; g0(x,y,z) = sin(x)^2

    results in::

        sage: f0 + g0
        (x, y, z) |--> cos(x)^2 + sin(x)^2

    To get `1`, one has to call
    :meth:`~sage.symbolic.expression.Expression.simplify_trig`::

        sage: (f0 + g0).simplify_trig()
        (x, y, z) |--> 1

    On the contrary, the sum of the corresponding :class:`CoordFunctionSymb`
    instances is automatically simplified (see
    :func:`~sage.manifolds.utilities.simplify_chain_real` and
    :func:`~sage.manifolds.utilities.simplify_chain_generic` for details)::

        sage: f = X.function(cos(x)^2) ; g = X.function(sin(x)^2)
        sage: f + g
        1

    Another difference regards the display of partial derivatives:
    for callable symbolic functions, it relies on Pynac notation
    ``D[0]``, ``D[1]``, etc.::

        sage: g = function('g')(x, y)
        sage: f0(x,y) = diff(g, x) + diff(g, y)
        sage: f0
        (x, y) |--> D[0](g)(x, y) + D[1](g)(x, y)

    while for coordinate functions, the display is more "textbook" like::

        sage: f = X.function(diff(g, x) + diff(g, y))
        sage: f
        d(g)/dx + d(g)/dy

    The difference is even more dramatic on LaTeX outputs::

        sage: latex(f0)
        \left( x, y \right) \ {\mapsto} \ D[0]\left(g\right)\left(x, y\right)
         + D[1]\left(g\right)\left(x, y\right)
        sage: latex(f)
        \frac{\partial\,g}{\partial x} + \frac{\partial\,g}{\partial y}

    Note that this regards only the display of coordinate functions:
    internally, the Pynac notation is still used, as we can check by asking
    for the symbolic expression stored in ``f``::

        sage: f.expr()
        D[0](g)(x, y) + D[1](g)(x, y)

    One can switch to Pynac notation by changing the options::

        sage: Manifold.options.textbook_output=False
        sage: latex(f)
        D[0]\left(g\right)\left(x, y\right) + D[1]\left(g\right)\left(x, y\right)
        sage: Manifold.options._reset()
        sage: latex(f)
        \frac{\partial\,g}{\partial x} + \frac{\partial\,g}{\partial y}

    Another difference between :class:`CoordFunctionSymb` and
    callable symbolic expression is the possibility to switch off the display
    of the arguments of unspecified functions. Consider for instance::

        sage: f = X.function(function('u')(x, y) * function('v')(x, y))
        sage: f
        u(x, y)*v(x, y)
        sage: f0(x,y) = function('u')(x, y) * function('v')(x, y)
        sage: f0
        (x, y) |--> u(x, y)*v(x, y)

    If there is a clear understanding that `u` and `v` are functions of
    `(x,y)`, the explicit mention of the latter can be cumbersome in lengthy
    tensor expressions. We can switch it off by::

        sage: Manifold.options.omit_function_arguments=True
        sage: f
        u*v

    Note that neither the callable symbolic expression ``f0`` nor the internal
    expression of ``f`` is affected by the above command::

        sage: f0
        (x, y) |--> u(x, y)*v(x, y)
        sage: f.expr()
        u(x, y)*v(x, y)

    We revert to the default behavior by::

        sage: Manifold.options._reset()
        sage: f
        u(x, y)*v(x, y)

    """
    def __init__(self, parent, expression):
        r"""
        Initialize ``self``.

        TESTS:

        Coordinate function on a real manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(1+x*y); f
            x*y + 1
            sage: type(f)
            <class 'sage.manifolds.coord_func_symb.CoordFunctionSymbRing_with_category.element_class'>
            sage: TestSuite(f).run()

        Coordinate function on a complex manifold::

            sage: N = Manifold(2, 'N', structure='topological', field='complex')
            sage: Y.<z,w> = N.chart()
            sage: g = Y.function(i*z + 2*w); g
            2*w + I*z
            sage: TestSuite(g).run()

        """
        CoordFunction.__init__(self, parent)
        self._express = SR(expression)  # symbolic expression enforced
        # Definition of the simplification chain to be applied in
        # symbolic calculus:
        if self.parent()._chart.manifold().base_field_type() == 'real':
            self._simplify = simplify_chain_real
        else:
            self._simplify = simplify_chain_generic
        # Derived quantities:
        self._der = None  # list of partial derivatives (to be set by diff()
                          # and unset by del_derived())

    # -------------------------------------------------------------
    # Methods to be implemented by derived classes of CoordFunction
    # -------------------------------------------------------------

    def _repr_(self):
        r"""
        String representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(1+x*y)
            sage: f._repr_()
            'x*y + 1'
            sage: repr(f)  # indirect doctest
            'x*y + 1'
            sage: f  # indirect doctest
            x*y + 1

        """
        if self.parent()._chart.manifold().options.textbook_output:
            return str(ExpressionNice(self._express))
        else:
            return str(self._express)

    def _latex_(self):
        r"""
        LaTeX representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(cos(x*y/2))
            sage: f._latex_()
            \cos\left(\frac{1}{2} \, x y\right)
            sage: latex(f)  # indirect doctest
            \cos\left(\frac{1}{2} \, x y\right)

        """
        from sage.misc.latex import latex
        if self.parent()._chart.manifold().options.textbook_output:
            return latex(ExpressionNice(self._express))
        else:
            return latex(self._express)

    def display(self):
        r"""
        Display ``self`` in arrow notation.

        The output is either text-formatted (console mode) or
        LaTeX-formatted (notebook mode).

        EXAMPLES:

        Coordinate function on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(cos(x*y/2))
            sage: f.display()
            (x, y) |--> cos(1/2*x*y)
            sage: latex(f.display())
            \left(x, y\right) \mapsto \cos\left(\frac{1}{2} \, x y\right)

        A shortcut is ``disp()``::

            sage: f.disp()
            (x, y) |--> cos(1/2*x*y)

        Display of the zero function::

            sage: X.zero_function().display()
            (x, y) |--> 0

        """
        from sage.tensor.modules.format_utilities import FormattedExpansion
        from sage.misc.latex import latex
        resu_txt = str(self.parent()._chart[:]) + ' |--> ' + \
                   str(ExpressionNice(self._express))
        resu_latex = latex(self.parent()._chart[:]) + r' \mapsto' + \
                     latex(ExpressionNice(self._express))
        return FormattedExpansion(resu_txt, resu_latex)

    disp = display

    def expr(self):
        r"""
        Return the symbolic expression representing the image of ``self``.

        OUTPUT:

        - :class:`symbolic expression <sage.symbolic.expression.Expression>`
          involving the chart coordinates

        EXAMPLES:

        Coordinate function of a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x^2+3*y+1)
            sage: f.expr()
            x^2 + 3*y + 1
            sage: type(f.expr())
            <type 'sage.symbolic.expression.Expression'>

        For a symbolic coordinate function, one shall always have::

            sage: bool( f.expr() == f(*(f.chart()[:])) )
            True

        The method :meth:`expr` is useful for accessing to all the
        symbolic expression functionalities in Sage; for instance::

            sage: var('a')
            a
            sage: f = X.function(a*x*y); f.display()
            (x, y) |--> a*x*y
            sage: f.expr()
            a*x*y
            sage: f.expr().subs(a=2)
            2*x*y

        Note that for substituting the value of a coordinate, the function
        call can be used as well::

            sage: f(x,3)
            3*a*x
            sage: bool( f(x,3) == f.expr().subs(y=3) )
            True

        """
        return self._express

    def __call__(self, *coords, **options):
        r"""
        Compute the value of the function at specified coordinates.

        INPUT:

        - ``*coords`` -- list of coordinates `(x^1, \ldots, x^n)`,
          where the function `f` is to be evaluated
        - ``**options`` -- allows to pass ``simplify=False`` to disable the
          call of the simplification chain on the result

        OUTPUT:

        - the value `f(x^1, \ldots, x^n)`, where `f` is the current
          coordinate function

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(sin(x*y))
            sage: f.__call__(-2, 3)
            -sin(6)
            sage: f(-2, 3)
            -sin(6)
            sage: var('a b')
            (a, b)
            sage: f.__call__(a, b)
            sin(a*b)
            sage: f(a,b)
            sin(a*b)
            sage: f.__call__(pi, 1)
            0
            sage: f.__call__(pi, 1/2)
            1

        """
        if len(coords) != self._nc:
            raise ValueError("bad number of coordinates")
        substitutions = dict(zip(self.parent()._chart._xx, coords))
        resu = self._express.subs(substitutions)
        if 'simplify' in options:
            if options['simplify']:
                return self._simplify(resu)
            else:
                return resu
        else:
            return self._simplify(resu)

    def is_zero(self):
        r"""
        Return ``True`` if the function is zero and ``False`` otherwise.

        EXAMPLES:

        Coordinate functions associated to a 2-dimensional chart::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x^2+3*y+1)
            sage: f.is_zero()
            False
            sage: f == 0
            False
            sage: g = X.function(0)
            sage: g.is_zero()
            True
            sage: g == 0
            True
            sage: X.zero_function().is_zero()
            True
            sage: X.zero_function() == 0
            True

        """
        return self._express.is_zero()

    def copy(self):
        r"""
        Return an exact copy of the object.

        OUTPUT:

        - a :class:`CoordFunctionSymb`

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x+y^2)
            sage: g = f.copy(); g
            y^2 + x

        By construction, ``g`` is identical to ``f``::

            sage: type(g) == type(f)
            True
            sage: g == f
            True

        but it is not the same object::

            sage: g is f
            False

        """
        return type(self)(self.parent(), self._express)

    def diff(self, coord):
        r"""
        Partial derivative with respect to a coordinate.

        INPUT:

        - ``coord`` -- either the coordinate `x^i` with respect
          to which the derivative of the coordinate function `f` is to be
          taken, or the index `i` labelling this coordinate (with the
          index convention defined on the chart domain via the parameter
          ``start_index``)

        OUTPUT:

        - a :class:`CoordFunctionSymb` representing the partial
          derivative `\frac{\partial f}{\partial x^i}`

        EXAMPLES:

        Partial derivatives of a 2-dimensional coordinate function::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x^2+3*y+1); f
            x^2 + 3*y + 1
            sage: f.diff(x)
            2*x
            sage: f.diff(y)
            3

        Each partial derivatives is itself a coordinate function::

            sage: type(f.diff(x))
            <class 'sage.manifolds.coord_func_symb.CoordFunctionSymbRing_with_category.element_class'>

        An index can be used instead of the coordinate symbol::

            sage: f.diff(0)
            2*x
            sage: f.diff(1)
            3

        The index range depends on the convention used on the chart's domain::

            sage: M = Manifold(2, 'M', structure='topological', start_index=1)
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x^2+3*y+1)
            sage: f.diff(0)
            Traceback (most recent call last):
            ...
            ValueError: coordinate index out of range
            sage: f.diff(1)
            2*x
            sage: f.diff(2)
            3

        """
        from sage.calculus.functional import diff
        from sage.rings.integer import Integer
        if self._der is None:
            # the list of partial derivatives has to be updated
            self._der = [type(self)(self.parent(),
                                    self._simplify(diff(self._express, xx)))
                         for xx in self.parent()._chart[:]]
        if isinstance(coord, (int, Integer)):
            # NB: for efficiency, we access directly to the "private" attributes
            # of other classes. A more conventional OOP writing would be
            # coordsi = coord - self.parent()._chart.domain().start_index()
            coordsi = coord - self.parent()._chart._domain._sindex
            if coordsi < 0 or coordsi >= self._nc:
                raise ValueError("coordinate index out of range")
            return self._der[coordsi]
        else:
            return self._der[self.parent()._chart[:].index(coord)]

    def __eq__(self, other):
        r"""
        Comparison (equality) operator.

        INPUT:

        - ``other`` -- a :class:`CoordFunction` or a value

        OUTPUT:

        - ``True`` if ``self`` is equal to ``other``,  or ``False`` otherwise

        TESTS:

        Coordinate functions associated to a 2-dimensional chart::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x+y^2)
            sage: g = X.function(x+y^2)
            sage: f == g
            True
            sage: f == 1
            False
            sage: h = X.function(1)
            sage: h == 1
            True
            sage: h == f
            False
            sage: h == 0
            False
            sage: X.function(0) == 0
            True
            sage: X.zero_function() == 0
            True

        """
        if other is self:
            return True
        if isinstance(other, CoordFunctionSymb):
            if other.parent() != self.parent():
                return False
            else:
                return bool(other._express == self._express)
        else:
            return bool(self._express == other)

    def __neg__(self):
        r"""
        Unary minus operator.

        OUTPUT:

        - the opposite of ``self``

        TESTS:

        Coordinate functions associated to a 2-dimensional chart::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x+y^2)
            sage: g = -f; g
            -y^2 - x
            sage: type(g)
            <class 'sage.manifolds.coord_func_symb.CoordFunctionSymbRing_with_category.element_class'>
            sage: -g == f
            True

        """
        return type(self)(self.parent(), self._simplify(-self._express))

    def __invert__(self):
        r"""
        Inverse operator.

        If `f` denotes the current coordinate function and `K` the topological
        field over which the manifold is defined, the *inverse* of `f` is the
        coordinate function `1/f`, where `1` of the multiplicative identity
        of `K`.

        OUTPUT:

        - the inverse of ``self``

        TESTS:

        Coordinate functions associated to a 2-dimensional chart::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(1+x^2+y^2)
            sage: g = f.__invert__(); g
            1/(x^2 + y^2 + 1)
            sage: type(g)
            <class 'sage.manifolds.coord_func_symb.CoordFunctionSymbRing_with_category.element_class'>
            sage: g == ~f
            True
            sage: g.__invert__() == f
            True

        """
        return type(self)(self.parent(),
                          self._simplify(SR.one() / self._express))
        # NB: self._express.__invert__() would return 1/self._express
        # (cf. the code of __invert__ in src/sage/symbolic/expression.pyx)
        # Here we prefer SR(1)/self._express

    def _add_(self, other):
        r"""
        Addition operator.

        INPUT:

        - ``other`` -- a :class:`CoordFunction` or a value

        OUTPUT:

        - coordinate function resulting from the addition of ``self``
          and ``other``

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x+y^2)
            sage: g = X.function(x+1)
            sage: s = f + g; s.display()
            (x, y) |--> y^2 + 2*x + 1
            sage: type(s)
            <class 'sage.manifolds.coord_func_symb.CoordFunctionSymbRing_with_category.element_class'>
            sage: (f + 0).display()
            (x, y) |--> y^2 + x
            sage: (f + X.zero_function()).display()
            (x, y) |--> y^2 + x
            sage: (f + 1).display()
            (x, y) |--> y^2 + x + 1
            sage: (f + pi).display()
            (x, y) |--> pi + y^2 + x
            sage: (f + x).display()
            (x, y) |--> y^2 + 2*x
            sage: (f + -f).display()
            (x, y) |--> 0

        """
        res = self._simplify(self._express + other._express)
        if res == 0:
            return self.parent().zero()
        else:
            return type(self)(self.parent(), res)

    def _sub_(self, other):
        r"""
        Subtraction operator.

        INPUT:

        - ``other`` -- a :class:`CoordFunction` or a value

        OUTPUT:

        - coordinate function resulting from the subtraction of ``other``
          from ``self``

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x+y^2)
            sage: g = X.function(x+1)
            sage: s = f - g; s.display()
            (x, y) |--> y^2 - 1
            sage: type(s)
            <class 'sage.manifolds.coord_func_symb.CoordFunctionSymbRing_with_category.element_class'>
            sage: (f - 0).display()
            (x, y) |--> y^2 + x
            sage: (f - X.zero_function()).display()
            (x, y) |--> y^2 + x
            sage: (f - 1).display()
            (x, y) |--> y^2 + x - 1
            sage: (f - x).display()
            (x, y) |--> y^2
            sage: (f - pi).display()
            (x, y) |--> -pi + y^2 + x
            sage: (f - f).display()
            (x, y) |--> 0
            sage: (f - g) == -(g - f)
            True
        """
        res = self._simplify(self._express - other._express)
        if res == 0:
            return self.parent().zero()
        else:
            return type(self)(self.parent(), res)

    def _mul_(self, other):
        r"""
        Multiplication operator.

        INPUT:

        - ``other`` -- a :class:`CoordFunction` or a value

        OUTPUT:

        - coordinate function resulting from the multiplication of ``self``
          by ``other``

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x+y)
            sage: g = X.function(x-y)
            sage: s = f._mul_(g); s.display()
            (x, y) |--> x^2 - y^2
            sage: type(s)
            <class 'sage.manifolds.coord_func_symb.CoordFunctionSymbRing_with_category.element_class'>
            sage: (f * 0).display()
            (x, y) |--> 0
            sage: (f * X.zero_function()).display()
            (x, y) |--> 0
            sage: (f * (1/f)).display()
            (x, y) |--> 1

        """
        res = self._simplify(self._express * other._express)
        if res == 0:
            return self.parent().zero()
        else:
            return type(self)(self.parent(), res)

    def _rmul_(self, other):
        """
        Return ``other * self``.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: one = X.function_ring().one()
            sage: 2 * one
            2
            sage: f = X.function(x+y)
            sage: (f * pi).display()
            (x, y) |--> pi*(x + y)
            sage: (x * f).display()
            (x, y) |--> (x + y)*x
        """
        try:
            other = SR(other)
        except (TypeError, ValueError):
            return
        return type(self)(self.parent(), other * self._express)

    def _lmul_(self, other):
        """
        Return ``self * other``.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: one = X.function_ring().one()
            sage: one * 2
            2
            sage: f = X.function(x+y)
            sage: (f * 2).display()
            (x, y) |--> 2*x + 2*y
            sage: (f * pi).display()
            (x, y) |--> pi*(x + y)
        """
        try:
            other = SR(other)
        except (TypeError, ValueError):
            return
        return type(self)(self.parent(), self._express * other)

    def _div_(self, other):
        r"""
        Division operator.

        INPUT:

        - ``other`` -- a :class:`CoordFunction` or a value

        OUTPUT:

        - coordinate function resulting from the division of ``self``
          by ``other``

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x+y)
            sage: g = X.function(1+x^2+y^2)
            sage: s = f._div_(g); s.display()
            (x, y) |--> (x + y)/(x^2 + y^2 + 1)
            sage: type(s)
            <class 'sage.manifolds.coord_func_symb.CoordFunctionSymbRing_with_category.element_class'>
            sage: f / X.zero_function()
            Traceback (most recent call last):
            ...
            ZeroDivisionError: division of a coordinate function by zero
            sage: (f / 1).display()
            (x, y) |--> x + y
            sage: (f / 2).display()
            (x, y) |--> 1/2*x + 1/2*y
            sage: (f / pi).display()
            (x, y) |--> (x + y)/pi
            sage: (f / (1+x^2)).display()
            (x, y) |--> (x + y)/(x^2 + 1)
            sage: (f / (1+x^2)).display()
            (x, y) |--> (x + y)/(x^2 + 1)
            sage: (f / g) == ~(g / f)
            True

        """
        if other._express.is_zero():
            raise ZeroDivisionError("division of a coordinate function by zero")
        res = self._simplify(self._express / SR(other))
        if res == 0:
            return self.parent().zero()
        else:
            return type(self)(self.parent(), res)

    def exp(self):
        r"""
        Exponential of ``self``.

        OUTPUT:

        - coordinate function `\exp(f)`, where `f` is the current
          coordinate function

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x+y)
            sage: f.exp()
            e^(x + y)
            sage: exp(f) # equivalent to f.exp()
            e^(x + y)
            sage: exp(f).display()
            (x, y) |--> e^(x + y)
            sage: exp(X.zero_function())
            1

        """
        return type(self)(self.parent(), self._simplify(self._express.exp()))

    def log(self, base=None):
        r"""
        Logarithm of ``self``.

        INPUT:

        - ``base`` -- (default: ``None``) base of the logarithm; if ``None``,
          the natural logarithm (i.e. logarithm to base `e`) is returned

        OUTPUT:

        - coordinate function `\log_a(f)`, where `f` is the current coordinate
          function and `a` is the base

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x+y)
            sage: f.log()
            log(x + y)
            sage: log(f) # equivalent to f.log()
            log(x + y)
            sage: log(f).display()
            (x, y) |--> log(x + y)
            sage: f.log(2)
            log(x + y)/log(2)
            sage: log(f, 2)
            log(x + y)/log(2)

        """
        return type(self)(self.parent(), self._simplify(self._express.log(base)))

    def __pow__(self, exponent):
        r"""
        Power of ``self``.

        INPUT:

        - ``exponent`` -- the exponent

        OUTPUT:

        - coordinate function `f^a`, where `f` is the current coordinate
          function and `a` is the exponent

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x+y)
            sage: f.__pow__(3)
            x^3 + 3*x^2*y + 3*x*y^2 + y^3
            sage: f^3  # equivalent to f.__pow__(3)
            x^3 + 3*x^2*y + 3*x*y^2 + y^3
            sage: f.__pow__(3).display()
            (x, y) |--> x^3 + 3*x^2*y + 3*x*y^2 + y^3
            sage: pow(f,3).display()
            (x, y) |--> x^3 + 3*x^2*y + 3*x*y^2 + y^3
            sage: (f^3).display()
            (x, y) |--> x^3 + 3*x^2*y + 3*x*y^2 + y^3
            sage: pow(X.zero_function(), 3).display()
            (x, y) |--> 0

        """
        return type(self)(self.parent(),
                          self._simplify(pow(self._express, exponent)))

    def sqrt(self):
        r"""
        Square root of ``self``.

        OUTPUT:

        - coordinate function `\sqrt{f}`, where `f` is the current
          coordinate function

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x+y)
            sage: f.sqrt()
            sqrt(x + y)
            sage: sqrt(f)  # equivalent to f.sqrt()
            sqrt(x + y)
            sage: sqrt(f).display()
            (x, y) |--> sqrt(x + y)
            sage: sqrt(X.zero_function()).display()
            (x, y) |--> 0

        """
        return type(self)(self.parent(),
                          self._simplify(self._express.sqrt()))

    def cos(self):
        r"""
        Cosine of ``self``.

        OUTPUT:

        - coordinate function `\cos(f)`, where `f` is the current
          coordinate function

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x*y)
            sage: f.cos()
            cos(x*y)
            sage: cos(f)  # equivalent to f.cos()
            cos(x*y)
            sage: cos(f).display()
            (x, y) |--> cos(x*y)
            sage: cos(X.zero_function()).display()
            (x, y) |--> 1

        """
        return type(self)(self.parent(),
                          self._simplify(self._express.cos()))

    def sin(self):
        r"""
        Sine of ``self``.

        OUTPUT:

        - coordinate function `\sin(f)`, where `f` is the current
          coordinate function

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x*y)
            sage: f.sin()
            sin(x*y)
            sage: sin(f)  # equivalent to f.sin()
            sin(x*y)
            sage: sin(f).display()
            (x, y) |--> sin(x*y)
            sage: sin(X.zero_function()) == X.zero_function()
            True

        """
        return type(self)(self.parent(),
                          self._simplify(self._express.sin()))

    def tan(self):
        r"""
        Tangent of ``self``.

        OUTPUT:

        - coordinate function `\tan(f)`, where `f` is the current
          coordinate function

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x*y)
            sage: f.tan()
            sin(x*y)/cos(x*y)
            sage: tan(f)  # equivalent to f.tan()
            sin(x*y)/cos(x*y)
            sage: tan(f).display()
            (x, y) |--> sin(x*y)/cos(x*y)
            sage: tan(X.zero_function()) == X.zero_function()
            True

        """
        return type(self)(self.parent(),
                          self._simplify(self._express.tan()))

    def arccos(self):
        r"""
        Arc cosine of ``self``.

        OUTPUT:

        - coordinate function `\arccos(f)`, where `f` is the current
          coordinate function

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x*y)
            sage: f.arccos()
            arccos(x*y)
            sage: arccos(f)  # equivalent to f.arccos()
            arccos(x*y)
            sage: acos(f)  # equivalent to f.arccos()
            arccos(x*y)
            sage: arccos(f).display()
            (x, y) |--> arccos(x*y)
            sage: arccos(X.zero_function()).display()
            (x, y) |--> 1/2*pi

        """
        return type(self)(self.parent(),
                          self._simplify(self._express.arccos()))

    def arcsin(self):
        r"""
        Arc sine of ``self``.

        OUTPUT:

        - coordinate function `\arcsin(f)`, where `f` is the current
          coordinate function

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x*y)
            sage: f.arcsin()
            arcsin(x*y)
            sage: arcsin(f)  # equivalent to f.arcsin()
            arcsin(x*y)
            sage: asin(f)  # equivalent to f.arcsin()
            arcsin(x*y)
            sage: arcsin(f).display()
            (x, y) |--> arcsin(x*y)
            sage: arcsin(X.zero_function()) == X.zero_function()
            True

        """
        return type(self)(self.parent(),
                          self._simplify(self._express.arcsin()))

    def arctan(self):
        r"""
        Arc tangent of ``self``.

        OUTPUT:

        - coordinate function `\arctan(f)`, where `f` is the current
          coordinate function

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x*y)
            sage: f.arctan()
            arctan(x*y)
            sage: arctan(f)  # equivalent to f.arctan()
            arctan(x*y)
            sage: atan(f)  # equivalent to f.arctan()
            arctan(x*y)
            sage: arctan(f).display()
            (x, y) |--> arctan(x*y)
            sage: arctan(X.zero_function()) == X.zero_function()
            True

        """
        return type(self)(self.parent(),
                          self._simplify(self._express.arctan()))

    def cosh(self):
        r"""
        Hyperbolic cosine of ``self``.

        OUTPUT:

        - coordinate function `\cosh(f)`, where `f` is the current
          coordinate function

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x*y)
            sage: f.cosh()
            cosh(x*y)
            sage: cosh(f)  # equivalent to f.cosh()
            cosh(x*y)
            sage: cosh(f).display()
            (x, y) |--> cosh(x*y)
            sage: cosh(X.zero_function()).display()
            (x, y) |--> 1

        """
        return type(self)(self.parent(),
                          self._simplify(self._express.cosh()))

    def sinh(self):
        r"""
        Hyperbolic sine of ``self``.

        OUTPUT:

        - coordinate function `\sinh(f)`, where `f` is the current
          coordinate function

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x*y)
            sage: f.sinh()
            sinh(x*y)
            sage: sinh(f)  # equivalent to f.sinh()
            sinh(x*y)
            sage: sinh(f).display()
            (x, y) |--> sinh(x*y)
            sage: sinh(X.zero_function()) == X.zero_function()
            True

        """
        return type(self)(self.parent(),
                          self._simplify(self._express.sinh()))

    def tanh(self):
        r"""
        Hyperbolic tangent of ``self``.

        OUTPUT:

        - coordinate function `\tanh(f)`, where `f` is the current
          coordinate function

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x*y)
            sage: f.tanh()
            sinh(x*y)/cosh(x*y)
            sage: tanh(f)  # equivalent to f.tanh()
            sinh(x*y)/cosh(x*y)
            sage: tanh(f).display()
            (x, y) |--> sinh(x*y)/cosh(x*y)
            sage: tanh(X.zero_function()) == X.zero_function()
            True

        """
        return type(self)(self.parent(),
                          self._simplify(self._express.tanh()))

    def arccosh(self):
        r"""
        Inverse hyperbolic cosine of ``self``.

        OUTPUT:

        - coordinate function `\mathrm{arcosh}(f)`, where `f` is the current
          coordinate function

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x*y)
            sage: f.arccosh()
            arccosh(x*y)
            sage: arccosh(f)  # equivalent to f.arccosh()
            arccosh(x*y)
            sage: acosh(f)  # equivalent to f.arccosh()
            arccosh(x*y)
            sage: arccosh(f).display()
            (x, y) |--> arccosh(x*y)
            sage: arccosh(X.function(1)) == X.zero_function()
            True

        """
        return type(self)(self.parent(),
                          self._simplify(self._express.arccosh()))

    def arcsinh(self):
        r"""
        Inverse hyperbolic sine of ``self``.

        OUTPUT:

        - coordinate function `\mathrm{arsinh}(f)`, where `f` is the current
          coordinate function

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x*y)
            sage: f.arcsinh()
            arcsinh(x*y)
            sage: arcsinh(f)  # equivalent to f.arcsinh()
            arcsinh(x*y)
            sage: asinh(f)  # equivalent to f.arcsinh()
            arcsinh(x*y)
            sage: arcsinh(f).display()
            (x, y) |--> arcsinh(x*y)
            sage: arcsinh(X.zero_function()) == X.zero_function()
            True

        """
        return type(self)(self.parent(),
                          self._simplify(self._express.arcsinh()))

    def arctanh(self):
        r"""
        Inverse hyperbolic tangent of ``self``.

        OUTPUT:

        - coordinate function `\mathrm{artanh}(f)`, where `f` is the
          current coordinate function

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x*y)
            sage: f.arctanh()
            arctanh(x*y)
            sage: arctanh(f)  # equivalent to f.arctanh()
            arctanh(x*y)
            sage: atanh(f)  # equivalent to f.arctanh()
            arctanh(x*y)
            sage: arctanh(f).display()
            (x, y) |--> arctanh(x*y)
            sage: arctanh(X.zero_function()) == X.zero_function()
            True

        """
        return type(self)(self.parent(),
                          self._simplify(self._express.arctanh()))

    # -------------------------------------
    # Methods specific to CoordFunctionSymb
    # -------------------------------------

    def _del_derived(self):
        r"""
        Delete the derived quantities.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(cos(x*y))
            sage: f._der
            sage: f.diff(x)
            -y*sin(x*y)
            sage: f._der
            [-y*sin(x*y), -x*sin(x*y)]
            sage: f._del_derived()
            sage: f._der

        """
        self._der = None  # reset of the partial derivatives

    def simplify(self):
        r"""
        Simplify the coordinate expression of ``self``.

        For details about the employed chain of simplifications, see
        :func:`~sage.manifolds.utilities.simplify_chain_real` for coordinate
        functions on real manifolds and
        :func:`~sage.manifolds.utilities.simplify_chain_generic` for the
        generic case.

        OUTPUT:

        - ``self`` with its expression simplified

        EXAMPLES:

        Simplification of a 2-dimension coordinate function::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(cos(x)^2+sin(x)^2 + sqrt(x^2))
            sage: f.display()
            (x, y) |--> cos(x)^2 + sin(x)^2 + sqrt(x^2)
            sage: f.simplify()
            abs(x) + 1

        The method ``simplify()`` has changed the expression of ``f``::

            sage: f.display()
            (x, y) |--> abs(x) + 1

        Another example::

            sage: f = X.function((x^2-1)/(x+1))
            sage: f
            (x^2 - 1)/(x + 1)
            sage: f.simplify()
            x - 1

        Examples taking into account the declared range of a coordinate::

            sage: M =  Manifold(2, 'M_1', structure='topological')
            sage: X.<x,y> = M.chart('x:(0,+oo) y')
            sage: f = X.function(sqrt(x^2))
            sage: f
            x
            sage: f.simplify()
            x

        ::

            sage: forget()  # to clear the previous assumption on x
            sage: M =  Manifold(2, 'M_2', structure='topological')
            sage: X.<x,y> = M.chart('x:(-oo,0) y')
            sage: f = X.function(sqrt(x^2))
            sage: f
            sqrt(x^2)
            sage: f.simplify()
            -x

        """
        self._express = self._simplify(self._express)
        self._del_derived()
        return self

    def factor(self):
        r"""
        Factorize the coordinate expression of ``self``.

        OUTPUT:

        - ``self`` with its expression factorized

        EXAMPLES:

        Factorization of a 2-dimensional coordinate function::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x^2 + 2*x*y + y^2)
            sage: f.display()
            (x, y) |--> x^2 + 2*x*y + y^2
            sage: f.factor()
            (x + y)^2

        The method ``factor()`` has changed the expression of ``f``::

            sage: f.display()
            (x, y) |--> (x + y)^2

        """
        self._express = self._express.factor()
        self._del_derived()
        return self

    def expand(self):
        r"""
        Expand the coordinate expression of ``self``.

        OUTPUT:

        - ``self`` with its expression expanded

        EXAMPLES:

        Expanding a 2-dimensional coordinate function::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function((x - y)^2)
            sage: f.display()
            (x, y) |--> (x - y)^2
            sage: f.expand()
            x^2 - 2*x*y + y^2

        The method ``expand()`` has changed the expression of ``f``::

            sage: f.display()
            (x, y) |--> x^2 - 2*x*y + y^2

        """
        self._express = self._express.expand()
        self._del_derived()
        return self

    def collect(self, s):
        r"""
        Collect the coefficients of `s` in the expression of ``self``
        into a group.

        INPUT:

        - ``s`` -- the symbol whose coefficients will be collected

        OUTPUT:

        - ``self`` with the coefficients of ``s`` grouped in
          its expression

        EXAMPLES:

        Action on a 2-dimensional coordinate function::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x^2*y + x*y + (x*y)^2)
            sage: f.display()
            (x, y) |--> x^2*y^2 + x^2*y + x*y
            sage: f.collect(y)
            x^2*y^2 + (x^2 + x)*y

        The method ``collect()`` has changed the expression of ``f``::

            sage: f.display()
            (x, y) |--> x^2*y^2 + (x^2 + x)*y

        """
        self._express = self._express.collect(s)
        self._del_derived()
        return self

    def collect_common_factors(self):
        r"""
        Collect common factors in the expression of ``self``.

        This method does not perform a full factorization but only looks
        for factors which are already explicitly present.

        OUTPUT:

        - ``self`` with the common factors collected in
          its expression

        EXAMPLES:

        Action on a 2-dimensional coordinate function::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = X.function(x/(x^2*y + x*y))
            sage: f.display()
            (x, y) |--> x/(x^2*y + x*y)
            sage: f.collect_common_factors()
            1/((x + 1)*y)


        The method ``collect_common_factors()`` has changed the expression
        of ``f``::

            sage: f.display()
            (x, y) |--> 1/((x + 1)*y)

        """
        self._express = self._express.collect_common_factors()
        self._del_derived()
        return self


class CoordFunctionSymbRing(Parent, UniqueRepresentation):
    """
    Ring of all symbolic coordinate functions on a chart.

    INPUT:

    - ``chart`` -- a coordinate chart, as an instance of class
      :class:`~sage.manifolds.chart.Chart`

    EXAMPLES::

        sage: M = Manifold(2, 'M', structure='topological')
        sage: X.<x,y> = M.chart()
        sage: FR = X.function_ring(); FR
        Ring of coordinate functions on Chart (M, (x, y))
        sage: type(FR)
        <class 'sage.manifolds.coord_func_symb.CoordFunctionSymbRing_with_category'>
        sage: FR.category()
        Category of commutative algebras over Symbolic Ring

    """
    def __init__(self, chart):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: FR = X.function_ring()
            sage: TestSuite(FR).run()
        """
        self._chart = chart
        Parent.__init__(self, base=SR, category=CommutativeAlgebras(SR))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: X.function_ring()
            Ring of coordinate functions on Chart (M, (x, y))
        """
        return "Ring of coordinate functions on {}".format(self._chart)

    @cached_method
    def zero(self):
        """
        Return the constant function `0` in ``self``.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: FR = X.function_ring()
            sage: FR.zero()
            0

            sage: M = Manifold(2, 'M', structure='topological', field=Qp(5))
            sage: X.<x,y> = M.chart()
            sage: X.function_ring().zero()
            0
        """
        if self._chart.manifold().base_field_type() in ['real', 'complex']:
            elt = SR.zero()
        else:
            elt = self._chart.manifold().base_field().zero()
        return self.element_class(self, elt)

    @cached_method
    def one(self):
        """
        Return the constant function `1` in ``self``.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: FR = X.function_ring()
            sage: FR.one()
            1

            sage: M = Manifold(2, 'M', structure='topological', field=Qp(5))
            sage: X.<x,y> = M.chart()
            sage: X.function_ring().one()
            1 + O(5^20)
        """
        if self._chart.manifold().base_field_type() in ['real', 'complex']:
            elt = SR.one()
        else:
            elt = self._chart.manifold().base_field().one()
        return self.element_class(self, elt)

    def from_base_ring(self, r):
        """
        Return the canonical embedding of ``r`` into ``self``.

        INPUT:

        - ``r`` -- an element of ``self.base_ring()``

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: FR = X.function_ring()
            sage: f = FR.from_base_ring(x*y)
            sage: f.display()
            (x, y) |--> x*y

        """
        return self.element_class(self, r)

    def is_integral_domain(self):
        """
        Return ``False`` as ``self`` is not an integral domain.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: FR = X.function_ring()
            sage: FR.is_integral_domain()
            False
            sage: FR.is_field()
            False
        """
        return False

    is_field = is_integral_domain

    Element = CoordFunctionSymb

