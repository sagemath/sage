r"""
Symbolic Series

Symbolic series are special kinds of symbolic expressions that are
constructed via the
:meth:`Expression.series <sage.symbolic.expression.Expression.series>`
method.
They usually have an ``Order()`` term unless the series representation
is exact, see
:meth:`~sage.symbolic.expression.SymbolicSeries.is_terminating_series`.

For series over general rings see
:class:`power series <sage.rings.power_series_poly.PowerSeries_poly>`
and
:class:`Laurent series<sage.rings.laurent_series_ring_element.LaurentSeries>`.

EXAMPLES:

We expand a polynomial in `x` about `0`, about `1`, and also truncate
it back to a polynomial::

    sage: var('x,y')
    (x, y)
    sage: f = (x^3 - sin(y)*x^2 - 5*x + 3); f
    x^3 - x^2*sin(y) - 5*x + 3
    sage: g = f.series(x, 4); g
    3 + (-5)*x + (-sin(y))*x^2 + 1*x^3 + Order(x^4)
    sage: g.truncate()
    x^3 - x^2*sin(y) - 5*x + 3
    sage: g = f.series(x==1, oo); g
    (-sin(y) - 1) + (-2*sin(y) - 2)*(x - 1) + (-sin(y) + 3)*(x - 1)^2 + 1*(x - 1)^3
    sage: h = g.truncate(); h
    (x - 1)^3 - (x - 1)^2*(sin(y) - 3) - 2*(x - 1)*(sin(y) + 1) - sin(y) - 1
    sage: h.expand()
    x^3 - x^2*sin(y) - 5*x + 3

We compute another series expansion of an analytic function::

    sage: f = sin(x)/x^2
    sage: f.series(x,7)
    1*x^(-1) + (-1/6)*x + 1/120*x^3 + (-1/5040)*x^5 + Order(x^7)
    sage: f.series(x==1,3)
    (sin(1)) + (cos(1) - 2*sin(1))*(x - 1) + (-2*cos(1) + 5/2*sin(1))*(x - 1)^2 + Order((x - 1)^3)
    sage: f.series(x==1,3).truncate().expand()
    -2*x^2*cos(1) + 5/2*x^2*sin(1) + 5*x*cos(1) - 7*x*sin(1) - 3*cos(1) + 11/2*sin(1)

Following the GiNaC tutorial, we use John Machin's amazing
formula `\pi = 16 \mathrm{tan}^{-1}(1/5) - 4 \mathrm{tan}^{-1}(1/239)`
to compute digits of `\pi`. We expand the arc tangent around 0 and insert
the fractions 1/5 and 1/239.

::

    sage: x = var('x')
    sage: f = atan(x).series(x, 10); f
    1*x + (-1/3)*x^3 + 1/5*x^5 + (-1/7)*x^7 + 1/9*x^9 + Order(x^10)
    sage: (16*f.subs(x==1/5) - 4*f.subs(x==1/239)).n()
    3.14159268240440

Note: The result of an operation or function of series is not automatically
expanded to a series. This must be explicitly done by the user::

    sage: ex1 = sin(x).series(x, 4); ex1
    1*x + (-1/6)*x^3 + Order(x^4)
    sage: ex2 = cos(x).series(x, 4); ex2
    1 + (-1/2)*x^2 + Order(x^4)
    sage: ex1 + ex2
    (1 + (-1/2)*x^2 + Order(x^4)) + (1*x + (-1/6)*x^3 + Order(x^4))
    sage: (ex1 + ex2).series(x,4)
    1 + 1*x + (-1/2)*x^2 + (-1/6)*x^3 + Order(x^4)
    sage: x*ex1
    x*(1*x + (-1/6)*x^3 + Order(x^4))
    sage: (x*ex1).series(x,5)
    1*x^2 + (-1/6)*x^4 + Order(x^5)
    sage: sin(ex1)
    sin(1*x + (-1/6)*x^3 + Order(x^4))
    sage: sin(ex1).series(x,9)
    1*x + (-1/3)*x^3 + 11/120*x^5 + (-53/2520)*x^7 + Order(x^9)
    sage: (sin(x^2)^(-5)).series(x,3)
    1*x^(-10) + 5/6*x^(-6) + 3/8*x^(-2) + 367/3024*x^2 + Order(x^3)
    sage: (cot(x)^(-3)).series(x,3)
    Order(x^3)
    sage: (cot(x)^(-3)).series(x,4)
    1*x^3 + Order(x^4)

TESTS:

Check that :issue:`20088` is fixed::

    sage: ((1+x).series(x)^pi).series(x,3)
    1 + pi*x + (-1/2*pi + 1/2*pi^2)*x^2 + Order(x^3)

Check that :issue:`14878` is fixed, this should take only microseconds::

    sage: sin(x*sin(x*sin(x*sin(x)))).series(x,8)
    1*x^4 + (-1/6)*x^6 + Order(x^8)
    sage: sin(x*sin(x*sin(x*sin(x)))).series(x,12)
    1*x^4 + (-1/6)*x^6 + (-19/120)*x^8 + (-421/5040)*x^10 + Order(x^12)

Check that :issue:`22959` is fixed::

    sage: (x/(1-x^2)).series(x==0, 10)
    1*x + 1*x^3 + 1*x^5 + 1*x^7 + 1*x^9 + Order(x^10)
    sage: (x/(1-x^2)).series(x==0, 11)
    1*x + 1*x^3 + 1*x^5 + 1*x^7 + 1*x^9 + Order(x^11)
    sage: (x^2/(1-x^2)).series(x==0, 10)
    1*x^2 + 1*x^4 + 1*x^6 + 1*x^8 + Order(x^10)
    sage: (x^2/(1-x^2)).series(x==0, 11)
    1*x^2 + 1*x^4 + 1*x^6 + 1*x^8 + 1*x^10 + Order(x^11)

Check that :issue:`22733` is fixed::

    sage: _ = var('z')
    sage: z.series(x)
    (z)
"""

# ****************************************************************************
#       Copyright (C) 2015 Ralf Stephan <ralf@ark.in-berlin.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


cpdef bint _is_order(Expression expr):
    """
    Return whether ``expr`` has the form ``Order(...)``.

    TESTS::

        sage: from sage.symbolic.expression import _is_order
        sage: s = (x^2).series(x,1)
        sage: expr = Expression.coefficients(s)[0][0]; expr
        Order(x)
        sage: _is_order(expr)
        True
        sage: expr = Expression.coefficients(x.series(x, 0))[0][0]; expr
        Order(1)
        sage: _is_order(expr)
        True
        sage: expr = Expression.coefficients((x^5).series(x, 3))[0][0]; expr
        Order(x^3)
        sage: _is_order(expr)
        True
        sage: from sage.functions.other import Order; _is_order(Order(x))
        True
        sage: _is_order(Order(x) + 1)
        False
        sage: _is_order(Order(x) + 1 - 1)
        True
        sage: var("y")
        y
        sage: _is_order(Order(y))
        True
    """
    from sage.functions.other import Order
    from sage.symbolic.operators import add_vararg
    if expr.operator() is add_vararg and len(expr.operands()) == 1:
        expr = expr.operands()[0]
    return expr.operator() is Order


cdef class SymbolicSeries(Expression):
    def __init__(self, SR):
        """
        Trivial constructor.

        EXAMPLES::

            sage: loads(dumps((x+x^3).series(x,2)))
            1*x + Order(x^2)
        """
        Expression.__init__(self, SR, 0)
        self._parent = SR

    def is_terminating_series(self):
        """
        Return ``True`` if the series is without order term.

        A series is terminating if it can be represented exactly,
        without requiring an order term. You can explicitly
        request terminating series by setting the order to
        positive infinity.

        OUTPUT: boolean; ``True`` if the series has no order term

        EXAMPLES::

            sage: (x^5+x^2+1).series(x, +oo)
            1 + 1*x^2 + 1*x^5
            sage: (x^5+x^2+1).series(x,+oo).is_terminating_series()
            True
            sage: SR(5).is_terminating_series()
            False
            sage: exp(x).series(x,10).is_terminating_series()
            False
        """
        return g_is_a_terminating_series((<Expression>self)._gobj)

    def truncate(self):
        """
        Given a power series or expression, return the corresponding
        expression without the big oh.

        OUTPUT: a symbolic expression

        EXAMPLES::

            sage: f = sin(x)/x^2
            sage: f.truncate()
            sin(x)/x^2
            sage: f.series(x,7)
            1*x^(-1) + (-1/6)*x + 1/120*x^3 + (-1/5040)*x^5 + Order(x^7)
            sage: f.series(x,7).truncate()
            -1/5040*x^5 + 1/120*x^3 - 1/6*x + 1/x
            sage: f.series(x==1,3).truncate().expand()
            -2*x^2*cos(1) + 5/2*x^2*sin(1) + 5*x*cos(1) - 7*x*sin(1) - 3*cos(1) + 11/2*sin(1)
        """
        return new_Expression_from_GEx(self._parent, series_to_poly(self._gobj))

    def default_variable(self):
        """
        Return the expansion variable of this symbolic series.

        EXAMPLES::

            sage: s = (1/(1-x)).series(x,3); s
            1 + 1*x + 1*x^2 + Order(x^3)
            sage: s.default_variable()
            x
        """
        cdef GEx x = g_series_var(self._gobj)
        cdef Expression ex = new_Expression_from_GEx(self._parent, x)
        return ex

    def coefficients(self, x=None, sparse=True):
        r"""
        Return the coefficients of this symbolic series. See :meth:`Expression.coefficients` for more details.

        When ``sparse=False``, unlike :meth:`Expression.coefficients`, the dense list of coefficients
        are padded to the degree of the asymptotic term of this series.
        Also, this correctly handles the case where there is no term below the asymptotic term.

        EXAMPLES::

            sage: s = (1/(1-x)).series(x,6); s
            1 + 1*x + 1*x^2 + 1*x^3 + 1*x^4 + 1*x^5 + Order(x^6)
            sage: s.coefficients()
            [[1, 0], [1, 1], [1, 2], [1, 3], [1, 4], [1, 5]]
            sage: s.coefficients(x, sparse=False)
            [1, 1, 1, 1, 1, 1]
            sage: x,y = var("x,y")
            sage: s = (1/(1-y*x-x)).series(x,3); s
            1 + (y + 1)*x + ((y + 1)^2)*x^2 + Order(x^3)
            sage: s.coefficients(x, sparse=False)
            [1, y + 1, (y + 1)^2]

        TESTS::

            sage: s = (x^-1 + x^2).series(x,6)
            sage: s.coefficients()
            [[1, -1], [1, 2]]
            sage: s.coefficients(sparse=False)
            Traceback (most recent call last):
            ...
            ValueError: cannot return dense coefficient list with negative valuation
            sage: Expression.coefficients(s)
            [[1, -1], [1, 2]]
            sage: Expression.coefficients(s, sparse=False)
            Traceback (most recent call last):
            ...
            ValueError: cannot return dense coefficient list with negative valuation

            sage: s = (x+x^3).series(x,6)
            sage: s.coefficients()
            [[1, 1], [1, 3]]
            sage: s.coefficients(sparse=False)
            [0, 1, 0, 1, 0, 0]
            sage: Expression.coefficients(s)
            [[1, 1], [1, 3]]
            sage: Expression.coefficients(s, sparse=False)
            [0, 1, 0, 1]

            sage: s = (x^5).series(x,1)
            sage: s.coefficients()
            []
            sage: s.coefficients(sparse=False)
            [0]
            sage: Expression.coefficients(s)  # incorrect result here, but the user will not see this
            [[Order(x), 0]]
            sage: Expression.coefficients(s, sparse=False)
            [Order(x)]

            sage: s = (x^5).series(x,4)
            sage: s.coefficients()
            []
            sage: s.coefficients(sparse=False)
            [0, 0, 0, 0]
            sage: Expression.coefficients(s)  # incorrect result here, but the user will not see this
            [[Order(x^4), 0]]
            sage: Expression.coefficients(s, sparse=False)
            [Order(x^4)]
        """
        if x is None:
            x = self.default_variable()
        result = super().coefficients(x, sparse)
        if sparse and len(result) == 1 and ZZ(result[0][1]) == 0 and _is_order(result[0][0]):
            result = []
        if sparse:
            return result
        if len(result) == 1 and _is_order(result[0]):
            result = []
        return result + [ZZ.zero()] * (self.degree(x) - len(result))

    def power_series(self, base_ring):
        """
        Return the algebraic power series associated to this symbolic series.

        The coefficients must be coercible to the base ring.

        EXAMPLES::

            sage: ex = (gamma(1-x)).series(x,3); ex
            1 + euler_gamma*x + (1/2*euler_gamma^2 + 1/12*pi^2)*x^2 + Order(x^3)
            sage: g = ex.power_series(SR); g
            1 + euler_gamma*x + (1/2*euler_gamma^2 + 1/12*pi^2)*x^2 + O(x^3)
            sage: g.parent()
            Power Series Ring in x over Symbolic Ring
        """
        from sage.rings.power_series_ring import PowerSeriesRing
        R = PowerSeriesRing(base_ring, names=str(self.default_variable()))
        return R(self.list(), self.degree(self.default_variable()))
