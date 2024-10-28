r"""
Utilities for Calculus

This module defines helper functions which are used for simplifications
and display of symbolic expressions.

AUTHORS:

- Michal Bejger (2015) : class :class:`ExpressionNice`
- Eric Gourgoulhon (2015, 2017) : simplification functions
- Travis Scrimshaw (2016): review tweaks
- Marius Gerbershagen (2022) : skip simplification of expressions with a single
  number or symbolic variable
"""

# *****************************************************************************
#
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#       Copyright (C) 2015, 2017 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2016 Travis Scrimshaw <tscrimsh@umn.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from operator import pow as _pow
from sage.symbolic.expression import Expression
from sage.symbolic.expression_conversions import ExpressionTreeWalker
from sage.symbolic.ring import SR
from sage.symbolic.constants import pi
from sage.functions.other import abs_symbolic
from sage.misc.functional import sqrt
from sage.functions.trig import cos, sin
from sage.rings.rational import Rational


class SimplifySqrtReal(ExpressionTreeWalker):
    r"""
    Class for simplifying square roots in the real domain, by walking the
    expression tree.

    The end user interface is the function :func:`simplify_sqrt_real`.

    INPUT:

    - ``ex`` -- a symbolic expression

    EXAMPLES:

    Let us consider the square root of an exact square under some assumption::

        sage: assume(x<1)
        sage: a = sqrt(x^2-2*x+1)

    The method :meth:`~sage.symbolic.expression.Expression.simplify_full()`
    is ineffective on such an expression::

        sage: a.simplify_full()
        sqrt(x^2 - 2*x + 1)

    and the more aggressive method :meth:`~sage.symbolic.expression.Expression.canonicalize_radical()`
    yields a wrong result, given that `x<1`::

        sage: a.canonicalize_radical()  # wrong output!
        x - 1

    We construct a :class:`SimplifySqrtReal` object ``s`` from the symbolic
    expression ``a``::

        sage: from sage.manifolds.utilities import SimplifySqrtReal
        sage: s = SimplifySqrtReal(a)

    We use the ``__call__`` method to walk the expression tree and produce a
    correctly simplified expression::

        sage: s()
        -x + 1

    Calling the simplifier ``s`` with an expression actually simplifies this
    expression::

        sage: s(a)  # same as s() since s is built from a
        -x + 1
        sage: s(sqrt(x^2))
        abs(x)
        sage: s(sqrt(1+sqrt(x^2-2*x+1)))  # nested sqrt's
        sqrt(-x + 2)

    Another example where both
    :meth:`~sage.symbolic.expression.Expression.simplify_full()` and
    :meth:`~sage.symbolic.expression.Expression.canonicalize_radical()`
    fail::

        sage: b = sqrt((x-1)/(x-2))*sqrt(1-x)
        sage: b.simplify_full()  # does not simplify
        sqrt(-x + 1)*sqrt((x - 1)/(x - 2))
        sage: b.canonicalize_radical()  # wrong output, given that x<1
        (I*x - I)/sqrt(x - 2)
        sage: SimplifySqrtReal(b)()  # OK, given that x<1
        -(x - 1)/sqrt(-x + 2)

    TESTS:

    We check that the inverse of a square root is well simplified; this is a
    a non-trivial test since ``1/sqrt(x)`` is represented by ``pow(x,-1/2)``
    in the expression tree::

        sage: SimplifySqrtReal(1/sqrt(x^2-4*x+4))()
        -1/(x - 2)
        sage: SimplifySqrtReal(sqrt((x-2)/((x-3)*(x^2-2*x+1))))()
        -sqrt(-x + 2)/((x - 1)*sqrt(-x + 3))
        sage: forget()  # for doctests below

    .. SEEALSO::

        :func:`simplify_sqrt_real` for more examples with
        :class:`SimplifySqrtReal` at work.
    """
    def arithmetic(self, ex, operator):
        r"""
        This is the only method of the base class
        :class:`~sage.symbolic.expression_conversions.ExpressionTreeWalker`
        that is reimplemented, since square roots are considered as
        arithmetic operations with ``operator`` = ``pow`` and
        ``ex.operands()[1]`` = ``1/2`` or ``-1/2``.

        INPUT:

        - ``ex`` -- a symbolic expression
        - ``operator`` -- an arithmetic operator

        OUTPUT:

        - a symbolic expression, equivalent to ``ex`` with square roots
          simplified

        EXAMPLES::

            sage: from sage.manifolds.utilities import SimplifySqrtReal
            sage: a = sqrt(x^2+2*x+1)
            sage: s = SimplifySqrtReal(a)
            sage: a.operator()
            <built-in function pow>
            sage: s.arithmetic(a, a.operator())
            abs(x + 1)

        ::

            sage: a = x + 1  # no square root
            sage: s.arithmetic(a, a.operator())
            x + 1

        ::

            sage: a = x + 1 + sqrt(function('f')(x)^2)
            sage: s.arithmetic(a, a.operator())
            x + abs(f(x)) + 1
        """
        if operator is _pow:
            operands = ex.operands()
            power = operands[1]
            one_half = Rational((1,2))
            minus_one_half = -one_half
            if (power == one_half) or (power == minus_one_half):
                # This is a square root or the inverse of a square root
                w0 = SR.wild(0)
                w1 = SR.wild(1)
                sqrt_pattern = w0**one_half
                inv_sqrt_pattern = w0**minus_one_half
                sqrt_ratio_pattern1 = w0**one_half * w1**minus_one_half
                sqrt_ratio_pattern2 = w0**minus_one_half * w1**one_half
                argum = operands[0]  # the argument of sqrt
                if argum.has(sqrt_pattern) or argum.has(inv_sqrt_pattern):
                    argum = self(argum)  # treatment of nested sqrt's
                den = argum.denominator()
                if not (den == 1):  # the argument of sqrt is a fraction
                    # NB: after #19312 (integrated in Sage 6.10.beta7), the
                    # above test cannot be written as "if den != 1:"
                    num = argum.numerator()
                    if num < 0 or den < 0:
                        ex = sqrt(-num) / sqrt(-den)
                    else:
                        ex = sqrt(argum)
                else:
                    ex = sqrt(argum)
                simpl = SR(ex._maxima_().radcan())
                if (not simpl.match(sqrt_pattern) and
                    not simpl.match(inv_sqrt_pattern) and
                    not simpl.match(sqrt_ratio_pattern1) and
                    not simpl.match(sqrt_ratio_pattern2)):
                    # radcan transformed substantially the expression,
                    # possibly getting rid of some sqrt; in order to ensure a
                    # positive result, the absolute value of radcan's output
                    # is taken, the call to simplify() taking care of possible
                    # assumptions regarding signs of subexpression of simpl:
                    simpl = abs(simpl).simplify()
                if power == minus_one_half:
                    simpl = SR(1)/simpl
                return simpl
        # If operator is not a square root, we default to ExpressionTreeWalker:
        return super().arithmetic(ex, operator)


class SimplifyAbsTrig(ExpressionTreeWalker):
    r"""
    Class for simplifying absolute values of cosines or sines (in the real
    domain), by walking the expression tree.

    The end user interface is the function :func:`simplify_abs_trig`.

    INPUT:

    - ``ex`` -- a symbolic expression

    EXAMPLES:

    Let us consider the following symbolic expression with some assumption
    on the range of the variable `x`::

        sage: assume(pi/2<x, x<pi)
        sage: a = abs(cos(x)) + abs(sin(x))

    The method :meth:`~sage.symbolic.expression.Expression.simplify_full()`
    is ineffective on such an expression::

        sage: a.simplify_full()
        abs(cos(x)) + abs(sin(x))

    We construct a :class:`SimplifyAbsTrig` object ``s`` from the symbolic
    expression ``a``::

        sage: from sage.manifolds.utilities import SimplifyAbsTrig
        sage: s = SimplifyAbsTrig(a)

    We use the ``__call__`` method to walk the expression tree and produce a
    correctly simplified expression, given that `x\in(\pi/2, \pi)`::

        sage: s()
        -cos(x) + sin(x)

    Calling the simplifier ``s`` with an expression actually simplifies this
    expression::

        sage: s(a)  # same as s() since s is built from a
        -cos(x) + sin(x)
        sage: s(abs(cos(x/2)) + abs(sin(x/2)))  #  pi/4 < x/2 < pi/2
        cos(1/2*x) + sin(1/2*x)
        sage: s(abs(cos(2*x)) + abs(sin(2*x)))  #  pi < 2 x < 2*pi
        abs(cos(2*x)) - sin(2*x)
        sage: s(abs(sin(2+abs(cos(x)))))  # nested abs(sin_or_cos(...))
        sin(-cos(x) + 2)

    TESTS::

        sage: forget()  # for doctests below

    .. SEEALSO::

        :func:`simplify_abs_trig` for more examples with
        :class:`SimplifyAbsTrig` at work.
    """
    def composition(self, ex, operator):
        r"""
        This is the only method of the base class
        :class:`~sage.symbolic.expression_conversions.ExpressionTreeWalker`
        that is reimplemented, since it manages the composition of
        ``abs`` with ``cos`` or ``sin``.

        INPUT:

        - ``ex`` -- a symbolic expression
        - ``operator`` -- an operator

        OUTPUT:

        - a symbolic expression, equivalent to ``ex`` with ``abs(cos(...))``
          and ``abs(sin(...))`` simplified, according to the range of their
          argument.

        EXAMPLES::

            sage: from sage.manifolds.utilities import SimplifyAbsTrig
            sage: assume(-pi/2 < x, x<0)
            sage: a = abs(sin(x))
            sage: s = SimplifyAbsTrig(a)
            sage: a.operator()
            abs
            sage: s.composition(a, a.operator())
            sin(-x)

        ::

            sage: a = exp(function('f')(x))  # no abs(sin_or_cos(...))
            sage: a.operator()
            exp
            sage: s.composition(a, a.operator())
            e^f(x)

        ::

            sage: forget()  # no longer any assumption on x
            sage: a = abs(cos(sin(x)))  # simplifiable since -1 <= sin(x) <= 1
            sage: s.composition(a, a.operator())
            cos(sin(x))
            sage: a = abs(sin(cos(x)))  # not simplifiable
            sage: s.composition(a, a.operator())
            abs(sin(cos(x)))
        """
        if operator is abs_symbolic:
            argum = ex.operands()[0]  # argument of abs
            if argum.operator() is sin:
                # Case of abs(sin(...))
                x = argum.operands()[0]  # argument of sin
                w0 = SR.wild()
                if x.has(abs_symbolic(sin(w0))) or x.has(abs_symbolic(cos(w0))):
                    x = self(x)  # treatment of nested abs(sin_or_cos(...))
                # Simplifications for values of x in the range [-pi, 2*pi]:
                if x >= 0 and x <= pi:
                    ex = sin(x)
                elif (x > pi and x <= 2*pi) or (x >= -pi and x < 0):
                    ex = -sin(x)
                return ex
            if argum.operator() is cos:
                # Case of abs(cos(...))
                x = argum.operands()[0]  # argument of cos
                w0 = SR.wild()
                if x.has(abs_symbolic(sin(w0))) or x.has(abs_symbolic(cos(w0))):
                    x = self(x)  # treatment of nested abs(sin_or_cos(...))
                # Simplifications for values of x in the range [-pi, 2*pi]:
                if (x >= -pi/2 and x <= pi/2) or (x >= 3*pi/2 and x <= 2*pi):
                    ex = cos(x)
                elif (x > pi/2 and x <= 3*pi/2) or (x >= -pi and x < -pi/2):
                    ex = -cos(x)
                return ex
        # If no pattern is found, we default to ExpressionTreeWalker:
        return super().composition(ex, operator)


def simplify_sqrt_real(expr):
    r"""
    Simplify ``sqrt`` in symbolic expressions in the real domain.

    EXAMPLES:

    Simplifications of basic expressions::

        sage: from sage.manifolds.utilities import simplify_sqrt_real
        sage: simplify_sqrt_real( sqrt(x^2) )
        abs(x)
        sage: assume(x<0)
        sage: simplify_sqrt_real( sqrt(x^2) )
        -x
        sage: simplify_sqrt_real( sqrt(x^2-2*x+1) )
        -x + 1
        sage: simplify_sqrt_real( sqrt(x^2) + sqrt(x^2-2*x+1) )
        -2*x + 1

    This improves over
    :meth:`~sage.symbolic.expression.Expression.canonicalize_radical`,
    which yields incorrect results when ``x < 0``::

        sage: forget()  # removes the assumption x<0
        sage: sqrt(x^2).canonicalize_radical()
        x
        sage: assume(x<0)
        sage: sqrt(x^2).canonicalize_radical()
        -x
        sage: sqrt(x^2-2*x+1).canonicalize_radical() # wrong output
        x - 1
        sage: ( sqrt(x^2) + sqrt(x^2-2*x+1) ).canonicalize_radical() # wrong output
        -1

    Simplification of nested ``sqrt``'s::

        sage: forget()  # removes the assumption x<0
        sage: simplify_sqrt_real( sqrt(1 + sqrt(x^2)) )
        sqrt(abs(x) + 1)
        sage: assume(x<0)
        sage: simplify_sqrt_real( sqrt(1 + sqrt(x^2)) )
        sqrt(-x + 1)
        sage: simplify_sqrt_real( sqrt(x^2 + sqrt(4*x^2) + 1) )
        -x + 1

    Again, :meth:`~sage.symbolic.expression.Expression.canonicalize_radical`
    fails on the last one::

        sage: (sqrt(x^2 + sqrt(4*x^2) + 1)).canonicalize_radical()
        x - 1

    TESTS:

    Simplification of expressions involving some symbolic derivatives::

        sage: f = function('f')
        sage: simplify_sqrt_real( diff(f(x), x)/sqrt(x^2-2*x+1) )  # x<0 => x-1<0
        -diff(f(x), x)/(x - 1)
        sage: g = function('g')
        sage: simplify_sqrt_real( sqrt(x^3*diff(f(g(x)), x)^2) )  # x<0
        (-x)^(3/2)*abs(D[0](f)(g(x)))*abs(diff(g(x), x))
        sage: forget()  # for doctests below
    """
    w0 = SR.wild()
    one_half = Rational((1,2))
    if expr.has(w0**one_half) or expr.has(w0**(-one_half)):
        return SimplifySqrtReal(expr)()
    return expr


def simplify_abs_trig(expr):
    r"""
    Simplify ``abs(sin(...))`` and ``abs(cos(...))`` in symbolic expressions.

    EXAMPLES::

        sage: M = Manifold(3, 'M', structure='topological')
        sage: X.<x,y,z> = M.chart(r'x y:(0,pi) z:(-pi/3,0)')
        sage: X.coord_range()
        x: (-oo, +oo); y: (0, pi); z: (-1/3*pi, 0)

    Since `x` spans all `\RR`, no simplification of ``abs(sin(x))``
    occurs, while ``abs(sin(y))`` and ``abs(sin(3*z))`` are correctly
    simplified, given that `y \in (0,\pi)` and `z \in (-\pi/3,0)`::

        sage: from sage.manifolds.utilities import simplify_abs_trig
        sage: simplify_abs_trig( abs(sin(x)) + abs(sin(y)) + abs(sin(3*z)) )
        abs(sin(x)) + sin(y) + sin(-3*z)

    Note that neither
    :meth:`~sage.symbolic.expression.Expression.simplify_trig` nor
    :meth:`~sage.symbolic.expression.Expression.simplify_full`
    works in this case::

        sage: s = abs(sin(x)) + abs(sin(y)) + abs(sin(3*z))
        sage: s.simplify_trig()
        abs(4*cos(-z)^2 - 1)*abs(sin(-z)) + abs(sin(x)) + abs(sin(y))
        sage: s.simplify_full()
        abs(4*cos(-z)^2 - 1)*abs(sin(-z)) + abs(sin(x)) + abs(sin(y))

    despite the following assumptions hold::

        sage: assumptions()
        [x is real, y is real, y > 0, y < pi, z is real, z > -1/3*pi, z < 0]

    Additional checks are::

        sage: simplify_abs_trig( abs(sin(y/2)) )  # shall simplify
        sin(1/2*y)
        sage: simplify_abs_trig( abs(sin(2*y)) )  # must not simplify
        abs(sin(2*y))
        sage: simplify_abs_trig( abs(sin(z/2)) )  # shall simplify
        sin(-1/2*z)
        sage: simplify_abs_trig( abs(sin(4*z)) )  # must not simplify
        abs(sin(-4*z))

    Simplification of ``abs(cos(...))``::

        sage: forget()
        sage: M = Manifold(3, 'M', structure='topological')
        sage: X.<x,y,z> = M.chart(r'x y:(0,pi/2) z:(pi/4,3*pi/4)')
        sage: X.coord_range()
        x: (-oo, +oo); y: (0, 1/2*pi); z: (1/4*pi, 3/4*pi)
        sage: simplify_abs_trig( abs(cos(x)) + abs(cos(y)) + abs(cos(2*z)) )
        abs(cos(x)) + cos(y) - cos(2*z)

    Additional tests::

        sage: simplify_abs_trig(abs(cos(y-pi/2)))  # shall simplify
        cos(-1/2*pi + y)
        sage: simplify_abs_trig(abs(cos(y+pi/2)))  # shall simplify
        -cos(1/2*pi + y)
        sage: simplify_abs_trig(abs(cos(y-pi)))  # shall simplify
        -cos(-pi + y)
        sage: simplify_abs_trig(abs(cos(2*y)))  # must not simplify
        abs(cos(2*y))
        sage: simplify_abs_trig(abs(cos(y/2)) * abs(sin(z)))  # shall simplify
        cos(1/2*y)*sin(z)

    TESTS:

    Simplification of expressions involving some symbolic derivatives::

        sage: f = function('f')
        sage: s = abs(cos(x)) + abs(cos(y))*diff(f(x),x) + abs(cos(2*z))
        sage: simplify_abs_trig(s)
        cos(y)*diff(f(x), x) + abs(cos(x)) - cos(2*z)
        sage: s = abs(sin(x))*diff(f(x),x).subs(x=y^2) + abs(cos(y))
        sage: simplify_abs_trig(s)
        abs(sin(x))*D[0](f)(y^2) + cos(y)
        sage: forget()  # for doctests below
    """
    w0 = SR.wild()
    if expr.has(abs_symbolic(sin(w0))) or expr.has(abs_symbolic(cos(w0))):
        return SimplifyAbsTrig(expr)()
    return expr


def simplify_chain_real(expr):
    r"""
    Apply a chain of simplifications to a symbolic expression, assuming the
    real domain.

    This is the simplification chain used in calculus involving coordinate
    functions on real manifolds, as implemented in
    :class:`~sage.manifolds.chart_func.ChartFunction`.

    The chain is formed by the following functions, called
    successively:

    #. :meth:`~sage.symbolic.expression.Expression.simplify_factorial`
    #. :meth:`~sage.symbolic.expression.Expression.simplify_trig`
    #. :meth:`~sage.symbolic.expression.Expression.simplify_rational`
    #. :func:`simplify_sqrt_real`
    #. :func:`simplify_abs_trig`
    #. :meth:`~sage.symbolic.expression.Expression.canonicalize_radical`
    #. :meth:`~sage.symbolic.expression.Expression.simplify_log`
    #. :meth:`~sage.symbolic.expression.Expression.simplify_rational`
    #. :meth:`~sage.symbolic.expression.Expression.simplify_trig`

    EXAMPLES:

    We consider variables that are coordinates of a chart on a real manifold::

        sage: M = Manifold(2, 'M', structure='topological')
        sage: X.<x,y> = M.chart('x:(0,1) y')

    The following assumptions then hold::

        sage: assumptions()
        [x is real, x > 0, x < 1, y is real]

    and we have::

        sage: from sage.manifolds.utilities import simplify_chain_real
        sage: s = sqrt(y^2)
        sage: simplify_chain_real(s)
        abs(y)

    The above result is correct since ``y`` is real. It is obtained by
    :meth:`~sage.symbolic.expression.Expression.simplify_real` as well::

        sage: s.simplify_real()
        abs(y)
        sage: s.simplify_full()
        abs(y)

    Furthermore, we have::

        sage: s = sqrt(x^2-2*x+1)
        sage: simplify_chain_real(s)
        -x + 1

    which is correct since `x \in (0,1)`. On this example, neither
    :meth:`~sage.symbolic.expression.Expression.simplify_real`
    nor :meth:`~sage.symbolic.expression.Expression.simplify_full`,
    nor :meth:`~sage.symbolic.expression.Expression.canonicalize_radical`
    give satisfactory results::

        sage: s.simplify_real()  # unsimplified output
        sqrt(x^2 - 2*x + 1)
        sage: s.simplify_full()  # unsimplified output
        sqrt(x^2 - 2*x + 1)
        sage: s.canonicalize_radical()  # wrong output since x in (0,1)
        x - 1

    Other simplifications::

        sage: s = abs(sin(pi*x))
        sage: simplify_chain_real(s)  # correct output since x in (0,1)
        sin(pi*x)
        sage: s.simplify_real()  # unsimplified output
        abs(sin(pi*x))
        sage: s.simplify_full()  # unsimplified output
        abs(sin(pi*x))

    ::

        sage: s = cos(y)^2 + sin(y)^2
        sage: simplify_chain_real(s)
        1
        sage: s.simplify_real()  # unsimplified output
        cos(y)^2 + sin(y)^2
        sage: s.simplify_full()  # OK
        1

    TESTS::

        sage: forget()  # for doctests below
    """
    if expr.number_of_operands() == 0:
        return expr
    expr = expr.simplify_factorial()
    expr = expr.simplify_trig()
    expr = expr.simplify_rational()
    expr = simplify_sqrt_real(expr)
    expr = simplify_abs_trig(expr)
    expr = expr.canonicalize_radical()
    expr = expr.simplify_log('one')
    expr = expr.simplify_rational()
    expr = expr.simplify_trig()
    return expr


def simplify_chain_generic(expr):
    r"""
    Apply a chain of simplifications to a symbolic expression.

    This is the simplification chain used in calculus involving coordinate
    functions on manifolds over fields different from `\RR`, as implemented in
    :class:`~sage.manifolds.chart_func.ChartFunction`.

    The chain is formed by the following functions, called
    successively:

    #. :meth:`~sage.symbolic.expression.Expression.simplify_factorial`
    #. :meth:`~sage.symbolic.expression.Expression.simplify_rectform`
    #. :meth:`~sage.symbolic.expression.Expression.simplify_trig`
    #. :meth:`~sage.symbolic.expression.Expression.simplify_rational`
    #. :meth:`~sage.symbolic.expression.Expression.expand_sum`

    NB: for the time being, this is identical to
    :meth:`~sage.symbolic.expression.Expression.simplify_full`.

    EXAMPLES:

    We consider variables that are coordinates of a chart on a complex
    manifold::

        sage: M = Manifold(2, 'M', structure='topological', field='complex')
        sage: X.<x,y> = M.chart()

    Then neither ``x`` nor ``y`` is assumed to be real::

        sage: assumptions()
        []

    Accordingly, ``simplify_chain_generic`` does not simplify
    ``sqrt(x^2)`` to ``abs(x)``::

        sage: from sage.manifolds.utilities import simplify_chain_generic
        sage: s = sqrt(x^2)
        sage: simplify_chain_generic(s)
        sqrt(x^2)

    This contrasts with the behavior of
    :func:`~sage.manifolds.utilities.simplify_chain_real`.

    Other simplifications::

        sage: s = (x+y)^2 - x^2 -2*x*y - y^2
        sage: simplify_chain_generic(s)
        0
        sage: s = (x^2 - 2*x + 1) / (x^2 -1)
        sage: simplify_chain_generic(s)
        (x - 1)/(x + 1)
        sage: s = cos(2*x) - 2*cos(x)^2 + 1
        sage: simplify_chain_generic(s)
        0

    TESTS::

        sage: forget()  # for doctests below
    """
    if expr.number_of_operands() == 0:
        return expr
    expr = expr.simplify_factorial()
    expr = expr.simplify_rectform()
    expr = expr.simplify_trig()
    expr = expr.simplify_rational()
    expr = expr.expand_sum()
    return expr


def simplify_chain_generic_sympy(expr):
    r"""
    Apply a chain of simplifications to a sympy expression.

    This is the simplification chain used in calculus involving coordinate
    functions on manifolds over fields different from `\RR`, as implemented in
    :class:`~sage.manifolds.chart_func.ChartFunction`.

    The chain is formed by the following functions, called
    successively:

    #. :meth:`~sympy.simplify.combsimp`
    #. :meth:`~sympy.simplify.trigsimp`
    #. :meth:`~sympy.core.expand`
    #. :meth:`~sympy.simplify.simplify`

    EXAMPLES:

    We consider variables that are coordinates of a chart on a complex
    manifold::

        sage: forget()  # for doctest only
        sage: M = Manifold(2, 'M', structure='topological', field='complex', calc_method='sympy')
        sage: X.<x,y> = M.chart()

    Then neither ``x`` nor ``y`` is assumed to be real::

        sage: assumptions()
        []

    Accordingly, ``simplify_chain_generic_sympy`` does not simplify
    ``sqrt(x^2)`` to ``abs(x)``::

        sage: from sage.manifolds.utilities import simplify_chain_generic_sympy
        sage: s = (sqrt(x^2))._sympy_()
        sage: simplify_chain_generic_sympy(s)
        sqrt(x**2)

    This contrasts with the behavior of
    :func:`~sage.manifolds.utilities.simplify_chain_real_sympy`.

    Other simplifications::

        sage: s = ((x+y)^2 - x^2 -2*x*y - y^2)._sympy_()
        sage: simplify_chain_generic_sympy(s)
        0
        sage: s = ((x^2 - 2*x + 1) / (x^2 -1))._sympy_()
        sage: simplify_chain_generic_sympy(s)
        (x - 1)/(x + 1)
        sage: s = (cos(2*x) - 2*cos(x)^2 + 1)._sympy_()
        sage: simplify_chain_generic_sympy(s)
        0
    """
    expr = expr.combsimp()
    expr = expr.trigsimp()
    expr = expr.expand()
    expr = expr.simplify()
    return expr


def simplify_chain_real_sympy(expr):
    r"""
    Apply a chain of simplifications to a sympy expression, assuming the
    real domain.

    This is the simplification chain used in calculus involving coordinate
    functions on real manifolds, as implemented in
    :class:`~sage.manifolds.chart_func.ChartFunction`.

    The chain is formed by the following functions, called
    successively:

    #. :meth:`~sympy.simplify.combsimp`
    #. :meth:`~sympy.simplify.trigsimp`
    #. :func:`simplify_sqrt_real`
    #. :func:`simplify_abs_trig`
    #. :meth:`~sympy.core.expand`
    #. :meth:`~sympy.simplify.simplify`

    EXAMPLES:

    We consider variables that are coordinates of a chart on a real manifold::

        sage: forget()  # for doctest only
        sage: M = Manifold(2, 'M', structure='topological',calc_method='sympy')
        sage: X.<x,y> = M.chart('x:(0,1) y')

    The following assumptions then hold::

        sage: assumptions()
        [x is real, x > 0, x < 1, y is real]

    and we have::

        sage: from sage.manifolds.utilities import simplify_chain_real_sympy
        sage: s = (sqrt(y^2))._sympy_()
        sage: simplify_chain_real_sympy(s)
        Abs(y)

    Furthermore, we have::

        sage: s = (sqrt(x^2-2*x+1))._sympy_()
        sage: simplify_chain_real_sympy(s)
        1 - x

    Other simplifications::

        sage: s = (abs(sin(pi*x)))._sympy_()
        sage: simplify_chain_real_sympy(s)  # correct output since x in (0,1)
        sin(pi*x)

    ::

        sage: s = (cos(y)^2 + sin(y)^2)._sympy_()
        sage: simplify_chain_real_sympy(s)
        1
    """
    # TODO: introduce pure SymPy functions instead of simplify_sqrt_real and
    #       simplify_abs_trig
    if 'sqrt(' in str(expr):
        expr = simplify_sqrt_real(expr._sage_())._sympy_()
    expr = expr.combsimp()
    expr = expr.trigsimp()
    if 'sqrt(' in str(expr):
        expr = simplify_sqrt_real(expr._sage_())._sympy_()
    if 'Abs(sin(' in str(expr):
        expr = simplify_abs_trig(expr._sage_())._sympy_()
    expr = expr.expand()
    expr = expr.simplify()
    return expr

#******************************************************************************


class ExpressionNice(Expression):
    r"""
    Subclass of :class:`~sage.symbolic.expression.Expression` for a
    "human-friendly" display of partial derivatives and the possibility to
    shorten the display by skipping the arguments of symbolic functions.

    INPUT:

    - ``ex`` -- symbolic expression

    EXAMPLES:

    An expression formed with callable symbolic expressions::

        sage: var('x y z')
        (x, y, z)
        sage: f = function('f')(x, y)
        sage: g = f.diff(y).diff(x)
        sage: h = function('h')(y, z)
        sage: k = h.diff(z)
        sage: fun = x*g + y*(k-z)^2

    The standard Pynac display of partial derivatives::

        sage: fun
        y*(z - diff(h(y, z), z))^2 + x*diff(f(x, y), x, y)
        sage: latex(fun)
        y {\left(z - \frac{\partial}{\partial z}h\left(y, z\right)\right)}^{2} + x \frac{\partial^{2}}{\partial x\partial y}f\left(x, y\right)

    With :class:`ExpressionNice`, the Pynac notation ``D[...]`` is replaced
    by textbook-like notation::

        sage: from sage.manifolds.utilities import ExpressionNice
        sage: ExpressionNice(fun)
        y*(z - d(h)/dz)^2 + x*d^2(f)/dxdy
        sage: latex(ExpressionNice(fun))
        y {\left(z - \frac{\partial\,h}{\partial z}\right)}^{2}
         + x \frac{\partial^2\,f}{\partial x\partial y}

    An example when function variables are themselves functions::

        sage: f = function('f')(x, y)
        sage: g = function('g')(x, f)  # the second variable is the function f
        sage: fun = (g.diff(x))*x - x^2*f.diff(x,y)
        sage: fun
        -x^2*diff(f(x, y), x, y) + (diff(f(x, y), x)*D[1](g)(x, f(x, y)) + D[0](g)(x, f(x, y)))*x
        sage: ExpressionNice(fun)
        -x^2*d^2(f)/dxdy + (d(f)/dx*d(g)/d(f(x, y)) + d(g)/dx)*x
        sage: latex(ExpressionNice(fun))
        -x^{2} \frac{\partial^2\,f}{\partial x\partial y}
         + {\left(\frac{\partial\,f}{\partial x}
           \frac{\partial\,g}{\partial \left( f\left(x, y\right) \right)}
         + \frac{\partial\,g}{\partial x}\right)} x

    Note that ``D[1](g)(x, f(x,y))`` is rendered as ``d(g)/d(f(x, y))``.

    An example with multiple differentiations::

        sage: fun = f.diff(x,x,y,y,x)*x
        sage: fun
        x*diff(f(x, y), x, x, x, y, y)
        sage: ExpressionNice(fun)
        x*d^5(f)/dx^3dy^2
        sage: latex(ExpressionNice(fun))
        x \frac{\partial^5\,f}{\partial x ^ 3\partial y ^ 2}

    Parentheses are added around powers of partial derivatives to avoid any
    confusion::

        sage: fun = f.diff(y)^2
        sage: fun
        diff(f(x, y), y)^2
        sage: ExpressionNice(fun)
        (d(f)/dy)^2
        sage: latex(ExpressionNice(fun))
        \left(\frac{\partial\,f}{\partial y}\right)^{2}

    The explicit mention of function arguments can be omitted for the sake of
    brevity::

        sage: fun = fun*f
        sage: ExpressionNice(fun)
        f(x, y)*(d(f)/dy)^2
        sage: Manifold.options.omit_function_arguments=True
        sage: ExpressionNice(fun)
        f*(d(f)/dy)^2
        sage: latex(ExpressionNice(fun))
        f \left(\frac{\partial\,f}{\partial y}\right)^{2}
        sage: Manifold.options._reset()
        sage: ExpressionNice(fun)
        f(x, y)*(d(f)/dy)^2
        sage: latex(ExpressionNice(fun))
        f\left(x, y\right) \left(\frac{\partial\,f}{\partial y}\right)^{2}
    """
    def __init__(self, ex):
        r"""
        Initialize ``self``.

        TESTS::

            sage: f = function('f')(x)
            sage: df = f.diff(x)
            sage: df
            diff(f(x), x)
            sage: from sage.manifolds.utilities import ExpressionNice
            sage: df_nice = ExpressionNice(df)
            sage: df_nice
            d(f)/dx
        """
        from sage.symbolic.ring import SR
        self._parent = SR
        Expression.__init__(self, SR, x=ex)

    def _repr_(self):
        r"""
        String representation of the object.

        EXAMPLES::

            sage: var('x y z')
            (x, y, z)
            sage: f = function('f')(x, y)
            sage: g = f.diff(y).diff(x)
            sage: h = function('h')(y, z)
            sage: k = h.diff(z)
            sage: fun = x*g + y*(k-z)^2
            sage: fun
            y*(z - diff(h(y, z), z))^2 + x*diff(f(x, y), x, y)
            sage: from sage.manifolds.utilities import ExpressionNice
            sage: ExpressionNice(fun)
            y*(z - d(h)/dz)^2 + x*d^2(f)/dxdy

        Check that :issue:`33399` is fixed::

            sage: ExpressionNice(function('f')(x+y, x-y).diff(y))
            d(f)/d(x + y) - d(f)/d(x - y)
        """
        d = self._parent._repr_element_(self)

        # find all occurrences of diff
        list_d = []
        _list_derivatives(self, list_d)

        # process the list
        for m in list_d:
            funcname = m[1]
            diffargs = m[3]
            numargs = len(diffargs)

            if numargs > 1:
                numargs = "^" + str(numargs)
            else:
                numargs = ""

            variables = m[4]
            strv = [str(v) for v in variables]

            # checking if the variable is composite
            comp_chars = ['+', '-', '*', '/', '^', '(']
            for i, sv in enumerate(strv):
                if any(c in sv for c in comp_chars):
                    strv[i] = "(" + sv + ")"

            # dictionary to group multiple occurrences of differentiation: d/dxdx -> d/dx^2 etc.
            occ = dict((i, strv[i] + "^" + str(diffargs.count(i))
                       if (diffargs.count(i) > 1) else strv[i])
                       for i in diffargs)

            res = "d" + str(numargs) + "(" + str(funcname) + ")/d" + "d".join(
                occ.values())

            # str representation of the operator
            s = self._parent._repr_element_(m[0])

            # if diff operator is raised to some power (m[5]), put brackets around
            if m[5]:
                res = "(" + res + ")^" + str(m[5])
                o = s + "^" + str(m[5])
            else:
                o = s

            d = d.replace(o, res)

        import re
        from sage.manifolds.manifold import TopologicalManifold
        if TopologicalManifold.options.omit_function_arguments:
            list_f = []
            _list_functions(self, list_f)

            for m in list_f:
                d = re.sub(m[1] + r'\([^)]+\)', m[1], d)

        return d

    def _latex_(self):
        r"""
        LaTeX representation of the object.

        EXAMPLES::

            sage: var('x y z')
            (x, y, z)
            sage: f = function('f')(x, y)
            sage: g = f.diff(y).diff(x)
            sage: h = function('h')(y, z)
            sage: k = h.diff(z)
            sage: fun = x*g + y*(k-z)^2
            sage: fun
            y*(z - diff(h(y, z), z))^2 + x*diff(f(x, y), x, y)
            sage: from sage.manifolds.utilities import ExpressionNice
            sage: ExpressionNice(fun)
            y*(z - d(h)/dz)^2 + x*d^2(f)/dxdy
            sage: latex(ExpressionNice(fun))
            y {\left(z - \frac{\partial\,h}{\partial z}\right)}^{2} + x \frac{\partial^2\,f}{\partial x\partial y}

        Testing the behavior if no latex_name of the function is given::

            sage: f = function('f_x')(x, y)
            sage: fun = f.diff(y)
            sage: latex(ExpressionNice(fun))
            \frac{\partial\,f_{x}}{\partial y}

        If latex_name, it should be used in LaTeX output::

            sage: f = function('f_x', latex_name=r"{\cal F}")(x,y)
            sage: fun = f.diff(y)
            sage: latex(ExpressionNice(fun))
            \frac{\partial\,{\cal F}}{\partial y}

        Check that :issue:`33399` is fixed::

            sage: latex(ExpressionNice(function('f')(x+y, x-y).diff(y)))
            \frac{\partial\,f}{\partial \left( x + y \right)}
             - \frac{\partial\,f}{\partial \left( x - y \right)}
        """
        from sage.misc.latex import latex

        d = self._parent._latex_element_(self)

        # find all occurrences of diff
        list_d = []
        _list_derivatives(self, list_d)

        for m in list_d:
            if str(m[1]) == str(m[2]):
                funcname = str(m[1])
            else:
                funcname = str(m[2])

            diffargs = m[3]
            numargs = len(diffargs)

            if numargs > 1:
                numargs = "^" + str(numargs)
            else:
                numargs = ""

            variables = m[4]

            strv = [str(v) for v in variables]
            latv = [latex(v) for v in variables]

            # checking if the variable is composite
            comp_chars = ['+', '-', '*', '/', '^', '(']
            for i, sv in enumerate(strv):
                if any(c in sv for c in comp_chars):
                    latv[i] = r"\left(" + latv[i] + r"\right)"

            # dictionary to group multiple occurrences of differentiation: d/dxdx -> d/dx^2 etc.
            occ = {i: (latv[i] + "^" + latex(diffargs.count(i))
                       if diffargs.count(i) > 1 else latv[i])
                   for i in diffargs}

            res = r"\frac{\partial" + numargs + r"\," + funcname + \
                  r"}{\partial " + r"\partial ".join(i for i in occ.values()) + "}"

            # representation of the operator
            s = self._parent._latex_element_(m[0])

            # if diff operator is raised to some power (m[5]), put brackets around
            if m[5]:
                res = r"\left(" + res + r"\right)^{" + str(m[5]) + "}"
                o = s + "^{" + str(m[5]) + "}"
            else:
                o = s

            d = d.replace(o, res)

        from sage.manifolds.manifold import TopologicalManifold
        if TopologicalManifold.options.omit_function_arguments:
            list_f = []
            _list_functions(self, list_f)

            for m in list_f:
                d = d.replace(str(m[3]) + str(m[4]), str(m[3]))

        return d


def _list_derivatives(ex, list_d, exponent=0):
    r"""
    Function to find the occurrences of ``FDerivativeOperator`` in a symbolic
    expression; inspired by
    http://ask.sagemath.org/question/10256/how-can-extract-different-terms-from-a-symbolic-expression/?answer=26136#post-id-26136

    INPUT:

    - ``ex`` -- symbolic expression to be analyzed
    - ``exponent`` -- (optional) exponent of ``FDerivativeOperator``,
      passed to a next level in the expression tree

    OUTPUT:

    - ``list_d`` -- tuple containing the details of ``FDerivativeOperator``
      found, in the following order:

      1. operator
      2. function name
      3. LaTeX function name
      4. parameter set
      5. operands
      6. exponent (if found, else 0)

    TESTS::

        sage: f = function('f_x', latex_name=r"{\cal F}")(x)
        sage: df = f.diff(x)^2
        sage: from sage.manifolds.utilities import _list_derivatives
        sage: list_d = []
        sage: _list_derivatives(df, list_d)
        sage: list_d
        [(diff(f_x(x), x), 'f_x', {\cal F}, [0], [x], 2)]
    """
    op = ex.operator()
    operands = ex.operands()

    import operator
    from sage.misc.latex import latex, latex_variable_name
    from sage.symbolic.operators import FDerivativeOperator

    if op:
        if op is operator.pow:
            if isinstance(operands[0].operator(), FDerivativeOperator):
                exponent = operands[1]

        if isinstance(op, FDerivativeOperator):
            parameter_set = op.parameter_set()
            function = repr(op.function())
            latex_function = latex(op.function())

            # case when no latex_name given
            if function == latex_function:
                latex_function = latex_variable_name(str(op.function()))

            list_d.append((ex, function, latex_function, parameter_set,
                           operands, exponent))

        for operand in operands:
            _list_derivatives(operand, list_d, exponent)


def _list_functions(ex, list_f):
    r"""
    Function to find the occurrences of symbolic functions in a symbolic
    expression.

    INPUT:

    - ``ex`` -- symbolic expression to be analyzed

    OUTPUT:

    - ``list_f`` -- tuple containing the details of a symbolic function found,
      in the following order:

      1. operator
      2. function name
      3. arguments
      4. LaTeX version of function name
      5. LaTeX version of arguments

    TESTS::

        sage: var('x y z')
        (x, y, z)
        sage: f = function('f', latex_name=r"{\cal F}")(x, y)
        sage: g = function('g_x')(x, y)
        sage: d = sin(x)*g.diff(x)*x*f - x^2*f.diff(x,y)/g
        sage: from sage.manifolds.utilities import _list_functions
        sage: list_f = []
        sage: _list_functions(d, list_f)
        sage: list_f
        [(f, 'f', '(x, y)', {\cal F}, \left(x, y\right)),
         (g_x, 'g_x', '(x, y)', 'g_{x}', \left(x, y\right))]
    """
    op = ex.operator()
    operands = ex.operands()

    from sage.misc.latex import latex, latex_variable_name

    if op:
        # FIXME: This hack is needed because the NewSymbolicFunction is
        #   a class defined inside of the *function* function_factory().
        if "NewSymbolicFunction" in str(type(op)):
            repr_function = repr(op)
            latex_function = latex(op)

            # case when no latex_name given
            if repr_function == latex_function:
                latex_function = latex_variable_name(str(op))

            repr_args = repr(ex.arguments())
            # remove comma in case of singleton
            if len(ex.arguments()) == 1:
                repr_args = repr_args.replace(",","")

            latex_args = latex(ex.arguments())

            list_f.append((op, repr_function, repr_args, latex_function, latex_args))

        for operand in operands:
            _list_functions(operand, list_f)

#******************************************************************************


def set_axes_labels(graph, xlabel, ylabel, zlabel, **kwds):
    r"""
    Set axes labels for a 3D graphics object ``graph``.

    This is a workaround for the lack of axes labels in 3D plots.
    This sets the labels as :func:`~sage.plot.plot3d.shapes2.text3d`
    objects at locations determined from the bounding box of the
    graphic object ``graph``.

    INPUT:

    - ``graph`` -- :class:`~sage.plot.plot3d.base.Graphics3d`;
      a 3D graphic object
    - ``xlabel`` -- string for the x-axis label
    - ``ylabel`` -- string for the y-axis label
    - ``zlabel`` -- string for the z-axis label
    - ``**kwds`` -- options (e.g. color) for text3d

    OUTPUT: the 3D graphic object with text3d labels added

    EXAMPLES::

        sage: # needs sage.plot
        sage: g = sphere()
        sage: g.all
        [Graphics3d Object]
        sage: from sage.manifolds.utilities import set_axes_labels
        sage: ga = set_axes_labels(g, 'X', 'Y', 'Z', color='red')
        sage: ga.all  # the 3D frame has now axes labels
        [Graphics3d Object, Graphics3d Object,
         Graphics3d Object, Graphics3d Object]
    """
    from sage.plot.plot3d.shapes2 import text3d
    xmin, ymin, zmin = graph.bounding_box()[0]
    xmax, ymax, zmax = graph.bounding_box()[1]
    dx = xmax - xmin
    dy = ymax - ymin
    dz = zmax - zmin
    x1 = xmin + dx / 2
    y1 = ymin + dy / 2
    z1 = zmin + dz / 2
    xmin1 = xmin - dx / 20
    ymin1 = ymin - dy / 20
    zmin1 = zmin - dz / 20
    graph += text3d('  ' + xlabel, (x1, ymin1, zmin1), **kwds)
    graph += text3d('  ' + ylabel, (xmin1, y1, zmin1), **kwds)
    graph += text3d('  ' + zlabel, (xmin1, ymin1, z1), **kwds)
    return graph


def exterior_derivative(form):
    r"""
    Exterior derivative of a differential form.

    INPUT:

    - ``form`` -- a differential form; this must an instance of either

      * :class:`~sage.manifolds.differentiable.scalarfield.DiffScalarField`
        for a 0-form (scalar field)
      * :class:`~sage.manifolds.differentiable.diff_form.DiffFormParal` for
        a `p`-form (`p\geq 1`) on a parallelizable manifold
      * :class:`~sage.manifolds.differentiable.diff_form.DiffForm` for a
        a `p`-form (`p\geq 1`) on a non-parallelizable manifold

    OUTPUT:

    - the `(p+1)`-form that is the exterior derivative of ``form``

    EXAMPLES:

    Exterior derivative of a scalar field (0-form)::

        sage: from sage.manifolds.utilities import exterior_derivative
        sage: M = Manifold(3, 'M')
        sage: X.<x,y,z> = M.chart()
        sage: f = M.scalar_field({X: x+y^2+z^3}, name='f')
        sage: df = exterior_derivative(f); df
        1-form df on the 3-dimensional differentiable manifold M
        sage: df.display()
        df = dx + 2*y dy + 3*z^2 dz

    An alias is ``xder``::

        sage: from sage.manifolds.utilities import xder
        sage: df == xder(f)
        True

    Exterior derivative of a 1-form::

        sage: a = M.one_form(name='a')
        sage: a[:] = [x+y*z, x-y*z, x*y*z]
        sage: da = xder(a); da
        2-form da on the 3-dimensional differentiable manifold M
        sage: da.display()
        da = (-z + 1) dx∧dy + (y*z - y) dx∧dz + (x*z + y) dy∧dz
        sage: dda = xder(da); dda
        3-form dda on the 3-dimensional differentiable manifold M
        sage: dda.display()
        dda = 0

    .. SEEALSO::

        :class:`sage.manifolds.differentiable.diff_form.DiffFormParal.exterior_derivative`
        or :class:`sage.manifolds.differentiable.diff_form.DiffForm.exterior_derivative`
        for more examples.
    """
    return form.exterior_derivative()


xder = exterior_derivative
