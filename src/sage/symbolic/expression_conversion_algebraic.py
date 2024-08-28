r"""
Conversion of symbolic expressions to algebraic numbers
"""
# ****************************************************************************
#       Copyright (C) 2009-2012 Mike Hansen
#                     2015-2018 Ralf Stephan
#                     2015      Nils Bruin
#                     2017      Jeroen Demeyer
#                     2019-2022 Frédéric Chapoton
#                     2021      Dave Witte Morris
#                     2023      Vincent Delecroix
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from operator import eq, ne, gt, lt, ge, le, mul, pow, neg, add, truediv
from functools import reduce

import sage.rings.abc

from sage.functions.all import exp
from sage.symbolic.expression_conversions import Converter
from sage.symbolic.operators import add_vararg, mul_vararg
from sage.symbolic.ring import SR
from sage.rings.qqbar import AlgebraicRealField


#############
# Algebraic #
#############
class AlgebraicConverter(Converter):
    def __init__(self, field):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import AlgebraicConverter
            sage: a = AlgebraicConverter(QQbar)
            sage: a.field
            Algebraic Field
            sage: a.reciprocal_trig_functions['cot']
            tan
        """
        self.field = field

        from sage.functions.all import reciprocal_trig_functions
        self.reciprocal_trig_functions = reciprocal_trig_functions

    def pyobject(self, ex, obj):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import AlgebraicConverter
            sage: a = AlgebraicConverter(QQbar)
            sage: f = SR(2)
            sage: a.pyobject(f, f.pyobject())
            2
            sage: _.parent()
            Algebraic Field
        """
        return self.field(obj)

    def arithmetic(self, ex, operator):
        """
        Convert a symbolic expression to an algebraic number.

        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import AlgebraicConverter
            sage: f = 2^(1/2)
            sage: a = AlgebraicConverter(QQbar)
            sage: a.arithmetic(f, f.operator())
            1.414213562373095?

        Note that converting an expression where an odd root is taken will take the real root
        if the target is ``AA``, but this behavior will change in the future::

            sage: AA((-1)^(2/3))
            doctest:warning...
            DeprecationWarning: Taking the root of an algebraic real number will yield the principal root in the future.
            See https://github.com/sagemath/sage/issues/38362 for details.
            1

        TESTS::

            sage: f = pi^6
            sage: a = AlgebraicConverter(QQbar)
            sage: a.arithmetic(f, f.operator())
            Traceback (most recent call last):
            ...
            TypeError: unable to convert pi^6 to Algebraic Field

        Test that :issue:`14602` is fixed::

            sage: K = QuadraticField(3)
            sage: K(sqrt(3)).parent() is K
            True
            sage: sqrt(K(3)).parent() is K
            True
            sage: (K(3)^(1/2)).parent()
            Symbolic Ring
            sage: bool(K.gen() == K(3)^(1/2) == sqrt(K(3)) == K(sqrt(3)) == sqrt(3))
            True
            sage: L = QuadraticField(3, embedding=-AA(3).sqrt())
            sage: bool(L.gen() == -sqrt(3))
            True
        """
        # We try to avoid simplifying, because maxima's simplify command
        # can change the value of a radical expression (by changing which
        # root is selected).
        try:
            if operator is pow:
                from sage.rings.rational import Rational
                base, expt = ex.operands()
                base = self.field(base)
                expt = Rational(expt)
                return self.field(base**expt)
            else:
                if operator is add_vararg:
                    operator = add
                elif operator is mul_vararg:
                    operator = mul
                return reduce(operator, map(self, ex.operands()))
        except TypeError:
            pass

        if operator is pow:
            from sage.symbolic.constants import e, pi, I
            from sage.rings.rational_field import QQ

            base, expt = ex.operands()
            if base == e and expt / (pi * I) in QQ:
                return exp(expt)._algebraic_(self.field)

        raise TypeError("unable to convert %r to %s" % (ex, self.field))

    def composition(self, ex, operator):
        """
        Coerce to an algebraic number.

        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import AlgebraicConverter
            sage: a = AlgebraicConverter(QQbar)
            sage: a.composition(exp(I*pi/3, hold=True), exp)
            0.500000000000000? + 0.866025403784439?*I
            sage: a.composition(sin(pi/7), sin)
            0.4338837391175581? + 0.?e-18*I

            sage: x = SR.var('x')
            sage: a.composition(complex_root_of(x^3 - x^2 - x - 1, 0), complex_root_of)
            1.839286755214161?
            sage: a.composition(complex_root_of(x^5 - 1, 3), complex_root_of)
            0.3090169943749474? - 0.9510565162951536?*I
            sage: a.composition(complex_root_of(x^2 + 1, 0), complex_root_of)
            1.?e-683 - 0.9999999999999999?*I
            sage: a.composition(complex_root_of(x^2 + 1, 1), complex_root_of)
            1.?e-683 + 0.9999999999999999?*I

        TESTS::

            sage: QQbar(zeta(7))
            Traceback (most recent call last):
            ...
            TypeError: unable to convert zeta(7) to Algebraic Field

        Test :issue:`22571`::

            sage: a.composition(exp(0, hold=True), exp)
            1
            sage: a.composition(exp(1, hold=True), exp)
            Traceback (most recent call last):
            ...
            ValueError: unable to represent as an algebraic number
            sage: a.composition(exp(pi*I*RR(1), hold=True), exp)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert e^(1.00000000000000*I*pi) to Algebraic Field
            sage: a.composition(exp(pi*CC.gen(), hold=True), exp)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert e^(1.00000000000000*I*pi) to Algebraic Field
            sage: bool(sin(pi*RR("0.7000000000000002")) > 0)
            True

        Check that :issue:`24440` is fixed::

            sage: QQbar(tanh(pi + 0.1))
            Traceback (most recent call last):
            ...
            ValueError: unable to represent as an algebraic number
            sage: QQbar(sin(I*pi/7))
            Traceback (most recent call last):
            ...
            ValueError: unable to represent as an algebraic number
            sage: QQbar(sin(I*pi/7, hold=True))
            Traceback (most recent call last):
            ...
            ValueError: unable to represent as an algebraic number
        """
        func = operator
        operands = ex.operands()
        if len(operands) == 1:
            operand = operands[0]
        else:
            operand = None

        if isinstance(self.field, sage.rings.abc.UniversalCyclotomicField):
            QQbar = self.field
            hold = True
        else:
            QQbar = self.field.algebraic_closure()
            hold = False

        zeta = QQbar.zeta
        # Note that comparing functions themselves goes via maxima, and is SLOW
        func_name = repr(func)
        if func_name == 'exp':
            if operand.is_trivial_zero():
                return self.field.one()
            if not (SR(-1).sqrt() * operand).is_real():
                raise ValueError("unable to represent as an algebraic number")
            # Coerce (not convert, see #22571) arg to a rational
            from sage.rings.rational_field import QQ
            arg = operand.imag()/(2*ex.parent().pi())
            try:
                rat_arg = QQ.coerce(arg.pyobject())
            except TypeError:
                raise TypeError("unable to convert %r to %s" % (ex, self.field))
            res = zeta(rat_arg.denom())**rat_arg.numer()
            return self.field(res)
        elif func_name in ['sin', 'cos', 'tan']:
            exp_ia = exp(SR(-1).sqrt() * operand, hold=hold)._algebraic_(QQbar)
            if func_name == 'sin':
                res = (exp_ia - ~exp_ia) / (2 * zeta(4))
            elif func_name == 'cos':
                res = (exp_ia + ~exp_ia) / 2
            else:
                res = -zeta(4) * (exp_ia - ~exp_ia) / (exp_ia + ~exp_ia)
            return self.field(res)
        elif func_name in ['sinh', 'cosh', 'tanh']:
            if not (SR(-1).sqrt()*operand).is_real():
                raise ValueError("unable to represent as an algebraic number")
            exp_a = exp(operand, hold=hold)._algebraic_(QQbar)
            if func_name == 'sinh':
                res = (exp_a - ~exp_a) / 2
            elif func_name == 'cosh':
                res = (exp_a + ~exp_a) / 2
            else:
                res = (exp_a - ~exp_a) / (exp_a + ~exp_a)
            return self.field(res)
        elif func_name in self.reciprocal_trig_functions:
            res = ~self.reciprocal_trig_functions[func_name](operand)._algebraic_(QQbar)
            return self.field(res)
        elif func_name == 'complex_root_of':
            cr = ex._sympy_()
            poly = cr.poly._sage_()
            interval = cr._get_interval()._sage_()
            return self.field.polynomial_root(poly, interval)
        elif operand is not None:
            res = func(operand._algebraic_(self.field))
            # We have to handle the case where we get the same symbolic
            # expression back.  For example, QQbar(zeta(7)).  See
            # issue #12665.
            if (res - ex).is_trivial_zero():
                raise TypeError("unable to convert %r to %s" % (ex, self.field))
            return self.field(res)

        raise ValueError("unable to represent as an algebraic number")


def algebraic(ex, field):
    """
    Return the symbolic expression ``ex`` as a element of the algebraic
    field ``field``.

    EXAMPLES::

        sage: a = SR(5/6)
        sage: AA(a)
        5/6
        sage: type(AA(a))
        <class 'sage.rings.qqbar.AlgebraicReal'>
        sage: QQbar(a)
        5/6
        sage: type(QQbar(a))
        <class 'sage.rings.qqbar.AlgebraicNumber'>
        sage: QQbar(i)
        I
        sage: AA(golden_ratio)
        1.618033988749895?
        sage: QQbar(golden_ratio)
        1.618033988749895?
        sage: QQbar(sin(pi/3))
        0.866025403784439?

        sage: QQbar(sqrt(2) + sqrt(8))
        4.242640687119285?
        sage: AA(sqrt(2) ^ 4) == 4
        True
        sage: AA(-golden_ratio)
        -1.618033988749895?
        sage: QQbar((2*SR(I))^(1/2))
        1 + 1*I
        sage: QQbar(e^(pi*I/3))
        0.50000000000000000? + 0.866025403784439?*I

        sage: AA(x*sin(0))
        0
        sage: QQbar(x*sin(0))
        0
    """
    return AlgebraicConverter(field)(ex)
