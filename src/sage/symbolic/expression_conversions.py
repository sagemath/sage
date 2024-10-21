# sage_setup: distribution = sagemath-symbolics
"""
Conversion of symbolic expressions to other types

This module provides routines for converting new symbolic expressions
to other types.  Primarily, it provides a class :class:`Converter`
which will walk the expression tree and make calls to methods
overridden by subclasses.
"""
###############################################################################
#   Sage: Open Source Mathematical Software
#       Copyright (C) 2009 Mike Hansen <mhansen@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  https://www.gnu.org/licenses/
###############################################################################

from operator import eq, ne, gt, lt, ge, le, mul, pow, neg, add, truediv
from functools import reduce

import sage.rings.abc

from sage.misc.lazy_import import lazy_import
from sage.symbolic.ring import SR
from sage.structure.element import Expression
from sage.functions.all import exp
from sage.symbolic.operators import arithmetic_operators, relation_operators, FDerivativeOperator, add_vararg, mul_vararg
from sage.rings.number_field.number_field_element_base import NumberFieldElement_base

lazy_import('sage.symbolic.expression_conversion_sympy', ['SympyConverter', 'sympy_converter'])
lazy_import('sage.symbolic.expression_conversion_algebraic', ['AlgebraicConverter', 'algebraic'])

lazy_import('sage.symbolic.expression_conversion_sympy', ['SympyConverter', 'sympy_converter'])
lazy_import('sage.symbolic.expression_conversion_algebraic', ['AlgebraicConverter', 'algebraic'])


class FakeExpression:
    r"""
    Pynac represents `x/y` as `xy^{-1}`.  Often, tree-walkers would prefer
    to see divisions instead of multiplications and negative exponents.
    To allow for this (since Pynac internally doesn't have division at all),
    there is a possibility to pass use_fake_div=True; this will rewrite
    an Expression into a mixture of Expression and FakeExpression nodes,
    where the FakeExpression nodes are used to represent divisions.
    These nodes are intended to act sufficiently like Expression nodes
    that tree-walkers won't care about the difference.
    """

    def __init__(self, operands, operator):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import FakeExpression
            sage: import operator; x,y = var('x,y')
            sage: FakeExpression([x, y], operator.truediv)
            FakeExpression([x, y], <built-in function truediv>)
        """
        self._operands = operands
        self._operator = operator

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import FakeExpression
            sage: import operator; x,y = var('x,y')
            sage: FakeExpression([x, y], operator.truediv)
            FakeExpression([x, y], <built-in function truediv>)
        """
        return "FakeExpression(%r, %r)" % (self._operands, self._operator)

    def pyobject(self):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import FakeExpression
            sage: import operator; x,y = var('x,y')
            sage: f = FakeExpression([x, y], operator.truediv)
            sage: f.pyobject()
            Traceback (most recent call last):
            ...
            TypeError: self must be a numeric expression
        """
        raise TypeError('self must be a numeric expression')

    def operands(self):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import FakeExpression
            sage: import operator; x,y = var('x,y')
            sage: f = FakeExpression([x, y], operator.truediv)
            sage: f.operands()
            [x, y]
        """
        return self._operands

    def __getitem__(self, i):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import FakeExpression
            sage: import operator; x,y = var('x,y')
            sage: f = FakeExpression([x, y], operator.truediv)
            sage: f[0]
            x
        """
        return self._operands[i]

    def operator(self):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import FakeExpression
            sage: import operator; x,y = var('x,y')
            sage: f = FakeExpression([x, y], operator.truediv)
            sage: f.operator()
            <built-in function truediv>
        """
        return self._operator

    def _fast_callable_(self, etb):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import FakeExpression
            sage: import operator; x,y = var('x,y')
            sage: f = FakeExpression([x, y], operator.truediv)
            sage: fast_callable(f, vars=['x','y']).op_list()
            [('load_arg', 0), ('load_arg', 1), 'div', 'return']
        """
        return fast_callable(self, etb)


class Converter:
    def __init__(self, use_fake_div=False):
        """
        If use_fake_div is set to True, then the converter will try to
        replace expressions whose operator is operator.mul with the
        corresponding expression whose operator is operator.truediv.

        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import Converter
            sage: c = Converter(use_fake_div=True)
            sage: c.use_fake_div
            True
        """
        self.use_fake_div = use_fake_div

    def __call__(self, ex=None):
        """
        .. NOTE::

            If this object does not have an attribute ``ex``, then an argument
            must be passed into :meth:`__call__`.

        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import Converter
            sage: c = Converter(use_fake_div=True)
            sage: c(SR(2))
            Traceback (most recent call last):
            ...
            NotImplementedError: pyobject
            sage: c(x+2)
            Traceback (most recent call last):
            ...
            NotImplementedError: arithmetic
            sage: c(x)
            Traceback (most recent call last):
            ...
            NotImplementedError: symbol
            sage: c(x==2)
            Traceback (most recent call last):
            ...
            NotImplementedError: relation
            sage: c(sin(x))
            Traceback (most recent call last):
            ...
            NotImplementedError: composition
            sage: c(function('f')(x).diff(x))
            Traceback (most recent call last):
            ...
            NotImplementedError: derivative

        We can set a default value for the argument by setting
        the ``ex`` attribute::

            sage: c.ex = SR(2)
            sage: c()
            Traceback (most recent call last):
            ...
            NotImplementedError: pyobject
        """
        if ex is None:
            ex = self.ex

        try:
            obj = ex.pyobject()
            return self.pyobject(ex, obj)
        except TypeError as err:
            if 'self must be a numeric expression' not in err.args:
                raise err

        operator = ex.operator()
        if operator is None:
            return self.symbol(ex)

        if operator in arithmetic_operators:
            if getattr(self, 'use_fake_div', False) and (operator is mul or operator is mul_vararg):
                div = self.get_fake_div(ex)
                return self.arithmetic(div, div.operator())
            return self.arithmetic(ex, operator)
        elif operator in relation_operators:
            return self.relation(ex, operator)
        elif isinstance(operator, FDerivativeOperator):
            return self.derivative(ex, operator)
        elif operator == tuple:
            return self.tuple(ex)
        else:
            return self.composition(ex, operator)

    def get_fake_div(self, ex):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import Converter
            sage: c = Converter(use_fake_div=True)
            sage: c.get_fake_div(sin(x)/x)
            FakeExpression([sin(x), x], <built-in function truediv>)
            sage: c.get_fake_div(-1*sin(x))
            FakeExpression([sin(x)], <built-in function neg>)
            sage: c.get_fake_div(-x)
            FakeExpression([x], <built-in function neg>)
            sage: c.get_fake_div((2*x^3+2*x-1)/((x-2)*(x+1)))
            FakeExpression([2*x^3 + 2*x - 1, FakeExpression([x + 1, x - 2], <built-in function mul>)], <built-in function truediv>)

        Check if :issue:`8056` is fixed, i.e., if numerator is 1.::

            sage: c.get_fake_div(1/pi/x)
            FakeExpression([1, FakeExpression([pi, x], <built-in function mul>)], <built-in function truediv>)
        """
        d = []
        n = []
        for arg in ex.operands():
            ops = arg.operands()
            try:
                if arg.operator() is pow and repr(ops[1]) == '-1':
                    d.append(ops[0])
                else:
                    n.append(arg)
            except TypeError:
                n.append(arg)

        len_d = len(d)
        if len_d == 0:
            repr_n = [repr(_) for _ in n]
            if len(n) == 2 and "-1" in repr_n:
                a = n[0] if repr_n[1] == "-1" else n[1]
                return FakeExpression([a], neg)
            else:
                return ex
        elif len_d == 1:
            d = d[0]
        else:
            d = FakeExpression(d, mul)

        if len(n) == 0:
            return FakeExpression([SR.one(), d], truediv)
        elif len(n) == 1:
            n = n[0]
        else:
            n = FakeExpression(n, mul)

        return FakeExpression([n, d], truediv)

    def pyobject(self, ex, obj):
        """
        The input to this method is the result of calling
        :meth:`pyobject` on a symbolic expression.

        .. NOTE::

           Note that if a constant such as ``pi`` is encountered in
           the expression tree, its corresponding pyobject which is an
           instance of :class:`sage.symbolic.constants.Pi` will be
           passed into this method.  One cannot do arithmetic using
           such an object.

        TESTS::

            sage: from sage.symbolic.expression_conversions import Converter
            sage: f = SR(1)
            sage: Converter().pyobject(f, f.pyobject())
            Traceback (most recent call last):
            ...
            NotImplementedError: pyobject
        """
        raise NotImplementedError("pyobject")

    def symbol(self, ex):
        """
        The input to this method is a symbolic expression which
        corresponds to a single variable.  For example, this method
        could be used to return a generator for a polynomial ring.

        TESTS::

            sage: from sage.symbolic.expression_conversions import Converter
            sage: Converter().symbol(x)
            Traceback (most recent call last):
            ...
            NotImplementedError: symbol
        """
        raise NotImplementedError("symbol")

    def relation(self, ex, operator):
        """
        The input to this method is a symbolic expression which
        corresponds to a relation.

        TESTS::

            sage: from sage.symbolic.expression_conversions import Converter
            sage: import operator
            sage: Converter().relation(x==3, operator.eq)
            Traceback (most recent call last):
            ...
            NotImplementedError: relation
            sage: Converter().relation(x==3, operator.lt)
            Traceback (most recent call last):
            ...
            NotImplementedError: relation
        """
        raise NotImplementedError("relation")

    def derivative(self, ex, operator):
        """
        The input to this method is a symbolic expression which
        corresponds to a relation.

        TESTS::

            sage: from sage.symbolic.expression_conversions import Converter
            sage: a = function('f')(x).diff(x); a
            diff(f(x), x)
            sage: Converter().derivative(a, a.operator())
            Traceback (most recent call last):
            ...
            NotImplementedError: derivative
        """
        raise NotImplementedError("derivative")

    def arithmetic(self, ex, operator):
        """
        The input to this method is a symbolic expression and the
        infix operator corresponding to that expression. Typically,
        one will convert all of the arguments and then perform the
        operation afterward.

        TESTS::

            sage: from sage.symbolic.expression_conversions import Converter
            sage: f = x + 2
            sage: Converter().arithmetic(f, f.operator())
            Traceback (most recent call last):
            ...
            NotImplementedError: arithmetic
        """
        raise NotImplementedError("arithmetic")

    def composition(self, ex, operator):
        """
        The input to this method is a symbolic expression and its
        operator.  This method will get called when you have a symbolic
        function application.

        TESTS::

            sage: from sage.symbolic.expression_conversions import Converter
            sage: f = sin(2)
            sage: Converter().composition(f, f.operator())
            Traceback (most recent call last):
            ...
            NotImplementedError: composition
        """
        raise NotImplementedError("composition")


class InterfaceInit(Converter):
    def __init__(self, interface):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import InterfaceInit
            sage: m = InterfaceInit(maxima)
            sage: a = pi + 2
            sage: m(a)
            '(%pi)+(2)'
            sage: m(sin(a))
            'sin((%pi)+(2))'
            sage: m(exp(x^2) + pi + 2)
            '(%pi)+(exp((_SAGE_VAR_x)^(2)))+(2)'
        """
        self.name_init = "_%s_init_" % interface.name()
        self.interface = interface
        self.relation_symbols = interface._relation_symbols()

    def symbol(self, ex):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import InterfaceInit
            sage: m = InterfaceInit(maxima)
            sage: m.symbol(x)
            '_SAGE_VAR_x'
            sage: f(x) = x
            sage: m.symbol(f)
            '_SAGE_VAR_x'
            sage: ii = InterfaceInit(gp)
            sage: ii.symbol(x)
            'x'
            sage: g = InterfaceInit(giac)
            sage: g.symbol(x)
            'sageVARx'
        """
        if self.interface.name() == 'maxima':
            return '_SAGE_VAR_' + repr(SR(ex))
        if self.interface.name() == 'giac':
            return 'sageVAR' + repr(SR(ex))
        return repr(SR(ex))

    def pyobject(self, ex, obj):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import InterfaceInit
            sage: ii = InterfaceInit(gp)
            sage: f = 2+SR(I)
            sage: ii.pyobject(f, f.pyobject())
            'I + 2'

            sage: ii.pyobject(SR(2), 2)
            '2'

            sage: ii.pyobject(pi, pi.pyobject())
            'Pi'
        """
        if (self.interface.name() in ['pari', 'gp'] and isinstance(obj, NumberFieldElement_base)):
            from sage.rings.number_field.number_field_element_quadratic import NumberFieldElement_gaussian
            if isinstance(obj, NumberFieldElement_gaussian):
                return repr(obj)
        try:
            return getattr(obj, self.name_init)()
        except AttributeError:
            return repr(obj)

    def relation(self, ex, operator):
        """
        EXAMPLES::

            sage: import operator
            sage: from sage.symbolic.expression_conversions import InterfaceInit
            sage: m = InterfaceInit(maxima)
            sage: m.relation(x==3, operator.eq)
            '_SAGE_VAR_x = 3'
            sage: m.relation(x==3, operator.lt)
            '_SAGE_VAR_x < 3'
        """
        return "%s %s %s" % (self(ex.lhs()), self.relation_symbols[operator],
                             self(ex.rhs()))

    def tuple(self, ex):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import InterfaceInit
            sage: m = InterfaceInit(maxima)
            sage: t = SR._force_pyobject((3, 4, e^x))
            sage: m.tuple(t)
            '[3,4,exp(_SAGE_VAR_x)]'
        """
        x = map(self, ex.operands())
        X = ','.join(x)
        return str(self.interface._left_list_delim()) + X + str(self.interface._right_list_delim())

    def derivative(self, ex, operator):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import InterfaceInit
            sage: m = InterfaceInit(maxima)
            sage: f = function('f')
            sage: a = f(x).diff(x); a
            diff(f(x), x)
            sage: print(m.derivative(a, a.operator()))
            diff('f(_SAGE_VAR_x), _SAGE_VAR_x, 1)
            sage: b = f(x).diff(x, x)
            sage: print(m.derivative(b, b.operator()))
            diff('f(_SAGE_VAR_x), _SAGE_VAR_x, 2)

        We can also convert expressions where the argument is not just a
        variable, but the result is an "at" expression using temporary
        variables::

            sage: y = var('y')
            sage: t = (f(x*y).diff(x))/y
            sage: t
            D[0](f)(x*y)
            sage: m.derivative(t, t.operator())
            "at(diff('f(_SAGE_VAR__symbol0), _SAGE_VAR__symbol0, 1), [_SAGE_VAR__symbol0 = (_SAGE_VAR_x)*(_SAGE_VAR_y)])"

        TESTS:

        Most of these confirm that :issue:`7401` was fixed::

            sage: t = var('t'); f = function('f')(t)
            sage: a = 2^e^t * f.subs(t=e^t) * diff(f, t).subs(t=e^t) + 2*t
            sage: solve(a == 0, diff(f, t).subs(t=e^t))
            [D[0](f)(e^t) == -2^(-e^t + 1)*t/f(e^t)]

        ::

            sage: f = function('f')(x)
            sage: df = f.diff(x); df
            diff(f(x), x)
            sage: maxima(df)
            'diff('f(_SAGE_VAR_x),_SAGE_VAR_x,1)

        ::

            sage: a = df.subs(x=exp(x)); a
            D[0](f)(e^x)
            sage: b = maxima(a); b
            %at('diff('f(_SAGE_VAR__symbol0),_SAGE_VAR__symbol0,1), _SAGE_VAR__symbol0 = %e^_SAGE_VAR_x)
            sage: bool(b.sage() == a)
            True

        ::

            sage: a = df.subs(x=4); a
            D[0](f)(4)
            sage: b = maxima(a); b
            %at('diff('f(_SAGE_VAR__symbol0),_SAGE_VAR__symbol0,1), _SAGE_VAR__symbol0 = 4)
            sage: bool(b.sage() == a)
            True

        It also works with more than one variable.  Note the preferred
        syntax ``function('f')(x, y)`` to create a general symbolic
        function of more than one variable::

            sage: x, y = var('x y')
            sage: f = function('f')(x, y)
            sage: f_x = f.diff(x); f_x
            diff(f(x, y), x)
            sage: maxima(f_x)
            'diff('f(_SAGE_VAR_x,_SAGE_VAR_y),_SAGE_VAR_x,1)

        ::

            sage: a = f_x.subs(x=4); a
            D[0](f)(4, y)
            sage: b = maxima(a); b
            %at('diff('f(_SAGE_VAR__symbol0,_SAGE_VAR_y),_SAGE_VAR__symbol0,1), _SAGE_VAR__symbol0 = 4)
            sage: bool(b.sage() == a)
            True

        ::

            sage: a = f_x.subs(x=4).subs(y=8); a
            D[0](f)(4, 8)
            sage: b = maxima(a); b
            %at('diff('f(_SAGE_VAR__symbol0,8),_SAGE_VAR__symbol0,1), _SAGE_VAR__symbol0 = 4)
            sage: bool(b.sage() == a)
            True

        Test a special case (:issue:`16697`)::

            sage: x,y = var('x,y')
            sage: (gamma_inc(x,y).diff(x))
            diff(gamma(x, y), x)
            sage: (gamma_inc(x,x+1).diff(x)).simplify()
            -(x + 1)^(x - 1)*e^(-x - 1) + D[0](gamma)(x, x + 1)
        """
        # This code should probably be moved into the interface
        # object in a nice way.
        if self.name_init != "_maxima_init_":
            raise NotImplementedError
        args = ex.operands()
        if (not all(isinstance(v, Expression) and v.is_symbol() for v in args) or
                len(args) != len(set(args))):
            # An evaluated derivative of the form f'(1) is not a
            # symbolic variable, yet we would like to treat it like
            # one. So, we replace the argument `1` with a temporary
            # variable e.g. `_symbol0` and then evaluate the
            # derivative f'(_symbol0) symbolically at _symbol0=1. See
            # trac #12796. Note that we cannot use SR.temp_var here
            # since two conversions of the same expression have to be
            # equal.
            temp_args = [SR.symbol("_symbol%s" % i) for i in range(len(args))]
            f = operator.function()(*temp_args)
            params = operator.parameter_set()
            params = ["%s, %s" % (temp_args[i]._maxima_init_(), params.count(i)) for i in set(params)]
            subs = ["%s = %s" % (t._maxima_init_(), a._maxima_init_())
                    for t, a in zip(temp_args, args)]
            outstr = "at(diff(%s, %s), [%s])" % (f._maxima_init_(),
                                                 ", ".join(params),
                                                 ", ".join(subs))
        else:
            f = operator.function()(*args)
            params = operator.parameter_set()
            params = ["%s, %s" % (args[i]._maxima_init_(), params.count(i))
                      for i in set(params)]
            outstr = "diff(%s, %s)" % (f._maxima_init_(),
                                       ", ".join(params))
        return outstr

    def arithmetic(self, ex, operator):
        """
        EXAMPLES::

            sage: import operator
            sage: from sage.symbolic.expression_conversions import InterfaceInit
            sage: m = InterfaceInit(maxima)
            sage: m.arithmetic(x+2, sage.symbolic.operators.add_vararg)
            '(_SAGE_VAR_x)+(2)'
        """
        args = ["(%s)" % self(op) for op in ex.operands()]
        return arithmetic_operators[operator].join(args)

    def composition(self, ex, operator):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import InterfaceInit
            sage: m = InterfaceInit(maxima)
            sage: m.composition(sin(x), sin)
            'sin(_SAGE_VAR_x)'
            sage: m.composition(ceil(x), ceil)
            'ceiling(_SAGE_VAR_x)'

            sage: m = InterfaceInit(mathematica)
            sage: m.composition(sin(x), sin)
            'Sin[x]'
        """
        ops = ex.operands()
        # FIXME: consider stripping pyobjects() in ops
        if hasattr(operator, self.name_init + "evaled_"):
            return getattr(operator, self.name_init + "evaled_")(*ops)
        else:
            ops = [self(_) for _ in ops]
        try:
            op = getattr(operator, self.name_init)()
        except (TypeError, AttributeError):
            op = repr(operator)

        return self.interface._function_call_string(op, ops, [])


##########
# FriCAS #
##########
class FriCASConverter(InterfaceInit):
    """
    Convert any expression to FriCAS.

    EXAMPLES::

        sage: var('x,y')
        (x, y)
        sage: f = exp(x^2) - arcsin(pi+x)/y
        sage: f._fricas_()                                                      # optional - fricas
             2
            x
        y %e   - asin(x + %pi)
        ----------------------
                   y
    """
    def __init__(self):
        import sage.interfaces.fricas
        super().__init__(sage.interfaces.fricas.fricas)

    def pyobject(self, ex, obj):
        r"""
        Return a string which, when evaluated by FriCAS, returns the
        object as an expression.

        We explicitly add the coercion to the FriCAS domains
        `Expression Integer` and `Expression Complex Integer` to make
        sure that elements of the symbolic ring are translated to
        these.  In particular, this is needed for integration, see
        :issue:`28641` and :issue:`28647`.

        EXAMPLES::

            sage: 2._fricas_().domainOf()                                       # optional - fricas
            PositiveInteger...

            sage: (-1/2)._fricas_().domainOf()                                  # optional - fricas
            Fraction(Integer...)

            sage: SR(2)._fricas_().domainOf()                                   # optional - fricas
            Expression(Integer...)

            sage: (sqrt(2))._fricas_().domainOf()                               # optional - fricas
            Expression(Integer...)

            sage: pi._fricas_().domainOf()                                      # optional - fricas
            Pi...

            sage: asin(pi)._fricas_()                                           # optional - fricas
            asin(%pi)

            sage: I._fricas_().domainOf()                                   # optional - fricas
            Complex(Integer...)

            sage: SR(I)._fricas_().domainOf()                                   # optional - fricas
            Expression(Complex(Integer...))

            sage: ex = (I+sqrt(2)+2)
            sage: ex._fricas_().domainOf()                                      # optional - fricas
            Expression(Complex(Integer...))

            sage: ex._fricas_()^2                                               # optional - fricas
                       +-+
            (4 + 2 %i)\|2  + 5 + 4 %i

            sage: (ex^2)._fricas_()                                             # optional - fricas
                       +-+
            (4 + 2 %i)\|2  + 5 + 4 %i
        """
        try:
            result = getattr(obj, self.name_init)()
        except AttributeError:
            result = repr(obj)
        else:
            if isinstance(obj, NumberFieldElement_base):
                from sage.rings.number_field.number_field_element_quadratic import NumberFieldElement_gaussian
                if isinstance(obj, NumberFieldElement_gaussian):
                    return "((%s)::EXPR COMPLEX INT)" % result
        return "((%s)::EXPR INT)" % result

    def symbol(self, ex):
        """
        Convert the argument, which is a symbol, to FriCAS.

        In this case, we do not return an `Expression Integer`,
        because FriCAS frequently requires elements of domain
        `Symbol` or `Variable` as arguments, for example to
        `integrate`.  Moreover, FriCAS is able to do the conversion
        itself, whenever the argument should be interpreted as a
        symbolic expression.

        EXAMPLES::

            sage: x._fricas_().domainOf()                                       # optional - fricas
            Variable(x)

            sage: (x^2)._fricas_().domainOf()                                   # optional - fricas
            Expression(Integer...)

            sage: (2*x)._fricas_().integrate(x)                                 # optional - fricas
             2
            x
        """
        return repr(ex)

    def derivative(self, ex, operator):
        """
        Convert the derivative of ``self`` in FriCAS.

        INPUT:

        - ``ex`` -- a symbolic expression

        - ``operator`` -- operator

        Note that ``ex.operator() == operator``.

        EXAMPLES::

            sage: var('x,y,z')
            (x, y, z)
            sage: f = function("F")
            sage: f(x)._fricas_()                                               # optional - fricas
            F(x)
            sage: diff(f(x,y,z), x, z, x)._fricas_()                            # optional - fricas
            F      (x,y,z)
             ,1,1,3

        Check that :issue:`25838` is fixed::

            sage: var('x')
            x
            sage: F = function('F')
            sage: integrate(F(x), x, algorithm='fricas')                        # optional - fricas
            integral(F(x), x)

            sage: integrate(diff(F(x), x)*sin(F(x)), x, algorithm='fricas')     # optional - fricas
            -cos(F(x))

        Check that :issue:`27310` is fixed::

            sage: f = function("F")
            sage: var("y")
            y
            sage: ex = (diff(f(x,y), x, x, y)).subs(y=x+y); ex
            D[0, 0, 1](F)(x, x + y)
            sage: fricas(ex)                                                    # optional - fricas
            F      (x,y + x)
             ,1,1,2
        """
        args = ex.operands()  # the arguments the derivative is evaluated at
        params = operator.parameter_set()
        params_set = set(params)
        mult = ",".join(str(params.count(i)) for i in params_set)
        if (not all(isinstance(v, Expression) and v.is_symbol() for v in args) or
                len(args) != len(set(args))):
            # An evaluated derivative of the form f'(1) is not a
            # symbolic variable, yet we would like to treat it like
            # one. So, we replace the argument `1` with a temporary
            # variable e.g. `_symbol0` and then evaluate the
            # derivative f'(_symbol0) symbolically at _symbol0=1. See
            # trac #12796. Note that we cannot use SR.temp_var here
            # since two conversions of the same expression have to be
            # equal.
            temp_args = [SR.symbol("_symbol%s" % i) for i in range(len(args))]
            f = operator.function()(*temp_args)
            vars = ",".join(temp_args[i]._fricas_init_() for i in params_set)
            subs = ",".join("%s = %s" % (t._fricas_init_(), a._fricas_init_())
                            for t, a in zip(temp_args, args))
            outstr = "eval(D(%s, [%s], [%s]), [%s])" % (f._fricas_init_(), vars, mult, subs)
        else:
            f = operator.function()(*args)
            vars = ",".join(args[i]._fricas_init_() for i in params_set)
            outstr = "D(%s, [%s], [%s])" % (f._fricas_init_(), vars, mult)

        return outstr


fricas_converter = FriCASConverter()


##############
# Polynomial #
##############
class PolynomialConverter(Converter):
    def __init__(self, ex, base_ring=None, ring=None):
        """
        A converter from symbolic expressions to polynomials.

        See :func:`polynomial` for details.

        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import PolynomialConverter
            sage: x, y = var('x,y')
            sage: p = PolynomialConverter(x+y, base_ring=QQ)
            sage: p.base_ring
            Rational Field
            sage: p.ring
            Multivariate Polynomial Ring in x, y over Rational Field

            sage: p = PolynomialConverter(x, base_ring=QQ)
            sage: p.base_ring
            Rational Field
            sage: p.ring
            Univariate Polynomial Ring in x over Rational Field

            sage: p = PolynomialConverter(x, ring=QQ['x,y'])
            sage: p.base_ring
            Rational Field
            sage: p.ring
            Multivariate Polynomial Ring in x, y over Rational Field

            sage: p = PolynomialConverter(x+y, ring=QQ['x'])
            Traceback (most recent call last):
            ...
            TypeError: y is not a variable of Univariate Polynomial Ring in x over Rational Field

        TESTS::

            sage: t, x, z = SR.var('t,x,z')
            sage: QQ[i]['x,y,z,t'](4*I*t + 2*x -12*z + 2)
            2*x - 12*z + (4*I)*t + 2
        """
        if not (ring is None or base_ring is None):
            raise TypeError("either base_ring or ring must be specified, but not both")
        self.ex = ex

        if ring is not None:
            base_ring = ring.base_ring()
            self.varnames = ring.variable_names_recursive()
            for v in ex.variables():
                if repr(v) not in self.varnames and v not in base_ring:
                    raise TypeError("%s is not a variable of %s" % (v, ring))
            self.ring = ring
            self.base_ring = base_ring
        elif base_ring is not None:
            self.base_ring = base_ring
            vars = self.ex.variables()
            if len(vars) == 0:
                vars = ['x']
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            self.ring = PolynomialRing(self.base_ring, names=vars)
            self.varnames = self.ring.variable_names()
        else:
            raise TypeError("either a ring or base ring must be specified")

    def symbol(self, ex):
        """
        Return a variable in the polynomial ring.

        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import PolynomialConverter
            sage: p = PolynomialConverter(x, base_ring=QQ)
            sage: p.symbol(x)
            x
            sage: _.parent()
            Univariate Polynomial Ring in x over Rational Field
            sage: y = var('y')
            sage: p = PolynomialConverter(x*y, ring=SR['x'])
            sage: p.symbol(y)
            y
        """
        try:
            # The symbol is one of the polynomial generators
            return self.ring(repr(ex))
        except TypeError:
            # The symbol should go into the base ring
            return self.base_ring(repr(ex))

    def pyobject(self, ex, obj):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import PolynomialConverter
            sage: p = PolynomialConverter(x, base_ring=QQ)
            sage: f = SR(2)
            sage: p.pyobject(f, f.pyobject())
            2
            sage: _.parent()
            Rational Field
        """
        return self.base_ring(obj)

    def composition(self, ex, operator):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import PolynomialConverter
            sage: a = sin(2)
            sage: p = PolynomialConverter(a*x, base_ring=RR)
            sage: p.composition(a, a.operator())
            0.909297426825682
        """
        return self.base_ring(ex)

    def relation(self, ex, op):
        """
        EXAMPLES::

            sage: import operator
            sage: from sage.symbolic.expression_conversions import PolynomialConverter

            sage: x, y = var('x, y')
            sage: p = PolynomialConverter(x, base_ring=RR)

            sage: p.relation(x==3, operator.eq)
            x - 3.00000000000000
            sage: p.relation(x==3, operator.lt)
            Traceback (most recent call last):
            ...
            ValueError: Unable to represent as a polynomial

            sage: p = PolynomialConverter(x - y, base_ring=QQ)
            sage: p.relation(x^2 - y^3 + 1 == x^3, operator.eq)
            -x^3 - y^3 + x^2 + 1
        """
        import operator
        if op == operator.eq:
            return self(ex.lhs()) - self(ex.rhs())
        raise ValueError("Unable to represent as a polynomial")

    def arithmetic(self, ex, operator):
        """
        EXAMPLES::

            sage: import operator
            sage: from sage.symbolic.expression_conversions import PolynomialConverter

            sage: x, y = var('x, y')
            sage: p = PolynomialConverter(x, base_ring=RR)
            sage: p.arithmetic(pi+e, operator.add)
            5.85987448204884
            sage: p.arithmetic(x^2, operator.pow)
            x^2

            sage: p = PolynomialConverter(x+y, base_ring=RR)
            sage: p.arithmetic(x*y+y^2, operator.add)
            x*y + y^2

            sage: p = PolynomialConverter(y^(3/2), ring=SR['x'])
            sage: p.arithmetic(y^(3/2), operator.pow)
            y^(3/2)
            sage: _.parent()
            Symbolic Ring
        """
        if not any(repr(v) in self.varnames for v in ex.variables()):
            return self.base_ring(ex)
        elif operator == pow:
            from sage.rings.integer import Integer
            base, exp = ex.operands()
            return self(base)**Integer(exp)
        if operator == add_vararg:
            operator = add
        elif operator == mul_vararg:
            operator = mul
        ops = [self(a) for a in ex.operands()]
        return reduce(operator, ops)


def polynomial(ex, base_ring=None, ring=None):
    """
    Return a polynomial from the symbolic expression ``ex``.

    INPUT:

    - ``ex`` -- a symbolic expression

    - ``base_ring``, ``ring`` -- either a
      ``base_ring`` or a polynomial ``ring`` can be
      specified for the parent of result.
      If just a ``base_ring`` is given, then the variables
      of the ``base_ring`` will be the variables of the expression ``ex``.

    OUTPUT: a polynomial

    EXAMPLES::

         sage: from sage.symbolic.expression_conversions import polynomial
         sage: f = x^2 + 2
         sage: polynomial(f, base_ring=QQ)
         x^2 + 2
         sage: _.parent()
         Univariate Polynomial Ring in x over Rational Field

         sage: polynomial(f, ring=QQ['x,y'])
         x^2 + 2
         sage: _.parent()
         Multivariate Polynomial Ring in x, y over Rational Field

         sage: x, y = var('x, y')
         sage: polynomial(x + y^2, ring=QQ['x,y'])
         y^2 + x
         sage: _.parent()
         Multivariate Polynomial Ring in x, y over Rational Field

         sage: s,t = var('s,t')
         sage: expr = t^2-2*s*t+1
         sage: expr.polynomial(None,ring=SR['t'])
         t^2 - 2*s*t + 1
         sage: _.parent()
         Univariate Polynomial Ring in t over Symbolic Ring

         sage: polynomial(x*y, ring=SR['x'])
         y*x

         sage: polynomial(y - sqrt(x), ring=SR['y'])
         y - sqrt(x)
         sage: _.list()
         [-sqrt(x), 1]

    The polynomials can have arbitrary (constant) coefficients so long as
    they coerce into the base ring::

         sage: polynomial(2^sin(2)*x^2 + exp(3), base_ring=RR)
         1.87813065119873*x^2 + 20.0855369231877
    """
    converter = PolynomialConverter(ex, base_ring=base_ring, ring=ring)
    res = converter()
    return converter.ring(res)


class LaurentPolynomialConverter(PolynomialConverter):
    def __init__(self, ex, base_ring=None, ring=None):
        """
        A converter from symbolic expressions to Laurent polynomials.

        See :func:`laurent_polynomial` for details.

        TESTS::

            sage: from sage.symbolic.expression_conversions import LaurentPolynomialConverter
            sage: x, y = var('x,y')
            sage: p = LaurentPolynomialConverter(x+1/y, base_ring=QQ)
            sage: p.base_ring
            Rational Field
            sage: p.ring
            Multivariate Laurent Polynomial Ring in x, y over Rational Field
        """
        super().__init__(ex, base_ring, ring)

        if ring is None and base_ring is not None:
            from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
            self.ring = LaurentPolynomialRing(self.base_ring,
                                              names=self.varnames)


def laurent_polynomial(ex, base_ring=None, ring=None):
    """
    Return a Laurent polynomial from the symbolic expression ``ex``.

    INPUT:

    - ``ex`` -- a symbolic expression

    - ``base_ring``, ``ring`` -- either a
      ``base_ring`` or a Laurent polynomial ``ring`` can be
      specified for the parent of result.
      If just a ``base_ring`` is given, then the variables
      of the ``base_ring`` will be the variables of the expression ``ex``.

    OUTPUT: a Laurent polynomial

    EXAMPLES::

         sage: from sage.symbolic.expression_conversions import laurent_polynomial
         sage: f = x^2 + 2/x
         sage: laurent_polynomial(f, base_ring=QQ)
         2*x^-1 + x^2
         sage: _.parent()
         Univariate Laurent Polynomial Ring in x over Rational Field

         sage: laurent_polynomial(f, ring=LaurentPolynomialRing(QQ, 'x, y'))
         x^2 + 2*x^-1
         sage: _.parent()
         Multivariate Laurent Polynomial Ring in x, y over Rational Field

         sage: x, y = var('x, y')
         sage: laurent_polynomial(x + 1/y^2, ring=LaurentPolynomialRing(QQ, 'x, y'))
         x + y^-2
         sage: _.parent()
         Multivariate Laurent Polynomial Ring in x, y over Rational Field
    """
    converter = LaurentPolynomialConverter(ex, base_ring=base_ring, ring=ring)
    res = converter()
    return converter.ring(res)


#################
# Fast Callable #
#################

class FastCallableConverter(Converter):
    def __init__(self, ex, etb):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import FastCallableConverter
            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=['x'])
            sage: f = FastCallableConverter(x+2, etb)
            sage: f.ex
            x + 2
            sage: f.etb
            <sage.ext.fast_callable.ExpressionTreeBuilder object at 0x...>
            sage: f.use_fake_div
            True
        """
        self.ex = ex
        self.etb = etb
        Converter.__init__(self, use_fake_div=True)

    def pyobject(self, ex, obj):
        r"""
        EXAMPLES::

            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=['x'])
            sage: pi._fast_callable_(etb)
            pi
            sage: etb = ExpressionTreeBuilder(vars=['x'], domain=RDF)
            sage: pi._fast_callable_(etb)
            3.141592653589793
        """
        from sage.symbolic.constants import Constant
        if isinstance(obj, Constant):
            obj = obj.expression()
        return self.etb.constant(obj)

    def relation(self, ex, operator):
        """
        EXAMPLES::

            sage: ff = fast_callable(x == 2, vars=['x'])
            sage: ff(2)
            0
            sage: ff(4)
            2
            sage: ff = fast_callable(x < 2, vars=['x'])
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if operator is not eq:
            raise NotImplementedError
        return self(ex.lhs() - ex.rhs())

    def arithmetic(self, ex, operator):
        r"""
        EXAMPLES::

            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=['x','y'])
            sage: var('x,y')
            (x, y)
            sage: (x+y)._fast_callable_(etb)
            add(v_0, v_1)
            sage: (-x)._fast_callable_(etb)
            neg(v_0)
            sage: (x+y+x^2)._fast_callable_(etb)
            add(add(ipow(v_0, 2), v_0), v_1)

        TESTS:

        Check if rational functions with numerator 1 can
        be converted. (:issue:`8056`)::

            sage: (1/pi/x)._fast_callable_(etb)
            div(1, mul(pi, v_0))

            sage: etb = ExpressionTreeBuilder(vars=['x'], domain=RDF)
            sage: (x^7)._fast_callable_(etb)
            ipow(v_0, 7)
            sage: f(x)=1/pi/x; plot(f,2,3)
            Graphics object consisting of 1 graphics primitive
        """
        # This used to convert the operands first.  Doing it this way
        # instead gives a chance to notice powers with an integer
        # exponent before the exponent gets (potentially) converted
        # to another type.
        operands = ex.operands()
        if operator is pow:
            exponent = operands[1]
            if exponent == -1:
                return self.etb.call(truediv, 1, operands[0])
            elif exponent == 0.5:
                from sage.misc.functional import sqrt
                return self.etb.call(sqrt, operands[0])
            elif exponent == -0.5:
                from sage.misc.functional import sqrt
                return self.etb.call(truediv, 1, self.etb.call(sqrt, operands[0]))
        elif operator is neg:
            return self.etb.call(operator, operands[0])
        if operator == add_vararg:
            operator = add
        elif operator == mul_vararg:
            operator = mul
        return reduce(lambda x, y: self.etb.call(operator, x, y), operands)

    def symbol(self, ex):
        r"""
        Given an ExpressionTreeBuilder, return an Expression representing
        this value.

        EXAMPLES::

            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=['x','y'])
            sage: x, y, z = var('x,y,z')
            sage: x._fast_callable_(etb)
            v_0
            sage: y._fast_callable_(etb)
            v_1
            sage: z._fast_callable_(etb)
            Traceback (most recent call last):
            ...
            ValueError: Variable 'z' not found...
        """
        return self.etb.var(SR(ex))

    def composition(self, ex, function):
        r"""
        Given an ExpressionTreeBuilder, return an Expression representing
        this value.

        EXAMPLES::

            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=['x','y'])
            sage: x,y = var('x,y')
            sage: sin(sqrt(x+y))._fast_callable_(etb)
            sin(sqrt(add(v_0, v_1)))
            sage: arctan2(x,y)._fast_callable_(etb)
            {arctan2}(v_0, v_1)
        """
        return self.etb.call(function, *ex.operands())

    def tuple(self, ex):
        r"""
        Given a symbolic tuple, return its elements as a Python list.

        EXAMPLES::

            sage: from sage.ext.fast_callable import ExpressionTreeBuilder
            sage: etb = ExpressionTreeBuilder(vars=['x'])
            sage: SR._force_pyobject((2, 3, x^2))._fast_callable_(etb)
            [2, 3, x^2]
        """
        return ex.operands()


def fast_callable(ex, etb):
    """
    Given an ExpressionTreeBuilder *etb*, return an Expression representing
    the symbolic expression *ex*.

    EXAMPLES::

        sage: from sage.ext.fast_callable import ExpressionTreeBuilder
        sage: etb = ExpressionTreeBuilder(vars=['x','y'])
        sage: x,y = var('x,y')
        sage: f = y+2*x^2
        sage: f._fast_callable_(etb)
        add(mul(ipow(v_0, 2), 2), v_1)

        sage: f = (2*x^3+2*x-1)/((x-2)*(x+1))
        sage: f._fast_callable_(etb)
        div(add(add(mul(ipow(v_0, 3), 2), mul(v_0, 2)), -1), mul(add(v_0, 1), add(v_0, -2)))
    """
    return FastCallableConverter(ex, etb)()


class RingConverter(Converter):
    def __init__(self, R, subs_dict=None):
        """
        A class to convert expressions to other rings.

        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import RingConverter
            sage: R = RingConverter(RIF, subs_dict={x:2})
            sage: R.ring
            Real Interval Field with 53 bits of precision
            sage: R.subs_dict
            {x: 2}
            sage: R(pi+e)
            5.85987448204884?
            sage: loads(dumps(R))
            <sage.symbolic.expression_conversions.RingConverter object at 0x...>
        """
        self.subs_dict = {} if subs_dict is None else subs_dict
        self.ring = R

    def symbol(self, ex):
        """
        All symbols appearing in the expression must either appear in
        *subs_dict* or be convertible by the ring's element
        constructor in order for the conversion to be successful.

        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import RingConverter
            sage: R = RingConverter(RIF, subs_dict={x:2})
            sage: R(x+pi)
            5.141592653589794?

            sage: R = RingConverter(RIF)
            sage: R(x+pi)
            Traceback (most recent call last):
            ...
            TypeError: unable to simplify to a real interval approximation

            sage: R = RingConverter(QQ['x'])
            sage: R(x^2+x)
            x^2 + x
            sage: R(x^2+x).parent()
            Univariate Polynomial Ring in x over Rational Field
        """
        try:
            return self.ring(self.subs_dict[ex])
        except KeyError:
            return self.ring(ex)

    def pyobject(self, ex, obj):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import RingConverter
            sage: R = RingConverter(RIF)
            sage: R(SR(5/2))
            2.5000000000000000?
        """
        return self.ring(obj)

    def arithmetic(self, ex, operator):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import RingConverter
            sage: P.<z> = ZZ[]
            sage: R = RingConverter(P, subs_dict={x:z})
            sage: a = 2*x^2 + x + 3
            sage: R(a)
            2*z^2 + z + 3
        """
        if operator not in [pow, add_vararg, mul_vararg]:
            raise TypeError

        operands = ex.operands()
        if operator is pow:
            from sage.rings.integer import Integer
            from sage.rings.rational import Rational
            base, expt = operands

            if expt == Rational((1, 2)):
                from sage.misc.functional import sqrt
                return sqrt(self(base))
            try:
                expt = Integer(expt)
            except TypeError:
                pass

            base = self(base)
            return base ** expt

        if operator == add_vararg:
            operator = add
        elif operator == mul_vararg:
            operator = mul
        return reduce(operator, map(self, operands))

    def composition(self, ex, operator):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import RingConverter
            sage: R = RingConverter(RIF)
            sage: R(cos(2))
            -0.4161468365471424?
        """
        res = operator(*[self(op) for op in ex.operands()])
        if res.parent() is not self.ring:
            raise TypeError
        else:
            return res


class ExpressionTreeWalker(Converter):
    def __init__(self, ex):
        """
        A class that walks the tree. Mainly for subclassing.

        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import ExpressionTreeWalker
            sage: from sage.symbolic.random_tests import random_expr
            sage: ex = sin(atan(0,hold=True)+hypergeometric((1,),(1,),x))
            sage: s = ExpressionTreeWalker(ex)
            sage: bool(s() == ex)
            True
            sage: set_random_seed(0)  # random_expr is unstable
            sage: foo = random_expr(20, nvars=2)
            sage: s = ExpressionTreeWalker(foo)
            sage: bool(s() == foo)
            True
        """
        self.ex = ex

    def symbol(self, ex):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import ExpressionTreeWalker
            sage: s = ExpressionTreeWalker(x)
            sage: bool(s.symbol(x) == x)
            True
        """
        return ex

    def pyobject(self, ex, obj):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import ExpressionTreeWalker
            sage: f = SR(2)
            sage: s = ExpressionTreeWalker(f)
            sage: bool(s.pyobject(f, f.pyobject()) == f.pyobject())
            True
        """
        return ex

    def relation(self, ex, operator):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import ExpressionTreeWalker
            sage: foo = function('foo')
            sage: eq = foo(x) == x
            sage: s = ExpressionTreeWalker(eq)
            sage: s.relation(eq, eq.operator()) == eq
            True
        """
        return operator(self(ex.lhs()), self(ex.rhs()))

    def arithmetic(self, ex, operator):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import ExpressionTreeWalker
            sage: foo = function('foo')
            sage: f = x*foo(x) + pi/foo(x)
            sage: s = ExpressionTreeWalker(f)
            sage: bool(s.arithmetic(f, f.operator()) == f)
            True
        """
        return reduce(operator, map(self, ex.operands()))

    def composition(self, ex, operator):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import ExpressionTreeWalker
            sage: foo = function('foo')
            sage: f = foo(atan2(0, 0, hold=True))
            sage: s = ExpressionTreeWalker(f)
            sage: bool(s.composition(f, f.operator()) == f)
            True
        """
        from sage.symbolic.function import Function
        if isinstance(operator, Function):
            return operator(*map(self, ex.operands()), hold=True)
        else:
            return operator(*map(self, ex.operands()))

    def derivative(self, ex, operator):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import ExpressionTreeWalker
            sage: foo = function('foo')
            sage: f = foo(x).diff(x)
            sage: s = ExpressionTreeWalker(f)
            sage: bool(s.derivative(f, f.operator()) == f)
            True
        """
        return operator(*map(self, ex.operands()))

    def tuple(self, ex):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import ExpressionTreeWalker
            sage: foo = function('foo')
            sage: f = hypergeometric((1,2,3,),(x,),x)
            sage: s = ExpressionTreeWalker(f)
            sage: bool(s() == f)
            True
        """
        return ex.operands()


class SubstituteFunction(ExpressionTreeWalker):
    def __init__(self, ex, *args):
        """
        A class that walks the tree and replaces occurrences of a
        function with another.

        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import SubstituteFunction
            sage: foo = function('foo'); bar = function('bar')
            sage: s = SubstituteFunction(foo(x), {foo: bar})
            sage: s(1/foo(foo(x)) + foo(2))
            1/bar(bar(x)) + bar(2)

        TESTS:

        Check that the old syntax still works::

            sage: s = SubstituteFunction(foo(x), foo, bar)
            sage: s(1/foo(foo(x)) + foo(2))
            1/bar(bar(x)) + bar(2)
        """
        if len(args) == 2:
            self.substitutions = {args[0]: args[1]}
        elif len(args) == 1:
            self.substitutions = args[0]
        else:
            raise TypeError('SubstituteFunction takes either one or two arguments.')
        self.ex = ex

    def composition(self, ex, operator):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import SubstituteFunction
            sage: foo = function('foo'); bar = function('bar')
            sage: s = SubstituteFunction(foo(x), {foo: bar})
            sage: f = foo(x)
            sage: s.composition(f, f.operator())
            bar(x)
            sage: f = foo(foo(x))
            sage: s.composition(f, f.operator())
            bar(bar(x))
            sage: f = sin(foo(x))
            sage: s.composition(f, f.operator())
            sin(bar(x))
            sage: f = foo(sin(x))
            sage: s.composition(f, f.operator())
            bar(sin(x))
        """
        new = self.substitutions.get(operator)
        if new is not None:
            return new(*[self(_) for _ in ex.operands()])
        else:
            return super().composition(ex, operator)

    def derivative(self, ex, operator):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import SubstituteFunction
            sage: foo = function('foo'); bar = function('bar')
            sage: s = SubstituteFunction(foo(x), {foo: bar})
            sage: f = foo(x).diff(x)
            sage: s.derivative(f, f.operator())
            diff(bar(x), x)

        TESTS:

        We can substitute functions under a derivative operator,
        :issue:`12801`::

            sage: f = function('f')
            sage: g = function('g')
            sage: f(g(x)).diff(x).substitute_function({g: sin})
            cos(x)*D[0](f)(sin(x))
        """
        new = self.substitutions.get(operator.function())
        if new is not None:
            return operator.change_function(new)(*[self(_) for _ in ex.operands()])
        else:
            return operator(*[self(_) for _ in ex.operands()])


class Exponentialize(ExpressionTreeWalker):
    # Implementation note: this code is executed once at first
    # reference in the code using it, therefore avoiding rebuilding
    # the same canned results dictionary at each call.
    from sage.calculus.var import function
    from sage.functions.hyperbolic import sinh, cosh, sech, csch, tanh, coth
    from sage.functions.log import exp
    from sage.functions.trig import sin, cos, sec, csc, tan, cot
    from sage.rings.integer import Integer
    from sage.symbolic.constants import e, I
    from sage.symbolic.ring import SR
    half = Integer(1) / Integer(2)
    two = Integer(2)
    x = SR.var("x")
    CircDict = {
        sin: (-half*I*exp(I*x) + half*I*exp(-I*x)).function(x),
        cos: (half*exp(I*x) + half*exp(-I*x)).function(x),
        sec: (two/(exp(I*x) + exp(-I*x))).function(x),
        csc: (two*I/(exp(I*x) - exp(-I*x))).function(x),
        tan: (-I*(exp(I*x) - exp(-I*x))/(exp(I*x) + exp(-I*x))).function(x),
        cot: (I*(exp(I*x) + exp(-I*x))/(exp(I*x) - exp(-I*x))).function(x),
        sinh: (-half*exp(-x) + half*exp(x)).function(x),
        cosh: (half*exp(-x) + half*exp(x)).function(x),
        sech: (two/(exp(-x) + exp(x))).function(x),
        csch: (-two/(exp(-x) - exp(x))).function(x),
        tanh: (-(exp(-x) - exp(x))/(exp(x) + exp(-x))).function(x),
        coth: (-(exp(-x) + exp(x))/(exp(-x) - exp(x))).function(x)
    }
    Circs = list(CircDict)

    def __init__(self, ex):
        """
        A class that walks a symbolic expression tree and replace circular
        and hyperbolic functions by their respective exponential
        expressions.

        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import Exponentialize
            sage: d = Exponentialize(sin(x))
            sage: d(sin(x))
            -1/2*I*e^(I*x) + 1/2*I*e^(-I*x)
            sage: d(cosh(x))
            1/2*e^(-x) + 1/2*e^x
        """
        self.ex = ex

    def composition(self, ex, op):
        r"""
        Return the composition of ``self`` with ``ex`` by ``op``.

        EXAMPLES::

            sage: x = SR.var("x")
            sage: from sage.symbolic.expression_conversions import Exponentialize
            sage: p = x
            sage: s = Exponentialize(p)
            sage: q = sin(x)
            sage: s.composition(q, q.operator())
            -1/2*I*e^(I*x) + 1/2*I*e^(-I*x)
        """
        if op in self.Circs:
            return self.CircDict.get(op)(*[self(oper)
                                           for oper in ex.operands()])
        return super().composition(ex, op)


class DeMoivre(ExpressionTreeWalker):
    def __init__(self, ex, force=False):
        r"""
        A class that walks a symbolic expression tree and replaces
        occurences of complex exponentials (optionally, all
        exponentials) by their respective trigonometric expressions.

        INPUT:

        - ``force`` -- boolean (default: ``False``); replace `\exp(x)`
          with `\cosh(x) + \sinh(x)`

        EXAMPLES::

            sage: a, b = SR.var("a, b")
            sage: from sage.symbolic.expression_conversions import DeMoivre
            sage: d = DeMoivre(e^a)
            sage: d(e^(a+I*b))
            (cos(b) + I*sin(b))*e^a
        """
        self.ex = ex
        self.force = force

    def composition(self, ex, op):
        """
        Return the composition of ``self`` with ``ex`` by ``op``.

        EXAMPLES::

            sage: x, a, b = SR.var('x, a, b')
            sage: from sage.symbolic.expression_conversions import DeMoivre
            sage: p = exp(x)
            sage: s = DeMoivre(p)
            sage: q = exp(a+I*b)
            sage: s.composition(q, q.operator())
            (cos(b) + I*sin(b))*e^a
        """
        from sage.functions.log import exp
        if op is not exp:
            # return super().composition(ex, op)
            return op(*[self(oper) for oper in ex.operands()])

        from sage.rings.imaginary_unit import I
        from sage.symbolic.ring import SR
        from sage.functions.hyperbolic import sinh, cosh
        from sage.functions.trig import sin, cos
        arg = self(ex.operands()[0])()
        w0, w1 = (SR.wild(u) for u in range(2))
        D = arg.match(w0 + I*w1)
        if D is not None:
            A = D.get(w1)
            return exp(D.get(w0))*(cos(A) + I*sin(A))
        D = arg.match(I*w0)
        if D is not None:
            A = D.get(w0)
            return cos(A) + I*sin(A)
        if self.force:
            return cosh(arg) + sinh(arg)
        return exp(arg)


# Half_angle transformation. Sometimes useful in integration

class HalfAngle(ExpressionTreeWalker):
    """
    A class that walks a symbolic expression tree, replacing each
    occurrence of a trigonometric or hyperbolic function by its
    expression as a rational fraction in the (hyperbolic) tangent
    of half the original argument.
    """
    # Code executed once at first class reference: create canned formulae.
    from sage.calculus.var import function
    from sage.functions.hyperbolic import sinh, cosh, sech, csch, tanh, coth
    from sage.functions.trig import sin, cos, sec, csc, tan, cot
    from sage.rings.integer import Integer
    from sage.symbolic.ring import SR
    x = SR.var("x")
    one = Integer(1)
    two = Integer(2)
    half = one / two
    halfx = half * x
    HalvesDict = {
        sin: two * tan(halfx) / (tan(halfx)**2 + one).function(x),
        cos: -(tan(halfx)**2 - one) / (tan(halfx)**2 + one).function(x),
        tan: -two * tan(halfx) / (tan(halfx)**2 - one).function(x),
        csc: half * (tan(halfx)**2 + one) / tan(halfx).function(x),
        sec: -(tan(halfx)**2 + one) / (tan(halfx)**2 - one).function(x),
        cot: -half * (tan(halfx)**2 - one) / tan(halfx).function(x),
        sinh: -two * tanh(halfx) / (tanh(halfx)**2 - one).function(x),
        cosh: -(tanh(halfx)**2 + one) / (tanh(halfx)**2 - one).function(x),
        tanh: two * tanh(halfx) / (tanh(halfx)**2 + one).function(x),
        csch: -half * (tanh(halfx)**2 - one) / tanh(halfx).function(x),
        sech: -(tanh(halfx)**2 - one) / (tanh(halfx)**2 + one).function(x),
        coth: half * (tanh(halfx)**2 + one) / tanh(halfx).function(x)
    }
    Halves = list(HalvesDict)

    def __init__(self, ex):
        """
        A class that walks a symbolic expression tree, replacing each
        occurrence of a trigonometric or hyperbolic function by its
        expression as a rational fraction in the (hyperbolic) tangent
        of half the original argument.

        EXAMPLES::

            sage: a, b = SR.var("a, b")
            sage: from sage.symbolic.expression_conversions import HalfAngle
            sage: HalfAngle(tan(a))(tan(a)+4)
            -2*tan(1/2*a)/(tan(1/2*a)^2 - 1) + 4
        """
        self.ex = ex

    def composition(self, ex, op):
        """
        Compose.

        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import HalfAngle
            sage: x, t = SR.var("x, t")
            sage: a = HalfAngle(cos(3*x)/(4-cos(x)).trig_expand())()
            sage: a.subs(tan(x/2) == t).simplify_full()
            (2*(t^2 + 1)*cos(3/2*x)^2 - t^2 - 1)/(5*t^2 + 3)
        """
        if op in self.Halves:
            return self.HalvesDict.get(op)(*[self(x) for x in ex.operands()])
        return super().composition(ex, op)


class HoldRemover(ExpressionTreeWalker):
    def __init__(self, ex, exclude=None):
        """
        A class that walks the tree and evaluates every operator
        that is not in a given list of exceptions.

        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import HoldRemover
            sage: ex = sin(pi*cos(0, hold=True), hold=True); ex
            sin(pi*cos(0))
            sage: h = HoldRemover(ex)
            sage: h()
            0
            sage: h = HoldRemover(ex, [sin])
            sage: h()
            sin(pi)
            sage: h = HoldRemover(ex, [cos])
            sage: h()
            sin(pi*cos(0))
            sage: ex = atan2(0, 0, hold=True) + hypergeometric([1,2], [3,4], 0, hold=True)
            sage: h = HoldRemover(ex, [atan2])
            sage: h()
            arctan2(0, 0) + 1
            sage: h = HoldRemover(ex, [hypergeometric])
            sage: h()
            NaN + hypergeometric((1, 2), (3, 4), 0)
        """
        self.ex = ex
        if exclude is None:
            exclude = []
        self._exclude = exclude

    def composition(self, ex, operator):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import HoldRemover
            sage: ex = sin(pi*cos(0, hold=True), hold=True); ex
            sin(pi*cos(0))
            sage: h = HoldRemover(ex)
            sage: h()
            0
        """
        from sage.calculus.calculus import symbolic_sum, symbolic_product
        from sage.functions.other import Function_sum, Function_prod
        if not operator:
            return self
        if isinstance(operator, Function_sum):
            return symbolic_sum(*map(self, ex.operands()))
        if isinstance(operator, Function_prod):
            return symbolic_product(*map(self, ex.operands()))
        if operator in self._exclude:
            return operator(*map(self, ex.operands()), hold=True)
        else:
            return operator(*map(self, ex.operands()))
