r"""
Conversion of symbolic expressions to SymPy
"""

from operator import eq, ne, gt, lt, ge, le, mul, pow, neg, add, truediv

from sage.structure.element import Expression
from sage.symbolic.expression_conversions import Converter
from sage.symbolic.operators import arithmetic_operators

#########
# Sympy #
#########
class SympyConverter(Converter):
    """
    Converts any expression to SymPy.

    EXAMPLES::

        sage: import sympy                                                              # optional - sympy
        sage: var('x,y')
        (x, y)
        sage: f = exp(x^2) - arcsin(pi+x)/y
        sage: f._sympy_()                                                               # optional - sympy
        exp(x**2) - asin(x + pi)/y
        sage: _._sage_()                                                                # optional - sympy
        -arcsin(pi + x)/y + e^(x^2)

        sage: sympy.sympify(x) # indirect doctest                                       # optional - sympy
        x

    TESTS:

    Make sure we can convert I (:trac:`6424`)::

        sage: bool(I._sympy_() == I)                                                    # optional - sympy
        True
        sage: (x+I)._sympy_()                                                           # optional - sympy
        x + I

    """
    def __init__(self):
        """
        TESTS::

            sage: from sage.symbolic.expression_conversions import SympyConverter       # optional - sympy
            sage: s = SympyConverter()  # indirect doctest                              # optional - sympy
            sage: TestSuite(s).run(skip="_test_pickling")                               # optional - sympy
        """
        from sage.interfaces.sympy import sympy_init
        sympy_init()

    def __call__(self, ex=None):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import SympyConverter       # optional - sympy
            sage: s = SympyConverter()                                                  # optional - sympy
            sage: f(x, y) = x^2 + y^2; f
            (x, y) |--> x^2 + y^2
            sage: s(f)                                                                  # optional - sympy
            Lambda((x, y), x**2 + y**2)
        """
        if isinstance(ex, Expression) and ex.is_callable():
            from sympy import Symbol, Lambda
            return Lambda(tuple(Symbol(str(arg)) for arg in ex.arguments()),
                          super().__call__(ex))
        return super().__call__(ex)

    def pyobject(self, ex, obj):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import SympyConverter       # optional - sympy
            sage: s = SympyConverter()                                                  # optional - sympy
            sage: f = SR(2)
            sage: s.pyobject(f, f.pyobject())                                           # optional - sympy
            2
            sage: type(_)                                                               # optional - sympy
            <class 'sympy.core.numbers.Integer'>
        """
        try:
            return obj._sympy_()
        except AttributeError:
            return obj

    def arithmetic(self, ex, operator):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import SympyConverter       # optional - sympy
            sage: s = SympyConverter()                                                  # optional - sympy
            sage: f = x + 2
            sage: s.arithmetic(f, f.operator())                                         # optional - sympy
            x + 2
        """
        import sympy
        operator = arithmetic_operators[operator]
        ops = [sympy.sympify(self(a), evaluate=False) for a in ex.operands()]
        if operator == "+":
            return sympy.Add(*ops)
        elif operator == "*":
            return sympy.Mul(*ops)
        elif operator == "-":
            return sympy.Sub(*ops)
        elif operator == "/":
            return sympy.Div(*ops)
        elif operator == "^":
            return sympy.Pow(*ops)
        else:
            raise NotImplementedError

    def symbol(self, ex):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import SympyConverter       # optional - sympy
            sage: s = SympyConverter()                                                  # optional - sympy
            sage: s.symbol(x)                                                           # optional - sympy
            x
            sage: type(_)                                                               # optional - sympy
            <class 'sympy.core.symbol.Symbol'>
        """
        import sympy
        return sympy.symbols(repr(ex))

    def relation(self, ex, op):
        """
        EXAMPLES::

            sage: import operator                                                       # optional - sympy
            sage: from sage.symbolic.expression_conversions import SympyConverter       # optional - sympy
            sage: s = SympyConverter()                                                  # optional - sympy
            sage: s.relation(x == 3, operator.eq)                                       # optional - sympy
            Eq(x, 3)
            sage: s.relation(pi < 3, operator.lt)                                       # optional - sympy
            pi < 3
            sage: s.relation(x != pi, operator.ne)                                      # optional - sympy
            Ne(x, pi)
            sage: s.relation(x > 0, operator.gt)                                        # optional - sympy
            x > 0
        """
        from sympy import Eq, Ne, Gt, Lt, Ge, Le
        ops = {eq: Eq, ne: Ne, gt: Gt, lt: Lt, ge: Ge, le: Le}
        return ops.get(op)(self(ex.lhs()), self(ex.rhs()), evaluate=False)

    def composition(self, ex, operator):
        """
        EXAMPLES::

            sage: from sage.symbolic.expression_conversions import SympyConverter       # optional - sympy
            sage: s = SympyConverter()                                                  # optional - sympy
            sage: f = sin(2)
            sage: s.composition(f, f.operator())                                        # optional - sympy
            sin(2)
            sage: type(_)                                                               # optional - sympy
            sin
            sage: f = arcsin(2)
            sage: s.composition(f, f.operator())                                        # optional - sympy
            asin(2)
        """
        g = ex.operands()
        try:
            return operator._sympy_(*g)
        except (AttributeError, TypeError):
            pass
        f = operator._sympy_init_()
        import sympy

        f_sympy = getattr(sympy, f, None)
        if f_sympy:
            return f_sympy(*sympy.sympify(g, evaluate=False))
        else:
            return sympy.Function(str(f))(*g, evaluate=False)

    def tuple(self, ex):
        """
        Conversion of tuples.

        EXAMPLES::

            sage: t = SR._force_pyobject((3, 4, e^x))
            sage: t._sympy_()                                                           # optional - sympy
            (3, 4, e^x)
            sage: t = SR._force_pyobject((cos(x),))
            sage: t._sympy_()                                                           # optional - sympy
            (cos(x),)

        TESTS::

            sage: from sage.symbolic.expression_conversions import sympy_converter      # optional - sympy
            sage: F = hypergeometric([1/3,2/3],[1,1],x)
            sage: F._sympy_()                                                           # optional - sympy
            hyper((1/3, 2/3), (1, 1), x)

            sage: F = hypergeometric([1/3,2/3],[1],x)
            sage: F._sympy_()                                                           # optional - sympy
            hyper((1/3, 2/3), (1,), x)

            sage: var('a,b,c,d')
            (a, b, c, d)
            sage: hypergeometric((a,b,),(c,),d)._sympy_()                               # optional - sympy
            hyper((a, b), (c,), d)
        """
        return tuple(ex.operands())

    def derivative(self, ex, operator):
        """
        Convert the derivative of ``self`` in sympy.

        INPUT:

        - ``ex`` -- a symbolic expression

        - ``operator`` -- operator

        TESTS::

            sage: var('x','y')
            (x, y)

            sage: f_sage = function('f_sage')(x, y)
            sage: f_sympy = f_sage._sympy_()                                            # optional - sympy

            sage: df_sage = f_sage.diff(x, 2, y, 1); df_sage
            diff(f_sage(x, y), x, x, y)
            sage: df_sympy = df_sage._sympy_(); df_sympy                                # optional - sympy
            Derivative(f_sage(x, y), (x, 2), y)
            sage: df_sympy == f_sympy.diff(x, 2, y, 1)                                  # optional - sympy
            True

        Check that :trac:`28964` is fixed::

            sage: f = function('f')
            sage: _ = var('x,t')
            sage: diff(f(x, t), x)._sympy_(), diff(f(x, t), t)._sympy_()                # optional - sympy
            (Derivative(f(x, t), x), Derivative(f(x, t), t))

        Check differentiating by variables with multiple occurrences
        (:trac:`28964`)::

            sage: f = function('f')
            sage: _ = var('x1,x2,x3,x,t')
            sage: f(x, x, t).diff(x)._sympy_()._sage_()                                 # optional - sympy
            D[0](f)(x, x, t) + D[1](f)(x, x, t)

            sage: g = f(x1, x2, x3, t).diff(x1, 2, x2).subs(x1==x, x2==x, x3==x); g
            D[0, 0, 1](f)(x, x, x, t)
            sage: g._sympy_()                                                           # optional - sympy
            Subs(Derivative(f(_xi_1, _xi_2, x, t), (_xi_1, 2), _xi_2),
                 (_xi_1, _xi_2), (x, x))
            sage: assert g._sympy_()._sage_() == g                                      # optional - sympy

        Check that the use of dummy variables does not cause a collision::

            sage: f = function('f')
            sage: _ = var('x1,x2,x,xi_1')
            sage: g = f(x1, x2, xi_1).diff(x1).subs(x1==x, x2==x); g
            D[0](f)(x, x, xi_1)
            sage: assert g._sympy_()._sage_() == g                                      # optional - sympy
        """
        import sympy

        # retrieve derivated function
        f = operator.function()

        # retrieve order
        order = operator._parameter_set
        # arguments
        _args = [a._sympy_() for a in ex.operands()]

        # when differentiating by a variable that occurs multiple times,
        # substitute it by a dummy variable
        subs_new = []
        subs_old = []
        sympy_arg = []
        for idx in order:
            a = _args[idx]
            if _args.count(a) > 1:
                D = sympy.Dummy('xi_%i' % (idx + 1))
                # to avoid collisions with ordinary symbols when converting
                # back to Sage, we pick an unused variable name for the dummy
                while D._sage_() in ex.variables():
                    D = sympy.Dummy(D.name + '_0')
                subs_old.append(a)
                subs_new.append(D)
                _args[idx] = D
                sympy_arg.append(D)
            else:
                sympy_arg.append(a)

        f_sympy = f._sympy_()(*_args)
        result = f_sympy.diff(*sympy_arg)
        if subs_new:
            return sympy.Subs(result, subs_new, subs_old)
        else:
            return result


sympy_converter = SympyConverter()
