"Operators"

import operator
from sage.symbolic.ring import is_SymbolicVariable, SR

def add_vararg(first,*rest):
    r"""
    Addition of a variable number of arguments.

    INPUT:

    - ``first``, ``rest`` - arguments to add

    OUTPUT: sum of arguments

    EXAMPLES::

        sage: from sage.symbolic.operators import add_vararg
        sage: add_vararg(1,2,3,4,5,6,7)
        28
        sage: F=(1+x+x^2)
        sage: bool(F.operator()(*F.operands()) == F)
        True
    """

    for r in rest:
        first = first + r
    return first

def mul_vararg(first,*rest):
    r"""
    Multiplication of a variable number of arguments.

    INPUT:

    - ``args`` - arguments to multiply

    OUTPUT: product of arguments

    EXAMPLES::

        sage: from sage.symbolic.operators import mul_vararg
        sage: mul_vararg(9,8,7,6,5,4)
        60480
        sage: G=x*cos(x)*sin(x)
        sage: bool(G.operator()(*G.operands())==G)
        True
    """

    for r in rest:
        first = first * r
    return first

arithmetic_operators = {add_vararg: '+',
                        mul_vararg: '*',
                        operator.add: '+',
                        operator.sub: '-',
                        operator.mul: '*',
                        operator.truediv: '/',
                        operator.floordiv: '//',
                        operator.pow: '^'}

relation_operators = {operator.eq:'==',
                      operator.lt:'<',
                      operator.gt:'>',
                      operator.ne:'!=',
                      operator.le:'<=',
                      operator.ge:'>='}

class FDerivativeOperator(object):
    def __init__(self, function, parameter_set):
        """
        EXAMPLES::

            sage: from sage.symbolic.operators import FDerivativeOperator
            sage: f = function('foo')
            sage: op = FDerivativeOperator(f, [0,1])
            sage: loads(dumps(op))
            D[0, 1](foo)
        """
        self._f = function
        self._parameter_set = [int(_) for _ in parameter_set]

    def __call__(self, *args):
        """
        EXAMPLES::

            sage: from sage.symbolic.operators import FDerivativeOperator
            sage: x,y = var('x,y')
            sage: f = function('foo')
            sage: op = FDerivativeOperator(f, [0,1])
            sage: op(x,y)
            diff(foo(x, y), x, y)
            sage: op(x,x^2)
            D[0, 1](foo)(x, x^2)

        TESTS:

        We should be able to operate on functions evaluated at a
        point, not just a symbolic variable, :trac:`12796`::

           sage: from sage.symbolic.operators import FDerivativeOperator
           sage: f = function('f')
           sage: op = FDerivativeOperator(f, [0])
           sage: op(1)
           D[0](f)(1)

        """
        if (not all(is_SymbolicVariable(x) for x in args) or
                len(args) != len(set(args))):
            # An evaluated derivative of the form f'(1) is not a
            # symbolic variable, yet we would like to treat it
            # like one. So, we replace the argument `1` with a
            # temporary variable e.g. `t0` and then evaluate the
            # derivative f'(t0) symbolically at t0=1. See trac
            # #12796.
            temp_args=SR.temp_var(n=len(args))
            vars=[temp_args[i] for i in self._parameter_set]
            return self._f(*temp_args).diff(*vars).function(*temp_args)(*args)
        vars = [args[i] for i in self._parameter_set]
        return self._f(*args).diff(*vars)

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.symbolic.operators import FDerivativeOperator
            sage: f = function('foo')
            sage: op = FDerivativeOperator(f, [0,1]); op
            D[0, 1](foo)
        """
        return "D[%s](%s)"%(", ".join(map(repr, self._parameter_set)), self._f)

    def function(self):
        """
        EXAMPLES::

            sage: from sage.symbolic.operators import FDerivativeOperator
            sage: f = function('foo')
            sage: op = FDerivativeOperator(f, [0,1])
            sage: op.function()
            foo
        """
        return self._f

    def change_function(self, new):
        """
        Returns a new FDerivativeOperator with the same parameter set
        for a new function.

            sage: from sage.symbolic.operators import FDerivativeOperator
            sage: f = function('foo')
            sage: b = function('bar')
            sage: op = FDerivativeOperator(f, [0,1])
            sage: op.change_function(bar)
            D[0, 1](bar)
        """
        return FDerivativeOperator(new, self._parameter_set)

    def parameter_set(self):
        """
        EXAMPLES::

            sage: from sage.symbolic.operators import FDerivativeOperator
            sage: f = function('foo')
            sage: op = FDerivativeOperator(f, [0,1])
            sage: op.parameter_set()
            [0, 1]
        """
        return self._parameter_set
