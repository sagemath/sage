# sage_setup: distribution = sagemath-modules
r"""
Constructor for symbolic expressions
"""

from sage.misc.lazy_import import lazy_import
from sage.structure.element import Matrix
from sage.symbolic.expression import Expression
from sage.symbolic.ring import SR

lazy_import('sage.modules.free_module_element', ['vector', 'FreeModuleElement'])
lazy_import('sage.matrix.constructor', 'matrix')


def symbolic_expression(x):
    """
    Create a symbolic expression or vector of symbolic expressions from x.

    INPUT:

    - ``x`` - an object

    OUTPUT: a symbolic expression

    EXAMPLES::

        sage: a = symbolic_expression(3/2); a
        3/2
        sage: type(a)
        <class 'sage.symbolic.expression.Expression'>
        sage: R.<x> = QQ[]; type(x)
        <class 'sage.rings.polynomial.polynomial_rational_flint.Polynomial_rational_flint'>
        sage: a = symbolic_expression(2*x^2 + 3); a
        2*x^2 + 3
        sage: type(a)
        <class 'sage.symbolic.expression.Expression'>
        sage: from sage.structure.element import Expression
        sage: isinstance(a, Expression)
        True
        sage: a in SR
        True
        sage: a.parent()
        Symbolic Ring

    Note that equations exist in the symbolic ring::

        sage: # needs sage.schemes
        sage: E = EllipticCurve('15a'); E
        Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 10*x - 10 over Rational Field
        sage: symbolic_expression(E)
        x*y + y^2 + y == x^3 + x^2 - 10*x - 10
        sage: symbolic_expression(E) in SR
        True

    If ``x`` is a list or tuple, create a vector of symbolic expressions::

        sage: v = symbolic_expression([x,1]); v
        (x, 1)
        sage: v.base_ring()
        Symbolic Ring
        sage: v = symbolic_expression((x,1)); v
        (x, 1)
        sage: v.base_ring()
        Symbolic Ring
        sage: v = symbolic_expression((3,1)); v
        (3, 1)
        sage: v.base_ring()
        Symbolic Ring
        sage: E = EllipticCurve('15a'); E
        Elliptic Curve defined by y^2 + x*y + y = x^3 + x^2 - 10*x - 10 over Rational Field
        sage: v = symbolic_expression([E,E]); v
        (x*y + y^2 + y == x^3 + x^2 - 10*x - 10, x*y + y^2 + y == x^3 + x^2 - 10*x - 10)
        sage: v.base_ring()
        Symbolic Ring

    Likewise, if ``x`` is a vector, create a vector of symbolic expressions::

        sage: u = vector([1, 2, 3])
        sage: v = symbolic_expression(u); v
        (1, 2, 3)
        sage: v.parent()
        Vector space of dimension 3 over Symbolic Ring

    If ``x`` is a list or tuple of lists/tuples/vectors, create a matrix of symbolic expressions::

        sage: M = symbolic_expression([[1, x, x^2], (x, x^2, x^3), vector([x^2, x^3, x^4])]); M
        [  1   x x^2]
        [  x x^2 x^3]
        [x^2 x^3 x^4]
        sage: M.parent()
        Full MatrixSpace of 3 by 3 dense matrices over Symbolic Ring

    If ``x`` is a matrix, create a matrix of symbolic expressions::

        sage: A = matrix([[1, 2, 3], [4, 5, 6]])
        sage: B = symbolic_expression(A); B
        [1 2 3]
        [4 5 6]
        sage: B.parent()
        Full MatrixSpace of 2 by 3 dense matrices over Symbolic Ring

    If ``x`` is a function, for example defined by a ``lambda`` expression, create a
    symbolic function::

        sage: f = symbolic_expression(lambda z: z^2 + 1); f
        z |--> z^2 + 1
        sage: f.parent()
        Callable function ring with argument z
        sage: f(7)
        50

    If ``x`` is a list or tuple of functions, or if ``x`` is a function that returns a list
    or tuple, create a callable symbolic vector::

        sage: symbolic_expression([lambda mu, nu: mu^2 + nu^2, lambda mu, nu: mu^2 - nu^2])
        (mu, nu) |--> (mu^2 + nu^2, mu^2 - nu^2)
        sage: f = symbolic_expression(lambda uwu: [1, uwu, uwu^2]); f
        uwu |--> (1, uwu, uwu^2)
        sage: f.parent()
        Vector space of dimension 3 over Callable function ring with argument uwu
        sage: f(5)
        (1, 5, 25)
        sage: f(5).parent()
        Vector space of dimension 3 over Symbolic Ring

    TESTS:

    Lists, tuples, and vectors of length 0 become vectors over a symbolic ring::

        sage: symbolic_expression([]).parent()
        Vector space of dimension 0 over Symbolic Ring
        sage: symbolic_expression(()).parent()
        Vector space of dimension 0 over Symbolic Ring
        sage: symbolic_expression(vector(QQ, 0)).parent()
        Vector space of dimension 0 over Symbolic Ring

    If a matrix has dimension 0, the result is still a matrix over a symbolic ring::

        sage: symbolic_expression(matrix(QQ, 2, 0)).parent()
        Full MatrixSpace of 2 by 0 dense matrices over Symbolic Ring
        sage: symbolic_expression(matrix(QQ, 0, 3)).parent()
        Full MatrixSpace of 0 by 3 dense matrices over Symbolic Ring

    Also functions defined using ``def`` can be used, but we do not advertise it as a use case::

        sage: def sos(x, y):
        ....:     return x^2 + y^2
        sage: symbolic_expression(sos)
        (x, y) |--> x^2 + y^2

    Functions that take a varying number of arguments or keyword-only arguments are not accepted::

        sage: def variadic(x, *y):
        ....:     return x
        sage: symbolic_expression(variadic)
        Traceback (most recent call last):
        ...
        TypeError: unable to convert <function variadic at 0x...> to a symbolic expression

        sage: def function_with_keyword_only_arg(x, *, sign=1):
        ....:     return sign * x
        sage: symbolic_expression(function_with_keyword_only_arg)
        Traceback (most recent call last):
        ...
        TypeError: unable to convert <function function_with_keyword_only_arg at 0x...>
        to a symbolic expression
    """
    if isinstance(x, Expression):
        return x
    elif hasattr(x, '_symbolic_'):
        return x._symbolic_(SR)
    elif isinstance(x, (tuple, list, FreeModuleElement)):
        expressions = [symbolic_expression(item) for item in x]
        if not expressions:
            # Make sure it is symbolic also when length is 0
            return vector(SR, 0)
        if isinstance(expressions[0], FreeModuleElement):
            return matrix(expressions)
        return vector(expressions)
    elif isinstance(x, Matrix):
        if not x.nrows() or not x.ncols():
            # Make sure it is symbolic and of correct dimensions
            # also when a matrix dimension is 0
            return matrix(SR, x.nrows(), x.ncols())
        rows = [symbolic_expression(row) for row in x.rows()]
        return matrix(rows)
    elif callable(x):
        from inspect import signature, Parameter
        try:
            s = signature(x)
        except ValueError:
            pass
        else:
            if all(param.kind in (Parameter.POSITIONAL_ONLY, Parameter.POSITIONAL_OR_KEYWORD)
                   for param in s.parameters.values()):
                vars = [SR.var(name) for name in s.parameters.keys()]
                result = x(*vars)
                if isinstance(result, (tuple, list)):
                    return vector(SR, result).function(*vars)
                else:
                    return SR(result).function(*vars)
    return SR(x)
