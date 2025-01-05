# sage.doctest: needs sage.symbolic
r"""
Calculus functions
"""
from sage.misc.lazy_import import lazy_import
from sage.structure.element import Matrix, Vector, Expression

lazy_import('sage.matrix.constructor', 'matrix')

from .functional import diff


def wronskian(*args):
    """
    Return the Wronskian of the provided functions, differentiating with
    respect to the given variable.

    If no variable is provided, diff(f) is called for each function f.

    wronskian(f1,...,fn, x) returns the Wronskian of f1,...,fn, with
    derivatives taken with respect to x.

    wronskian(f1,...,fn) returns the Wronskian of f1,...,fn where
    k'th derivatives are computed by doing ``.derivative(k)`` on each
    function.

    The Wronskian of a list of functions is a determinant of derivatives.
    The `n`-th row (starting from 0) is a list of the `n`-th derivatives of the
    given functions.

    For two functions::

                              | f   g  |
                 W(f, g) = det|        | = f*g' - g*f'.
                              | f'  g' |

    EXAMPLES::

        sage: wronskian(e^x, x^2)
        -x^2*e^x + 2*x*e^x

        sage: x,y = var('x, y')
        sage: wronskian(x*y, log(x), x)
        -y*log(x) + y

    If your functions are in a list, you can use `*' to turn them into
    arguments to :func:`wronskian`::

        sage: wronskian(*[x^k for k in range(1, 5)])
        12*x^4

    If you want to use 'x' as one of the functions in the Wronskian,
    you can't put it last or it will be interpreted as the variable
    with respect to which we differentiate. There are several ways to
    get around this.

    Two-by-two Wronskian of sin(x) and e^x::

        sage: wronskian(sin(x), e^x, x)
        -cos(x)*e^x + e^x*sin(x)

    Or don't put x last::

        sage: wronskian(x, sin(x), e^x)
        (cos(x)*e^x + e^x*sin(x))*x - 2*e^x*sin(x)

    Example where one of the functions is constant::

        sage: wronskian(1, e^(-x), e^(2*x))
        -6*e^x

    REFERENCES:

    - :wikipedia:`Wronskian`
    - http://planetmath.org/encyclopedia/WronskianDeterminant.html

    AUTHORS:

    - Dan Drake (2008-03-12)
    """
    if not args:
        raise TypeError('wronskian() takes at least one argument (0 given)')
    elif len(args) == 1:
        # a 1x1 Wronskian is just its argument
        return args[0]
    else:
        if isinstance(args[-1], Expression) and args[-1].is_symbol():
            # if last argument is a variable, peel it off and
            # differentiate the other args
            v = args[-1]
            fs = args[0:-1]

            def row(n):
                return [diff(f, v, n) for f in fs]
        else:
            # if the last argument is not a variable, just run
            # .derivative on everything
            fs = args

            def row(n):
                return [diff(f, n) for f in fs]
        # NOTE: I rewrote the below as two lines to avoid a possible subtle
        # memory management problem on some platforms (only VMware as far
        # as we know?).  See trac #2990.
        # There may still be a real problem that this is just hiding for now.
        return matrix([row(r) for r in range(len(fs))]).determinant()


def jacobian(functions, variables, hessian=False, sparse=False, rank=False):
    """
    Return the Jacobian matrix, which is the matrix of partial
    derivatives in which the i,j entry of the Jacobian matrix is the
    partial derivative diff(functions[i], variables[j]).

    EXAMPLES::

        sage: x,y = var('x,y')
        sage: g=x^2-2*x*y
        sage: jacobian(g, (x,y))
        [2*x - 2*y      -2*x]

    The Jacobian of the Jacobian should give us the "second derivative", which is the Hessian matrix::

        sage: jacobian(jacobian(g, (x,y)), (x,y))
        [ 2 -2]
        [-2  0]
        sage: g.hessian()
        [ 2 -2]
        [-2  0]

        sage: f=(x^3*sin(y), cos(x)*sin(y), exp(x))
        sage: jacobian(f, (x,y))
        [  3*x^2*sin(y)     x^3*cos(y)]
        [-sin(x)*sin(y)  cos(x)*cos(y)]
        [           e^x              0]
        sage: jacobian(f, (y,x))
        [    x^3*cos(y)   3*x^2*sin(y)]
        [ cos(x)*cos(y) -sin(x)*sin(y)]
        [             0            e^x]
    """
    if isinstance(functions, Matrix) and (functions.nrows() == 1
                                          or functions.ncols() == 1):
        functions = functions.list()
    elif not isinstance(functions, (tuple, list, Vector)):
        functions = [functions]

    if not isinstance(variables, (tuple, list, Vector)):
        variables = [variables]

    jacobian_matrix([[diff(f, v) for v in variables] for f in functions])
    
    if sparse:
        jacobian_matrix = jacobian_matrix.sparse_matrix()
    if hessian:
        if len(functions) > 1:
            raise ValueError("Hessian can only be computed for a single function")
        hessian_matrix = matrix([[diff(jacobian_matrix[0, j], v) for j, v in enumerate(variables)] for v in variables])
        return hessian_matrix
    if rank:
        return jacobian_matrix, jacobian_matrix.rank()
    
    return jacobian_matrix
