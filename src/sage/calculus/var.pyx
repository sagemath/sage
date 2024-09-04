# cython: old_style_globals=True
"""
Symbolic variables
"""

from sage.symbolic.function_factory import function as new_function
from sage.symbolic.ring import SR


def var(*args, **kwds):
    r"""
    Create a symbolic variable with the name *s*.

    INPUT:

    - ``args`` -- a single string ``var('x y')``, a list of strings
      ``var(['x','y'])``, or multiple strings ``var('x', 'y')``. A
      single string can be either a single variable name, or a space
      or comma separated list of variable names. In a list or tuple of
      strings, each entry is one variable. If multiple arguments are
      specified, each argument is taken to be one variable. Spaces
      before or after variable names are ignored.

    - ``kwds`` -- keyword arguments can be given to specify domain and
      custom latex_name for variables; see EXAMPLES for usage

    .. NOTE::

        The new variable is both returned and automatically injected
        into the global namespace. If you need a symbolic variable in
        library code, you must use either ``SR.var()``
        or ``SR.symbol()``.

    OUTPUT:

    If a single symbolic variable was created, the variable
    itself. Otherwise, a tuple of symbolic variables. The variable
    names are checked to be valid Python identifiers and a
    :exc:`ValueError` is raised otherwise.

    EXAMPLES:

    Here are the different ways to define three variables ``x``, ``y``,
    and ``z`` in a single line::

        sage: var('x y z')
        (x, y, z)
        sage: var('x, y, z')
        (x, y, z)
        sage: var(['x', 'y', 'z'])
        (x, y, z)
        sage: var('x', 'y', 'z')
        (x, y, z)
        sage: var('x'), var('y'), var(z)
        (x, y, z)

    We define some symbolic variables::

        sage: var('n xx yy zz')
        (n, xx, yy, zz)

    Then we make an algebraic expression out of them::

        sage: f = xx^n + yy^n + zz^n; f
        xx^n + yy^n + zz^n

    By default, var returns a complex variable. To define real or positive
    variables we can specify the domain as::

        sage: x = var('x', domain=RR); x; x.conjugate()
        x
        x
        sage: y = var('y', domain='real'); y.conjugate()
        y
        sage: y = var('y', domain='positive'); y.abs()
        y

    Custom latex expression can be assigned to variable::

        sage: x = var('sui', latex_name='s_{u,i}'); x._latex_()
        '{s_{u,i}}'

    In notebook, we can also colorize latex expression::

        sage: x = var('sui', latex_name='\\color{red}{s_{u,i}}'); x._latex_()
        '{\\color{red}{s_{u,i}}}'

    We can substitute a new variable name for n::

        sage: f(n = var('sigma'))
        xx^sigma + yy^sigma + zz^sigma

    If you make an important built-in variable into a symbolic variable,
    you can get back the original value using restore::

        sage: var('QQ RR')
        (QQ, RR)
        sage: QQ
        QQ
        sage: restore('QQ')
        sage: QQ
        Rational Field

    We make two new variables separated by commas::

        sage: var('theta, gamma')
        (theta, gamma)
        sage: theta^2 + gamma^3
        gamma^3 + theta^2

    The new variables are of type Expression, and belong
    to the symbolic expression ring::

        sage: type(theta)
        <class 'sage.symbolic.expression.Expression'>
        sage: parent(theta)
        Symbolic Ring
    """
    if len(args) == 1:
        name = args[0]
    else:
        name = args
    G = globals()  # this is the reason the code must be in Cython.
    v = SR.var(name, **kwds)
    if isinstance(v, tuple):
        for x in v:
            G[repr(x)] = x
    else:
        G[repr(v)] = v
    return v


def function(s, **kwds):
    r"""
    Create a formal symbolic function with the name *s*.

    INPUT:

    - ``nargs=0`` -- number of arguments the function accepts, defaults to
      variable number of arguments, or 0
    - ``latex_name`` -- name used when printing in latex mode
    - ``conversions`` -- dictionary specifying names of this function in
      other systems, this is used by the interfaces internally during conversion
    - ``eval_func`` -- method used for automatic evaluation
    - ``evalf_func`` -- method used for numeric evaluation
    - ``evalf_params_first`` -- boolean to indicate if parameters should be
      evaluated numerically before calling the custom evalf function
    - ``conjugate_func`` -- method used for complex conjugation
    - ``real_part_func`` -- method used when taking real parts
    - ``imag_part_func`` -- method used when taking imaginary parts
    - ``derivative_func`` -- method to be used for (partial) derivation
      This method should take a keyword argument deriv_param specifying
      the index of the argument to differentiate w.r.t
    - ``tderivative_func`` -- method to be used for derivatives
    - ``power_func`` -- method used when taking powers
      This method should take a keyword argument power_param specifying
      the exponent
    - ``series_func`` -- method used for series expansion
      This method should expect keyword arguments
      - ``order`` -- order for the expansion to be computed
      - ``var`` -- variable to expand w.r.t.
      - ``at`` -- expand at this value
    - ``print_func`` -- method for custom printing
    - ``print_latex_func`` -- method for custom printing in latex mode

    Note that custom methods must be instance methods, i.e., expect the instance
    of the symbolic function as the first argument.

    .. NOTE::

        The new function is both returned and automatically injected
        into the global namespace.  If you use this function in library
        code, it is better to use sage.symbolic.function_factory.function,
        since it will not touch the global namespace.

    EXAMPLES:

    We create a formal function called supersin ::

        sage: function('supersin')
        supersin

    We can immediately use supersin in symbolic expressions::

        sage: y, z, A = var('y z A')
        sage: supersin(y+z) + A^3
        A^3 + supersin(y + z)

    We can define other functions in terms of supersin::

        sage: g(x,y) = supersin(x)^2 + sin(y/2)
        sage: g
        (x, y) |--> supersin(x)^2 + sin(1/2*y)
        sage: g.diff(y)
        (x, y) |--> 1/2*cos(1/2*y)
        sage: k = g.diff(x); k
        (x, y) |--> 2*supersin(x)*diff(supersin(x), x)

    We create a formal function of one variable, write down
    an expression that involves first and second derivatives,
    and extract off coefficients::

        sage: r, kappa = var('r,kappa')
        sage: psi = function('psi', nargs=1)(r); psi
        psi(r)
        sage: g = 1/r^2*(2*r*psi.derivative(r,1) + r^2*psi.derivative(r,2)); g
        (r^2*diff(psi(r), r, r) + 2*r*diff(psi(r), r))/r^2
        sage: g.expand()
        2*diff(psi(r), r)/r + diff(psi(r), r, r)
        sage: g.coefficient(psi.derivative(r,2))
        1
        sage: g.coefficient(psi.derivative(r,1))
        2/r

    Custom typesetting of symbolic functions in LaTeX, either using latex_name
    keyword::

        sage: function('riemann', latex_name='\\mathcal{R}')
        riemann
        sage: latex(riemann(x))
        \mathcal{R}\left(x\right)

    or passing a custom callable function that returns a latex expression::

        sage: mu,nu = var('mu,nu')
        sage: def my_latex_print(self, *args): return "\\psi_{%s}"%(', '.join(map(latex, args)))
        sage: function('psi', print_latex_func=my_latex_print)
        psi
        sage: latex(psi(mu,nu))
        \psi_{\mu, \nu}

    Defining custom methods for automatic or numeric evaluation, derivation,
    conjugation, etc. is supported::

        sage: def ev(self, x): return 2*x
        sage: foo = function("foo", nargs=1, eval_func=ev)
        sage: foo(x)
        2*x
        sage: foo = function("foo", nargs=1, eval_func=lambda self, x: 5)
        sage: foo(x)
        5
        sage: def ef(self, x): pass
        sage: bar = function("bar", nargs=1, eval_func=ef)
        sage: bar(x)
        bar(x)

        sage: def evalf_f(self, x, parent=None, algorithm=None): return 6
        sage: foo = function("foo", nargs=1, evalf_func=evalf_f)
        sage: foo(x)
        foo(x)
        sage: foo(x).n()
        6

        sage: foo = function("foo", nargs=1, conjugate_func=ev)
        sage: foo(x).conjugate()
        2*x

        sage: def deriv(self, *args,**kwds): print("{} {}".format(args, kwds)); return args[kwds['diff_param']]^2
        sage: foo = function("foo", nargs=2, derivative_func=deriv)
        sage: foo(x,y).derivative(y)
        (x, y) {'diff_param': 1}
        y^2

        sage: def pow(self, x, power_param=None): print("{} {}".format(x, power_param)); return x*power_param
        sage: foo = function("foo", nargs=1, power_func=pow)
        sage: foo(y)^(x+y)
        y x + y
        (x + y)*y

        sage: from pprint import pformat
        sage: def expand(self, *args, **kwds):
        ....:     print("{} {}".format(args, pformat(kwds)))
        ....:     return sum(args[0]^i for i in range(kwds['order']))
        sage: foo = function("foo", nargs=1, series_func=expand)
        sage: foo(y).series(y, 5)
        (y,) {'at': 0, 'options': 0, 'order': 5, 'var': y}
        y^4 + y^3 + y^2 + y + 1

        sage: def my_print(self, *args):
        ....:     return "my args are: " + ', '.join(map(repr, args))
        sage: foo = function('t', nargs=2, print_func=my_print)
        sage: foo(x,y^z)
        my args are: x, y^z

        sage: latex(foo(x,y^z))
        t\left(x, y^{z}\right)
        sage: foo = function('t', nargs=2, print_latex_func=my_print)
        sage: foo(x,y^z)
        t(x, y^z)
        sage: latex(foo(x,y^z))
        my args are: x, y^z
        sage: foo = function('t', nargs=2, latex_name='foo')
        sage: latex(foo(x,y^z))
        foo\left(x, y^{z}\right)

    Chain rule::

        sage: def print_args(self, *args, **kwds): print("args: {}".format(args)); print("kwds: {}".format(kwds)); return args[0]
        sage: foo = function('t', nargs=2, tderivative_func=print_args)
        sage: foo(x,x).derivative(x)
        args: (x, x)
        kwds: {'diff_param': x}
        x
        sage: foo = function('t', nargs=2, derivative_func=print_args)
        sage: foo(x,x).derivative(x)
        args: (x, x)
        kwds: {'diff_param': 0}
        args: (x, x)
        kwds: {'diff_param': 1}
        2*x

    Since Sage 4.0, basic arithmetic with unevaluated functions is no
    longer supported::

        sage: x = var('x')
        sage: f = function('f')
        sage: 2*f
        Traceback (most recent call last):
        ...
        TypeError: unsupported operand parent(s) for *: 'Integer Ring' and '<class 'sage.symbolic.function_factory...NewSymbolicFunction'>'

    You now need to evaluate the function in order to do the arithmetic::

        sage: 2*f(x)
        2*f(x)

    Since Sage 4.0, you need to use :meth:`substitute_function` to
    replace all occurrences of a function with another::

        sage: var('a, b')
        (a, b)
        sage: cr = function('cr')
        sage: f = cr(a)
        sage: g = f.diff(a).integral(b)
        sage: g
        b*diff(cr(a), a)
        sage: g.substitute_function(cr, cos)
        -b*sin(a)

        sage: g.substitute_function(cr, (sin(x) + cos(x)).function(x))
        b*(cos(a) - sin(a))

    TESTS:

    Make sure that :issue:`15860` is fixed and whitespaces are removed::

        sage: function('A, B')
        (A, B)
        sage: B
        B
    """
    G = globals()  # this is the reason the code must be in Cython.
    v = new_function(s, **kwds)
    if isinstance(v, tuple):
        for x in v:
            G[repr(x)] = x
    else:
        G[repr(v)] = v
    return v


def clear_vars():
    """
    Delete all 1-letter symbolic variables that are predefined at
    startup of Sage.

    Any one-letter global variables that are not symbolic variables
    are not cleared.

    EXAMPLES::

        sage: var('x y z')
        (x, y, z)
        sage: (x+y)^z
        (x + y)^z
        sage: k = 15
        sage: clear_vars()
        sage: (x+y)^z
        Traceback (most recent call last):
        ...
        NameError: name 'x' is not defined
        sage: expand((e + i)^2)
        e^2 + 2*I*e - 1
        sage: k
        15
    """
    from sage.structure.element import Expression

    G = globals()
    for i in list(range(65, 65 + 26)) + list(range(97, 97 + 26)):
        chr_i = chr(i)
        if chr_i in G and isinstance(G[chr_i], Expression) and G[chr_i].is_symbol():
            # We check to see if there is a corresponding pyobject
            # associated with the expression.  This will work for
            # constants which we want to keep, but will fail for
            # variables that we want to delete.
            try:
                G[chr(i)].pyobject()
            except TypeError:
                del G[chr(i)]
