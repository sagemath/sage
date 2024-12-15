r"""
Symbolic Computation

AUTHORS:

- Bobby Moretti and William Stein (2006-2007)

- Robert Bradshaw (2007-10): minpoly(), numerical algorithm

- Robert Bradshaw (2008-10): minpoly(), algebraic algorithm

- Golam Mortuza Hossain (2009-06-15): _limit_latex()

- Golam Mortuza Hossain (2009-06-22): _laplace_latex(), _inverse_laplace_latex()

- Tom Coates (2010-06-11): fixed :issue:`9217`

EXAMPLES:

The basic units of the calculus package are symbolic expressions which
are elements of the symbolic expression ring (SR). To create a
symbolic variable object in Sage, use the :func:`var` function, whose
argument is the text of that variable. Note that Sage is intelligent
about LaTeXing variable names.

::

    sage: x1 = var('x1'); x1
    x1
    sage: latex(x1)
    x_{1}
    sage: theta = var('theta'); theta
    theta
    sage: latex(theta)
    \theta

Sage predefines ``x`` to be a global indeterminate.
Thus the following works::

    sage: x^2
    x^2
    sage: type(x)
    <class 'sage.symbolic.expression.Expression'>

More complicated expressions in Sage can be built up using ordinary
arithmetic. The following are valid, and follow the rules of Python
arithmetic: (The '=' operator represents assignment, and not
equality)

::

    sage: var('x,y,z')
    (x, y, z)
    sage: f = x + y + z/(2*sin(y*z/55))
    sage: g = f^f; g
    (x + y + 1/2*z/sin(1/55*y*z))^(x + y + 1/2*z/sin(1/55*y*z))

Differentiation and integration are available, but behind the
scenes through Maxima::

    sage: f = sin(x)/cos(2*y)
    sage: f.derivative(y)
    2*sin(x)*sin(2*y)/cos(2*y)^2
    sage: g = f.integral(x); g
    -cos(x)/cos(2*y)

Note that these methods usually require an explicit variable name. If none
is given, Sage will try to find one for you.

::

    sage: f = sin(x); f.derivative()
    cos(x)

If the expression is a callable symbolic expression (i.e., the
variable order is specified), then Sage can calculate the matrix
derivative (i.e., the gradient, Jacobian matrix, etc.) if no variables
are specified.  In the example below, we use the second derivative
test to determine that there is a saddle point at (0,-1/2).

::

    sage: f(x,y) = x^2*y + y^2 + y
    sage: f.diff()  # gradient
    (x, y) |--> (2*x*y, x^2 + 2*y + 1)
    sage: solve(list(f.diff()), [x,y])
    [[x == -I, y == 0], [x == I, y == 0], [x == 0, y == (-1/2)]]
    sage: H=f.diff(2); H  # Hessian matrix
    [(x, y) |--> 2*y (x, y) |--> 2*x]
    [(x, y) |--> 2*x   (x, y) |--> 2]
    sage: H(x=0, y=-1/2)
    [-1  0]
    [ 0  2]
    sage: H(x=0, y=-1/2).eigenvalues()
    [-1, 2]

Here we calculate the Jacobian for the polar coordinate transformation::

    sage: T(r,theta) = [r*cos(theta),r*sin(theta)]
    sage: T
    (r, theta) |--> (r*cos(theta), r*sin(theta))
    sage: T.diff() # Jacobian matrix
    [   (r, theta) |--> cos(theta) (r, theta) |--> -r*sin(theta)]
    [   (r, theta) |--> sin(theta)  (r, theta) |--> r*cos(theta)]
    sage: diff(T) # Jacobian matrix
    [   (r, theta) |--> cos(theta) (r, theta) |--> -r*sin(theta)]
    [   (r, theta) |--> sin(theta)  (r, theta) |--> r*cos(theta)]
    sage: T.diff().det() # Jacobian
    (r, theta) |--> r*cos(theta)^2 + r*sin(theta)^2

When the order of variables is ambiguous, Sage will raise an
exception when differentiating::

    sage: f = sin(x+y); f.derivative()
    Traceback (most recent call last):
    ...
    ValueError: No differentiation variable specified.

Simplifying symbolic sums is also possible, using the
:func:`sum` command, which also uses Maxima in the background::

    sage: k, m = var('k, m')
    sage: sum(1/k^4, k, 1, oo)
    1/90*pi^4
    sage: sum(binomial(m,k), k, 0, m)
    2^m

Symbolic matrices can be used as well in various ways,
including exponentiation::

    sage: M = matrix([[x,x^2],[1/x,x]])
    sage: M^2
    [x^2 + x   2*x^3]
    [      2 x^2 + x]
    sage: e^M
    [          1/2*(e^(2*sqrt(x)) + 1)*e^(x - sqrt(x)) 1/2*(x*e^(2*sqrt(x)) - x)*sqrt(x)*e^(x - sqrt(x))]
    [  1/2*(e^(2*sqrt(x)) - 1)*e^(x - sqrt(x))/x^(3/2)           1/2*(e^(2*sqrt(x)) + 1)*e^(x - sqrt(x))]

Complex exponentiation works, but may require a patched version of
maxima (:issue:`32898`) for now::

    sage: M = i*matrix([[pi]])
    sage: e^M  # not tested, requires patched maxima
    [-1]
    sage: M = i*matrix([[pi,0],[0,2*pi]])
    sage: e^M
    [-1  0]
    [ 0  1]
    sage: M = matrix([[0,pi],[-pi,0]])
    sage: e^M
    [-1  0]
    [ 0 -1]

Substitution works similarly. We can substitute with a python
dict::

    sage: f = sin(x*y - z)
    sage: f({x: var('t'), y: z})
    sin(t*z - z)

Also we can substitute with keywords::

    sage: f = sin(x*y - z)
    sage: f(x=t, y=z)
    sin(t*z - z)

Another example::

    sage: f = sin(2*pi*x/y)
    sage: f(x=4)
    sin(8*pi/y)

It is no longer allowed to call expressions with positional arguments::

    sage: f = sin(x)
    sage: f(y)
    Traceback (most recent call last):
    ...
    TypeError: Substitution using function-call syntax and unnamed
    arguments has been removed. You can use named arguments instead, like
    EXPR(x=..., y=...)
    sage: f(x=pi)
    0

We can also make a :class:`CallableSymbolicExpression`,
which is a :class:`SymbolicExpression` that is a function of
specified variables in a fixed order. Each
:class:`SymbolicExpression` has a
``function(...)`` method that is used to create a
:class:`CallableSymbolicExpression`, as illustrated below::

    sage: u = log((2-x)/(y+5))
    sage: f = u.function(x, y); f
    (x, y) |--> log(-(x - 2)/(y + 5))

There is an easier way of creating a
:class:`CallableSymbolicExpression`, which relies on the
Sage preparser.

::

    sage: f(x,y) = log(x)*cos(y); f
    (x, y) |--> cos(y)*log(x)

Then we have fixed an order of variables and there is no ambiguity
substituting or evaluating::

    sage: f(x,y) = log((2-x)/(y+5))
    sage: f(7,t)
    log(-5/(t + 5))

Some further examples::

    sage: f = 5*sin(x)
    sage: f
    5*sin(x)
    sage: f(x=2)
    5*sin(2)
    sage: f(x=pi)
    0
    sage: float(f(x=pi))
    0.0

Another example::

    sage: f = integrate(1/sqrt(9+x^2), x); f
    arcsinh(1/3*x)
    sage: f(x=3)
    arcsinh(1)
    sage: f.derivative(x)
    1/sqrt(x^2 + 9)

We compute the length of the parabola from 0 to 2::

    sage: x = var('x')
    sage: y = x^2
    sage: dy = derivative(y,x)
    sage: z = integral(sqrt(1 + dy^2), x, 0, 2)
    sage: z
    sqrt(17) + 1/4*arcsinh(4)
    sage: n(z,200)
    4.6467837624329358733826155674904591885104869874232887508703
    sage: float(z)
    4.646783762432936

We test pickling::

    sage: x, y = var('x,y')
    sage: f = -sqrt(pi)*(x^3 + sin(x/cos(y)))
    sage: bool(loads(dumps(f)) == f)
    True

Coercion examples:

We coerce various symbolic expressions into the complex numbers::

    sage: CC(I)
    1.00000000000000*I
    sage: CC(2*I)
    2.00000000000000*I
    sage: ComplexField(200)(2*I)
    2.0000000000000000000000000000000000000000000000000000000000*I
    sage: ComplexField(200)(sin(I))
    1.1752011936438014568823818505956008151557179813340958702296*I
    sage: f = sin(I) + cos(I/2); f
    cosh(1/2) + I*sinh(1)
    sage: CC(f)
    1.12762596520638 + 1.17520119364380*I
    sage: ComplexField(200)(f)
    1.1276259652063807852262251614026720125478471180986674836290
     + 1.1752011936438014568823818505956008151557179813340958702296*I
    sage: ComplexField(100)(f)
    1.1276259652063807852262251614 + 1.1752011936438014568823818506*I

We illustrate construction of an inverse sum where each denominator
has a new variable name::

    sage: f = sum(1/var('n%s'%i)^i for i in range(10))
    sage: f
    1/n1 + 1/n2^2 + 1/n3^3 + 1/n4^4 + 1/n5^5 + 1/n6^6 + 1/n7^7 + 1/n8^8 + 1/n9^9 + 1

Note that after calling var, the variables are immediately
available for use::

    sage: (n1 + n2)^5
    (n1 + n2)^5

We can, of course, substitute::

    sage: f(n9=9, n7=n6)
    1/n1 + 1/n2^2 + 1/n3^3 + 1/n4^4 + 1/n5^5 + 1/n6^6 + 1/n6^7 + 1/n8^8
     + 387420490/387420489

TESTS:

Substitution::

    sage: f = x
    sage: f(x=5)
    5

Simplifying expressions involving scientific notation::

    sage: k = var('k')
    sage: a0 = 2e-06; a1 = 12
    sage: c = a1 + a0*k; c
    (2.00000000000000e-6)*k + 12
    sage: sqrt(c)
    sqrt((2.00000000000000e-6)*k + 12)
    sage: sqrt(c^3)
    sqrt(((2.00000000000000e-6)*k + 12)^3)

The symbolic calculus package uses its own copy of Maxima for
simplification, etc., which is separate from the default
system-wide version::

    sage: maxima.eval('[x,y]: [1,2]')
    '[1,2]'
    sage: maxima.eval('expand((x+y)^3)')
    '27'

If the copy of Maxima used by the symbolic calculus package were
the same as the default one, then the following would return 27,
which would be very confusing indeed!

::

    sage: x, y = var('x,y')
    sage: expand((x+y)^3)
    x^3 + 3*x^2*y + 3*x*y^2 + y^3

Set x to be 5 in maxima::

    sage: maxima('x: 5')
    5
    sage: maxima('x + x + %pi')
    %pi+10

Simplifications like these are now done using Pynac::

    sage: x + x + pi
    pi + 2*x

But this still uses Maxima::

    sage: (x + x + pi).simplify()
    pi + 2*x

Note that ``x`` is still ``x``, since the
maxima used by the calculus package is different than the one in
the interactive interpreter.

Check to see that the problem with the variables method mentioned
in :issue:`3779` is actually fixed::

    sage: f = function('F')(x)
    sage: diff(f*SR(1),x)
    diff(F(x), x)

Doubly ensure that :issue:`7479` is working::

    sage: f(x)=x
    sage: integrate(f,x,0,1)
    1/2

Check that the problem with Taylor expansions of the gamma function
(:issue:`9217`) is fixed::

    sage: taylor(gamma(1/3+x),x,0,3)
    -1/432*((72*euler_gamma^3 + 36*euler_gamma^2*(sqrt(3)*pi + 9*log(3)) + ...
    sage: [f[0].n() for f in _.coefficients()]  # numerical coefficients to make comparison easier; Maple 12 gives same answer
    [2.6789385347..., -8.3905259853..., 26.662447494..., -80.683148377...]

Ensure that :issue:`8582` is fixed::

    sage: k = var("k")
    sage: sum(1/(1+k^2), k, -oo, oo)
    -1/2*I*psi(I + 1) + 1/2*I*psi(-I + 1) - 1/2*I*psi(I) + 1/2*I*psi(-I)

Ensure that :issue:`8624` is fixed::

    sage: integrate(abs(cos(x)) * sin(x), x, pi/2, pi)
    1/2
    sage: integrate(sqrt(cos(x)^2 + sin(x)^2), x, 0, 2*pi)
    2*pi

Ensure that :issue:`25626` is fixed. As the form of the answer is dependent of
the giac version, we simplify it (see :issue:`34037`)::

    sage: # needs sage.libs.giac
    sage: t = SR.var('t')
    sage: integrate(exp(t)/(t + 1)^2, t, algorithm='giac').full_simplify()
    ((t + 1)*Ei(t + 1) - e^(t + 1))/(t*e + e)

Check if maxima has redundant variables defined after initialization,
see :issue:`9538`::

    sage: maxima = sage.interfaces.maxima.maxima
    sage: maxima('f1')
    f1
    sage: sage.calculus.calculus.maxima('f1')
    f1

To check that :issue:`14821` is fixed::

    sage: H = exp(-1.0 * x)
    sage: H.integral(x, 0, 1)
    0.6321205588285577
    sage: result = integral(exp(-300.0/(-0.064*x+14.0)),x,0.0,120.0)
    ...
    sage: result  # abs tol 1e-10
    4.62770039817000e-9

To check that :issue:`27092` is fixed::

    sage: n = var('n')
    sage: sum(binomial(1, n), n, 0, oo)
    2
"""

import math
import re
from sage.arith.misc import algdep
from sage.rings.complex_arb import ComplexBallField
from sage.rings.complex_interval_field import ComplexIntervalField_class
from sage.rings.complex_mpfr import ComplexField_class
from sage.rings.integer import Integer
from sage.rings.rational_field import QQ
from sage.rings.real_arb import RealBallField
from sage.rings.real_mpfi import RealIntervalField_class
from sage.rings.real_double import RealDoubleElement
from sage.rings.real_mpfr import RR, create_RealNumber, RealField_class
from sage.rings.cc import CC
from sage.structure.element import parent

from sage.misc.latex import latex
from sage.misc.parser import Parser, LookupNameMaker
from sage.structure.element import Expression
from sage.symbolic.ring import var, SR
from sage.symbolic.symbols import symbol_table
from sage.symbolic.function import Function
from sage.symbolic.function_factory import function_factory
from sage.symbolic.integration.integral import (indefinite_integral,
        definite_integral)

from sage.misc.lazy_import import lazy_import
lazy_import('sage.interfaces.maxima_lib', 'maxima')
from types import FunctionType


########################################################
def symbolic_sum(expression, v, a, b, algorithm='maxima', hold=False):
    r"""
    Return the symbolic sum `\sum_{v = a}^b expression` with respect
    to the variable `v` with endpoints `a` and `b`.

    INPUT:

    - ``expression`` -- a symbolic expression

    - ``v`` -- a variable or variable name

    - ``a`` -- lower endpoint of the sum

    - ``b`` -- upper endpoint of the sum

    - ``algorithm`` -- (default: ``'maxima'``)  one of

      - ``'maxima'`` -- use Maxima (the default)

      - ``'maple'`` -- (optional) use Maple

      - ``'mathematica'`` -- (optional) use Mathematica

      - ``'giac'`` -- (optional) use Giac

      - ``'sympy'`` -- use SymPy

    - ``hold`` -- boolean (default: ``False``); if ``True``, don't evaluate

    EXAMPLES::

        sage: k, n = var('k,n')
        sage: from sage.calculus.calculus import symbolic_sum
        sage: symbolic_sum(k, k, 1, n).factor()
        1/2*(n + 1)*n

    ::

        sage: symbolic_sum(1/k^4, k, 1, oo)
        1/90*pi^4

    ::

        sage: symbolic_sum(1/k^5, k, 1, oo)
        zeta(5)

    A well known binomial identity::

        sage: symbolic_sum(binomial(n,k), k, 0, n)
        2^n

    And some truncations thereof::

        sage: assume(n>1)
        sage: symbolic_sum(binomial(n,k), k, 1, n)
        2^n - 1
        sage: symbolic_sum(binomial(n,k), k, 2, n)
        2^n - n - 1
        sage: symbolic_sum(binomial(n,k), k, 0, n-1)
        2^n - 1
        sage: symbolic_sum(binomial(n,k), k, 1, n-1)
        2^n - 2

    The binomial theorem::

        sage: x, y = var('x, y')
        sage: symbolic_sum(binomial(n,k) * x^k * y^(n-k), k, 0, n)
        (x + y)^n

    ::

        sage: symbolic_sum(k * binomial(n, k), k, 1, n)
        2^(n - 1)*n

    ::

        sage: symbolic_sum((-1)^k*binomial(n,k), k, 0, n)
        0

    ::

        sage: symbolic_sum(2^(-k)/(k*(k+1)), k, 1, oo)
        -log(2) + 1

    Summing a hypergeometric term::

        sage: symbolic_sum(binomial(n, k) * factorial(k) / factorial(n+1+k), k, 0, n)
        1/2*sqrt(pi)/factorial(n + 1/2)

    We check a well known identity::

        sage: bool(symbolic_sum(k^3, k, 1, n) == symbolic_sum(k, k, 1, n)^2)
        True

    A geometric sum::

        sage: a, q = var('a, q')
        sage: symbolic_sum(a*q^k, k, 0, n)
        (a*q^(n + 1) - a)/(q - 1)

    For the geometric series, we will have to assume
    the right values for the sum to converge::

        sage: assume(abs(q) < 1)
        sage: symbolic_sum(a*q^k, k, 0, oo)
        -a/(q - 1)

    A divergent geometric series.  Don't forget
    to forget your assumptions::

        sage: forget()
        sage: assume(q > 1)
        sage: symbolic_sum(a*q^k, k, 0, oo)
        Traceback (most recent call last):
        ...
        ValueError: Sum is divergent.
        sage: forget()
        sage: assumptions()  # check the assumptions were really forgotten
        []

    A summation performed by Mathematica::

        sage: symbolic_sum(1/(1+k^2), k, -oo, oo, algorithm='mathematica')     # optional - mathematica
        pi*coth(pi)

    An example of this summation with Giac::

        sage: # needs giac
        sage: symbolic_sum(1/(1+k^2), k, -oo, oo, algorithm='giac').factor()
        pi*(e^(2*pi) + 1)/((e^pi + 1)*(e^pi - 1))

    The same summation is solved by SymPy::

        sage: symbolic_sum(1/(1+k^2), k, -oo, oo, algorithm='sympy')
        pi/tanh(pi)

    SymPy and Maxima 5.39.0 can do the following (see
    :issue:`22005`)::

        sage: sum(1/((2*n+1)^2-4)^2, n, 0, Infinity, algorithm='sympy')
        1/64*pi^2
        sage: sum(1/((2*n+1)^2-4)^2, n, 0, Infinity)
        1/64*pi^2

    Use Maple as a backend for summation::

        sage: symbolic_sum(binomial(n,k)*x^k, k, 0, n, algorithm='maple')      # optional - maple
        (x + 1)^n

    If you don't want to evaluate immediately give the ``hold`` keyword::

        sage: s = sum(n, n, 1, k, hold=True); s
        sum(n, n, 1, k)
        sage: s.unhold()
        1/2*k^2 + 1/2*k
        sage: s.subs(k == 10)
        sum(n, n, 1, 10)
        sage: s.subs(k == 10).unhold()
        55
        sage: s.subs(k == 10).n()
        55.0000000000000

    TESTS:

    :issue:`10564` is fixed::

        sage: sum (n^3 * x^n, n, 0, infinity)
        (x^3 + 4*x^2 + x)/(x^4 - 4*x^3 + 6*x^2 - 4*x + 1)

    .. NOTE::

       Sage can currently only understand a subset of the output of Maxima,
       Maple and Mathematica, so even if the chosen backend can perform
       the summation the result might not be convertible into a Sage
       expression.
    """
    if not (isinstance(v, Expression) and v.is_symbol()):
        if isinstance(v, str):
            v = var(v)
        else:
            raise TypeError("need a summation variable")

    if v in SR(a).variables() or v in SR(b).variables():
        raise ValueError("summation limits must not depend on the summation variable")

    if hold:
        from sage.functions.other import symbolic_sum as ssum
        return ssum(expression, v, a, b)

    if algorithm == 'maxima':
        return maxima.sr_sum(expression,v,a,b)

    elif algorithm == 'mathematica':
        try:
            sum = "Sum[%s, {%s, %s, %s}]" % tuple([repr(expr._mathematica_()) for expr in (expression, v, a, b)])
        except TypeError:
            raise ValueError("Mathematica cannot make sense of input")
        from sage.interfaces.mathematica import mathematica
        try:
            result = mathematica(sum)
        except TypeError:
            raise ValueError("Mathematica cannot make sense of: %s" % sum)
        return result.sage()

    elif algorithm == 'maple':
        sum = "sum(%s, %s=%s..%s)" % tuple([repr(expr._maple_()) for expr in (expression, v, a, b)])
        from sage.interfaces.maple import maple
        try:
            result = maple(sum).simplify()
        except TypeError:
            raise ValueError("Maple cannot make sense of: %s" % sum)
        return result.sage()

    elif algorithm == 'giac':
        sum = "sum(%s, %s, %s, %s)" % tuple([repr(expr._giac_()) for expr in (expression, v, a, b)])
        from sage.interfaces.giac import giac
        try:
            result = giac(sum)
        except TypeError:
            raise ValueError("Giac cannot make sense of: %s" % sum)
        return result.sage()

    elif algorithm == 'sympy':
        expression,v,a,b = (expr._sympy_() for expr in (expression, v, a, b))
        from sympy import summation
        from sage.interfaces.sympy import sympy_init
        sympy_init()
        result = summation(expression, (v, a, b))
        try:
            return result._sage_()
        except AttributeError:
            raise AttributeError("Unable to convert SymPy result (={}) into"
                    " Sage".format(result))

    else:
        raise ValueError("unknown algorithm: %s" % algorithm)


def nintegral(ex, x, a, b,
              desired_relative_error='1e-8',
              maximum_num_subintervals=200):
    r"""
    Return a floating point machine precision numerical approximation
    to the integral of ``self`` from `a` to
    `b`, computed using floating point arithmetic via maxima.

    INPUT:

    - ``x`` -- variable to integrate with respect to

    - ``a`` -- lower endpoint of integration

    - ``b`` -- upper endpoint of integration

    - ``desired_relative_error`` -- (default: ``1e-8``) the
      desired relative error

    - ``maximum_num_subintervals`` -- (default: 200)
      maximal number of subintervals

    OUTPUT: float; approximation to the integral

    - float: estimated absolute error of the
      approximation

    - the number of integrand evaluations

    - an error code:

      - ``0`` -- no problems were encountered

      - ``1`` -- too many subintervals were done

      - ``2`` -- excessive roundoff error

      - ``3`` -- extremely bad integrand behavior

      - ``4`` -- failed to converge

      - ``5`` -- integral is probably divergent or slowly
        convergent

      - ``6`` -- the input is invalid; this includes the case of
        ``desired_relative_error`` being too small to be achieved

    ALIAS: :func:`nintegrate` is the same as :func:`nintegral`

    REMARK: There is also a function
    :func:`numerical_integral` that implements numerical
    integration using the GSL C library. It is potentially much faster
    and applies to arbitrary user defined functions.

    Also, there are limits to the precision to which Maxima can compute
    the integral due to limitations in quadpack.
    In the following example, remark that the last value of the returned
    tuple is ``6``, indicating that the input was invalid, in this case
    because of a too high desired precision.

    ::

        sage: f = x
        sage: f.nintegral(x, 0, 1, 1e-14)
        (0.0, 0.0, 0, 6)

    EXAMPLES::

        sage: f(x) = exp(-sqrt(x))
        sage: f.nintegral(x, 0, 1)
        (0.5284822353142306, 4.163...e-11, 231, 0)

    We can also use the :func:`numerical_integral` function,
    which calls the GSL C library.

    ::

        sage: numerical_integral(f, 0, 1)
        (0.528482232253147, 6.83928460...e-07)

    Note that in exotic cases where floating point evaluation of the
    expression leads to the wrong value, then the output can be
    completely wrong::

        sage: f = exp(pi*sqrt(163)) - 262537412640768744

    Despite appearance, `f` is really very close to 0, but one
    gets a nonzero value since the definition of
    ``float(f)`` is that it makes all constants inside the
    expression floats, then evaluates each function and each arithmetic
    operation using float arithmetic::

        sage: float(f)
        -480.0

    Computing to higher precision we see the truth::

        sage: f.n(200)
        -7.4992740280181431112064614366622348652078895136533593355718e-13
        sage: f.n(300)
        -7.49927402801814311120646143662663009137292462589621789352095066181709095575681963967103004e-13

    Now numerically integrating, we see why the answer is wrong::

        sage: f.nintegrate(x,0,1)
        (-480.000000000000..., 5.32907051820075...e-12, 21, 0)

    It is just because every floating point evaluation of `f` returns `-480.0`
    in floating point.

    Important note: using PARI/GP one can compute numerical integrals
    to high precision::

        sage: gp.eval('intnum(x=17,42,exp(-x^2)*log(x))')
        '2.565728500561051474934096410 E-127'            # 32-bit
        '2.5657285005610514829176211363206621657 E-127'  # 64-bit
        sage: old_prec = gp.set_real_precision(50)
        sage: gp.eval('intnum(x=17,42,exp(-x^2)*log(x))')
        '2.5657285005610514829173563961304957417746108003917 E-127'
        sage: gp.set_real_precision(old_prec)
        57

    Note that the input function above is a string in PARI syntax.
    """
    try:
        v = ex._maxima_().quad_qags(x, a, b,
                                    epsrel=desired_relative_error,
                                    limit=maximum_num_subintervals)
    except TypeError as err:
        if "ERROR" in str(err):
            raise ValueError("Maxima (via quadpack) cannot compute the integral")
        else:
            raise TypeError(err)

    # Maxima returns unevaluated expressions when the underlying library fails
    # to perform numerical integration. See:
    # http://www.math.utexas.edu/pipermail/maxima/2008/012975.html
    if 'quad_qags' in str(v):
        raise ValueError("Maxima (via quadpack) cannot compute the integral")

    return float(v[0]), float(v[1]), Integer(v[2]), Integer(v[3])


nintegrate = nintegral


def symbolic_product(expression, v, a, b, algorithm='maxima', hold=False):
    r"""
    Return the symbolic product `\prod_{v = a}^b expression` with respect
    to the variable `v` with endpoints `a` and `b`.

    INPUT:

    - ``expression`` -- a symbolic expression

    - ``v`` -- a variable or variable name

    - ``a`` -- lower endpoint of the product

    - ``b`` -- upper endpoint of the prduct

    - ``algorithm`` -- (default: ``'maxima'``)  one of

      - ``'maxima'`` -- use Maxima (the default)

      - ``'giac'`` -- use Giac (optional)

      - ``'sympy'`` -- use SymPy

      - ``'mathematica'`` -- (optional) use Mathematica

    - ``hold`` -- boolean (default: ``False``); if ``True``, don't evaluate

    EXAMPLES::

        sage: i, k, n = var('i,k,n')
        sage: from sage.calculus.calculus import symbolic_product
        sage: symbolic_product(k, k, 1, n)
        factorial(n)
        sage: symbolic_product(x + i*(i+1)/2, i, 1, 4)
        x^4 + 20*x^3 + 127*x^2 + 288*x + 180
        sage: symbolic_product(i^2, i, 1, 7)
        25401600
        sage: f = function('f')
        sage: symbolic_product(f(i), i, 1, 7)
        f(7)*f(6)*f(5)*f(4)*f(3)*f(2)*f(1)
        sage: symbolic_product(f(i), i, 1, n)
        product(f(i), i, 1, n)
        sage: assume(k>0)
        sage: symbolic_product(integrate (x^k, x, 0, 1), k, 1, n)
        1/factorial(n + 1)
        sage: symbolic_product(f(i), i, 1, n).log().log_expand()
        sum(log(f(i)), i, 1, n)

    TESTS:

    Verify that :issue:`30520` is fixed::

        sage: symbolic_product(-x^2,x,1,n)
        (-1)^n*factorial(n)^2
    """
    if not (isinstance(v, Expression) and v.is_symbol()):
        if isinstance(v, str):
            v = var(v)
        else:
            raise TypeError("need a multiplication variable")

    if v in SR(a).variables() or v in SR(b).variables():
        raise ValueError("product limits must not depend on the multiplication variable")

    if hold:
        from sage.functions.other import symbolic_product as sprod
        return sprod(expression, v, a, b)

    if algorithm == 'maxima':
        return maxima.sr_prod(expression,v,a,b)

    elif algorithm == 'mathematica':
        try:
            prod = "Product[%s, {%s, %s, %s}]" % tuple([repr(expr._mathematica_()) for expr in (expression, v, a, b)])
        except TypeError:
            raise ValueError("Mathematica cannot make sense of input")
        from sage.interfaces.mathematica import mathematica
        try:
            result = mathematica(prod)
        except TypeError:
            raise ValueError("Mathematica cannot make sense of: %s" % sum)
        return result.sage()

    elif algorithm == 'giac':
        prod = "product(%s, %s, %s, %s)" % tuple([repr(expr._giac_()) for expr in (expression, v, a, b)])
        from sage.interfaces.giac import giac
        try:
            result = giac(prod)
        except TypeError:
            raise ValueError("Giac cannot make sense of: %s" % sum)
        return result.sage()

    elif algorithm == 'sympy':
        expression,v,a,b = (expr._sympy_() for expr in (expression, v, a, b))
        from sympy import product as sproduct
        from sage.interfaces.sympy import sympy_init
        sympy_init()
        result = sproduct(expression, (v, a, b))
        try:
            return result._sage_()
        except AttributeError:
            raise AttributeError("Unable to convert SymPy result (={}) into"
                    " Sage".format(result))

    else:
        raise ValueError("unknown algorithm: %s" % algorithm)


def minpoly(ex, var='x', algorithm=None, bits=None, degree=None, epsilon=0):
    r"""
    Return the minimal polynomial of ``self``, if possible.

    INPUT:

    - ``var`` -- polynomial variable name (default: ``'x'``)

    - ``algorithm`` -- ``'algebraic'`` or ``'numerical'`` (default
      both, but with numerical first)

    - ``bits`` -- the number of bits to use in numerical
      approx

    - ``degree`` -- the expected algebraic degree

    - ``epsilon`` -- return without error as long as
      f(self) epsilon, in the case that the result cannot be proven

      All of the above parameters are optional, with epsilon=0, ``bits`` and
      ``degree`` tested up to 1000 and 24 by default respectively. The
      numerical algorithm will be faster if bits and/or degree are given
      explicitly. The algebraic algorithm ignores the last three
      parameters.

    OUTPUT: the minimal polynomial of ``self``. If the numerical algorithm
    is used, then it is proved symbolically when ``epsilon=0`` (default).

    If the minimal polynomial could not be found, two distinct kinds of
    errors are raised. If no reasonable candidate was found with the
    given ``bits``/``degree`` parameters, a :exc:`ValueError` will be
    raised. If a reasonable candidate was found but (perhaps due to
    limits in the underlying symbolic package) was unable to be proved
    correct, a :exc:`NotImplementedError` will be raised.

    ALGORITHM: Two distinct algorithms are used, depending on the
    algorithm parameter. By default, the numerical algorithm is
    attempted first, then the algebraic one.

    Algebraic: Attempt to evaluate this expression in ``QQbar``, using
    cyclotomic fields to resolve exponential and trig functions at
    rational multiples of `\pi`, field extensions to handle roots and
    rational exponents, and computing compositums to represent the full
    expression as an element of a number field where the minimal
    polynomial can be computed exactly. The ``bits``, ``degree``, and ``epsilon``
    parameters are ignored.

    Numerical: Computes a numerical approximation of
    ``self`` and use PARI's :pari:`algdep` to get a candidate
    minpoly `f`. If `f(\mathtt{self})`,
    evaluated to a higher precision, is close enough to 0 then evaluate
    `f(\mathtt{self})` symbolically, attempting to prove
    vanishing. If this fails, and ``epsilon`` is nonzero,
    return `f` if and only if
    `f(\mathtt{self}) < \mathtt{epsilon}`.
    Otherwise raise a :exc:`ValueError` (if no suitable
    candidate was found) or a :exc:`NotImplementedError` (if a
    likely candidate was found but could not be proved correct).

    EXAMPLES: First some simple examples::

        sage: sqrt(2).minpoly()
        x^2 - 2
        sage: minpoly(2^(1/3))
        x^3 - 2
        sage: minpoly(sqrt(2) + sqrt(-1))
        x^4 - 2*x^2 + 9
        sage: minpoly(sqrt(2)-3^(1/3))
        x^6 - 6*x^4 + 6*x^3 + 12*x^2 + 36*x + 1


    Works with trig and exponential functions too.

    ::

        sage: sin(pi/3).minpoly()
        x^2 - 3/4
        sage: sin(pi/7).minpoly()
        x^6 - 7/4*x^4 + 7/8*x^2 - 7/64
        sage: minpoly(exp(I*pi/17))
        x^16 - x^15 + x^14 - x^13 + x^12 - x^11 + x^10 - x^9 + x^8
         - x^7 + x^6 - x^5 + x^4 - x^3 + x^2 - x + 1

    Here we verify it gives the same result as the abstract number
    field.

    ::

        sage: (sqrt(2) + sqrt(3) + sqrt(6)).minpoly()
        x^4 - 22*x^2 - 48*x - 23
        sage: K.<a,b> = NumberField([x^2-2, x^2-3])
        sage: (a+b+a*b).absolute_minpoly()
        x^4 - 22*x^2 - 48*x - 23

    The :func:`minpoly` function is used implicitly when creating
    number fields::

        sage: x = var('x')
        sage: eqn =  x^3 + sqrt(2)*x + 5 == 0
        sage: a = solve(eqn, x)[0].rhs()
        sage: QQ[a]
        Number Field in a with defining polynomial x^6 + 10*x^3 - 2*x^2 + 25
         with a = 0.7185272465828846? - 1.721353471724806?*I

    Here we solve a cubic and then recover it from its complicated
    radical expansion.

    ::

        sage: f = x^3 - x + 1
        sage: a = f.solve(x)[0].rhs(); a
        -1/2*(1/18*sqrt(23)*sqrt(3) - 1/2)^(1/3)*(I*sqrt(3) + 1)
         - 1/6*(-I*sqrt(3) + 1)/(1/18*sqrt(23)*sqrt(3) - 1/2)^(1/3)
        sage: a.minpoly()
        x^3 - x + 1

    Note that simplification may be necessary to see that the minimal
    polynomial is correct.

    ::

        sage: a = sqrt(2)+sqrt(3)+sqrt(5)
        sage: f = a.minpoly(); f
        x^8 - 40*x^6 + 352*x^4 - 960*x^2 + 576
        sage: f(a)
        (sqrt(5) + sqrt(3) + sqrt(2))^8 - 40*(sqrt(5) + sqrt(3) + sqrt(2))^6
         + 352*(sqrt(5) + sqrt(3) + sqrt(2))^4 - 960*(sqrt(5) + sqrt(3) + sqrt(2))^2
         + 576
        sage: f(a).expand()
        0

    ::

        sage: a = sin(pi/7)
        sage: f = a.minpoly(algorithm='numerical'); f
        x^6 - 7/4*x^4 + 7/8*x^2 - 7/64
        sage: f(a).horner(a).numerical_approx(100)
        0.00000000000000000000000000000

    The degree must be high enough (default tops out at 24).

    ::

        sage: a = sqrt(3) + sqrt(2)
        sage: a.minpoly(algorithm='numerical', bits=100, degree=3)
        Traceback (most recent call last):
        ...
        ValueError: Could not find minimal polynomial (100 bits, degree 3).
        sage: a.minpoly(algorithm='numerical', bits=100, degree=10)
        x^4 - 10*x^2 + 1

    ::

        sage: cos(pi/33).minpoly(algorithm='algebraic')
        x^10 + 1/2*x^9 - 5/2*x^8 - 5/4*x^7 + 17/8*x^6 + 17/16*x^5
         - 43/64*x^4 - 43/128*x^3 + 3/64*x^2 + 3/128*x + 1/1024
        sage: cos(pi/33).minpoly(algorithm='numerical')
        x^10 + 1/2*x^9 - 5/2*x^8 - 5/4*x^7 + 17/8*x^6 + 17/16*x^5
         - 43/64*x^4 - 43/128*x^3 + 3/64*x^2 + 3/128*x + 1/1024

    Sometimes it fails, as it must given that some numbers aren't algebraic::

        sage: sin(1).minpoly(algorithm='numerical')
        Traceback (most recent call last):
        ...
        ValueError: Could not find minimal polynomial (1000 bits, degree 24).

    .. NOTE::

       Of course, failure to produce a minimal polynomial does not
       necessarily indicate that this number is transcendental.
    """
    if isinstance(ex, Expression):
        try:
            ex = ex.pyobject()
        except TypeError:
            pass

    if algorithm is None or algorithm.startswith('numeric'):
        def estimated_bit_length(f):
            return sum(x.numer().bit_length() + x.denom().bit_length() + 1 for x in f)

        if parent(ex) is float or isinstance(parent(ex),
                                             (RealBallField, ComplexBallField, RealIntervalField_class,
                                              ComplexIntervalField_class, RealField_class, ComplexField_class)):
            best_f = None
            cur_degree = degree or 2  # ``algdep(CC(I), 1)`` just return ``x``
            prec = parent(ex).prec()
            while True:
                if best_f is not None and cur_degree > estimated_bit_length(best_f) / 4:
                    break
                if cur_degree > prec / 6:
                    # experimentally algdep(RealField(500)(2^(1/100)), 100) returns the wrong answer
                    break
                try:
                    f = algdep(ex, cur_degree)
                except ValueError:
                    pass
                else:
                    f = QQ[var](f)
                    if best_f is None or estimated_bit_length(f) < estimated_bit_length(best_f):
                        best_f = f
                    if degree:
                        # user-specified, do not try any other degree value
                        break
                cur_degree = math.ceil(cur_degree*1.2)

            if best_f is not None:
                return best_f
            raise ValueError("Could not find minimal polynomial")

        bits_list = [bits] if bits else [100,200,500,1000]
        degree_list = [degree] if degree else [2,4,8,12,24]

        for bits in bits_list:
            a = ex.numerical_approx(bits)
            check_bits = int(1.25 * bits + 80)
            aa = ex.numerical_approx(check_bits)

            for degree in degree_list:

                f = QQ[var](algdep(a, degree))  # TODO: use the known_bits parameter?
                # If indeed we have found a minimal polynomial,
                # it should be accurate to a much higher precision.
                error = abs(f(aa))
                dx = ~RR(Integer(1) << (check_bits - degree - 2))
                expected_error = abs(f.derivative()(CC(aa))) * dx

                if error < expected_error:
                    # Degree might have been an over-estimate,
                    # factor because we want (irreducible) minpoly.
                    ff = f.factor()
                    for g, e in ff:
                        lead = g.leading_coefficient()
                        if lead != 1:
                            g = g / lead
                        expected_error = abs(g.derivative()(CC(aa))) * dx
                        error = abs(g(aa))
                        if error < expected_error:
                            # See if we can prove equality exactly
                            if g(ex).simplify_trig().canonicalize_radical() == 0:
                                return g
                            # Otherwise fall back to numerical guess
                            elif epsilon and error < epsilon:
                                return g
                            elif algorithm is not None:
                                raise NotImplementedError("Could not prove minimal polynomial %s (epsilon %s)" % (g, RR(error).str(no_sci=False)))

        if algorithm is not None:
            raise ValueError("Could not find minimal polynomial (%s bits, degree %s)." % (bits, degree))

    if algorithm is None or algorithm == 'algebraic':
        from sage.rings.qqbar import QQbar
        return QQ[var](QQbar(ex).minpoly())

    raise ValueError("Unknown algorithm: %s" % algorithm)


###################################################################
# limits
###################################################################
def limit(ex, dir=None, taylor=False, algorithm='maxima', **argv):
    r"""
    Return the limit as the variable `v` approaches `a`
    from the given direction.

    ::

       expr.limit(x = a)
       expr.limit(x = a, dir='+')

    INPUT:

    - ``dir`` -- (default: ``None``) may have the value
      ``'plus'`` (or ``'+'`` or ``'right'`` or ``'above'``) for a limit from above,
      ``'minus'`` (or ``'-'`` or ``'left'`` or ``'below'``) for a limit from below, or may be omitted
      (implying a two-sided limit is to be computed).

    - ``taylor`` -- (default: ``False``) if ``True``, use Taylor
      series, which allows more limits to be computed (but may also
      crash in some obscure cases due to bugs in Maxima).

    - ``**argv`` -- 1 named parameter

    .. NOTE::

        The output may also use ``und`` (undefined), ``ind``
        (indefinite but bounded), and ``infinity`` (complex
        infinity).

    EXAMPLES::

        sage: x = var('x')
        sage: f = (1 + 1/x)^x
        sage: f.limit(x=oo)
        e
        sage: f.limit(x=5)
        7776/3125

    Domain to real, a regression in 5.46.0, see https://sf.net/p/maxima/bugs/4138 ::

        sage: maxima_calculus.eval("domain:real")
        ...
        sage: f.limit(x=1.2).n()
        2.06961575467...
        sage: maxima_calculus.eval("domain:complex");
        ...

    Otherwise, it works ::

        sage: f.limit(x=I, taylor=True)
        (-I + 1)^I
        sage: f(x=1.2)
        2.0696157546720...
        sage: f(x=I)
        (-I + 1)^I
        sage: CDF(f(x=I))
        2.0628722350809046 + 0.7450070621797239*I
        sage: CDF(f.limit(x=I))
        2.0628722350809046 + 0.7450070621797239*I

    Notice that Maxima may ask for more information::

        sage: var('a')
        a
        sage: limit(x^a,x=0)
        Traceback (most recent call last):
        ...
        ValueError: Computation failed since Maxima requested additional
        constraints; using the 'assume' command before evaluation
        *may* help (example of legal syntax is 'assume(a>0)', see
        `assume?` for more details)
        Is a positive, negative or zero?

    With this example, Maxima is looking for a LOT of information::

        sage: assume(a>0)
        sage: limit(x^a,x=0)  # random - maxima 5.46.0 does not need extra assumption
        Traceback (most recent call last):
        ...
        ValueError: Computation failed since Maxima requested additional
        constraints; using the 'assume' command before evaluation *may* help
        (example of legal syntax is 'assume(a>0)', see `assume?` for
         more details)
        Is a an integer?
        sage: assume(a,'integer')
        sage: limit(x^a, x=0)  # random - maxima 5.46.0 does not need extra assumption
        Traceback (most recent call last):
        ...
        ValueError: Computation failed since Maxima requested additional
        constraints; using the 'assume' command before evaluation *may* help
        (example of legal syntax is 'assume(a>0)', see `assume?` for
         more details)
        Is a an even number?
        sage: assume(a, 'even')
        sage: limit(x^a, x=0)
        0
        sage: forget()

    More examples::

        sage: limit(x*log(x), x=0, dir='+')
        0
        sage: lim((x+1)^(1/x), x=0)
        e
        sage: lim(e^x/x, x=oo)
        +Infinity
        sage: lim(e^x/x, x=-oo)
        0
        sage: lim(-e^x/x, x=oo)
        -Infinity
        sage: lim((cos(x))/(x^2), x=0)
        +Infinity
        sage: lim(sqrt(x^2+1) - x, x=oo)
        0
        sage: lim(x^2/(sec(x)-1), x=0)
        2
        sage: lim(cos(x)/(cos(x)-1), x=0)
        -Infinity
        sage: lim(x*sin(1/x), x=0)
        0
        sage: limit(e^(-1/x), x=0, dir='right')
        0
        sage: limit(e^(-1/x), x=0, dir='left')
        +Infinity

    ::

        sage: f = log(log(x)) / log(x)
        sage: forget(); assume(x < -2); lim(f, x=0, taylor=True)
        0
        sage: forget()

    Here ind means "indefinite but bounded"::

        sage: lim(sin(1/x), x = 0)
        ind

    We can use other packages than maxima, namely "sympy", "giac", "fricas".

    With the standard package Giac::

        sage: # needs sage.libs.giac
        sage: (exp(-x)/(2+sin(x))).limit(x=oo, algorithm='giac')
        0
        sage: limit(e^(-1/x), x=0, dir='right', algorithm='giac')
        0
        sage: limit(e^(-1/x), x=0, dir='left', algorithm='giac')
        +Infinity
        sage: (x / (x+2^x+cos(x))).limit(x=-infinity, algorithm='giac')
        1

    With the optional package FriCAS::

        sage: (x / (x+2^x+cos(x))).limit(x=-infinity, algorithm='fricas')       # optional - fricas
        1
        sage: limit(e^(-1/x), x=0, dir='right', algorithm='fricas')             # optional - fricas
        0
        sage: limit(e^(-1/x), x=0, dir='left', algorithm='fricas')              # optional - fricas
        +Infinity

    One can also call Mathematica's online interface::

        sage: limit(pi+log(x)/x,x=oo, algorithm='mathematica_free') # optional - internet
        pi

    TESTS::

        sage: lim(x^2, x=2, dir='nugget')
        Traceback (most recent call last):
        ...
        ValueError: dir must be one of None, 'plus', '+', 'above', 'right',
        'minus', '-', 'below', 'left'

        sage: x.limit(x=3, algorithm='nugget')
        Traceback (most recent call last):
        ...
        ValueError: Unknown algorithm: nugget

    We check that :issue:`3718` is fixed, so that
    Maxima gives correct limits for the floor function::

        sage: limit(floor(x), x=0, dir='-')
        -1
        sage: limit(floor(x), x=0, dir='+')
        0
        sage: limit(floor(x), x=0)
        ...nd

    Maxima gives the right answer here, too, showing
    that :issue:`4142` is fixed::

        sage: f = sqrt(1 - x^2)
        sage: g = diff(f, x); g
        -x/sqrt(-x^2 + 1)
        sage: limit(g, x=1, dir='-')
        -Infinity

    ::

        sage: limit(1/x, x=0)
        Infinity
        sage: limit(1/x, x=0, dir='+')
        +Infinity
        sage: limit(1/x, x=0, dir='-')
        -Infinity

    Check that :issue:`8942` is fixed::

        sage: f(x) = (cos(pi/4 - x) - tan(x)) / (1 - sin(pi/4 + x))
        sage: limit(f(x), x=pi/4, dir='minus')
        +Infinity
        sage: limit(f(x), x=pi/4, dir='plus')
        -Infinity
        sage: limit(f(x), x=pi/4)
        Infinity

    Check that :issue:`12708` is fixed::

        sage: limit(tanh(x), x=0)
        0

    Check that :issue:`15386` is fixed::

        sage: n = var('n')
        sage: assume(n>0)
        sage: sequence = -(3*n^2 + 1)*(-1)^n / sqrt(n^5 + 8*n^3 + 8)
        sage: limit(sequence, n=infinity)
        0

    Check if :issue:`23048` is fixed::

        sage: (1/(x-3)).limit(x=3, dir='below')
        -Infinity

    From :issue:`14677`::

        sage: f = (x^x - sin(x)^sin(x)) / (x^3*log(x))
        sage: limit(f, x=0, algorithm='fricas')                                 # optional - fricas
        und

        sage: limit(f, x=0, dir='right', algorithm='fricas')                    # optional - fricas
        1/6

    From :issue:`26497`::

        sage: mu, y, sigma = var("mu, y, sigma")
        sage: f = 1/2*sqrt(2)*e^(-1/2*(mu - log(y))^2/sigma^2)/(sqrt(pi)*sigma*y)
        sage: limit(f, y=0, algorithm='fricas')                                 # optional - fricas
        0

    From :issue:`26060`::

        sage: limit(x / (x + 2^x + cos(x)), x=-infinity)
        1
    """
    if not isinstance(ex, Expression):
        ex = SR(ex)

    if len(argv) != 1:
        raise ValueError("call the limit function like this, e.g. limit(expr, x=2).")
    else:
        k, = argv.keys()
        v = var(k)
        a = argv[k]

    if taylor and algorithm == 'maxima':
        algorithm = 'maxima_taylor'

    dir_plus = ['plus', '+', 'above', 'right']
    dir_minus = ['minus', '-', 'below', 'left']
    dir_both = [None] + dir_plus + dir_minus
    if dir not in dir_both:
        raise ValueError("dir must be one of " + ", ".join(map(repr, dir_both)))

    if algorithm == 'maxima':
        if dir is None:
            l = maxima.sr_limit(ex, v, a)
        elif dir in dir_plus:
            l = maxima.sr_limit(ex, v, a, 'plus')
        elif dir in dir_minus:
            l = maxima.sr_limit(ex, v, a, 'minus')
    elif algorithm == 'maxima_taylor':
        if dir is None:
            l = maxima.sr_tlimit(ex, v, a)
        elif dir in dir_plus:
            l = maxima.sr_tlimit(ex, v, a, 'plus')
        elif dir in dir_minus:
            l = maxima.sr_tlimit(ex, v, a, 'minus')
    elif algorithm == 'sympy':
        import sympy
        if dir is None:
            l = sympy.limit(ex._sympy_(), v._sympy_(), a._sympy_())
        elif dir in dir_plus:
            l = sympy.limit(ex._sympy_(), v._sympy_(), a._sympy_(), dir='+')
        elif dir in dir_minus:
            l = sympy.limit(ex._sympy_(), v._sympy_(), a._sympy_(), dir='-')
    elif algorithm == 'fricas':
        from sage.interfaces.fricas import fricas
        eq = fricas.equation(v._fricas_(), a._fricas_())
        f = ex._fricas_()
        if dir is None:
            l = fricas.limit(f, eq).sage()
            if isinstance(l, dict):
                l = maxima("und")
        elif dir in dir_plus:
            l = fricas.limit(f, eq, '"right"').sage()
        elif dir in dir_minus:
            l = fricas.limit(f, eq, '"left"').sage()
    elif algorithm == 'giac':
        from sage.libs.giac.giac import libgiac
        v = v._giac_init_()
        a = a._giac_init_()
        if dir is None:
            l = libgiac.limit(ex, v, a).sage()
        elif dir in dir_plus:
            l = libgiac.limit(ex, v, a, 1).sage()
        elif dir in dir_minus:
            l = libgiac.limit(ex, v, a, -1).sage()
    elif algorithm == 'mathematica_free':
        return mma_free_limit(ex, v, a, dir)
    else:
        raise ValueError("Unknown algorithm: %s" % algorithm)
    return ex.parent()(l)


# lim is alias for limit
lim = limit


def mma_free_limit(expression, v, a, dir=None):
    """
    Limit using Mathematica's online interface.

    INPUT:

    - ``expression`` -- symbolic expression
    - ``v`` -- variable
    - ``a`` -- value where the variable goes to
    - ``dir`` -- ``'+'``, ``'-'`` or ``None`` (default: ``None``)

    EXAMPLES::

        sage: from sage.calculus.calculus import mma_free_limit
        sage: mma_free_limit(sin(x)/x, x, a=0) # optional - internet
        1

    Another simple limit::

        sage: mma_free_limit(e^(-x), x, a=oo) # optional - internet
        0
    """
    from sage.interfaces.mathematica import request_wolfram_alpha, parse_moutput_from_json, symbolic_expression_from_mathematica_string
    dir_plus = ['plus', '+', 'above', 'right']
    dir_minus = ['minus', '-', 'below', 'left']
    math_expr = expression._mathematica_init_()
    variable = v._mathematica_init_()
    a = a._mathematica_init_()
    if dir is None:
        input = "Limit[{},{} -> {}]".format(math_expr, variable, a)
    elif dir in dir_plus:
        dir = 'Direction -> "FromAbove"'
        input = "Limit[{}, {} -> {}, {}]".format(math_expr, variable, a, dir)
    elif dir in dir_minus:
        dir = 'Direction -> "FromBelow"'
        input = "Limit[{}, {} -> {}, {}]".format(math_expr, variable, a, dir)
    else:
        raise ValueError('wrong input for limit')
    json_page_data = request_wolfram_alpha(input)
    all_outputs = parse_moutput_from_json(json_page_data)
    if not all_outputs:
        raise ValueError("no outputs found in the answer from Wolfram Alpha")
    first_output = all_outputs[0]
    return symbolic_expression_from_mathematica_string(first_output)


###################################################################
# Laplace transform
###################################################################
def laplace(ex, t, s, algorithm='maxima'):
    r"""
    Return the Laplace transform with respect to the variable `t` and
    transform parameter `s`, if possible.

    If this function cannot find a solution, a formal function is returned.
    The function that is returned may be viewed as a function of `s`.

    DEFINITION:

    The Laplace transform of a function `f(t)`, defined for all real numbers
    `t \geq 0`, is the function `F(s)` defined by

    .. MATH::

                      F(s) = \int_{0}^{\infty} e^{-st} f(t) dt.

    INPUT:

    - ``ex`` -- a symbolic expression

    - ``t`` -- independent variable

    - ``s`` -- transform parameter

    - ``algorithm`` -- (default: ``'maxima'``)  one of

      - ``'maxima'`` -- use Maxima (the default)

      - ``'sympy'`` -- use SymPy

      - ``'giac'`` -- use Giac (optional)

    .. NOTE::

        The ``'sympy'`` algorithm returns the tuple (`F`, `a`, ``cond``)
        where `F` is the Laplace transform of `f(t)`,
        `Re(s)>a` is the half-plane of convergence, and ``cond`` are
        auxiliary convergence conditions.

    .. SEEALSO::

        :func:`inverse_laplace`

    EXAMPLES:

    We compute a few Laplace transforms::

        sage: var('x, s, z, t, t0')
        (x, s, z, t, t0)
        sage: sin(x).laplace(x, s)
        1/(s^2 + 1)
        sage: (z + exp(x)).laplace(x, s)
        z/s + 1/(s - 1)
        sage: log(t/t0).laplace(t, s)
        -(euler_gamma + log(s) + log(t0))/s

    We do a formal calculation::

        sage: f = function('f')(x)
        sage: g = f.diff(x); g
        diff(f(x), x)
        sage: g.laplace(x, s)
        s*laplace(f(x), x, s) - f(0)

    A BATTLE BETWEEN the X-women and the Y-men (by David
    Joyner): Solve

    .. MATH::

                   x' = -16y, x(0)=270,  y' = -x + 1, y(0) = 90.

    This models a fight between two sides, the "X-women" and the
    "Y-men", where the X-women have 270 initially and the Y-men have
    90, but the Y-men are better at fighting, because of the higher
    factor of "-16" vs "-1", and also get an occasional reinforcement,
    because of the "+1" term.

    ::

        sage: var('t')
        t
        sage: t = var('t')
        sage: x = function('x')(t)
        sage: y = function('y')(t)
        sage: de1 = x.diff(t) + 16*y
        sage: de2 = y.diff(t) + x - 1
        sage: de1.laplace(t, s)
        s*laplace(x(t), t, s) + 16*laplace(y(t), t, s) - x(0)
        sage: de2.laplace(t, s)
        s*laplace(y(t), t, s) - 1/s + laplace(x(t), t, s) - y(0)

    Next we form the augmented matrix of the above system::

        sage: A = matrix([[s, 16, 270], [1, s, 90+1/s]])
        sage: E = A.echelon_form()
        sage: xt = E[0,2].inverse_laplace(s,t)
        sage: yt = E[1,2].inverse_laplace(s,t)
        sage: xt
        -91/2*e^(4*t) + 629/2*e^(-4*t) + 1
        sage: yt
        91/8*e^(4*t) + 629/8*e^(-4*t)
        sage: p1 = plot(xt, 0, 1/2, rgbcolor=(1,0,0))                                   # needs sage.plot
        sage: p2 = plot(yt, 0, 1/2, rgbcolor=(0,1,0))                                   # needs sage.plot
        sage: import tempfile
        sage: with tempfile.NamedTemporaryFile(suffix='.png') as f:                     # needs sage.plot
        ....:     (p1 + p2).save(f.name)

    Another example::

        sage: var('a,s,t')
        (a, s, t)
        sage: f = exp (2*t + a) * sin(t) * t; f
        t*e^(a + 2*t)*sin(t)
        sage: L = laplace(f, t, s); L
        2*(s - 2)*e^a/(s^2 - 4*s + 5)^2
        sage: inverse_laplace(L, s, t)
        t*e^(a + 2*t)*sin(t)

    The Laplace transform of the exponential function::

        sage: laplace(exp(x), x, s)
        1/(s - 1)

    Dirac's delta distribution is handled (the output of SymPy is
    related to a choice that has to be made when defining Laplace
    transforms of distributions)::

        sage: laplace(dirac_delta(t), t, s)
        1
        sage: F, a, cond = laplace(dirac_delta(t), t, s, algorithm='sympy')
        sage: a, cond  # random - sympy <1.10 gives (-oo, True)
        (0, True)
        sage: F        # random - sympy <1.9 includes undefined heaviside(0) in answer
        1
        sage: laplace(dirac_delta(t), t, s, algorithm='giac')  # needs giac
        1

    Heaviside step function can be handled with different interfaces.
    Try with Maxima::

        sage: laplace(heaviside(t-1), t, s)
        e^(-s)/s

    Try with giac, if it is installed::

        sage: # needs giac
        sage: laplace(heaviside(t-1), t, s, algorithm='giac')
        e^(-s)/s

    Try with SymPy::

        sage: laplace(heaviside(t-1), t, s, algorithm='sympy')
        (e^(-s)/s, 0, True)

    TESTS:

    Testing Giac::

        sage: # needs giac
        sage: var('t, s')
        (t, s)
        sage: laplace(5*cos(3*t-2)*heaviside(t-2), t, s, algorithm='giac')
        5*(s*cos(4)*e^(-2*s) - 3*e^(-2*s)*sin(4))/(s^2 + 9)

    Check unevaluated expression from Giac (it is locale-dependent, see
    :issue:`22833`)::

        sage: # needs giac
        sage: n = SR.var('n')
        sage: laplace(t^n, t, s, algorithm='giac')
        laplace(t^n, t, s)

    Testing SymPy::

        sage: n = SR.var('n')
        sage: F, a, cond = laplace(t^n, t, s, algorithm='sympy')
        sage: a, cond
        (0, re(n) > -1)
        sage: F.simplify()
        s^(-n - 1)*gamma(n + 1)

    Testing Maxima::

        sage: laplace(t^n, t, s, algorithm='maxima')
        s^(-n - 1)*gamma(n + 1)

    Check that :issue:`24212` is fixed::

        sage: F, a, cond = laplace(cos(t^2), t, s, algorithm='sympy')
        sage: a, cond
        (0, True)
        sage: F._sympy_().simplify()
        sqrt(pi)*(sqrt(2)*sin(s**2/4)*fresnelc(sqrt(2)*s/(2*sqrt(pi))) -
        sqrt(2)*cos(s**2/4)*fresnels(sqrt(2)*s/(2*sqrt(pi))) + cos(s**2/4 + pi/4))/2

    Testing result from SymPy that Sage doesn't know how to handle::

        sage: laplace(cos(-1/t), t, s, algorithm='sympy')
        Traceback (most recent call last):
        ...
        AttributeError: Unable to convert SymPy result (=meijerg(((), ()), ((-1/2, 0, 1/2), (0,)), ...)/4) into Sage
    """
    if not isinstance(ex, (Expression, Function)):
        ex = SR(ex)

    if algorithm == 'maxima':
        return ex.parent()(ex._maxima_().laplace(var(t), var(s)))

    elif algorithm == 'sympy':
        ex_sy, t, s = (expr._sympy_() for expr in (ex, t, s))
        from sympy import laplace_transform
        from sage.interfaces.sympy import sympy_init
        sympy_init()
        result = laplace_transform(ex_sy, t, s)
        if isinstance(result, tuple):
            try:
                (result, a, cond) = result
                return result._sage_(), a, cond
            except AttributeError:
                raise AttributeError("Unable to convert SymPy result (={}) into"
                        " Sage".format(result))
        elif 'LaplaceTransform' in format(result):
            return dummy_laplace(ex, t, s)
        else:
            return result

    elif algorithm == 'giac':
        from sage.interfaces.giac import giac
        try:
            result = giac.laplace(ex, t, s)
        except TypeError:
            raise ValueError("Giac cannot make sense of: %s" % ex)
        if 'integrate' in format(result) or 'integration' in format(result):
            return dummy_laplace(ex, t, s)
        else:
            return result.sage()

    else:
        raise ValueError("Unknown algorithm: %s" % algorithm)


def inverse_laplace(ex, s, t, algorithm='maxima'):
    r"""
    Return the inverse Laplace transform with respect to the variable `t` and
    transform parameter `s`, if possible.

    If this function cannot find a solution, a formal function is returned.
    The function that is returned may be viewed as a function of `t`.

    DEFINITION:

    The inverse Laplace transform of a function `F(s)` is the function
    `f(t)`, defined by

    .. MATH::

                      f(t) = \frac{1}{2\pi i} \int_{\gamma-i\infty}^{\gamma + i\infty} e^{st} F(s) ds,

    where `\gamma` is chosen so that the contour path of
    integration is in the region of convergence of `F(s)`.

    INPUT:

    - ``ex`` -- a symbolic expression

    - ``s`` -- transform parameter

    - ``t`` -- independent variable

    - ``algorithm`` -- (default: ``'maxima'``)  one of

      - ``'maxima'`` -- use Maxima (the default)

      - ``'sympy'`` -- use SymPy

      - ``'giac'`` -- use Giac (optional)

    .. SEEALSO::

        :func:`laplace`

    EXAMPLES::

        sage: var('w, m')
        (w, m)
        sage: f = (1/(w^2+10)).inverse_laplace(w, m); f
        1/10*sqrt(10)*sin(sqrt(10)*m)
        sage: laplace(f, m, w)
        1/(w^2 + 10)

        sage: f(t) = t*cos(t)
        sage: s = var('s')
        sage: L = laplace(f, t, s); L
        t |--> 2*s^2/(s^2 + 1)^2 - 1/(s^2 + 1)
        sage: inverse_laplace(L, s, t)
        t |--> t*cos(t)
        sage: inverse_laplace(1/(s^3+1), s, t)
        1/3*(sqrt(3)*sin(1/2*sqrt(3)*t) - cos(1/2*sqrt(3)*t))*e^(1/2*t) + 1/3*e^(-t)

    No explicit inverse Laplace transform, so one is returned formally a
    function ``ilt``::

        sage: inverse_laplace(cos(s), s, t)
        ilt(cos(s), s, t)

    Transform an expression involving a time-shift, via SymPy::

        sage: inverse_laplace(1/s^2*exp(-s), s, t, algorithm='sympy').simplify()
        (t - 1)*heaviside(t - 1)

    The same instance with Giac::

        sage: # needs giac
        sage: inverse_laplace(1/s^2*exp(-s), s, t, algorithm='giac')
        (t - 1)*heaviside(t - 1)

    Transform a rational expression::

        sage: # needs giac
        sage: inverse_laplace((2*s^2*exp(-2*s) - exp(-s))/(s^3+1), s, t,
        ....:                 algorithm='giac')
        -1/3*(sqrt(3)*e^(1/2*t - 1/2)*sin(1/2*sqrt(3)*(t - 1))
         - cos(1/2*sqrt(3)*(t - 1))*e^(1/2*t - 1/2) + e^(-t + 1))*heaviside(t - 1)
         + 2/3*(2*cos(1/2*sqrt(3)*(t - 2))*e^(1/2*t - 1) + e^(-t + 2))*heaviside(t - 2)

        sage: inverse_laplace(1/(s - 1), s, x)
        e^x

    The inverse Laplace transform of a constant is a delta
    distribution::

        sage: inverse_laplace(1, s, t)
        dirac_delta(t)
        sage: inverse_laplace(1, s, t, algorithm='sympy')
        dirac_delta(t)
        sage: inverse_laplace(1, s, t, algorithm='giac')  # needs giac
        dirac_delta(t)

    TESTS:

    Testing unevaluated expression from Maxima::

        sage: var('t, s')
        (t, s)
        sage: inverse_laplace(exp(-s)/s, s, t)
        ilt(e^(-s)/s, s, t)

    Testing Giac::

        sage: # needs giac
        sage: inverse_laplace(exp(-s)/s, s, t, algorithm='giac')
        heaviside(t - 1)

    Testing SymPy::

        sage: inverse_laplace(exp(-s)/s, s, t, algorithm='sympy')
        heaviside(t - 1)

    Testing unevaluated expression from Giac::

        sage: # needs giac
        sage: n = SR.var('n')
        sage: inverse_laplace(1/s^n, s, t, algorithm='giac')
        ilt(1/(s^n), t, s)

    Try with Maxima::

        sage: n = SR.var('n')
        sage: inverse_laplace(1/s^n, s, t, algorithm='maxima')
        ilt(1/(s^n), s, t)

    Try with SymPy::

        sage: inverse_laplace(1/s^n, s, t, algorithm='sympy')
        t^(n - 1)*heaviside(t)/gamma(n)

    Testing unevaluated expression from SymPy::

        sage: inverse_laplace(cos(s), s, t, algorithm='sympy')
        ilt(cos(s), t, s)

    Testing the same with Giac::

        sage: # needs giac
        sage: inverse_laplace(cos(s), s, t, algorithm='giac')
        ilt(cos(s), t, s)
    """
    if not isinstance(ex, Expression):
        ex = SR(ex)

    if algorithm == 'maxima':
        return ex.parent()(ex._maxima_().ilt(var(s), var(t)))

    elif algorithm == 'sympy':
        ex_sy, s, t = (expr._sympy_() for expr in (ex, s, t))
        from sympy import inverse_laplace_transform
        from sage.interfaces.sympy import sympy_init
        sympy_init()
        result = inverse_laplace_transform(ex_sy, s, t)
        try:
            return result._sage_()
        except AttributeError:
            if 'InverseLaplaceTransform' in format(result):
                return dummy_inverse_laplace(ex, t, s)
            else:
                raise AttributeError("Unable to convert SymPy result (={}) into"
                                    " Sage".format(result))

    elif algorithm == 'giac':
        from sage.interfaces.giac import giac
        try:
            result = giac.invlaplace(ex, s, t)
        except TypeError:
            raise ValueError("Giac cannot make sense of: %s" % ex)
        if 'ilaplace' in format(result):
            return dummy_inverse_laplace(ex, t, s)
        else:
            return result.sage()

    else:
        raise ValueError("Unknown algorithm: %s" % algorithm)


###################################################################
# symbolic evaluation "at" a point
###################################################################
def at(ex, *args, **kwds):
    """
    Parses ``at`` formulations from other systems, such as Maxima.
    Replaces evaluation 'at' a point with substitution method of
    a symbolic expression.

    EXAMPLES:

    We do not import ``at`` at the top level, but we can use it
    as a synonym for substitution if we import it::

        sage: g = x^3 - 3
        sage: from sage.calculus.calculus import at
        sage: at(g, x=1)
        -2
        sage: g.subs(x=1)
        -2

    We find a formal Taylor expansion::

        sage: h,x = var('h,x')
        sage: u = function('u')
        sage: u(x + h)
        u(h + x)
        sage: diff(u(x+h), x)
        D[0](u)(h + x)
        sage: taylor(u(x+h), h, 0, 4)
        1/24*h^4*diff(u(x), x, x, x, x) + 1/6*h^3*diff(u(x), x, x, x)
         + 1/2*h^2*diff(u(x), x, x) + h*diff(u(x), x) + u(x)

    We compute a Laplace transform::

        sage: var('s,t')
        (s, t)
        sage: f = function('f')(t)
        sage: f.diff(t, 2)
        diff(f(t), t, t)
        sage: f.diff(t,2).laplace(t,s)
        s^2*laplace(f(t), t, s) - s*f(0) - D[0](f)(0)

    We can also accept a non-keyword list of expression substitutions,
    like Maxima does (:issue:`12796`)::

        sage: from sage.calculus.calculus import at
        sage: f = function('f')
        sage: at(f(x), [x == 1])
        f(1)

    TESTS:

    Our one non-keyword argument must be a list::

        sage: from sage.calculus.calculus import at
        sage: f = function('f')
        sage: at(f(x), x == 1)
        Traceback (most recent call last):
        ...
        TypeError: at can take at most one argument, which must be a list

    We should convert our first argument to a symbolic expression::

        sage: from sage.calculus.calculus import at
        sage: at(int(1), x=1)
        1
    """
    if not isinstance(ex, (Expression, Function)):
        ex = SR(ex)
    kwds = {(k[10:] if k[:10] == "_SAGE_VAR_" else k): v
            for k, v in kwds.items()}
    if len(args) == 1 and isinstance(args[0], list):
        for c in args[0]:
            kwds[str(c.lhs())] = c.rhs()
    elif args:
        raise TypeError("at can take at most one argument, which must be a list")

    return ex.subs(**kwds)


def dummy_diff(*args):
    """
    This function is called when 'diff' appears in a Maxima string.

    EXAMPLES::

        sage: from sage.calculus.calculus import dummy_diff
        sage: x,y = var('x,y')
        sage: dummy_diff(sin(x*y), x, SR(2), y, SR(1))
        -x*y^2*cos(x*y) - 2*y*sin(x*y)

    Here the function is used implicitly::

        sage: a = var('a')
        sage: f = function('cr')(a)
        sage: g = f.diff(a); g
        diff(cr(a), a)
    """
    f = args[0]
    args = list(args[1:])
    for i in range(1, len(args), 2):
        args[i] = Integer(args[i])
    return f.diff(*args)


def dummy_integrate(*args):
    """
    This function is called to create formal wrappers of integrals that
    Maxima can't compute:

    EXAMPLES::

        sage: from sage.calculus.calculus import dummy_integrate
        sage: f = function('f')
        sage: dummy_integrate(f(x), x)
        integrate(f(x), x)
        sage: a,b = var('a,b')
        sage: dummy_integrate(f(x), x, a, b)
        integrate(f(x), x, a, b)
    """
    if len(args) == 4:
        return definite_integral(*args, hold=True)
    else:
        return indefinite_integral(*args, hold=True)


def dummy_laplace(*args):
    """
    This function is called to create formal wrappers of laplace transforms
    that Maxima can't compute:

    EXAMPLES::

        sage: from sage.calculus.calculus import dummy_laplace
        sage: s,t = var('s,t')
        sage: f = function('f')
        sage: dummy_laplace(f(t), t, s)
        laplace(f(t), t, s)
    """
    return _laplace(args[0], var(repr(args[1])), var(repr(args[2])))


def dummy_inverse_laplace(*args):
    """
    This function is called to create formal wrappers of inverse Laplace
    transforms that Maxima can't compute:

    EXAMPLES::

        sage: from sage.calculus.calculus import dummy_inverse_laplace
        sage: s,t = var('s,t')
        sage: F = function('F')
        sage: dummy_inverse_laplace(F(s), s, t)
        ilt(F(s), s, t)
    """
    return _inverse_laplace(args[0], var(repr(args[1])), var(repr(args[2])))


def dummy_pochhammer(*args):
    """
    This function is called to create formal wrappers of Pochhammer symbols.

    EXAMPLES::

        sage: from sage.calculus.calculus import dummy_pochhammer
        sage: s,t = var('s,t')
        sage: dummy_pochhammer(s, t)
        gamma(s + t)/gamma(s)
    """
    x, y = args
    from sage.functions.gamma import gamma
    return gamma(x + y) / gamma(x)


#######################################################
#
# Helper functions for printing latex expression
#
#######################################################

def _laplace_latex_(self, *args):
    r"""
    Return LaTeX expression for Laplace transform of a symbolic function.

    EXAMPLES::

        sage: from sage.calculus.calculus import _laplace_latex_
        sage: var('s,t')
        (s, t)
        sage: f = function('f')(t)
        sage: _laplace_latex_(0,f,t,s)
        '\\mathcal{L}\\left(f\\left(t\\right), t, s\\right)'
        sage: latex(laplace(f, t, s))
        \mathcal{L}\left(f\left(t\right), t, s\right)
    """
    return "\\mathcal{L}\\left(%s\\right)" % (', '.join(latex(x) for x in args))


def _inverse_laplace_latex_(self, *args):
    r"""
    Return LaTeX expression for inverse Laplace transform
    of a symbolic function.

    EXAMPLES::

        sage: from sage.calculus.calculus import _inverse_laplace_latex_
        sage: var('s,t')
        (s, t)
        sage: F = function('F')(s)
        sage: _inverse_laplace_latex_(0,F,s,t)
        '\\mathcal{L}^{-1}\\left(F\\left(s\\right), s, t\\right)'
        sage: latex(inverse_laplace(F,s,t))
        \mathcal{L}^{-1}\left(F\left(s\right), s, t\right)
    """
    return "\\mathcal{L}^{-1}\\left(%s\\right)" % (', '.join(latex(x) for x in args))


# Return un-evaluated expression as instances of NewSymbolicFunction
_laplace = function_factory('laplace', print_latex_func=_laplace_latex_)
_inverse_laplace = function_factory('ilt',
        print_latex_func=_inverse_laplace_latex_)

######################################i################


# Conversion dict for special maxima objects
# c,k1,k2 are from ode2()
symtable = {'%pi': 'pi', '%e': 'e', '%i': 'I',
            '%gamma': 'euler_gamma',
            '%c': '_C', '%k1': '_K1', '%k2': '_K2',
            'e': '_e', 'i': '_i', 'I': '_I'}


maxima_qp = re.compile(r"\?\%[\w]*")  # e.g., ?%jacobi_cd

maxima_var = re.compile(r"[\w\%]*")  # e.g., %jacobi_cd

sci_not = re.compile(r"(-?(?:0|[1-9]\d*))(\.\d+)?([eE][-+]\d+)")

polylog_ex = re.compile(r'li\[([^\[\]]*)\]\(')

maxima_polygamma = re.compile(r"psi\[([^\[\]]*)\]\(")  # matches psi[n]( where n is a number

maxima_hyper = re.compile(r"\%f\[\d+,\d+\]")  # matches %f[m,n]


def _is_function(v):
    r"""
    Return whether a symbolic element is a function, not a variable.

    TESTS::

        sage: from sage.calculus.calculus import _is_function
        sage: _is_function(x)
        False
        sage: _is_function(sin)
        True

    Check that :issue:`31756` is fixed::

        sage: from sage.symbolic.expression import symbol_table
        sage: _is_function(symbol_table['mathematica'][('Gamma', 1)])
        True

        sage: from sage.symbolic.expression import register_symbol
        sage: foo = lambda x: x^2 + 1
        sage: register_symbol(foo, dict(mathematica='Foo'))  # optional - mathematica
        sage: mathematica('Foo[x]').sage()                   # optional - mathematica
        x^2 + 1
    """
    # note that Sage variables are callable, so we only check the type
    return isinstance(v, (Function, FunctionType))


def symbolic_expression_from_maxima_string(x, equals_sub=False, maxima=maxima):
    r"""
    Given a string representation of a Maxima expression, parse it and
    return the corresponding Sage symbolic expression.

    INPUT:

    - ``x`` -- string

    - ``equals_sub`` -- boolean (default: ``False``); if ``True``, replace
      '=' by '==' in self

    - ``maxima`` -- (default: the calculus package's copy of
      Maxima) the Maxima interpreter to use

    EXAMPLES::

        sage: from sage.calculus.calculus import symbolic_expression_from_maxima_string as sefms
        sage: sefms('x^%e + %e^%pi + %i + sin(0)')
        x^e + e^pi + I
        sage: f = function('f')(x)
        sage: sefms('?%at(f(x),x=2)#1')
        f(2) != 1
        sage: a = sage.calculus.calculus.maxima("x#0"); a
        x # 0
        sage: a.sage()
        x != 0

    TESTS:

    :issue:`8459` fixed::

        sage: maxima('3*li[2](u)+8*li[33](exp(u))').sage()
        3*dilog(u) + 8*polylog(33, e^u)

    Check if :issue:`8345` is fixed::

        sage: assume(x,'complex')
        sage: t = x.conjugate()
        sage: latex(t)
        \overline{x}
        sage: latex(t._maxima_()._sage_())
        \overline{x}

    Check that we can understand maxima's not-equals (:issue:`8969`)::

        sage: from sage.calculus.calculus import symbolic_expression_from_maxima_string as sefms
        sage: sefms("x!=3") == (factorial(x) == 3)
        True
        sage: sefms("x # 3") == SR(x != 3)
        True
        sage: solve([x != 5], x) in [[[x - 5 != 0]], [[x < 5], [5 < x]]]
        True
        sage: solve([2*x==3, x != 5], x)
        [[x == (3/2)...

    Make sure that we don't accidentally pick up variables in the maxima namespace (:issue:`8734`)::

        sage: sage.calculus.calculus.maxima('my_new_var : 2')
        2
        sage: var('my_new_var').full_simplify()
        my_new_var

    ODE solution constants are treated differently (:issue:`16007`)::

        sage: from sage.calculus.calculus import symbolic_expression_from_maxima_string as sefms
        sage: sefms('%k1*x + %k2*y + %c')
        _K1*x + _K2*y + _C

    Check that some hypothetical variables don't end up as special constants (:issue:`6882`)::

        sage: from sage.calculus.calculus import symbolic_expression_from_maxima_string as sefms
        sage: sefms('%i')^2
        -1
        sage: ln(sefms('%e'))
        1
        sage: sefms('i')^2
        _i^2
        sage: sefms('I')^2
        _I^2
        sage: sefms('ln(e)')
        ln(_e)
        sage: sefms('%inf')
        +Infinity
    """
    var_syms = {k[0]: v for k, v in symbol_table.get('maxima', {}).items()
                if not _is_function(v)}
    function_syms = {k[0]: v for k, v in symbol_table.get('maxima', {}).items()
                     if _is_function(v)}

    if not x:
        raise RuntimeError("invalid symbolic expression -- ''")
    maxima.set('_tmp_', x)

    # This is inefficient since it so rarely is needed:
    #r = maxima._eval_line('listofvars(_tmp_);')[1:-1]

    s = maxima._eval_line('_tmp_;')

    # We don't actually implement a parser for maxima expressions.
    # Instead we simply transform the string until it is a valid
    # sagemath expression and parse that.

    # Remove ticks in front of symbolic functions. You might think
    # there is a potential very subtle bug if 'foo is in a string
    # literal -- but string literals should *never* ever be part of a
    # symbolic expression.
    s = s.replace("'","")

    delayed_functions = maxima_qp.findall(s)
    if len(delayed_functions):
        for X in delayed_functions:
            if X == '?%at':  # we will replace Maxima's "at" with symbolic evaluation, not a SymbolicFunction
                pass
            else:
                function_syms[X[2:]] = function_factory(X[2:])
        s = s.replace("?%", "")

    s = maxima_hyper.sub('hypergeometric', s)

    # Look up every variable in the symtable keys and fill a replacement list.
    cursor = 0
    l = []
    for m in maxima_var.finditer(s):
        if m.group(0) in symtable:
            l.append(s[cursor:m.start()])
            l.append(symtable.get(m.group(0)))
            cursor = m.end()
    if cursor > 0:
        l.append(s[cursor:])
        s = "".join(l)

    s = s.replace("%","")

    s = s.replace("#","!=")  # a lot of this code should be refactored somewhere...
    #we apply the square-bracket replacing patterns repeatedly
    #to ensure that nested brackets get handled (from inside to out)
    while True:
        olds = s
        s = polylog_ex.sub('polylog(\\1,', s)
        s = maxima_polygamma.sub(r'psi(\g<1>,', s)  # this replaces psi[n](foo) with psi(n,foo), ensuring that derivatives of the digamma function are parsed properly below
        if s == olds:
            break

    if equals_sub:
        s = s.replace('=', '==')
        # unfortunately, this will turn != into !==, which we correct
        s = s.replace("!==", "!=")

    #replace %union from to_poly_solve with a list
    if s[0:5] == 'union':
        s = s[5:]
        s = s[s.find("(") + 1:s.rfind(")")]
        s = "[" + s + "]"  # turn it into a string that looks like a list

    # replace %solve from to_poly_solve with the expressions
    if s[0:5] == 'solve':
        s = s[5:]
        s = s[s.find("(") + 1:s.find("]") + 1]

    # replace all instances of Maxima's scientific notation
    # with regular notation
    search = sci_not.search(s)
    while search is not None:
        (start, end) = search.span()
        r = create_RealNumber(s[start:end]).str(no_sci=2, truncate=True)
        s = s.replace(s[start:end], r)
        search = sci_not.search(s)

    function_syms['diff'] = dummy_diff
    function_syms['integrate'] = dummy_integrate
    function_syms['laplace'] = dummy_laplace
    function_syms['ilt'] = dummy_inverse_laplace
    function_syms['at'] = at
    function_syms['pochhammer'] = dummy_pochhammer

    global is_simplified
    try:
        # use a global flag so all expressions obtained via
        # evaluation of maxima code are assumed pre-simplified
        is_simplified = True
        SRM_parser._variable_constructor().set_names(var_syms)
        SRM_parser._callable_constructor().set_names(function_syms)
        return SRM_parser.parse_sequence(s)
    except SyntaxError:
        raise TypeError("unable to make sense of Maxima expression '%s' in Sage" % s)
    finally:
        is_simplified = False


# Comma format options for Maxima
def mapped_opts(v):
    """
    Used internally when creating a string of options to pass to
    Maxima.

    INPUT:

    - ``v`` -- an object

    OUTPUT: string

    The main use of this is to turn Python bools into lower case
    strings.

    EXAMPLES::

        sage: sage.calculus.calculus.mapped_opts(True)
        'true'
        sage: sage.calculus.calculus.mapped_opts(False)
        'false'
        sage: sage.calculus.calculus.mapped_opts('bar')
        'bar'
    """
    if isinstance(v, bool):
        return str(v).lower()
    return str(v)


def maxima_options(**kwds):
    """
    Used internally to create a string of options to pass to Maxima.

    EXAMPLES::

        sage: sage.calculus.calculus.maxima_options(an_option=True, another=False, foo='bar')
        'an_option=true,another=false,foo=bar'
    """
    return ','.join('%s=%s' % (key, mapped_opts(val))
                    for key, val in sorted(kwds.items()))


# Parser for symbolic ring elements

# We keep two dictionaries syms_cur and syms_default to keep the current symbol
# table and the state of the table at startup respectively. These are used by
# the restore() function (see sage.misc.reset).

syms_cur = symbol_table.get('functions', {})
syms_default = dict(syms_cur)


def _find_var(name, interface=None):
    """
    Function to pass to Parser for constructing
    variables from strings.  For internal use.

    EXAMPLES::

        sage: y = SR.var('y')
        sage: sage.calculus.calculus._find_var('y')
        y
        sage: sage.calculus.calculus._find_var('I')
        I
        sage: sage.calculus.calculus._find_var(repr(maxima(y)), interface='maxima')
        y

    ::

        sage: # needs giac
        sage: y = SR.var('y')
        sage: sage.calculus.calculus._find_var(repr(giac(y)), interface='giac')
        y
    """
    if interface == 'maxima':
        if name.startswith("_SAGE_VAR_"):
            return var(name[10:])
    elif interface == 'giac':
        if name.startswith('sageVAR'):
            return var(name[7:])
    else:
        v = SR.symbols.get(name)
        if v is not None:
            return v

    # try to find the name in the global namespace
    # needed for identifiers like 'e', etc.
    import sage.all
    try:
        return SR(sage.all.__dict__[name])
    except (KeyError, TypeError):
        return var(name)


def _find_func(name, create_when_missing=True):
    """
    Function to pass to Parser for constructing
    functions from strings.  For internal use.

    EXAMPLES::

        sage: sage.calculus.calculus._find_func('limit')
        limit
        sage: sage.calculus.calculus._find_func('zeta_zeros')
        zeta_zeros
        sage: f(x)=sin(x)
        sage: sage.calculus.calculus._find_func('f')
        f
        sage: sage.calculus.calculus._find_func('g', create_when_missing=False)
        sage: s = sage.calculus.calculus._find_func('sin')
        sage: s(0)
        0
    """
    f = symbol_table['functions'].get(name)
    if f is not None:
        return f

    import sage.all
    try:
        f = SR(sage.all.__dict__[name])
        if not isinstance(f, Expression):
            return f
    except (KeyError, TypeError):
        if create_when_missing:
            return function_factory(name)
        else:
            return None


parser_make_var = LookupNameMaker({}, fallback=_find_var)
parser_make_function = LookupNameMaker({}, fallback=_find_func)

SR_parser = Parser(make_int=lambda x: SR(Integer(x)),
                   make_float=lambda x: SR(create_RealNumber(x)),
                   make_var=parser_make_var,
                   make_function=parser_make_function)


def symbolic_expression_from_string(s, syms=None, accept_sequence=False, *, parser=None):
    """
    Given a string, (attempt to) parse it and return the
    corresponding Sage symbolic expression.  Normally used
    to return Maxima output to the user.

    INPUT:

    - ``s`` -- string

    - ``syms`` -- (default: ``{}``) dictionary of
      strings to be regarded as symbols or functions;
      keys are pairs (string, number of arguments)

    - ``accept_sequence`` -- boolean (default: ``False``); controls whether
      to allow a (possibly nested) set of lists and tuples
      as input

    - ``parser`` -- (default: ``SR_parser``) parser for internal use

    EXAMPLES::

        sage: from sage.calculus.calculus import symbolic_expression_from_string
        sage: y = var('y')
        sage: symbolic_expression_from_string('[sin(0)*x^2,3*spam+e^pi]',
        ....:                                 syms={('spam',0): y}, accept_sequence=True)
        [0, 3*y + e^pi]

    TESTS:

    Check that the precision is preserved (:issue:`28814`)::

        sage: symbolic_expression_from_string(str(RealField(100)(1/3)))
        0.3333333333333333333333333333
        sage: symbolic_expression_from_string(str(RealField(100)(10^-500/3)))
        3.333333333333333333333333333e-501

    The Giac interface uses a different parser (:issue:`30133`)::

        sage: # needs giac
        sage: from sage.calculus.calculus import SR_parser_giac
        sage: symbolic_expression_from_string(repr(giac(SR.var('e'))), parser=SR_parser_giac)
        e
    """
    if syms is None:
        syms = {}
    if parser is None:
        parser = SR_parser
    parse_func = parser.parse_sequence if accept_sequence else parser.parse_expression
    # this assumes that the parser has constructors of type `LookupNameMaker`
    parser._variable_constructor().set_names({k[0]: v for k, v in syms.items()
                                              if not _is_function(v)})
    parser._callable_constructor().set_names({k[0]: v for k, v in syms.items()
                                              if _is_function(v)})
    return parse_func(s)


parser_make_Mvar = LookupNameMaker({}, fallback=lambda x: _find_var(x, interface='maxima'))

SRM_parser = Parser(make_int=lambda x: SR(Integer(x)),
                    make_float=lambda x: SR(RealDoubleElement(x)),
                    make_var=parser_make_Mvar,
                    make_function=parser_make_function)

SR_parser_giac = Parser(make_int=lambda x: SR(Integer(x)),
                        make_float=lambda x: SR(create_RealNumber(x)),
                        make_var=LookupNameMaker({}, fallback=lambda x: _find_var(x, interface='giac')),
                        make_function=parser_make_function)
