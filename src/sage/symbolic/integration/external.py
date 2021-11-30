"""Symbolic Integration via External Software

TESTS::

    sage: from sage.symbolic.integration.external import sympy_integrator
    sage: sympy_integrator(sin(x), x)
    -cos(x)
"""
from sage.symbolic.expression import Expression
from sage.symbolic.ring import SR


def maxima_integrator(expression, v, a=None, b=None):
    """
    Integration using Maxima

    EXAMPLES::

        sage: from sage.symbolic.integration.external import maxima_integrator
        sage: maxima_integrator(sin(x), x)
        -cos(x)
        sage: maxima_integrator(cos(x), x)
        sin(x)
        sage: f(x) = function('f')(x)
        sage: maxima_integrator(f(x), x)
        integrate(f(x), x)

    TESTS:

    Check that :trac:`25817` is fixed::

        sage: maxima_integrator(log(e^x*log(x)*sin(x))/x^2, x)
        1/2*(x*(Ei(-log(x)) + conjugate(Ei(-log(x))))
        - 2*x*integrate(sin(x)/(x*cos(x)^2 + x*sin(x)^2
        + 2*x*cos(x) + x), x) + 2*x*integrate(sin(x)/(x*cos(x)^2
        + x*sin(x)^2 - 2*x*cos(x) + x), x) + 2*x*log(x) + 2*log(2)
        - log(cos(x)^2 + sin(x)^2 + 2*cos(x) + 1) - log(cos(x)^2
        + sin(x)^2 - 2*cos(x) + 1) - 2*log(log(x)))/x
    """
    from sage.calculus.calculus import maxima
    if not isinstance(expression, Expression):
        expression = SR(expression)
    if a is None:
        result = maxima.sr_integral(expression, v)
    else:
        result = maxima.sr_integral(expression, v, a, b)
    return result._sage_()


def sympy_integrator(expression, v, a=None, b=None):
    """
    Integration using SymPy

    EXAMPLES::

        sage: from sage.symbolic.integration.external import sympy_integrator
        sage: sympy_integrator(sin(x), x)
        -cos(x)
        sage: sympy_integrator(cos(x), x)
        sin(x)
    """
    import sympy
    ex = expression._sympy_()
    v = v._sympy_()
    if a is None:
        result = sympy.integrate(ex, v)
    else:
        result = sympy.integrate(ex, (v, a._sympy_(), b._sympy_()))
    return result._sage_()


def mma_free_integrator(expression, v, a=None, b=None):
    """
    Integration using Mathematica's online integrator

    EXAMPLES::

        sage: from sage.symbolic.integration.external import mma_free_integrator
        sage: mma_free_integrator(sin(x), x) # optional - internet
        -cos(x)

    A definite integral::

        sage: mma_free_integrator(e^(-x), x, a=0, b=oo) # optional - internet
        1

    TESTS:

    Check that :trac:`18212` is resolved::

        sage: var('y')   # optional - internet
        y
        sage: result = integral(sin(y)^2, y, algorithm='mathematica_free') # optional - internet
        sage: result.simplify_trig()               # optional - internet
        -1/2*cos(y)*sin(y) + 1/2*y

    ::

    Check that :trac:`14764` is resolved::

        sage: integrate(x^2, x, 0, 1, algorithm="mathematica_free") # optional - internet
        1/3
        sage: integrate(sin(x), [x, 0, pi], algorithm="mathematica_free") # optional - internet
        2
        sage: integrate(sqrt(x), (x, 0, 1), algorithm="mathematica_free") # optional - internet
        2/3

    ::

        sage: mma_free_integrator(exp(-x^2)*log(x), x) # optional - internet
        1/2*sqrt(pi)*erf(x)*log(x) - x*hypergeometric((1/2, 1/2), (3/2, 3/2), -x^2)


    """
    from sage.interfaces.mathematica import request_wolfram_alpha, parse_moutput_from_json, symbolic_expression_from_mathematica_string
    math_expr = expression._mathematica_init_()
    variable = v._mathematica_init_()
    if a is None and b is None:
        input = "Integrate[{},{}]".format(math_expr, variable)
    elif a is not None and b is not None:
        input = "Integrate[{},{{{},{},{}}}]".format(math_expr, variable,
                    a._mathematica_init_(), b._mathematica_init_())
    else:
        raise ValueError('a(={}) and b(={}) should be both None'
                         ' or both defined'.format(a, b))
    json_page_data = request_wolfram_alpha(input)
    all_outputs = parse_moutput_from_json(json_page_data)
    if not all_outputs:
        raise ValueError("no outputs found in the answer from Wolfram Alpha")
    first_output = all_outputs[0]
    return symbolic_expression_from_mathematica_string(first_output)


def fricas_integrator(expression, v, a=None, b=None, noPole=True):
    """
    Integration using FriCAS

    EXAMPLES::

        sage: from sage.symbolic.integration.external import fricas_integrator  # optional - fricas
        sage: fricas_integrator(sin(x), x)                                      # optional - fricas
        -cos(x)
        sage: fricas_integrator(cos(x), x)                                      # optional - fricas
        sin(x)
        sage: fricas_integrator(1/(x^2-2), x, 0, 1)                             # optional - fricas
        -1/8*sqrt(2)*(log(2) - log(-24*sqrt(2) + 34))
        sage: fricas_integrator(1/(x^2+6), x, -oo, oo)                          # optional - fricas
        1/6*sqrt(6)*pi

    TESTS:

    Check that :trac:`25220` is fixed::

        sage: integral(sqrt(1-cos(x)), x, 0, 2*pi, algorithm="fricas")          # optional - fricas
        4*sqrt(2)

    Check that in case of failure one gets unevaluated integral::

        sage: integral(cos(ln(cos(x))), x, 0, pi/8, algorithm='fricas')         # optional - fricas
        integrate(cos(log(cos(x))), x, 0, 1/8*pi)

        sage: integral(cos(ln(cos(x))), x, algorithm='fricas')                  # optional - fricas
        integral(cos(log(cos(x))), x)

    Check that :trac:`28641` is fixed::

        sage: integrate(sqrt(2)*x^2 + 2*x, x, algorithm="fricas")               # optional - fricas
        1/3*sqrt(2)*x^3 + x^2

        sage: integrate(sqrt(2), x, algorithm="fricas")                         # optional - fricas
        sqrt(2)*x

        sage: integrate(1, x, algorithm="fricas")                               # optional - fricas
        x

    Check that :trac:`28630` is fixed::

        sage: f = polylog(3, x)
        sage: f.integral(x, algorithm='fricas')                                 # optional - fricas
        -x*dilog(x) - (x - 1)*log(-x + 1) + x*polylog(3, x) + x

    Check that :trac:`29043` is fixed::

        sage: var("a c d"); f = (I*a*tan(d*x + c) + a)*sec(d*x + c)^10
        (a, c, d)
        sage: ii = integrate(f, x, algorithm="fricas")                               # optional - fricas
        sage: 1/315*(64512*I*a*e^(10*I*d*x + 10*I*c) + 53760*I*a*e^(8*I*d*x + 8*I*c) + 30720*I*a*e^(6*I*d*x + 6*I*c) + 11520*I*a*e^(4*I*d*x + 4*I*c) + 2560*I*a*e^(2*I*d*x + 2*I*c) + 256*I*a)/(d*e^(20*I*d*x + 20*I*c) + 10*d*e^(18*I*d*x + 18*I*c) + 45*d*e^(16*I*d*x + 16*I*c) + 120*d*e^(14*I*d*x + 14*I*c) + 210*d*e^(12*I*d*x + 12*I*c) + 252*d*e^(10*I*d*x + 10*I*c) + 210*d*e^(8*I*d*x + 8*I*c) + 120*d*e^(6*I*d*x + 6*I*c) + 45*d*e^(4*I*d*x + 4*I*c) + 10*d*e^(2*I*d*x + 2*I*c) + d) - ii                        # optional - fricas
        0
    """
    if not isinstance(expression, Expression):
        expression = SR(expression)

    from sage.interfaces.fricas import fricas
    e_fricas = fricas(expression)
    v_fricas = fricas(v)

    if a is None:
        result = e_fricas.integrate(v_fricas)
    else:
        seg = fricas.equation(v_fricas, fricas.segment(a, b))

        if noPole:
            result = e_fricas.integrate(seg, '"noPole"')
        else:
            result = e_fricas.integrate(seg)

    result = result.sage()

    if result == "failed":
        result = expression.integrate(v, a, b, hold=True)

    elif result == "potentialPole":
        raise ValueError("The integrand has a potential pole"
                         " in the integration interval")

    return result


def giac_integrator(expression, v, a=None, b=None):
    r"""
    Integration using Giac

    EXAMPLES::

        sage: from sage.symbolic.integration.external import giac_integrator
        sage: giac_integrator(sin(x), x)
        -cos(x)
        sage: giac_integrator(1/(x^2+6), x, -oo, oo)
        1/6*sqrt(6)*pi

    TESTS::

        sage: giac_integrator(e^(-x^2)*log(x), x)
        integrate(e^(-x^2)*log(x), x)

    Check that :trac:`30133` is fixed::

        sage: ee = SR.var('e')
        sage: giac_integrator(ee^x, x)
        e^x/log(e)
        sage: y = SR.var('π')
        sage: giac_integrator(cos(y), y)
        sin(π)

    Check that :trac:`29966` is fixed::

        sage: giac_integrator(sqrt(x + sqrt(x)), x)
        1/12*(2*sqrt(x)*(4*sqrt(x) + 1) - 3)*sqrt(x + sqrt(x))...
    """
    ex = expression._giac_()
    if a is None:
        result = ex.integrate(v._giac_())
    else:
        result = ex.integrate(v._giac_(), a._giac_(), b._giac_())
    if 'integrate' in format(result) or 'integration' in format(result):
        return expression.integrate(v, a, b, hold=True)
    else:
        return result._sage_()

def libgiac_integrator(expression, v, a=None, b=None):
    r"""
    Integration using libgiac

    EXAMPLES::

        sage: import sage.libs.giac
        ...
        sage: from sage.symbolic.integration.external import libgiac_integrator
        sage: libgiac_integrator(sin(x), x)
        -cos(x)
        sage: libgiac_integrator(1/(x^2+6), x, -oo, oo)
        No checks were made for singular points of antiderivative...
        1/6*sqrt(6)*pi

    TESTS::

        sage: libgiac_integrator(e^(-x^2)*log(x), x)
        integrate(e^(-x^2)*log(x), x)

    The following integral fails with the Giac Pexpect interface, but works
    with libgiac (:trac:`31873`)::

        sage: a, x = var('a,x')
        sage: f = sec(2*a*x)
        sage: F = libgiac_integrator(f, x)
        ...
        sage: (F.derivative(x) - f).simplify_trig()
        0
    """
    from sage.libs.giac import libgiac
    from sage.libs.giac.giac import Pygen
    # We call Pygen on first argument because otherwise some expressions
    # involving derivatives result in doctest failures in interfaces/sympy.py
    # -- related to the fixme note in sage.libs.giac.giac.GiacFunction.__call__
    # regarding conversion of lists.
    if a is None:
        result = libgiac.integrate(Pygen(expression), v)
    else:
        result = libgiac.integrate(Pygen(expression), v, a, b)
    if 'integrate' in format(result) or 'integration' in format(result):
        return expression.integrate(v, a, b, hold=True)
    else:
        return result.sage()
