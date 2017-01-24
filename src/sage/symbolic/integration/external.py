"Symbolic integration via external software"

from sage.symbolic.expression import Expression
from sage.symbolic.ring import SR


def maxima_integrator(expression, v, a=None, b=None):
    """
        sage: from sage.symbolic.integration.external import maxima_integrator
        sage: maxima_integrator(sin(x), x)
        -cos(x)
        sage: maxima_integrator(cos(x), x)
        sin(x)
        sage: f(x) = function('f')(x)
        sage: maxima_integrator(f(x), x)
        integrate(f(x), x)
    """
    from sage.calculus.calculus import maxima
    if not isinstance(expression, Expression):
        expression = SR(expression)
    if a is None:
        result = maxima.sr_integral(expression,v)
    else:
        result = maxima.sr_integral(expression, v, a, b)
    return result._sage_()

def sympy_integrator(expression, v, a=None, b=None):
    """
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
        sage: from sage.symbolic.integration.external import mma_free_integrator
        sage: mma_free_integrator(sin(x), x) # optional - internet
        -cos(x)
    """
    import re
    # import compatible with py2 and py3
    from six.moves.urllib.request import urlopen
    from six.moves.urllib.parse import urlencode
    # We need to integrate against x
    vars = [str(x) for x in expression.variables()]
    if any(len(x)>1 for x in vars):
        raise NotImplementedError("Mathematica online integrator can only handle single letter variables.")
    x = SR.var('x')
    if repr(v) != 'x':
        for i in range(ord('a'), ord('z')+1):
            if chr(i) not in vars:
                shadow_x = SR.var(chr(i))
                break
        expression = expression.subs({x:shadow_x}).subs({dvar: x})
    params = urlencode({'expr': expression._mathematica_init_(), 'random': 'false'})
    page = urlopen("http://integrals.wolfram.com/index.jsp", params).read()
    page = page[page.index('"inputForm"'):page.index('"outputForm"')]
    page = re.sub("\s", "", page)
    mexpr = re.match(r".*Integrate.*==</em><br/>(.*)</p>", page).groups()[0]
    try:
        ans = SR(mexpr.lower().replace('[', '(').replace(']', ')'))
        if repr(v) != 'x':
            ans = ans.subs({x:v}).subs({shadow_x:x})
        return ans
    except TypeError:
        raise ValueError("Unable to parse: %s" % mexpr)


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
        1/4*sqrt(2)*(log(3*sqrt(2) - 4) - log(sqrt(2)))
        sage: fricas_integrator(1/(x^2+6), x, -oo, oo)                          # optional - fricas
        1/6*sqrt(6)*pi
    """
    if not isinstance(expression, Expression):
        expression = SR(expression)
    if a is None:
        result = expression._fricas_().integrate(v)
    else:
        import sage.rings.infinity
        if a == sage.rings.infinity.PlusInfinity():
            a = "%plusInfinity"
        elif a == sage.rings.infinity.MinusInfinity():
            a = "%minusInfinity"
        if b == sage.rings.infinity.PlusInfinity():
            b = "%plusInfinity"
        elif b == sage.rings.infinity.MinusInfinity():
            b = "%minusInfinity"

        if noPole:
            result = expression._fricas_().integrate("{}={}..{}".format(v, a, b), '"noPole"')
        else:
            result = expression._fricas_().integrate("{}={}..{}".format(v, a, b))

    locals = {str(v): v for v in expression.variables()}
    if str(result) == "potentialPole":
        raise ValueError("The integrand has a potential pole"
                         " in the integration interval")

    return result.sage()
