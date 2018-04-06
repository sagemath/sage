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

    TESTS:

    Check that :trac:`18212` is resolved::

        sage: var('y')   # optional - internet
        y
        sage: integral(sin(y)^2, y, algorithm='mathematica_free') # optional - internet
        -1/2*cos(y)*sin(y) + 1/2*y

        sage: mma_free_integrator(exp(-x^2)*log(x), x) # optional - internet
        1/2*sqrt(pi)*erf(x)*log(x) - x*hypergeometric((1/2, 1/2), (3/2, 3/2), -x^2)
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
        expression = expression.subs({x:shadow_x}).subs({v: x})
    params = urlencode({'expr': expression._mathematica_init_(), 'random': 'false'})
    page = urlopen("http://integrals.wolfram.com/home.jsp", params).read()
    page = page[page.index('"inputForm"'):page.index('"outputForm"')]
    page = re.sub("\s", "", page)
    mexpr = re.match(r".*Integrate.*==</em><br/>(.*)</p>", page).groups()[0]
    try:
        from sage.libs.pynac.pynac import symbol_table
        from sage.interfaces.mathematica import _un_camel as un_camel
        from sage.symbolic.constants import constants_name_table as constants
        from sage.calculus.calculus import symbolic_expression_from_string
        from sage.calculus.calculus import _find_func as find_func

        expr = mexpr.replace('\n',' ').replace('\r', '')
        expr = expr.replace('[', '(').replace(']', ')')
        expr = expr.replace('{', '[').replace('}', ']')
        lsymbols = symbol_table['mathematica'].copy()
        autotrans = [str.lower,      # Try it in lower case
                     un_camel,      # Convert `CamelCase` to `camel_case`
                     lambda x: x     # Try the original name
                    ]
        # Find the MMA funcs/vars/constants - they start with a letter.
        # Exclude exponents (e.g. 'e8' from 4.e8)
        p = re.compile('(?<!\.)[a-zA-Z]\w*')

        for m in p.finditer(expr):
            # If the function, variable or constant is already in the
            # translation dictionary, then just move on.
            if m.group() in lsymbols:
                pass
            # Now try to translate all other functions -- try each strategy
            # in `autotrans` and check if the function exists in Sage
            elif m.end() < len(expr) and expr[m.end()] == '(':
                for t in autotrans:
                    f = find_func(t(m.group()), create_when_missing = False)
                    if f is not None:
                        lsymbols[m.group()] = f
                        break
                else:
                    raise NotImplementedError("Don't know a Sage equivalent for Mathematica function '%s'." % m.group())
            # Check if Sage has an equivalent constant
            else:
                for t in autotrans:
                    if t(m.group()) in constants:
                        lsymbols[m.group()] = constants[t(m.group())]
                        break
        ans = symbolic_expression_from_string(expr, lsymbols, accept_sequence=True)
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

def giac_integrator(expression, v, a=None, b=None):
    """
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
