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


def request_wolfram_alpha(input, verbose=False):
    r"""
    Request Wolfram Alpha website.

    INPUT:

    - ``input`` -- string
    - ``verbose`` -- bool (default: ``False``)

    OUTPUT:

    json

    EXAMPLES::

        sage: from sage.symbolic.integration.external import request_wolfram_alpha
        sage: page_data = request_wolfram_alpha('integrate Sin[x]')      # optional internet
        sage: [str(a) for a in sorted(page_data.keys())]                 # optional internet
        ['queryresult']
        sage: [str(a) for a in sorted(page_data['queryresult'].keys())]  # optional internet
        ['datatypes',
         'encryptedEvaluatedExpression',
         'encryptedParsedExpression',
         'error',
         'host',
         'id',
         'numpods',
         'parsetimedout',
         'parsetiming',
         'pods',
         'recalculate',
         'related',
         'server',
         'sponsorCategories',
         'success',
         'timedout',
         'timedoutpods',
         'timing',
         'version']

    """
    # import compatible with py2 and py3
    from urllib.parse import urlencode
    from urllib.request import Request, build_opener, HTTPCookieProcessor
    import json
    from http.cookiejar import CookieJar

    # we need cookies for this...
    cj = CookieJar()
    opener = build_opener(HTTPCookieProcessor(cj))
    # build initial query for code
    req = Request("http://www.wolframalpha.com/input/api/v1/code")
    resp = opener.open(req)
    # the website returns JSON containing the code
    page_data = json.loads(resp.read().decode("utf-8"))
    if not ("code" in page_data.keys()):
        raise ValueError("Wolfram did not return a code")
    proxy_code = page_data['code']
    if verbose:
        print("Code: {}".format(proxy_code))
        print("Cookies: {}".format(cj))
    # now we can make a request
    # some parameters documented here:
    #   https://products.wolframalpha.com/api/documentation/#parameter-reference
    # the following are the parameters used by the website
    params = {
        'assumptionsversion': '2',
        'async': 'true',
        'banners': 'raw',
        'debuggingdata': 'false',
        'format': 'image,plaintext,imagemap,sound,minput,moutput',
        'formattimeout': '8',
        'input': input,
        'output': 'JSON',
        'parsetimeout': '5',
        'podinfosasync': 'true',
        'proxycode': proxy_code,
        'recalcscheme': 'parallel',
        'sbsdetails': 'true',
        'scantimeout': '0.5',
        'sponsorcategories': 'true',
        'statemethod': 'deploybutton',
        'storesubpodexprs': 'true'}
    # # we can also change some parameters
    # params = {
    #     'assumptionsversion': '2',
    #     'banners': 'raw',
    #     'format': 'minput,moutput',
    #     'formattimeout': '8',
    #     'input': input,
    #     'output': 'JSON',
    #     'parsetimeout': '5',
    #     'proxycode': proxy_code,
    #     'scantimeout': '0.5',
    #     'storesubpodexprs': 'true'
    # }
    params = urlencode(params)
    url = "https://www.wolframalpha.com/input/json.jsp?%s" % params
    req = Request(url)
    req.add_header('Referer', "https://www.wolframalpha.com/input/")  # seems important
    resp = opener.open(req)
    # the website returns JSON containing the code
    return json.loads(resp.read().decode("utf-8"))


def parse_moutput_from_json(page_data, verbose=False):
    r"""
    Return the list of outputs found in the json (with key ``'moutput'``)

    INPUT:

    - ``page_data`` -- json obtained from Wolfram Alpha
    - ``verbose`` -- bool (default: ``False``)

    OUTPUT:

    list of unicode strings

    EXAMPLES::

        sage: from sage.symbolic.integration.external import request_wolfram_alpha
        sage: from sage.symbolic.integration.external import parse_moutput_from_json
        sage: page_data = request_wolfram_alpha('integrate Sin[x]') # optional internet
        sage: parse_moutput_from_json(page_data)                    # optional internet
        [u'-Cos[x]']

    ::

        sage: page_data = request_wolfram_alpha('Sin[x]')           # optional internet
        sage: L = parse_moutput_from_json(page_data)                # optional internet
        sage: sorted(L)                                             # optional internet
        ['-Cos[x]', '{{x == 0}}', '{{x == Pi C[1], Element[C[1], Integers]}}']

    TESTS::

        sage: page_data = request_wolfram_alpha('Integrate(Sin[z], y)')  # optional internet
        sage: parse_moutput_from_json(page_data)                         # optional internet
        Traceback (most recent call last):
        ...
        ValueError: asking wolframalpha.com was not successful

    """
    queryresult = page_data['queryresult']
    if not queryresult['success']:
        raise ValueError('asking wolframalpha.com was not successful')
    if 'pods' not in queryresult:
        raise ValueError('json object contains no pods')
    pods = queryresult['pods']
    if verbose:
        print("  Query successful: {}".format(queryresult['success']))
        print("  Number of results: {}".format(len(pods)))
    L = []
    for i, result in enumerate(pods):
        if verbose:
            print("  Result #{}".format(i))
            print("    Title: {}".format(result['title']))
        if 'subpods' not in result:
            continue
        subpods = result[u'subpods']
        for j, subpod in enumerate(subpods):
            if verbose:
                print("    Subpod #{}".format(j))
                if 'minput' in subpod.keys():
                    print("      MInput: {}".format(subpod['minput']))
                if 'moutput' in subpod.keys():
                    print("      MOutput: {}".format(subpod['moutput']))
            if 'moutput' in subpod.keys():
                L.append(subpod['moutput'])
    return L


def symbolic_expression_from_mathematica_string(mexpr):
    r"""
    Translate a mathematica string into a symbolic expression

    INPUT:

    - ``mexpr`` -- string

    OUTPUT:

    symbolic expression

    EXAMPLES::

        sage: from sage.symbolic.integration.external import symbolic_expression_from_mathematica_string
        sage: symbolic_expression_from_mathematica_string(u'-Cos[x]')
        -cos(x)

    """
    import re
    from sage.libs.pynac.pynac import symbol_table
    from sage.interfaces.mathematica import _un_camel as un_camel
    from sage.symbolic.constants import constants_name_table as constants
    from sage.calculus.calculus import symbolic_expression_from_string
    from sage.calculus.calculus import _find_func as find_func

    expr = mexpr.replace('\n', ' ').replace('\r', '')
    expr = expr.replace('[', '(').replace(']', ')')
    expr = expr.replace('{', '[').replace('}', ']')
    lsymbols = symbol_table['mathematica'].copy()
    autotrans = [lambda x:x.lower(),      # Try it in lower case
                 un_camel,      # Convert `CamelCase` to `camel_case`
                 lambda x: x]     # Try the original name
    # Find the MMA funcs/vars/constants - they start with a letter.
    # Exclude exponents (e.g. 'e8' from 4.e8)
    p = re.compile(r'(?<!\.)[a-zA-Z]\w*')

    for m in p.finditer(expr):
        # If the function, variable or constant is already in the
        # translation dictionary, then just move on.
        if m.group() in lsymbols:
            pass
        # Now try to translate all other functions -- try each strategy
        # in `autotrans` and check if the function exists in Sage
        elif m.end() < len(expr) and expr[m.end()] == '(':
            for t in autotrans:
                f = find_func(t(m.group()), create_when_missing=False)
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
    return symbolic_expression_from_string(expr, lsymbols, accept_sequence=True)


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

    TESTS:

    Check that :trac:`25220` is fixed::

        sage: integral(sqrt(1-cos(x)), x, 0, 2*pi, algorithm="fricas")          # optional - fricas
        4*sqrt(2)

    Check that in case of failure one gets unevaluated integral::

        sage: integral(cos(ln(cos(x))), x, 0, pi/8, algorithm='fricas')   # optional - fricas
        integrate(cos(log(cos(x))), x, 0, 1/8*pi)
        sage: integral(cos(ln(cos(x))), x, algorithm='fricas')   # optional - fricas
        integral(cos(log(cos(x))), x)
    """
    if not isinstance(expression, Expression):
        expression = SR(expression)

    from sage.interfaces.fricas import fricas
    ex = fricas(expression)

    if a is None:
        result = ex.integrate(v)
    else:
        seg = fricas.equation(v, fricas.segment(a, b))

        if noPole:
            result = ex.integrate(seg, '"noPole"')
        else:
            result = ex.integrate(seg)

    result = result.sage()

    if result == "failed":
        return expression.integrate(v, a, b, hold=True)

    if result == "potentialPole":
        raise ValueError("The integrand has a potential pole"
                         " in the integration interval")

    return result


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
