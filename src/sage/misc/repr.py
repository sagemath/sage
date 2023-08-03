"""
Repr formatting support
"""


def coeff_repr(c, is_latex=False):
    r"""
    String representing coefficients in a linear combination.

    INPUT:

    - ``c`` -- a coefficient (i.e., an element of a ring)

    OUTPUT:

    A string

    EXAMPLES::

        sage: from sage.misc.repr import coeff_repr
        sage: coeff_repr(QQ(1/2))
        '1/2'
        sage: coeff_repr(-x^2)
        '(-x^2)'
        sage: coeff_repr(QQ(1/2), is_latex=True)
        '\\frac{1}{2}'
        sage: coeff_repr(-x^2, is_latex=True)
        '\\left(-x^{2}\\right)'
    """
    if not is_latex:
        try:
            return c._coeff_repr()
        except AttributeError:
            pass
    if isinstance(c, (int, float)):
        return str(c)
    if is_latex and hasattr(c, '_latex_'):
        s = c._latex_()
    else:
        s = str(c).replace(' ', '')
    if s.find("+") != -1 or s.find("-") != -1:
        if is_latex:
            return "\\left(%s\\right)" % s
        else:
            return "(%s)" % s
    return s


def repr_lincomb(terms, is_latex=False, scalar_mult="*", strip_one=False,
                 repr_monomial=None, latex_scalar_mult=None):
    """
    Compute a string representation of a linear combination of some
    formal symbols.

    INPUT:

    - ``terms`` -- list of terms, as pairs (support, coefficient)
    - ``is_latex`` -- whether to produce latex (default: ``False``)
    - ``scalar_mult`` -- string representing the multiplication (default:``'*'``)
    - ``latex_scalar_mult`` -- latex string representing the multiplication
      (default: a space if ``scalar_mult`` is ``'*'``; otherwise ``scalar_mult``)
    - ``coeffs`` -- for backward compatibility

    OUTPUT:

    -  ``str`` - a string

    EXAMPLES::

        sage: repr_lincomb([('a',1), ('b',-2), ('c',3)])
        'a - 2*b + 3*c'
        sage: repr_lincomb([('a',0), ('b',-2), ('c',3)])
        '-2*b + 3*c'
        sage: repr_lincomb([('a',0), ('b',2), ('c',3)])
        '2*b + 3*c'
        sage: repr_lincomb([('a',1), ('b',0), ('c',3)])
        'a + 3*c'
        sage: repr_lincomb([('a',-1), ('b','2+3*x'), ('c',3)])
        '-a + (2+3*x)*b + 3*c'
        sage: repr_lincomb([('a', '1+x^2'), ('b', '2+3*x'), ('c', 3)])
        '(1+x^2)*a + (2+3*x)*b + 3*c'
        sage: repr_lincomb([('a', '1+x^2'), ('b', '-2+3*x'), ('c', 3)])
        '(1+x^2)*a + (-2+3*x)*b + 3*c'
        sage: repr_lincomb([('a', 1), ('b', -2), ('c', -3)])
        'a - 2*b - 3*c'
        sage: t = PolynomialRing(RationalField(),'t').gen()
        sage: repr_lincomb([('a', -t), ('s', t - 2), ('', t^2 + 2)])
        '-t*a + (t-2)*s + (t^2+2)'

    Examples for ``scalar_mult``::

        sage: repr_lincomb([('a',1), ('b',2), ('c',3)], scalar_mult='*')
        'a + 2*b + 3*c'
        sage: repr_lincomb([('a',2), ('b',0), ('c',-3)], scalar_mult='**')
        '2**a - 3**c'
        sage: repr_lincomb([('a',-1), ('b',2), ('c',3)], scalar_mult='**')
        '-a + 2**b + 3**c'

    Examples for ``scalar_mult`` and ``is_latex``::

        sage: repr_lincomb([('a',-1), ('b',2), ('c',3)], is_latex=True)
        '-a + 2 b + 3 c'
        sage: repr_lincomb([('a',-1), ('b',-1), ('c',3)], is_latex=True, scalar_mult='*')
        '-a - b + 3 c'
        sage: repr_lincomb([('a',-1), ('b',2), ('c',-3)], is_latex=True, scalar_mult='**')
        '-a + 2**b - 3**c'
        sage: repr_lincomb([('a',-2), ('b',-1), ('c',-3)], is_latex=True, latex_scalar_mult='*')
        '-2*a - b - 3*c'
        sage: repr_lincomb([('a',-2), ('b',-1), ('c',-3)], is_latex=True, latex_scalar_mult='')
        '-2a - b - 3c'

    Examples for ``strip_one``::

        sage: repr_lincomb([ ('a',1), (1,-2), ('3',3) ])
        'a - 2*1 + 3*3'
        sage: repr_lincomb([ ('a',-1), (1,1), ('3',3) ])
        '-a + 1 + 3*3'
        sage: repr_lincomb([ ('a',1), (1,-2), ('3',3) ], strip_one = True)
        'a - 2 + 3*3'
        sage: repr_lincomb([ ('a',-1), (1,1), ('3',3) ], strip_one = True)
        '-a + 1 + 3*3'
        sage: repr_lincomb([ ('a',1), (1,-1), ('3',3) ], strip_one = True)
        'a - 1 + 3*3'

    Examples for ``repr_monomial``::

        sage: repr_lincomb([('a',1), ('b',2), ('c',3)], repr_monomial = lambda s: s+"1")
        'a1 + 2*b1 + 3*c1'

    TESTS:

    Verify that :trac:`31672` is fixed::

        sage: alpha = var("alpha")
        sage: repr_lincomb([(x, alpha)], is_latex=True)
        '\\alpha x'
        sage: A.<psi> = PolynomialRing(QQ)
        sage: B.<t> = FreeAlgebra(A)
        sage: (psi * t)._latex_()
        '\\psi t'
    """
    # Setting scalar_mult: symbol used for scalar multiplication
    if is_latex:
        if latex_scalar_mult is not None:
            scalar_mult = latex_scalar_mult
        elif scalar_mult == "*":
            scalar_mult = " "

    if repr_monomial is None:
        if is_latex:

            def repr_monomial(monomial):
                return monomial._latex_() if hasattr(monomial, '_latex_') else str(monomial)
        else:
            repr_monomial = str

    s = ""
    first = True

    if scalar_mult is None:
        scalar_mult = "" if is_latex else "*"

    for (monomial, c) in terms:
        if c != 0:
            coeff = coeff_repr(c)
            negative = False
            if len(coeff) and coeff[0] == "-":
                negative = True
            try:
                if c < 0:
                    negative = True
            except (NotImplementedError, TypeError):
                # comparisons may not be implemented for some coefficients
                pass
            if negative:
                coeff = coeff_repr(-c, is_latex)
            else:
                coeff = coeff_repr(c, is_latex)
            if coeff == "1":
                coeff = ""
            if coeff != "0":
                if negative:
                    if first:
                        sign = "-"  # add trailing space?
                    else:
                        sign = " - "
                else:
                    if first:
                        sign = ""
                    else:
                        sign = " + "
                b = repr_monomial(monomial)
                if len(b):
                    if coeff != "":
                        if b == "1" and strip_one:
                            b = ""
                        else:
                            b = scalar_mult + b
                s += "%s%s%s" % (sign, coeff, b)
                first = False
    if first:
        return "0"
        # this can happen only if are only terms with coeff_repr(c) == "0"
    # elif s == "":
        # return "1"  # is empty string representation invalid?
    else:
        return s
