r"""
RootSum -- sums over the roots of a polynomial.

A :class:`RootSumFunction` represents a symbolic expression of the form

.. MATH::

    F(x) = \sum_{r : P(r) = 0} f(r, x)

where `P` is a univariate polynomial, the sum runs over all roots `r`
of `P` (counted with multiplicity), and `f` is a symbolic expression
that may involve `r` and any number of free variables.

Such expressions arise, for example, as antiderivatives of rational
functions whose denominator does not factor over the coefficient field.
They are produced by Maxima, FriCAS, Wolfram and SymPy.  The classical
example (see also R. Fateman, "Simplifying RootSum Expressions") is

.. MATH::

    \int \frac{dx}{x^3 + a x + 1}
    = \sum_{r : r^3 + a r + 1 = 0} \frac{\log(x - r)}{a + 3 r^2}.

The root variable is bound.  Substitution acts on free occurrences only,
and alpha-renames the root variable when needed to avoid capture.  Bound
names are normalized at construction, so alpha-equivalent root sums have
identical symbolic representations.  As with other symbolic expressions,
``variables()`` lists all occurring symbols; use ``free_variables()`` to
obtain the unbound ones.

INPUT:

- ``poly`` -- a univariate polynomial, given as a symbolic expression,
  for example ``r^3 + a*r + 1``
- ``root_var`` -- the symbolic variable that is summed over; it must
  occur in ``poly``
- ``summand`` -- a symbolic expression in ``root_var`` (and possibly
  other free variables) giving the term to sum

OUTPUT:

A symbolic expression representing `\sum_{r:P(r)=0} f(r)`.

EXAMPLES::

    sage: from sage.symbolic.rootsum import root_sum
    sage: var('x a r')
    (x, a, r)
    sage: P = r^3 + a*r + 1
    sage: rs = root_sum(P, r, log(x - r)/(a + 3*r^2))
    sage: rs
    root_sum(r^3 + a*r + 1, r, log(-r + x)/(3*r^2 + a))

The three operands are accessible directly::

    sage: rs.operands()[0]
    r^3 + a*r + 1
    sage: rs.operands()[1]
    r
    sage: rs.operands()[2]
    log(-r + x)/(3*r^2 + a)

Substituting a free variable preserves the expression::

    sage: rs.subs(a == 0)
    root_sum(r^3 + 1, r, 1/3*log(-r + x)/r^2)

A low-degree polynomial evaluates explicitly when its roots are available::

    sage: root_sum(r^2 - 1, r, r^2)
    2

A quintic stays unevaluated::

    sage: root_sum(r^5 - r + 1, r, sin(r))
    root_sum(r^5 - r + 1, r, sin(r))

The derivative of a root sum is again a root sum, taken under the sum::

    sage: rs.diff(x)
    1/(x^3 + a*x + 1)

The motivating application is antidifferentiation of rational
functions.  The derivative of the returned root sum equals the
integrand::

    sage: # needs sympy
    sage: F = integrate(1/(x^3 + a*x + 1), x, algorithm='sympy')
    sage: D = F.diff(x)
    sage: val = root_sum.evaluate(D.subs(a == 1, x == 1), prec=100)
    sage: abs(val - 1/3) < 1e-25
    True
"""

from sage.symbolic.function import BuiltinFunction
from sage.symbolic.ring import SR
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing


def _canonical_root_variable(variables):
    """
    Choose a canonical name without free-variable or assumption collisions.

    Binding must not import the domain or assumptions of an unrelated
    user variable with the same printed name.

    EXAMPLES::

        sage: from sage.symbolic.rootsum import _canonical_root_variable
        sage: r = var('r')
        sage: _canonical_root_variable([r])
        r1
    """
    from sage.symbolic.assumptions import assumptions

    blocked = {str(v) for v in variables}
    for fact in assumptions():
        if hasattr(fact, 'variables'):
            blocked.update(str(v) for v in fact.variables())
        elif hasattr(fact, '_var'):
            blocked.add(str(fact._var))
    index = 0
    while True:
        name = 'r' if index == 0 else 'r%d' % index
        index += 1
        if name in blocked:
            continue
        variable = SR.var(name)
        if not variable.is_real():
            return variable


def _rational_root_sum(p, root_var, summand):
    """
    Return a rational root sum using modular inversion and Newton identities.

    Return ``None`` when the summand is not rational in ``root_var`` or
    its denominator is not invertible modulo ``p``.

    EXAMPLES::

        sage: from sage.symbolic.rootsum import _rational_root_sum
        sage: r = var('r')
        sage: R = PolynomialRing(SR, r)
        sage: p = (r^5-r+1).polynomial(ring=R)
        sage: _rational_root_sum(p, r, r^2)
        0
        sage: _rational_root_sum(p, r, sin(r)) is None
        True
        sage: p = ((r-1)^2).polynomial(ring=R)
        sage: _rational_root_sum(p, r, 1/(r-1)) is None
        True
    """
    numerator, denominator = summand.numerator_denominator()
    if not (numerator.is_polynomial(root_var)
            and denominator.is_polynomial(root_var)):
        return None
    R = p.parent()
    numerator = numerator.polynomial(ring=R)
    denominator = denominator.polynomial(ring=R)

    def normalize(q):
        return R([SR(c).normalize() for c in q.list()])

    # Extended Euclid, normalizing coefficients at each step to control
    # expression growth over the symbolic coefficient field.
    previous, current = p, denominator
    previous_inverse, inverse = R(0), R(1)
    while current:
        quotient, remainder = previous.quo_rem(current)
        previous, current = current, normalize(remainder)
        previous_inverse, inverse = inverse, normalize(
            previous_inverse - quotient*inverse)
    if previous.degree() != 0:
        return None
    inverse = previous_inverse / previous[0]
    reduced = normalize((numerator*inverse) % p)
    degree = p.degree()
    monic = normalize(p / p.leading_coefficient())
    powers = [SR(degree)]
    for k in range(1, reduced.degree()+1):
        powers.append((-sum(monic[degree-j]*powers[k-j]
                             for j in range(1, k))
                       - k*monic[degree-k]).normalize())
    return sum((c*powers[k] for k, c in enumerate(reduced)), SR(0)).normalize()


class RootSumFunction(BuiltinFunction):
    r"""
    The symbolic function ``root_sum``.

    See the module docstring for details.
    """

    def __call__(self, poly, root_var, summand, **kwds):
        """
        Validate and normalize the bound variable, including held sums.

        EXAMPLES::

            sage: r, s, x = var('r s x')
            sage: bool(root_sum(r^5-r+1, r, sin(r)) == root_sum(s^5-s+1, s, sin(s)))
            True
            sage: root_sum(r^5-r+1, 0, sin(r))
            Traceback (most recent call last):
            ...
            ValueError: root_var must be a symbolic variable
            sage: root_sum(0, r, sin(r))
            Traceback (most recent call last):
            ...
            ValueError: the zero polynomial does not define a finite root sum
        """
        poly, root_var, summand = map(SR, (poly, root_var, summand))
        if not root_var.is_symbol():
            raise ValueError("root_var must be a symbolic variable")
        if poly == 0:
            raise ValueError("the zero polynomial does not define a finite root sum")
        if not poly.is_polynomial(root_var):
            raise ValueError("poly must be polynomial in root_var")
        free = set(poly.free_variables()) | set(summand.free_variables())
        free.discard(root_var)
        canonical = _canonical_root_variable(free)
        if canonical != root_var:
            poly = poly.subs({root_var: canonical})
            summand = summand.subs({root_var: canonical})
        return BuiltinFunction.__call__(self, poly, canonical, summand, **kwds)

    def __init__(self):
        r"""
        Initialize the ``root_sum`` function.

        TESTS::

            sage: from sage.symbolic.rootsum import root_sum
            sage: root_sum
            root_sum
        """
        BuiltinFunction.__init__(self, "root_sum", nargs=3,
                                 evalf_params_first=False)

    def _eval_(self, poly, root_var, summand):
        r"""
        Simplify rational sums using modular inversion and Newton sums.

        No explicit algebraic roots are computed.  Nonrational sums are
        retained except for linear and quadratic defining polynomials.
        A denominator sharing a root with the defining polynomial is not
        inverted.  Identities with symbolic coefficients hold generically,
        where degree and denominator invertibility are preserved.

        EXAMPLES::

            sage: r, a, x = var('r a x')
            sage: root_sum(r^5-r+1, r, r^2)
            0
            sage: root_sum(r^5-r+1, r, 1/r)
            1
            sage: root_sum((r-1)^5, r, r)
            5
            sage: root_sum(r^3+a*r+1, r, 1/((3*r^2+a)*(x-r)))
            1/(x^3 + a*x + 1)
            sage: root_sum(r^5-r+1, r, sin(r))
            root_sum(r^5 - r + 1, r, sin(r))
        """
        R = PolynomialRing(SR, root_var)
        p = SR(poly).polynomial(ring=R)
        if p.is_zero():
            raise ValueError("the zero polynomial does not define a finite root sum")
        if p.degree() == 0:
            return SR(0)
        value = _rational_root_sum(p, root_var, SR(summand))
        if value is not None:
            return value
        if p.degree() > 2:
            return None
        try:
            roots = p.roots(SR)
        except (TypeError, NotImplementedError):
            return None
        if sum(m for _, m in roots) != p.degree():
            return None
        return sum((mult * summand.subs({root_var: r})
                    for r, mult in roots), SR(0))

    def _subs_(self, subs_map, options, poly, root_var, summand):
        """
        Substitute free occurrences, avoiding capture of the root variable.

        EXAMPLES::

            sage: r, x = var('r x')
            sage: rs = root_sum(r^5-r+1, r, log(x-r))
            sage: bool(rs.subs(r=0) == rs)
            True
            sage: rs.subs(x=r)
            root_sum(r1^5 - r1 + 1, r1, log(r - r1))
            sage: bool((r + rs).subs(r=0) == rs)
            True
            sage: sin(rs).subs(x=r)
            sin(root_sum(r1^5 - r1 + 1, r1, log(r - r1)))
            sage: rs.subs({rs: 7})
            7
        """
        original = BuiltinFunction.__call__(
            self, poly, root_var, summand, hold=True)
        # Whole-expression replacement must remain available even though
        # substitution into the bound operands is handled specially.
        for key, value in subs_map.items():
            if key.is_trivially_equal(original):
                return value
        with SR.temp_var() as dummy:
            rename = {root_var: dummy}
            new_poly = subs_map.apply_to(poly.subs(rename), options)
            new_summand = subs_map.apply_to(summand.subs(rename), options)
            free = set(new_poly.variables()) | set(new_summand.variables())
            free.discard(dummy)
            new_root = _canonical_root_variable(free)
            restore = {dummy: new_root}
            return self(new_poly.subs(restore), new_root,
                        new_summand.subs(restore))

    def _tderivative_(self, poly, root_var, summand, *args, **kwargs):
        r"""
        Differentiate through roots using implicit differentiation.

        For simple roots, ``dr/dx = -P_x/P_r``.  Repeated factors are
        squarefree-decomposed before applying this formula, and their
        multiplicities are retained.  The result is valid locally away
        from collisions, degree drops, and poles of the summand.

        EXAMPLES::

            sage: r, a, x = var('r a x')
            sage: F = root_sum(r^5+a*r+1, r, log(x-r)/(5*r^4+a))
            sage: F.diff(x)
            1/(x^5 + a*x + 1)
            sage: root_sum((r^3+a*r+1)^2, r, sin(r)).diff(a)
            2*root_sum(r^3 + a*r + 1, r, -r*cos(r)/(3*r^2 + a))
            sage: F.diff(r)
            0
        """
        variable = kwargs['diff_param']
        if variable == root_var:
            return SR(0)
        partial = SR(summand).diff(variable)
        if not SR(poly).has(variable):
            return root_sum(poly, root_var, partial)
        R = PolynomialRing(SR, root_var)
        p = SR(poly).polynomial(ring=R)
        factors = p.squarefree_decomposition()
        result = SR(0)
        for factor, multiplicity in factors:
            factor = R([SR(c).normalize() for c in factor.list()])
            q = factor(root_var)
            velocity = -q.diff(variable) / q.diff(root_var)
            result += multiplicity * root_sum(
                q, root_var, partial + SR(summand).diff(root_var)*velocity)
        return result

    def evaluate(self, expr, prec=53):
        r"""
        Numerically evaluate an expression containing root sums.

        Exact squarefree decomposition preserves multiplicities before
        numerical root finding.  ``prec`` is working bit precision; it
        does not certify the relative error of a nearly cancelling sum.

        EXAMPLES::

            sage: r = var('r')
            sage: rs = root_sum(r^5-r+1, r, sin(r))
            sage: root_sum.evaluate(rs, prec=100)
            -0.041691470337868009879662475806
            sage: root_sum.evaluate(root_sum((r-1)^5, r, r))
            5.00000000000000
            sage: root_sum.evaluate(1 + rs, prec=100)
            0.95830852966213199012033752419
        """
        from sage.rings.complex_mpfr import ComplexField
        from sage.symbolic.expression_conversions import ExpressionTreeWalker

        CF = ComplexField(prec)
        function = self

        class EvaluateRootSums(ExpressionTreeWalker):
            def composition(self, ex, operator):
                if operator is function:
                    poly, root_var, summand = ex.operands()
                    R = PolynomialRing(SR, root_var)
                    p = poly.polynomial(ring=R)
                    if p.is_zero():
                        raise ValueError("the zero polynomial does not define a finite root sum")
                    result = CF(0)
                    for factor, mult in p.squarefree_decomposition():
                        roots = factor.change_ring(CF).roots(multiplicities=False)
                        for root in roots:
                            term = self(summand.subs({root_var: root}))
                            result += mult * CF(term)
                    return SR(result)
                return super().composition(ex, operator)

        return CF(EvaluateRootSums(SR(expr))())

    def _evalf_(self, poly, root_var, summand, parent=None, algorithm=None):
        """
        Evaluate numerically without first evaluating the bound operands.

        EXAMPLES::

            sage: r = var('r')
            sage: rs = root_sum(r^5-r+1, r, sin(r))
            sage: abs(rs.n(100) - root_sum.evaluate(rs, prec=100)) < 1e-28
            True
            sage: abs((1 + rs).n(100) - root_sum.evaluate(1 + rs, prec=100)) < 1e-28
            True
        """
        precision = parent.precision() if parent is not None else 53
        expr = BuiltinFunction.__call__(self, poly, root_var, summand, hold=True)
        return self.evaluate(expr, prec=precision)

    def _print_latex_(self, poly, root_var, summand):
        r"""
        Return a LaTeX representation.

        EXAMPLES::

            sage: from sage.symbolic.rootsum import root_sum
            sage: var('x a r')
            (x, a, r)
            sage: latex(root_sum(r^3 + a*r + 1, r, log(x - r)/(a + 3*r^2)))
            \sum_{r : r^{3} + a r + 1 = 0} \frac{\log\left(-r + x\right)}{3 \, r^{2} + a}
        """
        from sage.misc.latex import latex
        return (r"\sum_{%s : %s = 0} %s"
                % (latex(root_var), latex(poly), latex(summand)))

    def _sympy_(self, poly, root_var, summand):
        r"""
        Convert to a SymPy ``RootSum``.

        EXAMPLES::

            sage: from sage.symbolic.rootsum import root_sum     # needs sympy
            sage: var('r')                                        # needs sympy
            r
            sage: root_sum(r^5 - r + 1, r, sin(r))._sympy_()      # needs sympy
            RootSum(r**5 - r + 1, Lambda(r, sin(r)))

        Coefficients may contain free parameters::

            sage: a = var('a')
            sage: rs = root_sum(r^5+a*r+1, r, sin(r))
            sage: rs._sympy_()._sage_().is_trivially_equal(rs)
            True
        """
        import sympy
        variable = root_var._sympy_()
        return sympy.RootSum(
            sympy.Poly(poly._sympy_(), variable),
            sympy.Lambda(variable, summand._sympy_()))


root_sum = RootSumFunction()

__all__ = ['RootSumFunction', 'root_sum']
