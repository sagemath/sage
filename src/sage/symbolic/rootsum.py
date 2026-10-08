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

.. WARNING::

    ``root_var`` is a bound variable, but substitution does not respect
    the binding.  In particular ``rs.subs(root_var == 0)`` substitutes
    into the polynomial and the summand.  Renaming and free-variable
    substitution behave as expected.  This matches the current behavior
    of :class:`~sage.calculus.calculus.Sum` and
    :class:`~sage.calculus.calculus.Integral`.  For example::

        sage: var('r x')
        (r, x)
        sage: rs = root_sum(r^5 - r + 1, r, log(x - r))
        sage: rs.variables()
        (r, x)
        sage: rs.subs(r == 0)
        root_sum(1, 0, log(x))

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
    root_sum(r^3 + a*r + 1, r, -1/((3*r^2 + a)*(r - x)))

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


class RootSumFunction(BuiltinFunction):
    r"""
    The symbolic function ``root_sum``.

    See the module docstring for details.
    """

    def __init__(self):
        r"""
        Initialize the ``root_sum`` function.

        TESTS::

            sage: from sage.symbolic.rootsum import root_sum
            sage: root_sum
            root_sum
        """
        BuiltinFunction.__init__(self, "root_sum", nargs=3)

    def _eval_(self, poly, root_var, summand):
        r"""
        Evaluate the sum explicitly when the roots are available.

        Returns ``None`` (leave unevaluated) unless the polynomial has
        degree at most 2 and all of its roots lie in the symbolic ring.

        EXAMPLES::

            sage: from sage.symbolic.rootsum import root_sum
            sage: var('r a')
            (r, a)
            sage: root_sum(r^2 - 1, r, r^2)
            2
            sage: root_sum(r^5 - r + 1, r, r^2)
            root_sum(r^5 - r + 1, r, r^2)
            sage: root_sum(r^3 + a*r + 1, r, r)
            root_sum(r^3 + a*r + 1, r, r)
            sage: root_sum(r^2 - 1, r, 1)
            2
            sage: root_sum(r - 1, r, r)
            1
        """
        try:
            R = PolynomialRing(SR, root_var)
            p = SR(poly).polynomial(ring=R)
        except (TypeError, ValueError, AttributeError):
            return None
        if p.is_zero():
            return None
        if p.degree() == 0:
            # Degree 0: either a genuine constant (empty sum = 0)
            # or an expression that is not polynomial in root_var
            # (leave unevaluated).
            if SR(poly).is_constant():
                return SR(0)
            return None
        if p.degree() > 2:
            return None
        try:
            roots = p.roots(SR)
        except (TypeError, NotImplementedError):
            return None
        if sum(m for _, m in roots) != p.degree():
            return None
        result = SR(0)
        for r, mult in roots:
            result += mult * summand.subs({root_var: r})
        return result

    def _tderivative_(self, poly, root_var, summand, *args, **kwargs):
        r"""
        Differentiate a ``root_sum`` expression under the sum.

        If the differentiation variable occurs only in the summand, the
        derivative acts under the sum:

        .. MATH::

            \frac{\partial}{\partial x}
            \sum_{r : P(r) = 0} f(r, x)
            = \sum_{r : P(r) = 0} \frac{\partial f}{\partial x}(r, x).

        When the differentiation variable occurs in the polynomial, the
        roots themselves depend on it, so the derivative involves an
        implicit-differentiation term and is not implemented here;
        a :class:`NotImplementedError` is raised.

        EXAMPLES::

            sage: from sage.symbolic.rootsum import root_sum
            sage: var('x a r')
            (x, a, r)
            sage: rs = root_sum(r^3 + a*r + 1, r, log(x - r)/(a + 3*r^2))
            sage: rs.diff(x)
            root_sum(r^3 + a*r + 1, r, -1/((3*r^2 + a)*(r - x)))
            sage: rs.diff(r)
            0
            sage: root_sum(r^5 - r + 1, r, 1).diff(x)
            0

        Differentiating with respect to a parameter that occurs in the
        polynomial is not supported::

            sage: rs.diff(a)
            Traceback (most recent call last):
            ...
            NotImplementedError: derivative of root_sum with respect to a
            parameter occurring in the polynomial is not implemented
        """
        diff_param = kwargs.get('diff_param')
        if diff_param is None:
            return None
        if diff_param == root_var:
            return SR(0)
        if SR(poly).has(diff_param):
            raise NotImplementedError(
                "derivative of root_sum with respect to a parameter "
                "occurring in the polynomial is not implemented")
        if not hasattr(summand, 'diff'):
            # Summand is a constant (e.g. Integer(1)).
            return SR(0)
        new_summand = summand.diff(diff_param)
        if new_summand == 0:
            return SR(0)
        return root_sum(poly, root_var, new_summand)

    def evaluate(self, expr, prec=53):
        r"""
        Numerically evaluate a ``root_sum`` expression.

        The roots of the polynomial are computed numerically in the
        complex field of precision ``prec`` and the summand is summed
        over them, preserving multiplicities.

        INPUT:

        - ``expr`` -- a ``root_sum`` expression
        - ``prec`` -- bit precision (default: 53)

        EXAMPLES::

            sage: from sage.symbolic.rootsum import root_sum
            sage: var('r')
            r
            sage: rs = root_sum(r^5 - r + 1, r, sin(r))
            sage: root_sum.evaluate(rs, prec=100)
            -0.041691470337868009879662475806

        Multiplicities are respected::

            sage: root_sum.evaluate(root_sum((r - 1)^5, r, r))
            5.00000000000000
        """
        from sage.rings.complex_mpfr import ComplexField
        ops = expr.operands()
        if len(ops) != 3:
            return None
        poly, root_var, summand = ops
        try:
            R = PolynomialRing(SR, root_var)
            p = SR(poly).polynomial(ring=R)
        except (TypeError, ValueError, AttributeError):
            return None
        if p.degree() == 0 and not SR(poly).is_constant():
            return None

        CF = ComplexField(prec)

        try:
            sfd = p.squarefree_decomposition()
        except (TypeError, NotImplementedError):
            sfd = None

        result = CF(0)
        if sfd is not None:
            for factor, mult in sfd:
                factor_cc = factor.change_ring(CF)
                roots = factor_cc.roots(multiplicities=False)
                for r in roots:
                    result += mult * CF(summand.subs({root_var: r}))
            return result

        try:
            roots = p.roots(ring=CF, multiplicities=False)
        except (TypeError, NotImplementedError):
            return None
        if not roots:
            return None
        for r in roots:
            result += CF(summand.subs({root_var: r}))
        return result

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
        """
        import sympy
        return sympy.RootSum(
            poly._sympy_(),
            sympy.Lambda(root_var._sympy_(), summand._sympy_()))


root_sum = RootSumFunction()

__all__ = ['RootSumFunction', 'root_sum']
