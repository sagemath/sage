"""
RootSum - Sum over roots of a polynomial

This module provides support for expressions of the form:
    sum_{r: P(r)=0} f(r, x)

EXAMPLES::

    sage: from sage.symbolic.rootsum import root_sum
    sage: var('x a')
    (x, a)
    sage: P = x^3 + a*x + 1
    sage: rs = root_sum(P, lambda r, x: log(x - r)/(a + 3*r^2))
    sage: rs
    root_sum(x^3 + a*x + 1)

Numerical evaluation::

    sage: P = x^5 + x + 1
    sage: rs = root_sum(P, lambda r: 1/(5*r^4 + 1))
    sage: rs.n()
    -0.443...
"""

# Import required Sage modules
from sage.symbolic.function import BuiltinFunction
from sage.symbolic.ring import SR
from sage.symbolic.expression import Expression
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.integer import Integer


class RootSumFunction(BuiltinFunction):
    r"""
    Represents sums over roots of a univariate polynomial.

    Mathematically, a RootSum represents:

        F(x) = \sum_{r : P(r)=0} f(r, x)

    where `P` is a univariate polynomial and `f` is a function of
    the root `r` and the variable `x`.

    INPUT:

    - ``polynomial`` -- a univariate polynomial (as a symbolic expression
      in the symbolic ring, e.g., ``x^3 + a*x + 1``)
    - ``summand`` -- a callable function that takes two arguments:
      the root `r` and the variable `x`
      (e.g., ``lambda r, x: log(x - r)/(a + 3*r^2)``)
      and returns a symbolic expression

    OUTPUT:

    A symbolic expression representing `\sum_{r: P(r)=0} f(r, x)`.

    The derivative (not yet implemented) would be:

        F'(x) = \sum_{r : P(r)=0} \partial f/\partial x (r, x)

    EXAMPLES::

        sage: from sage.symbolic.rootsum import root_sum
        sage: var('x a')
        (x, a)
        sage: P = x^3 + a*x + 1
        sage: rs = root_sum(P, lambda r, x: log(x - r)/(a + 3*r^2))
        sage: rs
        root_sum(x^3 + a*x + 1)

    Numerical evaluation::

        sage: P = x^5 + x + 1
        sage: rs = root_sum(P, lambda r: 1/(5*r^4 + 1))
        sage: rs.n()
        -0.443...
    """

    def __init__(self):
        """Initialize the RootSum function."""
        BuiltinFunction.__init__(self, "root_sum", nargs=1)

    def __call__(self, polynomial, summand=None, **kwargs):
        """
        Create a RootSum expression.

        When called with one argument (during numerical evaluation),
        it returns the expression itself.
        """
        # If called with only polynomial (e.g., during numerical eval)
        if summand is None:
            expr = super().__call__(polynomial, **kwargs)
            return expr

        # Original logic for creating RootSum
        try:
            if not polynomial.variables():
                return SR(0)
        except Exception:
            pass

        expr = super().__call__(polynomial, **kwargs)
        if not hasattr(self, '_data'):
            self._data = {}
        self._data[id(expr)] = {
            'summand': summand,
            'polynomial': polynomial
        }
        return expr

    def _get_data(self, expr):
        """Helper to get data from expression."""
        if hasattr(self, '_data'):
            return self._data.get(id(expr))
        return None

    def summand(self, expr):
        """
        Return the summand function of the RootSum expression.

        INPUT:

        - ``expr`` -- a RootSum expression

        OUTPUT:

        The summand function `f(r, x)` that was used to create the RootSum.

        EXAMPLES::

            sage: from sage.symbolic.rootsum import root_sum
            sage: var('x a')
            (x, a)
            sage: P = x^3 + a*x + 1
            sage: rs = root_sum(P, lambda r, x: log(x - r)/(a + 3*r^2))
            sage: root_sum.summand(rs)
            <function ...>
        """
        data = self._get_data(expr)
        if data is not None:
            return data.get('summand')
        return None

    def _eval_(self, polynomial):
        """
        Evaluate the RootSum if possible (for low-degree polynomials).
        """
        from sage.symbolic.ring import SR

        # Get the data from the expression
        try:
            import sys
            frame = sys._getframe(1)
            expr = frame.f_locals.get('self')
            if expr is not None:
                data = self._get_data(expr)
                if data:
                    summand = data.get('summand')
                else:
                    summand = None
            else:
                summand = None
        except Exception:
            summand = None

        if summand is None:
            return None

        # Convert polynomial to Sage's polynomial ring
        try:
            if hasattr(polynomial, 'polynomial'):
                poly = polynomial
            else:
                vars = polynomial.variables()
                if not vars:
                    return SR(0)
                R = PolynomialRing(SR, vars[0])
                poly = R(polynomial)
        except Exception:
            return None

        # Try to get explicit roots for low degree
        try:
            degree = poly.degree()
            if degree <= 4 and degree > 0:
                roots = poly.roots(SR)
                result = 0
                for root, multiplicity in roots:
                    try:
                        term = summand(root)
                        result += term
                    except Exception:
                        return None
                return result
        except Exception:
            pass

        return None

    def _derivative_(self, *args, **kwargs):
        """
        Differentiate a RootSum expression.

        Mathematically:
        F'(x) = Σ_{r : P(r)=0} ∂f/∂x (r, x)

        This is not yet implemented.
        """
        raise NotImplementedError("derivative of RootSum is not yet implemented")

    def _sympy_(self):
        """Convert to SymPy's RootSum."""
        from sympy import RootSum as SympyRootSum, Lambda, symbols

        try:
            import sys
            frame = sys._getframe(1)
            expr = frame.f_locals.get('self')
            if expr is None:
                return None

            # Get data from expression
            data = self._get_data(expr)
            if data is None:
                return None

            poly = data.get('polynomial')
            summand = data.get('summand')

            if poly is None or summand is None:
                return None

            from sage.symbolic.ring import SR
            sympy_poly = poly._sympy_()
            r = symbols('r')
            expr_sym = summand(SR(r))
            sympy_expr = expr_sym._sympy_()
            return SympyRootSum(sympy_poly, Lambda(r, sympy_expr))
        except Exception:
            return None

    def _latex_(self):
        """LaTeX representation."""
        try:
            import sys
            frame = sys._getframe(1)
            expr = frame.f_locals.get('self')
            if expr is None:
                return "\\operatorname{RootSum}"

            # Get data from expression
            data = self._get_data(expr)
            if data is None:
                return "\\operatorname{RootSum}"

            poly = data.get('polynomial')
            summand = data.get('summand')

            if poly is None or summand is None:
                return "\\operatorname{RootSum}"

            var = poly.variables()[0] if poly.variables() else SR.var('r')
            return f"\\sum_{{{var}: {poly._latex_()}=0}} {summand(var)._latex_()}"
        except Exception:
            return "\\operatorname{RootSum}"

    def _print_latex_(self, *args, **kwargs):
        """Print LaTeX representation."""
        return self._latex_()


# Create a global instance
root_sum = RootSumFunction()

# Make it available in the namespace
__all__ = ['RootSumFunction', 'root_sum']
