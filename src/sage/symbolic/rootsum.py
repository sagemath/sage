r"""
Sum over the roots of a polynomial

This module implements :class:`Function_root_sum`, a symbolic function
representing a sum

.. MATH::

    \sum_{r \,:\, p(r) = 0} f(r)

over the roots of a polynomial ``p`` in a dummy variable ``r``.  This
matches Maxima's ``lsum(f(r), r, rootsof(p, r))`` representation, which
Maxima emits when ``integrate_use_rootsof:true`` is set, and SymPy's
``RootSum(p, Lambda(r, f(r)))``.

Concretely, the antiderivative

.. MATH::

    \int \frac{1}{x^3 + a x + 1} \, dx
        = \sum_{r \,:\, r^3 + a r + 1 = 0} \frac{\log(x - r)}{3 r^2 + a}

is represented in Sage as
``root_sum(log(x - r)/(3*r^2 + a), r, r^3 + a*r + 1)``.

See :issue:`40356` for the original feature request and the converter in
:mod:`sage.interfaces.maxima_lib`.
"""
# ****************************************************************************
#       Copyright (C) 2026
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.symbolic.function import BuiltinFunction


class Function_root_sum(BuiltinFunction):
    r"""
    Symbolic representation of `\sum_{r \,:\, p(r) = 0} f(r)`, the sum of
    an expression ``f(r)`` over the roots of a polynomial ``p`` in a
    dummy variable ``r``.

    The three arguments are ``(body, dummy, poly)``, matching Maxima's
    ``lsum(body, dummy, rootsof(poly, dummy))`` argument order:

    - ``body`` -- the symbolic expression in the dummy variable
    - ``dummy`` -- the symbolic variable representing a generic root
    - ``poly`` -- the polynomial whose roots index the sum (expressed
      in the dummy variable)

    EXAMPLES::

        sage: from sage.symbolic.rootsum import root_sum
        sage: var('x a r')
        (x, a, r)
        sage: root_sum(log(x - r) / (a + 3*r^2), r, r^3 + a*r + 1)
        root_sum(log(x - r)/(3*r^2 + a), r, r^3 + a*r + 1)

    LaTeX renders the canonical math notation::

        sage: latex(root_sum(log(x - r) / (3*r^2 + a), r, r^3 + a*r + 1))
        \sum_{r \,:\, r^{3} + a r + 1 = 0} \frac{\log\left(x - r\right)}{3 \, r^{2} + a}

    The expression is intentionally inert.  No automatic simplification or
    evaluation is performed -- even when the polynomial has explicit
    rational roots, the sum stays symbolic::

        sage: root_sum(r^2, r, r^2 - 1)
        root_sum(r^2, r, r^2 - 1)

    Sage receives ``root_sum`` expressions from Maxima when
    ``integrate_use_rootsof:true`` is active; see :issue:`40356`.

    .. NOTE::

        Bridges for SymPy's ``RootSum`` and FriCAS's ``RootSum``, as well
        as ``_derivative_`` and ``_evalf_``, are intentionally not
        implemented in this initial version and are tracked as follow-up
        work.
    """
    def __init__(self):
        r"""
        EXAMPLES::

            sage: from sage.symbolic.rootsum import root_sum
            sage: loads(dumps(root_sum))
            root_sum
        """
        BuiltinFunction.__init__(self, "root_sum", nargs=3)

    def _eval_(self, body, dummy, poly):
        r"""
        Symbolic evaluation. Returns ``None`` to keep the expression
        unevaluated.

        TESTS::

            sage: from sage.symbolic.rootsum import root_sum
            sage: var('r')
            r
            sage: rs = root_sum(r^2, r, r^2 - 1)
            sage: rs.operator() is root_sum
            True
        """
        return None

    def _print_latex_(self, body, dummy, poly):
        r"""
        Render as ``\sum_{r \,:\, p(r) = 0} body``.

        EXAMPLES::

            sage: from sage.symbolic.rootsum import root_sum
            sage: var('a r')
            (a, r)
            sage: latex(root_sum(1/r, r, r^3 - a))
            \sum_{r \,:\, r^{3} - a = 0} \frac{1}{r}
        """
        from sage.misc.latex import latex
        return r"\sum_{%s \,:\, %s = 0} %s" % (latex(dummy),
                                               latex(poly),
                                               latex(body))


root_sum = Function_root_sum()
