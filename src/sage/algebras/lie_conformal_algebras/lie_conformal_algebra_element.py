# sage.doctest: needs sage.combinat sage.modules
"""
Lie Conformal Algebra Element

AUTHORS:

- Reimundo Heluani (2019-08-09): Initial implementation.
"""
# *****************************************************************************
#       Copyright (C) 2019 Reimundo Heluani <heluani@potuz.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.arith.misc import factorial
from sage.misc.misc_c import prod
from sage.misc.repr import repr_lincomb
from sage.misc.latex import latex
from sage.modules.with_basis.indexed_element import IndexedFreeModuleElement


class LCAWithGeneratorsElement(IndexedFreeModuleElement):
    """
    The element class of a Lie conformal algebra with a
    preferred set of generators.
    """
    def T(self, n=1):
        r"""
        The `n`-th derivative of this element.

        INPUT:

        - ``n`` -- nonnegative integer (default: `1`); how many
          times to apply `T` to this element

        We use the *divided powers* notation
        `T^{(j)} = \frac{T^j}{j!}`.

        EXAMPLES::

            sage: Vir = lie_conformal_algebras.Virasoro(QQ)
            sage: Vir.inject_variables()
            Defining L, C
            sage: L.T()
            TL
            sage: L.T(3)
            6*T^(3)L
            sage: C.T()
            0

            sage: R = lie_conformal_algebras.NeveuSchwarz(QQbar); R.inject_variables()
            Defining L, G, C
            sage: (L + 2*G.T() + 4*C).T(2)
            2*T^(2)L + 12*T^(3)G
        """
        from sage.rings.integer_ring import ZZ
        if n not in ZZ or n < 0:
            raise ValueError("n must be a nonnegative Integer")
        if n == 0 or self.is_zero():
            return self
        # it's faster to sum than to use recursion
        if self.is_monomial():
            p = self.parent()
            a, m = self.index()
            coef = self._monomial_coefficients[(a, m)]
            if (a, m + n) in p._indices:
                return coef * prod(range(m + 1, m + n + 1))\
                    * p.monomial((a, m + n))
            else:
                return p.zero()
        return sum(mon.T(n) for mon in self.terms())

    def is_monomial(self):
        """
        Whether this element is a monomial.

        EXAMPLES::

            sage: Vir = lie_conformal_algebras.Virasoro(QQ); L = Vir.0
            sage: (L + L.T()).is_monomial()
            False
            sage: L.T().is_monomial()
            True
        """
        return len(self._monomial_coefficients) == 1 or self.is_zero()


class LCAStructureCoefficientsElement(LCAWithGeneratorsElement):
    """
    An element of a Lie conformal algebra given by structure
    coefficients.
    """
    def _bracket_(self, right):
        """
        The lambda bracket of these two elements.

        The result is a dictionary with nonnegative integer keys.
        The value corresponding to the entry `j` is ``self_{(j)}right``.

        EXAMPLES::

            sage: Vir = lie_conformal_algebras.Virasoro(QQ); L = Vir.0
            sage: L.bracket(L)
            {0: TL, 1: 2*L, 3: 1/2*C}
            sage: L.T().bracket(L)
            {1: -TL, 2: -4*L, 4: -2*C}

            sage: R = lie_conformal_algebras.Affine(QQbar, 'A1', names=('e','h','f')); R
            The affine Lie conformal algebra of type ['A', 1] over Algebraic Field
            sage: R.inject_variables()
            Defining e, h, f, K
            sage: e.bracket(f)
            {0: h, 1: K}
            sage: h.bracket(h.T())
            {2: 4*K}
        """
        p = self.parent()
        if self.is_monomial() and right.is_monomial():
            if self.is_zero() or right.is_zero():
                return {}
            s_coeff = p._s_coeff
            a, k = self.index()
            coefa = self.monomial_coefficients()[(a, k)]
            b, m = right.index()
            coefb = right.monomial_coefficients()[(b, m)]
            try:
                mbr = dict(s_coeff[(a, b)])
            except KeyError:
                return {}
            pole = max(mbr.keys())
            ret = {l: coefa * coefb * (-1)**k / factorial(k) *
                   sum(factorial(l) / factorial(m + k + j - l)
                       / factorial(l - k - j) / factorial(j)
                       * mbr[j].T(m + k + j - l)
                       for j in mbr if l - m - k <= j <= l - k)
                   for l in range(m + k + pole + 1)}
            return {k: v for k, v in ret.items() if v}

        diclist = [i._bracket_(j) for i in self.terms()
                   for j in right.terms()]
        ret = {}
        pz = p.zero()
        for d in diclist:
            for k in d:
                ret[k] = ret.get(k, pz) + d[k]
        return {k: v for k, v in ret.items() if v}

    def _repr_(self):
        r"""
        A visual representation of this element.

        For a free generator `L`, the element `\frac{T^{j}}{j!}L` is
        denoted by ``T^(j)L``.

        EXAMPLES::

            sage: V = lie_conformal_algebras.Virasoro(QQ); V.inject_variables()
            Defining L, C
            sage: v = L.T(5).nproduct(L,6); v
            -1440*L
            sage: L.T(2) + L + C
            2*T^(2)L + L + C
            sage: L.T(4)
            24*T^(4)L

            sage: R = lie_conformal_algebras.Affine(QQ, 'B3')
            sage: R.2.T()+3*R.3
            TB[alpha[1]] + 3*B[alpha[2] + alpha[3]]  # 64-bit
            ...  # 32-bit
        """
        if self.is_zero():
            return "0"
        p = self.parent()
        if p._names:
            terms = [("T^({}){}".format(k1, p._names[p._index_to_pos[k0]]), v) if k1 > 1
                     else ("T{}".format(p._names[p._index_to_pos[k0]]), v) if k1 == 1
                     else ("{}".format(p._names[p._index_to_pos[k0]]), v)
                     for (k0, k1), v in self.monomial_coefficients().items()]
        else:
            terms = [("T^({}){}".format(k1, p._repr_generator(k0)), v) if k1 > 1
                     else ("T{}".format(p._repr_generator(k0)), v) if k1 == 1
                     else ("{}".format(p._repr_generator(k0)), v)
                     for (k0, k1), v in self.monomial_coefficients().items()]

        return repr_lincomb(terms, strip_one=True)

    def _latex_(self):
        r"""
        A visual representation of this element.

        For a free generator `L`, the element `\frac{T^{j}}{j!}L` is
        denoted by ``T^(j)L``.

        EXAMPLES::

            sage: V = lie_conformal_algebras.Virasoro(QQ); V.inject_variables()
            Defining L, C
            sage: latex(L.T(2))
            2 T^{(2)}L

            sage: R = lie_conformal_algebras.Affine(QQbar, 'A1', names=('e','h','f')); R.inject_variables()
            Defining e, h, f, K
            sage: latex(e.bracket(f))
            \left\{0 : h, 1 : K\right\}
            sage: latex(e.T(3))
            6 T^{(3)}e

            sage: R = lie_conformal_algebras.Affine(QQbar, 'A1')
            sage: latex(R.0.bracket(R.2))
            \left\{0 : \alpha^\vee_{1}, 1 : \text{\texttt{K}}\right\}

            sage: R = lie_conformal_algebras.Affine(QQ, 'A1'); latex(R.0.T(3))
            6 T^{(3)}\alpha_{1}
        """
        if self.is_zero():
            return "0"
        p = self.parent()
        try:
            names = p.latex_variable_names()
        except ValueError:
            names = None
        if names:
            terms = [("T^{{({})}}{}".format(k1, names[p._index_to_pos[k0]]), v) if k1 > 1
                     else ("T{}".format(names[p._index_to_pos[k0]]), v) if k1 == 1
                     else ("{}".format(names[p._index_to_pos[k0]]), v)
                     for (k0, k1), v in self.monomial_coefficients().items()]
        else:
            terms = [("T^{{({})}}{}".format(k1, latex(k0)), v) if k1 > 1
                     else ("T{}".format(latex(k0)), v) if k1 == 1
                     else ("{}".format(latex(k0)), v)
                     for (k0, k1), v in self.monomial_coefficients().items()]

        return repr_lincomb(terms, is_latex=True, strip_one=True)
