# sage_setup: distribution = sagemath-schemes
r"""
Eta-products on modular curves `X_0(N)`

This package provides a class for representing eta-products, which
are meromorphic functions on modular curves of the form

.. MATH::

    \prod_{d | N} \eta(q^d)^{r_d}

where `\eta(q)` is Dirichlet's eta function

.. MATH::

    q^{1/24} \prod_{n = 1}^\infty(1-q^n) .

These are useful for obtaining explicit models of modular curves.

See :issue:`3934` for background.

AUTHOR:

- David Loeffler (2008-08-22): initial version
"""

# ***************************************************************************
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#                     2008 David Loeffler <d.loeffler.01@cantab.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ***************************************************************************
from __future__ import annotations

from sage.arith.misc import divisors, prime_divisors, euler_phi, is_square, gcd
from sage.categories.groups import Groups
from sage.matrix.constructor import matrix
from sage.modules.free_module import FreeModule
from sage.rings.finite_rings.integer_mod import Mod
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.rational_field import QQ
from sage.structure.element import Element
from sage.structure.formal_sum import FormalSum
from sage.structure.parent import Parent
from sage.structure.richcmp import richcmp, richcmp_method, op_EQ, op_NE
from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation

import weakref

_cache = {}


def EtaGroup(level):
    r"""
    Create the group of eta products of the given level.

    EXAMPLES::

        sage: EtaGroup(12)
        Group of eta products on X_0(12)
        sage: EtaGroup(1/2)
        Traceback (most recent call last):
        ...
        TypeError: Level (=1/2) must be a positive integer
        sage: EtaGroup(0)
        Traceback (most recent call last):
        ...
        ValueError: Level (=0) must be a positive integer
    """
    if level in _cache:
        G = _cache[level]()
        if G is not None:
            return G
    G = EtaGroup_class(level)
    _cache[level] = weakref.ref(G)
    return G


class EtaGroupElement(Element):

    def __init__(self, parent, rdict):
        r"""
        Create an eta product object. Usually called implicitly via
        EtaGroup_class.__call__ or the EtaProduct factory function.

        EXAMPLES::

            sage: EtaProduct(8, {1:24, 2:-24})
            Eta product of level 8 : (eta_1)^24 (eta_2)^-24
            sage: g = _; g == loads(dumps(g))
            True
            sage: TestSuite(g).run()
        """
        self._N = parent._N
        N = self._N

        if isinstance(rdict, EtaGroupElement):
            rdict = rdict._rdict
            # Note: This is needed because the "x in G" test tries to call G(x)
            # and see if it returns an error. So sometimes this will be getting
            # called with rdict being an eta product, not a dictionary.

        if rdict == 1:
            rdict = {}
        # Check Ligozat criteria
        sumR = sumDR = sumNoverDr = 0
        prod = 1

        for d in list(rdict):
            if N % d:
                raise ValueError("%s does not divide %s" % (d, N))

            if rdict[d] == 0:
                del rdict[d]
                continue
            sumR += rdict[d]
            sumDR += rdict[d] * d
            sumNoverDr += rdict[d] * N / d
            prod *= (N // d)**rdict[d]

        if sumR != 0:
            raise ValueError("sum r_d (=%s) is not 0" % sumR)
        if sumDR % 24:
            raise ValueError("sum d r_d (=%s) is not 0 mod 24" % sumDR)
        if sumNoverDr % 24:
            raise ValueError("sum (N/d) r_d (=%s) is not 0 mod 24" % sumNoverDr)
        if not is_square(prod):
            raise ValueError("product (N/d)^(r_d) (=%s) is not a square" % prod)

        self._sumDR = ZZ(sumDR)  # this is useful to have around
        self._rdict = rdict

        Element.__init__(self, parent)

    def _mul_(self, other):
        r"""
        Return the product of ``self`` and ``other``.

        EXAMPLES::

            sage: eta1, eta2 = EtaGroup(4).basis() # indirect doctest
            sage: eta1 * eta2
            Eta product of level 4 : (eta_1)^24 (eta_2)^-48 (eta_4)^24
        """
        newdict = {d: self._rdict.get(d, 0) + other._rdict.get(d, 0)
                   for d in set(self._rdict).union(other._rdict)}
        P = self.parent()
        return P.element_class(P, newdict)

    def _div_(self, other):
        r"""
        Return `self * other^{-1}`.

        EXAMPLES::

            sage: eta1, eta2 = EtaGroup(4).basis()
            sage: eta1 / eta2 # indirect doctest
            Eta product of level 4 : (eta_1)^-8 (eta_4)^8
            sage: (eta1 / eta2) * eta2 == eta1
            True
        """
        newdict = {d: self._rdict.get(d, 0) - other._rdict.get(d, 0)
                   for d in set(self._rdict).union(other._rdict)}
        P = self.parent()
        return P.element_class(P, newdict)

    def __invert__(self):
        r"""
        Return the inverse of ``self``.

        EXAMPLES::

            sage: eta1, eta2 = EtaGroup(4).basis()
            sage: ~eta2  # indirect doctest
            Eta product of level 4 : (eta_1)^-16 (eta_2)^24 (eta_4)^-8
        """
        newdict = {d: -self._rdict[d] for d in self._rdict}
        P = self.parent()
        return P.element_class(P, newdict)

    def is_one(self) -> bool:
        r"""
        Return whether ``self`` is the one of the monoid.

        EXAMPLES::

            sage: e = EtaProduct(3, {3:12, 1:-12})
            sage: e.is_one()
            False
            sage: e.parent().one().is_one()
            True
            sage: ep = EtaProduct(5, {})
            sage: ep.is_one()
            True
            sage: ep.parent().one() == ep
            True
        """
        return not self._rdict

    def _richcmp_(self, other, op):
        r"""
        Compare ``self`` to ``other``.

        Eta products are compared according to their rdicts.

        EXAMPLES::

            sage: EtaProduct(2, {2:24,1:-24}) == 1
            False
            sage: EtaProduct(6, {1:-24, 2:24}) == EtaProduct(6, {1:-24, 2:24})
            True
            sage: EtaProduct(6, {1:-24, 2:24}) == EtaProduct(6, {1:24, 2:-24})
            False
            sage: EtaProduct(6, {1:-24, 2:24}) < EtaProduct(6, {1:-24, 2:24, 3:24, 6:-24})
            True
            sage: EtaProduct(6, {1:-24, 2:24, 3:24, 6:-24}) < EtaProduct(6, {1:-24, 2:24})
            False
        """
        if op in [op_EQ, op_NE]:
            test = (self._N == other._N and
                    self._rdict == other._rdict)
            return test == (op == op_EQ)
        return richcmp((self._N, sorted(self._rdict.items())),
                       (other._N, sorted(other._rdict.items())), op)

    def _short_repr(self) -> str:
        r"""
        A short string representation of ``self``, which does not specify the
        level.

        EXAMPLES::

            sage: EtaProduct(3, {3:12, 1:-12})._short_repr()
            '(eta_1)^-12 (eta_3)^12'
        """
        if self.degree() == 0:
            return "1"
        return " ".join("(eta_%s)^%s" % (d, exp)
                        for d, exp in sorted(self._rdict.items()))

    def _repr_(self) -> str:
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: EtaProduct(3, {3:12, 1:-12})._repr_()
            'Eta product of level 3 : (eta_1)^-12 (eta_3)^12'
        """
        return "Eta product of level %s : " % self.level() + self._short_repr()

    def level(self) -> Integer:
        r"""
        Return the level of this eta product.

        EXAMPLES::

            sage: e = EtaProduct(3, {3:12, 1:-12})
            sage: e.level()
            3
            sage: EtaProduct(12, {6:6, 2:-6}).level() # not the lcm of the d's
            12
            sage: EtaProduct(36, {6:6, 2:-6}).level() # not minimal
            36
        """
        return self._N

    def q_expansion(self, n):
        r"""
        Return the `q`-expansion of ``self`` at the cusp at infinity.

        INPUT:

        - ``n`` -- integer; number of terms to calculate

        OUTPUT:

        A power series over `\ZZ` in the variable `q`, with a *relative*
        precision of `1 + O(q^n)`.

        ALGORITHM: Calculates eta to (n/m) terms, where m is the smallest
        integer dividing self.level() such that self.r(m) != 0. Then
        multiplies.

        EXAMPLES::

            sage: EtaProduct(36, {6:6, 2:-6}).q_expansion(10)
            q + 6*q^3 + 27*q^5 + 92*q^7 + 279*q^9 + O(q^11)
            sage: R.<q> = ZZ[[]]
            sage: EtaProduct(2,{2:24,1:-24}).q_expansion(100) == delta_qexp(101)(q^2)/delta_qexp(101)(q)
            True
        """
        R, q = PowerSeriesRing(ZZ, 'q').objgen()
        pr = R.one().O(n)
        if not self._rdict:  # if self.is_one():
            return pr
        # self.r(d) should always be nonzero since we filtered out the 0s
        eta_n = max(n // d for d in self._rdict)  # if self.r(d)
        eta = qexp_eta(R, eta_n)
        for d in self._rdict:
            rd = self._rdict[d]
            if rd:
                pr *= eta(q ** d) ** ZZ(rd)
        return pr * q**(self._sumDR // 24)

    def qexp(self, n):
        """
        Alias for ``self.q_expansion()``.

        EXAMPLES::

            sage: e = EtaProduct(36, {6:8, 3:-8})
            sage: e.qexp(10)
            q + 8*q^4 + 36*q^7 + O(q^10)
            sage: e.qexp(30) == e.q_expansion(30)
            True
        """
        return self.q_expansion(n)

    def order_at_cusp(self, cusp: CuspFamily) -> Integer:
        r"""
        Return the order of vanishing of ``self`` at the given cusp.

        INPUT:

        - ``cusp`` -- a :class:`CuspFamily` object

        OUTPUT: integer

        EXAMPLES::

            sage: e = EtaProduct(2, {2:24, 1:-24})
            sage: e.order_at_cusp(CuspFamily(2, 1)) # cusp at infinity
            1
            sage: e.order_at_cusp(CuspFamily(2, 2)) # cusp 0
            -1
        """
        if not isinstance(cusp, CuspFamily):
            raise TypeError("argument (=%s) should be a CuspFamily" % cusp)
        if cusp.level() != self._N:
            raise ValueError("cusp not on right curve")
        sigma = sum(ell * self._rdict[ell] / cusp.width() *
                    (gcd(cusp.width(), self._N // ell))**2
                    for ell in self._rdict)
        return sigma / ZZ(24) / gcd(cusp.width(), self._N // cusp.width())

    def divisor(self):
        r"""
        Return the divisor of ``self``, as a formal sum of CuspFamily objects.

        EXAMPLES::

            sage: e = EtaProduct(12, {1:-336, 2:576, 3:696, 4:-216, 6:-576, 12:-144})
            sage: e.divisor() # FormalSum seems to print things in a random order?
            -131*(Inf) - 50*(c_{2}) + 11*(0) + 50*(c_{6}) + 169*(c_{4}) - 49*(c_{3})
            sage: e = EtaProduct(2^8, {8:1,32:-1})
            sage: e.divisor() # random
            -(c_{2}) - (Inf) - (c_{8,2}) - (c_{8,3}) - (c_{8,4}) - (c_{4,2})
             - (c_{8,1}) - (c_{4,1}) + (c_{32,4}) + (c_{32,3}) + (c_{64,1})
             + (0) + (c_{32,2}) + (c_{64,2}) + (c_{128}) + (c_{32,1})
        """
        return FormalSum([(self.order_at_cusp(c), c)
                          for c in AllCusps(self.level())])

    def degree(self) -> Integer:
        r"""
        Return the degree of ``self`` as a map `X_0(N) \to \mathbb{P}^1`.

        This is the sum of all the positive coefficients in the divisor
        of ``self``.

        EXAMPLES::

            sage: e = EtaProduct(12, {1:-336, 2:576, 3:696, 4:-216, 6:-576, 12:-144})
            sage: e.degree()
            230
        """
        return sum(self.order_at_cusp(c)
                   for c in AllCusps(self._N)
                   if self.order_at_cusp(c) > 0)

    def r(self, d) -> Integer:
        r"""
        Return the exponent `r_d` of `\eta(q^d)` in ``self``.

        EXAMPLES::

            sage: e = EtaProduct(12, {2:24, 3:-24})
            sage: e.r(3)
            -24
            sage: e.r(4)
            0
        """
        return self._rdict.get(d, 0)


class EtaGroup_class(UniqueRepresentation, Parent):
    r"""
    The group of eta products of a given level under multiplication.

    TESTS::

        sage: TestSuite(EtaGroup(12)).run()

        sage: EtaGroup(12) == EtaGroup(12)
        True
        sage: EtaGroup(12) == EtaGroup(13)
        False

        sage: EtaGroup(12) != EtaGroup(12)
        False
        sage: EtaGroup(12) != EtaGroup(13)
        True
    """

    def __init__(self, level):
        r"""
        Create the group of eta products of a given level, which must be a
        positive integer.

        EXAMPLES::

            sage: G = EtaGroup(12); G # indirect doctest
            Group of eta products on X_0(12)
            sage: TestSuite(G).run()
        """
        try:
            level = ZZ(level)
        except TypeError:
            raise TypeError("Level (=%s) must be a positive integer" % level)
        if level < 1:
            raise ValueError("Level (=%s) must be a positive integer" % level)
        self._N = level
        Parent.__init__(self, category=Groups().Commutative())

    def _repr_(self) -> str:
        r"""
        String representation of ``self``.

        EXAMPLES::

            sage: EtaGroup(12)._repr_()
            'Group of eta products on X_0(12)'
        """
        return "Group of eta products on X_0(%s)" % self.level()

    def one(self) -> EtaGroupElement:
        r"""
        Return the identity element of ``self``.

        EXAMPLES::

            sage: EtaGroup(12).one()
            Eta product of level 12 : 1
        """
        return self({})

    def _element_constructor_(self, dic):
        r"""
        Create an element of this group (an eta product object) with
        exponents from the given dictionary.

        INPUT:

        - ``dic`` -- dictionary

        See the docstring of :func:`EtaProduct` for how ``dic`` is used.

        EXAMPLES::

            sage: EtaGroup(2).__call__({1:24, 2:-24})
            Eta product of level 2 : (eta_1)^24 (eta_2)^-24
        """
        return self.element_class(self, dic)

    def level(self) -> Integer:
        r"""
        Return the level of ``self``.

        EXAMPLES::

            sage: EtaGroup(10).level()
            10
        """
        return self._N

    def basis(self, reduce=True) -> list:
        r"""
        Produce a basis for the free abelian group of eta-products of level
        N (under multiplication), attempting to find basis vectors of the
        smallest possible degree.

        INPUT:

        - ``reduce`` -- boolean (default: ``True``); whether or not to apply
          LLL-reduction to the calculated basis

        EXAMPLES::

            sage: EtaGroup(5).basis()
            [Eta product of level 5 : (eta_1)^6 (eta_5)^-6]
            sage: EtaGroup(12).basis()
            [Eta product of level 12 : (eta_1)^-3 (eta_2)^2 (eta_3)^1 (eta_4)^-1 (eta_6)^-2 (eta_12)^3,
             Eta product of level 12 : (eta_1)^-4 (eta_2)^2 (eta_3)^4 (eta_6)^-2,
             Eta product of level 12 : (eta_1)^6 (eta_2)^-9 (eta_3)^-2 (eta_4)^3 (eta_6)^3 (eta_12)^-1,
             Eta product of level 12 : (eta_1)^-1 (eta_2)^3 (eta_3)^3 (eta_4)^-2 (eta_6)^-9 (eta_12)^6,
             Eta product of level 12 : (eta_1)^3 (eta_3)^-1 (eta_4)^-3 (eta_12)^1]
            sage: EtaGroup(12).basis(reduce=False) # much bigger coefficients
            [Eta product of level 12 : (eta_1)^384 (eta_2)^-576 (eta_3)^-696 (eta_4)^216 (eta_6)^576 (eta_12)^96,
             Eta product of level 12 : (eta_2)^24 (eta_12)^-24,
             Eta product of level 12 : (eta_1)^-40 (eta_2)^116 (eta_3)^96 (eta_4)^-30 (eta_6)^-80 (eta_12)^-62,
             Eta product of level 12 : (eta_1)^-4 (eta_2)^-33 (eta_3)^-4 (eta_4)^1 (eta_6)^3 (eta_12)^37,
             Eta product of level 12 : (eta_1)^15 (eta_2)^-24 (eta_3)^-29 (eta_4)^9 (eta_6)^24 (eta_12)^5]

        ALGORITHM: An eta product of level `N` is uniquely
        determined by the integers `r_d` for `d | N` with
        `d < N`, since `\sum_{d | N} r_d = 0`. The valid
        `r_d` are those that satisfy two congruences modulo 24,
        and one congruence modulo 2 for every prime divisor of N. We beef
        up the congruences modulo 2 to congruences modulo 24 by multiplying
        by 12. To calculate the kernel of the ensuing map
        `\ZZ^m \to (\ZZ/24\ZZ)^n`
        we lift it arbitrarily to an integer matrix and calculate its Smith
        normal form. This gives a basis for the lattice.

        This lattice typically contains "large" elements, so by default we
        pass it to the reduce_basis() function which performs
        LLL-reduction to give a more manageable basis.
        """
        N = self.level()
        divs = divisors(N)[:-1]
        s = len(divs)
        primedivs = prime_divisors(N)

        rows = []
        for di in divs:
            # generate a row of relation matrix
            row = [Mod(di, 24) - Mod(N, 24), Mod(N // di, 24) - Mod(1, 24)]
            for p in primedivs:
                row.append(Mod(12 * (N // di).valuation(p), 24))
            rows.append(row)

        M = matrix(IntegerModRing(24), rows)
        Mlift = M.change_ring(ZZ)
        # now we compute elementary factors of Mlift
        S, U, V = Mlift.smith_form()
        good_vects = []
        for i in range(U.nrows()):
            vect = U.row(i)
            nf = (i < S.ncols() and S[i, i]) or 0  # ?
            good_vects.append((vect * 24 / gcd(nf, 24)).list())
        for v in good_vects:
            v.append(-sum(list(v)))
        dicts = []
        for v in good_vects:
            dicts.append({})
            for i in range(s):
                dicts[-1][divs[i]] = v[i]
            dicts[-1][N] = v[-1]
        if reduce:
            return self.reduce_basis([self(d) for d in dicts])
        else:
            return [self(d) for d in dicts]

    def reduce_basis(self, long_etas) -> list:
        r"""
        Produce a more manageable basis via LLL-reduction.

        INPUT:

        - ``long_etas`` -- a list of EtaGroupElement objects (which
          should all be of the same level)

        OUTPUT:

        - a new list of EtaGroupElement objects having
          hopefully smaller norm

        ALGORITHM: We define the norm of an eta-product to be the
        `L^2` norm of its divisor (as an element of the free
        `\ZZ`-module with the cusps as basis and the
        standard inner product). Applying LLL-reduction to this gives a
        basis of hopefully more tractable elements. Of course we'd like to
        use the `L^1` norm as this is just twice the degree, which
        is a much more natural invariant, but `L^2` norm is easier
        to work with!

        EXAMPLES::

            sage: EtaGroup(4).reduce_basis([ EtaProduct(4, {1:8,2:24,4:-32}), EtaProduct(4, {1:8, 4:-8})])
            [Eta product of level 4 : (eta_1)^8 (eta_4)^-8,
             Eta product of level 4 : (eta_1)^-8 (eta_2)^24 (eta_4)^-16]
        """
        N = self.level()
        cusps = AllCusps(N)
        r = matrix(ZZ, [[et.order_at_cusp(c) for c in cusps] for et in long_etas])
        V = FreeModule(ZZ, r.ncols())
        A = V.submodule_with_basis([V(rw) for rw in r.rows()])
        rred = r.LLL()
        short_etas = []
        for shortvect in rred.rows():
            bv = A.coordinates(shortvect)
            dic = {d: sum(bv[i] * long_etas[i].r(d)
                          for i in range(r.nrows()))
                   for d in divisors(N)}
            short_etas.append(self(dic))
        return short_etas

    Element = EtaGroupElement


def EtaProduct(level, dic) -> EtaGroupElement:
    r"""
    Create an :class:`EtaGroupElement` object representing the function
    `\prod_{d | N} \eta(q^d)^{r_d}`.

    This checks the criteria of Ligozat to ensure that this product
    really is the `q`-expansion of a meromorphic function on `X_0(N)`.

    INPUT:

    - ``level`` -- integer; the N such that this eta
      product is a function on `X_0(N)`

    - ``dic`` -- a dictionary indexed by divisors of N such that the
      coefficient of `\eta(q^d)` is r[d]. Only nonzero coefficients need be
      specified. If Ligozat's criteria are not satisfied, a :exc:`ValueError`
      will be raised.

    OUTPUT:

    An EtaGroupElement object, whose parent is the EtaGroup of level N and
    whose coefficients are the given dictionary.

    .. NOTE::

        The dictionary ``dic`` does not uniquely specify N. It is
        possible for two EtaGroupElements with different `N`'s to
        be created with the same dictionary, and these represent different
        objects (although they will have the same `q`-expansion at
        the cusp `\infty`).

    EXAMPLES::

        sage: EtaProduct(3, {3:12, 1:-12})
        Eta product of level 3 : (eta_1)^-12 (eta_3)^12
        sage: EtaProduct(3, {3:6, 1:-6})
        Traceback (most recent call last):
        ...
        ValueError: sum d r_d (=12) is not 0 mod 24
        sage: EtaProduct(3, {4:6, 1:-6})
        Traceback (most recent call last):
        ...
        ValueError: 4 does not divide 3
    """
    return EtaGroup(level)(dic)


def num_cusps_of_width(N, d) -> Integer:
    r"""
    Return the number of cusps on `X_0(N)` of width ``d``.

    INPUT:

    - ``N`` -- integer; the level

    - ``d`` -- integer; an integer dividing `N`, the cusp width

    EXAMPLES::

        sage: from sage.modular.etaproducts import num_cusps_of_width
        sage: [num_cusps_of_width(18,d) for d in divisors(18)]
        [1, 1, 2, 2, 1, 1]
        sage: num_cusps_of_width(4,8)
        Traceback (most recent call last):
        ...
        ValueError: N and d must be positive integers with d|N
    """
    N = ZZ(N)
    d = ZZ(d)
    if N <= 0 or d <= 0 or N % d:
        raise ValueError("N and d must be positive integers with d|N")

    return euler_phi(d.gcd(N // d))


def AllCusps(N) -> list:
    r"""
    Return a list of CuspFamily objects corresponding to the cusps of
    `X_0(N)`.

    INPUT:

    - ``N`` -- integer; the level

    EXAMPLES::

        sage: AllCusps(18)
        [(Inf), (c_{2}), (c_{3,1}), (c_{3,2}), (c_{6,1}), (c_{6,2}), (c_{9}), (0)]
        sage: AllCusps(0)
        Traceback (most recent call last):
        ...
        ValueError: N must be positive
    """
    N = ZZ(N)
    if N <= 0:
        raise ValueError("N must be positive")

    c = []
    for d in divisors(N):
        n = num_cusps_of_width(N, d)
        if n == 1:
            c.append(CuspFamily(N, d))
        elif n > 1:
            for i in range(n):
                c.append(CuspFamily(N, d, label=str(i + 1)))
    return c


@richcmp_method
class CuspFamily(SageObject):
    r"""
    A family of elliptic curves parametrising a region of `X_0(N)`.
    """
    def __init__(self, N, width, label=None):
        r"""
        Create the cusp of width d on X_0(N) corresponding to the family
        of Tate curves `(\CC_p/q^d, \langle \zeta q\rangle)`.

        Here `\zeta` is a primitive root of unity of order `r` with
        `\mathrm{lcm}(r,d) = N`. The cusp does not store zeta, so we
        store an arbitrary label instead.

        EXAMPLES::

            sage: CuspFamily(8, 4)
            (c_{4})
            sage: CuspFamily(16, 4, '1')
            (c_{4,1})
        """
        N = ZZ(N)
        if N <= 0:
            raise ValueError("N must be positive")
        self._N = N
        self._width = width
        if N % width:
            raise ValueError("bad width")
        if num_cusps_of_width(N, width) > 1 and label is None:
            raise ValueError("there are %s > 1 cusps of width %s on X_0(%s): specify a label" % (num_cusps_of_width(N, width), width, N))
        if num_cusps_of_width(N, width) == 1 and label is not None:
            raise ValueError("there is only one cusp of width %s on X_0(%s): no need to specify a label" % (width, N))
        self.label = label

    @property
    def __tuple(self):
        """
        The defining data of this ``CuspFamily`` as tuple, used for
        comparisons.
        """
        return (self._N, self._width, self.label)

    def __richcmp__(self, other, op) -> bool:
        """
        EXAMPLES::

            sage: a = CuspFamily(16, 4, "1"); a
            (c_{4,1})
            sage: b = CuspFamily(16, 4, "2"); b
            (c_{4,2})
            sage: c = CuspFamily(8, 8); c
            (0)
            sage: a == a
            True
            sage: a == b
            False
            sage: a != b
            True
            sage: a == c
            False
            sage: a < c
            False
            sage: a > c
            True
            sage: a != "foo"
            True
        """
        if not isinstance(other, CuspFamily):
            return NotImplemented
        return richcmp(self.__tuple, other.__tuple, op)

    def __hash__(self):
        """
        EXAMPLES::

            sage: hash(CuspFamily(10, 1))  # random
            -4769758480201659164
        """
        return hash(self.__tuple)

    def width(self) -> Integer:
        r"""
        Return the width of this cusp.

        EXAMPLES::

            sage: e = CuspFamily(10, 1)
            sage: e.width()
            1
        """
        return self._width

    def level(self) -> Integer:
        r"""
        Return the level of this cusp.

        EXAMPLES::

            sage: e = CuspFamily(10, 1)
            sage: e.level()
            10
        """
        return self._N

    def sage_cusp(self):
        r"""
        Return the corresponding element of `\mathbb{P}^1(\QQ)`.

        EXAMPLES::

            sage: CuspFamily(10, 1).sage_cusp() # not implemented
            Infinity
        """
        raise NotImplementedError

    def _repr_(self) -> str:
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: CuspFamily(16, 4, "1")._repr_()
            '(c_{4,1})'
        """
        if self.width() == 1:
            return "(Inf)"
        elif self.width() == self.level():
            return "(0)"
        return "(c_{%s%s})" % (self.width(), ((self.label and ("," + self.label)) or ""))


def qexp_eta(ps_ring, prec):
    r"""
    Return the `q`-expansion of `\eta(q) / q^{1/24}`.

    Here `\eta(q)` is Dedekind's function

    .. MATH::

        \eta(q) = q^{1/24}\prod_{n=1}^\infty (1-q^n).

    The result is an element of ``ps_ring``, with precision ``prec``.

    INPUT:

    - ``ps_ring`` -- PowerSeriesRing; a power series ring

    - ``prec`` -- integer; the number of terms to compute

    OUTPUT: an element of ``ps_ring`` which is the `q`-expansion of
    `\eta(q)/q^{1/24}` truncated to prec terms.

    ALGORITHM: We use the Euler identity

    .. MATH::

         \eta(q) = q^{1/24}( 1 + \sum_{n \ge 1} (-1)^n (q^{n(3n+1)/2} + q^{n(3n-1)/2})

    to compute the expansion.

    EXAMPLES::

        sage: from sage.modular.etaproducts import qexp_eta
        sage: qexp_eta(ZZ[['q']], 100)
        1 - q - q^2 + q^5 + q^7 - q^12 - q^15 + q^22 + q^26 - q^35 - q^40 + q^51 + q^57 - q^70 - q^77 + q^92 + O(q^100)
    """
    prec = Integer(prec)
    if not prec > 0:
        raise ValueError("prec must be a positive integer")
    v = [Integer(0)] * prec
    pm = Integer(1)
    v[0] = pm
    try:
        n = 1
        while True:
            pm = -pm
            v[n * (3 * n - 1) // 2] = pm
            v[n * (3 * n + 1) // 2] = pm
            n += 1
    except IndexError:
        pass
    return ps_ring(v, prec=prec)


def eta_poly_relations(eta_elements, degree, labels=['x1', 'x2'],
                       verbose=False):
    r"""
    Find polynomial relations between eta products.

    INPUT:

    - ``eta_elements`` -- list; a list of EtaGroupElement objects.
      Not implemented unless this list has precisely two elements. degree

    - ``degree`` -- integer; the maximal degree of polynomial to look for

    - ``labels`` -- list of strings; labels to use for the polynomial returned

    - ``verbose`` -- boolean (default: ``False``); if ``True``, prints information as
      it goes

    OUTPUT: list of polynomials which is a Groebner basis for the
    part of the ideal of relations between eta_elements which is
    generated by elements up to the given degree; or None, if no
    relations were found.

    ALGORITHM: An expression of the form
    `\sum_{0 \le i,j \le d} a_{ij} x^i y^j` is zero if and
    only if it vanishes at the cusp infinity to degree at least
    `v = d(deg(x) + deg(y))`. For all terms up to `q^v`
    in the `q`-expansion of this expression to be zero is a
    system of `v + k` linear equations in `d^2`
    coefficients, where `k` is the number of nonzero negative
    coefficients that can appear.

    Solving these equations and calculating a basis for the solution
    space gives us a set of polynomial relations, but this is generally
    far from a minimal generating set for the ideal, so we calculate a
    Groebner basis.

    As a test, we calculate five extra terms of `q`-expansion
    and check that this doesn't change the answer.

    EXAMPLES::

        sage: from sage.modular.etaproducts import eta_poly_relations
        sage: t = EtaProduct(26, {2:2,13:2,26:-2,1:-2})
        sage: u = EtaProduct(26, {2:4,13:2,26:-4,1:-2})
        sage: eta_poly_relations([t, u], 3)
        sage: eta_poly_relations([t, u], 4)
        [x1^3*x2 - 13*x1^3 - 4*x1^2*x2 - 4*x1*x2 - x2^2 + x2]

    Use ``verbose=True`` to see the details of the computation::

        sage: eta_poly_relations([t, u], 3, verbose=True)
        Trying to find a relation of degree 3
        Lowest order of a term at infinity = -12
        Highest possible degree of a term = 15
        Trying all coefficients from q^-12 to q^15 inclusive
        No polynomial relation of order 3 valid for 28 terms
        Check:
        Trying all coefficients from q^-12 to q^20 inclusive
        No polynomial relation of order 3 valid for 33 terms

    ::

        sage: eta_poly_relations([t, u], 4, verbose=True)
        Trying to find a relation of degree 4
        Lowest order of a term at infinity = -16
        Highest possible degree of a term = 20
        Trying all coefficients from q^-16 to q^20 inclusive
        Check:
        Trying all coefficients from q^-16 to q^25 inclusive
        [x1^3*x2 - 13*x1^3 - 4*x1^2*x2 - 4*x1*x2 - x2^2 + x2]
    """
    if len(eta_elements) > 2:
        raise NotImplementedError("do not know how to find relations between more than two elements")

    eta1, eta2 = eta_elements

    if verbose:
        print("Trying to find a relation of degree %s" % degree)
    inf = CuspFamily(eta1.level(), 1)
    loterm = -(min(0, eta1.order_at_cusp(inf)) + min(0, eta2.order_at_cusp(inf))) * degree
    if verbose:
        print("Lowest order of a term at infinity = %s" % -loterm)

    maxdeg = sum([eta1.degree(), eta2.degree()]) * degree
    if verbose:
        print("Highest possible degree of a term = %s" % maxdeg)
    m = loterm + maxdeg + 1
    oldgrob = _eta_relations_helper(eta1, eta2, degree, m, labels, verbose)
    if verbose:
        print("Check:")
    newgrob = _eta_relations_helper(eta1, eta2, degree, m + 5, labels, verbose)
    if oldgrob != newgrob:
        raise ArithmeticError("Check: answers different!")
    return newgrob


def _eta_relations_helper(eta1, eta2, degree, qexp_terms, labels, verbose):
    r"""
    Helper function used by eta_poly_relations. Finds a basis for the
    space of linear relations between the first qexp_terms of the
    `q`-expansions of the monomials
    `\eta_1^i * \eta_2^j` for `0 \le i,j < degree`,
    and calculates a Groebner basis for the ideal generated by these
    relations.

    Liable to return meaningless results if qexp_terms isn't at least
    `1 + d*(m_1,m_2)` where

    .. MATH::

       m_i = min(0, {\text degree of the pole of $\eta_i$ at $\infty$})

    as then 1 will be in the ideal.

    EXAMPLES::

        sage: from sage.modular.etaproducts import _eta_relations_helper
        sage: r,s = EtaGroup(4).basis()
        sage: _eta_relations_helper(r,s,4,100,['a','b'],False)
        [a + 1/16*b - 1/16]
        sage: _eta_relations_helper(EtaProduct(26, {2:2,13:2,26:-2,1:-2}),EtaProduct(26, {2:4,13:2,26:-4,1:-2}),3,12,['a','b'],False) # not enough terms, will return rubbish
        [1]
    """
    indices = [(i, j) for j in range(degree) for i in range(degree)]
    inf = CuspFamily(eta1.level(), 1)

    pole_at_infinity = -(min(0, eta1.order_at_cusp(inf)) + min(0, eta2.order_at_cusp(inf))) * degree
    if verbose:
        print("Trying all coefficients from q^%s to q^%s inclusive" % (-pole_at_infinity, -pole_at_infinity + qexp_terms - 1))

    rows = []
    for j in range(qexp_terms):
        rows.append([])
    for i in indices:
        func = (eta1**i[0] * eta2**i[1]).qexp(qexp_terms)
        for j in range(qexp_terms):
            rows[j].append(func[j - pole_at_infinity])
    M = matrix(rows)
    V = M.right_kernel()
    if V.dimension() == 0:
        if verbose:
            print("No polynomial relation of order %s valid for %s terms" % (degree, qexp_terms))
        return None
    if V.dimension() >= 1:
        R = PolynomialRing(QQ, 2, labels)
        x, y = R.gens()
        relations = [sum([c[v] * x**indices[v][0] * y**indices[v][1]
                          for v in range(len(indices))])
                     for c in V.basis()]
        id = R.ideal(relations)
        return id.groebner_basis()
