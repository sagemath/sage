# -*- coding: utf-8 -*-
r"""
Element class for Pollack-Stevens' modular symbols

This is the class of elements in the spaces of Pollack-Steven's modular symbols as described in [PS2011]_.

EXAMPLES::

    sage: E = EllipticCurve('11a')
    sage: phi = E.pollack_stevens_modular_symbol(); phi
    Modular symbol of level 11 with values in Sym^0 Q^2
    sage: phi.weight() # Note that weight k=2 of a modular form corresponds here to weight 0
    0
    sage: phi.values()
    [-1/5, 1, 0]
    sage: phi.is_ordinary(11)
    True
    sage: phi_lift = phi.lift(11, 5, eigensymbol = True) # long time
    sage: phi_lift.padic_lseries().series(5) # long time
    O(11^5) + (10 + 3*11 + 6*11^2 + 9*11^3 + O(11^4))*T + (6 + 3*11 + 2*11^2 + O(11^3))*T^2 + (2 + 2*11 + O(11^2))*T^3 + (5 + O(11))*T^4 + O(T^5)

::

    sage: A = ModularSymbols(Gamma1(8),4).decomposition()[0].plus_submodule().new_subspace()
    sage: from sage.modular.pollack_stevens.space import ps_modsym_from_simple_modsym_space
    sage: phi = ps_modsym_from_simple_modsym_space(A)
    sage: phi.values()
    [(-1, 0, 0), (1, 0, 0), (-9, -6, -4)]

"""
# ****************************************************************************
#        Copyright (C) 2012 Robert Pollack <rpollack@math.bu.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

import operator

from sage.arith.misc import next_prime, gcd, kronecker
from sage.categories.action import Action
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import
from sage.misc.verbose import verbose
from sage.rings.integer_ring import ZZ
from sage.rings.padics.precision_error import PrecisionError
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.structure.element import ModuleElement
from sage.structure.richcmp import op_EQ, op_NE

lazy_import('sage.rings.padics.factory', 'Qp')
lazy_import('sage.rings.padics.padic_generic', 'pAdicGeneric')

from .manin_map import ManinMap
from .sigma0 import Sigma0
from .fund_domain import M2Z


minusproj = [1, 0, 0, -1]


def _iterate_Up(Phi, p, M, ap, q, aq, check):
    r"""
    Return an overconvergent Hecke-eigensymbol lifting self -- self must be a
    `p`-ordinary eigensymbol

    INPUT:

    - ``p`` -- prime

    - ``M`` -- integer equal to the number of moments

    - ``ap`` -- Hecke eigenvalue at `p`

    - ``q`` -- prime

    - ``aq`` -- Hecke eigenvalue at `q`

    OUTPUT:

    - Hecke-eigenvalue overconvergent modular symbol lifting self.

    EXAMPLES::

        sage: E = EllipticCurve('57a')
        sage: p = 5
        sage: prec = 4
        sage: phi = E.pollack_stevens_modular_symbol()
        sage: phi_stabilized = phi.p_stabilize(p,M = prec)
        sage: Phi = phi_stabilized.lift(p,prec) # indirect doctest
    """
    if ap.valuation(p) > 0:
        raise ValueError("Lifting non-ordinary eigensymbols not implemented (issue #20)")

    # Act by Hecke to ensure values are in D and not D^dag after solving difference equation
    verbose("Applying Hecke", level=2)

    apinv = ~ap
    Phi = apinv * Phi.hecke(p)

    # Killing eisenstein part
    verbose("Killing eisenstein part with q = %s" % q, level=2)
    k = Phi.parent().weight()
    Phi = ((q ** (k + 1) + 1) * Phi - Phi.hecke(q))

    # Iterating U_p
    verbose("Iterating U_p", level=2)
    Psi = apinv * Phi.hecke(p)

    for attempts in range(M-1):
        verbose("%s attempt (val = %s/%s)" % (attempts + 1,(Phi-Psi).valuation(),M), level=2)
        Phi = Psi
        Psi = apinv * Phi.hecke(p)
        Psi._normalize()
    Phi = ~(q ** (k + 1) + 1 - aq) * Phi
    return Phi

class PSModSymAction(Action):
    def __init__(self, actor, MSspace):
        r"""
        Create the action

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: phi = E.pollack_stevens_modular_symbol()
            sage: g = phi._map._codomain._act._Sigma0(matrix(ZZ,2,2,[1,2,3,4]))
            sage: phi * g # indirect doctest
            Modular symbol of level 11 with values in Sym^0 Q^2
        """

        Action.__init__(self, actor, MSspace, False, operator.mul)

    def _act_(self, g, sym):
        r"""
        Return the result of sym * g

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: phi = E.pollack_stevens_modular_symbol()
            sage: g = phi._map._codomain._act._Sigma0(matrix(ZZ,2,2,[2,1,5,-1]))
            sage: phi * g # indirect doctest
            Modular symbol of level 11 with values in Sym^0 Q^2
        """

        return sym.__class__(sym._map * g, sym.parent(), construct=True)


class PSModularSymbolElement(ModuleElement):
    def __init__(self, map_data, parent, construct=False):
        r"""
        Initialize a modular symbol

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: phi = E.pollack_stevens_modular_symbol()
        """
        ModuleElement.__init__(self, parent)
        if construct:
            self._map = map_data
        else:
            self._map = ManinMap(parent._coefficients, parent._source, map_data)

    def _repr_(self):
        r"""
        Return the string representation of the symbol.

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: phi = E.pollack_stevens_modular_symbol()
            sage: phi._repr_()
            'Modular symbol of level 11 with values in Sym^0 Q^2'
        """
        return "Modular symbol of level %s with values in %s" % (self.parent().level(), self.parent().coefficient_module())

    def dict(self):
        r"""
        Return dictionary on the modular symbol self, where keys are generators and values are the corresponding values of self on generators

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: phi = E.pollack_stevens_modular_symbol()
            sage: Set([x.moment(0) for x in phi.dict().values()]) == Set([-1/5, 1, 0])
            True
        """
        D = {}
        for g in self.parent().source().gens():
            D[g] = self._map[g]
        return D

    def weight(self):
        r"""
        Return the weight of this Pollack-Stevens modular symbol.

        This is `k-2`, where `k` is the usual notion of weight for modular
        forms!

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: phi = E.pollack_stevens_modular_symbol()
            sage: phi.weight()
            0
        """
        return self.parent().weight()

    def values(self):
        r"""
        Return the values of the symbol ``self`` on our chosen generators.

        The generators are listed in ``self.dict()``.

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: phi = E.pollack_stevens_modular_symbol()
            sage: phi.values()
            [-1/5, 1, 0]
            sage: sorted(phi.dict())
            [
            [-1 -1]  [ 0 -1]  [1 0]
            [ 3  2], [ 1  3], [0 1]
            ]
            sage: sorted(phi.values()) == sorted(phi.dict().values())
            True
        """
        return [self._map[g] for g in self.parent().source().gens()]

    def _normalize(self, **kwds):
        """
        Normalize all of the values of the symbol ``self``.

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: phi = E.pollack_stevens_modular_symbol()
            sage: phi._normalize()
            Modular symbol of level 11 with values in Sym^0 Q^2
            sage: phi._normalize().values()
            [-1/5, 1, 0]
        """
        for val in self._map:
            val.normalize(**kwds)
        return self

    def _richcmp_(self, other, op):
        """
        Check if self == other.

        Here self and other have the same parent.

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: phi = E.pollack_stevens_modular_symbol()
            sage: phi == phi
            True
            sage: phi == 2*phi
            False
            sage: psi = EllipticCurve('37a').pollack_stevens_modular_symbol()
            sage: psi == phi
            False
        """
        if op not in [op_EQ, op_NE]:
            return NotImplemented

        b = all(self._map[g] == other._map[g]
                for g in self.parent().source().gens())

        return b == (op == op_EQ)

    def _add_(self, right):
        """
        Return self + right

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: phi = E.pollack_stevens_modular_symbol()
            sage: phi.values()
            [-1/5, 1, 0]
            sage: phi + phi
            Modular symbol of level 11 with values in Sym^0 Q^2
            sage: (phi + phi).values()
            [-2/5, 2, 0]
        """
        return self.__class__(self._map + right._map, self.parent(), construct=True)

    def _lmul_(self, right):
        """
        Return self * right

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: phi = E.pollack_stevens_modular_symbol()
            sage: phi.values()
            [-1/5, 1, 0]
            sage: 2*phi
            Modular symbol of level 11 with values in Sym^0 Q^2
            sage: (2*phi).values()
            [-2/5, 2, 0]
        """
        return self.__class__(self._map * right, self.parent(), construct=True)

    def _rmul_(self, right):
        """
        Return self * right

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: phi = E.pollack_stevens_modular_symbol()
            sage: phi.values()
            [-1/5, 1, 0]
            sage: phi*2
            Modular symbol of level 11 with values in Sym^0 Q^2
            sage: (phi*2).values()
            [-2/5, 2, 0]
        """
        return self.__class__(self._map * right, self.parent(), construct=True)

    def _sub_(self, right):
        """
        Return self - right

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: phi = E.pollack_stevens_modular_symbol()
            sage: phi.values()
            [-1/5, 1, 0]
            sage: phi - phi
            Modular symbol of level 11 with values in Sym^0 Q^2
            sage: (phi - phi).values()
            [0, 0, 0]
        """
        return self.__class__(self._map - right._map, self.parent(), construct=True)

    def _get_prime(self, p=None, alpha=None, allow_none=False):
        """
        Combine a prime specified by the user with the prime from the parent.

        INPUT:

        - ``p`` -- an integer or None (default None); if specified
          needs to match the prime of the parent.

        - ``alpha`` -- an element or None (default None); if p-adic
          can contribute a prime.

        - ``allow_none`` -- boolean (default False); whether to allow
          no prime to be specified.

        OUTPUT:

        - a prime or None.  If ``allow_none`` is False then a
          ``ValueError`` will be raised rather than returning None if no
          prime can be determined.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Symk
            sage: D = OverconvergentDistributions(0, 5, 10)
            sage: M = PollackStevensModularSymbols(Gamma0(5), coefficients=D)
            sage: f = M(1); f._get_prime()
            5
            sage: f._get_prime(5)
            5
            sage: f._get_prime(7)
            Traceback (most recent call last):
            ...
            ValueError: inconsistent prime
            sage: f._get_prime(alpha=Qp(5)(1))
            5
            sage: D = Symk(0)
            sage: M = PollackStevensModularSymbols(Gamma0(2), coefficients=D)
            sage: f = M(1); f._get_prime(allow_none=True) is None
            True
            sage: f._get_prime(alpha=Qp(7)(1))
            7
            sage: f._get_prime(7,alpha=Qp(7)(1))
            7
            sage: f._get_prime()
            Traceback (most recent call last):
            ...
            ValueError: you must specify a prime
        """
        pp = self.parent().prime()
        ppp = ((alpha is not None) and hasattr(alpha.parent(), 'prime')
               and alpha.parent().prime()) or None
        p = ZZ(p) or pp or ppp
        if not p:
            if not allow_none:
                raise ValueError("you must specify a prime")
        elif (pp and p != pp) or (ppp and p != ppp):
            raise ValueError("inconsistent prime")
        return p

    def plus_part(self):
        r"""
        Return the plus part of self -- i.e. ``self + self | [1,0,0,-1]``.

        Note that we haven't divided by 2.  Is this a problem?

        OUTPUT:

        - self + self | [1,0,0,-1]

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: phi = E.pollack_stevens_modular_symbol()
            sage: phi.values()
            [-1/5, 1, 0]
            sage: (phi.plus_part()+phi.minus_part()) == 2 * phi
            True
        """
        S0N = Sigma0(self.parent().level())
        return self + self * S0N(minusproj)

    def minus_part(self):
        r"""
        Return the minus part of self -- i.e. self - self | [1,0,0,-1]

        Note that we haven't divided by 2.  Is this a problem?

        OUTPUT:

        - self -- self | [1,0,0,-1]

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: phi = E.pollack_stevens_modular_symbol()
            sage: phi.values()
            [-1/5, 1, 0]
            sage: (phi.plus_part()+phi.minus_part()) == phi * 2
            True
        """
        S0N = Sigma0(self.parent().level())
        return self - self * S0N(minusproj)

    def hecke(self, ell, algorithm="prep"):
        r"""
        Return self | `T_{\ell}` by making use of the precomputations in
        self.prep_hecke()

        INPUT:

        - ``ell`` -- a prime

        - ``algorithm`` -- a string, either 'prep' (default) or
          'naive'

        OUTPUT:

        - The image of this element under the Hecke operator
          `T_{\ell}`

        ALGORITHMS:

        - If ``algorithm == 'prep'``, precomputes a list of matrices
          that only depend on the level, then uses them to speed up
          the action.

        - If ``algorithm == 'naive'``, just acts by the matrices
          defining the Hecke operator.  That is, it computes
          sum_a self | [1,a,0,ell] + self | [ell,0,0,1],
          the last term occurring only if the level is prime to ell.

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: phi = E.pollack_stevens_modular_symbol()
            sage: phi.values()
            [-1/5, 1, 0]
            sage: phi.hecke(2) == phi * E.ap(2)
            True
            sage: phi.hecke(3) == phi * E.ap(3)
            True
            sage: phi.hecke(5) == phi * E.ap(5)
            True
            sage: phi.hecke(101) == phi * E.ap(101)
            True

            sage: all(phi.hecke(p, algorithm='naive') == phi * E.ap(p) for p in [2,3,5,101]) # long time
            True
        """
        return self.__class__(self._map.hecke(ell, algorithm),
                                  self.parent(), construct=True)

    def valuation(self, p=None):
        r"""
        Return the valuation of ``self`` at `p`.

        Here the valuation is the minimum of the valuations of the
        values of ``self``.

        INPUT:

        - ``p`` - prime

        OUTPUT:

        - The valuation of ``self`` at `p`

        EXAMPLES::

           sage: E = EllipticCurve('11a')
           sage: phi = E.pollack_stevens_modular_symbol()
           sage: phi.values()
           [-1/5, 1, 0]
           sage: phi.valuation(2)
           0
           sage: phi.valuation(3)
           0
           sage: phi.valuation(5)
           -1
           sage: phi.valuation(7)
           0
           sage: phi.valuation()
           Traceback (most recent call last):
           ...
           ValueError: you must specify a prime

           sage: phi2 = phi.lift(11, M=2)
           sage: phi2.valuation()
           0
           sage: phi2.valuation(3)
           Traceback (most recent call last):
           ...
           ValueError: inconsistent prime
           sage: phi2.valuation(11)
           0
        """
        q = self._get_prime(p)
        return min([val.valuation(q) for val in self._map])

    def diagonal_valuation(self, p):
        """
        Return the minimum of the diagonal valuation on the values of self

        INPUT:

        - ``p`` -- a positive integral prime

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: phi = E.pollack_stevens_modular_symbol()
            sage: phi.values()
            [-1/5, 1, 0]
            sage: phi.diagonal_valuation(2)
            0
            sage: phi.diagonal_valuation(3)
            0
            sage: phi.diagonal_valuation(5)
            -1
            sage: phi.diagonal_valuation(7)
            0
        """
        return min([val.diagonal_valuation(p) for val in self._map])

    @cached_method
    def is_Tq_eigensymbol(self, q, p=None, M=None):
        r"""
        Determine if self is an eigenvector for `T_q` modulo `p^M`

        INPUT:

        - ``q`` -- prime of the Hecke operator

        - ``p`` -- prime we are working modulo

        - ``M`` -- degree of accuracy of approximation

        OUTPUT:

        - True/False

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: phi = E.pollack_stevens_modular_symbol()
            sage: phi.values()
            [-1/5, 1, 0]
            sage: phi_ord = phi.p_stabilize(p = 3, ap = E.ap(3), M = 10, ordinary = True)
            sage: phi_ord.is_Tq_eigensymbol(2,3,10)
            True
            sage: phi_ord.is_Tq_eigensymbol(2,3,100)
            False
            sage: phi_ord.is_Tq_eigensymbol(2,3,1000)
            False
            sage: phi_ord.is_Tq_eigensymbol(3,3,10)
            True
            sage: phi_ord.is_Tq_eigensymbol(3,3,100)
            False
        """
        try:
            self.Tq_eigenvalue(q, p, M)
            return True
        except ValueError:
            return False

    # what happens if a cached method raises an error?  Is it
    # recomputed each time?
    @cached_method
    def Tq_eigenvalue(self, q, p=None, M=None, check=True):
        r"""
        Eigenvalue of `T_q` modulo `p^M`

        INPUT:

        - ``q`` -- prime of the Hecke operator

        - ``p`` -- prime we are working modulo (default: None)

        - ``M`` -- degree of accuracy of approximation (default: None)

        - ``check`` -- check that ``self`` is an eigensymbol

        OUTPUT:

        - Constant `c` such that `self|T_q - c * self` has valuation greater than
          or equal to `M` (if it exists), otherwise raises ValueError

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: phi = E.pollack_stevens_modular_symbol()
            sage: phi.values()
            [-1/5, 1, 0]
            sage: phi_ord = phi.p_stabilize(p = 3, ap = E.ap(3), M = 10, ordinary = True)
            sage: phi_ord.Tq_eigenvalue(2,3,10) + 2
            O(3^10)

            sage: phi_ord.Tq_eigenvalue(3,3,10)
            2 + 3^2 + 2*3^3 + 2*3^4 + 2*3^6 + 3^8 + 2*3^9 + O(3^10)
            sage: phi_ord.Tq_eigenvalue(3,3,100)
            Traceback (most recent call last):
            ...
            ValueError: result not determined to high enough precision
        """
        qhecke = self.hecke(q)
        gens = self.parent().source().gens()
        if p is None:
            p = self.parent().prime()
        i = 0

        g = gens[i]
        verbose("Computing eigenvalue", level=2)
        while self._map[g].moment(0).is_zero():
            if not qhecke._map[g].moment(0).is_zero():
                raise ValueError("not a scalar multiple")
            i += 1
            try:
                g = gens[i]
            except IndexError:
                raise ValueError("self is zero")
        aq = self.parent().base_ring()(self._map[g].find_scalar_from_zeroth_moment(qhecke._map[g], p, M, check))

        verbose("Found eigenvalues of %s" % aq, level=2)
        if check:
            verbose("Checking that this is actually an eigensymbol", level=2)
            if p is None or M is None or not ZZ(p).is_prime():
                for g in gens[1:]:
                    try:
                        if not (qhecke._map[g] - aq * self._map[g]).is_zero():
                            # using != did not work
                            raise ValueError("not a scalar multiple")
                    except PrecisionError:
                        if qhecke._map[g] != aq * self._map[g]:
                            raise ValueError("not a scalar multiple")
            else:
                verbose('p = %s, M = %s' % (p, M), level=2)
                if qhecke != aq * self:
                    raise ValueError("not a scalar multiple")
        # if not aq.parent().is_exact() and M is not None:
        #     aq.add_bigoh(M)
        return aq

    def is_ordinary(self, p=None, P=None):
        r"""
        Return true if the `p`-th eigenvalue is a `p`-adic unit.

        INPUT:

        - ``p`` - a positive integral prime, or None (default None)
        - ``P`` - a prime of the base ring above `p`, or None. This is ignored
          unless the base ring is a number field.

        OUTPUT:

        - True/False

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: phi = E.pollack_stevens_modular_symbol()
            sage: phi.is_ordinary(2)
            False
            sage: E.ap(2)
            -2
            sage: phi.is_ordinary(3)
            True
            sage: E.ap(3)
            -1
            sage: phip = phi.p_stabilize(3,20)
            sage: phip.is_ordinary()
            True

        A number field example. Here there are multiple primes above `p`, and
        `\phi` is ordinary at one but not the other.::

            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_simple_modsym_space
            sage: f = Newforms(32, 8, names='a')[1]
            sage: phi = ps_modsym_from_simple_modsym_space(f.modular_symbols(1))
            sage: (p1, _), (p2, _) = phi.base_ring().ideal(3).factor()
            sage: phi.is_ordinary(p1) != phi.is_ordinary(p2)
            True
            sage: phi.is_ordinary(3)
            Traceback (most recent call last):
            ...
            TypeError: P must be an ideal
        """
        # q is the prime below p, if base is a number field; q = p otherwise
        if p is None:
            if self.parent().prime() == 0:
                raise ValueError("need to specify a prime")
            q = p = self.parent().prime()
        elif p in ZZ:
            q = p
        else:
            q = p.smallest_integer()
        if not q.is_prime():
            raise ValueError("p is not prime")
        if (self.parent().prime() != q) and (self.parent().prime() != 0):
            raise ValueError("prime does not match coefficient module's prime")
        aq = self.Tq_eigenvalue(q)
        return aq.valuation(p) == 0

    def evaluate_twisted(self, a, chi):
        r"""
        Return `\Phi_{\chi}(\{a/p\}-\{\infty\})` where `\Phi` is ``self`` and
        `\chi` is a quadratic character

        INPUT:

        - ``a`` -- integer in the range range(p)
        - ``chi`` -- the modulus of a quadratic character.

        OUTPUT:

        The distribution `\Phi_{\chi}(\{a/p\}-\{\infty\})`.

        EXAMPLES::

            sage: E = EllipticCurve('17a1')
            sage: L = E.padic_lseries(5, implementation="pollackstevens", precision=4) #long time
            sage: D = L.quadratic_twist()          # long time
            sage: L.symbol().evaluate_twisted(1,D) # long time
            (1 + 5 + 3*5^2 + 5^3 + O(5^4), 5^2 + O(5^3), 1 + O(5^2), 2 + O(5))

            sage: E = EllipticCurve('40a4')
            sage: L = E.padic_lseries(7, implementation="pollackstevens", precision=4) #long time
            sage: D = L.quadratic_twist()          # long time
            sage: L.symbol().evaluate_twisted(1,D) # long time
            (4 + 6*7 + 3*7^2 + O(7^4), 6*7 + 6*7^2 + O(7^3), 6 + O(7^2), 1 + O(7))

        TESTS:

        Check for :trac:`32878`::

            sage: E = EllipticCurve('11a1')
            sage: L = E.padic_lseries(3, implementation="pollackstevens", precision=4)
            sage: D = 5
            sage: L.symbol().evaluate_twisted(1, D)
            (2 + 3 + 2*3^2 + O(3^4), 2 + 3 + O(3^3), 2 + 3 + O(3^2), 2 + O(3))
        """
        p = self.parent().prime()
        S0p = Sigma0(p)
        Dists = self.parent().coefficient_module()
        M = Dists.precision_cap()
        p = Dists.prime()
        twisted_dist = Dists.zero()
        m_map = self._map
        for b in range(1, abs(chi) + 1):
            if gcd(b, chi) == 1:
                M1 = S0p([1, (b / abs(chi)) % p**M, 0, 1])
                new_dist = m_map(M2Z([a * abs(chi) + p * b,
                                      1, p * abs(chi), 0])) * M1
                new_dist = new_dist.scale(kronecker(chi, b)).normalize()
                twisted_dist += new_dist
        return twisted_dist.normalize()

    def _consistency_check(self):
        """
        Check that the map really does satisfy the Manin relations loop (for debugging).
        The two and three torsion relations are checked and it is checked that the symbol
        adds up correctly around the fundamental domain

        EXAMPLES::

            sage: E = EllipticCurve('37a1')
            sage: phi = E.pollack_stevens_modular_symbol()
            sage: phi._consistency_check()
            This modular symbol satisfies the Manin relations
        """
        f = self._map
        MR = self._map._manin
        # Test two torsion relations
        for g in MR.reps_with_two_torsion():
            gamg = MR.two_torsion_matrix(g)
            if not (f[g] * gamg + f[g]).is_zero():
                raise ValueError("Two torsion relation failed with", g)

        # Test three torsion relations
        for g in MR.reps_with_three_torsion():
            gamg = MR.three_torsion_matrix(g)
            if not (f[g] * (gamg ** 2) + f[g] * gamg + f[g]).is_zero():
                raise ValueError("Three torsion relation failed with", g)

        # Test that the symbol adds to 0 around the boundary of the
        # fundamental domain
        t = self.parent().coefficient_module().zero()
        for g in MR.gens()[1:]:
            if not (g in MR.reps_with_two_torsion()
                    or g in MR.reps_with_three_torsion()):
                t += f[g] * MR.gammas[g] - f[g]
            else:
                if g in MR.reps_with_two_torsion():
                    t -= f[g]
                else:
                    t -= f[g]   # what ?? same thing ??

        id = MR.gens()[0]
        if f[id] * MR.gammas[id] - f[id] != -t:
            print(t)
            print(f[id] * MR.gammas[id] - f[id])
            raise ValueError("does not add up correctly around loop")

        print("This modular symbol satisfies the Manin relations")


class PSModularSymbolElement_symk(PSModularSymbolElement):
    def _find_alpha(self, p, k, M=None, ap=None, new_base_ring=None, ordinary=True, check=True, find_extraprec=True):
        r"""
        Find `\alpha`, a `U_p` eigenvalue, which is found as a root of
        the polynomial `x^2 - a_p * x + p^{k+1} \chi(p)`.

        INPUT:

        - ``p`` -- prime

        - ``k`` -- Pollack-Stevens weight

        - ``M`` -- precision (default: None) of `\QQ_p`

        - ``ap`` -- Hecke eigenvalue at `p` (default: None)

        - ``new_base_ring`` -- field of definition of `\alpha` (default: None)

        - ``ordinary`` -- True if the prime is ordinary (default: True)

        - ``check`` -- check to see if the prime is ordinary (default: True)

        - ``find_extraprec`` -- setting this to True finds extra precision (default: True)

        OUTPUT:

        The output is a tuple (``alpha``, ``new_base_ring``,
        ``newM``, ``eisenloss``,``q``,``aq``), with

        - ``alpha`` --  `U_p` eigenvalue

        - ``new_base_ring`` -- field of definition of `\alpha` with precision at least ``newM``

        - ``newM`` -- new precision

        - ``eisenloss`` -- loss of precision

        - ``q`` -- a prime not equal to `p` which was used to find extra precision

        - ``aq`` -- the Hecke eigenvalue `a_q` corresponding to `q`

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: p = 5
            sage: M = 10
            sage: k = 0
            sage: phi = E.pollack_stevens_modular_symbol()
            sage: phi._find_alpha(p,k,M)
            (1 + 4*5 + 3*5^2 + 2*5^3 + 4*5^4 + 4*5^5 + 4*5^6 + 3*5^7 + 2*5^8 + 3*5^9 + 3*5^10 + 3*5^12 + 2*5^13 + O(5^14), 5-adic Field with capped relative precision 14, 13, 1, 2, -2)
        """
        if ap is None:
            ap = self.Tq_eigenvalue(p, check=check)
        if check and ap.valuation(p) > 0:
            raise ValueError("p is not ordinary")

        chi = self._map._codomain._character
        if chi is not None:
            eps = chi(p)
        else:
            eps = 1
        poly = PolynomialRing(ap.parent(), 'x')([p ** (k + 1) * eps, -ap, 1])
        if new_base_ring is None:
            # These should actually be completions of disc.parent()
            if p == 2:
                # is this the right precision adjustment for p=2?
                new_base_ring = Qp(2, M + 1)
            else:
                new_base_ring = Qp(p, M)
            set_padicbase = True
        else:
            set_padicbase = False
        try:
            verbose("finding alpha: rooting %s in %s" % (poly, new_base_ring), level=2)
            poly = poly.change_ring(new_base_ring)
            (v0, e0), (v1, e1) = poly.roots()
        except (TypeError, ValueError):
            raise ValueError("new base ring must contain a root of x^2 - ap * x + p^(k+1)")
        if v0.valuation(p) > 0:
            v0, v1 = v1, v0
        if ordinary:
            alpha = v0
        else:
            alpha = v1
        if find_extraprec:
            newM, eisenloss, q, aq = self._find_extraprec(p, M, alpha, check)
        else:
            newM, eisenloss, q, aq = M, None, None, None
        if set_padicbase:
            # We want to ensure that the relative precision of alpha
            # and (alpha-1) are both at least *newM*, where newM is
            # obtained from self._find_extraprec
            prec_cap = None
            verbose("testing prec_rel: newM = %s, alpha = %s" % (newM, alpha),
                    level=2)
            if alpha.precision_relative() < newM:
                prec_cap = newM + alpha.valuation(p) + (1 if p == 2 else 0)
            if ordinary:
                a1val = (alpha - 1).valuation(p)
                verbose("a1val = %s" % a1val, level=2)
                if a1val > 0 and ap != 1 + p ** (k + 1):
                    # if ap = 1 + p**(k+1) then alpha=1 and we need to give up.
                    if prec_cap is None:
                        prec_cap = newM + a1val + (1 if p == 2 else 0)
                    else:
                        prec_cap = max(prec_cap, newM + a1val + (1 if p == 2 else 0))
            verbose("prec_cap = %s" % prec_cap, level=2)
            if prec_cap is not None:
                new_base_ring = Qp(p, prec_cap)
                return self._find_alpha(p=p, k=k, M=M, ap=ap, new_base_ring=new_base_ring, ordinary=ordinary, check=False, find_extraprec=find_extraprec)
        return alpha, new_base_ring, newM, eisenloss, q, aq

    def p_stabilize(self, p=None, M=20, alpha=None, ap=None, new_base_ring=None, ordinary=True, check=True):
        r"""
        Return the `p`-stabilization of self to level `N p` on which `U_p` acts by `\alpha`.

        Note that since `\alpha` is `p`-adic, the resulting symbol
        is just an approximation to the true `p`-stabilization
        (depending on how well `\alpha` is approximated).

        INPUT:

        - ``p`` -- prime not dividing the level of self

        - ``M`` -- (default: 20) precision of `\QQ_p`

        - ``alpha`` -- `U_p` eigenvalue

        - ``ap`` -- Hecke eigenvalue

        - ``new_base_ring`` -- change of base ring

        - ``ordinary`` -- (default: True) whether to return the ordinary
                          (at ``p``) eigensymbol.

        - ``check`` -- (default: True) whether to perform extra sanity checks

        OUTPUT:

        A modular symbol with the same Hecke eigenvalues as
        self away from `p` and eigenvalue `\alpha` at `p`.
        The eigenvalue `\alpha` depends on the parameter ``ordinary``.

        If ``ordinary`` == True: the unique modular symbol of level
        `N p` with the same Hecke eigenvalues as self away from
        `p` and unit eigenvalue at `p`; else  the unique modular
        symbol of level `N p` with the same Hecke eigenvalues as
        self away from `p` and non-unit eigenvalue at `p`.

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: p = 5
            sage: prec = 4
            sage: phi = E.pollack_stevens_modular_symbol()
            sage: phis = phi.p_stabilize(p,M = prec)
            sage: phis
            Modular symbol of level 55 with values in Sym^0 Q_5^2
            sage: phis.hecke(7) == phis*E.ap(7)
            True
            sage: phis.hecke(5) == phis*E.ap(5)
            False
            sage: phis.hecke(3) == phis*E.ap(3)
            True
            sage: phis.Tq_eigenvalue(5)
            1 + 4*5 + 3*5^2 + 2*5^3 + O(5^4)
            sage: phis.Tq_eigenvalue(5,M = 3)
            1 + 4*5 + 3*5^2 + O(5^3)

            sage: phis = phi.p_stabilize(p,M = prec,ordinary=False)
            sage: phis.Tq_eigenvalue(5)
            5 + 5^2 + 2*5^3 + O(5^5)

        A complicated example (with nontrivial character)::

            sage: chi = DirichletGroup(24)([-1, -1, -1])
            sage: f = Newforms(chi,names='a')[0]
            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_simple_modsym_space
            sage: phi = ps_modsym_from_simple_modsym_space(f.modular_symbols(1))
            sage: phi11, h11 = phi.completions(11,20)[0]
            sage: phi11s = phi11.p_stabilize()
            sage: phi11s.is_Tq_eigensymbol(11) # long time
            True
        """
        if check:
            p = self._get_prime(p, alpha)
        k = self.parent().weight()
        M = ZZ(M)
        verbose("p stabilizing: M = %s" % M, level=2)
        if alpha is None:
            alpha, new_base_ring, newM, eisenloss, q, aq = self._find_alpha(p, k, M + 1, ap, new_base_ring, ordinary, check, find_extraprec=False)
            new_base_ring = Qp(p, M) if p != 2 else Qp(p, M + 1)
        else:
            if new_base_ring is None:
                new_base_ring = alpha.parent()
            if check:
                if ap is None:
                    ap = self.base_ring()(alpha + p ** (k + 1) / alpha)
                elif alpha ** 2 - ap * alpha + p ** (k + 1) != 0:
                    raise ValueError("alpha must be a root of x^2 - a_p*x + p^(k+1)")
                if self.hecke(p) != ap * self:
                    raise ValueError("alpha must be a root of x^2 - a_p*x + p^(k+1)")
        verbose("found alpha = %s" % alpha, level=2)

        V = self.parent()._p_stabilize_parent_space(p, new_base_ring)
        return self.__class__(self._map.p_stabilize(p, alpha, V), V, construct=True)

    def completions(self, p, M):
        r"""
        If `K` is the base_ring of self, this function takes all maps
        `K\to \QQ_p` and applies them to self return a list of
        (modular symbol,map: `K\to \QQ_p`) as map varies over all such maps.

        .. NOTE::

            This only returns all completions when `p` splits completely in `K`

        INPUT:

        - ``p`` -- prime

        - ``M`` -- precision

        OUTPUT:

        - A list of tuples (modular symbol,map: `K\to \QQ_p`) as map varies over all such maps

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_simple_modsym_space
            sage: D = ModularSymbols(67,2,1).cuspidal_submodule().new_subspace().decomposition()[1]
            sage: f = ps_modsym_from_simple_modsym_space(D)
            sage: S = f.completions(41,10); S
            [(Modular symbol of level 67 with values in Sym^0 Q_41^2, Ring morphism:
              From: Number Field in alpha with defining polynomial x^2 + 3*x + 1
              To:   41-adic Field with capped relative precision 10
              Defn: alpha |--> 5 + 22*41 + 19*41^2 + 10*41^3 + 28*41^4 + 22*41^5 + 9*41^6 + 25*41^7 + 40*41^8 + 8*41^9 + O(41^10)), (Modular symbol of level 67 with values in Sym^0 Q_41^2, Ring morphism:
              From: Number Field in alpha with defining polynomial x^2 + 3*x + 1
              To:   41-adic Field with capped relative precision 10
              Defn: alpha |--> 33 + 18*41 + 21*41^2 + 30*41^3 + 12*41^4 + 18*41^5 + 31*41^6 + 15*41^7 + 32*41^9 + O(41^10))]
            sage: TestSuite(S[0][0]).run(skip=['_test_category'])
        """
        K = self.base_ring()
        R = Qp(p, M + 10)['x']
        x = R.gen()
        if K == QQ:
            f = x - 1
        else:
            f = K.defining_polynomial()
        v = R(f).roots()
        if len(v) == 0:
            L = Qp(p, M).extension(f, names='a')
            # a = L.gen()
            V = self.parent().change_ring(L)
            Dist = V.coefficient_module()
            psi = K.hom([K.gen()], L)
            embedded_sym = self.parent().element_class(self._map.apply(psi, codomain=Dist, to_moments=True), V, construct=True)
            ans = [embedded_sym, psi]
            return ans
        else:
            roots = [r[0] for r in v]
            ans = []
            V = self.parent().change_ring(Qp(p, M))
            Dist = V.coefficient_module()
            for r in roots:
                psi = K.hom([r], Qp(p, M))
                embedded_sym = self.parent().element_class(self._map.apply(psi, codomain=Dist, to_moments=True), V, construct=True)
                ans.append((embedded_sym, psi))
            return ans

    def lift(self, p=None, M=None, alpha=None, new_base_ring=None,
             algorithm=None, eigensymbol=False, check=True):
        r"""
        Return a (`p`-adic) overconvergent modular symbol with
        `M` moments which lifts self up to an Eisenstein error

        Here the Eisenstein error is a symbol whose system of Hecke
        eigenvalues equals `\ell+1` for `T_\ell` when `\ell`
        does not divide `Np` and 1 for `U_q` when `q` divides `Np`.

        INPUT:

        - ``p`` -- prime

        - ``M`` -- integer equal to the number of moments

        - ``alpha`` -- `U_p` eigenvalue

        - ``new_base_ring`` -- change of base ring

        - ``algorithm`` -- 'stevens' or 'greenberg' (default 'stevens')

        - ``eigensymbol`` -- if True, lifts to Hecke eigensymbol (self must
          be a `p`-ordinary eigensymbol)

        (Note: ``eigensymbol = True`` does *not* just indicate to the code that
        self is an eigensymbol; it solves a wholly different problem, lifting
        an eigensymbol to an eigensymbol.)

        OUTPUT:

        An overconvergent modular symbol whose specialization equals self, up
        to some Eisenstein error if ``eigensymbol`` is False. If ``eigensymbol
        = True`` then the output will be an overconvergent Hecke eigensymbol
        (and it will lift the input exactly, the Eisenstein error disappears).

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: f = E.pollack_stevens_modular_symbol()
            sage: g = f.lift(11,4,algorithm='stevens',eigensymbol=True)
            sage: g.is_Tq_eigensymbol(2)
            True
            sage: g.Tq_eigenvalue(3)
            10 + 10*11 + 10*11^2 + 10*11^3 + O(11^4)
            sage: g.Tq_eigenvalue(11)
            1 + O(11^4)

        We check that lifting and then specializing gives back the original symbol::

            sage: g.specialize() == f
            True

        Another example, which showed precision loss in an earlier version of the code::

            sage: E = EllipticCurve('37a')
            sage: p = 5
            sage: prec = 4
            sage: phi = E.pollack_stevens_modular_symbol()
            sage: Phi = phi.p_stabilize_and_lift(p,prec, algorithm='stevens', eigensymbol=True) # long time
            sage: Phi.Tq_eigenvalue(5,M = 4) # long time
            3 + 2*5 + 4*5^2 + 2*5^3 + O(5^4)

        Another example::

            sage: from sage.modular.pollack_stevens.padic_lseries import pAdicLseries
            sage: E = EllipticCurve('37a')
            sage: p = 5
            sage: prec = 6
            sage: phi = E.pollack_stevens_modular_symbol()
            sage: Phi = phi.p_stabilize_and_lift(p=p,M=prec,alpha=None,algorithm='stevens',eigensymbol=True) #long time
            sage: L = pAdicLseries(Phi)          # long time
            sage: L.symbol() is Phi              # long time
            True

        Examples using Greenberg's algorithm::

            sage: E = EllipticCurve('11a')
            sage: phi = E.pollack_stevens_modular_symbol()
            sage: Phi = phi.lift(11,8,algorithm='greenberg',eigensymbol=True)
            sage: Phi2 = phi.lift(11,8,algorithm='stevens',eigensymbol=True)
            sage: Phi == Phi2
            True

        An example in higher weight::

            sage: from sage.modular.pollack_stevens.space import ps_modsym_from_simple_modsym_space
            sage: f = ps_modsym_from_simple_modsym_space(Newforms(7, 4)[0].modular_symbols(1))
            sage: fs = f.p_stabilize(5)
            sage: FsG = fs.lift(M=6, eigensymbol=True,algorithm='greenberg') # long time
            sage: FsG.values()[0]                                            # long time
            5^-1 * (2*5 + 5^2 + 3*5^3 + 4*5^4 + O(5^7), O(5^6), 2*5^2 + 3*5^3 + O(5^5), O(5^4), 5^2 + O(5^3), O(5^2))
            sage: FsS = fs.lift(M=6, eigensymbol=True,algorithm='stevens')   # long time
            sage: FsS == FsG                                                 # long time
            True
        """
        if p is None:
            p = self.parent().prime()
            if p == 0:
                raise ValueError("must specify a prime")
        elif (self.parent().prime() != 0) and p != self.parent().prime():
            raise ValueError("inconsistent prime")
        if M is None:
            M = self.parent().precision_cap() + 1
        elif M <= 1:
            raise ValueError("M must be at least 2")
        else:
            M = ZZ(M)
        if new_base_ring is None:
            if isinstance(self.parent().base_ring(), pAdicGeneric):
                new_base_ring = self.parent().base_ring()
            else:
                # We may need extra precision in solving the difference equation
                if algorithm == 'greenberg':
                    extraprec = 0
                else:
                    extraprec = (M - 1).exact_log(p)  # DEBUG: was M-1
                # should eventually be a completion
                new_base_ring = Qp(p, M + extraprec)
        if algorithm is None:
            # The default algorithm is Greenberg's, if possible.
            algorithm = 'greenberg' if eigensymbol else 'stevens'
        elif algorithm == 'greenberg':
            if not eigensymbol:
                raise ValueError("Greenberg's algorithm only works"
                                " for eigensymbols. Try 'stevens'")
        elif algorithm != 'stevens':
            raise ValueError("algorithm %s not recognized" % algorithm)
        if eigensymbol:
            # We need some extra precision due to the fact that solving
            # the difference equation can give denominators.
            if alpha is None:
                verbose('Finding alpha with M = %s' % M, level=2)
                alpha = self.Tq_eigenvalue(p, M=M + 1, check=check)
            newM, eisenloss, q, aq = self._find_extraprec(p, M + 1, alpha, check)
            Phi = self._lift_to_OMS(p, newM, new_base_ring, algorithm)
            Phi = _iterate_Up(Phi, p, newM, alpha, q, aq, check)
            Phi = Phi.reduce_precision(M)
            return Phi._normalize(include_zeroth_moment=True)
        else:
            return self._lift_to_OMS(p, M, new_base_ring, algorithm)

    def _lift_to_OMS(self, p, M, new_base_ring, algorithm='greenberg'):
        r"""
        Return a (`p`-adic) overconvergent modular symbol with
        `M` moments which lifts self up to an Eisenstein error

        Here the Eisenstein error is a symbol whose system of Hecke
        eigenvalues equals `\ell+1` for `T_\ell` when `\ell`
        does not divide `Np` and 1 for `U_q` when `q` divides `Np`.

        INPUT:

        - ``p`` -- prime

        - ``M`` -- integer equal to the number of moments

        - ``new_base_ring`` -- new base ring

        - ``algorithm`` -- (default: 'greenberg') a string, either 'greenberg'
          or 'stevens', specifying whether to use
          the lifting algorithm of M.Greenberg or that of Pollack--Stevens.
          The latter one solves the difference equation, which is not needed. The
          option to use Pollack--Stevens' algorithm here is just for historical reasons.

        OUTPUT:

        - An overconvergent modular symbol whose specialization
          equals self up to some Eisenstein error.

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: f = E.pollack_stevens_modular_symbol()
            sage: f._lift_to_OMS(11,4,Qp(11,4))
            Modular symbol of level 11 with values in Space of 11-adic distributions with k=0 action and precision cap 4
        """
        D = {}
        manin = self.parent().source()
        MSS = self.parent()._lift_parent_space(p, M, new_base_ring)
        if algorithm == 'greenberg':
            for g in manin.gens():
                D[g] = self._map[g].lift(p, M, new_base_ring)
        elif algorithm == 'stevens':
            half = ZZ(1) / ZZ(2)
            for g in manin.gens()[1:]:
                twotor = g in manin.reps_with_two_torsion()
                threetor = g in manin.reps_with_three_torsion()
                if twotor:
                    # See [PS2011] section 4.1
                    gam = manin.two_torsion_matrix(g)
                    mu = self._map[g].lift(p, M, new_base_ring)
                    D[g] = (mu - mu * gam) * half
                elif threetor:
                    # See [PS2011] section 4.1
                    gam = manin.three_torsion_matrix(g)
                    mu = self._map[g].lift(p, M, new_base_ring)
                    D[g] = (2 * mu - mu * gam - mu * (gam ** 2)) * half
                else:
                    # no two or three torsion
                    D[g] = self._map[g].lift(p, M, new_base_ring)

            t = self.parent().coefficient_module().lift(p, M, new_base_ring).zero()
            # This loops adds up around the boundary of fundamental
            # domain except the two vertical lines
            for g in manin.gens()[1:]:
                twotor = g in manin.reps_with_two_torsion()
                threetor = g in manin.reps_with_three_torsion()
                if twotor or threetor:
                    t = t - D[g]
                else:
                    t += D[g] * manin.gammas[g] - D[g]
            # t now should be sum Phi(D_i) | (gamma_i - 1) - sum
            # Phi(D'_i) - sum Phi(D''_i)

            # (Here I'm using the opposite sign convention of [PS2011]
            # regarding D'_i and D''_i)

            D[manin.gen(0)] = -t.solve_difference_equation()  # Check this!
        else:
            raise NotImplementedError

        return MSS(D)

    def _find_aq(self, p, M, check):
        r"""
        Helper function for finding Hecke eigenvalue `aq` for a prime `q`
        not equal to `p`. This is called in the case when `alpha = 1 (mod p^M)`
        (with `alpha` a `U_p`-eigenvalue), which creates the need to use
        other Hecke eigenvalues (and `alpha`s), because of division by `(alpha - 1)`.

        INPUT:

        - ``p`` -- working prime

        - ``M`` -- precision

        - ``check`` -- checks that ``self`` is a `T_q` eigensymbol

        OUTPUT:

        Tuple ``(q, aq, eisenloss)``, with

        - ``q`` -- a prime not equal to `p`

        - ``aq`` -- Hecke eigenvalue at `q`

        - ``eisenloss`` -- the `p`-adic valuation of `a_q - q^{k+1} - 1`

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: f = E.pollack_stevens_modular_symbol()
            sage: f._find_aq(5,10,True)
            (2, -2, 1)
        """
        N = self.parent().level()
        q = ZZ(2)
        k = self.parent().weight()
        aq = self.Tq_eigenvalue(q, check=check)
        eisenloss = (aq - q ** (k + 1) - 1).valuation(p)
        while ((q == p) or (N % q == 0) or (eisenloss >= M)) and (q < 50):
            q = next_prime(q)
            aq = self.Tq_eigenvalue(q, check=check)
            if q != p:
                eisenloss = (aq - q ** (k + 1) - 1).valuation(p)
            else:
                eisenloss = (aq - 1).valuation(p)
        if q >= 50:
            raise ValueError("The symbol appears to be eisenstein -- "
                             "not implemented yet")
        return q, aq, eisenloss

    def _find_extraprec(self, p, M, alpha, check):
        r"""
        Find the extra precision needed to account for:

        1) The denominators in the Hecke eigenvalue
        2) the denominators appearing when solving the difference equation,
        3) those denominators who might be also present in ``self``.

        INPUT:

        - ``p`` -- working prime
        - ``M`` -- precision
        - ``alpha`` -- the Up-eigenvalue
        - ``check`` -- whether to check that ``self`` is a `T_q` eigensymbol

        OUTPUT:

        A tuple (newM, eisenloss, q, aq), where ``newM`` is the new precision, `q` is
        a prime different from `p`, and ``aq`` is the eigenvalue of `T_q` of the eigensymbol.
        The value ``eisenloss`` is the loss of precision accounted for in the denominators of the Hecke
        eigenvalue.

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: p = 5
            sage: M = 10
            sage: k = 0
            sage: phi = E.pollack_stevens_modular_symbol()
            sage: alpha = phi.Tq_eigenvalue(p)
            sage: phi._find_extraprec(p,M,alpha,True)
            (13, 1, 2, -2)
        """
        q, aq, eisenloss = self._find_aq(p, M, check)
        newM = M + eisenloss
        # We also need to add precision to account for denominators appearing while solving the difference equation.
        eplog = (newM - 1).exact_log(p)
        while eplog < (newM + eplog).exact_log(p):
            eplog = (newM + eplog).exact_log(p)
            verbose("M = %s, newM = %s, eplog=%s" % (M, newM, eplog), level=2)
        newM += eplog

        # We also need to add precision to account for denominators that might be present in self
        s = self.valuation(p)
        if s < 0:
            newM += -s
        return newM, eisenloss, q, aq

    def p_stabilize_and_lift(self, p, M, alpha=None, ap=None,
                             new_base_ring=None,
                             ordinary=True, algorithm='greenberg', eigensymbol=False,
                             check=True):
        """
        `p`-stabilize and lift self

        INPUT:

        - ``p`` -- prime, not dividing the level of self

        - ``M`` -- precision

        - ``alpha`` -- (default: None) the `U_p` eigenvalue, if known

        - ``ap`` -- (default: None) the Hecke eigenvalue at p (before stabilizing), if known

        - ``new_base_ring`` -- (default: None) if specified, force the resulting eigensymbol to take values in the given ring

        - ``ordinary`` -- (default: True) whether to return the ordinary
                          (at ``p``) eigensymbol.

        - ``algorithm`` -- (default: 'greenberg') a string, either 'greenberg'
          or 'stevens', specifying whether to use
          the lifting algorithm of M.Greenberg or that of Pollack--Stevens.
          The latter one solves the difference equation, which is not needed. The
          option to use Pollack--Stevens' algorithm here is just for historical reasons.

        - ``eigensymbol`` -- (default: False) if True, return an overconvergent eigensymbol. Otherwise just perform a naive lift

        - ``check`` -- (default: True) whether to perform extra sanity checks

        OUTPUT:

        `p`-stabilized and lifted version of self.

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: f = E.pollack_stevens_modular_symbol()
            sage: g = f.p_stabilize_and_lift(3,10)  # long time
            sage: g.Tq_eigenvalue(5)                # long time
            1 + O(3^10)
            sage: g.Tq_eigenvalue(7)                # long time
            1 + 2*3 + 2*3^2 + 2*3^3 + 2*3^4 + 2*3^5 + 2*3^6 + 2*3^7 + 2*3^8 + 2*3^9 + O(3^10)
            sage: g.Tq_eigenvalue(3)                # long time
            2 + 3^2 + 2*3^3 + 2*3^4 + 2*3^6 + 3^8 + 2*3^9 + O(3^10)
        """
        if check:
            p = self._get_prime(p, alpha)
        k = self.parent().weight()
        M = ZZ(M)
        # alpha will be the eigenvalue of Up
        M0 = M + 1
        if alpha is None:
            alpha, new_base_ring, newM, eisenloss, q, aq = self._find_alpha(p, k, M0, ap, new_base_ring, ordinary, check)
        if new_base_ring is None:
            new_base_ring = alpha.parent()
        newM, eisenloss, q, aq = self._find_extraprec(p, M0, alpha, check)
        if hasattr(new_base_ring, 'precision_cap') and newM > new_base_ring.precision_cap():
            raise ValueError("Not enough precision in new base ring")

        # Now we can stabilize
        self = self.p_stabilize(p=p, alpha=alpha, ap=ap, M=newM,
                                new_base_ring=new_base_ring, check=check)
        # And use the standard lifting function for eigensymbols
        Phi = self._lift_to_OMS(p, newM, new_base_ring, algorithm)
        Phi = _iterate_Up(Phi, p=p, M=newM, ap=alpha, q=q, aq=aq, check=check)
        Phi = Phi.reduce_precision(M)
        return Phi._normalize(include_zeroth_moment=True)


class PSModularSymbolElement_dist(PSModularSymbolElement):

    def reduce_precision(self, M):
        r"""
        Only hold on to `M` moments of each value of self

        EXAMPLES::

            sage: D = OverconvergentDistributions(0, 5, 10)
            sage: M = PollackStevensModularSymbols(Gamma0(5), coefficients=D)
            sage: f = M(1)
            sage: f.reduce_precision(1)
            Modular symbol of level 5 with values in Space of 5-adic distributions with k=0 action and precision cap 10
        """
        return self.__class__(self._map.reduce_precision(M), self.parent(),
                              construct=True)

    def precision_relative(self):
        r"""
        Return the number of moments of each value of ``self``.

        EXAMPLES::

            sage: D = OverconvergentDistributions(0, 5, 10)
            sage: M = PollackStevensModularSymbols(Gamma0(5), coefficients=D)
            sage: f = M(1)
            sage: f.precision_relative()
            1
        """
        return min(len(a._moments) for a in self._map)

    def specialize(self, new_base_ring=None):
        r"""
        Return the underlying classical symbol of weight `k`.

        Namely, this applies the canonical map `D_k \to Sym^k` to all
        values of ``self``.

        EXAMPLES::

            sage: D = OverconvergentDistributions(0, 5, 10);  M = PollackStevensModularSymbols(Gamma0(5), coefficients=D); M
            Space of overconvergent modular symbols for Congruence Subgroup Gamma0(5) with sign 0
            and values in Space of 5-adic distributions with k=0 action and precision cap 10
            sage: f = M(1)
            sage: f.specialize()
            Modular symbol of level 5 with values in Sym^0 Z_5^2
            sage: f.specialize().values()
            [1 + O(5), 1 + O(5), 1 + O(5)]
            sage: f.values()
            [1 + O(5), 1 + O(5), 1 + O(5)]
            sage: f.specialize().parent()
            Space of modular symbols for Congruence Subgroup Gamma0(5) with sign 0 and values in Sym^0 Z_5^2
            sage: f.specialize().parent().coefficient_module()
            Sym^0 Z_5^2
            sage: f.specialize().parent().coefficient_module().is_symk()
            True
            sage: f.specialize(Qp(5,20))
            Modular symbol of level 5 with values in Sym^0 Q_5^2
        """
        if new_base_ring is None:
            new_base_ring = self.base_ring()
        return self.__class__(self._map.specialize(new_base_ring),
                              self.parent()._specialize_parent_space(new_base_ring), construct=True)

    def padic_lseries(self,*args, **kwds):
        """
        Return the `p`-adic L-series of this modular symbol.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: phi = E.pollack_stevens_modular_symbol()
            sage: L = phi.lift(37, M=6, eigensymbol=True).padic_lseries(); L  # long time
            37-adic L-series of Modular symbol of level 37 with values in Space of 37-adic distributions with k=0 action and precision cap 7
            sage: L.series(2) # long time
            O(37^6) + (4 + 37 + 36*37^2 + 19*37^3 + 21*37^4 + O(37^5))*T + O(T^2)
        """
        from sage.modular.pollack_stevens.padic_lseries import pAdicLseries
        return pAdicLseries(self, *args, **kwds)
