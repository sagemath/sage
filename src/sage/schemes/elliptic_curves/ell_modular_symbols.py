# -*- coding: utf-8 -*-
r"""
Modular symbols attached to an elliptic curve

To an elliptic curves `E` over the rational numbers one can associate
a space of modular symbols of level `N`,
equal to the conductor of `E`; because `E` is known to be modular.
The space is two-dimensional and contains a subspace on which complex
conjugation acts like `+1` and one on which it acts by `-1`.

There are two implementations of modular symbols, one within ``sage``
and the other as part of Cremona's ``eclib``. One can choose here which
one is used.

Associated to `E` there is a canonical generator in each space. They are maps
`[.]^+` and `[.]^{-}`, both `\QQ \to\QQ`. They are normalised such that

.. math::

   [r]^{+} \Omega^{+} + [r]^{-}\Omega^{-}  = \int_{\infty}^r 2\pi i f(z) dz

where `f` is the newform associated to the isogeny class of `E` and
`\Omega^{+}` is the smallest positive period of the Néron differential
of `E` and `\Omega^{-}` is the smallest positive purely imaginary
period. Note that it depends on `E` rather than on its
isogeny class.

The computation of the space provides first generators, but they are not
necessarily correctly normalised. There are two methods that try to
find the correct scaling factor.

Modular symbols are used to compute `p`-adic `L`-functions.

EXAMPLES::

    sage: E = EllipticCurve("19a1")
    sage: m = E.modular_symbol()
    sage: m(0)
    1/3
    sage: m(1/17)
    -2/3
    sage: m2 = E.modular_symbol(-1, implementation="sage")
    sage: m2(0)
    0
    sage: m2(1/5)
    1/2

    sage: V = E.modular_symbol_space()
    sage: V
    Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 2 for Gamma_0(19) of weight 2 with sign 1 over Rational Field
    sage: V.q_eigenform(30)
    q - 2*q^3 - 2*q^4 + 3*q^5 - q^7 + q^9 + 3*q^11 + 4*q^12 - 4*q^13 - 6*q^15 + 4*q^16 - 3*q^17 + q^19 - 6*q^20 + 2*q^21 + 4*q^25 + 4*q^27 + 2*q^28 + 6*q^29 + O(q^30)

For more details on modular symbols consult the following

REFERENCES:

- [MTT] B. Mazur, J. Tate, and J. Teitelbaum,
  On `p`-adic analogues of the conjectures of Birch and
  Swinnerton-Dyer, Inventiones mathematicae 84, (1986), 1-48.

- [Cre] John Cremona, Algorithms for modular elliptic curves,
  Cambridge University Press, 1997.

- [SW] William Stein and Christian Wuthrich, Algorithms for the
  Arithmetic of Elliptic Curves using Iwasawa Theory, Mathematics
  of Computation 82 (2013), 1757-1792.

AUTHORS:

- William Stein (2007): first version

- Chris Wuthrich (2008): add scaling and reference to eclib

"""

#*****************************************************************************
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function

from sage.structure.sage_object import SageObject
from sage.modular.modsym.all import ModularSymbols
from sage.libs.eclib.newforms import ECModularSymbol
from sage.databases.cremona import parse_cremona_label

from sage.arith.all import next_prime, kronecker_symbol, prime_divisors, valuation
from sage.rings.infinity import unsigned_infinity as infinity
from sage.rings.integer import Integer
from sage.modular.cusps import Cusps
from sage.rings.integer_ring import   ZZ
from sage.rings.rational_field import QQ
from sage.misc.all import verbose

from sage.schemes.elliptic_curves.constructor import EllipticCurve

oo = Cusps(infinity)
zero = Integer(0)

def modular_symbol_space(E, sign, base_ring, bound=None):
    r"""
    Creates the space of modular symbols of a given sign over a give base_ring,
    attached to the isogeny class of elliptic curves.

    INPUT:

    - ``E`` - an elliptic curve over `\QQ`
    - ``sign`` - integer, -1, 0, or 1
    - ``base_ring`` - ring
    - ``bound`` - (default: None) maximum number of Hecke operators to
      use to cut out modular symbols factor.  If None, use
      enough to provably get the correct answer.

    OUTPUT: a space of modular symbols

    EXAMPLES::

        sage: import sage.schemes.elliptic_curves.ell_modular_symbols
        sage: E=EllipticCurve('11a1')
        sage: M=sage.schemes.elliptic_curves.ell_modular_symbols.modular_symbol_space(E,-1,GF(37))
        sage: M
        Modular Symbols space of dimension 1 for Gamma_0(11) of weight 2 with sign -1 over Finite Field of size 37

    """
    _sign = int(sign)
    if _sign != sign:
        raise TypeError('sign must be an integer')
    if not (_sign in [-1,0,1]):
        raise TypeError('sign must -1, 0, or 1')
    N = E.conductor()
    M = ModularSymbols(N, sign=sign, base_ring=base_ring)
    if bound is None:
        bound = M.hecke_bound() + 10
    V = M
    p = 2
    while p <= bound and V.dimension() > 1:
        t = V.T(p)
        ap = E.ap(p)
        V = (t - ap).kernel()
        p = next_prime(p)

    return V

class ModularSymbol(SageObject):
    r"""
    A modular symbol attached to an elliptic curve, which is the map
    `\QQ\to \QQ` obtained by sending `r` to the normalized
    symmetrized (or anti-symmetrized) integral from `r` to `\infty`.

    This is as defined in [MTT], but normalized
    to depend on the curve and not only its isogeny class as in [SW].

    See the documentation of ``E.modular_symbol()`` in
    Elliptic curves over the rational numbers
    for help.

    REFERENCES:

    - [MTT] B. Mazur, J. Tate, and J. Teitelbaum,
      On `p`-adic analogues of the conjectures of Birch and
      Swinnerton-Dyer, Inventiones mathematicae 84, (1986), 1-48.

    - [SW] William Stein and Christian Wuthrich, Computations About Tate-Shafarevich Groups
      using Iwasawa theory, preprint 2009.

    """

    def sign(self):
        r"""
        Return the sign of this elliptic curve modular symbol.

        EXAMPLES::

            sage: m = EllipticCurve('11a1').modular_symbol()
            sage: m.sign()
            1
            sage: m = EllipticCurve('11a1').modular_symbol(sign=-1, implementation="sage")
            sage: m.sign()
            -1
        """

        return self._sign

    def elliptic_curve(self):
        r"""
        Return the elliptic curve of this modular symbol.

        EXAMPLES::

            sage: m = EllipticCurve('11a1').modular_symbol()
            sage: m.elliptic_curve()
            Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field

        """
        return self._E

    def base_ring(self):
        r"""
        Return the base ring for this modular symbol.

        EXAMPLES::

            sage: m = EllipticCurve('11a1').modular_symbol()
            sage: m.base_ring()
            Rational Field
        """
        return self._base_ring

    def _repr_(self):
        r"""
        String representation of modular symbols.

        EXAMPLES::

            sage: m = EllipticCurve('11a1').modular_symbol()
            sage: m
            Modular symbol with sign 1 over Rational Field attached to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
            sage: m = EllipticCurve('43a1').modular_symbol(sign=-1, implementation="sage")
            sage: m
            Modular symbol with sign -1 over Rational Field attached to Elliptic Curve defined by y^2 + y = x^3 + x^2 over Rational Field
        """
        return "Modular symbol with sign %s over %s attached to %s"%(
            self._sign, self._base_ring, self._E)

    def _find_scaling_L_ratio(self):
        r"""
        This function is use to set ``_scaling``, the factor used to adjust the
        scalar multiple of the modular symbol.
        If `[0]`, the modular symbol evaluated at 0, is non-zero, we can just scale
        it with respect to the approximation of the L-value. It is known that
        the quotient is a rational number with small denominator.
        Otherwise we try to scale using quadratic twists.

        ``_scaling`` will be set to a rational non-zero multiple if we succeed and to 1 otherwise.
        Even if we fail we scale at least to make up the difference between the periods
        of the `X_0`-optimal curve and our given curve `E` in the isogeny class.

        EXAMPLES::

            sage : m = EllipticCurve('11a1').modular_symbol()
            sage : m._scaling
            1
            sage: m = EllipticCurve('11a2').modular_symbol()
            sage: m._scaling
            5/2
            sage: m = EllipticCurve('11a3').modular_symbol()
            sage: m._scaling
            1/10
            sage: m = EllipticCurve('11a1').modular_symbol(implementation="sage")
            sage: m._scaling
            1/5
            sage: m = EllipticCurve('11a2').modular_symbol(implementation="sage")
            sage: m._scaling
            1
            sage: m = EllipticCurve('11a3').modular_symbol(implementation="sage")
            sage: m._scaling
            1/25
            sage: m = EllipticCurve('37a1').modular_symbol(implementation="sage")
            sage: m._scaling
            1
            sage: m = EllipticCurve('37a1').modular_symbol()
            sage: m._scaling
            -1
            sage: m = EllipticCurve('389a1').modular_symbol()
            sage: m._scaling
            -1/2
            sage: m = EllipticCurve('389a1').modular_symbol(implementation="sage")
            sage: m._scaling
            2
            sage: m = EllipticCurve('196a1').modular_symbol(implementation="sage")
            sage: m._scaling
            1/2

        Some harder cases fail::

            sage: m = EllipticCurve('121b1').modular_symbol(implementation="sage")
            Warning : Could not normalize the modular symbols, maybe all further results will be multiplied by -1, 2 or -2.
            sage: m._scaling
            1

        TESTS::

            sage: rk0 = ['11a1', '11a2', '15a1', '27a1', '37b1']
            sage: for la in rk0:  # long time (3s on sage.math, 2011)
            ....:          E = EllipticCurve(la)
            ....:          me = E.modular_symbol(implementation="eclib")
            ....:          ms = E.modular_symbol(implementation="sage")
            ....:          print("{} {} {}".format(E.lseries().L_ratio()*E.real_components(), me(0), ms(0)))
            1/5 1/5 1/5
            1 1 1
            1/4 1/4 1/4
            1/3 1/3 1/3
            2/3 2/3 2/3

            sage: rk1 = ['37a1','43a1','53a1', '91b1','91b2','91b3']
            sage: [EllipticCurve(la).modular_symbol()(0) for la in rk1]  # long time (1s on sage.math, 2011)
            [0, 0, 0, 0, 0, 0]
            sage: for la in rk1:  # long time (8s on sage.math, 2011)
            ....:       E = EllipticCurve(la)
            ....:       m = E.modular_symbol()
            ....:       lp = E.padic_lseries(5)
            ....:       for D in [5,17,12,8]:
            ....:           ED = E.quadratic_twist(D)
            ....:           md = sum([kronecker(D,u)*m(ZZ(u)/D) for u in range(D)])
            ....:           etaD = lp._quotient_of_periods_to_twist(D)
            ....:           assert ED.lseries().L_ratio()*ED.real_components() * etaD == md

        """
        E = self._E
        self._scaling = 1 # by now.
        self._failed_to_scale = False

        if self._sign == 1 :
            at0 = self(0)
            if at0 != 0 :
                l1 = self.__lalg__(1)
                if at0 != l1:
                    verbose('scale modular symbols by %s'%(l1/at0))
                    self._scaling = l1/at0
            else :
                # if [0] = 0, we can still hope to scale it correctly by considering twists of E
                Dlist = [5,8,12,13,17,21,24,28,29, 33, 37, 40, 41, 44, 53, 56, 57, 60, 61, 65, 69, 73, 76, 77, 85, 88, 89, 92, 93, 97]  # a list of positive fundamental discriminants
                j = 0
                at0 = 0
                # computes [0]+ for the twist of E by D until one value is non-zero
                while j < 30 and at0 == 0 :
                    D = Dlist[j]
                    # the following line checks if the twist of the newform of E by D is a newform
                    # this is to avoid that we 'twist back'
                    if all( valuation(E.conductor(),ell)<= valuation(D,ell) for ell in prime_divisors(D) ) :
                        at0 = sum([kronecker_symbol(D,u) * self(ZZ(u)/D) for u in range(1,abs(D))])
                    j += 1
                if j == 30 and at0 == 0: # curves like "121b1", "225a1", "225e1", "256a1", "256b1", "289a1", "361a1", "400a1", "400c1", "400h1", "441b1", "441c1", "441d1", "441f1 .. will arrive here
                    self._failed_to_scale = True
                    self.__scale_by_periods_only__()
                else :
                    l1 = self.__lalg__(D)
                    if at0 != l1:
                        verbose('scale modular symbols by %s found at D=%s '%(l1/at0,D), level=2)
                        self._scaling = l1/at0

        else : # that is when sign = -1
            Dlist = [-3,-4,-7,-8,-11,-15,-19,-20,-23,-24, -31, -35, -39, -40, -43, -47, -51, -52, -55, -56, -59, -67, -68, -71, -79, -83, -84, -87, -88, -91]  # a list of negative fundamental discriminants
            j = 0
            at0 = 0
            while j < 30 and at0 == 0 :
                # computes [0]+ for the twist of E by D until one value is non-zero
                D = Dlist[j]
                if all( valuation(E.conductor(),ell)<= valuation(D,ell) for ell in prime_divisors(D) ) :
                    at0 = - sum([kronecker_symbol(D,u) * self(ZZ(u)/D) for u in range(1,abs(D))])
                j += 1
            if j == 30 and at0 == 0: # no more hope for a normalization
                # we do at least a scaling with the quotient of the periods
                self._failed_to_scale = True
                self.__scale_by_periods_only__()
            else :
                l1 = self.__lalg__(D)
                if at0 != l1:
                    verbose('scale modular symbols by %s'%(l1/at0))
                    self._scaling = l1/at0


    def __lalg__(self,D):
        r"""
        For positive `D`, this function evaluates the quotient
        `L(E_D,1)\cdot \sqrt(D)/\Omega_E` where `E_D` is the twist of
        `E` by `D`, `\Omega_E` is the least positive period of `E`.
        For negative `E`, it is the quotient
        `L(E_D,1)\cdot \sqrt(-D)/\Omega^{-}_E`
        where `\Omega^{-}_E` is the least positive imaginary part of a
        non-real period of `E`.

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: m = E.modular_symbol(sign=+1)
            sage: m.__lalg__(1)
            1/5
            sage: m.__lalg__(3)
            5/2

        """
        from sage.functions.all import sqrt
        # the computation of the L-value could take a lot of time,
        # but then the conductor is so large
        # that the computation of modular symbols for E took even longer

        E = self._E
        ED = E.quadratic_twist(D)
        lv = ED.lseries().L_ratio() # this is L(ED,1) divided by the Neron period omD of ED
        lv *= ED.real_components() # now it is by the least positive period
        omD = ED.period_lattice().basis()[0]
        if D > 0 :
            om = E.period_lattice().basis()[0]
            q = sqrt(D)*omD/om * 8
        else :
            om = E.period_lattice().basis()[1].imag()
            if E.real_components() == 1:
                om *= 2
            q = sqrt(-D)*omD/om*8

        # see padic_lseries.pAdicLeries._quotient_of_periods_to_twist
        # for the explanation of the second factor
        verbose('real approximation is %s'%q)
        return lv/8 * QQ(int(round(q)))

    def __scale_by_periods_only__(self):
        r"""
        If we fail to scale with ``_find_scaling_L_ratio``, we drop here
        to try and find the scaling by the quotient of the
        periods to the `X_0`-optimal curve. The resulting ``_scaling``
        is not guaranteed to be correct, but could well be.

        EXAMPLES::

            sage: E = EllipticCurve('19a1')
            sage: m = E.modular_symbol(sign=+1)
            sage: m.__scale_by_periods_only__()
            Warning : Could not normalize the modular symbols, maybe all further results will be multiplied by -1, 2 or -2.
            sage: m._scaling
            1

            sage: E = EllipticCurve('19a2')
            sage: m = E.modular_symbol(sign=+1)
            sage: m._scaling
            3/2
            sage: m.__scale_by_periods_only__()
            Warning : Could not normalize the modular symbols, maybe all further results will be multiplied by -1, 2 or -2.
            sage: m._scaling
            3
        """
        # we only do this inside the cremona-tables.
        try :
            crla = parse_cremona_label(self._E.label())
        except RuntimeError: # raised when curve is outside of the table
            print("Warning : Could not normalize the modular symbols, maybe all further results will be multiplied by a rational number.")
            self._scaling = 1
        else :
            print("Warning : Could not normalize the modular symbols, maybe all further results will be multiplied by -1, 2 or -2.")
            cr0 = Integer(crla[0]).str() + crla[1] + '1'
            E0 = EllipticCurve(cr0)
            if self._sign == 1:
                q = E0.period_lattice().basis()[0]/self._E.period_lattice().basis()[0]
            else:
                q = E0.period_lattice().basis()[1].imag()/self._E.period_lattice().basis()[1].imag()
                if E0.real_components() == 1:
                    q *= 2
                if self._E.real_components() == 1:
                    q /= 2
            q = ZZ(int(round(q*200)))/200
            verbose('scale modular symbols by %s'%q)
            self._scaling = q


class ModularSymbolECLIB(ModularSymbol):
    def __init__(self, E, sign, normalize="L_ratio"):
        r"""
        Modular symbols attached to `E` using ``eclib``.

        INPUT:

        - ``E`` - an elliptic curve
        - ``sign`` - an integer, -1 or 1
        - ``normalize`` - either 'L_ratio' (default) or 'none';
          For 'L_ratio', the modular symbol is correctly normalized
          by comparing it to the quotient of `L(E,1)` by the least
          positive period for the curve and some small twists.
          For 'none', the modular symbol is almost certainly
          not correctly normalized, i.e. all values will be a
          fixed scalar multiple of what they should be.

        EXAMPLES::

            sage: import sage.schemes.elliptic_curves.ell_modular_symbols
            sage: E=EllipticCurve('11a1')
            sage: M=sage.schemes.elliptic_curves.ell_modular_symbols.ModularSymbolECLIB(E,+1)
            sage: M
            Modular symbol with sign 1 over Rational Field attached to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
            sage: M(0)
            1/5
            sage: E=EllipticCurve('11a2')
            sage: M=sage.schemes.elliptic_curves.ell_modular_symbols.ModularSymbolECLIB(E,+1)
            sage: M(0)
            1

        This is a rank 1 case with vanishing positive twists.
        The modular symbol can not be adjusted::

            sage: E=EllipticCurve('121b1')
            sage: M=sage.schemes.elliptic_curves.ell_modular_symbols.ModularSymbolECLIB(E,+1)
            Warning : Could not normalize the modular symbols, maybe all further results will be multiplied by -1, 2 or -2.
            sage: M(0)
            0
            sage: M(1/7)
            -2

            sage: M = EllipticCurve('121d1').modular_symbol(implementation="eclib")
            sage: M(0)
            2
            sage: M = EllipticCurve('121d1').modular_symbol(implementation="eclib",normalize='none')
            sage: M(0)
            8

            sage: E = EllipticCurve('15a1')
            sage: [C.modular_symbol(implementation="eclib",normalize='L_ratio')(0) for C in E.isogeny_class()]
            [1/4, 1/8, 1/4, 1/2, 1/8, 1/16, 1/2, 1]
            sage: [C.modular_symbol(implementation="eclib",normalize='none')(0) for C in E.isogeny_class()]
            [1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4]

        Currently, the interface for negative modular symbols in eclib is not yet written::

            sage: E.modular_symbol(implementation="eclib",sign=-1)
            Traceback (most recent call last):
            ...
            NotImplementedError: Despite that eclib has now -1 modular symbols the interface to them is not yet written.

        TESTS (for :trac:`10236`)::

            sage: E = EllipticCurve('11a1')
            sage: m = E.modular_symbol(implementation="eclib")
            sage: m(1/7)
            7/10
            sage: m(0)
            1/5
        """
        self._sign = ZZ(sign)
        if self._sign != sign:
            raise TypeError('sign must be an integer')
        if self._sign != -1 and self._sign != 1:
            raise TypeError('sign must -1 or 1')
        if self._sign == -1:
            raise NotImplementedError("Despite that eclib has now -1 modular symbols the interface to them is not yet written.")
        self._E = E
        self._use_eclib = True
        self._base_ring = QQ
        self._normalize = normalize
        self._modsym = ECModularSymbol(E)
        p = ZZ(2)
        while not E.is_good(p):
            p = p.next_prime()
        # this computes {0,oo} using the Hecke-operator at p
        self._atzero = sum([self._modsym(ZZ(a)/p) for a in range(p)])/E.Np(p)

        if normalize == "L_ratio":
            self._find_scaling_L_ratio()
        elif normalize == "none":
            self._scaling = ZZ(1)
        else :
            raise ValueError("no normalization '%s' known for modular symbols using John Cremona's eclib"%normalize)


    def _call_with_caching(self, r):
        r"""
        Evaluates the modular symbol at `r`, caching the computed value.

        EXAMPLES::

            sage: m = EllipticCurve('11a1').modular_symbol(implementation="eclib")
            sage: m._call_with_caching(0)
            1/5
        """
        try:
            return self.__cache[r]
        except AttributeError:
            self.__cache = {}
        except KeyError:
            pass
        c = (self._atzero - self._modsym(r))*self._scaling
        self.__cache[r] = c
        return c

    def __call__(self, r):
        r"""
        Evaluates the modular symbol at `r`.

        EXAMPLES::

            sage: m = EllipticCurve('11a1').modular_symbol(implementation="eclib")
            sage: m(0)
            1/5

        """
        # this computes {0,oo} - {0,r} = {r,oo}
        from sage.rings.rational import Rational
        if r != oo:
            r = Rational(r)
            r = r.numer() % r.denom() / r.denom()
        return (self._atzero - self._modsym(r))*self._scaling


class ModularSymbolSage(ModularSymbol):
    def __init__(self, E, sign, normalize="L_ratio"):
        """
        Modular symbols attached to `E` using ``sage``.

        INPUT:

        - ``E`` -- an elliptic curve
        - ``sign`` -- an integer, -1 or 1
        - ``normalize`` -- either 'L_ratio' (default), 'period', or 'none';
          For 'L_ratio', the modular symbol is correctly normalized
          by comparing it to the quotient of `L(E,1)` by the least
          positive period for the curve and some small twists.
          The normalization 'period' uses the integral_period_map
          for modular symbols and is known to be equal to the above
          normalization up to the sign and a possible power of 2.
          For 'none', the modular symbol is almost certainly
          not correctly normalized, i.e. all values will be a
          fixed scalar multiple of what they should be.  But
          the initial computation of the modular symbol is
          much faster, though evaluation of
          it after computing it won't be any faster.

        EXAMPLES::

            sage: E=EllipticCurve('11a1')
            sage: import sage.schemes.elliptic_curves.ell_modular_symbols
            sage: M=sage.schemes.elliptic_curves.ell_modular_symbols.ModularSymbolSage(E,+1)
            sage: M
            Modular symbol with sign 1 over Rational Field attached to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
            sage: M(0)
            1/5
            sage: E=EllipticCurve('11a2')
            sage: M=sage.schemes.elliptic_curves.ell_modular_symbols.ModularSymbolSage(E,+1)
            sage: M(0)
            1
            sage: M=sage.schemes.elliptic_curves.ell_modular_symbols.ModularSymbolSage(E,-1)
            sage: M(1/3)
            1/2

        This is a rank 1 case with vanishing positive twists.
        The modular symbol is adjusted by -2::

            sage: E=EllipticCurve('121b1')
            sage: M=sage.schemes.elliptic_curves.ell_modular_symbols.ModularSymbolSage(E,-1,normalize='L_ratio')
            sage: M(1/3)
            1
            sage: M._scaling
            -1

            sage: M = EllipticCurve('121d1').modular_symbol(implementation="sage")
            sage: M(0)
            2
            sage: M = EllipticCurve('121d1').modular_symbol(implementation="sage", normalize='none')
            sage: M(0)
            1

            sage: E = EllipticCurve('15a1')
            sage: [C.modular_symbol(implementation="sage", normalize='L_ratio')(0) for C in E.isogeny_class()]
            [1/4, 1/8, 1/4, 1/2, 1/8, 1/16, 1/2, 1]
            sage: [C.modular_symbol(implementation="sage", normalize='period')(0) for C in E.isogeny_class()]
            [1/8, 1/16, 1/8, 1/4, 1/16, 1/32, 1/4, 1/2]
            sage: [C.modular_symbol(implementation="sage", normalize='none')(0) for C in E.isogeny_class()]
            [1, 1, 1, 1, 1, 1, 1, 1]

        """
        self._sign = ZZ(sign)
        if self._sign != sign:
            raise TypeError('sign must be an integer')
        if self._sign != -1 and self._sign != 1:
            raise TypeError('sign must -1 or 1')
        self._E = E
        self._use_eclib = False
        self._normalize = normalize
        self._modsym = E.modular_symbol_space(sign=self._sign)
        self._base_ring = self._modsym.base_ring()
        self._ambient_modsym = self._modsym.ambient_module()

        if normalize == "L_ratio":
            self._e = self._modsym.dual_eigenvector()
            self._find_scaling_L_ratio()
            if self._failed_to_scale:
                self._find_scaling_period()  # will reset _e and _scaling
            else:
                self._e  *= self._scaling
        elif normalize == "period" :
            self._find_scaling_period()      # this will set _e and _scaling
        elif normalize == "none":
            self._scaling = 1
            self._e = self._modsym.dual_eigenvector()
        else :
            raise ValueError("no normalization %s known for modular symbols"%normalize)


    def _find_scaling_period(self):
        r"""
        Uses the integral period map of the modular symbol implementation in sage
        in order to determine the scaling. The resulting modular symbol is correct
        only for the `X_0`-optimal curve, at least up to a possible factor +-1 or +-2.

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: m = sage.schemes.elliptic_curves.ell_modular_symbols.ModularSymbolSage(E,+1,normalize='period')
            sage: m._e
            (1/5, 1)
            sage: E = EllipticCurve('11a2')
            sage: m = sage.schemes.elliptic_curves.ell_modular_symbols.ModularSymbolSage(E,+1,normalize='period')
            sage: m._e
            (1, 5)
            sage: E = EllipticCurve('121b2')
            sage: m = sage.schemes.elliptic_curves.ell_modular_symbols.ModularSymbolSage(E,+1,normalize='period')
            sage: m._e
            (0, 11/2, 0, 11/2, 11/2, 0, 0, -3, 2, 1/2, -1, 3/2)

        """

        P = self._modsym.integral_period_mapping()
        self._e = P.matrix().transpose().row(0)
        self._e /= 2
        E = self._E
        try :
            crla = parse_cremona_label(E.label())
        except RuntimeError: # raised when curve is outside of the table
            self._scaling = 1
        else :
            cr0 = Integer(crla[0]).str() + crla[1] + '1'
            E0 = EllipticCurve(cr0)
            if self._sign == 1:
                q = E0.period_lattice().basis()[0]/E.period_lattice().basis()[0]
            else:
                q = E0.period_lattice().basis()[1].imag()/E.period_lattice().basis()[1].imag()
                if E0.real_components() == 1:
                    q *= 2
                if E.real_components() == 1:
                    q /= 2
            q = QQ(int(round(q*200)))/200
            verbose('scale modular symbols by %s'%q)
            self._scaling = q
        self._e *= self._scaling

    def _call_with_caching(self, r):
        r"""
        Evaluates the modular symbol at `r`, caching the computed value.

        EXAMPLES::

            sage: m = EllipticCurve('11a1').modular_symbol(implementation="sage")
            sage: m._call_with_caching(0)
            1/5
        """
        try:
            return self.__cache[r]
        except AttributeError:
            self.__cache = {}
        except KeyError:
            pass
        w = self._ambient_modsym([oo,r]).element()
        c = (self._e).dot_product(w)
        self.__cache[r] = c
        return c

    def __call__(self, r):
        r"""
        Evaluates the modular symbol at `r`.

        EXAMPLES::

            sage: m = EllipticCurve('11a1').modular_symbol(implementation="sage")
            sage: m(0)
            1/5

        """
        # this next line takes most of the time
        w = self._ambient_modsym.modular_symbol([zero, oo, Cusps(r)], check=False)

        return (self._e).dot_product(w.element())
