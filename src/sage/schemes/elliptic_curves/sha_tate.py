# sage_setup: distribution = sagemath-schemes
r"""
Tate-Shafarevich group

If `E` is an elliptic curve over a global field `K`, the Tate-Shafarevich
group is the subgroup of elements in `H^1(K,E)` which map to zero under every
global-to-local restriction map `H^1(K,E) \to H^1(K_v,E)`, one for each place
`v` of `K`.

The group is usually denoted by the Russian letter Sha (Ш), in this document
it will be denoted by `Sha`.

`Sha` is known to be an abelian torsion group. It is conjectured that the
Tate-Shafarevich group is finite for any elliptic curve over a global field.
But it is not known in general.

A theorem of Kolyvagin and Gross-Zagier using Heegner points shows that if the
`L`-series of an elliptic curve `E/\QQ` does not vanish at 1 or has a simple
zero there, then `Sha` is finite.

A theorem of Kato, together with theorems from Iwasawa theory, allows for
certain primes `p` to show that the `p`-primary part of `Sha` is finite and
gives an effective upper bound for it.

The (`p`-adic) conjecture of Birch and Swinnerton-Dyer predicts the order of
`Sha` from the leading term of the (`p`-adic) `L`-series of the elliptic curve.

Sage can compute a few things about `Sha`. The commands ``an``,
``an_numerical`` and ``an_padic`` compute the conjectural order of `Sha` as a
real or `p`-adic number. With ``p_primary_bound`` one can find an upper bound
of the size of the `p`-primary part of `Sha`. Finally, if the analytic rank is
at most 1, then ``bound_kato`` and ``bound_kolyvagin`` find all primes for
which the theorems of Kato and Kolyvagin respectively do not prove the
triviality the `p`-primary part of `Sha`.

EXAMPLES::

    sage: E = EllipticCurve('11a1')
    sage: S = E.sha()
    sage: S.bound_kato()
    [2]
    sage: S.bound_kolyvagin()
    ([2, 5], 1)
    sage: S.an_padic(7,3)
    1 + O(7^5)
    sage: S.an()
    1
    sage: S.an_numerical()
    1.00000000000000

    sage: E = EllipticCurve('389a')
    sage: S = E.sha(); S
    Tate-Shafarevich group for the
     Elliptic Curve defined by y^2 + y = x^3 + x^2 - 2*x over Rational Field
    sage: S.an_numerical()
    1.00000000000000
    sage: S.p_primary_bound(5)
    0
    sage: S.an_padic(5)
    1 + O(5)
    sage: S.an_padic(5,prec=4)  # long time (2s on sage.math, 2011)
    1 + O(5^3)

AUTHORS:

- William Stein (2007) -- initial version

- Chris Wuthrich (April 2009) -- reformat docstrings

- Aly Deines, Chris Wuthrich, Jeaninne Van Order (2016-03): Added
  functionality that tests the Skinner-Urban condition.
"""
# ****************************************************************************
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from math import sqrt

import sage.arith.all as arith

from sage.misc.lazy_import import lazy_import
from sage.misc.verbose import verbose
from sage.modules.free_module_element import vector
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import Q
from sage.rings.real_mpfi import RIF
from sage.rings.real_mpfr import RealField
from sage.structure.sage_object import SageObject

lazy_import("sage.functions.log", "log")
lazy_import("sage.rings.padics.factory", "Qp")

factor = arith.factor
valuation = arith.valuation


class Sha(SageObject):
    r"""
    The Tate-Shafarevich group associated to an elliptic curve.

    If `E` is an elliptic curve over a global field `K`, the Tate-Shafarevich
    group is the subgroup of elements in `H^1(K,E)` which map to zero under
    every global-to-local restriction map `H^1(K,E) \to H^1(K_v,E)`, one for
    each place `v` of `K`.

    EXAMPLES::

        sage: E = EllipticCurve('571a1')
        sage: E._set_gens([])   # curve has rank 0, but non-trivial Sha[2]
        sage: S = E.sha()
        sage: S.bound_kato()
        [2]
        sage: S.bound_kolyvagin()
        ([2], 1)
        sage: S.an_padic(7,3)
        4 + O(7^5)
        sage: S.an()
        4
        sage: S.an_numerical()
        4.00000000000000

        sage: E = EllipticCurve('389a')
        sage: S = E.sha(); S
        Tate-Shafarevich group for the
         Elliptic Curve defined by y^2 + y = x^3 + x^2 - 2*x over Rational Field
        sage: S.an_numerical()
        1.00000000000000
        sage: S.p_primary_bound(5)  # long time
        0
        sage: S.an_padic(5)         # long time
        1 + O(5)
        sage: S.an_padic(5,prec=4)  # very long time
        1 + O(5^3)
    """
    def __init__(self, E):
        r"""
        The Tate-Shafarevich group associated to an elliptic curve.

        INPUT:

        - ``E`` -- an elliptic curve over `\QQ`

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: S = E.sha()
            sage: S
            Tate-Shafarevich group for the
             Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field

            sage: S == loads(dumps(S))
            True
        """
        self.E = E
        self.Emin = E.minimal_model() if not E.is_minimal() else E

    def __eq__(self, other):
        r"""
        Compare two Tate-Shafarevich groups by simply comparing the
        elliptic curves.

        EXAMPLES::

            sage: E = EllipticCurve('37a1')
            sage: S = E.sha()
            sage: S == S
            True
        """
        if not isinstance(other, Sha):
            return False
        return self.E == other.E

    def __ne__(self, other):
        """
        Check whether ``self`` is not equal to ``other``.

        EXAMPLES::

            sage: E = EllipticCurve('37a1')
            sage: S = E.sha()
            sage: S != S
            False
        """
        return not (self == other)

    def __repr__(self):
        r"""
        String representation of the Tate-Shafarevich group.

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: S = E.sha()
            sage: S.__repr__()
            'Tate-Shafarevich group for the Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field'
        """
        return "Tate-Shafarevich group for the " + repr(self.E)

    ########################################################################
    # Functions related to the BSD conjecture.
    ########################################################################

    def an_numerical(self, prec=None,
                     use_database=True, proof=None):
        r"""
        Return the numerical analytic order of `Sha`, which is
        a floating point number in all cases.

        INPUT:

        - ``prec`` -- integer (default: 53); bits precision -- used
          for the `L`-series computation, period,  regulator, etc.
        - ``use_database`` -- whether the rank and generators should
          be looked up in the database if possible. Default is ``True``
        - ``proof`` -- boolean or ``None`` (default: ``None``, see proof.[tab] or
          sage.structure.proof) proof option passed
          onto regulator and rank computation.

        .. NOTE::

            See also the :meth:`an` command, which will return a
            provably correct integer when the rank is 0 or 1.

        .. WARNING::

            If the curve's generators are not known, computing
            them may be very time-consuming.  Also, computation of the
            `L`-series derivative will be time-consuming for large rank and
            large conductor, and the computation time for this may
            increase substantially at greater precision.  However, use of
            very low precision less than about 10 can cause the underlying
            PARI library functions to fail.

        EXAMPLES::

            sage: EllipticCurve('11a').sha().an_numerical()
            1.00000000000000
            sage: EllipticCurve('37a').sha().an_numerical()
            1.00000000000000
            sage: EllipticCurve('389a').sha().an_numerical()
            1.00000000000000
            sage: EllipticCurve('66b3').sha().an_numerical()
            4.00000000000000
            sage: EllipticCurve('5077a').sha().an_numerical()
            1.00000000000000

        A rank 4 curve::

            sage: EllipticCurve([1, -1, 0, -79, 289]).sha().an_numerical()  # long time (3s on sage.math, 2011)
            1.00000000000000

        A rank 5 curve::

            sage: EllipticCurve([0, 0, 1, -79, 342]).sha().an_numerical(prec=10, proof=False)  # long time (22s on sage.math, 2011)
            1.0

        See :issue:`1115`::

            sage: sha = EllipticCurve('37a1').sha()
            sage: [sha.an_numerical(prec) for prec in range(40,100,10)]  # long time (3s on sage.math, 2013)
            [1.0000000000,
             1.0000000000000,
             1.0000000000000000,
             1.0000000000000000000,
             1.0000000000000000000000,
             1.0000000000000000000000000]
        """
        if prec is None:
            prec = RealField().precision()
        RR = RealField(prec)
        prec2 = prec + 2
        RR2 = RealField(prec2)
        try:
            an = self.__an_numerical
            if an.parent().precision() >= prec:
                return RR(an)
            else:  # cached precision too low
                pass
        except AttributeError:
            pass
        # it's critical to switch to the minimal model.
        E = self.Emin
        r = Integer(E.rank(use_database=use_database, proof=proof))
        L = E.lseries().dokchitser(prec=prec2)
        Lr = RR2(L.derivative(1, r))  # L.derivative() returns a Complex
        Om = RR2(E.period_lattice().omega(prec2))
        Reg = E.regulator(use_database=use_database, proof=proof, precision=prec2)
        T = E.torsion_order()
        cp = E.tamagawa_product()
        Sha = RR((Lr * T * T) / (r.factorial() * Om * cp * Reg))
        self.__an_numerical = Sha
        return Sha

    def an(self, use_database=False, descent_second_limit=12):
        r"""
        Return the Birch and Swinnerton-Dyer conjectural order of `Sha`
        as a provably correct integer, unless the analytic rank is > 1,
        in which case this function returns a numerical value.

        INPUT:

        - ``use_database`` -- boolean (default: ``False``); if ``True``, try
          to use any databases installed to lookup the analytic order of
          `Sha`, if possible.  The order of `Sha` is computed if it cannot
          be looked up.

        - ``descent_second_limit`` -- integer (default: 12); limit to use on
          point searching for the quartic twist in the hard case

        This result is proved correct if the order of vanishing is 0
        and the Manin constant is <= 2.

        If the optional parameter ``use_database`` is ``True`` (default:
        ``False``), this function returns the analytic order of `Sha` as
        listed in Cremona's tables, if this curve appears in Cremona's
        tables.

        NOTE:

        If you come across the following error::

            sage: E = EllipticCurve([0, 0, 1, -34874, -2506691])
            sage: E.sha().an()
            Traceback (most recent call last):
            ...
            RuntimeError: Unable to compute the rank, hence generators, with certainty
            (lower bound=0, generators found=[]).  This could be because Sha(E/Q)[2] is
            nontrivial. Try increasing descent_second_limit then trying this command again.

        You can increase the ``descent_second_limit`` (in the above example,
        set to the default, 12) option to try again::

            sage: E.sha().an(descent_second_limit=16)  # long time (2s on sage.math, 2011)
            1

        EXAMPLES::

            sage: E = EllipticCurve([0, -1, 1, -10, -20])   # 11A  = X_0(11)
            sage: E.sha().an()
            1
            sage: E = EllipticCurve([0, -1, 1, 0, 0])       # X_1(11)
            sage: E.sha().an()
            1

            sage: EllipticCurve('14a4').sha().an()
            1
            sage: EllipticCurve('14a4').sha().an(use_database=True)   # will be faster if you have large Cremona database installed
            1

        The smallest conductor curve with nontrivial `Sha`::

            sage: E = EllipticCurve([1,1,1,-352,-2689])     # 66b3
            sage: E.sha().an()
            4

        The four optimal quotients with nontrivial `Sha` and conductor <= 1000::

            sage: E = EllipticCurve([0, -1, 1, -929, -10595])       # 571A
            sage: E.sha().an()
            4
            sage: E = EllipticCurve([1, 1, 0, -1154, -15345])       # 681B
            sage: E.sha().an()
            9
            sage: E = EllipticCurve([0, -1, 0, -900, -10098])       # 960D
            sage: E.sha().an()
            4
            sage: E = EllipticCurve([0, 1, 0, -20, -42])            # 960N
            sage: E.sha().an()
            4

        The smallest conductor curve of rank > 1::

            sage: E = EllipticCurve([0, 1, 1, -2, 0])       # 389A (rank 2)
            sage: E.sha().an()
            1.00000000000000

        The following are examples that require computation of the Mordell-
        Weil group and regulator::

            sage: E = EllipticCurve([0, 0, 1, -1, 0])                     # 37A  (rank 1)
            sage: E.sha().an()
            1

            sage: E = EllipticCurve("1610f3")
            sage: E.sha().an()
            4

        In this case the input curve is not minimal, and if this function did
        not transform it to be minimal, it would give nonsense::

            sage: E = EllipticCurve([0, -432*6^2])
            sage: E.sha().an()
            1

        See :issue:`10096`: this used to give the wrong result 6.0000
        before since the minimal model was not used::

            sage: E = EllipticCurve([1215*1216, 0]) # non-minimal model
            sage: E.sha().an()  # long time (2s on sage.math, 2011)
            1.00000000000000
            sage: E.minimal_model().sha().an()  # long time (1s on sage.math, 2011)
            1.00000000000000
        """
        if hasattr(self, '__an'):
            return self.__an
        if use_database:
            d = self.Emin.database_curve()
            if hasattr(d, 'db_extra'):
                self.__an = Integer(round(float(d.db_extra[4])))
                return self.__an

        # it's critical to switch to the minimal model.
        E = self.Emin
        eps = E.root_number()
        if eps == 1:
            L1_over_omega = E.lseries().L_ratio()
            if L1_over_omega == 0:  # order of vanishing is at least 2
                return self.an_numerical(use_database=use_database)
            T = E.torsion_subgroup().order()
            Sha = (L1_over_omega * T * T) / Q(E.tamagawa_product())
            try:
                Sha = Integer(Sha)
            except ValueError:
                raise RuntimeError("There is a bug in an, since the computed conjectural order of Sha is %s, which is not an integer." % Sha)
            if not arith.is_square(Sha):
                raise RuntimeError("There is a bug in an, since the computed conjectural order of Sha is %s, which is not a square." % Sha)
            E.__an = Sha
            self.__an = Sha
            return Sha

        else:  # rank > 0  (Not provably correct)
            L1, error_bound = E.lseries().deriv_at1(10 * sqrt(E.conductor()) + 10)
            if abs(L1) < error_bound:
                s = self.an_numerical()
                E.__an = s
                self.__an = s
                return s

            regulator = E.regulator(use_database=use_database, descent_second_limit=descent_second_limit)
            T = E.torsion_subgroup().order()
            omega = E.period_lattice().omega()
            Sha = ((L1 * T * T) / (E.tamagawa_product() * regulator * omega)).round()
            try:
                Sha = Integer(Sha)
            except ValueError:
                raise RuntimeError("There is a bug in an, since the computed conjectural order of Sha is %s, which is not an integer." % Sha)
            if not arith.is_square(Sha):
                raise RuntimeError("There is a bug in an, since the computed conjectural order of Sha is %s, which is not a square." % Sha)
            E.__an = Sha
            self.__an = Sha
            return Sha

    def an_padic(self, p, prec=0, use_twists=True):
        r"""
        Return the conjectural order of `Sha(E/\QQ)`,
        according to the `p`-adic analogue of the Birch
        and Swinnerton-Dyer conjecture as formulated
        in [MTT1986]_ and [BP1993]_.

        INPUT:

        - ``p`` -- a prime > 3

        - ``prec`` -- (optional) the precision used in the computation of the
          `p`-adic L-Series

        - ``use_twists`` -- boolean (default: ``True``); if ``True`` the algorithm may
          change to a quadratic twist with minimal conductor to do the modular
          symbol computations rather than using the modular symbols of the
          curve itself. If ``False`` it forces the computation using the
          modular symbols of the curve itself.

        OUTPUT: `p`-adic number - that conjecturally equals `\# Sha(E/\QQ)`.

        If ``prec`` is set to zero (default) then the precision is set so that
        at least the first `p`-adic digit of conjectural `\# Sha(E/\QQ)` is
        determined.

        EXAMPLES:

        Good ordinary examples::

            sage: EllipticCurve('11a1').sha().an_padic(5)    # rank 0
            1 + O(5^22)
            sage: EllipticCurve('43a1').sha().an_padic(5)    # rank 1
            1 + O(5)
            sage: EllipticCurve('389a1').sha().an_padic(5,4) # rank 2, long time (2s on sage.math, 2011)
            1 + O(5^3)
            sage: EllipticCurve('858k2').sha().an_padic(7)   # rank 0, non trivial sha, long time (10s on sage.math, 2011)
            7^2 + O(7^24)
            sage: EllipticCurve('300b2').sha().an_padic(3)   # 9 elements in sha, long time (2s on sage.math, 2011)
            3^2 + O(3^24)
            sage: EllipticCurve('300b2').sha().an_padic(7, prec=6)  # long time
            2 + 7 + O(7^8)

        Exceptional cases::

            sage: EllipticCurve('11a1').sha().an_padic(11) # rank 0
            1 + O(11^22)
            sage: EllipticCurve('130a1').sha().an_padic(5) # rank 1
            1 + O(5)

        Non-split, but rank 0 case (:issue:`7331`)::

            sage: EllipticCurve('270b1').sha().an_padic(5) # rank 0, long time (2s on sage.math, 2011)
            1 + O(5^22)

        The output has the correct sign::

            sage: EllipticCurve('123a1').sha().an_padic(41) # rank 1, long time (3s on sage.math, 2011)
            1 + O(41)

        Supersingular cases::

            sage: EllipticCurve('34a1').sha().an_padic(5) # rank 0
            1 + O(5^22)
            sage: EllipticCurve('53a1').sha().an_padic(5) # rank 1, long time (11s on sage.math, 2011)
            1 + O(5)

        Cases that use a twist to a lower conductor::

            sage: EllipticCurve('99a1').sha().an_padic(5)
            1 + O(5)
            sage: EllipticCurve('240d3').sha().an_padic(5)  # sha has 4 elements here
            4 + O(5)
            sage: EllipticCurve('448c5').sha().an_padic(7, prec=4, use_twists=False)  # long time (2s on sage.math, 2011)
            2 + 7 + O(7^6)
            sage: EllipticCurve([-19,34]).sha().an_padic(5)  # see trac #6455, long time (4s on sage.math, 2011)
            1 + O(5)

        Test for :issue:`15737`::

            sage: E = EllipticCurve([-100,0])
            sage: s = E.sha()
            sage: s.an_padic(13)
            1 + O(13^20)
        """
        try:
            return self.__an_padic[(p, prec)]
        except AttributeError:
            self.__an_padic = {}
        except KeyError:
            pass

        E = self.Emin
        tam = E.tamagawa_product()
        tors = E.torsion_order()**2
        r = E.rank()
        if r > 0:
            reg = E.padic_regulator(p)
        else:
            if E.is_supersingular(p):
                reg = vector([Qp(p, 20)(1), 0])
            else:
                reg = Qp(p, 20)(1)

        if use_twists and p > 2:
            Et, D = E.minimal_quadratic_twist()
            # trac 6455 : we have to assure that the twist back is allowed
            D = ZZ(D)
            if not D % p:
                D = D // p
            for ell in D.prime_divisors():
                if ell % 2:
                    if not Et.conductor() % ell**2:
                        D = D // ell
            ve = valuation(D, 2)
            de = (D >> ve).abs()
            if de % 4 == 3:
                de = -de
            Et = E.quadratic_twist(de)
            # now check individually if we can twist by -1 or 2 or -2
            Nmin = Et.conductor()
            Dmax = de
            for DD in [-4 * de, 8 * de, -8 * de]:
                Et = E.quadratic_twist(DD)
                if Et.conductor() < Nmin and valuation(Et.conductor(), 2) <= valuation(DD, 2):
                    Nmin = Et.conductor()
                    Dmax = DD
            D = Dmax
            Et = E.quadratic_twist(D)
            lp = Et.padic_lseries(p)
        else:
            lp = E.padic_lseries(p)
            D = 1

        if r == 0 and D == 1:
            # short cut for rank 0 curves, we do not
            # to compute the p-adic L-function, the leading
            # term will be the L-value divided by the Neron
            # period.
            ms = E.modular_symbol(sign=+1, normalize='L_ratio')
            lstar = ms(0) / E.real_components()
            bsd = tam / tors
            if prec == 0:
                # prec = valuation(lstar/bsd, p)
                prec = 20
            shan = Qp(p, prec=prec + 2)(lstar / bsd)

        elif E.is_ordinary(p):
            K = reg.parent()
            lg = log(K(1 + p))

            if (E.is_good(p) or E.ap(p) == -1):
                if not E.is_good(p):
                    eps = 2
                else:
                    eps = (1 - arith.kronecker_symbol(D, p) / lp.alpha())**2
                # according to the p-adic BSD this should be equal to the leading term of the p-adic L-series divided by sha:
                bsdp = tam * reg * eps / tors / lg**r
            else:
                r += 1   # exceptional zero
                eq = E.tate_curve(p)
                Li = eq.L_invariant()

                # according to the p-adic BSD (Mazur-Tate-Teitelbaum)
                # this should be equal to the leading term of the p-adic L-series divided by sha:
                bsdp = tam * reg * Li / tors / lg**r

            v = bsdp.valuation()
            if v > 0:
                verbose("the prime is irregular for this curve.")

            # determine how much prec we need to prove at least the
            # triviality of the p-primary part of Sha

            if prec == 0:
                n = max(v, 2)
                bounds = lp._prec_bounds(n, r + 1)
                while bounds[r] <= v:
                    n += 1
                    bounds = lp._prec_bounds(n, r + 1)
                verbose("set precision to %s" % n)
            else:
                n = max(2, prec)

            not_yet_enough_prec = True
            while not_yet_enough_prec:
                lps = lp.series(n, quadratic_twist=D, prec=r + 1)
                lstar = lps[r]
                if (lstar != 0) or (prec != 0):
                    not_yet_enough_prec = False
                else:
                    n += 1
                    verbose("increased precision to %s" % n)

            shan = lstar / bsdp

        elif E.is_supersingular(p):
            K = reg[0].parent()
            lg = log(K(1 + p))

            # according to the p-adic BSD this should be equal to the leading term of the D_p - valued
            # L-series :
            bsdp = tam / tors / lg**r * reg
            # note this is an element in Q_p^2

            verbose("the algebraic leading terms : %s" % bsdp)

            v = [bsdp[0].valuation(), bsdp[1].valuation()]

            if prec == 0:
                n = max(min(v) + 2, 3)
            else:
                n = max(3, prec)

            verbose("...computing the p-adic L-series")
            not_yet_enough_prec = True
            while not_yet_enough_prec:
                lps = lp.Dp_valued_series(n, quadratic_twist=D, prec=r + 1)
                lstar = [lps[0][r], lps[1][r]]
                verbose("the leading terms : %s" % lstar)
                if (lstar[0] != 0 or lstar[1] != 0) or (prec != 0):
                    not_yet_enough_prec = False
                else:
                    n += 1
                    verbose("increased precision to %s" % n)

            verbose("...putting things together")
            if bsdp[0] != 0:
                shan0 = lstar[0] / bsdp[0]
            else:
                shan0 = 0   # this should actually never happen
            if bsdp[1] != 0:
                shan1 = lstar[1] / bsdp[1]
            else:
                shan1 = 0   # this should conjecturally only happen when the rank is 0
            verbose(f"the two values for Sha : {shan0}, {shan1}")

            # check consistency (the first two are only here to avoid a bug in the p-adic L-series
            # (namely the coefficients of zero-relative precision are treated as zero)
            if shan0 != 0 and shan1 != 0 and shan0 - shan1 != 0:
                raise RuntimeError("There must be a bug in the supersingular routines for the p-adic BSD.")

            # take the better
            if shan1 == 0 or shan0.precision_relative() > shan1.precision_relative():
                shan = shan0
            else:
                shan = shan1

        else:
            raise ValueError("The curve has to have semi-stable reduction at p.")

        self.__an_padic[(p, prec)] = shan
        return shan

    def p_primary_order(self, p):
        r"""
        Return the order of the `p`-primary part of the Tate-Shafarevich
        group.

        This uses the result of Skinner and Urban [SU2014]_ on the
        main conjecture in Iwasawa theory. In particular the elliptic
        curve must have good ordinary reduction at `p`, the residual
        Galois representation must be surjective. Furthermore there must
        be an auxiliary prime `\ell` dividing the conductor of the curve
        exactly once such that the residual representation is ramified
        at `p`.

        INPUT:

        - ``p`` -- an odd prime

        OUTPUT:

        - ``e`` -- nonnegative integer such that `p^e` is the
          order of the `p`-primary order if the conditions are satisfied
          and raises a :exc:`ValueError` otherwise.

        EXAMPLES::

            sage: E = EllipticCurve("389a1")  # rank 2
            sage: E.sha().p_primary_order(5)
            0
            sage: E = EllipticCurve("11a1")
            sage: E.sha().p_primary_order(7)
            0
            sage: E.sha().p_primary_order(5)
            Traceback (most recent call last):
            ...
            ValueError: The order is not provably known using Skinner-Urban.
            Try running p_primary_bound to get a bound.
        """
        E = self.E
        # does not work if p = 2
        if p == 2:
            raise ValueError(f"{p} is not an odd prime")
        if (E.is_ordinary(p) and E.conductor() % p and
                E.galois_representation().is_surjective(p)):
            N = E.conductor()
            fac = N.factor()
            # the auxiliary prime will be one dividing the conductor
            if all(E.tate_curve(ell).parameter().valuation() % p == 0
                   for ell, e in fac if e == 1):
                raise ValueError("The order is not provably known using Skinner-Urban.\n" +
                                 "Try running p_primary_bound to get a bound.")
        else:
            raise ValueError("The order is not provably known using Skinner-Urban.\n" +
                             "Try running p_primary_bound to get a bound.")
        return self.p_primary_bound(p)

    def p_primary_bound(self, p):
        r"""
        Return a provable upper bound for the order of the
        `p`-primary part `Sha(E)(p)` of the Tate-Shafarevich group.

        INPUT:

        - ``p`` -- a prime > 2

        OUTPUT:

        - ``e`` -- nonnegative integer such that `p^e` is an upper
          bound for the order of `Sha(E)(p)`

        In particular, if this algorithm does not fail, then it proves
        that the `p`-primary part of `Sha` is finite. This works also
        for curves of rank > 1.

        Note also that this bound is sharp if one assumes the main conjecture
        of Iwasawa theory of elliptic curves. One may use the method
        ``p_primary_order`` for checking if the extra conditions hold under
        which the main conjecture is known by the work of Skinner and Urban.
        This then returns the provable `p`-primary part of the Tate-Shafarevich
        group.

        Currently the algorithm is only implemented when the following
        conditions are verified:

        - The `p`-adic Galois representation must be surjective or
          must have its image contained in a Borel subgroup.

        - The reduction at `p` is not allowed to be additive.

        - If the reduction at `p` is non-split multiplicative, then
          the rank must be 0.

        - If `p = 3`, then the reduction at 3 must be good ordinary or
          split multiplicative, and the rank must be 0.

        ALGORITHM:

        The algorithm is described in [SW2013]_. The results for the
        reducible case can be found in [Wu2004]_. The main ingredient is
        Kato's result on the main conjecture in Iwasawa theory.

        EXAMPLES::

            sage: e = EllipticCurve('11a3')
            sage: e.sha().p_primary_bound(3)
            0
            sage: e.sha().p_primary_bound(5)
            0
            sage: e.sha().p_primary_bound(7)
            0
            sage: e.sha().p_primary_bound(11)
            0
            sage: e.sha().p_primary_bound(13)
            0

            sage: e = EllipticCurve('389a1')
            sage: e.sha().p_primary_bound(5)
            0
            sage: e.sha().p_primary_bound(7)
            0
            sage: e.sha().p_primary_bound(11)
            0
            sage: e.sha().p_primary_bound(13)
            0

            sage: e = EllipticCurve('858k2')
            sage: e.sha().p_primary_bound(3)  # long time (10s on sage.math, 2011)
            0

        Some checks for :issue:`6406` and :issue:`16959`::

            sage: e.sha().p_primary_bound(7)  # long time
            2

            sage: E = EllipticCurve('608b1')
            sage: E.sha().p_primary_bound(5)
            Traceback (most recent call last):
            ...
            ValueError: The p-adic Galois representation is not surjective or reducible.
            Current knowledge about Euler systems does not provide an upper bound
            in this case. Try an_padic for a conjectural bound.

            sage: E.sha().an_padic(5)           # long time
            1 + O(5^22)

            sage: E = EllipticCurve("5040bi1")
            sage: E.sha().p_primary_bound(5)    # long time
            0
        """
        p = Integer(p)
        if p == 2:
            raise ValueError("The prime p must be odd.")
        E = self.Emin
        if not (E.is_ordinary(p) or E.is_good(p)):
            raise ValueError("The curve has to have semi-stable reduction at p.")
        rho = E.galois_representation()
        su = rho.is_surjective(p)
        re = rho.is_reducible(p)
        if not su and not re:
            raise ValueError("The p-adic Galois representation is not surjective or reducible. Current knowledge about Euler systems does not provide an upper bound in this case. Try an_padic for a conjectural bound.")
        shan = self.an_padic(p, prec=0, use_twists=True)
        if shan == 0:
            raise RuntimeError("There is a bug in an_padic.")
        return shan.valuation()

    def two_selmer_bound(self):
        r"""
        Return the 2-rank, i.e. the `\GF{2}`-dimension
        of the 2-torsion part of `Sha`, provided we can determine the
        rank of `E`.

        EXAMPLES::

            sage: sh = EllipticCurve('571a1').sha()
            sage: sh.two_selmer_bound()
            2
            sage: sh.an()
            4

            sage: sh = EllipticCurve('66a1').sha()
            sage: sh.two_selmer_bound()
            0
            sage: sh.an()
            1

            sage: sh = EllipticCurve('960d1').sha()
            sage: sh.two_selmer_bound()
            2
            sage: sh.an()
            4
        """
        E = self.Emin
        S = E.selmer_rank()
        r = E.rank()
        t = E.two_torsion_rank()
        b = S - r - t
        return 0 if b < 0 else b

    def bound_kolyvagin(self, D=0, regulator=None,
                        ignore_nonsurj_hypothesis=False):
        r"""
        Given a fundamental discriminant `D \neq -3,-4` that satisfies the
        Heegner hypothesis for `E`, return a list of primes so that
        Kolyvagin's theorem (as in Gross's paper) implies that any
        prime divisor of `Sha` is in this list.

        INPUT:

        - ``D`` -- (optional) a fundamental discriminant < -4 that satisfies
          the Heegner hypothesis for `E`; if not given, use the first such `D`
        - ``regulator`` -- (optional) regulator of `E(K)`; if not given, will
          be computed (which could take a long time)
        - ``ignore_nonsurj_hypothesis`` -- (default: ``False``);
          if ``True``, then gives the bound coming from Heegner point
          index, but without any hypothesis on surjectivity
          of the mod-`p` representation

        OUTPUT:

        - ``list`` -- list of primes such that if `p` divides `Sha(E/K)`, then
          `p` is in this list, unless `E/K` has complex multiplication or
          analytic rank greater than 2 (in which case we return 0)

        - ``index`` -- the odd part of the index of the Heegner point in the full
          group of `K`-rational points on `E` (if `E` has CM, returns 0)

        REMARKS:

        1)      We do not have to assume that the Manin constant is 1
                (or a power of 2).  If the Manin constant were
                divisible by a prime, that prime would get included in
                the list of bad primes.

        2)      We assume the Gross-Zagier theorem is true under the
                hypothesis that `gcd(N,D) = 1`, instead of the stronger
                hypothesis `gcd(2\cdot N,D)=1` that is in the original
                Gross-Zagier paper.  That Gross-Zagier is true when
                `gcd(N,D)=1` is "well-known" to the experts, but does not
                seem to written up well in the literature.

        3)      Correctness of the computation is guaranteed using
                interval arithmetic, under the assumption that the
                regulator, square root, and period lattice are
                computed to precision at least `10^{-10}`, i.e., they are
                correct up to addition or a real number with absolute
                value less than `10^{-10}`.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: E.sha().bound_kolyvagin()
            ([2], 1)
            sage: E = EllipticCurve('141a')
            sage: E.sha().an()
            1
            sage: E.sha().bound_kolyvagin()
            ([2, 7], 49)

        We get no information when the curve has rank 2.::

            sage: E = EllipticCurve('389a')
            sage: E.sha().bound_kolyvagin()
            (0, 0)
            sage: E = EllipticCurve('681b')
            sage: E.sha().an()
            9
            sage: E.sha().bound_kolyvagin()
            ([2, 3], 9)
        """
        E = self.Emin
        if E.has_cm():
            return 0, 0

        if D == 0:
            D = -5
            while not E.satisfies_heegner_hypothesis(D):
                D -= 1

        if not E.satisfies_heegner_hypothesis(D):
            raise ArithmeticError("Discriminant (=%s) must be a fundamental discriminant that satisfies the Heegner hypothesis." % D)
        if D == -3 or D == -4:
            raise ArithmeticError("Discriminant (=%s) must not be -3 or -4." % D)
        eps = E.root_number()
        L1_vanishes = E.lseries().L1_vanishes()
        if eps == 1 and L1_vanishes:
            return 0, 0        # rank even hence >= 2, so Kolyvagin gives nothing.
        alpha = sqrt(abs(D)) / (2 * E.period_lattice().complex_area())
        F = E.quadratic_twist(D)
        k_E = 2 * sqrt(E.conductor()) + 10
        k_F = 2 * sqrt(F.conductor()) + 10
        # k_E = 2
        # k_F = 2

        MIN_ERR = 1e-10
        # we assume that regulator and
        # discriminant, etc., computed to this accuracy.

        tries = 0
        while True:
            tries += 1
            if tries >= 6:
                raise RuntimeError("Too many precision increases in bound_kolyvagin")
            if eps == 1:   # E has even rank
                verbose("Conductor of twist = %s" % F.conductor())
                LF1, err_F = F.lseries().deriv_at1(k_F)
                LE1, err_E = E.lseries().at1(k_E)
                err_F = max(err_F, MIN_ERR)
                err_E = max(err_E, MIN_ERR)
                if regulator is not None:
                    hZ = regulator / 2
                else:
                    hZ = F.regulator(use_database=True) / 2
                I = RIF(alpha) * RIF(LE1-err_E, LE1+err_E) * RIF(LF1-err_F, LF1+err_F) / RIF(hZ)

            else:          # E has odd rank

                if regulator is not None:
                    hZ = regulator / 2
                else:
                    hZ = E.regulator(use_database=True) / 2
                LE1, err_E = E.lseries().deriv_at1(k_E)
                LF1, err_F = F.lseries().at1(k_F)
                err_F = max(err_F, MIN_ERR)
                err_E = max(err_E, MIN_ERR)
                # I = alpha * LE1 * LF1 / hZ

                I = RIF(alpha) * RIF(LE1-err_E, LE1+err_E) * RIF(LF1-err_F, LF1+err_F) / RIF(hZ)

            verbose('interval = %s' % I)
            t, n = I.is_int()
            if t:
                break
            elif I.absolute_diameter() < 1:
                raise RuntimeError("Problem in bound_kolyvagin; square of index is not an integer -- D={}, I={}.".format(D, I))
            verbose("Doubling bounds")
            k_E *= 2
            k_F *= 2
        # end while

        # We include 2 since Kolyvagin (in Gross) says nothing there
        if n == 0:
            return 0, 0  # no bound
        B = [2]
        for p, e in n.factor():
            if p > 2:
                if e % 2:
                    raise RuntimeError("Problem in bound_kolyvagin; square of index is not a perfect square!  D={}, I={}, n={}, e={}.".format(D, I, n, e))
                B.append(p)
            else:
                n >>= e  # replace n by its odd part
        if not ignore_nonsurj_hypothesis:
            B.extend(E.galois_representation().non_surjective())
        return sorted({int(x) for x in B}), n

    def bound_kato(self):
        r"""
        Return a list of primes `p` such that the theorems of Kato's [Kat2004]_
        and others (e.g., as explained in a thesis of Grigor Grigorov [Gri2005]_)
        imply that if `p` divides  the order of `Sha(E/\QQ)` then `p` is in
        the list.

        If `L(E,1) = 0`, then this function gives no information, so
        it returns ``False``.

        THEOREM: Suppose `L(E,1) \neq 0` and `p \neq 2` is a prime such
        that

        - `E` does not have additive reduction at `p`,
        - either the `p`-adic representation is surjective or has its
          image contained in a Borel subgroup.

        Then `{ord}_p(\#Sha(E))` is bounded from above by the `p`-adic valuation of  `L(E,1)\cdot\#E(\QQ)_{tor}^2 / (\Omega_E \cdot \prod c_v)`.

        If the `L`-series vanishes, the method ``p_primary_bound`` can be used instead.

        EXAMPLES::

            sage: E = EllipticCurve([0, -1, 1, -10, -20])   # 11A  = X_0(11)
            sage: E.sha().bound_kato()
            [2]
            sage: E = EllipticCurve([0, -1, 1, 0, 0])       # X_1(11)
            sage: E.sha().bound_kato()
            [2]
            sage: E = EllipticCurve([1,1,1,-352,-2689])     # 66B3
            sage: E.sha().bound_kato()
            [2]

        For the following curve one really has that 25 divides the
        order of `Sha` (by [GJPST2009]_)::

            sage: E = EllipticCurve([1, -1, 0, -332311, -73733731])   # 1058D1
            sage: E.sha().bound_kato()                        # long time (about 1 second)
            [2, 5, 23]
            sage: E.galois_representation().non_surjective()  # long time (about 1 second)
            []

        For this one, `Sha` is divisible by 7::

            sage: E = EllipticCurve([0, 0, 0, -4062871, -3152083138])   # 3364C1
            sage: E.sha().bound_kato()                        # long time (< 10 seconds)
            [2, 7, 29]

        No information about curves of rank > 0::

            sage: E = EllipticCurve([0, 0, 1, -1, 0])       # 37A  (rank 1)
            sage: E.sha().bound_kato()
            False
        """
        E = self.Emin
        if E.has_cm():
            return False
        if E.lseries().L1_vanishes():
            return False
        B = [2]
        rho = E.galois_representation()
        B.extend(p for p in rho.non_surjective()
                 if p > 2 and p not in rho.reducible_primes())
        B.extend(p for p in E.conductor().prime_divisors()
                 if E.has_additive_reduction(p))

        # The only other p that might divide B are those that divide
        # the integer 2*#E(Q)_tor^2 * L(E,1)/omega.  So we compute
        # that to sufficient precision to determine it.  Note that
        # we have to assume the Manin constant is <=2 in order to provably
        # compute L(E,1)/omega.
        for p, n in self.an().factor():
            if n >= 2:    # use parity of Sha
                B.append(int(p))
        return sorted(set(B))

    def bound(self):
        r"""
        Compute a provably correct bound on the order of the Tate-Shafarevich
        group of this curve.

        The bound is either ``False`` (no bound) or a list ``B`` of primes
        such that any prime divisor of the order of `Sha` is in this list.

        EXAMPLES::

            sage: EllipticCurve('37a').sha().bound()
            ([2], 1)
        """
        if self.Emin.lseries().L1_vanishes():
            return self.bound_kolyvagin()
        return self.bound_kato()
