# -*- coding: utf-8 -*-
r"""
Formal groups of elliptic curves

AUTHORS:

- William Stein: original implementations

- David Harvey: improved asymptotics of some methods

- Nick Alexander: separation from ell_generic.py, bugfixes and
  docstrings
"""

from sage.structure.sage_object import SageObject

import sage.misc.misc as misc
import sage.rings.all as rings
from sage.rings.all import O


class EllipticCurveFormalGroup(SageObject):
    r"""
    The formal group associated to an elliptic curve.
    """
    def __init__(self, E):
        """
        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: F = E.formal_group(); F
            Formal Group associated to the Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
            sage: F == loads(dumps(F))
            True
        """
        self.__E = E

    def __eq__(self, other):
        """
        Check whether ``self`` is equal to ``other``.

        TESTS::

            sage: E = EllipticCurve('35a')
            sage: F1 = E.formal_group()
            sage: F2 = E.formal_group()
            sage: F1 == F2
            True
        """
        if not isinstance(other, EllipticCurveFormalGroup):
            return False

        return self.__E == other.__E

    def __ne__(self, other):
        """
        Check whether ``self`` is not equal to ``other``.

        EXAMPLES::

            sage: E = EllipticCurve('35a')
            sage: F1 = E.formal_group()
            sage: F2 = E.formal_group()
            sage: F1 != F2
            False
        """
        return not (self == other)

    def _repr_(self):
        """
        Return a string representation.

        EXAMPLES::

            sage: E = EllipticCurve('43a')
            sage: F = E.formal_group()
            sage: F._repr_()
            'Formal Group associated to the Elliptic Curve defined by y^2 + y = x^3 + x^2 over Rational Field'
        """
        return "Formal Group associated to the %s" % self.__E

    def curve(self):
        r"""
        Return the elliptic curve this formal group is associated to.

        EXAMPLES::

            sage: E = EllipticCurve("37a")
            sage: F = E.formal_group()
            sage: F.curve()
            Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
        """
        return self.__E

    def w(self, prec=20):
        r"""
        Return the formal group power series `w`.

        INPUT:


        -  ``prec`` - integer (default 20)


        OUTPUT: a power series with given precision

        Return the formal power series

        .. MATH::

                w(t) = t^3 + a_1 t^4 + (a_2 + a_1^2) t^5 + \cdots


        to precision `O(t^{prec})` of Proposition IV.1.1 of
        [Sil2009]_. This is the formal expansion of
        `w = -1/y` about the formal parameter `t = -x/y` at `\infty`.

        The result is cached, and a cached version is returned if
        possible.

        .. warning::

           The resulting power series will have precision prec, but
           its parent PowerSeriesRing will have default precision 20
           (or whatever the default default is).

        ALGORITHM: Uses Newton's method to solve the elliptic curve
        equation at the origin. Complexity is roughly `O(M(n))`
        where `n` is the precision and `M(n)` is the time
        required to multiply polynomials of length `n` over the
        coefficient ring of `E`.

        AUTHOR:

        - David Harvey (2006-09-09): modified to use Newton's
          method instead of a recurrence formula.

        EXAMPLES::

            sage: e = EllipticCurve([0, 0, 1, -1, 0])
            sage: e.formal_group().w(10)
             t^3 + t^6 - t^7 + 2*t^9 + O(t^10)

        Check that caching works::

            sage: e = EllipticCurve([3, 2, -4, -2, 5])
            sage: e.formal_group().w(20)
             t^3 + 3*t^4 + 11*t^5 + 35*t^6 + 101*t^7 + 237*t^8 + 312*t^9 - 949*t^10 - 10389*t^11 - 57087*t^12 - 244092*t^13 - 865333*t^14 - 2455206*t^15 - 4366196*t^16 + 6136610*t^17 + 109938783*t^18 + 688672497*t^19 + O(t^20)
            sage: e.formal_group().w(7)
             t^3 + 3*t^4 + 11*t^5 + 35*t^6 + O(t^7)
            sage: e.formal_group().w(35)
             t^3 + 3*t^4 + 11*t^5 + 35*t^6 + 101*t^7 + 237*t^8 + 312*t^9 - 949*t^10 - 10389*t^11 - 57087*t^12 - 244092*t^13 - 865333*t^14 - 2455206*t^15 - 4366196*t^16 + 6136610*t^17 + 109938783*t^18 + 688672497*t^19 + 3219525807*t^20 + 12337076504*t^21 + 38106669615*t^22 + 79452618700*t^23 - 33430470002*t^24 - 1522228110356*t^25 - 10561222329021*t^26 - 52449326572178*t^27 - 211701726058446*t^28 - 693522772940043*t^29 - 1613471639599050*t^30 - 421817906421378*t^31 + 23651687753515182*t^32 + 181817896829144595*t^33 + 950887648021211163*t^34 + O(t^35)
        """
        prec = max(prec, 0)
        k = self.curve().base_ring()

        try:
            # Try cached version
            w = self.__w
            cached_prec = w.prec()
            R = w.parent()
        except AttributeError:
            # No cached version available
            R = rings.PowerSeriesRing(k, "t")
            w = R([k(0), k(0), k(0), k(1)], 4)
            cached_prec = 4
            self.__w = w

        if prec < cached_prec:
            return R(w, prec)

        # We use the following iteration, which doubles the precision
        # at each step:
        #
        #              z^3 - a_3 w^2 - a_4 z w^2 - 2 a_6 w^3
        # w' = -----------------------------------------------------
        #      1 - a_1 z - a_2 z^2 - 2 a_3 w - 2 a_4 z w - 3 a_6 w^2

        a1, a2, a3, a4, a6 = self.curve().ainvs()
        current_prec = cached_prec
        w = w.truncate()   # work with polynomials instead of power series

        numerator_const = w.parent()([0, 0, 0, 1])      # z^3
        denominator_const = w.parent()([1, -a1, -a2])   # 1 - a_1 z - a_2 z^2

        last_prec = 0
        for next_prec in misc.newton_method_sizes(prec):
            if next_prec > current_prec:
                if w.degree() - 1 > last_prec:
                    # Here it's best to throw away some precision to get us
                    # in sync with the sizes recommended by
                    # newton_method_sizes(). This is especially counter-
                    # intuitive when we throw away almost half of our
                    # cached data!

                    # todo: this might not actually be true, depending on
                    # the overhead of truncate(), which is currently very
                    # high e.g. for NTL based polynomials (but this is
                    # slated to be fixed...)

                    w = w.truncate(last_prec)

                w_squared = w.square()
                w_cubed = (w_squared * w).truncate(next_prec)

                numerator = numerator_const                \
                            -  a3 * w_squared              \
                            -  a4 * w_squared.shift(1)     \
                            -  (2*a6) * w_cubed

                denominator = denominator_const           \
                              - (2*a3) * w                \
                              - (2*a4) * w.shift(1)       \
                              - (3*a6) * w_squared

                # todo: this is quite inefficient, because it gets
                # converted to a power series, then the power series
                # inversion works with polynomials again, and then
                # it gets converted *back* to a power series, and
                # then we convert it to a polynomial again! That's four
                # useless conversions!!!

                inverse = ~R(denominator, prec=next_prec)
                inverse = inverse.truncate(next_prec)

                w = (numerator * inverse).truncate(next_prec)

            last_prec = next_prec

        # convert back to power series
        w = R(w, prec)
        self.__w = (prec, w)
        return self.__w[1]

    def x(self, prec=20):
        r"""
        Return the formal series `x(t) = t/w(t)` in terms of the
        local parameter `t = -x/y` at infinity.

        INPUT:


        -  ``prec`` - integer (default 20)


        OUTPUT: a Laurent series with given precision

        Return the formal series

        .. MATH::

                x(t) = t^{-2} - a_1 t^{-1} - a_2 - a_3 t - \cdots


        to precision `O(t^{prec})` of page 113 of [Sil2009]_.

        .. warning::

           The resulting series will have precision prec, but its
           parent PowerSeriesRing will have default precision 20 (or
           whatever the default default is).

        EXAMPLES::

            sage: EllipticCurve([0, 0, 1, -1, 0]).formal_group().x(10)
             t^-2 - t + t^2 - t^4 + 2*t^5 - t^6 - 2*t^7 + 6*t^8 - 6*t^9 + O(t^10)
        """
        prec = max(prec, 0)
        y = self.y(prec)
        t = y.parent().gen()
        return -t*y + O(t**prec)

    def y(self, prec=20):
        r"""
        Return the formal series `y(t) = -1/w(t)` in terms of the
        local parameter `t = -x/y` at infinity.

        INPUT:


        -  ``prec`` - integer (default 20)


        OUTPUT: a Laurent series with given precision

        Return the formal series

        .. MATH::

                y(t) = - t^{-3} + a_1 t^{-2} + a_2 t + a_3 + \cdots


        to precision `O(t^{prec})` of page 113 of [Sil2009]_.

        The result is cached, and a cached version is returned if
        possible.

        .. warning::

           The resulting series will have precision prec, but its
           parent PowerSeriesRing will have default precision 20 (or
           whatever the default default is).

        EXAMPLES::

            sage: EllipticCurve([0, 0, 1, -1, 0]).formal_group().y(10)
             -t^-3 + 1 - t + t^3 - 2*t^4 + t^5 + 2*t^6 - 6*t^7 + 6*t^8 + 3*t^9 + O(t^10)
        """
        prec = max(prec,0)
        try:
            pr, y = self.__y
        except AttributeError:
            pr = -1
        if prec <= pr:
            t = y.parent().gen()
            return y + O(t**prec)
        w = self.w(prec+6) # XXX why 6?
        t = w.parent().gen()
        y = -(w**(-1)) + O(t**prec)
        self.__y = (prec, y)
        return self.__y[1]

    def differential(self, prec=20):
        r"""
        Return the power series `f(t) = 1 + \cdots` such that
        `f(t) dt` is the usual invariant differential
        `dx/(2y + a_1 x + a_3)`.

        INPUT:


        -  ``prec`` - nonnegative integer (default 20), answer
           will be returned `O(t^{\mathrm{prec}})`


        OUTPUT: a power series with given precision

        Return the formal series

        .. MATH::

                f(t) = 1 + a_1 t + ({a_1}^2 + a_2) t^2 + \cdots


        to precision `O(t^{prec})` of page 113 of [Sil2009]_.

        The result is cached, and a cached version is returned if
        possible.

        .. warning::

           The resulting series will have precision prec, but its
           parent PowerSeriesRing will have default precision 20 (or
           whatever the default default is).

        EXAMPLES::

            sage: EllipticCurve([-1, 1/4]).formal_group().differential(15)
             1 - 2*t^4 + 3/4*t^6 + 6*t^8 - 5*t^10 - 305/16*t^12 + 105/4*t^14 + O(t^15)
            sage: EllipticCurve(Integers(53), [-1, 1/4]).formal_group().differential(15)
             1 + 51*t^4 + 14*t^6 + 6*t^8 + 48*t^10 + 24*t^12 + 13*t^14 + O(t^15)

        AUTHOR:

        - David Harvey (2006-09-10): factored out of log
        """
        prec = max(prec,0)
        try:
            cached_prec, omega = self.__omega
        except AttributeError:
            cached_prec = -1
        if prec <= cached_prec:
            return omega.add_bigoh(prec)

        a = self.curve().ainvs()
        x = self.x(prec+1)
        y = self.y(prec+1)
        xprime = x.derivative()
        g = xprime / (2*y + a[0]*x + a[2])
        self.__omega = (prec, g.power_series().add_bigoh(prec))
        return self.__omega[1]

    def log(self, prec=20):
        r"""
        Return the power series `f(t) = t + \cdots` which is an
        isomorphism to the additive formal group.

        Generally this only makes sense in characteristic zero, although
        the terms before `t^p` may work in characteristic
        `p`.

        INPUT:


        -  ``prec`` - nonnegative integer (default 20)


        OUTPUT: a power series with given precision

        EXAMPLES::

            sage: EllipticCurve([-1, 1/4]).formal_group().log(15)
             t - 2/5*t^5 + 3/28*t^7 + 2/3*t^9 - 5/11*t^11 - 305/208*t^13 + O(t^15)

        AUTHORS:

        - David Harvey (2006-09-10): rewrote to use differential
        """
        return self.differential(prec-1).integral().add_bigoh(prec)

    def inverse(self, prec=20):
        r"""
        Return the formal group inverse law i(t), which satisfies F(t, i(t)) = 0.

        INPUT:


        -  ``prec`` - integer (default 20)


        OUTPUT: a power series with given precision

        Return the formal power series

        .. MATH::

                i(t) = - t + a_1 t^2 + \cdots


        to precision `O(t^{prec})` of page 114 of [Sil2009]_.

        The result is cached, and a cached version is returned if
        possible.

        .. warning::

           The resulting power series will have precision prec, but
           its parent PowerSeriesRing will have default precision 20
           (or whatever the default default is).

        EXAMPLES::

            sage: P.<a1, a2, a3, a4, a6> = ZZ[]
            sage: E = EllipticCurve(list(P.gens()))
            sage: i = E.formal_group().inverse(6); i
            -t - a1*t^2 - a1^2*t^3 + (-a1^3 - a3)*t^4 + (-a1^4 - 3*a1*a3)*t^5 + O(t^6)
            sage: F = E.formal_group().group_law(6)
            sage: F(i.parent().gen(), i)
            O(t^6)
        """
        prec = max(prec,0)
        try:
            pr, inv = self.__inverse
        except AttributeError:
            pr = -1
        if prec <= pr:
            t = inv.parent().gen()
            return inv + O(t**prec)
        x = self.x(prec)
        y = self.y(prec)
        a1, _, a3, _, _ = self.curve().ainvs()
        inv = x / ( y + a1*x + a3)          # page 114 of Silverman, AEC I
        inv = inv.power_series().add_bigoh(prec)
        self.__inverse = (prec, inv)
        return inv

    def group_law(self, prec=10):
        r"""
        Return the formal group law.

        INPUT:


        -  ``prec`` - integer (default 10)


        OUTPUT: a power series with given precision in R[['t1','t2']], where
        the curve is defined over R.

        Return the formal power series

        .. MATH::

           F(t_1, t_2) = t_1 + t_2 - a_1 t_1 t_2 - \cdots


        to precision `O(t1,t2)^{prec}` of page 115 of [Sil2009]_.

        The result is cached, and a cached version is returned if possible.


        AUTHORS:

        - Nick Alexander: minor fixes, docstring

        - Francis Clarke (2012-08): modified to use two-variable power series ring

        EXAMPLES::

            sage: e = EllipticCurve([1, 2])
            sage: e.formal_group().group_law(6)
            t1 + t2 - 2*t1^4*t2 - 4*t1^3*t2^2 - 4*t1^2*t2^3 - 2*t1*t2^4 + O(t1, t2)^6

            sage: e = EllipticCurve('14a1')
            sage: ehat = e.formal()
            sage: ehat.group_law(3)
            t1 + t2 - t1*t2 + O(t1, t2)^3
            sage: ehat.group_law(5)
            t1 + t2 - t1*t2 - 2*t1^3*t2 - 3*t1^2*t2^2 - 2*t1*t2^3 + O(t1, t2)^5

            sage: e = EllipticCurve(GF(7), [3, 4])
            sage: ehat = e.formal()
            sage: ehat.group_law(3)
            t1 + t2 + O(t1, t2)^3
            sage: F = ehat.group_law(7); F
            t1 + t2 + t1^4*t2 + 2*t1^3*t2^2 + 2*t1^2*t2^3 + t1*t2^4 + O(t1, t2)^7

        TESTS::

            sage: R.<x,y,z> = GF(7)[[]]
            sage: F(x, ehat.inverse()(x))
            0 + O(x, y, z)^7
            sage: F(x, y) == F(y, x)
            True
            sage: F(x, F(y, z)) == F(F(x, y), z)
            True

        Let's ensure caching with changed precision is working::

            sage: e.formal_group().group_law(4)
            t1 + t2 + O(t1, t2)^4

        Test for :trac:`9646`::

            sage: P.<a1, a2, a3, a4, a6> = PolynomialRing(ZZ, 5)
            sage: E = EllipticCurve(list(P.gens()))
            sage: F = E.formal().group_law(prec=5)
            sage: t1, t2 = F.parent().gens()
            sage: F(t1, 0)
            t1 + O(t1, t2)^5
            sage: F(0, t2)
            t2 + O(t1, t2)^5
            sage: F.coefficients()[t1*t2^2]
            -a2

        """
        prec = max(prec,0)
        if prec <= 0:
            raise ValueError("The precision must be positive.")

        R = rings.PowerSeriesRing(self.curve().base_ring(), 2, 't1,t2')
        t1, t2 = R.gens()

        if prec == 1:
            return R(0)
        elif prec == 2:
            return t1 + t2 - self.curve().a1()*t1*t2

        try:
            pr, F = self.__group_law
            if prec <= pr:
                return F.add_bigoh(prec)
        except AttributeError:
            pass

        w = self.w(prec+1)
        lam = sum([w[n]*sum(t2**m * t1**(n-m-1) for m in range(n)) for n in range(3, prec+1)])
        lam = lam.add_bigoh(prec)
        nu = w(t1) - lam*t1
        a1, a2, a3, a4, a6 = self.curve().ainvs()
        lam2 = lam*lam
        lam3 = lam2*lam
        # note that the following formula differs from the one in Silverman page 119.
        # See trac ticket 9646 for the explanation and justification.
        t3 = -t1 - t2 - \
             (a1*lam + a3*lam2 + a2*nu + 2*a4*lam*nu + 3*a6*lam2*nu)/  \
             (1 + a2*lam + a4*lam2 + a6*lam3)
        inv = self.inverse(prec)

        F = inv(t3).add_bigoh(prec)
        self.__group_law = (prec, F)
        return F

    def mult_by_n(self, n, prec=10):
        r"""
        Return the formal 'multiplication by n' endomorphism `[n]`.

        INPUT:


        -  ``prec`` - integer (default 10)


        OUTPUT: a power series with given precision

        Return the formal power series

        .. MATH::

                                [n](t) = n t + \cdots


        to precision `O(t^{prec})` of Proposition 2.3 of [Sil2009]_.

        .. warning::

           The resulting power series will have precision prec, but
           its parent PowerSeriesRing will have default precision 20
           (or whatever the default default is).

        AUTHORS:

        - Nick Alexander: minor fixes, docstring

        - David Harvey (2007-03): faster algorithm for char 0 field
          case

        - Hamish Ivey-Law (2009-06): double-and-add algorithm for
          non char 0 field case.

        - Tom Boothby (2009-06): slight improvement to double-and-add

        - Francis Clarke (2012-08): adjustments and simplifications using group_law
          code as modified to yield a two-variable power series.

        EXAMPLES::

            sage: e = EllipticCurve([1, 2, 3, 4, 6])
            sage: e.formal_group().mult_by_n(0, 5)
             O(t^5)
            sage: e.formal_group().mult_by_n(1, 5)
             t + O(t^5)

        We verify an identity of low degree::

            sage: none = e.formal_group().mult_by_n(-1, 5)
            sage: two = e.formal_group().mult_by_n(2, 5)
            sage: ntwo = e.formal_group().mult_by_n(-2, 5)
            sage: ntwo - none(two)
             O(t^5)
            sage: ntwo - two(none)
             O(t^5)

        It's quite fast::

            sage: E = EllipticCurve("37a"); F = E.formal_group()
            sage: F.mult_by_n(100, 20)
            100*t - 49999950*t^4 + 3999999960*t^5 + 14285614285800*t^7 - 2999989920000150*t^8 + 133333325333333400*t^9 - 3571378571674999800*t^10 + 1402585362624965454000*t^11 - 146666057066712847999500*t^12 + 5336978000014213190385000*t^13 - 519472790950932256570002000*t^14 + 93851927683683567270392002800*t^15 - 6673787211563812368630730325175*t^16 + 320129060335050875009191524993000*t^17 - 45670288869783478472872833214986000*t^18 + 5302464956134111125466184947310391600*t^19 + O(t^20)

        TESTS::

            sage: F = EllipticCurve(GF(17), [1, 1]).formal_group()
            sage: F.mult_by_n(10, 50)  # long time (13s on sage.math, 2011)
            10*t + 5*t^5 + 7*t^7 + 13*t^9 + t^11 + 16*t^13 + 13*t^15 + 9*t^17 + 16*t^19 + 15*t^23 + 15*t^25 + 2*t^27 + 10*t^29 + 8*t^31 + 15*t^33 + 6*t^35 + 7*t^37 + 9*t^39 + 10*t^41 + 5*t^43 + 4*t^45 + 6*t^47 + 13*t^49 + O(t^50)

            sage: F = EllipticCurve(GF(101), [1, 1]).formal_group()
            sage: F.mult_by_n(100, 20)
            100*t + O(t^20)

            sage: P.<a1, a2, a3, a4, a6> = PolynomialRing(ZZ, 5)
            sage: E = EllipticCurve(list(P.gens()))
            sage: E.formal().mult_by_n(2,prec=5)
            2*t - a1*t^2 - 2*a2*t^3 + (a1*a2 - 7*a3)*t^4 + O(t^5)

            sage: E = EllipticCurve(QQ, [1,2,3,4,6])
            sage: E.formal().mult_by_n(2,prec=5)
            2*t - t^2 - 4*t^3 - 19*t^4 + O(t^5)
        """
        if self.curve().base_ring().is_field() and self.curve().base_ring().characteristic() == 0 and n != 0:
            # The following algorithm only works over a field of
            # characteristic zero. I don't know whether something similar
            # can be done over a general ring. It would be nice if it did,
            # since it's much faster than using the formal group law.
            # -- dmharvey

            # Create a "formal point" on the original curve E.
            # Our answer only needs prec-1 coefficients (since lowest term
            # is t^1), and x(t) = t^(-2) + ... and y(t) = t^(-3) + ...,
            # so we only need x(t) mod t^(prec-3) and y(t) mod t^(prec-4)
            x = self.x(prec-3)
            y = self.y(prec-4)
            R = x.parent()    # the Laurent series ring over the base ring
            X = self.curve().change_ring(R)
            P = X(x, y)

            # and multiply it by n, using the group law on E
            Q = n*P

            # express it in terms of the formal parameter
            return -Q[0] / Q[1]


        # Now the general case, not necessarily over a field.

        R = rings.PowerSeriesRing(self.curve().base_ring(), "t")
        t = R.gen()

        if n == 1:
            return t.add_bigoh(prec)

        if n == 0:
            return R(0).add_bigoh(prec)

        if n == -1:
            return R(self.inverse(prec))

        if n < 0:
            return self.inverse(prec)(self.mult_by_n(-n, prec))

        F = self.group_law(prec)

        result = t
        if n < 4:
            for m in range(n - 1):
                result = F(result, t)
            return result

        # Double and add is faster than the naive method when n >= 4.
        g = t
        if n & 1:
            result = g
        else:
            result = 0
        n = n >> 1

        while n > 0:
            g = F(g, g)
            if n & 1:
                result = F(result, g)
            n = n >> 1

        return result

    def sigma(self, prec=10):
        r"""
        Return the Weierstrass sigma function as a formal power series
        solution to the differential equation

        .. MATH::

            \frac{d^2 \log \sigma}{dz^2} = - \wp(z)

        with initial conditions `\sigma(O)=0` and `\sigma'(O)=1`,
        expressed in the variable `t=\log_E(z)` of the formal group.

        INPUT:

        -  ``prec`` - integer (default 10)

        OUTPUT: a power series with given precision

        Other solutions can be obtained by multiplication with
        a function of the form `\exp(c z^2)`.
        If the curve has good ordinary reduction at a prime `p`
        then there is a canonical choice of `c` that produces
        the canonical `p`-adic sigma function.
        To obtain that,  please use ``E.padic_sigma(p)`` instead.
        See :meth:`~sage.schemes.elliptic_curves.ell_rational_field.EllipticCurve_rational_field.padic_sigma`

        EXAMPLES::

            sage: E = EllipticCurve('14a')
            sage: F = E.formal_group()
            sage: F.sigma(5)
            t + 1/2*t^2 + 1/3*t^3 + 3/4*t^4 + O(t^5)
        """
        a1,a2,a3,a4,a6 = self.curve().ainvs()

        k = self.curve().base_ring()
        fl = self.log(prec)
        F = fl.reverse()

        S = rings.LaurentSeriesRing(k,'z')
        z = S.gen()
        F = F(z + O(z**prec))
        wp = self.x()(F) + (a1**2 + 4*a2)/12
        g = (1/z**2 - wp).power_series()
        h = g.integral().integral()
        sigma_of_z = z.power_series() * h.exp()

        T = rings.PowerSeriesRing(k,'t')
        fl = fl(T.gen()+O(T.gen()**prec))
        sigma_of_t = sigma_of_z(fl)
        return sigma_of_t
