# sage_setup: distribution = sagemath-schemes
r"""
Modular parametrization of elliptic curves over `\QQ`

By the work of Taylor--Wiles et al. it is known that there
is a surjective morphism

.. MATH::

    \phi_E: X_0(N) \rightarrow E.

from the modular curve `X_0(N)`, where `N` is the conductor of `E`.
The map sends the cusp `\infty` to the origin of `E`.

EXAMPLES::

    sage: phi = EllipticCurve('11a1').modular_parametrization()
    sage: phi
    Modular parameterization
     from the upper half plane
       to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
    sage: phi(0.5+CDF(I))
    (285684.320516... + 7.0...e-11*I : 1.526964169...e8 + 5.6...e-8*I : 1.00000000000000)
    sage: phi.power_series(prec = 7)
    (q^-2 + 2*q^-1 + 4 + 5*q + 8*q^2 + q^3 + 7*q^4 + O(q^5),
     -q^-3 - 3*q^-2 - 7*q^-1 - 13 - 17*q - 26*q^2 - 19*q^3 + O(q^4))

AUTHORS:

- Chris Wuthrich (02/10): moved from ell_rational_field.py.
"""
######################################################################
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
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
######################################################################

from . import heegner

from sage.rings.laurent_series_ring import LaurentSeriesRing
from sage.rings.rational_field import RationalField
from sage.rings.complex_mpfr import ComplexField
from sage.rings.rational_field import QQ


class ModularParameterization:
    r"""
    This class represents the modular parametrization of an elliptic curve

    .. MATH::

        \phi_E: X_0(N) \rightarrow E.

    Evaluation is done by passing through the lattice representation of `E`.

    EXAMPLES::

        sage: phi = EllipticCurve('11a1').modular_parametrization()
        sage: phi
        Modular parameterization
         from the upper half plane
           to Elliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20
              over Rational Field
    """
    def __init__(self, E):
        r"""
        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.ell_rational_field import ModularParameterization
            sage: phi = ModularParameterization(EllipticCurve('389a'))
            sage: phi(CC.0/5)
            (27.1965586309057 : -144.727322178982 : 1.00000000000000)

            sage: phi == loads(dumps(phi))
            True
        """
        self._E = E

    def curve(self):
        """
        Return the curve associated to this modular parametrization.

        EXAMPLES::

            sage: E = EllipticCurve('15a')
            sage: phi = E.modular_parametrization()
            sage: phi.curve() is E
            True
        """
        return self._E

    def __repr__(self):
        """
        TESTS::

            sage: E = EllipticCurve('37a')
            sage: phi = E.modular_parametrization()
            sage: phi
            Modular parameterization from the upper half plane to Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field
            sage: phi.__repr__()
            'Modular parameterization from the upper half plane to Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field'
        """
        return "Modular parameterization from the upper half plane to %s" % self._E

    def __eq__(self, other):
        r"""
        Compare two modular parametrizations by simply comparing the elliptic curves.

        EXAMPLES::

            sage: E = EllipticCurve('37a1')
            sage: phi = E.modular_parametrization()
            sage: phi == phi
            True
        """
        if not isinstance(other, ModularParameterization):
            return False

        return self._E == other._E

    def __ne__(self, other):
        """
        Check whether ``self`` is not equal to ``other``.

        EXAMPLES::

            sage: E = EllipticCurve('37a1')
            sage: phi = E.modular_parametrization()
            sage: phi != phi
            False
        """
        return not (self == other)

    def __call__(self, z, prec=None):
        r"""
        Evaluate ``self`` at a point `z \in X_0(N)` where `z` is given by a
        representative in the upper half plane.

        All computations are done with ``prec``
        bits of precision. If ``prec`` is not given, use the precision of `z`.

        EXAMPLES::

            sage: E = EllipticCurve('37a')
            sage: phi = E.modular_parametrization()
            sage: phi((sqrt(7)*I - 17)/74, 53)                                          # needs sage.symbolic
            (...e-16 - ...e-16*I : ...e-16 + ...e-16*I : 1.00000000000000)

        Verify that the mapping is invariant under the action of `\Gamma_0(N)`
        on the upper half plane::

            sage: E = EllipticCurve('11a')
            sage: phi = E.modular_parametrization()
            sage: tau = CC((1+1j)/5)
            sage: phi(tau)
            (-3.92181329652811 - 12.2578555525366*I : 44.9649874434872 + 14.3257120944681*I : 1.00000000000000)
            sage: phi(tau+1)
            (-3.92181329652810 - 12.2578555525366*I : 44.9649874434872 + 14.3257120944681*I : 1.00000000000000)
            sage: phi((6*tau+1) / (11*tau+2))
            (-3.9218132965285... - 12.2578555525369*I : 44.964987443489... + 14.325712094467...*I : 1.00000000000000)

        We can also apply the modular parametrization to a Heegner point on `X_0(N)`::

            sage: H = heegner_points(389,-7,5); H
            All Heegner points of conductor 5 on X_0(389) associated to QQ[sqrt(-7)]
            sage: x = H[0]; x
            Heegner point 5/778*sqrt(-7) - 147/778 of discriminant -7 and conductor 5 on X_0(389)
            sage: E = EllipticCurve('389a'); phi = E.modular_parametrization()
            sage: phi(x)
            Heegner point of discriminant -7 and conductor 5 on elliptic curve of conductor 389
            sage: phi(x).quadratic_form()
            389*x^2 + 147*x*y + 14*y^2

        ALGORITHM:

        Integrate the modular form attached to this elliptic curve from
        `z` to `\infty` to get a point on the lattice representation of
        `E`, then use the Weierstrass `\wp` function to map it to the
        curve itself.
        """
        from sage.misc.verbose import verbose
        if isinstance(z, heegner.HeegnerPointOnX0N):
            return z.map_to_curve(self.curve())
        # Map to the CC of CC/PeriodLattice.
        tm = verbose("Evaluating modular parameterization to precision %s bits" % prec)
        w = self.map_to_complex_numbers(z, prec=prec)
        # Map to E via Weierstrass P
        z = self._E.elliptic_exponential(w)
        verbose("Finished evaluating modular parameterization", tm)
        return z

    def map_to_complex_numbers(self, z, prec=None):
        r"""
        Evaluate ``self`` at a point `z \in X_0(N)` where `z` is given by
        a representative in the upper half plane, returning a point in
        the complex numbers.

        All computations are done with ``prec`` bits
        of precision.  If ``prec`` is not given, use the precision of `z`.
        Use self(z) to compute the image of z on the Weierstrass equation
        of the curve.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: E = EllipticCurve('37a'); phi = E.modular_parametrization()
            sage: x = polygen(ZZ, 'x')
            sage: tau = (sqrt(7)*I - 17)/74
            sage: z = phi.map_to_complex_numbers(tau); z
            0.929592715285395 - 1.22569469099340*I
            sage: E.elliptic_exponential(z)
            (...e-16 - ...e-16*I : ...e-16 + ...e-16*I : 1.00000000000000)
            sage: phi(tau)
            (...e-16 - ...e-16*I : ...e-16 + ...e-16*I : 1.00000000000000)
        """
        if prec is None:
            try:
                prec = z.parent().prec()
            except AttributeError:
                prec = 53
        CC = ComplexField(prec)
        if z in QQ:
            raise NotImplementedError
        z = CC(z)
        if z.imag() <= 0:
            raise ValueError("Point must be in the upper half plane")
        # TODO: for very small imaginary part, maybe try to transform under
        # \Gamma_0(N) to a better representative?
        q = (2 * CC.gen() * CC.pi() * z).exp()
        #  nterms'th term is less than 2**-(prec+10) (c.f. eclib code)
        nterms = (-(prec + 10) / q.abs().log2()).ceil()
        # Use Horner's rule to sum the integral of the form
        enumerated_an = list(enumerate(self._E.anlist(nterms)))[1:]
        lattice_point = 0
        for n, an in reversed(enumerated_an):
            lattice_point += an / n
            lattice_point *= q
        return lattice_point

    def power_series(self, prec=20):
        r"""
        Return the power series of this modular parametrization.

        The curve must be a minimal model.  The prec parameter determines
        the number of significant terms.  This means that X will be given up
        to O(q^(prec-2)) and Y will be given up to O(q^(prec-3)).

        OUTPUT:

        A list of two Laurent series [`X(x)`,`Y(x)`] of degrees `-2`, `-3`,
        respectively, which satisfy the equation of the elliptic curve.
        There are modular functions on `\Gamma_0(N)` where `N` is the
        conductor.

        The series should satisfy the differential equation

        .. MATH::

            \frac{\mathrm{d}X}{2Y + a_1 X + a_3} = \frac{f(q)\, \mathrm{d}q}{q}

        where `f` is ``self.curve().q_expansion()``.

        EXAMPLES::

            sage: E = EllipticCurve('389a1')
            sage: phi = E.modular_parametrization()
            sage: X, Y = phi.power_series(prec=10)
            sage: X
            q^-2 + 2*q^-1 + 4 + 7*q + 13*q^2 + 18*q^3 + 31*q^4 + 49*q^5 + 74*q^6 + 111*q^7 + O(q^8)
            sage: Y
            -q^-3 - 3*q^-2 - 8*q^-1 - 17 - 33*q - 61*q^2 - 110*q^3 - 186*q^4 - 320*q^5 - 528*q^6 + O(q^7)
            sage: X,Y = phi.power_series()
            sage: X
            q^-2 + 2*q^-1 + 4 + 7*q + 13*q^2 + 18*q^3 + 31*q^4 + 49*q^5 + 74*q^6 + 111*q^7 + 173*q^8 + 251*q^9 + 379*q^10 + 560*q^11 + 824*q^12 + 1199*q^13 + 1773*q^14 + 2548*q^15 + 3722*q^16 + 5374*q^17 + O(q^18)
            sage: Y
            -q^-3 - 3*q^-2 - 8*q^-1 - 17 - 33*q - 61*q^2 - 110*q^3 - 186*q^4 - 320*q^5 - 528*q^6 - 861*q^7 - 1383*q^8 - 2218*q^9 - 3472*q^10 - 5451*q^11 - 8447*q^12 - 13020*q^13 - 19923*q^14 - 30403*q^15 - 46003*q^16 + O(q^17)

        The following should give 0, but only approximately::

            sage: q = X.parent().gen()
            sage: E.defining_polynomial()(X,Y,1) + O(q^11) == 0
            True

        Note that below we have to change variable from `x` to `q`::

            sage: a1,_,a3,_,_ = E.a_invariants()
            sage: f = E.q_expansion(17)
            sage: q = f.parent().gen()
            sage: f/q == (X.derivative()/(2*Y+a1*X+a3))
            True
        """
        R = LaurentSeriesRing(RationalField(), 'q')
        if not self._E.is_minimal():
            raise NotImplementedError("only implemented for minimal curves")
        XY = self._E.pari_mincurve().elltaniyama(prec - 1)
        return R(XY[0]), R(XY[1])
