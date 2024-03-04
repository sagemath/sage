r"""
Scalar-multiplication morphisms of elliptic curves

This class provides an :class:`EllipticCurveHom` instantiation for
multiplication-by-`m` maps on elliptic curves.

EXAMPLES:

We can construct and evaluate scalar multiplications::

    sage: from sage.schemes.elliptic_curves.hom_scalar import EllipticCurveHom_scalar
    sage: E = EllipticCurve('77a1')
    sage: phi = E.scalar_multiplication(5); phi
    Scalar-multiplication endomorphism [5] of Elliptic Curve defined by y^2 + y = x^3 + 2*x over Rational Field
    sage: P = E(2,3)
    sage: phi(P)
    (30 : 164 : 1)

The usual :class:`EllipticCurveHom` methods are supported::

    sage: phi.degree()
    25
    sage: phi.kernel_polynomial()
    x^12 + 124/5*x^10 + 19*x^9 - 84*x^8 + 24*x^7 - 483*x^6 - 696/5*x^5 - 448*x^4 - 37*x^3 - 332*x^2 - 84*x + 47/5
    sage: phi.rational_maps()
    ((x^25 - 200*x^23 - 520*x^22 + 9000*x^21 + ... + 1377010*x^3 + 20360*x^2 - 39480*x + 2209),
     (10*x^36*y - 620*x^36 + 3240*x^34*y - 44880*x^34 + ... + 424927560*x*y + 226380480*x + 42986410*y + 20974090)/(1250*x^36 + 93000*x^34 + 71250*x^33 + 1991400*x^32 + ... + 1212964050*x^3 + 138715800*x^2 - 27833400*x + 1038230))
    sage: phi.dual()
    Scalar-multiplication endomorphism [5] of Elliptic Curve defined by y^2 + y = x^3 + 2*x over Rational Field
    sage: phi.dual() is phi
    True
    sage: phi.formal()
    5*t - 310*t^4 - 2496*t^5 + 10540*t^7 + ... - 38140146674516*t^20 - 46800256902400*t^21 + 522178541079910*t^22 + O(t^23)
    sage: phi.is_normalized()
    False
    sage: phi.is_separable()
    True
    sage: phi.is_injective()
    False
    sage: phi.is_surjective()
    True

Contrary to constructing an :class:`EllipticCurveIsogeny` from
the division polynomial, :class:`EllipticCurveHom_scalar` can
deal with huge scalars very quickly::

    sage: E = EllipticCurve(GF(2^127-1), [1,2,3,4,5])
    sage: phi = E.scalar_multiplication(9^99); phi
    Scalar-multiplication endomorphism [29512665430652752148753480226197736314359272517043832886063884637676943433478020332709411004889] of Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over Finite Field of size 170141183460469231731687303715884105727
    sage: phi(E(1,2))
    (82124533143060719620799539030695848450 : 17016022038624814655722682134021402379 : 1)

Composition of scalar multiplications results in another scalar
multiplication::

    sage: E = EllipticCurve(GF(19), [4,4])
    sage: phi = E.scalar_multiplication(-3); phi
    Scalar-multiplication endomorphism [-3] of Elliptic Curve defined by y^2 = x^3 + 4*x + 4 over Finite Field of size 19
    sage: psi = E.scalar_multiplication(7); psi
    Scalar-multiplication endomorphism [7] of Elliptic Curve defined by y^2 = x^3 + 4*x + 4 over Finite Field of size 19
    sage: phi * psi
    Scalar-multiplication endomorphism [-21] of Elliptic Curve defined by y^2 = x^3 + 4*x + 4 over Finite Field of size 19
    sage: psi * phi
    Scalar-multiplication endomorphism [-21] of Elliptic Curve defined by y^2 = x^3 + 4*x + 4 over Finite Field of size 19
    sage: phi * psi == psi * phi
    True
    sage: -phi == E.scalar_multiplication(-1) * phi
    True

The zero endomorphism `[0]` is supported::

    sage: E = EllipticCurve(GF(71), [1,1])
    sage: zero = E.scalar_multiplication(0); zero
    Scalar-multiplication endomorphism [0] of Elliptic Curve defined by y^2 = x^3 + x + 1 over Finite Field of size 71
    sage: zero.is_zero()
    True
    sage: zero.is_injective()
    False
    sage: zero.is_surjective()
    False
    sage: zero(E.random_point())
    (0 : 1 : 0)

Retrieving multiplication-by-`m` maps when `m` is divisible by the
characteristic also works (since :issue:`37096`)::

    sage: E = EllipticCurve(GF(7), [1,0])
    sage: phi = E.scalar_multiplication(7); phi
    Scalar-multiplication endomorphism [7] of Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 7
    sage: phi.rational_maps()
    (x^49, -y^49)
    sage: phi.x_rational_map()
    x^49
    sage: psi = E.scalar_multiplication(-2); psi
    Scalar-multiplication endomorphism [-2] of Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 7
    sage: chi = E.scalar_multiplication(-14); chi
    Scalar-multiplication endomorphism [-14] of Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 7
    sage: chi == psi * phi
    True
    sage: chi.rational_maps()
    ((x^196 - 2*x^98 + 1)/(-3*x^147 - 3*x^49),
     (-x^294*y^49 + 2*x^196*y^49 - 2*x^98*y^49 + y^49)/(-x^294 - 2*x^196 - x^98))
    sage: chi.x_rational_map()
    (2*x^196 + 3*x^98 + 2)/(x^147 + x^49)
    sage: chi.rational_maps() == tuple(f(*phi.rational_maps()) for f in psi.rational_maps())
    True
    sage: chi.x_rational_map() == psi.x_rational_map()(phi.x_rational_map())
    True

::

    sage: E = EllipticCurve(GF(7), [0,1])
    sage: phi = E.scalar_multiplication(7); phi
    Scalar-multiplication endomorphism [7] of Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 7
    sage: phi.rational_maps()   # known bug -- #6413
    ((-3*x^49 - x^28 - x^7)/(x^42 - x^21 + 2),
     (-x^72*y - 3*x^69*y - 3*x^66*y - x^63*y + 3*x^51*y + 2*x^48*y + 2*x^45*y + 3*x^42*y - x^9*y - 3*x^6*y - 3*x^3*y - y)/(x^63 + 2*x^42 - x^21 - 1))
    sage: phi.x_rational_map()
    (4*x^49 + 6*x^28 + 6*x^7)/(x^42 + 6*x^21 + 2)

TESTS::

    sage: E = EllipticCurve(j = GF(65537^3).random_element())
    sage: m = randrange(-2^99, +2^99)
    sage: phi = E.scalar_multiplication(m)
    sage: phi.degree() == m**2
    True
    sage: P = E.random_point()
    sage: phi(P) == m*P
    True

AUTHORS:

- Lorenz Panny (2021): implement :class:`EllipticCurveHom_scalar`
"""

from sage.misc.cachefunc import cached_method
from sage.structure.richcmp import richcmp

from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

from sage.schemes.elliptic_curves.ell_generic import EllipticCurve_generic
from sage.schemes.elliptic_curves.weierstrass_morphism import negation_morphism
from sage.schemes.elliptic_curves.hom import EllipticCurveHom


class EllipticCurveHom_scalar(EllipticCurveHom):

    def __init__(self, E, m):
        """
        Construct a scalar-multiplication map on an elliptic curve.

        TESTS::

            sage: from sage.schemes.elliptic_curves.hom_scalar import EllipticCurveHom_scalar
            sage: E = EllipticCurve([1,1])
            sage: EllipticCurveHom_scalar(E, 123)
            Scalar-multiplication endomorphism [123] of Elliptic Curve defined by y^2 = x^3 + x + 1 over Rational Field
        """
        if not isinstance(E, EllipticCurve_generic):
            raise ValueError(f'not an elliptic curve: {E}')

        self._m = ZZ(m)

        self._degree = self._m**2
        self._domain = self._codomain = E

        EllipticCurveHom.__init__(self, self._domain, self._codomain)

        # TODO: should probably be in EllipticCurveHom?
        self._base_ring = self._domain.base_ring()
        self._poly_ring = PolynomialRing(self._base_ring, ['x'])
        self._mpoly_ring = PolynomialRing(self._base_ring, ['x','y'])

        self._rational_maps = None

    def _call_(self, P):
        """
        Evaluate this scalar-multiplication map `[m]` at a point `P`,
        i.e., return `[m]P`.

        TESTS::

            sage: p = random_prime(2^22)
            sage: q = p^randrange(1,5)
            sage: E = EllipticCurve_from_j(GF(q).random_element())
            sage: m = randrange(-9^99, 9^99)
            sage: phi = E.scalar_multiplication(m)
            sage: P = E.random_point()
            sage: phi(P) == m*P
            True
        """
        if P not in self._domain:
            raise ValueError(f'{P} is not a point on {self._domain}')
        return self._m * P

    def _eval(self, P):
        """
        Less strict evaluation method for internal use.

        In particular, this can be used to evaluate ``self`` at a
        point defined over an extension field.

        INPUT: a sequence of 3 coordinates defining a point on ``self``

        OUTPUT: the result of evaluating ``self`` at the given point

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_scalar import EllipticCurveHom_scalar
            sage: E = EllipticCurve(j=GF(419)(1728))
            sage: psi = EllipticCurveHom_scalar(E, 13)
            sage: P = E.change_ring(GF(419**2)).lift_x(5)
            sage: P = min({P, -P})  # fix choice of y
            sage: Q = psi._eval(P); Q
            (134 : 210*z2 + 314 : 1)
            sage: Q.curve()
            Elliptic Curve defined by y^2 = x^3 + x over Finite Field in z2 of size 419^2
        """
        if self._domain.defining_polynomial()(*P):
            raise ValueError(f'{P} not on {self._domain}')
        return self._m * P

    def _repr_(self):
        """
        Return basic facts about this scalar multiplication as a string.

        TESTS::

            sage: E = EllipticCurve([i,i])
            sage: E.scalar_multiplication(777)
            Scalar-multiplication endomorphism [777] of Elliptic Curve defined by y^2 = x^3 + I*x + I over Number Field in I with defining polynomial x^2 + 1 with I = 1*I
        """
        return f'Scalar-multiplication endomorphism [{self._m}] of {self._domain}'

    # EllipticCurveHom methods

    @staticmethod
    def _composition_impl(self, other):
        """
        Helper method to compose other elliptic-curve morphisms with
        :class:`EllipticCurveHom_scalar` objects. Called by
        :meth:`EllipticCurveHom._composition_`.

        This method only handles composing two scalar multiplications;
        all other cases are dealt with elsewhere.

        TESTS::

            sage: E = EllipticCurve([1,2,3,4,5])
            sage: phi = E.scalar_multiplication(5)
            sage: psi = E.scalar_multiplication(-7)
            sage: phi * psi     # implicit doctest
            Scalar-multiplication endomorphism [-35] of Elliptic Curve defined by y^2 + x*y + 3*y = x^3 + 2*x^2 + 4*x + 5 over Rational Field

        ::

            sage: phi._composition_impl(phi, E.automorphisms()[0])
            NotImplemented
        """
        if isinstance(self, EllipticCurveHom_scalar) and isinstance(other, EllipticCurveHom_scalar):
            assert self._domain == other._domain
            return EllipticCurveHom_scalar(self._domain, self._m * other._m)
        return NotImplemented

    def degree(self):
        """
        Return the degree of this scalar-multiplication morphism.

        The map `[m]` has degree `m^2`.

        EXAMPLES::

            sage: E = EllipticCurve(GF(23), [0,1])
            sage: phi = E.scalar_multiplication(1111111)
            sage: phi.degree()
            1234567654321

        TESTS:

        The degree is still `m^2` even in the inseparable case::

            sage: E = EllipticCurve(GF(23), [1,1])
            sage: E.scalar_multiplication(23).degree()
            529
            sage: E = EllipticCurve(GF(23), [0,1])
            sage: E.scalar_multiplication(23).degree()
            529
        """
        return self._degree

    def _richcmp_(self, other, op):
        """
        Compare this scalar multiplication to another elliptic-curve morphism.

        .. WARNING::

            This method sometimes calls :meth:`EllipticCurveHom._richcmp_`,
            which sometimes compares :meth:`rational_maps`. Therefore, the
            complexity is at least quadratic in `m` in the worst case.

        EXAMPLES::

            sage: E = EllipticCurve([i,i])
            sage: phi = E.scalar_multiplication(-5)
            sage: psi = E.scalar_multiplication(5)
            sage: phi == -psi
            True

        TESTS::

            sage: from sage.schemes.elliptic_curves.weierstrass_morphism import negation_morphism
            sage: neg = negation_morphism(E)
            sage: phi == neg*psi == psi*neg == -psi
            True
        """
        if isinstance(other, EllipticCurveHom_scalar):
            return richcmp((self._domain, self._m), (other._domain, other._m), op)
        return EllipticCurveHom._richcmp_(self, other, op)

    def rational_maps(self):
        """
        Return the pair of explicit rational maps defining this scalar
        multiplication.

        ALGORITHM: :meth:`EllipticCurve_generic.multiplication_by_m`

        EXAMPLES::

            sage: E = EllipticCurve('77a1')
            sage: phi = E.scalar_multiplication(5)
            sage: phi.rational_maps()
            ((x^25 - 200*x^23 - 520*x^22 + ... + 368660*x^2 + 163195*x + 16456)/(25*x^24 + 1240*x^22 + 950*x^21 + ... + 20360*x^2 - 39480*x + 2209),
             (10*x^36*y - 620*x^36 + 3240*x^34*y - ... + 226380480*x + 42986410*y + 20974090)/(1250*x^36 + 93000*x^34 + 71250*x^33 + ... + 138715800*x^2 - 27833400*x + 1038230))
            sage: P = (2,3)
            sage: Q = tuple(r(P) for r in phi.rational_maps()); Q
            (30, 164)
            sage: E(Q) == 5*E(P)
            True

        TESTS::

            sage: {r.parent() for r in phi.rational_maps()}
            {Fraction Field of Multivariate Polynomial Ring in x, y over Rational Field}
        """
        if not self._rational_maps or None in self._rational_maps:
            if not self._m:
                raise ValueError('[0] is not expressible in (x,y) coordinates')
            self._rational_maps = self._domain.multiplication_by_m(self._m)
        return self._rational_maps

    def x_rational_map(self):
        """
        Return the `x`-coordinate rational map of this scalar
        multiplication.

        ALGORITHM: :meth:`EllipticCurve_generic.multiplication_by_m`

        EXAMPLES::

            sage: E = EllipticCurve(GF(65537), [1,2,3,4,5])
            sage: phi = E.scalar_multiplication(7)
            sage: phi.x_rational_map() == phi.rational_maps()[0]
            True

        TESTS::

            sage: phi.x_rational_map().parent()
            Fraction Field of Univariate Polynomial Ring in x over Finite Field of size 65537
        """
        if not self._rational_maps:
            if not self._m:
                raise ValueError('[0] is not expressible in (x,y) coordinates')
            h = self._domain.multiplication_by_m(self._m, x_only=True)
            self._rational_maps = (self._mpoly_ring.fraction_field()(h), None)
        f,g = map(self._poly_ring, (self._rational_maps[0].numerator(),
                                    self._rational_maps[0].denominator()))
        return f / g

    def scaling_factor(self):
        r"""
        Return the Weierstrass scaling factor associated to this
        scalar multiplication.

        The scaling factor is the constant `u` (in the base field)
        such that `\varphi^* \omega_2 = u \omega_1`, where
        `\varphi: E_1\to E_2` is this morphism and `\omega_i` are
        the standard Weierstrass differentials on `E_i` defined by
        `\mathrm dx/(2y+a_1x+a_3)`.

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: phi = E.scalar_multiplication(5)
            sage: u = phi.scaling_factor()
            sage: u == phi.formal()[1]
            True
            sage: u == E.multiplication_by_m_isogeny(5).scaling_factor()
            doctest:warning ... DeprecationWarning: ...
            True

        The scaling factor lives in the base ring::

            sage: E = EllipticCurve(GF(101^2), [5,5])
            sage: phi = E.scalar_multiplication(123)
            sage: phi.scaling_factor()
            22
            sage: phi.scaling_factor().parent()
            Finite Field in z2 of size 101^2

        ALGORITHM: The scaling factor equals the scalar that is being
        multiplied by.
        """
        return self._base_ring(self._m)

    @cached_method
    def kernel_polynomial(self):
        r"""
        Return the kernel polynomial of this scalar-multiplication map.
        (When `m=0`, return `0`.)

        EXAMPLES::

            sage: E = EllipticCurve(GF(997), [7,7,7,7,7])
            sage: phi = E.scalar_multiplication(5)
            sage: phi.kernel_polynomial()
            x^12 + 77*x^11 + 380*x^10 + 198*x^9 + 840*x^8 + 376*x^7 + 946*x^6 + 848*x^5 + 246*x^4 + 778*x^3 + 77*x^2 + 518*x + 28

        ::

            sage: E = EllipticCurve(GF(997), [5,6,7,8,9])
            sage: phi = E.scalar_multiplication(11)
            sage: phi.kernel_polynomial()
            x^60 + 245*x^59 + 353*x^58 + 693*x^57 + 499*x^56 + 462*x^55 + 820*x^54 + 962*x^53 + ... + 736*x^7 + 939*x^6 + 429*x^5 + 267*x^4 + 116*x^3 + 770*x^2 + 491*x + 519

        TESTS::

            sage: E = EllipticCurve(j = GF(997^6).random_element())
            sage: m = choice([+1,-1]) * randrange(1,8)
            sage: phi = E.scalar_multiplication(m)
            sage: phi.kernel_polynomial() == phi.x_rational_map().denominator().monic().radical()
            True

        ::

            sage: E.scalar_multiplication(randint(-10,+10)).kernel_polynomial().parent()
            Univariate Polynomial Ring in x over Finite Field in z6 of size 997^6
        """
        if not self._m:
            return self._poly_ring(0)
        # TODO: inseparable case should be consistent with Frobenius' .kernel_polynomial()
        return self._domain.division_polynomial(self._m.abs()).monic().radical()

    def dual(self):
        """
        Return the dual isogeny of this scalar-multiplication map.

        This method simply returns ``self`` as scalars are self-dual.

        EXAMPLES::

            sage: E = EllipticCurve([5,5])
            sage: phi = E.scalar_multiplication(5)
            sage: phi.dual() is phi
            True
        """
        return self

    def inseparable_degree(self):
        r"""
        Return the inseparable degree of this scalar-multiplication map.

        EXAMPLES::

            sage: E = EllipticCurve(GF(7), [0,1])
            sage: E.is_supersingular()
            False
            sage: E.scalar_multiplication(4).inseparable_degree()
            1
            sage: E.scalar_multiplication(-7).inseparable_degree()
            7

        ::

            sage: E = EllipticCurve(GF(7), [1,0])
            sage: E.is_supersingular()
            True
            sage: E.scalar_multiplication(4).inseparable_degree()
            1
            sage: E.scalar_multiplication(-7).inseparable_degree()
            49
        """
        p = self.base_ring().characteristic()
        if not p:
            return ZZ.one()
        v = self._m.valuation(p)
        if not v:
            return ZZ.one()
        rk = 1 + self._domain.is_supersingular()
        return p**(rk*v)

    def __neg__(self):
        """
        Negate this scalar-multiplication map, i.e., return `[-m]`
        when this morphism equals `[m]`.

        If rational maps have been computed already, they will be
        reused for the negated morphism.

        EXAMPLES::

            sage: E = EllipticCurve(GF(2^8), [1,0,1,0,1])
            sage: phi = E.scalar_multiplication(23)
            sage: -phi
            Scalar-multiplication endomorphism [-23] of Elliptic Curve defined by y^2 + x*y + y = x^3 + 1 over Finite Field in z8 of size 2^8

        TESTS::

            sage: E = EllipticCurve(GF(79), [7,7])
            sage: phi = E.scalar_multiplication(5)
            sage: _ = phi.rational_maps()
            sage: (-phi)._rational_maps
            ((x^25 + 11*x^23 - 24*x^22 - ... - 7*x^2 + 34*x + 21)/(25*x^24 - 5*x^22 - 23*x^21 - ... - 11*x^2 + 36*x + 21),
             (29*x^36*y + 22*x^34*y - 27*x^33*y - ... + 14*x^2*y - 33*x*y + 37*y)/(9*x^36 + 21*x^34 - 14*x^33 + ... - 26*x^2 + 18*x + 7))
        """
        result = EllipticCurveHom_scalar(self._domain, -self._m)
        if self._rational_maps is not None:
            w = negation_morphism(self._domain).rational_maps()
            result._rational_maps = tuple(f(*w) if f is not None else None for f in self._rational_maps)
        return result
