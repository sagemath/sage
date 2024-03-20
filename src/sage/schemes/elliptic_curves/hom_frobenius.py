# sage_setup: distribution = sagemath-schemes
# sage.doctest: needs sage.rings.finite_rings
r"""
Frobenius isogenies of elliptic curves

Frobenius isogenies only exist in positive characteristic `p`. They
are given by `\pi_n:(x,y)\mapsto (x^{p^n},y^{p^n})`.

This class implements `\pi_n` for `n \geq 0`. Together with existing
tools for composing isogenies (see :class:`EllipticCurveHom_composite`),
we can therefore represent arbitrary inseparable isogenies in Sage.

EXAMPLES:

Constructing a Frobenius isogeny is straightforward::

    sage: from sage.schemes.elliptic_curves.hom_frobenius import EllipticCurveHom_frobenius
    sage: z5, = GF(17^5).gens()
    sage: E = EllipticCurve([z5,1])
    sage: pi = EllipticCurveHom_frobenius(E); pi
    Frobenius isogeny of degree 17:
      From: Elliptic Curve defined by y^2 = x^3 + z5*x + 1
            over Finite Field in z5 of size 17^5
      To:   Elliptic Curve defined by y^2 = x^3 + (9*z5^4+7*z5^3+10*z5^2+z5+14)*x + 1
            over Finite Field in z5 of size 17^5

By passing `n`, we can also construct higher-power Frobenius maps,
such as the Frobenius *endo*\morphism::

    sage: z5, = GF(7^5).gens()
    sage: E = EllipticCurve([z5,1])
    sage: pi = EllipticCurveHom_frobenius(E, 5); pi
    Frobenius endomorphism of degree 16807 = 7^5:
      From: Elliptic Curve defined by y^2 = x^3 + z5*x + 1
            over Finite Field in z5 of size 7^5
      To:   Elliptic Curve defined by y^2 = x^3 + z5*x + 1
            over Finite Field in z5 of size 7^5

The usual :class:`EllipticCurveHom` methods are supported::

    sage: z5, = GF(7^5).gens()
    sage: E = EllipticCurve([z5,1])
    sage: pi = EllipticCurveHom_frobenius(E,5)
    sage: pi.degree()
    16807
    sage: pi.rational_maps()
    (x^16807, y^16807)
    sage: pi.formal()                   # known bug
    ...
    sage: pi.is_normalized()            # known bug
    ...
    sage: pi.is_separable()
    False
    sage: pi.is_injective()
    True
    sage: pi.is_surjective()
    True

Computing the dual of Frobenius is supported as well::

    sage: E = EllipticCurve([GF(17^6).gen(), 0])
    sage: pi = EllipticCurveHom_frobenius(E)
    sage: pihat = pi.dual(); pihat
    Isogeny of degree 17
     from Elliptic Curve defined by y^2 = x^3 + (15*z6^5+5*z6^4+8*z6^3+12*z6^2+11*z6+7)*x
          over Finite Field in z6 of size 17^6
       to Elliptic Curve defined by y^2 = x^3 + z6*x
          over Finite Field in z6 of size 17^6
    sage: pihat.is_separable()
    True
    sage: pihat * pi == EllipticCurveHom_scalar(E,17)   # known bug -- #6413
    True

A supersingular example (with purely inseparable dual)::

    sage: E = EllipticCurve([0, GF(17^6).gen()])
    sage: E.is_supersingular()
    True
    sage: pi1 = EllipticCurveHom_frobenius(E)
    sage: pi1hat = pi1.dual(); pi1hat
    Composite morphism of degree 17 = 17*1:
      From: Elliptic Curve defined by y^2 = x^3 + (15*z6^5+5*z6^4+8*z6^3+12*z6^2+11*z6+7)
            over Finite Field in z6 of size 17^6
      To:   Elliptic Curve defined by y^2 = x^3 + z6
            over Finite Field in z6 of size 17^6
    sage: pi6 = EllipticCurveHom_frobenius(E,6)
    sage: pi6hat = pi6.dual(); pi6hat
    Composite morphism of degree 24137569 = 24137569*1:
      From: Elliptic Curve defined by y^2 = x^3 + z6
            over Finite Field in z6 of size 17^6
      To:   Elliptic Curve defined by y^2 = x^3 + z6
            over Finite Field in z6 of size 17^6
    sage: pi6hat.factors()
    (Frobenius endomorphism of degree 24137569 = 17^6:
       From: Elliptic Curve defined by y^2 = x^3 + z6
             over Finite Field in z6 of size 17^6
       To:   Elliptic Curve defined by y^2 = x^3 + z6
             over Finite Field in z6 of size 17^6,
     Elliptic-curve endomorphism of
      Elliptic Curve defined by y^2 = x^3 + z6 over Finite Field in z6 of size 17^6
       Via:  (u,r,s,t) = (2*z6^5 + 10*z6^3 + z6^2 + 8, 0, 0, 0))


TESTS::

    sage: z5, = GF(17^5).gens()
    sage: E = EllipticCurve([z5,1])
    sage: fs = [EllipticCurveHom_frobenius(E)]
    sage: while fs[-1].codomain() != E:
    ....:     fs.append(EllipticCurveHom_frobenius(fs[-1].codomain()))
    sage: fs
    [Frobenius isogeny of degree 17:
       From: Elliptic Curve defined by y^2 = x^3 + z5*x + 1 over Finite Field in z5 of size 17^5
       To:   Elliptic Curve defined by y^2 = x^3 + (9*z5^4+7*z5^3+10*z5^2+z5+14)*x + 1 over Finite Field in z5 of size 17^5,
     Frobenius isogeny of degree 17:
       From: Elliptic Curve defined by y^2 = x^3 + (9*z5^4+7*z5^3+10*z5^2+z5+14)*x + 1 over Finite Field in z5 of size 17^5
       To:   Elliptic Curve defined by y^2 = x^3 + (14*z5^4+7*z5^3+16*z5^2+14*z5+1)*x + 1 over Finite Field in z5 of size 17^5,
     Frobenius isogeny of degree 17:
       From: Elliptic Curve defined by y^2 = x^3 + (14*z5^4+7*z5^3+16*z5^2+14*z5+1)*x + 1 over Finite Field in z5 of size 17^5
       To:   Elliptic Curve defined by y^2 = x^3 + (16*z5^4+6*z5^3+7*z5^2+14*z5+6)*x + 1 over Finite Field in z5 of size 17^5,
     Frobenius isogeny of degree 17:
       From: Elliptic Curve defined by y^2 = x^3 + (16*z5^4+6*z5^3+7*z5^2+14*z5+6)*x + 1 over Finite Field in z5 of size 17^5
       To:   Elliptic Curve defined by y^2 = x^3 + (12*z5^4+14*z5^3+z5^2+4*z5+13)*x + 1 over Finite Field in z5 of size 17^5,
     Frobenius isogeny of degree 17:
       From: Elliptic Curve defined by y^2 = x^3 + (12*z5^4+14*z5^3+z5^2+4*z5+13)*x + 1 over Finite Field in z5 of size 17^5
       To:   Elliptic Curve defined by y^2 = x^3 + z5*x + 1 over Finite Field in z5 of size 17^5]
    sage: prod(fs[::-1])
    Composite morphism of degree 1419857 = 17^5:
      From: Elliptic Curve defined by y^2 = x^3 + z5*x + 1 over Finite Field in z5 of size 17^5
      To:   Elliptic Curve defined by y^2 = x^3 + z5*x + 1 over Finite Field in z5 of size 17^5

::

    sage: EllipticCurveHom_frobenius(EllipticCurve(GF(5),[1,1]), -1)
    Traceback (most recent call last):
    ...
    ValueError: negative powers of Frobenius are not isogenies

::

    sage: EllipticCurveHom_frobenius(EllipticCurve('11a1'))
    Traceback (most recent call last):
    ...
    ValueError: Frobenius isogenies do not exist in characteristic zero

AUTHORS:

- Lorenz Panny (2021): implement :class:`EllipticCurveHom_frobenius`
- Mickaël Montessinos (2021): computing the dual of a Frobenius isogeny
"""

from sage.misc.cachefunc import cached_method
from sage.structure.sequence import Sequence

from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

from sage.rings.finite_rings.finite_field_base import FiniteField

from sage.schemes.elliptic_curves.ell_generic import EllipticCurve_generic
from sage.schemes.elliptic_curves.constructor import EllipticCurve

from sage.schemes.elliptic_curves.hom import EllipticCurveHom, find_post_isomorphism
from sage.schemes.elliptic_curves.ell_curve_isogeny import EllipticCurveIsogeny
from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
from sage.schemes.elliptic_curves.hom_scalar import EllipticCurveHom_scalar


class EllipticCurveHom_frobenius(EllipticCurveHom):

    _degree = None

    def __init__(self, E, power=1):
        r"""
        Construct a Frobenius isogeny on a given curve with a given
        power of the base-ring characteristic.

        Writing `n` for the parameter ``power`` (default: `1`), the
        isogeny is defined by `(x,y) \to (x^{p^n}, y^{p^n})` where
        `p` is the characteristic of the base ring.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_frobenius import EllipticCurveHom_frobenius
            sage: E = EllipticCurve(j=GF(11^2).gen())
            sage: EllipticCurveHom_frobenius(E)
            Frobenius isogeny of degree 11:
              From: Elliptic Curve defined by y^2 = x^3 + (2*z2+6)*x + (8*z2+8) over Finite Field in z2 of size 11^2
              To:   Elliptic Curve defined by y^2 = x^3 + (9*z2+3)*x + (3*z2+7) over Finite Field in z2 of size 11^2
            sage: EllipticCurveHom_frobenius(E, 2)
            Frobenius endomorphism of degree 121 = 11^2:
              From: Elliptic Curve defined by y^2 = x^3 + (2*z2+6)*x + (8*z2+8) over Finite Field in z2 of size 11^2
              To:   Elliptic Curve defined by y^2 = x^3 + (2*z2+6)*x + (8*z2+8) over Finite Field in z2 of size 11^2

        TESTS::

            sage: EllipticCurveHom_frobenius(EllipticCurve('11a1'))
            Traceback (most recent call last):
            ...
            ValueError: Frobenius isogenies do not exist in characteristic zero

        ::

            sage: EllipticCurveHom_frobenius(E, -1)
            Traceback (most recent call last):
            ...
            ValueError: negative powers of Frobenius are not isogenies
        """
        if not isinstance(E, EllipticCurve_generic):
            raise ValueError(f'not an elliptic curve: {E}')

        self._base_ring = E.base_ring()

        self._p = self._base_ring.characteristic()
        if self._p == 0:
            raise ValueError('Frobenius isogenies do not exist in characteristic zero')

        self._n = ZZ(power)
        if self._n < 0:
            raise ValueError('negative powers of Frobenius are not isogenies')

        self._degree = self._p ** self._n

        self._domain = E
        as_ = [a**self._degree for a in self._domain.a_invariants()]
        self._codomain = EllipticCurve(as_)

        EllipticCurveHom.__init__(self, self._domain, self._codomain)

        # over finite fields, isogenous curves have the same number of points
        # (see #32786)
        if isinstance(self._base_ring, FiniteField):
            self._domain._fetch_cached_order(self._codomain)
            self._codomain._fetch_cached_order(self._domain)

        self._poly_ring = PolynomialRing(self._base_ring, ['x'], sparse=True)
        self._mpoly_ring = PolynomialRing(self._base_ring, ['x','y'], sparse=True)
        self._xfield = self._poly_ring.fraction_field()
        self._xyfield = self._mpoly_ring.fraction_field()

    def _call_(self, P):
        """
        Evaluate this Frobenius isogeny at a point `P`.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_frobenius import EllipticCurveHom_frobenius
            sage: z2 = GF(11^2).gen()
            sage: E = EllipticCurve(j=z2)
            sage: pi = EllipticCurveHom_frobenius(E)
            sage: P = E(7, 9*z2+4)
            sage: pi(P)     # implicit doctest
            (7 : 2*z2 + 7 : 1)
        """
        return self._codomain(*(c**self._degree for c in P))

    def _eval(self, P):
        """
        Less strict evaluation method for internal use.

        In particular, this can be used to evaluate ``self`` at a
        point defined over an extension field.

        INPUT: a sequence of 3 coordinates defining a point on ``self``

        OUTPUT: the result of evaluating ``self`` at the given point

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_frobenius import EllipticCurveHom_frobenius
            sage: E = EllipticCurve(GF(11), [1,1])
            sage: pi = EllipticCurveHom_frobenius(E)
            sage: Ebar = E.change_ring(GF(11^6))
            sage: z6 = GF(11^6).gen()
            sage: P = Ebar.lift_x(GF(11^3).gen())
            sage: p = Ebar(6*z6^5 + 8*z6^4 + 8*z6^3 + 6*z6^2 + 10*z6 + 5, 2*z6^5 + 2*z6^4 + 2*z6^3 + 4*z6 + 6, 1)
            sage: Q = pi._eval(P)
            sage: q = Ebar(z6^5 + 3*z6^4 + 3*z6^3 + 6*z6^2 + 9, z6^5 + 10*z6^4 + 10*z6^3 + 5*z6^2 + 4*z6 + 8, 1)
            sage: (P == p and Q == q) or (P == -p and Q == -q)
            True
        """
        if self._domain.defining_polynomial()(*P):
            raise ValueError(f'{P} not on {self._domain}')
        k = Sequence(P).universe()
        return self._codomain.base_extend(k)(*(c**self._degree for c in P))

    def _repr_(self):
        """
        Return basic facts about this Frobenius isogeny as a string.

        TESTS::

            sage: from sage.schemes.elliptic_curves.hom_frobenius import EllipticCurveHom_frobenius
            sage: z2 = GF(11^2).gen()
            sage: E = EllipticCurve(j=z2)
            sage: EllipticCurveHom_frobenius(E)
            Frobenius isogeny of degree 11:
              From: Elliptic Curve defined by y^2 = x^3 + (2*z2+6)*x + (8*z2+8) over Finite Field in z2 of size 11^2
              To:   Elliptic Curve defined by y^2 = x^3 + (9*z2+3)*x + (3*z2+7) over Finite Field in z2 of size 11^2
            sage: EllipticCurveHom_frobenius(E, E.base_field().degree())
            Frobenius endomorphism of degree 121 = 11^2:
              From: Elliptic Curve defined by y^2 = x^3 + (2*z2+6)*x + (8*z2+8) over Finite Field in z2 of size 11^2
              To:   Elliptic Curve defined by y^2 = x^3 + (2*z2+6)*x + (8*z2+8) over Finite Field in z2 of size 11^2
        """
        kind = 'endomorphism' if self._codomain == self._domain else 'isogeny'
        degs_str = '' if self._n == 1 else f' = {self._p}^{self._n}'
        return f'Frobenius {kind} of degree {self._degree}{degs_str}:' \
                f'\n  From: {self._domain}' \
                f'\n  To:   {self._codomain}'

    # EllipticCurveHom methods

    def rational_maps(self):
        """
        Return the explicit rational maps defining this Frobenius
        isogeny as (sparse) bivariate rational maps in `x` and `y`.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_frobenius import EllipticCurveHom_frobenius
            sage: E = EllipticCurve(GF(11), [1,1])
            sage: pi = EllipticCurveHom_frobenius(E, 4)
            sage: pi.rational_maps()
            (x^14641, y^14641)

        TESTS:

        See :issue:`34811`::

            sage: pi.rational_maps()[0].parent()
            Fraction Field of Multivariate Polynomial Ring in x, y over Finite Field of size 11
            sage: pi.rational_maps()[1].parent()
            Fraction Field of Multivariate Polynomial Ring in x, y over Finite Field of size 11
        """
        x,y = self._xyfield.gens()
        return (x**self._degree, y**self._degree)

    def x_rational_map(self):
        """
        Return the `x`-coordinate rational map of this Frobenius
        isogeny as a (sparse) univariate rational map in `x`.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_frobenius import EllipticCurveHom_frobenius
            sage: E = EllipticCurve(GF(11), [1,1])
            sage: pi = EllipticCurveHom_frobenius(E, 4)
            sage: pi.x_rational_map()
            x^14641

        TESTS:

        See :issue:`34811`::

            sage: pi.x_rational_map().parent()
            Fraction Field of Sparse Univariate Polynomial Ring in x over Finite Field of size 11
        """
        x, = self._xfield.gens()
        return x**self._degree

    def scaling_factor(self):
        r"""
        Return the Weierstrass scaling factor associated to this
        Frobenius morphism.

        The scaling factor is the constant `u` (in the base field)
        such that `\varphi^* \omega_2 = u \omega_1`, where
        `\varphi: E_1\to E_2` is this morphism and `\omega_i` are
        the standard Weierstrass differentials on `E_i` defined by
        `\mathrm dx/(2y+a_1x+a_3)`.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_frobenius import EllipticCurveHom_frobenius
            sage: E = EllipticCurve(GF(11), [1,1])
            sage: pi = EllipticCurveHom_frobenius(E)
            sage: pi.formal()
            t^11 + O(t^33)
            sage: pi.scaling_factor()
            0
            sage: pi = EllipticCurveHom_frobenius(E, 3)
            sage: pi.formal()
            t^1331 + O(t^1353)
            sage: pi.scaling_factor()
            0
            sage: pi = EllipticCurveHom_frobenius(E, 0)
            sage: pi == E.scalar_multiplication(1)
            True
            sage: pi.scaling_factor()
            1

        The scaling factor lives in the base ring::

            sage: pi.scaling_factor().parent()
            Finite Field of size 11

        ALGORITHM: Inseparable isogenies of degree `>1` have scaling
        factor `0`.
        """
        if self._degree == 1:
            return self._base_ring.one()
        return self._base_ring.zero()

    def kernel_polynomial(self):
        """
        Return the kernel polynomial of this Frobenius isogeny
        as a polynomial in `x`. This method always returns `1`.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_frobenius import EllipticCurveHom_frobenius
            sage: E = EllipticCurve(GF(11), [1,1])
            sage: pi = EllipticCurveHom_frobenius(E, 5)
            sage: pi.kernel_polynomial()
            1
        """
        return self._poly_ring(1)

    @cached_method
    def dual(self):
        """
        Compute the dual of this Frobenius isogeny.

        This method returns an :class:`EllipticCurveHom` object.

        EXAMPLES:

        An ordinary example::

            sage: from sage.schemes.elliptic_curves.hom_scalar import EllipticCurveHom_scalar
            sage: from sage.schemes.elliptic_curves.hom_frobenius import EllipticCurveHom_frobenius
            sage: E = EllipticCurve(GF(31), [0,1])
            sage: f = EllipticCurveHom_frobenius(E)
            sage: f.dual() * f == EllipticCurveHom_scalar(f.domain(), 31)
            True
            sage: f * f.dual() == EllipticCurveHom_scalar(f.codomain(), 31)
            True

        A supersingular example::

            sage: E = EllipticCurve(GF(31), [1,0])
            sage: f = EllipticCurveHom_frobenius(E)
            sage: f.dual() * f == EllipticCurveHom_scalar(f.domain(), 31)
            True
            sage: f * f.dual() == EllipticCurveHom_scalar(f.codomain(), 31)
            True

        TESTS:

        Some random testing (including small characteristic)::

            sage: p = random_prime(50)
            sage: q = p**randrange(1,10)
            sage: n = randrange(20)
            sage: while True:
            ....:     try:
            ....:         E = EllipticCurve([GF(q).random_element() for _ in range(5)])
            ....:         break
            ....:     except ArithmeticError:
            ....:         pass
            sage: f = EllipticCurveHom_frobenius(E, n)
            sage: f.dual() * f == EllipticCurveHom_scalar(E, p**n)
            True
            sage: f * f.dual() == EllipticCurveHom_scalar(f.codomain(), p**n)
            True
            sage: f.dual().dual() == f  # known bug -- broken in characteristic 2,3
            True
            sage: p in (2,3) or f.dual().dual() == f
            True

        ALGORITHM:

        - For supersingular curves, the dual of Frobenius is again purely
          inseparable, so we start out with a Frobenius isogeny of equal
          degree in the opposite direction.

        - For ordinary curves, we immediately reduce to the case of prime
          degree. The kernel of the dual is the unique subgroup of size `p`,
          which we compute from the `p`-division polynomial.

        In both cases, we then search for the correct post-isomorphism
        using :meth:`find_post_isomorphism`.
        """
        if self._degree == 1:
            return self

        if self._domain.is_supersingular():
            Phi = EllipticCurveHom_frobenius(self._codomain, self._n)

        else:
            E = self._domain
            poly = self._domain.division_polynomial(self._p)
            ker = self._poly_ring(list(poly)[::self._p]).monic()
            Phis = []
            for _ in range(self._n):
                Ep = EllipticCurve([a**self._p for a in E.a_invariants()])
                Phis.append(EllipticCurveIsogeny(Ep, ker, codomain=E))
                E, ker = Ep, ker.map_coefficients(lambda c: c**self._p)
            Phi = EllipticCurveHom_composite.from_factors(Phis[::-1], self._codomain)

        scalar_mul = EllipticCurveHom_scalar(self._domain, self._degree)
        iso = find_post_isomorphism(Phi * self, scalar_mul)
        return iso * Phi

    def inseparable_degree(self):
        """
        Return the inseparable degree of this Frobenius isogeny.

        Since this class implements only purely inseparable isogenies,
        the inseparable degree equals the degree.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_frobenius import EllipticCurveHom_frobenius
            sage: E = EllipticCurve(GF(11), [1,1])
            sage: pi = EllipticCurveHom_frobenius(E, 4)
            sage: pi.inseparable_degree()
            14641
            sage: pi.inseparable_degree() == pi.degree()
            True
        """
        return self._degree
