# sage_setup: distribution = sagemath-schemes
r"""
Sums of morphisms of elliptic curves

The set `\mathrm{Hom}(E,E')` of morphisms between two elliptic curves
forms an abelian group under pointwise addition. An important special
case is the endomorphism ring `\mathrm{End}(E) = \mathrm{Hom}(E,E)`.
However, it is not immediately obvious how to compute some properties
of the sum `\varphi+\psi` of two isogenies, even when both are given
explicitly. This class provides functionality for representing sums of
elliptic-curve morphisms (in particular, isogenies and endomorphisms)
formally, and explicitly computing important properties (such as the
degree or the kernel polynomial) from the formal representation.

EXAMPLES::

    sage: E = EllipticCurve(GF(101), [5,5])
    sage: phi = E.isogenies_prime_degree(7)[0]
    sage: phi + phi
    Sum morphism:
      From: Elliptic Curve defined by y^2 = x^3 + 5*x + 5 over Finite Field of size 101
      To:   Elliptic Curve defined by y^2 = x^3 + 12*x + 98 over Finite Field of size 101
      Via:  (Isogeny of degree 7 from Elliptic Curve defined by y^2 = x^3 + 5*x + 5 over Finite Field of size 101 to Elliptic Curve defined by y^2 = x^3 + 12*x + 98 over Finite Field of size 101, Isogeny of degree 7 from Elliptic Curve defined by y^2 = x^3 + 5*x + 5 over Finite Field of size 101 to Elliptic Curve defined by y^2 = x^3 + 12*x + 98 over Finite Field of size 101)
    sage: phi + phi == phi * E.scalar_multiplication(2)
    True
    sage: phi + phi + phi == phi * E.scalar_multiplication(3)
    True

An example of computing with a supersingular endomorphism ring::

    sage: E = EllipticCurve(GF(419^2), [1,0])
    sage: i = E.automorphisms()[-1]
    sage: j = E.frobenius_isogeny()
    sage: i * j == - j * i                          # i,j anticommute
    True
    sage: (i + j) * i == i^2 - i*j                  # distributive law
    True
    sage: (j - E.scalar_multiplication(1)).degree() # point counting!
    420

AUTHORS:

- Lorenz Panny (2023)
"""

from sage.misc.cachefunc import cached_method
from sage.structure.sequence import Sequence

from sage.arith.misc import gcd

from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring import polygen

from sage.sets.primes import Primes

from sage.schemes.elliptic_curves.ell_field import point_of_order
from sage.groups.generic import discrete_log, order_from_multiple

from sage.schemes.elliptic_curves.hom import EllipticCurveHom, compare_via_evaluation


class EllipticCurveHom_sum(EllipticCurveHom):

    _degree = None
    _phis = None

    def __init__(self, phis, domain=None, codomain=None):
        r"""
        Construct a sum morphism of elliptic curves from its summands.
        (For empty sums, the domain and codomain curves must be given.)

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.hom_sum import EllipticCurveHom_sum
            sage: E = EllipticCurve(GF(101), [5,5])
            sage: phi = E.isogenies_prime_degree(7)[0]
            sage: EllipticCurveHom_sum([phi, phi])
            Sum morphism:
              From: Elliptic Curve defined by y^2 = x^3 + 5*x + 5 over Finite Field of size 101
              To:   Elliptic Curve defined by y^2 = x^3 + 12*x + 98 over Finite Field of size 101
              Via:  (Isogeny of degree 7 from Elliptic Curve defined by y^2 = x^3 + 5*x + 5 over Finite Field of size 101 to Elliptic Curve defined by y^2 = x^3 + 12*x + 98 over Finite Field of size 101, Isogeny of degree 7 from Elliptic Curve defined by y^2 = x^3 + 5*x + 5 over Finite Field of size 101 to Elliptic Curve defined by y^2 = x^3 + 12*x + 98 over Finite Field of size 101)

        The zero morphism can be constructed even between non-isogenous curves::

            sage: E1 = EllipticCurve(GF(101), [5,5])
            sage: E2 = EllipticCurve(GF(101), [7,7])
            sage: E1.is_isogenous(E2)
            False
            sage: EllipticCurveHom_sum([], E1, E2)
            Sum morphism:
              From: Elliptic Curve defined by y^2 = x^3 + 5*x + 5 over Finite Field of size 101
              To:   Elliptic Curve defined by y^2 = x^3 + 7*x + 7 over Finite Field of size 101
              Via:  ()
        """
        phis = tuple(phis)

        if not phis and (domain is None or codomain is None):
            raise ValueError('need either phis or both domain and codomain')

        for phi in phis:
            if not isinstance(phi, EllipticCurveHom):
                raise ValueError(f'not an elliptic-curve morphism: {phi}')

        if domain is None:
            domain = phis[0].domain()
        if codomain is None:
            codomain = phis[0].codomain()
        for phi in phis:
            if phi.domain() != domain:
                raise ValueError(f'summand {phi} has incorrect domain (need {domain})')
            if phi.codomain() != codomain:
                raise ValueError(f'summand {phi} has incorrect codomain (need {codomain})')

        self._phis = phis
        self._domain = domain
        self._codomain = codomain

        # We temporarily overwrite the _degree attribute here to prevent the
        # EllipticCurveHom constructor from attempting to compute the degree.
        self._degree = 0
        EllipticCurveHom.__init__(self, self._domain, self._codomain)
        self._degree = None

    def _call_(self, P):
        r"""
        Evaluate this sum morphism at a point.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [5,5])
            sage: phi = E.isogenies_prime_degree(7)[0]
            sage: P = E.lift_x(0)
            sage: (phi + phi)(P)
            (72 : 56 : 1)
            sage: (phi - phi)(P)
            (0 : 1 : 0)
        """
        return sum((phi(P) for phi in self._phis), self._codomain(0))

    def _eval(self, P):
        r"""
        Less strict evaluation method for internal use.

        In particular, this can be used to evaluate ``self`` at a
        point defined over an extension field.

        INPUT: a sequence of 3 coordinates defining a point on ``self``

        OUTPUT: the result of evaluating ``self`` at the given point

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [5,5])
            sage: phi = E.isogenies_prime_degree(7)[0]
            sage: P = E.change_ring(GF(101^2)).lift_x(1)
            sage: (phi + phi)._eval(P)
            (11 : 15*z2 + 71 : 1)
            sage: (phi - phi)._eval(P)
            (0 : 1 : 0)
        """
        if self._domain.defining_polynomial()(*P):
            raise ValueError(f'{P} not on {self._domain}')
        k = Sequence(P).universe()
        return sum((phi._eval(P) for phi in self._phis), self._codomain.base_extend(k)(0))

    def _repr_(self):
        r"""
        Return basic facts about this sum morphism as a string.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [5,5])
            sage: phi = E.isogenies_prime_degree(7)[0]
            sage: phi + phi  # indirect doctest
            Sum morphism:
              From: Elliptic Curve defined by y^2 = x^3 + 5*x + 5 over Finite Field of size 101
              To:   Elliptic Curve defined by y^2 = x^3 + 12*x + 98 over Finite Field of size 101
              Via:  (Isogeny of degree 7 from Elliptic Curve defined by y^2 = x^3 + 5*x + 5 over Finite Field of size 101 to Elliptic Curve defined by y^2 = x^3 + 12*x + 98 over Finite Field of size 101, Isogeny of degree 7 from Elliptic Curve defined by y^2 = x^3 + 5*x + 5 over Finite Field of size 101 to Elliptic Curve defined by y^2 = x^3 + 12*x + 98 over Finite Field of size 101)
        """
        return f'Sum morphism:' \
                f'\n  From: {self._domain}' \
                f'\n  To:   {self._codomain}' \
                f'\n  Via:  {self._phis}'

    def summands(self):
        r"""
        Return the individual summands making up this sum morphism.

        EXAMPLES::

            sage: E = EllipticCurve(j=5)
            sage: m2 = E.scalar_multiplication(2)
            sage: m3 = E.scalar_multiplication(3)
            sage: m2 + m3
            Sum morphism:
              From: Elliptic Curve defined by y^2 + x*y = x^3 + x^2 + 180*x + 17255 over Rational Field
              To:   Elliptic Curve defined by y^2 + x*y = x^3 + x^2 + 180*x + 17255 over Rational Field
              Via:  (Scalar-multiplication endomorphism [2] of Elliptic Curve defined by y^2 + x*y = x^3 + x^2 + 180*x + 17255 over Rational Field, Scalar-multiplication endomorphism [3] of Elliptic Curve defined by y^2 + x*y = x^3 + x^2 + 180*x + 17255 over Rational Field)
        """
        return self._phis

    @cached_method
    def to_isogeny_chain(self):
        r"""
        Convert this formal sum of elliptic-curve morphisms into a
        :class:`~sage.schemes.elliptic_curves.hom_composite.EllipticCurveHom_composite`
        object representing the same morphism.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [5,5])
            sage: phi = E.isogenies_prime_degree(7)[0]
            sage: (phi + phi).to_isogeny_chain()
            Composite morphism of degree 28 = 4*1*7:
              From: Elliptic Curve defined by y^2 = x^3 + 5*x + 5 over Finite Field of size 101
              To:   Elliptic Curve defined by y^2 = x^3 + 12*x + 98 over Finite Field of size 101

        ::

            sage: p = 419
            sage: E = EllipticCurve(GF(p^2), [1,0])
            sage: iota = E.automorphisms()[2]   # sqrt(-1)
            sage: pi = E.frobenius_isogeny()    # sqrt(-p)
            sage: endo = iota + pi
            sage: endo.degree()
            420
            sage: endo.to_isogeny_chain()
            Composite morphism of degree 420 = 4*1*3*5*7:
              From: Elliptic Curve defined by y^2 = x^3 + x over Finite Field in z2 of size 419^2
              To:   Elliptic Curve defined by y^2 = x^3 + x over Finite Field in z2 of size 419^2

        The decomposition is impossible for the constant zero map::

             sage: endo = iota*pi + pi*iota
             sage: endo.degree()
             0
             sage: endo.to_isogeny_chain()
             Traceback (most recent call last):
             ...
             ValueError: zero morphism cannot be written as a composition of isogenies

        Isomorphisms are supported as well::

            sage: E = EllipticCurve(j=5); E
            Elliptic Curve defined by y^2 + x*y = x^3 + x^2 + 180*x + 17255 over Rational Field
            sage: m2 = E.scalar_multiplication(2)
            sage: m3 = E.scalar_multiplication(3)
            sage: (m2 - m3).to_isogeny_chain()
            Composite morphism of degree 1 = 1^2:
              From: Elliptic Curve defined by y^2 + x*y = x^3 + x^2 + 180*x + 17255 over Rational Field
              To:   Elliptic Curve defined by y^2 + x*y = x^3 + x^2 + 180*x + 17255 over Rational Field
            sage: (m2 - m3).rational_maps()
            (x, -x - y)
        """
        deg = self.degree()
        if deg.is_zero():
            raise ValueError('zero morphism cannot be written as a composition of isogenies')

        p = self.base_ring().characteristic()
        insep = self.inseparable_degree().valuation(p) if p else 0

        scalar = 1  #TODO Can we detect scalar factors earlier to save some extensions below?

        ker = []
        for l,m in deg.factor():
            if l == p:  # possibly inseparable
                if insep < m:
                    # kernel of the separable p-power part is unique
                    P = point_of_order(self.domain(), p**(m-insep))
                    ker.append(P)
                continue

#            F = self.domain().division_field(l**m)  #FIXME this can be used once #35936 is done; workaround below
            F = self.domain().division_polynomial(l**m).splitting_field('X').extension(2,'Y')

            P,Q = self.domain().change_ring(F).torsion_basis(l**m)
            if self.is_endomorphism():
                R,S = P,Q
            else:
                R,S = self.codomain().change_ring(F).torsion_basis(l**m)
            M = self.matrix_on_subgroup((P,Q), (R,S))
            g = ZZ(gcd(M.list())).p_primary_part(l)
            if g > 1:
                scalar *= g
                M = (M.change_ring(ZZ) / g).change_ring(M.base_ring())
            K = M.left_kernel_matrix()
            for row in K:
                u,v = map(ZZ, row)
                pt = u*P + v*Q
                pt.set_order(row.additive_order())
                ker.append(pt)

        from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
        phi = EllipticCurveHom_composite(self.domain(), [])

        if scalar != 1:
            phi *= phi.codomain().scalar_multiplication(scalar)

        while ker:
            K = ker.pop(0)

            (l,e), = K.order().factor()
            for i in reversed(range(e)):
                Kl = l**i * K
                Kl.set_order(l)

                from sage.groups.generic import multiples
                from sage.misc.misc_c import prod
                x = polygen(Kl.base_ring())
                poly = prod(x - T.x() for T in multiples(Kl, l//2, Kl))
                poly = poly.change_ring(self.base_ring())

                psi = phi.codomain().isogeny(poly)
                phi = psi * phi
                K = psi._eval(K)
                ker = [psi._eval(P) for P in ker]

        if insep:
            frob = phi.codomain().frobenius_isogeny(insep)
            phi = frob * phi

        from sage.schemes.elliptic_curves.hom import find_post_isomorphism
        iso = find_post_isomorphism(phi, self)
        return iso * phi

    # EllipticCurveHom methods

    def _degree_bounds(self):
        r"""
        Return a lower and upper bound on the degree of this sum morphism.

        EXAMPLES::

            sage: E = EllipticCurve(GF(307), [5,5])
            sage: phi = E.isogenies_prime_degree(3)[0]; phi
            Isogeny of degree 3 from Elliptic Curve defined by y^2 = x^3 + 5*x + 5 over Finite Field of size 307 to Elliptic Curve defined by y^2 = x^3 + 227*x + 163 over Finite Field of size 307
            sage: psi = next(iso*psi for psi in E.isogenies_prime_degree(43)
            ....:                    for iso in psi.codomain().isomorphisms(phi.codomain())); psi
            Isogeny of degree 43 from Elliptic Curve defined by y^2 = x^3 + 5*x + 5 over Finite Field of size 307 to Elliptic Curve defined by y^2 = x^3 + 227*x + 163 over Finite Field of size 307
            sage: (phi + psi)._degree_bounds()
            (24, 68)
            sage: (phi + psi).degree()
            31
            sage: (phi - phi)._degree_bounds()
            (0, 12)
            sage: (phi - phi).degree()
            0

        ::

            sage: E = EllipticCurve(GF(443), [1,1])
            sage: pi = E.frobenius_endomorphism()
            sage: m1 = E.scalar_multiplication(1)
            sage: (pi - m1)._degree_bounds()
            (402, 486)
            sage: (pi - m1)._degree_bounds() == Hasse_bounds(443)
            True
            sage: (pi - m1).degree()
            433

        ALGORITHM: Repeated application of the Cauchy-Schwarz inequality,
        here in the form
        `|\deg(f+g) - \deg(f) - \deg(g)| \leq 2\sqrt{\deg(f)\cdot\deg(g)}`.
        See for instance Lemma V.1.2 of [Sil2009]_.
        """
        lo, hi = ZZ.zero(), ZZ.zero()
        for phi in self._phis:
            m = (hi * phi.degree()).isqrt()
            hi += phi.degree() + 2*m
            lo += phi.degree() - 2*m
            lo = max(lo, 0)
        return lo, hi

    def _compute_degree(self):
        r"""
        Internal method to compute and cache the degree of this sum morphism
        (and its dual).

        ALGORITHM: Evaluate the composition with the dual on points of small
        order and solve logarithms to eventually recover the degree using CRT.
        (This is essentially Schoof's algorithm, applied to a scalar.)

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [5,5])
            sage: phi = E.isogenies_prime_degree(7)[0]
            sage: isog = phi + phi
            sage: print(isog._degree)
            None
            sage: isog._compute_degree()
            sage: isog._degree
            28

        ::

            sage: E = EllipticCurve(GF(443), [1,1])
            sage: pi = E.frobenius_endomorphism()
            sage: m1 = E.scalar_multiplication(1)
            sage: endo = pi - m1
            sage: print(endo._degree)
            None
            sage: endo._compute_degree()
            sage: endo._degree
            433
            sage: endo.dual()._degree
            433
        """
        if self._degree is not None:
            return
        if not self._phis:
            self._degree = ZZ.zero()
        elif len(self._phis) == 1:
            self._degree = self._phis[0].degree()
        else:
            #TODO In some cases it would probably be faster to simply
            # compute the kernel polynomial using the addition formulas?
            from sage.rings.finite_rings.integer_mod import Mod

            lo, hi = self._degree_bounds()
            M = hi - lo + 1
            rem = Mod(0,1)
            for l in Primes():
                if rem.modulus() >= M:
                    break
                try:
                    P = point_of_order(self._domain, l)
                except ValueError:
                    continue   # supersingular and l == p

                Q = self.dual()._eval(self._eval(P))
                d = discrete_log(Q, P, ord=l, operation='+')
                rem = rem.crt(Mod(d-lo, l))

            self._degree = lo + rem.lift()
            self.dual()._degree = self._degree

    @staticmethod
    def _comparison_impl(left, right, op):
        r"""
        Compare a sum morphism to another elliptic-curve morphism.

        Called by :meth:`EllipticCurveHom._richcmp_`.

        If possible, we use
        :func:`~sage.schemes.elliptic_curves.hom.compare_via_evaluation`.
        The complexity in that case is polynomial in the logarithm of
        the degree.

        TESTS::

            sage: from sage.schemes.elliptic_curves.hom_sum import EllipticCurveHom_sum
            sage: E = EllipticCurve(GF(419^2), [1,0])
            sage: i = E.automorphisms()[-1]
            sage: j = E.frobenius_isogeny()
            sage: i + j == j + i
            True
        """
        from sage.structure.richcmp import op_EQ
        if op != op_EQ:
            return NotImplemented
        try:
            return compare_via_evaluation(left, right)
        except NotImplementedError:
            return NotImplemented

    def degree(self):
        r"""
        Return the degree of this sum morphism.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [5,5])
            sage: phi = E.isogenies_prime_degree(7)[0]
            sage: (phi + phi).degree()
            28

        This method yields a simple toy point-counting algorithm::

            sage: E = EllipticCurve(GF(101), [5,5])
            sage: m1 = E.scalar_multiplication(1)
            sage: pi = E.frobenius_endomorphism()
            sage: (pi - m1).degree()
            119
            sage: E.count_points()
            119

        ALGORITHM: Essentially Schoof's algorithm; see :meth:`_compute_degree`.
        """
        if self._degree is None:
            self._compute_degree()
        return self._degree

    def rational_maps(self):
        r"""
        Return the rational maps of this sum morphism.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [5,5])
            sage: phi = E.isogenies_prime_degree(7)[0]
            sage: (phi + phi).rational_maps()
            ((5*x^28 + 43*x^27 + 26*x^26 - ... + 7*x^2 - 23*x + 38)/(23*x^27 + 16*x^26 + 9*x^25 + ... - 43*x^2 - 22*x + 37),
             (42*x^42*y - 44*x^41*y - 22*x^40*y + ... - 26*x^2*y - 50*x*y - 18*y)/(-24*x^42 - 47*x^41 - 12*x^40 + ... + 18*x^2 - 48*x + 18))

        ALGORITHM: :meth:`to_isogeny_chain`.
        """
        #TODO In some cases it would probably be faster to compute this
        # directly using the addition formulas?
        return self.to_isogeny_chain().rational_maps()

    def x_rational_map(self):
        r"""
        Return the `x`-coordinate rational map of this sum morphism.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [5,5])
            sage: phi = E.isogenies_prime_degree(7)[0]
            sage: (phi + phi).x_rational_map()
            (9*x^28 + 37*x^27 + 67*x^26 + ... + 53*x^2 + 100*x + 28)/(x^27 + 49*x^26 + 97*x^25 + ... + 64*x^2 + 21*x + 6)

        ALGORITHM: :meth:`to_isogeny_chain`.
        """
        #TODO In some cases it would probably be faster to compute this
        # directly using the addition formulas?
        return self.to_isogeny_chain().x_rational_map()

    def kernel_polynomial(self):
        r"""
        Return the kernel polynomial of this sum morphism.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [5,5])
            sage: phi = E.isogenies_prime_degree(7)[0]
            sage: (phi + phi).kernel_polynomial()
            x^15 + 75*x^14 + 16*x^13 + 59*x^12 + 28*x^11 + 60*x^10 + 69*x^9 + 79*x^8 + 79*x^7 + 52*x^6 + 35*x^5 + 11*x^4 + 37*x^3 + 69*x^2 + 66*x + 63

        ::

            sage: E = EllipticCurve(GF(11), [5,5])
            sage: pi = E.frobenius_endomorphism()
            sage: m1 = E.scalar_multiplication(1)
            sage: (pi - m1).kernel_polynomial()
            x^9 + 7*x^8 + 2*x^7 + 4*x^6 + 10*x^4 + 4*x^3 + 9*x^2 + 7*x

        ALGORITHM: :meth:`to_isogeny_chain`.
        """
        #TODO In some cases it would probably be faster to compute this
        # directly using the addition formulas?
        return self.to_isogeny_chain().kernel_polynomial()

    def scaling_factor(self):
        r"""
        Return the Weierstrass scaling factor associated to this
        sum morphism.

        The scaling factor is the constant `u` (in the base field)
        such that `\varphi^* \omega_2 = u \omega_1`, where
        `\varphi: E_1\to E_2` is this morphism and `\omega_i` are
        the standard Weierstrass differentials on `E_i` defined by
        `\mathrm dx/(2y+a_1x+a_3)`.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [5,5])
            sage: phi = E.isogenies_prime_degree(7)[0]
            sage: phi.scaling_factor()
            84
            sage: (phi + phi).scaling_factor()
            67

        ALGORITHM: The scaling factor is additive under addition
        of elliptic-curve morphisms, so we simply add together the
        scaling factors of the :meth:`summands`.
        """
        return sum(phi.scaling_factor() for phi in self._phis)

    @cached_method
    def dual(self):
        r"""
        Return the dual of this sum morphism.

        EXAMPLES::

            sage: E = EllipticCurve(GF(101), [5,5])
            sage: phi = E.isogenies_prime_degree(7)[0]
            sage: (phi + phi).dual()
            Sum morphism:
              From: Elliptic Curve defined by y^2 = x^3 + 12*x + 98 over Finite Field of size 101
              To:   Elliptic Curve defined by y^2 = x^3 + 5*x + 5 over Finite Field of size 101
              Via:  (Isogeny of degree 7 from Elliptic Curve defined by y^2 = x^3 + 12*x + 98 over Finite Field of size 101 to Elliptic Curve defined by y^2 = x^3 + 5*x + 5 over Finite Field of size 101, Isogeny of degree 7 from Elliptic Curve defined by y^2 = x^3 + 12*x + 98 over Finite Field of size 101 to Elliptic Curve defined by y^2 = x^3 + 5*x + 5 over Finite Field of size 101)
            sage: (phi + phi).dual() == phi.dual() + phi.dual()
            True

        ::

            sage: E = EllipticCurve(GF(431^2), [1,0])
            sage: iota = E.automorphisms()[2]
            sage: m2 = E.scalar_multiplication(2)
            sage: endo = m2 + iota
            sage: endo.dual()
            Sum morphism:
              From: Elliptic Curve defined by y^2 = x^3 + x over Finite Field in z2 of size 431^2
              To:   Elliptic Curve defined by y^2 = x^3 + x over Finite Field in z2 of size 431^2
              Via:  (Scalar-multiplication endomorphism [2] of Elliptic Curve defined by y^2 = x^3 + x over Finite Field in z2 of size 431^2, Elliptic-curve endomorphism of Elliptic Curve defined by y^2 = x^3 + x over Finite Field in z2 of size 431^2
              Via:  (u,r,s,t) = (8*z2 + 427, 0, 0, 0))
            sage: endo.dual() == (m2 - iota)
            True

        ALGORITHM: Taking the dual distributes over addition.
        """
        psi = EllipticCurveHom_sum((phi.dual() for phi in self._phis),
                                   domain=self._codomain, codomain=self._domain)
        psi._degree = self._degree
        if self.trace.is_in_cache():
            psi.trace.set_cache(-self.trace.cache)
        psi.dual.set_cache(self)
        return psi

    @cached_method
    def inseparable_degree(self):
        r"""
        Compute the inseparable degree of this sum morphism.

        EXAMPLES::

            sage: E = EllipticCurve(GF(7), [0,1])
            sage: m3 = E.scalar_multiplication(3)
            sage: m3.inseparable_degree()
            1
            sage: m4 = E.scalar_multiplication(4)
            sage: m7 = m3 + m4; m7
            Sum morphism:
              From: Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 7
              To:   Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 7
              Via:  (Scalar-multiplication endomorphism [3] of Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 7, Scalar-multiplication endomorphism [4] of Elliptic Curve defined by y^2 = x^3 + 1 over Finite Field of size 7)
            sage: m7.degree()
            49
            sage: m7.inseparable_degree()
            7

        A supersingular example::

            sage: E = EllipticCurve(GF(7), [1,0])
            sage: m3 = E.scalar_multiplication(3)
            sage: m3.inseparable_degree()
            1
            sage: m4 = E.scalar_multiplication(4)
            sage: m7 = m3 + m4; m7
            Sum morphism:
              From: Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 7
              To:   Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 7
              Via:  (Scalar-multiplication endomorphism [3] of Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 7, Scalar-multiplication endomorphism [4] of Elliptic Curve defined by y^2 = x^3 + x over Finite Field of size 7)
            sage: m7.inseparable_degree()
            49
        """
        if self.is_zero():
            raise ValueError('zero morphism is not an isogeny')

        p = self.base_ring().characteristic()
        if not p:
            return ZZ.one()

        m = self.degree().valuation(p)
        if not m:
            return ZZ.one()

        try:
            P = point_of_order(self.domain(), p**m)
        except ValueError:
            # supersingular; every p-isogeny is purely inseparable
            return p**m

        Q = self._eval(P)
        return order_from_multiple(Q, p**m)
