# sage_setup: distribution = sagemath-categories
"""
Ring homomorphisms from a polynomial ring to another ring

This module currently implements the canonical ring homomorphism from
`A[x]` to `B[x]` induced by a ring homomorphism from `A` to `B`.

.. TODO::

    Implement homomorphisms from `A[x]` to an arbitrary ring `R`,
    given by a ring homomorphism from `A` to `R` and the image of `x`
    in `R`.

AUTHORS:

- Peter Bruin (March 2014): initial version
"""

from sage.rings.morphism cimport RingHomomorphism_from_base
from sage.structure.element cimport Element

cdef class PolynomialRingHomomorphism_from_base(RingHomomorphism_from_base):
    """
    The canonical ring homomorphism from `R[x]` to `S[x]` induced by a
    ring homomorphism from `R` to `S`.

    EXAMPLES::

        sage: QQ['x'].coerce_map_from(ZZ['x'])
        Ring morphism:
          From: Univariate Polynomial Ring in x over Integer Ring
          To:   Univariate Polynomial Ring in x over Rational Field
          Defn: Induced from base ring by
                Natural morphism:
                  From: Integer Ring
                  To:   Rational Field
    """
    cpdef Element _call_(self, x):
        """
        Evaluate the homomorphism ``self`` at ``x``.

        TESTS::

            sage: from sage.rings.polynomial.polynomial_ring_homomorphism import PolynomialRingHomomorphism_from_base
            sage: R.<x> = ZZ[]
            sage: S = QQ['x']
            sage: f = ZZ.hom(QQ)
            sage: F = PolynomialRingHomomorphism_from_base(R.Hom(S), f)
            sage: F(2*x)
            2*x

            sage: A = PolynomialRing(QQ, 'x', sparse=True)
            sage: B = PolynomialRing(RR, 'x', sparse=True)
            sage: g = QQ.hom(RR)
            sage: G = PolynomialRingHomomorphism_from_base(A.Hom(B), g)
            sage: G(A.gen()^1000000)
            1.0...*x^1000000
        """
        P = self.codomain()
        f = self.underlying_map()
        if P.is_sparse():
            return P({a: f(b) for a, b in x.monomial_coefficients().items()})
        else:
            return P([f(b) for b in x])

    cpdef Element _call_with_args(self, x, args=(), kwds={}):
        """
        Evaluate ``self`` at ``x`` with additional (keyword) arguments.

        TESTS::

            sage: from sage.rings.polynomial.polynomial_ring_homomorphism import PolynomialRingHomomorphism_from_base
            sage: R.<x> = ZZ[]
            sage: S = GF(5)['x']
            sage: f = ZZ.hom(GF(5))
            sage: F = PolynomialRingHomomorphism_from_base(R.Hom(S), f)
            sage: F(2 * x, check=True)
            2*x

            sage: # needs sage.rings.finite_rings
            sage: k = GF(49, 'z')
            sage: A = PolynomialRing(GF(7), 'x', sparse=True)
            sage: B = PolynomialRing(k, 'x', sparse=True)
            sage: g = GF(7).hom(k)
            sage: G = PolynomialRingHomomorphism_from_base(A.Hom(B), g)
            sage: G(A.gen()^1000000, True, construct=False)
            x^1000000
        """
        P = self.codomain()
        f = self.underlying_map()
        if P.is_sparse():
            return P({a: f(b) for a, b in x.monomial_coefficients().items()},
                     *args, **kwds)
        else:
            return P([f(b) for b in x], *args, **kwds)

    def is_injective(self):
        r"""
        Return whether this morphism is injective.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: S.<x> = QQ[]
            sage: R.hom(S).is_injective()
            True
        """
        return self.underlying_map().is_injective()

    def is_surjective(self):
        r"""
        Return whether this morphism is surjective.

        EXAMPLES::

            sage: R.<x> = ZZ[]
            sage: S.<x> = Zmod(2)[]
            sage: R.hom(S).is_surjective()
            True
        """
        return self.underlying_map().is_surjective()
