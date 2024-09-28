"""
Space of homomorphisms between two rings
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.homset import HomsetWithBase
from sage.categories.rings import Rings

_Rings = Rings()

from sage.rings import morphism, quotient_ring


def is_RingHomset(H):
    """
    Return ``True`` if ``H`` is a space of homomorphisms between two rings.

    EXAMPLES::

        sage: from sage.rings.homset import is_RingHomset as is_RH
        sage: is_RH(Hom(ZZ, QQ))
        doctest:warning...
        DeprecationWarning: the function is_RingHomset is deprecated;
        use 'isinstance(..., RingHomset_generic)' instead
        See https://github.com/sagemath/sage/issues/37922 for details.
        True
        sage: is_RH(ZZ)
        False
        sage: is_RH(Hom(RR, CC))                                                        # needs sage.rings.real_mpfr
        True
        sage: is_RH(Hom(FreeModule(ZZ,1), FreeModule(QQ,1)))                            # needs sage.modules
        False
    """
    from sage.misc.superseded import deprecation
    deprecation(37922, "the function is_RingHomset is deprecated; use 'isinstance(..., RingHomset_generic)' instead")
    return isinstance(H, RingHomset_generic)


def RingHomset(R, S, category=None):
    """
    Construct a space of homomorphisms between the rings ``R`` and ``S``.

    For more on homsets, see :func:`Hom()`.

    EXAMPLES::

        sage: Hom(ZZ, QQ) # indirect doctest
        Set of Homomorphisms from Integer Ring to Rational Field
    """
    if isinstance(R, quotient_ring.QuotientRing_nc):
        from .polynomial.polynomial_quotient_ring import PolynomialQuotientRing_generic
        if not isinstance(R, PolynomialQuotientRing_generic):  # backwards compatibility
            return RingHomset_quo_ring(R, S, category=category)
    return RingHomset_generic(R, S, category=category)


class RingHomset_generic(HomsetWithBase):
    """
    A generic space of homomorphisms between two rings.

    EXAMPLES::

        sage: Hom(ZZ, QQ)
        Set of Homomorphisms from Integer Ring to Rational Field
        sage: QQ.Hom(ZZ)
        Set of Homomorphisms from Rational Field to Integer Ring
    """

    Element = morphism.RingHomomorphism

    def __init__(self, R, S, category=None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: Hom(ZZ, QQ)
            Set of Homomorphisms from Integer Ring to Rational Field
        """
        if category is None:
            category = _Rings
        HomsetWithBase.__init__(self, R, S, category)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: Hom(ZZ, QQ) # indirect doctest
            Set of Homomorphisms from Integer Ring to Rational Field
        """
        return "Set of Homomorphisms from %s to %s" % (self.domain(), self.codomain())

    def has_coerce_map_from(self, x):
        """
        The default for coercion maps between ring homomorphism spaces is
        very restrictive (until more implementation work is done).

        Currently this checks if the domains and the codomains are equal.

        EXAMPLES::

            sage: H = Hom(ZZ, QQ)
            sage: H2 = Hom(QQ, ZZ)
            sage: H.has_coerce_map_from(H2)
            False
        """
        return (x.domain() == self.domain() and x.codomain() == self.codomain())

    def _element_constructor_(self, x, check=True, base_map=None):
        """
        Construct an element of ``self`` from ``x``.

        EXAMPLES::

            sage: H = Hom(ZZ, QQ)
            sage: phi = H([1]); phi
            Ring morphism:
              From: Integer Ring
              To:   Rational Field
              Defn: 1 |--> 1
            sage: H2 = Hom(QQ, QQ)
            sage: phi2 = H2(phi); phi2
            Ring endomorphism of Rational Field
              Defn: 1 |--> 1
            sage: H(phi2)
            Ring morphism:
              From: Integer Ring
              To:   Rational Field
              Defn: 1 |--> 1

        You can provide a morphism on the base::

            sage: # needs sage.rings.finite_rings
            sage: k = GF(9)
            sage: z2 = k.gen()
            sage: cc = k.frobenius_endomorphism()
            sage: R.<x> = k[]
            sage: H = Hom(R, R)
            sage: phi = H([x^2], base_map=cc); phi
            Ring endomorphism of Univariate Polynomial Ring in x
             over Finite Field in z2 of size 3^2
              Defn: x |--> x^2
                    with map of base ring
            sage: phi(z2 * x) == z2^3 * x^2
            True

            sage: R.<x> = ZZ[]
            sage: K.<a> = GF(7^2)                                                       # needs sage.rings.finite_rings
            sage: L.<u> = K.extension(x^3 - 3)                                          # needs sage.rings.finite_rings
            sage: phi = L.hom([u^7], base_map=K.frobenius_endomorphism())               # needs sage.rings.finite_rings
            sage: phi(u) == u^7                                                         # needs sage.rings.finite_rings
            True
            sage: phi(a) == a^7                                                         # needs sage.rings.finite_rings
            True

        TESTS::

            sage: H = Hom(ZZ, QQ)
            sage: H == loads(dumps(H))
            True
        """
        from sage.categories.map import Map
        # Case 0: the homomorphism is given by images of generators
        if not (isinstance(x, Map) and x.category_for().is_subcategory(Rings())):
            return morphism.RingHomomorphism_im_gens(self, x, base_map=base_map, check=check)
        if base_map is not None:
            raise ValueError("cannot specify base_map when providing a map")
        # Case 1: the parent fits
        if x.parent() == self:
            if isinstance(x, morphism.RingHomomorphism_im_gens):
                return morphism.RingHomomorphism_im_gens(self, x.im_gens())
            elif isinstance(x, morphism.RingHomomorphism_cover):
                return morphism.RingHomomorphism_cover(self)
            elif isinstance(x, morphism.RingHomomorphism_from_base):
                return morphism.RingHomomorphism_from_base(self, x.underlying_map())
        # Case 2: unique extension via fraction field
        try:
            if (isinstance(x, morphism.RingHomomorphism_im_gens)
                and x.domain().fraction_field().has_coerce_map_from(self.domain())):
                return morphism.RingHomomorphism_im_gens(self, x.im_gens())
        except (TypeError, ValueError):
            pass
        # Case 3: the homomorphism can be extended by coercion
        try:
            return x.extend_codomain(self.codomain()).extend_domain(self.domain())
        except (TypeError, ValueError):
            pass
        # Case 4: the homomorphism is induced from the base ring
        if (self.domain() != self.domain().base()
            or self.codomain() != self.codomain().base()):
            x = self.domain().base().Hom(self.codomain().base())(x)
            return morphism.RingHomomorphism_from_base(self, x)
        raise ValueError('cannot convert {} to an element of {}'.format(x, self))

    def natural_map(self):
        """
        Return the natural map from the domain to the codomain.

        The natural map is the coercion map from the domain ring to the
        codomain ring.

        EXAMPLES::

            sage: H = Hom(ZZ, QQ)
            sage: H.natural_map()
            Natural morphism:
              From: Integer Ring
              To:   Rational Field
        """
        f = self.codomain().coerce_map_from(self.domain())
        if f is None:
            raise TypeError("natural coercion morphism from %s to %s not defined" % (self.domain(), self.codomain()))
        return f

    def zero(self):
        r"""
        Return the zero element of this homset.

        EXAMPLES:

        Since a ring homomorphism maps 1 to 1, there can only be a zero
        morphism when mapping to the trivial ring::

            sage: Hom(ZZ, Zmod(1)).zero()
            Ring morphism:
              From: Integer Ring
              To:   Ring of integers modulo 1
              Defn: 1 |--> 0
            sage: Hom(ZZ, Zmod(2)).zero()
            Traceback (most recent call last):
            ...
            ValueError: homset has no zero element
        """
        if not self.codomain().is_zero():
            raise ValueError("homset has no zero element")
        # there is only one map in this homset
        return self.an_element()


class RingHomset_quo_ring(RingHomset_generic):
    """
    Space of ring homomorphisms where the domain is a (formal) quotient
    ring.

    EXAMPLES::

        sage: R.<x,y> = PolynomialRing(QQ, 2)
        sage: S.<a,b> = R.quotient(x^2 + y^2)                                           # needs sage.libs.singular
        sage: phi = S.hom([b,a]); phi                                                   # needs sage.libs.singular
        Ring endomorphism of Quotient of Multivariate Polynomial Ring in x, y
         over Rational Field by the ideal (x^2 + y^2)
          Defn: a |--> b
                b |--> a
        sage: phi(a)                                                                    # needs sage.libs.singular
        b
        sage: phi(b)                                                                    # needs sage.libs.singular
        a

    TESTS:

    We test pickling of a homset from a quotient.

    ::

        sage: # needs sage.libs.singular
        sage: R.<x,y> = PolynomialRing(QQ, 2)
        sage: S.<a,b> = R.quotient(x^2 + y^2)
        sage: H = S.Hom(R)
        sage: H == loads(dumps(H))
        True

    We test pickling of actual homomorphisms in a quotient::

        sage: phi = S.hom([b,a])                                                        # needs sage.libs.singular
        sage: phi == loads(dumps(phi))                                                  # needs sage.libs.singular
        True
    """

    Element = morphism.RingHomomorphism_from_quotient

    def _element_constructor_(self, x, base_map=None, check=True):
        """
        Construct an element of ``self`` from ``x``.

        EXAMPLES::

            sage: # needs sage.libs.singular
            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: S.<a,b> = R.quotient(x^2 + y^2)
            sage: H = S.Hom(R)
            sage: phi = H([b, a]); phi
            Ring morphism:
              From: Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2 + y^2)
              To:   Multivariate Polynomial Ring in x, y over Rational Field
              Defn: a |--> b
                    b |--> a
            sage: R2.<x,y> = PolynomialRing(ZZ, 2)
            sage: H2 = Hom(R2, S)
            sage: H2(phi)
            Composite map:
              From: Multivariate Polynomial Ring in x, y over Integer Ring
              To:   Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2 + y^2)
              Defn:   Coercion map:
                      From: Multivariate Polynomial Ring in x, y over Integer Ring
                      To:   Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2 + y^2)
                    then
                      Ring morphism:
                      From: Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2 + y^2)
                      To:   Multivariate Polynomial Ring in x, y over Rational Field
                      Defn: a |--> b
                            b |--> a
                    then
                      Coercion map:
                      From: Multivariate Polynomial Ring in x, y over Rational Field
                      To:   Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2 + y^2)
        """
        if isinstance(x, morphism.RingHomomorphism_from_quotient):
            phi = x._phi()
        else:
            pi = self.domain().cover()
            phi = pi.domain().hom(x, base_map=base_map, check=check)
        return self.element_class(self, phi)
