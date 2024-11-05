# sage.doctest: needs sage.modules sage.rings.finite_rings
r"""
Galois groups of Finite Fields
"""

from sage.groups.abelian_gps.abelian_group_element import AbelianGroupElement
from sage.groups.galois_group import GaloisGroup_cyc
from sage.rings.integer_ring import ZZ
from sage.rings.finite_rings.hom_finite_field import FiniteFieldHomomorphism_generic, FrobeniusEndomorphism_finite_field


class GaloisGroup_GFElement(AbelianGroupElement):
    def as_hom(self):
        r"""
        Return the automorphism of the finite field corresponding to this element.

        EXAMPLES::

            sage: GF(3^6).galois_group()([4]).as_hom()
            Frobenius endomorphism z6 |--> z6^(3^4) on Finite Field in z6 of size 3^6
        """
        n = self.exponents()[0]
        return self.parent()._field.frobenius_endomorphism(n)

    def __call__(self, x):
        r"""
        Return the action of this automorphism on an element `x` of the finite field.

        EXAMPLES::

            sage: k.<a> = GF(3^6)
            sage: g = k.galois_group()([4])
            sage: g(a) == a^(3^4)
            True
        """
        return self.as_hom()(x)

    def fixed_field(self):
        r"""
        The fixed field of this automorphism.

        EXAMPLES::

            sage: k.<a> = GF(3^12)
            sage: g = k.galois_group()([8])
            sage: k0, embed = g.fixed_field()
            sage: k0.cardinality()
            81
            sage: embed.domain() is k0
            True
            sage: embed.codomain() is k
            True
        """
        return self.as_hom().fixed_field()


class GaloisGroup_GF(GaloisGroup_cyc):
    r"""
    The Galois group of a finite field.
    """
    Element = GaloisGroup_GFElement

    def __init__(self, field):
        r"""
        Create a Galois group.

        TESTS::

            sage: TestSuite(GF(9).galois_group()).run()
        """
        GaloisGroup_cyc.__init__(self, field, (field.degree(),), gen_names='Frob')

    def _repr_(self):
        r"""
        String representation of this Galois group.

        EXAMPLES::

            sage: GF(9).galois_group()
            Galois group C2 of GF(3^2)
        """
        return "Galois group C{0} of GF({1}^{0})".format(self._field.degree(), self._field.characteristic())

    def _element_constructor_(self, x, check=True):
        r"""
        Create an element of this Galois group from ``x``.

        INPUT:

        - ``x`` -- one of the following (`G` is this Galois group):

          - the integer 1, denoting the identity of `G`;

          - an element of `G`;

          - a list of length 1, giving the exponent of Frobenius

          - a permutation of the right length that defines an element of `G`,
            or anything that coerces into such a permutation;

          - an automorphism of the finite field.

        - ``check`` -- check that automorphisms have the correct domain and codomain

        EXAMPLES::

            sage: k = GF(3^3)
            sage: G = k.galois_group()
            sage: G(1)
            1
            sage: G([2])
            Frob^2
            sage: G(G.gens()[0])
            Frob
            sage: G([(1,3,2)])                                                          # needs sage.libs.gap
            Frob^2
            sage: G(k.hom(k.gen()^3, k))
            Frob
            sage: G(k.frobenius_endomorphism())
            Frob
        """
        if x == 1:
            return self.element_class(self, [0])

        k = self._field
        d = k.degree()
        n = None
        if isinstance(x, GaloisGroup_GFElement) and x.parent() is self:
            n = x.exponents()[0]
        elif isinstance(x, FiniteFieldHomomorphism_generic):
            if check and not (x.domain() is k and x.codomain() is k):
                raise ValueError("Not an automorphism of the correct finite field")
            a = k.gen()
            b = x(a)
            q = k.base_ring().cardinality()
            n = 0
            while n < d:
                if a == b:
                    break
                n += 1
                a = a**q
            else:
                raise RuntimeError("Automorphism was not a power of Frobenius")
        elif isinstance(x, FrobeniusEndomorphism_finite_field):
            if check and not x.domain() is k:
                raise ValueError("Not an automorphism of the correct finite field")
            n = x.power()
        elif isinstance(x, list) and len(x) == 1 and x[0] in ZZ:
            n = x[0]
        else:
            g = self.permutation_group()(x)
            n = g(1) - 1
        return self.element_class(self, [n])
