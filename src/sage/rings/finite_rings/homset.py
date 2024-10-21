# sage_setup: distribution = sagemath-categories
# sage.doctest: needs sage.rings.finite_rings
"""
Homset for finite fields

This is the set of all field homomorphisms between two finite fields.

EXAMPLES::

    sage: R.<t> = ZZ[]
    sage: E.<a> = GF(25, modulus = t^2 - 2)
    sage: F.<b> = GF(625)
    sage: H = Hom(E, F)
    sage: f = H([4*b^3 + 4*b^2 + 4*b]); f
    Ring morphism:
      From: Finite Field in a of size 5^2
      To:   Finite Field in b of size 5^4
      Defn: a |--> 4*b^3 + 4*b^2 + 4*b
    sage: f(2)
    2
    sage: f(a)
    4*b^3 + 4*b^2 + 4*b
    sage: len(H)
    2
    sage: [phi(2*a)^2 for phi in Hom(E, F)]
    [3, 3]

We can also create endomorphisms::

    sage: End(E)
    Automorphism group of Finite Field in a of size 5^2
    sage: End(GF(7))[0]
    Ring endomorphism of Finite Field of size 7
      Defn: 1 |--> 1
    sage: H = Hom(GF(7), GF(49, 'c'))
    sage: H[0](2)
    2
"""

from sage.rings.homset import RingHomset_generic
from sage.rings.finite_rings.hom_finite_field import FiniteFieldHomomorphism_generic
from sage.rings.finite_rings.finite_field_base import FiniteField
from sage.rings.integer import Integer
from sage.rings.morphism import RingHomomorphism_im_gens
from sage.structure.sequence import Sequence

class FiniteFieldHomset(RingHomset_generic):
    """
    Set of homomorphisms with domain a given finite field.
    """
#     def __init__(self, R, S, category=None):
#         if category is None:
#             from sage.categories.finite_fields import FiniteFields
#             category = FiniteFields()
#         RingHomset_generic.__init__(self, R, S, category)

    def __call__(self, im_gens, base_map=None, check=True):
        """
        Construct the homomorphism defined by ``im_gens``.

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: E.<a> = GF(25, modulus = t^2 - 2)
            sage: F.<b> = GF(625)
            sage: End(E)
            Automorphism group of Finite Field in a of size 5^2
            sage: list(Hom(E, F))
            [Ring morphism:
              From: Finite Field in a of size 5^2
              To:   Finite Field in b of size 5^4
              Defn: a |--> 4*b^3 + 4*b^2 + 4*b,
             Ring morphism:
              From: Finite Field in a of size 5^2
              To:   Finite Field in b of size 5^4
              Defn: a |--> b^3 + b^2 + b]
            sage: [phi(2*a)^2 for phi in Hom(E, F)]
            [3, 3]
            sage: End(GF(7))[0]
            Ring endomorphism of Finite Field of size 7
              Defn: 1 |--> 1
            sage: H = Hom(GF(7), GF(49, 'c'))
            sage: H[0](2)
            2
            sage: Hom(GF(49, 'c'), GF(7)).list()
            []
            sage: Hom(GF(49, 'c'), GF(81, 'd')).list()
            []
            sage: H = Hom(GF(9, 'a'), GF(81, 'b'))
            sage: H == loads(dumps(H))
            True
        """
        if isinstance(im_gens, FiniteFieldHomomorphism_generic):
            if base_map is not None:
                raise ValueError("Cannot specify base map when providing morphism")
            return self._coerce_impl(im_gens)
        try:
            if self.domain().degree() == 1:
                from sage.rings.finite_rings.hom_prime_finite_field import FiniteFieldHomomorphism_prime
                return FiniteFieldHomomorphism_prime(self, im_gens, base_map=base_map, check=check)
            if isinstance(self.codomain(), FiniteField):
                return FiniteFieldHomomorphism_generic(self, im_gens, base_map=base_map, check=check)
            # Currently, FiniteFieldHomomorphism_generic does not work if
            # the codomain is not derived from the finite field base class;
            # in that case, we have to fall back to the generic
            # implementation for rings
            else:
                return RingHomomorphism_im_gens(self, im_gens, base_map=base_map, check=check)
        except (NotImplementedError, ValueError):
            try:
                return self._coerce_impl(im_gens)
            except TypeError:
                raise TypeError("images do not define a valid homomorphism")

    def _coerce_impl(self, x):
        """
        Coercion of other morphisms.

        EXAMPLES::

            sage: k.<a> = GF(25)
            sage: l.<b> = GF(625)
            sage: H = Hom(k, l)
            sage: G = loads(dumps(H))
            sage: H is G
            True
            sage: G.coerce(list(H)[0]) # indirect doctest
            Ring morphism:
              From: Finite Field in a of size 5^2
              To:   Finite Field in b of size 5^4
              Defn: a |--> 4*b^3 + 4*b^2 + 4*b + 3
        """
        if not isinstance(x, FiniteFieldHomomorphism_generic):
            raise TypeError
        if x.parent() is self:
            return x
        if x.parent() == self:
            return FiniteFieldHomomorphism_generic(self, x.im_gens())
        raise TypeError

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: Hom(GF(4, 'a'), GF(16, 'b'))._repr_()
            'Set of field embeddings from Finite Field in a of size 2^2 to Finite Field in b of size 2^4'
            sage: Hom(GF(4, 'a'), GF(4, 'c'))._repr_()
            'Set of field embeddings from Finite Field in a of size 2^2 to Finite Field in c of size 2^2'
            sage: Hom(GF(4, 'a'), GF(4, 'a'))._repr_()
            'Automorphism group of Finite Field in a of size 2^2'
        """
        D = self.domain()
        C = self.codomain()
        if C == D:
            return "Automorphism group of %s" % D
        else:
            return "Set of field embeddings from %s to %s" % (D, C)

    def is_aut(self):
        """
        Check if ``self`` is an automorphism.

        EXAMPLES::

            sage: Hom(GF(4, 'a'), GF(16, 'b')).is_aut()
            False
            sage: Hom(GF(4, 'a'), GF(4, 'c')).is_aut()
            False
            sage: Hom(GF(4, 'a'), GF(4, 'a')).is_aut()
            True
        """
        return self.domain() == self.codomain()

    def order(self):
        """
        Return the order of this set of field homomorphisms.

        EXAMPLES::

            sage: K.<a> = GF(125)
            sage: End(K)
            Automorphism group of Finite Field in a of size 5^3
            sage: End(K).order()
            3
            sage: L.<b> = GF(25)
            sage: Hom(L, K).order() == Hom(K, L).order() == 0
            True
        """
        try:
            return self.__order
        except AttributeError:
            pass
        n = len(self.list())
        self.__order = n
        return n

    def __len__(self):
        """
        Return the number of elements of ``self``.

        EXAMPLES::

            sage: K.<a> = GF(25)
            sage: len(End(K))
            2
        """
        return self.order()

    def list(self):
        """
        Return a list of all the elements in this set of field homomorphisms.

        EXAMPLES::

            sage: K.<a> = GF(25)
            sage: End(K)
            Automorphism group of Finite Field in a of size 5^2
            sage: list(End(K))
            [Ring endomorphism of Finite Field in a of size 5^2
              Defn: a |--> 4*a + 1,
             Ring endomorphism of Finite Field in a of size 5^2
              Defn: a |--> a]
            sage: L.<z> = GF(7^6)
            sage: [g for g in End(L) if (g^3)(z) == z]
            [Ring endomorphism of Finite Field in z of size 7^6
              Defn: z |--> z,
             Ring endomorphism of Finite Field in z of size 7^6
              Defn: z |--> 5*z^4 + 5*z^3 + 4*z^2 + 3*z + 1,
             Ring endomorphism of Finite Field in z of size 7^6
              Defn: z |--> 3*z^5 + 5*z^4 + 5*z^2 + 2*z + 3]

        Between isomorphic fields with different moduli::

            sage: k1 = GF(1009)
            sage: k2 = GF(1009, modulus='primitive')
            sage: Hom(k1, k2).list()
            [
            Ring morphism:
              From: Finite Field of size 1009
              To:   Finite Field of size 1009
              Defn: 1 |--> 1
            ]
            sage: Hom(k2, k1).list()
            [
            Ring morphism:
              From: Finite Field of size 1009
              To:   Finite Field of size 1009
              Defn: 11 |--> 11
            ]

            sage: k1.<a> = GF(1009^2, modulus='first_lexicographic')
            sage: k2.<b> = GF(1009^2, modulus='conway')
            sage: Hom(k1, k2).list()
            [
            Ring morphism:
              From: Finite Field in a of size 1009^2
              To:   Finite Field in b of size 1009^2
              Defn: a |--> 290*b + 864,
            Ring morphism:
              From: Finite Field in a of size 1009^2
              To:   Finite Field in b of size 1009^2
              Defn: a |--> 719*b + 145
            ]

        TESTS:

        Check that :issue:`11390` is fixed::

            sage: K = GF(1<<16,'a'); L = GF(1<<32,'b')
            sage: K.Hom(L)[0]
            Ring morphism:
              From: Finite Field in a of size 2^16
              To:   Finite Field in b of size 2^32
              Defn: a |--> b^29 + b^27 + b^26 + b^23 + b^21 + b^19 + b^18 + b^16 + b^14 + b^13 + b^11 + b^10 + b^9 + b^8 + b^7 + b^6 + b^5 + b^2 + b
        """
        try:
            return self.__list
        except AttributeError:
            pass
        D = self.domain()
        C = self.codomain()
        if D.characteristic() == C.characteristic() and Integer(D.degree()).divides(Integer(C.degree())):
            f = D.modulus()
            g = C['x'](f)
            r = g.roots()
            v = [self(D.hom(a, C)) for a, _ in r]
            v = Sequence(v, immutable=True, cr=True)
        else:
            v = Sequence([], immutable=True, cr=False)
        self.__list = v
        return v

    def __getitem__(self, n):
        """
        EXAMPLES::

            sage: H = Hom(GF(32, 'a'), GF(1024, 'b'))
            sage: H[1]
            Ring morphism:
              From: Finite Field in a of size 2^5
              To:   Finite Field in b of size 2^10
              Defn: a |--> b^7 + b^5
            sage: H[2:4]
            [
            Ring morphism:
              From: Finite Field in a of size 2^5
              To:   Finite Field in b of size 2^10
              Defn: a |--> b^8 + b^6 + b^2,
            Ring morphism:
              From: Finite Field in a of size 2^5
              To:   Finite Field in b of size 2^10
              Defn: a |--> b^9 + b^7 + b^6 + b^5 + b^4
            ]
        """
        return self.list()[n]

    def index(self, item):
        """
        Return the index of ``self``.

        EXAMPLES::

            sage: K.<z> = GF(1024)
            sage: g = End(K)[3]
            sage: End(K).index(g) == 3
            True
        """
        return self.list().index(item)

    def _an_element_(self):
        """
        Return an element of ``self``.

        TESTS::

            sage: Hom(GF(3^3, 'a'), GF(3^6, 'b')).an_element()  # random
            Ring morphism:
              From: Finite Field in a of size 3^3
              To:   Finite Field in b of size 3^6
              Defn: a |--> 2*b^5 + 2*b^4

            sage: Hom(GF(3^3, 'a'), GF(3^2, 'c')).an_element()
            Traceback (most recent call last):
            ...
            EmptySetError: no homomorphisms from Finite Field in a of size 3^3 to Finite Field in c of size 3^2

        .. TODO::

            Use a more sophisticated algorithm; see also :issue:`8751`.
        """
        K = self.domain()
        L = self.codomain()
        if K.degree() == 1:
            return L.coerce_map_from(K)
        elif not K.degree().divides(L.degree()):
            from sage.categories.sets_cat import EmptySetError
            raise EmptySetError('no homomorphisms from %s to %s' % (K, L))
        return K.hom([K.modulus().any_root(L)])


from sage.misc.persist import register_unpickle_override
register_unpickle_override('sage.rings.finite_field_morphism', 'FiniteFieldHomset', FiniteFieldHomset)
