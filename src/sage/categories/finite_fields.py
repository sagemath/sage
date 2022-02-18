r"""
Finite fields
"""
# ****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.enumerated_sets import EnumeratedSets
from sage.misc.cachefunc import cached_method

class FiniteFields(CategoryWithAxiom):
    """
    The category of finite fields.

    EXAMPLES::

        sage: K = FiniteFields(); K
        Category of finite enumerated fields

    A finite field is a finite monoid with the structure of a field;
    it is currently assumed to be enumerated::

        sage: K.super_categories()
        [Category of fields,
         Category of finite commutative rings,
         Category of finite enumerated sets]

    Some examples of membership testing and coercion::

        sage: FiniteField(17) in K
        True
        sage: RationalField() in K
        False
        sage: K(RationalField())
        Traceback (most recent call last):
        ...
        TypeError: unable to canonically associate a finite field to Rational Field

    TESTS::

        sage: K is Fields().Finite()
        True
        sage: TestSuite(K).run()
    """

    def extra_super_categories(self):
        r"""
        Any finite field is assumed to be endowed with an enumeration.

        TESTS::

            sage: Fields().Finite().extra_super_categories()
            [Category of finite enumerated sets]
            sage: FiniteFields().is_subcategory(FiniteEnumeratedSets())
            True
        """
        return [EnumeratedSets().Finite()]

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: GF(4, "a") in FiniteFields()
            True
            sage: QQ in FiniteFields()
            False
            sage: IntegerModRing(4) in FiniteFields()
            False
        """
        from sage.categories.fields import Fields
        return x in Fields() and x.is_finite()

    # As is, this does no more than the usual __call__ of Category, but for the error message
    def _call_(self, x):
        """
        EXAMPLES::

            sage: FiniteFields()(GF(4, "a"))
            Finite Field in a of size 2^2
            sage: FiniteFields()(RationalField())   # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: unable to canonically associate a finite field to Rational Field
        """
        raise TypeError("unable to canonically associate a finite field to %s"%x)
        # TODO: local dvr ring?

    class ParentMethods:
        pass

    class ElementMethods:
        pass

    class MorphismMethods:
        @cached_method
        def _section(self):
            """
            Return the ``inverse`` of this embedding.

            It is a partially defined map whose domain is the codomain
            of the embedding, but which is only defined on the image of
            the embedding.

            EXAMPLES::

                sage: from sage.rings.finite_rings.hom_finite_field import FiniteFieldHomomorphism_generic
                sage: k.<t> = GF(3^7)
                sage: K.<T> = GF(3^21)
                sage: f = FiniteFieldHomomorphism_generic(Hom(k, K))
                sage: g = f.section(); g
                Section of Ring morphism:
                  From: Finite Field in t of size 3^7
                  To:   Finite Field in T of size 3^21
                  Defn: t |--> T^20 + 2*T^18 + T^16 + 2*T^13 + T^9 + 2*T^8 + T^7 + T^6 + T^5 + T^3 + 2*T^2 + T
                sage: g(f(t^3+t^2+1))
                t^3 + t^2 + 1
                sage: g(T)
                Traceback (most recent call last):
                ...
                ValueError: T is not in the image of Ring morphism:
                  From: Finite Field in t of size 3^7
                  To:   Finite Field in T of size 3^21
                  Defn: t |--> T^20 + 2*T^18 + T^16 + 2*T^13 + T^9 + 2*T^8 + T^7 + T^6 + T^5 + T^3 + 2*T^2 + T
            """
            if self.domain().is_prime_field():
                from sage.rings.finite_rings.hom_prime_finite_field import SectionFiniteFieldHomomorphism_prime as section_class
            else:
                from sage.rings.finite_rings.hom_finite_field import SectionFiniteFieldHomomorphism_generic as section_class
            return section_class(self)

        def _inverse_image_element(self, b):
            r"""
            Return the unique ``a`` such that ``self(a) = b`` if one such exists.

            This method is simply a shorthand for calling the map returned by
            ``self.section()`` on ``b``.

            EXAMPLES::

                sage: k.<t> = GF(3^7)
                sage: set_random_seed(0)
                sage: K.<T>, f = k.extension(3, absolute=True, map=True)
                sage: t.minpoly()(f(t))
                0
                sage: b = f(t^2); b
                2*T^19 + T^16 + T^14 + T^12 + 2*T^11 + T^10 + T^9 + 2*T^8 + T^7 + 2
                sage: f.inverse_image(b)
                t^2
                sage: f.inverse_image(T)
                Traceback (most recent call last):
                ...
                ValueError: T is not in the image of ...
            """
            return self.section()(b)

        def inverse_image(self, I):
            """
            Return the inverse image of an ideal or an element in the codomain.

            INPUT:

            - ``I`` -- an ideal or element in the codomain

            OUTPUT:

            For an ideal `I` in the codomain, this returns the largest ideal in the
            domain whose image is contained in `I`.

            Given an element `b` in the codomain, this returns the element
            `a` in the domain such that ``self(a) = b`` if one such exists.
            """
            from sage.categories.ring_ideals import RingIdeals
            B = self.codomain()
            if I in RingIdeals(B):
                if I.is_zero():
                    return self.domain().ideal(0)
                else:
                    return self.domain().ideal(1)
            elif I in B:
                return self.section()(I)
            else:
                raise ValueError("not an ideal or element in codomain %s" % B)
