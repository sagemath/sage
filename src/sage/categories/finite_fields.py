# sage_setup: distribution = sagemath-categories
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
from sage.rings.integer import Integer


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

            sage: GF(4, "a") in FiniteFields()                                          # needs sage.rings.finite_rings
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

            sage: FiniteFields()(GF(4, "a"))                                            # needs sage.rings.finite_rings
            Finite Field in a of size 2^2
            sage: FiniteFields()(RationalField())   # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: unable to canonically associate a finite field to Rational Field
        """
        raise TypeError("unable to canonically associate a finite field to %s" % x)
        # TODO: local dvr ring?

    class ParentMethods:
        def is_perfect(self):
            r"""
            Return whether this field is perfect, i.e., every element has a `p`-th
            root. Always returns ``True`` since finite fields are perfect.

            EXAMPLES::

                sage: GF(2).is_perfect()
                True
            """
            return True

        def zeta_order(self):
            """
            Return the order of the distinguished root of unity in ``self``.

            EXAMPLES::

                sage: GF(9,'a').zeta_order()
                8
                sage: GF(9,'a').zeta()
                a
                sage: GF(9,'a').zeta().multiplicative_order()
                8
            """
            return self.order() - 1

        def zeta(self, n=None):
            """
            Return an element of multiplicative order ``n`` in this finite
            field. If there is no such element, raise :exc:`ValueError`.

            .. WARNING::

                In general, this returns an arbitrary element of the correct
                order. There are no compatibility guarantees:
                ``F.zeta(9)^3`` may not be equal to ``F.zeta(3)``.

            EXAMPLES::

                sage: k = GF(7)
                sage: k.zeta()
                3
                sage: k.zeta().multiplicative_order()
                6
                sage: k.zeta(3)
                2
                sage: k.zeta(3).multiplicative_order()
                3
                sage: k = GF(49, 'a')
                sage: k.zeta().multiplicative_order()
                48
                sage: k.zeta(6)
                3
                sage: k.zeta(5)
                Traceback (most recent call last):
                ...
                ValueError: no 5th root of unity in Finite Field in a of size 7^2

            Even more examples::

                sage: GF(9,'a').zeta_order()
                8
                sage: GF(9,'a').zeta()
                a
                sage: GF(9,'a').zeta(4)
                a + 1
                sage: GF(9,'a').zeta()^2
                a + 1

            This works even in very large finite fields, provided that ``n``
            can be factored (see :issue:`25203`)::

                sage: k.<a> = GF(2^2000)
                sage: p = 8877945148742945001146041439025147034098690503591013177336356694416517527310181938001
                sage: z = k.zeta(p)
                sage: z
                a^1999 + a^1996 + a^1995 + a^1994 + ... + a^7 + a^5 + a^4 + 1
                sage: z ^ p
                1
            """
            if n is None:
                return self.multiplicative_generator()

            from sage.rings.integer import Integer
            n = Integer(n)
            grouporder = self.order() - 1
            co_order = grouporder // n
            if co_order * n != grouporder:
                raise ValueError("no {}th root of unity in {}".format(n, self))

            # If the co_order is small or we know a multiplicative
            # generator, use a multiplicative generator
            mg = self.multiplicative_generator
            if mg.cache is not None or co_order <= 500000:
                return mg() ** co_order
            return self._element_of_factored_order(n.factor())

        def _element_of_factored_order(self, F):
            """
            Return an element of ``self`` of order ``n`` where ``n`` is
            given in factored form.

            This is copied from the cython implementation in
            ``finite_field_base.pyx`` which is kept as it may be faster.

            INPUT:

            - ``F`` -- the factorization of the required order. The order
              must be a divisor of ``self.order() - 1`` but this is not
              checked.

            EXAMPLES::

                sage: k = Zmod(1913)
                sage: k in Fields()  # to let k be a finite field
                True
                sage: k._element_of_factored_order(factor(1912))
                3
            """
            n = Integer(1)
            primes = []
            for p, e in F:
                primes.append(p)
                n *= p**e

            N = self.order() - 1
            c = N // n

            # We check whether (x + g)^c has the required order, where
            # x runs through the finite field.
            # This has the advantage that g is the first element we try,
            # so if that was a chosen to be a multiplicative generator,
            # we are done immediately. Second, the PARI finite field
            # iterator gives all the constant elements first, so we try
            # (g+(constant))^c before anything else.
            g = self.gen()
            if g == self.one():
                # this allows to handle the ring Integers(prime)
                g = self.multiplicative_generator()
            for x in self:
                a = (g + x)**c
                if not a:
                    continue
                if all(a**(n // p) != 1 for p in primes):
                    return a
            raise AssertionError("no element found")

    class ElementMethods:
        pass
