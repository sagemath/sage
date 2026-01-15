r"""
Finite fields
"""
# ****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#                2025      Brian Heckel <heckelbri@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.enumerated_sets import EnumeratedSets
from sage.rings.integer import Integer
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

    def __contains__(self, x) -> bool:
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

        @cached_method
        def quadratic_nonresidue(self):
            r"""
            Return a random non square element of the finite field

            OUTPUT:
              A non-square element of the finite field; raises an error if
              the finite field is of even order.

            EXAMPLES::

                sage: k = GF((3, 10))
                sage: k.quadratic_nonresidue().is_square()
                False
                sage: k = GF((2, 10))
                sage: k in Fields()  # to let k be a finite field
                True
                sage: k.quadratic_nonresidue()
                Traceback (most recent call last):
                ...
                ValueError: there are no non-squares in finite fields of even order
            """
            # if the order is an even power of two
            # then every element is a square
            if self.characteristic() == 2:
                raise ValueError("there are no non-squares in finite fields of even order")
            for element in self:
                if not element.is_square():
                    return element

    class ElementMethods:
        def is_square(self) -> bool:
            r"""
            Test if the element is a square or has
            a square root element.

            OUTPUT:
              ``True`` if the element is a square ``False`` if not

            EXAMPLES::

                sage: S.<x> = GF(5)[]
                sage: f = S.irreducible_element(20)
                sage: k.<y> = S.quotient_ring(f)
                sage: k in Fields()
                True
                sage: k(2).is_square()
                True
                sage: k.quadratic_nonresidue().is_square()
                False
            """
            if self.is_zero():
                return True
            if self.parent().characteristic() == 2:
                return True
            q = self.parent().order()
            character = self**((q-1)//2)
            is_square = character == self.parent().one()
            return is_square

        def _tonelli(self):
            r"""
            Return a square root of the element if it exists
            using Tonelli's algorithm, only works for finite fields
            of odd characteristic.

            OUTPUT:
              A square root of the element; raises an error
              if the element is not a square

            EXAMPLES::

                sage: k.<a> = GF((5, 10))
                sage: k(2).is_square()
                True
                sage: k(2)._tonelli()^2 == k(2)
                True
                sage: k.quadratic_nonresidue()._tonelli()
                Traceback (most recent call last):
                ...
                ValueError: element is not a square
            """
            q = self.parent().cardinality()
            if not self.is_square():
                raise ValueError("element is not a square")
            g = self.parent().quadratic_nonresidue()
            even_exp, odd_order = (q - Integer(1)).val_unit(2)
            e = 0
            for i in range(2, even_exp+1):
                tmp = self * (pow(g, -e))

                condition = tmp**((q-1)//(2**i)) != self.parent().one()
                if condition:
                    e = 2**(i-1) + e
            h = self * (g**(-e))
            b = g**(e//2) * h**((odd_order+1)//2)
            return b

        def _cipolla(self):
            r"""
            Return a square root of the element if it exists
            using Cipolla's algorithm, more suited if order - 1
            is highly divisible by 2. Only works for finite fields
            of odd characteristic.

            OUTPUT:
              A square root of the element; raises an error
              if the element is not a square

            EXAMPLES::

                sage: k.<a> = GF((5, 10))
                sage: k(2).is_square()
                True
                sage: k(2)._cipolla()^2 == k(2)
                True
                sage: k.quadratic_nonresidue()._cipolla()
                Traceback (most recent call last):
                ...
                ValueError: element is not a square
            """
            parent = self.parent()
            q = parent.cardinality()
            if not self.is_square():
                raise ValueError("element is not a square")
            t = parent.random_element()
            root = t**2 - 4 * self
            while root.is_square():
                t = parent.random_element()
                root = t**2 - 4 * self
            from sage.rings.polynomial.polynomial_ring import polygen
            X = polygen(parent)
            f = X**2 - t*X + self
            b = pow(X, (q+1)//2, f)
            return b

        def sqrt(self, all: bool = False, algorithm: str = 'tonelli'):
            r"""
            Return the square root of the element if it exists.

            INPUT:

            - ``all`` -- boolean (default: ``False``); whether to return a list of
              all square roots or just a square root

            - ``algorithm`` -- string (default: 'tonelli'); the algorithm to use
              among ``'tonelli'``, ``'cipolla'``. Tonelli is typically faster but has
              a worse worst-case complexity than Cipolla. In particular, if the
              field cardinality minus 1 is highly divisible by 2 and has a large
              odd factor then Cipolla may perform better.

            OUTPUT:

            - if ``all=False``, a square root; raises an error if the element is not
              a square

            - if ``all=True``, a tuple of all distinct square roots. This tuple can have
              length 0, 1, or 2 depending on how many distinct square roots the
              element has.

            EXAMPLES::

                sage: S.<x> = GF(5)[]
                sage: f = S.irreducible_element(20)
                sage: k.<y> = S.quotient_ring(f)
                sage: k in Fields()
                True
                sage: k(2).is_square()
                True
                sage: k(2).sqrt()^2 == k(2)
                True
                sage: my_sqrts = k(4).sqrt(all=True)
                sage: len(k(4).sqrt(all=True))
                2
                sage: 2 in my_sqrts
                True
                sage: 3 in my_sqrts
                True
                sage: k.quadratic_nonresidue().sqrt()
                Traceback (most recent call last):
                ...
                ValueError: element is not a square
                sage: k.quadratic_nonresidue().sqrt(all=True)
                ()

            Here is an example where changing the algorithm results
            in a faster square root::

                sage: p = 141 * 2^141 + 1
                sage: S.<x> = GF(p)[]
                sage: f = S.irreducible_element(2)
                sage: k.<y> = S.quotient_ring(f)
                sage: k in Fields()
                True
                sage: k(2).sqrt(algorithm="cipolla")^2 == k(2)
                True

            ALGORITHM:

            The algorithms used come from chapter 7 of [BS1996]_.
            Let `q = p^n` be the order of the finite field, let `a` be the finite field element
            that we wish to find the square root of.

            - If `p = 2` then `a` is always a square, and the square root of `\sqrt{a} = a^{q / 2}`.
            - If `q \equiv 3 \pmod{4}` then if `a` is a square `\sqrt{a} = a^{\frac{q+1}{4}}`
            - For all other cases we use the algorithm given by the ``algorithm`` parameter.
            """
            cardinality = self.parent().order()
            if self.parent().characteristic() == 2:
                exponent = cardinality // 2
                square_root = self**exponent
                if all:
                    # we return a 1-tuple because the GF implementation does it
                    return (square_root,)
                else:
                    return square_root
            if not self.is_square():
                if all:
                    return ()
                else:
                    raise ValueError("element is not a square")
            if cardinality % 4 == 3:
                square_root = self**((cardinality+1)//4)
            elif algorithm == 'tonelli':
                square_root = self._tonelli()
            else:
                square_root = self._cipolla()
            if all:
                return (square_root, -square_root)
            return square_root

        square_root = sqrt
