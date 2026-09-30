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

from functools import lru_cache

from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.enumerated_sets import EnumeratedSets
from sage.misc.cachefunc import cached_method
from sage.rings.integer import Integer


def _pickleable_quotient_field_morphism(parent, mapping, codomain):
    r"""
    Rebuild ``mapping`` as a pickleable morphism defined by generator images.

    The isomorphism returned by polynomial quotient fields can contain local
    callables.  Rebuilding it from generator images also handles towers of
    quotient fields and makes the resulting map safe to store in extension
    parents::

        sage: from sage.categories.finite_fields import _pickleable_quotient_field_morphism
        sage: R.<x> = GF(5)[]
        sage: K.<a> = R.quotient(x**2 + 2)
        sage: _, mapping, F = K._isomorphic_ring()
        sage: morphism = _pickleable_quotient_field_morphism(K, mapping, F)
        sage: loads(dumps(morphism))(a) == morphism(a)
        True
    """
    from sage.rings.polynomial.polynomial_quotient_ring import (
        PolynomialQuotientRing_generic,
    )

    base = parent.base_ring()
    if isinstance(base, PolynomialQuotientRing_generic):
        base_map = _pickleable_quotient_field_morphism(base, mapping,
                                                       codomain)
    else:
        base_map = base.hom([mapping(gen) for gen in base.gens()], codomain,
                            check=False)
    return parent.hom([mapping(parent.gen())], codomain, base_map=base_map,
                      check=False)


@lru_cache(maxsize=256)
def _sqrt_extension(parent, name):
    r"""
    Return a reusable quadratic finite-field extension of ``parent``.

    The returned parent is a genuine finite field whenever Sage can construct
    an embedding of ``parent`` into its absolute representation.  This makes
    roots of different elements compatible with each other::

        sage: from sage.categories.finite_fields import _sqrt_extension
        sage: K.<a> = GF(7**3)
        sage: E, embedding = _sqrt_extension(K, 'sqrt_ext')
        sage: E in FiniteFields() and E.has_coerce_map_from(K)
        True
        sage: embedding(K.gen()).parent() is E
        True

    Polynomial quotient fields are connected to an isomorphic absolute
    finite field before the quadratic extension is formed::

        sage: R.<x> = GF(5)[]
        sage: L.<b> = R.quotient(x**2 + 2)
        sage: E, embedding = _sqrt_extension(L, 'sqrt_ext')
        sage: E in FiniteFields() and embedding(L.gen()).parent() is E
        True
    """
    from sage.rings.finite_rings.finite_field_base import (
        FiniteField as FiniteField_base,
    )

    if isinstance(parent, FiniteField_base):
        return parent.extension(2, name, map=True)

    from sage.rings.polynomial.polynomial_quotient_ring import (
        PolynomialQuotientRing_generic,
    )
    if (isinstance(parent, PolynomialQuotientRing_generic)
            and parent in FiniteFields()):
        try:
            _, to_field, field = parent._isomorphic_ring()
        except NotImplementedError:
            return None
        extension, embedding = field.extension(2, name, map=True)
        mapping = _pickleable_quotient_field_morphism(
            parent, embedding * to_field, extension
        )
        if not extension.has_coerce_map_from(parent):
            extension.register_coercion(mapping)
        return extension, extension.coerce_map_from(parent)

    if parent in FiniteFields():
        from sage.rings.finite_rings.finite_field_constructor import GF
        extension = GF(parent.order()**2, name)
        embedding = extension.coerce_map_from(parent)
        if embedding is not None:
            return extension, embedding

    return None


@lru_cache(maxsize=256)
def _sqrt_extension_nonresidue(parent, name):
    r"""
    Return a fixed nonsquare and its square root in a quadratic extension.

    Reusing this root reduces later square roots of nonsquares in ``parent``
    to a square root in ``parent`` and one multiplication in the extension::

        sage: from sage.categories.finite_fields import (
        ....:     _sqrt_extension, _sqrt_extension_nonresidue)
        sage: K = GF(1009)
        sage: nonsquare, radical = _sqrt_extension_nonresidue(K, 'sqrt_ext')
        sage: extension, embedding = _sqrt_extension(K, 'sqrt_ext')
        sage: radical**2 == embedding(nonsquare)
        True

    Return ``None`` when this reduction is not profitable.  In particular,
    Givaro computes the root directly faster, while polynomial quotient
    fields may have a comparatively expensive square-root implementation in
    the source parent.
    """
    extension_data = _sqrt_extension(parent, name)
    if extension_data is None:
        return None
    extension, embedding = extension_data
    from sage.rings.finite_rings.finite_field_base import (
        FiniteField as FiniteField_base,
    )
    from sage.rings.finite_rings.finite_field_givaro import FiniteField_givaro
    if (not isinstance(parent, FiniteField_base)
            or isinstance(extension, FiniteField_givaro)
            or parent.characteristic() == 2):
        return None
    nonsquare = parent.quadratic_nonresidue()
    radical = embedding(nonsquare).sqrt(extend=False)
    return nonsquare, radical


def _sqrt_in_extension(element, *, all_roots, name, algorithm=None):
    r"""
    Return square roots of ``element`` in a quadratic extension.

    Standard finite fields use a common genuine finite-field parent, so roots
    of different nonsquares can be combined::

        sage: K.<a> = GF(7**3)
        sage: from sage.categories.finite_fields import _sqrt_in_extension
        sage: r = _sqrt_in_extension(K(3), all_roots=False, name=None)
        sage: s = _sqrt_in_extension(K(5), all_roots=False, name=None)
        sage: r.parent() is s.parent() and (r + s).parent() is r.parent()
        True
        sage: from sage.rings.finite_rings.finite_field_base import FiniteField
        sage: isinstance(r.parent(), FiniteField)
        True
        sage: r.parent().variable_names()
        ('sqrt_ext',)
        sage: r.parent().multiplicative_generator().parent() is r.parent()
        True
        sage: r**2 == K(3) and s**2 == K(5)
        True

    For finite non-domains, repeated fallback calls reuse the polynomial
    quotient parent::

        sage: R.<x> = Zmod(4)[]
        sage: A.<a> = R.quotient(x**2)
        sage: r = _sqrt_in_extension(a, all_roots=False, name='w')
        sage: s = _sqrt_in_extension(a, all_roots=False, name='w')
        sage: r.parent() is s.parent() and r**2 == a
        True
    """
    from sage.rings.polynomial.polynomial_quotient_ring import PolynomialQuotientRing
    from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

    parent = element.parent()
    if name is None:
        name = 'sqrt_ext'
    extension_data = _sqrt_extension(parent, name)
    if extension_data is not None:
        extension, embedding = extension_data
        nonresidue_data = _sqrt_extension_nonresidue(parent, name)
        if nonresidue_data is not None and not element.is_square():
            nonsquare, radical = nonresidue_data
            # The quotient of two nonsquares in a finite field is a square.
            if element == nonsquare:
                square_root = radical
            else:
                ratio_root = (element / nonsquare).sqrt(
                    extend=False, algorithm=algorithm
                )
                square_root = embedding(ratio_root) * radical
            if all_roots:
                return [square_root, -square_root]
            return square_root
        return embedding(element).sqrt(extend=False, all=all_roots,
                                       algorithm=algorithm)

    polynomial_ring = PolynomialRing(parent, 'x')
    x = polynomial_ring.gen()
    extension = PolynomialQuotientRing(
        polynomial_ring, x**2 - polynomial_ring(element), names=name
    )
    square_root = extension.gen()
    if all_roots:
        if parent.characteristic() == 2:
            return [square_root]
        return [square_root, -square_root]
    return square_root


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

        def _tonelli(self, raise_on_failure=True):
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
                sage: k(2)._tonelli()**2 == k(2)
                True
                sage: k.quadratic_nonresidue()._tonelli()
                Traceback (most recent call last):
                ...
                ValueError: element is not a square
                sage: k.quadratic_nonresidue()._tonelli(
                ....:     raise_on_failure=False) is None
                True
            """
            parent = self.parent()
            if parent.characteristic() == 2:
                raise ValueError("Tonelli's algorithm requires odd characteristic")
            if self.is_zero():
                return self
            q = parent.cardinality()
            even_exp, odd_order = (q - 1).val_unit(2)
            one = parent.one()
            residue = self**odd_order
            character = residue
            for _ in range(even_exp - 1):
                character *= character
            if character != one:
                if raise_on_failure:
                    raise ValueError("element is not a square")
                return None
            generator = parent.quadratic_nonresidue()
            correction = generator**odd_order
            square_root = self**((odd_order + 1) // 2)
            exponent = even_exp

            while residue != one:
                i = 1
                power = residue * residue
                while i < exponent and power != one:
                    power *= power
                    i += 1
                if i == exponent:
                    raise ArithmeticError("Tonelli's algorithm failed")
                factor = correction**(2**(exponent - i - 1))
                square_root *= factor
                factor *= factor
                residue *= factor
                correction = factor
                exponent = i

            return square_root

        def _cipolla(self, check=True):
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
                sage: k(2)._cipolla()**2 == k(2)
                True
                sage: k.quadratic_nonresidue()._cipolla()
                Traceback (most recent call last):
                ...
                ValueError: element is not a square
                sage: k.quadratic_nonresidue()._cipolla(check=False)
                Traceback (most recent call last):
                ...
                ValueError: element is not a square
            """
            parent = self.parent()
            if parent.characteristic() == 2:
                raise ValueError("Cipolla's algorithm requires odd characteristic")
            if self.is_zero():
                return self
            q = parent.cardinality()
            if check and not self.is_square():
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
            square_root = b[0]
            if square_root * square_root != self:
                raise ValueError("element is not a square")
            return square_root

        def sqrt(self, *, extend=False, all=False, algorithm=None, name=None):
            r"""
            Return the square root of the element if it exists.

            INPUT:

            - ``extend`` -- boolean (default: ``False``); if ``True``, return
              roots in a quadratic extension when necessary

            - ``all`` -- boolean (default: ``False``); whether to return all
              square roots or just one

            - ``algorithm`` -- optional algorithm hint (default: ``None``).
              ``'cipolla'`` selects Cipolla's algorithm; ``'tonelli'``,
              ``None``, and unsupported hints select the backend default,
              Tonelli's algorithm. Tonelli is typically faster but has a worse
              worst-case complexity than Cipolla. In particular, if the field
              cardinality minus 1 is highly divisible by 2 and has a large odd
              factor then Cipolla may perform better.

            - ``name`` -- string (default: ``None``); name of the generator when
              a quadratic extension is created

            OUTPUT:

            - if ``all=False``, a square root in the parent or, when
              ``extend=True``, in a quadratic extension; raises an error if no
              root exists and extension is disabled

            - if ``all=True``, a list of all distinct square roots in the
              selected parent.  This list can have length 0, 1, or 2 depending
              on how many distinct square roots the element has.

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
                []

            TESTS:

            The common finite-field keyword interface is accepted, and roots
            returned by Cipolla's algorithm belong to the original field::

                sage: for method in ((y**2).sqrt, (y**2).square_root):
                ....:     r = method(extend=False, algorithm='cipolla')
                ....:     assert r.parent() is k and r**2 == y**2
                sage: for method in (k(0).sqrt, k(0).square_root):
                ....:     for algorithm in ('tonelli', 'cipolla'):
                ....:         assert method(algorithm=algorithm) == 0
                ....:         assert method(all=True,
                ....:                       algorithm=algorithm) == [k(0)]
                sage: for method in (k(1).sqrt, k(1).square_root):
                ....:     assert method(algorithm='backend-default')**2 == 1

            A nonsquare can be lifted to a quadratic extension::

                sage: a = k.quadratic_nonresidue()
                sage: for method in (a.sqrt, a.square_root):
                ....:     s = method(extend=True, name='s')
                ....:     assert s**2 == a and s.parent() in FiniteFields()

            Both method names implement the same keyword contract::

                sage: q = k(4)
                sage: for method in (q.sqrt, q.square_root):
                ....:     for algorithm in (None, 'tonelli', 'cipolla'):
                ....:         roots = method(extend=False, all=True,
                ....:                        algorithm=algorithm, name='s')
                ....:         assert isinstance(roots, list)
                ....:         assert len(roots) == 2
                ....:         assert all(r.parent() is k and r**2 == q
                ....:                    for r in roots)

            The contract is uniform across the concrete finite-field
            implementations::

                sage: from inspect import signature
                sage: fields = [GF(7),
                ....:           GF(next_prime(2^40)),
                ....:           GF(9, 'g', implementation='givaro'),
                ....:           GF(2^8, 'n', implementation='ntl'),
                ....:           GF(3^3, 'p', implementation='pari_ffelt')]
                sage: R.<u> = GF(5)[]
                sage: fields.append(R.quotient(u^2 + 2, 'q'))
                sage: signatures = {str(signature(method))
                ....:               for field in fields
                ....:               for method in (field.one().sqrt,
                ....:                              field.one().square_root)}
                sage: signatures
                {'(*, extend=False, all=False, algorithm=None, name=None)'}
                sage: for field in fields:
                ....:     value = field.gen()**2
                ....:     for method in (value.sqrt, value.square_root):
                ....:         roots = method(all=True,
                ....:                        algorithm='backend-default')
                ....:         assert isinstance(roots, list)
                ....:         assert roots and all(root**2 == value
                ....:                              for root in roots)

            Optional arguments are keyword-only, so old backend-specific
            positional orders cannot be confused with the common interface::

                sage: for method in (q.sqrt, q.square_root):
                ....:     try:
                ....:         method(True)
                ....:     except TypeError:
                ....:         pass
                ....:     else:
                ....:         raise AssertionError("optional argument was positional")

            The category implementation also handles characteristic two and
            the exponent shortcut for fields of order congruent to three
            modulo four::

                sage: R2.<z> = GF(2)[]
                sage: K2.<b> = R2.quotient(z^3 + z + 1)
                sage: b._cipolla()
                Traceback (most recent call last):
                ...
                ValueError: Cipolla's algorithm requires odd characteristic
                sage: for method in (b.sqrt, b.square_root):
                ....:     roots = method(extend=True, all=True,
                ....:                    algorithm='cipolla', name='w')
                ....:     assert isinstance(roots, list) and len(roots) == 1
                ....:     assert roots[0].parent() is K2 and roots[0]**2 == b
                sage: R3.<z> = GF(3)[]
                sage: K3.<b> = R3.quotient(z^3 - z + 1)
                sage: q3 = b**2
                sage: for method in (q3.sqrt, q3.square_root):
                ....:     roots = method(extend=False, all=True,
                ....:                    algorithm='cipolla')
                ....:     assert isinstance(roots, list) and len(roots) == 2
                ....:     assert all(r.parent() is K3 and r**2 == q3
                ....:                for r in roots)

            EXAMPLES:

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
            if self.is_zero():
                if all:
                    return [self]
                return self
            if self.parent().characteristic() == 2:
                exponent = cardinality // 2
                square_root = self**exponent
                if all:
                    return [square_root]
                return square_root
            is_square = True
            if cardinality % 4 == 3:
                square_root = self**((cardinality+1)//4)
                is_square = square_root * square_root == self
            elif algorithm == 'cipolla':
                is_square = self.is_square()
                if is_square:
                    square_root = self._cipolla(check=False)
            else:
                square_root = self._tonelli(raise_on_failure=False)
                is_square = square_root is not None
            if not is_square:
                if extend:
                    return _sqrt_in_extension(
                        self, all_roots=all, name=name, algorithm=algorithm
                    )
                if all:
                    return []
                raise ValueError("element is not a square")
            if all:
                return [square_root, -square_root]
            return square_root

        square_root = sqrt
