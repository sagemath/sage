r"""
Commutative rings
"""
# ****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2013 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.cartesian_product import CartesianProductsCategory


class CommutativeRings(CategoryWithAxiom):
    """
    The category of commutative rings

    commutative rings with unity, i.e. rings with commutative * and
    a multiplicative identity

    EXAMPLES::

         sage: C = CommutativeRings(); C
         Category of commutative rings
         sage: C.super_categories()
         [Category of rings, Category of commutative monoids]

    TESTS::

        sage: TestSuite(C).run()

        sage: QQ['x,y,z'] in CommutativeRings()
        True
        sage: GroupAlgebra(DihedralGroup(3), QQ) in CommutativeRings()
        False
        sage: MatrixSpace(QQ,2,2) in CommutativeRings()
        False

    GroupAlgebra should be fixed::

        sage: GroupAlgebra(CyclicPermutationGroup(3), QQ) in CommutativeRings() # todo: not implemented
        True

    """
    class ParentMethods:
        def _test_divides(self, **options):
            r"""
            Run generic tests on the method :meth:`divides`.

            EXAMPLES::

                sage: ZZ._test_divides()
            """
            tester = self._tester(**options)

            # 1. is there a divides method ?
            a = self.an_element()
            try:
                a.divides
            except AttributeError:
                return

            # 2. divisibility of 0 and 1
            z = self.zero()
            o = self.one()

            tester.assertTrue(z.divides(z))
            tester.assertTrue(o.divides(o))
            tester.assertTrue(o.divides(z))
            tester.assertIs(z.divides(o), self.is_zero())

            if not self.is_exact():
                return

            # 3. divisibility of some elements
            for a, b in tester.some_elements(repeat=2):
                try:
                    test = a.divides(a * b)
                except NotImplementedError:
                    pass
                else:
                    tester.assertTrue(test)

        def over(self, base=None, gen=None, gens=None, name=None, names=None):
            r"""
            Return this ring, considered as an extension of ``base``.

            INPUT:

            - ``base`` -- a commutative ring or a morphism or ``None``
              (default: ``None``); the base of this extension or its defining
              morphism

            - ``gen`` -- a generator of this extension (over its base) or ``None``
              (default: ``None``);

            - ``gens`` -- a list of generators of this extension (over its base)
              or ``None`` (default: ``None``);

            - ``name`` -- a variable name or ``None`` (default: ``None``)

            - ``names`` -- a list or a tuple of variable names or ``None``
              (default: ``None``)

            EXAMPLES:

            We construct an extension of finite fields::

                sage: F = GF(5^2)
                sage: k = GF(5^4)
                sage: z4 = k.gen()

                sage: K = k.over(F)
                sage: K
                Field in z4 with defining polynomial x^2 + (4*z2 + 3)*x + z2 over its base

            If not explicitly given, the default generator of the top ring
            (here k) is used and the same name is kept::

                sage: K.gen()
                z4
                sage: K(z4)
                z4

            However, it is possible to specify another generator and/or
            another name. For example::

                sage: Ka = k.over(F, name='a')
                sage: Ka
                Field in a with defining polynomial x^2 + (4*z2 + 3)*x + z2 over its base
                sage: Ka.gen()
                a
                sage: Ka(z4)
                a

                sage: Kb = k.over(F, gen=-z4+1, name='b')
                sage: Kb
                Field in b with defining polynomial x^2 + z2*x + 4 over its base
                sage: Kb.gen()
                b
                sage: Kb(-z4+1)
                b

            Note that the shortcut ``K.<a>`` is also available::

                sage: KKa.<a> = k.over(F)
                sage: KKa is Ka
                True

            Building an extension on top of another extension is allowed::

                sage: L = GF(5^12).over(K)
                sage: L
                Field in z12 with defining polynomial x^3 + (1 + (4*z2 + 2)*z4)*x^2 + (2 + 2*z4)*x - z4 over its base
                sage: L.base_ring()
                Field in z4 with defining polynomial x^2 + (4*z2 + 3)*x + z2 over its base

            The successive bases of an extension are accessible via the
            method :meth:`sage.rings.ring_extension.RingExtension_generic.bases`::

                sage: L.bases()
                [Field in z12 with defining polynomial x^3 + (1 + (4*z2 + 2)*z4)*x^2 + (2 + 2*z4)*x - z4 over its base,
                 Field in z4 with defining polynomial x^2 + (4*z2 + 3)*x + z2 over its base,
                 Finite Field in z2 of size 5^2]

            When ``base`` is omitted, the canonical base of the ring is used::

                sage: S.<x> = QQ[]
                sage: E = S.over()
                sage: E
                Univariate Polynomial Ring in x over Rational Field over its base
                sage: E.base_ring()
                Rational Field

            Here is an example where ``base`` is a defining morphism::

                sage: k.<a> = QQ.extension(x^2 - 2)
                sage: l.<b> = QQ.extension(x^4 - 2)
                sage: f = k.hom([b^2])
                sage: L = l.over(f)
                sage: L
                Field in b with defining polynomial x^2 - a over its base
                sage: L.base_ring()
                Number Field in a with defining polynomial x^2 - 2

            Similarly, one can create a tower of extensions::

                sage: K = k.over()
                sage: L = l.over(Hom(K,l)(f))
                sage: L
                Field in b with defining polynomial x^2 - a over its base
                sage: L.base_ring()
                Field in a with defining polynomial x^2 - 2 over its base
                sage: L.bases()
                [Field in b with defining polynomial x^2 - a over its base,
                 Field in a with defining polynomial x^2 - 2 over its base,
                 Rational Field]
            """
            from sage.rings.ring_extension import RingExtension
            if name is not None:
                if names is not None:
                    raise ValueError("keyword argument 'name' cannot be combined with 'names'")
                names = (name,)
            if gen is not None:
                if gens is not None:
                    raise ValueError("keyword argument 'gen' cannot be combined with 'gens'")
                gens = (gen,)
            return RingExtension(self, base, gens, names)

    class ElementMethods:
        pass

    class Finite(CategoryWithAxiom):
        r"""
        Check that Sage knows that Cartesian products of finite commutative
        rings is a finite commutative ring.

        EXAMPLES::

            sage: cartesian_product([Zmod(34), GF(5)]) in Rings().Commutative().Finite()
            True
        """
        class ParentMethods:
            def cyclotomic_cosets(self, q, cosets=None):
                r"""
                Return the (multiplicative) orbits of ``q`` in the ring.

                Let `R` be a finite commutative ring. The group of invertible
                elements `R^*` in `R` gives rise to a group action on `R` by
                multiplication.  An orbit of the subgroup generated by an
                invertible element `q` is called a `q`-*cyclotomic coset* (since
                in a finite ring, each invertible element is a root of unity).

                These cosets arise in the theory of minimal polynomials of
                finite fields, duadic codes and combinatorial designs. Fix a
                primitive element `z` of `GF(q^k)`. The minimal polynomial of
                `z^s` over `GF(q)` is given by

                .. MATH::

                         M_s(x) = \prod_{i \in C_s} (x - z^i),


                where `C_s` is the `q`-cyclotomic coset mod `n` containing `s`,
                `n = q^k - 1`.

                .. NOTE::

                    When `R = \ZZ / n \ZZ` the smallest element of each coset is
                    sometimes called a *coset leader*. This function returns
                    sorted lists so that the coset leader will always be the
                    first element of the coset.

                INPUT:

                - ``q`` -- an invertible element of the ring

                - ``cosets`` -- an optional lists of elements of ``self``. If
                  provided, the function only return the list of cosets that
                  contain some element from ``cosets``.

                OUTPUT:

                A list of lists.

                EXAMPLES::

                    sage: Zmod(11).cyclotomic_cosets(2)
                    [[0], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]]
                    sage: Zmod(15).cyclotomic_cosets(2)
                    [[0], [1, 2, 4, 8], [3, 6, 9, 12], [5, 10], [7, 11, 13, 14]]

                Since the group of invertible elements of a finite field is
                cyclic, the set of squares is a particular case of cyclotomic
                coset::

                    sage: K = GF(25,'z')
                    sage: a = K.multiplicative_generator()
                    sage: K.cyclotomic_cosets(a**2,cosets=[1])
                    [[1, 2, 3, 4, z + 1, z + 3,
                      2*z + 1, 2*z + 2, 3*z + 3,
                      3*z + 4, 4*z + 2, 4*z + 4]]
                    sage: sorted(b for b in K if not b.is_zero() and b.is_square())
                    [1, 2, 3, 4, z + 1, z + 3,
                     2*z + 1, 2*z + 2, 3*z + 3,
                     3*z + 4, 4*z + 2, 4*z + 4]

                We compute some examples of minimal polynomials::

                    sage: K = GF(27,'z')
                    sage: a = K.multiplicative_generator()
                    sage: R.<X> = PolynomialRing(K, 'X')
                    sage: a.minimal_polynomial('X')
                    X^3 + 2*X + 1
                    sage: cyc3 = Zmod(26).cyclotomic_cosets(3,cosets=[1]); cyc3
                    [[1, 3, 9]]
                    sage: prod(X - a**i for i in cyc3[0])
                    X^3 + 2*X + 1

                    sage: (a**7).minimal_polynomial('X')
                    X^3 + X^2 + 2*X + 1
                    sage: cyc7 = Zmod(26).cyclotomic_cosets(3,cosets=[7]); cyc7
                    [[7, 11, 21]]
                    sage: prod(X - a**i for i in cyc7[0])
                    X^3 + X^2 + 2*X + 1

                Cyclotomic cosets of fields are useful in combinatorial design
                theory to provide so called difference families (see
                :wikipedia:`Difference_set` and
                :mod:`~sage.combinat.designs.difference_family`). This is
                illustrated on the following examples::

                    sage: K = GF(5)
                    sage: a = K.multiplicative_generator()
                    sage: H = K.cyclotomic_cosets(a**2, cosets=[1,2]); H
                    [[1, 4], [2, 3]]
                    sage: sorted(x-y for D in H for x in D for y in D if x != y)
                    [1, 2, 3, 4]

                    sage: K = GF(37)
                    sage: a = K.multiplicative_generator()
                    sage: H = K.cyclotomic_cosets(a**4, cosets=[1]); H
                    [[1, 7, 9, 10, 12, 16, 26, 33, 34]]
                    sage: sorted(x-y for D in H for x in D for y in D if x != y)
                    [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, ..., 33, 34, 34, 35, 35, 36, 36]

                The method ``cyclotomic_cosets`` works on any finite commutative
                ring::

                    sage: R = cartesian_product([GF(7), Zmod(14)])
                    sage: a = R((3,5))
                    sage: R.cyclotomic_cosets((3,5), [(1,1)])
                    [[(1, 1), (2, 11), (3, 5), (4, 9), (5, 3), (6, 13)]]
                """
                q = self(q)

                try:
                    ~q
                except ZeroDivisionError:
                    raise ValueError("%s is not invertible in %s"%(q,self))

                if cosets is None:
                    rest = set(self)
                else:
                    rest = set(self(x) for x in cosets)

                orbits = []
                while rest:
                    x0 = rest.pop()
                    o = [x0]
                    x = q*x0
                    while x != x0:
                        o.append(x)
                        rest.discard(x)
                        x *= q
                    o.sort()
                    orbits.append(o)

                orbits.sort()
                return orbits

    class CartesianProducts(CartesianProductsCategory):
        def extra_super_categories(self):
            r"""
            Let Sage knows that Cartesian products of commutative rings is a
            commutative ring.

            EXAMPLES::

                sage: CommutativeRings().Commutative().CartesianProducts().extra_super_categories()
                [Category of commutative rings]
                sage: cartesian_product([ZZ, Zmod(34), QQ, GF(5)]) in CommutativeRings()
                True
            """
            return [CommutativeRings()]
