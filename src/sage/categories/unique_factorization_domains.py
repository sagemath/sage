# sage_setup: distribution = sagemath-categories
r"""
Unique factorization domains
"""
# ****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.misc.lazy_attribute import lazy_class_attribute
from sage.misc.misc_c import prod
from sage.categories.category_singleton import Category_singleton
from sage.categories.category_singleton import Category_contains_method_by_parent_class
from sage.categories.gcd_domains import GcdDomains


class UniqueFactorizationDomains(Category_singleton):
    """
    The category of (constructive) unique factorization domains.

    In a constructive unique factorization domain we can
    constructively factor members into a product of a finite number
    of irreducible elements.

    EXAMPLES::

        sage: UniqueFactorizationDomains()
        Category of unique factorization domains
        sage: UniqueFactorizationDomains().super_categories()
        [Category of gcd domains]

    TESTS::

        sage: TestSuite(UniqueFactorizationDomains()).run()
    """

    def super_categories(self):
        """
        EXAMPLES::

            sage: UniqueFactorizationDomains().super_categories()
            [Category of gcd domains]
        """
        return [GcdDomains()]

    def additional_structure(self):
        """
        Return whether ``self`` is a structure category.

        .. SEEALSO:: :meth:`Category.additional_structure`

        The category of unique factorization domains does not define
        additional structure: a ring morphism between unique factorization
        domains is a unique factorization domain morphism.

        EXAMPLES::

            sage: UniqueFactorizationDomains().additional_structure()
        """
        return None

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: GF(4, "a") in UniqueFactorizationDomains()                            # needs sage.rings.finite_rings
            True
            sage: QQ in UniqueFactorizationDomains()
            True
            sage: ZZ in UniqueFactorizationDomains()
            True
            sage: IntegerModRing(4) in UniqueFactorizationDomains()
            False
            sage: IntegerModRing(5) in UniqueFactorizationDomains()
            True

        This implementation will not be needed anymore once every
        field in Sage will be properly declared in the category
        :class:`UniqueFactorizationDomains() <UniqueFactorizationDomains>`.
        """
        try:
            return self._contains_helper(x) or x.is_unique_factorization_domain()
        except Exception:
            return False

    @lazy_class_attribute
    def _contains_helper(cls):
        """
        Helper for containment tests in the category of unique
        factorization domains.

        This helper just tests whether the given object's category
        is already known to be a sub-category of the category of
        unique factorization domains. There are, however, rings that
        are initialised as plain commutative rings and found out to be
        unique factorization domains only afterwards. Hence, this helper
        alone is not enough for a proper containment test.

        TESTS::

            sage: R = Zmod(7)
            sage: R.category()
            Join of Category of finite commutative rings
                and Category of subquotients of monoids
                and Category of quotients of semigroups
                and Category of finite enumerated sets
            sage: ID = UniqueFactorizationDomains()
            sage: ID._contains_helper(R)
            False
            sage: R in ID  # This changes the category!
            True
            sage: ID._contains_helper(R)
            True
        """
        return Category_contains_method_by_parent_class(cls())

    class ParentMethods:
        def is_unique_factorization_domain(self, proof=True):
            """
            Return True, since this in an object of the category of unique factorization domains.

            EXAMPLES::

                sage: UFD = UniqueFactorizationDomains()
                sage: Parent(QQ, category=UFD).is_unique_factorization_domain()
                True
            """
            return True

        def _gcd_univariate_polynomial(self, f, g):
            """
            Return the greatest common divisor of ``f`` and ``g``.

            INPUT:

            - ``f``, ``g`` -- two polynomials defined over this UFD

            .. NOTE::

                This is a helper method for
                :meth:`sage.rings.polynomial.polynomial_element.Polynomial.gcd`.

            ALGORITHM:

            Algorithm 3.3.1 in [Coh1993]_, based on pseudo-division.

            EXAMPLES::

                sage: R.<x> = PolynomialRing(ZZ, sparse=True)
                sage: S.<T> = R[]
                sage: p = (-3*x^2 - x)*T^3 - 3*x*T^2 + (x^2 - x)*T + 2*x^2 + 3*x - 2
                sage: q = (-x^2 - 4*x - 5)*T^2 + (6*x^2 + x + 1)*T + 2*x^2 - x
                sage: quo,rem=p.pseudo_quo_rem(q); quo,rem
                ((3*x^4 + 13*x^3 + 19*x^2 + 5*x)*T + 18*x^4 + 12*x^3 + 16*x^2 + 16*x,
                 (-113*x^6 - 106*x^5 - 133*x^4 - 101*x^3 - 42*x^2 - 41*x)*T
                         - 34*x^6 + 13*x^5 + 54*x^4 + 126*x^3 + 134*x^2 - 5*x - 50)
                sage: (-x^2 - 4*x - 5)^(3-2+1) * p == quo*q + rem
                True

            Check that :issue:`23620` has been resolved::

                sage: # needs sage.rings.padics
                sage: R.<x> = ZpFM(2)[]
                sage: f = 2*x + 2
                sage: g = 4*x + 2
                sage: f.gcd(g).parent() is R
                True
            """
            if f.degree() < g.degree():
                A,B = g, f
            else:
                A,B = f, g

            if B.is_zero():
                return A

            a = b = self.zero()
            for c in A.coefficients():
                a = a.gcd(c)
                if a.is_one():
                    break
            for c in B.coefficients():
                b = b.gcd(c)
                if b.is_one():
                    break

            d = a.gcd(b)
            A = A // a
            B = B // b
            g = h = 1

            delta = A.degree()-B.degree()
            _,R = A.pseudo_quo_rem(B)

            while R.degree() > 0:
                A = B
                B = R // (g*h**delta)
                g = A.leading_coefficient()
                h = h*g**delta // h**delta
                delta = A.degree() - B.degree()
                _, R = A.pseudo_quo_rem(B)

            if R.is_zero():
                b = self.zero()
                for c in B.coefficients():
                    b = b.gcd(c)
                    if b.is_one():
                        break

                return d*B // b

            return f.parent()(d)

    class ElementMethods:
        # prime?
        # squareFree
        # factor

        def radical(self, *args, **kwds):
            r"""
            Return the radical of this element, i.e. the product of its
            irreducible factors.

            This default implementation calls ``squarefree_decomposition`` if
            available, and ``factor`` otherwise.

            .. SEEALSO:: :meth:`squarefree_part`

            EXAMPLES::

                sage: Pol.<x> = QQ[]
                sage: (x^2*(x-1)^3).radical()
                x^2 - x
                sage: pol = 37 * (x-1)^3 * (x-2)^2 * (x-1/3)^7 * (x-3/7)
                sage: pol.radical()
                37*x^4 - 2923/21*x^3 + 1147/7*x^2 - 1517/21*x + 74/7

                sage: Integer(10).radical()
                10
                sage: Integer(-100).radical()
                10
                sage: Integer(0).radical()
                Traceback (most recent call last):
                ...
                ArithmeticError: radical of 0 is not defined

            The next example shows how to compute the radical of a number,
            assuming no prime > 100000 has exponent > 1 in the factorization::

                sage: n = 2^1000-1; n / radical(n, limit=100000)
                125

            TESTS::

                sage: radical(pol)
                37*x^4 - 2923/21*x^3 + 1147/7*x^2 - 1517/21*x + 74/7

                sage: Integer(20).radical()
                10
            """
            if self.is_zero():
                raise ArithmeticError("radical of 0 is not defined")
            try:
                decomp = self.squarefree_decomposition()
            except AttributeError:
                return self.factor(*args, **kwds).radical_value()
            else:
                return prod(fac for fac, mult in decomp)

        def squarefree_part(self):
            r"""
            Return the square-free part of this element, i.e. the product
            of its irreducible factors appearing with odd multiplicity.

            This default implementation calls ``squarefree_decomposition``.

            .. SEEALSO:: :meth:`radical`

            EXAMPLES::

                sage: Pol.<x> = QQ[]
                sage: (x^2*(x-1)^3).squarefree_part()
                x - 1
                sage: pol = 37 * (x-1)^3 * (x-2)^2 * (x-1/3)^7 * (x-3/7)
                sage: pol.squarefree_part()
                37*x^3 - 1369/21*x^2 + 703/21*x - 37/7

            TESTS::

                sage: squarefree_part(pol)
                37*x^3 - 1369/21*x^2 + 703/21*x - 37/7
            """
            decomp = self.squarefree_decomposition()
            return prod(fac for fac, mult in decomp if mult % 2 == 1)
