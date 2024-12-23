# sage_setup: distribution = sagemath-categories
r"""
Integral domains

TEST:

A few tests for the method ``is_integrally_closed``::

    sage: ZZ.is_integrally_closed()
    True
    sage: QQ.is_integrally_closed()
    True
    sage: QQbar.is_integrally_closed()                                          # needs sage.rings.number_field
    True
    sage: GF(5).is_integrally_closed()
    True
    sage: Z5 = Integers(5); Z5
    Ring of integers modulo 5
    sage: Z5.is_integrally_closed()
    Traceback (most recent call last):
    ...
    NotImplementedError

Note that this raises a :exc:`NotImplementedError` if the answer is not known.
"""
# ****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2012 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.misc.lazy_import import lazy_import
from sage.misc.lazy_attribute import lazy_class_attribute
from sage.categories.category_with_axiom import CategoryWithAxiom
from sage.categories.category_singleton import Category_contains_method_by_parent_class
from sage.categories.domains import Domains
lazy_import('sage.categories.fields', 'Fields')


class IntegralDomains(CategoryWithAxiom):
    """
    The category of integral domains.

    An integral domain is commutative ring with no zero divisors, or
    equivalently a commutative domain.

    EXAMPLES::

        sage: C = IntegralDomains(); C
        Category of integral domains
        sage: sorted(C.super_categories(), key=str)
        [Category of commutative rings, Category of domains]
        sage: C is Domains().Commutative()
        True
        sage: C is Rings().Commutative().NoZeroDivisors()
        True

    TESTS::

        sage: TestSuite(C).run()
    """
    _base_category_class_and_axiom = (Domains, "Commutative")

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: GF(4, "a") in IntegralDomains()                                       # needs sage.rings.finite_rings
            True
            sage: QQ in IntegralDomains()
            True
            sage: ZZ in IntegralDomains()
            True
            sage: IntegerModRing(4) in IntegralDomains()
            False
            sage: IntegerModRing(5) in IntegralDomains()
            True

        This implementation will not be needed anymore once every
        field in Sage will be properly declared in the category
        :class:`IntegralDomains() <IntegralDomains>`.
        """
        try:
            return self._contains_helper(x) or x.is_integral_domain()
        except Exception:
            return False

    @lazy_class_attribute
    def _contains_helper(cls):
        """
        Helper for containment tests in the category of integral domains.

        This helper just tests whether the given object's category
        is already known to be a sub-category of the category of
        integral domains. There are, however, rings that are initialised
        as plain commutative rings and found out to be integral domains
        only afterwards. Hence, this helper alone is not enough
        for a proper containment test.

        TESTS::

            sage: R = Zmod(7)
            sage: R.category()
            Join of Category of finite commutative rings
             and Category of subquotients of monoids
             and Category of quotients of semigroups
             and Category of finite enumerated sets
            sage: ID = IntegralDomains()
            sage: ID._contains_helper(R)
            False
            sage: R in ID  # This changes the category!
            True
            sage: ID._contains_helper(R)
            True
        """
        return Category_contains_method_by_parent_class(cls())

    class ParentMethods:
        def is_integral_domain(self, proof=True):
            r"""
            Return ``True``, since this in an object of the category
            of integral domains.

            EXAMPLES::

                sage: ZZ.is_integral_domain()
                True
                sage: QQ.is_integral_domain()
                True
                sage: Parent(QQ, category=IntegralDomains()).is_integral_domain()
                True

                sage: L.<z> = LazyLaurentSeriesRing(QQ)                                 # needs sage.combinat
                sage: L.is_integral_domain()                                            # needs sage.combinat
                True
                sage: L.is_integral_domain(proof=True)                                  # needs sage.combinat
                True

                sage: ZZ['x'].is_integral_domain()
                True
            """
            return True

        def is_field(self, proof=True):
            r"""
            Return ``True`` if this ring is a field.

            EXAMPLES::

                sage: ZZ['x'].is_field()
                False
            """
            if self.is_finite():
                return True
            if proof:
                raise NotImplementedError(f"unable to determine whether or not {self} is a field.")
            return False

        def localization(self, additional_units, names=None, normalize=True, category=None):
            """
            Return the localization of ``self`` at the given additional units.

            EXAMPLES::

                sage: R.<x, y> = GF(3)[]
                sage: R.localization((x*y, x**2 + y**2))                                    # needs sage.rings.finite_rings
                Multivariate Polynomial Ring in x, y over Finite Field of size 3
                 localized at (y, x, x^2 + y^2)
                sage: ~y in _                                                               # needs sage.rings.finite_rings
                True
            """
            from sage.rings.localization import Localization
            return Localization(self, additional_units, names=names, normalize=normalize, category=category)

        def _test_fraction_field(self, **options):
            r"""
            Test that the fraction field, if it is implemented, works
            correctly.

            EXAMPLES::

                sage: ZZ._test_fraction_field()
            """
            tester = self._tester(**options)
            try:
                fraction_field = self.fraction_field()
            except (AttributeError, ImportError):
                # some integral domains do not implement fraction_field() yet
                if self in Fields():
                    raise
                return

            for x in tester.some_elements():
                # check that we can coerce into the fraction field
                fraction_field.coerce(x)
                # and convert back from it
                z = self(x)
                tester.assertEqual(x, z)

    class ElementMethods:
        pass
