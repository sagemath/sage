# sage_setup: distribution = sagemath-categories
r"""
Principal ideal domains
"""
# ****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.categories.category_singleton import Category_singleton
from sage.categories.unique_factorization_domains import UniqueFactorizationDomains


class PrincipalIdealDomains(Category_singleton):
    """
    The category of (constructive) principal ideal domains.

    By constructive, we mean that a single generator can be
    constructively found for any ideal given by a finite set of
    generators. Note that this constructive definition only implies
    that finitely generated ideals are principal. It is not clear what
    we would mean by an infinitely generated ideal.

    EXAMPLES::

      sage: PrincipalIdealDomains()
      Category of principal ideal domains
      sage: PrincipalIdealDomains().super_categories()
      [Category of unique factorization domains]

    See also :wikipedia:`Principal_ideal_domain`

    TESTS::

        sage: TestSuite(PrincipalIdealDomains()).run()
    """

    def super_categories(self):
        """
        EXAMPLES::

            sage: PrincipalIdealDomains().super_categories()
            [Category of unique factorization domains]
        """
        return [UniqueFactorizationDomains()]

    def additional_structure(self):
        """
        Return ``None``.

        Indeed, the category of principal ideal domains defines no
        additional structure: a ring morphism between two principal
        ideal domains is a principal ideal domain morphism.

        EXAMPLES::

            sage: PrincipalIdealDomains().additional_structure()
        """
        return None

    class ParentMethods:
        def _test_gcd_vs_xgcd(self, **options):
            r"""
            Check that gcd and xgcd are compatible if implemented.

            This test will prevent things like :issue:`17671` to happen again.

            TESTS::

                sage: ZZ._test_gcd_vs_xgcd()
                sage: QQ._test_gcd_vs_xgcd()
                sage: QQ['x']._test_gcd_vs_xgcd()
                sage: QQbar['x']._test_gcd_vs_xgcd()                                    # needs sage.rings.number_field
                sage: RR._test_gcd_vs_xgcd()
                sage: RR['x']._test_gcd_vs_xgcd()

            A slightly more involved example of polynomial ring with a non UFD
            base ring::

                sage: # needs sage.rings.number_field
                sage: K = QuadraticField(5)
                sage: O = K.maximal_order()
                sage: O in UniqueFactorizationDomains()
                False
                sage: R = PolynomialRing(O, 'x')
                sage: F = R.fraction_field()
                sage: F in PrincipalIdealDomains()
                True
                sage: F._test_gcd_vs_xgcd()
            """
            tester = self._tester(**options)
            elts = list(tester.some_elements())

            # there are some strange things in Sage doctests... so it is better
            # to cut the list in order to avoid lists of size 531441.
            elts = elts[:10]
            pairs = [(x, y) for x in elts for y in elts]

            try:
                xgcds = [x.xgcd(y) for x, y in pairs]
            except (AttributeError, NotImplementedError):
                return

            has_gcd = True
            try:
                gcds = [x.gcd(y) for x, y in pairs]
            except (AttributeError, NotImplementedError):
                has_gcd = False

            tester.assertTrue(has_gcd,
                    "The ring {} provides a xgcd but no gcd".format(self))
            for (x, y), gcd, xgcd in zip(pairs, gcds, xgcds):
                tester.assertTrue(gcd.parent() == self,
                        "The parent of the gcd is {} for element of {}".format(
                            gcd.parent(), self))
                tester.assertTrue(xgcd[0].parent() == self and
                                  xgcd[1].parent() == self == xgcd[2].parent(),
                                  "The parent of output in xgcd is different from "
                                  "the parent of input for elements in {}".format(self))
                tester.assertTrue(gcd == xgcd[0],
                        "The methods gcd and xgcd disagree on {}:\n"
                        "  gcd({},{}) = {}\n"
                        " xgcd({},{}) = {}\n".format(self, x, y, gcd, x, y, xgcd))

        def is_noetherian(self) -> bool:
            """
            Every principal ideal domain is Noetherian, so we return ``True``.

            EXAMPLES::

                sage: Zp(5).is_noetherian()                                                 # needs sage.rings.padics
                True
            """
            return True

        def class_group(self):
            """
            Return the trivial group, since the class group of a PID is trivial.

            EXAMPLES::

                sage: QQ.class_group()                                                      # needs sage.groups
                Trivial Abelian group
            """
            from sage.groups.abelian_gps.abelian_group import AbelianGroup
            return AbelianGroup([])

        def gcd(self, x, y, coerce=True):
            r"""
            Return the greatest common divisor of ``x`` and ``y``, as elements
            of ``self``.

            EXAMPLES:

            The integers are a principal ideal domain and hence a GCD domain::

                sage: ZZ.gcd(42, 48)
                6
                sage: 42.factor(); 48.factor()
                2 * 3 * 7
                2^4 * 3
                sage: ZZ.gcd(2^4*7^2*11, 2^3*11*13)
                88
                sage: 88.factor()
                2^3 * 11

            In a field, any nonzero element is a GCD of any nonempty set
            of nonzero elements. In previous versions, Sage used to return
            1 in the case of the rational field. However, since :issue:`10771`,
            the rational field is considered as the
            *fraction field* of the integer ring. For the fraction field
            of an integral domain that provides both GCD and LCM, it is
            possible to pick a GCD that is compatible with the GCD of the
            base ring::

                sage: QQ.gcd(ZZ(42), ZZ(48)); type(QQ.gcd(ZZ(42), ZZ(48)))
                6
                <class 'sage.rings.rational.Rational'>
                sage: QQ.gcd(1/2, 1/3)
                1/6

            Polynomial rings over fields are GCD domains as well. Here is a simple
            example over the ring of polynomials over the rationals as well as
            over an extension ring. Note that ``gcd`` requires x and y to be
            coercible::

                sage: # needs sage.rings.number_field
                sage: R.<x> = PolynomialRing(QQ)
                sage: S.<a> = NumberField(x^2 - 2, 'a')
                sage: f = (x - a)*(x + a); g = (x - a)*(x^2 - 2)
                sage: print(f); print(g)
                x^2 - 2
                x^3 - a*x^2 - 2*x + 2*a
                sage: f in R
                True
                sage: g in R
                False
                sage: R.gcd(f, g)
                Traceback (most recent call last):
                ...
                TypeError: Unable to coerce 2*a to a rational
                sage: R.base_extend(S).gcd(f,g)
                x^2 - 2
                sage: R.base_extend(S).gcd(f, (x - a)*(x^2 - 3))
                x - a
            """
            if coerce:
                x = self(x)
                y = self(y)
            return x.gcd(y)

        def content(self, x, y, coerce=True):
            r"""
            Return the content of `x` and `y`.

            This is the unique element `c` of
            ``self`` such that `x/c` and `y/c`
            are coprime and integral.

            EXAMPLES::

                sage: QQ.content(ZZ(42), ZZ(48)); type(QQ.content(ZZ(42), ZZ(48)))
                6
                <class 'sage.rings.rational.Rational'>
                sage: QQ.content(1/2, 1/3)
                1/6
                sage: factor(1/2); factor(1/3); factor(1/6)
                2^-1
                3^-1
                2^-1 * 3^-1
                sage: a = (2*3)/(7*11); b = (13*17)/(19*23)
                sage: factor(a); factor(b); factor(QQ.content(a,b))
                2 * 3 * 7^-1 * 11^-1
                13 * 17 * 19^-1 * 23^-1
                7^-1 * 11^-1 * 19^-1 * 23^-1

            Note the changes to the second entry::

                sage: c = (2*3)/(7*11); d = (13*17)/(7*19*23)
                sage: factor(c); factor(d); factor(QQ.content(c,d))
                2 * 3 * 7^-1 * 11^-1
                7^-1 * 13 * 17 * 19^-1 * 23^-1
                7^-1 * 11^-1 * 19^-1 * 23^-1
                sage: e = (2*3)/(7*11); f = (13*17)/(7^3*19*23)
                sage: factor(e); factor(f); factor(QQ.content(e,f))
                2 * 3 * 7^-1 * 11^-1
                7^-3 * 13 * 17 * 19^-1 * 23^-1
                7^-3 * 11^-1 * 19^-1 * 23^-1
            """
            if coerce:
                x = self(x)
                y = self(y)
            return x.content(y)

        def _ideal_class_(self, n=0):
            """
            Ideals in PIDs have their own special class.

            EXAMPLES::

                sage: ZZ._ideal_class_()
                <class 'sage.rings.ideal.Ideal_pid'>
            """
            from sage.rings.ideal import Ideal_pid
            return Ideal_pid

    class ElementMethods:
        pass
