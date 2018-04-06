r"""
Principal ideal domains
"""
#*****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_singleton import Category_singleton
from sage.categories.unique_factorization_domains import UniqueFactorizationDomains

class PrincipalIdealDomains(Category_singleton):
    """
    The category of (constructive) principal ideal domains

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

            This test will prevent things like :trac:`17671` to happen again.

            TESTS::

                sage: ZZ._test_gcd_vs_xgcd()
                sage: QQ._test_gcd_vs_xgcd()
                sage: QQ['x']._test_gcd_vs_xgcd()
                sage: QQbar['x']._test_gcd_vs_xgcd()
                sage: RR._test_gcd_vs_xgcd()
                sage: RR['x']._test_gcd_vs_xgcd()

            A slightly more involved example of polynomial ring with a non UFD
            base ring::

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
            pairs = [(x,y) for x in elts for y in elts]

            try:
                xgcds = [x.xgcd(y) for x,y in pairs]
            except (AttributeError,NotImplementedError):
                return

            has_gcd = True
            try:
                gcds = [x.gcd(y) for x,y in pairs]
            except (AttributeError,NotImplementedError):
                has_gcd = False

            tester.assertTrue(has_gcd,
                    "The ring {} provides a xgcd but no gcd".format(self))
            for (x,y),gcd,xgcd in zip(pairs,gcds,xgcds):
                tester.assertTrue(gcd.parent()==self,
                        "The parent of the gcd is {} for element of {}".format(
                            gcd.parent(), self))
                tester.assertTrue(xgcd[0].parent()==self and
                        xgcd[1].parent()==self and xgcd[2].parent()==self,
                        "The parent of output in xgcd is different from "
                        "the parent of input for elements in {}".format(self))
                tester.assertTrue(gcd==xgcd[0],
                        "The methods gcd and xgcd disagree on {}:\n"
                        "  gcd({},{}) = {}\n"
                        " xgcd({},{}) = {}\n".format(self,x,y,gcd,x,y,xgcd))

    class ElementMethods:
        pass
