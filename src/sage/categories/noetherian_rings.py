# sage_setup: distribution = sagemath-categories
r"""
Noetherian rings

EXAMPLES::

    sage: from sage.categories.noetherian_rings import NoetherianRings
    sage: GF(4, "a") in NoetherianRings()                                       # needs sage.rings.finite_rings
    True
    sage: QQ in NoetherianRings()
    True
    sage: ZZ in NoetherianRings()
    True
    sage: IntegerModRing(4) in NoetherianRings()
    True
    sage: IntegerModRing(5) in NoetherianRings()
    True
"""
# ****************************************************************************
#  Copyright (C) 2008 Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2012 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.categories.category import Category
from sage.categories.commutative_rings import CommutativeRings


class NoetherianRings(Category):
    """
    The category of Noetherian rings.

    A Noetherian ring is a commutative ring in which
    every ideal is finitely generated.

    See :wikipedia:`Noetherian ring`

    EXAMPLES::

        sage: from sage.categories.noetherian_rings import NoetherianRings
        sage: C = NoetherianRings(); C
        Category of noetherian rings
        sage: sorted(C.super_categories(), key=str)
        [Category of commutative rings]

    TESTS::

        sage: TestSuite(C).run()
    """
    def super_categories(self):
        """
        EXAMPLES::

            sage: from sage.categories.noetherian_rings import NoetherianRings
            sage: NoetherianRings().super_categories()
            [Category of commutative rings]
        """
        return [CommutativeRings()]

    class ParentMethods:
        def is_noetherian(self, proof=True):
            r"""
            Return ``True``, since this in an object of the category
            of Noetherian rings.

            EXAMPLES::

                sage: ZZ.is_noetherian()
                True
                sage: QQ.is_noetherian()
                True
                sage: ZZ['x'].is_noetherian()
                True
                sage: R.<x> = PolynomialRing(QQ)
                sage: R.is_noetherian()
                True

                sage: L.<z> = LazyLaurentSeriesRing(QQ)                                 # needs sage.combinat
                sage: L.is_noetherian()                                            # needs sage.combinat
                True
            """
            return True

    class ElementMethods:
        pass
