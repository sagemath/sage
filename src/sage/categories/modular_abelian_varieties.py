r"""
Modular abelian varieties
"""
#*****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from .category_types import Category_over_base
from .category_with_axiom import CategoryWithAxiom
from .homsets import HomsetsCategory
from .rings import Rings
from .sets_cat import Sets

class ModularAbelianVarieties(Category_over_base):
    """
    The category of modular abelian varieties over a given field.

    EXAMPLES::

        sage: ModularAbelianVarieties(QQ)
        Category of modular abelian varieties over Rational Field
    """
    def __init__(self, Y):
        """
        TESTS::

            sage: C = ModularAbelianVarieties(QQ)
            sage: C
            Category of modular abelian varieties over Rational Field
            sage: TestSuite(C).run()

            sage: ModularAbelianVarieties(ZZ)
            Traceback (most recent call last):
            ...
              assert Y.is_field()
            AssertionError
        """
        assert Y.is_field()
        Category_over_base.__init__(self, Y)

    def base_field(self):
        """
        EXAMPLES::

            sage: ModularAbelianVarieties(QQ).base_field()
            Rational Field
        """
        return self.base()

    def super_categories(self):
        """
        EXAMPLES::

            sage: ModularAbelianVarieties(QQ).super_categories()
            [Category of sets]
        """
        return [Sets()] # FIXME

    class Homsets(HomsetsCategory):

        class Endset(CategoryWithAxiom):
            def extra_super_categories(self):
                """
                Implement the fact that an endset of modular abelian variety is a ring.

                EXAMPLES::

                    sage: ModularAbelianVarieties(QQ).Endsets().extra_super_categories()
                    [Category of rings]
                """
                return [Rings()]
