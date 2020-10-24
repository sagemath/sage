r"""
Function fields
"""
#*****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008      Teresa Gomez-Diaz (CNRS) <Teresa.Gomez-Diaz@univ-mlv.fr>
#                2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category import Category
from sage.misc.cachefunc import cached_method
from sage.categories.basic import Fields

class FunctionFields(Category):
    r"""
    The category of function fields.

    EXAMPLES:

    We create the category of function fields::

        sage: C = FunctionFields()
        sage: C
        Category of function fields

    TESTS::

        sage: TestSuite(FunctionFields()).run()
    """
    @cached_method
    def super_categories(self):
        """
        Returns the Category of which this is a direct sub-Category
        For a list off all super categories see all_super_categories

        EXAMPLES::

            sage: FunctionFields().super_categories()
            [Category of fields]
        """
        return[Fields()]

    def _call_(self, x):
        r"""
        Constructs an object in this category from the data in ``x``,
        or throws a TypeError.

        EXAMPLES::

            sage: C = FunctionFields()
            sage: K.<x>=FunctionField(QQ)
            sage: C(K)
            Rational function field in x over Rational Field
            sage: Ky.<y> = K[]
            sage: L = K.extension(y^2-x)
            sage: C(L)
            Function field in y defined by y^2 - x
            sage: C(L.equation_order())
            Function field in y defined by y^2 - x
        """
        try:
            return x.function_field()
        except AttributeError:
            raise  TypeError("unable to canonically associate a function field to %s"%x)

    class ParentMethods:
        pass

    class ElementMethods:
        pass
