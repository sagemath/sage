r"""
Number fields
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

from sage.categories.category_singleton import Category_singleton
from sage.categories.basic import Fields

class NumberFields(Category_singleton):
    r"""
    The category of number fields.

    EXAMPLES:

    We create the category of number fields::

        sage: C = NumberFields()
        sage: C
        Category of number fields

    By definition, it is infinite::

        sage: NumberFields().Infinite() is NumberFields()
        True

    Notice that the rational numbers `\QQ` *are* considered as
    an object in this category::

        sage: RationalField() in C
        True

    However, we can define a degree 1 extension of `\QQ`, which is of
    course also in this category::

        sage: x = PolynomialRing(RationalField(), 'x').gen()
        sage: K = NumberField(x - 1, 'a'); K
        Number Field in a with defining polynomial x - 1
        sage: K in C
        True

    Number fields all lie in this category, regardless of the name
    of the variable::

        sage: K = NumberField(x^2 + 1, 'a')
        sage: K in C
        True

    TESTS::

        sage: TestSuite(NumberFields()).run()
    """

    def super_categories(self):
        """
        EXAMPLES::

            sage: NumberFields().super_categories()
            [Category of infinite fields]
        """
        return [Fields().Infinite()]

    def __contains__(self, x):
        r"""
        Returns True if ``x`` is a number field.

        EXAMPLES::

            sage: NumberField(x^2+1,'a') in NumberFields()
            True
            sage: QuadraticField(-97,'theta') in NumberFields()
            True
            sage: CyclotomicField(97) in NumberFields()
            True

        Note that the rational numbers QQ are a number field::

            sage: QQ in NumberFields()
            True
            sage: ZZ in NumberFields()
            False
        """
        import sage.rings.number_field.number_field_base
        return sage.rings.number_field.number_field_base.is_NumberField(x)

    def _call_(self, x):
        r"""
        Constructs an object in this category from the data in ``x``,
        or throws a TypeError.

        EXAMPLES::

            sage: C = NumberFields()

            sage: C(QQ)
            Rational Field

            sage: C(NumberField(x^2+1,'a'))
            Number Field in a with defining polynomial x^2 + 1

            sage: C(UnitGroup(NumberField(x^2+1,'a')))  # indirect doctest
            Number Field in a with defining polynomial x^2 + 1

            sage: C(ZZ)
            Traceback (most recent call last):
            ...
            TypeError: unable to canonically associate a number field to Integer Ring
        """
        try:
            return x.number_field()
        except AttributeError:
            raise  TypeError("unable to canonically associate a number field to %s"%x)


    class ParentMethods:
        pass

    class ElementMethods:
        pass
