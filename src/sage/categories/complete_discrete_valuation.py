r"""
Complete Discrete Valuation Rings (CDVR) and Fields (CDVF)
"""
from __future__ import absolute_import
#**************************************************************************
#  Copyright (C) 2013 Xavier Caruso <xavier.caruso@normalesup.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#**************************************************************************


from sage.misc.abstract_method import abstract_method

from sage.categories.category_singleton import Category_singleton
from .discrete_valuation import DiscreteValuationRings, DiscreteValuationFields
#from sage.misc.cachefunc import cached_method

class CompleteDiscreteValuationRings(Category_singleton):
    """
    The category of complete discrete valuation rings

    EXAMPLES::

        sage: Zp(7) in CompleteDiscreteValuationRings()
        True
        sage: QQ in CompleteDiscreteValuationRings()
        False
        sage: QQ[['u']] in CompleteDiscreteValuationRings()
        True
        sage: Qp(7) in CompleteDiscreteValuationRings()
        False
        sage: TestSuite(CompleteDiscreteValuationRings()).run()
    """
    def super_categories(self):
        """
        EXAMPLES::

            sage: CompleteDiscreteValuationRings().super_categories()
            [Category of discrete valuation rings]
        """
        return [DiscreteValuationRings()]

    class ElementMethods:
        @abstract_method
        def precision_absolute(self):
            """
            Return the absolute precision of this element.

            EXAMPLES::

                sage: R = Zp(7)
                sage: x = R(7); x
                7 + O(7^21)
                sage: x.precision_absolute()
                21
            """

        @abstract_method
        def precision_relative(self):
            """
            Return the relative precision of this element.

            EXAMPLES::

                sage: R = Zp(7)
                sage: x = R(7); x
                7 + O(7^21)
                sage: x.precision_relative()
                20
            """

class CompleteDiscreteValuationFields(Category_singleton):
    """
    The category of complete discrete valuation fields

    EXAMPLES::

        sage: Zp(7) in CompleteDiscreteValuationFields()
        False
        sage: QQ in CompleteDiscreteValuationFields()
        False
        sage: LaurentSeriesRing(QQ,'u') in CompleteDiscreteValuationFields()
        True
        sage: Qp(7) in CompleteDiscreteValuationFields()
        True
        sage: TestSuite(CompleteDiscreteValuationFields()).run()
    """

    def super_categories(self):
        """
        EXAMPLES::

            sage: CompleteDiscreteValuationFields().super_categories()
            [Category of discrete valuation fields]
        """
        return [DiscreteValuationFields()]

    class ElementMethods:
        @abstract_method
        def precision_absolute(self):
            """
            Return the absolute precision of this element.

            EXAMPLES::

                sage: K = Qp(7)
                sage: x = K(7); x
                7 + O(7^21)
                sage: x.precision_absolute()
                21
            """

        @abstract_method
        def precision_relative(self):
            """
            Return the relative precision of this element.

            EXAMPLES::

                sage: K = Qp(7)
                sage: x = K(7); x
                7 + O(7^21)
                sage: x.precision_relative()
                20
            """
