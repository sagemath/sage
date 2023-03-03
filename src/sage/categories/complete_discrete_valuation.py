r"""
Complete Discrete Valuation Rings (CDVR) and Fields (CDVF)
"""
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

        sage: Zp(7) in CompleteDiscreteValuationRings()                                                                 # optional - sage.rings.padics
        True
        sage: QQ in CompleteDiscreteValuationRings()
        False
        sage: QQ[['u']] in CompleteDiscreteValuationRings()
        True
        sage: Qp(7) in CompleteDiscreteValuationRings()                                                                 # optional - sage.rings.padics
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
        def valuation(self):
            """
            Return the valuation of this element.

            EXAMPLES::

                sage: R = Zp(7)                                                                                         # optional - sage.rings.padics
                sage: x = R(7); x                                                                                       # optional - sage.rings.padics
                7 + O(7^21)
                sage: x.valuation()                                                                                     # optional - sage.rings.padics
                1
            """

        def denominator(self):
            """
            Return the denominator of this element normalized
            as a power of the uniformizer

            EXAMPLES::

                sage: K = Qp(7)                                                                                         # optional - sage.rings.padics
                sage: x = K(1/21)                                                                                       # optional - sage.rings.padics
                sage: x.denominator()                                                                                   # optional - sage.rings.padics
                7 + O(7^21)

                sage: x = K(7)                                                                                          # optional - sage.rings.padics
                sage: x.denominator()                                                                                   # optional - sage.rings.padics
                1 + O(7^20)

            Note that the denominator lives in the ring of integers::

                sage: x.denominator().parent()                                                                          # optional - sage.rings.padics
                7-adic Ring with capped relative precision 20

            When the denominator is indistinguishable from 0 and the
            precision on the input is `O(p^n)`, the return value is `1`
            if `n` is nonnegative and `p^(-n)` otherwise::

                sage: x = K(0, 5); x                                                                                    # optional - sage.rings.padics
                O(7^5)
                sage: x.denominator()                                                                                   # optional - sage.rings.padics
                1 + O(7^20)

                sage: x = K(0, -5); x                                                                                   # optional - sage.rings.padics
                O(7^-5)
                sage: x.denominator()                                                                                   # optional - sage.rings.padics
                7^5 + O(7^25)
            """
            return self.parent()(1)

        def numerator(self):
            """
            Return the numerator of this element, normalized in such a
            way that `x = x.numerator() / x.denominator()` always holds
            true.

            EXAMPLES::

                sage: K = Qp(7, 5)                                                                                      # optional - sage.rings.padics
                sage: x = K(1/21)                                                                                       # optional - sage.rings.padics
                sage: x.numerator()                                                                                     # optional - sage.rings.padics
                5 + 4*7 + 4*7^2 + 4*7^3 + 4*7^4 + O(7^5)

                sage: x == x.numerator() / x.denominator()                                                              # optional - sage.rings.padics
                True

            Note that the numerator lives in the ring of integers::

                sage: x.numerator().parent()                                                                            # optional - sage.rings.padics
                7-adic Ring with capped relative precision 5

            TESTS::

                sage: x = K(0, -5); x                                                                                   # optional - sage.rings.padics
                O(7^-5)
                sage: x.numerator()                                                                                     # optional - sage.rings.padics
                O(7^0)
                sage: x.denominator()                                                                                   # optional - sage.rings.padics
                7^5 + O(7^10)
            """
            return self

        @abstract_method
        def lift_to_precision(self, absprec=None):
            """
            Return another element of the same parent with absolute precision
            at least ``absprec``, congruent to this element modulo the
            precision of this element.

            INPUT:

            - ``absprec`` -- an integer or ``None`` (default: ``None``), the
              absolute precision of the result. If ``None``, lifts to the maximum
              precision allowed.

            .. NOTE::

                If setting ``absprec`` that high would violate the precision cap,
                raises a precision error.  Note that the new digits will not
                necessarily be zero.

            EXAMPLES::

                sage: R = ZpCA(17)                                                                                      # optional - sage.rings.padics
                sage: R(-1, 2).lift_to_precision(10)                                                                    # optional - sage.rings.padics
                16 + 16*17 + O(17^10)
                sage: R(1, 15).lift_to_precision(10)                                                                    # optional - sage.rings.padics
                1 + O(17^15)
                sage: R(1, 15).lift_to_precision(30)                                                                    # optional - sage.rings.padics
                Traceback (most recent call last):
                ...
                PrecisionError: precision higher than allowed by the precision cap
                sage: R(-1, 2).lift_to_precision().precision_absolute() == R.precision_cap()                            # optional - sage.rings.padics
                True

                sage: R = Zp(5); c = R(17, 3); c.lift_to_precision(8)                                                   # optional - sage.rings.padics
                2 + 3*5 + O(5^8)
                sage: c.lift_to_precision().precision_relative() == R.precision_cap()                                   # optional - sage.rings.padics
                True

            """

class CompleteDiscreteValuationFields(Category_singleton):
    """
    The category of complete discrete valuation fields

    EXAMPLES::

        sage: Zp(7) in CompleteDiscreteValuationFields()                                                                # optional - sage.rings.padics
        False
        sage: QQ in CompleteDiscreteValuationFields()
        False
        sage: LaurentSeriesRing(QQ, 'u') in CompleteDiscreteValuationFields()
        True
        sage: Qp(7) in CompleteDiscreteValuationFields()                                                                # optional - sage.rings.padics
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
        def valuation(self):
            """
            Return the valuation of this element.

            EXAMPLES::

                sage: K = Qp(7)                                                                                         # optional - sage.rings.padics
                sage: x = K(7); x                                                                                       # optional - sage.rings.padics
                7 + O(7^21)
                sage: x.valuation()                                                                                     # optional - sage.rings.padics
                1
            """

        def denominator(self):
            """
            Return the denominator of this element normalized
            as a power of the uniformizer

            EXAMPLES::

                sage: K = Qp(7)                                                                                         # optional - sage.rings.padics
                sage: x = K(1/21)                                                                                       # optional - sage.rings.padics
                sage: x.denominator()                                                                                   # optional - sage.rings.padics
                7 + O(7^21)

                sage: x = K(7)                                                                                          # optional - sage.rings.padics
                sage: x.denominator()                                                                                   # optional - sage.rings.padics
                1 + O(7^20)

            Note that the denominator lives in the ring of integers::

                sage: x.denominator().parent()                                                                          # optional - sage.rings.padics
                7-adic Ring with capped relative precision 20

            When the denominator is indistinguishable from 0 and the
            precision on the input is `O(p^n)`, the return value is `1`
            if `n` is nonnegative and `p^(-n)` otherwise::

                sage: x = K(0, 5); x                                                                                    # optional - sage.rings.padics
                O(7^5)
                sage: x.denominator()                                                                                   # optional - sage.rings.padics
                1 + O(7^20)

                sage: x = K(0, -5); x                                                                                   # optional - sage.rings.padics
                O(7^-5)
                sage: x.denominator()                                                                                   # optional - sage.rings.padics
                7^5 + O(7^25)
            """
            val = self.valuation()
            R = self.parent().integer_ring()
            if val >= 0:
                return R(1)
            else:
                return R(1) << (-val)

        def numerator(self):
            """
            Return the numerator of this element, normalized in such a
            way that `x = x.numerator() / x.denominator()` always holds
            true.

            EXAMPLES::

                sage: K = Qp(7, 5)                                                                                      # optional - sage.rings.padics
                sage: x = K(1/21)                                                                                       # optional - sage.rings.padics
                sage: x.numerator()                                                                                     # optional - sage.rings.padics
                5 + 4*7 + 4*7^2 + 4*7^3 + 4*7^4 + O(7^5)

                sage: x == x.numerator() / x.denominator()                                                              # optional - sage.rings.padics
                True

            Note that the numerator lives in the ring of integers::

                sage: x.numerator().parent()                                                                            # optional - sage.rings.padics
                7-adic Ring with capped relative precision 5

            TESTS::

                sage: x = K(0, -5); x                                                                                   # optional - sage.rings.padics
                O(7^-5)
                sage: x.numerator()                                                                                     # optional - sage.rings.padics
                O(7^0)
                sage: x.denominator()                                                                                   # optional - sage.rings.padics
                7^5 + O(7^10)
            """
            R = self.parent().integer_ring()
            return R(self * self.denominator())
