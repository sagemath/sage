r"""
Discrete Valuation Rings (DVR) and Fields (DVF)
"""
#**************************************************************************
#  Copyright (C) 2013 Xavier Caruso <xavier.caruso@normalesup.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#**************************************************************************


from sage.misc.abstract_method import abstract_method
from sage.categories.category_singleton import Category_singleton
from sage.categories.euclidean_domains import EuclideanDomains
from sage.categories.fields import Fields

class DiscreteValuationRings(Category_singleton):
    """
    The category of discrete valuation rings

    EXAMPLES::

        sage: GF(7)[['x']] in DiscreteValuationRings()
        True
        sage: TestSuite(DiscreteValuationRings()).run()
    """
    def super_categories(self):
        """
        EXAMPLES::

            sage: DiscreteValuationRings().super_categories()
            [Category of euclidean domains]
        """
        return [EuclideanDomains()]

    class ParentMethods:
        @abstract_method
        def uniformizer(self):
            """
            Return a uniformizer of this ring.

            EXAMPLES::

                sage: Zp(5).uniformizer()
                5 + O(5^21)

                sage: K.<u> = QQ[[]]
                sage: K.uniformizer()
                u
            """

        @abstract_method
        def residue_field(self):
            """
            Return the residue field of this ring.

            EXAMPLES::

                sage: Zp(5).residue_field()
                Finite Field of size 5

                sage: K.<u> = QQ[[]]
                sage: K.residue_field()
                Rational Field
            """

    class ElementMethods:
        @abstract_method
        def valuation(self):
            """
            Return the valuation of this element.

            EXAMPLES::

                sage: x = Zp(5)(50)
                sage: x.valuation()
                2
            """

        def euclidean_degree(self):
            """
            Return the Euclidean degree of this element.

            TESTS::

                sage: R.<q> = GF(5)[[]]
                sage: (q^3).euclidean_degree()
                3
                sage: R(0).euclidean_degree()
                Traceback (most recent call last):
                ...
                ValueError: Euclidean degree of the zero element not defined

            """
            if not self:
                raise ValueError("Euclidean degree of the zero element not defined")
            return self.valuation()

        def quo_rem(self, other):
            """
            Return the quotient and remainder for Euclidean division
            of ``self`` by ``other``.

            TESTS::

                sage: R.<q> = GF(5)[[]]
                sage: (q^2 + q).quo_rem(q)
                (1 + q, 0)
                sage: (q + 1).quo_rem(q^2)
                (0, 1 + q)
                sage: q.quo_rem(0)
                Traceback (most recent call last):
                ...
                ZeroDivisionError: Euclidean division by the zero element not defined

            """
            if not other:
                raise ZeroDivisionError("Euclidean division by the zero element not defined")
            P = self.parent()
            if self.valuation() >= other.valuation():
                return P(self / other), P.zero()
            else:
                return P.zero(), self

        def is_unit(self):
            """
            Return True if self is invertible.

            EXAMPLES::

                sage: x = Zp(5)(50)
                sage: x.is_unit()
                False

                sage: x = Zp(7)(50)
                sage: x.is_unit()
                True
            """
            return self.valuation() == 0

        def gcd(self,other):
            """
            Return the greatest common divisor of self and other,
            normalized so that it is a power of the distinguished
            uniformizer.
            """
            from sage.rings.infinity import Infinity
            val = min(self.valuation(), other.valuation())
            if val is Infinity:
                return self.parent()(0)
            else:
                return self.parent().uniformizer() ** val

        def lcm(self,other):
            """
            Return the least common multiple of self and other,
            normalized so that it is a power of the distinguished
            uniformizer.
            """
            from sage.rings.infinity import Infinity
            val = max(self.valuation(), other.valuation())
            if val is Infinity:
                return self.parent()(0)
            else:
                return self.parent().uniformizer() ** val


class DiscreteValuationFields(Category_singleton):
    """
    The category of discrete valuation fields

    EXAMPLES::

        sage: Qp(7) in DiscreteValuationFields()
        True
        sage: TestSuite(DiscreteValuationFields()).run()
    """

    def super_categories(self):
        """
        EXAMPLES::

            sage: DiscreteValuationFields().super_categories()
            [Category of fields]
        """
        return [Fields()]

    class ParentMethods:
        @abstract_method
        def uniformizer(self):
            """
            Return a uniformizer of this ring.

            EXAMPLES::

                sage: Qp(5).uniformizer()
                5 + O(5^21)
            """

        @abstract_method
        def residue_field(self):
            """
            Return the residue field of the ring of integers of 
            this discrete valuation field.

            EXAMPLES::

                sage: Qp(5).residue_field()
                Finite Field of size 5

                sage: K.<u> = LaurentSeriesRing(QQ)
                sage: K.residue_field()
                Rational Field
            """

    class ElementMethods:
        @abstract_method
        def valuation(self):
            """
            Return the valuation of this element.

            EXAMPLES::

                sage: x = Qp(5)(50)
                sage: x.valuation()
                2
            """
