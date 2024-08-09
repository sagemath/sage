r"""
Orders of function fields

An order of a function field is a subring that is, as a module over the base
maximal order, finitely generated and of maximal rank `n`, where `n` is the
extension degree of the function field. All orders are subrings of maximal
orders.

A rational function field has two maximal orders: maximal finite order `o` and
maximal infinite order `o_\infty`. The maximal order of a rational function
field over constant field `k` is just the polynomial ring `o=k[x]`. The
maximal infinite order is the set of rational functions whose denominator has
degree greater than or equal to that of the numerator.

EXAMPLES::

    sage: K.<x> = FunctionField(QQ)
    sage: O = K.maximal_order()
    sage: I = O.ideal(1/x); I
    Ideal (1/x) of Maximal order of Rational function field in x over Rational Field
    sage: 1/x in O
    False
    sage: Oinf = K.maximal_order_infinite()
    sage: 1/x in Oinf
    True

In an extension of a rational function field, an order over the maximal finite
order is called a finite order while an order over the maximal infinite order
is called an infinite order. Thus a function field has one maximal finite order
`O` and one maximal infinite order `O_\infty`. There are other non-maximal
orders such as equation orders::

    sage: # needs sage.rings.function_field
    sage: K.<x> = FunctionField(GF(3)); R.<y> = K[]
    sage: L.<y> = K.extension(y^3 - y - x)
    sage: O = L.equation_order()
    sage: 1/y in O
    False
    sage: x/y in O
    True

Sage provides an extensive functionality for computations in maximal orders of
function fields. For example, you can decompose a prime ideal of a rational
function field in an extension::

    sage: K.<x> = FunctionField(GF(2)); _.<t> = K[]
    sage: o = K.maximal_order()
    sage: p = o.ideal(x + 1)
    sage: p.is_prime()                                                                  # needs sage.libs.pari
    True

    sage: # needs sage.rings.function_field
    sage: F.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2)
    sage: O = F.maximal_order()
    sage: O.decomposition(p)
    [(Ideal (x + 1, y + 1) of Maximal order
     of Function field in y defined by y^3 + x^6 + x^4 + x^2, 1, 1),
     (Ideal (x + 1, (1/(x^3 + x^2 + x))*y^2 + y + 1) of Maximal order
     of Function field in y defined by y^3 + x^6 + x^4 + x^2, 2, 1)]

    sage: # needs sage.rings.function_field
    sage: p1, relative_degree,ramification_index = O.decomposition(p)[1]
    sage: p1.parent()
    Monoid of ideals of Maximal order of Function field in y
    defined by y^3 + x^6 + x^4 + x^2
    sage: relative_degree
    2
    sage: ramification_index
    1

When the base constant field is the algebraic field `\QQbar`, the only prime ideals
of the maximal order of the rational function field are linear polynomials. ::

    sage: # needs sage.rings.function_field sage.rings.number_field
    sage: K.<x> = FunctionField(QQbar)
    sage: R.<y> = K[]
    sage: L.<y> = K.extension(y^2 - (x^3-x^2))
    sage: p = K.maximal_order().ideal(x)
    sage: L.maximal_order().decomposition(p)
    [(Ideal (1/x*y - I) of Maximal order of Function field in y defined by y^2 - x^3 + x^2,
      1,
      1),
     (Ideal (1/x*y + I) of Maximal order of Function field in y defined by y^2 - x^3 + x^2,
      1,
      1)]


AUTHORS:

- William Stein (2010): initial version

- Maarten Derickx (2011-09-14): fixed ideal_with_gens_over_base() for rational function fields

- Julian Rueth (2011-09-14): added check in _element_constructor_

- Kwankyu Lee (2017-04-30): added maximal orders of global function fields

- Brent Baccala (2019-12-20): support orders in characteristic zero
"""

# ****************************************************************************
#       Copyright (C) 2010      William Stein <wstein@gmail.com>
#                     2011      Maarten Derickx <m.derickx.student@gmail.com>
#                     2011      Julian Rueth <julian.rueth@gmail.com>
#                     2017-2020 Kwankyu Lee
#                     2019      Brent Baccala
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.integral_domains import IntegralDomains
from sage.structure.parent import Parent
from sage.structure.unique_representation import CachedRepresentation, UniqueRepresentation

from .ideal import IdealMonoid, FunctionFieldIdeal


class FunctionFieldOrder_base(CachedRepresentation, Parent):
    """
    Base class for orders in function fields.

    INPUT:

    - ``field`` -- function field

    EXAMPLES::

        sage: F = FunctionField(QQ,'y')
        sage: F.maximal_order()
        Maximal order of Rational function field in y over Rational Field
    """
    def __init__(self, field, ideal_class=FunctionFieldIdeal, category=None):
        """
        Initialize.

        TESTS::

            sage: F = FunctionField(QQ,'y')
            sage: O = F.maximal_order()
            sage: TestSuite(O).run(skip='_test_gcd_vs_xgcd')
        """
        category = IntegralDomains().or_subcategory(category).Infinite()
        Parent.__init__(self, category=category, facade=field)

        self._ideal_class = ideal_class # element class for parent ideal monoid
        self._field = field

    def is_field(self, proof=True):
        """
        Return ``False`` since orders are never fields.

        EXAMPLES::

            sage: FunctionField(QQ,'y').maximal_order().is_field()
            False
        """
        return False

    def is_noetherian(self):
        """
        Return ``True`` since orders in function fields are Noetherian.

        EXAMPLES::

            sage: FunctionField(QQ,'y').maximal_order().is_noetherian()
            True
        """
        return True

    def function_field(self):
        """
        Return the function field to which the order belongs.

        EXAMPLES::

            sage: FunctionField(QQ,'y').maximal_order().function_field()
            Rational function field in y over Rational Field
        """
        return self._field

    fraction_field = function_field

    def is_subring(self, other):
        """
        Return ``True`` if the order is a subring of the other order.

        INPUT:

        - ``other`` -- order of the function field or the field itself

        EXAMPLES::

            sage: F = FunctionField(QQ,'y')
            sage: O = F.maximal_order()
            sage: O.is_subring(F)
            True
        """
        if other is self._field:
            return True
        else:
            raise NotImplementedError

    def ideal_monoid(self):
        """
        Return the monoid of ideals of the order.

        EXAMPLES::

            sage: FunctionField(QQ,'y').maximal_order().ideal_monoid()
            Monoid of ideals of Maximal order of Rational function field in y over Rational Field
        """
        return IdealMonoid(self)


class FunctionFieldOrder(FunctionFieldOrder_base):
    """
    Base class for orders in function fields.
    """
    def _repr_(self):
        """
        Return the string representation.

        EXAMPLES::

            sage: FunctionField(QQ,'y').maximal_order()
            Maximal order of Rational function field in y over Rational Field
        """
        return "Order in {}".format(self._field)


class FunctionFieldOrderInfinite(FunctionFieldOrder_base):
    """
    Base class for infinite orders in function fields.
    """
    def _repr_(self):
        """
        EXAMPLES::

            sage: FunctionField(QQ,'y').maximal_order_infinite()
            Maximal infinite order of Rational function field in y over Rational Field
        """
        return "Infinite order in {}".format(self.function_field())


class FunctionFieldMaximalOrder(UniqueRepresentation, FunctionFieldOrder):
    """
    Base class of maximal orders of function fields.
    """
    def _repr_(self):
        """
        Return the string representation of the order.

        EXAMPLES::

            sage: FunctionField(QQ,'y').maximal_order()._repr_()
            'Maximal order of Rational function field in y over Rational Field'
        """
        return "Maximal order of %s" % (self.function_field(),)


class FunctionFieldMaximalOrderInfinite(FunctionFieldMaximalOrder, FunctionFieldOrderInfinite):
    """
    Base class of maximal infinite orders of function fields.
    """
    def _repr_(self):
        """
        EXAMPLES::

            sage: FunctionField(QQ,'y').maximal_order_infinite()
            Maximal infinite order of Rational function field in y over Rational Field

            sage: K.<x> = FunctionField(GF(2)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3 - x^2*(x^2+x+1)^2)                            # needs sage.rings.function_field
            sage: F.maximal_order_infinite()                                            # needs sage.rings.function_field
            Maximal infinite order of Function field in y defined by y^3 + x^6 + x^4 + x^2
        """
        return "Maximal infinite order of %s" % (self.function_field(),)
