"""
Places of function fields

The places of a function field correspond, one-to-one, to valuation rings of
the function field, each of which defines a discrete valuation for the elements
of the function field. "Finite" places are in one-to-one correspondence with
the prime ideals of the finite maximal order while places "at infinity" are in
one-to-one correspondence with the prime ideals of the infinite maximal order.

EXAMPLES:

All rational places of a function field can be computed::

    sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
    sage: L.<y> = K.extension(Y^3 + x + x^3*Y)                                          # needs sage.rings.function_field
    sage: L.places()                                                                    # needs sage.rings.function_field
    [Place (1/x, 1/x^3*y^2 + 1/x),
     Place (1/x, 1/x^3*y^2 + 1/x^2*y + 1),
     Place (x, y)]

The residue field associated with a place is given as an extension of the
constant field::

    sage: F.<x> = FunctionField(GF(2))
    sage: O = F.maximal_order()
    sage: p = O.ideal(x^2 + x + 1).place()                                              # needs sage.libs.pari
    sage: k, fr_k, to_k = p.residue_field()                                             # needs sage.libs.pari sage.rings.function_field
    sage: k                                                                             # needs sage.libs.pari sage.rings.function_field
    Finite Field in z2 of size 2^2

The homomorphisms are between the valuation ring and the residue field::

    sage: fr_k                                                                          # needs sage.libs.pari sage.rings.function_field
    Ring morphism:
      From: Finite Field in z2 of size 2^2
      To:   Valuation ring at Place (x^2 + x + 1)
    sage: to_k                                                                          # needs sage.libs.pari sage.rings.function_field
    Ring morphism:
      From: Valuation ring at Place (x^2 + x + 1)
      To:   Finite Field in z2 of size 2^2

AUTHORS:

- Kwankyu Lee (2017-04-30): initial version

- Brent Baccala (2019-12-20): function fields of characteristic zero
"""

# ****************************************************************************
#       Copyright (C) 2016-2022 Kwankyu Lee <ekwankyu@gmail.com>
#                     2019      Brent Baccala
#                     2021      Jonathan Kliem
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.richcmp import richcmp
from sage.categories.sets_cat import Sets


class FunctionFieldPlace(Element):
    """
    Places of function fields.

    INPUT:

    - ``parent`` -- place set of a function field

    - ``prime`` -- prime ideal associated with the place

    EXAMPLES::

        sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
        sage: L.<y> = K.extension(Y^3 + x + x^3*Y)                                      # needs sage.rings.function_field
        sage: L.places_finite()[0]                                                      # needs sage.rings.function_field
        Place (x, y)
    """
    def __init__(self, parent, prime):
        """
        Initialize the place.

        TESTS::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)                                  # needs sage.rings.function_field
            sage: p = L.places_finite()[0]                                              # needs sage.rings.function_field
            sage: TestSuite(p).run()                                                    # needs sage.rings.function_field
        """
        Element.__init__(self, parent)

        self._prime = prime

    def __hash__(self):
        """
        Return the hash of the place.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)                                  # needs sage.rings.function_field
            sage: p = L.places_finite()[0]                                              # needs sage.rings.function_field
            sage: {p: 1}                                                                # needs sage.rings.function_field
            {Place (x, y): 1}
        """
        return hash((self.function_field(), self._prime))

    def _repr_(self):
        """
        Return the string representation of the place.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x^3*Y + x)                                  # needs sage.rings.function_field
            sage: p = L.places_finite()[0]                                              # needs sage.rings.function_field
            sage: p                                                                     # needs sage.rings.function_field
            Place (x, y)
        """
        try:
            gens = self._prime.gens_two()
        except AttributeError:
            gens = self._prime.gens()
        gens_str = ', '.join(repr(g) for g in gens)
        return "Place ({})".format(gens_str)

    def _latex_(self):
        r"""
        Return the LaTeX representation of the place.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)                                  # needs sage.rings.function_field
            sage: p = L.places_finite()[0]                                              # needs sage.rings.function_field
            sage: latex(p)                                                              # needs sage.rings.function_field
            \left(y\right)
        """
        return self._prime._latex_()

    def _richcmp_(self, other, op):
        """
        Compare the place with ``other`` place.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
            sage: p1, p2, p3 = L.places()[:3]
            sage: p1 < p2
            True
            sage: p2 < p1
            False
            sage: p1 == p3
            False
        """
        from sage.rings.function_field.order import FunctionFieldOrderInfinite

        # effect that places at infinity are ordered first
        s = not isinstance(self._prime.ring(), FunctionFieldOrderInfinite)
        o = not isinstance(other._prime.ring(), FunctionFieldOrderInfinite)
        return richcmp((s, self._prime), (o, other._prime), op)

    def _acted_upon_(self, other, self_on_left):
        """
        Define integer multiplication upon the prime divisor
        of the place on the left.

        The output is a divisor.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(GF(5)); R.<Y> = PolynomialRing(K)
            sage: F.<y> = K.extension(Y^2 - x^3 - 1)
            sage: O = F.maximal_order()
            sage: I = O.ideal(x + 1, y)
            sage: P = I.place()
            sage: -3*P + 5*P
            2*Place (x + 1, y)
        """
        if self_on_left:
            raise TypeError("only left multiplication by integers is allowed")
        return other * self.divisor()

    def _neg_(self):
        """
        Return the negative of the prime divisor of this place.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x^3*Y + x)
            sage: p1, p2, p3 = L.places()[:3]
            sage: -p1 + p2
            - Place (1/x, 1/x^3*y^2 + 1/x)
             + Place (1/x, 1/x^3*y^2 + 1/x^2*y + 1)
        """
        from .divisor import divisor
        return divisor(self.function_field(), {self: -1})

    def _add_(self, other):
        """
        Return the divisor that is the sum of the place and ``other``.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x^3*Y + x)
            sage: p1, p2, p3 = L.places()[:3]
            sage: p1 + p2 + p3
            Place (1/x, 1/x^3*y^2 + 1/x)
             + Place (1/x, 1/x^3*y^2 + 1/x^2*y + 1)
             + Place (x, y)
        """
        from .divisor import prime_divisor
        return prime_divisor(self.function_field(), self) + other

    def _sub_(self, other):
        """
        Return the divisor that is this place minus ``other``.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x^3*Y + x)
            sage: p1, p2 = L.places()[:2]
            sage: p1 - p2
            Place (1/x, 1/x^3*y^2 + 1/x)
             - Place (1/x, 1/x^3*y^2 + 1/x^2*y + 1)
        """
        from .divisor import prime_divisor
        return prime_divisor(self.function_field(), self) - other

    def __radd__(self, other):
        """
        Return the prime divisor of the place if ``other`` is zero.

        This is only to support the ``sum`` function, that adds
        the argument to initial (int) zero.

        EXAMPLES::

            sage: k.<a> = GF(2)
            sage: K.<x> = FunctionField(k)
            sage: sum(K.places_finite())                                                # needs sage.libs.pari sage.modules
            Place (x) + Place (x + 1)

        Note that this does not work, as wanted::

            sage: 0 + K.place_infinite()
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for +: ...

        The reason is that the ``0`` is a Sage integer, for which
        the coercion system applies.
        """
        if other == 0:
            from .divisor import prime_divisor
            return prime_divisor(self.function_field(), self)
        raise NotImplementedError

    def function_field(self):
        """
        Return the function field to which the place belongs.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x^3*Y + x)                                  # needs sage.rings.function_field
            sage: p = L.places()[0]                                                     # needs sage.rings.function_field
            sage: p.function_field() == L                                               # needs sage.rings.function_field
            True
        """
        return self.parent()._field

    def prime_ideal(self):
        """
        Return the prime ideal associated with the place.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x^3*Y + x)                                  # needs sage.rings.function_field
            sage: p = L.places()[0]                                                     # needs sage.rings.function_field
            sage: p.prime_ideal()                                                       # needs sage.rings.function_field
            Ideal (1/x^3*y^2 + 1/x) of Maximal infinite order of Function field
            in y defined by y^3 + x^3*y + x
        """
        return self._prime

    def divisor(self, multiplicity=1):
        """
        Return the prime divisor corresponding to the place.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(GF(5)); R.<Y> = PolynomialRing(K)
            sage: F.<y> = K.extension(Y^2 - x^3 - 1)
            sage: O = F.maximal_order()
            sage: I = O.ideal(x + 1, y)
            sage: P = I.place()
            sage: P.divisor()
            Place (x + 1, y)
        """
        from .divisor import prime_divisor
        return prime_divisor(self.function_field(), self, multiplicity)


class PlaceSet(UniqueRepresentation, Parent):
    """
    Sets of Places of function fields.

    INPUT:

    - ``field`` -- function field

    EXAMPLES::

        sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
        sage: L.<y> = K.extension(Y^3 + x^3*Y + x)                                      # needs sage.rings.function_field
        sage: L.place_set()                                                             # needs sage.rings.function_field
        Set of places of Function field in y defined by y^3 + x^3*y + x
    """
    Element = FunctionFieldPlace

    def __init__(self, field):
        """
        Initialize the set of places of the function ``field``.

        TESTS::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x^3*Y + x)                                  # needs sage.rings.function_field
            sage: places = L.place_set()                                                # needs sage.rings.function_field
            sage: TestSuite(places).run()                                               # needs sage.rings.function_field
        """
        self.Element = field._place_class
        Parent.__init__(self, category=Sets().Infinite())

        self._field = field

    def _repr_(self):
        """
        Return the string representation of the place.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x^3*Y + x)                                  # needs sage.rings.function_field
            sage: L.place_set()                                                         # needs sage.rings.function_field
            Set of places of Function field in y defined by y^3 + x^3*y + x
        """
        return "Set of places of {}".format(self._field)

    def _element_constructor_(self, x):
        """
        Create a place from ``x`` if ``x`` is a prime ideal.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x^3*Y + x)
            sage: places = L.place_set()
            sage: O = L.maximal_order()
            sage: places(O.ideal(x, y))
            Place (x, y)
        """
        from .ideal import FunctionFieldIdeal

        if isinstance(x, FunctionFieldIdeal) and x.is_prime():
            return self.element_class(self, x)
        else:
            raise ValueError("not a prime ideal")

    def _an_element_(self):
        """
        Return a place.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x^3*Y + x)
            sage: places = L.place_set()
            sage: places.an_element()  # random
            Ideal (x) of Maximal order of Rational function field in x
            over Finite Field of size 2
        """
        d = 1
        while True:
            try:
                p = self._field.places(d).pop()
            except IndexError:
                d = d + 1
            else:
                break
        return p

    def function_field(self):
        """
        Return the function field to which this place set belongs.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x^3*Y + x)
            sage: PS = L.place_set()
            sage: PS.function_field() == L
            True
        """
        return self._field
