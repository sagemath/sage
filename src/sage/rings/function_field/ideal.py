r"""
Ideals of function fields

Ideals of an order of a function field include all fractional ideals of the order.
Sage provides basic arithmetic with fractional ideals.

The fractional ideals of the maximal order of a global function field forms a multiplicative
monoid. Sage allows advanced arithmetic with the fractional ideals. For example, an ideal
of the maximal order can be factored into a product of prime ideals.

EXAMPLES:

Ideals in the maximal order of a rational function field::

    sage: K.<x> = FunctionField(QQ)
    sage: O = K.maximal_order()
    sage: I = O.ideal(x^3 + 1); I
    Ideal (x^3 + 1) of Maximal order of Rational function field in x over Rational Field
    sage: I^2
    Ideal (x^6 + 2*x^3 + 1) of Maximal order of Rational function field in x over Rational Field
    sage: ~I
    Ideal (1/(x^3 + 1)) of Maximal order of Rational function field in x over Rational Field
    sage: ~I * I
    Ideal (1) of Maximal order of Rational function field in x over Rational Field

Ideals in the equation order of an extension of a rational function field::

    sage: K.<x> = FunctionField(QQ); R.<y> = K[]
    sage: L.<y> = K.extension(y^2 - x^3 - 1)                                                                            # needs sage.rings.function_field
    sage: O = L.equation_order()                                                                                        # needs sage.rings.function_field
    sage: I = O.ideal(y); I                                                                                             # needs sage.rings.function_field
    Ideal (x^3 + 1, -y) of Order in Function field in y defined by y^2 - x^3 - 1
    sage: I^2                                                                                                           # needs sage.rings.function_field
    Ideal (x^3 + 1, (-x^3 - 1)*y) of Order in Function field in y defined by y^2 - x^3 - 1

Ideals in the maximal order of a global function field::

    sage: # needs sage.rings.finite_rings
    sage: K.<x> = FunctionField(GF(2)); R.<y> = K[]
    sage: L.<y> = K.extension(y^2 - x^3*y - x)                                                                          # needs sage.rings.function_field
    sage: O = L.maximal_order()                                                                                         # needs sage.rings.function_field
    sage: I = O.ideal(y)                                                                                                # needs sage.rings.function_field
    sage: I^2                                                                                                           # needs sage.rings.function_field
    Ideal (x) of Maximal order of Function field in y defined by y^2 + x^3*y + x
    sage: ~I                                                                                                            # needs sage.rings.function_field
    Ideal (1/x*y) of Maximal order of Function field in y defined by y^2 + x^3*y + x
    sage: ~I * I                                                                                                        # needs sage.rings.function_field
    Ideal (1) of Maximal order of Function field in y defined by y^2 + x^3*y + x

    sage: J = O.ideal(x + y) * I                                                                                        # needs sage.rings.finite_rings sage.rings.function_field
    sage: J.factor()                                                                                                    # needs sage.rings.finite_rings sage.rings.function_field
    (Ideal (y) of Maximal order of Function field in y defined by y^2 + x^3*y + x)^2 *
    (Ideal (x^3 + x + 1, y + x) of Maximal order of Function field in y defined by y^2 + x^3*y + x)

Ideals in the maximal infinite order of a global function field::

    sage: # needs sage.rings.finite_rings
    sage: K.<x> = FunctionField(GF(3^2)); R.<t> = K[]
    sage: F.<y> = K.extension(t^3 + t^2 - x^4)                                                                          # needs sage.rings.function_field
    sage: Oinf = F.maximal_order_infinite()                                                                             # needs sage.rings.function_field
    sage: I = Oinf.ideal(1/y)                                                                                           # needs sage.rings.function_field
    sage: I + I == I
    True
    sage: I^2                                                                                                           # needs sage.rings.function_field
    Ideal (1/x^4*y) of Maximal infinite order of Function field in y defined by y^3 + y^2 + 2*x^4
    sage: ~I                                                                                                            # needs sage.rings.function_field
    Ideal (y) of Maximal infinite order of Function field in y defined by y^3 + y^2 + 2*x^4
    sage: ~I * I                                                                                                        # needs sage.rings.function_field
    Ideal (1) of Maximal infinite order of Function field in y defined by y^3 + y^2 + 2*x^4
    sage: I.factor()                                                                                                    # needs sage.rings.function_field
    (Ideal (1/x^3*y^2) of Maximal infinite order of Function field in y defined by y^3 + y^2 + 2*x^4)^4

AUTHORS:

- William Stein (2010): initial version

- Maarten Derickx (2011-09-14): fixed ideal_with_gens_over_base()

- Kwankyu Lee (2017-04-30): added ideals for global function fields
"""

# ****************************************************************************
#       Copyright (C) 2010      William Stein <wstein@gmail.com>
#                     2011      Maarten Derickx <m.derickx.student@gmail.com>
#                     2017-2021 Kwankyu Lee
#                     2018      Frédéric Chapoton
#                     2019      Brent Baccala
#                     2021      Jonathan Kliem
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.latex import latex
from sage.combinat.subset import powerset
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.structure.richcmp import richcmp
from sage.structure.factorization import Factorization
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.monoids import Monoids
from sage.rings.ideal import Ideal_generic


class FunctionFieldIdeal(Element):
    """
    Base class of fractional ideals of function fields.

    INPUT:

    - ``ring`` -- ring of the ideal

    EXAMPLES::

        sage: K.<x> = FunctionField(GF(7))
        sage: O = K.equation_order()
        sage: O.ideal(x^3 + 1)
        Ideal (x^3 + 1) of Maximal order of Rational function field in x over Finite Field of size 7
    """
    def __init__(self, ring):
        """
        Initialize.

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(7))
            sage: O = K.equation_order()
            sage: I = O.ideal(x^3 + 1)
            sage: TestSuite(I).run()
        """
        Element.__init__(self, ring.ideal_monoid())
        self._ring = ring

    def _repr_short(self):
        """
        Return a string representation of this ideal that doesn't
        include the name of the ambient ring.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(3^2)); R.<t> = K[]
            sage: F.<y> = K.extension(t^3 + t^2 - x^4)                                                                  # needs sage.rings.function_field
            sage: Oinf = F.maximal_order_infinite()                                                                     # needs sage.rings.function_field
            sage: I = Oinf.ideal(1/y)                                                                                   # needs sage.rings.function_field
            sage: I._repr_short()                                                                                       # needs sage.rings.function_field
            '(1/x^4*y^2)'
        """
        if self.is_zero():
            return "(0)"

        return "(%s)" % (', '.join([repr(g) for g in self.gens_reduced()]), )

    def _repr_(self):
        """
        Return a string representation of this ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: O = K.maximal_order()
            sage: I = O.ideal(x, 1/(x+1)); I
            Ideal (1/(x + 1)) of Maximal order of Rational function field in x over Rational Field

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.equation_order()
            sage: O.ideal(x^2 + 1)
            Ideal (x^2 + 1) of Order in Function field in y defined by y^2 - x^3 - 1

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3*Y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y); I
            Ideal (y) of Maximal order of Function field in y defined by y^2 + x^3*y + x

            sage: # needs sage.rings.function_field
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y); I
            Ideal (y) of Maximal order of Function field in y
            defined by y^2 + y + (x^2 + 1)/x

            sage: K.<x> = FunctionField(GF(2))
            sage: Oinf = K.maximal_order_infinite()
            sage: I = Oinf.ideal(x/(x^2+1))
            sage: I
            Ideal (1/x) of Maximal infinite order of Rational function field
            in x over Finite Field of size 2

            sage: # needs sage.rings.finite_rings sage.rings.function_field
            sage: K.<x> = FunctionField(GF(3^2)); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3 + t^2 - x^4)
            sage: Oinf = F.maximal_order_infinite()
            sage: Oinf.ideal(1/y)
            Ideal (1/x^4*y^2) of Maximal infinite order of Function field
            in y defined by y^3 + y^2 + 2*x^4

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: Oinf.ideal(1/y)
            Ideal (1/x*y) of Maximal infinite order of Function field in y
            defined by y^2 + y + (x^2 + 1)/x
        """
        if self.is_zero():
            return "Zero ideal of %s" % (self._ring,)

        return "Ideal %s of %s" % (self._repr_short(), self.ring())

    def _latex_(self):
        r"""
        Return the LaTeX representation of the ideal.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 - x^3*Y - x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: latex(I)
            \left(y\right)
        """
        return '\\left(' + ', '.join(latex(g) for g in self.gens_reduced()) + '\\right)'

    def _div_(self, other):
        """
        Return the ideal divided by the ``other`` ideal.

        INPUT:

        - ``other`` -- ideal

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7))
            sage: O = K.equation_order()
            sage: I = O.ideal(x^3 + 1)
            sage: I / I
            Ideal (1) of Maximal order of Rational function field in x
            over Finite Field of size 7
        """
        return self * ~other

    def gens_reduced(self):
        r"""
        Return reduced generators.

        For now, this method just looks at the generators and sees if any
        can be removed without changing the ideal.  It prefers principal
        representations (a single generator) over all others, and otherwise
        picks the generator set with the shortest print representation.

        This method is provided so that ideals in function fields have
        the method :meth:`gens_reduced()`, just like ideals of number
        fields. Sage linear algebra machinery sometimes requires this.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7))
            sage: O = K.equation_order()
            sage: I = O.ideal(x, x^2, x^2 + x)
            sage: I.gens_reduced()
            (x,)
        """
        gens = self.gens()
        if len(gens) == 1:
            return gens
        candidate_gensets = []
        for genset in powerset(gens):
            if self.parent()(genset) == self:
                candidate_gensets.append(genset)
        candidate_gensets.sort(key=lambda item: (len(item), len(repr(item)), item))
        return candidate_gensets[0]

    def ring(self):
        """
        Return the ring to which this ideal belongs.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7))
            sage: O = K.equation_order()
            sage: I = O.ideal(x, x^2, x^2 + x)
            sage: I.ring()
            Maximal order of Rational function field in x over Finite Field of size 7
        """
        return self._ring

    def base_ring(self):
        r"""
        Return the base ring of this ideal.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.equation_order()
            sage: I = O.ideal(x^2 + 1)
            sage: I.base_ring()
            Order in Function field in y defined by y^2 - x^3 - 1
        """
        return self.ring()

    def place(self):
        """
        Return the place associated with this prime ideal.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(4))
            sage: O = K.maximal_order()
            sage: I = O.ideal(x^2 + x + 1)
            sage: I.place()
            Traceback (most recent call last):
            ...
            TypeError: not a prime ideal
            sage: I = O.ideal(x^3 + x + 1)
            sage: I.place()
            Place (x^3 + x + 1)

            sage: K.<x> = FunctionField(GF(2))
            sage: Oinf = K.maximal_order_infinite()
            sage: I = Oinf.ideal((x + 1)/(x^3 + 1))
            sage: p = I.factor()[0][0]
            sage: p.place()
            Place (1/x)

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(GF(2)); _.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3 - x^2*(x^2+x+1)^2)
            sage: O = F.maximal_order()
            sage: I = O.ideal(y)
            sage: [f.place() for f,_ in I.factor()]
            [Place (x, (1/(x^3 + x^2 + x))*y^2),
             Place (x^2 + x + 1, (1/(x^3 + x^2 + x))*y^2)]

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: [f.place() for f,_ in I.factor()]
            [Place (x, x*y), Place (x + 1, x*y)]

            sage: # needs sage.rings.finite_rings sage.rings.function_field
            sage: K.<x> = FunctionField(GF(3^2)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^3 + t^2 - x^4)
            sage: Oinf = F.maximal_order_infinite()
            sage: I = Oinf.ideal(1/x)
            sage: I.factor()
            (Ideal (1/x^3*y^2) of Maximal infinite order of Function field
            in y defined by y^3 + y^2 + 2*x^4)^3
            sage: J = I.factor()[0][0]
            sage: J.is_prime()
            True
            sage: J.place()
            Place (1/x, 1/x^3*y^2)

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: Oinf = L.maximal_order_infinite()
            sage: I = Oinf.ideal(1/x)
            sage: I.factor()
            (Ideal (1/x*y) of Maximal infinite order of Function field in y
            defined by y^2 + y + (x^2 + 1)/x)^2
            sage: J = I.factor()[0][0]
            sage: J.is_prime()
            True
            sage: J.place()
            Place (1/x, 1/x*y)
        """
        if not self.is_prime():
            raise TypeError("not a prime ideal")

        place_set = self.ring().fraction_field().place_set()
        return place_set.element_class(place_set, self)

    def factor(self):
        """
        Return the factorization of this ideal.

        Subclass of this class should define :meth:`_factor` method that
        returns a list of prime ideal and multiplicity pairs.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(4))
            sage: O = K.maximal_order()
            sage: I = O.ideal(x^3*(x + 1)^2)
            sage: I.factor()
            (Ideal (x) of Maximal order of Rational function field in x
            over Finite Field in z2 of size 2^2)^3 *
            (Ideal (x + 1) of Maximal order of Rational function field in x
            over Finite Field in z2 of size 2^2)^2

            sage: # needs sage.rings.finite_rings
            sage: Oinf = K.maximal_order_infinite()
            sage: I = Oinf.ideal((x + 1)/(x^3 + 1))
            sage: I.factor()
            (Ideal (1/x) of Maximal infinite order of Rational function field in x
            over Finite Field in z2 of size 2^2)^2

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(GF(2)); _.<T> = PolynomialRing(K)
            sage: F.<y> = K.extension(T^3 - x^2*(x^2 + x + 1)^2)
            sage: O = F.maximal_order()
            sage: I = O.ideal(y)
            sage: I == I.factor().prod()
            True

            sage: # needs sage.rings.function_field
            sage: Oinf = F.maximal_order_infinite()
            sage: f= 1/x
            sage: I = Oinf.ideal(f)
            sage: I.factor()
            (Ideal ((1/(x^4 + x^3 + x^2))*y^2 + 1/x^2*y + 1) of Maximal infinite order
            of Function field in y defined by y^3 + x^6 + x^4 + x^2) *
            (Ideal ((1/(x^4 + x^3 + x^2))*y^2 + 1) of Maximal infinite order
            of Function field in y defined by y^3 + x^6 + x^4 + x^2)

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ); _.<Y> = K[]
            sage: F.<y> = K.extension(Y^3 - x^2*(x^2 + x + 1)^2)
            sage: O = F.maximal_order()
            sage: I = O.ideal(y)
            sage: I == I.factor().prod()
            True

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: I == I.factor().prod()
            True
        """
        return Factorization(self._factor(), cr=True)

    def divisor(self):
        """
        Return the divisor corresponding to the ideal.

        EXAMPLES::

            sage: # needs sage.modules sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(4))
            sage: O = K.maximal_order()
            sage: I = O.ideal(x*(x + 1)^2/(x^2 + x + 1))
            sage: I.divisor()
            Place (x) + 2*Place (x + 1) - Place (x + z2) - Place (x + z2 + 1)

            sage: # needs sage.modules sage.rings.finite_rings
            sage: Oinf = K.maximal_order_infinite()
            sage: I = Oinf.ideal((x + 1)/(x^3 + 1))
            sage: I.divisor()
            2*Place (1/x)

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(GF(2)); _.<T> = PolynomialRing(K)
            sage: F.<y> = K.extension(T^3 - x^2*(x^2 + x + 1)^2)
            sage: O = F.maximal_order()
            sage: I = O.ideal(y)
            sage: I.divisor()
            2*Place (x, (1/(x^3 + x^2 + x))*y^2)
             + 2*Place (x^2 + x + 1, (1/(x^3 + x^2 + x))*y^2)

            sage: # needs sage.rings.function_field
            sage: Oinf = F.maximal_order_infinite()
            sage: I = Oinf.ideal(y)
            sage: I.divisor()
            -2*Place (1/x, 1/x^4*y^2 + 1/x^2*y + 1)
             - 2*Place (1/x, 1/x^2*y + 1)

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: I.divisor()
            - Place (x, x*y)
             + 2*Place (x + 1, x*y)

            sage: # needs sage.rings.function_field
            sage: Oinf = L.maximal_order_infinite()
            sage: I = Oinf.ideal(y)
            sage: I.divisor()
            - Place (1/x, 1/x*y)
        """
        from .divisor import divisor

        if self.is_zero():
            raise ValueError("not defined for zero ideal")

        F = self.ring().fraction_field()
        data = {prime.place(): multiplicity for prime, multiplicity in self._factor()}
        return divisor(F, data)

    def divisor_of_zeros(self):
        """
        Return the divisor of zeros corresponding to the ideal.

        EXAMPLES::

            sage: # needs sage.modules sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(4))
            sage: O = K.maximal_order()
            sage: I = O.ideal(x*(x + 1)^2/(x^2 + x + 1))
            sage: I.divisor_of_zeros()
            Place (x) + 2*Place (x + 1)

            sage: # needs sage.modules
            sage: K.<x> = FunctionField(GF(2))
            sage: Oinf = K.maximal_order_infinite()
            sage: I = Oinf.ideal((x + 1)/(x^3 + 1))
            sage: I.divisor_of_zeros()
            2*Place (1/x)

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: I.divisor_of_zeros()
            2*Place (x + 1, x*y)
        """
        from .divisor import divisor

        if self.is_zero():
            raise ValueError("not defined for zero ideal")

        F = self.ring().fraction_field()
        data = {prime.place(): multiplicity for prime, multiplicity in self._factor() if multiplicity > 0}
        return divisor(F, data)

    def divisor_of_poles(self):
        """
        Return the divisor of poles corresponding to the ideal.

        EXAMPLES::

            sage: # needs sage.modules sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(4))
            sage: O = K.maximal_order()
            sage: I = O.ideal(x*(x + 1)^2/(x^2 + x + 1))
            sage: I.divisor_of_poles()
            Place (x + z2) + Place (x + z2 + 1)

            sage: # needs sage.modules
            sage: K.<x> = FunctionField(GF(2))
            sage: Oinf = K.maximal_order_infinite()
            sage: I = Oinf.ideal((x + 1)/(x^3 + 1))
            sage: I.divisor_of_poles()
            0

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: O = L.maximal_order()
            sage: I = O.ideal(y)
            sage: I.divisor_of_poles()
            Place (x, x*y)
        """
        from .divisor import divisor

        if self.is_zero():
            raise ValueError("not defined for zero ideal")

        F = self.ring().fraction_field()
        data = {prime.place(): - multiplicity for prime, multiplicity in self._factor() if multiplicity < 0}
        return divisor(F, data)


class FunctionFieldIdeal_module(FunctionFieldIdeal, Ideal_generic):
    """
    A fractional ideal specified by a finitely generated module over
    the integers of the base field.

    INPUT:

    - ``ring`` -- an order in a function field

    - ``module`` -- a module of the order

    EXAMPLES:

    An ideal in an extension of a rational function field::

        sage: # needs sage.rings.function_field
        sage: K.<x> = FunctionField(QQ); R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x^3 - 1)
        sage: O = L.equation_order()
        sage: I = O.ideal(y)
        sage: I
        Ideal (x^3 + 1, -y) of Order in Function field in y defined by y^2 - x^3 - 1
        sage: I^2
        Ideal (x^3 + 1, (-x^3 - 1)*y) of Order in Function field in y defined by y^2 - x^3 - 1
    """
    def __init__(self, ring, module):
        """
        Initialize.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.equation_order()
            sage: I = O.ideal(y)
            sage: TestSuite(I).run()
        """
        FunctionFieldIdeal.__init__(self, ring)

        self._module = module
        self._structure = ring.fraction_field().vector_space()

        V, from_V, to_V = self._structure
        self._gens = tuple([from_V(a) for a in module.basis()])

        # module generators are still ideal generators
        Ideal_generic.__init__(self, ring, self._gens, coerce=False)

    def __contains__(self, x):
        """
        Return ``True`` if ``x`` is in this ideal.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.equation_order()
            sage: I = O.ideal(y); I
            Ideal (x^3 + 1, -y) of Order in Function field in y defined by y^2 - x^3 - 1
            sage: y in I
            True
            sage: y/x in I
            False
            sage: y^2 - 2 in I
            False
        """
        return self._structure[2](x) in self._module

    def __hash__(self):
        """
        Return the hash of this ideal.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.equation_order()
            sage: I = O.ideal(y)
            sage: d = {I: 1}  # indirect doctest
        """
        return hash((self._ring,self._module))

    def _richcmp_(self, other, op):
        """
        Compare this ideal with the ``other`` ideal with respect to ``op``.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.equation_order()
            sage: I = O.ideal(x, y);  J = O.ideal(y^2 - 2)
            sage: I + J == J + I  # indirect test
            True
            sage: I + I == I  # indirect doctest
            True
            sage: I == J
            False
            sage: I < J
            True
            sage: J < I
            False
        """
        return richcmp(self.module().basis(), other.module().basis(), op)

    def module(self):
        """
        Return the module over the maximal order of the base field that
        underlies this ideal.

        The formation of the module is compatible with the vector
        space corresponding to the function field.

        OUTPUT: a module over the maximal order of the base field of the ideal

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.equation_order();  O
            Order in Function field in y defined by y^2 - x^3 - 1
            sage: I = O.ideal(x^2 + 1)
            sage: I.gens()
            (x^2 + 1, (x^2 + 1)*y)
            sage: I.module()
            Free module of degree 2 and rank 2 over Maximal order of Rational function field in x over Rational Field
            Echelon basis matrix:
            [x^2 + 1       0]
            [      0 x^2 + 1]
            sage: V, from_V, to_V = L.vector_space(); V
            Vector space of dimension 2 over Rational function field in x over Rational Field
            sage: I.module().is_submodule(V)
            True
        """
        return self._module

    def gens(self):
        """
        Return a set of generators of this ideal.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.equation_order()
            sage: I = O.ideal(x^2 + 1)
            sage: I.gens()
            (x^2 + 1, (x^2 + 1)*y)
        """
        return self._gens

    def gen(self, i):
        """
        Return the ``i``-th generator in the current basis of this ideal.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.equation_order()
            sage: I = O.ideal(x^2 + 1)
            sage: I.gen(1)
            (x^2 + 1)*y
        """
        return self._gens[i]

    def ngens(self):
        """
        Return the number of generators in the basis.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.equation_order()
            sage: I = O.ideal(x^2 + 1)
            sage: I.ngens()
            2
        """
        return len(self._gens)

    def _add_(self, other):
        """
        Add this ideal with the ``other`` ideal.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.equation_order()
            sage: I = O.ideal(y)
            sage: J = O.ideal(x+y)
            sage: I + J
            Ideal ((-x^2 + x)*y + 1, -y) of Order in Function field in y defined by y^2 - x^3 - 1
        """
        return self.ring().ideal(self.gens() + other.gens())

    def _mul_(self, other):
        """
        Multiply this ideal with the ``other`` ideal.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.equation_order()
            sage: I = O.ideal(y)
            sage: J = O.ideal(x+y)
            sage: I * J
            Ideal ((-x^5 + x^4 - x^2 + x)*y + x^3 + 1, (x^3 - x^2 + 1)*y) of Order in Function field in y defined by y^2 - x^3 - 1
        """
        return self.ring().ideal([x*y for x in self.gens() for y in other.gens()])

    def _acted_upon_(self, other, on_left):
        """
        Multiply this ideal on the right with ``other``.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.equation_order()
            sage: I = O.ideal(y)
            sage: x * I
            Ideal (x^4 + x, -x*y) of Order in Function field in y defined by y^2 - x^3 - 1
        """
        return self.ring().ideal([other * x for x in self.gens()])

    def intersection(self, other):
        """
        Return the intersection of this ideal and ``other``.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.equation_order()
            sage: I = O.ideal(y^3); J = O.ideal(y^2)
            sage: Z = I.intersection(J); Z
            Ideal (x^6 + 2*x^3 + 1, (-x^3 - 1)*y) of Order in Function field in y defined by y^2 - x^3 - 1
            sage: y^2 in Z
            False
            sage: y^3 in Z
            True
        """
        if not isinstance(other, FunctionFieldIdeal):
            try:
                if self.ring().has_coerce_map_from(other):
                    return self
            except (TypeError,ArithmeticError,ValueError):
                pass
            other = self.ring().ideal(other)

        basis = self.module().intersection(other.module()).basis()

        V, from_V, to_V = self._structure
        return self.ring().ideal_with_gens_over_base([from_V(a) for a in basis])

    def __invert__(self):
        """
        Return the inverse of this ideal.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.equation_order()
            sage: I = O.ideal(y)
            sage: ~I
            Ideal (-1, (1/(x^3 + 1))*y) of Order in Function field in y defined by y^2 - x^3 - 1
            sage: I^-1
            Ideal (-1, (1/(x^3 + 1))*y) of Order in Function field in y defined by y^2 - x^3 - 1
            sage: ~I * I
            Ideal (1) of Order in Function field in y defined by y^2 - x^3 - 1
        """
        if len(self.gens()) == 0:
            raise ZeroDivisionError

        # NOTE: If  I = (g0, ..., gn), then {x : x*I is in R}
        # is the intersection over i of {x : x*gi is in R}
        # Thus (I + J)^(-1) = I^(-1) intersect J^(-1).

        G = self.gens()
        R = self.ring()
        inv = R.ideal(~G[0])
        for g in G[1:]:
            inv = inv.intersection(R.ideal(~g))
        return inv


class FunctionFieldIdealInfinite(FunctionFieldIdeal):
    """
    Base class of ideals of maximal infinite orders
    """
    pass


class FunctionFieldIdealInfinite_module(FunctionFieldIdealInfinite, Ideal_generic):
    """
    A fractional ideal specified by a finitely generated module over
    the integers of the base field.

    INPUT:

    - ``ring`` -- order in a function field

    - ``module`` -- module

    EXAMPLES::

        sage: # needs sage.rings.function_field
        sage: K.<x> = FunctionField(QQ); R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x^3 - 1)
        sage: O = L.equation_order()
        sage: O.ideal(y)
        Ideal (x^3 + 1, -y) of Order in Function field in y defined by y^2 - x^3 - 1
    """
    def __init__(self, ring, module):
        """
        Initialize.

        TESTS::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.equation_order()
            sage: I = O.ideal(y)
            sage: TestSuite(I).run()
        """
        FunctionFieldIdealInfinite.__init__(self, ring)

        self._module = module
        self._structure = ring.fraction_field().vector_space()

        V, from_V, to_V = self._structure
        gens = tuple([from_V(a) for a in module.basis()])
        self._gens = gens

        # module generators are still ideal generators
        Ideal_generic.__init__(self, ring, self._gens, coerce=False)

    def __contains__(self, x):
        """
        Return ``True`` if ``x`` is in this ideal.

        INPUT:

        - ``x`` -- element of the function field

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.equation_order()
            sage: I = O.ideal_with_gens_over_base([1, y]);  I
            Ideal (1) of Order in Function field in y defined by y^2 + 6*x^3 + 6
            sage: y in I
            True
            sage: y/x in I
            False
            sage: y^2 - 2 in I
            True
        """
        return self._structure[2](x) in self._module

    def __hash__(self):
        """
        Return the hash of this ideal.

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.equation_order()
            sage: I = O.ideal_with_gens_over_base([1, y])
            sage: d = {I: 2}  # indirect doctest
        """
        return hash((self._ring,self._module))

    def __eq__(self, other):
        """
        Test equality of this ideal with the ``other`` ideal.

        INPUT:

        - ``other`` -- ideal

        EXAMPLES::

            sage: # needs sage.rings.function_field
            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x^3 - 1)
            sage: O = L.equation_order()
            sage: I = O.ideal_with_gens_over_base([1, y])
            sage: I == I + I  # indirect doctest
            True
        """
        if not isinstance(other, FunctionFieldIdeal_module):
            other = self.ring().ideal(other)
        if self.ring() != other.ring():
            raise ValueError("rings must be the same")

        if (self.module().is_submodule(other.module()) and
            other.module().is_submodule(self.module())):
            return True
        else:
            return False

    def module(self):
        """
        Return the module over the maximal order of the base field that
        underlies this ideal.

        The formation of the module is compatible with the vector
        space corresponding to the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7))
            sage: O = K.maximal_order(); O
            Maximal order of Rational function field in x over Finite Field of size 7
            sage: K.polynomial_ring()
            Univariate Polynomial Ring in x over
             Rational function field in x over Finite Field of size 7
            sage: I = O.ideal([x^2 + 1, x*(x^2+1)])
            sage: I.gens()
            (x^2 + 1,)
            sage: I.module()                                                            # needs sage.modules
            Free module of degree 1 and rank 1 over
             Maximal order of Rational function field in x over Finite Field of size 7
            Echelon basis matrix:
            [x^2 + 1]
            sage: V, from_V, to_V = K.vector_space(); V                                 # needs sage.modules
            Vector space of dimension 1 over
             Rational function field in x over Finite Field of size 7
            sage: I.module().is_submodule(V)                                            # needs sage.modules
            True
        """
        return self._module


class IdealMonoid(UniqueRepresentation, Parent):
    r"""
    The monoid of ideals in orders of function fields.

    INPUT:

    - ``R`` -- order

    EXAMPLES::

        sage: K.<x> = FunctionField(GF(2))
        sage: O = K.maximal_order()
        sage: M = O.ideal_monoid(); M
        Monoid of ideals of Maximal order of Rational function field in x over Finite Field of size 2
    """

    def __init__(self, R):
        """
        Initialize the ideal monoid.

        TESTS::

            sage: K.<x> = FunctionField(GF(2))
            sage: O = K.maximal_order()
            sage: M = O.ideal_monoid()
            sage: TestSuite(M).run()
        """
        self.Element = R._ideal_class
        Parent.__init__(self, category=Monoids())

        self.__R = R
        self._populate_coercion_lists_()

    def _repr_(self):
        """
        Return the string representation of the ideal monoid.

        TESTS::

            sage: K.<x> = FunctionField(GF(2))
            sage: O = K.maximal_order()
            sage: M = O.ideal_monoid(); M._repr_()
            'Monoid of ideals of Maximal order of Rational function field in x over Finite Field of size 2'
        """
        return "Monoid of ideals of %s" % self.__R

    def ring(self):
        """
        Return the ring of which this is the ideal monoid.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: O = K.maximal_order()
            sage: M = O.ideal_monoid(); M.ring() is O
            True
        """
        return self.__R

    def _element_constructor_(self, x):
        """
        Create an ideal in the monoid from ``x``.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: O = K.maximal_order()
            sage: M = O.ideal_monoid()
            sage: M(x)
            Ideal (x) of Maximal order of Rational function field in x over Finite Field of size 2
            sage: M([x-4, 1/x])
            Ideal (1/x) of Maximal order of Rational function field in x over Finite Field of size 2
        """
        try: # x is an ideal
            x = x.gens()
        except AttributeError:
            pass
        return self.__R.ideal(x)

    def _coerce_map_from_(self, x):
        """
        Used by coercion framework.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: O = K.maximal_order()
            sage: M = O.ideal_monoid()
            sage: M.has_coerce_map_from(O) # indirect doctest
            True
            sage: M.has_coerce_map_from(O.ideal_monoid())
            True
        """
        if isinstance(x, IdealMonoid):
            return self.ring().has_coerce_map_from(x.ring())
        else:
            return self.ring().has_coerce_map_from(x)

    def _an_element_(self):
        """
        Return an element of the ideal monoid.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: O = K.maximal_order()
            sage: M = O.ideal_monoid()
            sage: M.an_element() # indirect doctest; random
            Ideal (x) of Maximal order of Rational function field in x
            over Finite Field of size 2
        """
        x = self.__R.an_element()
        return self.__R.ideal([x])
