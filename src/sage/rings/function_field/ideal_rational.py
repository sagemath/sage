r"""
Ideals of function fields: rational
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
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.structure.richcmp import richcmp
from sage.rings.infinity import infinity

from .ideal import FunctionFieldIdeal, FunctionFieldIdealInfinite


class FunctionFieldIdeal_rational(FunctionFieldIdeal):
    """
    Fractional ideals of the maximal order of a rational function field.

    INPUT:

    - ``ring`` -- the maximal order of the rational function field

    - ``gen`` -- generator of the ideal, an element of the function field

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ)
        sage: O = K.maximal_order()
        sage: I = O.ideal(1/(x^2+x)); I
        Ideal (1/(x^2 + x)) of Maximal order of Rational function field in x over Rational Field
    """
    def __init__(self, ring, gen):
        """
        Initialize.

        TESTS::

            sage: K.<x> = FunctionField(QQ)
            sage: O = K.maximal_order()
            sage: I = O.ideal(1/(x^2+x))
            sage: TestSuite(I).run()
        """
        FunctionFieldIdeal.__init__(self, ring)
        self._gen = gen

    def __hash__(self):
        """
        Return the hash computed from the data.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: O = K.maximal_order()
            sage: I = O.ideal(1/(x^2+x))
            sage: d = { I: 1, I^2: 2 }
        """
        return hash( (self._ring, self._gen) )

    def __contains__(self, element):
        """
        Test if ``element`` is in this ideal.

        INPUT:

        - ``element`` -- element of the function field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: O = K.maximal_order()
            sage: I = O.ideal(1/(x+1))
            sage: x in I
            True
        """
        return (element / self._gen) in self._ring

    def _richcmp_(self, other, op):
        """
        Compare the element with the other element with respect
        to the comparison operator.

        INPUT:

        - ``other`` -- element

        - ``op`` -- comparison operator

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: O = K.maximal_order()
            sage: I = O.ideal(x, x^2 + 1)
            sage: J = O.ideal(x^2 + x + 1, x)
            sage: I == J
            True
            sage: I = O.ideal(x)
            sage: J = O.ideal(x + 1)
            sage: I < J
            True
        """
        return richcmp(self._gen, other._gen, op)

    def _add_(self, other):
        """
        Add this ideal with the ``other`` ideal.

        INPUT:

        - ``other`` -- ideal

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: O = K.maximal_order()
            sage: I = O.ideal(x, x^2 + 1)
            sage: J = O.ideal(x^2 + x + 1, x)
            sage: I + J == J + I
            True
        """
        return self._ring.ideal([self._gen, other._gen])

    def _mul_(self, other):
        """
        Multiply this ideal with the ``other`` ideal.

        INPUT:

        - ``other`` -- ideal

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: O = K.maximal_order()
            sage: I = O.ideal(x, x^2 + x)
            sage: J = O.ideal(x^2, x)
            sage: I * J == J * I
            True
        """
        return self._ring.ideal([self._gen * other._gen])

    def _acted_upon_(self, other, on_left):
        """
        Multiply ``other`` with this ideal on the right.

        INPUT:

        - ``other`` -- ideal

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: O = K.maximal_order()
            sage: I = O.ideal(x^3 + x^2)
            sage: 2 * I
            Ideal (x^3 + x^2) of Maximal order of Rational function field in x over Rational Field
            sage: x * I
            Ideal (x^4 + x^3) of Maximal order of Rational function field in x over Rational Field
        """
        return self._ring.ideal([other * self._gen])

    def __invert__(self):
        """
        Return the ideal inverse of this fractional ideal.

        EXAMPLES::

            sage: F.<x> = FunctionField(QQ)
            sage: O = F.maximal_order()
            sage: I = O.ideal(x/(x^2+1))
            sage: ~I
            Ideal ((x^2 + 1)/x) of Maximal order of Rational function field
            in x over Rational Field
        """
        return self._ring.ideal([~(self._gen)])

    def denominator(self):
        """
        Return the denominator of this fractional ideal.

        EXAMPLES::

            sage: F.<x> = FunctionField(QQ)
            sage: O = F.maximal_order()
            sage: I = O.ideal(x/(x^2+1))
            sage: I.denominator()
            x^2 + 1
        """
        return self._gen.denominator()

    def is_prime(self):
        """
        Return ``True`` if this is a prime ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: O = K.maximal_order()
            sage: I = O.ideal(x^3 + x^2)
            sage: [f.is_prime() for f,m in I.factor()]                                  # needs sage.libs.pari
            [True, True]
        """
        return self._gen.denominator() == 1 and self._gen.numerator().is_prime()

    @cached_method
    def module(self):
        """
        Return the module, that is the ideal viewed as a module over the ring.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: O = K.maximal_order()
            sage: I = O.ideal(x^3 + x^2)
            sage: I.module()                                                                                            # needs sage.modules
            Free module of degree 1 and rank 1 over Maximal order of Rational
            function field in x over Rational Field
            Echelon basis matrix:
            [x^3 + x^2]
            sage: J = 0*I
            sage: J.module()                                                                                            # needs sage.modules
            Free module of degree 1 and rank 0 over Maximal order of Rational
            function field in x over Rational Field
            Echelon basis matrix:
            []
        """
        V, fr, to = self.ring().fraction_field().vector_space()
        return V.span([to(g) for g in self.gens()], base_ring=self.ring())

    def gen(self):
        """
        Return the unique generator of this ideal.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(4))
            sage: O = K.maximal_order()
            sage: I = O.ideal(x^2 + x)
            sage: I.gen()
            x^2 + x
        """
        return self._gen

    def gens(self):
        """
        Return the tuple of the unique generator of this ideal.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(4))
            sage: O = K.maximal_order()
            sage: I = O.ideal(x^2 + x)
            sage: I.gens()
            (x^2 + x,)
        """
        return (self._gen,)

    def gens_over_base(self):
        """
        Return the generator of this ideal as a rank one module over the maximal
        order.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(4))
            sage: O = K.maximal_order()
            sage: I = O.ideal(x^2 + x)
            sage: I.gens_over_base()
            (x^2 + x,)
        """
        return (self._gen,)

    def valuation(self, ideal):
        """
        Return the valuation of the ideal at this prime ideal.

        INPUT:

        - ``ideal`` -- fractional ideal

        EXAMPLES::

            sage: F.<x> = FunctionField(QQ)
            sage: O = F.maximal_order()
            sage: I = O.ideal(x^2*(x^2+x+1)^3)
            sage: [f.valuation(I) for f,_ in I.factor()]                                # needs sage.libs.pari
            [2, 3]
        """
        if not self.is_prime():
            raise TypeError("not a prime ideal")

        O = self.ring()
        d = ideal.denominator()
        return self._valuation(d*ideal) - self._valuation(O.ideal(d))

    def _valuation(self, ideal):
        """
        Return the valuation of the integral ideal at this prime ideal.

        INPUT:

        - ``ideal`` -- ideal

        EXAMPLES::

            sage: F.<x> = FunctionField(QQ)
            sage: O = F.maximal_order()
            sage: p = O.ideal(x)
            sage: p.valuation(O.ideal(x + 1))  # indirect doctest                       # needs sage.libs.pari
            0
            sage: p.valuation(O.ideal(x^2))  # indirect doctest                         # needs sage.libs.pari
            2
            sage: p.valuation(O.ideal(1/x^3))  # indirect doctest                       # needs sage.libs.pari
            -3
            sage: p.valuation(O.ideal(0))  # indirect doctest                           # needs sage.libs.pari
            +Infinity
        """
        return ideal.gen().valuation(self.gen())

    def _factor(self):
        """
        Return the list of prime and multiplicity pairs of the
        factorization of this ideal.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(4))
            sage: O = K.maximal_order()
            sage: I = O.ideal(x^3*(x+1)^2)
            sage: I.factor()  # indirect doctest
            (Ideal (x) of Maximal order of Rational function field in x
            over Finite Field in z2 of size 2^2)^3 *
            (Ideal (x + 1) of Maximal order of Rational function field in x
            over Finite Field in z2 of size 2^2)^2
        """
        return [(self.ring().ideal(f), m) for f, m in self._gen.factor()]


class FunctionFieldIdealInfinite_rational(FunctionFieldIdealInfinite):
    """
    Fractional ideal of the maximal order of rational function field.

    INPUT:

    - ``ring`` -- infinite maximal order

    - ``gen`` -- generator

    Note that the infinite maximal order is a principal ideal domain.

    EXAMPLES::

        sage: K.<x> = FunctionField(GF(2))
        sage: Oinf = K.maximal_order_infinite()
        sage: Oinf.ideal(x)
        Ideal (x) of Maximal infinite order of Rational function field in x over Finite Field of size 2
    """
    def __init__(self, ring, gen):
        """
        Initialize.

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2))
            sage: Oinf = K.maximal_order_infinite()
            sage: I = Oinf.ideal(x)
            sage: TestSuite(I).run()
        """
        FunctionFieldIdealInfinite.__init__(self, ring)
        self._gen = gen

    def __hash__(self):
        """
        Return the hash of this fractional ideal.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2))
            sage: Oinf = K.maximal_order_infinite()
            sage: I = Oinf.ideal(x)
            sage: J = Oinf.ideal(1/x)
            sage: d = { I: 1, J: 2 }
        """
        return hash( (self.ring(), self._gen) )

    def __contains__(self, element):
        """
        Test if ``element`` is in this ideal.

        INPUT:

        - ``element`` -- element of the function field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: O = K.maximal_order_infinite()
            sage: I = O.ideal(1/(x+1))
            sage: x in I
            False
            sage: 1/x in I
            True
            sage: x/(x+1) in I
            False
            sage: 1/(x*(x+1)) in I
            True
        """
        return (element / self._gen) in self._ring

    def _richcmp_(self, other, op):
        """
        Compare this ideal and ``other`` with respect to ``op``.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2))
            sage: Oinf = K.maximal_order_infinite()
            sage: I = Oinf.ideal(x + 1)
            sage: J = Oinf.ideal(x^2 + x)
            sage: I + J == J
            True
        """
        return richcmp(self._gen, other._gen, op)

    def _add_(self, other):
        """
        Add this ideal with the other ideal.

        INPUT:

        - ``other`` -- ideal

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2))
            sage: Oinf = K.maximal_order_infinite()
            sage: I = Oinf.ideal(x/(x^2+1))
            sage: J = Oinf.ideal(1/(x+1))
            sage: I + J
            Ideal (1/x) of Maximal infinite order of Rational function field
            in x over Finite Field of size 2
        """
        return self._ring.ideal([self._gen, other._gen])

    def _mul_(self, other):
        """
        Multiply this ideal with the ``other`` ideal.

        INPUT:

        - ``other`` -- ideal

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2))
            sage: Oinf = K.maximal_order_infinite()
            sage: I = Oinf.ideal(x/(x^2+1))
            sage: J = Oinf.ideal(1/(x+1))
            sage: I * J
            Ideal (1/x^2) of Maximal infinite order of Rational function field
            in x over Finite Field of size 2
        """
        return self._ring.ideal([self._gen * other._gen])

    def _acted_upon_(self, other, on_left):
        """
        Multiply this ideal with the ``other`` ideal.

        INPUT:

        - ``other`` -- ideal

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2))
            sage: Oinf = K.maximal_order_infinite()
            sage: I = Oinf.ideal(x/(x^2+1))
            sage: x * I
            Ideal (1) of Maximal infinite order of Rational function field
            in x over Finite Field of size 2
        """
        return self._ring.ideal([other * self._gen])

    def __invert__(self):
        """
        Return the multiplicative inverse of this ideal.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2))
            sage: Oinf = K.maximal_order_infinite()
            sage: I = Oinf.ideal(x/(x^2 + 1))
            sage: ~I  # indirect doctest
            Ideal (x) of Maximal infinite order of Rational function field in x
            over Finite Field of size 2
        """
        return self._ring.ideal([~self._gen])

    def is_prime(self):
        """
        Return ``True`` if this ideal is a prime ideal.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2))
            sage: Oinf = K.maximal_order_infinite()
            sage: I = Oinf.ideal(x/(x^2 + 1))
            sage: I.is_prime()
            True
        """
        x = self._ring.fraction_field().gen()
        return self._gen == 1/x

    def gen(self):
        """
        Return the generator of this principal ideal.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2))
            sage: Oinf = K.maximal_order_infinite()
            sage: I = Oinf.ideal((x+1)/(x^3+x), (x^2+1)/x^4)
            sage: I.gen()
            1/x^2
        """
        return self._gen

    def gens(self):
        """
        Return the generator of this principal ideal.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2))
            sage: Oinf = K.maximal_order_infinite()
            sage: I = Oinf.ideal((x+1)/(x^3+x), (x^2+1)/x^4)
            sage: I.gens()
            (1/x^2,)
        """
        return (self._gen,)

    def gens_over_base(self):
        """
        Return the generator of this ideal as a rank one module
        over the infinite maximal order.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2))
            sage: Oinf = K.maximal_order_infinite()
            sage: I = Oinf.ideal((x+1)/(x^3+x), (x^2+1)/x^4)
            sage: I.gens_over_base()
            (1/x^2,)
        """
        return (self._gen,)

    def valuation(self, ideal):
        """
        Return the valuation of ``ideal`` at this prime ideal.

        INPUT:

        - ``ideal`` -- fractional ideal

        EXAMPLES::

            sage: F.<x> = FunctionField(QQ)
            sage: O = F.maximal_order_infinite()
            sage: p = O.ideal(1/x)
            sage: p.valuation(O.ideal(x/(x+1)))
            0
            sage: p.valuation(O.ideal(0))
            +Infinity
        """
        if not self.is_prime():
            raise TypeError("not a prime ideal")

        f = ideal.gen()
        if f == 0:
            return infinity
        else:
            return f.denominator().degree() - f.numerator().degree()

    def _factor(self):
        """
        Return the factorization of this ideal into prime ideals.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(2))
            sage: Oinf = K.maximal_order_infinite()
            sage: I = Oinf.ideal((x+1)/(x^3+1))
            sage: I._factor()
            [(Ideal (1/x) of Maximal infinite order of Rational function field in x
            over Finite Field of size 2, 2)]
        """
        g = ~(self.ring().fraction_field().gen())
        m = self._gen.denominator().degree() - self._gen.numerator().degree()
        if m == 0:
            return []
        else:
            return [(self.ring().ideal(g), m)]
