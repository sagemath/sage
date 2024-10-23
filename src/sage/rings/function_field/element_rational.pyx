r"""
Elements of function fields: rational
"""

# *****************************************************************************
#       Copyright (C) 2010      William Stein <wstein@gmail.com>
#                     2010      Robert Bradshaw <robertwb@math.washington.edu>
#                     2011-2020 Julian Rueth <julian.rueth@gmail.com>
#                     2011      Maarten Derickx <m.derickx.student@gmail.com>
#                     2015      Nils Bruin
#                     2016      Frédéric Chapoton
#                     2017-2019 Kwankyu Lee
#                     2018-2020 Travis Scrimshaw
#                     2019      Brent Baccala
#                     2021      Saher Amasha
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.structure.richcmp cimport richcmp, richcmp_not_equal
from sage.structure.element cimport FieldElement

from sage.rings.function_field.element cimport FunctionFieldElement


cdef class FunctionFieldElement_rational(FunctionFieldElement):
    """
    Elements of a rational function field.

    EXAMPLES::

        sage: K.<t> = FunctionField(QQ); K
        Rational function field in t over Rational Field
        sage: t^2 + 3/2*t
        t^2 + 3/2*t
        sage: FunctionField(QQ,'t').gen()^3
        t^3
    """
    def __init__(self, parent, x, reduce=True):
        """
        Initialize.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: x = t^3
            sage: TestSuite(x).run()
        """
        FieldElement.__init__(self, parent)
        self._x = x

    def __pari__(self):
        r"""
        Coerce the element to PARI.

        EXAMPLES::

            sage: K.<a> = FunctionField(QQ)
            sage: ((a+1)/(a-1)).__pari__()                                              # needs sage.libs.pari
            (a + 1)/(a - 1)
        """
        return self.element().__pari__()

    def element(self):
        """
        Return the underlying fraction field element that represents the element.

        EXAMPLES::

            sage: K.<t> = FunctionField(GF(7))
            sage: t.element()
            t
            sage: type(t.element())                                                     # needs sage.libs.ntl
            <... 'sage.rings.fraction_field_FpT.FpTElement'>

            sage: # needs sage.rings.finite_rings
            sage: K.<t> = FunctionField(GF(131101))
            sage: t.element()
            t
            sage: type(t.element())
            <... 'sage.rings.fraction_field_element.FractionFieldElement_1poly_field'>
        """
        return self._x

    cpdef list list(self):
        """
        Return a list with just the element.

        The list represents the element when the rational function field is
        viewed as a (one-dimensional) vector space over itself.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: t.list()
            [t]
        """
        return [self]

    def _repr_(self):
        """
        Return the string representation of the element.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: t._repr_()
            't'
        """
        return repr(self._x)

    def __bool__(self):
        """
        Return ``True`` if the element is not zero.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: bool(t)
            True
            sage: bool(K(0))
            False
            sage: bool(K(1))
            True
        """
        return not not self._x

    def __hash__(self):
        """
        Return the hash of the element.

        TESTS:

        It would be nice if the following would produce a list of
        15 distinct hashes::

            sage: K.<t> = FunctionField(QQ)
            sage: len({hash(t^i+t^j) for i in [-2..2] for j in [i..2]}) >= 10
            True
        """
        return hash(self._x)

    cpdef _richcmp_(self, other, int op):
        """
        Compare the element with the other element with respect to ``op``.

        INPUT:

        - ``other`` -- element

        - ``op`` -- comparison operator

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: t > 0
            True
            sage: t < t^2
            True
        """
        cdef FunctionFieldElement left
        cdef FunctionFieldElement right
        try:
            left = <FunctionFieldElement?>self
            right = <FunctionFieldElement?>other
            lp = left._parent
            rp = right._parent
            if lp != rp:
                return richcmp_not_equal(lp, rp, op)
            return richcmp(left._x, right._x, op)
        except TypeError:
            return NotImplemented

    cpdef _add_(self, right):
        """
        Add the element with the other element.

        INPUT:

        - ``right`` -- element

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: t + (3*t^3)                      # indirect doctest
            3*t^3 + t
        """
        cdef FunctionFieldElement res = self._new_c()
        res._x = self._x + (<FunctionFieldElement>right)._x
        return res

    cpdef _sub_(self, right):
        """
        Subtract the other element from the element.

        INPUT:

        - ``right`` -- element

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: t - (3*t^3)                      # indirect doctest
            -3*t^3 + t
        """
        cdef FunctionFieldElement res = self._new_c()
        res._x = self._x - (<FunctionFieldElement>right)._x
        return res

    cpdef _mul_(self, right):
        """
        Multiply the element with the other element.

        INPUT:

        - ``right`` -- element

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: (t+1) * (t^2-1)                  # indirect doctest
            t^3 + t^2 - t - 1
        """
        cdef FunctionFieldElement res = self._new_c()
        res._x = self._x * (<FunctionFieldElement>right)._x
        return res

    cpdef _div_(self, right):
        """
        Divide the element with the other element.

        INPUT:

        - ``right`` -- element

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: (t+1) / (t^2 - 1)                # indirect doctest
            1/(t - 1)
        """
        cdef FunctionFieldElement res = self._new_c()
        res._parent = self._parent.fraction_field()
        res._x = self._x / (<FunctionFieldElement>right)._x
        return res

    def numerator(self):
        """
        Return the numerator of the rational function.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: f = (t+1) / (t^2 - 1/3); f
            (t + 1)/(t^2 - 1/3)
            sage: f.numerator()
            t + 1
        """
        return self._x.numerator()

    def denominator(self):
        """
        Return the denominator of the rational function.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: f = (t+1) / (t^2 - 1/3); f
            (t + 1)/(t^2 - 1/3)
            sage: f.denominator()
            t^2 - 1/3
        """
        return self._x.denominator()

    def valuation(self, place):
        """
        Return the valuation of the rational function at the place.

        Rational function field places are associated with irreducible
        polynomials.

        INPUT:

        - ``place`` -- a place or an irreducible polynomial

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: f = (t - 1)^2*(t + 1)/(t^2 - 1/3)^3
            sage: f.valuation(t - 1)
            2
            sage: f.valuation(t)
            0
            sage: f.valuation(t^2 - 1/3)
            -3

            sage: K.<x> = FunctionField(GF(2))
            sage: p = K.places_finite()[0]                                              # needs sage.libs.pari
            sage: (1/x^2).valuation(p)                                                  # needs sage.libs.pari
            -2
        """
        from sage.rings.function_field.place import FunctionFieldPlace

        if not isinstance(place, FunctionFieldPlace):
            # place is an irreducible polynomial
            R = self._parent._ring
            return self._x.valuation(R(self._parent(place)._x))

        prime = place.prime_ideal()
        ideal = prime.ring().ideal(self)
        return prime.valuation(ideal)

    def is_square(self):
        """
        Return whether the element is a square.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: t.is_square()
            False
            sage: (t^2/4).is_square()
            True
            sage: f = 9 * (t+1)^6 / (t^2 - 2*t + 1); f.is_square()
            True

            sage: K.<t> = FunctionField(GF(5))
            sage: (-t^2).is_square()                                                    # needs sage.libs.pari
            True
            sage: (-t^2).sqrt()                                                         # needs sage.libs.pari
            2*t
        """
        return self._x.is_square()

    def sqrt(self, all=False):
        """
        Return the square root of the rational function.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ)
            sage: f = t^2 - 2 + 1/t^2; f.sqrt()
            (t^2 - 1)/t
            sage: f = t^2; f.sqrt(all=True)
            [t, -t]

        TESTS::

            sage: K(4/9).sqrt()
            2/3
            sage: K(0).sqrt(all=True)
            [0]
        """
        if all:
            return [self._parent(r) for r in self._x.sqrt(all=True)]
        else:
            return self._parent(self._x.sqrt())

    cpdef bint is_nth_power(self, n) noexcept:
        r"""
        Return whether this element is an ``n``-th power in the rational
        function field.

        INPUT:

        - ``n`` -- integer

        OUTPUT:

        Returns ``True`` if there is an element `a` in the function field such
        that this element equals `a^n`.

        ALGORITHM:

        If ``n`` is a power of the characteristic of the field and the constant
        base field is perfect, then this uses the algorithm described in Lemma
        3 of [GiTr1996]_.

        .. SEEALSO::

            :meth:`nth_root`

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<x> = FunctionField(GF(3))
            sage: f = (x+1)/(x-1)
            sage: f.is_nth_power(1)
            True
            sage: f.is_nth_power(3)                                                     # needs sage.modules
            False
            sage: (f^3).is_nth_power(3)                                                 # needs sage.modules
            True
            sage: (f^9).is_nth_power(-9)                                                # needs sage.modules
            True
        """
        if n == 1:
            return True
        if n < 0:
            return (~self).is_nth_power(-n)

        p = self._parent.characteristic()
        if n == p:
            return self._parent.derivation()(self).is_zero()
        if p.divides(n):
            return self.is_nth_power(p) and self.nth_root(p).is_nth_power(n//p)
        if n == 2:
            return self.is_square()

        raise NotImplementedError("is_nth_power() not implemented for the given n")

    cpdef FunctionFieldElement nth_root(self, n):
        r"""
        Return an ``n``-th root of this element in the function field.

        INPUT:

        - ``n`` -- integer

        OUTPUT:

        Returns an element ``a`` in the rational function field such that this
        element equals `a^n`. Raises an error if no such element exists.

        ALGORITHM:

        If ``n`` is a power of the characteristic of the field and the constant
        base field is perfect, then this uses the algorithm described in
        Corollary 3 of [GiTr1996]_.

        .. SEEALSO::

            :meth:`is_nth_power`

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(3))
            sage: f = (x+1)/(x+2)
            sage: f.nth_root(1)
            (x + 1)/(x + 2)
            sage: f.nth_root(3)
            Traceback (most recent call last):
            ...
            ValueError: element is not an n-th power
            sage: (f^3).nth_root(3)                                                     # needs sage.modules
            (x + 1)/(x + 2)
            sage: (f^9).nth_root(-9)                                                    # needs sage.modules
            (x + 2)/(x + 1)
        """
        if n == 0:
            if not self.is_one():
                raise ValueError("element is not a 0-th power")
            return self
        if n == 1:
            return self
        if n < 0:
            return (~self).nth_root(-n)
        p = self._parent.characteristic()
        if p.divides(n):
            if not self.is_nth_power(p):
                raise ValueError("element is not an n-th power")
            return self._parent(self.numerator().nth_root(p) / self.denominator().nth_root(p)).nth_root(n//p)
        if n == 2:
            return self.sqrt()

        raise NotImplementedError("nth_root() not implemented for {}".format(n))

    def factor(self):
        """
        Factor the rational function.

        EXAMPLES::

            sage: # needs sage.libs.pari
            sage: K.<t> = FunctionField(QQ)
            sage: f = (t+1) / (t^2 - 1/3)
            sage: f.factor()
            (t + 1) * (t^2 - 1/3)^-1
            sage: (7*f).factor()
            (7) * (t + 1) * (t^2 - 1/3)^-1
            sage: ((7*f).factor()).unit()
            7
            sage: (f^3).factor()
            (t + 1)^3 * (t^2 - 1/3)^-3
        """
        P = self.parent()
        F = self._x.factor()
        from sage.structure.factorization import Factorization
        return Factorization([(P(a),e) for a,e in F], unit=F.unit())

    def inverse_mod(self, I):
        """
        Return an inverse of the element modulo the integral ideal `I`, if `I`
        and the element together generate the unit ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: O = K.maximal_order(); I = O.ideal(x^2 + 1)
            sage: t = O(x + 1).inverse_mod(I); t
            -1/2*x + 1/2
            sage: (t*(x+1) - 1) in I
            True
        """
        assert len(I.gens()) == 1
        f = I.gens()[0]._x
        assert f.denominator() == 1
        assert self._x.denominator() == 1
        return self.parent()(self._x.numerator().inverse_mod(f.numerator()))
