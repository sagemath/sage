r"""
Fraction field elements of Ore polynomial rings

AUTHOR:

- Xavier Caruso (2020-05)
"""

# ***************************************************************************
#    Copyright (C) 2020 Xavier Caruso <xavier.caruso@normalesup.org>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ***************************************************************************

from sage.structure.richcmp import richcmp, op_EQ, op_NE
from sage.misc.cachefunc import cached_method
from sage.misc.latex import latex

from sage.categories.map import Map
from sage.categories.morphism import Morphism
from sage.categories.homset import Hom
from sage.structure.element import AlgebraElement


class OreFunction(AlgebraElement):
    r"""
    An element in a Ore function field.
    """
    def __init__(self, parent, numerator, denominator=None, simplify=True):
        r"""
        Initialize this element.

        TESTS::

            sage: R.<t> = GF(5)[]
            sage: der = R.derivation()
            sage: A.<d> = R['d', der]
            sage: K = A.fraction_field()
            sage: f = K.random_element()

            sage: # TestSuite(f).run()
        """
        AlgebraElement.__init__(self, parent)
        ring = parent._ring
        numerator = ring(numerator)
        if denominator is None:
            denominator = ring.one()
        else:
            denominator = ring(denominator)
            if not denominator:
                raise ZeroDivisionError("denominator must be nonzero")
        # We normalize the fraction
        if numerator:
            if simplify and parent._simplification:
                D = numerator.left_gcd(denominator, monic=False)
                numerator, _ = numerator.left_quo_rem(D)
                denominator, _ = denominator.left_quo_rem(D)
            s = denominator.leading_coefficient()
            if s != 1:
                s = ~s
                numerator = s*numerator
                denominator = s*denominator
            self._numerator = numerator
            self._denominator = denominator
        else:
            self._numerator = ring.zero()
            self._denominator = ring.one()

    def _repr_(self):
        r"""
        Return a string representation of this element.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x', Frob]
            sage: f = (x+a)^(-1) * (x^2 + a^2)
            sage: f
            (x + a)^(-1) * (x^2 + a^2)

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: f = 1/x^3; f
            x^(-3)
            sage: f * x^5
            x^2
        """
        if not self._numerator:
            return "0"
        if self._denominator == 1:
            return str(self._numerator)
        if self._denominator.is_monomial():
            s = "%s^(-%s)" % (self.parent().variable_name(), self._denominator.degree())
        else:
            s = "(%s)^(-1)" % self._denominator
        if self._numerator == 1:
            return s
        if self._numerator._is_atomic():
            return "%s * %s" % (s, self._numerator)
        else:
            return "%s * (%s)" % (s, self._numerator)

    def _latex_(self):
        r"""
        Return a LaTeX representation of this element.

        EXAMPLES::

            sage: R.<t> = GF(5)[]
            sage: der = R.derivation()
            sage: S.<d> = R['d', der]
            sage: f = (d+t)^(-1) * (d^2 + t^2)
            sage: latex(f)
            \left(d + t\right)^{-1} \cdot \left(d^{2} + t^{2}\right)
        """
        if not self._numerator:
            return "0"
        if self._denominator == 1:
            return latex(self._numerator)
        if self._denominator.is_monomial():
            s = "%s^{-%s}" % (self.parent().latex_variable_names()[0], self._denominator.degree())
        else:
            s = "\\left(%s\\right)^{-1}" % self._denominator
        if self._numerator == 1:
            return s
        if self._numerator._is_atomic():
            return "%s \\cdot %s" % (s, latex(self._numerator))
        else:
            return "%s \\cdot \\left(%s\\right)" % (s, latex(self._numerator))

    def __hash__(self):
        r"""
        Return a hash of this element.

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x', Frob]
            sage: K = S.fraction_field()
            sage: f = K.random_element()
            sage: hash(f)  # random
            1700763101013238501
        """
        return hash((self._numerator, self._denominator))

    def _richcmp_(self, other, op):
        r"""
        Compare this element with ``other`` for the comparison
        operator ``op``.

        TESTS::

            sage: # needs sage.rings.function_field
            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: der = R.derivation(1, twist=sigma)
            sage: S.<delta> = R['delta', der]
            sage: K = S.fraction_field()
            sage: P = K.random_element()
            sage: Q = K.random_element()
            sage: D = K.random_element()
            sage: Q == 0 or D == 0 or (P*D) / (Q*D) == P/Q  # long time
            True
        """
        if self.parent()._simplification:
            return richcmp((self._numerator, self._denominator), (other._numerator, other._denominator), op)
        if op == op_EQ or op == op_NE:
            _, U, V = self._denominator.left_xlcm(other._denominator)
            return richcmp(U * self._numerator, V * other._numerator, op)
        return NotImplemented

    def left_denominator(self):
        r"""
        Return `s` if this element reads `s^{-1} t`.

        WARNING:

        When the twisting morphism is bijective, there is a unique
        irreducible fraction of the form `s^{-1} t` representing this
        element. Here irreducible means that `s` and `t` have no
        nontrivial common left divisor.
        Under this additional assumption, this method always returns
        this distinguished denominator `s`.

        On the contrary, when the twisting morphism is not bijective,
        this method returns the denominator of *some* fraction
        representing the input element.
        However, the software guarantees that the method :meth:`right_numerator`
        outputs the numerator of the *same* fraction.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x', Frob]
            sage: s = x + a
            sage: t = x^2 + a*x + a^2
            sage: f = s^(-1) * t
            sage: f.left_denominator()
            x + a

        In the example below, a simplification occurs::

            sage: # needs sage.rings.finite_rings
            sage: u = S.random_element(degree=2)
            sage: g = (u*s)^(-1) * (u*t)
            sage: g.left_denominator()
            x + a

        When the twisting morphism is not invertible, simplifications
        do not occur in general::

            sage: R.<z> = GF(11)[]
            sage: sigma = R.hom([z^2])
            sage: S.<x> = R['x', sigma]
            sage: s = (x + z)^2
            sage: t = (x + z) * (x^2 + z^2)
            sage: f = s^(-1) * t                                                        # needs sage.rings.function_field
            sage: f.left_denominator()                                                  # needs sage.rings.function_field
            x^2 + (z^2 + z)*x + z^2

        However, the following always holds true::

            sage: f == f.left_denominator()^(-1) * f.right_numerator()                  # needs sage.rings.function_field
            True

        .. SEEALSO::

            :meth:`right_numerator`, :meth:`left_numerator`, :meth:`right_denominator`
        """
        return self._denominator

    def right_numerator(self):
        r"""
        Return `t` if this element reads `s^{-1} t`.

        WARNING:

        When the twisting morphism is bijective, there is a unique
        irreducible fraction of the form `s^{-1} t` representing this
        element. Here irreducible means that `s` and `t` have no
        nontrivial common left divisor.
        Under this additional assumption, this method always returns
        this distinguished numerator `t`.

        On the contrary, when the twisting morphism is not bijective,
        this method returns the numerator of *some* fraction
        representing the input element.
        However, the software guarantees that the method :meth:`left_denominator`
        outputs the numerator of the *same* fraction.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x', Frob]
            sage: s = x + a
            sage: t = x^2 + a*x + a^2
            sage: f = s^(-1) * t
            sage: f.right_numerator()
            x^2 + a*x + a^2

        In the example below, a simplification occurs::

            sage: # needs sage.rings.finite_rings
            sage: u = S.random_element(degree=2)
            sage: g = (u*s)^(-1) * (u*t)
            sage: g.right_numerator()
            x^2 + a*x + a^2

        .. SEEALSO::

            :meth:`left_denominator`, :meth:`left_numerator`, :meth:`right_denominator`
        """
        return self._numerator

    @cached_method
    def _reverse_fraction(self):
        r"""
        Return the pair `(s,t)` if this element reads `t s^{-1}`.

        This is a helper function. Do not call it directly.

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(11^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: der = k.derivation(a+1, twist=Frob)
            sage: S.<x> = k['x', der]
            sage: P = S.random_element(degree=5)
            sage: Q = S.random_element(degree=5)
            sage: f = P / Q
            sage: f == f.left_numerator() / f.right_denominator()   # indirect doctest
            True
        """
        _, denominator, numerator = self._numerator.right_xlcm(self._denominator, monic=False)
        d = denominator.degree()
        s = ~(denominator.leading_coefficient())
        morphism = self.parent().twisting_morphism(-d)
        if morphism is not None:
            s = morphism(s)
        numerator = numerator * s
        denominator = denominator * s
        return numerator, denominator

    def right_denominator(self):
        r"""
        Return `s` if this element reads `t s^{-1}`.

        WARNING:

        When the twisting morphism is bijective, there is a unique
        irreducible fraction of the form `t s^{-1}` representing this
        element. Here irreducible means that `s` and `t` have no
        nontrivial common right divisor.
        Under this additional assumption, this method always returns
        this distinguished denominator `s`.

        On the contrary, when the twisting morphism is not bijective,
        the existence of the writing `t s^{-1}` is not guaranteed in
        general. In this case, this method raises an error.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x', Frob]
            sage: s = x + a
            sage: t = x^2 + a*x + a^2
            sage: f = t/s
            sage: f.right_denominator()
            x + a

        In the example below, a simplification occurs::

            sage: # needs sage.rings.finite_rings
            sage: u = S.random_element(degree=2)
            sage: g = (t*u) / (s*u)
            sage: g.right_denominator()
            x + a

        .. SEEALSO::

            :meth:`left_numerator`, :meth:`left_denominator`, :meth:`right_numerator`

        TESTS::

            sage: R.<z> = GF(11)[]
            sage: sigma = R.hom([z^2])
            sage: S.<x> = R['x', sigma]
            sage: f = (x + z) / (x - z)                                                 # needs sage.rings.function_field
            sage: f.right_denominator()                                                 # needs sage.rings.function_field
            Traceback (most recent call last):
            ...
            NotImplementedError: inversion of the twisting morphism Ring endomorphism
            of Fraction Field of Univariate Polynomial Ring in z over Finite Field of size 11
              Defn: z |--> z^2
        """
        return self._reverse_fraction()[1]

    def left_numerator(self):
        r"""
        Return `t` if this element reads `t s^{-1}`.

        WARNING:

        When the twisting morphism is bijective, there is a unique
        irreducible fraction of the form `t s^{-1}` representing this
        element. Here irreducible means that `s` and `t` have no
        nontrivial common right divisor.
        Under this additional assumption, this method always returns
        this distinguished numerator `t`.

        On the contrary, when the twisting morphism is not bijective,
        the existence of the writing `t s^{-1}` is not guaranteed in
        general. In this case, this method raises an error.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x', Frob]
            sage: s = x + a
            sage: t = x^2 + a*x + a^2
            sage: f = t/s
            sage: f.left_numerator()
            x^2 + a*x + a^2

        In the example below, a simplification occurs::

            sage: # needs sage.rings.finite_rings
            sage: u = S.random_element(degree=2)
            sage: g = (t*u) / (s*u)
            sage: g.left_numerator()
            x^2 + a*x + a^2
        """
        return self._reverse_fraction()[0]

    def is_zero(self):
        r"""
        Return ``True`` if this element is equal to zero.

        EXAMPLES::

            sage: R.<t> = GF(3)[]
            sage: der = R.derivation()
            sage: A.<d> = R['x', der]
            sage: f = t/d
            sage: f.is_zero()
            False
            sage: (f-f).is_zero()
            True
        """
        return self._numerator.is_zero()

    def _add_(self, other):
        r"""
        Return the sum of this element and ``other``.

        INPUT:

        - ``other`` -- a Ore function

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(5^2)
            sage: Frob = k.frobenius_endomorphism()
            sage: der = k.derivation(a, twist=Frob)
            sage: S.<x> = k['x', der]
            sage: K = S.fraction_field()
            sage: f = K.random_element()
            sage: g = K.random_element()
            sage: h = K.random_element()
            sage: f + g == g + f
            True
            sage: f + (g + h) == (f + g) + h
            True
        """
        denominator, U, V = self._denominator.left_xlcm(other._denominator)
        numerator = U * self._numerator + V * other._numerator
        return self.parent()(numerator, denominator)

    def _sub_(self, other):
        r"""
        Return the subtraction of this element and ``other``.

        INPUT:

        - ``other`` -- a Ore function

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(5^2)
            sage: Frob = k.frobenius_endomorphism()
            sage: der = k.derivation(a, twist=Frob)
            sage: S.<x> = k['x', der]
            sage: K = S.fraction_field()
            sage: f = K.random_element()
            sage: g = K.random_element()
            sage: h = K.random_element()
            sage: f - (g - h) == (f - g) + h
            True
        """
        denominator, U, V = self._denominator.left_xlcm(other._denominator)
        numerator = U * self._numerator - V * other._numerator
        return self.parent()(numerator, denominator)

    def _neg_(self):
        r"""
        Return the opposite of this element.

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(5^2)
            sage: Frob = k.frobenius_endomorphism()
            sage: der = k.derivation(a, twist=Frob)
            sage: S.<x> = k['x', der]
            sage: K = S.fraction_field()
            sage: f = K.random_element()
            sage: g = -f
            sage: (f+g).is_zero()
            True
        """
        return self.parent()(-self._numerator, self._denominator, simplify=False)

    def _mul_(self, other):
        r"""
        Return the product of this element and ``other``.

        INPUT:

        - ``other`` -- a Ore function

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(5^2)
            sage: Frob = k.frobenius_endomorphism()
            sage: der = k.derivation(a, twist=Frob)
            sage: S.<x> = k['x', der]
            sage: K = S.fraction_field()
            sage: f = K.random_element()
            sage: g = K.random_element()
            sage: h = K.random_element()
            sage: (f * g) * h == f * (g * h)
            True
            sage: f * (g + h) == f*g + f*h
            True
            sage: (f + g) * h == f*h + g*h
            True
        """
        if self.is_zero():
            return self
        if other.is_zero():
            return other
        L, U, V = self._numerator.left_xlcm(other._denominator, monic=False)
        denominator = U * self._denominator
        numerator = V * other._numerator
        return self.parent()(numerator, denominator)

    def _div_(self, other):
        r"""
        Return the division of this element by ``other``.

        INPUT:

        - ``other`` -- a Ore function

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(5^2)
            sage: Frob = k.frobenius_endomorphism()
            sage: der = k.derivation(a, twist=Frob)
            sage: S.<x> = k['x', der]
            sage: K = S.fraction_field()
            sage: f = K.random_element()
            sage: g = K.random_element()
            sage: h = K.random_element()
            sage: g == 0 or h == 0 or f / (g / h) == f*h / g
            True
            sage: 0/f
            0
            sage: f/0
            Traceback (most recent call last):
            ...
            ZeroDivisionError: cannot divide by zero

        We check that :issue:`32109` is fixed::

            sage: K(0)/K(0)                                                             # needs sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            ZeroDivisionError: cannot divide by zero
        """
        if not other._numerator:
            raise ZeroDivisionError("cannot divide by zero")
        if not self._numerator:
            return self
        L, U, V = self._numerator.left_xlcm(other._numerator, monic=False)
        denominator = U * self._denominator
        numerator = V * other._denominator
        return self.parent()(numerator, denominator)

    def __invert__(self):
        r"""
        Return the inverse of this element.

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(5^2)
            sage: Frob = k.frobenius_endomorphism()
            sage: der = k.derivation(a, twist=Frob)
            sage: S.<x> = k['x', der]
            sage: K = S.fraction_field()
            sage: f = K.random_element()
            sage: g = ~f
            sage: f * g
            1
            sage: ~K(0)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: cannot divide by zero
        """
        if not self._numerator:
            raise ZeroDivisionError("cannot divide by zero")
        return self.parent()(self._denominator, self._numerator)

    def hilbert_shift(self, s, var=None):
        r"""
        Return this Ore function with variable shifted by `s`,
        i.e. if this Ore function is `f(x)`, return `f(x+s)`.

        INPUT:

        - ``s`` -- an element in the base ring

        - ``var`` -- string; the variable name

        EXAMPLES::

            sage: R.<t> = GF(7)[]
            sage: der = R.derivation()
            sage: A.<d> = R['d', der]
            sage: K = A.fraction_field()

            sage: f = 1 / (d-t)
            sage: f.hilbert_shift(t)
            d^(-1)

        One can specify another variable name::

            sage: f.hilbert_shift(t, var='x')
            x^(-1)

        When the twisting morphism is not trivial, the output lies
        in a different Ore polynomial ring::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x', Frob]
            sage: K = S.fraction_field()
            sage: f = (x-a)^(-2)
            sage: g = f.hilbert_shift(a); g
            x^(-2)
            sage: g.parent()
            Ore Function Field in x over Finite Field in a of size 5^3
             twisted by a |--> a^5 and a*([a |--> a^5] - id)
            sage: g.parent() is S
            False

        This behavior ensures that the Hilbert shift by a fixed element
        defines a homomorphism of fields::

            sage: # needs sage.rings.finite_rings
            sage: U = K.random_element(degree=5)
            sage: V = K.random_element(degree=5)
            sage: s = k.random_element()
            sage: (U+V).hilbert_shift(s) == U.hilbert_shift(s) + V.hilbert_shift(s)
            True
            sage: (U*V).hilbert_shift(s) == U.hilbert_shift(s) * V.hilbert_shift(s)
            True
        """
        numerator = self._numerator.hilbert_shift(s, var)
        denominator = self._denominator.hilbert_shift(s, var)
        parent = numerator.parent().fraction_field()
        return parent(numerator, denominator, simplify=False)


class ConstantOreFunctionSection(Map):
    r"""
    Representation of the canonical homomorphism from the constants of a Ore
    function field to the base field.

    This class is needed by the coercion system.

    EXAMPLES::

        sage: # needs sage.rings.finite_rings
        sage: from sage.rings.polynomial.ore_polynomial_element import ConstantOrePolynomialSection
        sage: k.<a> = GF(5^3)
        sage: Frob = k.frobenius_endomorphism()
        sage: S.<x> = k['x', Frob]
        sage: K = S.fraction_field()
        sage: iota = K.coerce_map_from(k)
        sage: sigma = iota.section(); sigma
        Generic map:
          From: Ore Function Field in x over Finite Field in a of size 5^3
                twisted by a |--> a^5
          To:   Finite Field in a of size 5^3
    """
    def _call_(self, x):
        r"""
        Return `x` viewed in the base field,
        or raise an error if `x` is not a constant Ore function.

        TESTS::

            sage: R.<t> = QQ[]
            sage: F = R.fraction_field()
            sage: sigma = R.hom([t^2])
            sage: S.<x> = R['x', sigma]
            sage: P = S._random_nonzero_element()
            sage: f = (t*P) / P                                                         # needs sage.rings.function_field
            sage: F(f)                                                                  # needs sage.rings.function_field
            t
            sage: g = x / (x+t)                                                         # needs sage.rings.function_field
            sage: F(g)                                                                  # needs sage.rings.function_field
            Traceback (most recent call last):
            ...
            TypeError: (x + t^2)^(-1) * x is not a constant function
        """
        numerator = x._numerator
        denominator = x._denominator
        if numerator.degree() == denominator.degree() and denominator.right_divides(numerator):
            return numerator.leading_coefficient() / denominator.leading_coefficient()
        raise TypeError(f"{x} is not a constant function")


class OreFunctionBaseringInjection(Morphism):
    r"""
    Representation of the canonical homomorphism from a field `k` into a Ore
    function field over `k`.

    This class is needed by the coercion system.
    """
    def __init__(self, domain, codomain):
        r"""
        Initialize this morphism.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: K = S.fraction_field()                                                # needs sage.rings.function_field
            sage: K.coerce_map_from(K.base_ring())  # indirect doctest                  # needs sage.rings.function_field
            Ore Function base injection morphism:
              From: Fraction Field of Univariate Polynomial Ring in t over Rational Field
              To:   Ore Function Field in x over Fraction Field of Univariate Polynomial Ring in t over Rational Field twisted by t |--> t + 1
        """
        assert codomain.base_ring() is domain, \
            "the domain of the injection must be the base ring of the Ore function field"
        Morphism.__init__(self, Hom(domain,codomain))
        self._an_element = codomain.gen()
        self._repr_type_str = "Ore Function base injection"

    def an_element(self):
        r"""
        Return an element of the codomain of the ring homomorphism.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: K = S.fraction_field()
            sage: m = K.coerce_map_from(k)
            sage: m.an_element()
            x
        """
        return self._an_element

    def _call_(self, x):
        r"""
        Return the constant Ore function equal to `x`.

        INPUT:

        - ``e`` -- element in the base ring

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: K = S.fraction_field()
            sage: m = K.coerce_map_from(k)
            sage: m(t)
            t
            sage: m(t).parent()
            Ore Function Field in x over Finite Field in t of size 5^3 twisted by t |--> t^5
        """
        codomain = self.codomain()
        try:
            return codomain._element_constructor_(x)
        except AttributeError:
            return codomain(x)

    def section(self):
        r"""
        Return the canonical homomorphism from the constants of a Ore
        function filed to its base field.

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: K = S.fraction_field()
            sage: m = K.coerce_map_from(k)
            sage: m.section()
            Generic map:
              From: Ore Function Field in x over Finite Field in t of size 5^3
                    twisted by t |--> t^5
              To:   Finite Field in t of size 5^3
        """
        return ConstantOreFunctionSection(self.codomain(), self.domain())


# Ore functions over Ore function field with finite index center
################################################################

class OreFunction_with_large_center(OreFunction):
    r"""
    A special class for elements of Ore function fields whose
    center has finite index.

    TESTS::

        sage: # needs sage.rings.finite_rings
        sage: k.<a> = GF(5^3)
        sage: Frob = k.frobenius_endomorphism()
        sage: S.<x> = k['x', Frob]
        sage: K = S.fraction_field()
        sage: f = K.random_element()
        sage: from sage.rings.polynomial.ore_function_element import OreFunction_with_large_center
        sage: isinstance(f, OreFunction_with_large_center)
        True

        sage: # TestSuite(f).run()
    """
    def reduced_trace(self, var=None):
        r"""
        Return the reduced trace of this element.

        INPUT:

        - ``var`` -- string or ``None`` (default: ``None``);
          the name of the central variable

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x', Frob]
            sage: K = S.fraction_field()
            sage: a = 1 / (x^2 + t)
            sage: tr = a.reduced_trace(); tr
            3/(z^2 + 2)

        The reduced trace lies in the center of `S`, which is the fraction field
        of a univariate polynomial ring in the variable `z = x^3` over `GF(5)`::

            sage: # needs sage.rings.finite_rings
            sage: tr.parent()
            Fraction Field of Univariate Polynomial Ring in z over Finite Field of size 5
            sage: tr.parent() is K.center()
            True

        We can use explicit conversion to view ``tr`` as a Ore function::

            sage: K(tr)                                                                 # needs sage.rings.finite_rings
            (x^6 + 2)^(-1) * 3

        By default, the name of the central variable is usually ``z`` (see
        :meth:`sage.rings.polynomial.skew_polynomial_ring.OreFunctionField_with_large_center.center`
        for more details about this).
        However, the user can specify a different variable name if desired::

            sage: a.reduced_trace(var='u')                                              # needs sage.rings.finite_rings
            3/(u^2 + 2)

        TESTS:

        We check that the reduced trace is additive::

            sage: # needs sage.rings.finite_rings
            sage: a = K.random_element(degree=5)
            sage: b = K.random_element(degree=7)
            sage: a.reduced_trace() + b.reduced_trace() == (a+b).reduced_trace()
            True

        ::

            sage: (a*b).reduced_trace() == (b*a).reduced_trace()                        # needs sage.rings.finite_rings
            True
        """
        ring = self.parent()._ring
        denominator = self._denominator.reduced_norm(var)
        cofactor, _ = ring(denominator).right_quo_rem(self._denominator)
        numerator = (cofactor * self._numerator).reduced_trace(var)
        return numerator/denominator

    def reduced_norm(self, var=None):
        r"""
        Return the reduced norm of this Ore function.

        INPUT:

        - ``var`` -- string or ``None`` (default: ``None``);
          the name of the central variable

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x', Frob]
            sage: K = S.fraction_field()
            sage: a = (x + t) / (x^2 + t^2)
            sage: N = a.reduced_norm(); N
            (z + 2)/(z^2 + 4)

        The reduced norm lies in the center of `S`, which is the fraction field
        of a univariate polynomial ring in the variable `z = x^3` over `GF(5)`. ::

            sage: # needs sage.rings.finite_rings
            sage: N.parent()
            Fraction Field of Univariate Polynomial Ring in z over Finite Field of size 5
            sage: N.parent() is K.center()
            True

        We can use explicit conversion to view ``N`` as a skew polynomial::

            sage: K(N)                                                                  # needs sage.rings.finite_rings
            (x^6 + 4)^(-1) * (x^3 + 2)

        By default, the name of the central variable is usually ``z`` (see
        :meth:`sage.rings.polynomial.skew_polynomial_ring.SkewPolynomiaRing_finite_order.center`
        for more details about this).
        However, the user can specify a different variable name if desired::

            sage: a.reduced_norm(var='u')                                               # needs sage.rings.finite_rings
            (u + 2)/(u^2 + 4)

        TESTS:

        We check that the reduced norm is a multiplicative map::

            sage: # needs sage.rings.finite_rings
            sage: a = K.random_element()
            sage: b = K.random_element()
            sage: a.reduced_norm() * b.reduced_norm() == (a*b).reduced_norm()
            True
        """
        numerator = self._numerator.reduced_norm(var)
        denominator = self._denominator.reduced_norm(var)
        return numerator/denominator
