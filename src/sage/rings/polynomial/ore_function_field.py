r"""
Fraction fields of Ore polynomial rings

Sage provides support for building the fraction field of any Ore
polynomial ring and performing basic operations in it.
The fraction field is constructed by the method
:meth:`sage.rings.polynomial.ore_polynomial_ring.OrePolynomialRing.fraction_field`
as demonstrated below::

    sage: R.<t> = QQ[]
    sage: der = R.derivation()
    sage: A.<d> = R['d', der]
    sage: K = A.fraction_field()
    sage: K
    Ore Function Field in d over Fraction Field of Univariate Polynomial Ring in t
     over Rational Field twisted by d/dt

The simplest way to build elements in `K` is to use the division
operator over Ore polynomial rings::

    sage: f = 1/d
    sage: f
    d^(-1)
    sage: f.parent() is K
    True

REPRESENTATION OF ELEMENTS:

Elements in `K` are internally represented by fractions of the form `s^{-1} t`
with the denominator on the left. Notice that, because of noncommutativity,
this is not the same that fractions with denominator on the right.
For example, a fraction created by the division operator is usually displayed
with a different numerator and/or a different denominator::

    sage: g = t / d
    sage: g
    (d - 1/t)^(-1) * t

The left numerator and right denominator are accessible as follows::

    sage: g.left_numerator()
    t
    sage: g.right_denominator()
    d

Similarly the methods :meth:`OrePolynomial.left_denominator` and
:meth:`OrePolynomial.right_numerator` give access to the Ore polynomials
`s` and `t` in the representation `s^{-1} t`::

    sage: g.left_denominator()
    d - 1/t
    sage: g.right_numerator()
    t

We favored the writing `s^{-1} t` because it always exists.
On the contrary, the writing `s t^{-1}` is only guaranteed when the twisting
morphism defining the Ore polynomial ring is bijective. As a consequence, when
the latter assumption is not fulfilled (or actually if Sage cannot invert the
twisting morphism), computing the left numerator and the right denominator fails::

    sage: # needs sage.rings.function_field
    sage: sigma = R.hom([t^2])
    sage: S.<x> = R['x', sigma]
    sage: F = S.fraction_field()
    sage: f = F.random_element()
    sage: while not f:
    ....:     f = F.random_element()
    sage: f.left_numerator()
    Traceback (most recent call last):
    ...
    NotImplementedError: inversion of the twisting morphism
    Ring endomorphism of Fraction Field of Univariate Polynomial Ring in t over Rational Field
      Defn: t |--> t^2

On a related note, fractions are systematically simplified when the twisting
morphism is bijective but they are not otherwise. As an example, compare the
two following computations::

    sage: # needs sage.rings.function_field
    sage: P = d^2 + t*d + 1
    sage: Q = d + t^2
    sage: D = d^3 + t^2 + 1
    sage: f = P^(-1) * Q
    sage: f
    (d^2 + t*d + 1)^(-1) * (d + t^2)
    sage: g = (D*P)^(-1) * (D*Q)
    sage: g
    (d^2 + t*d + 1)^(-1) * (d + t^2)

    sage: # needs sage.rings.function_field
    sage: P = x^2 + t*x + 1
    sage: Q = x + t^2
    sage: D = x^3 + t^2 + 1
    sage: f = P^(-1) * Q
    sage: f
    (x^2 + t*x + 1)^(-1) * (x + t^2)
    sage: g = (D*P)^(-1) * (D*Q)
    sage: g
    (x^5 + t^8*x^4 + x^3 + (t^2 + 1)*x^2 + (t^3 + t)*x + t^2 + 1)^(-1)
    * (x^4 + t^16*x^3 + (t^2 + 1)*x + t^4 + t^2)
    sage: f == g
    True

OPERATIONS:

Basic arithmetical operations are available::

    sage: # needs sage.rings.function_field
    sage: f = 1 / d
    sage: g = 1 / (d + t)
    sage: u = f + g; u
    (d^2 + ((t^2 - 1)/t)*d)^(-1) * (2*d + (t^2 - 2)/t)
    sage: v = f - g; v
    (d^2 + ((t^2 - 1)/t)*d)^(-1) * t
    sage: u + v
    d^(-1) * 2

    sage: f * g
    (d^2 + t*d)^(-1)
    sage: f / g
    d^(-1) * (d + t)

Of course, multiplication remains noncommutative::

    sage: # needs sage.rings.function_field
    sage: g * f
    (d^2 + t*d + 1)^(-1)
    sage: g^(-1) * f
    (d - 1/t)^(-1) * (d + (t^2 - 1)/t)

TESTS:

The Ore function field is commutative if the twisting morphism is the
identity and the twisting derivation vanishes. ::

    sage: # needs sage.rings.finite_rings
    sage: k.<a> = GF(5^3)
    sage: Frob = k.frobenius_endomorphism()
    sage: S.<x> = k['x', Frob]
    sage: K = S.fraction_field()
    sage: K.is_commutative()
    False
    sage: T.<y> = k['y', Frob^3]
    sage: L = T.fraction_field()
    sage: L.is_commutative()
    True

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

from sage.categories.algebras import Algebras
from sage.categories.commutative_rings import CommutativeRings
from sage.categories.fields import Fields
from sage.categories.homset import Hom
from sage.categories.map import Section
from sage.rings.morphism import RingHomomorphism
from sage.rings.polynomial.ore_function_element import (
    OreFunction_with_large_center,
    OreFunctionBaseringInjection,
)
from sage.rings.polynomial.ore_polynomial_ring import OrePolynomialRing
from sage.structure.category_object import normalize_names
from sage.structure.parent import Parent
from sage.structure.richcmp import op_EQ
from sage.structure.unique_representation import UniqueRepresentation

WORKING_CENTER_MAX_TRIES = 1000


# Generic implementation of Ore function fields
###############################################

class OreFunctionField(Parent, UniqueRepresentation):
    r"""
    A class for fraction fields of Ore polynomial rings.
    """
    Element = None

    def __init__(self, ring, category=None):
        r"""
        Initialize this Ore function field.

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(11^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: der = k.derivation(a, twist=Frob)
            sage: S.<x> = k['x', der]
            sage: K = S.fraction_field()
            sage: TestSuite(K).run()
        """
        if self.Element is None:
            import sage.rings.polynomial.ore_function_element
            self.Element = sage.rings.polynomial.ore_function_element.OreFunction
        if not isinstance(ring, OrePolynomialRing):
            raise TypeError("not a Ore Polynomial Ring")
        if ring.base_ring() not in Fields():
            raise TypeError("the base ring must be a field")
        try:
            _ = ring.twisting_morphism(-1)
            self._simplification = True
        except (TypeError, ZeroDivisionError, NotImplementedError):
            self._simplification = False
        self._ring = ring
        base = ring.base_ring()
        if ring in CommutativeRings():
            category = Algebras(base).Commutative().or_subcategory(category)
        else:
            category = Algebras(base).or_subcategory(category)
        Parent.__init__(self, base=base, names=ring.variable_name(),
                        normalize=True, category=category)

    def _element_constructor_(self, *args, **kwds):
        r"""
        Construct an element in this Ore function field.

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x', Frob]
            sage: K = S.fraction_field()
            sage: f = K(x^2 + a, x + a^2)  # indirect doctest
            sage: f
            (x + a^2)^(-1) * (x^2 + a)
        """
        return self.element_class(self, *args, **kwds)

    def _coerce_map_from_base_ring(self):
        r"""
        Return a coercion map from the base ring of this Ore function field.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<t> = GF(11^2)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x', Frob]
            sage: K = S.fraction_field()
            sage: K.coerce_map_from(k)  # indirect doctest
            Ore Function base injection morphism:
              From: Finite Field in t of size 11^2
              To:   Ore Function Field in x over Finite Field in t of size 11^2
                    twisted by t |--> t^11
        """
        return OreFunctionBaseringInjection(self.base_ring(), self)

    def _coerce_map_from_(self, P):
        r"""
        Return ``True`` if `P` coerces to this Ore function field.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: der = R.derivation()
            sage: A.<d> = R['d', der]
            sage: K = A.fraction_field()

            sage: K.has_coerce_map_from(A)
            True
            sage: K.has_coerce_map_from(R)
            True
            sage: K.has_coerce_map_from(R.fraction_field())
            True
        """
        if isinstance(P, OreFunctionField):
            return P._ring.has_coerce_map_from(self._ring)
        if isinstance(P, Parent):
            return P.has_coerce_map_from(self._ring)

    def _repr_(self) -> str:
        r"""
        Return a string representation of this Ore function field.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = OrePolynomialRing(R, sigma)
            sage: S.fraction_field()                                                    # needs sage.rings.function_field
            Ore Function Field in x over Fraction Field of Univariate Polynomial Ring in t over Rational Field twisted by t |--> t + 1

            sage: der = R.derivation()
            sage: T.<d> = OrePolynomialRing(R, der)
            sage: T.fraction_field()                                                    # needs sage.rings.function_field
            Ore Function Field in d over Fraction Field of Univariate Polynomial Ring in t over Rational Field twisted by d/dt
        """
        s = "Ore Function Field in %s over %s twisted by " % (self.variable_name(), self.base_ring())
        morphism = self.twisting_morphism()
        derivation = self.twisting_derivation()
        if derivation is None:
            s += morphism._repr_short()
        else:
            if morphism is not None:
                s += "%s and " % morphism._repr_short()
            s += derivation._repr_()
        return s

    def _latex_(self):
        r"""
        Return a LaTeX representation of this Ore function field.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x', Frob]
            sage: K = S.fraction_field()
            sage: latex(K)  # indirect doctest
            \Bold{F}_{5^{3}}\left(x ; a \mapsto a^{5} \right)

            sage: R.<t> = QQ[]
            sage: der = R.derivation()
            sage: T.<delta> = R['delta', der]
            sage: L = T.fraction_field()
            sage: latex(L)  # indirect doctest
            \mathrm{Frac}(\Bold{Q}[t])\left(\delta ; \frac{d}{dt} \right)
        """
        from sage.misc.latex import latex
        s = "%s\\left(%s" % (latex(self.base_ring()), self.latex_variable_names()[0])
        sep = ";"
        if self.twisting_morphism() is not None:
            s += sep + latex(self.twisting_morphism())
            sep = ","
        if self.twisting_derivation() is not None:
            s += sep + latex(self.twisting_derivation())
        return s + "\\right)"

    def change_var(self, var):
        r"""
        Return the Ore function field in variable ``var`` with the same base
        ring, twisting morphism and twisting derivation as ``self``.

        INPUT:

        - ``var`` -- string representing the name of the new variable

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: R.<x> = OrePolynomialRing(k,Frob)
            sage: K = R.fraction_field()
            sage: K
            Ore Function Field in x over Finite Field in t of size 5^3 twisted by t |--> t^5
            sage: Ky = K.change_var('y'); Ky
            Ore Function Field in y over Finite Field in t of size 5^3 twisted by t |--> t^5
            sage: Ky is K.change_var('y')
            True
        """
        return OreFunctionField(self._ring.change_var(var))

    def characteristic(self):
        r"""
        Return the characteristic of this Ore function field.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S = R['x',sigma]
            sage: S.fraction_field().characteristic()                                   # needs sage.rings.function_field
            0

            sage: # needs sage.rings.finite_rings
            sage: k.<u> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S = k['y',Frob]
            sage: S.fraction_field().characteristic()                                   # needs sage.rings.function_field
            5
        """
        return self.base_ring().characteristic()

    def twisting_morphism(self, n=1):
        r"""
        Return the twisting endomorphism defining this Ore function field iterated ``n`` times
        or ``None`` if this Ore function field is not twisted by an endomorphism.

        INPUT:

        - ``n`` -- integer (default: 1)

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x', sigma]
            sage: K = S.fraction_field()                                                # needs sage.rings.function_field
            sage: K.twisting_morphism()                                                 # needs sage.rings.function_field
            Ring endomorphism of
             Fraction Field of Univariate Polynomial Ring in t over Rational Field
              Defn: t |--> t + 1

        When the Ore polynomial ring is only twisted by a derivation, this
        method returns nothing::

            sage: der = R.derivation()
            sage: A.<d> = R['x', der]
            sage: F = A.fraction_field()                                                # needs sage.rings.function_field
            sage: F.twisting_morphism()                                                 # needs sage.rings.function_field

        .. SEEALSO::

            :meth:`sage.rings.polynomial.ore_polynomial_element.OrePolynomial.twisting_morphism`,
            :meth:`twisting_derivation`
        """
        return self._ring.twisting_morphism(n)

    def twisting_derivation(self):
        r"""
        Return the twisting derivation defining this Ore function field
        or ``None`` if this Ore function field is not twisted by a derivation.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: der = R.derivation(); der
            d/dt
            sage: A.<d> = R['d', der]
            sage: F = A.fraction_field()
            sage: F.twisting_derivation()
            d/dt

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x', Frob]
            sage: K = S.fraction_field()
            sage: K.twisting_derivation()

        .. SEEALSO::

            :meth:`sage.rings.polynomial.ore_polynomial_element.OrePolynomial.twisting_derivation`,
            :meth:`twisting_morphism`
        """
        return self._ring.twisting_derivation()

    def gen(self, n=0):
        r"""
        Return the indeterminate generator of this Ore function field.

        INPUT:

        - ``n`` -- index of generator to return (default: 0); exists for
          compatibility with other polynomial rings

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(5^4)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x', Frob]
            sage: K = S.fraction_field()
            sage: K.gen()
            x
        """
        return self(self._ring.gen(n))

    def gens(self) -> tuple:
        """
        Return the tuple of generators of ``self``.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(5^4)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x', Frob]
            sage: K = S.fraction_field()
            sage: K.gens()
            (x,)
        """
        return (self(self._ring.gen(0)),)

    parameter = gen

    def gens_dict(self) -> dict:
        r"""
        Return a {name: variable} dictionary of the generators of
        this Ore function field.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: R.<t> = ZZ[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = OrePolynomialRing(R, sigma)
            sage: K = S.fraction_field()
            sage: K.gens_dict()
            {'x': x}
        """
        return dict(zip(self.variable_names(), self.gens()))

    def is_finite(self) -> bool:
        r"""
        Return ``False`` since Ore function field are not finite.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<t> = GF(5^3)
            sage: k.is_finite()
            True
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: K = S.fraction_field()
            sage: K.is_finite()
            False
        """
        return False

    def is_exact(self) -> bool:
        r"""
        Return ``True`` if elements of this Ore function field are exact.
        This happens if and only if elements of the base ring are exact.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x', Frob]
            sage: K = S.fraction_field()
            sage: K.is_exact()
            True

            sage: # needs sage.rings.padics
            sage: k.<u> = Qq(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x', Frob]
            sage: K = S.fraction_field()
            sage: K.is_exact()
            False
        """
        return self.base_ring().is_exact()

    def is_sparse(self) -> bool:
        r"""
        Return ``True`` if the elements of this Ore function field are sparsely
        represented.

        .. WARNING::

            Since sparse Ore polynomials are not yet implemented, this
            function always returns ``False``.

        EXAMPLES::

            sage: # needs sage.rings.function_field sage.rings.real_mpfr
            sage: R.<t> = RR[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x', sigma]
            sage: K = S.fraction_field()
            sage: K.is_sparse()
            False
        """
        return self._ring.is_sparse()

    def ngens(self) -> int:
        r"""
        Return the number of generators of this Ore function field,
        which is `1`.

        EXAMPLES::

            sage: # needs sage.rings.function_field sage.rings.real_mpfr
            sage: R.<t> = RR[]
            sage: sigma = R.hom([t+1])
            sage: S.<x> = R['x',sigma]
            sage: K = S.fraction_field()
            sage: K.ngens()
            1
        """
        return 1

    def random_element(self, degree=2, monic=False, *args, **kwds):
        r"""
        Return a random Ore function in this field.

        INPUT:

        - ``degree`` -- (default: 2) an integer or a list of
          two integers; the degrees of the denominator and numerator

        - ``monic`` -- boolean (default: ``False``); if ``True``, return a monic
          Ore function with monic numerator and denominator

        - ``*args``, ``**kwds`` -- passed in to the :meth:`random_element`
          method for the base ring

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x', Frob]
            sage: K = S.fraction_field()
            sage: K.random_element()              # random
            (x^2 + (2*t^2 + t + 1)*x + 2*t^2 + 2*t + 3)^(-1)
            * ((2*t^2 + 3)*x^2 + (4*t^2 + t + 4)*x + 2*t^2 + 2)
            sage: K.random_element(monic=True)    # random
            (x^2 + (4*t^2 + 3*t + 4)*x + 4*t^2 + t)^(-1)
            * (x^2 + (2*t^2 + t + 3)*x + 3*t^2 + t + 2)
            sage: K.random_element(degree=3)      # random
            (x^3 + (2*t^2 + 3)*x^2 + (2*t^2 + 4)*x + t + 3)^(-1)
            * ((t + 4)*x^3 + (4*t^2 + 2*t + 2)*x^2 + (2*t^2 + 3*t + 3)*x + 3*t^2 + 3*t + 1)
            sage: K.random_element(degree=[2,5])  # random
            (x^2 + (4*t^2 + 2*t + 2)*x + 4*t^2 + t + 2)^(-1)
            * ((3*t^2 + t + 1)*x^5 + (2*t^2 + 2*t)*x^4 + (t^2 + 2*t + 4)*x^3
               + (3*t^2 + 2*t)*x^2 + (t^2 + t + 4)*x)
        """
        if isinstance(degree, list):
            degdenom, degnum = degree
        else:
            degdenom = degnum = degree
        numerator = self._ring.random_element(degnum, monic, *args, **kwds)
        denominator = self._ring.random_element(degdenom, True, *args, **kwds)
        return self(numerator, denominator)

    def is_field(self, proof=False) -> bool:
        r"""
        Return always ``True`` since Ore function field are (skew) fields.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x', Frob]
            sage: K = S.fraction_field()
            sage: S.is_field()
            False
            sage: K.is_field()
            True

        TESTS:

        We check that :issue:`31470` is fixed::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(5^3)
            sage: S.<x> = k['x', k.frobenius_endomorphism()]
            sage: K = S.fraction_field()
            sage: zero_matrix(K, 2).row(0)
            ...
            (0, 0)
        """
        return True

    def fraction_field(self):
        r"""
        Return the fraction field of this Ore function field,
        i.e. this Ore function field itself.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: der = R.derivation()
            sage: A.<d> = R['d', der]
            sage: K = A.fraction_field(); K
            Ore Function Field in d
             over Fraction Field of Univariate Polynomial Ring in t over Rational Field
             twisted by d/dt
            sage: K.fraction_field()
            Ore Function Field in d
             over Fraction Field of Univariate Polynomial Ring in t over Rational Field
             twisted by d/dt
            sage: K.fraction_field() is K
            True
        """
        return self


# Special classes for twisting morphisms with finite order
##########################################################

class SectionOreFunctionCenterInjection(Section):
    r"""
    Section of the canonical injection of the center of a Ore
    function field into this field
    """
    def __init__(self, embed):
        r"""
        Initialize this map.

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(5^3)
            sage: S.<x> = OrePolynomialRing(k, k.frobenius_endomorphism())
            sage: K = S.fraction_field()
            sage: Z = K.center()
            sage: iota = K.coerce_map_from(Z)
            sage: sigma = iota.section()
            sage: TestSuite(sigma).run(skip=['_test_category'])
        """
        Section.__init__(self, embed)
        self._ringsection = embed._ringembed.section()
        self._simplify = not embed._codomain._simplification

    def _call_(self, x):
        r"""
        Return `x` viewed as an element of the center.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(5^3)
            sage: S.<x> = SkewPolynomialRing(k, k.frobenius_endomorphism())
            sage: K = S.fraction_field()
            sage: Z = K.center()
            sage: iota = K.coerce_map_from(Z)
            sage: sigma = iota.section()
            sage: sigma(1/x^3)
            1/z
            sage: sigma(1/x^2)
            Traceback (most recent call last):
            ...
            ValueError: x^(-2) is not in the center
        """
        numerator = x._numerator
        denominator = x._denominator
        if self._simplify:
            D = numerator.right_gcd(denominator)
            numerator, _ = numerator.right_quo_rem(D)
            denominator, _ = denominator.right_quo_rem(D)
        try:
            return self._ringsection(numerator) / self._ringsection(denominator)
        except ValueError:
            raise ValueError("%s is not in the center" % x)

    def _richcmp_(self, right, op):
        r"""
        Compare this morphism with ``right``.

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(5^3)
            sage: S.<x> = SkewPolynomialRing(k, k.frobenius_endomorphism())
            sage: K = S.fraction_field()
            sage: Z = K.center()
            sage: iota = K.coerce_map_from(Z)
            sage: sigma = iota.section()
            sage: s = loads(dumps(sigma))
            sage: s == sigma
            True
            sage: s is sigma
            False
        """
        if op == op_EQ:
            return (self.domain() is right.domain()) and (self.codomain() is right.codomain())
        return NotImplemented


class OreFunctionCenterInjection(RingHomomorphism):
    r"""
    Canonical injection of the center of a Ore function field
    into this field.
    """
    def __init__(self, domain, codomain, ringembed):
        r"""
        Initialize this morphism.

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(5^3)
            sage: S.<x> = SkewPolynomialRing(k, k.frobenius_endomorphism())
            sage: K = S.fraction_field()
            sage: Z = K.center()
            sage: iota = K.coerce_map_from(Z)
            sage: TestSuite(iota).run(skip=['_test_category'])
        """
        RingHomomorphism.__init__(self, Hom(domain, codomain))
        self._codomain = codomain
        self._ringembed = ringembed
        self._section = SectionOreFunctionCenterInjection(self)

    def _repr_(self):
        r"""
        Return a string representation of this morphism.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(5^3)
            sage: S.<x> = SkewPolynomialRing(k, k.frobenius_endomorphism())
            sage: K = S.fraction_field()
            sage: Z = K.center()
            sage: iota = K.coerce_map_from(Z)
            sage: iota
            Embedding of the center of Ore Function Field in x over Finite Field in a of size 5^3 twisted by a |--> a^5 into this field
            sage: iota._repr_()
            'Embedding of the center of Ore Function Field in x over Finite Field in a of size 5^3 twisted by a |--> a^5 into this field'
        """
        return "Embedding of the center of %s into this field" % self._codomain

    def _call_(self, x):
        r"""
        Return the image of `x` by this morphism.

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(5^3)
            sage: S.<x> = SkewPolynomialRing(k, k.frobenius_endomorphism())
            sage: K = S.fraction_field()
            sage: Z.<z> = K.center()
            sage: iota = K.coerce_map_from(Z)
            sage: iota(1/(z+1))
            (x^3 + 1)^(-1)
        """
        numerator = self._ringembed(x.numerator())
        denominator = self._ringembed(x.denominator())
        return self._codomain(numerator, denominator, simplify=False)

    def _richcmp_(self, right, op):
        r"""
        Compare this morphism with ``right``.

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(5^3)
            sage: S.<x> = SkewPolynomialRing(k, k.frobenius_endomorphism())
            sage: K = S.fraction_field()
            sage: Z = K.center()
            sage: iota = K.coerce_map_from(Z)
            sage: i = loads(dumps(iota))
            sage: i == iota
            True
            sage: i is iota
            False
        """
        if op == op_EQ:
            return (self.domain() is right.domain()) and (self.codomain() is right.codomain())
        return NotImplemented

    def section(self):
        r"""
        Return a section of this morphism.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(5^3)
            sage: S.<x> = SkewPolynomialRing(k, k.frobenius_endomorphism())
            sage: K = S.fraction_field()
            sage: Z = K.center()
            sage: iota = K.coerce_map_from(Z)
            sage: sigma = iota.section()
            sage: sigma(x^3 / (x^6 + 1))
            z/(z^2 + 1)
        """
        return self._section


class OreFunctionField_with_large_center(OreFunctionField):
    """
    A specialized class for Ore polynomial fields whose center has finite index.
    """
    def __init__(self, ring, category=None):
        r"""
        Initialize this Ore function field.

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x', Frob]
            sage: K = S.fraction_field()
            sage: TestSuite(K).run()
        """
        if self.Element is None:
            self.Element = OreFunction_with_large_center
        OreFunctionField.__init__(self, ring, category)
        self._center = {}
        self._center_variable_name = 'z'
        for i in range(WORKING_CENTER_MAX_TRIES):
            try:
                self._working_center = self.center()
                self._center_variable_name = None
                break
            except ValueError:
                self._center_variable_name = "z%s_" % i
        if self._center_variable_name is not None:
            raise NotImplementedError("unable to create the center")

    def center(self, name=None, names=None, default=False):
        r"""
        Return the center of this Ore function field.

        .. NOTE::

            One can prove that the center is a field of rational functions
            over a subfield of the base ring of this Ore function field.

        INPUT:

        - ``name`` -- string or ``None`` (default: ``None``);
          the name for the central variable

        - ``default`` -- boolean (default: ``False``); if ``True``,
          set the default variable name for the center to ``name``

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: S.<x> = k['x',Frob]
            sage: K = S.fraction_field()
            sage: Z = K.center(); Z
            Fraction Field of Univariate Polynomial Ring in z over Finite Field of size 5

        We can pass in another variable name::

            sage: K.center(name='y')                                                    # needs sage.rings.finite_rings
            Fraction Field of Univariate Polynomial Ring in y over Finite Field of size 5

        or use the bracket notation::

            sage: Zy.<y> = K.center(); Zy                                               # needs sage.rings.finite_rings
            Fraction Field of Univariate Polynomial Ring in y over Finite Field of size 5

        A coercion map from the center to the Ore function field is set::

            sage: K.has_coerce_map_from(Zy)                                             # needs sage.rings.finite_rings
            True

        and pushout works::

            sage: # needs sage.rings.finite_rings
            sage: x.parent()
            Ore Polynomial Ring in x over Finite Field in t of size 5^3 twisted by t |--> t^5
            sage: y.parent()
            Fraction Field of Univariate Polynomial Ring in y over Finite Field of size 5
            sage: P = x + y; P
            x^3 + x
            sage: P.parent()
            Ore Function Field in x over Finite Field in t of size 5^3 twisted by t |--> t^5

        A conversion map in the reverse direction is also set::

            sage: # needs sage.rings.finite_rings
            sage: Zy(x^(-6) + 2)
            (2*y^2 + 1)/y^2
            sage: Zy(1/x^2)
            Traceback (most recent call last):
            ...
            ValueError: x^(-2) is not in the center

        ABOUT THE DEFAULT NAME OF THE CENTRAL VARIABLE:

        A priori, the default is ``z``.

        However, a variable name is given the first time this method is
        called, the given name become the default for the next calls::

            sage: # needs sage.rings.finite_rings
            sage: k.<t> = GF(11^3)
            sage: phi = k.frobenius_endomorphism()
            sage: S.<X> = k['X', phi]
            sage: K = S.fraction_field()
            sage: C.<u> = K.center()  # first call
            sage: C
            Fraction Field of Univariate Polynomial Ring in u over Finite Field of size 11
            sage: K.center()  # second call: the variable name is still u
            Fraction Field of Univariate Polynomial Ring in u over Finite Field of size 11

        We can update the default variable name by passing in the argument
        ``default=True``::

            sage: # needs sage.rings.finite_rings
            sage: D.<v> = K.center(default=True)
            sage: D
            Fraction Field of Univariate Polynomial Ring in v over Finite Field of size 11
            sage: K.center()
            Fraction Field of Univariate Polynomial Ring in v over Finite Field of size 11
        """
        if name is not None and names is not None:
            raise ValueError("you must specify the name of the variable")
        if names is None:
            if name is None:
                name = self._center_variable_name
            if name is None:
                name = 'z'
            names = (name,)
        names = normalize_names(1, names)
        name = names[0]
        if name in self._center:
            center = self._center[name]
        else:
            ring = self._ring
            ringcenter = ring.center(name, default=False)
            ringembed = ring.coerce_map_from(ringcenter)
            center = ringcenter.fraction_field()
            embed = OreFunctionCenterInjection(center, self, ringembed)
            try:
                assert not self.has_coerce_map_from(center)
                self.register_coercion(embed)
                center.register_conversion(embed.section())
            except AssertionError:
                raise ValueError("creation of coercion map fails; consider using another variable name")
            self._center[name] = center
        if default or (self._center_variable_name is None):
            self._center_variable_name = name
        return center
