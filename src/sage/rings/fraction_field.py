r"""
Fraction Field of Integral Domains

AUTHORS:

- William Stein (with input from David Joyner, David Kohel, and Joe
  Wetherell)

- Burcin Erocal

- Julian Rüth (2017-06-27): embedding into the field of fractions and its
  section

EXAMPLES:

Quotienting is a constructor for an element of the fraction field::

    sage: R.<x> = QQ[]
    sage: (x^2-1)/(x+1)
    x - 1
    sage: parent((x^2-1)/(x+1))
    Fraction Field of Univariate Polynomial Ring in x over Rational Field

The GCD is not taken (since it doesn't converge sometimes) in the
inexact case::

    sage: # needs sage.rings.real_mpfr
    sage: Z.<z> = CC[]
    sage: I = CC.gen()
    sage: (1+I+z)/(z+0.1*I)
    (z + 1.00000000000000 + I)/(z + 0.100000000000000*I)
    sage: (1+I*z)/(z+1.1)
    (I*z + 1.00000000000000)/(z + 1.10000000000000)

TESTS::

    sage: F = FractionField(IntegerRing())
    sage: F == loads(dumps(F))
    True

::

    sage: F = FractionField(PolynomialRing(RationalField(),'x'))
    sage: F == loads(dumps(F))
    True

::

    sage: F = FractionField(PolynomialRing(IntegerRing(),'x'))
    sage: F == loads(dumps(F))
    True

::

    sage: F = FractionField(PolynomialRing(RationalField(),2,'x'))
    sage: F == loads(dumps(F))
    True

Test that :issue:`15971` is fixed::

    sage: for B in [QQ['t'], QQ['s, t'], ZZ['t'], ZZ['s, t']]:
    ....:     F = B.fraction_field()
    ....:     R = F['x, y']
    ....:     x = R.gen(0)
    ....:     print(x / x)
    1
    1
    1
    1
"""
# ****************************************************************************
#
#   Sage: Open Source Mathematical Software
#
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#                     2017 Julian Rüth <julian.rueth@fsfe.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import sage.misc.latex as latex
from sage.categories.basic import QuotientFields, Rings
from sage.categories.map import Section
from sage.misc.cachefunc import cached_method
from sage.rings import fraction_field_element, ring
from sage.rings.integer_ring import ZZ
from sage.structure.coerce import py_scalar_to_element
from sage.structure.coerce_maps import CallableConvertMap, DefaultConvertMap_unique
from sage.structure.element import parent
from sage.structure.parent import Parent
from sage.structure.richcmp import richcmp


def FractionField(R, names=None):
    """
    Create the fraction field of the integral domain ``R``.

    INPUT:

    - ``R`` -- an integral domain

    - ``names`` -- ignored

    EXAMPLES:

    We create some example fraction fields::

        sage: FractionField(IntegerRing())
        Rational Field
        sage: FractionField(PolynomialRing(RationalField(),'x'))
        Fraction Field of Univariate Polynomial Ring in x over Rational Field
        sage: FractionField(PolynomialRing(IntegerRing(),'x'))
        Fraction Field of Univariate Polynomial Ring in x over Integer Ring
        sage: FractionField(PolynomialRing(RationalField(),2,'x'))
        Fraction Field of Multivariate Polynomial Ring in x0, x1 over Rational Field

    Dividing elements often implicitly creates elements of the fraction
    field::

        sage: x = PolynomialRing(RationalField(), 'x').gen()
        sage: f = x/(x+1)
        sage: g = x**3/(x+1)
        sage: f/g
        1/x^2
        sage: g/f
        x^2

    The input must be an integral domain::

        sage: Frac(Integers(4))
        Traceback (most recent call last):
        ...
        TypeError: R must be an integral domain
    """
    if R not in Rings():
        raise TypeError("R must be a ring")
    if not R.is_integral_domain():
        raise TypeError("R must be an integral domain")
    return R.fraction_field()


def is_FractionField(x) -> bool:
    """
    Test whether or not ``x`` inherits from :class:`FractionField_generic`.

    EXAMPLES::

        sage: from sage.rings.fraction_field import is_FractionField
        sage: is_FractionField(Frac(ZZ['x']))
        doctest:warning...
        DeprecationWarning: The function is_FractionField is deprecated;
        use 'isinstance(..., FractionField_generic)' instead.
        See https://github.com/sagemath/sage/issues/38128 for details.
        True
        sage: is_FractionField(QQ)
        False
    """
    from sage.misc.superseded import deprecation
    deprecation(38128,
                "The function is_FractionField is deprecated; "
                "use 'isinstance(..., FractionField_generic)' instead.")
    return isinstance(x, FractionField_generic)


class FractionField_generic(ring.Field):
    """
    The fraction field of an integral domain.
    """
    def __init__(self, R,
                 element_class=fraction_field_element.FractionFieldElement,
                 category=QuotientFields()):
        """
        Create the fraction field of the integral domain ``R``.

        INPUT:

        - ``R`` -- an integral domain

        EXAMPLES::

            sage: Frac(QQ['x'])
            Fraction Field of Univariate Polynomial Ring in x over Rational Field
            sage: Frac(QQ['x,y']).variable_names()
            ('x', 'y')
            sage: category(Frac(QQ['x']))
            Category of infinite quotient fields

        TESTS::

            sage: F = FractionField(QQ['x'])
            sage: F.cardinality()
            +Infinity
        """
        self._R = R
        self._element_class = element_class
        cat = category
        if R in Rings().Infinite():
            cat = cat.Infinite()
        elif R in Rings().Finite():
            cat = cat.Finite()
        Parent.__init__(self, base=R, names=R._names, category=cat)

    def __reduce__(self):
        """
        For pickling.

        TESTS::

            sage: K = Frac(QQ['x'])
            sage: loads(dumps(K)) is K
            True
        """
        return FractionField, (self._R,)

    def _coerce_map_from_(self, S):
        """
        Return ``True`` if elements of ``S`` can be coerced into this
        fraction field.

        This fraction field has coercions from:

        - itself
        - any fraction field where the base ring coerces to the base
          ring of this fraction field
        - any ring that coerces to the base ring of this fraction field

        EXAMPLES::

            sage: F = QQ['x,y'].fraction_field()
            sage: F.has_coerce_map_from(F) # indirect doctest
            True

        ::

            sage: F.has_coerce_map_from(ZZ['x,y'].fraction_field())
            True

        ::

            sage: F.has_coerce_map_from(ZZ['x,y,z'].fraction_field())
            False

        ::

            sage: F.has_coerce_map_from(ZZ)
            True

        Test coercions::

            sage: F.coerce(1)
            1
            sage: F.coerce(int(1))
            1
            sage: F.coerce(1/2)
            1/2

        ::

            sage: K = ZZ['x,y'].fraction_field()
            sage: x,y = K.gens()
            sage: F.coerce(F.gen())
            x
            sage: F.coerce(x)
            x
            sage: F.coerce(x/y)
            x/y
            sage: L = ZZ['x'].fraction_field()
            sage: K.coerce(L.gen())
            x

        We demonstrate that :issue:`7958` is resolved in the case of
        number fields::

            sage: # needs sage.rings.number_field
            sage: _.<x> = ZZ[]
            sage: K.<a> = NumberField(x^5 - 3*x^4 + 2424*x^3 + 2*x - 232)
            sage: R = K.ring_of_integers()
            sage: S.<y> = R[]
            sage: F = FractionField(S)
            sage: F(1/a)
            (a^4 - 3*a^3 + 2424*a^2 + 2)/232

        Some corner cases have been known to fail in the past (:issue:`5917`)::

            sage: F1 = FractionField( QQ['a'] )
            sage: R12 = F1['x','y']
            sage: R12('a')
            a
            sage: F1(R12(F1('a')))
            a

            sage: F2 = FractionField( QQ['a','b'] )
            sage: R22 = F2['x','y']
            sage: R22('a')
            a
            sage: F2(R22(F2('a')))
            a

        Coercion from Laurent polynomials now works (:issue:`15345`)::

            sage: R = LaurentPolynomialRing(ZZ, 'x')
            sage: T = PolynomialRing(ZZ, 'x')
            sage: R.gen() + FractionField(T).gen()
            2*x
            sage: 1/(R.gen() + 1)
            1/(x + 1)

            sage: # needs sage.modules
            sage: R = LaurentPolynomialRing(ZZ, 'x,y')
            sage: FF = FractionField(PolynomialRing(ZZ, 'x,y'))
            sage: prod(R.gens()) + prod(FF.gens())
            2*x*y
            sage: 1/(R.gen(0) + R.gen(1))
            1/(x + y)

        Test for :issue:`31320`::

            sage: FQt = Frac(QQ['t'])
            sage: LCt = LaurentPolynomialRing(CC,'t')                                   # needs sage.rings.real_mpfr
            sage: coercion_model.common_parent(FQt, LCt)                                # needs sage.rings.real_mpfr
            Fraction Field of Univariate Polynomial Ring in t
            over Complex Field with 53 bits of precision

        Coercion from a localization::

            sage: # needs sage.libs.pari
            sage: R.<x> = ZZ[]
            sage: L = Localization(R, (x**2 + 1,7))
            sage: F = L.fraction_field()
            sage: f = F.coerce_map_from(L); f
            Coercion map:
              From: Univariate Polynomial Ring in x over Integer Ring localized at (7, x^2 + 1)
              To:   Fraction Field of Univariate Polynomial Ring in x over Integer Ring
            sage: f(L(1/7)) == 1/7
            True
        """
        from sage.rings.number_field.number_field_base import NumberField
        from sage.rings.polynomial.laurent_polynomial_ring_base import (
            LaurentPolynomialRing_generic,
        )
        from sage.rings.rational_field import QQ

        if S is self._R:
            parent = self._R.Hom(self)
            return parent.__make_element_class__(FractionFieldEmbedding)(self._R, self, category=parent.homset_category())

        def wrapper(x):
            return self._element_class(self, x.numerator(), x.denominator())

        # The case ``S`` being `\QQ` requires special handling since `\QQ` is
        # not implemented as a ``FractionField_generic``.
        if S is QQ and self._R.has_coerce_map_from(ZZ):
            return CallableConvertMap(S, self, wrapper, parent_as_first_arg=False)

        # special treatment for localizations
        from sage.rings.localization import Localization
        if isinstance(S, Localization):
            parent = S.Hom(self)
            return parent.__make_element_class__(FractionFieldEmbedding)(S, self, category=parent.homset_category())

        # Number fields also need to be handled separately.
        if isinstance(S, NumberField):
            return CallableConvertMap(S, self,
                                      self._number_field_to_frac_of_ring_of_integers,
                                      parent_as_first_arg=False)

        # special treatment for LaurentPolynomialRings
        if (isinstance(S, LaurentPolynomialRing_generic) and
                self._R.fraction_field().has_coerce_map_from(S.base_ring())):

            def converter(x, y=None):
                if y is None:
                    return self._element_class(self, *x._fraction_pair())
                xnum, xden = x._fraction_pair()
                ynum, yden = y._fraction_pair()
                return self._element_class(self, xnum * yden, xden * ynum)
            return CallableConvertMap(S, self, converter, parent_as_first_arg=False)

        if (isinstance(S, FractionField_generic) and
                self._R.has_coerce_map_from(S.ring())):
            return CallableConvertMap(S, self, wrapper, parent_as_first_arg=False)

        if self._R.has_coerce_map_from(S):
            return CallableConvertMap(S, self, self._element_class,
                                      parent_as_first_arg=True)

        return None

    def _number_field_to_frac_of_ring_of_integers(self, x):
        r"""
        Return the number field element ``x`` as an element of ``self``,
        explicitly treating the numerator of ``x``  as an element of the ring
        of integers and the denominator as an integer.

        INPUT:

        - ``x`` -- number field element

        OUTPUT: Element of ``self``

        TESTS:

        We demonstrate that :issue:`7958` is resolved in the case of
        number fields::

            sage: # needs sage.rings.number_field
            sage: _.<x> = ZZ[]
            sage: K.<a> = NumberField(x^5 - 3*x^4 + 2424*x^3 + 2*x - 232)
            sage: R = K.ring_of_integers()
            sage: S.<y> = R[]
            sage: F = FractionField(S)  # indirect doctest
            sage: F(1/a)
            (a^4 - 3*a^3 + 2424*a^2 + 2)/232
        """
        f = x.polynomial()   # Polynomial over QQ
        d = f.denominator()  # Integer
        return self._element_class(self, numerator=d * x, denominator=d)

    def is_field(self, proof=True):
        """
        Return ``True``, since the fraction field is a field.

        EXAMPLES::

            sage: Frac(ZZ).is_field()
            True
        """
        return True

    def is_finite(self):
        """
        Tells whether this fraction field is finite.

        .. NOTE::

            A fraction field is finite if and only if the associated
            integral domain is finite.

        EXAMPLES::

            sage: Frac(QQ['a','b','c']).is_finite()
            False
        """
        return self._R.is_finite()

    def base_ring(self):
        """
        Return the base ring of ``self``.

        This is the base ring of the ring
        which this fraction field is the fraction field of.

        EXAMPLES::

            sage: R = Frac(ZZ['t'])
            sage: R.base_ring()
            Integer Ring
        """
        return self._R.base_ring()

    def characteristic(self):
        """
        Return the characteristic of this fraction field.

        EXAMPLES::

            sage: R = Frac(ZZ['t'])
            sage: R.base_ring()
            Integer Ring
            sage: R = Frac(ZZ['t']); R.characteristic()
            0
            sage: R = Frac(GF(5)['w']); R.characteristic()
            5
        """
        return self._R.characteristic()

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: Frac(ZZ['x']) # indirect doctest
            Fraction Field of Univariate Polynomial Ring in x over Integer Ring
        """
        return "Fraction Field of %s" % self._R

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: latex(Frac(GF(7)['x,y,z'])) # indirect doctest
            \mathrm{Frac}(\Bold{F}_{7}[x, y, z])
        """
        return "\\mathrm{Frac}(%s)" % latex.latex(self._R)

    def _magma_init_(self, magma):
        """
        Return a string representation of ``self`` in the given magma instance.

        EXAMPLES::

            sage: QQ['x'].fraction_field()._magma_init_(magma)            # optional - magma
            'SageCreateWithNames(FieldOfFractions(SageCreateWithNames(PolynomialRing(_sage_ref...),["x"])),["x"])'
            sage: GF(9,'a')['x,y,z'].fraction_field()._magma_init_(magma)       # optional - magma, needs sage.rings.finite_rings
            'SageCreateWithNames(FieldOfFractions(SageCreateWithNames(PolynomialRing(_sage_ref...,3,"grevlex"),["x","y","z"])),["x","y","z"])'

        ``_magma_init_`` gets called implicitly below::

            sage: magma(QQ['x,y'].fraction_field())                  # optional - magma
            Multivariate rational function field of rank 2 over Rational Field
            Variables: x, y
            sage: magma(ZZ['x'].fraction_field())                    # optional - magma
            Univariate rational function field over Integer Ring
            Variables: x

        Verify that conversion is being properly cached::

            sage: k = Frac(QQ['x,z'])                                # optional - magma
            sage: magma(k) is magma(k)                               # optional - magma
            True
        """
        s = 'FieldOfFractions(%s)' % self.ring()._magma_init_(magma)
        return magma._with_names(s, self.variable_names())

    def ring(self):
        """
        Return the ring that this is the fraction field of.

        EXAMPLES::

            sage: R = Frac(QQ['x,y'])
            sage: R
            Fraction Field of Multivariate Polynomial Ring in x, y over Rational Field
            sage: R.ring()
            Multivariate Polynomial Ring in x, y over Rational Field
        """
        return self._R

    @cached_method
    def is_exact(self):
        """
        Return if ``self`` is exact which is if the underlying ring is exact.

        EXAMPLES::

            sage: Frac(ZZ['x']).is_exact()
            True
            sage: Frac(CDF['x']).is_exact()                                             # needs sage.rings.complex_double
            False
        """
        return self.ring().is_exact()

    def _element_constructor_(self, x, y=None, coerce=True):
        """
        Construct an element of this fraction field.

        EXAMPLES::

            sage: F = QQ['x,y'].fraction_field()
            sage: F._element_constructor_(1)
            1
            sage: F._element_constructor_(F.gen(0)/F.gen(1))
            x/y
            sage: F._element_constructor_('1 + x/y')
            (x + y)/y

        ::

            sage: K = ZZ['x,y'].fraction_field()
            sage: x,y = K.gens()

        ::

            sage: F._element_constructor_(x/y)
            x/y

        TESTS:

        The next example failed before :issue:`4376`::

            sage: K(pari((x + 1)/(x^2 + x + 1)))                                        # needs sage.libs.pari
            (x + 1)/(x^2 + x + 1)

        These examples failed before :issue:`11368`::

            sage: R.<x, y, z> = PolynomialRing(QQ)
            sage: S = R.fraction_field()
            sage: S(pari((x + y)/y))                                                    # needs sage.libs.pari
            (x + y)/y

            sage: S(pari(x + y + 1/z))                                                  # needs sage.libs.pari
            (x*z + y*z + 1)/z

        This example failed before :issue:`23664`::

            sage: P0.<x> = ZZ[]
            sage: P1.<y> = Frac(P0)[]
            sage: frac = (x/(x^2 + 1))*y + 1/(x^3 + 1)
            sage: Frac(ZZ['x,y'])(frac)
            (x^4*y + x^2 + x*y + 1)/(x^5 + x^3 + x^2 + 1)

        Test conversions where `y` is a string but `x` not::

            sage: K = ZZ['x,y'].fraction_field()
            sage: K._element_constructor_(2, 'x+y')
            2/(x + y)
            sage: K._element_constructor_(1, 'z')
            Traceback (most recent call last):
            ...
            TypeError: unable to evaluate 'z' in Fraction Field of Multivariate Polynomial Ring in x, y over Integer Ring

        Check that :issue:`17971` is fixed::

            sage: A.<a,c> = Frac(PolynomialRing(QQ,'a,c'))
            sage: B.<d,e> = PolynomialRing(A,'d,e')
            sage: R.<x> = PolynomialRing(B,'x')
            sage: (a*d*x^2+a+e+1).resultant(-4*c^2*x+1)                                 # needs sage.libs.pari
            a*d + (16*c^4)*e + (16*a*c^4 + 16*c^4)

        Check that :issue:`24539` is fixed::

            sage: tau = polygen(QQ, 'tau')
            sage: PolynomialRing(CyclotomicField(2), 'z').fraction_field()(tau/(1+tau))             # needs sage.rings.number_field
            z/(z + 1)

        Check that :issue:`26150` is fixed::

            sage: z = SR.var('z')                                                       # needs sage.symbolic
            sage: CyclotomicField(2)['z'].fraction_field()(2*(4*z + 5)/((z + 1)*(z - 1)^4))         # needs sage.rings.number_field sage.symbolic
            (8*z + 10)/(z^5 - 3*z^4 + 2*z^3 + 2*z^2 - 3*z + 1)

        ::

            sage: T.<t> = ZZ[]
            sage: S.<s> = ZZ[]
            sage: S.fraction_field()(s/(s+1), (t-1)/(t+2))
            (s^2 + 2*s)/(s^2 - 1)

        Check that :issue:`29713` is fixed::

            sage: F = FractionField(QQ['a'])
            sage: a = F.gen()
            sage: R = PolynomialRing(F, 'x')
            sage: FF = FractionField(R)
            sage: elt = F(-1/2/(a^2+a))
            sage: x = FF(elt)
            sage: F(x)
            -1/2/(a^2 + a)
        """
        if isinstance(x, (list, tuple)) and len(x) == 1:
            x = x[0]
        if y is None:
            if parent(x) is self:
                return x
            ring_one = self.ring().one()
            try:
                return self._element_class(self, x, ring_one, coerce=coerce)
            except (TypeError, ValueError):
                pass
            y = self._element_class(self, ring_one, ring_one,
                                    coerce=False, reduce=False)
        else:
            if parent(x) is self:
                y = self(y)
                x, y = x.numerator() * y.denominator(), y.numerator() * x.denominator()
            try:
                return self._element_class(self, x, y, coerce=coerce)
            except (TypeError, ValueError):
                pass

        if isinstance(x, str):
            from sage.misc.sage_eval import sage_eval
            try:
                x = sage_eval(x, self.gens_dict_recursive())
            except NameError:
                raise TypeError("unable to evaluate {!r} in {}".format(x, self))
        if isinstance(y, str):
            from sage.misc.sage_eval import sage_eval
            try:
                y = sage_eval(y, self.gens_dict_recursive())
            except NameError:
                raise TypeError("unable to evaluate {!r} in {}".format(y, self))

        x = py_scalar_to_element(x)
        y = py_scalar_to_element(y)

        try:
            from cypari2.gen import Gen as pari_gen
        except ImportError:
            pari_gen = ()

        if isinstance(x, pari_gen) and x.type() == 't_POL':
            # This recursive approach is needed because PARI
            # represents multivariate polynomials as iterated
            # univariate polynomials (see the above examples).
            # Below, v is the variable with highest priority,
            # and the x[i] are rational functions in the
            # remaining variables.
            d = x.poldegree()
            if d.type() == 't_INFINITY':
                return self.zero()
            v = self._element_class(self, x.variable(), 1)
            x = sum(self(x[i]) * v**i for i in range(d + 1))

        def resolve_fractions(x, y):
            xn = x.numerator()
            xd = x.denominator()
            yn = y.numerator()
            yd = y.denominator()
            try:
                return (xn * yd, yn * xd)
            except (AttributeError, TypeError, ValueError):
                pass
            try:
                P = parent(yd)
                return (P(xn) * yd, yn * P(xd))
            except (AttributeError, TypeError, ValueError):
                pass
            try:
                P = parent(xd)
                return (xn * P(yd), P(yn) * xd)
            except (AttributeError, TypeError, ValueError):
                pass
            raise TypeError

        while True:
            x0, y0 = x, y
            try:
                x, y = resolve_fractions(x0, y0)
            except (AttributeError, TypeError):
                raise TypeError("cannot convert {!r}/{!r} to an element of {}".format(
                                x0, y0, self))
            try:
                return self._element_class(self, x, y, coerce=coerce)
            except TypeError:
                if parent(x) is parent(x0):
                    raise

    def construction(self):
        """
        EXAMPLES::

            sage: Frac(ZZ['x']).construction()
            (FractionField, Univariate Polynomial Ring in x over Integer Ring)
            sage: K = Frac(GF(3)['t'])
            sage: f, R = K.construction()
            sage: f(R)
            Fraction Field of Univariate Polynomial Ring in t
             over Finite Field of size 3
            sage: f(R) == K
            True
        """
        from sage.categories.pushout import FractionField
        return FractionField(), self.ring()

    def __eq__(self, other):
        """
        Check whether ``self`` is equal to ``other``.

        EXAMPLES::

            sage: Frac(ZZ['x']) == Frac(ZZ['x'])
            True
            sage: Frac(ZZ['x']) == Frac(QQ['x'])
            False
            sage: Frac(ZZ['x']) == Frac(ZZ['y'])
            False
            sage: Frac(ZZ['x']) == QQ['x']
            False
        """
        if not isinstance(other, FractionField_generic):
            return False
        return self._R == other._R

    def __ne__(self, other):
        """
        Check whether ``self`` is not equal to ``other``.

        EXAMPLES::

            sage: Frac(ZZ['x']) != Frac(ZZ['x'])
            False
            sage: Frac(ZZ['x']) != Frac(QQ['x'])
            True
            sage: Frac(ZZ['x']) != Frac(ZZ['y'])
            True
            sage: Frac(ZZ['x']) != QQ['x']
            True
        """
        return not (self == other)

    def __hash__(self):
        """
        Compute the hash of ``self``.

        EXAMPLES::

            sage: h0 = hash(Frac(ZZ['x']))
            sage: h1 = hash(Frac(ZZ['x']))
            sage: h2 = hash(Frac(QQ['x']))
            sage: h3 = hash(ZZ['x'])
            sage: h0 == h1 and h1 != h2 and h1 != h3
            True
        """
        # to avoid having exactly the same hash as the base ring,
        # we change this hash using a random number
        return hash(self._R) ^ 147068341996611

    def ngens(self):
        """
        This is the same as for the parent object.

        EXAMPLES::

            sage: R = Frac(PolynomialRing(QQ,'z',10)); R
            Fraction Field of Multivariate Polynomial Ring
             in z0, z1, z2, z3, z4, z5, z6, z7, z8, z9 over Rational Field
            sage: R.ngens()
            10
        """
        return self._R.ngens()

    def gen(self, i=0):
        """
        Return the ``i``-th generator of ``self``.

        EXAMPLES::

            sage: R = Frac(PolynomialRing(QQ,'z',10)); R
            Fraction Field of Multivariate Polynomial Ring
             in z0, z1, z2, z3, z4, z5, z6, z7, z8, z9 over Rational Field
            sage: R.0
            z0
            sage: R.gen(3)
            z3
            sage: R.3
            z3
        """
        x = self._R.gen(i)
        one = self._R.one()
        r = self._element_class(self, x, one, coerce=False, reduce=False)
        return r

    def _is_valid_homomorphism_(self, codomain, im_gens, base_map=None):
        """
        Check if the homomorphism defined by sending generators of this
        fraction field to ``im_gens`` in ``codomain`` is valid.

        EXAMPLES::

            sage: F = QQ['x,y'].fraction_field()
            sage: x,y = F.gens()
            sage: F._is_valid_homomorphism_(F, [y,x])
            True
            sage: R = ZZ['x']; x = R.gen()
            sage: F._is_valid_homomorphism_(R, [x, x])
            False

        TESTS::

            sage: F._is_valid_homomorphism_(ZZ, [])
            False

        Test homomorphisms::

            sage: phi = F.hom([2*y, x])
            sage: phi(x+y)
            x + 2*y
            sage: phi(x/y)
            2*y/x
        """
        if len(im_gens) != self.ngens():
            return False
        # It is very difficult to check that the image of any element
        # is invertible.  Checking that the image of each generator
        # is a unit is not sufficient.  So we just give up and check
        # that elements of the base ring coerce to the codomain
        if base_map is None and not codomain.has_coerce_map_from(self.base_ring()):
            return False
        return True

    def random_element(self, *args, **kwds):
        """
        Return a random element in this fraction field.

        The arguments are passed to the random generator of the underlying ring.

        EXAMPLES::

            sage: F = ZZ['x'].fraction_field()
            sage: F.random_element()  # random
            (2*x - 8)/(-x^2 + x)

        ::

            sage: f = F.random_element(degree=5)
            sage: f.numerator().degree() == f.denominator().degree()
            True
            sage: f.denominator().degree() <= 5
            True
            sage: while f.numerator().degree() != 5:
            ....:      f = F.random_element(degree=5)
        """
        return self._element_class(self, self._R.random_element(*args, **kwds),
                                   self._R._random_nonzero_element(*args, **kwds),
                                   coerce=False, reduce=True)

    def some_elements(self):
        r"""
        Return some elements in this field.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: R.fraction_field().some_elements()
            [0,
             1,
             x,
             2*x,
             x/(x^2 + 2*x + 1),
             1/x^2,
             ...
             (2*x^2 + 2)/(x^2 + 2*x + 1),
             (2*x^2 + 2)/x^3,
             (2*x^2 + 2)/(x^2 - 1),
             2]
        """
        ret = [self.zero(), self.one()]
        ret.extend(self(a) / self(b)
                   for a in self._R.some_elements()
                   for b in self._R.some_elements()
                   if a != b and self(a) and self(b))
        return ret

    def _gcd_univariate_polynomial(self, f, g):
        r"""
        Helper method used to compute polynomial gcds over this field.

        See :meth:`sage.rings.polynomial.polynomial_element.Polynomial.gcd`.

        TESTS::

            sage: A.<x,y> = ZZ[]
            sage: C.<z> = Frac(A)[]
            sage: c = (2*y^2 - 11*x - 2*y + 1)/(-x^2 + x*y - 2*y^2)
            sage: p = (c*z^2 + x^10*z + 1)^6
            sage: q = (z^2 + c*x^10*z + 1)^6
            sage: g = p.gcd(q)
            sage: g
            1
            sage: g.parent() is p.parent()
            True
            sage: (p*(z-x)).gcd(q*(z-x))
            z - x
            sage: C.zero().gcd(2*z)
            z
            sage: (x*z).gcd(0)
            z
            sage: C.zero().gcd(0)
            0
        """
        if g.is_zero():
            if f.is_zero():
                return f
            else:
                return f.monic()
        Pol = f.parent()
        Num = Pol.change_ring(self.base())
        f1 = Num(f.numerator())
        g1 = Num(g.numerator())
        return Pol(f1.gcd(g1)).monic()


class FractionField_1poly_field(FractionField_generic):
    """
    The fraction field of a univariate polynomial ring over a field.

    Many of the functions here are included for coherence with number fields.
    """
    def __init__(self, R,
                 element_class=fraction_field_element.FractionFieldElement_1poly_field):
        """
        Just change the default for ``element_class``.

        EXAMPLES::

            sage: R.<t> = QQ[]; K = R.fraction_field()
            sage: K._element_class
            <class 'sage.rings.fraction_field_element.FractionFieldElement_1poly_field'>
        """
        FractionField_generic.__init__(self, R, element_class)

    def ring_of_integers(self):
        """
        Return the ring of integers in this fraction field.

        EXAMPLES::

            sage: K = FractionField(GF(5)['t'])
            sage: K.ring_of_integers()
            Univariate Polynomial Ring in t over Finite Field of size 5
        """
        return self._R

    def maximal_order(self):
        """
        Return the maximal order in this fraction field.

        EXAMPLES::

            sage: K = FractionField(GF(5)['t'])
            sage: K.maximal_order()
            Univariate Polynomial Ring in t over Finite Field of size 5
        """
        return self._R

    def class_number(self):
        """
        Here for compatibility with number fields and function fields.

        EXAMPLES::

            sage: R.<t> = GF(5)[]; K = R.fraction_field()
            sage: K.class_number()
            1
        """
        return 1

    def _factor_univariate_polynomial(self, f):
        r"""
        Return the factorization of ``f`` over this field.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: k.<a> = GF(9)
            sage: K = k['t'].fraction_field()
            sage: R.<x> = K[]
            sage: f = x^3 + a
            sage: f.factor()
            (x + 2*a + 1)^3
        """
        # The default implementation would try to convert this element to singular and factor there.
        # This fails silently over some base fields, see #23642, so we convert
        # to the function field and factor there.
        return f.change_ring(self.function_field()).factor().base_change(f.parent())

    def function_field(self):
        r"""
        Return the isomorphic function field.

        EXAMPLES::

            sage: R.<t> = GF(5)[]
            sage: K = R.fraction_field()
            sage: K.function_field()
            Rational function field in t over Finite Field of size 5

        .. SEEALSO::

            :meth:`sage.rings.function_field.RationalFunctionField.field`
        """
        from sage.rings.function_field.constructor import FunctionField
        return FunctionField(self.base_ring(), names=self.variable_name())

    def _coerce_map_from_(self, R):
        r"""
        Return a coerce map from ``R`` to this field.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: R.<t> = GF(5)[]
            sage: K = R.fraction_field()
            sage: L = K.function_field()
            sage: f = K.coerce_map_from(L); f # indirect doctest
            Isomorphism:
              From: Rational function field in t over Finite Field of size 5
              To:   Fraction Field of Univariate Polynomial Ring in t
                     over Finite Field of size 5
            sage: f(~L.gen())
            1/t
        """
        from sage.rings.function_field.function_field_rational import (
            RationalFunctionField,
        )
        if isinstance(R, RationalFunctionField) and self.variable_name() == R.variable_name() and self.base_ring() is R.constant_base_field():
            from sage.categories.homset import Hom
            parent = Hom(R, self)
            from sage.rings.function_field.maps import FunctionFieldToFractionField
            return parent.__make_element_class__(FunctionFieldToFractionField)(parent)

        return super()._coerce_map_from_(R)


class FractionFieldEmbedding(DefaultConvertMap_unique):
    r"""
    The embedding of an integral domain into its field of fractions.

    EXAMPLES::

        sage: R.<x> = QQ[]
        sage: f = R.fraction_field().coerce_map_from(R); f
        Coercion map:
          From: Univariate Polynomial Ring in x over Rational Field
          To:   Fraction Field of Univariate Polynomial Ring in x over Rational Field

    TESTS::

        sage: from sage.rings.fraction_field import FractionFieldEmbedding
        sage: isinstance(f, FractionFieldEmbedding)
        True
        sage: TestSuite(f).run()

    Check that :issue:`23185` has been resolved::

        sage: R.<x> = QQ[]
        sage: K.<x> = FunctionField(QQ)
        sage: R.is_subring(K)
        True
        sage: R.is_subring(R.fraction_field())
        True
    """
    def is_surjective(self):
        r"""
        Return whether this map is surjective.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: R.fraction_field().coerce_map_from(R).is_surjective()
            False
        """
        return self.domain().is_field()

    def is_injective(self):
        r"""
        Return whether this map is injective.

        EXAMPLES:

        The map from an integral domain to its fraction field is always
        injective::

            sage: R.<x> = QQ[]
            sage: R.fraction_field().coerce_map_from(R).is_injective()
            True
        """
        return True

    def section(self):
        r"""
        Return a section of this map.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: R.fraction_field().coerce_map_from(R).section()
            Section map:
              From: Fraction Field of Univariate Polynomial Ring in x over Rational Field
              To:   Univariate Polynomial Ring in x over Rational Field
        """
        from sage.categories.homset import Hom
        from sage.categories.sets_with_partial_maps import SetsWithPartialMaps
        parent = Hom(self.codomain(), self.domain(), SetsWithPartialMaps())
        return parent.__make_element_class__(FractionFieldEmbeddingSection)(self)

    def _richcmp_(self, other, op):
        r"""
        Compare this element to ``other`` with respect to ``op``.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: f = R.fraction_field().coerce_map_from(R)
            sage: S.<y> = GF(2)[]
            sage: g = S.fraction_field().coerce_map_from(S)

            sage: f == g # indirect doctest
            False
            sage: f == f
            True
        """
        if type(self) is not type(other):
            return NotImplemented
        return richcmp((self.domain(), self.codomain()), (other.domain(), other.codomain()), op)

    def __hash__(self):
        r"""
        Return a hash value for this embedding.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: hash(R.fraction_field().coerce_map_from(R)) == hash(R.fraction_field().coerce_map_from(R))
            True
        """
        return hash((type(self), self.domain()))


class FractionFieldEmbeddingSection(Section):
    r"""
    The section of the embedding of an integral domain into its field of
    fractions.

    EXAMPLES::

        sage: R.<x> = QQ[]
        sage: f = R.fraction_field().coerce_map_from(R).section(); f
        Section map:
          From: Fraction Field of Univariate Polynomial Ring in x over Rational Field
          To:   Univariate Polynomial Ring in x over Rational Field

    TESTS::

        sage: from sage.rings.fraction_field import FractionFieldEmbeddingSection
        sage: isinstance(f, FractionFieldEmbeddingSection)
        True
        sage: TestSuite(f).run()
    """
    def _call_(self, x, check=True):
        r"""
        Evaluate this map at ``x``.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: K = R.fraction_field()
            sage: x = K.gen()
            sage: f = K.coerce_map_from(R).section()
            sage: f(x)
            x
            sage: f(1/x)
            Traceback (most recent call last):
            ...
            TypeError: fraction must have unit denominator

        TESTS:

        Over inexact rings, we have to take the precision of the denominators
        into account::

            sage: # needs sage.rings.padics
            sage: R = ZpCR(2)
            sage: S.<x> = R[]
            sage: f = x/S(R(3, absprec=2))
            sage: S(f)
            (1 + 2 + O(2^2))*x

        Test for Localization::

            sage: R.<x> = ZZ[]
            sage: L = Localization(R, x**2 + 2*x + 1)                                   # needs sage.libs.pari
            sage: 1/(x + 1) in L                        # indirect doctest              # needs sage.libs.pari
            True
            sage: 1/(x + 2) in L                        # indirect doctest              # needs sage.libs.pari
            False
        """
        codom = self.codomain()
        if self.domain()._R is codom:
            num = x.numerator()
            den = x.denominator()
        else:
            # codomain may different from the fraction fields base ring
            # for example for localizations
            num = codom(x.numerator())
            den = codom(x.denominator())

        if codom.is_exact() and den.is_one():
            return num
        if check and not den.is_unit():
            # This should probably be a ValueError.
            # However, too much existing code is expecting this to throw a
            # TypeError, so we decided to keep it for the time being.
            raise TypeError("fraction must have unit denominator")
        return num * den.inverse_of_unit()

    def _call_with_args(self, x, args=(), kwds={}):
        r"""
        Evaluate this map at ``x``.

        INPUT:

        - ``check`` -- whether or not to check

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: K = R.fraction_field()
            sage: R(K.gen(), check=True)
            x
        """
        check = kwds.get('check', True)
        if args or any(key != 'check' for key in kwds):
            raise NotImplementedError("__call__ cannot be called with additional arguments other than check=True/False")
        return self._call_(x, check=check)

    def _richcmp_(self, other, op):
        r"""
        Compare this element to ``other`` with respect to ``op``.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: f = R.fraction_field().coerce_map_from(R).section()
            sage: S.<y> = GF(2)[]
            sage: g = S.fraction_field().coerce_map_from(S).section()

            sage: f == g # indirect doctest
            False
            sage: f == f
            True
        """
        if type(self) is not type(other):
            return NotImplemented
        return richcmp((self.domain(), self.codomain()), (other.domain(), other.codomain()), op)

    def __hash__(self):
        r"""
        Return a hash value for this section.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: hash(R.fraction_field().coerce_map_from(R).section()) == hash(R.fraction_field().coerce_map_from(R).section())
            True
        """
        return hash((type(self), self.codomain()))
