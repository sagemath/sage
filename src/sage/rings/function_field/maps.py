# -*- coding: utf-8 -*-
r"""
Morphisms of function fields

Maps and morphisms useful for computations with function fields.

EXAMPLES::

    sage: K.<x> = FunctionField(QQ); R.<y> = K[]
    sage: K.hom(1/x)
    Function Field endomorphism of Rational function field in x over Rational Field
      Defn: x |--> 1/x
    sage: L.<y> = K.extension(y^2 - x)
    sage: K.hom(y)
    Function Field morphism:
      From: Rational function field in x over Rational Field
      To:   Function field in y defined by y^2 - x
      Defn: x |--> y
    sage: L.hom([y,x])
    Function Field endomorphism of Function field in y defined by y^2 - x
      Defn: y |--> y
            x |--> x
    sage: L.hom([x,y])
    Traceback (most recent call last):
    ...
    ValueError: invalid morphism

For global function fields, which have positive characteristics, the higher
derivation is available::

    sage: K.<x> = FunctionField(GF(2)); _.<Y>=K[]
    sage: L.<y> = K.extension(Y^3+x+x^3*Y)
    sage: h = L.higher_derivation()
    sage: h(y^2, 2)
    ((x^7 + 1)/x^2)*y^2 + x^3*y

AUTHORS:

- William Stein (2010): initial version

- Julian Rüth (2011-09-14, 2014-06-23, 2017-08-21): refactored class hierarchy; added
  derivation classes; morphisms to/from fraction fields

- Kwankyu Lee (2017-04-30): added higher derivations and completions

"""
# ****************************************************************************
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
#       Copyright (C) 2011-2017 Julian Rüth <julian.rueth@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.morphism import Morphism, SetMorphism
from sage.categories.map import Map
from sage.categories.homset import Hom
from sage.categories.sets_cat import Sets

from sage.rings.infinity import infinity
from sage.rings.morphism import RingHomomorphism

from sage.modules.free_module_element import vector

from sage.functions.other import binomial

from sage.matrix.constructor import matrix


class FunctionFieldDerivation(Map):
    r"""
    Base class for derivations on function fields.

    A derivation on `R` is a map `R \to R` with
    `D(\alpha+\beta)=D(\alpha)+D(\beta)` and `D(\alpha\beta)=\beta
    D(\alpha)+\alpha D(\beta)` for all `\alpha,\beta\in R`.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ)
        sage: d = K.derivation()
        sage: d
        Derivation map:
          From: Rational function field in x over Rational Field
          To:   Rational function field in x over Rational Field
          Defn: x |--> 1
    """
    def __init__(self, K):
        r"""
        Initialize a derivation from `K` to `K`.

        INPUT:

        - ``K`` -- function field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: d = K.derivation()
            sage: TestSuite(d).run(skip=['_test_category', '_test_pickling'])

        .. TODO::

            Make the caching done at the map by subclassing
            ``UniqueRepresentation``, which will then implement a
            valid equality check. Then this will pass the pickling test.
        """
        from .function_field import is_FunctionField
        if not is_FunctionField(K):
            raise ValueError("K must be a function field")
        self.__field = K
        Map.__init__(self, Hom(K,K,Sets()))

    def _repr_type(self):
        r"""
        Return the type of this map (a derivation), for the purposes of printing out self.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: d = K.derivation()
            sage: d._repr_type()
            'Derivation'

        """
        return "Derivation"

    def is_injective(self):
        r"""
        Return ``False`` since a derivation is never injective.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: d = K.derivation()
            sage: d.is_injective()
            False
        """
        return False


class FunctionFieldDerivation_rational(FunctionFieldDerivation):
    """
    Derivations on rational function fields.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ)
        sage: K.derivation()
        Derivation map:
          From: Rational function field in x over Rational Field
          To:   Rational function field in x over Rational Field
          Defn: x |--> 1
    """
    def __init__(self, K, u):
        """
        Initialize a derivation of ``K`` which sends the generator of ``K`` to ``u``.

        INPUT:

        - ``K`` -- rational function field

        - ``u`` -- element of ``K``; the image of the generator of K under the
          derivation

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: d = K.derivation()
            sage: TestSuite(d).run(skip=["_test_category", "_test_pickling"])

        See the comment about the test suite run in
        ``FunctionFieldDerivation.__init__``.
        """
        FunctionFieldDerivation.__init__(self, K)

        self._u = u

    def _call_(self, x):
        """
        Compute the derivation of ``x``.

        INPUT:

        - ``x`` -- element of the rational function field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: d = K.derivation()
            sage: d(x)  # indirect doctest
            1
            sage: d(x^3)
            3*x^2
            sage: d(1/x)
            -1/x^2
        """
        f = x.numerator()
        g = x.denominator()

        numerator = f.derivative() * g - f * g.derivative()
        if numerator.is_zero():
            return self.codomain().zero()
        else:
            return self._u * self.codomain()(numerator / g**2)

    def _repr_defn(self):
        r"""
        Helper method to print this map.

        TESTS::

            sage: K.<x> = FunctionField(QQ)
            sage: K.derivation()  # indirect doctest
            Derivation map:
              From: Rational function field in x over Rational Field
              To:   Rational function field in x over Rational Field
              Defn: x |--> 1

        """
        return "%s |--> %s"%(self.domain().variable_name(), self(self.domain().gen()))


class FunctionFieldDerivation_separable(FunctionFieldDerivation):
    """
    Derivations of separable extensions.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ)
        sage: R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x)
        sage: L.derivation()
        Derivation map:
          From: Function field in y defined by y^2 - x
          To:   Function field in y defined by y^2 - x
          Defn: y |--> 1/2/x*y
                x |--> 1
    """
    def __init__(self, L, d):
        """
        Initialize.

        INPUT:

        - ``L`` -- function field; a separable extension of the domain of ``d``

        - ``d`` -- derivation on the base function field of ``L``

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: d = L.derivation()
            sage: TestSuite(d).run(skip=["_test_category", "_test_pickling"])

        See the comment about the test suite run in
        ``FunctionFieldDerivation.__init__``.
        """
        FunctionFieldDerivation.__init__(self, L)

        self._d = d
        if not L.is_separable():
            raise ValueError("L must be a separable extension of its base field.")

        x = self.domain().gen()
        f = L.polynomial()
        self._d = d
        self._gen_image = - f.map_coefficients(d)(x) / f.derivative()(x)

    def _call_(self, x):
        r"""
        Evaluate the derivation on ``x``.

        INPUT:

        - ``x`` -- element of the function field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: d = L.derivation()
            sage: d(x) # indirect doctest
            1
            sage: d(y)
            1/2/x*y
            sage: d(y^2)
            1
        """
        if x.is_zero():
            return self.codomain().zero()

        x = x._x
        y = self.domain().gen()

        return x.map_coefficients(self._d) + x.derivative()(y) * self._gen_image

    def _repr_defn(self):
        """
        Return the string representation of the map.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: L.derivation() # indirect doctest
            Derivation map:
              From: Function field in y defined by y^2 - x
              To:   Function field in y defined by y^2 - x
              Defn: y |--> 1/2/x*y
                    x |--> 1

            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^2 - y)
            sage: M.derivation()
            Derivation map:
              From: Function field in z defined by z^2 - y
              To:   Function field in z defined by z^2 - y
              Defn: z |--> 1/4/x*z
                    y |--> 1/2/x*y
                    x |--> 1

        """
        base = self._d._repr_defn()
        ret = "%s |--> %s"%(self.domain().variable_name(),self._gen_image)
        if base:
            return ret + "\n" + base
        else:
            return ret


class FunctionFieldDerivation_inseparable(FunctionFieldDerivation):
    r"""
    A generator of the space of derivations on ``L``.

    INPUT:

    - ``L`` -- a function field which is an inseparable extension of its base
      field.

    EXAMPLES::

        sage: K.<x> = FunctionField(GF(2))
        sage: R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x)
        sage: d = L.derivation()

    This also works for iterated non-monic extensions::

        sage: K.<x> = FunctionField(GF(2))
        sage: R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - 1/x)
        sage: R.<z> = L[]
        sage: M.<z> = L.extension(z^2*y - x^3)
        sage: M.derivation()
        Derivation map:
          From: Function field in z defined by y*z^2 + x^3
          To:   Function field in z defined by y*z^2 + x^3
          Defn: z |--> 1
                y |--> 0
                x |--> 0

    """
    def __init__(self, L):
        r"""
        Initialization.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(3))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^3 - x)
            sage: d = L.derivation() # indirect doctest
            sage: type(d)
            <class 'sage.rings.function_field.maps.FunctionFieldDerivation_inseparable'>

        """
        from .function_field import is_FunctionField
        if not is_FunctionField(L):
            raise TypeError("L must be a function field")
        FunctionFieldDerivation.__init__(self, L)

        if L.is_separable():
            raise ValueError("L must be an inseparable extension of its base field.")
        M, self._f, self._t = L.separable_model()
        self._d = M.derivation()

    def _call_(self, x):
        r"""
        Evaluate the derivation on ``x``.

        INPUT:

        - ``x`` -- an element of the function field

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: d = L.derivation()
            sage: d(x) # indirect doctest
            0
            sage: d(y)
            1
            sage: d(y^2)
            0

        """
        if x.is_zero():
            return self.codomain().zero()
        return self._f(self._d(self._t(x)))

    def _repr_defn(self):
        r"""
        Helper method to print this map.

        TESTS::

            sage: K.<x> = FunctionField(GF(2))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: L.derivation() # indirect doctest
            Derivation map:
              From: Function field in y defined by y^2 + x
              To:   Function field in y defined by y^2 + x
              Defn: y |--> 1
                    x |--> 0
            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^2 - y)
            sage: M.derivation()
            Derivation map:
              From: Function field in z defined by z^2 + y
              To:   Function field in z defined by z^2 + y
              Defn: z |--> 1
                    y |--> 0
                    x |--> 0

        """
        ret = ["%s |--> %s"%(k.variable_name(), self(k.gen())) for k in self.domain()._intermediate_fields(self.domain().rational_function_field())]
        return "\n".join(ret)


class FunctionFieldHigherDerivation(Map):
    """
    Base class of higher derivations on function fields.

    INPUT:

    - ``field`` -- function field on which the derivation operates

    EXAMPLES::

        sage: F.<x> = FunctionField(GF(2))
        sage: F.higher_derivation()
        Higher derivation map:
          From: Rational function field in x over Finite Field of size 2
          To:   Rational function field in x over Finite Field of size 2
    """
    def __init__(self, field):
        """
        Initialize.

        TESTS::

            sage: F.<x> = FunctionField(GF(4))
            sage: h = F.higher_derivation()
            sage: TestSuite(h).run(skip='_test_category')
        """
        Map.__init__(self, Hom(field, field, Sets()))

        self._field = field

        # elements of a prime finite field do not have pth_root method
        if field.constant_base_field().is_prime_field():
            self._pth_root_func = _pth_root_in_prime_field
        else:
            self._pth_root_func = _pth_root_in_finite_field

    def _repr_type(self):
        """
        Return a string containing the type of the map.

        EXAMPLES::

            sage: F.<x> = FunctionField(GF(2))
            sage: h = F.higher_derivation()
            sage: h  # indirect doctest
            Higher derivation map:
              From: Rational function field in x over Finite Field of size 2
              To:   Rational function field in x over Finite Field of size 2
        """
        return 'Higher derivation'

    def __eq__(self, other):
        """
        Test if ``self`` equals ``other``.

        TESTS::

            sage: F.<x> = FunctionField(GF(2))
            sage: h = F.higher_derivation()
            sage: loads(dumps(h)) == h
            True
        """
        if isinstance(other, FunctionFieldHigherDerivation):
            return self._field == other._field
        return False


def _pth_root_in_prime_field(e):
    """
    Return the `p`-th root of element ``e`` in a prime finite field.

    TESTS::

        sage: from sage.rings.function_field.maps import _pth_root_in_prime_field
        sage: p = 5
        sage: F.<a> = GF(p)
        sage: e = F.random_element()
        sage: _pth_root_in_prime_field(e)^p == e
        True
    """
    return e

def _pth_root_in_finite_field(e):
    """
    Return the `p`-th root of element ``e`` in a finite field.

    TESTS::

        sage: from sage.rings.function_field.maps import _pth_root_in_finite_field
        sage: p = 3
        sage: F.<a> = GF(p^2)
        sage: e = F.random_element()
        sage: _pth_root_in_finite_field(e)^p == e
        True
    """
    return e.pth_root()


class RationalFunctionFieldHigherDerivation_global(FunctionFieldHigherDerivation):
    """
    Higher derivations of rational function fields over finite fields.

    INPUT:

    - ``field`` -- function field on which the derivation operates

    EXAMPLES::

        sage: F.<x> = FunctionField(GF(2))
        sage: h = F.higher_derivation()
        sage: h
        Higher derivation map:
          From: Rational function field in x over Finite Field of size 2
          To:   Rational function field in x over Finite Field of size 2
        sage: h(x^2,2)
        1
    """
    def __init__(self, field):
        """
        Initialize.

        TESTS::

            sage: F.<x> = FunctionField(GF(2))
            sage: h = F.higher_derivation()
            sage: TestSuite(h).run(skip='_test_category')
        """
        FunctionFieldHigherDerivation.__init__(self, field)

        self._p = field.characteristic()
        self._separating_element = field.gen()

    def _call_with_args(self, f, args=(), kwds={}):
        """
        Call the derivation with args and kwds.

        EXAMPLES::

            sage: F.<x> = FunctionField(GF(2))
            sage: h = F.higher_derivation()
            sage: h(x^2,2)  # indirect doctest
            1
        """
        return self._derive(f, *args, **kwds)

    def _derive(self, f, i, separating_element=None):
        """
        Return the `i`-th derivative of ``f`` with respect to the
        separating element.

        This implements Hess' Algorithm 26 in [Hes2002b]_.

        EXAMPLES::

            sage: F.<x> = FunctionField(GF(2))
            sage: h = F.higher_derivation()
            sage: h._derive(x^3,0)
            x^3
            sage: h._derive(x^3,1)
            x^2
            sage: h._derive(x^3,2)
            x
            sage: h._derive(x^3,3)
            1
            sage: h._derive(x^3,4)
            0
        """
        F = self._field
        p = self._p

        if separating_element is None:
            x = self._separating_element
            derivative = lambda f: f.derivative()
        else:
            x = separating_element
            xderinv = ~(x.derivative())
            derivative = lambda f: xderinv * f.derivative()

        prime_power_representation = self._prime_power_representation

        def derive(f, i):
            # Step 1: zero-th derivative
            if i == 0:
                return f
            # Step 2:
            s = i % p
            r = i // p
            # Step 3:
            e = f
            while s > 0:
                e = derivative(e) / F(s)
                s -= 1
            # Step 4:
            if r == 0:
                return e
            else:
                # Step 5:
                lambdas = prime_power_representation(e, x)
                # Step 6 and 7:
                der = 0
                for i in range(p):
                    mu = derive(lambdas[i], r)
                    der += mu**p * x**i
                return der

        return derive(f, i)

    def _prime_power_representation(self, f, separating_element=None):
        """
        Return `p`-th power representation of the element ``f``.

        Here `p` is the characteristic of the function field.

        This implements Hess' Algorithm 25.

        EXAMPLES::

            sage: F.<x> = FunctionField(GF(2))
            sage: h = F.higher_derivation()
            sage: h._prime_power_representation(x^2 + x + 1)
            [x + 1, 1]
            sage: x^2 + x + 1 == _[0]^2 + _[1]^2 * x
            True
        """
        F = self._field
        p = self._p

        if separating_element is None:
            x = self._separating_element
            derivative = lambda f: f.derivative()
        else:
            x = separating_element
            xderinv = ~(x.derivative())
            derivative = lambda f: xderinv * f.derivative()

        # Step 1:
        a = [f]
        aprev = f
        j = 1
        while j < p:
            aprev = derivative(aprev) / F(j)
            a.append(aprev)
            j += 1
        # Step 2:
        b = a
        j = p - 2
        while j >= 0:
            b[j] -= sum(binomial(i,j) * b[i] * x**(i-j) for i in range(j+1, p))
            j -= 1
        # Step 3
        return [self._pth_root(c) for c in b]

    def _pth_root(self, c):
        """
        Return the `p`-th root of the rational function ``c``.

        INPUT:

        - ``c`` -- rational function

        EXAMPLES::

            sage: F.<x> = FunctionField(GF(2))
            sage: h = F.higher_derivation()
            sage: h._pth_root((x^2+1)^2)
            x^2 + 1
        """
        K = self._field
        p = self._p

        R = K._field.ring()

        poly = c.numerator()
        num = R([self._pth_root_func(poly[i]) for i in range(0, poly.degree()+1, p)])
        poly = c.denominator()
        den = R([self._pth_root_func(poly[i]) for i in range(0, poly.degree()+1, p)])
        return K.element_class(K, num/den)


class FunctionFieldHigherDerivation_global(FunctionFieldHigherDerivation):
    """
    Higher derivations of global function fields.

    INPUT:

    - ``field`` -- function field on which the derivation operates

    EXAMPLES::

        sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
        sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
        sage: h = L.higher_derivation()
        sage: h
        Higher derivation map:
          From: Function field in y defined by y^3 + x^3*y + x
          To:   Function field in y defined by y^3 + x^3*y + x
        sage: h(y^2, 2)
        ((x^7 + 1)/x^2)*y^2 + x^3*y
    """

    def __init__(self, field):
        """
        Initialize.

        TESTS::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
            sage: h = L.higher_derivation()
            sage: TestSuite(h).run(skip=['_test_category'])
        """
        FunctionFieldHigherDerivation.__init__(self, field)

        self._p = field.characteristic()
        self._separating_element = field(field.base_field().gen())

        p = field.characteristic()
        y = field.gen()

        # matrix for pth power map; used in _prime_power_representation method
        self.__pth_root_matrix = matrix([(y**(i*p)).list()
                                         for i in range(field.degree())]).transpose()

        # cache computed higher derivatives to speed up later computations
        self._cache = {}

    def _call_with_args(self, f, args, kwds):
        """
        Call the derivation with ``args`` and ``kwds``.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
            sage: h = L.higher_derivation()
            sage: h(y^2,2)  # indirect doctest
            ((x^7 + 1)/x^2)*y^2 + x^3*y
        """
        return self._derive(f, *args, **kwds)

    def _derive(self, f, i, separating_element=None):
        """
        Return ``i``-th derivative of ``f`` with respect to the separating
        element.

        This implements Hess' Algorithm 26 in [Hes2002b]_.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
            sage: h = L.higher_derivation()
            sage: y^3
            x^3*y + x
            sage: h._derive(y^3,0)
            x^3*y + x
            sage: h._derive(y^3,1)
            x^4*y^2 + 1
            sage: h._derive(y^3,2)
            x^10*y^2 + (x^8 + x)*y
            sage: h._derive(y^3,3)
            (x^9 + x^2)*y^2 + x^7*y
            sage: h._derive(y^3,4)
            (x^22 + x)*y^2 + ((x^21 + x^14 + x^7 + 1)/x)*y
        """
        F = self._field
        p = self._p
        frob = F.frobenius_endomorphism() # p-th power map

        if separating_element is None:
            x = self._separating_element
            derivative = lambda f: f.derivative()
        else:
            x = separating_element
            xderinv = ~(x.derivative())
            derivative = lambda f: xderinv * f.derivative()

        try:
            cache = self._cache[separating_element]
        except KeyError:
            cache = self._cache[separating_element] = {}

        def derive(f, i):
            # Step 1: zero-th derivative
            if i == 0:
                return f

            # Step 1.5: use cached result if available
            try:
                return cache[f,i]
            except KeyError:
                pass

            # Step 2:
            s = i % p
            r = i // p
            # Step 3:
            e = f
            while s > 0:
                e = derivative(e) / F(s)
                s -= 1
            # Step 4:
            if r == 0:
                der = e
            else:
                # Step 5: inlined self._prime_power_representation
                a = [e]
                aprev = e
                j = 1
                while j < p:
                    aprev = derivative(aprev) / F(j)
                    a.append(aprev)
                    j += 1
                b = a
                j = p - 2
                while j >= 0:
                    b[j] -= sum(binomial(k,j) * b[k] * x**(k-j) for k in range(j+1, p))
                    j -= 1
                lambdas = [self._pth_root(c) for c in b]

                # Step 6 and 7:
                der = 0
                xpow = 1
                for k in range(p):
                    mu = derive(lambdas[k], r)
                    der += frob(mu) * xpow
                    xpow *= x

            cache[f,i] = der
            return der

        return derive(f, i)

    def _prime_power_representation(self, f, separating_element=None):
        """
        Return `p`-th power representation of the element ``f``.

        Here `p` is the characteristic of the function field.

        This implements Hess' Algorithm 25 in [Hes2002b]_.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
            sage: h = L.higher_derivation()
            sage: b = h._prime_power_representation(y)
            sage: y == b[0]^2 + b[1]^2 * x
            True
        """
        F = self._field
        p = self._p

        if separating_element is None:
            x = self._separating_element
            derivative = lambda f: f.derivative()
        else:
            x = separating_element
            xderinv = ~(x.derivative())
            derivative = lambda f: xderinv * f.derivative()

        # Step 1:
        a = [f]
        aprev = f
        j = 1
        while j < p:
            aprev = derivative(aprev) / F(j)
            a.append(aprev)
            j += 1
        # Step 2:
        b = a
        j = p - 2
        while j >= 0:
            b[j] -= sum(binomial(i,j) * b[i] * x**(i-j) for i in range(j+1,p))
            j -= 1
        # Step 3
        return [self._pth_root(c) for c in b]

    def _pth_root(self, c):
        """
        Return the `p`-th root of function field element ``c``.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
            sage: h = L.higher_derivation()
            sage: h._pth_root((x^2 + y^2)^2)
            y^2 + x^2
        """
        K = self._field.base_field() # rational function field
        p = self._p

        coeffs = []
        for d in self.__pth_root_matrix.solve_right(vector(c.list())):
            poly = d.numerator()
            num = K([self._pth_root_func(poly[i]) for i in range(0, poly.degree()+1, p)])
            poly = d.denominator()
            den = K([self._pth_root_func(poly[i]) for i in range(0, poly.degree()+1, p)])
            coeffs.append(num / den)
        return self._field(coeffs)


class FunctionFieldHigherDerivation_char_zero(FunctionFieldHigherDerivation):
    """
    Higher derivations of function fields of characteristic zero.

    INPUT:

    - ``field`` -- function field on which the derivation operates

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ); _.<Y> = K[]
        sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
        sage: h = L.higher_derivation()
        sage: h
        Higher derivation map:
          From: Function field in y defined by y^3 + x^3*y + x
          To:   Function field in y defined by y^3 + x^3*y + x
        sage: h(y,1) == -(3*x^2*y+1)/(3*y^2+x^3)
        True
        sage: h(y^2,1) == -2*y*(3*x^2*y+1)/(3*y^2+x^3)
        True
        sage: e = L.random_element()
        sage: h(h(e,1),1) == 2*h(e,2)
        True
        sage: h(h(h(e,1),1),1) == 3*2*h(e,3)
        True
    """
    def __init__(self, field):
        """
        Initialize.

        TESTS::

            sage: K.<x> = FunctionField(QQ); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
            sage: h = L.higher_derivation()
            sage: TestSuite(h).run(skip=['_test_category'])
        """
        FunctionFieldHigherDerivation.__init__(self, field)

        self._separating_element = field(field.base_field().gen())

        # cache computed higher derivatives to speed up later computations
        self._cache = {}

    def _call_with_args(self, f, args, kwds):
        """
        Call the derivation with ``args`` and ``kwds``.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
            sage: h = L.higher_derivation()
            sage: e = L.random_element()
            sage: h(h(e,1),1) == 2*h(e,2)  # indirect doctest
            True
        """
        return self._derive(f, *args, **kwds)

    def _derive(self, f, i, separating_element=None):
        """
        Return ``i``-th derivative of ``f`` with respect to the separating
        element.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^3 + x + x^3*Y)
            sage: h = L.higher_derivation()
            sage: y^3
            -x^3*y - x
            sage: h._derive(y^3,0)
            -x^3*y - x
            sage: h._derive(y^3,1)
            (-21/4*x^4/(x^7 + 27/4))*y^2 + ((-9/2*x^9 - 45/2*x^2)/(x^7 + 27/4))*y + (-9/2*x^7 - 27/4)/(x^7 + 27/4)
        """
        F = self._field

        if separating_element is None:
            x = self._separating_element
            xderinv = 1
        else:
            x = separating_element
            xderinv = ~(x.derivative())

        try:
            cache = self._cache[separating_element]
        except KeyError:
            cache = self._cache[separating_element] = {}

        if i == 0:
            return f

        try:
            return cache[f,i]
        except KeyError:
            pass

        s = i
        e = f
        while s > 0:
            e = xderinv * e.derivative() / F(s)
            s -= 1

        der = e

        cache[f,i] = der
        return der


class FunctionFieldVectorSpaceIsomorphism(Morphism):
    r"""
    Base class for isomorphisms between function fields and vector spaces.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ); R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
        sage: V, f, t = L.vector_space()
        sage: isinstance(f, sage.rings.function_field.maps.FunctionFieldVectorSpaceIsomorphism)
        True
    """
    def _repr_(self):
        """
        Return the string representation of this isomorphism.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()
            sage: f
            Isomorphism:
              From: Vector space of dimension 2 over Rational function field in x over Rational Field
              To:   Function field in y defined by y^2 - x*y + 4*x^3
            sage: t
            Isomorphism:
              From: Function field in y defined by y^2 - x*y + 4*x^3
              To:   Vector space of dimension 2 over Rational function field in x over Rational Field
        """
        s = "Isomorphism:"
        s += "\n  From: {}".format(self.domain())
        s += "\n  To:   {}".format(self.codomain())
        return s

    def is_injective(self):
        """
        Return ``True``, since the isomorphism is injective.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()
            sage: f.is_injective()
            True
        """
        return True

    def is_surjective(self):
        """
        Return ``True``, since the isomorphism is surjective.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()
            sage: f.is_surjective()
            True
        """
        return True

    def _richcmp_(self, other, op):
        r"""
        Compare this map to ``other``.

        .. NOTE::

            This implementation assumes that this isomorphism is defined by its
            domain and codomain. Isomorphisms for which this is not true must
            override this implementation.

        EXAMPLES::

            sage: K = QQ['x'].fraction_field()
            sage: L = K.function_field()
            sage: f = K.coerce_map_from(L)

            sage: K = QQbar['x'].fraction_field()
            sage: L = K.function_field()
            sage: g = K.coerce_map_from(L)

            sage: f == g
            False
            sage: f == f
            True

        """
        if type(self) != type(other):
            return NotImplemented

        from sage.structure.richcmp import richcmp
        return richcmp((self.domain(),self.codomain()), (other.domain(),other.codomain()), op)

    def __hash__(self):
        r"""
        Return a hash value of this map.

        This implementation assumes that this isomorphism is defined by its
        domain and codomain. Isomorphisms for which this is not true should
        override this implementation.

        EXAMPLES::

            sage: K = QQ['x'].fraction_field()
            sage: L = K.function_field()
            sage: f = K.coerce_map_from(L)
            sage: hash(f) == hash(f)
            True

        """
        return hash((self.domain(), self.codomain()))


class MapVectorSpaceToFunctionField(FunctionFieldVectorSpaceIsomorphism):
    """
    Isomorphism from a vector space to a function field.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ); R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
        sage: V, f, t = L.vector_space(); f
        Isomorphism:
          From: Vector space of dimension 2 over Rational function field in x over Rational Field
          To:   Function field in y defined by y^2 - x*y + 4*x^3
    """
    def __init__(self, V, K):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space(); type(f)
            <class 'sage.rings.function_field.maps.MapVectorSpaceToFunctionField'>
        """
        self._V = V
        self._K = K
        self._R = K.polynomial_ring()
        FunctionFieldVectorSpaceIsomorphism.__init__(self, Hom(V, K))

    def _call_(self, v):
        """
        Map ``v`` to the function field.

        INPUT:

        - ``v`` -- element of the vector space

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()
            sage: f(x*V.0 + (1/x^3)*V.1) # indirect doctest
            1/x^3*y + x

        TESTS:

        Test that this map is a bijection for some random inputs::

            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^3 - y - x)
            sage: for F in [K,L,M]:
            ....:     for base in F._intermediate_fields(K):
            ....:         V, f, t = F.vector_space(base)
            ....:         for i in range(100):
            ....:             a = F.random_element()
            ....:             assert(f(t(a)) == a)

        """
        fields = self._K._intermediate_fields(self._V.base_field())
        fields.pop()
        degrees = [k.degree() for k in fields]
        gens = [k.gen() for k in fields]

        # construct the basis composed of powers of the generators of all the
        # intermediate fields, i.e., 1, x, y, x*y, ...
        from sage.misc.misc_c import prod
        from itertools import product
        exponents = product(*[range(d) for d in degrees])
        basis = [prod(g**e for g,e in zip(gens,es)) for es in exponents]

        # multiply the entries of v with the values in basis
        coefficients = self._V(v).list()
        ret = sum([c*b for (c,b) in zip(coefficients,basis)])
        return self._K(ret)

    def domain(self):
        """
        Return the vector space which is the domain of the isomorphism.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()
            sage: f.domain()
            Vector space of dimension 2 over Rational function field in x over Rational Field
        """
        return self._V

    def codomain(self):
        """
        Return the function field which is the codomain of the isomorphism.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()
            sage: f.codomain()
            Function field in y defined by y^2 - x*y + 4*x^3
        """
        return self._K


class MapFunctionFieldToVectorSpace(FunctionFieldVectorSpaceIsomorphism):
    """
    Isomorphism from a function field to a vector space.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ); R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
        sage: V, f, t = L.vector_space(); t
        Isomorphism:
          From: Function field in y defined by y^2 - x*y + 4*x^3
          To:   Vector space of dimension 2 over Rational function field in x over Rational Field
    """
    def __init__(self, K, V):
        """
        Initialize.

        INPUT:

        - ``K`` -- function field

        - ``V`` -- vector space isomorphic to the function field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()
            sage: TestSuite(t).run(skip="_test_category")
        """
        self._V = V
        self._K = K
        self._zero = K.base_ring()(0)
        self._n = K.degree()
        FunctionFieldVectorSpaceIsomorphism.__init__(self, Hom(K, V))

    def _call_(self, x):
        """
        Map ``x`` to the vector space.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ); R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x*y + 4*x^3)
            sage: V, f, t = L.vector_space()
            sage: t(x + (1/x^3)*y) # indirect doctest
            (x, 1/x^3)

        TESTS:

        Test that this map is a bijection for some random inputs::

            sage: R.<z> = L[]
            sage: M.<z> = L.extension(z^3 - y - x)
            sage: for F in [K,L,M]:
            ....:     for base in F._intermediate_fields(K):
            ....:         V, f, t = F.vector_space(base)
            ....:         for i in range(100):
            ....:             a = V.random_element()
            ....:             assert(t(f(a)) == a)

        """
        ret = [x]
        fields = self._K._intermediate_fields(self._V.base_field())
        fields.pop()
        from itertools import chain
        for k in fields:
            ret = chain.from_iterable([y.list() for y in ret])
        ret = list(ret)
        assert all(t.parent() is self._V.base_field() for t in ret)
        return self._V(ret)


class FunctionFieldMorphism(RingHomomorphism):
    """
    Base class for morphisms between function fields.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ)
        sage: f = K.hom(1/x); f
        Function Field endomorphism of Rational function field in x over Rational Field
          Defn: x |--> 1/x
    """
    def __init__(self, parent, im_gen, base_morphism):
        """
        Initialize.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: f = K.hom(1/x); f
            Function Field endomorphism of Rational function field in x over Rational Field
              Defn: x |--> 1/x
            sage: TestSuite(f).run(skip="_test_category")
        """
        RingHomomorphism.__init__(self, parent)

        self._im_gen = im_gen
        self._base_morphism = base_morphism

    def _repr_type(self):
        r"""
        Return the type of the morphism for the purpose of printing.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^3 + 6*x^3 + x)
            sage: f = L.hom(y*2)
            sage: f._repr_type()
            'Function Field'

        """
        return "Function Field"

    def _repr_defn(self):
        """
        Return the string containing the definition of the morphism.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^3 + 6*x^3 + x)
            sage: f = L.hom(y*2)
            sage: f._repr_defn()
            'y |--> 2*y'
        """
        a = '%s |--> %s'%(self.domain().variable_name(), self._im_gen)
        if self._base_morphism is not None:
            a += '\n' + self._base_morphism._repr_defn()
        return a


class FunctionFieldMorphism_polymod(FunctionFieldMorphism):
    """
    Morphism from a finite extension of a function field to a function field.

    EXAMPLES::

        sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
        sage: L.<y> = K.extension(y^3 + 6*x^3 + x)
        sage: f = L.hom(y*2); f
        Function Field endomorphism of Function field in y defined by y^3 + 6*x^3 + x
          Defn: y |--> 2*y
        sage: factor(L.polynomial())
        y^3 + 6*x^3 + x
        sage: f(y).charpoly('y')
        y^3 + 6*x^3 + x
    """
    def __init__(self, parent, im_gen, base_morphism):
        """
        Initialize.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^3 + 6*x^3 + x)
            sage: f = L.hom(y*2)
            sage: TestSuite(f).run(skip="_test_category")
        """
        FunctionFieldMorphism.__init__(self, parent, im_gen, base_morphism)
        # Verify that the morphism is valid:
        R = self.codomain()['X']
        v = parent.domain().polynomial().list()
        if base_morphism is not None:
            v = [base_morphism(a) for a in v]
        f = R(v)
        if f(im_gen):
            raise ValueError("invalid morphism")

    def _call_(self, x):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7)); R.<y> = K[]
            sage: L.<y> = K.extension(y^3 + 6*x^3 + x); f = L.hom(y*2)
            sage: f(y/x + x^2/(x+1))            # indirect doctest
            2/x*y + x^2/(x + 1)
            sage: f(y)
            2*y
        """
        v = x.list()
        if self._base_morphism is not None:
            v = [self._base_morphism(a) for a in v]
        f = v[0].parent()['X'](v)
        return f(self._im_gen)


class FunctionFieldMorphism_rational(FunctionFieldMorphism):
    """
    Morphism from a rational function field to a function field.
    """
    def __init__(self, parent, im_gen, base_morphism):
        """
        Initialize.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7))
            sage: f = K.hom(1/x); f
            Function Field endomorphism of Rational function field in x over Finite Field of size 7
              Defn: x |--> 1/x
        """
        FunctionFieldMorphism.__init__(self, parent, im_gen, base_morphism)

    def _call_(self, x):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(GF(7))
            sage: f = K.hom(1/x); f
            Function Field endomorphism of Rational function field in x over Finite Field of size 7
              Defn: x |--> 1/x
            sage: f(x+1)                          # indirect doctest
            (x + 1)/x
            sage: 1/x + 1
            (x + 1)/x

        You can specify a morphism on the base ring::

            sage: Qi = GaussianIntegers().fraction_field()
            sage: i = Qi.gen()
            sage: K.<x> = FunctionField(Qi)
            sage: phi1 = Qi.hom([CC.gen()])
            sage: phi2 = Qi.hom([-CC.gen()])
            sage: f = K.hom(CC.pi(),phi1)
            sage: f(1+i+x)
            4.14159265358979 + 1.00000000000000*I
            sage: g = K.hom(CC.pi(),phi2)
            sage: g(1+i+x)
            4.14159265358979 - 1.00000000000000*I
        """
        a = x.element()
        if self._base_morphism is None:
            return a.subs({a.parent().gen():self._im_gen})
        else:
            f = self._base_morphism
            num = a.numerator()
            den = a.denominator()
            R = self._im_gen.parent()['X']
            num = R([f(c) for c in num.list()])
            den = R([f(c) for c in den.list()])
            return num.subs(self._im_gen) / den.subs(self._im_gen)


class FunctionFieldConversionToConstantBaseField(Map):
    r"""
    Conversion map from the function field to its constant base field.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ)
        sage: QQ.convert_map_from(K)
        Conversion map:
          From: Rational function field in x over Rational Field
          To:   Rational Field
    """
    def __init__(self, parent):
        """
        Initialize.

        TESTS::

            sage: K.<x> = FunctionField(QQ)
            sage: f = QQ.convert_map_from(K)
            sage: from sage.rings.function_field.maps import FunctionFieldConversionToConstantBaseField
            sage: isinstance(f, FunctionFieldConversionToConstantBaseField)
            True
        """
        Map.__init__(self, parent)

    def _repr_type(self):
        r"""
        Return the type of this map (a conversion), for the purposes of printing out self.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: QQ.convert_map_from(K) # indirect doctest
            Conversion map:
              From: Rational function field in x over Rational Field
              To:   Rational Field

        """
        return "Conversion"

    def _call_(self, x):
        """
        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: QQ(K(1)) # indirect doctest
            1

        """
        return x.parent()._to_constant_base_field(x)


class FunctionFieldToFractionField(FunctionFieldVectorSpaceIsomorphism):
    r"""
    Isomorphism from rational function field to the isomorphic fraction
    field of a polynomial ring.

    EXAMPLES::

        sage: K = QQ['x'].fraction_field()
        sage: L = K.function_field()
        sage: f = K.coerce_map_from(L); f
        Isomorphism:
          From: Rational function field in x over Rational Field
          To:   Fraction Field of Univariate Polynomial Ring in x over Rational Field

    .. SEEALSO::

        :class:`FractionFieldToFunctionField`

    TESTS::

        sage: from sage.rings.function_field.maps import FunctionFieldToFractionField
        sage: isinstance(f, FunctionFieldToFractionField)
        True
        sage: TestSuite(f).run()

    """
    def _call_(self, f):
        r"""
        Return the value of this map at ``f``.

        EXAMPLES::

            sage: K = QQ['x'].fraction_field()
            sage: L = K.function_field()
            sage: f = K.coerce_map_from(L)
            sage: f(~L.gen())
            1/x

        """
        return self.codomain()(f.numerator(), f.denominator())

    def section(self):
        r"""
        Return the inverse map of this isomorphism.

        EXAMPLES::

            sage: K = QQ['x'].fraction_field()
            sage: L = K.function_field()
            sage: f = K.coerce_map_from(L)
            sage: f.section()
            Isomorphism:
                From: Fraction Field of Univariate Polynomial Ring in x over Rational Field
                To:   Rational function field in x over Rational Field


        """
        parent = Hom(self.codomain(), self.domain())
        return parent.__make_element_class__(FractionFieldToFunctionField)(parent.domain(), parent.codomain())


class FractionFieldToFunctionField(FunctionFieldVectorSpaceIsomorphism):
    r"""
    Isomorphism from a fraction field of a polynomial ring to the isomorphic
    function field.

    EXAMPLES::

        sage: K = QQ['x'].fraction_field()
        sage: L = K.function_field()
        sage: f = L.coerce_map_from(K); f
        Isomorphism:
            From: Fraction Field of Univariate Polynomial Ring in x over Rational Field
            To:   Rational function field in x over Rational Field

    .. SEEALSO::

        :class:`FunctionFieldToFractionField`

    TESTS::

        sage: from sage.rings.function_field.maps import FractionFieldToFunctionField
        sage: isinstance(f, FractionFieldToFunctionField)
        True
        sage: TestSuite(f).run()

    """
    def _call_(self, f):
        r"""
        Return the value of this morphism at ``f``.

        EXAMPLES::

            sage: K = QQ['x'].fraction_field()
            sage: L = K.function_field()
            sage: f = L.coerce_map_from(K)
            sage: f(~K.gen())
            1/x
        """
        return self.codomain()._element_constructor_(f)

    def section(self):
        r"""
        Return the inverse map of this isomorphism.

        EXAMPLES::

            sage: K = QQ['x'].fraction_field()
            sage: L = K.function_field()
            sage: f = L.coerce_map_from(K)
            sage: f.section()
            Isomorphism:
                From: Rational function field in x over Rational Field
                To:   Fraction Field of Univariate Polynomial Ring in x over Rational Field

        """
        parent = Hom(self.codomain(), self.domain())
        return parent.__make_element_class__(FunctionFieldToFractionField)(parent)


class FunctionFieldCompletion(Map):
    """
    Completions on function fields.

    INPUT:

    - ``field`` -- function field

    - ``place`` -- place of the function field

    - ``name`` -- string for the name of the series variable

    - ``prec`` -- positive integer; default precision

    - ``gen_name`` -- string; name of the generator of the residue
      field; used only when place is non-rational

    EXAMPLES::

        sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
        sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
        sage: p = L.places_finite()[0]
        sage: m = L.completion(p)
        sage: m
        Completion map:
          From: Function field in y defined by y^2 + y + (x^2 + 1)/x
          To:   Laurent Series Ring in s over Finite Field of size 2
        sage: m(x)
        s^2 + s^3 + s^4 + s^5 + s^7 + s^8 + s^9 + s^10 + s^12 + s^13
        + s^15 + s^16 + s^17 + s^19 + O(s^22)
        sage: m(y)
        s^-1 + 1 + s^3 + s^5 + s^7 + s^9 + s^13 + s^15 + s^17 + O(s^19)
        sage: m(x*y) == m(x) * m(y)
        True
        sage: m(x+y) == m(x) + m(y)
        True

    The variable name of the series can be supplied. If the place is not
    rational such that the residue field is a proper extension of the constant
    field, you can also specify the generator name of the extension::

        sage: p2 = L.places_finite(2)[0]
        sage: p2
        Place (x^2 + x + 1, x*y + 1)
        sage: m2 = L.completion(p2, 't', gen_name='b')
        sage: m2(x)
        (b + 1) + t + t^2 + t^4 + t^8 + t^16 + O(t^20)
        sage: m2(y)
        b + b*t + b*t^3 + b*t^4 + (b + 1)*t^5 + (b + 1)*t^7 + b*t^9 + b*t^11
        + b*t^12 + b*t^13 + b*t^15 + b*t^16 + (b + 1)*t^17 + (b + 1)*t^19 + O(t^20)
    """
    def __init__(self, field, place, name=None, prec=None, gen_name=None):
        """
        Initialize.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: p = L.places_finite()[0]
            sage: m = L.completion(p)
            sage: m
            Completion map:
              From: Function field in y defined by y^2 + y + (x^2 + 1)/x
              To:   Laurent Series Ring in s over Finite Field of size 2
        """
        if name is None:
            name = 's' # default

        if gen_name is None:
            gen_name = 'a' # default

        k, from_k, to_k = place.residue_field(name=gen_name)

        self._place = place
        self._gen_name = gen_name

        if prec == infinity:
            from sage.rings.lazy_series_ring import LazyLaurentSeriesRing
            codomain = LazyLaurentSeriesRing(k, name)
            self._precision = infinity
        else: # prec < infinity:
            # if prec is None, the Laurent series ring provides default precision
            from sage.rings.laurent_series_ring import LaurentSeriesRing
            codomain = LaurentSeriesRing(k, name=name, default_prec=prec)
            self._precision = codomain.default_prec()

        Map.__init__(self, field, codomain)

    def _repr_type(self):
        """
        Return a string containing the type of the map.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: p = L.places_finite()[0]
            sage: m = L.completion(p)
            sage: m  # indirect doctest
            Completion map:
              From: Function field in y defined by y^2 + y + (x^2 + 1)/x
              To:   Laurent Series Ring in s over Finite Field of size 2
        """
        return 'Completion'

    def _call_(self, f):
        """
        Call the completion for f

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: p = L.places_finite()[0]
            sage: m = L.completion(p)
            sage: m(y)
            s^-1 + 1 + s^3 + s^5 + s^7 + s^9 + s^13 + s^15 + s^17 + O(s^19)
        """
        if self._precision == infinity:
            return self._expand_lazy(f)
        else:
            return self._expand(f, prec=None)

    def _call_with_args(self, f, args, kwds):
        """
        Call the completion with ``args`` and ``kwds``.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: p = L.places_finite()[0]
            sage: m = L.completion(p)
            sage: m(x+y, 10)  # indirect doctest
            s^-1 + 1 + s^2 + s^4 + s^8 + O(s^9)
        """
        if self._precision == infinity:
            return self._expand_lazy(f, *args, **kwds)
        else:
            return self._expand(f, *args, **kwds)

    def _expand(self, f, prec=None):
        """
        Return the laurent series expansion of f with precision ``prec``.

        INPUT:

        - ``f`` -- element of the function field

        - ``prec`` -- positive integer; relative precision of the series

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: p = L.places_finite()[0]
            sage: m = L.completion(p)
            sage: m(x, prec=20)  # indirect doctest
            s^2 + s^3 + s^4 + s^5 + s^7 + s^8 + s^9 + s^10 + s^12 + s^13 + s^15
            + s^16 + s^17 + s^19 + O(s^22)
        """
        if prec is None:
            prec = self._precision

        place = self._place
        F = place.function_field()
        der = F.higher_derivation()

        k,from_k,to_k = place.residue_field(name=self._gen_name)
        sep = place.local_uniformizer()

        val = f.valuation(place)
        e = f * sep **(-val)

        coeffs = [to_k(der._derive(e, i, sep)) for i in range(prec)]
        return self.codomain()(coeffs, val).add_bigoh(prec + val)

    def _expand_lazy(self, f):
        """
        Return the lazy laurent series expansion of ``f``.

        INPUT:

        - ``f`` -- element of the function field

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: p = L.places_finite()[0]
            sage: m = L.completion(p, prec=infinity)
            sage: e = m(x); e
            s^2 + s^3 + s^4 + s^5 + s^7 + s^8 + ...
            sage: e.coefficient(99)  # indirect doctest
            0
            sage: e.coefficient(100)
            1
        """
        place = self._place
        F = place.function_field()
        der = F.higher_derivation()

        k, from_k, to_k = place.residue_field(name=self._gen_name)
        sep = place.local_uniformizer()

        val = f.valuation(place)
        e = f * sep**(-val)

        def coeff(s, n):
            return to_k(der._derive(e, n - val, sep))

        return self.codomain().series(coeff, valuation=val)

    def default_precision(self):
        """
        Return the default precision.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: p = L.places_finite()[0]
            sage: m = L.completion(p)
            sage: m.default_precision()
            20
        """
        return self._precision


class FunctionFieldRingMorphism(SetMorphism):
    """
    Ring homomorphism.
    """
    def _repr_(self):
        """
        Return the string representation of the map.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: p = L.places_finite()[0]
            sage: R = p.valuation_ring()
            sage: k, fr_k, to_k = R.residue_field()
            sage: k
            Finite Field of size 2
            sage: fr_k
            Ring morphism:
              From: Finite Field of size 2
              To:   Valuation ring at Place (x, x*y)
        """
        s = "Ring morphism:"
        s += "\n  From: {}".format(self.domain())
        s += "\n  To:   {}".format(self.codomain())
        return s


class FunctionFieldLinearMap(SetMorphism):
    """
    Linear map to function fields.
    """
    def _repr_(self):
        """
        Return the string representation of the map.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(5)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^2-x^3-1)
            sage: O = F.maximal_order()
            sage: I = O.ideal(x-2)
            sage: D = I.divisor()
            sage: V, from_V, to_V = D.function_space()
            sage: from_V
            Linear map:
              From: Vector space of dimension 2 over Finite Field of size 5
              To:   Function field in y defined by y^2 + 4*x^3 + 4
        """
        s = "Linear map:"
        s += "\n  From: {}".format(self.domain())
        s += "\n  To:   {}".format(self.codomain())
        return s


class FunctionFieldLinearMapSection(SetMorphism):
    """
    Section of linear map from function fields.
    """
    def _repr_(self):
        """
        Return the string representation of the map.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(5)); R.<t> = PolynomialRing(K)
            sage: F.<y> = K.extension(t^2-x^3-1)
            sage: O = F.maximal_order()
            sage: I = O.ideal(x-2)
            sage: D = I.divisor()
            sage: V, from_V, to_V = D.function_space()
            sage: to_V
            Section of linear map:
              From: Function field in y defined by y^2 + 4*x^3 + 4
              To:   Vector space of dimension 2 over Finite Field of size 5
        """
        s = "Section of linear map:"
        s += "\n  From: {}".format(self.domain())
        s += "\n  To:   {}".format(self.codomain())
        return s

