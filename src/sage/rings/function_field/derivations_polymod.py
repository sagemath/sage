r"""
Derivations of function fields: extension
"""

# ****************************************************************************
#       Copyright (C) 2010      William Stein <wstein@gmail.com>
#                     2011-2017 Julian Rüth <julian.rueth@gmail.com>
#                     2017      Alyson Deines
#                     2017-2019 Kwankyu Lee
#                     2018-2019 Travis Scrimshaw
#                     2019      Brent Baccala
#                     2022      Xavier Caruso
#                     2022      Frédéric Chapoton
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.arith.misc import binomial
from sage.categories.homset import Hom
from sage.categories.map import Map
from sage.categories.sets_cat import Sets
from sage.rings.derivation import RingDerivationWithoutTwist

from .derivations import FunctionFieldDerivation


class FunctionFieldDerivation_separable(FunctionFieldDerivation):
    """
    Derivations of separable extensions.

    EXAMPLES::

        sage: K.<x> = FunctionField(QQ)
        sage: R.<y> = K[]
        sage: L.<y> = K.extension(y^2 - x)
        sage: L.derivation()
        d/dx
    """
    def __init__(self, parent, d):
        """
        Initialize a derivation.

        INPUT:

        - ``parent`` -- the parent of this derivation

        - ``d`` -- a variable name or a derivation over
          the base field (or something capable to create
          such a derivation)

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: d = L.derivation()
            sage: TestSuite(d).run()

            sage: L.derivation(y)  # d/dy
            2*y*d/dx

            sage: dK = K.derivation([x]); dK
            x*d/dx
            sage: L.derivation(dK)
            x*d/dx
        """
        FunctionFieldDerivation.__init__(self, parent)
        L = parent.domain()
        C = parent.codomain()
        dm = parent._defining_morphism
        u = L.gen()
        if d == L.gen():
            d = parent._base_derivation(None)
            f = L.polynomial().change_ring(L)
            coeff = -f.derivative()(u) / f.map_coefficients(d)(u)
            if dm is not None:
                coeff = dm(coeff)
            self._d = parent._base_derivation([coeff])
            self._gen_image = C.one()
        else:
            if isinstance(d, RingDerivationWithoutTwist) and d.domain() is L.base_ring():
                self._d = d
            else:
                self._d = d = parent._base_derivation(d)
            f = L.polynomial()
            if dm is None:
                denom = f.derivative()(u)
            else:
                u = dm(u)
                denom = f.derivative().map_coefficients(dm, new_base_ring=C)(u)
            num = f.map_coefficients(d, new_base_ring=C)(u)
            self._gen_image = -num / denom

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
        parent = self.parent()
        if x.is_zero():
            return parent.codomain().zero()
        x = x._x
        y = parent.domain().gen()
        dm = parent._defining_morphism
        tmp1 = x.map_coefficients(self._d, new_base_ring=parent.codomain())
        tmp2 = x.derivative()(y)
        if dm is not None:
            tmp2 = dm(tmp2)
            y = dm(y)
        return tmp1(y) + tmp2 * self._gen_image

    def _add_(self, other):
        """
        Return the sum of this derivation and ``other``.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: d = L.derivation()
            sage: d
            d/dx
            sage: d + d
            2*d/dx
        """
        return type(self)(self.parent(), self._d + other._d)

    def _lmul_(self, factor):
        """
        Return the product of this derivation by the scalar ``factor``.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ)
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: d = L.derivation()
            sage: d
            d/dx
            sage: y * d
            y*d/dx
        """
        return type(self)(self.parent(), factor * self._d)


class FunctionFieldDerivation_inseparable(FunctionFieldDerivation):
    def __init__(self, parent, u=None):
        r"""
        Initialize this derivation.

        INPUT:

        - ``parent`` -- the parent of this derivation

        - ``u`` -- a parameter describing the derivation

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
            d/dz

        We can also create a multiple of the canonical derivation::

            sage: M.derivation([x])
            x*d/dz
        """
        FunctionFieldDerivation.__init__(self, parent)
        if u is None:
            self._u = parent.codomain().one()
        elif u == 0 or isinstance(u, (list, tuple)):
            if u == 0 or len(u) == 0:
                self._u = parent.codomain().zero()
            elif len(u) == 1:
                self._u = parent.codomain()(u[0])
            else:
                raise ValueError("the length does not match")
        else:
            raise ValueError("you must pass in either a name of a variable or a list of coefficients")

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
        parent = self.parent()
        return self._u * parent._d(parent._t(x))

    def _add_(self, other):
        """
        Return the sum of this derivation and ``other``.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(3))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^3 - x)
            sage: d = L.derivation()
            sage: d
            d/dy
            sage: d + d
            2*d/dy
        """
        return type(self)(self.parent(), [self._u + other._u])

    def _lmul_(self, factor):
        """
        Return the product of this derivation by the scalar ``factor``.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2))
            sage: R.<y> = K[]
            sage: L.<y> = K.extension(y^2 - x)
            sage: d = L.derivation()
            sage: d
            d/dy
            sage: y * d
            y*d/dy
        """
        return type(self)(self.parent(), [factor * self._u])


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

            sage: F.<x> = FunctionField(GF(4))                                                      # needs sage.rings.finite_rings
            sage: h = F.higher_derivation()                                                         # needs sage.rings.finite_rings
            sage: TestSuite(h).run(skip='_test_category')                                           # needs sage.rings.finite_rings
        """
        Map.__init__(self, Hom(field, field, Sets()))
        self._field = field
        # elements of a prime finite field do not have pth_root method
        if field.constant_base_field().is_prime_field():
            self._pth_root_func = _pth_root_in_prime_field
        else:
            self._pth_root_func = _pth_root_in_finite_field

    def _repr_type(self) -> str:
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

    def __eq__(self, other) -> bool:
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

        sage: from sage.rings.function_field.derivations_polymod import _pth_root_in_prime_field
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

        sage: from sage.rings.function_field.derivations_polymod import _pth_root_in_finite_field
        sage: p = 3
        sage: F.<a> = GF(p^2)                                                                       # needs sage.rings.finite_rings
        sage: e = F.random_element()                                                                # needs sage.rings.finite_rings
        sage: _pth_root_in_finite_field(e)^p == e                                                   # needs sage.rings.finite_rings
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
        sage: h(x^2, 2)
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
            sage: h(x^2, 2)  # indirect doctest
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
            sage: h._derive(x^3, 0)
            x^3
            sage: h._derive(x^3, 1)
            x^2
            sage: h._derive(x^3, 2)
            x
            sage: h._derive(x^3, 3)
            1
            sage: h._derive(x^3, 4)
            0
        """
        F = self._field
        p = self._p

        if separating_element is None:
            x = self._separating_element

            def derivative(f):
                return f.derivative()
        else:
            x = separating_element
            xderinv = ~(x.derivative())

            def derivative(f):
                return xderinv * f.derivative()

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

            def derivative(f):
                return f.derivative()
        else:
            x = separating_element
            xderinv = ~(x.derivative())

            def derivative(f):
                return xderinv * f.derivative()

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
            b[j] -= sum(binomial(i, j) * b[i] * x**(i - j)
                        for i in range(j + 1, p))
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
        num = R([self._pth_root_func(poly[i])
                 for i in range(0, poly.degree() + 1, p)])
        poly = c.denominator()
        den = R([self._pth_root_func(poly[i])
                 for i in range(0, poly.degree() + 1, p)])
        return K.element_class(K, num / den)


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
        from sage.matrix.constructor import matrix

        FunctionFieldHigherDerivation.__init__(self, field)

        self._p = field.characteristic()
        self._separating_element = field(field.base_field().gen())

        p = field.characteristic()
        y = field.gen()

        # matrix for pth power map; used in _prime_power_representation method
        self.__pth_root_matrix = matrix([(y**(i * p)).list()
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
            sage: h(y^2, 2)  # indirect doctest
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
            sage: h._derive(y^3, 0)
            x^3*y + x
            sage: h._derive(y^3, 1)
            x^4*y^2 + 1
            sage: h._derive(y^3, 2)
            x^10*y^2 + (x^8 + x)*y
            sage: h._derive(y^3, 3)
            (x^9 + x^2)*y^2 + x^7*y
            sage: h._derive(y^3, 4)
            (x^22 + x)*y^2 + ((x^21 + x^14 + x^7 + 1)/x)*y
        """
        F = self._field
        p = self._p
        frob = F.frobenius_endomorphism()  # p-th power map

        if separating_element is None:
            x = self._separating_element

            def derivative(f):
                return f.derivative()
        else:
            x = separating_element
            xderinv = ~(x.derivative())

            def derivative(f):
                return xderinv * f.derivative()

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
                return cache[f, i]
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
                    b[j] -= sum(binomial(k, j) * b[k] * x**(k - j)
                                for k in range(j + 1, p))
                    j -= 1
                lambdas = [self._pth_root(c) for c in b]

                # Step 6 and 7:
                der = 0
                xpow = 1
                for k in range(p):
                    mu = derive(lambdas[k], r)
                    der += frob(mu) * xpow
                    xpow *= x

            cache[f, i] = der
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

            def derivative(f):
                return f.derivative()
        else:
            x = separating_element
            xderinv = ~(x.derivative())

            def derivative(f):
                return xderinv * f.derivative()

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
            b[j] -= sum(binomial(i, j) * b[i] * x**(i - j)
                        for i in range(j + 1, p))
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
        from sage.modules.free_module_element import vector

        K = self._field.base_field()  # rational function field
        p = self._p

        coeffs = []
        for d in self.__pth_root_matrix.solve_right(vector(c.list())):
            poly = d.numerator()
            num = K([self._pth_root_func(poly[i])
                     for i in range(0, poly.degree() + 1, p)])
            poly = d.denominator()
            den = K([self._pth_root_func(poly[i])
                     for i in range(0, poly.degree() + 1, p)])
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
            sage: h._derive(y^3, 0)
            -x^3*y - x
            sage: h._derive(y^3, 1)
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
            return cache[f, i]
        except KeyError:
            pass

        s = i
        e = f
        while s > 0:
            e = xderinv * e.derivative() / F(s)
            s -= 1

        der = e

        cache[f, i] = der
        return der
