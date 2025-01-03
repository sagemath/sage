# sage.doctest: optional - sage.rings.finite_rings
r"""
Drinfeld modules over rings of characteristic zero

This module provides the class
:class:`sage.rings.function_fields.drinfeld_module.charzero_drinfeld_module.DrinfeldModule_charzero`,
which inherits
:class:`sage.rings.function_fields.drinfeld_module.drinfeld_module.DrinfeldModule`.

AUTHORS:

- David Ayotte (2023-09)
"""

# *****************************************************************************
#        Copyright (C) 2022 David Ayotte <david.ayotte@mail.concordia.ca>
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#                   http://www.gnu.org/licenses/
# *****************************************************************************

from .drinfeld_module import DrinfeldModule

from sage.rings.integer_ring import ZZ
from sage.rings.infinity import Infinity

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import

lazy_import('sage.rings.lazy_series_ring', 'LazyPowerSeriesRing')
lazy_import('sage.rings.power_series_ring', 'PowerSeriesRing')


class DrinfeldModule_charzero(DrinfeldModule):
    r"""
    This class implements Drinfeld `\mathbb{F}_q[T]`-modules defined
    over fields of `\mathbb{F}_q[T]`-characteristic zero.

    Recall that the `\mathbb{F}_q[T]`-*characteristic* is defined as the
    kernel of the underlying structure morphism. For general definitions
    and help on Drinfeld modules, see class
    :class:`sage.rings.function_fields.drinfeld_module.drinfeld_module.DrinfeldModule`.

    .. RUBRIC:: Construction:

    The user does not ever need to directly call
    ``DrinfeldModule_charzero`` --- the metaclass ``DrinfeldModule`` is
    responsible for instantiating the right class depending on the
    input::

        sage: A = GF(3)['T']
        sage: K.<T> = Frac(A)
        sage: phi = DrinfeldModule(A, [T, 1])
        sage: phi
        Drinfeld module defined by T |--> t + T

    ::

        sage: isinstance(phi, DrinfeldModule)
        True
        sage: from sage.rings.function_field.drinfeld_modules.charzero_drinfeld_module import DrinfeldModule_charzero
        sage: isinstance(phi, DrinfeldModule_charzero)
        True

    .. RUBRIC:: Logarithm and exponential

    It is possible to calculate the logarithm and the exponential of
    any Drinfeld modules of characteristic zero::

        sage: A = GF(2)['T']
        sage: K.<T> = Frac(A)
        sage: phi = DrinfeldModule(A, [T, 1])
        sage: phi.exponential()
        z + ((1/(T^2+T))*z^2) + ((1/(T^8+T^6+T^5+T^3))*z^4) + O(z^8)
        sage: phi.logarithm()
        z + ((1/(T^2+T))*z^2) + ((1/(T^6+T^5+T^3+T^2))*z^4) + O(z^8)

    .. RUBRIC:: Goss polynomials

    Goss polynomials are a sequence of polynomials related with the
    analytic theory of Drinfeld module. They provide a function field
    analogue of certain classical trigonometric functions::

        sage: A = GF(2)['T']
        sage: K.<T> = Frac(A)
        sage: phi = DrinfeldModule(A, [T, 1])
        sage: phi.goss_polynomial(1)
        X
        sage: phi.goss_polynomial(2)
        X^2
        sage: phi.goss_polynomial(3)
        X^3 + (1/(T^2 + T))*X^2

    .. RUBRIC:: Base fields of `\mathbb{F}_q[T]`-characteristic zero

    The base fields need not only be fraction fields of polynomials
    ring. In the following example, we construct a Drinfeld module over
    `\mathbb{F}_q((1/T))`, the completion of the rational function field
    at the place `1/T`::

        sage: A.<T> = GF(2)[]
        sage: L.<s> = LaurentSeriesRing(GF(2))  # s = 1/T
        sage: phi = DrinfeldModule(A, [1/s, s + s^2 + s^5 + O(s^6), 1+1/s])
        sage: phi(T)
        (s^-1 + 1)*t^2 + (s + s^2 + s^5 + O(s^6))*t + s^-1

    One can also construct Drinfeld modules over SageMath's global
    function fields::

        sage: A.<T> = GF(5)[]
        sage: K.<z> = FunctionField(GF(5))  # z = T
        sage: phi = DrinfeldModule(A, [z, 1, z^2])
        sage: phi(T)
        z^2*t^2 + t + z
    """
    @cached_method
    def _compute_coefficient_exp(self, k):
        r"""
        Return the `q^k`-th coefficient of the exponential of this
        Drinfeld module.

        INPUT:

        - ``k`` -- integer; the index of the coefficient

        TESTS::

            sage: A = GF(2)['T']
            sage: K.<T> = Frac(A)
            sage: phi = DrinfeldModule(A, [T, 1])
            sage: q = A.base_ring().cardinality()
            sage: phi._compute_coefficient_exp(0)
            1
            sage: phi._compute_coefficient_exp(1)
            1/(T^2 + T)
            sage: phi._compute_coefficient_exp(2)
            1/(T^8 + T^6 + T^5 + T^3)
            sage: phi._compute_coefficient_exp(3)
            1/(T^24 + T^20 + T^18 + T^17 + T^14 + T^13 + T^11 + T^7)
        """
        k = ZZ(k)
        if k.is_zero():
            return self._base.one()
        q = self._Fq.cardinality()
        c = self._base.zero()
        for i in range(k):
            j = k - i
            c += self._compute_coefficient_exp(i)*self._compute_coefficient_log(j)**(q**i)
        return -c

    def exponential(self, prec=Infinity, name='z'):
        r"""
        Return the exponential of this Drinfeld module.

        Note that the exponential is only defined when the
        `\mathbb{F}_q[T]`-characteristic is zero.

        INPUT:

        - ``prec`` -- an integer or ``Infinity`` (default: ``Infinity``);
          the precision at which the series is returned; if ``Infinity``,
          a lazy power series in returned

        - ``name`` -- string (default: ``'z'``); the name of the
          generator of the lazy power series ring

        EXAMPLES::

            sage: A = GF(2)['T']
            sage: K.<T> = Frac(A)
            sage: phi = DrinfeldModule(A, [T, 1])
            sage: q = A.base_ring().cardinality()

        When ``prec`` is ``Infinity`` (which is the default),
        the exponential is returned as a lazy power series, meaning
        that any of its coefficients can be computed on demands::

            sage: exp = phi.exponential(); exp
            z + ((1/(T^2+T))*z^2) + ((1/(T^8+T^6+T^5+T^3))*z^4) + O(z^8)
            sage: exp[2^4]
            1/(T^64 + T^56 + T^52 + ... + T^27 + T^23 + T^15)
            sage: exp[2^5]
            1/(T^160 + T^144 + T^136 + ... + T^55 + T^47 + T^31)

        On the contrary, when ``prec`` is a finite number, all the
        required coefficients are computed at once::

            sage: phi.exponential(prec=10)
            z + (1/(T^2 + T))*z^2 + (1/(T^8 + T^6 + T^5 + T^3))*z^4 + (1/(T^24 + T^20 + T^18 + T^17 + T^14 + T^13 + T^11 + T^7))*z^8 + O(z^10)

        Example in higher rank::

            sage: A = GF(5)['T']
            sage: K.<T> = Frac(A)
            sage: phi = DrinfeldModule(A, [T, T^2, T + T^2 + T^4, 1])
            sage: exp = phi.exponential(); exp
            z + ((T/(T^4+4))*z^5) + O(z^8)

        The exponential is the compositional inverse of the logarithm
        (see :meth:`logarithm`)::

            sage: log = phi.logarithm(); log
            z + ((4*T/(T^4+4))*z^5) + O(z^8)
            sage: exp.compose(log)
            z + O(z^8)
            sage: log.compose(exp)
            z + O(z^8)

        TESTS::

            sage: A = GF(2)['T']
            sage: K.<T> = Frac(A)
            sage: phi = DrinfeldModule(A, [T, 1])
            sage: exp = phi.exponential()
            sage: exp[2] == 1/(T**q - T)  # expected value
            True
            sage: exp[2^2] == 1/((T**(q**2) - T)*(T**q - T)**q)  # expected value
            True
            sage: exp[2^3] == 1/((T**(q**3) - T)*(T**(q**2) - T)**q*(T**q - T)**(q**2))  # expected value
            True

        REFERENCE:

        See section 4.6 of [Gos1998]_ for the definition of the
        exponential.
        """
        zero = self._base.zero()
        q = self._Fq.cardinality()

        def coeff_exp(k):
            # Return the k-th coefficient of the exponential.
            k = ZZ(k)
            v, u = k.val_unit(q)
            if u == 1:
                return self._compute_coefficient_exp(v)
            else:
                return zero

        if prec is Infinity:
            L = LazyPowerSeriesRing(self._base, name)
            return L(coeff_exp, valuation=1)
        L = PowerSeriesRing(self._base, name, default_prec=prec)
        return L([0] + [coeff_exp(i) for i in range(1,prec)], prec=prec)

    @cached_method
    def _compute_coefficient_log(self, k):
        r"""
        Return the `q^k`-th coefficient of the logarithm of this
        Drinfeld module.

        TESTS::

            sage: A = GF(2)['T']
            sage: K.<T> = Frac(A)
            sage: phi = DrinfeldModule(A, [T, 1])
            sage: q = A.base_ring().cardinality()
            sage: phi._compute_coefficient_log(0)
            1
            sage: phi._compute_coefficient_log(1)
            1/(T^2 + T)
            sage: phi._compute_coefficient_log(2)
            1/(T^6 + T^5 + T^3 + T^2)
            sage: phi._compute_coefficient_log(3)
            1/(T^14 + T^13 + T^11 + T^10 + T^7 + T^6 + T^4 + T^3)
        """
        k = ZZ(k)
        if k.is_zero():
            return self._base.one()
        r = self._gen.degree()
        T = self._gen[0]
        q = self._Fq.cardinality()
        c = self._base.zero()
        for i in range(k):
            j = k - i
            if j < r + 1:
                c += self._compute_coefficient_log(i)*self._gen[j]**(q**i)
        return c/(T - T**(q**k))

    def logarithm(self, prec=Infinity, name='z'):
        r"""
        Return the logarithm of the given Drinfeld module.

        By definition, the logarithm is the compositional inverse of the
        exponential (see :meth:`exponential`). Note that the logarithm
        is only defined when the `\mathbb{F}_q[T]`-characteristic is
        zero.

        INPUT:

        - ``prec`` -- an integer or ``Infinity`` (default: ``Infinity``);
          the precision at which the series is returned; if ``Infinity``,
          a lazy power series in returned

        - ``name`` -- string (default: ``'z'``); the name of the
          generator of the lazy power series ring

        EXAMPLES::

            sage: A = GF(2)['T']
            sage: K.<T> = Frac(A)
            sage: phi = DrinfeldModule(A, [T, 1])

        When ``prec`` is ``Infinity`` (which is the default),
        the logarithm is returned as a lazy power series, meaning
        that any of its coefficients can be computed on demands::

            sage: log = phi.logarithm(); log
            z + ((1/(T^2+T))*z^2) + ((1/(T^6+T^5+T^3+T^2))*z^4) + O(z^8)
            sage: log[2^4]
            1/(T^30 + T^29 + T^27 + ... + T^7 + T^5 + T^4)
            sage: log[2^5]
            1/(T^62 + T^61 + T^59 + ... + T^8 + T^6 + T^5)

        If ``prec`` is a finite number, all the
        required coefficients are computed at once::

            sage: phi.logarithm(prec=10)
            z + (1/(T^2 + T))*z^2 + (1/(T^6 + T^5 + T^3 + T^2))*z^4 + (1/(T^14 + T^13 + T^11 + T^10 + T^7 + T^6 + T^4 + T^3))*z^8 + O(z^10)

        Example in higher rank::

            sage: A = GF(5)['T']
            sage: K.<T> = Frac(A)
            sage: phi = DrinfeldModule(A, [T, T^2, T + T^2 + T^4, 1])
            sage: phi.logarithm()
            z + ((4*T/(T^4+4))*z^5) + O(z^8)

        TESTS::

            sage: A = GF(2)['T']
            sage: K.<T> = Frac(A)
            sage: phi = DrinfeldModule(A, [T, 1])
            sage: q = 2
            sage: log[2] == -1/((T**q - T))  # expected value
            True
            sage: log[2**2] == 1/((T**q - T)*(T**(q**2) - T))  # expected value
            True
            sage: log[2**3] == -1/((T**q - T)*(T**(q**2) - T)*(T**(q**3) - T))  # expected value
            True
        """
        q = self._Fq.cardinality()

        def coeff_log(k):
            # Return the k-th coefficient of the logarithm
            k = ZZ(k)
            v, u = k.val_unit(q)
            if u == 1:
                return self._compute_coefficient_log(v)
            else:
                return self._base.zero()

        if prec is Infinity:
            L = LazyPowerSeriesRing(self._base, name)
            return L(coeff_log, valuation=1)
        L = PowerSeriesRing(self._base, name, default_prec=prec)
        return L([0] + [coeff_log(i) for i in range(1, prec)], prec=prec)

    @cached_method
    def _compute_goss_polynomial(self, n, q, poly_ring, X):
        r"""
        Utility function for computing the n-th Goss polynomial.

        The user should not call this method directly, but
        :meth:`goss_polynomial` instead.

        TESTS::

            sage: A = GF(2^2)['T']
            sage: K.<T> = Frac(A)
            sage: phi = DrinfeldModule(A, [T, T+1, T^2, 1])
            sage: poly_ring = phi.base()['X']
            sage: X = poly_ring.gen()
            sage: phi._compute_goss_polynomial(0, 2^2, poly_ring, X)
            0
            sage: phi._compute_goss_polynomial(3, 2^2, poly_ring, X)
            X^3
            sage: phi._compute_goss_polynomial(4*3, 2^2, poly_ring, X)
            X^12
            sage: phi._compute_goss_polynomial(9, 2^2, poly_ring, X)
            X^9 + (1/(T^3 + T^2 + T))*X^6 + (1/(T^6 + T^4 + T^2))*X^3
        """
        # Trivial cases
        if n.is_zero():
            return poly_ring.zero()
        if n <= q - 1:
            return X**n
        if n % q == 0:
            return self.goss_polynomial(n // q)**q
        # General case
        pol = poly_ring.zero()
        m = q
        i = 1
        while m < n:
            pol += self._compute_coefficient_exp(i) * self._compute_goss_polynomial(n - m, q, poly_ring, X)
            m *= q
            i += 1
        return X*(self._compute_goss_polynomial(n - 1, q, poly_ring, X) + pol)

    def goss_polynomial(self, n, var='X'):
        r"""
        Return the `n`-th Goss polynomial of the Drinfeld module.

        Note that Goss polynomials are only defined for Drinfeld modules
        of characteristic zero.

        INPUT:

        - ``n`` -- integer; the index of the Goss polynomial

        - ``var``-- string (default: ``'X'``); the name of polynomial
          variable

        OUTPUT: a univariate polynomial in ``var`` over the base `A`-field

        EXAMPLES::

            sage: A = GF(3)['T']
            sage: K.<T> = Frac(A)
            sage: phi = DrinfeldModule(A, [T, 1])  # Carlitz module
            sage: phi.goss_polynomial(1)
            X
            sage: phi.goss_polynomial(2)
            X^2
            sage: phi.goss_polynomial(4)
            X^4 + (1/(T^3 + 2*T))*X^2
            sage: phi.goss_polynomial(5)
            X^5 + (2/(T^3 + 2*T))*X^3
            sage: phi.goss_polynomial(10)
            X^10 + (1/(T^3 + 2*T))*X^8 + (1/(T^6 + T^4 + T^2))*X^6 + (1/(T^9 + 2*T^3))*X^4 + (1/(T^18 + 2*T^12 + 2*T^10 + T^4))*X^2

        REFERENCE:

        Section 3 of [Gek1988]_ provides an exposition of Goss
        polynomials.
        """
        n = ZZ(n)
        K = self.base()
        poly_ring = K[var]
        X = poly_ring.gen()
        q = self._Fq.cardinality()
        return self._compute_goss_polynomial(n, q, poly_ring, X)
