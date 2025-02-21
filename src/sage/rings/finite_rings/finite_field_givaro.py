"""
Givaro finite fields

Finite fields that are implemented using Zech logs and the
cardinality must be less than `2^{16}`. By default, Conway polynomials are
used as minimal polynomial.
"""

#*****************************************************************************
#       Copyright (C) 2010-2012 David Roe
#       Copyright (C) 2012 Travis Scrimshaw
#       Copyright (C) 2013 Peter Bruin
#       Copyright (C) 2014 Jeroen Demeyer
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.finite_rings.finite_field_base import FiniteField
from sage.rings.integer import Integer
from sage.rings.finite_rings.element_givaro import Cache_givaro
from sage.libs.pari import pari


class FiniteField_givaro(FiniteField):
    """
    Finite field implemented using Zech logs and the cardinality must be
    less than `2^{16}`. By default, Conway polynomials are used as minimal
    polynomials.

    INPUT:

    - ``q`` -- `p^n` (must be prime power)

    - ``name`` -- (default: ``'a'``) variable used for
      :meth:`~sage.rings.finite_rings.element_givaro.FiniteField_givaroElement.poly_repr()`

    - ``modulus`` -- a minimal polynomial to use for reduction

    - ``repr`` -- (default: ``'poly'``) controls the way elements are printed
      to the user:

      - 'log': repr is
        :meth:`~sage.rings.finite_rings.element_givaro.FiniteField_givaroElement.log_repr()`
      - 'int': repr is
        :meth:`~sage.rings.finite_rings.element_givaro.FiniteField_givaroElement.int_repr()`
      - 'poly': repr is
        :meth:`~sage.rings.finite_rings.element_givaro.FiniteField_givaroElement.poly_repr()`

    - ``cache`` -- boolean (default: ``False``); if ``True`` a cache of all
      elements of this field is created. Thus, arithmetic does not create new
      elements which speeds calculations up. Also, if many elements are needed
      during a calculation this cache reduces the memory requirement as at most
      :meth:`order` elements are created.

    OUTPUT: Givaro finite field with characteristic `p` and cardinality `p^n`

    EXAMPLES:

    By default, Conway polynomials are used for extension fields::

        sage: k.<a> = GF(2**8)
        sage: -a ^ k.degree()
        a^4 + a^3 + a^2 + 1
        sage: f = k.modulus(); f
        x^8 + x^4 + x^3 + x^2 + 1

    You may enforce a modulus::

        sage: P.<x> = PolynomialRing(GF(2))
        sage: f = x^8 + x^4 + x^3 + x + 1 # Rijndael Polynomial
        sage: k.<a> = GF(2^8, modulus=f)
        sage: k.modulus()
        x^8 + x^4 + x^3 + x + 1
        sage: a^(2^8)
        a

    You may enforce a random modulus::

        sage: k = GF(3**5, 'a', modulus='random')
        sage: k.modulus() # random polynomial
        x^5 + 2*x^4 + 2*x^3 + x^2 + 2

    Three different representations are possible::

        sage: FiniteField(9, 'a', impl='givaro', repr='poly').gen()
        a
        sage: FiniteField(9, 'a', impl='givaro', repr='int').gen()
        3
        sage: FiniteField(9, 'a', impl='givaro', repr='log').gen()
        1

    For prime fields, the default modulus is the polynomial `x - 1`,
    but you can ask for a different modulus::

        sage: GF(1009, impl='givaro').modulus()
        x + 1008
        sage: GF(1009, impl='givaro', modulus='conway').modulus()
        x + 998
    """
    def __init__(self, q, name='a', modulus=None, repr='poly', cache=False):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: k.<a> = GF(2^3)
            sage: j.<b> = GF(3^4)
            sage: k == j
            False

            sage: GF(2^3,'a') == copy(GF(2^3,'a'))
            True
            sage: TestSuite(GF(2^3, 'a')).run()
        """
        if repr not in ['int', 'log', 'poly']:
            raise ValueError("Unknown representation %s" % repr)

        q = Integer(q)
        if q < 2:
            raise ValueError("q must be a prime power")
        F = q.factor()
        if len(F) > 1:
            raise ValueError("q must be a prime power")
        p = F[0][0]
        k = F[0][1]

        if q >= 1 << 16:
            raise ValueError("q must be < 2^16")

        from .finite_field_constructor import GF
        FiniteField.__init__(self, GF(p), name, normalize=False)

        from sage.rings.polynomial.polynomial_element import Polynomial
        if not isinstance(modulus, Polynomial):
            raise TypeError("modulus must be a polynomial")

        self._cache = Cache_givaro(self, p, k, modulus, repr, cache)
        self._modulus = modulus

    def characteristic(self):
        """
        Return the characteristic of this field.

        EXAMPLES::

            sage: p = GF(19^5,'a').characteristic(); p
            19
            sage: type(p)
            <class 'sage.rings.integer.Integer'>
        """
        return Integer(self._cache.characteristic())

    def order(self):
        """
        Return the cardinality of this field.

        OUTPUT: integer; the number of elements in ``self``

        EXAMPLES::

            sage: n = GF(19^5,'a').order(); n
            2476099
            sage: type(n)
            <class 'sage.rings.integer.Integer'>
        """
        return self._cache.order()

    def degree(self):
        r"""
        If the cardinality of ``self`` is `p^n`, then this returns `n`.

        OUTPUT: integer; the degree

        EXAMPLES::

            sage: GF(3^4,'a').degree()
            4
        """
        return Integer(self._cache.exponent())

    def _repr_option(self, key):
        """
        Metadata about the :meth:`_repr_` output.

        See :meth:`sage.structure.parent._repr_option` for details.

        EXAMPLES::

            sage: GF(23**3, 'a', repr='log')._repr_option('element_is_atomic')
            True
            sage: GF(23**3, 'a', repr='int')._repr_option('element_is_atomic')
            True
            sage: GF(23**3, 'a', repr='poly')._repr_option('element_is_atomic')
            False
        """
        if key == 'element_is_atomic':
            return self._cache.repr != 0   # 0 means repr='poly'
        return super()._repr_option(key)

    def random_element(self, *args, **kwds):
        """
        Return a random element of ``self``.

        EXAMPLES::

            sage: k = GF(23**3, 'a')
            sage: e = k.random_element()
            sage: e.parent() is k
            True
            sage: type(e)
            <class 'sage.rings.finite_rings.element_givaro.FiniteField_givaroElement'>

            sage: P.<x> = PowerSeriesRing(GF(3^3, 'a'))
            sage: P.random_element(5).parent() is P
            True
        """
        return self._cache.random_element()

    def _element_constructor_(self, e):
        """
        Coerces several data types to ``self``.

        INPUT:

        - ``e`` -- data to coerce

        EXAMPLES:

        :class:`FiniteField_givaroElement` are accepted where the parent
        is either ``self``, equals ``self`` or is the prime subfield::

            sage: k = GF(2**8, 'a')
            sage: k.gen() == k(k.gen())
            True

        Floats, ints, longs, Integer are interpreted modulo characteristic::

            sage: k(2) # indirect doctest
            0

        Floats are converted like integers::

            sage: k(float(2.0))
            0

        Rational are interpreted as ``self(numerator)/self(denominator)``.
        Both may not be greater than :meth:`characteristic`.
        ::

            sage: k = GF(3**8, 'a')
            sage: k(1/2) == k(1)/k(2)
            True

        Free module elements over :meth:`prime_subfield()` are interpreted
        'little endian'::

            sage: k = GF(2**8, 'a')
            sage: e = k.vector_space(map=False).gen(1); e
            (0, 1, 0, 0, 0, 0, 0, 0)
            sage: k(e)
            a

        ``None`` yields zero::

            sage: k(None)
            0

        Strings are evaluated as polynomial representation of elements in
        ``self``::

            sage: k('a^2+1')
            a^2 + 1

        Univariate polynomials coerce into finite fields by evaluating
        the polynomial at the field's generator::

            sage: R.<x> = QQ[]
            sage: k.<a> = FiniteField(5^2, 'a', impl='givaro')
            sage: k(R(2/3))
            4
            sage: k(x^2)
            a + 3
            sage: R.<x> = GF(5)[]
            sage: k(x^3-2*x+1)
            2*a + 4

            sage: x = polygen(QQ)
            sage: k(x^25)
            a

            sage: Q.<q> = FiniteField(5^3, 'q', impl='givaro')
            sage: L = GF(5)
            sage: LL.<xx> = L[]
            sage: Q(xx^2 + 2*xx + 4)
            q^2 + 2*q + 4

        Multivariate polynomials only coerce if constant::

            sage: R = k['x,y,z']; R
            Multivariate Polynomial Ring in x, y, z over Finite Field in a of size 5^2
            sage: k(R(2))
            2
            sage: R = QQ['x,y,z']
            sage: k(R(1/5))
            Traceback (most recent call last):
            ...
            ZeroDivisionError: division by zero in finite field

        PARI elements are interpreted as finite field elements; this PARI
        flexibility is (absurdly!) liberal::

            sage: k.<a> = GF(2^8)
            sage: k(pari('Mod(1,2)'))
            1
            sage: k(pari('Mod(2,3)'))
            a
            sage: k(pari('Mod(1,3)*a^20'))
            a^7 + a^5 + a^4 + a^2
            sage: k(pari('O(x)'))
            Traceback (most recent call last):
            ...
            TypeError: unable to convert PARI t_SER to Finite Field in a of size 2^8

        We can coerce from PARI finite field implementations::

            sage: K.<a> = GF(3^10, impl='givaro')
            sage: a^20
            2*a^9 + 2*a^8 + a^7 + 2*a^5 + 2*a^4 + 2*a^3 + 1
            sage: M.<c> = GF(3^10, impl='pari_ffelt')
            sage: K(c^20)
            2*a^9 + 2*a^8 + a^7 + 2*a^5 + 2*a^4 + 2*a^3 + 1

        GAP elements need to be finite field elements::

            sage: x = gap('Z(13)')
            sage: F = FiniteField(13, impl='givaro')
            sage: F(x)
            2
            sage: F(gap('0*Z(13)'))
            0
            sage: F = FiniteField(13^2, 'a', impl='givaro')
            sage: x = gap('Z(13)')
            sage: F(x)
            2
            sage: x = gap('Z(13^2)^3')
            sage: F(x)
            12*a + 11
            sage: F.multiplicative_generator()^3
            12*a + 11

            sage: k.<a> = GF(29^3)
            sage: k(48771/1225)
            28

            sage: F9 = FiniteField(9, impl='givaro', prefix='a')
            sage: F81 = FiniteField(81, impl='givaro', prefix='a')
            sage: F81(F9.gen())
            2*a4^3 + 2*a4^2 + 1
        """
        return self._cache.element_from_data(e)

    def gen(self, n=0):
        r"""
        Return a generator of ``self`` over its prime field, which is a
        root of ``self.modulus()``.

        INPUT:

        - ``n`` -- must be 0

        OUTPUT:

        An element `a` of ``self`` such that ``self.modulus()(a) == 0``.

        .. WARNING::

            This generator is not guaranteed to be a generator for the
            multiplicative group.  To obtain the latter, use
            :meth:`~sage.rings.finite_rings.finite_field_base.FiniteFields.multiplicative_generator()`
            or use the ``modulus="primitive"`` option when constructing
            the field.

        EXAMPLES::

            sage: k = GF(3^4, 'b'); k.gen()
            b
            sage: k.gen(1)
            Traceback (most recent call last):
            ...
            IndexError: only one generator
            sage: F = FiniteField(31, impl='givaro')
            sage: F.gen()
            1
        """
        if n:
            raise IndexError("only one generator")
        return self._cache.gen()

    def prime_subfield(self):
        r"""
        Return the prime subfield `\GF{p}` of ``self`` if ``self`` is `\GF{p^n}`.

        EXAMPLES::

            sage: GF(3^4, 'b').prime_subfield()
            Finite Field of size 3

            sage: S.<b> = GF(5^2); S
            Finite Field in b of size 5^2
            sage: S.prime_subfield()
            Finite Field of size 5
            sage: type(S.prime_subfield())
            <class 'sage.rings.finite_rings.finite_field_prime_modn.FiniteField_prime_modn_with_category'>
        """
        try:
            return self._prime_subfield
        except AttributeError:
            from .finite_field_constructor import GF
            self._prime_subfield = GF(self.characteristic())
            return self._prime_subfield

    def log_to_int(self, n):
        r"""
        Given an integer `n` this method returns ``i`` where ``i``
        satisfies `g^n = i` where `g` is the generator of ``self``; the
        result is interpreted as an integer.

        INPUT:

        - ``n`` -- log representation of a finite field element

        OUTPUT: integer representation of a finite field element

        EXAMPLES::

            sage: k = GF(2**8, 'a')
            sage: k.log_to_int(4)
            16
            sage: k.log_to_int(20)
            180
        """
        return self._cache.log_to_int(n)

    def int_to_log(self, n):
        r"""
        Given an integer `n` this method returns `i` where `i` satisfies
        `g^i = n \mod p` where `g` is the generator and `p` is the
        characteristic of ``self``.

        INPUT:

        - ``n`` -- integer representation of a finite field element

        OUTPUT: log representation of ``n``

        EXAMPLES::

            sage: k = GF(7**3, 'a')
            sage: k.int_to_log(4)
            228
            sage: k.int_to_log(3)
            57
            sage: k.gen()^57
            3
        """
        return self._cache.int_to_log(n)

    def from_integer(self, n):
        r"""
        Given an integer `n` return a finite field element in ``self``
        which equals `n` under the condition that :meth:`gen()` is set to
        :meth:`characteristic()`.

        EXAMPLES::

            sage: k.<a> = GF(2^8)
            sage: k.from_integer(8)
            a^3
            sage: e = k.from_integer(151); e
            a^7 + a^4 + a^2 + a + 1
            sage: 2^7 + 2^4 + 2^2 + 2 + 1
            151
        """
        return self._cache.fetch_int(n)

    def _pari_modulus(self):
        """
        Return the modulus of ``self`` in a format for PARI.

        EXAMPLES::

            sage: GF(3^4,'a')._pari_modulus()
            Mod(1, 3)*a^4 + Mod(2, 3)*a^3 + Mod(2, 3)
        """
        f = pari(str(self.modulus()))
        return f.subst('x', 'a') * pari("Mod(1,%s)" % self.characteristic())

    def __iter__(self):
        """
        Finite fields may be iterated over.

        EXAMPLES::

            sage: list(GF(2**2, 'a'))
            [0, a, a + 1, 1]
        """
        from .element_givaro import FiniteField_givaro_iterator
        return FiniteField_givaro_iterator(self._cache)

    def a_times_b_plus_c(self, a, b, c):
        """
        Return ``a*b + c``. This is faster than multiplying ``a`` and ``b``
        first and adding ``c`` to the result.

        INPUT:

        - ``a``, ``b``, ``c`` -- :class:`~~sage.rings.finite_rings.element_givaro.FiniteField_givaroElement`

        EXAMPLES::

            sage: k.<a> = GF(2**8)
            sage: k.a_times_b_plus_c(a,a,k(1))
            a^2 + 1
        """
        return self._cache.a_times_b_plus_c(a, b, c)

    def a_times_b_minus_c(self, a, b, c):
        """
        Return ``a*b - c``.

        INPUT:

        - ``a``, ``b``, ``c`` -- :class:`~sage.rings.finite_rings.element_givaro.FiniteField_givaroElement`

        EXAMPLES::

            sage: k.<a> = GF(3**3)
            sage: k.a_times_b_minus_c(a,a,k(1))
            a^2 + 2
        """
        return self._cache.a_times_b_minus_c(a, b, c)

    def c_minus_a_times_b(self, a, b, c):
        """
        Return ``c - a*b``.

        INPUT:

        - ``a``, ``b``, ``c`` -- :class:`~sage.rings.finite_rings.element_givaro.FiniteField_givaroElement`

        EXAMPLES::

            sage: k.<a> = GF(3**3)
            sage: k.c_minus_a_times_b(a,a,k(1))
            2*a^2 + 1
        """
        return self._cache.c_minus_a_times_b(a, b, c)

    def frobenius_endomorphism(self, n=1):
        """
        INPUT:

        - ``n`` -- integer (default: 1)

        OUTPUT:

        The `n`-th power of the absolute arithmetic Frobenius
        endomorphism on this finite field.

        EXAMPLES::

            sage: k.<t> = GF(3^5)
            sage: Frob = k.frobenius_endomorphism(); Frob
            Frobenius endomorphism t |--> t^3 on Finite Field in t of size 3^5

            sage: a = k.random_element()
            sage: Frob(a) == a^3
            True

        We can specify a power::

            sage: k.frobenius_endomorphism(2)
            Frobenius endomorphism t |--> t^(3^2) on Finite Field in t of size 3^5

        The result is simplified if possible::

            sage: k.frobenius_endomorphism(6)
            Frobenius endomorphism t |--> t^3 on Finite Field in t of size 3^5
            sage: k.frobenius_endomorphism(5)
            Identity endomorphism of Finite Field in t of size 3^5

        Comparisons work::

            sage: k.frobenius_endomorphism(6) == Frob
            True
            sage: from sage.categories.morphism import IdentityMorphism
            sage: k.frobenius_endomorphism(5) == IdentityMorphism(k)
            True

        AUTHOR:

        - Xavier Caruso (2012-06-29)
        """
        from sage.rings.finite_rings.hom_finite_field_givaro import FrobeniusEndomorphism_givaro
        return FrobeniusEndomorphism_givaro(self, n)
