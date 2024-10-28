"""
Finite fields of characteristic 2
"""

#*****************************************************************************
#       Copyright (C) 2011 David Roe
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
from sage.libs.pari.all import pari
from sage.rings.integer import Integer
from sage.misc.superseded import deprecated_function_alias


def late_import():
    """
    Imports various modules after startup.

    EXAMPLES::

       sage: sage.rings.finite_rings.finite_field_ntl_gf2e.late_import()
       sage: sage.rings.finite_rings.finite_field_ntl_gf2e.GF2 is None # indirect doctest
       False
    """
    if "GF2" in globals():
        return
    global Cache_ntl_gf2e, GF, GF2

    import sage.rings.finite_rings.element_ntl_gf2e
    Cache_ntl_gf2e = sage.rings.finite_rings.element_ntl_gf2e.Cache_ntl_gf2e

    import sage.rings.finite_rings.finite_field_constructor
    GF = sage.rings.finite_rings.finite_field_constructor.GF
    GF2 = GF(2)


class FiniteField_ntl_gf2e(FiniteField):
    """
    Finite Field of characteristic 2 and order `2^n`.

    INPUT:

    - ``q`` -- `2^n` (must be 2 power)

    - ``names`` -- variable used for poly_repr (default: ``'a'``)

    - ``modulus`` -- a minimal polynomial to use for reduction

    - ``repr`` -- controls the way elements are printed to the user:
                 (default: ``'poly'``)

      - ``'poly'`` -- polynomial representation

    OUTPUT: finite field with characteristic 2 and cardinality `2^n`

    EXAMPLES::

        sage: k.<a> = GF(2^16)
        sage: type(k)
        <class 'sage.rings.finite_rings.finite_field_ntl_gf2e.FiniteField_ntl_gf2e_with_category'>
        sage: k.<a> = GF(2^1024)
        sage: k.modulus()
        x^1024 + x^19 + x^6 + x + 1
        sage: set_random_seed(6397)
        sage: k.<a> = GF(2^17, modulus='random')
        sage: k.modulus()
        x^17 + x^16 + x^15 + x^10 + x^8 + x^6 + x^4 + x^3 + x^2 + x + 1
        sage: k.modulus().is_irreducible()
        True
        sage: k.<a> = GF(2^211, modulus='minimal_weight')
        sage: k.modulus()
        x^211 + x^11 + x^10 + x^8 + 1
        sage: k.<a> = GF(2^211, modulus='conway')
        sage: k.modulus()
        x^211 + x^9 + x^6 + x^5 + x^3 + x + 1
        sage: k.<a> = GF(2^23, modulus='conway')
        sage: a.multiplicative_order() == k.order() - 1
        True
    """

    def __init__(self, q, names='a', modulus=None, repr='poly'):
        """
        Initialize ``self``.

        TESTS::

            sage: k.<a> = GF(2^100, modulus='strangeinput')
            Traceback (most recent call last):
            ...
            ValueError: no such algorithm for finding an irreducible polynomial: strangeinput
            sage: k.<a> = GF(2^20) ; type(k)
            <class 'sage.rings.finite_rings.finite_field_ntl_gf2e.FiniteField_ntl_gf2e_with_category'>
            sage: loads(dumps(k)) is k
            True
            sage: k1.<a> = GF(2^16)
            sage: k2.<a> = GF(2^17)
            sage: k1 == k2
            False
            sage: k3.<a> = GF(2^16, impl='pari_ffelt')
            sage: k1 == k3
            False

            sage: TestSuite(k).run()

            sage: k.<a> = GF(2^64)
            sage: k._repr_option('element_is_atomic')
            False
            sage: P.<x> = PolynomialRing(k)
            sage: (a+1)*x # indirect doctest
            (a + 1)*x
        """
        late_import()
        q = Integer(q)
        if q < 2:
            raise ValueError("q must be a 2-power")
        k = q.exact_log(2)
        if q != 1 << k:
            raise ValueError("q must be a 2-power")
        FiniteField.__init__(self, GF2, names, normalize=True)

        from sage.rings.polynomial.polynomial_element import Polynomial
        if not isinstance(modulus, Polynomial):
            raise TypeError("modulus must be a polynomial")

        self._cache = Cache_ntl_gf2e(self, k, modulus)
        self._modulus = modulus

    def characteristic(self):
        """
        Return the characteristic of ``self`` which is 2.

        EXAMPLES::

            sage: k.<a> = GF(2^16,modulus='random')
            sage: k.characteristic()
            2
        """
        return Integer(2)

    def order(self):
        """
        Return the cardinality of this field.

        EXAMPLES::

            sage: k.<a> = GF(2^64)
            sage: k.order()
            18446744073709551616
        """
        return self._cache.order()

    def degree(self):
        r"""
        If this field has cardinality `2^n` this method returns `n`.

        EXAMPLES::

            sage: k.<a> = GF(2^64)
            sage: k.degree()
            64
        """
        return self._cache.degree()

    def _element_constructor_(self, e):
        """
        Coerces several data types to ``self``.

        INPUT:

        - ``e`` -- data to coerce

        EXAMPLES::

            sage: k.<a> = GF(2^20)
            sage: k(1) # indirect doctest
            1
            sage: k(int(2))
            0

            sage: k('a+1')
            a + 1
            sage: k('b+1')
            Traceback (most recent call last):
            ...
            NameError: name 'b' is not defined

            sage: R.<x>=GF(2)[]
            sage: k(1+x+x^10+x^55)
            a^19 + a^17 + a^16 + a^15 + a^12 + a^11 + a^8 + a^6 + a^4 + a^2 + 1

            sage: V = k.vector_space(map=False)
            sage: v = V.random_element()
            sage: k(v) == sum(a^i if v[i] else 0 for i in range(len(v)))
            True
            sage: vector(k(v)) == v
            True

            sage: k(pari('Mod(1,2)*a^20'))
            a^10 + a^9 + a^7 + a^6 + a^5 + a^4 + a + 1
        """
        return self._cache.import_data(e)

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

            sage: k.<a> = GF(2^19)
            sage: k.gen() == a
            True
            sage: a
            a

        TESTS::

            sage: GF(2, impl='ntl').gen()
            1
            sage: GF(2, impl='ntl', modulus=polygen(GF(2)) ).gen()
            0
            sage: GF(2^19, 'a').gen(1)
            Traceback (most recent call last):
            ...
            IndexError: only one generator
        """
        if n:
            raise IndexError("only one generator")
        return self._cache._gen

    def prime_subfield(self):
        r"""
        Return the prime subfield `\GF{p}` of ``self`` if ``self`` is
        `\GF{p^n}`.

        EXAMPLES::

            sage: F.<a> = GF(2^16)
            sage: F.prime_subfield()
            Finite Field of size 2
        """
        return GF2

    def from_integer(self, number):
        r"""
        Given an integer `n` less than :meth:`cardinality` with base `2`
        representation `a_0 + 2 \cdot a_1 + \cdots + 2^k a_k`, returns
        `a_0 + a_1 \cdot x + \cdots + a_k x^k`, where `x` is the
        generator of this finite field.

        INPUT:

        - ``number`` -- integer

        EXAMPLES::

            sage: k.<a> = GF(2^48)
            sage: k.from_integer(2^43 + 2^15 + 1)
            a^43 + a^15 + 1
            sage: k.from_integer(33793)
            a^15 + a^10 + 1
            sage: 33793.digits(2) # little endian
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
        """
        return self._cache.fetch_int(number)

    fetch_int = deprecated_function_alias(33941, from_integer)

    def _pari_modulus(self):
        """
        Return PARI object which is equivalent to the
        polynomial/modulus of ``self``.

        EXAMPLES::

            sage: k1.<a> = GF(2^16)
            sage: k1._pari_modulus()
            Mod(1, 2)*a^16 + Mod(1, 2)*a^5 + Mod(1, 2)*a^3 + Mod(1, 2)*a^2 + Mod(1, 2)
        """
        f = pari(str(self.modulus()))
        return f.subst('x', 'a') * pari("Mod(1,%s)" % self.characteristic())
