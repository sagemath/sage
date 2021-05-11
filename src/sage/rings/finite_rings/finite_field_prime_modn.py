"""
Finite Prime Fields

AUTHORS:

- William Stein: initial version

- Martin Albrecht (2008-01): refactoring

TESTS::

    sage: k = GF(3)
    sage: TestSuite(k).run()
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#       Copyright (C) 2008 Martin Albrecht <malb@informatik.uni-bremen.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************`

from sage.rings.finite_rings.finite_field_base import FiniteField as FiniteField_generic
from sage.categories.finite_fields import FiniteFields
_FiniteFields = FiniteFields()

import sage.rings.finite_rings.integer_mod_ring as integer_mod_ring
from sage.rings.integer import Integer
import sage.rings.finite_rings.integer_mod as integer_mod

from sage.rings.integer_ring import ZZ
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing_generic

from sage.misc.persist import register_unpickle_override


class FiniteField_prime_modn(FiniteField_generic, integer_mod_ring.IntegerModRing_generic):
    r"""
    Finite field of order `p` where `p` is prime.

    EXAMPLES::

        sage: FiniteField(3)
        Finite Field of size 3

        sage: FiniteField(next_prime(1000))
        Finite Field of size 1009
    """
    def __init__(self, p, check=True, modulus=None):
        """
        Return a new finite field of order `p` where `p` is prime.

        INPUT:

        - ``p`` -- an integer at least 2

        - ``check`` -- bool (default: ``True``); if ``False``, do not
          check ``p`` for primality

        EXAMPLES::

            sage: F = FiniteField(3); F
            Finite Field of size 3
        """
        p = Integer(p)
        if check and not p.is_prime():
            raise ArithmeticError("p must be prime")
        self.__char = p
        # FiniteField_generic does nothing more than IntegerModRing_generic, and
        # it saves a non trivial overhead
        integer_mod_ring.IntegerModRing_generic.__init__(self, p, category=_FiniteFields)

        # If modulus is None, it will be created on demand as x-1
        # by the modulus() method.
        if modulus is not None:
            self._modulus = modulus

    def __reduce__(self):
        """
        For pickling.

        EXAMPLES::

            sage: k = FiniteField(5); type(k)
            <class 'sage.rings.finite_rings.finite_field_prime_modn.FiniteField_prime_modn_with_category'>
            sage: k is loads(dumps(k))
            True
        """
        return self._factory_data[0].reduce_data(self)

    def _coerce_map_from_(self, S):
        """
        This is called implicitly by arithmetic methods.

        EXAMPLES::

            sage: k = GF(7)
            sage: e = k(6)
            sage: e * 2 # indirect doctest
            5
            sage: 12 % 7
            5
            sage: ZZ.residue_field(7).hom(GF(7))(1)  # See trac 11319
            1
            sage: K.<w> = QuadraticField(337)  # See trac 11319
            sage: pp = K.ideal(13).factor()[0][0]
            sage: RF13 = K.residue_field(pp)
            sage: RF13.hom([GF(13)(1)])
            Ring morphism:
             From: Residue field of Fractional ideal (w + 18)
             To:   Finite Field of size 13
             Defn: 1 |--> 1

        Check that :trac:`19573` is resolved::

            sage: Integers(9).hom(GF(3))
            Natural morphism:
              From: Ring of integers modulo 9
              To:   Finite Field of size 3

            sage: Integers(9).hom(GF(5))
            Traceback (most recent call last):
            ...
            TypeError: natural coercion morphism from Ring of integers modulo 9 to Finite Field of size 5 not defined

        There is no coercion from a `p`-adic ring to its residue field::

            sage: GF(3).has_coerce_map_from(Zp(3))
            False
        """
        if S is int:
            return integer_mod.Int_to_IntegerMod(self)
        elif S is ZZ:
            return integer_mod.Integer_to_IntegerMod(self)
        elif isinstance(S, IntegerModRing_generic):
            from .residue_field import ResidueField_generic
            if (S.characteristic() % self.characteristic() == 0 and
                    (not isinstance(S, ResidueField_generic) or
                     S.degree() == 1)):
                try:
                    return integer_mod.IntegerMod_to_IntegerMod(S, self)
                except TypeError:
                    pass
        to_ZZ = ZZ._internal_coerce_map_from(S)
        if to_ZZ is not None:
            return integer_mod.Integer_to_IntegerMod(self) * to_ZZ

    def _convert_map_from_(self, R):
        """
        Conversion from p-adic fields.

        EXAMPLES::

            sage: GF(3).convert_map_from(Qp(3))
            Reduction morphism:
              From: 3-adic Field with capped relative precision 20
              To:   Finite Field of size 3
            sage: GF(3).convert_map_from(Zp(3))
            Reduction morphism:
              From: 3-adic Ring with capped relative precision 20
              To:   Finite Field of size 3
        """
        from sage.rings.padics.padic_generic import pAdicGeneric, ResidueReductionMap
        if isinstance(R, pAdicGeneric) and R.residue_field() is self:
            return ResidueReductionMap._create_(R, self)

    def construction(self):
        """
        Returns the construction of this finite field (for use by sage.categories.pushout)

        EXAMPLES::

            sage: GF(3).construction()
            (QuotientFunctor, Integer Ring)
        """
        return integer_mod_ring.IntegerModRing_generic.construction(self)

    def characteristic(self):
        r"""
        Return the characteristic of \code{self}.

        EXAMPLES::

            sage: k = GF(7)
            sage: k.characteristic()
            7
        """
        return self.__char

    def is_prime_field(self):
        """
        Return ``True`` since this is a prime field.

        EXAMPLES::

            sage: k.<a> = GF(3)
            sage: k.is_prime_field()
            True

            sage: k.<a> = GF(3^2)
            sage: k.is_prime_field()
            False
        """
        return True

    def polynomial(self, name=None):
        """
        Returns the polynomial ``name``.

        EXAMPLES::

            sage: k.<a> = GF(3)
            sage: k.polynomial()
            x
        """
        if name is None:
            name = self.variable_name()
        try:
            return self.__polynomial[name]
        except  AttributeError:
            f = self[name]([0, 1])
            try:
                self.__polynomial[name] = f
            except (KeyError, AttributeError):
                self.__polynomial = {}
                self.__polynomial[name] = f
            return f

    def order(self):
        """
        Return the order of this finite field.

        EXAMPLES::

            sage: k = GF(5)
            sage: k.order()
            5
        """
        return self.__char

    def gen(self, n=0):
        """
        Return a generator of ``self`` over its prime field, which is a
        root of ``self.modulus()``.

        Unless a custom modulus was given when constructing this prime
        field, this returns `1`.

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

            sage: k = GF(13)
            sage: k.gen()
            1
            sage: k = GF(1009, modulus="primitive")
            sage: k.gen()  # this gives a primitive element
            11
            sage: k.gen(1)
            Traceback (most recent call last):
            ...
            IndexError: only one generator
        """
        if n:
            raise IndexError("only one generator")
        try:
            return self.__gen
        except AttributeError:
            pass

        try:
            self.__gen = -(self._modulus[0])
        except AttributeError:
            self.__gen = self.one()
        return self.__gen

    def __iter__(self):
        """
        Return an iterator over ``self``.

        EXAMPLES::

            sage: list(GF(7))
            [0, 1, 2, 3, 4, 5, 6]

        We can even start iterating over something that would be too big
        to actually enumerate::

            sage: K = GF(next_prime(2^256))
            sage: all = iter(K)
            sage: next(all)
            0
            sage: next(all)
            1
            sage: next(all)
            2
        """
        yield self(0)
        i = one = self(1)
        while i:
            yield i
            i += one

    def degree(self):
        """
        Return the degree of ``self`` over its prime field.

        This always returns 1.

        EXAMPLES::

            sage: FiniteField(3).degree()
            1
        """
        return Integer(1)

register_unpickle_override(
    'sage.rings.finite_field_prime_modn', 'FiniteField_prime_modn',
    FiniteField_prime_modn)
