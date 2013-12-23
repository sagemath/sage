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
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.rings.finite_rings.finite_field_base import FiniteField as FiniteField_generic
from sage.categories.finite_fields import FiniteFields
_FiniteFields = FiniteFields()

import sage.rings.finite_rings.integer_mod_ring as integer_mod_ring
import sage.rings.integer as integer
import sage.rings.finite_rings.integer_mod as integer_mod
import sage.rings.arith as arith

from sage.rings.integer_ring import ZZ
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing_generic
import sage.structure.factorization as factorization
from sage.structure.parent import Parent

class FiniteField_prime_modn(FiniteField_generic, integer_mod_ring.IntegerModRing_generic):
    r"""
    Finite field of order `p` where `p` is prime.

    EXAMPLES::

        sage: FiniteField(3)
        Finite Field of size 3

        sage: FiniteField(next_prime(1000))
        Finite Field of size 1009
    """
    def __init__(self, p, name=None, check=True):
        """
        Return a new finite field of order `p` where `p` is prime.

        INPUT:

        - ``p`` -- an integer at least 2

        - ``name`` -- ignored

        - ``check`` -- bool (default: ``True``); if ``False``, do not
          check ``p`` for primality

        EXAMPLES::

            sage: F = FiniteField(3); F
            Finite Field of size 3
        """
        p = integer.Integer(p)
        if check and not arith.is_prime(p):
            raise ArithmeticError, "p must be prime"
        self.__char = p
        self._kwargs = {}
        # FiniteField_generic does nothing more than IntegerModRing_generic, and
        # it saves a non trivial overhead
        integer_mod_ring.IntegerModRing_generic.__init__(self, p, category = _FiniteFields)

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

    def __cmp__(self, other):
        r"""
        Compare ``self`` with ``other``.

        Two finite prime fields are considered equal if their characteristic
        is equal.

        EXAMPLES::

            sage: K = FiniteField(3)
            sage: copy(K) == K
            True
        """
        if not isinstance(other, FiniteField_prime_modn):
            return cmp(type(self), type(other))
#        elif other.__class__ != FiniteField_prime_modn:
#            return -cmp(other, self)
        return cmp(self.__char, other.__char)

    def __richcmp__(left, right, op):
        r"""
        Compare ``self`` with ``right``.

        EXAMPLES::

            sage: k = GF(2)
            sage: j = GF(3)
            sage: k == j
            False

            sage: GF(2) == copy(GF(2))
            True
        """
        return left._richcmp_helper(right, op)

    def _is_valid_homomorphism_(self, codomain, im_gens):
        """
        This is called implicitly by the ``hom`` constructor.

        EXAMPLES::

            sage: k = GF(73^2,'a')
            sage: f = k.modulus()
            sage: r = f.change_ring(k).roots()
            sage: k.hom([r[0][0]]) # indirect doctest
            Ring endomorphism of Finite Field in a of size 73^2
              Defn: a |--> 72*a + 3
        """
        try:
            return im_gens[0] == codomain._coerce_(self.gen(0))
        except TypeError:
            return False

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
             From: Residue field of Fractional ideal (w - 18)
             To:   Finite Field of size 13
             Defn: 1 |--> 1
        """
        if S is int:
            return integer_mod.Int_to_IntegerMod(self)
        elif S is ZZ:
            return integer_mod.Integer_to_IntegerMod(self)
        elif isinstance(S, IntegerModRing_generic):
            from sage.rings.residue_field import ResidueField_generic
            if S.characteristic() == self.characteristic() and \
               (not isinstance(S, ResidueField_generic) or S.degree() == 1):
                try:
                    return integer_mod.IntegerMod_to_IntegerMod(S, self)
                except TypeError:
                    pass
        to_ZZ = ZZ.coerce_map_from(S)
        if to_ZZ is not None:
            return integer_mod.Integer_to_IntegerMod(self) * to_ZZ

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

    def modulus(self):
        """
        Return the minimal polynomial of ``self``, which is always `x - 1`.

        EXAMPLES::

            sage: k = GF(199)
            sage: k.modulus()
            x + 198
        """
        try:
            return self.__modulus
        except AttributeError:
            x = self['x'].gen()
            self.__modulus = x - 1
        return self.__modulus

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
            from sage.rings.finite_rings.constructor import FiniteField
            R = FiniteField(self.characteristic())[name]
            f = self[name]([0,1])
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
        Return generator of this finite field as an extension of its
        prime field.

        .. NOTE::

            If you want a primitive element for this finite field
            instead, use :meth:`multiplicative_generator()`.

        EXAMPLES::

            sage: k = GF(13)
            sage: k.gen()
            1
            sage: k.gen(1)
            Traceback (most recent call last):
            ...
            IndexError: only one generator
        """
        if n != 0:
            raise IndexError, "only one generator"
        return self(1)

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
            sage: all.next()
            0
            sage: all.next()
            1
            sage: all.next()
            2
        """
        yield self(0)
        i = one = self(1)
        while i:
            yield i
            i += one

    def degree(self):
        """
        Returns the degree of the finite field, which is a positive
        integer.

        EXAMPLES::

            sage: FiniteField(3).degree()
            1
        """
        return integer.Integer(1)
