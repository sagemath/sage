"""
`p`-Adic Generic Nodes

This file contains a bunch of intermediate classes for the `p`-adic
parents, allowing a function to be implemented at the right level of
generality.

AUTHORS:

- David Roe
"""

#*****************************************************************************
#       Copyright (C) 2007-2013 David Roe <roed.math@gmail.com>
#                               William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.padics.local_generic import LocalGeneric
from sage.rings.padics.padic_generic import pAdicGeneric
from sage.rings.ring import EuclideanDomain, Field
from sage.rings.padics.padic_base_generic import pAdicBaseGeneric
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

class CappedAbsoluteGeneric(LocalGeneric):
    def is_capped_absolute(self):
        """
        Returns whether this `p`-adic ring bounds precision in a
        capped absolute fashion.

        The absolute precision of an element is the power of `p` modulo
        which that element is defined.  In a capped absolute ring, the
        absolute precision of elements are bounded by a constant
        depending on the ring.

        EXAMPLES::

            sage: R = ZpCA(5, 15)
            sage: R.is_capped_absolute()
            True
            sage: R(5^7)
            5^7 + O(5^15)
            sage: S = Zp(5, 15)
            sage: S.is_capped_absolute()
            False
            sage: S(5^7)
            5^7 + O(5^22)
        """
        return True

    def _prec_type(self):
        """
        Returns the precision handling type.

        EXAMPLES::

            sage: ZpCA(5)._prec_type()
            'capped-abs'
        """
        return 'capped-abs'

class CappedRelativeGeneric(LocalGeneric):
    def is_capped_relative(self):
        """
        Returns whether this `p`-adic ring bounds precision in a capped
        relative fashion.

        The relative precision of an element is the power of p modulo
        which the unit part of that element is defined.  In a capped
        relative ring, the relative precision of elements are bounded
        by a constant depending on the ring.

        EXAMPLES::

            sage: R = ZpCA(5, 15)
            sage: R.is_capped_relative()
            False
            sage: R(5^7)
            5^7 + O(5^15)
            sage: S = Zp(5, 15)
            sage: S.is_capped_relative()
            True
            sage: S(5^7)
            5^7 + O(5^22)
        """
        return True

    def _prec_type(self):
        """
        Returns the precision handling type.

        EXAMPLES::

            sage: Zp(5)._prec_type()
            'capped-rel'
        """
        return 'capped-rel'

class FixedModGeneric(LocalGeneric):
    def is_fixed_mod(self):
        """
        Returns whether this `p`-adic ring bounds precision in a fixed
        modulus fashion.

        The absolute precision of an element is the power of p modulo
        which that element is defined.  In a fixed modulus ring, the
        absolute precision of every element is defined to be the
        precision cap of the parent.  This means that some operations,
        such as division by `p`, don't return a well defined answer.

        EXAMPLES::

            sage: R = ZpFM(5,15)
            sage: R.is_fixed_mod()
            True
            sage: R(5^7,absprec=9)
            5^7 + O(5^15)
            sage: S = ZpCA(5, 15)
            sage: S.is_fixed_mod()
            False
            sage: S(5^7,absprec=9)
            5^7 + O(5^9)
        """
        return True

    def _prec_type(self):
        """
        Returns the precision handling type.

        EXAMPLES::

            sage: ZpFM(5)._prec_type()
            'fixed-mod'
        """
        return 'fixed-mod'

class CappedRelativeRingGeneric(CappedRelativeGeneric):
    pass
class CappedRelativeFieldGeneric(CappedRelativeGeneric):#, sage.rings.ring.Field):
    pass


def is_pAdicRing(R):
    """
    Returns ``True`` if and only if ``R`` is a `p`-adic ring (not a
    field).

    EXAMPLES::

        sage: is_pAdicRing(Zp(5))
        True
        sage: is_pAdicRing(RR)
        False
    """
    return isinstance(R, pAdicRingGeneric)

class pAdicRingGeneric(pAdicGeneric, EuclideanDomain):
    def is_field(self, proof = True):
        """
        Returns whether this ring is actually a field, ie ``False``.

        EXAMPLES::

            sage: Zp(5).is_field()
            False
        """
        return False


    def krull_dimension(self):
        r"""
        Returns the Krull dimension of self, i.e. 1

        INPUT:

        - self -- a `p`-adic ring

        OUTPUT:

        - the Krull dimension of self.  Since self is a `p`-adic ring,
          this is 1.

        EXAMPLES::

            sage: Zp(5).krull_dimension()
            1
        """
        return 1

def is_pAdicField(R):
    """
    Returns ``True`` if and only if ``R`` is a `p`-adic field.

    EXAMPLES::

        sage: is_pAdicField(Zp(17))
        False
        sage: is_pAdicField(Qp(17))
        True
    """
    return isinstance(R, pAdicFieldGeneric)

class pAdicFieldGeneric(pAdicGeneric, Field):
    pass

    #def class_field(self, group=None, map=None, generators=None):
    #    raise NotImplementedError

    #def composite(self, subfield1, subfield2):
    #    raise NotImplementedError

    #def norm_equation(self):
    #    raise NotImplementedError

    #def norm_group(self):
    #    raise NotImplementedError

    #def norm_group_discriminant(self, group=None, map=None, generators=None):
    #    raise NotImplementedError

    #def number_of_extensions(self, degree, discriminant=None, e=None, f=None):
    #    raise NotImplementedError

    #def list_of_extensions(self, degree, discriminant=None, e=None, f=None):
    #    raise NotImplementedError

    #def subfield(self, list):
    #    raise NotImplementedError

    #def subfield_lattice(self):
    #    raise NotImplementedError

    #def subfields_of_degree(self, n):
    #    raise NotImplementedError

class pAdicFixedModRingGeneric(pAdicRingGeneric, FixedModGeneric):
    pass
class pAdicCappedAbsoluteRingGeneric(pAdicRingGeneric, CappedAbsoluteGeneric):
    pass
class pAdicCappedRelativeRingGeneric(pAdicRingGeneric, CappedRelativeRingGeneric):
    pass
class pAdicCappedRelativeFieldGeneric(pAdicFieldGeneric, CappedRelativeFieldGeneric):
    pass

class pAdicRingBaseGeneric(pAdicBaseGeneric, pAdicRingGeneric):
    def construction(self):
        """
        Returns the functorial construction of self, namely,
        completion of the rational numbers with respect a given prime.

        Also preserves other information that makes this field unique
        (e.g. precision, rounding, print mode).

        EXAMPLE::

            sage: K = Zp(17, 8, print_mode='val-unit', print_sep='&')
            sage: c, L = K.construction(); L
            Integer Ring
            sage: c(L)
            17-adic Ring with capped relative precision 8
            sage: K == c(L)
            True
        """
        from sage.categories.pushout import CompletionFunctor
        return (CompletionFunctor(self.prime(),
                                  self.precision_cap(),
                                  {'print_mode':self._printer.dict(), 'type':self._prec_type(), 'names':self._names}),
                ZZ)

    def random_element(self, algorithm='default'):
        r"""
        Returns a random element of self, optionally using the
        algorithm argument to decide how it generates the
        element. Algorithms currently implemented:

        - default: Choose `a_i`, `i >= 0`, randomly between `0` and
          `p-1` until a nonzero choice is made. Then continue choosing
          `a_i` randomly between `0` and `p-1` until we reach
          precision_cap, and return `\sum a_i p^i`.

        EXAMPLES::

            sage: Zp(5,6).random_element()
            3 + 3*5 + 2*5^2 + 3*5^3 + 2*5^4 + 5^5 + O(5^6)
            sage: ZpCA(5,6).random_element()
            4*5^2 + 5^3 + O(5^6)
            sage: ZpFM(5,6).random_element()
            2 + 4*5^2 + 2*5^4 + 5^5 + O(5^6)
        """
        if (algorithm == 'default'):
            if self.is_capped_relative():
                i = 0
                a_i = ZZ.random_element(self.prime())
                while a_i.is_zero():
                    i += 1
                    a_i = ZZ.random_element(self.prime())
                return self((self.prime()**i)*(a_i + self.prime()*ZZ.random_element(self.prime_pow.pow_Integer_Integer(self.precision_cap()-1))))
            else:
                return self(ZZ.random_element(self.prime_pow.pow_Integer_Integer(self.precision_cap())))
        else:
            raise NotImplementedError("Don't know %s algorithm"%algorithm)

    #def unit_group(self):
    #    raise NotImplementedError

    #def unit_group_gens(self):
    #    raise NotImplementedError

    #def principal_unit_group(self):
    #    raise NotImplementedError

class pAdicFieldBaseGeneric(pAdicBaseGeneric, pAdicFieldGeneric):
    def composite(self, subfield1, subfield2):
        r"""
        Returns the composite of two subfields of self, i.e., the
        largest subfield containing both

        INPUT:

        - ``self`` -- a `p`-adic field
        - ``subfield1`` -- a subfield
        - ``subfield2`` -- a subfield

        OUTPUT:

        - the composite of subfield1 and subfield2

        EXAMPLES::

            sage: K = Qp(17); K.composite(K, K) is K
            True
        """
        #should be overridden for extension fields
        if (subfield1 is self) and (subfield2 is self):
            return self
        raise ValueError("Arguments must be subfields of self.")

    def subfields_of_degree(self, n):
        r"""
        Returns the number of subfields of self of degree `n`

        INPUT:

        - ``self`` -- a `p`-adic field
        - ``n`` -- an integer

        OUTPUT:

        - integer -- the number of subfields of degree ``n`` over self.base_ring()

        EXAMPLES::

            sage: K = Qp(17)
            sage: K.subfields_of_degree(1)
            1
        """
        if n == 1:
            return 1
        else:
            return 0

    def subfield(self, list):
        r"""
        Returns the subfield generated by the elements in list

        INPUT:

        - ``self`` -- a `p`-adic field
        - ``list`` -- a list of elements of ``self``

        OUTPUT:

        - the subfield of ``self`` generated by the elements of list

        EXAMPLES::

            sage: K = Qp(17); K.subfield([K(17), K(1827)]) is K
            True
        """
        for x in list:
            if not self.__contains__(x):
                raise TypeError("Members of the list of generators must be elements of self.")
        return self

    def construction(self):
        """
        Returns the functorial construction of ``self``, namely,
        completion of the rational numbers with respect a given prime.

        Also preserves other information that makes this field unique
        (e.g. precision, rounding, print mode).

        EXAMPLE::

            sage: K = Qp(17, 8, print_mode='val-unit', print_sep='&')
            sage: c, L = K.construction(); L
            Rational Field
            sage: c(L)
            17-adic Field with capped relative precision 8
            sage: K == c(L)
            True
        """
        from sage.categories.pushout import CompletionFunctor
        return (CompletionFunctor(self.prime(),
                                  self.precision_cap(),
                                  {'print_mode':self._printer.dict(), 'type':self._prec_type(), 'names':self._names}),
                QQ)
