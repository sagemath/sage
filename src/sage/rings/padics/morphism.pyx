# sage.doctest: needs sage.libs.ntl
"""
Frobenius endomorphisms on `p`-adic fields
"""
# ****************************************************************************
#       Copyright (C) 2013 Xavier Caruso <xavier.caruso@normalesup.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.integer cimport Integer
from sage.rings.infinity import Infinity

from sage.categories.homset import Hom
from sage.structure.element cimport Element
from sage.structure.richcmp cimport (richcmp, rich_to_bool,
                                     richcmp_not_equal)

from sage.rings.morphism cimport RingHomomorphism
from sage.rings.padics.padic_generic import pAdicGeneric


cdef class FrobeniusEndomorphism_padics(RingHomomorphism):
    """
    A class implementing Frobenius endomorphisms on `p`-adic fields.
    """
    def __init__ (self,domain,n=1):
        """
        INPUT:

        - ``domain`` -- an unramified `p`-adic field

        - ``n`` -- integer (default: 1)

        .. NOTE::

            `n` may be negative.

        OUTPUT:

        The `n`-th power of the absolute (arithmetic) Frobenius
        endomorphism on ``domain``

        TESTS::

            sage: from sage.rings.padics.morphism import FrobeniusEndomorphism_padics
            sage: K.<a> = Qq(5^3)
            sage: FrobeniusEndomorphism_padics(K)
            Frobenius endomorphism on 5-adic Unramified Extension ... lifting a |--> a^5 on the residue field
            sage: FrobeniusEndomorphism_padics(K,2)
            Frobenius endomorphism on 5-adic Unramified Extension ... lifting a |--> a^(5^2) on the residue field

            sage: FrobeniusEndomorphism_padics(K,a)
            Traceback (most recent call last):
            ...
            TypeError: n (=a + O(5^20)) is not an integer

            sage: K = Qp(5)
            sage: x = polygen(ZZ, 'x')
            sage: L.<pi> = K.extension(x^2 - 5)
            sage: FrobeniusEndomorphism_padics(L)
            Traceback (most recent call last):
            ...
            TypeError: The domain must be unramified

            sage: FrobeniusEndomorphism_padics(QQ)
            Traceback (most recent call last):
            ...
            TypeError: The domain must be an instance of pAdicGeneric
        """
        if not isinstance(domain, pAdicGeneric):
            raise TypeError("The domain must be an instance of pAdicGeneric")
        if domain.absolute_e() != 1:
            raise TypeError("The domain must be unramified")
        try:
            n = Integer(n)
        except (ValueError, TypeError):
            raise TypeError("n (=%s) is not an integer" % n)

        self._degree = domain.absolute_f()
        self._power = n % self._degree
        self._order = self._degree / domain.absolute_f().gcd(self._power)
        RingHomomorphism.__init__(self, Hom(domain, domain))

    cdef dict _extra_slots(self):
        """
        Helper for copying and pickling.

        TESTS::

            sage: R.<x> = QQ[]
            sage: U.<a> = Qp(2).extension(x^2 + x + 1)
            sage: F = U.frobenius_endomorphism(); F
            Frobenius endomorphism on 2-adic Unramified Extension Field in a defined by x^2 + x + 1 lifting a |--> a^2 on the residue field
            sage: copy(F)
            Frobenius endomorphism on 2-adic Unramified Extension Field in a defined by x^2 + x + 1 lifting a |--> a^2 on the residue field
        """
        slots = RingHomomorphism._extra_slots(self)
        slots['_degree'] = self._degree
        slots['_power'] = self._power
        slots['_order'] = self._order
        return slots

    cdef _update_slots(self, dict slots):
        """
        Helper for copying and pickling.

        TESTS::

            sage: R.<x> = ZZ[]
            sage: U.<a> = Zp(2).extension(x^2 + x + 1)
            sage: F = U.frobenius_endomorphism(); F
            Frobenius endomorphism on 2-adic Unramified Extension Ring in a defined by x^2 + x + 1 lifting a |--> a^2 on the residue field
            sage: loads(dumps(F)) == F
            True
        """
        self._degree = slots['_degree']
        self._power = slots['_power']
        self._order = slots['_order']
        RingHomomorphism._update_slots(self, slots)

    def _repr_(self):
        """
        Return a string representation of this endomorphism.

        EXAMPLES::

            sage: K.<a> = Qq(5^3)
            sage: Frob = K.frobenius_endomorphism(); Frob
            Frobenius endomorphism on 5-adic Unramified Extension ... lifting a |--> a^5 on the residue field

            sage: Frob._repr_()
            'Frobenius endomorphism on 5-adic Unramified Extension ... lifting a |--> a^5 on the residue field'
        """
        name = self.domain().variable_name()
        if self._power == 0:
            s = "Identity endomorphism of %s" % self.domain()
        elif self._power == 1:
            s = "Frobenius endomorphism on %s lifting %s |--> %s^%s on the residue field" % (self.domain(), name, name, self.domain().residue_characteristic())
        else:
            s = "Frobenius endomorphism on %s lifting %s |--> %s^(%s^%s) on the residue field" % (self.domain(), name, name, self.domain().residue_characteristic(), self._power)
        return s

    def _repr_short(self):
        """
        Return a short string representation of this endomorphism.

        EXAMPLES::

            sage: K.<a> = Qq(5^3)
            sage: Frob = K.frobenius_endomorphism(); Frob
            Frobenius endomorphism on 5-adic Unramified Extension ... lifting a |--> a^5 on the residue field

            sage: Frob._repr_short()
            'Frob'
        """
        if self._power == 0:
            s = "Identity"
        elif self._power == 1:
            s = "Frob"
        else:
            s = "Frob^%s" % self._power
        return s

    cpdef Element _call_ (self, x):
        """
        TESTS::

            sage: K.<a> = Qq(5^3)
            sage: Frob = K.frobenius_endomorphism()
            sage: Frob(a) == a.frobenius()
            True
        """
        res = x
        for i in range(self._power):
            res = res.frobenius()
        return res

    def order(self):
        """
        Return the order of this endomorphism.

        EXAMPLES::

            sage: K.<a> = Qq(5^12)
            sage: Frob = K.frobenius_endomorphism()
            sage: Frob.order()
            12
            sage: (Frob^2).order()
            6
            sage: (Frob^9).order()
            4
        """
        if self._order == 0:
            return Infinity
        else:
            return Integer(self._order)

    def power(self):
        """
        Return the smallest integer `n` such that this endomorphism
        is the `n`-th power of the absolute (arithmetic) Frobenius.

        EXAMPLES::

            sage: K.<a> = Qq(5^12)
            sage: Frob = K.frobenius_endomorphism()
            sage: Frob.power()
            1
            sage: (Frob^9).power()
            9
            sage: (Frob^13).power()
            1
        """
        return self._power

    def __pow__(self,n,modulus):
        """
        Return the `n`-th iterate of this endomorphism.

        EXAMPLES::

            sage: K.<a> = Qq(5^12)
            sage: Frob = K.frobenius_endomorphism()
            sage: Frob^2
            Frobenius endomorphism on 5-adic Unramified Extension ... lifting a |--> a^(5^2) on the residue field

        The result is simplified if possible::

            sage: Frob^15
            Frobenius endomorphism on 5-adic Unramified Extension ... lifting a |--> a^(5^3) on the residue field
            sage: Frob^36
            Identity endomorphism of 5-adic Unramified Extension ...
        """
        return self.__class__(self.domain(), self.power()*n)

    def _composition(self, right):
        """
        Return ``self`` o ``right``.

        EXAMPLES::

            sage: K.<a> = Qq(5^12)
            sage: f = K.frobenius_endomorphism(); f
            Frobenius endomorphism on 5-adic Unramified Extension ... lifting a |--> a^5 on the residue field
            sage: g = K.frobenius_endomorphism(2); g
            Frobenius endomorphism on 5-adic Unramified Extension ... lifting a |--> a^(5^2) on the residue field
            sage: f * g
            Frobenius endomorphism on 5-adic Unramified Extension ... lifting a |--> a^(5^3) on the residue field

        The result is simplified if possible::

            sage: f = K.frobenius_endomorphism(9)
            sage: g = K.frobenius_endomorphism(10)
            sage: f * g
            Frobenius endomorphism on 5-adic Unramified Extension ... lifting a |--> a^(5^7) on the residue field
        """
        if isinstance(right,FrobeniusEndomorphism_padics):
            return self.__class__(self.domain(), self._power+right.power())
        else:
            return RingHomomorphism._composition(self,right)

    def is_injective(self):
        """
        Return ``True`` since any power of the Frobenius endomorphism
        over an unramified `p`-adic field is always injective.

        EXAMPLES::

            sage: K.<a> = Qq(5^3)
            sage: Frob = K.frobenius_endomorphism()
            sage: Frob.is_injective()
            True
        """
        return True

    def is_surjective(self):
        """
        Return ``True`` since any power of the Frobenius endomorphism
        over an unramified `p`-adic field is always surjective.

        EXAMPLES::

            sage: K.<a> = Qq(5^3)
            sage: Frob = K.frobenius_endomorphism()
            sage: Frob.is_surjective()
            True
        """
        return True

    def is_identity(self):
        """
        Return ``True`` if this morphism is the identity morphism.

        EXAMPLES::

            sage: K.<a> = Qq(5^3)
            sage: Frob = K.frobenius_endomorphism()
            sage: Frob.is_identity()
            False
            sage: (Frob^3).is_identity()
            True
        """
        return self.power() == 0

    def __hash__(self):
        r"""
        Return a hash of this morphism.

        It is the hash of ``(domain, codomain, ('Frob', power)``
        where ``power`` is the smallest integer `n` such that
        this morphism acts by `x \mapsto x^(p^n)` on the residue field.

        EXAMPLES::

            sage: K.<a> = Qq(5^3)
            sage: Frob = K.frobenius_endomorphism()
            sage: hash(Frob)  # indirect doctest, random
            2818440606874670810
        """
        domain = self.domain()
        codomain = self.codomain()
        return hash((domain, codomain, ('Frob', self._power)))

    cpdef _richcmp_(left, right, int op):
        """
        Compare ``left`` and ``right``.

        EXAMPLES::

            sage: K.<a> = Qq(5^3)
            sage: F = K.frobenius_endomorphism()
            sage: G = K.frobenius_endomorphism(4)

            sage: F == G
            True
        """
        if left is right:
            return rich_to_bool(op, 0)
        l_domain = left.domain()
        r_domain = right.domain()
        if l_domain != r_domain:
            return richcmp_not_equal(l_domain, r_domain, op)

        l_codomain = left.codomain()
        r_codomain = right.codomain()
        if l_codomain != r_codomain:
            return richcmp_not_equal(l_codomain, r_codomain, op)

        if isinstance(right, FrobeniusEndomorphism_padics):
            return richcmp(left._power, (<FrobeniusEndomorphism_padics>right)._power, op)

        try:
            for x in l_domain.gens():
                lx = left(x)
                rx = right(x)
                if lx != rx:
                    return richcmp_not_equal(lx, rx, op)
            return rich_to_bool(op, 0)
        except (AttributeError, NotImplementedError):
            raise NotImplementedError
