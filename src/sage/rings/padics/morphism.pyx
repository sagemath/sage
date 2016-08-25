"""
Frobenius endomorphisms on p-adic fields
"""

#*****************************************************************************
#       Copyright (C) 2013 Xavier Caruso <xavier.caruso@normalesup.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.integer cimport Integer
from sage.rings.infinity import Infinity

from sage.rings.ring import CommutativeRing
from sage.categories.homset import Hom
from sage.structure.element cimport Element

from sage.rings.morphism cimport RingHomomorphism
from padic_generic import pAdicGeneric

from sage.categories.morphism cimport Morphism


cdef class FrobeniusEndomorphism_padics(RingHomomorphism):
    """
    A class implementing Frobenius endomorphisms on padic fields.
    """    
    def __init__ (self,domain,n=1):
        """
        INPUT:
            
        -  ``domain`` -- an unramified padic field
            
        -  ``n`` -- an integer (default: 1)
            
        .. NOTE::
           
            `n` may be negative.
            
        OUTPUT:
            
        The `n`-th power of the absolute (arithmetic) Frobenius
        endomorphism on ``domain``
    
        TESTS::
    
            sage: from sage.rings.padics.morphism import FrobeniusEndomorphism_padics
            sage: K.<a> = Qq(5^3)
            sage: FrobeniusEndomorphism_padics(K)
            Frobenius endomorphism on Unramified Extension of 5-adic Field ... lifting a |--> a^5 on the residue field
            sage: FrobeniusEndomorphism_padics(K,2)
            Frobenius endomorphism on Unramified Extension of 5-adic Field ... lifting a |--> a^(5^2) on the residue field

            sage: FrobeniusEndomorphism_padics(K,a)
            Traceback (most recent call last):
            ...
            TypeError: n (=a + O(5^20)) is not an integer

            sage: K = Qp(5)
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
        if domain.e() != 1:
            raise TypeError("The domain must be unramified")
        try:
            n = Integer(n)
        except (ValueError, TypeError):
            raise TypeError("n (=%s) is not an integer" % n)

        self._degree = domain.f()
        self._power = n % self._degree
        self._order = self._degree / domain.degree().gcd(self._power)
        RingHomomorphism.__init__(self, Hom(domain, domain))


    def _repr_(self):
        """
        Return a string representation of this endomorphism.
            
        EXAMPLES::
            
            sage: K.<a> = Qq(5^3)
            sage: Frob = K.frobenius_endomorphism(); Frob
            Frobenius endomorphism on Unramified Extension of 5-adic Field ... lifting a |--> a^5 on the residue field
        
            sage: Frob._repr_()
            'Frobenius endomorphism on Unramified Extension of 5-adic Field ... lifting a |--> a^5 on the residue field'
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
            Frobenius endomorphism on Unramified Extension of 5-adic Field ... lifting a |--> a^5 on the residue field

            sage: Frob._repr_short()
            'Frob'
        """
        name = self.domain().variable_name()
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
            sage: Frob = K.frobenius_endomorphism();
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
        Return the smallest integer `n` such that this endormorphism
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
            Frobenius endomorphism on Unramified Extension of 5-adic Field ... lifting a |--> a^(5^2) on the residue field

        The result is simplified if possible::

            sage: Frob^15
            Frobenius endomorphism on Unramified Extension of 5-adic Field ... lifting a |--> a^(5^3) on the residue field
            sage: Frob^36
            Identity endomorphism of Unramified Extension of 5-adic Field ...
        """
        return self.__class__(self.domain(), self.power()*n)


    def _composition(self,right):
        """
        Return self o right.

        EXAMPLES::

            sage: K.<a> = Qq(5^12)
            sage: f = K.frobenius_endomorphism(); f
            Frobenius endomorphism on Unramified Extension of 5-adic Field ... lifting a |--> a^5 on the residue field
            sage: g = K.frobenius_endomorphism(2); g
            Frobenius endomorphism on Unramified Extension of 5-adic Field ... lifting a |--> a^(5^2) on the residue field
            sage: f * g
            Frobenius endomorphism on Unramified Extension of 5-adic Field ... lifting a |--> a^(5^3) on the residue field

        The result is simplified if possible::

            sage: f = K.frobenius_endomorphism(9)
            sage: g = K.frobenius_endomorphism(10)
            sage: f * g
            Frobenius endomorphism on Unramified Extension of 5-adic Field ... lifting a |--> a^(5^7) on the residue field
        """
        if isinstance(right,FrobeniusEndomorphism_padics):
            return self.__class__(self.domain(), self._power+right.power())
        else:
            return RingHomomorphism._composition(self,right)

    def is_injective(self):
        """
        Return true since any power of the Frobenius endomorphism
        over an unramified padic field is always injective.

        EXAMPLES::

            sage: K.<a> = Qq(5^3)
            sage: Frob = K.frobenius_endomorphism()
            sage: Frob.is_injective()
            True
        """
        return True


    def is_surjective(self):
        """
        Return true since any power of the Frobenius endomorphism
        over an unramified padic field is always surjective.

        EXAMPLES::

            sage: K.<a> = Qq(5^3)
            sage: Frob = K.frobenius_endomorphism()
            sage: Frob.is_surjective()
            True
        """
        return True


    def is_identity(self):
        """
        Return true if this morphism is the identity morphism.

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
        """
        Return a hash of this morphism.

        It is the hash of ``(domain, codomain, ('Frob', power)``
        where ``power`` is the smalles integer `n` such that 
        this morphism acts by `x \mapsto x^(p^n)` on the
        residue field
        """
        domain = self.domain()
        codomain = self.codomain()
        return hash((domain,codomain,('Frob',self._power)))

    cpdef int _cmp_(left, right) except -2:
        """
        Compare left and right
        """ 
        if left is right: return 0
        domain = left.domain()
        c = cmp(domain,right.domain())
        if c: return c
        c = cmp(left.codomain(),right.codomain())
        if c: return c
        if isinstance(right, FrobeniusEndomorphism_padics):
            return cmp(left._power, (<FrobeniusEndomorphism_padics>right)._power)
        try:
            gens = domain.gens()
            for x in gens:
                c = cmp(left(x),right(x))
                if c: return c
        except (AttributeError, NotImplementedError):
            raise NotImplementedError
