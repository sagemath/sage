"""
Finite field morphisms for prime fields

Special implementation for prime finite field of:

- embeddings of such field into general finite fields

- Frobenius endomorphisms (= identity with our assumptions)

.. SEEALSO::

    :mod:`sage.rings.finite_rings.hom_finite_field`

AUTHOR:

- Xavier Caruso (2012-06-29)
"""

#############################################################################
#    Copyright (C) 2012 Xavier Caruso <xavier.caruso@normalesup.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#****************************************************************************

from sage.categories.homset import Hom
from sage.structure.element cimport Element

from sage.rings.finite_rings.hom_finite_field cimport SectionFiniteFieldHomomorphism_generic
from sage.rings.finite_rings.hom_finite_field cimport FiniteFieldHomomorphism_generic
from sage.rings.finite_rings.hom_finite_field cimport FrobeniusEndomorphism_finite_field

from sage.rings.finite_rings.finite_field_base import FiniteField


cdef class SectionFiniteFieldHomomorphism_prime(SectionFiniteFieldHomomorphism_generic):
    cpdef Element _call_(self, x):
        try:
            return self._codomain._element_constructor_(x)
        except TypeError:
            raise ValueError("%s is not in the image of %s" % (x, self._inverse))


cdef class FiniteFieldHomomorphism_prime(FiniteFieldHomomorphism_generic):
    """
    A class implementing embeddings of prime finite fields into
    general finite fields.
    """
    def __init__(self, parent, im_gens=None, base_map=None, check=False, section_class=None):
        """
        TESTS::

            sage: from sage.rings.finite_rings.hom_prime_finite_field import FiniteFieldHomomorphism_prime
            sage: k = GF(3)
            sage: K.<T> = GF(3^4)                                                       # needs sage.rings.finite_rings
            sage: f = FiniteFieldHomomorphism_prime(Hom(k, K)); f                       # needs sage.rings.finite_rings
            Ring morphism:
              From: Finite Field of size 3
              To:   Finite Field in T of size 3^4
              Defn: 1 |--> 1

            sage: k.<t> = GF(3^2)                                                       # needs sage.rings.finite_rings
            sage: K.<T> = GF(3^4)                                                       # needs sage.rings.finite_rings
            sage: f = FiniteFieldHomomorphism_prime(Hom(k, K)); f                       # needs sage.rings.finite_rings
            Traceback (most recent call last):
            ...
            TypeError: The domain is not a finite prime field
        """
        domain = parent.domain()
        if not isinstance(domain, FiniteField) or not domain.is_prime_field():
            raise TypeError("The domain is not a finite prime field")
        if section_class is None:
            section_class = SectionFiniteFieldHomomorphism_prime
        FiniteFieldHomomorphism_generic.__init__(self, parent, im_gens, base_map=base_map,
                                                 check=check, section_class=section_class)

    cpdef Element _call_(self, x):
        """
        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: from sage.rings.finite_rings.hom_prime_finite_field import FiniteFieldHomomorphism_prime
            sage: k = GF(3)
            sage: K.<t> = GF(3^5)
            sage: f = FiniteFieldHomomorphism_prime(Hom(k, K))
            sage: a = f(4); a
            1
            sage: a.parent()
            Finite Field in t of size 3^5
        """
        return self._codomain._element_constructor_(x)


cdef class FrobeniusEndomorphism_prime(FrobeniusEndomorphism_finite_field):
    """
    A class implementing Frobenius endomorphism on prime finite
    fields (i.e. identity map :-).
    """
    def __init__(self, domain, power=1):
        if not isinstance(domain, FiniteField) or not domain.is_prime_field():
            raise TypeError("The domain is not a finite prime field")
        FrobeniusEndomorphism_finite_field.__init__(self, Hom(domain, domain))
        self._order = 1
        self._power = 0

    cpdef Element _call_(self, x):
        """
        TESTS::

            sage: k = GF(5)
            sage: Frob = k.frobenius_endomorphism()
            sage: Frob(2)
            2
        """
        return x

    def _composition(self, right):
        """
        Return ``self`` o ``right``.

        It is always ``right``, since ``self`` is always identity because
        the domain is a prime field.
        """
        return right

    def __pow__(self, n, modulus):
        """
        Return the `n`-th iterate of this endomorphism
        (that is the identity since the domain is a prime
        field).
        """
        return self

    def fixed_field(self):
        """
        Return the fixed field of ``self``.

        OUTPUT:

        - a tuple `(K, e)`, where `K` is the subfield of the domain
          consisting of elements fixed by ``self`` and `e` is an
          embedding of `K` into the domain.

        .. NOTE::

            Since here the domain is a prime field, the subfield
            is the same prime field and the embedding is necessarily
            the identity map.

        EXAMPLES::

            sage: k.<t> = GF(5)
            sage: f = k.frobenius_endomorphism(2); f
            Identity endomorphism of Finite Field of size 5
            sage: kfixed, embed = f.fixed_field()

            sage: kfixed == k
            True
            sage: [ embed(x) == x for x in kfixed ]
            [True, True, True, True, True]
        """
        return self.domain(), self.domain().hom(self.domain())
