# sage_setup: distribution = sagemath-linbox
# distutils: libraries = givaro NTL_LIBRARIES gmp m
# distutils: extra_compile_args = NTL_CFLAGS
# distutils: include_dirs = NTL_INCDIR
# distutils: library_dirs = NTL_LIBDIR
# distutils: extra_link_args = NTL_LIBEXTRA
# distutils: language = c++
"""
Givaro finite field morphisms

Special implementation for givaro finite fields of:

- embeddings between finite fields

- frobenius endomorphisms

SEEALSO::

    :mod:`sage.rings.finite_rings.hom_finite_field`

AUTHOR:

- Xavier Caruso (2012-06-29)
"""

# ###########################################################################
#    Copyright (C) 2012 Xavier Caruso <xavier.caruso@normalesup.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
# ***************************************************************************


from sage.rings.finite_rings.finite_field_constructor import FiniteField

from sage.rings.finite_rings.hom_finite_field cimport SectionFiniteFieldHomomorphism_generic
from sage.rings.finite_rings.hom_finite_field cimport FiniteFieldHomomorphism_generic
from sage.rings.finite_rings.hom_finite_field cimport FrobeniusEndomorphism_finite_field

from sage.rings.finite_rings.hom_prime_finite_field cimport FiniteFieldHomomorphism_prime

from sage.categories.homset import Hom
from sage.structure.element cimport Element

from sage.rings.finite_rings.finite_field_givaro import FiniteField_givaro
from sage.rings.finite_rings.element_givaro cimport FiniteField_givaroElement
#from element_givaro cimport make_FiniteField_givaroElement

from sage.structure.parent cimport Parent
from sage.rings.finite_rings.element_givaro cimport Cache_givaro


cdef class SectionFiniteFieldHomomorphism_givaro(SectionFiniteFieldHomomorphism_generic):
    def __init__(self, inverse):
        """
        TESTS::

            sage: from sage.rings.finite_rings.hom_finite_field_givaro import FiniteFieldHomomorphism_givaro
            sage: k.<t> = GF(3^2)
            sage: K.<T> = GF(3^4)
            sage: f = FiniteFieldHomomorphism_givaro(Hom(k, K))
            sage: g = f.section(); g # random
            Section of Ring morphism:
              From: Finite Field in t of size 3^2
              To:   Finite Field in T of size 3^4
              Defn: t |--> 2*T^3 + 2*T^2 + 1
            sage: a = k.random_element()
            sage: b = k.random_element()
            sage: g(f(a) + f(b)) == g(f(a)) + g(f(b)) == a + b
            True
        """
        if not isinstance(inverse, FiniteFieldHomomorphism_givaro):
            raise TypeError("The given map is not an instance of FiniteFieldHomomorphism_givaro")
        SectionFiniteFieldHomomorphism_generic.__init__(self, inverse)

        cdef long inverse_power = (<FiniteFieldHomomorphism_givaro>inverse)._power
        cdef long order = self.domain().cardinality() - 1
        self._order_codomain = self.codomain().cardinality() - 1

        # Compute a = gcd(inverse_power, order)
        # and solve inverse_power*x = a (mod order)
        cdef long a = inverse_power, b = order
        cdef unsigned long q
        cdef long x = 1, y = 0
        cdef long sb, sy
        while b != 0:
            q = a // b
            sb = b
            b = a - q * b
            a = sb
            sy = y
            y = x - q * y
            x = sy

        self._gcd = a
        if x < 0:
            x += order
        self._power = x % self._order_codomain

        self._codomain_cache = (<FiniteField_givaroElement>(self._codomain.gen()))._cache

    cpdef Element _call_(self, x):
        """
        TESTS::

            sage: from sage.rings.finite_rings.hom_finite_field_givaro import FiniteFieldHomomorphism_givaro
            sage: k.<t> = GF(3^2)
            sage: K.<T> = GF(3^4)
            sage: f = FiniteFieldHomomorphism_givaro(Hom(k, K))
            sage: g = f.section()
            sage: g(f(t+1))
            t + 1

            sage: g(T)
            Traceback (most recent call last):
            ...
            ValueError: T is not in the image of Ring morphism:
              From: Finite Field in t of size 3^2
              To:   Finite Field in T of size 3^4
              Defn: t |--> ...
        """
        if x.parent() != self.domain():
            raise TypeError("%s is not in %s" % (x, self.domain()))
        cdef FiniteField_givaroElement y = <FiniteField_givaroElement?>x
        if y._cache.objectptr.isZero(y.element):
            return make_FiniteField_givaroElement(self._codomain_cache, self._codomain_cache.objectptr.zero)
        if y._cache.objectptr.isOne(y.element):
            return make_FiniteField_givaroElement(self._codomain_cache, self._codomain_cache.objectptr.one)
        cdef int log = y.element
        cdef int q = log / self._gcd
        if log == q*self._gcd:
            q = (q*self._power) % self._order_codomain
            return make_FiniteField_givaroElement(self._codomain_cache, q)
        else:
            raise ValueError("%s is not in the image of %s" % (x, self._inverse))


cdef class FiniteFieldHomomorphism_givaro(FiniteFieldHomomorphism_generic):
    def __init__(self, parent, im_gens=None, check=False):
        """
        TESTS::

            sage: from sage.rings.finite_rings.hom_finite_field_givaro import FiniteFieldHomomorphism_givaro
            sage: k.<t> = GF(3^2)
            sage: K.<T> = GF(3^4)
            sage: f = FiniteFieldHomomorphism_givaro(Hom(k, K)); f # random
            Ring morphism:
              From: Finite Field in t of size 3^2
              To:   Finite Field in T of size 3^4
              Defn: t |--> 2*T^3 + 2*T^2 + 1
            sage: a = k.random_element()
            sage: b = k.random_element()
            sage: f(a) + f(b) == f(a + b)
            True

            sage: k.<t> = GF(3^10)
            sage: K.<T> = GF(3^20)
            sage: f = FiniteFieldHomomorphism_givaro(Hom(k, K)); f
            Traceback (most recent call last):
            ...
            TypeError: The codomain is not an instance of FiniteField_givaro
        """
        domain = parent.domain()
        codomain = parent.codomain()
        if not isinstance(domain, FiniteField_givaro):
            raise TypeError("The domain is not an instance of FiniteField_givaro")
        if not isinstance(codomain, FiniteField_givaro):
            raise TypeError("The codomain is not an instance of FiniteField_givaro")

        FiniteFieldHomomorphism_generic.__init__(self, parent, im_gens, check=check,
                                                 section_class=SectionFiniteFieldHomomorphism_givaro)

        cdef Cache_givaro domain_cache = (<FiniteField_givaroElement>(domain.gen()))._cache
        self._codomain_cache = (<FiniteField_givaroElement>(codomain.gen()))._cache

        cdef FiniteField_givaroElement g = make_FiniteField_givaroElement(domain_cache, 1)
        cdef FiniteField_givaroElement G = FiniteFieldHomomorphism_generic._call_(self, g)
        self._power = G.element

        self._order_domain = domain.cardinality() - 1
        self._order_codomain = codomain.cardinality() - 1

    cpdef Element _call_(self, x):
        """
        TESTS::

            sage: from sage.rings.finite_rings.hom_finite_field_givaro import FiniteFieldHomomorphism_givaro
            sage: k.<t> = GF(3^2)
            sage: K.<T> = GF(3^4)
            sage: f = FiniteFieldHomomorphism_givaro(Hom(k, K))
            sage: f(t) # random
            2*T^3 + 2*T^2 + 1
            sage: f(t) == f.im_gens()[0]
            True
        """
        if x.parent() != self.domain():
            raise TypeError("%s is not in %s" % (x, self.domain()))
        cdef FiniteField_givaroElement y = <FiniteField_givaroElement?>x
        if y._cache.objectptr.isZero(y.element):
            return make_FiniteField_givaroElement(self._codomain_cache, self._codomain_cache.objectptr.zero)
        if y._cache.objectptr.isOne(y.element):
            return make_FiniteField_givaroElement(self._codomain_cache, self._codomain_cache.objectptr.one)
        cdef int log = y.element
        log = (log*self._power) % self._order_codomain
        return make_FiniteField_givaroElement(self._codomain_cache, log)


cdef class FrobeniusEndomorphism_givaro(FrobeniusEndomorphism_finite_field):
    def __init__(self, domain, power=1):
        """
        TESTS::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism(); Frob
            Frobenius endomorphism t |--> t^5 on Finite Field in t of size 5^3
            sage: type(Frob)
            <class 'sage.rings.finite_rings.hom_finite_field_givaro.FrobeniusEndomorphism_givaro'>

            sage: k.<t> = GF(5^20)
            sage: Frob = k.frobenius_endomorphism(); Frob
            Frobenius endomorphism t |--> t^5 on Finite Field in t of size 5^20
            sage: type(Frob)
            <class 'sage.rings.finite_rings.hom_finite_field.FrobeniusEndomorphism_finite_field'>
        """
        if not isinstance(domain, FiniteField_givaro):
            raise TypeError("The domain is not an instance of FiniteField_givaro")
        FrobeniusEndomorphism_finite_field.__init__(self, domain, power)

    def fixed_field(self):
        """
        Return the fixed field of ``self``.

        OUTPUT:

        - a tuple `(K, e)`, where `K` is the subfield of the domain
          consisting of elements fixed by ``self`` and `e` is an
          embedding of `K` into the domain.

        .. NOTE::

            The name of the variable used for the subfield (if it
            is not a prime subfield) is suffixed by ``_fixed``.

        EXAMPLES::

            sage: k.<t> = GF(5^6)
            sage: f = k.frobenius_endomorphism(2)
            sage: kfixed, embed = f.fixed_field()
            sage: kfixed
            Finite Field in t_fixed of size 5^2
            sage: embed # random
            Ring morphism:
              From: Finite Field in t_fixed of size 5^2
              To:   Finite Field in t of size 5^6
              Defn: t_fixed |--> 4*t^5 + 2*t^4 + 4*t^2 + t

            sage: tfixed = kfixed.gen()
            sage: embed(tfixed) # random
            4*t^5 + 2*t^4 + 4*t^2 + t
        """
        if self._degree_fixed == 1:
            k = FiniteField(self.domain().characteristic())
            f = FiniteFieldHomomorphism_prime(Hom(k, self.domain()))
        else:
            k = FiniteField(self.domain().characteristic()**self._degree_fixed,
                            name=self.domain().variable_name() + "_fixed")
            f = FiniteFieldHomomorphism_givaro(Hom(k, self.domain()))
        return k, f


# copied from element_givaro.pyx
cdef inline FiniteField_givaroElement make_FiniteField_givaroElement(Cache_givaro cache, int x):
    cdef FiniteField_givaroElement y

    if cache._has_array:
        return <FiniteField_givaroElement>cache._array[x]
    else:
        y = FiniteField_givaroElement.__new__(FiniteField_givaroElement)
        y._parent = <Parent> cache.parent
        y._cache = cache
        y.element = x
        return y
