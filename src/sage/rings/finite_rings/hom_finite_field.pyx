# sage.doctest: needs sage.rings.finite_rings
"""
Finite field morphisms

This file provides several classes implementing:

- embeddings between finite fields

- Frobenius isomorphism on finite fields

EXAMPLES::

    sage: from sage.rings.finite_rings.hom_finite_field import FiniteFieldHomomorphism_generic

Construction of an embedding::

    sage: k.<t> = GF(3^7)
    sage: K.<T> = GF(3^21)
    sage: f = FiniteFieldHomomorphism_generic(Hom(k, K)); f # random
    Ring morphism:
      From: Finite Field in t of size 3^7
      To:   Finite Field in T of size 3^21
      Defn: t |--> T^20 + 2*T^18 + T^16 + 2*T^13 + T^9 + 2*T^8 + T^7 + T^6 + T^5 + T^3 + 2*T^2 + T

    sage: f(t) # random
    T^20 + 2*T^18 + T^16 + 2*T^13 + T^9 + 2*T^8 + T^7 + T^6 + T^5 + T^3 + 2*T^2 + T
    sage: f(t) == f.im_gens()[0]
    True

The map `f` has a method ``section`` which returns a partially defined
map which is the inverse of `f` on the image of `f`::

    sage: g = f.section(); g # random
    Section of Ring morphism:
      From: Finite Field in t of size 3^7
      To:   Finite Field in T of size 3^21
      Defn: t |--> T^20 + 2*T^18 + T^16 + 2*T^13 + T^9 + 2*T^8 + T^7 + T^6 + T^5 + T^3 + 2*T^2 + T
    sage: a = k.random_element()
    sage: g(f(a)) == a
    True
    sage: g(T)
    Traceback (most recent call last):
    ...
    ValueError: T is not in the image of Ring morphism:
      From: Finite Field in t of size 3^7
      To:   Finite Field in T of size 3^21
      Defn: ...

There is no embedding of `GF(5^6)` into `GF(5^11)`::

    sage: k.<t> = GF(5^6)
    sage: K.<T> = GF(5^11)
    sage: FiniteFieldHomomorphism_generic(Hom(k, K))
    Traceback (most recent call last):
    ...
    ValueError: No embedding of Finite Field in t of size 5^6 into Finite Field in T of size 5^11


Construction of Frobenius endomorphisms::

    sage: k.<t> = GF(7^14)
    sage: Frob = k.frobenius_endomorphism(); Frob
    Frobenius endomorphism t |--> t^7 on Finite Field in t of size 7^14
    sage: Frob(t)
    t^7

Some basic arithmetics is supported::

    sage: Frob^2
    Frobenius endomorphism t |--> t^(7^2) on Finite Field in t of size 7^14
    sage: f = k.frobenius_endomorphism(7); f
    Frobenius endomorphism t |--> t^(7^7) on Finite Field in t of size 7^14
    sage: f*Frob
    Frobenius endomorphism t |--> t^(7^8) on Finite Field in t of size 7^14

    sage: Frob.order()
    14
    sage: f.order()
    2

Note that simplifications are made automatically::

    sage: Frob^16
    Frobenius endomorphism t |--> t^(7^2) on Finite Field in t of size 7^14
    sage: Frob^28
    Identity endomorphism of Finite Field in t of size 7^14

And that comparisons work::

    sage: Frob == Frob^15
    True
    sage: Frob^14 == Hom(k, k).identity()
    True

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

from sage.rings.integer cimport Integer

from sage.categories.homset import Hom
from sage.structure.element cimport Element

from sage.rings.finite_rings.finite_field_base import FiniteField as FiniteField_base
from sage.rings.morphism cimport RingHomomorphism, RingHomomorphism_im_gens, FrobeniusEndomorphism_generic
from sage.rings.finite_rings.finite_field_constructor import FiniteField

from sage.categories.map cimport Section
from sage.categories.morphism cimport Morphism

from sage.misc.cachefunc import cached_method


cdef class SectionFiniteFieldHomomorphism_generic(Section):
    """
    A class implementing sections of embeddings between finite fields.
    """
    cpdef Element _call_(self, x):  # Not optimized
        """
        TESTS::

            sage: from sage.rings.finite_rings.hom_finite_field import FiniteFieldHomomorphism_generic
            sage: k.<t> = GF(3^7)
            sage: K.<T> = GF(3^21)
            sage: f = FiniteFieldHomomorphism_generic(Hom(k, K))
            sage: g = f.section()
            sage: a = k.random_element()
            sage: g(f(a)) == a
            True

            sage: g(T)
            Traceback (most recent call last):
            ...
            ValueError: T is not in the image of Ring morphism:
              From: Finite Field in t of size 3^7
              To:   Finite Field in T of size 3^21
              Defn: t |--> ...
        """
        for root, _ in x.minimal_polynomial().roots(ring=self.codomain()):
            if self._inverse(root) == x:
                return root
        raise ValueError("%s is not in the image of %s" % (x, self._inverse))

    def _repr_(self):
        """
        Return a string representation of this section.

        EXAMPLES::

            sage: from sage.rings.finite_rings.hom_finite_field import FiniteFieldHomomorphism_generic
            sage: k.<t> = GF(3^7)
            sage: K.<T> = GF(3^21)
            sage: f = FiniteFieldHomomorphism_generic(Hom(k, K))
            sage: g = f.section()
            sage: g._repr_() # random
            'Section of Ring morphism:\n  From: Finite Field in t of size 3^7\n  To:   Finite Field in T of size 3^21\n  Defn: t |--> T^20 + 2*T^18 + T^16 + 2*T^13 + T^9 + 2*T^8 + T^7 + T^6 + T^5 + T^3 + 2*T^2 + T'
        """
        return "Section of %s" % self._inverse

    def _latex_(self):
        r"""
        Return a latex representation of this section.

        EXAMPLES::

            sage: from sage.rings.finite_rings.hom_finite_field import FiniteFieldHomomorphism_generic
            sage: k.<t> = GF(3^7)
            sage: K.<T> = GF(3^21)
            sage: f = FiniteFieldHomomorphism_generic(Hom(k, K))
            sage: g = f.section()
            sage: g._latex_()
            '\\verb"Section of "\\Bold{F}_{3^{7}} \\hookrightarrow \\Bold{F}_{3^{21}}'
        """
        return '\\verb"Section of "' + self._inverse._latex_()


cdef class FiniteFieldHomomorphism_generic(RingHomomorphism_im_gens):
    """
    A class implementing embeddings between finite fields.

    TESTS::

        sage: from sage.rings.finite_rings.hom_finite_field import FiniteFieldHomomorphism_generic
        sage: k.<t> = GF(3^7)
        sage: K.<T> = GF(3^21)
        sage: f = FiniteFieldHomomorphism_generic(Hom(k, K))
        sage: TestSuite(f).run()
    """
    def __init__(self, parent, im_gens=None, base_map=None, check=True, section_class=None):
        """
        TESTS::

            sage: from sage.rings.finite_rings.hom_finite_field import FiniteFieldHomomorphism_generic
            sage: k.<t> = GF(3^7)
            sage: K.<T> = GF(3^21)
            sage: f = FiniteFieldHomomorphism_generic(Hom(k, K)); f # random
            Ring morphism:
              From: Finite Field in t of size 3^7
              To:   Finite Field in T of size 3^21
              Defn: t |--> T^20 + 2*T^18 + T^16 + 2*T^13 + T^9 + 2*T^8 + T^7 + T^6 + T^5 + T^3 + 2*T^2 + T
            sage: a = k.random_element()
            sage: b = k.random_element()
            sage: f(a) + f(b) == f(a + b)
            True

            sage: k.<t> = GF(3^6)
            sage: K.<t> = GF(3^9)
            sage: FiniteFieldHomomorphism_generic(Hom(k, K))
            Traceback (most recent call last):
            ...
            ValueError: No embedding of Finite Field in t of size 3^6 into Finite Field in t of size 3^9

            sage: FiniteFieldHomomorphism_generic(Hom(ZZ, QQ))
            Traceback (most recent call last):
            ...
            TypeError: The domain is not a finite field or does not provide the required interface for finite fields

            sage: R.<x> = k[]
            sage: FiniteFieldHomomorphism_generic(Hom(k, R))
            Traceback (most recent call last):
            ...
            TypeError: The codomain is not a finite field or does not provide the required interface for finite fields
        """
        domain = parent.domain()
        codomain = parent.codomain()
        if not isinstance(domain, FiniteField_base):
            raise TypeError("The domain is not a finite field or does not provide the required interface for finite fields")
        if not isinstance(codomain, FiniteField_base):
            raise TypeError("The codomain is not a finite field or does not provide the required interface for finite fields")
        if domain.characteristic() != codomain.characteristic() or codomain.degree() % domain.degree() != 0:
            raise ValueError("No embedding of %s into %s" % (domain, codomain))
        if im_gens is None:
            im_gens = domain.modulus().any_root(codomain)
            check=False
        RingHomomorphism_im_gens.__init__(self, parent, im_gens, base_map=base_map, check=check)
        if section_class is None:
            self._section_class = SectionFiniteFieldHomomorphism_generic
        else:
            self._section_class = section_class

    def __copy__(self):
        """
        Return a copy of this map.

        TESTS::

            sage: from sage.rings.finite_rings.hom_finite_field import FiniteFieldHomomorphism_generic
            sage: k.<t> = GF(3^7)
            sage: K.<T> = GF(3^21)
            sage: f = FiniteFieldHomomorphism_generic(Hom(k, K))
            sage: g = copy(f)
            sage: g.section()(g(t)) == f.section()(f(t))
            True

        ::

            sage: F = GF(2)
            sage: E = GF(4)
            sage: phi = E.coerce_map_from(F); phi
            Ring morphism:
              From: Finite Field of size 2
              To:   Finite Field in z2 of size 2^2
              Defn: 1 |--> 1
            sage: phi.section()
            Section of Ring morphism:
              From: Finite Field of size 2
              To:   Finite Field in z2 of size 2^2
              Defn: 1 |--> 1
        """
        cdef FiniteFieldHomomorphism_generic out = super().__copy__()
        out._section_class = self._section_class
        return out

    def _latex_(self):
        r"""
        Return a latex representation of this embedding.

        EXAMPLES::

            sage: from sage.rings.finite_rings.hom_finite_field import FiniteFieldHomomorphism_generic
            sage: k.<t> = GF(3^7)
            sage: K.<T> = GF(3^21)
            sage: f = FiniteFieldHomomorphism_generic(Hom(k, K))
            sage: f._latex_()
            '\\Bold{F}_{3^{7}} \\hookrightarrow \\Bold{F}_{3^{21}}'
        """
        return self.domain()._latex_() + " \\hookrightarrow " + self.codomain()._latex_()

    cpdef Element _call_(self, x):
        """
        TESTS::

            sage: from sage.rings.finite_rings.hom_finite_field import FiniteFieldHomomorphism_generic
            sage: k.<t> = GF(3^3)
            sage: K.<T> = GF(3^9)
            sage: f = FiniteFieldHomomorphism_generic(Hom(k, K))
            sage: f(t) # random
            2*T^6 + 2*T^4 + T^2 + T

            sage: a = k.random_element()
            sage: b = k.random_element()
            sage: f(a+b) == f(a) + f(b)
            True
            sage: f(a*b) == f(a) * f(b)
            True
        """
        f = x.polynomial()
        bm = self.base_map()
        if bm is not None:
            f = f.map_coefficients(bm)
        return f(self.im_gens()[0])

    def is_injective(self):
        """
        Return ``True`` since a embedding between finite fields is
        always injective.

        EXAMPLES::

            sage: from sage.rings.finite_rings.hom_finite_field import FiniteFieldHomomorphism_generic
            sage: k.<t> = GF(3^3)
            sage: K.<T> = GF(3^9)
            sage: f = FiniteFieldHomomorphism_generic(Hom(k, K))
            sage: f.is_injective()
            True
        """
        return True

    def is_surjective(self):
        """
        Return ``True`` if this embedding is surjective (and hence an
        isomorphism.

        EXAMPLES::

            sage: from sage.rings.finite_rings.hom_finite_field import FiniteFieldHomomorphism_generic
            sage: k.<t> = GF(3^3)
            sage: K.<T> = GF(3^9)
            sage: f = FiniteFieldHomomorphism_generic(Hom(k, K))
            sage: f.is_surjective()
            False
            sage: g = FiniteFieldHomomorphism_generic(Hom(k, k))
            sage: g.is_surjective()
            True
        """
        return self.domain().cardinality() == self.codomain().cardinality()

    @cached_method
    def section(self):
        """
        Return the ``inverse`` of this embedding.

        It is a partially defined map whose domain is the codomain
        of the embedding, but which is only defined on the image of
        the embedding.

        EXAMPLES::

            sage: from sage.rings.finite_rings.hom_finite_field import FiniteFieldHomomorphism_generic
            sage: k.<t> = GF(3^7)
            sage: K.<T> = GF(3^21)
            sage: f = FiniteFieldHomomorphism_generic(Hom(k, K))
            sage: g = f.section(); g # random
            Section of Ring morphism:
              From: Finite Field in t of size 3^7
              To:   Finite Field in T of size 3^21
              Defn: t |--> T^20 + 2*T^18 + T^16 + 2*T^13 + T^9 + 2*T^8 + T^7 + T^6 + T^5 + T^3 + 2*T^2 + T
            sage: a = k.random_element()
            sage: b = k.random_element()
            sage: g(f(a) + f(b)) == a + b
            True
            sage: g(T)
            Traceback (most recent call last):
            ...
            ValueError: T is not in the image of Ring morphism:
              From: Finite Field in t of size 3^7
              To:   Finite Field in T of size 3^21
              Defn: ...
        """
        if self.base_map() is not None:
            raise NotImplementedError
        return self._section_class(self)

    def _inverse_image_element(self, b):
        """
        Return the unique ``a`` such that ``self(a) = b`` if one such exists.

        This method is simply a shorthand for calling the map returned by
        ``self.section()`` on ``b``.

        EXAMPLES::

            sage: k.<t> = GF(3^7)
            sage: K.<T>, f = k.extension(3, map=True)
            sage: b = f(t^2); b
            2*T^20 + 2*T^19 + T^18 + T^15 + 2*T^14 + 2*T^13 + 2*T^12 + T^8 + 2*T^6 + T^5 + 2*T^4 + T^3 + 2*T^2 + T
            sage: f.inverse_image(b)
            t^2
            sage: f.inverse_image(T)
            Traceback (most recent call last):
            ...
            ValueError: T is not in the image of Ring morphism:
            From: Finite Field in t of size 3^7
            To:   Finite Field in T of size 3^21
            Defn: t |--> T^20 + 2*T^18 + T^16 + 2*T^13 + T^9 + 2*T^8 + T^7 + T^6 + T^5 + T^3 + 2*T^2 + T
        """
        return self.section()(b)

    def __hash__(self):
        r"""
        Return a hash of this morphism.

        TESTS::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: embed = Frob.fixed_field()[1]
            sage: hash(embed) # random
            -2441354824160407762
        """
        return Morphism.__hash__(self)

    cdef dict _extra_slots(self):
        r"""
        Helper function for copying and pickling.

        TESTS::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: embed = Frob.fixed_field()[1]
            sage: embed.__reduce__()  # indirect doctest
            (<built-in function unpickle_map>,
             (<class 'sage.rings.finite_rings.hom_prime_finite_field.FiniteFieldHomomorphism_prime'>,
              Set of field embeddings from Finite Field of size 5 to Finite Field in t of size 5^3,
              {},
              {'__im_gens': [1],
               '_base_map': None,
               '_codomain': Finite Field in t of size 5^3,
               '_domain': Finite Field of size 5,
               '_is_coercion': False,
               '_lift': None,
               '_repr_type_str': None,
               '_section_class': <class 'sage.rings.finite_rings.hom_prime_finite_field.SectionFiniteFieldHomomorphism_prime'>}))
        """
        cdef dict slots
        slots = RingHomomorphism_im_gens._extra_slots(self)
        slots['_section_class'] = self._section_class
        return slots

    cdef _update_slots(self, dict slots):
        r"""
        Helper function for copying and pickling.

        TESTS::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: embed = Frob.fixed_field()[1]
            sage: f = loads(dumps(embed))
            sage: f == embed
            True
            sage: f.section()
            Section of Ring morphism:
              From: Finite Field of size 5
              To:   Finite Field in t of size 5^3
              Defn: 1 |--> 1
        """
        RingHomomorphism_im_gens._update_slots(self, slots)
        self._section_class = slots['_section_class']


cdef class FrobeniusEndomorphism_finite_field(FrobeniusEndomorphism_generic):
    """
    A class implementing Frobenius endomorphisms on finite fields.

    TESTS::

        sage: k.<a> = GF(7^11)
        sage: Frob = k.frobenius_endomorphism(5)
        sage: TestSuite(Frob).run()
    """
    def __init__(self, domain, n=1):
        """
        INPUT:

        - ``domain`` -- a finite field

        - ``n`` -- integer (default: 1)

        .. NOTE::

            `n` may be negative.

        OUTPUT:

        The `n`-th power of the absolute (arithmetic) Frobenius
        endomorphism on ``domain``

        TESTS::

            sage: from sage.rings.finite_rings.hom_finite_field import FrobeniusEndomorphism_finite_field
            sage: k.<t> = GF(5^3)
            sage: FrobeniusEndomorphism_finite_field(k)
            Frobenius endomorphism t |--> t^5 on Finite Field in t of size 5^3
            sage: FrobeniusEndomorphism_finite_field(k, 2)
            Frobenius endomorphism t |--> t^(5^2) on Finite Field in t of size 5^3

            sage: FrobeniusEndomorphism_finite_field(k, t)
            Traceback (most recent call last):
            ...
            TypeError: n (=t) is not an integer

            sage: FrobeniusEndomorphism_finite_field(k['x'])
            Traceback (most recent call last):
            ...
            TypeError: The domain is not a finite field or does not provide the required interface for finite fields
        """
        if not isinstance(domain, FiniteField_base):
            raise TypeError("The domain is not a finite field or does not provide the required interface for finite fields")
        try:
            n = Integer(n)
        except TypeError:
            raise TypeError("n (=%s) is not an integer" % n)

        self._degree = domain.degree()
        self._power = n % self._degree
        self._degree_fixed = domain.degree().gcd(self._power)
        self._order = self._degree / self._degree_fixed
        self._q = domain.characteristic() ** self._power
        RingHomomorphism.__init__(self, Hom(domain, domain))

    def _repr_(self):
        """
        Return a string representation of this endomorphism.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism(); Frob
            Frobenius endomorphism t |--> t^5 on Finite Field in t of size 5^3

            sage: Frob._repr_()
            'Frobenius endomorphism t |--> t^5 on Finite Field in t of size 5^3'
        """
        name = self.domain().variable_name()
        if self._power == 0:
            s = "Identity endomorphism of"
        elif self._power == 1:
            s = "Frobenius endomorphism %s |--> %s^%s on" % (name, name, self.domain().characteristic())
        else:
            s = "Frobenius endomorphism %s |--> %s^(%s^%s) on" % (name, name, self.domain().characteristic(), self._power)
        s += " %s" % self.domain()
        return s

    def _repr_short(self):
        """
        Return a short string representation of this endomorphism.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism(); Frob
            Frobenius endomorphism t |--> t^5 on Finite Field in t of size 5^3

            sage: Frob._repr_short()
            't |--> t^5'
        """
        name = self.domain().variable_name()
        if self._power == 0:
            s = "Identity"
        elif self._power == 1:
            s = "%s |--> %s^%s" % (name, name, self.domain().characteristic())
        else:
            s = "%s |--> %s^(%s^%s)" % (name, name, self.domain().characteristic(), self._power)
        return s

    def _latex_(self):
        r"""
        Return a latex representation of this endomorphism.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: Frob._latex_()
            't \\mapsto t^{5}'
        """
        try:
            name = self.domain().latex_variable_names()[0]
        except IndexError:
            name = "x"
        if self._power == 0:
            s = '\\verb"id"'
        elif self._power == 1:
            s = "%s \\mapsto %s^{%s}" % (name, name, self.domain().characteristic())
        else:
            s = "%s \\mapsto %s^{%s^{%s}}" % (name, name, self.domain().characteristic(), self._power)
        return s

    cpdef Element _call_(self, x):
        """
        TESTS::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: Frob(t)
            2*t^2 + 4*t + 4
            sage: Frob(t) == t^5
            True
        """
        if self.is_identity():
            return x
        else:
            return x.pth_power(self._power)

    def order(self):
        """
        Return the order of this endomorphism.

        EXAMPLES::

            sage: k.<t> = GF(5^12)
            sage: Frob = k.frobenius_endomorphism()
            sage: Frob.order()
            12
            sage: (Frob^2).order()
            6
            sage: (Frob^9).order()
            4
        """
        if self._order == 0:
            from sage.rings.infinity import Infinity
            return Infinity
        else:
            return Integer(self._order)

    def power(self):
        """
        Return an integer `n` such that this endomorphism
        is the `n`-th power of the absolute (arithmetic)
        Frobenius.

        EXAMPLES::

            sage: k.<t> = GF(5^12)
            sage: Frob = k.frobenius_endomorphism()
            sage: Frob.power()
            1
            sage: (Frob^9).power()
            9
            sage: (Frob^13).power()
            1
        """
        return self._power

    def __pow__(self, n, modulus):
        """
        Return the `n`-th iterate of this endomorphism.

        EXAMPLES::

            sage: k.<t> = GF(5^12)
            sage: Frob = k.frobenius_endomorphism(); Frob
            Frobenius endomorphism t |--> t^5 on Finite Field in t of size 5^12
            sage: Frob^2
            Frobenius endomorphism t |--> t^(5^2) on Finite Field in t of size 5^12

        The result is simplified if possible::

            sage: Frob^15
            Frobenius endomorphism t |--> t^(5^3) on Finite Field in t of size 5^12
            sage: Frob^36
            Identity endomorphism of Finite Field in t of size 5^12
        """
        return self.__class__(self.domain(), self.power()*n)

    @cached_method
    def inverse(self):
        """
        Return the inverse of this Frobenius endomorphism.

        EXAMPLES::

            sage: k.<a> = GF(7^11)
            sage: f = k.frobenius_endomorphism(5)
            sage: (f.inverse() * f).is_identity()
            True
        """
        return self.__class__(self.domain(), -self.power())

    def _composition(self, right):
        """
        Return ``self`` o ``right``.

        EXAMPLES::

            sage: k.<t> = GF(5^12)
            sage: f = k.frobenius_endomorphism(); f
            Frobenius endomorphism t |--> t^5 on Finite Field in t of size 5^12
            sage: g = k.frobenius_endomorphism(2); g
            Frobenius endomorphism t |--> t^(5^2) on Finite Field in t of size 5^12
            sage: f * g
            Frobenius endomorphism t |--> t^(5^3) on Finite Field in t of size 5^12

        The result is simplified if possible::

            sage: f = k.frobenius_endomorphism(9)
            sage: g = k.frobenius_endomorphism(10)
            sage: f * g
            Frobenius endomorphism t |--> t^(5^7) on Finite Field in t of size 5^12
        """
        if isinstance(right, FrobeniusEndomorphism_finite_field):
            return self.__class__(self.domain(), self._power + right.power())
        else:
            return RingHomomorphism._composition(self, right)

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
            sage: embed(tfixed) == embed.im_gens()[0]
            True
        """
        if self._degree_fixed == 1:
            k = FiniteField(self.domain().characteristic())
            from sage.rings.finite_rings.hom_prime_finite_field import FiniteFieldHomomorphism_prime
            f = FiniteFieldHomomorphism_prime(Hom(k, self.domain()))
        else:
            k = FiniteField(self.domain().characteristic()**self._degree_fixed,
                            name=self.domain().variable_name() + "_fixed")
            f = FiniteFieldHomomorphism_generic(Hom(k, self.domain()))
        return k, f

    def is_injective(self):
        """
        Return ``True`` since any power of the Frobenius endomorphism
        over a finite field is always injective.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: Frob.is_injective()
            True
        """
        return True

    def is_surjective(self):
        """
        Return ``True`` since any power of the Frobenius endomorphism
        over a finite field is always surjective.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: Frob.is_surjective()
            True
        """
        return True

    def is_identity(self):
        """
        Return ``True`` if this morphism is the identity morphism.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: Frob.is_identity()
            False
            sage: (Frob^3).is_identity()
            True
        """
        return self.power() == 0

    def __hash__(self):
        r"""
        Return a hash of this morphism.

        EXAMPLES::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism()
            sage: hash(Frob)  # random
            383183030479672104
        """
        return Morphism.__hash__(self)

    cdef _update_slots(self, dict slots):
        r"""
        Helper function for copying and pickling.

        TESTS::

            sage: k.<t> = GF(5^3)
            sage: Frob = k.frobenius_endomorphism(2)
            sage: Frob
            Frobenius endomorphism t |--> t^(5^2) on Finite Field in t of size 5^3

            sage: phi = copy(Frob)
            sage: phi
            Frobenius endomorphism t |--> t^(5^2) on Finite Field in t of size 5^3
            sage: Frob == phi
            True
        """
        FrobeniusEndomorphism_generic._update_slots(self, slots)
        domain = self.domain()
        self._degree = domain.degree()
        self._degree_fixed = domain.degree().gcd(self._power)
        self._order = self._degree / self._degree_fixed
        self._q = domain.characteristic() ** self._power

from sage.misc.persist import register_unpickle_override
register_unpickle_override('sage.rings.finite_field_morphism', 'FiniteFieldHomomorphism_generic', FiniteFieldHomomorphism_generic)
