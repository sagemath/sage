# sage.doctest: needs sage.rings.number_field
r"""
Embeddings into ambient fields

This module provides classes to handle embeddings of number fields into ambient
fields (generally `\RR` or `\CC`).
"""
# ****************************************************************************
#      Copyright (C) 2008 Robert Bradshaw <robertwb@math.washington.edu>
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
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import sage.rings.complex_double

from sage.structure.element cimport Element
from sage.categories.morphism cimport Morphism
from sage.categories.map cimport Map
from sage.categories.pushout import pushout

from sage.rings.real_lazy import RLF, CLF, LazyField, LazyAlgebraic


cdef class NumberFieldEmbedding(Morphism):

    cdef _gen_image

    def __init__(self, K, R, gen_embedding):
        """
        If R is a lazy field, the closest root to gen_embedding will be chosen.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: from sage.rings.number_field.number_field_morphisms import NumberFieldEmbedding
            sage: K.<a> = NumberField(x^3-2)
            sage: f = NumberFieldEmbedding(K, RLF, 1)
            sage: f(a)^3
            2.00000000000000?
            sage: RealField(200)(f(a)^3)
            2.0000000000000000000000000000000000000000000000000000000000

            sage: sigma_a = K.polynomial().change_ring(CC).roots()[1][0]; sigma_a
            -0.62996052494743... - 1.09112363597172*I
            sage: g = NumberFieldEmbedding(K, CC, sigma_a)
            sage: g(a+1)
            0.37003947505256... - 1.09112363597172*I
        """
        from sage.categories.homset import Hom
        Morphism.__init__(self, Hom(K, R))
        if isinstance(R, LazyField) and not isinstance(gen_embedding.parent(), LazyField):
            self._gen_image = LazyAlgebraic(R, K.polynomial(), gen_embedding, prec=0)
        else:
            self._gen_image = R(gen_embedding)

    cdef dict _extra_slots(self):
        """
        A helper for pickling and copying.

        INPUT:

        - ``_slots`` -- dictionary

        OUTPUT: the given dictionary, with the generator image added

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: from sage.rings.number_field.number_field_morphisms import NumberFieldEmbedding
            sage: K.<a> = NumberField(x^3-2)
            sage: f = NumberFieldEmbedding(K, RLF, 1)
            sage: g = copy(f)    # indirect doctest
            sage: g
            Generic morphism:
              From: Number Field in a with defining polynomial x^3 - 2
              To:   Real Lazy Field
              Defn: a -> 1.259921049894873?
            sage: g(a)^3
            2.00000000000000?
        """
        slots = Morphism._extra_slots(self)
        slots['_gen_image'] = self._gen_image
        return slots

    cdef _update_slots(self, dict _slots):
        """
        A helper for unpickling and copying.

        INPUT:

        - ``_slots`` -- dictionary providing values for the c(p)def slots of ``self``

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: from sage.rings.number_field.number_field_morphisms import NumberFieldEmbedding
            sage: K.<a> = NumberField(x^3-2)
            sage: f = NumberFieldEmbedding(K, RLF, 1)
            sage: g = copy(f)    # indirect doctest
            sage: g
            Generic morphism:
              From: Number Field in a with defining polynomial x^3 - 2
              To:   Real Lazy Field
              Defn: a -> 1.259921049894873?
            sage: g(a)^3
            2.00000000000000?
        """
        Morphism._update_slots(self, _slots)
        self._gen_image = _slots['_gen_image']

    cpdef Element _call_(self, x):
        """
        EXAMPLES::

            sage: x = polygen(QQ)
            sage: from sage.rings.number_field.number_field_morphisms import NumberFieldEmbedding
            sage: K.<a> = NumberField(x^2-2)
            sage: f = NumberFieldEmbedding(K, RLF, 1.4)
            sage: f(a) # indirect doctest
            1.414213562373095?
        """
        return x.polynomial()(self._gen_image)

    def _repr_defn(self):
        """
        EXAMPLES::

            sage: from sage.rings.number_field.number_field_morphisms import NumberFieldEmbedding
            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^2 - 2)
            sage: f = NumberFieldEmbedding(K, RLF, 1.4)
            sage: f # indirect doctest
            Generic morphism:
              From: Number Field in a with defining polynomial x^2 - 2
              To:   Real Lazy Field
              Defn: a -> 1.414213562373095?
        """
        return "{} -> {}".format(self.domain().variable_name(), self._gen_image)

    def gen_image(self):
        """
        Return the image of the generator under this embedding.

        EXAMPLES::

            sage: f = QuadraticField(7, 'a', embedding=2).coerce_embedding()
            sage: f.gen_image()
            2.645751311064591?
        """
        return self._gen_image


cdef class EmbeddedNumberFieldMorphism(NumberFieldEmbedding):
    r"""
    This allows one to go from one number field in another consistently,
    assuming they both have specified embeddings into an ambient field.

    If no ambient field is supplied, then the following ambient fields are
    tried:

    * the pushout of the fields where the number fields are embedded;

    * the algebraic closure of the previous pushout;

    * `\CC`.

    EXAMPLES::

        sage: x = polygen(ZZ, 'x')
        sage: K.<i> = NumberField(x^2 + 1, embedding=QQbar(I))
        sage: L.<i> = NumberField(x^2 + 1, embedding=-QQbar(I))
        sage: from sage.rings.number_field.number_field_morphisms import EmbeddedNumberFieldMorphism
        sage: EmbeddedNumberFieldMorphism(K, L, CDF)
        Generic morphism:
          From: Number Field in i with defining polynomial x^2 + 1 with i = I
          To:   Number Field in i with defining polynomial x^2 + 1 with i = -I
          Defn: i -> -i
        sage: EmbeddedNumberFieldMorphism(K, L, QQbar)
        Generic morphism:
          From: Number Field in i with defining polynomial x^2 + 1 with i = I
          To:   Number Field in i with defining polynomial x^2 + 1 with i = -I
          Defn: i -> -i
    """
    cdef readonly ambient_field

    def __init__(self, K, L, ambient_field=None):
        """
        EXAMPLES::

            sage: from sage.rings.number_field.number_field_morphisms import EmbeddedNumberFieldMorphism
            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^2 - 17, embedding=4.1)
            sage: L.<b> = NumberField(x^4 - 17, embedding=2.0)
            sage: f = EmbeddedNumberFieldMorphism(K, L)
            sage: f(a)
            b^2

            sage: K.<zeta12> = CyclotomicField(12)
            sage: L.<zeta36> = CyclotomicField(36)
            sage: f = EmbeddedNumberFieldMorphism(K, L)
            sage: f(zeta12)
            zeta36^3
            sage: f(zeta12^5-zeta12+1)
            zeta36^9 - 2*zeta36^3 + 1
            sage: f
            Generic morphism:
              From: Cyclotomic Field of order 12 and degree 4
              To:   Cyclotomic Field of order 36 and degree 12
              Defn: zeta12 -> zeta36^3

        The embeddings must be compatible::

            sage: F1 = NumberField(x^3 + 2, 'a', embedding=2)
            sage: F2 = NumberField(x^3 + 2, 'a', embedding=CC.0)
            sage: F1.gen() + F2.gen()
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for +:
            'Number Field in a with defining polynomial x^3 + 2 with a = -1.259921049894873?' and
            'Number Field in a with defining polynomial x^3 + 2 with a = 0.6299605249474365? + 1.091123635971722?*I'

        The following was fixed to raise a :exc:`TypeError` in :issue:`15331`::

            sage: L.<i> = NumberField(x^2 + 1)
            sage: K = NumberField(L(i/2+3).minpoly(), names=('i0',), embedding=L(i/2+3))
            sage: EmbeddedNumberFieldMorphism(K, L)
            Traceback (most recent call last):
            ...
            TypeError: No embedding available for Number Field in i with defining polynomial x^2 + 1
        """
        if ambient_field is None:
            if K.coerce_embedding() is None:
                raise TypeError("No embedding available for %s"%K)
            Kemb = K
            while Kemb.coerce_embedding() is not None:
                Kemb = Kemb.coerce_embedding().codomain()
            if L.coerce_embedding() is None:
                raise TypeError("No embedding available for %s"%L)
            Lemb = L
            while Lemb.coerce_embedding() is not None:
                Lemb = Lemb.coerce_embedding().codomain()
            ambient_field = pushout(Kemb, Lemb)
            candidate_ambient_fields = [ambient_field]
            try:
                candidate_ambient_fields.append(ambient_field.algebraic_closure())
            except NotImplementedError:
                pass
            candidate_ambient_fields.append(sage.rings.complex_double.CDF)
        else:
            candidate_ambient_fields = [ambient_field]

        for ambient_field in candidate_ambient_fields:
            gen_image = matching_root(K.polynomial().change_ring(L), K.gen(), ambient_field=ambient_field, margin=2)
            if gen_image is not None:
                NumberFieldEmbedding.__init__(self, K, L, gen_image)
                self.ambient_field = ambient_field
                return
        else:
            raise ValueError("No consistent embedding of all of %s into %s." % (K, L))

    def section(self):
        """
        EXAMPLES::

            sage: from sage.rings.number_field.number_field_morphisms import EmbeddedNumberFieldMorphism
            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^2 - 700, embedding=25)
            sage: L.<b> = NumberField(x^6 - 700, embedding=3)
            sage: f = EmbeddedNumberFieldMorphism(K, L)
            sage: f(2*a - 1)
            2*b^3 - 1
            sage: g = f.section()
            sage: g(2*b^3 - 1)
            2*a - 1
        """
        return EmbeddedNumberFieldConversion(self.codomain(), self.domain(), self.ambient_field)


cdef class EmbeddedNumberFieldConversion(Map):
    r"""
    This allows one to cast one number field in another consistently,
    assuming they both have specified embeddings into an ambient field
    (by default it looks for an embedding into `\CC`).

    This is done by factoring the minimal polynomial of the input
    in the number field of the codomain. This may fail if the element is
    not actually in the given field.
    """
    cdef _gen_image
    cdef readonly ambient_field

    def __init__(self, K, L, ambient_field=None):
        """
        EXAMPLES::

            sage: from sage.rings.number_field.number_field_morphisms import EmbeddedNumberFieldConversion
            sage: x = polygen(ZZ, 'x')
            sage: K.<a> = NumberField(x^2 - 17, embedding=4.1)
            sage: L.<b> = NumberField(x^4 - 17, embedding=2.0)
            sage: f = EmbeddedNumberFieldConversion(K, L)
            sage: f(a)
            b^2
            sage: f(K(b^2/2-11))
            1/2*b^2 - 11
        """
        if ambient_field is None:
            from sage.rings.complex_double import CDF
            ambient_field = CDF
        self.ambient_field = ambient_field
        Map.__init__(self, K, L)

    cpdef Element _call_(self, x):
        """
        EXAMPLES::

            sage: from sage.rings.number_field.number_field_morphisms import EmbeddedNumberFieldConversion
            sage: K.<zeta12> = CyclotomicField(12)
            sage: L.<zeta15> = CyclotomicField(15)
            sage: f = EmbeddedNumberFieldConversion(K, L)
            sage: f(zeta12^4) # indirect doctest
            zeta15^5
            sage: f(zeta12)
            Traceback (most recent call last):
            ...
            ValueError: No consistent embedding of Cyclotomic Field of order 12 and degree 4 into Cyclotomic Field of order 15 and degree 8.
        """
        minpoly = x.minpoly()
        gen_image = matching_root(minpoly.change_ring(self.codomain()), x, self.ambient_field, 4)
        if gen_image is None:
            raise ValueError("No consistent embedding of {} into {}.".format(self.domain(), self.codomain()))
        return gen_image


cpdef matching_root(poly, target, ambient_field=None, margin=1, max_prec=None):
    """
    Given a polynomial and a ``target``, choose the root that
    ``target`` best approximates as compared in ``ambient_field``.

    If the parent of ``target`` is exact, the equality is required, otherwise
    find closest root (with respect to the ``abs`` function) in the
    ambient field to the ``target``, and return the root of ``poly`` (if any) that
    approximates it best.

    EXAMPLES::

        sage: from sage.rings.number_field.number_field_morphisms import matching_root
        sage: R.<x> = CC[]
        sage: matching_root(x^2-2, 1.5)
        1.41421356237310
        sage: matching_root(x^2-2, -100.0)
        -1.41421356237310
        sage: matching_root(x^2-2, .00000001)
        1.41421356237310
        sage: matching_root(x^3-1, CDF.0)
        -0.50000000000000... + 0.86602540378443...*I
        sage: matching_root(x^3-x, 2, ambient_field=RR)
        1.00000000000000
    """
    if isinstance(poly, list):
        roots = poly
    else:
        roots = poly.roots()
        if len(roots) == 0:
            return None
        elif isinstance(roots[0], tuple): # as returned from the roots method
            roots = [r for r, e in roots]

    if ambient_field is None:
        ambient_field = target.parent()

    if ambient_field.is_exact():
        target_approx = ambient_field(target)
        for r in roots:
            if ambient_field(r) == target_approx:
                return r
    else:
        # since things are inexact, try and pick the closest one
        # -- unless the ambient field is inexact and has no prec(),
        # which holds, e.g., for the symbolic ring
        if not hasattr(ambient_field,'prec'):
            return None
        if max_prec is None:
            max_prec = ambient_field.prec() * 32
        while ambient_field.prec() < max_prec:
            if isinstance(poly, list):
                ambient_roots = [ambient_field(r) for r in poly]
            else:
                ambient_roots = [r for r, e in poly.change_ring(ambient_field).roots()]
            target_root = closest(ambient_field(target), ambient_roots, margin)
            if target_root is not None:
                for r in roots:
                    if closest(ambient_field(r), ambient_roots, margin) is target_root:
                        return r
            ambient_field = ambient_field.to_prec(ambient_field.prec() * 2)


cpdef closest(target, values, margin=1):
    """
    This is a utility function that returns the item in ``values`` closest to
    target (with respect to the ``abs`` function). If ``margin`` is greater
    than 1, and `x` and `y` are the first and second closest elements to ``target``,
    then only return `x` if `x` is ``margin`` times closer to ``target`` than `y`, i.e.
    ``margin * abs(target-x) < abs(target-y)``.

    TESTS::

        sage: from sage.rings.number_field.number_field_morphisms import closest
        sage: closest(1.2, [0,1,2,3,4])
        1
        sage: closest(1.7, [0,1,2,3,4])
        2
        sage: closest(1.7, [0,1,2,3,4], margin=5)
        sage: closest(1.9, [0,1,2,3,4], margin=5)
        2
        sage: closest(.2, [-1, 1, CDF.0, -CDF.0])
        1
    """
    cdef int i
    if len(values) == 0:
        raise ValueError
    elif len(values) == 1:
        return values[0]
    else:
        dists = [abs(target - r) for r in values]
        sdists = sorted(dists)
        min_dist = sdists[0]
        if margin*min_dist < sdists[1]:
            for i in range(len(values)):
                if dists[i] is min_dist:
                    return values[i]
        else:
            return None


def root_from_approx(f, a):
    """
    Return an exact root of the polynomial `f` closest to `a`.

    INPUT:

    - ``f`` -- polynomial with rational coefficients

    - ``a`` -- element of a ring

    OUTPUT:

    A root of ``f`` in the parent of ``a`` or, if ``a`` is not already
    an exact root of ``f``, in the corresponding lazy field.  The root
    is taken to be closest to ``a`` among all roots of ``f``.

    EXAMPLES::

        sage: from sage.rings.number_field.number_field_morphisms import root_from_approx
        sage: R.<x> = QQ[]

        sage: root_from_approx(x^2 - 1, -1)
        -1
        sage: root_from_approx(x^2 - 2, 1)
        1.414213562373095?
        sage: root_from_approx(x^3 - x - 1, RR(1))
        1.324717957244746?
        sage: root_from_approx(x^3 - x - 1, CC.gen())
        -0.6623589786223730? + 0.5622795120623013?*I

        sage: root_from_approx(x^2 + 1, 0)
        Traceback (most recent call last):
        ...
        ValueError: x^2 + 1 has no real roots
        sage: root_from_approx(x^2 + 1, CC(0))
        -1*I

        sage: root_from_approx(x^2 - 2, sqrt(2))                                        # needs sage.symbolic
        sqrt(2)
        sage: root_from_approx(x^2 - 2, sqrt(3))                                        # needs sage.symbolic
        Traceback (most recent call last):
        ...
        ValueError: sqrt(3) is not a root of x^2 - 2
    """
    P = a.parent()
    if P.is_exact() and not f(a):
        return a
    elif P._is_real_numerical():
        return LazyAlgebraic(RLF, f, a, prec=0)
    elif P._is_numerical():
        return LazyAlgebraic(CLF, f, a, prec=0)
    # p-adic lazy, when implemented, would go here
    else:
        from sage.symbolic.relation import test_relation_maxima
        rel = (f(a) != 0)
        if (rel is True
            or (not isinstance(rel, bool) and test_relation_maxima(rel))):
            raise ValueError("{} is not a root of {}".format(a, f))
        return a


def create_embedding_from_approx(K, gen_image):
    """
    Return an embedding of ``K`` determined by ``gen_image``.

    The codomain of the embedding is the parent of ``gen_image`` or,
    if ``gen_image`` is not already an exact root of the defining
    polynomial of ``K``, the corresponding lazy field.  The embedding
    maps the generator of ``K`` to a root of the defining polynomial
    of ``K`` closest to ``gen_image``.

    EXAMPLES::

        sage: from sage.rings.number_field.number_field_morphisms import create_embedding_from_approx
        sage: x = polygen(ZZ, 'x')
        sage: K.<a> = NumberField(x^3 - x + 1/10)
        sage: create_embedding_from_approx(K, 1)
        Generic morphism:
          From: Number Field in a with defining polynomial x^3 - x + 1/10
          To:   Real Lazy Field
          Defn: a -> 0.9456492739235915?
        sage: create_embedding_from_approx(K, 0)
        Generic morphism:
          From: Number Field in a with defining polynomial x^3 - x + 1/10
          To:   Real Lazy Field
          Defn: a -> 0.10103125788101081?
        sage: create_embedding_from_approx(K, -1)
        Generic morphism:
          From: Number Field in a with defining polynomial x^3 - x + 1/10
          To:   Real Lazy Field
          Defn: a -> -1.046680531804603?

    We can define embeddings from one number field to another::

        sage: L.<b> = NumberField(x^6-x^2+1/10)
        sage: create_embedding_from_approx(K, b^2)
        Generic morphism:
          From: Number Field in a with defining polynomial x^3 - x + 1/10
          To:   Number Field in b with defining polynomial x^6 - x^2 + 1/10
          Defn: a -> b^2

    If the embedding is exact, it must be valid::

        sage: create_embedding_from_approx(K, b)
        Traceback (most recent call last):
        ...
        ValueError: b is not a root of x^3 - x + 1/10
    """
    if gen_image is None:
        return None
    elif isinstance(gen_image, Map):
        return gen_image
    elif isinstance(gen_image, Element):
        x = root_from_approx(K.defining_polynomial(), gen_image)
        return NumberFieldEmbedding(K, x.parent(), x)
    else:
        raise TypeError("Embedding (type %s) must be a morphism or element." % type(gen_image))


cdef class CyclotomicFieldEmbedding(NumberFieldEmbedding):
    """
    Specialized class for converting cyclotomic field elements into a
    cyclotomic field of higher order. All the real work is done by
    :meth:`_lift_cyclotomic_element`.
    """

    cdef ratio

    def __init__(self, K, L):
        """
        Check and cache the parameters.

        EXAMPLES::

            sage: from sage.rings.number_field.number_field_morphisms import CyclotomicFieldEmbedding
            sage: CyclotomicFieldEmbedding(CyclotomicField(7), CyclotomicField(21))
            Generic morphism:
              From: Cyclotomic Field of order 7 and degree 6
              To:   Cyclotomic Field of order 21 and degree 12
              Defn: zeta7 -> zeta21^3

        Note that this only handles the easy case of cyclotomic fields where
        the order of the smaller dividing the order of the larger, regardless
        of whether or not there is an actual coercion::

            sage: CyclotomicFieldEmbedding(CyclotomicField(3), QuadraticField(-3, 'a'))
            Traceback (most recent call last):
            ...
            TypeError: CyclotomicFieldEmbedding only valid for cyclotomic fields.
            sage: CyclotomicFieldEmbedding(CyclotomicField(14), CyclotomicField(21))
            Traceback (most recent call last):
            ...
            TypeError: The zeta_order of the new field must be a multiple of the zeta_order of the original.

        Check that :issue:`13765` is fixed::

            sage: z3=(CC(-1)^(1/3))^2
            sage: Ka.<a>=CyclotomicField(3,embedding=z3)
            sage: Kb.<b>=CyclotomicField(3,embedding=z3^2)
            sage: CyclotomicFieldEmbedding(Ka, Kb)
            Generic morphism:
              From: Cyclotomic Field of order 3 and degree 2
              To:   Cyclotomic Field of order 3 and degree 2
              Defn: a -> -b - 1
            sage: Ka(b)
            -a - 1
            sage: a + b
            -1
            sage: b + a
            -1
        """
        Morphism.__init__(self, K, L)
        from sage.rings.number_field.number_field import NumberField_cyclotomic
        if not isinstance(K, NumberField_cyclotomic) or not isinstance(L, NumberField_cyclotomic):
            raise TypeError("CyclotomicFieldEmbedding only valid for cyclotomic fields.")
        Kn = K._n()
        Ln = L._n()
        if not Kn.divides(Ln):
            raise TypeError("The zeta_order of the new field must be a multiple of the zeta_order of the original.")
        self.ratio = L._log_gen(K.coerce_embedding()(K.gen()))
        self._gen_image = L.gen() ** self.ratio

    cdef dict _extra_slots(self):
        """
        A helper for pickling and copying.

        INPUT:

        - ``_slots`` -- dictionary

        OUTPUT: the given dictionary, with _gen_image and ratio added

        EXAMPLES::

            sage: from sage.rings.number_field.number_field_morphisms import CyclotomicFieldEmbedding
            sage: cf6 = CyclotomicField(6)
            sage: cf12 = CyclotomicField(12)
            sage: f = CyclotomicFieldEmbedding(cf6, cf12)
            sage: g = copy(f) # indirect doctest
            sage: g
            Generic morphism:
              From: Cyclotomic Field of order 6 and degree 2
              To:   Cyclotomic Field of order 12 and degree 4
              Defn: zeta6 -> zeta12^2
            sage: g(cf6.0)
            zeta12^2
        """
        slots = NumberFieldEmbedding._extra_slots(self)
        slots['ratio'] = self.ratio
        return slots

    cdef _update_slots(self, dict _slots):
        """
        A helper for unpickling and copying.

        INPUT:

        - ``_slots`` -- dictionary providing values for the c(p)def slots of ``self``

        EXAMPLES::

            sage: from sage.rings.number_field.number_field_morphisms import CyclotomicFieldEmbedding
            sage: cf6 = CyclotomicField(6)
            sage: cf12 = CyclotomicField(12)
            sage: f = CyclotomicFieldEmbedding(cf6, cf12)
            sage: g = copy(f) # indirect doctest
            sage: g
            Generic morphism:
              From: Cyclotomic Field of order 6 and degree 2
              To:   Cyclotomic Field of order 12 and degree 4
              Defn: zeta6 -> zeta12^2
            sage: g(cf6.0)
            zeta12^2
        """
        Morphism._update_slots(self, _slots)
        self._gen_image = _slots['_gen_image']
        self.ratio = _slots['ratio']

    cpdef Element _call_(self, x):
        """
        EXAMPLES::

            sage: from sage.rings.number_field.number_field_morphisms import CyclotomicFieldEmbedding
            sage: K = CyclotomicField(7)
            sage: L = CyclotomicField(21)
            sage: f = CyclotomicFieldEmbedding(K, L)
            sage: f(K.gen()) # indirect doctest
            zeta21^3
            sage: f(K.gen()^2 + 3) # indirect doctest
            zeta21^6 + 3
        """
        return x._lift_cyclotomic_element(self.codomain(), False, self.ratio)

    def section(self):
        """
        Return the section of ``self``.

        EXAMPLES::

            sage: from sage.rings.number_field.number_field_morphisms import CyclotomicFieldEmbedding
            sage: K = CyclotomicField(7)
            sage: L = CyclotomicField(21)
            sage: f = CyclotomicFieldEmbedding(K, L)
            sage: h = f.section()
            sage: h(f(K.gen())) # indirect doctest
            zeta7
        """
        return CyclotomicFieldConversion(self.codomain(), self.domain())

cdef class CyclotomicFieldConversion(Map):
    r"""
    This allows one to cast one cyclotomic field in another consistently.

    EXAMPLES::

        sage: from sage.rings.number_field.number_field_morphisms import CyclotomicFieldConversion
        sage: K1.<z1> = CyclotomicField(12)
        sage: K2.<z2> = CyclotomicField(18)
        sage: f = CyclotomicFieldConversion(K1, K2)
        sage: f(z1^2)
        z2^3
        sage: f(z1)
        Traceback (most recent call last):
        ...
        ValueError: Element z1 has no image in the codomain

    Tests from :issue:`29511`::

        sage: K.<z> = CyclotomicField(12)
        sage: K1.<z1> = CyclotomicField(3)
        sage: K(2) in K1 # indirect doctest
        True
        sage: K1(K(2)) # indirect doctest
        2
    """
    cdef ambient_field
    cdef phi

    def __init__(self, K, L):
        """
        Construct a conversion map between cyclotomic fields.

        EXAMPLES::

            sage: from sage.rings.number_field.number_field_morphisms import CyclotomicFieldEmbedding
            sage: K.<a> = CyclotomicField(7)
            sage: L.<b> = CyclotomicField(21)
            sage: K(b^3) # indirect doctest
            a
        """
        from sage.rings.number_field.number_field import CyclotomicField
        n1 = K._n()
        n2 = L._n()
        n3 = n1.lcm(n2)
        M = CyclotomicField(n3)
        self.ambient_field = M
        self.phi = L.hom([M.gen()**(n3//n2)])
        Map.__init__(self, K, L)

    cpdef Element _call_(self, x):
        """
        Call a conversion map between cyclotomic fields.

        EXAMPLES::

            sage: from sage.rings.number_field.number_field_morphisms import CyclotomicFieldEmbedding
            sage: K.<a> = CyclotomicField(7)
            sage: L.<b> = CyclotomicField(21)
            sage: K(b^3) # indirect doctest
            a
        """
        M = self.ambient_field
        try:
            return self.phi.preimage(M(x))
        except ValueError:
            raise ValueError('Element {} has no image in the codomain'.format(x))
