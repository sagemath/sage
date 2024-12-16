r"""
Morphisms between extension of rings

AUTHOR:

- Xavier Caruso (2019)
"""

#############################################################################
#    Copyright (C) 2019 Xavier Caruso <xavier.caruso@normalesup.org>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  http://www.gnu.org/licenses/
#****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.structure.richcmp import op_EQ, op_NE

from sage.structure.element cimport Element
from sage.categories.map import Map
from sage.structure.parent cimport Parent
from sage.rings.morphism cimport RingMap
from sage.rings.ring_extension cimport RingExtension_generic, RingExtensionWithBasis
from sage.rings.ring_extension_element cimport RingExtensionElement
from sage.rings.ring_extension_conversion cimport backend_parent, backend_element, backend_morphism


# I don't trust the operator ==
cdef are_equal_morphisms(f, g):
    r"""
    Return ``True`` if ``f`` and ``g`` coincide on the
    generators of the domain, ``False`` otherwise.

    INPUT:

    - ``f`` -- a ring homomorphism or ``None``; if ``None``,
      we consider that ``f`` is a coercion map

    - ``g`` -- a ring homomorphism or ``None``

    TESTS::

        sage: S.<x> = QQ[]
        sage: T.<y> = S[]
        sage: TT = T.over(QQ)
        sage: H = End(TT)

        sage: cc = S.hom([-x])
        sage: f = T.hom([x^2 + y^2], base_map=cc)
        sage: g = T.hom([x^2 + y^2])

        sage: H(f) == H(g)      # indirect doctest
        False
        sage: H(f^2) == H(g^2)  # indirect doctest
        True
    """
    cdef Parent b
    cdef tuple gens
    if f is None and g is None:
        return True
    if f is None:
        f, g = g, f
    gens = tuple()
    b = f.domain()
    while b is not b._base:
        gens += b.gens()
        b = b._base
    if g is None:
        for x in gens:
            if f(x) != x: return False
    else:
        for x in gens:
            if f(x) != g(x): return False
    return True


cdef class RingExtensionHomomorphism(RingMap):
    r"""
    A class for ring homomorphisms between extensions.

    TESTS::

        sage: # needs sage.rings.finite_rings
        sage: K.<a> = GF(5^2).over()
        sage: L.<b> = GF(5^4).over(K)
        sage: phi = L.hom([b^5, a^5]); phi
        Ring endomorphism of Field in b
         with defining polynomial x^2 + (3 - a)*x + a over its base
          Defn: b |--> (2 + a) + 2*b
                with map on base ring:
                a |--> 1 - a
        sage: type(phi)
        <class 'sage.rings.ring_extension_morphism.RingExtensionHomomorphism'>
        sage: TestSuite(phi).run()
    """
    def __init__(self, parent, defn, base_map=None, check=True):
        r"""
        Initialize this morphism.

        INPUT:

        - ``defn`` -- the definition of the morphism (either a map or images of generators)

        - ``base_map`` -- a ring homomorphism or ``None`` (default: ``None``);
          the action of this morphism on one of the bases of the domain;
          if ``None``, a coercion map is used

        - ``check`` -- boolean (default: ``True``); whether to check if
          the given data define a valid homomorphism

        TESTS::

            sage: S.<x> = QQ[]
            sage: T.<x,y> = QQ[]
            sage: f = T.hom([x^2, y^2]); f
            Ring endomorphism of Multivariate Polynomial Ring in x, y over Rational Field
              Defn: x |--> x^2
                    y |--> y^2

            sage: TT = T.over(QQ)
            sage: End(TT)(f)
            Ring endomorphism of Multivariate Polynomial Ring in x, y over Rational Field over its base
              Defn: x |--> x^2
                    y |--> y^2

            sage: TT = T.over(S)
            sage: End(TT)(f)
            Ring endomorphism of Multivariate Polynomial Ring in x, y over Rational Field over its base
              Defn: y |--> y^2
                    with map on base ring:
                    x |--> x^2
        """
        RingMap.__init__(self, parent)
        domain = self.domain()
        backend_domain = backend_parent(domain)
        codomain = self.codomain()
        backend_codomain = backend_parent(codomain)
        # We construct the backend morphism
        if isinstance(defn, Map):
            if base_map is not None:
                raise ValueError("base_map cannot be set when passing in the backend morphism")
            backend = backend_morphism(defn)
            if backend.domain() is not backend_domain:
                raise TypeError("the domain of the backend morphism is not correct")
            if backend.codomain() is not backend_codomain:
                raise TypeError("the codomain of the backend morphism is not correct")
            self._backend = backend
            self._im_gens = None
            self._base_map_construction = False
        elif isinstance(defn, (list, tuple)):
            # We figure out what is the base
            if base_map is not None:
                base = base_map.domain()
                gens = domain.gens(base)
            else:
                base = domain
                gens = tuple()
                while True:
                    if len(gens) == len(defn):
                        break
                    if len(gens) > len(defn) or base is base.base_ring():
                        raise ValueError("the number of images does not match the number of generators")
                    gens += base.gens()
                    base = base.base_ring()
            # We construct the backend morphism
            im_gens = [ codomain(x) for x in defn ]
            backend_bases = [ backend_domain ]
            b = backend_domain.base_ring()
            while b is not b.base_ring():
                backend_bases.append(b)
                b = b.base_ring()
            backend_bases.reverse()
            current_morphism = None
            for current_domain in backend_bases:
                current_im_gens = [ ]
                for x in current_domain.gens():
                    pol = domain(backend_domain(x)).polynomial(base)
                    if base_map is not None:
                        pol = pol.map_coefficients(base_map)
                    y = pol(im_gens)
                    current_im_gens.append(backend_element(y))
                current_morphism = current_domain.hom(current_im_gens, base_map=current_morphism, check=check)
            # We check that everything went well
            if check:
                for i in range(len(gens)):
                    x = backend_element(domain(gens[i]))
                    y = backend_element(im_gens[i])
                    if current_morphism(x) != y:
                        raise ValueError("images do not define a valid homomorphism")
                coercion_morphism = backend_morphism(domain.defining_morphism(base))
                if base_map is None:
                    backend_base_map = coercion_morphism
                else:
                    backend_base_map = backend_morphism(base_map)
                restriction_current_morphism = current_morphism * coercion_morphism
                if not are_equal_morphisms(restriction_current_morphism, backend_base_map):
                    raise ValueError("images do not define a valid homomorphism")
            self._backend = current_morphism
            self._im_gens = im_gens[:domain.ngens()]
            if base is domain.base_ring():
                self._base_map_construction = base_map
            else:
                self._base_map_construction = {
                    'im_gens': defn[domain.ngens():],
                    'base_map': base_map,
                    'check': False
                }
        else:
            raise TypeError

    def _repr_type(self):
        r"""
        Return a string that describes the type of this morphism.

        EXAMPLES::

            sage: K.<a> = GF(5^2).over()  # over GF(5)                                  # needs sage.rings.finite_rings
            sage: f = K.hom([a^5]); f                                                   # needs sage.rings.finite_rings
            Ring endomorphism of Field in a with defining polynomial x^2 + 4*x + 2 over its base
              Defn: a |--> 1 - a

            sage: f._repr_type()                                                        # needs sage.rings.finite_rings
            'Ring'
        """
        return "Ring"

    cpdef Element _call_(self, x):
        r"""
        Return the image of ``x`` under this morphism.

        INPUT:

        - ``x`` -- an element in the domain of this morphism

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: x = polygen(QQ, 'x')
            sage: A.<sqrt2> = QQ.extension(x^2 - 2)
            sage: K.<sqrt2> = A.over()
            sage: f = K.hom([-sqrt2]); f
            Ring endomorphism of Field in sqrt2 with defining polynomial x^2 - 2 over its base
              Defn: sqrt2 |--> -sqrt2
            sage: f(sqrt2)
            -sqrt2

        TESTS::

            sage: a = QQ.random_element()
            sage: b = QQ.random_element()
            sage: f(a + b*sqrt2) == a - b*sqrt2                                         # needs sage.rings.number_field
            True
        """
        y = self._backend(backend_element(x))
        if isinstance(self.codomain(), RingExtension_generic):
            y = self._codomain(y)
        return y

    @cached_method
    def base_map(self):
        r"""
        Return the base map of this morphism
        or just ``None`` if the base map is a coercion map.

        EXAMPLES::

            sage: F = GF(5)
            sage: K.<a> = GF(5^2).over(F)                                               # needs sage.rings.finite_rings
            sage: L.<b> = GF(5^6).over(K)                                               # needs sage.rings.finite_rings

        We define the absolute Frobenius of L::

            sage: FrobL = L.hom([b^5, a^5]); FrobL                                      # needs sage.rings.finite_rings
            Ring endomorphism of
             Field in b with defining polynomial x^3 + (2 + 2*a)*x - a over its base
              Defn: b |--> (-1 + a) + (1 + 2*a)*b + a*b^2
                    with map on base ring:
                    a |--> 1 - a
            sage: FrobL.base_map()                                                      # needs sage.rings.finite_rings
            Ring morphism:
              From: Field in a with defining polynomial x^2 + 4*x + 2 over its base
              To:   Field in b with defining polynomial x^3 + (2 + 2*a)*x - a over its base
              Defn: a |--> 1 - a

        The square of ``FrobL`` acts trivially on K; in other words, it has
        a trivial base map::

            sage: phi = FrobL^2; phi                                                    # needs sage.rings.finite_rings
            Ring endomorphism of
             Field in b with defining polynomial x^3 + (2 + 2*a)*x - a over its base
              Defn: b |--> 2 + 2*a*b + (2 - a)*b^2
            sage: phi.base_map()                                                        # needs sage.rings.finite_rings
        """
        domain = self.domain()
        codomain = self.codomain()
        base = domain.base_ring()
        if base is base.base_ring():
            return None
        base_map = self._base_map_construction
        if base_map is False:
            if base is domain:
                base_map = None
            else:
                base_map = self * domain.coerce_map_from(base)
        elif isinstance(base_map, dict):
            base_map = base.hom(**self._base_map_construction)
        if base_map is None:
            return None
        if (codomain.has_coerce_map_from(base) and
            are_equal_morphisms(backend_morphism(base_map),
                                backend_morphism(codomain.coerce_map_from(base)))):
            return None
        if base_map.codomain() is not self.codomain():
            base_map = base_map.extend_codomain(self.codomain())
        return base_map

    cpdef _richcmp_(self, other, int op):
        r"""
        Compare this element with ``other`` according to
        the rich comparison operator ``op``.

        INPUT:

        - ``other`` -- a morphism with the same codomain and codomain

        - ``op`` -- the comparison operator

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: K.<a> = GF(5^2).over()   # over GF(5)
            sage: L.<b> = GF(5^6).over(K)
            sage: FrobK = K.hom([a^5])
            sage: FrobL = L.hom([b^5], base_map=FrobK)
            sage: FrobK^2 == End(K).identity()
            True
            sage: FrobL^6 == End(L).identity()
            True
        """
        eq = are_equal_morphisms(self._backend, backend_morphism(other))
        if op == op_EQ:
            return eq
        if op == op_NE:
            return not eq
        return NotImplemented

    def is_identity(self):
        r"""
        Return whether this morphism is the identity.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<a> = GF(5^2).over()   # over GF(5)
            sage: FrobK = K.hom([a^5])
            sage: FrobK.is_identity()
            False
            sage: (FrobK^2).is_identity()
            True

        Coercion maps are not considered as identity morphisms::

            sage: # needs sage.rings.finite_rings
            sage: L.<b> = GF(5^6).over(K)
            sage: iota = L.defining_morphism(); iota
            Ring morphism:
              From: Field in a with defining polynomial x^2 + 4*x + 2 over its base
              To:   Field in b with defining polynomial x^3 + (2 + 2*a)*x - a over its base
              Defn: a |--> a
            sage: iota.is_identity()
            False
        """
        if self.domain() is not self.codomain():
            return False
        return are_equal_morphisms(self._backend, None)

    def is_injective(self):
        r"""
        Return whether this morphism is injective.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K = GF(5^10).over(GF(5^5))
            sage: iota = K.defining_morphism(); iota
            Ring morphism:
              From: Finite Field in z5 of size 5^5
              To:   Field in z10 with defining polynomial
                    x^2 + (2*z5^3 + 2*z5^2 + 4*z5 + 4)*x + z5 over its base
              Defn: z5 |--> z5
            sage: iota.is_injective()
            True

            sage: K = GF(7).over(ZZ)
            sage: iota = K.defining_morphism(); iota
            Ring morphism:
              From: Integer Ring
              To:   Finite Field of size 7 over its base
              Defn: 1 |--> 1
            sage: iota.is_injective()
            False
        """
        return self._backend.is_injective()

    def is_surjective(self):
        r"""
        Return whether this morphism is surjective.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K = GF(5^10).over(GF(5^5))
            sage: iota = K.defining_morphism(); iota
            Ring morphism:
              From: Finite Field in z5 of size 5^5
              To:   Field in z10 with defining polynomial
                    x^2 + (2*z5^3 + 2*z5^2 + 4*z5 + 4)*x + z5 over its base
              Defn: z5 |--> z5
            sage: iota.is_surjective()
            False

            sage: K = GF(7).over(ZZ)
            sage: iota = K.defining_morphism(); iota
            Ring morphism:
              From: Integer Ring
              To:   Finite Field of size 7 over its base
              Defn: 1 |--> 1
            sage: iota.is_surjective()
            True
        """
        return self._backend.is_surjective()

    def _repr_defn(self):
        r"""
        Return a string definition of this morphism.

        By default, we show the action of the morphism on the
        generators of the domain.

        EXAMPLES::

            sage: # needs sage.rings.finite_rings
            sage: K.<a> = GF(5^2).over()   # over GF(5)
            sage: L.<b> = GF(5^6).over(K)
            sage: FrobL = L.hom([b^5, a^5])  # absolute Frobenius
            sage: print(FrobL._repr_defn())
            b |--> (-1 + a) + (1 + 2*a)*b + a*b^2
            with map on base ring:
            a |--> 1 - a
        """
        import re
        s = ""
        gens = self.domain().gens()
        if self._im_gens is None:
            self._im_gens = [ self(x) for x in gens ]
        for i in range(len(gens)):
            s += "%s |--> %s\n" % (gens[i], self._im_gens[i])
        if self.base_map() is not None:
            s += "with map on base ring"
            ss = self.base_map()._repr_defn()
            ss = re.sub('\nwith map on base ring:?$', '', ss, count=0, flags=re.MULTILINE)
            if ss != "": s += ":\n" + ss
        if s != "" and s[-1] == "\n":
            s = s[:-1]
        return s

    def _composition(self, right):
        r"""
        Return the composite ``self o right``.

        TESTS::

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: A.<sqrt5> = QQ.extension(x^2 - 5)
            sage: K.<sqrt5> = A.over()
            sage: f = K.hom([-sqrt5]); f
            Ring endomorphism of Field in sqrt5 with defining polynomial x^2 - 5 over its base
              Defn: sqrt5 |--> -sqrt5
            sage: f^2  # indirect doctest
            Ring endomorphism of Field in sqrt5 with defining polynomial x^2 - 5 over its base
              Defn: sqrt5 |--> sqrt5
        """
        domain = right.domain()
        codomain = self.codomain()
        backend_right = backend_morphism(right)
        backend = self._backend * backend_right
        if isinstance(domain, RingExtension_generic) or isinstance(codomain, RingExtension_generic):
            return RingExtensionHomomorphism(domain.Hom(codomain), backend)
        else:
            return backend

    cdef _update_slots(self, dict _slots):
        """
        Helper function for copying and pickling.

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: K.<a> = GF(5^2).over()   # over GF(5)
            sage: f = K.hom([a^5])
            sage: g = copy(f)    # indirect doctest
            sage: f == g
            True
            sage: f is g
            False
        """
        self._backend = _slots['_backend']
        RingMap._update_slots(self, _slots)

    cdef dict _extra_slots(self):
        """
        Helper function for copying and pickling.

        TESTS::

            sage: K.<a> = GF(5^2).over()   # over GF(5)                                 # needs sage.rings.finite_rings
            sage: f = K.hom([a^5])                                                      # needs sage.rings.finite_rings
            sage: loads(dumps(f)) == f                                                  # needs sage.rings.finite_rings
            True
        """
        slots = RingMap._extra_slots(self)
        slots['_backend'] = self._backend
        return slots


cdef class RingExtensionBackendIsomorphism(RingExtensionHomomorphism):
    r"""
    A class for implementating isomorphisms taking an element of the
    backend to its ring extension.

    TESTS::

        sage: # needs sage.rings.finite_rings
        sage: K = GF(11^9).over(GF(11^3))
        sage: f = K.coerce_map_from(GF(11^9)); f
        Coercion morphism:
          From: Finite Field in z9 of size 11^9
          To:   Field in z9 with defining polynomial
                x^3 + (9*z3^2 + 5*z3 + 1)*x^2 + (4*z3 + 3)*x + 10*z3 over its base
        sage: type(f)
        <class 'sage.rings.ring_extension_morphism.RingExtensionBackendIsomorphism'>
        sage: TestSuite(f).run()
    """
    def __init__(self, parent):
        r"""
        Initialize this morphism.

        TESTS::

            sage: x = polygen(ZZ, 'x')
            sage: A.<a> = QQ.extension(x^2 - 5)                                         # needs sage.rings.number_field
            sage: K = A.over()                                                          # needs sage.rings.number_field
            sage: K.coerce_map_from(A)                                                  # needs sage.rings.number_field
            Coercion morphism:
              From: Number Field in a with defining polynomial x^2 - 5
              To:   Field in a with defining polynomial x^2 - 5 over its base
        """
        RingMap.__init__(self, parent)
        domain = self.domain()
        self._backend = domain.Hom(domain).identity()

    def _repr_type(self):
        r"""
        Return a string that describes the type of this morphism.

        EXAMPLES::

            sage: K.<a> = GF(5^2).over()  # over GF(5)                                  # needs sage.rings.finite_rings
            sage: f = K.coerce_map_from(GF(5^2)); f                                     # needs sage.rings.finite_rings
            Coercion morphism:
              From: Finite Field in z2 of size 5^2
              To:   Field in a with defining polynomial x^2 + 4*x + 2 over its base

            sage: f._repr_type()                                                        # needs sage.rings.finite_rings
            'Coercion'
        """
        return "Coercion"

    def _repr_defn(self):
        r"""
        Return the empty string since this morphism is canonical.

        EXAMPLES::

            sage: K.<a> = GF(5^2).over()  # over GF(5)                                  # needs sage.rings.finite_rings
            sage: f = K.coerce_map_from(GF(5^2)); f                                     # needs sage.rings.finite_rings
            Coercion morphism:
              From: Finite Field in z2 of size 5^2
              To:   Field in a with defining polynomial x^2 + 4*x + 2 over its base

            sage: f._repr_defn()                                                        # needs sage.rings.finite_rings
            ''
        """
        return ""

    cpdef Element _call_(self, x):
        r"""
        Return the image of ``x`` under this morphism.

        INPUT:

        - ``x`` -- an element in the domain of this morphism

        EXAMPLES::

            sage: K.<a> = GF(5^2).over()  # over GF(5)                                  # needs sage.rings.finite_rings
            sage: f = K.coerce_map_from(GF(5^2))                                        # needs sage.rings.finite_rings
            sage: f(GF(5^2).gen())                                                      # needs sage.rings.finite_rings
            a
        """
        codomain = self.codomain()
        return codomain.element_class(codomain, x)


cdef class RingExtensionBackendReverseIsomorphism(RingExtensionHomomorphism):
    r"""
    A class for implementating isomorphisms from a ring extension to
    its backend.

    TESTS::

        sage: # needs sage.rings.finite_rings
        sage: K = GF(11^9).over(GF(11^3))
        sage: f = GF(11^9).convert_map_from(K); f
        Canonical morphism:
          From: Field in z9 with defining polynomial
                x^3 + (9*z3^2 + 5*z3 + 1)*x^2 + (4*z3 + 3)*x + 10*z3 over its base
          To:   Finite Field in z9 of size 11^9
        sage: type(f)
        <class 'sage.rings.ring_extension_morphism.RingExtensionBackendReverseIsomorphism'>
        sage: TestSuite(f).run()
    """
    def __init__(self, parent):
        r"""
        Initialize this morphism.

        TESTS::

            sage: x = polygen(ZZ, 'x')
            sage: A.<a> = QQ.extension(x^2 - 5)                                         # needs sage.rings.number_field
            sage: K = A.over()                                                          # needs sage.rings.number_field
            sage: A.convert_map_from(K)                                                 # needs sage.rings.number_field
            Canonical morphism:
              From: Field in a with defining polynomial x^2 - 5 over its base
              To:   Number Field in a with defining polynomial x^2 - 5
        """
        RingMap.__init__(self, parent)
        codomain = self.codomain()
        self._backend = codomain.Hom(codomain).identity()

    def _repr_type(self):
        r"""
        Return a string that describes the type of this morphism.

        EXAMPLES::

            sage: K.<a> = GF(5^2).over()  # over GF(5)                                  # needs sage.rings.finite_rings
            sage: f = GF(5^2).convert_map_from(K); f                                    # needs sage.rings.finite_rings
            Canonical morphism:
              From: Field in a with defining polynomial x^2 + 4*x + 2 over its base
              To:   Finite Field in z2 of size 5^2

            sage: f._repr_type()                                                        # needs sage.rings.finite_rings
            'Canonical'
        """
        return "Canonical"

    def _repr_defn(self):
        r"""
        Return the empty string since this morphism is canonical.

        EXAMPLES::

            sage: K.<a> = GF(5^2).over()  # over GF(5)                                  # needs sage.rings.finite_rings
            sage: f = GF(5^2).convert_map_from(K); f                                    # needs sage.rings.finite_rings
            Canonical morphism:
              From: Field in a with defining polynomial x^2 + 4*x + 2 over its base
              To:   Finite Field in z2 of size 5^2

            sage: f._repr_defn()                                                        # needs sage.rings.finite_rings
            ''
        """
        return ""

    cpdef Element _call_(self, x):
        r"""
        Return the image of ``x`` under this morphism.

        INPUT:

        - ``x`` -- an element in the domain of this morphism

        EXAMPLES::

            sage: K.<a> = GF(5^2).over()  # over GF(5)                                  # needs sage.rings.finite_rings
            sage: f = GF(5^2).convert_map_from(K)                                       # needs sage.rings.finite_rings
            sage: f(a)                                                                  # needs sage.rings.finite_rings
            z2
        """
        return (<RingExtensionElement>x)._backend


cdef class MapFreeModuleToRelativeRing(Map):
    """
    Base class of the module isomorphism between a ring extension
    and a free module over one of its bases.

    TESTS::

        sage: K = GF(5^2).over()                                                        # needs sage.rings.finite_rings
        sage: V, i, j = K.free_module()                                                 # needs sage.rings.finite_rings
        sage: type(i)                                                                   # needs sage.rings.finite_rings
        <class 'sage.rings.ring_extension_morphism.MapFreeModuleToRelativeRing'>
    """
    def __init__(self, E, K):
        r"""
        Initialize this morphism.

        INPUT:

        - ``E`` -- a ring extension

        - ``K`` -- a commutative ring; one base of ``E``

        TESTS::

            sage: K = GF(11^6).over(GF(11^3))                                           # needs sage.rings.finite_rings
            sage: V, i, j = K.free_module()                                             # needs sage.rings.finite_rings
            sage: i                                                                     # needs sage.rings.finite_rings
            Generic map:
              From: Vector space of dimension 2 over Finite Field in z3 of size 11^3
              To:   Field in z6 with defining polynomial x^2 + (10*z3^2 + z3 + 6)*x + z3 over its base
        """
        self._degree = E.degree(K)
        self._basis = [ (<RingExtensionElement>x)._backend for x in E.basis_over(K) ]
        self._f = backend_morphism(E.defining_morphism(K), forget='codomain')
        domain = K ** self._degree
        parent = domain.Hom(E)
        Map.__init__(self, parent)

    def is_injective(self):
        r"""
        Return whether this morphism is injective.

        EXAMPLES::

            sage: K = GF(11^6).over(GF(11^3))                                           # needs sage.rings.finite_rings
            sage: V, i, j = K.free_module()                                             # needs sage.rings.finite_rings
            sage: i.is_injective()                                                      # needs sage.rings.finite_rings
            True
        """
        return True

    def is_surjective(self):
        r"""
        Return whether this morphism is surjective.

        EXAMPLES::

            sage: K = GF(11^6).over(GF(11^3))                                           # needs sage.rings.finite_rings
            sage: V, i, j = K.free_module()                                             # needs sage.rings.finite_rings
            sage: i.is_surjective()                                                     # needs sage.rings.finite_rings
            True
        """
        return True

    cpdef Element _call_(self, v):
        r"""
        Return the image of ``x`` under this morphism.

        INPUT:

        - ``x`` -- an element in the domain of this morphism

        EXAMPLES::

            sage: K.<a> = GF(11^6).over(GF(11^3))                                       # needs sage.rings.finite_rings
            sage: V, i, j = K.free_module()                                             # needs sage.rings.finite_rings
            sage: i((0,1))                                                              # needs sage.rings.finite_rings
            a
        """
        cdef Element elt
        elt = self._f(v[0]) * self._basis[0]
        for i in range(1, self._degree):
            elt += self._f(v[i]) * self._basis[i]
        return self.codomain()(elt)


cdef class MapRelativeRingToFreeModule(Map):
    """
    Base class of the module isomorphism between a ring extension
    and a free module over one of its bases.

    TESTS::

        sage: K = GF(5^2).over()                                                        # needs sage.rings.finite_rings
        sage: V, i, j = K.free_module()                                                 # needs sage.rings.finite_rings
        sage: type(j)                                                                   # needs sage.rings.finite_rings
        <class 'sage.rings.ring_extension_morphism.MapRelativeRingToFreeModule'>
    """
    def __init__(self, E, K):
        r"""
        Initialize this morphism.

        INPUT:

        - ``E`` -- a ring extension

        - ``K`` -- a commutative ring; one base of ``E``

        TESTS::

            sage: K = GF(11^6).over(GF(11^3))                                           # needs sage.rings.finite_rings
            sage: V, i, j = K.free_module()                                             # needs sage.rings.finite_rings
            sage: j                                                                     # needs sage.rings.finite_rings
            Generic map:
              From: Field in z6 with defining polynomial x^2 + (10*z3^2 + z3 + 6)*x + z3 over its base
              To:   Vector space of dimension 2 over Finite Field in z3 of size 11^3
        """
        cdef Parent L, base

        self._degree = (<RingExtensionWithBasis>E)._degree_over(K)
        self._basis = [ (<RingExtensionElement>x)._backend for x in E.basis_over(K) ]
        f = backend_morphism(E.defining_morphism(K), forget='codomain')
        codomain = K ** self._degree
        Map.__init__(self, E.Hom(codomain))

        K = backend_parent(K)
        L = (<RingExtensionWithBasis>E)._backend

        # We compute the matrix of our isomorphism (over base)
        from sage.rings.ring_extension import common_base
        base = common_base(K, L, False)
        EK, iK, jK = K.free_module(base, map=True)
        EL, iL, jL = L.free_module(base, map=True)

        self._dimK = EK.dimension()
        self._iK = iK
        self._jL = jL

        M = [ ]
        for x in self._basis:
            for v in EK.basis():
                y = x * f(iK(v))
                M.append(jL(y))
        from sage.matrix.matrix_space import MatrixSpace
        self._matrix = MatrixSpace(base,len(M))(M).inverse_of_unit()

    def is_injective(self):
        r"""
        Return whether this morphism is injective.

        EXAMPLES::

            sage: K = GF(11^6).over(GF(11^3))                                           # needs sage.rings.finite_rings
            sage: V, i, j = K.free_module()                                             # needs sage.rings.finite_rings
            sage: j.is_injective()                                                      # needs sage.rings.finite_rings
            True
        """
        return True

    def is_surjective(self):
        r"""
        Return whether this morphism is injective.

        EXAMPLES::

            sage: K = GF(11^6).over(GF(11^3))                                           # needs sage.rings.finite_rings
            sage: V, i, j = K.free_module()                                             # needs sage.rings.finite_rings
            sage: j.is_surjective()                                                     # needs sage.rings.finite_rings
            True
        """
        return True

    cpdef Element _call_(self, x):
        r"""
        Return the image of ``x`` under this morphism.

        INPUT:

        - ``x`` -- an element in the domain of this morphism

        EXAMPLES::

            sage: K.<a> = GF(11^6).over(GF(11^3))                                       # needs sage.rings.finite_rings
            sage: V, i, j = K.free_module()                                             # needs sage.rings.finite_rings
            sage: j(a)                                                                  # needs sage.rings.finite_rings
            (0, 1)
        """
        coeffs = self.backend_coefficients(x)
        return self.codomain()(coeffs)

    cdef list backend_coefficients(self, RingExtensionElement x):
        r"""
        Return the coordinates of the image of ``x``
        as elements of the backend ring.

        INPUT:

        - ``x`` -- an element in the domain of this morphism

        TESTS::

            sage: K.<a> = GF(11^9).over(GF(11^3))                                       # needs sage.rings.finite_rings
            sage: V, i, j = K.free_module()                                             # needs sage.rings.finite_rings
            sage: j(a + 2*a^2)   # indirect doctest                                     # needs sage.rings.finite_rings
            (0, 1, 2)
        """
        cdef list coeffs = [ ]
        dK = self._dimK
        w = (self._jL(x._backend) * self._matrix).list()
        for i in range(self._degree):
            coeff = self._iK(w[i*dK:(i+1)*dK])
            coeffs.append(coeff)
        return coeffs
