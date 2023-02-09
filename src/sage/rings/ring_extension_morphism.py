r"""
Morphisms between extension of rings

AUTHOR:

- Xavier Caruso (2019)
"""

#############################################################################
#    Copyright (C) 2019 Xavier Caruso <xavier.caruso@normalesup.org>
#                  2022 Julian Rüth <julian.rueth@fsfe.org>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  https://www.gnu.org/licenses/
#****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.structure.richcmp import op_EQ, op_NE

from sage.categories.map import Map
from sage.rings.morphism import RingMap
from sage.rings.ring_extension_conversion import backend_parent

# I don't trust the operator ==
def are_different_morphisms(f, g):
    r"""
    Return a list of triples encoding how ``f`` and ``g`` differ.

    If they have different domains, return ``[("domain", domain(f), domain(g))]``
    Otherwise, if they have different codomains, return ``[("codomain", codomain(f), codomain(g))]``.
    Otherwise, return the list of triples ``(x, f(x), g(x))``
    where `x` varies over the generators of the domain so that ``f(x) != g(x)``.

    INPUT:

    - ``f`` - a ring homomorphism or ``None``; if ``None``,
      we consider that ``f`` is a coercion map

    - ``g`` - a ring homomorphism or ``None``

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
    if f is None and g is None:
        return []
    elif f is None:
        b = g.domain()
        f = g.codomain().coerce_map_from(b)
        if f is None:
            return [("no coercion", b, g.codomain())]
    else:
        b = f.domain()
        if g is None:
            g = f.codomain().coerce_map_from(b)
            if g is None:
                return [("no coercion", b, f.codomain())]
        else:
            if b is not g.domain():
                return [("domain", f.domain(), g.domain())]
            elif f.codomain() is not g.codomain():
                return [("codomain", f.codomain(), g.codomain())]
    gens = tuple()
    while b is not b._base:
        gens += b.gens()
        b = b._base
    fvalues = [f(x) for x in gens]
    gvalues = [g(x) for x in gens]
    return [(x, y, z) for x, y, z in zip(gens, fvalues, gvalues) if y != z]

class RingExtensionHomomorphism(RingMap):
    r"""
    A class for ring homomorphisms between extensions.

    TESTS::

        sage: K.<a> = GF(5^2).over()
        sage: L.<b> = GF(5^4).over(K)
        sage: phi = L.hom([b^5, a^5])
        sage: phi
        Ring endomorphism of Field in b with defining polynomial x^2 + (3 - a)*x + a over its base
          Defn: b |--> (2 + a) + 2*b
                with map on base ring:
                a |--> 1 - a

        sage: type(phi)
        <class 'sage.rings.ring_extension_morphism.RingExtensionHomomorphism_with_category'>

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

        - ``check`` -- a boolean (default: ``True``); whether to check if
          the given data define a valid homomorphism

        TESTS::

            sage: S.<x> = QQ[]
            sage: T.<x,y> = QQ[]
            sage: f = T.hom([x^2, y^2])
            sage: f
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

        domain, codomain = self.domain(), self.codomain()
        backend_domain, from_backend_domain, _ = backend_parent(domain, map=True)

        # We construct the backend morphism
        if isinstance(defn, Map):
            ddom, dcodom = defn.domain(), defn.codomain()
            if base_map is not None:
                raise ValueError("base_map cannot be set when passing in the backend morphism")
            if ddom is not backend_domain:
                raise TypeError(f"the domain of the backend morphism is not correct, expected {backend_domain} but found {ddom}")
            from sage.rings.ring_extension import RingExtension_generic
            if isinstance(codomain, RingExtension_generic) and dcodom is codomain._backend:
                defn = codomain._from_backend_morphism * defn
            elif dcodom is not codomain:
                raise TypeError(f"the codomain of the backend morphism is not correct, expected {codomain} but found {dcodom}")
            self._backend = defn
            self._im_gens = None
            self._base_map_construction = False
            self._base_map = None
        elif isinstance(defn, (list, tuple)):
            # We figure out what is the base
            if base_map is not None:
                base = base_map.domain()
                gens = domain.gens(base)
            else:
                base = domain
                gens = tuple()
                while len(gens) != len(defn):
                    if len(gens) > len(defn) or base is base.base_ring():
                        raise ValueError("the number of images does not match the number of generators")
                    gens += base.gens()
                    base = base.base_ring()
            im_gens = []
            for x in backend_domain.gens():
                pol = from_backend_domain(x).polynomial(base=base)
                if base_map is not None:
                    pol = pol.map_coefficients(base_map)
                elif check and not codomain.has_coerce_map_from(pol.base_ring()):
                    raise ValueError("There is no coercion from base ring to codomain; try specifying base_map")
                im_gens.append(pol(defn))
            # There is no base map on the backend by assumption that the backends in a tower all have the same base and that base coerces into the absolute base of the tower
            self._backend = backend_domain.hom(im_gens, codomain=codomain, check=check)
            if check:
                for x, y in zip(domain.gens(base=base), defn):
                    if self._backend(x._backend) != y:
                        raise ValueError("images do not define a valid homomorphism")
            self._im_gens = defn
            if base_map is None:
                if base is domain or base is domain.base_ring():
                    self._base_map_construction = None
                else:
                    self._base_map_construction = (None,)
            else:
                if base is domain:
                    # This is really the same as the isinstance(defn, Map) case
                    self._base_map_construction = False
                else:
                    # At least one generator provided
                    self._base_map_construction = (base_map,)
        else:
            raise TypeError(f"Unsupported type for hom: {type(im_gens)}")

    def _repr_type(self):
        r"""
        Return a string that describes the type of this morphism.

        EXAMPLES::

            sage: K.<a> = GF(5^2).over()  # over GF(5)
            sage: f = K.hom([a^5])
            sage: f
            Ring endomorphism of Field in a with defining polynomial x^2 + 4*x + 2 over its base
              Defn: a |--> 1 - a

            sage: f._repr_type()
            'Ring'
        """
        return "Ring"

    def _call_(self, x):
        r"""
        Return the image of ``x`` under this morphism.

        INPUT:

        - ``x`` -- an element in the domain of this morphism

        EXAMPLES::

            sage: A.<sqrt2> = QQ.extension(x^2 - 2)
            sage: K.<sqrt2> = A.over()
            sage: f = K.hom([-sqrt2])
            sage: f
            Ring endomorphism of Field in sqrt2 with defining polynomial x^2 - 2 over its base
              Defn: sqrt2 |--> -sqrt2
            sage: f(sqrt2)
            -sqrt2

        TESTS::

            sage: a = QQ.random_element()
            sage: b = QQ.random_element()
            sage: f(a + b*sqrt2) == a - b*sqrt2
            True
        """
        # Using the _backend morphism we can map from the domain to the codomain.
        #
        # codomain ← backend
        #               ↑ (_backend)
        # domain   → backend
        #
        # Note that the domain might not actually be a ring extension but just
        # the base ring. This class is also used to implement maps into the
        # backend, i.e., codomain might not be a ring extension either.
        return self._backend(x._backend)
        # domain = self.domain()
        # from sage.rings.ring_extension import RingExtension_generic
        # if isinstance(domain, RingExtension_generic):
        #     x = domain._to_backend_morphism(x)

        # x = self._backend(x)

        # if isinstance(self.codomain(), RingExtension_generic):
        #     x = self.codomain()._from_backend_morphism(x)

        # return x

    @cached_method
    def base_map(self):
        r"""
        Return the base map of this morphism
        or just ``None`` if the base map is a coercion map.

        EXAMPLES::

            sage: F = GF(5)
            sage: K.<a> = GF(5^2).over(F)
            sage: L.<b> = GF(5^6).over(K)

        We define the absolute Frobenius of L::

            sage: FrobL = L.hom([b^5, a^5])
            sage: FrobL
            Ring endomorphism of Field in b with defining polynomial x^3 + (2 + 2*a)*x - a over its base
              Defn: b |--> (-1 + a) + (1 + 2*a)*b + a*b^2
                    with map on base ring:
                    a |--> 1 - a
            sage: FrobL.base_map()
            Ring morphism:
              From: Field in a with defining polynomial x^2 + 4*x + 2 over its base
              To:   Field in b with defining polynomial x^3 + (2 + 2*a)*x - a over its base
              Defn: a |--> 1 - a

        The square of ``FrobL`` acts trivially on K; in other words, it has
        a trivial base map::

            sage: phi = FrobL^2
            sage: phi
            Ring endomorphism of Field in b with defining polynomial x^3 + (2 + 2*a)*x - a over its base
              Defn: b |--> 2 + 2*a*b + (2 - a)*b^2
            sage: phi.base_map()

        """
        domain = self.domain()
        codomain = self.codomain()
        base = domain.base_ring()
        base_map = self._base_map_construction
        if base_map is None:
            # User provided images of generators and just relied on coercion for mapping the base
            return None
        elif base_map is False:
            # User provided an explicit morphism
            if base is domain:
                return None
            else:
                base_map = self * domain.coerce_map_from(base)
                if (codomain.has_coerce_map_from(base) and not are_different_morphisms(base_map, None)):
                    return None
        else:
            if base is base.base_ring():
                # User provided images of generators, and there's no way it could have gone below base
                return None
            # User provided images of generators, and they do affect the map on the base ring
            # (otherwise _base_map_construction would be set to None)
            n = domain.ngens()
            base_map = base_map[0]
            if len(self._im_gens) > n:
                base_map = base.hom(self._im_gens[n:], codomain=codomain, base_map=self._base_map_construction[0], category=self.category_for())
            if (codomain.has_coerce_map_from(base) and
                not are_different_morphisms(base_map._backend, None)):
                return None
        return base_map

    def _richcmp_(self, other, op):
        r"""
        Compare this element with ``other`` according to
        the rich comparison operator ``op``.

        INPUT:

        - ``other`` -- a morphism with the same codomain and codomain

        - ``op`` -- the comparison operator

        TESTS::

            sage: K.<a> = GF(5^2).over()   # over GF(5)
            sage: L.<b> = GF(5^6).over(K)

            sage: FrobK = K.hom([a^5])
            sage: FrobL = L.hom([b^5], base_map=FrobK)

            sage: FrobK^2 == End(K).identity()
            True
            sage: FrobL^6 == End(L).identity()
            True
        """
        eq = not are_different_morphisms(self, other)
        if op == op_EQ:
            return eq
        if op == op_NE:
            return not eq
        return NotImplemented

    def is_identity(self):
        r"""
        Return whether this morphism is the identity.

        EXAMPLES::

            sage: K.<a> = GF(5^2).over()   # over GF(5)
            sage: FrobK = K.hom([a^5])
            sage: FrobK.is_identity()
            False
            sage: (FrobK^2).is_identity()
            True
        """
        if self.domain() is not self.codomain():
            return False
        return not are_different_morphisms(self._backend, None)

    def is_injective(self):
        r"""
        Return whether this morphism is injective.

        EXAMPLES::

            sage: K = GF(5^5).over()
            sage: f = K.hom([K.gen()^5])
            sage: f.is_injective()
            True
        """
        return self._backend.is_injective()

    def is_surjective(self):
        r"""
        Return whether this morphism is surjective.

        EXAMPLES::

            sage: K = GF(5^10).over(GF(5^5))
            sage: iota = K.defining_morphism()
            sage: iota
            Base injection morphism:
              From: Finite Field in z5 of size 5^5
              To:   Field in z10 with defining polynomial x^2 + (2*z5^3 + 2*z5^2 + 4*z5 + 4)*x + z5 over its base
            sage: iota.is_surjective()
            False

            sage: K = GF(7).over(ZZ)
            sage: iota = K.defining_morphism()
            sage: iota
            Base injection morphism:
              From: Integer Ring
              To:   Finite Field of size 7 over its base
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
        for x, y in zip(gens, self._im_gens):
            s += f"{x} |--> {y}\n"
        if self.base_map() is not None:
            # Give images of generators as long as they're relevant
            base = self.domain().base_ring()
            codomain = self.codomain()
            gens = []
            imgs = []
            coimgs = []
            while base is not base._base:
                bgens = base.gens()
                gens.append(bgens)
                imgs.append([self(x) for x in bgens])
                co = codomain.coerce_map_from(base)
                if co is None:
                    coimgs.append([None] * len(bgens))
                else:
                    coimgs.append([co(x) for x in bgens])
                base = base._base
            i = len(gens)
            while i > 0 and imgs[i-1] == coimgs[i-1]:
                i -= 1
            if i > 0: # should always happen since base_map wasn't None
                s += "with map on base ring:"
                for X, Y in zip(gens[:i], imgs[:i]):
                    for x, y in zip(X, Y):
                        s += f"\n{x} |--> {y}"
        if s != "" and s[-1] == "\n":
            s = s[:-1]
        return s

    def _composition(self, right):
        r"""
        Return the composite ``self o right``.

        TESTS::

            sage: A.<sqrt5> = QQ.extension(x^2 - 5)
            sage: K.<sqrt5> = A.over()
            sage: f = K.hom([-sqrt5])
            sage: f
            Ring endomorphism of Field in sqrt5 with defining polynomial x^2 - 5 over its base
              Defn: sqrt5 |--> -sqrt5

            sage: f^2  # indirect doctest
            Ring endomorphism of Field in sqrt5 with defining polynomial x^2 - 5 over its base
              Defn: sqrt5 |--> sqrt5
        """
        domain = right.domain()
        middle = self.domain()
        codomain = self.codomain()
        from sage.rings.ring_extension import RingExtension_generic
        if isinstance(right, RingExtensionHomomorphism):
            backend_right = right._backend
        elif isinstance(domain, RingExtension_generic):
            backend_right = right * domain._from_backend_morphism
        else:
            backend_right = right
        backend_middle, from_backend_middle, to_backend_middle = backend_parent(middle, map=True)
        backend = self._backend * to_backend_middle * backend_right
        if isinstance(domain, RingExtension_generic):
            return RingExtensionHomomorphism(domain.Hom(codomain), backend)
        else:
            return backend

class RingExtensionHomomorphism_baseinclusion(RingMap):
    """
    Class for the maps from the base ring of a ring extension to the top ring.

    EXAMPLES::

        sage: K.<a> = GF(5^2).over()  # over GF(5)
        sage: TestSuite(K.defining_morphism()).run()
    """
    def _repr_type(self):
        """
        Return a string that describes the type of this morphism.

        EXAMPLES::

            sage: K.<a> = GF(5^2).over()  # over GF(5)
            sage: K.defining_morphism()._repr_type()
            'Base injection'
        """
        return "Base injection"

    def _repr_defn(self):
        r"""
        The definintion of this morphism, either empty (canonical case) or
        given by its action on generators (noncanonical case)

        EXAMPLES::

            sage: K.<a> = GF(5^2).over()  # over GF(5)
            sage: K.defining_morphism()._repr_defn()
            ''
            sage: K.<x> = FunctionField(GF(2)); _.<Y> = K[]
            sage: F.<y> = K.extension(Y^2 + Y + x + 1/x)
            sage: E = F.extension_constant_field(GF(2^4))
            sage: E.defining_morphism()._repr_defn()
            'y |--> y\nx |--> x\n1 |--> 1'
        """
        if self.codomain()._canonical_backend:
            return ""
        gens = self.domain().gens()
        s = ""
        for x in gens:
            s += f"{x} |--> {self(x)}\n"
        return s

    def _call_(self, x):
        """
        Return the image of ``x`` under this morphism.

        INPUT:

        - ``x`` -- an element in the domain of this morphism

        EXAMPLES::

            sage: K.<a> = GF(5^2).over()  # over GF(5)
            sage: L.<b> = GF(5^4).over(K)
            sage: f = L.defining_morphism()
            sage: f(a).parent() is L
            True
        """
        codomain = self.codomain()
        y = codomain._backend_defining_morphism(x)
        return codomain.element_class(codomain, y)

    def is_injective(self):
        r"""
        Return whether this morphism is injective.

        EXAMPLES::

            sage: K = GF(5^10).over(GF(5^5))
            sage: iota = K.defining_morphism()
            sage: iota
            Base injection morphism:
              From: Finite Field in z5 of size 5^5
              To:   Field in z10 with defining polynomial x^2 + (2*z5^3 + 2*z5^2 + 4*z5 + 4)*x + z5 over its base
            sage: iota.is_injective()
            True

            sage: K = GF(7).over(ZZ)
            sage: iota = K.defining_morphism()
            sage: iota
            Base injection morphism:
              From: Integer Ring
              To:   Finite Field of size 7 over its base
            sage: iota.is_injective()
            False
        """
        return self.codomain()._backend_defining_morphism.is_injective()

    def is_surjective(self):
        r"""
        Return whether this morphism is surjective.

        EXAMPLES::

            sage: K = GF(5^10).over(GF(5^5))
            sage: iota = K.defining_morphism()
            sage: iota
            Base injection morphism:
              From: Finite Field in z5 of size 5^5
              To:   Field in z10 with defining polynomial x^2 + (2*z5^3 + 2*z5^2 + 4*z5 + 4)*x + z5 over its base
            sage: iota.is_surjective()
            False

            sage: K = GF(7).over(ZZ)
            sage: iota = K.defining_morphism()
            sage: iota
            Base injection morphism:
              From: Integer Ring
              To:   Finite Field of size 7 over its base
            sage: iota.is_surjective()
            True
        """
        return self.codomain()._backend_defining_morphism.is_surjective()

class RingExtensionBackendIsomorphism(RingMap):
    r"""
    The isomorphism taking an element of the backend to its ring extension.

    TESTS::

        sage: K = GF(11^9).over(GF(11^3))
        sage: f = K.coerce_map_from(GF(11^9))
        sage: f
        Coercion morphism:
          From: Finite Field in z9 of size 11^9
          To:   Field in z9 with defining polynomial x^3 + (9*z3^2 + 5*z3 + 1)*x^2 + (4*z3 + 3)*x + 10*z3 over its base

        sage: type(f)
        <class 'sage.rings.ring_extension_morphism.RingExtensionBackendIsomorphism_with_category'>

        sage: TestSuite(f).run()
    """
    def __init__(self, parent):
        r"""
        Initialize this morphism.

        TESTS::

            sage: A.<a> = QQ.extension(x^2 - 5)
            sage: K = A.over()
            sage: K.coerce_map_from(A)
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

            sage: K.<a> = GF(5^2).over()  # over GF(5)
            sage: f = K.coerce_map_from(GF(5^2))
            sage: f
            Coercion morphism:
              From: Finite Field in z2 of size 5^2
              To:   Field in a with defining polynomial x^2 + 4*x + 2 over its base

            sage: f._repr_type()
            'Coercion'
        """
        return "Coercion"

    def _repr_defn(self):
        r"""
        Return the empty string since this morphism is canonical.

        EXAMPLES::

            sage: K.<a> = GF(5^2).over()  # over GF(5)
            sage: f = K.coerce_map_from(GF(5^2))
            sage: f
            Coercion morphism:
              From: Finite Field in z2 of size 5^2
              To:   Field in a with defining polynomial x^2 + 4*x + 2 over its base

            sage: f._repr_defn()
            ''
        """
        return ""

    def _call_(self, x):
        r"""
        Return the image of ``x`` under this morphism.

        INPUT:

        - ``x`` -- an element in the domain of this morphism

        EXAMPLES::

            sage: K.<a> = GF(5^2).over()  # over GF(5)
            sage: f = K.coerce_map_from(GF(5^2))
            sage: f(GF(5^2).gen())
            a
        """
        codomain = self.codomain()
        return codomain.element_class(codomain, x)


class RingExtensionBackendReverseIsomorphism(RingMap):
    r"""
    The isomorphism from a ring extension to its backend.

    TESTS::

        sage: K = GF(11^9).over(GF(11^3))
        sage: f = GF(11^9).convert_map_from(K)
        sage: f
        Canonical morphism:
          From: Field in z9 with defining polynomial x^3 + (9*z3^2 + 5*z3 + 1)*x^2 + (4*z3 + 3)*x + 10*z3 over its base
          To:   Finite Field in z9 of size 11^9

        sage: type(f)
        <class 'sage.rings.ring_extension_morphism.RingExtensionBackendReverseIsomorphism_with_category'>

        sage: TestSuite(f).run()

        sage: A.<a> = QQ.extension(x^2 - 5)
        sage: K = A.over()
        sage: A.convert_map_from(K)
        Canonical morphism:
          From: Field in a with defining polynomial x^2 - 5 over its base
          To:   Number Field in a with defining polynomial x^2 - 5

    """
    def _repr_type(self):
        r"""
        Return a string that describes the type of this morphism.

        EXAMPLES::

            sage: K.<a> = GF(5^2).over()  # over GF(5)
            sage: f = GF(5^2).convert_map_from(K)
            sage: f
            Canonical morphism:
              From: Field in a with defining polynomial x^2 + 4*x + 2 over its base
              To:   Finite Field in z2 of size 5^2

            sage: f._repr_type()
            'Canonical'
        """
        return "Canonical"

    def _repr_defn(self):
        r"""
        Return the empty string since this morphism is canonical.

        EXAMPLES::

            sage: K.<a> = GF(5^2).over()  # over GF(5)
            sage: f = GF(5^2).convert_map_from(K)
            sage: f
            Canonical morphism:
              From: Field in a with defining polynomial x^2 + 4*x + 2 over its base
              To:   Finite Field in z2 of size 5^2

            sage: f._repr_defn()
            ''
        """
        return ""

    def _call_(self, x):
        r"""
        Return the image of ``x`` under this morphism.

        INPUT:

        - ``x`` -- an element in the domain of this morphism

        EXAMPLES::

            sage: K.<a> = GF(5^2).over()  # over GF(5)
            sage: f = GF(5^2).convert_map_from(K)
            sage: f(a)
            z2
        """
        return x._backend

class MapFreeModuleToRelativeRing(Map):
    """
    Base class of the module isomorphism between a ring extension
    and a free module over one of its bases.

    TESTS::

        sage: K = GF(5^2).over()
        sage: V, i, j = K.free_module()
        sage: type(i)
        <class 'sage.rings.ring_extension_morphism.MapFreeModuleToRelativeRing_with_category'>
        sage: TestSuite(i).run()

    """
    def __init__(self, parent):
        r"""
        Initialize this morphism.

        INPUT:

        - ``parent`` -- the homset to a ring extension ``E`` from a free module over a base of ``E``.

        TESTS::

            sage: K = GF(11^6).over(GF(11^3))
            sage: V, i, j = K.free_module()
            sage: i
            Generic map:
              From: Vector space of dimension 2 over Finite Field in z3 of size 11^3
              To:   Field in z6 with defining polynomial x^2 + (10*z3^2 + z3 + 6)*x + z3 over its base
        """
        K = parent.domain().base_ring()
        E = parent.codomain()
        self._degree = parent.domain().dimension()
        self._basis = [ x._backend for x in E.basis_over(K) ]
        self._f = E._to_backend_morphism * E.defining_morphism(K)
        Map.__init__(self, parent)

    def is_injective(self):
        r"""
        Return whether this morphism is injective.

        EXAMPLES::

            sage: K = GF(11^6).over(GF(11^3))
            sage: V, i, j = K.free_module()
            sage: i.is_injective()
            True
        """
        return True

    def is_surjective(self):
        r"""
        Return whether this morphism is surjective.

        EXAMPLES::

            sage: K = GF(11^6).over(GF(11^3))
            sage: V, i, j = K.free_module()
            sage: i.is_surjective()
            True
        """
        return True

    def _call_(self, v):
        r"""
        Return the image of ``x`` under this morphism.

        INPUT:

        - ``x`` -- an element in the domain of this morphism

        EXAMPLES::

            sage: K.<a> = GF(11^6).over(GF(11^3))
            sage: V, i, j = K.free_module()
            sage: i((0,1))
            a
        """
        elt = self._f(v[0]) * self._basis[0]
        for i in range(1, self._degree):
            elt += self._f(v[i]) * self._basis[i]

        return self.codomain()._from_backend_morphism(elt)

    def _richcmp_(self, other, op):
        r"""
        Compare this element with ``other`` according to
        the rich comparison operator ``op``.

        INPUT:

        - ``other`` -- a morphism with the same codomain and codomain

        - ``op`` -- the comparison operator

        TESTS::

            sage: K = GF(5^2).over()
            sage: V, i, j = K.free_module()
            sage: i == i
            True
            sage: i == j
            False
        """
        eq = type(other) == type(self)
        if op == op_EQ:
            return eq
        if op == op_NE:
            return not eq
        return NotImplemented

class MapRelativeRingToFreeModule(Map):
    """
    Base class of the module isomorphism between a ring extension
    and a free module over one of its bases.

    TESTS::

        sage: K = GF(5^2).over()
        sage: V, i, j = K.free_module()
        sage: type(j)
        <class 'sage.rings.ring_extension_morphism.MapRelativeRingToFreeModule_with_category'>
        sage: TestSuite(j).run()
    """
    def __init__(self, parent):
        r"""
        Initialize this morphism.

        INPUT:

        - ``parent`` -- the homset from a ring extension ``E`` to a free module over a base of ``E``.

        TESTS::

            sage: K = GF(11^6).over(GF(11^3))
            sage: V, i, j = K.free_module()
            sage: j
            Generic map:
              From: Field in z6 with defining polynomial x^2 + (10*z3^2 + z3 + 6)*x + z3 over its base
              To:   Vector space of dimension 2 over Finite Field in z3 of size 11^3
        """
        E = parent.domain()
        K = parent.codomain().base_ring()
        self._degree = parent.codomain().dimension()
        self._basis = [ x._backend for x in E.basis_over(K) ]
        defining_morphism = E._to_backend_morphism * E.defining_morphism(K)
        Map.__init__(self, parent)

        K_backend, K_backend_to_K, K_to_K_backend = backend_parent(K, map=True)
        E_backend = E._backend

        # We compute the matrix of our isomorphism (over base)
        from sage.rings.ring_extension import common_base
        base = common_base(K_backend, E_backend, False)
        EK, iK, jK = K_backend.free_module(base, map=True)
        _, _, jL = E_backend.free_module(base, map=True)

        self._dimK = EK.dimension()
        self._iK = iK
        self._jL = jL

        M = [ ]
        for x in self._basis:
            for v in EK.basis():
                y = x * defining_morphism(K_backend_to_K(iK(v)))
                M.append(jL(y))
        from sage.matrix.matrix_space import MatrixSpace
        self._matrix = MatrixSpace(base, len(M))(M).inverse_of_unit()

    def is_injective(self):
        r"""
        Return whether this morphism is injective.

        EXAMPLES::

            sage: K = GF(11^6).over(GF(11^3))
            sage: V, i, j = K.free_module()
            sage: j.is_injective()
            True
        """
        return True

    def is_surjective(self):
        r"""
        Return whether this morphism is injective.

        EXAMPLES::

            sage: K = GF(11^6).over(GF(11^3))
            sage: V, i, j = K.free_module()
            sage: j.is_surjective()
            True
        """
        return True

    def _call_(self, x):
        r"""
        Return the image of ``x`` under this morphism.

        INPUT:

        - ``x`` -- an element in the domain of this morphism

        EXAMPLES::

            sage: K.<a> = GF(11^6).over(GF(11^3))
            sage: V, i, j = K.free_module()
            sage: j(a)
            (0, 1)
        """
        coeffs = self.backend_coefficients(x)

        _, from_backend, _ = backend_parent(self.codomain().base_ring(), map=True)
        coeffs = [from_backend(c) for c in coeffs]

        return self.codomain()(coeffs)

    def backend_coefficients(self, x):
        r"""
        Return the coordinates of the image of ``x``
        as elements of the base's backend ring.

        INPUT:

        - ``x`` -- an element in the domain of this morphism

        TESTS::

            sage: K.<a> = GF(11^9).over(GF(11^3))
            sage: V, i, j = K.free_module()
            sage: j(a + 2*a^2)   # indirect doctest
            (0, 1, 2)
        """
        coeffs = []
        dK = self._dimK
        w = (self._jL(x._backend) * self._matrix).list()
        for i in range(self._degree):
            coeff = self._iK(w[i*dK:(i+1)*dK])
            coeffs.append(coeff)
        return coeffs

    def _richcmp_(self, other, op):
        r"""
        Compare this element with ``other`` according to
        the rich comparison operator ``op``.

        INPUT:

        - ``other`` -- a morphism with the same codomain and codomain

        - ``op`` -- the comparison operator

        TESTS::

            sage: K = GF(5^2).over()
            sage: V, i, j = K.free_module()
            sage: j == j
            True
            sage: i == j
            False
        """
        eq = type(other) == type(self)
        if op == op_EQ:
            return eq
        if op == op_NE:
            return not eq
        return NotImplemented
