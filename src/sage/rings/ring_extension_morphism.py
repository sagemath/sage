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
            self._base_map = None
        elif isinstance(defn, (list, tuple)):
            im_gens = []
            for x in backend_domain.gens():
                pol = from_backend_domain(x).polynomial()
                if base_map is not None:
                    pol = pol.map_coefficients(base_map)
                elif check and not codomain.has_coerce_map_from(pol.base_ring()):
                    raise ValueError("There is no coercion from base ring to codomain; try specifying base_map")
                im_gens.append(pol(defn))
            # There is no base map on the backend by assumption that the backends in a tower all have the same base and that base coerces into the absolute base of the tower
            self._backend = backend_domain.hom(im_gens, codomain=codomain, check=check)
            if check:
                for x, y in zip(domain.gens(), defn):
                    if self._backend(x._backend) != y:
                        raise ValueError("images do not define a valid homomorphism")
            self._im_gens = defn
            self._base_map = base_map
        else:
            raise TypeError(f"Unsupported type for hom: {type(im_gens)}")

            # We figure out what is the base
        #     if base_map is not None:
        #         base = base_map.domain()
        #         gens = domain.gens(base=base)
        #     else:
        #         base = domain
        #         gens = tuple([])
        #         while True:
        #             if len(gens) == len(defn):
        #                 break
        #             if len(gens) > len(defn) or base is base.base_ring():
        #                 raise ValueError("the number of images does not match the number of generators")
        #             gens += base.gens()
        #             base = base.base_ring()
        #     # We construct the backend morphism
        #     im_gens = [ codomain(x) for x in defn ]
        #     backend_bases = [ backend_domain ]
        #     b = backend_domain.base_ring()
        #     while b is not b.base_ring():
        #         backend_bases.append(b)
        #         b = b.base_ring()
        #     backend_bases.reverse()
        #     current_morphism = None
        #     for current_domain in backend_bases:
        #         current_im_gens = [ ]
        #         for x in current_domain.gens():
        #             pol = from_backend_domain(backend_domain(x)).polynomial(base=base)
        #             if base_map is not None:
        #                 pol = pol.map_coefficients(base_map)
        #             y = pol(im_gens)
        #             # Multivariate polynomials can have the wrong parent
        #             assert all(g.parent() is y.parent() for g in im_gens)
        #             current_im_gens.append(backend_element(y))
        #         current_morphism = current_domain.hom(current_im_gens, base_map=current_morphism, check=check)
        #     # We check that everything went well
        #     if check:
        #         for i in range(len(gens)):
        #             x = backend_element(domain(gens[i]))
        #             y = backend_element(im_gens[i])
        #             if current_morphism(x) != y:
        #                 raise ValueError("images do not define a valid homomorphism")
        #         coercion_morphism = backend_morphism(domain.defining_morphism(base))
        #         restriction_current_morphism = current_morphism * coercion_morphism
        #         if base_map is None:
        #             backend_base_map = coercion_morphism
        #         else:
        #             backend_base_map = backend_morphism(base_map)
        #             # the base map might be an automorphism of the base
        #             if backend_base_map.codomain() is coercion_morphism.domain():
        #                 backend_base_map = coercion_morphism * backend_base_map
        #             if backend_base_map.domain() is not restriction_current_morphism.domain():
        #                 phi = backend_base_map.domain().coerce_map_from(restriction_current_morphism.domain())
        #                 if phi is None:
        #                     msg = "Cannot coerce base map into correct domain:\n"
        #                     msg += f" Domain is {backend_base_map.domain()}\n"
        #                     msg += f" Needs to be {restriction_current_morphism.domain()}"
        #                     raise ValueError(msg)
        #                 backend_base_map = backend_base_map * phi
        #             if backend_base_map.codomain() is not restriction_current_morphism.codomain():
        #                 R = backend_base_map.codomain()
        #                 phi = restriction_current_morphism.codomain().coerce_map_from(R)
        #                 if phi is None:
        #                     # Try into the backend
        #                     back, from_back, to_back = backend_parent(R, map=True)
        #                     if back is not R and to_back is not None and restriction_current_morphism.codomain().has_coerce_map_from(back):
        #                         phi = restriction_current_morphism.codomain().coerce_map_from(back) * to_back
        #                 if phi is None:
        #                     msg = "Cannot coerce base map into correct codomain:\n"
        #                     msg += f" Codomain is {backend_base_map.codomain()}\n"
        #                     msg += f" Needs to be {restriction_current_morphism.codomain()}"
        #                     raise ValueError(msg)
        #                 backend_base_map = phi * backend_base_map
        #         differing = are_different_morphisms(restriction_current_morphism, backend_base_map)
        #         if differing:
        #             msg = "images do not define a valid homomorphism:\n"
        #             for x, y, z in differing:
        #                 if isinstance(x, str):
        #                     msg += f" different {x}:\n  {y}\n  {z}"
        #                 else:
        #                     msg += f" f({x}) = {y}\n g({x}) = {z}\n"
        #             raise ValueError(msg)
        #     self._backend = current_morphism
        #     self._im_gens = im_gens[:domain.ngens()]
        #     if base is domain.base_ring():
        #         self._base_map_construction = base_map
        #     else:
        #         self._base_map_construction = {
        #             'im_gens': defn[domain.ngens():],
        #             'base_map': base_map,
        #             'check': False
        #         }
        # else:
        #     raise TypeError("%s has type %s" % (defn, type(defn)))

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
        return self._base_map

        # domain = self.domain()
        # codomain = self.codomain()
        # base = domain.base_ring()
        # if base is base.base_ring():
        #     return None
        # base_map = self._base_map_construction
        # if base_map is False:
        #     if base is domain:
        #         base_map = None
        #     else:
        #         base_map = self * domain.coerce_map_from(base)
        # elif isinstance(base_map, dict):
        #     base_map = base.hom(**self._base_map_construction)
        # if base_map is None:
        #     return None
        # if (codomain.has_coerce_map_from(base) and
        #     not are_different_morphisms(backend_morphism(base_map),
        #                                 backend_morphism(codomain.coerce_map_from(base)))):
        #     return None
        # if base_map.codomain() is not self.codomain():
        #     base_map = base_map.extend_codomain(self.codomain())
        # return base_map

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
        if isinstance(other, RingExtensionHomomorphism):
            other = other._backend
        eq = not are_different_morphisms(self._backend, other)
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
            s += "with map on base ring"
            ss = self.base_map()._repr_defn()
            ss = re.sub('\nwith map on base ring:?$', '', ss, 0, re.MULTILINE)
            if ss != "":
                s += ":\n" + ss
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
        if isinstance(right, RingExtensionHomomorphism):
            backend_right = right._backend
        else:
            backend_right = right
        backend_middle, from_backend_middle, to_backend_middle = backend_parent(middle, map=True)
        backend = self._backend * to_backend_middle * backend_right
        from sage.rings.ring_extension import RingExtension_generic
        if isinstance(domain, RingExtension_generic):
            return RingExtensionHomomorphism(domain.Hom(codomain), backend)
        else:
            return backend

class RingExtensionHomomorphism_baseinclusion(RingMap):
    def _repr_type(self):
        return "Base injection"

    def _repr_defn(self):
        if self.codomain()._canonical_backend:
            return ""
        gens = self.domain().gens()
        s = ""
        for x in gens:
            s += f"{x} |--> {self(x)}\n"
        return s

    def _call_(self, x):
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
        <class 'sage.rings.ring_extension_morphism.RingExtensionBackendIsomorphism'>

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
        <class 'sage.rings.ring_extension_morphism.RingExtensionBackendReverseIsomorphism'>

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
        <class 'sage.rings.ring_extension_morphism.MapFreeModuleToRelativeRing'>

    """
    def __init__(self, E, K):
        r"""
        Initialize this morphism.

        INPUT:

        - ``E`` -- a ring extension

        - ``K`` -- a commutative ring; one base of ``E``

        TESTS::

            sage: K = GF(11^6).over(GF(11^3))
            sage: V, i, j = K.free_module()
            sage: i
            Generic map:
              From: Vector space of dimension 2 over Finite Field in z3 of size 11^3
              To:   Field in z6 with defining polynomial x^2 + (10*z3^2 + z3 + 6)*x + z3 over its base
        """
        self._degree = E.degree(K)
        self._basis = [ x._backend for x in E.basis_over(K) ]
        self._f = E._to_backend_morphism * E.defining_morphism(K)
        domain = K ** self._degree
        parent = domain.Hom(E)
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


class MapRelativeRingToFreeModule(Map):
    """
    Base class of the module isomorphism between a ring extension
    and a free module over one of its bases.

    TESTS::

        sage: K = GF(5^2).over()
        sage: V, i, j = K.free_module()
        sage: type(j)
        <class 'sage.rings.ring_extension_morphism.MapRelativeRingToFreeModule'>

    """
    def __init__(self, E, K):
        r"""
        Initialize this morphism.

        INPUT:

        - ``E`` -- a ring extension

        - ``K`` -- a commutative ring; one base of ``E``

        TESTS::

            sage: K = GF(11^6).over(GF(11^3))
            sage: V, i, j = K.free_module()
            sage: j
            Generic map:
              From: Field in z6 with defining polynomial x^2 + (10*z3^2 + z3 + 6)*x + z3 over its base
              To:   Vector space of dimension 2 over Finite Field in z3 of size 11^3
        """

        self._degree = E._degree_over(K)
        self._basis = [ x._backend for x in E.basis_over(K) ]
        defining_morphism = E._to_backend_morphism * E.defining_morphism(K)
        codomain = K ** self._degree
        Map.__init__(self, E.Hom(codomain))

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
