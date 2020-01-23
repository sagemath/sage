r"""
Extension of rings

See :func:`RingExtension` for documentation.

AUTHOR:

- Xavier Caruso (2016)
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

from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.commutative_rings import CommutativeRings
from sage.categories.pushout import pushout
from sage.categories.commutative_algebras import CommutativeAlgebras
from sage.rings.integer_ring import IntegerRing, ZZ
from sage.rings.infinity import Infinity
from sage.rings.ring import CommutativeRing, CommutativeAlgebra

from sage.rings.ring_extension_element import RingExtensionElement
from sage.rings.ring_extension_element import RingExtensionWithBasisElement


# Helper functions

def _common_base(K, L, degree=False):
    def tower_bases(ring):
        bases = [ ]
        base = ring
        deg = 1
        degrees = [ ]
        while True:
            bases.append(base)
            if degree:
                degrees.append(deg)
                try:
                    d = base.relative_degree()
                except AttributeError:
                    try:
                        d = base.degree()
                    except AttributeError:
                        break
                if d is Infinity: break
                deg *= d
            newbase = base.base_ring()
            if newbase is base: break
            base = newbase
        return bases, degrees
    bases_K, degrees_K = tower_bases(K)
    bases_L, degrees_L = tower_bases(L)
    base = None
    for iL in range(len(bases_L)):
        try:
            iK = bases_K.index(bases_L[iL])
            base = bases_L[iL]
            break
        except ValueError:
            pass
    if base is None:
        raise NotImplementedError("unable to find a common base")
    if degree:
        return base, degrees_K[iK], degrees_L[iL]
    else:
        return base


# General extensions
####################

class RingExtension_class(CommutativeAlgebra):
    r"""
    Create a ring extension

    Do not call this function directly
    Use :func:`RingExtension` instead

    INPUT:

    - ``defining_morphism`` -- a ring homomorphism

    - ``coerce`` -- boolean (specify whether defining_morphism
    is a coercion map or not)

    OUTPUT:

    The extension defined by ``defining_morphism``

    EXAMPLES::

        sage: K = GF(5^2)
        sage: L = GF(5^4)
        sage: E = RingExtension(L,K); E
        Finite Field in z4 of size 5^4 viewed as an algebra over its base

        sage: from sage.rings.ring_extension import RingExtension_class
        sage: isinstance(E, RingExtension_class)
        True

    See :func:`RingExtension` for more documentation

    AUTHOR:

    - Xavier Caruso (2016)
    """
    Element = RingExtensionElement

    def __init__(self, defining_morphism, coerce=False):
        r"""
        TESTS::

        The attribute _coerce indicates if the defining morphism is
        a coercion map::

            sage: K = GF(5^2)
            sage: E = RingExtension(K, K, K.frobenius_endomorphism())
            sage: E._coerce
            False

            sage: L = GF(5^4)
            sage: E1 = RingExtension(L,K)
            sage: E1._coerce
            True

        """
        from sage.rings.ring_extension_morphism import RingExtensionHomomorphism
        base = defining_morphism.domain()
        ring = defining_morphism.codomain()

        if isinstance(ring, RingExtension_class):
            from sage.rings.ring_extension_morphism import backend_morphism
            defining_morphism = backend_morphism(defining_morphism, forget="codomain")
            coerce &= ring._coerce
            ring = ring._backend()

        CommutativeAlgebra.__init__(self, ZZ, category=CommutativeAlgebras(base))
        self._base = base
        self._ring = ring
        self._coerce = coerce
        self._defining_morphism = defining_morphism

        self._unset_coercions_used()
        f = RingExtensionHomomorphism(self._base.Hom(self), defining_morphism)
        self.register_coercion(f)
        if coerce:
            self._populate_coercion_lists_(embedding = RingExtensionHomomorphism(self.Hom(ring), ring.Hom(ring).identity()))

    def from_base_ring(self, r):
        r"""
        Return the canonical embedding of ``r`` into ``self``.

        INPUT:

        - ``r`` -- an element of the base of the ring of this extension

        EXAMPLES::

            sage: k = GF(5)
            sage: K = GF(5^2); z2 = K.gen()
            sage: L = GF(5^4); z4 = L.gen()
            sage: E = RingExtension(L,K); E
            Finite Field in z4 of size 5^4 viewed as an algebra over its base

            sage: E.base()
            Finite Field in z2 of size 5^2
            sage: E._backend().base()
            Finite Field of size 5

            sage: x = E.from_base_ring(k(2)); x
            2
            sage: x.parent()
            Finite Field in z4 of size 5^4 viewed as an algebra over its base

            sage: x = E.from_base_ring(z2); x
            z4^3 + z4^2 + z4 + 3

        TESTS::

            sage: R.<t> = K[]
            sage: F = RingExtension(R, K, K.frobenius_endomorphism())
            sage: F.from_base_ring(z2)
            4*z2 + 1
            sage: F.defining_morphism()(z2)
            4*z2 + 1
        """
        if r not in self._base:
            raise TypeError("%s is not an element of the base of %s (= %s)" % (r, self._ring, self._base))
        return self.element_class(self, r)

    def _Hom_(self, other, category):
        from sage.rings.ring_extension_homset import RingExtensionHomset
        return RingExtensionHomset(self, other, category)

    def hom(self, im_gens, codomain=None, base_map=None, check=True, category=None):
        from sage.rings.ring_extension_homset import RingExtensionHomset
        from sage.rings.ring_extension_morphism import RingExtensionHomomorphism
        if codomain is None:
            from sage.structure.sequence import Sequence
            codomain = Sequence(im_gens).universe()
        parent = RingExtensionHomset(self, codomain, category)
        return RingExtensionHomomorphism(parent, im_gens, base_map, check)

    def _pushout_(self, other):
        r"""
        Construct the pushout of this extension and ``other``
        This is called by :func:`sage.categories.pushout.pushout`.

        ALGORITHM:

        If ``other`` coerces to the base of ``self``, we return ``self``
        Similarly, if ``self`` coerces to the base of ``other``, we return
        ``other``

        Otherwise, if ``self`` or ``other`` are extensions defined by
        coercion morphisms, we compute the pushout of ``self._backend()``
        and ``other._backend()`` and call it ``ring``. We then check
        whether there exists a coercion map between ``self.base()`` and
        ``other.base()`` and, if so, we return the extension ``ring/base``
        where ``base`` the domain of the coercion map.

        In all other cases, return None

        .. TODO::

            Handle the case where defining morphisms are not
            coercion maps

            If no coercion map is discovered between ``self.base()``
            and ``other.base()``, compute the greatest base ``base``
            on which ``self.base()`` and ``other.base()`` are defined
            and return the extension ``ring/base``

        TESTS::

            sage: K1 = GF(3^4); L1 = GF(3^8);  E1 = RingExtension(L1,K1)
            sage: K2 = GF(3^2); L2 = GF(3^16); E2 = RingExtension(L2,K2)
            sage: E2.has_coerce_map_from(E1)
            True
            sage: E2._pushout_(E1) is E2
            True
            sage: E1._pushout_(E2) is E2
            True

            sage: E1._pushout_(L2)

        """
        if self._base.has_coerce_map_from(other):
            return self
        if isinstance(other,RingExtension_class):
            if other._backend().has_coerce_map_from(self):
                return other
            if self._coerce and other._coerce:
                ring = pushout(self._ring, other._backend())
                sbase = self._base
                obase = other._base
                base = None
                if sbase.has_coerce_map_from(obase):
                    base = obase
                elif obase.has_coerce_map_from(sbase):
                    base = sbase
                if base is not None:
                    from sage.rings.ring_extension_constructor import RingExtension
                    return RingExtension(ring, base)

    def _repr_(self):
        r"""
        Return a string representation of this extension.

        EXAMPLES::

            sage: K = GF(5^2)
            sage: L = GF(5^4)
            sage: E = RingExtension(L,K)
            sage: E
            Finite Field in z4 of size 5^4 viewed as an algebra over its base

            sage: E._repr_()
            'Finite Field in z4 of size 5^4 viewed as an algebra over its base'
        """
        return "%s viewed as an algebra over its base" % self._ring

    def _latex_(self):
        r"""
        Return a latex representation of this extension.

        EXAMPLES::

            sage: K = GF(5^2)
            sage: L = GF(5^4)
            sage: E = RingExtension(L,K)
            sage: latex(E)
            \Bold{F}_{5^{4}}/\Bold{F}_{5^{2}}
        """
        from sage.misc.latex import latex
        return "%s/%s" % (latex(self._ring), latex(self._base))

    def _coerce_map_from_(self, other):
        r"""
        An extension L1/K1 coerces to another extension L2/K2 when L1 coerces
        to L2 and K2 coerces to K1 (be careful at the direction) and the
        defining morphisms are coercion maps.

        .. TODO::

            Relax the condition on defining morphisms (and only require
            that they agree on K2)

        TESTS::

            sage: K1 = GF(3^4); L1 = GF(3^8);  E1 = RingExtension(L1,K1)
            sage: K2 = GF(3^2); L2 = GF(3^16); E2 = RingExtension(L2,K2)
            sage: E2.has_coerce_map_from(E1)  # indirect doctest
            True

        In the next example, the defining morphisms are not coercion maps.
        So there is no coercion between the extensions (through it probably
        should).

            sage: E1p = RingExtension(L1, K1, K1.frobenius_endomorphism())
            sage: E2p = RingExtension(L2, K2, K2.frobenius_endomorphism())
            sage: E2p.has_coerce_map_from(E1p)  # indirect doctest
            False
        """
        if isinstance(other, RingExtension_class):
            if self._coerce and other._coerce:
                return other.base().has_coerce_map_from(self._base) and self._ring.has_coerce_map_from(other._backend())

    def defining_morphism(self, base=None):
        r"""
        Return the defining morphism of this extension

        EXAMPLES::

            sage: K = GF(5^2); z2 = K.gen()
            sage: L = GF(5^4); z4 = L.gen()
            sage: E = RingExtension(L,K); E
            Finite Field in z4 of size 5^4 viewed as an algebra over its base
            sage: E.defining_morphism()
            Ring morphism:
              From: Finite Field in z2 of size 5^2
              To:   Finite Field in z4 of size 5^4
              Defn: z2 |--> z4^3 + z4^2 + z4 + 3

            sage: E2 = RingExtension(L, K, K.frobenius_endomorphism()); E2
            Finite Field in z4 of size 5^4 viewed as an algebra over its base
            sage: E2.defining_morphism()
            Composite map:
              From: Finite Field in z2 of size 5^2
              To:   Finite Field in z4 of size 5^4
              Defn:   Frobenius endomorphism z2 |--> z2^5 on Finite Field in z2 of size 5^2
                    then
                      Ring morphism:
                      From: Finite Field in z2 of size 5^2
                      To:   Finite Field in z4 of size 5^4
                      Defn: z2 |--> z4^3 + z4^2 + z4 + 3
        """
        from sage.rings.ring_extension_morphism import RingExtensionHomomorphism
        from sage.rings.ring_extension_morphism import backend_morphism

        if base is None or base is self._base:
            return self._defining_morphism
        ring = self._ring
        f = ring.Hom(ring).identity()
        b = self
        while b is not base:
            f = f * backend_morphism(b.defining_morphism())
            if b is b.base_ring():
                raise ValueError("(%s) is not defined over (%s)" % (self, base))
            b = b.base_ring()
        if isinstance(base, RingExtension_class):
            f = RingExtensionHomomorphism(base.Hom(ring), f)
        return f

    def base(self):
        r"""
        Return the base of this extension

        EXAMPLES::

            sage: K = GF(5^2); z2 = K.gen()
            sage: L = GF(5^4); z4 = L.gen()
            sage: E = RingExtension(L,K); E
            Finite Field in z4 of size 5^4 viewed as an algebra over its base
            sage: E.base()
            Finite Field in z2 of size 5^2
            sage: E.base() is K
            True

        Note that the base of the extension is generally different
        from the base of the ring of the extension::

            sage: E._backend().base()
            Finite Field of size 5
        """
        return self._base

    def _backend(self):
        r"""
        Return the top ring of this extension

        EXAMPLES::

            sage: K = GF(5^2); z2 = K.gen()
            sage: L = GF(5^4); z4 = L.gen()
            sage: E = RingExtension(L,K); E
            Finite Field in z4 of size 5^4 viewed as an algebra over its base
            sage: E._backend()
            Finite Field in z4 of size 5^4
            sage: E._backend() is L
            True
        """
        return self._ring

    def _an_element_(self):
        r"""
        Return an element of this extension

        This element is constructed by converting an element of
        the ring of this extension to the extension itself

        TESTS::

            sage: K = GF(5^2)
            sage: L = GF(5^4)
            sage: E = RingExtension(L,K); E
            Finite Field in z4 of size 5^4 viewed as an algebra over its base

            sage: x = E.an_element()  # indirect doctest
            sage: x
            0
            sage: x.parent()
            Finite Field in z4 of size 5^4 viewed as an algebra over its base

            sage: a = E._backend().an_element()
            sage: a == x
            True
        """
        elt = self._ring.an_element()
        return self.element_class(self, elt)

    def gens(self, base=None):
        return tuple([ self(x) for x in self._ring.gens() ])

    def ngens(self, base=None):
        return len(self.gens(base))

    def gen(self):
        r"""
        Return a generator of this extension

        This element is constructed by converting a generator of
        the ring of this extension to the extension itself

        TESTS::

            sage: K = GF(5^2)
            sage: L = GF(5^4)
            sage: E = RingExtension(L,K); E
            Finite Field in z4 of size 5^4 viewed as an algebra over its base

            sage: x = E.gen()
            sage: x
            z4
            sage: x.parent()
            Finite Field in z4 of size 5^4 viewed as an algebra over its base

            sage: a = E._backend().gen()
            sage: a == x
            True
        """
        return self.gens()[0]

    def random_element(self):
        r"""
        Return a random element of this extension

        This element is constructed by converting a random element of
        the ring of this extension to the extension itself

        TESTS::

            sage: K = GF(5^2)
            sage: L = GF(5^4)
            sage: E = RingExtension(L,K); E
            Finite Field in z4 of size 5^4 viewed as an algebra over its base

            sage: x = E.random_element(); x  # random
            0
            sage: x.parent()
            Finite Field in z4 of size 5^4 viewed as an algebra over its base
        """
        elt = self._ring.random_element()
        return self.element_class(self, elt)

    def is_field(self, proof=False):
        return self._ring.is_field(proof=proof)

    def scalar_restriction(self, newbase):
        r"""
        Return the scalar restriction of this extension to
        the new base ``newbase``.

        INPUT:

        - ``newbase`` -- it can be either
          * a ring which coerces to the base of ``self`` or
          * an extension whose ring coerces to the base of ``self`` or
          * a morphism whose codomain coerces to the base of ``self``

        OUTPUT:

        The scalar restriction of ``self`` to the ``newbase``.

        It is defined as follows.
        In each above three cases, we get a ring homorphism ``f``
        from a domain ``domain`` to the base of ``self``.
        The scalar restriction of ``self`` to ``newbase`` is the
        extension ``self._backend()/domain`` with defining morphism the
        compositum of ``f`` with the defining morphism of ``self``

        EXAMPLES::

            sage: k = GF(5^2)
            sage: K = GF(5^4)
            sage: L = GF(5^8)
            sage: E = RingExtension(L,K); E
            Finite Field in z8 of size 5^8 viewed as an algebra over its base

            sage: E1 = E.scalar_restriction(k); E1
            Finite Field in z8 of size 5^8 viewed as an algebra over its base
            sage: E1.defining_morphism()
            Ring morphism:
              From: Finite Field in z2 of size 5^2
              To:   Finite Field in z8 of size 5^8
              Defn: z2 |--> z8^7 + z8^6 + z8^5 + z8^4 + 3*z8^3 + 3*z8^2 + z8 + 3

            sage: Frob = k.frobenius_endomorphism()
            sage: E2 = E.scalar_restriction(Frob); E2
            Finite Field in z8 of size 5^8 viewed as an algebra over its base
            sage: E2.defining_morphism()
            Composite map:
              From: Finite Field in z2 of size 5^2
              To:   Finite Field in z8 of size 5^8
              Defn:   Frobenius endomorphism z2 |--> z2^5 on Finite Field in z2 of size 5^2
                    then
                      Ring morphism:
                      From: Finite Field in z2 of size 5^2
                      To:   Finite Field in z4 of size 5^4
                      Defn: z2 |--> z4^3 + z4^2 + z4 + 3
                    then
                      Ring morphism:
                      From: Finite Field in z4 of size 5^4
                      To:   Finite Field in z8 of size 5^8
                      Defn: z4 |--> 4*z8^7 + z8^6 + 3*z8^4 + z8^3 + z8^2 + 4

            sage: F = RingExtension(K, k, Frob)
            sage: E3 = E.scalar_restriction(F); E3
            Finite Field in z8 of size 5^8 viewed as an algebra over its base
            sage: E3.defining_morphism()
            Composite map:
              From: Finite Field in z2 of size 5^2
              To:   Finite Field in z8 of size 5^8
              Defn:   Frobenius endomorphism z2 |--> z2^5 on Finite Field in z2 of size 5^2
                    then
                      Ring morphism:
                      From: Finite Field in z2 of size 5^2
                      To:   Finite Field in z4 of size 5^4
                      Defn: z2 |--> z4^3 + z4^2 + z4 + 3
                    then
                      Ring morphism:
                      From: Finite Field in z4 of size 5^4
                      To:   Finite Field in z8 of size 5^8
                      Defn: z4 |--> 4*z8^7 + z8^6 + 3*z8^4 + z8^3 + z8^2 + 4

            sage: E2 is E3
            True
        """
        coerce = self._coerce
        if isinstance(newbase, RingExtension_class):
            if coerce:
                coerce = newbase._coerce
            newbase = newbase.defining_morphism()
        elif isinstance(newbase, CommutativeRing):
            if self._base.has_coerce_map_from(newbase):
                newbase = self._base.coerce_map_from(newbase)
            else:
                raise TypeError("No coercion map from %s to %s" % (newbase, self._base))
        else:
            domain = newbase.domain()
            codomain = newbase.codomain()
            if coerce:
                coerce = (newbase == codomain.coerce_map_from(domain))
            if self._base.has_coerce_map_from(codomain):
                newbase = newbase.post_compose(self._base.coerce_map_from(codomain))
            else:
                raise TypeError("No coercion map from %s to %s" % (codomain, self._base))
        defining_morphism = self._defining_morphism.pre_compose(newbase)
        return RingExtension_class(defining_morphism, coerce)

    def intermediate_rings(self):
        L = [ self ]
        while isinstance(L[-1], RingExtension_class):
            L.append(L[-1].base())
        return L


# Finite free extensions
########################

class RingExtensionWithBasis(RingExtension_class):
    Element = RingExtensionWithBasisElement

    def __init__(self, defining_morphism, basis, names=None, coerce=False, check=True):
        RingExtension_class.__init__(self, defining_morphism, coerce)
        self._basis = [ self(b) for b in basis ]
        if names is None:
            names = [ ]
            for b in self._basis:
                b = b._backend()
                if b == 1:
                    names.append("")
                if b._is_atomic():
                    names.append(str(b))
                else:
                    names.append("(%s)" % b)
        else:
            if len(names) != len(self._basis):
                raise ValueError("the number of names does not match the cardinality of the basis")
        self._names = names
        if check:
            try:
                _ = self.vector_space()
            except ZeroDivisionError:
                raise ValueError("the given family is not a basis")

    def degree(self, base):
        if base is self:
            return ZZ(1)
        elif base is self._base:
            return len(self._basis)
        else:
            try:
                deg = self._base.degree(base)
            except TypeError:
                if base is not self._base.base_ring():
                    raise NotImplementedError
                deg = self._base.relative_degree()
            return len(self._basis) * deg

    def relative_degree(self):
        return len(self._basis)

    def absolute_degree(self):
        return self.relative_degree() * self._base.absolute_degree()

    def basis(self, base=None):
        if base is self:
            return [ self.one() ]
        elif base is None or base is self._base:
            return self._basis[:]
        else:
            b = self._base.basis(base)
            return [ x*y for y in b for x in self._basis ]

    def is_finite(self):
        return True

    def is_free(self):
        return True

    def dimension(self, base):
        return self.degree(base)

    @cached_method
    def vector_space(self, base=None, map=True):
        if base is None:
            base = self._base
        d = self.degree(base)
        if map:
            from sage.rings.ring_extension_morphism import MapVectorSpaceToRelativeField, MapRelativeFieldToVectorSpace
            return base**d, MapVectorSpaceToRelativeField(self, base), MapRelativeFieldToVectorSpace(self, base)
        else:
            return base**d


class RingExtensionWithGen(RingExtensionWithBasis):
    def __init__(self, defining_morphism, gen, name, coerce=False, check=True):
        self._name = str(name)
        _, deg_domain, deg_codomain = _common_base(defining_morphism.domain(), defining_morphism.codomain(), degree=True)
        degree = ZZ(deg_codomain / deg_domain)
        names = [ "", self._name ] + [ "%s^%s" % (self._name, i) for i in range(2,degree) ]
        basis = [ gen ** i for i in range(degree) ]
        try:
            RingExtensionWithBasis.__init__(self, defining_morphism, basis, names, coerce, check)
        except ValueError:
            raise ValueError("the given element is not a generator")
        self._gen = self(gen)._backend()
        self._type = "Ring"
        # if self._ring in Fields():
        #     self._type = "Field"

    def modulus(self, var='x'):
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        _, _, j = self.vector_space()
        d = self.relative_degree()
        coeffs = [ -c for c in j(self._gen**d) ] + [ 1 ]
        S = PolynomialRing(self._base, name=var)
        return S(coeffs)

    def gens(self, base=None):
        if base is None:
            return (self(self._gen),)
        gens = tuple([])
        b = self
        while b is not base:
            gens += b.gens()
            if b is b.base_ring():
                raise ValueError("(%s) is not defined over (%s)" % (self, base))
            b = b.base_ring()
        return gens

    def _repr_(self, base=None):
        return "%s in %s with defining polynomial %s over its base" % (self._type, self._name, self.modulus())
