r"""
Extension of rings

See :func:`RingExtension` for documentation.

AUTHOR:

- Xavier Caruso (2016)
"""

#############################################################################
#    Copyright (C) 2016 Xavier Caruso <xavier.caruso@normalesup.org>
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
from sage.categories.map import Map
from sage.rings.integer_ring import IntegerRing, ZZ
from sage.rings.ring import CommutativeRing, CommutativeAlgebra
from sage.rings.morphism import AlgebraFromMorphismHomomorphism
from sage.rings.infinity import Infinity

from sage.rings.algebra_from_morphism_element import AlgebraFMElement
from sage.rings.algebra_from_morphism_element import RingExtensionWithBasisElement


class AlgebraFromMorphism(CommutativeAlgebra, UniqueRepresentation):
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

        sage: from sage.rings.algebra_from_morphism import AlgebraFromMorphism
        sage: isinstance(E, AlgebraFromMorphism)
        True

    See :func:`RingExtension` for more documentation

    AUTHOR:

    - Xavier Caruso (2016)
    """
    Element = AlgebraFMElement

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
        base = defining_morphism.domain()
        ring = defining_morphism.codomain()

        if isinstance(ring, AlgebraFromMorphism):
            from sage.rings.morphism import backend_morphism
            defining_morphism = backend_morphism(defining_morphism, forget="codomain")
            coerce &= ring._coerce
            ring = ring._backend()

        CommutativeAlgebra.__init__(self, ZZ, category=CommutativeAlgebras(base))
        self._base = base
        self._ring = ring
        self._coerce = coerce
        self._defining_morphism = defining_morphism

        self._unset_coercions_used()
        f = AlgebraFromMorphismHomomorphism(self._base.Hom(self), defining_morphism)
        self.register_coercion(f)
        if coerce:
            self._populate_coercion_lists_(embedding = AlgebraFromMorphismHomomorphism(self.Hom(ring), ring.Hom(ring).identity()))

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
        if isinstance(self, AlgebraFromMorphism) or isinstance(other, AlgebraFromMorphism):
            from sage.rings.homset import AlgebraFromMorphismHomset
            return AlgebraFromMorphismHomset(self, other, category)

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
        if isinstance(other,AlgebraFromMorphism):
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
                    from sage.rings.algebra_from_morphism_constructor import RingExtension
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
        if isinstance(other, AlgebraFromMorphism):
            if self._coerce and other._coerce:
                return other.base().has_coerce_map_from(self._base) and self._ring.has_coerce_map_from(other._backend())

    def defining_morphism(self):
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
        return self._defining_morphism

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

    def gens(self):
        return tuple([ self(x) for x in self._ring.gens() ])

    def ngens(self):
        return len(self.gens())

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
        if isinstance(newbase, AlgebraFromMorphism):
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
        return AlgebraFromMorphism(defining_morphism, coerce)

    def intermediate_rings(self):
        L = [ self ]
        while isinstance(L[-1], AlgebraFromMorphism):
            L.append(L[-1].base())
        return L


# Finite free extensions (with a given basis)
#############################################


class RingExtensionWithBasis(AlgebraFromMorphism):
    Element = RingExtensionWithBasisElement

    def __init__(self, defining_morphism, basis, names=None, coerce=False):
        AlgebraFromMorphism.__init__(self, defining_morphism, coerce)
        self._basis = [ self(b)._backend() for b in basis ]
        if names is None:
            names = [ ]
            for b in self._basis:
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

    def relative_degree(self):
        return len(self._basis)

    def absolute_degree(self):
        return self.relative_degree() * self._base.absolute_degree()

    def basis(self):
        return [ self(x) for x in self._basis ]

    def is_finite(self):
        return True

    def is_free(self):
        return True

    def dimension(self):
        return len(self._basis)

    @cached_method
    def vector_space(self, K=None, map=True):
        if K is None:
            K = self._base
        if K is not self._base:
            raise NotImplemetedError
        d = self.relative_degree()
        if map:
            return K**d, MapVectorSpaceToRelativeField(self), MapRelativeFieldToVectorSpace(self)
        else:
            return K**d


class RingExtensionWithGen(RingExtensionWithBasis):
    def __init__(self, defining_morphism, gen, name, degree, coerce=False):
        self._name = str(name)
        names = [ "", self._name ] + [ "%s^%s" % (self._name, i) for i in range(2,degree) ]
        basis = [ gen ** i for i in range(degree) ]
        RingExtensionWithBasis.__init__(self, defining_morphism, basis, names, coerce)
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

    def gens(self):
        return (self(self._gen),)

    def _repr_(self):
        return "%s in %s with defining polynomial %s over its base" % (self._type, self._name, self.modulus())


class MapVectorSpaceToRelativeField(Map):
    def __init__(self, E):
        if not isinstance(E, RingExtensionWithBasis):
            raise TypeError("you must pass in a RingExtensionWithBasis")
        self._degree = E.relative_degree()
        self._basis = E._basis
        self._f = E.defining_morphism()
        domain = E.base_ring() ** self._degree
        parent = domain.Hom(E)
        Map.__init__(self, parent)

    def is_injective(self):
        return True

    def is_surjective(self):
        return True

    def _call_(self, v):
        return sum(self._f(v[i]) * self._basis[i] for i in range(self._degree))


class MapRelativeFieldToVectorSpace(Map):
    def __init__(self, E):
        if not isinstance(E, RingExtensionWithBasis):
            raise TypeError("you must pass in a RingExtensionWithBasis")
        self._degree = E.relative_degree()
        self._basis = E._basis
        K = E.base_ring()
        codomain = K ** self._degree
        parent = E.Hom(codomain)
        Map.__init__(self, parent)

        # We now perform some precomputations
        if isinstance(K, AlgebraFromMorphism):
            K = K._backend()
        L = E._backend()

        # We compute a common base for K and L
        def tower_finite_bases(ring):
            bases = [ ]
            base = ring
            while True:
                bases.append(base)
                try:
                    deg = base.relative_degree()
                except AttributeError:
                    try:
                        deg = base.degree()
                    except AttributeError:
                        break
                if deg is Infinity: break
                newbase = base.base_ring()
                if newbase is base: break
                base = newbase
            return bases
        bases_K = tower_finite_bases(K)
        bases_L = tower_finite_bases(L)
        base = None
        for b in bases_L:
            if b in bases_K:
                base = b
                break
        if base is None:
            raise NotImplementedError("unable to find a common base")

        # We compute the matrix of our isomorphism (over base)
        # EK, iK, jK = K.vector_space(base, with_maps=True)
        # EL, iL, jL = L.vector_space(base, with_maps=True)
        EK, iK, jK = K.absolute_vector_space()
        EL, iL, jL = L.absolute_vector_space()

        self._dimK = EK.dimension()
        self._iK = iK
        self._jL = jL

        M = [ ]
        f = E.defining_morphism()
        for x in self._basis:
            for v in EK.basis():
                y = x * f(iK(v))
                M.append(jL(y))
        from sage.matrix.matrix_space import MatrixSpace
        self._matrix = MatrixSpace(base,len(M))(M).inverse()

    def is_injective(self):
        return True

    def is_surjective(self):
        return True

    def _call_(self, x):
        dK = self._dimK
        w = (self._jL(x._backend()) * self._matrix).list()
        coeffs = [ ]
        for i in range(self._degree):
            coeff = self._iK(w[i*dK:(i+1)*dK])
            coeffs.append(coeff)
        return self.codomain()(coeffs)
