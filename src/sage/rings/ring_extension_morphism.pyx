#############################################################################
#    Copyright (C) 2019 Xavier Caruso <xavier.caruso@normalesup.org>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  http://www.gnu.org/licenses/
#****************************************************************************

from sage.structure.element cimport Element
from sage.categories.map import Map
from sage.rings.morphism cimport RingHomomorphism
from sage.rings.ring_extension import RingExtension_class, RingExtensionWithBasis
from sage.rings.ring_extension import _common_base


# Helper functions

def _backend_morphism(f):
    from sage.categories.map import FormalCompositeMap
    if not isinstance(f.domain(), RingExtension_class) and not isinstance(f.codomain(), RingExtension_class):
        return f
    elif isinstance(f, RingExtensionHomomorphism):
        return f._backend()
    elif isinstance(f, FormalCompositeMap):
        return _backend_morphism(f.then()) * _backend_morphism(f.first())
    else:
        raise NotImplementedError

def backend_morphism(f, forget="all"):
    try:
        g = _backend_morphism(f)
        if forget is None and (isinstance(f.domain(), RingExtension_class) or isinstance(f.codomain(), RingExtension_class)):
            g = RingExtensionHomomorphism(f.domain().Hom(f.codomain()), g)
        if forget == "domain" and isinstance(f.codomain(), RingExtension_class):
            g = RingExtensionHomomorphism(g.domain().Hom(f.codomain()), g)
        if forget == "codomain" and isinstance(f.domain(), RingExtension_class):
            g = RingExtensionHomomorphism(f.domain().Hom(g.codomain()), g)
    except NotImplementedError:
        g = f
        if (forget == "all" or forget == "domain") and isinstance(f.domain(), RingExtension_class):
            ring = f.domain()._backend()
            g = g * RingExtensionHomomorphism(ring.Hom(f.domain()), ring.Hom(ring).identity())
        if (forget == "all" or forget == "codomain") and isinstance(f.codomain(), RingExtension_class):
            ring = f.codomain()._backend()
            g = RingExtensionHomomorphism(f.codomain().Hom(ring), ring.Hom(ring).identity()) * g
    return g


# Classes

cdef class RingExtensionHomomorphism(RingHomomorphism):
    r"""
    Homomorphisms between extensions

    EXAMPLES:

        sage: F = GF(5^2)
        sage: K = GF(5^4)
        sage: L = GF(5^8)
        sage: E1 = RingExtension(K,F)
        sage: E2 = RingExtension(L,K)

    """
    def __init__(self, parent, backend):
        RingHomomorphism.__init__(self, parent)
        backend_domain = self.domain()
        if isinstance(backend_domain, RingExtension_class):
            backend_domain = backend_domain._backend()
        backend_codomain = self.codomain()
        if isinstance(backend_codomain, RingExtension_class):
            backend_codomain = backend_codomain._backend()
        backend = backend_morphism(backend)
        if backend.domain() is not backend_domain:
            raise TypeError("the domain of the backend morphism is not correct")
        if backend.codomain() is not backend_codomain:
            raise TypeError("the codomain of the backend morphism is not correct")
        self._backend_morphism = backend
        # We should probably allow for more general constructions but
        #   self._backend_morphism = backend_domain.Hom(backend_codomain)(*args, **kwargs)
        # does not work currently

    cpdef Element _call_(self, x):
        if isinstance(self.domain(), RingExtension_class):
            x = x._backend()
        y = self._backend_morphism(x)
        if isinstance(self.codomain(), RingExtension_class):
            y = self._codomain(y)
        return y

    def _backend(self):
        return self._backend_morphism

    cdef _update_slots(self, dict _slots):
        self._backend_morphism = _slots['_backend_morphism']
        RingHomomorphism._update_slots(self, _slots)

    cdef dict _extra_slots(self):
        slots = RingHomomorphism._extra_slots(self)
        slots['_backend_morphism'] = self._backend_morphism
        return slots


class MapVectorSpaceToRelativeField(Map):
    def __init__(self, E, K):
        if not isinstance(E, RingExtensionWithBasis):
            raise TypeError("you must pass in a RingExtensionWithBasis")
        self._degree = E.degree(K)
        self._basis = [ x._backend() for x in E.basis(K) ]
        self._f = E.defining_morphism(K)
        domain = K ** self._degree
        parent = domain.Hom(E)
        Map.__init__(self, parent)

    def is_injective(self):
        return True

    def is_surjective(self):
        return True

    def _call_(self, v):
        return sum(self._f(v[i]) * self._basis[i] for i in range(self._degree))


class MapRelativeFieldToVectorSpace(Map):
    def __init__(self, E, K):
        if not isinstance(E, RingExtensionWithBasis):
            raise TypeError("you must pass in a RingExtensionWithBasis")
        self._degree = E.degree(K)
        self._basis = [ x._backend() for x in E.basis(K) ]
        f = E.defining_morphism(K)
        codomain = K ** self._degree
        parent = E.Hom(codomain)
        Map.__init__(self, parent)

        if isinstance(K, RingExtension_class):
            K = K._backend()
        L = E._backend()

        # We compute the matrix of our isomorphism (over base)
        # EK, iK, jK = K.vector_space(base, map=True)
        # EL, iL, jL = L.vector_space(base, map=True)
        base = _common_base(K,L)
        EK, iK, jK = K.absolute_vector_space()
        EL, iL, jL = L.absolute_vector_space()

        self._dimK = EK.dimension()
        self._iK = iK
        self._jL = jL

        M = [ ]
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
