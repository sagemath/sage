"""
General extensions of p-adic rings and fields; the base ring may also be an extension.

These are implemented as proxy parents, backed by an absolute extension.
"""

from sage.misc.cachefunc import cached_method
from .extension_element import pAdicGenericExtensionElement
from .padic_extension_generic import pAdicExtensionGeneric

# NotImplementedError: AlgebraFromMorphism shouldn't inherit from UniqueRepresentation
class pAdicGeneralExtension(AlgebraFromMorphism, pAdicExtensionGeneric):
    def __init__(self, exact_modulus, poly, prec, print_mode, shift_seed, names, implementation='FLINT'):
        self._exact_modulus = exact_modulus
        self._shift_seed = shift_seed
        self._implementation = 'General'
        defining_morphism = None # NotImplementedError
        AlgebraFromMorphism.__init__(self, defining_morphism, False)
        pAdicGeneric.__init__(self, poly, prec, print_mode, names, pAdicGenericExtensionElement)
        # Fix the following
        self._gen = None

    def teichmuller(self, x, prec=None):
        R = self._backend()
        x = R(x) if prec is None else R(x, prec)
        return self(R.teichmuller(x))

    def _prec_type(self):
        return self._backend()._prec_type()

    def is_field(self):
        return self._backend().is_field()

    def random_element(self, **kwds):
        return self(self._backend().random_element(**kwds))

    def residue_ring(self, n):
        raise NotImplementedError

    def residue_class_field(self):
        raise NotImplementedError

    def inertia_subring(self):
        raise NotImplementedError

    def gen(self, n=0):
        if n != 0:
            raise IndexError("only one generator")
        return self._gen

    def uniformizer(self):
        return self(self._backend().uniformizer())

    def uniformizer_pow(self, n):
        return self(self._backend().uniformizer_pow(n))

    def _uniformizer_print(self):
        return self._backend()._uniformizer_print()

    def _unram_print(self):
        return self._backend()._unram_print()

    def has_pth_root(self):
        return self._backend().has_pth_root()

    def has_root_of_unity(self, n):
        return self._backend().has_root_of_unity(self, n)

    
