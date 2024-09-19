from sage.modules.free_module import FreeModule_ambient
from sage.modules.free_module.ore_module_element import OreModule_element


class OreModule(FreeModule_ambient):
    Element = OreModule_element

    def __init__(self, f, twist=None, names=None):
        base = matrix.base_ring()
        rank = matrix.nrows()
        if matrix.ncols() != rank:
            raise ValueError("matrix must be square")
        if names is not None:
            if isinstance(names, (list, tuple)):
                if rank != len(names):
                    raise ValueError
            elif isinstance(names, str):
                names = [ names + str(i) for i in range(rank) ]
        self._names = names
        FreeModule_ambient.__init__(self, base, rank)
        self._pseudohom = FreeModule_ambient.pseudohom(self, f, twist)
        self._twist = twist

    def _repr_(self):
        return "Ore module of rank %s over %s" % (self.rank(), self.base_ring())

    def twisting_morphism(self):
        return self._pseudohom.twisting_morphism()

    def twisting_derivation(self):
        return self._pseudhom.twisting_derivation()

    def _element_constructor_(self, *args, **kwds):
        return self.element_class(self, *args, **kwds)

    def pseudohom(self):
        return self._pseudohom

    def basis(self):
        rank = self.rank()
        zero = self.base_ring().zero()
        one = self.base_ring().one()
        coeffs = [zero] * rank
        B = [ ]
        for i in range(rank):
            coeffs[i] = one
            B.append(self.element_class(self, coeffs))
            coeffs[i] = zero
        return B
