from sage.misc.lazy_attribute import lazy_attribute

from sage.categories.modules import Modules
from sage.categories.category_types import Category_over_base_ring
from sage.categories.homsets import Homsets

from sage.rings.polynomial.ore_polynomial_ring import OrePolynomialRing

class OreModules(Category_over_base_ring):
    @staticmethod
    def __classcall_private__(cls, ring, twist):
        if isinstance(twist, OrePolynomialRing):
            ore = twist.change_var('x')
            if ore.base_ring() is not ring:
                raise ValueError("base rings do not match")
        else:
            ore = OrePolynomialRing(ring, twist, names='x')
        return cls.__classcall__(cls, ore)

    def __init__(self, ore):
        base = ore.base_ring()
        Category_over_base_ring.__init__(self, base)
        self._ore = ore

    def super_categories(self):
        return [Modules(self.base())]

    def _repr_object_names(self):
        return "Ore modules over %s %s" % (self.base_ring(), self._ore._repr_twist())

    def Homsets(self):
        return Homsets()

    class SubcategoryMethods:

        @lazy_attribute
        def _ore(self):
            for cat in self.super_categories():
                if isinstance(cat, OreModules):
                    return cat._ore

        def ore_ring(self, var='x'):
            return self._ore.change_var(var)

        def twisting_morphism(self):
            return self._ore.twisting_morphism()

        def twisting_derivation(self):
            return self._ore.twisting_derivation()

    class ParentMethods:

        def twisting_morphism(self):
            return self.category().twisting_morphism()

        def twisting_derivation(self):
            return self.category().twisting_derivation()
