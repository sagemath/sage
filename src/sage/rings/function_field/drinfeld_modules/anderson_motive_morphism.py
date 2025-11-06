from sage.categories.map import Map

from sage.modules.ore_module_homspace import OreModule_homspace
from sage.modules.ore_module_morphism import OreModuleMorphism

from sage.rings.function_field.drinfeld_modules.morphism import DrinfeldModuleMorphism


# Morphisms between Anderson modules

class AndersonMotiveMorphism(OreModuleMorphism):
    def _repr_type(self):
        return "Anderson motive"

    def __init__(self, parent, im_gens, check=True):
        from sage.rings.function_field.drinfeld_modules.anderson_motive import AndersonMotive_drinfeld
        if isinstance(im_gens, DrinfeldModuleMorphism):
            domain = parent.domain()
            codomain = parent.codomain()
            if not isinstance(domain, AndersonMotive_drinfeld)\
            or domain.drinfeld_module() is not im_gens.codomain():
                raise ValueError("the domain must be the Anderson module of the codomain of the isogeny")
            if not isinstance(codomain, AndersonMotive_drinfeld)\
            or codomain.drinfeld_module() is not im_gens.domain():
                raise ValueError("the codomain must be the Anderson module of the domain of the isogeny")
            u = im_gens._ore_polynomial
            im_gens = {codomain.gen(0): u*domain.gen(0)}
            check = False
        OreModuleMorphism.__init__(self, parent, im_gens, check)

    def characteristic_polynomial(self, var='X'):
        chi = OreModuleMorphism.characteristic_polynomial(self, var)
        A = self.domain().function_ring()
        return chi.change_ring(A)

    charpoly = characteristic_polynomial


class AndersonMotive_homspace(OreModule_homspace):
    Element = AndersonMotiveMorphism


# Coercion maps

class DrinfeldToAnderson(Map):
    def __init__(self, parent, phi):
        Map.__init__(self, parent)
        self._phi = phi
        self._motive = parent.codomain()
        self._AK = self._motive.base_combined()

    def _call_(self, f):
        phi = self._phi
        r = phi.rank()
        phiT = phi.gen()
        coords = []
        for _ in range(r):
            coords.append([])
        while f:
            f, rem = f.right_quo_rem(phiT)
            for i in range(r):
                coords[i].append(rem[i])
        coords = [self._AK(c) for c in coords]
        return self._motive(coords)

class AndersonToDrinfeld(Map):
    def __init__(self, parent, phi):
        Map.__init__(self, parent)
        self._phi = phi
        self._Ktau = parent.codomain()

    def _call_(self, x):
        phi = self._phi
        r = phi.rank()
        phiT = phi.gen()
        S = self._Ktau
        xs = []
        for i in range(r):
            if x[i].denominator() != 1:
                raise ValueError("not in the Anderson motive")
            xs.append(x[i].numerator())
        ans = S.zero()
        d = max(xi.degree() for xi in xs)
        for j in range(d, -1, -1):
            ans = ans*phiT + S([xs[i][j] for i in range(r)])
        return ans
