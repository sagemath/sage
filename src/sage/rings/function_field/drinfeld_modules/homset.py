from sage.categories.homset import Homset, Hom

class DrinfeldModuleHomset(Homset):
    def __init__(self, X, Y, base=None, check=True):
        if X.category() != Y.category() \
                and not isinstance(X.category(), DrinfeldModules):
            raise NotImplementedError('Drinfeld modules must be in the same category')
        Homset.__init__(self, X, Y, base=base, check=check)

    def __contains__(self, candidate):
        phi = self.domain()
        psi = self.codomain()
        if candidate not in phi.ore_polring():
            raise TypeError('morphism must be in the Ore polynomial ring')
        return candidate * phi.gen() == psi.gen() * candidate
