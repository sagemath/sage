from sage.rings.polynomial.ore_polynomial_ring import OrePolynomialRing

class DifferentialPolynomialRing(OrePolynomialRing):
    def __init__(self, base_ring, morphism, derivation, name, sparse, category=None):
        if morphism is not None:
            raise NotImplementedError
        if self.Element is None:
            import sage.rings.polynomial.differential_polynomial_element
            self.Element = sage.rings.polynomial.differential_polynomial_element.DifferentialPolynomial_generic_dense
        OrePolynomialRing.__init__(self, base_ring, None, derivation, name, sparse, category)
