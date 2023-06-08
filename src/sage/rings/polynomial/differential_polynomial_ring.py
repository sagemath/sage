from sage.rings.polynomial.ore_polynomial_ring import OrePolynomialRing

class DifferentialPolynomialRing(OrePolynomialRing):
    def __init__(self, base_ring, morphism, derivation, name, sparse, category=None):
        r"""
        Initialize ``self``.

        INPUT:

        - ``base_ring`` -- a commutative ring

        - ``derivation`` -- a derivation of the base ring

        - ``name`` -- string or list of strings representing the name of
          the variables of ring

        - ``sparse`` -- boolean (default: ``False``)

        - ``category`` -- a category

        EXAMPLES::

            sage: R.<t> = ZZ[]
            sage: derivation = t*R.derivation()
            sage: from sage.rings.polynomial.differential_polynomial_ring import DifferentialPolynomialRing
            sage: S.<D> = DifferentialPolynomialRing(R,derivation)
            sage: S.category()
            Category of algebras over Univariate Polynomial Ring in t over Integer Ring
            sage: D*t
            t*D + t
        """
        if morphism is not None:
            raise NotImplementedError
        if self.Element is None:
            import sage.rings.polynomial.differential_polynomial_element
            self.Element = sage.rings.polynomial.differential_polynomial_element.DifferentialPolynomial_generic_dense
        OrePolynomialRing.__init__(self, base_ring, None, derivation, name, sparse, category)



