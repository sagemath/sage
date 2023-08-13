from .multi_polynomial_libsingular import MPolynomialRing_libsingular, MPolynomial_libsingular
from .polynomial_ring_constructor import PolynomialRing
from .multi_polynomial_subring_element import MPolynomial_subring_element

class MPolynomial_subring(MPolynomialRing_libsingular):

    def __init__(self, parent_ring, gens):
        # assert parent ring is immutable ?
        self._hom=PolynomialRing(parent_ring.base_ring(), len(gens), "a").hom(gens)
        self.generators=[self._element_constructor_(gen) for gen in gens]
        self.parent_ring=parent_ring        
        self._ideal=self._hom.kernel()         
        self._zero_element = self._element_constructor_(0)         
        self._one_element = self._element_constructor_(1)
        self._homdomain = self._hom.domain()   




    #########################
    ##!!!!!! DISCUSS !!!!!!##
    #########################
    def _element_constructor_(self, element):
        assert element in self
        return MPolynomial_subring_element(self, element)


    def __repr__(self):
        """
        sage: R.<a,b>=QQ[]; S=R.subring_generated_by([a**2+b])
        sage: S
        Subring of Multivariate Polynomial Ring in a, b over Rational Field generated by [a^2 + b]
        """
        return f"Subring of {repr(self.parent_ring)} generated by {self.gens()}"

    def gens(self):
        return self.generators
    
    def ngens(self):
        return len(self.generators)
    
    def hilbert_series(self, algorithm="sage"):
        """
        sage: G.<x,y,z> = GF(5)[]
        sage: SG=G.subring_generated_by([x**2+y**3,z**3+y**2])
        sage: SG.hilbert_series()
        1/(t^6 - 2*t^3 + 1)
        """
        grading=[s.degree() for s in self.generators]
        return self._ideal.hilbert_series(grading=grading, algorithm=algorithm)
    
    def __contains__(self, other):
        """
        sage: R.<a,b,c>=ZZ[]
        sage: S=R.subring_generated_by([a**2, b**3+c])
        sage: b in S
        False
        sage: b**3 in S
        False
        sage: b**3 + c in S
        True
        sage: S2=R.subring_generated_by([a**2, b**3+c, c])
        sage: b**3 in S2
        True
        """
        try:
            self._hom.inverse_image(other)
            return True
        except ValueError as e:
            return False
    
    def random_element(self):
        return self._element_constructor_(self._hom(self._homdomain.random_element()))
    
    def _an_element_(self):
        return self._element_constructor_(self._hom(self._homdomain._an_element_()))



    



