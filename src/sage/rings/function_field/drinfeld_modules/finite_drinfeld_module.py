from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.function_field.drinfeld_modules.drinfeld_module import DrinfeldModule

class FiniteDrinfeldModule(DrinfeldModule):

    def __init__(self, gen, category):

        # NOTE: There used to be no __init__ here (which was fine). I
        # added one to ensure that FiniteDrinfeldModule would always
        # have _frobenius_norm and _frobenius_trace attributes.

        super().__init__(gen, category)
        self._frobenius_norm = None
        self._frobenius_trace = None

    def frobenius_endomorphism(self):
        t = self.ore_variable()
        L = self._base_ring
        Fq = self._function_ring.base_ring()
        deg = L.over(Fq).degree(Fq)
        return self._Hom_(self, category=self.category())(t**deg)

    def frobenius_charpoly(self, var='T'):
        A = self._function_ring  # Fq[X]
        S = PolynomialRing(A, name=var)  # Fq[X][T]
        # Does not work when Fq is not a prime field...
        # chi = self._gen.reduced_charpoly()
        # return -chi(A.gen(), S.gen())
        return S([self.frobenius_norm(), -self.frobenius_trace(), 1])

    def frobenius_norm(self):
        self._check_rank_two()
        # Notations from Schost-Musleh:
        if self._frobenius_norm is None:
            n = self._base_ring.over(self._Fq).degree_over(self._Fq)
            d = self.characteristic().degree()
            m = n // d
            delta = self._gen[2]
            norm = self._base_ring.over(self._Fq)(delta).norm()
            self._frobenius_norm = ((-1)**n) * (self.characteristic()**m) / norm
        return self._frobenius_norm

    def frobenius_trace(self):
        self._check_rank_two()
        # Notations from Schost-Musleh:
        if self._frobenius_trace is None:
            n = self._base_ring.over(self._Fq).degree_over(self._Fq)
            B = self.frobenius_norm()
            t = self.ore_polring().gen()
            self._frobenius_trace = self.invert(t**n + (self(B) // t**n))
        return self._frobenius_trace

    def is_ordinary(self):
        self._check_rank_two()
        return not self.is_supersingular()

    def is_supersingular(self):
        self._check_rank_two()
        return self.characteristic().divides(self.frobenius_trace())
