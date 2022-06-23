from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.function_field.drinfeld_modules.drinfeld_module import DrinfeldModule

class FiniteDrinfeldModule(DrinfeldModule):
    def frobenius_charpoly(self, var='x'):
        # Does not work when Fq is not a prime field...
        chi = self._gen.reduced_charpoly()
        A = self._polring
        S = PolynomialRing(A, name=var)
        return -chi(A.gen(), S.gen())

    def characteristic_polynomial(self):
        self._test_rank_two()
        FqXT = PolynomialRing(self.functions_ring(), 'T')
        return FqXT([self.frobenius_norm(), -self.frobenius_trace(), 1])

    def frobenius_norm(self):
        self._test_rank_two()
        # Notations from Schost-Musleh:
        n = self._L.over(self._Fq).degree_over(self._Fq)
        d = self.characteristic().degree()
        m = n // d
        norm = self._L.over(self._Fq)(self.delta()).norm()
        return ((-1)**n) * (self.characteristic()**m) / norm

    def frobenius_trace(self):
        self._test_rank_two()
        # Notations from Schost-Musleh:
        n = self._L.over(self._Fq).degree_over(self._Fq)
        B = self.frobenius_norm()
        t = self.ore_polring().gen()
        return self.invert(t**n + (self(B) // t**n))

    def is_ordinary(self):
        self._test_rank_two()
        return not self.is_supersingular()

    def is_supersingular(self):
        self._test_rank_two()
        return self.characteristic().divides(self.frobenius_trace())
