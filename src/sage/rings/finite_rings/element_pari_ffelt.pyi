from sage.rings.polynomial.polynomial_element import Polynomial

from sage.rings.finite_rings.element_base import FinitePolyExtElement
from sage.rings.integer import Integer

class FiniteFieldElement_pari_ffelt(FinitePolyExtElement):
    def lift(self) -> Integer: ...

    def log(self, base, order=None, check=False) -> Integer: ...

    def multiplicative_order(self) -> Integer: ...
