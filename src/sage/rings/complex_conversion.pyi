from sage.categories.map import Map
from sage.rings.complex_double import ComplexDoubleElement
from sage.rings.complex_mpfr import ComplexNumber

class CCtoCDF(Map[ComplexNumber, ComplexDoubleElement]):
    def _call_(self, x: ComplexNumber) -> ComplexDoubleElement: ...
