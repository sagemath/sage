from cypari2.gen cimport Gen
from sage.rings.complex_double cimport ComplexDoubleElement

cpdef ComplexDoubleElement pari_to_cdf(Gen g)

cpdef Gen new_gen_from_complex_double_element(ComplexDoubleElement self)
