# sage_setup: distribution = sagemath-pari

from cypari2.gen cimport Gen
from sage.rings.real_double cimport RealDoubleElement

cpdef Gen new_gen_from_real_double_element(RealDoubleElement self)
