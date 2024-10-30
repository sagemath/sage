from cypari2.gen cimport Gen

from sage.rings.real_mpfr cimport RealNumber

cpdef Gen new_gen_from_real_mpfr_element(RealNumber self)
cpdef bint set_real_mpfr_element_from_gen(RealNumber self, Gen x) noexcept
