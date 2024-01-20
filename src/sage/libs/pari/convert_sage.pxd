from cypari2.gen cimport Gen
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational

cpdef gen_to_sage(Gen z, locals=*) noexcept

cpdef set_integer_from_gen(Integer self, Gen x) noexcept
cpdef Gen new_gen_from_integer(Integer self) noexcept
cpdef set_rational_from_gen(Rational self, Gen x) noexcept
cpdef Gen new_gen_from_rational(Rational self) noexcept

cpdef pari_is_prime(Integer p) noexcept
cpdef pari_is_prime_power(Integer q, bint get_data) noexcept
cpdef unsigned long pari_maxprime() noexcept
cpdef list pari_prime_range(long c_start, long c_stop, bint py_ints=*) noexcept
