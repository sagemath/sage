from cypari2.gen cimport Gen
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational

cpdef gen_to_sage(Gen z, locals=*)

cpdef set_integer_from_gen(Integer self, Gen x)
cpdef Gen new_gen_from_integer(Integer self)
cpdef set_rational_from_gen(Rational self, Gen x)
cpdef Gen new_gen_from_rational(Rational self)

cpdef pari_is_prime(Integer p)
cpdef pari_is_prime_power(Integer q, bint get_data)
cpdef unsigned long pari_maxprime() noexcept
cpdef list pari_prime_range(long c_start, long c_stop, bint py_ints=*)
