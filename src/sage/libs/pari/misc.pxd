# sage_setup: distribution = sagemath-pari
from cypari2.gen cimport Gen

cdef Gen new_t_POL_from_int_star(int* vals, unsigned long length, long varnum)
