
from sage.rings.integer cimport Integer
from sage.rings.padics.pow_computer cimport PowComputer_class

cdef long get_ordp(x, PowComputer_class prime_pow) except? -10000
cdef long get_preccap(x, PowComputer_class prime_pow) except? -10000
cdef long comb_prec(iprec, long prec) except? -10000
cdef int _process_args_and_kwds(long *aprec, long *rprec, args, kwds, bint absolute, PowComputer_class prime_pow) except -1

