# sage_setup: distribution = sagemath-ntl
from sage.libs.flint.types cimport fmpz, fmpz_t, fmpz_poly_t, flint_rand_t

ctypedef fmpz_t cdigit
ctypedef fmpz* cdigit_ptr
ctypedef fmpz_poly_t celement
ctypedef flint_rand_t randgen


include "relaxed_template_header.pxi"


cdef class pAdicRelaxedElement(RelaxedElement):
    pass

cdef class pAdicRelaxedElement_zero(RelaxedElement_zero):
    pass

cdef class pAdicRelaxedElement_one(RelaxedElement_one):
    pass

cdef class pAdicRelaxedElement_bound(RelaxedElement_bound):
    pass

cdef class pAdicRelaxedElement_value(RelaxedElement_value):
    pass

cdef class pAdicRelaxedElement_random(RelaxedElement_random):
    pass

cdef class pAdicRelaxedElement_slice(RelaxedElement_slice):
    pass

cdef class pAdicRelaxedElement_add(RelaxedElement_add):
    pass

cdef class pAdicRelaxedElement_sub(RelaxedElement_sub):
    pass

cdef class pAdicRelaxedElement_mul(RelaxedElement_mul):
    pass

cdef class pAdicRelaxedElement_muldigit(RelaxedElement_muldigit):
    pass

cdef class pAdicRelaxedElement_div(RelaxedElement_div):
    pass

cdef class pAdicRelaxedElement_sqrt(RelaxedElement_sqrt):
    pass

cdef class pAdicRelaxedElement_teichmuller(RelaxedElement_teichmuller):
    pass

cdef class pAdicRelaxedElement_unknown(RelaxedElement_unknown):
    pass
