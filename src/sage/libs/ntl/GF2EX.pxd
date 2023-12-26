# sage_setup: distribution = sagemath-ntl
from sage.libs.ntl.types cimport GF2EX_c


cdef extern from "ntlwrap.h":
    void GF2EX_add "add"( GF2EX_c x, GF2EX_c a, GF2EX_c b)
    void GF2EX_sub "sub"( GF2EX_c x, GF2EX_c a, GF2EX_c b)
    void GF2EX_mul "mul"( GF2EX_c x, GF2EX_c a, GF2EX_c b)
    void GF2EX_negate "NTL::negate"(GF2EX_c x, GF2EX_c a)
    void GF2EX_power "NTL::power"(GF2EX_c t, GF2EX_c x, long e)
    int GF2EX_IsOne "IsOne"(GF2EX_c x)
    int GF2EX_IsZero "IsZero"(GF2EX_c x)
