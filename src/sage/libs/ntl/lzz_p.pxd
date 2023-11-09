# distutils: depends = NTL/ZZ.h

from sage.libs.ntl.types cimport zz_p_c

cdef extern from "ntlwrap.h":
    long zz_p_rep "rep"(zz_p_c x)
    long zz_p_isZero "IsZero"(zz_p_c x)
    void zz_p_add "add"(zz_p_c x, zz_p_c a, zz_p_c b)
    void zz_p_sub "sub"(zz_p_c x, zz_p_c a, zz_p_c b)
    void zz_p_negate "NTL::negate"(zz_p_c x, zz_p_c a)
    void zz_p_mul "mul"(zz_p_c x, zz_p_c a, zz_p_c b)
    void zz_p_div "div"(zz_p_c x, zz_p_c a, zz_p_c b)
    void zz_p_inv "inv"(zz_p_c x, zz_p_c a)
    void zz_p_power "NTL::power"(zz_p_c x, zz_p_c a, long e)
    void zz_p_sqr "sqr"(zz_p_c x, zz_p_c a)
    void zz_p_clear "clear"(zz_p_c x)
    void zz_p_set_zero "clear"(zz_p_c x)
    void zz_p_set_one "set"(zz_p_c x)
    void zz_p_swap "swap"(zz_p_c x, zz_p_c y)
