# sage_setup: distribution = sagemath-ntl
from sage.libs.ntl.types cimport vec_GF2_c, GF2_c


cdef extern from "ntlwrap.h":
    int vec_GF2_IsZero "IsZero"(vec_GF2_c x)

    void vec_GF2_append_GF2 "append"(vec_GF2_c v, GF2_c a)
    void vec_GF2_append_vec "append"(vec_GF2_c v,  vec_GF2_c a)
    void vec_GF2_VectorCopy "VectorCopy"(vec_GF2_c x, vec_GF2_c a, long n)

    void vec_GF2_clear "clear"(vec_GF2_c x)
    void vec_GF2_shift "shift"(vec_GF2_c x,  vec_GF2_c a, long n)
    void vec_GF2_reverse "reverse"(vec_GF2_c x,  vec_GF2_c a)
    long vec_GF2_weight "weight"(vec_GF2_c a)

    void vec_GF2_random "random"(vec_GF2_c x, long n)

    void vec_GF2_add "add"( vec_GF2_c x, vec_GF2_c a, vec_GF2_c b)
    void vec_GF2_sub "sub"( vec_GF2_c x, vec_GF2_c a, vec_GF2_c b)
    void vec_GF2_mul "mul"( vec_GF2_c x, vec_GF2_c a, vec_GF2_c b)
    void vec_GF2_negate "NTL::negate"(vec_GF2_c x, vec_GF2_c a)
    void vec_GF2_InnerProduct "InnerProduct"(GF2_c x, vec_GF2_c a, vec_GF2_c b)

    void vec_GF2_conv_long "conv"(vec_GF2_c x, long i)
    long vec_GF2_conv_to_long "rep"(vec_GF2_c x)
