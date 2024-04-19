# sage_setup: distribution = sagemath-ntl
from sage.libs.ntl.types cimport mat_ZZ_c, ZZ_c, ZZX_c


cdef extern from "ntlwrap.h":
    void mat_ZZ_mul "mul"( mat_ZZ_c x, mat_ZZ_c a, mat_ZZ_c b)
    void mat_ZZ_add "add"( mat_ZZ_c x, mat_ZZ_c a, mat_ZZ_c b)
    void mat_ZZ_sub "sub"( mat_ZZ_c x, mat_ZZ_c a, mat_ZZ_c b)
    void mat_ZZ_power "NTL::power"( mat_ZZ_c x, mat_ZZ_c a, long e)
    void mat_ZZ_CharPoly "CharPoly"(ZZX_c r, mat_ZZ_c m)

    cdef long mat_ZZ_LLL_FP   "LLL_FP"(mat_ZZ_c B, double delta, int deep, int check , int verbose)
    cdef long mat_ZZ_LLL_FP_U "LLL_FP"(mat_ZZ_c B, mat_ZZ_c U, double delta, int deep, int check , int verbose)
    cdef long mat_ZZ_LLL_QP   "LLL_QP"(mat_ZZ_c B, double delta, int deep, int check , int verbose)
    cdef long mat_ZZ_LLL_QP_U "LLL_QP"(mat_ZZ_c B, mat_ZZ_c U, double delta, int deep, int check , int verbose)
    cdef long mat_ZZ_LLL_XD   "LLL_XD"(mat_ZZ_c B, double delta, int deep, int check , int verbose)
    cdef long mat_ZZ_LLL_XD_U "LLL_XD"(mat_ZZ_c B, mat_ZZ_c U, double delta, int deep, int check , int verbose)
    cdef long mat_ZZ_LLL_RR   "LLL_RR"(mat_ZZ_c B, double delta, int deep, int check , int verbose)
    cdef long mat_ZZ_LLL_RR_U "LLL_RR"(mat_ZZ_c B, mat_ZZ_c U, double delta, int deep, int check , int verbose)

    cdef long mat_ZZ_G_LLL_FP   "G_LLL_FP"(mat_ZZ_c B, double delta, int deep, int check , int verbose)
    cdef long mat_ZZ_G_LLL_FP_U "G_LLL_FP"(mat_ZZ_c B, mat_ZZ_c U, double delta, int deep, int check , int verbose)
    cdef long mat_ZZ_G_LLL_QP   "G_LLL_QP"(mat_ZZ_c B, double delta, int deep, int check , int verbose)
    cdef long mat_ZZ_G_LLL_QP_U "G_LLL_QP"(mat_ZZ_c B, mat_ZZ_c U, double delta, int deep, int check , int verbose)
    cdef long mat_ZZ_G_LLL_XD   "G_LLL_XD"(mat_ZZ_c B, double delta, int deep, int check , int verbose)
    cdef long mat_ZZ_G_LLL_XD_U "G_LLL_XD"(mat_ZZ_c B, mat_ZZ_c U, double delta, int deep, int check , int verbose)
    cdef long mat_ZZ_G_LLL_RR   "G_LLL_RR"(mat_ZZ_c B, double delta, int deep, int check , int verbose)
    cdef long mat_ZZ_G_LLL_RR_U "G_LLL_RR"(mat_ZZ_c B, mat_ZZ_c U, double delta, int deep, int check , int verbose)

    cdef long mat_ZZ_BKZ_FP     "BKZ_FP"(mat_ZZ_c B, double delta, long BlockSize, long prune, int check, long verbose)
    cdef long mat_ZZ_BKZ_FP_U   "BKZ_FP"(mat_ZZ_c B, mat_ZZ_c U, double delta, long BlockSize, long prune, int check, long verbose)
    cdef long mat_ZZ_BKZ_XD     "BKZ_XD"(mat_ZZ_c B, double delta, long BlockSize, long prune, int check, long verbose)
    cdef long mat_ZZ_BKZ_XD_U   "BKZ_XD"(mat_ZZ_c B, mat_ZZ_c U, double delta, long BlockSize, long prune, int check, long verbose)
    cdef long mat_ZZ_BKZ_QP     "BKZ_QP"(mat_ZZ_c B, double delta, long BlockSize, long prune, int check, long verbose)
    cdef long mat_ZZ_BKZ_QP_U   "BKZ_QP"(mat_ZZ_c B, mat_ZZ_c U, double delta, long BlockSize, long prune, int check, long verbose)
    cdef long mat_ZZ_BKZ_QP1    "BKZ_QP1"(mat_ZZ_c B, double delta, long BlockSize, long prune, int check, long verbose)
    cdef long mat_ZZ_BKZ_QP1_U  "BKZ_QP1"(mat_ZZ_c B, mat_ZZ_c U, double delta, long BlockSize, long prune, int check, long verbose)
    cdef long mat_ZZ_BKZ_RR     "BKZ_RR"(mat_ZZ_c B, double delta, long BlockSize, long prune, int check, long verbose)
    cdef long mat_ZZ_BKZ_RR_U   "BKZ_RR"(mat_ZZ_c B, mat_ZZ_c U, double delta, long BlockSize, long prune, int check, long verbose)

    cdef long mat_ZZ_G_BKZ_FP     "G_BKZ_FP"(mat_ZZ_c B, double delta, long BlockSize, long prune, int check, long verbose)
    cdef long mat_ZZ_G_BKZ_FP_U   "G_BKZ_FP"(mat_ZZ_c B, mat_ZZ_c U, double delta, long BlockSize, long prune, int check, long verbose)
    cdef long mat_ZZ_G_BKZ_XD     "G_BKZ_XD"(mat_ZZ_c B, double delta, long BlockSize, long prune, int check, long verbose)
    cdef long mat_ZZ_G_BKZ_XD_U   "G_BKZ_XD"(mat_ZZ_c B, mat_ZZ_c U, double delta, long BlockSize, long prune, int check, long verbose)
    cdef long mat_ZZ_G_BKZ_QP     "G_BKZ_QP"(mat_ZZ_c B, double delta, long BlockSize, long prune, int check, long verbose)
    cdef long mat_ZZ_G_BKZ_QP_U   "G_BKZ_QP"(mat_ZZ_c B, mat_ZZ_c U, double delta, long BlockSize, long prune, int check, long verbose)
    cdef long mat_ZZ_G_BKZ_QP1    "G_BKZ_QP1"(mat_ZZ_c B, double delta, long BlockSize, long prune, int check, long verbose)
    cdef long mat_ZZ_G_BKZ_QP1_U  "G_BKZ_QP1"(mat_ZZ_c B, mat_ZZ_c U, double delta, long BlockSize, long prune, int check, long verbose)
    cdef long mat_ZZ_G_BKZ_RR     "G_BKZ_RR"(mat_ZZ_c B, double delta, long BlockSize, long prune, int check, long verbose)
    cdef long mat_ZZ_G_BKZ_RR_U   "G_BKZ_RR"(mat_ZZ_c B, mat_ZZ_c U, double delta, long BlockSize, long prune, int check, long verbose)


cdef extern from "ntlwrap_impl.h":
    void mat_ZZ_setitem(mat_ZZ_c* x, int i, int j, ZZ_c* z)
    ZZ_c* mat_ZZ_getitem(mat_ZZ_c* x, int i, int j)
    ZZ_c* mat_ZZ_determinant(mat_ZZ_c* x, long deterministic)
    mat_ZZ_c* mat_ZZ_HNF(mat_ZZ_c* A, ZZ_c* D)
    cdef long mat_ZZ_LLL(ZZ_c **det, mat_ZZ_c *x, long a, long b, long verbose)
    cdef long mat_ZZ_LLL_U(ZZ_c **det, mat_ZZ_c *x, mat_ZZ_c *U, long a, long b, long verbose)
