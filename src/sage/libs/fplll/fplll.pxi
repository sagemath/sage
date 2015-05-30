#
# general include
#

from libcpp.vector cimport vector

#
# integers
#

cdef extern from "fplll/nr.h" namespace "fplll":
    cdef cppclass Z_NR[T]:
        T& getData()
        void set(mpz_t d)

#
# matrices over the integers
#

cdef extern from "fplll/matrix.h" namespace "fplll":
    cdef cppclass MatrixRow[T]:
        Z_NR[T]& operator[](int i)

    cdef cppclass ZZ_mat[T]:

        ZZ_mat()
        ZZ_mat(int r, int c)

        int getRows()
        int getCols()

        T& operator()(int i, int j)
        MatrixRow[T] operator[](int i)

        void gen_intrel(int bits)
        void gen_simdioph(int bits, int bits2)
        void gen_uniform(int bits)
        void gen_ntrulike(int bits, int q)
        void gen_ntrulike2(int bits, int q)
        void gen_ajtai(double alpha)


cdef extern from "fplll/defs.h" namespace "fplll":

    cdef enum RedStatus:
        RED_SUCCESS
        RED_GSO_FAILURE
        RED_BABAI_FAILURE
        RED_LLL_FAILURE
        RED_ENUM_FAILURE
        RED_BKZ_FAILURE
        RED_BKZ_TIME_LIMIT
        RED_BKZ_LOOPS_LIMIT
        RED_STATUS_MAX

    cdef enum LLLFlags:
        LLL_VERBOSE
        LLL_EARLY_RED
        LLL_SIEGEL
        LLL_DEFAULT

    cdef enum BKZFlags:
        BKZ_DEFAULT
        BKZ_VERBOSE
        BKZ_NO_LLL
        BKZ_MAX_LOOPS
        BKZ_MAX_TIME
        BKZ_BOUNDED_LLL
        BKZ_AUTO_ABORT

    cdef enum LLLMethod:
        LM_WRAPPER
        LM_PROVED
        LM_HEURISTIC
        LM_FAST

    cdef enum IntType:
        ZT_MPZ
        ZT_LONG
        ZT_DOUBLE

    cdef enum FloatType:
        FT_DEFAULT
        FT_DOUBLE
        FT_LONG_DOUBLE
        FT_DPE
        FT_MPFR

    cdef enum SVPMethod:
        SVPM_FAST
        SVPM_PROVED

    cdef double LLL_DEF_DELTA
    cdef double LLL_DEF_ETA

cdef extern from "fplll/fplll.h" namespace "fplll":

    int lllReduction(ZZ_mat[mpz_t] b, double delta, double eta,
                     LLLMethod method, FloatType floatType,
                     int precision, int flags)
    int lllReduction(ZZ_mat[mpz_t] b, ZZ_mat[mpz_t] u,
                     double delta, double eta,
                     LLLMethod method, FloatType floatType,
                     int precision, int flags)

    cdef cppclass BKZParam:
         ZZ_mat[mpz_t]* b
         ZZ_mat[mpz_t]* u
         int blockSize
         double delta
         FloatType floatType
         int precision
         int flags
         int maxLoops
         double maxTime
         vector[double] pruning

    int bkzReduction(BKZParam &param)
    int bkzReduction(ZZ_mat[mpz_t] b, int blockSize)

    int hkzReduction(ZZ_mat[mpz_t] b)
    int shortestVector(ZZ_mat[mpz_t] b,
                       vector[Z_NR[mpz_t]] &solCoord,
                       SVPMethod method)
    const char* getRedStatusStr (int status)

cdef extern from "fplll/util.h" namespace "fplll":
    void vectMatrixProduct(vector[Z_NR[mpz_t]] &result,
                           vector[Z_NR[mpz_t]] &x,
                           const ZZ_mat[mpz_t] &m)


#
# fpLLL 3.x interface
#

#
# fastest LLL
#

cdef  extern from "fplll/fplllv31.h":
    ctypedef struct fast_double "fplll::fast<mpz_t,double>":
      int (*LLL)()
      ZZ_mat* (*GetBase)()

    fast_double *fast_double_new "new fplll::fast<mpz_t,double>"(ZZ_mat[mpz_t] *B,int precision, double eta, double delta)
    void fast_double_delete "delete "(fast_double *mem)

#
# fastest LLL with early reduction
#

cdef  extern from "fplll/fplllv31.h":
    ctypedef struct fast_early_red_double "fplll::fast_early_red<mpz_t,double>":
      int (*LLL)()
      ZZ_mat* (*GetBase)()

    fast_early_red_double *fast_early_red_double_new "new fplll::fast_early_red<mpz_t,double>"(ZZ_mat[mpz_t] *B,int precision, double eta, double delta)
    void fast_early_red_double_delete "delete "(fast_early_red_double *mem)

#
# heuristic
#

cdef  extern from "fplll/fplllv31.h":
    ctypedef struct heuristic_mpfr "fplll::heuristic<mpz_t,mpfr_t>":
      int (*LLL)()
      ZZ_mat* (*GetBase)()

    heuristic_mpfr *heuristic_mpfr_new "new fplll::heuristic<mpz_t,mpfr_t>"(ZZ_mat[mpz_t] *B,int precision, double eta, double delta)
    void heuristic_mpfr_delete "delete "(heuristic_mpfr *mem)

    ctypedef struct heuristic_dpe "fplll::heuristic<mpz_t,dpe_t>":
      int (*LLL)()
      ZZ_mat* (*GetBase)()

    heuristic_dpe *heuristic_dpe_new "new fplll::heuristic<mpz_t,dpe_t>"(ZZ_mat[mpz_t] *B,int precision, double eta, double delta)
    void heuristic_dpe_delete "delete "(heuristic_dpe *mem)

    ctypedef struct heuristic_double "fplll::heuristic<mpz_t,double>":
      int (*LLL)()
      ZZ_mat* (*GetBase)()

    heuristic_double *heuristic_double_new "new fplll::heuristic<mpz_t,double>"(ZZ_mat[mpz_t] *B,int precision, double eta, double delta)
    void heuristic_double_delete "delete "(heuristic_double *mem)

#
# heuristic with early reduction
#

cdef  extern from "fplll/fplllv31.h":
    ctypedef struct heuristic_early_red_mpfr "fplll::heuristic_early_red<mpz_t,mpfr_t>":
      int (*LLL)()
      ZZ_mat* (*GetBase)()

    heuristic_early_red_mpfr *heuristic_early_red_mpfr_new "new fplll::heuristic_early_red<mpz_t,mpfr_t>"(ZZ_mat[mpz_t] *B,int precision, double eta, double delta)
    void heuristic_early_red_mpfr_delete "delete "(heuristic_early_red_mpfr *mem)

    ctypedef struct heuristic_early_red_dpe "fplll::heuristic_early_red<mpz_t,dpe_t>":
      int (*LLL)()
      ZZ_mat* (*GetBase)()

    heuristic_early_red_dpe *heuristic_early_red_dpe_new "new fplll::heuristic_early_red<mpz_t,dpe_t>"(ZZ_mat[mpz_t] *B,int precision, double eta, double delta)
    void heuristic_early_red_dpe_delete "delete "(heuristic_early_red_dpe *mem)

    ctypedef struct heuristic_early_red_double "fplll::heuristic_early_red<mpz_t,double>":
      int (*LLL)()
      ZZ_mat* (*GetBase)()

    heuristic_early_red_double *heuristic_early_red_double_new "new fplll::heuristic_early_red<mpz_t,double>"(ZZ_mat[mpz_t] *B,int precision, double eta, double delta)
    void heuristic_early_red_double_delete "delete "(heuristic_early_red_double *mem)



#
# provable LLL
#

cdef  extern from "fplll/fplllv31.h":
    ctypedef struct proved_mpfr "fplll::proved<mpz_t,mpfr_t>":
      int (*LLL)()
      ZZ_mat* (*GetBase)()

    proved_mpfr *proved_mpfr_new "new fplll::proved<mpz_t,mpfr_t>"(
        ZZ_mat[mpz_t] *B,int precision, double eta, double delta)
    void proved_mpfr_delete "delete "(proved_mpfr *mem)

    ctypedef struct proved_dpe "fplll::proved<mpz_t,dpe_t>":
      int (*LLL)()
      ZZ_mat* (*GetBase)()

    proved_dpe *proved_dpe_new "new fplll::proved<mpz_t,dpe_t>"(
        ZZ_mat[mpz_t] *B,int precision, double eta, double delta)
    void proved_dpe_delete "delete "(proved_dpe *mem)

    ctypedef struct proved_double "fplll::proved<mpz_t,double>":
      int (*LLL)()
      ZZ_mat* (*GetBase)()

    proved_double *proved_double_new "new fplll::proved<mpz_t,double>"(
        ZZ_mat[mpz_t] *B,int precision, double eta, double delta)
    void proved_double_delete "delete "(proved_double *mem)

#
# wrapper code which chooses a LLL sequence automatically
#

cdef  extern from "fplll/fplllv31.h":
    ctypedef struct wrapper "fplll::wrapper":
      int (*LLL)()
      ZZ_mat* (*GetBase)()

    wrapper *wrapper_new "new fplll::wrapper"(
        ZZ_mat[mpz_t] *B, int precision, double eta, double delta)
    void wrapper_delete "delete "(wrapper *mem)
