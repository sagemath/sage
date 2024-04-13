# sage_setup: distribution = sagemath-ntl
# distutils: depends = NTL/ZZ.h

cdef extern from "ntlwrap.h":
    long NTL_OVFBND, NTL_SP_BOUND

    cdef cppclass ZZ_c "ZZ":
        bint operator==(ZZ_c)
        bint operator!=(ZZ_c)
        bint operator<(ZZ_c)
        bint operator<=(ZZ_c)
        bint operator>(ZZ_c)
        bint operator>=(ZZ_c)

    cdef cppclass zz_p_c "zz_p":
        zz_p_c operator=(long)
        bint operator==(zz_p_c)
        bint operator!=(zz_p_c)
        void *rep

    cdef cppclass ZZ_p_c "ZZ_p":
        bint operator==(ZZ_p_c)
        bint operator!=(ZZ_p_c)

    cdef cppclass ZZ_pE_c "ZZ_pE":
        bint operator==(ZZ_pE_c)
        bint operator!=(ZZ_pE_c)

    cdef cppclass GF2_c "GF2":
        bint operator==(GF2_c)
        bint operator!=(GF2_c)

    cdef cppclass GF2E_c "GF2E":
        bint operator==(GF2E_c)
        bint operator!=(GF2E_c)

    cdef cppclass vec_ZZ_c "vec_ZZ":
        ZZ_c RawGet(long i)
        ZZ_c *elts()
        long length()

    cdef cppclass vec_ZZ_p_c "vec_ZZ_p":
        pass

    cdef cppclass vec_ZZ_pE_c "vec_ZZ_pE":
        pass

    cdef cppclass vec_GF2_c "vec_GF2":
        void SetLength(long n)
        void SetMaxLength(long n)
        long length()
        long MaxLength()
        long allocated()

        GF2_c (*get)(long i)
        void  (*put_GF2 "put")(long i, GF2_c a)
        void  (*put_long "put")(long i, long a)

    cdef cppclass vec_GF2E_c "vec_GF2E":
        pass

    cdef cppclass ZZX_c "ZZX":
        bint operator==(ZZX_c)
        bint operator!=(ZZX_c)
        vec_ZZ_c rep

    cdef cppclass zz_pX_c "zz_pX":
        bint operator==(zz_pX_c)
        bint operator!=(zz_pX_c)
        void *rep
        void SetMaxLength(long n)

    cdef cppclass zz_pX_Modulus_c "zz_pXModulus":
        zz_pX_c val()

    cdef cppclass ZZ_pX_c "ZZ_pX":
        ZZ_pX_c()
        bint operator==(ZZ_pX_c)
        bint operator!=(ZZ_pX_c)
        void *rep
        void SetMaxLength(long n)

    cdef cppclass ZZ_pX_Modulus_c "ZZ_pXModulus":
        ZZ_pX_c val()

    cdef cppclass ZZ_pX_Multiplier_c "ZZ_pXMultiplier":
        ZZ_pX_c val()

    cdef cppclass ZZ_pEX_c "ZZ_pEX":
        bint operator==(ZZ_pEX_c)
        bint operator!=(ZZ_pEX_c)
        void *rep
        void (* SetMaxLength)(long n)

    cdef cppclass ZZ_pEX_Modulus_c "ZZ_pEXModulus":
        ZZ_pEX_c val()

    cdef cppclass GF2X_c "GF2X":
        bint operator==(GF2X_c)
        bint operator!=(GF2X_c)

    cdef cppclass GF2XModulus_c "GF2XModulus":
        pass

    cdef cppclass GF2EX_c "GF2EX":
        bint operator==(GF2EX_c)
        bint operator!=(GF2EX_c)

    cdef cppclass mat_ZZ_c "mat_ZZ":
        long NumRows()
        long NumCols()
        void SetDims(long, long)

    cdef cppclass mat_ZZ_p_c "mat_ZZ_p":
        bint operator==(mat_ZZ_p_c)
        bint operator!=(mat_ZZ_p_c)
        void SetDims(long nrows, long ncols)
        long NumRows()
        long NumCols()
        ZZ_p_c get(long i, long j)
        void put(long i, long j, ZZ_p_c &a)

    cdef cppclass mat_GF2_c "mat_GF2":
        bint operator==(mat_GF2_c)
        bint operator!=(mat_GF2_c)
        void SetDims(long nrows, long ncols)
        long NumRows()
        long NumCols()
        GF2_c (*get "operator()") (long i, long j)

    cdef cppclass mat_GF2E_c "mat_GF2E":
        bint operator==(mat_GF2E_c)
        bint operator!=(mat_GF2E_c)
        void SetDims(long nrows, long ncols)
        long NumRows()
        long NumCols()
        GF2E_c (*get "operator()") (long i, long j)

    cdef cppclass zz_pContext_c "zz_pContext":
        zz_pContext_c()
        zz_pContext_c(long)
        void restore()

    cdef cppclass ZZ_pContext_c "ZZ_pContext":
        ZZ_pContext_c()
        ZZ_pContext_c(ZZ_c)
        void restore()

    cdef cppclass ZZ_pEContext_c "ZZ_pEContext":
        ZZ_pEContext_c()
        ZZ_pEContext_c(ZZ_pX_c)
        void restore()

    cdef cppclass GF2EContext_c "GF2EContext":
        GF2EContext_c()
        GF2EContext_c(GF2X_c)
        void restore()
