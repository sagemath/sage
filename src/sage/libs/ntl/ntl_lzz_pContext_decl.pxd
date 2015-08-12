# distutils: depends = NTL/ZZ.h

cdef extern from "sage/libs/ntl/ntlwrap.cpp":
    ctypedef struct zz_pContext_c "struct zz_pContext":
        pass

    zz_pContext_c* zz_pContext_new "New<zz_pContext>"()
    zz_pContext_c* zz_pContext_construct "Construct<zz_pContext>"(void *mem)
    zz_pContext_c* zz_pContext_new_long "zz_pContext_new"(long p)
    zz_pContext_c* zz_pContext_construct_long "zz_pContext_construct"(void *mem, long p)
    void zz_pContext_destruct "Destruct<zz_pContext>"(zz_pContext_c *mem)
    void zz_pContext_delete "Delete<zz_pContext>"(zz_pContext_c *mem)

    void zz_pContext_restore(zz_pContext_c *c)

