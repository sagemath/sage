# distutils: depends = NTL/ZZ.h

from sage.libs.ntl.ntl_ZZ_decl cimport ZZ_c

cdef extern from "sage/libs/ntl/ntlwrap.cpp":
    ctypedef struct ZZ_pContext_c "struct ZZ_pContext":
        void (*restore)()

    ZZ_pContext_c* ZZ_pContext_new "New<ZZ_pContext>"()
    ZZ_pContext_c* ZZ_pContext_construct "Construct<ZZ_pContext>"(void *mem)
    ZZ_pContext_c* ZZ_pContext_new_ZZ "ZZ_pContext_new"(ZZ_c* p)
    ZZ_pContext_c* ZZ_pContext_construct_ZZ "ZZ_pContext_construct"(void *mem, ZZ_c* p)
    void ZZ_pContext_destruct "Destruct<ZZ_pContext>"(ZZ_pContext_c *mem)
    void ZZ_pContext_delete "Delete<ZZ_pContext>"(ZZ_pContext_c *mem)
