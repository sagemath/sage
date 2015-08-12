# distutils: depends = NTL/ZZ.h

from sage.libs.ntl.ntl_ZZ_pX_decl cimport ZZ_pX_c

cdef extern from "sage/libs/ntl/ntlwrap.cpp":
    ctypedef struct ZZ_pEContext_c "struct ZZ_pEContext":
        void (*restore)()

    ZZ_pEContext_c* ZZ_pEContext_new "New<ZZ_pEContext>"()
    ZZ_pEContext_c* ZZ_pEContext_construct "Construct<ZZ_pEContext>"(void *mem)
    ZZ_pEContext_c* ZZ_pEContext_new_ZZ_pX "ZZ_pEContext_new"(ZZ_pX_c* p)
    ZZ_pEContext_c* ZZ_pEContext_construct_ZZ_pX "ZZ_pEContext_construct"(void *mem, ZZ_pX_c* p)
    void ZZ_pEContext_destruct "Destruct<ZZ_pEContext>"(ZZ_pEContext_c *mem)
    void ZZ_pEContext_delete "Delete<ZZ_pEContext>"(ZZ_pEContext_c *mem)

    void ZZ_pEContext_restore(ZZ_pEContext_c *c)

