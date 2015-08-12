# distutils: depends = NTL/ZZ.h

cdef extern from "sage/libs/ntl/ntlwrap.h":
    ctypedef struct vec_ZZ_pE_c "struct vec_ZZ_pE":
        pass
