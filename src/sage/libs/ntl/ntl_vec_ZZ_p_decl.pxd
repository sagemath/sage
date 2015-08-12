# distutils: depends = NTL/ZZ.h

cdef extern from "sage/libs/ntl/ntlwrap.h":
    #### vec_ZZ_p_c
    ctypedef struct vec_ZZ_p_c "struct vec_ZZ_p":
        pass

