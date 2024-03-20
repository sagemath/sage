# sage_setup: distribution = sagemath-ntl
cdef extern from "NTL/tools.h" namespace "NTL":
    void (*ErrorMsgCallback)(const char *) except *
