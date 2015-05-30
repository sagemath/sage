cdef extern from "NTL/tools.h" namespace "NTL":
    void SetErrorCallbackFunction(void (*func)(const char *s, void *context), void *context)
