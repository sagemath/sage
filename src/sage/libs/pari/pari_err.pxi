cdef extern from 'interrupt.h':
    int pari_catch_sig_on "sig_on"() except 0
    int pari_catch_sig_str "sig_str"(char *) except 0
    void pari_catch_sig_off "sig_off"()
