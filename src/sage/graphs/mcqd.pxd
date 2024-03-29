# sage_setup: distribution = sagemath-mcqd

from libcpp cimport bool

cdef extern from "mcqd.h":
    cdef cppclass Maxclique:
        Maxclique()
        Maxclique(bool **, int n)
        void mcqdyn(int * maxclique, int& size)
