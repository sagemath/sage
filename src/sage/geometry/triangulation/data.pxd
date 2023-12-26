# sage_setup: distribution = sagemath-polyhedra
cdef extern from "data.h":
    cdef cppclass compact_simplices():
        void push_back(int encoded_simplex)
