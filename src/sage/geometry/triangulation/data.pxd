cdef extern from "data.h":
    cdef cppclass compact_simplices():
        void push_back(int encoded_simplex)
