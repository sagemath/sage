cdef extern from "data.h":
    cdef cppclass compact_simplices(object):
        void push_back(int encoded_simplex)
