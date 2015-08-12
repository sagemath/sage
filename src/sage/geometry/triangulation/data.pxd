cdef extern from "sage/geometry/triangulation/data.h":
    cdef cppclass compact_simplices(object):
        void push_back(int encoded_simplex)
