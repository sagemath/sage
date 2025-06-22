
cdef extern from "sage/graphs/path_enumeration.h":
    cdef struct PathCandidate:
        double cost
        int path_idx
        int dev_idx
        bint is_simple
