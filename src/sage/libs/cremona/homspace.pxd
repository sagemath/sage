cdef extern from "eclib/moddata.h":
    pass
cdef extern from "eclib/symb.h":
    pass
cdef extern from "eclib/cusp.h":
    pass


cdef extern from "eclib/homspace.h":
    # From mat.h
    ctypedef struct mat "mat":
        pass

    # From homspace.h
    ctypedef struct homspace "homspace":
        # attributes
        long modulus
        int plusflag # sign
        int cuspidal

        # member functions
        int (* h1dim)()
        int (* h1cuspdim)()
        mat (* heckeop)(long p, int dual, int display)
        long (* h1ncusps)()

    # Constructors
    homspace *new_homspace "new homspace" (long n, int hp, int hcusp, int verbose)

    # General C++ stuff
    void delete_homspace "delete "(homspace* H)

cdef class ModularSymbols:
    cdef homspace* H
##     cpdef long level(self)
##     cpdef int dimension(self)
##     cpdef int sign(self)
##     cpdef bint is_cuspidal(self)

