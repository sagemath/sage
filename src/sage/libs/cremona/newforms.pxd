include "sage/libs/ntl/decl.pxi"
include "defs.pxi"

cdef extern from "eclib/moddata.h":
    pass
cdef extern from "eclib/symb.h":
    pass
cdef extern from "eclib/cusp.h":
    pass
cdef extern from "eclib/xsplit.h":
    pass
cdef extern from "eclib/method.h":
    pass
cdef extern from "eclib/oldforms.h":
    pass
cdef extern from "eclib/homspace.h":
    pass
cdef extern from "eclib/cperiods.h":
    pass

cdef extern from "eclib/newforms.h":
    ctypedef struct nfvec "std::vector<newform>":
        pass

    # From newforms.h
    ctypedef struct newforms "newforms":
        #attributes
        nfvec nflist
        int verbose
        int basisflag
        long maxdepth
        long cuspidal
        long sign

        #member functions
        void (* createfromcurve)(int sign, CurveRed CR)
        void (* display)()
        rational (* plus_modular_symbol)(rational r)
        rational (* minus_modular_symbol)(rational r)

    ctypedef struct newform "newform":
        #attributes
        newforms *nf

    # Constructors
    newforms *new_newforms "new newforms" (long n, int disp)
    # General C++ stuff
    void delete_newforms "delete "(newforms* nfs)


cdef class ECModularSymbol:
    cdef newforms* nfs
    cdef int n
    cdef object _E
    cdef int sign
