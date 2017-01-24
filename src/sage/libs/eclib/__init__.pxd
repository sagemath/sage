# distutils: language = c++
# distutils: extra_compile_args = -DNTL_ALL
# distutils: libraries = ec ntl pari gmp m

# NOTE: eclib includes have specific dependencies and must be included
# in a specific order. That explains the various empty
# "cdef extern from" blocks below.

cdef extern from "eclib/matrix.h":
    ctypedef int scalar

    cdef cppclass mat:
        mat()
        mat(mat m)
        scalar* get_entries()
        scalar sub(long, long)
        long nrows()
        long ncols()
        long rank()

    mat addscalar(mat M, scalar)

from sage.libs.ntl.types cimport ZZ_c
ctypedef ZZ_c bigint

cdef extern from "eclib/interface.h":
    int I2int(bigint)

cdef extern from "eclib/rat.h":
    cdef cppclass rational:
        rational()
        rational(long num, long den)
    cdef long rational_num "num"(rational q)
    cdef long rational_den "den"(rational q)

cdef extern from "eclib/bigrat.h":
    cdef cppclass bigrational:
        pass
    cdef bigint bigrational_num "num"(bigrational q)
    cdef bigint bigrational_den "den"(bigrational q)

cdef extern from "eclib/isogs.h":
    pass

cdef extern from "eclib/curve.h":
    cdef cppclass Curve:
        Curve(bigint aa1, bigint aa2, bigint aa3, bigint aa4, bigint aa6)
        void getai(bigint a1, bigint a2, bigint a3, bigint a4, bigint a6)

    cdef cppclass Curvedata:
        Curvedata(Curve C, int m)
        Curvedata(bigint a1, bigint a2, bigint a3, bigint a4, bigint a6,
                int min_on_init)
        void getai(bigint a1, bigint a2, bigint a3, bigint a4, bigint a6)

    cdef cppclass CurveRed:
        CurveRed(Curvedata CD)

    bigint getconductor(CurveRed CR)

cdef extern from "eclib/descent.h":
    cdef cppclass two_descent:
        two_descent(Curvedata* curve,
                    int verb, int sel,
                    long firstlim, long secondlim,
                    long n_aux, int second_descent)

    cdef cppclass mw:
        mw(Curvedata* curve, int verb, int pp, int maxr)

cdef extern from "eclib/egr.h":
    pass
cdef extern from "eclib/htconst.h":
    pass
cdef extern from "eclib/moddata.h":
    pass
cdef extern from "eclib/symb.h":
    pass

cdef extern from "eclib/homspace.h":
    cdef cppclass homspace:
        long modulus
        int plusflag # sign
        int cuspidal

        homspace(long n, int hp, int hcusp, int verbose)

        int h1dim()
        int h1cuspdim()
        mat heckeop(long p, int dual, int display)
        long h1ncusps()

cdef extern from "eclib/oldforms.h":
    pass

from libcpp.vector cimport vector

cdef extern from "eclib/newforms.h":
    cdef cppclass newforms:
        vector[newform] nflist
        int verbose
        int basisflag
        long maxdepth
        long cuspidal
        long sign

        newforms(long n, int disp)

        void createfromcurve(int sign, CurveRed CR)
        void display()
        rational plus_modular_symbol(rational r)
        rational minus_modular_symbol(rational r)

    cdef cppclass newform:
        newforms* nf
