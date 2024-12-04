# distutils: language = c++
# distutils: libraries = ec NTL_LIBRARIES pari gmp m
# distutils: extra_compile_args = NTL_CFLAGS
# distutils: include_dirs = NTL_INCDIR
# distutils: library_dirs = NTL_LIBDIR
# distutils: extra_link_args = NTL_LIBEXTRA


from libcpp.map cimport map
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from sage.libs.ntl.types cimport ZZ_c


# NOTE: eclib used to have specific dependencies, so that they had to
# be included in a specific order. Although this is no longer the
# case, we start by listing all relevant include files in the correct
# order.

cdef extern from "eclib/vector.h": pass
cdef extern from "eclib/xmod.h": pass
cdef extern from "eclib/svector.h": pass
cdef extern from "eclib/matrix.h": pass
cdef extern from "eclib/smatrix.h": pass
cdef extern from "eclib/interface.h": pass
cdef extern from "eclib/rat.h": pass
cdef extern from "eclib/bigrat.h": pass
cdef extern from "eclib/isogs.h": pass
cdef extern from "eclib/curve.h": pass
cdef extern from "eclib/descent.h": pass
cdef extern from "eclib/egr.h": pass
cdef extern from "eclib/htconst.h": pass
cdef extern from "eclib/moddata.h": pass
cdef extern from "eclib/symb.h": pass
cdef extern from "eclib/homspace.h": pass
cdef extern from "eclib/oldforms.h": pass
cdef extern from "eclib/newforms.h": pass


# Now the actual declarations, where the order no longer matters
cdef extern from "eclib/vector.h":
    ctypedef int scalar

    cdef cppclass vec:
        scalar operator[](long)

cdef extern from "eclib/svector.h":
    cdef cppclass svec:
        vec as_vec()
        scalar elem(int) const
        map[int,scalar].iterator begin()
        map[int,scalar].iterator end()

cdef extern from "eclib/matrix.h":
    cdef cppclass mat:
        mat()
        mat(mat m)
        scalar sub(long, long)
        long nrows()
        long ncols()
        long rank()

    mat addscalar(mat M, scalar)

cdef extern from "eclib/smatrix.h":
    cdef cppclass smat:
        smat()
        smat(smat m)
        scalar sub(long, long)
        long nrows()
        long ncols()
        long rank()
        scalar elem(int, int) const
        svec row(int) const

cdef extern from "eclib/interface.h":
    ctypedef ZZ_c bigint
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

cdef extern from "eclib/curve.h":
    cdef cppclass Curve:
        Curve()
        Curve(bigint aa1, bigint aa2, bigint aa3, bigint aa4, bigint aa6)
        void getai(bigint a1, bigint a2, bigint a3, bigint a4, bigint a6)

    cdef cppclass Curvedata:
        Curvedata()
        Curvedata(Curve C, int m)
        Curvedata(bigint a1, bigint a2, bigint a3, bigint a4, bigint a6,
                int min_on_init)
        void getai(bigint a1, bigint a2, bigint a3, bigint a4, bigint a6)

    cdef cppclass CurveRed:
        CurveRed()
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

cdef extern from "eclib/homspace.h":
    cdef cppclass homspace:
        long modulus
        int plusflag # sign
        int cuspidal

        homspace(long n, int hp, int hcusp, int verbose)

        int h1dim()
        int h1cuspdim()
        mat heckeop(long p, int dual, int display)
        vec heckeop_col(long p, int j, int display)
        smat s_heckeop(long p, int dual, int display)
        long h1ncusps()

cdef extern from "eclib/newforms.h":
    cdef cppclass newforms:
        vector[newform] nflist
        int verbose
        int basisflag
        long maxdepth
        long cuspidal
        long sign

        newforms(long n, int disp)

        void createfromcurve(int sign, CurveRed CR, int nap)
        void display()
        # Here i is the index of the relevant newform in the space,
        # which for us will always be 0:
        rational plus_modular_symbol(rational r, int i, int base_at_infinity)
        rational minus_modular_symbol(rational r, int i, int base_at_infinity)
        pair[rational,rational] full_modular_symbol(rational r, int i, int base_at_infinity)

    cdef cppclass newform:
        newforms* nf
