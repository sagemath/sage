# distutils: depends = pari/pari.h
"""
Declarations for types used by PARI

This includes both the C types as well as the PARI types (and a few
macros for dealing with those).

It is important that the functionality in this file does not call any
PARI library functions. The reason is that we want to allow just using
these types (for example, to define a Cython extension type) without
linking to PARI. This file should consist only of typedefs and macros
from PARI's include files.
"""

cdef extern from "sage/libs/pari/parisage.h":
    ctypedef unsigned long ulong "pari_ulong"

    ctypedef long* GEN
    ctypedef char* byteptr
    ctypedef unsigned long pari_sp

    # PARI types: these are actually an enum type, but that doesn't
    # matter for Cython.
    int t_INT, t_REAL, t_INTMOD, t_FRAC, t_FFELT, t_COMPLEX, t_PADIC, \
        t_QUAD, t_POLMOD, t_POL, t_SER, t_RFRAC, t_QFR, t_QFI, t_VEC, \
        t_COL, t_MAT, t_LIST, t_STR, t_VECSMALL, t_CLOSURE, t_ERROR, \
        t_INFINITY

    int BITS_IN_LONG
    int DEFAULTPREC       #  64 bits precision
    int MEDDEFAULTPREC    # 128 bits precision
    int BIGDEFAULTPREC    # 192 bits precision

    long    typ(GEN x)
    long    settyp(GEN x, long s)
    long    isclone(GEN x)
    long    setisclone(GEN x)
    long    unsetisclone(GEN x)
    long    lg(GEN x)
    long    setlg(GEN x, long s)
    long    signe(GEN x)
    long    setsigne(GEN x, long s)
    long    lgefint(GEN x)
    long    setlgefint(GEN x, long s)
    long    expo(GEN x)
    long    setexpo(GEN x, long s)
    long    valp(GEN x)
    long    setvalp(GEN x, long s)
    long    precp(GEN x)
    long    setprecp(GEN x, long s)
    long    varn(GEN x)
    long    setvarn(GEN x, long s)
    long    evaltyp(long x)
    long    evallg(long x)
    long    evalvarn(long x)
    long    evalsigne(long x)
    long    evalprecp(long x)
    long    evalvalp(long x)
    long    evalexpo(long x)
    long    evallgefint(long x)

    # Various structures that we don't use in Sage but which need to be
    # declared, such that Cython understands the declarations of
    # functions using these types.
    struct bb_group
    struct bb_field
    struct bb_ring
    struct bb_algebra
    struct qfr_data
    struct nfmaxord_t
    struct forcomposite_t
    struct forpart_t
    struct forprime_t
    struct forvec_t
    struct entree
    struct gp_context
    struct nfbasic_t
    struct pariFILE
    struct pari_mt
    struct pari_stack
    struct pari_thread
    struct pari_timer
    struct GENbin
    struct hashentry
    struct hashtable
