# distutils: include_dirs = SINGULAR_INCDIR
# distutils: extra_compile_args = SINGULAR_CFLAGS
# distutils: libraries = SINGULAR_LIBRARIES
# distutils: library_dirs = SINGULAR_LIBDIR
# distutils: language = c++
# distutils: extra_compile_args = -std=c++11
"""
Declarations of Singular's C/C++ Functions

.. NOTE::

    Our ``ring``, ``poly``,... types are not the Singular ``ring``,
    ``poly``,... types. They are dereferences. So a Singular ``ring`` is
    a pointer to a ``ring`` from Sage.

AUTHOR:

- Martin Albrecht (2009-07): initial implementation
"""
#*****************************************************************************
#       Copyright (C) 2009 Martin Albrecht <malb@informatik.uni-bremen.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.libs.gmp.types cimport mpz_t, mpz_ptr

cdef extern from "factory/factory.h":

    #
    # CF OPTIONS
    #

    void On( int )
    void Off( int )
    int isOn( int )

    cdef int SW_USE_CHINREM_GCD
    cdef int SW_USE_EZGCD
    cdef int SW_USE_NTL
    cdef int SW_USE_NTL_GCD_0
    cdef int SW_USE_NTL_GCD_P
    cdef int SW_USE_NTL_SORT

cdef extern from "singular/Singular/libsingular.h":
    """
    // compatibility for singular 4.3.2p10 and before
    #if SINGULAR_VERSION <= 4330
    #define ringorder_ip ringorder_rp
    #define BIGINTVEC_CMD INTVEC_CMD
    #endif
    """

    #
    # OPTIONS
    #

    cdef unsigned int singular_options "si_opt_1"           # previously 'test'
    cdef unsigned int singular_verbose_options "si_opt_2"   # previously 'verbose'

    # actual options
    cdef int OPT_PROT
    cdef int OPT_REDSB
    cdef int OPT_NOT_BUCKETS
    cdef int OPT_NOT_SUGAR
    cdef int OPT_INTERRUPT
    cdef int OPT_SUGARCRIT
    cdef int OPT_DEBUG
    cdef int OPT_REDTHROUGH
    cdef int OPT_RETURN_SB
    cdef int OPT_FASTHC
    cdef int OPT_OLDSTD
    cdef int OPT_KEEPVARS
    cdef int OPT_STAIRCASEBOUND
    cdef int OPT_MULTBOUND
    cdef int OPT_DEGBOUND
    cdef int OPT_REDTAIL
    cdef int OPT_INTSTRATEGY
    cdef int OPT_INFREDTAIL
    cdef int OPT_SB_1
    cdef int OPT_NOTREGULARITY
    cdef int OPT_WEIGHTM

    # getter/setter functions
    int Sy_bit(int)
    int Sy_inset(int x,int s)
    int BTEST1(int)
    int BVERBOSE(int)

    # ideal flags
    cdef int  FLAG_STD
    cdef int  FLAG_TWOSTD

    #
    # STRUCTS
    #

    # rational numbers
    ctypedef struct number "snumber":
        mpz_t z
        mpz_t n
        int s

    # See singular/libpolys/coeffs/coeffs.h for documentation
    cdef enum n_coeffType:
        n_unknown
        n_Zp
        n_Q
        n_R
        n_GF
        n_long_R
        n_algExt
        n_transExt
        n_long_C
        n_Z
        n_Zn
        n_Znm
        n_Z2m
        n_CF

    ctypedef struct ring "ip_sring"
    ctypedef struct AlgExtInfo

    ctypedef struct n_Procs_s:

        number* cfDiv(number *, number *, const n_Procs_s* r)
        number* cfAdd(number *, number *, const n_Procs_s* r)  # algebraic number addition
        number* cfSub(number *, number *, const n_Procs_s* r)
        number* cfMult(number *, number *, const n_Procs_s* r)  # algebraic number multiplication

        number*  (*cfInit)(int i, const n_Procs_s* r ) # algebraic number from int
        number*  (*cfInitMPZ)(mpz_t i, const n_Procs_s* r)
        number*  (*cfParameter)(int i, const n_Procs_s* r)
        int     (*cfParDeg)(number* n, const n_Procs_s* r)
        int     (*cfSize)(number* n, const n_Procs_s* r)
        int     (*cfInt)(number* n, const n_Procs_s* r)
        int     (*cdDivComp)(number* a,number* b, const n_Procs_s* r)
        number*  (*cfGetUnit)(number* a, const n_Procs_s* r)
        number*  (*cfExtGcd)(number* a, number* b, number* *s, number* *t , const n_Procs_s* r)

        void (*cfDelete)(number **, const n_Procs_s*)

        number*  (*cfInpNeg)(number* a,  const n_Procs_s* r)
        number*  (*cfInvers)(number* a,  const n_Procs_s* r)
        number*  (*cfCopy)(number* a,  const n_Procs_s* r) # deep copy of algebraic number
        number*  (*cfRePart)(number* a, const n_Procs_s* cf)
        number*  (*cfImPart)(number* a, const n_Procs_s* cf)
        void    (*cfWrite)(number* a, const n_Procs_s* r)
        void    (*cfNormalize)(number* a,  const n_Procs_s* r)

        bint (*cfDivBy)(number* a, number* b, const n_Procs_s* r)
        bint (*cfEqual)(number* a,number* b, const n_Procs_s* )
        bint (*cfIsZero)(number* a, const n_Procs_s* ) # algebraic number comparison with zero
        bint (*cfIsOne)(number* a, const n_Procs_s* )  # algebraic number comparison with one
        bint (*cfIsMOne)(number* a, const n_Procs_s* )
        bint (*cfGreaterZero)(number* a, const n_Procs_s* )
        void (*cfPower)(number* a, int i, number* * result,  const n_Procs_s* r) # algebraic number power

        ring *extRing
        int ch
        mpz_ptr modBase
        unsigned long modExponent

        #n_coeffType type
        int type

    # polynomials

    const char ** n_ParameterNames(const n_Procs_s* r)

    int n_NumberOfParameters(const n_Procs_s* r)

    ctypedef struct poly "spolyrec":
        poly *next
        number *coef
        unsigned long exp[1]

    # ideals

    ctypedef struct ideal "sip_sideal":
        poly **m # gens array
        long rank # rank of module, 1 for ideals
        int nrows # always 1
        int ncols # number of gens

    # polynomial procs
    ctypedef struct p_Procs_s "p_Procs_s":
        pass
    # rings

    ctypedef struct ring "ip_sring":
        int  *order  # array of orderings
        int  *block0 # starting pos
        int  *block1 # ending pos
        int  **wvhdl # weight vectors
        int  OrdSgn  # 1 for polynomial rings
        int  ShortOut # control printing
        int  CanShortOut # control printing capabilities
        number *minpoly # minpoly for base extension field
        char **names # variable names
        p_Procs_s *p_Procs #polxnomial procs
        ideal *qideal #quotient ideal

        short N # number of variables

        int pCompIndex # index of components
        unsigned long bitmask # mask for getting single exponents

        n_Procs_s*    cf # coefficient field/ring
        int ref

        # return total degree of p

        long (*pLDeg)(poly *p, int *l, ring *r)
        long (*pLDegOrig)(poly *p, int *l, ring *r)
        long (*pFDeg)(poly *p, ring *r)
        long (*pFDegOrig)(poly *p, ring *r)

    # a bunch of functions for degree

    long p_FDeg(poly *p, ring *r)
    long p_LDeg(poly *p, int *l, ring *r)
    long p_Deg(poly *p, ring *r)
    long p_WTotaldegree(poly *p, ring *r)
    long p_Totaldegree(poly *p, ring *r)
    long p_WDegree(poly *p, ring *r)

    # available ring orders

    ctypedef struct AlgExtInfo:
        ring * r

    cdef enum rRingOrder_t:
        ringorder_no
        ringorder_a
        ringorder_a64 # for int64 weights
        ringorder_c
        ringorder_C
        ringorder_M
        ringorder_S
        ringorder_s
        ringorder_lp
        ringorder_dp
        ringorder_ip
        ringorder_Dp
        ringorder_wp
        ringorder_Wp
        ringorder_ls
        ringorder_ds
        ringorder_Ds
        ringorder_ws
        ringorder_Ws
        ringorder_L

    # groebner basis options

    cdef enum tHomog:
        isNotHomog
        isHomog
        testHomog

    # ideals

    ctypedef struct ideal "sip_sideal":
        poly **m # gens array
        long rank # rank of module, 1 for ideals
        int nrows # always 1
        int ncols # number of gens

    # dense matrices

    ctypedef struct matrix "ip_smatrix":
        poly **m # gens array
        long rank # rank of module, 1 for normal matrices
        int nrows # number of rows
        int ncols # number of columns

    # integer vectors

    cdef cppclass intvec:
        intvec()
        intvec(int, int, int)
        int *(*ivGetVec)() # this is internal actually
        int (*rows)()
        int (*cols)()
        int (*length)()
        int (*resize)(size_t s)
        int (*get "operator[]")(int i)
        int row
        int col

    cdef cppclass bigintmat:
        int (*length)()
        number* (*get)(int i)

    # omalloc bins

    ctypedef struct omBin "omBin_s"

    # pairs and polys for Groebner Strategy objects

    ctypedef struct LObject:
        pass

    ctypedef struct TObject:
        pass

    # Groebner Strategy objects

    ctypedef struct skStrategy:
        int ak # free module rank
        int sl # last index (inclusive)
        int (*red)(LObject * L,skStrategy *strat)
        void (*initEcart)(LObject *L)
        int (*posInT)(TObject *T,int tl,LObject h) # return the required position in T
        int (*posInL)(LObject *set, int length, LObject* L, skStrategy *strat) # return the required position in L
        void (*enterS)(LObject h, int pos, skStrategy *strat, int atR) # enter into S
        void (*initEcartPair)(LObject * h, poly *f, poly *g, int ecartF, int ecartG)
        int (*posInLOld)(LObject * Ls, int Ll, LObject* Lo, skStrategy *strat)
        LObject P
        ideal *Shdl
        ideal *D
        ideal *M
        poly ** S
        int * ecartS
        int * lenS
        int * fromQ
        unsigned long* sevS
        unsigned long* sevT
        TObject *T
        LObject *L
        LObject *B
        poly*    kHEdge
        poly*    kNoether
        poly*    t_kHEdge
        poly*    kNoetherTail()
        poly*    t_kNoether
        bint *NotUsedAxis
        bint *pairtest
        void *R
        int *S_2_R
        bint noTailReduction

    ctypedef struct data_union:
        ring *uring

    # interpreter objects

    ctypedef struct attr "sattr":
        void (*Init)()
        char *  name
        void *  data
        attr *  next
        int     atyp # the type of the attribute, describes the data field

        void (*Print)()
        attr *(*Copy)() # copy all arguments
        void *(*CopyA)()# copy the data of this attribute
        attr *(*set)(char * s, void * data, int t)
        attr *(*get)(char * s)
        void (*kill)()
        void (*killAll)()

    ctypedef struct idhdl "idrec":
        idhdl *next
        char *id
        int typ
        short lev
        short ref
        data_union data

    ctypedef struct leftv "sleftv":
        leftv *next
        char  *id
        void* data
        #data is some union, so this might be very dangerous, but I am lazy now
        attr *attribute
        void (* Copy)(leftv*)
        void (* Init)()
        void (* CleanUp)(ring *r)
        int  rtyp

    ctypedef struct package "ip_package":
        int language
        idhdl *idroot

    #
    # GLOBAL VARIABLES
    #

    # current ring

    cdef ring *currRing
    cdef ideal *currQuotient
    # omalloc bin for numbers

    cdef omBin *rnumber_bin

    # omalloc bin for rings

    cdef omBin *sip_sring_bin

    # omalloc bin for lists

    cdef omBin *slists_bin

    # integer conversion constant

    cdef long SR_INT

    cdef package *basePack
    cdef package *currPack

    cdef idhdl *basePackHdl
    cdef idhdl *currPackHdl
    cdef idhdl *currRingHdl

    cdef int errorreported
    cdef int si_opt_2    #  previously 'verbose'
    cdef void * currentVoice
    cdef int myynest
    cdef int printlevel

    ctypedef char * const_char_ptr "const char *"
    cdef extern void (*WerrorS_callback)(const_char_ptr)

    #
    # FUNCTIONS
    #

    # singular init

    int siInit(char *)

    ctypedef short (*cfInitCharProc)(coeffs, void *)

    n_coeffType nRegister(n_coeffType n, cfInitCharProc p)

    # external resource init

    void feInitResources(char *name)

    # external resource query

    char* feGetResource(const char id)

    void *omAlloc(size_t size)

    # calloc

    void *omAlloc0(size_t size)

    # typed malloc from bin

    void *omAllocBin(omBin *bin)
    void omFreeBin(void* foo, void* bin)

    # typed calloc from bin

    void *omAlloc0Bin(omBin *bin)

    # strdup

    char *omStrDup(char *)

    # free

    void omFree(void *)

    void omfree(void *)

    # construct ring with characteristic, number of vars and names

    ring *rDefault(int char , int nvars, char **names)
    ring *rDefault(const n_Procs_s* cf, int nvars, char **names)
    ring *rDefault(int ch             , int nvars, char **names,int ord_size, rRingOrder_t *ord, int *block0, int *block1, int **wvhdl)
    ring *rDefault(const n_Procs_s* cf, int nvars, char **names,int ord_size, rRingOrder_t *ord, int *block0, int *block1, int **wvhdl)

    # see coeffs.h
    ctypedef struct  GFInfo:
        int GFChar
        int GFDegree
        const char* GFPar_name

    # parameter is pointer to gGFInfo
    #
    n_Procs_s* nInitChar(n_coeffType t, void * parameter)

    # ring destructor

    void rDelete(ring *r)

    # return characteristic of r

    int rChar(ring *r)

    # return name of the i-th variable of ring r

    char* rRingVar(short i, ring *r)

    # before changing a ring struct, call this

    void rUnComplete(ring *r)

    # after changing a ring struct, call this, computes internal
    # representation for monomials etc.
    void rComplete(ring *r, int force)

    # deep copy of ring

    ring *rCopy0(ring *)

    # change current ring to r

    void rChangeCurrRing(ring *r)

    # return True if ring has components

    int rRing_has_Comp(ring *r)

    int rHasGlobalOrdering(ring *r)
    int rHasLocalOrMixedOrdering(ring *r)

    # return new empty monomial

    poly *p_Init(ring *r)

    # Returns new monomial with exp vector copied from Lm(p)
    poly *p_LmInit(poly *p, ring *r)

    # return constant polynomial from int

    poly *p_ISet(int i, ring *r)

    # return constant polynomial from number

    poly *p_NSet(number *n,ring *r)

    # destructor for polynomials

    void p_Delete(poly **p, ring *r)

    # set the coefficient n for the current list element p in r

    int p_SetCoeff(poly *p, number *n, ring *r)

    # set the coefficient n for the current list element p in r
    # without deleting the current coefficient first

    int p_SetCoeff0(poly *p, number *n, ring *r)

    # get the coefficient of the current list element p in r

    number *p_GetCoeff(poly *p, ring *r)

    # get the module component
    int p_GetComp(poly *p, ring *r)

    # set the module component
    void p_SetComp(poly *p, int v, ring *r)

    # set the exponent e at index v for the monomial p, v starts at 1

    int p_SetExp(poly *p, int v, int e, ring *r)

    # get the exponent at index v of the monomial p in r, v starts at 1

    int p_GetExp(poly *p, int v, ring *r)

    # get the maximal exponent in p

    unsigned long p_GetMaxExp(poly *p, ring *r)

    # if SetExp is called on p, p_Setm needs to be called afterwards to finalize the change.

    void p_Setm(poly *p, ring *r)

    # gets a component out of a polynomial vector

    poly *pTakeOutComp(poly **, int)

    # deep copy p

    poly *p_Copy(poly *p, ring *r)

    # homogenizes p by multiplying certain powers of the varnum-th variable

    poly *p_Homogen (poly *p, int varnum, ring *r)

    # return whether a polynomial is homogeneous

    int p_IsHomogeneous(poly *p, const  ring *r)

    # return string representation of p

    char *p_String(poly *p, ring *r, ring *r)

    # normalize p, needed e.g. for polynomials over the rationals

    void p_Normalize(poly *p, ring *r)

    # makes it so that leading coefficient == 1
    void pNorm(poly *p)

    # return -p, p is destroyed

    poly *p_Neg(poly *p, ring *r)

    # return p*n, p is const (i.e. copied)

    poly *pp_Mult_nn(poly *p, number *n, ring *r)

    # return p*m, does neither destroy p nor m

    poly *pp_Mult_mm(poly *p, poly *m, ring *r)

    # return p+q, destroys p and q

    poly *p_Add_q(poly *p, poly *q, ring *r)

    # return p - m*q, destroys p; const: q,m

    poly *p_Minus_mm_Mult_qq(poly *p, poly *m, poly *q, ring *r)

    # return p + m*q destroys p, const: q, m

    poly *p_Plus_mm_Mult_qq(poly *p, poly *m, poly *q, ring *r)

    # return p*q, does neither destroy p nor q
    poly *pp_Mult_qq(poly *p, poly *q, ring *r)

    # return p*q, destroys p and q
    poly *p_Mult_q(poly *p, poly *q, ring *r)

    # polynomial division, ignoring the remainder
    # via singclap_pdivide resp. idLift, destroys p,q
    poly *p_Divide(poly *p, poly *q, ring *r)

    # divide monomial p by monomial q, p,q const

    poly *pMDivide(poly *p,poly *q)
    poly *p_MDivide(poly *p, poly *q, ring *r)

    # return the i-th power of p; p destroyed, requires global ring

    poly *pPower(poly *p, int exp)

    # return new copy of lm(p), coefficient copied, next=NULL, p may be NULL

    poly *p_Head(poly *p, ring *r)

    # return TRUE, if leading monom of a divides leading monom of b
    # i.e., if there exists a expvector c > 0, s.t. b = a + c

    int p_DivisibleBy(poly *a, poly *b, ring *r)

    # like pDivisibleBy, except that it is assumed that a!=NULL, b!=NULL

    int p_LmDivisibleBy(poly *a, poly *b, ring *r)

    # least common multiplies for monomials only, result is written to m
    # p_Setm must be called on m afterwards.

    void pLcm(poly *a, poly *b, poly *m)

    # total degree of p

    long p_Totaldegree(poly *p, ring *r)

    long pLDeg1_Totaldegree(poly * p,int *l, ring * r)

    # iterate through the monomials of p

    poly *pNext(poly *p)

    # Returns the number of monomials in the poly
    int pLength(poly *a)

    # compare l and r

    int p_Cmp(poly *l, poly *r, ring *r)

    # compare exponent vectors only

    int p_ExpVectorEqual(poly *p, poly *m, ring *r)

    # TRUE if poly is constant

    int p_IsConstant(poly *, ring *)

    # like p_IsConstant but p must be !=NULL

    int p_LmIsConstant(poly *p, ring *)

    # TRUE if poly is unit

    int p_IsUnit(poly *p, ring *r)

    # TRUE if poly  is one

    bint p_IsOne(const poly *p, const ring *r)

    # substitute monomial for variable given by varidx in poly

    poly *pSubst(poly *p, int varidx, poly *value)

    # inverse of poly, if possible

    poly *p_Series(int n, poly *, poly *, intvec *, ring *r)

    # gcd of f and g

    poly *singclap_gcd ( poly *f, poly *g, ring * r )

    # resultant of f and g in x

    poly *singclap_resultant ( poly *f, poly *g , poly *x, ring * r)

    # extended gcd of f and g

    int singclap_extgcd( poly *f, poly *g, poly *res, poly *pa, poly *pb, ring * r )

    # full polynomial division (as opposed to monomial division)

    poly *singclap_pdivide ( poly *f, poly *g, ring * r )

    # factorization

    ideal *singclap_factorize ( poly *f, intvec ** v , int with_exps, ring * r)

    # return determinant of i
    poly *singclap_det(matrix *i, ring * r)

    # normal form calculation of p with respect to i, q is quotient
    # ring.

    poly *kNF(ideal *i, ideal *q, poly *p)
    poly *kNF(ideal *i, ideal *q, poly *p, int syzComp, int lazyReduce)

    # derive p with respect to i-th variable

    poly *pDiff(poly *p, int i)

    # TRUE if p is a vector

    int pIsVector(poly *p)

    # return current component level

    int pGetComp(poly *p)

    # general number constructor

    number *n_Init(int n, n_Procs_s *cf)

    # general number destructor

    void n_Delete(number **n, n_Procs_s *cf)

    # Copy this number
    number *n_Copy(number *n, n_Procs_s *cf)

    # Invert this number
    int n_IsUnit(number *n, const n_Procs_s *cf)
    number *n_Invers(number *n, const n_Procs_s *cf)

    # Characteristic of coefficient domain
    int n_GetChar(const n_Procs_s *cf)

    # rational number from int

    number *nlInit(int)

    # rational number from int

    number *nlRInit(int)

    # rational number from numerator and denominator

    number *nlInit2(int i, int j,const n_Procs_s* cf)

    # simplify rational number (cancel common factors)

    number *nlNormalize(number *)

    # copy a number

    number *nlCopy(number *)

    # number to integer handle

    long SR_TO_INT(number *)

    # mpz_t to integer handle

    long SR_HDL(number *)

    # get C int from ZmodN
    int nrnInt(number *v)

    # ideal constructor

    ideal *idInit(int size, int rank)

    # ideal destructor

    void id_Delete(ideal **, ring *)

    # lifting

    ideal *idLift(ideal *mod, ideal *submod, ideal **rest, int goodShape, int isSB, int divide)

    # number of generators of ideal

    int IDELEMS(ideal *i)

    # remove zero entries

    void idSkipZeroes (ideal *ide)

    # rank of free module for m

    long id_RankFreeModule(ideal *m, ring *r)

    # buchberger's algorithm

    ideal *kStd(ideal *i, ideal *q, tHomog h, intvec *w)

    # slimgb algorithm

    ideal *t_rep_gb(ring *r,ideal *arg_I, int syz_comp, int F4_mode)

    # interreduction

    ideal *kInterRed(ideal *i, ideal *q)

    # TRUE if ideal is module

    int idIsModule(ideal *m, ring *r)

    # convert module to matrix

    matrix *idModule2Matrix(ideal *i)

    # convert matrix to module

    ideal * idMatrix2Module(matrix *m)

    # dense matrix constructor

    matrix *mpNew(int i, int j)

    # gauss-bareiss algorithm

    void smCallBareiss(ideal *, int, int, ideal *, intvec**)

    # determinant
    poly *smCallDet(ideal *)

    # check which det algorithm to choose
    int smCheckDet(ideal *, int, int)

    # get element at row i and column j

    poly *MATELEM(matrix *, int i, int j)

    # number columns of matrix

    int MATCOLS(matrix *)

    # number for rows of matrix

    int MATROWS(matrix *)

    # Groebner Strategy functions

    void initBuchMoraCrit(skStrategy *strat)
    void initEcartNormal (LObject* h)
    void initEcartBBA (LObject* h)
    void initEcartPairBba (LObject* Lp,poly *f,poly *g,int ecartF,int ecartG)

    # init strat with F and Q
    void initS (ideal *F, ideal *Q,skStrategy *strat)

    void enterL (LObject *set,int *length, int *LSetmax, LObject p,int at)
    void enterSBba (LObject p,int atS,skStrategy *strat, int atR)

    skStrategy * new_skStrategy "new skStrategy"()
    void * delete_skStrategy "delete "(skStrategy *doomed) # needs currRing

    # test validity of strat
    void kTest(skStrategy *strat)

    # head reduction
    poly *redNF(poly *p, int index, int nonorm, skStrategy *strat)
    # tail reduction
    poly *redtailBba(poly *p, int index, skStrategy *strat)

    cdef int CMD_1
    cdef int CMD_2
    cdef int CMD_12
    cdef int CMD_3
    cdef int CMD_13
    cdef int CMD_23
    cdef int CMD_123
    cdef int CMD_M
    cdef int ROOT_DECL
    cdef int ROOT_DECL_LIST
    cdef int RING_DECL
    cdef int RING_DECL_LIST
    cdef int INT_CMD
    cdef int INTMAT_CMD
    cdef int POLY_CMD
    cdef int PROC_CMD
    cdef int RING_CMD
    cdef int QRING_CMD
    cdef int STRING_CMD
    cdef int VECTOR_CMD
    cdef int IDEAL_CMD
    cdef int MODUL_CMD
    cdef int NUMBER_CMD
    cdef int MATRIX_CMD
    cdef int LIST_CMD
    cdef int INTVEC_CMD
    cdef int BIGINTVEC_CMD
    cdef int NONE
    cdef int RESOLUTION_CMD
    cdef int PACKAGE_CMD

    cdef int V_SHOW_MEM
    cdef int V_YACC
    cdef int V_REDEFINE
    cdef int V_READING
    cdef int V_LOAD_LIB
    cdef int V_DEBUG_LIB
    cdef int V_LOAD_PROC
    cdef int V_DEF_RES
    cdef int V_SHOW_USE
    cdef int V_IMAP
    cdef int V_PROMPT
    cdef int V_NSB
    cdef int V_CONTENTSB
    cdef int V_CANCELUNIT
    cdef int V_DEG_STOP

    #
    # INTERPRETER
    #
    leftv iiRETURNEXPR

    cdef omBin* sleftv_bin

    int     IsCmd(char *n, int  tok)

    idhdl* ggetid(char *n)

    bint iiMake_proc(idhdl *pn, package *pack, leftv *sl)

    bint iiExprArith1(leftv *res, leftv* a, int op)
    bint iiExprArith2(leftv *res, leftv* a, int op, leftv *b, bint proc_call)
    bint iiExprArith3(leftv *res, int op, leftv *a, leftv *b, leftv *c)
    bint iiExprArithM(leftv *res, leftv* a, int op)

    bint iiLibCmd(char *, bint autoexport, bint tellerror, bint force)

    package *IDPACKAGE(idhdl *)

    cdef idhdl *IDROOT

    cdef idhdl *IDNEXT(idhdl *)

    cdef int IDTYP(idhdl *)

    idhdl *enterid(char * a, int lev, int t, idhdl** root, bint init)

    void atSet(leftv *root, char * name, void * data, int typ)
    void *atGet(leftv *root, char *name, int typ)

    int hasFlag(leftv *A, int F)
    void setFlag(leftv *A, int F)
    void resetFlag(leftv *A, int F)

    ctypedef number* (*nMapFunc)(number *c, const n_Procs_s* src, const n_Procs_s* dst)

cdef extern from "singular/coeffs/coeffs.h":

    number *ndCopyMap(number *, const n_Procs_s* src,const n_Procs_s* dst)

cdef extern from "singular/coeffs/rmodulo2m.h":

    #init 2^m from a long
    number *nr2mMapZp(number *,const n_Procs_s* src,const n_Procs_s* dst)

cdef extern from "singular/kernel/maps/gen_maps.h":

    # mapping from p in r1 by i2 to r2

    poly *maMapPoly(poly *p, ring *r1, ideal *i2, ring *r2, const nMapFunc nMap)

cdef extern from "singular/kernel/maps/fast_maps.h":

    # mapping from ideal i1 in r1 by i2 to r2

    ideal *fast_map_common_subexp(ideal *i1, ring *r1, ideal *i2, ring *r2)

cdef extern from "singular/polys/ext_fields/algext.h":

    naInitChar(n_Procs_s* cf, void * infoStruct)

    nMapFunc naSetMap(const n_Procs_s* src, const n_Procs_s* dst)

cdef extern from "singular/coeffs/rmodulon.h":

    # init ZmodN from GMP
    number *nrnMapGMP(number *v,const n_Procs_s* src,const n_Procs_s* dst)

    nMapFunc nrnSetMap(const n_Procs_s* src,const n_Procs_s* dst)

cdef extern from "singular/coeffs/rmodulon.h":
    # see rmodulon.h

    ctypedef struct ZnmInfo:
       mpz_ptr base
       unsigned long exp

cdef extern from "singular/coeffs/rintegers.h":

    # init integer
    number *nrzInit(int i, const n_Procs_s* cf)

cdef extern from "singular/polys/weight.h":

    double wFunctionalBuch(int *degw, int *lpol, int npol, double *rel, double wx, double wNsqr)

cdef extern from "singular/polys/prCopy.h":
    poly *prCopyR_NoSort(poly *p, ring *r, ring *dest_r)
    poly *prCopyR(poly *p, ring *r, ring *dest_r)

    cdef int LANG_TOP

cdef extern from "singular/polys/sbuckets.h":
    #sBucket is actually a class, but we handle it opaquely, so we call it a "struct" here.
    ctypedef struct sBucket:
        pass

    #create an sBucket
    sBucket *sBucketCreate(ring *r)

    #destroy an sBucket (note: pointer to pointer)
    void sBucketDestroy(sBucket **bucket)

    #merge contents of sBucket into polynomial and clear bucket
    #(use when monomials are distinct).
    #assumes length <= 0 || pLength(p) == length
    void sBucketClearMerge(sBucket *bucket, poly **p, int *length)

    #add contents of sBucket into polynomial and clear bucket
    #(can handle repeated monomials)
    void sBucketClearAdd(sBucket *bucket, poly **p, int *length)

    #inline versions that in addition clear the pointer bucket afterwards
    void sBucketDestroyMerge(sBucket *bucket, poly **p, int *length)
    void sBucketDestroyAdd(sBucket *bucket, poly *p, int *length)

    #delete bucket constant and clear pointer
    void sBucketDeleteAndDestroy(sBucket **bucket_pt)

    #merge p into bucket (distinct monomials assumed)
    #destroys poly in the process
    void sBucket_Merge_p(sBucket *bucket, poly *p, int lp)

    #merge p into bucket  (distinct monomials assumed)
    #destroys poly in the process
    void sBucket_Merge_m(sBucket *bucket, poly *p)

    #adds p into bucket (distinct monomials assumed)
    #destroys poly in the process
    void sBucket_Add_p(sBucket *bucket, poly *p, int lp)

    #adds p into bucket (distinct monomials assumed)
    #destroys poly in the process
    void sBucket_Add_m(sBucket *bucket, poly *p)

    #sorts p with bucketSort: assumes all monomials of p are different
    poly *sBucketSortMerge(poly *p, const ring *r)

    #sorts p with bucketSort: p may have equal monomials
    poly *sBucketSortAdd(poly *p, const ring *r)

cdef extern from "singular/polys/nc/nc.h":
    # Non-commutative functions
    ctypedef enum nc_type:
        nc_error
        nc_general
        nc_skew
        nc_comm
        nc_lie
        nc_undef
        nc_exterior

    void ncRingType(ring *, nc_type)
    nc_type ncRingType_get "ncRingType" (ring *)
    int nc_CallPlural(matrix* CC, matrix* DD, poly* CN, poly* DN, ring* r)
    bint nc_SetupQuotient(ring *, ring *, bint)

cdef extern from "singular/coeffs/longrat.h":

    # get numerator

    number *nlGetNumerator(number *n, const n_Procs_s* cf)

    # get denominator

    number *nlGetDenom(number *n, const n_Procs_s* cf)

    # rational number from numerator and denominator

    number *nlInit2gmp(mpz_t n, mpz_t d,const n_Procs_s* cf)

    # delete rational number

    void nlDelete(number **n, const n_Procs_s* cf)

cdef extern from "singular/polys/nc/sca.h":
    void sca_p_ProcsSet(ring *, p_Procs_s *)
    void scaFirstAltVar(ring *, int)
    void scaLastAltVar(ring *, int)

cdef extern from "singular/polys/monomials/ring.h":
    bint rIsPluralRing(ring* r)
    void rPrint "rWrite"(ring* r)
    char* rOrderingString "rOrdStr"(ring* r)
    void pDebugPrint "p_DebugPrint" (poly*p, ring* r)

cdef extern from "singular/kernel/combinatorics/stairc.h":
    # Computes the monomial basis for R[x]/I
    ideal *scKBase(int deg, ideal *s, ideal *Q)

cdef extern from "singular/Singular/lists.h":
    ctypedef struct lists "slists":
        int    nr
        leftv  *m
        void (*Init)(int n)

cdef extern from "singular/kernel/GBEngine/kstd1.h":
    cdef extern int Kstd1_deg   # degBound, default 0
    cdef extern int Kstd1_mu    # multBound, default 0

cdef extern from "singular/kernel/ideals.h":
    ctypedef ideal * resolvente

cdef extern from "singular/kernel/GBEngine/syz.h":
    ctypedef struct syStrategy "ssyStrategy":
        resolvente fullres;
        resolvente minres;
        int length;
        int regularity;
        short list_length;
        short references

cdef extern from "singular/polys/ext_fields/transext.h":
    ctypedef struct TransExtInfo:
        ring * r
