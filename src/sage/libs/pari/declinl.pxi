r"""
Declarations for inline functions from PARI.

This file contains all declarations from headers/pariinl.h from the
PARI distribution.  All these functions are simple inline functions.
This file is included by sage/libs/pari/decl.pxi


AUTHORS:

 - Jeroen Demeyer (2010-08-15): initial version (#9898)

"""

cdef extern from "sage/libs/pari/parisage.h":

    ###################################################################
    #                                                                 #
    #                          CONSTRUCTORS                           #
    #                                                                 #
    ###################################################################

    GEN     mkintmod(GEN x, GEN y)
    GEN     mkintmodu(ulong x, ulong y)
    GEN     mkpolmod(GEN x, GEN y)
    GEN     mkfrac(GEN x, GEN y)
    GEN     mkfraccopy(GEN x, GEN y)
    GEN     mkrfrac(GEN x, GEN y)
    GEN     mkcomplex(GEN x, GEN y)
    GEN     gen_I()
    GEN     cgetc(long l)
    GEN     mkquad(GEN n, GEN x, GEN y)
    GEN     mkvecsmall(long x)
    GEN     mkvecsmall2(long x,long y)
    GEN     mkvecsmall3(long x,long y,long z)
    GEN     mkvecsmall4(long x,long y,long z,long t)
    GEN     mkvec(GEN x)
    GEN     mkvec2(GEN x, GEN y)
    GEN     mkvec3(GEN x, GEN y, GEN z)
    GEN     mkvec4(GEN x, GEN y, GEN z, GEN t)
    GEN     mkvec5(GEN x, GEN y, GEN z, GEN t, GEN u)
    GEN     mkvecs(long x)
    GEN     mkvec2s(long x, long y)
    GEN     mkvec3s(long x, long y, long z)
    GEN     mkveccopy(GEN x)
    GEN     mkvec2copy(GEN x, GEN y)
    GEN     mkcol(GEN x)
    GEN     mkcol2(GEN x, GEN y)
    GEN     mkcolcopy(GEN x)
    GEN     mkmat(GEN x)
    GEN     mkmat2(GEN x, GEN y)
    GEN     mkmatcopy(GEN x)
    GEN     pol_x(long v)
    GEN     pol_1(long v)
    GEN     const_vec(long n, GEN x)
    GEN     const_col(long n, GEN x)
    GEN     const_vecsmall(long n, long c)

    ### Zero ###
    GEN     zeropadic(GEN p, long e)
    GEN     zeroser(long v, long e)
    GEN     zeropol(long v)
    GEN     zerocol(long n)
    GEN     zerovec(long n)
    GEN     zeromat(long m, long n)
    GEN     zero_Flx(long sv)
    GEN     zero_Flv(long n)
    GEN     zero_Flm(long m, long n)
    GEN     zero_F2v(long m)
    GEN     zero_F2m(long m, long n)
    GEN     zero_F2m_copy(long m, long n)
    GEN     zeromatcopy(long m, long n)

    GEN     col_ei(long n, long i)
    GEN     vec_ei(long n, long i)
    GEN     vecsmall_ei(long n, long i)
    GEN     Rg_col_ei(GEN x, long n, long i)
    GEN     shallowcopy(GEN x)
    GEN     vecsmall_copy(GEN x)
    GEN     vectrunc_init(long l)
    void    vectrunc_append(GEN x, GEN t)
    GEN     vecsmalltrunc_init(long l)
    void    vecsmalltrunc_append(GEN x, long t)

    ###################################################################
    #                                                                 #
    #                        VEC / COL / VECSMALL                     #
    #                                                                 #
    ###################################################################

    GEN     vec_shorten(GEN v, long n)
    GEN     vec_lengthen(GEN v, long n)
    GEN     vec_setconst(GEN v, GEN x)
    GEN     vecsmall_shorten(GEN v, long n)
    GEN     vecsmall_lengthen(GEN v, long n)
    GEN     vec_to_vecsmall(GEN z)
    GEN     vecsmall_to_vec(GEN z)
    GEN     vecsmall_to_col(GEN z)
    int     vecsmall_lexcmp(GEN x, GEN y)
    int     vecsmall_prefixcmp(GEN x, GEN y)
    GEN     vecsmall_prepend(GEN V, long s)
    GEN     vecsmall_append(GEN V, long s)
    GEN     vecsmall_concat(GEN u, GEN v)
    long    vecsmall_coincidence(GEN u, GEN v)
    long    vecsmall_isin(GEN v, long x)
    long    vecsmall_pack(GEN V, long base, long mod)
    long    vecsmall_max(GEN x)
    long    vecsmall_min(GEN x)
    int     ZV_isscalar(GEN x)
    int     QV_isscalar(GEN x)
    int     RgV_isscalar(GEN x)
    int     RgX_isscalar(GEN x)
    int     RgX_is_rational(GEN x)
    int     RgX_is_ZX(GEN x)
    int     RgX_is_monomial(GEN x)
    int     RgM_is_ZM(GEN x)

    ###################################################################
    #                                                                 #
    #            Dynamic arrays implementation                        #
    #                                                                 #
    ###################################################################

    # Omitted

    ###################################################################
    #                                                                 #
    #                            EXTRACT                              #
    #                                                                 #
    ###################################################################

    GEN     vecslice(GEN A, long y1, long y2)
    GEN     vecslicepermute(GEN A, GEN p, long y1, long y2)
    GEN     rowslicepermute(GEN A, GEN p, long x1, long x2)
    GEN     rowslice(GEN A, long x1, long x2)
    GEN     row(GEN A, long x0)
    GEN     row_Flm(GEN A, long x0)
    GEN     rowcopy(GEN A, long x0)
    GEN     row_i(GEN A, long x0, long x1, long x2)
    GEN     vecreverse(GEN A)
    GEN     vecpermute(GEN A, GEN p)
    GEN     rowpermute(GEN A, GEN p)
    void    vecselect_p(GEN A, GEN B, GEN p, long init, long lB)
    void    rowselect_p(GEN A, GEN B, GEN p, long init)

    ###################################################################
    #                                                                 #
    #                          PERMUTATIONS                           #
    #                                                                 #
    ###################################################################

    GEN     identity_perm(long n)
    GEN     cyclic_perm(long n, long d)
    GEN     perm_mul(GEN s, GEN t)
    GEN     perm_inv(GEN x)
    GEN     perm_conj(GEN s, GEN t)

    ###################################################################
    #                                                                 #
    #                      MALLOC/FREE WRAPPERS                       #
    #                                                                 #
    ###################################################################

    void    pari_free(void *pointer)
    void*   pari_malloc(size_t size)
    void*   pari_realloc(void *pointer, size_t size)
    void*   pari_calloc(size_t size)
    GEN     cgetalloc(long t, size_t l)

    ###################################################################
    #                                                                 #
    #                       GARBAGE COLLECTION                        #
    #                                                                 #
    ###################################################################

    GEN     icopy_avma(GEN x, pari_sp av)
    GEN     gerepileuptoleaf(pari_sp av, GEN x)
    GEN     gerepileuptoint(pari_sp av, GEN x)
    GEN     gerepileupto(pari_sp av, GEN x)
    GEN     gerepilecopy(pari_sp av, GEN x)
    void    gerepilemany(pari_sp av, GEN* gptr[], int n)
    void    gerepileall(pari_sp av, int n, ...)
    void    gerepilecoeffs(pari_sp av, GEN x, int n)
    void    gerepilecoeffs2(pari_sp av, GEN x, int n, GEN y, int o)
    void    cgiv(GEN x)
    void    killblock(GEN x)
    int     is_universal_constant(GEN x)

    ###################################################################
    #                                                                 #
    #                    CONVERSION / ASSIGNMENT                      #
    #                                                                 #
    ###################################################################

    GEN     cxcompotor(GEN z, long prec)
    GEN     cxtofp(GEN x, long prec)
    double  gtodouble(GEN x)
    long    gtos(GEN x)
    GEN     absfrac(GEN x)
    GEN     Q_abs(GEN x)
    GEN     gtofp(GEN z, long prec)
    GEN     RgX_gtofp(GEN x, long prec)
    GEN     RgC_gtofp(GEN x, long prec)
    GEN     RgM_gtofp(GEN x, long prec)
    GEN     RgX_fpnorml2(GEN x, long prec)
    GEN     RgC_fpnorml2(GEN x, long prec)
    GEN     RgM_fpnorml2(GEN x, long prec)
    void    affgr(GEN x, GEN y)
    GEN     affc_fixlg(GEN x, GEN res)
    GEN     trunc_safe(GEN x)

    ###################################################################
    #                                                                 #
    #                          LENGTH CONVERSIONS                     #
    #                                                                 #
    ###################################################################

    # Omitted

    ###################################################################
    #                                                                 #
    #                      OPERATIONS MODULO m                        #
    #                                                                 #
    ###################################################################

    GEN     Fp_red(GEN a, GEN m)
    GEN     Fp_add(GEN a, GEN b, GEN m)
    GEN     Fp_sub(GEN a, GEN b, GEN m)
    GEN     Fp_neg(GEN b, GEN m)
    GEN     Fp_center(GEN u, GEN p, GEN ps2)
    GEN     Fp_mul(GEN a, GEN b, GEN m)
    GEN     Fp_sqr(GEN a, GEN m)
    GEN     Fp_mulu(GEN a, ulong b, GEN m)
    GEN     Fp_inv(GEN a, GEN m)
    GEN     Fp_invsafe(GEN a, GEN m)
    GEN     Fp_div(GEN a, GEN b, GEN m)

    ###################################################################
    #                                                                 #
    #                          GEN SUBTYPES                           #
    #                                                                 #
    ###################################################################

    int     is_const_t(long t)
    int     is_extscalar_t(long t)
    int     is_intreal_t(long t)
    int     is_matvec_t(long t)
    int     is_noncalc_t(long tx)
    int     is_rational_t(long t)
    int     is_recursive_t(long t)
    int     is_scalar_t(long t)
    int     is_vec_t(long t)

    ###################################################################
    #                                                                 #
    #                         TRANSCENDENTAL                          #
    #                                                                 #
    ###################################################################

    GEN     sqrtr(GEN x)
    GEN     sqrtnr(GEN x, long n)

    ###################################################################
    #                                                                 #
    #                         MISCELLANEOUS                           #
    #                                                                 #
    ###################################################################

    int     isintzero(GEN x)
    int     isint1(GEN x)
    int     isintm1(GEN x)
    int     equali1(GEN n)
    int     equalim1(GEN n)
    int     is_pm1(GEN n)
    int     is_bigint(GEN n)

    # Many functions omitted

    ### POLYNOMIALS
    GEN     constant_term(GEN x)
    GEN     leading_term(GEN x)
    long    degpol(GEN x)
    long    lgpol(GEN x)
    GEN     truecoeff(GEN x, long n)

    ###################################################################
    #                                                                 #
    #                             ASSIGNMENTS                         #
    #                                                                 #
    ###################################################################

    # Omitted

    ###################################################################
    #                                                                 #
    #                       ELLIPTIC CURVES                           #
    #                                                                 #
    ###################################################################

    GEN     ell_get_a1(GEN e)
    GEN     ell_get_a2(GEN e)
    GEN     ell_get_a3(GEN e)
    GEN     ell_get_a4(GEN e)
    GEN     ell_get_a6(GEN e)
    GEN     ell_get_b2(GEN e)
    GEN     ell_get_b4(GEN e)
    GEN     ell_get_b6(GEN e)
    GEN     ell_get_b8(GEN e)
    GEN     ell_get_c4(GEN e)
    GEN     ell_get_c6(GEN e)
    GEN     ell_get_disc(GEN e)
    GEN     ell_get_j(GEN e)
    GEN     ell_get_roots(GEN e)

    int     ell_is_inf(GEN z)
    int     ell_is_padic(GEN x)
    int     ell_is_real(GEN x)

    ###################################################################
    #                                                                 #
    #                    ALGEBRAIC NUMBER THEORY                      #
    #                                                                 #
    ###################################################################

    GEN     pr_get_p(GEN pr)
    GEN     pr_get_gen(GEN pr)
    long    pr_get_e(GEN pr)
    long    pr_get_f(GEN pr)
    GEN     pr_get_tau(GEN pr)
    int     pr_is_inert(GEN P)
    GEN     pr_norm(GEN pr)

    long    nf_get_varn(GEN nf)
    GEN     nf_get_pol(GEN nf)
    long    nf_get_degree(GEN nf)
    long    nf_get_r1(GEN nf)
    long    nf_get_r2(GEN nf)
    GEN     nf_get_disc(GEN nf)
    GEN     nf_get_index(GEN nf)
    GEN     nf_get_M(GEN nf)
    GEN     nf_get_G(GEN nf)
    GEN     nf_get_roundG(GEN nf)
    GEN     nf_get_Tr(GEN nf)
    GEN     nf_get_TrInv(GEN nf)
    GEN     nf_get_roots(GEN nf)
    GEN     nf_get_zk(GEN nf)
    GEN     nf_get_invzk(GEN nf)
    void    nf_get_sign(GEN nf, long *r1, long *r2)

    GEN     bnf_get_nf(GEN bnf)
    GEN     bnf_get_clgp(GEN bnf)
    GEN     bnf_get_no(GEN bnf)
    GEN     bnf_get_cyc(GEN bnf)
    GEN     bnf_get_gen(GEN bnf)
    GEN     bnf_get_reg(GEN bnf)
    GEN     bnf_get_logfu(GEN bnf)
    GEN     bnf_get_tuU(GEN bnf)
    long    bnf_get_tuN(GEN bnf)
    GEN     bnf_get_fu(GEN bnf)
    GEN     bnf_get_fu_nocheck(GEN bnf)

    GEN     bnr_get_bnf(GEN bnr)
    GEN     bnr_get_bid(GEN bnr)
    GEN     bnr_get_mod(GEN bnr)
    GEN     bnr_get_nf(GEN bnr)
    GEN     bnr_get_no(GEN bnr)
    GEN     bnr_get_cyc(GEN bnr)
    GEN     bnr_get_gen_nocheck(GEN bnr)
    GEN     bnr_get_gen(GEN bnr)

    GEN     bid_get_mod(GEN bid)
    GEN     bid_get_ideal(GEN bid)
    GEN     bid_get_arch(GEN bid)
    GEN     bid_get_cyc(GEN bid)
    GEN     bid_get_gen_nocheck(GEN bid)
    GEN     bid_get_gen(GEN bid)

    GEN     gal_get_pol(GEN gal)
    GEN     gal_get_p(GEN gal)
    GEN     gal_get_e(GEN gal)
    GEN     gal_get_mod(GEN gal)
    GEN     gal_get_roots(GEN gal)
    GEN     gal_get_invvdm(GEN gal)
    GEN     gal_get_den(GEN gal)
    GEN     gal_get_group(GEN gal)
    GEN     gal_get_gen(GEN gal)
    GEN     gal_get_orders(GEN gal)

    long    rnf_get_degree(GEN rnf)

    GEN     idealpseudomin(GEN I, GEN G)
    GEN     idealpseudomin_nonscalar(GEN I, GEN G)
    GEN     idealred_elt(GEN nf, GEN I)
    GEN     idealred(GEN nf, GEN I)
