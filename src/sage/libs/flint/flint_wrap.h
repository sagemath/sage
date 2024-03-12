/* sage_setup: distribution = sagemath-flint
 */
/* WARNING: src/sage/libs/flint/flint_wrap.h is generated from
 * src/sage_setup/autogen/flint/templates/flint_wrap.h.template
 * please make sure that you are modifying the correct file! */

#ifndef SAGE_FLINT_WRAP_H
#define SAGE_FLINT_WRAP_H
/* Using flint headers together in the same module as headers from
 * some other libraries (pari, possibly others) as it defines the
 * macros ulong and slong all over the place.
 *
 * What's worse is they are defined to types from GMP (mp_limb_t and
 * mp_limb_signed_t respectively) which themselves can have system-dependent
 * typedefs, so there is no guarantee that all these 'ulong' definitions from
 * different libraries' headers will be compatible.
 *
 * When including flint headers in Sage it should be done through this wrapper
 * to prevent confusion.  We rename flint's ulong and slong to fulong and
 * fslong.  This is consistent with flint's other f-prefixed typedefs.
 */

#include <gmp.h>
#include <mpfr.h>

/* Save previous definition of ulong if any, as pari also uses it */
/* Should work on GCC, clang, MSVC */
#pragma push_macro("ulong")
#undef ulong

#include <flint/flint.h>

/* If flint was already previously included via another header (e.g.
 * arb_wrap.h) then it may be necessary to redefine ulong and slong again */

#ifndef ulong
#define ulong mp_limb_t
#define slong mp_limb_signed_t
#endif

#include <flint/acb.h>
#include <flint/acb_calc.h>
#include <flint/acb_dft.h>
#include <flint/acb_dirichlet.h>
#include <flint/acb_elliptic.h>
#include <flint/acb_hypgeom.h>
#include <flint/acb_mat.h>
#include <flint/acb_modular.h>
#include <flint/acb_poly.h>
#include <flint/acf.h>
#include <flint/aprcl.h>
#include <flint/arb.h>
#include <flint/arb_calc.h>
#include <flint/arb_fmpz_poly.h>
#include <flint/arb_fpwrap.h>
#include <flint/arb_hypgeom.h>
#include <flint/arb_mat.h>
#include <flint/arb_poly.h>
#include <flint/arf.h>
#include <flint/arith.h>
#include <flint/bernoulli.h>
#include <flint/bool_mat.h>
#include <flint/ca.h>
#include <flint/ca_ext.h>
#include <flint/ca_field.h>
#include <flint/ca_mat.h>
#include <flint/ca_poly.h>
#include <flint/ca_vec.h>
#include <flint/calcium.h>
#include <flint/d_mat.h>
#include <flint/d_vec.h>
#include <flint/dirichlet.h>
#include <flint/dlog.h>
#include <flint/double_extras.h>
#include <flint/double_interval.h>
#include <flint/fexpr.h>
#include <flint/fexpr_builtin.h>
#include <flint/fft.h>
#include <flint/flint.h>
#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>
#include <flint/fmpq_mpoly.h>
#include <flint/fmpq_mpoly_factor.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpq_vec.h>
#include <flint/fmpz.h>
#include <flint/fmpz_extras.h>
#include <flint/fmpz_factor.h>
#include <flint/fmpz_lll.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_mod.h>
#include <flint/fmpz_mod_mat.h>
#include <flint/fmpz_mod_mpoly.h>
#include <flint/fmpz_mod_mpoly_factor.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/fmpz_mod_poly_factor.h>
#include <flint/fmpz_mod_vec.h>
#include <flint/fmpz_mpoly.h>
#include <flint/fmpz_mpoly_factor.h>
#include <flint/fmpz_mpoly_q.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_factor.h>
#include <flint/fmpz_poly_mat.h>
#include <flint/fmpz_poly_q.h>
#include <flint/fmpz_vec.h>
#include <flint/fmpzi.h>
#include <flint/fq.h>
#include <flint/fq_default.h>
#include <flint/fq_default_mat.h>
#include <flint/fq_default_poly.h>
#include <flint/fq_default_poly_factor.h>
#include <flint/fq_embed.h>
#include <flint/fq_mat.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_embed.h>
#include <flint/fq_nmod_mat.h>
#include <flint/fq_nmod_mpoly.h>
#include <flint/fq_nmod_mpoly_factor.h>
#include <flint/fq_nmod_poly.h>
#include <flint/fq_nmod_poly_factor.h>
#include <flint/fq_nmod_vec.h>
#include <flint/fq_poly.h>
#include <flint/fq_poly_factor.h>
#include <flint/fq_vec.h>
#include <flint/fq_zech.h>
#include <flint/fq_zech_embed.h>
#include <flint/fq_zech_mat.h>
#include <flint/fq_zech_poly.h>
#include <flint/fq_zech_poly_factor.h>
#include <flint/fq_zech_vec.h>
#include <flint/gr.h>
#include <flint/gr_generic.h>
#include <flint/gr_mat.h>
#include <flint/gr_mpoly.h>
#include <flint/gr_poly.h>
#include <flint/gr_special.h>
#include <flint/gr_vec.h>
#include <flint/hypgeom.h>
#include <flint/long_extras.h>
#include <flint/mag.h>
#include <flint/mpfr_mat.h>
#include <flint/mpfr_vec.h>
#include <flint/mpn_extras.h>
#include <flint/mpoly.h>
#include <flint/nf.h>
#include <flint/nf_elem.h>
#include <flint/nmod.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_mpoly.h>
#include <flint/nmod_mpoly_factor.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_factor.h>
#include <flint/nmod_poly_mat.h>
#include <flint/nmod_types.h>
#include <flint/nmod_vec.h>
#include <flint/padic.h>
#include <flint/padic_mat.h>
#include <flint/padic_poly.h>
#include <flint/partitions.h>
#include <flint/perm.h>
#include <flint/profiler.h>
#include <flint/qadic.h>
#include <flint/qfb.h>
#include <flint/qqbar.h>
#include <flint/qsieve.h>
#include <flint/thread_pool.h>
#include <flint/ulong_extras.h>

#undef ulong
#undef slong
#undef mp_bitcnt_t

#pragma pop_macro("ulong")

/* CPU_SIZE_1 and SIZE_RED_FAILURE_THRESH are defined as macros in flint/fmpz_lll.h
 * and as variables in fplll/defs.h, which breaks build if linbox is compiled with fplll */

#undef CPU_SIZE_1
#undef SIZE_RED_FAILURE_THRESH

#endif
