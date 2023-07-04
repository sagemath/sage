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

#include <flint/arith.h>
#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpz_factor.h>
#include <flint/fmpz_mod.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_poly_mat.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_q.h>
#include <flint/fmpz_vec.h>
#include <flint/fq.h>
#include <flint/fq_nmod.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_vec.h>
#include <flint/padic.h>
#include <flint/padic_poly.h>
#include <flint/qadic.h>
#include <flint/qsieve.h>
#include <flint/ulong_extras.h>

#undef ulong
#undef slong
#undef mp_bitcnt_t

#pragma pop_macro("ulong")

#endif
