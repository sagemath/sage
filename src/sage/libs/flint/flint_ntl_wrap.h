/* sage_setup: distribution = sagemath-flint
 */
#ifndef SAGE_FLINT_NTL_WRAP_H
#define SAGE_FLINT_NTL_WRAP_H
/*
 * Similar to flint_wrap.h but specifically for wrapping the flint-NTL
 * interface.  It is separate from flint_wrap.h so that not every module
 * which uses flint has to pull in NTL headers as well.
 */

#include <gmp.h>

/* Save previous definition of ulong if any, as pari also uses it */
#pragma push_macro("ulong")
#undef ulong

#include <flint/flint.h>

/* If flint was already previously included via another header (e.g.
 * arb_wrap.h) then it may be necessary to redefine ulong and slong again */

#ifndef ulong
#define ulong mp_limb_t
#define slong mp_limb_signed_t
#endif

#include <flint/NTL-interface.h>

#undef ulong
#undef slong
#undef mp_bitcnt_t

#pragma pop_macro("ulong")

#endif
