#ifndef SAGE_ARB_WRAP_H
#define SAGE_ARB_WRAP_H
/*
 * Similar to flint_wrap.h but specifically for wrapping the headers supplied
 * by arb, most of which rely on flint's ulong and slong defines.
 */

#include <mpfr.h>

#undef ulong
#undef slong

#define ulong mp_limb_t
#define slong mp_limb_signed_t

#include <flint/acb.h>
#include <flint/acb_calc.h>
#include <flint/acb_elliptic.h>
#include <flint/acb_hypgeom.h>
#include <flint/acb_mat.h>
#include <flint/acb_modular.h>
#include <flint/acb_poly.h>
#include <flint/arb.h>
#include <flint/arb_fmpz_poly.h>
#include <flint/arb_hypgeom.h>
#include <flint/arf.h>
#include <flint/bernoulli.h>
#include <flint/mag.h>

#undef ulong
#undef slong
#endif
