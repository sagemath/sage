#ifndef SAGE_ARB_WRAP_H
#define SAGE_ARB_WRAP_H
/*
 * Similar to flint_wrap.h but specifically for wrapping the headers supplied
 * by arb, most of which rely on flint's ulong and slong defines.
 */

#undef ulong
#undef slong

#define ulong mp_limb_t
#define slong mp_limb_signed_t

#include <acb.h>
#include <acb_calc.h>
#include <acb_elliptic.h>
#include <acb_hypgeom.h>
#include <acb_mat.h>
#include <acb_modular.h>
#include <acb_poly.h>
#include <arb.h>
#include <arb_fmpz_poly.h>
#include <arb_hypgeom.h>
#include <arf.h>
#include <bernoulli.h>
#include <mag.h>

#undef ulong
#undef slong
#endif
