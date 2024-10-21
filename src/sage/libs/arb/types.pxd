# sage_setup: distribution = sagemath-flint

# Deprecated header file; use sage/libs/flint/types.pxd instead
# See https://github.com/sagemath/sage/pull/36449

from sage.libs.flint.types cimport (
    mag_struct,
    mag_t,
    mag_ptr,
    mag_srcptr,
    arf_struct,
    arf_t,
    arf_ptr,
    arf_srcptr,
    arf_rnd_t,
    ARF_RND_DOWN,
    ARF_RND_UP,
    ARF_RND_FLOOR,
    ARF_RND_CEIL,
    ARF_RND_NEAR,
    ARF_PREC_EXACT,
    arb_struct,
    arb_t,
    arb_ptr,
    arb_srcptr,
    acb_struct,
    acb_t,
    acb_ptr,
    acb_srcptr,
    acb_mat_struct,
    acb_mat_t,
    acb_poly_struct,
    acb_poly_t,
    acb_calc_integrate_opt_struct,
    acb_calc_integrate_opt_t,
    acb_calc_func_t,
    arb_poly_struct,
    arb_poly_t)
