# sage_setup: distribution = sagemath-flint

# Deprecated header file; use sage/libs/flint/acb_calc.pxd instead
# See https://github.com/sagemath/sage/pull/36449

from sage.libs.flint.types cimport acb_calc_integrate_opt_t, acb_calc_func_t

from sage.libs.flint.acb_calc cimport (
    acb_calc_integrate,
    acb_calc_integrate_opt_init)
