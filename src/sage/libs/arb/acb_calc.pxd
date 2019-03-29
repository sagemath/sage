# distutils: libraries = gmp flint arb

from sage.libs.arb.types cimport *
from sage.libs.flint.types cimport fmpz_t, fmpq_t

cdef extern from "acb_calc.h":

    void acb_calc_integrate_opt_init(acb_calc_integrate_opt_t options)

    int acb_calc_integrate(
             acb_t res,
             acb_calc_func_t func,
             void * param,
             const acb_t a,
             const acb_t b,
             long rel_goal,
             const mag_t abs_tol,
             const acb_calc_integrate_opt_t options,
             long prec)
