# sage_setup: distribution = sagemath-flint

# Deprecated header file; use sage/libs/flint/acb_elliptic.pxd instead
# See https://github.com/sagemath/sage/pull/36449

from sage.libs.flint.acb_elliptic cimport (
    acb_elliptic_k,
    acb_elliptic_k_jet,
    acb_elliptic_k_series,
    acb_elliptic_e,
    acb_elliptic_rf,
    acb_elliptic_rj,
    acb_elliptic_rg,
    acb_elliptic_rc1,
    acb_elliptic_f,
    acb_elliptic_e_inc,
    acb_elliptic_pi,
    acb_elliptic_pi_inc,
    acb_elliptic_p,
    acb_elliptic_p_jet,
    acb_elliptic_p_series,
    acb_elliptic_zeta,
    acb_elliptic_sigma,
    acb_elliptic_roots,
    acb_elliptic_invariants,
    acb_elliptic_inv_p)
