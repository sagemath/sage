# Deprecated header file; use sage/libs/flint/acb_modular.pxd instead
# See https://github.com/sagemath/sage/pull/36449

from sage.libs.flint.acb_modular cimport (
    acb_modular_theta,
    acb_modular_j,
    acb_modular_eta,
    acb_modular_lambda,
    acb_modular_delta,
    acb_modular_eisenstein,
    acb_modular_elliptic_p,
    acb_modular_elliptic_p_zpx,
    acb_modular_elliptic_k,
    acb_modular_elliptic_k_cpx,
    acb_modular_elliptic_e,
    acb_modular_hilbert_class_poly)
