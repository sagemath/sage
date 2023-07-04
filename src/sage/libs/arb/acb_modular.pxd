# distutils: libraries = gmp flint ARB_LIBRARY
# distutils: depends = acb_modular.h

from sage.libs.arb.types cimport *
from sage.libs.flint.types cimport fmpz_poly_t

# acb_modular.h
cdef extern from "arb_wrap.h":
    void acb_modular_theta(acb_t theta1, acb_t theta2, acb_t theta3, acb_t theta4, const acb_t z, const acb_t tau, long prec)
    void acb_modular_j(acb_t z, const acb_t tau, long prec)
    void acb_modular_eta(acb_t z, const acb_t tau, long prec)
    void acb_modular_lambda(acb_t r, const acb_t tau, long prec)
    void acb_modular_delta(acb_t r, const acb_t tau, long prec)
    void acb_modular_eisenstein(acb_ptr r, const acb_t tau, long len, long prec)
    void acb_modular_elliptic_p(acb_t r, const acb_t z, const acb_t tau, long prec)
    void acb_modular_elliptic_p_zpx(acb_ptr r, const acb_t z, const acb_t tau, long len, long prec)
    void acb_modular_elliptic_k(acb_t k, const acb_t m, long prec)
    void acb_modular_elliptic_k_cpx(acb_ptr w, const acb_t m, long len, long prec)
    void acb_modular_elliptic_e(acb_t res, const acb_t m, long prec)
    void acb_modular_hilbert_class_poly(fmpz_poly_t res, long D)

