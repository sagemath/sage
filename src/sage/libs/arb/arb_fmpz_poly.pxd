# sage_setup: distribution = sagemath-flint

# Deprecated header file; use sage/libs/flint/arb_fmpz_poly.pxd instead
# See https://github.com/sagemath/sage/pull/36449

from sage.libs.flint.arb_fmpz_poly cimport (
    _arb_fmpz_poly_evaluate_arb_horner,
    arb_fmpz_poly_evaluate_arb_horner,
    _arb_fmpz_poly_evaluate_arb_rectangular,
    arb_fmpz_poly_evaluate_arb_rectangular,
    _arb_fmpz_poly_evaluate_arb,
    arb_fmpz_poly_evaluate_arb,
    _arb_fmpz_poly_evaluate_acb_horner,
    arb_fmpz_poly_evaluate_acb_horner,
    _arb_fmpz_poly_evaluate_acb_rectangular,
    arb_fmpz_poly_evaluate_acb_rectangular,
    _arb_fmpz_poly_evaluate_acb,
    arb_fmpz_poly_evaluate_acb,
    arb_fmpz_poly_deflation,
    arb_fmpz_poly_deflate,
    arb_fmpz_poly_complex_roots,
    arb_fmpz_poly_gauss_period_minpoly)
