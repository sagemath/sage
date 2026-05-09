from sage.algebras.fusion_rings.shm_managers cimport KSHandler
from sage.rings.number_field.number_field_element cimport NumberFieldElement_absolute
from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomial_libsingular, MPolynomialRing_libsingular
from sage.rings.polynomial.polydict cimport ETuple

cpdef tuple poly_to_tup(MPolynomial_libsingular poly)
cpdef MPolynomial_libsingular _tup_to_poly(tuple eq_tup, MPolynomialRing_libsingular parent)
cpdef tuple resize(tuple eq_tup, dict idx_map, int nvars)
cpdef list get_variables_degrees(list eqns, int nvars)
cpdef list variables(tuple eq_tup)
cpdef constant_coeff(tuple eq_tup, field)
cpdef tuple apply_coeff_map(tuple eq_tup, coeff_map)
# cpdef bint tup_fixes_sq(tuple eq_tup)
cdef bint tup_fixes_sq(tuple eq_tup) noexcept
cdef dict subs_squares(dict eq_dict, KSHandler known_sq)
cpdef dict compute_known_powers(max_degs, dict val_dict, one)
cdef dict subs(tuple poly_tup, dict known_powers, one)
cpdef tup_to_univ_poly(tuple eq_tup, univ_poly_ring)
cpdef tuple poly_tup_sortkey(tuple eq_tup)
cdef tuple reduce_poly_dict(dict eq_dict, ETuple nonz, KSHandler known_sq, NumberFieldElement_absolute one)
cdef tuple _flatten_coeffs(tuple eq_tup)
cpdef tuple _unflatten_coeffs(field, tuple eq_tup)
cdef int has_appropriate_linear_term(tuple eq_tup) noexcept
