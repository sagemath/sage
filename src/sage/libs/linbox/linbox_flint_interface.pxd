# sage_setup: distribution = sagemath-linbox
# distutils: libraries = LINBOX_LIBRARIES
# distutils: library_dirs = LINBOX_LIBDIR
# distutils: extra_link_args = LINBOX_LIBEXTRA

from sage.libs.flint.types cimport fmpz_t, fmpz_mat_t, fmpz_poly_t

# set C <- A * B
cdef void linbox_fmpz_mat_mul(fmpz_mat_t C, fmpz_mat_t A, fmpz_mat_t B) noexcept

# set cp to the characteristic polynomial of A
cdef void linbox_fmpz_mat_charpoly(fmpz_poly_t cp, fmpz_mat_t A) noexcept

# set mp to the minimal polynomial of A
cdef void linbox_fmpz_mat_minpoly(fmpz_poly_t mp, fmpz_mat_t A) noexcept

# return the rank of A
cdef size_t linbox_fmpz_mat_rank(fmpz_mat_t A) noexcept

# set det to the determinant of A
cdef void linbox_fmpz_mat_det(fmpz_t det, fmpz_mat_t A) noexcept
