#############################################################
#
#    Sparse Vector over mpz_t (the GMP integers)
#
#############################################################

from sage.libs.gmp.types cimport mpz_t

cdef struct mpz_vector:
    mpz_t *entries      # array of nonzero entries
    Py_ssize_t   *positions    # positions of those nonzero entries, starting at 0
    Py_ssize_t    degree       # the degree of this sparse vector
    Py_ssize_t    num_nonzero  # the number of nonzero entries of this vector.

