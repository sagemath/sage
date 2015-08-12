#############################################################
#
#    Sparse Vector over mpq_t (the GMP rationals)
#
#############################################################

from sage.libs.gmp.types cimport mpq_t

cdef struct mpq_vector:
    mpq_t *entries      # array of nonzero entries
    Py_ssize_t   *positions    # positions of those nonzero entries, starting at 0
    Py_ssize_t    degree       # the degree of this sparse vector
    Py_ssize_t    num_nonzero  # the number of nonzero entries of this vector.

