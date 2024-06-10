r"""
Misc matrix algorithms using MPFR
"""

from cysignals.signals cimport sig_check

cimport sage.rings.abc

from sage.libs.mpfr cimport *
from sage.rings.real_mpfr cimport RealNumber

from sage.matrix.matrix0 cimport Matrix


def hadamard_row_bound_mpfr(Matrix A):
    r"""
    Given a matrix `A` with entries that coerce to ``RR``, compute the row
    Hadamard bound on the determinant.

    INPUT:

    - ``A`` -- a matrix over ``RR``

    OUTPUT:

    integer -- an integer n such that the absolute value of the
    determinant of this matrix is at most `10^n`.

    EXAMPLES:

    We create a very large matrix, compute the row Hadamard bound,
    and also compute the row Hadamard bound of the transpose, which
    happens to be sharp. ::

        sage: a = matrix(ZZ, 2, [2^10000, 3^10000, 2^50, 3^19292])
        sage: from sage.matrix.misc_mpfr import hadamard_row_bound_mpfr
        sage: hadamard_row_bound_mpfr(a.change_ring(RR))
        13976
        sage: len(str(a.det()))
        12215
        sage: hadamard_row_bound_mpfr(a.transpose().change_ring(RR))
        12215

    Note that in the above example using RDF would overflow::

        sage: b = a.change_ring(RDF)
        sage: b._hadamard_row_bound()
        Traceback (most recent call last):
        ...
        OverflowError: cannot convert float infinity to integer
    """
    if not isinstance(A.base_ring(), sage.rings.abc.RealField):
        raise TypeError("A must have base field an mpfr real field.")

    cdef RealNumber a, b
    cdef mpfr_t s, d, pr
    cdef Py_ssize_t i, j

    mpfr_init(s)
    mpfr_init(d)
    mpfr_init(pr)
    mpfr_set_si(d, 0, MPFR_RNDU)

    for i in range(A._nrows):
        mpfr_set_si(s, 0, MPFR_RNDU)
        for j in range(A._ncols):
            sig_check()
            a = A.get_unsafe(i, j)
            mpfr_mul(pr, a.value, a.value, MPFR_RNDU)
            mpfr_add(s, s, pr, MPFR_RNDU)
        mpfr_log10(s, s, MPFR_RNDU)
        mpfr_add(d, d, s, MPFR_RNDU)
    b = a._new()
    mpfr_set(b.value, d, MPFR_RNDU)
    b /= 2
    mpfr_clear(s)
    mpfr_clear(d)
    mpfr_clear(pr)
    return b.ceil()
