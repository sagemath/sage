r"""
Misc matrix algorithms using FLINT
"""

from cysignals.signals cimport sig_check

from sage.arith.rational_reconstruction cimport mpq_rational_reconstruction
from sage.libs.gmp.mpq cimport *
from sage.libs.gmp.mpz cimport *
from sage.libs.flint.fmpq cimport fmpq_set_mpq, fmpq_canonicalise
from sage.libs.flint.fmpq_mat cimport fmpq_mat_entry_num, fmpq_mat_entry_den, fmpq_mat_entry
from sage.libs.flint.fmpz cimport fmpz_set_mpz, fmpz_one
from sage.rings.integer cimport Integer
from sage.rings.rational_field import QQ

from sage.matrix.matrix_integer_dense cimport Matrix_integer_dense
from sage.matrix.matrix_rational_dense cimport Matrix_rational_dense


def matrix_integer_dense_rational_reconstruction(Matrix_integer_dense A, Integer N):
    r"""
    Given a matrix over the integers and an integer modulus, do
    rational reconstruction on all entries of the matrix, viewed as
    numbers mod `N`.  This is done efficiently by assuming there is a
    large common factor dividing the denominators.

    INPUT:

    - ``A`` -- matrix
    - ``N`` -- integer

    EXAMPLES::

        sage: B = ((matrix(ZZ, 3,4, [1,2,3,-4,7,2,18,3,4,3,4,5])/3)%500).change_ring(ZZ)
        sage: from sage.matrix.misc_flint import matrix_integer_dense_rational_reconstruction
        sage: matrix_integer_dense_rational_reconstruction(B, 500)
        [ 1/3  2/3    1 -4/3]
        [ 7/3  2/3    6    1]
        [ 4/3    1  4/3  5/3]

    TESTS:

    Check that :issue:`9345` is fixed::

        sage: A = random_matrix(ZZ, 3)
        sage: matrix_integer_dense_rational_reconstruction(A, 0)
        Traceback (most recent call last):
        ...
        ZeroDivisionError: The modulus cannot be zero
    """
    if not N:
        raise ZeroDivisionError("The modulus cannot be zero")
    cdef Matrix_rational_dense R
    R = Matrix_rational_dense.__new__(Matrix_rational_dense,
                                      A.parent().change_ring(QQ), 0,0,0)

    cdef mpz_t a, bnd, other_bnd, denom, tmp
    cdef mpq_t qtmp
    cdef Integer _bnd
    cdef Py_ssize_t i, j
    cdef int do_it

    mpz_init_set_si(denom, 1)
    mpz_init(a)
    mpz_init(tmp)
    mpz_init(other_bnd)
    mpq_init(qtmp)

    _bnd = (N//2).isqrt()
    mpz_init_set(bnd, _bnd.value)
    mpz_sub(other_bnd, N.value, bnd)

    for i in range(A._nrows):
        for j in range(A._ncols):
            sig_check()
            A.get_unsafe_mpz(i, j, a)
            if mpz_cmp_ui(denom, 1) != 0:
                mpz_mul(a, a, denom)
            mpz_fdiv_r(a, a, N.value)
            do_it = 0
            if mpz_cmp(a, bnd) <= 0:
                do_it = 1
            elif mpz_cmp(a, other_bnd) >= 0:
                mpz_sub(a, a, N.value)
                do_it = 1
            if do_it:
                fmpz_set_mpz(fmpq_mat_entry_num(R._matrix, i, j), a)
                if mpz_cmp_ui(denom, 1) != 0:
                    fmpz_set_mpz(fmpq_mat_entry_den(R._matrix, i, j), denom)
                    fmpq_canonicalise(fmpq_mat_entry(R._matrix, i, j))
                else:
                    fmpz_one(fmpq_mat_entry_den(R._matrix, i, j))
            else:
                # Otherwise have to do it the hard way
                A.get_unsafe_mpz(i, j, tmp)
                mpq_rational_reconstruction(qtmp, tmp, N.value)
                mpz_lcm(denom, denom, mpq_denref(qtmp))
                fmpq_set_mpq(fmpq_mat_entry(R._matrix, i, j), qtmp)

    mpz_clear(denom)
    mpz_clear(a)
    mpz_clear(tmp)
    mpz_clear(other_bnd)
    mpz_clear(bnd)
    mpq_clear(qtmp)

    return R
