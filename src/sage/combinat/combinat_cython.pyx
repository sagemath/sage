"""
Fast computation of combinatorial functions (Cython + mpz).

Currently implemented:
- Stirling numbers of the second kind

AUTHORS:
- Fredrik Johansson (2010-10): Stirling numbers of second kind

"""

include "sage/ext/stdsage.pxi"


from sage.libs.gmp.all cimport *
from sage.rings.integer cimport Integer

cdef void mpz_addmul_alt(mpz_t s, mpz_t t, mpz_t u, unsigned long parity):
    """
    Set s = s + t*u * (-1)^parity
    """
    if parity & 1:
        mpz_submul(s, t, u)
    else:
        mpz_addmul(s, t, u)


cdef mpz_stirling_s2(mpz_t s, unsigned long n, unsigned long k):
    """
    Set s = S(n,k) where S(n,k) denotes a Stirling number of the
    second kind.

    Algorithm: S(n,k) = (sum_{j=0}^k (-1)^(k-j) C(k,j) j^n) / k!

    TODO: compute S(n,k) efficiently for large n when n-k is small
    (e.g. when k > 20 and n-k < 20)
    """
    cdef mpz_t t, u
    cdef mpz_t *bc
    cdef unsigned long j, max_bc
    # Some important special cases
    if k+1 >= n:
        # Upper triangle of n\k table
        if k > n:
            mpz_set_ui(s, 0)
        elif n == k:
            mpz_set_ui(s, 1)
        elif k+1 == n:
            # S(n,n-1) = C(n,2)
            mpz_set_ui(s, n)
            mpz_mul_ui(s, s, n-1)
            mpz_tdiv_q_2exp(s, s, 1)
    elif k <= 2:
        # Leftmost three columns of n\k table
        if k == 0:
            mpz_set_ui(s, 0)
        elif k == 1:
            mpz_set_ui(s, 1)
        elif k == 2:
            # 2^(n-1)-1
            mpz_set_ui(s, 1)
            mpz_mul_2exp(s, s, n-1)
            mpz_sub_ui(s, s, 1)
    # Direct sequential evaluation of the sum
    elif n < 200:
        mpz_init(t)
        mpz_init(u)
        mpz_set_ui(t, 1)
        mpz_set_ui(s, 0)
        for j in range(1, k//2+1):
            mpz_mul_ui(t, t, k+1-j)
            mpz_tdiv_q_ui(t, t, j)
            mpz_set_ui(u, j)
            mpz_pow_ui(u, u, n)
            mpz_addmul_alt(s, t, u, k+j)
            if 2*j != k:
                # Use the fact that C(k,j) = C(k,k-j)
                mpz_set_ui(u, k-j)
                mpz_pow_ui(u, u, n)
                mpz_addmul_alt(s, t, u, j)
        # Last term not included because loop starts from 1
        mpz_set_ui(u, k)
        mpz_pow_ui(u, u, n)
        mpz_add(s, s, u)
        mpz_fac_ui(t, k)
        mpz_tdiv_q(s, s, t)
        mpz_clear(t)
        mpz_clear(u)
    # Only compute odd powers, saving about half of the time for large n.
    # We need to precompute binomial coefficients since they will be accessed
    # out of order, adding overhead that makes this slower for small n.
    else:
        mpz_init(t)
        mpz_init(u)
        max_bc = (k+1)//2
        bc = <mpz_t*> sage_malloc((max_bc+1) * sizeof(mpz_t))
        mpz_init_set_ui(bc[0], 1)
        for j in range(1, max_bc+1):
            mpz_init_set(bc[j], bc[j-1])
            mpz_mul_ui(bc[j], bc[j], k+1-j)
            mpz_tdiv_q_ui(bc[j], bc[j], j)
        mpz_set_ui(s, 0)
        for j in range(1, k+1, 2):
            mpz_set_ui(u, j)
            mpz_pow_ui(u, u, n)
            # Process each 2^p * j, where j is odd
            while True:
                if j > max_bc:
                    mpz_addmul_alt(s, bc[k-j], u, k+j)
                else:
                    mpz_addmul_alt(s, bc[j], u, k+j)
                j *= 2
                if j > k:
                    break
                mpz_mul_2exp(u, u, n)
        for j in range(max_bc+1):   # careful: 0 ... max_bc
            mpz_clear(bc[j])
        sage_free(bc)
        mpz_fac_ui(t, k)
        mpz_tdiv_q(s, s, t)
        mpz_clear(t)
        mpz_clear(u)

def _stirling_number2(n, k):
    """
    Python wrapper of mpz_stirling_s2.

        sage: from sage.combinat.combinat_cython import _stirling_number2
        sage: _stirling_number2(3, 2)
        3

    This is wrapped again by stirling_number2 in combinat.py.
    """
    cdef Integer s
    s = PY_NEW(Integer)
    mpz_stirling_s2(s.value, n, k)
    return s
