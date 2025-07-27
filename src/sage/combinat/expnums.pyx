r"""
Bell and Uppuluri-Carpenter numbers

AUTHORS:

- Nick Alexander
"""

from cysignals.memory cimport check_allocarray, sig_free
from sage.libs.gmp.mpz cimport *
from sage.rings.integer cimport Integer


def expnums(int n, int aa):
    r"""
    Compute the first `n` exponential numbers around
    `aa`, starting with the zero-th.

    INPUT:

    - ``n`` -- C machine int

    - ``aa`` -- C machine int

    OUTPUT: list of length `n`

    ALGORITHM: We use the same integer addition algorithm as GAP. This
    is an extension of Bell's triangle to the general case of
    exponential numbers. The recursion performs `O(n^2)`
    additions, but the running time is dominated by the cost of the
    last integer addition, because the growth of the integer results of
    partial computations is exponential in `n`. The algorithm
    stores `O(n)` integers, but each is exponential in
    `n`.

    EXAMPLES::

        sage: expnums(10, 1)
        [1, 1, 2, 5, 15, 52, 203, 877, 4140, 21147]

    ::

        sage: expnums(10, -1)
        [1, -1, 0, 1, 1, -2, -9, -9, 50, 267]

    ::

        sage: expnums(1, 1)
        [1]
        sage: expnums(0, 1)
        []
        sage: expnums(-1, 0)
        []

    AUTHORS:

    - Nick Alexander
    """
    if n < 1:
        return []

    cdef Integer z
    if n == 1:
        z = Integer.__new__(Integer)
        mpz_set_si(z.value, 1)
        return [z]

    n = n - 1

    r = []
    z = Integer.__new__(Integer)
    mpz_set_si(z.value, 1)
    r.append(z)

    z = Integer.__new__(Integer)
    mpz_set_si(z.value, aa)
    r.append(z)

    cdef mpz_t *bell
    bell = <mpz_t *>check_allocarray(n + 1, sizeof(mpz_t))
    cdef mpz_t a
    mpz_init_set_si(a, aa)
    mpz_init_set_si(bell[1], aa)
    cdef int i
    cdef int k
    for i in range(1, n):
        mpz_init(bell[i+1])
        mpz_mul(bell[i+1], a, bell[1])
        for k in range(i):
            mpz_add(bell[i-k], bell[i-k], bell[i-k+1])

        z = Integer.__new__(Integer)
        mpz_set(z.value, bell[1])
        r.append(z)

    for i in range(1, n + 1):
        mpz_clear(bell[i])
    sig_free(bell)

    return r

# The following code is from GAP 4.4.9.
###################################################
# InstallGlobalFunction(Bell,function ( n )
#     local   bell, k, i;
#     bell := [ 1 ];
#     for i  in [1..n-1]  do
#         bell[i+1] := bell[1];
#         for k  in [0..i-1]  do
#             bell[i-k] := bell[i-k] + bell[i-k+1];
#         od;
#     od;
#     return bell[1];
# end);


def expnums2(n, aa):
    r"""
    A vanilla python (but compiled via Cython) implementation of
    expnums.

    We Compute the first `n` exponential numbers around
    `aa`, starting with the zero-th.

    EXAMPLES::

        sage: from sage.combinat.expnums import expnums2
        sage: expnums2(10, 1)
        [1, 1, 2, 5, 15, 52, 203, 877, 4140, 21147]
    """
    if n < 1:
        return []
    if n == 1:
        return [Integer(1)]

    n = n - 1
    r = [Integer(1), Integer(aa)]

    bell = [Integer(0)] * (n+1)
    a = aa
    bell[1] = aa
    for i in range(1, n):
        bell[i+1] = a * bell[1]
        for k in range(i):
            bell[i-k] = bell[i-k] + bell[i-k+1]
        r.append(bell[1])
    return r
