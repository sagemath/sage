"""
Eisenstein series, optimized
"""

from cysignals.memory cimport check_allocarray, sig_free
from cysignals.signals cimport sig_check

from sage.arith.misc import primes, bernoulli
from sage.rings.integer cimport Integer
from sage.rings.fast_arith cimport prime_range

from cpython.list cimport PyList_GET_ITEM
from sage.libs.flint.fmpz_poly cimport *
from sage.libs.gmp.mpz cimport *
from sage.libs.flint.fmpz_poly_sage cimport Fmpz_poly, fmpz_poly_set_coeff_mpz, fmpz_poly_scalar_mul_mpz

cpdef Ek_ZZ(int k, int prec=10):
    """
    Return list of prec integer coefficients of the weight k
    Eisenstein series of level 1, normalized so the coefficient of q
    is 1, except that the 0th coefficient is set to 1 instead of its
    actual value.

    INPUT:

    - ``k`` -- integer
    - ``prec`` -- integer

    OUTPUT: list of integers

    EXAMPLES::

        sage: from sage.modular.modform.eis_series_cython import Ek_ZZ
        sage: Ek_ZZ(4,10)
        [1, 1, 9, 28, 73, 126, 252, 344, 585, 757]
        sage: [sigma(n,3) for n in [1..9]]
        [1, 9, 28, 73, 126, 252, 344, 585, 757]
        sage: Ek_ZZ(10,10^3) == [1] + [sigma(n,9) for n in range(1,10^3)]
        True
    """
    cdef Integer pp
    cdef mpz_t q, current_sum, q_plus_one

    cdef unsigned long p
    cdef Py_ssize_t i, ind
    cdef bint continue_flag

    cdef list power_sum_ls

    cdef unsigned long max_power_sum, temp_index
    cdef unsigned long remainder, prev_index
    cdef unsigned long additional_p_powers

    mpz_init(q)
    mpz_init(current_sum)
    mpz_init(q_plus_one)

    # allocate the list for the result
    cdef list val = []
    for i in range(prec):
        pp = Integer.__new__(Integer)
        mpz_set_si(pp.value, 1)
        val.append(pp)

    # no need to reallocate this list every time -- just reuse the
    # Integers in it
    power_sum_ls = []
    max_power_sum = prec
    while max_power_sum:
        max_power_sum >>= 1
        pp = Integer.__new__(Integer)
        mpz_set_si(pp.value, 1)
        power_sum_ls.append(pp)

    for pp in prime_range(prec):
        p = mpz_get_ui((<Integer>pp).value)
        mpz_ui_pow_ui(q, p, k - 1)
        mpz_add_ui(q_plus_one, q, 1)
        mpz_set(current_sum, q_plus_one)

        # NOTE: I (wstein) did benchmarks, and the use of
        # PyList_GET_ITEM in the code below is worth it since it leads
        # to a significant speedup by not doing bounds checking.

        # set power_sum_ls[1] = q+1
        mpz_set((<Integer>(PyList_GET_ITEM(power_sum_ls, 1))).value,
                current_sum)
        max_power_sum = 1

        ind = 0
        while True:
            continue_flag = 0
            # do the first p-1
            for i in range(1, p):
                ind += p
                if (ind >= prec):
                    continue_flag = 1
                    break
                mpz_mul((<Integer>(PyList_GET_ITEM(val, ind))).value,
                        (<Integer>(PyList_GET_ITEM(val, ind))).value,
                        q_plus_one)
            ind += p
            if (ind >= prec or continue_flag):
                break

            # now do the pth one, which is harder.

            # compute the valuation of n at p
            additional_p_powers = 0
            temp_index = ind // p
            remainder = 0
            while not remainder:
                additional_p_powers += 1
                prev_index = temp_index
                temp_index = temp_index // p
                remainder = prev_index - p*temp_index

            # if we need a new sum, it has to be the next uncomputed one.
            if (additional_p_powers > max_power_sum):
                mpz_mul(current_sum, current_sum, q)
                mpz_add_ui(current_sum, current_sum, 1)
                max_power_sum += 1

                mpz_set((<Integer>(PyList_GET_ITEM(power_sum_ls,
                                                   max_power_sum))).value,
                        current_sum)

            # finally, set the coefficient
            mpz_mul((<Integer>(PyList_GET_ITEM(val, ind))).value,
                    (<Integer>(PyList_GET_ITEM(val, ind))).value,
                    (<Integer>(PyList_GET_ITEM(power_sum_ls,
                                               additional_p_powers))).value)

    mpz_clear(q)
    mpz_clear(current_sum)
    mpz_clear(q_plus_one)

    return val


cpdef eisenstein_series_poly(int k, int prec=10):
    r"""
    Return the `q`-expansion up to precision ``prec`` of the weight `k`
    Eisenstein series, as a FLINT :class:`~sage.libs.flint.fmpz_poly.Fmpz_poly`
    object, normalised so the coefficients are integers with no common factor.

    Used internally by the functions
    :func:`~sage.modular.modform.eis_series.eisenstein_series_qexp` and
    :func:`~sage.modular.modform.vm_basis.victor_miller_basis`; see the
    docstring of the former for further details.

    EXAMPLES::

        sage: from sage.modular.modform.eis_series_cython import eisenstein_series_poly
        sage: eisenstein_series_poly(12, prec=5)
        5  691 65520 134250480 11606736960 274945048560
    """
    cdef mpz_t *val = <mpz_t *>check_allocarray(prec, sizeof(mpz_t))
    cdef mpz_t one, mult, term, last, term_m1, last_m1
    cdef long ind
    cdef int i
    cdef Fmpz_poly res = Fmpz_poly.__new__(Fmpz_poly)

    if k % 2 or k < 2:
        raise ValueError("k (=%s) must be an even positive integer" % k)
    if prec < 0:
        raise ValueError("prec (=%s) must be an even nonnegative integer" % prec)
    if (prec == 0):
        return Fmpz_poly.__new__(Fmpz_poly)

    mpz_init(one)
    mpz_init(term)
    mpz_init(last)
    mpz_init(mult)
    mpz_init(term_m1)
    mpz_init(last_m1)

    for i in range(prec):
        mpz_init_set_si(val[i], 1)

    mpz_set_si(one, 1)

    cdef unsigned long expt = k - 1
    a0 = -bernoulli(k) / (2*k)

    cdef long p, ppow
    for p in primes(1, prec) :
        ppow = p

        mpz_set_si(mult, p)
        mpz_pow_ui(mult, mult, expt)
        mpz_mul(term, mult, mult)
        mpz_set(last, mult)

        while ppow < prec:
            sig_check()
            ind = ppow
            mpz_sub(term_m1, term, one)
            mpz_sub(last_m1, last, one)
            while ind < prec:
                mpz_mul(val[ind], val[ind], term_m1)
                mpz_fdiv_q(val[ind], val[ind], last_m1)
                ind += ppow
            ppow *= p
            mpz_set(last, term)
            mpz_mul(term, term, mult)

    mpz_clear(one)
    mpz_clear(term)
    mpz_clear(last)
    mpz_clear(mult)
    mpz_clear(term_m1)
    mpz_clear(last_m1)

    fmpz_poly_set_coeff_mpz(res.poly, prec-1, val[prec-1])
    for i in range(1, prec - 1):
        fmpz_poly_set_coeff_mpz(res.poly, i, val[i])

    fmpz_poly_scalar_mul_mpz(res.poly, res.poly, (<Integer>(a0.denominator())).value)
    fmpz_poly_set_coeff_mpz(res.poly, 0, (<Integer>(a0.numerator())).value)

    sig_free(val)

    return res
