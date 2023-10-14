# distutils: libraries = flint
# distutils: depends = flint/fmpz_factor.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void fmpz_factor_init(fmpz_factor_t factor)
    # Initialises an ``fmpz_factor_t`` structure.

    void fmpz_factor_clear(fmpz_factor_t factor)
    # Clears an ``fmpz_factor_t`` structure.

    void _fmpz_factor_append_ui(fmpz_factor_t factor, mp_limb_t p, ulong exp)
    # Append a factor `p` to the given exponent to the
    # ``fmpz_factor_t`` structure ``factor``.

    void _fmpz_factor_append(fmpz_factor_t factor, const fmpz_t p, ulong exp)
    # Append a factor `p` to the given exponent to the
    # ``fmpz_factor_t`` structure ``factor``.

    void fmpz_factor(fmpz_factor_t factor, const fmpz_t n)
    # Factors `n` into prime numbers. If `n` is zero or negative, the
    # sign field of the ``factor`` object will be set accordingly.

    int fmpz_factor_smooth(fmpz_factor_t factor, const fmpz_t n, slong bits, int proved)
    # Factors `n` into prime numbers up to approximately the given number of
    # bits and possibly one additional cofactor, which may or may not be prime.
    # If the number is definitely factored fully, the return value is `1`,
    # otherwise the final factor (which may have exponent greater than `1`)
    # is composite and needs to be factored further.
    # If the number has a factor of around the given number of bits, there is
    # at least a two-thirds chance of finding it. Smaller factors should be
    # found with a much higher probability.
    # The amount of time spent factoring can be controlled by lowering or
    # increasing ``bits``. However, the quadratic sieve may be faster if
    # ``bits`` is set to more than one third of the number of bits of `n`.
    # The function uses trial factoring up to ``bits = 15``, followed by
    # a primality test and a perfect power test to check if the factorisation
    # is complete. If ``bits`` is at least 16, it proceeds to use the
    # elliptic curve method to look for larger factors.
    # The behavior of primality testing is determined by the ``proved``
    # parameter:
    # If ``proved`` is set to `1` the function will prove all factors prime
    # (other than the last factor, if the return value is `0`).
    # If ``proved`` is set to `0`, the function will only check that factors are
    # probable primes.
    # If ``proved`` is set to `-1`, the function will not test primality
    # after performing trial division. A perfect power test is still performed.
    # As an exception to the rules stated above, this function will call
    # ``n_factor`` internally if `n` or the remainder after trial division
    # is smaller than one word, guaranteeing a complete factorisation.

    void fmpz_factor_si(fmpz_factor_t factor, slong n)
    # Like ``fmpz_factor``, but takes a machine integer `n` as input.

    int fmpz_factor_trial_range(fmpz_factor_t factor, const fmpz_t n, ulong start, ulong num_primes)
    # Factors `n` into prime factors using trial division. If `n` is
    # zero or negative, the sign field of the ``factor`` object will be
    # set accordingly.
    # The algorithm starts with the given start index in the ``flint_primes``
    # table and uses at most ``num_primes`` primes from that point.
    # The function returns 1 if `n` is completely factored, otherwise it returns
    # `0`.

    int fmpz_factor_trial(fmpz_factor_t factor, const fmpz_t n, slong num_primes)
    # Factors `n` into prime factors using trial division. If `n` is
    # zero or negative, the sign field of the ``factor`` object will be
    # set accordingly.
    # The algorithm uses the given number of primes, which must be in the range
    # `[0, 3512]`. An exception is raised if a number outside this range is
    # passed.
    # The function returns 1 if `n` is completely factored, otherwise it returns
    # `0`.
    # The final entry in the factor struct is set to the cofactor after removing
    # prime factors, if this is not `1`.

    void fmpz_factor_refine(fmpz_factor_t res, const fmpz_factor_t f)
    # Attempts to improve a partial factorization of an integer by "refining"
    # the factorization ``f`` to a more complete factorization ``res``
    # whose bases are pairwise relatively prime.
    # This function does not require its input to be in canonical form,
    # nor does it guarantee that the resulting factorization will be canonical.

    void fmpz_factor_expand_iterative(fmpz_t n, const fmpz_factor_t factor)
    # Evaluates an integer in factored form back to an ``fmpz_t``.
    # This currently exponentiates the bases separately and multiplies
    # them together one by one, although much more efficient algorithms
    # exist.

    int fmpz_factor_pp1(fmpz_t factor, const fmpz_t n, ulong B1, ulong B2_sqrt, ulong c)
    # Use Williams' `p + 1` method to factor `n`, using a prime bound in
    # stage 1 of ``B1`` and a prime limit in stage 2 of at least the square
    # of ``B2_sqrt``. If a factor is found, the function returns `1` and
    # ``factor`` is set to the factor that is found. Otherwise, the function
    # returns `0`.
    # The value `c` should be a random value greater than `2`. Successive
    # calls to the function with different values of `c` give additional
    # chances to factor `n` with roughly exponentially decaying probability
    # of finding a factor which has been missed (if `p+1` or `p-1` is not
    # smooth for any prime factors `p` of `n` then the function will
    # not ever succeed).

    int fmpz_factor_pollard_brent_single(fmpz_t p_factor, fmpz_t n_in, fmpz_t yi, fmpz_t ai, mp_limb_t max_iters)
    # Pollard Rho algorithm for integer factorization. Assumes that the `n` is
    # not prime. ``factor`` is set as the factor if found. Takes as input the initial
    # value `y`, to start polynomial evaluation, and `a`, the constant of the polynomial
    # used. It is not assured that the factor found will be prime. Does not compute
    # the complete factorization, just one factor. Returns the number of limbs of
    # factor if factorization is successful (non trivial factor is found), else returns 0.
    # ``max_iters`` is the number of iterations tried in process of finding the cycle.
    # If the algorithm fails to find a non trivial factor in one call, it tries again
    # (this time with a different set of random values).

    int fmpz_factor_pollard_brent(fmpz_t factor, flint_rand_t state, fmpz_t n, mp_limb_t max_tries, mp_limb_t max_iters)
    # Pollard Rho algorithm for integer factorization. Assumes that the `n` is
    # not prime. ``factor`` is set as the factor if found. It is not assured that the
    # factor found will be prime. Does not compute the complete factorization,
    # just one factor. Returns the number of limbs of factor if factorization is
    # successful (non trivial factor is found), else returns 0.
    # ``max_iters`` is the number of iterations tried in process of finding the cycle.
    # If the algorithm fails to find a non trivial factor in one call, it tries again
    # (this time with a different set of random values). This process is repeated a
    # maximum of ``max_tries`` times.
    # The algorithm used is a modification of the original Pollard Rho algorithm,
    # suggested by Richard Brent. It can be found in the paper available at
    # https://maths-people.anu.edu.au/~brent/pd/rpb051i.pdf

    void fmpz_factor_ecm_init(ecm_t ecm_inf, mp_limb_t sz)
    # Initializes the ``ecm_t`` struct. This is needed in some functions
    # and carries data between subsequent calls.

    void fmpz_factor_ecm_clear(ecm_t ecm_inf)
    # Clears the ``ecm_t`` struct.

    void fmpz_factor_ecm_addmod(mp_ptr a, mp_ptr b, mp_ptr c, mp_ptr n, mp_limb_t n_size)
    # Sets `a` to `(b + c)` ``%`` `n`. This is not a normal add mod function,
    # it assumes `n` is normalized (highest bit set) and `b` and `c` are reduced
    # modulo `n`.
    # Used for arithmetic operations in ``fmpz_factor_ecm``.

    void fmpz_factor_ecm_submod(mp_ptr x, mp_ptr a, mp_ptr b, mp_ptr n, mp_limb_t n_size)
    # Sets `x` to `(a - b)` ``%`` `n`. This is not a normal subtract mod
    # function, it assumes `n` is normalized (highest bit set)
    # and `b` and `c` are reduced modulo `n`.
    # Used for arithmetic operations in ``fmpz_factor_ecm``.

    void fmpz_factor_ecm_double(mp_ptr x, mp_ptr z, mp_ptr x0, mp_ptr z0, mp_ptr n, ecm_t ecm_inf)
    # Sets the point `(x : z)` to two times `(x_0 : z_0)` modulo `n` according
    # to the formula
    # .. math ::
    # x = (x_0 + z_0)^2 \cdot (x_0 - z_0)^2 \mod n,
    # .. math ::
    # z = 4 x_0 z_0 \left((x_0 - z_0)^2 + 4a_{24}x_0z_0\right) \mod n.
    # ``ecm_inf`` is used just to use temporary ``mp_ptr``'s in the
    # structure. This group doubling is valid only for points expressed in
    # Montgomery projective coordinates.

    void fmpz_factor_ecm_add(mp_ptr x, mp_ptr z, mp_ptr x1, mp_ptr z1, mp_ptr x2, mp_ptr z2, mp_ptr x0, mp_ptr z0, mp_ptr n, ecm_t ecm_inf)
    # Sets the point `(x : z)` to the sum of `(x_1 : z_1)` and `(x_2 : z_2)`
    # modulo `n`, given the difference `(x_0 : z_0)` according to the formula
    # .. math ::
    # x = 4z_0(x_1x_2 - z_1z_2)^2 \mod n, \\ z = 4x_0(x_2z_1 - x_1z_2)^2 \mod n.
    # ``ecm_inf`` is used just to use temporary ``mp_ptr``'s in the
    # structure. This group addition is valid only for points expressed in
    # Montgomery projective coordinates.

    void fmpz_factor_ecm_mul_montgomery_ladder(mp_ptr x, mp_ptr z, mp_ptr x0, mp_ptr z0, mp_limb_t k, mp_ptr n, ecm_t ecm_inf)
    # Montgomery ladder algorithm for scalar multiplication of elliptic points.
    # Sets the point `(x : z)` to `k(x_0 : z_0)` modulo `n`.
    # ``ecm_inf`` is used just to use temporary ``mp_ptr``'s in the
    # structure. Valid only for points expressed in Montgomery projective
    # coordinates.

    int fmpz_factor_ecm_select_curve(mp_ptr f, mp_ptr sigma, mp_ptr n, ecm_t ecm_inf)
    # Selects a random elliptic curve given a random integer ``sigma``,
    # according to Suyama's parameterization. If the factor is found while
    # selecting the curve, the number of limbs required to store the factor
    # is returned, otherwise `0`.
    # It could be possible that the selected curve is unsuitable for further
    # computations, in such a case, `-1` is returned.
    # Also selects the initial point `x_0`, and the value of `(a + 2)/4`, where `a`
    # is a curve parameter. Sets `z_0` as `1`. All these are stored in the
    # ``ecm_t`` struct.
    # The curve selected is of Montgomery form, the points selected satisfy the
    # curve and are projective coordinates.

    int fmpz_factor_ecm_stage_I(mp_ptr f, const mp_limb_t *prime_array, mp_limb_t num, mp_limb_t B1, mp_ptr n, ecm_t ecm_inf)
    # Stage I implementation of the ECM algorithm.
    # ``f`` is set as the factor if found. ``num`` is number of prime numbers
    # `\le` the bound ``B1``. ``prime_array`` is an array of first ``B1``
    # primes. `n` is the number being factored.
    # If the factor is found, number of words required to store the factor is
    # returned, otherwise `0`.

    int fmpz_factor_ecm_stage_II(mp_ptr f, mp_limb_t B1, mp_limb_t B2, mp_limb_t P, mp_ptr n, ecm_t ecm_inf)
    # Stage II implementation of the ECM algorithm.
    # ``f`` is set as the factor if found. ``B1``, ``B2`` are the two
    # bounds. ``P`` is the primorial (approximately equal to `\sqrt{B2}`).
    # `n` is the number being factored.
    # If the factor is found, number of words required to store the factor is
    # returned, otherwise `0`.

    int fmpz_factor_ecm(fmpz_t f, mp_limb_t curves, mp_limb_t B1, mp_limb_t B2, flint_rand_t state, const fmpz_t n_in)
    # Outer wrapper function for the ECM algorithm. In case ``f`` can fit
    # in a single unsigned word, a call to ``n_factor_ecm`` is made.
    # The function calls stage I and II, and
    # the precomputations (builds ``prime_array`` for stage I,
    # ``GCD_table`` and ``prime_table`` for stage II).
    # ``f`` is set as the factor if found. ``curves`` is the number of
    # random curves being tried. ``B1``, ``B2`` are the two bounds or
    # stage I and stage II. `n` is the number being factored.
    # If a factor is found in stage I, `1` is returned.
    # If a factor is found in stage II, `2` is returned.
    # If a factor is found while selecting the curve, `-1` is returned.
    # Otherwise `0` is returned.
