# distutils: libraries = flint
# distutils: depends = flint/ulong_extras.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    ulong n_randlimb(flint_rand_t state)
    # Returns a uniformly pseudo random limb.
    # The algorithm generates two random half limbs `s_j`, `j = 0, 1`,
    # by iterating respectively `v_{i+1} = (v_i a + b) \bmod{p_j}` for
    # some initial seed `v_0`, randomly chosen values `a` and `b` and
    # ``p_0 = 4294967311 = nextprime(2^32)`` on a 64-bit machine
    # and ``p_0 = nextprime(2^16)`` on a 32-bit machine and
    # ``p_1 = nextprime(p_0)``.

    ulong n_randbits(flint_rand_t state, unsigned int bits)
    # Returns a uniformly pseudo random number with the given number of
    # bits. The most significant bit is always set, unless zero is passed,
    # in which case zero is returned.

    ulong n_randtest_bits(flint_rand_t state, int bits)
    # Returns a uniformly pseudo random number with the given number of
    # bits. The most significant bit is always set, unless zero is passed,
    # in which case zero is returned. The probability of a value with a
    # sparse binary representation being returned is increased. This
    # function is intended for use in test code.

    ulong n_randint(flint_rand_t state, ulong limit)
    # Returns a uniformly pseudo random number up to but not including
    # the given limit. If zero is passed as a parameter, an entire random
    # limb is returned.

    ulong n_urandint(flint_rand_t state, ulong limit)
    # Returns a uniformly pseudo random number up to but not including
    # the given limit. If zero is passed as a parameter, an entire
    # random limb is returned. This function provides somewhat better
    # randomness as compared to :func:`n_randint`, especially for larger
    # values of limit.

    ulong n_randtest(flint_rand_t state)
    # Returns a pseudo random number with a random number of bits,
    # from `0` to ``FLINT_BITS``.  The probability of the special
    # values `0`, `1`, ``COEFF_MAX`` and ``WORD_MAX`` is increased
    # as is the probability of a value with sparse binary representation.
    # This random function is mainly used for testing purposes.
    # This function is intended for use in test code.

    ulong n_randtest_not_zero(flint_rand_t state)
    # As for :func:`n_randtest`, but does not return `0`.
    # This function is intended for use in test code.

    ulong n_randprime(flint_rand_t state, ulong bits, int proved)
    # Returns a random prime number ``(proved = 1)`` or probable prime
    # ``(proved = 0)``
    # with ``bits`` bits, where ``bits`` must be at least 2 and
    # at most ``FLINT_BITS``.

    ulong n_randtest_prime(flint_rand_t state, int proved)
    # Returns a random prime number ``(proved = 1)`` or probable
    # prime ``(proved = 0)``
    # with size randomly chosen between 2 and ``FLINT_BITS`` bits.
    # This function is intended for use in test code.

    ulong n_pow(ulong n, ulong exp)
    # Returns ``n^exp``. No checking is done for overflow. The exponent
    # may be zero. We define `0^0 = 1`.
    # The algorithm simply uses a for loop. Repeated squaring is
    # unlikely to speed up this algorithm.

    ulong n_flog(ulong n, ulong b)
    # Returns `\lfloor\log_b n\rfloor`.
    # Assumes that `n \geq 1` and `b \geq 2`.

    ulong n_clog(ulong n, ulong b)
    # Returns `\lceil\log_b n\rceil`.
    # Assumes that `n \geq 1` and `b \geq 2`.

    ulong n_clog_2exp(ulong n, ulong b)
    # Returns `\lceil\log_b 2^n\rceil`.
    # Assumes that `b \geq 2`.

    ulong n_revbin(ulong n, ulong b)
    # Returns the binary reverse of `n`, assuming it is `b` bits in length,
    # e.g. ``n_revbin(10110, 6)`` will return ``110100``.

    int n_sizeinbase(ulong n, int base)
    # Returns the exact number of digits needed to represent `n` as a
    # string in base ``base`` assumed to be between 2 and 36.
    # Returns 1 when `n = 0`.

    ulong n_preinvert_limb_prenorm(ulong n)
    # Computes an approximate inverse ``invxl`` of the limb ``xl``,
    # with an implicit leading~`1`. More formally it computes::
    # invxl = (B^2 - B*x - 1)/x = (B^2 - 1)/x - B
    # Note that `x` must be normalised, i.e. with msb set. This inverse
    # makes use of the following theorem of Torbjorn Granlund and Peter
    # Montgomery~[Lemma~8.1][GraMon1994]_:
    # Let `d` be normalised, `d < B`, i.e. it fits in a word, and suppose
    # that `m d < B^2 \leq (m+1) d`. Let `0 \leq n \leq B d - 1`.  Write
    # `n = n_2 B + n_1 B/2 + n_0` with `n_1 = 0` or `1` and `n_0 < B/2`.
    # Suppose `q_1 B + q_0 = n_2 B + (n_2 + n_1) (m - B) + n_1 (d-B/2) + n_0`
    # and `0 \leq q_0 < B`. Then `0 \leq q_1 < B` and `0 \leq n - q_1 d < 2 d`.
    # In the theorem, `m` is the inverse of `d`. If we let
    # ``m = invxl + B`` and `d = x` we have `m d = B^2 - 1 < B^2` and
    # `(m+1) x = B^2 + d - 1 \geq B^2`.
    # The theorem is often applied as follows: note that `n_0` and `n_1 (d-B/2)`
    # are both less than `B/2`. Also note that `n_1 (m-B) < B`. Thus the sum of
    # all these terms contributes at most `1` to `q_1`. We are left with
    # `n_2 B + n_2 (m-B)`. But note that `(m-B)` is precisely our precomputed
    # inverse ``invxl``. If we write `q_1 B + q_0 = n_2 B + n_2 (m-B)`,
    # then from the theorem, we have `0 \leq n - q_1 d < 3 d`, i.e. the
    # quotient is out by at most `2` and is always either correct or too small.

    ulong n_preinvert_limb(ulong n)
    # Returns a precomputed inverse of `n`, as defined in [GraMol2010]_.
    # This precomputed inverse can be used with all of the functions that
    # take a precomputed inverse whose names are suffixed by ``_preinv``.
    # We require `n > 0`.

    double n_precompute_inverse(ulong n)
    # Returns a precomputed inverse of `n` with double precision value `1/n`.
    # This precomputed inverse can be used with all of the functions that
    # take a precomputed inverse whose names are suffixed by ``_precomp``.
    # We require `n > 0`.

    ulong n_mod_precomp(ulong a, ulong n, double ninv)
    # Returns `a \bmod{n}` given a precomputed inverse of `n` computed
    # by :func:`n_precompute_inverse`. We require ``n < 2^FLINT_D_BITS``
    # and ``a < 2^(FLINT_BITS-1)`` and `0 \leq a < n^2`.
    # We assume the processor is in the standard round to nearest
    # mode. Thus ``ninv`` is correct to `53` binary bits, the least
    # significant bit of which we shall call a place, and can be at most
    # half a place out. When `a` is multiplied by `ninv`, the binary
    # representation of `a` is exact and the mantissa is less than `2`, thus we
    # see that ``a * ninv`` can be at most one out in the mantissa. We now
    # truncate ``a * ninv`` to the nearest integer, which is always a round
    # down. Either we already have an integer, or we need to make a change down
    # of at least `1` in the last place. In the latter case we either get
    # precisely the exact quotient or below it as when we rounded the
    # product to the nearest place we changed by at most half a place.
    # In the case that truncating to an integer takes us below the
    # exact quotient, we have rounded down by less than `1` plus half a
    # place. But as the product is less than `n` and `n` is less than `2^{53}`,
    # half a place is less than `1`, thus we are out by less than `2` from
    # the exact quotient, i.e. the quotient we have computed is the
    # quotient we are after or one too small. That leaves only the case
    # where we had to round up to the nearest place which happened to
    # be an integer, so that truncating to an integer didn't change
    # anything. But this implies that the exact quotient `a/n` is less
    # than `2^{-54}` from an integer. We deal with this rare case by
    # subtracting 1 from the quotient. Then the quotient we have computed is
    # either exactly what we are after, or one too small.

    ulong n_mod2_precomp(ulong a, ulong n, double ninv)
    # Returns `a \bmod{n}` given a precomputed inverse of `n` computed by
    # :func:`n_precompute_inverse`. There are no restrictions on `a` or
    # on `n`.
    # As for :func:`n_mod_precomp` for `n < 2^{53}` and `a < n^2` the
    # computed quotient is either what we are after or one too large or small.
    # We deal with these cases. Otherwise we can be sure that the
    # top `52` bits of the quotient are computed correctly. We take
    # the remainder and adjust the quotient by multiplying the
    # remainder by ``ninv`` to compute another approximate quotient as
    # per :func:`mod_precomp`. Now the remainder may be either
    # negative or positive, so the quotient we compute may be one
    # out in either direction.

    ulong n_divrem2_preinv(ulong * q, ulong a, ulong n, ulong ninv)
    # Returns `a \bmod{n}` and sets `q` to the quotient of `a` by `n`, given a
    # precomputed inverse of `n` computed by :func:`n_preinvert_limb()`. There are
    # no restrictions on `a` and the only restriction on `n` is that it be
    # nonzero.
    # This uses the algorithm of Granlund and Möller [GraMol2010]_. First
    # `n` is normalised and `a` is shifted into two limbs to compensate. Then
    # their algorithm is applied verbatim and the remainder shifted back.

    ulong n_div2_preinv(ulong a, ulong n, ulong ninv)
    # Returns the Euclidean quotient of `a` by `n` given a precomputed inverse of
    # `n` computed by :func:`n_preinvert_limb`. There are no restrictions on `a`
    # and the only restriction on `n` is that it be nonzero.
    # This uses the algorithm of Granlund and Möller [GraMol2010]_. First
    # `n` is normalised and `a` is shifted into two limbs to compensate. Then
    # their algorithm is applied verbatim.

    ulong n_mod2_preinv(ulong a, ulong n, ulong ninv)
    # Returns `a \bmod{n}` given a precomputed inverse of `n` computed by
    # :func:`n_preinvert_limb()`. There are no restrictions on `a` and the only
    # restriction on `n` is that it be nonzero.
    # This uses the algorithm of Granlund and Möller [GraMol2010]_. First
    # `n` is normalised and `a` is shifted into two limbs to compensate. Then
    # their algorithm is applied verbatim and the result shifted back.

    ulong n_divrem2_precomp(ulong * q, ulong a, ulong n, double npre)
    # Returns `a \bmod{n}` given a precomputed inverse of `n` computed by
    # :func:`n_precompute_inverse` and sets `q` to the quotient. There
    # are no restrictions on `a` or on `n`.
    # This is as for :func:`n_mod2_precomp` with some additional care taken
    # to retain the quotient information. There are also special
    # cases to deal with the case where `a` is already reduced modulo
    # `n` and where `n` is `64` bits and `a` is not reduced modulo `n`.

    ulong n_ll_mod_preinv(ulong a_hi, ulong a_lo, ulong n, ulong ninv)
    # Returns `a \bmod{n}` given a precomputed inverse of `n` computed by
    # :func:`n_preinvert_limb`. There are no restrictions on `a`, which
    # will be two limbs ``(a_hi, a_lo)``, or on `n`.
    # The old version of this function merely reduced the top limb
    # ``a_hi`` modulo `n` so that :func:`udiv_qrnnd_preinv()` could
    # be used.
    # The new version reduces the top limb modulo `n` as per
    # :func:`n_mod2_preinv` and then the algorithm of Granlund and
    # Möller [GraMol2010]_ is used again to reduce modulo `n`.

    ulong n_lll_mod_preinv(ulong a_hi, ulong a_mi, ulong a_lo, ulong n, ulong ninv)
    # Returns `a \bmod{n}`, where `a` has three limbs ``(a_hi, a_mi, a_lo)``,
    # given a precomputed inverse of `n` computed by :func:`n_preinvert_limb`.
    # It is assumed that ``a_hi`` is reduced modulo `n`. There are no
    # restrictions on `n`.
    # This function uses the algorithm of Granlund and
    # Möller [GraMol2010]_ to first reduce the top two limbs
    # modulo `n`, then does the same on the bottom two limbs.

    ulong n_mulmod_precomp(ulong a, ulong b, ulong n, double ninv)
    # Returns `a b \bmod{n}` given a precomputed inverse of `n`
    # computed by :func:`n_precompute_inverse`. We require
    # ``n < 2^FLINT_D_BITS`` and `0 \leq a, b < n`.
    # We assume the processor is in the standard round to nearest
    # mode. Thus ``ninv`` is correct to `53` binary bits, the least
    # significant bit of which we shall call a place, and can be at most half
    # a place out. The product of `a` and `b` is computed with error at most
    # half a place. When ``a * b`` is multiplied by `ninv` we find that the
    # exact quotient and computed quotient differ by less than two places. As
    # the quotient is less than `n` this means that the exact quotient is at
    # most `1` away from the computed quotient. We truncate this quotient to
    # an integer which reduces the value by less than `1`. We end up with a
    # value which can be no more than two above the quotient we are after and
    # no less than two below. However an argument similar to that for
    # :func:`n_mod_precomp` shows that the truncated computed quotient cannot
    # be two smaller than the truncated exact quotient. In other words the
    # computed integer quotient is at most two above and one below the quotient
    # we are after.

    ulong n_mulmod2_preinv(ulong a, ulong b, ulong n, ulong ninv)
    # Returns `a b \bmod{n}` given a precomputed inverse of `n` computed by
    # :func:`n_preinvert_limb`. There are no restrictions on `a`, `b` or
    # on `n`. This is implemented by multiplying using :func:`umul_ppmm` and
    # then reducing using :func:`n_ll_mod_preinv`.

    ulong n_mulmod2(ulong a, ulong b, ulong n)
    # Returns `a b \bmod{n}`. There are no restrictions on `a`, `b` or
    # on `n`. This is implemented by multiplying using :func:`umul_ppmm` and
    # then reducing using :func:`n_ll_mod_preinv` after computing a precomputed
    # inverse.

    ulong n_mulmod_preinv(ulong a, ulong b, ulong n, ulong ninv, ulong norm)
    # Returns `a b \pmod{n}` given a precomputed inverse of `n` computed by
    # :func:`n_preinvert_limb`, assuming `a` and `b` are reduced modulo `n`
    # and `n` is normalised, i.e. with most significant bit set. There are
    # no other restrictions on `a`, `b` or `n`.
    # The value ``norm`` is provided for convenience. As `n` is required
    # to be normalised, it may be that `a` and `b` have been shifted to the
    # left by ``norm`` bits before calling the function. Their product
    # then has an extra factor of `2^\text{norm}`. Specifying a nonzero
    # ``norm`` will shift the product right by this many bits before
    # reducing it.
    # The algorithm used is that of Granlund and Möller [GraMol2010]_.

    ulong n_gcd(ulong x, ulong y)
    # Returns the greatest common divisor `g` of `x` and `y`. No assumptions
    # are made about the values `x` and `y`.
    # This function wraps GMP's ``mpn_gcd_1``.

    ulong n_gcdinv(ulong * a, ulong x, ulong y)
    # Returns the greatest common divisor `g` of `x` and `y` and computes
    # `a` such that `0 \leq a < y` and `a x = \gcd(x, y) \bmod{y}`, when
    # this is defined. We require `x < y`.
    # When `y = 1` the greatest common divisor is set to `1` and `a` is
    # set to `0`.
    # This is merely an adaption of the extended Euclidean algorithm
    # computing just one cofactor and reducing it modulo `y`.

    ulong n_xgcd(ulong * a, ulong * b, ulong x, ulong y)
    # Returns the greatest common divisor `g` of `x` and `y` and unsigned
    # values `a` and `b` such that `a x - b y = g`. We require `x \geq y`.
    # We claim that computing the extended greatest common divisor via the
    # Euclidean algorithm always results in cofactor `\lvert a \rvert < x/2`,
    # `\lvert b\rvert < x/2`, with perhaps some small degenerate exceptions.
    # We proceed by induction.
    # Suppose we are at some step of the algorithm, with `x_n = q y_n + r`
    # with `r \geq 1`, and suppose `1 = s y_n - t r` with
    # `s < y_n / 2`, `t < y_n / 2` by hypothesis.
    # Write `1 = s y_n - t (x_n - q y_n) = (s + t q) y_n - t x_n`.
    # It suffices to show that `(s + t q) < x_n / 2` as `t < y_n / 2 < x_n / 2`,
    # which will complete the induction step.
    # But at the previous step in the backsubstitution we would have had
    # `1 = s r - c d` with `s < r/2` and `c < r/2`.
    # Then `s + t q < r/2 + y_n / 2 q = (r + q y_n)/2 = x_n / 2`.
    # See the documentation of :func:`n_gcd` for a description of the
    # branching in the algorithm, which is faster than using division.

    int n_jacobi(mp_limb_signed_t x, ulong y)
    # Computes the Jacobi symbol `\left(\frac{x}{y}\right)` for any `x` and odd `y`.

    int n_jacobi_unsigned(ulong x, ulong y)
    # Computes the Jacobi symbol, allowing `x` to go up to a full limb.

    ulong n_addmod(ulong a, ulong b, ulong n)
    # Returns `(a + b) \bmod{n}`.

    ulong n_submod(ulong a, ulong b, ulong n)
    # Returns `(a - b) \bmod{n}`.

    ulong n_invmod(ulong x, ulong y)
    # Returns the inverse of `x` modulo `y`, if it exists. Otherwise an exception
    # is thrown.
    # This is merely an adaption of the extended Euclidean algorithm
    # with appropriate normalisation.

    ulong n_powmod_precomp(ulong a, mp_limb_signed_t exp, ulong n, double npre)
    # Returns ``a^exp`` modulo `n` given a precomputed inverse of `n`
    # computed by :func:`n_precompute_inverse`. We require `n < 2^{53}`
    # and `0 \leq a < n`. There are no restrictions on ``exp``, i.e.
    # it can be negative.
    # This is implemented as a standard binary powering algorithm using
    # repeated squaring and reducing modulo `n` at each step.

    ulong n_powmod_ui_precomp(ulong a, ulong exp, ulong n, double npre)
    # Returns ``a^exp`` modulo `n` given a precomputed inverse of `n`
    # computed by :func:`n_precompute_inverse`. We require `n < 2^{53}`
    # and `0 \leq a < n`. The exponent ``exp`` is unsigned and so
    # can be larger than allowed by :func:`n_powmod_precomp`.
    # This is implemented as a standard binary powering algorithm using
    # repeated squaring and reducing modulo `n` at each step.

    ulong n_powmod(ulong a, mp_limb_signed_t exp, ulong n)
    # Returns ``a^exp`` modulo `n`. We require ``n < 2^FLINT_D_BITS``
    # and `0 \leq a < n`. There are no restrictions on ``exp``, i.e.
    # it can be negative.
    # This is implemented by precomputing an inverse and calling the
    # ``precomp`` version of this function.

    ulong n_powmod2_preinv(ulong a, mp_limb_signed_t exp, ulong n, ulong ninv)
    # Returns ``(a^exp) % n`` given a precomputed inverse of `n` computed
    # by :func:`n_preinvert_limb`. We require `0 \leq a < n`, but there are no
    # restrictions on `n` or on ``exp``, i.e. it can be negative.
    # This is implemented as a standard binary powering algorithm using
    # repeated squaring and reducing modulo `n` at each step.
    # If ``exp`` is negative but `a` is not invertible modulo `n`, an
    # exception is raised.

    ulong n_powmod2(ulong a, mp_limb_signed_t exp, ulong n)
    # Returns ``(a^exp) % n``. We require `0 \leq a < n`, but there are
    # no restrictions on `n` or on ``exp``, i.e. it can be negative.
    # This is implemented by precomputing an inverse limb and calling the
    # ``preinv`` version of this function.
    # If ``exp`` is negative but `a` is not invertible modulo `n`, an
    # exception is raised.

    ulong n_powmod2_ui_preinv(ulong a, ulong exp, ulong n, ulong ninv)
    # Returns ``(a^exp) % n`` given a precomputed inverse of `n` computed
    # by :func:`n_preinvert_limb`. We require `0 \leq a < n`, but there are no
    # restrictions on `n`. The exponent ``exp`` is unsigned and so can be
    # larger than allowed by :func:`n_powmod2_preinv`.
    # This is implemented as a standard binary powering algorithm using
    # repeated squaring and reducing modulo `n` at each step.

    ulong n_powmod2_fmpz_preinv(ulong a, const fmpz_t exp, ulong n, ulong ninv)
    # Returns ``(a^exp) % n`` given a precomputed inverse of `n` computed
    # by :func:`n_preinvert_limb`. We require `0 \leq a < n`, but there are no
    # restrictions on `n`. The exponent ``exp`` must not be negative.
    # This is implemented as a standard binary powering algorithm using
    # repeated squaring and reducing modulo `n` at each step.

    ulong n_sqrtmod(ulong a, ulong p)
    # If `p` is prime, compute a square root of `a` modulo `p` if `a` is a
    # quadratic residue modulo `p`, otherwise return `0`.
    # If `p` is not prime the result is with high probability `0`, indicating
    # that `p` is not prime, or `a` is not a square modulo `p`. Otherwise the
    # result is meaningless.
    # Assumes that `a` is reduced modulo `p`.

    slong n_sqrtmod_2pow(ulong ** sqrt, ulong a, slong exp)
    # Computes all the square roots of ``a`` modulo ``2^exp``. The roots
    # are stored in an array which is created and whose address is stored in
    # the location pointed to by ``sqrt``. The array of roots is allocated
    # by the function but must be cleaned up by the user by calling
    # ``flint_free``. The number of roots is returned by the function. If
    # ``a`` is not a quadratic residue modulo ``2^exp`` then 0 is
    # returned by the function and the location ``sqrt`` points to is set to
    # NULL.

    slong n_sqrtmod_primepow(ulong ** sqrt, ulong a, ulong p, slong exp)
    # Computes all the square roots of ``a`` modulo ``p^exp``. The roots
    # are stored in an array which is created and whose address is stored in
    # the location pointed to by ``sqrt``. The array of roots is allocated
    # by the function but must be cleaned up by the user by calling
    # ``flint_free``. The number of roots is returned by the function. If
    # ``a`` is not a quadratic residue modulo ``p^exp`` then 0 is
    # returned by the function and the location ``sqrt`` points to is set to
    # NULL.

    slong n_sqrtmodn(ulong ** sqrt, ulong a, n_factor_t * fac)
    # Computes all the square roots of ``a`` modulo ``m`` given the
    # factorisation of ``m`` in ``fac``. The roots are stored in an array
    # which is created and whose address is stored in the location pointed to by
    # ``sqrt``. The array of roots is allocated by the function but must be
    # cleaned up by the user by calling :func:`flint_free`. The number of roots
    # is returned by the function. If ``a`` is not a quadratic residue modulo
    # ``m`` then 0 is returned by the function and the location ``sqrt``
    # points to is set to NULL.

    mp_limb_t n_mulmod_shoup(mp_limb_t w, mp_limb_t t, mp_limb_t w_precomp, mp_limb_t p)
    # Returns `w t \bmod{p}` given a precomputed scaled approximation of `w / p`
    # computed by :func:`n_mulmod_precomp_shoup`. The value of `p` should be
    # less than `2^{\mathtt{FLINT\_BITS} - 1}`. `w` and `t` should be less than `p`.
    # Works faster than :func:`n_mulmod2_preinv` if `w` fixed and `t` from array
    # (for example, scalar multiplication of vector).

    mp_limb_t n_mulmod_precomp_shoup(mp_limb_t w, mp_limb_t p)
    # Returns `w'`, scaled approximation of `w / p`. `w'`  is equal to the integer
    # part of `w \cdot 2^{\mathtt{FLINT\_BITS}} / p`.

    int n_divides(mp_limb_t * q, mp_limb_t n, mp_limb_t p)

    void n_primes_init(n_primes_t it)
    # Initialises the prime number iterator ``iter`` for use.

    void n_primes_clear(n_primes_t it)
    # Clears memory allocated by the prime number iterator ``iter``.

    ulong n_primes_next(n_primes_t it)
    # Returns the next prime number and advances the state of ``iter``.
    # The first call returns 2.
    # Small primes are looked up from ``flint_small_primes``.
    # When this table is exhausted, primes are generated in blocks
    # by calling :func:`n_primes_sieve_range`.

    void n_primes_jump_after(n_primes_t it, ulong n)
    # Changes the state of ``iter`` to start generating primes
    # after `n` (excluding `n` itself).

    void n_primes_extend_small(n_primes_t it, ulong bound)
    # Extends the table of small primes in ``iter`` to contain
    # at least two primes larger than or equal to ``bound``.

    void n_primes_sieve_range(n_primes_t it, ulong a, ulong b)
    # Sets the block endpoints of ``iter`` to the smallest and
    # largest odd numbers between `a` and `b` inclusive, and
    # sieves to mark all odd primes in this range.
    # The iterator state is changed to point to the first
    # number in the sieved range.

    void n_compute_primes(ulong num_primes)
    # Precomputes at least ``num_primes`` primes and their ``double``
    # precomputed inverses and stores them in an internal cache.
    # Assuming that FLINT has been built with support for thread-local storage,
    # each thread has its own cache.

    const ulong * n_primes_arr_readonly(ulong num_primes)
    # Returns a pointer to a read-only array of the first ``num_primes``
    # prime numbers. The computed primes are cached for repeated calls.
    # The pointer is valid until the user calls :func:`n_cleanup_primes`
    # in the same thread.

    const double * n_prime_inverses_arr_readonly(ulong n)
    # Returns a pointer to a read-only array of inverses of the first
    # ``num_primes`` prime numbers. The computed primes are cached for
    # repeated calls. The pointer is valid until the user calls
    # :func:`n_cleanup_primes` in the same thread.

    void n_cleanup_primes()
    # Frees the internal cache of prime numbers used by the current thread.
    # This will invalidate any pointers returned by
    # :func:`n_primes_arr_readonly` or :func:`n_prime_inverses_arr_readonly`.

    ulong n_nextprime(ulong n, int proved)
    # Returns the next prime after `n`. Assumes the result will fit in an
    # ``ulong``. If proved is `0`, i.e. false, the prime is not
    # proven prime, otherwise it is.

    ulong n_prime_pi(ulong n)
    # Returns the value of the prime counting function `\pi(n)`, i.e. the
    # number of primes less than or equal to `n`. The invariant
    # ``n_prime_pi(n_nth_prime(n)) == n``.
    # Currently, this function simply extends the table of cached primes up to
    # an upper limit and then performs a binary search.

    void n_prime_pi_bounds(ulong *lo, ulong *hi, ulong n)
    # Calculates lower and upper bounds for the value of the prime counting
    # function ``lo <= pi(n) <= hi``. If ``lo`` and ``hi`` point to
    # the same location, the high value will be stored.
    # This does a table lookup for small values, then switches over to some
    # proven bounds.
    # The upper approximation is `1.25506 n / \ln n`, and the
    # lower is `n / \ln n`.  These bounds are due to Rosser and
    # Schoenfeld [RosSch1962]_ and valid for `n \geq 17`.
    # We use the number of bits in `n` (or one less) to form an
    # approximation to `\ln n`, taking care to use a value too
    # small or too large to maintain the inequality.

    ulong n_nth_prime(ulong n)
    # Returns the `n`\th prime number `p_n`, using the mathematical indexing
    # convention `p_1 = 2, p_2 = 3, \dotsc`.
    # This function simply ensures that the table of cached primes is large
    # enough and then looks up the entry.

    void n_nth_prime_bounds(ulong *lo, ulong *hi, ulong n)
    # Calculates lower and upper bounds for the  `n`\th prime number `p_n` ,
    # ``lo <= p_n <= hi``. If ``lo`` and ``hi`` point to the same
    # location, the high value will be stored. Note that this function will
    # overflow for sufficiently large `n`.
    # We use the following estimates, valid for `n > 5` :
    # .. math ::
    # p_n  & >  n (\ln n + \ln \ln n - 1) \\
    # p_n  & <  n (\ln n + \ln \ln n) \\
    # p_n  & <  n (\ln n + \ln \ln n - 0.9427) \quad (n \geq 15985)
    # The first inequality was proved by Dusart [Dus1999]_, and the last
    # is due to Massias and Robin [MasRob1996]_.  For a further overview,
    # see http://primes.utm.edu/howmany.shtml .
    # We bound `\ln n` using the number of bits in `n` as in
    # ``n_prime_pi_bounds()``, and estimate `\ln \ln n` to the nearest
    # integer; this function is nearly constant.

    bint n_is_oddprime_small(ulong n)
    # Returns `1` if `n` is an odd prime smaller than
    # ``FLINT_ODDPRIME_SMALL_CUTOFF``. Expects `n`
    # to be odd and smaller than the cutoff.
    # This function merely uses a lookup table with one bit allocated for each
    # odd number up to the cutoff.

    bint n_is_oddprime_binary(ulong n)
    # This function performs a simple binary search through
    # the table of cached primes for `n`. If it exists in the array it returns
    # `1`, otherwise `0`. For the algorithm to operate correctly
    # `n` should be odd and at least `17`.
    # Lower and upper bounds are computed with :func:`n_prime_pi_bounds`.
    # Once we have bounds on where to look in the table, we
    # refine our search with a simple binary algorithm, taking
    # the top or bottom of the current interval as necessary.

    bint n_is_prime_pocklington(ulong n, ulong iterations)
    # Tests if `n` is a prime using the Pocklington--Lehmer primality
    # test. If `1` is returned `n` has been proved prime. If `0` is returned
    # `n` is composite. However `-1` may be returned if nothing was proved
    # either way due to the number of iterations being too small.
    # The most time consuming part of the algorithm is factoring
    # `n - 1`. For this reason :func:`n_factor_partial` is used,
    # which uses a combination of trial factoring and Hart's one
    # line factor algorithm [Har2012]_ to try to quickly factor `n - 1`.
    # Additionally if the cofactor is less than the square root of
    # `n - 1` the algorithm can still proceed.
    # One can also specify a number of iterations if less time
    # should be taken. Simply set this to ``WORD(0)`` if this is irrelevant.
    # In most cases a greater number of iterations will not
    # significantly affect timings as most of the time is spent
    # factoring.
    # See
    # https://mathworld.wolfram.com/PocklingtonsTheorem.html
    # for a description of the algorithm.

    bint n_is_prime_pseudosquare(ulong n)
    # Tests if `n` is a prime according to Theorem 2.7 [LukPatWil1996]_.
    # We first factor `N` using trial division up to some limit `B`.
    # In fact, the number of primes used in the trial factoring is at
    # most ``FLINT_PSEUDOSQUARES_CUTOFF``.
    # Next we compute `N/B` and find the next pseudosquare `L_p` above
    # this value, using a static table as per
    # https://oeis.org/A002189/b002189.txt .
    # As noted in the text, if `p` is prime then Step 3 will pass. This
    # test rejects many composites, and so by this time we suspect
    # that `p` is prime. If `N` is `3` or `7` modulo `8`, we are done,
    # and `N` is prime.
    # We now run a probable prime test, for which no known
    # counterexamples are known, to reject any composites. We then
    # proceed to prove `N` prime by executing Step 4. In the case that
    # `N` is `1` modulo `8`, if Step 4 fails, we extend the number of primes
    # `p_i` at Step 3 and hope to find one which passes Step 4. We take
    # the test one past the largest `p` for which we have pseudosquares
    # `L_p` tabulated, as this already corresponds to the next `L_p` which
    # is bigger than `2^{64}` and hence larger than any prime we might be
    # testing.
    # As explained in the text, Condition 4 cannot fail if `N` is prime.
    # The possibility exists that the probable prime test declares a
    # composite prime. However in that case an error is printed, as
    # that would be of independent interest.

    bint n_is_prime(ulong n)
    # Tests if `n` is a prime. This first sieves for small prime factors,
    # then simply calls :func:`n_is_probabprime`. This has been checked
    # against the tables of Feitsma and Galway
    # http://www.cecm.sfu.ca/Pseudoprimes/index-2-to-64.html and thus
    # constitutes a check for primality (rather than just pseudoprimality)
    # up to `2^{64}`.
    # In future, this test may produce and check a certificate of
    # primality. This is likely to be significantly slower for prime
    # inputs.

    bint n_is_strong_probabprime_precomp(ulong n, double npre, ulong a, ulong d)
    # Tests if `n` is a strong probable prime to the base `a`. We
    # require that `d` is set to the largest odd factor of `n - 1` and
    # ``npre`` is a precomputed inverse of `n` computed with
    # :func:`n_precompute_inverse`.  We also require that `n < 2^{53}`,
    # `a` to be reduced modulo `n` and not `0` and `n` to be odd.
    # If we write `n - 1 = 2^s d` where `d` is odd then `n` is a strong
    # probable prime to the base `a`, i.e. an `a`-SPRP, if either
    # `a^d = 1 \pmod n` or `(a^d)^{2^r} = -1 \pmod n` for some `r` less
    # than `s`.
    # A description of strong probable primes is given here:
    # https://mathworld.wolfram.com/StrongPseudoprime.html

    bint n_is_strong_probabprime2_preinv(ulong n, ulong ninv, ulong a, ulong d)
    # Tests if `n` is a strong probable prime to the base `a`. We require
    # that `d` is set to the largest odd factor of `n - 1` and ``npre``
    # is a precomputed inverse of `n` computed with :func:`n_preinvert_limb`.
    # We require a to be reduced modulo `n` and not `0` and `n` to be odd.
    # If we write `n - 1 = 2^s d` where `d` is odd then `n` is a strong
    # probable prime to the base `a` (an `a`-SPRP) if either `a^d = 1 \pmod n`
    # or `(a^d)^{2^r} = -1 \pmod n` for some `r` less than `s`.
    # A description of strong probable primes is given here:
    # https://mathworld.wolfram.com/StrongPseudoprime.html

    bint n_is_probabprime_fermat(ulong n, ulong i)
    # Returns `1` if `n` is a base `i` Fermat probable prime. Requires
    # `1 < i < n` and that `i` does not divide `n`.
    # By Fermat's Little Theorem if `i^{n-1}` is not congruent to `1`
    # then `n` is not prime.

    bint n_is_probabprime_fibonacci(ulong n)
    # Let `F_j` be the `j`\th element of the Fibonacci sequence
    # `0, 1, 1, 2, 3, 5, \dotsc`, starting at `j = 0`. Then if `n` is prime
    # we have `F_{n - (n/5)} = 0 \pmod n`, where `(n/5)` is the Jacobi
    # symbol.
    # For further details, see  pp. 142 [CraPom2005]_.
    # We require that `n` is not divisible by `2` or `5`.

    bint n_is_probabprime_BPSW(ulong n)
    # Implements a Baillie--Pomerance--Selfridge--Wagstaff probable primality
    # test. This is a variant of the usual BPSW test (which only uses strong
    # base-2 probable prime and Lucas-Selfridge tests, see Baillie and
    # Wagstaff [BaiWag1980]_).
    # This implementation makes use of a weakening of the usual Baillie-PSW
    # test given in  [Chen2003]_, namely replacing the Lucas test with a
    # Fibonacci test when `n \equiv 2, 3 \pmod{5}` (see also the comment on
    # page 143 of [CraPom2005]_), regarding Fibonacci pseudoprimes.
    # There are no known counterexamples to this being a primality test.
    # Up to `2^{64}` the test we use has been checked against tables of
    # pseudoprimes. Thus it is a primality test up to this limit.

    bint n_is_probabprime_lucas(ulong n)
    # For details on Lucas pseudoprimes, see [pp. 143] [CraPom2005]_.
    # We implement a variant of the Lucas pseudoprime test similar to that
    # described by Baillie and Wagstaff [BaiWag1980]_.

    bint n_is_probabprime(ulong n)
    # Tests if `n` is a probable prime. Up to ``FLINT_ODDPRIME_SMALL_CUTOFF``
    # this algorithm uses :func:`n_is_oddprime_small` which uses a lookup table.
    # Next it calls :func:`n_compute_primes` with the maximum table size and
    # uses this table to perform a binary search for `n` up to the table limit.
    # Then up to `1050535501` it uses a number of strong probable prime tests,
    # :func:`n_is_strong_probabprime_preinv`, etc., for various bases. The
    # output of the algorithm is guaranteed to be correct up to this bound due
    # to exhaustive tables, described at
    # http://uucode.com/obf/dalbec/alg.html .
    # Beyond that point the BPSW probabilistic primality test is used, by
    # calling the function :func:`n_is_probabprime_BPSW`. There are no known
    # counterexamples, and it has been checked against the tables of Feitsma
    # and Galway and up to the accuracy of those tables, this is an exhaustive
    # check up to `2^{64}`, i.e. there are no counterexamples.

    ulong n_CRT(ulong r1, ulong m1, ulong r2, ulong m2)
    # Use the Chinese Remainder Theorem to return the unique value
    # `0 \le x < M` congruent to `r_1` modulo `m_1` and `r_2` modulo `m_2`,
    # where `M = m_1 \times m_2` is assumed to fit a ulong.
    # It is assumed that `m_1` and `m_2` are positive integers greater
    # than `1` and coprime. It is assumed that `0 \le r_1 < m_1` and `0 \le r_2 < m_2`.

    ulong n_sqrt(ulong a)
    # Computes the integer truncation of the square root of `a`.
    # The implementation uses a call to the IEEE floating point sqrt function.
    # The integer itself is represented by the nearest double and its square
    # root is computed to the nearest place. If `a` is one below a square, the
    # rounding may be up, whereas if it is one above a square, the rounding
    # will be down. Thus the square root may be one too large in some
    # instances which we then adjust by checking if we have the right value.
    # We also have to be careful when the square of this too large
    # value causes an overflow. The same assumptions hold for a single
    # precision float provided the square root itself can be represented
    # in a single float, i.e. for `a < 281474976710656 = 2^{46}`.

    ulong n_sqrtrem(ulong * r, ulong a)
    # Computes the integer truncation of the square root of `a`.
    # The integer itself is represented by the nearest double and its square
    # root is computed to the nearest place. If `a` is one below a square, the
    # rounding may be up, whereas if it is one above a square, the rounding
    # will be down. Thus the square root may be one too large in some
    # instances which we then adjust by checking if we have the right value.
    # We also have to be careful when the square of this too
    # large value causes an overflow. The same assumptions hold for a
    # single precision float provided the square root itself can be
    # represented in a single float, i.e. for \
    # `a < 281474976710656 = 2^{46}`.
    # The remainder is computed by subtracting the square of the computed square
    # root from `a`.

    bint n_is_square(ulong x)
    # Returns `1` if `x` is a square, otherwise `0`.
    # This code first checks if `x` is a square modulo `64`,
    # `63 = 3 \times 3 \times 7` and `65 = 5 \times 13`, using lookup tables,
    # and if so it then takes a square root and checks that the square of this
    # equals the original value.

    bint n_is_perfect_power235(ulong n)
    # Returns `1` if `n` is a perfect square, cube or fifth power.
    # This function uses a series of modular tests to reject most
    # non 235-powers. Each modular test returns a value from 0 to 7
    # whose bits respectively indicate whether the value is a square,
    # cube or fifth power modulo the given modulus. When these are
    # logically ``AND``-ed together, this gives a powerful test which will
    # reject most non-235 powers.
    # If a bit remains set indicating it may be a square, a standard
    # square root test is performed. Similarly a cube root or fifth
    # root can be taken, if indicated, to determine whether the power
    # of that root is exactly equal to `n`.

    bint n_is_perfect_power(ulong * root, ulong n)
    # If `n = r^k`, return `k` and set ``root`` to `r`. Note that `0` and
    # `1` are considered squares. No guarantees are made about `r` or `k`
    # being the minimum possible value.

    ulong n_rootrem(ulong* remainder, ulong n, ulong root)
    # This function uses the Newton iteration method to calculate the nth root of
    # a number.
    # First approximation is calculated by an algorithm mentioned in this
    # article:  https://en.wikipedia.org/wiki/Fast_inverse_square_root .
    # Instead of the inverse square root, the nth root is calculated.
    # Returns the integer part of ``n ^ 1/root``. Remainder is set as
    # ``n - base^root``. In case `n < 1` or ``root < 1``, `0` is returned.

    ulong n_cbrt(ulong n)
    # This function returns the integer truncation of the cube root of `n`.
    # First approximation is calculated by an algorithm mentioned in this
    # article: https://en.wikipedia.org/wiki/Fast_inverse_square_root .
    # Instead of the inverse square root, the cube root is calculated.
    # This functions uses different algorithms to calculate the cube root,
    # depending upon the size of `n`. For numbers greater than `2^{46}`, it uses
    # :func:`n_cbrt_chebyshev_approx`. Otherwise, it makes use of the iteration,
    # `x \leftarrow x - (x\cdot x\cdot x - a)\cdot x/(2\cdot x\cdot x\cdot x + a)` for getting a good estimate,
    # as mentioned in the paper by W. Kahan [Kahan1991]_ .

    ulong n_cbrt_newton_iteration(ulong n)
    # This function returns the integer truncation of the cube root of `n`.
    # Makes use of Newton iterations to get a close value, and then adjusts the
    # estimate so as to get the correct value.

    ulong n_cbrt_binary_search(ulong n)
    # This function returns the integer truncation of the cube root of `n`.
    # Uses binary search to get the correct value.

    ulong n_cbrt_chebyshev_approx(ulong n)
    # This function returns the integer truncation of the cube root of `n`.
    # The number is first expressed in the form ``x * 2^exp``. This ensures
    # `x` is in the range [0.5, 1]. Cube root of x is calculated using
    # Chebyshev's approximation polynomial for the function `y = x^{1/3}`. The
    # values of the coefficient are calculated from the Python module mpmath,
    # https://mpmath.org, using the function chebyfit. x is multiplied
    # by ``2^exp`` and the cube root of 1, 2 or 4 (according to ``exp%3``).

    ulong n_cbrtrem(ulong* remainder, ulong n)
    # This function returns the integer truncation of the cube root of `n`.
    # Remainder is set as `n` minus the cube of the value returned.

    void n_factor_init(n_factor_t * factors)
    # Initializes factors.

    int n_remove(ulong * n, ulong p)
    # Removes the highest possible power of `p` from `n`, replacing
    # `n` with the quotient. The return value is the highest
    # power of `p` that divided `n`. Assumes `n` is not `0`.
    # For `p = 2` trailing zeroes are counted. For other primes
    # `p` is repeatedly squared and stored in a table of powers
    # with the current highest power of `p` removed at each step
    # until no higher power can be removed. The algorithm then
    # proceeds down the power tree again removing powers of `p`
    # until none remain.

    int n_remove2_precomp(ulong * n, ulong p, double ppre)
    # Removes the highest possible power of `p` from `n`, replacing
    # `n` with the quotient. The return value is the highest
    # power of `p` that divided `n`. Assumes `n` is not `0`. We require
    # ``ppre`` to be set to a precomputed inverse of `p` computed
    # with :func:`n_precompute_inverse`.
    # For `p = 2` trailing zeroes are counted. For other primes
    # `p` we make repeated use of :func:`n_divrem2_precomp` until division
    # by `p` is no longer possible.

    void n_factor_insert(n_factor_t * factors, ulong p, ulong exp)
    # Inserts the given prime power factor ``p^exp`` into
    # the ``n_factor_t`` ``factors``. See the documentation for
    # :func:`n_factor_trial` for a description of the ``n_factor_t`` type.
    # The algorithm performs a simple search to see if `p` already
    # exists as a prime factor in the structure. If so the exponent
    # there is increased by the supplied exponent. Otherwise a new
    # factor ``p^exp`` is added to the end of the structure.
    # There is no test code for this function other than its use by
    # the various factoring functions, which have test code.

    ulong n_factor_trial_range(n_factor_t * factors, ulong n, ulong start, ulong num_primes)
    # Trial factor `n` with the first ``num_primes`` primes, but
    # starting at the prime with index start (counting from zero).
    # One requires an initialised ``n_factor_t`` structure, but factors
    # will be added by default to an already used ``n_factor_t``. Use
    # the function :func:`n_factor_init` defined in ``ulong_extras`` if
    # initialisation has not already been completed on factors.
    # Once completed, ``num`` will contain the number of distinct
    # prime factors found. The field `p` is an array of ``ulong``\s
    # containing the distinct prime factors, ``exp`` an array
    # containing the corresponding exponents.
    # The return value is the unfactored cofactor after trial
    # factoring is done.
    # The function calls :func:`n_compute_primes` automatically. See
    # the documentation for that function regarding limits.
    # The algorithm stops when the current prime has a square
    # exceeding `n`, as no prime factor of `n` can exceed this
    # unless `n` is prime.
    # The precomputed inverses of all the primes computed by
    # :func:`n_compute_primes` are utilised with the :func:`n_remove2_precomp`
    # function.

    ulong n_factor_trial(n_factor_t * factors, ulong n, ulong num_primes)
    # This function calls :func:`n_factor_trial_range`, with the value of
    # `0` for ``start``. By default this adds factors to an already existing
    # ``n_factor_t`` or to a newly initialised one.

    ulong n_factor_power235(ulong *exp, ulong n)
    # Returns `0` if `n` is not a perfect square, cube or fifth power.
    # Otherwise it returns the root and sets ``exp`` to either `2`,
    # `3` or `5` appropriately.
    # This function uses a series of modular tests to reject most
    # non 235-powers. Each modular test returns a value from 0 to 7
    # whose bits respectively indicate whether the value is a square,
    # cube or fifth power modulo the given modulus. When these are
    # logically ``AND``-ed together, this gives a powerful test which will
    # reject most non-235 powers.
    # If a bit remains set indicating it may be a square, a standard
    # square root test is performed. Similarly a cube root or fifth
    # root can be taken, if indicated, to determine whether the power
    # of that root is exactly equal to `n`.

    ulong n_factor_one_line(ulong n, ulong iters)
    # This implements Bill Hart's one line factoring algorithm [Har2012]_.
    # It is a variant of Fermat's algorithm which cycles through a large number
    # of multipliers instead of incrementing the square root. It is faster than
    # SQUFOF for `n` less than about `2^{40}`.

    ulong n_factor_lehman(ulong n)
    # Lehman's factoring algorithm. Currently works up to `10^{16}`, but is
    # not particularly efficient and so is not used in the general factor
    # function. Always returns a factor of `n`.

    ulong n_factor_SQUFOF(ulong n, ulong iters)
    # Attempts to split `n` using the given number of iterations
    # of SQUFOF. Simply set ``iters`` to ``WORD(0)`` for maximum
    # persistence.
    # The version of SQUFOF implemented here is as described by Gower
    # and Wagstaff [GowWag2008]_.
    # We start by trying SQUFOF directly on `n`. If that fails we
    # multiply it by each of the primes in ``flint_primes_small`` in
    # turn. As this multiplication may result in a two limb value
    # we allow this in our implementation of SQUFOF. As SQUFOF
    # works with values about half the size of `n` it only needs
    # single limb arithmetic internally.
    # If SQUFOF fails to factor `n` we return `0`, however with
    # ``iters`` large enough this should never happen.

    void n_factor(n_factor_t * factors, ulong n, int proved)
    # Factors `n` with no restrictions on `n`. If the prime factors are
    # required to be checked with a primality test, one may set
    # ``proved`` to `1`, otherwise set it to `0`, and they will only be
    # probable primes. NB: at the present there is no difference because
    # the probable prime tests have been exhaustively tested up to `2^{64}`.
    # However, in future, this flag may produce and separately check
    # a primality certificate. This may be quite slow (and probably no
    # less reliable in practice).
    # For details on the ``n_factor_t`` structure, see
    # :func:`n_factor_trial`.
    # This function first tries trial factoring with a number of primes
    # specified by the constant ``FLINT_FACTOR_TRIAL_PRIMES``. If the
    # cofactor is `1` or prime the function returns with all the factors.
    # Otherwise, the cofactor is placed in the array ``factor_arr``. Whilst
    # there are factors remaining in there which have not been split, the
    # algorithm continues. At each step each factor is first checked to
    # determine if it is a perfect power. If so it is replaced by the power
    # that has been found. Next if the factor is small enough and composite,
    # in particular, less than ``FLINT_FACTOR_ONE_LINE_MAX`` then
    # :func:`n_factor_one_line` is called with
    # ``FLINT_FACTOR_ONE_LINE_ITERS`` to try and split the factor. If
    # that fails or the factor is too large for :func:`n_factor_one_line`
    # then :func:`n_factor_SQUFOF` is called, with
    # ``FLINT_FACTOR_SQUFOF_ITERS``. If that fails an error results and
    # the program aborts. However this should not happen in practice.

    ulong n_factor_trial_partial(n_factor_t * factors, ulong n, ulong * prod, ulong num_primes, ulong limit)
    # Attempts trial factoring of `n` with the first ``num_primes primes``,
    # but stops when the product of prime factors so far exceeds ``limit``.
    # One requires an initialised ``n_factor_t`` structure, but factors
    # will be added by default to an already used ``n_factor_t``. Use
    # the function :func:`n_factor_init` defined in ``ulong_extras`` if
    # initialisation has not already been completed on ``factors``.
    # Once completed, ``num`` will contain the number of distinct
    # prime factors found. The field `p` is an array of ``ulong``\s
    # containing the distinct prime factors, ``exp`` an array
    # containing the corresponding exponents.
    # The return value is the unfactored cofactor after trial
    # factoring is done. The value ``prod`` will be set to the product
    # of the factors found.
    # The function calls :func:`n_compute_primes` automatically. See
    # the documentation for that function regarding limits.
    # The algorithm stops when the current prime has a square
    # exceeding `n`, as no prime factor of `n` can exceed this
    # unless `n` is prime.
    # The precomputed inverses of all the primes computed by
    # :func:`n_compute_primes` are utilised with the :func:`n_remove2_precomp`
    # function.

    ulong n_factor_partial(n_factor_t * factors, ulong n, ulong limit, int proved)
    # Factors `n`, but stops when the product of prime factors so far
    # exceeds ``limit``.
    # One requires an initialised ``n_factor_t`` structure, but factors
    # will be added by default to an already used ``n_factor_t``. Use
    # the function ``n_factor_init()`` defined in ``ulong_extras`` if
    # initialisation has not already been completed on ``factors``.
    # On exit, ``num`` will contain the number of distinct prime factors
    # found. The field `p` is an array of ``ulong``\s containing the
    # distinct prime factors, ``exp`` an array containing the corresponding
    # exponents.
    # The return value is the unfactored cofactor after factoring is done.
    # The factors are proved prime if ``proved`` is `1`, otherwise
    # they are merely probably prime.

    ulong n_factor_pp1(ulong n, ulong B1, ulong c)
    # Factors `n` using Williams' `p + 1` factoring algorithm, with prime
    # limit set to `B1`. We require `c` to be set to a random value. Each
    # trial of the algorithm with a different value of `c` gives another
    # chance to factor `n`, with roughly exponentially decreasing chance
    # of finding a missing factor. If `p + 1` (or `p - 1`) is not smooth
    # for any factor `p` of `n`, the algorithm will never succeed. The
    # value `c` should be less than `n` and greater than `2`.
    # If the algorithm succeeds, it returns the factor, otherwise it
    # returns `0` or `1` (the trivial factors modulo `n`).

    ulong n_factor_pp1_wrapper(ulong n)
    # A simple wrapper around ``n_factor_pp1`` which works in the range
    # `31`-`64` bits. Below this point, trial factoring will always succeed.
    # This function mainly exists for ``n_factor`` and is tuned to minimise
    # the time for ``n_factor`` on numbers that reach the ``n_factor_pp1``
    # stage, i.e. after trial factoring and one line factoring.

    int n_factor_pollard_brent_single(mp_limb_t *factor, mp_limb_t n, mp_limb_t ninv, mp_limb_t ai, mp_limb_t xi, mp_limb_t normbits, mp_limb_t max_iters)
    # Pollard Rho algorithm (with Brent modification) for integer factorization.
    # Assumes that the `n` is not prime. `factor` is set as the factor if found.
    # It is not assured that the factor found will be prime. Does not compute the complete
    # factorization, just one factor. Returns 1 if factorization is successful
    # (non trivial factor is found), else returns 0. Assumes `n` is normalized
    # (shifted by normbits bits), and takes as input a precomputed inverse of `n` as
    # computed by :func:`n_preinvert_limb`. `ai` and `xi` should also be shifted
    # left by `normbits`.
    # `ai` is the constant of the polynomial used, `xi` is the initial value.
    # `max\_iters` is the number of iterations tried in process of finding the
    # cycle.
    # The algorithm used is a modification of the original Pollard Rho algorithm,
    # suggested by Richard Brent in the paper, available at
    # https://maths-people.anu.edu.au/~brent/pd/rpb051i.pdf

    int n_factor_pollard_brent(mp_limb_t *factor, flint_rand_t state, mp_limb_t n_in, mp_limb_t max_tries, mp_limb_t max_iters)
    # Pollard Rho algorithm, modified as suggested by Richard Brent. Makes a call to
    # :func:`n_factor_pollard_brent_single`. The input parameters ai and xi for
    # :func:`n_factor_pollard_brent_single` are selected at random.
    # If the algorithm fails to find a non trivial factor in one call, it tries again
    # (this time with a different set of random values). This process is repeated a
    # maximum of `max\_tries` times.
    # Assumes `n` is not prime. `factor` is set as the factor found, if factorization
    # is successful. In such a case, 1 is returned. Otherwise, 0 is returned. Factor
    # discovered is not necessarily prime.

    int n_moebius_mu(ulong n)
    # Computes the Moebius function `\mu(n)`, which is defined as `\mu(n) = 0`
    # if `n` has a prime factor of multiplicity greater than `1`, `\mu(n) = -1`
    # if `n` has an odd number of distinct prime factors, and `\mu(n) = 1` if
    # `n` has an even number of distinct prime factors. By convention,
    # `\mu(0) = 0`.
    # For even numbers, we use the identities `\mu(4n) = 0` and
    # `\mu(2n) = - \mu(n)`. Odd numbers up to a cutoff are then looked up from
    # a precomputed table storing `\mu(n) + 1` in groups of two bits.
    # For larger `n`, we first check if `n` is divisible by a small odd square
    # and otherwise call ``n_factor()`` and count the factors.

    void n_moebius_mu_vec(int * mu, ulong len)
    # Computes `\mu(n)` for ``n = 0, 1, ..., len - 1``. This
    # is done by sieving over each prime in the range, flipping the sign
    # of `\mu(n)` for every multiple of a prime `p` and setting `\mu(n) = 0`
    # for every multiple of `p^2`.

    bint n_is_squarefree(ulong n)
    # Returns `0` if `n` is divisible by some perfect square, and `1` otherwise.
    # This simply amounts to testing whether `\mu(n) \neq 0`. As special
    # cases, `1` is considered squarefree and `0` is not considered squarefree.

    ulong n_euler_phi(ulong n)
    # Computes the Euler totient function `\phi(n)`, counting the number of
    # positive integers less than or equal to `n` that are coprime to `n`.

    ulong n_factorial_fast_mod2_preinv(ulong n, ulong p, ulong pinv)
    # Returns `n! \bmod p` given a precomputed inverse of `p` as computed
    # by :func:`n_preinvert_limb`. `p` is not required to be a prime, but
    # no special optimisations are made for composite `p`.
    # Uses fast multipoint evaluation, running in about `O(n^{1/2})` time.

    ulong n_factorial_mod2_preinv(ulong n, ulong p, ulong pinv)
    # Returns `n! \bmod p` given a precomputed inverse of `p` as computed
    # by :func:`n_preinvert_limb`. `p` is not required to be a prime, but
    # no special optimisations are made for composite `p`.
    # Uses a lookup table for small `n`, otherwise computes the product
    # if `n` is not too large, and calls the fast algorithm for extremely
    # large `n`.

    ulong n_primitive_root_prime_prefactor(ulong p, n_factor_t * factors)
    # Returns a primitive root for the multiplicative subgroup of `\mathbb{Z}/p\mathbb{Z}`
    # where `p` is prime given the factorisation (``factors``) of `p - 1`.

    ulong n_primitive_root_prime(ulong p)
    # Returns a primitive root for the multiplicative subgroup of `\mathbb{Z}/p\mathbb{Z}`
    # where `p` is prime.

    ulong n_discrete_log_bsgs(ulong b, ulong a, ulong n)
    # Returns the discrete logarithm of `b` with  respect to `a` in the
    # multiplicative subgroup of `\mathbb{Z}/n\mathbb{Z}` when `\mathbb{Z}/n\mathbb{Z}`
    # is cyclic. That is,
    # it returns a number `x` such that `a^x = b \bmod n`.  The
    # multiplicative subgroup is only cyclic when `n` is `2`, `4`,
    # `p^k`, or `2p^k` where `p` is an odd prime and `k` is a positive
    # integer.

    void n_factor_ecm_double(mp_limb_t *x, mp_limb_t *z, mp_limb_t x0, mp_limb_t z0, mp_limb_t n, n_ecm_t n_ecm_inf)
    # Sets the point `(x : z)` to two times `(x_0 : z_0)` modulo `n` according
    # to the formula
    # `x = (x_0 + z_0)^2 \cdot (x_0 - z_0)^2 \mod n,`
    # `z = 4 x_0 z_0 \left((x_0 - z_0)^2 + 4a_{24}x_0z_0\right) \mod n.`
    # This group doubling is valid only for points expressed in
    # Montgomery projective coordinates.

    void n_factor_ecm_add(mp_limb_t *x, mp_limb_t *z, mp_limb_t x1, mp_limb_t z1, mp_limb_t x2, mp_limb_t z2, mp_limb_t x0, mp_limb_t z0, mp_limb_t n, n_ecm_t n_ecm_inf)
    # Sets the point `(x : z)` to the sum of `(x_1 : z_1)` and `(x_2 : z_2)`
    # modulo `n`, given the difference `(x_0 : z_0)` according to the formula
    # This group doubling is valid only for points expressed in
    # Montgomery projective coordinates.

    void n_factor_ecm_mul_montgomery_ladder(mp_limb_t *x, mp_limb_t *z, mp_limb_t x0, mp_limb_t z0, mp_limb_t k, mp_limb_t n, n_ecm_t n_ecm_inf)
    # Montgomery ladder algorithm for scalar multiplication of elliptic points.
    # Sets the point `(x : z)` to `k(x_0 : z_0)` modulo `n`.
    # Valid only for points expressed in Montgomery projective coordinates.

    int n_factor_ecm_select_curve(mp_limb_t *f, mp_limb_t sigma, mp_limb_t n, n_ecm_t n_ecm_inf)
    # Selects a random elliptic curve given a random integer ``sigma``,
    # according to Suyama's parameterization. If the factor is found while
    # selecting the curve, `1` is returned. In case the curve found is not
    # suitable, `0` is returned.
    # Also selects the initial point `x_0`, and the value of `(a + 2)/4`, where `a`
    # is a curve parameter. Sets `z_0` as `1` (shifted left by
    # ``n_ecm_inf->normbits``). All these are stored in the
    # ``n_ecm_t`` struct.
    # The curve selected is of Montgomery form, the points selected satisfy the
    # curve and are projective coordinates.

    int n_factor_ecm_stage_I(mp_limb_t *f, const mp_limb_t *prime_array, mp_limb_t num, mp_limb_t B1, mp_limb_t n, n_ecm_t n_ecm_inf)
    # Stage I implementation of the ECM algorithm.
    # ``f`` is set as the factor if found. ``num`` is number of prime numbers
    # `<=` the bound ``B1``. ``prime_array`` is an array of first ``B1``
    # primes. `n` is the number being factored.
    # If the factor is found, `1` is returned, otherwise `0`.

    int n_factor_ecm_stage_II(mp_limb_t *f, mp_limb_t B1, mp_limb_t B2, mp_limb_t P, mp_limb_t n, n_ecm_t n_ecm_inf)
    # Stage II implementation of the ECM algorithm.
    # ``f`` is set as the factor if found. ``B1``, ``B2`` are the two
    # bounds. ``P`` is the primorial (approximately equal to `\sqrt{B2}`).
    # `n` is the number being factored.
    # If the factor is found, `1` is returned, otherwise `0`.

    int n_factor_ecm(mp_limb_t *f, mp_limb_t curves, mp_limb_t B1, mp_limb_t B2, flint_rand_t state, mp_limb_t n)
    # Outer wrapper function for the ECM algorithm. It factors `n` which
    # must fit into a ``mp_limb_t``.
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
