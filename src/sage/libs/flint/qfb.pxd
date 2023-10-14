# distutils: libraries = flint
# distutils: depends = flint/qfb.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void qfb_init(qfb_t q)
    # Initialise a ``qfb_t`` `q` for use.

    void qfb_clear(qfb_t q)
    # Clear a ``qfb_t`` after use. This releases any memory allocated for
    # `q` back to flint.

    void qfb_array_clear(qfb ** forms, slong num)
    # Clean up an array of ``qfb`` structs allocated by a qfb function.
    # The parameter ``num`` must be set to the length of the array.

    qfb_hash_t * qfb_hash_init(slong depth)
    # Initialises a hash table of size `2^{depth}`.

    void qfb_hash_clear(qfb_hash_t * qhash, slong depth)
    # Frees all memory used by a hash table of size `2^{depth}`.

    void qfb_hash_insert(qfb_hash_t * qhash, qfb_t q, qfb_t q2, slong it, slong depth)
    # Insert the binary quadratic form ``q`` into the given hash table
    # of size `2^{depth}` in the field ``q`` of the hash structure.
    # Also store the second binary quadratic form ``q2`` (if not
    # ``NULL``) in the similarly named field and ``iter`` in the
    # similarly named field of the hash structure.

    slong qfb_hash_find(qfb_hash_t * qhash, qfb_t q, slong depth)
    # Search for the given binary quadratic form or its inverse in the
    # given hash table of size `2^{depth}`. If it is found, return
    # the index in the table (which is an array of ``qfb_hash_t``
    # structs), otherwise return ``-1``.

    void qfb_set(qfb_t f, qfb_t g)
    # Set the binary quadratic form `f` to be equal to `g`.

    bint qfb_equal(qfb_t f, qfb_t g)
    # Returns `1` if `f` and `g` are identical binary quadratic forms,
    # otherwise returns `0`.

    void qfb_print(qfb_t q)
    # Print a binary quadratic form `q` in the format `(a, b, c)` where
    # `a`, `b`, `c` are the entries of `q`.

    void qfb_discriminant(fmpz_t D, qfb_t f)
    # Set `D` to the discriminant of the binary quadratic form `f`, i.e. to
    # `b^2 - 4ac`, where `f = (a, b, c)`.

    void qfb_reduce(qfb_t r, qfb_t f, fmpz_t D)
    # Set `r` to a reduced form equivalent to the binary quadratic form `f`
    # of discriminant `D`.

    bint qfb_is_reduced(qfb_t r)
    # Returns `1` if `q` is a reduced binary quadratic form, otherwise
    # returns `0`. Note that this only tests for definite quadratic
    # forms, so a form `r = (a,b,c)` is reduced if and only if `|b| \le a \le
    # c` and if either inequality is an equality, then `b \ge 0`.

    slong qfb_reduced_forms(qfb ** forms, slong d)
    # Given a discriminant `d` (negative for negative definite forms), compute
    # all the reduced binary quadratic forms of that discriminant. The function
    # allocates space for these and returns it in the variable ``forms``
    # (the user is responsible for cleaning this up by a single call to
    # ``qfb_array_clear`` on ``forms``, after use.) The function returns
    # the number of forms generated (the form class number). The forms are
    # stored in an array of ``qfb`` structs, which contain fields
    # ``a, b, c`` corresponding to forms `(a, b, c)`.

    slong qfb_reduced_forms_large(qfb ** forms, slong d)
    # As for ``qfb_reduced_forms``. However, for small `|d|` it requires
    # fewer primes to be computed at a small cost in speed. It is called
    # automatically by ``qfb_reduced_forms`` for large `|d|` so that
    # ``flint_primes`` is not exhausted.

    void qfb_nucomp(qfb_t r, const qfb_t f, const qfb_t g, fmpz_t D, fmpz_t L)
    # Shanks' NUCOMP as described in [JvdP2002]_.
    # Computes the near reduced composition of forms `f` and `g` given
    # `L = \lfloor |D|^{1/4} \rfloor` where `D` is the common discriminant of
    # `f` and `g`. The result is returned in `r`.
    # We require that `f` is a primitive form.

    void qfb_nudupl(qfb_t r, const qfb_t f, fmpz_t D, fmpz_t L)
    # As for ``nucomp`` except that the form `f` is composed with itself.
    # We require that `f` is a primitive form.

    void qfb_pow_ui(qfb_t r, qfb_t f, fmpz_t D, ulong exp)
    # Compute the near reduced form `r` which is the result of composing the
    # principal form (identity) with `f` ``exp`` times.
    # We require `D` to be set to the discriminant of `f` and that `f` is a
    # primitive form.

    void qfb_pow(qfb_t r, qfb_t f, fmpz_t D, fmpz_t exp)
    # As per ``qfb_pow_ui``.

    void qfb_inverse(qfb_t r, qfb_t f)
    # Set `r` to the inverse of the binary quadratic form `f`.

    bint qfb_is_principal_form(qfb_t f, fmpz_t D)
    # Return `1` if `f` is the reduced principal form of discriminant `D`,
    # i.e. the identity in the form class group, else `0`.

    void qfb_principal_form(qfb_t f, fmpz_t D)
    # Set `f` to the principal form of discriminant `D`, i.e. the identity in
    # the form class group.

    bint qfb_is_primitive(qfb_t f)
    # Return `1` if `f` is primitive, i.e. the greatest common divisor of its
    # three coefficients is `1`. Otherwise the function returns `0`.

    void qfb_prime_form(qfb_t r, fmpz_t D, fmpz_t p)
    # Sets `r` to the unique prime `(p, b, c)` of discriminant `D`, i.e. with
    # `0 < b \leq p`. We require that `p` is a prime.

    int qfb_exponent_element(fmpz_t exponent, qfb_t f, fmpz_t n, ulong B1, ulong B2_sqrt)
    # Find the exponent of the element `f` in the form class group of forms of
    # discriminant `n`, doing a stage `1` with primes up to at least ``B1``
    # and a stage `2` for a single large prime up to at least the square of
    # ``B2_sqrt``. If the function fails to find the exponent it returns `0`,
    # otherwise the function returns `1` and ``exponent`` is set to the
    # exponent of `f`, i.e. the minimum power of `f` which gives the identity.
    # It is assumed that the form `f` is reduced. We require that ``iters``
    # is a power of `2` and that ``iters`` `\ge 1024`.
    # The function performs a stage `2` which stores up to `4\times`
    # ``iters`` binary quadratic forms, and `12\times` ``iters``
    # additional limbs of data in a hash table, where ``iters`` is the
    # square root of ``B2``.

    int qfb_exponent(fmpz_t exponent, fmpz_t n, ulong B1, ulong B2_sqrt, slong c)
    # Compute the exponent of the class group of discriminant `n`, doing
    # a stage `1` with primes up to at least ``B1`` and a stage `2` for
    # a single large prime up to at least the square of ``B2_sqrt``, and
    # with probability at least `1 - 2^{-c}`. If the prime limits are
    # exhausted without finding the exponent, the function returns `0`,
    # otherwise it returns `1` and ``exponent`` is set to the computed
    # exponent, i.e. the minimum power to which every element of the
    # class group has to be raised in order to get the identity.
    # The function performs a stage `2` which stores up to `4\times`
    # ``iters`` binary quadratic forms, and `12\times` ``iters``
    # additional limbs of data in a hash table, where ``iters`` is the
    # square root of ``B2``.
    # We use algorithm 8.1 of [Sut2007]_.

    int qfb_exponent_grh(fmpz_t exponent, fmpz_t n, ulong B1, ulong B2_sqrt)
    # Similar to ``qfb_exponent`` except that the bound ``c`` is
    # automatically generated such that the exponent is guaranteed to be
    # correct, if found, assuming the GRH, namely that the class group is
    # generated by primes less than `6\log^2(|n|)` as described in [BD1992]_.
