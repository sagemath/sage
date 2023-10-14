# distutils: libraries = flint
# distutils: depends = flint/nmod.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void nmod_init(nmod_t * mod, mp_limb_t n)
    # Initialises the given ``nmod_t`` structure for reduction modulo `n`
    # with a precomputed inverse.

    mp_limb_t _nmod_add(mp_limb_t a, mp_limb_t b, nmod_t mod)
    # Returns `a + b` modulo ``mod.n``. It is assumed that ``mod`` is
    # no more than ``FLINT_BITS - 1`` bits. It is assumed that `a` and `b`
    # are already reduced modulo ``mod.n``.

    mp_limb_t nmod_add(mp_limb_t a, mp_limb_t b, nmod_t mod)
    # Returns `a + b` modulo ``mod.n``. No assumptions are made about
    # ``mod.n``. It is assumed that `a` and `b` are already reduced
    # modulo ``mod.n``.

    mp_limb_t _nmod_sub(mp_limb_t a, mp_limb_t b, nmod_t mod)
    # Returns `a - b` modulo ``mod.n``. It is assumed that ``mod``
    # is no more than ``FLINT_BITS - 1`` bits. It is assumed that
    # `a` and `b` are already reduced modulo ``mod.n``.

    mp_limb_t nmod_sub(mp_limb_t a, mp_limb_t b, nmod_t mod)
    # Returns `a - b` modulo ``mod.n``. No assumptions are made about
    # ``mod.n``. It is assumed that `a` and `b` are already reduced
    # modulo ``mod.n``.

    mp_limb_t nmod_neg(mp_limb_t a, nmod_t mod)
    # Returns `-a` modulo ``mod.n``. It is assumed that `a` is already
    # reduced modulo ``mod.n``, but no assumptions are made about the
    # latter.

    mp_limb_t nmod_mul(mp_limb_t a, mp_limb_t b, nmod_t mod)
    # Returns `ab` modulo ``mod.n``. No assumptions are made about
    # ``mod.n``. It is assumed that `a` and `b` are already reduced
    # modulo ``mod.n``.

    mp_limb_t _nmod_mul_fullword(mp_limb_t a, mp_limb_t b, nmod_t mod)
    # Returns `ab` modulo ``mod.n``. Requires that ``mod.n`` is exactly
    # ``FLINT_BITS`` large. It is assumed that `a` and `b` are already
    # reduced modulo ``mod.n``.

    mp_limb_t nmod_inv(mp_limb_t a, nmod_t mod)
    # Returns `a^{-1}` modulo ``mod.n``. The inverse is assumed to exist.

    mp_limb_t nmod_div(mp_limb_t a, mp_limb_t b, nmod_t mod)
    # Returns `ab^{-1}` modulo ``mod.n``. The inverse of `b` is assumed to
    # exist. It is assumed that `a` is already reduced modulo ``mod.n``.

    mp_limb_t nmod_pow_ui(mp_limb_t a, ulong e, nmod_t mod)
    # Returns `a^e` modulo ``mod.n``. No assumptions are made about
    # ``mod.n``. It is assumed that `a` is already reduced
    # modulo ``mod.n``.

    mp_limb_t nmod_pow_fmpz(mp_limb_t a, const fmpz_t e, nmod_t mod)
    # Returns `a^e` modulo ``mod.n``. No assumptions are made about
    # ``mod.n``. It is assumed that `a` is already reduced
    # modulo ``mod.n`` and that `e` is not negative.

    void nmod_discrete_log_pohlig_hellman_init(nmod_discrete_log_pohlig_hellman_t L)
    # Initialize ``L``. Upon initialization ``L`` is not ready for computation.

    void nmod_discrete_log_pohlig_hellman_clear(nmod_discrete_log_pohlig_hellman_t L)
    # Free any space used by ``L``.

    double nmod_discrete_log_pohlig_hellman_precompute_prime(nmod_discrete_log_pohlig_hellman_t L, mp_limb_t p)
    # Configure ``L`` for discrete logarithms modulo ``p`` to an internally chosen base. It is assumed that ``p`` is prime.
    # The return is an estimate on the number of multiplications needed for one run.

    mp_limb_t nmod_discrete_log_pohlig_hellman_primitive_root(const nmod_discrete_log_pohlig_hellman_t L)
    # Return the internally stored base.

    ulong nmod_discrete_log_pohlig_hellman_run(const nmod_discrete_log_pohlig_hellman_t L, mp_limb_t y)
    # Return the logarithm of ``y`` with respect to the internally stored base. ``y`` is expected to be reduced modulo the ``p``.
    # The function is undefined if the logarithm does not exist.
