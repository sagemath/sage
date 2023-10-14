# distutils: libraries = flint
# distutils: depends = flint/nmod_vec.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    mp_ptr _nmod_vec_init(slong len)
    # Returns a vector of the given length. The entries are not necessarily
    # zero.

    void _nmod_vec_clear(mp_ptr vec)
    # Frees the memory used by the given vector.

    void _nmod_vec_randtest(mp_ptr vec, flint_rand_t state, slong len, nmod_t mod)
    # Sets ``vec`` to a random vector of the given length with entries
    # reduced modulo ``mod.n``.

    void _nmod_vec_set(mp_ptr res, mp_srcptr vec, slong len)
    # Copies ``len`` entries from the vector ``vec`` to ``res``.

    void _nmod_vec_zero(mp_ptr vec, slong len)
    # Zeros the given vector of the given length.

    void _nmod_vec_swap(mp_ptr a, mp_ptr b, slong length)
    # Swaps the vectors ``a`` and ``b`` of length `n` by actually
    # swapping the entries.

    void _nmod_vec_reduce(mp_ptr res, mp_srcptr vec, slong len, nmod_t mod)
    # Reduces the entries of ``(vec, len)`` modulo ``mod.n`` and set
    # ``res`` to the result.

    flint_bitcnt_t _nmod_vec_max_bits(mp_srcptr vec, slong len)
    # Returns the maximum number of bits of any entry in the vector.

    bint _nmod_vec_equal(mp_srcptr vec, mp_srcptr vec2, slong len)
    # Returns~`1` if ``(vec, len)`` is equal to ``(vec2, len)``,
    # otherwise returns~`0`.

    void _nmod_vec_print_pretty(mp_srcptr vec, slong len, nmod_t mod)
    # Pretty-prints ``vec`` to ``stdout``. A header is printed followed by the
    # vector enclosed in brackets. Each entry is right-aligned to the width of
    # the modulus written in decimal, and the entries are separated by spaces.
    # For example::
    # <length-12 integer vector mod 197>
    # [ 33 181 107  61  32  11  80 138  34 171  86 156]

    int _nmod_vec_fprint_pretty(FILE * file, mp_srcptr vec, slong len, nmod_t mod)
    # Same as ``_nmod_vec_print_pretty`` but printing to ``file``.

    int _nmod_vec_print(mp_srcptr vec, slong len, nmod_t mod)
    # Currently, same as ``_nmod_vec_print_pretty``.

    int _nmod_vec_fprint(FILE * f, mp_srcptr vec, slong len, nmod_t mod)
    # Currently, same as ``_nmod_vec_fprint_pretty``.

    void _nmod_vec_add(mp_ptr res, mp_srcptr vec1, mp_srcptr vec2, slong len, nmod_t mod)
    # Sets ``(res, len)`` to the sum of ``(vec1, len)``
    # and ``(vec2, len)``.

    void _nmod_vec_sub(mp_ptr res, mp_srcptr vec1, mp_srcptr vec2, slong len, nmod_t mod)
    # Sets ``(res, len)`` to the difference of ``(vec1, len)``
    # and ``(vec2, len)``.

    void _nmod_vec_neg(mp_ptr res, mp_srcptr vec, slong len, nmod_t mod)
    # Sets ``(res, len)`` to the negation of ``(vec, len)``.

    void _nmod_vec_scalar_mul_nmod(mp_ptr res, mp_srcptr vec, slong len, mp_limb_t c, nmod_t mod)
    # Sets ``(res, len)`` to ``(vec, len)`` multiplied by `c`. The element
    # `c` and all elements of `vec` are assumed to be less than `mod.n`.

    void _nmod_vec_scalar_mul_nmod_shoup(mp_ptr res, mp_srcptr vec, slong len, mp_limb_t c, nmod_t mod)
    # Sets ``(res, len)`` to ``(vec, len)`` multiplied by `c` using
    # :func:`n_mulmod_shoup`. `mod.n` should be less than `2^{\mathtt{FLINT\_BITS} - 1}`. `c`
    # and all elements of `vec` should be less than `mod.n`.

    void _nmod_vec_scalar_addmul_nmod(mp_ptr res, mp_srcptr vec, slong len, mp_limb_t c, nmod_t mod)
    # Adds ``(vec, len)`` times `c` to the vector ``(res, len)``. The element
    # `c` and all elements of `vec` are assumed to be less than `mod.n`.

    int _nmod_vec_dot_bound_limbs(slong len, nmod_t mod)
    # Returns the number of limbs (0, 1, 2 or 3) needed to represent the
    # unreduced dot product of two vectors of length ``len`` having entries
    # modulo ``mod.n``, assuming that ``len`` is nonnegative and that
    # ``mod.n`` is nonzero. The computed bound is tight. In other words,
    # this function returns the precise limb size of ``len`` times
    # ``(mod.n - 1) ^ 2``.

    mp_limb_t _nmod_vec_dot(mp_srcptr vec1, mp_srcptr vec2, slong len, nmod_t mod, int nlimbs)
    # Returns the dot product of (``vec1``, ``len``) and
    # (``vec2``, ``len``). The ``nlimbs`` parameter should be
    # 0, 1, 2 or 3, specifying the number of limbs needed to represent the
    # unreduced result.

    mp_limb_t _nmod_vec_dot_rev(mp_srcptr vec1, mp_srcptr vec2, slong len, nmod_t mod, int nlimbs)
    # The same as ``_nmod_vec_dot``, but reverses ``vec2``.

    mp_limb_t _nmod_vec_dot_ptr(mp_srcptr vec1, const mp_ptr * vec2, slong offset, slong len, nmod_t mod, int nlimbs)
    # Returns the dot product of (``vec1``, ``len``) and the values at
    # ``vec2[i][offset]``. The ``nlimbs`` parameter should be
    # 0, 1, 2 or 3, specifying the number of limbs needed to represent the
    # unreduced result.
