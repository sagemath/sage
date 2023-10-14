# distutils: libraries = flint
# distutils: depends = flint/fq_vec.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    fq_struct * _fq_vec_init(slong len, const fq_ctx_t ctx)
    # Returns an initialised vector of ``fq``'s of given length.

    void _fq_vec_clear(fq_struct * vec, slong len, const fq_ctx_t ctx)
    # Clears the entries of ``(vec, len)`` and frees the space allocated
    # for ``vec``.

    void _fq_vec_randtest(fq_struct * f, flint_rand_t state, slong len, const fq_ctx_t ctx)
    # Sets the entries of a vector of the given length to elements of
    # the finite field.

    int _fq_vec_fprint(FILE * file, const fq_struct * vec, slong len, const fq_ctx_t ctx)
    # Prints the vector of given length to the stream ``file``. The
    # format is the length followed by two spaces, then a space separated
    # list of coefficients. If the length is zero, only `0` is printed.
    # In case of success, returns a positive value.  In case of failure,
    # returns a non-positive value.

    int _fq_vec_print(const fq_struct * vec, slong len, const fq_ctx_t ctx)
    # Prints the vector of given length to ``stdout``.
    # For further details, see ``_fq_vec_fprint()``.

    void _fq_vec_set(fq_struct * vec1, const fq_struct * vec2, slong len2, const fq_ctx_t ctx)
    # Makes a copy of ``(vec2, len2)`` into ``vec1``.

    void _fq_vec_swap(fq_struct * vec1, fq_struct * vec2, slong len2, const fq_ctx_t ctx)
    # Swaps the elements in ``(vec1, len2)`` and ``(vec2, len2)``.

    void _fq_vec_zero(fq_struct * vec, slong len, const fq_ctx_t ctx)
    # Zeros the entries of ``(vec, len)``.

    void _fq_vec_neg(fq_struct * vec1, const fq_struct * vec2, slong len2, const fq_ctx_t ctx)
    # Negates ``(vec2, len2)`` and places it into ``vec1``.

    bint _fq_vec_equal(const fq_struct * vec1, const fq_struct * vec2, slong len, const fq_ctx_t ctx)
    # Compares two vectors of the given length and returns `1` if they are
    # equal, otherwise returns `0`.

    bint _fq_vec_is_zero(const fq_struct * vec, slong len, const fq_ctx_t ctx)
    # Returns `1` if ``(vec, len)`` is zero, and `0` otherwise.

    void _fq_vec_add(fq_struct * res, const fq_struct * vec1, const fq_struct * vec2, slong len2, const fq_ctx_t ctx)
    # Sets ``(res, len2)`` to the sum of ``(vec1, len2)``
    # and ``(vec2, len2)``.

    void _fq_vec_sub(fq_struct * res, const fq_struct * vec1, const fq_struct * vec2, slong len2, const fq_ctx_t ctx)
    # Sets ``(res, len2)`` to ``(vec1, len2)`` minus ``(vec2, len2)``.

    void _fq_vec_scalar_addmul_fq(fq_struct * vec1, const fq_struct * vec2, slong len2, const fq_t c, const fq_ctx_t ctx)
    # Adds ``(vec2, len2)`` times `c` to ``(vec1, len2)``, where
    # `c` is a ``fq_t``.

    void _fq_vec_scalar_submul_fq(fq_struct * vec1, const fq_struct * vec2, slong len2, const fq_t c, const fq_ctx_t ctx)
    # Subtracts ``(vec2, len2)`` times `c` from ``(vec1, len2)``,
    # where `c` is a ``fq_t``.

    void _fq_vec_dot(fq_t res, const fq_struct * vec1, const fq_struct * vec2, slong len2, const fq_ctx_t ctx)
    # Sets ``res`` to the dot product of (``vec1``, ``len``)
    # and (``vec2``, ``len``).
