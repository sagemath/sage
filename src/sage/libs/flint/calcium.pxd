# distutils: libraries = flint
# distutils: depends = flint/calcium.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    const char * calcium_version()
    # Returns a pointer to the version of the library as a string ``X.Y.Z``.

    ulong calcium_fmpz_hash(const fmpz_t x)
    # Hash function for integers. The algorithm may change;
    # presently, this simply extracts the low word (with sign).

    void calcium_stream_init_file(calcium_stream_t out, FILE * fp)
    # Initializes the stream *out* for writing to the file *fp*.
    # The file can be *stdout*, *stderr*, or any file opened for writing
    # by the user.

    void calcium_stream_init_str(calcium_stream_t out)
    # Initializes the stream *out* for writing to a string in memory.
    # When finished, the user should free the string (the *s* member
    # of *out* with ``flint_free()``).

    void calcium_write(calcium_stream_t out, const char * s)
    # Writes the string *s* to *out*.

    void calcium_write_free(calcium_stream_t out, char * s)
    # Writes *s* to *out* and then frees *s* by calling ``flint_free()``.

    void calcium_write_si(calcium_stream_t out, slong x)
    void calcium_write_fmpz(calcium_stream_t out, const fmpz_t x)
    # Writes the integer *x* to *out*.

    void calcium_write_arb(calcium_stream_t out, const arb_t z, slong digits, ulong flags)
    void calcium_write_acb(calcium_stream_t out, const acb_t z, slong digits, ulong flags)
    # Writes the Arb number *z* to *out*, showing *digits*
    # digits and with the display style specified by *flags*
    # (``ARB_STR_NO_RADIUS``, etc.).
