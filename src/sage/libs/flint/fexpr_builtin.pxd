# distutils: libraries = flint
# distutils: depends = flint/fexpr_builtin.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    slong fexpr_builtin_lookup(const char * s)
    # Returns the internal index used to encode the builtin symbol
    # with name *s* in expressions. If *s* is not the name of a builtin
    # symbol, returns -1.

    const char * fexpr_builtin_name(slong n)
    # Returns a read-only pointer for a string giving the name of the
    # builtin symbol with index *n*.

    slong fexpr_builtin_length()
    # Returns the number of builtin symbols.
