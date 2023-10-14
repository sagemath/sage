# distutils: libraries = flint
# distutils: depends = flint/nf.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void nf_init(nf_t nf, const fmpq_poly_t pol)
    # Perform basic initialisation of a number field (for element arithmetic)
    # given a defining polynomial over `\mathbb{Q}`.

    void nf_clear(nf_t nf)
    # Release resources used by a number field object. The object will need
    # initialisation again before it can be used.
