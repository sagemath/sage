# distutils: libraries = flint
# distutils: depends = flint/ca_ext.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void ca_ext_init_qqbar(ca_ext_t res, const qqbar_t x, ca_ctx_t ctx)
    # Initializes *res* and sets it to the algebraic constant *x*.

    void ca_ext_init_const(ca_ext_t res, calcium_func_code func, ca_ctx_t ctx)
    # Initializes *res* and sets it to the constant defined by *func*
    # (example: *func* = *CA_Pi* for `x = \pi`).

    void ca_ext_init_fx(ca_ext_t res, calcium_func_code func, const ca_t x, ca_ctx_t ctx)
    # Initializes *res* and sets it to the univariate function value `f(x)`
    # where *f* is defined by *func*  (example: *func* = *CA_Exp* for `e^x`).

    void ca_ext_init_fxy(ca_ext_t res, calcium_func_code func, const ca_t x, const ca_t y, ca_ctx_t ctx)
    # Initializes *res* and sets it to the bivariate function value `f(x, y)`
    # where *f* is defined by *func*  (example: *func* = *CA_Pow* for `x^y`).

    void ca_ext_init_fxn(ca_ext_t res, calcium_func_code func, ca_srcptr x, long nargs, ca_ctx_t ctx)
    # Initializes *res* and sets it to the multivariate function value
    # `f(x_1, \ldots, x_n)` where *f* is defined by *func* and *n* is
    # given by *nargs*.

    void ca_ext_init_set(ca_ext_t res, const ca_ext_t x, ca_ctx_t ctx)
    # Initializes *res* and sets it to a copy of *x*.

    void ca_ext_clear(ca_ext_t res, ca_ctx_t ctx)
    # Clears *res*.

    long ca_ext_nargs(const ca_ext_t x, ca_ctx_t ctx)
    # Returns the number of function arguments of *x*.
    # The return value is 0 for any algebraic constant and for any built-in
    # symbolic constant such as `\pi`.

    void ca_ext_get_arg(ca_t res, const ca_ext_t x, long i, ca_ctx_t ctx)
    # Sets *res* to argument *i* (indexed from zero) of *x*.
    # This calls *flint_abort* if *i* is out of range.

    unsigned long ca_ext_hash(const ca_ext_t x, ca_ctx_t ctx)
    # Returns a hash of the structural representation of *x*.

    int ca_ext_equal_repr(const ca_ext_t x, const ca_ext_t y, ca_ctx_t ctx)
    # Tests *x* and *y* for structural equality, returning 0 (false) or 1 (true).

    int ca_ext_cmp_repr(const ca_ext_t x, const ca_ext_t y, ca_ctx_t ctx)
    # Compares the representations of *x* and *y* in a canonical sort order,
    # returning -1, 0 or 1. This only performs a structural comparison
    # of the symbolic representations; the return value does not say
    # anything meaningful about the numbers represented by *x* and *y*.

    void ca_ext_print(const ca_ext_t x, ca_ctx_t ctx)
    # Prints a description of *x* to standard output.

    void ca_ext_get_acb_raw(acb_t res, ca_ext_t x, long prec, ca_ctx_t ctx)
    # Sets *res* to an enclosure of the numerical value of *x*.
    # A working precision of *prec* bits is used for the evaluation,
    # without adaptive refinement.

    void ca_ext_cache_init(ca_ext_cache_t cache, ca_ctx_t ctx)
    # Initializes *cache* for use.

    void ca_ext_cache_clear(ca_ext_cache_t cache, ca_ctx_t ctx)
    # Clears *cache*, freeing the memory allocated internally.

    ca_ext_ptr ca_ext_cache_insert(ca_ext_cache_t cache, const ca_ext_t x, ca_ctx_t ctx)
    # Adds *x* to *cache* without duplication. If a structurally identical
    # instance already exists in *cache*, a pointer to that instance is returned.
    # Otherwise, a copy of *x* is inserted into *cache* and a pointer to that new
    # instance is returned.
