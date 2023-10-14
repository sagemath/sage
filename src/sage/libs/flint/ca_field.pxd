# distutils: libraries = flint
# distutils: depends = flint/ca_field.h

################################################################################
# This file is auto-generated. Do not modify by hand
################################################################################

from libc.stdio cimport FILE
from sage.libs.gmp.types cimport *
from sage.libs.mpfr.types cimport *
from sage.libs.flint.types cimport *

cdef extern from "flint_wrap.h":

    void ca_field_init_qq(ca_field_t K, ca_ctx_t ctx)
    # Initializes *K* to represent the trivial field `\mathbb{Q}`.

    void ca_field_init_nf(ca_field_t K, const qqbar_t x, ca_ctx_t ctx)
    # Initializes *K* to represent the algebraic number field `\mathbb{Q}(x)`.

    void ca_field_init_const(ca_field_t K, calcium_func_code func, ca_ctx_t ctx)
    # Initializes *K* to represent the field
    # `\mathbb{Q}(x)` where *x* is a builtin constant defined by
    # *func* (example: *func* = *CA_Pi* for `x = \pi`).

    void ca_field_init_fx(ca_field_t K, calcium_func_code func, const ca_t x, ca_ctx_t ctx)
    # Initializes *K* to represent the field
    # `\mathbb{Q}(a)` where `a = f(x)`, given a number *x* and a builtin
    # univariate function *func* (example: *func* = *CA_Exp* for `e^x`).

    void ca_field_init_fxy(ca_field_t K, calcium_func_code func, const ca_t x, const ca_t y, ca_ctx_t ctx)
    # Initializes *K* to represent the field
    # `\mathbb{Q}(a,b)` where `a = f(x, y)`.

    void ca_field_init_multi(ca_field_t K, slong len, ca_ctx_t ctx)
    # Initializes *K* to represent a multivariate field
    # `\mathbb{Q}(a_1, \ldots, a_n)` in *n*
    # extension numbers. The extension numbers must subsequently be
    # assigned one by one using :func:`ca_field_set_ext`.

    void ca_field_set_ext(ca_field_t K, slong i, ca_ext_srcptr x_index, ca_ctx_t ctx)
    # Sets the extension number at position *i* (here indexed from 0) of *K*
    # to the generator of the field with index *x_index* in *ctx*.
    # (It is assumed that the generating field is a univariate field.)
    # This only inserts a shallow reference: the field at index *x_index* must
    # be kept alive until *K* has been cleared.

    void ca_field_clear(ca_field_t K, ca_ctx_t ctx)
    # Clears the field *K*. This does not clear the individual extension
    # numbers, which are only held as references.

    void ca_field_print(const ca_field_t K, ca_ctx_t ctx)
    # Prints a description of the field *K* to standard output.

    void ca_field_build_ideal(ca_field_t K, ca_ctx_t ctx)
    # Given *K* with assigned extension numbers,
    # builds the reduction ideal in-place.

    void ca_field_build_ideal_erf(ca_field_t K, ca_ctx_t ctx)
    # Builds relations for error functions present among the extension
    # numbers in *K*. This heuristic adds relations that are consequences
    # of the functional equations
    # `\operatorname{erf}(x) = -\operatorname{erf}(-x)`,
    # `\operatorname{erfc}(x) = 1-\operatorname{erf}(x)`,
    # `\operatorname{erfi}(x) = -i\operatorname{erf}(ix)`.

    int ca_field_cmp(const ca_field_t K1, const ca_field_t K2, ca_ctx_t ctx)
    # Compares the field objects *K1* and *K2* in a canonical sort order,
    # returning -1, 0 or 1. This only performs a lexicographic comparison
    # of the representations of *K1* and *K2*; the return value does not say
    # anything meaningful about the relative structures of *K1* and *K2*
    # as mathematical fields.

    void ca_field_cache_init(ca_field_cache_t cache, ca_ctx_t ctx)
    # Initializes *cache* for use.

    void ca_field_cache_clear(ca_field_cache_t cache, ca_ctx_t ctx)
    # Clears *cache*, freeing the memory allocated internally.
    # This does not clear the individual extension
    # numbers, which are only held as references.

    ca_field_ptr ca_field_cache_insert_ext(ca_field_cache_t cache, ca_ext_struct ** x, slong len, ca_ctx_t ctx)
    # Adds the field defined by the length-*len* list of extension numbers *x*
    # to *cache* without duplication. If such a field already exists in *cache*,
    # a pointer to that instance is returned. Otherwise, a field with
    # extension numbers *x* is inserted into *cache* and a pointer to that
    # new instance is returned. Upon insertion of a new field, the
    # reduction ideal is constructed via :func:`ca_field_build_ideal`.
