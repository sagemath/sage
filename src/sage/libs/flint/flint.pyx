"""
Flint imports

TESTS:

Import this module::

    sage: import sage.libs.flint.flint

We verify that :trac:`6919` is correctly fixed::

    sage: R.<x> = PolynomialRing(ZZ)
    sage: A = 2^(2^17+2^15)
    sage: a = A * x^31
    sage: b = (A * x) * x^30
    sage: a == b
    True
"""

# cimport all .pxd files to make sure they compile
cimport sage.libs.flint.arith
cimport sage.libs.flint.fmpq_poly
cimport sage.libs.flint.fmpq
cimport sage.libs.flint.fmpz_mat
cimport sage.libs.flint.fmpz_mod_poly
cimport sage.libs.flint.fmpz_poly
cimport sage.libs.flint.fmpz
cimport sage.libs.flint.fmpz_vec
cimport sage.libs.flint.fq_nmod
cimport sage.libs.flint.fq
cimport sage.libs.flint.nmod_poly
cimport sage.libs.flint.nmod_vec
cimport sage.libs.flint.padic
cimport sage.libs.flint.types
cimport sage.libs.flint.ulong_extras


def free_flint_stack():
    _fmpz_cleanup_mpz_content()
