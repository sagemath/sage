r"""
Fast binary operations for basic types
"""

from .types cimport mpz_t, mpq_t
from .mpz cimport mpz_set, mpz_add, mpz_mul
from .mpq cimport mpq_canonicalize, mpq_numref, mpq_denref, mpq_add

cdef inline void mpq_add_z(mpq_t res, mpq_t op1, mpz_t op2):
    mpz_mul(mpq_numref(res), mpq_denref(op1), op2)
    mpz_add(mpq_numref(res), mpq_numref(res), mpq_numref(op1))
    mpz_set(mpq_denref(res), mpq_denref(op1))

cdef inline void mpq_div_zz(mpq_t res, mpz_t op1, mpz_t op2):
    mpz_set(mpq_numref(res), op1)
    mpz_set(mpq_denref(res), op2)
    mpq_canonicalize(res)

cdef inline void mpq_mul_z(mpq_t res, mpq_t op1, mpz_t op2):
    # (A/B) * C = (A/C) * B
    mpq_div_zz(res, op2, mpq_denref(op1))
    mpz_mul(mpq_numref(res), mpq_numref(res), mpq_numref(op1))
