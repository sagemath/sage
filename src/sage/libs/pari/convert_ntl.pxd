"""
Convert PARI objects to/from NTL objects

AUTHORS:

- Vincent Delecroix (2024)
"""
#*****************************************************************************
#       Copyright (C) 2024 Vincent Delecroix <20100.delecroix@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from cysignals.signals cimport sig_on

from cypari2.types cimport GEN
from cypari2.gen cimport Gen
from cypari2.paridecl cimport *
from cypari2.stack cimport new_gen

from .convert_gmp cimport _new_GEN_from_mpz_t

from sage.libs.gmp.types cimport mpz_t
from sage.libs.gmp.mpz cimport mpz_init, mpz_clear
from sage.libs.pari.convert_gmp cimport INT_to_mpz
from sage.libs.ntl.types cimport *
from sage.libs.ntl.convert cimport mpz_to_ZZ, ZZ_to_mpz


cdef Gen new_gen_from_ZZ(ZZ_c * x)


cdef inline void INT_to_ZZ(ZZ_c * x, GEN g):
    r"""
    Set ``x`` to the value of the pari INT ``g``
    """
    # TODO: we should not go through mpz here but implement a direct conversion
    cdef mpz_t tmp
    mpz_init(tmp)
    INT_to_mpz(tmp, g)
    mpz_to_ZZ(x, tmp)
    mpz_clear(tmp)


cdef inline GEN _new_GEN_from_ZZ(ZZ_c * x):
    r"""
    Create a new PARI ``t_INT`` from a NTL `ZZ``
    """
    # TODO: we should not go through mpz here but implement a direct conversion
    cdef mpz_t tmp
    mpz_init(tmp)
    ZZ_to_mpz(tmp, x)
    cdef GEN g = _new_GEN_from_mpz_t(tmp)
    mpz_clear(tmp)
    return g

# TODO
# cdef inline void FFp_to_ZZ_p(ZZ_p x, ZZ_pContext_c * p_context, GEN g):
#     r"""
#     Set ``x`` to the value of the pari Fp element ``g`` given the NTL context ``p_context``
#    """
#    pass

# TODO
# cdef inline void FFELT_to_ZZ_pE(ZZ_pE_c x, ZZ_pEContext_c * p_context, ZZ_pEContext_c * pE_context, GEN g):
#     r"""
#     Set ``x`` to the value of the pari FFELT ``g`` given the NTL contexts ``p_context`` and ``pE_context``.
#     """
#     pass
