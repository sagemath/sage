r"""
Interface between flint matrices and linbox

This module only contains C++ code (and the interface is fully C
compatible). It basically contains what used to be in the LinBox
source code under interfaces/sage/linbox-sage.C written by M. Albrecht
and C. Pernet. The functions available are:

- ``void  linbox_fmpz_mat_mul(fmpz_mat_t C, fmpz_mat_t A, fmpz_mat_t B)``: set
  ``C`` to be the result of the multiplication ``A * B``

- ``void linbox_fmpz_mat_charpoly(fmpz_poly_t cp, fmpz_mat_t A)``: set ``cp``
  to be the characteristic polynomial of the square matrix ``A``

- ``void  linbox_fmpz_mat_minpoly(fmpz_poly_t mp, fmpz_mat_t A)``: set ``mp``
  to be the minimal polynomial of the square matrix ``A``

- ``size_t linbox_fmpz_mat_rank(fmpz_mat_t A)``: return the rank of the
  matrix ``A``

- ``void linbox_fmpz_mat_det(fmpz_t det, fmpz_mat_t A)``: set ``det`` to the
  determinant of the square matrix ``A``
"""
# ****************************************************************************
#       Copyright (C) 2007 Martin Albrecht
#       Copyright (C) 2008 Clement Pernet
#       Copyright (C) 2017-2018 Vincent Delecroix
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.libs.gmp.types cimport mpz_t
from sage.libs.flint.types cimport fmpz_t
from sage.libs.flint.fmpz cimport fmpz_get_mpz, fmpz_set_mpz
from sage.libs.flint.fmpz_mat cimport fmpz_mat_entry, fmpz_mat_nrows, fmpz_mat_ncols
from sage.libs.flint.fmpz_poly cimport fmpz_poly_fit_length, _fmpz_poly_set_length
from sage.libs.flint.fmpz_poly_sage cimport fmpz_poly_set_coeff_mpz

cimport sage.libs.linbox.givaro as givaro
cimport sage.libs.linbox.linbox as linbox
from sage.libs.linbox.linbox cimport PolynomialRing_integer


cdef void fmpz_mat_get_linbox(linbox.DenseMatrix_integer& A, fmpz_mat_t m) noexcept:
    r"""
    Set the entries of A from m (no allocation performed).

    NOTE: this function does not appear in the Cython header
    (the .pxd file) in order to keep the header C-compatible
    """
    cdef size_t i,j
    cdef givaro.Integer t

    for i in range(fmpz_mat_nrows(m)):
        for j in range(fmpz_mat_ncols(m)):
            fmpz_get_mpz(t.get_mpz(), fmpz_mat_entry(m, i, j))
            A.setEntry(i, j, t)


cdef void fmpz_mat_set_linbox(fmpz_mat_t m, linbox.DenseMatrix_integer& A) noexcept:
    r"""
    Set the entries of m from A (no allocation performed).

    NOTE: this function does not appear in the Cython header
    (the .pxd file) in order to keep the header C-compatible
    """
    cdef size_t i,j
    for i in range(A.rowdim()):
        for j in range(A.coldim()):
            fmpz_set_mpz(fmpz_mat_entry(m, i, j), A.getEntry(i, j).get_mpz_const())


cdef void fmpz_poly_set_linbox(fmpz_poly_t p, PolynomialRing_integer.Element& q) noexcept:
    r"""
    Set the entries of the polynomial p from q (no allocation performed).

    NOTE: this function does not appear in the Cython header
    (the .pxd file) in order to keep the header C-compatible
    """
    cdef size_t i
    cdef mpz_t tmp

    fmpz_poly_fit_length(p, q.size())

    for i in range(q.size()):
        tmp = q[i].get_mpz_const()
        fmpz_poly_set_coeff_mpz(p, i, tmp)

    _fmpz_poly_set_length(p, q.size())


cdef void linbox_fmpz_mat_mul(fmpz_mat_t C, fmpz_mat_t A, fmpz_mat_t B) noexcept:
    r"""
    Set C to be A * B.
    """
    cdef givaro.ZRing ZZ
    cdef linbox.DenseMatrix_integer *LBA
    cdef linbox.DenseMatrix_integer *LBB
    cdef linbox.DenseMatrix_integer *LBC
    cdef linbox.MatrixDomain_integer * MD

    LBA = new linbox.DenseMatrix_integer(ZZ, fmpz_mat_nrows(A), fmpz_mat_ncols(A))
    fmpz_mat_get_linbox(LBA[0], A)

    LBB = new linbox.DenseMatrix_integer(ZZ, fmpz_mat_nrows(B), fmpz_mat_ncols(B))
    fmpz_mat_get_linbox(LBB[0], B)

    LBC = new linbox.DenseMatrix_integer(ZZ, fmpz_mat_nrows(A), fmpz_mat_ncols(B))

    MD = new linbox.MatrixDomain_integer(ZZ)
    MD.mul(LBC[0], LBA[0], LBB[0])

    del MD

    fmpz_mat_set_linbox(C, LBC[0])


cdef void linbox_fmpz_mat_charpoly(fmpz_poly_t cp, fmpz_mat_t A) noexcept:
    r"""
    Set cp to the characteristic polynomial of A.
    """
    cdef givaro.ZRing ZZ
    cdef linbox.DenseMatrix_integer * LBA
    cdef linbox.DensePolynomial_integer * m_A

    LBA = new linbox.DenseMatrix_integer(ZZ, fmpz_mat_nrows(A), fmpz_mat_ncols(A))
    fmpz_mat_get_linbox(LBA[0], A)
    m_A = new linbox.DensePolynomial_integer(ZZ, fmpz_mat_nrows(A))
    linbox.charpoly(m_A[0], LBA[0])
    fmpz_poly_set_linbox(cp, m_A[0])

    del LBA
    del m_A


cdef void linbox_fmpz_mat_minpoly(fmpz_poly_t mp, fmpz_mat_t A) noexcept:
    r"""
    Set mp to the minimal polynomial of A.
    """
    cdef givaro.ZRing ZZ
    cdef linbox.DenseMatrix_integer * LBA
    cdef linbox.DensePolynomial_integer * m_A

    LBA = new linbox.DenseMatrix_integer(ZZ, fmpz_mat_nrows(A), fmpz_mat_ncols(A))
    m_A = new linbox.DensePolynomial_integer(ZZ)
    fmpz_mat_get_linbox(LBA[0], A)
    linbox.minpoly(m_A[0], LBA[0])
    fmpz_poly_set_linbox(mp, m_A[0])

    del LBA
    del m_A


cdef size_t linbox_fmpz_mat_rank(fmpz_mat_t A) noexcept:
    r"""
    Return the rank of A
    """
    cdef givaro.ZRing ZZ
    cdef linbox.DenseMatrix_integer * LBA
    cdef size_t r = 0

    LBA = new linbox.DenseMatrix_integer(ZZ, fmpz_mat_nrows(A), fmpz_mat_ncols(A))
    fmpz_mat_get_linbox(LBA[0], A)
    linbox.rank(r, LBA[0])

    del LBA

    return r


cdef void linbox_fmpz_mat_det(fmpz_t det, fmpz_mat_t A) noexcept:
    r"""
    Set det to the determinant of A.
    """
    cdef givaro.ZRing ZZ
    cdef linbox.DenseMatrix_integer * LBA
    cdef givaro.Integer d

    LBA = new linbox.DenseMatrix_integer(ZZ, fmpz_mat_nrows(A), fmpz_mat_ncols(A))
    fmpz_mat_get_linbox(LBA[0], A)
    linbox.det(d, LBA[0])
    fmpz_set_mpz(det, d.get_mpz_const())

    del LBA
