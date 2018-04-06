# distutils: extra_compile_args = LINBOX_CFLAGS
# distutils: libraries = LINBOX_LIBRARIES
# distutils: library_dirs = LINBOX_LIBDIR
# distutils: language = c++
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

- ``unsigned long linbox_fmpz_mat_rank(fmpz_mat_t A)``: return the rank of the
  matrix ``A``

- ``void linbox_fmpz_mat_det(fmpz_t det, fmpz_mat_t A)``: set ``det`` to the
  determinat of the square matrix ``A``
"""
#*****************************************************************************
#       Copyright (C) 2007 Martin Albrecht
#       Copyright (C) 2008 Clement Pernet
#       Copyright (C) 2017 Vincent Delecroix
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.libs.gmp.types cimport mpz_t, mpz_srcptr, mpz_ptr
from sage.libs.gmp.mpz cimport mpz_set
from sage.libs.flint.types cimport fmpz, fmpz_t
from sage.libs.flint.fmpz cimport fmpz_get_mpz, fmpz_set_mpz
from sage.libs.flint.fmpz_mat cimport fmpz_mat_entry, fmpz_mat_nrows, fmpz_mat_ncols
from sage.libs.flint.fmpz_poly cimport fmpz_poly_set_coeff_mpz, fmpz_poly_fit_length, _fmpz_poly_set_length, fmpz_poly_one

cdef extern from "givaro/givconfig.h":
    pass
cdef extern from "linbox/linbox-config.h":
    pass

cdef extern from "gmp++/gmp++.h":
    cdef cppclass GivaroInteger "Givaro::Integer":
        mpz_ptr get_mpz()
        mpz_srcptr get_mpz_const()

cdef extern from "linbox/matrix/dense-matrix.h":
    ## template <class _Field>
    ## using DenseMatrix = BlasMatrix<_Field> ;
    ##
    # indirect include from linbox/matrix/densematrix/blas-matrix.h
    ##
    ## template <class _Field, class _Storage>
    ## class BlasMatrix
    cdef cppclass LinBoxIntegerDenseMatrix "LinBox::DenseMatrix<Givaro::ZRing<Givaro::Integer>>":
        ctypedef GivaroIntegerRing Field
        ctypedef GivaroInteger Element
        size_t rowdim()
        size_t coldim()
        LinBoxIntegerDenseMatrix(Field &F, size_t m, size_t n)
        void setEntry(size_t i, size_t j, Element &a)
        Element &getEntry(size_t i, size_t j)

cdef extern from "givaro/zring.h":
    cdef cppclass GivaroIntegerRing "Givaro::ZRing<Givaro::Integer>":
        ctypedef GivaroInteger Element


cdef extern from "givaro/givpoly1.h":
    ## template < typename T, typename A=std::allocator<T> >
    ## class givvector : public __GIV_STANDARD_VECTOR<T,A>
    cdef cppclass givvector "Givaro::givvector" [T,ALLOCATOR=*]:
        T& operator[](size_t i)
        size_t size()

cdef extern from "linbox/ring/givaro-polynomial.h":
    ## template <class Domain, class StorageTag= Givaro::Dense>
    ## class GivPolynomialRing : public Givaro::Poly1FactorDom< Domain,StorageTag>
    cdef cppclass LinBoxIntegerPolynomialRing "LinBox::GivPolynomialRing<Givaro::ZRing<Givaro::Integer>, Givaro::Dense>":
        ctypedef givvector[GivaroInteger] Element
        ctypedef givvector[GivaroInteger] Polynomial

cdef extern from "linbox/matrix/matrix-domain.h":
    ## template <class Field_ >
    ## class MatrixDomain : public MVProductDomain<Field_> {
    cdef cppclass LinBoxIntegerDenseMatrixDomain "LinBox::MatrixDomain<Givaro::ZRing<Givaro::Integer>>":
        LinBoxIntegerDenseMatrixDomain(GivaroIntegerRing&)
        LinBoxIntegerDenseMatrix& mul(LinBoxIntegerDenseMatrix ans,
                                      LinBoxIntegerDenseMatrix left,
                                      LinBoxIntegerDenseMatrix right)

cdef extern from "linbox/solutions/charpoly.h":
    ## template<class Blackbox, class Polynomial>
    ## Polynomial &charpoly (Polynomial & P, const Blackbox & A)
    LinBoxIntegerPolynomialRing.Element& LinBoxIntegerDense_charpoly "LinBox::charpoly" (LinBoxIntegerPolynomialRing.Element&, LinBoxIntegerDenseMatrix&)

cdef extern from "linbox/solutions/minpoly.h":
    ## template<class Polynomial, class Blackbox>
    ## Polynomial &minpoly (Polynomial & P, const Blackbox & A)
    LinBoxIntegerPolynomialRing.Element& LinBoxIntegerDense_minpoly "LinBox::minpoly" (LinBoxIntegerPolynomialRing.Element&, LinBoxIntegerDenseMatrix&)

cdef extern from "linbox/solutions/rank.h":
    ## template <class Blackbox, class Method, class DomainCategory>
    ## inline unsigned long &rank (unsigned long &r, const Blackbox &A,
    ##                             const DomainCategory &tag, const Method &M);
    unsigned long & LinBoxIntegerDense_rank "LinBox::rank" (unsigned long &, LinBoxIntegerDenseMatrix)

cdef extern from "linbox/solutions/det.h":
    GivaroInteger& LinBoxIntegerDense_det "LinBox::det" (GivaroInteger&, LinBoxIntegerDenseMatrix)



###############################################################################
# end of LinBox declarations -- begining of the code                          #
###############################################################################


# set the entries of A from m (no allocation performed)
# NOTE: this function is not part of the interface (ie the .pxd file) to keep the
# module C-compatible
cdef void fmpz_mat_get_linbox(LinBoxIntegerDenseMatrix& A, fmpz_mat_t m):
    cdef size_t i,j
    cdef GivaroInteger t

    for i in range(fmpz_mat_nrows(m)):
        for j in range(fmpz_mat_ncols(m)):
            fmpz_get_mpz(t.get_mpz(), fmpz_mat_entry(m, i, j))
            A.setEntry(i, j, t)


# set the entries of m from A (no allocation performed)
# NOTE: this function is not part of the interface (ie the .pxd file) to keep the
# module C-compatible
cdef void fmpz_mat_set_linbox(fmpz_mat_t m, LinBoxIntegerDenseMatrix& A):
    cdef size_t i,j
    for i in range(A.rowdim()):
        for j in range(A.coldim()):
            fmpz_set_mpz(fmpz_mat_entry(m, i, j), A.getEntry(i, j).get_mpz_const())


# set the entries of the polynomial p from q (no allocation performed)
# NOTE: this function is not part of the interface (ie the .pxd file) to keep the
# module C-compatible
cdef void fmpz_poly_set_linbox(fmpz_poly_t p, LinBoxIntegerPolynomialRing.Element& q):
    cdef size_t i

    fmpz_poly_fit_length(p, q.size())

    for i in range(q.size()):
        fmpz_poly_set_coeff_mpz(p, i, q[i].get_mpz_const())

    _fmpz_poly_set_length(p, q.size())


# set C <- A * B
cdef void linbox_fmpz_mat_mul(fmpz_mat_t C, fmpz_mat_t A, fmpz_mat_t B):
    cdef GivaroIntegerRing ZZ
    cdef LinBoxIntegerDenseMatrix *LBA
    cdef LinBoxIntegerDenseMatrix *LBB
    cdef LinBoxIntegerDenseMatrix *LBC
    cdef LinBoxIntegerDenseMatrixDomain * MD

    LBA = new LinBoxIntegerDenseMatrix(ZZ, fmpz_mat_nrows(A), fmpz_mat_ncols(A))
    fmpz_mat_get_linbox(LBA[0], A)

    LBB = new LinBoxIntegerDenseMatrix(ZZ, fmpz_mat_nrows(B), fmpz_mat_ncols(B))
    fmpz_mat_get_linbox(LBB[0], B)

    LBC = new LinBoxIntegerDenseMatrix(ZZ, fmpz_mat_nrows(A), fmpz_mat_ncols(B))

    MD = new LinBoxIntegerDenseMatrixDomain(ZZ)
    MD.mul(LBC[0], LBA[0], LBB[0])

    del MD

    fmpz_mat_set_linbox(C, LBC[0])


# set cp to the characteristic polynomial of A
cdef void linbox_fmpz_mat_charpoly(fmpz_poly_t cp, fmpz_mat_t A):
    cdef GivaroIntegerRing ZZ
    cdef LinBoxIntegerDenseMatrix * LBA
    cdef LinBoxIntegerPolynomialRing.Element m_A

    # FIXME: bug in LinBox
    # see https://github.com/linbox-team/linbox/issues/51
    if fmpz_mat_nrows(A) == 0:
        fmpz_poly_one(cp)
        return

    LBA = new LinBoxIntegerDenseMatrix(ZZ, fmpz_mat_nrows(A), fmpz_mat_ncols(A))
    fmpz_mat_get_linbox(LBA[0], A)
    LinBoxIntegerDense_charpoly(m_A, LBA[0])
    fmpz_poly_set_linbox(cp, m_A)

    del LBA


# set mp to the minimal polynomial of A
cdef void linbox_fmpz_mat_minpoly(fmpz_poly_t mp, fmpz_mat_t A):
    cdef GivaroIntegerRing ZZ
    cdef LinBoxIntegerDenseMatrix * LBA
    cdef LinBoxIntegerPolynomialRing.Element m_A

    # FIXME: bug in LinBox
    # see https://github.com/linbox-team/linbox/issues/51
    if fmpz_mat_nrows(A) == 0:
        fmpz_poly_one(mp)
        return

    LBA = new LinBoxIntegerDenseMatrix(ZZ, fmpz_mat_nrows(A), fmpz_mat_ncols(A))
    fmpz_mat_get_linbox(LBA[0], A)
    LinBoxIntegerDense_minpoly(m_A, LBA[0])
    fmpz_poly_set_linbox(mp, m_A)

    del LBA


# return the rank of A
cdef unsigned long linbox_fmpz_mat_rank(fmpz_mat_t A):
    cdef GivaroIntegerRing ZZ
    cdef LinBoxIntegerDenseMatrix * LBA
    cdef unsigned long r = 0

    LBA = new LinBoxIntegerDenseMatrix(ZZ, fmpz_mat_nrows(A), fmpz_mat_ncols(A))
    fmpz_mat_get_linbox(LBA[0], A)
    LinBoxIntegerDense_rank(r, LBA[0])

    del LBA

    return r


# set det to the determinant of A
cdef void linbox_fmpz_mat_det(fmpz_t det, fmpz_mat_t A):
    cdef GivaroIntegerRing ZZ
    cdef LinBoxIntegerDenseMatrix * LBA
    cdef GivaroInteger d

    LBA = new LinBoxIntegerDenseMatrix(ZZ, fmpz_mat_nrows(A), fmpz_mat_ncols(A))
    fmpz_mat_get_linbox(LBA[0], A)
    LinBoxIntegerDense_det(d, LBA[0])
    fmpz_set_mpz(det, d.get_mpz_const())

    del LBA
