# sage_setup: distribution = sagemath-linbox
# distutils: extra_compile_args = LINBOX_CFLAGS
# distutils: include_dirs = LINBOX_INCDIR
# distutils: libraries = LINBOX_LIBRARIES
# distutils: library_dirs = LINBOX_LIBDIR
# distutils: extra_link_args = LINBOX_LIBEXTRA
# distutils: language = c++

from libc.stdint cimport uint32_t, uint64_t
from libcpp.vector cimport vector as cppvector

from sage.libs.linbox.givaro cimport *

cdef extern from "linbox/matrix/dense-matrix.h":
    ## template <class _Field, class _blasRep=typename Vector<_Field>::Dense >
    ## class DenseMatrix ;
    ##
    ## template <class _Field>
    ## using DenseMatrix = DenseMatrix<_Field> ;
    cdef cppclass DenseMatrix_integer "LinBox::DenseMatrix<Givaro::ZRing<Givaro::Integer>>":
        ctypedef ZRing Field
        ctypedef Integer Element
        DenseMatrix_integer(Field &F, size_t m, size_t n)
        size_t rowdim()
        size_t coldim()
        void setEntry(size_t i, size_t j, Element &a)
        Element &getEntry(size_t i, size_t j)
        Field& field()

        ostream& write(ostream&)

    cdef cppclass DenseMatrix_Modular_double "LinBox::DenseMatrix<Givaro::Modular<double>>":
        ctypedef Modular_double Field
        ctypedef double Element
        DenseMatrix_Modular_double(Field F, size_t m, size_t n)
        void setEntry(size_t i, size_t j, Element& a)
        Element &getEntry(size_t i, size_t j)

        ostream& write(ostream&)

    cdef cppclass DenseMatrix_Modular_float "LinBox::DenseMatrix<Givaro::Modular<float>>":
        ctypedef Modular_float Field
        ctypedef float Element
        DenseMatrix_Modular_float(Field F, size_t m, size_t n)
        void setEntry(size_t i, size_t j, Element& a)
        Element &getEntry(size_t i, size_t j)

        ostream& write(ostream&)

cdef extern from "linbox/matrix/sparse-matrix.h":
    ## template<class _Field, class _Storage = SparseMatrixFormat::SparseSeq >
    ## class SparseMatrix ;
    cdef cppclass SparseMatrix_integer "LinBox::SparseMatrix<Givaro::ZRing<Givaro::Integer>>":
        ctypedef ZRing Field
        ctypedef Integer Element
        SparseMatrix_integer(Field &F, size_t m, size_t n)
        size_t rowdim()
        size_t coldim()
        void setEntry(size_t i, size_t j, Element &a)
        Element &getEntry(size_t i, size_t j)
        Field& field()

        ostream& write(ostream&)

    cdef cppclass SparseMatrix_Modular_uint64 "LinBox::SparseMatrix<Givaro::Modular<uint64_t>, LinBox::SparseMatrixFormat::SparseSeq>":
        ctypedef Modular_uint64 Field
        ctypedef uint64_t Element
        SparseMatrix_Modular_uint64(Field &F, size_t m, size_t n)
        SparseMatrix_Modular_uint64(SparseMatrix_Modular_uint64&)
        size_t rowdim()
        size_t coldim()
        void setEntry(size_t i, size_t j, Element &a)
        Element &getEntry(size_t i, size_t j)
        Field& field()

        ostream& write(ostream&)

cdef extern from "linbox/polynomial/dense-polynomial.h":
    ## template<class Field>
    ## class DensePolynomial : public Givaro::Poly1FactorDom<Field, Givaro::Dense>::Element
    cdef cppclass DensePolynomial_integer "LinBox::DensePolynomial<Givaro::ZRing<Givaro::Integer>>":
        ctypedef ZRing BaseRing
        ctypedef Integer BaseRingElement
        DensePolynomial_integer(BaseRing &F)
        DensePolynomial_integer(BaseRing &F, size_t s)
        BaseRingElement& operator[](size_t i)
        size_t size()

cdef extern from "linbox/ring/polynomial-ring.h":
    ## template <class BaseRing, class StorageTag= Givaro::Dense>
    ## class PolynomialRing : public Givaro::Poly1FactorDom<BaseRing,StorageTag>
    cdef cppclass PolynomialRing_integer "LinBox::PolynomialRing<Givaro::ZRing<Givaro::Integer>, Givaro::Dense>":
        ctypedef DensePolynomial_integer Element
        ctypedef DensePolynomial_integer Polynomial

cdef extern from "linbox/vector/vector.h":
    cdef cppclass DenseVector_integer "LinBox::DenseVector<Givaro::ZRing<Givaro::Integer>>":
        ctypedef ZRing Field
        ctypedef Integer Element
        DenseVector_integer (Field &F)
        DenseVector_integer (Field &F, long& m)
        DenseVector_integer (Field &F, cppvector[Integer]&)
        size_t size()
        void resize(size_t)
        void resize(size_t n, const Element&)
        void push_back(Element&)
        void clear()
        void reserve(const size_t&)
        void setEntry(size_t i, Element&)
        Element& getEntry(size_t i)

cdef extern from "linbox/matrix/matrix-domain.h":
    ## template <class Field_ >
    ## class MatrixDomain : public MVProductDomain<Field_> {
    cdef cppclass MatrixDomain_integer "LinBox::MatrixDomain<Givaro::ZRing<Givaro::Integer>>":
        MatrixDomain_integer(ZRing&)
        DenseMatrix_integer& mul(DenseMatrix_integer ans,
                                      DenseMatrix_integer left,
                                      DenseMatrix_integer right)

cdef extern from "linbox/solutions/methods.h" namespace "LinBox":
    cdef struct HybridSpecifier:
        pass
    cdef struct BlackboxSpecifier:
        pass
    cdef struct EliminationSpecifier:
        pass
    cdef struct WiedemannTraits:
        pass
    cdef struct DenseEliminationTraits:
        pass
    cdef struct SparseEliminationTraits:
        pass

    cdef cppclass Method:
        ctypedef HybridSpecifier Hybrid
        ctypedef BlackboxSpecifier Blackbox
        ctypedef EliminationSpecifier Elimination
        ctypedef WiedemannTraits Wiedemann
        ctypedef DenseEliminationTraits DenseElimination
        ctypedef SparseEliminationTraits SparseElimination

cdef extern from "linbox/solutions/charpoly.h" namespace "LinBox":
    PolynomialRing_integer.Element& charpoly (PolynomialRing_integer.Element&, DenseMatrix_integer&)
    PolynomialRing_integer.Element& charpoly (PolynomialRing_integer.Element&, SparseMatrix_integer&)

cdef extern from "linbox/solutions/minpoly.h" namespace "LinBox":
    PolynomialRing_integer.Element& minpoly (PolynomialRing_integer.Element&, DenseMatrix_integer&)
    PolynomialRing_integer.Element& minpoly (PolynomialRing_integer.Element&, SparseMatrix_integer&)

cdef extern from "linbox/algorithms/gauss.h":
    cdef cppclass GaussDomain_Modular_uint64 "LinBox::GaussDomain<Givaro::Modular<uint64_t>>":
        ctypedef Modular_uint64 Field
        ctypedef uint64_t Element
        GaussDomain_Modular_uint64(Field &)
        unsigned long& InPlaceLinearPivoting(unsigned long &rank,
                                             Element& determinant,
                                             SparseMatrix_Modular_uint64 &A,
                                             unsigned long Ni,
                                             unsigned long Nj)

cdef extern from "linbox/solutions/echelon.h" namespace "LinBox":
    size_t rowEchelon (DenseMatrix_Modular_float&, const DenseMatrix_Modular_float&)
    size_t rowEchelonize (DenseMatrix_Modular_float&)
    size_t reducedRowEchelon (DenseMatrix_Modular_float&, const DenseMatrix_Modular_float&)
    size_t reducedRowEchelonize (DenseMatrix_Modular_float&)
    size_t colEchelon (DenseMatrix_Modular_float&, const DenseMatrix_Modular_float&)
    size_t colEchelonize (DenseMatrix_Modular_float&)
    size_t reducedColEchelon (DenseMatrix_Modular_float&, const DenseMatrix_Modular_float&)
    size_t reducedColEchelonize (DenseMatrix_Modular_float&)

    size_t rowEchelon (DenseMatrix_Modular_double&, const DenseMatrix_Modular_double&)
    size_t rowEchelonize (DenseMatrix_Modular_double&)
    size_t reducedRowEchelon (DenseMatrix_Modular_double&, const DenseMatrix_Modular_double&)
    size_t reducedRowEchelonize (DenseMatrix_Modular_double&)
    size_t colEchelon (DenseMatrix_Modular_double&, const DenseMatrix_Modular_double&)
    size_t colEchelonize (DenseMatrix_Modular_double&)
    size_t reducedColEchelon (DenseMatrix_Modular_double&, const DenseMatrix_Modular_double&)
    size_t reducedColEchelonize (DenseMatrix_Modular_double&)

cdef extern from "linbox/solutions/rank.h" namespace "LinBox":
    unsigned long & rank (unsigned long&, DenseMatrix_integer)
    unsigned long & rank (unsigned long&, SparseMatrix_integer)

cdef extern from "linbox/solutions/det.h" namespace "LinBox":
    Integer& det (Integer&, DenseMatrix_integer)
    Integer& det (Integer&, SparseMatrix_integer)

cdef extern from "linbox/solutions/solve.h" namespace "LinBox":
    # integer solve

    DenseVector_integer& solve (DenseVector_integer &,
                                Integer &,
                                SparseMatrix_integer &,
                                DenseVector_integer &,
                                Method.DenseElimination) except +

    DenseVector_integer& solve (DenseVector_integer &,
                                Integer &,
                                SparseMatrix_integer &,
                                DenseVector_integer &,
                                Method.SparseElimination) except +

    DenseVector_integer& solve (DenseVector_integer &,
                                Integer &,
                                SparseMatrix_integer &,
                                DenseVector_integer &,
                                Method.Blackbox) except +

    DenseVector_integer& solve (DenseVector_integer &,
                                Integer &,
                                SparseMatrix_integer &,
                                DenseVector_integer &,
                                Method.Wiedemann) except +

    DenseVector_integer& solve (DenseVector_integer &,
                                Integer &,
                                SparseMatrix_integer &,
                                DenseVector_integer &)
