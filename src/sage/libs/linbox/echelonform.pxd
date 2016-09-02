# distutils: extra_compile_args = LINBOX_CFLAGS
# distutils: libraries = LINBOX_LIBRARIES
# distutils: library_dirs = LINBOX_LIBDIR
# distutils: language = c++

from sage.libs.linbox.modular cimport ModDoubleField, ModDoubleFieldElement, ModFloatField, ModFloatFieldElement
from libc.stdint cimport uint64_t

cdef extern from "linbox/matrix/dense-matrix.h" namespace "LinBox":
    cdef cppclass BlasMatrixDouble "LinBox::DenseMatrix<Givaro::Modular<double> >":
        BlasMatrixDouble(ModDoubleField F, uint64_t nrows, uint64_t ncols)
        void setEntry(size_t i, size_t j, ModDoubleFieldElement t)
        ModDoubleFieldElement &getEntry(size_t i, size_t j)

    cdef cppclass BlasMatrixFloat "LinBox::DenseMatrix<Givaro::Modular<float> >":
        BlasMatrixFloat(ModFloatField F, uint64_t nrows, uint64_t ncols)
        void setEntry(size_t i, size_t j, ModFloatFieldElement t)
        ModFloatFieldElement &getEntry(size_t i, size_t j)

cdef extern from "linbox/algorithms/echelon-form.h":
    cdef cppclass EchelonFormDomainDouble "LinBox::EchelonFormDomain<Givaro::Modular<double> >":
        EchelonFormDomainDouble(ModDoubleField)
        int rowReducedEchelon(BlasMatrixDouble, BlasMatrixDouble)

    cdef cppclass EchelonFormDomainFloat "LinBox::EchelonFormDomain<Givaro::Modular<float> >":
        EchelonFormDomainFloat(ModFloatField)
        int rowReducedEchelon(BlasMatrixFloat, BlasMatrixFloat)

