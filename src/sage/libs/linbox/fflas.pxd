from .modular cimport ModDoubleField, ModFloatField, ModDoubleFieldElement, ModFloatFieldElement

cdef extern from "fflas-ffpack/fflas-ffpack.h" namespace "std":
    cdef cppclass vector[T]:
        cppclass iterator:
            T operator*()
            iterator operator++()
            bint operator==(iterator)
            bint operator!=(iterator)
        vector()
        void push_back(T&)
        T& operator[](int)
        T& at(int)
        iterator begin()
        iterator end()
        size_t size()

    cdef cppclass list[T]:
        cppclass iterator:
            T operator*()
            iterator operator++()
            bint operator==(iterator)
            bint operator!=(iterator)
        void push_back(T&)
        void pop_front()
        T& front()
        iterator begin()
        iterator end()
        size_t size()
        void clear()

cdef extern from "fflas-ffpack/fflas-ffpack.h":
    ctypedef enum fflas_trans_enum "FFLAS::FFLAS_TRANSPOSE":
        fflas_no_trans  "FFLAS::FflasNoTrans"
        fflas_trans  "FFLAS::FflasTrans"

    ctypedef enum fflas_side_enum "FFLAS::FFLAS_SIDE":
        fflas_right  "FFLAS::FflasRight"

    # double
    void ModDouble_fgemv "FFLAS::fgemv" \
            (ModDoubleField F, fflas_trans_enum transA,
             size_t nrows, size_t ncols,
             ModDoubleFieldElement alpha, ModDoubleFieldElement* A,
             size_t lda, ModDoubleFieldElement* X, size_t incX,
             ModDoubleFieldElement beta, ModDoubleFieldElement* Y,
             size_t incY)

    ModDoubleFieldElement* ModDouble_fgemm "FFLAS::fgemm" \
            (ModDoubleField F,
             fflas_trans_enum transA, fflas_trans_enum transB,
             size_t nrowsA, size_t ncolsB, size_t ncolsA,
             ModDoubleFieldElement alpha, ModDoubleFieldElement* A,
             size_t A_stride, ModDoubleFieldElement* B, int B_stride,
             ModDoubleFieldElement beta, ModDoubleFieldElement* C,
             size_t C_stride)


    # float
    void ModFloat_fgemv "FFLAS::fgemv" \
            (ModFloatField F, fflas_trans_enum transA,
             size_t nrows, size_t ncols,
             ModFloatFieldElement alpha, ModFloatFieldElement* A,
             size_t lda, ModFloatFieldElement* X, size_t incX,
             ModFloatFieldElement beta, ModFloatFieldElement* Y,
             size_t incY)

    ModFloatFieldElement* ModFloat_fgemm "FFLAS::fgemm" \
            (ModFloatField F,
             fflas_trans_enum transA, fflas_trans_enum transB,
             size_t nrowsA, size_t ncolsB, size_t ncolsA,
             ModFloatFieldElement alpha, ModFloatFieldElement* A,
             size_t A_stride, ModFloatFieldElement* B, int B_stride,
             ModFloatFieldElement beta, ModFloatFieldElement* C,
             size_t C_stride)

cdef extern from "fflas-ffpack/fflas-ffpack.h":
    # double
    bint ModDouble_is_singular "FFPACK::IsSingular" (ModDoubleField F,
                                                     size_t nrows, size_t ncols, ModDoubleFieldElement* A,
                                                     size_t A_stride)

    ModDoubleFieldElement* ModDouble_invert_in_place "FFPACK::Invert"  (ModDoubleField F, size_t order,
                                                                        ModDoubleFieldElement* A, size_t A_stride, int nullity)

    ModDoubleFieldElement ModDoubleDet "FFPACK::Det" (ModDoubleField F,
                                                      size_t nrows, size_t ncols,
                                                      ModDoubleFieldElement* A, size_t A_stride)

    int ModDoubleRank "FFPACK::Rank" (ModDoubleField,
                                      size_t nrows, size_t ncols,
                                      ModDoubleFieldElement *A, size_t lda)

    size_t ModDouble_echelon "FFPACK::ReducedRowEchelonForm" (ModDoubleField F, size_t a, size_t b,
                                                              ModDoubleFieldElement* matrix,
                                                              size_t s, size_t* P, size_t* Q)

    void ModDouble_applyp "FFPACK::applyP" (ModDoubleField F,
                                            fflas_side_enum s, fflas_trans_enum tr,
                                            size_t nr, size_t foo, size_t r,
                                            ModDoubleFieldElement* matrix, size_t nc, size_t* Q)

    void ModDouble_MinPoly "FFPACK::MinPoly" ( ModDoubleField F,
                                               vector[ModDoubleFieldElement] minP, size_t N,
                                               ModDoubleFieldElement* A, size_t lda,
                                               ModDoubleFieldElement* X, size_t ldx, size_t* P)

    void ModDouble_CharPoly "FFPACK::CharPoly" ( ModDoubleField F,
                                                 list[vector[ModDoubleFieldElement]] charp, size_t N,
                                                 ModDoubleFieldElement* A, size_t lda)

    # float

    bint ModFloat_is_singular "FFPACK::IsSingular" (ModFloatField F,
                                                    size_t nrows, size_t ncols, ModFloatFieldElement* A,
                                                    size_t A_stride)

    ModFloatFieldElement* ModFloat_invert_in_place "FFPACK::Invert" (ModFloatField F, size_t order,
                                                                     ModFloatFieldElement* A, size_t A_stride, int nullity)

    ModFloatFieldElement ModFloatDet "FFPACK::Det" (ModFloatField F,
                                                    size_t nrows, size_t ncols,
                                                    ModFloatFieldElement* A, size_t A_stride)

    int ModFloatRank "FFPACK::Rank" (ModFloatField,
                                     size_t nrows, size_t ncols,
                                     ModFloatFieldElement *A, size_t lda)

    size_t ModFloat_echelon "FFPACK::ReducedRowEchelonForm" (ModFloatField F, size_t a, size_t b,
                                                             ModFloatFieldElement* matrix,
                                                             size_t s, size_t* P, size_t* Q)

    void ModFloat_applyp "FFPACK::applyP" (ModFloatField F,
                                           fflas_side_enum s, fflas_trans_enum tr,
                                           size_t nr, size_t foo, size_t r,
                                           ModFloatFieldElement* matrix, size_t nc, size_t* Q)

    void ModFloat_MinPoly "FFPACK::MinPoly" ( ModFloatField F,
                                              vector[ModFloatFieldElement] minP, size_t N,
                                              ModFloatFieldElement* A, size_t lda,
                                              ModFloatFieldElement* X, size_t ldx, size_t* P)

    void ModFloat_CharPoly "FFPACK::CharPoly" ( ModFloatField F,
                                                list[vector[ModFloatFieldElement]] charp, size_t N,
                                                ModFloatFieldElement* A, size_t lda )
