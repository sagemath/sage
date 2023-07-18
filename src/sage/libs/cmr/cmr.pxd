# -*- python -*-
# distutils: libraries = cmr

# (progn (replace-regexp "/[*]\\(.\\|\n\\)*?[*]/" "") (replace-regexp "[;{}]" "") (replace-regexp "CMR_EXPORT *" "") (replace-regexp "bool" "bint"))

cdef extern from "cmr/env.h":
    int CMR_OKAY
    int CMR_ERROR_INPUT
    int CMR_ERROR_OUTPUT
    int CMR_ERROR_MEMORY
    int CMR_ERROR_INVALID
    int CMR_ERROR_OVERFLOW
    int CMR_ERROR_TIMEOUT

    ctypedef int CMR_ERROR

    ctypedef struct CMR

    int CMRcreateEnvironment(CMR** pcmr)
    int CMRfreeEnvironment(CMR** pcmr)

    char* CMRgetErrorMessage(CMR* cmr)
    void CMRclearErrorMessage(CMR* cmr)

cdef extern from "cmr/matrix.h":

    ctypedef struct CMR_SUBMAT:
        size_t numRows
        size_t* rows
        size_t numColumns
        size_t* columns

    CMR_ERROR CMRsubmatCreate(CMR* cmr, size_t numRows, size_t numColumns, CMR_SUBMAT** psubmatrix)
    CMR_ERROR CMRsubmatCreate1x1(CMR* cmr, size_t row, size_t column, CMR_SUBMAT** psubmatrix)
    CMR_ERROR CMRsubmatFree(CMR* cmr, CMR_SUBMAT** psubmatrix)

    ctypedef struct CMR_CHRMAT:
        size_t numRows
        size_t numColumns
        size_t numNonzeros
        size_t* rowSlice
        size_t* entryColumns
        char* entryValues

    CMR_ERROR CMRchrmatCreate(CMR* cmr, CMR_CHRMAT** presult, int numRows, int numColumns, int numNonzeros)
    CMR_ERROR CMRchrmatSortNonzeros(CMR* cmr, CMR_CHRMAT* matrix)
    # CMR_ERROR CMRchrmatPrintDense(CMR* cmr, CMR_CHRMAT* matrix, FILE* stream, char zeroChar, bint header)
    CMR_ERROR CMRchrmatFindEntry(CMR_CHRMAT* matrix, size_t row, size_t column, size_t* pentry)
    CMR_ERROR CMRchrmatZoomSubmat(CMR* cmr, CMR_CHRMAT* matrix, CMR_SUBMAT* submatrix, CMR_CHRMAT** presult)
    CMR_ERROR CMRchrmatFree(CMR* cmr, CMR_CHRMAT** pmatrix)

cdef extern from "cmr/k_modular.h":

    CMR_ERROR CMRtestUnimodularity(CMR* cmr, CMR_CHRMAT* matrix, bint* pisUnimodular)
    CMR_ERROR CMRtestStrongUnimodularity(CMR* cmr, CMR_CHRMAT* matrix, bint* pisStronglyUnimodular)
    CMR_ERROR CMRtestKmodularity(CMR* cmr, CMR_CHRMAT* matrix, bint* pisKmodular, size_t* pk)
    CMR_ERROR CMRtestStrongKmodularity(CMR* cmr, CMR_CHRMAT* matrix, bint* pisStronglyKmodular, size_t* pk)

cdef extern from "cmr/camion.h":

    ctypedef struct CMR_CAMION_STATISTICS:
        size_t totalCount
        double totalTime

    CMR_ERROR CMRstatsCamionInit(CMR_CAMION_STATISTICS* stats)
    # CMR_ERROR CMRstatsCamionPrint(FILE* stream, CMR_CAMION_STATISTICS* stats, const char* prefix)
    CMR_ERROR CMRtestCamionSigned(CMR* cmr, CMR_CHRMAT* matrix, bint* pisCamionSigned, CMR_SUBMAT** psubmatrix, CMR_CAMION_STATISTICS* stats, double timeLimit)
    CMR_ERROR CMRcomputeCamionSigned(CMR* cmr, CMR_CHRMAT* matrix, bint* pwasCamionSigned, CMR_SUBMAT** psubmatrix, CMR_CAMION_STATISTICS* stats, double timeLimit)

cdef extern from "cmr/matroid.h":

    ctypedef struct CMR_MINOR:
        size_t numPivots
        size_t* pivotRows
        size_t* pivotColumns
        CMR_SUBMAT* remainingSubmatrix

    CMR_ERROR CMRminorCreate(CMR* cmr, CMR_MINOR** pminor, size_t numPivots, CMR_SUBMAT* submatrix)
    CMR_ERROR CMRminorFree(CMR* cmr, CMR_MINOR** pminor)
