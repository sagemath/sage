# sage_setup: distribution = sagemath-meataxe

#*****************************************************************************
#       Copyright (C) 2015 Simon King <simon.king@uni-jena.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

#
# Import SOME features from meataxe.h
# (most types are not needed, but listed here
# in the comments, for completeness)
#
cdef extern from "meataxe.h":
    # general ctype emulations
    # ctypedef int size_t   # size_t should be a standard type!
    # ctypedef unsigned long Ulong
    # ctypedef unsigned short Ushort
    # ctypedef unsigned char Uchar
    ctypedef unsigned char FEL
    ctypedef FEL *PTR

    # global constants
    cdef extern int FfOrder             # Current field order
    cdef extern int FfChar              # Current characteristic
    cdef extern FEL FfGen               # Generator
    cdef extern int FfNoc               # Number of columns for row ops
    cdef extern int MPB                 # Number of marks stored in a single byte
    cdef extern size_t FfCurrentRowSize # The byte size of a single row in memory,
                                        # always a multiple of sizeof(long)
    cdef extern size_t FfCurrentRowSizeIo # The number of bytes actually used in a row
    cdef extern char MtxLibDir[1024]    # Where to search/create multiplication tables

    # we only wrap MeatAxe for small fields (size < 255)
    cdef extern FEL mtx_tmult[256][256]
    cdef extern FEL mtx_tadd[256][256]
    cdef extern FEL mtx_taddinv[256]
    cdef extern FEL mtx_tmultinv[256]
    cdef extern FEL mtx_tinsert[8][256]
    cdef extern FEL mtx_textract[8][256]
    cdef extern FEL FF_ONE, FF_ZERO

#########################################
# function prototypes
    ## global parameters
    size_t FfRowSize(int noc)
    size_t FfTrueRowSize(int noc) # Difference to FfRowSize: Doesn't count padding bytes
    int FfSetField(int field) except -1
    int FfSetNoc(int ncols) except -1

    ## Finite Fields
    FEL FfAdd(FEL a, FEL b)
    FEL FfSub(FEL a, FEL b)
    FEL FfNeg(FEL a)
    FEL FfMul(FEL a, FEL b)
    FEL FfDiv(FEL a, FEL b)
    FEL FfInv(FEL a)
    # FEL FfEmbed(FEL a, int subfield) except 255
    # FEL FfRestrict(FEL a, int subfield) except 255
    FEL FfFromInt(int l)
    int FfToInt(FEL f)

    ## Rows
    void FfMulRow(PTR row, FEL mark)
    void FfAddMulRow(PTR dest, PTR src, FEL f)
    void FfAddMulRowPartial(PTR dest, PTR src, FEL f, int first, int len)
    PTR FfAddRow(PTR dest, PTR src)
    PTR FfAddRowPartial(PTR dest, PTR src, int first, int len)
    PTR FfSubRow(PTR dest, PTR src)
    PTR FfSubRowPartial(PTR dest, PTR src, int first, int len)
    FEL FfExtract(PTR row, int col)
    void FfInsert(PTR row, int col, FEL mark)
    int FfFindPivot(PTR row, FEL *mark)
    FEL FfScalarProduct(PTR a, PTR b)
    void FfSwapRows(PTR dest, PTR src)
    void FfPermRow(PTR row, long *perm, PTR result)
    int FfCmpRows(PTR p1, PTR p2)

    ## multiple rows
    PTR FfAlloc(int nor) except NULL
    void FfExtractColumn(PTR mat,int nor,int col,PTR result)
    int FfStepPtr(PTR *x)  # Advance to next row
    PTR FfGetPtr(PTR base, int row)  # Advance to "row" rows after base
    void FfInsert(PTR row, int col, FEL mark)
    void FfMapRow(PTR row, PTR matrix, int nor, PTR result)

    ############
    ## Skip: Application, error handling, i/o

    ############
    ## Matrices
    ############
    ctypedef struct Matrix_t:
        int Field, Nor, Noc        # Field, #rows, #columns
        PTR Data                   # Pointer to data area
        int RowSize                # Size (in bytes) of one row
        int *PivotTable            # Pivot table (if matrix is in echelon form)

    ## Basic memory operations
    Matrix_t *MatAlloc(int field, int nor, int noc) except NULL
    int MatFree(Matrix_t *mat)
    PTR MatGetPtr(Matrix_t *mat, int row) except NULL
    int MatCompare(Matrix_t *a, Matrix_t *b) except -2
    int MatCopyRegion(Matrix_t *dest, int destrow, int destcol, Matrix_t *src, int row1, int col1, int nrows, int ncols) except -1
    Matrix_t *MatCut(Matrix_t *src, int row1, int col1, int nrows, int ncols) except NULL
    Matrix_t *MatCutRows(Matrix_t *src, int row1, int nrows) except NULL
    Matrix_t *MatDup(Matrix_t *src) except NULL
    Matrix_t *MatId(int fl, int nor) except NULL
    Matrix_t *MatLoad(char *fn) except? NULL
    int MatSave(Matrix_t *mat, char *fn) except -1

    ## Basic Arithmetic  ## general rule: dest is changed, src/mat are unchanged!
    Matrix_t *MatTransposed(Matrix_t *src) except NULL
    Matrix_t *MatAdd(Matrix_t *dest, Matrix_t *src) except NULL
    Matrix_t *MatAddMul(Matrix_t *dest, Matrix_t *src, FEL coeff) except NULL
    Matrix_t *MatMul(Matrix_t *dest, Matrix_t *src) except NULL
    Matrix_t *MatMulScalar(Matrix_t *dest, FEL coeff) except NULL
    Matrix_t *MatPower(Matrix_t *mat, long n) except NULL
    int StablePower(Matrix_t *mat, int *pwr, Matrix_t **ker) except -1
    FEL MatTrace(Matrix_t *mat) except 255
    Matrix_t *MatMulStrassen(Matrix_t *dest, Matrix_t *A, Matrix_t *B) except NULL
    void StrassenSetCutoff(size_t size)

    ## "Higher" Arithmetic
    Matrix_t *MatTensor(Matrix_t *m1, Matrix_t *m2) except NULL
    Matrix_t *TensorMap(Matrix_t *vec, Matrix_t *a, Matrix_t *b) except NULL

    int MatClean(Matrix_t *mat, Matrix_t *sub) except -1
    int MatEchelonize(Matrix_t *mat) except -1
    int MatOrder(Matrix_t *mat) except? -1
    long MatNullity(Matrix_t *mat)
    Matrix_t *MatInverse(Matrix_t *src) except NULL
    Matrix_t *MatNullSpace(Matrix_t *mat) except NULL
    Matrix_t *MatNullSpace__(Matrix_t *mat) except NULL

    ## Error handling
    cdef extern int MTX_ERR_NOMEM, MTX_ERR_GAME_OVER, MTX_ERR_DIV0, MTX_ERR_FILEFMT, MTX_ERR_BADARG
    cdef extern int MTX_ERR_RANGE, MTX_ERR_NOTECH, MTX_ERR_NOTSQUARE, MTX_ERR_INCOMPAT
    cdef extern int MTX_ERR_BADUSAGE, MTX_ERR_OPTION, MTX_ERR_NARGS, MTX_ERR_NOTMATRIX, MTX_ERR_NOTPERM
    ctypedef struct MtxFileInfo_t:
        char *Name
        char *BaseName

    ctypedef struct MtxErrorRecord_t:
        MtxFileInfo_t *FileInfo
        int LineNo
        char *Text

    ctypedef void MtxErrorHandler_t(MtxErrorRecord_t*)
    MtxErrorHandler_t *MtxSetErrorHandler(MtxErrorHandler_t *h)

cdef Matrix_t *rawMatrix(int Field, list entries) except NULL
