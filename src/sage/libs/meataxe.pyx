# sage_setup: distribution = sagemath-meataxe
# distutils: libraries = mtx

# ****************************************************************************
#       Copyright (C) 2017 Simon King <simon.king@uni-jena.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cpython.exc cimport PyErr_SetObject
from cysignals.signals cimport sig_block, sig_unblock
cdef dict ErrMsg = {
    "Not enough memory": MemoryError,
    "Time limit exceeded": RuntimeError,
    "Division by zero": ZeroDivisionError,
    "Bad file format": OSError,
    "No such file or directory": OSError,
    "Bad argument": ValueError,
    "Argument out of range": IndexError,

    "Matrix not in echelon form": ValueError,
    "Matrix not square": ArithmeticError,
    "Incompatible objects": TypeError,

    "Bad syntax, try `-help'": SyntaxError,
    "Bad usage of option, try `-help'": ValueError,
    "Bad number of arguments, try `-help'": ValueError,

    "Not a matrix": TypeError,
    "Not a permutation": TypeError
}


cdef Matrix_t *rawMatrix(int Field, list entries) except NULL:
    """
    Return a meataxe matrix.

    INPUT:

    - ``Field`` -- integer; the field size
    - ``entries`` -- list of lists, the entries of the matrix, also
      defining the matrix dimensions. It is *not* tested that all rows
      in ``entries`` have the same length, and it is assumed that both
      the number of rows and the number of columns is positive. All
      elements are given by ints, they are automatically interpreted as
      field elements.
    """
    cdef Matrix_t *M = MatAlloc(Field, len(entries), len(entries[0]))
    cdef PTR x = M.Data
    cdef int idx, i, j
    cdef list dt_i
    for i in range(M.Nor):
        idx = 0
        dt_i = entries[i]
        for j in range(M.Noc):
            FfInsert(x, j, FfFromInt(dt_i[j]))
        FfStepPtr(&(x))
    return M

###############################################################
## It is needed to do some initialisation. Since the meataxe
## version used in Sage (SharedMeatAxe) is a dynamic (shared)
## library, it sufficed to do this initialisation only once.
## For convenience, the meataxe_init() function is called in
## this module. Hence, `import sage.libs.meataxe` is enough
## to make sure that MeatAxe is initialised.

from sage.cpython.string cimport str_to_bytes, char_to_str

cdef void sage_meataxe_error_handler(const MtxErrorRecord_t *err) noexcept:
    sig_block()
    ErrText  = char_to_str(err.Text)
    BaseName = char_to_str(err.FileInfo.BaseName)
    LineNo   = err.LineNo
    PyErr_SetObject(ErrMsg.get(ErrText.split(': ')[-1], RuntimeError), f"{ErrText} in file {BaseName} (line {LineNo})")
    sig_unblock()

cdef inline meataxe_init() noexcept:
    ## Assign to a variable that enables MeatAxe to find
    ## its multiplication tables.
    global MtxLibDir
    from sage import env
    mtxdir = str_to_bytes(env.MTXLIB)
    if len(mtxdir) >= 1024:
        raise RuntimeError(f"the path for the meataxe library {mtxdir!r} is too long, it needs to be of length < 1024")
    MtxLibDir[:len(mtxdir)] = mtxdir
    MtxLibDir[len(mtxdir)] = c'\0'
    ## Error handling for MeatAxe, to prevent immediate exit of the program
    MtxSetErrorHandler(sage_meataxe_error_handler)

meataxe_init()
