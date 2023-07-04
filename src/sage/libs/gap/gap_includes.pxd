# distutils: libraries = gap gmp m
###############################################################################
#       Copyright (C) 2009, William Stein <wstein@gmail.com>
#       Copyright (C) 2012, Volker Braun <vbraun.name@gmail.com>
#
#   Distributed under the terms of the GNU General Public License (GPL)
#   as published by the Free Software Foundation; either version 2 of
#   the License, or (at your option) any later version.
#                   http://www.gnu.org/licenses/
###############################################################################

from libc.stdint cimport uintptr_t, uint8_t, uint16_t, uint32_t, uint64_t

cdef extern from "gap/system.h" nogil:
    ctypedef char Char
    ctypedef int Int
    ctypedef uintptr_t UInt
    ctypedef uint8_t  UInt1
    ctypedef uint16_t UInt2
    ctypedef uint32_t UInt4
    ctypedef uint64_t UInt8
    ctypedef void* Obj


cdef extern from "gap/ariths.h" nogil:
    Obj SUM(Obj, Obj)
    Obj DIFF(Obj, Obj)
    Obj PROD(Obj, Obj)
    Obj QUO(Obj, Obj)
    Obj POW(Obj, Obj)
    Obj MOD(Obj, Obj)
    bint EQ(Obj opL, Obj opR)
    bint LT(Obj opL, Obj opR)


cdef extern from "gap/bool.h" nogil:
    cdef Obj GAP_True "True"
    cdef Obj GAP_False "False"


cdef extern from "gap/calls.h" nogil:
    bint IS_FUNC(Obj)
    Obj CALL_0ARGS(Obj f)              # 0 arguments
    Obj CALL_1ARGS(Obj f, Obj a1)      # 1 argument
    Obj CALL_2ARGS(Obj f, Obj a1, Obj a2)
    Obj CALL_3ARGS(Obj f, Obj a1, Obj a2, Obj a3)
    Obj CALL_4ARGS(Obj f, Obj a1, Obj a2, Obj a3, Obj a4)
    Obj CALL_5ARGS(Obj f, Obj a1, Obj a2, Obj a3, Obj a4, Obj a5)
    Obj CALL_6ARGS(Obj f, Obj a1, Obj a2, Obj a3, Obj a4, Obj a5, Obj a6)
    Obj CALL_XARGS(Obj f, Obj args)   # more than 6 arguments


cdef extern from "gap/gasman.h" nogil:
    void MarkBag(Obj bag)
    UInt CollectBags(UInt size, UInt full)


cdef extern from "gap/integer.h" nogil:
    Int IS_INT(Obj)


cdef extern from "gap/intobj.h" nogil:
    bint IS_INTOBJ(Obj obj)
    Obj INTOBJ_INT(Int)
    Int INT_INTOBJ(Obj)


cdef extern from "gap/io.h" nogil:
    ctypedef struct TypOutputFile:
        pass
    UInt OpenOutputStream(TypOutputFile* output, Obj stream)
    UInt CloseOutput(TypOutputFile* output)


cdef extern from "gap/libgap-api.h" nogil:
    """
    #define sig_GAP_Enter()  {int t = GAP_Enter(); if (!t) sig_error();}
    """
    ctypedef void (*GAP_CallbackFunc)()
    void GAP_Initialize(int argc, char ** argv,
            GAP_CallbackFunc markBagsCallback, GAP_CallbackFunc errorCallback,
            int handleSignals)
    Obj GAP_EvalString(const char *) except *
    Obj GAP_EvalStringNoExcept "GAP_EvalString"(const char *)
    Obj GAP_ValueGlobalVariable(const char *)
    cdef void GAP_EnterStack()
    cdef void GAP_LeaveStack()
    cdef int GAP_Enter() except 0
    cdef void sig_GAP_Enter()
    cdef void GAP_Leave()
    cdef int GAP_Error_Setjmp() except 0


cdef extern from "gap/lists.h" nogil:
    bint IS_LIST(Obj lst)
    int LEN_LIST(Obj lst)
    Obj ELM_LIST(Obj lst, int pos)
    Obj ELM0_LIST(Obj lst, int pos)
    void ASS_LIST(Obj lst, int pos, Obj elt)


cdef extern from "gap/listfunc.h" nogil:
    void AddList(Obj list, Obj obj)


cdef extern from "gap/macfloat.h" nogil:
    double VAL_MACFLOAT(Obj obj)


cdef extern from "gap/objects.h" nogil:
    bint IS_MUTABLE_OBJ(Obj obj)
    Obj SHALLOW_COPY_OBJ(Obj obj)
    Obj CopyObj(Obj obj, int mut)

    UInt TNUM_OBJ(Obj obj)
    char* TNAM_OBJ(Obj obj)

    cdef enum TNUM:
        T_RAT
        T_CYC
        T_FFE
        T_MACFLOAT
        T_PERM2
        T_PERM4
        T_BOOL
        T_CHAR
        T_FUNCTION
        T_PLIST
        T_COMOBJ
        T_POSOBJ


cdef extern from "gap/permutat.h" nogil:
    UInt DEG_PERM2(Obj)
    UInt DEG_PERM4(Obj)
    Obj NEW_PERM2(UInt)
    Obj NEW_PERM4(UInt)
    UInt2* ADDR_PERM2(Obj)
    UInt4* ADDR_PERM4(Obj)
    const UInt2* CONST_ADDR_PERM2(Obj)
    const UInt4* CONST_ADDR_PERM4(Obj)


cdef extern from "gap/precord.h" nogil:
    Obj NEW_PREC(int len)
    int LEN_PREC(Obj rec)
    int GET_RNAM_PREC(Obj rec, int i)
    Obj GET_ELM_PREC(Obj rec, int i)
    void AssPRec(Obj rec, UInt rnam, Obj val)


cdef extern from "gap/records.h" nogil:
    char* NAME_RNAM(UInt rnam)
    bint IS_REC(Obj obj)
    Obj ELM_REC(Obj rec, UInt rnam)
    UInt RNamName(Char* name)


cdef extern from "gap/stringobj.h" nogil:
    char* CSTR_STRING(Obj list)
    bint IS_STRING(Obj obj)
    bint IsStringConv(Obj obj)
    Obj NEW_STRING(Int)
    Obj MakeStringWithLen(const char* buf, size_t len)
