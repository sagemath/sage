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

from libc.stdint cimport intptr_t, uintptr_t, uint8_t, uint16_t, uint32_t, uint64_t

cdef extern from "gap/system.h" nogil:
    ctypedef char Char
    ctypedef intptr_t Int
    ctypedef uintptr_t UInt
    ctypedef uint8_t  UInt1
    ctypedef uint16_t UInt2
    ctypedef uint32_t UInt4
    ctypedef uint64_t UInt8
    ctypedef void* Obj


cdef extern from "gap/calls.h" nogil:
    bint IS_FUNC(Obj)


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

    void GAP_MarkBag(Obj bag)
    void GAP_CollectBags(UInt full)

    Obj GAP_SUM(Obj, Obj)
    Obj GAP_DIFF(Obj, Obj)
    Obj GAP_PROD(Obj, Obj)
    Obj GAP_QUO(Obj, Obj)
    Obj GAP_POW(Obj, Obj)
    Obj GAP_MOD(Obj, Obj)
    bint GAP_EQ(Obj opL, Obj opR)
    bint GAP_LT(Obj opL, Obj opR)
    bint GAP_IN(Obj opL, Obj opR)

    cdef Obj GAP_True
    cdef Obj GAP_False

    Obj GAP_CallFuncList(Obj func, Obj args);
    Obj GAP_CallFuncArray(Obj func, UInt narg, Obj * args);
    Obj GAP_CallFunc0Args(Obj func);
    Obj GAP_CallFunc1Args(Obj func, Obj a1);
    Obj GAP_CallFunc2Args(Obj func, Obj a1, Obj a2);
    Obj GAP_CallFunc3Args(Obj func, Obj a1, Obj a2, Obj a3);

    bint GAP_IsMacFloat(Obj obj)
    double GAP_ValueMacFloat(Obj obj)

    bint GAP_IsInt(Obj)
    bint GAP_IsSmallInt(Obj)
    Obj GAP_NewObjIntFromInt(Int val)
    Int GAP_ValueInt(Obj)
    Int GAP_SizeInt(Obj)
    UInt* GAP_AddrInt(Obj)

    bint GAP_IsList(Obj lst)
    UInt GAP_LenList(Obj lst)
    void GAP_AssList(Obj lst, UInt pos, Obj val)
    Obj GAP_ElmList(Obj lst, UInt pos)
    Obj GAP_NewPlist(Int capacity)

    bint GAP_IsRecord(Obj obj)
    Obj GAP_NewPrecord(Int capacity)

    bint GAP_IsString(Obj obj)
    UInt GAP_LenString(Obj string)
    char* GAP_CSTR_STRING(Obj list)
    Obj GAP_MakeStringWithLen(const char* buf, UInt len)

    Int GAP_ValueOfChar(Obj obj)


cdef extern from "gap/lists.h" nogil:
    Obj ELM_LIST(Obj lst, int pos)


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
        T_PERM2
        T_PERM4
        T_BOOL
        T_FUNCTION
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
    int LEN_PREC(Obj rec)
    int GET_RNAM_PREC(Obj rec, int i)
    Obj GET_ELM_PREC(Obj rec, int i)
    void AssPRec(Obj rec, UInt rnam, Obj val)


cdef extern from "gap/records.h" nogil:
    char* NAME_RNAM(UInt rnam)
    Obj ELM_REC(Obj rec, UInt rnam)
    UInt RNamName(Char* name)


cdef extern from "gap/stringobj.h" nogil:
    bint IS_STRING(Obj obj)
    bint IsStringConv(Obj obj)
    Obj NEW_STRING(Int)


cdef extern from "<structmember.h>" nogil:
    """
    /* Hack: Cython 3.0 automatically includes <structmember.h>, which
     * defines several macros that collides with enum definitions in
     * gap/objects.h. We need to include the header explicitly and
     * undefine these macros.
     */
    #undef T_INT
    #undef T_STRING
    #undef T_CHAR
    #undef T_BOOL
    """
    pass
