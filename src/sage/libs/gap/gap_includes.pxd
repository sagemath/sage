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


cdef extern from "<gap/system.h>":
    ctypedef char Char
    ctypedef int Int
    ctypedef unsigned int UInt
    ctypedef void* Obj


cdef extern from "<gap/ariths.h>":
    Obj SUM(Obj, Obj)
    Obj DIFF(Obj, Obj)
    Obj PROD(Obj, Obj)
    Obj QUO(Obj, Obj)
    Obj POW(Obj, Obj)
    Obj MOD(Obj, Obj)
    bint EQ(Obj opL, Obj opR)
    bint LT(Obj opL, Obj opR)


cdef extern from "<gap/bool.h>":
    cdef Obj GAP_True "True"
    cdef Obj GAP_False "False"


cdef extern from "<gap/calls.h>":
    bint IS_FUNC(Obj)
    Obj CALL_0ARGS(Obj f)              # 0 arguments
    Obj CALL_1ARGS(Obj f, Obj a1)      # 1 argument
    Obj CALL_2ARGS(Obj f, Obj a1, Obj a2)
    Obj CALL_3ARGS(Obj f, Obj a1, Obj a2, Obj a3)
    Obj CALL_4ARGS(Obj f, Obj a1, Obj a2, Obj a3, Obj a4)
    Obj CALL_5ARGS(Obj f, Obj a1, Obj a2, Obj a3, Obj a4, Obj a5)
    Obj CALL_6ARGS(Obj f, Obj a1, Obj a2, Obj a3, Obj a4, Obj a5, Obj a6)
    Obj CALL_XARGS(Obj f, Obj args)   # more than 6 arguments


cdef extern from "<gap/gasman.h>":
    Obj NewBag "NewBag"(UInt type, UInt size)
    void MarkBag(Obj bag)
    UInt CollectBags(UInt size, UInt full)


cdef extern from "<gap/gasman_intern.h>":
    void CallbackForAllBags(void (*func)(Obj))


cdef extern from "<gap/gvars.h>":
    UInt GVarName "GVarName"(char* name)
    void AssGVar "AssGVar"(UInt gvar, Obj val)


cdef extern from "<gap/integer.h>":
    Int IS_INT(Obj)


cdef extern from "<gap/intobj.h>":
    bint IS_INTOBJ(Obj obj)
    Obj INTOBJ_INT(Int)
    Int INT_INTOBJ(Obj)


cdef extern from "<gap/io.h>":
    UInt OpenOutputStream(Obj stream)
    UInt CloseOutput()


cdef extern from "<gap/libgap-api.h>":
    ctypedef void (*CallbackFunc)()
    void GAP_Initialize(int argc, char ** argv, char ** env,
        CallbackFunc, CallbackFunc)
    Obj GAP_EvalString(const char *) except *
    Obj GAP_EvalStringNoExcept "GAP_EvalString"(const char *)
    Obj GAP_ValueGlobalVariable(const char *)


cdef extern from "<gap/libgap-api.h>" nogil:
    """
    #define sig_GAP_Enter()  {int t = GAP_Enter(); if (!t) sig_error();}
    """
    cdef void GAP_EnterStack()
    cdef void GAP_LeaveStack()
    cdef int GAP_Enter() except 0
    cdef void sig_GAP_Enter()
    cdef void GAP_Leave()
    cdef int GAP_Error_Setjmp() except 0


cdef extern from "<gap/lists.h>":
    bint IS_LIST(Obj lst)
    int LEN_LIST(Obj lst)
    Obj ELM_LIST(Obj lst, int pos)
    Obj ELM0_LIST(Obj lst, int pos)
    void ASS_LIST(Obj lst, int pos, Obj elt)


cdef extern from "<gap/listfunc.h>":
    void AddList(Obj list, Obj obj)


cdef extern from "<gap/macfloat.h>":
    double VAL_MACFLOAT(Obj obj)


cdef extern from "<gap/objects.h>":
    bint IS_MUTABLE_OBJ(Obj obj)
    Obj SHALLOW_COPY_OBJ(Obj obj)
    Obj CopyObj(Obj obj, int mut)

    UInt SIZE_OBJ(Obj obj)
    UInt TNUM_OBJ(Obj obj)
    char* TNAM_OBJ(Obj obj)

    cdef int T_INTPOS
    cdef int T_INTNEG
    cdef int T_RAT
    cdef int T_CYC
    cdef int T_FFE
    cdef int T_PERM2
    cdef int T_PERM4
    cdef int T_BOOL
    cdef int T_CHAR
    cdef int T_FUNCTION
    cdef int T_MACFLOAT
    cdef int T_PLIST
    cdef int T_PLIST_CYC
    cdef int T_BLIST
    cdef int T_STRING
    cdef int T_COMOBJ
    cdef int T_POSOBJ
    cdef int T_DATOBJ
    cdef int T_WPOBJ


cdef extern from "<gap/precord.h>":
    Obj NEW_PREC(int len)
    int LEN_PREC(Obj rec)
    int GET_RNAM_PREC(Obj rec, int i)
    Obj GET_ELM_PREC(Obj rec, int i)
    void AssPRec(Obj rec, UInt rnam, Obj val)


cdef extern from "<gap/records.h>":
    char* NAME_RNAM(UInt rnam)
    bint IS_REC(Obj obj)
    Obj ELM_REC(Obj rec, UInt rnam)
    UInt RNamName(Char* name)


cdef extern from "<gap/stringobj.h>":
    char* CSTR_STRING(Obj list)
    bint IS_STRING(Obj obj)
    bint IsStringConv(Obj obj)
    Obj NEW_STRING(Int)
    void C_NEW_STRING(Obj new_gap_string, int length, char* c_string)
