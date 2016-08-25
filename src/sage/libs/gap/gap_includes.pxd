###############################################################################
#       Copyright (C) 2009, William Stein <wstein@gmail.com>
#       Copyright (C) 2012, Volker Braun <vbraun.name@gmail.com>
#
#   Distributed under the terms of the GNU General Public License (GPL)
#   as published by the Free Software Foundation; either version 2 of
#   the License, or (at your option) any later version.
#                   http://www.gnu.org/licenses/
###############################################################################


cdef extern from "gap/system.h":
    ctypedef char libGAP_Char
    ctypedef int libGAP_Int
    ctypedef unsigned char libGAP_UChar

cdef extern from "gap/libgap.h":
    void libgap_initialize(int argc, char** argv)
    ctypedef void(*libgap_gasman_callback_ptr)()
    void libgap_set_gasman_callback(libgap_gasman_callback_ptr callback)
    ctypedef void(*libgap_error_func_ptr)(char* msg)
    void libgap_set_error_handler(libgap_error_func_ptr error_handler)
    void libgap_call_error_handler()
    void libgap_finalize()
    void libgap_start_interaction(char* inputline)
    char* libgap_get_output()
    char* libgap_get_error()
    void libgap_finish_interaction()
    void libgap_mark_stack_bottom()
    void libgap_enter()
    void libgap_exit()

cdef extern from "gap/code.h":
    ctypedef unsigned int libGAP_Stat
    ctypedef libGAP_Stat* libGAP_PtrBody

cdef extern from "gap/gap.h":
    ctypedef unsigned int libGAP_UInt
    ctypedef void* libGAP_ExecStatus
    void libGAP_ViewObjHandler(void*)
    void libGAP_InitializeGap(int*, char** argv)
    void libGAP_set_system_variables(char**, char**)
    cdef libGAP_UInt libGAP_Last
    cdef libGAP_UInt libGAP_Last2
    cdef libGAP_UInt libGAP_Last3
    cdef libGAP_ExecStatus libGAP_STATUS_END
    cdef libGAP_ExecStatus libGAP_STATUS_RETURN_VAL
    cdef libGAP_ExecStatus libGAP_STATUS_RETURN_VOID
    cdef libGAP_ExecStatus libGAP_STATUS_TNM
    cdef libGAP_ExecStatus libGAP_STATUS_QUIT
    cdef libGAP_ExecStatus libGAP_STATUS_EOF
    cdef libGAP_ExecStatus libGAP_STATUS_ERROR
    cdef libGAP_ExecStatus libGAP_STATUS_QQUIT

cdef extern from "gap/objects.h":
    ctypedef void* libGAP_Obj
    libGAP_Obj libGAP_SHALLOW_COPY_OBJ(libGAP_Obj obj)
    bint libGAP_IS_INTOBJ(libGAP_Obj obj)
    libGAP_Obj libGAP_INTOBJ_INT(libGAP_Int)
    libGAP_Int libGAP_INT_INTOBJ(libGAP_Obj)
    libGAP_UInt libGAP_TNUM_OBJ(libGAP_Obj obj)
    char* libGAP_TNAM_OBJ(libGAP_Obj obj)
    cdef int libGAP_FIRST_REAL_TNUM
    cdef int libGAP_FIRST_CONSTANT_TNUM
    cdef int libGAP_T_INT
    cdef int libGAP_T_INTPOS
    cdef int libGAP_T_INTNEG
    cdef int libGAP_T_RAT
    cdef int libGAP_T_CYC
    cdef int libGAP_T_FFE
    cdef int libGAP_T_PERM2
    cdef int libGAP_T_PERM4
    cdef int libGAP_T_BOOL
    cdef int libGAP_T_CHAR
    cdef int libGAP_T_FUNCTION
    cdef int libGAP_T_FLAGS
    cdef int libGAP_T_MACFLOAT
    cdef int libGAP_T_RESERVED_BY_GAP
    cdef int libGAP_LAST_CONSTANT_TNUM
    cdef int libGAP_IMMUTABLE
    cdef int libGAP_FIRST_IMM_MUT_TNUM
    cdef int libGAP_FIRST_RECORD_TNUM
    cdef int libGAP_T_PREC
    cdef int libGAP_LAST_RECORD_TNUM
    cdef int libGAP_FIRST_LIST_TNUM
    cdef int libGAP_FIRST_PLIST_TNUM
    cdef int libGAP_T_PLIST
    cdef int libGAP_T_PLIST_NDENSE
    cdef int libGAP_T_PLIST_DENSE
    cdef int libGAP_T_PLIST_DENSE_NHOM
    cdef int libGAP_T_PLIST_DENSE_NHOM_SSORT
    cdef int libGAP_T_PLIST_DENSE_NHOM_NSORT
    cdef int libGAP_T_PLIST_EMPTY
    cdef int libGAP_T_PLIST_HOM
    cdef int libGAP_T_PLIST_HOM_NSORT
    cdef int libGAP_T_PLIST_HOM_SSORT
    cdef int libGAP_T_PLIST_TAB
    cdef int libGAP_T_PLIST_TAB_NSORT
    cdef int libGAP_T_PLIST_TAB_SSORT
    cdef int libGAP_T_PLIST_TAB_RECT
    cdef int libGAP_T_PLIST_TAB_RECT_NSORT
    cdef int libGAP_T_PLIST_TAB_RECT_SSORT
    cdef int libGAP_T_PLIST_CYC
    cdef int libGAP_T_PLIST_CYC_NSORT
    cdef int libGAP_T_PLIST_CYC_SSORT
    cdef int libGAP_T_PLIST_FFE
    cdef int libGAP_LAST_PLIST_TNUM
    cdef int libGAP_T_RANGE_NSORT
    cdef int libGAP_T_RANGE_SSORT
    cdef int libGAP_T_BLIST
    cdef int libGAP_T_BLIST_NSORT
    cdef int libGAP_T_BLIST_SSORT
    cdef int libGAP_T_STRING
    cdef int libGAP_T_STRING_NSORT
    cdef int libGAP_T_STRING_SSORT
    cdef int libGAP_LAST_LIST_TNUM
    cdef int libGAP_LAST_IMM_MUT_TNUM
    cdef int libGAP_FIRST_EXTERNAL_TNUM
    cdef int libGAP_T_COMOBJ
    cdef int libGAP_T_POSOBJ
    cdef int libGAP_T_DATOBJ
    cdef int libGAP_T_WPOBJ
    cdef int libGAP_LAST_EXTERNAL_TNUM
    cdef int libGAP_LAST_REAL_TNUM
    cdef int libGAP_LAST_VIRTUAL_TNUM
    cdef int libGAP_FIRST_COPYING_TNUM
    cdef int libGAP_COPYING
    cdef int libGAP_LAST_COPYING_TNUM
    cdef int libGAP_FIRST_TESTING_TNUM
    cdef int libGAP_TESTING
    cdef int libGAP_LAST_TESTING_TNUM

cdef extern from "gap/read.h":
    void* libGAP_ReadEvalCommand(libGAP_Obj context, libGAP_UInt *dualSemicolon)
    void* libGAP_ReadEvalFile()
    void* libGAP_ReadEvalResult
    bint libGAP_READ_ERROR()

cdef extern from "gap/scanner.h":
    void libGAP_ClearError()
    libGAP_UInt libGAP_NrError
    libGAP_UInt libGAP_Symbol
    void libGAP_GetSymbol()
    void libGAP_Match (libGAP_UInt symbol, char* msg, libGAP_UInt skipto)
    int libGAP_S_ILLEGAL
    int libGAP_S_IDENT
    int libGAP_S_UNBIND
    int libGAP_S_ISBOUND
    int libGAP_S_TRYNEXT
    int libGAP_S_INFO
    int libGAP_S_ASSERT
    int libGAP_S_SAVEWS
    int libGAP_S_LOADWS
    int libGAP_S_LBRACK
    int libGAP_S_LBRACE
    int libGAP_S_BLBRACK
    int libGAP_S_BLBRACE
    int libGAP_S_RBRACK
    int libGAP_S_RBRACE
    int libGAP_S_DOT
    int libGAP_S_BDOT
    int libGAP_S_LPAREN
    int libGAP_S_RPAREN
    int libGAP_S_COMMA
    int libGAP_S_DOTDOT
    int libGAP_S_COLON
    int libGAP_S_PARTIALINT
    int libGAP_S_INT
    int libGAP_S_TRUE
    int libGAP_S_FALSE
    int libGAP_S_CHAR
    int libGAP_S_STRING
    int libGAP_S_PARTIALSTRING
    int libGAP_S_REC
    int libGAP_S_FUNCTION
    int libGAP_S_LOCAL
    int libGAP_S_END
    int libGAP_S_MAPTO
    int libGAP_S_MULT
    int libGAP_S_DIV
    int libGAP_S_MOD
    int libGAP_S_POW
    int libGAP_S_PLUS
    int libGAP_S_MINUS
    int libGAP_S_EQ
    int libGAP_S_LT
    int libGAP_S_GT
    int libGAP_S_NE
    int libGAP_S_LE
    int libGAP_S_GE
    int libGAP_S_IN
    int libGAP_S_NOT
    int libGAP_S_AND
    int libGAP_S_OR
    int libGAP_S_ASSIGN
    int libGAP_S_IF
    int libGAP_S_FOR
    int libGAP_S_WHILE
    int libGAP_S_REPEAT
    int libGAP_S_THEN
    int libGAP_S_ELIF
    int libGAP_S_ELSE
    int libGAP_S_FI
    int libGAP_S_DO
    int libGAP_S_OD
    int libGAP_S_UNTIL
    int libGAP_S_BREAK
    int libGAP_S_RETURN
    int libGAP_S_QUIT
    int libGAP_S_QQUIT
    int libGAP_S_CONTINUE
    int libGAP_S_SEMICOLON
    int libGAP_S_EOF

cdef extern from "gap/gvars.h":
    libGAP_UInt libGAP_GVarName(char* name)
    void libGAP_AssGVar(libGAP_UInt gvar, libGAP_Obj val)
    libGAP_Obj libGAP_VAL_GVAR(libGAP_UInt gvar)

cdef extern from "gap/string.h":
    char* libGAP_CSTR_STRING(libGAP_Obj list)
    int libGAP_GET_LEN_STRING(libGAP_Obj list)
    bint libGAP_IS_STRING(libGAP_Obj obj)
    bint libGAP_IsStringConv(libGAP_Obj obj)
    bint libGAP_ConvString(libGAP_Obj obj)
    void libGAP_C_NEW_STRING(libGAP_Obj new_gap_string, int length, char* c_string)

cdef extern from "gap/gasman.h":
    void libGAP_InitGlobalBag(libGAP_Obj* addr, char* cookie)
    libGAP_Obj libGAP_NewBag(libGAP_UInt type, libGAP_UInt size)
    void libGAP_CHANGED_BAG(libGAP_Obj bag)
    void libGAP_MARK_BAG(libGAP_Obj bag)
    bint libGAP_IS_MARKED_ALIVE(libGAP_Obj bag)
    bint libGAP_IS_MARKED_DEAD(libGAP_Obj bag)
    bint libGAP_IS_MARKED_HALFDEAD(libGAP_Obj bag)
    cdef libGAP_UInt libGAP_NrAllBags
    cdef libGAP_UInt libGAP_SizeAllBags
    cdef libGAP_UInt libGAP_NrLiveBags
    cdef libGAP_UInt libGAP_SizeLiveBags
    cdef libGAP_UInt libGAP_NrDeadBags
    cdef libGAP_UInt libGAP_SizeDeadBags
    cdef libGAP_UInt libGAP_NrHalfDeadBags
    libGAP_UInt libGAP_CollectBags(libGAP_UInt size, libGAP_UInt full)
    void libGAP_CallbackForAllBags(void (*func)(libGAP_Obj))
    char* libGAP_TNAM_BAG(libGAP_Obj obj)
    libGAP_UInt libGAP_TNUM_BAG(libGAP_Obj)
    libGAP_UInt libGAP_SIZE_BAG(libGAP_Obj)
    void libGAP_CheckMasterPointers()
    libGAP_Obj* libGAP_MptrBags
    libGAP_Obj* libGAP_YoungBags
    libGAP_Obj* libGAP_OldBags
    libGAP_Obj* libGAP_AllocBags
    libGAP_Obj* libGAP_MarkedBags
    libGAP_Obj* libGAP_ChangedBags

# in gasman.c but not declared in gasman.h
cdef extern libGAP_Obj* libGAP_StopBags
cdef extern libGAP_Obj* libGAP_EndBags

cdef extern from "gap/ariths.h":
    libGAP_Obj libGAP_SUM (libGAP_Obj, libGAP_Obj)
    libGAP_Obj libGAP_DIFF(libGAP_Obj, libGAP_Obj)
    libGAP_Obj libGAP_PROD(libGAP_Obj, libGAP_Obj)
    libGAP_Obj libGAP_QUO(libGAP_Obj, libGAP_Obj)
    libGAP_Obj libGAP_POW(libGAP_Obj, libGAP_Obj)
    libGAP_Obj libGAP_MOD(libGAP_Obj, libGAP_Obj)
    libGAP_Obj libGAP_CALL_0ARGS(libGAP_Obj f)              # 0 arguments
    libGAP_Obj libGAP_CALL_1ARGS(libGAP_Obj f, libGAP_Obj a1)      # 1 argument
    libGAP_Obj libGAP_CALL_2ARGS(libGAP_Obj f, libGAP_Obj a1, libGAP_Obj a2)
    libGAP_Obj libGAP_CALL_3ARGS(libGAP_Obj f, libGAP_Obj a1, libGAP_Obj a2, libGAP_Obj a3)
    libGAP_Obj libGAP_CALL_4ARGS(libGAP_Obj f, libGAP_Obj a1, libGAP_Obj a2, libGAP_Obj a3,
                                 libGAP_Obj a4)
    libGAP_Obj libGAP_CALL_5ARGS(libGAP_Obj f, libGAP_Obj a1, libGAP_Obj a2, libGAP_Obj a3,
                                 libGAP_Obj a4, libGAP_Obj a5)
    libGAP_Obj libGAP_CALL_6ARGS(libGAP_Obj f, libGAP_Obj a1, libGAP_Obj a2, libGAP_Obj a3,
                                 libGAP_Obj a4, libGAP_Obj a5, libGAP_Obj a6)
    libGAP_Obj libGAP_CALL_XARGS(libGAP_Obj f, libGAP_Obj args)   # more than 6 arguments
    bint libGAP_EQ(libGAP_Obj opL, libGAP_Obj opR)
    bint libGAP_LT(libGAP_Obj opL, libGAP_Obj opR)

cdef extern from "gap/calls.h":
    bint libGAP_IS_FUNC(libGAP_Obj)

cdef extern from "gap/plist.h":
    libGAP_Obj libGAP_NEW_PLIST(int type, int len)
    bint libGAP_IS_PLIST(libGAP_Obj lst)
    int libGAP_LEN_PLIST(libGAP_Obj lst)
    libGAP_Obj libGAP_ELM_PLIST(libGAP_Obj lst, int pos)

cdef extern from "gap/lists.h":
    void libGAP_UNB_LIST(libGAP_Obj list, int pos)
    bint libGAP_IS_LIST(libGAP_Obj lst)
    int libGAP_LEN_LIST(libGAP_Obj lst)
    libGAP_Obj libGAP_ELM_LIST(libGAP_Obj lst, int pos)

cdef extern from "gap/listfunc.h":
    void libGAP_AddList(libGAP_Obj list, libGAP_Obj obj)
    void libGAP_AddPlist(libGAP_Obj list, libGAP_Obj obj)

cdef extern from "gap/records.h":
    char* libGAP_NAME_RNAM(libGAP_UInt rnam)
    libGAP_UInt libGAP_RNamIntg(int i)
    bint libGAP_IS_REC(libGAP_Obj obj)
    libGAP_Obj libGAP_ELM_REC(libGAP_Obj rec, libGAP_UInt rnam)
    libGAP_UInt libGAP_RNamName(libGAP_Char* name)

cdef extern from "gap/precord.h":
    libGAP_Obj libGAP_NEW_PREC(int len)
    int libGAP_LEN_PREC(libGAP_Obj rec)
    int libGAP_GET_RNAM_PREC(libGAP_Obj rec, int i)
    libGAP_Obj libGAP_GET_ELM_PREC(libGAP_Obj rec, int i)
    void libGAP_AssPRec(libGAP_Obj rec, libGAP_UInt rnam, libGAP_Obj val)
    void libGAP_UnbPRec(libGAP_Obj rec, libGAP_UInt rnam)
    bint libGAP_IsbPRec(libGAP_Obj rec, libGAP_UInt rnam)
    libGAP_Obj libGAP_ElmPRec(libGAP_Obj rec, libGAP_UInt rnam)

cdef extern from "gap/cyclotom.h":
    pass

cdef extern from "gap/bool.h":
    cdef libGAP_Obj libGAP_True
    cdef libGAP_Obj libGAP_False

cdef extern from "gap/vars.h":
     cdef int libGAP_T_LVARS
     libGAP_Obj libGAP_BottomLVars





