#*****************************************************************************
#       Copyright (C) 2009-2013 Mike Hansen <mhansen@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "sage/ext/cdefs.pxi"

cdef extern from "ccobject.h":
    #We do this to get access to Construct_p, etc.
    pass

cdef extern from "globals.h":
    ctypedef unsigned long Ulong

#cdef extern from "klsupport.h":
#ctypedef unsigned short KLCoeff

########
# Bits #
########
cdef extern from "bits.h":
    ctypedef Ulong LFlags

cdef extern from "coxtypes.h":
    ctypedef unsigned short Rank
    ctypedef Ulong CoxSize         # should hold at least 32 bits */
    ctypedef Ulong BettiNbr        # should hold at least 32 bits */
    ctypedef unsigned CoxNbr         # should fit into a CoxSize */
    ctypedef CoxNbr ParSize          # this should not be changed */
    ctypedef unsigned short ParNbr   # should fit into a CoxNbr */
    ctypedef ParNbr *CoxArr
    ctypedef unsigned char CoxLetter # for string representations */
    ctypedef CoxLetter Generator     # internal representation of generators*/
    ctypedef unsigned short Length
    ctypedef Ulong StarOp          # for numbering star operations */


    #################
    #    CoxWord    #
    #################
    ctypedef struct c_CoxWord "coxtypes::CoxWord":
        CoxLetter get_index "operator[]"(Length j)
        Length length()
        void setLength(Length n)
        c_CoxWord append_letter "append" (CoxLetter a)
        c_CoxWord reset()
        c_CoxWord set "operator="(c_CoxWord w)

    void CoxWord_destruct "Destruct<coxtypes::CoxWord>"(c_CoxWord *mem)
    c_CoxWord* CoxWord_construct "Construct<coxtypes::CoxWord>" (void *mem)


################
#    String    #
################
cdef extern from "io.h":
    ctypedef struct c_String "io::String":
        void setLength(Ulong n)
        Ulong length()
        char* ptr()
        void setData(char* source, Ulong r)

    c_String* String_new "New<io::String>"()
    c_String* String_construct_str "Construct_p<io::String, char*>"(void *mem, char* s)
    void String_destruct "Destruct<io::String>"(c_String *mem)


##############
#    Type    #
##############
cdef extern from "type.h":
    ctypedef struct c_Type "type::Type":
        c_String name()
    void Type_destruct "Destruct<type::Type>"(c_Type *mem)
    c_Type* Type_construct "Construct_p<type::Type, char*>" (void *mem, char* s)
    c_Type* Type_construct_str "Construct_p<type::Type, char*>" (void *mem, char* s)

##############
#    Bits    #
##############
cdef extern from "bits.h":
    Generator firstBit "bits::firstBit"(Ulong n)

##################
#    CoxGraph    #
##################
cdef extern from "graph.h":
    ctypedef unsigned short CoxEntry
    ctypedef struct c_CoxGraph "graph::CoxGraph":
        pass

###############
#    KLPol    #
###############
cdef extern from "kl.h":
    ctypedef struct c_KLPol "kl::KLPol":
        pass

cdef extern from "polynomials.h":
    c_String klpoly_append "polynomials::append"(c_String str, c_KLPol, char* x)

##################
#    List    #
##################
cdef extern from "list.h":
    ctypedef struct c_List_CoxWord "list::List<coxtypes::CoxWord>":
        c_CoxWord get_index "operator[]"(Length j)
        Ulong size()
    c_List_CoxWord  c_List_CoxWord_factory "list::List<coxtypes::CoxWord>"(Ulong len)
    void List_CoxWord_destruct "Destruct< list::List<coxtypes::CoxWord> >"(c_List_CoxWord *mem)


###################
#     CoxGroup    #
###################
cdef extern from "coxgroup.h":
    ctypedef struct c_CoxGroup "coxgroup::CoxGroup":
        Rank rank()
        c_Type type()
        CoxSize order()
        c_CoxWord reduced(c_CoxWord res, c_CoxWord w)
        c_CoxWord normalForm(c_CoxWord w)
        c_CoxWord prod(c_CoxWord g, c_CoxWord h)
        int prod_nbr "prod"(c_CoxWord g, CoxNbr x)

        c_CoxGraph graph()
        CoxEntry M(Generator s, Generator t)
        CoxNbr extendContext(c_CoxWord w)
        c_KLPol klPol(CoxNbr x, CoxNbr y)
        bint inOrder_word "inOrder"(c_CoxWord u, c_CoxWord w)
        bint inOrder_nbr "inOrder"(CoxNbr x, CoxNbr y)

        LFlags descent(c_CoxWord w)
        LFlags ldescent(c_CoxWord w)
        LFlags rdescent(c_CoxWord w)
        bint isDescent(c_CoxWord w, Generator s)

        void coatoms(c_List_CoxWord, c_CoxWord u)

        unsigned short mu(CoxNbr x, CoxNbr y) #KLCoeff

    c_CoxGroup* CoxGroup_construct_type_rank "Construct_pp<coxgroup::CoxGroup, type::Type, Rank>" (void *mem, c_Type t, Rank r)
    void CoxGroup_delete "Delete<coxgroup::CoxGroup>"(c_CoxGroup *mem)


#####################
#    Interactive    #
#####################
cdef extern from "interactive.h":
    c_CoxGroup* coxeterGroup "interactive::coxeterGroup"(c_Type x, Rank l)

cdef extern from "constants.h":
    void initConstants()

###############################
#    Finite Coxeter groups    #
###############################
cdef extern from "fcoxgroup.h":
    ctypedef struct c_FiniteCoxGroup "fcoxgroup::FiniteCoxGroup":
        void fullContext()
        bint isFullContext()
        Length maxLength()
        c_CoxWord longest_coxword()

    bint isFiniteType "fcoxgroup::isFiniteType"(c_CoxGroup *W)


######################
#    Sage specific   #
######################
cdef extern from "sage.h":
    void interval "sage::interval"(c_List_CoxWord, c_CoxGroup W, c_CoxWord g, c_CoxWord h)
