# sage_setup: distribution = sagemath-coxeter3

#*****************************************************************************
#       Copyright (C) 2009-2013 Mike Hansen <mhansen@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

cdef extern from "coxeter/globals.h":
    ctypedef unsigned long Ulong


########
# Bits #
########
cdef extern from "coxeter/bits.h":
    ctypedef Ulong LFlags

cdef extern from "coxeter/coxtypes.h" namespace "coxtypes":
    ctypedef unsigned short Rank
    ctypedef Ulong CoxSize            # should hold at least 32 bits
    ctypedef Ulong BettiNbr           # should hold at least 32 bits
    ctypedef unsigned CoxNbr          # should fit into a CoxSize
    ctypedef CoxNbr ParSize           # this should not be changed
    ctypedef unsigned short ParNbr    # should fit into a CoxNbr
    ctypedef ParNbr *CoxArr
    ctypedef unsigned char CoxLetter  # for string representations
    ctypedef CoxLetter Generator      # internal representation of generators
    ctypedef unsigned short Length
    ctypedef Ulong StarOp             # for numbering star operations

    #################
    #    CoxWord    #
    #################
    cdef cppclass c_CoxWord "coxtypes::CoxWord":
        CoxLetter operator[](Length j)
        Length length()
        void setLength(Length n)
        c_CoxWord append(CoxLetter& a)
        c_CoxWord reset()


################
#    String    #
################
cdef extern from "coxeter/io.h" namespace "io":
    cdef cppclass c_String "io::String":
        c_String()
        c_String(char* s)
        void setLength(Ulong& n)
        Ulong length()
        char* ptr()
        void setData(char* source, Ulong r)


##############
#    Type    #
##############
cdef extern from "coxeter/type.h":
    cdef cppclass c_Type "coxeter::Type":
        c_Type()
        c_Type(char* s)
        c_String name()


##############
#    Bits    #
##############
cdef extern from "coxeter/bits.h" namespace "bits":
    Generator firstBit(Ulong n)

##################
#    CoxGraph    #
##################
cdef extern from "coxeter/graph.h" namespace "graph":
    ctypedef unsigned short CoxEntry
    ctypedef struct c_CoxGraph "graph::CoxGraph":
        pass

###############
#    KLPol    #
###############
cdef extern from "coxeter/kl.h" namespace "kl":
    cdef cppclass c_KLPol "kl::KLPol":
        const unsigned short& operator[](Ulong j)
        unsigned long deg()
        int isZero()

cdef extern from "coxeter/polynomials.h" namespace "polynomials":
    c_String klpoly_append "polynomials::append"(c_String str, c_KLPol, char* x)

##################
#    List    #
##################
cdef extern from "coxeter/list.h" namespace "list":
    cdef cppclass c_List_CoxWord "list::List<coxtypes::CoxWord> ":
        c_List_CoxWord()
        c_List_CoxWord(Ulong len)
        c_CoxWord operator[](Ulong j)
        Ulong size()


###################
#     CoxGroup    #
###################
cdef extern from "coxeter/coxgroup.h":
    cdef cppclass c_CoxGroup "coxeter::CoxGroup":
        c_CoxGroup()
        c_CoxGroup(c_Type t, Rank r)

        Rank rank()
        c_Type type()
        CoxSize order()
        c_CoxWord reduced(c_CoxWord& res, c_CoxWord& w)
        c_CoxWord normalForm(c_CoxWord& w)
        c_CoxWord prod(c_CoxWord& g, c_CoxWord& h)
        int prod_nbr "prod"(c_CoxWord& g, CoxNbr& x)

        c_CoxGraph graph()
        CoxEntry M(Generator s, Generator t)
        CoxNbr extendContext(c_CoxWord& w)
        c_KLPol klPol(CoxNbr& x, CoxNbr& y)
        bint inOrder(c_CoxWord& u, c_CoxWord& w)
        bint inOrder(CoxNbr x, CoxNbr y)

        LFlags descent(c_CoxWord& w)
        LFlags ldescent(c_CoxWord& w)
        LFlags rdescent(c_CoxWord& w)
        bint isDescent(c_CoxWord& w, Generator& s)

        void coatoms(c_List_CoxWord& l, c_CoxWord& u)

        unsigned short mu(CoxNbr& x, CoxNbr& y)


#####################
#    Interactive    #
#####################
cdef extern from "coxeter/interactive.h" namespace "interactive":
    c_CoxGroup* coxeterGroup(c_Type x, Rank l)

cdef extern from "coxeter/constants.h":
    void initConstants()

###############################
#    Finite Coxeter groups    #
###############################
cdef extern from "coxeter/fcoxgroup.h" namespace "fcoxgroup":
    ctypedef struct c_FiniteCoxGroup "fcoxgroup::FiniteCoxGroup":
        void fullContext()
        bint isFullContext()
        Length maxLength()
        c_CoxWord longest_coxword()

    bint isFiniteType(c_CoxGroup *W)


######################
#    Sage specific   #
######################
cdef extern from "coxeter/sage.h" namespace "sage":
    void interval(c_List_CoxWord& l, c_CoxGroup& W, c_CoxWord& g, c_CoxWord& h)
