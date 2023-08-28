# sage_setup: distribution = sagemath-glpk
#*****************************************************************************
#       Copyright (C) 2010 Nathann Cohen <nathann.cohen@gmail.com>
#       Copyright (C) 2015 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

cdef extern from "glpk.h":
    # constants for smcp control
    int GLP_MSG_OFF
    int GLP_MSG_ERR
    int GLP_MSG_ON
    int GLP_MSG_ALL

    int GLP_PRIMAL
    int GLP_DUALP
    int GLP_DUAL

    int GLP_PT_STD
    int GLP_PT_PSE

    int GLP_RT_STD
    int GLP_RT_HAR

    int GLP_ON
    int GLP_OFF

    # constants for iocp control, not already in simplex
    int GLP_BR_FFV
    int GLP_BR_LFV
    int GLP_BR_MFV
    int GLP_BR_DTH
    int GLP_BR_PCH

    int GLP_BT_DFS
    int GLP_BT_BFS
    int GLP_BT_BLB
    int GLP_BT_BPH

    int GLP_PP_NONE
    int GLP_PP_ROOT
    int GLP_PP_ALL

    # error codes
    int GLP_EBADB
    int GLP_ESING
    int GLP_ECOND
    int GLP_EBOUND
    int GLP_EFAIL
    int GLP_EOBJLL
    int GLP_EOBJUL
    int GLP_EITLIM
    int GLP_ETMLIM
    int GLP_ENOPFS
    int GLP_ENODFS
    int GLP_EROOT
    int GLP_ESTOP
    int GLP_EMIPGAP
    int GLP_EDATA
    int GLP_ERANGE

    int GLP_UNDEF
    int GLP_OPT
    int GLP_FEAS
    int GLP_NOFEAS
    int GLP_INFEAS
    int GLP_UNBND

    # other constants
    int GLP_MAX
    int GLP_MIN
    int GLP_UP
    int GLP_FR
    int GLP_DB
    int GLP_FX
    int GLP_LO
    int GLP_CV
    int GLP_IV
    int GLP_BV
    int GLP_MPS_DECK
    int GLP_MPS_FILE

    int GLP_MSG_DBG

    int GLP_BS      # basic variable
    int GLP_NL      # non-basic variable on lower bound
    int GLP_NU      # non-basic variable on upper bound
    int GLP_NF      # non-basic free (unbounded) variable
    int GLP_NS      # non-basic fixed variable
