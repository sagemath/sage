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

from libc.stdio cimport FILE

cdef extern from "glpk.h":
    ctypedef struct glp_tree

    ctypedef struct glp_prob

    ctypedef struct glp_iocp:
        int msg_lev
        int br_tech
        int bt_tech
        int pp_tech
        int fp_heur
        int gmi_cuts
        int mir_cuts
        int cov_cuts
        int clq_cuts
        double tol_int
        double tol_obj
        double mip_gap
        int tm_lim
        int out_frq
        int out_dly
        int presolve
        int binarize
        void (*cb_func)(glp_tree *T, void *info) # callback function
        void *cb_info                            # callback function input

    ctypedef struct glp_smcp:
        int msg_lev
        int meth
        int pricing
        int r_test
        double tol_bnd
        double tol_dj
        double tol_piv
        double obj_ll
        double obj_ul
        int it_lim
        int tm_lim
        int out_frq
        int out_dly
        int presolve

    # Graph structure
    ctypedef struct glp_graph:
        void *pool
        char *name
        int nv_max
        int nv
        int na
        glp_vertex **v
        void *index
        int v_size
        int a_size

    # Arc structure
    ctypedef struct glp_arc:
        glp_vertex *tail
        glp_vertex *head
        void *data
        void *temp
        glp_arc *t_prev
        glp_arc *t_next
        glp_arc *h_prev
        glp_arc *h_next

    # Vertex structure
    ctypedef struct glp_vertex:
        int i
        char *name
        void *entry
        void *data
        void *temp
        glp_arc *out
