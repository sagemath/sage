# sage_setup: distribution = sagemath-glpk
# distutils: libraries = glpk z gmp

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

from sage.libs.glpk.types cimport *

cdef extern from "glpk.h":
    void glp_init_iocp(glp_iocp *)
    void glp_init_smcp(glp_smcp *)
    glp_prob * glp_create_prob()
    void glp_set_prob_name(glp_prob *, char *)
    const char* glp_get_prob_name(glp_prob *)
    void glp_set_obj_dir(glp_prob *, int)
    void glp_add_rows(glp_prob *, int)
    void glp_add_cols(glp_prob *, int)
    void glp_del_rows(glp_prob *, int, int *)
    void glp_set_row_name(glp_prob *, int, char *)
    void glp_set_col_name(glp_prob *, int, char *)
    void glp_set_row_bnds(glp_prob *, int, int, double, double)
    void glp_set_col_bnds(glp_prob *, int, int, double, double)
    void glp_set_obj_coef(glp_prob *, int, double)
    void glp_set_row_stat(glp_prob *, int, int)
    void glp_set_col_stat(glp_prob *, int, int)
    int glp_warm_up(glp_prob *)
    void glp_load_matrix(glp_prob *, int, int *, int *, double *)
    int glp_simplex(glp_prob *, glp_smcp *)
    int glp_exact(glp_prob *, glp_smcp *)
    int glp_intopt(glp_prob *, glp_iocp *)
    int lpx_intopt(glp_prob *)
    void glp_std_basis(glp_prob *)
    void glp_delete_prob(glp_prob *)
    double glp_get_col_prim(glp_prob *, int)
    double glp_get_obj_val(glp_prob *)
    double glp_get_col_dual(glp_prob *, int)
    double glp_get_row_prim(glp_prob *, int)
    double glp_get_row_dual(glp_prob *, int)
    int glp_get_col_stat(glp_prob *, int)
    int glp_get_row_stat(glp_prob *, int)
    int glp_print_ranges(glp_prob *lp, int,int, int, char *fname)
    int glp_get_num_rows(glp_prob *)
    int glp_get_num_cols(glp_prob *)
    double glp_mip_col_val(glp_prob *, int)
    double glp_mip_obj_val(glp_prob *)
    void glp_set_col_kind(glp_prob *, int, int)
    int glp_write_mps(glp_prob *lp, int fmt, void *parm, char *fname)
    int glp_write_lp(glp_prob *lp, void *parm, char *fname)
    void glp_set_prob_name(glp_prob *lp, char *name)
    void glp_set_obj_name(glp_prob *lp, char *name)
    void glp_set_row_name(glp_prob *lp, int i, char *name)
    void glp_set_col_name(glp_prob *lp, int i, char *name)

    double glp_get_row_ub(glp_prob *lp, int i)
    double glp_get_row_lb(glp_prob *lp, int i)

    double glp_get_col_ub(glp_prob *lp, int i)
    double glp_get_col_lb(glp_prob *lp, int i)

    void glp_set_col_ub(glp_prob *lp, int i, double value)
    void glp_set_col_lb(glp_prob *lp, int i, double value)

    int glp_eval_tab_row(glp_prob *lp, int k, int ind[], double val[])
    int glp_eval_tab_col(glp_prob *lp, int k, int ind[], double val[])

    const char* glp_get_row_name(glp_prob *lp, int i)
    const char* glp_get_col_name(glp_prob *lp, int i)

    void glp_create_index(glp_prob *lp)

    int glp_get_prim_stat(glp_prob *lp)
    int glp_get_status(glp_prob *lp)
    int glp_mip_status(glp_prob *lp)
    int glp_set_mat_row(glp_prob *lp, int, int, int *, double * )
    int glp_set_mat_col(glp_prob *lp, int, int, int *, double * )
    int glp_get_mat_row(glp_prob *lp, int, int *, double * )
    double glp_get_row_ub(glp_prob *lp, int)
    double glp_get_row_lb(glp_prob *lp, int)
    int glp_get_col_kind(glp_prob *lp, int)
    double glp_get_obj_coef(glp_prob *lp, int)
    int glp_get_obj_dir(glp_prob *lp)
    void glp_copy_prob(glp_prob *dst, glp_prob *src, int names)
    double glp_ios_mip_gap(glp_tree *T)
    int glp_ios_best_node(glp_tree *tree)
    double glp_ios_node_bound(glp_tree *T, int p)

    int glp_bf_exists(glp_prob *lp)
