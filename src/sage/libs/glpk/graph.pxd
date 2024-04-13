# sage_setup: distribution = sagemath-glpk
# distutils: libraries = glpk z gmp

#*****************************************************************************
#       Copyright (C) 2012 Christian Kuper <christian.kuper@t-online.de>
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
    glp_graph *glp_create_graph(int v_size, int a_size)
    void glp_set_graph_name(glp_graph *G, char *name)
    int glp_add_vertices(glp_graph *G, int nv)
    void glp_set_vertex_name(glp_graph *G, int i, char *name)
    glp_arc *glp_add_arc(glp_graph *G, int i, int j)
    void glp_del_vertices(glp_graph *G, int ndel, int num[])
    void glp_del_arc(glp_graph *G, glp_arc *a)
    void glp_delete_graph(glp_graph *G)
    void glp_create_v_index(glp_graph *G)
    int glp_find_vertex(glp_graph *G, char *name)
    int glp_read_graph(glp_graph *G, char *fname)
    int glp_write_graph(glp_graph *G, char *fname)
    int glp_read_ccdata(glp_graph *G, int v_wgt, char *fname)
    int glp_write_ccdata(glp_graph *G, int v_wgt, char *fname)
    int glp_read_mincost(glp_graph *G, int v_rhs, int a_low,
                         int a_cap, int a_cost, char *fname)
    int glp_write_mincost(glp_graph *G, int v_rhs, int a_low,
                          int a_cap, int a_cost, char *fname)
    int glp_mincost_okalg(glp_graph *G, int v_rhs, int a_low, int a_cap,
                          int a_cost, double *sol, int a_x, int v_pi)
    int glp_read_maxflow(glp_graph *G, int *s, int *t,
                         int a_cap, char *fname)
    int glp_write_maxflow(glp_graph *G, int s, int t, int a_cap, char *fname)
    int glp_maxflow_ffalg(glp_graph *G, int s, int t, int a_cap,
                          double *sol, int a_x, int v_cut)
    double glp_cpp(glp_graph *G, int v_t, int v_es, int v_ls)
