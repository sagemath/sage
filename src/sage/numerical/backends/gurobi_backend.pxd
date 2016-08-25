##############################################################################
#       Copyright (C) 2010 Nathann Cohen <nathann.cohen@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

from .generic_backend cimport GenericBackend


cdef extern from "gurobi_c.h":
     ctypedef struct GRBmodel:
         pass
     ctypedef struct GRBenv:
         pass

     int GRBloadenv(GRBenv **, char *)
     int GRBnewmodel(GRBenv *env, GRBmodel **modelP, char *Pname, int numvars, double *obj, double *lb, double *ub, char *vtype, char **varnames)
     GRBmodel * GRBcopymodel (GRBmodel *model)

     int GRBaddvar(GRBmodel *model, int numnz, int *vind, double *vval, double obj, double lb, double ub, char vtype, char *varname)
     int GRBaddvars (GRBmodel*model, intnumvars, intnumnz, int*vbeg, int*vind, double*vval, double*obj, double*lb, double*ub, char*vtype, char** varnames )



     int GRBaddconstr(GRBmodel *model, int numnz, int *cind, double *cval, char sense, double rhs, char *constrnames)

     int GRBdelconstrs (GRBmodel *model, int numdel, int * ind )
     int GRBaddrangeconstr(GRBmodel *model, int numnz, int *cind, double *cval, double lower, double upper, char *constrnames)


     int GRBfreemodel(GRBmodel *model)
     void GRBfreeenv(GRBenv *env)
     int GRBupdatemodel(GRBmodel *model)
     int GRBoptimize(GRBmodel *model)
     char * GRBgeterrormsg(GRBenv *env)
     int GRBgetintattr(GRBmodel *model, char * attrname, int * result)
     int GRBsetintattr(GRBmodel *model, char * attrname, int value)

     int GRBgetdblattr(GRBmodel *model, char * attrname, double * result)
     int GRBsetdblattr(GRBmodel *model, char * attrname, double value)

     int GRBgetdblattrelement (GRBmodel*model, char* attrname, int element, double* valueP )
     int GRBsetdblattrelement (GRBmodel*model, char* attrname, int element, double valueP )

     int GRBgetcharattrelement (GRBmodel*model, char* attrname, int element, char* valueP )
     int GRBsetcharattrelement (GRBmodel*model, char* attrname, int element, char valueP )
     int GRBsetstrattr (GRBmodel*model, char* attrname, char * valueP )
     int GRBgetstrattr (GRBmodel*model, char* attrname, char ** valueP )

     #int GRBsetstrattrelement (GRBmodel*model, char* attrname, int index, char * valueP )
     int GRBgetstrattrelement (GRBmodel*model, char* attrname, int index, char ** valueP )


     int GRBwrite (GRBmodel*model, char* filename)

     int GRBgetdblparam(GRBenv *env, char * attrname, double * value)
     int GRBgetintparam(GRBenv *env, char * attrname, int * value)
     int GRBgetstrparam(GRBenv *env, char * attrname, char * value)
     int GRBsetdblparam(GRBenv *env, char * attrname, double value)
     int GRBsetintparam(GRBenv *env, char * attrname, int value)
     int GRBsetstrparam(GRBenv *env, char * attrname, char * value)

     GRBenv * GRBgetenv (GRBmodel * model )

     int GRBgetconstrs (GRBmodel * model, int * numnzP, int * cbeg, int * cind, double * cval, int start, int len )

     int GRB_BINARY
     int GRB_CONTINUOUS
     int GRB_INTEGER
     double GRB_INFINITY

     char GRB_LESS_EQUAL
     char GRB_GREATER_EQUAL
     char GRB_EQUAL

     int GRB_LOADED
     int GRB_OPTIMAL
     int GRB_INFEASIBLE
     int GRB_INF_OR_UNBD
     int GRB_UNBOUNDED
     int GRB_CUTOFF
     int GRB_ITERATION_LIMIT
     int GRB_NODE_LIMIT
     int GRB_TIME_LIMIT
     int GRB_SOLUTION_LIMIT
     int GRB_INTERRUPTED
     int GRB_NUMERIC
     int GRB_SUBOPTIMAL




cdef class GurobiBackend(GenericBackend):

    cdef GRBenv * env
    cdef GRBenv * env_master
    cdef GRBmodel * model
    cpdef __copy__(self)

    cdef int num_vars


