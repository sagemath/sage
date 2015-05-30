cdef extern from *:
    ctypedef double* const_double_ptr "const double*"
    ctypedef char * const_char_ptr "const char*"

cdef extern from "float.h":
    cdef double DBL_MAX

cdef extern from "glpk.h":

     # Graph structure
     ctypedef struct _glp_graph "glp_graph":
         void *pool
         char *name
         int nv_max
         int nv
         int na
         _glp_vertex **v
         void *index
         int v_size
         int a_size

     # Arc structure
     ctypedef struct _glp_arc "glp_arc":
         _glp_vertex *tail
         _glp_vertex *head
         void *data
         void *temp
         _glp_arc *t_prev
         _glp_arc *t_next
         _glp_arc *h_prev
         _glp_arc *h_next

     # Vertex structure
     ctypedef struct _glp_vertex "glp_vertex":
         int i
         char *name
         void *entry
         void *data
         void *temp
         #_glp_arc *in
         _glp_arc *out


     _glp_graph *glp_create_graph(int v_size, int a_size)
     void glp_set_graph_name(_glp_graph *G, char *name)
     int glp_add_vertices(_glp_graph *G, int nv)
     void glp_set_vertex_name(_glp_graph *G, int i, char *name)
     _glp_arc *glp_add_arc(_glp_graph *G, int i, int j)
     void glp_del_vertices(_glp_graph *G, int ndel, int num[])
     void glp_del_arc(_glp_graph *G, _glp_arc *a)
     void glp_delete_graph(_glp_graph *G)
     void glp_create_v_index(_glp_graph *G)
     int glp_find_vertex(_glp_graph *G, char *name)
     int glp_read_graph(_glp_graph *G, char *fname)
     int glp_write_graph(_glp_graph *G, char *fname)
     int glp_read_ccdata(_glp_graph *G, int v_wgt, char *fname)
     int glp_write_ccdata(_glp_graph *G, int v_wgt, char *fname)
     int glp_read_mincost(_glp_graph *G, int v_rhs, int a_low,
                          int a_cap, int a_cost, char *fname)
     int glp_write_mincost(_glp_graph *G, int v_rhs, int a_low,
                           int a_cap, int a_cost, char *fname)
     int glp_mincost_okalg(_glp_graph *G, int v_rhs, int a_low, int a_cap,
                           int a_cost, double *sol, int a_x, int v_pi)
     int glp_read_maxflow(_glp_graph *G, int *s, int *t,
                          int a_cap, char *fname)
     int glp_write_maxflow(_glp_graph *G, int s, int t, int a_cap, char *fname)
     int glp_maxflow_ffalg(_glp_graph *G, int s, int t, int a_cap,
                           double *sol, int a_x, int v_cut)
     double glp_cpp(_glp_graph *G, int v_t, int v_es, int v_ls)

     int GLP_ENOPFS
     int GLP_EDATA
     int GLP_ERANGE
     int GLP_EFAIL

ctypedef struct c_v_data:
         double rhs
         double pi
         double es
         double ls
         long cut

ctypedef struct c_a_data:
         double low
         double cap
         double cost
         double x


cdef class GLPKGraphBackend(object):

    cdef _glp_graph * graph
    cpdef add_vertex(self, char* name = ?)
    cpdef list add_vertices(self, vertices)
    cpdef __add_vertices_sage(self, g)
    cpdef dict get_vertex(self, char* vertex)
    cpdef dict get_vertices(self, verts)
    cpdef set_vertex_demand(self, char* vertex, param)
    cpdef set_vertices_demand(self, list pairs)
    cpdef list vertices(self)
    cpdef add_edge(self, char* u, char* v, dict params = ?)
    cpdef __add_edges_sage(self, g)
    cpdef list add_edges(self, edges)
    cpdef delete_edge(self, char* u, char* v, dict params = ?)
    cpdef tuple get_edge(self, char* u, char* v)
    cpdef list edges(self)
    cpdef delete_vertex(self, char* vert)
    cpdef delete_vertices(self, list verts)
    cpdef int _find_vertex(self, char *)
    cpdef int write_graph(self, char *fname)
    cpdef int write_ccdata(self, char *fname)
    cpdef int write_mincost(self, char *fname)
    cpdef double mincost_okalg(self) except -1
    cdef int s
    cdef int t
    cpdef int write_maxflow(self, char *fname) except -1
    cpdef double maxflow_ffalg(self, u = ?, v = ?) except -1
    cpdef double cpp(self)




