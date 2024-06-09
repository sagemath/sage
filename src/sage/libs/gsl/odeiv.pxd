# sage_setup: distribution = sagemath-modules
# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
cdef extern from "gsl/gsl_odeiv.h":

  ctypedef struct gsl_odeiv_system:
    int (* function) (double t,  double y[], double dydt[], void * params)
    int (* jacobian) (double t,  double y[], double * dfdy, double dfdt[], void * params)
    size_t dimension
    void * params

  # no:
  #define GSL_ODEIV_FN_EVAL(S,t,y,f)  (*((S)->function))(t,y,f,(S)->params)
  #define GSL_ODEIV_JA_EVAL(S,t,y,dfdy,dfdt)  (*((S)->jacobian))(t,y,dfdy,dfdt,(S)->params)


  ctypedef struct gsl_odeiv_step_type

  ctypedef struct gsl_odeiv_step

  gsl_odeiv_step_type *gsl_odeiv_step_rk2
  gsl_odeiv_step_type *gsl_odeiv_step_rk4
  gsl_odeiv_step_type *gsl_odeiv_step_rkf45
  gsl_odeiv_step_type *gsl_odeiv_step_rkck
  gsl_odeiv_step_type *gsl_odeiv_step_rk8pd
  gsl_odeiv_step_type *gsl_odeiv_step_rk2imp
  gsl_odeiv_step_type *gsl_odeiv_step_rk4imp
  gsl_odeiv_step_type *gsl_odeiv_step_bsimp
  gsl_odeiv_step_type *gsl_odeiv_step_gear1
  gsl_odeiv_step_type *gsl_odeiv_step_gear2


  gsl_odeiv_step * gsl_odeiv_step_alloc( gsl_odeiv_step_type * T, size_t dim)
  int  gsl_odeiv_step_reset(gsl_odeiv_step * s)
  void gsl_odeiv_step_free(gsl_odeiv_step * s)

  char * gsl_odeiv_step_name( gsl_odeiv_step *)
  unsigned int gsl_odeiv_step_order( gsl_odeiv_step * s)

  int  gsl_odeiv_step_apply(gsl_odeiv_step *, double t, double h, double y[], double yerr[],  double dydt_in[], double dydt_out[],  gsl_odeiv_system * dydt)


  ctypedef struct gsl_odeiv_control_type

  ctypedef struct gsl_odeiv_control
  cdef enum:
    GSL_ODEIV_HADJ_DEC = -1
    GSL_ODEIV_HADJ_NIL = 0
    GSL_ODEIV_HADJ_INC = 1

  gsl_odeiv_control * gsl_odeiv_control_alloc( gsl_odeiv_control_type * T)
  int gsl_odeiv_control_init(gsl_odeiv_control * c, double eps_abs, double eps_rel, double a_y, double a_dydt)
  void gsl_odeiv_control_free(gsl_odeiv_control * c)
  int gsl_odeiv_control_hadjust (gsl_odeiv_control * c, gsl_odeiv_step * s,  double y0[],  double yerr[],  double dydt[], double * h)
  char * gsl_odeiv_control_name( gsl_odeiv_control * c)


  gsl_odeiv_control * gsl_odeiv_control_standard_new(double eps_abs, double eps_rel, double a_y, double a_dydt)
  gsl_odeiv_control * gsl_odeiv_control_y_new(double eps_abs, double eps_rel)
  gsl_odeiv_control * gsl_odeiv_control_yp_new(double eps_abs, double eps_rel)

  gsl_odeiv_control * gsl_odeiv_control_scaled_new(double eps_abs, double eps_rel, double a_y, double a_dydt,  double scale_abs[], size_t dim)

  ctypedef struct gsl_odeiv_evolve

  gsl_odeiv_evolve * gsl_odeiv_evolve_alloc(size_t dim)
  int gsl_odeiv_evolve_apply(gsl_odeiv_evolve *, gsl_odeiv_control * con, gsl_odeiv_step * step,  gsl_odeiv_system * dydt, double * t, double t1, double * h, double y[])
  int gsl_odeiv_evolve_reset(gsl_odeiv_evolve *)
  void gsl_odeiv_evolve_free(gsl_odeiv_evolve *)
