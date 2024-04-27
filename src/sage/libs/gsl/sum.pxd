# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR

cdef extern from "gsl/gsl_sum.h":

  ctypedef struct gsl_sum_levin_u_workspace:
    size_t terms_used
    double sum_plain

  gsl_sum_levin_u_workspace *gsl_sum_levin_u_alloc (size_t n)
  void gsl_sum_levin_u_free (gsl_sum_levin_u_workspace * w)


  int gsl_sum_levin_u_accel ( double *array,
                              size_t n,
                             gsl_sum_levin_u_workspace * w,
                             double *sum_accel, double *abserr)


  int gsl_sum_levin_u_minmax ( double *array,
                               size_t n,
                               size_t min_terms,
                               size_t max_terms,
                              gsl_sum_levin_u_workspace * w,
                              double *sum_accel, double *abserr)


  int gsl_sum_levin_u_step ( double term,
                         size_t n,
                         size_t nmax,
                        gsl_sum_levin_u_workspace * w,
                        double *sum_accel)


  ctypedef struct gsl_sum_levin_utrunc_workspace

  gsl_sum_levin_utrunc_workspace *gsl_sum_levin_utrunc_alloc (size_t n)
  void gsl_sum_levin_utrunc_free (gsl_sum_levin_utrunc_workspace * w)

  int gsl_sum_levin_utrunc_accel ( double *array,  size_t n,
                                  gsl_sum_levin_utrunc_workspace * w,
                                  double *sum_accel, double *abserr_trunc)

  int gsl_sum_levin_utrunc_minmax ( double *array,  size_t n,
                                    size_t min_terms,
                                    size_t max_terms,
                                   gsl_sum_levin_utrunc_workspace * w,
                                   double *sum_accel, double *abserr_trunc)

  int gsl_sum_levin_utrunc_step ( double term,  size_t n,
                                 gsl_sum_levin_utrunc_workspace * w,
                                 double *sum_accel)

