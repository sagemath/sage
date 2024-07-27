# sage_setup: distribution = sagemath-modules
# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR

cdef extern from "gsl/gsl_fit.h":


  int gsl_fit_linear ( double * x,  size_t xstride,
                       double * y,  size_t ystride,
                       size_t n,
                      double * c0, double * c1,
                      double * cov00, double * cov01, double * cov11,
                      double * sumsq)


  int gsl_fit_wlinear ( double * x,  size_t xstride,
                        double * w,  size_t wstride,
                        double * y,  size_t ystride,
                        size_t n, double * c0, double * c1,
                       double * cov00, double * cov01, double * cov11,
                       double * chisq)

  int gsl_fit_linear_est ( double x, double c0,  double c1,
                       double c00,  double c01,  double c11,
                      double *y, double *y_err)


  int gsl_fit_mul ( double * x,  size_t xstride,
                    double * y,  size_t ystride,
                    size_t n, double * c1, double * cov11, double * sumsq)

  int gsl_fit_wmul ( double * x,  size_t xstride,
                     double * w,  size_t wstride,
                     double * y,  size_t ystride,
                   size_t n, double * c1, double * cov11, double * sumsq)


  int gsl_fit_mul_est ( double x, double c1, double c11,
                   double *y, double *y_err)


  int gsl_fit_poly ( double * x, double * w, double * y, size_t n,
                    double * c, size_t m, double * chisq)

  int gsl_fit_fns ( double * A, double * w, double * y,
                   size_t n, double * c, size_t m, double * chisq)

  int gsl_fit_linear_nd (double * m, double * y, double * w)

