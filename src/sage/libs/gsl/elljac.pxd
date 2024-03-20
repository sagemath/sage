# sage_setup: distribution = sagemath-modules
# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
cdef extern from "gsl/gsl_sf_elljac.h":

  int  gsl_sf_elljac_e(double u, double m, double * sn, double * cn, double * dn)

