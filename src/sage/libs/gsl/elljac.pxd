# distutils: libraries = GSL_LIBRARIES
cdef extern from "gsl/gsl_sf_elljac.h":

  int  gsl_sf_elljac_e(double u, double m, double * sn, double * cn, double * dn)

