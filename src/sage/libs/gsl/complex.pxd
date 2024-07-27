# sage_setup: distribution = sagemath-modules
# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from sage.libs.gsl.types cimport *

cdef extern from "gsl/gsl_complex.h":
  double GSL_REAL(gsl_complex z)

  double GSL_IMAG(gsl_complex z)

  int GSL_COMPLEX_EQ(gsl_complex z1,gsl_complex z2)

  void GSL_SET_COMPLEX(gsl_complex * zp, double  x, double  y)

  void GSL_SET_REAL(gsl_complex * zp, double x)

  void GSL_SET_IMAG(gsl_complex * zp, double y)


cdef extern from "gsl/gsl_complex_math.h":

# Complex numbers
  gsl_complex  gsl_complex_rect(double x, double y)

  gsl_complex  gsl_complex_polar(double r, double theta)

  gsl_complex GSL_COMPLEX_ONE
  gsl_complex GSL_COMPLEX_ZERO
  gsl_complex GSL_COMPLEX_NEGONE

  # Properties of complex numbers
  double  gsl_complex_arg(gsl_complex z)

  double  gsl_complex_abs(gsl_complex z)

  double  gsl_complex_abs2(gsl_complex z)

  double  gsl_complex_logabs(gsl_complex z)

  #Complex arithmetic operators
  gsl_complex  gsl_complex_add(gsl_complex a, gsl_complex b)

  gsl_complex  gsl_complex_sub(gsl_complex a, gsl_complex b)

  gsl_complex  gsl_complex_mul(gsl_complex a, gsl_complex b)

  gsl_complex  gsl_complex_div(gsl_complex a, gsl_complex b)

  gsl_complex  gsl_complex_add_real(gsl_complex a, double x)

  gsl_complex  gsl_complex_sub_real(gsl_complex a, double x)

  gsl_complex  gsl_complex_mul_real(gsl_complex a, double x)

  gsl_complex  gsl_complex_div_real(gsl_complex a, double x)

  gsl_complex  gsl_complex_add_imag(gsl_complex a, double y)

  gsl_complex  gsl_complex_sub_imag(gsl_complex a, double y)

  gsl_complex  gsl_complex_mul_imag(gsl_complex a, double y)

  gsl_complex  gsl_complex_div_imag(gsl_complex a, double y)

  gsl_complex  gsl_complex_conjugate(gsl_complex z)

  gsl_complex  gsl_complex_inverse(gsl_complex z)

  gsl_complex  gsl_complex_negative(gsl_complex z)

  #Elementary Complex Functions
  gsl_complex  gsl_complex_sqrt(gsl_complex z)

  gsl_complex  gsl_complex_sqrt_real(double x)

  gsl_complex  gsl_complex_pow(gsl_complex z, gsl_complex a)

  gsl_complex  gsl_complex_pow_real(gsl_complex z, double x)

  gsl_complex  gsl_complex_exp(gsl_complex z)

  gsl_complex  gsl_complex_log(gsl_complex z)

  gsl_complex  gsl_complex_log10(gsl_complex z)

  gsl_complex  gsl_complex_log_b(gsl_complex z, gsl_complex b)

  #Complex Trigonometric Functions
  gsl_complex  gsl_complex_sin(gsl_complex z)

  gsl_complex  gsl_complex_cos(gsl_complex z)

  gsl_complex  gsl_complex_tan(gsl_complex z)

  gsl_complex  gsl_complex_sec(gsl_complex z)

  gsl_complex  gsl_complex_csc(gsl_complex z)

  gsl_complex  gsl_complex_cot(gsl_complex z)

  #Inverse Complex Trigonometric Functions
  gsl_complex  gsl_complex_arcsin(gsl_complex z)

  gsl_complex  gsl_complex_arcsin_real(double z)

  gsl_complex  gsl_complex_arccos(gsl_complex z)

  gsl_complex  gsl_complex_arccos_real(double z)

  gsl_complex  gsl_complex_arctan(gsl_complex z)

  gsl_complex  gsl_complex_arcsec(gsl_complex z)

  gsl_complex  gsl_complex_arcsec_real(double z)

  gsl_complex  gsl_complex_arccsc(gsl_complex z)

  gsl_complex  gsl_complex_arccsc_real(double z)

  gsl_complex  gsl_complex_arccot(gsl_complex z)

  #Complex Hyperbolic Functions
  gsl_complex  gsl_complex_sinh(gsl_complex z)

  gsl_complex  gsl_complex_cosh(gsl_complex z)

  gsl_complex  gsl_complex_tanh(gsl_complex z)

  gsl_complex  gsl_complex_sech(gsl_complex z)

  gsl_complex  gsl_complex_csch(gsl_complex z)

  gsl_complex  gsl_complex_coth(gsl_complex z)

  #Inverse Complex Hyperbolic Functions
  gsl_complex  gsl_complex_arcsinh(gsl_complex z)

  gsl_complex  gsl_complex_arccosh(gsl_complex z)

  gsl_complex  gsl_complex_arccosh_real(double z)

  gsl_complex  gsl_complex_arctanh(gsl_complex z)

  gsl_complex  gsl_complex_arctanh_real(double z)

  gsl_complex  gsl_complex_arcsech(gsl_complex z)

  gsl_complex  gsl_complex_arccsch(gsl_complex z)

  gsl_complex  gsl_complex_arccoth(gsl_complex z)
