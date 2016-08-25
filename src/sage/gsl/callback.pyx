# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from sage.libs.gsl.types cimport gsl_function

cdef extern from "gsl/gsl_diff.h":
  int gsl_diff_central ( gsl_function *f, double x,
                        double *result, double *abserr)


from sage.misc.superseded import deprecation
deprecation(20135, "the module sage.gsl.callback is deprecated")


foo = None
cdef double f(double x, void* p):
    return foo(x)

def diff(g, double x=1.0):
    cdef double result, abserr
    cdef gsl_function F

    global foo
    foo = g
    F.function = f
    F.params = NULL

    gsl_diff_central(&F, x, &result, &abserr)

    return result, abserr
