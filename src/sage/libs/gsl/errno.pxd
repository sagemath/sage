# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = GSL_LIBDIR
# distutils: include_dirs = GSL_INCDIR
from libc.stdio cimport FILE

cdef extern from "gsl/gsl_errno.h":

    ctypedef void gsl_error_handler_t (char * reason, char * file,int line, int gsl_errno)

    ctypedef void gsl_stream_handler_t (char * label, char * file,int line, char * reason)


    void gsl_error (char * reason,char * file, int line, int gsl_errno)

    void gsl_stream_printf (char *label, char *file,int line, char *reason)

    char * gsl_strerror (int gsl_errno)

    gsl_error_handler_t * gsl_set_error_handler (gsl_error_handler_t * new_handler)

    gsl_error_handler_t * gsl_set_error_handler_off()

    gsl_stream_handler_t * gsl_set_stream_handler (gsl_stream_handler_t * new_handler)

    FILE * gsl_set_stream (FILE * new_stream)
