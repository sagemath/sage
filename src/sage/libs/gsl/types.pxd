# distutils: include_dirs = GSL_INCDIR
from libc.stdio cimport FILE

cdef enum:
    GSL_SUCCESS

from sage.libs.gsl.blas_types cimport *

cdef extern from "gsl/gsl_mode.h":
    ctypedef unsigned int gsl_mode_t
    double GSL_PREC_DOUBLE
    double GSL_PREC_SINGLE
    double GSL_PREC_APPROX


cdef extern from "gsl/gsl_sf_result.h":
    ctypedef struct gsl_sf_result:
        double val
        double err

    ctypedef struct gsl_sf_result_e10:
        double val
        double err
        int    e10


cdef extern from "gsl/gsl_math.h":
    # Definition of an arbitrary function with parameters
    ctypedef struct gsl_function:
        double (* function) (double x, void * params)
        void * params

    # Definition of an arbitrary function returning two values, r1, r2
    ctypedef struct gsl_function_fdf:
        double (* f) (double x, void * params)
        double (* df) (double x, void * params)
        void (* fdf) (double x, void * params, double * f, double * df)
        void * params


cdef extern from "gsl/gsl_complex.h":
    ctypedef double * gsl_complex_packed_array
    ctypedef double * gsl_complex_packed_ptr

    ctypedef struct gsl_complex:
        pass


cdef extern from "gsl/gsl_block_double.h":
    ctypedef struct gsl_block:
        size_t size
        double * data


cdef extern from "gsl/gsl_block_complex_double.h":
    ctypedef struct gsl_block_complex:
        size_t size
        double * data


cdef extern from "gsl/gsl_vector.h":
    ctypedef struct gsl_vector:
        size_t size
        size_t stride
        double *data
        gsl_block *block
        int owner

    ctypedef struct gsl_vector_view:
        gsl_vector vector

    ctypedef struct gsl_vector_const_view:
        gsl_vector vector


cdef extern from "gsl/gsl_vector_complex_double.h":
    ctypedef struct gsl_vector_complex:
        size_t size
        size_t stride
        double *data
        gsl_block_complex *block
        int owner

    ctypedef struct gsl_vector_complex_view:
        gsl_vector_complex vector

    ctypedef struct gsl_vector_complex_const_view:
        gsl_vector_complex vector_complex


cdef extern from "gsl/gsl_matrix_double.h":
    ctypedef struct gsl_matrix:
        size_t size1
        size_t size2
        size_t tda
        double * data
        gsl_block * block
        int owner

    ctypedef struct gsl_matrix_view:
        gsl_matrix matrix

    ctypedef struct gsl_matrix_const_view:
        gsl_matrix matrix


cdef extern from "gsl/gsl_matrix_complex_double.h":
    ctypedef struct gsl_matrix_complex:
        size_t size1
        size_t size2
        size_t tda
        double * data
        gsl_block_complex * block
        int owner

    ctypedef struct gsl_matrix_complex_view:
        gsl_matrix_complex matrix

    ctypedef struct gsl_matrix_complex_const_view:
        gsl_matrix_complex matrix


cdef extern from "gsl/gsl_rng.h":
    ctypedef struct gsl_rng_type
    ctypedef struct gsl_rng


cdef extern from "gsl/gsl_permutation.h":
    ctypedef struct gsl_permutation:
        size_t size
        size_t *data


cdef extern from "gsl/gsl_histogram.h":
    ctypedef struct gsl_histogram
    ctypedef struct gsl_histogram_pdf
