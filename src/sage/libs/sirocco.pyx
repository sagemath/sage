#cython: boundscheck=False, wraparound=False
# distutils: libraries = sirocco
# distutils: language = c++
# sage_setup: distribution = sage-sirocco

r"""
Cython wrapper for sirocco library

This is used to call the sirocco library directly from Python.

AUTHORS:

- Miguel Marco (2016-07-19): initial version.
"""

from cysignals.signals cimport sig_on, sig_off
from cysignals.memory cimport check_allocarray, sig_free as free

from sage.libs.mpfr cimport *
from sage.rings.real_mpfr cimport RealNumber
from sage.rings.real_mpfr import RealField

cdef extern from "sirocco.h":
    mpfr_t* homotopyPath_mp(int degree, mpfr_t *_coef, mpfr_t _y0R, mpfr_t _y0I, int prec)
    double* homotopyPath(int degree, double *_coef, double _y0R, double _y0I)



cpdef list[list] contpath_mp(int deg, list values, RealNumber y0r, RealNumber y0i, int prec):
    """
    Mimics :func:`contpath`, but with the following differences:

    - The floating point numbers can be arbitrary precision RealNumbers.

    - A extra argument is needed, indicating the bits of precision used
      in the computations.
    """
    cdef mpfr_t* cvalues = <mpfr_t*> check_allocarray(len(values), sizeof(mpfr_t))
    cdef mpfr_t* rop
    cdef int i, j
    cdef mpfr_t y0R
    cdef mpfr_t y0I

    for j in range(len(values)):
        mpfr_init2(cvalues[j], prec)
        mpfr_set(cvalues[j], (<RealNumber>values[j]).value, MPFR_RNDN)

    sig_on()
    mpfr_init2(y0R, prec)
    mpfr_set(y0R, (<RealNumber>y0r).value, MPFR_RNDN)
    mpfr_init2(y0I, prec)
    mpfr_set(y0I, (<RealNumber>y0i).value, MPFR_RNDN)
    rop = homotopyPath_mp(deg, cvalues, y0R, y0I, prec)
    sig_off()

    for j in range(len(values)):
        mpfr_clear(cvalues[j])
    free(cvalues)

    if rop == NULL:
        raise ValueError("libsirocco could not guarantee one step")

    cdef int n = mpfr_get_si(rop[0], MPFR_RNDN)
    cdef list l = []
    cdef list inner
    cdef RealNumber RN
    field = RealField(prec)
    for i in range(n):
        inner = []
        for j in range(3*i+1, 3*(i+1)+1):
            RN = <RealNumber> RealNumber.__new__(RealNumber, field)
            mpfr_set(RN.value, rop[j], MPFR_RNDN)
            mpfr_clear(rop[j])
            inner.append(RN)
        l.append(tuple(inner))
    free(rop)
    return l

cpdef list[list] contpath(int deg, list values, double y0r, double y0i):
    """
    INPUT:

    - An integer, representing the degree of the polynomial

    - A list of floating point numbers. Each four consecutive elements
      of this list represent the interval corresponding to a coefficient.
      Coefficients are listed in increasing deglex order, and inside each
      coefficients, the four numbers represent the lower real, upper real,
      lower imaginary and real imaginary limits of the interval.

    - A float representing the real part of the initial root approximation

    - A float representing the imaginary part of the initial root.

    OUTPUT:

    A list of tuples. Each tuple represents the `x` value (between 0 and 1)
    and the real and imaginary parts of the `y` value of a vertex in
    the piecewise linear approximation of the path tracked by the root.
    """
    cdef double* rop
    cdef double* c_values = <double*> check_allocarray(len(values), sizeof(double))
    cdef int clen = <int> len(values)
    cdef int i
    for i,v in enumerate(values):
        c_values[i] = values[i]
    cdef double y0R = y0r
    cdef double y0I = y0i
    sig_on()
    rop = homotopyPath(deg, c_values, y0R, y0I)
    sig_off()
    if rop == NULL:
        raise ValueError("libsirocco could not guarantee one step")
    cdef int n = int(rop[0])
    cdef list l = [0] * n
    for i in range(n):
        l[i] = (rop[3*i+1], rop[3*i+2], rop[3*i+3])
    free(rop)
    free(c_values)
    return l

