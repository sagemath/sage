# cython: boundscheck=False, wraparound=False
# distutils: libraries = sirocco
# distutils: language = c++
# sage_setup: distribution = sagemath-sirocco

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
    mpfr_t* homotopyPath_mp_comps(int degree, mpfr_t *_coef, mpfr_t _y0R, mpfr_t _y0I, int prec, int nothercomps, int *degreescomps, mpfr_t *_coefscomps)
    double* homotopyPath_comps(int degree, double *_coef, double _y0R, double _y0I, int nothercomps, int *degreescomps, double *_coefscomps)


cpdef list[list] contpath_mp(int deg, list values, RealNumber y0r, RealNumber y0i, int prec) noexcept:
    """
    Mimics :func:`contpath`, but with the following differences:

    - The floating point numbers can be arbitrary precision RealNumbers.

    - A extra argument is needed, indicating the bits of precision used
      in the computations.

    EXAMPLES::

        sage: from sage.libs.sirocco import contpath_mp   # optional - sirocco
        sage: from sage.rings.real_mpfr import RR
        sage: pol = list(map(RR, [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]))
        sage: contpath_mp(2, pol, RR(0), RR(0), 53)   # optional - sirocco # abs tol 1e-15
        [(0.000000000000000, 0.000000000000000, 0.000000000000000),
         (0.500000000000000, -0.250000000000000, 0.000000000000000),
         (1.00000000000000, -1.00000000000000, 0.000000000000000)]
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

cpdef list[list] contpath_mp_comps(int deg, list values, RealNumber y0r, RealNumber y0i, int prec, list otherdegs, list othercoefs) noexcept:
    """
    Mimics :func:`contpath`, but with the following differences:

    - The floating point numbers can be arbitrary precision RealNumbers.

    - A extra argument is needed, indicating the bits of precision used
      in the computations.

    EXAMPLES::

        sage: from sage.libs.sirocco import contpath_mp_comps   # optional - sirocco
        sage: from sage.rings.real_mpfr import RR
        sage: pol = list(map(RR,[0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]))
        sage: fac = list(map(RR,[0, 0, 0.1, 0.2, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]))
        sage: contpath_mp_comps(2, pol, RR(0), RR(0), 53, [2], fac)   # optional - sirocco # abs tol 1e-15
        [(0.000000000000000, 0.000000000000000, 0.000000000000000),
         (0.125000000000000, -0.0156250000000000, 0.000000000000000),
         (0.250000000000000, -0.0625000000000000, 0.000000000000000),
         (0.375000000000000, -0.140625000000000, 0.000000000000000),
         (0.500000000000000, -0.250000000000000, 0.000000000000000),
         (0.625000000000000, -0.390625000000000, 0.000000000000000),
         (0.750000000000000, -0.562500000000000, 0.000000000000000),
         (0.875000000000000, -0.765625000000000, 0.000000000000000),
         (1.00000000000000, -1.00000000000000, 0.000000000000000)]
    """

    cdef mpfr_t* cvalues = <mpfr_t*> check_allocarray(len(values), sizeof(mpfr_t))
    cdef mpfr_t* cothercoefs = <mpfr_t*> check_allocarray(len(othercoefs), sizeof(mpfr_t))
    cdef int* cotherdegs = <int*> check_allocarray(len(otherdegs), sizeof(int))
    cdef mpfr_t* rop
    cdef int i, j
    cdef mpfr_t y0R
    cdef mpfr_t y0I

    for j in range(len(values)):
        mpfr_init2(cvalues[j], prec)
        mpfr_set(cvalues[j], (<RealNumber>values[j]).value, MPFR_RNDN)

    for j in range(len(othercoefs)):
        mpfr_init2(cothercoefs[j], prec)
        mpfr_set(cothercoefs[j], (<RealNumber>othercoefs[j]).value, MPFR_RNDN)

    for j in range(len(otherdegs)):
        cotherdegs[j] = int(otherdegs[j])
    sig_on()
    mpfr_init2(y0R, prec)
    mpfr_set(y0R, (<RealNumber>y0r).value, MPFR_RNDN)
    mpfr_init2(y0I, prec)
    mpfr_set(y0I, (<RealNumber>y0i).value, MPFR_RNDN)
    rop = homotopyPath_mp_comps(deg, cvalues, y0R, y0I, prec, int(len(otherdegs)), cotherdegs, cothercoefs)
    sig_off()
    for j in range(len(values)):
        mpfr_clear(cvalues[j])
    free(cvalues)
    for j in range(len(othercoefs)):
        mpfr_clear(cothercoefs[j])
    free(cothercoefs)
    free(cotherdegs)
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


cpdef list[list] contpath(int deg, list values, double y0r, double y0i) noexcept:
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

    EXAMPLES::

        sage: from sage.libs.sirocco import contpath   # optional - sirocco
        sage: from sage.rings.real_mpfr import RR
        sage: pol = list(map(RR,[0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]))
        sage: contpath(2, pol, RR(0), RR(0))  # optional - sirocco # abs tol 1e-15
        [(0.0, 0.0, 0.0),
         (0.3535533905932738, -0.12500000000000003, 0.0),
         (0.7071067811865476, -0.5000000000000001, 0.0),
         (1.0, -1.0, 0.0)]
    """
    cdef double* rop
    cdef double* c_values = <double*> check_allocarray(len(values), sizeof(double))
    cdef int clen = <int> len(values)
    cdef int i
    for i, v in enumerate(values):
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

cpdef list[list] contpath_comps(int deg, list values, double y0r, double y0i, list otherdegrees, list othercoefs) noexcept:
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

    EXAMPLES::

        sage: from sage.libs.sirocco import contpath_comps   # optional - sirocco
        sage: from sage.rings.real_mpfr import RR
        sage: pol = list(map(RR,[0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]))
        sage: fac = list(map(RR,[0, 0, 0.1, 0.2, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]))
        sage: contpath_comps(2, pol, RR(0), RR(0), [2], fac)   # optional - sirocco # abs tol 1e-15
        [(0.0, 0.0, 0.0),
         (0.125, -0.015625, 0.0),
         (0.25, -0.0625, 0.0),
         (0.375, -0.140625, 0.0),
         (0.5, -0.25, 0.0),
         (0.625, -0.390625, 0.0),
         (0.75, -0.5625, 0.0),
         (0.875, -0.765625, 0.0),
         (1.0, -1.0, 0.0)]
    """
    cdef double* rop
    cdef double* c_values = <double*> check_allocarray(len(values), sizeof(double))
    cdef int* c_otherdegrees = <int*> check_allocarray(len(otherdegrees), sizeof(int))
    cdef double* c_othercoefs = <double*> check_allocarray(len(othercoefs), sizeof(double))
    cdef int clen = <int> len(values)
    cdef int i
    for i, v in enumerate(values):
        c_values[i] = values[i]

    for i, v in enumerate(otherdegrees):
        c_otherdegrees[i] = otherdegrees[i]
    for i, v in enumerate(othercoefs):
        c_othercoefs[i] = othercoefs[i]

    cdef double y0R = y0r
    cdef double y0I = y0i
    sig_on()
    rop = homotopyPath_comps(deg, c_values, y0R, y0I, int(len(otherdegrees)), c_otherdegrees, c_othercoefs)
    sig_off()
    if rop == NULL:
        raise ValueError("libsirocco could not guarantee one step")
    cdef int n = int(rop[0])
    cdef list l = [0] * n
    for i in range(n):
        l[i] = (rop[3*i+1], rop[3*i+2], rop[3*i+3])
    free(rop)
    free(c_values)
    free(c_otherdegrees)
    free(c_othercoefs)
    return l
