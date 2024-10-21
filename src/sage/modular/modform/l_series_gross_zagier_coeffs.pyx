# sage_setup: distribution = sagemath-schemes
from cysignals.memory cimport check_allocarray, check_calloc, sig_free
from cysignals.signals cimport sig_check, sig_on, sig_off

from sage.rings.fast_arith cimport arith_llong
cdef arith_llong arith = arith_llong()

from sage.arith.misc import kronecker_symbol
from sage.rings.integer_ring import ZZ
from sage.rings.power_series_ring import PowerSeriesRing

from libc.math cimport ceil, floor, sqrt
from libc.string cimport memcpy


cpdef to_series(L, var):
    """
    Create a power series element out of a list ``L`` in the variable`` var``.

    EXAMPLES::

        sage: from sage.modular.modform.l_series_gross_zagier_coeffs import to_series
        sage: to_series([1,10,100], 't')
        1 + 10*t + 100*t^2 + O(t^3)
        sage: to_series([0..5], CDF[['z']].0)
        0.0 + 1.0*z + 2.0*z^2 + 3.0*z^3 + 4.0*z^4 + 5.0*z^5 + O(z^6)
    """
    if var is None:
        return L
    if isinstance(var, str):
        R = PowerSeriesRing(ZZ, var)
    else:
        R = var.parent()
    return R(L).O(len(L))


# TODO, when quadratic form code stabilizes, add this there.
def bqf_theta_series(Q, long bound, var=None):
    r"""
    Return the theta series associated to a positive definite quadratic form.

    For a given form `f = ax^2 + bxy + cy^2` this is the sum

    .. MATH::

        \sum_{(x,y) \in \Z^2} q^{f(x,y)} = \sum_{n=-\infty}^{\infty} r(n)q^n

    where `r(n)` give the number of way `n` is represented by `f`.

    INPUT:

    - ``Q`` -- a positive definite quadratic form
    - ``bound`` -- how many terms to compute
    - ``var`` -- (optional) the variable in which to express this power series

    OUTPUT: a power series in ``var``, or list of ints if ``var`` is unspecified

    EXAMPLES::

        sage: from sage.modular.modform.l_series_gross_zagier_coeffs import bqf_theta_series
        sage: bqf_theta_series([2,1,5], 10)
        [1, 0, 2, 0, 0, 2, 2, 0, 4, 0, 0]
        sage: Q = BinaryQF([41,1,1])
        sage: bqf_theta_series(Q, 50, ZZ[['q']].gen())
        1 + 2*q + 2*q^4 + 2*q^9 + 2*q^16 + 2*q^25 + 2*q^36 + 4*q^41 + 4*q^43 + 4*q^47 + 2*q^49 + O(q^51)
    """
    cdef long a, b, c
    a, b, c = Q
    cdef long* terms = bqf_theta_series_c(NULL, bound, a, b, c)
    L = [terms[i] for i in range(bound + 1)]
    sig_free(terms)
    return to_series(L, var)


cdef long* bqf_theta_series_c(long* terms, long bound, long a, long b, long c) except NULL:
    cdef long x, y, yD
    cdef long xmax, ymin, ymax
    cdef double sqrt_yD

    if a < 0 or 4 * a * c - b * b < 0:
        raise ValueError("Not positive definite.")
    xmax = <long>ceil(2 * sqrt((c * bound) / <double>(4 * a * c - b * b)))
    if terms == NULL:
        terms = <long*>check_calloc(1 + bound, sizeof(long))

    sig_on()
    for x in range(-xmax, xmax + 1):
        yD = b * b * x * x - 4 * c * (a * x * x - bound)
        if yD > 0:
            sqrt_yD = sqrt(yD)
            ymin = <long>ceil((-b * x - sqrt_yD) / (2 * c))
            ymax = <long>floor((-b * x + sqrt_yD) / (2 * c))
            for y in range(ymin, ymax + 1):
                terms[a * x * x + b * x * y + c * y * y] += 1
    sig_off()
    return terms


def gross_zagier_L_series(an_list, Q, long N, long u, var=None):
    """
    Compute the coefficients of the Gross-Zagier `L`-series.

    INPUT:

    - ``an_list`` -- list of coefficients of the `L`-series of an elliptic curve
    - ``Q`` -- a positive definite quadratic form
    - ``N`` -- conductor of the elliptic curve
    - ``u`` -- number of roots of unity in the field associated with ``Q``
    - ``var`` -- (optional) the variable in which to express this power series

    OUTPUT: a power series in ``var``, or list of ints if ``var`` is unspecified

    The number of terms is the length of the input ``an_list``.

    EXAMPLES::

        sage: from sage.modular.modform.l_series_gross_zagier_coeffs import gross_zagier_L_series
        sage: E = EllipticCurve('37a')
        sage: N = 37
        sage: an = E.anlist(60); len(an)
        61
        sage: K.<a> = QuadraticField(-40)
        sage: A = K.class_group().gen(0)
        sage: Q = A.ideal().quadratic_form().reduced_form()
        sage: Q2 = (A**2).ideal().quadratic_form().reduced_form()
        sage: u = K.zeta_order()
        sage: t = PowerSeriesRing(ZZ,'t').gen()
        sage: LA = gross_zagier_L_series(an,Q,N,u,t); LA
        -2*t^2 - 2*t^5 - 2*t^7 - 4*t^13 - 6*t^18 - 4*t^20 + 20*t^22 + 4*t^23
        - 4*t^28 + 8*t^32 - 2*t^37 - 6*t^45 - 18*t^47 + 2*t^50 - 8*t^52
        + 2*t^53 + 20*t^55 + O(t^61)
        sage: len(gross_zagier_L_series(an,Q,N,u))
        61

        sage: LA + gross_zagier_L_series(an,Q2,N,u,t)
        t - 2*t^2 + 2*t^4 - 2*t^5 - 2*t^7 + 3*t^9 + 4*t^10 - 10*t^11 - 4*t^13
        + 4*t^14 - 4*t^16 - 6*t^18 - 4*t^20 + 20*t^22 + 4*t^23 - t^25 + 8*t^26
        - 4*t^28 + 8*t^32 + 4*t^35 + 6*t^36 - 2*t^37 - 18*t^41 - 20*t^44
        - 6*t^45 - 8*t^46 - 18*t^47 - 11*t^49 + 2*t^50 - 8*t^52 + 2*t^53
        + 20*t^55 + 16*t^59 + O(t^61)
    """
    cdef long bound = len(an_list) + 1  # one more term for safety
    cdef long a, b, c
    a, b, c = Q
    cdef long D = b * b - 4 * a * c
    cdef long i, m, n, me, j
    cdef long* con_terms = bqf_theta_series_c(NULL, bound - 1, a, b, c)
    cdef long* terms = NULL
    try:
        terms = <long*>check_allocarray(bound, sizeof(long))
    except MemoryError:
        sig_free(con_terms)
        raise
    i = 0
    for an in an_list:
        sig_check()
        con_terms[i] = con_terms[i] / u * an
        i += 1
    sig_on()
    memcpy(terms, con_terms, sizeof(long) * bound)  # m = 1
    for m in range(2, <long>sqrt(bound) + 1):
        if arith.c_gcd_longlong(D * N, m) == 1:
            me = m * kronecker_symbol(D, m)
            j = 0
            n = 0
            while j < bound:
                terms[j] += me * con_terms[n]
                j += m * m
                n += 1
    sig_off()
    L = [terms[i] for i in range(bound - 1)]
    sig_free(con_terms)
    sig_free(terms)
    return to_series(L, var)
