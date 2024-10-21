# sage_setup: distribution = sagemath-mpmath
"""
Faster versions of some key functions in mpmath.libmp
"""
from sage.libs.mpmath.ext_impl cimport *
from sage.libs.gmp.all cimport *

# the next line is used by mpmath
from sage.libs.mpmath.ext_impl import exp_fixed, cos_sin_fixed, log_int_fixed

# Note: not thread-safe
cdef MPF tmp1
cdef MPF tmp2
MPF_init(&tmp1)
MPF_init(&tmp2)


def mpf_add(tuple x, tuple y, int prec=0, str rnd='d'):
    cdef MPopts opts
    MPF_set_tuple(&tmp1, x)
    MPF_set_tuple(&tmp2, y)
    opts.rounding = rndmode_from_python(rnd)
    opts.prec = prec
    MPF_add(&tmp1, &tmp1, &tmp2, opts)
    return MPF_to_tuple(&tmp1)


def mpf_sub(tuple x, tuple y, int prec=0, str rnd='d'):
    cdef MPopts opts
    MPF_set_tuple(&tmp1, x)
    MPF_set_tuple(&tmp2, y)
    opts.rounding = rndmode_from_python(rnd)
    opts.prec = prec
    MPF_sub(&tmp1, &tmp1, &tmp2, opts)
    return MPF_to_tuple(&tmp1)


def mpf_mul(tuple x, tuple y, int prec=0, str rnd='d'):
    cdef MPopts opts
    MPF_set_tuple(&tmp1, x)
    MPF_set_tuple(&tmp2, y)
    opts.rounding = rndmode_from_python(rnd)
    opts.prec = prec
    MPF_mul(&tmp1, &tmp1, &tmp2, opts)
    return MPF_to_tuple(&tmp1)


def mpf_div(tuple x, tuple y, int prec, str rnd='d'):
    cdef MPopts opts
    MPF_set_tuple(&tmp1, x)
    MPF_set_tuple(&tmp2, y)
    opts.rounding = rndmode_from_python(rnd)
    opts.prec = prec
    MPF_div(&tmp1, &tmp1, &tmp2, opts)
    return MPF_to_tuple(&tmp1)


def mpf_sqrt(tuple x, int prec, str rnd='d'):
    """
    Compute sqrt(x) with mpf value tuples.

    EXAMPLES::

        sage: from sage.libs.mpmath._vendor.mpmath.libmp import mpf_sqrt, from_float, to_float
        sage: x = from_float(2)
        sage: y = mpf_sqrt(x, 53, 'n')
        sage: to_float(y)
        1.4142135623730951
    """
    if x[0]:
        import sage.libs.mpmath._vendor.mpmath.libmp as libmp
        raise libmp.ComplexResult("square root of a negative number")
    cdef MPopts opts
    MPF_set_tuple(&tmp1, x)
    opts.rounding = rndmode_from_python(rnd)
    opts.prec = prec
    MPF_sqrt(&tmp1, &tmp1, opts)
    return MPF_to_tuple(&tmp1)


def mpf_log(tuple x, int prec, str rnd='d'):
    """
    Compute log(x) with mpf value tuples.

    EXAMPLES::

        sage: from sage.libs.mpmath._vendor.mpmath.libmp import mpf_log, from_float, to_float
        sage: x = from_float(2)
        sage: y = mpf_log(x, 53, 'n')
        sage: to_float(y)
        0.6931471805599453
    """
    if x[0]:
        import sage.libs.mpmath._vendor.mpmath.libmp as libmp
        raise libmp.ComplexResult("logarithm of a negative number")
    cdef MPopts opts
    MPF_set_tuple(&tmp1, x)
    opts.rounding = rndmode_from_python(rnd)
    opts.prec = prec
    MPF_log(&tmp1, &tmp1, opts)
    return MPF_to_tuple(&tmp1)


def mpf_exp(tuple x, int prec, str rnd='d'):
    """
    Compute exp(x) with mpf value tuples.

    EXAMPLES::

        sage: from sage.libs.mpmath._vendor.mpmath.libmp import mpf_exp, from_float, to_float
        sage: x = from_float(2)
        sage: z = mpf_exp(x, 53, 'n')
        sage: to_float(z)
        7.38905609893065
    """
    cdef MPopts opts
    MPF_set_tuple(&tmp1, x)
    opts.rounding = rndmode_from_python(rnd)
    opts.prec = prec
    MPF_exp(&tmp1, &tmp1, opts)
    return MPF_to_tuple(&tmp1)


def mpf_cos(tuple x, int prec, str rnd='d'):
    """
    Compute cos(x) with mpf value tuples.

    EXAMPLES::

        sage: from sage.libs.mpmath._vendor.mpmath.libmp import mpf_cos, from_float, to_float
        sage: x = from_float(1)
        sage: y = mpf_cos(x, 53, 'n')
        sage: to_float(y)
        0.5403023058681398
    """
    cdef MPopts opts
    MPF_set_tuple(&tmp1, x)
    opts.rounding = rndmode_from_python(rnd)
    opts.prec = prec
    MPF_cos(&tmp1, &tmp1, opts)
    return MPF_to_tuple(&tmp1)


def mpf_sin(tuple x, int prec, str rnd='d'):
    """
    Compute sin(x) with mpf value tuples.

    EXAMPLES::

        sage: from sage.libs.mpmath._vendor.mpmath.libmp import mpf_sin, from_float, to_float
        sage: x = from_float(1)
        sage: y = mpf_sin(x, 53, 'n')
        sage: to_float(y)
        0.8414709848078965
    """
    cdef MPopts opts
    MPF_set_tuple(&tmp1, x)
    opts.rounding = rndmode_from_python(rnd)
    opts.prec = prec
    MPF_sin(&tmp1, &tmp1, opts)
    return MPF_to_tuple(&tmp1)


def mpc_sqrt(tuple z, int prec, str rnd='d'):
    """
    Compute sqrt(z) with mpc value tuples.

    EXAMPLES::

        sage: from sage.libs.mpmath._vendor.mpmath.libmp import mpc_sqrt, from_float, to_float
        sage: z = from_float(-2), from_float(0)
        sage: re, im = mpc_sqrt(z, 53, 'n')
        sage: to_float(re), to_float(im)
        (0.0, 1.4142135623730951)
    """
    cdef tuple a, b
    cdef MPopts opts
    a, b = z
    MPF_set_tuple(&tmp1, a)
    MPF_set_tuple(&tmp2, b)
    opts.rounding = rndmode_from_python(rnd)
    opts.prec = prec
    MPF_complex_sqrt(&tmp1, &tmp2, &tmp1, &tmp2, opts)
    return MPF_to_tuple(&tmp1), MPF_to_tuple(&tmp2)


def mpc_exp(tuple z, int prec, str rnd='d'):
    """
    Compute exp(z) with mpc value tuples.

    EXAMPLES::

        sage: from sage.libs.mpmath._vendor.mpmath.libmp import mpc_exp, from_float, to_float
        sage: z = from_float(0), from_float(1)
        sage: re, im = mpc_exp(z, 53, 'n')
        sage: to_float(re), to_float(im)
        (0.5403023058681398, 0.8414709848078965)
    """
    cdef tuple a, b
    cdef MPopts opts
    a, b = z
    MPF_set_tuple(&tmp1, a)
    MPF_set_tuple(&tmp2, b)
    opts.rounding = rndmode_from_python(rnd)
    opts.prec = prec
    MPF_complex_exp(&tmp1, &tmp2, &tmp1, &tmp2, opts)
    return MPF_to_tuple(&tmp1), MPF_to_tuple(&tmp2)


def mpf_pow(tuple x, tuple y, int prec, str rnd='d'):
    """
    Compute x ^ y with mpf value tuples.

    EXAMPLES::

        sage: from sage.libs.mpmath._vendor.mpmath.libmp import mpf_pow, from_float, to_float
        sage: x = from_float(2)
        sage: y = from_float(3)
        sage: z = mpf_pow(x, y, 53, 'n')
        sage: to_float(z)
        8.0
    """
    cdef MPopts opts
    MPF_set_tuple(&tmp1, x)
    MPF_set_tuple(&tmp2, y)
    opts.rounding = rndmode_from_python(rnd)
    opts.prec = prec
    if MPF_pow(&tmp1, &tmp1, &tmp2, opts):
        import sage.libs.mpmath._vendor.mpmath.libmp as libmp
        raise libmp.ComplexResult("negative number raised to a fractional power")
    return MPF_to_tuple(&tmp1)
