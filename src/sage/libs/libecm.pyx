# distutils: libraries = ecm
# distutils: extra_link_args = LINUX_NOEXECSTACK
r"""
The Elliptic Curve Method for Integer Factorization (ECM)

Sage includes GMP-ECM, which is a highly optimized implementation of
Lenstra's elliptic curve factorization method.
See https://gitlab.inria.fr/zimmerma/ecm for more about GMP-ECM.
This file provides a Cython interface to the GMP-ECM library.

AUTHORS:

- Robert L Miller (2008-01-21): library interface (clone of ecmfactor.c)

- Jeroen Demeyer (2012-03-29): signal handling, documentation
- Paul Zimmermann (2011-05-22) -- added input/output of sigma

EXAMPLES::

    sage: from sage.libs.libecm import ecmfactor
    sage: result = ecmfactor(999, 0.00)
    sage: result[0]
    True
    sage: result[1] in [3, 9, 27, 37, 111, 333, 999] or result[1]
    True
    sage: result = ecmfactor(999, 0.00, verbose=True)
    Performing one curve with B1=0
    Found factor in step 1: ...
    sage: result[0]
    True
    sage: result[1] in [3, 9, 27, 37, 111, 333, 999] or result[1]
    True
    sage: ecmfactor(2^128+1,1000,sigma=227140902)
    (True, 5704689200685129054721, 227140902)
"""

#*****************************************************************************
#       Copyright (C) 2008 Robert Miller
#       Copyright (C) 2012 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cysignals.signals cimport sig_on, sig_off

from sage.libs.gmp.mpz cimport *
from sage.rings.integer cimport Integer

cdef extern from "ecm.h":
    ctypedef struct __ecm_param_struct:
        int method
        mpz_t x
        mpz_t sigma
    ctypedef __ecm_param_struct ecm_params[1]
    int ecm_factor (mpz_t, mpz_t, double, ecm_params)
    void ecm_init (ecm_params)
    void ecm_clear (ecm_params)
    int ECM_NO_FACTOR_FOUND


def ecmfactor(number, double B1, verbose=False, sigma=0):
    """
    Try to find a factor of a positive integer using ECM (Elliptic Curve Method).
    This function tries one elliptic curve.

    INPUT:

    - ``number`` -- positive integer to be factored

    - ``B1`` -- bound for step 1 of ECM

    - ``verbose`` -- boolean (default: ``False``); print some debugging
      information

    OUTPUT:

    Either ``(False, None)`` if no factor was found, or ``(True, f)``
    if the factor `f` was found.

    EXAMPLES::

        sage: from sage.libs.libecm import ecmfactor

    This number has a small factor which is easy to find for ECM::

        sage: N = 2^167 - 1
        sage: factor(N)
        2349023 * 79638304766856507377778616296087448490695649
        sage: ecmfactor(N, 2e5)
        (True, 2349023, ...)

    If a factor was found, we can reproduce the factorization with the same
    sigma value::

        sage: N = 2^167 - 1
        sage: ecmfactor(N, 2e5, sigma=1473308225)
        (True, 2349023, 1473308225)

    With a smaller B1 bound, we may or may not succeed::

        sage: ecmfactor(N, 1e2)  # random
        (False, None)

    The following number is a Mersenne prime, so we don't expect to
    find any factors (there is an extremely small chance that we get
    the input number back as factorization)::

        sage: N = 2^127 - 1
        sage: N.is_prime()
        True
        sage: ecmfactor(N, 1e3)
        (False, None)

    If we have several small prime factors, it is possible to find a
    product of primes as factor::

        sage: N = 2^179 - 1
        sage: factor(N)
        359 * 1433 * 1489459109360039866456940197095433721664951999121
        sage: ecmfactor(N, 1e3)  # random
        (True, 514447, 3475102204)

    We can ask for verbose output::

        sage: N = 12^97 - 1
        sage: factor(N)
        11 * 43570062353753446053455610056679740005056966111842089407838902783209959981593077811330507328327968191581
        sage: ecmfactor(N, 100, verbose=True)
        Performing one curve with B1=100
        Found factor in step 1: 11
        (True, 11, ...)
        sage: ecmfactor(N/11, 100, verbose=True)
        Performing one curve with B1=100
        Found no factor.
        (False, None)

    TESTS:

    Check that ``ecmfactor`` can be interrupted (factoring a large
    prime number)::

        sage: alarm(0.5); ecmfactor(2^521-1, 1e7)
        Traceback (most recent call last):
        ...
        AlarmInterrupt

    Some special cases::

        sage: ecmfactor(1, 100)
        Traceback (most recent call last):
        ...
        ValueError: Input number (1) must be greater than 1
        sage: ecmfactor(0, 100)
        Traceback (most recent call last):
        ...
        ValueError: Input number (0) must be greater than 1
    """
    cdef mpz_t n, f
    cdef int res
    cdef Integer sage_int_f, sage_int_number, sage_int_sigma
    cdef ecm_params q

    sage_int_f = Integer(0)
    sage_int_number = Integer(number)
    sage_int_sigma = Integer(sigma)

    if number <= 1:
        raise ValueError("Input number (%s) must be greater than 1"%number)

    if verbose:
        print("Performing one curve with B1=%1.0f" % B1)

    sig_on()
    mpz_init(n)
    mpz_set(n, sage_int_number.value)
    mpz_init(f) # For potential factor
    ecm_init(q)
    mpz_set(q.sigma,sage_int_sigma.value)

    res = ecm_factor(f, n, B1, q)

    if res > 0:
        mpz_set(sage_int_f.value, f)
        mpz_set(sage_int_sigma.value, q.sigma)

    mpz_clear(f)
    mpz_clear(n)
    ecm_clear(q)
    sig_off()

    if res > 0:
        if verbose:
            print("Found factor in step %d: %d" % (res,sage_int_f))
        return (True, sage_int_f, sage_int_sigma)
    elif res == ECM_NO_FACTOR_FOUND:
        if verbose:
            print("Found no factor.")
        return (False, None)
    else:
        raise RuntimeError( "ECM lib error" )
