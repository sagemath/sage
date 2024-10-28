"""
Convert Sage/Python objects to real/complex intervals
"""
#*****************************************************************************
#       Copyright (C) 2018 Jeroen Demeyer <J.Demeyer@UGent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cpython.float cimport PyFloat_AS_DOUBLE
from cpython.complex cimport PyComplex_RealAsDouble, PyComplex_ImagAsDouble

from sage.libs.mpfr cimport *
from sage.libs.mpfi cimport *
from sage.libs.gsl.complex cimport *

from sage.arith.long cimport integer_check_long
from sage.cpython.string cimport bytes_to_str
from sage.structure.element cimport Element

import sage.rings.abc
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from sage.rings.real_mpfi cimport RealIntervalFieldElement, RealIntervalField_class
from sage.rings.real_mpfr cimport RealNumber
from sage.rings.real_double cimport RealDoubleElement
from sage.rings.complex_mpfr cimport ComplexNumber
from sage.rings.complex_interval cimport ComplexIntervalFieldElement
from sage.rings.complex_double cimport ComplexDoubleElement

from cypari2.gen cimport Gen


cdef inline int return_real(mpfi_ptr im) noexcept:
    """
    Called by ``mpfi_set_sage`` on the imaginary part when converting
    a real number.
    """
    if im is not NULL:
        mpfi_set_ui(im, 0)
    return 0


cdef int mpfi_set_sage(mpfi_ptr re, mpfi_ptr im, x, field, int base) except -1:
    """
    Convert any object ``x`` to an MPFI interval or a pair of
    real and complex intervals.

    INPUT:

    - ``re`` -- a pre-initialized MPFI interval

    - ``im`` -- a pre-initialized MPFI interval or NULL

    - ``x`` -- any Sage or Python object to be converted to an interval

    - ``field`` -- a ``RealIntervalField`` or ``ComplexIntervalField``
      of the right precision (real or complex doesn't matter)

    - ``base`` -- base to use for string conversion

    OUTPUT:

    - if conversion is possible: set ``re`` and ``im`` (if applicable)
      and return 0.

    - if ``x`` is complex but ``im`` is ``NULL``: convert only if the
      imaginary component is 0.

    - in all other cases: raise an exception.
    """
    cdef RealIntervalFieldElement ri
    cdef ComplexIntervalFieldElement zi
    cdef ComplexNumber zn
    cdef ComplexDoubleElement zd
    cdef bytes s

    if im is not NULL and isinstance(x, tuple):
        # For complex numbers, interpret tuples as real/imag parts
        if len(x) != 2:
            raise TypeError("tuple defining a complex number must have length 2")
        mpfi_set_sage(re, NULL, x[0], field, base)
        mpfi_set_sage(im, NULL, x[1], field, base)
        return 0
    if isinstance(x, (list, tuple)):
        # Interpret entries in x as endpoints of interval
        if len(x) != 2:
            raise TypeError("list defining an interval must have length 2")
        return mpfi_interv_sage(re, im, x[0], x[1], field, base)

    cdef long value
    cdef int err

    # Check for known types. First check for Element to reduce the
    # number of checks below.
    if isinstance(x, Element):
        # Real
        if isinstance(x, RealIntervalFieldElement):
            mpfi_set(re, (<RealIntervalFieldElement>x).value)
            return return_real(im)
        if isinstance(x, RealNumber):
            mpfi_set_fr(re, (<RealNumber>x).value)
            return return_real(im)
        if isinstance(x, Rational):
            mpfi_set_q(re, (<Rational>x).value)
            return return_real(im)
        if isinstance(x, Integer):
            mpfi_set_z(re, (<Integer>x).value)
            return return_real(im)
        if isinstance(x, RealDoubleElement):
            mpfi_set_d(re, (<RealDoubleElement>x)._value)
            return return_real(im)

        # Complex
        if isinstance(x, ComplexIntervalFieldElement):
            zi = <ComplexIntervalFieldElement>x
            if im is NULL:
                if not mpfi_is_zero(zi.__im):
                    raise TypeError(f"unable to convert complex interval {x!r} to real interval")
            else:
                mpfi_set(im, zi.__im)
            mpfi_set(re, zi.__re)
            return 0
        if isinstance(x, ComplexNumber):
            zn = <ComplexNumber>x
            if im is NULL:
                if mpfr_sgn(zn.__im) != 0:
                    raise TypeError(f"unable to convert complex number {x!r} to real interval")
            else:
                mpfi_set_fr(im, zn.__im)
            mpfi_set_fr(re, zn.__re)
            return 0
        if isinstance(x, ComplexDoubleElement):
            zd = <ComplexDoubleElement>x
            if im is NULL:
                if GSL_IMAG(zd._complex):
                    raise TypeError(f"unable to convert complex number {x!r} to real interval")
            else:
                mpfi_set_d(im, GSL_IMAG(zd._complex))
            mpfi_set_d(re, GSL_REAL(zd._complex))
            return 0
    else:  # not a Sage Element
        # Real
        if isinstance(x, float):
            mpfi_set_d(re, PyFloat_AS_DOUBLE(x))
            return return_real(im)
        if integer_check_long(x, &value, &err):
            if err == 0:
                mpfi_set_si(re, value)
            else:
                mpfi_set_via_RR(re, x, field)
            return return_real(im)
        if isinstance(x, unicode):
            x = x.encode("ascii")
        if isinstance(x, bytes):
            s = (<bytes>x).replace(b'..', b',').replace(b' ', b'').replace(b'+infinity', b'@inf@').replace(b'-infinity', b'-@inf@')
            if mpfi_set_str(re, s, base):
                x = bytes_to_str(x)
                raise TypeError(f"unable to convert {x!r} to real interval")
            return return_real(im)

        # Complex
        if isinstance(x, Gen):
            imag = x.imag()
            if im is NULL:
                if imag:
                    raise TypeError(f"unable to convert complex PARI/GP element {x!r} to real interval")
            else:
                mpfi_set_via_RR(im, imag, field)
            mpfi_set_via_RR(re, x.real(), field)
            return 0
        if isinstance(x, complex):
            imag = PyComplex_ImagAsDouble(x)
            if im is NULL:
                if imag:
                    raise TypeError(f"unable to convert complex number {x!r} to real interval")
            else:
                mpfi_set_d(im, imag)
            mpfi_set_d(re, PyComplex_RealAsDouble(x))
            return 0

    # No known type, try _real_mpfi_ or _complex_mpfi_ methods
    if im is not NULL:
        try:
            m = x._complex_mpfi_
        except AttributeError:
            pass
        else:
            if not isinstance(field, sage.rings.abc.ComplexIntervalField):
                field = field.complex_field()
            e = <ComplexIntervalFieldElement?>m(field)
            mpfi_swap(re, e.__re)
            mpfi_swap(im, e.__im)
            return 0

    try:
        m = x._real_mpfi_
    except AttributeError:
        pass
    else:
        if not isinstance(field, RealIntervalField_class):
            field = field.real_field()
        ri = <RealIntervalFieldElement?>m(field)
        mpfi_swap(re, ri.value)
        return return_real(im)

    # Finally, try converting via the corresponding RealField
    mpfi_set_via_RR(re, x, field)
    return return_real(im)


cdef int mpfi_interv_sage(mpfi_ptr re, mpfi_ptr im, x, y, field, int base) except -1:
    """
    Like ``mpfi_set_sage`` but construct the interval around ``x`` and
    ``y``. It is not required that ``x <= y`` or that ``x`` and ``y``
    are of the same type.

    INPUT: see ``mpfi_set_sage``
    """
    cdef long valx, valy
    cdef int err
    if type(x) is type(y):
        if isinstance(x, RealNumber):
            mpfi_interv_fr(re, (<RealNumber>x).value, (<RealNumber>y).value)
            return return_real(im)
        if isinstance(x, RealDoubleElement):
            mpfi_interv_d(re, (<RealDoubleElement>x)._value, (<RealDoubleElement>y)._value)
            return return_real(im)
        if isinstance(x, Integer):
            mpfi_interv_z(re, (<Integer>x).value, (<Integer>y).value)
            return return_real(im)
        if isinstance(x, Rational):
            mpfi_interv_q(re, (<Rational>x).value, (<Rational>y).value)
            return return_real(im)
        if isinstance(x, float):
            mpfi_interv_d(re, PyFloat_AS_DOUBLE(x), PyFloat_AS_DOUBLE(y))
            return return_real(im)
        # General check for C long
        integer_check_long(x, &valx, &err)
        if err == 0:
            integer_check_long(y, &valy, &err)
            if err == 0:
                mpfi_interv_si(re, valx, valy)
                return return_real(im)

    # General case: convert both x and y to an interval and take the
    # union

    # First handle x
    mpfi_set_sage(re, im, x, field, base)

    # Now handle y, which requires temporary mpfi variables.
    cdef mpfi_t tmp1, tmp2
    cdef mpfr_prec_t prec = mpfi_get_prec(re)

    mpfi_init2(tmp1, prec)
    cdef mpfi_ptr tmpim = NULL
    if im is not NULL:
        mpfi_init2(tmp2, prec)
        tmpim = tmp2

    try:
        mpfi_set_sage(tmp1, tmpim, y, field, base)
        mpfi_union(re, re, tmp1)
        if im is not NULL:
            mpfi_union(im, im, tmpim)
    finally:
        mpfi_clear(tmp1)
        if tmpim is not NULL:
            mpfi_clear(tmpim)


cdef int mpfi_set_via_RR(mpfi_ptr re, x, field) except -1:
    """
    Convert ``x`` to an MPFI interval by converting ``x`` to the
    appropriate real fields.

    INPUT: see ``mpfi_set_sage``
    """
    cdef RealIntervalField_class RIF
    if isinstance(field, RealIntervalField_class):
        RIF = <RealIntervalField_class>field
    else:
        RIF = <RealIntervalField_class?>field.real_field()

    try:
        ra = RIF.__lower_field(x)
        rb = RIF.__upper_field(x)
    except TypeError:
        raise TypeError(f"unable to convert {x!r} to real interval")
    mpfi_interv_fr(re, (<RealNumber>ra).value, (<RealNumber>rb).value)
