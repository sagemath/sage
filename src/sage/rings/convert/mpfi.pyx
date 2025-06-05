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

import re

from cpython.float cimport PyFloat_AS_DOUBLE
from cpython.complex cimport PyComplex_RealAsDouble, PyComplex_ImagAsDouble

from libc.stdio cimport printf

from sage.libs.mpfr cimport *
from sage.libs.mpfi cimport *
from sage.libs.gmp.mpz cimport *
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


NUMBER = re.compile(rb'([+-]?(0[XxBb])?[0-9A-Za-z]+)\.([0-9A-Za-z]*)\?([0-9]*)(?:([EePp@])([+-]?[0-9]+))?')
# example: -0xABC.DEF?12@5
# match groups: (-0xABC) (0x) (DEF) (12) (@) (5)

cdef int _from_str_question_style(mpfi_ptr x, bytes s, int base) except -1:
    """
    Convert a string in question style to an MPFI interval.

    INPUT:

    - ``x`` -- a pre-initialized MPFI interval

    - ``s`` -- the string to convert

    - ``base`` -- base to use for string conversion

    OUTPUT:

    - if conversion is possible: set ``x`` and return 0.

    - in all other cases: return some nonzero value, or raise an exception.

    TESTS:

    Double check that ``ZZ``, ``RR`` and ``RIF`` follows the string
    conversion rule for base different from `10` (except ``ZZ``
    which only allows base up to `36`)::

        sage: ZZ("0x123", base=0)
        291
        sage: RR("0x123.e1", base=0)  # rel tol 1e-12
        291.878906250000
        sage: RR("0x123.@1", base=0)  # rel tol 1e-12
        4656.00000000000
        sage: RIF("0x123.4@1", base=0)
        4660
        sage: ZZ("1Xx", base=36)  # case insensitive
        2517
        sage: ZZ("1Xx", base=62)
        Traceback (most recent call last):
        ...
        ValueError: base (=62) must be 0 or between 2 and 36
        sage: RR("1Xx", base=36)  # rel tol 1e-12
        2517.00000000000
        sage: RR("0x123", base=36)  # rel tol 1e-12
        1.54101900000000e6
        sage: RR("-1Xx@-1", base=62)  # rel tol 1e-12
        -95.9516129032258
        sage: RIF("1Xx@-1", base=62)  # rel tol 1e-12
        95.95161290322580?
        sage: RIF("1aE1", base=11)
        Traceback (most recent call last):
        ...
        TypeError: unable to convert '1aE1' to real interval
        sage: RIF("1aE1", base=11)
        Traceback (most recent call last):
        ...
        TypeError: unable to convert '1aE1' to real interval

    General checks::

        sage: RIF("123456.?2").endpoints()  # rel tol 1e-12
        (123454.0, 123458.0)
        sage: RIF("1234.56?2").endpoints()  # rel tol 1e-12
        (1234.54, 1234.58)
        sage: RIF("1234.56?2e2").endpoints()  # rel tol 1e-12
        (123454.0, 123458.0)
        sage: x = RIF("-1234.56?2e2"); x.endpoints()  # rel tol 1e-12
        (-123458.0, -123454.0)
        sage: x
        -1.2346?e5
        sage: x.str(style="question", error_digits=1)
        '-123456.?2'
        sage: RIF("1.?100").endpoints()  # rel tol 1e-12
        (-99.0, 101.0)
        sage: RIF("1.?100").str(style="question", error_digits=3)
        '1.?100'

    Large exponent (ensure precision is not lost)::

        sage: x = RIF("1.123456?2e100000000"); x
        1.12346?e100000000
        sage: x.str(style="question", error_digits=3)
        '1.12345600?201e100000000'

    Large precision::

        sage: F = RealIntervalField(1000)
        sage: x = F(sqrt(2)); x.endpoints()  # rel tol 1e-290
        (1.41421356237309504880168872420969807856967187537694807317667973799073247846210703885038753432764157273501384623091229702492483605585073721264412149709993583141322266592750559275579995050115278206057147010955997160597027453459686201472851741864088919860955232923048430871432145083976260362799525140798,
         1.41421356237309504880168872420969807856967187537694807317667973799073247846210703885038753432764157273501384623091229702492483605585073721264412149709993583141322266592750559275579995050115278206057147010955997160597027453459686201472851741864088919860955232923048430871432145083976260362799525140799)
        sage: x in F(x.str(style="question", error_digits=3))
        True
        sage: x in F(x.str(style="question", error_digits=0))
        True
        sage: F("1.123456789123456789123456789123456789123456789123456789123456789123456789?987654321987654321987654321e500").endpoints()  # rel tol 1e-290
        (1.123456789123456789123456789123456789123456788135802467135802467135802468e500,
         1.12345678912345678912345678912345678912345679011111111111111111111111111e500)

    Stress test::

        sage: for F in [RealIntervalField(15), RIF, RealIntervalField(100), RealIntervalField(1000)]:
        ....:     for i in range(1000):
        ....:         a, b = randint(-10^9, 10^9), randint(0, 50)
        ....:         c, d = randint(-2^b, 2^b), randint(2, 5)
        ....:         x = a * F(d)^c
        ....:         assert x in F(x.str(style="question", error_digits=3)), (x, a, c, d)
        ....:         assert x in F(x.str(style="question", error_digits=0)), (x, a, c, d)

    Base different from `10` (note that the error and exponent are specified in decimal)::

        sage: RIF("10000.?0", base=2).endpoints()  # rel tol 1e-12
        (16.0, 16.0)
        sage: RIF("10000.?0e10", base=2).endpoints()  # rel tol 1e-12
        (16384.0, 16384.0)
        sage: x = RIF("10000.?10", base=2); x.endpoints()  # rel tol 1e-12
        (6.0, 26.0)
        sage: x.str(base=2, style="question", error_digits=2)
        '10000.000?80'
        sage: x = RIF("10000.000?80", base=2); x.endpoints()  # rel tol 1e-12
        (6.0, 26.0)
        sage: x = RIF("12a.?", base=16); x.endpoints()  # rel tol 1e-12
        (297.0, 299.0)
        sage: x = RIF("12a.BcDeF?", base=16); x.endpoints()  # rel tol 1e-12
        (298.737775802611, 298.737777709962)
        sage: x = RIF("12a.BcDeF?@10", base=16); x.endpoints()  # rel tol 1e-12
        (3.28465658150911e14, 3.28465660248065e14)
        sage: x = RIF("12a.BcDeF?p10", base=16); x.endpoints()  # rel tol 1e-12
        (305907.482421875, 305907.484375000)
        sage: x = RIF("0x12a.BcDeF?p10", base=0); x.endpoints()  # rel tol 1e-12
        (305907.482421875, 305907.484375000)

    Space is allowed::

        sage: RIF("-1234.56?2").endpoints()  # rel tol 1e-12
        (-1234.58, -1234.54)
        sage: RIF("- 1234.56 ?2").endpoints()  # rel tol 1e-12
        (-1234.58, -1234.54)

    Erroneous input::

        sage: RIF("1234.56?2e2.3")
        Traceback (most recent call last):
        ...
        TypeError: unable to convert '1234.56?2e2.3' to real interval
        sage: RIF("1234?2")  # decimal point required
        Traceback (most recent call last):
        ...
        TypeError: unable to convert '1234?2' to real interval
        sage: RIF("1234.?2e")
        Traceback (most recent call last):
        ...
        TypeError: unable to convert '1234.?2e' to real interval
        sage: RIF("1.?e999999999999999999999999")
        [-infinity .. +infinity]
        sage: RIF("0X1.?", base=33)  # X is not valid digit in base 33
        Traceback (most recent call last):
        ...
        TypeError: unable to convert '0X1.?' to real interval
        sage: RIF("1.a?1e10", base=12)
        Traceback (most recent call last):
        ...
        TypeError: unable to convert '1.a?1e10' to real interval
        sage: RIF("1.1?a@10", base=12)
        Traceback (most recent call last):
        ...
        TypeError: unable to convert '1.1?a@10' to real interval
        sage: RIF("0x1?2e1", base=0)  # e is not allowed in base > 10, use @ instead
        Traceback (most recent call last):
        ...
        TypeError: unable to convert '0x1?2e1' to real interval
        sage: RIF("0x1?2p1", base=36)
        Traceback (most recent call last):
        ...
        TypeError: unable to convert '0x1?2p1' to real interval
    """
    cdef mpz_t error_part
    cdef mpfi_t error
    cdef mpfr_t radius, neg_radius
    cdef bytes int_part_string, base_prefix, frac_part_string, error_string, e, sci_expo_string, optional_expo, tmp

    match = NUMBER.fullmatch(s)
    if match is None:
        return 1
    int_part_string, base_prefix, frac_part_string, error_string, e, sci_expo_string = match.groups()

    if (base > 10 or (base == 0 and base_prefix in (b'0X', b'0X'))) and e in (b'e', b'E'):
        return 1
    if base > 16 and e in (b'p', b'P'):
        return 1
    if base > 16 or not base_prefix:
        base_prefix = b''

    if error_string:
        if mpz_init_set_str(error_part, error_string, 10):
            mpz_clear(error_part)
            return 1
    else:
        mpz_init_set_ui(error_part, 1)

    optional_expo = e + sci_expo_string if e else b''
    if mpfi_set_str(x, int_part_string + b'.' + frac_part_string + optional_expo, base):
        mpz_clear(error_part)
        return 1

    mpfr_init2(radius, mpfi_get_prec(x))
    tmp = base_prefix + (
            b'0.' + b'0'*(len(frac_part_string)-1) + b'1' + optional_expo
            if frac_part_string else
            b'1.' + optional_expo)
    # if base = 0:
    #     when s = '-0x123.456@7', tmp = '0x0.001@7'
    #     when s = '-0x123.@7', tmp = '0x1.@7'
    # if base = 36:
    #     when s = '-0x123.456@7', tmp = '0.001@7'
    if mpfr_set_str(radius, tmp, base, MPFR_RNDU):
        mpfr_clear(radius)
        mpz_clear(error_part)
        return 1

    mpfr_mul_z(radius, radius, error_part, MPFR_RNDU)
    mpz_clear(error_part)

    mpfr_init2(neg_radius, mpfi_get_prec(x))
    mpfr_neg(neg_radius, radius, MPFR_RNDD)

    mpfi_init2(error, mpfi_get_prec(x))
    mpfi_interv_fr(error, neg_radius, radius)
    mpfr_clear(radius)
    mpfr_clear(neg_radius)

    mpfi_add(x, x, error)
    mpfi_clear(error)

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

    TESTS::

        sage: RIF('0xabc')
        Traceback (most recent call last):
        ...
        TypeError: unable to convert '0xabc' to real interval
        sage: RIF("0x123.e1", base=0)  # rel tol 1e-12
        291.87890625000000?
        sage: RIF("0x123.@1", base=0)  # rel tol 1e-12
        4656
        sage: RIF("1Xx", base=36)  # rel tol 1e-12
        2517
        sage: RIF("-1Xx@-10", base=62)  # rel tol 1e-12
        -7.088054920481391?e-15
        sage: RIF("1", base=1)
        Traceback (most recent call last):
        ...
        ValueError: base (=1) must be 0 or between 2 and 62
        sage: RIF("1", base=-1)
        Traceback (most recent call last):
        ...
        ValueError: base (=-1) must be 0 or between 2 and 62
        sage: RIF("1", base=63)
        Traceback (most recent call last):
        ...
        ValueError: base (=63) must be 0 or between 2 and 62
    """
    cdef RealIntervalFieldElement ri
    cdef ComplexIntervalFieldElement zi
    cdef ComplexNumber zn
    cdef ComplexDoubleElement zd
    cdef bytes s

    if base != 0 and (base < 2 or base > 62):
        raise ValueError(f"base (={base}) must be 0 or between 2 and 62")
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
            if b"?" in x:
                if _from_str_question_style(re, (<bytes>x).replace(b' ', b''), base):
                    x = bytes_to_str(x)
                    raise TypeError(f"unable to convert {x!r} to real interval")
                return return_real(im)
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
