"""
Generic implementation of powering

This implements powering of arbitrary objects using a
square-and-multiply algorithm.
"""

# ****************************************************************************
#       Copyright (C) 2017 Jeroen Demeyer <J.Demeyer@UGent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.signals cimport sig_check

from .long cimport integer_check_long


cpdef generic_power(a, n):
    """
    Return `a^n`.

    If `n` is negative, return `(1/a)^{-n}`.

    INPUT:

    - ``a`` -- any object supporting multiplication
      (and division if n < 0)

    - ``n`` -- any integer (in the duck typing sense)

    EXAMPLES::

        sage: from sage.arith.power import generic_power
        sage: generic_power(int(12), int(0))
        1
        sage: generic_power(int(0), int(100))
        0
        sage: generic_power(Integer(10), Integer(0))
        1
        sage: generic_power(Integer(0), Integer(23))
        0
        sage: sum([generic_power(2,i) for i in range(17)]) #test all 4-bit combinations
        131071
        sage: F = Zmod(5)
        sage: a = generic_power(F(2), 5); a
        2
        sage: a.parent() is F
        True
        sage: a = generic_power(F(1), 2)
        sage: a.parent() is F
        True

        sage: generic_power(int(5), 0)
        1
        sage: generic_power(2, 5/4)
        Traceback (most recent call last):
        ...
        NotImplementedError: non-integral exponents not supported

    ::

        sage: class SymbolicMul(str):
        ....:     def __mul__(self, other):
        ....:         s = "({}*{})".format(self, other)
        ....:         return type(self)(s)
        sage: x = SymbolicMul("x")
        sage: print(generic_power(x, 7))
        (((x*x)*(x*x))*((x*x)*x))
    """
    if not n:
        return one(a)

    cdef long value = 0
    cdef int err
    if not integer_check_long(n, &value, &err):
        raise NotImplementedError("non-integral exponents not supported")
    if not err:
        return generic_power_long(a, value)

    if n < 0:
        n = -n
        a = invert(a)
    return generic_power_pos(a, n)


cdef generic_power_long(a, long n):
    """
    As ``generic_power`` but where ``n`` is a C long.
    """
    if not n:
        return one(a)

    cdef unsigned long u = <unsigned long>n
    if n < 0:
        u = -u
        a = invert(a)
    return generic_power_pos(a, u)


cdef generic_power_pos(a, ulong_or_object n):
    """
    Return `a^n` where `n > 0`.
    """
    # Find least significant set bit as starting point
    apow = a
    while not (n & 1):
        sig_check()
        apow *= apow
        n >>= 1

    # Now multiply together the correct factors a^(2^i)
    res = apow
    n >>= 1
    while n:
        sig_check()
        apow *= apow
        if n & 1:
            res = apow * res
        n >>= 1

    return res
