# cython: cdivision = True
"""
Convert PARI objects to/from Python integers

PARI integers are stored as an array of limbs of type ``pari_ulong``
(which are 32-bit or 64-bit integers). Depending on the kernel
(GMP or native), this array is stored little-endian or big-endian.
This is encapsulated in macros like ``int_W()``:
see section 4.5.1 of the
`PARI library manual <http://pari.math.u-bordeaux.fr/pub/pari/manuals/2.7.0/libpari.pdf>`_.

Python integers of type ``int`` are just C longs. Python integers of
type ``long`` are stored as a little-endian array of type ``digit``
with 15 or 30 bits used per digit. The internal format of a ``long`` is
not documented, but there is some information in
`longintrepr.h <https://github.com/python-git/python/blob/master/Include/longintrepr.h>`_.

Because of this difference in bit lengths, converting integers involves
some bit shuffling.
"""

#*****************************************************************************
#       Copyright (C) 2016 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function

include "cysignals/signals.pxi"
from cpython.int cimport PyInt_AS_LONG
from libc.limits cimport LONG_MIN, LONG_MAX
from .paridecl cimport *
from .pari_instance cimport pari_instance as P

cdef extern from "longintrepr.h":
    cdef _PyLong_New(Py_ssize_t s)
    ctypedef unsigned int digit
    ctypedef struct PyLongObject:
        Py_ssize_t ob_size
        digit* ob_digit

    cdef long PyLong_SHIFT
    cdef digit PyLong_MASK


cpdef integer_to_gen(x):
    """
    Convert a Python ``int`` or ``long`` to a PARI ``gen`` of type
    ``t_INT``.

    EXAMPLES::

        sage: from sage.libs.pari.convert import integer_to_gen
        sage: a = integer_to_gen(int(12345)); a; type(a)
        12345
        <type 'sage.libs.pari.gen.gen'>
        sage: a = integer_to_gen(long(12345)); a; type(a)
        12345
        <type 'sage.libs.pari.gen.gen'>
        sage: integer_to_gen(float(12345))
        Traceback (most recent call last):
        ...
        TypeError: integer_to_gen() needs an int or long argument, not float

    TESTS::

        sage: for i in range(10000):
        ....:     x = 3**i
        ....:     if pari(long(x)) != pari(x):
        ....:         print(x)
    """
    if isinstance(x, int):
        sig_on()
        return P.new_gen(stoi(PyInt_AS_LONG(x)))
    elif isinstance(x, long):
        sig_on()
        return P.new_gen(PyLong_AsGEN(x))
    raise TypeError("integer_to_gen() needs an int or long argument, not {}".format(type(x).__name__))


cpdef gen_to_integer(gen x):
    """
    Convert a PARI ``gen`` to a Python ``int`` or ``long``.

    INPUT:

    - ``x`` -- a PARI ``t_INT``, ``t_FRAC``, ``t_REAL``, a purely
      real ``t_COMPLEX``, a ``t_INTMOD`` or ``t_PADIC`` (which are
      lifted).

    EXAMPLES::

        sage: from sage.libs.pari.convert import gen_to_integer
        sage: a = gen_to_integer(pari("12345")); a; type(a)
        12345
        <type 'int'>
        sage: a = gen_to_integer(pari("10^30")); a; type(a)
        1000000000000000000000000000000L
        <type 'long'>
        sage: gen_to_integer(pari("19/5"))
        3
        sage: gen_to_integer(pari("1 + 0.0*I"))
        1
        sage: gen_to_integer(pari("3/2 + 0.0*I"))
        1
        sage: gen_to_integer(pari("Mod(-1, 11)"))
        10
        sage: gen_to_integer(pari("5 + O(5^10)"))
        5
        sage: gen_to_integer(pari("Pol(42)"))
        42
        sage: gen_to_integer(pari("x"))
        Traceback (most recent call last):
        ...
        TypeError: unable to convert PARI object x of type t_POL to an integer
        sage: gen_to_integer(pari("x + O(x^2)"))
        Traceback (most recent call last):
        ...
        TypeError: unable to convert PARI object x + O(x^2) of type t_SER to an integer
        sage: gen_to_integer(pari("1 + I"))
        Traceback (most recent call last):
        ...
        TypeError: unable to convert PARI object 1 + I of type t_COMPLEX to an integer

    TESTS::

        sage: for i in range(10000):
        ....:     x = 3**i
        ....:     if long(pari(x)) != long(x):
        ....:         print(x)
        sage: gen_to_integer(pari("1.0 - 2^64"))
        -18446744073709551615L
        sage: gen_to_integer(pari("1 - 2^64"))
        -18446744073709551615L

    Check some corner cases::

        sage: for s in [1, -1]:
        ....:     for a in [1, 2^31, 2^32, 2^63, 2^64]:
        ....:         for b in [-1, 0, 1]:
        ....:             Nstr = str(s * (a + b))
        ....:             N1 = gen_to_integer(pari(Nstr))  # Convert via PARI
        ....:             N2 = int(Nstr)                   # Convert via Python
        ....:             if N1 != N2:
        ....:                 print(Nstr, N1, N2)
        ....:             if type(N1) is not type(N2):
        ....:                 print(N1, type(N1), N2, type(N2))
    """
    # First convert the input to a t_INT
    cdef GEN g = gtoi(x.g)

    if not signe(g):
        return 0

    # Try converting to a C long first. Note that we cannot use itos()
    # from PARI since that does not deal with LONG_MIN correctly.
    cdef ulong u
    if lgefint(g) == 3:  # abs(x) fits in a ulong
        u = g[2]         # u = abs(x)
        # Check that <long>(u) or <long>(-u) does not overflow
        if signe(g) >= 0:
            if u <= <ulong>LONG_MAX:
                return <long>(u)
        else:
            if u <= -<ulong>LONG_MIN:
                return <long>(-u)

    # Result does not fit in a C long
    return PyLong_FromGEN(g)


cdef GEN gtoi(GEN g0) except NULL:
    """
    Convert a PARI object to a PARI integer.

    This function is shallow and not stack-clean.
    """
    if typ(g0) == t_INT:
        return g0
    cdef GEN g
    try:
        sig_on()
        g = simplify_shallow(g0)
        if typ(g) == t_COMPLEX:
            if gequal0(gel(g,2)):
                g = gel(g,1)
        if typ(g) == t_INTMOD:
            g = gel(g,2)
        g = trunc_safe(g)
        if typ(g) != t_INT:
            sig_error()
        sig_off()
    except RuntimeError:
        raise TypeError(stack_sprintf(
            "unable to convert PARI object %Ps of type %s to an integer",
            g0, type_name(typ(g0))))
    return g


cdef GEN PyLong_AsGEN(x):
    cdef PyLongObject* L = <PyLongObject*>(x)
    cdef const digit* D = L.ob_digit

    # Size of the input
    cdef size_t sizedigits
    cdef long sgn
    if L.ob_size == 0:
        return gen_0
    elif L.ob_size > 0:
        sizedigits = L.ob_size
        sgn = evalsigne(1)
    else:
        sizedigits = -L.ob_size
        sgn = evalsigne(-1)

    # Size of the output, in bits and in words
    cdef size_t sizebits = sizedigits * PyLong_SHIFT
    cdef size_t sizewords = (sizebits + BITS_IN_LONG - 1) // BITS_IN_LONG

    # Compute the most significant word of the output.
    # This is a special case because we need to be careful not to
    # overflow the ob_digit array. We also need to check for zero,
    # in which case we need to decrease sizewords.
    # See the loop below for an explanation of this code.
    cdef size_t bit = (sizewords - 1) * BITS_IN_LONG
    cdef size_t dgt = bit // PyLong_SHIFT
    bit = bit % PyLong_SHIFT

    cdef ulong w = <ulong>(D[dgt]) >> bit
    if 1*PyLong_SHIFT - bit < BITS_IN_LONG and dgt+1 < sizedigits:
        w += <ulong>(D[dgt+1]) << (1*PyLong_SHIFT - bit)
    if 2*PyLong_SHIFT - bit < BITS_IN_LONG and dgt+2 < sizedigits:
        w += <ulong>(D[dgt+2]) << (2*PyLong_SHIFT - bit)
    if 3*PyLong_SHIFT - bit < BITS_IN_LONG and dgt+3 < sizedigits:
        w += <ulong>(D[dgt+3]) << (3*PyLong_SHIFT - bit)
    if 4*PyLong_SHIFT - bit < BITS_IN_LONG and dgt+4 < sizedigits:
        w += <ulong>(D[dgt+4]) << (4*PyLong_SHIFT - bit)
    if 5*PyLong_SHIFT - bit < BITS_IN_LONG and dgt+5 < sizedigits:
        w += <ulong>(D[dgt+5]) << (5*PyLong_SHIFT - bit)

    # Effective size in words plus 2 special codewords
    cdef pariwords = sizewords+2 if w else sizewords+1
    cdef GEN g = cgeti(pariwords)
    g[1] = sgn + evallgefint(pariwords)

    if w:
        int_MSW(g)[0] = w

    # Fill all words
    cdef GEN ptr = int_LSW(g)
    cdef size_t i
    for i in range(sizewords - 1):
        # The least significant bit of word number i of the output
        # integer is bit number "bit" of Python digit "dgt".
        bit = i * BITS_IN_LONG
        dgt = bit // PyLong_SHIFT
        bit = bit % PyLong_SHIFT

        # Now construct the output word from the Python digits:
        # 6 digits are enough assuming that PyLong_SHIFT >= 15 and
        # BITS_IN_LONG <= 76.  The compiler should optimize away all
        # but one of the "if" statements below.
        w = <ulong>(D[dgt]) >> bit
        if 1*PyLong_SHIFT - bit < BITS_IN_LONG:
            w += <ulong>(D[dgt+1]) << (1*PyLong_SHIFT - bit)
        if 2*PyLong_SHIFT - bit < BITS_IN_LONG:
            w += <ulong>(D[dgt+2]) << (2*PyLong_SHIFT - bit)
        if 3*PyLong_SHIFT - bit < BITS_IN_LONG:
            w += <ulong>(D[dgt+3]) << (3*PyLong_SHIFT - bit)
        if 4*PyLong_SHIFT - bit < BITS_IN_LONG:
            w += <ulong>(D[dgt+4]) << (4*PyLong_SHIFT - bit)
        if 5*PyLong_SHIFT - bit < BITS_IN_LONG:
            w += <ulong>(D[dgt+5]) << (5*PyLong_SHIFT - bit)

        ptr[0] = w
        ptr = int_nextW(ptr)

    return g


cdef PyLong_FromGEN(GEN g):
    # Size of input in words, bits and Python digits. The size in
    # digits might be a small over-estimation, but that is not a
    # problem.
    cdef size_t sizewords = (lgefint(g) - 2)
    cdef size_t sizebits = sizewords * BITS_IN_LONG
    cdef size_t sizedigits = (sizebits + PyLong_SHIFT - 1) // PyLong_SHIFT

    # Actual correct computed size
    cdef Py_ssize_t sizedigits_final = 0

    x = _PyLong_New(sizedigits)
    cdef PyLongObject* L = <PyLongObject*>(x)
    cdef digit* D = L.ob_digit

    cdef digit d
    cdef ulong w
    cdef size_t i, j, bit
    for i in range(sizedigits):
        # The least significant bit of digit number i of the output
        # integer is bit number "bit" of word "j".
        bit = i * PyLong_SHIFT
        j = bit // BITS_IN_LONG
        bit = bit % BITS_IN_LONG

        w = int_W(g, j)[0]
        d = w >> bit

        # Do we need bits from the next word too?
        if BITS_IN_LONG - bit < PyLong_SHIFT and j+1 < sizewords:
            w = int_W(g, j+1)[0]
            d += w << (BITS_IN_LONG - bit)

        d = d & PyLong_MASK
        D[i] = d

        # Keep track of last non-zero digit
        if d:
            sizedigits_final = i+1

    # Set correct size
    if signe(g) > 0:
        L.ob_size = sizedigits_final
    else:
        L.ob_size = -sizedigits_final

    return x
