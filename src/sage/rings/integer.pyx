r"""
Elements of the ring `\ZZ` of integers

AUTHORS:

- William Stein (2005): initial version

- Gonzalo Tornaria (2006-03-02): vastly improved python/GMP
  conversion; hashing

- Didier Deshommes (2006-03-06): numerous examples
  and docstrings

- William Stein (2006-03-31): changes to reflect GMP bug fixes

- William Stein (2006-04-14): added GMP factorial method (since it's
  now very fast).

- David Harvey (2006-09-15): added nth_root, exact_log

- David Harvey (2006-09-16): attempt to optimise Integer constructor

- Rishikesh (2007-02-25): changed quo_rem so that the rem is positive

- David Harvey, Martin Albrecht, Robert Bradshaw (2007-03-01):
  optimized Integer constructor and pool

- Pablo De Napoli (2007-04-01): multiplicative_order should return
  +infinity for non zero numbers

- Robert Bradshaw (2007-04-12): is_perfect_power, Jacobi symbol (with
  Kronecker extension).  Convert some methods to use GMP directly
  rather than PARI, Integer(), PY_NEW(Integer)

- David Roe (2007-03-21): sped up valuation and is_square, added
  val_unit, is_power, is_power_of and divide_knowing_divisible_by

- Robert Bradshaw (2008-03-26): gamma function, multifactorials

- Robert Bradshaw (2008-10-02): bounded squarefree part

- David Loeffler (2011-01-15): fixed bug #10625 (inverse_mod should accept an ideal as argument)

- Vincent Delecroix (2010-12-28): added unicode in Integer.__init__

- David Roe (2012-03): deprecate :meth:`~sage.rings.integer.Integer.is_power`
  in favour of :meth:`~sage.rings.integer.Integer.is_perfect_power` (see
  :trac:`12116`)

EXAMPLES:

Add 2 integers::

    sage: a = Integer(3) ; b = Integer(4)
    sage: a + b == 7
    True

Add an integer and a real number::

    sage: a + 4.0
    7.00000000000000

Add an integer and a rational number::

    sage: a + Rational(2)/5
    17/5

Add an integer and a complex number::

    sage: b = ComplexField().0 + 1.5
    sage: loads((a+b).dumps()) == a+b
    True

    sage: z = 32
    sage: -z
    -32
    sage: z = 0; -z
    0
    sage: z = -0; -z
    0
    sage: z = -1; -z
    1

Multiplication::

    sage: a = Integer(3) ; b = Integer(4)
    sage: a * b == 12
    True
    sage: loads((a * 4.0).dumps()) == a*b
    True
    sage: a * Rational(2)/5
    6/5

::

    sage: list([2,3]) * 4
    [2, 3, 2, 3, 2, 3, 2, 3]

::

    sage: 'sage'*Integer(3)
    'sagesagesage'

COERCIONS:

Returns version of this integer in the multi-precision floating
real field R::

    sage: n = 9390823
    sage: RR = RealField(200)
    sage: RR(n)
    9.3908230000000000000000000000000000000000000000000000000000e6

"""
#*****************************************************************************
#       Copyright (C) 2004,2006 William Stein <wstein@gmail.com>
#       Copyright (C) 2006 Gonzalo Tornaria <tornaria@math.utexas.edu>
#       Copyright (C) 2006 Didier Deshommes <dfdeshom@gmail.com>
#       Copyright (C) 2007 David Harvey <dmharvey@math.harvard.edu>
#       Copyright (C) 2007 Martin Albrecht <malb@informatik.uni-bremen.de>
#       Copyright (C) 2007,2008 Robert Bradshaw <robertwb@math.washington.edu>
#       Copyright (C) 2007 David Roe <roed314@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

# Do not create any Integer, especially non cdef'ed ones, before the hooked
# creation and deletion are setup by the call to hook_fast_tp_functions

doc="""
Integers
"""

import cython

import operator

import sys

include "sage/ext/gmp.pxi"
include "sage/ext/interrupt.pxi"  # ctrl-c interrupt block support
include "sage/ext/stdsage.pxi"
from cpython.list cimport *
from cpython.number cimport *
from cpython.int cimport *
include "sage/ext/python_debug.pxi"
include "../structure/coerce.pxi"   # for parent_c
include "sage/libs/pari/decl.pxi"

cdef extern from "limits.h":
    long LONG_MAX

cdef extern from "math.h":
    double sqrt_double "sqrt"(double)
    cdef double log_c "log" (double)
    cdef double ceil_c "ceil" (double)
    int isnan(double)

cdef extern from "mpz_pylong.h":
    cdef mpz_get_pylong(mpz_t src)
    cdef mpz_get_pyintlong(mpz_t src)
    cdef int mpz_set_pylong(mpz_t dst, src) except -1
    cdef long mpz_pythonhash(mpz_t src)

cdef extern from "mpz_longlong.h":
    cdef void mpz_set_longlong(mpz_t dst, long long src)
    cdef void mpz_set_ulonglong(mpz_t dst, unsigned long long src)

cdef extern from "convert.h":
    cdef void t_INT_to_ZZ( mpz_t value, long *g )

from sage.libs.pari.gen cimport gen as pari_gen
from sage.libs.pari.pari_instance cimport PariInstance
from sage.libs.flint.ulong_extras cimport *

import sage.rings.infinity

import sage.libs.pari.pari_instance
cdef PariInstance pari = sage.libs.pari.pari_instance.pari

from sage.structure.element import canonical_coercion
from sage.misc.superseded import deprecated_function_alias

cdef object numpy_long_interface = {'typestr': '=i4' if sizeof(long) == 4 else '=i8' }
cdef object numpy_int64_interface = {'typestr': '=i8'}
cdef object numpy_object_interface = {'typestr': '|O'}

cdef mpz_t mpz_tmp
mpz_init(mpz_tmp)

def init_mpz_globals():
    init_mpz_globals_c()

def clear_mpz_globals():
    clear_mpz_globals_c()

cdef int set_mpz(Integer self, mpz_t value):
    mpz_set(self.value, value)

cdef set_from_Integer(Integer self, Integer other):
    mpz_set(self.value, other.value)

cdef set_from_int(Integer self, int other):
    mpz_set_si(self.value, other)

cdef mpz_t* get_value(Integer self):
    return &self.value


def _test_mpz_set_longlong(long long v):
    """
    TESTS::

        sage: from sage.rings.integer import _test_mpz_set_longlong
        sage: _test_mpz_set_longlong(1)
        1
        sage: _test_mpz_set_longlong(-1)
        -1
        sage: _test_mpz_set_longlong(100)
        100
        sage: _test_mpz_set_longlong(100000000000)
        100000000000
        sage: _test_mpz_set_longlong(2^32) == 2^32
        True
        sage: _test_mpz_set_longlong(-2^32) == -2^32
        True
        sage: all([_test_mpz_set_longlong(2^n) == 2^n for n in range(63)])
        True
    """
    cdef Integer z = PY_NEW(Integer)
    mpz_set_longlong(z.value, v)
    return z

cdef _digits_naive(mpz_t v,l,int offset,Integer base,digits):
    """
    This method fills in digit entries in the list, l, using the most
    basic digit algorithm -- repeat division by base.

    INPUT:

    - ``v`` - the value whose digits we want to put into the list

    - ``l`` - the list to file

    - ``offset`` - offset from the beginning of the list that we want
      to fill at

    - ``base`` -- the base to which we finding digits

    - ``digits`` - a python sequence type with objects to use for digits
                 note that python negative index semantics are relied upon

    AUTHORS:

    - Joel B. Mohler (2009-01-16)
    """
    cdef mpz_t mpz_value
    cdef mpz_t mpz_res  # used on one side of the 'if'
    cdef Integer z      # used on the other side of the 'if'

    mpz_init(mpz_value)
    mpz_set(mpz_value, v)

    # we aim to avoid sage Integer creation if possible
    if digits is None:
        while mpz_cmp_si(mpz_value,0):
            z = PY_NEW(Integer)
            mpz_tdiv_qr(mpz_value, z.value, mpz_value, base.value)
            l[offset] = z
            offset += 1
    else:
        mpz_init(mpz_res)
        while mpz_cmp_si(mpz_value,0):
            mpz_tdiv_qr(mpz_value, mpz_res, mpz_value, base.value)
            l[offset] = digits[mpz_get_si(mpz_res)]
            offset += 1
        mpz_clear(mpz_res)

    mpz_clear(mpz_value)

cdef _digits_internal(mpz_t v,l,int offset,int power_index,power_list,digits):
    """
    INPUT:

    - ``v`` - the value whose digits we want to put into the list

    - ``l`` - the list to file

    - ``offset`` - offset from the beginning of the list that we want
      to fill at

    - ``power_index`` - a measure of size to fill and index to
      power_list we're filling 1 << (power_index+1) digits

    - ``power_list`` - a list of powers of the base, precomputed in
      method digits digits - a python sequence type with objects to
      use for digits note that python negative index semantics are
      relied upon

    AUTHORS:

    - Joel B. Mohler (2008-03-13)
    """
    cdef mpz_t mpz_res
    cdef mpz_t mpz_quot
    cdef Integer temp
    cdef int v_int
    if power_index < 5:
        # It turns out that simple repeated division is very fast for
        # relatively few digits.  I don't think this is a real algorithmic
        # statement, it's an annoyance introduced by memory allocation.
        # I think that manual memory management with mpn_* would make the
        # divide & conquer approach even faster, but the code would be much
        # more complicated.
        _digits_naive(v,l,offset,power_list[0],digits)
    else:
        mpz_init(mpz_quot)
        mpz_init(mpz_res)
        temp = power_list[power_index]
        mpz_tdiv_qr(mpz_quot, mpz_res, v, temp.value)
        if mpz_sgn(mpz_res) != 0:
            _digits_internal(mpz_res,l,offset,power_index-1,power_list,digits)
        if mpz_sgn(mpz_quot) != 0:
            _digits_internal(mpz_quot,l,offset+(1<<power_index),power_index-1,power_list,digits)
        mpz_clear(mpz_quot)
        mpz_clear(mpz_res)

arith = None
cdef void late_import():
    global arith
    if arith is None:
        import sage.rings.arith
        arith = sage.rings.arith

MAX_UNSIGNED_LONG = 2 * sys.maxint

from sage.structure.sage_object cimport SageObject
from sage.structure.element cimport EuclideanDomainElement, ModuleElement, Element
from sage.structure.element import  bin_op
from sage.structure.coerce_exceptions import CoercionException

import integer_ring
the_integer_ring = integer_ring.ZZ

# The documentation for the ispseudoprime() function in the PARI
# manual states that its result is always prime up to this 10^13.
cdef mpz_t PARI_PSEUDOPRIME_LIMIT
mpz_init(PARI_PSEUDOPRIME_LIMIT)
mpz_set_str(PARI_PSEUDOPRIME_LIMIT, "10000000000000", 10)

def is_Integer(x):
    """
    Return true if x is of the Sage integer type.

    EXAMPLES::

        sage: from sage.rings.integer import is_Integer
        sage: is_Integer(2)
        True
        sage: is_Integer(2/1)
        False
        sage: is_Integer(int(2))
        False
        sage: is_Integer(long(2))
        False
        sage: is_Integer('5')
        False
    """
    return PY_TYPE_CHECK(x, Integer)

cdef inline Integer as_Integer(x):
    if PY_TYPE_CHECK(x, Integer):
        return <Integer>x
    else:
        return Integer(x)

cdef class IntegerWrapper(Integer):
    r"""
    Rationale for the ``IntegerWrapper`` class:

    With ``Integers``, the allocation/deallocation function slots are
    hijacked with custom functions that stick already allocated
    ``Integers`` (with initialized ``parent`` and ``mpz_t`` fields)
    into a pool on "deallocation" and then pull them out whenever a
    new one is needed. Because ``Integers`` are so common, this is
    actually a significant savings. However , this does cause issues
    with subclassing a Python class directly from ``Integer`` (but
    that's ok for a Cython class).

    As a workaround, one can instead derive a class from the
    intermediate class ``IntegerWrapper``, which sets statically its
    alloc/dealloc methods to the *original* ``Integer`` alloc/dealloc
    methods, before they are swapped manually for the custom ones.

    The constructor of ``IntegerWrapper`` further allows for
    specifying an alternative parent to ``IntegerRing()``.
    """

    def __init__(self, parent=None, x=None, unsigned int base=0):
        """
        We illustrate how to create integers with parents different
        from ``IntegerRing()``::

            sage: from sage.rings.integer import IntegerWrapper

            sage: n = IntegerWrapper(Primes(), 3) # indirect doctest
            sage: n
            3
            sage: n.parent()
            Set of all prime numbers: 2, 3, 5, 7, ...

        Pickling seems to work now (as of #10314)

            sage: nn = loads(dumps(n))
            sage: nn
            3
            sage: nn.parent()
            Integer Ring

            sage: TestSuite(n).run()
        """
        if parent is not None:
            Element.__init__(self, parent=parent)
        Integer.__init__(self, x, base=base)

cdef class Integer(sage.structure.element.EuclideanDomainElement):
    r"""
    The ``Integer`` class represents arbitrary precision
    integers. It derives from the ``Element`` class, so
    integers can be used as ring elements anywhere in Sage.

    Integer() interprets numbers and strings that begin with 0 as octal
    numbers, and numbers and strings that begin with 0x as hexadecimal
    numbers.

        The class ``Integer`` is implemented in Cython, as a
        wrapper of the GMP ``mpz_t`` integer type.

    EXAMPLES::

        sage: Integer(010)
        8
        sage: Integer(0x10)
        16
        sage: Integer(10)
        10
        sage: Integer('0x12')
        18
        sage: Integer('012')
        10

    Conversion from PARI::

        sage: Integer(pari('-10380104371593008048799446356441519384'))
        -10380104371593008048799446356441519384
        sage: Integer(pari('Pol([-3])'))
        -3
    """

    def __cinit__(self):
        mpz_init(self.value)
        self._parent = <SageObject>the_integer_ring

    def __init__(self, x=None, base=0):
        """
        EXAMPLES::

            sage: a = long(-901824309821093821093812093810928309183091832091)
            sage: b = ZZ(a); b
            -901824309821093821093812093810928309183091832091
            sage: ZZ(b)
            -901824309821093821093812093810928309183091832091
            sage: ZZ('-901824309821093821093812093810928309183091832091')
            -901824309821093821093812093810928309183091832091
            sage: ZZ(int(-93820984323))
            -93820984323
            sage: ZZ(ZZ(-901824309821093821093812093810928309183091832091))
            -901824309821093821093812093810928309183091832091
            sage: ZZ(QQ(-901824309821093821093812093810928309183091832091))
            -901824309821093821093812093810928309183091832091
            sage: ZZ(RR(2.0)^80)
            1208925819614629174706176
            sage: ZZ(QQbar(sqrt(28-10*sqrt(3)) + sqrt(3)))
            5
            sage: ZZ(AA(32).nth_root(5))
            2
            sage: ZZ(pari('Mod(-3,7)'))
            4
            sage: ZZ('sage')
            Traceback (most recent call last):
            ...
            TypeError: unable to convert x (=sage) to an integer
            sage: Integer('zz',36).str(36)
            'zz'
            sage: ZZ('0x3b').str(16)
            '3b'
            sage: ZZ( ZZ(5).digits(3) , 3)
            5
            sage: import numpy
            sage: ZZ(numpy.int64(7^7))
            823543
            sage: ZZ(numpy.ubyte(-7))
            249
            sage: ZZ(True)
            1
            sage: ZZ(False)
            0
            sage: ZZ(1==0)
            0
            sage: ZZ('+10')
            10

        ::

            sage: k = GF(2)
            sage: ZZ( (k(0),k(1)), 2)
            2

        ::

            sage: ZZ(float(2.0))
            2
            sage: ZZ(float(1.0/0.0))
            Traceback (most recent call last):
            ...
            OverflowError: cannot convert float infinity to integer
            sage: ZZ(float(0.0/0.0))
            Traceback (most recent call last):
            ...
            ValueError: cannot convert float NaN to integer

        ::

            sage: class MyInt(int):
            ...       pass
            sage: class MyLong(long):
            ...       pass
            sage: class MyFloat(float):
            ...       pass
            sage: ZZ(MyInt(3))
            3
            sage: ZZ(MyLong(4))
            4
            sage: ZZ(MyFloat(5))
            5

        ::

            sage: Integer(u'0')
            0
            sage: Integer(u'0X2AEEF')
            175855

        Test conversion from PARI (#11685)::

            sage: ZZ(pari(-3))
            -3
            sage: ZZ(pari("-3.0"))
            -3
            sage: ZZ(pari("-3.5"))
            Traceback (most recent call last):
            ...
            TypeError: Attempt to coerce non-integral real number to an Integer
            sage: ZZ(pari("1e100"))
            Traceback (most recent call last):
            ...
            PariError: precision too low in truncr (precision loss in truncation)
            sage: ZZ(pari("10^50"))
            100000000000000000000000000000000000000000000000000
            sage: ZZ(pari("Pol(3)"))
            3
            sage: ZZ(GF(3^20,'t')(1))
            1
            sage: ZZ(pari(GF(3^20,'t')(1)))
            1
            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^2+3)
            sage: ZZ(a^2)
            -3
            sage: ZZ(pari(a)^2)
            -3
            sage: ZZ(pari("Mod(x, x^3+x+1)"))   # Note error message refers to lifted element
            Traceback (most recent call last):
            ...
            TypeError: Unable to coerce PARI x to an Integer

        Test coercion of p-adic with negative valuation::

            sage: ZZ(pari(Qp(11)(11^-7)))
            Traceback (most recent call last):
            ...
            TypeError: Cannot convert p-adic with negative valuation to an integer

        Test converting a list with a very large base::

            sage: a=ZZ(randint(0,2^128-1))
            sage: L = a.digits(2^64)
            sage: a == sum([x * 2^(64*i) for i,x in enumerate(L)])
            True
            sage: a == ZZ(L,base=2^64)
            True
        """
        # TODO: All the code below should somehow be in an external
        # cdef'd function.  Then e.g., if a matrix or vector or
        # polynomial is getting filled by mpz_t's, it can use the
        # rules below to do the fill construction of mpz_t's, but
        # without the overhead of creating any Python objects at all.
        # The cdef's function should be of the form
        #     mpz_init_set_sage(mpz_t y, object x)
        # Then this function becomes the one liner:
        #     mpz_init_set_sage(self.value, x)

        cdef Integer tmp
        cdef char* xs
        cdef int paritype
        cdef Py_ssize_t j
        cdef object otmp
        cdef unsigned int ibase

        cdef Element lift

        if x is None:
            if mpz_sgn(self.value) != 0:
                mpz_set_si(self.value, 0)

        else:
            # First do all the type-check versions; these are fast.

            if PY_TYPE_CHECK(x, Integer):
                set_from_Integer(self, <Integer>x)

            elif PY_TYPE_CHECK(x, bool):
                mpz_set_si(self.value, PyInt_AS_LONG(x))

            elif PyInt_Check(x):
                mpz_set_si(self.value, PyInt_AS_LONG(x))

            elif PyLong_Check(x):
                mpz_set_pylong(self.value, x)

            elif PyFloat_Check(x):
                n = long(x)
                if n == x:
                    mpz_set_pylong(self.value, n)
                else:
                    raise TypeError, "Cannot convert non-integral float to integer"

            elif PY_TYPE_CHECK(x, pari_gen):
                # Simplify and lift until we get an integer
                while typ((<pari_gen>x).g) != t_INT:
                    x = x.simplify()
                    paritype = typ((<pari_gen>x).g)
                    if paritype == t_INT:
                        break
                    elif paritype == t_REAL:
                        # Check that the fractional part is zero
                        if not x.frac().gequal0():
                            raise TypeError, "Attempt to coerce non-integral real number to an Integer"
                        # floor yields an integer
                        x = x.floor()
                        break
                    elif paritype == t_PADIC:
                        if x._valp() < 0:
                            raise TypeError("Cannot convert p-adic with negative valuation to an integer")
                        # Lifting a PADIC yields an integer
                        x = x.lift()
                        break
                    elif paritype == t_INTMOD:
                        # Lifting an INTMOD yields an integer
                        x = x.lift()
                        break
                    elif paritype == t_POLMOD:
                        x = x.lift()
                    elif paritype == t_FFELT:
                        # x = (f modulo defining polynomial of finite field);
                        # we extract f.
                        sig_on()
                        x = pari.new_gen(FF_to_FpXQ_i((<pari_gen>x).g))
                    else:
                        raise TypeError, "Unable to coerce PARI %s to an Integer"%x

                # Now we have a true PARI integer, convert it to Sage
                t_INT_to_ZZ(self.value, (<pari_gen>x).g)

            elif PyString_Check(x) or PY_TYPE_CHECK(x,unicode):
                if base < 0 or base > 36:
                    raise ValueError, "`base` (=%s) must be between 2 and 36."%base
                ibase = base
                xs = x
                if xs[0] == c'+':
                    xs += 1
                if mpz_set_str(self.value, xs, ibase) != 0:
                    raise TypeError, "unable to convert x (=%s) to an integer"%x

            elif PyObject_HasAttrString(x, "_integer_"):
                # TODO: Note that PyObject_GetAttrString returns NULL if
                # the attribute was not found. If we could test for this,
                # we could skip the double lookup. Unfortunately Cython doesn't
                # seem to let us do this; it flags an error if the function
                # returns NULL, because it can't construct an "object" object
                # out of the NULL pointer. This really sucks. Perhaps we could
                # make the function prototype have return type void*, but
                # then how do we make Cython handle the reference counting?
                set_from_Integer(self, (<object> PyObject_GetAttrString(x, "_integer_"))(the_integer_ring))

            elif (PY_TYPE_CHECK(x, list) or PY_TYPE_CHECK(x, tuple)) and base > 1:
                b = the_integer_ring(base)
                if b == 2: # we use a faster method
                    for j from 0 <= j < len(x):
                        otmp = x[j]
                        if not PY_TYPE_CHECK(otmp, Integer):
                            # should probably also have fast code for Python ints...
                            otmp = Integer(otmp)
                        if mpz_cmp_si((<Integer>otmp).value, 1) == 0:
                            mpz_setbit(self.value, j)
                        elif mpz_sgn((<Integer>otmp).value) != 0:
                            # one of the entries was something other than 0 or 1.
                            break
                    else:
                        return
                tmp = the_integer_ring(0)
                for i in range(len(x)):
                    tmp += the_integer_ring(x[i])*b**i
                mpz_set(self.value, tmp.value)

            else:
                import numpy
                if isinstance(x, numpy.integer):
                    mpz_set_pylong(self.value, x.__long__())
                    return

                elif PY_TYPE_CHECK(x, Element):
                    try:
                        lift = x.lift()
                        if lift._parent != (<Element>x)._parent:
                            tmp = the_integer_ring(lift)
                            mpz_swap(tmp.value, self.value)
                            return
                    except AttributeError:
                        pass

                raise TypeError, "unable to coerce %s to an integer" % type(x)

    def __reduce__(self):
        """
        This is used when pickling integers.

        EXAMPLES::

            sage: n = 5
            sage: t = n.__reduce__(); t
            (<built-in function make_integer>, ('5',))
            sage: t[0](*t[1])
            5
            sage: loads(dumps(n)) == n
            True
        """
        # This single line below took me HOURS to figure out.
        # It is the *trick* needed to pickle Cython extension types.
        # The trick is that you must put a pure Python function
        # as the first argument, and that function must return
        # the result of unpickling with the argument in the second
        # tuple as input. All kinds of problems happen
        # if we don't do this.
        return sage.rings.integer.make_integer, (self.str(32),)

    cdef _reduce_set(self, s):
        """
        Set this integer from a string in base 32.

        .. note::

           Integers are supposed to be immutable, so you should not
           use this function.
        """
        mpz_set_str(self.value, s, 32)

    def __index__(self):
        """
        Needed so integers can be used as list indices.

        EXAMPLES::

            sage: v = [1,2,3,4,5]
            sage: v[Integer(3)]
            4
            sage: v[Integer(2):Integer(4)]
            [3, 4]
        """
        return mpz_get_pyintlong(self.value)

    def _im_gens_(self, codomain, im_gens):
        """
        Return the image of self under the map that sends the generators of
        the parent to im_gens. Since ZZ maps canonically in the category
        of rings, this is just the natural coercion.

        EXAMPLES::

            sage: n = -10
            sage: R = GF(17)
            sage: n._im_gens_(R, [R(1)])
            7
        """
        return codomain._coerce_(self)

    cdef _xor(Integer self, Integer other):
        cdef Integer x
        x = PY_NEW(Integer)
        mpz_xor(x.value, self.value, other.value)
        return x

    def __xor__(x, y):
        """
        Compute the exclusive or of x and y.

        EXAMPLES::

            sage: n = ZZ(2); m = ZZ(3)
            sage: n.__xor__(m)
            1
        """
        if PY_TYPE_CHECK(x, Integer) and PY_TYPE_CHECK(y, Integer):
            return (<Integer>x)._xor(y)
        return bin_op(x, y, operator.xor)

    def __richcmp__(left, right, int op):
        """
        cmp for integers

        EXAMPLES::

            sage: 2 < 3
            True
            sage: 2 > 3
            False
            sage: 2 == 3
            False
            sage: 3 > 2
            True
            sage: 3 < 2
            False

            sage: 1000000000000000000000000000000000000000000000000000.0r==1000000000000000000000000000000000000000000000000000
            False
            sage: 1000000000000000000000000000000000000000000000000000.1r==1000000000000000000000000000000000000000000000000000
            False

        Canonical coercions are used but non-canonical ones are not.

        ::

            sage: 4 == 4/1
            True
            sage: 4 == '4'
            False

        TESTS::

            sage: 3 < 4r
            True
            sage: 3r < 4
            True
            sage: 3 >= 4r
            False
            sage: 4r <= 3
            False
            sage: 12345678901234567890123456789r == 12345678901234567890123456789
            True
            sage: 12345678901234567890123456788 < 12345678901234567890123456789r
            True
            sage: 2 < 2.7r
            True
            sage: 4 < 3.1r
            False
            sage: -1r < 1
            True
            sage: -1.5r < 3
            True
            sage: Ilist = [-2,-1,0,1,2,12345678901234567890123456788]
            sage: ilist = [-4r,-1r,0r,1r,2r,5r]
            sage: llist = [-12345678901234567890123456788r, 12345678901234567890123456788r, 12345678901234567890123456900r]
            sage: flist = [-21.8r, -1.2r, -.000005r, 0.0r, .999999r, 1000000000000.0r]
            sage: all([(a < b) == (RR(a) < RR(b)) for (a, b) in zip(Ilist, ilist)])
            True
            sage: all([(a > b) == (RR(a) > RR(b)) for (a, b) in zip(Ilist, ilist)])
            True
            sage: all([(a == b) == (RR(a) == RR(b)) for (a, b) in zip(Ilist, ilist)])
            True
            sage: all([(a <= b) == (RR(a) <= RR(b)) for (a, b) in zip(Ilist, ilist)])
            True
            sage: all([(a >= b) == (RR(a) >= RR(b)) for (a, b) in zip(Ilist, ilist)])
            True
            sage: all([(a != b) == (RR(a) != RR(b)) for (a, b) in zip(Ilist, ilist)])
            True
            sage: all([(a < b) == (RR(a) < RR(b)) for (a, b) in zip(Ilist, llist)])
            True
            sage: all([(a > b) == (RR(a) > RR(b)) for (a, b) in zip(Ilist, llist)])
            True
            sage: all([(a < b) == (RR(a) < RR(b)) for (a, b) in zip(Ilist, flist)])
            True
            sage: all([(a > b) == (RR(a) > RR(b)) for (a, b) in zip(Ilist, flist)])
            True

        Verify that trac 12149 was fixed (and the fix is consistent
        with Python ints)::

            sage: a = int(1); b = 1; n = float('nan')
            sage: a == n
            False
            sage: a == n, b == n
            (False, False)
            sage: a != n, b != n, n != b
            (True, True, True)
            sage: a < n, b < n, n > b
            (False, False, False)
            sage: a > n, b > n, n < b
            (False, False, False)
            sage: a <= n, b <= n, n >= b
            (False, False, False)
            sage: a >= n, b >= n, n <= b
            (False, False, False)
        """
        cdef int c
        cdef double d
        if PY_TYPE_CHECK(left, Integer):
            if PY_TYPE_CHECK(right, Integer):
                c = mpz_cmp((<Integer>left).value, (<Integer>right).value)
            elif PyInt_CheckExact(right):
                c = mpz_cmp_si((<Integer>left).value, PyInt_AS_LONG(right))
            elif PyLong_CheckExact(right):
                mpz_set_pylong(mpz_tmp, right)
                c = mpz_cmp((<Integer>left).value, mpz_tmp)
            elif PyFloat_CheckExact(right):
                d = PyFloat_AsDouble(right)
                if isnan(d): return op == 3
                c = mpz_cmp_d((<Integer>left).value, d)
            else:
                return (<sage.structure.element.Element>left)._richcmp(right, op)

        else: # right is an Integer
            if PyInt_CheckExact(left):
                c = -mpz_cmp_si((<Integer>right).value, PyInt_AS_LONG(left))
            elif PyLong_CheckExact(left):
                mpz_set_pylong(mpz_tmp, left)
                c = mpz_cmp(mpz_tmp, (<Integer>right).value)
            elif PyFloat_CheckExact(left):
                d = PyFloat_AsDouble(left)
                if isnan(d): return op == 3
                c = -mpz_cmp_d((<Integer>right).value, d)
            else:
                return (<sage.structure.element.Element>left)._richcmp(right, op)
        return (<sage.structure.element.Element>left)._rich_to_bool(op, c)

    cdef int _cmp_c_impl(left, sage.structure.element.Element right) except -2:
        cdef int i
        i = mpz_cmp((<Integer>left).value, (<Integer>right).value)
        if i < 0: return -1
        elif i == 0: return 0
        else: return 1

    def __copy__(self):
        """
        Return a copy of the integer.

        EXAMPLES::

            sage: n = 2
            sage: copy(n)
            2
            sage: copy(n) is n
            False
        """
        cdef Integer z
        z = PY_NEW(Integer)
        set_mpz(z,self.value)
        return z

    def list(self):
        """
        Return a list with this integer in it, to be compatible with the
        method for number fields.

        EXAMPLES::

            sage: m = 5
            sage: m.list()
            [5]
        """
        return [ self ]

    def  __dealloc__(self):
        mpz_clear(self.value)

    def __repr__(self):
        """
        Return string representation of this integer.

        EXAMPLES::

            sage: n = -5; n.__repr__()
            '-5'
        """
        return self.str()

    def _latex_(self):
        """
        Return latex representation of this integer. This is just the
        underlying string representation and nothing more. This is called
        by the latex function.

        EXAMPLES::

            sage: n = -5; n._latex_()
            '-5'
            sage: latex(n)
            -5
        """
        return self.str()

    def _sympy_(self):
        """
        Convert Sage Integer() to SymPy Integer.

        EXAMPLES::

            sage: n = 5; n._sympy_()
            5
            sage: n = -5; n._sympy_()
            -5
        """
        import sympy
        return sympy.sympify(int(self))

    def _mathml_(self):
        """
        Return mathml representation of this integer.

        EXAMPLES::

            sage: mathml(-45)
            <mn>-45</mn>
            sage: (-45)._mathml_()
            '<mn>-45</mn>'
        """
        return '<mn>%s</mn>'%self

    def str(self, int base=10):
        r"""
        Return the string representation of ``self`` in the
        given base.

        EXAMPLES::

            sage: Integer(2^10).str(2)
            '10000000000'
            sage: Integer(2^10).str(17)
            '394'

        ::

            sage: two=Integer(2)
            sage: two.str(1)
            Traceback (most recent call last):
            ...
            ValueError: base (=1) must be between 2 and 36

        ::

            sage: two.str(37)
            Traceback (most recent call last):
            ...
            ValueError: base (=37) must be between 2 and 36

        ::

            sage: big = 10^5000000
            sage: s = big.str()       # long time (> 20 seconds)
            sage: len(s)              # long time (depends on above defn of s)
            5000001
            sage: s[:10]              # long time (depends on above defn of s)
            '1000000000'
        """
        if base < 2 or base > 36:
            raise ValueError, "base (=%s) must be between 2 and 36"%base
        cdef size_t n
        cdef char *s
        n = mpz_sizeinbase(self.value, base) + 2
        s = <char *>PyMem_Malloc(n)
        if s == NULL:
            raise MemoryError, "Unable to allocate enough memory for the string representation of an integer."
        sig_on()
        mpz_get_str(s, base, self.value)
        sig_off()
        k = <object> PyString_FromString(s)
        PyMem_Free(s)
        return k

    def __format__(self, *args, **kwargs):
        """
        Returns a string representation using Python's Format protocol.
        Valid format descriptions are exactly those for Python integers.

        EXAMPLES::

            sage: "{0:#x}; {0:#b}; {0:+05d}".format(ZZ(17))
            '0x11; 0b10001; +0017'

        """
        return int(self).__format__(*args,**kwargs)

    def ordinal_str(self):
        """
        Returns a string representation of the ordinal associated to self.

        EXAMPLES::

            sage: [ZZ(n).ordinal_str() for n in range(25)]
            ['0th',
            '1st',
            '2nd',
            '3rd',
            '4th',
            ...
            '10th',
            '11th',
            '12th',
            '13th',
            '14th',
            ...
            '20th',
            '21st',
            '22nd',
            '23rd',
            '24th']

            sage: ZZ(1001).ordinal_str()
            '1001st'

            sage: ZZ(113).ordinal_str()
            '113th'
            sage: ZZ(112).ordinal_str()
            '112th'
            sage: ZZ(111).ordinal_str()
            '111th'

        """
        if self<0:
            raise ValueError, "Negative integers are not ordinals."
        n = self.abs()
        if ((n%100)!=11 and n%10==1):
            th = 'st'
        elif ((n%100)!=12 and n%10==2):
            th = 'nd'
        elif ((n%100)!=13 and n%10==3):
            th = 'rd'
        else:
            th = 'th'
        return n.str()+th

    def __hex__(self):
        r"""
        Return the hexadecimal digits of self in lower case.

        .. note::

           '0x' is *not* prepended to the result like is done by the
           corresponding Python function on int or long. This is for
           efficiency sake--adding and stripping the string wastes
           time; since this function is used for conversions from
           integers to other C-library structures, it is important
           that it be fast.

        EXAMPLES::

            sage: print hex(Integer(15))
            f
            sage: print hex(Integer(16))
            10
            sage: print hex(Integer(16938402384092843092843098243))
            36bb1e3929d1a8fe2802f083
            sage: print hex(long(16938402384092843092843098243))
            0x36bb1e3929d1a8fe2802f083L
        """
        return self.str(16)

    def __oct__(self):
        r"""
        Return the digits of self in base 8.

        .. note::

           '0' is *not* prepended to the result like is done by the
           corresponding Python function on int or long. This is for
           efficiency sake--adding and stripping the string wastes
           time; since this function is used for conversions from
           integers to other C-library structures, it is important
           that it be fast.

        EXAMPLES::

            sage: print oct(Integer(800))
            1440
            sage: print oct(Integer(8))
            10
            sage: print oct(Integer(-50))
            -62
            sage: print oct(Integer(-899))
            -1603
            sage: print oct(Integer(16938402384092843092843098243))
            15535436162247215217705000570203

        Behavior of Sage integers vs. Python integers::

            sage: oct(Integer(10))
            '12'
            sage: oct(int(10))
            '012'
            sage: oct(Integer(-23))
            '-27'
            sage: oct(int(-23))
            '-027'
        """
        return self.str(8)

    def binary(self):
        """
        Return the binary digits of self as a string.

        EXAMPLES::

            sage: print Integer(15).binary()
            1111
            sage: print Integer(16).binary()
            10000
            sage: print Integer(16938402384092843092843098243).binary()
            1101101011101100011110001110010010100111010001101010001111111000101000000000101111000010000011
        """
        return self.str(2)

    def bits(self):
        """
        Return the bits in self as a list, least significant first. The
        result satisfies the identity

        ::

            x == sum(b*2^e for e, b in enumerate(x.bits()))

        Negative numbers will have negative "bits". (So, strictly
        speaking, the entries of the returned list are not really
        members of `\ZZ/2\ZZ`.)

        This method just calls :func:`digits` with ``base=2``.

        SEE ALSO:

        :func:`nbits` (number of bits; a faster way to compute
        ``len(x.bits())``; and :func:`binary`, which returns a string in
        more-familiar notation.

        EXAMPLES::

            sage: 500.bits()
            [0, 0, 1, 0, 1, 1, 1, 1, 1]
            sage: 11.bits()
            [1, 1, 0, 1]
            sage: (-99).bits()
            [-1, -1, 0, 0, 0, -1, -1]
        """
        return self.digits(base=2)

    def nbits(self):
        """
        Return the number of bits in self.

        EXAMPLES::

            sage: 500.nbits()
            9
            sage: 5.nbits()
            3
            sage: 0.nbits() == len(0.bits()) == 0.ndigits(base=2)
            True
            sage: 12345.nbits() == len(12345.binary())
            True
        """
        # mpz_sizeinbase(0,2) always returns 1
        if mpz_cmp_si(self.value,0) == 0:
           return int(0)
        else:
           return int(mpz_sizeinbase(self.value, 2))

    def trailing_zero_bits(self):
        """
        Return the number of trailing zero bits in self, i.e.
        the exponent of the largest power of 2 dividing self.

        EXAMPLES::

            sage: 11.trailing_zero_bits()
            0
            sage: (-11).trailing_zero_bits()
            0
            sage: (11<<5).trailing_zero_bits()
            5
            sage: (-11<<5).trailing_zero_bits()
            5
            sage: 0.trailing_zero_bits()
            0

        """
        if mpz_sgn(self.value) == 0:
            return int(0)
        return int(mpz_scan1(self.value, 0))

    def digits(self, base=10, digits=None, padto=0):
        r"""
        Return a list of digits for ``self`` in the given base in little
        endian order.

        The returned value is unspecified if self is a negative number
        and the digits are given.

        INPUT:

        -  ``base`` - integer (default: 10)

        -  ``digits`` - optional indexable object as source for
           the digits

        -  ``padto`` - the minimal length of the returned list,
           sufficient number of zeros are added to make the list minimum that
           length (default: 0)

        As a shorthand for ``digits(2)``, you can use :meth:`.bits`.

        Also see :meth:`ndigits`.

        EXAMPLES::

            sage: 17.digits()
            [7, 1]
            sage: 5.digits(base=2, digits=["zero","one"])
            ['one', 'zero', 'one']
            sage: 5.digits(3)
            [2, 1]
            sage: 0.digits(base=10)  # 0 has 0 digits
            []
            sage: 0.digits(base=2)  # 0 has 0 digits
            []
            sage: 10.digits(16,'0123456789abcdef')
            ['a']
            sage: 0.digits(16,'0123456789abcdef')
            []
            sage: 0.digits(16,'0123456789abcdef',padto=1)
            ['0']
            sage: 123.digits(base=10,padto=5)
            [3, 2, 1, 0, 0]
            sage: 123.digits(base=2,padto=3)       # padto is the minimal length
            [1, 1, 0, 1, 1, 1, 1]
            sage: 123.digits(base=2,padto=10,digits=(1,-1))
            [-1, -1, 1, -1, -1, -1, -1, 1, 1, 1]
            sage: a=9939082340; a.digits(10)
            [0, 4, 3, 2, 8, 0, 9, 3, 9, 9]
            sage: a.digits(512)
            [100, 302, 26, 74]
            sage: (-12).digits(10)
            [-2, -1]
            sage: (-12).digits(2)
            [0, 0, -1, -1]

        We support large bases.

        ::

            sage: n=2^6000
            sage: n.digits(2^3000)
            [0, 0, 1]

        ::

            sage: base=3; n=25
            sage: l=n.digits(base)
            sage: # the next relationship should hold for all n,base
            sage: sum(base^i*l[i] for i in range(len(l)))==n
            True
            sage: base=3; n=-30; l=n.digits(base); sum(base^i*l[i] for i in range(len(l)))==n
            True

        The inverse of this method -- constructing an integer from a
        list of digits and a base -- can be done using the above method
        or by simply using :class:`ZZ()
        <sage.rings.integer_ring.IntegerRing_class>` with a base::

            sage: x = 123; ZZ(x.digits(), 10)
            123
            sage: x == ZZ(x.digits(6), 6)
            True
            sage: x == ZZ(x.digits(25), 25)
            True

        Using :func:`sum` and :func:`enumerate` to do the same thing is
        slightly faster in many cases (and
        :func:`~sage.misc.misc_c.balanced_sum` may be faster yet). Of
        course it gives the same result::

            sage: base = 4
            sage: sum(digit * base^i for i, digit in enumerate(x.digits(base))) == ZZ(x.digits(base), base)
            True

        Note: In some cases it is faster to give a digits collection. This
        would be particularly true for computing the digits of a series of
        small numbers. In these cases, the code is careful to allocate as
        few python objects as reasonably possible.

        ::

            sage: digits = range(15)
            sage: l=[ZZ(i).digits(15,digits) for i in range(100)]
            sage: l[16]
            [1, 1]

        This function is comparable to ``str`` for speed.

        ::

            sage: n=3^100000
            sage: n.digits(base=10)[-1]  # slightly slower than str
            1
            sage: n=10^10000
            sage: n.digits(base=10)[-1]  # slightly faster than str
            1

        AUTHORS:

        - Joel B. Mohler (2008-03-02):  significantly rewrote this entire function
        """
        cdef Integer _base
        cdef Integer self_abs = self
        cdef int power_index = 0
        cdef list power_list
        cdef list l
        cdef int i
        cdef size_t s

        if PY_TYPE_CHECK(base, Integer):
            _base = <Integer>base
        else:
            _base = Integer(base)

        if mpz_cmp_si(_base.value,2) < 0:
            raise ValueError, "base must be >= 2"

        if mpz_sgn(self.value) < 0:
            self_abs = -self

        cdef bint do_sig_on
        if mpz_sgn(self.value) == 0:
            l = [the_integer_ring._zero_element if digits is None else digits[0]]*padto
        elif mpz_cmp_si(_base.value,2) == 0:
            s = mpz_sizeinbase(self.value, 2)
            if digits:
                o = digits[1]
                z = digits[0]
            else:
                if mpz_sgn(self.value) == 1:
                    o = the_integer_ring._one_element
                else:
                    o = -the_integer_ring._one_element
                z = the_integer_ring._zero_element
            l = [z]*(s if s >= padto else padto)
            for i from 0<= i < s:
                # mpz_tstbit seems to return 0 for the high-order bit of
                # negative numbers?!
                if mpz_tstbit(self_abs.value,i):
                    l[i] = o
        else:
            s = mpz_sizeinbase(self.value, 2)
            do_sig_on = (s > 256)
            if do_sig_on: sig_on()

            # We use a divide and conquer approach (suggested by the prior
            # author, malb?, of the digits method) here: for base b, compute
            # b^2, b^4, b^8, ... (repeated squaring) until you get larger
            # than your number; then compute (n // b^256, n % b^256)
            # (if b^512 > number) to split the number in half and recurse

            # Pre-computing the exact number of digits up-front is actually
            # faster (especially for large values of self) than trimming off
            # trailing zeros after the fact.  It also seems that it would
            # avoid duplicating the list in memory with a list-slice.
            z = the_integer_ring._zero_element if digits is None else digits[0]
            s = self_abs.exact_log(_base)
            l = [z]*(s+1 if s+1 >= padto else padto)

            # set up digits for optimal access once we get inside the worker
            # functions
            if not digits is None:
                # list objects have fastest access in the innermost loop
                if not PyList_CheckExact(digits):
                    digits = [digits[i] for i in range(_base)]
            elif mpz_cmp_ui(_base.value,s) < 0 and mpz_cmp_ui(_base.value,10000):
                # We can get a speed boost by pre-allocating digit values in
                # big cases.
                # We do this we have more digits than the base and the base
                # is not too extremely large (currently, "extremely" means
                # larger than 10000 -- that's very arbitrary.)
                if mpz_sgn(self.value) > 0:
                    digits = [Integer(i) for i in range(_base)]
                else:
                    # All the digits will be negated in the recursive function.
                    # we'll just compensate for python index semantics
                    digits = [Integer(i) for i in range(-_base,0)]
                    digits[0] = the_integer_ring._zero_element

            if s < 40:
                _digits_naive(self.value,l,0,_base,digits)
            else:
                # count the bits of s
                i = 0
                while s != 0:
                    s >>= 1
                    i += 1

                power_list = [_base]*i
                for power_index from 1 <= power_index < i:
                    power_list[power_index] = power_list[power_index-1]**2

                # Note that it may appear that the recursive calls to
                # _digit_internal would be assigning list elements i in l for
                # anywhere from 0<=i<(1<<power_index).  However, this is not
                # the case due to the optimization of skipping assigns
                # assigning zero.
                _digits_internal(self.value,l,0,i-1,power_list,digits)

            if do_sig_on: sig_off()

        # padding should be taken care of with-in the function
        # all we need to do is return
        return l

    def ndigits(self, base=10):
        """
        Return the number of digits of self expressed in the given base.

        INPUT:

        -  ``base`` - integer (default: 10)

        EXAMPLES::

            sage: n = 52
            sage: n.ndigits()
            2
            sage: n = -10003
            sage: n.ndigits()
            5
            sage: n = 15
            sage: n.ndigits(2)
            4
            sage: n = 1000**1000000+1
            sage: n.ndigits()
            3000001
            sage: n = 1000**1000000-1
            sage: n.ndigits()
            3000000
            sage: n = 10**10000000-10**9999990
            sage: n.ndigits()
            10000000
        """
        cdef Integer temp

        if mpz_sgn(self.value) == 0:
            temp = PY_NEW(Integer)
            mpz_set_ui(temp.value, 0)
            return temp

        if mpz_sgn(self.value) > 0:
            temp = self.exact_log(base)
            mpz_add_ui(temp.value, temp.value, 1)
            return temp
        else:
            return self.abs().exact_log(base) + 1

    cdef void set_from_mpz(Integer self, mpz_t value):
        mpz_set(self.value, value)

    cdef mpz_t* get_value(Integer self):
        return &self.value

    cdef void _to_ZZ(self, ZZ_c *z):
        sig_on()
        mpz_to_ZZ(z, &self.value)
        sig_off()

    cpdef ModuleElement _add_(self, ModuleElement right):
        """
        Integer addition.

        TESTS::

            sage: Integer(32) + Integer(23)
            55
            sage: sum(Integer(i) for i in [1..100])
            5050
            sage: a = ZZ.random_element(10^50000)
            sage: b = ZZ.random_element(10^50000)
            sage: a+b == b+a
            True
        """
        # self and right are guaranteed to be Integers
        cdef Integer x = <Integer>PY_NEW(Integer)
        mpz_add(x.value, self.value, (<Integer>right).value)
        return x

    cpdef ModuleElement _iadd_(self, ModuleElement right):
        # self and right are guaranteed to be Integers, self safe to mutate
        mpz_add(self.value, self.value, (<Integer>right).value)
        return self

    cdef RingElement _add_long(self, long n):
        """
        Fast path for adding a C long.

        TESTS::
            sage: int(10) + Integer(100)
            110
            sage: Integer(100) + int(10)
            110
            sage: Integer(10^100) + int(10)
            10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010

        Also called for subtraction::

            sage: Integer(100) - int(10)
            90
            sage: Integer(10^100) - int(10)
            9999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999990

        Make sure it works when -<long>n would overflow::

            sage: most_neg_long = int(-sys.maxint - 1)
            sage: type(most_neg_long), type(-most_neg_long)
            (<type 'int'>, <type 'long'>)
            sage: 0 + most_neg_long == most_neg_long
            True
            sage: 0 - most_neg_long == -most_neg_long
            True
        """
        cdef Integer x = <Integer>PY_NEW(Integer)
        if n > 0:
            mpz_add_ui(x.value, self.value, n)
        else:
            # Note that 0-<unsigned long>n is always -n as an unsigned
            # long (whereas -n may overflow).
            mpz_sub_ui(x.value, self.value, 0 - <unsigned long>n)
        return x

    cpdef ModuleElement _sub_(self, ModuleElement right):
        """
        Integer subtraction.

        TESTS::

            sage: Integer(32) - Integer(23)
            9
            sage: Integer(10^100) - Integer(1)
            9999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999
            sage: Integer(1) - Integer(10^100)
            -9999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999
            sage: a = ZZ.random_element(10^50000)
            sage: b = ZZ.random_element(10^50000)
            sage: a-b == -(b-a) == a + -b
            True
        """
        # self and right are guaranteed to be Integers
        cdef Integer x = <Integer>PY_NEW(Integer)
        mpz_sub(x.value, self.value, (<Integer>right).value)
        return x

    cpdef ModuleElement _isub_(self, ModuleElement right):
        mpz_sub(self.value, self.value, (<Integer>right).value)
        return self

    def __neg__(self):
        """
        TESTS::

            sage: a = Integer(3)
            sage: -a
            -3
            sage: a = Integer(3^100); a
            515377520732011331036461129765621272702107522001
            sage: -a
            -515377520732011331036461129765621272702107522001
        """
        cdef Integer x = <Integer>PY_NEW(Integer)
        mpz_neg(x.value, self.value)
        return x

    cpdef ModuleElement _neg_(self):
        cdef Integer x = <Integer>PY_NEW(Integer)
        mpz_neg(x.value, self.value)
        return x

    cpdef _act_on_(self, s, bint self_on_left):
        """
        EXAMPLES::

            sage: 8 * [0] #indirect doctest
            [0, 0, 0, 0, 0, 0, 0, 0]
            sage: 'hi' * 8
            'hihihihihihihihi'
        """
        if isinstance(s, (list, tuple, basestring)):
            if mpz_fits_slong_p(self.value):
                return s * mpz_get_si(self.value)
            else:
                return s * int(self) # will raise the appropriate exception

    cdef ModuleElement _mul_long(self, long n):
        """
        Fast path for multiplying a C long.

        TESTS::

            sage: Integer(25) * int(4)
            100
            sage: int(4) * Integer(25)
            100
            sage: Integer(10^100) * int(4)
            40000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
        """
        cdef Integer x = <Integer>PY_NEW(Integer)
        if mpz_size(self.value) > 100000:
            sig_on()
            mpz_mul_si(x.value, self.value, n)
            sig_off()
        else:
            mpz_mul_si(x.value, self.value, n)
        return x

    cpdef RingElement _mul_(self, RingElement right):
        """
        Integer multiplication.

            sage: Integer(25) * Integer(4)
            100
            sage: Integer(5^100) * Integer(2^100)
            10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
            sage: a = ZZ.random_element(10^50000)
            sage: b = ZZ.random_element(10^50000)
            sage: a*b == b*a
            True
        """
        # self and right are guaranteed to be Integers
        cdef Integer x = <Integer>PY_NEW(Integer)
        if mpz_size(self.value) + mpz_size((<Integer>right).value) > 100000:
            # We only use the signal handler (to enable ctrl-c out) when the
            # product might take a while to compute
            sig_on()
            mpz_mul(x.value, self.value, (<Integer>right).value)
            sig_off()
        else:
            mpz_mul(x.value, self.value, (<Integer>right).value)
        return x

    cpdef RingElement _imul_(self, RingElement right):
        if mpz_size(self.value) + mpz_size((<Integer>right).value) > 100000:
            # We only use the signal handler (to enable ctrl-c out) when the
            # product might take a while to compute
            sig_on()
            mpz_mul(self.value, self.value, (<Integer>right).value)
            sig_off()
        else:
            mpz_mul(self.value, self.value, (<Integer>right).value)
        return self

    cpdef RingElement _div_(self, RingElement right):
        r"""
        Computes `\frac{a}{b}`

        EXAMPLES::

            sage: a = Integer(3) ; b = Integer(4)
            sage: a / b == Rational(3) / 4
            True
            sage: Integer(32) / Integer(32)
            1
        """
        # This is vastly faster than doing it here, since here
        # we can't cimport rationals.
        return the_integer_ring._div(self, right)

    def __floordiv__(x, y):
        r"""
        Computes the whole part of `\frac{x}{y}`.

        EXAMPLES::

            sage: a = Integer(321) ; b = Integer(10)
            sage: a // b
            32
            sage: z = Integer(-231)
            sage: z // 2
            -116
            sage: z = Integer(231)
            sage: z // 2
            115
            sage: z // -2
            -116
            sage: z // 0
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Integer division by zero
            sage: 101 // int(5)
            20
            sage: 100 // int(-3)
            -34

        TESTS::

            sage: signs = [(11,5), (11,-5), (-11,5), (-11,-5)]
            sage: control = [int(a) // int(b) for a, b in signs]
            sage: [a // b for a,b in signs] == control
            True
            sage: [a // int(b) for a,b in signs] == control
            True
            sage: [int(a) // b for a,b in signs] == control
            True
        """
        cdef Integer z = <Integer>PY_NEW(Integer)
        cdef long yy, res
        if PY_TYPE(x) is PY_TYPE(y):
            if not mpz_sgn((<Integer>y).value):
                raise ZeroDivisionError, "Integer division by zero"
            if mpz_size((<Integer>x).value) > 100000:
                sig_on()
                mpz_fdiv_q(z.value, (<Integer>x).value, (<Integer>y).value)
                sig_off()
            else:
                mpz_fdiv_q(z.value, (<Integer>x).value, (<Integer>y).value)
            return z

        elif PyInt_CheckExact(y):
            yy = PyInt_AS_LONG(y)
            if yy > 0:
                mpz_fdiv_q_ui(z.value, (<Integer>x).value, yy)
            elif yy == 0:
                raise ZeroDivisionError, "Integer division by zero"
            else:
                res = mpz_fdiv_q_ui(z.value, (<Integer>x).value, -yy)
                mpz_neg(z.value, z.value)
                if res:
                    mpz_sub_ui(z.value, z.value, 1)
            return z

        else:
            return bin_op(x, y, operator.floordiv)

    def __pow__(self, n, modulus):
        r"""
        Computes `\text{self}^n`

        EXAMPLES::

            sage: 2^-6
            1/64
            sage: 2^6
            64
            sage: 2^0
            1
            sage: 2^-0
            1
            sage: (-1)^(1/3)
            (-1)^(1/3)

        For consistency with Python and MPFR, 0^0 is defined to be 1 in
        Sage::

            sage: 0^0
            1

        The base need not be an integer (it can be a builtin Python type).

        ::

            sage: int(2)^10
            1024
            sage: float(2.5)^10
            9536.7431640625
            sage: 'sage'^3
            'sagesagesage'

        The exponent must fit in a long unless the base is -1, 0, or 1.

        ::

            sage: x = 2^100000000000000000000000
            Traceback (most recent call last):
            ...
            RuntimeError: exponent must be at most 2147483647  # 32-bit
            RuntimeError: exponent must be at most 9223372036854775807 # 64-bit
            sage: (-1)^100000000000000000000000
            1

        We raise 2 to various interesting exponents::

            sage: 2^x                # symbolic x
            2^x
            sage: 2^1.5              # real number
            2.82842712474619
            sage: 2^float(1.5)       # python float
            2.8284271247461903
            sage: 2^I                # complex number
            2^I
            sage: f = 2^(sin(x)-cos(x)); f
            2^(-cos(x) + sin(x))
            sage: f(x=3)
            2^(-cos(3) + sin(3))

        A symbolic sum::

            sage: x,y,z = var('x,y,z')
            sage: 2^(x+y+z)
            2^(x + y + z)
            sage: 2^(1/2)
            sqrt(2)
            sage: 2^(-1/2)
            1/2*sqrt(2)

        TESTS::

            sage: complex(0,1)^2
            (-1+0j)
            sage: R.<t> = QQ[]
            sage: 2^t
            Traceback (most recent call last):
            ...
            TypeError: non-integral exponents not supported
            sage: int(3)^-3
            1/27
            sage: type(int(3)^2)
            <type 'sage.rings.integer.Integer'>
            sage: type(int(3)^int(2))
            <type 'int'>

        Check that a ``MemoryError`` is thrown if the resulting number
        would be ridiculously large, see :trac:`15363`::

            sage: 2^(2^63-2)
            Traceback (most recent call last):
            ...
            RuntimeError: exponent must be at most 2147483647          # 32-bit
            MemoryError: failed to allocate 1152921504606847008 bytes  # 64-bit
        """
        if modulus is not None:
            from sage.rings.finite_rings.integer_mod import Mod
            return Mod(self, modulus) ** n

        if not PY_TYPE_CHECK(self, Integer):
            if isinstance(self, str):
                return self * n
            if not PY_TYPE_CHECK(self, int):
                return self ** int(n)
            else:
                self = Integer(self)        #convert from int to Integer
        cdef Integer _self = <Integer>self
        cdef long nn

        try:
            nn = PyNumber_Index(n)
        except TypeError:
            s = parent_c(n)(self)
            return s**n
        except OverflowError:
            if mpz_cmp_si(_self.value, 1) == 0:
                return self
            elif mpz_cmp_si(_self.value, 0) == 0:
                return self
            elif mpz_cmp_si(_self.value, -1) == 0:
                return self if n % 2 else -self
            raise RuntimeError, "exponent must be at most %s" % sys.maxint

        if nn == 0:
            return one

        cdef Integer x = PY_NEW(Integer)

        sig_on()
        mpz_pow_ui(x.value, (<Integer>self).value, nn if nn > 0 else -nn)
        sig_off()

        if nn < 0:
            return ~x
        else:
            return x

    def nth_root(self, int n, bint truncate_mode=0):
        r"""
        Returns the (possibly truncated) n'th root of self.

        INPUT:

        -  ``n`` - integer >= 1 (must fit in C int type).

        -  ``truncate_mode`` - boolean, whether to allow truncation if
           self is not an n'th power.

        OUTPUT:

        If truncate_mode is 0 (default), then returns the exact n'th root
        if self is an n'th power, or raises a ValueError if it is not.

        If truncate_mode is 1, then if either n is odd or self is
        positive, returns a pair (root, exact_flag) where root is the
        truncated nth root (rounded towards zero) and exact_flag is a
        boolean indicating whether the root extraction was exact;
        otherwise raises a ValueError.

        AUTHORS:

        - David Harvey (2006-09-15)
        - Interface changed by John Cremona (2009-04-04)

        EXAMPLES::

            sage: Integer(125).nth_root(3)
            5
            sage: Integer(124).nth_root(3)
            Traceback (most recent call last):
            ...
            ValueError: 124 is not a 3rd power
            sage: Integer(124).nth_root(3, truncate_mode=1)
            (4, False)
            sage: Integer(125).nth_root(3, truncate_mode=1)
            (5, True)
            sage: Integer(126).nth_root(3, truncate_mode=1)
            (5, False)

        ::

            sage: Integer(-125).nth_root(3)
            -5
            sage: Integer(-125).nth_root(3,truncate_mode=1)
            (-5, True)
            sage: Integer(-124).nth_root(3,truncate_mode=1)
            (-4, False)
            sage: Integer(-126).nth_root(3,truncate_mode=1)
            (-5, False)

        ::

            sage: Integer(125).nth_root(2, True)
            (11, False)
            sage: Integer(125).nth_root(3, True)
            (5, True)

        ::

            sage: Integer(125).nth_root(-5)
            Traceback (most recent call last):
            ...
            ValueError: n (=-5) must be positive

        ::

            sage: Integer(-25).nth_root(2)
            Traceback (most recent call last):
            ...
            ValueError: cannot take even root of negative number

        ::

            sage: a=9
            sage: a.nth_root(3)
            Traceback (most recent call last):
            ...
            ValueError: 9 is not a 3rd power

            sage: a.nth_root(22)
            Traceback (most recent call last):
            ...
            ValueError: 9 is not a 22nd power

            sage: ZZ(2^20).nth_root(21)
            Traceback (most recent call last):
            ...
            ValueError: 1048576 is not a 21st power

            sage: ZZ(2^20).nth_root(21, truncate_mode=1)
            (1, False)

        """
        if n < 1:
            raise ValueError, "n (=%s) must be positive" % n
        if (mpz_sgn(self.value) < 0) and not (n & 1):
            raise ValueError, "cannot take even root of negative number"
        cdef Integer x
        cdef bint is_exact
        x = PY_NEW(Integer)
        sig_on()
        is_exact = mpz_root(x.value, self.value, n)
        sig_off()

        if truncate_mode:
            return x, is_exact
        else:
            if is_exact:
                return x
            else:
                raise ValueError, "%s is not a %s power"%(self,integer_ring.ZZ(n).ordinal_str())

    cpdef size_t _exact_log_log2_iter(self,Integer m):
        """
        This is only for internal use only.  You should expect it to crash
        and burn for negative or other malformed input.  In particular, if
        the base `2 \leq m < 4` the log2 approximation of m is 1 and certain
        input causes endless loops.  Along these lines, it is clear that
        this function is most useful for m with a relatively large number
        of bits.

        For ``small`` values (which I'll leave quite ambiguous), this function
        is a fast path for exact log computations.  Any integer division with
        such input tends to dominate the runtime.  Thus we avoid division
        entirely in this function.

        AUTHOR::

        - Joel B. Mohler (2009-04-10)

        EXAMPLES::

            sage: Integer(125)._exact_log_log2_iter(4)
            3
            sage: Integer(5^150)._exact_log_log2_iter(5)
            150
        """
        cdef size_t n_log2
        cdef size_t m_log2
        cdef size_t l_min
        cdef size_t l_max
        cdef size_t l
        cdef Integer result
        cdef mpz_t accum
        cdef mpz_t temp_exp

        if mpz_cmp_si(m.value,4) < 0:
            raise ValueError, "This is undefined or possibly non-convergent with this algorithm."

        n_log2=mpz_sizeinbase(self.value,2)-1
        m_log2=mpz_sizeinbase(m.value,2)-1
        l_min=n_log2/(m_log2+1)
        l_max=n_log2/m_log2
        if l_min != l_max:
            sig_on()
            mpz_init(accum)
            mpz_init(temp_exp)
            mpz_set_ui(accum,1)
            l = 0
            while l_min != l_max:
                # print "self=...",m,l_min,l_max
                if l_min + 1 == l_max:
                    mpz_pow_ui(temp_exp,m.value,l_min+1-l)
                    # This might over-shoot and make accum > self, but
                    # we'll know that it's only over by a factor of m^1.
                    mpz_mul(accum,accum,temp_exp)
                    if mpz_cmp(self.value,accum) >= 0:
                        l_min += 1
                    break
                mpz_pow_ui(temp_exp,m.value,l_min-l)
                mpz_mul(accum,accum,temp_exp)
                l = l_min

                # Let x=n_log2-(mpz_sizeinbase(accum,2)-1) and y=m_log2.
                # Now, with x>0 and y>0, we have the following observation.
                # If floor((x-1)/(y+1))=0, then x-1<y+1 which implies that
                # x/y<1+2/y.
                # So long as y>=2, this means that floor(x/y)<=1.  This shows
                # that this iteration is forced to converge for input m >= 4.
                # If m=3, we can find input so that floor((x-1)/(y+1))=0 and
                # floor(x/y)=2 which results in non-convergence.

                # We need the additional '-1' in the l_min computation
                # because mpz_sizeinbase(accum,2)-1 is smaller than the
                # true log_2(accum)
                l_min=l+(n_log2-(mpz_sizeinbase(accum,2)-1)-1)/(m_log2+1)
                l_max=l+(n_log2-(mpz_sizeinbase(accum,2)-1))/m_log2
            mpz_clear(temp_exp)
            mpz_clear(accum)
            sig_off()
        return l_min

    cpdef size_t _exact_log_mpfi_log(self,m):
        """
        This is only for internal use only.  You should expect it to crash
        and burn for negative or other malformed input.

        I avoid using this function until the input is large.  The overhead
        associated with computing the floating point log entirely dominates
        the runtime for small values.  Note that this is most definitely not
        an artifact of format conversion.  Tricks with log2 approximations
        and using exact integer arithmetic are much better for small input.

        AUTHOR::

        - Joel B. Mohler (2009-04-10)

        EXAMPLES::

            sage: Integer(125)._exact_log_mpfi_log(3)
            4
            sage: Integer(5^150)._exact_log_mpfi_log(5)
            150
        """
        cdef int i
        cdef list pow_2_things
        cdef int pow_2
        cdef size_t upper,lower,middle

        import real_mpfi
        R=real_mpfi.RIF

        rif_self = R(self)

        sig_on()
        rif_m = R(m)
        rif_log = rif_self.log()/rif_m.log()
        # upper is *greater* than the answer
        try:
            upper = rif_log.upper().ceiling()
        except StandardError:
            # ceiling is probably Infinity
            # I'm not sure what to do now
            upper = 0
        lower = rif_log.lower().floor()
        # since the log function is monotonic increasing, lower
        # and upper bracket our desired answer

        # if upper - lower == 1: "we are done"
        if upper - lower == 2:
            # You could test it by checking rif_m**(lower+1), but I think
            # that's a waste of time since it won't be conclusive.
            # We must test with exact integer arithmetic which takes all
            # the bits of self into account.
            sig_off()
            if self >= m**(lower+1):
                return lower + 1
            else:
                return lower
        elif upper - lower > 2:
            # this case would only happen in cases with extremely large 'self'
            rif_m = R(m)
            min_power = rif_m**lower
            middle = upper-lower
            pow_2 = 0
            while middle != 0:
                middle >>= 1
                pow_2 += 1
            # if middle was an exact power of 2, adjust down
            if (1 << (pow_2-1)) == upper-lower:
                pow_2 -= 1
            #print upper, lower, pow_2
            pow_2_things = [rif_m]*pow_2
            for i from 1<=i<pow_2:
                pow_2_things[i] = pow_2_things[i-1]**2
            for i from pow_2>i>=0:
                middle = lower + int(2)**i
                #print "Upper:  %i;  Lower:  %i;  Middle:  %i" % (upper,lower, middle)
                exp = min_power*pow_2_things[i]
                if exp > rif_self:
                    upper = middle
                elif exp < rif_self:
                    lower = middle
                    min_power = exp
                else:
                    sig_off()
                    if m**middle <= self:
                        return middle
                    else:
                        return lower
        sig_off()

        if upper == 0:
            raise ValueError, "The input for exact_log is too large and support is not implemented."

        return lower

    def exact_log(self, m):
        r"""
        Returns the largest integer `k` such that `m^k \leq \text{self}`,
        i.e., the floor of `\log_m(\text{self})`.

        This is guaranteed to return the correct answer even when the usual
        log function doesn't have sufficient precision.

        INPUT:

        -  ``m`` - integer >= 2

        AUTHORS:

        - David Harvey (2006-09-15)
        - Joel B. Mohler (2009-04-08) -- rewrote this to handle small cases
               and/or easy cases up to 100x faster..

        EXAMPLES::

            sage: Integer(125).exact_log(5)
            3
            sage: Integer(124).exact_log(5)
            2
            sage: Integer(126).exact_log(5)
            3
            sage: Integer(3).exact_log(5)
            0
            sage: Integer(1).exact_log(5)
            0
            sage: Integer(178^1700).exact_log(178)
            1700
            sage: Integer(178^1700-1).exact_log(178)
            1699
            sage: Integer(178^1700+1).exact_log(178)
            1700
            sage: # we need to exercise the large base code path too
            sage: Integer(1780^1700-1).exact_log(1780)
            1699

            sage: # The following are very very fast.
            sage: # Note that for base m a perfect power of 2, we get the exact log by counting bits.
            sage: n=2983579823750185701375109835; m=32
            sage: n.exact_log(m)
            18
            sage: # The next is a favorite of mine.  The log2 approximate is exact and immediately provable.
            sage: n=90153710570912709517902579010793251709257901270941709247901209742124;m=213509721309572
            sage: n.exact_log(m)
            4

        ::

            sage: x = 3^100000
            sage: RR(log(RR(x), 3))
            100000.000000000
            sage: RR(log(RR(x + 100000), 3))
            100000.000000000

        ::

            sage: x.exact_log(3)
            100000
            sage: (x+1).exact_log(3)
            100000
            sage: (x-1).exact_log(3)
            99999

        ::

            sage: x.exact_log(2.5)
            Traceback (most recent call last):
            ...
            TypeError: Attempt to coerce non-integral RealNumber to Integer
        """
        cdef Integer _m
        cdef Integer result
        cdef size_t n_log2
        cdef size_t m_log2
        cdef size_t guess # this will contain the final answer
        cdef bint guess_filled = 0  # this variable is only used in one branch below
        cdef mpz_t z
        if PY_TYPE_CHECK(m,Integer):
            _m=<Integer>m
        else:
            _m=<Integer>Integer(m)

        if mpz_sgn(self.value) <= 0 or mpz_sgn(_m.value) <= 0:
            raise ValueError, "both self and m must be positive"
        if mpz_cmp_si(_m.value,2) < 0:
            raise ValueError, "m must be at least 2"

        n_log2=mpz_sizeinbase(self.value,2)-1
        m_log2=mpz_sizeinbase(_m.value,2)-1
        if mpz_divisible_2exp_p(_m.value,m_log2):
            # Here, m is a power of 2 and the correct answer is found
            # by a log 2 approximation.
            guess = n_log2/m_log2  # truncating division
        elif n_log2/(m_log2+1) == n_log2/m_log2:
            # In this case, we have an upper bound and lower bound which
            # give the same answer, thus, the correct answer.
            guess = n_log2/m_log2
        elif m_log2 < 8:  # i.e. m<256
            # if the base m is at most 256, we can use mpz_sizeinbase
            # to get the following guess which is either the exact
            # log, or 1+ the exact log
            guess = mpz_sizeinbase(self.value, mpz_get_si(_m.value)) - 1

            # we've already excluded the case when m is an exact power of 2

            if n_log2/m_log2 > 8000:
                # If we have a very large number of digits, it can be a nice
                # shortcut to test the guess using interval arithmetic.
                # (suggested by David Harvey and Carl Witty)
                # "for randomly distributed integers, the chance of this
                # interval-based comparison failing is absurdly low"
                import real_mpfi
                approx_compare = real_mpfi.RIF(m)**guess
                if self > approx_compare:
                    guess_filled = 1
                elif self < approx_compare:
                    guess_filled = 1
                    guess =  guess - 1
            if not guess_filled:
                # At this point, either
                #  1)  self is close enough to a perfect power of m that we
                #      need an exact comparison, or
                #  2)  the numbers are small enough that converting to the
                #      interval field is more work than the exact comparison.
                compare = _m**guess
                if self < compare:
                    guess = guess - 1
        elif n_log2 < 5000:
            # for input with small exact log, it's very fast to work in exact
            # integer arithmetic starting from log2 approximations
            guess = self._exact_log_log2_iter(_m)
        else:
            # finally, we are out of easy cases this subroutine uses interval
            # arithmetic to guess and check the exact log.
            guess = self._exact_log_mpfi_log(_m)

        result = PY_NEW(Integer)
        mpz_set_ui(result.value,guess)
        return result

    def log(self, m=None, prec=None):
        r"""
        Returns symbolic log by default, unless the logarithm is exact (for
        an integer base). When precision is given, the RealField
        approximation to that bit precision is used.

        This function is provided primarily so that Sage integers may be
        treated in the same manner as real numbers when convenient. Direct
        use of exact_log is probably best for arithmetic log computation.

        INPUT:

        -  ``m`` - default: natural log base e

        -  ``prec`` - integer (default: None): if None, returns
           symbolic, else to given bits of precision as in RealField

        EXAMPLES::

            sage: Integer(124).log(5)
            log(124)/log(5)
            sage: Integer(124).log(5,100)
            2.9950093311241087454822446806
            sage: Integer(125).log(5)
            3
            sage: Integer(125).log(5,prec=53)
            3.00000000000000
            sage: log(Integer(125))
            log(125)

        For extremely large numbers, this works::

            sage: x = 3^100000
            sage: log(x,3)
            100000

        With the new Pynac symbolic backend, log(x) also
        works in a reasonable amount of time for this x::

            sage: x = 3^100000
            sage: log(x)
            log(1334971414230...5522000001)

        But approximations are probably more useful in this
        case, and work to as high a precision as we desire::

            sage: x.log(3,53) # default precision for RealField
            100000.000000000
            sage: (x+1).log(3,53)
            100000.000000000
            sage: (x+1).log(3,1000)
            100000.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

        We can use non-integer bases, with default e::

            sage: x.log(2.5,prec=53)
            119897.784671579

        We also get logarithms of negative integers, via the
        symbolic ring, using the branch from `-pi` to `pi`::

            sage: log(-1)
            I*pi

        The logarithm of zero is done likewise::

            sage: log(0)
            -Infinity
        """
        if mpz_sgn(self.value) <= 0:
            from sage.symbolic.all import SR
            return SR(self).log()
        if m <= 0 and m != None:
            raise ValueError, "m must be positive"
        if prec:
            from sage.rings.real_mpfr import RealField
            if m is None:
                return RealField(prec)(self).log()
            return RealField(prec)(self).log(m)
        if type(m)==Integer and type(self)==Integer and m**(self.exact_log(m))==self:
            return self.exact_log(m)

        from sage.symbolic.all import SR
        from sage.functions.log import function_log
        if m is None:
            return function_log(self,dont_call_method_on_arg=True)
        return function_log(self,dont_call_method_on_arg=True)/\
                function_log(m,dont_call_method_on_arg=True)

    def exp(self, prec=None):
        r"""
        Returns the exponential function of self as a real number.

        This function is provided only so that Sage integers may be treated
        in the same manner as real numbers when convenient.

        INPUT:


        -  ``prec`` - integer (default: None): if None, returns
           symbolic, else to given bits of precision as in RealField


        EXAMPLES::

            sage: Integer(8).exp()
            e^8
            sage: Integer(8).exp(prec=100)
            2980.9579870417282747435920995
            sage: exp(Integer(8))
            e^8

        For even fairly large numbers, this may not be useful.

        ::

            sage: y=Integer(145^145)
            sage: y.exp()
            e^25024207011349079210459585279553675697932183658421565260323592409432707306554163224876110094014450895759296242775250476115682350821522931225499163750010280453185147546962559031653355159703678703793369785727108337766011928747055351280379806937944746847277089168867282654496776717056860661614337004721164703369140625
            sage: y.exp(prec=53) # default RealField precision
            +infinity
        """
        from sage.functions.all import exp
        res = exp(self, dont_call_method_on_arg=True)
        if prec:
            return res.n(prec=prec)
        return res

    def prime_to_m_part(self, m):
        """
        Returns the prime-to-m part of self, i.e., the largest divisor of
        self that is coprime to m.

        INPUT:

        -  ``m`` - Integer

        OUTPUT: Integer

        EXAMPLES::

            sage: z = 43434
            sage: z.prime_to_m_part(20)
            21717
        """
        late_import()
        return sage.rings.arith.prime_to_m_part(self, m)

    def prime_divisors(self):
        """
        The prime divisors of self, sorted in increasing order. If n is
        negative, we do *not* include -1 among the prime divisors, since
        -1 is not a prime number.

        EXAMPLES::

            sage: a = 1; a.prime_divisors()
            []
            sage: a = 100; a.prime_divisors()
            [2, 5]
            sage: a = -100; a.prime_divisors()
            [2, 5]
            sage: a = 2004; a.prime_divisors()
            [2, 3, 167]
        """
        late_import()
        return sage.rings.arith.prime_divisors(self)

    prime_factors = prime_divisors

    def divisors(self):
        """
        Returns a list of all positive integer divisors of the integer
        self.

        EXAMPLES::

            sage: a = -3; a.divisors()
            [1, 3]
            sage: a = 6; a.divisors()
            [1, 2, 3, 6]
            sage: a = 28; a.divisors()
            [1, 2, 4, 7, 14, 28]
            sage: a = 2^5; a.divisors()
            [1, 2, 4, 8, 16, 32]
            sage: a = 100; a.divisors()
            [1, 2, 4, 5, 10, 20, 25, 50, 100]
            sage: a = 1; a.divisors()
            [1]
            sage: a = 0; a.divisors()
            Traceback (most recent call last):
            ...
            ValueError: n must be nonzero
            sage: a = 2^3 * 3^2 * 17; a.divisors()
            [1, 2, 3, 4, 6, 8, 9, 12, 17, 18, 24, 34, 36, 51, 68, 72, 102, 136, 153, 204, 306, 408, 612, 1224]
            sage: a = odd_part(factorial(31))
            sage: v = a.divisors(); len(v)
            172800
            sage: prod(e+1 for p,e in factor(a))
            172800
            sage: all([t.divides(a) for t in v])
            True

        .. note::

           If one first computes all the divisors and then sorts it,
           the sorting step can easily dominate the runtime. Note,
           however, that (non-negative) multiplication on the left
           preserves relative order. One can leverage this fact to
           keep the list in order as one computes it using a process
           similar to that of the merge sort when adding new elements.
        """
        cdef list all, prev, sorted
        cdef long tip, top
        cdef long i, j, e, ee
        cdef Integer apn, p, pn, z, all_tip

        if not self:
            raise ValueError, "n must be nonzero"
        f = self.factor()

        # All of the declarations below are for optimizing the word-sized
        # case.  Operations are performed in c as far as possible without
        # overflow before moving to python objects.
        cdef long long p_c, pn_c, apn_c
        cdef long all_len, sorted_len, prev_len
        cdef long long* ptr
        cdef long long* empty_c
        cdef long long* swap_tmp
        cdef long long* all_c
        cdef long long* sorted_c
        cdef long long* prev_c

        cdef long divisor_count = 1
        for p,e in f: divisor_count *= (1+e)
        ptr = <long long*>sage_malloc(sizeof(long long) * 3 * divisor_count)
        if ptr == NULL:
            raise MemoryError
        all_c = ptr
        sorted_c = ptr + divisor_count
        prev_c = ptr + (2*divisor_count)

        # These are used to keep track of whether or not we are able to
        # perform the operations in machine words. A factor of two safety
        # margin is added to cover any floating-point rounding issues.
        cdef bint fits_c = True
        cdef double cur_max = 1
        cdef double fits_max = 2.0**(sizeof(long long)*8-2)

        sorted_c[0] = 1
        sorted_len = 1

        for p, e in f:

            cur_max *= (<double>p)**e
            if fits_c and cur_max > fits_max:
                sorted = []
                for i from 0 <= i < sorted_len:
                    z = <Integer>PY_NEW(Integer)
                    mpz_set_longlong(z.value, sorted_c[i])
                    sorted.append(z)
                sage_free(ptr)
                fits_c = False

            # The two cases below are essentially the same algorithm, one
            # operating on Integers in Python lists, the other on long longs.
            if fits_c:

                pn_c = p_c = p

                swap_tmp = sorted_c
                sorted_c = prev_c
                prev_c = swap_tmp
                prev_len = sorted_len
                sorted_len = 0

                tip = 0
                prev_c[prev_len] = prev_c[prev_len-1] * pn_c
                for i from 0 <= i < prev_len:
                    apn_c = prev_c[i] * pn_c
                    while prev_c[tip] < apn_c:
                        sorted_c[sorted_len] = prev_c[tip]
                        sorted_len += 1
                        tip += 1
                    sorted_c[sorted_len] = apn_c
                    sorted_len += 1

                for ee in range(1, e):

                    swap_tmp = all_c
                    all_c = sorted_c
                    sorted_c = swap_tmp
                    all_len = sorted_len
                    sorted_len = 0

                    pn_c *= p_c
                    tip = 0
                    all_c[all_len] = prev_c[prev_len-1] * pn_c
                    for i from 0 <= i < prev_len:
                        apn_c = prev_c[i] * pn_c
                        while all_c[tip] < apn_c:
                            sorted_c[sorted_len] = all_c[tip]
                            sorted_len += 1
                            tip += 1
                        sorted_c[sorted_len] = apn_c
                        sorted_len += 1

            else:
                prev = sorted
                pn = <Integer>PY_NEW(Integer)
                mpz_set_ui(pn.value, 1)
                for ee in range(e):
                    all = sorted
                    sorted = []
                    tip = 0
                    top = len(all)
                    mpz_mul(pn.value, pn.value, p.value) # pn *= p
                    for a in prev:
                        # apn = a*pn
                        apn = <Integer>PY_NEW(Integer)
                        mpz_mul(apn.value, (<Integer>a).value, pn.value)
                        while tip < top:
                            all_tip = <Integer>all[tip]
                            if mpz_cmp(all_tip.value, apn.value) > 0:
                                break
                            sorted.append(all_tip)
                            tip += 1
                        sorted.append(apn)

        if fits_c:
            # all the data is in sorted_c
            sorted = []
            for i from 0 <= i < sorted_len:
                z = <Integer>PY_NEW(Integer)
                mpz_set_longlong(z.value, sorted_c[i])
                sorted.append(z)
            sage_free(ptr)

        return sorted


    def __pos__(self):
        """
        EXAMPLES::

            sage: z=43434
            sage: z.__pos__()
            43434
        """
        return self

    def __abs__(self):
        """
        Computes `|self|`

        EXAMPLES::

            sage: z = -1
            sage: abs(z)
            1
            sage: abs(z) == abs(1)
            True
        """
        cdef Integer x = PY_NEW(Integer)
        mpz_abs(x.value, self.value)
        return x

    def sign(self):
        """
        Returns the sign of this integer, which is -1, 0, or 1
        depending on whether this number is negative, zero, or positive
        respectively.

        OUTPUT: Integer

        EXAMPLES::

            sage: 500.sign()
            1
            sage: 0.sign()
            0
            sage: (-10^43).sign()
            -1
        """
        return smallInteger(mpz_sgn(self.value))

    def __mod__(x, y):
        r"""
         Returns x modulo y.

         EXAMPLES::

             sage: z = 43
             sage: z % 2
             1
             sage: z % 0
             Traceback (most recent call last):
             ...
             ZeroDivisionError: Integer modulo by zero
             sage: -5 % 7
             2
             sage: -5 % -7
             -5
             sage: 5 % -7
             -2

        TESTS::

            sage: signs = [(11,5), (11,-5), (-11,5), (-11,-5)]
            sage: control = [int(a) % int(b) for a, b in signs]
            sage: [a % b for a,b in signs] == control
            True
            sage: [a % int(b) for a,b in signs] == control
            True
            sage: [int(a) % b for a,b in signs] == control
            True

        This example caused trouble in trac #6083::

            sage: a = next_prime(2**31)
            sage: b = Integers(a)(100)
            sage: a % b
            59
         """
        cdef Integer z = PY_NEW(Integer)
        cdef long yy, res

        # first case: Integer % Integer
        if PY_TYPE(x) is PY_TYPE(y):
            if not mpz_sgn((<Integer>y).value):
                raise ZeroDivisionError, "Integer modulo by zero"
            if mpz_size((<Integer>x).value) > 100000:
                sig_on()
                mpz_fdiv_r(z.value, (<Integer>x).value, (<Integer>y).value)
                sig_off()
            else:
                mpz_fdiv_r(z.value, (<Integer>x).value, (<Integer>y).value)
            return z

        # next: Integer % python int
        elif PyInt_CheckExact(y):
            yy = PyInt_AS_LONG(y)
            if yy > 0:
                mpz_fdiv_r_ui(z.value, (<Integer>x).value, yy)
            elif yy == 0:
                raise ZeroDivisionError, "Integer modulo by zero"
            else:
                res = mpz_fdiv_r_ui(z.value, (<Integer>x).value, -yy)
                if res:
                    mpz_sub_ui(z.value, z.value, -yy)
            return z

        # all other cases
        else:
            try:
                # we explicitly try coercing both to ZZ here to
                # avoid infinite loops in some cases (such as
                # Integers and Integers(n)), see trac #6083
                x = integer(x)
                y = integer(y)
                return x % y
            except ValueError:
                return bin_op(x, y, operator.mod)

    def quo_rem(Integer self, other):
        """
        Returns the quotient and the remainder of self divided by other.
        Note that the remainder returned is always either zero or of the
        same sign as other.

        INPUT:

        -  ``other`` - the divisor

        OUTPUT:

        -  ``q`` - the quotient of self/other

        -  ``r`` - the remainder of self/other

        EXAMPLES::

            sage: z = Integer(231)
            sage: z.quo_rem(2)
            (115, 1)
            sage: z.quo_rem(-2)
            (-116, -1)
            sage: z.quo_rem(0)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Integer division by zero

            sage: a = ZZ.random_element(10**50)
            sage: b = ZZ.random_element(10**15)
            sage: q, r = a.quo_rem(b)
            sage: q*b + r == a
            True

            sage: 3.quo_rem(ZZ['x'].0)
            (0, 3)

        TESTS:

        The divisor can be rational as well, although the remainder
        will always be zero (trac #7965)::

            sage: 5.quo_rem(QQ(2))
            (5/2, 0)
            sage: 5.quo_rem(2/3)
            (15/2, 0)

        """
        cdef Integer q = PY_NEW(Integer)
        cdef Integer r = PY_NEW(Integer)
        cdef long d, res

        if PyInt_CheckExact(other):
            d = PyInt_AS_LONG(other)
            if d > 0:
                mpz_fdiv_qr_ui(q.value, r.value, self.value, d)
            elif d == 0:
                raise ZeroDivisionError, "Integer division by zero"
            else:
                res = mpz_fdiv_qr_ui(q.value, r.value, self.value, -d)
                mpz_neg(q.value, q.value)
                if res:
                    mpz_sub_ui(q.value, q.value, 1)
                    mpz_sub_ui(r.value, r.value, -d)

        elif PY_TYPE_CHECK_EXACT(other, Integer):
            if mpz_sgn((<Integer>other).value) == 0:
                raise ZeroDivisionError, "Integer division by zero"
            if mpz_size((<Integer>x).value) > 100000:
                sig_on()
                mpz_fdiv_qr(q.value, r.value, self.value, (<Integer>other).value)
                sig_off()
            else:
                mpz_fdiv_qr(q.value, r.value, self.value, (<Integer>other).value)

        else:
            left, right = canonical_coercion(self, other)
            return left.quo_rem(right)

        return q, r

    def powermod(self, exp, mod):
        """
        Compute self\*\*exp modulo mod.

        EXAMPLES::

            sage: z = 2
            sage: z.powermod(31,31)
            2
            sage: z.powermod(0,31)
            1
            sage: z.powermod(-31,31) == 2^-31 % 31
            True

        As expected, the following is invalid::

            sage: z.powermod(31,0)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: cannot raise to a power modulo 0
        """
        cdef Integer x, _exp, _mod
        _exp = Integer(exp); _mod = Integer(mod)
        if mpz_cmp_si(_mod.value,0) == 0:
            raise ZeroDivisionError, "cannot raise to a power modulo 0"

        x = PY_NEW(Integer)

        sig_on()
        mpz_powm(x.value, self.value, _exp.value, _mod.value)
        sig_off()

        return x

    def rational_reconstruction(self, Integer m):
        """
        Return the rational reconstruction of this integer modulo m, i.e.,
        the unique (if it exists) rational number that reduces to self
        modulo m and whose numerator and denominator is bounded by
        sqrt(m/2).

        EXAMPLES::

            sage: (3/7)%100
            29
            sage: (29).rational_reconstruction(100)
            3/7

        TEST:

        Check that ticket #9345 is fixed::

            sage: ZZ(1).rational_reconstruction(0)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: The modulus cannot be zero
            sage: m = ZZ.random_element(-10^6,10^6)
            sage: m.rational_reconstruction(0)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: The modulus cannot be zero
        """
        import rational
        return rational.pyrex_rational_reconstruction(self, m)

    def powermodm_ui(self, exp, mod):
        r"""
        Computes self\*\*exp modulo mod, where exp is an unsigned long
        integer.

        EXAMPLES::

            sage: z = 32
            sage: z.powermodm_ui(2, 4)
            0
            sage: z.powermodm_ui(2, 14)
            2
            sage: z.powermodm_ui(2^32-2, 14)
            2
            sage: z.powermodm_ui(2^32-1, 14)
            Traceback (most recent call last):                              # 32-bit
            ...                                                             # 32-bit
            OverflowError: exp (=4294967295) must be <= 4294967294          # 32-bit
            8              # 64-bit
            sage: z.powermodm_ui(2^65, 14)
            Traceback (most recent call last):
            ...
            OverflowError: exp (=36893488147419103232) must be <= 4294967294  # 32-bit
            OverflowError: exp (=36893488147419103232) must be <= 18446744073709551614     # 64-bit
        """
        if exp < 0:
            raise ValueError, "exp (=%s) must be nonnegative"%exp
        elif exp > MAX_UNSIGNED_LONG:
            raise OverflowError, "exp (=%s) must be <= %s"%(exp, MAX_UNSIGNED_LONG)
        cdef Integer x, _mod
        _mod = Integer(mod)
        x = PY_NEW(Integer)

        sig_on()
        mpz_powm_ui(x.value, self.value, exp, _mod.value)
        sig_off()

        return x

    def __int__(self):
        """
        Return the Python int (or long) corresponding to this Sage
        integer.

        EXAMPLES::

            sage: n = 920938; n
            920938
            sage: int(n)
            920938
            sage: type(n.__int__())
            <type 'int'>
            sage: n = 99028390823409823904823098490238409823490820938; n
            99028390823409823904823098490238409823490820938
            sage: type(n.__int__())
            <type 'long'>
        """
        # TODO -- this crashes on sage.math, since it is evidently written incorrectly.
        return mpz_get_pyintlong(self.value)
        #return int(mpz_get_pylong(self.value))

    def __long__(self):
        """
        Return long integer corresponding to this Sage integer.

        EXAMPLES::

            sage: n = 9023408290348092849023849820934820938490234290; n
            9023408290348092849023849820934820938490234290
            sage: long(n)
            9023408290348092849023849820934820938490234290L
            sage: n = 920938; n
            920938
            sage: long(n)
            920938L
            sage: n.__long__()
            920938L
        """
        return mpz_get_pylong(self.value)

    def __float__(self):
        """
        Return double precision floating point representation of this
        integer.

        EXAMPLES::

            sage: n = Integer(17); float(n)
            17.0
            sage: n = Integer(902834098234908209348209834092834098); float(n)
            9.028340982349081e+35
            sage: n = Integer(-57); float(n)
            -57.0
            sage: n.__float__()
            -57.0
            sage: type(n.__float__())
            <type 'float'>
        """
        return mpz_get_d(self.value)

    def _rpy_(self):
        """
        Returns int(self) so that rpy can convert self into an object it
        knows how to work with.

        EXAMPLES::

            sage: n = 100
            sage: n._rpy_()
            100
            sage: type(n._rpy_())
            <type 'int'>
        """
        return self.__int__()

    def __hash__(self):
        """
        Return the hash of this integer.

        This agrees with the Python hash of the corresponding Python int or
        long.

        EXAMPLES::

            sage: n = -920384; n.__hash__()
            -920384
            sage: hash(int(n))
            -920384
            sage: n = -920390823904823094890238490238484; n.__hash__()
            -873977844            # 32-bit
            6874330978542788722   # 64-bit
            sage: hash(long(n))
            -873977844            # 32-bit
            6874330978542788722   # 64-bit

        TESTS::

            sage: hash(0)
            0
            sage: hash(-1)
            -2
            sage: n = 2^31 + 2^63 + 2^95 + 2^127 + 2^128*(2^32-2)
            sage: hash(n) == hash(long(n))
            True
            sage: hash(n-1) == hash(long(n-1))
            True
            sage: hash(-n) == hash(long(-n))
            True
            sage: hash(1-n) == hash(long(1-n))
            True
            sage: n = 2^63 + 2^127 + 2^191 + 2^255 + 2^256*(2^64-2)
            sage: hash(n) == hash(long(n))
            True
            sage: hash(n-1) == hash(long(n-1))
            True
            sage: hash(-n) == hash(long(-n))
            True
            sage: hash(1-n) == hash(long(1-n))
            True

        These tests come from Trac #4957:

            sage: n = 2^31 + 2^13
            sage: hash(n)
            -2147475456               # 32-bit
            2147491840                # 64-bit
            sage: hash(n) == hash(int(n))
            True
            sage: n = 2^63 + 2^13
            sage: hash(n)
            -2147475456               # 32-bit
            -9223372036854767616      # 64-bit
            sage: hash(n) == hash(int(n))
            True
        """
        return mpz_pythonhash(self.value)

    cdef hash_c(self):
        """
        A C version of the __hash__ function.
        """
        return mpz_pythonhash(self.value)

    def trial_division(self, long bound=LONG_MAX, long start=2):
        """
        Return smallest prime divisor of self up to bound, beginning
        checking at start, or abs(self) if no such divisor is found.

        INPUT:

            - ``bound`` -- a positive integer that fits in a C signed long
            - ``start`` -- a positive integer that fits in a C signed long

        OUTPUT:

            - a positive integer

        EXAMPLES::

            sage: n = next_prime(10^6)*next_prime(10^7); n.trial_division()
            1000003
            sage: (-n).trial_division()
            1000003
            sage: n.trial_division(bound=100)
            10000049000057
            sage: n.trial_division(bound=-10)
            Traceback (most recent call last):
            ...
            ValueError: bound must be positive
            sage: n.trial_division(bound=0)
            Traceback (most recent call last):
            ...
            ValueError: bound must be positive
            sage: ZZ(0).trial_division()
            Traceback (most recent call last):
            ...
            ValueError: self must be nonzero

            sage: n = next_prime(10^5) * next_prime(10^40); n.trial_division()
            100003
            sage: n.trial_division(bound=10^4)
            1000030000000000000000000000000000000012100363
            sage: (-n).trial_division(bound=10^4)
            1000030000000000000000000000000000000012100363
            sage: (-n).trial_division()
            100003
            sage: n = 2 * next_prime(10^40); n.trial_division()
            2
            sage: n = 3 * next_prime(10^40); n.trial_division()
            3
            sage: n = 5 * next_prime(10^40); n.trial_division()
            5
            sage: n = 2 * next_prime(10^4); n.trial_division()
            2
            sage: n = 3 * next_prime(10^4); n.trial_division()
            3
            sage: n = 5 * next_prime(10^4); n.trial_division()
            5

        You can specify a starting point::

            sage: n = 3*5*101*103
            sage: n.trial_division(start=50)
            101
        """
        if bound <= 0:
            raise ValueError, "bound must be positive"
        if mpz_sgn(self.value) == 0:
            raise ValueError, "self must be nonzero"
        cdef unsigned long n, m=7, i=1, limit, dif[8]
        if start > 7:
            # We need to find i.
            m = start % 30
            if 0 <= m <= 1:
                i = 0; m = start + (1-m)
            elif 1 < m <= 7:
                i = 1; m = start + (7-m)
            elif 7 < m <= 11:
                i = 2; m = start + (11-m)
            elif 11 < m <= 13:
                i = 3; m = start + (13-m)
            elif 13 < m <= 17:
                i = 4; m = start + (17-m)
            elif 17 < m <= 19:
                i = 5; m = start + (19-m)
            elif 19 < m <= 23:
                i = 6; m = start + (23-m)
            elif 23 < m <= 29:
                i = 7; m = start + (29-m)
        dif[0]=6;dif[1]=4;dif[2]=2;dif[3]=4;dif[4]=2;dif[5]=4;dif[6]=6;dif[7]=2
        cdef Integer x = PY_NEW(Integer)
        if mpz_fits_ulong_p(self.value):
            n = mpz_get_ui(self.value)   # ignores the sign automatically
            if n == 1: return one
            if start <= 2 and n%2==0:
                mpz_set_ui(x.value,2); return x
            if start <= 3 and n%3==0:
                mpz_set_ui(x.value,3); return x
            if start <= 5 and n%5==0:
                mpz_set_ui(x.value,5); return x
            limit = <unsigned long> sqrt_double(<double> n)
            if bound < limit: limit = bound
            # Algorithm: only trial divide by numbers that
            # are congruent to 1,7,11,13,17,19,23,29 mod 30=2*3*5.
            while m <= limit:
                if n%m == 0:
                    mpz_set_ui(x.value, m); return x
                m += dif[i%8]
                i += 1
            mpz_abs(x.value, self.value)
            return x
        else:
            # self is big -- it doesn't fit in unsigned long.
            if start <= 2 and mpz_even_p(self.value):
                mpz_set_ui(x.value,2); return x
            if start <= 3 and mpz_divisible_ui_p(self.value,3):
                mpz_set_ui(x.value,3); return x
            if start <= 5 and mpz_divisible_ui_p(self.value,5):
                mpz_set_ui(x.value,5); return x

            # x.value = floor(sqrt(self.value))
            sig_on()
            mpz_abs(x.value, self.value)
            mpz_sqrt(x.value, x.value)
            if mpz_cmp_si(x.value, bound) < 0:
                limit = mpz_get_ui(x.value)
            else:
                limit = bound
            while m <= limit:
                if  mpz_divisible_ui_p(self.value, m):
                    mpz_set_ui(x.value, m)
                    sig_off()
                    return x
                m += dif[i%8]
                i += 1
            mpz_abs(x.value, self.value)
            sig_off()
            return x

    def factor(self, algorithm='pari', proof=None, limit=None, int_=False,
                     verbose=0):
        """
        Return the prime factorization of this integer as a
        formal Factorization object.

        INPUT:

        -  ``algorithm`` - string

           - ``'pari'`` - (default) use the PARI library

           - ``'kash'`` - use the KASH computer algebra system (requires
             the optional kash package)

        -  ``proof`` - bool (default: True) whether or not to
           prove primality of each factor (only applicable for PARI).

        -  ``limit`` - int or None (default: None) if limit is
           given it must fit in a signed int, and the factorization is done
           using trial division and primes up to limit.

        OUTPUT:

        -  a Factorization object containing the prime factors and
           their multiplicities

        EXAMPLES::

            sage: n = 2^100 - 1; n.factor()
            3 * 5^3 * 11 * 31 * 41 * 101 * 251 * 601 * 1801 * 4051 * 8101 * 268501

        This factorization can be converted into a list of pairs `(p,
        e)`, where `p` is prime and `e` is a positive integer.  Each
        pair can also be accessed directly by its index (ordered by
        increasing size of the prime)::

            sage: f = 60.factor()
            sage: list(f)
            [(2, 2), (3, 1), (5, 1)]
            sage: f[2]
            (5, 1)

        Similarly, the factorization can be converted to a dictionary
        so the exponent can be extracted for each prime::

            sage: f = (3^6).factor()
            sage: dict(f)
            {3: 6}
            sage: dict(f)[3]
            6

        We use proof=False, which doesn't prove correctness of the primes
        that appear in the factorization::

            sage: n = 920384092842390423848290348203948092384082349082
            sage: n.factor(proof=False)
            2 * 11 * 1531 * 4402903 * 10023679 * 619162955472170540533894518173
            sage: n.factor(proof=True)
            2 * 11 * 1531 * 4402903 * 10023679 * 619162955472170540533894518173

        We factor using trial division only::

            sage: n.factor(limit=1000)
            2 * 11 * 41835640583745019265831379463815822381094652231

        TESTS::

            sage: n.factor(algorithm='foobar')
            Traceback (most recent call last):
            ...
            ValueError: Algorithm is not known
        """
        from sage.structure.factorization import Factorization
        from sage.structure.factorization_integer import IntegerFactorization

        cdef Integer n, p, unit
        cdef int i
        cdef n_factor_t f

        if mpz_sgn(self.value) == 0:
            raise ArithmeticError, "Prime factorization of 0 not defined."

        if mpz_sgn(self.value) > 0:
            n    = self
            unit = one
        else:
            n    = PY_NEW(Integer)
            unit = PY_NEW(Integer)
            mpz_neg(n.value, self.value)
            mpz_set_si(unit.value, -1)

        if mpz_cmpabs_ui(n.value, 1) == 0:
            return IntegerFactorization([], unit=unit, unsafe=True,
                                            sort=False, simplify=False)

        if limit is not None:
            from sage.rings.factorint import factor_trial_division
            return factor_trial_division(self, limit)

        if mpz_fits_slong_p(n.value):
            if proof is None:
                from sage.structure.proof.proof import get_flag
                proof = get_flag(proof, "arithmetic")
            n_factor_init(&f)
            n_factor(&f, mpz_get_ui(n.value), proof)
            F = [(Integer(f.p[i]), int(f.exp[i])) for i from 0 <= i < f.num]
            F.sort()
            return IntegerFactorization(F, unit=unit, unsafe=True,
                                           sort=False, simplify=False)

        if mpz_sizeinbase(n.value, 2) < 40:
            from sage.rings.factorint import factor_trial_division
            return factor_trial_division(self)

        if algorithm == 'pari':
            from sage.rings.factorint import factor_using_pari
            F = factor_using_pari(n, int_=int_, debug_level=verbose, proof=proof)
            F.sort()
            return IntegerFactorization(F, unit=unit, unsafe=True,
                                           sort=False, simplify=False)
        elif algorithm in ['kash', 'magma']:
            if algorithm == 'kash':
                from sage.interfaces.all import kash as I
            else:
                from sage.interfaces.all import magma as I
            str_res = I.eval('Factorization(%s)'%n)
            # The result looks like "[ <n1, p1>, <p2, e2>, ... ]
            str_res = str_res.replace(']', '').replace('[', '').replace('>', '').replace('<', '').split(',')
            res = [int(s.strip()) for s in str_res]
            exp_type = int if int_ else Integer
            F = [(Integer(p), exp_type(e)) for p,e in zip(res[0::2], res[1::2])]
            return Factorization(F, unit)
        else:
            raise ValueError, "Algorithm is not known"

    def support(self):
        """
        Return a sorted list of the primes dividing this integer.

        OUTPUT: The sorted list of primes appearing in the factorization of
        this rational with positive exponent.

        EXAMPLES::

            sage: factorial(10).support()
            [2, 3, 5, 7]
            sage: (-999).support()
            [3, 37]

        Trying to find the support of 0 gives an arithmetic error::

            sage: 0.support()
            Traceback (most recent call last):
            ...
            ArithmeticError: Support of 0 not defined.
        """
        if self.is_zero():
            raise ArithmeticError, "Support of 0 not defined."
        return sage.rings.arith.prime_factors(self)

    def coprime_integers(self, m):
        """
        Return the positive integers `< m` that are coprime to
        self.

        EXAMPLES::

            sage: n = 8
            sage: n.coprime_integers(8)
            [1, 3, 5, 7]
            sage: n.coprime_integers(11)
            [1, 3, 5, 7, 9]
            sage: n = 5; n.coprime_integers(10)
            [1, 2, 3, 4, 6, 7, 8, 9]
            sage: n.coprime_integers(5)
            [1, 2, 3, 4]
            sage: n = 99; n.coprime_integers(99)
            [1, 2, 4, 5, 7, 8, 10, 13, 14, 16, 17, 19, 20, 23, 25, 26, 28, 29, 31, 32, 34, 35, 37, 38, 40, 41, 43, 46, 47, 49, 50, 52, 53, 56, 58, 59, 61, 62, 64, 65, 67, 68, 70, 71, 73, 74, 76, 79, 80, 82, 83, 85, 86, 89, 91, 92, 94, 95, 97, 98]

        AUTHORS:

        - Naqi Jaffery (2006-01-24): examples

        ALGORITHM: Naive - compute lots of GCD's. If this isn't good enough
        for you, please code something better and submit a patch.
        """
        # TODO -- make VASTLY faster
        v = []
        for n in range(1,m):
            if self.gcd(n) == 1:
                v.append(Integer(n))
        return v

    def divides(self, n):
        """
        Return True if self divides n.

        EXAMPLES::

            sage: Z = IntegerRing()
            sage: Z(5).divides(Z(10))
            True
            sage: Z(0).divides(Z(5))
            False
            sage: Z(10).divides(Z(5))
            False
        """
        cdef bint t
        cdef Integer _n
        _n = Integer(n)
        if mpz_sgn(self.value) == 0:
            return mpz_sgn(_n.value) == 0
        sig_on()
        t = mpz_divisible_p(_n.value, self.value)
        sig_off()
        return t

    cpdef RingElement _valuation(Integer self, Integer p):
        r"""
        Return the p-adic valuation of self.

        We do not require that p be prime, but it must be at least 2. For
        more documentation see ``valuation``

        AUTHORS:

        - David Roe (3/31/07)
        """
        if mpz_sgn(self.value) == 0:
            return sage.rings.infinity.infinity
        if mpz_cmp_ui(p.value, 2) < 0:
            raise ValueError, "You can only compute the valuation with respect to a integer larger than 1."

        cdef Integer v = PY_NEW(Integer)
        cdef mpz_t u
        mpz_init(u)
        sig_on()
        mpz_set_ui(v.value, mpz_remove(u, self.value, p.value))
        sig_off()
        mpz_clear(u)
        return v

    cdef object _val_unit(Integer self, Integer p):
        r"""
        Returns a pair: the p-adic valuation of self, and the p-adic unit
        of self.

        We do not require the p be prime, but it must be at least 2. For
        more documentation see ``val_unit``

        AUTHORS:

        - David Roe (2007-03-31)
        """
        cdef Integer v, u
        if mpz_cmp_ui(p.value, 2) < 0:
            raise ValueError, "You can only compute the valuation with respect to a integer larger than 1."
        if self == 0:
            u = one
            return (sage.rings.infinity.infinity, u)
        v = PY_NEW(Integer)
        u = PY_NEW(Integer)
        sig_on()
        mpz_set_ui(v.value, mpz_remove(u.value, self.value, p.value))
        sig_off()
        return (v, u)

    def valuation(self, p):
        """
        Return the p-adic valuation of self.

        INPUT:

        -  ``p`` - an integer at least 2.

        EXAMPLE::

            sage: n = 60
            sage: n.valuation(2)
            2
            sage: n.valuation(3)
            1
            sage: n.valuation(7)
            0
            sage: n.valuation(1)
            Traceback (most recent call last):
            ...
            ValueError: You can only compute the valuation with respect to a integer larger than 1.

        We do not require that p is a prime::

            sage: (2^11).valuation(4)
            5
        """
        return self._valuation(Integer(p))

    # Alias for valuation
    ord = valuation

    def val_unit(self, p):
        r"""
        Returns a pair: the p-adic valuation of self, and the p-adic unit
        of self.

        INPUT:

        -  ``p`` - an integer at least 2.

        OUTPUT:

        -  ``v_p(self)`` - the p-adic valuation of ``self``

        -  ``u_p(self)`` - ``self`` / `p^{v_p(\mathrm{self})}`

        EXAMPLE::

            sage: n = 60
            sage: n.val_unit(2)
            (2, 15)
            sage: n.val_unit(3)
            (1, 20)
            sage: n.val_unit(7)
            (0, 60)
            sage: (2^11).val_unit(4)
            (5, 2)
            sage: 0.val_unit(2)
            (+Infinity, 1)
        """
        return self._val_unit(Integer(p))

    def odd_part(self):
        r"""
        The odd part of the integer `n`. This is `n / 2^v`,
        where `v = \mathrm{valuation}(n,2)`.

        IMPLEMENTATION:

        Currently returns 0 when self is 0.  This behaviour is fairly arbitrary,
        and in Sage 4.6 this special case was not handled at all, eventually
        propagating a TypeError.  The caller should not rely on the behaviour
        in case self is 0.

        EXAMPLES::

            sage: odd_part(5)
            5
            sage: odd_part(4)
            1
            sage: odd_part(factorial(31))
            122529844256906551386796875
        """
        cdef Integer odd
        cdef unsigned long bits

        if mpz_cmpabs_ui(self.value, 1) <= 0:
            return self

        odd  = PY_NEW(Integer)
        bits = mpz_scan1(self.value, 0)
        mpz_tdiv_q_2exp(odd.value, self.value, bits)
        return odd

    cdef Integer _divide_knowing_divisible_by(Integer self, Integer right):
        r"""
        Returns the integer self / right when self is divisible by right.

        If self is not divisible by right, the return value is undefined,
        and may not even be close to self/right. For more documentation see
        ``divide_knowing_divisible_by``

        AUTHORS:

        - David Roe (2007-03-31)
        """
        if mpz_cmp_ui(right.value, 0) == 0:
            raise ZeroDivisionError
        cdef Integer x
        x = PY_NEW(Integer)
        if mpz_size(self.value) + mpz_size((<Integer>right).value) > 100000:
            # Only use the signal handler (to enable ctrl-c out) when the
            # quotient might take a while to compute
            sig_on()
            mpz_divexact(x.value, self.value, right.value)
            sig_off()
        else:
            mpz_divexact(x.value, self.value, right.value)
        return x

    def divide_knowing_divisible_by(self, right):
        r"""
        Returns the integer self / right when self is divisible by right.

        If self is not divisible by right, the return value is undefined,
        and may not even be close to self/right for multi-word integers.

        EXAMPLES::

            sage: a = 8; b = 4
            sage: a.divide_knowing_divisible_by(b)
            2
            sage: (100000).divide_knowing_divisible_by(25)
            4000
            sage: (100000).divide_knowing_divisible_by(26) # close (random)
            3846

        However, often it's way off.

        ::

            sage: a = 2^70; a
            1180591620717411303424
            sage: a // 11  # floor divide
            107326510974310118493
            sage: a.divide_knowing_divisible_by(11) # way off and possibly random
            43215361478743422388970455040
        """
        return self._divide_knowing_divisible_by(right)

    def _lcm(self, Integer n):
        """
        Returns the least common multiple of self and `n`.

        EXAMPLES::

            sage: n = 60
            sage: n._lcm(150)
            300
        """
        cdef Integer z = PY_NEW(Integer)
        sig_on()
        mpz_lcm(z.value, self.value, n.value)
        sig_off()
        return z

    def denominator(self):
        """
        Return the denominator of this integer, which of course is
        always 1.

        EXAMPLES::

            sage: x = 5
            sage: x.denominator()
            1
            sage: x = 0
            sage: x.denominator()
            1
        """
        return one

    def numerator(self):
        """
        Return the numerator of this integer.

        EXAMPLES::

            sage: x = 5
            sage: x.numerator()
            5

        ::

            sage: x = 0
            sage: x.numerator()
            0
        """
        return self

    def factorial(self):
        r"""
        Return the factorial `n! = 1 \cdot 2 \cdot 3 \cdots n`.

        If the input does not fit in an ``unsigned long int`` a symbolic
        expression is returned.

        EXAMPLES::

            sage: for n in srange(7):
            ...    print n, n.factorial()
            0 1
            1 1
            2 2
            3 6
            4 24
            5 120
            6 720
            sage: 234234209384023842034.factorial()
            factorial(234234209384023842034)
        """
        if mpz_sgn(self.value) < 0:
            raise ValueError, "factorial -- self = (%s) must be nonnegative"%self

        if not mpz_fits_uint_p(self.value):
            from sage.functions.all import factorial
            return factorial(self, hold=True)

        cdef Integer z = PY_NEW(Integer)

        sig_on()
        mpz_fac_ui(z.value, mpz_get_ui(self.value))
        sig_off()

        return z

    @cython.cdivision(True)
    def multifactorial(self, int k):
        r"""
        Computes the k-th factorial `n!^{(k)}` of self. For k=1
        this is the standard factorial, and for k greater than one it is
        the product of every k-th terms down from self to k. The recursive
        definition is used to extend this function to the negative
        integers.

        EXAMPLES::

            sage: 5.multifactorial(1)
            120
            sage: 5.multifactorial(2)
            15
            sage: 23.multifactorial(2)
            316234143225
            sage: prod([1..23, step=2])
            316234143225
            sage: (-29).multifactorial(7)
            1/2640
        """
        if k <= 0:
            raise ValueError, "multifactorial only defined for positive values of k"

        if not mpz_fits_sint_p(self.value):
            raise ValueError, "multifactorial not implemented for n >= 2^32.\nThis is probably OK, since the answer would have billions of digits."

        cdef int n = mpz_get_si(self.value)

        # base case
        if 0 < n < k:
            return one

        # easy to calculate
        elif n % k == 0:
            factorial = Integer(n/k).factorial()
            if k == 2:
                return factorial << (n/k)
            else:
                return factorial * Integer(k)**(n/k)

        # negative base case
        elif -k < n < 0:
            return one / (self+k)

        # reflection case
        elif n < -k:
            if (n/k) % 2:
                sign = -one
            else:
                sign = one
            return sign / Integer(-k-n).multifactorial(k)

        # compute the actual product, optimizing the number of large
        # multiplications
        cdef int i,j

        # we need (at most) log_2(#factors) concurrent sub-products
        cdef int prod_count = <int>ceil_c(log_c(n/k+1)/log_c(2))
        cdef mpz_t* sub_prods = <mpz_t*>sage_malloc(prod_count * sizeof(mpz_t))
        if sub_prods == NULL:
            raise MemoryError
        for i from 0 <= i < prod_count:
            mpz_init(sub_prods[i])

        sig_on()

        cdef residue = n % k
        cdef int tip = 0
        for i from 1 <= i <= n//k:
            mpz_set_ui(sub_prods[tip], k*i + residue)
            # for the i-th terms we use the bits of i to calculate how many
            # times we need to multiply "up" the stack of sub-products
            for j from 0 <= j < 32:
                if i & (1 << j):
                    break
                tip -= 1
                mpz_mul(sub_prods[tip], sub_prods[tip], sub_prods[tip+1])
            tip += 1
        cdef int last = tip-1
        for tip from last > tip >= 0:
            mpz_mul(sub_prods[tip], sub_prods[tip], sub_prods[tip+1])

        sig_off()

        cdef Integer z = PY_NEW(Integer)
        mpz_swap(z.value, sub_prods[0])

        for i from 0 <= i < prod_count:
            mpz_clear(sub_prods[i])
        sage_free(sub_prods)

        return z

    def gamma(self):
        r"""
        The gamma function on integers is the factorial function (shifted by
        one) on positive integers, and `\pm \infty` on non-positive integers.

        EXAMPLES::

            sage: gamma(5)
            24
            sage: gamma(0)
            Infinity
            sage: gamma(-1)
            Infinity
            sage: gamma(-2^150)
            Infinity
        """
        if mpz_sgn(self.value) > 0:
            return (self-one).factorial()
        else:
            from sage.rings.infinity import unsigned_infinity
            return unsigned_infinity

    def floor(self):
        """
        Return the floor of self, which is just self since self is an
        integer.

        EXAMPLES::

            sage: n = 6
            sage: n.floor()
            6
        """
        return self

    def ceil(self):
        """
        Return the ceiling of self, which is self since self is an
        integer.

        EXAMPLES::

            sage: n = 6
            sage: n.ceil()
            6
        """
        return self

    def real(self):
        """
        Returns the real part of self, which is self.

        EXAMPLES::

            sage: Integer(-4).real()
            -4
        """
        return self

    def imag(self):
        """
        Returns the imaginary part of self, which is zero.

        EXAMPLES::

            sage: Integer(9).imag()
            0
        """
        return zero

    def is_one(self):
        r"""
        Returns ``True`` if the integer is `1`, otherwise ``False``.

        EXAMPLES::

            sage: Integer(1).is_one()
            True
            sage: Integer(0).is_one()
            False
        """
        return mpz_cmp_si(self.value, 1) == 0

    def __nonzero__(self):
        r"""
        Returns ``True`` if the integer is not `0`, otherwise ``False``.

        EXAMPLES::

            sage: Integer(1).is_zero()
            False
            sage: Integer(0).is_zero()
            True
        """
        return mpz_sgn(self.value) != 0

    def is_integral(self):
        """
        Return ``True`` since integers are integral, i.e.,
        satisfy a monic polynomial with integer coefficients.

        EXAMPLES::

            sage: Integer(3).is_integral()
            True
        """
        return True

    def is_unit(self):
        r"""
        Returns ``true`` if this integer is a unit, i.e., 1 or `-1`.

        EXAMPLES::

            sage: for n in srange(-2,3):
            ...    print n, n.is_unit()
            -2 False
            -1 True
            0 False
            1 True
            2 False
        """
        return mpz_cmpabs_ui(self.value, 1) == 0

    def is_square(self):
        r"""
        Returns ``True`` if self is a perfect square.

        EXAMPLES::

            sage: Integer(4).is_square()
            True
            sage: Integer(41).is_square()
            False
        """
        return mpz_perfect_square_p(self.value)

    is_power = deprecated_function_alias(12116, is_perfect_power)

    def perfect_power(self):
        r"""
        Returns ``(a, b)``, where this integer is `a^b` and `b` is maximal.

        If called on `-1`, `0` or `1`, `b` will be `1`, since there is no
        maximal value of `b`.

        .. seealso::

            - :meth:`is_perfect_power`: testing whether an integer is a perfect
              power is usually faster than finding `a` and `b`.
            - :meth:`is_prime_power`: checks whether the base is prime.
            - :meth:`is_power_of`: if you know the base already, this method is
              the fastest option.

        EXAMPLES::

            sage: 144.perfect_power()
            (12, 2)
            sage: 1.perfect_power()
            (1, 1)
            sage: 0.perfect_power()
            (0, 1)
            sage: (-1).perfect_power()
            (-1, 1)
            sage: (-8).perfect_power()
            (-2, 3)
            sage: (-4).perfect_power()
            (-4, 1)
            sage: (101^29).perfect_power()
            (101, 29)
            sage: (-243).perfect_power()
            (-3, 5)
            sage: (-64).perfect_power()
            (-4, 3)
        """
        parians = self._pari_().ispower()
        return Integer(parians[1]), Integer(parians[0])

    def global_height(self, prec=None):
        r"""
        Returns the absolute logarithmic height of this rational integer.

        INPUT:

        - ``prec`` (int) -- desired floating point precision (default:
          default RealField precision).

        OUTPUT:

        (real) The absolute logarithmic height of this rational integer.

        ALGORITHM:

        The height of the integer `n` is `\log |n|`.

        EXAMPLES::

            sage: ZZ(5).global_height()
            1.60943791243410
            sage: ZZ(-2).global_height(prec=100)
            0.69314718055994530941723212146
            sage: exp(_)
            2.0000000000000000000000000000
        """
        from sage.rings.real_mpfr import RealField
        if prec is None:
            R = RealField()
        else:
            R = RealField(prec)
        if self.is_zero():
            return R.zero_element()
        return R(self).abs().log()

    cdef bint _is_power_of(Integer self, Integer n):
        r"""
        Returns a non-zero int if there is an integer b with
        `\mathtt{self} = n^b`.

        For more documentation see ``is_power_of``.

        AUTHORS:

        - David Roe (2007-03-31)
        """
        cdef int a
        cdef unsigned long b, c
        cdef mpz_t u, sabs, nabs
        a = mpz_cmp_ui(n.value, 2)
        if a <= 0: # n <= 2
            if a == 0: # n == 2
                if mpz_popcount(self.value) == 1: #number of bits set in self == 1
                    return 1
                else:
                    return 0
            a = mpz_cmp_si(n.value, -2)
            if a >= 0: # -2 <= n < 2:
                a = mpz_get_si(n.value)
                if a == 1: # n == 1
                    if mpz_cmp_ui(self.value, 1) == 0: # Only 1 is a power of 1
                        return 1
                    else:
                        return 0
                elif a == 0: # n == 0
                    if mpz_cmp_ui(self.value, 0) == 0 or mpz_cmp_ui(self.value, 1) == 0: # 0^0 = 1, 0^x = 0
                        return 1
                    else:
                        return 0
                elif a == -1: # n == -1
                    if mpz_cmp_ui(self.value, 1) == 0 or mpz_cmp_si(self.value, -1) == 0: # 1 and -1 are powers of -1
                        return 1
                    else:
                        return 0
                elif a == -2: # n == -2
                    mpz_init(sabs)
                    mpz_abs(sabs, self.value)
                    if mpz_popcount(sabs) == 1: # number of bits set in |self| == 1
                        b = mpz_scan1(sabs, 0) % 2 # b == 1 if |self| is an odd power of 2, 0 if |self| is an even power
                        mpz_clear(sabs)
                        if (b == 1 and mpz_cmp_ui(self.value, 0) < 0) or (b == 0 and mpz_cmp_ui(self.value, 0) > 0):
                            # An odd power of -2 is negative, an even power must be positive.
                            return 1
                        else: # number of bits set in |self| is not 1, so self cannot be a power of -2
                            return 0
                    else: # |self| is not a power of 2, so self cannot be a power of -2
                        return 0
            else: # n < -2
                mpz_init(nabs)
                mpz_neg(nabs, n.value)
                if mpz_popcount(nabs) == 1: # |n| = 2^k for k >= 2.  We special case this for speed
                    mpz_init(sabs)
                    mpz_abs(sabs, self.value)
                    if mpz_popcount(sabs) == 1: # |self| = 2^L for some L >= 0.
                        b = mpz_scan1(sabs, 0) # the bit that self is set at
                        c = mpz_scan1(nabs, 0) # the bit that n is set at
                        # Having obtained b and c, we're done with nabs and sabs (on this branch anyway)
                        mpz_clear(nabs)
                        mpz_clear(sabs)
                        if b % c == 0: # Now we know that |self| is a power of |n|
                            b = (b // c) % 2 # Whether b // c is even or odd determines whether (-2^c)^(b // c) is positive or negative
                            a = mpz_cmp_ui(self.value, 0)
                            if b == 0 and a > 0 or b == 1 and a < 0:
                                # These two cases are that b // c is even and self positive, or b // c is odd and self negative
                                return 1
                            else: # The sign of self is wrong
                                return 0
                        else: # Since |self| is not a power of |n|, self cannot be a power of n
                            return 0
                    else: # self is not a power of 2, and thus cannot be a power of n, which is a power of 2.
                        mpz_clear(nabs)
                        mpz_clear(sabs)
                        return 0
                else: # |n| is not a power of 2, so we use mpz_remove
                    mpz_init(u)
                    sig_on()
                    b = mpz_remove(u, self.value, nabs)
                    sig_off()
                    # Having obtained b and u, we're done with nabs
                    mpz_clear(nabs)
                    if mpz_cmp_ui(u, 1) == 0: # self is a power of |n|
                        mpz_clear(u)
                        if b % 2 == 0: # an even power of |n|, and since self > 0, this means that self is a power of n
                            return 1
                        else:
                            return 0
                    elif mpz_cmp_si(u, -1) == 0: # -self is a power of |n|
                        mpz_clear(u)
                        if b % 2 == 1: # an odd power of |n|, and thus self is a power of n
                            return 1
                        else:
                            return 0
                    else: # |self| is not a power of |n|, so self cannot be a power of n
                        mpz_clear(u)
                        return 0
        elif mpz_popcount(n.value) == 1: # n > 2 and in fact n = 2^k for k >= 2
            if mpz_popcount(self.value) == 1: # since n is a power of 2, so must self be.
                if mpz_scan1(self.value, 0) % mpz_scan1(n.value, 0) == 0: # log_2(self) is divisible by log_2(n)
                    return 1
                else:
                    return 0
            else: # self is not a power of 2, and thus not a power of n
                return 0
        else: # n > 2, but not a power of 2, so we use mpz_remove
            mpz_init(u)
            sig_on()
            mpz_remove(u, self.value, n.value)
            sig_off()
            a = mpz_cmp_ui(u, 1)
            mpz_clear(u)
            if a == 0:
                return 1
            else:
                return 0

    def is_power_of(Integer self, n):
        r"""
        Returns ``True`` if there is an integer b with
        `\mathtt{self} = n^b`.

        .. seealso::

            - :meth:`perfect_power`: Finds the minimal base for which this
              integer is a perfect power.
            - :meth:`is_perfect_power`: If you don't know the base but just
              want to know if this integer is a perfect power, use this
              function.
            - :meth:`is_prime_power`: Checks whether the base is prime.

        EXAMPLES::

            sage: Integer(64).is_power_of(4)
            True
            sage: Integer(64).is_power_of(16)
            False

        TESTS::

            sage: Integer(-64).is_power_of(-4)
            True
            sage: Integer(-32).is_power_of(-2)
            True
            sage: Integer(1).is_power_of(1)
            True
            sage: Integer(-1).is_power_of(-1)
            True
            sage: Integer(0).is_power_of(1)
            False
            sage: Integer(0).is_power_of(0)
            True
            sage: Integer(1).is_power_of(0)
            True
            sage: Integer(1).is_power_of(8)
            True
            sage: Integer(-8).is_power_of(2)
            False
            sage: Integer(-81).is_power_of(-3)
            False

        .. note::

           For large integers self, is_power_of() is faster than
           is_perfect_power(). The following examples gives some indication of
           how much faster.

        ::

            sage: b = lcm(range(1,10000))
            sage: b.exact_log(2)
            14446
            sage: t=cputime()
            sage: for a in range(2, 1000): k = b.is_perfect_power()
            sage: cputime(t)      # random
            0.53203299999999976
            sage: t=cputime()
            sage: for a in range(2, 1000): k = b.is_power_of(2)
            sage: cputime(t)      # random
            0.0
            sage: t=cputime()
            sage: for a in range(2, 1000): k = b.is_power_of(3)
            sage: cputime(t)      # random
            0.032002000000000308

        ::

            sage: b = lcm(range(1, 1000))
            sage: b.exact_log(2)
            1437
            sage: t=cputime()
            sage: for a in range(2, 10000): k = b.is_perfect_power() # note that we change the range from the example above
            sage: cputime(t)      # random
            0.17201100000000036
            sage: t=cputime(); TWO=int(2)
            sage: for a in range(2, 10000): k = b.is_power_of(TWO)
            sage: cputime(t)      # random
            0.0040000000000000036
            sage: t=cputime()
            sage: for a in range(2, 10000): k = b.is_power_of(3)
            sage: cputime(t)      # random
            0.040003000000000011
            sage: t=cputime()
            sage: for a in range(2, 10000): k = b.is_power_of(a)
            sage: cputime(t)      # random
            0.02800199999999986
        """
        if not PY_TYPE_CHECK(n, Integer):
            n = Integer(n)
        return self._is_power_of(n)

    def is_prime_power(self, flag=0):
        r"""
        Returns True if this integer is a prime power, and False otherwise.

        INPUT:

        - ``flag`` (for primality testing) - int. Values are:

          - ``0`` (default): use a combination of algorithms.
          - ``1``: certify primality using the Pocklington-Lehmer test.
          - ``2``: certify primality using the APRCL test.

        .. seealso::

            - :meth:`perfect_power`: Finds the minimal base for which integer
              is a perfect power.
            - :meth:`is_perfect_power`: Doesn't test whether the base is prime.
            - :meth:`is_power_of`: If you know the base already this method is
              the fastest option.

        EXAMPLES::

            sage: (-10).is_prime_power()
            False
            sage: (10).is_prime_power()
            False
            sage: (64).is_prime_power()
            True
            sage: (3^10000).is_prime_power()
            True
            sage: (10000).is_prime_power(flag=1)
            False

        We check that :trac:`4777` is fixed::

            sage: n = 150607571^14
            sage: n.is_prime_power()
            True
        """
        if self.is_zero():
            return False
        elif self.is_one():
            return True
        elif mpz_sgn(self.value) < 0:
            return False
        if self.is_prime():
            return True
        if not self.is_perfect_power():
            return False
        k, g = self._pari_().ispower()
        if k == 1:
            raise RuntimeError, "Inconsistent results between GMP and PARI"
        return g.isprime(flag=flag)

    def is_prime(self, proof=None):
        r"""
        Returns ``True`` if ``self`` is prime.

        INPUT:

        - ``proof`` -- If False, use a strong pseudo-primality test.
          If True, use a provable primality test.  If unset, use the
          default arithmetic proof flag.

        .. note::

           Integer primes are by definition *positive*! This is
           different than Magma, but the same as in PARI. See also the
           :meth:`is_irreducible()` method.

        EXAMPLES::

            sage: z = 2^31 - 1
            sage: z.is_prime()
            True
            sage: z = 2^31
            sage: z.is_prime()
            False
            sage: z = 7
            sage: z.is_prime()
            True
            sage: z = -7
            sage: z.is_prime()
            False
            sage: z.is_irreducible()
            True

        ::

            sage: z = 10^80 + 129
            sage: z.is_prime(proof=False)
            True
            sage: z.is_prime(proof=True)
            True

        IMPLEMENTATION: Calls the PARI ``isprime`` function.
        """
        from sage.structure.proof.proof import get_flag
        proof = get_flag(proof, "arithmetic")
        if proof:
            return bool(self._pari_().isprime())
        else:
            return bool(self._pari_().ispseudoprime())

    def is_irreducible(self):
        r"""
        Returns ``True`` if self is irreducible, i.e. +/-
        prime

        EXAMPLES::

            sage: z = 2^31 - 1
            sage: z.is_irreducible()
            True
            sage: z = 2^31
            sage: z.is_irreducible()
            False
            sage: z = 7
            sage: z.is_irreducible()
            True
            sage: z = -7
            sage: z.is_irreducible()
            True
        """
        n = self if self >= 0 else -self
        return bool(n._pari_().isprime())

    def is_pseudoprime(self):
        r"""
        Returns ``True`` if self is a pseudoprime

        EXAMPLES::

            sage: z = 2^31 - 1
            sage: z.is_pseudoprime()
            True
            sage: z = 2^31
            sage: z.is_pseudoprime()
            False
        """
        return bool(self._pari_().ispseudoprime())

    def is_perfect_power(self):
        r"""
        Returns ``True`` if ``self`` is a perfect power, ie if there exist integers
        `a` and `b`, `b > 1` with ``self`` `= a^b`.

        .. seealso::

            - :meth:`perfect_power`: Finds the minimal base for which this
              integer is a perfect power.
            - :meth:`is_power_of`: If you know the base already this method is
              the fastest option.
            - :meth:`is_prime_power`: Checks whether the base is prime.

        EXAMPLES::

            sage: Integer(-27).is_perfect_power()
            True
            sage: Integer(12).is_perfect_power()
            False

            sage: z = 8
            sage: z.is_perfect_power()
            True
            sage: 144.is_perfect_power()
            True
            sage: 10.is_perfect_power()
            False
            sage: (-8).is_perfect_power()
            True
            sage: (-4).is_perfect_power()
            False

            sage: (4).is_power()
            doctest:...: DeprecationWarning: is_power is deprecated. Please use is_perfect_power instead.
            See http://trac.sagemath.org/12116 for details.
            True

        TESTS:

        This is a test to make sure we work around a bug in GMP, see
        :trac:`4612`.

        ::

            sage: [ -a for a in srange(100) if not (-a^3).is_perfect_power() ]
            []
        """
        cdef mpz_t tmp
        cdef int res
        if mpz_sgn(self.value) < 0:
            if mpz_cmp_si(self.value, -1) == 0:
                return True
            mpz_init(tmp)
            mpz_neg(tmp, self.value)
            while mpz_perfect_square_p(tmp):
                mpz_sqrt(tmp, tmp)
            res = mpz_perfect_power_p(tmp)
            mpz_clear(tmp)
            if res:
                return True
            else:
                return False
        return mpz_perfect_power_p(self.value)

    def is_norm(self, K, element=False, proof=True):
        r"""
        See ``QQ(self).is_norm()``.

        EXAMPLES::

            sage: K = NumberField(x^2 - 2, 'beta')
            sage: n = 4
            sage: n.is_norm(K)
            True
            sage: 5.is_norm(K)
            False
            sage: 7.is_norm(QQ)
            True
            sage: n.is_norm(K, element=True)
            (True, -4*beta + 6)
            sage: n.is_norm(K, element=True)[1].norm()
            4
            sage: n = 5
            sage: n.is_norm(K, element=True)
            (False, None)
            sage: n = 7
            sage: n.is_norm(QQ, element=True)
            (True, 7)

        """
        from sage.rings.rational_field import QQ
        return QQ(self).is_norm(K, element=element, proof=proof)

    def _bnfisnorm(self, K, proof=True, extra_primes=0):
        r"""
        See ``QQ(self)._bnfisnorm()``.

        EXAMPLES::

            sage: 3._bnfisnorm(QuadraticField(-1, 'i'))
            (1, 3)
            sage: 7._bnfisnorm(CyclotomicField(7))
            (-zeta7 + 1, 1)            # 64-bit
            (-zeta7^5 + zeta7^4, 1)    # 32-bit
        """
        from sage.rings.rational_field import QQ
        return QQ(self)._bnfisnorm(K, proof=proof, extra_primes=extra_primes)


    def jacobi(self, b):
        r"""
        Calculate the Jacobi symbol `\left(\frac{self}{b}\right)`.

        EXAMPLES::

            sage: z = -1
            sage: z.jacobi(17)
            1
            sage: z.jacobi(19)
            -1
            sage: z.jacobi(17*19)
            -1
            sage: (2).jacobi(17)
            1
            sage: (3).jacobi(19)
            -1
            sage: (6).jacobi(17*19)
            -1
            sage: (6).jacobi(33)
            0
            sage: a = 3; b = 7
            sage: a.jacobi(b) == -b.jacobi(a)
            True
        """
        cdef long tmp
        if PY_TYPE_CHECK(b, int):
            tmp = b
            if (tmp & 1) == 0:
                raise ValueError, "Jacobi symbol not defined for even b."
            return mpz_kronecker_si(self.value, tmp)
        if not PY_TYPE_CHECK(b, Integer):
            b = Integer(b)
        if mpz_even_p((<Integer>b).value):
            raise ValueError, "Jacobi symbol not defined for even b."
        return mpz_jacobi(self.value, (<Integer>b).value)

    def kronecker(self, b):
        r"""
        Calculate the Kronecker symbol `\left(\frac{self}{b}\right)`
        with the Kronecker extension `(self/2)=(2/self)` when `self` is odd,
        or `(self/2)=0` when `self` is even.

        EXAMPLES::

            sage: z = 5
            sage: z.kronecker(41)
            1
            sage: z.kronecker(43)
            -1
            sage: z.kronecker(8)
            -1
            sage: z.kronecker(15)
            0
            sage: a = 2; b = 5
            sage: a.kronecker(b) == b.kronecker(a)
            True
        """
        if PY_TYPE_CHECK(b, int):
            return mpz_kronecker_si(self.value, b)
        if not PY_TYPE_CHECK(b, Integer):
            b = Integer(b)
        return mpz_kronecker(self.value, (<Integer>b).value)

    def radical(self, *args, **kwds):
        r"""
        Return the product of the prime divisors of self. Computing
        the radical of zero gives an error.

        EXAMPLES::

            sage: Integer(10).radical()
            10
            sage: Integer(20).radical()
            10
            sage: Integer(-100).radical()
            10
            sage: Integer(0).radical()
            Traceback (most recent call last):
            ...
            ArithmeticError: Radical of 0 not defined.
        """
        if self.is_zero():
            raise ArithmeticError, "Radical of 0 not defined."
        return self.factor(*args, **kwds).radical_value()

    def squarefree_part(self, long bound=-1):
        r"""
        Return the square free part of `x` (=self), i.e., the unique integer
        `z` that `x = z y^2`, with `y^2` a perfect square and `z` square-free.

        Use ``self.radical()`` for the product of the primes that divide self.

        If self is 0, just returns 0.

        EXAMPLES::

            sage: squarefree_part(100)
            1
            sage: squarefree_part(12)
            3
            sage: squarefree_part(17*37*37)
            17
            sage: squarefree_part(-17*32)
            -34
            sage: squarefree_part(1)
            1
            sage: squarefree_part(-1)
            -1
            sage: squarefree_part(-2)
            -2
            sage: squarefree_part(-4)
            -1

        ::

            sage: a = 8 * 5^6 * 101^2
            sage: a.squarefree_part(bound=2).factor()
            2 * 5^6 * 101^2
            sage: a.squarefree_part(bound=5).factor()
            2 * 101^2
            sage: a.squarefree_part(bound=1000)
            2
            sage: a.squarefree_part(bound=2**14)
            2
            sage: a = 7^3 * next_prime(2^100)^2 * next_prime(2^200)
            sage: a / a.squarefree_part(bound=1000)
            49
        """
        cdef Integer z
        cdef long even_part, p, p2
        cdef char switch_p
        if mpz_sgn(self.value) == 0:
            return self
        if 0 <= bound < 2:
            return self
        elif 2 <= bound <= 10000:
            z = PY_NEW(Integer)
            even_part = mpz_scan1(self.value, 0)
            mpz_fdiv_q_2exp(z.value, self.value, even_part ^ (even_part&1))
            sig_on()
            if bound >= 3:
                while mpz_divisible_ui_p(z.value, 9):
                    mpz_divexact_ui(z.value, z.value, 9)
            if bound >= 5:
                while mpz_divisible_ui_p(z.value, 25):
                    mpz_divexact_ui(z.value, z.value, 25)
            for p from 7 <= p <= bound by 2:
                switch_p = p % 30
                if switch_p in [1, 7, 11, 13, 17, 19, 23, 29]:
                    p2 = p*p
                    while mpz_divisible_ui_p(z.value, p2):
                        mpz_divexact_ui(z.value, z.value, p2)
            sig_off()
            return z
        else:
            if bound == -1:
                F = self.factor()
            else:
                from sage.rings.factorint import factor_trial_division
                F = factor_trial_division(self,bound)
            n = one
            for pp, e in F:
                if e % 2 != 0:
                    n = n * pp
            return n * F.unit()

    def next_probable_prime(self):
        """
        Returns the next probable prime after self, as determined by PARI.

        EXAMPLES::

            sage: (-37).next_probable_prime()
            2
            sage: (100).next_probable_prime()
            101
            sage: (2^512).next_probable_prime()
            13407807929942597099574024998205846127479365820592393377723561443721764030073546976801874298166903427690031858186486050853753882811946569946433649006084171
            sage: 0.next_probable_prime()
            2
            sage: 126.next_probable_prime()
            127
            sage: 144168.next_probable_prime()
            144169
        """
        return Integer( self._pari_().nextprime(True) )

    def next_prime(self, proof=None):
        r"""
        Returns the next prime after self.

        INPUT:

        -  ``proof`` - bool or None (default: None, see
           proof.arithmetic or sage.structure.proof) Note that the global Sage
           default is proof=True

        EXAMPLES::

            sage: Integer(100).next_prime()
            101

        Use Proof = False, which is way faster::

            sage: b = (2^1024).next_prime(proof=False)

        ::

            sage: Integer(0).next_prime()
            2
            sage: Integer(1001).next_prime()
            1009
        """
        if mpz_cmp(self.value, PARI_PSEUDOPRIME_LIMIT) < 0:
            return Integer( self._pari_().nextprime(True) )
        if proof is None:
            from sage.structure.proof.proof import get_flag
            proof = get_flag(proof, "arithmetic")
        if not proof:
            return Integer( self._pari_().nextprime(True) )
        n = self
        if n % 2 == 0:
            n += 1
        else:
            n += 2
        while not n.is_prime():  # PARI's isprime() is provably correct
            n += 2
        return integer_ring.ZZ(n)

    def additive_order(self):
        """
        Return the additive order of self.

        EXAMPLES::

            sage: ZZ(0).additive_order()
            1
            sage: ZZ(1).additive_order()
            +Infinity
        """
        import sage.rings.infinity
        if self.is_zero():
            return one
        else:
            return sage.rings.infinity.infinity

    def multiplicative_order(self):
        r"""
        Return the multiplicative order of self.

        EXAMPLES::

            sage: ZZ(1).multiplicative_order()
            1
            sage: ZZ(-1).multiplicative_order()
            2
            sage: ZZ(0).multiplicative_order()
            +Infinity
            sage: ZZ(2).multiplicative_order()
            +Infinity
        """
        import sage.rings.infinity
        if  mpz_cmp_si(self.value, 1) == 0:
                return one
        elif mpz_cmp_si(self.value, -1) == 0:
                return Integer(2)
        else:
                return sage.rings.infinity.infinity

    def is_squarefree(self):
        """
        Returns True if this integer is not divisible by the square of any
        prime and False otherwise.

        EXAMPLES::

            sage: Integer(100).is_squarefree()
            False
            sage: Integer(102).is_squarefree()
            True
            sage: Integer(0).is_squarefree()
            False
        """
        return self._pari_().issquarefree()

    def _pari_(self):
        """
        Returns the PARI version of this integer.

        EXAMPLES::

            sage: n = 9390823
            sage: m = n._pari_(); m
            9390823
            sage: type(m)
            <type 'sage.libs.pari.gen.gen'>

        TESTS::

            sage: n = 10^10000000
            sage: m = n._pari_() ## crash from trac 875
            sage: m % 1234567
            1041334
        """
        return self._pari_c()

    cdef _pari_c(self):
        return pari.new_gen_from_mpz_t(self.value)

    def _interface_init_(self, I=None):
        """
        Return canonical string to coerce this integer to any other math
        software, i.e., just the string representation of this integer in
        base 10.

        EXAMPLES::

            sage: n = 9390823
            sage: n._interface_init_()
            '9390823'
        """
        return str(self)

    property __array_interface__:
        def __get__(self):
            """
            Used for NumPy conversion.

            EXAMPLES::

                sage: import numpy
                sage: numpy.array([1, 2, 3])
                array([1, 2, 3])
                sage: numpy.array([1, 2, 3]).dtype
                dtype('int32')                         # 32-bit
                dtype('int64')                         # 64-bit

                sage: numpy.array(2**40).dtype
                dtype('int64')
                sage: numpy.array(2**400).dtype
                dtype('O')

                sage: numpy.array([1,2,3,0.1]).dtype
                dtype('float64')
            """
            if mpz_fits_slong_p(self.value):
                return numpy_long_interface
            elif sizeof(long) == 4 and mpz_sizeinbase(self.value, 2) <= 63:
                return numpy_int64_interface
            else:
                return numpy_object_interface

    def _magma_init_(self, magma):
        """
        Return string that evaluates in Magma to this element.

        For small integers we just use base 10.  For large integers we use
        base 16, but use Magma's StringToInteger command, which (for no
        good reason) is much faster than 0x[string literal].  We only use
        base 16 for integers with at least 10000 binary digits, since e.g.,
        for a large list of small integers the overhead of calling
        StringToInteger can be a killer.

        EXAMPLES:
            sage: (117)._magma_init_(magma)           # optional - magma
            '117'

        Large integers use hex:
            sage: m = 3^(2^20)                        # optional - magma
            sage: s = m._magma_init_(magma)           # optional - magma
            sage: 'StringToInteger' in s              # optional - magma
            True
            sage: magma(m).sage() == m                # optional - magma
            True
        """
        if self.ndigits(2) > 10000:
            return 'StringToInteger("%s",16)'%self.str(16)
        else:
            return str(self)

    def _sage_input_(self, sib, coerced):
        r"""
        Produce an expression which will reproduce this value when
        evaluated.

        EXAMPLES::

            sage: sage_input(1, verify=True)
            # Verified
            1
            sage: sage_input(1, preparse=False)
            ZZ(1)
            sage: sage_input(-12435, verify=True)
            # Verified
            -12435
            sage: sage_input(0, verify=True)
            # Verified
            0
            sage: sage_input(-3^70, verify=True)
            # Verified
            -2503155504993241601315571986085849
            sage: sage_input(-37, preparse=False)
            -ZZ(37)
            sage: sage_input(-37 * polygen(ZZ), preparse=False)
            R = ZZ['x']
            x = R.gen()
            -37*x
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: (-314159)._sage_input_(SageInputBuilder(preparse=False), False)
            {unop:- {call: {atomic:ZZ}({atomic:314159})}}
            sage: (314159)._sage_input_(SageInputBuilder(preparse=False), True)
            {atomic:314159}
        """
        if coerced or sib.preparse():
            return sib.int(self)
        else:
            if self < 0:
                return -sib.name('ZZ')(sib.int(-self))
            else:
                return sib.name('ZZ')(sib.int(self))

    def sqrtrem(self):
        r"""
        Return (s, r) where s is the integer square root of self and
        r is the remainder such that `\text{self} = s^2 + r`.
        Raises ``ValueError`` if self is negative.

        EXAMPLES::

            sage: 25.sqrtrem()
            (5, 0)
            sage: 27.sqrtrem()
            (5, 2)
            sage: 0.sqrtrem()
            (0, 0)

        ::

            sage: Integer(-102).sqrtrem()
            Traceback (most recent call last):
            ...
            ValueError: square root of negative integer not defined.

        """
        if mpz_sgn(self.value) < 0:
            raise ValueError, "square root of negative integer not defined."
        cdef Integer s = PY_NEW(Integer)
        cdef Integer r  = PY_NEW(Integer)
        mpz_sqrtrem(s.value, r.value, self.value)
        return s, r

    def isqrt(self):
        r"""
        Returns the integer floor of the square root of self, or raises an
        ``ValueError`` if self is negative.

        EXAMPLE::

            sage: a = Integer(5)
            sage: a.isqrt()
            2

        ::

            sage: Integer(-102).isqrt()
            Traceback (most recent call last):
            ...
            ValueError: square root of negative integer not defined.
        """
        if mpz_sgn(self.value) < 0:
            raise ValueError, "square root of negative integer not defined."

        cdef Integer x = PY_NEW(Integer)

        sig_on()
        mpz_sqrt(x.value, self.value)
        sig_off()

        return x

    def sqrt_approx(self, prec=None, all=False):
        """
        EXAMPLES::

            sage: 5.sqrt_approx(prec=200)
            doctest:...: DeprecationWarning: This function is deprecated.  Use sqrt with a given number of bits of precision instead.
            See http://trac.sagemath.org/10107 for details.
            2.2360679774997896964091736687312762354406183596115257242709
            sage: 5.sqrt_approx()
            2.23606797749979
            sage: 4.sqrt_approx()
            2
        """
        from sage.misc.superseded import deprecation
        deprecation(10107, "This function is deprecated.  Use sqrt with a given number of bits of precision instead.")
        try:
            return self.sqrt(extend=False,all=all)
        except ValueError:
            pass
        if prec is None:
            prec = max(53, 2*(mpz_sizeinbase(self.value, 2)+2))
        return self.sqrt(prec=prec, all=all)

    def sqrt(self, prec=None, extend=True, all=False):
        """
        The square root function.

        INPUT:

        -  ``prec`` - integer (default: None): if None, returns
           an exact square root; otherwise returns a numerical square root if
           necessary, to the given bits of precision.

        -  ``extend`` - bool (default: True); if True, return a
           square root in an extension ring, if necessary. Otherwise, raise a
           ValueError if the square is not in the base ring. Ignored if prec
           is not None.

        -  ``all`` - bool (default: False); if True, return all
           square roots of self, instead of just one.

        EXAMPLES::

            sage: Integer(144).sqrt()
            12
            sage: sqrt(Integer(144))
            12
            sage: Integer(102).sqrt()
            sqrt(102)

        ::

            sage: n = 2
            sage: n.sqrt(all=True)
            [sqrt(2), -sqrt(2)]
            sage: n.sqrt(prec=10)
            1.4
            sage: n.sqrt(prec=100)
            1.4142135623730950488016887242
            sage: n.sqrt(prec=100,all=True)
            [1.4142135623730950488016887242, -1.4142135623730950488016887242]
            sage: n.sqrt(extend=False)
            Traceback (most recent call last):
            ...
            ValueError: square root of 2 not an integer
            sage: Integer(144).sqrt(all=True)
            [12, -12]
            sage: Integer(0).sqrt(all=True)
            [0]
            sage: type(Integer(5).sqrt())
            <type 'sage.symbolic.expression.Expression'>
            sage: type(Integer(5).sqrt(prec=53))
            <type 'sage.rings.real_mpfr.RealNumber'>
            sage: type(Integer(-5).sqrt(prec=53))
            <type 'sage.rings.complex_number.ComplexNumber'>

        TESTS:

        Check that :trac:`9466` is fixed::

            sage: 3.sqrt(extend=False, all=True)
            []
        """
        if mpz_sgn(self.value) == 0:
            return [self] if all else self

        if mpz_sgn(self.value) < 0:
            if not extend:
                raise ValueError, "square root of negative number not an integer"
            from sage.functions.other import _do_sqrt
            return _do_sqrt(self, prec=prec, all=all)

        cdef int non_square
        cdef Integer z = PY_NEW(Integer)
        cdef mpz_t tmp
        sig_on()
        mpz_init(tmp)
        mpz_sqrtrem(z.value, tmp, self.value)
        non_square = mpz_sgn(tmp) != 0
        mpz_clear(tmp)
        sig_off()

        if non_square:
            if not extend:
                if not all:
                   raise ValueError, "square root of %s not an integer"%self
                else:
                    return []
            from sage.functions.other import _do_sqrt
            return _do_sqrt(self, prec=prec, all=all)

        if prec:
            from sage.functions.other import _do_sqrt
            return _do_sqrt(self, prec=prec, all=all)

        if all:
           return [z, -z]
        return z

    def _xgcd(self, Integer n, bint minimal=0):
        r"""
        Computes extended gcd of self and `n`.

        INPUT:

        -  ``self`` - integer

        -  ``n`` - integer

        -  ``minimal`` - boolean (default false), whether to
           compute minimal cofactors (see below)

        OUTPUT:

        A triple (g, s, t), where `g` is the non-negative gcd of self
        and `n`, and `s` and `t` are cofactors satisfying the Bezout identity

        .. math::

             g = s \cdot \mbox{\rm self} + t \cdot n.

        .. note::

           If minimal is False, then there is no guarantee that the returned
           cofactors will be minimal in any sense; the only guarantee is that
           the Bezout identity will be satisfied (see examples below).

           If minimal is True, the cofactors will satisfy the following
           conditions. If either self or `n` are zero, the trivial solution
           is returned. If both self and `n` are nonzero, the function returns
           the unique solution such that `0 \leq s < |n|/g` (which then must
           also satisfy `0 \leq |t| \leq |\mbox{\rm self}|/g`).

        EXAMPLES::

            sage: Integer(5)._xgcd(7)
            (1, 3, -2)
            sage: 5*3 + 7*-2
            1
            sage: g,s,t = Integer(58526524056)._xgcd(101294172798)
            sage: g
            22544886
            sage: 58526524056 * s + 101294172798 * t
            22544886

        Without minimality guarantees, weird things can happen::

            sage: Integer(3)._xgcd(21)
            (3, 1, 0)
            sage: Integer(3)._xgcd(24)
            (3, 1, 0)
            sage: Integer(3)._xgcd(48)
            (3, 1, 0)

        Weirdness disappears with minimal option::

            sage: Integer(3)._xgcd(21, minimal=True)
            (3, 1, 0)
            sage: Integer(3)._xgcd(24, minimal=True)
            (3, 1, 0)
            sage: Integer(3)._xgcd(48, minimal=True)
            (3, 1, 0)
            sage: Integer(21)._xgcd(3, minimal=True)
            (3, 0, 1)

        Try minimal option with various edge cases::

            sage: Integer(5)._xgcd(0, minimal=True)
            (5, 1, 0)
            sage: Integer(-5)._xgcd(0, minimal=True)
            (5, -1, 0)
            sage: Integer(0)._xgcd(5, minimal=True)
            (5, 0, 1)
            sage: Integer(0)._xgcd(-5, minimal=True)
            (5, 0, -1)
            sage: Integer(0)._xgcd(0, minimal=True)
            (0, 1, 0)

        Exhaustive tests, checking minimality conditions::

            sage: for a in srange(-20, 20):
            ...     for b in srange(-20, 20):
            ...       if a == 0 or b == 0: continue
            ...       g, s, t = a._xgcd(b)
            ...       assert g > 0
            ...       assert a % g == 0 and b % g == 0
            ...       assert a*s + b*t == g
            ...       g, s, t = a._xgcd(b, minimal=True)
            ...       assert g > 0
            ...       assert a % g == 0 and b % g == 0
            ...       assert a*s + b*t == g
            ...       assert s >= 0 and s < abs(b)/g
            ...       assert abs(t) <= abs(a)/g

        AUTHORS:

        - David Harvey (2007-12-26): added minimality option
        """
        cdef Integer g = PY_NEW(Integer)
        cdef Integer s = PY_NEW(Integer)
        cdef Integer t = PY_NEW(Integer)

        sig_on()
        mpz_gcdext(g.value, s.value, t.value, self.value, n.value)
        sig_off()

        # Note: the GMP documentation for mpz_gcdext (or mpn_gcdext for that
        # matter) makes absolutely no claims about any minimality conditions
        # satisfied by the returned cofactors. They guarantee a non-negative
        # gcd, but that's it. So we have to do some work ourselves.

        if not minimal:
            return g, s, t

        # handle degenerate cases n == 0 and self == 0

        if not mpz_sgn(n.value):
            mpz_set_ui(t.value, 0)
            mpz_abs(g.value, self.value)
            mpz_set_si(s.value, 1 if mpz_sgn(self.value) >= 0 else -1)
            return g, s, t

        if not mpz_sgn(self.value):
            mpz_set_ui(s.value, 0)
            mpz_abs(g.value, n.value)
            mpz_set_si(t.value, 1 if mpz_sgn(n.value) >= 0 else -1)
            return g, s, t

        # both n and self are nonzero, so we need to do a division and
        # make the appropriate adjustment

        cdef mpz_t u1, u2
        mpz_init(u1)
        mpz_init(u2)
        mpz_divexact(u1, n.value, g.value)
        mpz_divexact(u2, self.value, g.value)
        if mpz_sgn(u1) > 0:
            mpz_fdiv_qr(u1, s.value, s.value, u1)
        else:
            mpz_cdiv_qr(u1, s.value, s.value, u1)
        mpz_addmul(t.value, u1, u2)
        mpz_clear(u2)
        mpz_clear(u1)

        return g, s, t

    cpdef _shift_helper(Integer self, y, int sign):
        """
        Function used to compute left and right shifts of integers.
        Shifts self y bits to the left if sign is 1, and to the right
        if sign is -1.

        WARNING: This function does no error checking. In particular,
        it assumes that sign is either 1 or -1,

        EXAMPLES::

            sage: n = 1234
            sage: factor(n)
            2 * 617
            sage: n._shift_helper(1, 1)
            2468
            sage: n._shift_helper(1, -1)
            617
            sage: n._shift_helper(100, 1)
            1564280840681635081446931755433984
            sage: n._shift_helper(100, -1)
            0
            sage: n._shift_helper(-100, 1)
            0
            sage: n._shift_helper(-100, -1)
            1564280840681635081446931755433984
        """
        cdef long n

        if PyInt_CheckExact(y):
            # For a Python int, we can just use the Python/C API.
            n = PyInt_AS_LONG(y)
        else:
            # If it's not already an Integer, try to convert it.
            if not PY_TYPE_CHECK(y, Integer):
                try:
                    y = Integer(y)
                except TypeError:
                    raise TypeError, "unsupported operands for %s: %s, %s"%(("<<" if sign == 1 else ">>"), self, y)
                except ValueError:
                    return bin_op(self, y, operator.lshift if sign == 1 else operator.rshift)

            # If y wasn't a Python int, it's now an Integer, so set n
            # accordingly.
            if mpz_fits_slong_p((<Integer>y).value):
                n = mpz_get_si((<Integer>y).value)
            elif sign * mpz_sgn((<Integer>y).value) < 0:
                # Doesn't fit in a long so shifting to the right by
                # this much will be 0.
                return PY_NEW(Integer)
            else:
                # Doesn't fit in a long so shifting to the left by
                # this much will raise appropriate overflow error
                n = y

        # Decide which way we're shifting
        n *= sign

        # Now finally call into MPIR to do the shifting.
        cdef Integer z = PY_NEW(Integer)
        if n < 0:
            mpz_fdiv_q_2exp(z.value, self.value, -n)
        else:
            mpz_mul_2exp(z.value, self.value, n)
        return z

    def __lshift__(x, y):
        """
        Shift x to the left by y bits.

        EXAMPLES::

            sage: 32 << 2
            128
            sage: 32 << int(2)
            128
            sage: int(32) << 2
            128
            sage: 1 << 2.5
            Traceback (most recent call last):
            ...
            TypeError: unsupported operands for <<: 1, 2.5000...

            sage: 32 << (4/2)
            128

        A negative shift to the left is treated as a right shift::

            sage: 128 << -2
            32
            sage: 128 << (-2^100)
            0
        """
        # note that x need not be self -- int(3) << ZZ(2) will
        # dispatch this function
        if not PY_TYPE_CHECK(x, Integer):
            return x << int(y)
        return (<Integer>x)._shift_helper(y, 1)

    def __rshift__(x, y):
        """
        Shift x to the right by y bits.

        EXAMPLES::

            sage: 32 >> 2
            8
            sage: 32 >> int(2)
            8
            sage: int(32) >> 2
            8
            sage: 1 >> 2.5
            Traceback (most recent call last):
            ...
            TypeError: unsupported operands for >>: 1, 2.5000...
            sage: 10^5 >> 10^100
            0

        A negative shift to the right is treated as a left shift::

            sage: 8 >> -2
            32
        """
        # note that x need not be self -- int(3) >> ZZ(2) will
        # dispatch this function
        if not PY_TYPE_CHECK(x, Integer):
            return x >> int(y)
        return (<Integer>x)._shift_helper(y, -1)

    cdef _and(Integer self, Integer other):
        cdef Integer x = PY_NEW(Integer)
        mpz_and(x.value, self.value, other.value)
        return x

    def __and__(x, y):
        """
        Return the bitwise and two integers.

        EXAMPLES::

            sage: n = Integer(6);  m = Integer(2)
            sage: n & m
            2
            sage: n.__and__(m)
            2
        """
        if PY_TYPE_CHECK(x, Integer) and PY_TYPE_CHECK(y, Integer):
            return (<Integer>x)._and(y)
        return bin_op(x, y, operator.and_)

    cdef _or(Integer self, Integer other):
        cdef Integer x = PY_NEW(Integer)
        mpz_ior(x.value, self.value, other.value)
        return x

    def __or__(x, y):
        """
        Return the bitwise or of the integers x and y.

        EXAMPLES::

            sage: n = 8; m = 4
            sage: n.__or__(m)
            12
        """
        if PY_TYPE_CHECK(x, Integer) and PY_TYPE_CHECK(y, Integer):
            return (<Integer>x)._or(y)
        return bin_op(x, y, operator.or_)

    def __invert__(self):
        """
        Return the multiplicative inverse of self, as a rational number.

        EXAMPLE::

            sage: n = 10
            sage: 1/n
            1/10
            sage: n.__invert__()
            1/10
        """
        return one / self

    def inverse_of_unit(self):
        """
        Return inverse of self if self is a unit in the integers, i.e.,
        self is -1 or 1. Otherwise, raise a ZeroDivisionError.

        EXAMPLES::

            sage: (1).inverse_of_unit()
            1
            sage: (-1).inverse_of_unit()
            -1
            sage: 5.inverse_of_unit()
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Inverse does not exist.
            sage: 0.inverse_of_unit()
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Inverse does not exist.
        """
        if mpz_cmpabs_ui(self.value, 1) == 0:
            return self
        else:
            raise ZeroDivisionError, "Inverse does not exist."

    def inverse_mod(self, n):
        """
        Returns the inverse of self modulo `n`, if this inverse exists.
        Otherwise, raises a ``ZeroDivisionError`` exception.

        INPUT:

        -  ``self`` - Integer

        -  ``n`` - Integer, or ideal of integer ring

        OUTPUT:

        -  ``x`` - Integer such that x\*self = 1 (mod m), or
           raises ZeroDivisionError.

        IMPLEMENTATION:

        Call the mpz_invert GMP library function.

        EXAMPLES::

            sage: a = Integer(189)
            sage: a.inverse_mod(10000)
            4709
            sage: a.inverse_mod(-10000)
            4709
            sage: a.inverse_mod(1890)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Inverse does not exist.
            sage: a = Integer(19)**100000
            sage: b = a*a
            sage: c = a.inverse_mod(b)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Inverse does not exist.

        We check that #10625 is fixed::

            sage: ZZ(2).inverse_mod(ZZ.ideal(3))
            2

        We check that #9955 is fixed::

            sage: Rational(3)%Rational(-1)
            0
        """
        cdef int r
        if isinstance(n, sage.rings.ideal.Ideal_pid) and n.ring() == the_integer_ring:
            n = n.gen()
        cdef Integer m = as_Integer(n)
        cdef Integer ans = <Integer>PY_NEW(Integer)
        if mpz_cmpabs_ui(m.value, 1) == 0:
            return zero
        sig_on()
        r = mpz_invert(ans.value, self.value, m.value)
        sig_off()
        if r == 0:
            raise ZeroDivisionError, "Inverse does not exist."
        return ans

    def gcd(self, n):
        """
        Return the greatest common divisor of self and `n`.

        EXAMPLE::

            sage: gcd(-1,1)
            1
            sage: gcd(0,1)
            1
            sage: gcd(0,0)
            0
            sage: gcd(2,2^6)
            2
            sage: gcd(21,2^6)
            1
        """
        if not isinstance(n, Integer) and not isinstance(n, int):
            left, right = canonical_coercion(self, n)
            return left.gcd(right)
        cdef Integer m = as_Integer(n)
        cdef Integer g = PY_NEW(Integer)
        sig_on()
        mpz_gcd(g.value, self.value, m.value)
        sig_off()
        return g

    def crt(self, y, m, n):
        """
        Return the unique integer between `0` and `mn` that is congruent to
        the integer modulo `m` and to `y` modulo `n`. We assume that `m` and
        `n` are coprime.

        EXAMPLES::

            sage: n = 17
            sage: m = n.crt(5, 23, 11); m
            247
            sage: m%23
            17
            sage: m%11
            5
        """
        cdef object g, s, t
        cdef Integer _y, _m, _n
        _y = Integer(y); _m = Integer(m); _n = Integer(n)
        g, s, t = _m.xgcd(_n)
        if not g.is_one():
            raise ArithmeticError, "CRT requires that gcd of moduli is 1."
        # Now s*m + t*n = 1, so the answer is x + (y-x)*s*m, where x=self.
        return (self + (_y - self) * s * _m) % (_m * _n)

    def test_bit(self, long index):
        r"""
        Return the bit at ``index``.

        If the index is negative, returns 0.

        Although internally a sign-magnitude representation is used
        for integers, this method pretends to use a two's complement
        representation.  This is illustrated with a negative integer
        below.

        EXAMPLES::

            sage: w = 6
            sage: w.str(2)
            '110'
            sage: w.test_bit(2)
            1
            sage: w.test_bit(-1)
            0
            sage: x = -20
            sage: x.str(2)
            '-10100'
            sage: x.test_bit(4)
            0
            sage: x.test_bit(5)
            1
            sage: x.test_bit(6)
            1
        """
        if index < 0:
            return 0
        else:
            return mpz_tstbit(self.value, index)

    def popcount(self):
        """
        Return the number of 1 bits in the binary representation.
        If self<0, we return Infinity.

        EXAMPLES::

            sage: n = 123
            sage: n.str(2)
            '1111011'
            sage: n.popcount()
            6

            sage: n = -17
            sage: n.popcount()
            +Infinity
        """
        if self<0:
            return sage.rings.infinity.Infinity
        return int(mpz_popcount(self.value))


    def conjugate(self):
        """
        Return the complex conjugate of this integer, which is the
        integer itself.

        EXAMPLES:
            sage: n = 205
            sage: n.conjugate()
            205
        """
        return self

    def binomial(self, m, algorithm='mpir'):
        """
        Return the binomial coefficient "self choose m".

        INPUT:

            - ``m`` -- an integer

            - ``algorithm`` -- ``'mpir'`` (default) or ``'pari'``; ``'mpir'`` is
              faster for small ``m``, and ``'pari'`` tends to be faster for
              large ``m``

        OUTPUT:

            - integer

        EXAMPLES::

            sage: 10.binomial(2)
            45
            sage: 10.binomial(2, algorithm='pari')
            45
            sage: 10.binomial(-2)
            0
        """
        cdef Integer x
        if m < 0:
            return the_integer_ring.zero()
        if algorithm == 'mpir':
            x = PY_NEW(Integer)
            sig_on()
            mpz_bin_ui(x.value, self.value, m)
            sig_off()
            return x
        elif algorithm == 'pari':
            return the_integer_ring(self._pari_().binomial(m))
        else:
            raise ValueError("algorithm must be one of: 'pari', 'mpir'")


cpdef LCM_list(v):
    """
    Return the LCM of a list v of integers. Elements of v are converted
    to Sage integers if they aren't already.

    This function is used, e.g., by rings/arith.py

    INPUT:

    -  ``v`` - list or tuple

    OUTPUT: integer

    EXAMPLES::

        sage: from sage.rings.integer import LCM_list
        sage: w = LCM_list([3,9,30]); w
        90
        sage: type(w)
        <type 'sage.rings.integer.Integer'>

    The inputs are converted to Sage integers.

    ::

        sage: w = LCM_list([int(3), int(9), '30']); w
        90
        sage: type(w)
        <type 'sage.rings.integer.Integer'>
    """
    cdef int i, n = len(v)
    cdef Integer z = <Integer>PY_NEW(Integer)

    for i from 0 <= i < n:
        if not isinstance(v[i], Integer):
            if not isinstance(v, list):
                v = list(v)
            v[i] = Integer(v[i])

    if n == 0:
        return one
    elif n == 1:
        return v[0].abs()

    sig_on()
    mpz_lcm(z.value, (<Integer>v[0]).value, (<Integer>v[1]).value)
    for i from 2 <= i < n:
        mpz_lcm(z.value, z.value, (<Integer>v[i]).value)
    sig_off()

    return z

def GCD_list(v):
    r"""
    Return the GCD of a list v of integers. Elements of v are converted
    to Sage integers if they aren't already.

    This function is used, e.g., by rings/arith.py

    INPUT:

    -  ``v`` - list or tuple

    OUTPUT: integer

    EXAMPLES::

        sage: from sage.rings.integer import GCD_list
        sage: w = GCD_list([3,9,30]); w
        3
        sage: type(w)
        <type 'sage.rings.integer.Integer'>

    Check that the bug reported in trac #3118 has been fixed::

        sage: sage.rings.integer.GCD_list([2,2,3])
        1

    The inputs are converted to Sage integers.

    ::

        sage: w = GCD_list([int(3), int(9), '30']); w
        3
        sage: type(w)
        <type 'sage.rings.integer.Integer'>
    """
    cdef int i, n = len(v)
    cdef Integer z = <Integer>PY_NEW(Integer)

    for i from 0 <= i < n:
        if not isinstance(v[i], Integer):
            if not isinstance(v, list):
                v = list(v)
            v[i] = Integer(v[i])

    if n == 0:
        return one
    elif n == 1:
        return v[0].abs()

    sig_on()
    mpz_gcd(z.value, (<Integer>v[0]).value, (<Integer>v[1]).value)
    for i from 2 <= i < n:
        if mpz_cmp_ui(z.value, 1) == 0:
            break
        mpz_gcd(z.value, z.value, (<Integer>v[i]).value)
    sig_off()

    return z

def make_integer(s):
    """
    Create a Sage integer from the base-32 Python *string* s. This is
    used in unpickling integers.

    EXAMPLES::

        sage: from sage.rings.integer import make_integer
        sage: make_integer('-29')
        -73
        sage: make_integer(29)
        Traceback (most recent call last):
        ...
        TypeError: expected string or Unicode object, sage.rings.integer.Integer found
    """
    cdef Integer r = PY_NEW(Integer)
    r._reduce_set(s)
    return r

cdef class int_to_Z(Morphism):
    def __init__(self):
        import integer_ring
        import sage.categories.homset
        from sage.structure.parent import Set_PythonType
        Morphism.__init__(self, sage.categories.homset.Hom(Set_PythonType(int), integer_ring.ZZ))
    cpdef Element _call_(self, a):
        cdef Integer r
        r = <Integer>PY_NEW(Integer)
        mpz_set_si(r.value, PyInt_AS_LONG(a))
        return r
    def _repr_type(self):
        return "Native"

cdef class long_to_Z(Morphism):
    """
    EXAMPLES::

        sage: f = ZZ.coerce_map_from(long); f
        Native morphism:
          From: Set of Python objects of type 'long'
          To:   Integer Ring
        sage: f(1rL)
        1
        sage: f(-10000000000000000000001r)
        -10000000000000000000001
    """
    def __init__(self):
        import integer_ring
        import sage.categories.homset
        from sage.structure.parent import Set_PythonType
        Morphism.__init__(self, sage.categories.homset.Hom(Set_PythonType(long), integer_ring.ZZ))
    cpdef Element _call_(self, a):
        cdef Integer r
        r = <Integer>PY_NEW(Integer)
        mpz_set_pylong(r.value, a)
        return r
    def _repr_type(self):
        return "Native"

############### INTEGER CREATION CODE #####################

# We need a couple of internal GMP datatypes.

# This may be potentially very dangerous as it reaches
# deeply into the internal structure of GMP which may not
# be consistent across future versions of GMP.
# See extensive note in the fast_tp_new() function below.

include "sage/ext/python_rich_object.pxi"
cdef extern from "gmp.h":
    ctypedef long mp_limb_t
    ctypedef mp_limb_t* mp_ptr #"mp_ptr"

    # We allocate _mp_d directly (mpz_t is typedef of this in GMP)
    ctypedef struct __mpz_struct "__mpz_struct":
        int _mp_alloc
        int _mp_size
        mp_ptr _mp_d

    # sets the three free, alloc, and realloc function pointers to the
    # memory management functions set in GMP. Accepts NULL pointer.
    # Potentially dangerous if changed by calling
    # mp_set_memory_functions again after we initialized this module.
    void mp_get_memory_functions (void *(**alloc) (size_t), void *(**realloc)(void *, size_t, size_t), void (**free) (void *, size_t))

    # GMP's configuration of how many Bits are stuffed into a limb
    cdef int GMP_LIMB_BITS

# This variable holds the size of any Integer object in bytes.
cdef int sizeof_Integer

# We use a global Integer element to steal all the references
# from.  DO NOT INITIALIZE IT AGAIN and DO NOT REFERENCE IT!
cdef Integer global_dummy_Integer
global_dummy_Integer = Integer()

# Accessing the .value attribute of an Integer object causes Cython to
# refcount it. This is problematic, because that causes overhead and
# more importantly an infinite loop in the destructor. If you refcount
# in the destructor and the refcount reaches zero (which is true
# every time) the destructor is called.
#
# To avoid this we calculate the byte offset of the value member and
# remember it in this variable.
#
# Eventually this may be rendered obsolete by a change in Cython allowing
# non-reference counted extension types.
cdef long mpz_t_offset
mpz_t_offset_python = None


# stores the GMP alloc function
cdef void * (* mpz_alloc)(size_t)

# stores the GMP free function
cdef void (* mpz_free)(void *, size_t)

# A global  pool for performance when integers are rapidly created and destroyed.
# It operates on the following principles:
#
# - The pool starts out empty.
# - When an new integer is needed, one from the pool is returned
#   if available, otherwise a new Integer object is created
# - When an integer is collected, it will add it to the pool
#   if there is room, otherwise it will be deallocated.


cdef int integer_pool_size = 100

cdef PyObject** integer_pool
cdef int integer_pool_count = 0

# used for profiling the pool
cdef int total_alloc = 0
cdef int use_pool = 0


cdef PyObject* fast_tp_new(PyTypeObject *t, PyObject *a, PyObject *k):

    global integer_pool, integer_pool_count, total_alloc, use_pool

    cdef PyObject* new

    # for profiling pool usage
    # total_alloc += 1

    # If there is a ready integer in the pool, we will
    # decrement the counter and return that.

    if integer_pool_count > 0:

        # for profiling pool usage
        # use_pool += 1

        integer_pool_count -= 1
        new = <PyObject *> integer_pool[integer_pool_count]

    # Otherwise, we have to create one.

    else:

        # allocate enough room for the Integer, sizeof_Integer is
        # sizeof(Integer). The use of PyObject_MALLOC directly
        # assumes that Integers are not garbage collected, i.e.
        # they do not possess references to other Python
        # objects (as indicated by the Py_TPFLAGS_HAVE_GC flag).
        # See below for a more detailed description.
        new = PyObject_MALLOC( sizeof_Integer )

        # Now set every member as set in z, the global dummy Integer
        # created before this tp_new started to operate.

        memcpy(new, (<void*>global_dummy_Integer), sizeof_Integer )

        # We take the address 'new' and move mpz_t_offset bytes (chars)
        # to the address of 'value'. We treat that address as a pointer
        # to a mpz_t struct and allocate memory for the _mp_d element of
        # that struct. We allocate one limb.
        #
        # What is done here is potentially very dangerous as it reaches
        # deeply into the internal structure of GMP. Consequently things
        # may break if a new release of GMP changes some internals. To
        # emphasize this, this is what the GMP manual has to say about
        # the documentation for the struct we are using:
        #
        #  "This chapter is provided only for informational purposes and the
        #  various internals described here may change in future GMP releases.
        #  Applications expecting to be compatible with future releases should use
        #  only the documented interfaces described in previous chapters."
        #
        # If this line is used Sage is not such an application.
        #
        # The clean version of the following line is:
        #
        #  mpz_init( <mpz_t>(<char *>new + mpz_t_offset) )
        #
        # We save time both by avoiding an extra function call and
        # because the rest of the mpz struct was already initialized
        # fully using the memcpy above.

        (<__mpz_struct *>( <char *>new + mpz_t_offset) )._mp_d = <mp_ptr>mpz_alloc(GMP_LIMB_BITS >> 3)

    # This line is only needed if Python is compiled in debugging mode
    # './configure --with-pydebug' or SAGE_DEBUG=yes. If that is the
    # case a Python object has a bunch of debugging fields which are
    # initialized with this macro.

    if_Py_TRACE_REFS_then_PyObject_INIT\
        (new, (<PyObject*>global_dummy_Integer).ob_type)

    # The global_dummy_Integer may have a reference count larger than
    # one, but it is expected that newly created objects have a
    # reference count of one. This is potentially unneeded if
    # everybody plays nice, because the gobal_dummy_Integer has only
    # one reference in that case.

    # Objects from the pool have reference count zero, so this
    # needs to be set in this case.

    new.ob_refcnt = 1

    return new

cdef void fast_tp_dealloc(PyObject* o):

    # If there is room in the pool for a used integer object,
    # then put it in rather than deallocating it.

    global integer_pool, integer_pool_count

    if integer_pool_count < integer_pool_size:

        # Here we free any extra memory used by the mpz_t by
        # setting it to a single limb.
        if (<__mpz_struct *>( <char *>o + mpz_t_offset))._mp_alloc > 10:
            _mpz_realloc(<mpz_t *>( <char *>o + mpz_t_offset), 10)

        # It's cheap to zero out an integer, so do it here.
        (<__mpz_struct *>( <char *>o + mpz_t_offset))._mp_size = 0

        # And add it to the pool.
        integer_pool[integer_pool_count] = o
        integer_pool_count += 1
        return

    # Again, we move to the mpz_t and clear it. See above, why this is evil.
    # The clean version of this line would be:
    #   mpz_clear(<mpz_t>(<char *>o + mpz_t_offset))

    mpz_free((<__mpz_struct *>( <char *>o + mpz_t_offset) )._mp_d, 0)

    # Free the object. This assumes that Py_TPFLAGS_HAVE_GC is not
    # set. If it was set another free function would need to be
    # called.

    PyObject_FREE(o)

from sage.misc.allocator cimport hook_tp_functions
cdef hook_fast_tp_functions():
    """
    Initialize the fast integer creation functions.
    """
    global global_dummy_Integer, mpz_t_offset, sizeof_Integer, integer_pool

    integer_pool = <PyObject**>sage_malloc(integer_pool_size * sizeof(PyObject*))

    cdef PyObject* o
    o = <PyObject *>global_dummy_Integer

    # calculate the offset of the GMP mpz_t to avoid casting to/from
    # an Integer which includes reference counting. Reference counting
    # is bad in constructors and destructors as it potentially calls
    # the destructor.
    # Eventually this may be rendered obsolete by a change in Cython allowing
    # non-reference counted extension types.
    mpz_t_offset = <char *>(&global_dummy_Integer.value) - <char *>o
    global mpz_t_offset_python
    mpz_t_offset_python = mpz_t_offset

    # store how much memory needs to be allocated for an Integer.
    sizeof_Integer = (<RichPyTypeObject *>o.ob_type).tp_basicsize

    # get the functions to do memory management for the GMP elements
    # WARNING: if the memory management functions are changed after
    # this initialisation, we are/you are doomed.

    mp_get_memory_functions(&mpz_alloc, NULL, &mpz_free)

    # Finally replace the functions called when an Integer needs
    # to be constructed/destructed.
    hook_tp_functions(global_dummy_Integer, NULL, &fast_tp_new, NULL, &fast_tp_dealloc, False)

cdef integer(x):
    if PY_TYPE_CHECK(x, Integer):
        return x
    return Integer(x)


def free_integer_pool():
    cdef int i
    cdef PyObject *o

    global integer_pool_count, integer_pool_size

    for i from 0 <= i < integer_pool_count:
        o = integer_pool[i]
        mpz_clear(<__mpz_struct *>(<char*>o + mpz_t_offset))
        # Free the object. This assumes that Py_TPFLAGS_HAVE_GC is not
        # set. If it was set another free function would need to be
        # called.
        PyObject_FREE(o)

    integer_pool_size = 0
    integer_pool_count = 0
    sage_free(integer_pool)

# Replace default allocation and deletion with faster custom ones
hook_fast_tp_functions()

# zero and one initialization
initialized = False
cdef set_zero_one_elements():
    global the_integer_ring, initialized
    if initialized: return
    the_integer_ring._zero_element = Integer(0)
    the_integer_ring._one_element = Integer(1)
    init_mpz_globals()
    initialized = True
set_zero_one_elements()

cdef Integer zero = the_integer_ring._zero_element
cdef Integer one = the_integer_ring._one_element

# pool of small integer for fast sign computation
cdef long small_pool_min = -1
cdef long small_pool_max = 1
# we could use the above zero and one here
cdef list small_pool = [Integer(k) for k in range(small_pool_min, small_pool_max+1)]
cdef inline Integer smallInteger(long value):
    cdef Integer z
    if small_pool_min <= value <= small_pool_max:
        return <Integer>small_pool[value - small_pool_min]
    else:
        z = PY_NEW(Integer)
        mpz_set_si(z.value, value)
        return z
