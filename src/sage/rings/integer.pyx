r"""
Elements of the ring `\ZZ` of integers

Sage has highly optimized and extensive functionality for arithmetic with integers
and the ring of integers.

EXAMPLES:

Add 2 integers::

    sage: a = Integer(3); b = Integer(4)
    sage: a + b == 7
    True

Add an integer and a real number::

    sage: a + 4.0                                                                       # needs sage.rings.real_mpfr
    7.00000000000000

Add an integer and a rational number::

    sage: a + Rational(2)/5
    17/5

Add an integer and a complex number::

    sage: # needs sage.rings.real_mpfr
    sage: b = ComplexField().0 + 1.5
    sage: loads((a + b).dumps()) == a + b
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

    sage: a = Integer(3); b = Integer(4)
    sage: a * b == 12
    True
    sage: loads((a * 4.0).dumps()) == a*b
    True
    sage: a * Rational(2)/5
    6/5

::

    sage: [2,3] * 4
    [2, 3, 2, 3, 2, 3, 2, 3]

::

    sage: 'sage' * Integer(3)
    'sagesagesage'

COERCIONS:

Return version of this integer in the multi-precision floating
real field `\RR`::

    sage: n = 9390823
    sage: RR = RealField(200)                                                           # needs sage.rings.real_mpfr
    sage: RR(n)                                                                         # needs sage.rings.real_mpfr
    9.3908230000000000000000000000000000000000000000000000000000e6

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
  :issue:`12116`)

- Vincent Delecroix (2017-05-03): faster integer-rational comparisons

- Vincent Klein (2017-05-11): add __mpz__() to class Integer

- Vincent Klein (2017-05-22): Integer constructor support gmpy2.mpz parameter

- Samuel Lelièvre (2018-08-02): document that divisors are sorted (:issue:`25983`)
"""
# ****************************************************************************
#       Copyright (C) 2004, 2006 William Stein <wstein@gmail.com>
#       Copyright (C) 2006 Gonzalo Tornaria <tornaria@math.utexas.edu>
#       Copyright (C) 2006 Didier Deshommes <dfdeshom@gmail.com>
#       Copyright (C) 2007 David Harvey <dmharvey@math.harvard.edu>
#       Copyright (C) 2007 Martin Albrecht <malb@informatik.uni-bremen.de>
#       Copyright (C) 2007, 2008 Robert Bradshaw <robertwb@math.washington.edu>
#       Copyright (C) 2007 David Roe <roed314@gmail.com>
#       Copyright (C) 2017 Vincent Delecroix <20100.delecroix@gmail.com>
#       Copyright (C) 2018 Samuel Lelièvre <samuel.lelievre@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

# Do not create any Integer, especially non cdef'ed ones, before the hooked
# creation and deletion are setup by the call to hook_fast_tp_functions

cimport cython
from libc.math cimport (ldexp, sqrt as sqrt_double, isnan)
from libc.string cimport memcpy
from libc.limits cimport LONG_MAX

from cysignals.memory cimport check_allocarray, check_malloc, sig_free
from cysignals.signals cimport sig_on, sig_off, sig_check, sig_occurred

import operator

from sage.ext.stdsage cimport PY_NEW
from sage.cpython.python_debug cimport if_Py_TRACE_REFS_then_PyObject_INIT

from sage.libs.gmp.mpz cimport *
from sage.libs.gmp.mpq cimport *
from sage.cpython.string cimport char_to_str, str_to_bytes
from sage.arith.long cimport (integer_check_long,
                              integer_check_long_py, is_small_python_int)

from cpython.list cimport *
from cpython.number cimport *
from cpython.long cimport *
from cpython.object cimport *
from libc.stdint cimport uint64_t
cimport sage.structure.element
from sage.structure.coerce cimport coercion_model
from sage.structure.element cimport (Element, parent)
from sage.structure.parent cimport Parent
from sage.rings.rational cimport Rational
from sage.arith.rational_reconstruction cimport mpq_rational_reconstruction
from sage.libs.gmp.pylong cimport *
from sage.libs.gmp.binop cimport mpq_add_z, mpq_mul_z, mpq_div_zz

import sage.rings.infinity

from sage.structure.coerce cimport is_numpy_type
from sage.structure.element import coerce_binop

from sage.structure.richcmp cimport rich_to_bool_sgn

from sage.rings import integer_ring

cimport gmpy2
gmpy2.import_gmpy2()


try:
    from cypari2.gen import Gen as pari_gen
except ImportError:
    pari_gen = ()


set_integer_from_gen = None
pari_divisors_small = None
n_factor_to_list = None
pari_is_prime_power = None
pari_is_prime = None
objtogen = None
new_gen_from_integer = None


cdef extern from *:
    int unlikely(int) nogil  # Defined by Cython

cdef object numpy_long_interface = {'typestr': '=i4' if sizeof(long) == 4 else '=i8'}
cdef object numpy_int64_interface = {'typestr': '=i8'}
cdef object numpy_object_interface = {'typestr': '|O'}

cdef set_from_Integer(Integer self, Integer other):
    mpz_set(self.value, other.value)


cdef _digits_naive(mpz_t v, l, int offset, Integer base, digits):
    """
    This method fills in digit entries in the list, l, using the most
    basic digit algorithm -- repeat division by base.

    INPUT:

    - ``v`` -- the value whose digits we want to put into the list

    - ``l`` -- the list to file

    - ``offset`` -- offset from the beginning of the list that we want
      to fill at

    - ``base`` -- the base to which we finding digits

    - ``digits`` -- a python sequence type with objects to use for digits
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

cdef _digits_internal(mpz_t v, l, int offset, int power_index, power_list, digits):
    """
    INPUT:

    - ``v`` -- the value whose digits we want to put into the list

    - ``l`` -- the list to file

    - ``offset`` -- offset from the beginning of the list that we want
      to fill at

    - ``power_index`` -- a measure of size to fill and index to
      power_list we're filling ``1 << (power_index+1)`` digits

    - ``power_list`` -- list of powers of the base, precomputed in
      method digits digits - a python sequence type with objects to
      use for digits note that python negative index semantics are
      relied upon

    AUTHORS:

    - Joel B. Mohler (2008-03-13)
    """
    cdef mpz_t mpz_res
    cdef mpz_t mpz_quot
    cdef Integer temp
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


cdef Parent the_integer_ring = integer_ring.ZZ

# The documentation for the ispseudoprime() function in the PARI
# manual states that its result is always prime up to this 2^64.
cdef mpz_t PARI_PSEUDOPRIME_LIMIT
mpz_init(PARI_PSEUDOPRIME_LIMIT)
mpz_ui_pow_ui(PARI_PSEUDOPRIME_LIMIT, 2, 64)


def is_Integer(x):
    """
    Return ``True`` if ``x`` is of the Sage :class:`Integer` type.

    EXAMPLES::

        sage: from sage.rings.integer import is_Integer
        sage: is_Integer(2)
        doctest:warning...
        DeprecationWarning: The function is_Integer is deprecated;
        use 'isinstance(..., Integer)' instead.
        See https://github.com/sagemath/sage/issues/38128 for details.
        True
        sage: is_Integer(2/1)
        False
        sage: is_Integer(int(2))
        False
        sage: is_Integer('5')
        False
    """
    from sage.misc.superseded import deprecation_cython
    deprecation_cython(38128,
                       "The function is_Integer is deprecated; "
                       "use 'isinstance(..., Integer)' instead.")
    return isinstance(x, Integer)


cdef inline Integer as_Integer(x):
    if isinstance(x, Integer):
        return <Integer>x
    else:
        return Integer(x)


cdef class IntegerWrapper(Integer):
    r"""
    Rationale for the :class:`IntegerWrapper` class:

    With :class:`Integer` objects, the allocation/deallocation function slots are
    hijacked with custom functions that stick already allocated
    :class:`Integer` objects (with initialized ``parent`` and ``mpz_t`` fields)
    into a pool on "deallocation" and then pull them out whenever a
    new one is needed. Because :class:`Integers` objects are so common, this is
    actually a significant savings. However, this does cause issues
    with subclassing a Python class directly from :class:`Integer` (but
    that's ok for a Cython class).

    As a workaround, one can instead derive a class from the
    intermediate class :class:`IntegerWrapper`, which sets statically its
    alloc/dealloc methods to the *original* :class:`Integer` alloc/dealloc
    methods, before they are swapped manually for the custom ones.

    The constructor of :class:`IntegerWrapper` further allows for
    specifying an alternative parent to :class:`IntegerRing`.
    """

    def __init__(self, parent=None, x=None, unsigned int base=0):
        """
        We illustrate how to create integers with parents different
        from :class:`IntegerRing`::

            sage: from sage.rings.integer import IntegerWrapper

            sage: n = IntegerWrapper(Primes(), 3) # indirect doctest
            sage: n
            3
            sage: n.parent()
            Set of all prime numbers: 2, 3, 5, 7, ...

        Pickling seems to work now (as of :issue:`10314`)::

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
    The :class:`Integer` class represents arbitrary precision
    integers. It derives from the :class:`Element` class, so
    integers can be used as ring elements anywhere in Sage.

    The constructor of :class:`Integer` interprets strings that begin with ``0o`` as octal numbers,
    strings that begin with ``0x`` as hexadecimal numbers and strings
    that begin with ``0b`` as binary numbers.

    The class :class:`Integer` is implemented in Cython, as a wrapper of the
    GMP ``mpz_t`` integer type.

    EXAMPLES::

        sage: Integer(123)
        123
        sage: Integer("123")
        123

    Sage Integers support :pep:`3127` literals::

        sage: Integer('0x12')
        18
        sage: Integer('-0o12')
        -10
        sage: Integer('+0b101010')
        42

    Conversion from PARI::

        sage: Integer(pari('-10380104371593008048799446356441519384'))                  # needs sage.libs.pari
        -10380104371593008048799446356441519384
        sage: Integer(pari('Pol([-3])'))                                                # needs sage.libs.pari
        -3

    Conversion from gmpy2::

        sage: from gmpy2 import mpz
        sage: Integer(mpz(3))
        3

    .. automethod:: __pow__
    """

    def __cinit__(self):
        global the_integer_ring
        mpz_init(self.value)
        self._parent = the_integer_ring

    def __init__(self, x=None, base=0):
        """
        EXAMPLES::

            sage: a = int(-901824309821093821093812093810928309183091832091)
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
            sage: ZZ(QQbar(sqrt(28-10*sqrt(3)) + sqrt(3)))                              # needs sage.rings.number_field sage.symbolic
            5
            sage: ZZ(AA(32).nth_root(5))                                                # needs sage.rings.number_field
            2
            sage: ZZ(pari('Mod(-3,7)'))                                                 # needs sage.libs.pari
            4
            sage: ZZ('sage')
            Traceback (most recent call last):
            ...
            TypeError: unable to convert 'sage' to an integer
            sage: Integer('zz',36).str(36)
            'zz'
            sage: ZZ('0x3b').str(16)
            '3b'
            sage: ZZ( ZZ(5).digits(3) , 3)
            5
            sage: import numpy                                                          # needs numpy
            sage: ZZ(numpy.int64(7^7))                                                  # needs numpy
            823543
            sage: ZZ(numpy.ubyte(-7))                                                   # needs numpy
            249
            sage: ZZ(True)
            1
            sage: ZZ(False)
            0
            sage: ZZ(1==0)
            0
            sage: ZZ('+10')
            10
            sage: from gmpy2 import mpz
            sage: ZZ(mpz(42))
            42

        ::

            sage: k = GF(2)
            sage: ZZ((k(0),k(1)), 2)
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
            ....:     pass
            sage: class MyFloat(float):
            ....:     pass
            sage: ZZ(MyInt(3))
            3
            sage: ZZ(MyFloat(5))
            5

        ::

            sage: Integer('0')
            0
            sage: Integer('0X2AEEF')
            175855

        Test conversion from PARI (:issue:`11685`)::

            sage: # needs sage.libs.pari
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
            sage: ZZ(GF(3^20,'t')(1))                                                   # needs sage.rings.finite_rings
            1
            sage: ZZ(pari(GF(3^20,'t')(1)))                                             # needs sage.rings.finite_rings
            1
            sage: x = polygen(QQ)
            sage: K.<a> = NumberField(x^2 + 3)                                          # needs sage.rings.number_field
            sage: ZZ(a^2)                                                               # needs sage.rings.number_field
            -3
            sage: ZZ(pari(a)^2)                                                         # needs sage.rings.number_field
            -3
            sage: ZZ(pari("Mod(x, x^3+x+1)"))   # Note error message refers to lifted element
            Traceback (most recent call last):
            ...
            TypeError: Unable to coerce PARI x to an Integer

        Test coercion of `p`-adic with negative valuation::

            sage: ZZ(pari(Qp(11)(11^-7)))                                               # needs sage.libs.pari sage.rings.padics
            Traceback (most recent call last):
            ...
            TypeError: cannot convert p-adic with negative valuation to an integer

        Test converting a list with a very large base::

            sage: a = ZZ(randint(0, 2^128 - 1))
            sage: L = a.digits(2^64)
            sage: a == sum([x * 2^(64*i) for i,x in enumerate(L)])
            True
            sage: a == ZZ(L, base=2^64)
            True

        Test comparisons with numpy types (see :issue:`13386` and :issue:`18076`)::

            sage: # needs numpy
            sage: import numpy
            sage: if int(numpy.version.short_version[0]) > 1:
            ....:     numpy.set_printoptions(legacy="1.25")
            sage: numpy.int8('12') == 12
            True
            sage: 12 == numpy.int8('12')
            True

            sage: float('15') == 15
            True
            sage: 15 == float('15')
            True

        Test underscores as digit separators (PEP 515,
        https://www.python.org/dev/peps/pep-0515/)::

            sage: Integer('1_3')
            13
            sage: Integer(b'1_3')
            13
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
        cdef Py_ssize_t j
        cdef object otmp

        cdef Element lift

        if x is None:
            if mpz_sgn(self.value) != 0:
                mpz_set_si(self.value, 0)

        else:
            # First do all the type-check versions (these are fast to test),
            # except those for which the conversion itself will be slow.

            if isinstance(x, Integer):
                set_from_Integer(self, <Integer>x)

            elif isinstance(x, int):
                mpz_set_pylong(self.value, x)

            elif isinstance(x, float):
                n = int(x)
                if n == x:
                    mpz_set_pylong(self.value, n)
                else:
                    raise TypeError("cannot convert non-integral float to integer")

            elif isinstance(x, pari_gen):
                global set_integer_from_gen
                if set_integer_from_gen is None:
                    from sage.libs.pari.convert_sage import set_integer_from_gen
                set_integer_from_gen(self, x)

            else:

                otmp = getattr(x, "_integer_", None)
                if otmp is not None:
                    set_from_Integer(self, otmp(the_integer_ring))
                    return

                if isinstance(x, Element):
                    try:
                        lift = x.lift()
                        if lift._parent is the_integer_ring:
                            set_from_Integer(self, lift)
                            return
                    except AttributeError:
                        pass

                elif isinstance(x, bytes):
                    if b'_' in x:
                        x = x.replace(b'_', b'')
                    mpz_set_str_python(self.value, x, base)
                    return
                elif isinstance(x, unicode):
                    if '_' in x:
                        x = x.replace('_', '')
                    mpz_set_str_python(self.value, str_to_bytes(x), base)
                    return

                elif isinstance(x, (list, tuple)) and base > 1:
                    b = the_integer_ring(base)
                    if b == 2:  # we use a faster method
                        for j in range(len(x)):
                            otmp = x[j]
                            if isinstance(otmp, int):
                                if (<long> otmp) == 1:
                                    mpz_setbit(self.value, j)
                                if (<long> otmp) != 0:
                                    break
                            else:
                                if not isinstance(otmp, Integer):
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
                    return

                elif is_numpy_type(type(x)):
                    import numpy
                    if isinstance(x, numpy.integer):
                        mpz_set_pylong(self.value, int(x))
                        return

                elif type(x) is gmpy2.mpz:
                    mpz_set(self.value, (<gmpy2.mpz>x).z)
                    return

                raise TypeError("unable to coerce %s to an integer" % type(x))

    def __reduce__(self):
        """
        This is used when pickling integers.

        EXAMPLES::

            sage: n = 5
            sage: t = n.__reduce__(); t
            (<cyfunction make_integer at ...>, ('5',))
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

    def __index__(self):
        """
        Needed so integers can be used as list indices.

        EXAMPLES::

            sage: v = [1,2,3,4,5]
            sage: v[Integer(3)]
            4
            sage: v[Integer(2):Integer(4)]
            [3, 4]

        See :issue:`20750`::

            sage: import re
            sage: p = re.compile('(a)b')
            sage: m = p.match('ab')
            sage: m.group(Integer(0))
            'ab'
            sage: m.group(Integer(1))
            'a'
        """
        return mpz_get_pyintlong(self.value)

    def _im_gens_(self, codomain, im_gens, base_map=None):
        """
        Return the image of ``self`` under the map that sends the generators of
        the parent to im_gens. Since ZZ maps canonically in the category
        of rings, this is just the natural coercion.

        EXAMPLES::

            sage: n = -10
            sage: R = GF(17)
            sage: n._im_gens_(R, [R(1)])
            7
        """
        return codomain.coerce(self)

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
        if isinstance(x, Integer) and isinstance(y, Integer):
            return (<Integer>x)._xor(y)
        return coercion_model.bin_op(x, y, operator.xor)

    def __richcmp__(left, right, int op):
        """
        ``cmp`` for integers.

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
            sage: all((a < b) == (RR(a) < RR(b)) for (a, b) in zip(Ilist, ilist))
            True
            sage: all((a > b) == (RR(a) > RR(b)) for (a, b) in zip(Ilist, ilist))
            True
            sage: all((a == b) == (RR(a) == RR(b)) for (a, b) in zip(Ilist, ilist))
            True
            sage: all((a <= b) == (RR(a) <= RR(b)) for (a, b) in zip(Ilist, ilist))
            True
            sage: all((a >= b) == (RR(a) >= RR(b)) for (a, b) in zip(Ilist, ilist))
            True
            sage: all((a != b) == (RR(a) != RR(b)) for (a, b) in zip(Ilist, ilist))
            True
            sage: all((a < b) == (RR(a) < RR(b)) for (a, b) in zip(Ilist, llist))
            True
            sage: all((a > b) == (RR(a) > RR(b)) for (a, b) in zip(Ilist, llist))
            True
            sage: all((a < b) == (RR(a) < RR(b)) for (a, b) in zip(Ilist, flist))
            True
            sage: all((a > b) == (RR(a) > RR(b)) for (a, b) in zip(Ilist, flist))
            True

        Verify that :issue:`12149` was fixed (and the fix is consistent
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
        cdef mpz_t mpz_tmp

        assert isinstance(left, Integer)

        if isinstance(right, Integer):
            c = mpz_cmp((<Integer>left).value, (<Integer>right).value)
        elif isinstance(right, Rational):
            c = -mpq_cmp_z((<Rational>right).value, (<Integer>left).value)
        elif isinstance(right, int):
            mpz_init(mpz_tmp)
            mpz_set_pylong(mpz_tmp, right)
            c = mpz_cmp((<Integer>left).value, mpz_tmp)
            mpz_clear(mpz_tmp)
        elif isinstance(right, float):
            d = right
            if isnan(d):
                return op == Py_NE
            c = mpz_cmp_d((<Integer>left).value, d)
        else:
            return coercion_model.richcmp(left, right, op)

        return rich_to_bool_sgn(op, c)

    cpdef _richcmp_(left, right, int op):
        r"""
        EXAMPLES::

            sage: from sage.structure.richcmp import op_EQ, op_NE, op_LT, op_GE
            sage: 1._richcmp_(2, op_LT)
            True
            sage: 0._richcmp_(0, op_EQ)
            True
            sage: (-4)._richcmp_(-4, op_NE)
            False
            sage: (-3**10 + 1)._richcmp_(-3**10, op_GE)
            True
        """
        cdef int c
        c = mpz_cmp((<Integer>left).value, (<Integer>right).value)
        return rich_to_bool_sgn(op, c)

    def __copy__(self):
        """
        EXAMPLES::

            sage: n = 2
            sage: copy(n)
            2
            sage: copy(n) is n
            True
        """
        # integers are immutable
        return self

    def __deepcopy__(self, memo):
        """
        EXAMPLES::

            sage: n = 2
            sage: deepcopy(n) is n
            True
        """
        # integers are immutable
        return self

    def list(self):
        """
        Return a list with this integer in it, to be compatible with the
        method for number fields.

        EXAMPLES::

            sage: m = 5
            sage: m.list()
            [5]
        """
        return [self]

    def __dealloc__(self):
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

    def _symbolic_(self, sring):
        """
        Return this integer as symbolic expression.

        EXAMPLES::

            sage: ex = SR(ZZ(7)); ex                                                    # needs sage.symbolic
            7
            sage: parent(ex)                                                            # needs sage.symbolic
            Symbolic Ring
        """
        return sring._force_pyobject(self, force=True)

    def _sympy_(self):
        """
        Convert Sage Integer() to SymPy Integer.

        EXAMPLES::

            sage: n = 5; n._sympy_()                                                    # needs sympy
            5
            sage: n = -5; n._sympy_()                                                   # needs sympy
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
        return '<mn>%s</mn>' % self

    def __mpz__(self):
        """
        Return a gmpy2 integer.

        EXAMPLES::

            sage: a = 5
            sage: a.__mpz__()
            mpz(5)
            sage: from gmpy2 import mpz
            sage: mpz(a)
            mpz(5)

        TESTS::

            sage: a.__mpz__(); raise NotImplementedError("gmpy2 is not installed")
            Traceback (most recent call last):
            ...
            NotImplementedError: gmpy2 is not installed
        """
        return gmpy2.GMPy_MPZ_From_mpz(self.value)

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

            sage: two = Integer(2)
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
            sage: s = big.str()       # long time (2s on sage.math, 2014)
            sage: len(s)              # long time (depends on above defn of s)
            5000001
            sage: s[:10]              # long time (depends on above defn of s)
            '1000000000'
        """
        if base < 2 or base > 36:
            raise ValueError(f"base (={base}) must be between 2 and 36")
        cdef size_t n
        cdef char *s
        n = mpz_sizeinbase(self.value, base) + 2
        s = <char*>check_malloc(n)
        sig_on()
        mpz_get_str(s, base, self.value)
        sig_off()
        k = char_to_str(s)
        sig_free(s)
        return k

    def __format__(self, *args, **kwargs):
        """
        Return a string representation using Python's Format protocol.
        Valid format descriptions are exactly those for Python integers.

        EXAMPLES::

            sage: "{0:#x}; {0:#b}; {0:+05d}".format(ZZ(17))
            '0x11; 0b10001; +0017'
        """
        return int(self).__format__(*args, **kwargs)

    def ordinal_str(self):
        """
        Return a string representation of the ordinal associated to ``self``.

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
        if self < 0:
            raise ValueError("Negative integers are not ordinals.")
        n = self.abs()
        if (n % 100) != 11 and n % 10 == 1:
            th = 'st'
        elif (n % 100) != 12 and n % 10 == 2:
            th = 'nd'
        elif (n % 100) != 13 and n % 10 == 3:
            th = 'rd'
        else:
            th = 'th'
        return n.str() + th

    def hex(self):
        r"""
        Return the hexadecimal digits of ``self`` in lower case.

        .. NOTE::

           '0x' is *not* prepended to the result like is done by the
           corresponding Python function on ``int``. This is for
           efficiency sake--adding and stripping the string wastes
           time; since this function is used for conversions from
           integers to other C-library structures, it is important
           that it be fast.

        EXAMPLES::

            sage: print(Integer(15).hex())
            f
            sage: print(Integer(16).hex())
            10
            sage: print(Integer(16938402384092843092843098243).hex())
            36bb1e3929d1a8fe2802f083
        """
        return self.str(16)

    def oct(self):
        r"""
        Return the digits of ``self`` in base 8.

        .. NOTE::

           '0' (or '0o') is *not* prepended to the result like is done by the
           corresponding Python function on ``int``. This is for
           efficiency sake--adding and stripping the string wastes
           time; since this function is used for conversions from
           integers to other C-library structures, it is important
           that it be fast.

        EXAMPLES::

            sage: print(Integer(800).oct())
            1440
            sage: print(Integer(8).oct())
            10
            sage: print(Integer(-50).oct())
            -62
            sage: print(Integer(-899).oct())
            -1603
            sage: print(Integer(16938402384092843092843098243).oct())
            15535436162247215217705000570203

        Behavior of Sage integers vs. Python integers::

            sage: Integer(10).oct()
            '12'
            sage: oct(int(10))
            '0o12'

            sage: Integer(-23).oct()
            '-27'
            sage: oct(int(-23))
            '-0o27'
        """
        return self.str(8)

    def binary(self):
        """
        Return the binary digits of ``self`` as a string.

        EXAMPLES::

            sage: print(Integer(15).binary())
            1111
            sage: print(Integer(16).binary())
            10000
            sage: print(Integer(16938402384092843092843098243).binary())
            1101101011101100011110001110010010100111010001101010001111111000101000000000101111000010000011
        """
        return self.str(2)

    def bits(self):
        r"""
        Return the bits in ``self`` as a list, least significant first. The
        result satisfies the identity

        ::

            x == sum(b*2^e for e, b in enumerate(x.bits()))

        Negative numbers will have negative "bits". (So, strictly
        speaking, the entries of the returned list are not really
        members of `\ZZ/2\ZZ`.)

        This method just calls :func:`digits` with ``base=2``.

        .. SEEALSO::

            - :meth:`bit_length`, a faster way to compute ``len(x.bits())``
            - :meth:`binary`, which returns a string in perhaps more familiar notation

        EXAMPLES::

            sage: 500.bits()
            [0, 0, 1, 0, 1, 1, 1, 1, 1]
            sage: 11.bits()
            [1, 1, 0, 1]
            sage: (-99).bits()
            [-1, -1, 0, 0, 0, -1, -1]
        """
        return self.digits(base=2)

    def bit_length(self):
        """
        Return the number of bits required to represent this integer.

        Identical to :meth:`int.bit_length`.

        EXAMPLES::

            sage: 500.bit_length()
            9
            sage: 5.bit_length()
            3
            sage: 0.bit_length() == len(0.bits()) == 0.ndigits(base=2)
            True
            sage: 12345.bit_length() == len(12345.binary())
            True
            sage: 1023.bit_length()
            10
            sage: 1024.bit_length()
            11

        TESTS::

            sage: {ZZ(n).bit_length() == int(n).bit_length() for n in range(-9999, 9999)}
            {True}
            sage: n = randrange(-2^99, 2^99)
            sage: ZZ(n).bit_length() == int(n).bit_length()
            True
            sage: n = randrange(-2^999, 2^999)
            sage: ZZ(n).bit_length() == int(n).bit_length()
            True
            sage: n = randrange(-2^9999, 2^9999)
            sage: ZZ(n).bit_length() == int(n).bit_length()
            True
        """
        # mpz_sizeinbase(0, 2) == 1, but int(0).bit_length() == 0
        if mpz_sgn(self.value) == 0:
            return int(0)

        return int(mpz_sizeinbase(self.value, 2))

    def nbits(self):
        r"""
        Alias for :meth:`bit_length`.

        TESTS::

            sage: {ZZ(n).nbits() == ZZ(n).bit_length() for n in range(-9999, 9999)}
            {True}
        """
        return self.bit_length()

    def trailing_zero_bits(self):
        """
        Return the number of trailing zero bits in ``self``, i.e.
        the exponent of the largest power of 2 dividing ``self``.

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

        The returned value is unspecified if ``self`` is a negative number
        and the digits are given.

        INPUT:

        - ``base`` -- integer (default: 10)

        - ``digits`` -- (optional) indexable object as source for
          the digits

        - ``padto`` -- the minimal length of the returned list, sufficient
          number of zeros are added to make the list minimum that length
          (default: 0)

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

            sage: digits = list(range(15))
            sage: l = [ZZ(i).digits(15,digits) for i in range(100)]
            sage: l[16]
            [1, 1]

        This function is comparable to :func:`str` for speed.

        ::

            sage: n=3^100000
            sage: n.digits(base=10)[-1]  # slightly slower than str                     # needs sage.rings.real_interval_field
            1
            sage: n=10^10000
            sage: n.digits(base=10)[-1]  # slightly faster than str                     # needs sage.rings.real_interval_field
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

        if isinstance(base, Integer):
            _base = <Integer>base
        else:
            _base = Integer(base)

        if mpz_cmp_si(_base.value,2) < 0:
            raise ValueError("base must be >= 2")

        if mpz_sgn(self.value) < 0:
            self_abs = -self

        cdef bint do_sig_on
        if mpz_sgn(self.value) == 0:
            l = [zero if digits is None else digits[0]]*padto
        elif mpz_cmp_si(_base.value,2) == 0:
            s = mpz_sizeinbase(self.value, 2)
            if digits:
                o = digits[1]
                z = digits[0]
            else:
                if mpz_sgn(self.value) == 1:
                    o = one
                else:
                    o = -one
                z = zero
            l = [z]*(s if s >= padto else padto)
            for i in range(s):
                # mpz_tstbit seems to return 0 for the high-order bit of
                # negative numbers?!
                if mpz_tstbit(self_abs.value,i):
                    l[i] = o
        else:
            s = mpz_sizeinbase(self.value, 2)
            do_sig_on = (s > 256)
            if do_sig_on:
                sig_on()

            # We use a divide and conquer approach (suggested by the prior
            # author, malb?, of the digits method) here: for base b, compute
            # b^2, b^4, b^8, ... (repeated squaring) until you get larger
            # than your number; then compute (n // b^256, n % b^256)
            # (if b^512 > number) to split the number in half and recurse

            # Pre-computing the exact number of digits up-front is actually
            # faster (especially for large values of self) than trimming off
            # trailing zeros after the fact.  It also seems that it would
            # avoid duplicating the list in memory with a list-slice.
            z = zero if digits is None else digits[0]
            s = self_abs.exact_log(_base)
            l = [z]*(s+1 if s+1 >= padto else padto)

            # set up digits for optimal access once we get inside the worker
            # functions
            if digits is not None:
                # list objects have fastest access in the innermost loop
                if type(digits) is not list:
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

            if do_sig_on:
                sig_off()

        # padding should be taken care of with-in the function
        # all we need to do is return
        return l

    def balanced_digits(self, base=10, positive_shift=True):
        r'''
        Return the list of balanced digits for ``self`` in the given base.

        The balanced base ``b`` uses ``b`` digits centered around zero. Thus
        if ``b`` is odd, there is only one possibility, namely digits
        between ``-b//2`` and ``b//2`` (both included). For instance in base 9,
        one uses digits from -4 to 4. If ``b`` is even, one has to choose
        between digits from ``-b//2`` to ``b//2 - 1`` or ``-b//2 + 1`` to ``b//2``
        (base 10 for instance: either `-5` to `4` or `-4` to `5`), and this is
        defined by the value of ``positive_shift``.

        INPUT:

        - ``base`` -- integer (default: 10); when ``base`` is 2, only the
          nonnegative or the nonpositive integers can be represented by
          ``balanced_digits``. Thus we say base must be greater than 2.

        - ``positive_shift`` -- boolean (default: ``True``); for even bases, the
          representation uses digits from ``-b//2 + 1`` to ``b//2`` if set to
          ``True``, and from ``-b//2`` to ``b//2 - 1`` otherwise. This has no
          effect for odd bases.

        EXAMPLES::

            sage: 8.balanced_digits(3)
            [-1, 0, 1]
            sage: (-15).balanced_digits(5)
            [0, 2, -1]
            sage: 17.balanced_digits(6)
            [-1, 3]
            sage: 17.balanced_digits(6, positive_shift=False)
            [-1, -3, 1]
            sage: (-46).balanced_digits()
            [4, 5, -1]
            sage: (-46).balanced_digits(positive_shift=False)
            [4, -5]
            sage: (-23).balanced_digits(12)
            [1, -2]
            sage: (-23).balanced_digits(12, positive_shift=False)
            [1, -2]
            sage: 0.balanced_digits(7)
            []
            sage: 14.balanced_digits(5.8)
            Traceback (most recent call last):
            ...
            ValueError: base must be an integer
            sage: 14.balanced_digits(2)
            Traceback (most recent call last):
            ...
            ValueError: base must be > 2

        TESTS::

            sage: base = 5; n = 39
            sage: l = n.balanced_digits(base)
            sage: sum(l[i]*base^i for i in range(len(l))) == n
            True
            sage: base = 12; n = -52
            sage: l = n.balanced_digits(base)
            sage: sum(l[i]*base^i for i in range(len(l))) == n
            True
            sage: base = 8; n = 37
            sage: l = n.balanced_digits(base)
            sage: sum(l[i]*base^i for i in range(len(l))) == n
            True
            sage: base = 8; n = 37
            sage: l = n.balanced_digits(base, positive_shift=False)
            sage: sum(l[i]*base^i for i in range(len(l))) == n
            True

        .. SEEALSO::

            :func:`digits <digits>`
        '''
        if not isinstance(base, Integer):
            try:
                base = Integer(base)
            except TypeError:
                raise ValueError('base must be an integer')
        if base <= 2:
            raise ValueError('base must be > 2')

        neg = False
        if self < 0:
            neg = True
            positive_shift = not positive_shift

        if positive_shift or base % 2 == 1:
            m = base//2
        else:
            m = base//2 - 1
        digits = abs(self).digits(base)

        for i in range(len(digits)):
            if digits[i] > m:
                digits[i] = digits[i] - base
                try:
                    digits[i+1] += 1
                except IndexError:
                    if neg:
                        digits.append(-1)
                    else:
                        digits.append(1)
            if neg:
                digits[i] = -digits[i]
        return digits

    def ndigits(self, base=10):
        """
        Return the number of digits of ``self`` expressed in the given base.

        INPUT:

        - ``base`` -- integer (default: 10)

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
            sage: n.ndigits()                                                           # needs sage.rings.real_interval_field
            3000001
            sage: n = 1000**1000000-1
            sage: n.ndigits()                                                           # needs sage.rings.real_interval_field
            3000000
            sage: n = 10**10000000-10**9999990
            sage: n.ndigits()                                                           # needs sage.rings.real_interval_field
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

    cdef void set_from_mpz(Integer self, mpz_t value) noexcept:
        mpz_set(self.value, value)

    def __add__(left, right):
        r"""
        TESTS::

            sage: 1 + 2
            3
            sage: sum(Integer(i) for i in [1..100])
            5050
            sage: 1 + 2/3
            5/3
            sage: 1 + (-2/3)
            1/3
        """
        cdef Integer x
        cdef Rational y
        if type(left) is type(right):
            x = <Integer>PY_NEW(Integer)
            mpz_add(x.value, (<Integer>left).value, (<Integer>right).value)
            return x
        elif type(right) is Rational:
            y = <Rational> Rational.__new__(Rational)
            mpq_add_z(y.value, (<Rational>right).value, (<Integer>left).value)
            return y

        return coercion_model.bin_op(left, right, operator.add)

    cpdef _add_(self, right):
        """
        Integer addition.

        TESTS::

            sage: 32._add_(23)
            55
            sage: a = ZZ.random_element(10^50000)
            sage: b = ZZ.random_element(10^50000)
            sage: a._add_(b) == b._add_(a)
            True
        """
        # self and right are guaranteed to be Integers
        cdef Integer x = <Integer>PY_NEW(Integer)
        mpz_add(x.value, self.value, (<Integer>right).value)
        return x

    cdef _add_long(self, long n):
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

            sage: most_neg_long = int(-sys.maxsize - 1)
            sage: type(most_neg_long), type(-most_neg_long)
            (<class 'int'>, <class 'int'>)
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

    def __sub__(left, right):
        r"""
        TESTS::

            sage: 1 - 2
            -1
            sage: 1 - 2/3
            1/3
            sage: 1 - (-2/3)
            5/3
            sage: (-1) - (-5/4)
            1/4
        """
        cdef Integer x
        cdef Rational y
        if type(left) is type(right):
            x = <Integer>PY_NEW(Integer)
            mpz_sub(x.value, (<Integer>left).value, (<Integer>right).value)
            return x
        elif type(right) is Rational:
            y = <Rational> Rational.__new__(Rational)
            mpz_mul(mpq_numref(y.value), (<Integer>left).value,
                    mpq_denref((<Rational>right).value))
            mpz_sub(mpq_numref(y.value), mpq_numref(y.value),
                    mpq_numref((<Rational>right).value))
            mpz_set(mpq_denref(y.value), mpq_denref((<Rational>right).value))
            return y

        return coercion_model.bin_op(left, right, operator.sub)

    cpdef _sub_(self, right):
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

    cpdef _neg_(self):
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
            sage: b'hi' * 8 == b'hihihihihihihihi'
            True
        """
        if isinstance(s, (list, tuple, str, bytes)):
            if mpz_fits_slong_p(self.value):
                return s * mpz_get_si(self.value)
            else:
                return s * int(self)  # will raise the appropriate exception

    cdef _mul_long(self, long n):
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

    def __mul__(left, right):
        r"""
        TESTS::

            sage: 3 * 2
            6
            sage: 5 * QQ((2,3))
            10/3
            sage: 3 * (-5/6)
            -5/2
            sage: (-2) * (-5/4)
            5/2
        """
        cdef Integer x
        cdef Rational y
        if type(left) is type(right):
            x = <Integer>PY_NEW(Integer)
            mpz_mul(x.value, (<Integer>left).value, (<Integer>right).value)
            return x
        elif type(right) is Rational:
            y = <Rational> Rational.__new__(Rational)
            mpq_mul_z(y.value, (<Rational>right).value, (<Integer>left).value)
            return y

        return coercion_model.bin_op(left, right, operator.mul)

    cpdef _mul_(self, right):
        """
        Integer multiplication.

            sage: 25._mul_(4)
            100
            sage: (5^100)._mul_(2^100)
            10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
            sage: a = ZZ.random_element(10^50000)
            sage: b = ZZ.random_element(10^50000)
            sage: a._mul_(b) == b._mul_(a)
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

    def __truediv__(left, right):
        r"""
        TESTS::

            sage: 3 / 2
            3/2
            sage: 5 / QQ((10,3))
            3/2
            sage: 3 / (-5/6)
            -18/5
            sage: (-2) / (-5/4)
            8/5
            sage: 3 / polygen(ZZ)
            3/x

            sage: 3 / 0
            Traceback (most recent call last):
            ...
            ZeroDivisionError: rational division by zero
            sage: 3 / QQ.zero()
            Traceback (most recent call last):
            ...
            ZeroDivisionError: rational division by zero
            sage: 3 / QQbar.zero()                                                      # needs sage.rings.number_field
            Traceback (most recent call last):
            ...
            ZeroDivisionError: division by zero in algebraic field
        """
        cdef Rational x
        if type(left) is type(right):
            if mpz_sgn((<Integer>right).value) == 0:
                raise ZeroDivisionError("rational division by zero")
            x = <Rational> Rational.__new__(Rational)
            mpq_div_zz(x.value, (<Integer>left).value, (<Integer>right).value)
            return x
        elif type(right) is Rational:
            if mpq_sgn((<Rational>right).value) == 0:
                raise ZeroDivisionError("rational division by zero")
            # left * den(right) / num(right)
            y = <Rational> Rational.__new__(Rational)
            mpq_div_zz(y.value, (<Integer>left).value,
                       mpq_numref((<Rational>right).value))
            mpz_mul(mpq_numref(y.value), mpq_numref(y.value),
                    mpq_denref((<Rational>right).value))
            return y

        return coercion_model.bin_op(left, right, operator.truediv)

    cpdef _div_(self, right):
        r"""
        Compute `\frac{a}{b}`.

        EXAMPLES::

            sage: 3._div_(4)
            3/4
            sage: (-32)._div_(-32)
            1
        """
        if mpz_sgn((<Integer>right).value) == 0:
            raise ZeroDivisionError("rational division by zero")
        x = <Rational> Rational.__new__(Rational)
        mpq_div_zz(x.value, self.value, (<Integer>right).value)
        return x

    cpdef _floordiv_(self, right):
        r"""
        Compute the whole part of `\frac{x}{y}`.

        EXAMPLES::

            sage: a = Integer(321); b = Integer(10)
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
        if not mpz_sgn((<Integer>right).value):
            raise ZeroDivisionError("Integer division by zero")

        cdef Integer z = <Integer>PY_NEW(Integer)
        if mpz_size(self.value) > 1000:
            sig_on()
            mpz_fdiv_q(z.value, self.value, (<Integer>right).value)
            sig_off()
        else:
            mpz_fdiv_q(z.value, self.value, (<Integer>right).value)
        return z

    def __pow__(left, right, modulus):
        r"""
        Return ``(left ^ right) % modulus``.

        EXAMPLES::

            sage: 2^-6
            1/64
            sage: 2^6
            64
            sage: 2^0
            1
            sage: 2^-0
            1
            sage: (-1)^(1/3)                                                            # needs sage.symbolic
            (-1)^(1/3)

        For consistency with Python and MPFR, 0^0 is defined to be 1 in
        Sage::

            sage: 0^0
            1

        See also `<http://www.faqs.org/faqs/sci-math-faq/0to0/>`_ and
        `<https://math.stackexchange.com/questions/11150/zero-to-the-zero-power-is-00-1>`_.

        The base need not be a Sage integer. If it is a Python type, the
        result is a Python type too::

            sage: r = int(2) ^ 10; r; type(r)
            1024
            <... 'int'>
            sage: r = int(3) ^ -3; r; type(r)
            0.037037037037037035
            <... 'float'>
            sage: r = float(2.5) ^ 10; r; type(r)
            9536.7431640625
            <... 'float'>

        We raise 2 to various interesting exponents::

            sage: 2^x                # symbolic x                                       # needs sage.symbolic
            2^x
            sage: 2^1.5              # real number                                      # needs sage.rings.real_mpfr
            2.82842712474619
            sage: 2^float(1.5)       # python float  abs tol 3e-16
            2.8284271247461903
            sage: 2^I                # complex number                                   # needs sage.symbolic
            2^I
            sage: r = 2 ^ int(-3); r; type(r)
            1/8
            <class 'sage.rings.rational.Rational'>
            sage: f = 2^(sin(x)-cos(x)); f                                              # needs sage.symbolic
            2^(-cos(x) + sin(x))
            sage: f(x=3)                                                                # needs sage.symbolic
            2^(-cos(3) + sin(3))

        A symbolic sum::

            sage: # needs sage.symbolic
            sage: x, y, z = var('x,y,z')
            sage: 2^(x + y + z)
            2^(x + y + z)
            sage: 2^(1/2)
            sqrt(2)
            sage: 2^(-1/2)
            1/2*sqrt(2)

        TESTS::

            sage: R.<t> = QQ[]
            sage: 2^t
            Traceback (most recent call last):
            ...
            TypeError: no canonical coercion from Univariate Polynomial
            Ring in t over Rational Field to Rational Field

        Test for :issue:`34143`::

            sage: pow(5,7,13).parent()
            Integer Ring
        """
        if modulus is not None:
            from sage.rings.finite_rings.integer_mod import Mod
            return (Mod(left, modulus) ** right).lift()

        if type(left) is type(right):
            return (<Integer>left)._pow_(right)
        elif isinstance(left, Element):
            return coercion_model.bin_op(left, right, operator.pow)
        # left is a non-Element: do the powering with a Python int
        return left ** int(right)

    cpdef _pow_(self, other):
        """
        Integer powering.

        TESTS::

            sage: 2._pow_(3)
            8
            sage: (-2)._pow_(3)
            -8
            sage: 2._pow_(-3)
            1/8
            sage: (-2)._pow_(-3)
            -1/8
            sage: 2._pow_(4)
            16
            sage: (-2)._pow_(4)
            16
            sage: 2._pow_(-4)
            1/16
            sage: (-2)._pow_(-4)
            1/16
            sage: 0._pow_(3)
            0
            sage: 0._pow_(-3)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: rational division by zero

        The exponent must fit in a ``long`` unless the base is `-1`, `0`, or `1`::

            sage: 2 ^ 100000000000000000000000
            Traceback (most recent call last):
            ...
            OverflowError: exponent must be at most 2147483647           # 32-bit
            OverflowError: exponent must be at most 9223372036854775807  # 64-bit
            sage: 1 ^ 100000000000000000000000
            1
            sage: 1 ^ -100000000000000000000000
            1
            sage: 0 ^ 100000000000000000000000
            0
            sage: 0 ^ -100000000000000000000000
            Traceback (most recent call last):
            ...
            ZeroDivisionError: rational division by zero
            sage: (-1) ^ 100000000000000000000000
            1
            sage: (-1) ^ 100000000000000000000001
            -1
            sage: (-1) ^ -100000000000000000000000
            1
            sage: (-1) ^ -100000000000000000000001
            -1
        """
        cdef mpz_ptr exp = (<Integer>other).value

        if mpz_fits_slong_p(exp):
            return self._pow_long(mpz_get_si(exp))

        # Raising to an exponent which doesn't fit in a long overflows
        # except if the base is -1, 0 or 1.
        cdef long s = LONG_MAX
        if mpz_fits_slong_p(self.value):
            s = mpz_get_si(self.value)

        if s == 0 or s == 1:
            r = self
        elif s == -1:
            if mpz_odd_p(exp):
                r = self
            else:
                r = smallInteger(1)
        else:
            raise OverflowError(f"exponent must be at most {LONG_MAX}")
        if mpz_sgn(exp) >= 0:
            return r
        else:
            return ~r

    cdef _pow_long(self, long n):
        if n == 0:
            return smallInteger(1)
        elif n == 1:
            return self

        cdef Integer x
        cdef Rational q
        if n > 0:
            x = PY_NEW(Integer)
            sig_on()
            mpz_pow_ui(x.value, self.value, n)
            sig_off()
            return x
        else:
            if mpz_sgn(self.value) == 0:
                raise ZeroDivisionError("rational division by zero")
            q = Rational.__new__(Rational)
            sig_on()
            mpz_pow_ui(mpq_denref(q.value), self.value, -n)
            if mpz_sgn(mpq_denref(q.value)) > 0:
                mpz_set_ui(mpq_numref(q.value), 1)
            else:
                # If the denominator was negative, change sign and set
                # numerator to -1
                mpz_set_si(mpq_numref(q.value), -1)
                mpz_abs(mpq_denref(q.value), mpq_denref(q.value))
            sig_off()
            return q

    cpdef _pow_int(self, n):
        """
        Integer powering to an integer exponent.

        TESTS::

            sage: 2._pow_int(int(20))
            1048576
            sage: 1._pow_int(int(2^100))
            1
        """
        return self._pow_(Integer(n))

    def nth_root(self, int n, bint truncate_mode=0):
        r"""
        Return the (possibly truncated) ``n``-th root of ``self``.

        INPUT:

        - ``n`` -- integer `\geq 1` (must fit in the C ``int`` type)

        - ``truncate_mode`` -- boolean, whether to allow truncation if
          ``self`` is not an ``n``-th power

        OUTPUT:

        If ``truncate_mode`` is 0 (default), then returns the exact n'th root
        if ``self`` is an n'th power, or raises a :exc:`ValueError`
        if it is not.

        If ``truncate_mode`` is 1, then if either ``n`` is odd or ``self`` is
        positive, returns a pair ``(root, exact_flag)`` where ``root`` is the
        truncated ``n``-th root (rounded towards zero) and ``exact_flag`` is a
        boolean indicating whether the root extraction was exact;
        otherwise raises a :exc:`ValueError`.

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
            raise ValueError("n (=%s) must be positive" % n)
        if (mpz_sgn(self.value) < 0) and not (n & 1):
            raise ValueError("cannot take even root of negative number")
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
                raise ValueError("%s is not a %s power" % (self,
                                                           integer_ring.ZZ(n).ordinal_str()))

    cpdef size_t _exact_log_log2_iter(self, Integer m) noexcept:
        r"""
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
        cdef mpz_t accum
        cdef mpz_t temp_exp

        if mpz_cmp_si(m.value,4) < 0:
            raise ValueError("This is undefined or possibly non-convergent with this algorithm.")

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

    cpdef size_t _exact_log_mpfi_log(self, m) noexcept:
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

            sage: Integer(125)._exact_log_mpfi_log(3)                                   # needs sage.rings.real_interval_field
            4
            sage: Integer(5^150)._exact_log_mpfi_log(5)                                 # needs sage.rings.real_interval_field
            150
        """
        cdef int i
        cdef list pow_2_things
        cdef int pow_2
        cdef size_t upper,lower,middle

        from sage.rings.real_mpfi import RIF as R

        rif_self = R(self)

        sig_on()
        rif_m = R(m)
        rif_log = rif_self.log()/rif_m.log()
        # upper is *greater* than the answer
        try:
            upper = rif_log.upper().ceiling()
        except Exception:
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
            pow_2_things = [rif_m]*pow_2
            for i in range(pow_2):
                pow_2_things[i] = pow_2_things[i-1]**2
            for i from pow_2>i>=0:
                middle = lower + int(2)**i
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
            raise ValueError("The input for exact_log is too large and support is not implemented.")

        return lower

    def exact_log(self, m):
        r"""
        Return the largest integer `k` such that `m^k \leq \text{self}`,
        i.e., the floor of `\log_m(\text{self})`.

        This is guaranteed to return the correct answer even when the usual
        log function doesn't have sufficient precision.

        INPUT:

        - ``m`` -- integer `\geq 2`

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
            sage: Integer(1780^1700-1).exact_log(1780)                                  # needs sage.rings.real_interval_field
            1699

            sage: # The following are very very fast.
            sage: # Note that for base m a perfect power of 2, we get the exact log by counting bits.
            sage: n = 2983579823750185701375109835; m = 32
            sage: n.exact_log(m)
            18
            sage: # The next is a favorite of mine.  The log2 approximate is exact and immediately provable.
            sage: n = 90153710570912709517902579010793251709257901270941709247901209742124
            sage: m = 213509721309572
            sage: n.exact_log(m)
            4

        ::

            sage: # needs sage.rings.real_mpfr
            sage: x = 3^100000
            sage: RR(log(RR(x), 3))
            100000.000000000
            sage: RR(log(RR(x + 100000), 3))
            100000.000000000

        ::

            sage: # needs sage.rings.real_mpfr
            sage: x.exact_log(3)
            100000
            sage: (x + 1).exact_log(3)
            100000
            sage: (x - 1).exact_log(3)
            99999

        ::

            sage: # needs sage.rings.real_mpfr
            sage: x.exact_log(2.5)
            Traceback (most recent call last):
            ...
            TypeError: Attempt to coerce non-integral RealNumber to Integer
        """
        cdef Integer _m
        cdef Integer result
        cdef size_t n_log2
        cdef size_t m_log2
        cdef size_t guess  # this will contain the final answer
        cdef bint guess_filled = 0  # this variable is only used in one branch below
        if isinstance(m, Integer):
            _m = <Integer>m
        else:
            _m = <Integer>Integer(m)

        self_sgn = mpz_sgn(self.value)
        if self_sgn == 0:
            return -sage.rings.infinity.infinity
        if self_sgn < 0 or mpz_sgn(_m.value) <= 0:
            raise ValueError("self must be nonnegative and m must be positive")
        if mpz_cmp_si(_m.value,2) < 0:
            raise ValueError("m must be at least 2")

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

            if n_log2 / m_log2 > 8000:
                # If we have a very large number of digits, it can be a nice
                # shortcut to test the guess using interval arithmetic.
                # (suggested by David Harvey and Carl Witty)
                # "for randomly distributed integers, the chance of this
                # interval-based comparison failing is absurdly low"
                from sage.rings.real_mpfi import RIF
                approx_compare = RIF(m)**guess
                if self > approx_compare:
                    guess_filled = 1
                elif self < approx_compare:
                    guess_filled = 1
                    guess = guess - 1
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
        Return symbolic log by default, unless the logarithm is exact (for
        an integer argument). When ``prec`` is given, the :class:`RealField`
        approximation to that bit precision is used.

        This function is provided primarily so that Sage integers may be
        treated in the same manner as real numbers when convenient. Direct
        use of :meth:`exact_log` is probably best for arithmetic log computation.

        INPUT:

        - ``m`` -- (default: natural) log base e

        - ``prec`` -- integer (default: ``None``); if ``None``, returns
          symbolic, else to given bits of precision as in :class:`RealField`

        EXAMPLES::

            sage: Integer(124).log(5)                                                   # needs sage.symbolic
            log(124)/log(5)
            sage: Integer(124).log(5, 100)                                              # needs sage.rings.real_mpfr
            2.9950093311241087454822446806
            sage: Integer(125).log(5)
            3
            sage: Integer(125).log(5, prec=53)                                          # needs sage.rings.real_mpfr
            3.00000000000000
            sage: log(Integer(125))                                                     # needs sage.symbolic
            3*log(5)

        For extremely large numbers, this works::

            sage: x = 3^100000
            sage: log(x, 3)                                                             # needs sage.rings.real_interval_field
            100000

        Also ``log(x)``, giving a symbolic output,
        works in a reasonable amount of time for this ``x``::

            sage: x = 3^100000
            sage: log(x)                                                                # needs sage.symbolic
            log(1334971414230...5522000001)

        But approximations are probably more useful in this
        case, and work to as high a precision as we desire::

            sage: x.log(3, 53)  # default precision for RealField                       # needs sage.rings.real_mpfr
            100000.000000000
            sage: (x + 1).log(3, 53)                                                    # needs sage.rings.real_mpfr
            100000.000000000
            sage: (x + 1).log(3, 1000)                                                  # needs sage.rings.real_mpfr
            100000.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

        We can use non-integer bases, with default e::

            sage: x.log(2.5, prec=53)                                                   # needs sage.rings.real_mpfr
            119897.784671579

        We also get logarithms of negative integers, via the
        symbolic ring, using the branch from `-\pi` to `\pi`::

            sage: log(-1)                                                               # needs sage.symbolic
            I*pi

        The logarithm of zero is done likewise::

            sage: log(0)                                                                # needs sage.symbolic
            -Infinity

        Some rational bases yield integer logarithms (:issue:`21517`)::

            sage: ZZ(8).log(1/2)
            -3

        Check that Python ints are accepted (:issue:`21518`)::

            sage: ZZ(8).log(int(2))
            3

        TESTS::

            sage: (-2).log(3)                                                           # needs sage.symbolic
            (I*pi + log(2))/log(3)
        """
        cdef int self_sgn
        if m is not None and m <= 0:
            raise ValueError("log base must be positive")
        self_sgn = mpz_sgn(self.value)
        if self_sgn < 0 and prec is None:
            from sage.symbolic.ring import SR
            return SR(self).log(m)
        if prec:
            if self_sgn >= 0:
                from sage.rings.real_mpfr import RealField
                return RealField(prec)(self).log(m)
            else:
                from sage.rings.complex_mpfr import ComplexField
                return ComplexField(prec)(self).log(m)

        if m is None:
            from sage.functions.log import function_log
            return function_log(self,dont_call_method_on_arg=True)
        try:
            m = Integer(m)
        except (ValueError, TypeError):
            pass

        if isinstance(m, Integer):
            elog = self.exact_log(m)
            if elog == -sage.rings.infinity.infinity or m**elog == self:
                return elog

        if isinstance(m, Rational) and m.numer() == 1:
            elog = -self.exact_log(m.denom())
            if m**elog == self:
                return elog

        from sage.functions.log import function_log
        return function_log(self, dont_call_method_on_arg=True)/\
            function_log(m, dont_call_method_on_arg=True)

    def exp(self, prec=None):
        r"""
        Return the exponential function of ``self`` as a real number.

        This function is provided only so that Sage integers may be treated
        in the same manner as real numbers when convenient.

        INPUT:

        - ``prec`` -- integer (default: ``None``); if ``None``, returns
          symbolic, else to given bits of precision as in :class:`RealField`

        EXAMPLES::

            sage: Integer(8).exp()                                                      # needs sage.symbolic
            e^8
            sage: Integer(8).exp(prec=100)                                              # needs sage.symbolic
            2980.9579870417282747435920995
            sage: exp(Integer(8))                                                       # needs sage.symbolic
            e^8

        For even fairly large numbers, this may not be useful.

        ::

            sage: y = Integer(145^145)
            sage: y.exp()                                                               # needs sage.symbolic
            e^25024207011349079210459585279553675697932183658421565260323592409432707306554163224876110094014450895759296242775250476115682350821522931225499163750010280453185147546962559031653355159703678703793369785727108337766011928747055351280379806937944746847277089168867282654496776717056860661614337004721164703369140625
            sage: y.exp(prec=53)  # default RealField precision                         # needs sage.symbolic
            +infinity
        """
        from sage.functions.all import exp
        res = exp(self, dont_call_method_on_arg=True)
        if prec:
            return res.n(prec=prec)
        return res

    def prime_to_m_part(self, m):
        """
        Return the prime-to-`m` part of ``self``, i.e., the largest divisor of
        ``self`` that is coprime to `m`.

        INPUT:

        - ``m`` -- integer

        OUTPUT: integer

        EXAMPLES::

            sage: 43434.prime_to_m_part(20)
            21717
            sage: 2048.prime_to_m_part(2)
            1
            sage: 2048.prime_to_m_part(3)
            2048

            sage: 0.prime_to_m_part(2)
            Traceback (most recent call last):
            ...
            ArithmeticError: self must be nonzero
        """
        cdef Integer mm = Integer(m)

        if not self:
            raise ArithmeticError("self must be nonzero")
        if not mm:
            return one

        cdef Integer n = Integer(self)  # need a copy as it is modified below

        sig_on()
        while mpz_cmp_ui(mm.value, 1):
            mpz_gcd(mm.value, n.value, mm.value)
            mpz_divexact(n.value, n.value, mm.value)
        sig_off()

        return n

    def prime_divisors(self, *args, **kwds):
        """
        Return the prime divisors of this integer, sorted in increasing order.

        If this integer is negative, we do *not* include `-1` among
        its prime divisors, since `-1` is not a prime number.

        INPUT:

        - ``limit`` -- (integer, optional keyword argument)
          Return only prime divisors up to this bound, and the factorization
          is done by checking primes up to ``limit`` using trial division.

        Any additional arguments are passed on to the :meth:`factor` method.

        EXAMPLES::

            sage: a = 1; a.prime_divisors()
            []
            sage: a = 100; a.prime_divisors()
            [2, 5]
            sage: a = -100; a.prime_divisors()
            [2, 5]
            sage: a = 2004; a.prime_divisors()
            [2, 3, 167]

        Setting the optional ``limit`` argument works as expected::

            sage: a = 10^100 + 1
            sage: a.prime_divisors()                                                    # needs sage.libs.pari
            [73, 137, 401, 1201, 1601, 1676321, 5964848081,
             129694419029057750551385771184564274499075700947656757821537291527196801]
            sage: a.prime_divisors(limit=10^3)
            [73, 137, 401]
            sage: a.prime_divisors(limit=10^7)
            [73, 137, 401, 1201, 1601, 1676321]
        """
        res = [r[0] for r in self.factor(*args, **kwds)]
        limit = kwds.get('limit')
        if limit is not None:
            res = [r for r in res if r <= limit]
        return res

    prime_factors = prime_divisors

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def divisors(self, method=None):
        """
        Return the list of all positive integer divisors of this integer,
        sorted in increasing order.

        EXAMPLES:

        ::

            sage: (-3).divisors()
            [1, 3]
            sage: 6.divisors()
            [1, 2, 3, 6]
            sage: 28.divisors()
            [1, 2, 4, 7, 14, 28]
            sage: (2^5).divisors()
            [1, 2, 4, 8, 16, 32]
            sage: 100.divisors()
            [1, 2, 4, 5, 10, 20, 25, 50, 100]
            sage: 1.divisors()
            [1]
            sage: 0.divisors()
            Traceback (most recent call last):
            ...
            ValueError: n must be nonzero
            sage: (2^3 * 3^2 * 17).divisors()
            [1, 2, 3, 4, 6, 8, 9, 12, 17, 18, 24, 34, 36, 51, 68, 72,
            102, 136, 153, 204, 306, 408, 612, 1224]
            sage: a = odd_part(factorial(31))
            sage: v = a.divisors()                                                      # needs sage.libs.pari
            sage: len(v)                                                                # needs sage.libs.pari
            172800
            sage: prod(e + 1 for p, e in factor(a))                                     # needs sage.libs.pari
            172800
            sage: all(t.divides(a) for t in v)                                          # needs sage.libs.pari
            True

        ::

            sage: n = 2^551 - 1
            sage: L = n.divisors()                                                      # needs sage.libs.pari
            sage: len(L)                                                                # needs sage.libs.pari
            256
            sage: L[-1] == n                                                            # needs sage.libs.pari
            True

        TESTS:

        Overflow::

            sage: prod(primes_first_n(64)).divisors()                                   # needs sage.libs.pari
            Traceback (most recent call last):
            ...
            OverflowError: value too large
            sage: prod(primes_first_n(58)).divisors()                                   # needs sage.libs.pari
            Traceback (most recent call last):
            ...
            OverflowError: value too large                                 # 32-bit
            MemoryError: failed to allocate 288230376151711744 * 24 bytes  # 64-bit

        Check for memory leaks and ability to interrupt
        (the ``divisors`` call below allocates about 800 MB every time,
        so a memory leak will not go unnoticed)::

            sage: n = prod(primes_first_n(25))                                          # needs sage.libs.pari
            sage: for i in range(20):           # long time                             # needs sage.libs.pari
            ....:     try:
            ....:         alarm(RDF.random_element(1e-3, 0.5))
            ....:         _ = n.divisors()
            ....:         cancel_alarm()  # we never get here
            ....:     except AlarmInterrupt:
            ....:         pass

        Test a strange method::

            sage: 100.divisors(method='hey')
            Traceback (most recent call last):
            ...
            ValueError: method must be 'pari' or 'sage'


        .. NOTE::

           If one first computes all the divisors and then sorts it,
           the sorting step can easily dominate the runtime. Note,
           however, that (nonnegative) multiplication on the left
           preserves relative order. One can leverage this fact to
           keep the list in order as one computes it using a process
           similar to that of the merge sort algorithm.
        """
        if mpz_cmp_ui(self.value, 0) == 0:
            raise ValueError("n must be nonzero")

        if (method is None or method == 'pari') and mpz_fits_slong_p(self.value):
            global pari_divisors_small
            if pari_divisors_small is None:
                try:
                    from sage.libs.pari.convert_sage import pari_divisors_small
                except ImportError:
                    if method == 'pari':
                        raise ImportError("method `pari` requested, but cypari2 not present")
            if pari_divisors_small is not None:
                if mpz_sgn(self.value) > 0:
                    return pari_divisors_small(self)
                else:
                    return pari_divisors_small(-self)
        elif method is not None and method != 'sage':
            raise ValueError("method must be 'pari' or 'sage'")

        cdef list all, prev, sorted
        cdef Py_ssize_t tip, top
        cdef Py_ssize_t i, e, ee
        cdef Integer apn, p, pn, z, all_tip

        f = self.factor()

        # All of the declarations below are for optimizing the unsigned long-sized
        # case.  Operations are performed in C as far as possible without
        # overflow before moving to Python objects.
        cdef unsigned long p_c, pn_c, apn_c
        cdef Py_ssize_t all_len, sorted_len, prev_len
        cdef unsigned long* ptr
        cdef unsigned long* swap_tmp
        cdef unsigned long* all_c
        cdef unsigned long* sorted_c
        cdef unsigned long* prev_c

        # These are used to keep track of whether or not we are able to
        # perform the operations in machine words. A factor of 0.999
        # safety margin is added to cover any floating-point rounding
        # issues.
        cdef bint fits_c = True
        cdef double cur_max = 1
        cdef double fits_max = 0.999 * 2.0 ** (8*sizeof(unsigned long))

        cdef Py_ssize_t divisor_count = 1
        with cython.overflowcheck(True):
            for p, e in f:
                # Using *= does not work, see
                # https://github.com/cython/cython/issues/1381
                divisor_count = divisor_count * (1 + e)

        ptr = <unsigned long*>check_allocarray(divisor_count, 3 * sizeof(unsigned long))
        all_c = ptr
        sorted_c = ptr + divisor_count
        prev_c = sorted_c + divisor_count

        try:
            sorted_c[0] = 1
            sorted_len = 1

            for p, e in f:
                cur_max *= (<double>p)**e
                if fits_c and cur_max > fits_max:
                    sorted = []
                    for i in range(sorted_len):
                        z = <Integer>PY_NEW(Integer)
                        mpz_set_ui(z.value, sorted_c[i])
                        sorted.append(z)
                    fits_c = False
                    sig_free(ptr)
                    ptr = NULL

                # The two cases below are essentially the same algorithm, one
                # operating on Integers in Python lists, the other on unsigned long's.
                if fits_c:
                    sig_on()

                    pn_c = p_c = p

                    swap_tmp = sorted_c
                    sorted_c = prev_c
                    prev_c = swap_tmp
                    prev_len = sorted_len
                    sorted_len = 0

                    tip = 0
                    prev_c[prev_len] = prev_c[prev_len-1] * pn_c
                    for i in range(prev_len):
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
                        for i in range(prev_len):
                            apn_c = prev_c[i] * pn_c
                            while all_c[tip] < apn_c:
                                sorted_c[sorted_len] = all_c[tip]
                                sorted_len += 1
                                tip += 1
                            sorted_c[sorted_len] = apn_c
                            sorted_len += 1

                    sig_off()

                else:
                    # fits_c is False: use mpz integers
                    prev = sorted
                    pn = <Integer>PY_NEW(Integer)
                    mpz_set_ui(pn.value, 1)
                    for ee in range(e):
                        all = sorted
                        sorted = []
                        tip = 0
                        top = len(all)
                        mpz_mul(pn.value, pn.value, p.value)  # pn *= p
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
                for i in range(sorted_len):
                    z = <Integer>PY_NEW(Integer)
                    mpz_set_ui(z.value, sorted_c[i])
                    sorted.append(z)
        finally:
            sig_free(ptr)

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
        Compute ``|self|``.

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

    def euclidean_degree(self):
        r"""
        Return the degree of this element as an element of a Euclidean domain.

        If this is an element in the ring of integers, this is simply its
        absolute value.

        EXAMPLES::

            sage: ZZ(1).euclidean_degree()
            1
        """
        from sage.rings.integer_ring import ZZ
        if self.parent() is ZZ:
            return abs(self)
        raise NotImplementedError

    def sign(self):
        """
        Return the sign of this integer, which is `-1`, `0`, or `1`
        depending on whether this number is negative, zero, or positive
        respectively.

        OUTPUT: integer

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
        Return x modulo y.

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
            sage: 5 % int(-7)
            -2
            sage: int(5) % -7
            -2
            sage: int(5) % int(-7)
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

        This example caused trouble in :issue:`6083`::

            sage: a = next_prime(2**31)                                                 # needs sage.libs.pari
            sage: b = Integers(a)(100)                                                  # needs sage.libs.pari
            sage: a % b                                                                 # needs sage.libs.pari
            Traceback (most recent call last):
            ...
            ArithmeticError: reduction modulo 100 not defined
        """
        cdef Integer z

        # First case: Integer % Integer
        if type(x) is type(y):
            if not mpz_sgn((<Integer>y).value):
                raise ZeroDivisionError("Integer modulo by zero")
            z = <Integer>PY_NEW(Integer)
            if mpz_size((<Integer>x).value) > 100000:
                sig_on()
                mpz_fdiv_r(z.value, (<Integer>x).value, (<Integer>y).value)
                sig_off()
            else:
                mpz_fdiv_r(z.value, (<Integer>x).value, (<Integer>y).value)
            return z

        # Next: Integer % C long
        cdef long yy = 0
        cdef int err = 0
        if not isinstance(y, Element):
            # x must be an Integer in this case
            if not integer_check_long(y, &yy, &err):
                # y cannot be converted to an integer
                return NotImplemented
            if err:
                # y is some kind of integer,
                # but too large for a C long
                return x % Integer(y)

            if yy == 0:
                raise ZeroDivisionError("Integer modulo by zero")
            z = <Integer>PY_NEW(Integer)
            if yy > 0:
                mpz_fdiv_r_ui(z.value, (<Integer>x).value, yy)
            else:
                mpz_cdiv_r_ui(z.value, (<Integer>x).value, -<unsigned long>yy)
            return z

        # Use the coercion model
        return coercion_model.bin_op(x, y, operator.mod)

    def quo_rem(Integer self, other):
        """
        Return the quotient and the remainder of ``self`` divided by ``other``.
        Note that the remainder returned is always either zero or of the
        same sign as ``other``.

        INPUT:

        - ``other`` -- the divisor

        OUTPUT:

        - ``q`` -- the quotient of ``self/other``

        - ``r`` -- the remainder of ``self/other``

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
        will always be zero (:issue:`7965`)::

            sage: 5.quo_rem(QQ(2))
            (5/2, 0)
            sage: 5.quo_rem(2/3)
            (15/2, 0)

        Check that :issue:`29009` is fixed:

            sage: divmod(1, sys.maxsize+1r)  # should not raise OverflowError: Python int too large to convert to C long
            (0, 1)

            sage: # needs mpmath
            sage: import mpmath
            sage: mpmath.mp.prec = 1000
            sage: root = mpmath.findroot(lambda x: x^2 - 3, 2)
            sage: len(str(root))
            301
        """
        cdef Integer q = PY_NEW(Integer)
        cdef Integer r = PY_NEW(Integer)
        cdef long d, res

        if is_small_python_int(other):
            d = PyLong_AsLong(other)
            if d > 0:
                mpz_fdiv_qr_ui(q.value, r.value, self.value, d)
            elif d == 0:
                raise ZeroDivisionError("Integer division by zero")
            else:
                res = mpz_fdiv_qr_ui(q.value, r.value, self.value, -d)
                mpz_neg(q.value, q.value)
                if res:
                    mpz_sub_ui(q.value, q.value, 1)
                    mpz_sub_ui(r.value, r.value, -d)

        elif type(other) is Integer:
            if mpz_sgn((<Integer>other).value) == 0:
                raise ZeroDivisionError("Integer division by zero")
            if mpz_size(self.value) > 100000:
                sig_on()
                mpz_fdiv_qr(q.value, r.value, self.value, (<Integer>other).value)
                sig_off()
            else:
                mpz_fdiv_qr(q.value, r.value, self.value, (<Integer>other).value)

        else:
            left, right = coercion_model.canonical_coercion(self, other)
            return left.quo_rem(right)

        return q, r

    def powermod(self, exp, mod):
        r"""
        Compute ``self**exp`` modulo ``mod``.

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
        _exp = Integer(exp)
        _mod = Integer(mod)
        if mpz_cmp_si(_mod.value,0) == 0:
            raise ZeroDivisionError("cannot raise to a power modulo 0")

        x = PY_NEW(Integer)

        sig_on()
        mpz_powm(x.value, self.value, _exp.value, _mod.value)
        sig_off()

        return x

    def rational_reconstruction(self, Integer m):
        r"""
        Return the rational reconstruction of this integer modulo `m`, i.e.,
        the unique (if it exists) rational number that reduces to ``self``
        modulo m and whose numerator and denominator is bounded by
        `\sqrt{m/2}`.

        INPUT:

        - ``self`` -- integer

        - ``m`` -- integer

        OUTPUT: a :class:`Rational`

        EXAMPLES::

            sage: (3/7)%100
            29
            sage: (29).rational_reconstruction(100)
            3/7

        TESTS:

        Check that :issue:`9345` is fixed::

            sage: 0.rational_reconstruction(0)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: rational reconstruction with zero modulus
            sage: ZZ.random_element(-10^6, 10^6).rational_reconstruction(0)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: rational reconstruction with zero modulus
        """
        cdef Integer a
        cdef Rational x = <Rational>Rational.__new__(Rational)
        try:
            mpq_rational_reconstruction(x.value, self.value, m.value)
        except ValueError:
            a = self % m
            raise ArithmeticError("rational reconstruction of %s (mod %s) does not exist" % (a, m))
        return x

    def __int__(self):
        """
        Return the Python int corresponding to this Sage integer.

        EXAMPLES::

            sage: n = 920938
            sage: int(n)
            920938
            sage: int(-n)
            -920938
            sage: type(n.__int__())
            <... 'int'>
            sage: n = 99028390823409823904823098490238409823490820938
            sage: int(n)
            99028390823409823904823098490238409823490820938
            sage: int(-n)
            -99028390823409823904823098490238409823490820938
            sage: type(n.__int__())
            <class 'int'>
            sage: int(-1), int(0), int(1)
            (-1, 0, 1)
        """
        return mpz_get_pyintlong(self.value)

    def __float__(self):
        """
        Return double precision floating point representation of this
        integer.

        EXAMPLES::

            sage: n = Integer(17); float(n)
            17.0
            sage: n = Integer(902834098234908209348209834092834098); float(n)
            9.028340982349083e+35
            sage: n = Integer(-57); float(n)
            -57.0
            sage: n.__float__()
            -57.0
            sage: type(n.__float__())
            <... 'float'>
        """
        return mpz_get_d_nearest(self.value)

    def _rpy_(self):
        """
        Return int(self) so that rpy can convert ``self`` into an object it
        knows how to work with.

        EXAMPLES::

            sage: n = 100
            sage: n._rpy_()
            100
            sage: type(n._rpy_())
            <... 'int'>
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
            sage: n = -920390823904823094890238490238484
            sage: n.__hash__()    # random
            -43547310504077801
            sage: n.__hash__() == hash(int(n))
            True

        TESTS::

            sage: hash(-1), hash(0), hash(1)
            (-2, 0, 1)
            sage: n = 2^31 + 2^63 + 2^95 + 2^127 + 2^128*(2^32-2)
            sage: hash(n) == hash(int(n))
            True
            sage: hash(n-1) == hash(int(n-1))
            True
            sage: hash(-n) == hash(int(-n))
            True
            sage: hash(1-n) == hash(int(1-n))
            True
            sage: n = 2^63 + 2^127 + 2^191 + 2^255 + 2^256*(2^64-2)
            sage: hash(n) == hash(int(n))
            True
            sage: hash(n-1) == hash(int(n-1))
            True
            sage: hash(-n) == hash(int(-n))
            True
            sage: hash(1-n) == hash(int(1-n))
            True

        These tests come from :issue:`4957`::

            sage: n = 2^31 + 2^13
            sage: hash(n)             # random
            2147491840
            sage: hash(n) == hash(int(n))
            True
            sage: n = 2^63 + 2^13
            sage: hash(n)             # random
            8196
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
        Return smallest prime divisor of ``self`` up to bound, beginning
        checking at ``start``, or ``abs(self)`` if no such divisor is found.

        INPUT:

        - ``bound`` -- positive integer that fits in a C ``signed long``
        - ``start`` -- positive integer that fits in a C ``signed long``

        OUTPUT: a positive integer

        EXAMPLES::

            sage: # needs sage.libs.pari
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

            sage: # needs sage.libs.pari
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
            raise ValueError("bound must be positive")
        if mpz_sgn(self.value) == 0:
            raise ValueError("self must be nonzero")
        cdef unsigned long n, m=7, i=1, limit
        cdef unsigned long dif[8]
        if start > 7:
            # We need to find i.
            m = start % 30
            if 0 <= m <= 1:
                i = 0
                m = start + (1-m)
            elif 1 < m <= 7:
                i = 1
                m = start + (7-m)
            elif 7 < m <= 11:
                i = 2
                m = start + (11-m)
            elif 11 < m <= 13:
                i = 3
                m = start + (13-m)
            elif 13 < m <= 17:
                i = 4
                m = start + (17-m)
            elif 17 < m <= 19:
                i = 5
                m = start + (19-m)
            elif 19 < m <= 23:
                i = 6
                m = start + (23-m)
            elif 23 < m <= 29:
                i = 7
                m = start + (29-m)
        dif[0] = 6
        dif[1] = 4
        dif[2] = 2
        dif[3] = 4
        dif[4] = 2
        dif[5] = 4
        dif[6] = 6
        dif[7] = 2
        cdef Integer x = PY_NEW(Integer)
        if mpz_fits_ulong_p(self.value):
            n = mpz_get_ui(self.value)   # ignores the sign automatically
            if n == 1:
                return one
            if start <= 2 and n % 2 == 0:
                mpz_set_ui(x.value,2)
                return x
            if start <= 3 and n % 3 == 0:
                mpz_set_ui(x.value,3)
                return x
            if start <= 5 and n % 5 == 0:
                mpz_set_ui(x.value,5)
                return x
            limit = <unsigned long> sqrt_double(<double> n)
            if bound < limit:
                limit = bound
            # Algorithm: only trial divide by numbers that
            # are congruent to 1,7,11,13,17,19,23,29 mod 30=2*3*5.
            while m <= limit:
                if n % m == 0:
                    mpz_set_ui(x.value, m)
                    return x
                m += dif[i % 8]
                i += 1
            mpz_abs(x.value, self.value)
            return x
        else:
            # self is big -- it doesn't fit in unsigned long.
            if start <= 2 and mpz_even_p(self.value):
                mpz_set_ui(x.value, 2)
                return x
            if start <= 3 and mpz_divisible_ui_p(self.value, 3):
                mpz_set_ui(x.value, 3)
                return x
            if start <= 5 and mpz_divisible_ui_p(self.value, 5):
                mpz_set_ui(x.value, 5)
                return x

            # x.value = floor(sqrt(self.value))
            sig_on()
            mpz_abs(x.value, self.value)
            mpz_sqrt(x.value, x.value)
            if mpz_cmp_si(x.value, bound) < 0:
                limit = mpz_get_ui(x.value)
            else:
                limit = bound
            while m <= limit:
                if mpz_divisible_ui_p(self.value, m):
                    mpz_set_ui(x.value, m)
                    sig_off()
                    return x
                m += dif[i % 8]
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

        - ``algorithm`` -- string; one of

          - ``'pari'`` -- (default) use the PARI library

          - ``'flint'`` -- use the FLINT library

          - ``'kash'`` -- use the KASH computer algebra system (requires
            kash)

          - ``'magma'`` -- use the MAGMA computer algebra system (requires
            an installation of MAGMA)

          - ``'qsieve'`` -- use Bill Hart's quadratic sieve code;
            WARNING: this may not work as expected, see qsieve? for
            more information

          - ``'ecm'`` -- use ECM-GMP, an implementation of Hendrik
            Lenstra's elliptic curve method

        - ``proof`` -- boolean (default: ``True``); whether or not to prove
          primality of each factor (only applicable for ``'pari'`` and ``'ecm'``)

        - ``limit`` -- integer or ``None`` (default: ``None``); if limit is
          given it must fit in a ``signed int``, and the factorization is done
          using trial division and primes up to limit

        OUTPUT: a Factorization object containing the prime factors and
        their multiplicities

        EXAMPLES::

            sage: n = 2^100 - 1; n.factor()                                             # needs sage.libs.pari
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

        We use ``proof=False``, which doesn't prove correctness of the primes
        that appear in the factorization::

            sage: n = 920384092842390423848290348203948092384082349082
            sage: n.factor(proof=False)                                                 # needs sage.libs.pari
            2 * 11 * 1531 * 4402903 * 10023679 * 619162955472170540533894518173
            sage: n.factor(proof=True)                                                  # needs sage.libs.pari
            2 * 11 * 1531 * 4402903 * 10023679 * 619162955472170540533894518173

        We factor using trial division only::

            sage: n.factor(limit=1000)
            2 * 11 * 41835640583745019265831379463815822381094652231

        An example where FLINT is used::

            sage: n = 82862385732327628428164127822
            sage: n.factor(algorithm='flint')                                           # needs sage.libs.flint
            2 * 3 * 11 * 13 * 41 * 73 * 22650083 * 1424602265462161

        We factor using a quadratic sieve algorithm::

            sage: # needs sage.libs.pari
            sage: p = next_prime(10^20)
            sage: q = next_prime(10^21)
            sage: n = p * q
            sage: n.factor(algorithm='qsieve')                                          # needs sage.libs.flint
            doctest:... RuntimeWarning: the factorization returned
            by qsieve may be incomplete (the factors may not be prime)
            or even wrong; see qsieve? for details
            100000000000000000039 * 1000000000000000000117

        We factor using the elliptic curve method::

            sage: # needs sage.libs.pari
            sage: p = next_prime(10^15)
            sage: q = next_prime(10^21)
            sage: n = p * q
            sage: n.factor(algorithm='ecm')
            1000000000000037 * 1000000000000000000117

        TESTS::

            sage: n = 42
            sage: n.factor(algorithm='foobar')
            Traceback (most recent call last):
            ...
            ValueError: Algorithm is not known
        """
        from sage.structure.factorization import Factorization
        from sage.structure.factorization_integer import IntegerFactorization

        if algorithm not in ['pari', 'flint', 'kash', 'magma', 'qsieve', 'ecm']:
            raise ValueError("Algorithm is not known")

        cdef Integer n, p, unit

        if mpz_sgn(self.value) == 0:
            raise ArithmeticError("factorization of 0 is not defined")

        if mpz_sgn(self.value) > 0:
            n = self
            unit = one
        else:
            n = PY_NEW(Integer)
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
            global n_factor_to_list
            if n_factor_to_list is None:
                try:
                    from sage.libs.flint.ulong_extras_sage import n_factor_to_list
                except ImportError:
                    pass
            if n_factor_to_list is not None:
                if proof is None:
                    from sage.structure.proof.proof import get_flag
                    proof = get_flag(proof, "arithmetic")
                F = n_factor_to_list(mpz_get_ui(n.value), proof)
                F = [(smallInteger(a), smallInteger(b)) for a, b in F]
                F.sort()
                return IntegerFactorization(F, unit=unit, unsafe=True,
                                            sort=False, simplify=False)

        if mpz_sizeinbase(n.value, 2) < 40:
            from sage.rings.factorint import factor_trial_division
            return factor_trial_division(self)

        if algorithm == 'pari':
            from sage.rings.factorint_pari import factor_using_pari
            F = factor_using_pari(n, int_=int_, debug_level=verbose, proof=proof)
            F.sort()
            return IntegerFactorization(F, unit=unit, unsafe=True,
                                        sort=False, simplify=False)
        elif algorithm == 'flint':
            from sage.rings.factorint_flint import factor_using_flint
            F = factor_using_flint(n)
            F.sort()
            return IntegerFactorization(F, unit=unit, unsafe=True,
                                        sort=False, simplify=False)

        elif algorithm in ['kash', 'magma']:
            if algorithm == 'kash':
                from sage.interfaces.kash import kash as I
            else:
                from sage.interfaces.magma import magma as I
            str_res = I.eval('Factorization(%s)' % n)
            # The result looks like "[ <n1, p1>, <p2, e2>, ... ]
            str_res = str_res.replace(']', '').replace('[', '').replace('>', '').replace('<', '').split(',')
            res = [int(s.strip()) for s in str_res]
            exp_type = int if int_ else Integer
            F = [(Integer(p), exp_type(e)) for p,e in zip(res[0::2], res[1::2])]
            return Factorization(F, unit)
        elif algorithm == 'qsieve':
            message = "the factorization returned by qsieve may be incomplete (the factors may not be prime) or even wrong; see qsieve? for details"
            from warnings import warn
            warn(message, RuntimeWarning, stacklevel=5)
            from sage.libs.flint.qsieve_sage import qsieve
            F = qsieve(n)
            F.sort()
            return IntegerFactorization(F, unit=unit, unsafe=True,
                                        sort=False, simplify=False)
        else:
            from sage.interfaces.ecm import ecm
            res = [(p, 1) for p in ecm.factor(n, proof=proof)]
            F = IntegerFactorization(res, unit)
            return F

    def support(self):
        """
        Return a sorted list of the primes dividing this integer.

        OUTPUT: the sorted list of primes appearing in the factorization of
        this rational with positive exponent

        EXAMPLES::

            sage: factorial(10).support()
            [2, 3, 5, 7]
            sage: (-999).support()
            [3, 37]

        Trying to find the support of 0 raises an :exc:`ArithmeticError`::

            sage: 0.support()
            Traceback (most recent call last):
            ...
            ArithmeticError: support of 0 not defined
        """
        if self.is_zero():
            raise ArithmeticError("support of 0 not defined")
        from sage.arith.misc import prime_factors

        return prime_factors(self)

    def coprime_integers(self, m):
        """
        Return the nonnegative integers `< m` that are coprime to
        this integer.

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

        TESTS::

            sage: 0.coprime_integers(10^100)
            [1]
            sage: 1.coprime_integers(10^100)
            Traceback (most recent call last):
            ...
            OverflowError: bound is too large
            sage: for n in srange(-6, 7):
            ....:     for m in range(-1, 10):
            ....:         assert n.coprime_integers(m) == [k for k in srange(0, m) if gcd(k, n) == 1]

        AUTHORS:

        - Naqi Jaffery (2006-01-24): examples

        - David Roe (2017-10-02): Use sieving

        - Jeroen Demeyer (2018-06-25): allow returning zero (only relevant for 1.coprime_integers(n))

        ALGORITHM:

        Create an integer with `m` bits and set bits at every multiple
        of a prime `p` that divides this integer and is less than `m`.
        Then return a list of integers corresponding to the unset bits.
        """
        cdef Integer sieve, p, slf, mInteger = Integer(m)
        cdef long k
        cdef unsigned long ilong, plong

        # Trivial case m <= 0
        if mpz_sgn(mInteger.value) <= 0:
            return []

        # Handle 0.coprime_integers(n) first because it's the only case
        # where very large n are allowed
        if mpz_sgn(self.value) == 0:
            if mpz_cmp_ui(mInteger.value, 1) <= 0:
                return []
            return [one]

        if mpz_fits_slong_p(mInteger.value) == 0:
            raise OverflowError("bound is too large")
        cdef long mlong = mpz_get_si(mInteger.value)
        if mpz_cmpabs_ui(self.value, 1) == 0:
            return [smallInteger(k) for k in range(mlong)]
        if (mpz_cmpabs(self.value, mInteger.value) >= 0 and
            (mpz_sgn(self.value) > 0 and self.is_prime() or
             mpz_sgn(self.value) < 0 and (-self).is_prime())):
            return [smallInteger(k) for k in range(1, mlong)]
        sieve = PY_NEW(Integer)
        slf = PY_NEW(Integer)
        mpz_set(slf.value, self.value)
        p = one

        while True:
            sig_check()
            p = slf.trial_division(mlong, mpz_get_si(p.value)+1)
            if mpz_cmp_si(p.value, mlong) >= 0:
                # p is larger than m, so no more primes are needed.
                break
            ilong = plong = mpz_get_ui(p.value)
            while ilong < <unsigned long>mlong:
                # Set bits in sieve at each multiple of p
                mpz_setbit(sieve.value, ilong)
                ilong += plong
            # Now divide by p until no ps remain
            mpz_divexact_ui(slf.value, slf.value, plong)
            while mpz_divisible_ui_p(slf.value, plong):
                mpz_divexact_ui(slf.value, slf.value, plong)
            # If we have found all factors, we break
            if mpz_cmpabs_ui(slf.value, 1) == 0:
                break
        return [smallInteger(k) for k in range(1, mlong)
                if mpz_tstbit(sieve.value, k) == 0]

    def divides(self, n):
        """
        Return ``True`` if ``self`` divides ``n``.

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
        Return the `p`-adic valuation of ``self``.

        We do not require that p be prime, but it must be at least 2. For
        more documentation see ``valuation``

        AUTHORS:

        - David Roe (3/31/07)
        """
        if mpz_sgn(self.value) == 0:
            return sage.rings.infinity.infinity
        if mpz_cmp_ui(p.value, 2) < 0:
            raise ValueError("You can only compute the valuation with respect to a integer larger than 1.")

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
        Return a pair: the `p`-adic valuation of ``self``, and the `p`-adic unit
        of ``self``.

        We do not require the p be prime, but it must be at least 2. For
        more documentation see ``val_unit``

        AUTHORS:

        - David Roe (2007-03-31)
        """
        cdef Integer v, u
        if mpz_cmp_ui(p.value, 2) < 0:
            raise ValueError("You can only compute the valuation with respect to a integer larger than 1.")
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
        Return the `p`-adic valuation of ``self``.

        INPUT:

        - ``p`` -- integer at least 2

        EXAMPLES::

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

        We do not require that ``p`` is a prime::

            sage: (2^11).valuation(4)
            5
        """
        return self._valuation(Integer(p))

    # Alias for valuation
    ord = valuation

    def p_primary_part(self, p):
        """
        Return the ``p``-primary part of ``self``.

        INPUT:

        - ``p`` -- prime integer

        OUTPUT: largest power of ``p`` dividing ``self``

        EXAMPLES::

            sage: n = 40
            sage: n.p_primary_part(2)
            8
            sage: n.p_primary_part(5)
            5
            sage: n.p_primary_part(7)
            1
            sage: n.p_primary_part(6)
            Traceback (most recent call last):
            ...
            ValueError: 6 is not a prime number
        """
        p = smallInteger(p)
        if not p.is_prime():
            raise ValueError("{} is not a prime number".format(p))
        return p**self._valuation(p)

    def val_unit(self, p):
        r"""
        Return a pair: the `p`-adic valuation of ``self``, and th`p`-adicic unit
        of ``self``.

        INPUT:

        - ``p`` -- integer at least 2

        OUTPUT:

        - ``v_p(self)`` -- the `p`-adic valuation of ``self``

        - ``u_p(self)`` -- ``self`` / `p^{v_p(\mathrm{self})}`

        EXAMPLES::

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

        Currently returns 0 when ``self`` is 0.  This behaviour is fairly arbitrary,
        and in Sage 4.6 this special case was not handled at all, eventually
        propagating a :exc:`TypeError`.  The caller should not rely on the behaviour
        in case ``self`` is 0.

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

        odd = PY_NEW(Integer)
        bits = mpz_scan1(self.value, 0)
        mpz_tdiv_q_2exp(odd.value, self.value, bits)
        return odd

    cdef Integer _divide_knowing_divisible_by(Integer self, Integer right):
        r"""
        Return the integer ``self`` / ``right`` when ``self`` is divisible by right.

        If ``self`` is not divisible by right, the return value is undefined,
        and may not even be close to ``self`` / ``right``. For more documentation see
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
        Return the integer ``self`` / ``right`` when ``self`` is divisible by ``right``.

        If ``self`` is not divisible by right, the return value is undefined,
        and may not even be close to ``self`` / ``right`` for multi-word integers.

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
        Return the least common multiple of ``self`` and `n`.

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

    def _gcd(self, Integer n):
        """
        Return the greatest common divisor of ``self`` and `n`.

        EXAMPLES::

            sage: 1._gcd(-1)
            1
            sage: 0._gcd(1)
            1
            sage: 0._gcd(0)
            0
            sage: 2._gcd(2^6)
            2
            sage: 21._gcd(2^6)
            1
        """
        cdef Integer z = PY_NEW(Integer)
        sig_on()
        mpz_gcd(z.value, self.value, n.value)
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

    def as_integer_ratio(self):
        """
        Return the pair ``(self.numerator(), self.denominator())``,
        which is ``(self, 1)``.

        EXAMPLES::

            sage: x = -12
            sage: x.as_integer_ratio()
            (-12, 1)
        """
        return (self, one)

    def factorial(self):
        r"""
        Return the factorial `n! = 1 \cdot 2 \cdot 3 \cdots n`.

        If the input does not fit in an ``unsigned long int``, an :exc:`OverflowError`
        is raised.

        EXAMPLES::

            sage: for n in srange(7):
            ....:     print("{} {}".format(n, n.factorial()))
            0 1
            1 1
            2 2
            3 6
            4 24
            5 120
            6 720

        Large integers raise an :exc:`OverflowError`::

            sage: (2**64).factorial()
            Traceback (most recent call last):
            ...
            OverflowError: argument too large for factorial

        And negative ones a :exc:`ValueError`::

            sage: (-1).factorial()
            Traceback (most recent call last):
            ...
            ValueError: factorial only defined for nonnegative integers
        """
        if mpz_sgn(self.value) < 0:
            raise ValueError("factorial only defined for nonnegative integers")

        if not mpz_fits_ulong_p(self.value):
            raise OverflowError("argument too large for factorial")

        cdef Integer z = PY_NEW(Integer)

        sig_on()
        mpz_fac_ui(z.value, mpz_get_ui(self.value))
        sig_off()

        return z

    def multifactorial(self, long k):
        r"""
        Compute the `k`-th factorial `n!^{(k)}` of ``self``.

        The multifactorial number `n!^{(k)}` is defined for nonnegative
        integers `n` as follows. For `k=1` this is the standard factorial,
        and for `k` greater than `1` it is the product of every `k`-th
        terms down from `n` to `1`. The recursive definition is used to
        extend this function to the negative integers `n`.

        This function uses direct call to GMP if `k` and `n` are nonnegative
        and uses simple transformation for other cases.

        EXAMPLES::

            sage: 5.multifactorial(1)
            120
            sage: 5.multifactorial(2)
            15
            sage: 5.multifactorial(3)
            10

            sage: 23.multifactorial(2)
            316234143225
            sage: prod([1..23, step=2])
            316234143225

            sage: (-29).multifactorial(7)
            1/2640
            sage: (-3).multifactorial(5)
            1/2
            sage: (-9).multifactorial(3)
            Traceback (most recent call last):
            ...
            ValueError: multifactorial undefined

        When entries are too large an :exc:`OverflowError` is raised::

            sage: (2**64).multifactorial(2)
            Traceback (most recent call last):
            ...
            OverflowError: argument too large for multifactorial
        """
        if k <= 0:
            raise ValueError("multifactorial only defined for nonpositive k")

        if not mpz_fits_slong_p(self.value):
            raise OverflowError("argument too large for multifactorial")

        cdef long n = mpz_get_si(self.value)

        cdef Integer z

        if n >= 0:
            # nonnegative n: call native GMP functions
            z = PY_NEW(Integer)
            if k == 1:
                mpz_fac_ui(z.value, n)
            elif k == 2:
                mpz_2fac_ui(z.value, n)
            else:
                mpz_mfac_uiui(z.value, n, k)

            return z

        elif n % k == 0:
            # undefined negative case
            raise ValueError("multifactorial undefined")

        elif -k < n < 0:
            # negative base case
            return one / (self+k)

        # reflection case
        elif n < -k:
            if (n/k) % 2:
                sign = -one
            else:
                sign = one
            return sign / Integer(-k-n).multifactorial(k)

    def gamma(self):
        r"""
        The gamma function on integers is the factorial function (shifted by
        one) on positive integers, and `\pm \infty` on nonpositive integers.

        EXAMPLES::

            sage: # needs sage.symbolic
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
            return sage.rings.infinity.unsigned_infinity

    def floor(self):
        """
        Return the floor of ``self``, which is just ``self`` since ``self`` is an
        integer.

        EXAMPLES::

            sage: n = 6
            sage: n.floor()
            6
        """
        return self

    def ceil(self):
        """
        Return the ceiling of ``self``, which is ``self`` since ``self`` is an
        integer.

        EXAMPLES::

            sage: n = 6
            sage: n.ceil()
            6
        """
        return self

    def trunc(self):
        """
        Round this number to the nearest integer, which is ``self`` since
        ``self`` is an integer.

        EXAMPLES::

            sage: n = 6
            sage: n.trunc()
            6
        """
        return self

    def round(Integer self, mode='away'):
        """
        Return the nearest integer to ``self``, which is ``self`` since
        ``self`` is an integer.

        EXAMPLES:

        This example addresses :issue:`23502`::

            sage: n = 6
            sage: n.round()
            6
        """
        return self

    def real(self):
        """
        Return the real part of ``self``, which is ``self``.

        EXAMPLES::

            sage: Integer(-4).real()
            -4
        """
        return self

    def imag(self):
        """
        Return the imaginary part of ``self``, which is zero.

        EXAMPLES::

            sage: Integer(9).imag()
            0
        """
        return zero

    def is_one(self):
        r"""
        Return ``True`` if the integer is `1`, otherwise ``False``.

        EXAMPLES::

            sage: Integer(1).is_one()
            True
            sage: Integer(0).is_one()
            False
        """
        return mpz_cmp_si(self.value, 1) == 0

    def __bool__(self):
        r"""
        Return ``True`` if the integer is not `0`, otherwise ``False``.

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

    def is_rational(self):
        r"""
        Return ``True`` as an integer is a rational number.

        EXAMPLES::

            sage: 5.is_rational()
            True
        """
        return True

    def is_integer(self):
        """
        Return ``True`` as they are integers.

        EXAMPLES::

            sage: sqrt(4).is_integer()
            True
        """
        return True

    def is_unit(self):
        r"""
        Return ``True`` if this integer is a unit, i.e., `1` or `-1`.

        EXAMPLES::

            sage: for n in srange(-2,3):
            ....:     print("{} {}".format(n, n.is_unit()))
            -2 False
            -1 True
            0 False
            1 True
            2 False
        """
        return mpz_cmpabs_ui(self.value, 1) == 0

    def is_square(self):
        r"""
        Return ``True`` if ``self`` is a perfect square.

        EXAMPLES::

            sage: Integer(4).is_square()
            True
            sage: Integer(41).is_square()
            False
        """
        return mpz_perfect_square_p(self.value)

    def perfect_power(self):
        r"""
        Return ``(a, b)``, where this integer is `a^b` and `b` is maximal.

        If called on `-1`, `0` or `1`, `b` will be `1`, since there is no
        maximal value of `b`.

        .. SEEALSO::

            - :meth:`is_perfect_power`: testing whether an integer is a perfect
              power is usually faster than finding `a` and `b`.
            - :meth:`is_prime_power`: checks whether the base is prime.
            - :meth:`is_power_of`: if you know the base already, this method is
              the fastest option.

        EXAMPLES::

            sage: 144.perfect_power()                                                   # needs sage.libs.pari
            (12, 2)
            sage: 1.perfect_power()
            (1, 1)
            sage: 0.perfect_power()
            (0, 1)
            sage: (-1).perfect_power()
            (-1, 1)
            sage: (-8).perfect_power()                                                  # needs sage.libs.pari
            (-2, 3)
            sage: (-4).perfect_power()
            (-4, 1)
            sage: (101^29).perfect_power()                                              # needs sage.libs.pari
            (101, 29)
            sage: (-243).perfect_power()                                                # needs sage.libs.pari
            (-3, 5)
            sage: (-64).perfect_power()                                                 # needs sage.libs.pari
            (-4, 3)

        TESTS::

            sage: 4.perfect_power()
            (2, 2)
            sage: 256.perfect_power()
            (2, 8)
        """
        cdef long n
        # Fast PARI-free path
        if mpz_fits_slong_p(self.value):
            n = mpz_get_si(self.value)
            if -8 < n < 4:
                return self, one
            if n >= 4:
                if not (n & 1):
                    if mpz_popcount(self.value) == 1:
                        return smallInteger(2), smallInteger(mpz_sizeinbase(self.value, 2) - 1)
                if n < 1000:
                    if _small_primes_table[n >> 1]:
                        return self, one

        parians = self.__pari__().ispower()
        return Integer(parians[1]), Integer(parians[0])

    def global_height(self, prec=None):
        r"""
        Return the absolute logarithmic height of this rational integer.

        INPUT:

        - ``prec`` -- integer; desired floating point precision (default:
          default RealField precision)

        OUTPUT:

        (real) The absolute logarithmic height of this rational integer.

        ALGORITHM:

        The height of the integer `n` is `\log |n|`.

        EXAMPLES::

            sage: # needs sage.rings.real_mpfr
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
            return R.zero()
        return R(self).abs().log()

    cdef bint _is_power_of(Integer self, Integer n) noexcept:
        r"""
        Return a nonzero int if there is an integer b with
        `\mathtt{self} = n^b`.

        For more documentation see :meth:`is_power_of`.

        AUTHORS:

        - David Roe (2007-03-31)
        """
        cdef int a
        cdef unsigned long b, c
        cdef mpz_t u, sabs, nabs
        a = mpz_cmp_ui(n.value, 2)
        if a <= 0:  # n <= 2
            if a == 0:  # n == 2
                if mpz_popcount(self.value) == 1:  # number of bits set in self == 1
                    return 1
                else:
                    return 0
            a = mpz_cmp_si(n.value, -2)
            if a >= 0:  # -2 <= n < 2:
                a = mpz_get_si(n.value)
                if a == 1:  # n == 1
                    if mpz_cmp_ui(self.value, 1) == 0:  # Only 1 is a power of 1
                        return 1
                    else:
                        return 0
                elif a == 0:  # n == 0
                    if mpz_cmp_ui(self.value, 0) == 0 or mpz_cmp_ui(self.value, 1) == 0:  # 0^0 = 1, 0^x = 0
                        return 1
                    else:
                        return 0
                elif a == -1:  # n == -1
                    if mpz_cmp_ui(self.value, 1) == 0 or mpz_cmp_si(self.value, -1) == 0:  # 1 and -1 are powers of -1
                        return 1
                    else:
                        return 0
                elif a == -2:  # n == -2
                    mpz_init(sabs)
                    mpz_abs(sabs, self.value)
                    if mpz_popcount(sabs) == 1:  # number of bits set in |self| == 1
                        b = mpz_scan1(sabs, 0) % 2  # b == 1 if |self| is an odd power of 2, 0 if |self| is an even power
                        mpz_clear(sabs)
                        if (b == 1 and mpz_cmp_ui(self.value, 0) < 0) or (b == 0 and mpz_cmp_ui(self.value, 0) > 0):
                            # An odd power of -2 is negative, an even power must be positive.
                            return 1
                        else:  # number of bits set in |self| is not 1, so self cannot be a power of -2
                            return 0
                    else:  # |self| is not a power of 2, so self cannot be a power of -2
                        return 0
            else:  # n < -2
                mpz_init(nabs)
                mpz_neg(nabs, n.value)
                if mpz_popcount(nabs) == 1:  # |n| = 2^k for k >= 2.  We special case this for speed
                    mpz_init(sabs)
                    mpz_abs(sabs, self.value)
                    if mpz_popcount(sabs) == 1:  # |self| = 2^L for some L >= 0.
                        b = mpz_scan1(sabs, 0)  # the bit that self is set at
                        c = mpz_scan1(nabs, 0)  # the bit that n is set at
                        # Having obtained b and c, we're done with nabs and sabs (on this branch anyway)
                        mpz_clear(nabs)
                        mpz_clear(sabs)
                        if b % c == 0:  # Now we know that |self| is a power of |n|
                            b = (b // c) % 2  # Whether b // c is even or odd determines whether (-2^c)^(b // c) is positive or negative
                            a = mpz_cmp_ui(self.value, 0)
                            if b == 0 and a > 0 or b == 1 and a < 0:
                                # These two cases are that b // c is even and self positive, or b // c is odd and self negative
                                return 1
                            else:  # The sign of self is wrong
                                return 0
                        else:  # Since |self| is not a power of |n|, self cannot be a power of n
                            return 0
                    else:  # self is not a power of 2, and thus cannot be a power of n, which is a power of 2.
                        mpz_clear(nabs)
                        mpz_clear(sabs)
                        return 0
                else:  # |n| is not a power of 2, so we use mpz_remove
                    mpz_init(u)
                    sig_on()
                    b = mpz_remove(u, self.value, nabs)
                    sig_off()
                    # Having obtained b and u, we're done with nabs
                    mpz_clear(nabs)
                    if mpz_cmp_ui(u, 1) == 0:  # self is a power of |n|
                        mpz_clear(u)
                        if b % 2 == 0:  # an even power of |n|, and since self > 0, this means that self is a power of n
                            return 1
                        else:
                            return 0
                    elif mpz_cmp_si(u, -1) == 0:  # -self is a power of |n|
                        mpz_clear(u)
                        if b % 2 == 1:  # an odd power of |n|, and thus self is a power of n
                            return 1
                        else:
                            return 0
                    else:  # |self| is not a power of |n|, so self cannot be a power of n
                        mpz_clear(u)
                        return 0
        elif mpz_popcount(n.value) == 1:  # n > 2 and in fact n = 2^k for k >= 2
            if mpz_popcount(self.value) == 1:  # since n is a power of 2, so must self be.
                if mpz_scan1(self.value, 0) % mpz_scan1(n.value, 0) == 0:  # log_2(self) is divisible by log_2(n)
                    return 1
                else:
                    return 0
            else:  # self is not a power of 2, and thus not a power of n
                return 0
        else:  # n > 2, but not a power of 2, so we use mpz_remove
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
        Return ``True`` if there is an integer `b` with
        `\mathtt{self} = n^b`.

        .. SEEALSO::

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

        .. NOTE::

           For large integers ``self``, :meth:`is_power_of` is faster than
           :meth:`is_perfect_power`. The following examples give some indication of
           how much faster.

        ::

            sage: b = lcm(range(1,10000))
            sage: b.exact_log(2)
            14446
            sage: t = cputime()
            sage: for a in range(2, 1000): k = b.is_perfect_power()
            sage: cputime(t)      # random
            0.53203299999999976
            sage: t = cputime()
            sage: for a in range(2, 1000): k = b.is_power_of(2)
            sage: cputime(t)      # random
            0.0
            sage: t = cputime()
            sage: for a in range(2, 1000): k = b.is_power_of(3)
            sage: cputime(t)      # random
            0.032002000000000308

        ::

            sage: b = lcm(range(1, 1000))
            sage: b.exact_log(2)
            1437
            sage: t = cputime()
            sage: for a in range(2, 10000):  # note: changed range from the example above
            ....:     k = b.is_perfect_power()
            sage: cputime(t)      # random
            0.17201100000000036
            sage: t = cputime(); TWO = int(2)
            sage: for a in range(2, 10000): k = b.is_power_of(TWO)
            sage: cputime(t)      # random
            0.0040000000000000036
            sage: t = cputime()
            sage: for a in range(2, 10000): k = b.is_power_of(3)
            sage: cputime(t)      # random
            0.040003000000000011
            sage: t = cputime()
            sage: for a in range(2, 10000): k = b.is_power_of(a)
            sage: cputime(t)      # random
            0.02800199999999986
        """
        if not isinstance(n, Integer):
            n = Integer(n)
        return self._is_power_of(n)

    def is_prime_power(self, *, proof=None, bint get_data=False):
        r"""
        Return ``True`` if this integer is a prime power, and ``False`` otherwise.

        A prime power is a prime number raised to a positive power. Hence `1` is
        not a prime power.

        For a method that uses a pseudoprimality test instead see
        :meth:`is_pseudoprime_power`.

        INPUT:

        - ``proof`` -- boolean or ``None`` (default). If ``False``, use a strong
          pseudo-primality test (see :meth:`is_pseudoprime`).  If ``True``, use
          a provable primality test. If unset, use the default arithmetic proof
          flag.

        - ``get_data`` -- (default: ``False``), if ``True`` return a pair
          ``(p,k)`` such that this integer equals ``p^k`` with ``p`` a prime
          and ``k`` a positive integer or the pair ``(self,0)`` otherwise.

        .. SEEALSO::

            - :meth:`perfect_power`: Finds the minimal base for which integer
              is a perfect power.
            - :meth:`is_perfect_power`: Doesn't test whether the base is prime.
            - :meth:`is_power_of`: If you know the base already this method is
              the fastest option.
            - :meth:`is_pseudoprime_power`: If the entry is very large.

        EXAMPLES::

            sage: # needs sage.libs.pari
            sage: 17.is_prime_power()
            True
            sage: 10.is_prime_power()
            False
            sage: 64.is_prime_power()
            True
            sage: (3^10000).is_prime_power()
            True
            sage: (10000).is_prime_power()
            False
            sage: (-3).is_prime_power()
            False
            sage: 0.is_prime_power()
            False
            sage: 1.is_prime_power()
            False
            sage: p = next_prime(10^20); p
            100000000000000000039
            sage: p.is_prime_power()
            True
            sage: (p^97).is_prime_power()
            True
            sage: (p + 1).is_prime_power()
            False

        With the ``get_data`` keyword set to ``True``::

            sage: # needs sage.libs.pari
            sage: (3^100).is_prime_power(get_data=True)
            (3, 100)
            sage: 12.is_prime_power(get_data=True)
            (12, 0)
            sage: (p^97).is_prime_power(get_data=True)
            (100000000000000000039, 97)
            sage: q = p.next_prime(); q
            100000000000000000129
            sage: (p*q).is_prime_power(get_data=True)
            (10000000000000000016800000000000000005031, 0)

        The method works for large entries when ``proof=False``::

            sage: proof.arithmetic(False)
            sage: ((10^500 + 961)^4).is_prime_power()                                   # needs sage.libs.pari
            True
            sage: proof.arithmetic(True)

        We check that :issue:`4777` is fixed::

            sage: n = 150607571^14
            sage: n.is_prime_power()                                                    # needs sage.libs.pari
            True

        TESTS::

            sage: 2.is_prime_power(get_data=True)
            (2, 1)
            sage: 4.is_prime_power(get_data=True)
            (2, 2)
            sage: 512.is_prime_power(get_data=True)
            (2, 9)
        """
        cdef long n

        if mpz_sgn(self.value) <= 0:
            return (self, zero) if get_data else False

        if mpz_fits_slong_p(self.value):
            # Fast PARI-free path
            n = mpz_get_si(self.value)
            if not (n & 1):
                if mpz_popcount(self.value) != 1:
                    return (self, zero) if get_data else False
                return (smallInteger(2), smallInteger(mpz_sizeinbase(self.value, 2) - 1)) if get_data else True
            if n < 1000:
                if _small_primes_table[n >> 1]:
                    return (self, one) if get_data else True

            global pari_is_prime_power
            if pari_is_prime_power is None:
                try:
                    from sage.libs.pari.convert_sage import pari_is_prime_power
                except ImportError:
                    pass
            if pari_is_prime_power is not None:
                return pari_is_prime_power(self, get_data)

        if proof is None:
            from sage.structure.proof.proof import get_flag
            proof = get_flag(proof, "arithmetic")

        if proof:
            n, pari_p = self.__pari__().isprimepower()
        else:
            n, pari_p = self.__pari__().ispseudoprimepower()

        if n:
            return (Integer(pari_p), smallInteger(n)) if get_data else True
        else:
            return (self, zero) if get_data else False

    _small_primes_table[:] = [
        0,1,1,1,0,1,1,0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0,1,0,  # 1,3,...,49
        0,1,0,0,1,1,0,0,1,0,1,1,0,0,1,0,1,0,0,1,0,0,0,1,0,  # 51,53,...,99
        1,1,0,1,1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,1,0,0,0,0,1,  # 101,103,...,149
        1,0,0,1,0,0,1,0,1,0,0,1,0,0,1,1,0,0,0,0,1,1,0,1,1,  # 151,153,...,199
        0,0,0,0,0,1,0,0,0,0,0,1,0,1,1,0,1,0,0,1,1,0,0,0,0,  # 201,203,...,249
        1,0,0,1,0,0,1,0,0,1,1,0,0,1,0,1,1,0,0,0,0,1,0,0,0,  # 251,253,...,299
        0,0,0,1,0,1,1,0,1,0,0,0,0,0,0,1,0,0,1,0,0,0,0,1,1,  # 301,303,...,349
        0,1,0,0,1,0,0,0,1,0,0,1,0,0,1,0,1,0,0,1,0,0,0,1,0,  # 351,353,...,399
        1,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,1,0,0,1,0,1,0,0,1,  # 401,403,...,449
        0,0,0,1,0,1,1,0,1,0,0,0,0,0,1,0,0,0,1,0,1,0,0,0,1,  # 451,453,...,499
        0,1,0,0,1,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,1,0,  # 501,503,...,549
        0,0,0,1,0,0,1,0,0,1,1,0,0,1,0,0,0,0,1,0,0,1,0,0,1,  # 551,553,...,599
        1,0,0,1,0,0,1,0,1,1,0,0,0,0,0,1,0,0,0,0,1,1,0,1,0,  # 601,603,...,649
        0,1,0,0,1,1,0,0,0,0,0,1,0,1,0,0,1,0,0,0,1,0,0,0,0,  # 651,653,...,699
        1,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,1,0,0,1,0,1,0,0,0,  # 701,703,...,749
        1,0,0,1,0,1,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,0,1,0,  # 751,753,...,799
        0,0,0,0,1,1,0,0,0,0,1,1,0,1,1,0,0,0,0,1,0,0,0,0,0,  # 801,803,...,849
        0,1,0,1,1,0,1,0,0,0,0,0,0,1,0,1,1,0,1,0,0,0,0,0,0,  # 851,853,...,899
        0,0,0,1,0,1,0,0,0,1,0,0,0,0,1,0,0,0,1,0,1,0,0,1,0,  # 901,903,...,949
        0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,0,1,0,0,0,1,0,0,1,0,  # 951,953,...,999
    ]

    def is_prime(self, proof=None):
        r"""
        Test whether ``self`` is prime.

        INPUT:

        - ``proof`` -- boolean or ``None`` (default). If ``False``, use a
          strong pseudo-primality test (see :meth:`is_pseudoprime`).
          If ``True``, use a provable primality test.  If unset, use the
          :mod:`default arithmetic proof flag <sage.structure.proof.proof>`.

        .. NOTE::

           Integer primes are by definition *positive*! This is
           different than Magma, but the same as in PARI. See also the
           :meth:`is_irreducible()` method.

        EXAMPLES::

            sage: z = 2^31 - 1
            sage: z.is_prime()                                                          # needs sage.libs.pari
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
            sage: z.is_prime(proof=False)                                               # needs sage.libs.pari
            True
            sage: z.is_prime(proof=True)                                                # needs sage.libs.pari
            True

        When starting Sage the arithmetic proof flag is True. We can change
        it to False as follows::

            sage: proof.arithmetic()
            True
            sage: n = 10^100 + 267
            sage: timeit("n.is_prime()")        # not tested                            # needs sage.libs.pari
            5 loops, best of 3: 163 ms per loop
            sage: proof.arithmetic(False)
            sage: proof.arithmetic()
            False
            sage: timeit("n.is_prime()")        # not tested                            # needs sage.libs.pari
            1000 loops, best of 3: 573 us per loop

        ALGORITHM:

        Calls the PARI function :pari:`isprime`.

        TESTS:

        We compare the output of this method to a straightforward sieve::

            sage: size = 10000
            sage: tab = [0,0] + [1] * (size-2)
            sage: for i in range(size):
            ....:     if tab[i]:
            ....:         for j in range(2*i, size, i):
            ....:             tab[j] = 0
            sage: all(ZZ(i).is_prime() == b for i,b in enumerate(tab))                  # needs sage.libs.pari
            True
        """
        if mpz_sgn(self.value) <= 0:
            return False

        cdef unsigned long u
        if mpz_fits_ulong_p(self.value):
            u = mpz_get_ui(self.value)
            if not (u & 1):
                return u == 2
            if u < 1000:
                return _small_primes_table[u >> 1]

            global pari_is_prime
            if pari_is_prime is None:
                try:
                    from sage.libs.pari.convert_sage import pari_is_prime
                except ImportError:
                    pass
            if pari_is_prime is not None:
                return pari_is_prime(self)

        if proof is None:
            from sage.structure.proof.proof import get_flag
            proof = get_flag(proof, "arithmetic")
        if proof:
            return self.__pari__().isprime()
        else:
            return self.__pari__().ispseudoprime()

    cdef bint _pseudoprime_is_prime(self, proof) except -1:
        """
        Given a pseudoprime, return ``self.is_prime(proof)``.

        INPUT:

        - ``self`` -- a PARI pseudoprime

        - ``proof`` -- mandatory proof flag (``True``, ``False`` or ``None``)

        OUTPUT:

        - The result of ``self.is_prime(proof)`` but faster
        """
        if mpz_cmp(self.value, PARI_PSEUDOPRIME_LIMIT) < 0:
            return True
        if proof is None:
            from sage.structure.proof.proof import get_flag
            proof = get_flag(proof, "arithmetic")
        if proof:
            return self.__pari__().isprime()
        else:
            return True

    def is_irreducible(self):
        r"""
        Return ``True`` if ``self`` is irreducible, i.e. +/-
        prime

        EXAMPLES::

            sage: z = 2^31 - 1
            sage: z.is_irreducible()                                                    # needs sage.libs.pari
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
        cdef Integer n = self if self >= 0 else -self
        return n.is_prime(proof=True)

    def is_pseudoprime(self):
        r"""
        Test whether ``self`` is a pseudoprime.

        This uses PARI's Baillie-PSW probabilistic primality
        test. Currently, there are no known pseudoprimes for
        Baillie-PSW that are not actually prime. However, it is
        conjectured that there are infinitely many.

        See :wikipedia:`Baillie-PSW_primality_test`

        EXAMPLES::

            sage: z = 2^31 - 1
            sage: z.is_pseudoprime()                                                    # needs sage.libs.pari
            True
            sage: z = 2^31
            sage: z.is_pseudoprime()                                                    # needs sage.libs.pari
            False
        """
        return self.__pari__().ispseudoprime()

    def is_pseudoprime_power(self, get_data=False):
        r"""
        Test if this number is a power of a pseudoprime number.

        For large numbers, this method might be faster than
        :meth:`is_prime_power`.

        INPUT:

        - ``get_data`` -- (default: ``False``) if ``True`` return a pair `(p,k)`
          such that this number equals `p^k` with `p` a pseudoprime and `k` a
          positive integer or the pair ``(self,0)`` otherwise

        EXAMPLES::

            sage: # needs sage.libs.pari
            sage: x = 10^200 + 357
            sage: x.is_pseudoprime()
            True
            sage: (x^12).is_pseudoprime_power()
            True
            sage: (x^12).is_pseudoprime_power(get_data=True)
            (1000...000357, 12)
            sage: (997^100).is_pseudoprime_power()
            True
            sage: (998^100).is_pseudoprime_power()
            False
            sage: ((10^1000 + 453)^2).is_pseudoprime_power()
            True

        TESTS::

            sage: 0.is_pseudoprime_power()
            False
            sage: (-1).is_pseudoprime_power()
            False
            sage: 1.is_pseudoprime_power()                                              # needs sage.libs.pari
            False
        """
        return self.is_prime_power(proof=False, get_data=get_data)

    def is_perfect_power(self):
        r"""
        Return ``True`` if ``self`` is a perfect power, ie if there exist integers
        `a` and `b`, `b > 1` with ``self`` `= a^b`.

        .. SEEALSO::

            - :meth:`perfect_power`: Finds the minimal base for which this
              integer is a perfect power.
            - :meth:`is_power_of`: If you know the base already, this method is
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

        TESTS:

        This is a test to make sure we work around a bug in GMP, see
        :issue:`4612`.

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
            return res != 0
        return mpz_perfect_power_p(self.value)

    def is_norm(self, K, element=False, proof=True):
        r"""
        See ``QQ(self).is_norm()``.

        EXAMPLES::

            sage: n = 7
            sage: n.is_norm(QQ)
            True
            sage: n.is_norm(QQ, element=True)
            (True, 7)

            sage: # needs sage.rings.number_field
            sage: x = polygen(ZZ, 'x')
            sage: K = NumberField(x^2 - 2, 'beta')
            sage: n = 4
            sage: n.is_norm(K)
            True
            sage: 5.is_norm(K)
            False
            sage: n.is_norm(K, element=True)
            (True, -4*beta + 6)
            sage: n.is_norm(K, element=True)[1].norm()
            4
            sage: n = 5
            sage: n.is_norm(K, element=True)
            (False, None)
        """
        from sage.rings.rational_field import QQ
        return QQ(self).is_norm(K, element=element, proof=proof)

    def _bnfisnorm(self, K, proof=True, extra_primes=0):
        r"""
        See ``QQ(self)._bnfisnorm()``.

        EXAMPLES::

            sage: 3._bnfisnorm(QuadraticField(-1, 'i'))                                 # needs sage.rings.number_field
            (1, 3)
            sage: 7._bnfisnorm(CyclotomicField(7))                                      # needs sage.rings.number_field
            (zeta7^5 - zeta7^2, 1)
        """
        from sage.rings.rational_field import QQ
        return QQ(self)._bnfisnorm(K, proof=proof, extra_primes=extra_primes)

    def jacobi(self, b):
        r"""
        Calculate the Jacobi symbol `\left(\frac{\text{self}}{b}\right)`.

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
        if is_small_python_int(b):
            tmp = b
            if (tmp & 1) == 0:
                raise ValueError("Jacobi symbol not defined for even b.")
            return mpz_kronecker_si(self.value, tmp)
        if not isinstance(b, Integer):
            b = Integer(b)
        if mpz_even_p((<Integer>b).value):
            raise ValueError("Jacobi symbol not defined for even b.")
        return mpz_jacobi(self.value, (<Integer>b).value)

    def kronecker(self, b):
        r"""
        Calculate the Kronecker symbol `\left(\frac{\text{self}}{b}\right)`
        with the Kronecker extension `(\text{self}/2)=(2/\text{self})` when ``self`` is odd,
        or `(\text{self}/2)=0` when ``self`` is even.

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
        if is_small_python_int(b):
            return mpz_kronecker_si(self.value, b)
        if not isinstance(b, Integer):
            b = Integer(b)
        return mpz_kronecker(self.value, (<Integer>b).value)

    def class_number(self, proof=True):
        r"""
        Return the class number of the quadratic order with this discriminant.

        INPUT:

        - ``self`` -- integer congruent to `0` or `1` mod `4` which is
          not a square

        - ``proof`` -- boolean (default: ``True``); if ``False``, then
          for negative discriminants a faster algorithm is used by
          the PARI library which is known to give incorrect results
          when the class group has many cyclic factors.  However, the
          results are correct for discriminants `D` with `|D|\le 2\cdot10^{10}`.

        OUTPUT:

        (integer) the class number of the quadratic order with this
        discriminant.

        .. NOTE::

           For positive `D`, this is not always equal to the number of classes of
           primitive binary quadratic forms of discriminant `D`, which
           is equal to the narrow class number. The two notions are
           the same when `D<0`, or `D>0` and the fundamental unit of
           the order has negative norm; otherwise the number of
           classes of forms is twice this class number.

        EXAMPLES::

            sage: (-163).class_number()                                                 # needs sage.libs.pari
            1
            sage: (-104).class_number()                                                 # needs sage.libs.pari
            6
            sage: [((4*n + 1), (4*n + 1).class_number()) for n in [21..29]]             # needs sage.libs.pari
            [(85, 2),
            (89, 1),
            (93, 1),
            (97, 1),
            (101, 1),
            (105, 2),
            (109, 1),
            (113, 1),
            (117, 1)]

        TESTS:

        The integer must not be a square, or an error is raised::

           sage: 100.class_number()
           Traceback (most recent call last):
           ...
           ValueError: class_number not defined for square integers


        The integer must be 0 or 1 mod 4, or an error is raised::

           sage: 10.class_number()
           Traceback (most recent call last):
           ...
           ValueError: class_number only defined for integers congruent to 0 or 1 modulo 4
           sage: 3.class_number()
           Traceback (most recent call last):
           ...
           ValueError: class_number only defined for integers congruent to 0 or 1 modulo 4
        """
        if self.is_square():
            raise ValueError("class_number not defined for square integers")
        if self % 4 not in [0, 1]:
            raise ValueError("class_number only defined for integers congruent to 0 or 1 modulo 4")

        global objtogen
        if objtogen is None:
            from cypari2.gen import objtogen
        flag = self < 0 and proof
        return objtogen(self).qfbclassno(flag).sage()

    def squarefree_part(self, long bound=-1):
        r"""
        Return the square free part of `x` (=``self``), i.e., the unique integer
        `z` that `x = z y^2`, with `y^2` a perfect square and `z` square-free.

        Use ``self.radical()`` for the product of the primes that divide ``self``.

        If ``self`` is 0, just returns 0.

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
            sage: a = 7^3 * next_prime(2^100)^2 * next_prime(2^200)                     # needs sage.libs.pari
            sage: a / a.squarefree_part(bound=1000)                                     # needs sage.libs.pari
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
                if e % 2:
                    n = n * pp
            return n * F.unit()

    def next_probable_prime(self):
        """
        Return the next probable prime after ``self``, as determined by PARI.

        EXAMPLES::

            sage: # needs sage.libs.pari
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
        return Integer(self.__pari__().nextprime(True))

    def next_prime(self, proof=None):
        r"""
        Return the next prime after ``self``.

        This method calls the PARI function :pari:`nextprime`.

        INPUT:

        - ``proof`` -- boolean or ``None`` (default: ``None``, see
          ``proof.arithmetic`` or :mod:`sage.structure.proof`); note that the
          global Sage default is ``proof=True``

        EXAMPLES::

            sage: 100.next_prime()                                                      # needs sage.libs.pari
            101
            sage: (10^50).next_prime()                                                  # needs sage.libs.pari
            100000000000000000000000000000000000000000000000151

        Use ``proof=False``, which is way faster since it does not need
        a primality proof::

            sage: b = (2^1024).next_prime(proof=False)                                  # needs sage.libs.pari
            sage: b - 2^1024                                                            # needs sage.libs.pari
            643

        ::

            sage: Integer(0).next_prime()                                               # needs sage.libs.pari
            2
            sage: Integer(1001).next_prime()                                            # needs sage.libs.pari
            1009
        """
        # Use PARI to compute the next *pseudo*-prime
        p = Integer(self.__pari__().nextprime(True))
        while not p._pseudoprime_is_prime(proof):
            p = Integer(p.__pari__().nextprime(True))
        return p

    def previous_prime(self, proof=None):
        r"""
        Return the previous prime before ``self``.

        This method calls the PARI function :pari:`precprime`.

        INPUT:

        - ``proof`` -- if ``True`` ensure that the returned value is the next
          prime power and if set to ``False`` uses probabilistic methods
          (i.e. the result is not guaranteed). By default it uses global
          configuration variables to determine which alternative to use (see
          :mod:`proof.arithmetic` or :mod:`sage.structure.proof`).

        .. SEEALSO::

            - :meth:`next_prime`

        EXAMPLES::

            sage: 10.previous_prime()                                                   # needs sage.libs.pari
            7
            sage: 7.previous_prime()                                                    # needs sage.libs.pari
            5
            sage: 14376485.previous_prime()                                             # needs sage.libs.pari
            14376463

            sage: 2.previous_prime()
            Traceback (most recent call last):
            ...
            ValueError: no prime less than 2

        An example using ``proof=False``, which is way faster since it does not
        need a primality proof::

            sage: b = (2^1024).previous_prime(proof=False)                              # needs sage.libs.pari
            sage: 2^1024 - b                                                            # needs sage.libs.pari
            105
        """
        if mpz_cmp_ui(self.value, 2) <= 0:
            raise ValueError("no prime less than 2")
        cdef Integer p = self-1
        p = Integer(p.__pari__().precprime())
        while not p._pseudoprime_is_prime(proof):
            mpz_sub_ui(p.value, p.value, 1)
            p = Integer(p.__pari__().precprime())
        return p

    def next_prime_power(self, proof=None):
        r"""
        Return the next prime power after ``self``.

        INPUT:

        - ``proof`` -- if ``True`` ensure that the returned value is the next
          prime power and if set to ``False`` uses probabilistic methods
          (i.e. the result is not guaranteed). By default it uses global
          configuration variables to determine which alternative to use (see
          :mod:`proof.arithmetic` or :mod:`sage.structure.proof`).

        ALGORITHM:

        The algorithm is naive. It computes the next power of 2 and goes through
        the odd numbers calling :meth:`is_prime_power`.

        .. SEEALSO::

            - :meth:`previous_prime_power`
            - :meth:`is_prime_power`
            - :meth:`next_prime`
            - :meth:`previous_prime`

        EXAMPLES::

            sage: (-1).next_prime_power()
            2
            sage: 2.next_prime_power()
            3
            sage: 103.next_prime_power()                                                # needs sage.libs.pari
            107
            sage: 107.next_prime_power()
            109
            sage: 2044.next_prime_power()                                               # needs sage.libs.pari
            2048

        TESTS::

            sage: [(2**k - 1).next_prime_power() for k in range(1,10)]
            [2, 4, 8, 16, 32, 64, 128, 256, 512]
            sage: [(2**k).next_prime_power() for k in range(10)]                        # needs sage.libs.pari
            [2, 3, 5, 9, 17, 37, 67, 131, 257, 521]

            sage: for _ in range(10):                                                   # needs sage.libs.pari
            ....:     n = ZZ.random_element(2**256).next_prime_power()
            ....:     m = n.next_prime_power().previous_prime_power()
            ....:     assert m == n, "problem with n = {}".format(n)
        """
        if mpz_cmp_ui(self.value, 2) < 0:
            return smallInteger(2)

        cdef mp_bitcnt_t bit_index = mpz_sizeinbase(self.value, 2)
        cdef Integer n = PY_NEW(Integer)

        mpz_add_ui(n.value, self.value, 1 if mpz_even_p(self.value) else 2)

        while not mpz_tstbit(n.value, bit_index):
            if n.is_prime_power(proof=proof):
                return n
            mpz_add_ui(n.value, n.value, 2)

        # return the power of 2 we just skipped
        mpz_sub_ui(n.value, n.value, 1)
        return n

    def previous_prime_power(self, proof=None):
        r"""
        Return the previous prime power before ``self``.

        INPUT:

        - ``proof`` -- if ``True`` ensure that the returned value is the next
          prime power and if set to ``False`` uses probabilistic methods
          (i.e. the result is not guaranteed). By default it uses global
          configuration variables to determine which alternative to use (see
          :mod:`proof.arithmetic` or :mod:`sage.structure.proof`).

        ALGORITHM:

        The algorithm is naive. It computes the previous power of 2 and goes
        through the odd numbers calling the method :meth:`is_prime_power`.

        .. SEEALSO::

            - :meth:`next_prime_power`
            - :meth:`is_prime_power`
            - :meth:`previous_prime`
            - :meth:`next_prime`

        EXAMPLES::

            sage: # needs sage.libs.pari
            sage: 3.previous_prime_power()
            2
            sage: 103.previous_prime_power()
            101
            sage: 107.previous_prime_power()
            103
            sage: 2044.previous_prime_power()
            2039

            sage: 2.previous_prime_power()
            Traceback (most recent call last):
            ...
            ValueError: no prime power less than 2

        TESTS::

            sage: [(2**k + 1).previous_prime_power() for k in range(1,10)]
            [2, 4, 8, 16, 32, 64, 128, 256, 512]
            sage: [(2**k).previous_prime_power() for k in range(2, 10)]                 # needs sage.libs.pari
            [3, 7, 13, 31, 61, 127, 251, 509]

            sage: for _ in range(10):                                                   # needs sage.libs.pari
            ....:     n = ZZ.random_element(3,2**256).previous_prime_power()
            ....:     m = n.previous_prime_power().next_prime_power()
            ....:     assert m == n, "problem with n = {}".format(n)
        """
        if mpz_cmp_ui(self.value, 2) <= 0:
            raise ValueError("no prime power less than 2")

        cdef Integer n = PY_NEW(Integer)

        mpz_sub_ui(n.value, self.value, 1)
        cdef mp_bitcnt_t bit_index = mpz_sizeinbase(n.value, 2)-1
        if mpz_even_p(n.value):
            mpz_sub_ui(n.value, n.value, 1)

        while mpz_tstbit(n.value, bit_index):
            if n.is_prime_power(proof=proof):
                return n
            mpz_sub_ui(n.value, n.value, 2)

        # return the power of 2 we just skipped
        mpz_add_ui(n.value, n.value, 1)
        return n

    def additive_order(self):
        """
        Return the additive order of ``self``.

        EXAMPLES::

            sage: ZZ(0).additive_order()
            1
            sage: ZZ(1).additive_order()
            +Infinity
        """
        if mpz_sgn(self.value) == 0:
            return one
        else:
            return sage.rings.infinity.infinity

    def multiplicative_order(self):
        r"""
        Return the multiplicative order of ``self``.

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
        if mpz_cmp_si(self.value, 1) == 0:
            return one
        elif mpz_cmp_si(self.value, -1) == 0:
            return smallInteger(2)
        else:
            return sage.rings.infinity.infinity

    def is_squarefree(self):
        """
        Return ``True`` if this integer is not divisible by the square of any
        prime and ``False`` otherwise.

        EXAMPLES::

            sage: 100.is_squarefree()                                                   # needs sage.libs.pari
            False
            sage: 102.is_squarefree()                                                   # needs sage.libs.pari
            True
            sage: 0.is_squarefree()                                                     # needs sage.libs.pari
            False
        """
        return self.__pari__().issquarefree()

    def is_discriminant(self):
        """
        Return ``True`` if this integer is a discriminant.

        .. NOTE::

            A discriminant is an integer congruent to 0 or 1 modulo 4.

        EXAMPLES::

            sage: (-1).is_discriminant()
            False
            sage: (-4).is_discriminant()
            True
            sage: 100.is_discriminant()
            True
            sage: 101.is_discriminant()
            True

        TESTS::

            sage: 0.is_discriminant()
            True
            sage: 1.is_discriminant()
            True
            sage: len([D for D in srange(-100,100) if D.is_discriminant()])
            100
        """
        return self % 4 in [0, 1]

    def is_fundamental_discriminant(self):
        """
        Return ``True`` if this integer is a fundamental discriminant.

        .. NOTE::

            A fundamental discriminant is a discrimimant, not 0 or 1 and not a square multiple of a smaller discriminant.

        EXAMPLES::

            sage: (-4).is_fundamental_discriminant()                                    # needs sage.libs.pari
            True
            sage: (-12).is_fundamental_discriminant()
            False
            sage: 101.is_fundamental_discriminant()                                     # needs sage.libs.pari
            True

        TESTS::

            sage: 0.is_fundamental_discriminant()
            False
            sage: 1.is_fundamental_discriminant()
            False
            sage: len([D for D in srange(-100,100)                                      # needs sage.libs.pari
            ....:      if D.is_fundamental_discriminant()])
            61
        """
        if self in [0, 1]:
            return False
        Dmod4 = self % 4
        if Dmod4 in [2, 3]:
            return False
        if Dmod4 == 1:
            return self.is_squarefree()
        d = self // 4
        return d % 4 in [2, 3] and d.is_squarefree()

    cpdef __pari__(self):
        """
        Return the PARI version of this integer.

        EXAMPLES::

            sage: n = 9390823
            sage: m = n.__pari__(); m                                                   # needs sage.libs.pari
            9390823
            sage: type(m)                                                               # needs sage.libs.pari
            <class 'cypari2.gen.Gen'>

        TESTS::

            sage: n = 10^10000000
            sage: m = n.__pari__()  # crash from trac 875                               # needs sage.libs.pari
            sage: m % 1234567                                                           # needs sage.libs.pari
            1041334
        """
        global new_gen_from_integer
        if new_gen_from_integer is None:
            from sage.libs.pari.convert_sage import new_gen_from_integer
        return new_gen_from_integer(self)

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

    @property
    def __array_interface__(self):
        """
        Used for NumPy conversion.

        EXAMPLES::

            sage: # needs numpy
            sage: import numpy
            sage: numpy.array([1, 2, 3])
            array([1, 2, 3])
            sage: numpy.array([1, 2, 3]).dtype
            dtype('int32')                         # 32-bit
            dtype('int64')                         # 64-bit

            sage: # needs numpy (this has to be repeated until #36099 is fixed)
            sage: import numpy
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

        EXAMPLES::

            sage: (117)._magma_init_(magma)           # optional - magma
            '117'

        Large integers use hex:
            sage: # optional - magma
            sage: m = 3^(2^20)
            sage: s = m._magma_init_(magma)
            sage: 'StringToInteger' in s
            True
            sage: magma(m).sage() == m
            True
        """
        if self.ndigits(2) > 10000:
            return 'StringToInteger("%s",16)' % self.str(16)
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
        Return `(s, r)` where `s` is the integer square root of ``self`` and
        `r` is the remainder such that `\text{self} = s^2 + r`.
        Raises :exc:`ValueError` if ``self`` is negative.

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
            ValueError: square root of negative integer not defined
        """
        if mpz_sgn(self.value) < 0:
            raise ValueError("square root of negative integer not defined")
        cdef Integer s = PY_NEW(Integer)
        cdef Integer r = PY_NEW(Integer)
        mpz_sqrtrem(s.value, r.value, self.value)
        return s, r

    def isqrt(self):
        r"""
        Return the integer floor of the square root of ``self``, or raises an
        :exc:`ValueError` if ``self`` is negative.

        EXAMPLES::

            sage: a = Integer(5)
            sage: a.isqrt()
            2

        ::

            sage: Integer(-102).isqrt()
            Traceback (most recent call last):
            ...
            ValueError: square root of negative integer not defined
        """
        if mpz_sgn(self.value) < 0:
            raise ValueError("square root of negative integer not defined")

        cdef Integer x = PY_NEW(Integer)

        sig_on()
        mpz_sqrt(x.value, self.value)
        sig_off()

        return x

    def sqrt(self, prec=None, extend=True, all=False):
        """
        The square root function.

        INPUT:

        - ``prec`` -- integer (default: ``None``); if ``None``, return an exact
          square root; otherwise return a numerical square root, to the
          given bits of precision.

        - ``extend`` -- boolean (default: ``True``); if ``True``, return a
          square root in an extension ring, if necessary. Otherwise, raise a
          :exc:`ValueError` if the square is not in the base ring. Ignored if
          ``prec`` is not ``None``.

        - ``all`` -- boolean (default: ``False``); if ``True``, return all
          square roots of ``self`` (a list of length 0, 1, or 2)

        EXAMPLES::

            sage: Integer(144).sqrt()
            12
            sage: sqrt(Integer(144))
            12
            sage: Integer(102).sqrt()                                                   # needs sage.symbolic
            sqrt(102)

        ::

            sage: n = 2
            sage: n.sqrt(all=True)                                                      # needs sage.symbolic
            [sqrt(2), -sqrt(2)]
            sage: n.sqrt(prec=10)                                                       # needs sage.rings.real_mpfr
            1.4
            sage: n.sqrt(prec=100)                                                      # needs sage.rings.real_mpfr
            1.4142135623730950488016887242
            sage: n.sqrt(prec=100, all=True)                                            # needs sage.rings.real_mpfr
            [1.4142135623730950488016887242, -1.4142135623730950488016887242]
            sage: n.sqrt(extend=False)
            Traceback (most recent call last):
            ...
            ArithmeticError: square root of 2 is not an integer
            sage: (-1).sqrt(extend=False)
            Traceback (most recent call last):
            ...
            ArithmeticError: square root of -1 is not an integer
            sage: Integer(144).sqrt(all=True)
            [12, -12]
            sage: Integer(0).sqrt(all=True)
            [0]

        TESTS::

            sage: type(5.sqrt())                                                        # needs sage.symbolic
            <class 'sage.symbolic.expression.Expression'>
            sage: type(5.sqrt(prec=53))                                                 # needs sage.rings.real_mpfr
            <class 'sage.rings.real_mpfr.RealNumber'>
            sage: type((-5).sqrt(prec=53))                                              # needs sage.rings.real_mpfr
            <class 'sage.rings.complex_mpfr.ComplexNumber'>
            sage: type(0.sqrt(prec=53))                                                 # needs sage.rings.real_mpfr
            <class 'sage.rings.real_mpfr.RealNumber'>

        Check that :issue:`9466` and :issue:`26509` are fixed::

            sage: 3.sqrt(extend=False, all=True)
            []
            sage: (-1).sqrt(extend=False, all=True)
            []
        """
        if prec is not None:
            from sage.misc.functional import _do_sqrt
            return _do_sqrt(self, prec=prec, all=all)

        if mpz_sgn(self.value) == 0:
            return [self] if all else self

        cdef bint is_square
        cdef Integer z
        cdef mpz_t tmp
        if mpz_sgn(self.value) < 0:
            is_square = False
        else:
            sig_on()
            mpz_init(tmp)
            z = PY_NEW(Integer)
            mpz_sqrtrem(z.value, tmp, self.value)
            is_square = (mpz_sgn(tmp) == 0)
            mpz_clear(tmp)
            sig_off()

        if not is_square:
            if extend:
                from sage.misc.functional import _do_sqrt
                return _do_sqrt(self, all=all)
            if all:
                return []
            raise ArithmeticError(f"square root of {self} is not an integer")

        if all:
            return [z, -z]
        return z

    @coerce_binop
    def xgcd(self, Integer n):
        r"""
        Return the extended gcd of this element and ``n``.

        INPUT:

        - ``n`` -- integer

        OUTPUT:

        A triple ``(g, s, t)`` such that ``g`` is the nonnegative gcd of
        ``self`` and ``n``, and ``s`` and ``t`` are cofactors satisfying the
        Bezout identity

        .. MATH::

            g = s \cdot \mathrm{self} + t \cdot n.

        .. NOTE::

            There is no guarantee that the cofactors will be minimal. If you
            need the cofactors to be minimal use :meth:`_xgcd`. Also, using
            :meth:`_xgcd` directly might be faster in some cases, see
            :issue:`13628`.

        EXAMPLES::

            sage: 6.xgcd(4)
            (2, 1, -1)
        """
        return self._xgcd(n)

    def _xgcd(self, Integer n, bint minimal=0):
        r"""
        Return the extended gcd of ``self`` and ``n``.

        INPUT:

        - ``n`` -- integer
        - ``minimal`` -- boolean (default: ``False``); whether to compute
          minimal cofactors (see below)

        OUTPUT:

        A triple ``(g, s, t)`` such that ``g`` is the nonnegative gcd of
        ``self`` and ``n``, and ``s`` and ``t`` are cofactors satisfying the
        Bezout identity.

        .. MATH::

            g = s \cdot \mathrm{self} + t \cdot n.

        .. NOTE::

            If ``minimal`` is ``False``, then there is no guarantee that the
            returned cofactors will be minimal in any sense; the only guarantee
            is that the Bezout identity will be satisfied (see examples below).

            If ``minimal`` is ``True``, the cofactors will satisfy the following
            conditions. If either ``self`` or ``n`` are zero, the trivial
            solution is returned. If both ``self`` and ``n`` are nonzero, the
            function returns the unique solution such that `0 \leq s < |n|/g`
            (which then must also satisfy
            `0 \leq |t| \leq |\mbox{\rm self}|/g`).

        EXAMPLES::

            sage: 5._xgcd(7)
            (1, 3, -2)
            sage: 5*3 + 7*-2
            1
            sage: g,s,t = 58526524056._xgcd(101294172798)
            sage: g
            22544886
            sage: 58526524056 * s + 101294172798 * t
            22544886

        Try ``minimal`` option with various edge cases::

            sage: 5._xgcd(0, minimal=True)
            (5, 1, 0)
            sage: (-5)._xgcd(0, minimal=True)
            (5, -1, 0)
            sage: 0._xgcd(5, minimal=True)
            (5, 0, 1)
            sage: 0._xgcd(-5, minimal=True)
            (5, 0, -1)
            sage: 0._xgcd(0, minimal=True)
            (0, 1, 0)

        Output may differ with and without the ``minimal`` option::

            sage: 5._xgcd(6)
            (1, -1, 1)
            sage: 5._xgcd(6, minimal=True)
            (1, 5, -4)

        Exhaustive tests, checking minimality conditions::

            sage: for a in srange(-20, 20):
            ....:   for b in srange(-20, 20):
            ....:     if a == 0 or b == 0: continue
            ....:     g, s, t = a._xgcd(b)
            ....:     assert g > 0
            ....:     assert a % g == 0 and b % g == 0
            ....:     assert a*s + b*t == g
            ....:     g, s, t = a._xgcd(b, minimal=True)
            ....:     assert g > 0
            ....:     assert a % g == 0 and b % g == 0
            ....:     assert a*s + b*t == g
            ....:     assert s >= 0 and s < abs(b)/g
            ....:     assert abs(t) <= abs(a)/g

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
        # satisfied by the returned cofactors. They guarantee a nonnegative
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
        Compute left and right shifts of integers.
        Shifts ``self`` ``y`` bits to the left if ``sign`` is `1`, and to the right
        if ``sign`` is `-1`.

        WARNING: This function does no error checking. In particular,
        it assumes that ``sign`` is either `1` or `-1`.

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

        TESTS::

            sage: try:
            ....:     print('Possible error output from gmp', flush=True)
            ....:     1 << (2^60)
            ....: except (MemoryError, OverflowError, RuntimeError, FloatingPointError):
            ....:     pass
            ....: else:
            ....:     print("Failed to raise exception")
            Possible error output from gmp...
        """
        cdef long n

        if type(y) is int:
            # For a Python int, we can just use the Python/C API.
            n = PyLong_AsLong(y)
        else:
            # If it's not already an Integer, try to convert it.
            if not isinstance(y, Integer):
                try:
                    y = Integer(y)
                except TypeError:
                    raise TypeError("unsupported operands for %s: %s, %s" % (("<<" if sign == 1 else ">>"), self, y))
                except ValueError:
                    return coercion_model.bin_op(self, y, operator.lshift if sign == 1 else operator.rshift)

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
        sig_on()
        if n < 0:
            mpz_fdiv_q_2exp(z.value, self.value, -n)
        else:
            mpz_mul_2exp(z.value, self.value, n)
        sig_off()
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
            sage: 1 << 2.5                                                              # needs sage.rings.real_mpfr
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
        if not isinstance(x, Integer):
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
            sage: 1 >> 2.5                                                              # needs sage.rings.real_mpfr
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
        if not isinstance(x, Integer):
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
        if isinstance(x, Integer) and isinstance(y, Integer):
            return (<Integer>x)._and(y)
        return coercion_model.bin_op(x, y, operator.and_)

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
        if isinstance(x, Integer) and isinstance(y, Integer):
            return (<Integer>x)._or(y)
        return coercion_model.bin_op(x, y, operator.or_)

    def __invert__(self):
        """
        Return the multiplicative inverse of ``self``, as a rational number.

        EXAMPLES::

            sage: n = 10
            sage: 1/n
            1/10
            sage: n.__invert__()
            1/10
            sage: n = -3
            sage: ~n
            -1/3
        """
        if mpz_sgn(self.value) == 0:
            raise ZeroDivisionError("rational division by zero")
        cdef Rational x
        x = <Rational> Rational.__new__(Rational)
        mpz_set_ui(mpq_numref(x.value), 1)
        mpz_set(mpq_denref(x.value), self.value)
        if mpz_sgn(self.value) == -1:
            mpz_neg(mpq_numref(x.value), mpq_numref(x.value))
            mpz_neg(mpq_denref(x.value), mpq_denref(x.value))
        return x

    def inverse_of_unit(self):
        """
        Return inverse of ``self`` if ``self`` is a unit in the integers, i.e.,
        ``self`` is `-1` or `1`. Otherwise, raise a :exc:`ZeroDivisionError`.

        EXAMPLES::

            sage: (1).inverse_of_unit()
            1
            sage: (-1).inverse_of_unit()
            -1
            sage: 5.inverse_of_unit()
            Traceback (most recent call last):
            ...
            ArithmeticError: inverse does not exist
            sage: 0.inverse_of_unit()
            Traceback (most recent call last):
            ...
            ArithmeticError: inverse does not exist
        """
        if mpz_cmpabs_ui(self.value, 1) == 0:
            return self
        else:
            raise ArithmeticError("inverse does not exist")

    def inverse_mod(self, n):
        r"""
        Return the inverse of ``self`` modulo `n`, if this inverse exists.

        Otherwise, raise a :exc:`ZeroDivisionError` exception.

        INPUT:

        - ``self`` -- integer

        - ``n`` -- integer or ideal of integer ring

        OUTPUT:

        - ``x`` -- integer such that x * ``self`` = 1 (mod m), or
          raises :exc:`ZeroDivisionError`

        IMPLEMENTATION:

        Call the ``mpz_invert`` GMP library function.

        EXAMPLES::

            sage: a = Integer(189)
            sage: a.inverse_mod(10000)
            4709
            sage: a.inverse_mod(-10000)
            4709
            sage: a.inverse_mod(1890)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: inverse of Mod(189, 1890) does not exist
            sage: a = Integer(19)**100000  # long time
            sage: c = a.inverse_mod(a*a)   # long time
            Traceback (most recent call last):
            ...
            ZeroDivisionError: inverse of Mod(..., ...) does not exist

        We check that :issue:`10625` is fixed::

            sage: ZZ(2).inverse_mod(ZZ.ideal(3))
            2

        We check that :issue:`9955` is fixed::

            sage: Rational(3) % Rational(-1)
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
            raise ZeroDivisionError(f"inverse of Mod({self}, {m}) does not exist")
        return ans

    def crt(self, y, m, n):
        """
        Return the unique integer between `0` and `mn` that is congruent to
        the integer modulo `m` and to `y` modulo `n`.

        We assume that `m` and `n` are coprime.

        EXAMPLES::

            sage: n = 17
            sage: m = n.crt(5, 23, 11); m
            247
            sage: m%23
            17
            sage: m%11
            5
        """
        cdef object g, s
        cdef Integer _y, _m, _n
        _y = Integer(y)
        _m = Integer(m)
        _n = Integer(n)
        g, s, _ = _m.xgcd(_n)
        if not g.is_one():
            raise ArithmeticError("CRT requires that gcd of moduli is 1.")
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
        If ``self`` < 0, we return Infinity.

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
        if mpz_sgn(self.value) < 0:
            return sage.rings.infinity.Infinity
        return smallInteger(mpz_popcount(self.value))

    def conjugate(self):
        """
        Return the complex conjugate of this integer, which is the
        integer itself.

        EXAMPLES::

            sage: n = 205
            sage: n.conjugate()
            205
        """
        return self

    def binomial(self, m, algorithm='gmp'):
        """
        Return the binomial coefficient "``self`` choose ``m``".

        INPUT:

        - ``m`` -- integer

        - ``algorithm`` -- ``'gmp'`` (default), ``'mpir'`` (an alias for
          ``gmp``), or ``'pari'``; ``'gmp'`` is faster for small ``m``,
          and ``'pari'`` tends to be faster for large ``m``

        OUTPUT: integer

        EXAMPLES::

            sage: 10.binomial(2)
            45
            sage: 10.binomial(2, algorithm='pari')                                      # needs sage.libs.pari
            45
            sage: 10.binomial(-2)
            0
            sage: (-2).binomial(3)
            -4
            sage: (-3).binomial(0)
            1

        The argument ``m`` or (``self - m``) must fit into an ``unsigned long``::

            sage: (2**256).binomial(2**256)
            1
            sage: (2**256).binomial(2**256 - 1)
            115792089237316195423570985008687907853269984665640564039457584007913129639936
            sage: (2**256).binomial(2**128)
            Traceback (most recent call last):
            ...
            OverflowError: m must fit in an unsigned long

        TESTS::

            sage: 0.binomial(0)
            1
            sage: 0.binomial(1)
            0
            sage: 0.binomial(-1)
            0
            sage: 13.binomial(2r)
            78

        Check that it can be interrupted (:issue:`17852`)::

            sage: from sage.doctest.util import ensure_interruptible_after
            sage: with ensure_interruptible_after(0.5): (2^100).binomial(2^22, algorithm='mpir')

        For PARI, we try 10 interrupts with increasing intervals to
        check for reliable interrupting, see :issue:`18919`::

            sage: from cysignals import AlarmInterrupt
            sage: for i in [1..10]:             # long time (5s)                        # needs sage.libs.pari
            ....:     with ensure_interruptible_after(i/11):
            ....:         (2^100).binomial(2^22, algorithm='pari')
            doctest:...: RuntimeWarning: cypari2 leaked ... bytes on the PARI stack...
        """
        cdef Integer x
        cdef Integer mm

        if isinstance(m, Integer):
            mm = m
        else:
            mm = Integer(m)

        # trivial cases and potential simplification binom(n,x) -> binom(n,n-x)
        if self == zero or mm < zero or mm > self > zero:
            return one if mm == zero else zero

        if 2*mm > self > zero:
            mm = self - mm

        if mm == zero:
            return one
        if mm == one:
            return self

        # now call the various backend
        if algorithm == 'gmp' or algorithm == 'mpir':
            x = PY_NEW(Integer)
            if mpz_fits_ulong_p(mm.value):
                sig_on()
                mpz_bin_ui(x.value, self.value, mpz_get_ui(mm.value))
                sig_off()
            else:
                raise OverflowError("m must fit in an unsigned long")
            return x
        elif algorithm == 'pari':
            return the_integer_ring(self.__pari__().binomial(mm))
        else:
            raise ValueError("algorithm must be one of: 'pari' or 'gmp' (alias: 'mpir')")

    def to_bytes(self, length=1, byteorder='big', is_signed=False):
        r"""
        Return an array of bytes representing an integer.

        Internally relies on the python ``int.to_bytes()`` method.

        INPUT:

        - ``length`` -- positive integer (default: 1); integer represented
          in ``length`` bytes
        - ``byteorder`` -- string (default: ``'big'``); determines the byte
          order of the output (can only be ``'big'`` or ``'little'``)
        - ``is_signed`` -- boolean (default: ``False``); determines whether to use two's
          compliment to represent the integer

        .. TODO::

            It should be possible to convert straight from the gmp type in cython.
            This could be significantly faster, but I am unsure of the fastest and cleanest
            way to do this.

        EXAMPLES::

            sage: (1024).to_bytes(2, byteorder='big')
            b'\x04\x00'
            sage: (1024).to_bytes(10, byteorder='big')
            b'\x00\x00\x00\x00\x00\x00\x00\x00\x04\x00'
            sage: (-1024).to_bytes(10, byteorder='big', is_signed=True)
            b'\xff\xff\xff\xff\xff\xff\xff\xff\xfc\x00'
            sage: x = 1000
            sage: x.to_bytes((x.bit_length() + 7) // 8, byteorder='little')
            b'\xe8\x03'
        """
        return int(self).to_bytes(length=length, byteorder=byteorder, signed=is_signed)

cdef int mpz_set_str_python(mpz_ptr z, char* s, int base) except -1:
    """
    Wrapper around ``mpz_set_str()`` which supports :pep:`3127`
    literals.

    If the string is invalid, a :exc:`TypeError` will be raised.

    INPUT:

    - ``z`` -- a pre-allocated ``mpz_t`` where the result will be stored

    - ``s`` -- string to be converted to an ``mpz_t``

    - ``base`` -- either 0 or a base between 2 and 36: a base to use
      for the string conversion. 0 means auto-detect using prefixes

    EXAMPLES::

        sage: Integer('12345')
        12345
        sage: Integer('   -      1  2   3  4   5  ')
        -12345
        sage: Integer('  -  0x  1  2   3  4   5  ')
        -74565
        sage: Integer('-0012345', 16)
        -74565
        sage: Integer('+0x12345')
        74565
        sage: Integer('0X12345', 16)
        Traceback (most recent call last):
        ...
        TypeError: unable to convert '0X12345' to an integer
        sage: Integer('0x12345', 1000)
        Traceback (most recent call last):
        ...
        ValueError: base (=1000) must be 0 or between 2 and 36
        sage: Integer('0x00DeadBeef')
        3735928559
        sage: Integer('0x0x12345')
        Traceback (most recent call last):
        ...
        TypeError: unable to convert '0x0x12345' to an integer
        sage: Integer('-0B100')
        -4
        sage: Integer('-0B100', 16)
        -45312
        sage: Integer('0B12345')
        Traceback (most recent call last):
        ...
        TypeError: unable to convert '0B12345' to an integer

    Test zeros::

        sage: Integer('')
        Traceback (most recent call last):
        ...
        TypeError: unable to convert '' to an integer
        sage: Integer("0")
        0
        sage: Integer("  0O  0  ")  # second character is the letter O
        0
        sage: Integer("-00")
        0
        sage: Integer("+00000", 4)
        0

    For octals, the old leading-zero style is no longer available (unless an
    explicit base is given)::

        sage: Integer('0o12')
        10
        sage: Integer('012', 8)
        10
        sage: Integer('012')
        12

    We disallow signs in unexpected places::

        sage: Integer('+ -0')
        Traceback (most recent call last):
        ...
        TypeError: unable to convert '+ -0' to an integer
        sage: Integer('0o-0')
        Traceback (most recent call last):
        ...
        TypeError: unable to convert '0o-0' to an integer
    """
    cdef int sign = 1
    cdef char* x = s

    if base != 0 and (base < 2 or base > 36):
        raise ValueError("base (=%s) must be 0 or between 2 and 36" % base)

    while x[0] == c' ':
        x += 1  # Strip spaces

    # Check for signs
    if x[0] == c'-':
        sign = -1
        x += 1
    elif x[0] == c'+':
        x += 1

    while x[0] == c' ':
        x += 1  # Strip spaces

    # If no base was given, check for PEP 3127 prefixes
    if base == 0:
        if x[0] != c'0':
            base = 10
        else:
            # String starts with "0"
            if x[1] == c'b' or x[1] == c'B':
                x += 2
                base = 2
            elif x[1] == c'o' or x[1] == c'O':
                x += 2
                base = 8
            elif x[1] == c'x' or x[1] == c'X':
                x += 2
                base = 16
            else:
                base = 10  # no longer giving an octal

    while x[0] == c' ':
        x += 1  # Strip spaces

    # Disallow a sign here
    if x[0] == c'-' or x[0] == c'+':
        x = ""  # Force an error below

    assert base >= 2
    if mpz_set_str(z, x, base) != 0:
        raise TypeError("unable to convert %r to an integer" % char_to_str(s))
    if sign < 0:
        mpz_neg(z, z)


def GCD_list(v):
    r"""
    Return the greatest common divisor of a list of integers.

    INPUT:

    - ``v`` -- list or tuple

    Elements of `v` are converted to Sage integers.  An empty list has
    GCD zero.

    This function is used, for example, by ``rings/arith.py``.

    EXAMPLES::

        sage: from sage.rings.integer import GCD_list
        sage: w = GCD_list([3,9,30]); w
        3
        sage: type(w)
        <class 'sage.rings.integer.Integer'>

    Check that the bug reported in :issue:`3118` has been fixed::

        sage: sage.rings.integer.GCD_list([2,2,3])
        1

    The inputs are converted to Sage integers.

    ::

        sage: w = GCD_list([int(3), int(9), '30']); w
        3
        sage: type(w)
        <class 'sage.rings.integer.Integer'>

    Check that the GCD of the empty list is zero (:issue:`17257`)::

        sage: GCD_list([])
        0
    """
    cdef int i, n = len(v)
    cdef Integer z = <Integer>PY_NEW(Integer)

    for i in range(n):
        if not isinstance(v[i], Integer):
            if not isinstance(v, list):
                v = list(v)
            v[i] = Integer(v[i])

    if n == 0:
        return zero
    elif n == 1:
        return v[0].abs()

    sig_on()
    mpz_gcd(z.value, (<Integer>v[0]).value, (<Integer>v[1]).value)
    for i in range(2, n):
        if mpz_cmp_ui(z.value, 1) == 0:
            break
        mpz_gcd(z.value, z.value, (<Integer>v[i]).value)
    sig_off()

    return z


@cython.binding(True)
def make_integer(s):
    """
    Create a Sage integer from the base-32 Python *string* ``s``. This is
    used in unpickling integers.

    EXAMPLES::

        sage: from sage.rings.integer import make_integer
        sage: make_integer('-29')
        -73
        sage: make_integer(29)
        Traceback (most recent call last):
        ...
        TypeError: expected str...Integer found
    """
    cdef Integer r = PY_NEW(Integer)
    mpz_set_str(r.value, str_to_bytes(s), 32)
    return r


cdef class int_to_Z(Morphism):
    """
    EXAMPLES::

        sage: f = ZZ.coerce_map_from(int)
        sage: f
        Native morphism:
          From: Set of Python objects of class 'int'
          To:   Integer Ring
        sage: f(1rL)
        1
    """
    def __init__(self):
        import sage.categories.homset
        from sage.sets.pythonclass import Set_PythonType
        Morphism.__init__(self, sage.categories.homset.Hom(Set_PythonType(int), integer_ring.ZZ))

    cpdef Element _call_(self, a):
        cdef Integer r
        cdef long l
        cdef int err = 0

        integer_check_long_py(a, &l, &err)
        if not err:
            return smallInteger(l)

        r = <Integer>PY_NEW(Integer)
        mpz_set_pylong(r.value, a)
        return r

    def _repr_type(self):
        return "Native"

# ############## INTEGER CREATION CODE #####################

# This variable holds the size of any Integer object in bytes.
cdef int sizeof_Integer

# We use a global Integer element to steal all the references
# from. DO NOT INITIALIZE IT AGAIN and DO NOT REFERENCE IT!
cdef Integer global_dummy_Integer
global_dummy_Integer = Integer()
# Reallocate to one limb to fix :issue:`31340` and :issue:`33081`
_mpz_realloc(global_dummy_Integer.value, 1)


def _check_global_dummy_Integer():
    """
    Return ``True`` if the global dummy :class:`Integer` is ok.

    TESTS::

        sage: from sage.rings.integer import _check_global_dummy_Integer
        sage: _check_global_dummy_Integer()
        True
    """
    # Check that it has exactly one limb allocated
    # This is assumed later in fast_tp_new() :issue:`33081`
    cdef mpz_ptr dummy = global_dummy_Integer.value
    if dummy._mp_alloc == 1 and dummy._mp_size == 0:
        return True

    raise AssertionError(
      "global dummy Integer is corrupt (_mp_alloc = %d, _mp_size = %d)"
      % (dummy._mp_alloc, dummy._mp_size))


# A global pool for performance when integers are rapidly created and destroyed.
# It operates on the following principles:
#
# - The pool starts out empty.
# - When a new integer is needed, one from the pool is returned
#   if available, otherwise a new Integer object is created
# - When an integer is collected, it will add it to the pool
#   if there is room, otherwise it will be deallocated.
cdef int integer_pool_size = 100

cdef PyObject** integer_pool
cdef int integer_pool_count = 0

# used for profiling the pool
cdef int total_alloc = 0
cdef int use_pool = 0


cdef PyObject* fast_tp_new(type t, args, kwds) except NULL:
    global integer_pool, integer_pool_count, total_alloc, use_pool

    cdef PyObject* new
    cdef mpz_ptr new_mpz

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
        # sizeof(Integer). The use of PyObject_Malloc directly
        # assumes that Integers are not garbage collected, i.e.
        # they do not possess references to other Python
        # objects (as indicated by the Py_TPFLAGS_HAVE_GC flag).
        # See below for a more detailed description.
        new = <PyObject*>PyObject_Malloc(sizeof_Integer)
        if unlikely(new == NULL):
            raise MemoryError

        # Now set every member as set in z, the global dummy Integer
        # created before this tp_new started to operate.
        memcpy(new, (<void*>global_dummy_Integer), sizeof_Integer)

        # We allocate memory for the _mp_d element of the value of this
        # new Integer. We allocate one limb. Normally, one would use
        # mpz_init() for this, but we allocate the memory directly.
        # This saves time both by avoiding extra function calls and
        # because the rest of the mpz struct was already initialized
        # fully using the memcpy above.
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
        # NOTE: This assumes global_dummy_Integer.value._mp_alloc == 1
        new_mpz = <mpz_ptr>((<Integer>new).value)
        new_mpz._mp_d = <mp_ptr>check_malloc(GMP_LIMB_BITS >> 3)

    # This line is only needed if Python is compiled in debugging mode
    # './configure --with-pydebug' or SAGE_DEBUG=yes. If that is the
    # case a Python object has a bunch of debugging fields which are
    # initialized with this macro.

    if_Py_TRACE_REFS_then_PyObject_INIT(
            new, Py_TYPE(global_dummy_Integer))

    # The global_dummy_Integer may have a reference count larger than
    # one, but it is expected that newly created objects have a
    # reference count of one. This is potentially unneeded if
    # everybody plays nice, because the gobal_dummy_Integer has only
    # one reference in that case.

    # Objects from the pool have reference count zero, so this
    # needs to be set in this case.

    new.ob_refcnt = 1

    return new


cdef void fast_tp_dealloc(PyObject* o) noexcept:
    # If there is room in the pool for a used integer object,
    # then put it in rather than deallocating it.
    global integer_pool, integer_pool_count

    cdef mpz_ptr o_mpz = <mpz_ptr>((<Integer>o).value)

    # If we are recovering from an interrupt, throw the mpz_t away
    # without recycling or freeing it because it might be in an
    # inconsistent state (see Issue #24986).
    if sig_occurred() is NULL:
        if integer_pool_count < integer_pool_size:
            # Here we free any extra memory used by the mpz_t by
            # setting it to a single limb.
            if o_mpz._mp_alloc > 10:
                _mpz_realloc(o_mpz, 1)

            # It's cheap to zero out an integer, so do it here.
            o_mpz._mp_size = 0

            # And add it to the pool.
            integer_pool[integer_pool_count] = o
            integer_pool_count += 1
            return

        # No space in the pool, so just free the mpz_t.
        sig_free(o_mpz._mp_d)

    # Free the object. This assumes that Py_TPFLAGS_HAVE_GC is not
    # set. If it was set another free function would need to be
    # called.
    PyObject_Free(o)


from sage.misc.allocator cimport hook_tp_functions
cdef hook_fast_tp_functions():
    """
    Initialize the fast integer creation functions.
    """
    global global_dummy_Integer, sizeof_Integer, integer_pool

    integer_pool = <PyObject**>check_allocarray(integer_pool_size, sizeof(PyObject*))

    cdef PyObject* o
    o = <PyObject *>global_dummy_Integer

    # store how much memory needs to be allocated for an Integer.
    sizeof_Integer = o.ob_type.tp_basicsize

    # Finally replace the functions called when an Integer needs
    # to be constructed/destructed.
    hook_tp_functions(global_dummy_Integer, <newfunc>(&fast_tp_new), <destructor>(&fast_tp_dealloc), False)

cdef integer(x):
    if isinstance(x, Integer):
        return x
    return Integer(x)


def free_integer_pool():
    cdef int i
    cdef PyObject *o

    global integer_pool_count, integer_pool_size

    for i in range(integer_pool_count):
        o = integer_pool[i]
        mpz_clear((<Integer>o).value)
        # Free the object. This assumes that Py_TPFLAGS_HAVE_GC is not
        # set. If it was set another free function would need to be
        # called.
        PyObject_Free(o)

    integer_pool_size = 0
    integer_pool_count = 0
    sig_free(integer_pool)


# Replace default allocation and deletion with faster custom ones
hook_fast_tp_functions()

# zero and one initialization
initialized = False
cdef set_zero_one_elements():
    global the_integer_ring, initialized
    if initialized:
        return
    the_integer_ring._zero_element = Integer(0)
    the_integer_ring._one_element = Integer(1)
    initialized = True
set_zero_one_elements()

cdef Integer zero = the_integer_ring._zero_element
cdef Integer one = the_integer_ring._one_element

# pool of small integer for fast sign computation
# Use the same defaults as Python3 documented at
# https://docs.python.org/3/c-api/long.html#c.PyLong_FromLong
DEF small_pool_min = -5
DEF small_pool_max = 256
# we could use the above zero and one here
cdef list small_pool = [Integer(k) for k in range(small_pool_min, small_pool_max+1)]

cdef inline Integer smallInteger(long value):
    """
    This is the fastest way to create a (likely) small Integer.
    """
    cdef Integer z
    if small_pool_min <= value <= small_pool_max:
        return <Integer>small_pool[value - small_pool_min]
    else:
        z = PY_NEW(Integer)
        mpz_set_si(z.value, value)
        return z


# The except value is just some random double, it doesn't matter what it is.
cdef double mpz_get_d_nearest(mpz_t x) except? -648555075988944.5:
    """
    Convert a ``mpz_t`` to a ``double``, with round-to-nearest-even.
    This differs from ``mpz_get_d()`` which does round-to-zero.

    TESTS::

        sage: x = ZZ(); float(x)
        0.0
        sage: x = 2^54 - 1
        sage: float(x)
        1.8014398509481984e+16
        sage: float(-x)
        -1.8014398509481984e+16
        sage: x = 2^10000; float(x)
        inf
        sage: float(-x)
        -inf

    ::

        sage: x = (2^53 - 1) * 2^971; float(x)  # Largest double
        1.7976931348623157e+308
        sage: float(-x)
        -1.7976931348623157e+308
        sage: x = (2^53) * 2^971; float(x)
        inf
        sage: float(-x)
        -inf
        sage: x = ZZ((2^53 - 1/2) * 2^971); float(x)
        inf
        sage: float(-x)
        -inf
        sage: x = ZZ((2^53 - 3/4) * 2^971); float(x)
        1.7976931348623157e+308
        sage: float(-x)
        -1.7976931348623157e+308

    AUTHORS:

    - Jeroen Demeyer (:issue:`16385`, based on :issue:`14416`)
    """
    cdef mp_bitcnt_t sx = mpz_sizeinbase(x, 2)

    # Easy case: x is exactly representable as double.
    if sx <= 53:
        return mpz_get_d(x)

    cdef int resultsign = mpz_sgn(x)

    # Check for overflow
    if sx > 1024:
        if resultsign < 0:
            return float('-inf')
        else:
            return float('inf')

    # General case

    # We should shift x right by this amount in order
    # to have 54 bits remaining.
    cdef mp_bitcnt_t shift = sx - 54

    # Compute q = trunc(x / 2^shift) and let remainder_is_zero be True
    # if and only if no truncation occurred.
    cdef int remainder_is_zero
    remainder_is_zero = mpz_divisible_2exp_p(x, shift)

    sig_on()

    cdef mpz_t q
    mpz_init(q)
    mpz_tdiv_q_2exp(q, x, shift)

    # Convert abs(q) to a 64-bit integer.
    cdef mp_limb_t* q_limbs = (<mpz_ptr>q)._mp_d
    cdef uint64_t q64
    if sizeof(mp_limb_t) >= 8:
        q64 = q_limbs[0]
    else:
        assert sizeof(mp_limb_t) == 4
        q64 = q_limbs[1]
        q64 = (q64 << 32) + q_limbs[0]

    mpz_clear(q)
    sig_off()

    # Round q from 54 to 53 bits of precision.
    if ((q64 & 1) == 0):
        # Round towards zero
        pass
    else:
        if not remainder_is_zero:
            # Remainder is nonzero: round away from zero
            q64 += 1
        else:
            # Halfway case: round to even
            q64 += (q64 & 2) - 1

    # The conversion of q64 to double is *exact*.
    # This is because q64 is even and satisfies 2^53 <= q64 <= 2^54.
    cdef double d = <double>q64
    if resultsign < 0:
        d = -d
    return ldexp(d, shift)


# Support Python's numbers abstract base class
import numbers
numbers.Integral.register(Integer)

# Free the memory used by the integer pool when sage exits. This is
# not strictly necessary because the OS should immediately reclaim
# these resources when sage terminates. However, it may aid valgrind
# or similar tools, and can help expose bugs in other code.
import atexit
atexit.register(free_integer_pool)
