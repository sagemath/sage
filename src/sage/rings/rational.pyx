# distutils: libraries = NTL_LIBRARIES
# distutils: extra_compile_args = NTL_CFLAGS
# distutils: include_dirs = NTL_INCDIR
# distutils: library_dirs = NTL_LIBDIR
# distutils: extra_link_args = NTL_LIBEXTRA
# distutils: language = c++
r"""
Rational Numbers

AUTHORS:

- William Stein (2005): first version

- William Stein (2006-02-22): floor and ceil (pure fast GMP versions).

- Gonzalo Tornaria and William Stein (2006-03-02): greatly improved
  python/GMP conversion; hashing

- William Stein and Naqi Jaffery (2006-03-06): height, sqrt examples,
  and improve behavior of sqrt.

- David Harvey (2006-09-15): added nth_root

- Pablo De Napoli (2007-04-01): corrected the implementations of
  multiplicative_order, is_one; optimized __bool__ ; documented:
  lcm,gcd

- John Cremona (2009-05-15): added support for local and global
  logarithmic heights.

- Travis Scrimshaw (2012-10-18): Added doctests for full coverage.

- Vincent Delecroix (2013): continued fraction

- Vincent Delecroix (2017-05-03): faster integer-rational comparison

- Vincent Klein (2017-05-11): add __mpq__() to class Rational

- Vincent Klein (2017-05-22): Rational constructor support gmpy2.mpq
  or gmpy2.mpz parameter. Add __mpz__ to class Rational.

TESTS::

    sage: a = -2/3
    sage: a == loads(dumps(a))
    True
"""

# ****************************************************************************
#       Copyright (C) 2004, 2006 William Stein <wstein@gmail.com>
#       Copyright (C) 2017 Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

cimport cython
from cpython cimport *
from cpython.object cimport Py_EQ, Py_NE

from cysignals.signals cimport sig_on, sig_off

import operator
import fractions

from sage.rings.rational_field import QQ

from sage.arith.long cimport integer_check_long_py
from sage.categories.morphism cimport Morphism
from sage.categories.map cimport Map
from sage.cpython.string cimport char_to_str, str_to_bytes
from sage.libs.gmp.pylong cimport mpz_set_pylong
from sage.rings.integer cimport Integer, smallInteger
from sage.structure.coerce cimport coercion_model, is_numpy_type
from sage.structure.element cimport Element
from sage.structure.parent cimport Parent
from sage.structure.richcmp cimport rich_to_bool_sgn


RealNumber_classes = ()


def _register_real_number_class(cls):
    r"""
    Register ``cls``.

    This is called by ``sage.rings.real_mpfr``, to avoid a cyclic import.
    """
    global RealNumber_classes
    RealNumber_classes += (cls,)


RealDouble_classes = (float,)
try:
    from sage.rings.real_double import RealDoubleElement
    RealDouble_classes += (RealDoubleElement,)
except ImportError:
    pass


from libc.stdint cimport uint64_t
from sage.libs.gmp.binop cimport mpq_add_z, mpq_mul_z, mpq_div_zz

cimport sage.rings.fast_arith
import sage.rings.fast_arith


try:
    from cypari2.gen import Gen as pari_gen
except ImportError:
    pari_gen = ()


set_rational_from_gen = None
new_gen_from_rational = None


cdef sage.rings.fast_arith.arith_int ai
ai = sage.rings.fast_arith.arith_int()

cdef object numpy_long_interface = {'typestr': '=i4' if sizeof(long) == 4 else '=i8' }
cdef object numpy_int64_interface = {'typestr': '=i8'}
cdef object numpy_object_interface = {'typestr': '|O'}
cdef object numpy_double_interface = {'typestr': '=f8'}

from libc.math cimport ldexp
from sage.libs.gmp.all cimport *

cimport gmpy2
gmpy2.import_gmpy2()


cdef class Rational(sage.structure.element.FieldElement)

cdef inline void set_from_mpq(Rational self, mpq_t value) noexcept:
    mpq_set(self.value, value)

cdef inline void set_from_Rational(Rational self, Rational other) noexcept:
    mpq_set(self.value, other.value)

cdef inline void set_from_Integer(Rational self, integer.Integer other) noexcept:
    mpq_set_z(self.value, other.value)

cdef object Rational_mul_(Rational a, Rational b):
    cdef Rational x
    x = <Rational> Rational.__new__(Rational)

    sig_on()
    mpq_mul(x.value, a.value, b.value)
    sig_off()

    return x

cdef object Rational_div_(Rational a, Rational b):
    cdef Rational x
    x = <Rational> Rational.__new__(Rational)

    sig_on()
    mpq_div(x.value, a.value, b.value)
    sig_off()

    return x

cdef Rational_add_(Rational self, Rational other):
    cdef Rational x
    x = <Rational> Rational.__new__(Rational)
    sig_on()
    mpq_add(x.value, self.value, other.value)
    sig_off()
    return x

cdef Rational_sub_(Rational self, Rational other):
    cdef Rational x
    x = <Rational> Rational.__new__(Rational)

    sig_on()
    mpq_sub(x.value, self.value, other.value)
    sig_off()

    return x

cpdef Integer integer_rational_power(Integer a, Rational b):
    """
    Compute `a^b` as an integer, if it is integral, or return ``None``.

    The nonnegative real root is taken for even denominators.

    INPUT:

    - ``a`` -- an ``Integer``
    - ``b`` -- a nonnegative ``Rational``

    OUTPUT: `a^b` as an ``Integer`` or ``None``

    EXAMPLES::

        sage: from sage.rings.rational import integer_rational_power
        sage: integer_rational_power(49, 1/2)
        7
        sage: integer_rational_power(27, 1/3)
        3
        sage: integer_rational_power(-27, 1/3) is None
        True
        sage: integer_rational_power(-27, 2/3) is None
        True
        sage: integer_rational_power(512, 7/9)
        128

        sage: integer_rational_power(27, 1/4) is None
        True
        sage: integer_rational_power(-16, 1/4) is None
        True

        sage: integer_rational_power(0, 7/9)
        0
        sage: integer_rational_power(1, 7/9)
        1
        sage: integer_rational_power(-1, 7/9) is None
        True
        sage: integer_rational_power(-1, 8/9) is None
        True
        sage: integer_rational_power(-1, 9/8) is None
        True

    TESTS (:issue:`11228`)::

        sage: integer_rational_power(-10, QQ(2))
        100
        sage: integer_rational_power(0, QQ(0))
        1
    """
    cdef Integer z = Integer.__new__(Integer)
    if mpz_sgn(mpq_numref(b.value)) < 0:
        raise ValueError("Only positive exponents supported.")
    cdef int sgn = mpz_sgn(a.value)
    cdef bint exact
    if (mpz_cmp_ui(a.value, 1) == 0 or
          mpz_cmp_ui(mpq_numref(b.value), 0) == 0):
        mpz_set_ui(z.value, 1)
    elif sgn == 0:
        pass # z is 0
    elif sgn < 0 and mpz_cmp_ui(mpq_denref(b.value), 1):
        return None
    else:
        if (not mpz_fits_ulong_p(mpq_numref(b.value))
            or not mpz_fits_ulong_p(mpq_denref(b.value))):
            # too big to take roots/powers
            return None
        elif mpz_cmp_ui(mpq_denref(b.value), 2) == 0:
            if mpz_perfect_square_p(a.value):
                mpz_sqrt(z.value, a.value)
            else:
                return None
        else:
            exact = mpz_root(z.value, a.value, mpz_get_ui(mpq_denref(b.value)))
            if not exact:
                return None
        mpz_pow_ui(z.value, z.value, mpz_get_ui(mpq_numref(b.value)))
    return z


cpdef rational_power_parts(a, Rational b, factor_limit=10**5):
    """
    Compute rationals or integers `c` and `d` such that `a^b = c*d^b`
    with `d` small. This is used for simplifying radicals.

    INPUT:

    - ``a`` -- a rational or integer
    - ``b`` -- a rational
    - ``factor_limit`` -- the limit used in factoring ``a``

    EXAMPLES::

        sage: from sage.rings.rational import rational_power_parts
        sage: rational_power_parts(27, 1/2)
        (3, 3)
        sage: rational_power_parts(-128, 3/4)
        (8, -8)
        sage: rational_power_parts(-4, 1/2)
        (2, -1)
        sage: rational_power_parts(-4, 1/3)
        (1, -4)
        sage: rational_power_parts(9/1000, 1/2)
        (3/10, 1/10)

    TESTS:

    Check if :issue:`8540` is fixed::

        sage: rational_power_parts(3/4, -1/2)
        (2, 3)
        sage: t = (3/4)^(-1/2); t                                                       # needs sage.symbolic
        2/3*sqrt(3)
        sage: t^2                                                                       # needs sage.symbolic
        4/3

    Check if :issue:`15605` is fixed::

        sage: rational_power_parts(-1, -1/3)
        (1, -1)
        sage: rational_power_parts(-1, 2/3)
        (1, -1)
        sage: all(rational_power_parts(-1, i/77) == (1,-1) for i in range(1,9))
        True

        sage: # needs sage.symbolic
        sage: (-1)^(-1/3)
        -(-1)^(2/3)
        sage: 1 / ((-1)^(1/3))
        -(-1)^(2/3)
        sage: (-1)^(2/3)
        (-1)^(2/3)
        sage: (-1)^(1/3)*(-1)^(1/5)
        (-1)^(8/15)
        sage: bool((-1)^(2/3) == -1/2 + sqrt(3)/2*I)
        True
        sage: all((-1)^(p/q) == cos(p*pi/q) + I * sin(p*pi/q)
        ....:     for p in srange(1, 6) for q in srange(1, 6))
        True

    A few more tests added in :issue:`26414`::

        sage: rational_power_parts(-1, 2/1)
        (1, 1)
        sage: rational_power_parts(-8, 2/3)
        (4, -1)
        sage: all(isinstance(z, Integer) for z in rational_power_parts(-1, 1/1))
        True
        sage: all(isinstance(z, Integer) for z in rational_power_parts(-1, 2/3))
        True
    """
    cdef bint b_negative = (b < 0)
    if b_negative:
        b = -b
        a = ~a

    if isinstance(a, Integer):
        pass
    elif isinstance(a, Rational):
        c1, d1 = rational_power_parts(a.numerator(), b)
        c2, d2 = rational_power_parts(a.denominator(), b)
        return (c1/c2, d1/d2) if not b_negative else (c1/c2, d2/d1)
    else:
        a = Integer(a)

    c = integer_rational_power(a, b)
    if c is not None:
        return c, smallInteger(1)

    numer, denom = b.numerator(), b.denominator()
    if a < factor_limit*factor_limit:
        f = a.factor()
    else:
        from sage.rings.factorint import factor_trial_division
        f = factor_trial_division(a, factor_limit)
    c = smallInteger(1)
    # The sign is not handled by the loop below. We don't want to
    # simplify (-1)^(2/3) to 1 (see Issue #15605), so we always move
    # the sign over to d. Note that the case (-1)^2 is already
    # handled by integer_rational_power() above.
    if a >= 0:
        # d = 1
        d = c
    else:
        # d = -1
        d = smallInteger(-1)
    for p, e in f:
        c *= p**((e // denom)*numer)
        d *= p**(e % denom)
    return (c, d) if not b_negative else (c, ~d)


def is_Rational(x):
    """
    Return ``True`` if ``x`` is of the Sage :class:`Rational` type.

    EXAMPLES::

        sage: from sage.rings.rational import is_Rational
        sage: is_Rational(2)
        doctest:warning...
        DeprecationWarning: The function is_Rational is deprecated;
        use 'isinstance(..., Rational)' instead.
        See https://github.com/sagemath/sage/issues/38128 for details.
        False
        sage: is_Rational(2/1)
        True
        sage: is_Rational(int(2))
        False
        sage: is_Rational('5')
        False
    """
    from sage.misc.superseded import deprecation_cython
    deprecation_cython(38128,
                       "The function is_Rational is deprecated; "
                       "use 'isinstance(..., Rational)' instead.")
    return isinstance(x, Rational)


cdef class Rational(sage.structure.element.FieldElement):
    """
    A rational number.

    Rational numbers are implemented using the GMP C library.

    EXAMPLES::

        sage: a = -2/3
        sage: type(a)
        <class 'sage.rings.rational.Rational'>
        sage: parent(a)
        Rational Field
        sage: Rational('1/0')
        Traceback (most recent call last):
        ...
        TypeError: unable to convert '1/0' to a rational
        sage: Rational(1.5)
        3/2
        sage: Rational('9/6')
        3/2
        sage: Rational((2^99,2^100))
        1/2
        sage: Rational(("2", "10"), 16)
        1/8
        sage: Rational(QQbar(125/8).nth_root(3))                                        # needs sage.rings.number_field
        5/2
        sage: Rational(AA(209735/343 - 17910/49*golden_ratio).nth_root(3)               # needs sage.rings.number_field sage.symbolic
        ....:          + 3*AA(golden_ratio))
        53/7
        sage: QQ(float(1.5))
        3/2
        sage: QQ(RDF(1.2))
        6/5

    Conversion from fractions::

        sage: import fractions
        sage: f = fractions.Fraction(1r, 2r)
        sage: Rational(f)
        1/2

    Conversion from PARI::

        sage: Rational(pari('-939082/3992923'))                                         # needs sage.libs.pari
        -939082/3992923
        sage: Rational(pari('Pol([-1/2])'))  #9595                                      # needs sage.libs.pari
        -1/2

    Conversions from numpy::

        sage: # needs numpy
        sage: import numpy as np
        sage: QQ(np.int8('-15'))
        -15
        sage: QQ(np.int16('-32'))
        -32
        sage: QQ(np.int32('-19'))
        -19
        sage: QQ(np.uint32('1412'))
        1412

        sage: QQ(np.float16('12'))                                                      # needs numpy
        12

    Conversions from gmpy2::

        sage: from gmpy2 import *
        sage: QQ(mpq('3/4'))
        3/4
        sage: QQ(mpz(42))
        42
        sage: Rational(mpq(2/3))
        2/3
        sage: Rational(mpz(5))
        5

    TESTS:

    Check that :issue:`28321` is fixed::

        sage: QQ((2r^100r, 3r^100r))
        1267650600228229401496703205376/515377520732011331036461129765621272702107522001
        sage: QQ((-2r^100r, -3r^100r))
        1267650600228229401496703205376/515377520732011331036461129765621272702107522001
    """
    def __cinit__(self):
        r"""
        Initialize ``self`` as an element of `\QQ`.

        EXAMPLES::

            sage: p = Rational(3) # indirect doctest
            sage: p.parent()
            Rational Field
        """
        mpq_init(self.value)
        self._parent = QQ

    def __init__(self, x=None, unsigned int base=0):
        """
        Create a new rational number.

        INPUT:

        - ``x`` -- object (default: ``None``)

        - ``base`` -- base if ``x`` is a string

        EXAMPLES::

            sage: a = Rational()
            sage: a.__init__(7); a
            7
            sage: a.__init__('70', base=8); a
            56
            sage: a.__init__(pari('2/3')); a                                            # needs sage.libs.pari
            2/3
            sage: a.__init__('-h/3ki', 32); a
            -17/3730
            sage: from gmpy2 import mpq
            sage: a.__init__(mpq('3/5')); a
            3/5

        TESTS:

        Check that :issue:`19835` is fixed::

            sage: QQ((0r,-1r))
            0
            sage: QQ((-1r,-1r))
            1

        .. NOTE::

           This is for doctesting purposes only.  Rationals are defined
           to be immutable.
        """
        if x is not None:
            self.__set_value(x, base)

    def __reduce__(self):
        """
        Used in pickling rational numbers.

        EXAMPLES::

            sage: a = 3/5
            sage: a.__reduce__()
            (<cyfunction make_rational at ...>, ('3/5',))
        """
        return sage.rings.rational.make_rational, (self.str(32),)

    def __index__(self):
        """
        Needed so integers can be used as list indices.

        EXAMPLES::

            sage: v = [1,2,3,4,5]
            sage: v[3/1]
            4
            sage: v[3/2]
            Traceback (most recent call last):
            ...
            TypeError: unable to convert rational 3/2 to an integer
        """
        if self.denominator() == 1:
            return int(self)

        raise TypeError(f"unable to convert rational {self} to an integer")

    cdef __set_value(self, x, unsigned int base):
        cdef int n
        cdef Rational temp_rational
        cdef integer.Integer a, b

        if isinstance(x, Rational):
            set_from_Rational(self, x)

        elif isinstance(x, int):
            mpz_set_pylong(mpq_numref(self.value), x)

        elif isinstance(x, integer.Integer):
            set_from_Integer(self, x)

        elif isinstance(x, sage.rings.real_mpfr.RealNumber):

            if x == 0:
                mpq_set_si(self.value, 0, 1)
                return
            if not base:
                set_from_Rational(self, x.simplest_rational())
            else:
                # Truncate in base 10 to match repr(x).
                # See https://github.com/sagemath/sage/issues/21124
                xstr = x.str(base, truncate=(base == 10))
                if '.' in xstr:
                    exp = (len(xstr) - (xstr.index('.') +1))
                    p = base**exp
                    pstr = '1'+'0'*exp
                    s = xstr.replace('.','') +'/'+pstr
                    n = mpq_set_str(self.value, str_to_bytes(s), base)
                    if n or mpz_cmp_si(mpq_denref(self.value), 0) == 0:
                        raise TypeError("unable to convert {!r} to a rational".format(x))
                    mpq_canonicalize(self.value)
                else:
                    n = mpq_set_str(self.value, xstr, base)
                    if n or mpz_cmp_si(mpq_denref(self.value), 0) == 0:
                        raise TypeError("unable to convert {!r} to a rational".format(x))
                    mpq_canonicalize(self.value)
        elif isinstance(x, bytes):
            n = mpq_set_str(self.value, x, base)
            if n or mpz_cmp_si(mpq_denref(self.value), 0) == 0:
                raise TypeError("unable to convert {!r} to a rational".format(x))
            mpq_canonicalize(self.value)
        elif isinstance(x, unicode):
            n = mpq_set_str(self.value, str_to_bytes(x), base)
            if n or mpz_cmp_si(mpq_denref(self.value), 0) == 0:
                raise TypeError("unable to convert {!r} to a rational".format(x))
            mpq_canonicalize(self.value)

        elif hasattr(x, "_rational_"):
            set_from_Rational(self, x._rational_())

        elif isinstance(x, tuple) and len(x) == 2:
            num = x[0]
            denom = x[1]

            if isinstance(num, int):
                mpz_set_pylong(mpq_numref(self.value), num)
            else:
                if not isinstance(num, integer.Integer):
                    num = integer.Integer(num, base)
                mpz_set(mpq_numref(self.value), (<integer.Integer>num).value)

            if isinstance(denom, int):
                mpz_set_pylong(mpq_denref(self.value), denom)
            else:
                if not isinstance(denom, integer.Integer):
                    denom = integer.Integer(denom, base)
                mpz_set(mpq_denref(self.value), (<integer.Integer>denom).value)

            if mpz_sgn(mpq_denref(self.value)) == 0:
                raise ValueError("denominator must not be 0")

            mpq_canonicalize(self.value)

        elif isinstance(x, pari_gen):
            global set_rational_from_gen
            if set_rational_from_gen is None:
                from sage.libs.pari.convert_sage import set_rational_from_gen
            set_rational_from_gen(self, x)

        elif isinstance(x, list) and len(x) == 1:
            self.__set_value(x[0], base)

        elif hasattr(x, 'rational_reconstruction'):
            temp_rational = x.rational_reconstruction()
            mpq_set(self.value, temp_rational.value)

        elif isinstance(x, (float, sage.rings.real_double.RealDoubleElement)):
            self.__set_value(sage.rings.real_mpfr.RealNumber(sage.rings.real_mpfr.RR, x), base)

        elif is_numpy_type(type(x)):
            import numpy
            if isinstance(x, numpy.integer):
                self.__set_value(integer.Integer(x), base)
            elif isinstance(x, numpy.floating):
                self.__set_value(sage.rings.real_mpfr.RR(x), base)
            else:
                raise TypeError("unable to convert {!r} to a rational".format(x))

        elif isinstance(x, fractions.Fraction):
            mpz_set(mpq_numref(self.value), (<integer.Integer> integer.Integer(x.numerator)).value)
            mpz_set(mpq_denref(self.value), (<integer.Integer> integer.Integer(x.denominator)).value)

        elif type(x) is gmpy2.mpq:
            mpq_set(self.value, (<gmpy2.mpq>x).q)

        elif type(x) is gmpy2.mpz:
            mpq_set_z(self.value, (<gmpy2.mpz>x).z)

        else:
            raise TypeError("unable to convert {!r} to a rational".format(x))

    cdef void set_from_mpq(Rational self, mpq_t value) noexcept:
        mpq_set(self.value, value)

    def list(self):
        """
        Return a list with the rational element in it, to be compatible
        with the method for number fields.

        OUTPUT: the list ``[self]``

        EXAMPLES::

            sage: m = 5/3
            sage: m.list()
            [5/3]
        """
        return [ self ]

    def continued_fraction_list(self, type='std'):
        r"""
        Return the list of partial quotients of this rational number.

        INPUT:

        - ``type`` -- either ``'std'`` (the default) for the standard continued
          fractions or ``'hj'`` for the Hirzebruch-Jung ones

        EXAMPLES::

            sage: (13/9).continued_fraction_list()
            [1, 2, 4]
            sage: 1 + 1/(2 + 1/4)
            13/9

            sage: (225/157).continued_fraction_list()
            [1, 2, 3, 4,  5]
            sage: 1 + 1/(2 + 1/(3 + 1/(4 + 1/5)))
            225/157

            sage: (fibonacci(20)/fibonacci(19)).continued_fraction_list()               # needs sage.libs.pari
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2]

            sage: (-1/3).continued_fraction_list()
            [-1, 1, 2]

        Check that the partial quotients of an integer ``n`` is simply ``[n]``::

            sage: QQ(1).continued_fraction_list()
            [1]
            sage: QQ(0).continued_fraction_list()
            [0]
            sage: QQ(-1).continued_fraction_list()
            [-1]

        Hirzebruch-Jung continued fractions::

            sage: (11/19).continued_fraction_list("hj")
            [1, 3, 2, 3, 2]
            sage: 1 - 1/(3 - 1/(2 - 1/(3 - 1/2)))
            11/19

            sage: (225/137).continued_fraction_list("hj")
            [2, 3, 5, 10]
            sage: 2 - 1/(3 - 1/(5 - 1/10))
            225/137

            sage: (-23/19).continued_fraction_list("hj")
            [-1, 5, 4]
            sage: -1 - 1/(5 - 1/4)
            -23/19
        """
        cdef Integer z
        cdef mpz_t p,q,tmp
        cdef list res = []

        mpz_init(tmp)
        mpz_init(p)
        mpz_init(q)
        mpz_set(p, mpq_numref(self.value))
        mpz_set(q, mpq_denref(self.value))

        if type == "std":
            while mpz_sgn(q) != 0:
                z = Integer.__new__(Integer)
                mpz_fdiv_qr(z.value,tmp,p,q)
                mpz_set(p,q)
                mpz_set(q,tmp)
                res.append(z)
        elif type == "hj":
            while mpz_sgn(q) != 0:
                z = Integer.__new__(Integer)
                mpz_cdiv_qr(z.value,tmp,p,q)
                mpz_set(p,q)
                mpz_set(q,tmp)
                res.append(z)
                if mpz_sgn(q) == 0:
                    break
                z = Integer.__new__(Integer)
                mpz_fdiv_qr(z.value,tmp,p,q)
                mpz_set(p,q)
                mpz_set(q,tmp)
                mpz_neg(z.value,z.value)
                res.append(z)
        else:
            mpz_clear(p)
            mpz_clear(q)
            mpz_clear(tmp)
            raise ValueError("the type must be one of 'floor', 'hj'")

        mpz_clear(p)
        mpz_clear(q)
        mpz_clear(tmp)

        return res

    def continued_fraction(self):
        r"""
        Return the continued fraction of that rational.

        EXAMPLES::

            sage: (641/472).continued_fraction()
            [1; 2, 1, 3, 1, 4, 1, 5]

            sage: a = (355/113).continued_fraction(); a
            [3; 7, 16]
            sage: a.n(digits=10)                                                        # needs sage.rings.real_mpfr
            3.141592920
            sage: pi.n(digits=10)                                                       # needs sage.rings.real_mpfr sage.symbolic
            3.141592654

        It's almost pi!
        """
        #TODO: do better
        from sage.rings.continued_fraction import ContinuedFraction_periodic
        l = self.continued_fraction_list()
        return ContinuedFraction_periodic(l)

    cpdef _richcmp_(left, right, int op):
        """
        Compare two rational numbers.

        INPUT:

        - ``left``, ``right`` -- objects

        - ``op`` -- integer

        EXAMPLES::

            sage: 1/3 < 2/3
            True
            sage: 2/3 < 1/3
            False
            sage: 4/5 < 2.0
            True
            sage: 4/5 < 0.8
            False

            sage: ones = [1, 1r, 1l, 1/1, 1.0r, 1.0]
            sage: twos = [2, 2r, 2l, 2/1, 2.0r, 2.0]
            sage: threes = [3, 3r, 3l, 3/1, 3.0r, 3.0]
            sage: from itertools import product
            sage: for one,two,three in product(ones,twos,threes):
            ....:     assert one < two < three
            ....:     assert one <= two <= three
            ....:     assert three > two > one
            ....:     assert three >= two >= one
            ....:     assert one != two and one != three and two != three
            sage: for one1, one2 in product(ones,repeat=2):
            ....:     assert (one1 == one2) is True
            ....:     assert (one1 <= one2) is True
            ....:     assert (one1 >= one2) is True

        Comparisons with gmpy2 values (:issue:`28394`)::

            sage: import gmpy2
            sage: values = [(-2,5),(-1,3),(0,1),(2,9),(1,1),(73,2)]
            sage: for num1, den1 in values:
            ....:     for num2, den2 in values:
            ....:         a1 = QQ((num1, den1))
            ....:         a2 = QQ((num2, den2))
            ....:         b1 = gmpy2.mpq(num1, den1)
            ....:         b2 = gmpy2.mpq(num2, den2)
            ....:         assert a1 == b1 and b1 == a1 and a2 == b2 and b2 == a2
            ....:         assert (a1 == a2) == (b1 == b2) == (a1 == b2) == (b1 == a2)
            ....:         assert (a1 != a2) == (b1 != b2) == (a1 != b2) == (b1 != a2)
            ....:         assert (a1 <  a2) == (b1 <  b2) == (a1 <  b2) == (b1 <  a2)
            ....:         assert (a1 <= a2) == (b1 <= b2) == (a1 <= b2) == (b1 <= a2)
            ....:         assert (a1 >  a2) == (b1 >  b2) == (a1 >  b2) == (b1 >  a2)
            ....:         assert (a1 >= a2) == (b1 >= b2) == (a1 >= b2) == (b1 >= a2)
        """
        cdef int c
        if op == Py_EQ:
            return <bint> mpq_equal((<Rational>left).value,
                                    (<Rational>right).value)
        elif op == Py_NE:
            return not mpq_equal((<Rational>left).value,
                                 (<Rational>right).value)
        c = mpq_cmp((<Rational>left).value, (<Rational>right).value)
        return rich_to_bool_sgn(op, c)

    def __copy__(self):
        """
        EXAMPLES::

            sage: a = -17/37
            sage: copy(a) is a
            True

        Coercion does not make a new copy::

            sage: QQ(a) is a
            True

        Calling the constructor directly makes a new copy::

            sage: Rational(a) is a
            False
        """
        # immutable
        return self

    def __deepcopy__(self, memo):
        """
        EXAMPLES::

            sage: a = -17/37
            sage: deepcopy(a) is a
            True
        """
        # immutable
        return self

    def __dealloc__(self):
        """
        Free memory occupied by this rational number.

        EXAMPLES::

            sage: a = -17/37
            sage: del a          # indirect test
        """
        mpq_clear(self.value)

    def __repr__(self):
        """
        Return string representation of this rational number.

        EXAMPLES::

            sage: a = -17/37; a.__repr__()
            '-17/37'
        """
        return self.str()

    def _latex_(self):
        """
        Return Latex representation of this rational number.

        EXAMPLES::

            sage: a = -17/37
            sage: a._latex_()
            '-\\frac{17}{37}'
        """
        if self.denom() == 1:
            return str(self.numer())
        if self < 0:
            return "-\\frac{%s}{%s}" % (-self.numer(), self.denom())
        else:
            return "\\frac{%s}{%s}" % (self.numer(), self.denom())

    def _symbolic_(self, sring):
        """
        Return this rational as symbolic expression.

        EXAMPLES::

            sage: ex = SR(QQ(7)/3); ex                                                  # needs sage.symbolic
            7/3
            sage: parent(ex)                                                            # needs sage.symbolic
            Symbolic Ring
        """
        return sring._force_pyobject(self, force=True)

    def _sympy_(self):
        """
        Convert Sage ``Rational`` to SymPy ``Rational``.

        EXAMPLES::

            sage: # needs sympy
            sage: n = 1/2; n._sympy_()
            1/2
            sage: n = -1/5; n._sympy_()
            -1/5
            sage: from sympy import Symbol
            sage: QQ(1) + Symbol('x')*QQ(2)
            2*x + 1
        """
        import sympy
        return sympy.Rational(int(self.numerator()), int(self.denominator()))

    def __mpz__(self):
        """
        Return a gmpy2 ``mpz`` if this Rational is an integer.

        EXAMPLES::

            sage: q = 6/2
            sage: q.__mpz__()
            mpz(3)
            sage: q = 1/4
            sage: q.__mpz__()
            Traceback (most recent call last):
            ...
            TypeError: unable to convert rational 1/4 to an integer

        TESTS::

            sage: QQ().__mpz__(); raise NotImplementedError("gmpy2 is not installed")
            Traceback (most recent call last):
            ...
            NotImplementedError: gmpy2 is not installed
        """
        if self.denominator() != 1:
            raise TypeError(f"unable to convert rational {self} to an integer")
        return self.numerator().__mpz__()

    def __mpq__(self):
        """
        Convert Sage ``Rational`` to gmpy2 ``Rational``.

        EXAMPLES::

            sage: r = 5/3
            sage: r.__mpq__()
            mpq(5,3)
            sage: from gmpy2 import mpq
            sage: mpq(r)
            mpq(5,3)
        """
        return gmpy2.GMPy_MPQ_From_mpq(self.value)

    def _magma_init_(self, magma):
        """
        Return the magma representation of ``self``.

        EXAMPLES::

            sage: n = -485/82847
            sage: n._magma_init_(magma)                         # optional - magma
            '-485/82847'
        """
        return self.numerator()._magma_init_(magma) + '/' + self.denominator()._magma_init_(magma)

    @property
    def __array_interface__(self):
        """
        Used for NumPy conversion. If ``self`` is integral, it converts to
        an ``Integer``. Otherwise it converts to a double floating point
        value.

        EXAMPLES::

            sage: # needs numpy
            sage: import numpy
            sage: numpy.array([1, 2, 3/1])
            array([1, 2, 3])
            sage: numpy.array(QQ(2**40)).dtype
            dtype('int64')
            sage: numpy.array(QQ(2**400)).dtype
            dtype('O')
            sage: numpy.array([1, 1/2, 3/4])
            array([1.  , 0.5 , 0.75])
        """
        if mpz_cmp_ui(mpq_denref(self.value), 1) == 0:
            if mpz_fits_slong_p(mpq_numref(self.value)):
                return numpy_long_interface
            elif sizeof(long) == 4 and mpz_sizeinbase(mpq_numref(self.value), 2) <= 63:
                return numpy_int64_interface
            else:
                return numpy_object_interface
        else:
            return numpy_double_interface

    def _mathml_(self):
        """
        Return mathml representation of this rational number.

        EXAMPLES::

            sage: a = -17/37; a._mathml_()
            '<mo>-</mo><mfrac><mrow><mn>17</mn></mrow><mrow><mn>37</mn></mrow></mfrac>'
        """
        if self.denom() == 1:
            return '<mn>%s</mn>'%(self.numer())
        else:
            from sage.misc.mathml import mathml
            t = ''
            if self < 0:
                t = t + '<mo>-</mo>'
            t = t + '<mfrac><mrow>%s</mrow><mrow>%s</mrow></mfrac>'%(
                mathml(abs(self.numer())), mathml(self.denom()))
            return t

    def _im_gens_(self, codomain, im_gens, base_map=None):
        """
        Return the image of ``self`` under the homomorphism from the rational
        field to ``codomain``.

        This always just returns ``self`` coerced into the ``codomain``.

        INPUT:

        - ``codomain`` -- object (usually a ring)

        - ``im_gens`` -- list of elements of ``codomain``

        EXAMPLES::

            sage: a = -17/37
            sage: a._im_gens_(QQ, [1/1])
            -17/37
        """
        return codomain.coerce(self)

    def content(self, other):
        """
        Return the content of ``self`` and ``other``, i.e., the unique positive
        rational number `c` such that ``self/c`` and ``other/c`` are coprime
        integers.

        ``other`` can be a rational number or a list of rational numbers.

        EXAMPLES::

            sage: a = 2/3
            sage: a.content(2/3)
            2/3
            sage: a.content(1/5)
            1/15
            sage: a.content([2/5, 4/9])
            2/45
        """
        from sage.structure.sequence import Sequence
        seq = Sequence(other)
        seq.append(self)
        nums = [x.numerator() for x in seq]
        denoms = [x.denominator() for x in seq]
        from sage.arith.misc import GCD as gcd
        from sage.arith.functions import lcm
        return gcd(nums) / lcm(denoms)

    def valuation(self, p):
        r"""
        Return the power of ``p`` in the factorization of ``self``.

        INPUT:

        - ``p`` -- a prime number

        OUTPUT:

        (integer or infinity) ``Infinity`` if ``self`` is zero, otherwise the
        (positive or negative) integer `e` such that ``self`` = `m*p^e`
        with `m` coprime to `p`.

        .. NOTE::

           See also :meth:`val_unit()` which returns the pair `(e,m)`. The
           function :meth:`ord()` is an alias for :meth:`valuation()`.

        EXAMPLES::

            sage: x = -5/9
            sage: x.valuation(5)
            1
            sage: x.ord(5)
            1
            sage: x.valuation(3)
            -2
            sage: x.valuation(2)
            0

        Some edge cases::

            sage: (0/1).valuation(4)
            +Infinity
            sage: (7/16).valuation(4)
            -2
        """
        return self.numerator().valuation(p) - self.denominator().valuation(p)

    ord = valuation

    def local_height(self, p, prec=None):
        r"""
        Return the local height of this rational number at the prime `p`.

        INPUT:

        - ``p`` -- a prime number

        - ``prec`` -- integer (default: default :class:`RealField` precision);
          desired floating point precision

        OUTPUT:

        (real) The local height of this rational number at the
        prime `p`.

        EXAMPLES::

            sage: a = QQ(25/6)
            sage: a.local_height(2)                                                     # needs sage.rings.real_mpfr
            0.693147180559945
            sage: a.local_height(3)                                                     # needs sage.rings.real_mpfr
            1.09861228866811
            sage: a.local_height(5)                                                     # needs sage.rings.real_mpfr
            0.000000000000000
        """
        from sage.rings.real_mpfr import RealField
        if prec is None:
            R = RealField()
        else:
            R = RealField(prec)
        if self.is_zero():
            return R.zero()
        val = self.valuation(p)
        if val >= 0:
            return R.zero()
        return -val * R(p).log()

    def local_height_arch(self, prec=None):
        r"""
        Return the Archimedean local height of this rational number at the
        infinite place.

        INPUT:

        - ``prec`` -- integer (default: default :class:`RealField` precision);
          desired floating point precision

        OUTPUT:

        (real) The local height of this rational number `x` at the
        unique infinite place of `\QQ`, which is
        `\max(\log(|x|),0)`.

        EXAMPLES::

            sage: a = QQ(6/25)
            sage: a.local_height_arch()                                                 # needs sage.rings.real_mpfr
            0.000000000000000
            sage: (1/a).local_height_arch()                                             # needs sage.rings.real_mpfr
            1.42711635564015
            sage: (1/a).local_height_arch(100)                                          # needs sage.rings.real_mpfr
            1.4271163556401457483890413081
        """
        from sage.rings.real_mpfr import RealField
        if prec is None:
            R = RealField()
        else:
            R = RealField(prec)
        a = self.abs()
        if a <= 1:
            return R.zero()
        return R(a).log()

    def global_height_non_arch(self, prec=None):
        r"""
        Return the total non-archimedean component of the height of this
        rational number.

        INPUT:

        - ``prec`` -- integer (default: default :class:`RealField` precision);
          desired floating point precision

        OUTPUT:

        (real) The total non-archimedean component of the height of
        this rational number.

        ALGORITHM:

        This is the sum of the local heights at all primes `p`, which
        may be computed without factorization as the log of the
        denominator.

        EXAMPLES::

            sage: a = QQ(5/6)
            sage: a.support()
            [2, 3, 5]
            sage: a.global_height_non_arch()                                            # needs sage.rings.real_mpfr
            1.79175946922805
            sage: [a.local_height(p) for p in a.support()]                              # needs sage.rings.real_mpfr
            [0.693147180559945, 1.09861228866811, 0.000000000000000]
            sage: sum([a.local_height(p) for p in a.support()])                         # needs sage.rings.real_mpfr
            1.79175946922805
        """
        from sage.rings.real_mpfr import RealField
        if prec is None:
            R = RealField()
        else:
            R = RealField(prec)
        d = self.denominator()
        if d.is_one():
            return R.zero()
        return R(d).log()

    def global_height_arch(self, prec=None):
        r"""
        Return the total archimedean component of the height of this rational
        number.

        INPUT:

        - ``prec`` -- integer (default: default :class:`RealField` precision);
          desired floating point precision

        OUTPUT:

        (real) The total archimedean component of the height of
        this rational number.

        ALGORITHM:

        Since `\QQ` has only one infinite place this is just the value
        of the local height at that place.  This separate function is
        included for compatibility with number fields.

        EXAMPLES::

            sage: a = QQ(6/25)
            sage: a.global_height_arch()                                                # needs sage.rings.real_mpfr
            0.000000000000000
            sage: (1/a).global_height_arch()                                            # needs sage.rings.real_mpfr
            1.42711635564015
            sage: (1/a).global_height_arch(100)                                         # needs sage.rings.real_mpfr
            1.4271163556401457483890413081
        """
        return self.local_height_arch(prec)

    def global_height(self, prec=None):
        r"""
        Return the absolute logarithmic height of this rational number.

        INPUT:

        - ``prec`` -- integer (default: default :class:`RealField` precision);
          desired floating point precision

        OUTPUT:

        (real) The absolute logarithmic height of this rational number.

        ALGORITHM:

        The height is the sum of the total archimedean and
        non-archimedean components, which is equal to
        `\max(\log(n),\log(d))` where `n,d` are the numerator and
        denominator of the rational number.

        EXAMPLES::

            sage: # needs sage.rings.real_mpfr
            sage: a = QQ(6/25)
            sage: a.global_height_arch() + a.global_height_non_arch()
            3.21887582486820
            sage: a.global_height()
            3.21887582486820
            sage: (1/a).global_height()
            3.21887582486820
            sage: QQ(0).global_height()
            0.000000000000000
            sage: QQ(1).global_height()
            0.000000000000000
        """
        from sage.rings.real_mpfr import RealField
        if prec is None:
            R = RealField()
        else:
            R = RealField(prec)
        return R(max(self.numerator().abs(),self.denominator())).log()

    def is_square(self):
        """
        Return whether or not this rational number is a square.

        OUTPUT: boolean

        EXAMPLES::

            sage: x = 9/4
            sage: x.is_square()
            True
            sage: x = (7/53)^100
            sage: x.is_square()
            True
            sage: x = 4/3
            sage: x.is_square()
            False
            sage: x = -1/4
            sage: x.is_square()
            False
        """
        return mpq_sgn(self.value) >= 0 and mpz_perfect_square_p(mpq_numref(self.value)) and mpz_perfect_square_p(mpq_denref(self.value))

    def is_norm(self, L, element=False, proof=True):
        r"""
        Determine whether ``self`` is the norm of an element of ``L``.

        INPUT:

        - ``L`` -- a number field
        - ``element`` -- boolean (default: ``False``); whether to also output
          an element of which ``self`` is a norm
        - ``proof`` -- if ``True``, then the output is correct unconditionally;
          if ``False``, then the output assumes GRH

        OUTPUT:

        If element is ``False``, then the output is a boolean ``B``, which is
        ``True`` if and only if ``self`` is the norm of an element of ``L``.
        If ``element`` is ``False``, then the output is a pair ``(B, x)``,
        where ``B`` is as above. If ``B`` is ``True``, then ``x`` an element of
        ``L`` such that ``self == x.norm()``. Otherwise, ``x is None``.

        ALGORITHM:

        Uses the PARI function :pari:`bnfisnorm`. See :meth:`_bnfisnorm()`.

        EXAMPLES::

            sage: # needs sage.rings.number_field
            sage: x = polygen(QQ, 'x')
            sage: K = NumberField(x^2 - 2, 'beta')
            sage: (1/7).is_norm(K)
            True
            sage: (1/10).is_norm(K)
            False
            sage: 0.is_norm(K)
            True
            sage: (1/7).is_norm(K, element=True)
            (True, 1/7*beta + 3/7)
            sage: (1/10).is_norm(K, element=True)
            (False, None)
            sage: (1/691).is_norm(QQ, element=True)
            (True, 1/691)

        The number field doesn't have to be defined by an
        integral polynomial::

            sage: B, e = (1/5).is_norm(QuadraticField(5/4, 'a'), element=True)          # needs sage.rings.number_field
            sage: B                                                                     # needs sage.rings.number_field
            True
            sage: e.norm()                                                              # needs sage.rings.number_field
            1/5

        A non-Galois number field::

            sage: # needs sage.rings.number_field
            sage: K.<a> = NumberField(x^3 - 2)
            sage: B, e = (3/5).is_norm(K, element=True); B
            True
            sage: e.norm()
            3/5
            sage: 7.is_norm(K)                                                          # needs sage.groups
            Traceback (most recent call last):
            ...
            NotImplementedError: is_norm is not implemented unconditionally
             for norms from non-Galois number fields
            sage: 7.is_norm(K, proof=False)
            False

        AUTHORS:

        - Craig Citro (2008-04-05)

        - Marco Streng (2010-12-03)
        """
        if not element:
            return self.is_norm(L, element=True, proof=proof)[0]

        from sage.rings.number_field.number_field_base import NumberField
        if not isinstance(L, NumberField):
            raise ValueError("L (=%s) must be a NumberField in is_norm" % L)
        if L.degree() == 1 or self.is_zero():
            return True, L(self)
        d = L.polynomial().denominator()
        if not d == 1:
            M, M_to_L = L.subfield(L.gen()*d)
            b, x = self.is_norm(M, element=True, proof=proof)
            if b:
                x = M_to_L(x)
            return b, x
        a, b = self._bnfisnorm(L, proof=proof)
        if b == 1:
            assert a.norm() == self
            return True, a
        if L.is_galois():
            return False, None
        M = L.galois_closure('a')
        from sage.functions.log import log
        from sage.functions.other import floor
        extra_primes = floor(12*log(abs(M.discriminant()))**2)
        a, b = self._bnfisnorm(L, proof=proof, extra_primes=extra_primes)
        if b == 1:
            assert a.norm() == self
            return True, a
        if proof:
            raise NotImplementedError("is_norm is not implemented unconditionally for norms from non-Galois number fields")
        return False, None

    def _bnfisnorm(self, K, proof=True, extra_primes=0):
        r"""
        Return the output of the PARI function :pari:`bnfisnorm`.

        Tries to tell whether the rational number ``self`` is the norm of some
        element `y` in ``K``. Returns a pair `(a, b)` where
        ``self = Norm(a)*b``. Looks for a solution that is an `S`-unit, with
        `S` a certain set of prime ideals containing (among others) all primes
        dividing ``self``.

        If `K` is known to be Galois, set ``extra_primes = 0`` (in this case,
        ``self`` is a norm iff `b = 1`).

        If ``extra_primes`` is nonzero, the program adds to `S` the following
        prime ideals, depending on the sign of extra_primes.
        If ``extra_primes > 0``, the ideals of norm less than ``extra_primes``.
        And if ``extra_primes < 0``, the ideals dividing ``extra_primes``.

        Assuming GRH, the answer is guaranteed (i.e., ``self`` is a norm
        iff `b = 1`), if `S` contains all primes less than
        `12\log(\disc(L))^2`,
        where `L` is the Galois closure of `K`.

        INPUT:

        - ``K`` -- a number field
        - ``proof`` -- whether to certify the output of bnfinit;
          if ``False``, then correctness of the output depends on GRH
        - ``extra_primes`` -- integer as explained above

        OUTPUT:

        A pair `(a, b)` with `a` in `K` and `b` in `\QQ` such that
        ``self == Norm(a)*b`` as explained above.

        ALGORITHM:

        Uses PARI's bnfisnorm.

        EXAMPLES::

            sage: QQ(2)._bnfisnorm(QuadraticField(-1, 'i'))                             # needs sage.rings.number_field
            (i + 1, 1)
            sage: x = polygen(QQ, 'x')
            sage: 7._bnfisnorm(NumberField(x^3 - 2, 'b'))                               # needs sage.rings.number_field
            (1, 7)

        AUTHORS:

        - Craig Citro (2008-04-05)

        - Marco Streng (2010-12-03)
        """
        from sage.rings.number_field.number_field_base import NumberField
        if not isinstance(K, NumberField):
            raise ValueError("K must be a NumberField in bnfisnorm")

        a, b = K.pari_bnf(proof=proof).bnfisnorm(self, flag=extra_primes)
        return K(a), Rational(b)

    def is_perfect_power(self, expected_value=False):
        r"""
        Return ``True`` if ``self`` is a perfect power.

        INPUT:

        - ``expected_value`` -- boolean; whether or not this rational is
          expected to be a perfect power. This does not affect the correctness
          of the output, only the runtime.

        If ``expected_value`` is ``False`` (default) it will check the
        smallest of the numerator and denominator is a perfect power
        as a first step, which is often faster than checking if the
        quotient is a perfect power.

        EXAMPLES::

            sage: (4/9).is_perfect_power()
            True
            sage: (144/1).is_perfect_power()
            True
            sage: (4/3).is_perfect_power()
            False
            sage: (2/27).is_perfect_power()
            False
            sage: (4/27).is_perfect_power()
            False
            sage: (-1/25).is_perfect_power()
            False
            sage: (-1/27).is_perfect_power()
            True
            sage: (0/1).is_perfect_power()
            True

        The second parameter does not change the result, but may
        change the runtime.

        ::

            sage: (-1/27).is_perfect_power(True)
            True
            sage: (-1/25).is_perfect_power(True)
            False
            sage: (2/27).is_perfect_power(True)
            False
            sage: (144/1).is_perfect_power(True)
            True

        This test makes sure we workaround a bug in GMP (see :issue:`4612`)::

            sage: [-a for a in srange(100) if not QQ(-a^3).is_perfect_power()]
            []
            sage: [-a for a in srange(100) if not QQ(-a^3).is_perfect_power(True)]
            []
        """
        cdef int s

        if (mpz_cmp_ui(mpq_numref(self.value), 0) == 0):
            return True
        elif (mpz_cmp_ui(mpq_numref(self.value), 1) == 0):
            return mpz_perfect_power_p(mpq_denref(self.value))

        cdef mpz_t prod
        cdef bint res

        # We should be able to run the code in the sign == 1 case
        # below for both cases. However, we need to do extra work to
        # avoid a bug in GMP's mpz_perfect_power_p; see Issue #4612 for
        # more details.
        #
        # The code in the case of sign == -1 could definitely be
        # cleaned up, but it will be removed shortly, since both GMP
        # and eMPIRe have fixes for the mpz_perfect_power_p bug.

        s = mpz_sgn(mpq_numref(self.value))
        if s == 1: # self is positive

            if (mpz_cmp_ui(mpq_denref(self.value), 1) == 0):
                return mpz_perfect_power_p(mpq_numref(self.value))
            if not expected_value:
                # A necessary condition is that both the numerator and denominator
                # be perfect powers, which can be faster to disprove than the full
                # product (especially if both have a large prime factor).
                if mpz_cmpabs(mpq_numref(self.value), mpq_denref(self.value)) < 0:
                    if not mpz_perfect_power_p(mpq_numref(self.value)):
                        return False
                else:
                    if not mpz_perfect_power_p(mpq_denref(self.value)):
                        return False
            mpz_init(prod)
            mpz_mul(prod, mpq_numref(self.value), mpq_denref(self.value))
            res = mpz_perfect_power_p(prod)
            mpz_clear(prod)
            return res == 1

        else: # self is negative

            if (mpz_cmp_ui(mpq_denref(self.value), 1) == 0):
                if (mpz_cmp_si(mpq_numref(self.value), -1) == 0):
                    return True
                mpz_init(prod)
                mpz_mul_si(prod, mpq_numref(self.value), -1)
                while mpz_perfect_square_p(prod):
                    mpz_sqrt(prod, prod)
                s = mpz_perfect_power_p(prod)
                mpz_clear(prod)
                return s == 1

            if not expected_value:
                if mpz_cmpabs(mpq_numref(self.value), mpq_denref(self.value)) < 0:
                    mpz_init(prod)
                    mpz_mul_si(prod, mpq_numref(self.value), -1)
                    if mpz_cmp_ui(prod, 1) != 0:
                        while mpz_perfect_square_p(prod):
                            mpz_sqrt(prod, prod)
                        if not mpz_perfect_power_p(prod):
                            mpz_clear(prod)
                            return False
                else:
                    if not mpz_perfect_power_p(mpq_denref(self.value)):
                        return False
                    mpz_init(prod)
            else:
                mpz_init(prod)

            mpz_mul(prod, mpq_numref(self.value), mpq_denref(self.value))
            mpz_mul_si(prod, prod, -1)
            while mpz_perfect_square_p(prod):
                mpz_sqrt(prod, prod)
            res = mpz_perfect_power_p(prod)
            mpz_clear(prod)
            return res == 1

    def squarefree_part(self):
        """
        Return the square free part of `x`, i.e., an integer `z` such
        that `x = z y^2`, for a perfect square `y^2`.

        EXAMPLES::

            sage: a = 1/2
            sage: a.squarefree_part()
            2
            sage: b = a/a.squarefree_part()
            sage: b, b.is_square()
            (1/4, True)
            sage: a = 24/5
            sage: a.squarefree_part()
            30
        """
        return self.numer().squarefree_part() * self.denom().squarefree_part()

    def is_padic_square(self, p, check=True):
        r"""
        Determines whether this rational number is a square in `\QQ_p` (or in
        `R` when ``p = infinity``).

        INPUT:

        - ``p`` -- a prime number, or ``infinity``

        - ``check`` -- boolean (default: ``True``); check if `p` is prime

        EXAMPLES::

            sage: QQ(2).is_padic_square(7)
            True
            sage: QQ(98).is_padic_square(7)
            True
            sage: QQ(2).is_padic_square(5)
            False

        TESTS::

            sage: QQ(5/7).is_padic_square(int(2))
            False
        """
        ## Special case when self is zero
        if self.is_zero():
            return True

        ## Deal with p = infinity (i.e. the real numbers)
        import sage.rings.infinity
        if p == sage.rings.infinity.infinity:
            return (self > 0)

        ## Check that p is prime
        from sage.rings.integer_ring import ZZ
        p = ZZ(p)
        if check and not p.is_prime():
            raise ValueError('p must be "infinity" or a positive prime number.')

        ## Deal with finite primes
        e, m = self.val_unit(p)

        if e % 2:
            return False

        if p == 2:
            return ((m % 8) == 1)

        from sage.arith.misc import kronecker as kronecker_symbol
        return (kronecker_symbol(m, p) == 1)

    def val_unit(self, p):
        r"""
        Return a pair: the `p`-adic valuation of ``self``, and the `p`-adic
        unit of ``self``, as a :class:`Rational`.

        We do not require the `p` be prime, but it must be at least 2. For
        more documentation see :meth:`Integer.val_unit()`.

        INPUT:

        - ``p`` -- a prime

        OUTPUT:

        - integer; the `p`-adic valuation of this rational

        - ``Rational``; `p`-adic unit part of ``self``

        EXAMPLES::

            sage: (-4/17).val_unit(2)
            (2, -1/17)
            sage: (-4/17).val_unit(17)
            (-1, -4)
            sage: (0/1).val_unit(17)
            (+Infinity, 1)

        AUTHORS:

        - David Roe (2007-04-12)
        """
        return self._val_unit(p)

    # TODO -- change to use cpdef?  If so, must fix
    # code in padics, etc.  Do search_src('_val_unit').
    cdef _val_unit(Rational self, integer.Integer p):
        """
        This is called by :meth:`val_unit()`.

        EXAMPLES::

            sage: (-4/17).val_unit(2) # indirect doctest
            (2, -1/17)
        """
        cdef Rational u
        if mpz_cmp_ui(p.value, 2) < 0:
            raise ValueError("p must be at least 2.")
        if mpq_sgn(self.value) == 0:
            import sage.rings.infinity
            u = Rational.__new__(Rational)
            mpq_set_ui(u.value, 1, 1)
            return (sage.rings.infinity.infinity, u)
        cdef Integer v = Integer.__new__(Integer)
        u = Rational.__new__(Rational)
        sig_on()
        mpz_set_ui(v.value, mpz_remove(mpq_numref(u.value), mpq_numref(self.value), p.value))
        sig_off()
        if mpz_sgn(v.value) != 0:
            mpz_set(mpq_denref(u.value), mpq_denref(self.value))
        else:
            sig_on()
            mpz_set_ui(v.value, mpz_remove(mpq_denref(u.value), mpq_denref(self.value), p.value))
            sig_off()
            mpz_neg(v.value, v.value)
        return (v, u)

    def prime_to_S_part(self, S=[]):
        r"""
        Return ``self`` with all powers of all primes in ``S`` removed.

        INPUT:

        - ``S`` -- list or tuple of primes

        OUTPUT: rational

        .. NOTE::

           Primality of the entries in `S` is not checked.

        EXAMPLES::

            sage: QQ(3/4).prime_to_S_part()
            3/4
            sage: QQ(3/4).prime_to_S_part([2])
            3
            sage: QQ(-3/4).prime_to_S_part([3])
            -1/4
            sage: QQ(700/99).prime_to_S_part([2,3,5])
            7/11
            sage: QQ(-700/99).prime_to_S_part([2,3,5])
            -7/11
            sage: QQ(0).prime_to_S_part([2,3,5])
            0
            sage: QQ(-700/99).prime_to_S_part([])
            -700/99
        """
        if self.is_zero():
            return self
        a = self
        for p in S:
            e, a = a.val_unit(p)
        return a

    def sqrt(self, prec=None, extend=True, all=False):
        r"""
        The square root function.

        INPUT:

        - ``prec`` -- integer (default: ``None``); if ``None``, returns
          an exact square root; otherwise returns a numerical square root if
          necessary, to the given bits of precision.

        - ``extend`` -- boolean (default: ``True``); if ``True``, return a
          square root in an extension ring, if necessary. Otherwise, raise a
          :exc:`ValueError` if the square is not in the base ring. Ignored if
          ``prec`` is not ``None``.

        - ``all`` -- boolean (default: ``False``); if ``True``, return all
          square roots of ``self`` (a list of length 0, 1, or 2)

        EXAMPLES::

            sage: x = 25/9
            sage: x.sqrt()
            5/3
            sage: sqrt(x)
            5/3
            sage: x = 64/4
            sage: x.sqrt()
            4
            sage: x = 100/1
            sage: x.sqrt()
            10
            sage: x.sqrt(all=True)
            [10, -10]
            sage: x = 81/5
            sage: x.sqrt()                                                              # needs sage.symbolic
            9*sqrt(1/5)
            sage: x = -81/3
            sage: x.sqrt()                                                              # needs sage.symbolic
            3*sqrt(-3)

        ::

            sage: n = 2/3
            sage: n.sqrt()                                                              # needs sage.symbolic
            sqrt(2/3)

            sage: # needs sage.rings.real_mpfr
            sage: n.sqrt(prec=10)
            0.82
            sage: n.sqrt(prec=100)
            0.81649658092772603273242802490
            sage: n.sqrt(prec=100)^2
            0.66666666666666666666666666667
            sage: n.sqrt(prec=53, all=True)
            [0.816496580927726, -0.816496580927726]
            sage: sqrt(-2/3, prec=53)
            0.816496580927726*I
            sage: sqrt(-2/3, prec=53, all=True)
            [0.816496580927726*I, -0.816496580927726*I]

            sage: n.sqrt(extend=False)
            Traceback (most recent call last):
            ...
            ValueError: square root of 2/3 not a rational number
            sage: n.sqrt(extend=False, all=True)
            []
            sage: sqrt(-2/3, all=True)                                                  # needs sage.symbolic
            [sqrt(-2/3), -sqrt(-2/3)]

        TESTS:

        Ensure that :issue:`37153` is fixed, so that behaviour aligns
        with other rings and fields.
        See :issue:`9466` and :issue:`26509` for context::

            sage: QQ(3).sqrt(extend=False, all=True)
            []
            sage: QQ(-1).sqrt(extend=False, all=True)
            []

        AUTHORS:

        - Naqi Jaffery (2006-03-05): some examples
        """
        if prec is not None:
            from sage.misc.functional import _do_sqrt
            return _do_sqrt(self, prec=prec, all=all)

        if mpq_sgn(self.value) == 0:
            return [self] if all else self

        if mpq_sgn(self.value) < 0:
            if extend:
                from sage.misc.functional import _do_sqrt
                return _do_sqrt(self, prec=prec, all=all)
            if all:
                return []
            raise ValueError("square root of negative number not rational")

        cdef Rational z = <Rational> Rational.__new__(Rational)
        cdef mpz_t tmp
        cdef int non_square = 0

        sig_on()
        mpz_init(tmp)
        mpz_sqrtrem(mpq_numref(z.value), tmp, mpq_numref(self.value))
        if mpz_sgn(tmp) != 0:
            non_square = 1
        else:
            mpz_sqrtrem(mpq_denref(z.value), tmp, mpq_denref(self.value))
            if mpz_sgn(tmp) != 0:
                non_square = 1
        mpz_clear(tmp)
        sig_off()

        if non_square:
            if extend:
                from sage.misc.functional import _do_sqrt
                return _do_sqrt(self, prec=prec, all=all)
            if all:
                return []
            raise ValueError("square root of %s not a rational number" % self)

        if all:
            return [z, -z]
        return z

    def period(self):
        r"""
        Return the period of the repeating part of the decimal expansion of
        this rational number.

        ALGORITHM:

        When a rational number `n/d` with `(n,d)=1` is
        expanded, the period begins after `s` terms and has length
        `t`, where `s` and `t` are the smallest numbers satisfying
        `10^s=10^{s+t} \mod d`. In general if `d=2^a 5^b m` where `m`
        is coprime to 10, then `s=\max(a,b)` and `t` is the order of
        10 modulo `m`.

        EXAMPLES::

            sage: (1/7).period()                                                        # needs sage.libs.pari
            6
            sage: RR(1/7)                                                               # needs sage.rings.real_mpfr
            0.142857142857143
            sage: (1/8).period()                                                        # needs sage.libs.pari
            1
            sage: RR(1/8)                                                               # needs sage.rings.real_mpfr
            0.125000000000000
            sage: RR(1/6)                                                               # needs sage.rings.real_mpfr
            0.166666666666667
            sage: (1/6).period()                                                        # needs sage.libs.pari
            1
            sage: x = 333/106
            sage: x.period()                                                            # needs sage.libs.pari
            13
            sage: RealField(200)(x)                                                     # needs sage.rings.real_mpfr
            3.1415094339622641509433962264150943396226415094339622641509
        """
        cdef unsigned int alpha, beta
        d = self.denominator()
        alpha, d = d.val_unit(2)
        beta, d  = d.val_unit(5)
        from sage.rings.finite_rings.integer_mod import Mod
        return Mod(10, d).multiplicative_order()

    def nth_root(self, int n):
        r"""
        Compute the `n`-th root of ``self``, or raises a
        :exc:`ValueError` if ``self`` is not a perfect `n`-th power.

        INPUT:

        - ``n`` -- integer (must fit in C ``int`` type)

        AUTHORS:

        - David Harvey (2006-09-15)

        EXAMPLES::

            sage: (25/4).nth_root(2)
            5/2
            sage: (125/8).nth_root(3)
            5/2
            sage: (-125/8).nth_root(3)
            -5/2
            sage: (25/4).nth_root(-2)
            2/5

        ::

            sage: (9/2).nth_root(2)
            Traceback (most recent call last):
            ...
            ValueError: not a perfect 2nd power

        ::

            sage: (-25/4).nth_root(2)
            Traceback (most recent call last):
            ...
            ValueError: cannot take even root of negative number
        """
        # TODO -- this could be quicker, by using GMP directly.
        cdef integer.Integer num
        cdef integer.Integer den
        cdef int negative

        if n > 0:
            negative = 0
        elif n < 0:
            n = -n
            negative = 1
        else:
            raise ValueError("n cannot be zero")

        num, exact = self.numerator().nth_root(n, 1)
        if not exact:
            from sage.rings.integer_ring import ZZ

            raise ValueError("not a perfect %s power" % ZZ(n).ordinal_str())

        den, exact = self.denominator().nth_root(n, 1)
        if not exact:
            from sage.rings.integer_ring import ZZ

            raise ValueError("not a perfect %s power" % ZZ(n).ordinal_str())

        if negative:
            return den / num
        else:
            return num / den

    def is_nth_power(self, int n):
        r"""
        Return ``True`` if ``self`` is an `n`-th power, else ``False``.

        INPUT:

        - ``n`` -- integer (must fit in C ``int`` type)

        .. NOTE::

           Use this function when you need to test if a rational
           number is an `n`-th power, but do not need to know the value
           of its `n`-th root.  If the value is needed, use :meth:`nth_root()`.

        AUTHORS:

        - John Cremona (2009-04-04)

        EXAMPLES::

            sage: QQ(25/4).is_nth_power(2)
            True
            sage: QQ(125/8).is_nth_power(3)
            True
            sage: QQ(-125/8).is_nth_power(3)
            True
            sage: QQ(25/4).is_nth_power(-2)
            True

            sage: QQ(9/2).is_nth_power(2)
            False
            sage: QQ(-25).is_nth_power(2)
            False
        """
        if n == 0:
            raise ValueError("n cannot be zero")
        if n < 0:
            n = -n
        if not n % 2 and self < 0:
            return False
        return self.numerator().nth_root(n, 1)[1]\
               and self.denominator().nth_root(n, 1)[1]

    def str(self, int base=10):
        """
        Return a string representation of ``self`` in the given ``base``.

        INPUT:

        - ``base`` -- integer (default: 10); base must be between 2 and 36

        OUTPUT: string

        EXAMPLES::

            sage: (-4/17).str()
            '-4/17'
            sage: (-4/17).str(2)
            '-100/10001'

        Note that the base must be at most 36.

        ::

            sage: (-4/17).str(40)
            Traceback (most recent call last):
            ...
            ValueError: base (=40) must be between 2 and 36
            sage: (-4/17).str(1)
            Traceback (most recent call last):
            ...
            ValueError: base (=1) must be between 2 and 36
        """
        if base < 2 or base > 36:
            raise ValueError("base (=%s) must be between 2 and 36" % base)
        cdef size_t n
        cdef char *s

        n = mpz_sizeinbase (mpq_numref(self.value), base) \
            + mpz_sizeinbase (mpq_denref(self.value), base) + 3
        s = <char *>PyMem_Malloc(n)
        if s == NULL:
            raise MemoryError("Unable to allocate enough memory for the string representation of an integer.")

        sig_on()
        mpq_get_str(s, base, self.value)
        sig_off()
        k = char_to_str(s)
        PyMem_Free(s)
        return k

    def __float__(self):
        """
        Return floating point approximation to ``self`` as a Python float.

        OUTPUT: float

        EXAMPLES::

            sage: (-4/17).__float__()
            -0.23529411764705882
            sage: float(-4/17)
            -0.23529411764705882
            sage: float(1/3)
            0.3333333333333333
            sage: float(1/10)
            0.1
            sage: n = QQ(902834098234908209348209834092834098); float(n)
            9.028340982349083e+35

        TESTS:

        Test that conversion agrees with `RR`::

            sage: Q = [a/b for a in [-99..99] for b in [1..99]]
            sage: all(RDF(q) == RR(q) for q in Q)
            True

        Test that the conversion has correct rounding on simple rationals::

            sage: for p in [-100..100]:                                                 # needs sage.rings.real_mpfr
            ....:   for q in [1..100]:
            ....:       r = RDF(p/q)
            ....:       assert (RR(r).exact_rational() - p/q) <= r.ulp()/2

        Test larger rationals::

            sage: Q = continued_fraction(pi).convergents()[:100]                        # needs sage.symbolic
            sage: all(RDF(q) == RR(q) for q in Q)
            True

        At some point, the continued fraction and direct conversion
        to ``RDF`` should agree::

            sage: RDFpi = RDF(pi)                                                       # needs sage.symbolic
            sage: all(RDF(q) == RDFpi for q in Q[20:])                                  # needs sage.symbolic
            True
        """
        return mpq_get_d_nearest(self.value)

    def __hash__(self):
        """
        Return hash of ``self``.

        OUTPUT: integer

        EXAMPLES::

            sage: QQ(42).__hash__()
            42
            sage: QQ(1/42).__hash__()
            1488680910            # 32-bit
            -7658195599476688946  # 64-bit
            sage: n = ZZ.random_element(10^100)
            sage: hash(n) == hash(QQ(n)) or n
            True
            sage: hash(-n) == hash(-QQ(n)) or n
            True
            sage: hash(-4/17)
            -47583156            # 32-bit
            8709371129873690700  # 64-bit
        """
        cdef Py_hash_t n = mpz_pythonhash(mpq_numref(self.value))
        cdef Py_hash_t d = mpz_pythonhash(mpq_denref(self.value))
        # The constant below is (1 + sqrt(5)) << 61
        return n + (d - 1) * <Py_hash_t>(7461864723258187525)

    def __getitem__(self, int n):
        """
        Return ``n``-th element of ``self``, viewed as a list. This is for
        consistency with how number field elements work.

        INPUT:

        - ``n`` -- integer (error if not 0 or -1)

        OUTPUT: rational

        EXAMPLES::

            sage: (-4/17)[0]
            -4/17
            sage: (-4/17)[1]
            Traceback (most recent call last):
            ...
            IndexError: index n (=1) out of range; it must be 0
            sage: (-4/17)[-1]   # indexing from the right
            -4/17
        """
        if n == 0 or n == -1:
            return self
        raise IndexError(f"index n (={n}) out of range; it must be 0")

    ################################################################
    # Optimized arithmetic
    ################################################################
    def __add__(left, right):
        """
        Return ``left`` plus ``right``.

        EXAMPLES::

            sage: (2/3) + (1/6)
            5/6
            sage: (1/3) + (1/2)
            5/6
            sage: (1/3) + 2
            7/3
        """
        cdef Rational x
        if type(left) is type(right):
            x = <Rational> Rational.__new__(Rational)
            mpq_add(x.value, (<Rational>left).value, (<Rational>right).value)
            return x
        elif type(right) is Integer:
            x = <Rational> Rational.__new__(Rational)
            mpq_add_z(x.value, (<Rational>left).value, (<Integer>right).value)
            return x

        return coercion_model.bin_op(left, right, operator.add)

    cpdef _add_(self, right):
        """
        Return ``right`` plus ``self``.

        EXAMPLES::

            sage: (2/3)._add_(1/6)
            5/6
            sage: (1/3)._add_(1/2)
            5/6
        """
        cdef Rational x
        x = <Rational> Rational.__new__(Rational)
        mpq_add(x.value, self.value, (<Rational>right).value)
        return x

    def __sub__(left, right):
        """
        Return ``left`` minus ``right``.

        EXAMPLES::

            sage: 11/3 - 5/4
            29/12

            sage: (2/3) - 2
            -4/3
            sage: (-2/3) - 1
            -5/3
            sage: (2/3) - (-3)
            11/3
            sage: (-2/3) - (-3)
            7/3
            sage: 2/3 - polygen(QQ)
            -x + 2/3
        """
        cdef Rational x
        if type(left) is type(right):
            x = <Rational> Rational.__new__(Rational)
            mpq_sub(x.value, (<Rational>left).value, (<Rational>right).value)
            return x
        elif type(right) is Integer:
            x = <Rational> Rational.__new__(Rational)
            mpz_mul(mpq_numref(x.value), mpq_denref((<Rational>left).value),
                    (<Integer>right).value)
            mpz_sub(mpq_numref(x.value), mpq_numref((<Rational>left).value),
                    mpq_numref(x.value))
            mpz_set(mpq_denref(x.value), mpq_denref((<Rational>left).value))
            return x

        return coercion_model.bin_op(left, right, operator.sub)

    cpdef _sub_(self, right):
        """
        Return ``self`` minus ``right``.

        EXAMPLES::

            sage: (2/3)._sub_(1/6)
            1/2
        """
        cdef Rational x
        x = <Rational> Rational.__new__(Rational)
        mpq_sub(x.value, self.value, (<Rational>right).value)
        return x

    cpdef _neg_(self):
        """
        Negate ``self``.

        EXAMPLES::

            sage: -(2/3) # indirect doctest
            -2/3
        """
        cdef Rational x
        x = <Rational> Rational.__new__(Rational)
        mpq_neg(x.value, self.value)
        return x

    def __mul__(left, right):
        """
        Return ``left`` times ``right``.

        EXAMPLES::

            sage: (3/14) * 2/3
            1/7
            sage: (3/14) * 10
            15/7
            sage: 3/14 * polygen(QQ)
            3/14*x
        """
        cdef Rational x
        if type(left) is type(right):
            x = <Rational> Rational.__new__(Rational)
            mpq_mul(x.value, (<Rational>left).value, (<Rational>right).value)
            return x
        elif type(right) is Integer:
            x = <Rational> Rational.__new__(Rational)
            mpq_mul_z(x.value, (<Rational>left).value, (<Integer>right).value)
            return x

        return coercion_model.bin_op(left, right, operator.mul)

    cpdef _mul_(self, right):
        """
        Return ``self`` times ``right``.

        EXAMPLES::

            sage: (3/14)._mul_(2/3)
            1/7
        """
        cdef Rational x
        x = <Rational> Rational.__new__(Rational)
        if mpz_sizeinbase (mpq_numref(self.value), 2)  > 100000 or \
             mpz_sizeinbase (mpq_denref(self.value), 2) > 100000:
            # We only use the signal handler (to enable ctrl-c out) in case
            # self is huge, so the product might actually take a while to compute.
            sig_on()
            mpq_mul(x.value, self.value, (<Rational>right).value)
            sig_off()
        else:
            mpq_mul(x.value, self.value, (<Rational>right).value)
        return x

    def __truediv__(left, right):
        """
        Return ``left`` divided by ``right``.

        EXAMPLES::

            sage: QQ((2,3)) / QQ((-5,4))
            -8/15
            sage: QQ((22,3)) / 4
            11/6
            sage: QQ((-2,3)) / (-4)
            1/6
            sage: QQ((2,3)) / QQ.zero()
            Traceback (most recent call last):
            ...
            ZeroDivisionError: rational division by zero
        """
        cdef Rational x
        if type(left) is type(right):
            if mpq_cmp_si((<Rational> right).value, 0, 1) == 0:
                raise ZeroDivisionError('rational division by zero')
            x = <Rational> Rational.__new__(Rational)
            mpq_div(x.value, (<Rational>left).value, (<Rational>right).value)
            return x
        elif type(right) is Integer:
            if mpz_cmp_si((<Integer> right).value, 0) == 0:
                raise ZeroDivisionError('rational division by zero')
            x = <Rational> Rational.__new__(Rational)
            mpq_div_zz(x.value, mpq_numref((<Rational>left).value), (<Integer>right).value)
            mpz_mul(mpq_denref(x.value), mpq_denref(x.value),
                    mpq_denref((<Rational>left).value))
            return x

        return coercion_model.bin_op(left, right, operator.truediv)

    cpdef _div_(self, right):
        """
        Return ``self`` divided by ``right``.

        EXAMPLES::

            sage: 2/3 # indirect doctest
            2/3
            sage: 3/0 # indirect doctest
            Traceback (most recent call last):
            ...
            ZeroDivisionError: rational division by zero
        """
        if mpq_cmp_si((<Rational> right).value, 0, 1) == 0:
            raise ZeroDivisionError('rational division by zero')
        cdef Rational x
        x = <Rational> Rational.__new__(Rational)
        mpq_div(x.value, self.value, (<Rational>right).value)
        return x

    ################################################################
    # Other arithmetic operations.
    ################################################################

    def __invert__(self):
        """
        Return the multiplicative inverse of ``self``.

        OUTPUT: rational

        EXAMPLES::

            sage: (-4/17).__invert__()
            -17/4
            sage: ~(-4/17)
            -17/4
        """
        if self.is_zero():
            raise ZeroDivisionError('rational division by zero')
        cdef Rational x
        x = <Rational> Rational.__new__(Rational)
        mpq_inv(x.value, self.value)
        return x

    cpdef _pow_(self, other):
        """
        Raise ``self`` to the rational power ``other``.

        EXAMPLES::

            sage: (2/3)^5
            32/243
            sage: (-1/1)^(1/3)                                                          # needs sage.symbolic
            (-1)^(1/3)
            sage: (2/3)^(3/4)                                                           # needs sage.symbolic
            (2/3)^(3/4)
            sage: (-1/3)^0
            1
            sage: a = (0/1)^(0/1); a
            1
            sage: type(a)
            <class 'sage.rings.rational.Rational'>

        If the result is rational, it is returned as a rational::

            sage: a = (4/9)^(1/2); a
            2/3
            sage: parent(a)
            Rational Field
            sage: (-27/125)^(1/3)                                                       # needs sage.symbolic
            3/5*(-1)^(1/3)
            sage: (-27/125)^(1/2)                                                       # needs sage.symbolic
            3/5*sqrt(-3/5)

        The result is normalized to have the rational power in the numerator::

            sage: 2^(-1/2)                                                              # needs sage.symbolic
            1/2*sqrt(2)
            sage: 8^(-1/5)                                                              # needs sage.symbolic
            1/8*8^(4/5)
            sage: 3^(-3/2)                                                              # needs sage.symbolic
            1/9*sqrt(3)

        TESTS::

            sage: QQ(0)^(-1)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: rational division by zero

        This works even if the base is a Python integer::

            sage: int(2)^(1/2)                                                          # needs sage.symbolic
            sqrt(2)
            sage: a = int(2)^(3/1); a
            8
            sage: type(a)
            <class 'sage.rings.rational.Rational'>

        The exponent must fit in a ``long`` unless the base is -1, 0, or 1::

            sage: (1/2)^(2^100)
            Traceback (most recent call last):
            ...
            OverflowError: exponent must be at most 2147483647           # 32-bit
            OverflowError: exponent must be at most 9223372036854775807  # 64-bit
            sage: (1/2)^(-2^100)
            Traceback (most recent call last):
            ...
            OverflowError: exponent must be at most 2147483647           # 32-bit
            OverflowError: exponent must be at most 9223372036854775807  # 64-bit
            sage: QQ(-1)^(2^100)                                                        # needs sage.symbolic
            1
        """
        n = <Rational?>other

        # If the exponent happens to be an integer, we can use optimized code
        if mpz_cmp_ui(mpq_denref(n.value), 1) == 0:
            if mpz_fits_slong_p(mpq_numref(n.value)):
                return self._pow_long(mpz_get_si(mpq_numref(n.value)))

        # Perhaps it can be done exactly
        c, d = rational_power_parts(self, n)
        if d == 1:
            # It was an exact power
            return c
        elif d == -1 and n.denominator() == 2:
            # Exact rational times a power of I
            from sage.rings.imaginary_unit import I
            return c * I ** (n.numerator() % 4)

        # Result is c * d^n but we cannot simplify d^n further:
        # return a symbolic expression.
        # We use the hold=True keyword argument to prevent the
        # symbolics library from trying to simplify this expression
        # again. This would lead to infinite loops otherwise.
        from sage.symbolic.ring import SR
        return SR(c) * SR(d).power(n, hold=True)

    cpdef _pow_int(self, n):
        """
        Raise ``self`` to the integer power ``n``.

        TESTS::

            sage: a = QQ(3)._pow_int(2); a
            9
            sage: parent(a)
            Rational Field
            sage: a = (1/3)._pow_int(-2); a
            9
            sage: parent(a)
            Rational Field
        """
        # This is only called whenever n does not fit in a long
        # (otherwise _pow_long would be called). So we don't have to
        # be particularly efficient here.
        return self._pow_(Rational(n))

    cdef _pow_long(self, long n):
        """
        TESTS::

            sage: QQ(2/5) ^ 0
            1
            sage: QQ(2/5) ^ int(3)
            8/125
            sage: QQ(2/5) ^ (-3)
            125/8
            sage: QQ(2/5) ^ QQ(-3)
            125/8
            sage: QQ(-2/5) ^ int(0)
            1
            sage: QQ(-2/5) ^ 3
            -8/125
            sage: QQ(-2/5) ^ int(-3)
            -125/8
            sage: QQ(-2/5) ^ QQ(3)
            -8/125
        """
        cdef Rational x = <Rational>Rational.__new__(Rational)

        if n == 0:
            mpq_set_si(x.value, 1, 1)
        elif n > 0:
            sig_on()
            mpz_pow_ui(mpq_numref(x.value), mpq_numref(self.value), n)
            mpz_pow_ui(mpq_denref(x.value), mpq_denref(self.value), n)
            sig_off()
        else:  # n < 0:
            if mpz_sgn(mpq_numref(self.value)) == 0:
                raise ZeroDivisionError("rational division by zero")
            sig_on()
            mpz_pow_ui(mpq_denref(x.value), mpq_numref(self.value), -<unsigned long>n)
            mpz_pow_ui(mpq_numref(x.value), mpq_denref(self.value), -<unsigned long>n)
            # Fix signs
            if mpz_sgn(mpq_denref(x.value)) < 0:
                mpz_neg(mpq_numref(x.value), mpq_numref(x.value))
                mpz_neg(mpq_denref(x.value), mpq_denref(x.value))
            sig_off()
        return x

    def __pos__(self):
        """
        Return ``self``.

        OUTPUT: rational

        EXAMPLES::

            sage: (-4/17).__pos__()
            -4/17
            sage: +(-4/17)
            -4/17
        """
        return self

    def __neg__(self):
        """
        Return the negative of ``self``.

        OUTPUT: rational

        EXAMPLES::

            sage: (-4/17).__neg__()
            4/17
            sage: - (-4/17)
            4/17
        """
        cdef Rational x
        x = <Rational> Rational.__new__(Rational)
        mpq_neg(x.value, self.value)
        return x

    def __bool__(self):
        """
        Return ``True`` if this rational number is nonzero.

        OUTPUT: boolean

        EXAMPLES::

            sage: bool(0/5)
            False
            sage: bool(-4/17)
            True
        """
        # A rational number is zero iff its numerator is zero.
        return mpq_sgn(self.value) != 0

    def __abs__(self):
        """
        Return the absolute value of this rational number.

        OUTPUT: rational

        EXAMPLES::

            sage: (-4/17).__abs__()
            4/17
            sage: abs(-4/17)
            4/17
        """
        cdef Rational x
        x = <Rational> Rational.__new__(Rational)
        mpq_abs(x.value, self.value)
        return x

    def sign(self):
        """
        Return the sign of this rational number, which is -1, 0, or 1
        depending on whether this number is negative, zero, or positive
        respectively.

        OUTPUT: integer

        EXAMPLES::

            sage: (2/3).sign()
            1
            sage: (0/3).sign()
            0
            sage: (-1/6).sign()
            -1
        """
        return smallInteger(mpq_sgn(self.value))

    def mod_ui(Rational self, unsigned long int n):
        """
        Return the remainder upon division of ``self`` by the unsigned long
        integer ``n``.

        INPUT:

        - ``n`` -- an unsigned long integer

        OUTPUT: integer

        EXAMPLES::

            sage: (-4/17).mod_ui(3)
            1
            sage: (-4/17).mod_ui(17)
            Traceback (most recent call last):
            ...
            ArithmeticError: The inverse of 0 modulo 17 is not defined.
        """
        cdef unsigned int num, den, a

        # Documentation from GMP manual:
        # "For the ui variants the return value is the remainder, and
        # in fact returning the remainder is all the div_ui functions do."
        sig_on()
        num = mpz_fdiv_ui(mpq_numref(self.value), n)
        den = mpz_fdiv_ui(mpq_denref(self.value), n)
        sig_off()
        return int((num * ai.inverse_mod_int(den, n)) % n)

    def __mod__(x, y):
        """
        Return the remainder of division of ``x`` by ``y``, where ``y`` is
        something that can be coerced to an integer.

        INPUT:

        - ``other`` -- object that coerces to an integer

        OUTPUT: integer

        EXAMPLES::

            sage: (-4/17).__mod__(3/1)
            1

        TESTS:

        Check that :issue:`14870` is fixed::

            sage: int(4) % QQ(3)
            1
        """
        cdef Rational rat
        if not isinstance(x, Rational):
            rat = Rational(x)
        else:
            rat = x
        cdef other = integer.Integer(y)
        if not other:
            raise ZeroDivisionError("Rational modulo by zero")
        n = rat.numer() % other
        d = rat.denom() % other
        d = d.inverse_mod(other)
        return (n * d) % other

    def norm(self):
        r"""
        Return the norm from `\QQ` to `\QQ` of `x` (which is just `x`). This
        was added for compatibility with :class:`NumberField`.

        OUTPUT: ``Rational`` -- reference to ``self``

        EXAMPLES::

            sage: (1/3).norm()
             1/3

        AUTHORS:

        - Craig Citro
        """
        return self

    def relative_norm(self):
        """
        Return the norm from Q to Q of x (which is just x). This was added for
        compatibility with NumberFields.

        EXAMPLES::

            sage: (6/5).relative_norm()
            6/5

            sage: QQ(7/5).relative_norm()
            7/5
        """
        return self

    def absolute_norm(self):
        """
        Return the norm from Q to Q of x (which is just x). This was added for
        compatibility with NumberFields.

        EXAMPLES::

            sage: (6/5).absolute_norm()
            6/5

            sage: QQ(7/5).absolute_norm()
            7/5
        """
        return self

    def trace(self):
        r"""
        Return the trace from `\QQ` to `\QQ` of `x` (which is just `x`). This
        was added for compatibility with :class:`NumberFields`.

        OUTPUT: ``Rational`` -- reference to ``self``

        EXAMPLES::

            sage: (1/3).trace()
             1/3

        AUTHORS:

        - Craig Citro
        """
        return self

    def charpoly(self, var='x'):
        """
        Return the characteristic polynomial of this rational number. This
        will always be just ``var - self``; this is really here so that code
        written for number fields won't crash when applied to rational
        numbers.

        INPUT:

        - ``var`` -- string

        OUTPUT: polynomial

        EXAMPLES::

            sage: (1/3).charpoly('x')
             x - 1/3

        The default is ``var='x'``. (:issue:`20967`)::

            sage: a = QQ(2); a.charpoly('x')
            x - 2


        AUTHORS:

        - Craig Citro
        """
        QQ = self.parent()
        return QQ[var]([-self,1])

    def minpoly(self, var='x'):
        """
        Return the minimal polynomial of this rational number. This will
        always be just ``x - self``; this is really here so that code written
        for number fields won't crash when applied to rational numbers.

        INPUT:

        - ``var`` -- string

        OUTPUT: polynomial

        EXAMPLES::

            sage: (1/3).minpoly()
            x - 1/3
            sage: (1/3).minpoly('y')
            y - 1/3

        AUTHORS:

        - Craig Citro
        """
        QQ = self.parent()
        return QQ[var]([-self,1])

    def _integer_(self, Z=None):
        """
        Return ``self`` coerced to an integer. Of course this rational number
        must have a denominator of 1.

        OUTPUT: integer

        EXAMPLES::

            sage: (-4/17)._integer_()
            Traceback (most recent call last):
            ...
            TypeError: no conversion of this rational to integer
            sage: (-4/1)._integer_()
            -4
        """
        if not mpz_cmp_si(mpq_denref(self.value), 1) == 0:
            raise TypeError("no conversion of this rational to integer")
        cdef Integer n = Integer.__new__(Integer)
        n.set_from_mpz(mpq_numref(self.value))
        return n

    def numerator(self):
        """
        Return the numerator of this rational number.
        :meth:`numer` is an alias of :meth:`numerator`.

        EXAMPLES::

            sage: x = 5/11
            sage: x.numerator()
            5

            sage: x = 9/3
            sage: x.numerator()
            3

            sage: x = -5/11
            sage: x.numer()
            -5
        """
        cdef Integer n = Integer.__new__(Integer)
        n.set_from_mpz(mpq_numref(self.value))
        return n

    #Define an alias for numerator
    numer = numerator

    def __int__(self):
        """
        Convert this rational to a Python ``int``.

        This truncates ``self`` if ``self`` has a denominator (which is
        consistent with Python's ``int(floats)``).

        EXAMPLES::

            sage: int(7/1)
            7
            sage: int(7/2)
            3
        """
        cdef mpz_t x
        if mpz_cmp_si(mpq_denref(self.value),1) != 0:
            mpz_init(x)
            mpz_tdiv_q(x, mpq_numref(self.value), mpq_denref(self.value))
            n = mpz_get_pylong(x)
            mpz_clear(x)
            return n
        else:
            return mpz_get_pylong(mpq_numref(self.value))

    def denominator(self):
        """
        Return the denominator of this rational number.
        :meth:`denom` is an alias of :meth:`denominator`.

        EXAMPLES::

            sage: x = -5/11
            sage: x.denominator()
            11

            sage: x = 9/3
            sage: x.denominator()
            1

            sage: x = 5/13
            sage: x.denom()
            13
        """
        cdef Integer n = Integer.__new__(Integer)
        n.set_from_mpz(mpq_denref(self.value))
        return n

    #Define an alias for denominator
    denom = denominator

    def as_integer_ratio(self):
        """
        Return the pair ``(self.numerator(), self.denominator())``.

        EXAMPLES::

            sage: x = -12/29
            sage: x.as_integer_ratio()
            (-12, 29)
        """
        n = <Integer>Integer.__new__(Integer)
        d = <Integer>Integer.__new__(Integer)
        n.set_from_mpz(mpq_numref(self.value))
        d.set_from_mpz(mpq_denref(self.value))
        return (n, d)

    def factor(self):
        """
        Return the factorization of this rational number.

        OUTPUT: factorization

        EXAMPLES::

            sage: (-4/17).factor()
            -1 * 2^2 * 17^-1

        Trying to factor 0 gives an arithmetic error::

            sage: (0/1).factor()
            Traceback (most recent call last):
            ...
            ArithmeticError: factorization of 0 is not defined
        """
        from sage.structure.factorization import Factorization
        return self.numerator().factor() * \
           Factorization([(p, -e) for p, e in self.denominator().factor()])

    def support(self):
        """
        Return a sorted list of the primes where this rational number has
        nonzero valuation.

        OUTPUT: the set of primes appearing in the factorization of this
        rational with nonzero exponent, as a sorted list.

        EXAMPLES::

            sage: (-4/17).support()
            [2, 17]

        Trying to find the support of 0 gives an arithmetic error::

            sage: (0/1).support()
            Traceback (most recent call last):
            ...
            ArithmeticError: Support of 0 not defined.
        """
        if self.is_zero():
            raise ArithmeticError("Support of 0 not defined.")
        from sage.arith.misc import prime_factors

        return prime_factors(self)

    def log(self, m=None, prec=None):
        r"""
        Return the log of ``self``.

        INPUT:

        - ``m`` -- the base (default: natural log base e)

        - ``prec`` -- integer (optional); the precision in bits

        OUTPUT:

        When ``prec`` is not given, the log as an element in symbolic
        ring unless the logarithm is exact. Otherwise the log is a
        :class:`RealField` approximation to ``prec`` bit precision.

        EXAMPLES::

            sage: (124/345).log(5)                                                      # needs sage.symbolic
            log(124/345)/log(5)
            sage: (124/345).log(5, 100)                                                 # needs sage.rings.real_mpfr
            -0.63578895682825611710391773754
            sage: log(QQ(125))                                                          # needs sage.symbolic
            3*log(5)
            sage: log(QQ(125), 5)
            3
            sage: log(QQ(125), 3)                                                       # needs sage.symbolic
            3*log(5)/log(3)
            sage: QQ(8).log(1/2)
            -3
            sage: (1/8).log(1/2)
            3
            sage: (1/2).log(1/8)
            1/3
            sage: (1/2).log(8)
            -1/3
            sage: (16/81).log(8/27)                                                     # needs sage.libs.pari
            4/3
            sage: (8/27).log(16/81)                                                     # needs sage.libs.pari
            3/4
            sage: log(27/8, 16/81)                                                      # needs sage.libs.pari
            -3/4
            sage: log(16/81, 27/8)                                                      # needs sage.libs.pari
            -4/3
            sage: (125/8).log(5/2)                                                      # needs sage.libs.pari
            3
            sage: (125/8).log(5/2, prec=53)                                             # needs sage.rings.real_mpfr
            3.00000000000000

        TESTS::

            sage: (25/2).log(5/2)                                                       # needs sage.symbolic
            log(25/2)/log(5/2)
            sage: (-1/2).log(3)                                                         # needs sage.symbolic
            (I*pi + log(1/2))/log(3)
        """
        from sage.rings.integer_ring import ZZ

        cdef int self_sgn
        if self.denom().is_one():
            return ZZ(self.numer()).log(m, prec)
        if m is not None and m <= 0:
            raise ValueError("log base must be positive")
        self_sgn = mpz_sgn(mpq_numref(self.value))
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

        from sage.functions.log import function_log
        if m is None:
            return function_log(self, dont_call_method_on_arg=True)

        anum = self.numer()
        aden = self.denom()
        mrat = Rational(m)
        bnum = mrat.numer()
        bden = mrat.denom()

        anp = anum.perfect_power()
        bnp = bnum.perfect_power()
        adp = aden.perfect_power()
        bdp = bden.perfect_power()

        if anum.is_one():
            a_exp=adp[1]
            a_base=1/adp[0]
        # we already know that aden!=0
        else:
            a_exp=anp[1].gcd(adp[1])
            a_base=(anp[0]**(anp[1]//a_exp))/(adp[0]**(adp[1]//a_exp))

        if bnum.is_one():
            b_exp=bdp[1]
            b_base=1/bdp[0]
        elif bden.is_one():
            b_exp=bnp[1]
            b_base=bnp[0]
        else:
            b_exp=bnp[1].gcd(bdp[1])
            b_base=(bnp[0]**(bnp[1]//b_exp))/(bdp[0]**(bdp[1]//b_exp))

        if a_base == b_base:
            return a_exp/b_exp
        elif a_base*b_base == 1:
            return -a_exp/b_exp

        return (function_log(self, dont_call_method_on_arg=True) /
                function_log(m, dont_call_method_on_arg=True))

    def gamma(self, *, prec=None):
        """
        Return the gamma function evaluated at ``self``. This value is exact
        for integers and half-integers, and returns a symbolic value
        otherwise.  For a numerical approximation, use keyword ``prec``.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: gamma(1/2)
            sqrt(pi)
            sage: gamma(7/2)
            15/8*sqrt(pi)
            sage: gamma(-3/2)
            4/3*sqrt(pi)
            sage: gamma(6/1)
            120
            sage: gamma(1/3)
            gamma(1/3)

        This function accepts an optional precision argument::

            sage: (1/3).gamma(prec=100)                                                 # needs sage.rings.real_mpfr
            2.6789385347077476336556929410
            sage: (1/2).gamma(prec=100)                                                 # needs sage.rings.real_mpfr
            1.7724538509055160272981674833

        TESTS:

        This is not the incomplete gamma function! ::

            sage: (1/2).gamma(5)
            Traceback (most recent call last):
            ...
            TypeError: ...gamma() takes exactly 0 positional arguments (1 given)
        """
        if prec:
            return self.n(prec).gamma()
        else:
            if mpz_cmp_ui(mpq_denref(self.value), 1) == 0:
                return integer.Integer(self).gamma()
            elif mpz_cmp_ui(mpq_denref(self.value), 2) == 0:
                numer = self.numer()
                rat_part = Rational((numer-2).multifactorial(2)) >> ((numer-1)//2)
                from sage.symbolic.constants import pi
                from sage.misc.functional import sqrt
                return sqrt(pi) * rat_part
            else:
                from sage.symbolic.ring import SR
                return SR(self).gamma()

    def floor(self):
        """
        Return the floor of this rational number as an integer.

        OUTPUT: integer

        EXAMPLES::

            sage: n = 5/3; n.floor()
            1
            sage: n = -17/19; n.floor()
            -1
            sage: n = -7/2; n.floor()
            -4
            sage: n = 7/2; n.floor()
            3
            sage: n = 10/2; n.floor()
            5
        """
        cdef integer.Integer n
        n = integer.Integer()
        mpz_fdiv_q(n.value, mpq_numref(self.value), mpq_denref(self.value))
        return n

    def ceil(self):
        """
        Return the ceiling of this rational number.

        OUTPUT: integer

        If this rational number is an integer, this returns this number,
        otherwise it returns the floor of this number +1.

        EXAMPLES::

            sage: n = 5/3; n.ceil()
            2
            sage: n = -17/19; n.ceil()
            0
            sage: n = -7/2; n.ceil()
            -3
            sage: n = 7/2; n.ceil()
            4
            sage: n = 10/2; n.ceil()
            5
        """
        cdef integer.Integer n
        n = integer.Integer()
        mpz_cdiv_q(n.value, mpq_numref(self.value), mpq_denref(self.value))
        return n

    def trunc(self):
        """
        Round this rational number to the nearest integer toward zero.

        EXAMPLES::

            sage: (5/3).trunc()
            1
            sage: (-5/3).trunc()
            -1
            sage: QQ(42).trunc()
            42
            sage: QQ(-42).trunc()
            -42
        """
        cdef integer.Integer n
        n = integer.Integer()
        mpz_tdiv_q(n.value, mpq_numref(self.value), mpq_denref(self.value))
        return n

    def round(Rational self, mode=None):
        """
        Return the nearest integer to ``self``, rounding away by default.
        Deprecation: in the future the default will be changed to rounding to
        even, for consistency with the builtin Python :func:`round`.

        INPUT:

        - ``self`` -- a rational number

        - ``mode`` -- a rounding mode for half integers:

           - ``'toward'`` rounds toward zero
           - ``'away'`` (default) rounds away from zero
           - ``'up'`` rounds up
           - ``'down'`` rounds down
           - ``'even'`` rounds toward the even integer
           - ``'odd'`` rounds toward the odd integer

        OUTPUT: integer

        EXAMPLES::

            sage: (9/2).round()
            doctest:...: DeprecationWarning: the default rounding for rationals, currently `away`, will be changed to `even`.
            See https://github.com/sagemath/sage/issues/35473 for details.
            5
            sage: n = 4/3; n.round()
            1
            sage: n = -17/4; n.round()
            -4
            sage: n = -5/2; n.round()
            -3
            sage: n.round("away")
            -3
            sage: n.round("up")
            -2
            sage: n.round("down")
            -3
            sage: n.round("even")
            -2
            sage: n.round("odd")
            -3
        """
        if mode is None:
            if self.denominator() == 2:
                from sage.misc.superseded import deprecation
                deprecation(35473,
                    "the default rounding for rationals, currently `away`, will be changed to `even`.")
            mode = "away"
        if not (mode in ['toward', 'away', 'up', 'down', 'even', 'odd']):
            raise ValueError("rounding mode must be one of 'toward', 'away', 'up', 'down', 'even', or 'odd'")
        if self.denominator() == 1:
            from sage.rings.integer import Integer
            return Integer(self)
        if self.denominator() == 2:
            # round down:
            if (mode == "down") or \
                   (mode == "toward" and self > 0) or \
                   (mode == "away" and self < 0) or \
                   (mode == "even" and self.numerator() % 4 == 1) or \
                   (mode == "odd" and self.numerator() % 4 == 3):
                return self.numerator() // self.denominator()
            else:
                return self.numerator() // self.denominator() + 1
        else:
            q, r = self.numerator().quo_rem(self.denominator())
            if r < self.denominator() / 2:
                return q
            else:
                return q+1

    __round__ = round

    def real(self):
        """
        Return the real part of ``self``, which is ``self``.

        EXAMPLES::

            sage: (1/2).real()
            1/2
        """
        return self

    def imag(self):
        """
        Return the imaginary part of ``self``, which is zero.

        EXAMPLES::

            sage: (1/239).imag()
            0
        """
        return self._parent(0)

    def height(self):
        """
        The max absolute value of the numerator and denominator of ``self``, as
        an :class:`Integer`.

        OUTPUT: integer

        EXAMPLES::

            sage: a = 2/3
            sage: a.height()
            3
            sage: a = 34/3
            sage: a.height()
            34
            sage: a = -97/4
            sage: a.height()
            97

        AUTHORS:

        - Naqi Jaffery (2006-03-05): examples

        .. NOTE::

           For the logarithmic height, use :meth:`global_height()`.
        """
        x = abs(self.numer())
        if x > self.denom():
            return x
        return self.denom()

    def _lcm(self, Rational other):
        """
        Return the least common multiple, in the rational numbers, of ``self``
        and ``other``. This function returns either 0 or 1 (as a rational
        number).

        INPUT:

        - ``other`` -- Rational

        OUTPUT: ``Rational`` -- 0 or 1

        EXAMPLES::

            sage: (2/3)._lcm(3/5)
            1
            sage: (0/1)._lcm(0/1)
            0
            sage: type((2/3)._lcm(3/5))
            <class 'sage.rings.rational.Rational'>
        """
        if mpz_cmp_si(mpq_numref(self.value), 0) == 0 and \
               mpz_cmp_si(mpq_numref(other.value), 0) == 0:
            return Rational(0)
        return Rational(1)

    def additive_order(self):
        """
        Return the additive order of ``self``.

        OUTPUT: integer or infinity

        EXAMPLES::

            sage: QQ(0).additive_order()
            1
            sage: QQ(1).additive_order()
            +Infinity
        """
        import sage.rings.infinity
        if self.is_zero():
            return integer.Integer(1)
        else:
            return sage.rings.infinity.infinity

    def multiplicative_order(self):
        """
        Return the multiplicative order of ``self``.

        OUTPUT: integer or ``infinity``

        EXAMPLES::

            sage: QQ(1).multiplicative_order()
            1
            sage: QQ('1/-1').multiplicative_order()
            2
            sage: QQ(0).multiplicative_order()
            +Infinity
            sage: QQ('2/3').multiplicative_order()
            +Infinity
            sage: QQ('1/2').multiplicative_order()
            +Infinity
        """
        import sage.rings.infinity
        if self.is_one():
            return integer.Integer(1)
        elif mpz_cmpabs(mpq_numref(self.value),mpq_denref(self.value)) == 0:
            # if the numerator and the denominator are equal in absolute value,
            # then the rational number is -1
            return integer.Integer(2)
        else:
            return sage.rings.infinity.infinity

    def is_one(self):
        r"""
        Determine if a rational number is one.

        OUTPUT: boolean

        EXAMPLES::

            sage: QQ(1/2).is_one()
            False
            sage: QQ(4/4).is_one()
            True
        """
        # A rational number is equal to 1 iff its numerator and denominator are equal
        return mpz_cmp(mpq_numref(self.value),mpq_denref(self.value)) == 0

    def is_integral(self):
        r"""
        Determine if a rational number is integral (i.e., is in
        `\ZZ`).

        OUTPUT: boolean

        EXAMPLES::

            sage: QQ(1/2).is_integral()
            False
            sage: QQ(4/4).is_integral()
            True
        """
        return mpz_cmp_si(mpq_denref(self.value), 1) == 0

    def is_rational(self):
        r"""
        Return ``True`` since this is a rational number.

        EXAMPLES::

            sage: (3/4).is_rational()
            True
        """
        return True

    # Function alias for checking if the number is a integer.  Added to solve issue 15500
    is_integer = is_integral

    def is_S_integral(self, S=[]):
        r"""
        Determine if the rational number is ``S``-integral.

        ``x`` is ``S``-integral if ``x.valuation(p)>=0`` for all ``p`` not in
        ``S``, i.e., the denominator of ``x`` is divisible only by the primes
        in ``S``.

        INPUT:

        - ``S`` -- list or tuple of primes

        OUTPUT: boolean

        .. NOTE::

           Primality of the entries in ``S`` is not checked.

        EXAMPLES::

            sage: QQ(1/2).is_S_integral()
            False
            sage: QQ(1/2).is_S_integral([2])
            True
            sage: [a for a in range(1,11) if QQ(101/a).is_S_integral([2,5])]
            [1, 2, 4, 5, 8, 10]
        """
        if self.is_integral():
            return True
        return self.prime_to_S_part(S).is_integral()

    def is_S_unit(self, S=None):
        r"""
        Determine if the rational number is an ``S``-unit.

        ``x`` is an ``S``-unit if ``x.valuation(p)==0`` for all ``p`` not in
        ``S``, i.e., the numerator and denominator of ``x`` are divisible only
        by the primes in `S`.

        INPUT:

        - ``S`` -- list or tuple of primes

        OUTPUT: boolean

        .. NOTE::

           Primality of the entries in ``S`` is not checked.

        EXAMPLES::

            sage: QQ(1/2).is_S_unit()
            False
            sage: QQ(1/2).is_S_unit([2])
            True
            sage: [a for a in range(1,11) if QQ(10/a).is_S_unit([2,5])]
            [1, 2, 4, 5, 8, 10]
        """
        a = self.abs()
        if a == 1:
            return True
        if S is None:
            return False
        return a.prime_to_S_part(S) == 1

    cdef _lshift(self, long int exp):
        r"""
        Return ``self * 2^exp``.
        """
        cdef Rational x
        x = <Rational> Rational.__new__(Rational)
        sig_on()
        if exp < 0:
            mpq_div_2exp(x.value,self.value,-exp)
        else:
            mpq_mul_2exp(x.value,self.value,exp)
        sig_off()
        return x

    def __lshift__(x, y):
        """
        Left shift operator ``x << y``.

        INPUT:

        - ``x``, ``y`` -- integer or rational

        OUTPUT: rational

        EXAMPLES::

            sage: (2/3).__lshift__(4/1)
            32/3
            sage: (2/3).__lshift__(4/7)
            Traceback (most recent call last):
            ...
            ValueError: denominator must be 1
            sage: (2).__lshift__(4/1)
            32
            sage: (2/3).__lshift__(4)
            32/3
            sage: (2/3) << (4/1)
            32/3
        """
        if isinstance(x, Rational):
            if isinstance(y, (int, integer.Integer)):
                return (<Rational>x)._lshift(y)
            if isinstance(y, Rational):
                if mpz_cmp_si(mpq_denref((<Rational>y).value), 1) != 0:
                    raise ValueError("denominator must be 1")
                return (<Rational>x)._lshift(y)
        return coercion_model.bin_op(x, y, operator.lshift)

    cdef _rshift(self, long int exp):
        r"""
        Return ``self / 2^exp``.
        """
        cdef Rational x
        x = <Rational> Rational.__new__(Rational)
        sig_on()
        if exp < 0:
            mpq_mul_2exp(x.value,self.value,-exp)
        else:
            mpq_div_2exp(x.value,self.value,exp)
        sig_off()
        return x

    def __rshift__(x, y):
        """
        Right shift operator ``x >> y``.

        INPUT:

        - ``x``, ``y`` -- integer or rational

        OUTPUT: rational

        EXAMPLES::

            sage: (2/3).__rshift__(4/1)
            1/24
            sage: (2/3).__rshift__(4/7)
            Traceback (most recent call last):
            ...
            ValueError: denominator must be 1
            sage: (2).__rshift__(4/1)
            0
            sage: (2/1).__rshift__(4)
            1/8
            sage: (2/1) >>(4/1)
            1/8
        """
        if isinstance(x, Rational):
            if isinstance(y, (int, integer.Integer)):
                return (<Rational>x)._rshift(y)
            if isinstance(y, Rational):
                if mpz_cmp_si(mpq_denref((<Rational>y).value), 1) != 0:
                    raise ValueError("denominator must be 1")
                return (<Rational>x)._rshift(y)
        return coercion_model.bin_op(x, y, operator.rshift)

    def conjugate(self):
        """
        Return the complex conjugate of this rational number, which is
        the number itself.

        EXAMPLES::

            sage: n = 23/11
            sage: n.conjugate()
            23/11
        """
        return self

    ##################################################
    # Support for interfaces
    ##################################################

    def __pari__(self):
        """
        Return the PARI version of this rational number.

        EXAMPLES::

            sage: n = 9390823/17
            sage: m = n.__pari__(); m                                                   # needs sage.libs.pari
            9390823/17
            sage: type(m)                                                               # needs sage.libs.pari
            <class 'cypari2.gen.Gen'>
            sage: m.type()                                                              # needs sage.libs.pari
            't_FRAC'
        """
        global new_gen_from_rational
        if new_gen_from_rational is None:
            from sage.libs.pari.convert_sage import new_gen_from_rational
        return new_gen_from_rational(self)

    def _interface_init_(self, I=None):
        """
        Return representation of this rational suitable for coercing into
        almost any computer algebra system.

        OUTPUT: string

        EXAMPLES::

            sage: (2/3)._interface_init_()
            '2/3'
            sage: kash(3/1).Type()              # optional - kash
            elt-fld^rat
            sage: magma(3/1).Type()             # optional - magma
            FldRatElt
        """
        return '%s/%s'%(self.numerator(), self.denominator())

    def _sage_input_(self, sib, coerced):
        r"""
        Produce an expression which will reproduce this value when evaluated.

        EXAMPLES::

            sage: sage_input(QQ(1), verify=True)
            # Verified
            QQ(1)
            sage: sage_input(-22/7, verify=True)
            # Verified
            -22/7
            sage: sage_input(-22/7, preparse=False)
            -ZZ(22)/7
            sage: sage_input(10^-50, verify=True)
            # Verified
            1/100000000000000000000000000000000000000000000000000
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: (-2/37)._sage_input_(SageInputBuilder(preparse=False), False)
            {unop:- {binop:/ {call: {atomic:ZZ}({atomic:2})} {atomic:37}}}
            sage: QQ(5)._sage_input_(SageInputBuilder(preparse=False), True)
            {atomic:5}
        """

        # This code is extensively described in the docstring
        # for sage_input.py.

        num = self.numerator()
        neg = (num < 0)
        if neg: num = -num
        if self.denominator() == 1:
            if coerced:
                v = sib.int(num)
            else:
                v = sib.name('QQ')(sib.int(num))
        else:
            v = sib(num)/sib.int(self.denominator())
        if neg: v = -v
        return v


# The except value is just some random double, it doesn't matter what it is.
cdef double mpq_get_d_nearest(mpq_t x) except? -648555075988944.5:
    """
    Convert a ``mpq_t`` to a ``double``, with round-to-nearest-even.
    This differs from ``mpq_get_d()`` which does round-to-zero.

    TESTS::

        sage: q = QQ(); float(q)
        0.0
        sage: q = 2^-10000; float(q)
        0.0
        sage: float(-q)
        -0.0
        sage: q = 2^10000/1; float(q)
        inf
        sage: float(-q)
        -inf

    ::

        sage: q = 2^-1075; float(q)
        0.0
        sage: float(-q)
        -0.0
        sage: q = 2^52 / 2^1074; float(q)  # Smallest normal double
        2.2250738585072014e-308
        sage: float(-q)
        -2.2250738585072014e-308
        sage: q = (2^52 + 1/2) / 2^1074; float(q)
        2.2250738585072014e-308
        sage: float(-q)
        -2.2250738585072014e-308
        sage: q = (2^52 + 1) / 2^1074; float(q)  # Next normal double
        2.225073858507202e-308
        sage: float(-q)
        -2.225073858507202e-308
        sage: q = (2^52 - 1) / 2^1074; float(q)  # Largest denormal double
        2.225073858507201e-308
        sage: float(-q)
        -2.225073858507201e-308
        sage: q = 1 / 2^1074; float(q)  # Smallest denormal double
        5e-324
        sage: float(-q)
        -5e-324
        sage: q = (1/2) / 2^1074; float(q)
        0.0
        sage: float(-q)
        -0.0
        sage: q = (3/2) / 2^1074; float(q)
        1e-323
        sage: float(-q)
        -1e-323
        sage: q = (2/3) / 2^1074; float(q)
        5e-324
        sage: float(-q)
        -5e-324
        sage: q = (1/3) / 2^1074; float(q)
        0.0
        sage: float(-q)
        -0.0
        sage: q = (2^53 - 1) * 2^971/1; float(q)  # Largest double
        1.7976931348623157e+308
        sage: float(-q)
        -1.7976931348623157e+308
        sage: q = (2^53) * 2^971/1; float(q)
        inf
        sage: float(-q)
        -inf
        sage: q = (2^53 - 1/2) * 2^971/1; float(q)
        inf
        sage: float(-q)
        -inf
        sage: q = (2^53 - 2/3) * 2^971/1; float(q)
        1.7976931348623157e+308
        sage: float(-q)
        -1.7976931348623157e+308

    AUTHORS:

    - Paul Zimmermann, Jeroen Demeyer (:issue:`14416`)
    """
    cdef mpz_ptr a = mpq_numref(x)
    cdef mpz_ptr b = mpq_denref(x)
    cdef int resultsign = mpz_sgn(a)

    if resultsign == 0:
        return 0.0

    cdef Py_ssize_t sa = mpz_sizeinbase(a, 2)
    cdef Py_ssize_t sb = mpz_sizeinbase(b, 2)

    # Easy case: both numerator and denominator are exactly
    # representable as doubles.
    if sa <= 53 and sb <= 53:
        return mpz_get_d(a) / mpz_get_d(b)

    # General case

    # We should shift a right by this amount
    cdef Py_ssize_t shift = sa - sb - 54

    # At this point, we know that q0 = a/b / 2^shift satisfies
    # 2^53 < q0 < 2^55.
    # The end result d = q0 * 2^shift (rounded).

    # Check for obvious overflow/underflow before shifting
    if shift <= -1130:  # |d| < 2^-1075
        if resultsign < 0:
            return -0.0
        else:
            return 0.0
    elif shift >= 971:  # |d| > 2^1024
        if resultsign < 0:
            return float('-inf')
        else:
            return float('inf')

    sig_on()

    # Compute q = trunc(a / 2^shift) and let remainder_is_zero be True
    # if and only if no truncation occurred.
    cdef mpz_t q, r
    mpz_init(q)
    mpz_init(r)
    cdef int remainder_is_zero
    if shift > 0:
        remainder_is_zero = mpz_divisible_2exp_p(a, shift)
        mpz_tdiv_q_2exp(q, a, shift)
    else:
        mpz_mul_2exp(q, a, -shift)
        remainder_is_zero = True

    # Now divide by b to get q = trunc(a/b / 2^shift).
    # remainder_is_zero is True if and only if no truncation occurred
    # (in neither division).
    mpz_tdiv_qr(q, r, q, b)
    if remainder_is_zero:
        remainder_is_zero = (mpz_cmp_ui(r, 0) == 0)

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
    mpz_clear(r)
    sig_off()

    # The quotient q64 has 54 or 55 bits, but we need exactly 54.
    # Shift it down by 1 one if needed.
    cdef Py_ssize_t add_shift
    if q64 < (1ULL << 54):
        add_shift = 0
    else:
        add_shift = 1

    if (shift + add_shift) < -1075:
        # The result will be denormal, ensure the final shift is -1075
        # to avoid a double rounding.
        add_shift = -1075 - shift

    # Add add_shift to shift and let q = trunc(a/b / 2^shift)
    # for the new shift value.
    cdef uint64_t mask
    if add_shift:
        assert add_shift > 0
        assert add_shift < 64
        shift += add_shift
        # We do an additional division of q by 2^add_shift.
        if remainder_is_zero:
            mask = ((1ULL << add_shift)-1)
            remainder_is_zero = ((q64 & mask) == 0)
        q64 = q64 >> add_shift

    # Round q64 from 54 to 53 bits of precision.
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
    # This is because q64 is even and satisfies q64 <= 2^54,
    # (with 2^53 <= q64 <= 2^54 unless in the denormal case).
    cdef double d = <double>q64
    if resultsign < 0:
        d = -d
    return ldexp(d, shift)


@cython.binding(True)
def make_rational(s):
    """
    Make a rational number from ``s`` (a string in base 32).

    INPUT:

    - ``s`` -- string in base 32

    OUTPUT: rational

    EXAMPLES::

        sage: (-7/15).str(32)
        '-7/f'
        sage: sage.rings.rational.make_rational('-7/f')
        -7/15
    """
    r = Rational()
    mpq_set_str(r.value, str_to_bytes(s), 32)
    return r


cdef class Z_to_Q(Morphism):
    r"""
    A morphism from `\ZZ` to `\QQ`.
    """

    def __init__(self):
        """
        Create morphism from integers to rationals.

        EXAMPLES::

            sage: sage.rings.rational.Z_to_Q()
            Natural morphism:
              From: Integer Ring
              To:   Rational Field
        """
        from sage.rings.integer_ring import ZZ

        import sage.categories.homset
        Morphism.__init__(self, sage.categories.homset.Hom(ZZ, QQ))

    cpdef Element _call_(self, x):
        """
        Return the image of the morphism on ``x``.

        EXAMPLES::

            sage: sage.rings.rational.Z_to_Q()(2) # indirect doctest
            2
        """
        cdef Rational rat
        rat = <Rational> Rational.__new__(Rational)
        mpq_set_z(rat.value, (<integer.Integer>x).value)
        return rat

    def _repr_type(self):
        """
        Return string that describes the type of morphism.

        EXAMPLES::

            sage: sage.rings.rational.Z_to_Q()._repr_type()
            'Natural'
        """
        return "Natural"

    def section(self):
        """
        Return a section of this morphism.

        EXAMPLES::

            sage: f = QQ.coerce_map_from(ZZ).section(); f
            Generic map:
              From: Rational Field
              To:   Integer Ring

        This map is a morphism in the category of sets with partial
        maps (see :issue:`15618`)::

            sage: f.parent()
            Set of Morphisms from Rational Field to Integer Ring
             in Category of sets with partial maps
        """
        from sage.categories.sets_with_partial_maps import SetsWithPartialMaps
        return Q_to_Z(self._codomain.Hom(self.domain(), category=SetsWithPartialMaps()))

    def is_surjective(self):
        r"""
        Return whether this morphism is surjective.

        EXAMPLES::

            sage: QQ.coerce_map_from(ZZ).is_surjective()
            False
        """
        return False

cdef class Q_to_Z(Map):
    r"""
    A morphism from `\QQ` to `\ZZ`.

    TESTS::

        sage: type(ZZ.convert_map_from(QQ))
        <class 'sage.rings.rational.Q_to_Z'>
    """
    cpdef Element _call_(self, x):
        """
        A fast map from the rationals to the integers.

        EXAMPLES::

            sage: f = sage.rings.rational.Q_to_Z(QQ, ZZ)
            sage: f(1/2) # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: no conversion of this rational to integer
            sage: f(4/2) # indirect doctest
            2
        """
        if not mpz_cmp_si(mpq_denref((<Rational>x).value), 1) == 0:
            raise TypeError("no conversion of this rational to integer")
        cdef Integer n = Integer.__new__(Integer)
        n.set_from_mpz(mpq_numref((<Rational>x).value))
        return n

    def section(self):
        """
        Return a section of this morphism.

        EXAMPLES::

            sage: sage.rings.rational.Q_to_Z(QQ, ZZ).section()
            Natural morphism:
              From: Integer Ring
              To:   Rational Field
        """
        return Z_to_Q()


cdef class int_to_Q(Morphism):
    r"""
    A morphism from Python 3 ``int`` to `\QQ`.
    """
    def __init__(self):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: sage.rings.rational.int_to_Q()
            Native morphism:
              From: Set of Python objects of class 'int'
              To:   Rational Field
        """
        from sage.rings import rational_field
        import sage.categories.homset
        from sage.sets.pythonclass import Set_PythonType
        Morphism.__init__(self, sage.categories.homset.Hom(
            Set_PythonType(int), rational_field.QQ))

    cpdef Element _call_(self, a):
        """
        Return the image of the morphism on ``a``.

        EXAMPLES::

            sage: f = sage.rings.rational.int_to_Q()
            sage: f(4^100)
            1606938044258990275541962092341162602522202993782792835301376
        """

        cdef Rational rat
        cdef long a_long
        cdef int err = 0

        rat = <Rational> Rational.__new__(Rational)

        integer_check_long_py(a, &a_long, &err)

        if not err:
            mpq_set_si(rat.value, a_long, 1)
        else:
            mpz_set_pylong(mpq_numref(rat.value), a)

        return rat

    def _repr_type(self):
        """
        Return string that describes the type of morphism.

        EXAMPLES::

            sage: sage.rings.rational.int_to_Q()._repr_type()
            'Native'
        """
        return "Native"


# Support Python's numbers abstract base class
import numbers
numbers.Rational.register(Rational)
