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
  multiplicative_order, is_one; optimized __nonzero__ ; documented:
  lcm,gcd

- John Cremona (2009-05-15): added support for local and global
  logarithmic heights.

- Travis Scrimshaw (2012-10-18): Added doctests for full coverage.

TESTS::

    sage: a = -2/3
    sage: a == loads(dumps(a))
    True
"""

#*****************************************************************************
#       Copyright (C) 2004, 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


include "sage/ext/interrupt.pxi"  # ctrl-c interrupt block support
include "sage/ext/gmp.pxi"
include "sage/ext/stdsage.pxi"
include "sage/ext/python.pxi"
include "sage/libs/pari/decl.pxi"


import sys
import operator

from sage.misc.mathml import mathml

import sage.misc.misc as misc
import sage.rings.rational_field

cimport integer
import integer
from integer cimport Integer

import sage.libs.pari.pari_instance
from sage.libs.pari.gen cimport gen as pari_gen
from sage.libs.pari.pari_instance cimport PariInstance

from integer_ring import ZZ

from sage.structure.element cimport Element, RingElement, ModuleElement
from sage.structure.element import bin_op
from sage.categories.morphism cimport Morphism
from sage.categories.map cimport Map

import sage.structure.factorization

import sage.rings.real_mpfr
import sage.rings.real_double
from libc.stdint cimport uint64_t

cimport sage.rings.fast_arith
import  sage.rings.fast_arith

cdef sage.rings.fast_arith.arith_int ai
ai = sage.rings.fast_arith.arith_int()

cdef object numpy_long_interface = {'typestr': '=i4' if sizeof(long) == 4 else '=i8' }
cdef object numpy_int64_interface = {'typestr': '=i8'}
cdef object numpy_object_interface = {'typestr': '|O'}
cdef object numpy_double_interface = {'typestr': '=f8'}

cdef extern from "convert.h":
    ctypedef long* GEN
    void t_FRAC_to_QQ ( mpq_t value, GEN g )
    void t_INT_to_ZZ ( mpz_t value, GEN g )


cdef extern from "mpz_pylong.h":
    cdef mpz_get_pylong(mpz_t src)
    cdef int mpz_set_pylong(mpz_t dst, src) except -1
    cdef long mpz_pythonhash(mpz_t src)

cdef class Rational(sage.structure.element.FieldElement)

cdef inline void set_from_mpq(Rational self, mpq_t value):
    mpq_set(self.value, value)

cdef inline void set_from_Rational(Rational self, Rational other):
    mpq_set(self.value, other.value)

cdef inline void set_from_Integer(Rational self, integer.Integer other):
    mpq_set_z(self.value, other.value)

cdef object Rational_mul_(Rational a, Rational b):
    cdef Rational x
    x = <Rational> PY_NEW(Rational)

    sig_on()
    mpq_mul(x.value, a.value, b.value)
    sig_off()

    return x

cdef object Rational_div_(Rational a, Rational b):
    cdef Rational x
    x = <Rational> PY_NEW(Rational)

    sig_on()
    mpq_div(x.value, a.value, b.value)
    sig_off()

    return x

cdef Rational_add_(Rational self, Rational other):
    cdef Rational x
    x = <Rational> PY_NEW(Rational)
    sig_on()
    mpq_add(x.value, self.value, other.value)
    sig_off()
    return x

cdef Rational_sub_(Rational self, Rational other):
    cdef Rational x
    x = <Rational> PY_NEW(Rational)

    sig_on()
    mpq_sub(x.value, self.value, other.value)
    sig_off()

    return x

cdef object the_rational_ring
the_rational_ring = sage.rings.rational_field.Q

# make sure zero/one elements are set
cdef set_zero_one_elements():
    global the_rational_ring
    the_rational_ring._zero_element = Rational(0)
    the_rational_ring._one_element = Rational(1)

set_zero_one_elements()

cpdef Integer integer_rational_power(Integer a, Rational b):
    """
    Compute `a^b` as an integer, if it is integral, or return None.
    The positive real root is taken for even denominators.

    INPUT::

        a -- an Integer
        b -- a positive Rational

    OUTPUT::

        `a^b` as an ``Integer`` or ``None``

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
    """
    cdef Integer z = <Integer>PY_NEW(Integer)
    if mpz_sgn(mpq_numref(b.value)) < 0:
        raise ValueError, "Only positive exponents supported."
    cdef int sgn = mpz_sgn(a.value)
    cdef bint exact
    if sgn == 0:
        pass # z is 0
    elif sgn < 0:
        return None
    elif mpz_cmp_ui(a.value, 1) == 0:
        mpz_set_ui(z.value, 1)
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

cpdef rational_power_parts(a, b, factor_limit=10**5):
    """
    Compute rationals or integers `c` and `d` such that `a^b = c*d^b`
    with `d` small. This is used for simplifying radicals.

    INPUT::

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

    Check if :trac:`8540` is fixed::

        sage: rational_power_parts(3/4, -1/2)
        (2, 3)
        sage: t = (3/4)^(-1/2); t
        2/3*sqrt(3)
        sage: t^2
        4/3
    """
    b_negative=False
    if b < 0:
        b_negative = True
        b = -b
        a = ~a
    if isinstance(a, Rational):
        c1, d1 = rational_power_parts(a.numerator(), b)
        c2, d2 = rational_power_parts(a.denominator(), b)
        return (c1/c2, d1/d2) if not b_negative else (c1/c2, d2/d1)
    elif not isinstance(a, Integer):
        a = Integer(a)
    c = integer_rational_power(a, b)
    if c is not None:
        return c, 1
    numer, denom = b.numerator(), b.denominator()
    if a < factor_limit*factor_limit:
        f = a.factor()
    else:
        from sage.rings.factorint import factor_trial_division
        f = factor_trial_division(a,factor_limit)
    c = 1
    d = 1
    for p, e in f:
        c *= p**((e // denom)*numer)
        d *= p**(e % denom)
    if a < 0 and numer & 1:
        d = -d
    return (c, d) if not b_negative else (c, ~d)


def is_Rational(x):
    """
    Return true if x is of the Sage rational number type.

    EXAMPLES::

        sage: from sage.rings.rational import is_Rational
        sage: is_Rational(2)
        False
        sage: is_Rational(2/1)
        True
        sage: is_Rational(int(2))
        False
        sage: is_Rational(long(2))
        False
        sage: is_Rational('5')
        False
    """
    return isinstance(x, Rational)


cdef class Rational(sage.structure.element.FieldElement):
    """
    A rational number.

    Rational numbers are implemented using the GMP C library.

    EXAMPLES::

        sage: a = -2/3
        sage: type(a)
        <type 'sage.rings.rational.Rational'>
        sage: parent(a)
        Rational Field
        sage: Rational('1/0')
        Traceback (most recent call last):
        ...
        TypeError: unable to convert 1/0 to a rational
        sage: Rational(1.5)
        3/2
        sage: Rational('9/6')
        3/2
        sage: Rational((2^99,2^100))
        1/2
        sage: Rational(("2", "10"), 16)
        1/8
        sage: Rational(QQbar(125/8).nth_root(3))
        5/2
        sage: Rational(AA(209735/343 - 17910/49*golden_ratio).nth_root(3) + 3*AA(golden_ratio))
        53/7
        sage: QQ(float(1.5))
        3/2
        sage: QQ(RDF(1.2))
        6/5

    Conversion from PARI::

        sage: Rational(pari('-939082/3992923'))
        -939082/3992923
        sage: Rational(pari('Pol([-1/2])'))  #9595
        -1/2
    """
    def __cinit__(self):
        r"""
        Initialize ``self`` as an element of `\QQ`.

        EXAMPLES::

            sage: p = Rational(3) # indirect doctest
            sage: p.parent()
            Rational Field
        """
        global the_rational_ring
        mpq_init(self.value)
        self._parent = the_rational_ring

    def __init__(self, x=None, unsigned int base=0):
        """
        Create a new rational number.

        INPUT:

        -  ``x`` - object (default: ``None``)

        -  ``base`` - base if ``x`` is a string

        EXAMPLES::

            sage: a = Rational()
            sage: a.__init__(7); a
            7
            sage: a.__init__('70', base=8); a
            56
            sage: a.__init__(pari('2/3')); a
            2/3
            sage: a.__init__('-h/3ki', 32); a
            -17/3730

        .. NOTE::

           This is for demonstration purposes only, mutating rationals
           is almost always the wrong thing to do.
        """
        if x is not None:
            self.__set_value(x, base)

    def __reduce__(self):
        """
        Used in pickling rational numbers.

        EXAMPLES::

            sage: a = 3/5
            sage: a.__reduce__()
            (<built-in function make_rational>, ('3/5',))
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
            TypeError: rational is not an integer
        """
        if self.denominator() == 1:
            return int(self)
        raise TypeError, "rational is not an integer"

    def _reduce_set(self, s):
        """
        Used in setting a rational number when unpickling. Do not call this
        from external code since it violates immutability.

        INPUT:

        -  ``s`` - string representation of rational in base 32

        EXAMPLES::

            sage: a = -17/3730; _, (s,) = a.__reduce__(); s
            '-h/3ki'
            sage: b = 2/3; b._reduce_set('-h/3ki'); b
            -17/3730

            sage: Rational(pari(-345/7687))
            -345/7687
            sage: Rational(pari(-345))
            -345
            sage: Rational(pari('Mod(2,3)'))
            2
            sage: Rational(pari('x'))
            Traceback (most recent call last):
            ...
            TypeError: Unable to coerce PARI x to an Integer
        """
        mpq_set_str(self.value, s, 32)

    cdef __set_value(self, x, unsigned int base):
        cdef int n
        cdef Rational temp_rational
        cdef integer.Integer a, b

        if isinstance(x, Rational):
            set_from_Rational(self, x)

        elif isinstance(x, int):
            i = x
            mpq_set_si(self.value, i, 1)

        elif isinstance(x, long):
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
                xstr = x.str(base)
                if '.' in xstr:
                    exp = (len(xstr) - (xstr.index('.') +1))
                    p = base**exp
                    pstr = '1'+'0'*exp
                    s = xstr.replace('.','') +'/'+pstr
                    n = mpq_set_str( self.value, s, base)
                    if n or mpz_cmp_si(mpq_denref(self.value), 0) == 0:
                        raise TypeError, "unable to convert %s to a rational"%x
                    mpq_canonicalize(self.value)
                else:
                    n = mpq_set_str(self.value, xstr, base)
                    if n or mpz_cmp_si(mpq_denref(self.value), 0) == 0:
                        raise TypeError, "unable to convert %s to a rational"%x
                    mpq_canonicalize(self.value)

        elif isinstance(x, str):
            n = mpq_set_str(self.value, x, base)
            if n or mpz_cmp_si(mpq_denref(self.value), 0) == 0:
                raise TypeError, "unable to convert %s to a rational"%x
            mpq_canonicalize(self.value)

        elif hasattr(x, "_rational_"):
            set_from_Rational(self, x._rational_())

        elif isinstance(x, tuple) and len(x) == 2:
            num = x[0]
            denom = x[1]
            if isinstance(num, int) and isinstance(denom, int):
                mpq_set_si(self.value, num, denom)
            else:
                if not isinstance(num, integer.Integer):
                    num = integer.Integer(num, base)
                if not isinstance(denom, integer.Integer):
                    denom = integer.Integer(denom, base)
                mpz_set(mpq_numref(self.value), (<integer.Integer>num).value)
                mpz_set(mpq_denref(self.value), (<integer.Integer>denom).value)
            if mpz_sgn(mpq_denref(self.value)) == 0:
                raise ValueError, "denominator must not be 0"
            mpq_canonicalize(self.value)

        elif isinstance(x, pari_gen):
            x = x.simplify()
            if typ((<pari_gen>x).g) == t_FRAC:
                t_FRAC_to_QQ(self.value, (<pari_gen>x).g)
            elif typ((<pari_gen>x).g) == t_INT:
                t_INT_to_ZZ(mpq_numref(self.value), (<pari_gen>x).g)
                mpz_set_si(mpq_denref(self.value), 1)
            else:
                a = integer.Integer(x)
                mpz_set(mpq_numref(self.value), a.value)
                mpz_set_si(mpq_denref(self.value), 1)

        elif isinstance(x, list) and len(x) == 1:
            self.__set_value(x[0], base)

        elif hasattr(x, 'rational_reconstruction'):
            temp_rational = x.rational_reconstruction()
            mpq_set(self.value, temp_rational.value)

        elif isinstance(x, (float, sage.rings.real_double.RealDoubleElement)):
            self.__set_value(sage.rings.real_mpfr.RealNumber(sage.rings.real_mpfr.RR, x), base)

        else:

            raise TypeError, "Unable to coerce %s (%s) to Rational"%(x,type(x))

    cdef void set_from_mpq(Rational self, mpq_t value):
        mpq_set(self.value, value)

    def list(self):
        """
        Return a list with the rational element in it, to be compatible
        with the method for number fields.

        OUTPUT:

        -  ``list`` - the list ``[self]``

        EXAMPLES::

            sage: m = 5/3
            sage: m.list()
            [5/3]
        """
        return [ self ]


    def __richcmp__(left, right, int op):
        """
        Rich comparison between two rational numbers.

        INPUT:

        -  ``left, right`` -- objects

        -  ``op`` -- integer

        EXAMPLES::

            sage: 1/3 < 2/3
            True
            sage: 2/3 < 1/3
            False
            sage: 4/5 < 2.0
            True
            sage: 4/5 < 0.8
            False
        """
        return (<sage.structure.element.Element>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, sage.structure.element.Element right) except -2:
        cdef int i
        i = mpq_cmp((<Rational>left).value, (<Rational>right).value)
        if i < 0: return -1
        elif i == 0: return 0
        else: return 1

    def __copy__(self):
        """
        Return a copy of ``self``.

        OUTPUT: Rational

        EXAMPLES::

            sage: a = -17/37
            sage: copy(a) is a
            False

        Coercion does not make a new copy::

            sage: QQ(a) is a
            True

        The constructor also makes a new copy::

            sage: Rational(a) is a
            False
        """
        cdef Rational z
        z = <Rational> PY_NEW(Rational)
        mpq_set(z.value, self.value)
        return z

    def  __dealloc__(self):
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
        else:
            if self < 0:
                return "-\\frac{%s}{%s}"%(-self.numer(), self.denom())
            else:
               return "\\frac{%s}{%s}"%(self.numer(), self.denom())

    def _sympy_(self):
        """
        Convert Sage ``Rational`` to SymPy ``Rational``.

        EXAMPLES::

            sage: n = 1/2; n._sympy_()
            1/2
            sage: n = -1/5; n._sympy_()
            -1/5
            sage: from sympy import Symbol
            sage: QQ(1)+Symbol('x')*QQ(2)
            2*x + 1
        """
        import sympy
        return sympy.Rational(int(self.numerator()), int(self.denominator()))

    def _magma_init_(self, magma):
        """
        Return the magma representation of ``self``.

        EXAMPLES::

            sage: n = -485/82847
            sage: n._magma_init_(magma) # optional - magma
            '-485/82847'
        """
        return self.numerator()._magma_init_(magma) + '/' + self.denominator()._magma_init_(magma)

    property __array_interface__:
        def __get__(self):
            """
            Used for NumPy conversion. If ``self`` is integral, it converts to
            an ``Integer``. Otherwise it converts to a double floating point
            value.

            EXAMPLES::

                sage: import numpy
                sage: numpy.array([1, 2, 3/1])
                array([1, 2, 3])

                sage: numpy.array(QQ(2**40)).dtype
                dtype('int64')
                sage: numpy.array(QQ(2**400)).dtype
                dtype('O')

                sage: numpy.array([1, 1/2, 3/4])
                array([ 1.  ,  0.5 ,  0.75])
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
            t = ''
            if self < 0:
                t = t + '<mo>-</mo>'
            t = t + '<mfrac><mrow>%s</mrow><mrow>%s</mrow></mfrac>'%(
                mathml(abs(self.numer())), mathml(self.denom()))
            return t

    def _im_gens_(self, codomain, im_gens):
        """
        Return the image of ``self`` under the homomorphism from the rational
        field to ``codomain``.

        This always just returns ``self`` coerced into the ``codomain``.

        INPUT:

        -  ``codomain`` -- object (usually a ring)

        -  ``im_gens`` -- list of elements of ``codomain``

        EXAMPLES::

            sage: a = -17/37
            sage: a._im_gens_(QQ, [1/1])
            -17/37
        """
        return codomain._coerce_(self)

    def lcm(self, Rational other):
        r"""
        Return the least common multiple of self and other.

        One way to define this notion is the following:

        Note that each rational positive rational number can be written as
        a product of primes with integer (positive or negative) powers in a
        unique way.

        Then, the LCM of two rational numbers `x,y` can be defined by
        specifying that the exponent of every prime `p` in ``lcm(x,y)`` is the
        supremum of the exponents of `p` in `x`, and the exponent of `p` in `y`
        (The primes that does not appear in the decomposition of `x` or `y` are
        considered to have exponent zero).

        This definition is consistent with the definition of the LCM in the
        rational integers. Our hopefully interesting notion of LCM for
        rational numbers is illustrated in the examples below.

        EXAMPLES::

            sage: lcm(2/3,1/5)
            2

        This is consistent with the definition above, since:

        .. math::

                        2/3 = 2^1 * 3^{-1}*5^0

        .. math::

                        1/5 = 2^0 * 3^0   *5^{-1}

        and hence,

        .. math::

                        lcm(2/3,1/5)= 2^1*3^0*5^0 = 2.

        ::

            sage: lcm(2/3,7/5)
            14

        In this example:

        .. math::

                        2/3 = 2^1*3^{-1}*5^0    * 7^0

        .. math::

                        7/5 = 2^0*3^0   *5^{-1} * 7^1

        .. math::

                        lcm(2/3,7/5) = 2^1*3^0*5^0*7^1 = 14

        ::

            sage: lcm(1/3,1/5)
            1

        In this example:

        .. math::

                        1/3 = 3^{-1}*5^0

        .. math::

                        1/5 = 3^0 * 5^{-1}

        .. math::

                        lcm(1/3,1/5)=3^0*5^0=1

        ::

            sage: lcm(1/3,1/6)
            1/3

        In this example:

        .. math::

                        1/3 = 2^0*3^{-1}

        .. math::

                        1/6 = 2^{-1}*3^{-1}

        .. math::

                        lcm(1/3,1/6)=2^0*3^{-1}=1/3
        """
        d = self.denom()*other.denom()
        self_d = self.numer()*other.denom()
        other_d = other.numer()*self.denom()
        return self_d.lcm(other_d) / d

    def content(self, other):
        """
        Return the content of ``self`` and ``other``, i.e. the unique positive
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
        from sage.rings.arith import gcd, lcm
        return gcd(nums) / lcm(denoms)

#    def gcd_rational(self, other, **kwds):
#        """
#        Return a gcd of the rational numbers self and other.
#
#        If self = other = 0, this is by convention 0.  In all other
#        cases it can (mathematically) be any nonzero rational number,
#        but for simplicity we choose to always return 1.
#
#        EXAMPLES::
#
#            sage: (1/3).gcd_rational(2/1)
#            1
#            sage: (1/1).gcd_rational(0/1)
#            1
#            sage: (0/1).gcd_rational(0/1)
#            0
#        """
#        if self == 0 and other == 0:
#            return Rational(0)
#        else:
#            return Rational(1)

    def valuation(self, p):
        r"""
        Return the power of ``p`` in the factorization of self.

        INPUT:


        -  ``p`` - a prime number

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
        Returns the local height of this rational number at the prime `p`.

        INPUT:

        -  ``p`` -- a prime number

        - ``prec`` (int) -- desired floating point precision (default:
          default RealField precision).

        OUTPUT:

        (real) The local height of this rational number at the
        prime `p`.

        EXAMPLES::

            sage: a = QQ(25/6)
            sage: a.local_height(2)
            0.693147180559945
            sage: a.local_height(3)
            1.09861228866811
            sage: a.local_height(5)
            0.000000000000000
        """
        from sage.rings.real_mpfr import RealField
        if prec is None:
            R = RealField()
        else:
            R = RealField(prec)
        if self.is_zero():
            return R.zero_element()
        val = self.valuation(p)
        if val >= 0:
            return R.zero_element()
        return -val * R(p).log()

    def local_height_arch(self, prec=None):
        r"""
        Returns the Archimedean local height of this rational number at the
        infinite place.

        INPUT:

        - ``prec`` (int) -- desired floating point precision (default:
          default RealField precision).

        OUTPUT:

        (real) The local height of this rational number `x` at the
        unique infinite place of `\QQ`, which is
        `\max(\log(|x|),0)`.

        EXAMPLES::

            sage: a = QQ(6/25)
            sage: a.local_height_arch()
            0.000000000000000
            sage: (1/a).local_height_arch()
            1.42711635564015
            sage: (1/a).local_height_arch(100)
            1.4271163556401457483890413081
        """
        from sage.rings.real_mpfr import RealField
        if prec is None:
            R = RealField()
        else:
            R = RealField(prec)
        a = self.abs()
        if a <= 1:
            return R.zero_element()
        return R(a).log()

    def global_height_non_arch(self, prec=None):
        r"""
        Returns the total non-archimedean component of the height of this
        rational number.

        INPUT:

        - ``prec`` (int) -- desired floating point precision (default:
          default RealField precision).

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
            sage: a.global_height_non_arch()
            1.79175946922805
            sage: [a.local_height(p) for p in a.support()]
            [0.693147180559945, 1.09861228866811, 0.000000000000000]
            sage: sum([a.local_height(p) for p in a.support()])
            1.79175946922805
        """
        from sage.rings.real_mpfr import RealField
        if prec is None:
            R = RealField()
        else:
            R = RealField(prec)
        d = self.denominator()
        if d.is_one():
            return R.zero_element()
        return R(d).log()

    def global_height_arch(self, prec=None):
        r"""
        Returns the total archimedean component of the height of this rational
        number.

        INPUT:

        - ``prec`` (int) -- desired floating point precision (default:
          default RealField precision).

        OUTPUT:

        (real) The total archimedean component of the height of
        this rational number.

        ALGORITHM:

        Since `\QQ` has only one infinite place this is just the value
        of the local height at that place.  This separate function is
        included for compatibility with number fields.

        EXAMPLES::

            sage: a = QQ(6/25)
            sage: a.global_height_arch()
            0.000000000000000
            sage: (1/a).global_height_arch()
            1.42711635564015
            sage: (1/a).global_height_arch(100)
            1.4271163556401457483890413081
        """
        return self.local_height_arch(prec)

    def global_height(self, prec=None):
        r"""
        Returns the absolute logarithmic height of this rational number.

        INPUT:

        - ``prec`` (int) -- desired floating point precision (default:
          default RealField precision).

        OUTPUT:

        (real) The absolute logarithmic height of this rational number.

        ALGORITHM:

        The height is the sum of the total archimedean and
        non-archimedean components, which is equal to
        `\max(\log(n),\log(d))` where `n,d` are the numerator and
        denominator of the rational number.

        EXAMPLES::

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

        OUTPUT: bool

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
         - ``element`` -- (default: ``False``) boolean whether to also output
           an element of which ``self`` is a norm
         - proof -- If ``True``, then the output is correct unconditionally.
           If ``False``, then the output assumes GRH.

        OUTPUT:

        If element is ``False``, then the output is a boolean ``B``, which is
        ``True`` if and only if ``self`` is the norm of an element of ``L``.
        If ``element`` is ``False``, then the output is a pair ``(B, x)``,
        where ``B`` is as above. If ``B`` is ``True``, then ``x`` an element of
        ``L`` such that ``self == x.norm()``. Otherwise, ``x is None``.

        ALGORITHM:

        Uses PARI's bnfisnorm. See ``_bnfisnorm()``.

        EXAMPLES::

            sage: K = NumberField(x^2 - 2, 'beta')
            sage: (1/7).is_norm(K)
            True
            sage: (1/10).is_norm(K)
            False
            sage: 0.is_norm(K)
            True
            sage: (1/7).is_norm(K, element=True)
            (True, -1/7*beta - 3/7)
            sage: (1/10).is_norm(K, element=True)
            (False, None)
            sage: (1/691).is_norm(QQ, element=True)
            (True, 1/691)

        The number field doesn't have to be defined by an
        integral polynomial::

            sage: B, e = (1/5).is_norm(QuadraticField(5/4, 'a'), element=True)
            sage: B
            True
            sage: e.norm()
            1/5

        A non-Galois number field::

            sage: K.<a> = NumberField(x^3-2)
            sage: B, e = (3/5).is_norm(K, element=True); B
            True
            sage: e.norm()
            3/5

            sage: 7.is_norm(K)
            Traceback (most recent call last):
            ...
            NotImplementedError: is_norm is not implemented unconditionally for norms from non-Galois number fields
            sage: 7.is_norm(K, proof=False)
            False

        AUTHORS:

        - Craig Citro (2008-04-05)

        - Marco Streng (2010-12-03)
        """
        if not element:
            return self.is_norm(L, element=True, proof=proof)[0]

        from sage.rings.number_field.all import is_NumberField
        if not is_NumberField(L):
            raise ValueError, "L (=%s) must be a NumberField in is_norm" % L
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
            raise NotImplementedError, "is_norm is not implemented unconditionally for norms from non-Galois number fields"
        return False, None

    def _bnfisnorm(self, K, proof=True, extra_primes=0):
        r"""
        This gives the output of the PARI function bnfisnorm.

        Tries to tell whether the rational number ``self`` is the norm of some
        element `y` in ``K``. Returns a pair `(a, b)` where
        ``self = Norm(a)*b``. Looks for a solution that is an `S`-unit, with
        `S` a certain set of prime ideals containing (among others) all primes
        dividing ``self``.

        If `K` is known to be Galois, set ``extra_primes = 0`` (in this case,
        ``self`` is a norm iff `b = 1`).

        If ``extra_primes`` is non-zero, the program adds to `S` the following
        prime ideals, depending on the sign of extra_primes.
        If ``extra_primes > 0``, the ideals of norm less than ``extra_primes``.
        And if ``extra_primes < 0``, the ideals dividing ``extra_primes``.

        Assuming GRH, the answer is guaranteed (i.e., ``self`` is a norm
        iff `b = 1`), if `S` contains all primes less than
        `12\log(\disc(L))^2`,
        where `L` is the Galois closure of `K`.

        INPUT:

         - ``K`` -- a number field
         - ``proof`` -- whether to certify the output of bnfinit.
           If ``False``, then correctness of the output depends on GRH.
         - ``extra_primes`` -- an integer as explained above

        OUTPUT:

        A pair `(a, b)` with `a` in `K` and `b` in `\QQ` such that
        ``self == Norm(a)*b`` as explained above.

        ALGORITHM:

        Uses PARI's bnfisnorm.

        EXAMPLES::

            sage: QQ(2)._bnfisnorm(QuadraticField(-1, 'i'))
            (i + 1, 1)
            sage: 7._bnfisnorm(NumberField(x^3-2, 'b'))
            (1, 7)

        AUTHORS:

        - Craig Citro (2008-04-05)

        - Marco Streng (2010-12-03)
        """
        from sage.rings.number_field.all import is_NumberField
        if not is_NumberField(K):
            raise ValueError, "K must be a NumberField in bnfisnorm"

        a, b = K.pari_bnf(proof=proof).bnfisnorm(self, flag=extra_primes)
        return K(a), Rational(b)


    def is_perfect_power(self, expected_value=False):
        r"""
        Returns ``True`` if ``self`` is a perfect power.

        INPUT:

        - ``expected_value`` - (bool) whether or not this rational is expected
          be a perfect power. This does not affect the  correctness of the
          output, only the runtime.

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

        This test makes sure we workaround a bug in GMP (see :trac:`4612`)::

            sage: [ -a for a in srange(100) if not QQ(-a^3).is_perfect_power() ]
            []
            sage: [ -a for a in srange(100) if not QQ(-a^3).is_perfect_power(True) ]
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
        # avoid a bug in GMP's mpz_perfect_power_p; see trac #4612 for
        # more details.
        #
        # The code in the case of sign == -1 could definitely be
        # cleaned up, but it will be removed shortly, since both GMP
        # and eMPIRe have fixes for the mpz_perfect_power_p bug.

        s = mpz_sgn(mpq_numref(self.value))
        if s == 1: # self is positive

            if (mpz_cmp_ui(mpq_denref(self.value), 1) == 0):
                return mpz_perfect_power_p(mpq_numref(self.value))
            if expected_value == False:
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

            if expected_value == False:
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
        Return the square free part of `x`, i.e., an integer z such
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

    def sqrt_approx(self, prec=None, all=False):
        """
        Return numerical approximation with given number of bits of
        precision to this rational number. If all is given, return both
        approximations.

        INPUT:

        -  ``prec`` -- integer

        -  ``all`` -- bool

        EXAMPLES::

            sage: (5/3).sqrt_approx()
            doctest:...: DeprecationWarning: This function is deprecated.  Use sqrt with a given number of bits of precision instead.
            See http://trac.sagemath.org/9859 for details.
            1.29099444873581
            sage: (990829038092384908234098239048230984/4).sqrt_approx()
            4.9770197862083713747374920870362581922510725585130996993055116540856385e17
            sage: (5/3).sqrt_approx(prec=200)
            1.2909944487358056283930884665941332036109739017638636088625
            sage: (9/4).sqrt_approx()
            3/2
        """
        from sage.misc.superseded import deprecation
        deprecation(9859, "This function is deprecated.  Use sqrt with a given number of bits of precision instead.")
        try:
            return self.sqrt(extend=False,all=all)
        except ValueError:
            pass
        if prec is None:
            prec = max(max(53, 2*(mpz_sizeinbase(mpq_numref(self.value), 2)+2)),
                   2*(mpz_sizeinbase(mpq_denref(self.value), 2)+2))
        return self.sqrt(prec=prec, all=all)

    def is_padic_square(self, p):
        """
        Determines whether this rational number is a square in `\QQ_p` (or in
        `R` when ``p = infinity``).

        INPUT:

        -  ``p`` - a prime number, or ``infinity``

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
        p = ZZ(p)
        if not p.is_prime():
            raise ValueError, 'p must be "infinity" or a positive prime number.'

        ## Deal with finite primes
        e, m = self.val_unit(p)

        if e % 2 == 1:
            return False

        if p == 2:
            return ((m % 8) == 1)

        from sage.rings.arith import kronecker_symbol
        return (kronecker_symbol(m, p) == 1)

    def val_unit(self, p):
        r"""
        Returns a pair: the `p`-adic valuation of ``self``, and the `p`-adic
        unit of ``self``, as a :class:`Rational`.

        We do not require the `p` be prime, but it must be at least 2. For
        more documentation see :meth:`Integer.val_unit()`.

        INPUT:

        -  ``p`` - a prime

        OUTPUT:

        -  ``int`` - the `p`-adic valuation of this rational

        -  ``Rational`` - `p`-adic unit part of ``self``

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
        cdef integer.Integer v
        cdef Rational u
        if mpz_cmp_ui(p.value, 2) < 0:
            raise ValueError, "p must be at least 2."
        if mpq_sgn(self.value) == 0:
            import sage.rings.infinity
            u = PY_NEW(Rational)
            mpq_set_ui(u.value, 1, 1)
            return (sage.rings.infinity.infinity, u)
        v = PY_NEW(integer.Integer)
        u = PY_NEW(Rational)
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
        Returns ``self`` with all powers of all primes in ``S`` removed.

        INPUT:

        -  ``S`` - list or tuple of primes.

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

        -  ``prec`` -- integer (default: ``None``): if ``None``, returns
           an exact square root; otherwise returns a numerical square root if
           necessary, to the given bits of precision.

        -  ``extend`` -- bool (default: ``True``); if ``True``, return a
           square root in an extension ring, if necessary. Otherwise, raise a
           ``ValueError`` if the square is not in the base ring.

        -  ``all`` -- bool (default: ``False``); if ``True``, return all
           square roots of self, instead of just one.

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
            sage: x.sqrt()
            9*sqrt(1/5)
            sage: x = -81/3
            sage: x.sqrt()
            3*sqrt(-3)

        ::

            sage: n = 2/3
            sage: n.sqrt()
            sqrt(2/3)
            sage: n.sqrt(prec=10)
            0.82
            sage: n.sqrt(prec=100)
            0.81649658092772603273242802490
            sage: n.sqrt(prec=100)^2
            0.66666666666666666666666666667
            sage: n.sqrt(prec=53, all=True)
            [0.816496580927726, -0.816496580927726]
            sage: n.sqrt(extend=False, all=True)
            Traceback (most recent call last):
            ...
            ValueError: square root of 2/3 not a rational number
            sage: sqrt(-2/3, all=True)
            [sqrt(-2/3), -sqrt(-2/3)]
            sage: sqrt(-2/3, prec=53)
            0.816496580927726*I
            sage: sqrt(-2/3, prec=53, all=True)
            [0.816496580927726*I, -0.816496580927726*I]

        AUTHORS:

        - Naqi Jaffery (2006-03-05): some examples
        """
        if mpq_sgn(self.value) == 0:
            return [self] if all else self

        if mpq_sgn(self.value) < 0:
            if not extend:
                raise ValueError, "square root of negative number not rational"
            from sage.functions.other import _do_sqrt
            return _do_sqrt(self, prec=prec, all=all)

        cdef Rational z = <Rational> PY_NEW(Rational)
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
            if not extend:
                raise ValueError, "square root of %s not a rational number"%self
            from sage.functions.other import _do_sqrt
            return _do_sqrt(self, prec=prec, all=all)

        if prec:
            from sage.functions.other import _do_sqrt
            return _do_sqrt(self, prec=prec, all=all)

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
        `10^s=10^{s+t} \mod d`. In general if `d=2^a3^bm` where `m`
        is coprime to 10, then `s=\max(a,b)` and `t` is the order of
        10 modulo `d`.

        EXAMPLES::

            sage: (1/7).period()
            6
            sage: RR(1/7)
            0.142857142857143
            sage: (1/8).period()
            1
            sage: RR(1/8)
            0.125000000000000
            sage: RR(1/6)
            0.166666666666667
            sage: (1/6).period()
            1
            sage: x = 333/106
            sage: x.period()
            13
            sage: RealField(200)(x)
            3.1415094339622641509433962264150943396226415094339622641509
        """
        cdef unsigned int alpha, beta
        d = self.denominator()
        alpha, d = d.val_unit(2)
        beta, d  = d.val_unit(5)
        from sage.rings.finite_rings.integer_mod import Mod
        return Mod(ZZ(10),d).multiplicative_order()

    def nth_root(self, int n):
        r"""
        Computes the `n`-th root of ``self``, or raises a
        ``ValueError`` if ``self`` is not a perfect `n`-th power.

        INPUT:

        -  ``n`` - integer (must fit in C int type)

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
            raise ValueError, "n cannot be zero"

        num, exact = self.numerator().nth_root(n, 1)
        if not exact:
            raise ValueError, "not a perfect %s power"%ZZ(n).ordinal_str()

        den, exact = self.denominator().nth_root(n, 1)
        if not exact:
            raise ValueError, "not a perfect %s power"%ZZ(n).ordinal_str()

        if negative:
            return den / num
        else:
            return num / den

    def is_nth_power(self, int n):
        r"""
        Returns ``True`` if self is an `n`-th power, else ``False``.

        INPUT:

        -  ``n`` - integer (must fit in C int type)

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
            raise ValueError, "n cannot be zero"
        if n<0:
            n = -n
        if n%2==0 and self<0:
            return False
        return self.numerator().nth_root(n, 1)[1]\
               and self.denominator().nth_root(n, 1)[1]

    def str(self, int base=10):
        """
        Return a string representation of ``self`` in the given ``base``.

        INPUT:

        -  ``base`` -- integer (default: 10); base must be between 2 and 36.

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
            raise ValueError, "base (=%s) must be between 2 and 36"%base
        cdef size_t n
        cdef char *s

        n = mpz_sizeinbase (mpq_numref(self.value), base) \
            + mpz_sizeinbase (mpq_denref(self.value), base) + 3
        s = <char *>PyMem_Malloc(n)
        if s == NULL:
            raise MemoryError, "Unable to allocate enough memory for the string representation of an integer."

        sig_on()
        mpq_get_str(s, base, self.value)
        sig_off()
        k = <object> PyString_FromString(s)
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

        TESTS:

        Test that conversion agrees with `RR`::

            sage: Q = [a/b for a in [-99..99] for b in [1..99]]
            sage: all([RDF(q) == RR(q) for q  in Q])
            True

        Test that the conversion has correct rounding on simple rationals::

            sage: for p in [-100..100]:
            ....:   for q in [1..100]:
            ....:       r = RDF(p/q)
            ....:       assert (RR(r).exact_rational() - p/q) <= r.ulp()/2

        Test larger rationals::

            sage: Q = continued_fraction(pi, bits=3000).convergents()
            sage: all([RDF(q) == RR(q) for q  in Q])
            True

        At some point, the continued fraction and direct conversion
        to ``RDF`` should agree::

            sage: RDFpi = RDF(pi)
            sage: all([RDF(q) == RDFpi for q  in Q[20:]])
            True
        """
        return mpq_get_d_nearest(self.value)

    def __hash__(self):
        """
        Return hash of ``self``.

        OUTPUT: integer

        EXAMPLES::

            sage: (-4/17).__hash__()
            -19
            sage: (-5/1).__hash__()
            -5
        """
        cdef long n, d
        n = mpz_pythonhash(mpq_numref(self.value))
        d = mpz_pythonhash(mpq_denref(self.value))
        if d == 1:
            return n
        n = n ^ d
        if n == -1:
            return -2
        return n

    def __getitem__(self, int n):
        """
        Return ``n``-th element of ``self``, viewed as a list. This is for
        consistency with how number field elements work.

        INPUT:

        -  ``n`` - an integer (error if not 0 or -1)

        OUTPUT: Rational

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
        raise IndexError, "index n (=%s) out of range; it must be 0"%n

    ################################################################
    # Optimized arithmetic
    ################################################################
    cpdef ModuleElement _add_(self, ModuleElement right):
        """
        Return ``right`` plus ``self``.

        EXAMPLES::

            sage: (2/3) + (1/6) # indirect doctest
            5/6
            sage: (1/3) + (1/2) # indirect doctest
            5/6
        """
        cdef Rational x
        x = <Rational> PY_NEW(Rational)
        mpq_add(x.value, self.value, (<Rational>right).value)
        return x

    cpdef ModuleElement _iadd_(self, ModuleElement right):
        """
        Add ``right`` to ``self`` (and store in ``self``).

            sage: p = 2/3
            sage: p += 1/6; p # indirect doctest
            5/6
        """
        mpq_add(self.value, self.value, (<Rational>right).value)
        return self

    cpdef ModuleElement _sub_(self, ModuleElement right):
        """
        Return ``self`` minus ``right``.

        EXAMPLES::

            sage: (2/3) - (1/6) # indirect doctest
            1/2
        """
        # self and right are guaranteed to be Integers
        cdef Rational x
        x = <Rational> PY_NEW(Rational)
        mpq_sub(x.value, self.value, (<Rational>right).value)
        return x

    cpdef ModuleElement _isub_(self, ModuleElement right):
        """
        Subtract ``right`` from ``self`` (and store in ``self``).

        EXAMPLES::

            sage: p = (2/3)
            sage: p -= 1/6; p # indirect doctest
            1/2
        """
        mpq_sub(self.value, self.value, (<Rational>right).value)
        return self

    cpdef ModuleElement _neg_(self):
        """
        Negate ``self``.

        EXAMPLES::

            sage: -(2/3) # indirect doctest
            -2/3
        """
        cdef Rational x
        x = <Rational> PY_NEW(Rational)
        mpq_neg(x.value, self.value)
        return x

    cpdef RingElement _mul_(self, RingElement right):
        """
        Return ``self`` times ``right``.

        EXAMPLES::

            sage: (3/14) * 2 # indirect doctest
            3/7
        """
        cdef Rational x
        x = <Rational> PY_NEW(Rational)
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

    cpdef RingElement _imul_(self, RingElement right):
        """
        Multiply ``right`` with ``self`` (and store in ``self``).

        EXAMPLES::
            sage: p = 3/14
            sage: p *= 2; p # indirect doctest
            3/7
        """
        if mpz_sizeinbase (mpq_numref(self.value), 2)  > 100000 or \
             mpz_sizeinbase (mpq_denref(self.value), 2) > 100000:
            # We only use the signal handler (to enable ctrl-c out) in case
            # self is huge, so the product might actually take a while to compute.
            sig_on()
            mpq_mul(self.value, self.value, (<Rational>right).value)
            sig_off()
        else:
            mpq_mul(self.value, self.value, (<Rational>right).value)
        return self

    cpdef RingElement _div_(self, RingElement right):
        """
        Return ``self`` divided by ``right``.

        EXAMPLES::

            sage: 2/3 # indirect doctest
            2/3
            sage: 3/0 # indirect doctest
            Traceback (most recent call last):
            ...
            ZeroDivisionError: Rational division by zero
        """
        if mpq_cmp_si((<Rational> right).value, 0, 1) == 0:
            raise ZeroDivisionError, "Rational division by zero"
        cdef Rational x
        x = <Rational> PY_NEW(Rational)
        mpq_div(x.value, self.value, (<Rational>right).value)
        return x

    cpdef RingElement _idiv_(self, RingElement right):
        """
        Divide ``self`` by ``right`` (and store in ``self``).

        EXAMPLES::

            sage: p = 2/3
            sage: p /= 2; p # indirect doctest
            1/3
        """
        if mpq_cmp_si((<Rational> right).value, 0, 1) == 0:
            raise ZeroDivisionError, "Rational division by zero"
        mpq_div(self.value, self.value, (<Rational>right).value)
        return self


    ################################################################
    # Other arithmetic operations.
    ################################################################

    def __invert__(self):
        """
        Return the multiplicative inverse of ``self``.

        OUTPUT: Rational

        EXAMPLES::

            sage: (-4/17).__invert__()
            -17/4
            sage: ~(-4/17)
            -17/4
        """
        if self.is_zero():
            raise ZeroDivisionError, "rational division by zero"
        cdef Rational x
        x = <Rational> PY_NEW(Rational)
        mpq_inv(x.value, self.value)
        return x

    def __pow__(self, n, dummy):
        """
        Raise ``self`` to the integer power ``n``.

        EXAMPLES::

            sage: (2/3)^5
            32/243
            sage: (-1/1)^(1/3)
            (-1)^(1/3)

        We raise to some interesting powers::

            sage: (2/3)^I
            (2/3)^I
            sage: (2/3)^sqrt(2)
            (2/3)^sqrt(2)
            sage: x,y,z,n = var('x,y,z,n')
            sage: (2/3)^(x^n + y^n + z^n)
            (2/3)^(x^n + y^n + z^n)
            sage: (-7/11)^(tan(x)+exp(x))
            (-7/11)^(e^x + tan(x))
            sage: (2/3)^(3/4)
            (2/3)^(3/4)
            sage: (-1/3)^0
            1
            sage: a = (0/1)^(0/1); a
            1
            sage: type(a)
            <type 'sage.rings.rational.Rational'>

        The exponent must fit in a long unless the base is -1, 0, or 1.

        ::

            sage: s = (1/2)^(2^100)
            Traceback (most recent call last):
            ...
            RuntimeError: exponent must be at most 2147483647  # 32-bit
            RuntimeError: exponent must be at most 9223372036854775807 # 64-bit
            sage: s = (1/2)^(-2^100)
            Traceback (most recent call last):
            ...
            RuntimeError: exponent must be at most 2147483647  # 32-bit
            RuntimeError: exponent must be at most 9223372036854775807 # 64-bit
            sage: (-3/3)^(2^100)
            1

        This works even if the base is a float or Python complex or other
        type::

            sage: float(1.2)**(1/2)
            1.0954451150103321
            sage: complex(1,2)**(1/2)
            (1.272019649514069+0.786151377757423...j)
            sage: int(2)^(1/2)
            sqrt(2)
            sage: a = int(2)^(3/1); a
            8
            sage: type(a)
            <type 'sage.rings.rational.Rational'>

        If the result is rational, it is returned as a rational::

            sage: (4/9)^(1/2)
            2/3
            sage: parent((4/9)^(1/2))
            Rational Field
            sage: (-27/125)^(1/3)
            3/5*(-1)^(1/3)
            sage: (-27/125)^(1/2)
            3/5*sqrt(-3/5)

        The result is normalized to have the rational power in the numerator::

            sage: 2^(-1/2)
            1/2*sqrt(2)
            sage: 8^(-1/5)
            1/8*8^(4/5)
            sage: 3^(-3/2)
            1/9*sqrt(3)
        """
        if dummy is not None:
            raise ValueError, "__pow__ dummy variable not used"

        if not PY_TYPE_CHECK(self, Rational):
            # If the base is not a rational, e.g., it is an int, complex, float, user-defined type, etc.
            try:
                self_coerced = n.parent().coerce(self)
            except TypeError:
                n_coerced = type(self)(n)
                if n != n_coerced:
                    # dangerous coercion -- don't use -- try symbolic result
                    from sage.calculus.calculus import SR
                    return SR(self)**SR(n)
                return self.__pow__(n_coerced)
            return self_coerced.__pow__(n)

        cdef Rational _self = <Rational>self
        cdef long nn

        try:
            nn = PyNumber_Index(n)
        except TypeError:
            if PY_TYPE_CHECK(n, Rational):
                # Perhaps it can be done exactly
                c, d = rational_power_parts(self, n)
                if d == 1:
                    # It was an exact power
                    return c
                elif d == -1 and n.denominator() == 2:
                    from sage.symbolic.all import I
                    return c * I.pyobject() ** (n.numerator() % 4)
                elif c != 1:
                    return c * d**n
                # Can't simplify the power, return a symbolic expression.
                # We use the hold=True keyword argument to prevent the
                # symbolics library from trying to simplify this expression
                # again. This would lead to infinite loops otherwise.
                from sage.symbolic.ring import SR
                if n < 0:
                    int_exp = -(n.floor())
                    return SR(1/self**int_exp)*\
                            SR(self).power(int_exp+n, hold=True)
                else:
                    return SR(self).power(n, hold=True)

            if PY_TYPE_CHECK(n, Element):
                return (<Element>n)._parent(self)**n
            try:
                return n.parent()(self)**n
            except AttributeError:
                try:
                    return type(n)(self)**n
                except StandardError:
                    raise TypeError, "exponent (=%s) must be an integer.\nCoerce your numbers to real or complex numbers first."%n

        except OverflowError:
            if mpz_cmp_si(mpq_denref(_self.value), 1) == 0:
                if mpz_cmp_si(mpq_numref(_self.value), 1) == 0:
                    return self
                elif mpz_cmp_si(mpq_numref(_self.value), 0) == 0:
                    return self
                elif mpz_cmp_si(mpq_numref(_self.value), -1) == 0:
                    return self if n % 2 else -self
            raise RuntimeError, "exponent must be at most %s" % sys.maxint

        cdef Rational x = <Rational> PY_NEW(Rational)

        if nn == 0:
            mpq_set_si(x.value, 1, 1)
            return x

        if nn < 0:
            sig_on()
            # mpz_pow_ui(mpq_denref(x.value), mpq_numref(_self.value), <unsigned long int>(-nn))
            # mpz_pow_ui(mpq_numref(x.value), mpq_denref(_self.value), <unsigned long int>(-nn))
            # The above causes segfaults, so swap after instead...
            mpz_pow_ui(mpq_numref(x.value), mpq_numref(_self.value), -nn)
            mpz_pow_ui(mpq_denref(x.value), mpq_denref(_self.value), -nn)
            # mpz_swap(mpq_numref(x.value), mpq_denref(x.value)) # still a segfault
            mpq_inv(x.value, x.value)
            sig_off()
            return x
        elif nn > 0:
            sig_on()
            mpz_pow_ui(mpq_numref(x.value), mpq_numref(_self.value), nn)
            mpz_pow_ui(mpq_denref(x.value), mpq_denref(_self.value), nn)
            sig_off()
            return x


        if n < 0:  # this doesn't make sense unless n is an integer.
            x = _self**(-n)
            return x.__invert__()

        cdef mpz_t num, den

        sig_on()
        mpz_init(num)
        mpz_init(den)
        mpz_pow_ui(num, mpq_numref(_self.value), nn)
        mpz_pow_ui(den, mpq_denref(_self.value), nn)
        mpq_set_num(x.value, num)
        mpq_set_den(x.value, den)
        mpz_clear(num)
        mpz_clear(den)
        sig_off()

        return x

    def __pos__(self):
        """
        Return ``self``.

        OUTPUT: Rational

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

        OUTPUT: Rational

        EXAMPLES::

            sage: (-4/17).__neg__()
            4/17
            sage: - (-4/17)
            4/17
        """
        cdef Rational x
        x = <Rational> PY_NEW(Rational)
        mpq_neg(x.value, self.value)
        return x

    def __nonzero__(self):
        """
        Return ``True`` if this rational number is nonzero.

        OUTPUT: bool

        EXAMPLES::

            sage: (-4/17).__nonzero__()
            True
            sage: (0/5).__nonzero__()
            False
            sage: bool(-4/17)
            True
        """
        # A rational number is zero iff its numerator is zero.
        return mpq_sgn(self.value) != 0

    def __abs__(self):
        """
        Return the absolute value of this rational number.

        OUTPUT: Rational

        EXAMPLES::

            sage: (-4/17).__abs__()
            4/17
            sage: abs(-4/17)
            4/17
        """
        cdef Rational x
        x = <Rational> PY_NEW(Rational)
        mpq_abs(x.value, self.value)
        return x

    def sign(self):
        """
        Returns the sign of this rational number, which is -1, 0, or 1
        depending on whether this number is negative, zero, or positive
        respectively.

        OUTPUT: Integer

        EXAMPLES::

            sage: (2/3).sign()
            1
            sage: (0/3).sign()
            0
            sage: (-1/6).sign()
            -1
        """
        return integer.smallInteger(mpq_sgn(self.value))

    def mod_ui(Rational self, unsigned long int n):
        """
        Return the remainder upon division of ``self`` by the unsigned long
        integer ``n``.

        INPUT:

        -  ``n`` - an unsigned long integer

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

        -  ``other`` - object that coerces to an integer.

        OUTPUT: integer

        EXAMPLES::

            sage: (-4/17).__mod__(3/1)
            1

        TESTS:

        Check that :trac:`14870` is fixed::

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
        Returns the norm from `\QQ` to `\QQ` of `x` (which is just `x`). This
        was added for compatibility with :class:`NumberFields`.

        OUTPUT:

        -  ``Rational`` - reference to ``self``

        EXAMPLES::

            sage: (1/3).norm()
             1/3

        AUTHORS:

        - Craig Citro
        """
        return self

    def relative_norm(self):
        """
        Returns the norm from Q to Q of x (which is just x). This was added for compatibility with NumberFields

        EXAMPLES::

            sage: (6/5).relative_norm()
            6/5

            sage: QQ(7/5).relative_norm()
            7/5
        """
        return self

    def absolute_norm(self):
        """
        Returns the norm from Q to Q of x (which is just x). This was added for compatibility with NumberFields

        EXAMPLES::

            sage: (6/5).absolute_norm()
            6/5

            sage: QQ(7/5).absolute_norm()
            7/5
        """
        return self

    def trace(self):
        r"""
        Returns the trace from `\QQ` to `\QQ` of `x` (which is just `x`). This
        was added for compatibility with :class:`NumberFields`.

        OUTPUT:

        -  ``Rational`` - reference to self

        EXAMPLES::

            sage: (1/3).trace()
             1/3

        AUTHORS:

        - Craig Citro
        """
        return self

    def charpoly(self, var):
        """
        Return the characteristic polynomial of this rational number. This
        will always be just ``var - self``; this is really here so that code
        written for number fields won't crash when applied to rational
        numbers.

        INPUT:

        -  ``var`` - a string

        OUTPUT: Polynomial

        EXAMPLES::

            sage: (1/3).charpoly('x')
             x - 1/3

        AUTHORS:

        - Craig Citro
        """
        QQ = self.parent()
        return QQ[var]([-self,1])

    def minpoly(self, var):
        """
        Return the minimal polynomial of this rational number. This will
        always be just ``x - self``; this is really here so that code written
        for number fields won't crash when applied to rational numbers.

        INPUT:

        -  ``var`` - a string

        OUTPUT: Polynomial

        EXAMPLES::

            sage: (1/3).minpoly('x')
             x - 1/3

        AUTHORS:

        - Craig Citro
        """
        QQ = self.parent()
        return QQ[var]([-self,1])

    def _integer_(self, Z=None):
        """
        Return ``self`` coerced to an integer. Of course this rational number
        must have a denominator of 1.

        OUTPUT: Integer

        EXAMPLES::

            sage: (-4/17)._integer_()
            Traceback (most recent call last):
            ...
            TypeError: no conversion of this rational to integer
            sage: (-4/1)._integer_()
            -4
        """
        if not mpz_cmp_si(mpq_denref(self.value), 1) == 0:
            raise TypeError, "no conversion of this rational to integer"
        cdef integer.Integer n
        n = PY_NEW(integer.Integer)
        n.set_from_mpz(mpq_numref(self.value))
        return n

    # TODO -- this should be deprecated
    def numer(self):
        """
        Return the numerator of this rational number.

        EXAMPLE::

            sage: x = -5/11
            sage: x.numer()
            -5
        """
        cdef integer.Integer n
        n = PY_NEW(integer.Integer)
        n.set_from_mpz(mpq_numref(self.value))
        return n

    def numerator(self):
        """
        Return the numerator of this rational number.

        EXAMPLE::

            sage: x = 5/11
            sage: x.numerator()
            5

        ::

            sage: x = 9/3
            sage: x.numerator()
            3
        """
        cdef integer.Integer n
        n = PY_NEW(integer.Integer)
        n.set_from_mpz(mpq_numref(self.value))
        return n

    def __int__(self):
        """
        Return coercion of ``self`` to Python ``int``.

        This takes the floor of ``self`` if ``self`` has a denominator (which
        is consistent with Python's ``long(floats)``).

        EXAMPLES::

            sage: int(7/3)
            2
            sage: int(-7/3)
            -3
        """
        return int(self.__long__())

    def __long__(self):
        """
        Return coercion of ``self`` to Python ``long``.

        This takes the floor of ``self`` if ``self`` has a denominator (which
        is consistent with Python's ``long(floats)``).

        EXAMPLES::

            sage: long(7/3)
            2L
            sage: long(-7/3)
            -3L
        """
        cdef mpz_t x
        if mpz_cmp_si(mpq_denref(self.value),1) != 0:
            mpz_init(x)
            mpz_fdiv_q(x, mpq_numref(self.value), mpq_denref(self.value))
            n = mpz_get_pylong(x)
            mpz_clear(x)
            return n
        else:
            return mpz_get_pylong(mpq_numref(self.value))

    def denom(self):
        """
        Returns the denominator of this rational number.

        EXAMPLES::

            sage: x = 5/13
            sage: x.denom()
            13
            sage: x = -9/3
            sage: x.denom()
            1
        """
        cdef integer.Integer n
        n = PY_NEW(integer.Integer)
        n.set_from_mpz(mpq_denref(self.value))
        return n

    def denominator(self):
        """
        Returns the denominator of this rational number.

        EXAMPLES::

            sage: x = -5/11
            sage: x.denominator()
            11
            sage: x = 9/3
            sage: x.denominator()
            1
        """
        cdef integer.Integer n
        n = PY_NEW(integer.Integer)
        n.set_from_mpz(mpq_denref(self.value))
        return n

    def factor(self):
        """
        Return the factorization of this rational number.

        OUTPUT: Factorization

        EXAMPLES::

            sage: (-4/17).factor()
            -1 * 2^2 * 17^-1

        Trying to factor 0 gives an arithmetic error::

            sage: (0/1).factor()
            Traceback (most recent call last):
            ...
            ArithmeticError: Prime factorization of 0 not defined.
        """
        return self.numerator().factor() * \
           sage.structure.factorization.Factorization([(p,-e) for p, e in self.denominator().factor()])

    def support(self):
        """
        Return a sorted list of the primes where this rational number has
        non-zero valuation.

        OUTPUT: The set of primes appearing in the factorization of this
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
            raise ArithmeticError, "Support of 0 not defined."
        return sage.rings.arith.prime_factors(self)

    def gamma(self, prec=None):
        """
        Return the gamma function evaluated at ``self``. This value is exact
        for integers and half-integers, and returns a symbolic value
        otherwise.  For a numerical approximation, use keyword ``prec``.

        EXAMPLES::

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

            sage: (1/3).gamma(prec=100)
            2.6789385347077476336556929410
            sage: (1/2).gamma(prec=100)
            1.7724538509055160272981674833
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
                from sage.functions.all import sqrt
                return sqrt(pi) * rat_part
            else:
                from sage.symbolic.all import SR
                return SR(self).gamma()

    def floor(self):
        """
        Return the floor of this rational number as an integer.

        OUTPUT: Integer

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

        OUTPUT: Integer

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

    def round(Rational self, mode="away"):
        """
        Returns the nearest integer to ``self``, rounding away from 0 by
        default, for consistency with the builtin Python round.

        INPUT:

        -  ``self`` - a rational number

        -  ``mode`` - a rounding mode for half integers:

           - 'toward' rounds toward zero
           - 'away' (default) rounds away from zero
           - 'up' rounds up
           - 'down' rounds down
           - 'even' rounds toward the even integer
           - 'odd' rounds toward the odd integer


        OUTPUT: Integer

        EXAMPLES::

            sage: (9/2).round()
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

    def real(self):
        """
        Returns the real part of ``self``, which is ``self``.

        EXAMPLES::

            sage: (1/2).real()
            1/2
        """
        return self

    def imag(self):
        """
        Returns the imaginary part of ``self``, which is zero.

        EXAMPLES::

            sage: (1/239).imag()
            0
        """
        return self._parent(0)

    def height(self):
        """
        The max absolute value of the numerator and denominator of ``self``, as
        an :class:`Integer`.

        OUTPUT: Integer

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
        Returns the least common multiple, in the rational numbers, of ``self``
        and ``other``. This function returns either 0 or 1 (as a rational
        number).

        INPUT:

        -  ``other`` - Rational

        OUTPUT:

        -  ``Rational`` - 0 or 1

        EXAMPLES::

            sage: (2/3)._lcm(3/5)
            1
            sage: (0/1)._lcm(0/1)
            0
            sage: type((2/3)._lcm(3/5))
            <type 'sage.rings.rational.Rational'>
        """
        if mpz_cmp_si(mpq_numref(self.value), 0) == 0 and \
               mpz_cmp_si(mpq_numref(other.value), 0) == 0:
            return Rational(0)
        return Rational(1)

    def _gcd(self, Rational other):
        """
        Returns the least common multiple, in the rational numbers, of ``self``
        and ``other``. This function returns either 0 or 1 (as a rational
        number).

        INPUT:

        -  ``other`` - Rational

        OUTPUT:

        -  ``Rational`` - 0 or 1

        EXAMPLES::

            sage: (2/3)._gcd(3/5)
            1
            sage: (0/1)._gcd(0/1)
            0
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

        OUTPUT: Integer or ``infinity``

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
        elif mpz_cmpabs(mpq_numref(self.value),mpq_denref(self.value))==0:
            # if the numerator and the denominator are equal in absolute value,
            # then the rational number is -1
            return integer.Integer(2)
        else:
            return sage.rings.infinity.infinity

    def is_one(self):
        r"""
        Determine if a rational number is one.

        OUTPUT: bool

        EXAMPLES::

            sage: QQ(1/2).is_one()
            False
            sage: QQ(4/4).is_one()
            True
        """
        # A rational number is equal to 1 iff its numerator and denominator are equal
        return mpz_cmp(mpq_numref(self.value),mpq_denref(self.value))==0

    def is_integral(self):
        r"""
        Determine if a rational number is integral (i.e is in
        `\ZZ`).

        OUTPUT: bool

        EXAMPLES::

            sage: QQ(1/2).is_integral()
            False
            sage: QQ(4/4).is_integral()
            True
        """
        return mpz_cmp_si(mpq_denref(self.value), 1) == 0

    def is_S_integral(self, S=[]):
        r"""
        Determine if the rational number is ``S``-integral.

        ``x`` is ``S``-integral if ``x.valuation(p)>=0`` for all ``p`` not in
        ``S``, i.e., the denominator of ``x`` is divisible only by the primes
        in ``S``.

        INPUT:

        -  ``S`` -- list or tuple of primes.

        OUTPUT: bool

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

        -  ``S`` -- list or tuple of primes.

        OUTPUT: bool

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
        if a==1:
            return True
        if S is None:
            return False
        return a.prime_to_S_part(S) == 1

    cdef _lshift(self, long int exp):
        r"""
        Return ``self * 2^exp``.
        """
        cdef Rational x
        x = <Rational> PY_NEW(Rational)
        sig_on()
        if exp < 0:
            mpq_div_2exp(x.value,self.value,-exp)
        else:
            mpq_mul_2exp(x.value,self.value,exp)
        sig_off()
        return x

    def __lshift__(x,y):
        """
        Left shift operator ``x << y``.

        INPUT:

        -  ``x, y`` -- integer or rational

        OUTPUT: Rational

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
        if PY_TYPE_CHECK(x, Rational):
            if isinstance(y, (int, long, integer.Integer)):
                return (<Rational>x)._lshift(y)
            if PY_TYPE_CHECK(y, Rational):
                if mpz_cmp_si(mpq_denref((<Rational>y).value), 1) != 0:
                    raise ValueError, "denominator must be 1"
                return (<Rational>x)._lshift(y)
        return bin_op(x, y, operator.lshift)

    cdef _rshift(self, long int exp):
        r"""
        Return ``self / 2^exp``.
        """
        cdef Rational x
        x = <Rational> PY_NEW(Rational)
        sig_on()
        if exp < 0:
            mpq_mul_2exp(x.value,self.value,-exp)
        else:
            mpq_div_2exp(x.value,self.value,exp)
        sig_off()
        return x

    def __rshift__(x,y):
        """
        Right shift operator ``x >> y``.

        INPUT:

        -  ``x, y`` -- integer or rational

        OUTPUT: Rational

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
        if PY_TYPE_CHECK(x, Rational):
            if isinstance(y, (int, long, integer.Integer)):
                return (<Rational>x)._rshift(y)
            if PY_TYPE_CHECK(y, Rational):
                if mpz_cmp_si(mpq_denref((<Rational>y).value), 1) != 0:
                    raise ValueError, "denominator must be 1"
                return (<Rational>x)._rshift(y)
        return bin_op(x, y, operator.rshift)

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

    def _pari_(self):
        """
        Returns the PARI version of this rational number.

        EXAMPLES::

            sage: n = 9390823/17
            sage: m = n._pari_(); m
            9390823/17
            sage: type(m)
            <type 'sage.libs.pari.gen.gen'>
            sage: m.type()
            't_FRAC'
        """
        cdef PariInstance P = sage.libs.pari.pari_instance.pari
        return P.new_gen_from_mpq_t(self.value)

    def _interface_init_(self, I=None):
        """
        Return representation of this rational suitable for coercing into
        almost any computer algebra system.

        OUTPUT: string

        EXAMPLES::

            sage: (2/3)._interface_init_()
            '2/3'
            sage: kash(3/1).Type()              # optional
            elt-fld^rat
            sage: magma(3/1).Type()             # optional
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

        sage: q= QQ(); float(q)
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

    - Paul Zimmermann, Jeroen Demeyer (:trac:`14416`)
    """
    cdef __mpz_struct* a = <__mpz_struct*>(mpq_numref(x))
    cdef __mpz_struct* b = <__mpz_struct*>(mpq_denref(x))
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
            return -1.0/0.0
        else:
            return 1.0/0.0

    sig_on()

    # Compute q = trunc(a / 2^shift) and let remainder_is_zero be True
    # if and only if no truncation occurred.
    cdef mpz_t q, r
    mpz_init(q)
    mpz_init(r)
    cdef bint remainder_is_zero
    if shift > 0:
        mpz_tdiv_r_2exp(r, a, shift)
        remainder_is_zero = (mpz_cmp_ui(r, 0) == 0)
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

    # Convert q to a 64-bit integer.
    cdef mp_limb_t* q_limbs = (<__mpz_struct*>q)._mp_d
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

    if ((q64 & 1) == 0):
        # Round towards zero
        pass
    else:
        if not remainder_is_zero:
            # Remainder is non-zero: round away from zero
            q64 += 1
        else:
            q64 += (q64 & 2) - 1

    # The conversion of q64 to double is *exact*.
    # This is because q64 is even and satisfies q64 <= 2^54,
    # (with 2^53 <= q64 <= 2^54 unless in the denormal case).
    cdef double d = <double>q64
    if resultsign < 0:
        d = -d
    return ldexp(d, shift)


def pyrex_rational_reconstruction(integer.Integer a, integer.Integer m):
    """
    Find the rational reconstruction of ``a mod m``, if it exists.

    INPUT:

    -  ``a`` - Integer

    -  ``m`` - Integer

    OUTPUT:

    -  ``x`` - rings.rational.Rational

    EXAMPLES::

        sage: Integers(100)(2/3)
        34
        sage: sage.rings.rational.pyrex_rational_reconstruction(34, 100)
        2/3

    TEST:

    Check that ticket #9345 is fixed::

        sage: sage.rings.rational.pyrex_rational_reconstruction(0,0)
        Traceback (most recent call last):
        ...
        ZeroDivisionError: The modulus cannot be zero
        sage: sage.rings.rational.pyrex_rational_reconstruction(ZZ.random_element(-10^6, 10^6), 0)
        Traceback (most recent call last):
        ...
        ZeroDivisionError: The modulus cannot be zero
    """
    if not m.__nonzero__():
        raise ZeroDivisionError("The modulus cannot be zero")
    cdef Rational x
    x = <Rational> PY_NEW(Rational)
    mpq_rational_reconstruction(x.value, a.value, m.value)
    return x

def make_rational(s):
    """
    Make a rational number from ``s`` (a string in base 32)

    INPUT:

    -  ``s`` - string in base 32

    OUTPUT: Rational

    EXAMPLES::

        sage: (-7/15).str(32)
        '-7/f'
        sage: sage.rings.rational.make_rational('-7/f')
        -7/15
    """
    r = Rational()
    r._reduce_set(s)
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
        import integer_ring
        import rational_field
        import sage.categories.homset
        Morphism.__init__(self, sage.categories.homset.Hom(integer_ring.ZZ, rational_field.QQ))

    cpdef Element _call_(self, x):
        """
        Return the image of the morphism on ``x``.

        EXAMPLES::

            sage: sage.rings.rational.Z_to_Q()(2) # indirect doctest
            2
        """
        cdef Rational rat
        rat = <Rational> PY_NEW(Rational)
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

            sage: QQ.coerce_map_from(ZZ).section()
            Generic map:
              From: Rational Field
              To:   Integer Ring
        """
        return Q_to_Z(self._codomain, self._domain)

cdef class Q_to_Z(Map):
    r"""
    A morphism from `\QQ` to `\ZZ`.

    TESTS::

        sage: type(ZZ.convert_map_from(QQ))
        <type 'sage.rings.rational.Q_to_Z'>
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
            raise TypeError, "no conversion of this rational to integer"
        cdef integer.Integer n
        n = <integer.Integer>PY_NEW(integer.Integer)
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
    A morphism from ``int`` to `\QQ`.
    """
    def __init__(self):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: sage.rings.rational.int_to_Q()
            Native morphism:
              From: Set of Python objects of type 'int'
              To:   Rational Field
        """
        import rational_field
        import sage.categories.homset
        from sage.structure.parent import Set_PythonType
        Morphism.__init__(self, sage.categories.homset.Hom(Set_PythonType(int), rational_field.QQ))

    cpdef Element _call_(self, a):
        """
        Return the image of the morphism on ``a``.

        EXAMPLES::

            sage: f = sage.rings.rational.int_to_Q()
            sage: f(int(4)) # indirect doctest
            4
            sage: f(int(4^100))   # random garbage, not an int
            14
        """
        cdef Rational rat
        rat = <Rational> PY_NEW(Rational)
        mpq_set_si(rat.value, PyInt_AS_LONG(a), 1)
        return rat

    def _repr_type(self):
        """
        Return string that describes the type of morphism.

        EXAMPLES::

            sage: sage.rings.rational.int_to_Q()._repr_type()
            'Native'
        """
        return "Native"

