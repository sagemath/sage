"""
Arbitrary precision floating point complex numbers using GNU MPC

This is a binding for the arbitrary-precision floating point library
:ref:`GNU MPC <spkg_mpc>`.

We define a class :class:`MPComplexField`, each instance of which
specifies a field of floating-point complex numbers with
a specified precision shared by the real and imaginary part and a rounding
mode stating the rounding mode directions specific to real and imaginary
parts.

Individual floating-point numbers are of class :class:`MPComplexNumber`.

For floating-point representation and rounding mode description see the
documentation for the :mod:`sage.rings.real_mpfr`.

AUTHORS:

- Philippe Theveny (2008-10-13): initial version, adapted from
  :mod:`sage.rings.real_mpfr` and :mod:`sage.rings.complex_mpfr`.

- Alex Ghitza (2008-11): cache, generators, random element, and many doctests.

- Yann Laigle-Chapuy (2010-01): improves compatibility with CC, updates.

- Jeroen Demeyer (2012-02): reformat documentation, make MPC a standard
  package.

- Travis Scrimshaw (2012-10-18): Added doctests for full coverage.

- Vincent Klein (2017-11-15) : add __mpc__() to class MPComplexNumber.
  MPComplexNumber constructor support gmpy2.mpz, gmpy2.mpq, gmpy2.mpfr
  and gmpy2.mpc parameters.

EXAMPLES::

    sage: MPC = MPComplexField(42)
    sage: a = MPC(12, '15.64E+32'); a
    12.0000000000 + 1.56400000000e33*I
    sage: a *a *a *a
    5.98338564121e132 - 1.83633318912e101*I
    sage: a + 1
    13.0000000000 + 1.56400000000e33*I
    sage: a / 3
    4.00000000000 + 5.21333333333e32*I
    sage: MPC("infinity + NaN *I")
    +infinity + NaN*I
"""
# ****************************************************************************
#       Copyright (C) 2008 Philippe Theveny <thevenyp@loria.fr>
#                     2008 Alex Ghitza
#                     2010 Yann Laigle-Chapuy
#                     2012 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import re
from sage.rings import real_mpfr
import weakref
from cpython.object cimport Py_NE

import sage
from sage.cpython.string cimport str_to_bytes
from sage.libs.mpfr cimport *
from sage.libs.mpc cimport *
from sage.structure.parent cimport Parent
from sage.structure.element cimport Element
from sage.structure.richcmp cimport rich_to_bool
from sage.categories.map cimport Map

try:
    from sage.libs.pari import pari
    from cypari2.handle_error import PariError
    from cypari2.gen import Gen as pari_gen
except ImportError:
    pari_gen = PariError = ()

from sage.rings.integer cimport Integer
from sage.rings.complex_mpfr cimport ComplexNumber
from sage.rings.complex_mpfr import ComplexField_class

from sage.misc.randstate cimport randstate, current_randstate
from sage.rings.real_mpfr cimport RealField_class, RealNumber
from sage.rings.real_mpfr import mpfr_prec_min, mpfr_prec_max
from sage.structure.richcmp cimport rich_to_bool, richcmp
from sage.categories.fields import Fields

cimport gmpy2
gmpy2.import_gmpy2()

AA = None
QQbar = None
CDF = CLF = RLF = None


def late_import():
    """
    Import the objects/modules after build (when needed).

    TESTS::

        sage: sage.rings.complex_mpc.late_import()
    """
    global AA, QQbar
    global CLF, RLF, CDF
    if AA is None:
        import sage.rings.qqbar
        AA = sage.rings.qqbar.AA
        QQbar = sage.rings.qqbar.QQbar
        from sage.rings.real_lazy import CLF, RLF
        from sage.rings.complex_double import CDF


_mpfr_rounding_modes = ['RNDN', 'RNDZ', 'RNDU', 'RNDD']

_mpc_rounding_modes = [ 'RNDNN', 'RNDZN', 'RNDUN', 'RNDDN',
                        '', '', '', '', '', '', '', '', '', '', '', '',
                        'RNDNZ', 'RNDZZ', 'RNDUZ', 'RNDDZ',
                        '', '', '', '', '', '', '', '', '', '', '', '',
                        'RNDUN', 'RNDZU', 'RNDUU', 'RNDDU',
                        '', '', '', '', '', '', '', '', '', '', '', '',
                        'RNDDN', 'RNDZD', 'RNDUD', 'RNDDD' ]

cdef inline mpfr_rnd_t rnd_re(mpc_rnd_t rnd) noexcept:
    """
    Return the numeric value of the real part rounding mode.  This
    is an internal function.
    """
    return <mpfr_rnd_t>(rnd & 3)

cdef inline mpfr_rnd_t rnd_im(mpc_rnd_t rnd) noexcept:
    """
    Return the numeric value of the imaginary part rounding mode.
    This is an internal function.
    """
    return <mpfr_rnd_t>(rnd >> 4)

sign = '[+-]'
digit_ten = '[0123456789]'
exponent_ten = '[e@]' + sign + '?[0123456789]+'
number_ten = r'inf(?:inity)?|@inf@|nan(?:\([0-9A-Z_]*\))?|@nan@(?:\([0-9A-Z_]*\))?'\
    '|(?:' + digit_ten + r'*\.' + digit_ten + '+|' + digit_ten + r'+\.?)(?:' + exponent_ten + ')?'
imaginary_ten = r'i(?:\s*\*\s*(?:' + number_ten + '))?|(?:' + number_ten + r')\s*\*\s*i'
complex_ten = '(?P<im_first>(?P<im_first_im_sign>' + sign + r')?\s*(?P<im_first_im_abs>' + imaginary_ten + r')' \
    r'(\s*(?P<im_first_re_sign>' + sign + r')\s*(?P<im_first_re_abs>' + number_ten + '))?)' \
    '|' \
    '(?P<re_first>(?P<re_first_re_sign>' + sign + r')?\s*(?P<re_first_re_abs>' + number_ten + r')' \
    r'(\s*(?P<re_first_im_sign>' + sign + r')\s*(?P<re_first_im_abs>' + imaginary_ten + '))?)'
re_complex_ten = re.compile(r'^\s*(?:' + complex_ten + r')\s*$', re.I)

cpdef inline split_complex_string(string, int base=10):
    """
    Split and return in that order the real and imaginary parts
    of a complex in a string.

    This is an internal function.

    EXAMPLES::

        sage: sage.rings.complex_mpc.split_complex_string('123.456e789')
        ('123.456e789', None)
        sage: sage.rings.complex_mpc.split_complex_string('123.456e789*I')
        (None, '123.456e789')
        sage: sage.rings.complex_mpc.split_complex_string('123.+456e789*I')
        ('123.', '+456e789')
        sage: sage.rings.complex_mpc.split_complex_string('123.456e789', base=2)
        (None, None)
    """
    if base == 10:
        number = number_ten
        z = re_complex_ten.match(string)
    else:
        all_digits = "0123456789abcdefghijklmnopqrstuvwxyz"
        digit = '[' + all_digits[0:base] + ']'

        # In MPFR, '1e42'-> 10^42, '1p42'->2^42, '1@42'->base^42
        if base == 2:
            exponent = '[e@p]'
        elif base <= 10:
            exponent = '[e@]'
        elif base == 16:
            exponent = '[@p]'
        else:
            exponent = '@'
            exponent +=  sign + '?' + digit + '+'

        # Warning: number, imaginary, and complex should be enclosed in parentheses
        # when used as regexp because of alternatives '|'
        number = r'@nan@(?:\([0-9A-Z_]*\))?|@inf@|(?:' + digit + r'*\.' + digit + '+|' + digit + r'+\.?)(?:' + exponent + ')?'
        if base <= 10:
            number = r'nan(?:\([0-9A-Z_]*\))?|inf(?:inity)?|' + number
        imaginary = r'i(?:\s*\*\s*(?:' + number + '))?|(?:' + number + r')\s*\*\s*i'
        complex = '(?P<im_first>(?P<im_first_im_sign>' + sign + r')?\s*(?P<im_first_im_abs>' + imaginary + r')' \
            r'(\s*(?P<im_first_re_sign>' + sign + r')\s*(?P<im_first_re_abs>' + number + '))?)' \
            '|' \
            '(?P<re_first>(?P<re_first_re_sign>' + sign + r')?\s*(?P<re_first_re_abs>' + number + r')' \
            r'(\s*(?P<re_first_im_sign>' + sign + r')\s*(?P<re_first_im_abs>' + imaginary + '))?)'

        z = re.match(r'^\s*(?:' + complex + r')\s*$', string, re.I)

    x, y = None, None
    if z is not None:
        if z.group('im_first') is not None:
            prefix = 'im_first'
        elif z.group('re_first') is not None:
            prefix = 're_first'
        else:
            return None

        if z.group(prefix + '_re_abs') is not None:
            x = z.expand(r'\g<' + prefix + '_re_abs>')
            if z.group(prefix + '_re_sign') is not None:
                x = z.expand(r'\g<' + prefix + '_re_sign>') + x

        if z.group(prefix + '_im_abs') is not None:
            y = re.search('(?P<im_part>' + number + ')', z.expand(r'\g<' + prefix + '_im_abs>'), re.I)
            if y is None:
                y = '1'
            else:
                y = y.expand(r'\g<im_part>')
            if z.group(prefix + '_im_sign') is not None:
                y = z.expand(r'\g<' + prefix + '_im_sign>') + y

    return x, y

#*****************************************************************************
#
#       MPComplex Field
#
#*****************************************************************************
# The complex field is in Cython, so mpc elements will have access to
# their parent via direct C calls, which will be faster.

cache = {}


def MPComplexField(prec=53, rnd='RNDNN', names=None):
    """
    Return the complex field with real and imaginary parts having
    prec *bits* of precision.

    EXAMPLES::

        sage: MPComplexField()
        Complex Field with 53 bits of precision
        sage: MPComplexField(100)
        Complex Field with 100 bits of precision
        sage: MPComplexField(100).base_ring()
        Real Field with 100 bits of precision
        sage: i = MPComplexField(200).gen()
        sage: i^2
        -1.0000000000000000000000000000000000000000000000000000000000
    """
    global cache
    mykey = (prec, rnd)
    if mykey in cache:
        X = cache[mykey]
        C = X()
        if C is not None:
            return C
    C = MPComplexField_class(prec, rnd)
    cache[mykey] = weakref.ref(C)
    return C


cdef class MPComplexField_class(Field):
    def __init__(self, int prec=53, rnd="RNDNN"):
        """
        Initialize ``self``.

        INPUT:

        - ``prec`` -- integer (default: 53); precision

          prec is the number of bits used to represent the mantissa of
          both the real and imaginary part of complex floating-point number.

        - ``rnd`` -- string (default: ``'RNDNN'``); the rounding mode

          Rounding mode is of the form ``'RNDxy'`` where ``x`` and ``y`` are
          the rounding mode for respectively the real and imaginary parts and
          are one of:

          - ``'N'`` for rounding to nearest
          - ``'Z'`` for rounding towards zero
          - ``'U'`` for rounding towards plus infinity
          - ``'D'`` for rounding towards minus infinity

          For example, ``'RNDZU'`` indicates to round the real part towards
          zero, and the imaginary part towards plus infinity.

        EXAMPLES::

            sage: MPComplexField(17)
            Complex Field with 17 bits of precision
            sage: MPComplexField()
            Complex Field with 53 bits of precision
            sage: MPComplexField(1042,'RNDDZ')
            Complex Field with 1042 bits of precision and rounding RNDDZ

        ALGORITHMS: Computations are done using the MPC library.

        TESTS::

            sage: TestSuite(MPComplexField(17)).run()

            sage: MPComplexField(17).is_finite()
            False
        """
        if prec < mpfr_prec_min() or prec > mpfr_prec_max():
            raise ValueError("prec (=%s) must be >= %s and <= %s." % (
                prec, mpfr_prec_min(), mpfr_prec_max()))
        self._prec = prec
        if not isinstance(rnd, str):
            raise TypeError("rnd must be a string")
        try:
            n = _mpc_rounding_modes.index(rnd)
        except ValueError:
            raise ValueError("rnd (=%s) must be of the form RNDxy"
                             "where x and y are one of N, Z, U, D" % rnd)
        self.__rnd = n
        self.__rnd_str = rnd

        self.__real_field = real_mpfr.RealField(prec, rnd=_mpfr_rounding_modes[rnd_re(n)])
        self.__imag_field = real_mpfr.RealField(prec, rnd=_mpfr_rounding_modes[rnd_im(n)])

        Parent.__init__(self, self._real_field(), names=('I',), normalize=False,
                        category=Fields().Infinite())
        self._populate_coercion_lists_(coerce_list=[MPFRtoMPC(self._real_field(), self)])

    cdef MPComplexNumber _new(self):
        """
        Return a new complex number with parent ``self``.
        """
        cdef MPComplexNumber z
        z = MPComplexNumber.__new__(MPComplexNumber)
        z._parent = self
        mpc_init2(z.value, self._prec)
        z.init = 1
        return z

    def _repr_ (self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: MPComplexField(200, 'RNDDU') # indirect doctest
            Complex Field with 200 bits of precision and rounding RNDDU
        """
        s = "Complex Field with %s bits of precision"%self._prec
        if self.__rnd != MPC_RNDNN:
            s = s + " and rounding %s"%(self.__rnd_str)
        return s

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: MPC = MPComplexField(10)
            sage: latex(MPC) # indirect doctest
            \C
        """
        return "\\C"

    def __call__(self, x, im=None):
        """
        Create a floating-point complex using ``x`` and optionally an imaginary
        part ``im``.

        EXAMPLES::

            sage: MPC = MPComplexField()
            sage: MPC(2) # indirect doctest
            2.00000000000000
            sage: MPC(0, 1) # indirect doctest
            1.00000000000000*I
            sage: MPC(1, 1)
            1.00000000000000 + 1.00000000000000*I
            sage: MPC(2, 3)
            2.00000000000000 + 3.00000000000000*I
        """
        if x is None:
            return self.zero()
        # We implement __call__ to gracefully accept the second argument.
        if im is not None:
            x = x, im
        return Parent.__call__(self, x)

    def _element_constructor_(self, z):
        """
        Coerce `z` into this complex field.

        EXAMPLES::

            sage: C20 = MPComplexField(20) # indirect doctest

        The value can be set with a couple of reals::

            sage: a = C20(1.5625, 17.42); a
            1.5625 + 17.420*I
            sage: a.str(2)
            '1.1001000000000000000 + 10001.011010111000011*I'
            sage: C20(0, 2)
            2.0000*I

        Complex number can be coerced into MPComplexNumber::

            sage: C20(14.7+0.35*I)
            14.700 + 0.35000*I
            sage: C20(i*4, 7)
            Traceback (most recent call last):
            ...
            TypeError: unable to coerce to a ComplexNumber: <class 'sage.rings.number_field.number_field_element_quadratic.NumberFieldElement_gaussian'>

        Each part can be set with strings (written in base ten)::

            sage: C20('1.234', '56.789')
            1.2340 + 56.789*I

        The string can represent the whole complex value::

            sage: C20('42 + I * 100')
            42.000 + 100.00*I
            sage: C20('-42 * I')
            - 42.000*I

        The imaginary part can be written first::

            sage: C20('100*i+42')
            42.000 + 100.00*I

        Use ``'inf'`` for infinity and ``'nan'`` for Not a Number::

            sage: C20('nan+inf*i')
            NaN + +infinity*I
        """
        cdef MPComplexNumber zz
        zz = self._new()
        zz._set(z)
        return zz

    cpdef _coerce_map_from_(self, S):
        """
        Canonical coercion of `z` to this mpc complex field.

        The rings that canonically coerce to this mpc complex field are:

        - any mpc complex field with precision that is as large as this one
        - anything that canonically coerces to the mpfr real
          field with this prec and the rounding mode of real part.

        EXAMPLES::

            sage: MPComplexField(100)(17, '4.2') + MPComplexField(20)('6.0', -23) # indirect doctest
            23.000 - 18.800*I
            sage: a = MPComplexField(100)(17, '4.2') + MPComplexField(20)('6.0', -23)
            sage: a.parent()
            Complex Field with 20 bits of precision
        """
        if isinstance(S, RealField_class):
            return MPFRtoMPC(S, self)
        if isinstance(S, sage.rings.integer_ring.IntegerRing_class):
            return INTEGERtoMPC(S, self)

        RR = self.__real_field
        if RR.has_coerce_map_from(S):
            return self._coerce_map_via([RR], S)

        if isinstance(S, MPComplexField_class) and S.prec() >= self._prec:
            #FIXME: What map when rounding modes differ but prec is the same ?
            #       How to provide commutativity of morphisms ?
            #       Change _cmp_ when done
            return MPCtoMPC(S, self)

        if isinstance(S, ComplexField_class) and S.prec() >= self._prec:
            return CCtoMPC(S, self)

        late_import()
        if S in [AA, QQbar, CLF, RLF] or (S == CDF and self._prec <= 53):
            return self._generic_coerce_map(S)

        return self._coerce_map_via([CLF], S)

    def __reduce__(self):
        """
        For pickling.

        EXAMPLES::

            sage: C = MPComplexField(prec=200, rnd='RNDDZ')
            sage: loads(dumps(C)) == C
            True
        """
        return __create__MPComplexField_version0, (self._prec, self.__rnd_str)

    def __richcmp__(left, right, int op):
        """
        Compare ``self`` and ``other``, ignoring the rounding mode.

        EXAMPLES::

            sage: MPComplexField(10) == MPComplexField(11) # indirect doctest
            False
            sage: MPComplexField(10) == MPComplexField(10)
            True
            sage: MPComplexField(10,rnd='RNDZN') == MPComplexField(10,rnd='RNDZU')
            True
        """
        if left is right:
            return rich_to_bool(op, 0)

        if not isinstance(right, MPComplexField_class):
            return op == Py_NE

        cdef MPComplexField_class s = <MPComplexField_class>left
        cdef MPComplexField_class o = <MPComplexField_class>right
        return richcmp(s._prec, o._prec, op)

    def gen(self, n=0):
        """
        Return the generator of this complex field over its real subfield.

        EXAMPLES::

            sage: MPComplexField(34).gen()
            1.00000000*I
        """
        if n != 0:
            raise IndexError("n must be 0")
        return self(0, 1)

    def ngens(self):
        """
        Return 1, the number of generators of this complex field over its real
        subfield.

        EXAMPLES::

            sage: MPComplexField(34).ngens()
            1
        """
        return 1

    cpdef _an_element_(self):
        """
        Return an element of this complex field.

        EXAMPLES::

            sage: MPC = MPComplexField(20)
            sage: MPC._an_element_()
            1.0000*I
        """
        return self(0, 1)

    def random_element(self, min=0, max=1):
        """
        Return a random complex number, uniformly distributed with
        real and imaginary parts between min and max (default 0 to 1).

        EXAMPLES::

            sage: MPComplexField(100).random_element(-5, 10)  # random
            1.9305310520925994224072377281 + 0.94745292506956219710477444855*I
            sage: MPComplexField(10).random_element()  # random
            0.12 + 0.23*I
        """
        cdef MPComplexNumber z
        z = self._new()
        cdef randstate rstate = current_randstate()
        mpc_urandom(z.value, rstate.gmp_state)
        if min == 0 and max == 1:
            return z
        else:
            return (max-min)*z + min*self(1,1)

    cpdef bint is_exact(self) except -2:
        """
        Return whether or not this field is exact, which is always ``False``.

        EXAMPLES::

            sage: MPComplexField(42).is_exact()
            False
        """
        return False

    def characteristic(self):
        """
        Return 0, since the field of complex numbers has characteristic 0.

        EXAMPLES::

            sage: MPComplexField(42).characteristic()
            0
        """
        return Integer(0)

    def name(self):
        """
        Return the name of the complex field.

        EXAMPLES::

            sage: C = MPComplexField(10, 'RNDNZ'); C.name()
            'MPComplexField10_RNDNZ'
        """
        return "MPComplexField%s_%s"%(self._prec, self.__rnd_str)

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: MPC = MPComplexField()
            sage: hash(MPC) % 2^32 == hash(MPC.name()) % 2^32
            True
        """
        return hash(self.name())

    def prec(self):
        """
        Return the precision of this field of complex numbers.

        EXAMPLES::

            sage: MPComplexField().prec()
            53
            sage: MPComplexField(22).prec()
            22
        """
        return self._prec

    def rounding_mode(self):
        """
        Return rounding modes used for each part of a complex number.

        EXAMPLES::

            sage: MPComplexField().rounding_mode()
            'RNDNN'
            sage: MPComplexField(rnd='RNDZU').rounding_mode()
            'RNDZU'
        """
        return self.__rnd_str

    def rounding_mode_real(self):
        """
        Return rounding mode used for the real part of complex number.

        EXAMPLES::

            sage: MPComplexField(rnd='RNDZU').rounding_mode_real()
            'RNDZ'
        """
        return _mpfr_rounding_modes[rnd_re(self.__rnd)]

    def rounding_mode_imag(self):
        """
        Return rounding mode used for the imaginary part of complex number.

        EXAMPLES::

            sage: MPComplexField(rnd='RNDZU').rounding_mode_imag()
            'RNDU'
        """
        return _mpfr_rounding_modes[rnd_im(self.__rnd)]

    def _real_field(self):
        """
        Return real field for the real part.

        EXAMPLES::

            sage: MPComplexField()._real_field()
            Real Field with 53 bits of precision
        """
        return self.__real_field

    def _imag_field(self):
        """
        Return real field for the imaginary part.

        EXAMPLES::

            sage: MPComplexField(prec=100)._imag_field()
            Real Field with 100 bits of precision
        """
        return self.__imag_field


#*****************************************************************************
#
#       MPComplex Number -- element of MPComplex Field
#
#*****************************************************************************

cdef class MPComplexNumber(sage.structure.element.FieldElement):
    """
    A floating point approximation to a complex number using any specified
    precision common to both real and imaginary part.
    """
    cdef MPComplexNumber _new(self):
        """
        Return a new complex number with same parent as ``self``.
        """
        cdef MPComplexNumber z
        z = MPComplexNumber.__new__(MPComplexNumber)
        z._parent = self._parent
        mpc_init2(z.value, (<MPComplexField_class>self._parent)._prec)
        z.init = 1
        return z

    def __init__(self, MPComplexField_class parent, x, y=None, int base=10):
        """
        Create a complex number.

        INPUT:

        - ``x`` -- real part or the complex value in a string

        - ``y`` -- imaginary part

        - ``base`` -- when ``x`` or ``y`` is a string, base in which the
          number is written

        A :class:`MPComplexNumber` should be called by first creating a
        :class:`MPComplexField`, as illustrated in the examples.

        EXAMPLES::

            sage: C200 = MPComplexField(200)
            sage: C200(1/3, '0.6789')
            0.33333333333333333333333333333333333333333333333333333333333 + 0.67890000000000000000000000000000000000000000000000000000000*I
            sage: C3 = MPComplexField(3)
            sage: C3('1.2345', '0.6789')
            1.2 + 0.62*I
            sage: C3(3.14159)
            3.0

        Rounding modes::

            sage: w = C3(5/2, 7/2); w.str(2)
            '10.1 + 11.1*I'
            sage: MPComplexField(2, rnd="RNDZN")(w).str(2)
            '10. + 100.*I'
            sage: MPComplexField(2, rnd="RNDDU")(w).str(2)
            '10. + 100.*I'
            sage: MPComplexField(2, rnd="RNDUD")(w).str(2)
            '11. + 11.*I'
            sage: MPComplexField(2, rnd="RNDNZ")(w).str(2)
            '10. + 11.*I'

        TESTS::

            sage: MPComplexField(42)._repr_option('element_is_atomic')
            False
        """
        self.init = 0
        if parent is None:
            raise TypeError
        self._parent = parent
        mpc_init2(self.value, parent._prec)
        self.init = 1
        if x is None: return
        self._set(x, y, base)

    def _set(self, z, y=None, base=10):
        """
        EXAMPLES::

            sage: MPC = MPComplexField(100)
            sage: r = RealField(100).pi()
            sage: z = MPC(r); z # indirect doctest
            3.1415926535897932384626433833
            sage: MPComplexField(10, rnd='RNDDD')(z)
            3.1
            sage: c = ComplexField(53)(1, r)
            sage: MPC(c)
            1.0000000000000000000000000000 + 3.1415926535897931159979634685*I
            sage: MPC(I)
            1.0000000000000000000000000000*I
            sage: MPC('-0 +i')
            1.0000000000000000000000000000*I
            sage: MPC(1+i)
            1.0000000000000000000000000000 + 1.0000000000000000000000000000*I
            sage: MPC(1/3)
            0.33333333333333333333333333333

        ::

            sage: MPC(1, r/3)
            1.0000000000000000000000000000 + 1.0471975511965977461542144611*I
            sage: MPC(3, 2)
            3.0000000000000000000000000000 + 2.0000000000000000000000000000*I
            sage: MPC(0, r)
            3.1415926535897932384626433833*I
            sage: MPC('0.625e-26', '0.0000001')
            6.2500000000000000000000000000e-27 + 1.0000000000000000000000000000e-7*I

        Conversion from gmpy2 numbers::

            sage: from gmpy2 import *
            sage: MPC(mpc(int(2),int(1)))
            2.0000000000000000000000000000 + 1.0000000000000000000000000000*I
            sage: MPC(mpfr(2.5))
            2.5000000000000000000000000000
            sage: MPC(mpq('3/2'))
            1.5000000000000000000000000000
            sage: MPC(mpz(int(5)))
            5.0000000000000000000000000000
            sage: re = mpfr('2.5')
            sage: im = mpz(int(2))
            sage: MPC(re, im)
            2.5000000000000000000000000000 + 2.0000000000000000000000000000*I
        """
        # This should not be called except when the number is being created.
        # Complex Numbers are supposed to be immutable.
        cdef mpc_rnd_t rnd
        rnd =(<MPComplexField_class>self._parent).__rnd
        if y is None:
            if z is None: return
            if isinstance(z, MPComplexNumber):
                mpc_set(self.value, (<MPComplexNumber>z).value, rnd)
                return
            elif isinstance(z, str):
                a, b = split_complex_string(z, base)
                if a is None and b is None:
                    raise TypeError("unable to convert {!r} to a MPComplexNumber".format(z))
                # set real part
                if a is None:
                    mpfr_set_ui(self.value.re, 0, MPFR_RNDN)
                else:
                    mpfr_set_str(self.value.re, str_to_bytes(a),
                                 base, rnd_re(rnd))
                # set imag part
                if b is None:
                    mpfr_set_ui(self.value.im, 0, MPFR_RNDN)
                else:
                    mpfr_set_str(self.value.im, str_to_bytes(b),
                                 base, rnd_im(rnd))
                return
            elif isinstance(z, ComplexNumber):
                mpc_set_fr_fr(self.value, (<ComplexNumber>z).__re, (<ComplexNumber>z).__im, rnd)
                return
            elif isinstance(z, pari_gen):
                real, imag = z.real(), z.imag()
            elif isinstance(z, (list, tuple)):
                real, imag = z
            elif isinstance(z, complex):
                real, imag = z.real, z.imag
            elif type(z) is gmpy2.mpc:
                mpc_set(self.value, (<gmpy2.mpc>z).c, rnd)
                return
            elif type(z) is gmpy2.mpfr:
                mpc_set_fr(self.value, (<gmpy2.mpfr>z).f, rnd)
                return
            elif type(z) is gmpy2.mpq:
                mpc_set_q(self.value, (<gmpy2.mpq>z).q, rnd)
                return
            elif type(z) is gmpy2.mpz:
                mpc_set_z(self.value, (<gmpy2.mpz>z).z, rnd)
                return
            elif isinstance(z, RealNumber):
                zz = sage.rings.real_mpfr.RealField(self._parent.prec())(z)
                mpc_set_fr(self.value, (<RealNumber>zz).value, rnd)
                return
            elif isinstance(z, Integer):
                mpc_set_z(self.value, (<Integer>z).value, rnd)
                return
            elif isinstance(z, int):
                mpc_set_si(self.value, z, rnd)
                return
            else:
                C = sage.rings.complex_mpfr.ComplexField(self._parent.prec())
                if C.has_coerce_map_from(z.parent()) or isinstance(z, sage.symbolic.expression.Expression):
                    self._set(C(z))
                    return
                else:
                    real = z
                    imag = 0
        else:
            real = z
            imag = y

        cdef RealField_class R = self._parent._real_field()
        try:
            rr = R(real)
            ii = R(imag)
            mpc_set_fr_fr(self.value, (<RealNumber>rr).value, (<RealNumber>ii).value, rnd)

        except TypeError:
            raise TypeError("unable to coerce to a ComplexNumber: %s" % type(real))

    def __reduce__(self):
        """
        For pickling.

        EXAMPLES::

            sage: C = MPComplexField(prec=200, rnd='RNDUU')
            sage: b = C(393.39203845902384098234098230948209384028340)
            sage: loads(dumps(b)) == b
            True
            sage: b = C(-1).sqrt(); b
            1.0000000000000000000000000000000000000000000000000000000000*I
            sage: loads(dumps(b)) == b
            True

        Some tests with ``NaN``, which cannot be compared to anything::

            sage: b = C(1)/C(0); b
            NaN + NaN*I
            sage: loads(dumps(b))
            NaN + NaN*I
            sage: b = C(-1)/C(0.); b
            NaN + NaN*I
            sage: loads(dumps(b))
            NaN + NaN*I
        """
        s = self.str(32)
        return (__create_MPComplexNumber_version0, (self._parent, s, 32))

    def __dealloc__(self):
        if self.init:
            mpc_clear(self.value)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: MPComplexField()(2, -3) # indirect doctest
            2.00000000000000 - 3.00000000000000*I
        """
        return self.str(truncate=True)

    def _latex_(self):
        """
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: latex(MPComplexField()(2, -3)) # indirect doctest
            2.00000000000000 - 3.00000000000000i
        """
        import re
        s = repr(self).replace('*I', 'i')
        return re.sub(r"e(-?\d+)", r" \\times 10^{\1}", s)

    def __hash__(self):
        """
        Return the hash of ``self``, which coincides with the python
        complex and float (and often int) types.

        This has the drawback that two very close high precision
        numbers will have the same hash, but allows them to play
        nicely with other real types.

        EXAMPLES::

            sage: hash(MPComplexField()('1.2', 33)) == hash(complex(1.2, 33))
            True
        """
        return hash(complex(self))

    def __getitem__(self, i):
        r"""
        Return either the real or imaginary component of ``self``
        depending on the choice of ``i``: real (``i``=0), imaginary (``i``=1).

        INPUT:

        - ``i`` -- 0 or 1

          - ``0`` -- will return the real component of ``self``
          - ``1`` -- will return the imaginary component of ``self``

        EXAMPLES::

            sage: MPC = MPComplexField()
            sage: a = MPC(2,1)
            sage: a.__getitem__(0)
            2.00000000000000
            sage: a.__getitem__(1)
            1.00000000000000

        ::

            sage: b = MPC(42,0)
            sage: b
            42.0000000000000
            sage: b.__getitem__(1)
            0.000000000000000
        """
        if i == 0:
            return self.real()
        elif i == 1:
            return self.imag()
        raise IndexError("i must be between 0 and 1.")

    def prec(self):
        """
        Return precision of this complex number.

        EXAMPLES::

            sage: i = MPComplexField(2000).0
            sage: i.prec()
            2000
        """
        return <MPComplexField_class>(self._parent)._prec

    def real(self):
        """
        Return the real part of ``self``.

        EXAMPLES::

            sage: C = MPComplexField(100)
            sage: z = C(2, 3)
            sage: x = z.real(); x
            2.0000000000000000000000000000
            sage: x.parent()
            Real Field with 100 bits of precision
        """
        cdef RealNumber x
        x = RealNumber(self._parent._real_field())
        mpfr_set (x.value, self.value.re, (<RealField_class>x._parent).rnd)
        return x

    def imag(self):
        """
        Return imaginary part of ``self``.

        EXAMPLES::

            sage: C = MPComplexField(100)
            sage: z = C(2, 3)
            sage: x = z.imag(); x
            3.0000000000000000000000000000
            sage: x.parent()
            Real Field with 100 bits of precision
        """
        cdef RealNumber y
        y = RealNumber(self._parent._imag_field())
        mpfr_set (y.value, self.value.im, (<RealField_class>y._parent).rnd)
        return y

    def str(self, base=10, **kwds):
        """
        Return a string of ``self``.

        INPUT:

        - ``base`` -- (default: 10) base for output

        - ``**kwds`` -- other arguments to pass to the ``str()``
          method of the real numbers in the real and imaginary parts

        EXAMPLES::

            sage: MPC = MPComplexField(64)
            sage: z = MPC(-4, 3)/7
            sage: z.str()
            '-0.571428571428571428564 + 0.428571428571428571436*I'
            sage: z.str(16)
            '-0.92492492492492490 + 0.6db6db6db6db6db70*I'
            sage: z.str(truncate=True)
            '-0.571428571428571429 + 0.428571428571428571*I'
            sage: z.str(2)
            '-0.1001001001001001001001001001001001001001001001001001001001001001 + 0.01101101101101101101101101101101101101101101101101101101101101110*I'
        """
        s = ""
        if self.real() != 0:
            s = self.real().str(base, **kwds)
        if self.imag() != 0:
            if mpfr_signbit(self.value.im):
                s += ' - ' + (-self.imag()).str(base, **kwds) + '*I'
            else:
                if s:
                    s += ' + '
                s += self.imag().str(base, **kwds) + '*I'
        if not s:
            return "0"
        return s

    def __copy__(self):
        """
        Return copy of ``self``.

        Since ``self`` is immutable, we just return ``self`` again.

        EXAMPLES::

            sage: a = MPComplexField()(3.5, 3)
            sage: copy(a) is  a
            True
        """
        return self    # since object is immutable.

    def __deepcopy__(self, memo):
        """
        EXAMPLES::

            sage: a = MPComplexField()(3.5, 3)
            sage: deepcopy(a) is  a
            True
        """
        return self    # since object is immutable.

    def __int__(self):
        r"""
        Method for converting ``self`` to type ``int``.

        Called by the ``int`` function. Note that calling this method returns
        an error since, in general, complex numbers cannot be coerced into
        integers.

        EXAMPLES::

            sage: MPC = MPComplexField()
            sage: a = MPC(2,1)
            sage: int(a)
            Traceback (most recent call last):
            ...
            TypeError: can...t convert complex to int; use int(abs(z))
            sage: a.__int__()
            Traceback (most recent call last):
            ...
            TypeError: can...t convert complex to int; use int(abs(z))
        """
        raise TypeError("can't convert complex to int; use int(abs(z))")

    def __float__(self):
        r"""
        Method for converting ``self`` to type ``float``.

        Called by the ``float`` function. Note that calling this method returns
        an error since if the imaginary part of the number is not zero.

        EXAMPLES::

            sage: MPC = MPComplexField()
            sage: a = MPC(1, 0)
            sage: float(a)
            1.0
            sage: a = MPC(2,1)
            sage: float(a)
            Traceback (most recent call last):
            ...
            TypeError: can...t convert complex to float; use abs(z)
            sage: a.__float__()
            Traceback (most recent call last):
            ...
            TypeError: can...t convert complex to float; use abs(z)
        """
        if mpfr_zero_p(self.value.im):
            return mpfr_get_d(self.value.re,
                              rnd_re((<MPComplexField_class>self._parent).__rnd))
        else:
            raise TypeError("can't convert complex to float; use abs(z)")

    def __complex__(self):
        r"""
        Method for converting ``self`` to type ``complex``.

        Called by the ``complex`` function.

        EXAMPLES::

            sage: MPC = MPComplexField()
            sage: a = MPC(2,1)
            sage: complex(a)
            (2+1j)
            sage: type(complex(a))
            <... 'complex'>
            sage: a.__complex__()
            (2+1j)
        """
        # Fixme: is it the right choice for rounding modes ?
        cdef mpc_rnd_t rnd
        rnd = (<MPComplexField_class>self._parent).__rnd
        return complex(mpfr_get_d(self.value.re, rnd_re(rnd)), mpfr_get_d(self.value.im, rnd_im(rnd)))

    def __pari__(self):
        r"""
        Convert ``self`` to a PARI object.

        OUTPUT: a PARI ``t_COMPLEX`` object if the input is not purely
        real. If the input is real, a ``t_REAL`` is returned.

        EXAMPLES::

            sage: MPC = MPComplexField()
            sage: a = MPC(2,1)
            sage: a.__pari__()                                                          # needs sage.libs.pari
            2.00000000000000 + 1.00000000000000*I
            sage: pari(a)                                                               # needs sage.libs.pari
            2.00000000000000 + 1.00000000000000*I
            sage: pari(a).type()                                                        # needs sage.libs.pari
            't_COMPLEX'
            sage: a = MPC(pi)                                                           # needs sage.libs.pari sage.symbolic
            sage: pari(a)                                                               # needs sage.libs.pari sage.symbolic
            3.14159265358979
            sage: pari(a).type()                                                        # needs sage.libs.pari sage.symbolic
            't_REAL'
            sage: a = MPC(-2).sqrt()
            sage: pari(a)                                                               # needs sage.libs.pari
            1.41421356237310*I

        The precision is preserved, rounded up to the wordsize::

            sage: MPC = MPComplexField(250)
            sage: MPC(1,2).__pari__().bitprecision()                                    # needs sage.libs.pari
            256
            sage: MPC(pi).__pari__().bitprecision()                                     # needs sage.libs.pari
            256
        """
        if mpfr_zero_p(self.value.re):
            re = pari.PARI_ZERO
        else:
            re = self.real().__pari__()
        if mpfr_zero_p(self.value.im):
            return re
        im = self.imag().__pari__()
        return pari.complex(re, im)

    def __mpc__(self):
        """
        Convert Sage ``MPComplexNumber`` to gmpy2 ``mpc``.

        EXAMPLES::

            sage: MPC = MPComplexField()
            sage: c = MPC(2,1)
            sage: c.__mpc__()
            mpc('2.0+1.0j')
            sage: from gmpy2 import mpc
            sage: mpc(c)
            mpc('2.0+1.0j')
            sage: MPCF = MPComplexField(42)
            sage: mpc(MPCF(12, 12)).precision
            (42, 42)
            sage: MPCF = MPComplexField(236)
            sage: mpc(MPCF(12, 12)).precision
            (236, 236)
            sage: MPCF = MPComplexField(63)
            sage: x = MPCF('15.64E+128', '15.64E+128')
            sage: y = mpc(x)
            sage: y.precision
            (63, 63)
            sage: MPCF(y) == x
            True
            sage: x = mpc('1.324+4e50j', precision=(70,70))
            sage: y = MPComplexField(70)(x)
            sage: mpc(y) == x
            True
        """
        return gmpy2.GMPy_MPC_From_mpfr(self.value.re, self.value.im)

    cpdef _richcmp_(self, other, int op):
        r"""
        Compare ``self`` to ``other``.

        EXAMPLES::

            sage: MPC = MPComplexField()
            sage: a = MPC(2,1)
            sage: b = MPC(1,2)
            sage: a < b
            False
            sage: a > b
            True
        """
        cdef MPComplexNumber z = <MPComplexNumber>other
        # NaN should compare to nothing
        if mpfr_nan_p(self.value.re) or mpfr_nan_p(self.value.im) or mpfr_nan_p(z.value.re) or mpfr_nan_p(z.value.im):
            return False
        cdef int c = mpc_cmp(self.value, z.value)
        cdef int cre = MPC_INEX_RE(c)
        cdef int cim
        if cre:
            return rich_to_bool(op, cre)
        else:
            cim = MPC_INEX_IM(c)
            return rich_to_bool(op, cim)

    def __bool__(self):
        """
        Return ``True`` if ``self`` is not zero.

        This is an internal function; use :meth:`is_zero()` instead.

        EXAMPLES::

            sage: MPC = MPComplexField()
            sage: z = 1 + MPC(I)
            sage: z.is_zero()
            False
        """
        return not (mpfr_zero_p(self.value.re) and mpfr_zero_p(self.value.im))

    def is_square(self):
        r"""
        This function always returns true as `\CC` is algebraically closed.

        EXAMPLES::

            sage: C200 = MPComplexField(200)
            sage: a = C200(2,1)
            sage: a.is_square()
            True

        `\CC` is algebraically closed, hence every element is a square::

            sage: b = C200(5)
            sage: b.is_square()
            True
        """
        return True

    def is_real(self):
        """
        Return ``True`` if ``self`` is real, i.e. has imaginary part zero.

        EXAMPLES::

            sage: C200 = MPComplexField(200)
            sage: C200(1.23).is_real()
            True
            sage: C200(1+i).is_real()
            False
        """
        return (mpfr_zero_p(self.value.im) != 0)

    def is_imaginary(self):
        """
        Return ``True`` if ``self`` is imaginary, i.e. has real part zero.

        EXAMPLES::

            sage: C200 = MPComplexField(200)
            sage: C200(1.23*i).is_imaginary()
            True
            sage: C200(1+i).is_imaginary()
            False
        """
        return (mpfr_zero_p(self.value.re) != 0)

    def algebraic_dependency(self, n, **kwds):
        """
        Return an irreducible polynomial of degree at most `n` which is
        approximately satisfied by this complex number.

        ALGORITHM: Uses the PARI C-library :pari:`algdep` command.

        INPUT: Type ``algdep?`` at the top level prompt. All additional
        parameters are passed onto the top-level algdep command.

        EXAMPLES::

            sage: MPC = MPComplexField()
            sage: z = (1/2)*(1 + sqrt(3.0) * MPC.0); z
            0.500000000000000 + 0.866025403784439*I
            sage: p = z.algebraic_dependency(5)
            sage: p
            x^2 - x + 1
            sage: p(z)
            1.11022302462516e-16
        """
        from sage.arith.misc import algdep
        return algdep(self, n, **kwds)

    ################################
    # Basic Arithmetic
    ################################

    cpdef _add_(self, right):
        """
        Add two complex numbers with the same parent.

        EXAMPLES::

           sage: MPC = MPComplexField(30)
           sage: MPC(-1.5, 2) + MPC(0.2, 1) # indirect doctest
           -1.3000000 + 3.0000000*I
        """
        cdef MPComplexNumber z
        z = self._new()
        mpc_add(z.value, self.value, (<MPComplexNumber>right).value, (<MPComplexField_class>self._parent).__rnd)
        return z

    cpdef _sub_(self, right):
        """
        Subtract two complex numbers with the same parent.

        EXAMPLES::

           sage: MPC = MPComplexField(30)
           sage: MPC(-1.5, 2) - MPC(0.2, 1) # indirect doctest
           -1.7000000 + 1.0000000*I
        """
        cdef MPComplexNumber z
        z = self._new()
        mpc_sub(z.value, self.value, (<MPComplexNumber>right).value, (<MPComplexField_class>self._parent).__rnd)
        return z

    cpdef _mul_(self, right):
        """
        Multiply two complex numbers with the same parent.

        EXAMPLES::

           sage: MPC = MPComplexField(30)
           sage: MPC(-1.5, 2) * MPC(0.2, 1) # indirect doctest
           -2.3000000 - 1.1000000*I
        """
        cdef MPComplexNumber z
        z = self._new()
        mpc_mul(z.value, self.value, (<MPComplexNumber>right).value, (<MPComplexField_class>self._parent).__rnd)
        return z

    cpdef _div_(self, right):
        """
        Divide two complex numbers with the same parent.

        EXAMPLES::

           sage: MPC = MPComplexField(30)
           sage: MPC(-1.5, 2) / MPC(0.2, 1) # indirect doctest
           1.6346154 + 1.8269231*I
           sage: MPC(-1, 1) / MPC(0)
           NaN + NaN*I
        """
        cdef MPComplexNumber z, x
        z = self._new()
        x = <MPComplexNumber>right
        if not x.is_zero():
            mpc_div(z.value, self.value, x.value, (<MPComplexField_class>self._parent).__rnd)
        return z

    cpdef _neg_(self):
        """
        Return the negative of this complex number.

        EXAMPLES::

            sage: MPC = MPComplexField(30)
            sage: - MPC(-1.5, 2) # indirect doctest
            1.5000000 - 2.0000000*I
            sage: - MPC(0)
            0
        """
        cdef MPComplexNumber z
        z = self._new()
        mpc_neg(z.value, self.value, (<MPComplexField_class>self._parent).__rnd)
        return z

    def __invert__(self):
        """
        Return the multiplicative inverse.

        EXAMPLES::

            sage: C = MPComplexField()
            sage: a = ~C(5, 1)
            sage: a * C(5, 1)
            1.00000000000000
        """
        cdef MPComplexNumber z
        z = self._new()
        if not self.is_zero():
            mpc_ui_div (z.value, 1, self.value, (<MPComplexField_class>self._parent).__rnd)
        return z

    def __neg__(self):
        r"""
        Return the negative of this complex number.

            -(a + ib) = -a -i b

        EXAMPLES::

            sage: MPC = MPComplexField()
            sage: a = MPC(2,1)
            sage: -a
            -2.00000000000000 - 1.00000000000000*I
            sage: a.__neg__()
            -2.00000000000000 - 1.00000000000000*I
        """
        cdef MPComplexNumber z
        z = self._new()
        mpc_neg(z.value, self.value, (<MPComplexField_class>self._parent).__rnd)
        return z

    def __abs__(self):
        r"""
        Absolute value or modulus of this complex number,
        rounded with the rounding mode of the real part:

        .. MATH::

            |a + ib| = \sqrt(a^2 + b^2).

        OUTPUT:

        A floating-point number in the real field of the real part
        (same precision, same rounding mode).

        EXAMPLES:

        Note that the absolute value of a complex number with imaginary
        component equal to zero is the absolute value of the real component::

            sage: MPC = MPComplexField()
            sage: a = MPC(2,1)
            sage: abs(a)
            2.23606797749979
            sage: a.__abs__()
            2.23606797749979
            sage: float(sqrt(2^2 + 1^1))                                                # needs sage.symbolic
            2.23606797749979

            sage: b = MPC(42,0)
            sage: abs(b)
            42.0000000000000
            sage: b.__abs__()
            42.0000000000000
            sage: b
            42.0000000000000
        """
        cdef RealNumber x
        x = RealNumber(self._parent._real_field())
        mpc_abs (x.value, self.value, (<RealField_class>x._parent).rnd)
        return x

    def norm(self):
        r"""
        Return the norm of a complex number, rounded with the rounding
        mode of the real part.  The norm is the square of the absolute
        value:

        .. MATH::

            \mathrm{norm}(a + ib) = a^2 + b^2.

        OUTPUT:

        A floating-point number in the real field of the real part
        (same precision, same rounding mode).

        EXAMPLES:

        This indeed acts as the square function when the imaginary
        component of ``self`` is equal to zero::

            sage: MPC = MPComplexField()
            sage: a = MPC(2,1)
            sage: a.norm()
            5.00000000000000
            sage: b = MPC(4.2,0)
            sage: b.norm()
            17.6400000000000
            sage: b^2
            17.6400000000000
        """
        cdef RealNumber x
        x = RealNumber(self._parent._real_field())
        mpc_norm(x.value, self.value, (<RealField_class>x._parent).rnd)
        return x

    def __pow__(self, right, modulus):
        """
        Compute ``self`` raised to the power of exponent, rounded in
        the direction specified by the parent of ``self``.

        .. TODO::

            FIXME: Branch cut

        EXAMPLES::

            sage: MPC.<i> = MPComplexField(20)
            sage: a = i^2; a
            -1.0000
            sage: a.parent()
            Complex Field with 20 bits of precision
            sage: a = (1+i)^i; a
            0.42883 + 0.15487*I
            sage: (1+i)^(1+i)
            0.27396 + 0.58370*I
            sage: a.parent()
            Complex Field with 20 bits of precision
            sage: i^i
            0.20788
            sage: (2+i)^(0.5)
            1.4553 + 0.34356*I
        """
        cdef MPComplexNumber z, x, p
        x = <MPComplexNumber>self
        z = x._new()

        if isinstance(right, int):
            mpc_pow_si(z.value, x.value, right, (<MPComplexField_class>x._parent).__rnd)
        elif isinstance(right, Integer):
            mpc_pow_z(z.value, x.value, (<Integer>right).value, (<MPComplexField_class>x._parent).__rnd)
        else:
            try:
                p = (<MPComplexField_class>x._parent)(right)
            except Exception:
                raise ValueError
            mpc_pow(z.value, x.value, p.value, (<MPComplexField_class>x._parent).__rnd)

        return z

    ################################
    # Trigonometric & hyperbolic functions
    ################################

    def cos(self):
        r"""
        Return the cosine of this complex number:

        .. MATH::

            \cos(a + ib) = \cos a \cosh b -i \sin a \sinh b.

        EXAMPLES::

            sage: MPC = MPComplexField()
            sage: u = MPC(2, 4)
            sage: cos(u)
            -11.3642347064011 - 24.8146514856342*I
        """
        cdef MPComplexNumber z
        z = self._new()
        mpc_cos (z.value, self.value, (<MPComplexField_class>self._parent).__rnd)
        return z

    def sin(self):
        r"""
        Return the sine of this complex number:

        .. MATH::

            \sin(a + ib) = \sin a \cosh b + i \cos x \sinh b.

        EXAMPLES::

            sage: MPC = MPComplexField()
            sage: u = MPC(2, 4)
            sage: sin(u)
            24.8313058489464 - 11.3566127112182*I
        """
        cdef MPComplexNumber z
        z = self._new()
        mpc_sin (z.value, self.value, (<MPComplexField_class>self._parent).__rnd)
        return z

    def tan(self):
        r"""
        Return the tangent of this complex number:

        .. MATH::

            \tan(a + ib) = (\sin 2a + i \sinh 2b)/(\cos 2a + \cosh 2b).

        EXAMPLES::

            sage: MPC = MPComplexField()
            sage: u = MPC(-2, 4)
            sage: tan(u)
            0.000507980623470039 + 1.00043851320205*I
        """
        cdef MPComplexNumber z
        z = self._new()
        mpc_tan (z.value, self.value, (<MPComplexField_class>self._parent).__rnd)
        return z

    def cosh(self):
        r"""
        Return the hyperbolic cosine of this complex number:

        .. MATH::

            \cosh(a + ib) = \cosh a \cos b + i \sinh a \sin b.

        EXAMPLES::

            sage: MPC = MPComplexField()
            sage: u = MPC(2, 4)
            sage: cosh(u)
            -2.45913521391738 - 2.74481700679215*I
        """
        cdef MPComplexNumber z
        z = self._new()
        mpc_cosh (z.value, self.value, (<MPComplexField_class>self._parent).__rnd)
        return z

    def sinh(self):
        r"""
        Return the hyperbolic sine of this complex number:

        .. MATH::

            \sinh(a + ib) = \sinh a \cos b + i \cosh a \sin b.

        EXAMPLES::

            sage: MPC = MPComplexField()
            sage: u = MPC(2, 4)
            sage: sinh(u)
            -2.37067416935200 - 2.84723908684883*I
        """
        cdef MPComplexNumber z
        z = self._new()
        mpc_sinh (z.value, self.value, (<MPComplexField_class>self._parent).__rnd)
        return z

    def tanh(self):
        r"""
        Return the hyperbolic tangent of this complex number:

        .. MATH::

            \tanh(a + ib) = (\sinh 2a + i \sin 2b)/(\cosh 2a + \cos 2b).

        EXAMPLES::

            sage: MPC = MPComplexField()
            sage: u = MPC(2, 4)
            sage: tanh(u)
            1.00468231219024 + 0.0364233692474037*I
        """
        cdef MPComplexNumber z
        z = self._new()
        mpc_tanh (z.value, self.value, (<MPComplexField_class>self._parent).__rnd)
        return z

    def arccos(self):
        """
        Return the arccosine of this complex number.

        EXAMPLES::

            sage: MPC = MPComplexField()
            sage: u = MPC(2, 4)
            sage: arccos(u)
            1.11692611683177 - 2.19857302792094*I
        """
        cdef MPComplexNumber z
        z = self._new()
        mpc_acos (z.value, self.value, (<MPComplexField_class>self._parent).__rnd)
        return z

    def arcsin(self):
        """
        Return the arcsine of this complex number.

        EXAMPLES::

            sage: MPC = MPComplexField()
            sage: u = MPC(2, 4)
            sage: arcsin(u)
            0.453870209963122 + 2.19857302792094*I
        """
        cdef MPComplexNumber z
        z = self._new()
        mpc_asin (z.value, self.value, (<MPComplexField_class>self._parent).__rnd)
        return z

    def arctan(self):
        """
        Return the arctangent of this complex number.

        EXAMPLES::

            sage: MPC = MPComplexField()
            sage: u = MPC(-2, 4)
            sage: arctan(u)
            -1.46704821357730 + 0.200586618131234*I
        """
        cdef MPComplexNumber z
        z = self._new()
        mpc_atan (z.value, self.value, (<MPComplexField_class>self._parent).__rnd)
        return z

    def arccosh(self):
        """
        Return the hyperbolic arccos of this complex number.

        EXAMPLES::

            sage: MPC = MPComplexField()
            sage: u = MPC(2, 4)
            sage: arccosh(u)
            2.19857302792094 + 1.11692611683177*I
        """
        cdef MPComplexNumber z
        z = self._new()
        mpc_acosh (z.value, self.value, (<MPComplexField_class>self._parent).__rnd)
        return z

    def arcsinh(self):
        """
        Return the hyperbolic arcsine of this complex number.

        EXAMPLES::

            sage: MPC = MPComplexField()
            sage: u = MPC(2, 4)
            sage: arcsinh(u)
            2.18358521656456 + 1.09692154883014*I
        """
        cdef MPComplexNumber z
        z = self._new()
        mpc_asinh (z.value, self.value, (<MPComplexField_class>self._parent).__rnd)
        return z

    def arctanh(self):
        """
        Return the hyperbolic arctangent of this complex number.

        EXAMPLES::

            sage: MPC = MPComplexField()
            sage: u = MPC(2, 4)
            sage: arctanh(u)
            0.0964156202029962 + 1.37153510396169*I
        """
        cdef MPComplexNumber z
        z = self._new()
        mpc_atanh (z.value, self.value, (<MPComplexField_class>self._parent).__rnd)
        return z

    def coth(self):
        """
        Return the hyperbolic cotangent of this complex number.

        EXAMPLES::

            sage: MPC = MPComplexField(100)
            sage: MPC(1,1).coth()
            0.86801414289592494863584920892 - 0.21762156185440268136513424361*I
        """
        return ~(self.tanh())

    def arccoth(self):
        """
        Return the hyperbolic arccotangent of this complex number.

        EXAMPLES::

            sage: MPC = MPComplexField(100)
            sage: MPC(1,1).arccoth()
            0.40235947810852509365018983331 - 0.55357435889704525150853273009*I
        """
        return (~self).arctanh()

    def csc(self):
        """
        Return the cosecant of this complex number.

        EXAMPLES::

            sage: MPC = MPComplexField(100)
            sage: MPC(1,1).csc()
            0.62151801717042842123490780586 - 0.30393100162842645033448560451*I
        """
        return ~(self.sin())

    def csch(self):
        """
        Return the hyperbolic cosecant of this complex number.

        EXAMPLES::

            sage: MPC = MPComplexField(100)
            sage: MPC(1,1).csch()
            0.30393100162842645033448560451 - 0.62151801717042842123490780586*I
        """
        return ~(self.sinh())

    def arccsch(self):
        """
        Return the hyperbolic arcsine of this complex number.

        EXAMPLES::

            sage: MPC = MPComplexField(100)
            sage: MPC(1,1).arccsch()
            0.53063753095251782601650945811 - 0.45227844715119068206365839783*I
        """
        return (~self).arcsinh()

    def sec(self):
        """
        Return the secant of this complex number.

        EXAMPLES::

            sage: MPC = MPComplexField(100)
            sage: MPC(1,1).sec()
            0.49833703055518678521380589177 + 0.59108384172104504805039169297*I
        """
        return ~(self.cos())

    def sech(self):
        """
        Return the hyperbolic secant of this complex number.

        EXAMPLES::

            sage: MPC = MPComplexField(100)
            sage: MPC(1,1).sech()
            0.49833703055518678521380589177 - 0.59108384172104504805039169297*I
        """
        return ~(self.cosh())

    def arcsech(self):
        """
        Return the hyperbolic arcsecant of this complex number.

        EXAMPLES::

            sage: MPC = MPComplexField(100)
            sage: MPC(1,1).arcsech()
            0.53063753095251782601650945811 - 1.1185178796437059371676632938*I
        """
        return (~self).arccosh()

    def cot(self):
        """
        Return the cotangent of this complex number.

        EXAMPLES::

            sage: MPC = MPComplexField(53)
            sage: (1+MPC(I)).cot()
            0.217621561854403 - 0.868014142895925*I
            sage: i = MPComplexField(200).0
            sage: (1+i).cot()
            0.21762156185440268136513424360523807352075436916785404091068 - 0.86801414289592494863584920891627388827343874994609327121115*I
            sage: i = MPComplexField(220).0
            sage: (1+i).cot()
            0.21762156185440268136513424360523807352075436916785404091068124239 - 0.86801414289592494863584920891627388827343874994609327121115071646*I
        """
        return ~(self.tan())

    ################################
    # Other functions
    ################################

    def argument(self):
        r"""
        The argument (angle) of the complex number, normalized so that
        `-\pi < \theta \leq \pi`.

        EXAMPLES::

            sage: MPC = MPComplexField()
            sage: i = MPC.0
            sage: (i^2).argument()
            3.14159265358979
            sage: (1+i).argument()
            0.785398163397448
            sage: i.argument()
            1.57079632679490
            sage: (-i).argument()
            -1.57079632679490
            sage: (RR('-0.001') - i).argument()
            -1.57179632646156
        """
        cdef RealNumber x
        x = RealNumber(self._parent._real_field())
        mpc_arg(x.value, self.value, (<RealField_class>x._parent).rnd)
        return x

    def conjugate(self):
        r"""
        Return the complex conjugate of this complex number:

        .. MATH::

            \mathrm{conjugate}(a + ib) = a - ib.

        EXAMPLES::

            sage: MPC = MPComplexField()
            sage: i = MPC(0, 1)
            sage: (1+i).conjugate()
            1.00000000000000 - 1.00000000000000*I
        """
        cdef MPComplexNumber z
        z = self._new()
        mpc_conj(z.value, self.value, (<MPComplexField_class>self._parent).__rnd)
        return z

    def sqr(self):
        """
        Return the square of a complex number:

        .. MATH::

            (a + ib)^2 = (a^2 - b^2) + 2iab.

        EXAMPLES::

            sage: C = MPComplexField()
            sage: a = C(5, 1)
            sage: a.sqr()
            24.0000000000000 + 10.0000000000000*I
        """
        cdef MPComplexNumber z
        z = self._new()
        mpc_sqr (z.value, self.value, (<MPComplexField_class>self._parent).__rnd)
        return z

    def sqrt(self):
        r"""
        Return the square root, taking the branch cut to be the negative real
        axis:

        .. MATH::

            \sqrt z = \sqrt{|z|}(\cos(\arg(z)/2) + i \sin(\arg(z)/2)).

        EXAMPLES::

            sage: C = MPComplexField()
            sage: a = C(24, 10)
            sage: a.sqrt()
            5.00000000000000 + 1.00000000000000*I
        """
        cdef MPComplexNumber z
        z = self._new()
        mpc_sqrt(z.value, self.value, (<MPComplexField_class>self._parent).__rnd)
        return z

    def exp(self):
        r"""
        Return the exponential of this complex number:

        .. MATH::

            \exp(a + ib) = \exp(a) (\cos b + i \sin b).

        EXAMPLES::

            sage: MPC = MPComplexField()
            sage: u = MPC(2, 4)
            sage: exp(u)
            -4.82980938326939 - 5.59205609364098*I
        """
        cdef MPComplexNumber z
        z = self._new()
        mpc_exp(z.value, self.value, (<MPComplexField_class>self._parent).__rnd)
        return z

    def log(self):
        r"""
        Return the logarithm of this complex number with the branch
        cut on the negative real axis:

        .. MATH::

            \log(z) = \log |z| + i \arg(z).

        EXAMPLES::

            sage: MPC = MPComplexField()
            sage: u = MPC(2, 4)
            sage: log(u)
            1.49786613677700 + 1.10714871779409*I
        """
        cdef MPComplexNumber z
        z = self._new()
        mpc_log(z.value, self.value, (<MPComplexField_class>self._parent).__rnd)
        return z

    def __lshift__(self, n):
        """
        Fast multiplication by `2^n`.

        EXAMPLES::

            sage: MPC = MPComplexField()
            sage: u = MPC(2, 4)
            sage: u<<2
            8.00000000000000 + 16.0000000000000*I
            sage: u<<(-1)
            1.00000000000000 + 2.00000000000000*I
        """
        cdef MPComplexNumber z, x
        x = <MPComplexNumber>self
        z = x._new()
        mpc_mul_2si(z.value, x.value, n, (<MPComplexField_class>x._parent).__rnd)
        return z

    def __rshift__(self, n):
        """
        Fast division by `2^n`.

        EXAMPLES::

            sage: MPC = MPComplexField()
            sage: u = MPC(2, 4)
            sage: u>>2
            0.500000000000000 + 1.00000000000000*I
            sage: u>>(-1)
            4.00000000000000 + 8.00000000000000*I
        """
        cdef MPComplexNumber z, x
        x = <MPComplexNumber>self
        z = x._new()
        mpc_div_2si(z.value, x.value, n, (<MPComplexField_class>x._parent).__rnd)
        return z

    def nth_root(self, n, all=False):
        """
        The `n`-th root function.

        INPUT:

        - ``all`` -- boolean (default: ``False``); if ``True``, return a
          list of all `n`-th roots

        EXAMPLES::

            sage: MPC = MPComplexField()
            sage: a = MPC(27)
            sage: a.nth_root(3)
            3.00000000000000
            sage: a.nth_root(3, all=True)
            [3.00000000000000, -1.50000000000000 + 2.59807621135332*I, -1.50000000000000 - 2.59807621135332*I]
        """
        if self.is_zero():
            return [self] if all else self

        cdef RealField_class R = self._parent._real_field()
        cdef mpfr_rnd_t rrnd = R.rnd
        cdef mpc_rnd_t crnd = (<MPComplexField_class>(self._parent)).__rnd

        cdef RealNumber a,r
        a = self.argument()/n
        r = self.abs()
        mpfr_rootn_ui(r.value, r.value, n, rrnd)

        cdef MPComplexNumber z
        z = self._new()
        mpfr_sin_cos(z.value.im, z.value.re, a.value, rrnd)
        mpc_mul_fr(z.value, z.value, r.value, crnd)

        if not all:
            return z

        cdef RealNumber t, tt
        t = R.pi()*2/n
        tt= RealNumber(R)

        cdef list zlist = [z]
        for i in range(1,n):
            z = self._new()
            mpfr_mul_ui(tt.value, t.value, i, rrnd)
            mpfr_add(tt.value, tt.value, a.value, rrnd)
            mpfr_sin_cos(z.value.im, z.value.re, tt.value, rrnd)
            mpc_mul_fr(z.value, z.value, r.value, crnd)
            zlist.append(z)

        return zlist

    def dilog(self):
        r"""
        Return the complex dilogarithm of ``self``.

        The complex dilogarithm, or Spence's function, is defined by

        .. MATH::

            Li_2(z) = - \int_0^z \frac{\log|1-\zeta|}{\zeta} d(\zeta)
            = \sum_{k=1}^\infty \frac{z^k}{k^2}.

        Note that the series definition can only be used for `|z| < 1`.

        EXAMPLES::

            sage: MPC = MPComplexField()
            sage: a = MPC(1,0)
            sage: a.dilog()
            1.64493406684823
            sage: float(pi^2/6)                                                         # needs sage.symbolic
            1.6449340668482262

        ::

            sage: b = MPC(0,1)
            sage: b.dilog()
            -0.205616758356028 + 0.915965594177219*I

        ::

            sage: c = MPC(0,0)
            sage: c.dilog()
            0
        """
        return self._parent(self.__pari__().dilog())

    def eta(self, omit_frac=False):
        r"""
        Return the value of the Dedekind `\eta` function on ``self``,
        intelligently computed using `\mathbb{SL}(2,\ZZ)` transformations.

        The `\eta` function is

        .. MATH::

            \eta(z) = e^{\pi i z / 12} \prod_{n=1}^{\infty}(1-e^{2\pi inz})

        INPUT:

        - ``self`` -- element of the upper half plane (if not,
          raises a :exc:`ValueError`)

        - ``omit_frac`` -- -- boolean (default: ``False``); if ``True``,
          omit the `e^{\pi i z / 12}` factor

        OUTPUT: a complex number

        ALGORITHM: Uses the PARI C library.

        EXAMPLES::

            sage: MPC = MPComplexField()
            sage: i = MPC.0
            sage: z = 1+i; z.eta()
            0.742048775836565 + 0.198831370229911*I
        """
        try:
            return self._parent(self.__pari__().eta(not omit_frac))
        except PariError:
            raise ValueError("value must be in the upper half plane")

    def gamma(self):
        """
        Return the Gamma function evaluated at this complex number.

        EXAMPLES::

            sage: MPC = MPComplexField(30)
            sage: i = MPC.0
            sage: (1+i).gamma()
            0.49801567 - 0.15494983*I

        TESTS::

            sage: MPC(0).gamma()
            Infinity

        ::

            sage: MPC(-1).gamma()
            Infinity
        """
        try:
            return self._parent(self.__pari__().gamma())
        except PariError:
            from sage.rings.infinity import UnsignedInfinityRing
            return UnsignedInfinityRing.gen()

    def gamma_inc(self, t):
        """
        Return the incomplete Gamma function evaluated at this complex number.

        EXAMPLES::

            sage: C, i = MPComplexField(30).objgen()
            sage: (1+i).gamma_inc(2 + 3*i)  # abs tol 2e-10
            0.0020969149 - 0.059981914*I
            sage: (1+i).gamma_inc(5)
            -0.0013781309 + 0.0065198200*I
            sage: C(2).gamma_inc(1 + i)
            0.70709210 - 0.42035364*I
        """
        return self._parent(self.__pari__().incgam(t, precision=self.prec()))

    def zeta(self):
        """
        Return the Riemann zeta function evaluated at this complex number.

        EXAMPLES::

            sage: i = MPComplexField(30).gen()
            sage: z = 1 + i
            sage: z.zeta()
            0.58215806 - 0.92684856*I
        """
        return self._parent(self.__pari__().zeta())

    def agm(self, right, algorithm='optimal'):
        """
        Return the algebro-geometric mean of ``self`` and ``right``.

        EXAMPLES::

            sage: MPC = MPComplexField()
            sage: u = MPC(1, 4)
            sage: v = MPC(-2,5)
            sage: u.agm(v, algorithm='pari')
            -0.410522769709397 + 4.60061063922097*I
            sage: u.agm(v, algorithm='principal')
            1.24010691168158 - 0.472193567796433*I
            sage: u.agm(v, algorithm='optimal')
            -0.410522769709397 + 4.60061063922097*I
        """
        if algorithm == "pari":
            t = self._parent(right).__pari__()
            return self._parent(self.__pari__().agm(t))

        cdef MPComplexNumber a, b, d, s, res
        cdef mp_exp_t rel_prec
        cdef bint optimal = algorithm == "optimal"

        cdef mpc_rnd_t rnd = (<MPComplexField_class>(self._parent)).__rnd

        cdef int prec = self._parent._prec

        if optimal or algorithm == "principal":
            if not isinstance(right, MPComplexNumber) or (<MPComplexNumber>right)._parent is not self._parent:
                right = self._parent(right)

            res = self._new()

            if self.is_zero():
                return self
            elif (<MPComplexNumber>right).is_zero():
                return right
            elif (mpfr_cmpabs(self.value.re, (<MPComplexNumber>right).value.re) == 0 and
                  mpfr_cmpabs(self.value.im, (<MPComplexNumber>right).value.im) == 0 and
                  mpfr_cmp(self.value.re, (<MPComplexNumber>right).value.re) != 0 and
                  mpfr_cmp(self.value.im, (<MPComplexNumber>right).value.im) != 0):
                # self = -right
                mpc_set_ui(res.value, 0, rnd)
                return res
            # Do the computations to a bit higher precision so rounding error
            # won't obscure the termination condition.
            a = MPComplexNumber(MPComplexField(prec+5), None)
            b = a._new()
            d = a._new()
            s = a._new()

            # Make copies so we don't mutate self or right.
            mpc_set(a.value, self.value, rnd)
            mpc_set(b.value, (<MPComplexNumber>right).value, rnd)

            if optimal:
                mpc_add(s.value, a.value, b.value, rnd)

            while True:
                # s = a+b
                if not optimal:
                    mpc_add(s.value, a.value, b.value, rnd)

                # b = sqrt(a*b)
                mpc_mul(d.value, a.value, b.value, rnd)
                mpc_sqrt(b.value, d.value, rnd)

                # a = s/2
                mpc_div_2ui(a.value, s.value, 1, rnd)

                # d = a - b
                mpc_sub(d.value, a.value, b.value, rnd)
                if d.is_zero():
                    mpc_set(res.value, a.value, rnd)
                    return res

                if optimal:
                    # s = a+b
                    mpc_add(s.value, a.value, b.value, rnd)
                    if s.is_zero():
                        mpc_set(res.value, a.value, rnd)
                        return res

                    # |s| < |d|
                    if s.norm() < d.norm():
                        mpc_swap(d.value, s.value)
                        mpc_neg(b.value, b.value, rnd)

                rel_prec = min_exp_t(max_exp(a), max_exp(b)) - max_exp(d)
                if rel_prec > prec:
                    mpc_set(res.value, a.value, rnd)
                    return res

cdef inline mp_exp_t min_exp_t(mp_exp_t a, mp_exp_t b) noexcept:
    return a if a < b else b
cdef inline mp_exp_t max_exp_t(mp_exp_t a, mp_exp_t b) noexcept:
    return a if a > b else b

cdef inline mp_exp_t max_exp(MPComplexNumber z) noexcept:
    """
    Quickly return the maximum exponent of the real and complex parts of ``z``,
    which is useful for estimating its magnitude.
    """
    if mpfr_zero_p(z.value.im):
        return mpfr_get_exp(z.value.re)
    elif mpfr_zero_p(z.value.re):
        return mpfr_get_exp(z.value.im)
    return max_exp_t(mpfr_get_exp(z.value.re), mpfr_get_exp(z.value.im))


def __create__MPComplexField_version0 (prec, rnd):
    """
    Create a :class:`MPComplexField`.

    EXAMPLES::

        sage: sage.rings.complex_mpc.__create__MPComplexField_version0(200, 'RNDDZ')
        Complex Field with 200 bits of precision and rounding RNDDZ
    """
    return MPComplexField(prec, rnd)


def __create__MPComplexNumber_version0 (parent, s, base=10):
    """
    Create a :class:`MPComplexNumber`.

    EXAMPLES::

        sage: C = MPComplexField(prec=20, rnd='RNDUU')
        sage: sage.rings.complex_mpc.__create__MPComplexNumber_version0(C, 3.2+2*i)
        3.2001 + 2.0000*I
        sage: C = MPComplexField(prec=20, rnd='RNDUU')
        sage: sage.rings.complex_mpc.__create__MPComplexNumber_version0(C,"33.",10)
        33.000
    """
    return MPComplexNumber(parent, s, base=base)


# original version of the file had this with only 1 underscore - TCS
__create_MPComplexNumber_version0 = __create__MPComplexNumber_version0

#*****************************************************************************
#
#       Morphisms
#
#*****************************************************************************

cdef class MPCtoMPC(Map):
    cpdef Element _call_(self, z):
        """
        EXAMPLES::

            sage: from sage.rings.complex_mpc import *
            sage: C10 = MPComplexField(10)
            sage: C100 = MPComplexField(100)
            sage: f = MPCtoMPC(C100, C10) # indirect doctest
            sage: a = C100(1.2, 24)
            sage: f(a)
            1.2 + 24.*I
            sage: f
            Generic map:
              From: Complex Field with 100 bits of precision
              To:   Complex Field with 10 bits of precision
        """
        cdef MPComplexNumber y
        y = (<MPComplexField_class>self.codomain())._new()
        y._set(z)
        return y

    def section(self):
        """
        EXAMPLES::

            sage: from sage.rings.complex_mpc import *
            sage: C10 = MPComplexField(10)
            sage: C100 = MPComplexField(100)
            sage: f = MPCtoMPC(C100, C10)
            sage: f.section()
            Generic map:
              From: Complex Field with 10 bits of precision
              To:   Complex Field with 100 bits of precision
        """
        return MPCtoMPC(self.codomain(), self.domain())

cdef class INTEGERtoMPC(Map):
    cpdef Element _call_(self, x):
        """
        EXAMPLES::

            sage: from sage.rings.complex_mpc import *
            sage: I = IntegerRing()
            sage: C100 = MPComplexField(100)
            sage: f = MPFRtoMPC(I, C100); f # indirect doctest
            Generic map:
              From: Integer Ring
              To:   Complex Field with 100 bits of precision
            sage: a = I(625)
            sage: f(a)
            625.00000000000000000000000000
        """
        cdef MPComplexNumber y
        cdef mpc_rnd_t rnd
        rnd =(<MPComplexField_class>self._parent).__rnd
        y = (<MPComplexField_class>self.codomain())._new()
        mpc_set_z(y.value, (<Integer>x).value, rnd)
        return y

cdef class MPFRtoMPC(Map):
    cpdef Element _call_(self, x):
        """
        EXAMPLES::

            sage: from sage.rings.complex_mpc import *
            sage: R10 = RealField(10)
            sage: C100 = MPComplexField(100)
            sage: f = MPFRtoMPC(R10, C100); f # indirect doctest
            Generic map:
              From: Real Field with 10 bits of precision
              To:   Complex Field with 100 bits of precision
            sage: a = R10(1.625)
            sage: f(a)
            1.6250000000000000000000000000
        """
        cdef MPComplexNumber y
#        cdef mpc_rnd_t rnd
#        rnd =(<MPComplexField_class>self._parent).__rnd
        y = (<MPComplexField_class>self.codomain())._new()
#        mpc_set_fr(y.value, (<RealNumber>x).value, rnd)
        y._set(x)
        return y

cdef class CCtoMPC(Map):
    cpdef Element _call_(self, z):
        """
        EXAMPLES::

            sage: from sage.rings.complex_mpc import *
            sage: C10 = ComplexField(10)
            sage: MPC100 = MPComplexField(100)
            sage: f = CCtoMPC(C10, MPC100); f # indirect doctest
            Generic map:
              From: Complex Field with 10 bits of precision
              To:   Complex Field with 100 bits of precision
            sage: a = C10(1.625, 42)
            sage: f(a)
            1.6250000000000000000000000000 + 42.000000000000000000000000000*I
        """
        cdef MPComplexNumber y
        cdef mpc_rnd_t rnd
        rnd =(<MPComplexField_class>self._parent).__rnd
        y = (<MPComplexField_class>self.codomain())._new()
        mpc_set_fr_fr(y.value, (<ComplexNumber>z).__re, (<ComplexNumber>z).__im, rnd)
        return y


# Support Python's numbers abstract base class
# import numbers
