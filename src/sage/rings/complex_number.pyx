"""
Arbitrary Precision Complex Numbers

AUTHORS:

- William Stein (2006-01-26): complete rewrite

- Joel B. Mohler (2006-12-16): naive rewrite into pyrex

- William Stein(2007-01): rewrite of Mohler's rewrite

- Vincent Delecroix (2010-01): plot function

- Travis Scrimshaw (2012-10-18): Added documentation for full coverage
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function

import math
import operator

from sage.structure.element cimport FieldElement, RingElement, Element, ModuleElement
from sage.categories.map cimport Map

from complex_double cimport ComplexDoubleElement
from real_mpfr cimport RealNumber

import sage.misc.misc
import integer
import infinity

from sage.libs.mpmath.utils cimport mpfr_to_mpfval
from sage.rings.integer_ring import ZZ


cdef object numpy_complex_interface = {'typestr': '=c16'}
cdef object numpy_object_interface = {'typestr': '|O'}

cdef mpfr_rnd_t rnd
rnd = MPFR_RNDN

cdef double LOG_TEN_TWO_PLUS_EPSILON = 3.321928094887363 # a small overestimate of log(10,2)

def set_global_complex_round_mode(n):
    """
    Set the global complex rounding mode.

    .. WARNING::

        Do not call this function explicitly. The default rounding mode is
        ``n = 0``.

    EXAMPLES::

        sage: sage.rings.complex_number.set_global_complex_round_mode(0)
    """
    global rnd
    rnd = n

#from sage.databases.odlyzko import zeta_zeroes

def is_ComplexNumber(x):
    r"""
    Returns ``True`` if ``x`` is a complex number. In particular, if ``x`` is
    of the :class:`ComplexNumber` type.

    EXAMPLES::

        sage: from sage.rings.complex_number import is_ComplexNumber
        sage: a = ComplexNumber(1,2); a
        1.00000000000000 + 2.00000000000000*I
        sage: is_ComplexNumber(a)
        True
        sage: b = ComplexNumber(1); b
        1.00000000000000
        sage: is_ComplexNumber(b)
        True

    Note that the global element ``I`` is of type :class:`SymbolicConstant`.
    However, elements of the class :class:`ComplexField_class` are of type
    :class:`ComplexNumber`::

        sage: c = 1 + 2*I
        sage: is_ComplexNumber(c)
        False
        sage: d = CC(1 + 2*I)
        sage: is_ComplexNumber(d)
        True
    """
    return isinstance(x, ComplexNumber)

cdef class ComplexNumber(sage.structure.element.FieldElement):
    """
    A floating point approximation to a complex number using any
    specified precision. Answers derived from calculations with such
    approximations may differ from what they would be if those
    calculations were performed with true complex numbers. This is due
    to the rounding errors inherent to finite precision calculations.

    EXAMPLES::

        sage: I = CC.0
        sage: b = 1.5 + 2.5*I
        sage: loads(b.dumps()) == b
        True
    """

    cdef ComplexNumber _new(self):
        """
        Quickly creates a new initialized complex number with the same
        parent as ``self``.
        """
        cdef ComplexNumber x
        x = ComplexNumber.__new__(ComplexNumber)
        x._parent = self._parent
        x._prec = self._prec
        x._multiplicative_order = None
        mpfr_init2(x.__re, self._prec)
        mpfr_init2(x.__im, self._prec)
        return x

    def __cinit__(self, parent=None, real=None, imag=None):
        """
        Cython initialize ``self``.

        EXAMPLES::

            sage: ComplexNumber(2,1) # indirect doctest
            2.00000000000000 + 1.00000000000000*I
        """
        self._prec = -1

    def __init__(self, parent, real, imag=None):
        r"""
        Initialize :class:`ComplexNumber` instance.

        EXAMPLES::

            sage: a = ComplexNumber(2,1)
            sage: a.__init__(CC,2,1)
            sage: a
            2.00000000000000 + 1.00000000000000*I
            sage: parent(a)
            Complex Field with 53 bits of precision
            sage: real(a)
            2.00000000000000
            sage: imag(a)
            1.00000000000000
        """
        cdef RealNumber rr, ii
        self._parent = parent
        self._prec = self._parent._prec
        self._multiplicative_order = None

        mpfr_init2(self.__re, self._prec)
        mpfr_init2(self.__im, self._prec)

        if imag is None:
            if real is None: return

            if isinstance(real, ComplexNumber):
                real, imag = (<ComplexNumber>real).real(), (<ComplexNumber>real).imag()
            elif isinstance(real, sage.libs.pari.all.pari_gen):
                real, imag = real.real(), real.imag()
            elif isinstance(real, list) or isinstance(real, tuple):
                re, imag = real
                real = re
            elif isinstance(real, complex):
                real, imag = real.real, real.imag
            else:
                imag = 0
        try:
            R = parent._real_field()
            rr = R(real)
            ii = R(imag)
            mpfr_set(self.__re, rr.value, rnd)
            mpfr_set(self.__im, ii.value, rnd)
        except TypeError:
            raise TypeError("unable to coerce to a ComplexNumber: %s" % type(real))


    def  __dealloc__(self):
        """
        TESTS:

        Check that :trac:`12038` is resolved::

            sage: from sage.rings.complex_number import ComplexNumber as CN
            sage: coerce(CN, 1+I)
            Traceback (most recent call last):
            ...
            TypeError: __init__() takes at least 2 positional arguments (1 given)
        """
        if self._prec != -1:
            mpfr_clear(self.__re)
            mpfr_clear(self.__im)

    def _interface_init_(self, I=None):
        """
        Returns ``self`` formatted as a string, suitable as input to another
        computer algebra system. (This is the default function used for
        exporting to other computer algebra systems.)

        EXAMPLES::

            sage: s1 = CC(exp(I)); s1
            0.540302305868140 + 0.841470984807897*I
            sage: s1._interface_init_()
            '0.54030230586813977 + 0.84147098480789650*I'
            sage: s1 == CC(gp(s1))
            True
        """
        return self.str(truncate=False)

    def _maxima_init_(self, I=None):
        """
        Return a string representation of this complex number in the syntax of
        Maxima. That is, use ``%i`` to represent the complex unit.

        EXAMPLES::

            sage: CC.0._maxima_init_()
            '1.0000000000000000*%i'
            sage: CC(.5 + I)._maxima_init_()
            '0.50000000000000000 + 1.0000000000000000*%i'
        """
        return self.str(truncate=False, istr='%i')

    @property
    def __array_interface__(self):
        """
        Used for NumPy conversion.

        EXAMPLES::

            sage: import numpy
            sage: numpy.array([1.0, 2.5j]).dtype
            dtype('complex128')
            sage: numpy.array([1.000000000000000000000000000000000000j]).dtype
            dtype('O')
        """
        if self._prec <= 53:
            return numpy_complex_interface
        else:
            return numpy_object_interface

    def _sage_input_(self, sib, coerced):
        r"""
        Produce an expression which will reproduce this value when evaluated.

        EXAMPLES::

            sage: for prec in (2, 53, 200):
            ...       fld = ComplexField(prec)
            ...       var = polygen(fld)
            ...       ins = [-20, 0, 1, -2^4000, 2^-4000] + [fld._real_field().random_element() for _ in range(3)]
            ...       for v1 in ins:
            ...           for v2 in ins:
            ...               v = fld(v1, v2)
            ...               _ = sage_input(fld(v), verify=True)
            ...               _ = sage_input(fld(v) * var, verify=True)
            sage: x = polygen(CC)
            sage: for v1 in [-2, 0, 2]:
            ....:     for v2 in [-2, -1, 0, 1, 2]:
            ....:         print(str(sage_input(x + CC(v1, v2))).splitlines()[1])
            x + CC(-2 - RR(2)*I)
            x + CC(-2 - RR(1)*I)
            x - 2
            x + CC(-2 + RR(1)*I)
            x + CC(-2 + RR(2)*I)
            x - CC(RR(2)*I)
            x - CC(RR(1)*I)
            x
            x + CC(RR(1)*I)
            x + CC(RR(2)*I)
            x + CC(2 - RR(2)*I)
            x + CC(2 - RR(1)*I)
            x + 2
            x + CC(2 + RR(1)*I)
            x + CC(2 + RR(2)*I)
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: sib = SageInputBuilder()
            sage: sib_np = SageInputBuilder(preparse=False)
            sage: CC(-infinity)._sage_input_(sib, True)
            {unop:- {call: {atomic:RR}({atomic:Infinity})}}
            sage: CC(0, infinity)._sage_input_(sib, True)
            {call: {atomic:CC}({call: {atomic:RR}({atomic:0})}, {call: {atomic:RR}({atomic:Infinity})})}
            sage: CC(NaN, 5)._sage_input_(sib, True)
            {call: {atomic:CC}({call: {atomic:RR}({atomic:NaN})}, {call: {atomic:RR}({atomic:5})})}
            sage: CC(5, NaN)._sage_input_(sib, True)
            {call: {atomic:CC}({call: {atomic:RR}({atomic:5})}, {call: {atomic:RR}({atomic:NaN})})}
            sage: CC(12345)._sage_input_(sib, True)
            {atomic:12345}
            sage: CC(-12345)._sage_input_(sib, False)
            {call: {atomic:CC}({binop:+ {unop:- {atomic:12345}} {binop:* {call: {atomic:RR}({atomic:0})} {atomic:I}}})}
            sage: CC(0, 12345)._sage_input_(sib, True)
            {call: {atomic:CC}({binop:* {call: {atomic:RR}({atomic:12345})} {atomic:I}})}
            sage: CC(0, -12345)._sage_input_(sib, False)
            {unop:- {call: {atomic:CC}({binop:* {call: {atomic:RR}({atomic:12345})} {atomic:I}})}}
            sage: CC(1.579)._sage_input_(sib, True)
            {atomic:1.579}
            sage: CC(1.579)._sage_input_(sib_np, True)
            {atomic:1.579}
            sage: ComplexField(150).zeta(37)._sage_input_(sib, True)
            {call: {call: {atomic:ComplexField}({atomic:150})}({binop:+ {atomic:0.98561591034770846226477029397621845736859851519} {binop:* {call: {call: {atomic:RealField}({atomic:150})}({atomic:0.16900082032184907409303555538443060626072476297})} {atomic:I}}})}
            sage: ComplexField(150).zeta(37)._sage_input_(sib_np, True)
            {call: {call: {atomic:ComplexField}({atomic:150})}({binop:+ {call: {call: {atomic:RealField}({atomic:150})}({atomic:'0.98561591034770846226477029397621845736859851519'})} {binop:* {call: {call: {atomic:RealField}({atomic:150})}({atomic:'0.16900082032184907409303555538443060626072476297'})} {atomic:I}}})}
        """
        if coerced and self.imag() == 0:
            return sib(self.real(), True)

        # The body will be coerced first to symbolics and then to CC.
        # This works fine if we produce integer or float literals, but
        # not for infinity or NaN.
        if not (mpfr_number_p(self.__re) and mpfr_number_p(self.__im)):
            return sib(self.parent())(self.real(), self.imag())

        # The following uses of .sum() and .prod() will simplify
        # 3+0*I to 3, 0+1*I to I, etc.
        real_part = sib(self.real(), 2)
        imag_part = sib.prod([sib(self.imag()), sib.name('I')],
                             simplify=True)
        sum = sib.sum([real_part, imag_part], simplify=True)
        if sum._sie_is_negation():
            return -sib(self.parent())(sum._sie_operand)
        else:
            return sib(self.parent())(sum)

        # The following is an (untested) implementation that produces
        # CC_I = CC.gen()
        # 2 + 3*CC_I
        # instead of CC(2 + 3*I)
#         cdef int prec

#         if self.real().is_zero() and self.imag() == 1:
#             v = sib(self.parent()).gen()
#             prec = self.prec()
#             if prec == 53:
#                 gen_name = 'CC_I'
#             else:
#                 gen_name = 'CC%d_I' % prec
#             sib.cache(self, v, gen_name)

#         real_part = sib(self.real())
#         imag_part = sib.prod([self.imag(), self.parent().gen()], simplify=True)
#         return sib.sum([real_part, imag_part], simplify=True)

    def _repr_(self):
        r"""
        Returns ``self`` formatted as a string.

        EXAMPLES::

            sage: a = ComplexNumber(2,1); a
            2.00000000000000 + 1.00000000000000*I
            sage: a._repr_()
            '2.00000000000000 + 1.00000000000000*I'
        """
        return self.str(10)

    def __hash__(self):
        """
        Returns the hash of ``self``, which coincides with the python complex
        and float (and often int) types.

        This has the drawback that two very close high precision numbers
        will have the same hash, but allows them to play nicely with other
        real types.

        EXAMPLES::

            sage: hash(CC(1.2, 33)) == hash(complex(1.2, 33))
            True
        """
        return hash(complex(self))

    def __getitem__(self, i):
        r"""
        Returns either the real or imaginary component of self depending on
        the choice of i: real (i=0), imaginary (i=1)

        INPUT:

        - ``i`` - 0 or 1
            - ``0`` -- will return the real component of ``self``
            - ``1`` -- will return the imaginary component of ``self``

        EXAMPLES::

            sage: a = ComplexNumber(2,1)
            sage: a.__getitem__(0)
            2.00000000000000
            sage: a.__getitem__(1)
            1.00000000000000

        ::

            sage: b = CC(42,0)
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

    def __reduce__( self ):
        """
        Pickling support

        EXAMPLES::

            sage: a = CC(1 + I)
            sage: loads(dumps(a)) == a
            True
        """
        # TODO: This is potentially slow -- make a 1 version that
        # is native and much faster -- doesn't use .real()/.imag()
        return (make_ComplexNumber0, (self._parent, self._multiplicative_order, self.real(), self.imag()))

    def _set_multiplicative_order(self, n):
        r"""
        Function for setting the attribute :meth:`multiplicative_order` of
        ``self``.

        .. WARNING::

            It is not advisable to explicitly call
            ``_set_multiplicative_order()`` for explicitly declared complex
            numbers.

        INPUT:

        - ``n`` -- an integer which will define the multiplicative order

        EXAMPLES::

            sage: a = ComplexNumber(2,1)
            sage: a.multiplicative_order()
            +Infinity
            sage: a._set_multiplicative_order(5)
            sage: a.multiplicative_order()
            5
            sage: a^5
            -38.0000000000000 + 41.0000000000000*I
        """
        self._multiplicative_order = integer.Integer(n)

    def str(self, base=10, truncate=True, istr='I'):
        r"""
        Return a string representation of ``self``.

        INPUT:

        - ``base`` --  (Default: 10) The base to use for printing

        - ``truncate`` -- (Default: ``True``) Whether to print fewer
          digits than are available, to mask errors in the last bits.

        - ``istr`` -- (Default: ``I``) String representation of the complex unit

        EXAMPLES::

            sage: a = CC(pi + I*e)
            sage: a.str()
            '3.14159265358979 + 2.71828182845905*I'
            sage: a.str(truncate=False)
            '3.1415926535897931 + 2.7182818284590451*I'
            sage: a.str(base=2)
            '11.001001000011111101101010100010001000010110100011000 + 10.101101111110000101010001011000101000101011101101001*I'
            sage: CC(0.5 + 0.625*I).str(base=2)
            '0.10000000000000000000000000000000000000000000000000000 + 0.10100000000000000000000000000000000000000000000000000*I'
            sage: a.str(base=16)
            '3.243f6a8885a30 + 2.b7e151628aed2*I'
            sage: a.str(base=36)
            '3.53i5ab8p5fc + 2.puw5nggjf8f*I'
            sage: CC(0)
            0.000000000000000
            sage: CC.0.str(istr='%i')
            '1.00000000000000*%i'
        """
        s = ""
        if self.real() != 0:
            s = self.real().str(base, truncate=truncate)
        if self.imag() != 0:
            y  =  self.imag()
            if s!="":
                if y < 0:
                    s = s+" - "
                    y = -y
                else:
                    s = s+" + "
            s = s+"{ystr}*{istr}".format(ystr=y.str(base, truncate=truncate),
                    istr=istr)
        if len(s) == 0:
            s = self.real().str(base, truncate=truncate)
        return s

    def _latex_(self):
        r"""
        Method for converting ``self`` to a string with latex formatting.
        Called by the global function ``latex``.

        EXAMPLES::

            sage: a = ComplexNumber(2,1)
            sage: a
            2.00000000000000 + 1.00000000000000*I
            sage: latex(a)
            2.00000000000000 + 1.00000000000000i
            sage: a._latex_()
            '2.00000000000000 + 1.00000000000000i'

        ::

            sage: b = ComplexNumber(7,4,min_prec=16)
            sage: b
            7.000 + 4.000*I
            sage: latex(b)
            7.000 + 4.000i
            sage: b._latex_()
            '7.000 + 4.000i'

        ::

            sage: ComplexNumber(0).log()._latex_()
            '-\\infty'
        """
        import re
        s = self.str().replace('*I', 'i').replace('infinity','\\infty')
        return re.sub(r"e(-?\d+)", r" \\times 10^{\1}", s)

    def _pari_(self):
        r"""
        Coerces ``self`` into a PARI ``t_COMPLEX`` object,
        or a ``t_REAL`` if ``self`` is real.

        EXAMPLES:

        Coerce the object using the ``pari`` function::

            sage: a = ComplexNumber(2,1)
            sage: pari(a)
            2.00000000000000 + 1.00000000000000*I
            sage: pari(a).type()
            't_COMPLEX'
            sage: type(pari(a))
            <type 'sage.libs.pari.gen.gen'>
            sage: a._pari_()
            2.00000000000000 + 1.00000000000000*I
            sage: type(a._pari_())
            <type 'sage.libs.pari.gen.gen'>
            sage: a = CC(pi)
            sage: pari(a)
            3.14159265358979
            sage: pari(a).type()
            't_REAL'
            sage: a = CC(-2).sqrt()
            sage: pari(a)
            1.41421356237310*I
        """
        if self.is_real():
            return self.real()._pari_()
        return sage.libs.pari.all.pari.complex(self.real() or 0, self.imag())

    def _mpmath_(self, prec=None, rounding=None):
        """
        Returns an mpmath version of ``self``.

        .. NOTE::

           Currently, the rounding mode is ignored.

        EXAMPLES::

            sage: CC(1,2)._mpmath_()
            mpc(real='1.0', imag='2.0')
        """
        if prec is not None:
            from .complex_field import ComplexField
            return ComplexField(prec)(self)._mpmath_()
        from sage.libs.mpmath.all import make_mpc
        re = mpfr_to_mpfval(self.__re)
        im = mpfr_to_mpfval(self.__im)
        return make_mpc((re, im))

    cpdef _add_(self, right):
        """
        Add ``self`` to ``right``.

        EXAMPLES::

            sage: CC(2, 1)._add_(CC(1, -2))
            3.00000000000000 - 1.00000000000000*I
        """
        cdef ComplexNumber x
        x = self._new()
        mpfr_add(x.__re, self.__re, (<ComplexNumber>right).__re, rnd)
        mpfr_add(x.__im, self.__im, (<ComplexNumber>right).__im, rnd)
        return x

    cpdef _sub_(self, right):
        """
        Subtract ``right`` from ``self``.

        EXAMPLES::

            sage: CC(2, 1)._sub_(CC(1, -2))
            1.00000000000000 + 3.00000000000000*I
        """
        cdef ComplexNumber x
        x = self._new()
        mpfr_sub(x.__re, self.__re, (<ComplexNumber>right).__re, rnd)
        mpfr_sub(x.__im, self.__im, (<ComplexNumber>right).__im, rnd)
        return x

    cpdef _mul_(self, right):
        """
        Multiply ``self`` by ``right``.

        EXAMPLES::

            sage: CC(2, 1)._mul_(CC(1, -2))
            4.00000000000000 - 3.00000000000000*I
        """
        cdef ComplexNumber x
        x = self._new()
        cdef mpfr_t t0, t1
        mpfr_init2(t0, self._prec)
        mpfr_init2(t1, self._prec)
        mpfr_mul(t0, self.__re, (<ComplexNumber>right).__re, rnd)
        mpfr_mul(t1, self.__im, (<ComplexNumber>right).__im, rnd)
        mpfr_sub(x.__re, t0, t1, rnd)
        mpfr_mul(t0, self.__re, (<ComplexNumber>right).__im, rnd)
        mpfr_mul(t1, self.__im, (<ComplexNumber>right).__re, rnd)
        mpfr_add(x.__im, t0, t1, rnd)
        mpfr_clear(t0)
        mpfr_clear(t1)
        return x

    def norm(self):
        r"""
        Returns the norm of this complex number.

        If `c = a + bi` is a complex number, then the norm of `c` is defined as
        the product of `c` and its complex conjugate:

        .. MATH::

            \text{norm}(c)
            =
            \text{norm}(a + bi)
            =
            c \cdot \overline{c}
            =
            a^2 + b^2.

        The norm of a complex number is different from its absolute value.
        The absolute value of a complex number is defined to be the square
        root of its norm. A typical use of the complex norm is in the
        integral domain `\ZZ[i]` of Gaussian integers, where the norm of
        each Gaussian integer `c = a + bi` is defined as its complex norm.

        .. SEEALSO::

            - :func:`sage.misc.functional.norm`

            - :meth:`sage.rings.complex_double.ComplexDoubleElement.norm`

        EXAMPLES:

        This indeed acts as the square function when the
        imaginary component of ``self`` is equal to zero::

            sage: a = ComplexNumber(2,1)
            sage: a.norm()
            5.00000000000000
            sage: b = ComplexNumber(4.2,0)
            sage: b.norm()
            17.6400000000000
            sage: b^2
            17.6400000000000
        """
        return self.norm_c()

    cdef RealNumber norm_c(ComplexNumber self):
        cdef RealNumber x
        x = RealNumber(self._parent._real_field(), None)

        cdef mpfr_t t0, t1
        mpfr_init2(t0, self._prec)
        mpfr_init2(t1, self._prec)

        mpfr_mul(t0, self.__re, self.__re, rnd)
        mpfr_mul(t1, self.__im, self.__im, rnd)

        mpfr_add(x.value, t0, t1, rnd)

        mpfr_clear(t0)
        mpfr_clear(t1)
        return x

    cdef RealNumber abs_c(ComplexNumber self):
        cdef RealNumber x
        x = RealNumber(self._parent._real_field(), None)

        cdef mpfr_t t0, t1
        mpfr_init2(t0, self._prec)
        mpfr_init2(t1, self._prec)

        mpfr_mul(t0, self.__re, self.__re, rnd)
        mpfr_mul(t1, self.__im, self.__im, rnd)

        mpfr_add(x.value, t0, t1, rnd)
        mpfr_sqrt(x.value, x.value, rnd)

        mpfr_clear(t0)
        mpfr_clear(t1)
        return x

    cpdef _div_(self, right):
        """
        Divide ``self`` by ``right``.

        EXAMPLES::

            sage: CC(2, 1)._div_(CC(1, -2))
            1.00000000000000*I
        """
        cdef ComplexNumber x
        x = self._new()
        cdef mpfr_t a, b, t0, t1, right_nm
        mpfr_init2(t0, self._prec)
        mpfr_init2(t1, self._prec)
        mpfr_init2(a, self._prec)
        mpfr_init2(b, self._prec)
        mpfr_init2(right_nm, self._prec)

        mpfr_mul(t0, (<ComplexNumber>right).__re, (<ComplexNumber>right).__re, rnd)
        mpfr_mul(t1, (<ComplexNumber>right).__im, (<ComplexNumber>right).__im, rnd)
        mpfr_add(right_nm, t0, t1, rnd)

        mpfr_div(a, (<ComplexNumber>right).__re, right_nm, rnd)
        mpfr_div(b, (<ComplexNumber>right).__im, right_nm, rnd)

        ## Do this: x.__re =  a * self.__re + b * self.__im
        mpfr_mul(t0, a, self.__re, rnd)
        mpfr_mul(t1, b, self.__im, rnd)
        mpfr_add(x.__re, t0, t1, rnd)

        ## Do this: x.__im =  a * self.__im - b * self.__re
        mpfr_mul(t0, a, self.__im, rnd)
        mpfr_mul(t1, b, self.__re, rnd)
        mpfr_sub(x.__im, t0, t1, rnd)
        mpfr_clear(t0)
        mpfr_clear(t1)
        mpfr_clear(a)
        mpfr_clear(b)
        mpfr_clear(right_nm)
        return x

    def __rdiv__(self, left):
        r"""
        Returns the quotient of left with ``self``, that is:

        ``left/self``

        as a complex number.

        INPUT:

        - ``left`` - a complex number to divide by ``self``

        EXAMPLES::

            sage: a = ComplexNumber(2,0)
            sage: a.__rdiv__(CC(1))
            0.500000000000000
            sage: CC(1)/a
            0.500000000000000
        """
        return ComplexNumber(self._parent, left)/self

    def __pow__(self, right, modulus):
        r"""
        Raise ``self`` to the ``right`` expontent.

        This takes `a^b` and compues `\exp(b \log(a))`.

        EXAMPLES::

            sage: C.<i> = ComplexField(20)
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
        if isinstance(right, (int, long, integer.Integer)):
            return RingElement.__pow__(self, right)

        try:
            return (self.log()*right).exp()
        except TypeError:
            pass

        try:
            self = right.parent()(self)
            return self**right
        except AttributeError:
            raise TypeError

    def _magma_init_(self, magma):
        r"""
        EXAMPLES::

            sage: magma(CC([1, 2])) # indirect doctest, optional - magma
            1.00000000000000 + 2.00000000000000*$.1
            sage: v = magma(CC([1, 2])).sage(); v # indirect, optional - magma
            1.00000000000000 + 2.00000000000000*I
            sage: v.parent() # optional - magma
            Complex Field with 53 bits of precision

            sage: i = ComplexField(200).gen()
            sage: sqrt(i)
            0.70710678118654752440084436210484903928483593768847403658834 + 0.70710678118654752440084436210484903928483593768847403658834*I
            sage: magma(sqrt(i)) # indirect, optional - magma
            0.707106781186547524400844362104849039284835937688474036588340 + 0.707106781186547524400844362104849039284835937688474036588340*$.1
            sage: magma(i).Sqrt() # indirect, optional - magma
            0.707106781186547524400844362104849039284835937688474036588340 + 0.707106781186547524400844362104849039284835937688474036588340*$.1

            sage: magma(ComplexField(200)(1/3)) # indirect, optional - magma
            0.333333333333333333333333333333333333333333333333333333333333
        """
        real_string = self.real().str(truncate=False)
        imag_string = self.imag().str(truncate=False)
        digit_precision_bound = len(real_string)
        return "%s![%sp%s, %sp%s]" % (self.parent()._magma_init_(magma),
                                      real_string, digit_precision_bound,
                                      imag_string, digit_precision_bound)

    def __nonzero__(self):
        """
        Return ``True`` if ``self`` is not zero. This is an internal function;
        use :meth:`is_zero()` instead.

        EXAMPLES::

            sage: z = 1 + CC(I)
            sage: z.is_zero()
            False
        """
        return not (mpfr_zero_p(self.__re) and mpfr_zero_p(self.__im))

    def prec(self):
        """
        Return precision of this complex number.

        EXAMPLES::

            sage: i = ComplexField(2000).0
            sage: i.prec()
            2000
        """
        return self._parent.prec()

    def real(self):
        """
        Return real part of ``self``.

        EXAMPLES::

            sage: i = ComplexField(100).0
            sage: z = 2 + 3*i
            sage: x = z.real(); x
            2.0000000000000000000000000000
            sage: x.parent()
            Real Field with 100 bits of precision
            sage: z.real_part()
            2.0000000000000000000000000000
        """
        cdef RealNumber x
        x = RealNumber(self._parent._real_field(), None)
        mpfr_set(x.value, self.__re, rnd)
        return x

    real_part = real

    def imag(self):
        """
        Return imaginary part of ``self``.

        EXAMPLES::

            sage: i = ComplexField(100).0
            sage: z = 2 + 3*i
            sage: x = z.imag(); x
            3.0000000000000000000000000000
            sage: x.parent()
            Real Field with 100 bits of precision
            sage: z.imag_part()
            3.0000000000000000000000000000
        """
        cdef RealNumber x
        x = RealNumber(self._parent._real_field(), None)
        mpfr_set(x.value, self.__im, rnd)
        return x

    imag_part = imag

    def __neg__(self):
        r"""
        Method for computing the negative of ``self``.

        .. MATH::

            -(a + bi) = -a - bi

        EXAMPLES::

            sage: a = ComplexNumber(2,1)
            sage: -a
            -2.00000000000000 - 1.00000000000000*I
            sage: a.__neg__()
            -2.00000000000000 - 1.00000000000000*I
        """
        cdef ComplexNumber x
        x = self._new()
        mpfr_neg(x.__re, self.__re, rnd)
        mpfr_neg(x.__im, self.__im, rnd)
        return x

    def __pos__(self):
        r"""
        Method for computing the "positive" of ``self``.

        EXAMPLES::

            sage: a = ComplexNumber(2,1)
            sage: +a
            2.00000000000000 + 1.00000000000000*I
            sage: a.__pos__()
            2.00000000000000 + 1.00000000000000*I
        """
        return self

    def __abs__(self):
        r"""
        Method for computing the absolute value or modulus of ``self``

        .. MATH::

            `|a + bi| = sqrt(a^2 + b^2)`

        EXAMPLES:

        Note that the absolute value of a complex number with imaginary
        component equal to zero is the absolute value of the real component.

        ::

            sage: a = ComplexNumber(2,1)
            sage: abs(a)
            2.23606797749979
            sage: a.__abs__()
            2.23606797749979
            sage: float(sqrt(2^2 + 1^1))
            2.23606797749979

        ::

            sage: b = ComplexNumber(42,0)
            sage: abs(b)
            42.0000000000000
            sage: b.__abs__()
            42.0000000000000
            sage: b
            42.0000000000000
        """
        return self.abs_c()

    def __invert__(self):
        """
        Return the multiplicative inverse.

        EXAMPLES::

            sage: I = CC.0
            sage: a = ~(5+I)
            sage: a * (5+I)
            1.00000000000000
        """
        cdef ComplexNumber x
        x = self._new()

        cdef mpfr_t t0, t1
        mpfr_init2(t0, self._prec)
        mpfr_init2(t1, self._prec)

        mpfr_mul(t0, self.__re, self.__re, rnd)
        mpfr_mul(t1, self.__im, self.__im, rnd)

        mpfr_add(t0, t0, t1, rnd)         # now t0 is the norm
        mpfr_div(x.__re, self.__re, t0, rnd)   #     x.__re = self.__re/norm

        mpfr_neg(t1, self.__im, rnd)
        mpfr_div(x.__im, t1, t0, rnd)  #     x.__im = -self.__im/norm

        mpfr_clear(t0)
        mpfr_clear(t1)

        return x

    def __int__(self):
        r"""
        Method for converting ``self`` to type ``int``.

        Called by the ``int`` function. Note that calling this method returns
        an error since, in general, complex numbers cannot be coerced into
        integers.

        EXAMPLES::

            sage: a = ComplexNumber(2,1)
            sage: int(a)
            Traceback (most recent call last):
            ...
            TypeError: can't convert complex to int; use int(abs(z))
            sage: a.__int__()
            Traceback (most recent call last):
            ...
            TypeError: can't convert complex to int; use int(abs(z))
        """
        raise TypeError("can't convert complex to int; use int(abs(z))")

    def __long__(self):
        r"""
        Method for converting ``self`` to type ``long``.

        Called by the ``long`` function. Note that calling this method
        returns an error since, in general, complex numbers cannot be
        coerced into integers.

        EXAMPLES::

            sage: a = ComplexNumber(2,1)
            sage: long(a)
            Traceback (most recent call last):
            ...
            TypeError: can't convert complex to long; use long(abs(z))
            sage: a.__long__()
            Traceback (most recent call last):
            ...
            TypeError: can't convert complex to long; use long(abs(z))
        """
        raise TypeError("can't convert complex to long; use long(abs(z))")

    def __float__(self):
        r"""
        Method for converting ``self`` to type ``float``.

        Called by the ``float`` function.  This conversion will throw an error
        if the number has a nonzero imaginary part.

        EXAMPLES::

            sage: a = ComplexNumber(1, 0)
            sage: float(a)
            1.0
            sage: a = ComplexNumber(2,1)
            sage: float(a)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert 2.00000000000000 + 1.00000000000000*I to float; use abs() or real_part() as desired
            sage: a.__float__()
            Traceback (most recent call last):
            ...
            TypeError: unable to convert 2.00000000000000 + 1.00000000000000*I to float; use abs() or real_part() as desired
            sage: float(abs(ComplexNumber(1,1)))
            1.4142135623730951
        """
        if mpfr_zero_p(self.__im):
            return mpfr_get_d(self.__re, rnd)
        else:
            raise TypeError("unable to convert {!r} to float; use abs() or real_part() as desired".format(self))

    def __complex__(self):
        r"""
        Method for converting ``self`` to type ``complex``.

        Called by the ``complex`` function.

        EXAMPLES::

            sage: a = ComplexNumber(2,1)
            sage: complex(a)
            (2+1j)
            sage: type(complex(a))
            <type 'complex'>
            sage: a.__complex__()
            (2+1j)
        """
        return complex(mpfr_get_d(self.__re, rnd),
                       mpfr_get_d(self.__im, rnd))

    cpdef int _cmp_(left, right) except -2:
        """
        Compare ``left`` and ``right``.

        EXAMPLES::

            sage: cmp(CC(2, 1), CC(-1, 2))
            1
            sage: cmp(CC(2, 1), CC(2, 1))
            0
        """
        cdef int a, b
        a = mpfr_nan_p(left.__re)
        b = mpfr_nan_p((<ComplexNumber>right).__re)
        if a != b:
            return -1

        cdef int i
        i = mpfr_cmp(left.__re, (<ComplexNumber>right).__re)
        if i < 0:
            return -1
        elif i > 0:
            return 1
        i = mpfr_cmp(left.__im, (<ComplexNumber>right).__im)
        if i < 0:
            return -1
        elif i > 0:
            return 1
        return 0

    def multiplicative_order(self):
        """
        Return the multiplicative order of this complex number, if known,
        or raise a ``NotImplementedError``.

        EXAMPLES::

            sage: C.<i> = ComplexField()
            sage: i.multiplicative_order()
            4
            sage: C(1).multiplicative_order()
            1
            sage: C(-1).multiplicative_order()
            2
            sage: C(i^2).multiplicative_order()
            2
            sage: C(-i).multiplicative_order()
            4
            sage: C(2).multiplicative_order()
            +Infinity
            sage: w = (1+sqrt(-3.0))/2; w
            0.500000000000000 + 0.866025403784439*I
            sage: abs(w)
            1.00000000000000
            sage: w.multiplicative_order()
            Traceback (most recent call last):
            ...
            NotImplementedError: order of element not known
        """
        if self == 1:
            return integer.Integer(1)
        elif self == -1:
            return integer.Integer(2)
        elif self == self._parent.gen():
            return integer.Integer(4)
        elif self == -self._parent.gen():
            return integer.Integer(4)
        elif not self._multiplicative_order is None:
            return integer.Integer(self._multiplicative_order)
        elif abs(abs(self) - 1) > 0.1:  # clearly not a root of unity
            return infinity.infinity
        raise NotImplementedError("order of element not known")


    ########################################################################
    # Plotting
    ########################################################################

    def plot(self, **kargs):
        """
        Plots this complex number as a point in the plane

        The accepted options are the ones of :meth:`~sage.plot.point.point2d`.
        Type ``point2d.options`` to see all options.

        .. NOTE::

            Just wraps the sage.plot.point.point2d method

        EXAMPLES:

        You can either use the indirect::

            sage: z = CC(0,1)
            sage: plot(z)
            Graphics object consisting of 1 graphics primitive

        or the more direct::

            sage: z = CC(0,1)
            sage: z.plot()
            Graphics object consisting of 1 graphics primitive
        """
        return sage.plot.point.point2d((self.real(), self.imag()), **kargs)

    ########################################################################
    # Transcendental (and other) functions
    ########################################################################

    # Trig functions
    def arccos(self):
        """
        Return the arccosine of ``self``.

        EXAMPLES::

            sage: (1+CC(I)).arccos()
            0.904556894302381 - 1.06127506190504*I
        """
        return self._parent(self._pari_().acos())

    def arccosh(self):
        """
        Return the hyperbolic arccosine of ``self``.

        EXAMPLES::

            sage: (1+CC(I)).arccosh()
            1.06127506190504 + 0.904556894302381*I
        """
        return self._parent(self._pari_().acosh())

    def arcsin(self):
        """
        Return the arcsine of ``self``.

        EXAMPLES::

            sage: (1+CC(I)).arcsin()
            0.666239432492515 + 1.06127506190504*I
        """
        return self._parent(self._pari_().asin())

    def arcsinh(self):
        """
        Return the hyperbolic arcsine of ``self``.

        EXAMPLES::

            sage: (1+CC(I)).arcsinh()
            1.06127506190504 + 0.666239432492515*I
        """
        return self._parent(self._pari_().asinh())

    def arctan(self):
        """
        Return the arctangent of ``self``.

        EXAMPLES::

            sage: (1+CC(I)).arctan()
            1.01722196789785 + 0.402359478108525*I
        """
        return self._parent(self._pari_().atan())

    def arctanh(self):
        """
        Return the hyperbolic arctangent of ``self``.

        EXAMPLES::

            sage: (1+CC(I)).arctanh()
            0.402359478108525 + 1.01722196789785*I
        """
        return self._parent(self._pari_().atanh())

    def coth(self):
        """
        Return the hyperbolic cotangent of ``self``.

        EXAMPLES::

            sage: ComplexField(100)(1,1).coth()
            0.86801414289592494863584920892 - 0.21762156185440268136513424361*I
        """
        return ~(self.tanh())

    def arccoth(self):
        """
        Return the hyperbolic arccotangent of ``self``.

        EXAMPLES::

            sage: ComplexField(100)(1,1).arccoth()
            0.40235947810852509365018983331 - 0.55357435889704525150853273009*I
        """
        return (~self).arctanh()

    def csc(self):
        """
        Return the cosecant of ``self``.

        EXAMPLES::

            sage: ComplexField(100)(1,1).csc()
            0.62151801717042842123490780586 - 0.30393100162842645033448560451*I
        """
        return ~(self.sin())

    def csch(self):
        """
        Return the hyperbolic cosecant of ``self``.

        EXAMPLES::

            sage: ComplexField(100)(1,1).csch()
            0.30393100162842645033448560451 - 0.62151801717042842123490780586*I
        """
        return ~(self.sinh())

    def arccsch(self):
        """
        Return the hyperbolic arccosecant of ``self``.

        EXAMPLES::

            sage: ComplexField(100)(1,1).arccsch()
            0.53063753095251782601650945811 - 0.45227844715119068206365839783*I
        """
        return (~self).arcsinh()

    def sec(self):
        """
        Return the secant of ``self``.

        EXAMPLES::

            sage: ComplexField(100)(1,1).sec()
            0.49833703055518678521380589177 + 0.59108384172104504805039169297*I
        """
        return ~(self.cos())

    def sech(self):
        """
        Return the hyperbolic secant of ``self``.

        EXAMPLES::

            sage: ComplexField(100)(1,1).sech()
            0.49833703055518678521380589177 - 0.59108384172104504805039169297*I
        """
        return ~(self.cosh())

    def arcsech(self):
        """
        Return the hyperbolic arcsecant of ``self``.

        EXAMPLES::

            sage: ComplexField(100)(1,1).arcsech()
            0.53063753095251782601650945811 - 1.1185178796437059371676632938*I
        """
        return (~self).arccosh()

    def cotan(self):
        """
        Return the cotangent of ``self``.

        EXAMPLES::

            sage: (1+CC(I)).cotan()
            0.217621561854403 - 0.868014142895925*I
            sage: i = ComplexField(200).0
            sage: (1+i).cotan()
            0.21762156185440268136513424360523807352075436916785404091068 - 0.86801414289592494863584920891627388827343874994609327121115*I
            sage: i = ComplexField(220).0
            sage: (1+i).cotan()
            0.21762156185440268136513424360523807352075436916785404091068124239 - 0.86801414289592494863584920891627388827343874994609327121115071646*I
        """
        return ~(self.tan())

    def cos(self):
        """
        Return the cosine of ``self``.

        EXAMPLES::

            sage: (1+CC(I)).cos()
            0.833730025131149 - 0.988897705762865*I
        """
        # write self = a + i*b, then
        # cos(self) = cosh(b)*cos(a) - i*sinh(b)*sin(a)
        cdef ComplexNumber z
        z = self._new()
        cdef mpfr_t ch, sh
        mpfr_init2(sh, self._prec)
        mpfr_sinh(sh, self.__im, rnd)
        mpfr_init2(ch, self._prec)
        mpfr_sqr(ch, sh, rnd)
        mpfr_add_ui(ch, ch, 1, rnd)
        mpfr_sqrt(ch, ch, rnd)
        mpfr_neg(sh, sh, rnd)
        mpfr_sin_cos(z.__im, z.__re, self.__re, rnd)
        mpfr_mul(z.__re, z.__re, ch, rnd)
        mpfr_mul(z.__im, z.__im, sh, rnd)
        mpfr_clear(sh)
        mpfr_clear(ch)
        return z

    def cosh(self):
        """
        Return the hyperbolic cosine of ``self``.

        EXAMPLES::

            sage: (1+CC(I)).cosh()
            0.833730025131149 + 0.988897705762865*I
        """
        # write self = a + i*b, then
        # cosh(self) = cosh(a)*cos(b) + i*sinh(a)*sin(b)
        cdef ComplexNumber z
        z = self._new()
        cdef mpfr_t ch, sh
        mpfr_init2(sh, self._prec)
        mpfr_sinh(sh, self.__re, rnd)
        mpfr_init2(ch, self._prec)
        mpfr_sqr(ch, sh, rnd)
        mpfr_add_ui(ch, ch, 1, rnd)
        mpfr_sqrt(ch, ch, rnd)
        mpfr_sin_cos(z.__im, z.__re, self.__im, rnd)
        mpfr_mul(z.__re, z.__re, ch, rnd)
        mpfr_mul(z.__im, z.__im, sh, rnd)
        mpfr_clear(sh)
        mpfr_clear(ch)
        return z



    def eta(self, omit_frac=False):
        r"""
        Return the value of the Dedekind `\eta` function on ``self``,
        intelligently computed using `\mathbb{SL}(2,\ZZ)`
        transformations.

        The `\eta` function is

        .. MATH::

            \eta(z) = e^{\pi i z / 12} \prod_{n=1}^{\infty}(1-e^{2\pi inz})

        INPUT:

        -  ``self`` -- element of the upper half plane (if not,
           raises a ``ValueError``).

        -  ``omit_frac`` -- (bool, default: ``False``), if ``True``,
           omit the `e^{\pi i z / 12}` factor.

        OUTPUT: a complex number

        ALGORITHM: Uses the PARI C library.

        EXAMPLES:

        First we compute `\eta(1+i)`::

            sage: i = CC.0
            sage: z = 1+i; z.eta()
            0.742048775836565 + 0.198831370229911*I

        We compute eta to low precision directly from the definition::

            sage: z = 1 + i; z.eta()
            0.742048775836565 + 0.198831370229911*I
            sage: pi = CC(pi)        # otherwise we will get a symbolic result.
            sage: exp(pi * i * z / 12) * prod([1-exp(2*pi*i*n*z) for n in range(1,10)])
            0.742048775836565 + 0.198831370229911*I

        The optional argument allows us to omit the fractional part::

            sage: z = 1 + i
            sage: z.eta(omit_frac=True)
            0.998129069925959
            sage: prod([1-exp(2*pi*i*n*z) for n in range(1,10)])
            0.998129069925958 + 4.59099857829247e-19*I

        We illustrate what happens when `z` is not in the upper
        half plane::

            sage: z = CC(1)
            sage: z.eta()
            Traceback (most recent call last):
            ...
            ValueError: value must be in the upper half plane

        You can also use functional notation::

            sage: eta(1+CC(I))
            0.742048775836565 + 0.198831370229911*I
        """
        try:
            return self._parent(self._pari_().eta(not omit_frac))
        except sage.libs.pari.all.PariError:
            raise ValueError("value must be in the upper half plane")


    def sin(self):
        """
        Return the sine of ``self``.

        EXAMPLES::

            sage: (1+CC(I)).sin()
            1.29845758141598 + 0.634963914784736*I
        """
        # write self = a + i*b, then
        # sin(self) = cosh(b)*sin(a) + i*sinh(b)*cos(a)
        cdef ComplexNumber z
        z = self._new()
        cdef mpfr_t ch, sh
        mpfr_init2(sh, self._prec)
        mpfr_sinh(sh, self.__im, rnd)
        mpfr_init2(ch, self._prec)
        mpfr_sqr(ch, sh, rnd)
        mpfr_add_ui(ch, ch, 1, rnd)
        mpfr_sqrt(ch, ch, rnd)
        mpfr_sin_cos(z.__re, z.__im, self.__re, rnd)
        mpfr_mul(z.__re, z.__re, ch, rnd)
        mpfr_mul(z.__im, z.__im, sh, rnd)
        mpfr_clear(sh)
        mpfr_clear(ch)
        return z

    def sinh(self):
        """
        Return the hyperbolic sine of ``self``.

        EXAMPLES::

            sage: (1+CC(I)).sinh()
            0.634963914784736 + 1.29845758141598*I
        """
        # write self = a + i*b, then
        # sinh(self) = sinh(a)*cos(b) + i*cosh(a)*sin(b)
        cdef ComplexNumber z
        z = self._new()
        cdef mpfr_t ch, sh
        mpfr_init2(sh, self._prec)
        mpfr_sinh(sh, self.__re, rnd)
        mpfr_init2(ch, self._prec)
        mpfr_sqr(ch, sh, rnd)
        mpfr_add_ui(ch, ch, 1, rnd)
        mpfr_sqrt(ch, ch, rnd)
        mpfr_sin_cos(z.__im, z.__re, self.__im, rnd)
        mpfr_mul(z.__re, z.__re, sh, rnd)
        mpfr_mul(z.__im, z.__im, ch, rnd)
        mpfr_clear(sh)
        mpfr_clear(ch)
        return z

    def tan(self):
        """
        Return the tangent of ``self``.

        EXAMPLES::

            sage: (1+CC(I)).tan()
            0.271752585319512 + 1.08392332733869*I
        """
        # write self = a + i*b, then
        # tan(self) = [cos(a)*sin(a) + i*cosh(b)*sinh(b)]/[sinh^2(b)+cos^2(a)]
        cdef ComplexNumber z
        z = self._new()
        cdef mpfr_t ch, sh, c, s, a, b
        mpfr_init2(sh, self._prec)
        mpfr_sinh(sh, self.__im, rnd)
        mpfr_init2(ch, self._prec)
        mpfr_init2(a, self._prec)
        mpfr_sqr(a, sh, rnd)
        mpfr_add_ui(ch, a, 1, rnd)
        mpfr_sqrt(ch, ch, rnd)
        mpfr_init2(c, self._prec)
        mpfr_init2(s, self._prec)
        mpfr_sin_cos(s, c, self.__re, rnd)
        mpfr_init2(b, self._prec)
        mpfr_sqr(b, c, rnd)
        mpfr_add(a, a, b, rnd)
        mpfr_mul(z.__re, c, s, rnd)
        mpfr_div(z.__re, z.__re, a, rnd)
        mpfr_mul(z.__im, ch, sh, rnd)
        mpfr_div(z.__im, z.__im, a, rnd)
        mpfr_clear(sh)
        mpfr_clear(ch)
        mpfr_clear(c)
        mpfr_clear(s)
        mpfr_clear(b)
        mpfr_clear(a)
        return z


    def tanh(self):
        """
        Return the hyperbolic tangent of ``self``.

        EXAMPLES::

            sage: (1+CC(I)).tanh()
            1.08392332733869 + 0.271752585319512*I
        """
        # write self = a + i*b, then
        # tanh(self) = [cosh(a)*sinh(a) + i*cos(b)*sin(b)]/[sinh^2(a)+cos^2(b)]
        cdef ComplexNumber z
        z = self._new()
        cdef mpfr_t ch, sh, c, s, a, b
        mpfr_init2(sh, self._prec)
        mpfr_sinh(sh, self.__re, rnd)
        mpfr_init2(ch, self._prec)
        mpfr_init2(a, self._prec)
        mpfr_sqr(a, sh, rnd)
        mpfr_add_ui(ch, a, 1, rnd)
        mpfr_sqrt(ch, ch, rnd)
        mpfr_init2(c, self._prec)
        mpfr_init2(s, self._prec)
        mpfr_sin_cos(s, c, self.__im, rnd)
        mpfr_init2(b, self._prec)
        mpfr_sqr(b, c, rnd)
        mpfr_add(a, a, b, rnd)
        mpfr_mul(z.__im, c, s, rnd)
        mpfr_div(z.__im, z.__im, a, rnd)
        mpfr_mul(z.__re, ch, sh, rnd)
        mpfr_div(z.__re, z.__re, a, rnd)
        mpfr_clear(sh)
        mpfr_clear(ch)
        mpfr_clear(c)
        mpfr_clear(s)
        mpfr_clear(b)
        mpfr_clear(a)
        return z

    # Other special functions
    def agm(self, right, algorithm="optimal"):
        """
        Return the Arithmetic-Geometric Mean (AGM) of ``self`` and ``right``.

        INPUT:

        - ``right`` (complex) -- another complex number

        - ``algorithm`` (string, default "optimal") -- the algorithm to use
          (see below).

        OUTPUT:

        (complex) A value of the AGM of ``self`` and ``right``.  Note that
        this is a multi-valued function, and the algorithm used
        affects the value returned, as follows:

        - "pari": Call the sgm function from the pari library.

        - "optimal": Use the AGM sequence such that at each stage
              `(a,b)` is replaced by `(a_1,b_1)=((a+b)/2,\pm\sqrt{ab})`
              where the sign is chosen so that `|a_1-b_1|\le|a_1+b_1|`, or
              equivalently `\Re(b_1/a_1)\ge 0`.  The resulting limit is
              maximal among all possible values.

        - "principal": Use the AGM sequence such that at each stage
              `(a,b)` is replaced by `(a_1,b_1)=((a+b)/2,\pm\sqrt{ab})`
              where the sign is chosen so that `\Re(b_1)\ge 0` (the
              so-called principal branch of the square root).

        The values `AGM(a,0)`, `AGM(0,a)`, and `AGM(a,-a)` are all taken to be 0.

        EXAMPLES::

            sage: a = CC(1,1)
            sage: b = CC(2,-1)
            sage: a.agm(b)
            1.62780548487271 + 0.136827548397369*I
            sage: a.agm(b, algorithm="optimal")
            1.62780548487271 + 0.136827548397369*I
            sage: a.agm(b, algorithm="principal")
            1.62780548487271 + 0.136827548397369*I
            sage: a.agm(b, algorithm="pari")
            1.62780548487271 + 0.136827548397369*I

        An example to show that the returned value depends on the algorithm
        parameter::

            sage: a = CC(-0.95,-0.65)
            sage: b = CC(0.683,0.747)
            sage: a.agm(b, algorithm="optimal")
            -0.371591652351761 + 0.319894660206830*I
            sage: a.agm(b, algorithm="principal")
            0.338175462986180 - 0.0135326969565405*I
            sage: a.agm(b, algorithm="pari")
            -0.371591652351761 + 0.319894660206830*I
            sage: a.agm(b, algorithm="optimal").abs()
            0.490319232466314
            sage: a.agm(b, algorithm="principal").abs()
            0.338446122230459
            sage: a.agm(b, algorithm="pari").abs()
            0.490319232466314

        TESTS:

        An example which came up in testing::

            sage: I = CC(I)
            sage: a =  0.501648970493109 + 1.11877240294744*I
            sage: b =  1.05946309435930 + 1.05946309435930*I
            sage: a.agm(b)
            0.774901870587681 + 1.10254945079875*I

            sage: a = CC(-0.32599972608379413, 0.60395514542928641)
            sage: b = CC( 0.6062314525690593,  0.1425693337776659)
            sage: a.agm(b)
            0.199246281325876 + 0.478401702759654*I
            sage: a.agm(-a)
            0.000000000000000
            sage: a.agm(0)
            0.000000000000000
            sage: CC(0).agm(a)
            0.000000000000000

        Consistency::

            sage: a = 1 + 0.5*I
            sage: b = 2 - 0.25*I
            sage: a.agm(b) - ComplexField(100)(a).agm(b)
            0.000000000000000
            sage: ComplexField(200)(a).agm(b) - ComplexField(500)(a).agm(b)
            0.00000000000000000000000000000000000000000000000000000000000
            sage: ComplexField(500)(a).agm(b) - ComplexField(1000)(a).agm(b)
            0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
        """
        if algorithm=="pari":
            t = self._parent(right)._pari_()
            return self._parent(self._pari_().agm(t))

        cdef ComplexNumber a, b, a1, b1, d, e, res
        cdef mp_exp_t rel_prec
        cdef bint optimal = algorithm == "optimal"

        if optimal or algorithm == "principal":

            if not isinstance(right, ComplexNumber) or (<ComplexNumber>right)._parent is not self._parent:
                right = self._parent(right)

            res = self._new()

            if mpfr_zero_p(self.__re) and mpfr_zero_p(self.__im):
                return self
            elif mpfr_zero_p((<ComplexNumber>right).__re) and mpfr_zero_p((<ComplexNumber>right).__im):
                return right
            elif (mpfr_cmpabs(self.__re, (<ComplexNumber>right).__re) == 0 and
                  mpfr_cmpabs(self.__im, (<ComplexNumber>right).__im) == 0 and
                  mpfr_cmp(self.__re, (<ComplexNumber>right).__re) != 0 and
                  mpfr_cmp(self.__im, (<ComplexNumber>right).__im) != 0):
                # self = -right
                mpfr_set_ui(res.__re, 0, rnd)
                mpfr_set_ui(res.__im, 0, rnd)
                return res

            # Do the computations to a bit higher precicion so rounding error
            # won't obscure the termination condition.
            a = ComplexNumber(self._parent.to_prec(self._prec+5), None)
            b = a._new()
            a1 = a._new()
            b1 = a._new()

            d = a._new()
            if optimal:
                e = a._new()

            # Make copies so we don't mutate self or right.
            mpfr_set(a.__re, self.__re, rnd)
            mpfr_set(a.__im, self.__im, rnd)
            mpfr_set(b.__re, (<ComplexNumber>right).__re, rnd)
            mpfr_set(b.__im, (<ComplexNumber>right).__im, rnd)

            if optimal:
                mpfr_add(e.__re, a.__re, b.__re, rnd)
                mpfr_add(e.__im, a.__im, b.__im, rnd)

            while True:

                # a1 = (a+b)/2
                if optimal:
                    mpfr_swap(a1.__re, e.__re)
                    mpfr_swap(a1.__im, e.__im)
                else:
                    mpfr_add(a1.__re, a.__re, b.__re, rnd)
                    mpfr_add(a1.__im, a.__im, b.__im, rnd)
                mpfr_mul_2si(a1.__re, a1.__re, -1, rnd)
                mpfr_mul_2si(a1.__im, a1.__im, -1, rnd)

                # b1 = sqrt(a*b)
                mpfr_mul(d.__re, a.__re, b.__re, rnd)
                mpfr_mul(d.__im, a.__im, b.__im, rnd)
                mpfr_sub(b1.__re, d.__re, d.__im, rnd)
                mpfr_mul(d.__re, a.__re, b.__im, rnd)
                mpfr_mul(d.__im, a.__im, b.__re, rnd)
                mpfr_add(b1.__im, d.__re, d.__im, rnd)
                b1 = b1.sqrt() # this would be a *lot* of code duplication

                # d = a1 - b1
                mpfr_sub(d.__re, a1.__re, b1.__re, rnd)
                mpfr_sub(d.__im, a1.__im, b1.__im, rnd)
                if mpfr_zero_p(d.__re) and mpfr_zero_p(d.__im):
                    mpfr_set(res.__re, a1.__re, rnd)
                    mpfr_set(res.__im, a1.__im, rnd)
                    return res

                if optimal:
                    # e = a1+b1
                    mpfr_add(e.__re, a1.__re, b1.__re, rnd)
                    mpfr_add(e.__im, a1.__im, b1.__im, rnd)
                    if mpfr_zero_p(e.__re) and mpfr_zero_p(e.__im):
                        mpfr_set(res.__re, a1.__re, rnd)
                        mpfr_set(res.__im, a1.__im, rnd)
                        return res

                    # |e| < |d|
                    if cmp_abs(e, d) < 0:
                        mpfr_swap(d.__re, e.__re)
                        mpfr_swap(d.__im, e.__im)
                        mpfr_neg(b1.__re, b1.__re, rnd)
                        mpfr_neg(b1.__im, b1.__im, rnd)

                rel_prec = min_exp_t(max_exp(a1), max_exp(b1)) - max_exp(d)
                if rel_prec > self._prec:
                    mpfr_set(res.__re, a1.__re, rnd)
                    mpfr_set(res.__im, a1.__im, rnd)
                    return res

                # a, b = a1, b1
                mpfr_swap(a.__re, a1.__re)
                mpfr_swap(a.__im, a1.__im)
                mpfr_swap(b.__re, b1.__re)
                mpfr_swap(b.__im, b1.__im)

        raise ValueError("agm algorithm must be one of 'pari', 'optimal', 'principal'")

    def argument(self):
        r"""
        The argument (angle) of the complex number, normalized so that
        `-\pi < \theta \leq \pi`.

        EXAMPLES::

            sage: i = CC.0
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
        x = RealNumber(self._parent._real_field(), None)
        mpfr_atan2(x.value, self.__im, self.__re, rnd)
        return x


    def arg(self):
        """
        See :meth:`argument()`.

        EXAMPLES::

            sage: i = CC.0
            sage: (i^2).arg()
            3.14159265358979
        """
        return self.argument()

    def conjugate(self):
        """
        Return the complex conjugate of this complex number.

        EXAMPLES::

            sage: i = CC.0
            sage: (1+i).conjugate()
            1.00000000000000 - 1.00000000000000*I
        """
        cdef ComplexNumber x
        x = self._new()

        cdef mpfr_t i
        mpfr_init2(i, self._prec)
        mpfr_neg(i, self.__im, rnd)
        mpfr_set(x.__re, self.__re, rnd)
        mpfr_set(x.__im, i, rnd)
        mpfr_clear(i)
        return x

    def dilog(self):
        r"""
        Returns the complex dilogarithm of ``self``.

        The complex dilogarithm, or Spence's function, is defined by

        .. MATH::

            Li_2(z) = - \int_0^z \frac{\log|1-\zeta|}{\zeta} d(\zeta)
            = \sum_{k=1}^\infty \frac{z^k}{k}

        Note that the series definition can only be used for `|z| < 1`.

        EXAMPLES::

            sage: a = ComplexNumber(1,0)
            sage: a.dilog()
            1.64493406684823
            sage: float(pi^2/6)
            1.6449340668482262

        ::

            sage: b = ComplexNumber(0,1)
            sage: b.dilog()
            -0.205616758356028 + 0.915965594177219*I

        ::

            sage: c = ComplexNumber(0,0)
            sage: c.dilog()
            0.000000000000000
        """
        return self._parent(self._pari_().dilog())

    def exp(ComplexNumber self):
        r"""
        Compute `e^z` or `\exp(z)`.

        EXAMPLES::

            sage: i = ComplexField(300).0
            sage: z = 1 + i
            sage: z.exp()
            1.46869393991588515713896759732660426132695673662900872279767567631093696585951213872272450 + 2.28735528717884239120817190670050180895558625666835568093865811410364716018934540926734485*I
        """
        # write self = a + i*b, then
        # exp(self) = exp(a)*(cos(b) + i*sin(b))
        cdef ComplexNumber z
        z = self._new()
        cdef mpfr_t r
        mpfr_init2(r, self._prec)
        mpfr_exp(r, self.__re, rnd)
        mpfr_sin_cos(z.__im, z.__re, self.__im, rnd)
        mpfr_mul(z.__re, z.__re, r, rnd)
        mpfr_mul(z.__im, z.__im, r, rnd)
        mpfr_clear(r)
        return z

    def gamma(self):
        """
        Return the Gamma function evaluated at this complex number.

        EXAMPLES::

            sage: i = ComplexField(30).0
            sage: (1+i).gamma()
            0.49801567 - 0.15494983*I

        TESTS::

            sage: CC(0).gamma()
            Infinity

        ::

            sage: CC(-1).gamma()
            Infinity
        """
        try:
            return self._parent(self._pari_().gamma())
        except sage.libs.pari.all.PariError:
            from sage.rings.infinity import UnsignedInfinityRing
            return UnsignedInfinityRing.gen()

    def gamma_inc(self, t):
        """
        Return the incomplete Gamma function evaluated at this complex
        number.

        EXAMPLES::

            sage: C, i = ComplexField(30).objgen()
            sage: (1+i).gamma_inc(2 + 3*i)  # abs tol 2e-10
            0.0020969149 - 0.059981914*I
            sage: (1+i).gamma_inc(5)
            -0.0013781309 + 0.0065198200*I
            sage: C(2).gamma_inc(1 + i)
            0.70709210 - 0.42035364*I
            sage: CC(2).gamma_inc(5)
            0.0404276819945128

        TESTS:

        Check that :trac:`7099` is fixed::

            sage: C = ComplexField(400)
            sage: C(2 + I).gamma_inc(C(3 + I))  # abs tol 1e-120
            0.121515644664508695525971545977439666159749344176962379708992904126499444842886620664991650378432544392118359044438541515 + 0.101533909079826033296475736021224621546966200987295663190553587086145836461236284668967411665020429964946098113930918850*I

        """
        return self._parent(self._pari_().incgam(t, precision=self.prec()))

    def log(self,base=None):
        r"""
        Complex logarithm of `z` with branch chosen as follows: Write
        `z = \rho e^{i \theta}` with `-\pi < \theta <= pi`. Then
        `\mathrm{log}(z) = \mathrm{log}(\rho) + i \theta`.

        .. WARNING::

           Currently the real log is computed using floats, so there
           is potential precision loss.

        EXAMPLES::

            sage: a = ComplexNumber(2,1)
            sage: a.log()
            0.804718956217050 + 0.463647609000806*I
            sage: log(a.abs())
            0.804718956217050
            sage: a.argument()
            0.463647609000806

        ::

            sage: b = ComplexNumber(float(exp(42)),0)
            sage: b.log()
            41.99999999999971

        ::

            sage: c = ComplexNumber(-1,0)
            sage: c.log()
            3.14159265358979*I

        The option of a base is included for compatibility with other logs::

            sage: c = ComplexNumber(-1,0)
            sage: c.log(2)
            4.53236014182719*I

        If either component (real or imaginary) of the complex number
        is NaN (not a number), log will return the complex NaN::

            sage: c = ComplexNumber(NaN,2)
            sage: c.log()
            NaN - NaN*I

        """
        if mpfr_nan_p(self.__re):
            return ComplexNumber(self._parent,self.real(),self.real())
        if mpfr_nan_p(self.__im):
            return ComplexNumber(self._parent,self.imag(),self.imag())
        theta = self.argument()
        rho = abs(self)
        if base is None:
            return ComplexNumber(self._parent, rho.log(), theta)
        else:
            from real_mpfr import RealField
            return ComplexNumber(self._parent, rho.log()/RealNumber(RealField(self.prec()),base).log(), theta/RealNumber(RealField(self.prec()),base).log())

    def additive_order(self):
        """
        Return the additive order of ``self``.

        EXAMPLES::

            sage: CC(0).additive_order()
            1
            sage: CC.gen().additive_order()
            +Infinity
        """
        if self == 0:
            return 1
        else:
            return infinity.infinity

    def sqrt(self, all=False):
        """
        The square root function, taking the branch cut to be the negative
        real axis.

        INPUT:

        -  ``all`` - bool (default: ``False``); if ``True``, return a
           list of all square roots.

        EXAMPLES::

            sage: C.<i> = ComplexField(30)
            sage: i.sqrt()
            0.70710678 + 0.70710678*I
            sage: (1+i).sqrt()
            1.0986841 + 0.45508986*I
            sage: (C(-1)).sqrt()
            1.0000000*I
            sage: (1 + 1e-100*i).sqrt()^2
            1.0000000 + 1.0000000e-100*I
            sage: i = ComplexField(200).0
            sage: i.sqrt()
            0.70710678118654752440084436210484903928483593768847403658834 + 0.70710678118654752440084436210484903928483593768847403658834*I
        """
        cdef ComplexNumber z = self._new()
        if mpfr_zero_p(self.__im):
            if mpfr_sgn(self.__re) >= 0:
                mpfr_set_ui(z.__im, 0, rnd)
                mpfr_sqrt(z.__re, self.__re, rnd)
            else:
                mpfr_set_ui(z.__re, 0, rnd)
                mpfr_neg(z.__im, self.__re, rnd)
                mpfr_sqrt(z.__im, z.__im, rnd)
            if all:
                return [z, -z] if z else [z]
            else:
                return z
        # self = x + yi = (a+bi)^2
        # expand, substitute, solve
        # a^2 = (x + sqrt(x^2+y^2))/2
        cdef bint avoid_branch = mpfr_sgn(self.__re) < 0 and mpfr_cmpabs(self.__im, self.__re) < 0
        cdef mpfr_t a2
        mpfr_init2(a2, self._prec)
        mpfr_hypot(a2, self.__re, self.__im, rnd)
        if avoid_branch:
            # x + sqrt(x^2+y^2) numerically unstable for x near negative real axis
            # so we compute sqrt of (-z) and shift by i at the end
            mpfr_sub(a2, a2, self.__re, rnd)
        else:
            mpfr_add(a2, a2, self.__re, rnd)
        mpfr_mul_2si(a2, a2, -1, rnd)
        # a = sqrt(a2)
        mpfr_sqrt(z.__re, a2, rnd)
        # b = y/(2a)
        mpfr_div(z.__im, self.__im, z.__re, rnd)
        mpfr_mul_2si(z.__im, z.__im, -1, rnd)
        mpfr_clear(a2)
        if avoid_branch:
            mpfr_swap(z.__re, z.__im)
            # Note that y (hence b) was never negated, so we have z=i*sqrt(self).
            # if we were below the branch cut, we want the other branch
            if mpfr_sgn(self.__im) < 0:
                mpfr_neg(z.__re, z.__re, rnd)
                mpfr_neg(z.__im, z.__im, rnd)
        if all:
            return [z, -z]
        else:
            return z

    def nth_root(self, n, all=False):
        """
        The `n`-th root function.

        INPUT:

        -  ``all`` - bool (default: ``False``); if ``True``, return a
           list of all `n`-th roots.

        EXAMPLES::

            sage: a = CC(27)
            sage: a.nth_root(3)
            3.00000000000000
            sage: a.nth_root(3, all=True)
            [3.00000000000000, -1.50000000000000 + 2.59807621135332*I, -1.50000000000000 - 2.59807621135332*I]
            sage: a = ComplexField(20)(2,1)
            sage: [r^7 for r in a.nth_root(7, all=True)]
            [2.0000 + 1.0000*I, 2.0000 + 1.0000*I, 2.0000 + 1.0000*I, 2.0000 + 1.0000*I, 2.0000 + 1.0000*I, 2.0000 + 1.0001*I, 2.0000 + 1.0001*I]
        """
        if self.is_zero():
            return [self] if all else self

        cdef ComplexNumber z
        z = self._new()

        cdef RealNumber arg, rho
        cdef mpfr_t r
        rho = abs(self)
        arg = self.argument() / n
        mpfr_init2(r, self._prec)
        mpfr_root(r, rho.value, n, rnd)

        mpfr_sin_cos(z.__im, z.__re, arg.value, rnd)
        mpfr_mul(z.__re, z.__re, r, rnd)
        mpfr_mul(z.__im, z.__im, r, rnd)

        if not all:
            mpfr_clear(r)
            return z

        R = self._parent._real_field()
        cdef RealNumber theta
        theta = R.pi()*2/n
        zlist = [z]
        for k in range(1, n):
            z = self._new()
            arg += theta
            mpfr_sin_cos(z.__im, z.__re, arg.value, rnd)
            mpfr_mul(z.__re, z.__re, r, rnd)
            mpfr_mul(z.__im, z.__im, r, rnd)
            zlist.append(z)

        mpfr_clear(r)
        return zlist


    def is_square(self):
        r"""
        This function always returns true as `\CC` is algebraically closed.

        EXAMPLES::

            sage: a = ComplexNumber(2,1)
            sage: a.is_square()
            True

        `\CC` is algebraically closed, hence every element
        is a square::

            sage: b = ComplexNumber(5)
            sage: b.is_square()
            True
        """
        return True

    def is_real(self):
        """
        Return ``True`` if ``self`` is real, i.e. has imaginary part zero.

        EXAMPLES::

            sage: CC(1.23).is_real()
            True
            sage: CC(1+i).is_real()
            False
        """
        return (mpfr_zero_p(self.__im) != 0)

    def is_imaginary(self):
        """
        Return ``True`` if ``self`` is imaginary, i.e. has real part zero.

        EXAMPLES::

            sage: CC(1.23*i).is_imaginary()
            True
            sage: CC(1+i).is_imaginary()
            False
        """
        return (mpfr_zero_p(self.__re) != 0)

    def is_integer(self):
        """
        Return ``True`` if ``self`` is a integer

        EXAMPLES::

            sage: CC(3).is_integer()
            True
            sage: CC(1,2).is_integer()
            False
        """
        return self.is_real() and self.real() in ZZ

    def is_positive_infinity(self):
        r"""
        Check if ``self`` is `+\infty`.

        EXAMPLES::

            sage: CC(1, 2).is_positive_infinity()
            False
            sage: CC(oo, 0).is_positive_infinity()
            True
            sage: CC(0, oo).is_positive_infinity()
            False
        """
        return self.real().is_positive_infinity() and self.imag().is_zero()

    def is_negative_infinity(self):
        r"""
        Check if ``self`` is `-\infty`.

        EXAMPLES::

            sage: CC(1, 2).is_negative_infinity()
            False
            sage: CC(-oo, 0).is_negative_infinity()
            True
            sage: CC(0, -oo).is_negative_infinity()
            False
        """
        return self.real().is_negative_infinity() and self.imag().is_zero()

    def is_infinity(self):
        r"""
        Check if ``self`` is `\infty`.

        EXAMPLES::

            sage: CC(1, 2).is_infinity()
            False
            sage: CC(0, oo).is_infinity()
            True
        """
        return self.real().is_infinity() or self.imag().is_infinity()

    def zeta(self):
        """
        Return the Riemann zeta function evaluated at this complex number.

        EXAMPLES::

            sage: i = ComplexField(30).gen()
            sage: z = 1 + i
            sage: z.zeta()
            0.58215806 - 0.92684856*I
            sage: zeta(z)
            0.58215806 - 0.92684856*I

            sage: CC(1).zeta()
            Infinity
        """
        if mpfr_zero_p(self.__im) and mpfr_cmp_ui(self.__re, 1) == 0:
            import infinity
            return infinity.unsigned_infinity
        return self._parent(self._pari_().zeta())

    def algdep(self, n, **kwds):
        """
        Returns a polynomial of degree at most `n` which is
        approximately satisfied by this complex number. Note that the
        returned polynomial need not be irreducible, and indeed usually
        won't be if `z` is a good approximation to an algebraic
        number of degree less than `n`.

        ALGORITHM: Uses the PARI C-library algdep command.

        INPUT: Type algdep? at the top level prompt. All additional
        parameters are passed onto the top-level algdep command.

        EXAMPLE::

            sage: C = ComplexField()
            sage: z = (1/2)*(1 + sqrt(3.0) *C.0); z
            0.500000000000000 + 0.866025403784439*I
            sage: p = z.algdep(5); p
            x^3 + 1
            sage: p.factor()
            (x + 1) * (x^2 - x + 1)
            sage: z^2 - z + 1
            1.11022302462516e-16
        """
        import sage.arith.all
        return sage.arith.all.algdep(self, n, **kwds)

    def algebraic_dependancy( self, n ):
        """
        Returns a polynomial of degree at most `n` which is
        approximately satisfied by this complex number. Note that the
        returned polynomial need not be irreducible, and indeed usually
        won't be if `z` is a good approximation to an algebraic
        number of degree less than `n`.

        ALGORITHM: Uses the PARI C-library algdep command.

        INPUT: Type algdep? at the top level prompt. All additional
        parameters are passed onto the top-level algdep command.

        EXAMPLE::

            sage: C = ComplexField()
            sage: z = (1/2)*(1 + sqrt(3.0) *C.0); z
            0.500000000000000 + 0.866025403784439*I
            sage: p = z.algebraic_dependancy(5); p
            x^3 + 1
            sage: p.factor()
            (x + 1) * (x^2 - x + 1)
            sage: z^2 - z + 1
            1.11022302462516e-16
        """
        return self.algdep( n )

def make_ComplexNumber0( fld, mult_order, re, im ):
    """
    Create a complex number for pickling.

    EXAMPLES::

        sage: a = CC(1 + I)
        sage: loads(dumps(a)) == a # indirect doctest
        True
    """
    x = ComplexNumber( fld, re, im )
    x._set_multiplicative_order( mult_order )
    return x



def create_ComplexNumber(s_real, s_imag=None, int pad=0, min_prec=53):
    r"""
    Return the complex number defined by the strings ``s_real`` and
    ``s_imag`` as an element of ``ComplexField(prec=n)``,
    where `n` potentially has slightly more (controlled by pad) bits than
    given by `s`.

    INPUT:

    - ``s_real`` -- a string that defines a real number
      (or something whose string representation defines a number)

    - ``s_imag`` -- a string that defines a real number
      (or something whose string representation defines a number)

    - ``pad`` -- an integer at least 0.

    - ``min_prec`` -- number will have at least this many bits of precision,
      no matter what.

    EXAMPLES::

        sage: ComplexNumber('2.3')
        2.30000000000000
        sage: ComplexNumber('2.3','1.1')
        2.30000000000000 + 1.10000000000000*I
        sage: ComplexNumber(10)
        10.0000000000000
        sage: ComplexNumber(10,10)
        10.0000000000000 + 10.0000000000000*I
        sage: ComplexNumber(1.000000000000000000000000000,2)
        1.00000000000000000000000000 + 2.00000000000000000000000000*I
        sage: ComplexNumber(1,2.000000000000000000000)
        1.00000000000000000000 + 2.00000000000000000000*I

    ::

        sage: sage.rings.complex_number.create_ComplexNumber(s_real=2,s_imag=1)
        2.00000000000000 + 1.00000000000000*I

    TESTS:

    Make sure we've rounded up ``log(10,2)`` enough to guarantee
    sufficient precision (:trac:`10164`)::

        sage: s = "1." + "0"*10**6 + "1"
        sage: sage.rings.complex_number.create_ComplexNumber(s,0).real()-1 == 0
        False
        sage: sage.rings.complex_number.create_ComplexNumber(0,s).imag()-1 == 0
        False

    """
    if s_imag is None:
        s_imag = 0

    if not isinstance(s_real, str):
        s_real = str(s_real).strip()
    if not isinstance(s_imag, str):
        s_imag = str(s_imag).strip()
    #if base == 10:
    bits = max(int(LOG_TEN_TWO_PLUS_EPSILON*len(s_real)),
               int(LOG_TEN_TWO_PLUS_EPSILON*len(s_imag)))
    #else:
    #    bits = max(int(math.log(base,2)*len(s_imag)),int(math.log(base,2)*len(s_imag)))
    from .complex_field import ComplexField
    C = ComplexField(prec=max(bits+pad, min_prec))

    return ComplexNumber(C, s_real, s_imag)


cdef class RRtoCC(Map):

    cdef ComplexNumber _zero

    def __init__(self, RR, CC):
        """
        EXAMPLES::

            sage: from sage.rings.complex_number import RRtoCC
            sage: RRtoCC(RR, CC)
            Natural map:
              From: Real Field with 53 bits of precision
              To:   Complex Field with 53 bits of precision
        """
        Map.__init__(self, RR, CC)
        self._zero = ComplexNumber(CC, 0)
        self._repr_type_str = "Natural"

    cdef dict _extra_slots(self, dict _slots):
        """
        A helper for pickling and copying.

        INPUT:

        ``_slots`` -- a dictionary

        OUTPUT:

        The given dictionary, with zero added.

        EXAMPLES::

            sage: from sage.rings.complex_number import RRtoCC
            sage: f = RRtoCC(RR, CC)
            sage: g = copy(f) # indirect doctest
            sage: g
            Natural map:
              From: Real Field with 53 bits of precision
              To:   Complex Field with 53 bits of precision
        """
        _slots['_zero'] = self._zero
        return Map._extra_slots(self, _slots)

    cdef _update_slots(self, dict _slots):
        """
        A helper for unpickling and copying.

        INPUT:

        ``_slots`` -- a dictionary providing values for the c(p)def slots of self.

        EXAMPLES::

            sage: from sage.rings.complex_number import RRtoCC
            sage: RRtoCC(RR, CC)
            Natural map:
              From: Real Field with 53 bits of precision
              To:   Complex Field with 53 bits of precision
        """
        Map._update_slots(self, _slots)
        self._zero = _slots['_zero']

    cpdef Element _call_(self, x):
        """
        EXAMPLES::

            sage: from sage.rings.complex_number import RRtoCC
            sage: f = RRtoCC(RealField(100), ComplexField(10)) # indirect doctest
            sage: f(1/3)
            0.33
        """
        cdef ComplexNumber z = self._zero._new()
        mpfr_set(z.__re, (<RealNumber>x).value, rnd)
        mpfr_set_ui(z.__im, 0, rnd)
        return z


cdef class CCtoCDF(Map):

    cpdef Element _call_(self, x):
        """
        EXAMPLES::

            sage: from sage.rings.complex_number import CCtoCDF
            sage: f = CCtoCDF(CC, CDF) # indirect doctest
            sage: f(CC.0)
            1.0*I
            sage: f(exp(pi*CC.0/4))
            0.7071067811865476 + 0.7071067811865475*I
        """
        cdef ComplexDoubleElement z = <ComplexDoubleElement>ComplexDoubleElement.__new__(ComplexDoubleElement)
        z._complex.dat[0] = mpfr_get_d((<ComplexNumber>x).__re, MPFR_RNDN)
        z._complex.dat[1] = mpfr_get_d((<ComplexNumber>x).__im, MPFR_RNDN)
        return z


cdef inline mp_exp_t min_exp_t(mp_exp_t a, mp_exp_t b):
    return a if a < b else b

cdef inline mp_exp_t max_exp_t(mp_exp_t a, mp_exp_t b):
    return a if a > b else b

cdef inline mp_exp_t max_exp(ComplexNumber z):
    """
    Quickly return the maximum exponent of the real and complex parts of z,
    which is useful for estimating its magnitude.
    """
    if mpfr_zero_p(z.__im):
        return mpfr_get_exp(z.__re)
    elif mpfr_zero_p(z.__re):
        return mpfr_get_exp(z.__im)
    return max_exp_t(mpfr_get_exp(z.__re), mpfr_get_exp(z.__im))

cpdef int cmp_abs(ComplexNumber a, ComplexNumber b):
    """
    Returns -1, 0, or 1 according to whether `|a|` is less than, equal to, or
    greater than `|b|`.

    Optimized for non-close numbers, where the ordering can be determined by
    examining exponents.

    EXAMPLES::

        sage: from sage.rings.complex_number import cmp_abs
        sage: cmp_abs(CC(5), CC(1))
        1
        sage: cmp_abs(CC(5), CC(4))
        1
        sage: cmp_abs(CC(5), CC(5))
        0
        sage: cmp_abs(CC(5), CC(6))
        -1
        sage: cmp_abs(CC(5), CC(100))
        -1
        sage: cmp_abs(CC(-100), CC(1))
        1
        sage: cmp_abs(CC(-100), CC(100))
        0
        sage: cmp_abs(CC(-100), CC(1000))
        -1
        sage: cmp_abs(CC(1,1), CC(1))
        1
        sage: cmp_abs(CC(1,1), CC(2))
        -1
        sage: cmp_abs(CC(1,1), CC(1,0.99999))
        1
        sage: cmp_abs(CC(1,1), CC(1,-1))
        0
        sage: cmp_abs(CC(0), CC(1))
        -1
        sage: cmp_abs(CC(1), CC(0))
        1
        sage: cmp_abs(CC(0), CC(0))
        0
        sage: cmp_abs(CC(2,1), CC(1,2))
        0
    """
    if mpfr_zero_p(b.__re) and mpfr_zero_p(b.__im):
        return not ((mpfr_zero_p(a.__re) and mpfr_zero_p(a.__im)))
    elif (mpfr_zero_p(a.__re) and mpfr_zero_p(a.__im)):
        return -1
    cdef mp_exp_t exp_diff = max_exp(a) - max_exp(b)
    if exp_diff <= -2:
        return -1
    elif exp_diff >= 2:
        return 1

    cdef int res
    cdef mpfr_t abs_a, abs_b, tmp
    mpfr_init2(abs_a, mpfr_get_prec(a.__re))
    mpfr_init2(abs_b, mpfr_get_prec(b.__re))
    mpfr_init2(tmp, mpfr_get_prec(a.__re))

    mpfr_sqr(abs_a, a.__re, rnd)
    mpfr_sqr(tmp, a.__im, rnd)
    mpfr_add(abs_a, abs_a, tmp, rnd)

    mpfr_sqr(abs_b, b.__re, rnd)
    mpfr_sqr(tmp, b.__im, rnd)
    mpfr_add(abs_b, abs_b, tmp, rnd)

    res = mpfr_cmpabs(abs_a, abs_b)

    mpfr_clear(abs_a)
    mpfr_clear(abs_b)
    mpfr_clear(tmp)

    return res
