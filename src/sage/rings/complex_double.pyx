r"""
Double Precision Complex Numbers

Sage supports arithmetic using double-precision complex numbers. A
double-precision complex number is a complex number ``x + I*y`` with
`x`, `y` 64-bit (8 byte) floating point numbers (double precision).

The field :class:`ComplexDoubleField` implements the field
of all double-precision complex numbers. You can refer to this
field by the shorthand CDF. Elements of this field are of type
:class:`ComplexDoubleElement`. If `x` and `y` are coercible to
doubles, you can create a complex double element using
``ComplexDoubleElement(x,y)``. You can coerce more
general objects `z` to complex doubles by typing either
``ComplexDoubleField(x)`` or ``CDF(x)``.

EXAMPLES::

    sage: ComplexDoubleField()
    Complex Double Field
    sage: CDF
    Complex Double Field
    sage: type(CDF.0)
    <type 'sage.rings.complex_double.ComplexDoubleElement'>
    sage: ComplexDoubleElement(sqrt(2),3)
    1.41421356237 + 3.0*I
    sage: parent(CDF(-2))
    Complex Double Field

::

    sage: CC == CDF
    False
    sage: CDF is ComplexDoubleField()     # CDF is the shorthand
    True
    sage: CDF == ComplexDoubleField()
    True

The underlying arithmetic of complex numbers is implemented using
functions and macros in GSL (the GNU Scientific Library), and
should be very fast. Also, all standard complex trig functions,
log, exponents, etc., are implemented using GSL, and are also
robust and fast. Several other special functions, e.g. eta, gamma,
incomplete gamma, etc., are implemented using the PARI C library.

AUTHORS:

- William Stein (2006-09): first version

- Travis Scrimshaw (2012-10-18): Added doctests to get full coverage

- Jeroen Demeyer (2013-02-27): fixed all PARI calls (:trac:`14082`)
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#       Copyright (C) 2013 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import operator

from sage.misc.randstate cimport randstate, current_randstate

include 'sage/ext/interrupt.pxi'
include 'sage/ext/stdsage.pxi'
include 'sage/libs/pari/pari_err.pxi'

cdef extern from "<complex.h>":
    double complex csqrt(double complex)
    double cabs(double complex)

cimport sage.rings.ring

cimport sage.rings.integer

from sage.structure.element cimport RingElement, Element, ModuleElement, FieldElement
from sage.structure.parent  cimport Parent

from sage.libs.pari.gen cimport gen as pari_gen
from sage.libs.pari.pari_instance cimport PariInstance

import sage.libs.pari.pari_instance
cdef PariInstance pari = sage.libs.pari.pari_instance.pari

import complex_number

cdef RR, CC, RDF
import complex_field
CC = complex_field.ComplexField()

import real_mpfr
RR = real_mpfr.RealField()

from real_double import RealDoubleElement, RDF


from sage.structure.parent_gens import ParentWithGens
from sage.categories.morphism cimport Morphism

def is_ComplexDoubleField(x):
    """
    Return ``True`` if ``x`` is the complex double field.

    EXAMPLE::

        sage: from sage.rings.complex_double import is_ComplexDoubleField
        sage: is_ComplexDoubleField(CDF)
        True
        sage: is_ComplexDoubleField(ComplexField(53))
        False
    """
    return PY_TYPE_CHECK(x, ComplexDoubleField_class)

cdef class ComplexDoubleField_class(sage.rings.ring.Field):
    """
    An approximation to the field of complex numbers using double
    precision floating point numbers. Answers derived from calculations
    in this approximation may differ from what they would be if those
    calculations were performed in the true field of complex numbers.
    This is due to the rounding errors inherent to finite precision
    calculations.

    ALGORITHM:

    Arithmetic is done using GSL (the GNU Scientific Library).
    """
    def __init__(self):
        r"""
        Construct field of complex double precision numbers.

        EXAMPLE::

            sage: from sage.rings.complex_double import ComplexDoubleField_class
            sage: CDF == ComplexDoubleField_class()
            True
            sage: TestSuite(CDF).run(skip = ["_test_prod"])

        .. WARNING:: due to rounding errors, one can have `x^2 != x*x`::

            sage: x = CDF.an_element()
            sage: x
            1.0*I
            sage: x*x, x**2, x*x == x**2
            (-1.0, -1.0 + 1.2246...e-16*I, False)
        """
        from sage.categories.fields import Fields
        ParentWithGens.__init__(self, self, ('I',), normalize=False, category = Fields())
        self._populate_coercion_lists_()

    def __reduce__(self):
        """
        For pickling.

        EXAMPLES::

            sage: loads(dumps(CDF)) is CDF
            True
        """
        return ComplexDoubleField, ()

    cpdef bint is_exact(self) except -2:
        """
        Returns whether or not this field is exact, which is always ``False``.

        EXAMPLES::

            sage: CDF.is_exact()
            False
        """
        return False

    def __richcmp__(left, right, int op):
        """
        Rich comparison of ``left`` against ``right``.

        EXAMPLES::

            sage: cmp(CDF, CDF)
            0
        """
        return (<Parent>left)._richcmp_helper(right, op)

    cdef int _cmp_c_impl(left, Parent right) except -2:
        # There is only one CDF.
        return cmp(type(left),type(right))

    def __hash__(self):
        """
        Return the hash for ``self``.

        TEST::

            sage: hash(CDF) % 2^32 == hash(str(CDF)) % 2^32
            True
        """
        return 561162115
        #return hash(self.str())

    def characteristic(self):
        """
        Return the characteristic of the complex double field, which is 0.

        EXAMPLES::

            sage: CDF.characteristic()
            0
        """
        import integer
        return integer.Integer(0)

    def random_element(self, double xmin=-1, double xmax=1, double ymin=-1, double ymax=1):
        """
        Return a random element of this complex double field with real and
        imaginary part bounded by ``xmin``, ``xmax``, ``ymin``, ``ymax``.

        EXAMPLES::

            sage: CDF.random_element()
            -0.436810529675 + 0.736945423566*I
            sage: CDF.random_element(-10,10,-10,10)
            -7.08874026302 - 9.54135400334*I
            sage: CDF.random_element(-10^20,10^20,-2,2)
            -7.58765473764e+19 + 0.925549022839*I
        """
        cdef randstate rstate = current_randstate()
        global _CDF
        cdef ComplexDoubleElement z
        cdef double imag = (ymax-ymin)*rstate.c_rand_double() + ymin
        cdef double real = (xmax-xmin)*rstate.c_rand_double() + xmin
        z = PY_NEW(ComplexDoubleElement)
        z._complex = gsl_complex_rect(real, imag)
        return z

    def _repr_(self):
        """
        Print out this complex double field.

        EXAMPLES::

            sage: ComplexDoubleField() # indirect doctest
            Complex Double Field
            sage: CDF # indirect doctest
            Complex Double Field
        """
        return "Complex Double Field"

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.

        OUTPUT:

        - a string.

        TESTS::

            sage: print CDF._latex_()
            \Bold{C}
        """
        return r"\Bold{C}"

    def _cmp_(self, x):
        """
        Compare ``x`` to ``self``.

        EXAMPLES::

            sage: CDF == 5 # indirect doctest
            False
            sage: loads(dumps(CDF)) == CDF # indirect doctest
            True
        """
        if PY_TYPE_CHECK(x, ComplexDoubleField_class):
            return 0
        return cmp(type(self), type(x))

    def __call__(self, x, im=None):
        """
        Create a complex double using ``x`` and optionally an imaginary part
        ``im``.

        EXAMPLES::

            sage: CDF(0,1) # indirect doctest
            1.0*I
            sage: CDF(2/3) # indirect doctest
            0.666666666667
            sage: CDF(5) # indirect doctest
            5.0
            sage: CDF('i') # indirect doctest
            1.0*I
            sage: CDF(complex(2,-3)) # indirect doctest
            2.0 - 3.0*I
            sage: CDF(4.5) # indirect doctest
            4.5
            sage: CDF(1+I) # indirect doctest
            1.0 + 1.0*I
            sage: CDF(pari(1))
            1.0
            sage: CDF(pari("I"))
            1.0*I
            sage: CDF(pari("x^2 + x + 1").polroots()[0])
            -0.5 - 0.866025403784*I

        A ``TypeError`` is raised if the coercion doesn't make sense::

            sage: CDF(QQ['x'].0)
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce nonconstant polynomial to float

        One can convert back and forth between double precision complex
        numbers and higher-precision ones, though of course there may be
        loss of precision::

            sage: a = ComplexField(200)(-2).sqrt(); a
            1.4142135623730950488016887242096980785696718753769480731767*I
            sage: b = CDF(a); b
            1.41421356237*I
            sage: a.parent()(b)
            1.4142135623730951454746218587388284504413604736328125000000*I
            sage: a.parent()(b) == b
            True
            sage: b == CC(a)
            True
        """
        # We implement __call__ to gracefully accept the second argument.
        if im is not None:
            x = x, im
        return Parent.__call__(self, x)

    def _element_constructor_(self, x):
        """
        See ``__call__()``.

        EXAMPLES::

            sage: CDF((1,2)) # indirect doctest
            1.0 + 2.0*I
        """
        cdef pari_gen g
        if PY_TYPE_CHECK(x, ComplexDoubleElement):
            return x
        elif PY_TYPE_CHECK(x, tuple):
            return ComplexDoubleElement(x[0], x[1])
        elif isinstance(x, (float, int, long)):
            return ComplexDoubleElement(x, 0)
        elif isinstance(x, complex):
            return ComplexDoubleElement(x.real, x.imag)
        elif isinstance(x, complex_number.ComplexNumber):
            return ComplexDoubleElement(x.real(), x.imag())
        elif isinstance(x, pari_gen):
            g = x
            if typ(g.g) == t_COMPLEX:
                return ComplexDoubleElement(gtodouble(gel(g.g, 1)), gtodouble(gel(g.g, 2)))
            else:
                return ComplexDoubleElement(gtodouble(g.g), 0.0)
        elif isinstance(x, str):
            t = cdf_parser.parse_expression(x)
            if isinstance(t, float):
                return ComplexDoubleElement(t, 0)
            else:
                return t
        elif hasattr(x, '_complex_double_'):
            return x._complex_double_(self)
        else:
            return ComplexDoubleElement(x, 0)

    cpdef _coerce_map_from_(self, S):
        """
        Return the canonical coerce of `x` into the complex double field, if
        it is defined, otherwise raise a ``TypeError``.

        The rings that canonically coerce to the complex double field are:

        - the complex double field itself
        - anything that canonically coerces to real double field.
        - mathematical constants
        - the 53-bit mpfr complex field

        EXAMPLES::

            sage: CDF._coerce_(5) # indirect doctest
            5.0
            sage: CDF._coerce_(RDF(3.4))
            3.4

        Thus the sum of a CDF and a symbolic object is symbolic::

            sage: a = pi + CDF.0; a
            pi + 1.0*I
            sage: parent(a)
            Symbolic Ring

        TESTS::

            sage: CDF(1) + RR(1)
            2.0
            sage: CDF.0 - CC(1) - long(1) - RR(1) - QQbar(1)
            -4.0 + 1.0*I
            sage: CDF.has_coerce_map_from(ComplexField(20))
            False
        """
        from integer_ring import ZZ
        from rational_field import QQ
        from real_lazy import RLF
        from real_mpfr import RR, RealField_class
        from complex_field import ComplexField, ComplexField_class
        CC = ComplexField()
        from complex_number import CCtoCDF
        if S in [int, float, ZZ, QQ, RDF, RLF] or isinstance(S, RealField_class) and S.prec() >= 53:
            return FloatToCDF(S)
        elif RR.has_coerce_map_from(S):
            return FloatToCDF(RR) * RR.coerce_map_from(S)
        elif isinstance(S, ComplexField_class) and S.prec() >= 53:
            return CCtoCDF(S, self)
        elif CC.has_coerce_map_from(S):
            return CCtoCDF(CC, self) * CC.coerce_map_from(S)

    def _magma_init_(self, magma):
        r"""
        Return a string representation of ``self`` in the Magma language.

        EXAMPLES::

            sage: CDF._magma_init_(magma) # optional - magma
            'ComplexField(53 : Bits := true)'
            sage: magma(CDF) # optional - magma
            Complex field of precision 15
            sage: floor(RR(log(2**53, 10)))
            15
            sage: magma(CDF).sage() # optional - magma
            Complex Field with 53 bits of precision
        """
        return "ComplexField(%s : Bits := true)" % self.prec()

    def prec(self):
        """
        Return the precision of this complex double field (to be more
        similar to :class:`ComplexField`). Always returns 53.

        EXAMPLES::

            sage: CDF.prec()
            53
        """
        return 53

    precision=prec

    def to_prec(self, prec):
        """
        Returns the complex field to the specified precision. As doubles
        have fixed precision, this will only return a complex double field
        if prec is exactly 53.

        EXAMPLES::

            sage: CDF.to_prec(53)
            Complex Double Field
            sage: CDF.to_prec(250)
            Complex Field with 250 bits of precision
        """
        if prec == 53:
            return self
        else:
            from complex_field import ComplexField
            return ComplexField(prec)


    def gen(self, n=0):
        """
        Return the generator of the complex double field.

        EXAMPLES::

            sage: CDF.0
            1.0*I
            sage: CDF.gen(0)
            1.0*I
        """
        if n != 0:
            raise ValueError, "only 1 generator"
        return I

    def ngens(self):
        r"""
        The number of generators of this complex field as an `\RR`-algebra.

        There is one generator, namely ``sqrt(-1)``.

        EXAMPLES::

            sage: CDF.ngens()
            1
        """
        return 1

    def algebraic_closure(self):
        r"""
        Returns the algebraic closure of ``self``, i.e., the complex double
        field.

        EXAMPLES::

            sage: CDF.algebraic_closure()
            Complex Double Field
        """
        return self

    def real_double_field(self):
        """
        The real double field, which you may view as a subfield of this
        complex double field.

        EXAMPLES::

            sage: CDF.real_double_field()
            Real Double Field
        """
        import real_double
        return real_double.RDF

    def pi(self):
        r"""
        Returns `\pi` as a double precision complex number.

        EXAMPLES::

            sage: CDF.pi()
            3.14159265359
        """
        return self(3.1415926535897932384626433832)

    def construction(self):
        """
        Returns the functorial construction of ``self``, namely, algebraic
        closure of the real double field.

        EXAMPLES::

            sage: c, S = CDF.construction(); S
            Real Double Field
            sage: CDF == c(S)
            True
        """
        from sage.categories.pushout import AlgebraicClosureFunctor
        return (AlgebraicClosureFunctor(), self.real_double_field())

    def zeta(self, n=2):
        r"""
        Return a primitive `n`-th root of unity in this CDF, for
        `n \geq 1`.

        INPUT:

        -  ``n`` -- a positive integer (default: 2)

        OUTPUT: a complex `n`-th root of unity.

        EXAMPLES::

            sage: CDF.zeta(7)
            0.623489801859 + 0.781831482468*I
            sage: CDF.zeta(1)
            1.0
            sage: CDF.zeta()
            -1.0
            sage: CDF.zeta() == CDF.zeta(2)
            True

        ::

            sage: CDF.zeta(0.5)
            Traceback (most recent call last):
            ...
            ValueError: n must be a positive integer
            sage: CDF.zeta(0)
            Traceback (most recent call last):
            ...
            ValueError: n must be a positive integer
            sage: CDF.zeta(-1)
            Traceback (most recent call last):
            ...
            ValueError: n must be a positive integer
        """
        from integer import Integer
        try:
           n = Integer(n)
        except TypeError:
           raise ValueError, "n must be a positive integer"

        if n<1:
           raise ValueError, "n must be a positive integer"

        if n == 1:
            x = self(1)
        elif n == 2:
            x = self(-1)
        elif n >= 3:
            # Use De Moivre
            # e^(2*pi*i/n) = cos(2pi/n) + i *sin(2pi/n)
            pi = RDF.pi()
            z = 2*pi/n
            x = CDF(z.cos(), z.sin())
#         x._set_multiplicative_order( n ) # not implemented for CDF
        return x


cdef ComplexDoubleElement new_ComplexDoubleElement():
    """
    Creates a new (empty) :class:`ComplexDoubleElement`.
    """
    cdef ComplexDoubleElement z
    z = PY_NEW(ComplexDoubleElement)
    return z

def is_ComplexDoubleElement(x):
    """
    Return ``True`` if ``x`` is a :class:`ComplexDoubleElement`.

    EXAMPLES::

        sage: from sage.rings.complex_double import is_ComplexDoubleElement
        sage: is_ComplexDoubleElement(0)
        False
        sage: is_ComplexDoubleElement(CDF(0))
        True
    """
    return PY_TYPE_CHECK(x, ComplexDoubleElement)

cdef class ComplexDoubleElement(FieldElement):
    """
    An approximation to a complex number using double precision
    floating point numbers. Answers derived from calculations with such
    approximations may differ from what they would be if those
    calculations were performed with true complex numbers. This is due
    to the rounding errors inherent to finite precision calculations.
    """

    __array_interface__ = {'typestr': '=c16'}

    def __cinit__(self):
        r"""
        Initialize ``self`` as an element of `\CC`.

        EXAMPLES::

            sage: ComplexDoubleElement(1,-2) # indirect doctest
            1.0 - 2.0*I
        """
        self._parent = _CDF

    def __init__(self, real, imag):
        """
        Constructs an element of a complex double field with specified real
        and imaginary values.

        EXAMPLES::

            sage: ComplexDoubleElement(1,-2)
            1.0 - 2.0*I
        """
        self._complex = gsl_complex_rect(real, imag)

    def __reduce__(self):
        """
        For pickling.

        EXAMPLES::

            sage: a = CDF(-2.7, -3)
            sage: loads(dumps(a)) == a
            True
        """
        return (ComplexDoubleElement,
                (self._complex.dat[0], self._complex.dat[1]))

    cdef ComplexDoubleElement _new_c(self, gsl_complex x):
        """
        C-level code for creating a :class:`ComplexDoubleElement` from a
        ``gsl_complex``.
        """
        cdef ComplexDoubleElement z = <ComplexDoubleElement>PY_NEW(ComplexDoubleElement)
        z._complex = x
        return z

    cdef _new_from_gen(self, pari_gen g):
        """
        C-level code for creating a :class:`ComplexDoubleElement` from a
        PARI gen.

        INPUT:

        -  ``g`` -- A PARI ``gen``.

        OUTPUT: A ``ComplexDoubleElement``.
        """
        cdef gsl_complex x
        if typ(g.g) == t_COMPLEX:
            x = gsl_complex_rect(gtodouble(gel(g.g, 1)), gtodouble(gel(g.g, 2)))
        else:
            x = gsl_complex_rect(gtodouble(g.g), 0.0)
        return self._new_c(x)

    def __hash__(self):
        """
        Returns the hash of ``self``, which coincides with the python ``float``
        and ``complex`` (and often ``int``) types for ``self``.

        EXAMPLES::

            sage: hash(CDF(1.2)) == hash(1.2r)
            True
            sage: hash(CDF(-1))
            -2
            sage: hash(CDF(1.2, 1.3)) == hash(complex(1.2r, 1.3r))
            True
        """
        return hash(complex(self))

    def __richcmp__(left, right, int op):
        """
        Rich comparison between ``left`` and ``right``.

        EXAMPLES::

            sage: cmp(CDF(1.2), CDF(i))
            1
            sage: cmp(CDF(1), CDF(2))
            -1
            sage: cmp(CDF(1 + i), CDF(-1 - i))
            1
        """
        return (<Element>left)._richcmp(right, op)

    cdef int _cmp_c_impl(left, Element right) except -2:
        """
        We order the complex numbers in dictionary order by real parts then
        imaginary parts.

        This order, of course, does not respect the field structure, though
        it agrees with the usual order on the real numbers.

        EXAMPLES::

            sage: CDF(2,3) < CDF(3,1)
            True
            sage: CDF(2,3) > CDF(3,1)
            False
            sage: CDF(2,-1) < CDF(2,3)
            True

        It's dictionary order, not absolute value::

            sage: CDF(-1,3) < CDF(-1,-20)
            False

        Numbers are coerced before comparison::

            sage: CDF(3,5) < 7
            True
            sage: 4.3 > CDF(5,1)
            False
        """
        if left._complex.dat[0] < (<ComplexDoubleElement>right)._complex.dat[0]:
            return -1
        if left._complex.dat[0] > (<ComplexDoubleElement>right)._complex.dat[0]:
            return 1
        if left._complex.dat[1] < (<ComplexDoubleElement>right)._complex.dat[1]:
            return -1
        if left._complex.dat[1] > (<ComplexDoubleElement>right)._complex.dat[1]:
            return 1
        return 0

    def __getitem__(self, n):
        """
        Returns the real or imaginary part of ``self``.

        INPUT:

        -  ``n`` -- integer (either 0 or 1)

        Raises an ``IndexError`` if ``n`` is not 0 or 1.

        EXAMPLES::

            sage: P = CDF(2,3)
            sage: P[0]
            2.0
            sage: P[1]
            3.0
            sage: P[3]
            Traceback (most recent call last):
            ...
            IndexError: index n must be 0 or 1
        """
        if n >= 0 and n <= 1:
            return self._complex.dat[n]
        raise IndexError, "index n must be 0 or 1"

    def _magma_init_(self, magma):
        r"""
        Return the magma representation of ``self``.

        EXAMPLES::

            sage: CDF((1.2, 0.3))._magma_init_(magma) # optional - magma
            'ComplexField(53 : Bits := true)![1.2, 0.3]'
            sage: magma(CDF(1.2, 0.3)) # optional - magma # indirect doctest
            1.20000000000000 + 0.300000000000000*$.1
            sage: s = magma(CDF(1.2, 0.3)).sage(); s # optional - magma # indirect doctest
            1.20000000000000 + 0.300000000000000*I
            sage: s.parent() # optional - magma
            Complex Field with 53 bits of precision
        """
        return "%s![%s, %s]" % (self.parent()._magma_init_(magma), self.real(), self.imag())

    def prec(self):
        """
        Returns the precision of this number (to be more similar to
        :class:`ComplexNumber`). Always returns 53.

        EXAMPLES::

            sage: CDF(0).prec()
            53
        """
        return 53

    #######################################################################
    # Coercions
    #######################################################################

    def __int__(self):
        """
        Convert ``self`` to an ``int``.

        EXAMPLES::

            sage: int(CDF(1,1))
            Traceback (most recent call last):
            ...
            TypeError: can't convert complex to int; use int(abs(z))
            sage: int(abs(CDF(1,1)))
            1
        """
        raise TypeError, "can't convert complex to int; use int(abs(z))"

    def __long__(self):
        """
        Convert ``self`` to a ``long``.

        EXAMPLES::

            sage: long(CDF(1,1))
            Traceback (most recent call last):
            ...
            TypeError: can't convert complex to long; use long(abs(z))
            sage: long(abs(CDF(1,1)))
            1L
        """
        raise TypeError, "can't convert complex to long; use long(abs(z))"

    def __float__(self):
        """
        Method for converting ``self`` to type ``float``. Called by the
        ``float`` function.  This conversion will throw an error if
        the number has a nonzero imaginary part.

        EXAMPLES::

            sage: a = CDF(1, 0)
            sage: float(a)
            1.0
            sage: a = CDF(2,1)
            sage: float(a)
            Traceback (most recent call last):
            ...
            TypeError: Unable to convert 2.0 + 1.0*I to float; use abs() or real_part() as desired
            sage: a.__float__()
            Traceback (most recent call last):
            ...
            TypeError: Unable to convert 2.0 + 1.0*I to float; use abs() or real_part() as desired
            sage: float(abs(CDF(1,1)))
            1.4142135623730951
        """
        if self._complex.dat[1]==0:
            return float(self._complex.dat[0])
        raise TypeError, "Unable to convert %s to float; use abs() or real_part() as desired"%self

    def __complex__(self):
        """
        Convert ``self`` to python's ``complex`` object.

        EXAMPLES::

            sage: a = complex(2303,-3939)
            sage: CDF(a)
            2303.0 - 3939.0*I
            sage: complex(CDF(a))
            (2303-3939j)
        """
        return complex(self._complex.dat[0], self._complex.dat[1])

    def _interface_init_(self, I=None):
        """
        Returns ``self`` formatted as a string, suitable as input to another
        computer algebra system. (This is the default function used for
        exporting to other computer algebra systems.)

        EXAMPLES::

            sage: s1 = CDF(exp(I)); s1
            0.540302305868 + 0.841470984808*I
            sage: s1._interface_init_()
            '0.54030230586813977 + 0.84147098480789650*I'
            sage: s1 == CDF(gp(s1))
            True
        """
        # Sending to another computer algebra system is slow anyway, right?
        return CC(self)._interface_init_(I)

    def _maxima_init_(self, I=None):
        """
        Return a string representation of this complex number in the syntax of
        Maxima. That is, use ``%i`` to represent the complex unit.

        EXAMPLES::

            sage: CDF.0._maxima_init_()
            '1.0000000000000000*%i'
            sage: CDF(.5 + I)._maxima_init_()
            '0.50000000000000000 + 1.0000000000000000*%i'
        """
        return CC(self)._maxima_init_(I)

    def __repr__(self):
        """
        Return print version of ``self``.

        EXAMPLES::

            sage: a = CDF(2,-3); a # indirect doctest
            2.0 - 3.0*I
            sage: a^2
            -5.0 - 12.0*I
            sage: (1/CDF(0,0)).__repr__()
            'NaN + NaN*I'
            sage: CDF(oo,1)
            +infinity + 1.0*I
            sage: CDF(1,oo)
            1.0 + +infinity*I
            sage: CDF(1,-oo)
            1.0 - +infinity*I
            sage: CC(CDF(1,-oo))
            1.00000000000000 - +infinity*I
            sage: CDF(oo,oo)
            +infinity + +infinity*I
            sage: CC(CDF(oo,oo))
            +infinity + +infinity*I
            sage: CDF(0)
            0.0
        """
        if self._complex.dat[0]:
            # real part is nonzero
            s = double_to_str(self._complex.dat[0])
        else:
            # real part is zero
            if self._complex.dat[1]:   # imag is nonzero
                s = ''
            else:
                return double_to_str(self._complex.dat[0]) # imag is zero

        cdef double y = self._complex.dat[1]
        if y:
            if s != "":
                if y < 0:
                    s = s+" - "
                    y = -y
                else:
                    s = s+" + "
            t = double_to_str(y)
            s += t + "*I"
        return s

    def _latex_(self):
        """
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: CDF(1, 2)._latex_()
            '1.0 + 2.0i'
            sage: z = CDF(1,2)^100
            sage: z._latex_()
            '-6.44316469099 \\times 10^{34} - 6.11324130776 \\times 10^{34}i'
        """
        import re
        s = str(self).replace('*I', 'i')
        return re.sub(r"e\+?(-?\d+)", r" \\times 10^{\1}", s)

    cdef GEN _gen(self):
        cdef GEN y
        y = cgetg(3, t_COMPLEX)    # allocate space for a complex number
        set_gel(y, 1, pari.double_to_GEN(self._complex.dat[0]))
        set_gel(y, 2, pari.double_to_GEN(self._complex.dat[1]))
        return y

    def _pari_(self):
        """
        Return PARI version of ``self``.

        EXAMPLES::

            sage: CDF(1,2)._pari_()
            1.00000000000000 + 2.00000000000000*I
            sage: pari(CDF(1,2))
            1.00000000000000 + 2.00000000000000*I
        """
        pari_catch_sig_on()
        return pari.new_gen(self._gen())

    #######################################################################
    # Arithmetic
    #######################################################################

    cpdef ModuleElement _add_(self, ModuleElement right):
        """
        Add ``self`` and ``right``.

        EXAMPLES::

            sage: CDF(2,-3)._add_(CDF(1,-2))
            3.0 - 5.0*I
        """
        return self._new_c(gsl_complex_add(self._complex,
                                           (<ComplexDoubleElement>right)._complex))

    cpdef ModuleElement _sub_(self, ModuleElement right):
        """
        Subtract ``self`` and ``right``.

        EXAMPLES::

            sage: CDF(2,-3)._sub_(CDF(1,-2))
            1.0 - 1.0*I
        """
        return self._new_c(gsl_complex_sub(self._complex,
                                (<ComplexDoubleElement>right)._complex))

    cpdef RingElement _mul_(self, RingElement right):
        """
        Multiply ``self`` and ``right``.

        EXAMPLES::

            sage: CDF(2,-3)._mul_(CDF(1,-2))
            -4.0 - 7.0*I
        """
        return self._new_c(gsl_complex_mul(self._complex,
                       (<ComplexDoubleElement>right)._complex))

    cpdef RingElement _div_(self, RingElement right):
        """
        Divide ``self`` by ``right``.

        EXAMPLES::

            sage: CDF(2,-3)._div_(CDF(1,-2))
            1.6 + 0.2*I
        """
        return self._new_c(gsl_complex_div(self._complex, (<ComplexDoubleElement>right)._complex))

    def __invert__(self):
        r"""
        This function returns the inverse, or reciprocal, of the complex
        number `z`:

        .. MATH::

            1/z = (x - i y)/(x^2 + y^2).

        EXAMPLES::

            sage: ~CDF(2,1)
            0.4 - 0.2*I
            sage: 1/CDF(2,1)
            0.4 - 0.2*I

        The inverse of 0 is ``NaN`` (it doesn't raise an exception)::

            sage: ~(0*CDF(0,1))
            NaN + NaN*I
        """
        return self._new_c(gsl_complex_inverse(self._complex))

    cpdef ModuleElement _neg_(self):
        """
        This function returns the negative of the complex number `z`:

        .. MATH::

            -z = (-x) + i(-y).

        EXAMPLES::

            sage: -CDF(2,1) # indirect doctest
            -2.0 - 1.0*I
        """
        return self._new_c(gsl_complex_negative(self._complex))

    def conjugate(self):
        r"""
        This function returns the complex conjugate of the complex number `z`:

        .. MATH::

            \overline{z} = x - i y.

        EXAMPLES::

            sage: z = CDF(2,3); z.conjugate()
            2.0 - 3.0*I
        """
        return self._new_c(gsl_complex_conjugate(self._complex))

    def conj(self):
        r"""
        This function returns the complex conjugate of the complex number `z`:

        .. MATH::

            \overline{z} = x - i y.

        EXAMPLES::

            sage: z = CDF(2,3); z.conj()
            2.0 - 3.0*I
        """
        return self._new_c(gsl_complex_conjugate(self._complex))

    #######################################################################
    # Properties of Complex Numbers
    #######################################################################

    def arg(self):
        r"""
        This function returns the argument of ``self``, the complex number
        `z`, denoted by `\arg(z)`, where `-\pi < \arg(z) <= \pi`.

        EXAMPLES::

            sage: CDF(1,0).arg()
            0.0
            sage: CDF(0,1).arg()
            1.57079632679
            sage: CDF(0,-1).arg()
            -1.57079632679
            sage: CDF(-1,0).arg()
            3.14159265359
        """
        return RealDoubleElement(gsl_complex_arg(self._complex))

    def __abs__(self):
        """
        This function returns the magnitude of the complex number `z`, `|z|`.

        EXAMPLES::

            sage: abs(CDF(1,2)) # indirect doctest
            2.2360679775
            sage: abs(CDF(1,0)) # indirect doctest
            1.0
            sage: abs(CDF(-2,3))   # slightly random-ish arch dependent output
            3.6055512754639891
        """
        return RealDoubleElement(gsl_complex_abs(self._complex))

    def abs(self):
        """
        This function returns the magnitude `|z|` of the complex number `z`.

        .. SEEALSO::

            - :meth:`norm`

        EXAMPLES::

            sage: CDF(2,3).abs()   # slightly random-ish arch dependent output
            3.6055512754639891
        """
        return RealDoubleElement(gsl_complex_abs(self._complex))

    def argument(self):
        r"""
        This function returns the argument of the ``self``, the complex number
        `z`, in the interval `-\pi < arg(z) \leq \pi`.

        EXAMPLES::

            sage: CDF(6).argument()
            0.0
            sage: CDF(i).argument()
            1.57079632679
            sage: CDF(-1).argument()
            3.14159265359
            sage: CDF(-1 - 0.000001*i).argument()
            -3.14159165359
        """
        return RealDoubleElement(gsl_complex_arg(self._complex))

    def abs2(self):
        """
        This function returns the squared magnitude `|z|^2` of the complex
        number `z`, otherwise known as the complex norm.

        .. SEEALSO::

            - :meth:`norm`

        EXAMPLES::

            sage: CDF(2,3).abs2()
            13.0
        """
        return RealDoubleElement(gsl_complex_abs2(self._complex))

    def norm(self):
        r"""
        This function returns the squared magnitude `|z|^2` of the complex
        number `z`, otherwise known as the complex norm. If `c = a + bi`
        is a complex number, then the norm of `c` is defined as the product of
        `c` and its complex conjugate:

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

            - :meth:`abs`

            - :meth:`abs2`

            - :func:`sage.misc.functional.norm`

            - :meth:`sage.rings.complex_number.ComplexNumber.norm`

        EXAMPLES::

            sage: CDF(2,3).norm()
            13.0
        """
        return RealDoubleElement(gsl_complex_abs2(self._complex))

    def logabs(self):
        r"""
        This function returns the natural logarithm of the magnitude of the
        complex number `z`, `\log|z|`.

        This allows for an accurate evaluation of `\log|z|` when `|z|` is
        close to `1`. The direct evaluation of ``log(abs(z))`` would lead
        to a loss of precision in this case.

        EXAMPLES::

            sage: CDF(1.1,0.1).logabs()
            0.0994254293726
            sage: log(abs(CDF(1.1,0.1)))
            0.0994254293726

        ::

            sage: log(abs(ComplexField(200)(1.1,0.1)))
            0.099425429372582595066319157757531449594489450091985182495705
        """
        return RealDoubleElement(gsl_complex_logabs(self._complex))

    def real(self):
        """
        Return the real part of this complex double.

        EXAMPLES::

            sage: a = CDF(3,-2)
            sage: a.real()
            3.0
            sage: a.real_part()
            3.0
        """
        return RealDoubleElement(self._complex.dat[0])

    real_part = real

    def imag(self):
        """
        Return the imaginary part of this complex double.

        EXAMPLES::

            sage: a = CDF(3,-2)
            sage: a.imag()
            -2.0
            sage: a.imag_part()
            -2.0
        """
        return RealDoubleElement(self._complex.dat[1])

    imag_part = imag

    def parent(self):
        """
        Return the complex double field, which is the parent of ``self``.

        EXAMPLES::

            sage: a = CDF(2,3)
            sage: a.parent()
            Complex Double Field
            sage: parent(a)
            Complex Double Field
        """
        return CDF


    #######################################################################
    # Elementary Complex Functions
    #######################################################################
    def sqrt(self, all=False, **kwds):
        r"""
        The square root function.

        INPUT:

        -  ``all`` - bool (default: ``False``); if ``True``, return a
           list of all square roots.

        If all is ``False``, the branch cut is the negative real axis. The
        result always lies in the right half of the complex plane.

        EXAMPLES:

        We compute several square roots::

            sage: a = CDF(2,3)
            sage: b = a.sqrt(); b
            1.67414922804 + 0.89597747613*I
            sage: b^2
            2.0 + 3.0*I
            sage: a^(1/2)
            1.67414922804 + 0.89597747613*I

        We compute the square root of -1::

            sage: a = CDF(-1)
            sage: a.sqrt()
            1.0*I

        We compute all square roots::

            sage: CDF(-2).sqrt(all=True)
            [1.41421356237*I, -1.41421356237*I]
            sage: CDF(0).sqrt(all=True)
            [0.0]
        """
        z = self._new_c(gsl_complex_sqrt(self._complex))
        if all:
             if z.is_zero():
                 return [z]
             else:
                 return [z, -z]
        return z

    def nth_root(self, n, all=False):
        """
        The ``n``-th root function.

        INPUT:

        -  ``all`` -- bool (default: ``False``); if ``True``, return a
           list of all ``n``-th roots.

        EXAMPLES::

            sage: a = CDF(125)
            sage: a.nth_root(3)
            5.0
            sage: a = CDF(10, 2)
            sage: [r^5 for r in a.nth_root(5, all=True)]
            [10.0 + 2.0*I, 10.0 + 2.0*I, 10.0 + 2.0*I, 10.0 + 2.0*I, 10.0 + 2.0*I]
            sage: abs(sum(a.nth_root(111, all=True))) # random but close to zero
            6.00659385991e-14
        """
        if not self:
            return [self] if all else self
        arg = self.argument() / n
        abs = self.abs().nth_root(n)
        z = ComplexDoubleElement(abs * arg.cos(), abs*arg.sin())
        if all:
            zeta = self._parent.zeta(n)
            return [z * zeta**k for k in range(n)]
        else:
            return z


    def is_square(self):
        r"""
        This function always returns ``True`` as `\CC` is algebraically closed.

        EXAMPLES::

            sage: CDF(-1).is_square()
            True
        """
        return True

    def _pow_(self, ComplexDoubleElement a):
        """
        The function returns the complex number `z` raised to the
        complex power `a`, `z^a`.

        INPUT:

        -  ``a`` - a :class:`ComplexDoubleElement`

        OUTPUT: :class:`ComplexDoubleElement`

        EXAMPLES::

            sage: a = CDF(1,1); b = CDF(2,3)
            sage: a._pow_(b)
            -0.163450932107 + 0.0960049836089*I
        """
        return self._new_c(gsl_complex_pow(self._complex, a._complex))

    def __pow__(z, a, dummy):
        r"""
        The function returns the complex number `z` raised to the
        complex power `a`, `z^a`.

        This is computed as `\exp(\log(z)*a)` using complex
        logarithms and complex exponentials.

        EXAMPLES::

            sage: a = CDF(1,1); b = CDF(2,3)
            sage: c = a^b; c # indirect doctest
            -0.163450932107 + 0.0960049836089*I
            sage: c^(1/b)
            1.0 + 1.0*I

        We compute the cube root of `-1` then cube it and observe a
        rounding error::

            sage: a = CDF(-1)^(1/3); a
            0.5 + 0.866025403784*I
            sage: a^3                  # slightly random-ish arch dependent output
            -1.0 + 1.22460635382e-16*I

        We raise to symbolic powers::

            sage: x, n = var('x, n')
            sage: CDF(1.2)^x
            1.2^x
            sage: CDF(1.2)^(x^n + n^x)
            1.2^(n^x + x^n)
        """
        try:
            return z._pow_(a)
        except AttributeError:
            # z is not a complex number
            return CDF(z)._pow_(a)
        except TypeError:
            # a is not a complex number
            try:
                return z._pow_(CDF(a))
            except TypeError:
                try:
                    return a.parent()(z)**a
                except AttributeError:
                    raise TypeError


    def exp(self):
        r"""
        This function returns the complex exponential of the complex number
        `z`, `\exp(z)`.

        EXAMPLES::

            sage: CDF(1,1).exp()
            1.46869393992 + 2.28735528718*I

        We numerically verify a famous identity to the precision of a double::

            sage: z = CDF(0, 2*pi); z
            6.28318530718*I
            sage: exp(z)         # somewhat random-ish output depending on platform
            1.0 - 2.44921270764e-16*I
        """
        return self._new_c(gsl_complex_exp(self._complex))

    def log(self, base=None):
        r"""
        This function returns the complex natural logarithm to the given
        base of the complex number `z`, `\log(z)`. The
        branch cut is the negative real axis.

        INPUT:

        -  ``base`` - default: `e`, the base of the natural logarithm

        EXAMPLES::

            sage: CDF(1,1).log()
            0.34657359028 + 0.785398163397*I

        This is the only example different from the GSL::

            sage: CDF(0,0).log()
            -infinity
        """
        if self == 0:
            return RDF(0).log()
        if base is None:
            return self._new_c(gsl_complex_log(self._complex))
        cdef ComplexDoubleElement z
        try:
            z = base
        except TypeError:
            z = CDF(base)
        return self._new_c(gsl_complex_log_b(self._complex, z._complex))

    def log10(self):
        r"""
        This function returns the complex base-10 logarithm of the complex
        number `z`, `\log_{10}(z)`.

        The branch cut is the negative real axis.

        EXAMPLES::

            sage: CDF(1,1).log10()
            0.150514997832 + 0.34109408846*I
        """
        if self == 0:
            return RDF(0).log()
        return self._new_c(gsl_complex_log10(self._complex))

    def log_b(self, b):
        r"""
        This function returns the complex base-`b` logarithm of the
        complex number `z`, `\log_b(z)`. This quantity is
        computed as the ratio `\log(z)/\log(b)`.

        The branch cut is the negative real axis.

        EXAMPLES::

            sage: CDF(1,1).log_b(10)
            0.150514997832 + 0.34109408846*I
        """
        cdef ComplexDoubleElement _b
        if self == 0:
            return RDF(0).log()
        try:
            _b = b
        except TypeError:
            _b = CDF(b)
        return self._new_c(gsl_complex_log_b(self._complex, _b._complex))

    #######################################################################
    # Complex Trigonometric Functions
    #######################################################################
    def sin(self):
        r"""
        This function returns the complex sine of the complex number `z`:

        .. MATH::

            \sin(z) = \frac{e^{iz} - e^{-iz}}{2i}.

        EXAMPLES::

            sage: CDF(1,1).sin()
            1.29845758142 + 0.634963914785*I
        """
        return self._new_c(gsl_complex_sin(self._complex))

    def cos(self):
        r"""
        This function returns the complex cosine of the complex number `z`:

        .. MATH::

            \cos(z) = \frac{e^{iz} + e^{-iz}}{2}

        EXAMPLES::

            sage: CDF(1,1).cos()
            0.833730025131 - 0.988897705763*I
        """
        return self._new_c(gsl_complex_cos(self._complex))

    def tan(self):
        r"""
        This function returns the complex tangent of the complex number `z`:

        .. MATH::

            \tan(z) = \frac{\sin(z)}{\cos(z)}.

        EXAMPLES::

            sage: CDF(1,1).tan()
            0.27175258532 + 1.08392332734*I
        """
        return self._new_c(gsl_complex_tan(self._complex))

    def sec(self):
        r"""
        This function returns the complex secant of the complex number `z`:

        .. MATH::

            {\rm sec}(z) = \frac{1}{\cos(z)}.

        EXAMPLES::

            sage: CDF(1,1).sec()
            0.498337030555 + 0.591083841721*I
        """
        return self._new_c(gsl_complex_sec(self._complex))

    def csc(self):
        r"""
        This function returns the complex cosecant of the complex number `z`:

        .. MATH::

            \csc(z) = \frac{1}{\sin(z)}.

        EXAMPLES::

            sage: CDF(1,1).csc()
            0.62151801717 - 0.303931001628*I
        """
        return self._new_c(gsl_complex_csc(self._complex))

    def cot(self):
        r"""
        This function returns the complex cotangent of the complex number `z`:

        .. MATH::

            \cot(z) = \frac{1}{\tan(z)}.

        EXAMPLES::

            sage: CDF(1,1).cot()
            0.217621561854 - 0.868014142896*I
        """
        return self._new_c(gsl_complex_cot(self._complex))

    #######################################################################
    # Inverse Complex Trigonometric Functions
    #######################################################################
    def arcsin(self):
        r"""
        This function returns the complex arcsine of the complex number
        `z`, `{\rm arcsin}(z)`. The branch cuts are on the
        real axis, less than -1 and greater than 1.

        EXAMPLES::

            sage: CDF(1,1).arcsin()
            0.666239432493 + 1.06127506191*I
        """
        return self._new_c(gsl_complex_arcsin(self._complex))

    def arccos(self):
        r"""
        This function returns the complex arccosine of the complex number
        `z`, `{\rm arccos}(z)`. The branch cuts are on the
        real axis, less than -1 and greater than 1.

        EXAMPLES::

            sage: CDF(1,1).arccos()
            0.904556894302 - 1.06127506191*I
        """
        return self._new_c(gsl_complex_arccos(self._complex))

    def arctan(self):
        r"""
        This function returns the complex arctangent of the complex number
        `z`, `{\rm arctan}(z)`. The branch cuts are on the
        imaginary axis, below `-i` and above `i`.

        EXAMPLES::

            sage: CDF(1,1).arctan()
            1.0172219679 + 0.402359478109*I
        """
        return self._new_c(gsl_complex_arctan(self._complex))

    def arccsc(self):
        r"""
        This function returns the complex arccosecant of the complex number
        `z`, `{\rm arccsc}(z) = {\rm arcsin}(1/z)`.

        EXAMPLES::

            sage: CDF(1,1).arccsc()
            0.452278447151 - 0.530637530953*I
        """
        return self._new_c(gsl_complex_arccsc(self._complex))

    def arccot(self):
        r"""
        This function returns the complex arccotangent of the complex
        number `z`, `{\rm arccot}(z) = {\rm arctan}(1/z).`

        EXAMPLES::

            sage: CDF(1,1).arccot()
            0.553574358897 - 0.402359478109*I
        """
        return self._new_c(gsl_complex_arccot(self._complex))


    #######################################################################
    # Complex Hyperbolic Functions
    #######################################################################
    def sinh(self):
        r"""
        This function returns the complex hyperbolic sine of the complex
        number `z`:

        .. MATH::

            \sinh(z) = \frac{e^z - e^{-z}}{2}.

        EXAMPLES::

            sage: CDF(1,1).sinh()
            0.634963914785 + 1.29845758142*I
        """
        return self._new_c(gsl_complex_sinh(self._complex))

    def cosh(self):
        r"""
        This function returns the complex hyperbolic cosine of the complex
        number `z`:

        .. MATH::

            \cosh(z) = \frac{e^z + e^{-z}}{2}.

        EXAMPLES::

            sage: CDF(1,1).cosh()
            0.833730025131 + 0.988897705763*I
        """
        return self._new_c(gsl_complex_cosh(self._complex))

    def tanh(self):
        r"""
        This function returns the complex hyperbolic tangent of the complex
        number `z`:

        .. MATH::

            \tanh(z) = \frac{\sinh(z)}{\cosh(z)}.

        EXAMPLES::

            sage: CDF(1,1).tanh()
            1.08392332734 + 0.27175258532*I
        """
        return self._new_c(gsl_complex_tanh(self._complex))


    def sech(self):
        r"""
        This function returns the complex hyperbolic secant of the complex
        number `z`:

        .. MATH::

            {\rm sech}(z) = \frac{1}{{\rm cosh}(z)}.

        EXAMPLES::

            sage: CDF(1,1).sech()
            0.498337030555 - 0.591083841721*I
        """
        return self._new_c(gsl_complex_sech(self._complex))

    def csch(self):
        r"""
        This function returns the complex hyperbolic cosecant of the
        complex number `z`:

        .. MATH::

            {\rm csch}(z) = \frac{1}{{\rm sinh}(z)}.

        EXAMPLES::

            sage: CDF(1,1).csch()
            0.303931001628 - 0.62151801717*I
        """
        return self._new_c(gsl_complex_csch(self._complex))

    def coth(self):
        r"""
        This function returns the complex hyperbolic cotangent of the
        complex number `z`:

        .. MATH::

            \coth(z) = \frac{1}{\tanh(z)}.

        EXAMPLES::

            sage: CDF(1,1).coth()
            0.868014142896 - 0.217621561854*I
        """
        return self._new_c(gsl_complex_coth(self._complex))

    #######################################################################
    # Inverse Complex Hyperbolic Functions
    #######################################################################
    def arcsinh(self):
        r"""
        This function returns the complex hyperbolic arcsine of the complex
        number `z`, `{\rm arcsinh}(z)`. The branch cuts are
        on the imaginary axis, below `-i` and above `i`.

        EXAMPLES::

            sage: CDF(1,1).arcsinh()
            1.06127506191 + 0.666239432493*I
        """
        return self._new_c(gsl_complex_arcsinh(self._complex))

    def arccosh(self):
        r"""
        This function returns the complex hyperbolic arccosine of the
        complex number `z`, `{\rm arccosh}(z)`. The branch
        cut is on the real axis, less than 1.

        EXAMPLES::

            sage: CDF(1,1).arccosh()
            1.06127506191 + 0.904556894302*I
        """
        return self._new_c(gsl_complex_arccosh(self._complex))

    def arctanh(self):
        r"""
        This function returns the complex hyperbolic arctangent of the
        complex number `z`, `{\rm arctanh} (z)`. The branch
        cuts are on the real axis, less than -1 and greater than 1.

        EXAMPLES::

            sage: CDF(1,1).arctanh()
            0.402359478109 + 1.0172219679*I
        """
        return self._new_c(gsl_complex_arctanh(self._complex))

    def arcsech(self):
        r"""
        This function returns the complex hyperbolic arcsecant of the
        complex number `z`, `{\rm arcsech}(z) = {\rm arccosh}(1/z)`.

        EXAMPLES::

            sage: CDF(1,1).arcsech()
            0.530637530953 - 1.11851787964*I
        """
        return self._new_c(gsl_complex_arcsech(self._complex))

    def arccsch(self):
        r"""
        This function returns the complex hyperbolic arccosecant of the
        complex number `z`, `{\rm arccsch}(z) = {\rm arcsin}(1/z)`.

        EXAMPLES::

            sage: CDF(1,1).arccsch()
            0.530637530953 - 0.452278447151*I
        """
        return self._new_c(gsl_complex_arccsch(self._complex))

    def arccoth(self):
        r"""
        This function returns the complex hyperbolic arccotangent of the
        complex number `z`, `{\rm arccoth}(z) = {\rm arctanh(1/z)}`.

        EXAMPLES::

            sage: CDF(1,1).arccoth()
            0.402359478109 - 0.553574358897*I
        """
        return self._new_c(gsl_complex_arccoth(self._complex))

    #######################################################################
    # Special Functions (from PARI)
    #######################################################################
    def eta(self, int omit_frac=0):
        r"""
        Return the value of the Dedekind `\eta` function on self.

        INPUT:

        -  ``self`` - element of the upper half plane (if not,
           raises a ValueError).

        -  ``omit_frac`` - (bool, default: ``False``), if ``True``,
           omit the `e^{\pi i z / 12}` factor.

        OUTPUT: a complex double number

        ALGORITHM: Uses the PARI C library.

        The `\eta` function is

        .. MATH::

            \eta(z) = e^{\pi i z / 12} \prod_{n=1}^{\infty} (1 - e^{2\pi inz})

        EXAMPLES:

        We compute a few values of :meth:`eta()`::

            sage: CDF(0,1).eta()
            0.768225422326
            sage: CDF(1,1).eta()
            0.742048775837 + 0.19883137023*I
            sage: CDF(25,1).eta()
            0.742048775837 + 0.19883137023*I

        :meth:`eta()` works even if the inputs are large::

            sage: CDF(0, 10^15).eta()
            0.0
            sage: CDF(10^15, 0.1).eta()  # abs tol 1e-10
            -0.115342592727 - 0.19977923088*I

        We compute a few values of :meth:`eta()`, but with the fractional power
        of `e` omitted::

            sage: CDF(0,1).eta(True)
            0.998129069926

        We compute :meth:`eta()` to low precision directly from the
        definition::

            sage: z = CDF(1,1); z.eta()
            0.742048775837 + 0.19883137023*I
            sage: i = CDF(0,1); pi = CDF(pi)
            sage: exp(pi * i * z / 12) * prod([1-exp(2*pi*i*n*z) for n in range(1,10)])
            0.742048775837 + 0.19883137023*I

        The optional argument allows us to omit the fractional part::

            sage: z.eta(omit_frac=True)  # abs tol 1e-12
            0.998129069926 - 8.12769318782e-22*I
            sage: pi = CDF(pi)
            sage: prod([1-exp(2*pi*i*n*z) for n in range(1,10)])  # abs tol 1e-12
            0.998129069926 + 4.59084695545e-19*I

        We illustrate what happens when `z` is not in the upper half plane::

            sage: z = CDF(1)
            sage: z.eta()
            Traceback (most recent call last):
            ...
            ValueError: value must be in the upper half plane

        You can also use functional notation::

            sage: z = CDF(1,1) ; eta(z)
            0.742048775837 + 0.19883137023*I
        """
        cdef GEN a, b, c, y, t

        if self._complex.dat[1] <= 0:
            raise ValueError, "value must be in the upper half plane"

        if self._complex.dat[1] > 100000 and not omit_frac:
            # To the precision of doubles for such large imaginary
            # part, the answer is automatically 0. If we don't do
            # this, PARI can easily underflow.
            return ComplexDoubleElement(0,0)

        cdef int flag = 0 if omit_frac else 1
        return self._new_from_gen(self._pari_().eta(flag))

    def agm(self, right, algorithm="optimal"):
        r"""
        Return the Arithmetic-Geometric Mean (AGM) of ``self`` and ``right``.

        INPUT:

        - ``right`` (complex) -- another complex number

        - ``algorithm`` (string, default ``"optimal"``) -- the algorithm to use
          (see below).

        OUTPUT:

        (complex) A value of the AGM of self and right.  Note that
        this is a multi-valued function, and the algorithm used
        affects the value returned, as follows:

        - ``'pari'``: Call the agm function from the pari library.

        - ``'optimal'``: Use the AGM sequence such that at each stage
          `(a,b)` is replaced by `(a_1,b_1)=((a+b)/2,\pm\sqrt{ab})`
          where the sign is chosen so that `|a_1-b_1| \leq |a_1+b_1|`, or
          equivalently `\Re(b_1/a_1) \geq 0`.  The resulting limit is
          maximal among all possible values.

        - ``'principal'``: Use the AGM sequence such that at each stage
          `(a,b)` is replaced by `(a_1,b_1)=((a+b)/2,\pm\sqrt{ab})`
          where the sign is chosen so that `\Re(b_1/a_1) \geq 0` (the
          so-called principal branch of the square root).

        EXAMPLES::

            sage: i = CDF(I)
            sage: (1+i).agm(2-i)
            1.62780548487 + 0.136827548397*I

        An example to show that the returned value depends on the algorithm
        parameter::

            sage: a = CDF(-0.95,-0.65)
            sage: b = CDF(0.683,0.747)
            sage: a.agm(b, algorithm='optimal')
            -0.371591652352 + 0.319894660207*I
            sage: a.agm(b, algorithm='principal')
            0.338175462986 - 0.0135326969565*I
            sage: a.agm(b, algorithm='pari')
            0.080689185076 + 0.239036532686*I

        Some degenerate cases::

            sage: CDF(0).agm(a)
            0.0
            sage: a.agm(0)
            0.0
            sage: a.agm(-a)
            0.0
        """
        cdef double complex a, b, a1, b1, r
        cdef double d, e, eps = 2.0**-51

        if algorithm == "pari":
            return self._new_from_gen(self._pari_().agm(right))

        if not isinstance(right, ComplexDoubleElement):
            right = CDF(right)

        a = extract_double_complex(self)
        b = extract_double_complex(<ComplexDoubleElement>right)

        if a == 0 or b == 0 or a == -b:
            return ComplexDoubleElement(0, 0)

        if algorithm=="optimal":
            while True:
                a1 = (a+b)/2
                b1 = csqrt(a*b)
                r = b1/a1
                d  = cabs(r-1)
                e  = cabs(r+1)
                if e < d:
                    b1=-b1
                    d = e
                if d < eps: return ComplexDoubleElement_from_doubles(a1.real, a1.imag)
                a, b = a1, b1

        elif algorithm=="principal":
            while True:
                a1 = (a+b)/2
                b1 = csqrt(a*b)
                if cabs((b1/a1)-1) < eps: return ComplexDoubleElement_from_doubles(a1.real, a1.imag)
                a, b = a1, b1

        else:
            raise ValueError, "agm algorithm must be one of 'pari', 'optimal', 'principal'"

    def dilog(self):
        r"""
        Returns the principal branch of the dilogarithm of `x`, i.e., analytic
        continuation of the power series

        .. MATH::

            \log_2(x) = \sum_{n \ge 1} x^n / n^2.

        EXAMPLES::

            sage: CDF(1,2).dilog()
            -0.0594747986738 + 2.07264797177*I
            sage: CDF(10000000,10000000).dilog()
            -134.411774491 + 38.793962999*I
        """
        return self._new_from_gen(self._pari_().dilog())

    def gamma(self):
        r"""
        Return the gamma function `\Gamma(z)` evaluated at ``self``, the
        complex number `z`.

        EXAMPLES::

            sage: CDF(5,0).gamma()
            24.0
            sage: CDF(1,1).gamma()
            0.498015668118 - 0.154949828302*I
            sage: CDF(0).gamma()
            Infinity
            sage: CDF(-1,0).gamma()
            Infinity
        """
        if self._complex.dat[1] == 0:
            if self._complex.dat[0] == 0:
                import infinity
                return infinity.unsigned_infinity
            try:
                from sage.rings.all import Integer, CC
                if Integer(self._complex.dat[0]) < 0:
                    return CC(self).gamma()
            except TypeError:
                pass
        return self._new_from_gen(self._pari_().gamma())

    def gamma_inc(self, t):
        r"""
        Return the incomplete gamma function evaluated at this complex number.

        EXAMPLES::

            sage: CDF(1,1).gamma_inc(CDF(2,3))
            0.00209691486365 - 0.0599819136554*I
            sage: CDF(1,1).gamma_inc(5)
            -0.00137813093622 + 0.00651982002312*I
            sage: CDF(2,0).gamma_inc(CDF(1,1))
            0.707092096346 - 0.42035364096*I
        """
        return self._new_from_gen(self._pari_().incgam(t))

    def zeta(self):
        """
        Return the Riemann zeta function evaluated at this complex number.

        EXAMPLES::

            sage: z = CDF(1, 1)
            sage: z.zeta()
            0.582158059752 - 0.926848564331*I
            sage: zeta(z)
            0.582158059752 - 0.926848564331*I
            sage: zeta(CDF(1))
            Infinity
        """
        if self._complex.dat[0] == 1 and self._complex.dat[1] == 0:
            import infinity
            return infinity.unsigned_infinity
        return self._new_from_gen(self._pari_().zeta())

    def algdep(self, long n):
        """
        Returns a polynomial of degree at most `n` which is
        approximately satisfied by this complex number. Note that the
        returned polynomial need not be irreducible, and indeed usually
        won't be if `z` is a good approximation to an algebraic
        number of degree less than `n`.

        ALGORITHM: Uses the PARI C-library algdep command.

        EXAMPLES::

            sage: z = (1/2)*(1 + RDF(sqrt(3)) *CDF.0); z
            0.5 + 0.866025403784*I
            sage: p = z.algdep(5); p
            x^3 + 1
            sage: p.factor()
            (x + 1) * (x^2 - x + 1)
            sage: abs(z^2 - z + 1) < 1e-14
            True

        ::

            sage: CDF(0,2).algdep(10)
            x^2 + 4
            sage: CDF(1,5).algdep(2)
            x^2 - 2*x + 26
        """
        pari_catch_sig_on()
        f = pari.new_gen(algdep0(self._gen(), n, 0))
        from polynomial.polynomial_ring_constructor import PolynomialRing
        from integer_ring import ZZ
        R = PolynomialRing(ZZ ,'x')
        return R(list(reversed(eval(str(f.Vec())))))


cdef class FloatToCDF(Morphism):
    """
    Fast morphism from anything with a ``__float__`` method to an RDF element.

    EXAMPLES::

        sage: f = CDF.coerce_map_from(ZZ); f
        Native morphism:
          From: Integer Ring
          To:   Complex Double Field
        sage: f(4)
        4.0
        sage: f = CDF.coerce_map_from(QQ); f
        Native morphism:
          From: Rational Field
          To:   Complex Double Field
        sage: f(1/2)
        0.5
        sage: f = CDF.coerce_map_from(int); f
        Native morphism:
          From: Set of Python objects of type 'int'
          To:   Complex Double Field
        sage: f(3r)
        3.0
        sage: f = CDF.coerce_map_from(float); f
        Native morphism:
          From: Set of Python objects of type 'float'
          To:   Complex Double Field
        sage: f(3.5)
        3.5
    """
    def __init__(self, R):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: f = CDF.coerce_map_from(ZZ); f
            Native morphism:
              From: Integer Ring
              To:   Complex Double Field
        """
        from sage.categories.homset import Hom
        if isinstance(R, type):
            from sage.structure.parent import Set_PythonType
            R = Set_PythonType(R)
        Morphism.__init__(self, Hom(R, CDF))

    cpdef Element _call_(self, x):
        """
        Create an :class:`ComplexDoubleElement`.

        EXAMPLES::

            sage: CDF((1,2)) # indirect doctest
            1.0 + 2.0*I
            sage: CDF('i') # indirect doctest
            1.0*I
            sage: CDF(2+i) # indirect doctest
            2.0 + 1.0*I
        """
        cdef ComplexDoubleElement z = <ComplexDoubleElement>PY_NEW(ComplexDoubleElement)
        z._complex = gsl_complex_rect(x, 0)
        return z

    def _repr_type(self):
        """
        Return string that describes the type of morphism.

        EXAMPLES::

            sage: sage.rings.complex_double.FloatToCDF(QQ)._repr_type()
            'Native'
        """
        return "Native"

#####################################################
# Create new ComplexDoubleElement from a
# gsl_complex C struct.
#####################################################

cdef GEN complex_gen(x):
    """
    Complex generator.

    INPUT:

    -  A Python object ``x``

    OUTPUT:

    -  A GEN of type t_COMPLEX, or raise a ``TypeError``.
    """
    cdef ComplexDoubleElement z
    global _CDF
    try:
        z = x
    except TypeError:
        z = _CDF(x)
    return z._gen()





#####################################################
# unique objects
#####################################################
cdef ComplexDoubleField_class _CDF
_CDF = ComplexDoubleField_class()
CDF = _CDF  # external interface
cdef ComplexDoubleElement I = ComplexDoubleElement(0,1)

def ComplexDoubleField():
    """
    Returns the field of double precision complex numbers.

    EXAMPLE::

        sage: ComplexDoubleField()
        Complex Double Field
        sage: ComplexDoubleField() is CDF
        True
    """
    return _CDF

from sage.misc.parser import Parser
cdef cdf_parser = Parser(float, float,  {"I" : _CDF.gen(), "i" : _CDF.gen()})

cdef extern from "math.h":
       int isfinite(double x)
       int isnan(double x)
       int isinf(double x)

cdef double_to_str(double x):
    """
    Convert a double to a string.
    """
    if isfinite(x):
        return str(x)
    if isnan(x):
        return "NaN"
    # C99 only guarantees that isinf() returns a nonzero value (actually: 1) if x is an
    # infinity (positive or negative). Modern Linux distros return -1 or +1 depending
    # on the sign of infinity, but that is not the case on OS X or Solaris
    if isinf(x) != 0 and x < 0:
        return '-infinity'
    elif isinf(x) != 0 and x > 0:
        return '+infinity'

cdef inline double complex extract_double_complex(ComplexDoubleElement x):
    """
    Return the value of ``x`` as a c99 complex double.
    """
    cdef double complex z
    z.real = x._complex.dat[0]
    z.imag = x._complex.dat[1]
    return z

cdef ComplexDoubleElement ComplexDoubleElement_from_doubles(double re, double im):
    """
    Create a new :class:`ComplexDoubleElement` with the specified real and
    imaginary parts.
    """
    cdef ComplexDoubleElement z = <ComplexDoubleElement>PY_NEW(ComplexDoubleElement)
    z._complex.dat[0] = re
    z._complex.dat[1] = im
    return z

#####
