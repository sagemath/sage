# distutils: extra_compile_args = -D_XPG6
# distutils: libraries = m
r"""
Double precision floating point complex numbers

Sage supports arithmetic using double-precision complex numbers. A
double-precision complex number is a complex number ``x + I*y`` with
`x`, `y` 64-bit (8 byte) floating point numbers (double precision).

The field :class:`ComplexDoubleField` implements the field
of all double-precision complex numbers. You can refer to this
field by the shorthand ``CDF``. Elements of this field are of type
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
    <class 'sage.rings.complex_double.ComplexDoubleElement'>
    sage: ComplexDoubleElement(sqrt(2), 3)                                              # needs sage.symbolic
    1.4142135623730951 + 3.0*I
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
functions and macros in :ref:`GSL <spkg_gsl>` (the GNU Scientific Library), and
should be very fast. Also, all standard complex trig functions,
log, exponents, etc., are implemented using GSL, and are also
robust and fast. Several other special functions, e.g. eta, gamma,
incomplete gamma, etc., are implemented using the :ref:`PARI <spkg_pari>`
C library.

AUTHORS:

- William Stein (2006-09): first version

- Travis Scrimshaw (2012-10-18): Added doctests to get full coverage

- Jeroen Demeyer (2013-02-27): fixed all PARI calls (:issue:`14082`)

- Vincent Klein (2017-11-15) : add __mpc__() to class ComplexDoubleElement.
  ComplexDoubleElement constructor support and gmpy2.mpc parameter.
"""

# ****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#       Copyright (C) 2013 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cpython.object cimport Py_NE

from sage.misc.randstate cimport randstate, current_randstate

from sage.libs.gsl.complex cimport *

cdef extern from "<complex.h>":
    double complex csqrt(double complex)
    double cabs(double complex)

import sage.rings.abc
cimport sage.rings.integer

from sage.structure.element cimport Element, FieldElement
from sage.structure.parent cimport Parent
from sage.structure.richcmp cimport rich_to_bool
from sage.categories.morphism cimport Morphism
from sage.structure.coerce cimport is_numpy_type

try:
    from cypari2.gen import Gen as pari_gen
    from sage.libs.pari.convert_sage_complex_double import pari_to_cdf

except ImportError:
    pari_gen = ()


new_gen_from_complex_double_element = None
complex_double_element_eta = None
complex_double_element_agm = None
complex_double_element_dilog = None
complex_double_element_gamma = None
complex_double_element_gamma_inc = None
complex_double_element_zeta = None

from sage.rings.complex_conversion cimport CCtoCDF

from sage.rings.real_double cimport RealDoubleElement, double_repr
from sage.rings.real_double import RDF
from sage.rings.integer_ring import ZZ
from sage.structure.richcmp cimport rich_to_bool

cimport gmpy2
gmpy2.import_gmpy2()


cdef class ComplexDoubleField_class(sage.rings.abc.ComplexDoubleField):
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

        EXAMPLES::

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
        Parent.__init__(self, self, names=('I',), normalize=False,
                        category=Fields().Infinite().Metric().Complete())
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
        Return whether or not this field is exact, which is always ``False``.

        EXAMPLES::

            sage: CDF.is_exact()
            False
        """
        return False

    def __richcmp__(left, right, int op):
        """
        Rich comparison of ``left`` against ``right``.

        EXAMPLES::

            sage: CDF == CDF
            True
        """
        if left is right:
            return rich_to_bool(op, 0)
        if isinstance(right, ComplexDoubleField_class):
            return rich_to_bool(op, 0)
        return op == Py_NE

    def __hash__(self):
        """
        Return the hash for ``self``.

        This class is intended for use as a singleton so any instance
        of it should be equivalent from a hashing perspective.

        TESTS::

            sage: from sage.rings.complex_double import ComplexDoubleField_class
            sage: hash(CDF) == hash(ComplexDoubleField_class())
            True
        """
        return 561162115

    def characteristic(self):
        """
        Return the characteristic of the complex double field, which is 0.

        EXAMPLES::

            sage: CDF.characteristic()
            0
        """
        from sage.rings.integer import Integer
        return Integer(0)

    def random_element(self, double xmin=-1, double xmax=1, double ymin=-1, double ymax=1):
        """
        Return a random element of this complex double field with real and
        imaginary part bounded by ``xmin``, ``xmax``, ``ymin``, ``ymax``.

        EXAMPLES::

            sage: CDF.random_element().parent() is CDF
            True
            sage: re, im = CDF.random_element()
            sage: -1 <= re <= 1, -1 <= im <= 1
            (True, True)
            sage: re, im = CDF.random_element(-10,10,-10,10)
            sage: -10 <= re <= 10, -10 <= im <= 10
            (True, True)
            sage: re, im = CDF.random_element(-10^20,10^20,-2,2)
            sage: -10^20 <= re <= 10^20, -2 <= im <= 2
            (True, True)
        """
        cdef randstate rstate = current_randstate()
        global _CDF
        cdef ComplexDoubleElement z
        cdef double imag = (ymax-ymin)*rstate.c_rand_double() + ymin
        cdef double real = (xmax-xmin)*rstate.c_rand_double() + xmin
        z = ComplexDoubleElement.__new__(ComplexDoubleElement)
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

        OUTPUT: string

        TESTS::

            sage: print(CDF._latex_())
            \Bold{C}
        """
        return r"\Bold{C}"

    def __call__(self, x=None, im=None):
        """
        Create a complex double using ``x`` and optionally an imaginary part
        ``im``.

        EXAMPLES::

            sage: CDF(0,1) # indirect doctest
            1.0*I
            sage: CDF(2/3) # indirect doctest
            0.6666666666666666
            sage: CDF(5) # indirect doctest
            5.0
            sage: CDF('i') # indirect doctest
            1.0*I
            sage: CDF(complex(2,-3)) # indirect doctest
            2.0 - 3.0*I
            sage: CDF(4.5) # indirect doctest
            4.5
            sage: CDF(1+I)  # indirect doctest                                          # needs sage.symbolic
            1.0 + 1.0*I
            sage: CDF(pari(1))                                                          # needs sage.libs.pari
            1.0
            sage: CDF(pari("I"))                                                        # needs sage.libs.pari
            1.0*I
            sage: CDF(pari("x^2 + x + 1").polroots()[0])                                # needs sage.libs.pari
            -0.5 - 0.8660254037844386*I
            sage: from gmpy2 import mpc
            sage: CDF(mpc('2.0+1.0j'))
            2.0 + 1.0*I

        A :exc:`TypeError` is raised if the coercion doesn't make sense::

            sage: CDF(QQ['x'].0)
            Traceback (most recent call last):
            ...
            TypeError: cannot convert nonconstant polynomial

        One can convert back and forth between double precision complex
        numbers and higher-precision ones, though of course there may be
        loss of precision::

            sage: # needs sage.rings.real_mpfr
            sage: a = ComplexField(200)(-2).sqrt(); a
            1.4142135623730950488016887242096980785696718753769480731767*I
            sage: b = CDF(a); b
            1.4142135623730951*I
            sage: a.parent()(b)
            1.4142135623730951454746218587388284504413604736328125000000*I
            sage: a.parent()(b) == b
            True
            sage: b == CC(a)
            True

        TESTS:

        Check that :issue:`31836` is fixed::

            sage: a = CDF(); a
            0.0
            sage: a.parent()
            Complex Double Field
        """
        # We implement __call__ to gracefully accept the second argument.
        if x is None:
            return self.zero()
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
        from sage.rings import complex_mpfr
        if isinstance(x, ComplexDoubleElement):
            return x
        elif isinstance(x, tuple):
            return ComplexDoubleElement(x[0], x[1])
        elif isinstance(x, (float, int)):
            return ComplexDoubleElement(x, 0)
        elif isinstance(x, complex):
            return ComplexDoubleElement(x.real, x.imag)
        elif isinstance(x, complex_mpfr.ComplexNumber):
            return ComplexDoubleElement(x.real(), x.imag())
        elif isinstance(x, pari_gen):
            return pari_to_cdf(x)
        elif type(x) is gmpy2.mpc:
            return ComplexDoubleElement((<gmpy2.mpc>x).real, (<gmpy2.mpc>x).imag)
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
        it is defined, otherwise raise a :exc:`TypeError`.

        The rings that canonically coerce to the complex double field are:

        - the complex double field itself
        - anything that canonically coerces to real double field.
        - mathematical constants
        - the 53-bit mpfr complex field

        EXAMPLES::

            sage: CDF.coerce(5) # indirect doctest
            5.0
            sage: CDF.coerce(RDF(3.4))
            3.4

        Thus the sum of a CDF and a symbolic object is symbolic::

            sage: a = pi + CDF.0; a                                                     # needs sage.symbolic
            pi + 1.0*I
            sage: parent(a)                                                             # needs sage.symbolic
            Symbolic Ring

        TESTS::

            sage: # needs sage.rings.real_mpfr
            sage: CDF(1) + RR(1)
            2.0
            sage: CDF.0 - CC(1) - int(1) - RR(1) - QQbar(1)                             # needs sage.rings.number_field
            -4.0 + 1.0*I
            sage: CDF.has_coerce_map_from(ComplexField(20))
            False
            sage: CDF.has_coerce_map_from(complex)
            True
        """
        if S is int or S is float:
            return FloatToCDF(S)
        from sage.rings.rational_field import QQ
        from sage.rings.real_lazy import RLF
        from sage.rings.real_mpfr import RR
        from sage.rings.cc import CC

        if S is ZZ or S is QQ or S is RDF or S is RLF:
            return FloatToCDF(S)
        if isinstance(S, sage.rings.abc.RealField):
            if S.prec() >= 53:
                return FloatToCDF(S)
            else:
                return None
        elif is_numpy_type(S):
            import numpy
            if issubclass(S, numpy.integer) or issubclass(S, numpy.floating):
                return FloatToCDF(S)
            elif issubclass(S, numpy.complexfloating):
                return ComplexToCDF(S)
            else:
                return None
        elif RR.has_coerce_map_from(S):
            return FloatToCDF(RR) * RR._internal_coerce_map_from(S)
        elif isinstance(S, sage.rings.abc.ComplexField) and S.prec() >= 53:
            return CCtoCDF(S, self)
        elif CC.has_coerce_map_from(S):
            return CCtoCDF(CC, self) * CC._internal_coerce_map_from(S)

    def _magma_init_(self, magma):
        r"""
        Return a string representation of ``self`` in the Magma language.

        EXAMPLES::

            sage: CDF._magma_init_(magma)                       # optional - magma
            'ComplexField(53 : Bits := true)'
            sage: magma(CDF)                                    # optional - magma
            Complex field of precision 15
            sage: floor(RR(log(2**53, 10)))                                             # needs sage.symbolic
            15
            sage: magma(CDF).sage()                             # optional - magma
            Complex Field with 53 bits of precision
        """
        return "ComplexField(%s : Bits := true)" % self.prec()

    def _fricas_init_(self):
        r"""
        Return a string representation of ``self`` in the FriCAS language.

        EXAMPLES::

            sage: fricas(CDF)                     # indirect doctest, optional - fricas
            Complex(DoubleFloat)
        """
        return "Complex DoubleFloat"

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
        Return the complex field to the specified precision. As doubles
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
            from sage.rings.complex_mpfr import ComplexField
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
            raise ValueError("only 1 generator")
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
        Return the algebraic closure of ``self``, i.e., the complex double
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
        return RDF

    def pi(self):
        r"""
        Return `\pi` as a double precision complex number.

        EXAMPLES::

            sage: CDF.pi()
            3.141592653589793
        """
        return self(3.1415926535897932384626433832)

    def construction(self):
        """
        Return the functorial construction of ``self``, namely, algebraic
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

        - ``n`` -- positive integer (default: 2)

        OUTPUT: a complex `n`-th root of unity

        EXAMPLES::

            sage: CDF.zeta(7)  # rel tol 1e-15
            0.6234898018587336 + 0.7818314824680298*I
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
        from sage.rings.integer import Integer
        try:
            n = Integer(n)
        except TypeError:
            raise ValueError("n must be a positive integer")

        if n<1:
            raise ValueError("n must be a positive integer")

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

    def _factor_univariate_polynomial(self, f):
        """
        Factor the univariate polynomial ``f``.

        INPUT:

        - ``f`` -- a univariate polynomial defined over the double precision
          complex numbers

        OUTPUT:

        - A factorization of ``f`` over the double precision complex numbers
          into a unit and monic irreducible factors

        .. NOTE::

            This is a helper method for
            :meth:`sage.rings.polynomial.polynomial_element.Polynomial.factor`.

        TESTS::

            sage: R.<x> = CDF[]
            sage: CDF._factor_univariate_polynomial(x)
            x
            sage: CDF._factor_univariate_polynomial(2*x)
            (2.0) * x
            sage: CDF._factor_univariate_polynomial(x^2)
            x^2
            sage: f = x^2 + 1
            sage: F = CDF._factor_univariate_polynomial(f)
            sage: [f(t[0][0]).abs() for t in F] # abs tol 1e-9
            [5.55111512313e-17, 6.66133814775e-16]
            sage: f = (x^2 + 2*R(I))^3
            sage: F = f.factor()
            sage: [f(t[0][0]).abs() for t in F] # abs tol 1e-9
            [1.979365054e-14, 1.97936298566e-14, 1.97936990747e-14, 3.6812407475e-14, 3.65211563729e-14, 3.65220890052e-14]
        """
        unit = f.leading_coefficient()
        f *= ~unit
        roots = f.roots()
        from sage.misc.flatten import flatten
        roots = flatten([[r]*m for r, m in roots])
        from sage.structure.factorization import Factorization
        x = f.parent().gen()
        return Factorization([(x - a, 1) for a in roots], unit)


cdef ComplexDoubleElement new_ComplexDoubleElement():
    """
    Create a new (empty) :class:`ComplexDoubleElement`.
    """
    cdef ComplexDoubleElement z
    z = ComplexDoubleElement.__new__(ComplexDoubleElement)
    return z


def is_ComplexDoubleElement(x):
    """
    Return ``True`` if ``x`` is a :class:`ComplexDoubleElement`.

    EXAMPLES::

        sage: from sage.rings.complex_double import is_ComplexDoubleElement
        sage: is_ComplexDoubleElement(0)
        doctest:warning...
        DeprecationWarning: The function is_ComplexDoubleElement is deprecated;
        use 'isinstance(..., ComplexDoubleElement)' instead.
        See https://github.com/sagemath/sage/issues/38128 for details.
        False
        sage: is_ComplexDoubleElement(CDF(0))
        True
    """
    from sage.misc.superseded import deprecation_cython
    deprecation_cython(38128,
                       "The function is_ComplexDoubleElement is deprecated; "
                       "use 'isinstance(..., ComplexDoubleElement)' instead.")
    return isinstance(x, ComplexDoubleElement)


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
        Construct an element of a complex double field with specified real
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
                (GSL_REAL(self._complex), GSL_IMAG(self._complex)))

    cdef ComplexDoubleElement _new_c(self, gsl_complex x):
        """
        C-level code for creating a :class:`ComplexDoubleElement` from a
        ``gsl_complex``.
        """
        cdef ComplexDoubleElement z = <ComplexDoubleElement>ComplexDoubleElement.__new__(ComplexDoubleElement)
        z._complex = x
        return z

    def __hash__(self):
        """
        Return the hash of ``self``, which coincides with the python ``float``
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

    cpdef _richcmp_(left, right, int op):
        """
        We order the complex numbers in dictionary order by real parts then
        imaginary parts.

        This order, of course, does not respect the field structure, though
        it agrees with the usual order on the real numbers.

        EXAMPLES::

            sage: CDF(1.2) > CDF(i)                                                     # needs sage.symbolic
            True
            sage: CDF(1) < CDF(2)
            True
            sage: CDF(1 + i) > CDF(-1 - i)                                              # needs sage.symbolic
            True

        ::

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
        if GSL_REAL(left._complex) < GSL_REAL((<ComplexDoubleElement>right)._complex):
            return rich_to_bool(op, -1)
        if GSL_REAL(left._complex) > GSL_REAL((<ComplexDoubleElement>right)._complex):
            return rich_to_bool(op, 1)
        if GSL_IMAG(left._complex) < GSL_IMAG((<ComplexDoubleElement>right)._complex):
            return rich_to_bool(op, -1)
        if GSL_IMAG(left._complex) > GSL_IMAG((<ComplexDoubleElement>right)._complex):
            return rich_to_bool(op, 1)
        return rich_to_bool(op, 0)

    def __getitem__(self, n):
        """
        Return the real or imaginary part of ``self``.

        INPUT:

        - ``n`` -- integer (either 0 or 1)

        Raises an :exc:`IndexError` if ``n`` is not 0 or 1.

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
        if n == 0:
            return GSL_REAL(self._complex)
        if n == 1:
            return GSL_IMAG(self._complex)
        raise IndexError("index n must be 0 or 1")

    def _magma_init_(self, magma):
        r"""
        Return the magma representation of ``self``.

        EXAMPLES::

            sage: # optional - magma
            sage: CDF((1.2, 0.3))._magma_init_(magma)
            'ComplexField(53 : Bits := true)![1.2, 0.3]'
            sage: magma(CDF(1.2, 0.3))                  # indirect doctest
            1.20000000000000 + 0.300000000000000*$.1
            sage: s = magma(CDF(1.2, 0.3)).sage(); s    # indirect doctest
            1.20000000000000 + 0.300000000000000*I
            sage: s.parent()
            Complex Field with 53 bits of precision
        """
        return "%s![%s, %s]" % (self.parent()._magma_init_(magma), self.real(), self.imag())

    def prec(self):
        """
        Return the precision of this number (to be more similar to
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
            TypeError: can...t convert complex to int; use int(abs(z))
            sage: int(abs(CDF(1,1)))
            1
        """
        raise TypeError("can't convert complex to int; use int(abs(z))")

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
            TypeError: unable to convert 2.0 + 1.0*I to float; use abs() or real_part() as desired
            sage: a.__float__()
            Traceback (most recent call last):
            ...
            TypeError: unable to convert 2.0 + 1.0*I to float; use abs() or real_part() as desired
            sage: float(abs(CDF(1,1)))
            1.4142135623730951
        """
        if GSL_IMAG(self._complex):
            raise TypeError(f"unable to convert {self} to float; use abs() or real_part() as desired")
        return GSL_REAL(self._complex)

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
        return complex(GSL_REAL(self._complex), GSL_IMAG(self._complex))

    def _interface_init_(self, I=None):
        """
        Return ``self`` formatted as a string, suitable as input to another
        computer algebra system. (This is the default function used for
        exporting to other computer algebra systems.)

        EXAMPLES::

            sage: s1 = CDF(exp(I)); s1                                                  # needs sage.symbolic
            0.5403023058681398 + 0.8414709848078965*I
            sage: s1._interface_init_()                                                 # needs sage.symbolic
            '0.54030230586813977 + 0.84147098480789650*I'
            sage: s1 == CDF(gp(s1))                                                     # needs sage.symbolic
            True
        """
        # Sending to another computer algebra system is slow anyway, right?
        from sage.rings.complex_mpfr import ComplexField
        return ComplexField()(self)._interface_init_(I)

    def _mathematica_init_(self):
        """
        TESTS:

        Check that :issue:`28814` is fixed::

            sage: mathematica(CDF(1e-25, 1e25))  # optional - mathematica
            1.*^-25 + 1.*^25*I
        """
        from sage.rings.complex_mpfr import ComplexField
        return ComplexField()(self)._mathematica_init_()

    def _maxima_init_(self, I=None):
        """
        Return a string representation of this complex number in the syntax of
        Maxima. That is, use ``%i`` to represent the complex unit.

        EXAMPLES::

            sage: CDF.0._maxima_init_()
            '1.0000000000000000*%i'
            sage: CDF(.5 + I)._maxima_init_()                                           # needs sage.symbolic
            '0.50000000000000000 + 1.0000000000000000*%i'
        """
        from sage.rings.complex_mpfr import ComplexField
        return ComplexField()(self)._maxima_init_(I)

    def _sympy_(self):
        """
        Convert this complex number to Sympy.

        EXAMPLES::

            sage: CDF(1, 0)._sympy_()                                                   # needs sympy
            1.00000000000000
            sage: CDF(1/3, 1)._sympy_()                                                 # needs sympy
            0.333333333333333 + 1.0*I
            sage: type(_)                                                               # needs sympy
            <class 'sympy.core.add.Add'>
        """
        x = GSL_REAL(self._complex)
        y = GSL_IMAG(self._complex)
        import sympy
        return sympy.Float(x) + sympy.Float(y) * sympy.I

    def _repr_(self):
        """
        Return the string representation of ``self``.

        EXAMPLES::

            sage: CDF(0, 2/3)
            0.6666666666666666*I
            sage: a = CDF(2,-3); a # indirect doctest
            2.0 - 3.0*I
            sage: a^2   # abs tol 4e-15
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
        x = GSL_REAL(self._complex)
        y = GSL_IMAG(self._complex)
        if x == 0:
            if y == 0:
                # Not sure what to do with the signs of the real and
                # imaginary zeros, let's not print any sign.
                return "0.0"
            s = ''
        else:
            s = double_repr(x)
            if y == 0:
                return s
            elif y < 0:
                s += " - "
                y = -y
            else:
                s += " + "

        return s + double_repr(y) + "*I"

    def __format__(self, format_spec):
        """
        Return a formatted string representation of this complex number.

        INPUT:

        - ``format_spec`` -- string; a floating point format specificier as
          defined by :python:`the format specification mini-language
          <library/string.html#formatspec>` in Python

        EXAMPLES::

            sage: format(CDF(32/3, 0), ' .4f')
            ' 10.6667 + 0.0000*I'
            sage: format(CDF(-2/3, -2/3), '.4e')
            '-6.6667e-01 - 6.6667e-01*I'
            sage: format(CDF(3, 0), '.4g')
            '3 + 0*I'
            sage: format(CDF(3, 0), '#.4g')
            '3.000 + 0.000*I'

        If the representation type character is absent, the output matches the
        string representation of the complex number. This has the effect that
        real and imaginary part are only shown if they are not zero::

            sage: format(CDF(0, 2/3), '.4')
            '0.6667*I'
            sage: format(CDF(2, 0), '.4')
            '2.0'
            sage: format(CDF(0, 0), '+#.4')
            '+0.000'
        """
        from sage.rings.complex_mpfr import _format_complex_number
        return _format_complex_number(GSL_REAL(self._complex),
                                                     GSL_IMAG(self._complex),
                                                     format_spec)

    def _latex_(self):
        """
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: CDF(1, 2)._latex_()
            '1.0 + 2.0i'
            sage: z = CDF(1,2)^100
            sage: z._latex_()
            '-6.443164690985956 \\times 10^{34} - 6.113241307762409 \\times 10^{34}i'
        """
        import re
        s = str(self).replace('*I', 'i')
        return re.sub(r"e\+?(-?\d+)", r" \\times 10^{\1}", s)

    def __pari__(self):
        """
        Return PARI version of ``self``, as ``t_COMPLEX`` or ``t_REAL``.

        EXAMPLES::

            sage: # needs sage.libs.pari
            sage: CDF(1,2).__pari__()
            1.00000000000000 + 2.00000000000000*I
            sage: pari(CDF(1,2))
            1.00000000000000 + 2.00000000000000*I
            sage: pari(CDF(2.0))
            2.00000000000000
            sage: pari(CDF(I))
            1.00000000000000*I
        """
        global new_gen_from_complex_double_element
        if new_gen_from_complex_double_element is None:
            from sage.libs.pari.convert_sage_complex_double import new_gen_from_complex_double_element
        return new_gen_from_complex_double_element(self)

    def __mpc__(self):
        """
        Convert Sage ``ComplexDoubleElement`` to gmpy2 ``mpc``.

        EXAMPLES::

            sage: c = CDF(2,1)
            sage: c.__mpc__()
            mpc('2.0+1.0j')
            sage: from gmpy2 import mpc
            sage: mpc(c)
            mpc('2.0+1.0j')
        """
        return gmpy2.mpc(GSL_REAL(self._complex), GSL_IMAG(self._complex))

    #######################################################################
    # Arithmetic
    #######################################################################

    cpdef _add_(self, right):
        """
        Add ``self`` and ``right``.

        EXAMPLES::

            sage: CDF(2,-3)._add_(CDF(1,-2))
            3.0 - 5.0*I
        """
        return self._new_c(gsl_complex_add(self._complex,
                                           (<ComplexDoubleElement>right)._complex))

    cpdef _sub_(self, right):
        """
        Subtract ``self`` and ``right``.

        EXAMPLES::

            sage: CDF(2,-3)._sub_(CDF(1,-2))
            1.0 - 1.0*I
        """
        return self._new_c(gsl_complex_sub(self._complex,
                                (<ComplexDoubleElement>right)._complex))

    cpdef _mul_(self, right):
        """
        Multiply ``self`` and ``right``.

        EXAMPLES::

            sage: CDF(2,-3)._mul_(CDF(1,-2))
            -4.0 - 7.0*I
        """
        return self._new_c(gsl_complex_mul(self._complex,
                       (<ComplexDoubleElement>right)._complex))

    cpdef _div_(self, right):
        """
        Divide ``self`` by ``right``.

        EXAMPLES::

            sage: CDF(2,-3)._div_(CDF(1,-2))  # rel tol 1e-15
            1.5999999999999999 + 0.19999999999999998*I
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
            0.39999999999999997 - 0.19999999999999998*I
            sage: 1/CDF(2,1)
            0.39999999999999997 - 0.19999999999999998*I

        The inverse of 0 is ``NaN`` (it doesn't raise an exception)::

            sage: ~(0*CDF(0,1))
            NaN + NaN*I
        """
        return self._new_c(gsl_complex_inverse(self._complex))

    cpdef _neg_(self):
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
            1.5707963267948966
            sage: CDF(0,-1).arg()
            -1.5707963267948966
            sage: CDF(-1,0).arg()
            3.141592653589793
        """
        return RealDoubleElement(gsl_complex_arg(self._complex))

    def __abs__(self):
        """
        This function returns the magnitude of the complex number `z`, `|z|`.

        EXAMPLES::

            sage: abs(CDF(1,2)) # indirect doctest
            2.23606797749979
            sage: abs(CDF(1,0)) # indirect doctest
            1.0
            sage: abs(CDF(-2,3))
            3.605551275463989
        """
        return RealDoubleElement(gsl_complex_abs(self._complex))

    def abs(self):
        """
        This function returns the magnitude `|z|` of the complex number `z`.

        .. SEEALSO::

            - :meth:`norm`

        EXAMPLES::

            sage: CDF(2,3).abs()
            3.605551275463989
        """
        return RealDoubleElement(gsl_complex_abs(self._complex))

    def argument(self):
        r"""
        This function returns the argument of the ``self``, the complex number
        `z`, in the interval `-\pi < arg(z) \leq \pi`.

        EXAMPLES::

            sage: CDF(6).argument()
            0.0
            sage: CDF(i).argument()                                                     # needs sage.symbolic
            1.5707963267948966
            sage: CDF(-1).argument()
            3.141592653589793
            sage: CDF(-1 - 0.000001*i).argument()                                       # needs sage.symbolic
            -3.1415916535897934
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

            - :meth:`sage.rings.complex_mpfr.ComplexNumber.norm`

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
            0.09942542937258267
            sage: log(abs(CDF(1.1,0.1)))
            0.09942542937258259

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
        return RealDoubleElement(GSL_REAL(self._complex))

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
        return RealDoubleElement(GSL_IMAG(self._complex))

    imag_part = imag

    #######################################################################
    # Elementary Complex Functions
    #######################################################################
    def sqrt(self, all=False, **kwds):
        r"""
        The square root function.

        INPUT:

        - ``all`` -- boolean (default: ``False``); if ``True``, return a
          list of all square roots

        If all is ``False``, the branch cut is the negative real axis. The
        result always lies in the right half of the complex plane.

        EXAMPLES:

        We compute several square roots::

            sage: a = CDF(2,3)
            sage: b = a.sqrt(); b  # rel tol 1e-15
            1.6741492280355401 + 0.8959774761298381*I
            sage: b^2  # rel tol 1e-15
            2.0 + 3.0*I
            sage: a^(1/2)   # abs tol 1e-16
            1.6741492280355401 + 0.895977476129838*I

        We compute the square root of -1::

            sage: a = CDF(-1)
            sage: a.sqrt()
            1.0*I

        We compute all square roots::

            sage: CDF(-2).sqrt(all=True)
            [1.4142135623730951*I, -1.4142135623730951*I]
            sage: CDF(0).sqrt(all=True)
            [0.0]
        """
        z = self._new_c(gsl_complex_sqrt(self._complex))
        if all:
            return [z] if z.is_zero() else [z, -z]
        return z

    def nth_root(self, n, all=False):
        """
        The ``n``-th root function.

        INPUT:

        - ``all`` -- boolean (default: ``False``); if ``True``, return a
          list of all ``n``-th roots

        EXAMPLES::

            sage: a = CDF(125)
            sage: a.nth_root(3)
            5.000000000000001
            sage: a = CDF(10, 2)
            sage: [r^5 for r in a.nth_root(5, all=True)]  # rel tol 1e-14
            [9.999999999999998 + 2.0*I, 9.999999999999993 + 2.000000000000002*I, 9.999999999999996 + 1.9999999999999907*I, 9.999999999999993 + 2.0000000000000004*I, 9.999999999999998 + 1.9999999999999802*I]
            sage: abs(sum(a.nth_root(111, all=True)))  # rel tol 0.1
            1.1057313523818259e-13
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

    def is_integer(self):
        """
        Return ``True`` if this number is a integer.

        EXAMPLES::

            sage: CDF(0.5).is_integer()
            False
            sage: CDF(I).is_integer()                                                   # needs sage.symbolic
            False
            sage: CDF(2).is_integer()
            True
        """
        return (self.real() in ZZ) and (self.imag()==0)

    def is_positive_infinity(self):
        r"""
        Check if ``self`` is `+\infty`.

        EXAMPLES::

            sage: CDF(1, 2).is_positive_infinity()
            False
            sage: CDF(oo, 0).is_positive_infinity()
            True
            sage: CDF(0, oo).is_positive_infinity()
            False
        """
        return self.real().is_positive_infinity() and self.imag().is_zero()

    def is_negative_infinity(self):
        r"""
        Check if ``self`` is `-\infty`.

        EXAMPLES::

            sage: CDF(1, 2).is_negative_infinity()
            False
            sage: CDF(-oo, 0).is_negative_infinity()
            True
            sage: CDF(0, -oo).is_negative_infinity()
            False
        """
        return self.real().is_negative_infinity() and self.imag().is_zero()

    def is_infinity(self):
        r"""
        Check if ``self`` is `\infty`.

        EXAMPLES::

            sage: CDF(1, 2).is_infinity()
            False
            sage: CDF(0, oo).is_infinity()
            True
        """
        return self.real().is_infinity() or self.imag().is_infinity()

    def is_NaN(self):
        r"""
        Check if ``self`` is not-a-number.

        EXAMPLES::

            sage: CDF(1, 2).is_NaN()
            False
            sage: CDF(NaN).is_NaN()                                                     # needs sage.symbolic
            True
            sage: (1/CDF(0, 0)).is_NaN()
            True
        """
        return self.real().is_NaN() or self.imag().is_NaN()

    cpdef _pow_(self, other):
        r"""
        The complex number ``self`` raised to the power ``other``.

        This is computed using complex logarithms and exponentials
        as `z ^ a = \exp(a \log(z))`, unless ``self`` is real and
        ``other`` is an exact integer: in that case, only real
        arithmetic is used.

        EXAMPLES::

            sage: a = CDF(1,1)
            sage: b = CDF(2,3)
            sage: c = a ^ b; c  # rel tol 5e-16
            -0.163450932107355 + 0.09600498360894891*I
            sage: c ^ (1/b)     # rel tol 2e-16
            1.0 + 1.0*I

        We compute the cube root of `-1` then cube it and observe a
        rounding error::

            sage: a = CDF(-1)^(1/3); a
            0.5000000000000001 + 0.8660254037844386*I
            sage: a^3  # rel tol 1e-4
            -1.0 + 1.2246467991473532e-16*I

        We raise to symbolic powers::

            sage: var('x, n')                                                           # needs sage.symbolic
            (x, n)
            sage: CDF(1.2)^x                                                            # needs sage.symbolic
            1.2^x
            sage: CDF(1.2)^(x^n + n^x)                                                  # needs sage.symbolic
            1.2^(n^x + x^n)

        A real number powered to an exact integer always yields a real
        number::

            sage: CDF(-0.26319743704743886) ^ 1
            -0.26319743704743886
            sage: CDF(-0.26319743704743886) ^ 0
            1.0
            sage: CDF(-0.26319743704743886) ^ -1
            -3.7994290948196405
            sage: CDF(-0.26319743704743886) ^ 2
            0.06927289086834054
            sage: CDF(-0.26319743704743886) ^ -2
            14.435661446561994
            sage: CDF(-1) ^ (2^100)
            1.0

        Note the difference with ::

            sage: CDF(-0.26319743704743886) ^ QQ(2)
            0.06927289086834053 - 1.6966964813919428e-17*I
            sage: CDF(-0.26319743704743886) ^ RDF(2)
            0.06927289086834053 - 1.6966964813919428e-17*I

        TESTS::

            sage: (-2) ^ CDF(-0.5)  # abs tol 1e-16
            4.329780281177467e-17 - 0.7071067811865476*I
        """
        z = <ComplexDoubleElement>other
        return self._new_c(gsl_complex_pow(self._complex, z._complex))

    cdef _pow_long(self, long other):
        if other == 1:
            return self
        elif other == 0:
            res = gsl_complex_rect(1, 0)
        elif other == -1:
            res = gsl_complex_inverse(self._complex)
        elif not GSL_IMAG(self._complex):
            # If self is real, the result should be real too
            real = GSL_REAL(self._complex) ** <double>other
            res = gsl_complex_rect(real, 0)
        else:
            # General case
            res = gsl_complex_pow_real(self._complex, other)
        return self._new_c(res)

    cpdef _pow_int(self, other):
        if not GSL_IMAG(self._complex):
            # If self is real, the result should be real too
            real = GSL_REAL(self._complex) ** <double>other
            res = gsl_complex_rect(real, 0)
        else:
            # General case
            res = gsl_complex_pow_real(self._complex, other)
        return self._new_c(res)

    def exp(self):
        r"""
        This function returns the complex exponential of the complex number
        `z`, `\exp(z)`.

        EXAMPLES::

            sage: CDF(1,1).exp()  # abs tol 4e-16
            1.4686939399158851 + 2.2873552871788423*I

        We numerically verify a famous identity to the precision of a double::

            sage: z = CDF(0, 2*pi); z                                                   # needs sage.symbolic
            6.283185307179586*I
            sage: exp(z)  # rel tol 1e-4                                                # needs sage.symbolic
            1.0 - 2.4492935982947064e-16*I
        """
        return self._new_c(gsl_complex_exp(self._complex))

    def log(self, base=None):
        r"""
        This function returns the complex natural logarithm to the given
        base of the complex number `z`, `\log(z)`. The
        branch cut is the negative real axis.

        INPUT:

        - ``base`` -- (default: `e`) the base of the natural logarithm

        EXAMPLES::

            sage: CDF(1,1).log()
            0.34657359027997264 + 0.7853981633974483*I

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
            0.15051499783199057 + 0.3410940884604603*I
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

            sage: CDF(1,1).log_b(10)  # rel tol 1e-15
            0.15051499783199057 + 0.3410940884604603*I
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
            1.2984575814159773 + 0.6349639147847361*I
        """
        return self._new_c(gsl_complex_sin(self._complex))

    def cos(self):
        r"""
        This function returns the complex cosine of the complex number `z`:

        .. MATH::

            \cos(z) = \frac{e^{iz} + e^{-iz}}{2}

        EXAMPLES::

            sage: CDF(1,1).cos()   # abs tol 1e-16
            0.8337300251311491 - 0.9888977057628651*I
        """
        return self._new_c(gsl_complex_cos(self._complex))

    def tan(self):
        r"""
        This function returns the complex tangent of the complex number `z`:

        .. MATH::

            \tan(z) = \frac{\sin(z)}{\cos(z)}.

        EXAMPLES::

            sage: CDF(1,1).tan()
            0.27175258531951174 + 1.0839233273386946*I
        """
        return self._new_c(gsl_complex_tan(self._complex))

    def sec(self):
        r"""
        This function returns the complex secant of the complex number `z`:

        .. MATH::

            {\rm sec}(z) = \frac{1}{\cos(z)}.

        EXAMPLES::

            sage: CDF(1,1).sec()  # rel tol 1e-15
            0.4983370305551868 + 0.591083841721045*I
        """
        return self._new_c(gsl_complex_sec(self._complex))

    def csc(self):
        r"""
        This function returns the complex cosecant of the complex number `z`:

        .. MATH::

            \csc(z) = \frac{1}{\sin(z)}.

        EXAMPLES::

            sage: CDF(1,1).csc()  # rel tol 1e-15
            0.6215180171704284 - 0.30393100162842646*I
        """
        return self._new_c(gsl_complex_csc(self._complex))

    def cot(self):
        r"""
        This function returns the complex cotangent of the complex number `z`:

        .. MATH::

            \cot(z) = \frac{1}{\tan(z)}.

        EXAMPLES::

            sage: CDF(1,1).cot()  # rel tol 1e-15
            0.21762156185440268 - 0.8680141428959249*I
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
            0.6662394324925152 + 1.0612750619050357*I
        """
        return self._new_c(gsl_complex_arcsin(self._complex))

    def arccos(self):
        r"""
        This function returns the complex arccosine of the complex number
        `z`, `{\rm arccos}(z)`. The branch cuts are on the
        real axis, less than -1 and greater than 1.

        EXAMPLES::

            sage: CDF(1,1).arccos()
            0.9045568943023814 - 1.0612750619050357*I
        """
        return self._new_c(gsl_complex_arccos(self._complex))

    def arctan(self):
        r"""
        This function returns the complex arctangent of the complex number
        `z`, `{\rm arctan}(z)`. The branch cuts are on the
        imaginary axis, below `-i` and above `i`.

        EXAMPLES::

            sage: CDF(1,1).arctan()
            1.0172219678978514 + 0.4023594781085251*I
        """
        return self._new_c(gsl_complex_arctan(self._complex))

    def arccsc(self):
        r"""
        This function returns the complex arccosecant of the complex number
        `z`, `{\rm arccsc}(z) = {\rm arcsin}(1/z)`.

        EXAMPLES::

            sage: CDF(1,1).arccsc()  # rel tol 1e-15
            0.45227844715119064 - 0.5306375309525178*I
        """
        return self._new_c(gsl_complex_arccsc(self._complex))

    def arccot(self):
        r"""
        This function returns the complex arccotangent of the complex
        number `z`, `{\rm arccot}(z) = {\rm arctan}(1/z).`

        EXAMPLES::

            sage: CDF(1,1).arccot()  # rel tol 1e-15
            0.5535743588970452 - 0.4023594781085251*I
        """
        return self._new_c(gsl_complex_arccot(self._complex))

    def arcsec(self):
        r"""
        This function returns the complex arcsecant of the complex number
        `z`, `{\rm arcsec}(z) = {\rm arccos}(1/z)`.

        EXAMPLES::

            sage: CDF(1,1).arcsec()  # rel tol 1e-15
            1.118517879643706 + 0.5306375309525178*I
        """
        return self._new_c(gsl_complex_arcsec(self._complex))

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
            0.6349639147847361 + 1.2984575814159773*I
        """
        return self._new_c(gsl_complex_sinh(self._complex))

    def cosh(self):
        r"""
        This function returns the complex hyperbolic cosine of the complex
        number `z`:

        .. MATH::

            \cosh(z) = \frac{e^z + e^{-z}}{2}.

        EXAMPLES::

            sage: CDF(1,1).cosh()  # abs tol 1e-16
            0.8337300251311491 + 0.9888977057628651*I
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
            1.0839233273386946 + 0.27175258531951174*I
        """
        return self._new_c(gsl_complex_tanh(self._complex))

    def sech(self):
        r"""
        This function returns the complex hyperbolic secant of the complex
        number `z`:

        .. MATH::

            {\rm sech}(z) = \frac{1}{{\rm cosh}(z)}.

        EXAMPLES::

            sage: CDF(1,1).sech()  # rel tol 1e-15
            0.4983370305551868 - 0.591083841721045*I
        """
        return self._new_c(gsl_complex_sech(self._complex))

    def csch(self):
        r"""
        This function returns the complex hyperbolic cosecant of the
        complex number `z`:

        .. MATH::

            {\rm csch}(z) = \frac{1}{{\rm sinh}(z)}.

        EXAMPLES::

            sage: CDF(1,1).csch()  # rel tol 1e-15
            0.30393100162842646 - 0.6215180171704284*I
        """
        return self._new_c(gsl_complex_csch(self._complex))

    def coth(self):
        r"""
        This function returns the complex hyperbolic cotangent of the
        complex number `z`:

        .. MATH::

            \coth(z) = \frac{1}{\tanh(z)}.

        EXAMPLES::

            sage: CDF(1,1).coth()  # rel tol 1e-15
            0.8680141428959249 - 0.21762156185440268*I
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
            1.0612750619050357 + 0.6662394324925152*I
        """
        return self._new_c(gsl_complex_arcsinh(self._complex))

    def arccosh(self):
        r"""
        This function returns the complex hyperbolic arccosine of the
        complex number `z`, `{\rm arccosh}(z)`. The branch
        cut is on the real axis, less than 1.

        EXAMPLES::

            sage: CDF(1,1).arccosh()
            1.0612750619050357 + 0.9045568943023814*I
        """
        return self._new_c(gsl_complex_arccosh(self._complex))

    def arctanh(self):
        r"""
        This function returns the complex hyperbolic arctangent of the
        complex number `z`, `{\rm arctanh} (z)`. The branch
        cuts are on the real axis, less than -1 and greater than 1.

        EXAMPLES::

            sage: CDF(1,1).arctanh()
            0.4023594781085251 + 1.0172219678978514*I
        """
        return self._new_c(gsl_complex_arctanh(self._complex))

    def arcsech(self):
        r"""
        This function returns the complex hyperbolic arcsecant of the
        complex number `z`, `{\rm arcsech}(z) = {\rm arccosh}(1/z)`.

        EXAMPLES::

            sage: CDF(1,1).arcsech()  # rel tol 1e-15
            0.5306375309525176 - 1.118517879643706*I
        """
        return self._new_c(gsl_complex_arcsech(self._complex))

    def arccsch(self):
        r"""
        This function returns the complex hyperbolic arccosecant of the
        complex number `z`, `{\rm arccsch}(z) = {\rm arcsin}(1/z)`.

        EXAMPLES::

            sage: CDF(1,1).arccsch()  # rel tol 1e-15
            0.5306375309525178 - 0.45227844715119064*I
        """
        return self._new_c(gsl_complex_arccsch(self._complex))

    def arccoth(self):
        r"""
        This function returns the complex hyperbolic arccotangent of the
        complex number `z`, `{\rm arccoth}(z) = {\rm arctanh(1/z)}`.

        EXAMPLES::

            sage: CDF(1,1).arccoth()  # rel tol 1e-15
            0.4023594781085251 - 0.5535743588970452*I
        """
        return self._new_c(gsl_complex_arccoth(self._complex))

    #######################################################################
    # Special Functions (from PARI)
    #######################################################################
    def eta(self, int omit_frac=0):
        r"""
        Return the value of the Dedekind `\eta` function on ``self``.

        INPUT:

        - ``self`` -- element of the upper half plane (if not,
          raises a ValueError)

        - ``omit_frac`` -- boolean (default: ``False``); if ``True``,
          omit the `e^{\pi i z / 12}` factor

        OUTPUT: a complex double number

        ALGORITHM: Uses the PARI C library.

        The `\eta` function is

        .. MATH::

            \eta(z) = e^{\pi i z / 12} \prod_{n=1}^{\infty} (1 - e^{2\pi inz})

        EXAMPLES:

        We compute a few values of :meth:`eta()`::

            sage: CDF(0,1).eta()                                                        # needs sage.libs.pari
            0.7682254223260566
            sage: CDF(1,1).eta()                                                        # needs sage.libs.pari
            0.7420487758365647 + 0.1988313702299107*I
            sage: CDF(25,1).eta()                                                       # needs sage.libs.pari
            0.7420487758365647 + 0.1988313702299107*I

        :meth:`eta()` works even if the inputs are large::

            sage: CDF(0, 10^15).eta()
            0.0
            sage: CDF(10^15, 0.1).eta()  # abs tol 1e-10                                # needs sage.libs.pari
            -0.115342592727 - 0.19977923088*I

        We compute a few values of :meth:`eta()`, but with the fractional power
        of `e` omitted::

            sage: CDF(0,1).eta(True)                                                    # needs sage.libs.pari
            0.9981290699259585

        We compute :meth:`eta()` to low precision directly from the
        definition::

            sage: z = CDF(1,1); z.eta()                                                 # needs sage.libs.pari
            0.7420487758365647 + 0.1988313702299107*I
            sage: i = CDF(0,1); pi = CDF(pi)                                            # needs sage.symbolic
            sage: exp(pi * i * z / 12) * prod(1 - exp(2*pi*i*n*z)                       # needs sage.libs.pari sage.symbolic
            ....:                             for n in range(1, 10))
            0.7420487758365647 + 0.19883137022991068*I

        The optional argument allows us to omit the fractional part::

            sage: z.eta(omit_frac=True)                                                 # needs sage.libs.pari
            0.9981290699259585
            sage: pi = CDF(pi)                                                          # needs sage.symbolic
            sage: prod(1 - exp(2*pi*i*n*z) for n in range(1,10))  # abs tol 1e-12       # needs sage.libs.pari sage.symbolic
            0.998129069926 + 4.59084695545e-19*I

        We illustrate what happens when `z` is not in the upper half plane::

            sage: z = CDF(1)
            sage: z.eta()
            Traceback (most recent call last):
            ...
            ValueError: value must be in the upper half plane

        You can also use functional notation::

            sage: z = CDF(1,1)
            sage: eta(z)                                                                # needs sage.libs.pari
            0.7420487758365647 + 0.1988313702299107*I
        """
        if GSL_IMAG(self._complex) <= 0:
            raise ValueError("value must be in the upper half plane")

        if GSL_IMAG(self._complex) > 100000 and not omit_frac:
            # To the precision of doubles for such large imaginary
            # part, the answer is automatically 0. If we don't do
            # this, PARI can easily underflow.
            return ComplexDoubleElement(0,0)

        global complex_double_element_eta
        if complex_double_element_eta is None:
            from sage.libs.pari.convert_sage_complex_double import complex_double_element_eta

        cdef int flag = 0 if omit_frac else 1
        return complex_double_element_eta(self, flag)

    def agm(self, right, algorithm='optimal'):
        r"""
        Return the Arithmetic-Geometric Mean (AGM) of ``self`` and ``right``.

        INPUT:

        - ``right`` -- complex; another complex number

        - ``algorithm`` -- string (default: ``'optimal'``); the algorithm to
          use (see below)

        OUTPUT:

        (complex) A value of the AGM of ``self`` and ``right``.  Note that
        this is a multi-valued function, and the algorithm used
        affects the value returned, as follows:

        - ``'pari'``: Call the :pari:`agm` function from the pari library.

        - ``'optimal'``: Use the AGM sequence such that at each stage
          `(a,b)` is replaced by `(a_1,b_1)=((a+b)/2,\pm\sqrt{ab})`
          where the sign is chosen so that `|a_1-b_1| \leq |a_1+b_1|`, or
          equivalently `\Re(b_1/a_1) \geq 0`.  The resulting limit is
          maximal among all possible values.

        - ``'principal'``: Use the AGM sequence such that at each stage
          `(a,b)` is replaced by `(a_1,b_1)=((a+b)/2,\pm\sqrt{ab})`
          where the sign is chosen so that `\Re(b_1/a_1) \geq 0` (the
          so-called principal branch of the square root).

        See :wikipedia:`Arithmetic-geometric mean`

        EXAMPLES::

            sage: i = CDF(I)                                                            # needs sage.symbolic
            sage: (1+i).agm(2-i)  # rel tol 1e-15                                       # needs sage.symbolic
            1.6278054848727064 + 0.1368275483973686*I

        An example to show that the returned value depends on the algorithm
        parameter::

            sage: a = CDF(-0.95,-0.65)
            sage: b = CDF(0.683,0.747)
            sage: a.agm(b, algorithm='optimal')
            -0.3715916523517613 + 0.31989466020683*I
            sage: a.agm(b, algorithm='principal')  # rel tol 1e-15
            0.33817546298618006 - 0.013532696956540503*I
            sage: a.agm(b, algorithm='pari')
            -0.37159165235176134 + 0.31989466020683005*I

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
            global complex_double_element_agm
            if complex_double_element_agm is None:
                from sage.libs.pari.convert_sage_complex_double import complex_double_element_agm
            return complex_double_element_agm(self, right)

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
            raise ValueError("agm algorithm must be one of 'pari', 'optimal', 'principal'")

    def dilog(self):
        r"""
        Return the principal branch of the dilogarithm of `x`, i.e., analytic
        continuation of the power series

        .. MATH::

            \log_2(x) = \sum_{n \ge 1} x^n / n^2.

        EXAMPLES::

            sage: CDF(1,2).dilog()                                                      # needs sage.libs.pari
            -0.059474798673809476 + 2.0726479717747566*I
            sage: CDF(10000000,10000000).dilog()                                        # needs sage.libs.pari
            -134.411774490731 + 38.79396299904504*I
        """
        global complex_double_element_dilog
        if complex_double_element_dilog is None:
            from sage.libs.pari.convert_sage_complex_double import complex_double_element_dilog
        return complex_double_element_dilog(self)

    def gamma(self):
        r"""
        Return the gamma function `\Gamma(z)` evaluated at ``self``, the
        complex number `z`.

        EXAMPLES::

            sage: # needs sage.libs.pari
            sage: CDF(5,0).gamma()
            24.0
            sage: CDF(1,1).gamma()
            0.49801566811835607 - 0.15494982830181067*I
            sage: CDF(0).gamma()
            Infinity
            sage: CDF(-1,0).gamma()
            Infinity
        """
        if not GSL_IMAG(self._complex):
            if GSL_REAL(self._complex) == 0:
                from sage.rings.infinity import unsigned_infinity
                return unsigned_infinity
            try:
                from sage.rings.integer import Integer
                if Integer(GSL_REAL(self._complex)) < 0:
                    from sage.rings.cc import CC
                    return CC(self).gamma()
            except TypeError:
                pass
        global complex_double_element_gamma
        if complex_double_element_gamma is None:
            from sage.libs.pari.convert_sage_complex_double import complex_double_element_gamma
        return complex_double_element_gamma(self)

    def gamma_inc(self, t):
        r"""
        Return the incomplete gamma function evaluated at this complex number.

        EXAMPLES::

            sage: CDF(1,1).gamma_inc(CDF(2,3))                                          # needs sage.libs.pari
            0.0020969148636468277 - 0.059981913655449706*I
            sage: CDF(1,1).gamma_inc(5)                                                 # needs sage.libs.pari
            -0.001378130936215849 + 0.006519820023119819*I
            sage: CDF(2,0).gamma_inc(CDF(1,1))                                          # needs sage.libs.pari
            0.7070920963459381 - 0.4203536409598115*I
        """
        global complex_double_element_gamma_inc
        if complex_double_element_gamma_inc is None:
            from sage.libs.pari.convert_sage_complex_double import complex_double_element_gamma_inc
        return complex_double_element_gamma_inc(self, t)

    def zeta(self):
        """
        Return the Riemann zeta function evaluated at this complex number.

        EXAMPLES::

            sage: z = CDF(1, 1)
            sage: z.zeta()                                                              # needs sage.libs.pari
            0.5821580597520036 - 0.9268485643308071*I
            sage: zeta(z)                                                               # needs sage.libs.pari
            0.5821580597520036 - 0.9268485643308071*I
            sage: zeta(CDF(1))                                                          # needs sage.libs.pari
            Infinity
        """
        if GSL_REAL(self._complex) == 1 and GSL_IMAG(self._complex) == 0:
            from sage.rings.infinity import unsigned_infinity
            return unsigned_infinity

        global complex_double_element_zeta
        if complex_double_element_zeta is None:
            from sage.libs.pari.convert_sage_complex_double import complex_double_element_zeta
        return complex_double_element_zeta(self)

    def algdep(self, long n):
        """
        Return a polynomial of degree at most `n` which is
        approximately satisfied by this complex number. Note that the
        returned polynomial need not be irreducible, and indeed usually
        won't be if `z` is a good approximation to an algebraic
        number of degree less than `n`.

        ALGORITHM: Uses the PARI C-library algdep command.

        EXAMPLES::

            sage: z = (1/2)*(1 + RDF(sqrt(3)) * CDF.0); z   # abs tol 1e-16             # needs sage.symbolic
            0.5 + 0.8660254037844387*I
            sage: p = z.algdep(5); p                                                    # needs sage.libs.pari sage.symbolic
            x^2 - x + 1
            sage: abs(z^2 - z + 1) < 1e-14                                              # needs sage.symbolic
            True

        ::

            sage: CDF(0,2).algdep(10)                                                   # needs sage.libs.pari
            x^2 + 4
            sage: CDF(1,5).algdep(2)                                                    # needs sage.libs.pari
            x^2 - 2*x + 26
        """
        from sage.arith.misc import algdep
        return algdep(self, n)

cdef class FloatToCDF(Morphism):
    """
    Fast morphism from anything with a ``__float__`` method to a CDF element.

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
          From: Set of Python objects of class 'int'
          To:   Complex Double Field
        sage: f(3r)
        3.0
        sage: f = CDF.coerce_map_from(float); f
        Native morphism:
          From: Set of Python objects of class 'float'
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
            from sage.sets.pythonclass import Set_PythonType
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
            sage: CDF(2+i)  # indirect doctest                                          # needs sage.symbolic
            2.0 + 1.0*I
        """
        cdef ComplexDoubleElement z = <ComplexDoubleElement>ComplexDoubleElement.__new__(ComplexDoubleElement)
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


cdef class ComplexToCDF(Morphism):
    r"""
    Fast morphism for anything such that the elements have attributes ``.real``
    and ``.imag`` (e.g. numpy complex types).

    EXAMPLES::

        sage: # needs numpy
        sage: import numpy
        sage: f = CDF.coerce_map_from(numpy.complex128)
        sage: f(numpy.complex128(I))
        1.0*I
        sage: f(numpy.complex128(I)).parent()
        Complex Double Field
    """
    def __init__(self, R):
        from sage.categories.homset import Hom
        if isinstance(R, type):
            from sage.sets.pythonclass import Set_PythonType
            R = Set_PythonType(R)
        Morphism.__init__(self, Hom(R, CDF))

    cpdef Element _call_(self, x):
        """
        Create an :class:`ComplexDoubleElement`.

        EXAMPLES::

            sage: import numpy                                                          # needs numpy
            sage: CDF(numpy.complex128(I))    # indirect doctest                          # needs numpy
            1.0*I
        """
        cdef ComplexDoubleElement z = <ComplexDoubleElement>ComplexDoubleElement.__new__(ComplexDoubleElement)
        z._complex = gsl_complex_rect(x.real, x.imag)
        return z

    def _repr_type(self):
        """
        Return string that describes the type of morphism.

        EXAMPLES::

            sage: import numpy                                                          # needs numpy
            sage: f = sage.rings.complex_double.ComplexToCDF(numpy.complex128)          # needs numpy
            sage: f._repr_type()                                                        # needs numpy
            'Native'
        """
        return "Native"


#####################################################
# unique objects
#####################################################
cdef ComplexDoubleField_class _CDF
_CDF = ComplexDoubleField_class()
CDF = _CDF  # external interface
cdef ComplexDoubleElement I = ComplexDoubleElement(0,1)


def ComplexDoubleField():
    """
    Return the field of double precision complex numbers.

    EXAMPLES::

        sage: ComplexDoubleField()
        Complex Double Field
        sage: ComplexDoubleField() is CDF
        True
    """
    return _CDF


from sage.misc.parser import Parser
cdef cdf_parser = Parser(float, float,  {"I" : _CDF.gen(), "i" : _CDF.gen()})

cdef inline double complex extract_double_complex(ComplexDoubleElement x) noexcept:
    """
    Return the value of ``x`` as a c99 complex double.
    """
    cdef double complex z
    z.real = GSL_REAL(x._complex)
    z.imag = GSL_IMAG(x._complex)
    return z


cdef inline ComplexDoubleElement ComplexDoubleElement_from_doubles(double re, double im):
    """
    Create a new :class:`ComplexDoubleElement` with the specified real and
    imaginary parts.
    """
    z = <ComplexDoubleElement>ComplexDoubleElement.__new__(ComplexDoubleElement)
    GSL_SET_COMPLEX(&z._complex, re, im)
    return z


# Support Python's numbers abstract base class
import numbers
numbers.Complex.register(ComplexDoubleElement)
