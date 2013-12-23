r"""
Field of Arbitrary Precision Complex Numbers

AUTHORS:

- William Stein (2006-01-26): complete rewrite

- Niles Johnson (2010-08): :trac:`3893`: ``random_element()`` should pass on
  ``*args`` and ``**kwds``.

- Travis Scrimshaw (2012-10-18): Added documentation for full coverage.
"""

#################################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import complex_number
import complex_double
import field
import integer
import real_mpfr
import weakref
from sage.misc.sage_eval import sage_eval

from sage.structure.parent import Parent
from sage.structure.parent_gens import ParentWithGens

NumberField_quadratic = None
NumberFieldElement_quadratic = None
AlgebraicNumber_base = None
AlgebraicNumber = None
AlgebraicReal = None
AA = None
QQbar = None
SR = None
CDF = CLF = RLF = None

def late_import():
    """
    Import the objects/modules after build (when needed).

    TESTS::

        sage: sage.rings.complex_field.late_import()
    """
    global NumberField_quadratic
    global NumberFieldElement_quadratic
    global AlgebraicNumber_base
    global AlgebraicNumber
    global AlgebraicReal
    global AA, QQbar, SR
    global CLF, RLF, CDF
    if NumberFieldElement_quadratic is None:
        import sage.rings.number_field.number_field
        import sage.rings.number_field.number_field_element_quadratic as nfeq
        NumberField_quadratic = sage.rings.number_field.number_field.NumberField_quadratic
        NumberFieldElement_quadratic = nfeq.NumberFieldElement_quadratic
        import sage.rings.qqbar
        AlgebraicNumber_base = sage.rings.qqbar.AlgebraicNumber_base
        AlgebraicNumber = sage.rings.qqbar.AlgebraicNumber
        AlgebraicReal = sage.rings.qqbar.AlgebraicReal
        AA = sage.rings.qqbar.AA
        QQbar = sage.rings.qqbar.QQbar
        import sage.symbolic.ring
        SR = sage.symbolic.ring.SR
        from real_lazy import CLF, RLF
        from complex_double import CDF

def is_ComplexField(x):
    """
    Check if ``x`` is a :class:`complex field <ComplexField_class>`.

    EXAMPLES::

        sage: from sage.rings.complex_field import is_ComplexField as is_CF
        sage: is_CF(ComplexField())
        True
        sage: is_CF(ComplexField(12))
        True
        sage: is_CF(CC)
        True
    """
    return isinstance(x, ComplexField_class)

cache = {}
def ComplexField(prec=53, names=None):
    """
    Return the complex field with real and imaginary parts having prec
    *bits* of precision.

    EXAMPLES::

        sage: ComplexField()
        Complex Field with 53 bits of precision
        sage: ComplexField(100)
        Complex Field with 100 bits of precision
        sage: ComplexField(100).base_ring()
        Real Field with 100 bits of precision
        sage: i = ComplexField(200).gen()
        sage: i^2
        -1.0000000000000000000000000000000000000000000000000000000000
    """
    global cache
    if cache.has_key(prec):
        X = cache[prec]
        C = X()
        if not C is None:
            return C
    C = ComplexField_class(prec)
    cache[prec] = weakref.ref(C)
    return C


class ComplexField_class(field.Field):
    """
    An approximation to the field of complex numbers using floating
    point numbers with any specified precision. Answers derived from
    calculations in this approximation may differ from what they would
    be if those calculations were performed in the true field of
    complex numbers. This is due to the rounding errors inherent to
    finite precision calculations.

    EXAMPLES::

        sage: C = ComplexField(); C
        Complex Field with 53 bits of precision
        sage: Q = RationalField()
        sage: C(1/3)
        0.333333333333333
        sage: C(1/3, 2)
        0.333333333333333 + 2.00000000000000*I
        sage: C(RR.pi())
        3.14159265358979
        sage: C(RR.log2(), RR.pi())
        0.693147180559945 + 3.14159265358979*I

    We can also coerce rational numbers and integers into C, but
    coercing a polynomial will raise an exception::

        sage: Q = RationalField()
        sage: C(1/3)
        0.333333333333333
        sage: S = PolynomialRing(Q, 'x')
        sage: C(S.gen())
        Traceback (most recent call last):
        ...
        TypeError: unable to coerce to a ComplexNumber: <type 'sage.rings.polynomial.polynomial_rational_flint.Polynomial_rational_flint'>

    This illustrates precision::

        sage: CC = ComplexField(10); CC(1/3, 2/3)
        0.33 + 0.67*I
        sage: CC
        Complex Field with 10 bits of precision
        sage: CC = ComplexField(100); CC
        Complex Field with 100 bits of precision
        sage: z = CC(1/3, 2/3); z
        0.33333333333333333333333333333 + 0.66666666666666666666666666667*I

    We can load and save complex numbers and the complex field::

        sage: loads(z.dumps()) == z
        True
        sage: loads(CC.dumps()) == CC
        True
        sage: k = ComplexField(100)
        sage: loads(dumps(k)) == k
        True

    This illustrates basic properties of a complex field::

        sage: CC = ComplexField(200)
        sage: CC.is_field()
        True
        sage: CC.characteristic()
        0
        sage: CC.precision()
        200
        sage: CC.variable_name()
        'I'
        sage: CC == ComplexField(200)
        True
        sage: CC == ComplexField(53)
        False
        sage: CC == 1.1
        False
    """
    def __init__(self, prec=53):
        """
        Initialize ``self``.

        TESTS::

            sage: C = ComplexField(200)
            sage: C.category()
            Category of fields
            sage: TestSuite(C).run()
        """
        self._prec = int(prec)
        from sage.categories.fields import Fields
        ParentWithGens.__init__(self, self._real_field(), ('I',), False, category = Fields())
#        self._populate_coercion_lists_()
        self._populate_coercion_lists_(coerce_list=[complex_number.RRtoCC(self._real_field(), self)])

    def __reduce__(self):
        """
        For pickling.

        EXAMPLES::

            sage: loads(dumps(ComplexField())) == ComplexField()
            True
        """
        return ComplexField, (self._prec, )

    def is_exact(self):
        """
        Return whether or not this field is exact, which is always ``False``.

        EXAMPLES::

            sage: ComplexField().is_exact()
            False
        """
        return False

    def prec(self):
        """
        Return the precision of this complex field.

        EXAMPLES::

            sage: ComplexField().prec()
            53
            sage: ComplexField(15).prec()
            15
        """
        return self._prec

    def _magma_init_(self, magma):
        r"""
        Return a string representation of ``self`` in the Magma language.

        EXAMPLES::

            sage: ComplexField()._magma_init_(magma) # optional - magma
            'ComplexField(53 : Bits := true)'
            sage: magma(ComplexField(200)) # optional - magma
            Complex field of precision 60
            sage: 10^60 < 2^200 < 10^61
            True
            sage: s = magma(ComplexField(200)).sage(); s # optional - magma
            Complex Field with 200 bits of precision
            sage: 2^199 < 10^60 < 2^200
            True
            sage: s is ComplexField(200) # optional - magma
            True
        """
        return "ComplexField(%s : Bits := true)" % self.prec()

    precision = prec

    def to_prec(self, prec):
        """
        Returns the complex field to the specified precision.

        EXAMPLES::

            sage: CC.to_prec(10)
            Complex Field with 10 bits of precision
            sage: CC.to_prec(100)
            Complex Field with 100 bits of precision
        """
        return ComplexField(prec)


    # very useful to cache this.
    def _real_field(self):
        """
        Return the underlying real field with the same precision.

        EXAMPLES::

            sage: RF = ComplexField(10)._real_field(); RF
            Real Field with 10 bits of precision
            sage: ComplexField(10)._real_field() is RF
            True
        """
        try:
            return self.__real_field
        except AttributeError:
            self.__real_field = real_mpfr.RealField(self._prec)
            return self.__real_field

    def __cmp__(self, other):
        """
        Compare ``self`` to ``other``.

        If ``other`` is not a :class:`ComplexField_class', then this compares
        by their types. Otherwise it compares by their precision.

        EXAMPLES::

            sage: cmp(ComplexField(), ComplexField())
            0
            sage: cmp(ComplexField(10), ComplexField(15))
            -1
        """
        if not isinstance(other, ComplexField_class):
            return cmp(type(self), type(other))
        return cmp(self._prec, other._prec)

    def __call__(self, x=None, im=None):
        """
        Create a complex number.

        EXAMPLES::

            sage: CC(2) # indirect doctest
            2.00000000000000
            sage: CC(CC.0)
            1.00000000000000*I
            sage: CC('1+I')
            1.00000000000000 + 1.00000000000000*I
            sage: CC(2,3)
            2.00000000000000 + 3.00000000000000*I
            sage: CC(QQ[I].gen())
            1.00000000000000*I
            sage: CC.gen() + QQ[I].gen()
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '+': 'Complex Field with 53 bits of precision' and 'Number Field in I with defining polynomial x^2 + 1'

        In the absence of arguments we return zero::

            sage: a = CC(); a
            0.000000000000000
            sage: a.parent()
            Complex Field with 53 bits of precision
        """
        if x is None:
            return self.zero_element()
        # we leave this here to handle the imaginary parameter
        if im is not None:
            x = x, im
        return Parent.__call__(self, x)

    def _element_constructor_(self, x):
        """
        Construct a complex number.

        EXAMPLES::

            sage: CC((1,2)) # indirect doctest
            1.00000000000000 + 2.00000000000000*I

        Check that :trac:`14989` is fixed::

            sage: QQi = NumberField(x^2+1, 'i', embedding=CC(0,1))
            sage: i = QQi.order(QQi.gen()).gen(1)
            sage: CC(i)
            1.00000000000000*I

        """
        if not isinstance(x, (real_mpfr.RealNumber, tuple)):
            if isinstance(x, complex_double.ComplexDoubleElement):
                return complex_number.ComplexNumber(self, x.real(), x.imag())
            elif isinstance(x, str):
                # TODO: this is probably not the best and most
                # efficient way to do this.  -- Martin Albrecht
                return complex_number.ComplexNumber(self,
                            sage_eval(x.replace(' ',''), locals={"I":self.gen(),"i":self.gen()}))

            late_import()
            if isinstance(x, NumberFieldElement_quadratic):
                if isinstance(x.parent(), NumberField_quadratic) and list(x.parent().polynomial()) == [1, 0, 1]:
                    (re, im) = list(x)
                    return complex_number.ComplexNumber(self, re, im)

            try:
                return self(x.sage())
            except (AttributeError, TypeError):
                pass
            try:
                return x._complex_mpfr_field_( self )
            except AttributeError:
                pass
        return complex_number.ComplexNumber(self, x)

    def _coerce_map_from_(self, S):
        """
        The rings that canonically coerce to the MPFR complex field are:

        - This MPFR complex field, or any other of higher precision

        - Anything that canonically coerces to the mpfr real field
          with this prec

        EXAMPLES::

            sage: ComplexField(200)(1) + RealField(90)(1) # indirect doctest
            2.0000000000000000000000000
            sage: parent(ComplexField(200)(1) + RealField(90)(1)) # indirect doctest
            Complex Field with 90 bits of precision
            sage: CC.0 + RLF(1/3) # indirect doctest
            0.333333333333333 + 1.00000000000000*I
            sage: ComplexField(20).has_coerce_map_from(CDF)
            True
            sage: ComplexField(200).has_coerce_map_from(CDF)
            False
        """
        RR = self._real_field()
        if RR.has_coerce_map_from(S):
            return complex_number.RRtoCC(RR, self) * RR.coerce_map_from(S)
        if is_ComplexField(S) and S._prec >= self._prec:
            return self._generic_convert_map(S)
        late_import()
        if S in [AA, QQbar, CLF, RLF] or (S == CDF and self._prec <= 53):
            return self._generic_convert_map(S)
        return self._coerce_map_via([CLF], S)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: ComplexField() # indirect doctest
            Complex Field with 53 bits of precision
            sage: ComplexField(15) # indirect doctest
            Complex Field with 15 bits of precision
        """
        return "Complex Field with %s bits of precision"%self._prec

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: latex(ComplexField()) # indirect doctest
            \Bold{C}
            sage: latex(ComplexField(15)) # indirect doctest
            \Bold{C}
        """
        return "\\Bold{C}"

    def _sage_input_(self, sib, coerce):
        r"""
        Produce an expression which will reproduce this value when evaluated.

        EXAMPLES::

            sage: sage_input(CC, verify=True)
            # Verified
            CC
            sage: sage_input(ComplexField(25), verify=True)
            # Verified
            ComplexField(25)
            sage: k = (CC, ComplexField(75))
            sage: sage_input(k, verify=True)
            # Verified
            (CC, ComplexField(75))
            sage: sage_input((k, k), verify=True)
            # Verified
            CC75 = ComplexField(75)
            ((CC, CC75), (CC, CC75))
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: ComplexField(99)._sage_input_(SageInputBuilder(), False)
            {call: {atomic:ComplexField}({atomic:99})}
        """
        if self.prec() == 53:
            return sib.name('CC')

        v = sib.name('ComplexField')(sib.int(self.prec()))

        name = 'CC%d' % (self.prec())
        sib.cache(self, v, name)
        return v

    def characteristic(self):
        r"""
        Return the characteristic of `\CC`, which is 0.

        EXAMPLES::

            sage: ComplexField().characteristic()
            0
        """
        return integer.Integer(0)

    def gen(self, n=0):
        """
        Return the generator of the complex field.

        EXAMPLES::

            sage: ComplexField().gen(0)
            1.00000000000000*I
        """
        if n != 0:
            raise IndexError, "n must be 0"
        return complex_number.ComplexNumber(self, 0, 1)

    def is_field(self, proof = True):
        """
        Return ``True`` since the complex numbers are a field.

        EXAMPLES::

            sage: CC.is_field()
            True
        """
        return True

    def is_finite(self):
        """
        Return ``False`` since there are infinite number of complex numbers.

        EXAMPLES::

            sage: CC.is_finite()
            False
        """
        return False

    def construction(self):
        """
        Returns the functorial construction of ``self``, namely the algebraic
        closure of the real field with the same precision.

        EXAMPLES::

            sage: c, S = CC.construction(); S
            Real Field with 53 bits of precision
            sage: CC == c(S)
            True
        """
        from sage.categories.pushout import AlgebraicClosureFunctor
        return (AlgebraicClosureFunctor(), self._real_field())


    def random_element(self, component_max=1, *args, **kwds):
        r"""
        Returns a uniformly distributed random number inside a square
        centered on the origin (by default, the square `[-1,1] \times [-1,1]`).

        Passes additional arguments and keywords to underlying real field.

        EXAMPLES::

            sage: [CC.random_element() for _ in range(5)]
            [0.153636193785613 - 0.502987375247518*I,
             0.609589964322241 - 0.948854594338216*I,
             0.968393085385764 - 0.148483595843485*I,
             -0.908976099636549 + 0.126219184235123*I,
             0.461226845462901 - 0.0420335212948924*I]
            sage: CC6 = ComplexField(6)
            sage: [CC6.random_element(2^-20) for _ in range(5)]
            [-5.4e-7 - 3.3e-7*I, 2.1e-7 + 8.0e-7*I, -4.8e-7 - 8.6e-7*I, -6.0e-8 + 2.7e-7*I, 6.0e-8 + 1.8e-7*I]
            sage: [CC6.random_element(pi^20) for _ in range(5)]
            [6.7e8 - 5.4e8*I, -9.4e8 + 5.0e9*I, 1.2e9 - 2.7e8*I, -2.3e9 - 4.0e9*I, 7.7e9 + 1.2e9*I]

        Passes extra positional or keyword arguments through::

            sage: [CC.random_element(distribution='1/n') for _ in range(5)]
            [-0.900931453455899 - 0.932172283929307*I,
             0.327862582226912 + 0.828104487111727*I,
             0.246299162813240 + 0.588214960163442*I,
             0.892970599589521 - 0.266744694790704*I,
             0.878458776600692 - 0.905641181799996*I]
        """
        size = self._real_field()(component_max)
        re = self._real_field().random_element(-size, size, *args, **kwds)
        im = self._real_field().random_element(-size, size, *args, **kwds)
        return self(re, im)

    def pi(self):
        r"""
        Returns `\pi` as a complex number.

        EXAMPLES::

            sage: ComplexField().pi()
            3.14159265358979
            sage: ComplexField(100).pi()
            3.1415926535897932384626433833
        """
        return self(self._real_field().pi())

    def ngens(self):
        r"""
        The number of generators of this complex field as an `\RR`-algebra.

        There is one generator, namely ``sqrt(-1)``.

        EXAMPLES::

            sage: ComplexField().ngens()
            1
        """
        return 1

    def zeta(self, n=2):
        """
        Return a primitive `n`-th root of unity.

        INPUT:

        -  ``n`` - an integer (default: 2)

        OUTPUT: a complex `n`-th root of unity.

        EXAMPLES::

            sage: C = ComplexField()
            sage: C.zeta(2)
            -1.00000000000000
            sage: C.zeta(5)
            0.309016994374947 + 0.951056516295154*I
        """
        from integer import Integer
        n = Integer(n)
        if n == 1:
            x = self(1)
        elif n == 2:
            x = self(-1)
        elif n >= 3:
            # Use De Moivre
            # e^(2*pi*i/n) = cos(2pi/n) + i *sin(2pi/n)
            RR = self._real_field()
            pi = RR.pi()
            z = 2*pi/n
            x = complex_number.ComplexNumber(self, z.cos(), z.sin())
        x._set_multiplicative_order( n )
        return x

    def scientific_notation(self, status=None):
        """
        Set or return the scientific notation printing flag.

        If this flag is ``True`` then complex numbers with this space as parent
        print using scientific notation.

        EXAMPLES::

            sage: C = ComplexField()
            sage: C((0.025, 2))
            0.0250000000000000 + 2.00000000000000*I
            sage: C.scientific_notation(True)
            sage: C((0.025, 2))
            2.50000000000000e-2 + 2.00000000000000e0*I
            sage: C.scientific_notation(False)
            sage: C((0.025, 2))
            0.0250000000000000 + 2.00000000000000*I
        """
        return self._real_field().scientific_notation(status)

    def algebraic_closure(self):
        """
        Return the algebraic closure of ``self`` (which is itself).

        EXAMPLES::

            sage: CC
            Complex Field with 53 bits of precision
            sage: CC.algebraic_closure()
            Complex Field with 53 bits of precision
            sage: CC = ComplexField(1000)
            sage: CC.algebraic_closure() is CC
            True
        """
        return self

