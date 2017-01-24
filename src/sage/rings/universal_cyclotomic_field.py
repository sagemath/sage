r"""
Universal cyclotomic field

The universal cyclotomic field is the smallest subfield of the complex field
containing all roots of unity. It is also the maximal Galois Abelian extension
of the rational numbers.

The implementation simply wraps GAP Cyclotomic. As mentioned in their
documentation: arithmetical operations are quite expensive, so the use of
internally represented cyclotomics is not recommended for doing arithmetic over
number fields, such as calculations with matrices of cyclotomics.

.. NOTE::

    There used to be a native Sage version of the universal cyclotomic field
    written by Christian Stump (see :trac:`8327`). It was slower on most
    operations and it was decided to use a version based on libGAP instead (see
    :trac:`18152`). One main difference in the design choices is that GAP stores
    dense vectors whereas the native ones used Python dictionaries (storing only
    nonzero coefficients). Most operations are faster with libGAP except some
    operation on very sparse elements. All details can be found in
    :trac:`18152`.

REFERENCES:

.. [Bre97] \T. Breuer "Integral bases for subfields of cyclotomic fields" AAECC 8, 279--289 (1997).

EXAMPLES::

    sage: UCF = UniversalCyclotomicField(); UCF
    Universal Cyclotomic Field

To generate cyclotomic elements::

    sage: UCF.gen(5)
    E(5)
    sage: UCF.gen(5,2)
    E(5)^2

    sage: E = UCF.gen

Equality and inequality checks::

    sage: E(6,2) == E(6)^2 == E(3)
    True

    sage: E(6)^2 != E(3)
    False

Addition and multiplication::

    sage: E(2) * E(3)
    -E(3)
    sage: f = E(2) + E(3); f
    2*E(3) + E(3)^2

Inverses::

    sage: f^-1
    1/3*E(3) + 2/3*E(3)^2
    sage: f.inverse()
    1/3*E(3) + 2/3*E(3)^2
    sage: f * f.inverse()
    1

Conjugation and Galois conjugates::

    sage: f.conjugate()
    E(3) + 2*E(3)^2

    sage: f.galois_conjugates()
    [2*E(3) + E(3)^2, E(3) + 2*E(3)^2]
    sage: f.norm_of_galois_extension()
    3

One can create matrices and polynomials::

    sage: m = matrix(2,[E(3),1,1,E(4)]); m
    [E(3)    1]
    [   1 E(4)]
    sage: m.parent()
    Full MatrixSpace of 2 by 2 dense matrices over Universal Cyclotomic Field
    sage: m**2
    [                       -E(3) E(12)^4 - E(12)^7 - E(12)^11]
    [E(12)^4 - E(12)^7 - E(12)^11                            0]

    sage: m.charpoly()
    x^2 + (-E(12)^4 + E(12)^7 + E(12)^11)*x + E(12)^4 + E(12)^7 + E(12)^8

    sage: m.echelon_form()
    [1 0]
    [0 1]

    sage: m.pivots()
    (0, 1)

    sage: m.rank()
    2

    sage: R.<x> = PolynomialRing(UniversalCyclotomicField(), 'x')
    sage: E(3) * x - 1
    E(3)*x - 1

TESTS::

    sage: UCF.one()
    1
    sage: UCF.zero()
    0
    sage: UCF.one().is_one()
    True
    sage: UCF.one().is_zero()
    False
    sage: UCF.zero().is_zero()
    True

Check that :trac:`14240` is fixed::

    sage: K.<rho> = CyclotomicField(245)
    sage: h = K.random_element()
    sage: h_rho = rho.coordinates_in_terms_of_powers()(h)
    sage: h_ucf = sum( c * E(245, i) for (i, c) in enumerate(h_rho) )
    sage: h_ucf**2  # random
    -169539876343/589714020*E(245) + 27815735177/20058300*E(245)^2  + ... + + 7828432097501/842448600*E(245)^244

Check that :trac:`16130` is fixed::

    sage: mat = matrix(UCF, 2, [-4, 2*E(7)^6, -5*E(13)^3 + 5*E(13)^8 - 4*E(13)^9, 0])
    sage: mat._echelon_classical()
    [1 0]
    [0 1]

Check that :trac:`16631` is fixed::

    sage: UCF.one() / 2
    1/2
    sage: UCF.one() / 2r
    1/2

Check that :trac:`17117` is fixed::

    sage: e3 = UCF.gen(3)
    sage: N(e3)
    -0.500000000000000 + 0.866025403784439*I
    sage: real(e3)
    -1/2
    sage: imag(e3)
    -1/2*E(12)^7 + 1/2*E(12)^11

AUTHORS:

- Christian Stump (2013): initial Sage version (see :trac:`8327`)

- Vincent Delecroix (2015): complete rewriting using libgap (see :trac:`18152`)
"""
from sage.misc.cachefunc import cached_method
from sage.misc.superseded import deprecated_function_alias

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import FieldElement, parent
from sage.structure.coerce import py_scalar_to_element
from sage.categories.morphism import Morphism
from sage.rings.ring import Field

from sage.rings.integer import Integer
from sage.rings.rational import Rational

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.infinity import Infinity

# Deprecations from the old universal cyclotomic field
from sage.misc.lazy_import import lazy_import
lazy_import("sage.rings", "universal_cyclotomic_field", deprecation=18152)

libgap = GapElement_Integer = GapElement_Rational = GapElement_Cyclotomic = None
gap = gap3 = None


def late_import():
    r"""
    This function avoids importing libgap on startup. It is called once through
    the constructor of :class:`UniversalCyclotomicField`.

    EXAMPLES::

        sage: import sage.rings.universal_cyclotomic_field as ucf
        sage: _ = UniversalCyclotomicField()   # indirect doctest
        sage: ucf.libgap is None               # indirect doctest
        False
    """
    global gap, gap3, libgap
    global GapElement_Integer, GapElement_Rational, GapElement_Cyclotomic
    from sage.libs.gap.libgap import libgap
    from sage.libs.gap.element import (GapElement_Integer,
                                       GapElement_Rational,
                                       GapElement_Cyclotomic)
    from sage.interfaces import (gap, gap3)


class UCFtoQQbar(Morphism):
    r"""
    Conversion to ``QQbar``.

    EXAMPLES::

        sage: UCF = UniversalCyclotomicField()
        sage: QQbar(UCF.gen(3))
        -0.500000000000000? + 0.866025403784439?*I

        sage: CC(UCF.gen(7,2) + UCF.gen(7,6))
        0.400968867902419 + 0.193096429713794*I

        sage: complex(E(7)+E(7,2))
        (0.40096886790241915+1.7567593946498534j)
        sage: complex(UCF.one()/2)
        (0.5+0j)
    """
    def __init__(self, UCF):
        r"""
        INPUT:

        - ``UCF`` -- a universal cyclotomic field

        TESTS::

            sage: UCF = UniversalCyclotomicField()
            sage: UCF.coerce_embedding()
            Generic morphism:
              From: Universal Cyclotomic Field
              To:   Algebraic Field
        """
        from sage.rings.qqbar import QQbar
        Morphism.__init__(self, UCF, QQbar)

    def _call_(self, x):
        r"""
        TESTS::

            sage: UCF = UniversalCyclotomicField()
            sage: UCFtoQQbar = UCF.coerce_embedding()
            sage: UCFtoQQbar(UCF.gen(3))  # indirect doctest
            -0.500000000000000? + 0.866025403784439?*I

        Test that the bug reported in :trac:`19912` has been fixed::

            sage: UCFtoQQbar(UCF.gen(4)+1)
            I + 1
        """
        obj = x._obj
        QQbar = self.codomain()
        if obj.IsRat():
            return QQbar(obj.sage())
        k = obj.Conductor().sage()
        coeffs = obj.CoeffsCyc(k).sage()
        zeta = QQbar.zeta(k)
        return QQbar(sum(coeffs[a] * zeta**a for a in range(k)))


class UniversalCyclotomicFieldElement(FieldElement):
    def __init__(self, parent, obj):
        r"""
        INPUT:

        - ``parent`` - a universal cyclotomic field

        - ``obj`` - a libgap element (either an integer, a rational or a
          cyclotomic)

        TESTS::

            sage: UCF = UniversalCyclotomicField()
            sage: a = UCF.an_element()
            sage: TestSuite(a).run()
        """
        self._obj = obj
        FieldElement.__init__(self, parent)

    def __bool__(self):
        r"""
        TESTS::

            sage: UCF = UniversalCyclotomicField()
            sage: map(bool, [UCF.zero(), UCF.one(), UCF.gen(3), UCF.gen(5) + UCF.gen(5,3)])
            [False, True, True, True]
        """
        return bool(self._obj)

    __nonzero__ = __bool__

    def __reduce__(self):
        r"""
        TESTS::

            sage: UCF = UniversalCyclotomicField()
            sage: a = UCF.zero()
            sage: loads(dumps(a))
            0
            sage: parent(_)
            Universal Cyclotomic Field

            sage: b = UCF.gen(5,1) - 3*UCF.gen(5,4)
            sage: c = loads(dumps(b))
            sage: c
            E(5) - 3*E(5)^4
            sage: c == b
            True
            sage: parent(c)
            Universal Cyclotomic Field
        """
        return self.parent(), (str(self),)

    def __eq__(self, other):
        r"""
        Equality test.

        EXAMPLES::

            sage: UCF = UniversalCyclotomicField()
            sage: UCF.one() == 1
            True
            sage: 1 == UCF.one()
            True

            sage: UCF(2/3) == 2/3
            True
            sage: 2/3 == UCF(2/3)
            True

            sage: UCF.gen(3) == UCF.gen(5)
            False
            sage: UCF.gen(5) + UCF.gen(3) == UCF.gen(3) + UCF.gen(5)
            True

            sage: UCF.zero() == None
            False

            sage: QQbar.zeta(5) == UCF.gen(5)
            True
            sage: UCF.gen(5) == QQbar.zeta(5)
            True
            sage: QQbar.zeta(5) == UCF.gen(5,2)
            False
            sage: UCF.gen(5,2) == QQbar.zeta(5)
            False
        """
        if parent(self) is not parent(other):
            from sage.structure.element import get_coercion_model
            cm = get_coercion_model()
            try:
                self, other = cm.canonical_coercion(self, other)
            except TypeError:
                return False
            return self == other
        return self._obj == other._obj

    def __ne__(self, other):
        r"""
        Difference test.

        EXAMPLES::

            sage: UCF = UniversalCyclotomicField()
            sage: UCF.one() != 1
            False
            sage: 1 != UCF.one()
            False

            sage: UCF(2/3) != 3/2
            True
            sage: 3/2 != UCF(2/3)
            True

            sage: UCF.gen(3) != UCF.gen(5)
            True
            sage: UCF.gen(3) + UCF.gen(5) != UCF.gen(5) + UCF.gen(3)
            False

            sage: UCF.gen(7) != QQbar.zeta(7)
            False
            sage: UCF.gen(7,2) != QQbar.zeta(7)
            True
        """
        return not self == other

    def real(self):
        r"""
        Return the real part of this element.

        EXAMPLES::

            sage: E(3).real()
            -1/2
            sage: E(5).real()
            1/2*E(5) + 1/2*E(5)^4

            sage: a = E(5) - 2*E(3)
            sage: AA(a.real()) == QQbar(a).real()
            True
        """
        P = self.parent()
        return P.element_class(P, self._obj.RealPart())

    real_part = real

    def imag(self):
        r"""
        Return the imaginary part of this element.

        EXAMPLES::

            sage: E(3).imag()
            -1/2*E(12)^7 + 1/2*E(12)^11
            sage: E(5).imag()
            1/2*E(20) - 1/2*E(20)^9

            sage: a = E(5) - 2*E(3)
            sage: AA(a.imag()) == QQbar(a).imag()
            True
        """
        P = self.parent()
        return P.element_class(P, self._obj.ImaginaryPart())

    imag_part = imag

    def is_real(self):
        r"""
        Test whether this element is real.

        EXAMPLES::

            sage: E(3).is_real()
            False
            sage: (E(3) + E(3,2)).is_real()
            True

            sage: a = E(3) - 2*E(7)
            sage: a.real_part().is_real()
            True
            sage: a.imag_part().is_real()
            True
        """
        return self._obj.RealPart() == self._obj

    def conductor(self):
        r"""
        Return the conductor of ``self``.

        EXAMPLES::

            sage: E(3).conductor()
            3
            sage: (E(5) + E(3)).conductor()
            15
        """
        return self._obj.Conductor().sage()

    field_order = deprecated_function_alias(18152, conductor)

    def _symbolic_(self, R):
        r"""
        TESTS::

            sage: SR(E(7))
            e^(2/7*I*pi)
            sage: SR(E(5) + 2*E(5,2) + 3*E(5,3))
            -sqrt(5) + 1/4*I*sqrt(2*sqrt(5) + 10) - 1/4*I*sqrt(-2*sqrt(5) + 10) - 3/2

        Test that the bug reported in :trac:`19912` has been fixed::

            sage: SR(1+E(4))
            I + 1
        """
        from sage.symbolic.constants import pi
        from sage.symbolic.all import i as I
        k = self._obj.Conductor().sage()
        coeffs = self._obj.CoeffsCyc(k).sage()
        s = R.zero()
        for a in range(k):
            if coeffs[a]:
                s += coeffs[a] * (2 * a * I * pi / k).exp()
        return s

    def to_cyclotomic_field(self, R=None):
        r"""
        Return this element as an element of a cyclotomic field.

        EXAMPLES::

            sage: UCF = UniversalCyclotomicField()

            sage: UCF.gen(3).to_cyclotomic_field()
            zeta3
            sage: UCF.gen(3,2).to_cyclotomic_field()
            -zeta3 - 1

            sage: CF = CyclotomicField(5)
            sage: CF(E(5)) # indirect doctest
            zeta5

            sage: CF = CyclotomicField(7)
            sage: CF(E(5)) # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: Cannot coerce zeta5 into Cyclotomic Field of order 7 and
            degree 6

            sage: CF = CyclotomicField(10)
            sage: CF(E(5)) # indirect doctest
            zeta10^2

        Matrices are correctly dealt with::

            sage: M = Matrix(UCF,2,[E(3),E(4),E(5),E(6)]); M
            [   E(3)    E(4)]
            [   E(5) -E(3)^2]

            sage: Matrix(CyclotomicField(60),M) # indirect doctest
            [zeta60^10 - 1     zeta60^15]
            [    zeta60^12     zeta60^10]

        Using a non-standard embedding::

            sage: CF = CyclotomicField(5,embedding=CC(exp(4*pi*i/5)))
            sage: x = E(5)
            sage: CC(x)
            0.309016994374947 + 0.951056516295154*I
            sage: CC(CF(x))
            0.309016994374947 + 0.951056516295154*I

        Test that the bug reported in :trac:`19912` has been fixed::

            sage: a = 1+E(4); a
            1 + E(4)
            sage: a.to_cyclotomic_field()
            zeta4 + 1
        """
        from sage.rings.number_field.number_field import CyclotomicField
        k = self._obj.Conductor().sage()
        Rcan = CyclotomicField(k)
        if R is None:
            R = Rcan
        obj = self._obj
        if obj.IsRat():
            return R(obj.sage())
        zeta = Rcan.gen()
        coeffs = obj.CoeffsCyc(k).sage()
        return R(sum(coeffs[a] * zeta**a for a in range(k)))

    def __hash__(self):
        r"""
        EXAMPLES::

            sage: UCF = UniversalCyclotomicField()
            sage: hash(UCF.zero())  # indirect doctest
            0
            sage: hash(UCF.gen(3,2))
            313156239               # 32-bit
            1524600308199219855     # 64-bit

        TESTS:

        See :trac:`19514`::

            sage: hash(UCF.one())
            1
        """
        k = self._obj.Conductor().sage()
        coeffs = self._obj.CoeffsCyc(k).sage()
        if k == 1:
            return hash(coeffs[0])
        else:
            return hash((k,) + tuple(coeffs))

    def _algebraic_(self, R):
        r"""
        TESTS::

            sage: UCF = UniversalCyclotomicField()
            sage: AA(UCF.gen(5) + UCF.gen(5,4))
            0.618033988749895?
            sage: AA(UCF.gen(5))
            Traceback (most recent call last):
            ...
            ValueError: Cannot coerce algebraic number with non-zero imaginary
            part to algebraic real
        """
        from sage.rings.qqbar import QQbar
        return R(QQbar(self))

    def __float__(self):
        r"""
        TESTS::

            sage: float(E(7) + E(7,6))
            1.246979603717467
        """
        from sage.rings.real_mpfr import RR
        return float(RR(self))

    def __complex__(self):
        r"""
        TESTS::

            sage: complex(E(3))
            (-0.5+0.8660254037844386j)
        """
        f = self.parent().coerce_embedding()
        return complex(f(self))

    def _eval_complex_(self, R):
        r"""
        Return a complex value of this element in ``R``.

        TESTS::

            sage: CC(E(3))
            -0.500000000000000 + 0.866025403784439*I

        Check that :trac:`19825` is fixed::

            sage: CIF(E(3))
            -0.500000000000000? + 0.866025403784439?*I
            sage: CIF(E(5))
            0.309016994374948? + 0.9510565162951536?*I
            sage: CIF(E(12))
            0.86602540378444? + 0.50000000000000?*I

        If the input is real, the imaginary part is exactly 0::

            sage: CIF(E(17,2) + E(17,15))
            1.47801783444132?
            sage: _.imag().is_zero()
            True
        """
        if self._obj.IsRat():
            return R(self._obj.sage())

        k = self._obj.Conductor().sage()
        coeffs = self._obj.CoeffsCyc(k).sage()
        zeta = R.zeta(k)
        s = sum(coeffs[i] * zeta ** i for i in range(1, k))
        if self.is_real():
            return R(s.real())
        return s

    _complex_mpfi_ = _eval_complex_
    _complex_mpfr_field_ = _eval_complex_

    def _eval_real_(self, R):
        r"""
        Return a real value of this element in ``R``.

        TESTS::

            sage: RR(E(7) + E(7,6))
            1.24697960371747
            sage: 2*cos(2*pi/7).n()
            1.24697960371747
        """
        if not self.is_real():
            raise TypeError("self is not real")

        if self._obj.IsRat():
            return R(self._obj.sage())

        k = self._obj.Conductor().sage()
        coeffs = self._obj.CoeffsCyc(k).sage()
        t = (2 * R.pi()) / k
        return sum(coeffs[i] * (i * t).cos() for i in range(1, k))

    _mpfr_ = _eval_real_

    def _cmp_(self, other):
        r"""
        Comparison (using the complex embedding).

        TESTS::

            sage: UCF = UniversalCyclotomicField()
            sage: l = [UCF.gen(3), UCF.gen(3)+1, UCF.gen(5), UCF.gen(5,2),
            ....:      UCF.gen(4), 2*UCF.gen(4), UCF.gen(5)-22/3]
            sage: lQQbar = map(QQbar,l)
            sage: lQQbar.sort()
            sage: l.sort()
            sage: lQQbar == map(QQbar,l)
            True

            sage: for i in range(len(l)):
            ....:     assert l[i] >= l[i] and l[i] <= l[i]
            ....:     for j in range(i):
            ....:         assert l[i] > l[j] and l[j] < l[i]

            sage: fibonacci(200)*(E(5)+E(5,4)) <= fibonacci(199)
            True
            sage: fibonacci(201)*(E(5)+E(5,4)) <= fibonacci(200)
            False
        """
        if self._obj == other._obj:
            return 0

        s = self.real_part()
        o = other.real_part()
        if s == o:
            s = self.imag_part()
            o = other.imag_part()

        from sage.rings.real_mpfi import RealIntervalField
        prec = 53
        R = RealIntervalField(prec)
        sa = s._eval_real_(R)
        oa = o._eval_real_(R)
        while sa.overlaps(oa):
            prec <<= 2
            R = RealIntervalField(prec)
            sa = s._eval_real_(R)
            oa = o._eval_real_(R)
        return sa._cmp_(oa)

    def denominator(self):
        r"""
        Return the denominator of this element.

        EXAMPLES::

            sage: a = E(5) + 1/2*E(5,2) + 1/3*E(5,3)
            sage: a
            E(5) + 1/2*E(5)^2 + 1/3*E(5)^3
            sage: a.denominator()
            6
            sage: parent(_)
            Integer Ring
        """
        return self._obj.DenominatorCyc().sage()

    def multiplicative_order(self):
        r"""
        The multiplicative order.

        EXAMPLES::

            sage: E(5).multiplicative_order()
            5
            sage: (E(5) + E(12)).multiplicative_order()
            +Infinity
            sage: UniversalCyclotomicField().zero().multiplicative_order()
            Traceback (most recent call last):
            ...
            ValueError: libGAP: Error, argument must be nonzero
        """
        return self._obj.Order().sage()

    def additive_order(self):
        r"""
        The additive order.

        EXAMPLES::

            sage: UCF = UniversalCyclotomicField()
            sage: UCF.zero().additive_order()
            0
            sage: UCF.one().additive_order()
            +Infinity
            sage: UCF.gen(3).additive_order()
            +Infinity
        """
        return Infinity if self else ZZ.zero()

    def is_rational(self):
        r"""
        Test whether this element is a rational number.

        EXAMPLES::

            sage: E(3).is_rational()
            False
            sage: (E(3) + E(3,2)).is_rational()
            True

        TESTS::

            sage: type(E(3).is_rational())
            <type 'bool'>
        """
        return self._obj.IsRat().sage()

    def _rational_(self):
        r"""
        TESTS::

            sage: UCF = UniversalCyclotomicField()
            sage: QQ(UCF.zero())         # indirect doctest
            0
            sage: parent(_)
            Rational Field

            sage: QQ(UCF.one())          # indirect doctest
            1
            sage: parent(_)
            Rational Field

            sage: QQ(E(3)/2 + E(3,2)/2)  # indirect doctest
            -1/2
        """
        if not self._obj.IsRat():
            raise TypeError("Unable to coerce to a rational")
        return Rational(self._obj.sage())

    def _repr_(self):
        r"""
        TESTS::

            sage: U1 = UniversalCyclotomicField(names='E')
            sage: U2 = UniversalCyclotomicField(names='UCF')
            sage: U1.gen(5,2)
            E(5)^2
            sage: U2.gen(5,2)
            E(5)^2
        """
        s = str(self._obj)
        first_char = s[0]
        s = s[1:].replace('+', ' + ').replace('-', ' - ')
        return first_char + s

    def _add_(self, other):
        r"""
        TESTS::

            sage: E(3) + E(5)
            -E(15)^2 - 2*E(15)^8 - E(15)^11 - E(15)^13 - E(15)^14
            sage: 1/2 + E(3)
            1/2*E(3) - 1/2*E(3)^2
        """
        P = self.parent()
        return P.element_class(P, self._obj + other._obj)

    def _sub_(self, other):
        r"""
        TESTS::

            sage: E(3) - E(5)
            -E(15)^2 - E(15)^11 + E(15)^13 - E(15)^14
        """
        P = self.parent()
        return P.element_class(P, self._obj - other._obj)

    def __neg__(self):
        r"""
        Return the inverse of ``self``.

        TESTS::

            sage: -E(5)
            -E(5)
        """
        P = self.parent()
        return P.element_class(P, -self._obj)

    def _mul_(self, other):
        r"""
        TESTS::

            sage: E(3) * E(4)
            E(12)^7
            sage: 3 * E(4)
            3*E(4)
            sage: E(4) * 3
            3*E(4)
        """
        P = self.parent()
        return P.element_class(P, self._obj * other._obj)

    def _div_(self, other):
        r"""
        TESTS::

            sage: E(3)/2
            1/2*E(3)
            sage: 2/E(3)
            2*E(3)^2
        """
        P = self.parent()
        try:
            return P.element_class(P, self._obj / other._obj)
        except ValueError:
            raise ZeroDivisionError("division by zero")

    def __invert__(self):
        r"""
        TESTS::

            sage: UCF = UniversalCyclotomicField()
            sage: ~(UCF.one())
            1
            sage: ~UCF.gen(4)
            -E(4)
        """
        P = self.parent()
        return P.element_class(P, ~self._obj)

    inverse = __invert__

    def conjugate(self):
        r"""
        Return the complex conjugate.

        EXAMPLES::

            sage: (E(7) + 3*E(7,2) - 5 * E(7,3)).conjugate()
            -5*E(7)^4 + 3*E(7)^5 + E(7)^6
        """
        P = self.parent()
        return P.element_class(P, self._obj.ComplexConjugate())

    def galois_conjugates(self, n=None):
        r"""
        Return the Galois conjugates of ``self``.

        INPUT:

        - ``n`` -- an optional integer. If provided, return the orbit of the
          Galois group of the ``n``-th cyclotomic field over `\QQ`. Note that
          ``n`` must be such that this element belongs to the ``n``-th
          cyclotomic field (in other words, it must be a multiple of the
          conductor).

        EXAMPLES::

            sage: E(6).galois_conjugates()
            [-E(3)^2, -E(3)]

            sage: E(6).galois_conjugates()
            [-E(3)^2, -E(3)]

            sage: (E(9,2) - E(9,4)).galois_conjugates()
            [E(9)^2 - E(9)^4,
             E(9)^2 + E(9)^4 + E(9)^5,
             -E(9)^2 - E(9)^5 - E(9)^7,
             -E(9)^2 - E(9)^4 - E(9)^7,
             E(9)^4 + E(9)^5 + E(9)^7,
             -E(9)^5 + E(9)^7]

            sage: zeta = E(5)
            sage: zeta.galois_conjugates(5)
            [E(5), E(5)^2, E(5)^3, E(5)^4]
            sage: zeta.galois_conjugates(10)
            [E(5), E(5)^3, E(5)^2, E(5)^4]
            sage: zeta.galois_conjugates(15)
            [E(5), E(5)^2, E(5)^4, E(5)^2, E(5)^3, E(5), E(5)^3, E(5)^4]

            sage: zeta.galois_conjugates(17)
            Traceback (most recent call last):
            ...
            ValueError: n = 17 must be a multiple of the conductor (5)
        """
        P = self.parent()
        obj = self._obj
        k = obj.Conductor().sage()
        n = k if n is None else ZZ(n)
        if not k.divides(n):
            raise ValueError("n = {} must be a multiple of the conductor ({})".format(n, k))
        return [P.element_class(P, obj.GaloisCyc(i)) for i in range(n) if n.gcd(i) == 1]

    def norm_of_galois_extension(self):
        r"""
        Return the norm as a Galois extension of `\QQ`, which is
        given by the product of all galois_conjugates.

        EXAMPLES::

            sage: E(3).norm_of_galois_extension()
            1
            sage: E(6).norm_of_galois_extension()
            1
            sage: (E(2) + E(3)).norm_of_galois_extension()
            3
            sage: parent(_)
            Integer Ring
        """
        obj = self._obj
        k = obj.Conductor().sage()
        return libgap.Product(libgap([obj.GaloisCyc(i) for i in range(k) if k.gcd(i) == 1])).sage()

    def minpoly(self, var='x'):
        r"""
        The minimal polynomial of ``self`` element over `\QQ`.

        INPUT:

        - ``var`` -- (optional, default 'x') the name of the variable to use.

        EXAMPLES::

            sage: UCF.<E> = UniversalCyclotomicField()

            sage: UCF(4).minpoly()
            x - 4

            sage: UCF(4).minpoly(var='y')
            y - 4

            sage: E(3).minpoly()
            x^2 + x + 1

            sage: E(3).minpoly(var='y')
            y^2 + y + 1

        TESTS::

            sage: for elt in UCF.some_elements():
            ....:     assert elt.minpoly() == elt.to_cyclotomic_field().minpoly()
            ....:     assert elt.minpoly(var='y') == elt.to_cyclotomic_field().minpoly(var='y')

        .. TODO::

            Polynomials with libgap currently does not implement a ``.sage()`` method
            (see :trac:`18266`). It would be faster/safer to not use string to
            construct the polynomial.
        """
        gap_p = libgap.MinimalPolynomial(libgap.eval("Rationals"), self._obj)
        return QQ[var](QQ['x_1'](str(gap_p)))


class UniversalCyclotomicField(UniqueRepresentation, Field):
    r"""
    The universal cyclotomic field.

    The universal cyclotomic field is the infinite algebraic extension of `\QQ`
    generated by the roots of unity. It is also the maximal Abelian extension of
    `\QQ` in the sense that any Abelian Galois extension of `\QQ` is also a
    subfield of the universal cyclotomic field.
    """
    Element = UniversalCyclotomicFieldElement

    @staticmethod
    def __classcall__(cls, names=None):
        r"""
        Just ignoring the argument ``names``.

        TESTS::

            sage: UCF.<E> = UniversalCyclotomicField()
            sage: E(3,1)
            E(3)
            sage: E(3,2)
            E(3)^2
        """
        return super(UniversalCyclotomicField, cls).__classcall__(cls, None)

    def __init__(self, names=None):
        r"""
        TESTS::

            sage: UCF = UniversalCyclotomicField()
            sage: TestSuite(UCF).run()
        """
        from sage.categories.fields import Fields
        Field.__init__(self, base_ring=QQ, category=Fields().Infinite())
        self._populate_coercion_lists_(embedding=UCFtoQQbar(self))
        late_import()

    def _first_ngens(self, n):
        r"""
        Return the function :meth:`gen` if ``n=1``, and raises an error otherwise.

        This method is needed to make the following work::

            sage: UCF.<E> = UniversalCyclotomicField() # indirect doctest
        """
        if n == 1:
            return (self.gen,)
        else:
            raise ValueError("This ring has only a single generator method.")

    def an_element(self):
        r"""
        Return an element.

        EXAMPLES::

            sage: UniversalCyclotomicField().an_element()
            E(5) - 3*E(5)^2
        """
        return self.gen(5, 1) - self(3) * self.gen(5, 2)

    def some_elements(self):
        r"""
        Return a tuple of some elements in the universal cyclotomic field.

        EXAMPLES::

            sage: UniversalCyclotomicField().some_elements()
            (0, 1, -1, E(3), E(7) - 2/3*E(7)^2)
            sage: all(parent(x) is UniversalCyclotomicField() for x in _)
            True
        """
        return (self.zero(), self.one(), -self.one(),
                self.gen(3, 1),
                self.gen(7, 1) - self(2) / self(3) * self.gen(7, 2))

    def is_finite(self):
        r"""
        Return ``True``.

        EXAMPLES::

            sage: UniversalCyclotomicField().is_finite()
            True

        .. TODO::

            this method should be provided by the category.
        """
        return True

    def _repr_(self):
        r"""
        TESTS::

            sage: UniversalCyclotomicField()  # indirect doctest
            Universal Cyclotomic Field
        """
        return "Universal Cyclotomic Field"

    def is_exact(self):
        r"""
        Return ``True`` as this is an exact ring (i.e. not numerical).

        EXAMPLES::

            sage: UniversalCyclotomicField().is_exact()
            True
        """
        return True

    @cached_method
    def zero(self):
        r"""
        Return zero.

        EXAMPLES::

            sage: UCF = UniversalCyclotomicField()
            sage: UCF.zero()
            0
            sage: parent(_)
            Universal Cyclotomic Field
        """
        return self.element_class(self, libgap.zero())

    @cached_method
    def one(self):
        r"""
        Return one.

        EXAMPLES::

            sage: UCF = UniversalCyclotomicField()
            sage: UCF.one()
            1
            sage: parent(_)
            Universal Cyclotomic Field
        """
        return self.element_class(self, libgap.one())

    def characteristic(self):
        r"""
        Return the characteristic.

        EXAMPLES::

            sage: UniversalCyclotomicField().characteristic()
            0
            sage: parent(_)
            Integer Ring
        """
        return ZZ.zero()

    def gen(self, n, k=1):
        r"""
        Return the standard primitive ``n``-th root of unity.

        If ``k`` is not ``None``, return the ``k``-th power of it.

        EXAMPLES::

            sage: UCF = UniversalCyclotomicField()
            sage: UCF.gen(15)
            E(15)
            sage: UCF.gen(7,3)
            E(7)^3
            sage: UCF.gen(4,2)
            -1

        There is an alias ``zeta`` also available::

            sage: UCF.zeta(6)
            -E(3)^2
        """
        return self.element_class(self, libgap.E(n)**k)

    zeta = gen

    def _element_constructor_(self, elt):
        r"""
        TESTS::

            sage: UCF = UniversalCyclotomicField()
            sage: UCF(3)
            3
            sage: UCF(3/2)
            3/2

            sage: C = CyclotomicField(13)
            sage: UCF(C.gen())
            E(13)
            sage: UCF(C.gen() - 3*C.gen()**2 + 5*C.gen()**5)
            E(13) - 3*E(13)^2 + 5*E(13)^5

            sage: C = CyclotomicField(12)
            sage: zeta12 = C.gen()
            sage: a = UCF(zeta12 - 3* zeta12**2)
            sage: a
            -E(12)^7 + 3*E(12)^8
            sage: C(_) == a
            True

            sage: UCF('[[0, 1], [0, 2]]')
            Traceback (most recent call last):
            ...
            TypeError: [ [ 0, 1 ], [ 0, 2 ] ] of type <type
            'sage.libs.gap.element.GapElement_List'> not valid to initialize an
            element of the universal cyclotomic field

        .. TODO::

            Implement conversion from QQbar (and as a consequence from the
            symbolic ring)
        """
        elt = py_scalar_to_element(elt)

        if isinstance(elt, (Integer, Rational)):
            return self.element_class(self, libgap(elt))
        elif isinstance(elt, (GapElement_Integer, GapElement_Rational, GapElement_Cyclotomic)):
            return self.element_class(self, elt)
        elif not elt:
            return self.zero()

        obj = None
        if isinstance(elt, gap.GapElement):
            obj = libgap(elt)
        elif isinstance(elt, gap3.GAP3Element):
            obj = libgap.eval(str(elt))
        elif isinstance(elt, str):
            obj = libgap.eval(elt)
        if obj is not None:
            if not isinstance(obj, (GapElement_Integer, GapElement_Rational, GapElement_Cyclotomic)):
                raise TypeError("{} of type {} not valid to initialize an element of the universal cyclotomic field".format(obj, type(obj)))
            return self.element_class(self, obj)

        # late import to avoid slowing down the above conversions
        from sage.rings.number_field.number_field_element import NumberFieldElement
        from sage.rings.number_field.number_field import NumberField_cyclotomic, CyclotomicField
        P = parent(elt)
        if isinstance(elt, NumberFieldElement) and isinstance(P, NumberField_cyclotomic):
            n = P.gen().multiplicative_order()
            elt = CyclotomicField(n)(elt)
            return sum(c * self.gen(n, i)
                       for i, c in enumerate(elt._coefficients()))
        else:
            raise TypeError("{} of type {} not valid to initialize an element of the universal cyclotomic field".format(elt, type(elt)))

    def _coerce_map_from_(self, other):
        r"""
        TESTS::

            sage: UCF = UniversalCyclotomicField()
            sage: UCF.has_coerce_map_from(ZZ)
            True
            sage: UCF.has_coerce_map_from(QQ)
            True
            sage: ZZ.has_coerce_map_from(UCF)
            False
            sage: QQ.has_coerce_map_from(UCF)
            False
            sage: QQbar.has_coerce_map_from(UCF)
            True
            sage: CC.has_coerce_map_from(UCF)
            True
        """
        if other is ZZ or other is QQ:
            return True
        from sage.rings.number_field.number_field import NumberField_cyclotomic
        if isinstance(other, NumberField_cyclotomic):
            return True

    def degree(self):
        r"""
        Return the *degree* of ``self`` as a field extension over the Rationals.

        EXAMPLES::

            sage: UCF = UniversalCyclotomicField()
            sage: UCF.degree()
            +Infinity
        """
        return Infinity

    def _gap_init_(self):
        r"""
        Return gap string representation of ``self``.

        EXAMPLES::

            sage: UCF = UniversalCyclotomicField()
            sage: UCF._gap_init_()
            'Cyclotomics'
        """
        return 'Cyclotomics'

    def algebraic_closure(self):
        r"""
        The algebraic closure.

        EXAMPLES::

            sage: UniversalCyclotomicField().algebraic_closure()
            Algebraic Field
        """
        from sage.rings.qqbar import QQbar
        return QQbar


def E(n, k=1):
    r"""
    Return the ``n``-th root of unity as an element of the universal cyclotomic
    field.

    EXAMPLES::

        sage: E(3)
        E(3)
        sage: E(3) + E(5)
        -E(15)^2 - 2*E(15)^8 - E(15)^11 - E(15)^13 - E(15)^14
    """
    return UniversalCyclotomicField().gen(n, k)
