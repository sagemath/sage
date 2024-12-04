r"""
Computation of Frobenius matrix on Monsky-Washnitzer cohomology

The most interesting functions to be exported here are
:func:`matrix_of_frobenius` and :func:`adjusted_prec`.

Currently this code is limited to the case `p \geq 5` (no
`GF(p^n)` for `n > 1`), and only handles the
elliptic curve case (not more general hyperelliptic curves).

REFERENCES:

- [Ked2001]_

- [Edix]_

AUTHORS:

- David Harvey and Robert Bradshaw: initial code developed at the 2006
  MSRI graduate workshop, working with Jennifer Balakrishnan and Liang
  Xiao

- David Harvey (2006-08): cleaned up, rewrote some chunks, lots more
  documentation, added Newton iteration method, added more complete
  'trace trick', integrated better into Sage.

- David Harvey (2007-02): added algorithm with sqrt(p) complexity
  (removed in May 2007 due to better C++ implementation)

- Robert Bradshaw (2007-03): keep track of exact form in reduction
  algorithms

- Robert Bradshaw (2007-04): generalization to hyperelliptic curves

- Julian Rueth (2014-05-09): improved caching
"""

# ****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#                     2006 Robert Bradshaw <robertwb@math.washington.edu>
#                     2006 David Harvey <dmharvey@math.harvard.edu>
#                     2014 Julian Rueth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.arith.misc import binomial
from sage.arith.misc import integer_ceil as ceil
from sage.categories.algebras import Algebras
from sage.categories.integral_domains import IntegralDomains
from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_import import lazy_import
from sage.misc.misc import newton_method_sizes
from sage.misc.profiler import Profiler
from sage.misc.repr import repr_lincomb
from sage.modules.free_module import FreeModule
from sage.modules.free_module_element import FreeModuleElement, vector
from sage.modules.module import Module
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing as Integers
from sage.rings.infinity import Infinity
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.laurent_series_ring import LaurentSeriesRing
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.rational import Rational
from sage.rings.rational_field import QQ
from sage.rings.rational_field import RationalField as Rationals
from sage.schemes.elliptic_curves.constructor import EllipticCurve
from sage.schemes.elliptic_curves.ell_generic import EllipticCurve_generic
from sage.schemes.hyperelliptic_curves_smooth_model.hyperelliptic_constructor import (
    HyperellipticCurveSmoothModel,
)
from sage.schemes.hyperelliptic_curves_smooth_model.hyperelliptic_generic import (
    HyperellipticCurveSmoothModel_generic,
)
from sage.structure.element import ModuleElement
from sage.structure.parent import Parent
from sage.structure.richcmp import richcmp
from sage.structure.unique_representation import UniqueRepresentation

# TODO: why are these lazy imports used?
lazy_import("sage.functions.log", "log")
lazy_import("sage.rings.lazy_series_ring", "LazyLaurentSeriesRing")
lazy_import("sage.rings.padics.factory", "Qp", as_="pAdicField")


class SpecialCubicQuotientRingElement(ModuleElement):
    """
    An element of a :class:`SpecialCubicQuotientRing`.
    """

    def __init__(self, parent, p0, p1, p2, check=True):
        """
        Construct the element `p_0 + p_1*x + p_2*x^2`, where
        the `p_i` are polynomials in `T`.

        INPUT:

        - ``parent`` -- a :class:`SpecialCubicQuotientRing`

        - ``p0``, ``p1``, ``p2`` -- coefficients; must be coercible
          into parent.poly_ring()

        - ``check`` -- boolean (default: ``True``); whether to carry
          out coercion

        EXAMPLES::

            sage: B.<t> = PolynomialRing(Integers(125))
            sage: from sage.schemes.hyperelliptic_curves_smooth_model.monsky_washnitzer import SpecialCubicQuotientRing, SpecialCubicQuotientRingElement
            sage: R = SpecialCubicQuotientRing(t^3 - t + B(1/4))
            sage: SpecialCubicQuotientRingElement(R, 2, 3, 4)
            (2) + (3)*x + (4)*x^2

        TESTS::

            sage: B.<t> = PolynomialRing(Integers(125))
            sage: from sage.schemes.hyperelliptic_curves_smooth_model.monsky_washnitzer import SpecialCubicQuotientRing
            sage: R = SpecialCubicQuotientRing(t^3 - t + B(1/4))
            sage: TestSuite(R).run()
            sage: p = R.create_element(t, t^2 - 2, 3)
            sage: -p
            (124*T) + (124*T^2 + 2)*x + (122)*x^2
        """
        if not isinstance(parent, SpecialCubicQuotientRing):
            raise TypeError(f"parent (={parent}) must be a SpecialCubicQuotientRing")

        ModuleElement.__init__(self, parent)

        if check:
            poly_ring = parent._poly_ring
            p0 = poly_ring(p0)
            p1 = poly_ring(p1)
            p2 = poly_ring(p2)

        self._triple = (p0, p1, p2)

    def coeffs(self):
        """
        Return list of three lists of coefficients, corresponding to the
        `x^0`, `x^1`, `x^2` coefficients.

        The lists are zero padded to the same length. The list entries
        belong to the base ring.

        EXAMPLES::

            sage: B.<t> = PolynomialRing(Integers(125))
            sage: from sage.schemes.hyperelliptic_curves_smooth_model.monsky_washnitzer import SpecialCubicQuotientRing
            sage: R = SpecialCubicQuotientRing(t^3 - t + B(1/4))
            sage: p = R.create_element(t, t^2 - 2, 3)
            sage: p.coeffs()
            [[0, 1, 0], [123, 0, 1], [3, 0, 0]]
        """
        coeffs = [column.coefficients(sparse=False) for column in self._triple]
        degree = max([len(x) for x in coeffs])
        base_ring = self.parent().base_ring()
        for column in coeffs:
            column.extend([base_ring(0)] * (degree - len(column)))
        return coeffs

    def __bool__(self) -> bool:
        """
        EXAMPLES::

            sage: B.<t> = PolynomialRing(Integers(125))
            sage: from sage.schemes.hyperelliptic_curves_smooth_model.monsky_washnitzer import SpecialCubicQuotientRing
            sage: R = SpecialCubicQuotientRing(t^3 - t + B(1/4))
            sage: x, T = R.gens()
            sage: not x
            False
            sage: not T
            False
            sage: not R.create_element(0, 0, 0)
            True
        """
        return bool(self._triple[0]) or bool(self._triple[1]) or bool(self._triple[2])

    def _richcmp_(self, other, op) -> bool:
        """
        EXAMPLES::

            sage: B.<t> = PolynomialRing(Integers(125))
            sage: x, t = monsky_washnitzer.SpecialCubicQuotientRing(t^3 - t + B(1/4)).gens()
            sage: x == t
            False
            sage: x == x
            True
            sage: x == x + x - x
            True
        """
        return richcmp(self._triple, other._triple, op)

    def _repr_(self) -> str:
        """
        EXAMPLES::

            sage: B.<t> = PolynomialRing(Integers(125))
            sage: from sage.schemes.hyperelliptic_curves_smooth_model.monsky_washnitzer import SpecialCubicQuotientRing
            sage: R = SpecialCubicQuotientRing(t^3 - t + B(1/4))
            sage: x, T = R.gens()
            sage: x + T*x - 2*T^2
            (123*T^2) + (T + 1)*x + (0)*x^2
        """
        return "({}) + ({})*x + ({})*x^2".format(*self._triple)

    def _latex_(self) -> str:
        """
        EXAMPLES::

            sage: B.<t> = PolynomialRing(Integers(125))
            sage: from sage.schemes.hyperelliptic_curves_smooth_model.monsky_washnitzer import SpecialCubicQuotientRing
            sage: R = SpecialCubicQuotientRing(t^3 - t + B(1/4))
            sage: x, T = R.gens()
            sage: f = x + T*x - 2*T^2
            sage: latex(f)
            (123 T^{2}) + (T + 1)x + (0)x^2
        """
        return "(%s) + (%s)x + (%s)x^2" % tuple(
            column._latex_() for column in self._triple
        )

    def _add_(self, other):
        """
        EXAMPLES::

            sage: B.<t> = PolynomialRing(Integers(125))
            sage: from sage.schemes.hyperelliptic_curves_smooth_model.monsky_washnitzer import SpecialCubicQuotientRing
            sage: R = SpecialCubicQuotientRing(t^3 - t + B(1/4))
            sage: f = R.create_element(2, t, t^2 - 3)
            sage: g = R.create_element(3 + t, -t, t)
            sage: f + g
            (T + 5) + (0)*x + (T^2 + T + 122)*x^2
        """
        P = self.parent()
        return P.element_class(
            P,
            self._triple[0] + other._triple[0],
            self._triple[1] + other._triple[1],
            self._triple[2] + other._triple[2],
            check=False,
        )

    def _sub_(self, other):
        """
        EXAMPLES::

            sage: B.<t> = PolynomialRing(Integers(125))
            sage: from sage.schemes.hyperelliptic_curves_smooth_model.monsky_washnitzer import SpecialCubicQuotientRing
            sage: R = SpecialCubicQuotientRing(t^3 - t + B(1/4))
            sage: f = R.create_element(2, t, t^2 - 3)
            sage: g = R.create_element(3 + t, -t, t)
            sage: f - g
            (124*T + 124) + (2*T)*x + (T^2 + 124*T + 122)*x^2
        """
        P = self.parent()
        return P.element_class(
            P,
            self._triple[0] - other._triple[0],
            self._triple[1] - other._triple[1],
            self._triple[2] - other._triple[2],
            check=False,
        )

    def shift(self, n):
        """
        Return this element multiplied by `T^n`.

        EXAMPLES::

            sage: B.<t> = PolynomialRing(Integers(125))
            sage: from sage.schemes.hyperelliptic_curves_smooth_model.monsky_washnitzer import SpecialCubicQuotientRing
            sage: R = SpecialCubicQuotientRing(t^3 - t + B(1/4))
            sage: f = R.create_element(2, t, t^2 - 3)
            sage: f
            (2) + (T)*x + (T^2 + 122)*x^2
            sage: f.shift(1)
            (2*T) + (T^2)*x + (T^3 + 122*T)*x^2
            sage: f.shift(2)
            (2*T^2) + (T^3)*x + (T^4 + 122*T^2)*x^2
        """
        P = self.parent()
        return P.element_class(
            P,
            self._triple[0].shift(n),
            self._triple[1].shift(n),
            self._triple[2].shift(n),
            check=False,
        )

    def scalar_multiply(self, scalar):
        """
        Multiply this element by a scalar, i.e. just multiply each
        coefficient of `x^j` by the scalar.

        INPUT:

        - ``scalar`` -- either an element of ``base_ring``, or an
          element of ``poly_ring``

        EXAMPLES::

            sage: B.<t> = PolynomialRing(Integers(125))
            sage: from sage.schemes.hyperelliptic_curves_smooth_model.monsky_washnitzer import SpecialCubicQuotientRing
            sage: R = SpecialCubicQuotientRing(t^3 - t + B(1/4))
            sage: x, T = R.gens()
            sage: f = R.create_element(2, t, t^2 - 3)
            sage: f
            (2) + (T)*x + (T^2 + 122)*x^2
            sage: f.scalar_multiply(2)
            (4) + (2*T)*x + (2*T^2 + 119)*x^2
            sage: f.scalar_multiply(t)
            (2*T) + (T^2)*x + (T^3 + 122*T)*x^2
        """
        P = self.parent()
        scalar = P._poly_ring(scalar)
        return P.element_class(
            P,
            scalar * self._triple[0],
            scalar * self._triple[1],
            scalar * self._triple[2],
            check=False,
        )

    def square(self):
        """
        Return the square of the element.

        EXAMPLES::

            sage: B.<t> = PolynomialRing(Integers(125))
            sage: from sage.schemes.hyperelliptic_curves_smooth_model.monsky_washnitzer import SpecialCubicQuotientRing
            sage: R = SpecialCubicQuotientRing(t^3 - t + B(1/4))
            sage: x, T = R.gens()

        ::

            sage: f = R.create_element(1 + 2*t + 3*t^2, 4 + 7*t + 9*t^2, 3 + 5*t + 11*t^2)
            sage: f.square()
            (73*T^5 + 16*T^4 + 38*T^3 + 39*T^2 + 70*T + 120)
            + (121*T^5 + 113*T^4 + 73*T^3 + 8*T^2 + 51*T + 61)*x
            + (18*T^4 + 60*T^3 + 22*T^2 + 108*T + 31)*x^2
        """
        return self * self

    def _mul_(self, other):
        """
        EXAMPLES::

            sage: B.<t> = PolynomialRing(Integers(125))
            sage: from sage.schemes.hyperelliptic_curves_smooth_model.monsky_washnitzer import SpecialCubicQuotientRing
            sage: R = SpecialCubicQuotientRing(t^3 - t + B(1/4))
            sage: x, T = R.gens()

        ::

            sage: f = R.create_element(1 + 2*t + 3*t^2, 4 + 7*t + 9*t^2, 3 + 5*t + 11*t^2)
            sage: g = R.create_element(4 + 3*t + 7*t^2, 2 + 3*t + t^2, 8 + 4*t + 6*t^2)
            sage: f * g
            (65*T^5 + 27*T^4 + 33*T^3 + 75*T^2 + 120*T + 57)
            + (66*T^5 + T^4 + 123*T^3 + 95*T^2 + 24*T + 50)*x
            + (45*T^4 + 75*T^3 + 37*T^2 + 2*T + 52)*x^2
        """
        # Here we do Toom-Cook three-way multiplication, which reduces
        # the naive 9 polynomial multiplications to only 5 polynomial
        # multiplications.

        a0, a1, a2 = self._triple
        b0, b1, b2 = other._triple
        M = self.parent()._speedup_matrix

        if self is other:
            # faster method if we are squaring
            p0 = a0 * a0
            temp = a0 + 2 * a1 + 4 * a2
            p1 = temp * temp
            temp = a0 + a1 + a2
            p2 = temp * temp
            temp = 4 * a0 + 2 * a1 + a2
            p3 = temp * temp
            p4 = a2 * a2

        else:
            p0 = a0 * b0
            p1 = (a0 + 2 * a1 + 4 * a2) * (b0 + 2 * b1 + 4 * b2)
            p2 = (a0 + a1 + a2) * (b0 + b1 + b2)
            p3 = (4 * a0 + 2 * a1 + a2) * (4 * b0 + 2 * b1 + b2)
            p4 = a2 * b2

        q1 = p1 - p0 - 16 * p4
        q2 = p2 - p0 - p4
        q3 = p3 - 16 * p0 - p4

        c0 = p0
        c1 = M[0] * q1 + M[1] * q2 + M[2] * q3
        c2 = M[3] * q1 + M[4] * q2 + M[5] * q3
        c3 = M[6] * q1 + M[7] * q2 + M[8] * q3
        c4 = p4

        # Now the product is c0 + c1 x + c2 x^2 + c3 x^3 + c4 x^4.
        # We need to reduce mod y = x^3 + ax + b and return result.

        parent = self.parent()
        T = parent._poly_generator
        b = parent._b
        a = parent._a

        # todo: These lines are necessary to get binop stuff working
        # for certain base rings, e.g. when we compute b*c3 in the
        # final line. They shouldn't be necessary. Need to fix this
        # somewhere else in Sage.
        a = parent._poly_ring(a)
        b = parent._poly_ring(b)

        return parent.element_class(
            parent,
            -b * c3 + c0 + c3 * T,
            -b * c4 - a * c3 + c1 + c4 * T,
            -a * c4 + c2,
            check=False,
        )


class SpecialCubicQuotientRing(UniqueRepresentation, Parent):
    r"""
    Specialised class for representing the quotient ring
    `R[x,T]/(T - x^3 - ax - b)`, where `R` is an
    arbitrary commutative base ring (in which 2 and 3 are invertible),
    `a` and `b` are elements of that ring.

    Polynomials are represented internally in the form
    `p_0 + p_1 x + p_2 x^2` where the `p_i` are
    polynomials in `T`. Multiplication of polynomials always
    reduces high powers of `x` (i.e. beyond `x^2`) to
    powers of `T`.

    Hopefully this ring is faster than a general quotient ring because
    it uses the special structure of this ring to speed multiplication
    (which is the dominant operation in the frobenius matrix
    calculation). I haven't actually tested this theory though...

    .. TODO::

        Eventually we will want to run this in characteristic 3, so we
        need to: (a) Allow `Q(x)` to contain an `x^2` term, and (b) Remove
        the requirement that 3 be invertible. Currently this is used in
        the Toom-Cook algorithm to speed multiplication.

    EXAMPLES::

        sage: B.<t> = PolynomialRing(Integers(125))
        sage: from sage.schemes.hyperelliptic_curves_smooth_model.monsky_washnitzer import SpecialCubicQuotientRing
        sage: R = SpecialCubicQuotientRing(t^3 - t + B(1/4))
        sage: R
        SpecialCubicQuotientRing over Ring of integers modulo 125
        with polynomial T = x^3 + 124*x + 94
        sage: TestSuite(R).run()

    Get generators::

        sage: x, T = R.gens()
        sage: x
        (0) + (1)*x + (0)*x^2
        sage: T
        (T) + (0)*x + (0)*x^2

    Coercions::

        sage: R(7)
        (7) + (0)*x + (0)*x^2

    Create elements directly from polynomials::

        sage: A = R.poly_ring()
        sage: A
        Univariate Polynomial Ring in T over Ring of integers modulo 125
        sage: z = A.gen()
        sage: R.create_element(z^2, z+1, 3)
        (T^2) + (T + 1)*x + (3)*x^2

    Some arithmetic::

        sage: x^3
        (T + 31) + (1)*x + (0)*x^2
        sage: 3 * x**15 * T**2 + x - T
        (3*T^7 + 90*T^6 + 110*T^5 + 20*T^4 + 58*T^3 + 26*T^2 + 124*T) +
        (15*T^6 + 110*T^5 + 35*T^4 + 63*T^2 + 1)*x +
        (30*T^5 + 40*T^4 + 8*T^3 + 38*T^2)*x^2

    Retrieve coefficients (output is zero-padded)::

        sage: x^10
        (3*T^2 + 61*T + 8) + (T^3 + 93*T^2 + 12*T + 40)*x + (3*T^2 + 61*T + 9)*x^2
        sage: (x^10).coeffs()
        [[8, 61, 3, 0], [40, 12, 93, 1], [9, 61, 3, 0]]

    .. TODO::

        write an example checking multiplication of these polynomials
        against Sage's ordinary quotient ring arithmetic. I cannot seem
        to get the quotient ring stuff happening right now...
    """

    def __init__(self, Q, laurent_series=False):
        """
        Constructor.

        INPUT:

        - ``Q`` -- a polynomial of the form
          `Q(x) = x^3 + ax + b`, where `a`, `b` belong to a ring in which
          2, 3 are invertible.

        - ``laurent_series`` -- boolean (default: ``False``); whether or not to allow
          negative powers of `T`

        EXAMPLES::

            sage: B.<t> = PolynomialRing(Integers(125))
            sage: from sage.schemes.hyperelliptic_curves_smooth_model.monsky_washnitzer import SpecialCubicQuotientRing
            sage: R = SpecialCubicQuotientRing(t^3 - t + B(1/4))
            sage: R
            SpecialCubicQuotientRing over Ring of integers modulo 125
            with polynomial T = x^3 + 124*x + 94

        ::

            sage: from sage.schemes.hyperelliptic_curves_smooth_model.monsky_washnitzer import SpecialCubicQuotientRing
            sage: R = SpecialCubicQuotientRing(t^3 + 2*t^2 - t + B(1/4))
            Traceback (most recent call last):
            ...
            ValueError: Q (=t^3 + 2*t^2 + 124*t + 94) must be of the form x^3 + ax + b

        ::

            sage: B.<t> = PolynomialRing(Integers(10))
            sage: from sage.schemes.hyperelliptic_curves_smooth_model.monsky_washnitzer import SpecialCubicQuotientRing
            sage: R = SpecialCubicQuotientRing(t^3 - t + 1)
            Traceback (most recent call last):
            ...
            ArithmeticError: 2 and 3 must be invertible in the coefficient ring
            (=Ring of integers modulo 10) of Q
        """
        if not isinstance(Q, Polynomial):
            raise TypeError("Q (=%s) must be a polynomial" % Q)

        if Q.degree() != 3 or not Q[2].is_zero():
            raise ValueError("Q (=%s) must be of the form x^3 + ax + b" % Q)

        base_ring = Q.parent().base_ring()

        if not base_ring(6).is_unit():
            raise ArithmeticError(
                "2 and 3 must be invertible in the "
                "coefficient ring (=%s) of Q" % base_ring
            )

        self._a = Q[1]
        self._b = Q[0]
        if laurent_series:
            self._poly_ring = LaurentSeriesRing(base_ring, "T")  # R[T]
        else:
            self._poly_ring = PolynomialRing(base_ring, "T")  # R[T]
        self._poly_generator = self._poly_ring.gen(0)  # the generator T
        Parent.__init__(
            self, base=base_ring, category=Algebras(base_ring).Commutative()
        )

        # Precompute a matrix that is used in the Toom-Cook multiplication.
        # This is where we need 2 and 3 invertible.

        # a good description of Toom-Cook is online at:
        # https://gmplib.org/manual/Multiplication-Algorithms
        m = matrix(QQ, [[1, -12, 2], [-3, 30, -3], [2, -12, 1]]) / 6
        self._speedup_matrix = m.change_ring(base_ring).list()

    def _repr_(self) -> str:
        """
        String representation.

        EXAMPLES::

            sage: B.<t> = PolynomialRing(Integers(125))
            sage: from sage.schemes.hyperelliptic_curves_smooth_model.monsky_washnitzer import SpecialCubicQuotientRing
            sage: R = SpecialCubicQuotientRing(t^3 - t + B(1/4))
            sage: R
            SpecialCubicQuotientRing over Ring of integers modulo 125
            with polynomial T = x^3 + 124*x + 94
        """
        return "SpecialCubicQuotientRing over %s with polynomial T = %s" % (
            self.base_ring(),
            PolynomialRing(self.base_ring(), "x")([self._b, self._a, 0, 1]),
        )

    def poly_ring(self):
        """
        Return the underlying polynomial ring in `T`.

        EXAMPLES::

            sage: B.<t> = PolynomialRing(Integers(125))
            sage: from sage.schemes.hyperelliptic_curves_smooth_model.monsky_washnitzer import SpecialCubicQuotientRing
            sage: R = SpecialCubicQuotientRing(t^3 - t + B(1/4))
            sage: R.poly_ring()
            Univariate Polynomial Ring in T over Ring of integers modulo 125
        """
        return self._poly_ring

    def gens(self):
        """
        Return a list [x, T] where x and T are the generators of the ring
        (as element *of this ring*).

        .. NOTE::

             I have no idea if this is compatible with the usual Sage
             'gens' interface.

        EXAMPLES::

            sage: B.<t> = PolynomialRing(Integers(125))
            sage: from sage.schemes.hyperelliptic_curves_smooth_model.monsky_washnitzer import SpecialCubicQuotientRing
            sage: R = SpecialCubicQuotientRing(t^3 - t + B(1/4))
            sage: x, T = R.gens()
            sage: x
            (0) + (1)*x + (0)*x^2
            sage: T
            (T) + (0)*x + (0)*x^2
        """
        zero = self._poly_ring.zero()
        one = self._poly_ring.one()
        return (
            self.element_class(self, zero, one, zero, check=False),
            self.element_class(self, self._poly_generator, zero, zero, check=False),
        )

    def _element_constructor_(self, *args, check=True):
        """
        Create the element `p_0 + p_1*x + p_2*x^2`, where the `p_i`
        are polynomials in `T`.

        INPUT:

        - ``p0``, ``p1``, ``p2`` -- coefficients; must be coercible
          into poly_ring()

        - ``check`` -- boolean (default: ``True``); whether to carry
          out coercion

        EXAMPLES::

            sage: B.<t> = PolynomialRing(Integers(125))
            sage: from sage.schemes.hyperelliptic_curves_smooth_model.monsky_washnitzer import SpecialCubicQuotientRing
            sage: R = SpecialCubicQuotientRing(t^3 - t + B(1/4))
            sage: A, z = R.poly_ring().objgen()
            sage: R.create_element(z^2, z+1, 3)  # indirect doctest
            (T^2) + (T + 1)*x + (3)*x^2
        """
        if len(args) == 1 and args[0] in self.base_ring():
            p0 = self._poly_ring.coerce(args[0])
            p1 = p2 = self._poly_ring.zero()
        else:
            p0, p1, p2 = args
        return self.element_class(self, p0, p1, p2, check=check)

    create_element = _element_constructor_

    @cached_method
    def one(self):
        """
        Return the unit of ``self``.

        EXAMPLES::

            sage: B.<t> = PolynomialRing(Integers(125))
            sage: from sage.schemes.hyperelliptic_curves_smooth_model.monsky_washnitzer import SpecialCubicQuotientRing
            sage: R = SpecialCubicQuotientRing(t^3 - t + B(1/4))
            sage: R.one()
            (1) + (0)*x + (0)*x^2
        """
        R = self._poly_ring
        return self.element_class(self, R.one(), R.zero(), R.zero(), check=False)

    def _coerce_map_from_(self, R):
        """
        Coercion system.

        EXAMPLES::

            sage: Z125 = Integers(125)
            sage: B.<t> = PolynomialRing(Z125)
            sage: from sage.schemes.hyperelliptic_curves_smooth_model.monsky_washnitzer import SpecialCubicQuotientRing
            sage: R = SpecialCubicQuotientRing(t^3 - t + B(1/4))
            sage: R.has_coerce_map_from(Z125)
            True
        """
        return self._poly_ring.has_coerce_map_from(R)

    Element = SpecialCubicQuotientRingElement


def transpose_list(input) -> list[list]:
    """
    INPUT:

    - ``input`` -- list of lists, each list of the same length

    OUTPUT: list of lists such that ``output[i][j] = input[j][i]``

    EXAMPLES::

        sage: from sage.schemes.hyperelliptic_curves_smooth_model.monsky_washnitzer import transpose_list
        sage: L = [[1, 2], [3, 4], [5, 6]]
        sage: transpose_list(L)
        [[1, 3, 5], [2, 4, 6]]
    """
    h = len(input)
    w = len(input[0])

    output = []
    for i in range(w):
        row = []
        for j in range(h):
            row.append(input[j][i])
        output.append(row)
    return output


def helper_matrix(Q):
    r"""
    Compute the (constant) matrix used to calculate the linear
    combinations of the `d(x^i y^j)` needed to eliminate the
    negative powers of `y` in the cohomology (i.e., in
    :func:`reduce_negative`).

    INPUT:

    - ``Q`` -- cubic polynomial

    EXAMPLES::

        sage: t = polygen(QQ,'t')
        sage: from sage.schemes.hyperelliptic_curves_smooth_model.monsky_washnitzer import helper_matrix
        sage: helper_matrix(t**3-4*t-691)
        [     64/12891731  -16584/12891731 4297329/12891731]
        [   6219/12891731     -32/12891731    8292/12891731]
        [    -24/12891731    6219/12891731     -32/12891731]
    """
    a = Q[1]
    b = Q[0]

    # Discriminant (should be invertible for a curve of good reduction)
    D = 4 * a**3 + 27 * b**2
    Dinv = D ** (-1)  # NB do not use 1/D

    # This is the inverse of the matrix
    #   [  a,  -3b,    0 ]
    #   [  0,  -2a,  -3b ]
    #   [  3,    0,  -2a ]

    return Dinv * matrix(
        [
            [4 * a**2, -6 * b * a, 9 * b**2],
            [-9 * b, -2 * a**2, 3 * b * a],
            [6 * a, -9 * b, -2 * a**2],
        ]
    )


def lift(x):
    r"""
    Try to call ``x.lift()``, presumably from the `p`-adics to `\ZZ`.

    If this fails, it assumes the input is a power series, and tries to
    lift it to a power series over `\QQ`.

    This function is just a very kludgy solution to the problem of
    trying to make the reduction code (below) work over both `\ZZ_p` and
    `\ZZ_p[[t]]`.

    EXAMPLES::

        sage: # needs sage.rings.padics
        sage: from sage.schemes.hyperelliptic_curves_smooth_model.monsky_washnitzer import lift
        sage: l = lift(Qp(13)(131)); l
        131
        sage: l.parent()
        Integer Ring
        sage: x = PowerSeriesRing(Qp(17),'x').gen()
        sage: l = lift(4 + 5*x + 17*x**6); l
        4 + 5*t + 17*t^6
        sage: l.parent()
        Power Series Ring in t over Rational Field
    """
    try:
        return x.lift()
    except AttributeError:
        return PowerSeriesRing(Rationals(), "t")(x.list(), x.prec())


def reduce_negative(Q, p, coeffs, offset, exact_form=None):
    """
    Apply cohomology relations to incorporate negative powers of
    `y` into the `y^0` term.

    INPUT:

    - ``p`` -- prime

    - ``Q`` -- cubic polynomial

    - ``coeffs`` -- list of length 3 lists. The
      `i`-th list ``[a, b, c]`` represents
      `y^{2(i - offset)} (a + bx + cx^2) dx/y`.

    - ``offset`` -- nonnegative integer

    OUTPUT:

    The reduction is performed in-place. The output is placed
    in coeffs[offset]. Note that coeffs[i] will be meaningless for i
    offset after this function is finished.

    EXAMPLES::

        sage: R.<x> = Integers(5^3)['x']
        sage: Q = x^3 - x + R(1/4)
        sage: coeffs = [[10, 15, 20], [1, 2, 3], [4, 5, 6], [7, 8, 9]]
        sage: coeffs = [[R.base_ring()(a) for a in row] for row in coeffs]
        sage: monsky_washnitzer.reduce_negative(Q, 5, coeffs, 3)
        sage: coeffs[3]
        [28, 52, 9]

    ::

        sage: R.<x> = Integers(7^3)['x']
        sage: Q = x^3 - x + R(1/4)
        sage: coeffs = [[7, 14, 21], [1, 2, 3], [4, 5, 6], [7, 8, 9]]
        sage: coeffs = [[R.base_ring()(a) for a in row] for row in coeffs]
        sage: monsky_washnitzer.reduce_negative(Q, 7, coeffs, 3)
        sage: coeffs[3]
        [245, 332, 9]
    """

    m = helper_matrix(Q).list()
    base_ring = Q.base_ring()
    next_a = coeffs[0]

    if exact_form is not None:
        x = exact_form.parent().gen(0)
        y = exact_form.parent()(exact_form.parent().base_ring().gen(0))

    try:
        three_j_plus_5 = 5 - base_ring(6 * offset)
        three_j_plus_7 = 7 - base_ring(6 * offset)
        six = base_ring(6)

        for i in range(offset):
            j = 2 * (i - offset)
            a = next_a
            next_a = coeffs[i + 1]

            # todo: the following divisions will sometimes involve
            # a division by (a power of) p. In all cases, we know (from
            # Kedlaya's estimates) that the answer should be p-integral.
            # However, since we're working over $Z/p^k Z$, we're not allowed
            # to "divide by p". So currently we lift to Q, divide, and coerce
            # back. Eventually, when pAdicInteger is implemented, and plays
            # nicely with pAdicField, we should reimplement this stuff
            # using pAdicInteger.

            if p.divides(j + 1):
                # need to lift here to perform the division
                a[0] = base_ring(lift(a[0]) / (j + 1))
                a[1] = base_ring(lift(a[1]) / (j + 1))
                a[2] = base_ring(lift(a[2]) / (j + 1))
            else:
                j_plus_1_inv = ~base_ring(j + 1)
                a[0] = a[0] * j_plus_1_inv
                a[1] = a[1] * j_plus_1_inv
                a[2] = a[2] * j_plus_1_inv

            c1 = m[3] * a[0] + m[4] * a[1] + m[5] * a[2]
            c2 = m[6] * a[0] + m[7] * a[1] + m[8] * a[2]
            next_a[0] = next_a[0] - three_j_plus_5 * c1
            next_a[1] = next_a[1] - three_j_plus_7 * c2

            three_j_plus_7 = three_j_plus_7 + six
            three_j_plus_5 = three_j_plus_5 + six

            if exact_form is not None:
                c0 = m[0] * a[0] + m[1] * a[1] + m[2] * a[2]
                exact_form += (c0 + c1 * x + c2 * x**2) * y ** (j + 1)

    except NotImplementedError:
        raise NotImplementedError(
            "It looks like you've found a "
            "non-integral matrix of Frobenius! "
            f"(Q={Q}, p={p})\nTime to write a paper."
        )

    coeffs[int(offset)] = next_a

    return exact_form


def reduce_positive(Q, p, coeffs, offset, exact_form=None):
    """
    Apply cohomology relations to incorporate positive powers of
    `y` into the `y^0` term.

    INPUT:

    - ``Q`` -- cubic polynomial

    - ``coeffs`` -- list of length 3 lists. The
      `i`-th list [a, b, c] represents
      `y^{2(i - offset)} (a + bx + cx^2) dx/y`.

    - ``offset`` -- nonnegative integer

    OUTPUT:

    The reduction is performed in-place. The output is placed
    in coeffs[offset]. Note that coeffs[i] will be meaningless for i
    offset after this function is finished.

    EXAMPLES::

        sage: R.<x> = Integers(5^3)['x']
        sage: Q = x^3 - x + R(1/4)

    ::

        sage: coeffs = [[1, 2, 3], [10, 15, 20]]
        sage: coeffs = [[R.base_ring()(a) for a in row] for row in coeffs]
        sage: monsky_washnitzer.reduce_positive(Q, 5, coeffs, 0)
        sage: coeffs[0]
        [16, 102, 88]

    ::

        sage: coeffs = [[9, 8, 7], [10, 15, 20]]
        sage: coeffs = [[R.base_ring()(a) for a in row] for row in coeffs]
        sage: monsky_washnitzer.reduce_positive(Q, 5, coeffs, 0)
        sage: coeffs[0]
        [24, 108, 92]
    """

    base_ring = Q.base_ring()
    next_a = coeffs[len(coeffs) - 1]

    Qa = Q[1]
    Qb = Q[0]

    A = 2 * Qa
    B = 3 * Qb

    offset = Integer(offset)

    if exact_form is not None:
        x = exact_form.parent().gen(0)
        y = exact_form.parent().base_ring().gen(0)
    #        y = exact_form.parent()(exact_form.parent().base_ring().gen(0))

    for i in range(len(coeffs) - 1, offset, -1):
        j = 2 * (i - offset) - 2
        a = next_a
        next_a = coeffs[i - 1]

        a[0] = a[0] - Qa * a[2] / 3  # subtract d(y^j + 3)
        if exact_form is not None:
            exact_form += Q.base_ring()(a[2].lift() / (3 * j + 9)) * y ** (j + 3)

        # todo: see comments about pAdicInteger in reduceNegative()

        # subtract off c1 of d(x y^j + 1), and
        if p.divides(3 * j + 5):
            c1 = base_ring(lift(a[0]) / (3 * j + 5))
        else:
            c1 = a[0] / (3 * j + 5)

        # subtract off c2 of d(x^2 y^j + 1)
        if p.divides(3 * j + 7):
            c2 = base_ring(lift(a[1]) / (3 * j + 7))
        else:
            c2 = a[1] / (3 * j + 7)

        next_a[0] = next_a[0] + B * c1 * (j + 1)
        next_a[1] = next_a[1] + A * c1 * (j + 1) + B * c2 * (j + 1)
        next_a[2] = next_a[2] + A * c2 * (j + 1)

        if exact_form is not None:
            exact_form += (c1 * x + c2 * x**2) * y ** (j + 1)

    coeffs[int(offset)] = next_a

    return exact_form


def reduce_zero(Q, coeffs, offset, exact_form=None):
    """
    Apply cohomology relation to incorporate `x^2 y^0` term
    into `x^0 y^0` and `x^1 y^0` terms.

    INPUT:

    - ``Q`` -- cubic polynomial

    - ``coeffs`` -- list of length 3 lists. The
      `i`-th list [a, b, c] represents
      `y^{2(i - offset)} (a + bx + cx^2) dx/y`.

    - ``offset`` -- nonnegative integer

    OUTPUT:

    The reduction is performed in-place. The output is placed
    in coeffs[offset]. This method completely ignores coeffs[i] for i
    != offset.

    EXAMPLES::

        sage: R.<x> = Integers(5^3)['x']
        sage: Q = x^3 - x + R(1/4)
        sage: coeffs = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        sage: coeffs = [[R.base_ring()(a) for a in row] for row in coeffs]
        sage: monsky_washnitzer.reduce_zero(Q, coeffs, 1)
        sage: coeffs[1]
        [6, 5, 0]
    """

    a = coeffs[int(offset)]
    if a[2] == 0:
        return exact_form

    Qa = Q[1]

    a[0] = a[0] - a[2] * Qa / 3  # $3x^2 dx/y = -a dx/y$

    coeffs[int(offset)] = a

    if exact_form is not None:
        y = exact_form.parent()(exact_form.parent().base_ring().gen(0))
        exact_form += Q.base_ring()(a[2] / 3) * y

    a[2] = 0

    coeffs[int(offset)] = a
    return exact_form


def reduce_all(Q, p, coeffs, offset, compute_exact_form=False):
    """
    Apply cohomology relations to reduce all terms to a linear
    combination of `dx/y` and `x dx/y`.

    INPUT:

    - ``Q`` -- cubic polynomial

    - ``coeffs`` -- list of length 3 lists. The
      `i`-th list [a, b, c] represents
      `y^{2(i - offset)} (a + bx + cx^2) dx/y`.

    - ``offset`` -- nonnegative integer

    OUTPUT:

    - ``A``, ``B`` -- pair such that the input differential is
      cohomologous to (A + Bx) dx/y

    .. NOTE::

        The algorithm operates in-place, so the data in coeffs is
        destroyed.

    EXAMPLES::

        sage: R.<x> = Integers(5^3)['x']
        sage: Q = x^3 - x + R(1/4)
        sage: coeffs = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        sage: coeffs = [[R.base_ring()(a) for a in row] for row in coeffs]
        sage: monsky_washnitzer.reduce_all(Q, 5, coeffs, 1)
        (21, 106)
    """

    R = Q.base_ring()

    if compute_exact_form:
        # exact_form = SpecialCubicQuotientRing(Q, laurent_series=True)(0)
        exact_form = PolynomialRing(LaurentSeriesRing(Q.base_ring(), "y"), "x").zero()
        # t = (Q.base_ring().order().factor())[0]
        # from sage.rings.padics.qp import pAdicField
        # exact_form = PolynomialRing(LaurentSeriesRing(pAdicField(p, t[1]), 'y'), 'x')(0)
    else:
        exact_form = None

    while len(coeffs) <= offset:
        coeffs.append([R(0), R(0), R(0)])

    exact_form = reduce_negative(Q, p, coeffs, offset, exact_form)
    exact_form = reduce_positive(Q, p, coeffs, offset, exact_form)
    exact_form = reduce_zero(Q, coeffs, offset, exact_form)

    if exact_form is None:
        return coeffs[int(offset)][0], coeffs[int(offset)][1]
    else:
        return (coeffs[int(offset)][0], coeffs[int(offset)][1]), exact_form


def frobenius_expansion_by_newton(Q, p, M):
    r"""
    Compute the action of Frobenius on `dx/y` and on
    `x dx/y`, using Newton's method (as suggested in Kedlaya's
    paper [Ked2001]_).

    (This function does *not* yet use the cohomology relations - that
    happens afterwards in the "reduction" step.)

    More specifically, it finds `F_0` and `F_1` in
    the quotient ring `R[x, T]/(T - Q(x))`, such that

    .. MATH::

       F(  dx/y) = T^{-r} F0 dx/y, \text{\ and\ } F(x dx/y) = T^{-r} F1 dx/y

    where

    .. MATH::

       r = ( (2M-3)p - 1 )/2.


    (Here `T` is `y^2 = z^{-2}`, and `R` is the
    coefficient ring of `Q`.)

    `F_0` and `F_1` are computed in the
    :class:`SpecialCubicQuotientRing` associated to `Q`, so all powers
    of `x^j` for `j \geq 3` are reduced to powers of
    `T`.

    INPUT:

    - ``Q`` -- cubic polynomial of the form
      `Q(x) = x^3 + ax + b`, whose coefficient ring is a
      `Z/(p^M)Z`-algebra

    - ``p`` -- residue characteristic of the `p`-adic field

    - ``M`` -- `p`-adic precision of the coefficient ring
      (this will be used to determine the number of Newton iterations)

    OUTPUT:

    - ``F0``, ``F1`` -- elements of
      ``SpecialCubicQuotientRing(Q)``, as described above

    - ``r`` -- nonnegative integer, as described above

    EXAMPLES::

        sage: from sage.schemes.hyperelliptic_curves_smooth_model.monsky_washnitzer import frobenius_expansion_by_newton
        sage: R.<x> = Integers(5^3)['x']
        sage: Q = x^3 - x + R(1/4)
        sage: frobenius_expansion_by_newton(Q,5,3)
        ((25*T^5 + 75*T^3 + 100*T^2 + 100*T + 100) + (5*T^6 + 80*T^5 + 100*T^3
        + 25*T + 50)*x + (55*T^5 + 50*T^4 + 75*T^3 + 25*T^2 + 25*T + 25)*x^2,
        (5*T^8 + 15*T^7 + 95*T^6 + 10*T^5 + 25*T^4 + 25*T^3 + 100*T^2 + 50)
        + (65*T^7 + 55*T^6 + 70*T^5 + 100*T^4 + 25*T^2 + 100*T)*x
        + (15*T^6 + 115*T^5 + 75*T^4 + 100*T^3 + 50*T^2 + 75*T + 75)*x^2, 7)
    """

    S = SpecialCubicQuotientRing(Q)
    x, _ = S.gens()  # T = y^2
    base_ring = S.base_ring()

    # When we compute Frob(1/y) we actually only need precision M-1, since
    # we're going to multiply by p at the end anyway.
    M = float(M - 1)

    # Kedlaya sets s = Q(x^p)/T^p = 1 + p T^{-p} E, where
    # E = (Q(x^p) - Q(x)^p) / p (has integral coefficients).
    # Then he computes s^{-1/2} in S, using Newton's method to find
    # successive approximations. We follow this plan, but we normalise our
    # approximations so that we only ever need positive powers of T.

    # Start by setting r = Q(x^p)/2 = 1/2 T^p s.
    # (The 1/2 is for convenience later on.)
    x_to_p_less_one = x ** (p - 1)
    x_to_p = x_to_p_less_one * x
    x_to_p_cubed = x_to_p.square() * x_to_p
    r = (base_ring(1) / base_ring(2)) * (x_to_p_cubed + Q[1] * x_to_p + S(Q[0]))

    # todo: this next loop would be clearer if it used the newton_method_sizes()
    # function

    # We will start with a hard-coded initial approximation, which we provide
    # up to precision 3. First work out what precision is best to start with.
    if M <= 3:
        initial_precision = M
    elif ceil(log(M / 2, 2)) == ceil(log(M / 3, 2)):
        # In this case there is no advantage to starting with precision three,
        # because we'll overshoot at the end. E.g. suppose the final precision
        # is 8. If we start with precision 2, we need two iterations to get us
        # to 8. If we start at precision 3, we will still need two iterations,
        # but we do more work along the way. So may as well start with only 2.
        initial_precision = 2
    else:
        initial_precision = 3

    # Now compute the first approximation. In the main loop below, X is the
    # normalised approximation, and k is the precision. More specifically,
    # X = T^{p(k-1)} x_i, where x_i is an approximation to s^{-1/2}, and the
    # approximation is correct mod p^k.
    if initial_precision == 1:
        k = 1
        X = S(1)
    elif initial_precision == 2:
        # approximation is 3/2 - 1/2 s
        k = 2
        X = S(base_ring(3) / base_ring(2)).shift(p) - r
    elif initial_precision == 3:
        # approximation is (15 - 10 s + 3 s^2) / 8
        k = 3
        X = (base_ring(1) / base_ring(8)) * (
            S(15).shift(2 * p)
            - (base_ring(20) * r).shift(p)
            + (base_ring(12) * r.square())
        )
    # The key to the following calculation is that the T^{-m} coefficient
    # of every x_i is divisible by p^(ceil(m/p)) (for m >= 0). Therefore if
    # we are only expecting an answer correct mod p^k, we can truncate
    # beyond the T^{-(k-1)p} term without any problems.

    # todo: what would be really nice is to be able to work in a lower
    # precision *coefficient ring* when we start the iteration, and move up to
    # higher precision rings as the iteration proceeds. This would be feasible
    # over Integers(p**n), but quite complicated (maybe impossible) over a more
    # general base ring. This might give a decent constant factor speedup;
    # or it might not, depending on how much the last iteration dominates the
    # whole runtime. My guess is that it isn't worth the effort.

    three_halves = base_ring(3) / base_ring(2)

    # Newton iteration loop
    while k < M:
        # target_k = k' = precision we want our answer to be after this iteration
        target_k = 2 * k

        # This prevents us overshooting. For example if the current precision
        # is 3 and we want to get to 10, we're better off going up to 5
        # instead of 6, because it is less work to get from 5 to 10 than it
        # is to get from 6 to 10.
        if ceil(log(M / target_k, 2)) == ceil(log(M / (target_k - 1), 2)):
            target_k -= 1

        # temp = T^{p(3k-2)} 1/2 s x_i^3
        temp = X.square() * (X * r)

        # We know that the final result is only going to be correct mod
        # p^(target_k), so we might as well truncate the extraneous terms now.
        # temp = T^{p(k'-1)} 1/2 s x_i^3
        temp = temp.shift(-p * (3 * k - target_k - 1))

        # X = T^{p(k'-1)} (3/2 x_i - 1/2 s x_i^3)
        #   = T^{p(k'-1)} x_{i+1}
        X = (three_halves * X).shift(p * (target_k - k)) - temp

        k = target_k

    # Now k should equal M, since we're up to the correct precision
    assert k == M, "Oops, something went wrong in the iteration"

    # We should have s^{-1/2} correct to precision M.
    # The following line can be uncommented to verify this.
    # (It is a slow verification though, can double the whole computation time.)

    # assert (p * X.square() * r * base_ring(2)).coeffs() == \
    #       R(p).shift(p*(2*M - 1)).coeffs()

    # Finally incorporate frobenius of dx and x dx, and choose offset that
    # compensates for our normalisations by powers of T.
    F0 = base_ring(p) * x_to_p_less_one * X
    F1 = F0 * x_to_p
    offset = ((2 * k - 1) * p - 1) / 2

    return F0, F1, offset


def frobenius_expansion_by_series(Q, p, M):
    r"""
    Compute the action of Frobenius on `dx/y` and on `x dx/y`, using a
    series expansion.

    (This function computes the same thing as
    frobenius_expansion_by_newton(), using a different method.
    Theoretically the Newton method should be asymptotically faster,
    when the precision gets large. However, in practice, this functions
    seems to be marginally faster for moderate precision, so I'm
    keeping it here until I figure out exactly why it is faster.)

    (This function does *not* yet use the cohomology relations - that
    happens afterwards in the "reduction" step.)

    More specifically, it finds F0 and F1 in the quotient ring
    `R[x, T]/(T - Q(x))`, such that
    `F(  dx/y) = T^{-r} F0 dx/y`, and
    `F(x dx/y) = T^{-r} F1 dx/y` where
    `r = ( (2M-3)p - 1 )/2`. (Here `T` is `y^2 = z^{-2}`,
    and `R` is the coefficient ring of `Q`.)

    `F_0` and `F_1` are computed in the
    :class:`SpecialCubicQuotientRing` associated to `Q`, so all powers
    of `x^j` for `j \geq 3` are reduced to powers of
    `T`.

    It uses the sum

    .. MATH::

         F0 = \sum_{k=0}^{M-2} \binom{-1/2}{k} p x^{p-1} E^k T^{(M-2-k)p}

    and

    .. MATH::

         F1 = x^p F0,

     where `E = Q(x^p) - Q(x)^p`.

    INPUT:

    - ``Q`` -- cubic polynomial of the form
      `Q(x) = x^3 + ax + b`, whose coefficient ring is a
      `\ZZ/(p^M)\ZZ` -algebra

    - ``p`` -- residue characteristic of the `p`-adic field

    - ``M`` -- `p`-adic precision of the coefficient ring
      (this will be used to determine the number of terms in the
      series)

    OUTPUT:

    - ``F0``, ``F1`` -- elements of
      ``SpecialCubicQuotientRing(Q)``, as described above

    - ``r`` -- nonnegative integer, as described above

    EXAMPLES::

        sage: from sage.schemes.hyperelliptic_curves_smooth_model.monsky_washnitzer import frobenius_expansion_by_series
        sage: R.<x> = Integers(5^3)['x']
        sage: Q = x^3 - x + R(1/4)
        sage: frobenius_expansion_by_series(Q,5,3)                                      # needs sage.libs.pari
        ((25*T^5 + 75*T^3 + 100*T^2 + 100*T + 100) + (5*T^6 + 80*T^5 + 100*T^3
        + 25*T + 50)*x + (55*T^5 + 50*T^4 + 75*T^3 + 25*T^2 + 25*T + 25)*x^2,
        (5*T^8 + 15*T^7 + 95*T^6 + 10*T^5 + 25*T^4 + 25*T^3 + 100*T^2 + 50)
        + (65*T^7 + 55*T^6 + 70*T^5 + 100*T^4 + 25*T^2 + 100*T)*x
        + (15*T^6 + 115*T^5 + 75*T^4 + 100*T^3 + 50*T^2 + 75*T + 75)*x^2, 7)
    """

    S = SpecialCubicQuotientRing(Q)
    x, _ = S.gens()
    base_ring = S.base_ring()

    x_to_p_less_1 = x ** (p - 1)
    x_to_p = x_to_p_less_1 * x

    # compute frobQ = Q(x^p)
    x_to_p_squared = x_to_p * x_to_p
    x_to_p_cubed = x_to_p_squared * x_to_p
    frobQ = x_to_p_cubed + Q[1] * x_to_p + Q[0] * S.one()
    # anticipating the day when p = 3 is supported:
    # frobQ = x_to_p_cubed + Q[2]*x_to_p_squared + Q[1]*x_to_p + Q[0]*S(1)

    E = frobQ - S.one().shift(p)  # E =  Q(x^p) - Q(x)^p

    offset = int(((2 * M - 3) * p - 1) / 2)
    term = p * x_to_p_less_1
    F0 = term.shift((M - 2) * p)

    # todo: Possible speedup idea, perhaps by a factor of 2, but
    # it requires a lot of work:
    # Note that p divides E, so p^k divides E^k. So when we are
    # working with high powers of E, we're doing a lot more work
    # in the multiplications than we need to. To take advantage of
    # this we would need some protocol for "lowering the precision"
    # of a SpecialCubicQuotientRing. This would be quite messy to
    # do properly over an arbitrary base ring. Perhaps it is
    # feasible to do for the most common case (i.e. Z/p^nZ).
    # (but it probably won't save much time unless p^n is very
    # large, because the machine word size is probably pretty
    # big anyway.)

    for k in range(1, int(M - 1)):
        term = term * E
        c = base_ring(binomial(QQ((-1, 2)), k))
        F0 += (term * c).shift((M - k - 2) * p)

    return F0, F0 * x_to_p, offset


def adjusted_prec(p, prec):
    r"""
    Compute how much precision is required in ``matrix_of_frobenius`` to
    get an answer correct to ``prec`` `p`-adic digits.

    The issue is that the algorithm used in
    :func:`matrix_of_frobenius` sometimes performs divisions by `p`,
    so precision is lost during the algorithm.

    The estimate returned by this function is based on Kedlaya's result
    (Lemmas 2 and 3 of [Ked2001]_),
    which implies that if we start with `M` `p`-adic
    digits, the total precision loss is at most
    `1 + \lfloor \log_p(2M-3) \rfloor` `p`-adic
    digits. (This estimate is somewhat less than the amount you would
    expect by naively counting the number of divisions by
    `p`.)

    INPUT:

    - ``p`` -- a prime ``p >= 5``

    - ``prec`` -- integer; desired output precision, ``prec >= 1``

    OUTPUT: adjusted precision (usually slightly more than ``prec``)

    EXAMPLES::

        sage: from sage.schemes.hyperelliptic_curves_smooth_model.monsky_washnitzer import adjusted_prec
        sage: adjusted_prec(5,2)
        3
    """
    # initial estimate:
    if prec <= 2:
        defect = 0
        adjusted = 2
    else:
        defect = Integer(2 * prec - 3).exact_log(p)
        adjusted = prec + defect - 1

    # increase it until we have enough
    while adjusted - defect - 1 < prec:
        adjusted += 1

    return adjusted


def matrix_of_frobenius(Q, p, M, trace=None, compute_exact_forms=False):
    r"""
    Compute the matrix of Frobenius on Monsky-Washnitzer cohomology,
    with respect to the basis `(dx/y, x dx/y)`.

    INPUT:

    - ``Q`` -- cubic polynomial `Q(x) = x^3 + ax + b`
      defining an elliptic curve `E` by
      `y^2 = Q(x)`. The coefficient ring of `Q` should be a
      `\ZZ/(p^M)\ZZ`-algebra in which the matrix of
      frobenius will be constructed.

    - ``p`` -- prime >= 5 for which E has good reduction

    - ``M`` -- integer >= 2; `p` -adic precision of the coefficient ring

    - ``trace`` -- (optional) the trace of the matrix, if
      known in advance. This is easy to compute because it is just the
      `a_p` of the curve. If the trace is supplied,
      matrix_of_frobenius will use it to speed the computation (i.e. we
      know the determinant is `p`, so we have two conditions, so
      really only column of the matrix needs to be computed. it is
      actually a little more complicated than that, but that's the basic
      idea.) If trace=None, then both columns will be computed
      independently, and you can get a strong indication of correctness
      by verifying the trace afterwards.

      .. WARNING::

          THE RESULT WILL NOT NECESSARILY BE CORRECT TO M p-ADIC
          DIGITS. If you want prec digits of precision, you need to use
          the function adjusted_prec(), and then you need to reduce the
          answer mod `p^{\mathrm{prec}}` at the end.

    OUTPUT:

    `2 \times 2` matrix of Frobenius acting on Monsky-Washnitzer cohomology,
    with entries in the coefficient ring of ``Q``.

    EXAMPLES:

    A simple example::

        sage: p = 5
        sage: prec = 3
        sage: M = monsky_washnitzer.adjusted_prec(p, prec); M
        4
        sage: R.<x> = PolynomialRing(Integers(p**M))
        sage: A = monsky_washnitzer.matrix_of_frobenius(x^3 - x + R(1/4), p, M)
        sage: A
        [340  62]
        [ 70 533]

    But the result is only accurate to ``prec`` digits::

        sage: B = A.change_ring(Integers(p**prec))
        sage: B
        [90 62]
        [70 33]

    Check trace (123 = -2 mod 125) and determinant::

        sage: B.det()
        5
        sage: B.trace()
        123
        sage: EllipticCurve([-1, 1/4]).ap(5)
        -2

    Try using the trace to speed up the calculation::

        sage: A = monsky_washnitzer.matrix_of_frobenius(x^3 - x + R(1/4),
        ....:                                           p, M, -2)
        sage: A
        [ 90  62]
        [320 533]

    Hmmm... it looks different, but that's because the trace of our
    first answer was only -2 modulo `5^3`, not -2 modulo
    `5^5`. So the right answer is::

        sage: A.change_ring(Integers(p**prec))
        [90 62]
        [70 33]

    Check it works with only one digit of precision::

        sage: p = 5
        sage: prec = 1
        sage: M = monsky_washnitzer.adjusted_prec(p, prec)
        sage: R.<x> = PolynomialRing(Integers(p**M))
        sage: A = monsky_washnitzer.matrix_of_frobenius(x^3 - x + R(1/4), p, M)
        sage: A.change_ring(Integers(p))
        [0 2]
        [0 3]

    Here is an example that is particularly badly conditioned for
    using the trace trick::

        sage: # needs sage.libs.pari
        sage: p = 11
        sage: prec = 3
        sage: M = monsky_washnitzer.adjusted_prec(p, prec)
        sage: R.<x> = PolynomialRing(Integers(p**M))
        sage: A = monsky_washnitzer.matrix_of_frobenius(x^3 + 7*x + 8, p, M)
        sage: A.change_ring(Integers(p**prec))
        [1144  176]
        [ 847  185]

    The problem here is that the top-right entry is divisible by 11,
    and the bottom-left entry is divisible by `11^2`. So when
    you apply the trace trick, neither `F(dx/y)` nor
    `F(x dx/y)` is enough to compute the whole matrix to the
    desired precision, even if you try increasing the target precision
    by one. Nevertheless, ``matrix_of_frobenius`` knows
    how to get the right answer by evaluating `F((x+1) dx/y)`
    instead::

        sage: A = monsky_washnitzer.matrix_of_frobenius(x^3 + 7*x + 8, p, M, -2)
        sage: A.change_ring(Integers(p**prec))
        [1144  176]
        [ 847  185]

    The running time is about ``O(p*prec**2)`` (times some logarithmic
    factors), so it is feasible to run on fairly large primes, or
    precision (or both?!?!)::

        sage: # long time, needs sage.libs.pari
        sage: p = 10007
        sage: prec = 2
        sage: M = monsky_washnitzer.adjusted_prec(p, prec)
        sage: R.<x> = PolynomialRing(Integers(p**M))
        sage: A = monsky_washnitzer.matrix_of_frobenius(x^3 - x + R(1/4), p, M)
        sage: B = A.change_ring(Integers(p**prec)); B
        [74311982 57996908]
        [95877067 25828133]
        sage: B.det()
        10007
        sage: B.trace()
        66
        sage: EllipticCurve([-1, 1/4]).ap(10007)
        66

    ::

        sage: # long time, needs sage.libs.pari
        sage: p = 5
        sage: prec = 300
        sage: M = monsky_washnitzer.adjusted_prec(p, prec)
        sage: R.<x> = PolynomialRing(Integers(p**M))
        sage: A = monsky_washnitzer.matrix_of_frobenius(x^3 - x + R(1/4), p, M)
        sage: B = A.change_ring(Integers(p**prec))
        sage: B.det()
        5
        sage: -B.trace()
        2
        sage: EllipticCurve([-1, 1/4]).ap(5)
        -2

    Let us check consistency of the results for a range of precisions::

        sage: # long time, needs sage.libs.pari
        sage: p = 5
        sage: max_prec = 60
        sage: M = monsky_washnitzer.adjusted_prec(p, max_prec)
        sage: R.<x> = PolynomialRing(Integers(p**M))
        sage: A = monsky_washnitzer.matrix_of_frobenius(x^3 - x + R(1/4), p, M)
        sage: A = A.change_ring(Integers(p**max_prec))
        sage: result = []
        sage: for prec in range(1, max_prec):
        ....:     M = monsky_washnitzer.adjusted_prec(p, prec)
        ....:     R.<x> = PolynomialRing(Integers(p^M),'x')
        ....:     B = monsky_washnitzer.matrix_of_frobenius(x^3 - x + R(1/4), p, M)
        ....:     B = B.change_ring(Integers(p**prec))
        ....:     result.append(B == A.change_ring(Integers(p**prec)))
        sage: result == [True] * (max_prec - 1)
        True

    The remaining examples discuss what happens when you take the
    coefficient ring to be a power series ring; i.e. in effect you're
    looking at a family of curves.

    The code does in fact work...

    ::

        sage: # needs sage.libs.pari
        sage: p = 11
        sage: prec = 3
        sage: M = monsky_washnitzer.adjusted_prec(p, prec)
        sage: S.<t> = PowerSeriesRing(Integers(p**M), default_prec=4)
        sage: a = 7 + t + 3*t^2
        sage: b = 8 - 6*t + 17*t^2
        sage: R.<x> = PolynomialRing(S)
        sage: Q = x**3 + a*x + b
        sage: A = monsky_washnitzer.matrix_of_frobenius(Q, p, M)            # long time
        sage: B = A.change_ring(PowerSeriesRing(Integers(p**prec), 't',     # long time
        ....:                                   default_prec=4)); B
        [1144 + 264*t + 841*t^2 + 1025*t^3 + O(t^4)  176 + 1052*t + 216*t^2 + 523*t^3 + O(t^4)]
        [   847 + 668*t + 81*t^2 + 424*t^3 + O(t^4)   185 + 341*t + 171*t^2 + 642*t^3 + O(t^4)]

    The trace trick should work for power series rings too, even in the
    badly-conditioned case. Unfortunately I do not know how to compute
    the trace in advance, so I am not sure exactly how this would help.
    Also, I suspect the running time will be dominated by the
    expansion, so the trace trick will not really speed things up anyway.
    Another problem is that the determinant is not always p::

        sage: B.det()                                               # long time
        11 + 484*t^2 + 451*t^3 + O(t^4)

    However, it appears that the determinant always has the property
    that if you substitute t - 11t, you do get the constant series p
    (mod p\*\*prec). Similarly for the trace. And since the parameter
    only really makes sense when it is divisible by p anyway, perhaps
    this is not a problem after all.
    """
    M = int(M)
    if M < 2:
        raise ValueError("M (=%s) must be at least 2" % M)

    base_ring = Q.base_ring()

    # Expand out frobenius of dx/y and x dx/y.
    # (You can substitute frobenius_expansion_by_series here, that will work
    # as well. See its docstring for some performance notes.)
    F0, F1, offset = frobenius_expansion_by_newton(Q, p, M)
    # F0, F1, offset = frobenius_expansion_by_series(Q, p, M)

    if compute_exact_forms:
        # we need to do all the work to get the exact expressions f such that F(x^i dx/y) = df + \sum a_i x^i dx/y
        F0_coeffs = transpose_list(F0.coeffs())
        F0_reduced, f_0 = reduce_all(Q, p, F0_coeffs, offset, True)

        F1_coeffs = transpose_list(F1.coeffs())
        F1_reduced, f_1 = reduce_all(Q, p, F1_coeffs, offset, True)

    elif M == 2:
        # This implies that only one digit of precision is valid, so we only need
        # to reduce the second column. Also, the trace doesn't help at all.

        F0_reduced = [base_ring(0), base_ring(0)]

        F1_coeffs = transpose_list(F1.coeffs())
        F1_reduced = reduce_all(Q, p, F1_coeffs, offset)

    elif trace is None:
        # No trace provided, just reduce F(dx/y) and F(x dx/y) separately.

        F0_coeffs = transpose_list(F0.coeffs())
        F0_reduced = reduce_all(Q, p, F0_coeffs, offset)

        F1_coeffs = transpose_list(F1.coeffs())
        F1_reduced = reduce_all(Q, p, F1_coeffs, offset)

    else:
        # Trace has been provided.

        # In most cases this can be used to quickly compute F(dx/y) from
        # F(x dx/y). However, if we're unlucky, the (dx/y)-component of
        # F(x dx/y) (i.e. the top-right corner of the matrix) may be divisible
        # by p, in which case there isn't enough information to get the
        # (x dx/y)-component of F(dx/y) to the desired precision. When this
        # happens, it turns out that F((x+1) dx/y) always *does* give enough
        # information (together with the trace) to get both columns to the
        # desired precision.

        # First however we need a quick way of telling whether the top-right
        # corner is divisible by p, i.e. we want to compute the second column
        # of the matrix mod p. We could do this by just running the entire
        # algorithm with M = 2 (which assures precision 1). Luckily, we've
        # already done most of the work by computing F1 to high precision; so
        # all we need to do is extract the coefficients that would correspond
        # to the first term of the series, and run the reduction on them.

        # todo: actually we only need to do this reduction step mod p^2, not
        # mod p^M, which is what the code currently does. If the base ring
        # is Integers(p^M), then it is easy. Otherwise it is tricky to construct
        # the right ring, I don't know how to do it.

        F1_coeffs = transpose_list(F1.coeffs())
        F1_modp_coeffs = F1_coeffs[int((M - 2) * p) :]
        # make a copy, because reduce_all will destroy the coefficients:
        F1_modp_coeffs = [list(row) for row in F1_modp_coeffs]
        F1_modp_offset = offset - (M - 2) * p
        F1_modp_reduced = reduce_all(Q, p, F1_modp_coeffs, F1_modp_offset)

        if F1_modp_reduced[0].is_unit():
            # If the first entry is invertible mod p, then F(x dx/y) is sufficient
            # to get the whole matrix.

            F1_reduced = reduce_all(Q, p, F1_coeffs, offset)

            F0_reduced = [base_ring(trace) - F1_reduced[1], None]
            # using that the determinant is p:
            F0_reduced[1] = (F0_reduced[0] * F1_reduced[1] - base_ring(p)) / F1_reduced[
                0
            ]

        else:
            # If the first entry is zero mod p, then F((x+1) dx/y) will be sufficient
            # to get the whole matrix. (Here we are using the fact that the second
            # entry *cannot* be zero mod p. This is guaranteed by some results in
            # section 3.2 of ``Computation of p-adic Heights and Log Convergence''
            # by Mazur, Stein, Tate. But let's quickly check it anyway :-))
            msg = "The second entry in the second column "
            msg += "should be invertible mod p!"
            assert F1_modp_reduced[1].is_unit(), msg

            G0_coeffs = transpose_list((F0 + F1).coeffs())
            G0_reduced = reduce_all(Q, p, G0_coeffs, offset)

            # Now G0_reduced expresses F((x+1) dx/y) in terms of dx/y and x dx/y.
            # Re-express this in terms of (x+1) dx/y and x dx/y.
            H0_reduced = [G0_reduced[0], G0_reduced[1] - G0_reduced[0]]

            # The thing we're about to divide by better be a unit.
            msg = "The second entry in this column "
            msg += "should be invertible mod p!"
            assert H0_reduced[1].is_unit(), msg

            # Figure out the second column using the trace...
            H1_reduced = [None, base_ring(trace) - H0_reduced[0]]
            # ... and using that the determinant is p:
            H1_reduced[0] = (H0_reduced[0] * H1_reduced[1] - base_ring(p)) / H0_reduced[
                1
            ]

            # Finally, change back to the usual basis (dx/y, x dx/y)
            F1_reduced = [H1_reduced[0], H1_reduced[0] + H1_reduced[1]]
            F0_reduced = [
                H0_reduced[0] - F1_reduced[0],
                H0_reduced[0] + H0_reduced[1] - F1_reduced[1],
            ]

        # One more sanity check: our final result should be congruent mod p
        # to the approximation we used earlier.
        msg = "The output matrix is not congruent mod p "
        msg += "to the approximation found earlier!"
        assert not (
            (F1_reduced[0] - F1_modp_reduced[0]).is_unit()
            or (F1_reduced[1] - F1_modp_reduced[1]).is_unit()
            or F0_reduced[0].is_unit()
            or F0_reduced[1].is_unit()
        ), msg

    if compute_exact_forms:
        return (
            matrix(
                base_ring,
                2,
                2,
                [F0_reduced[0], F1_reduced[0], F0_reduced[1], F1_reduced[1]],
            ),
            f_0,
            f_1,
        )
    else:
        return matrix(
            base_ring,
            2,
            2,
            [F0_reduced[0], F1_reduced[0], F0_reduced[1], F1_reduced[1]],
        )


# ****************************************************************************
# This is a generalization of the above functionality for hyperelliptic curves.
#
# THIS IS A WORK IN PROGRESS.
#
# I tried to embed must stuff into the rings themselves rather than
# just extract and manipulate lists of coefficients. Hence the implementations
# below are much less optimized, so are much slower, but should hopefully be
# easier to follow. (E.g. one can print/make sense of intermediate results.)
#
# AUTHOR:
#    -- Robert Bradshaw (2007-04)
#
# ****************************************************************************


def matrix_of_frobenius_hyperelliptic(Q, p=None, prec=None, M=None):
    r"""
    Compute the matrix of Frobenius on Monsky-Washnitzer cohomology,
    with respect to the basis `(dx/2y, x dx/2y, ...x^{d-2} dx/2y)`, where
    `d` is the degree of `Q`.

    INPUT:

    - ``Q`` -- monic polynomial `Q(x)`

    - ``p`` -- prime `\geq 5` for which `E` has good reduction

    - ``prec`` -- (optional) `p`-adic precision of the coefficient ring

    - ``M`` -- (optional) adjusted `p`-adic precision of the coefficient ring

    OUTPUT:

    `(d-1)` x `(d-1)` matrix `M` of Frobenius on Monsky-Washnitzer cohomology,
    and list of differentials \{f_i \} such that

    .. MATH::

        \phi^* (x^i dx/2y) = df_i + M[i]*vec(dx/2y, ..., x^{d-2} dx/2y)

    EXAMPLES::

        sage: # needs sage.rings.padics
        sage: p = 5
        sage: prec = 3
        sage: R.<x> = QQ['x']
        sage: A,f = monsky_washnitzer.matrix_of_frobenius_hyperelliptic(x^5 - 2*x + 3, p, prec)
        sage: A
        [            4*5 + O(5^3)       5 + 2*5^2 + O(5^3) 2 + 3*5 + 2*5^2 + O(5^3)     2 + 5 + 5^2 + O(5^3)]
        [      3*5 + 5^2 + O(5^3)             3*5 + O(5^3)             4*5 + O(5^3)         2 + 5^2 + O(5^3)]
        [    4*5 + 4*5^2 + O(5^3)     3*5 + 2*5^2 + O(5^3)       5 + 3*5^2 + O(5^3)     2*5 + 2*5^2 + O(5^3)]
        [            5^2 + O(5^3)       5 + 4*5^2 + O(5^3)     4*5 + 3*5^2 + O(5^3)             2*5 + O(5^3)]
    """
    prof = Profiler()
    prof("setup")
    if p is None:
        try:
            K = Q.base_ring()
            p = K.prime()
            prec = K.precision_cap()
        except AttributeError:
            raise ValueError(
                "p and prec must be specified if Q is not " "defined over a p-adic ring"
            )
    if M is None:
        M = adjusted_prec(p, prec)
    extra_prec_ring = Integers(p**M)
    #    extra_prec_ring = pAdicField(p, M) # SLOW!

    real_prec_ring = pAdicField(p, prec)  # pAdicField(p, prec) # To capped absolute?
    S = SpecialHyperellipticQuotientRing(Q, extra_prec_ring, True)
    MW = S.monsky_washnitzer()
    prof("frob basis elements")
    F = MW.frob_basis_elements(M, p)

    prof("rationalize")
    # do reduction over Q in case we have non-integral entries (and it is so much faster than padics)
    rational_S = S.change_ring(QQ)
    # this is a hack until pAdics are fast
    # (They are in the latest development bundle, but its not standard and I'd need to merge.
    # (it will periodically cast into this ring to reduce coefficient size)
    rational_S._prec_cap = p**M
    rational_S._p = p
    #    S._p = p
    #    rational_S(F[0]).reduce_fast()
    #    prof("reduce others")

    #    rational_S = S.change_ring(pAdicField(p, M))
    F = [rational_S(F_i) for F_i in F]

    prof("reduce")
    reduced = [F_i.reduce_fast(True) for F_i in F]

    # but the coeffs are WAY more precision than they need to be

    prof("make matrix")
    # now take care of precision capping
    M = matrix(real_prec_ring, [a for f, a in reduced])
    for i in range(M.ncols()):
        for j in range(M.nrows()):
            M[i, j] = M[i, j].add_bigoh(prec)
    return M.transpose(), [f for f, a in reduced]


class SpecialHyperellipticQuotientElement(ModuleElement):
    r"""
    Element in the Hyperelliptic quotient ring.

    EXAMPLES::

        sage: R.<x> = QQ['x']
        sage: E = HyperellipticCurveSmoothModel(x^5 - 36*x + 1)
        sage: x,y = E.monsky_washnitzer_gens()
        sage: MW = x.parent()
        sage: MW(x + x**2 + y - 77)
        -(77-y)*1 + x + x^2
    """

    def __init__(self, parent, val=0, offset=0, check=True):
        """
        Elements in the Hyperelliptic quotient ring.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 36*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: MW = x.parent()
            sage: elt = MW(x + x**2 + y - 77)
            sage: TestSuite(elt).run()
        """
        ModuleElement.__init__(self, parent)
        if not check:
            self._f = parent._poly_ring(val, check=False)
            return
        if isinstance(val, SpecialHyperellipticQuotientElement):
            R = parent.base_ring()
            self._f = parent._poly_ring([a.change_ring(R) for a in val._f])
            return
        if isinstance(val, tuple):
            val, offset = val
        if isinstance(val, list) and val and isinstance(val[0], FreeModuleElement):
            val = transpose_list(val)
        self._f = parent._poly_ring(val)
        if offset != 0:
            self._f = self._f.parent()([a << offset for a in self._f], check=False)

    def _richcmp_(self, other, op) -> bool:
        """
        Compare the elements.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 36*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: x == x
            True
            sage: x > y
            True
        """
        return richcmp(self._f, other._f, op)

    def change_ring(self, R):
        """
        Return the same element after changing the base ring to `R`.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 36*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: MW = x.parent()
            sage: z = MW(x + x**2 + y - 77)
            sage: z.change_ring(AA).parent()                                            # needs sage.rings.number_field
            SpecialHyperellipticQuotientRing K[x,y,y^-1] / (y^2 = x^5 - 36*x + 1)
            over Algebraic Real Field
        """
        return self.parent().change_ring(R)(self)

    def __call__(self, *x):
        """
        Evaluate ``self`` at given arguments.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 36*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: MW = x.parent()
            sage: z = MW(x + x**2 + y - 77); z
            -(77-y)*1 + x + x^2
            sage: z(66)
            4345 + y
            sage: z(5,4)
            -43
        """
        return self._f(*x)

    def __invert__(self):
        """
        Return the inverse of ``self``.

        The general element in our ring is not invertible, but `y` may
        be. We do not want to pass to the fraction field.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 36*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: MW = x.parent()
            sage: z = y**(-1)  # indirect doctest
            sage: z.parent()
            SpecialHyperellipticQuotientRing K[x,y,y^-1] / (y^2 = x^5 - 36*x + 1)
            over Rational Field

            sage: z = (x+y)**(-1)  # indirect doctest
            Traceback (most recent call last):
            ...
            ZeroDivisionError: element not invertible
        """
        if self._f.degree() == 0 and self._f[0].is_unit():
            P = self.parent()
            return P.element_class(P, ~self._f[0])
        else:
            raise ZeroDivisionError("element not invertible")

    def __bool__(self):
        """
        Return ``True`` iff ``self`` is not zero.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: bool(x)
            True
        """
        return bool(self._f)

    def __eq__(self, other):
        """
        Return ``True`` iff ``self`` is equal to ``other``.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: x == y  # indirect doctest
            False
        """
        if not isinstance(other, SpecialHyperellipticQuotientElement):
            other = self.parent()(other)
        return self._f == other._f

    def _add_(self, other):
        """
        Return the sum of two elements.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 36*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: x + y
            y*1 + x
        """
        P = self.parent()
        return P.element_class(P, self._f + other._f)

    def _sub_(self, other):
        """
        Return the difference of two elements.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 36*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: y - x
            y*1 - x
        """
        P = self.parent()
        return P.element_class(P, self._f - other._f)

    def _mul_(self, other):
        """
        Return the product of two elements.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 36*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: y*x
            y*x
        """
        # over Laurent series, addition and subtraction can be
        # expensive, and the degree of this poly is small enough that
        # Karatsuba actually hurts significantly in some cases
        if self._f[0].valuation() + other._f[0].valuation() > -200:
            prod = self._f._mul_generic(other._f)
        else:
            prod = self._f * other._f
        v = prod.list()
        parent = self.parent()
        Q_coeffs = parent._Q_coeffs
        n = len(Q_coeffs) - 1
        y2 = self.parent()._series_ring_y << 1
        for i in range(len(v) - 1, n - 1, -1):
            for j in range(n):
                v[i - n + j] -= Q_coeffs[j] * v[i]
            v[i - n] += y2 * v[i]
        P = self.parent()
        return P.element_class(P, v[0:n])

    def _rmul_(self, c):
        """
        Return the product of ``self`` with `c` on the left.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: x._rmul_(y)                                                           # needs sage.rings.real_interval_field
            y*1*x
        """
        P = self.parent()
        if not c:
            return P.zero()
        ret = [c * a for a in self._f.list(copy=False)]
        while ret and not ret[-1]:  # strip off trailing 0s
            ret.pop()
        return P.element_class(P, ret, check=False)

    def _lmul_(self, c):
        """
        Return the product of ``self`` with `c` on the right.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5-3*x+1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: x._lmul_(y)                                                           # needs sage.rings.real_interval_field
            y*1*x
        """
        P = self.parent()
        if not c:
            return P.zero()
        ret = [a * c for a in self._f.list(copy=False)]
        while ret and not ret[-1]:  # strip off trailing 0s
            ret.pop()
        return P.element_class(P, ret, check=False)

    def __lshift__(self, k):
        """
        Return the left shift of ``self`` by `k`.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: x.__lshift__(3)
            y^3*x
        """
        coeffs = self._f.list(copy=False)
        P = self.parent()
        return P.element_class(P, [a << k for a in coeffs], check=False)

    def __rshift__(self, k):
        """
        Return the right shift of ``self`` by `k`.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: y.__rshift__(3)
            (y^-2)*1
        """
        coeffs = self._f.list(copy=False)
        P = self.parent()
        return P.element_class(P, [a >> k for a in coeffs], check=False)

    def truncate_neg(self, n):
        """
        Return ``self`` minus its terms of degree less than `n` wrt `y`.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: (x + 3*y + 7*x*2*y**4).truncate_neg(1)
            3*y*1 + 14*y^4*x
        """
        coeffs = self._f.list(copy=False)
        P = self.parent()
        return P.element_class(P, [a.truncate_neg(n) for a in coeffs], check=False)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: (x + 3*y)._repr_()
            '3*y*1 + x'
        """
        x = PolynomialRing(QQ, "x").gen(0)
        coeffs = self._f.list()
        return repr_lincomb([(x**i, coeffs[i]) for i in range(len(coeffs))])

    def _latex_(self):
        r"""
        Return a LateX string for ``self``.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: (x + 3*y)._latex_()
            '3y 1 + x'
        """
        x = PolynomialRing(QQ, "x").gen(0)
        coeffs = self._f.list()
        return repr_lincomb(
            [(x**i, coeffs[i]) for i in range(len(coeffs))], is_latex=True
        )

    def diff(self):
        """
        Return the differential of ``self``.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: (x + 3*y).diff()
            (-(9-2*y)*1 + 15*x^4) dx/2y
        """
        #        try:
        #            return self._diff_x
        #        except AttributeError:
        #            pass

        # d(self) = A dx + B dy
        #         = (2y A + BQ') dx/2y
        P = self.parent()
        R = P.base_ring()
        v = self._f.list()
        n = len(v)
        A = P([R(i) * v[i] for i in range(1, n)])
        B = P([a.derivative() for a in v])
        dQ = P._dQ
        return P._monsky_washnitzer((R(2) * A << 1) + dQ * B)

    #        self._diff = self.parent()._monsky_washnitzer(two_y * A + dQ * B)
    #        return self._diff

    def extract_pow_y(self, k):
        r"""
        Return the coefficients of `y^k` in ``self`` as a list.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: (x + 3*y + 9*x*y).extract_pow_y(1)
            [3, 9, 0, 0, 0]
        """
        v = [a[k] for a in self._f.list()]
        v += [ZZ.zero()] * (self.parent()._n - len(v))
        return v

    def min_pow_y(self):
        r"""
        Return the minimal degree of ``self`` with respect to `y`.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: (x + 3*y).min_pow_y()
            0
        """
        if self._f.degree() == -1:
            return ZZ.zero()
        return min([a.valuation() for a in self._f.list()])

    def max_pow_y(self):
        r"""
        Return the maximal degree of ``self`` with respect to `y`.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: (x + 3*y).max_pow_y()
            1
        """
        if self._f.degree() == -1:
            return ZZ.zero()
        return max([a.degree() for a in self._f.list()])

    def coeffs(self, R=None):
        r"""
        Return the raw coefficients of this element.

        INPUT:

        - ``R`` -- an (optional) base-ring in which to cast the coefficients

        OUTPUT:

        - ``coeffs`` -- list of coefficients of powers of `x` for each power
          of `y`

        - ``n`` -- an offset indicating the power of `y` of the first list
          element

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: x.coeffs()
            ([(0, 1, 0, 0, 0)], 0)
            sage: y.coeffs()
            ([(0, 0, 0, 0, 0), (1, 0, 0, 0, 0)], 0)

            sage: a = sum(n*x^n for n in range(5)); a
            x + 2*x^2 + 3*x^3 + 4*x^4
            sage: a.coeffs()
            ([(0, 1, 2, 3, 4)], 0)
            sage: a.coeffs(Qp(7))                                                       # needs sage.rings.padics
            ([(0, 1 + O(7^20), 2 + O(7^20), 3 + O(7^20), 4 + O(7^20))], 0)
            sage: (a*y).coeffs()
            ([(0, 0, 0, 0, 0), (0, 1, 2, 3, 4)], 0)
            sage: (a*y^-2).coeffs()
            ([(0, 1, 2, 3, 4), (0, 0, 0, 0, 0), (0, 0, 0, 0, 0)], -2)

        Note that the coefficient list is transposed compared to how they
        are stored and printed::

            sage: a*y^-2
            (y^-2)*x + (2*y^-2)*x^2 + (3*y^-2)*x^3 + (4*y^-2)*x^4

        A more complicated example::

            sage: a = x^20*y^-3 - x^11*y^2; a
            (y^-3-4*y^-1+6*y-4*y^3+y^5)*1 - (12*y^-3-36*y^-1+36*y+y^2-12*y^3-2*y^4+y^6)*x
             + (54*y^-3-108*y^-1+54*y+6*y^2-6*y^4)*x^2 - (108*y^-3-108*y^-1+9*y^2)*x^3
             + (81*y^-3)*x^4
            sage: raw, offset = a.coeffs()
            sage: a.min_pow_y()
            -3
            sage: offset
            -3
            sage: raw
            [(1, -12, 54, -108, 81),
             (0, 0, 0, 0, 0),
             (-4, 36, -108, 108, 0),
             (0, 0, 0, 0, 0),
             (6, -36, 54, 0, 0),
             (0, -1, 6, -9, 0),
             (-4, 12, 0, 0, 0),
             (0, 2, -6, 0, 0),
             (1, 0, 0, 0, 0),
             (0, -1, 0, 0, 0)]
            sage: sum(c * x^i * y^(j+offset)
            ....:     for j, L in enumerate(raw) for i, c in enumerate(L)) == a
            True

        Can also be used to construct elements::

            sage: a.parent()(raw, offset) == a
            True
        """
        zero = self.base_ring()(0) if R is None else R(0)
        y_offset = min(self.min_pow_y(), 0)
        y_degree = max(self.max_pow_y(), 0)
        coeffs = []
        n = y_degree - y_offset + 1
        for a in self._f.list():
            k = a.valuation()
            if k is Infinity:
                k = 0
            k -= y_offset
            z = a.list()
            coeffs.append([zero] * k + z + [zero] * (n - len(z) - k))
        while len(coeffs) < self.parent().degree():
            coeffs.append([zero] * n)
        V = FreeModule(self.base_ring() if R is None else R, self.parent().degree())
        coeffs = transpose_list(coeffs)
        return [V(a) for a in coeffs], y_offset


class SpecialHyperellipticQuotientRing(UniqueRepresentation, Parent):
    """
    The special hyperelliptic quotient ring.
    """

    _p = None

    def __init__(self, Q, R=None, invert_y=True):
        r"""
        Initialize ``self``.

        TESTS::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: from sage.schemes.hyperelliptic_curves_smooth_model.monsky_washnitzer import SpecialHyperellipticQuotientRing
            sage: HQR = SpecialHyperellipticQuotientRing(E)
            sage: TestSuite(HQR).run()                                                  # needs sage.rings.real_interval_field

        Check that caching works::

            sage: HQR is SpecialHyperellipticQuotientRing(E)
            True
        """
        if R is None:
            R = Q.base_ring()

        x = PolynomialRing(R, "xx").gen()
        if isinstance(Q, EllipticCurve_generic):
            E = Q
            if E.a1() != 0 or E.a2() != 0:
                raise NotImplementedError("curve must be in Weierstrass normal form")
            Q = -E.change_ring(R).defining_polynomial()(x, 0, 1)
            self._curve = E

        elif isinstance(Q, HyperellipticCurveSmoothModel_generic):
            C = Q
            if C.hyperelliptic_polynomials()[1] != 0:
                raise NotImplementedError("curve must be of form y^2 = Q(x)")
            Q = C.hyperelliptic_polynomials()[0].change_ring(R)
            self._curve = C

        if isinstance(Q, Polynomial):
            self._Q = Q.change_ring(R)
            self._coeffs = self._Q.coefficients(sparse=False)
            if self._coeffs.pop() != 1:
                raise NotImplementedError("polynomial must be monic")
            if not hasattr(self, "_curve"):
                if self._Q.degree() == 3:
                    ainvs = [0, self._Q[2], 0, self._Q[1], self._Q[0]]
                    self._curve = EllipticCurve(ainvs, check_squarefree=R.is_field())
                else:
                    self._curve = HyperellipticCurveSmoothModel(
                        self._Q, check_squarefree=R.is_field()
                    )

        else:
            raise NotImplementedError(
                "must be an elliptic curve or polynomial "
                "Q for y^2 = Q(x)\n(Got element of %s)" % Q.parent()
            )

        self._n = int(Q.degree())
        self._series_ring = (LaurentSeriesRing if invert_y else PolynomialRing)(R, "y")
        self._series_ring_y = self._series_ring.gen(0)
        self._series_ring_0 = self._series_ring.zero()

        Parent.__init__(self, base=R, category=Algebras(R).Commutative())

        self._poly_ring = PolynomialRing(self._series_ring, "x")

        self._x = self.element_class(self, self._poly_ring.gen(0))
        self._y = self.element_class(self, self._series_ring.gen(0))

        self._Q_coeffs = Q.change_ring(self._series_ring).list()
        self._dQ = Q.derivative().change_ring(self)(self._x)
        self._monsky_washnitzer = MonskyWashnitzerDifferentialRing(self)

        self._monomial_diffs = {}
        self._monomial_diff_coeffs = {}

    def _repr_(self) -> str:
        r"""
        String representation.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: x.parent()  # indirect doctest
            SpecialHyperellipticQuotientRing K[x,y,y^-1] / (y^2 = x^5 - 3*x + 1) over Rational Field
        """
        y_inverse = (
            ",y^-1"
            if isinstance(self._series_ring, (LaurentSeriesRing, LazyLaurentSeriesRing))
            else ""
        )
        return "SpecialHyperellipticQuotientRing K[x,y%s] / (y^2 = %s) over %s" % (
            y_inverse,
            self._Q,
            self.base_ring(),
        )

    def base_extend(self, R):
        r"""
        Return the base extension of ``self`` to the ring ``R`` if possible.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: x.parent().base_extend(UniversalCyclotomicField())                    # needs sage.libs.gap
            SpecialHyperellipticQuotientRing K[x,y,y^-1] / (y^2 = x^5 - 3*x + 1)
            over Universal Cyclotomic Field
            sage: x.parent().base_extend(ZZ)
            Traceback (most recent call last):
            ...
            TypeError: no such base extension
        """
        if R.has_coerce_map_from(self.base_ring()):
            return self.change_ring(R)
        raise TypeError("no such base extension")

    def change_ring(self, R):
        r"""
        Return the analog of ``self`` over the ring ``R``.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: x.parent().change_ring(ZZ)
            SpecialHyperellipticQuotientRing K[x,y,y^-1] / (y^2 = x^5 - 3*x + 1) over Integer Ring
        """
        return SpecialHyperellipticQuotientRing(
            self._Q.change_ring(R),
            R,
            isinstance(self._series_ring, (LaurentSeriesRing, LazyLaurentSeriesRing)),
        )

    def _element_constructor_(self, val, offset=0, check=True):
        r"""
        Construct an element of ``self``.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: x.parent()(x^6)
            -(1-y^2)*x + 3*x^2
        """
        if (
            isinstance(val, SpecialHyperellipticQuotientElement)
            and val.parent() is self
        ):
            if offset == 0:
                return val
            else:
                return val << offset
        elif isinstance(val, MonskyWashnitzerDifferential):
            return self._monsky_washnitzer(val)
        return self.element_class(self, val, offset, check)

    @cached_method
    def one(self):
        """
        Return the unit of ``self``.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: x.parent().one()
            1
        """
        return self.element_class(self, self._poly_ring.one(), check=False)

    @cached_method
    def zero(self):
        """
        Return the zero of ``self``.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: x.parent().zero()
            0
        """
        return self.element_class(self, self._poly_ring.zero(), check=False)

    def gens(self):
        """
        Return the generators of ``self``.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: x.parent().gens()
            (x, y*1)
        """
        return (self._x, self._y)

    def x(self):
        r"""
        Return the generator `x` of ``self``.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: x.parent().x()
            x
        """
        return self._x

    def y(self):
        r"""
        Return the generator `y` of ``self``.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: x.parent().y()
            y*1
        """
        return self._y

    def monomial(self, i, j, b=None):
        """
        Return `b y^j x^i`, computed quickly.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: x.parent().monomial(4,5)
            y^5*x^4
        """
        i = int(i)
        j = int(j)

        if 0 < i and i < self._n:
            if b is None:
                by_to_j = self._series_ring_y << (j - 1)
            else:
                by_to_j = self._series_ring(b) << j
            v = [self._series_ring_0] * self._n
            v[i] = by_to_j
            return self.element_class(self, v)

        if b is not None:
            b = self.base_ring()(b)
        return (self._x**i) << j if b is None else b * (self._x**i) << j

    def monomial_diff_coeffs(self, i, j):
        r"""
        Compute coefficients of the basis representation of `d(x^iy^j)`.

        The key here is that the formula for `d(x^iy^j)` is messy
        in terms of `i`, but varies nicely with `j`.

        .. MATH::

             d(x^iy^j) = y^{j-1} (2ix^{i-1}y^2 + j (A_i(x) + B_i(x)y^2)) \frac{dx}{2y},

        where `A,B` have degree at most `n-1` for each
        `i`. Pre-compute `A_i, B_i` for each `i`
        the "hard" way, and the rest are easy.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: x.parent().monomial_diff_coeffs(2,3)
            ((0, -15, 36, 0, 0), (0, 19, 0, 0, 0))
        """
        try:
            return self._monomial_diff_coeffs[i, j]
        except KeyError:
            pass
        if i < self._n:
            try:
                A, B, two_i_x_to_i = self._precomputed_diff_coeffs[i]
            except AttributeError:
                self._precomputed_diff_coeffs = self._precompute_monomial_diffs()
                A, B, two_i_x_to_i = self._precomputed_diff_coeffs[i]
            if i == 0:
                return j * A, j * B
            else:
                return j * A, j * B + two_i_x_to_i
        else:
            dg = self.monomial(i, j).diff()
            coeffs = [dg.extract_pow_y(j - 1), dg.extract_pow_y(j + 1)]
            self._monomial_diff_coeffs[i, j] = coeffs
            return coeffs

    def monomial_diff_coeffs_matrices(self):
        r"""
        Compute tables of coefficients of the basis representation
        of `d(x^iy^j)` for small `i`, `j`.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: x.parent().monomial_diff_coeffs_matrices()
            (
            [0 5 0 0 0]  [0 2 0 0 0]
            [0 0 5 0 0]  [0 0 4 0 0]
            [0 0 0 5 0]  [0 0 0 6 0]
            [0 0 0 0 5]  [0 0 0 0 8]
            [0 0 0 0 0], [0 0 0 0 0]
            )
        """
        self.monomial_diff_coeffs(0, 0)  # precompute stuff
        R = self.base_ring()
        mat_1 = matrix(R, self._n, self._n)
        mat_2 = matrix(R, self._n, self._n)
        for i in range(self._n):
            mat_1[i] = self._precomputed_diff_coeffs[i][1]
            mat_2[i] = self._precomputed_diff_coeffs[i][2]
        return mat_1.transpose(), mat_2.transpose()

    def _precompute_monomial_diffs(self):
        r"""
        Precompute coefficients of the basis representation of `d(x^iy^j)`
        for small `i`, `j`.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: x.parent()._precompute_monomial_diffs()
            [((-3, 0, 0, 0, 5), (0, 0, 0, 0, 0), (0, 0, 0, 0, 0)),
            ((-5, 12, 0, 0, 0), (5, 0, 0, 0, 0), (2, 0, 0, 0, 0)),
            ((0, -5, 12, 0, 0), (0, 5, 0, 0, 0), (0, 4, 0, 0, 0)),
            ((0, 0, -5, 12, 0), (0, 0, 5, 0, 0), (0, 0, 6, 0, 0)),
            ((0, 0, 0, -5, 12), (0, 0, 0, 5, 0), (0, 0, 0, 8, 0))]
        """
        x, y = self.gens()
        R = self.base_ring()
        V = FreeModule(R, self.degree())
        As = []
        for i in range(self.degree()):
            dg = self.monomial(i, 1).diff()
            two_i_x_to_i = R(2 * i) * x ** (i - 1) * y * y if i > 0 else self(0)
            A = dg - self._monsky_washnitzer(two_i_x_to_i)
            As.append(
                (
                    V(A.extract_pow_y(0)),
                    V(A.extract_pow_y(2)),
                    V(two_i_x_to_i.extract_pow_y(2)),
                )
            )
        return As

    def Q(self):
        """
        Return the defining polynomial of the underlying hyperelliptic curve.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5-2*x+1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: x.parent().Q()
            x^5 - 2*x + 1
        """
        return self._Q

    def curve(self):
        """
        Return the underlying hyperelliptic curve.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: x.parent().curve()
            Hyperelliptic Curve over Rational Field defined by y^2 = x^5 - 3*x + 1
        """
        return self._curve

    def degree(self):
        """
        Return the degree of the underlying hyperelliptic curve.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: x.parent().degree()
            5
        """
        return Integer(self._n)

    def prime(self):
        """
        Return the stored prime number `p`.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: x.parent().prime() is None
            True
        """
        return self._p

    def monsky_washnitzer(self):
        """
        Return the stored Monsky-Washnitzer differential ring.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: type(x.parent().monsky_washnitzer())
            <class 'sage.schemes.hyperelliptic_curves_smooth_model.monsky_washnitzer.MonskyWashnitzerDifferentialRing_with_category'>
        """
        return self._monsky_washnitzer

    def is_field(self, proof=True) -> bool:
        """
        Return ``False`` as ``self`` is not a field.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: x,y = E.monsky_washnitzer_gens()
            sage: x.parent().is_field()
            False
        """
        return False

    Element = SpecialHyperellipticQuotientElement


SpecialHyperellipticQuotientRing_class = SpecialHyperellipticQuotientRing


class MonskyWashnitzerDifferential(ModuleElement):
    r"""
    An element of the Monsky-Washnitzer ring of differentials, of
    the form `F dx/2y`.

    EXAMPLES::

        sage: R.<x> = QQ['x']
        sage: C = HyperellipticCurveSmoothModel(x^5 - 4*x + 4)
        sage: x,y = C.monsky_washnitzer_gens()
        sage: MW = C.invariant_differential().parent()
        sage: MW(x)
        x dx/2y
        sage: MW(y)
        y*1 dx/2y
        sage: MW(x, 10)
        y^10*x dx/2y
    """

    def __init__(self, parent, val, offset=0):
        r"""
        Initialize ``self``.

        INPUT:

        - ``parent`` -- Monsky-Washnitzer differential ring (instance of class
          :class:`~MonskyWashnitzerDifferentialRing`

        - ``val`` -- element of the base ring, or list of coefficients

        - ``offset`` -- (default: 0) if nonzero, shift val by `y^\text{offset}`

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: C = HyperellipticCurveSmoothModel(x^5 - 4*x + 4)
            sage: x,y = C.monsky_washnitzer_gens()
            sage: MW = C.invariant_differential().parent()
            sage: elt = MW(x)
            sage: TestSuite(elt).run()
        """
        ModuleElement.__init__(self, parent)
        R = parent.base_ring()
        self._coeff = R._element_constructor_(val, offset)

    def _add_(self, other):
        r"""
        Return the sum of ``self`` and ``other``, both elements of the
        Monsky-Washnitzer ring of differentials.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: C = HyperellipticCurveSmoothModel(x^5 - 4*x + 4)
            sage: x,y = C.monsky_washnitzer_gens()
            sage: w = C.invariant_differential()
            sage: w + w
            2*1 dx/2y
            sage: x*w + w
            (1 + x) dx/2y
            sage: x*w + y*w
            (y*1 + x) dx/2y
        """
        P = self.parent()
        return P.element_class(P, self._coeff + other._coeff)

    def _sub_(self, other):
        r"""
        Return the difference of ``self`` and ``other``, both elements of the
        Monsky-Washnitzer ring of differentials.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: C = HyperellipticCurveSmoothModel(x^5 - 4*x + 4)
            sage: x,y = C.monsky_washnitzer_gens()
            sage: w = C.invariant_differential()
            sage: w-w
            0 dx/2y
            sage: x*w-w
            (-1 + x) dx/2y
            sage: w - x*w - y*w
            ((1-y)*1 - x) dx/2y
        """
        P = self.parent()
        return P.element_class(P, self._coeff - other._coeff)

    def __neg__(self):
        r"""
        Return the additive inverse of ``self``.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: C = HyperellipticCurveSmoothModel(x^5 - 4*x + 4)
            sage: x,y = C.monsky_washnitzer_gens()
            sage: w = C.invariant_differential()
            sage: -w
            -1 dx/2y
            sage: -((y-x)*w)
            (-y*1 + x) dx/2y
        """
        P = self.parent()
        return P.element_class(P, -self._coeff)

    def _lmul_(self, a):
        r"""
        Return `self * a`.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: C = HyperellipticCurveSmoothModel(x^5 - 4*x + 4)
            sage: x,y = C.monsky_washnitzer_gens()
            sage: w = C.invariant_differential()
            sage: w*x
            x dx/2y
            sage: (w*x)*2
            2*x dx/2y
            sage: w*y
            y*1 dx/2y
            sage: w*(x+y)
            (y*1 + x) dx/2y
        """
        P = self.parent()
        return P.element_class(P, self._coeff * a)

    def _rmul_(self, a):
        r"""
        Return `a * self`.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: C = HyperellipticCurveSmoothModel(x^5 - 4*x + 4)
            sage: x,y = C.monsky_washnitzer_gens()
            sage: w = C.invariant_differential()
            sage: x*w
            x dx/2y
            sage: 2*(x*w)
            2*x dx/2y
            sage: y*w
            y*1 dx/2y
            sage: (x+y)*w
            (y*1 + x) dx/2y
        """
        P = self.parent()
        return P.element_class(P, a * self._coeff)

    def coeff(self):
        r"""
        Return `A`, where this element is `A dx/2y`.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: C = HyperellipticCurveSmoothModel(x^5 - 4*x + 4)
            sage: x,y = C.monsky_washnitzer_gens()
            sage: w = C.invariant_differential()
            sage: w
            1 dx/2y
            sage: w.coeff()
            1
            sage: (x*y*w).coeff()
            y*x
        """
        return self._coeff

    def __bool__(self):
        r"""
        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: C = HyperellipticCurveSmoothModel(x^5 - 4*x + 4)
            sage: x,y = C.monsky_washnitzer_gens()
            sage: w = C.invariant_differential()
            sage: not w
            False
            sage: not 0*w
            True
            sage: not x*y*w
            False
        """
        return bool(self._coeff)

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: C = HyperellipticCurveSmoothModel(x^5 - 4*x + 4)
            sage: x,y = C.monsky_washnitzer_gens()
            sage: w = C.invariant_differential()
            sage: w
            1 dx/2y
            sage: (2*x+y)*w
            (y*1 + 2*x) dx/2y
        """
        s = self._coeff._repr_()
        if s.find("+") != -1 or s.find("-") > 0:
            s = "(%s)" % s
        return s + " dx/2y"

    def _latex_(self):
        r"""
        Return the latex representation of ``self``.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: C = HyperellipticCurveSmoothModel(x^5 - 4*x + 4)
            sage: x,y = C.monsky_washnitzer_gens()
            sage: w = C.invariant_differential()
            sage: latex(w)
            1 \frac{dx}{2y}
            sage: latex(x*w)
            x \frac{dx}{2y}
        """
        s = self._coeff._latex_()
        if s.find("+") != -1 or s.find("-") > 0:
            s = "\\left(%s\\right)" % s
        return s + " \\frac{dx}{2y}"

    def _richcmp_(self, other, op):
        """
        Rich comparison of ``self`` to ``other``.

        .. TODO::

            This does not compare elements by any reduction;
            it only compares the coefficients. The comparison
            should be done against a normal form or possibly
            after some reduction steps.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: C = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: x,y = C.monsky_washnitzer_gens()
            sage: (y^-1).diff() == (y^-1).diff()
            True

        This element is the zero differential in the ring, but it does
        not compare as equal because the representative is different::

            sage: MW = C.invariant_differential().parent()
            sage: (y^-1).diff().reduce_neg_y()
            ((y^-1)*1, 0 dx/2y)
            sage: (y^-1).diff() == MW.zero()
            False
        """
        return richcmp(self._coeff, other._coeff, op)

    def extract_pow_y(self, k):
        r"""
        Return the power of `y` in `A` where ``self`` is `A dx/2y`.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: C = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: x,y = C.monsky_washnitzer_gens()
            sage: A = y^5 - x*y^3
            sage: A.extract_pow_y(5)
            [1, 0, 0, 0, 0]
            sage: (A * C.invariant_differential()).extract_pow_y(5)
            [1, 0, 0, 0, 0]
        """
        return self._coeff.extract_pow_y(k)

    def min_pow_y(self):
        r"""
        Return the minimum power of `y` in `A` where ``self`` is `A dx/2y`.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: C = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: x,y = C.monsky_washnitzer_gens()
            sage: w = y^5 * C.invariant_differential()
            sage: w.min_pow_y()
            5
            sage: w = (x^2*y^4 + y^5) * C.invariant_differential()
            sage: w.min_pow_y()
            4
        """
        return self._coeff.min_pow_y()

    def max_pow_y(self):
        r"""
        Return the maximum power of `y` in `A` where ``self`` is `A dx/2y`.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: C = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: x,y = C.monsky_washnitzer_gens()
            sage: w = y^5 * C.invariant_differential()
            sage: w.max_pow_y()
            5
            sage: w = (x^2*y^4 + y^5) * C.invariant_differential()
            sage: w.max_pow_y()
            5
        """
        return self._coeff.max_pow_y()

    def reduce_neg_y(self):
        r"""
        Use homology relations to eliminate negative powers of `y`.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: C = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: x,y = C.monsky_washnitzer_gens()
            sage: (y^-1).diff().reduce_neg_y()
            ((y^-1)*1, 0 dx/2y)
            sage: (y^-5*x^2+y^-1*x).diff().reduce_neg_y()
            ((y^-1)*x + (y^-5)*x^2, 0 dx/2y)
        """
        S = self.parent().base_ring()
        R = S.base_ring()
        M = self.parent().helper_matrix()
        p = S._p
        f = S.zero()
        reduced = self
        for j in range(self.min_pow_y() + 1, 0):
            if p is not None and p.divides(j):
                cs = [a / j for a in reduced.extract_pow_y(j - 1)]
            else:
                j_inverse = ~R(j)
                cs = [a * j_inverse for a in reduced.extract_pow_y(j - 1)]
            lin_comb = M * vector(M.base_ring(), cs)
            if lin_comb.is_zero():
                continue
            g = S.sum(S.monomial(i, j, val) for i, val in enumerate(lin_comb) if val)
            if not g.is_zero():
                f += g
                reduced -= g.diff()

        return f, reduced

    def reduce_neg_y_fast(self, even_degree_only=False):
        r"""
        Use homology relations to eliminate negative powers of `y`.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: x, y = E.monsky_washnitzer_gens()
            sage: (y^-1).diff().reduce_neg_y_fast()
            ((y^-1)*1, 0 dx/2y)
            sage: (y^-5*x^2+y^-1*x).diff().reduce_neg_y_fast()
            ((y^-1)*x + (y^-5)*x^2, 0 dx/2y)

        It leaves nonnegative powers of `y` alone::

            sage: y.diff()
            (-3*1 + 5*x^4) dx/2y
            sage: y.diff().reduce_neg_y_fast()
            (0, (-3*1 + 5*x^4) dx/2y)
        """
        #        prof = Profiler()
        #        prof("reduce setup")
        S = self.parent().base_ring()
        R = S.base_ring()
        M = self.parent().helper_matrix()

        #        prof("extract coeffs")
        coeffs, offset = self.coeffs(R)
        V = coeffs[0].parent()

        if offset == 0:
            return S(0), self

        #        prof("loop %s"%self.min_pow_y())
        forms = []
        p = S._p
        for j in range(self.min_pow_y() + 1, 0):
            if (even_degree_only and j % 2 == 0) or coeffs[j - offset - 1].is_zero():
                forms.append(V(0))
            else:
                # this is a total hack to deal with the fact that we're using
                # rational numbers to approximate fixed precision p-adics
                if p is not None and j % 3 == 1:
                    try:
                        v = coeffs[j - offset - 1]
                        for kk in range(len(v)):
                            a = v[kk]
                            ppow = p ** max(-a.valuation(S._p), 0)
                            v[kk] = ((a * ppow) % S._prec_cap) / ppow
                    except AttributeError:
                        pass
                lin_comb = ~R(j) * (M * coeffs[j - offset - 1])
                forms.append(lin_comb)
                for i in lin_comb.nonzero_positions():
                    # g = lin_comb[i] x^i y^j
                    # self -= dg
                    coeffs[j - offset + 1] -= (
                        lin_comb[i] * S.monomial_diff_coeffs(i, j)[1]
                    )

        #        prof("recreate forms")
        f = S(forms, offset + 1)
        reduced = S._monsky_washnitzer(coeffs[-1 - offset :], -1)
        return f, reduced

    def reduce_neg_y_faster(self, even_degree_only=False):
        r"""
        Use homology relations to eliminate negative powers of `y`.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: C = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: x,y = C.monsky_washnitzer_gens()
            sage: (y^-1).diff().reduce_neg_y()
            ((y^-1)*1, 0 dx/2y)
            sage: (y^-5*x^2+y^-1*x).diff().reduce_neg_y_faster()
            ((y^-1)*x + (y^-5)*x^2, 0 dx/2y)
        """
        # Timings indicate that this is not any faster after all...

        S = self.parent().base_ring()
        R = S.base_ring()
        M = self.parent().helper_matrix()

        coeffs, offset = self.coeffs(R)
        V = coeffs[0].parent()
        zeroV = V(0)

        if offset == 0:
            return S.zero(), self

        # See monomial_diff_coeffs
        # this is the B_i and x_to_i contributions respectively for all i
        d_mat_1, d_mat_2 = S.monomial_diff_coeffs_matrices()

        forms = []
        for j in range(self.min_pow_y() + 1, 0):
            if coeffs[j - offset - 1].is_zero():
                forms.append(zeroV)
            else:
                # this is a total hack to deal with the fact that we're using
                # rational numbers to approximate fixed precision p-adics
                if j % 3 == 0:
                    try:
                        v = coeffs[j - offset - 1]
                        for kk in range(len(v)):
                            a = v[kk]
                            ppow = S._p ** max(-a.valuation(S._p), 0)
                            v[kk] = ((a * ppow) % S._prec_cap) / ppow
                    except AttributeError:
                        pass
                j_inverse = ~R(j)
                lin_comb = M * coeffs[j - offset - 1]
                forms.append(j_inverse * lin_comb)
                coeffs[j - offset + 1] -= (d_mat_1 + j_inverse * d_mat_2) * lin_comb

        f = S(forms, offset + 1)
        reduced = S._monsky_washnitzer(coeffs[-1 - offset :], -1)
        #        reduced = self - f.diff()
        return f, reduced

    def reduce_pos_y(self):
        r"""
        Use homology relations to eliminate positive powers of `y`.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: C = HyperellipticCurveSmoothModel(x^3-4*x+4)
            sage: x,y = C.monsky_washnitzer_gens()
            sage: (y^2).diff().reduce_pos_y()
            (y^2*1, 0 dx/2y)
            sage: (y^2*x).diff().reduce_pos_y()
            (y^2*x, 0 dx/2y)
            sage: (y^92*x).diff().reduce_pos_y()
            (y^92*x, 0 dx/2y)
            sage: w = (y^3 + x).diff()
            sage: w += w.parent()(x)
            sage: w.reduce_pos_y_fast()
            (y^3*1 + x, x dx/2y)
        """
        S = self.parent().base_ring()
        n = S.Q().degree()
        f = S.zero()
        reduced = self
        for j in range(self.max_pow_y(), 0, -1):
            for i in range(n - 1, -1, -1):
                c = reduced.extract_pow_y(j)[i]
                if c:
                    g = S.monomial(0, j + 1) if i == n - 1 else S.monomial(i + 1, j - 1)
                    dg = g.diff()
                    denom = dg.extract_pow_y(j)[i]
                    c /= denom
                    c = g.parent()(c)
                    f += c * g
                    reduced -= c * dg

        return f, reduced

    def reduce_pos_y_fast(self, even_degree_only=False):
        r"""
        Use homology relations to eliminate positive powers of `y`.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^3 - 4*x + 4)
            sage: x, y = E.monsky_washnitzer_gens()
            sage: y.diff().reduce_pos_y_fast()
            (y*1, 0 dx/2y)
            sage: (y^2).diff().reduce_pos_y_fast()
            (y^2*1, 0 dx/2y)
            sage: (y^2*x).diff().reduce_pos_y_fast()
            (y^2*x, 0 dx/2y)
            sage: (y^92*x).diff().reduce_pos_y_fast()
            (y^92*x, 0 dx/2y)
            sage: w = (y^3 + x).diff()
            sage: w += w.parent()(x)
            sage: w.reduce_pos_y_fast()
            (y^3*1 + x, x dx/2y)
        """
        S = self.parent().base_ring()
        R = S.base_ring()
        n = S.Q().degree()

        coeffs, offset = self.coeffs(R)
        V = coeffs[0].parent()
        zeroV = V(0)
        forms = [V(0), V(0)]

        for j in range(self.max_pow_y(), -1, -1):
            if (even_degree_only and j % 2) or (j > 0 and coeffs[j - offset].is_zero()):
                forms.append(zeroV)
                continue

            form = V(0)
            i = n - 1
            c = coeffs[j - offset][i]
            if c:
                dg_coeffs = S.monomial_diff_coeffs(0, j + 1)[0]
                c /= dg_coeffs[i]
                forms[len(forms) - 2][0] = c
                # self -= c d(y^{j+1})
                coeffs[j - offset] -= c * dg_coeffs

            if j == 0:
                # the others are basis elements
                break

            for i in range(n - 2, -1, -1):
                c = coeffs[j - offset][i]
                if c:
                    dg_coeffs = S.monomial_diff_coeffs(i + 1, j - 1)
                    denom = dg_coeffs[1][i]
                    c /= denom
                    form[i + 1] = c
                    # self -= c d(x^{i+1} y^{j-1})
                    coeffs[j - offset] -= c * dg_coeffs[1]
                    coeffs[j - offset - 2] -= c * dg_coeffs[0]
            forms.append(form)

        forms.reverse()
        f = S(forms)
        reduced = self.parent()(coeffs[: 1 - offset], offset)
        return f, reduced

    def reduce(self):
        r"""
        Use homology relations to find `a` and `f` such that this element is
        equal to `a + df`, where `a` is given in terms of the `x^i dx/2y`.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: C = HyperellipticCurveSmoothModel(x^5 - 4*x + 4)
            sage: x,y = C.monsky_washnitzer_gens()
            sage: w = (y*x).diff()
            sage: w.reduce()
            (y*x, 0 dx/2y)

            sage: w = x^4 * C.invariant_differential()
            sage: w.reduce()
            (1/5*y*1, 4/5*1 dx/2y)

            sage: w = sum(QQ.random_element() * x^i * y^j
            ....:         for i in [0..4] for j in [-3..3]) * C.invariant_differential()
            sage: f, a = w.reduce()
            sage: f.diff() + a - w
            0 dx/2y
        """
        n = self.parent().base_ring().Q().degree()
        f1, a = self.reduce_neg_y()
        f2, a = a.reduce_pos_y()
        f = f1 + f2

        c = a.extract_pow_y(0)[n - 1]
        if c:
            x, y = self.parent().base_ring().gens()
            g = y
            dg = g.diff()
            c = g.parent()(c / dg.extract_pow_y(0)[n - 1])
            f += c * g
            a -= c * dg

        return f, a

    def reduce_fast(self, even_degree_only=False):
        r"""
        Use homology relations to find `a` and `f` such that this element is
        equal to `a + df`, where `a` is given in terms of the `x^i dx/2y`.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^3 - 4*x + 4)
            sage: x, y = E.monsky_washnitzer_gens()
            sage: x.diff().reduce_fast()
            (x, (0, 0))
            sage: y.diff().reduce_fast()
            (y*1, (0, 0))
            sage: (y^-1).diff().reduce_fast()
            ((y^-1)*1, (0, 0))
            sage: (y^-11).diff().reduce_fast()
            ((y^-11)*1, (0, 0))
            sage: (x*y^2).diff().reduce_fast()
            (y^2*x, (0, 0))
        """
        f1, reduced = self.reduce_neg_y_fast(even_degree_only)
        f2, reduced = reduced.reduce_pos_y_fast(even_degree_only)
        #        f1, reduced = self.reduce_neg_y()
        #        f2, reduced = reduced.reduce_pos_y()
        v = reduced.extract_pow_y(0)
        v.pop()
        V = FreeModule(self.base_ring().base_ring(), len(v))
        return f1 + f2, V(v)

    def coeffs(self, R=None):
        """
        Used to obtain the raw coefficients of a differential, see
        :meth:`SpecialHyperellipticQuotientElement.coeffs`

        INPUT:

        - ``R`` -- an (optional) base ring in which to cast the coefficients

        OUTPUT: the raw coefficients of `A` where ``self`` is `A dx/2y`

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: C = HyperellipticCurveSmoothModel(x^5 - 4*x + 4)
            sage: x,y = C.monsky_washnitzer_gens()
            sage: w = C.invariant_differential()
            sage: w.coeffs()
            ([(1, 0, 0, 0, 0)], 0)
            sage: (x*w).coeffs()
            ([(0, 1, 0, 0, 0)], 0)
            sage: (y*w).coeffs()
            ([(0, 0, 0, 0, 0), (1, 0, 0, 0, 0)], 0)
            sage: (y^-2*w).coeffs()
            ([(1, 0, 0, 0, 0), (0, 0, 0, 0, 0), (0, 0, 0, 0, 0)], -2)
        """
        return self._coeff.coeffs(R)

    def coleman_integral(self, P, Q):
        r"""
        Compute the definite integral of ``self`` from `P` to `Q`.

        INPUT:

        - `P`, `Q` -- two points on the underlying curve

        OUTPUT: `\int_P^Q \text{self}`

        EXAMPLES::

            sage: K = pAdicField(5,7)
            sage: E = EllipticCurve(K,[-31/3,-2501/108]) #11a
            sage: P = E(K(14/3), K(11/2))
            sage: w = E.invariant_differential()
            sage: w.coleman_integral(P, 2*P)
            O(5^6)

            sage: Q = E([3,58332])
            sage: w.coleman_integral(P,Q)
            2*5 + 4*5^2 + 3*5^3 + 4*5^4 + 3*5^5 + O(5^6)
            sage: w.coleman_integral(2*P,Q)
            2*5 + 4*5^2 + 3*5^3 + 4*5^4 + 3*5^5 + O(5^6)
            sage: (2*w).coleman_integral(P, Q) == 2*(w.coleman_integral(P, Q))
            True
        """
        return self.parent().base_ring().curve().coleman_integral(self, P, Q)

    integrate = coleman_integral


class MonskyWashnitzerDifferentialRing(UniqueRepresentation, Module):
    r"""
    A ring of Monsky--Washnitzer differentials over ``base_ring``.
    """

    def __init__(self, base_ring):
        r"""
        Initialize ``self``.

        TESTS::

            sage: R.<x> = QQ['x']
            sage: E = HyperellipticCurveSmoothModel(x^5 - 3*x + 1)
            sage: from sage.schemes.hyperelliptic_curves_smooth_model.monsky_washnitzer import SpecialHyperellipticQuotientRing, MonskyWashnitzerDifferentialRing
            sage: S = SpecialHyperellipticQuotientRing(E)
            sage: DR = MonskyWashnitzerDifferentialRing(S)
            sage: TestSuite(DR).run()                                                   # needs sage.rings.real_interval_field

        Check that caching works::

            sage: DR is MonskyWashnitzerDifferentialRing(S)
            True
        """
        Module.__init__(self, base_ring)

    def invariant_differential(self):
        r"""
        Return `dx/2y` as an element of ``self``.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: C = HyperellipticCurveSmoothModel(x^5 - 4*x + 4)
            sage: MW = C.invariant_differential().parent()
            sage: MW.invariant_differential()
            1 dx/2y
        """
        return self.element_class(self, 1)

    def base_extend(self, R):
        """
        Return a new differential ring which is ``self`` base-extended to `R`.

        INPUT:

        - ``R`` -- ring

        OUTPUT:

        Self, base-extended to `R`.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: C = HyperellipticCurveSmoothModel(x^5 - 4*x + 4)
            sage: MW = C.invariant_differential().parent()
            sage: MW.base_ring()
            SpecialHyperellipticQuotientRing K[x,y,y^-1] / (y^2 = x^5 - 4*x + 4)
             over Rational Field
            sage: MW.base_extend(Qp(5,5)).base_ring()                                   # needs sage.rings.padics
            SpecialHyperellipticQuotientRing K[x,y,y^-1] / (y^2 = (1 + O(5^5))*x^5
                        + (1 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + O(5^5))*x + 4 + O(5^5))
             over 5-adic Field with capped relative precision 5
        """
        return MonskyWashnitzerDifferentialRing(self.base_ring().base_extend(R))

    def change_ring(self, R):
        """
        Return a new differential ring which is ``self`` with the coefficient
        ring changed to `R`.

        INPUT:

        - ``R`` -- ring of coefficients

        OUTPUT: ``self`` with the coefficient ring changed to `R`

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: C = HyperellipticCurveSmoothModel(x^5 - 4*x + 4)
            sage: MW = C.invariant_differential().parent()
            sage: MW.base_ring()
            SpecialHyperellipticQuotientRing K[x,y,y^-1] / (y^2 = x^5 - 4*x + 4)
             over Rational Field
            sage: MW.change_ring(Qp(5,5)).base_ring()                                   # needs sage.rings.padics
            SpecialHyperellipticQuotientRing K[x,y,y^-1] / (y^2 = (1 + O(5^5))*x^5
                        + (1 + 4*5 + 4*5^2 + 4*5^3 + 4*5^4 + O(5^5))*x + 4 + O(5^5))
             over 5-adic Field with capped relative precision 5
        """
        return MonskyWashnitzerDifferentialRing(self.base_ring().change_ring(R))

    def degree(self):
        """
        Return the degree of `Q(x)`, where the model of the underlying
        hyperelliptic curve of ``self`` is given by `y^2 = Q(x)`.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: C = HyperellipticCurveSmoothModel(x^5 - 4*x + 4)
            sage: MW = C.invariant_differential().parent()
            sage: MW.Q()
            x^5 - 4*x + 4
            sage: MW.degree()
            5
        """
        return self.base_ring().degree()

    def dimension(self):
        """
        Return the dimension of ``self``.

        EXAMPLES::

            sage: # needs sage.rings.padics
            sage: R.<x> = QQ['x']
            sage: C = HyperellipticCurveSmoothModel(x^5 - 4*x + 4)
            sage: K = Qp(7,5)
            sage: CK = C.change_ring(K)
            sage: MW = CK.invariant_differential().parent()
            sage: MW.dimension()
            4
        """
        return self.base_ring().degree() - 1

    def Q(self):
        """
        Return `Q(x)` where the model of the underlying hyperelliptic curve
        of ``self`` is given by `y^2 = Q(x)`.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: C = HyperellipticCurveSmoothModel(x^5 - 4*x + 4)
            sage: MW = C.invariant_differential().parent()
            sage: MW.Q()
            x^5 - 4*x + 4
        """
        return self.base_ring().Q()

    @cached_method
    def x_to_p(self, p):
        """
        Return and cache `x^p`, reduced via the relations coming from the
        defining polynomial of the hyperelliptic curve.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: C = HyperellipticCurveSmoothModel(x^5 - 4*x + 4)
            sage: MW = C.invariant_differential().parent()
            sage: MW.x_to_p(3)
            x^3
            sage: MW.x_to_p(5)
            -(4-y^2)*1 + 4*x
            sage: MW.x_to_p(101) is MW.x_to_p(101)
            True
        """
        return self.base_ring().x() ** p

    @cached_method
    def frob_Q(self, p):
        r"""
        Return and cache `Q(x^p)`, which is used in computing the image of
        `y` under a `p`-power lift of Frobenius to `A^{\dagger}`.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: C = HyperellipticCurveSmoothModel(x^5 - 4*x + 4)
            sage: MW = C.invariant_differential().parent()
            sage: MW.frob_Q(3)
            -(60-48*y^2+12*y^4-y^6)*1 + (192-96*y^2+12*y^4)*x - (192-48*y^2)*x^2 + 60*x^3
            sage: MW.Q()(MW.x_to_p(3))                                                  # needs sage.rings.real_interval_field
            -(60-48*y^2+12*y^4-y^6)*1 + (192-96*y^2+12*y^4)*x - (192-48*y^2)*x^2 + 60*x^3
            sage: MW.frob_Q(11) is MW.frob_Q(11)
            True
        """
        return self.base_ring()._Q.change_ring(self.base_ring())(self.x_to_p(p))

    def frob_invariant_differential(self, prec, p):
        r"""
        Kedlaya's algorithm allows us to calculate the action of Frobenius on
        the Monsky-Washnitzer cohomology. First we lift `\phi` to `A^{\dagger}`
        by setting

        .. MATH::

            \phi(x) = x^p,
            \qquad\qquad
            \phi(y) = y^p \sqrt{1  + \frac{Q(x^p) - Q(x)^p}{Q(x)^p}}.

        Pulling back the differential `dx/2y`, we get

        .. MATH::

           \phi^*(dx/2y) = px^{p-1} y(\phi(y))^{-1} dx/2y
           = px^{p-1} y^{1-p} \sqrt{1+ \frac{Q(x^p) - Q(x)^p}{Q(x)^p}} dx/2y.

        Use Newton's method to calculate the square root.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: C = HyperellipticCurveSmoothModel(x^5 - 4*x + 4)
            sage: prec = 2
            sage: p = 7
            sage: MW = C.invariant_differential().parent()
            sage: MW.frob_invariant_differential(prec, p)
            ((67894400*y^-20-81198880*y^-18+40140800*y^-16-10035200*y^-14+1254400*y^-12-62720*y^-10)*1
             - (119503944*y^-20-116064242*y^-18+43753472*y^-16-7426048*y^-14+514304*y^-12-12544*y^-10+1568*y^-8-70*y^-6-7*y^-4)*x
             + (78905288*y^-20-61014016*y^-18+16859136*y^-16-2207744*y^-14+250880*y^-12-37632*y^-10+3136*y^-8-70*y^-6)*x^2
             - (39452448*y^-20-26148752*y^-18+8085490*y^-16-2007040*y^-14+376320*y^-12-37632*y^-10+1568*y^-8)*x^3
             + (21102144*y^-20-18120592*y^-18+8028160*y^-16-2007040*y^-14+250880*y^-12-12544*y^-10)*x^4) dx/2y
        """
        prof = Profiler()
        prof("setup")
        # TODO, would it be useful to be able to take Frobenius of any element? Less efficient?
        x, y = self.base_ring().gens()
        prof("x_to_p")
        x_to_p_less_1 = x ** (p - 1)
        x_to_p = x * x_to_p_less_1

        # cache for future use
        self.x_to_p.set_cache(p, x_to_p)

        prof("frob_Q")
        a = self.frob_Q(p) >> 2 * p  # frobQ * y^{-2p}

        prof("sqrt")

        #        Q = self.base_ring()._Q
        #        three_halves = Q.parent().base_ring()(Rational((3,2)))
        #        one_half = Q.parent().base_ring()(Rational((1,2)))
        three_halves = self.base_ring()._series_ring.base_ring()(Rational((3, 2)))
        one_half = self.base_ring()._series_ring.base_ring()(Rational((1, 2)))
        half_a = a._rmul_(one_half)

        # We are solving for t = a^{-1/2} = (F_pQ y^{-p})^{-1/2}
        # Newton's method converges because we know the root is in the same residue class as 1.

        #        t = self.base_ring()(1)
        t = self.base_ring()(three_halves) - half_a
        # first iteration trivial, start with prec 2

        for cur_prec in newton_method_sizes(prec)[2:]:
            # newton_method_sizes = [1, 2, ...]
            y_prec = -(2 * cur_prec - 1) * p + 1
            # binomial expansion is  $\sum p^{k+1} y^{-(2k+1)p+1} f(x)$
            # so if we are only correct mod p^prec,
            # can ignore y powers less than y_prec
            t_cube = (t * t * t).truncate_neg(y_prec)
            t = t._rmul_(three_halves) - (half_a * t_cube).truncate_neg(y_prec)
            # t = (3/2) t - (1/2) a t^3

        prof("compose")
        F_dx_y = (p * x_to_p_less_1 * t) >> (p - 1)  # px^{p-1} sqrt(a) * y^{-p+1}

        prof("done")
        return MonskyWashnitzerDifferential(self, F_dx_y)

    def frob_basis_elements(self, prec, p):
        r"""
        Return the action of a `p`-power lift of Frobenius on the basis.

        .. MATH::

            \{ dx/2y, x dx/2y, ..., x^{d-2} dx/2y \},

        where `d` is the degree of the underlying hyperelliptic curve.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: C = HyperellipticCurveSmoothModel(x^5 - 4*x + 4)
            sage: prec = 1
            sage: p = 5
            sage: MW = C.invariant_differential().parent()
            sage: MW.frob_basis_elements(prec, p)
            [((92000*y^-14-74200*y^-12+32000*y^-10-8000*y^-8+1000*y^-6-50*y^-4)*1
              - (194400*y^-14-153600*y^-12+57600*y^-10-9600*y^-8+600*y^-6)*x
              + (204800*y^-14-153600*y^-12+38400*y^-10-3200*y^-8)*x^2
              - (153600*y^-14-76800*y^-12+9600*y^-10)*x^3
              + (63950*y^-14-18550*y^-12+1600*y^-10-400*y^-8+50*y^-6+5*y^-4)*x^4) dx/2y,
             (-(1391200*y^-14-941400*y^-12+302000*y^-10-76800*y^-8+14400*y^-6-1320*y^-4+30*y^-2)*1
              + (2168800*y^-14-1402400*y^-12+537600*y^-10-134400*y^-8+16800*y^-6-720*y^-4)*x
              - (1596800*y^-14-1433600*y^-12+537600*y^-10-89600*y^-8+5600*y^-6)*x^2
              + (1433600*y^-14-1075200*y^-12+268800*y^-10-22400*y^-8)*x^3
              - (870200*y^-14-445350*y^-12+63350*y^-10-3200*y^-8+600*y^-6-30*y^-4-5*y^-2)*x^4) dx/2y,
             ((19488000*y^-14-15763200*y^-12+4944400*y^-10-913800*y^-8+156800*y^-6-22560*y^-4+1480*y^-2-10)*1
              - (28163200*y^-14-18669600*y^-12+5774400*y^-10-1433600*y^-8+268800*y^-6-25440*y^-4+760*y^-2)*x
              + (15062400*y^-14-12940800*y^-12+5734400*y^-10-1433600*y^-8+179200*y^-6-8480*y^-4)*x^2
              - (12121600*y^-14-11468800*y^-12+4300800*y^-10-716800*y^-8+44800*y^-6)*x^3
              + (9215200*y^-14-6952400*y^-12+1773950*y^-10-165750*y^-8+5600*y^-6-720*y^-4+10*y^-2+5)*x^4) dx/2y,
             (-(225395200*y^-14-230640000*y^-12+91733600*y^-10-18347400*y^-8+2293600*y^-6-280960*y^-4+31520*y^-2-1480-10*y^2)*1
              + (338048000*y^-14-277132800*y^-12+89928000*y^-10-17816000*y^-8+3225600*y^-6-472320*y^-4+34560*y^-2-720)*x
              - (172902400*y^-14-141504000*y^-12+58976000*y^-10-17203200*y^-8+3225600*y^-6-314880*y^-4+11520*y^-2)*x^2
              + (108736000*y^-14-109760000*y^-12+51609600*y^-10-12902400*y^-8+1612800*y^-6-78720*y^-4)*x^3
              - (85347200*y^-14-82900000*y^-12+31251400*y^-10-5304150*y^-8+367350*y^-6-8480*y^-4+760*y^-2+10-5*y^2)*x^4) dx/2y]
        """
        F_i = self.frob_invariant_differential(prec, p)
        x_to_p = self.x_to_p(p)
        F = [F_i]
        for i in range(1, self.degree() - 1):
            F_i *= x_to_p
            F.append(F_i)
        return F

    def helper_matrix(self):
        r"""
        We use this to solve for the linear combination of
        `x^i y^j` needed to clear all terms with `y^{j-1}`.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: C = HyperellipticCurveSmoothModel(x^5 - 4*x + 4)
            sage: MW = C.invariant_differential().parent()
            sage: MW.helper_matrix()
            [ 256/2101  320/2101  400/2101  500/2101  625/2101]
            [-625/8404  -64/2101  -80/2101 -100/2101 -125/2101]
            [-125/2101 -625/8404  -64/2101  -80/2101 -100/2101]
            [-100/2101 -125/2101 -625/8404  -64/2101  -80/2101]
            [ -80/2101 -100/2101 -125/2101 -625/8404  -64/2101]
        """
        try:
            return self._helper_matrix
        except AttributeError:
            pass

        # The smallest y term of (1/j) d(x^i y^j) is constant for all j.
        x, y = self.base_ring().gens()
        n = self.degree()
        L = [(y * x**i).diff().extract_pow_y(0) for i in range(n)]
        A = matrix(L).transpose()
        if A.base_ring() not in IntegralDomains():
            # must be using integer_mod or something to approximate
            self._helper_matrix = (~A.change_ring(QQ)).change_ring(A.base_ring())
        else:
            self._helper_matrix = ~A
        return self._helper_matrix

    def _element_constructor_(self, val=0, offset=0):
        r"""
        Construct an element of ``self``.

        INPUT:

        - ``parent`` -- Monsky-Washnitzer differential ring (instance of class
          :class:`~MonskyWashnitzerDifferentialRing`
        - ``val`` -- element of the base ring, or list of coefficients
        - ``offset`` -- (default: 0) if nonzero, shift val by `y^\text{offset}`

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: C = HyperellipticCurveSmoothModel(x^5 - 4*x + 4)
            sage: MW = C.invariant_differential().parent()
            sage: MW(3)
            3*1 dx/2y
        """
        if isinstance(val, MonskyWashnitzerDifferential):
            val = val._coeff
        return self.element_class(self, val, offset)

    Element = MonskyWashnitzerDifferential


MonskyWashnitzerDifferentialRing_class = MonskyWashnitzerDifferentialRing
