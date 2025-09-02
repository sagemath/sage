"""
Arithmetic on the Jacobian

This module implements the group operation in the Picard group of a
hyperelliptic curve, represented as divisors in Mumford
representation, using Cantor's algorithm.
"""

from sage.groups.generic import order_from_multiple
from sage.misc.cachefunc import cached_method
from sage.rings.finite_rings.finite_field_base import FiniteField as FiniteField_generic
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.schemes.generic.morphism import SchemeMorphism
from sage.structure.element import AdditiveGroupElement
from sage.structure.richcmp import richcmp


class MumfordDivisorClassField(AdditiveGroupElement, SchemeMorphism):
    r"""
    An element of a Jacobian defined over a field, i.e. in
    `J(K) = \mathrm{Pic}^0_K(C)`.
    """

    def __init__(self, parent, u, v, check=True):
        """
        Create an element of the Jacobian of a hyperelliptic curve.

        EXAMPLES::

            sage: R.<x> = GF(19)[]
            sage: H = HyperellipticCurveSmoothModel(x^5 + x, x^2 + 1)
            sage: J = H.jacobian()
            sage: D = J(x^2 + 14*x + 16, 3*x + 4); D # indirect doctest
            (x^2 + 14*x + 16, 3*x + 4)
        """
        SchemeMorphism.__init__(self, parent)

        if not isinstance(u, Polynomial) or not isinstance(v, Polynomial):
            raise TypeError(f"arguments u={u} and v={v} must be polynomials")

        # TODO:
        # 1. allow elements of the base field as input
        #   (in particular something like (u,v) = (x-alpha, 0))

        # Ensure the divisor is valid
        if check:
            f, h = parent.curve().hyperelliptic_polynomials()
            assert (
                v**2 + v * h - f
            ) % u == 0, f"u={u}, v={v} do not define a divisor on the Jacobian"

            # TODO: should we automatically do reduction here if the degree of u is
            # too large?

        self._parent = parent
        self._u = u
        self._v = v

    def parent(self):
        """
        Return the parent of the divisor class.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurve(x^5 - 1)
            sage: J = H.jacobian()
            sage: D = J(x-1,0)
            sage: D.parent()
            Set of rational points of Jacobian of Hyperelliptic Curve over Rational Field defined by y^2 = x^5 - 1
        """
        return self._parent

    def scheme(self):
        r"""
        Return the scheme this morphism maps to; or, where this divisor lives.

        .. WARNING::

           Although a pointset is defined over a specific field, the
           scheme returned may be over a different (usually smaller)
           field.  The example below demonstrates this: the pointset
           is determined over a number field of absolute degree 2 but
           the scheme returned is defined over the rationals.

        EXAMPLES::

            sage: x = QQ['x'].gen()
            sage: f = x^5 + x
            sage: H = HyperellipticCurve(f)
            sage: F.<a> = NumberField(x^2 - 2, 'a')                                     # needs sage.rings.number_field
            sage: J = H.jacobian()(F); J                                                # needs sage.rings.number_field
            Set of rational points of Jacobian of Hyperelliptic Curve
             over Number Field in a with defining polynomial x^2 - 2
             defined by y^2 = x^5 + x

        ::

            sage: P = J(H.lift_x(F(1)))                                                 # needs sage.rings.number_field
            sage: P.scheme()                                                            # needs sage.rings.number_field
            Jacobian of Hyperelliptic Curve over Rational Field defined by y^2 = x^5 + x
        """
        return self.codomain()

    def _repr_(self) -> str:
        """
        Return the Mumford presentation of the divisor class.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: H = HyperellipticCurve(x^5 + 2*x^3 + x, x^3 + 1)
            sage: P = H([0,-1])
            sage: Q = H([0,0])
            sage: J = H.jacobian()
            sage: D = J(P,Q); D # indirect doctest
            (x^2, y + x + 1)
        """
        return f"({self._u}, {self._v})"

    def is_zero(self):
        """
        Return ``True`` if this point is zero.

        EXAMPLES::

            sage: x = polygen(GF(5))
            sage: H = HyperellipticCurveSmoothModel(x^5 + 3 * x + 1)
            sage: J = H.jacobian(); J
            Jacobian of Hyperelliptic Curve over Finite Field of size 5 defined by y^2 = x^5 + 3*x + 1
            sage: points = list(J)
            sage: [P for P in points if P.is_zero()]
            [(1, 0)]
            sage: [P for P in points if P == 0]
            [(1, 0)]
        """
        return self._u.is_one() and self._v.is_zero()

    def uv(self):
        """
        Return the `u` and `v` component of this Mumford divisor.

        EXAMPLES::

            sage: x = polygen(GF(1993))
            sage: H = HyperellipticCurveSmoothModel(x^7 + 3 * x + 1)
            sage: J = H.jacobian(); J
            Jacobian of Hyperelliptic Curve over Finite Field of size 1993 defined by y^2 = x^7 + 3*x + 1
            sage: u, v = x^3 + 1570*x^2 + 1930*x + 81, 368*x^2 + 1478*x + 256
            sage: P = J(u, v); P
            (x^3 + 1570*x^2 + 1930*x + 81, 368*x^2 + 1478*x + 256)
            sage: P.uv() == (u, v)
            True
        """
        return (self._u, self._v)

    def _richcmp_(self, other, op) -> bool:
        """
        Method for rich comparison.

        TESTS::

            sage: x = polygen(GF(23))
            sage: H = HyperellipticCurveSmoothModel(x^7 + x + 1)
            sage: J = H.jacobian(); J
            Jacobian of Hyperelliptic Curve over Finite Field of size 23 defined by y^2 = x^7 + x + 1
            sage: P = J.random_element()
            sage: P == P
            True
            sage: P != P
            False

            sage: Z = J(0); Z
            (1, 0)
            sage: Z == 0
            True
            sage: Z != 0
            False
        """
        # _richcmp_ is called after type unification/coercion
        assert isinstance(other, MumfordDivisorClassField)
        return richcmp(tuple(self), tuple(other), op)

    def __iter__(self):
        """
        TESTS:

        Indirect tests::

            sage: R.<x> = GF(13)[]
            sage: H = HyperellipticCurveSmoothModel(x^5 - x + 1)
            sage: J = H.jacobian()
            sage: P = J(x + 5, 12)
            sage: list(P)
            [x + 5, 12]
            sage: tuple(P)
            (x + 5, 12)
        """
        yield from [self._u, self._v]

    def __getitem__(self, n):
        """
        Return the n-th item in the Mumford presentation of the divisor.

        TESTS::

            sage: R.<x> = GF(23)[]
            sage: H = HyperellipticCurveSmoothModel(x^5 + x + 1, x^2 + x +  1)
            sage: J = H.jacobian()
            sage: D = J(x^2 + 21*x + 10, 14*x + 22)
            sage: D[0] # indirect doctest
            x^2 + 21*x + 10
            sage: D[1] # indirect doctest
            14*x + 22
        """
        return (self._u, self._v)[n]

    def __hash__(self):
        """
        Compute the hash value of this element.

        TESTS::

            sage: R.<x> = GF(23)[]
            sage: H = HyperellipticCurveSmoothModel(x^5 + x + 1, x^2 + x +  1)
            sage: J = H.jacobian()
            sage: D = J(x^2 + 21*x + 10, 14*x + 22)
            sage: hash(D) == hash(2*D)
            False
        """
        return hash(tuple(self))

    def __bool__(self):
        """
        Return "True" if this is not the zero element of the Jacobian.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurveSmoothModel(x^2 + x, x^3 + 1)
            sage: J = H.jacobian()
            sage: bool(J(x+1,0))
            True
            sage: bool(J.zero())
            False
        """
        return not self.is_zero()

    @cached_method
    def order(self):
        """
        Returns the order of self.
        This is only implemented over finite fields.

        EXAMPLES::

            sage: K = FiniteField(7)
            sage: R.<x> = K[]
            sage: H = HyperellipticCurveSmoothModel(x^6 + 3*x + 2)
            sage: JK = Jacobian(H)(K)
            sage: D = JK(x^2 + 5*x + 6, 6*x + 3)
            sage: D.order()
            38
        """
        if not isinstance(self.base_ring(), FiniteField_generic):
            raise NotImplementedError(
                "this is only implemented for Jacobians over a finite field"
            )
        n = self.parent().order()
        return order_from_multiple(self, n)

    def degree(self):
        """
        Returns the degree of the affine part of the divisor.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurveSmoothModel(6*x^5 + 9*x^4 - x^3 - 3*x^2, 1)
            sage: J = H.jacobian()
            sage: J.zero().degree()
            0
            sage: J(x, 0).degree()
            1
            sage: J(x^2 + 1/2*x, -x - 1).degree()
            2
        """
        return self._u.degree()

    def _add_(self, other):
        """
        Add `other` to the divisor `self`.

        TESTS::

            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurveSmoothModel(x^5 - x^4 + x^2 - x, 1)
            sage: J = H.jacobian()
            sage: D1 = J(x^2 + x, x)
            sage: D2 = J(x^2, x - 1)
            sage: D1 + D2
            (x^2 + x, -1)

            sage: Jp = J.change_ring(GF(13))
            sage: D1, D2, D3 = [Jp.random_element() for i in range(3)]
            sage: D1 + D2 + D3 == D3 + D1 + D2
            True
        """
        # `other` should be coerced automatically before reaching here
        assert isinstance(other, type(self))

        # Collect data from HyperellipticCurve
        H = self.parent().curve()
        g = H.genus()

        # Extract out mumford coordinates
        u1, v1 = self.uv()
        u2, v2 = other.uv()

        # Step one: cantor composition of the two divisors
        u3, v3 = self._parent.cantor_composition(u1, v1, u2, v2)

        # Step two: cantor reduction of the above to ensure
        # the degree of u is smaller than g + 1
        while u3.degree() > g:
            u3, v3 = self._parent.cantor_reduction(u3, v3)
        u3 = u3.monic()
        v3 = v3 % u3

        return self._parent(u3, v3, check=False)

    def _neg_(self):
        """
        Multiply the divisor by `[-1]`.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurveSmoothModel(x^5 + 2*x^4 + 3*x^3 + 2*x^2 + x, 1)
            sage: P = H(1,0,0)
            sage: Q = H(0,0)
            sage: J = H.jacobian()
            sage: D = J(P,Q); D
            (x, -1)
            sage: -D #indirect doctest
            (x, 0)
            sage: J(Q,P) == - D #indirect doctest
            True
        """
        _, h = self.parent().curve().hyperelliptic_polynomials()
        u0, v0 = self.uv()

        if h.is_zero():
            return self._parent(u0, -v0, check=False)

        v1 = (-h - v0) % u0
        return self._parent(u0, v1, check=False)


# =======================================================================
# Child classes for representation for ramified, inert and split models
# =======================================================================


class MumfordDivisorClassFieldRamified(MumfordDivisorClassField):
    def __init__(self, parent, u, v, check=True):
        """
            Create an element of the Jacobian of a ramified
            hyperelliptic curve.

            TESTS::

                sage: R.<x> = GF(5)[]
                sage: H = HyperellipticCurveSmoothModel(x^7 + 3*x + 2, x^2 + 1)
                sage: H.is_ramified()
                True
                sage: J = Jacobian(H)
                sage: D = J(x^3 + x^2 + 3*x + 1, 4*x^2 + 2*x + 3)
                sage: type(D)
                <class 'sage.schemes.hyperelliptic_curves_smooth_model.jacobian_morphism.MumfordDivisorClassFieldRamified'>
        """
        if not parent.curve().is_ramified():
            raise TypeError("hyperelliptic curve must be ramified")

        super().__init__(parent, u, v, check=check)


class MumfordDivisorClassFieldInert(MumfordDivisorClassField):
    def __init__(self, parent, u, v, check=True):
        """
            Create an element of the Jacobian of a ramified
            hyperelliptic curve.

            TESTS::

                sage: R.<x> = GF(5)[]
                sage: H = HyperellipticCurveSmoothModel(2*x^6 + 1)
                sage: H.is_inert()
                True
                sage: J = Jacobian(H)
                sage: D = J(x^2 + x + 4, 4*x)
                sage: type(D)
                <class 'sage.schemes.hyperelliptic_curves_smooth_model.jacobian_morphism.MumfordDivisorClassFieldInert'>
        """
        if not parent.curve().is_inert():
            raise TypeError("hyperelliptic curve must be inert")

        if u.degree() % 2:
            raise ValueError(f"mumford coordinate {u} must have even degree")

        g = parent.curve().genus()
        self._n = (g - u.degree()) // 2
        super().__init__(parent, u, v, check=check)

    def __repr__(self):
        """
            Return a representation of the element.

            TESTS::

                sage: R.<x> = GF(5)[]
                sage: H = HyperellipticCurveSmoothModel(2*x^6 + 1)
                sage: H.is_inert()
                True
                sage: J = Jacobian(H)
                sage: D = J(x^2 + x + 4, 4*x); D # indirect doctest
                (x^2 + x + 4, 4*x : 0)
        """
        return f"({self._u}, {self._v} : {self._n})"


class MumfordDivisorClassFieldSplit(MumfordDivisorClassField):
    def __init__(self, parent, u, v, n=0, check=True):
        """
            Create an element of the Jacobian of a split
            hyperelliptic curve.

            TESTS::

                sage: R.<x> = GF(5)[]
                sage: H = HyperellipticCurveSmoothModel(x^6 + 1)
                sage: H.is_split()
                True
                sage: J = Jacobian(H)
                sage: D = J(x^2 + 3, 3)
                sage: type(D)
                <class 'sage.schemes.hyperelliptic_curves_smooth_model.jacobian_morphism.MumfordDivisorClassFieldSplit'>
        """
        if not parent.curve().is_split():
            raise TypeError("hyperelliptic curve must be split")

        # Ensure the weight is set correctly
        g = parent.curve().genus()
        assert 0 <= n <= (g - u.degree())
        self._n = n
        self._m = g - u.degree() - n

        super().__init__(parent, u, v, check=check)

    def __repr__(self):
        """
            Return a representation of the element.

            TESTS::

                sage: R.<x> = GF(5)[]
                sage: H = HyperellipticCurveSmoothModel(x^6 + 1)
                sage: J = Jacobian(H)
                sage: D = J(x^2 + 3, 3); D
                (x^2 + 3, 3 : 0)
        """
        return f"({self._u}, {self._v} : {self._n})"

    def is_zero(self):
        """
            Return ``True`` if this element is zero.

            EXAMPLES:

            We consider a genus-3 hyperelliptic curve with two points
            at infinity ``P_0``, ``P_1``. Elements of the Jacobian are
            represented as ``[D - 2P_0 - P_1]``.
            In particular the zero element has Mumford presentation
            `(1, 0 : 2)`::

                sage: R.<x> = QQ[]
                sage: H = HyperellipticCurveSmoothModel(x^8  + 1)
                sage: J = Jacobian(H)
                sage: D = J(1,0,2); D
                (1, 0 : 2)
                sage: D.is_zero()
                True

            There are more elements that are only supported at infinity,
            but are non-zero::

                sage: [P0,P1] = H.points_at_infinity()
                sage: D = J(P0,P1); D
                (1, 0 : 3)
                sage: D.is_zero()
                False
        """
        g = self._parent.curve().genus()
        return self._u.is_one() and self._v.is_zero() and self._n == (g + 1) // 2

    def __iter__(self):
        """
        TESTS:

        Indirect tests::

            sage: R.<x> = GF(7)[]
            sage: H = HyperellipticCurveSmoothModel(x^6 - x + 1)
            sage: J = H.jacobian()
            sage: P = J(x + 2, 5)
            sage: list(P)
            [x + 2, 5, 1]
            sage: tuple(P)
            (x + 2, 5, 1)
            sage: Q = J(x + 2, 5, 0)
            sage: list(Q)
            [x + 2, 5, 0]
            sage: tuple(Q)
            (x + 2, 5, 0)
            sage: P == Q
            False
        """
        yield from [self._u, self._v, self._n]

    def __hash__(self):
        """
        Compute the hash value of this element.

        TESTS::

            sage: R.<x> = GF(23)[]
            sage: H = HyperellipticCurveSmoothModel(x^6 + 1)
            sage: J = H.jacobian()
            sage: D1 = J(x + 6, 17, 0)
            sage: D2 = J(x + 6, 17, 1)
            sage: hash(D1) == hash(D2)
            False
        """
        return hash(tuple(self))

    def _add_(self, other):
        r"""
        Return the sum of the two elements.

        Follows algorithm 3.7 of

        Efficient Arithmetic on Hyperelliptic Curves With Real Representation
        David J. Mireles Morales (2008)
        https://www.math.auckland.ac.nz/~sgal018/Dave-Mireles-Full.pdf

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurveSmoothModel(x^5 - x^4 + x^2 - x, x^3 + 1)
            sage: J = Jacobian(H)
            sage: D1 = J(x^2 + x, 0)
            sage: D2 = J(x^2, -x)
            sage: D1 + D2
            (x^2 - 3/2*x + 1/2, -5/2*x + 1/2 : 0)
            sage: D3 = J(x - 1, -2)
            sage: D1 + D3
            (x^2 - 1/2*x, 5/4*x - 1 : 0)
        """
        # Ensure we are adding two divisors
        assert isinstance(other, type(self))

        # Collect data from HyperellipticCurve
        H = self.parent().curve()
        g = H.genus()

        # Extract out mumford coordinates
        u1, v1 = self.uv()
        u2, v2 = other.uv()

        # Extract out integers for weights
        n1, n2 = self._n, other._n

        # Step one: cantor composition of the two divisors
        u3, v3, n3 = self._parent.cantor_composition(u1, v1, n1, u2, v2, n2)

        # Step two: cantor reduction of the above to ensure
        # the degree of u is smaller than g + 1
        while u3.degree() > (g + 1):
            u3, v3, n3 = self._parent.cantor_reduction(u3, v3, n3)
        u3 = u3.monic()

        # Step three: compose and then reduce at infinity to ensure
        # unique representation of D
        while n3 < 0 or n3 > g - u3.degree():
            if n3 < 0:
                u3, v3, n3 = self._parent.cantor_compose_at_infinity(
                    u3, v3, n3, plus=False
                )
            else:
                u3, v3, n3 = self._parent.cantor_compose_at_infinity(
                    u3, v3, n3, plus=True
                )

        return self._parent(u3, v3, n3, check=False)

    def _neg_(self):
        r"""
        Return the divisor multiplied by -1.

        Follows algorithm 3.8 of

        Efficient Arithmetic on Hyperelliptic Curves With Real Representation
        David J. Mireles Morales (2008)
        https://www.math.auckland.ac.nz/~sgal018/Dave-Mireles-Full.pdf

        TESTS::

            sage: R.<x> = GF(101)[]
            sage: H = HyperellipticCurveSmoothModel(x^5 + 1, x^3 + x)
            sage: assert H.is_split()
            sage: J = Jacobian(H)
            sage: [P0,P1] = H.points_at_infinity()
            sage: P = H.random_point()
            sage: Q = H.hyperelliptic_involution(P)
            sage: - J(P,P0) == J(Q,P1)
            True
        """
        # Collect data from HyperellipticCurve
        H = self.parent().curve()
        _, h = H.hyperelliptic_polynomials()
        g = H.genus()

        u0, v0 = self.uv()
        n0, m0 = self._n, self._m

        # Case for even genus
        if not (g % 2):
            v1 = (-h - v0) % u0
            u1 = u0
            n1 = m0
        # Odd genus, positive n0
        elif n0 > 0:
            v1 = (-h - v0) % u0
            # Note: I think this is a typo in the paper
            # n1 = g - m0 - u1.degree() + 1
            u1 = u0
            n1 = m0 + 1
        else:
            # Composition at infinity always with plus=True.
            # want to "substract" \infty_+ - \infty_-
            (u1, v1, n1) = self._parent.cantor_compose_at_infinity(
                u0, -h - v0, n0, plus=True
            )
            n1 = n1 - n0 + m0 + 1

        return self._parent(u1, v1, n1, check=False)
