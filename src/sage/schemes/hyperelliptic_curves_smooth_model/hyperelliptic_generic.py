from sage.functions.all import log
from sage.misc.cachefunc import cached_method
from sage.rings.big_oh import O
from sage.rings.laurent_series_ring import LaurentSeriesRing
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.real_mpfr import RR

# TODO: move this
from sage.schemes.hyperelliptic_curves_smooth_model.weighted_projective_curve import (
    WeightedProjectiveCurve,
)
from sage.schemes.hyperelliptic_curves_smooth_model.weighted_projective_space import (
    WeightedProjectiveSpace,
)


class HyperellipticCurveSmoothModel_generic(WeightedProjectiveCurve):
    def __init__(self, defining_polynomial, f, h, genus):
        self._genus = genus
        self._hyperelliptic_polynomials = (f, h)

        self._polynomial_ring = f.parent()
        self._base_ring = f.base_ring()

        # TODO: is this simply genus + 1
        self._d = max(h.degree(), (f.degree() + 1) // 2)

        # Initalise the underlying curve
        A = WeightedProjectiveSpace([1, self._genus + 1, 1], self._base_ring)
        WeightedProjectiveCurve.__init__(self, A, defining_polynomial)

    def _repr_(self):
        old_gen = str(self._polynomial_ring.gen())
        f, h = self._hyperelliptic_polynomials

        # TODO:
        # The old class has these weird internal gens and then
        # printing polynomial rings to change output.
        #
        # Will do something hacky here and we can talk about it.
        f_str, h_str = repr(f).replace(old_gen, "x"), repr(h).replace(old_gen, "x")

        if h:
            if h.is_one():
                curve = f"y^2 + y = {f_str}"
            else:
                curve = f"y^2 + ({h_str})*y = {f_str}"
        else:
            curve = f"y^2 = {f_str}"

        return f"Hyperelliptic Curve over {self.base_ring()} defined by {curve}"

    def genus(self):
        """
        Return the genus of the hyperelliptic curve.

        EXAMPLES:
            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurveSmoothModel(x^7 + 3*x + 2)
            sage: H.genus()
            3
            sage: H = HyperellipticCurveSmoothModel(x^3 + 2, x^5 + 1)
            sage: H.genus()
            4
        """
        return self._genus

    def base_ring(self):
        """
        Return the base ring of the hyperelliptic curve.

         EXAMPLES::

            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurveSmoothModel(x^7 + 3*x + 2)
            sage: H.base_ring()
            Rational Field

            sage: S.<x> = FiniteField(19)[]
            sage: H = HyperellipticCurveSmoothModel(x^5 - x + 2)
            sage: H.base_ring()
            Finite Field of size 19
        """
        return self._base_ring

    def change_ring(self, R):
        """
        Return this hyperelliptic curve over a new ring "R".

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurveSmoothModel(x^5 - 10*x + 9)
            sage: K = Qp(3, 5)                                                          # optional - sage.rings.padics
            sage: L.<a> = K.extension(x^30 - 3)                                         # optional - sage.rings.padics
            sage: HK = H.change_ring(K)                                                 # optional - sage.rings.padics
            sage: HL = HK.change_ring(L); HL                                            # optional - sage.rings.padics
            Hyperelliptic Curve over 3-adic Eisenstein Extension Field in a defined by x^30 - 3 defined by y^2 = x^5 + (2 + 2*a^30 + a^60 + 2*a^90 + 2*a^120 + O(a^150))*x + a^60 + O(a^210)

            sage: R.<x> = FiniteField(7)[]                                              # optional - sage.rings.finite_rings
            sage: H = HyperellipticCurveSmoothModel(x^8 + x + 5)                                   # optional - sage.rings.finite_rings
            sage: H.base_extend(FiniteField(7^2, 'a'))                                  # optional - sage.rings.finite_rings
            Hyperelliptic Curve over Finite Field in a of size 7^2 defined by y^2 = x^8 + x + 5

        It is also possible to compute the reduction at a prime by changing
        the base ring to the residue field::

            sage: R.<x> = PolynomialRing(QQ);
            sage: H = HyperellipticCurve(R([0, -1, 2, 0, -2]), R([0, 1, 0, 1])); #LMFDB label: 763.a.763.1
            sage: H.change_ring(FiniteField(2))
            Hyperelliptic Curve over Finite Field of size 2 defined by y^2 + (x^3 + x)*y = x
            sage: H.change_ring(FiniteField(3))
            Hyperelliptic Curve over Finite Field of size 3 defined by y^2 + (x^3 + x)*y = x^4 + 2*x^2 + 2*x

        Note that this only works when the curve has good reduction at p::

            sage: H.change_ring(FiniteField(7))
            Traceback (most recent call last):
            ...
            ValueError: not a hyperelliptic curve: singularity in the provided affine patch
        """
        from sage.schemes.hyperelliptic_curves_smooth_model.hyperelliptic_constructor import (
            HyperellipticCurveSmoothModel,
        )

        f, h = self._hyperelliptic_polynomials
        fR = f.change_ring(R)
        hR = h.change_ring(R)
        return HyperellipticCurveSmoothModel(fR, hR)

    base_extend = change_ring

    def polynomial_ring(self):
        """
        Returns the parent of f,h, where self is the hyperelliptic curve
        defined by y^2 + h*y = f.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurveSmoothModel(x^7 + 3*x + 2)
            sage: H.polynomial_ring()
            Univariate Polynomial Ring in x over Rational Field

            sage: # optional - sage.rings.padics
            sage: K = Qp(7, 10)
            sage: HK = H.change_ring(K)
            sage: HK.polynomial_ring()
            Univariate Polynomial Ring in x over 7-adic Field with capped relative precision 10
        """
        return self._polynomial_ring

    def hyperelliptic_polynomials(self):
        """
        Return the polynomials (f, h) such that
        C : y^2 + h*y = f

        EXAMPLES:

            sage: R.<x> = PolynomialRing(QQ)
            sage: H = HyperellipticCurveSmoothModel(R([0, 1, 2, 0, 0, 1]), R([1, 1, 0, 1]))
            sage: H.hyperelliptic_polynomials()
            (x^5 + 2*x^2 + x, x^3 + x + 1)

            sage: H = HyperellipticCurveSmoothModel(x^7 + x + 2)
            sage: H.hyperelliptic_polynomials()
            (x^7 + x + 2, 0)
        """
        return self._hyperelliptic_polynomials

    def roots_at_infinity(self):
        """
        Compute the roots of: Y^2 + h[d]Y - f[2d] = 0.
        When the curve is ramified, we expect one root, when
        the curve is inert or split we expect zero or two roots.
        """
        if hasattr(self, "_alphas"):
            return self._alphas

        f, h = self._hyperelliptic_polynomials
        x = f.parent().gen()
        d = self._d

        if h.is_zero():
            coeff = f[2 * d]
            # Handle the ramified case
            if coeff.is_zero():
                return [coeff]
            return f[2 * d].sqrt(all=True)

        self._alphas = (x**2 + x * h[d] - f[2 * d]).roots(multiplicities=False)
        return self._alphas

    def is_split(self):
        """
        Return True if the curve is split, i.e. there are two rational
        points at infinity.

        EXAMPLES:

            sage: R.<x> = PolynomialRing(QQ)
            sage: H = HyperellipticCurveSmoothModel(x^6+1, x^3+1)
            sage: H.is_split()
            False

            sage: HK = H.change_ring(FiniteField(19))
            sage: HK.is_split()
            True
        """
        return len(self.roots_at_infinity()) == 2

    def is_ramified(self):
        """
        Return True if the curve is ramified, i.e. there is one rational
        point at infinity.

        EXAMPLES:

            sage: R.<x> = PolynomialRing(QQ)
            sage: H = HyperellipticCurveSmoothModel(x^5+1)
            sage: H.is_ramified()
            True

            sage: H = HyperellipticCurveSmoothModel(x^5+1, x^3+1)
            sage: H.is_ramified()
            False
        """
        return len(self.roots_at_infinity()) == 1

    def is_inert(self):
        """
        Return True if the curve is inert, i.e. there are no rational
        points at infinity.

        EXAMPLES:

            sage: R.<x> = PolynomialRing(QQ)
            sage: H = HyperellipticCurveSmoothModel(x^6+1,-x^3+1)
            sage: H.is_inert()
            True

            sage: K.<a> = QQ.extension(x^2+x-1)
            sage: HK = H.change_ring(K)
            sage: HK.is_inert()
            False

            sage: HF = H.change_ring(FiniteField(29))
            sage: HF.is_inert()
            False
        """
        return len(self.roots_at_infinity()) == 0

    def infinite_polynomials(self):
        """
        TODO: the name of this function could be better?

        Computes G^±(x) for curves in the split degree model used for
        Cantor composition with points at infinity when performing
        arithmetic with the Jac(H).

        .. SEEALSO::

            :func:`~sage.schemes.hyperelliptic_curves_smooth_model.jacobian_homset_split.cantor_compose_at_infinity`
        """
        if hasattr(self, "_infinite_polynomials"):
            return self._infinite_polynomials

        alphas = self.roots_at_infinity()

        # This function only makes sense for the split model
        if not len(alphas) == 2:
            raise ValueError("hyperelliptic curve does not have the split model")

        f, h = self._hyperelliptic_polynomials
        alpha_plus, alpha_minus = alphas
        d = self._d

        # Construct G_plus from alpha_plus
        g = [None] * (d + 1)
        g[d] = alpha_plus
        for i in range(d - 1, -1, -1):
            # We need (g * (g + h))[x^(i + d)] to match f_{i + d}
            the_rest = g[d] * h[i] + sum(
                g[k] * (g[i + d - k] + h[i + d - k]) for k in range(i + 1, d)
            )
            g[i] = (f[i + d] - the_rest) / (2 * g[d] + h[d])

        G_plus = self._polynomial_ring(g)
        G_minus = -G_plus - h

        # Checks for the assumptions on G^±
        genus = self.genus()
        assert G_plus.degree() <= (genus + 1)
        assert (G_plus**2 + h * G_plus - f).degree() <= genus
        assert G_minus.leading_coefficient() == alpha_minus

        self._infinite_polynomials = G_plus, G_minus
        return self._infinite_polynomials

    def points_at_infinity(self):
        """
        Compute the points at infinity on the curve. Assumes we are using
        a weighted projective model for the curve
        """
        return [self.point([1, y, 0], check=False) for y in self.roots_at_infinity()]

    def is_x_coord(self, x):
        """
        Return True if ``x`` is the `x`-coordinate of a point on this curve.

        .. SEEALSO::

            See also :meth:`lift_x` to find the point(s) with a given
            `x`-coordinate.  This function may be useful in cases where
            testing an element of the base field for being a square is
            faster than finding its square root.

        INPUT:

        - ``x`` -- an element of the base ring of the curve

        OUTPUT:

        A bool stating whether or not `x` is a x-coordinate of a point on the curve

        EXAMPLES:

        When `x` is the `x`-coordinate of a rational point on the
        curve, we can request these::

            sage: R.<x> = PolynomialRing(QQ)
            sage: f = x^5 + x^3 + 1
            sage: H = HyperellipticCurveSmoothModel(f)
            sage: H.is_x_coord(0)
            True

        There are no rational points with `x`-coordinate 3::

            sage: H.is_x_coord(3)
            False

        The function also handles the case when `h(x)` is not zero::

            sage: R.<x> = PolynomialRing(QQ)
            sage: f = x^5 + x^3 + 1
            sage: h = x + 1
            sage: H = HyperellipticCurveSmoothModel(f, h)
            sage: H.is_x_coord(1)
            True

        We can perform these operations over finite fields too::

            sage: # needs sage.rings.finite_rings
            sage: R.<x> = PolynomialRing(GF(163))
            sage: f = x^7 + x + 1
            sage: H = HyperellipticCurveSmoothModel(f)
            sage: H.is_x_coord(13)
            True

        Including the case of characteristic two::

            sage: # needs sage.rings.finite_rings
            sage: F.<z4> = GF(2^4)
            sage: R.<x> = PolynomialRing(F)
            sage: f = x^7 + x^3 + 1
            sage: h = x + 1
            sage: H = HyperellipticCurveSmoothModel(f, h)
            sage: H.is_x_coord(z4^3 + z4^2 + z4)
            True

        AUTHORS:

        - Giacomo Pope (2024): adapted from :meth:`lift_x`

        TESTS:

        The `x`-coordinate must be defined over the base field of the curve::

            sage: p = 11
            sage: F = GF(11)
            sage: F_ext = GF(11^2)
            sage: R.<x> = PolynomialRing(F)
            sage: f = x^7 + x^3 + 1
            sage: H = HyperellipticCurveSmoothModel(f)
            sage: H.is_x_coord(F_ext.gen())
            Traceback (most recent call last):
            ...
            TypeError: x must be coercible into the base ring of the curve
        """
        f, h = self.hyperelliptic_polynomials()
        K = self.base_ring()
        try:
            x = K(x)
        except (ValueError, TypeError):
            raise TypeError("x must be coercible into the base ring of the curve")

        # When h is zero then x is a valid coordinate if y2 is square
        if not h:
            y2 = f(x)
            return y2.is_square()

        # Generic case for h != 0
        a = f(x)
        b = h(x)

        # Special case for char 2
        if K.characteristic() == 2:
            R = f.parent()  # Polynomial ring K[x]
            F = R([-a, b, 1])
            return bool(F.roots())

        # Otherwise x is a point on the curve if the discriminant is a square
        D = b * b + 4 * a
        return D.is_square()

    def lift_x(self, x, all=False):
        """
        Return one or all finite points with given `x`-coordinate.

        This method is deterministic: It returns the same data each
        time when called again with the same `x`.

        INPUT:

        - ``x`` -- an element of the base ring of the curve

        - ``all`` (bool, default ``False``) -- if ``True``, return a
          (possibly empty) list of all points; if ``False``, return
          just one point, or raise a :class:`ValueError` if there are none.

        OUTPUT:

        A point or list of up to two points on this curve.

        .. SEEALSO::

            :meth:`is_x_coord`

        AUTHORS:

        - Giacomo Pope (2024): Allowed for the case of characteristic two

        EXAMPLES:

        When `x` is the `x`-coordinate of a rational point on the
        curve, we can request these::

            sage: R.<x> = PolynomialRing(QQ)
            sage: f = x^5 + x^3 + 1
            sage: H = HyperellipticCurveSmoothModel(f)
            sage: H.lift_x(0)
            (0 : -1 : 1)
            sage: H.lift_x(4, all=True)
            [(4 : -33 : 1), (4 : 33 : 1)]

        There are no rational points with `x`-coordinate 3::

            sage: H.lift_x(3)
            Traceback (most recent call last):
            ...
            ValueError: No point with x-coordinate 3 on Hyperelliptic Curve over Rational Field defined by y^2 = x^5 + x^3 + 1

        An empty list is returned when there are no points and ``all=True``::

            sage: H.lift_x(3, all=True)
            []

        The function also handles the case when `h(x)` is not zero::

            sage: R.<x> = PolynomialRing(QQ)
            sage: f = x^5 + x^3 + 1
            sage: h = x + 1
            sage: H = HyperellipticCurveSmoothModel(f, h)
            sage: H.lift_x(1)
            (1 : -3 : 1)

        We can perform these operations over finite fields too::

            sage: # needs sage.rings.finite_rings
            sage: R.<x> = PolynomialRing(GF(163))
            sage: f = x^7 + x + 1
            sage: H = HyperellipticCurveSmoothModel(f)
            sage: H.lift_x(13)
            (13 : 41 : 1)

        Including the case of characteristic two::

            sage: # needs sage.rings.finite_rings
            sage: F.<z4> = GF(2^4)
            sage: R.<x> = PolynomialRing(F)
            sage: f = x^7 + x^3 + 1
            sage: h = x + 1
            sage: H = HyperellipticCurveSmoothModel(f, h)
            sage: H.lift_x(z4^3 + z4^2 + z4, all=True)
            [(z4^3 + z4^2 + z4 : z4^2 + z4 + 1 : 1), (z4^3 + z4^2 + z4 : z4^3 : 1)]

        Points at infinity are not included, as they do not have a unique x-coordinate::

            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurveSmoothModel(x^8 + 1)
            sage: H(1, -1, 0)
            (1 : -1 : 0)
            sage: H.lift_x(1, all=True)
            []

        TESTS::

            sage: # needs sage.rings.finite_rings
            sage: F1 = GF(11)
            sage: F2 = GF(13)
            sage: R.<x> = PolynomialRing(F1)
            sage: f = x^7 + x^3 + 1
            sage: H = HyperellipticCurveSmoothModel(f)
            sage: H.lift_x(F2.random_element())
            Traceback (most recent call last):
            ...
            ValueError: x must have a common parent with the base ring

        Ensure that :issue:`37097` is fixed::

            sage: # needs sage.rings.finite_rings
            sage: F.<z4> = GF(2^4)
            sage: R.<x> = PolynomialRing(F)
            sage: f = x^7 + x^3 + 1
            sage: h = x + 1
            sage: H = HyperellipticCurveSmoothModel(f, h)
            sage: H.lift_x(z4^3 + z4^2 + z4, all=True)
            [(z4^3 + z4^2 + z4 : z4^2 + z4 + 1 : 1), (z4^3 + z4^2 + z4 : z4^3 : 1)]
        """
        from sage.structure.element import get_coercion_model

        cm = get_coercion_model()

        f, h = self.hyperelliptic_polynomials()
        K = self.base_ring()

        # Compute the common parent between the base ring of the curve and
        # the parent of the input x-coordinate.
        try:
            L = cm.common_parent(x.parent(), K)
            x = L(x)
        except (TypeError, ValueError):
            raise ValueError("x must have a common parent with the base ring")

        # First we compute the y-coordinates the given x-coordinate
        ys = []
        one = L.one()

        # When h is zero we find all y-coordinates with a single sqrt
        if not h:
            y2 = f(x)
            # When y2 is not a square, ys will be an empty list
            ys = y2.sqrt(all=True, extend=False)
        # Otherwise we need roots of the discriminant
        else:
            a = f(x)
            b = h(x)
            # Special case for char 2
            if K.characteristic() == 2:
                R = f.parent()
                F = R([-a, b, 1])
                ys = F.roots(L, multiplicities=False)
            else:
                D = b * b + 4 * a
                # When D is not a square, ys will be an empty list
                ys = [(-b + d) / 2 for d in D.sqrt(all=True, extend=False)]

        if ys:
            ys.sort()  # Make lifting deterministic
            if all:
                return [self.point([x, y, one], check=False) for y in ys]
            else:
                return self.point([x, ys[0], one], check=False)

        if all:
            return []
        else:
            raise ValueError(f"No point with x-coordinate {x} on {self}")

    def affine_coordinates(self, P):
        """
        Returns the affine coordinates of a point P of self.
        That is for P = [X,Y,Z], the output is X/Z¸ Y/Z^(g+1).

        EXAMPLES::
            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurveSmoothModel(x^6 - 1)
            sage: P = H.point([2,0,2])
            sage: H.affine_coordinates(P)
            (1, 0)
            sage: Q = H.point([1,1,0])
            sage: H.affine_coordinates(Q)
            Traceback (most recent call last):
            ...
            ValueError: The point (1 : 1 : 0) is not an affine point of Hyperelliptic Curve over Rational Field defined by y^2 = x^6 - 1
        """
        if P[2] == 0:
            raise ValueError(f"The point {P} is not an affine point of {self}")

        g = self.genus()
        return P[0] / P[2], P[1] / P[2] ** (g + 1)

    def is_weierstrass_point(self, P):
        """
        Return True if P is a Weierstrass point of self.
        TODO: It would be better to define this function for points directly.

        EXAMPLES::

        We consider the hyperelliptic curve `y^2 = x^6 -1` over the rational numbers.
        For instance, the point P = (1,0) is a Weierstrass point,
        while the points at infinity are not::

            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurveSmoothModel(x^6 - 1)
            sage: P = H.point([1,0])
            sage: H.is_weierstrass_point(P)
            True
            sage: Q = H.point([1,1,0])
            sage: H.is_weierstrass_point(Q)
            False

        This also works for hyperelliptic curves with `h(x)` nonzero.
        Note that in this case the `y`-coordinate of a Weierstrass point
        is not necessarily zero::

            sage: R.<x> = FiniteField(17)[]
            sage: H = HyperellipticCurveSmoothModel(x^6 + 2, x^2 + 1)
            sage: P = H.point([15,6,1])
            sage: H.is_weierstrass_point(P)
            True
            sage: Q = H.point([3,0,1])
            sage: H.is_weierstrass_point(Q)
            False

        TESTS::

        Check that the examples from the p-adic file work.

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurveSmoothModel(x^3-10*x+9)
            sage: K = Qp(5,8)
            sage: HK = H.change_ring(K)
            sage: P = HK(0,3)
            sage: HK.is_weierstrass_point(P)
            False
            sage: Q = HK(1,0,0)
            sage: HK.is_weierstrass_point(Q)
            True
            sage: S = HK(1,0)
            sage: HK.is_weierstrass_point(S)
            True
            sage: T = HK.lift_x(1+3*5^2); T
            (1 + 3*5^2 + O(5^8) : 3*5 + 4*5^2 + 5^4 + 3*5^5 + 5^6 + O(5^7) : 1 + O(5^8))
            sage: HK.is_weierstrass_point(T)
            False
        """

        f, h = self.hyperelliptic_polynomials()
        if P[2] == 0:
            return self.is_ramified()
        else:
            x, y = self.affine_coordinates(P)
            return y == -y - h(x)

    def rational_weierstrass_points(self):
        """
        Return the rational Weierstrass points of the hyperelliptic curve.
        These are the points that are fixed by the hyperelliptic involution.

        EXAMPLES::

        When `h(x)` is zero, then the Weierstrass points are the points with
        `y`-coordinate equal to zero::

            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurveSmoothModel(x^5 - x)
            sage: H.rational_weierstrass_points()
            [(1 : 0 : 0), (1 : 0 : 1), (0 : 0 : 1), (-1 : 0 : 1)]

        The function also handles the case with `h(x)` nonzero::

            sage: R.<x> = FiniteField(17)[]
            sage: H = HyperellipticCurveSmoothModel(x^6 + 2, x^2 + 1)
            sage: H.rational_weierstrass_points()
            [(15 : 6 : 1), (2 : 6 : 1)]
            sage: P = H.point([15,6,1])
            sage: H.is_weierstrass_point(P)
            True
        """
        f, h = self.hyperelliptic_polynomials()

        F = h**2 + 4 * f
        affine_weierstrass_points = [
            self(r, -h(r) / 2) for r in F.roots(multiplicities=False)
        ]

        if self.is_ramified():  # the point at infinity is Weierstrass
            return self.points_at_infinity() + affine_weierstrass_points
        else:
            return affine_weierstrass_points

    def hyperelliptic_involution(self, P):
        """
        Returns the image of P under the hyperelliptic involution.

        EXAMPLES::

            sage: R.<x> = FiniteField(17)[]
            sage: H = HyperellipticCurveSmoothModel(x^6 + 2, x^2 + 1)
            sage: P = H.point([8,12])
            sage: P_inv = H.hyperelliptic_involution(P); P_inv
            (8 : 8 : 1)
            sage: H.hyperelliptic_involution(P_inv) == P
            True
            sage: Q = H.point([15,6])
            sage: H.is_weierstrass_point(Q)
            True
            sage: H.hyperelliptic_involution(Q) == Q
            True
        """
        [X, Y, Z] = P._coords
        f, h = self.hyperelliptic_polynomials()
        if Z == 0:
            if self.is_ramified():
                return P
            else:
                points = self.points_at_infinity()
                if P == points[0]:
                    return points[1]
                else:
                    return points[0]
        elif Z == 1:
            Y_inv = -Y - h(X)
            return self.point([X, Y_inv])
        else:
            raise ValueError("the point P has to be normalized")

    def distinguished_point(self):
        """
        Return the distinguished point of the hyperelliptic curve.
        By default, this is one of the points at infinity if possible.
        """
        if hasattr(self, "_distinguished_point"):
            return self._distinguished_point

        if not self.is_inert():
            # For the the split and ramified case, a point at infinity is chosen,
            self._distinguished_point = self.points_at_infinity()[0]
            return self._distinguished_point

        assert (
            self.base_ring().characteristic() > 0
        ), "in characteristic 0, a distinguished_point needs to be specified"

        # in the inert case we choose a point with minimal x-coordinate
        for x0 in self.base_ring():
            try:
                self._distinguished_point = self.lift_x(x0)
                return self._distinguished_point
            except ValueError:
                pass

        raise ValueError("distinguished point not found")

    def set_distinguished_point(self, P0):
        """
        Change the distinguished point of the hyperelliptic curve to P0.
        """
        try:
            P0 = self.point(P0)
        except ValueError:
            raise TypeError("P0 must be a point on the hyperelliptic curve")
        self._distinguished_point = P0

    @cached_method
    def jacobian(self):
        """
        Returns the Jacobian of the hyperelliptic curve.
        """
        from sage.schemes.hyperelliptic_curves_smooth_model.jacobian_generic import (
            HyperellipticJacobian_generic,
        )

        return HyperellipticJacobian_generic(self)

    @cached_method
    def projective_curve(self):
        """
        Returns a singular plane model of the hyperelliptic curve self.

        TODO: renaming to plane_model ?

        EXAMPLES::

        We consider the hyperelliptic curve with affine equation
         `y^2 = x^5  + x` ::

            sage: R.<x> = FiniteField(11)[]
            sage: H = HyperellipticCurveSmoothModel(x^6 + 2)
            sage: C = H.projective_curve(); C
            Projective Plane Curve over Finite Field of size 11 defined by -x^6 + y^2*z^4 - 2*z^6


        Note that the projective coordinates of points on H and their images in C
        are in general not the same::

            sage: P = H.point([9,4,2])
            sage: Q = C.point([9,4,2])
            Traceback (most recent call last):
            ...
            TypeError: Coordinates [10, 2, 1] do not define a point on Projective Plane Curve over Finite Field of size 11 defined by -x^6 + y^2*z^4 - 2*z^6


        However, the affine coordinates coincide::

            sage: H.affine_coordinates(P)
            (10, 6)
            sage: Q = C.point([10,6,1])

        The model C has one singular point at infinity, while H is non-singular
        and has two points at infinity.

            sage: H.points_at_infinity()
            [(1 : 1 : 0), (1 : 10 : 0)]
            sage: [P for P in C.rational_points() if P[2]==0]
            [(0 : 1 : 0)]
        """
        from sage.schemes.curves.constructor import Curve

        f, h = self._hyperelliptic_polynomials
        R, (_, y, z) = PolynomialRing(self.base_ring(), 3, "x, y, z").objgens()
        return Curve(R(y**2 + h * y - f).homogenize(z))

    def rational_points(self, **kwds):
        r"""
        Find rational points on the hyperelliptic curve. Arguments are passed
        on to :meth:`sage.schemes.generic.algebraic_scheme.rational_points`.

        ALGORITHM:

        We use :meth:`points_at_infinity` to compute the points at infinity, and
        `sage.schemes.generic.algebraic_scheme.rational_points` on this curve's
        :meth:`projective_curve` for the affine points.

        EXAMPLES::

        For the LMFDB genus 2 curve `932.a.3728.1 <https://www.lmfdb.org/Genus2Curve/Q/932/a/3728/1>`_::

            sage: R.<x> = PolynomialRing(QQ)
            sage: C = HyperellipticCurveSmoothModel(R([0, -1, 1, 0, 1, -2, 1]), R([1]))
            sage: C.rational_points(bound=8)
            [(1 : 1 : 0), (1 : -1 : 0), (-1 : -3 : 1), (-1 : 2 : 1), (0 : -1 : 1),
             (0 : 0 : 1), (1/2 : -5/8 : 1), (1/2 : -3/8 : 1), (1 : -1 : 1), (1 : 0 : 1)]

        Check that :issue:`29509` is fixed for the LMFDB genus 2 curve
        `169.a.169.1 <https://www.lmfdb.org/Genus2Curve/Q/169/a/169/1>`_::

            sage: C = HyperellipticCurveSmoothModel(R([0, 0, 0, 0, 1, 1]), R([1, 1, 0, 1]))
            sage: C.rational_points(bound=10) # long time (6s)
            [(1 : 0 : 0), (1 : -1 : 0), (-1 : 0 : 1), (-1 : 1 : 1), (0 : -1 : 1), (0 : 0 : 1)]

         An example over a number field::

            sage: R.<x> = PolynomialRing(QuadraticField(2))                             # needs sage.rings.number_field
            sage: C = HyperellipticCurveSmoothModel(R([1, 0, 0, 0, 0, 1]))              # needs sage.rings.number_field
            sage: C.rational_points(bound=2)
            [(1 : 0 : 0), (-1 : 0 : 1), (0 : -1 : 1), (0 : 1 : 1), (1 : -a : 1), (1 : a : 1)]

            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurveSmoothModel(x^8 + 33)
            sage: H.rational_points(bound=20)                                           # long time (6s)
            [(1 : 1 : 0), (1 : -1 : 0), (-2 : -17 : 1), (-2 : 17 : 1), (2 : -17 : 1), (2 : 17 : 1)]
        """
        proj_pts = self.projective_curve().rational_points(**kwds)
        return self.points_at_infinity() + [self(*P) for P in proj_pts if P[2] != 0]

    # -------------------------------------------
    # Hacky functions from old implementation.
    # -------------------------------------------

    def is_singular(self, *args, **kwargs):
        r"""
        Returns False, because hyperelliptic curves are smooth projective
        curves, as checked on construction.

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: H = HyperellipticCurveSmoothModel(x^5 + 1)
            sage: H.is_singular()
            False
        """
        return False

    def is_smooth(self):
        r"""
        Returns True, because hyperelliptic curves are smooth projective
        curves, as checked on construction.

        EXAMPLES::

            sage: R.<x> = GF(13)[]
            sage: H = HyperellipticCurveSmoothModel(x^8 + 1)
            sage: H.is_smooth()
            True
        """
        return True

    # -------------------------------------------
    # Odd degree model functions
    # -------------------------------------------

    def odd_degree_model(self):
        r"""
        Return an odd degree model of self, or raise ValueError if one does not
        exist over the field of definition. The term odd degree model refers to
        a model of the form y^2 = f(x) with deg(f) = 2*g + 1.

        EXAMPLES::

            sage: x = QQ['x'].gen()
            sage: H = HyperellipticCurveSmoothModel((x^2 + 2)*(x^2 + 3)*(x^2 + 5)); H
            Hyperelliptic Curve over Rational Field defined by y^2 = x^6 + 10*x^4 + 31*x^2 + 30
            sage: H.odd_degree_model()
            Traceback (most recent call last):
            ...
            ValueError: No odd degree model exists over field of definition

            sage: K2 = QuadraticField(-2, 'a')                                          # needs sage.rings.number_field
            sage: Hp2 = H.change_ring(K2).odd_degree_model(); Hp2                       # needs sage.rings.number_field
            Hyperelliptic Curve over Number Field in a
             with defining polynomial x^2 + 2 with a = 1.414213562373095?*I
             defined by y^2 = 6*a*x^5 - 29*x^4 - 20*x^2 + 6*a*x + 1

            sage: K3 = QuadraticField(-3, 'b')                                          # needs sage.rings.number_field
            sage: Hp3 = H.change_ring(QuadraticField(-3, 'b')).odd_degree_model(); Hp3  # needs sage.rings.number_field
            Hyperelliptic Curve over Number Field in b
             with defining polynomial x^2 + 3 with b = 1.732050807568878?*I
             defined by y^2 = -4*b*x^5 - 14*x^4 - 20*b*x^3 - 35*x^2 + 6*b*x + 1

            Of course, ``Hp2`` and ``Hp3`` are isomorphic over the composite
            extension.  One consequence of this is that odd degree models
            reduced over "different" fields should have the same number of
            points on their reductions.  43 and 67 split completely in the
            compositum, so when we reduce we find:

            sage: # needs sage.rings.number_field
            sage: P2 = K2.factor(43)[0][0]
            sage: P3 = K3.factor(43)[0][0]
            sage: Hp2.change_ring(K2.residue_field(P2)).frobenius_polynomial()
            x^4 - 16*x^3 + 134*x^2 - 688*x + 1849
            sage: Hp3.change_ring(K3.residue_field(P3)).frobenius_polynomial()
            x^4 - 16*x^3 + 134*x^2 - 688*x + 1849

            sage: H.change_ring(GF(43)).odd_degree_model().frobenius_polynomial()       # needs sage.rings.finite_rings
            x^4 - 16*x^3 + 134*x^2 - 688*x + 1849

            sage: # needs sage.rings.number_field
            sage: P2 = K2.factor(67)[0][0]
            sage: P3 = K3.factor(67)[0][0]
            sage: Hp2.change_ring(K2.residue_field(P2)).frobenius_polynomial()
            x^4 - 8*x^3 + 150*x^2 - 536*x + 4489
            sage: Hp3.change_ring(K3.residue_field(P3)).frobenius_polynomial()
            x^4 - 8*x^3 + 150*x^2 - 536*x + 4489

            sage: H.change_ring(GF(67)).odd_degree_model().frobenius_polynomial()       # needs sage.rings.finite_rings
            x^4 - 8*x^3 + 150*x^2 - 536*x + 4489

            The case where `h(x)` is nonzero is also supported::

            sage: HyperellipticCurveSmoothModel(x^5 + 1, 1).odd_degree_model()
            Hyperelliptic Curve over Rational Field defined by y^2 = 4*x^5 + 5
        """
        from sage.schemes.hyperelliptic_curves_smooth_model.hyperelliptic_constructor import (
            HyperellipticCurveSmoothModel,
        )

        f, h = self._hyperelliptic_polynomials
        if f.base_ring().characteristic() == 2:
            raise ValueError(
                "There are no odd degree models over a field with characteristic 2."
            )
        if h:
            f = 4 * f + h**2  # move h to the right side of the equation
        if f.degree() % 2:
            # already odd
            return HyperellipticCurveSmoothModel(f, 0)

        rts = f.roots(multiplicities=False)
        if not rts:
            raise ValueError("No odd degree model exists over field of definition")
        rt = rts[0]
        x = f.parent().gen()
        f_new = f((x * rt + 1) / x).numerator()  # move rt to "infinity"

        return HyperellipticCurveSmoothModel(f_new, 0)

    def has_odd_degree_model(self):
        r"""
        Return True if an odd degree model of self exists over the field of definition; False otherwise.

        Use ``odd_degree_model`` to calculate an odd degree model.

        EXAMPLES::

            sage: x = QQ['x'].0
            sage: HyperellipticCurveSmoothModel(x^5 + x).has_odd_degree_model()
            True
            sage: HyperellipticCurveSmoothModel(x^6 + x).has_odd_degree_model()
            True
            sage: HyperellipticCurveSmoothModel(x^6 + x + 1).has_odd_degree_model()
            False
        """
        try:
            return bool(self.odd_degree_model())
        except ValueError:
            return False

    # -------------------------------------------
    # Magma
    # -------------------------------------------

    def _magma_init_(self, magma):
        """
        Internal function. Returns a string to initialize this elliptic
        curve in the Magma subsystem.

        EXAMPLES::

            sage: # optional - magma
            sage: R.<x> = QQ[]; C = HyperellipticCurveSmoothModel(x^3 + x - 1, x); C
            Hyperelliptic Curve over Rational Field
            defined by y^2 + x*y = x^3 + x - 1
            sage: magma(C)
            Hyperelliptic Curve defined by y^2 + x*y = x^3 + x - 1 over Rational Field
            sage: R.<x> = GF(9,'a')[]; C = HyperellipticCurveSmoothModel(x^3 + x - 1, x^10); C     # needs sage.rings.finite_rings
            Hyperelliptic Curve over Finite Field in a of size 3^2
            defined by y^2 + x^10*y = x^3 + x + 2
            sage: D = magma(C); D                                                       # needs sage.rings.finite_rings
            Hyperelliptic Curve defined by y^2 + (x^10)*y = x^3 + x + 2 over GF(3^2)
            sage: D.sage()                                                              # needs sage.rings.finite_rings
            Hyperelliptic Curve over Finite Field in a of size 3^2
            defined by y^2 + x^10*y = x^3 + x + 2
        """
        f, h = self._hyperelliptic_polynomials
        return f"HyperellipticCurve({f._magma_init_(magma)}, {h._magma_init_(magma)})"

    # -------------------------------------------
    # monsky washnitzer things...
    # -------------------------------------------

    def monsky_washnitzer_gens(self):
        """
        Compute the generators of the special hyperelliptic quotient ring

        EXAMPLES::

            TODO
        """
        from sage.schemes.hyperelliptic_curves_smooth_model import monsky_washnitzer

        S = monsky_washnitzer.SpecialHyperellipticQuotientRing(self)
        return S.gens()

    def invariant_differential(self):
        """
        Returns `dx/2y`, as an element of the Monsky-Washnitzer cohomology
        of ``self``.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: C = HyperellipticCurveSmoothModel(x^5 - 4*x + 4)
            sage: C.invariant_differential()
            1 dx/2y
        """
        from sage.schemes.hyperelliptic_curves_smooth_model import (
            monsky_washnitzer as m_w,
        )

        S = m_w.SpecialHyperellipticQuotientRing(self)
        MW = m_w.MonskyWashnitzerDifferentialRing(S)
        return MW.invariant_differential()

    # -------------------------------------------
    # Local coordinates
    # -------------------------------------------

    def local_coordinates_at_nonweierstrass(self, P, prec=20, name="t"):
        """
        For a non-Weierstrass point `P = (a,b)` on the hyperelliptic
        curve `y^2 + y*h(x) = f(x)`, return `(x(t), y(t))` such that `(y(t))^2 = f(x(t))`,
        where `t = x - a` is the local parameter.

        INPUT:

        - ``P = (a, b)`` -- a non-Weierstrass point on self
        - ``prec`` --  desired precision of the local coordinates
        - ``name`` -- gen of the power series ring (default: ``t``)

        OUTPUT:

        `(x(t),y(t))` such that `y(t)^2 + y(t)*h(x(t)) = f(x(t))` and `t = x - a`
        is the local parameter at `P`.

        EXAMPLES::

        We compute the local coordinates of `H : y^2 = x^5 - 23*x^3 + 18*x^2 + 40*x` at
        the point `P = (1, 6)`::

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurveSmoothModel(x^5 - 23*x^3 + 18*x^2 + 40*x)
            sage: P = H(1, 6)
            sage: xt, yt = H.local_coordinates_at_nonweierstrass(P, prec=5)
            sage: (xt, yt)
            (1 + t + O(t^5),  6 + t - 7/2*t^2 - 1/2*t^3 - 25/48*t^4 + O(t^5))

        We verify that `y(t) = f(x(t))`::

            sage: f,_ = H.hyperelliptic_polynomials()
            sage: (yt^2 - f(xt)).is_zero()
            True

        We can also compute the local coordinates of points on a hyperelliptic curve
        with equation `y^2 + h(x)*y = f(x)`::

            sage: H = HyperellipticCurveSmoothModel(x^6+3*x^5+6*x^4+7*x^3+6*x^2+3*x+1, x^2+x) # 196.a.21952.1
            sage: P = H(-1,1)
            sage: xt, yt = H.local_coordinates_at_nonweierstrass(P, prec=5)
            sage: (xt, yt)
            (-1 + t + O(t^5), 1 - t + 3/2*t^2 - 3/4*t^3 + O(t^5))
            sage: f,h = H.hyperelliptic_polynomials()
            sage: (yt^2 + h(xt)*yt -f(xt)).is_zero()
            True

        AUTHOR:

        - Jennifer Balakrishnan (2007-12)
        - Sabrina Kunzweiler, Gareth Ma, Giacomo Pope (2024): adapt to smooth model
        """
        if P[2] == 0:
            raise TypeError(
                f"P = {P} is a point at infinity. Use local_coordinates_at_infinity instead"
            )
        if self.is_weierstrass_point(P):
            raise TypeError(
                f"P = {P} is a Weierstrass point. Use local_coordinates_at_weierstrass instead"
            )
        f, h = self.hyperelliptic_polynomials()
        a, b = self.affine_coordinates(P)

        L = PowerSeriesRing(self.base_ring(), name, default_prec=prec)
        t = L.gen()
        K = PowerSeriesRing(L, "x")
        f = K(f)
        h = K(h)

        ft = f(t + a)
        ht = h(t + a)
        for _ in range((RR(log(prec, 2))).ceil()):
            b = b - (b**2 + b * ht - ft) / (2 * b + ht)
        return t + a + O(t ** (prec)), b + O(t ** (prec))

    def local_coordinates_at_weierstrass(self, P, prec=20, name="t"):
        """
        For a finite Weierstrass point P = (a,b) on the hyperelliptic
        curve `y^2 + h(x)*y = f(x)`, returns `(x(t), y(t))` such that
        `y(t)^2 + h(x(t))*y(t) = f(x(t))`, where `t = y - b` is the local parameter.

        INPUT:

        - ``P`` -- a finite Weierstrass point on self
        - ``prec`` -- desired precision of the local coordinates
        - ``name`` -- gen of the power series ring (default: `t`)

        OUTPUT:

        `(x(t),y(t))` such that `y(t)^2 + h(x(t))*y(t) = f(x(t))` and `t = y - b`
        is the local parameter at `P = (a,b)`.

        EXAMPLES::

        We compute the local coordinates of the Weierstrass point `P = (4,0)`
        on the hyperelliptic curve `y^2 = x^5 - 23*x^3 + 18*x^2 + 40*x`.

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurveSmoothModel(x^5 - 23*x^3 + 18*x^2 + 40*x)
            sage: P = H(4, 0)
            sage: xt, yt = H.local_coordinates_at_weierstrass(P, prec=7)
            sage: xt
            4 + 1/360*t^2 - 191/23328000*t^4 + 7579/188956800000*t^6 + O(t^7)
            sage: yt
            t + O(t^7)

        We verify that `y(t) = f(x(t))`::

            sage: f,_ = H.hyperelliptic_polynomials()
            sage: (yt^2 - f(xt)).is_zero()
            True

        We compute the local coordinates at the Weierstrass point `(1,-1)`
        of the hyperelliptic curve `y^2 + (x^3 + 1)*y = -x^2.

            sage: H = HyperellipticCurveSmoothModel(-x^2, x^3+1)
            sage: P = H(1,-1)
            sage: xt,yt = H.local_coordinates_at_weierstrass(P)
            sage: f,h = H.hyperelliptic_polynomials()
            sage: (yt^2 + h(xt)*yt - f(xt)).is_zero()
            True

        AUTHOR:

            - Jennifer Balakrishnan (2007-12)
            - Francis Clarke (2012-08-26)
            - Sabrina Kunzweiler, Gareth Ma, Giacomo Pope (2024): adapt to smooth model
        """
        if P[2] == 0:
            raise TypeError(
                f"P = {P} is a point at infinity. Use local_coordinates_at_infinity instead"
            )
        if not self.is_weierstrass_point(P):
            raise TypeError(
                f"P = {P} is not a Weierstrass point. Use local_coordinates_at_nonweierstrass instead"
            )

        L = PowerSeriesRing(self.base_ring(), name, default_prec=prec)
        t = L.gen()
        f, h = self.hyperelliptic_polynomials()
        f_prime = f.derivative()
        h_prime = h.derivative()

        a, b = self.affine_coordinates(P)
        yt = (t + b).add_bigoh(prec)
        yt2 = yt**2
        for _ in range(int(log(prec, 2))):
            a = a - (yt2 + yt * h(a) - f(a)) / (yt * h_prime(a) - f_prime(a))
        return (a, yt)

    def local_coordinates_at_infinity_ramified(self, prec=20, name="t"):
        """
        For a hyperelliptic curve with ramified model `y^2 + h(x)*y = f(x)`,
        return `(x(t), y(t))` such that `(y(t))^2 = f(x(t))`, where
        `t = y/x^{g+1}` is the local parameter at the unique (Weierstrass) point
        at infinity.

        TODO/NOTE: In the previous implementation `t = x^g/y` was used.
        This is not a valid parameter on the smooth model, and the output is necessarily different.

        INPUT:

        - ``prec`` -- desired precision of the local coordinates
        - ``name`` -- generator of the power series ring (default: ``t``)

        OUTPUT:

        `(x(t),y(t))` such that `y(t)^2 = f(x(t))` and `t = y/x^{g+1}`
        is the local parameter at infinity

        EXAMPLES::

        We compute the local coordinates at the point at infinity of the
        hyperelliptic curve `y^2 = x^5 - 5*x^2 + 1`.

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurveSmoothModel(x^5 - 5*x^2 + 1)
            sage: xt, yt = H.local_coordinates_at_infinity_ramified(prec=10)
            sage: (xt,yt)
            (t^-2 - 5*t^4 + t^8 - 75*t^10 + O(t^12),
            t^-5 - 15*t + 3*t^5 - 150*t^7 + 90*t^11 + O(t^12))

        We verify that `y(t)^2 = f(x(t))` and `t = y(t)/x(t)^3`::

            sage: f,_ = H.hyperelliptic_polynomials()
            sage: (yt^2 - f(xt)).is_zero()
            True
            sage: yt/xt^3
            t + O(t^15)

        The method also works when `h` is nonzero. We compute the local
        coordinates of the point at infinity of the hyperelliptic curve
        `y^2 + y = x^5 - 9*x^4 + 14*x^3 - 19*x^2 + 11*x - 6`::

            sage: H = HyperellipticCurveSmoothModel(x^5-9*x^4+14*x^3-19*x^2+11*x-6,1)
            sage: f,h = H.hyperelliptic_polynomials()
            sage: xt,yt = H.local_coordinates_at_infinity_ramified()
            sage: (yt^2 + h(xt)*yt - f(xt)).is_zero()
            True

        Note that the point at infinity has to be a Weierstrass point.


        AUTHOR:

        - Jennifer Balakrishnan (2007-12)
        - Sabrina Kunzweiler, Gareth Ma, Giacomo Pope (2024): adapt to smooth model
        """

        if not self.is_ramified():
            raise TypeError(
                "The point at infinity is not a Weierstrass point. Use local_coordinates_at_infinity_split instead!"
            )

        g = self.genus()
        f, h = self.hyperelliptic_polynomials()
        K = LaurentSeriesRing(self.base_ring(), name, default_prec=prec + 2)
        t = K.gen()
        L = PolynomialRing(K, "x")
        x = L.gen()

        # note that P = H(1,0,0)
        yt = x ** (g + 1) * t
        w = yt**2 + h * yt - f
        wprime = w.derivative(x)
        xt = t**-2
        for _ in range((RR(log(prec + 2) / log(2))).ceil()):
            xt = xt - w(xt) / wprime(xt)
        yt = xt ** (g + 1) * t
        return xt + O(t ** (prec + 2)), yt + O(
            t ** (prec + 2)
        )  # TODO: Why the prec+2 ? Not sure if this is adapted in the correct way.

    def local_coordinates_at_infinity_split(self, P, prec=20, name="t"):
        """
        For a point at infinity P = (a:b:0) on a hyperelliptic curve with
        split model `y^2 + h(x)*y = f(x)`,
        return `(x(t), y(t))` such that `(y(t))^2 = f(x(t))`.
        Here  `t = a/x` is the local parameter at P.

        INPUT:

        - ``P`` -- a point at infinity of a self
        - ``prec`` -- desired precision of the local coordinates
        - ``name`` -- generator of the power series ring (default: ``t``)

        OUTPUT:

        `(x(t),y(t))` such that `y(t)^2 = f(x(t))` and `t = y/x^{g+1}`
        is the local parameter at infinity

        EXAMPLES::

        We compute the local coordinates at the point at infinity of the
        hyperelliptic curve ` `

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurveSmoothModel(x^6+4*x^4 + 4*x^2+1)
            sage: P1 = H(1,-1,0)
            sage: xt1,yt1 = H.local_coordinates_at_infinity_split(P1)

        Note the similarity to the local coordinates of the other point at infinity ::

            sage: P2 = H(1,1,0)
            sage: xt2,yt2 = H.local_coordinates_at_infinity_split(P2)
            sage: xt1 == xt2 and yt1 == - yt2
            True

        Similarly, if `h` is nonzero, the relation between the local coordinates at the
        points at infinity is obtained from the hyperelliptic involution::

            sage: H = HyperellipticCurveSmoothModel(-x^5, x^3+x+1)
            sage: f,h = H.hyperelliptic_polynomials()
            sage: P1 = H(1,0,0)
            sage: P2 = H(1,-1,0)
            sage: xt1,yt1 = H.local_coordinates_at_infinity_split(P1)
            sage: xt2,yt2 = H.local_coordinates_at_infinity_split(P2)
            sage: (yt1 + yt2 + h(xt1)).is_zero()
            True

        AUTHOR:

        - Jennifer Balakrishnan (2007-12)
        - Sabrina Kunzweiler, Gareth Ma, Giacomo Pope (2024): adapt to smooth model
        """

        if not self.is_split():
            raise TypeError(
                "The point at infinity is a Weierstrass point. Use local_coordinates_at_infinity_ramified instead!"
            )
        if not P[2] == 0:
            raise TypeError(
                f"P = {P} is not a point at infinity. Use local_coordinates_at_nonweierstrass or local_coordinates_at_weierstrass instead"
            )

        K = LaurentSeriesRing(self.base_ring(), name, default_prec=prec + 2)
        t = K.gen()

        # note that P = H(a,b,0)
        xt = P[0] / t
        f, h = self.hyperelliptic_polynomials()
        ft = f(xt)
        ht = h(xt)
        yt = P[1] / t**3
        for _ in range((RR(log(prec + 2) / log(2))).ceil()):
            yt = yt - (yt**2 + ht * yt - ft) / (2 * yt + ht)
        return xt + O(t ** (prec + 2)), yt + O(t ** (prec + 2))

    def local_coord(self, P, prec=20, name="t"):
        """
        For point `P = (a,b)` on the hyperelliptic curve
        `y^2 + y*h(x) = f(x)`, return `(x(t), y(t))` such that
        `(y(t))^2 + h(x(t))*y(t) = f(x(t))`, where `t` is the local parameter at
        that point.

        INPUT:

        - ``P`` -- a point on self
        - ``prec`` -- desired precision of the local coordinates
        - ``name`` -- generator of the power series ring (default: ``t``)

        OUTPUT:

        `(x(t),y(t))` such that `y(t)^2 + h(x(t))*y(t) = f(x(t))`, where `t`
        is the local parameter at `P`

        EXAMPLES::

        We compute the local coordinates of several points of the curve with
        defining equation `y^2 = x^5 - 23*x^3 + 18*x^2 + 40*x`. ::

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurveSmoothModel(x^5 - 23*x^3 + 18*x^2 + 40*x)
            sage: H.local_coord(H(1 ,6), prec=5)
            (1 + t + O(t^5), 6 + t - 7/2*t^2 - 1/2*t^3 - 25/48*t^4 + O(t^5))
            sage: H.local_coord(H(4, 0), prec=7)
            (4 + 1/360*t^2 - 191/23328000*t^4 + 7579/188956800000*t^6 + O(t^7), t + O(t^7))
            sage: H.local_coord(H(1, 0, 0), prec=5)
            (t^-2 - 23*t^2 + 18*t^4 - 1018*t^6 + O(t^7),
            t^-5 - 69*t^-1 + 54*t - 1467*t^3 + 3726*t^5 + O(t^6))

        We compute the local coordinates of several points of the curve with
        defining equation `y^2 + (x^3 + 1)*y = x^4 + 2*x^3 + x^2 - x`. ::

            sage: H = HyperellipticCurveSmoothModel(x^4+2*x^3+x^2-x,x^3+1)
            sage: H.local_coord(H(1,-3,1), prec=5)
            (1 + t + O(t^5), -3 - 5*t - 3*t^2 - 3/4*t^3 - 3/16*t^4 + O(t^5))
            sage: H.local_coord(H(1,-1,0), prec=5)
            (t^-1 + O(t^7), -t^-3 - t^-1 - 3 + 6*t^2 + 6*t^3 - 12*t^4 - 42*t^5 + O(t^6))

        AUTHOR:

        - Jennifer Balakrishnan (2007-12)
        - Sabrina Kunzweiler, Gareth Ma, Giacomo Pope (2024): adapt to smooth model
        """

        if P[2] == 0:
            if self.is_ramified():
                return self.local_coordinates_at_infinity_ramified(prec, name)
            else:
                return self.local_coordinates_at_infinity_split(P, prec, name)

        elif self.is_weierstrass_point(P):
            return self.local_coordinates_at_weierstrass(P, prec, name)
        else:
            return self.local_coordinates_at_nonweierstrass(P, prec, name)
