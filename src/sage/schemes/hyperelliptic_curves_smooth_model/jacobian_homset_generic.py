import itertools

from sage.misc.banner import SAGE_VERSION
from sage.misc.cachefunc import cached_method
from sage.misc.functional import symbolic_prod as product
from sage.misc.prandom import choice
from sage.rings.finite_rings.finite_field_base import FiniteField as FiniteField_generic
from sage.rings.integer import Integer
from sage.rings.polynomial.polynomial_ring import polygen
from sage.schemes.generic.homset import SchemeHomset_points

# TODO: move this
from sage.schemes.hyperelliptic_curves_smooth_model.weighted_projective_point import (
    SchemeMorphism_point_weighted_projective_ring,
)
from sage.structure.element import parent

assert SAGE_VERSION.startswith("10."), "please update to Sage 10"
assert (
    int(SAGE_VERSION.lstrip("10.").split(".")[0]) >= 4
), "please update to Sage 10.4+ (#37118)"


class HyperellipticJacobianHomset(SchemeHomset_points):
    def __init__(self, Y, X, **kwds):
        SchemeHomset_points.__init__(self, Y, X, **kwds)
        self._morphism_element = None

    def _repr_(self) -> str:
        return f"Abelian group of points on {self.codomain()}"

    def _morphism(self, *args, **kwds):
        return self._morphism_element(*args, **kwds)

    def curve(self):
        # It is unlike you will want to use this, as this returns the original curve H even when
        return self.codomain().curve()

    def extended_curve(self):
        # Code from schemes/generic/homset.py
        if "_extended_curve" in self.__dict__:
            return self._extended_curve
        R = self.domain().coordinate_ring()
        if R is not self.curve().base_ring():
            X = self.curve().base_extend(R)
        else:
            X = self.curve()
        self._extended_curve = X
        return X

    @cached_method
    def order(self):
        """
        TODO: currently using lazy methods by calling sage
        """
        return sum(self.extended_curve().frobenius_polynomial())

    @cached_method
    def _curve_frobenius_roots(self):
        r"""
        Return the roots of the charpoly of frobenius on the extended curve.

        EXAMPLES::
        """
        from sage.rings.qqbar import QQbar

        roots = self.extended_curve().frobenius_polynomial().roots(QQbar)
        return [r for r, e in roots for _ in range(e)]

    def cardinality(self, extension_degree=1):
        r"""
        Return `|Jac(C) / \F_{q^n}|`.

        EXAMPLES::

            sage: R.<x> = GF(5)[]
            sage: H = HyperellipticCurveSmoothModel(x^6 + x + 1)
            sage: J = H.jacobian()
            sage: J.count_points(10) == [J.change_ring(GF((5, k))).order() for k in range(1, 11)]
            True
            sage: J2 = J(GF((5, 2)))
            sage: J2.count_points(5) == J.count_points(10)[1::2]
            True
        """
        K = self.extended_curve().base_ring()
        if not isinstance(K, FiniteField_generic):
            raise NotImplementedError(
                "cardinality is only implemented for Jacobians over finite fields"
            )

        return Integer(
            product(1 - r**extension_degree for r in self._curve_frobenius_roots())
        )

    def count_points(self, n=1):
        try:
            n = Integer(n)
        except TypeError:
            raise TypeError("n must be a positive integer")

        if n < 1:
            raise ValueError("n must be a positive integer")

        if n == 1:
            return self.cardinality()

        return [self.cardinality(extension_degree=i) for i in range(1, n + 1)]

    def point_to_mumford_coordinates(self, P):
        """
        On input a point P, return the Mumford coordinates
        of (the affine part of) the divisor [P].
        """
        R, x = self.extended_curve().polynomial_ring().objgen()
        X, Y, Z = P._coords
        if Z == 0:
            return R.one(), R.zero()
        u = x - X
        v = R(Y)
        return u, v

    def __call__(self, *args, check=True):
        r"""
        Return a rational point in the abstract Homset `J(K)`, given:

        1. No arguments or the integer `0`; return `0 \in J`;

        2. A point `P` on `J = Jac(C)`, return `P`;

        3. A point `P` on the curve `H` such that `J = Jac(H)`;
           return `[P - P_0]`, where `P_0` is the distinguished point of `H`.
           By default, `P_0 = \infty`;

        4. Two points `P, Q` on the curve `H` such that `J = Jac(H)`;
           return `[P - Q]`;

        5. Polynomials `(u, v)` such that `v^2 + hv - f \equiv 0 \pmod u`;
           reutrn `[(u(x), y - v(x))]`.

        EXAMPLES:

        First consider a hyperelliptic curve with an odd-degree model,
        hence a unique point at infinity::

            sage: R.<x> = GF(13)[]
            sage: H = HyperellipticCurveSmoothModel(x^7 + x + 1)
            sage: J = Jacobian(H)
            sage: JH = J.point_homset()
            sage: P = H.lift_x(1)
            sage: Q = H.lift_x(2)
            sage: D1 = JH(P); D1
            (x + 12, 4)
            sage: D2 = JH(Q); D2
            (x + 11, 1)
            sage: D = JH(P,Q); D
            (x^2 + 10*x + 2, 8*x + 9)
            sage: D == D1 - D2
            True
            sage: JH(x^2+10*x+2, 8*x+9) == D
            True

        In general, a distinguished point is used to embed points of the curve
        in the Jacobian. This works for general models of hyperelliptic curves::

            sage: R.<x> = PolynomialRing(GF(13))
            sage: H = HyperellipticCurveSmoothModel(2*x^8 + x + 1)
            sage: H.is_inert()
            True
            sage: J = Jacobian(H)
            sage: JH = J.point_homset()
            sage: P = H.lift_x(1)
            sage: D1 = JH(P); D1
            (x^2 + 12*x, 3*x + 12 : 0)

        To understand this output, one needs to look at the distinguished
        point::

            sage: P0 = H.distinguished_point(); P0
            (0 : 1 : 1)
            sage: JH(P) == JH(P,P0)
            True
            sage: JH(P0)
            (1, 0 : 1)

        We may change the distinguished point. Of course, the divisor `[P - Q]`
        does not depend on the choice of the distinguished point::

            sage: Q = H.lift_x(3)
            sage: JH(P,Q)
            (x^2 + 9*x + 3, 4*x + 11 : 0)
            sage: H.set_distinguished_point(H.random_point())
            sage: JH(P)  # random
            (x^2 + 6*x + 6, 10*x + 5 : 0)
            sage: JH(P,Q)
            (x^2 + 9*x + 3, 4*x + 11 : 0)

        TESTS:

        Ensure that field elements are treated as mumford coordinates::

            sage: R.<x> = GF(7)[]
            sage: H = HyperellipticCurveSmoothModel(x^7 - x^2 - 1)
            sage: J = H.jacobian(); J
            Jacobian of Hyperelliptic Curve over Finite Field of size 7 defined by y^2 = x^7 + 6*x^2 + 6
            sage: J(H.lift_x(3))
            (x + 4, 0)
            sage: _ == J(x + 4, 0) == J(x + 4, R(0))
            True

        TODO:

        - Allow sending a field element corresponding to the x-coordinate of a point?

        - Use ``__classcall__`` to sanitise input?

        - Add doctest for base extension
        """
        R = self.extended_curve().polynomial_ring()

        if len(args) == 0 or (len(args) == 1 and args[0] == ()):
            # this returns the incorrect result because subclasses like the inert case may implement
            # .zero differently
            # return self._morphism_element(self, R.one(), R.zero(), check=check)
            return self.zero(check=check)

        if len(args) == 1 and isinstance(args[0], (list, tuple)):
            args = args[0]

        if len(args) == 1:
            P1 = args[0]
            if P1 == 0:
                return self.zero(check=check)
            elif isinstance(P1, self._morphism_element):
                return P1
            elif isinstance(P1, SchemeMorphism_point_weighted_projective_ring):
                args = args + (
                    self.extended_curve().distinguished_point(),
                )  # this case will now be handled below.
            else:
                raise ValueError(
                    "the input must consist of one or two points, or Mumford coordinates"
                )

        if len(args) == 2:
            P1 = args[0]
            P2 = args[1]
            if isinstance(P1, SchemeMorphism_point_weighted_projective_ring) and isinstance(
                P2, SchemeMorphism_point_weighted_projective_ring
            ):
                u1, v1 = self.point_to_mumford_coordinates(P1)
                P2_inv = self.extended_curve().hyperelliptic_involution(P2)
                u2, v2 = self.point_to_mumford_coordinates(P2_inv)
                u, v = self.cantor_composition(u1, v1, u2, v2)
            # This checks whether P1 and P2 can be interpreted as polynomials
            elif R.coerce_map_from(parent(P1)) and R.coerce_map_from(parent(P2)):
                u = R(P1)
                v = R(P2)
            else:
                raise ValueError(
                    "the input must consist of one or two points, or Mumford coordinates"
                )

        if len(args) > 2:
            raise ValueError("at most two arguments are allowed as input")

        return self._morphism_element(self, u, v, check=check)

    def zero(self, check=True):
        """
        Return the zero element of this jacobian homset.
        """
        H = self.extended_curve()
        R = H.polynomial_ring()
        return self._morphism_element(self, R.one(), R.zero(), check=check)

    def __cantor_double_generic(self, u1, v1):
        """
        Efficient cantor composition for doubling an affine divisor

        Returns the Cantor composition of (u1, v1) with (u1, v1)
        together with the degree of the polynomial ``s`` which is
        needed for computing weights for the split and inert models
        """
        f, h = self.extended_curve().hyperelliptic_polynomials()

        # New mumford coordinates
        if h.is_zero():
            s, _, e2 = u1.xgcd(v1 + v1)
            u3 = (u1 // s) ** 2
            v3 = v1 + e2 * (f - v1**2) // s
        else:
            s, _, e2 = u1.xgcd(v1 + v1 + h)
            u3 = (u1 // s) ** 2
            v3 = v1 + e2 * (f - v1 * h - v1**2) // s
        v3 = v3 % u3

        return u3, v3, s.degree()

    def _cantor_composition_generic(self, u1, v1, u2, v2):
        """
        Helper function for the Cantor composition algorithm.

        OUTPUT:

        The Cantor composition of ``(u1, v1)`` with ``(u2, v2)``,
        together with the degree of the polynomial ``s`` which is
        needed for computing the weights for the split and inert models.
        """
        # Collect data from HyperellipticCurve
        H = self.extended_curve()
        f, h = H.hyperelliptic_polynomials()
        g = H.genus()

        # Ensure D1 and D2 are semi-reduced divisors
        assert (
            v1.degree() < u1.degree() and v2.degree() < u2.degree()
        ), "The degree of bi must be smaller than ai"
        assert (
            u1.degree() <= 2 * g + 2 and u2.degree() <= 2 * g + 2
        ), f"The degree of ai must be smaller than 2g+2, {u1.degree()}, {u2.degree()}"

        # Special case: duplication law
        if u1 == u2 and v1 == v2:
            return self.__cantor_double_generic(u1, v1)

        # Step One
        s0, _, e2 = u1.xgcd(u2)
        v1_m_v2 = v1 - v2

        # Special case: when gcd(u0, u1) == 1 we can
        # avoid many expensive steps as we have s = 1
        if s0.is_one():
            u3 = u1 * u2
            v3 = v2 + e2 * u2 * v1_m_v2
            v3 = v3 % u3
            return u3, v3, 0

        # Step Two
        w0 = v1 + v2 + h

        # Another special case, when w0 is zero we skip
        # a xgcd and can return early
        if w0.is_zero():
            u3 = (u1 * u2) // (s0**2)
            v3 = v2 + e2 * v1_m_v2 * (u2 // s0)
            v3 = v3 % u3
            return u3, v3, s0.degree()

        # Step Three
        s, c1, c2 = s0.xgcd(w0)
        u3 = (u1 * u2) // (s**2)
        v3 = v2 + (c1 * e2 * v1_m_v2 * u2 + c2 * (f - h * v2 - v2**2)) // s
        v3 = v3 % u3
        return u3, v3, s.degree()

    def _cantor_reduction_generic(self, u0, v0):
        """
        Helper function for the Cantor composition algorithm.

        OUTPUT:

        The reduced divisor of ``(u0, v0)``.

        NOTE:

        The `u`-coordinate of the output is not necessarily monic. That step is
        delayed to :meth:`cantor_reduction` to save time.
        """
        # Collect data from HyperellipticCurve
        H = self.extended_curve()
        f, h = H.hyperelliptic_polynomials()

        # Compute u' and v'
        u1 = (v0**2 + h * v0 - f) // u0
        v1 = (-h - v0) % u1

        return u1, v1

    def cantor_composition(self, u1, v1, u2, v2):
        """
        Return the Cantor composition of ``(u1, v1)`` and ``(u2, v2)``.

        EXAMPLES::

        TODO
        """
        u3, v3, _ = self._cantor_composition_generic(u1, v1, u2, v2)
        return u3, v3

    def cantor_reduction(self, u0, v0):
        """
        Return the Cantor reduced ``(u0, v0)``.

        EXAMPLES::

        TODO
        """
        return self._cantor_reduction_generic(u0, v0)

    def lift_u(self, u, all=False):
        """
        Return one or all points with given `u`-coordinate.

        This method is deterministic: it returns the same data each time when
        called with the same `u`.

        Currently only implemented for Jacobians over a finite field.

        INPUT:

        * ``u`` -- an element of the base ring of the Jacobian

        * ``all`` -- boolean (default: ``False``); if ``True``, return a
          (possibly empty) list of all points with the given `u`-coordinate; if
          ``False``, return just one point, or raise a ``ValueError`` if there
          are none.

        OUTPUT:

        A point on this Jacobian.

        EXAMPLES::

            sage: R.<x> = GF(1993)[]
            sage: H = HyperellipticCurveSmoothModel(x^5 + x + 1)
            sage: J = H.jacobian()
            sage: P = J.lift_u(x^2 + 42*x + 270); P
            (x^2 + 42*x + 270, 1837*x + 838)
        """
        H = self.extended_curve()
        K = H.base_ring()
        g = H.genus()

        H_is_split = H.is_split()

        if u.is_zero():
            if H_is_split:
                if not all:
                    return self._morphism_element(self, R.one(), R.zero(), 0)
                return [
                    self._morphism_element(self, R.one(), R.zero(), n)
                    for n in range(g + 1)
                ]
            if not all:
                return self.zero()
            return [self.zero()]

        if not isinstance(K, FiniteField_generic) or K.degree() != 1:
            raise NotImplementedError(
                "lift_u is only implemented for Jacobians over a finite field of prime order"
            )

        R = H.polynomial_ring()
        f, h = H.hyperelliptic_polynomials()
        u1, v1 = R.one(), R.zero()

        u = R(u).monic()
        u_factors = u.factor()
        vss = []

        for x, e in u_factors:
            # Solve y^2 + hy - f = 0 mod x

            # TODO: is this the most efficient method? Maybe we should write
            # a helper function which computes y^2 + hy - f = 0 mod x which
            # properly handles trivial cases like when x is linear?
            K_ext = K.extension(modulus=x, names="a")
            y_ext = polygen(K_ext, "y_ext")
            h_ = K_ext(h % x)
            f_ = K_ext(f % x)
            vs = (y_ext**2 + h_ * y_ext - f_).roots(multiplicities=False)
            if len(vs) == 0 and not all:
                raise ValueError(f"no point with u-coordinate {u} on {self}")

            vss.append(vs)

        def postprocess_uv(u1, v1, n=None):
            # This function converts (u, v) into the point being returned, handling split cases
            if H.is_split():
                if n is None or not (0 <= n <= g - u1.degree()):
                    # this is an internal function
                    raise ValueError(
                        f"bug: n must be an integer between 0 and {g - u1.degree()}"
                    )
                return self._morphism_element(self, u1, v1, n, check=False)

            # We need to ensure the degree of u is even
            if H.is_inert():
                if u1.degree() % 2:
                    # TODO: make composition with distinguished_point its own function?
                    P0 = self.extended_curve().distinguished_point()
                    X0, Y0, _ = P0._coords
                    X = R.gen()  # TODO use better variable names in this function
                    _, h = self.extended_curve().hyperelliptic_polynomials()
                    u0 = X - X0
                    v0 = R(-Y0 - h(X0))
                    u1, v1, _ = self._cantor_composition_generic(u1, v1, u0, v0)
                assert not (u1.degree() % 2), f"{u1} must have even degree"

            return self._morphism_element(self, u1, v1, check=False)

        points = []
        for vv in itertools.product(*vss):
            u1, v1 = R.one(), R.zero()
            try:
                for (x, e), v in zip(u_factors, map(R, vv)):
                    for _ in range(e):
                        u1, v1, _ = self._cantor_composition_generic(u1, v1, x, v)
            # v is not rational, so we skip it
            except (ValueError, AttributeError):
                pass

            if not all:
                return postprocess_uv(u1, v1, n=0)

            if H_is_split:
                for n in range(g - u1.degree() + 1):
                    points.append(postprocess_uv(u1, v1, n=n))
            else:
                points.append(postprocess_uv(u1, v1))

        return points

    def _random_element_cover(self, degree=None):
        r"""
        Return a random element from the Jacobian.

        Distribution is not uniformly random, but returns the entire group.
        """
        H = self.extended_curve()
        R = H.polynomial_ring()
        g = H.genus()

        # For the inert case, the genus must be even
        if H.is_inert():
            assert not (g % 2)

        if degree is None:
            degree = (-1, g)

        while True:
            # TODO: i think we can skip this and simply ensure u
            #       is even degree with composition with the distinguished
            #       point?
            # if H.is_inert() and (u.degree() % 2) == 1:
            #     #TODO: better method to sample even degree polynomials
            #     continue
            u = R.random_element(degree=degree, monic=True)
            try:
                return choice(self.lift_u(u, all=True))
            # TODO: better handling rather than looping with try / except?
            except IndexError:
                pass

    def _random_element_rational(self):
        r"""
        Return a random element from the Jacobian. This algorithm is faster
        than :meth:`_random_element_cover` but is **NOT** surjective on the set
        of points.

        EXAMPLES::

            sage: R.<x> = GF(5)[]
            sage: f = x^5 + 2*x^4 + 4*x^3 + x^2 + 4*x + 3
            sage: H = HyperellipticCurveSmoothModel(f)
            sage: J = H.jacobian()
            sage: J.order()
            16
            sage: JH = J.point_homset()
            sage: len(set(JH._random_element_rational() for _ in range(300)))
            8
        """
        H = self.extended_curve()
        g = H.genus()

        # We randomly sample 2g + 1 points on the hyperelliptic curve
        points = [H.random_point() for _ in range(2 * g + 1)]

        # We create 2g + 1 divisors of the form (P) - infty
        divisors = [self(P) for P in points if P[2] != 0]

        # If we happened to only sample the point at infinity, we return this
        # Otherwise we compute the sum of all divisors.
        if not divisors:
            return self.zero()
        return sum(divisors, start=self.zero())

    def random_element(self, fast=False, *args, **kwargs):
        r"""
        Returns a random element from the Jacobian. Distribution is **NOT**
        uniformly random.

        INPUT:

        - ``fast`` -- (boolean, default ``True``) If set to ``True``, a fast
          algorithm is used, but the output is **NOT** guaranteed to cover the
          entire Jacobian. See examples below. If set to ``False``, a slower
          algorithm is used, but covers the entire Jacobian.

        EXAMPLES::

            sage: R.<x> = GF(5)[]
            sage: f = x^5 + 2*x^4 + 4*x^3 + x^2 + 4*x + 3
            sage: H = HyperellipticCurveSmoothModel(f)
            sage: J = H.jacobian()

        This example demonstrates that the ``fast`` algorithm is not
        necessarily surjective on the rational points of the Jacobian::

            sage: JH = J.point_homset()
            sage: len(set(JH.random_element(fast=True) for _ in range(300)))
            8
            sage: len(set(JH.random_element(fast=False) for _ in range(300)))
            16
            sage: J.order()
            16
        """
        if not isinstance(self.base_ring(), FiniteField_generic):
            raise NotImplementedError(
                "random element of Jacobian is only implemented over Finite Fields"
            )

        if fast:
            return self._random_element_rational(*args, **kwargs)
        return self._random_element_cover(*args, **kwargs)

    def points(self):
        """
        Return all points on this Jacobian `J(K)`.

        .. WARNING::

            This code is not efficient at all.
        """
        H = self.extended_curve()
        R = H.polynomial_ring()
        g = H.genus()

        # TODO: after `monic` argument is added, use it
        ss = []
        for u in R.polynomials(max_degree=g):
            if u.is_monic():
                for P in self.lift_u(u, all=True):
                    # UGLY HACK: it keeps overcounting 0 without this...
                    # failing example: y^2 = x^5 + x over F_5
                    if not u.is_one() and list(P)[:2] == [1, 0]:
                        continue
                    ss.append(P)

        if H.is_ramified() or H.is_split():
            assert len(ss) == self.order()
            return ss

        # TODO: remove this
        # failing example: y^2 = 2x^6 + 1 over F_5
        ss = sorted(set(ss))
        assert len(ss) == self.order()
        return ss

    rational_points = points

    # TODO: Remove this after #38566 is merged
    def __iter__(self):
        yield from self.points()
