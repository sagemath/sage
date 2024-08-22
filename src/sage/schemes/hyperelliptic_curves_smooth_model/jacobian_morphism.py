from sage.groups.generic import order_from_multiple
from sage.misc.cachefunc import cached_method
from sage.rings.integer import Integer
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.schemes.generic.morphism import SchemeMorphism
from sage.structure.element import AdditiveGroupElement

# TODO _richcmp_

class MumfordDivisorClassField(AdditiveGroupElement, SchemeMorphism):
    r"""
    An element of a Jacobian defined over a field, i.e. in
    `J(K) = \mathrm{Pic}^0_K(C)`.
    """
    def __init__(self, parent, u, v, check=True):
        SchemeMorphism.__init__(self, parent)

        if not isinstance(u, Polynomial) or not isinstance(v, Polynomial):
            raise TypeError(f"arguments {u = } and {v = } must be polynomials")

        # TODO:
        # 1. allow elements of the base field as input
        #   (in particular something like (u,v) = (x-alpha, 0))

        # Ensure the divisor is valid
        if check:
            f, h = parent.curve().hyperelliptic_polynomials()
            assert (
                v**2 + v * h - f
            ) % u == 0, f"{u = }, {v = } do not define a divisor on the Jacobian"

            # TODO: should we automatically do reduction here if the degree of u is
            # too large?

        self._parent = parent
        self._u = u
        self._v = v

    def parent(self):
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

    def __repr__(self) -> str:
        return f"({self._u}, {self._v})"

    def is_zero(self):
        return self._u.is_one() and self._v.is_zero()

    def uv(self):
        return (self._u, self._v)

    def __eq__(self, other):
        if not isinstance(other, MumfordDivisorClassField):
            return False

        u1, v1 = self.uv()
        u2, v2 = other.uv()

        return u1 == u2 and v1 == v2

    def __list__(self):
        return list(self._u, self._v)

    def __tuple__(self):
        return tuple(self._u, self._v)

    def __getitem__(self, n):
        return (self._u, self._v)[n]

    def __hash__(self):
        data = (self._u, self._v)
        return hash(data)

    def __bool__(self):
        return not self.is_zero()

    @cached_method
    def order(self):
        n = self.parent().order()
        return order_from_multiple(self, n)

    def degree(self):
        """
        Returns the degree of the affine part of the divisor.
        """
        return self._u.degree()

    def __add__(self, other):
        # Ensure we are adding two divisors
        if not isinstance(other, type(self)):
            raise ValueError("TODO")

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

    def __neg__(self):
        _, h = self.parent().curve().hyperelliptic_polynomials()
        u0, v0 = self.uv()

        if h.is_zero():
            return self._parent(u0, -v0, check=False)

        v1 = (-h - v0) % u0
        return self._parent(u0, v1, check=False)

    def __sub__(self, other):
        return self.__add__(-other)

    def __rsub__(self, other):
        return (-self).__add__(other)

    def __mul__(self, n):
        """ """
        # TODO: is there a better handlings for this?
        if not isinstance(n, (int, Integer)):
            raise ValueError

        if not n:
            return self._parent().zero()

        P = self

        # Handle negative scalars
        if n < 0:
            n = -n
            P = -P

        # Double and Add
        Q = P
        R = self.parent().zero()
        while n > 0:
            if n % 2 == 1:
                R = R + Q
            Q = Q + Q
            n = n // 2
        return R

    __rmul__ = __mul__



# =======================================================================
# Child classes for representation for ramified, inert and split models
# =======================================================================


class MumfordDivisorClassFieldRamified(MumfordDivisorClassField):
    def __init__(self, parent, u, v, check=True):
        if not parent.curve().is_ramified():
            raise TypeError("hyperelliptic curve must be ramified")

        super().__init__(parent, u, v, check=check)


class MumfordDivisorClassFieldInert(MumfordDivisorClassField):
    def __init__(self, parent, u, v, check=True):
        if not parent.curve().is_inert():
            raise TypeError("hyperelliptic curve must be inert")

        if u.degree() % 2:
            raise ValueError(f"mumford coordinate {u} must have even degree")

        g = parent.curve().genus()
        self._n = (g - u.degree()) // 2
        super().__init__(parent, u, v, check=check)

    def __repr__(self):
        return f"({self._u}, {self._v} : {self._n})"

    def __list__(self):
        return list(self._u, self._v, self._n)

    def __tuple__(self):
        return tuple(self._u, self._v, self._n)

    def __hash__(self):
        return hash(tuple(self))


class MumfordDivisorClassFieldSplit(MumfordDivisorClassField):
    def __init__(self, parent, u, v, n=0, check=True):
        if not parent.curve().is_split():
            raise TypeError("hyperelliptic curve must be split")

        # Ensure the weight is set correctly
        g = parent.curve().genus()
        assert 0 <= n <= (g - u.degree())
        self._n = n
        self._m = g - u.degree() - n

        super().__init__(parent, u, v, check=check)

    def __repr__(self):
        return f"({self._u}, {self._v} : {self._n})"

    def is_zero(self):
        g = self._parent.curve().genus()
        if self._n != (g / 2).ceil():
            return False
        return self._u.is_one() and self._v.is_zero()

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return False

        n1, n2 = self._n, other._n

        if n1 != n2:
            return False

        u1, v1 = self.uv()
        u2, v2 = other.uv()

        return u1 == u2 and v1 == v2

    def __list__(self):
        return list(self._u, self._v, self._n)

    def __tuple__(self):
        return tuple(self._u, self._v, self._n)

    def __hash__(self):
        return hash(tuple(self))

    def __add__(self, other):
        r"""
        Follows algorithm 3.7 of

        Efficient Arithmetic on Hyperelliptic Curves With Real Representation
        David J. Mireles Morales (2008)
        https://www.math.auckland.ac.nz/~sgal018/Dave-Mireles-Full.pdf
        """
        # Ensure we are adding two divisors
        if not isinstance(other, type(self)):
            raise ValueError("TODO")

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

    def __neg__(self):
        r"""
        Follows algorithm 3.8 of

        Efficient Arithmetic on Hyperelliptic Curves With Real Representation
        David J. Mireles Morales (2008)
        https://www.math.auckland.ac.nz/~sgal018/Dave-Mireles-Full.pdf
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
