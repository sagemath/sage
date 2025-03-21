r"""
Witt vectors

Implementation of the class :class:`WittVector` of truncated Witt vectors.

AUTHORS:

- Jacob Dennerlein (2022-11-28): initial version
- Rubén Muñoz-\-Bertrand (2025-02-13): major refactoring and clean-up
"""

# ****************************************************************************
#       Copyright (C) 2025 Rubén Muñoz--Bertrand
#                          <ruben.munoz--bertrand@univ-fcomte.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.integer import Integer
from sage.rings.padics.factory import Zp, Zq
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.structure.element import CommutativeRingElement
from sage.structure.richcmp import op_EQ, op_NE


def padic_to_vector(x, F):
    r"""
    Return the coordinates in `W(F)` of `x \in \mathbb Z_q`.

    INPUT:

    - ``x`` -- an element in an unramified `p`-adic ring

    - ``F`` -- (a field isomorphic to) the residue field

    EXAMPLES::

        sage: from sage.rings.padics.witt_vector import padic_to_vector
        sage: R.<a> = ZqFM(5^2, prec=5)
        sage: F.<b> = GF(5^2)
        sage: x = (3*a + 4)*5 + (2*a + 1)*5^2 + 5^4
        sage: padic_to_vector(x, F)
        [0, 2*b + 2, 2*b + 3, 4*b, 2*b + 2]
    """
    prec = x.precision_absolute()
    val = x.valuation()
    if x == 0:
        return prec * [F.zero()]
    p = x.parent().prime()
    E = list(x.teichmuller_expansion())
    return [F(E[i].residue().polynomial()) ** (p**i) for i in range(prec)]


def vector_to_padic(v, R):
    r"""
    Return the element `x \in R` corresponding to the element
    of `W(\mathbb F_q)` whose coordinates are `v`.

    INPUT:

    - ``v`` -- a tuple of elements in a finite field

    - ``R`` -- (a ring isomorphic to) the corresponding `p`-adic ring

    EXAMPLES::

        sage: from sage.rings.padics.witt_vector import vector_to_padic
        sage: R.<a> = ZqFM(5^2, prec=5)
        sage: F.<b> = GF(5^2)
        sage: v = [F(0), 2*b + 2, 2*b + 3, 4*b, 2*b + 2]
        sage: vector_to_padic(v, R)
        (3*a + 4)*5 + (2*a + 1)*5^2 + 5^4
    """
    p = R.prime()
    F = R.residue_field()
    v = [F(c.polynomial()) for c in v]
    return sum(R.teichmuller(v[i].nth_root(p**i)) << i for i in range(len(v)))


class WittVector(CommutativeRingElement):
    """
    Base class for truncated Witt vectors.

    EXAMPLES::

        sage: W = WittVectorRing(GF(25), p=5, prec=3)
        sage: W(12)
        (2, 1, 3)

        sage: W = WittVectorRing(Integers(6), p=3, prec=4)
        sage: w = W([1,2,3,4]) * W([4,5,0,0])
        sage: w
        (4, 1, 3, 4)

        sage: TestSuite(w).run()
    """
    def __init__(self, parent, vec=None):
        """
        Common class for all kinds of Witt vectors.

        EXAMPLES::

            sage: W = WittVectorRing(GF(3))
            sage: e = W.one(); e
            (1)
            sage: e^2
            (1)
            sage: -e
            (2)

            sage: W = WittVectorRing(GF(3), prec=4)
            sage: t = W([1,2,0,1])
            sage: t^2
            (1, 1, 0, 2)
            sage: -t
            (2, 1, 0, 2)
            sage: 1/t
            (1, 1, 1, 0)

            sage: W=WittVectorRing(ZZ, p=5, prec=2)
            sage: WW=WittVectorRing(ZZ, p=5, prec=2)
            sage: W((4,10)) * WW((-5,12,1))
            (-20, -18362)
            sage: WW((1,2,3)) + W((1,2))
            (2, -2)
        """
        self._prec = parent.precision()
        B = parent.base()
        if vec is not None:
            if isinstance(vec, int) or isinstance(vec, Integer):
                self._int_to_vector(vec, parent)
            elif isinstance(vec, WittVector):
                if vec.parent().precision() < self._prec:
                    raise ValueError(f'{vec} has not the correct length. '
                                     'Expected length has to be at least '
                                     f'{self._prec}.')
                if not B.has_coerce_map_from(vec.parent().base()):
                    raise ValueError('Cannot coerce an element of '
                                     f'{vec.base()} to an element of {B}.')
                self._vec = tuple(B(vec.vec()[i]) for i in range(self._prec))
            elif isinstance(vec, tuple) or isinstance(vec, list):
                if len(vec) < self._prec:
                    raise ValueError(f'{vec} has not the correct length. '
                                     'Expected length has to be at least '
                                     f'{self._prec}.')
                self._vec = tuple(B(vec[i]) for i in range(self._prec))
            else:
                raise ValueError(f'{vec} cannot be interpreted as a Witt '
                                 'vector.')
        else:
            self._vec = (B(0) for i in range(self._prec))
        CommutativeRingElement.__init__(self, parent)

    def __hash__(self) -> int:
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: W = WittVectorRing(GF(3), prec=4)
            sage: t = W([1,2,0,1])
            sage: hash(t)  # random
            -2438844084280889141
        """
        return hash(self._vec)

    def _richcmp_(self, other, op) -> bool:
        """
        Compare ``self`` and ``other``.

        EXAMPLES::

            sage: W = WittVectorRing(GF(3), prec=4)
            sage: t = W([1,2,0,1])
            sage: u = 1/t
            sage: u == t
            False
            sage: t != t
            False
        """
        if not isinstance(other, WittVector):
            return NotImplemented
        if op == op_EQ:
            return self._vec == other.vec()
        if op == op_NE:
            return self._vec != other.vec()
        return NotImplemented

    def _repr_(self) -> str:
        """
        Return a string representation.

        EXAMPLES::

            sage: W = WittVectorRing(GF(3), prec=4)
            sage: t = W([1,2,0,1]); t
            (1, 2, 0, 1)
        """
        return '(' + ', '.join(map(str, self._vec)) + ')'

    def _add_(self, other):
        """
        Return the sum of ``self`` and ``other``.

        EXAMPLES::

            sage: W = WittVectorRing(GF(3), prec=4)
            sage: t = W([1,2,0,1])
            sage: u = 1/t
            sage: u + t
            (2, 1, 1, 1)
        """
        P = self.parent()

        # As a slight optimization, we'll check for zero ahead of time.
        if other == P.zero():
            return self
        elif self == P.zero():
            return other

        alg = P.algorithm()
        if alg == 'standard':
            s = P.sum_polynomials()
            # note here this is tuple addition, i.e. concatenation
            sum_vec = tuple(s[i](*(self._vec + other.vec()))
                            for i in range(self._prec))

        elif alg == 'finotti':
            x = self._vec
            y = other.vec()
            prec = P.precision()
            G = []
            for n in range(prec):
                G_n = [x[n], y[n]]
                for i in range(n):
                    G_n.append(P._eta_bar(G[i], n - i))
                G.append(G_n)
            sum_vec = tuple(sum(G[i]) for i in range(self._prec))

        elif alg == 'Zq_isomorphism':
            x = vector_to_padic(self._vec, P._padic_ring)
            y = vector_to_padic(other.vec(), P._padic_ring)
            sum_vec = padic_to_vector(x + y, P.base_ring())

        elif alg == 'p_invertible':
            p = P.prime()  # we know p is a unit in this case!
            x = self._vec
            y = other.vec()
            sum_vec = [x[0] + y[0]]
            for n in range(1, self._prec):
                next_sum = x[n] + y[n] + \
                    sum((x[i]**(p**(n - i)) + y[i]**(p**(n - i))
                         - sum_vec[i]**(p**(n - i)))
                        / p**(n - i)
                        for i in range(n))
                sum_vec.append(next_sum)

        return P(sum_vec)

    def _mul_(self, other):
        """
        Return the product of ``self`` and ``other``.

        EXAMPLES::

            sage: W = WittVectorRing(GF(3), prec=4)
            sage: t = W([1,2,0,1])
            sage: u = 1/t + 1
            sage: u * t
            (2, 0, 0, 1)
        """
        P = self.parent()

        # As a slight optimization, we'll check for zero or one ahead of time.
        if self == P.zero() or other == P.zero():
            return P.zero()
        if other == P.one():
            return self
        if self == P.one():
            return other

        alg = P.algorithm()
        if alg == 'standard':
            p = P.prod_polynomials()
            # note here this is tuple addition, i.e. concatenation
            prod_vec = tuple(p[i](*(self._vec + other.vec()))
                             for i in range(self._prec))

        elif alg == 'finotti':
            from sage.rings.padics.witt_vector_ring import fast_char_p_power
            x = self._vec
            y = other.vec()
            prec = P.precision()
            p = P.prime()
            G = [[x[0] * y[0]]]
            for n in range(1, prec):
                G_n = [fast_char_p_power(x[0], p**n) * y[n],
                       fast_char_p_power(y[0], p**n) * x[n]]
                G_n.extend(fast_char_p_power(x[i], p**(n - i)) * fast_char_p_power(y[n - i], p**i)
                           for i in range(1, n))
                for i in range(n):
                    G_n.append(P._eta_bar(G[i], n - i))
                G.append(G_n)
            prod_vec = tuple(sum(G[i]) for i in range(prec))

        elif alg == 'Zq_isomorphism':
            x = vector_to_padic(self._vec, P._padic_ring)
            y = vector_to_padic(other.vec(), P._padic_ring)
            prod_vec = padic_to_vector(x * y, P.base_ring())

        elif alg == 'p_invertible':
            p = P.prime()  # we know p is a unit in this case!
            x = self._vec
            y = other.vec()
            prod_vec = [x[0] * y[0]]
            for n in range(1, self._prec):
                next_prod = (
                    sum(p**i * x[i]**(p**(n - i)) for i in range(n + 1)) *
                    sum(p**i * y[i]**(p**(n - i)) for i in range(n + 1)) -
                    sum(p**i * prod_vec[i]**(p**(n - i)) for i in range(n))
                ) / p**n
                prod_vec.append(next_prod)

        return P(prod_vec)

    def _neg_(self):
        """
        Return the opposite of ``self``.

        EXAMPLES::

            sage: W = WittVectorRing(GF(3), prec=4)
            sage: t = W([1,2,0,1])
            sage: -t
            (2, 1, 0, 2)
        """
        P = self.parent()
        # If p == 2, -1 == (-1, -1, -1, ...)
        # Otherwise, -1 == (-1, 0, 0, ...)
        if P.prime() == 2:
            all_ones = P(tuple(-1 for _ in range(self._prec)))
            return all_ones * self
        neg_vec = tuple(-self._vec[i] for i in range(self._prec))
        return P(neg_vec)

    def _div_(self, other):
        """
        Return the quotient of ``self`` and ``other``.

        EXAMPLES::

            sage: W = WittVectorRing(GF(3), prec=4)
            sage: t = W([1,2,0,1])
            sage: u = 1/t + 1
            sage: u / t
            (2, 1, 2, 1)
        """
        P = self.parent()
        # As a slight optimization, we'll check for one ahead of time.
        if other == P.one():
            return self
        elif self == P.one():
            return other._invert_()

        return self * other._invert_()

    def _invert_(self):
        """
        Return the inverse of ``self``.

        EXAMPLES::

            sage: W = WittVectorRing(GF(3), prec=4)
            sage: t = W([1,2,0,1])
            sage: ~t
            (1, 1, 1, 0)
        """
        if not self._vec[0].is_unit():
            raise ZeroDivisionError(f"Inverse of {self} does not exist.")
        P = self.parent()

        if self == P.one():
            return self
        if self._prec == 1:
            return P((self._vec[0]**-1, ))

        # Strategy: Multiply ``self`` by ``(Y_0, Y_1, ...)``, set equal
        # to (1, 0, 0, ...), and solve.
        var_names = [f'Y{i}' for i in range(1, self._prec)]
        poly_ring = PolynomialRing(P.base(), var_names)
        inv_vec = list((self._vec[0]**-1,) + poly_ring.gens())
        # We'll fill this in one-by-one

        from sage.rings.padics.witt_vector_ring import WittVectorRing
        W = WittVectorRing(poly_ring, p=P.prime(), prec=P.precision())
        prod_vec = (W(self._vec) * W(inv_vec)).vec()
        for i in range(1, self._prec):
            poly = prod_vec[i](inv_vec[1:])
            Y_i = poly.parent().gens()[i - 1]
            try:
                inv_vec[i] = (-poly.constant_coefficient()
                              / poly.monomial_coefficient(Y_i))
            except ZeroDivisionError:
                raise ZeroDivisionError(f"Inverse of {self} does not exist.")
            try:
                inv_vec[i] = P.base()(inv_vec[i])
            except ValueError:
                raise ZeroDivisionError(f"Inverse of {self} does not exist.")

        return P(inv_vec)

    def _int_to_vector(self, k, parent):
        """
        Return the image of ``k`` in ``self`` with coefficients in ``parent``.

        EXAMPLES::

            sage: W = WittVectorRing(ZZ, p=23, prec=2)
            sage: W(-123)
            (-123, 50826444131062300759362981690761165250849615528)
        """
        p = parent.prime()
        R = parent.base()

        if p == R.characteristic():
            Z = Zp(p, prec=self._prec + 1, type='fixed-mod')
            self._vec = padic_to_vector(Z(k), R)

        should_negate = False
        if k < 0:
            k = -k
            should_negate = True

        vec_k = [k]
        for n in range(1, self._prec):
            total = (
                k - k**(p**n)
                - sum(p**(n-i) * vec_k[n-i]**(p**i) for i in range(1, n))
            )
            total //= p**n
            vec_k.append(total)

        if should_negate:
            if p == 2:
                vec_k = (
                    parent(vec_k)
                    * parent((tuple(-1 for _ in range(self._prec))))
                ).vec()
            else:
                vec_k = (-x for x in vec_k)

        self._vec = tuple([R(x) for x in vec_k])

    def vec(self):
        """
        Return the underlying tuple of the truncated Witt vector.

        EXAMPLES::

            sage: W = WittVectorRing(GF(7), p=7, prec=3)
            sage: v = W([1,2,3])
            sage: v.vec()
            (1, 2, 3)
        """
        return self._vec


class WittVector_phantom(WittVector):
    def __init__(self, parent, vec=None, phantom=None):
        self._prec = prec = parent.precision()
        R = parent.base_ring()
        p = parent._prime
        mod = parent._mod
        lift = parent._lift
        if phantom is not None:
            self._phantom = phantom
            self._vec = [mod(phantom[0])]
            self._powers = [phantom[0]]
        elif vec is None:
            zero = R.zero()
            zerolift = lift(zero)
            self._vec = prec * [zero]
            self._phantom = prec * [zero]
        elif isinstance(vec, WittVector_phantom):
            self._vec = vec._vec
            self._phantom = vec._phantom
            self._powers = vec._powers
        elif isinstance(vec, int) or isinstance(vec, Integer):
            y = lift(vec)
            self._vec = [mod(y)]
            self._powers = [y]
            self._phantom = prec * [y]
        elif isinstance(vec, tuple) or isinstance(vec, list):
            if len(vec) < self._prec:
                raise ValueError(f'{vec} has not the correct length. '
                                 'Expected length has to be at least '
                                 f'{self._prec}.')
            # We compute the phantom components
            self._vec = [R(v) for v in vec]
            x = [lift(v) for v in self._vec]
            self._phantom = [x[0]]
            for n in range(1, prec):
                for i in range(n):
                    x[i] = x[i] ** p
                self._phantom.append(sum(x[i] * p**i for i in range(n+1)))
            self._powers = None
        else:
            raise ValueError(f'{vec} cannot be interpreted as a Witt vector.')
        CommutativeRingElement.__init__(self, parent)

    def _compute_vector(self, prec=None):
        maxprec = self._prec
        if prec is None:
            prec = maxprec
        else:
            prec = min(prec, maxprec)
        phantom = self._phantom
        powers = self._powers
        vec = self._vec
        p = self.parent()._prime
        mod = self.parent()._mod
        for n in range(len(vec), prec):
            for i in range(n):
                powers[i] = powers[i] ** p
            c = (phantom[n] - sum(powers[i] * p**i for i in range(n))) // p**n
            self._vec.append(mod(c))
            self._powers.append(c)

    def __repr__(self):
        self._compute_vector()
        return str(tuple(self._vec))

    def phantom(self):
        return self._phantom

    def __getitem__(self, i):
        if i < 0 or i >= self._prec:
            raise IndexError
        self._compute_vector(i+1)
        return self._vec[i]

    def _add_(self, other):
        phantom = [self._phantom[i] + other._phantom[i] for i in range(self._prec)]
        return self.__class__(self.parent(), phantom=phantom)

    def _sub_(self, other):
        phantom = [self._phantom[i] - other._phantom[i] for i in range(self._prec)]
        return self.__class__(self.parent(), phantom=phantom)

    def _neg_(self):
        phantom = [-v for v in self._phantom]
        return self.__class__(self.parent(), phantom=phantom)

    def _mul_(self, other):
        phantom = [self._phantom[i] * other._phantom[i] for i in range(self._prec)]
        return self.__class__(self.parent(), phantom=phantom)
