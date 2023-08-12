r"""
This module provides facilities for manipulation Anderson motives
and computing the associated L-functions.

    sage: from sage.rings.function_field.drinfeld_modules.anderson_motive import AndersonMotive


.. RUBRIC:: Construction of Anderson motives

We first need to define an Anderson motive.
This is achieved by passing in the matrix of `\tau` to the
constructor :class:`AndersonMotive`::

    sage: q = 5
    sage: Fq = GF(q)
    sage: A.<t> = Fq[]
    sage: R.<theta> = A[]
    sage: tau = matrix(2, 2, [t, 1, theta, 1])
    sage: M = AndersonMotive(tau)
    sage: M
    Anderson motive defined by the matrix
    [    t     1]
    [theta     1]

As a shortcut, the trivial Anderson motive can be defined by
just passing in the base ring::

    sage: F = AndersonMotive(R)
    sage: F
    Anderson motive defined by the matrix
    [1]

We can create (positive or negative) twists.
For example, the Carlitz module can be defined as follows::

    sage: C = F(1)
    sage: C
    Anderson motive defined by the matrix
    [4*theta + t]

Another twist::

    sage: M(-3)
    Anderson motive defined by the matrix
    [    t/(4*theta^3 + 3*t*theta^2 + 2*t^2*theta + t^3)     1/(4*theta^3 + 3*t*theta^2 + 2*t^2*theta + t^3)]
    [theta/(4*theta^3 + 3*t*theta^2 + 2*t^2*theta + t^3)     1/(4*theta^3 + 3*t*theta^2 + 2*t^2*theta + t^3)]


.. RUBRIC:: Tensorial constructions

Basic additive and tensorial constructions are implemented.

Direct sums::

    sage: N = M + F(-2)
    sage: N
    Anderson motive defined by the matrix
    [                            t                             1                             0]
    [                        theta                             1                             0]
    [                            0                             0 1/(theta^2 + 3*t*theta + t^2)]

Tensor products::

    sage: M * M
    Anderson motive defined by the matrix
    [    t^2       t       t       1]
    [t*theta       t   theta       1]
    [t*theta   theta       t       1]
    [theta^2   theta   theta       1]

Duals::

    sage: M.dual()
    Anderson motive defined by the matrix
    [      1/(4*theta + t)       4/(4*theta + t)]
    [4*theta/(4*theta + t)       t/(4*theta + t)]

Symmetric powers::

    sage: M.symmetric_power(2)
    Anderson motive defined by the matrix
    [      t^2       2*t         1]
    [  t*theta theta + t         1]
    [  theta^2   2*theta         1]

Exterior powers::

    sage: N.exterior_power(2)
    Anderson motive defined by the matrix
    [                      4*theta + t                                 0                                 0]
    [                                0     t/(theta^2 + 3*t*theta + t^2)     1/(theta^2 + 3*t*theta + t^2)]
    [                                0 theta/(theta^2 + 3*t*theta + t^2)     1/(theta^2 + 3*t*theta + t^2)]

As a shortcut, the method :meth:`determinant` computes the maximal
exterior power::

    sage: N.determinant()
    Anderson motive defined by the matrix
    [1/(4*theta + t)]


.. RUBRIC:: L-functions

The L-function assciated to an Anderson motive is computed
thanks to the method :meth:`Lfunction`.
This method takes as input a place of `\FF_q[t]` (encoded
either by ``infinity`` for the place at infinity, an element
of `\FF_q` for a rational place ot an irreducible polynomial
for a general finite place) and a precision::

    sage: F(-3).Lfunction(infinity, prec=20)
    (4*u^15 + 2*u^19 + O(u^20))*x + 1 + O(u^20)

    sage: F(-3).Lfunction(0, prec=20)
    (3*u^18 + O(u^20))*x^2 + (3*u + u^5 + u^17 + O(u^20))*x + 1 + O(u^20)

    sage: F(-3).Lfunction(t^2 + t + 1, prec=10)
    (u^2 + 2*u^3 + 3*u^4 + 4*u^6 + 3*u^8 + 3*u^9 + O(u^10))*x^2 + ((3*a + 4) + (3*a + 4)*u + (3*a + 4)*u^2 + (2*a + 1)*u^3 + (a + 3)*u^6 + (a + 3)*u^7 + (4*a + 2)*u^8 + O(u^10))*x + 1 + O(u^10)

In the output:

- the variable `u` corresponds to a uniformizer of the completion
  of `\mathbb F_q[t]` at the given place: when the place is infinity,
  we have `u = 1/t` whereas, when the place is finite given by 
  irreducible polynomial `\mathfrak p(t)`, we have `u = \mathfrak p(t)`,

- the variable `a` is the image of `t` in the residue field,

- the variable `x` is the variable of the L-function.

It is possible to evaluate the `L`-function at a given `x` by just
passing in ``x = value``::

    sage: F(-3).Lfunction(infinity, prec=100, x=1)
    1 + 4*u^15 + 2*u^19 + 4*u^23 + 4*u^35 + 2*u^39 + 4*u^43 + 4*u^55 + 2*u^59 + 4*u^63 + 4*u^75 + 2*u^79 + 4*u^83 + u^90 + 3*u^94 + 4*u^95 + u^98 + 2*u^99 + O(u^100)

We check that the L-function of a direct sum is the product of the
L-functions of the summands::

    sage: N = M(-1) + F(-3)
    sage: N.Lfunction(2, prec=20, x=1)
    1 + 3*u + 3*u^3 + 4*u^5 + 3*u^7 + 3*u^8 + 4*u^9 + 3*u^11 + 3*u^12 + 4*u^13 + 3*u^15 + 3*u^16 + 3*u^18 + 3*u^19 + O(u^20)
    sage: M(-1).Lfunction(2, prec=20, x=1) * F(-3).Lfunction(2, prec=20, x=1)
    1 + 3*u + 3*u^3 + 4*u^5 + 3*u^7 + 3*u^8 + 4*u^9 + 3*u^11 + 3*u^12 + 4*u^13 + 3*u^15 + 3*u^16 + 3*u^18 + 3*u^19 + O(u^20)


TESTS::

    sage: f = M(-7).Lfunction(t^2 + t + 1, prec=100, x=1)
    sage: for prec in [10, 20, 50, 80]:
    ....:     assert(f == M(-7).Lfunction(t^2 + t + 1, prec=prec, x=1))

::

    sage: f = F(-3).Lfunction(1, prec=50)
    sage: for m in range(1, 10):
    ....:     Fm = AndersonMotive(matrix(R, 1, 1, [(t - theta)^m]))
    ....:     assert(f == Fm(-3-m).Lfunction(1, prec=50))

::

    sage: for h in range(50):
    ....:     cond1 = h % (q - 1) == 0
    ....:     cond2 = F(-h).Lfunction(0, prec=100, x=1) == 0
    ....:     assert(cond1 == cond2)


AUTHORS:

- Xavier Caruso (2023-08): initial version

"""

# *****************************************************************************
#        Copyright (C) 2023 Xavier Caruso <xavier.caruso@normalesup.org>
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#                   http://www.gnu.org/licenses/
# *****************************************************************************


from sage.structure.sage_object import SageObject
from sage.categories.homset import Hom

from sage.functions.other import ceil
from sage.misc.mrange import mrange
from sage.misc.misc_c import prod

from sage.sets.set import Set

from sage.rings.ring import CommutativeRing
from sage.rings.morphism import RingHomomorphism
from sage.rings.integer_ring import ZZ
from sage.rings.infinity import Infinity
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
from sage.rings.laurent_series_ring import LaurentSeriesRing
from sage.rings.function_field.function_field import FunctionField

from sage.matrix.constructor import matrix
from sage.matrix.special import identity_matrix, block_diagonal_matrix


def normalize_place(A, place):
    if place is Infinity:
        return place
    R = place.parent()
    if A.base_ring().has_coerce_map_from(R):
        return A.gen() - place
    if A.has_coerce_map_from(R):
        place = A(place)
        if place.degree() == 0:
            return A.gen() - place
        if place.is_irreducible():
            return place.monic()
    raise ValueError("place must be Infinity or an irreducible polynomial")


class MorphismToCompletion(RingHomomorphism):
    def __init__(self, domain, place, default_prec):
        if not isinstance(domain, PolynomialRing_general):
            raise NotImplementedError
        k = base = domain.base_ring()
        place = self._place = normalize_place(domain, place)
        if place is Infinity:
            codomain = LaurentSeriesRing(base, name='u', default_prec=default_prec)
            image = codomain.one() >> 1
        elif place.degree() == 1:
            codomain = LaurentSeriesRing(base, name='u', default_prec=default_prec)
            image = codomain.gen() - place[0]
        else:
            k = base.extension(place, name='a')
            codomain = LaurentSeriesRing(k, name='u', default_prec=default_prec)
            image = codomain(k.gen()).add_bigoh(1)
        parent = Hom(domain, codomain)
        RingHomomorphism.__init__(self, parent)
        self._image = image
        self._k = k
        self._q = k.cardinality()

    def degree(self):
        if self._place is Infinity:
            return ZZ(1)
        else:
            return self._place.degree()

    def residue_field(self):
        return self._k

    def image_generator(self, prec):
        S = self.codomain()
        image = self._image
        current_prec = image.precision_absolute()
        if current_prec < prec:
            # We update the image
            u = S.gen()
            place = self._place
            der = place.derivative()
            while current_prec < prec:
                current_prec *= 2
                image = image.lift_to_precision(current_prec)
                image -= (place(image) - u) / der(image)
            self._image = image
        return image.add_bigoh(prec)

    def _repr_type(self):
        return "Completion"

    def place(self):
        return self._place

    def _call_with_args(self, P, args, kwds):
        if 'prec' in kwds:
            prec = kwds['prec']
        else:
            prec = self.codomain().default_prec()
        image = self.image_generator(prec)
        return P(image)

    def _call_(self, P):
        image = self._image
        if image.precision_absolute() is not Infinity:
            prec = self.codomain().default_prec()
            image = self.image_generator(prec)
        return P(image)

    def teichmuller_lift(self, x, prec):
        place = self._place
        q = self._q
        modulus = place
        Q = self._k(x).polynomial()
        if Q == 0:
            return Q
        Q = self.domain()(Q)
        current_prec = 1
        while current_prec < prec:
            modulus = modulus ** 2
            current_prec *= 2
            R = Q ** (q-2)
            G, D, _ = R.xgcd(modulus)
            if G != 1:
                raise RuntimeError
            Q = 2*Q - D
        return Q % (place**prec)


class AndersonMotive(SageObject):  # should be a Parent
    def __init__(self, tau, twist=0, check=True, normalize=True):
        self._twist = ZZ(twist)
        if isinstance(tau, CommutativeRing):
            tau = identity_matrix(tau, 1)
        if tau.nrows() != tau.ncols():
            raise ValueError("tau must be a square matrix")
        self._base = base = tau.base_ring()  # it's Fq[t][theta]
        if not isinstance(base, PolynomialRing_general):
            raise NotImplementedError("tau has to be defined over Fq[t][theta]")
        self._A = A = base.base_ring()       # it's Fq(t)
        if not isinstance(A, PolynomialRing_general):
            raise NotImplementedError("tau has to be defined over Fq[t][theta]")
        self._t = A.gen()
        self._Fq = Fq = A.base_ring()
        self._q = Fq.cardinality()
        if self._q is Infinity:
            raise ValueError("not defined over a finite field")
        self._theta = base.gen()
        self._theta_name = base.variable_name()

        if normalize:
            divisor = self._t - self._theta
            exponent = Infinity
            for entry in tau.list():
                if not entry:
                    continue
                e = 0
                while entry.degree() > 0 and e < exponent:
                    entry, R = entry.quo_rem(divisor)
                    if R:
                        break
                    e += 1
                exponent = e
                if exponent == 0:
                    break
            if exponent is not Infinity and exponent > 0:
                denom = divisor ** exponent
                tau = tau.parent()([entry // denom for entry in tau.list()])
                self._twist += exponent
        self._tau = tau

        self._det = None
        if check:
            det = self._dettau()
            if det == 0:
                raise ValueError("tau does not define an Anderson motive")
            while det.degree() > 0:
                det, R = det.quo_rem(divisor)
                if R:
                    raise ValueError("tau does not define an Anderson motive")
            if det[0].degree() > 0:
                raise NotImplementedError("bad reduction")

    def _dettau(self):
        if self._det is None:
            self._det = self._tau.det()
        return self._det

    def __repr__(self):
        mat = ((self._t - self._theta) ** self._twist) * self._tau
        return "Anderson motive defined by the matrix\n%s" % mat

    def rank(self):
        return self._tau.nrows()

    def twist(self, n):
        n = ZZ(n)
        return self.__class__(self._tau, self._twist + n, check=False, normalize=False)

    def __call__(self, n):
        return self.twist(n)

    def tensor_product(self, other):
        if not isinstance(other, AndersonMotive):
            raise TypeError("cannot compute the tensor product of an Anderson module with something else")
        if self._base is not other._base:
            raise TypeError("Anderson modules are not defined over the same base")
        if self.rank() == 0:
            return self
        if other.rank() == 0:
            return other
        tau = self._tau.tensor_product(other._tau)
        tau.subdivide()
        twist = self._twist + other._twist
        return self.__class__(tau, twist, check=False)

    def __mul__(self, other):
        return self.tensor_product(other)

    def direct_sum(self, other):
        if not isinstance(other, AndersonMotive):
            raise TypeError("cannot compute the direct sum of an Anderson module with something else")
        if self._base is not other._base:
            raise TypeError("Anderson modules are not defined over the same base")
        if self.rank() == 0:
            return other
        if other.rank() == 0:
            return self
        t = self._t
        theta = self._theta
        n = self._twist - other._twist
        if n > 0:
            tau = block_diagonal_matrix([(t - theta)**n * self._tau, other._tau])
            twist = other._twist
        else:
            tau = block_diagonal_matrix([self._tau, (t - theta)**(-n) * other._tau])
            twist = self._twist
        tau.subdivide()
        return self.__class__(tau, twist, check=False, normalize=(n == 0))

    def __add__(self, other):
        return self.direct_sum(other)

    def dual(self):
        # TODO: double-check the formula
        det = self._dettau()
        deg = det.degree()
        scalar = ~(det.leading_coefficient()[0])
        tau = (-1)**deg * scalar * self._tau.adjugate()
        twist = -det.degree() - self._twist
        return self.__class__(tau, twist, check=False, normalize=True)

    def symmetric_power(self, n):
        # Optimization:
        # instead of running over mrange([r] * n),
        # run over E and then consider all permutations
        n = ZZ(n)
        if n < 0:
            raise ValueError("exponent must be nonnegative")
        if n == 0:
            return self.__class__(self._base)
        R = self._base
        r = self.rank()
        S = Set(range(r+n-1))
        E = [ ]
        Ed = { }
        for s in S.subsets(r-1):
            es = (n + r - 2 - s[-1]) * (0,)
            for i in range(r-2, 0, -1):
                es += (s[i] - s[i-1] - 1) * (r-1-i,)
            es += s[0] * (r-1,)
            Ed[es] = len(E)
            E.append(es)
        rows = [ ]
        for es in E:
            row = [R.zero() for i in range(len(E))]
            for fs in mrange([r] * n):
                j = Ed[tuple(sorted(fs))]
                row[j] += prod(self._tau[es[i], fs[i]] for i in range(n))
            rows.append(row)
        tau = matrix(rows)
        twist = n * self._twist
        return self.__class__(tau, twist, check=False, normalize=True)

    def exterior_power(self, n):
        n = ZZ(n)
        if n < 0:
            raise ValueError("exponent must be nonnegative")
        if n == 0:
            return self.__class__(self._base)
        r = self.rank()
        if n > r:
            return self.__class__(identity_matrix(self._base, 0))
        I = Set(range(r))
        In = [ sorted(list(J)) for J in I.subsets(n) ]
        rows = [ ]
        for J1 in In:
            row = [ ]
            for J2 in In:
                M = self._tau.matrix_from_rows_and_columns(J1, J2)
                row.append(M.det())
            rows.append(row)
        tau = matrix(rows)
        twist = n * self._twist
        return self.__class__(tau, twist, check=False, normalize=True)

    def determinant(self):
        tau = matrix(self._base, 1, 1, [self._dettau()])
        twist = self.rank() * self._twist
        return self.__class__(tau, twist, check=False, normalize=True)

    def Lfunction(self, place, prec, x=None):
        n = self.rank()
        h = -self._twist
        q = self._q
        completion = MorphismToCompletion(self._A, place, prec)
        place = completion.place()
        C = completion.codomain()
        k = completion.residue_field()
        ktheta = PolynomialRing(k, self._theta_name)
        Ctheta = LaurentSeriesRing(ktheta, name='u')
        theta = Ctheta(ktheta.gen())

        val = 0
        if place is Infinity:
            for entry in self._tau.list():
                for c in entry.list():
                    val = max(val, c.degree())
            prectau = prec + n*val - min(0, h)
        else:
            prectau = prec

        # Computation of rho
        hi = h
        current_prec = 1
        if place is Infinity:
            v = Ctheta.gen()
            rho = Ctheta(1).add_bigoh(prectau - h)
            while current_prec < prectau - h:
                hj = ceil(hi / q)
                rho *= (1 - v*theta) ** (q*hj - hi)
                hi = hj
                v = v**q
                current_prec *= q
            rho <<= h
        else:
            d = completion.degree()
            v = completion(self._t)
            a = v[0]
            m = 0
            ev = -h
            rho = Ctheta(1).add_bigoh(prectau)
            while current_prec < prectau:
                hj = ceil(hi / q)
                e = q*hj - hi
                rho *= (v - theta) ** e
                ev -= e * q**m
                hi = hj
                v = v ** q
                current_prec *= q
                m = (m + 1) % d
            ev %= q**d - 1
            for _ in range(d):
                ev, e = ev.quo_rem(q)
                rho *= (a - theta) ** (e + q - 1)
                a = a ** q

        # Computation of tau*
        u = Ctheta.gen()
        B = [ ]
        kmax = 1
        for entry in self._tau.list():
            e = Ctheta(0)
            for i in range(entry.degree() + 1):
                e += completion(entry[i]) * theta**i
            e *= rho
            B.append(e)
            d = max([0] + [P.degree() for P in e.list()])
            kmax = max(kmax, d // (q - 1))
        rows = [ ]
        for i in range(n):
            for k in range(kmax):
                row = [ ]
                for j in range(n):
                    s = B[n*i + j]
                    for kp in range(kmax):
                        l = q*kp - k + q - 1
                        coeffs = [ s[e][l] for e in range(val, s.degree() + 1) ]
                        cs = (C(coeffs) << val).add_bigoh(prectau)
                        row.append(cs)
                rows.append(row)
        taudual = matrix(rows)

        # Computation of the L-function
        if x is None:
            S = PolynomialRing(C, name='x')
            L = (1 - S.gen()*taudual).det()
            L = S([l.add_bigoh(prec) for l in L])
        else:
            L = (1 - completion(x)*taudual).det()
            L = L.add_bigoh(prec)
        return L
