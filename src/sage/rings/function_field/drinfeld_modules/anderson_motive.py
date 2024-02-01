r"""
This module provides facilities for manipulation Anderson motives
and computing the associated L-series.

    sage: from sage.rings.function_field.drinfeld_modules.anderson_motive import AndersonMotive


.. RUBRIC:: Construction of Anderson motives

We first need to define an Anderson motive.
This is achieved by passing in the matrix of `\tau` to the
constructor :class:`AndersonMotive`::

    sage: q = 5
    sage: Fq = GF(q)
    sage: K.<theta> = Fq[]
    sage: R.<t> = K[]
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
    [t + 4*theta]

Another twist::

    sage: M(-3)
    Anderson motive defined by the matrix
    [    t/(t^3 + 2*theta*t^2 + 3*theta^2*t + 4*theta^3)     1/(t^3 + 2*theta*t^2 + 3*theta^2*t + 4*theta^3)]
    [theta/(t^3 + 2*theta*t^2 + 3*theta^2*t + 4*theta^3)     1/(t^3 + 2*theta*t^2 + 3*theta^2*t + 4*theta^3)]


.. RUBRIC:: Tensorial constructions

Basic additive and tensorial constructions are implemented.

Direct sums::

    sage: N = M + F(-2)
    sage: N
    Anderson motive defined by the matrix
    [                            t                             1                             0]
    [                        theta                             1                             0]
    [                            0                             0 1/(t^2 + 3*theta*t + theta^2)]

Tensor products::

    sage: M * M
    Anderson motive defined by the matrix
    [    t^2       t       t       1]
    [theta*t       t   theta       1]
    [theta*t   theta       t       1]
    [theta^2   theta   theta       1]

Duals::

    sage: M.dual()
    Anderson motive defined by the matrix
    [    4/(t + 4*theta)     1/(t + 4*theta)]
    [theta/(t + 4*theta)   4*t/(t + 4*theta)]

Symmetric powers::

    sage: M.symmetric_power(2)
    Anderson motive defined by the matrix
    [      t^2       2*t         1]
    [  theta*t t + theta         1]
    [  theta^2   2*theta         1]

Exterior powers::

    sage: N.exterior_power(2)
    Anderson motive defined by the matrix
    [                      t + 4*theta                                 0                                 0]
    [                                0     t/(t^2 + 3*theta*t + theta^2)     1/(t^2 + 3*theta*t + theta^2)]
    [                                0 theta/(t^2 + 3*theta*t + theta^2)     1/(t^2 + 3*theta*t + theta^2)]

As a shortcut, the method :meth:`determinant` computes the maximal
exterior power::

    sage: N.determinant()
    Anderson motive defined by the matrix
    [1/(t + 4*theta)]


.. RUBRIC:: L-series

The L-series assciated to an Anderson motive is computed
thanks to the method :meth:`Lseries`.
This method takes as input a place of `\FF_q[t]` (encoded
either by ``infinity`` for the place at infinity, an element
of `\FF_q` for a rational place ot an irreducible polynomial
for a general finite place) and a precision::

    sage: F(-3).Lseries(infinity, prec=20)
    (4*u^15 + 2*u^19 + O(u^20))*x + 1 + O(u^20)

    sage: F(-3).Lseries(0, prec=20)
    (3*u^18 + O(u^20))*x^2 + (3*u + u^5 + u^17 + O(u^20))*x + 1 + O(u^20)

    sage: F(-3).Lseries(t^2 + t + 1, prec=10)
    (u^2 + 2*u^3 + 3*u^4 + 4*u^6 + 3*u^8 + 3*u^9 + O(u^10))*x^2 + ((3*a + 4) + (3*a + 4)*u + (3*a + 4)*u^2 + (2*a + 1)*u^3 + (a + 3)*u^6 + (a + 3)*u^7 + (4*a + 2)*u^8 + O(u^10))*x + 1 + O(u^10)

In the output:

- the variable `u` corresponds to a uniformizer of the completion
  of `\mathbb F_q[t]` at the given place: when the place is infinity,
  we have `u = 1/t` whereas, when the place is finite given by
  irreducible polynomial `\mathfrak p(t)`, we have `u = \mathfrak p(t)`,

- the variable `a` is the image of `t` in the residue field,

- the variable `x` is the variable of the L-series.

It is possible to evaluate the `L`-series at a given `x` by just
passing in ``x = value``::

    sage: F(-3).Lseries(infinity, prec=100, x=1)
    1 + 4*u^15 + 2*u^19 + 4*u^23 + 4*u^35 + 2*u^39 + 4*u^43 + 4*u^55 + 2*u^59 + 4*u^63 + 4*u^75 + 2*u^79 + 4*u^83 + u^90 + 3*u^94 + 4*u^95 + u^98 + 2*u^99 + O(u^100)

We check that the L-series of a direct sum is the product of the
L-series of the summands::

    sage: N = M(-1) + F(-3)
    sage: N.Lseries(2, prec=20, x=1)
    1 + 3*u + 3*u^3 + 4*u^5 + 3*u^7 + 3*u^8 + 4*u^9 + 3*u^11 + 3*u^12 + 4*u^13 + 3*u^15 + 3*u^16 + 3*u^18 + 3*u^19 + O(u^20)
    sage: M(-1).Lseries(2, prec=20, x=1) * F(-3).Lseries(2, prec=20, x=1)
    1 + 3*u + 3*u^3 + 4*u^5 + 3*u^7 + 3*u^8 + 4*u^9 + 3*u^11 + 3*u^12 + 4*u^13 + 3*u^15 + 3*u^16 + 3*u^18 + 3*u^19 + O(u^20)


TESTS::

    sage: f = M(-7).Lseries(t^2 + t + 1, prec=100, x=1)
    sage: for prec in [10, 20, 50, 80]:
    ....:     assert(f == M(-7).Lseries(t^2 + t + 1, prec=prec, x=1))

::

    sage: f = F(-3).Lseries(1, prec=50)
    sage: for m in range(1, 10):
    ....:     Fm = AndersonMotive(matrix(R, 1, 1, [(t - theta)^m]))
    ....:     assert(f == Fm(-3-m).Lseries(1, prec=50))

::

    sage: for h in range(50):
    ....:     cond1 = h % (q - 1) == 0
    ....:     cond2 = F(-h).Lseries(0, prec=100, x=1) == 0
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

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute

from sage.categories.homset import Hom

from sage.functions.other import ceil
from sage.misc.functional import log
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



def unimodular_completion(F):
    # unimodular completion based on Zhou-Labahn ISSAC 2024
    # (by Vincent Neiger, Jan 2024)

    # input F: m x n matrix over K[x] with unimodular column bases, m <= n
    # output: (n-m) x n matrix C over K[x] such that [[F],[C]] is unimodular

    pring = F.base_ring()

    if F.nrows() == 0:
        return identity_matrix(pring, F.ncols())

    # find column degrees and reverse column-wise
    # (for zero columns put degree 0 instead of -1:
    # F.reverse requires degrees >=0 and zero columns won't be touched)
    cdeg = [max(d,0) for d in F.column_degrees()]
    Frev = F.reverse(row_wise=False, degree=cdeg)

    # compute right kernel basis
    K = Frev.minimal_kernel_basis(row_wise=False, shifts=cdeg)
    kernel_cdeg = K.column_degrees(shifts=cdeg)

    # compute left approximant basis
    approx_shift = [-c for c in cdeg]
    approx_order = [d+1 for d in kernel_cdeg]
    P = K.minimal_approximant_basis(order=approx_order, shifts=approx_shift)

    # select right rows from P
    approx_rdeg = P.row_degrees(shifts=approx_shift)
    subrows = [i for i in range(len(approx_rdeg)) if approx_rdeg[i] > 0]
    C = P.matrix_from_rows(subrows)
    rdeg = [approx_rdeg[i] for i in subrows]

    # reverse, left-multiply by diag(x**rdeg), and return
    # warning: cannot use reverse here, there are negative degrees...
    left_diag = matrix.diagonal([pring.monomial(d) for d in rdeg])
    right_diag = matrix.diagonal([pring.monomial(d) for d in cdeg])
    ring,var = F.base_ring().objgen()
    Crev = matrix(ring, left_diag * C(var**(-1)) * right_diag)

    return Crev


def normalize_place(A, place, infty=True):
    if infty and place is Infinity:
        return place
    if place in A.base_ring():
        return A.gen() - place
    elif place in A:
        place = A(place)
        if place.degree() == 0:
            return A.gen() - place
        if place.is_irreducible():
            return place.monic()
    if infty:
        raise ValueError("place must be Infinity or an irreducible polynomial")
    else:
        raise ValueError("place must an irreducible polynomial")


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
        self._base = base = tau.base_ring()  # it's Fq[theta][t]
        if not isinstance(base, PolynomialRing_general):
            raise NotImplementedError("tau has to be defined over Fq[theta][t]")
        self._t = base.gen()
        self._K = K = base.base_ring()       # it's Fq[theta]
        # TODO: allow coefficients in Fq(theta)
        if not isinstance(K, PolynomialRing_general):
            raise NotImplementedError("tau has to be defined over Fq[theta][t]")
        self._theta = K.gen()
        self._Fq = Fq = K.base_ring()
        self._q = Fq.cardinality()
        self._deg = ZZ(log(self._q, Fq.characteristic()))
        if self._q is Infinity:
            raise ValueError("not defined over a finite field")
        self._t_name = base.variable_name()
        self._theta_name = K.variable_name()
        self._A = PolynomialRing(Fq, self._t_name)

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

        if check:
            self._dettau

    @lazy_attribute
    def _dettau(self):
        det = self._tau.det()
        if det == 0:
            raise ValueError("tau does not define an Anderson motive")
        h = det.degree()
        disc, R = det.quo_rem((self._t - self._theta) ** h)
        if R:
            raise ValueError("tau does not define an Anderson motive")
        return disc[0], det.degree()

    @lazy_attribute
    def _bad_places(self):
        disc, _ = self._dettau
        return [F for F, _ in disc.factor()]

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

    @cached_method
    def _local_maximal_model(self, place):
        q = self._q
        Fq = self._Fq
        r = self.rank()

        F = Fq.extension(place, name='w')
        if place.degree() == 1:
            a = -place[0]
        else:
            a = F.gen()
        A = PolynomialRing(F, self._t_name)
        K = PolynomialRing(F, 'v')
        v = K.gen()
        S = PolynomialRing(K, self._t_name)
        psiA = A.hom([A.gen()], base_map = F.frobenius_endomorphism(-self._deg))

        tau = matrix(S, r, r,
                     [ S([c(v+a) for c in entry.list()]) for entry in self._tau.transpose().list() ])
        disc, _ = self._dettau
        val = val0 = disc.valuation(place)

        # 1st sequence: the N_i
        while val >= q - 1:
            rows = [ [ R([f[e] for f in c.list()])
                       for c in taurow.list() for e in range(q) ]
                     for taurow in tau.rows() ]
            M = matrix(A, rows)
            basis = M.minimal_kernel_basis()
            dim = basis.nrows()
            if dim == 0:
                break
            M = basis.stack(unimodular_completion(basis))
            N = M.inverse_of_unit()
            N = N.parent()([psiA(x) for x in N.list()])
            tau = M * tau * N
            for i in range(dim):
                for j in range(dim):
                    tau[i,j] = S([x >> q-1 for x in tau[i,j].list()])
                for j in range(dim, r):
                    tau[i,j] = S([x >> q for x in tau[i,j].list()])
                    tau[j,i] = S([x << 1 for x in tau[j,i].list()])
            val -= (q-1) * dim

        # 2nd sequence: the L_i
        if val >= q - 1:
            # The following can be improved
            while True:
                rows = [ [ R([f[e] for f in c.list()])
                           for c in taurow.list() for e in range(q-1) ]
                         for taurow in tau.rows() ]
                M = matrix(A, rows)
                if M.is_zero():
                    break
                basis = M.minimal_kernel_basis()
                dim = basis.nrows()
                M = basis.stack(unimodular_completion(basis))
                N = M.inverse_of_unit()
                N = N.parent()([psiA(x) for x in N.list()])
                tau = M * tau * N
                for i in range(dim, r):
                    for j in range(dim, r):
                        tau[i,j] = S([x << q-1 for x in tau[i,j].list()])
                    for j in range(dim):
                        tau[i,j] = S([x << q for x in tau[i,j].list()])
                        tau[j,i] = S([x >> 1 for x in tau[j,i].list()])
                val += (q-1) * (r - dim)
            for i in range(r):
                for j in range(r):
                    tau[i,j] = S([x >> q-1 for x in tau[i,j].list()])
            val -= (q-1) * r

        if val < val0:
            tau0 = matrix(A, r, r,
                          [ A([c(a) for c in entry.list()]) for entry in self._tau.transpose().list() ])
            tau1 = matrix(A, r, r,
                          [ A([c[0] for c in entry.list()]) for entry in tau.list() ])
            return val, (tau0, tau1)
        else:
            return val, None

    def discriminant(self):
        disc = self._K.one()
        for place in self._bad_places:
            val = self._local_maximal_model(place)[0]
            disc *= place ** val
        return disc

    def bad_reduction_places(self):
        return [place for place in self._bad_places if self._local_maximal_model(place)[0]]

    def has_good_reduction(self):
        return not bool(self.bad_reduction_places())

    def has_good_reduction_at(self, place):
        place = normalize_place(self._K, place, infty=False)
        return not bool(self._local_maximal_model(place)[0])

    def dual(self):
        disc, deg = self._dettau
        scalar = self._K(~disc)
        tau = (-1)**deg * scalar * self._tau.adjugate()
        twist = -deg - self._twist
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
        det, deg = self._dettau
        det *= (self._t - self._theta)**deg
        tau = matrix(self._base, 1, 1, [det])
        twist = self.rank() * self._twist
        return self.__class__(tau, twist, check=False, normalize=True)

    def Lseries(self, place, prec, x=None):
        n = self.rank()
        h = -self._twist
        q = self._q
        completion = MorphismToCompletion(self._A, place, prec + 100)
        place = completion.place()
        C = completion.codomain()
        k = completion.residue_field()
        t = completion(self._A.gen())   # t in C
        ktheta = PolynomialRing(k, self._theta_name)
        Ctheta = LaurentSeriesRing(ktheta, name='u')
        theta = Ctheta(ktheta.gen())

        if x is None:
            S = PolynomialRing(C, name='x')
            x = S.gen()
            charpoly = True
        else:
            S = C
            x = completion(x)
            charpoly = False

        # Correction at bad places
        if place is Infinity and self._bad_places:
            raise NotImplementedError("bad places")
        corr_num = S.one()
        corr_denom = S.one()
        for pl in self._bad_places:
            if pl == place(self._theta):
                continue
            _, taus = self._local_maximal_model(pl)
            if taus is None:
                continue
            tau0, tau1 = taus
            d = pl.degree()
            A = tau0.base_ring()
            F = A.base_ring()
            phiA = A.hom([A.gen()], base_map = F.frobenius_endomorphism(self._deg))
            T0 = tau0; T1 = tau1
            for _ in range(1, d):
                tau0 = tau0.parent()([phiA(y) for y in tau0.list()])
                T0 = tau0 * T0
                tau1 = tau1.parent()([phiA(y) for y in tau1.list()])
                T1 = tau1 * T1
            if pl.degree() > 1:
                ell = PolynomialRing(k, name='w').quo(pl)
                Aell = A.change_ring(ell)
                T0 = T0.change_ring(Aell)
                T1 = T1.change_ring(Aell)
                Cell = LaurentSeriesRing(ell, name='u')
                if charpoly:
                    Sell = PolynomialRing(Cell, name='x')
                else:
                    Sell = Cell
            else:
                Cell = C
            T0 = matrix(n, n, [f(Cell(t)) for f in T0.list()])
            T1 = matrix(n, n, [f(Cell(t)) for f in T1.list()])
            scalar = x**d * pl(t) ** (self._twist)
            chi0 = (1 - scalar*T0).det()
            chi1 = (1 - scalar*T1).det()
            if pl.degree() > 1:
                if charpoly:
                    chi0 = S([C([y.lift()[0] for y in z.list()], z.valuation(), z.precision_absolute())
                          for z in chi0.list()])
                    chi1 = S([C([y.lift()[0] for y in z.list()], z.valuation(), z.precision_absolute())
                              for z in chi1.list()])
                else:
                    chi0 = S([y.lift()[0] for y in chi0.list()], chi0.valuation(), chi0.precision_absolute())
                    chi1 = S([y.lift()[0] for y in chi1.list()], chi1.valuation(), chi1.precision_absolute())
            corr_num *= chi0
            corr_denom *= chi1

        # Figure out which precision is needed
        val = 0
        if place is Infinity:
            for entry in self._tau.list():
                val = max(val, entry.degree())  # t-degree
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
                e += entry[i](theta) * t**i
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

        # Computation of the L-series
        L = (1 - x*taudual).det() * corr_num
        if charpoly:
            L //= corr_denom
            L = S([l.add_bigoh(prec) for l in L])
        else:
            L /= corr_denom
            L = L.add_bigoh(prec)
        return L
