r"""
Normal forms for products of gamma functions

This module rewrites an expression that is a finite product and quotient of
gamma functions, factorials and binomials with rational-linear arguments
``a*n + b`` into a *normal form* built from factorials::

    C * d^n * prod (A_i*n + B_i)!^e_i * prod p_k(n)^eps_k * prod gamma(a_j*n + b_j)^delta_j

where ``C`` is a constant (a rational times rational powers of primes and of
``pi``, times gamma values at rational points in irreducible position),
``d`` is a nonzero rational, the factorial arguments are nonnegative
integers for all ``n >= 0``, and the ``p_k(n)`` are products of consecutive
linear factors.  The remaining gamma factors, if any, are in *irreducible position*:
their shifts have pairwise distinct fractional parts within each slope, are
normalized into the interval `[0, 1)`, and no complete "circle" is present.

For example::

    sage: from sage.symbolic.gamma_products import simplify_gamma
    sage: n = var('n')
    sage: simplify_gamma(gamma(n + 2/3)*gamma(n + 1/3)/gamma(2/3)/gamma(1/3), n)
    (1/27)^n*factorial(3*n)/factorial(n)

A term can be written without gamma functions if and only if its fractional
gamma shifts decompose into complete circles; this is decided by the algorithm
(see :func:`simplify_gamma`), and otherwise the irreducible gamma factors are
returned.

The same machinery evaluates finite products with a symbolic number of
factors exactly, see :func:`factorial_product`::

    sage: k = var('k')
    sage: from sage.symbolic.gamma_products import factorial_product
    sage: f = (54*k^3 + 81*k^2 + 39*k + 6)/(2*k^3 + 9*k^2 + 13*k + 6)
    sage: factorial_product(f, k, 0, n - 1)
    factorial(3*n)/((2*n + 1)*(n + 1)*factorial(n)^3)

Products of more general bodies, whose rational part is possibly hidden, for
example in ratios of factorials, are evaluated by :func:`exact_product`, which
is used by :func:`~sage.calculus.calculus.symbolic_product`::

    sage: from sage.symbolic.gamma_products import exact_product
    sage: exact_product(factorial(k)/factorial(k + 2)*2^k, k, 0, n)
    2^(1/2*n^2 + 1/2*n)/((n + 2)*factorial(n + 1)^2)

ALGORITHM:

The expression ``E`` is hypergeometric, so the ratio ``E(n+1)/E(n)``
reduces, by applying the recurrence of the gamma function, to a rational function.
Factoring the ratio into linear factors ``(A*n + B)^e``, each factor
contributes ``A^n * gamma(n + B/A)^e / gamma(B/A)^e``, so everything becomes
a finite multiset of gammas with slope 1 and rational shift, plus constants.

Fractional shifts are then eliminated with the Gauss multiplication formula:
a complete circle ``{gamma(a*n + b + j/m) : j = 0..m-1}`` collapses to
``gamma(m*a*n + m*b)``, where members with integer argument are factorials and
may be introduced freely.  Every circle modulus divides the common denominator
of the shifts present, and each move strictly decreases the number of
fractional gammas, so the procedure terminates without any search bound.  A
depth-first search with memoization explores all possible moves, so a
factorial form is found whenever one exists.

The reduction of the term ratio to a rational function and its
factorization into linear factors is the starting point of Gosper's
algorithm [Gos1978]_; see also [PWZ1997]_, Chapter 5, and [GG1999]_,
Section 23.4.  The products of consecutive linear factors arising in the
normal form correspond to greatest factorial factorizations [Paule1995]_,
see also [GG1999]_, Section 23.3.  That the three classical formulas
generate all identities between products of gamma functions is folklore.

AUTHORS:

- Martin Rubey (2026-10)
"""

# ****************************************************************************
#  Copyright (C) 2026 Martin Rubey <martin.rubey at tuwien.ac.at>
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import operator
from collections import Counter

from sage.arith.functions import lcm
from sage.arith.misc import divisors, factor
from sage.calculus.calculus import symbolic_sum
from sage.functions.gamma import gamma
from sage.functions.log import exp
from sage.functions.other import binomial, factorial, symbolic_product
from sage.misc.misc_c import prod
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.symbolic.constants import pi
from sage.symbolic.expression import Expression
from sage.symbolic.operators import mul_vararg
from sage.symbolic.ring import SR

# operators of the symbolic factorial and gamma functions; obtained from
# probe expressions, since the functions in sage.functions are wrappers
_t = SR.var('_gamma_products_probe_')
FAC_OP = factorial(_t).operator()
GAMMA_OP = gamma(_t).operator()
EXP_OP = exp(_t).operator()
E_CONST = exp(Integer(1))

# -----------------------------------------------------------------------
# stage 1: rewrite to gammas and reduce the term ratio to linear factors
# -----------------------------------------------------------------------
def _linear_coeff(exp, n):
    r"""
    Return ``(a, b)`` if ``exp == a*n + b`` with ``a, b`` rational, else ``None``.

    INPUT:

    - ``exp`` -- a symbolic expression

    - ``n`` -- the variable

    The name of the variable must be a valid polynomial ring generator name;
    the callers guarantee this.

    EXAMPLES::

        sage: from sage.symbolic.gamma_products import _linear_coeff
        sage: k = var('k')
        sage: _linear_coeff(2*k + 1, k)
        (2, 1)
        sage: _linear_coeff(k^2 - (k - 1)^2, k)
        (2, -1)
        sage: _linear_coeff(k^2, k) is None
        True
        sage: _linear_coeff((k - 1)^3, k) is None
        True
        sage: x = var('x')
        sage: _linear_coeff(x*k + 1, k) is None
        True
    """
    try:
        p = PolynomialRing(QQ, n)(exp)
    except TypeError:
        return None
    if p.degree() > 1:
        return None
    return p[1], p[0]


def _generator_name(v, *exprs):
    r"""
    Return a valid polynomial ring generator name for the variable ``v``,
    together with ``exprs`` with ``v`` substituted.

    The name of ``v`` itself is returned if it is valid; otherwise a fresh
    temporary variable is used, as happens for converted FriCAS expressions,
    whose variables are called like ``_B``.

    INPUT:

    - ``v`` -- a symbolic variable

    - ``*exprs`` -- symbolic expressions in ``v``

    EXAMPLES::

        sage: from sage.symbolic.gamma_products import _generator_name
        sage: k = var('k')
        sage: _generator_name(k, k + 1)
        (k, k + 1)
        sage: w, e = _generator_name(var('_B'), var('_B') + 1)
        sage: bool(w == var('_B')), bool(e == w + 1)
        (False, True)
    """
    try:
        PolynomialRing(QQ, v)
    except ValueError:
        w = SR.temp_var()
        return w, *(e.subs({v: w}) for e in exprs)
    return v, *exprs


def _strip_exponentials(E, n):
    r"""
    Replace each atom ``c^(a*n+b)`` with positive rational base by `1`.

    INPUT:

    - ``E`` -- a symbolic expression

    - ``n`` -- the variable

    OUTPUT: a pair ``(expression, counter)`` where ``counter`` records the
    constant contribution of the stripped atoms as prime exponents.

    The stripped atoms are independent of the variable if and only if the
    total ``n``-slope vanishes for every prime occurring in the bases, in
    other words, if the exponential part of the expression is a constant.
    Otherwise a :exc:`ValueError` is raised, since the ratio of consecutive
    terms of a hypergeometric term is rational, whereas the exponential part
    then grows nontrivially.

    EXAMPLES::

        sage: from sage.symbolic.gamma_products import _strip_exponentials
        sage: k = var('k')
        sage: _strip_exponentials((k + 1)*2^(k + 2)*2^(-k + 1), k)
        (k + 1, Counter({2: 3}))
        sage: _strip_exponentials((k + 1)*4^k*(1/4)^(k + 1), k)
        (k + 1, Counter({2: -2}))

        The slopes of the bases must cancel prime by prime::

            sage: _strip_exponentials((k + 1)*2^(k + 2)*3^(-k + 1), k)
            Traceback (most recent call last):
            ...
            ValueError: expression is not hypergeometric (exponential part
            grows like (2/3)^n)
    """
    slopes = Counter()
    const = Counter()

    def walk(ex, mult):
        op = ex.operator()
        if op is None:
            return ex
        args = ex.operands()
        if op == operator.pow:
            base, exp = args
            try:
                c = QQ(base)
            except TypeError:
                c = None
            if c is not None and c < 0:
                # (-b)^f with f not an integer is not a hypergeometric factor
                try:
                    integral = QQ(exp).is_integer()
                except TypeError:
                    integral = False
                if not integral:
                    raise ValueError("negative base with non-integer exponent: %s" % ex)
            if c is not None and c > 0:
                lin = _linear_coeff(exp, n)
                if lin is not None:
                    a, b = lin
                    for p, kp in factor(c):
                        slopes[p] += kp * a * mult
                        const[p] += kp * b * mult
                    return SR(1)
            # integer exponent: descend, tracking the multiplicity
            try:
                k = QQ(args[1])
                if k in ZZ:
                    return walk(args[0], mult * ZZ(k)) ** args[1]
            except TypeError:
                pass
        return op(*[walk(a, mult) for a in args])

    new = walk(E, 1)
    if any(slopes.values()):
        growth = prod(p ** s for p, s in slopes.items())
        raise ValueError("expression is not hypergeometric "
                         "(exponential part grows like (%s)^n)" % growth)
    return new, const


def _factor_ratio(R, n, constp=None):
    r"""
    Factor a rational function in ``n`` into a unit and linear factors.

    INPUT:

    - ``R`` -- a symbolic expression which is a rational function in ``n``
    - ``n`` -- the variable
    - ``constp`` -- (default: ``None``) a counter of prime exponents
      contributing to the unit

    OUTPUT: a pair ``(unit, factors)`` with ``R = unit * prod (A*n + B)^e``,
    where ``unit`` is a rational and the entries of ``factors`` are triples
    of integers ``A > 0``, ``B`` and ``e``.

    A :exc:`ValueError` is raised if ``R`` does not split into linear factors
    over ``QQ``, or if the unit is irrational.

    EXAMPLES::

        sage: from sage.symbolic.gamma_products import _factor_ratio
        sage: k = var('k')
        sage: _factor_ratio((k + 1)*(k + 2)/k, k)
        (1, [(1, 1, 1), (1, 2, 1), (1, 0, -1)])
    """
    P = PolynomialRing(QQ, n)
    # any failure of the conversion means that the ratio is not a rational
    # function of the variable
    try:
        fr = R.fraction(P)
    except Exception as exc:
        raise ValueError("expression is not hypergeometric (ratio %s)" % R) from exc
    num, den = fr.numerator(), fr.denominator()
    unitp = Counter(constp)
    sign = 1
    factors = []
    for poly, side in ((num, 1), (den, -1)):
        F = poly.factor()
        U = factor(F.unit())
        sign *= U.unit()
        for p, kp in U:
            unitp[p] += kp * side
        # flint factorizations have no constant factors, the content is
        # carried by the unit
        for f, e in F:
            if f.degree() != 1:
                raise ValueError("irrational factor %s in the ratio" % f)
            A, B = ZZ(f[1]), ZZ(f[0])
            if A < 0:
                A, B = -A, -B
            factors.append((A, B, e * side))
    for v in unitp.values():
        if v not in ZZ:
            raise ValueError("expression is not hypergeometric "
                             "(irrational ratio unit)")
    return sign * prod(p ** v for p, v in unitp.items() if v), factors


def _linear_factors(E, n):
    r"""
    Return ``(unit, factors)`` with ``E(n+1)/E(n) = unit * prod (A*n+B)^e``.

    INPUT:

    - ``E`` -- a symbolic expression, a hypergeometric term in ``n``

    - ``n`` -- the variable

    The ratio is reduced to a rational function by applying the recurrence
    of the gamma function (:meth:`~sage.symbolic.expression.Expression.gamma_normalize`)
    and stripping exponential atoms.

    EXAMPLES::

        sage: from sage.symbolic.gamma_products import _linear_factors
        sage: k = var('k')
        sage: _linear_factors(factorial(k + 1)/factorial(k), k)
        (1, [(1, 2, 1), (1, 1, -1)])

    Constants such as ``pi`` are kept unchanged by the conversion::

        sage: _linear_factors(gamma(k + 1/2)^2/pi, k)
        (1/4, [(2, 1, 2)])
        sage: _linear_factors(1/(pi*factorial(k + 2)), k)
        (1, [(1, 3, -1)])
    """
    # to_gamma also converts Pochhammer symbols
    Eg = E.to_gamma()
    R = (Eg.subs({n: n + 1}) / Eg).gamma_normalize()
    # any gamma or factorial atoms surviving the normalization are caught by
    # the conversion of the ratio to a rational function in _factor_ratio
    R, constp = _strip_exponentials(R, n)
    return _factor_ratio(R, n, constp)


def _parse_const_product(expression):
    r"""
    Split a constant symbolic expression into ``(unit, gammas, rest)``.

    INPUT:

    - ``expression`` -- a symbolic expression free of the product variable

    OUTPUT:

    - ``unit`` -- the rational part
    - ``gammas`` -- a counter mapping rational arguments to the exponent of
      their gamma value
    - ``rest`` -- the remaining factors, which are carried along unchanged

    EXAMPLES::

        sage: from sage.symbolic.gamma_products import _parse_const_product
        sage: _parse_const_product(6*gamma(1/3)^2*sqrt(2))
        (6, Counter({1/3: 2}), sqrt(2))

    Fractional powers of positive rationals are split into prime powers, so
    that they combine with the prime powers contributed elsewhere::

        sage: _parse_const_product(6^(3/2))
        (6, Counter(), sqrt(3)*sqrt(2))
    """
    unit = QQ(1)
    cts = Counter()
    rest = SR(1)
    for base, mult in _mult_split(expression):
        op = base.operator()
        if op is None:
            try:
                unit *= QQ(base) ** mult
                continue
            except TypeError:
                pass
        elif op == GAMMA_OP:
            arg = base.operands()[0]
            if arg.operator() is None:
                try:
                    cts[QQ(arg)] += mult
                    continue
                except TypeError:
                    pass
        elif op == operator.pow:
            # a positive rational raised to a nonintegral rational power is
            # split into prime powers, so that they combine with the prime
            # powers contributed elsewhere
            u, v = base.operands()
            try:
                b = QQ(u)
                if b > 0 and QQ(v) not in ZZ:
                    for p, kp in factor(b):
                        rest *= p ** (kp * QQ(v) * mult)
                    continue
            except TypeError:
                pass
        rest *= base ** mult
    return unit, cts, rest


# -----------------------------------------------------------------------
# stage 2: the solver
# -----------------------------------------------------------------------
class _GammaSolver:
    r"""
    Bookkeeping of the multiset of gammas ``gamma(a*n + x)^e`` occurring in a
    normal form, together with the corrections accumulated by the classical
    gamma identities.

    The state is a dictionary mapping pairs ``(a, x)``, consisting of a
    nonnegative integer slope ``a`` and a rational shift ``x``, to the exponent
    ``e`` of ``gamma(a*n + x)``.  Slope 0 entries are constants.
    """

    def __init__(self, n):
        r"""
        Initialize the solver for the variable ``n``.

        INPUT:

        - ``n`` -- the variable

        EXAMPLES::

            sage: from sage.symbolic.gamma_products import _GammaSolver
            sage: k = var('k')
            sage: _GammaSolver(k).state
            Counter()
        """
        self.n = n
        self.state = Counter()     # (slope, shift) -> exponent
        self.d = QQ(1)             # exponential base
        self.unit = QQ(1)           # rational part of the constant
        self.primes = Counter()     # prime -> exponent, for the constant
        self.pi_exp = QQ(0)         # exponent of pi
        self.fact = Counter()      # (A, B) -> exponent of factorial(A*n + B)
        self.polys = []            # (a, x0, k, e): product of (a*n+x0+j), j<k
        self.memo = set()          # states known to have no solution

    def feed(self, factors):
        r"""
        Feed the linear factors of the term ratio into the solver.

        INPUT:

        - ``factors`` -- a list of triples ``(A, B, e)`` of integers, as
          returned by :func:`_factor_ratio`

        EXAMPLES::

            sage: from sage.symbolic.gamma_products import _GammaSolver
            sage: k = var('k')
            sage: s = _GammaSolver(k)
            sage: s.feed([(1, 1, 1), (1, 0, -1)])
            sage: sorted(s.state.items())
            [((0, 0), 1), ((0, 1), -1), ((1, 0), -1), ((1, 1), 1)]
        """
        for A, B, e in factors:
            x = B/A
            self.d *= A ** e
            self.state[(1, x)] += e
            self.state[(0, x)] -= e

    def normalize(self):
        r"""
        Reduce the state within each (slope, residue) class of shifts.

        Each residue class of shifts is combined at its representative, the
        fractional part of the shifts: it survives as a gamma, the difference
        becomes a product of consecutive linear factors.  Integral classes
        become factorials, anchored at the smallest shift, which must be
        positive, otherwise the term has a pole.

        EXAMPLES::

            sage: from sage.symbolic.gamma_products import _GammaSolver
            sage: k = var('k')
            sage: s = _GammaSolver(k); s.feed([(2, 5, 1), (2, 1, 1)])
            sage: s.normalize()
            sage: sorted(s.state.items()), s.polys, s.unit, s.d
            ([((0, 1/2), -2), ((1, 1/2), 2)], [(1, 1/2, 2, 1)], 4/3, 4)

        Integral residue classes whose smallest member is not positive have
        a pole::

            sage: s = _GammaSolver(k); s.feed([(1, -2, 1), (1, 1, 1)])
            sage: s.normalize()
            Traceback (most recent call last):
            ...
            ValueError: gamma argument reaches a pole at -2 (term undefined at some n >= 0)
        """
        # split the entries into the residue classes of their shifts modulo
        # the integers
        groups = {}
        for (a, x), e in self.state.items():
            if not e:
                continue
            groups.setdefault((a, x - x.floor()), []).append((x, e))

        state = Counter()
        for (a, r), lst in groups.items():
            # the representative is the fractional part of the shifts for
            # fractional classes, and the smallest member for integral ones
            x0 = r if r else min(x for x, _ in lst)
            tot = 0
            if not a:
                # accumulate the correction exponents per offset, so that
                # poles and zeros of the members cancel before the value is
                # evaluated
                exps = Counter()
                for x, e in lst:
                    k = ZZ(x - x0)
                    for j in range(min(0, k), max(0, k)):
                        exps[j] += e if k > 0 else -e
                    tot += e
                for j, ej in exps.items():
                    if not ej:
                        continue
                    if ej < 0 and not x0 + j:
                        raise ValueError("gamma argument reaches a "
                                         "pole at %s (term undefined "
                                         "at some n >= 0)" % (x0 + j))
                    self.unit *= (x0 + j) ** ej
            else:
                for x, e in lst:
                    k = ZZ(x - x0)
                    if k > 0:
                        self.polys.append((a, x0, k, e))
                    elif k < 0:
                        self.polys.append((a, x, -k, -e))
                    tot += e
            if not tot:
                continue
            if not r:
                # integral shifts become factorials, anchored at the
                # smallest shift
                if x0 < 1:
                    raise ValueError("gamma argument reaches a pole at %s "
                                     "(term undefined at some n >= 0)" % x0)
                if a:
                    self.fact[(a, ZZ(x0) - 1)] += tot
                else:
                    self.unit *= factorial(ZZ(x0) - 1) ** tot
            else:
                state[(a, x0)] = tot
        self.state = state

    def circles(self, key):
        r"""
        Return all Gauss circles through the state entry ``key``.

        A circle is a tuple ``(m, b, e, consumed, introduced)``, meaning that
        ``prod_{j=0}^{m-1} gamma(a*n + b + j/m)`` is replaced by
        ``gamma(m*a*n + m*b)``, consuming the entries ``consumed`` with total
        exponent ``e`` and introducing the factorials ``factorial(a*n + y - 1)``
        for the integer members ``y``.

        INPUT:

        - ``key`` -- a key of ``self.state`` whose shifts are to be combined

        EXAMPLES::

            sage: from sage.symbolic.gamma_products import _GammaSolver
            sage: k = var('k')
            sage: s = _GammaSolver(k); s.feed([(3, 1, 1), (3, 2, 1)])
            sage: s.circles((1, 1/3))
            [(3, 1/3, 1, [(1, 1/3), (1, 2/3)], [1])]
        """
        a, x0 = key
        q = lcm(x.denominator() for _, x in self.state)
        out = []
        for m in divisors(q):
            if m < 2:
                continue
            for j0 in range(m):
                b = x0 - QQ((j0, m))
                members = [b + QQ((j, m)) for j in range(m)]
                consumed, introduced, ok = [], [], True
                for mem in members:
                    if mem in ZZ:
                        if mem < 1:
                            ok = False
                            break
                        introduced.append(ZZ(mem))
                    else:
                        if not self.state[(a, mem)]:
                            ok = False
                            break
                        consumed.append((a, mem))
                if not ok:
                    continue
                if key != min(consumed):
                    continue          # generated from its smallest member
                mb = m * b
                if mb in ZZ and mb < 1:
                    continue          # would land on a pole
                s0 = self.state[key]
                if all(self.state[c] * s0 > 0 for c in consumed):
                    # the admissible multiplicities have the common sign of
                    # the consumed entries, and magnitudes up to their minimum
                    sgn = 1 if s0 > 0 else -1
                    mag = min(abs(self.state[c]) for c in consumed)
                    out.extend((m, b, e, consumed, introduced)
                               for e in range(sgn * mag, 0, -sgn))
        return out

    def apply_circle(self, mv):
        r"""
        Apply the Gauss circle ``mv``, updating the corrections.

        INPUT:

        - ``mv`` -- a circle, as returned by :meth:`circles`

        EXAMPLES::

            sage: from sage.symbolic.gamma_products import _GammaSolver
            sage: k = var('k')
            sage: s = _GammaSolver(k); s.feed([(3, 1, 1), (3, 2, 1)])
            sage: s.apply_circle(s.circles((1, 1/3))[0])
            sage: sorted(s.state.items()), sorted(s.fact.items()), s.d
            ([((0, 1/3), -1), ((0, 2/3), -1)], [((1, 0), -1), ((3, 0), 1)], 1/3)
        """
        m, b, e, consumed, introduced = mv
        a = consumed[0][0]
        for c in consumed:
            self.state[c] -= e
            if not self.state[c]:
                del self.state[c]
        for y in introduced:
            if not a:
                self.unit *= factorial(y - 1) ** (-e)
            else:
                self.fact[(a, y - 1)] += -e
        mb = m * b
        if mb in ZZ:
            self.fact[(m * a, mb.numerator() - 1)] += e
        else:
            self.state[(m * a, mb)] += e
        # corrections:  (2*pi)^((m-1)/2) * m^(1/2 - m*(a*n+b))
        corr = (m - 1) * e / 2
        self.pi_exp += corr
        self.primes[2] += corr
        for p, k in factor(m):
            self.primes[p] += k * (QQ((1, 2)) - m * b) * e
        if a:
            self.d *= m ** (-m * a * e)

    def solve(self):
        r"""
        Eliminate all fractional gamma shifts, if possible.

        Return ``True`` if the state is free of fractional shifts at return,
        and ``False`` if no sequence of Gauss circles achieves this.

        EXAMPLES::

            sage: from sage.symbolic.gamma_products import _GammaSolver
            sage: k = var('k')
            sage: s = _GammaSolver(k); s.feed([(3, 1, 1), (3, 2, 1)])
            sage: s.solve(), sorted(s.state.items())
            (True, [])
        """
        self.normalize()
        if not self.state:
            return True
        key = tuple(sorted(self.state.items()))
        if key in self.memo:
            return False
        for key0 in sorted(self.state):
            if key0[1] in ZZ:
                continue
            for mv in self.circles(key0):
                saved = (Counter(self.state), Counter(self.fact), self.d,
                         self.unit, Counter(self.primes),
                         self.pi_exp, len(self.polys))
                self.apply_circle(mv)
                if self.solve():
                    return True
                (self.state, self.fact, self.d, self.unit, self.primes,
                 self.pi_exp, npolys) = saved
                del self.polys[npolys:]
        self.memo.add(key)
        return False

    def expression(self, rest_const):
        r"""
        Assemble the normal form from the current state.

        INPUT:

        - ``rest_const`` -- a constant expression multiplied into the result

        EXAMPLES::

            sage: from sage.symbolic.gamma_products import _GammaSolver
            sage: k = var('k')
            sage: s = _GammaSolver(k); s.feed([(2, 1, 1)]); s.solve()
            True
            sage: s.expression(SR(1))
            (1/2)^k*factorial(2*k)/factorial(k)
        """
        n = self.n
        res = rest_const * self.unit
        for p, f in self.primes.items():
            if f:
                res *= p ** f
        if self.pi_exp:
            res *= pi ** self.pi_exp
        if self.d != 1:
            res *= self.d ** n
        for (A, B), e in self.fact.items():
            if e:
                res *= factorial(A * n + B) ** e
        for (a, x0, k, e) in self.polys:
            res *= prod(a * n + x0 + j for j in range(k)) ** e
        for (a, x), e in self.state.items():
            if e:
                res *= gamma(a * n + x) ** e
        return res


# -----------------------------------------------------------------------
# stage 3: public functions
# -----------------------------------------------------------------------
def simplify_gamma(expression, var):
    r"""
    Rewrite a hypergeometric term as a product of factorials, in normal form.

    INPUT:

    - ``expression`` -- a symbolic expression which is a finite product and
      quotient of powers of rationals, of exponentials with rational base and
      rational-linear exponent, and of factorials, binomials, Pochhammer
      symbols and gamma functions with rational-linear argument.  It must be
      a hypergeometric term in ``var`` and must be defined and nonzero at
      ``var = 0``.

    - ``var`` -- the variable


    OUTPUT: an expression equivalent to ``expression`` of the form

    .. MATH::

        C d^n \prod_i (A_i n + B_i)!^{e_i} \prod_k p_k(n)^{\varepsilon_k}
        \prod_j \Gamma(a_j n + b_j)^{\delta_j}

    where `C` is a constant, `d` a nonzero rational, whose sign absorbs the
    sign of negative linear factors, `A_i, B_i` are nonnegative integers,
    `p_k(n)` are products of consecutive linear factors, and the residual
    gamma factors, if there are any, are in irreducible position: the shifts
    with the same slope `a_j` are distinct modulo the integers, and all
    shifts are normalized into the interval `[0, 1)`.

    A :exc:`ValueError` is raised if ``expression`` is not a hypergeometric
    term, does not split into linear factors over ``QQ``, or is zero or
    undefined at ``var = 0``.

    EXAMPLES::

        sage: from sage.symbolic.gamma_products import simplify_gamma
        sage: n = var('n')
        sage: simplify_gamma(binomial(3*n, n), n)
        factorial(3*n)/(factorial(2*n)*factorial(n))
        sage: simplify_gamma(gamma(n + 1/2), n)
        sqrt(pi)*(1/4)^n*factorial(2*n)/factorial(n)
        sage: simplify_gamma(factorial(n)*gamma(n + 1/2)/sqrt(pi), n)
        (1/4)^n*factorial(2*n)

    A lone gamma with fractional shift cannot be reduced::

        sage: simplify_gamma(gamma(n + 2/3)/gamma(2/3), n)
        gamma(n + 2/3)/gamma(2/3)

    but two of them combine, by the Gauss multiplication formula::

        sage: simplify_gamma(gamma(n + 1/3)*gamma(n + 2/3), n)
        2/3*sqrt(3)*pi*(1/27)^n*factorial(3*n)/factorial(n)

    The following term arises in the enumeration of certain trees; the
    denominator stems from a conversion of a product and was left in an
    unexpanded form::

        sage: E = (3^(3*n)*gamma(n + 2/3)*gamma(n + 1/3)
        ....:      /((2*n*gamma(2/3)*gamma(1/3) + gamma(2/3)*gamma(1/3))
        ....:        *factorial(n + 1)*factorial(n)))
        sage: simplify_gamma(E, n)
        factorial(3*n)/((2*n + 1)*(n + 1)*factorial(n)^3)

    TESTS::

        sage: simplify_gamma(4^n/binomial(2*n, n), n)
        4^n*factorial(n)^2/factorial(2*n)
        sage: simplify_gamma(rising_factorial(2*n + 1, n), n)
        factorial(3*n)/factorial(2*n)
        sage: B = SR.var('_B')
        sage: simplify_gamma(binomial(2*B, B), B)
        factorial(2*_B)/factorial(_B)^2
        sage: simplify_gamma(gamma(n - 1/2), n)
        2*sqrt(pi)*(1/4)^n*factorial(2*n)/((2*n - 1)*factorial(n))
        sage: simplify_gamma(factorial(n) + 1, n)
        Traceback (most recent call last):
        ...
        ValueError: expression is not hypergeometric (ratio ...)
        sage: simplify_gamma((-1)^n*factorial(n), n)
        Traceback (most recent call last):
        ...
        ValueError: negative base with non-integer exponent: ...
        sage: simplify_gamma(n*factorial(n), n)
        Traceback (most recent call last):
        ...
        ValueError: the term vanishes at the anchor var=0...
    """
    n = SR(var) if not isinstance(var, str) else SR.var(var)
    if not (isinstance(n, Expression) and n.is_symbol()):
        raise TypeError("the variable must be a symbol")
    E = SR(expression)

    # the variable is bound, so rename it if its name is not a valid
    # polynomial ring generator name, as happens for converted FriCAS
    # expressions, whose variables are called like ``_B``
    n0 = n
    n, E = _generator_name(n0, E)

    unit, factors = _linear_factors(E, n)
    solver = _GammaSolver(n)
    solver.feed(factors)
    solver.d *= unit

    # anchor the constants: split E(0) into a rational part, gamma values at
    # rational points, and the rest
    try:
        v = E.subs({n: 0})
    except Exception as ex:
        raise ValueError("cannot evaluate the term at var=0: %s" % ex)
    if hasattr(v, 'is_infinity') and v.is_infinity():
        raise ValueError("the term has a pole at var=0")
    try:
        zero = bool(v == 0)
    except Exception:
        zero = False
    if zero:
        raise ValueError("the term vanishes at the anchor var=0; "
                         "factor out the polynomial zero first")
    u0, cts, rest = _parse_const_product(v)
    solver.unit *= u0
    for x, e in cts.items():
        solver.state[(0, x)] = solver.state.get((0, x), 0) + e

    solver.solve()
    res = solver.expression(rest)
    return res if n is n0 else res.subs({n: n0})


def factorial_product(f, var, a, b):
    r"""
    Evaluate a finite product with a symbolic bound in factorial normal form.

    INPUT:

    - ``f`` -- a rational function of ``var`` with rational coefficients

    - ``var`` -- the product variable

    - ``a`` -- the lower bound, a rational number

    - ``b`` -- the upper bound, such that the number of factors
      ``b - a + 1`` is a nonnegative integer, or an integer-linear expression
      ``q*var + r`` in a single variable with `q` a positive integer and `r`
      a nonnegative integer


    OUTPUT: the product `\prod_{v=a}^{b} f(v)`, in the normal form described in
    :func:`simplify_gamma`.

    A :exc:`ValueError` is raised if the input is not of the required form.

    EXAMPLES::

        sage: from sage.symbolic.gamma_products import factorial_product
        sage: k, n = var('k n')
        sage: factorial_product(k, k, 1, n)
        factorial(n)
        sage: factorial_product(2*k + 1, k, 0, n - 1)
        (1/2)^n*factorial(2*n)/factorial(n)
        sage: factorial_product(1/(k + 1), k, 1, n)
        1/factorial(n + 1)
        sage: factorial_product(k^2, k, 1, 2*n)
        factorial(2*n)^2

        sage: f = (54*k^3 + 81*k^2 + 39*k + 6)/(2*k^3 + 9*k^2 + 13*k + 6)
        sage: factorial_product(f, k, 0, n - 1)
        factorial(3*n)/((2*n + 1)*(n + 1)*factorial(n)^3)

    Finite products are evaluated exactly::

        sage: factorial_product(k^2 - 1, k, 2, 6)
        302400

    TESTS::

        sage: factorial_product(2^k, k, 0, n - 1)
        Traceback (most recent call last):
        ...
        ValueError: 2^k is not a rational function of k ...
        sage: x = var('x')
        sage: factorial_product(x + k, k, 0, n)
        Traceback (most recent call last):
        ...
        ValueError: k + x is not a rational function of k with rational coefficients

        Converted FriCAS expressions may have variables with names like
        ``_B``, which are not valid polynomial ring generator names::

            sage: B = SR.var('_B')
            sage: factorial_product(k, k, 1, B)
            factorial(_B)
    """
    v = SR(var) if not isinstance(var, str) else SR.var(var)
    if not (isinstance(v, Expression) and v.is_symbol()):
        raise TypeError("the product variable must be a symbol")
    f = SR(f)
    aS, bS = SR(a), SR(b)
    if v in aS.variables() or v in bS.variables():
        raise ValueError("product bounds must not depend on the product variable")

    # the product variable is bound, so rename it if its name is not a valid
    # polynomial ring generator name, as happens for converted FriCAS
    # expressions, whose variables are called like ``_B``
    w, f = _generator_name(v, f)

    m = bS - aS + 1
    mvars = m.variables()

    # check that f is a rational function of the product variable with
    # rational coefficients; any failure of the conversion counts as a
    # rejection
    Pv = PolynomialRing(QQ, w)
    try:
        Pv.fraction_field()(f)
    except Exception as exc:
        raise ValueError("%s is not a rational function of %s with rational "
                         "coefficients" % (f, v)) from exc

    if not mvars:
        try:
            count = QQ(m)
        except TypeError:
            raise ValueError("the number of factors %s is not rational" % m)
        if count not in ZZ or count < 0:
            raise ValueError("the number of factors %s is not a nonnegative "
                             "integer" % count)
        try:
            aq = QQ(a)
        except (TypeError, ValueError):
            raise ValueError("the lower bound %s is not rational" % a)
        try:
            res = prod(QQ(f.subs({w: aq + k})) for k in range(ZZ(count)))
        except TypeError:
            raise ValueError("%s cannot be evaluated exactly over the given "
                             "range" % f)
        if not res:
            raise ValueError("the product vanishes")
        return SR(res)

    if len(mvars) > 1:
        raise ValueError("the number of factors %s is not linear in a single "
                         "variable" % m)
    # the count variable is bound, so rename it if its name is not a
    # valid polynomial ring generator name, as happens for converted
    # FriCAS expressions, whose variables are called like ``_B``
    n0 = mvars[0]
    n, m = _generator_name(n0, m)
    lin = _linear_coeff(m, n)
    if lin is None:
        raise ValueError("the number of factors %s is not linear in a single "
                         "variable with rational coefficients" % m)
    q, r0 = lin
    if q not in ZZ or q < 1:
        raise ValueError("the number of factors grows by %s per step" % q)
    if r0 not in ZZ or r0 < 0:
        raise ValueError("the number of factors is %s at the anchor" % r0)
    q, r0 = ZZ(q), ZZ(r0)
    try:
        aq = QQ(a)
    except (TypeError, ValueError):
        raise ValueError("the lower bound %s is not rational" % a)

    # the term ratio is the product of the last q factors,
    # E(n+1)/E(n) = prod_{j=0}^{q-1} f(a + m(n) + j)
    R = prod(f.subs({w: aq + q*n + r0 + j}) for j in range(q))
    unit, factors = _factor_ratio(R, n)

    solver = _GammaSolver(n)
    solver.feed(factors)
    solver.d *= unit

    # anchor: E(0) is the product of the first r0 factors
    try:
        E0 = prod(QQ(f.subs({w: aq + p})) for p in range(r0))
    except TypeError as ex:
        raise ValueError("cannot evaluate the product at the anchor: %s" % ex)
    if not E0:
        raise ValueError("the product vanishes at the anchor")
    solver.unit *= E0

    solver.solve()
    res = solver.expression(1)
    return res if n is n0 else res.subs({n: n0})


def _factorials_to_gamma(E, conv=None):
    r"""
    Return ``E`` with factorials rewritten as gamma functions.

    INPUT:

    - ``E`` -- a symbolic expression

    - ``conv`` -- (default: ``None``) a dictionary in which the replacements
      are recorded, so that they can be undone later with ``subs``

    EXAMPLES::

        sage: from sage.symbolic.gamma_products import _factorials_to_gamma
        sage: k = var('k')
        sage: conv = {}
        sage: _factorials_to_gamma(factorial(k + 1), conv)
        gamma(k + 2)
        sage: conv
        {gamma(k + 2): factorial(k + 1)}
    """
    op = E.operator()
    if op is None:
        return E
    args = [_factorials_to_gamma(a, conv) for a in E.operands()]
    if op == FAC_OP:
        g = gamma(args[0] + 1)
        if conv is not None:
            conv[g] = E
        return g
    return op(*args)


def _mult_split(E):
    r"""
    Split ``E`` multiplicatively into a list of pairs ``(base, exponent)``
    with integer exponents, such that ``E`` is the product of all
    ``base^exponent``.

    INPUT:

    - ``E`` -- a symbolic expression

    EXAMPLES::

        sage: from sage.symbolic.gamma_products import _mult_split
        sage: k = var('k')
        sage: _mult_split(2*k^2*(k + 1))
        [(k + 1, 1), (k, 2), (2, 1)]
    """
    op = E.operator()
    if op is None:
        return [(E, Integer(1))]
    if op in (mul_vararg, operator.mul):
        return [pair for a in E.operands() for pair in _mult_split(a)]
    if op == operator.pow:
        u, v = E.operands()
        try:
            vi = Integer(v)
        except TypeError:
            return [(E, Integer(1))]
        return [(g, vi * e) for g, e in _mult_split(u)]
    return [(E, Integer(1))]


def exact_product(f, var, a, b):
    r"""
    Evaluate a finite product of a symbolic body exactly, as far as possible.

    INPUT:

    - ``f`` -- a symbolic expression, the body of the product

    - ``var`` -- the product variable

    - ``a`` -- the lower bound of the product

    - ``b`` -- the upper bound of the product, such that the number of
      factors ``b - a + 1`` is a nonnegative integer, or an integer-linear
      expression in a single variable


    The body is first normalized to expose hidden rationality: factorials
    are rewritten as gamma functions, ratios of gamma functions whose
    arguments differ by integers are converted to polynomials, and
    binomials as well as falling and rising factorials with constant second
    argument are expanded.  For example, the body ``factorial(k)/factorial(k+2)``
    is rational in disguise and normalizes to ``1/((k+1)*(k+2))``.

    The normalized body is then split multiplicatively into

    - a rational function of the product variable, which is evaluated in
      the factorial normal form described in :func:`simplify_gamma`,

    - powers ``c^p`` with base ``c`` independent of the product variable and
      polynomial exponent ``p``, which are evaluated as
      ``c^symbolic_sum(p, var, a, b)``; for instance, ``2^k`` over ``k = 0..n``
      evaluates to ``2^(n*(n+1)/2)``,

    - linear factors ``A*var + B`` with rational slope ``A`` and ``B`` free of
      the product variable, which are evaluated in terms of gamma functions,
      for instance ``x + k`` over ``k = 0..n`` evaluates to
      ``gamma(x+n+1)/gamma(x)``,

    - factors independent of the product variable, which are evaluated as
      ``c^(b - a + 1)``,

    - a residual part, which is returned as an unevaluated product in the
      original variable and bounds.

    A :exc:`ValueError` is raised if nothing can be evaluated; in that case
    the caller may want to fall back to another method.

    EXAMPLES::

        sage: from sage.symbolic.gamma_products import exact_product
        sage: k, n = var('k n')
        sage: exact_product(2^k, k, 0, n)
        2^(1/2*n^2 + 1/2*n)
        sage: exact_product(factorial(k)/factorial(k + 2), k, 0, n)
        1/((n + 2)*factorial(n + 1)^2)
        sage: exact_product((3*k + 1)*(3*k + 2)*factorial(k), k, 0, n)
        (1/3)^n*(3*n + 2)*(3*n + 1)*factorial(3*n)*product(factorial(k), k, 0, n)/factorial(n)

        Powers with symbolic base and polynomial exponent, as well as
        factors independent of the product variable, are evaluated::

            sage: x = var('x')
            sage: exact_product(x^k*(k + 1), k, 0, n)
            x^(1/2*n^2 + 1/2*n)*factorial(n + 1)
            sage: exact_product(sqrt(2)^k, k, 0, n)
            2^(1/4*n^2 + 1/4*n)

        The residual part is returned as an unevaluated product::

            sage: exact_product((k + 1)/(k + 2)*cos(k), k, 0, n)
            product(cos(k), k, 0, n)/(n + 2)

        Linear factors with symbolic constant term are evaluated in terms of
        gamma functions::

            sage: exact_product((x + k)/(k + 1), k, 0, n)
            gamma(n + x + 1)/(factorial(n + 1)*gamma(x))

        Nothing can be done for genuinely non-rational bodies::

            sage: exact_product(cos(k), k, 0, n)
            Traceback (most recent call last):
            ...
            ValueError: the product cannot be evaluated
    """
    v = SR(var) if not isinstance(var, str) else SR.var(var)
    if not (isinstance(v, Expression) and v.is_symbol()):
        raise TypeError("the product variable must be a symbol")
    f = SR(f)
    aS, bS = SR(a), SR(b)
    if v in aS.variables() or v in bS.variables():
        raise ValueError("product bounds must not depend on the product variable")
    if aS.is_infinity() or bS.is_infinity():
        raise ValueError("product bounds must be finite")

    # the product variable is bound, so rename it if its name is not a valid
    # polynomial ring generator name, as happens for converted FriCAS
    # expressions, whose variables are called like ``_B``
    w, f = _generator_name(v, f)

    # normalize the body to expose hidden rationality
    conv = {}
    nf = _factorials_to_gamma(f, conv)
    nf = nf.gamma_normalize()

    # split the normalized body into a rational, an exponential, a constant
    # and a residual part
    F = PolynomialRing(QQ, w).fraction_field()
    R = SR(1)
    C = SR(1)
    expo = []
    param = []
    rest = []
    for g, e in _mult_split(nf):
        # the conversion to QQ(w) is used as a membership test of an arbitrary
        # symbolic expression in the rational function field; which of the
        # arithmetic exceptions it raises on failure is an implementation
        # detail of the conversion
        try:
            F(g)
        except (ValueError, TypeError, ArithmeticError):
            pass
        else:
            R *= g**e
            continue
        op = g.operator()
        if op is operator.pow:
            c, p = g.operands()
            if (w not in c.variables() and w in p.variables()
                    and p.is_polynomial(w)):
                expo.append((c, p * e))
                continue
        elif op == EXP_OP:
            p = g.operands()[0]
            if w in p.variables() and p.is_polynomial(w):
                expo.append((E_CONST, p * e))
                continue
        try:
            coeffs = {d: c for c, d in g.coefficients(w)}
        except TypeError:
            coeffs = None
        if coeffs is not None and set(coeffs) <= {0, 1}:
            A = coeffs.get(1, SR(0))
            B = coeffs.get(0, SR(0))
            if A != 0 and w not in B.variables():
                try:
                    A = QQ(A)
                except TypeError:
                    A = None
                if A:
                    # prod_{k=a}^{b} (A*k + B)^e is
                    # A^(e*(b-a+1)) * (gamma(B/A + b + 1)/gamma(B/A + a))^e
                    param.append((A, B/A, e))
                    continue
        if w not in g.variables():
            C *= g**e
        else:
            rest.append(g**e)

    parts = []
    if R != 1:
        # a failure of the exact evaluation of the rational part, for instance
        # because of symbolic bounds, falls back to the caller
        parts.append(factorial_product(R, w, a, b))
    for A, c, e in param:
        parts.append(A**(e*(bS - aS + 1)) * (gamma(c + bS + 1)/gamma(c + aS))**e)
    if expo:
        for c, p in expo:
            parts.append(c**symbolic_sum(p, w, a, b))
    if C != 1:
        parts.append(C**(bS - aS + 1))
    if not parts:
        if rest:
            raise ValueError("the product cannot be evaluated")
        # everything cancelled: the body is multiplicatively trivial
        return SR(1)
    if rest:
        h = prod(rest).subs(conv)
        if w is not v:
            h = h.subs({w: v})
        parts.append(symbolic_product(h, v, a, b))

    res = prod(parts)
    return res
