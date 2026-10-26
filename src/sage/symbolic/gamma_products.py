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

Infinite sums of hypergeometric terms are evaluated in hypergeometric
notation by :func:`hypergeometric_sum`::

    sage: k, a = var('k a')
    sage: from sage.symbolic.gamma_products import hypergeometric_sum
    sage: hypergeometric_sum(rising_factorial(a, k)*x^k/factorial(k), k, 0, oo)
    hypergeometric((a,), (), x)

The method :meth:`~sage.symbolic.expression.Expression.to_hypergeometric`
applies this to every unevaluated sum in an expression::

    sage: k = var('k')
    sage: sum(1/factorial(k), k, 0, oo, hold=True).to_hypergeometric()
    hypergeometric((), (), 1)

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
of the shifts present, and each move strictly decreases the total number of
fractional gamma factors, counted with multiplicity, so the depth-first
search with memoization terminates.  The search is bounded; when no factorial
form is found within the bound, or none exists, a normal form is committed
instead: a state to which no circle applies.

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
from sage.rings.infinity import infinity, minus_infinity
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


def _rational_atom(ex, n):
    r"""
    Return ``(c, a, b)`` if ``ex == c^(a*n + b)`` with positive rational
    ``c`` and rational ``a`` and ``b``, else ``None``.

    INPUT:

    - ``ex`` -- a symbolic expression

    - ``n`` -- the variable

    EXAMPLES::

        sage: from sage.symbolic.gamma_products import _rational_atom
        sage: k = var('k')
        sage: _rational_atom(4^(k + 1), k)
        (4, 1, 1)
        sage: _rational_atom((1/2)^k, k)
        (1/2, 1, 0)
        sage: _rational_atom((-4)^k, k) is None
        True
        sage: _rational_atom(2^(k^2), k) is None
        True
        sage: _rational_atom(k^2, k) is None
        True
    """
    op = ex.operator()
    if op != operator.pow:
        return None
    base, exp = ex.operands()
    try:
        c = QQ(base)
    except TypeError:
        return None
    if c <= 0:
        return None
    lin = _linear_coeff(exp, n)
    if lin is None:
        return None
    return c, lin[0], lin[1]


def _linear_symbolic(ex, n):
    r"""
    Return ``(a, b)`` if ``ex == a*n + b`` with ``a`` and ``b`` independent
    of ``n``, else ``None``.

    INPUT:

    - ``ex`` -- a symbolic expression

    - ``n`` -- the variable

    EXAMPLES::

        sage: from sage.symbolic.gamma_products import _linear_symbolic
        sage: k, x = var('k x')
        sage: _linear_symbolic(x*k + 2, k)
        (x, 2)
        sage: _linear_symbolic(k^2, k) is None
        True
        sage: _linear_symbolic(sin(x), k)
        (0, sin(x))
    """
    if not ex.is_polynomial(n):
        return None
    coeffs = {d: c for c, d in ex.coefficients(n)}
    if not set(coeffs) <= {0, 1}:
        return None
    return coeffs.get(1, SR(0)), coeffs.get(0, SR(0))


def _strip_exponentials(E, n):
    r"""
    Replace each atom ``c^(a*n+b)`` with positive rational base by `1`.

    INPUT:

    - ``E`` -- a symbolic expression

    - ``n`` -- the variable

    OUTPUT: a pair ``(expression, counter)`` where ``counter`` records the
    constant contribution of the stripped atoms as prime exponents.

    Only atoms which are factors of ``E`` are stripped: atoms occurring in a
    non-multiplicative context, for example under an addition, are kept, so
    that the conversion of ``E`` to a rational function fails, as it should.

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
        sage: _strip_exponentials(sqrt(4^(k + 1))/(k + 1)/sqrt(4^k), k)
        (1/(k + 1), Counter({2: 1}))

        Atoms which are not factors of the expression are kept::

            sage: _strip_exponentials(2^k + 1, k)
            (2^k + 1, Counter())

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
            # an atom, or a rational power of an atom, is again an atom
            atom = _rational_atom(ex, n)
            if atom is None:
                try:
                    r = QQ(exp)
                except TypeError:
                    r = None
                if r is not None and r != 0:
                    inner = _rational_atom(base, n)
                    if inner is not None:
                        c0, a, b = inner
                        atom = (c0, a * r, b * r)
            if atom is not None:
                c0, a, b = atom
                for p, kp in factor(c0):
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
        if op in (mul_vararg, operator.mul):
            return op(*[walk(a, mult) for a in args])
        # atoms in a non-multiplicative context cannot be factored out, so
        # they are kept and the conversion to a rational function must fail
        return ex

    new = walk(E, 1)
    if any(slopes.values()):
        growth = prod(p ** s for p, s in slopes.items())
        raise ValueError("expression is not hypergeometric "
                         "(exponential part grows like (%s)^n)" % growth)
    return new, const


def _split_envelope(E, n):
    r"""
    Split ``E`` into a part free of exponentials and its exponential envelope.

    INPUT:

    - ``E`` -- a symbolic expression

    - ``n`` -- the variable

    OUTPUT: a pair ``(h, z)`` of symbolic expressions such that
    ``E == h*z^n``, where every atom ``c^(a*n + b)`` of ``E``, with ``c``
    independent of ``n`` and ``a*n + b`` linear in ``n``, is replaced by
    the constant ``c^b``, and ``z`` is the product of the ``c^a``.

    EXAMPLES::

        sage: from sage.symbolic.gamma_products import _split_envelope
        sage: k, x = var('k x')
        sage: _split_envelope(2^(3*k + 1)*x^k*(k + 1), k)
        (2*k + 2, 8*x)
        sage: _split_envelope(e^(2*k)*sin(x), k)
        (sin(x), e^2)
        sage: _split_envelope((x^k + 1)^2, k)
        ((x^k + 1)^2, 1)
    """
    z = SR(1)

    def walk(ex, mult):
        nonlocal z
        op = ex.operator()
        if op is None:
            return ex
        args = ex.operands()
        if op == EXP_OP:
            lin = _linear_symbolic(args[0], n)
            if lin is None:
                return ex
            a, b = lin
            z *= E_CONST ** (a * mult)
            return E_CONST ** (b * mult)
        if op == operator.pow:
            base, expo = args
            if n not in base.variables():
                lin = _linear_symbolic(expo, n)
                if lin is not None:
                    a, b = lin
                    z *= base ** (a * mult)
                    return base ** (b * mult)
            # integer exponent: descend, tracking the multiplicity
            try:
                k = ZZ(expo)
            except (TypeError, ValueError):
                k = None
            if k is not None:
                return walk(base, mult * k) ** expo
            return ex
        if op in (mul_vararg, operator.mul):
            return op(*[walk(arg, mult) for arg in args])
        # atoms in a non-multiplicative context cannot be factored out
        return ex

    return walk(E, 1), z


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
        sage: _factor_ratio(3/2, k)
        (3/2, [])
        sage: _factor_ratio(sqrt(2), k)
        Traceback (most recent call last):
        ...
        ValueError: expression is not hypergeometric (irrational ratio unit sqrt(2))
    """
    def unit_from_primes(sign):
        # assemble the unit from the prime exponents, rejecting irrational
        # contributions such as sqrt(2)
        for v in unitp.values():
            if v not in ZZ:
                raise ValueError("expression is not hypergeometric "
                                 "(irrational ratio unit)")
        return sign * prod(p ** v for p, v in unitp.items() if v)

    R = SR(R)
    # a constant ratio contributes only to the unit; the conversion of a
    # constant to a fraction field is not supported
    if not R.has(n):
        unitp = Counter(constp)
        try:
            c = QQ(R)
        except (TypeError, ValueError):
            raise ValueError("expression is not hypergeometric "
                             "(irrational ratio unit %s)" % R) from None
        if not c:
            raise ValueError("expression is not hypergeometric (vanishing ratio)")
        F = factor(c)
        for p, kp in F:
            unitp[p] += kp
        return unit_from_primes(F.unit()), []
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
    return unit_from_primes(sign), factors


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
    try:
        R = (Eg.subs({n: n + 1}) / Eg).gamma_normalize()
    except ZeroDivisionError as exc:
        # the term is zero or has a pole, so that the ratio cannot be formed
        raise ValueError("the term is zero or undefined at var=0 or var=1"
                         ) from exc
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

# bound on the number of states explored in the search for a complete
# elimination of the fractional shifts; when it is exceeded, a normal form is
# committed instead of discarding the progress
_SEARCH_BUDGET = 2000


class _SearchBudget(Exception):
    r"""
    Raised when the search for a complete elimination exceeds its bound.
    """


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
        self._nodes = 0            # states explored by _search

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

        Constant classes whose members telescope across a pole are
        rejected, since the term vanishes from some point on::

            sage: s = _GammaSolver(k); s.feed([(1, -1, 1), (1, 11, -1)])
            sage: s.normalize()
            Traceback (most recent call last):
            ...
            ValueError: gamma argument reaches a pole at 0 (term undefined at some n >= 0)
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
                    if not x0 + j:
                        # the telescoping of the constant gamma ratio
                        # crosses a pole: the term vanishes from some
                        # point on, or has a pole
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

    def _circle_moves(self):
        r"""
        Return all applicable Gauss circles of the current state.

        A circle is a tuple ``(m, b, e, consumed, introduced)``, meaning that
        ``prod_{j=0}^{m-1} gamma(a*n + b + j/m)`` is replaced by
        ``gamma(m*a*n + m*b)``, consuming the entries ``consumed`` with total
        exponent ``e`` and introducing the factorials ``factorial(a*n + y - 1)``
        for the integer members ``y``.

        The circles are ordered by their smallest member, so that each is
        generated exactly once.

        EXAMPLES::

            sage: from sage.symbolic.gamma_products import _GammaSolver
            sage: k = var('k')
            sage: s = _GammaSolver(k); s.feed([(3, 1, 1), (3, 2, 1)])
            sage: s._circle_moves()
            [(3, 1/3, -1, [(0, 1/3), (0, 2/3)], [1]),
             (3, 1/3, 1, [(1, 1/3), (1, 2/3)], [1])]
        """
        state = self.state
        if not state:
            return []
        q = lcm(x.denominator() for _, x in state)
        # represent the shifts as integer multiples of 1/q, so that the
        # membership test in the circle enumeration is a set lookup
        ent = {}
        for (a, x), e in state.items():
            if e:
                ent.setdefault(a, {})[x * q] = e
        out = []
        for a in sorted(ent):
            ae = ent[a]
            for x0u in sorted(ae):
                # a circle of modulus m requires m members at this slope,
                # of which the member 1 needs no state entry
                for m in divisors(q):
                    if m < 2 or m > len(ae) + 1:
                        continue
                    step = q // m
                    r = x0u % step
                    # the members are the residue class of the shift modulo
                    # 1/m, as shifts in the interval (0, 1); the class
                    # representative 0 is replaced by 1, which contributes a
                    # factorial instead of a gamma
                    window = ([r + t*step if r else (t + 1)*step
                               for t in range(m)])
                    consumed, introduced, ok = [], [], True
                    for u in window:
                        if u == q:
                            introduced.append(ZZ(1))
                        elif u in ae:
                            consumed.append((a, QQ((u, q))))
                        else:
                            ok = False
                            break
                    if not ok or x0u != window[0]:
                        continue      # generated from its smallest member
                    b = QQ((r, q)) if r else QQ((1, m))
                    mb = m * b
                    if mb in ZZ and mb < 1:
                        continue      # would land on a pole
                    s0 = ae[x0u]
                    if all(ae[u] * s0 > 0 for u in window if u != q):
                        # the admissible multiplicities have the common sign
                        # of the consumed entries, and magnitudes up to their
                        # minimum
                        sgn = 1 if s0 > 0 else -1
                        mag = min(abs(ae[u]) for u in window if u != q)
                        out.extend((m, b, e, consumed, introduced)
                                   for e in range(sgn * mag, 0, -sgn))
        return out

    def apply_circle(self, mv):
        r"""
        Apply the Gauss circle ``mv``, updating the corrections.

        INPUT:

        - ``mv`` -- a circle, as returned by :meth:`_circle_moves`

        EXAMPLES::

            sage: from sage.symbolic.gamma_products import _GammaSolver
            sage: k = var('k')
            sage: s = _GammaSolver(k); s.feed([(3, 1, 1), (3, 2, 1)])
            sage: mv = [mv for mv in s._circle_moves() if min(mv[3]) == (1, 1/3)][0]
            sage: s.apply_circle(mv)
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
            self.fact[(m * a, ZZ(mb) - 1)] += e
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

    def _snapshot(self):
        r"""
        Return a deep copy of the solver state, for later :meth:`_restore`.

        EXAMPLES::

            sage: from sage.symbolic.gamma_products import _GammaSolver
            sage: k = var('k')
            sage: s = _GammaSolver(k); s.feed([(3, 1, 1)])
            sage: snap = s._snapshot()
            sage: s.state[(1, 1/3)] += 1; s.d *= 2
            sage: s._restore(snap)
            sage: sorted(s.state.items()), s.d
            ([((0, 1/3), -1), ((1, 1/3), 1)], 3)
        """
        return (Counter(self.state), Counter(self.fact), self.d, self.unit,
                Counter(self.primes), self.pi_exp, list(self.polys))

    def _restore(self, snap):
        r"""
        Restore a state returned by :meth:`_snapshot`.

        EXAMPLES::

            sage: from sage.symbolic.gamma_products import _GammaSolver
            sage: k = var('k')
            sage: s = _GammaSolver(k)
            sage: s._restore(s._snapshot())
            sage: sorted(s.state.items())
            []
        """
        (self.state, self.fact, self.d, self.unit, self.primes,
         self.pi_exp, polys) = snap
        self.polys = polys

    def _search(self, budget):
        r"""
        Depth-first search for a sequence of Gauss circles emptying the state.

        Raise :class:`_SearchBudget` when more than ``budget`` states have
        been explored.

        EXAMPLES::

            sage: from sage.symbolic.gamma_products import _GammaSolver, _SearchBudget
            sage: k = var('k')
            sage: s = _GammaSolver(k); s.feed([(3, 1, 1), (3, 2, 1)])
            sage: s.normalize(); s._search(10)
            True
            sage: sorted(s.state.items())
            []
            sage: s = _GammaSolver(k); s.feed([(3, 1, 1), (3, 2, 1)])
            sage: s.normalize(); s._search(1)
            Traceback (most recent call last):
            ...
            sage.symbolic.gamma_products._SearchBudget
        """
        self._nodes += 1
        if self._nodes > budget:
            raise _SearchBudget
        if not self.state:
            return True
        key = tuple(sorted(self.state.items()))
        if key in self.memo:
            return False
        for mv in self._circle_moves():
            saved = self._snapshot()
            self.apply_circle(mv)
            if self._search(budget):
                return True
            self._restore(saved)
        self.memo.add(key)
        return False

    def _descend(self):
        r"""
        Apply Gauss circles until the state is in normal form.

        Every application strictly decreases the total number of fractional
        gamma factors, counted with multiplicity, so this terminates; the
        final state admits no circle.

        EXAMPLES::

            sage: from sage.symbolic.gamma_products import _GammaSolver
            sage: k = var('k')
            sage: s = _GammaSolver(k); s.feed([(6, 1, 1), (6, 5, 1)])
            sage: s._descend(), sorted(s.state.items())
            (None, [((0, 1/6), -1), ((0, 5/6), -1), ((1, 1/6), 1), ((1, 5/6), 1)])
        """
        while True:
            moves = self._circle_moves()
            if not moves:
                return
            self.apply_circle(moves[0])

    def solve(self, budget=_SEARCH_BUDGET):
        r"""
        Eliminate all fractional gamma shifts, as far as possible.

        Return ``True`` if the state could be emptied completely.  Otherwise
        a normal form is installed: a state to which no Gauss circle applies,
        so that no progress is discarded.

        INPUT:

        - ``budget`` -- the number of states explored in the search for a
          complete elimination; when it is exceeded, the normal form is
          committed without knowing whether the state can be emptied

        EXAMPLES::

            sage: from sage.symbolic.gamma_products import _GammaSolver
            sage: k = var('k')
            sage: s = _GammaSolver(k); s.feed([(3, 1, 1), (3, 2, 1)])
            sage: s.solve(), sorted(s.state.items())
            (True, [])

            A state which admits no circle is a normal form; the residual
            gammas are irreducible::

            sage: s = _GammaSolver(k); s.feed([(6, 1, 1), (6, 5, 1)])
            sage: s.solve(), sorted(s.state.items())
            (False, [((0, 1/6), -1), ((0, 5/6), -1), ((1, 1/6), 1), ((1, 5/6), 1)])

            With a vanishing search budget, the normal form is found by
            greedy descent alone; here it empties the state::

            sage: s = _GammaSolver(k); s.feed([(3, 1, 1), (3, 2, 1)])
            sage: s.solve(budget=0), sorted(s.state.items())
            (True, [])
        """
        self.normalize()
        if not self.state:
            return True
        saved = self._snapshot()
        self._nodes = 0
        try:
            if self._search(budget):
                return True
        except _SearchBudget:
            pass
        # either the search failed after full backtracking, or it was cut
        # short: restore the initial state in both cases and commit the
        # progress of the greedy descent below
        self._restore(saved)
        self._descend()
        return not self.state

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
    gamma factors, if there are any, are in irreducible position: all shifts
    lie in the interval `[0, 1)`, shifts with the same slope `a_j` are
    distinct modulo the integers, and no application of the Gauss
    multiplication formula reduces the product further.

    A :exc:`ValueError` is raised if ``expression`` is not a hypergeometric
    term, does not split into linear factors over ``QQ``, or is zero or
    undefined at ``var = 0``.

    The normal form is not canonical: in particular, constants may be
    represented by gamma values at rational points or by evaluated powers
    of ``pi``, and simplifying an already simplified expression again may
    change the representation.

    EXAMPLES::

        sage: from sage.symbolic.gamma_products import simplify_gamma
        sage: n = var('n')
        sage: simplify_gamma(binomial(3*n, n), n)
        factorial(3*n)/(factorial(2*n)*factorial(n))
        sage: simplify_gamma(gamma(n + 1/2), n)
        sqrt(pi)*(1/4)^n*factorial(2*n)/factorial(n)
        sage: simplify_gamma(factorial(n)*gamma(n + 1/2)/sqrt(pi), n)
        (1/4)^n*factorial(2*n)

    Purely exponential terms are accepted, in particular those whose ratio
    is constant::

        sage: simplify_gamma(2^n, n)
        2^n
        sage: simplify_gamma((2/3)^n, n)
        (2/3)^n

    Exponential atoms which are not factors of the term are kept, so that
    the conversion of the ratio to a rational function fails::

        sage: simplify_gamma(2^n + 1, n)
        Traceback (most recent call last):
        ...
        ValueError: expression is not hypergeometric (ratio (2*2^n + 1)/(2^n + 1))

    A lone gamma with fractional shift cannot be reduced::

        sage: simplify_gamma(gamma(n + 2/3)/gamma(2/3), n)
        gamma(n + 2/3)/gamma(2/3)

    but two of them combine, by the Gauss multiplication formula::

        sage: simplify_gamma(gamma(n + 1/3)*gamma(n + 2/3), n)
        2/3*sqrt(3)*pi*(1/27)^n*factorial(3*n)/factorial(n)

        When a complete reduction is impossible, the progress is committed
        and the residual gamma factors are in irreducible position; here the
        constants gamma(1/6) and gamma(5/6) admit no multiplication formula::

            sage: simplify_gamma(gamma(2*n + 1/3)*gamma(2*n + 2/3), n)
            4/3*sqrt(3)*pi^2*(1/729)^n*factorial(6*n)/(factorial(2*n)*gamma(5/6)*gamma(1/6))

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
        sage: simplify_gamma(gamma(n - 1)/gamma(-1), n)
        Traceback (most recent call last):
        ...
        ValueError: the term is zero or undefined at var=0 or var=1
        sage: simplify_gamma(rising_factorial(-1, n), n)
        Traceback (most recent call last):
        ...
        ValueError: the term is zero or undefined at var=0 or var=1
    """
    n = SR(var) if not isinstance(var, str) else SR.var(var)
    if not (isinstance(n, Expression) and n.is_symbol()):
        raise TypeError("the variable must be a symbol")
    E = SR(expression)


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
        solver.state[(0, x)] += e

    solver.solve()
    res = solver.expression(rest)
    return res


def _factor_count(m):
    r"""
    Validate the number of factors of a finite product.

    INPUT:

    - ``m`` -- the expression ``b - a + 1`` for the bounds of the product

    OUTPUT: the pair ``(count, lin)``.  If the number of factors is a
    concrete constant, ``count`` is a nonnegative integer and ``lin`` is
    ``None``.  Otherwise ``count`` is ``None`` and ``lin`` is the pair
    ``((n, q, r0), n0)``, meaning that the number of factors equals
    ``q*n0 + r0`` with a positive integer slope ``q`` and a nonnegative
    integer ``r0``; the count variable ``n`` is a valid generator name,
    renamed if necessary.

    A :exc:`ValueError` is raised if the number of factors is neither of
    these.

    EXAMPLES::

        sage: from sage.symbolic.gamma_products import _factor_count
        sage: n = var('n')
        sage: _factor_count(3*n + 2)
        (None, ((n, 3, 2), n))
        sage: _factor_count(SR(5))
        (5, None)
        sage: _factor_count(n^2)
        Traceback (most recent call last):
        ...
        ValueError: the number of factors n^2 is not linear in a single variable with rational coefficients
        sage: _factor_count(1/2)
        Traceback (most recent call last):
        ...
        ValueError: the number of factors 1/2 is not a nonnegative integer
    """
    m = SR(m)
    m0 = m
    mvars = m.variables()
    if not mvars:
        try:
            count = QQ(m)
        except (TypeError, ValueError):
            raise ValueError("the number of factors %s is not rational" % m)
        if count not in ZZ or count < 0:
            raise ValueError("the number of factors %s is not a nonnegative "
                             "integer" % count)
        return ZZ(count), None
    if len(mvars) > 1:
        raise ValueError("the number of factors %s is not linear in a single "
                         "variable" % m)
    n = mvars[0]
    lin = _linear_coeff(m, n)
    if lin is None:
        # report the failure against the original expression, not the
        # renamed one
        raise ValueError("the number of factors %s is not linear in a single "
                         "variable with rational coefficients" % m0)
    q, r0 = lin
    if q not in ZZ or q < 1:
        raise ValueError("the number of factors grows by %s per step" % q)
    if r0 not in ZZ or r0 < 0:
        raise ValueError("the number of factors is %s at the anchor" % r0)
    return None, ((n, ZZ(q), ZZ(r0)), mvars[0])


def _rational_range_product(f, w, a, count):
    r"""
    Evaluate ``prod_{j=0}^{count-1} f(a + j)`` exactly over the rationals.

    A :exc:`TypeError` or :exc:`ValueError` is raised if some factor is
    irrational or has a pole.

    EXAMPLES::

        sage: from sage.symbolic.gamma_products import _rational_range_product
        sage: k = var('k')
        sage: _rational_range_product(k/(k + 1), k, 1, 4)
        1/5
        sage: _rational_range_product(sqrt(2)*k, k, 0, 3)
        Traceback (most recent call last):
        ...
        TypeError: unable to convert sqrt(2) to a rational
    """
    return prod(QQ(f.subs({w: a + j})) for j in range(count))


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

        Products which vanish from some point on, or have a pole in the
        range, are rejected::

            sage: factorial_product((k - 7)/(k + 5), k, 2, n + 5)
            Traceback (most recent call last):
            ...
            ValueError: the product vanishes from some point on: the body has a zero at 7
            sage: factorial_product((k - 4)/(k - 3), k, 0, n + 1)
            Traceback (most recent call last):
            ...
            ValueError: the product is undefined from some point on: the body has a pole at 3
            sage: factorial_product(1/(k - 3), k, 0, n + 3)
            Traceback (most recent call last):
            ...
            ValueError: cannot evaluate the product at the anchor: the product has a pole or an irrational factor

        Infinite bounds are rejected::

            sage: factorial_product(k, k, 1, oo)
            Traceback (most recent call last):
            ...
            ValueError: product bounds must be finite
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

    w = v
    f0 = f

    # check that f is a rational function of the product variable with
    # rational coefficients; any failure of the conversion counts as a
    # rejection
    Pv = PolynomialRing(QQ, w)
    try:
        fr = Pv.fraction_field()(f)
    except Exception as exc:
        raise ValueError("%s is not a rational function of %s with rational "
                         "coefficients" % (f0, v)) from exc

    count, lin = _factor_count(bS - aS + 1)
    if lin is None:
        try:
            aq = QQ(a)
        except (TypeError, ValueError):
            raise ValueError("the lower bound %s is not rational" % a)
        try:
            res = _rational_range_product(f, w, aq, count)
        except (TypeError, ValueError):
            raise ValueError("%s cannot be evaluated exactly over the given "
                             "range" % f)
        if not res:
            raise ValueError("the product vanishes")
        return SR(res)

    (n, q, r0), n0 = lin
    try:
        aq = QQ(a)
    except (TypeError, ValueError):
        raise ValueError("the lower bound %s is not rational" % a)

    # reject bodies with zeros or poles in the shifting part of the
    # range, which grows without bound: the product vanishes or is
    # undefined from some point on and has no normal form.  the anchor
    # part of the range is checked by direct evaluation below.  the
    # product variable runs over the lattice a + ZZ, so only roots in
    # that lattice can be hit
    for root, _ in fr.denominator().roots():
        if (root - aq) in ZZ and root > aq + r0 - 1:
            raise ValueError("the product is undefined from some point on: "
                             "the body has a pole at %s" % root)
    num = fr.numerator()
    if num.is_zero():
        raise ValueError("the product vanishes")
    for root, _ in num.roots():
        if (root - aq) in ZZ and root > aq + r0 - 1:
            raise ValueError("the product vanishes from some point on: "
                             "the body has a zero at %s" % root)

    # the term ratio is the product of the last q factors,
    # E(n+1)/E(n) = prod_{j=0}^{q-1} f(a + m(n) + j)
    R = prod(f.subs({w: aq + q*n + r0 + j}) for j in range(q))
    unit, factors = _factor_ratio(R, n)

    solver = _GammaSolver(n)
    solver.feed(factors)
    solver.d *= unit

    # anchor: E(0) is the product of the first r0 factors
    try:
        E0 = _rational_range_product(f, w, aq, r0)
    except (TypeError, ValueError):
        raise ValueError("cannot evaluate the product at the anchor: "
                         "the product has a pole or an irrational factor")
    if not E0:
        raise ValueError("the product vanishes at the anchor")
    solver.unit *= E0

    solver.solve()
    res = solver.expression(1)
    return res


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
      ``gamma(x+n+1)/gamma(x)``; over a concrete range the finite product is
      evaluated directly, since the gamma quotient has removable poles at the
      zeros of the factors, so that its size grows linearly with the number
      of factors,

    - factors independent of the product variable, which are evaluated as
      ``c^(b - a + 1)``,

    - a residual part, which is returned as an unevaluated product in the
      original variable and bounds.

    The number of factors ``b - a + 1`` must be a nonnegative integer, or an
    integer-linear expression ``q*n + r`` in a single variable ``n`` with
    ``q`` a positive integer and ``r`` a nonnegative integer; a
    :exc:`ValueError` is raised otherwise, and also if nothing can be
    evaluated.  In both cases the caller may want to fall back to another
    method.

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

        Parameter-linear factors over a concrete range are evaluated
        directly, so that removable gamma poles do not become evaluation
        errors::

            sage: exact_product(x + k, k, 0, 1)
            (x + 1)*x
            sage: exact_product(x + k, k, 0, 1).subs(x = -2)
            2

        For a symbolic number of factors the gamma quotient is returned,
        which then evaluates like any other gamma quotient, in particular
        failing at the removable poles::

            sage: exact_product(x + k, k, 0, n)
            gamma(n + x + 1)/gamma(x)
            sage: try:
            ....:     exact_product(x + k, k, 0, n).subs(x = -2, n = 1)
            ....: except RuntimeError as err:
            ....:     print('RuntimeError: %s' % err)
            RuntimeError: indeterminate expression: 0 * infinity encountered.

        The empty product is 1, whatever the body is::

            sage: exact_product(cos(k), k, 0, -1)
            1

        The number of factors must be a nonnegative integer, or an
        integer-linear expression in a single variable::

            sage: exact_product(x + k, k, 0, 1/2)
            Traceback (most recent call last):
            ...
            ValueError: the number of factors 3/2 is not a nonnegative integer
            sage: exact_product(x, k, 0, n^2)
            Traceback (most recent call last):
            ...
            ValueError: the number of factors n^2 + 1 is not linear in a single variable with rational coefficients

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

    # the formulas below are valid only if the number of factors is a
    # nonnegative integer, or an integer-linear expression in a single
    # variable with a positive slope; otherwise fall back to the caller
    count, lin = _factor_count(bS - aS + 1)
    if count == 0:
        # the empty product
        return SR(1)

    w = v

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
        if count is not None:
            # over a concrete range, evaluate the finite product directly:
            # the gamma quotient has removable poles at the zeros of the
            # factors
            parts.append(A**(e*count) * prod(c + aS + j for j in range(count))**e)
        else:
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
        parts.append(symbolic_product(h, v, a, b))

    res = prod(parts)
    return res


def hypergeometric_sum(f, var, a, b):
    r"""
    Evaluate the infinite sum of a hypergeometric term in hypergeometric
    notation.

    The summand must be hypergeometric in the summation variable: its term
    ratio, reduced by applying the recurrence of the gamma function, must
    factor into a constant times integer powers of linear polynomials.  The
    summation variable must not occur in the coefficients of these
    polynomials, but arbitrary other symbols may occur, both in the term and
    in its parameters.

    INPUT:

    - ``f`` -- a symbolic expression, the term of the sum

    - ``var`` -- the summation variable

    - ``a`` -- the lower bound, a rational number

    - ``b`` -- the upper bound, ``oo`` or ``-oo``

    As in Sage's symbolic sums, the summation variable runs over ``a``,
    ``a + 1``, ``a + 2``, ..., also for non-integral lower bounds.  Since
    a hypergeometric term is a product of gamma functions, it is defined
    for non-integral arguments, unless one of the gamma functions has a
    pole there.

    OUTPUT: the sum `\sum_{k=a}^{b} f(k)`, a constant times
    ``hypergeometric(A, B, z)``.  Identical parameters in ``A`` and ``B``
    are eliminated, in particular the parameter `1` accounting for the
    implicit `1/k!` of the hypergeometric series.

    See also
    :meth:`~sage.symbolic.expression.Expression.to_hypergeometric`, which
    rewrites every unevaluated symbolic sum in an expression in this form.

    A :exc:`ValueError` is raised if the term is not hypergeometric, if the
    term is not defined at some point of the range, or if the bounds are not
    of the required form.  In particular, a term which is zero at some point
    of the range, but not identically zero from then on, is not
    hypergeometric, and a terminating series must terminate via a negative
    integer numerator parameter.

    EXAMPLES::

        sage: from sage.symbolic.gamma_products import hypergeometric_sum
        sage: k, x, a = var('k x a')
        sage: hypergeometric_sum(rising_factorial(a, k)*x^k/factorial(k), k, 0, oo)
        hypergeometric((a,), (), x)

        The binomial series is a hypergeometric series with a single
        numerator parameter, which cancels the implicit parameter 1::

            sage: hypergeometric_sum(binomial(a, k)*x^k, k, 0, oo)
            hypergeometric((-a,), (), -x)

        The Bessel function of the first kind::

            sage: hypergeometric_sum((-1)^k*x^(2*k)/(4^k*factorial(k)^2), k, 0, oo)
            hypergeometric((), (1,), -1/4*x^2)
            sage: hypergeometric_sum((-1)^k*x^(2*k)/(4^k*factorial(k)^2), k, 0, oo).simplify_hypergeometric()
            bessel_J(0, x)

        Nonzero lower bounds shift the parameters, leading zeroes of the
        term are skipped::

            sage: hypergeometric_sum(k^2*x^k/factorial(k), k, 1, oo)
            x*hypergeometric((2,), (1,), x)

            sage: hypergeometric_sum(k^3*x^k/factorial(k), k, 0, oo)
            x*hypergeometric((2, 2), (1, 1), x)

        Gamma values with a symbolic shift contribute a constant prefactor::

            sage: hypergeometric_sum(gamma(k + a)/(gamma(a)*factorial(k))*x^k, k, 0, oo)
            hypergeometric((a,), (), x)

        Non-integral lower bounds make the summation variable run over
        ``a + Z``, exactly as in Sage's symbolic sums::

            sage: hypergeometric_sum(x^k, k, 1/2, oo).simplify_hypergeometric()
            -sqrt(x)/(x - 1)
            sage: sum(x^k, k, 1/2, 5/2)
            (x^2 + x + 1)*sqrt(x)

        Symbolic parameters may also occur in the denominator::

            sage: hypergeometric_sum(factorial(k)*x^k/rising_factorial(a, k), k, 0, oo)
            hypergeometric((1, 1), (a,), x)

        TESTS::

        The term ratio is verified numerically against direct summation::

            sage: f = binomial(2*k, k)^2*x^k/factorial(4*k)
            sage: r = hypergeometric_sum(f, k, 0, oo)
            sage: bool(abs(r.subs(x=1/8).n(50)
            ....:           - sum(f.subs({x: 1/8, k: m}).n(50) for m in range(30))) < 1e-12)
            True

        Sums to negative infinity negate the summation variable::

            sage: hypergeometric_sum(2^k, k, -2, -oo).simplify_hypergeometric()
            1/2

        Leading zeroes of the term may hide poles of the gamma function::

            sage: hypergeometric_sum(1/gamma(k - 2), k, 0, oo).simplify_hypergeometric()
            e

        A term which vanishes on the whole summation lattice sums to
        zero::

            sage: hypergeometric_sum(1/gamma(-k), k, 0, oo)
            0
            sage: hypergeometric_sum(sin(pi*k), k, 0, oo)
            0

        Terms whose term ratio does not factor into linear factors are
        rejected::

            sage: hypergeometric_sum(2^(k^2), k, 0, oo)
            Traceback (most recent call last):
            ...
            ValueError: the term 2^(k^2) is not hypergeometric in k
            sage: hypergeometric_sum(sin(k)*x^k, k, 0, oo)
            Traceback (most recent call last):
            ...
            ValueError: the term x^k*sin(k) is not hypergeometric in k

        So are terms with a pole in the range::

            sage: hypergeometric_sum(1/(k - 2), k, 0, oo)
            Traceback (most recent call last):
            ...
            ValueError: the term 1/(k - 2) is not defined in the range

        A term whose evaluation at a lattice point is indeterminate, for
        instance because a zero cancels a pole of a gamma factor, is
        rejected as well::

            sage: hypergeometric_sum(gamma(-k - 2)/gamma(-k - 3), k, 0, oo)
            Traceback (most recent call last):
            ...
            ValueError: the term cannot be evaluated at the lower bound 0

        A term may be defined everywhere, but still fail to be
        hypergeometric, because its term ratio has a pole inside the range::

            sage: hypergeometric_sum((k - 3)*x^k/factorial(k), k, 0, oo)
            Traceback (most recent call last):
            ...
            ValueError: the term (k - 3)*x^k/factorial(k) is not hypergeometric in k

        The bounds must have the required form::

            sage: n = var('n')
            sage: hypergeometric_sum(x^k, k, 0, n)
            Traceback (most recent call last):
            ...
            ValueError: the upper bound n is not infinite
            sage: hypergeometric_sum(x^k, k, x, oo)
            Traceback (most recent call last):
            ...
            ValueError: the lower bound x is not rational
    """
    if not (isinstance(var, Expression) and var.is_symbol()):
        raise TypeError("the variable must be a symbol")
    try:
        aq = QQ(a)
    except (TypeError, ValueError):
        raise ValueError("the lower bound %s is not rational" % a) from None
    if b == minus_infinity:
        return hypergeometric_sum(f.subs({var: -var}), var, -aq, infinity)
    if not b == infinity:
        # note that the comparison with a symbolic bound must not evaluate
        # to a relational
        raise ValueError("the upper bound %s is not infinite" % b)

    j = SR.temp_var()
    g = f.subs({var: aq + j})
    if g.is_trivial_zero():
        return SR(0)

    # the term is hypergeometric iff its term ratio, reduced by the
    # recurrence of the gamma function, factors into linear polynomials;
    # compute the ratio before skipping leading zeros, since its poles at
    # nonnegative integers are the only positions at which a vanishing
    # term can be followed by a nonzero one, and therefore bound the number
    # of leading zeros
    h0, _ = _split_envelope(g.to_gamma(), j)
    Rf0 = (h0.subs({j: j + 1}) / h0).gamma_normalize().factor()
    bound = 0
    for ex, side in ((Rf0.numerator(), 1), (Rf0.denominator(), -1)):
        for base, e in _mult_split(ex):
            lin = _linear_symbolic(base, j)
            if lin is None:
                raise ValueError("the term %s is not hypergeometric in %s"
                                 % (f, var))
            if side < 0:
                A, B = lin
                if A.is_trivial_zero():
                    # a constant factor has no zeros
                    continue
                try:
                    pole = ZZ(-B/A)
                except (TypeError, ValueError, ArithmeticError):
                    # the position of the pole depends on symbolic parameters
                    continue
                if pole >= 0:
                    bound = max(bound, pole + 1)

    try:
        skipped = 0
        while g.subs({j: 0}).is_trivial_zero():
            if skipped == bound:
                # beyond the last pole of the term ratio, a vanishing term
                # is followed by a vanishing one: every term on the
                # summation lattice vanishes
                return SR(0)
            g = g.subs({j: j + 1})
            skipped += 1
        g0 = g.subs({j: 0})
        undefined = g0.is_infinity() or not (g0 - g0).is_trivial_zero()
    except (ValueError, ArithmeticError, RuntimeError) as exc:
        raise ValueError("the term cannot be evaluated at the "
                        "lower bound %s" % a) from exc
    if undefined:
        raise ValueError("the term is not defined at the lower bound %s" % a)

    # split off the exponential envelope:  g = h * z^j
    h, z = _split_envelope(g.to_gamma(), j)
    R = (h.subs({j: j + 1}) / h).gamma_normalize()
    Rf = R.factor()

    ap, bq = [], []
    for ex, side in ((Rf.numerator(), 1), (Rf.denominator(), -1)):
        for base, e in _mult_split(ex):
            lin = _linear_symbolic(base, j)
            if lin is None:
                raise ValueError("the term %s is not hypergeometric in %s"
                                 % (f, var))
            A, B = lin
            if A.is_trivial_zero():
                z *= base ** (e * side)
            else:
                z *= A ** (e * side)
                (ap if side > 0 else bq).extend([B/A] * e)

    # the implicit 1/j! of the hypergeometric series, and pairwise
    # cancellation of identical parameters
    ap.append(SR(1))
    for q in list(bq):
        if q in ap:
            ap.remove(q)
            bq.remove(q)

    # a nonpositive integer denominator parameter is a pole of the term
    # ratio inside the range: either the term is undefined there, or the
    # term vanishes and the recurrence of the ratio breaks
    for q in bq:
        if q in ZZ and q <= 0:
            m = ZZ(-q)
            try:
                values = (g.subs({j: m}), g.subs({j: m + 1}))
                undefined = any(v.is_infinity() or not (v - v).is_trivial_zero()
                               for v in values)
            except (ValueError, ArithmeticError, TypeError):
                undefined = True
            if undefined:
                raise ValueError("the term %s is not defined in the range" % f)
            raise ValueError("the term %s is not hypergeometric in %s"
                             % (f, var))

    from sage.functions.hypergeometric import hypergeometric
    return g0 * hypergeometric(ap, bq, z)


def _to_hypergeometric(ex):
    r"""
    Rewrite every unevaluated symbolic sum in ``ex`` in hypergeometric
    notation.

    INPUT:

    - ``ex`` -- a symbolic expression

    Every subexpression which is an unevaluated symbolic sum is replaced by
    its hypergeometric form, see :func:`hypergeometric_sum`.  Sums which are
    not hypergeometric, or whose bounds are not of the required form, raise
    a :exc:`ValueError`.

    EXAMPLES::

        sage: from sage.symbolic.gamma_products import _to_hypergeometric
        sage: k, x = var('k x')
        sage: _to_hypergeometric(sum(k^2*x^k/factorial(k), k, 1, oo, hold=True))
        x*hypergeometric((2,), (1,), x)
        sage: _to_hypergeometric(sum(k^2*x^k/factorial(k), k, 1, oo, hold=True)
        ....:                      + 2*sum(k^3*x^k/factorial(k), k, 0, oo, hold=True))
        2*x*hypergeometric((2, 2), (1, 1), x) + x*hypergeometric((2,), (1,), x)

        Subexpressions which do not contain an unevaluated sum are kept, in
        particular unevaluated products are not evaluated::

            sage: from sage.calculus.calculus import symbolic_product
            sage: _to_hypergeometric(symbolic_product(k, k, 1, oo, hold=True))
            product(k, k, 1, +Infinity)
            sage: _to_hypergeometric(sin(sum(1/factorial(k), k, 0, oo, hold=True)))
            sin(hypergeometric((), (), 1))

        A sum which is not hypergeometric is rejected::

            sage: _to_hypergeometric(sum(sin(k)*x^k, k, 0, oo, hold=True))
            Traceback (most recent call last):
            ...
            ValueError: the term x^k*sin(k) is not hypergeometric in k
    """
    from sage.functions.other import Function_prod, Function_sum

    def walk(ex):
        op = ex.operator()
        if op is None:
            return ex
        args = ex.operands()
        new = [walk(a) for a in args]
        if isinstance(op, Function_sum):
            return hypergeometric_sum(*new)
        if all(n is a for n, a in zip(new, args)):
            return ex
        if isinstance(op, Function_prod):
            # keep unevaluated products unevaluated
            return op(*new, hold=True)
        return op(*new)

    return walk(ex)
