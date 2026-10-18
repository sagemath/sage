r"""
Small roots of multivariate modular and integer polynomials

This module implements Coppersmith's method for multivariate polynomials,
following the general strategy of Jochemsz and May [JM06]_.  Two variants
are provided: :func:`small_roots_modular` for polynomials with a root
modulo `N`, and :func:`small_roots_integer` for polynomials with a small
integer root (via Coron's reformulation).

The method is heuristic for more than one variable: LLL only guarantees
that the produced polynomials vanish over the integers under the
assumption that they are algebraically independent (Assumption 1 of
[JM06]_), which holds for all the examples included here.

AUTHORS:

- Xing Ji (starnotes-xj) (2026-10): clean-room implementation from
  [JM06]_

REFERENCES:

.. [JM06] \E. Jochemsz, A. May. A Strategy for Finding Roots of
   Multivariate Polynomials with New Applications in Attacking RSA
   Variants. ASIACRYPT 2006, LNCS 4284, pp. 267--282.
"""

# ****************************************************************************
#       Copyright (C) 2026 Xing Ji
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from functools import reduce

from sage.arith.misc import gcd
from sage.matrix.constructor import matrix
from sage.misc.mrange import prod
from sage.rings.finite_rings.integer_mod_ring import Zmod
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ


# ---------------------------------------------------------------------------
# monomial helpers (a monomial is a plain exponent tuple; index = variable)
# ---------------------------------------------------------------------------

def _terms(f):
    """[(exponent_tuple, coeff)] of ``f`` (handles univariate rings whose
    dict keys are plain ints)."""
    n = f.parent().ngens()
    return [((e,) if n == 1 else tuple(e), ZZ(c)) for e, c in f.dict().items()]


def _tadd(a, b):
    return tuple(x + y for x, y in zip(a, b))


def _tsub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def _tle(a, b):
    """componentwise a <= b"""
    return all(x <= y for x, y in zip(a, b))


def _closure(base, m):
    r"""Monomials of f^j for j = 0..m, given exponent tuples ``base`` of f.

    Returns ``levels`` with ``levels[j]`` the set of monomials of f^j.
    """
    n = len(base[0])
    levels = [{tuple([0] * n)}]
    for _ in range(m):
        nxt = set()
        for u in levels[-1]:
            for b in base:
                nxt.add(_tadd(u, b))
        levels.append(nxt)
    return levels


def _from_dict(P, d):
    out = P(0)
    for e, c in d.items():
        if c:
            out += c * prod(v**i for v, i in zip(P.gens(), e))
    return out


def _reduce_mod(f, N):
    P = f.parent()
    return _from_dict(P, {e: c % N for e, c in _terms(f) if c % N})


def _coeff(f, e):
    for et, c in _terms(f):
        if et == e:
            return c
    return ZZ(0)


def _leading_monomial(f):
    r"""A monomial `l` of `f` such that no other monomial of `f` is
    divisible by `l` (needed for a triangular basis; paper Section 2.1).
    Candidates are scanned from largest to smallest monomial."""
    exps = [e for e, _ in _terms(f)]
    for cand in sorted(exps, key=lambda e: (-sum(e), tuple(-x for x in e))):
        if not any(cand != e and _tle(cand, e) for e in exps):
            return cand
    raise ValueError("no valid leading monomial; try a change of variables")


def _scale_matrix(ring, shifts, bounds):
    """Build the (unscaled-row) lattice matrix.

    ``shifts``: list of (monomial_tuple, polynomial).  Columns are the union
    of the monomials of the shift polynomials, ordered by (degree, lex);
    column j is scaled by prod(bounds^e).  Returns (Matrix, cols, W).
    """
    n = ring.ngens()
    cols = sorted({u for u, _ in shifts}, key=lambda e: (-sum(e), e))
    W = [ZZ(prod(bounds[j]**e[j] for j in range(n))) for e in cols]
    idx = {e: i for i, e in enumerate(cols)}
    M = matrix(ZZ, len(shifts), len(cols))
    for r, (_, g) in enumerate(shifts):
        for e, c in _terms(g):
            j = idx[e]
            M[r, j] = c * W[j]
    return M, cols, W


def _solve_lattice(ring, shifts, bounds, hgbound, verify):
    """LLL + Howgrave-Graham filter + integer root extraction.

    ``shifts``: list of (monomial_tuple, polynomial) as in ``_scale_matrix``.
    ``hgbound``: a polynomial h with ||h(xX)|| < hgbound/sqrt(omega) holds
    over the integers at the root.
    """
    M, cols, W = _scale_matrix(ring, shifts, bounds)
    M = M.LLL()
    omega = len(cols)

    polys = []
    for r in range(M.nrows()):
        if M[r] == 0:
            continue
        if sum(ZZ(v)**2 for v in M[r]) * omega >= hgbound**2:
            continue
        polys.append(_from_dict(ring, {cols[j]: ZZ(M[r, j]) // W[j]
                                       for j in range(omega) if M[r, j]}))
    return _integer_points(polys, ring, bounds, verify)


def _integer_points(polys, ring, bounds, verify):
    """Integer points inside ``bounds`` satisfying ``verify``, that are
    common roots of ``polys``."""
    nz = [p for p in polys if p != 0]
    if not nz:
        return []

    def extract(polys_):
        QQring = ring.change_ring(QQ)
        gens = QQring.gens()
        pts = []
        try:
            J = QQring.ideal([QQring(p) for p in polys_])
            V = J.variety(ring=QQ)
            pts = [tuple(pt[v] for v in gens) for pt in V]
        except Exception:
            if ring.ngens() == 1:
                for p in polys_:
                    for r0, _ in p.roots(ring=ZZ):
                        pts.append((r0,))
            elif ring.ngens() == 2:
                xx, yy = ring.gens()
                for i in range(len(polys_)):
                    for j in range(i + 1, len(polys_)):
                        res = polys_[i].resultant(polys_[j], xx)
                        if res == 0:
                            continue
                        for y0, _ in res.univariate_polynomial() \
                                        .roots(ring=ZZ):
                            for h in polys_:
                                for x0, _ in h(y0).univariate_polynomial() \
                                                  .roots(ring=ZZ):
                                    pts.append((x0, y0))
            else:
                # Resultant cascade for n >= 3: eliminate the first
                # variable pairwise, eliminate the second variable between
                # two of the resultants, then back-substitute each integer
                # root of the final univariate and solve the remaining
                # bivariate system by resultant.  Only the shortest rows
                # are used: the cascade's cost is dominated by resultants
                # of resultants, which explode with degree.
                if len(polys_) > 4:
                    polys_ = polys_[:4]
                gens3 = ring.gens()
                for i in range(len(polys_)):
                    for j in range(i + 1, len(polys_)):
                        rij = polys_[i].resultant(polys_[j], gens3[0])
                        if rij == 0:
                            continue
                        for k2 in range(len(polys_)):
                            if k2 in (i, j):
                                continue
                            rk = polys_[k2].resultant(polys_[j], gens3[0])
                            if rk == 0:
                                continue
                            last = rij.resultant(rk, gens3[1])
                            if last == 0:
                                continue
                            lu = last.univariate_polynomial()
                            if lu.degree() < 1:
                                continue
                            for v0, _ in lu.roots(ring=ZZ):
                                sub = [h.subs({gens3[2]: ZZ(v0)})
                                       for h in polys_]
                                sub = [h for h in sub if h != 0]
                                if not sub:
                                    continue
                                for a in range(len(sub)):
                                    for b in range(a + 1, len(sub)):
                                        rs = sub[a].resultant(sub[b],
                                                              gens3[0])
                                        if rs == 0:
                                            continue
                                        rsub = rs.univariate_polynomial()
                                        if rsub.degree() < 1:
                                            continue
                                        for y0, _ in rsub.roots(ring=ZZ):
                                            for h in sub:
                                                hy = h.subs(
                                                    {gens3[1]: ZZ(y0)})
                                                if hy == 0:
                                                    pts.append((0, y0, v0))
                                                    continue
                                                if hy.degree() < 1:
                                                    continue
                                                for x0, _ in hy \
                                                        .univariate_polynomial() \
                                                        .roots(ring=ZZ):
                                                    pts.append((x0, y0, v0))
                raise ValueError("resultant cascade exhausted")
        return pts
        # (any extraction strategy may fail for heuristic reasons;
        # the caller falls back to the next strategy)

    points = []
    try:
        points = extract(nz)
    except Exception:
        pass
    g = reduce(gcd, nz)
    qs = [p // g for p in nz if p // g != 0] if g.degree() > 0 else []
    if qs:
        # The rows are LLL-ordered, so the quotients h/g that still vanish
        # at the root concentrate at the front of the list; growing prefixes
        # walk the generated ideal down toward dimension 0.
        # (a) prefix sweep over the quotients,
        # (b) g paired with each single quotient (g itself vanishes at the
        #     root whenever all rows do; one more vanishing equation then
        #     cuts the solution set down to finitely many points).
        n = ring.ngens()
        ks = list(range(n, min(len(qs), n + 5) + 1))
        if len(qs) not in ks:
            ks.append(len(qs))
        systems = [qs[:k] for k in ks] + [[g, q_] for q_ in qs]
        for system in systems:
            try:
                points += extract(system)
            except Exception:
                pass
    out = []
    for pt in points:
        if any(v.denominator() != 1 for v in pt):
            continue
        pt = tuple(ZZ(v) for v in pt)
        if all(abs(pt[j]) <= bounds[j] for j in range(len(pt))) and verify(pt):
            out.append(pt)
    if not out and ring.ngens() == 2 and len(nz) >= 1 and min(bounds) <= 2**20:
        # Last resort for two variables: the shortest polynomial that
        # vanishes at the root defines a curve through it; enumerate one
        # variable over its bound and solve for the other one.
        j = 0 if bounds[0] <= bounds[1] else 1
        vj = ring.gen(j)
        for t in range(-int(bounds[j]), int(bounds[j]) + 1):
            uni = nz[0].subs({vj: t})
            if uni.degree() < 1:
                continue
            for r0, _ in uni.univariate_polynomial().roots(ring=ZZ):
                pt = (t, r0) if j == 0 else (r0, t)
                if abs(pt[1 - j]) <= bounds[1 - j] and verify(pt):
                    out.append(pt)
    if not out and ring.ngens() == 3 and min(bounds) <= 2**14:
        # Last resort for three variables: two vanishing polynomials cut
        # the root down to a curve; enumerate the tightest-bounded
        # variable and eliminate another one by resultant per candidate.
        pair = [g, qs[0]] if qs and g.degree() > 0 else nz[:2]
        _curve_enum_3vars(pair, ring, bounds, verify, out)
    return sorted(set(out))


def _curve_enum_3vars(pair, ring, bounds, verify, out):
    j = min(range(3), key=lambda i: bounds[i])
    others = [i for i in range(3) if i != j]
    vj = ring.gen(j)
    for t in range(-int(bounds[j]), int(bounds[j]) + 1):
        pa = pair[0].subs({vj: ZZ(t)})
        pb = pair[1].subs({vj: ZZ(t)})
        if pa == 0 or pb == 0:
            continue
        for (ka, kb) in (tuple(others), tuple(others)[::-1]):
            # eliminate va = gen(ka); the resultant is univariate in
            # vb = gen(kb), so its roots are vb-values.
            va, vb = ring.gen(ka), ring.gen(kb)
            res = pa.resultant(pb, va)
            if res == 0:
                continue
            ru = res.univariate_polynomial()
            if ru.degree() < 1:
                continue
            for b0, _ in ru.roots(ring=ZZ):
                last = pa.subs({vb: ZZ(b0)})
                if last.degree() < 1:
                    last = pb.subs({vb: ZZ(b0)})
                    if last.degree() < 1:
                        continue
                for a0, _ in last.univariate_polynomial().roots(ring=ZZ):
                    coord = [None, None, None]
                    coord[j], coord[ka], coord[kb] = t, a0, b0
                    pt = tuple(coord)
                    if all(abs(pt[k]) <= bounds[k] for k in range(3)) \
                            and verify(pt):
                        out.append(pt)
            break


# ---------------------------------------------------------------------------
# lattice planning (analytic determinant via the triangular structure)
# ---------------------------------------------------------------------------

def _modular_lattice_plan(f, N, bounds, m, t):
    r"""Return ``(omega, det, hgbound)`` of the modular lattice for
    ``(m, t)`` *without* building any polynomial: the basis is triangular
    in the shift monomials (Jochemsz--May, Section 2.1), so the
    determinant is the product of the diagonal entries
    `N^{m-k} \cdot \prod_j X_j^{u_j}`.  Verified numerically against the
    built matrix on bivariate, trivariate and univariate instances."""
    ring = f.parent()
    n = ring.ngens()
    N = ZZ(N)
    l = _leading_monomial(f)
    base = [e for e, _ in _terms(f)]
    levels = _closure(base, m)
    lMk = []
    for k in range(m + 1):
        lk = tuple(k * c for c in l)
        Mk = {u for u in levels[m]
              if min(_tsub(u, lk)) >= 0 and _tsub(u, lk) in levels[m - k]}
        lMk.append(Mk)
    if t:
        lMk = [{tuple([u[0] + j] + list(u[1:]))
                for u in Mk for j in range(t + 1)} for Mk in lMk]
    lMk.append(set())
    omega = 0
    det = ZZ(1)
    for k in range(m + 1):
        for u in lMk[k] - lMk[k + 1]:
            omega += 1
            det *= N**(m - k) * prod(bounds[j]**u[j] for j in range(n))
    return omega, det, N**m


def _integer_lattice_plan(f, bounds, m, t):
    r"""Return ``(omega, det, hgbound)`` of the Coron-reformulated integer
    lattice for ``(m, t)``, analytically as in
    :func:`_modular_lattice_plan`."""
    ring = f.parent()
    n = ring.ngens()
    bounds = [ZZ(b) for b in bounds]
    dj = [f.degree(ring.gen(j)) for j in range(n)]
    W = max(abs(c) * prod(bounds[j]**e[j] for j in range(n))
            for e, c in _terms(f))
    a0 = ZZ(f.constant_coefficient())
    while gcd(a0, W) != 1:
        W += 1
    R = W * prod(bounds[j]**(dj[j] * (m - 1)) for j in range(n))
    base = [e for e, _ in _terms(f)]
    levels = _closure(base, m)
    S = set(levels[m - 1])
    if t:
        S = {tuple([u[0] + j] + list(u[1:])) for u in S for j in range(t + 1)}
    M = set()
    for u in S:
        for b in base:
            M.add(_tadd(u, b))
    l = [max(u[j] for u in S) for j in range(n)]
    det = ZZ(1)
    for u in S:
        for j in range(n):
            det *= bounds[j]**(l[j] - u[j])
    det *= R**len(M - S)
    return len(M), det, R


def _hg_feasible(omega, det, hgbound):
    r"""Sufficient LLL/Howgrave--Graham success condition
    `2^{(\omega-1)/4} \det(L)^{1/\omega} \sqrt{\omega} < hgbound`,
    compared in bits."""
    from math import log2
    lhs = (det.nbits() - 1)         + omega * (log2(omega) / 2 + (omega - 1) / 4)
    return lhs < omega * (hgbound.nbits() - 1)


# ---------------------------------------------------------------------------
# modular case (paper Section 2.1)
# ---------------------------------------------------------------------------

def small_roots_modular(f, N, bounds, m=2, t=0):
    r"""Return small modular roots of ``f`` modulo ``N`` inside ``bounds``.

    Implements the Jochemsz--May basic strategy (and the extended strategy
    for ``t > 0``, extra shifts of the first variable) for finding roots
    `(x^{(0)}_1,\dots,x^{(0)}_n)` of `f(x_1,\dots,x_n) \equiv 0 \pmod N`.

    INPUT:

    - ``f`` -- polynomial in ``n`` variables over `\ZZ`, not zero modulo
      ``N``, with a small root modulo ``N``
    - ``N`` -- positive integer, the modulus (need not be prime)
    - ``bounds`` -- list `[X_1, \dots, X_n]`, the root satisfies
      `|x^{(0)}_j| \le X_j`
    - ``m`` -- positive integer lattice parameter (default 2); larger
      ``m`` grows the lattice and approaches the asymptotic bound;
      ``None`` auto-selects: the lattice determinant is computed
      analytically from the triangular structure for each `m` (and, for
      more than one variable, a few `t`), and the smallest parameter set
      meeting the sufficient LLL condition is run first, escalating on
      failure
    - ``t`` -- nonnegative integer, number of extra shifts of the first
      variable (extended strategy, default 0)

    OUTPUT:

    A sorted list of tuples of integers `(x^{(0)}_1,\dots,x^{(0)}_n)` such
    that `f(x^{(0)}) \equiv 0 \pmod N` and `|x^{(0)}_j| \le X_j`.  For more
    than one variable the method is heuristic (Assumption 1 of [JM06]_);
    roots are found when the polynomials produced by LLL are algebraically
    independent.  For two variables, if the reduced polynomials share a
    common factor, the shortest one is used as a curve through the root
    and the variable with the smaller bound is enumerated (documented
    last-resort; keep `X_1` moderate in that case).

    For a single variable, the call is delegated to Sage's optimized
    univariate :meth:`~sage.rings.polynomial.polynomial_modn_dense_ntl.Polynomial_dense_modn_ntl_ZZ_fc.small_roots` (the ``m`` and ``t``
parameters are then ignored); the generic lattice path runs as a
fallback if that fails.

    EXAMPLES::

        sage: R = PolynomialRing(ZZ, 'x'); x = R.gen()
        sage: N = next_prime(2**30) * next_prime(2**31)
        sage: f = x**3 + 15*x + 8 - (999**3 + 15*999 + 8) % N
        sage: small_roots_modular(f, N, [2**12], m=4)
        [(999,)]

    The result agrees with the univariate builtin::

        sage: fb = PolynomialRing(Zmod(N), 'x').gen()
        sage: ZZ(fb(f).small_roots(X=2**12, beta=1.0)[0]) == 999
        True

    A bivariate polynomial of Boneh--Durfee type `1 + x(A + y)`, whose root
    `(k, -(p+q))` modulo `e` reveals the factorization of `p*q`::

        sage: p = next_prime(10**12); q = next_prime(10**12 + 500)
        sage: e = inverse_mod(257, (p - 1)*(q - 1))
        sage: R2 = PolynomialRing(ZZ, 'x, y'); x, y = R2.gens()
        sage: g = 1 + x*(p*q + 1 + y)
        sage: roots = small_roots_modular(g, e, [2**9, 3*isqrt(p*q)], m=2)
        sage: sorted(k for k, _ in roots) == [(257*e - 1)//((p-1)*(q-1))]
        True

    ``m = None`` escalates the lattice parameter automatically::

        sage: roots = small_roots_modular(g, e, [2**9, 3*isqrt(p*q)], m=None)
        sage: sorted(k for k, _ in roots) == [(257*e - 1)//((p-1)*(q-1))]
        True

    For three or more variables the method works the same way, but the
    algebraic-independence assumption (Assumption 1 of [JM06]) fails more
    often on small lattices: success then depends on the instance and may
    require larger ``m``, extra shifts (``t > 0``), or strategy tuning.
    Following the paper, experiments must be done for the specific
    polynomial at hand.

    REFERENCES:

    .. [JM06] E. Jochemsz, A. May, "A Strategy for Finding Roots of
       Multivariate Polynomials with New Applications in Attacking RSA
       Variants", ASIACRYPT 2006, LNCS 4284, pp. 267--282.
    """
    ring = f.parent()
    if ring.ngens() == 1:
        try:
            Rz = PolynomialRing(Zmod(ZZ(N)), ring.variable_names())
            roots = sorted({(ZZ(r),) for r in
                            Rz(f).small_roots(X=ZZ(bounds[0]), beta=1.0)})
            if roots:
                return roots
        except Exception:
            pass
    if m is None:
        # Analytic auto-tuning: run the lattice for the smallest (m, t)
        # whose triangular determinant meets the sufficient LLL condition,
        # escalating on failure (the condition is sufficient, not
        # necessary, and extraction can still fail on degenerate shapes).
        tried = False
        for mm in range(2, 7):
            for tt in sorted({0, mm // 4} if ring.ngens() > 1 else {0}):
                omega, det, hg = _modular_lattice_plan(f, N, bounds, mm, tt)
                if _hg_feasible(omega, det, hg):
                    tried = True
                    roots = _small_roots_modular_lattice(f, N, bounds, mm, tt)
                    if roots:
                        return roots
        if not tried:  # nothing feasible in range: best effort at the top
            roots = _small_roots_modular_lattice(
                f, N, bounds, 6, 1 if ring.ngens() > 1 else 0)
            if roots:
                return roots
        return []
    return _small_roots_modular_lattice(f, N, bounds, m, t)


def _small_roots_modular_lattice(f, N, bounds, m, t):
    r"""
    Jochemsz--May lattice path of :func:`small_roots_modular`; see there
    for the documentation and examples.
    """

    ring = f.parent()
    n = ring.ngens()
    N = ZZ(N)
    bounds = [ZZ(b) for b in bounds]
    l = _leading_monomial(f)
    al = _coeff(f, l)
    if gcd(al, N) != 1:
        raise ValueError("leading coefficient shares a factor with N; "
                         "a factor of N is %s" % gcd(al, N))
    inv = al.inverse_mod(N)
    fp = _reduce_mod(_from_dict(ring, {e: inv * c
                                       for e, c in _terms(f)}), N)

    base = [e for e, _ in _terms(fp)]
    levels = _closure(base, m)
    lMk = []
    for k in range(m + 1):
        lk = tuple(k * c for c in l)
        Mk = set()
        for u in levels[m]:
            d = _tsub(u, lk)
            if min(d) >= 0 and tuple(d) in levels[m - k]:
                Mk.add(u)
        lMk.append(Mk)
    if t:
        lMk = [{tuple([u[0] + j] + list(u[1:]))
                for u in Mk for j in range(t + 1)} for Mk in lMk]
    lMk.append(set())  # M_{m+1} = empty

    # NOTE: f_bar^k must be an honest power of the reduced f_bar; reducing
    # the power's coefficients mod N again would break the congruence
    # g(root) = 0 mod N^m (evaluation is multiplicative, reduction is not).
    shifts = []
    for k in range(m + 1):
        fk = fp**k if k else ring(1)
        for u in lMk[k] - lMk[k + 1]:
            e = _tsub(u, tuple(k * c for c in l))
            shifts.append((u, _from_dict(ring, {e: ZZ(1)}) * fk * N**(m - k)))

    return _solve_lattice(ring, shifts, bounds, N**m,
                          lambda pt: f(*pt) % N == 0)


# ---------------------------------------------------------------------------
# integer case (paper Section 2.2, Coron's reformulation)
# ---------------------------------------------------------------------------

def small_roots_integer(f, bounds, m=2, t=0):
    r"""Return small integer roots of ``f`` inside ``bounds``.

    Implements the Jochemsz--May basic strategy (extended strategy for
    ``t > 0``) for finding integer roots
    `(x^{(0)}_1,\dots,x^{(0)}_n)` of `f(x_1,\dots,x_n) = 0` over `\ZZ`,
    via Coron's reformulation.

    INPUT:

    - ``f`` -- polynomial in ``n`` variables over `\ZZ` with nonzero
      constant term, and a small integer root
    - ``bounds`` -- list `[X_1, \dots, X_n]`, the root satisfies
      `|x^{(0)}_j| \le X_j`
    - ``m`` -- positive integer lattice parameter (default 2); ``None``
      tries `m = 2, 3, 4` in order and returns the first nonempty result
    - ``t`` -- nonnegative integer, extra shifts of the first variable
      (extended strategy, default 0)

    OUTPUT:

    A sorted list of tuples of integers `(x^{(0)}_1,\dots,x^{(0)}_n)` with
    `f(x^{(0)}) = 0` and `|x^{(0)}_j| \le X_j`.  Heuristic for `n \ge 2`
    (Assumption 1 of [JM06]_).

    EXAMPLES:

    Factoring with known high bits: given `N = pq` with `p`, `q` of 40 and
    41 bits and their top 24/25 bits, recover `p` and `q`::

        sage: p = next_prime(2**39); q = next_prime(2**40)
        sage: N = p*q; X = Y = 2**16
        sage: p0 = p - 4242; q0 = q - 3417
        sage: R = PolynomialRing(ZZ, 'x, y'); x, y = R.gens()
        sage: f = (p0 + x)*(q0 + y) - N
        sage: roots = small_roots_integer(f, [X, Y], m=2)
        sage: p0 + roots[0][0] == p and q0 + roots[0][1] == q
        True

    Three variables: factoring `N = p q r` with known approximations of
    each factor (small residuals) works at `m = 2`, via the documented
    curve-enumeration fallback::

        sage: p3 = next_prime(2**40); q3 = next_prime(2**41)
        sage: r3 = next_prime(2**42); N3 = p3*q3*r3
        sage: R3 = PolynomialRing(ZZ, 'x, y, z'); x, y, z = R3.gens()
        sage: f3 = (p3 - 1234 + x)*(q3 - 1567 + y)*(r3 - 1891 + z) - N3
        sage: roots3 = small_roots_integer(f3, [2**11, 2**11, 2**11], m=2)
        sage: (1234, 1567, 1891) in roots3
        True

    REFERENCES:

    .. [JM06] E. Jochemsz, A. May, "A Strategy for Finding Roots of
       Multivariate Polynomials with New Applications in Attacking RSA
       Variants", ASIACRYPT 2006, LNCS 4284, pp. 267--282.
    """
    if m is None:
        tried = False
        for mm in range(2, 6):
            for tt in sorted({0, mm // 4} if f.parent().ngens() > 1 else {0}):
                omega, det, hg = _integer_lattice_plan(f, bounds, mm, tt)
                if _hg_feasible(omega, det, hg):
                    tried = True
                    roots = small_roots_integer(f, bounds, mm, tt)
                    if roots:
                        return roots
        if not tried:
            roots = small_roots_integer(
                f, bounds, 5, 1 if f.parent().ngens() > 1 else 0)
            if roots:
                return roots
        return []
    ring = f.parent()
    n = ring.ngens()
    bounds = [ZZ(b) for b in bounds]
    if f.constant_coefficient() == 0:
        raise NotImplementedError("zero constant term; shift the polynomial "
                                  "first (paper Section 2.2, Coron)")

    dj = [f.degree(ring.gen(j)) for j in range(n)]
    W = max(abs(c) * prod(bounds[j]**e[j] for j in range(n))
            for e, c in _terms(f))
    # W may be tweaked to make a0 invertible modulo R (paper Section 2.2)
    a0 = ZZ(f.constant_coefficient())
    while gcd(a0, W) != 1:
        W += 1
    R = W * prod(bounds[j]**(dj[j] * (m - 1)) for j in range(n))
    inv = a0.inverse_mod(R)
    fp = _reduce_mod(_from_dict(ring, {e: inv * c
                                       for e, c in _terms(f)}), R)

    base = [e for e, _ in _terms(fp)]
    levels = _closure(base, m)
    S = set(levels[m - 1])
    if t:
        S = {tuple([u[0] + j] + list(u[1:])) for u in S for j in range(t + 1)}
    M = set()
    for u in S:
        for b in base:
            M.add(_tadd(u, b))
    l = [max(u[j] for u in S) for j in range(n)]

    shifts = []
    for u in sorted(S):
        e = tuple(l[j] - u[j] for j in range(n))
        w = prod(bounds[j]**e[j] for j in range(n))
        shifts.append((u, _from_dict(ring, {u: ZZ(1)}) * fp * w))
    for u in sorted(M - S):
        shifts.append((u, _from_dict(ring, {u: ZZ(1)}) * R))

    return _solve_lattice(ring, shifts, bounds, R,
                          lambda pt: f(*pt) == 0)
