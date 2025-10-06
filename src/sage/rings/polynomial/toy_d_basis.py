r"""
Educational version of the `d`-Groebner basis algorithm over PIDs

No attempt was made to optimize this algorithm as the emphasis of this
implementation is a clean and easy presentation.

.. NOTE::

    The notion of 'term' and 'monomial' in [BW1993]_ is swapped from the
    notion of those words in Sage (or the other way around, however you
    prefer it). In Sage a term is a monomial multiplied by a
    coefficient, while in [BW1993]_ a monomial is a term multiplied by a
    coefficient. Also, what is called LM (the leading monomial) in
    Sage is called HT (the head term) in [BW1993]_.

EXAMPLES::

    sage: from sage.rings.polynomial.toy_d_basis import d_basis

First, consider an example from arithmetic geometry::

    sage: A.<x,y> = PolynomialRing(ZZ, 2)
    sage: B.<X,Y> = PolynomialRing(Rationals(), 2)
    sage: f = -y^2 - y + x^3 + 7*x + 1
    sage: fx = f.derivative(x)
    sage: fy = f.derivative(y)
    sage: I = B.ideal([B(f), B(fx), B(fy)])
    sage: I.groebner_basis()                                                            # needs sage.libs.singular
    [1]

Since the output is 1, we know that there are no generic
singularities.

To look at the singularities of the arithmetic surface, we need to do
the corresponding computation over `\ZZ`::

    sage: I = A.ideal([f, fx, fy])
    sage: gb = d_basis(I); gb                                                           # needs sage.libs.singular
    [x - 2020, y - 11313, 22627]

    sage: gb[-1].factor()                                                               # needs sage.libs.singular
    11^3 * 17

This Groebner Basis gives a lot of information.  First, the only
fibers (over `\ZZ`) that are not smooth are at 11 = 0, and 17 = 0.
Examining the Groebner Basis, we see that we have a simple node in
both the fiber at 11 and at 17.  From the factorization, we see that
the node at 17 is regular on the surface (an `I_1` node), but the node
at 11 is not.  After blowing up this non-regular point, we find that
it is an `I_3` node.

Another example. This one is from the Magma Handbook::

    sage: # needs sage.libs.singular
    sage: P.<x, y, z> = PolynomialRing(IntegerRing(), 3, order='lex')
    sage: I = ideal(x^2 - 1, y^2 - 1, 2*x*y - z)
    sage: I = Ideal(d_basis(I))
    sage: x.reduce(I)
    x
    sage: (2*x).reduce(I)
    y*z

To compute modulo 4, we can add the generator 4 to our basis.::

    sage: # needs sage.libs.singular
    sage: I = ideal(x^2 - 1, y^2 - 1, 2*x*y - z, 4)
    sage: gb = d_basis(I)
    sage: R = P.change_ring(IntegerModRing(4))
    sage: gb = [R(f) for f in gb if R(f)]; gb
    [x^2 - 1, x*z + 2*y, 2*x - y*z, y^2 - 1, z^2, 2*z]

A third example is also from the Magma Handbook.

This example shows how one can use Groebner bases over the integers to
find the primes modulo which a system of equations has a solution,
when the system has no solutions over the rationals.

We first form a certain ideal `I` in `\ZZ[x, y, z]`, and note that the
Groebner basis of `I` over `\QQ` contains 1, so there are no solutions
over `\QQ` or an algebraic closure of it (this is not surprising as
there are 4 equations in 3 unknowns). ::

    sage: P.<x, y, z> = PolynomialRing(IntegerRing(), 3, order='degneglex')
    sage: I = ideal(x^2 - 3*y, y^3 - x*y, z^3 - x, x^4 - y*z + 1)
    sage: I.change_ring(P.change_ring(RationalField())).groebner_basis()                # needs sage.libs.singular
    [1]

However, when we compute the Groebner basis of I (defined over `\ZZ`), we
note that there is a certain integer in the ideal which is not 1::

    sage: gb = d_basis(I); gb                                                           # needs sage.libs.singular
    [z ..., y ..., x ..., 282687803443]

Now for each prime `p` dividing this integer 282687803443, the Groebner
basis of I modulo `p` will be non-trivial and will thus give a solution
of the original system modulo `p`.::

    sage: factor(282687803443)
    101 * 103 * 27173681

    sage: I.change_ring(P.change_ring(GF(101))).groebner_basis()                        # needs sage.libs.singular
    [z - 33, y + 48, x + 19]

    sage: I.change_ring(P.change_ring(GF(103))).groebner_basis()                        # needs sage.libs.singular
    [z - 18, y + 8, x + 39]

    sage: I.change_ring(P.change_ring(GF(27173681))).groebner_basis()                   # needs sage.libs.singular sage.rings.finite_rings
    [z + 10380032, y + 3186055, x - 536027]

Of course, modulo any other prime the Groebner basis is trivial so
there are no other solutions. For example::

    sage: I.change_ring(P.change_ring(GF(3))).groebner_basis()                          # needs sage.libs.singular
    [1]

AUTHOR:

- Martin Albrecht (2008-08): initial version
"""
from sage.rings.integer_ring import ZZ
from sage.arith.functions import lcm
from sage.arith.misc import XGCD as xgcd, GCD as gcd
from sage.rings.polynomial.toy_buchberger import inter_reduction
from sage.structure.sequence import Sequence


def spol(g1, g2):
    """
    Return the S-Polynomial of ``g_1`` and ``g_2``.

    Let `a_i t_i` be `LT(g_i)`, `b_i = a/a_i` with `a = LCM(a_i,a_j)`,
    and `s_i = t/t_i` with `t = LCM(t_i,t_j)`. Then the S-Polynomial
    is defined as: `b_1s_1g_1 - b_2s_2g_2`.

    INPUT:

    - ``g1`` -- polynomial
    - ``g2`` -- polynomial

    EXAMPLES::

        sage: from sage.rings.polynomial.toy_d_basis import spol
        sage: P.<x, y, z> = PolynomialRing(IntegerRing(), 3, order='lex')
        sage: f = x^2 - 1
        sage: g = 2*x*y - z
        sage: spol(f,g)
        x*z - 2*y
    """
    a1, a2 = g1.lc(), g2.lc()
    a = a1.lcm(a2)
    b1, b2 = a // a1, a // a2

    t1, t2 = g1.lm(), g2.lm()
    t = t1.parent().monomial_lcm(t1, t2)
    s1, s2 = t // t1, t // t2

    return b1 * s1 * g1 - b2 * s2 * g2


def gpol(g1, g2):
    """
    Return the G-Polynomial of ``g_1`` and ``g_2``.

    Let `a_i t_i` be `LT(g_i)`, `a = a_i*c_i + a_j*c_j` with `a =
    GCD(a_i,a_j)`, and `s_i = t/t_i` with `t = LCM(t_i,t_j)`. Then the
    G-Polynomial is defined as: `c_1s_1g_1 - c_2s_2g_2`.

    INPUT:

    - ``g1`` -- polynomial
    - ``g2`` -- polynomial

    EXAMPLES::

        sage: from sage.rings.polynomial.toy_d_basis import gpol
        sage: P.<x, y, z> = PolynomialRing(IntegerRing(), 3, order='lex')
        sage: f = x^2 - 1
        sage: g = 2*x*y - z
        sage: gpol(f,g)
        x^2*y - y
    """
    a1, a2 = g1.lc(), g2.lc()
    a, c1, c2 = xgcd(a1, a2)

    t1, t2 = g1.lm(), g2.lm()
    t = t1.parent().monomial_lcm(t1, t2)
    s1, s2 = t // t1, t // t2

    return c1 * s1 * g1 + c2 * s2 * g2


def LM(f):
    return f.lm()


def LC(f):
    return f.lc()


def d_basis(F, strat=True):
    r"""
    Return the `d`-basis for the Ideal ``F`` as defined in [BW1993]_.

    INPUT:

    - ``F`` -- an ideal
    - ``strat`` -- use update strategy (default: ``True``)

    EXAMPLES::

        sage: from sage.rings.polynomial.toy_d_basis import d_basis
        sage: A.<x,y> = PolynomialRing(ZZ, 2)
        sage: f = -y^2 - y + x^3 + 7*x + 1
        sage: fx = f.derivative(x)
        sage: fy = f.derivative(y)
        sage: I = A.ideal([f,fx,fy])
        sage: gb = d_basis(I); gb                                                       # needs sage.libs.singular
        [x - 2020, y - 11313, 22627]
    """
    R = F.ring()

    G = set(inter_reduction(F.gens()))
    B = set((f1, f2) for f1 in G for f2 in G if f1 != f2)
    D = set()
    C = set(B)

    LCM = R.monomial_lcm
    divides = R.monomial_divides

    def divides_ZZ(x, y):
        return ZZ(x).divides(ZZ(y))

    while B:
        while C:
            f1, f2 = select(C)
            C.remove((f1, f2))
            lcm_lmf1_lmf2 = LCM(LM(f1), LM(f2))
            if not any(divides(LM(g), lcm_lmf1_lmf2) and
                       divides_ZZ(LC(g), LC(f1)) and
                       divides_ZZ(LC(g), LC(f2))
                       for g in G):
                h = gpol(f1, f2)
                h0 = h.reduce(G)
                if h0.lc() < 0:
                    h0 *= -1
                if not strat:
                    D = D.union([(g, h0) for g in G])
                    G.add(h0)
                else:
                    G, D = update(G, D, h0)
                G = inter_reduction(G)

        f1, f2 = select(B)
        B.remove((f1, f2))
        h = spol(f1, f2)
        h0 = h.reduce(G)
        if h0 != 0:
            if h0.lc() < 0:
                h0 *= -1
            if not strat:
                D = D.union([(g, h0) for g in G])
                G.add(h0)
            else:
                G, D = update(G, D, h0)

        B = B.union(D)
        C = D
        D = set()

    return Sequence(sorted(inter_reduction(G), reverse=True))


def select(P):
    """
    The normal selection strategy.

    INPUT:

    - ``P`` -- list of critical pairs

    OUTPUT: an element of P

    EXAMPLES::

        sage: from sage.rings.polynomial.toy_d_basis import select
        sage: A.<x,y> = PolynomialRing(ZZ, 2)
        sage: f = -y^2 - y + x^3 + 7*x + 1
        sage: fx = f.derivative(x)
        sage: fy = f.derivative(y)
        sage: G = [f, fx, fy]
        sage: B = set((f1, f2) for f1 in G for f2 in G if f1 != f2)
        sage: select(B)
        (-2*y - 1, 3*x^2 + 7)
    """
    min_d = 2**20
    min_pair = 0, 0
    for fi, fj in sorted(P):
        d = fi.parent().monomial_lcm(fi.lm(), fj.lm()).total_degree()
        if d < min_d:
            min_d = d
            min_pair = fi, fj
    return min_pair


def update(G, B, h):
    """
    Update ``G`` using the list of critical pairs ``B`` and the
    polynomial ``h`` as presented in [BW1993]_, page 230. For this,
    Buchberger's first and second criterion are tested.

    This function uses the Gebauer-Moeller Installation.

    INPUT:

    - ``G`` -- an intermediate Groebner basis
    - ``B`` -- list of critical pairs
    - ``h`` -- a polynomial

    OUTPUT: ``G,B`` where ``G`` and ``B`` are updated

    EXAMPLES::

        sage: from sage.rings.polynomial.toy_d_basis import update
        sage: A.<x,y> = PolynomialRing(ZZ, 2)
        sage: G = set([3*x^2 + 7, 2*y + 1, x^3 - y^2 + 7*x - y + 1])
        sage: B = set()
        sage: h = x^2*y - x^2 + y - 3
        sage: update(G,B,h)
        ({2*y + 1, 3*x^2 + 7, x^2*y - x^2 + y - 3, x^3 - y^2 + 7*x - y + 1},
         {(x^2*y - x^2 + y - 3, 2*y + 1),
          (x^2*y - x^2 + y - 3, 3*x^2 + 7),
          (x^2*y - x^2 + y - 3, x^3 - y^2 + 7*x - y + 1)})
    """
    R = h.parent()
    LCM = R.monomial_lcm

    def lt_divides(x, y):
        return R.monomial_divides(LM(h), LM(g)) and LC(h).divides(LC(g))

    def lt_pairwise_prime(x, y):
        return (R.monomial_pairwise_prime(LM(x), LM(y))
                and gcd(LC(x), LC(y)) == 1)

    def lcm_divides(f, g1, h):
        return (R.monomial_divides(LCM(LM(h), LM(f[1])), LCM(LM(h), LM(g1)))
                and lcm(LC(h), LC(f[1])).divides(lcm(LC(h), LC(g1))))

    C = set((h, g) for g in G)

    D = set()
    while C:
        (h, g1) = C.pop()

        if (lt_pairwise_prime(h, g1) or
                (not any(lcm_divides(f, g1, h) for f in C) and
                 not any(lcm_divides(f, g1, h) for f in D))):
            D.add((h, g1))

    E = set()

    while D:
        (h, g) = D.pop()
        if not lt_pairwise_prime(h, g):
            E.add((h, g))

    B_new = set()
    while B:
        g1, g2 = B.pop()

        lcm_12 = lcm(LC(g1), LC(g2)) * LCM(LM(g1), LM(g2))
        if (not lt_divides(lcm_12, h) or
                lcm(LC(g1), LC(h)) * R.monomial_lcm(LM(g1), LM(h)) == lcm_12 or
                lcm(LC(h), LC(g2)) * R.monomial_lcm(LM(h), LM(g2)) == lcm_12):
            B_new.add((g1, g2))

    B_new = B_new.union(E)

    G_new = set()
    while G:
        g = G.pop()
        if not lt_divides(g, h):
            G_new.add(g)

    G_new.add(h)

    return G_new, B_new
