r"""
Educational versions of Groebner basis algorithms

Following [BW1993]_, the original Buchberger algorithm (algorithm GROEBNER in
[BW1993]_) and an improved version of Buchberger's algorithm (algorithm
GROEBNERNEW2 in [BW1993]_) are implemented.

No attempt was made to optimize either algorithm as the emphasis of these
implementations is a clean and easy presentation. To compute a Groebner basis
most efficiently in Sage, use the :meth:`.MPolynomialIdeal.groebner_basis`
method on multivariate polynomial objects instead.

.. NOTE::

   The notion of 'term' and 'monomial' in [BW1993]_ is swapped from the
   notion of those words in Sage (or the other way around, however you
   prefer it). In Sage a term is a monomial multiplied by a
   coefficient, while in [BW1993]_ a monomial is a term multiplied by a
   coefficient. Also, what is called LM (the leading monomial) in
   Sage is called HT (the head term) in [BW1993]_.

EXAMPLES:

Consider Katsura-6 with respect to a ``degrevlex`` ordering. ::

    sage: from sage.rings.polynomial.toy_buchberger import *
    sage: P.<a,b,c,e,f,g,h,i,j,k> = PolynomialRing(GF(32003))
    sage: I = sage.rings.ideal.Katsura(P, 6)

    sage: g1 = buchberger(I)
    sage: g2 = buchberger_improved(I)
    sage: g3 = I.groebner_basis()

All algorithms actually compute a Groebner basis::

    sage: Ideal(g1).basis_is_groebner()
    True
    sage: Ideal(g2).basis_is_groebner()
    True
    sage: Ideal(g3).basis_is_groebner()
    True

The results are correct::

    sage: Ideal(g1) == Ideal(g2) == Ideal(g3)
    True

If ``get_verbose()`` is `\ge 1`, a protocol is provided::

    sage: from sage.misc.verbose import set_verbose
    sage: set_verbose(1)
    sage: P.<a,b,c> = PolynomialRing(GF(127))
    sage: I = sage.rings.ideal.Katsura(P)
    // sage... ideal

    sage: I
    Ideal (a + 2*b + 2*c - 1, a^2 + 2*b^2 + 2*c^2 - a, 2*a*b + 2*b*c - b) of Multivariate Polynomial Ring in a, b, c over Finite Field of size 127

    sage: buchberger(I)  # random
    (a + 2*b + 2*c - 1, a^2 + 2*b^2 + 2*c^2 - a) => -2*b^2 - 6*b*c - 6*c^2 + b + 2*c
    G: set([a + 2*b + 2*c - 1, 2*a*b + 2*b*c - b, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c])
    <BLANKLINE>
    (a^2 + 2*b^2 + 2*c^2 - a, a + 2*b + 2*c - 1) => 0
    G: set([a + 2*b + 2*c - 1, 2*a*b + 2*b*c - b, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c])
    <BLANKLINE>
    (a + 2*b + 2*c - 1, 2*a*b + 2*b*c - b) => -5*b*c - 6*c^2 - 63*b + 2*c
    G: set([a + 2*b + 2*c - 1, 2*a*b + 2*b*c - b, -5*b*c - 6*c^2 - 63*b + 2*c, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c])
    <BLANKLINE>
    (2*a*b + 2*b*c - b, a + 2*b + 2*c - 1) => 0
    G: set([a + 2*b + 2*c - 1, 2*a*b + 2*b*c - b, -5*b*c - 6*c^2 - 63*b + 2*c, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c])
    <BLANKLINE>
    (2*a*b + 2*b*c - b, -5*b*c - 6*c^2 - 63*b + 2*c) => -22*c^3 + 24*c^2 - 60*b - 62*c
    G: set([a + 2*b + 2*c - 1, -22*c^3 + 24*c^2 - 60*b - 62*c, 2*a*b + 2*b*c - b, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c, -5*b*c - 6*c^2 - 63*b + 2*c])
    <BLANKLINE>
    (2*a*b + 2*b*c - b, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c) => 0
    G: set([a + 2*b + 2*c - 1, -22*c^3 + 24*c^2 - 60*b - 62*c, 2*a*b + 2*b*c - b, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c, -5*b*c - 6*c^2 - 63*b + 2*c])
    <BLANKLINE>
    (2*a*b + 2*b*c - b, a^2 + 2*b^2 + 2*c^2 - a) => 0
    G: set([a + 2*b + 2*c - 1, -22*c^3 + 24*c^2 - 60*b - 62*c, 2*a*b + 2*b*c - b, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c, -5*b*c - 6*c^2 - 63*b + 2*c])
    <BLANKLINE>
    (a + 2*b + 2*c - 1, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c) => 0
    G: set([a + 2*b + 2*c - 1, -22*c^3 + 24*c^2 - 60*b - 62*c, 2*a*b + 2*b*c - b, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c, -5*b*c - 6*c^2 - 63*b + 2*c])
    <BLANKLINE>
    (a^2 + 2*b^2 + 2*c^2 - a, 2*a*b + 2*b*c - b) => 0
    G: set([a + 2*b + 2*c - 1, -22*c^3 + 24*c^2 - 60*b - 62*c, 2*a*b + 2*b*c - b, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c, -5*b*c - 6*c^2 - 63*b + 2*c])
    <BLANKLINE>
    (-2*b^2 - 6*b*c - 6*c^2 + b + 2*c, -5*b*c - 6*c^2 - 63*b + 2*c) => 0
    G: set([a + 2*b + 2*c - 1, -22*c^3 + 24*c^2 - 60*b - 62*c, 2*a*b + 2*b*c - b, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c, -5*b*c - 6*c^2 - 63*b + 2*c])
    <BLANKLINE>
    (a + 2*b + 2*c - 1, -5*b*c - 6*c^2 - 63*b + 2*c) => 0
    G: set([a + 2*b + 2*c - 1, -22*c^3 + 24*c^2 - 60*b - 62*c, 2*a*b + 2*b*c - b, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c, -5*b*c - 6*c^2 - 63*b + 2*c])
    <BLANKLINE>
    (a^2 + 2*b^2 + 2*c^2 - a, -5*b*c - 6*c^2 - 63*b + 2*c) => 0
    G: set([a + 2*b + 2*c - 1, -22*c^3 + 24*c^2 - 60*b - 62*c, 2*a*b + 2*b*c - b, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c, -5*b*c - 6*c^2 - 63*b + 2*c])
    <BLANKLINE>
    (-5*b*c - 6*c^2 - 63*b + 2*c, -22*c^3 + 24*c^2 - 60*b - 62*c) => 0
    G: set([a + 2*b + 2*c - 1, -22*c^3 + 24*c^2 - 60*b - 62*c, 2*a*b + 2*b*c - b, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c, -5*b*c - 6*c^2 - 63*b + 2*c])
    <BLANKLINE>
    (a + 2*b + 2*c - 1, -22*c^3 + 24*c^2 - 60*b - 62*c) => 0
    G: set([a + 2*b + 2*c - 1, -22*c^3 + 24*c^2 - 60*b - 62*c, 2*a*b + 2*b*c - b, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c, -5*b*c - 6*c^2 - 63*b + 2*c])
    <BLANKLINE>
    (a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c) => 0
    G: set([a + 2*b + 2*c - 1, -22*c^3 + 24*c^2 - 60*b - 62*c, 2*a*b + 2*b*c - b, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c, -5*b*c - 6*c^2 - 63*b + 2*c])
    <BLANKLINE>
    (-2*b^2 - 6*b*c - 6*c^2 + b + 2*c, -22*c^3 + 24*c^2 - 60*b - 62*c) => 0
    G: set([a + 2*b + 2*c - 1, -22*c^3 + 24*c^2 - 60*b - 62*c, 2*a*b + 2*b*c - b, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c, -5*b*c - 6*c^2 - 63*b + 2*c])
    <BLANKLINE>
    (2*a*b + 2*b*c - b, -22*c^3 + 24*c^2 - 60*b - 62*c) => 0
    G: set([a + 2*b + 2*c - 1, -22*c^3 + 24*c^2 - 60*b - 62*c, 2*a*b + 2*b*c - b, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c, -5*b*c - 6*c^2 - 63*b + 2*c])
    <BLANKLINE>
    (a^2 + 2*b^2 + 2*c^2 - a, -22*c^3 + 24*c^2 - 60*b - 62*c) => 0
    G: set([a + 2*b + 2*c - 1, -22*c^3 + 24*c^2 - 60*b - 62*c, 2*a*b + 2*b*c - b, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c, -5*b*c - 6*c^2 - 63*b + 2*c])
    <BLANKLINE>
    15 reductions to zero.
    [a + 2*b + 2*c - 1, -22*c^3 + 24*c^2 - 60*b - 62*c, 2*a*b + 2*b*c - b, a^2 + 2*b^2 + 2*c^2 - a, -2*b^2 - 6*b*c - 6*c^2 + b + 2*c, -5*b*c - 6*c^2 - 63*b + 2*c]

The original Buchberger algorithm performs 15 useless reductions to
zero for this example::

    sage: gb = buchberger(I)
    ...
    15 reductions to zero.

The 'improved' Buchberger algorithm in contrast only performs 1 reduction to
zero::

    sage: gb = buchberger_improved(I)
    ...
    1 reductions to zero.
    sage: sorted(gb)
    [a + 2*b + 2*c - 1, b*c + 52*c^2 + 38*b + 25*c, b^2 - 26*c^2 - 51*b + 51*c, c^3 + 22*c^2 - 55*b + 49*c]

AUTHORS:

- Martin Albrecht (2007-05-24): initial version

- Marshall Hampton (2009-07-08): some doctest additions

"""

from sage.misc.verbose import get_verbose
from sage.structure.sequence import Sequence

# some aliases that conform to Becker and Weispfenning's notation:
LCM = lambda f, g: f.parent().monomial_lcm(f, g)
LM = lambda f: f.lm()
LT = lambda f: f.lt()


def spol(f, g):
    """
    Compute the S-polynomial of f and g.

    INPUT:

    -  ``f, g`` -- polynomials

    OUTPUT: the S-polynomial of f and g

    EXAMPLES::

        sage: R.<x,y,z> = PolynomialRing(QQ)
        sage: from sage.rings.polynomial.toy_buchberger import spol
        sage: spol(x^2 - z - 1, z^2 - y - 1)
        x^2*y - z^3 + x^2 - z^2
    """
    fg_lcm = LCM(LM(f), LM(g))
    return fg_lcm//LT(f)*f - fg_lcm//LT(g)*g


def buchberger(F):
    """
    Compute a Groebner basis using the original version of Buchberger's
    algorithm as presented in [BW1993]_, page 214.

    INPUT:

    - ``F`` -- an ideal in a multivariate polynomial ring

    OUTPUT: a Groebner basis for F

    .. NOTE::

       The verbosity of this function may be controlled with a
       ``set_verbose()`` call. Any value >=1 will result in this
       function printing intermediate bases.

    EXAMPLES::

        sage: from sage.rings.polynomial.toy_buchberger import buchberger
        sage: R.<x,y,z> = PolynomialRing(QQ)
        sage: I = R.ideal([x^2 - z - 1, z^2 - y - 1, x*y^2 - x - 1])
        sage: set_verbose(0)
        sage: gb = buchberger(I)
        sage: gb.is_groebner()
        True
        sage: gb.ideal() == I
        True
    """
    G = set(F.gens())
    B = set((g1, g2) for g1 in G for g2 in G if g1 != g2)

    if get_verbose() >= 1:
        reductions_to_zero = 0

    while B:
        g1, g2 = select(B)
        B.remove((g1, g2))

        h = spol(g1, g2).reduce(G)
        if h != 0:
            B = B.union((g, h) for g in G)
            G.add(h)

        if get_verbose() >= 1:
            print("(%s, %s) => %s" % (g1, g2, h))
            print("G: %s\n" % G)
            if h == 0:
                reductions_to_zero += 1

    if get_verbose() >= 1:
        print("%d reductions to zero." % reductions_to_zero)

    return Sequence(G)


def buchberger_improved(F):
    """
    Compute a Groebner basis using an improved version of Buchberger's
    algorithm as presented in [BW1993]_, page 232.

    This variant uses the Gebauer-Moeller Installation to apply
    Buchberger's first and second criterion to avoid useless pairs.

    INPUT:

    - ``F`` -- an ideal in a multivariate polynomial ring

    OUTPUT: a Groebner basis for F

    .. NOTE::

       The verbosity of this function may be controlled with a
       ``set_verbose()`` call. Any value ``>=1`` will result in this
       function printing intermediate Groebner bases.

    EXAMPLES::

        sage: from sage.rings.polynomial.toy_buchberger import buchberger_improved
        sage: R.<x,y,z> = PolynomialRing(QQ)
        sage: set_verbose(0)
        sage: sorted(buchberger_improved(R.ideal([x^4 - y - z, x*y*z - 1])))
        [x*y*z - 1, x^3 - y^2*z - y*z^2, y^3*z^2 + y^2*z^3 - x^2]
    """
    F = inter_reduction(F.gens())

    G = set()
    B = set()

    if get_verbose() >= 1:
        reductions_to_zero = 0

    while F:
        f = min(F)
        F.remove(f)
        G, B = update(G, B, f)

    while B:

        g1, g2 = select(B)
        B.remove((g1, g2))
        h = spol(g1, g2).reduce(G)
        if h != 0:
            G, B = update(G, B, h)

        if get_verbose() >= 1:
            print("(%s, %s) => %s" % (g1, g2, h))
            print("G: %s\n" % G)
            if h == 0:
                reductions_to_zero += 1

    if get_verbose() >= 1:
        print("%d reductions to zero." % reductions_to_zero)

    return Sequence(inter_reduction(G))


def update(G, B, h):
    """
    Update ``G`` using the set of critical pairs ``B`` and the
    polynomial ``h`` as presented in [BW1993]_, page 230. For this,
    Buchberger's first and second criterion are tested.

    This function implements the Gebauer-Moeller Installation.

    INPUT:

    - ``G`` -- an intermediate Groebner basis

    - ``B`` -- a set of critical pairs

    - ``h`` -- a polynomial

    OUTPUT: a tuple of

    - an intermediate Groebner basis

    - a set of critical pairs

    EXAMPLES::

        sage: from sage.rings.polynomial.toy_buchberger import update
        sage: R.<x,y,z> = PolynomialRing(QQ)
        sage: set_verbose(0)
        sage: update(set(), set(), x*y*z)
        ({x*y*z}, set())
        sage: G, B = update(set(), set(), x*y*z - 1)
        sage: G, B = update(G, B, x*y^2 - 1)
        sage: G, B
        ({x*y*z - 1, x*y^2 - 1}, {(x*y^2 - 1, x*y*z - 1)})
    """
    R = h.parent()

    C = set((h, g) for g in G)
    D = set()

    while C:
        (h, g) = C.pop()

        lcm_divides = lambda rhs: R.monomial_divides(LCM(LM(h), LM(rhs[1])),
                                                     LCM(LM(h), LM(g)))

        if R.monomial_pairwise_prime(LM(h), LM(g)) or \
                (
                   not any(lcm_divides(f) for f in C)
                   and
                   not any(lcm_divides(f) for f in D)
                ):
            D.add((h, g))

    E = set()

    while D:
        (h, g) = D.pop()
        if not R.monomial_pairwise_prime(LM(h), LM(g)):
            E.add((h, g))

    B_new = set()

    while B:
        g1, g2 = B.pop()
        if not R.monomial_divides(LM(h), LCM(LM(g1), LM(g2))) or \
           R.monomial_lcm(LM(g1), LM(h)) == LCM(LM(g1), LM(g2)) or \
           R.monomial_lcm(LM(h), LM(g2)) == LCM(LM(g1), LM(g2)):
            B_new.add((g1, g2))

    B_new = B_new.union(E)

    G_new = set()

    while G:
        g = G.pop()
        if not R.monomial_divides(LM(h), LM(g)):
            G_new.add(g)

    G_new.add(h)

    return G_new, B_new


def select(P):
    """
    Select a polynomial using the normal selection strategy.

    INPUT:

    - ``P`` -- a list of critical pairs

    OUTPUT: an element of P

    EXAMPLES::

        sage: from sage.rings.polynomial.toy_buchberger import select
        sage: R.<x,y,z> = PolynomialRing(QQ, order='lex')
        sage: ps = [x^3 - z -1, z^3 - y - 1, x^5 - y - 2]
        sage: pairs = [[ps[i], ps[j]] for i in range(3) for j in range(i+1, 3)]
        sage: select(pairs)
        [x^3 - z - 1, -y + z^3 - 1]
    """
    return min(P, key=lambda fi_fj: LCM(LM(fi_fj[0]),
                                        LM(fi_fj[1])).total_degree())


def inter_reduction(Q):
    r"""
    Compute inter-reduced polynomials from a set of polynomials.

    INPUT:

    - ``Q`` -- a set of polynomials

    OUTPUT: if ``Q`` is the set `(f_1, ..., f_n)`, this method returns `(g_1,
    ..., g_s)` such that:

    - `<f_1,...,f_n> = <g_1,...,g_s>`
    - `LM(g_i) \neq LM(g_j)` for all `i \neq j`
    - `LM(g_i)` does not divide `m` for all monomials `m` of
      `\{g_1,...,g_{i-1}, g_{i+1},...,g_s\}`
    - `LC(g_i) = 1` for all `i`.

    EXAMPLES::

        sage: from sage.rings.polynomial.toy_buchberger import inter_reduction
        sage: inter_reduction(set())
        set()

    ::

        sage: P.<x,y> = QQ[]
        sage: reduced = inter_reduction(set([x^2 - 5*y^2, x^3]))
        sage: reduced == set([x*y^2, x^2-5*y^2])
        True
        sage: reduced == inter_reduction(set([2*(x^2 - 5*y^2), x^3]))
        True
    """
    if not Q:
        return Q  # if Q is empty we cannot get a base ring
    base_ring = next(iter(Q)).base_ring()

    Q = set(Q)
    while True:
        Qbar = set(Q)
        for p in sorted(Qbar):
            Q.remove(p)
            h = p.reduce(Q)
            if not h.is_zero():
                Q.add(h)
        if Qbar == Q:
            if base_ring.is_field():
                return set(f.lc()**(-1) * f for f in Qbar)
            else:
                return Qbar
