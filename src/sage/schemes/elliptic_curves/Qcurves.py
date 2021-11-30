# -*- coding: utf-8 -*-
r"""
Testing whether elliptic curves over number fields are `\QQ`-curves

AUTHORS:

- John Cremona (February 2021)

The code here implements the algorithm of Cremona and Najman presented
in [CrNa2020]_.
"""

##############################################################################
#       Copyright (C) 2020-2021 John Cremona <john.cremona@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
##############################################################################

from sage.rings.rational_field import QQ
from sage.rings.polynomial.polynomial_ring import polygen

def is_Q_curve(E, maxp=100, certificate=False, verbose=False):
    r"""
    Return whether ``E`` is a `\QQ`-curve, with optional certificate.

    INPUT:

    - ``E`` (elliptic curve) -- an elliptic curve over a number field.

    - ``maxp`` (int, default 100): bound on primes used for checking
      necessary local conditions.  The result will not depend on this,
      but using a larger value may return ``False`` faster.

    - ``certificate`` (bool, default ``False``): if ``True`` then a
      second value is returned giving a certificate for the
      `\QQ`-curve property.

    OUTPUT:

    If ``certificate`` is ``False``: either ``True`` (if `E` is a
    `\QQ`-curve), or ``False``.

    If ``certificate`` is ``True``: a tuple consisting of a boolean
    flag as before and a certificate, defined as follows:

    - when the flag is ``True``, so `E` is a `\QQ`-curve:

        - either {'CM':`D`} where `D` is a negative discriminant, when
          `E` has potential CM with discriminant `D`;

        - otherwise {'CM': `0`, 'core_poly': `f`, 'rho': `\rho`, 'r':
          `r`, 'N': `N`}, when `E` is a non-CM `\QQ`-curve, where the
          core polynomial `f` is an irreducible monic polynomial over
          `QQ` of degree `2^\rho`, all of whose roots are
          `j`-invariants of curves isogenous to `E`, the core level
          `N` is a square-free integer with `r` prime factors which is
          the LCM of the degrees of the isogenies between these
          conjugates.  For example, if there exists a curve `E'`
          isogenous to `E` with `j(E')=j\in\QQ`, then the certificate
          is {'CM':0, 'r':0, 'rho':0, 'core_poly': x-j, 'N':1}.

    - when the flag is ``False``, so `E` is not a `\QQ`-curve, the
      certificate is a prime `p` such that the reductions of `E` at
      the primes dividing `p` are inconsistent with the property of
      being a `\QQ`-curve.  See the ALGORITHM section for details.

    ALGORITHM:

    See [CrNa2020]_ for details.

    1. If `E` has rational `j`-invariant, or has CM, then return
    ``True``.

    2. Replace `E` by a curve defined over `K=\QQ(j(E))`. Let `N` be
    the conductor norm.

    3. For all primes `p\mid N` check that the valuations of `j` at
    all `P\mid p` are either all negative or all non-negative; if not,
    return ``False``.

    4. For `p\le maxp`, `p\not\mid N`, check that either `E` is
    ordinary mod `P` for all `P\mid p`, or `E` is supersingular mod
    `P` for all `P\mid p`; if neither, return ``False``.  If all are
    ordinary, check that the integers `a_P(E)^2-4N(P)` have the same
    square-free part; if not, return ``False``.

    5. Compute the `K`-isogeny class of `E` using the "heuristic"
    option (which is faster, but not guaranteed to be complete).
    Check whether the set of `j`-invariants of curves in the class of
    `2`-power degree contains a complete Galois orbit.  If so, return
    ``True``.

    6. Otherwise repeat step 4 for more primes, and if still
    undecided, repeat Step 5 without the "heuristic" option, to get
    the complete `K`-isogeny class (which will probably be no bigger
    than before).  Now return ``True`` if the set of `j`-invariants of
    curves in the class contains a complete Galois orbit, otherwise
    return ``False``.

    EXAMPLES:

    A non-CM curve over `\QQ` and a CM curve over `\QQ` are both
    trivially `\QQ`-curves::

        sage: from sage.schemes.elliptic_curves.Qcurves import is_Q_curve
        sage: E = EllipticCurve([1,2,3,4,5])
        sage: flag, cert = is_Q_curve(E, certificate=True)
        sage: flag
        True
        sage: cert
        {'CM': 0, 'N': 1, 'core_poly': x, 'r': 0, 'rho': 0}

        sage: E = EllipticCurve(j=8000)
        sage: flag, cert = is_Q_curve(E, certificate=True)
        sage: flag
        True
        sage: cert
        {'CM': -8}

    A non-`\QQ`-curve over a quartic field.  The local data at bad
    primes above `3` is inconsistent::

        sage: from sage.schemes.elliptic_curves.Qcurves import is_Q_curve
        sage: R.<x> = PolynomialRing(QQ)
        sage: K.<a> = NumberField(R([3, 0, -5, 0, 1]))
        sage: E = EllipticCurve([K([-3,-4,1,1]),K([4,-1,-1,0]),K([-2,0,1,0]),K([-621,778,138,-178]),K([9509,2046,-24728,10380])])
        sage: is_Q_curve(E, certificate=True, verbose=True)
        Checking whether Elliptic Curve defined by y^2 + (a^3+a^2-4*a-3)*x*y + (a^2-2)*y = x^3 + (-a^2-a+4)*x^2 + (-178*a^3+138*a^2+778*a-621)*x + (10380*a^3-24728*a^2+2046*a+9509) over Number Field in a with defining polynomial x^4 - 5*x^2 + 3 is a Q-curve
        No: inconsistency at the 2 primes dividing 3
        - potentially multiplicative: [True, False]
        (False, 3)

    A non-`\QQ`-curve over a quadratic field.  The local data at bad
    primes is consistent, but the local test at good primes above `13`
    is not::

        sage: K.<a> = NumberField(R([-10, 0, 1]))
        sage: E = EllipticCurve([K([0,1]),K([-1,-1]),K([0,0]),K([-236,40]),K([-1840,464])])
        sage: is_Q_curve(E, certificate=True, verbose=True)
        Checking whether Elliptic Curve defined by y^2 + a*x*y = x^3 + (-a-1)*x^2 + (40*a-236)*x + (464*a-1840) over Number Field in a with defining polynomial x^2 - 10 is a Q-curve
        Applying local tests at good primes above p<=100
        No: inconsistency at the 2 ordinary primes dividing 13
        - Frobenius discriminants mod squares: [-1, -3]
        No: local test at p=13 failed
        (False, 13)

    A quadratic `\QQ`-curve with CM discriminant `-15` (`j`-invariant not in `\QQ`)::

        sage: from sage.schemes.elliptic_curves.Qcurves import is_Q_curve
        sage: R.<x> = PolynomialRing(QQ)
        sage: K.<a> = NumberField(R([-1, -1, 1]))
        sage: E = EllipticCurve([K([1,0]),K([-1,0]),K([0,1]),K([0,-2]),K([0,1])])
        sage: is_Q_curve(E, certificate=True, verbose=True)
        Checking whether Elliptic Curve defined by y^2 + x*y + a*y = x^3 + (-1)*x^2 + (-2*a)*x + a over Number Field in a with defining polynomial x^2 - x - 1 is a Q-curve
        Yes: E is CM (discriminant -15)
        (True, {'CM': -15})

    An example over `\QQ(\sqrt{2},\sqrt{3})`.  The `j`-invariant is in
    `\QQ(\sqrt{6})`, so computations will be done over that field, and
    in fact there is an isogenous curve with rational `j`, so we have
    a so-called rational `\QQ`-curve::

        sage: K.<a> = NumberField(R([1, 0, -4, 0, 1]))
        sage: E = EllipticCurve([K([-2,-4,1,1]),K([0,1,0,0]),K([0,1,0,0]),K([-4780,9170,1265,-2463]),K([163923,-316598,-43876,84852])])
        sage: flag, cert = is_Q_curve(E, certificate=True)
        sage: flag
        True
        sage: cert
        {'CM': 0, 'N': 1, 'core_degs': [1], 'core_poly': x - 85184/3, 'r': 0, 'rho': 0}

    Over the same field, a so-called strict `\QQ`-curve which is not
    isogenous to one with rational `j`, but whose core field is
    quadratic. In fact the isogeny class over `K` consists of `6`
    curves, four with conjugate quartic `j`-invariants and `2` with
    quadratic conjugate `j`-invariants in `\QQ(\sqrt{3})` (but which
    are not base-changes from the quadratic subfield)::

        sage: E = EllipticCurve([K([0,-3,0,1]),K([1,4,0,-1]),K([0,0,0,0]),K([-2,-16,0,4]),K([-19,-32,4,8])])
        sage: flag, cert = is_Q_curve(E, certificate=True)
        sage: flag
        True
        sage: cert
        {'CM': 0,
        'N': 2,
        'core_degs': [1, 2],
        'core_poly': x^2 - 840064*x + 1593413632,
        'r': 1,
        'rho': 1}

    """
    from sage.rings.number_field.number_field_base import is_NumberField

    if verbose:
        print("Checking whether {} is a Q-curve".format(E))

    try:
        assert is_NumberField(E.base_field())
    except (AttributeError, AssertionError):
        raise TypeError("{} must be an elliptic curve defined over a number field in is_Q_curve()")

    from sage.rings.integer_ring import ZZ
    from sage.arith.functions import lcm
    from sage.libs.pari import pari
    from sage.rings.number_field.number_field import NumberField
    from sage.schemes.elliptic_curves.constructor import EllipticCurve
    from sage.schemes.elliptic_curves.cm import cm_j_invariants_and_orders, is_cm_j_invariant

    # Step 1

    # all curves with rational j-invariant are Q-curves:
    jE = E.j_invariant()
    if jE in QQ:
        if verbose:
            print("Yes: j(E) is in QQ")
        if certificate:
            # test for CM
            for d, f, j in cm_j_invariants_and_orders(QQ):
                if jE == j:
                    return True, {'CM': d*f**2}
            # else not CM
            return True, {'CM': ZZ(0), 'r': ZZ(0), 'rho': ZZ(0), 'N': ZZ(1), 'core_poly': polygen(QQ)}
        else:
            return True

    # CM curves are Q-curves:
    flag, df = is_cm_j_invariant(jE)
    if flag:
        d, f = df
        D = d*f**2
        if verbose:
            print("Yes: E is CM (discriminant {})".format(D))
        if certificate:
            return True, {'CM': D}
        else:
            return True

    # Step 2: replace E by a curve defined over Q(j(E)):

    K = E.base_field()
    jpoly = jE.minpoly()
    if jpoly.degree()<K.degree():
        if verbose:
            print("switching to smaller base field: j's minpoly is {}".format(jpoly))
        f = pari(jpoly).polredbest().sage({'x':jpoly.parent().gen()})
        K2 = NumberField(f, 'b')
        jE = jpoly.roots(K2)[0][0]
        if verbose:
            print("New j is {} over {}, with minpoly {}".format(jE, K2, jE.minpoly()))
        #assert jE.minpoly()==jpoly
        E = EllipticCurve(j=jE)
        K = K2
        if verbose:
            print("New test curve is {}".format(E))

    # Step 3: check primes of bad reduction

    NN = E.conductor().norm()
    for p in NN.support():
        Plist = K.primes_above(p)
        if len(Plist)<2:
            continue
        # pot_mult = potential multiplicative reduction
        pot_mult = [jE.valuation(P) < 0 for P in Plist]
        consistent = all(pot_mult) or not any(pot_mult)
        if not consistent:
            if verbose:
                print("No: inconsistency at the {} primes dividing {}".format(len(Plist),p))
                print("  - potentially multiplicative: {}".format(pot_mult))
            if certificate:
                return False, p
            else:
                return False

    # Step 4 check: primes P of good reduction above p<=B:

    if verbose:
        print("Applying local tests at good primes above p<={}".format(maxp))

    res4, p = Step4Test(E, B=maxp, oldB=0, verbose=verbose)
    if not res4:
        if verbose:
            print("No: local test at p={} failed".format(p))
        if certificate:
            return False, p
        else:
            return False

    if verbose:
        print("...all local tests pass for p<={}".format(maxp))

    # Step 5: compute the (partial) K-isogeny class of E and test the
    # set of j-invariants in the class:

    C = E.isogeny_class(algorithm='heuristic', minimal_models=False)
    jC = [E2.j_invariant() for E2 in C]
    centrejpols = conjugacy_test(jC, verbose=verbose)
    if centrejpols:
        if verbose:
            print("Yes: the isogeny class contains a complete conjugacy class of j-invariants")
        if certificate:
            for f in centrejpols:
                rho = f.degree().valuation(2)
                centre_indices = [i for i,j in enumerate(jC) if f(j) == 0]
                M = C.matrix()
                core_degs = [M[centre_indices[0], i] for i in centre_indices]
                level = lcm(core_degs)
                if level.is_squarefree():
                    r = len(level.prime_divisors())
                    cert = {'CM': ZZ(0), 'core_poly':f, 'rho':rho, 'r':r, 'N':level, 'core_degs':core_degs}
                    return True, cert
            print("No central curve found")
        else:
            return True

    # Now we are undecided.  This can happen if either (1) E is not a
    # Q-curve but we did not use enough primes in Step 4 to detect
    # this, or (2) E is a Q-curve but in Step 5 we did not compute the
    # complete isogeny class.  Case (2) is most unlikely since the
    # heuristic bound used in computing isogeny classes means that we
    # have all isogenous curves linked to E by an isogeny of degree
    # supported on primes<1000.

    # We first rerun Step 4 with a larger bound.

    xmaxp = 10 * maxp
    if verbose:
        print("Undecided after first round, so we apply more local tests, up to {}".format(xmaxp))

    res4, p = Step4Test(E, B=xmaxp, oldB=maxp, verbose=verbose)
    if not res4:
        if verbose:
            print("No: local test at p={} failed".format(p))
        if certificate:
            return False, p
        else:
            return False

    # Now we rerun Step 5 using a rigorous computation of the complete
    # isogeny class.  This will probably contain no more curves than
    # before, in which case -- since we already tested that the set of
    # j-invariants does not contain a complete Galois conjugacy class
    # -- we can deduce that E is not a Q-curve.

    if verbose:
        print("...all local tests pass for p<={}".format(xmaxp))
        print("We now compute the complete isogeny class...")

    Cfull = E.isogeny_class(minimal_models=False)
    jCfull = [E2.j_invariant() for E2 in Cfull]

    if len(jC) == len(jCfull):
        if verbose:
            print("...and find that we already had the complete class:so No")
        if certificate:
            return False, 0
        else:
            return False
    if verbose:
        print("...and find that the class contains {} curves, not just the {} we computed originally".format(len(jCfull), len(jC)))
    centrejpols = conjugacy_test(jCfull, verbose=verbose)
    if cert:
        if verbose:
            print("Yes: the isogeny class contains a complete conjugacy class of j-invariants")
        if certificate:
            return True, centrejpols
        else:
            return True
    if verbose:
        print("No: the isogeny class does *not* contain a complete conjugacy class of j-invariants")
    if certificate:
        return False, 0
    else:
        return False

def Step4Test(E, B, oldB=0, verbose=False):
    r"""
    Apply local Q-curve test to E at all primes up to B.

    INPUT:

    - `E` (elliptic curve): an elliptic curve defined over a number field

    - `B` (integer): upper bound on primes to test

    - `oldB` (integer, default 0): lower bound on primes to test

    - `verbose` (boolean, default ``False``): verbosity flag

    OUTPUT:

    Either (``False``, `p`), if the local test at `p` proves that `E`
    is not a `\QQ`-curve, or (``True``, `0`) if all local tests at
    primes between ``oldB`` and ``B`` fail to prove that `E` is not a
    `\QQ`-curve.

    ALGORITHM (see [CrNa2020]_ for details):

    This local test at `p` only applies if `E` has good reduction at
    all of the primes lying above `p` in the base field `K` of `E`.  It
    tests whether (1) `E` is either ordinary at all `P\mid p`, or
    supersingular at all; (2) if ordinary at all, it tests that the
    squarefree part of `a_P^2-4N(P)` is the same for all `P\mid p`.

    EXAMPLES:

    A non-`\QQ`-curve over a quartic field (with LMFDB label
    '4.4.8112.1-12.1-a1') fails this test at `p=13`::

        sage: from sage.schemes.elliptic_curves.Qcurves import Step4Test
        sage: R.<x> = PolynomialRing(QQ)
        sage: K.<a> = NumberField(R([3, 0, -5, 0, 1]))
        sage: E = EllipticCurve([K([-3,-4,1,1]),K([4,-1,-1,0]),K([-2,0,1,0]),K([-621,778,138,-178]),K([9509,2046,-24728,10380])])
        sage: Step4Test(E, 100, verbose=True)
        No: inconsistency at the 2 ordinary primes dividing 13
        - Frobenius discriminants mod squares: [-3, -1]
        (False, 13)

    A `\QQ`-curve over a sextic field (with LMFDB label
    '6.6.1259712.1-64.1-a6') passes this test for all `p<100`::

        sage: from sage.schemes.elliptic_curves.Qcurves import Step4Test
        sage: R.<x> = PolynomialRing(QQ)
        sage: K.<a> = NumberField(R([-3, 0, 9, 0, -6, 0, 1]))
        sage: E = EllipticCurve([K([1,-3,0,1,0,0]),K([5,-3,-6,1,1,0]),K([1,-3,0,1,0,0]),K([-139,-129,331,277,-76,-63]),K([2466,1898,-5916,-4582,1361,1055])])
        sage: Step4Test(E, 100, verbose=True)
        (True, 0)

    """
    from sage.arith.misc import primes
    K = E.base_field()
    NN = E.conductor().norm()
    for p in primes(B):
        if p <= oldB or p.divides(NN):
            continue
        Plist = K.primes_above(p)
        if len(Plist) < 2:
            continue

        EmodP = [E.reduction(P) for P in Plist]

        # (a) Check all are ordinary or all supersingular:
        ordinary = [Ei.is_ordinary() for Ei in EmodP]
        consistent = all(ordinary) or not any(ordinary)
        if not consistent:
            if verbose:
                print("No: inconsistency at the {} primes dividing {} ".format(len(Plist),p))
                print("  - ordinary: {}".format(ordinary))
            return False, p

        # (b) Skip if all are supersingular:
        if not ordinary[0]:
            continue

        # else compare a_P^2-4*N(P) which should have the same squarefree part:

        discs = [(Ei.trace_of_frobenius()**2 - 4 * P.norm()).squarefree_part() for P, Ei in zip(Plist, EmodP)]
        if any(d != discs[0] for d in discs[1:]):
            if verbose:
                print("No: inconsistency at the {} ordinary primes dividing {} ".format(len(Plist), p))
                print("  - Frobenius discriminants mod squares: {}".format(discs))
            return False, p
    # Now we have failed to prove that E is not a Q-curve
    return True, 0

def conjugacy_test(jlist, verbose=False):
    r"""
    Test whether a list of algebraic numbers contains a complete
    conjugacy class of 2-power degree.

    INPUT:

    - `jlist` (list): a list of algebraic numbers in the same field

    - `verbose` (boolean, default ``False``): verbosity flag

    OUTPUT:

    A possibly empty list of irreducible polynomials over `\QQ` of
    2-power degree all of whose roots are in the list.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.Qcurves import conjugacy_test
        sage: conjugacy_test([3])
        [x - 3]
        sage: K.<a> = QuadraticField(2)
        sage: conjugacy_test([K(3), a])
        [x - 3]
        sage: conjugacy_test([K(3), 3+a])
        [x - 3]
        sage: conjugacy_test([3+a])
        []
        sage: conjugacy_test([3+a, 3-a])
        [x^2 - 6*x + 7]
        sage: x = polygen(QQ)
        sage: f = x^3-3
        sage: K.<a> = f.splitting_field()
        sage: js = f.roots(K, multiplicities=False)
        sage: conjugacy_test(js)
        []
        sage: f = x^4-3
        sage: K.<a> = NumberField(f)
        sage: js = f.roots(K, multiplicities=False)
        sage: conjugacy_test(js)
        []
        sage: K.<a> = f.splitting_field()
        sage: js = f.roots(K, multiplicities=False)
        sage: conjugacy_test(js)
        [x^4 - 3]

    """
    from sage.sets.set import Set

    # First test to see if the list contains a rational

    jQ = next((j for j in jlist if j in QQ), None)
    if jQ:
        if verbose:
            print("Yes: an isogenous curve has rational j-invariant {}".format(jQ))
        x = polygen(QQ)
        return [x-jQ]

    # If the degree d is odd then we know that none of the
    # j-invariants in the class have 2-power degree, so we can exit.

    K = jlist[0].parent()
    if K.degree() % 2:
        if verbose:
            print("Odd-degree case: no rational j-invariant in the class {}".format(jlist))
        return []

    # If K has no quadratic subfields we can similarly conclude right
    # away.  This is one way of determining this.

    if K(1).descend_mod_power(QQ,2) == [1]:
        if verbose:
            print("No-quadratic-subfield case: no rational j-invariant in the class {}".format(jlist))
        return []

    # compute the minimum polynomials of the j-invariants in the class
    pols = [j.minpoly() for j in jlist]

    # pick out those of 2-power degree
    pols = [f for f in pols if f.degree().prime_to_m_part(2) == 1]

    # If none, there is nothing to do
    if not pols:
        return []

    # see if there's a poly of degree d appearing d times.  NB There
    # may be more than one of these, possibly including some conjugacy
    # classes defined over the core field but not central, so we
    # return all those with the minimal degree.

    mindeg = min([f.degree() for f in pols])
    minpols = [f for f in pols if f.degree() == mindeg]
    centrepols = list(Set([f for f in pols if f.degree() == minpols.count(f)]))
    if centrepols:
        if verbose:
            print("Yes: the isogeny class contains all j-invariants with min poly {}".format(centrepols))
        return centrepols
    if verbose:
        print("No complete conjugacy class of 2-power size found in {}".format(jlist))
    return []
