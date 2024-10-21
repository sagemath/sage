# sage_setup: distribution = sagemath-schemes
r"""
Elliptic curves with prescribed good reduction

Construction of elliptic curves with good reduction outside a finite
set of primes

A theorem of Shafarevich states that, over a number field `K`, given
any finite set `S` of primes of `K`, there are (up to isomorphism)
only a finite set of elliptic curves defined over `K` with good
reduction at all primes outside `S`.  An explicit form of the theorem
with an algorithm for finding this finite set was given in "Finding
all elliptic curves with good reduction outside a given set of primes"
by John Cremona and Mark Lingham, Experimental Mathematics 16 No.3
(2007), 303-312.  The method requires computation of the class and
unit groups of `K` as well as all the `S`-integral points on a
collection of auxiliary elliptic curves defined over `K`.

This implementation (April 2009) is only for the case `K=\QQ`, where in
many cases the determination of the necessary sets of `S`-integral
points is possible.  The main user-level function is
:func:`EllipticCurves_with_good_reduction_outside_S`, defined in
constructor.py.  Users should note carefully the following points:

(1) the number of auxiliary curves to be considered is exponential in
the size of `S` (specifically, `2.6^s` where `s=|S|`).

(2) For some of the auxiliary curves it is impossible at present to
provably find all the `S`-integral points using the current
algorithms, which rely on first finding a basis for their Mordell-Weil
groups using 2-descent.  A warning is output in cases where the set of
points (and hence the final output) is not guaranteed to be complete.
Using the ``proof=False`` flag suppresses these warnings.

EXAMPLES: We find all elliptic curves with good reduction outside 2,
listing the label of each::

    sage: [e.label() for e in EllipticCurves_with_good_reduction_outside_S([2])]  # long time (5s on sage.math, 2013)
    ['32a1',
    '32a2',
    '32a3',
    '32a4',
    '64a1',
    '64a2',
    '64a3',
    '64a4',
    '128a1',
    '128a2',
    '128b1',
    '128b2',
    '128c1',
    '128c2',
    '128d1',
    '128d2',
    '256a1',
    '256a2',
    '256b1',
    '256b2',
    '256c1',
    '256c2',
    '256d1',
    '256d2']

Secondly we try the same with `S={11}`; note that warning messages are
printed without ``proof=False`` (unless the optional database is
installed: two of the auxiliary curves whose Mordell-Weil bases are
required have conductors 13068 and 52272 so are in the database)::

    sage: [e.label() for e in EllipticCurves_with_good_reduction_outside_S([11], proof=False)]  # long time (13s on sage.math, 2011)
    ['11a1', '11a2', '11a3', '121a1', '121a2', '121b1', '121b2', '121c1', '121c2', '121d1', '121d2', '121d3']

AUTHORS:

- John Cremona (6 April 2009): initial version (over `\QQ` only).
"""

# ****************************************************************************
#   Copyright (C) 2009 John Cremona <john.cremona@gmail.com>
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
# ****************************************************************************

from sage.misc.mrange import xmrange
from sage.rings.rational_field import QQ
from .constructor import EllipticCurve, EllipticCurve_from_j


def is_possible_j(j, S=[]):
    r"""
    Test if the rational `j` is a possible `j`-invariant of an
    elliptic curve with good reduction outside `S`.

    .. NOTE::

        The condition used is necessary but not sufficient unless S
        contains both 2 and 3.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_egros import is_possible_j
        sage: is_possible_j(0,[])
        False
        sage: is_possible_j(1728,[])
        True
        sage: is_possible_j(-4096/11,[11])
        True
    """
    j = QQ(j)
    return (j.is_zero() and 3 in S) or (j == 1728) \
        or (j.is_S_integral(S) and j.prime_to_S_part(S).is_nth_power(3)
            and (j - 1728).prime_to_S_part(S).abs().is_square())


def curve_key(E1):
    r"""
    Comparison key for elliptic curves over `\QQ`.

    The key is a tuple:

    - if the curve is in the database: (conductor, 0, label, number)

    - otherwise: (conductor, 1, a_invariants)

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_egros import curve_key
        sage: E = EllipticCurve_from_j(1728)
        sage: curve_key(E)
        (32, 0, 0, 2)
        sage: E = EllipticCurve_from_j(1729)
        sage: curve_key(E)
        (2989441, 1, (1, 0, 0, -36, -1))
    """
    try:
        from sage.databases.cremona import parse_cremona_label, class_to_int
        N, l, k = parse_cremona_label(E1.label())
        return (N, 0, class_to_int(l), k)
    except LookupError:
        return (E1.conductor(), 1, E1.ainvs())


def egros_from_j_1728(S=[]):
    r"""
    Given a list of primes S, returns a list of elliptic curves over `\QQ`
    with j-invariant 1728 and good reduction outside S, by checking
    all relevant quartic twists.

    INPUT:

    - ``S`` -- list of primes (default: empty list)

    .. NOTE::

        Primality of elements of S is not checked, and the output
        is undefined if S is not a list or contains non-primes.

    OUTPUT:

    A sorted list of all elliptic curves defined over `\QQ` with
    `j`-invariant equal to `1728` and with good reduction at
    all primes outside the list ``S``.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_egros import egros_from_j_1728
        sage: egros_from_j_1728([])
        []
        sage: egros_from_j_1728([3])
        []
        sage: [e.cremona_label() for e in egros_from_j_1728([2])]
        ['32a1', '32a2', '64a1', '64a4', '256b1', '256b2', '256c1', '256c2']
    """
    Elist = []
    no2 = 2 not in S
    for ei in xmrange([2] + [4] * len(S)):
        u = QQ.prod(p**e for p, e in zip([-1] + S, ei))
        if no2:
            u *= 4  # make sure 12|val(D,2)
        Eu = EllipticCurve([0, 0, 0, u, 0]).minimal_model()
        if Eu.has_good_reduction_outside_S(S):
            Elist += [Eu]
    Elist.sort(key=curve_key)
    return Elist


def egros_from_j_0(S=[]):
    r"""
    Given a list of primes S, returns a list of elliptic curves over `\QQ`
    with j-invariant 0 and good reduction outside S, by checking all
    relevant sextic twists.

    INPUT:

    - ``S`` -- list of primes (default: empty list)

    .. NOTE::

        Primality of elements of S is not checked, and the output
        is undefined if S is not a list or contains non-primes.

    OUTPUT:

    A sorted list of all elliptic curves defined over `\QQ` with
    `j`-invariant equal to `0` and with good reduction at
    all primes outside the list ``S``.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_egros import egros_from_j_0
        sage: egros_from_j_0([])
        []
        sage: egros_from_j_0([2])
        []
        sage: [e.label() for e in egros_from_j_0([3])]
        ['27a1', '27a3', '243a1', '243a2', '243b1', '243b2']
        sage: len(egros_from_j_0([2,3,5]))  # long time (8s on sage.math, 2013)
        432
    """
    Elist = []
    if 3 not in S:
        return Elist
    no2 = 2 not in S
    for ei in xmrange([2] + [6] * len(S)):
        u = QQ.prod(p**e for p, e in zip([-1] + S, ei))
        if no2:
            u *= 16  # make sure 12|val(D,2)
        Eu = EllipticCurve([0, 0, 0, 0, u]).minimal_model()
        if Eu.has_good_reduction_outside_S(S):
            Elist += [Eu]
    Elist.sort(key=curve_key)
    return Elist


def egros_from_j(j, S=[]):
    r"""
    Given a rational j and a list of primes S, returns a list of
    elliptic curves over `\QQ` with j-invariant j and good reduction
    outside S, by checking all relevant quadratic twists.

    INPUT:

    - ``j`` -- a rational number

    - ``S`` -- list of primes (default: empty list)

    .. NOTE::

        Primality of elements of S is not checked, and the output
        is undefined if S is not a list or contains non-primes.

    OUTPUT:

    A sorted list of all elliptic curves defined over `\QQ` with
    `j`-invariant equal to `j` and with good reduction at
    all primes outside the list ``S``.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_egros import egros_from_j
        sage: [e.label() for e in egros_from_j(0,[3])]
        ['27a1', '27a3', '243a1', '243a2', '243b1', '243b2']
        sage: [e.label() for e in egros_from_j(1728,[2])]
        ['32a1', '32a2', '64a1', '64a4', '256b1', '256b2', '256c1', '256c2']
        sage: elist=egros_from_j(-4096/11,[11])
        sage: [e.label() for e in elist]
        ['11a3', '121d1']
    """
    if j == 1728:
        return egros_from_j_1728(S)

    if j == 0:
        return egros_from_j_0(S)

    # Now j != 0, 1728

    E = EllipticCurve_from_j(j)
    Elist = []

    for ei in xmrange([2] * (1 + len(S))):
        u = QQ.prod(p**e for p, e in zip(reversed([-1] + S), ei))
        Eu = E.quadratic_twist(u).minimal_model()
        if Eu.has_good_reduction_outside_S(S):
            Elist += [Eu]

    Elist.sort(key=curve_key)
    return Elist


def egros_from_jlist(jlist, S=[]):
    r"""
    Given a list of rational j and a list of primes S, returns a list
    of elliptic curves over `\QQ` with j-invariant in the list and good
    reduction outside S.

    INPUT:

    - ``j`` -- list of rational numbers

    - ``S`` -- list of primes (default: empty list)

    .. NOTE::

        Primality of elements of S is not checked, and the output
        is undefined if S is not a list or contains non-primes.

    OUTPUT:

    A sorted list of all elliptic curves defined over `\QQ` with
    `j`-invariant in the list ``jlist`` and with good reduction at
    all primes outside the list ``S``.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_egros import egros_get_j, egros_from_jlist
        sage: jlist=egros_get_j([3])
        sage: elist=egros_from_jlist(jlist,[3])
        sage: [e.label() for e in elist]
        ['27a1', '27a2', '27a3', '27a4', '243a1', '243a2', '243b1', '243b2']
        sage: [e.ainvs() for e in elist]
        [(0, 0, 1, 0, -7),
        (0, 0, 1, -270, -1708),
        (0, 0, 1, 0, 0),
        (0, 0, 1, -30, 63),
        (0, 0, 1, 0, -1),
        (0, 0, 1, 0, 20),
        (0, 0, 1, 0, 2),
        (0, 0, 1, 0, -61)]
    """
    elist = [e for j in jlist for e in egros_from_j(j, S)]
    elist.sort(key=curve_key)
    return elist


def egros_get_j(S=[], proof=None, verbose=False):
    r"""
    Return a list of rational `j` such that all elliptic curves
    defined over `\QQ` with good reduction outside `S` have
    `j`-invariant in the list, sorted by height.

    INPUT:

    - ``S`` -- list of primes (default: empty list)

    - ``proof`` -- boolean (default: ``True``); the MW basis for
      auxiliary curves will be computed with this proof flag

    - ``verbose`` -- boolean (default: ``False``); if ``True``, some
      details of the computation will be output

    .. NOTE::

        Proof flag: The algorithm used requires determining all
        S-integral points on several auxiliary curves, which in turn
        requires the computation of their generators.  This is not
        always possible (even in theory) using current knowledge.

        The value of this flag is passed to the function which
        computes generators of various auxiliary elliptic curves, in
        order to find their S-integral points.  Set to ``False`` if the
        default (``True``) causes warning messages, but note that you can
        then not rely on the set of invariants returned being
        complete.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.ell_egros import egros_get_j
        sage: egros_get_j([])
        [1728]
        sage: egros_get_j([2])  # long time (3s on sage.math, 2013)
        [128, 432, -864, 1728, 3375/2, -3456, 6912, 8000, 10976, -35937/4, 287496, -784446336, -189613868625/128]
        sage: egros_get_j([3])  # long time (3s on sage.math, 2013)
        [0, -576, 1536, 1728, -5184, -13824, 21952/9, -41472, 140608/3, -12288000]
        sage: jlist=egros_get_j([2,3]); len(jlist) # long time (30s)
        83
    """
    if not all(p.is_prime() for p in S):
        raise ValueError("Elements of S must be prime.")

    if proof is None:
        from sage.structure.proof.proof import get_flag
        proof = get_flag(proof, "elliptic_curve")
    else:
        proof = bool(proof)

    if verbose:
        import sys  # so we can flush stdout for debugging

    SS = [-1] + S

    jlist = []
    wcount = 0
    nw = 6**len(S) * 2

    if verbose:
        print("Finding possible j invariants for S = ", S)
        print("Using ", nw, " twists of base curve")
        sys.stdout.flush()

    for ei in xmrange([6] * len(S) + [2]):
        w = QQ.prod(p**e for p, e in zip(reversed(SS), ei))
        wcount += 1
        if verbose:
            print("Curve #", wcount, "/", nw, ":")
            print("w = ", w, "=", w.factor())
            sys.stdout.flush()
        a6 = -1728 * w
        E = EllipticCurve([0, 0, 0, 0, a6])
        # This curve may not be minimal at 2 or 3, but the
        # S-integral_points function requires minimality at primes in
        # S, so we find a new model which is p-minimal at both 2 and 3
        # if they are in S.  Note that the isomorphism between models
        # will preserve S-integrality of points.
        E2 = E.local_minimal_model(2) if 2 in S else E
        E23 = E2.local_minimal_model(3) if 3 in S else E2
        urst = E23.isomorphism_to(E)

        try:
            pts = E23.S_integral_points(S, proof=proof)
        except RuntimeError:
            pts = []
            print("Failed to find S-integral points on ", E23.ainvs())
            if proof:
                if verbose:
                    print("--trying again with proof=False")
                    sys.stdout.flush()
                pts = E23.S_integral_points(S, proof=False)
                if verbose:
                    print("--done")
        if verbose:
            print(len(pts), " S-integral points: ", pts)
            sys.stdout.flush()
        for P in pts:
            P = urst(P)
            x = P[0]
            y = P[1]
            j = x**3 / w
            assert j - 1728 == y**2 / w
            if is_possible_j(j, S):
                if j not in jlist:
                    if verbose:
                        print("Adding possible j = ", j)
                        sys.stdout.flush()
                    jlist += [j]
            else:
                if verbose:
                    print("Discarding illegal j = ", j)
                    sys.stdout.flush()
    return sorted(jlist, key=lambda j: j.height())
