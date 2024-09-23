"Birch and Swinnerton-Dyer formulas"

from sage.arith.misc import prime_divisors
from sage.misc.lazy_import import lazy_import
from sage.rings.infinity import Infinity
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ

lazy_import("sage.functions.other", "ceil")
lazy_import("sage.rings.number_field.number_field", "QuadraticField")


class BSD_data:
    """
    Helper class used to keep track of information in proving BSD.

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.BSD import BSD_data
        sage: D = BSD_data()
        sage: D.Sha is None
        True
        sage: D.curve=EllipticCurve('11a')
        sage: D.update()
        sage: D.Sha
        Tate-Shafarevich group for the Elliptic Curve
         defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
    """
    def __init__(self):
        self.curve = None
        self.two_tor_rk = None
        self.Sha = None
        self.sha_an = None
        self.N = None

        self.rank = None
        self.gens = None
        self.bounds = {}  # p : (low_bd, up_bd) bounds on ord_p(#sha)
        self.primes = None  # BSD(E,p) holds for odd primes p outside this set
        self.heegner_indexes = {}  # D : I_K, K = QQ(\sqrt(D))
        self.heegner_index_upper_bound = {}  # D : M, I_K <= M
        self.N_factorization = None
        self.proof = {}

    def update(self):
        """
        Update some properties from ``curve``.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.BSD import BSD_data
            sage: D = BSD_data()
            sage: D.Sha is None
            True
            sage: D.curve = EllipticCurve('11a')
            sage: D.update()
            sage: D.Sha
            Tate-Shafarevich group for the Elliptic Curve
             defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field
        """
        self.two_tor_rk = self.curve.two_torsion_rank()
        self.Sha = self.curve.sha()
        self.sha_an = self.Sha.an(use_database=True)
        self.N = self.curve.conductor()


def simon_two_descent_work(E, two_tor_rk):
    """
    Prepare the output from Simon two-descent.

    INPUT:

    - ``E`` -- an elliptic curve

    - ``two_tor_rk`` -- its two-torsion rank

    OUTPUT:

    - a lower bound on the rank

    - an upper bound on the rank

    - a lower bound on the rank of Sha[2]

    - an upper bound on the rank of Sha[2]

    - a list of the generators found

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.BSD import simon_two_descent_work
        sage: E = EllipticCurve('14a')
        sage: simon_two_descent_work(E, E.two_torsion_rank())
        doctest:warning
        ...
        DeprecationWarning: Use E.rank(algorithm="pari") instead, as this script has been ported over to pari.
        See https://github.com/sagemath/sage/issues/35621 for details.
        (0, 0, 0, 0, [])
        sage: E = EllipticCurve('37a')
        sage: simon_two_descent_work(E, E.two_torsion_rank())
        (1, 1, 0, 0, [(0 : 0 : 1)])
    """
    from sage.misc.superseded import deprecation
    deprecation(35621, 'Use the two-descent in pari instead, as this script has been ported over to pari.')
    rank_lower_bd, two_sel_rk, gens = E.simon_two_descent()
    rank_upper_bd = two_sel_rk - two_tor_rk
    gens = [P for P in gens if P.additive_order() == Infinity]
    return rank_lower_bd, rank_upper_bd, 0, rank_upper_bd - rank_lower_bd, gens


def mwrank_two_descent_work(E, two_tor_rk):
    """
    Prepare the output from mwrank two-descent.

    INPUT:

    - ``E`` -- an elliptic curve

    - ``two_tor_rk`` -- its two-torsion rank

    OUTPUT:

    - a lower bound on the rank

    - an upper bound on the rank

    - a lower bound on the rank of Sha[2]

    - an upper bound on the rank of Sha[2]

    - a list of the generators found

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.BSD import mwrank_two_descent_work
        sage: E = EllipticCurve('14a')
        sage: mwrank_two_descent_work(E, E.two_torsion_rank())
        (0, 0, 0, 0, [])
        sage: E = EllipticCurve('37a')
        sage: mwrank_two_descent_work(E, E.two_torsion_rank())
        (1, 1, 0, 0, [(0 : -1 : 1)])
    """
    MWRC = E.mwrank_curve()
    rank_upper_bd = MWRC.rank_bound()
    rank_lower_bd = MWRC.rank()
    gens = [E(P) for P in MWRC.gens()]
    sha2_lower_bd = MWRC.selmer_rank() - two_tor_rk - rank_upper_bd
    sha2_upper_bd = MWRC.selmer_rank() - two_tor_rk - rank_lower_bd
    return rank_lower_bd, rank_upper_bd, sha2_lower_bd, sha2_upper_bd, gens


def pari_two_descent_work(E):
    r"""
    Prepare the output from pari by two-isogeny.

    INPUT:

    - ``E`` -- an elliptic curve

    OUTPUT: a tuple of 5 elements with the first 4 being integers

    - a lower bound on the rank

    - an upper bound on the rank

    - a lower bound on the rank of Sha[2]

    - an upper bound on the rank of Sha[2]

    - a list of the generators found

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.BSD import pari_two_descent_work
        sage: E = EllipticCurve('14a')
        sage: pari_two_descent_work(E)
        (0, 0, 0, 0, [])
        sage: E = EllipticCurve('37a')
        sage: pari_two_descent_work(E) # random, up to sign
        (1, 1, 0, 0, [(0 : -1 : 1)])
        sage: E = EllipticCurve('210e7')
        sage: pari_two_descent_work(E)
        (0, 2, 0, 2, [])
        sage: E = EllipticCurve('66b3')
        sage: pari_two_descent_work(E)
        (0, 0, 2, 2, [])
    """
    ep = E.pari_curve()
    lower, rank_upper_bd, s, pts = ep.ellrank()
    gens = sorted([E.point([QQ(x[0]),QQ(x[1])], check=True) for x in pts])
    gens = E.saturation(gens)[0]
    # this is explained in the pari-gp documentation:
    # s is the dimension of Sha[2]/2Sha[4],
    # which is a lower bound for dim Sha[2]
    # dim Sha[2] = dim Sel2 - rank E(Q) - dim tors
    # rank_upper_bd = dim Sel_2 - dim tors - s
    sha_upper_bd = rank_upper_bd - len(gens) + s
    return len(gens), rank_upper_bd, s, sha_upper_bd, gens


def native_two_isogeny_descent_work(E, two_tor_rk):
    """
    Prepare the output from two-descent by two-isogeny.

    INPUT:

    - ``E`` -- an elliptic curve

    - ``two_tor_rk`` -- its two-torsion rank

    OUTPUT:

    - a lower bound on the rank

    - an upper bound on the rank

    - a lower bound on the rank of Sha[2]

    - an upper bound on the rank of Sha[2]

    - a list of the generators found (currently None, since we don't store them)

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.BSD import native_two_isogeny_descent_work
        sage: E = EllipticCurve('14a')
        sage: native_two_isogeny_descent_work(E, E.two_torsion_rank())
        (0, 0, 0, 0, None)
        sage: E = EllipticCurve('65a')
        sage: native_two_isogeny_descent_work(E, E.two_torsion_rank())
        (1, 1, 0, 0, None)
    """
    from sage.schemes.elliptic_curves.descent_two_isogeny import two_descent_by_two_isogeny
    n1, n2, n1p, n2p = two_descent_by_two_isogeny(E)
    # bring n1 and n1p up to the nearest power of two
    two = ZZ(2)  # otherwise "log" is symbolic >.<
    e1 = ceil(ZZ(n1).log(two))
    e1p = ceil(ZZ(n1p).log(two))
    e2 = ZZ(n2).log(two)
    e2p = ZZ(n2p).log(two)
    rank_lower_bd = e1 + e1p - 2
    rank_upper_bd = e2 + e2p - 2
    sha_upper_bd = e2 + e2p - e1 - e1p
    gens = None  # right now, we are not keeping track of them
    return rank_lower_bd, rank_upper_bd, 0, sha_upper_bd, gens


def heegner_index_work(E):
    """
    Prepare the input and output for computing the heegner index.

    INPUT:

    - ``E`` -- an elliptic curve

    OUTPUT:

    - a Heegner index

    - the discriminant used

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.BSD import heegner_index_work
        sage: heegner_index_work(EllipticCurve('14a'))
        (1, -31)
    """
    for D in E.heegner_discriminants_list(10):
        I = None
        while I is None:
            dsl = 15
            try:
                I = E.heegner_index(D, descent_second_limit=dsl)
            except RuntimeError as err:
                if err.args[0][-33:] == 'Generators not provably computed.':
                    dsl += 1
                else:
                    raise RuntimeError(err)
        J = I.is_int()
        if J[0] and J[1] > 0:
            I = J[1]
        else:
            J = (2 * I).is_int()
            if J[0] and J[1] > 0:
                I = J[1]
            else:
                I = None
        if I is not None:
            return I, D


def prove_BSD(E, verbosity=0, two_desc='mwrank', proof=None, secs_hi=5,
              return_BSD=False):
    r"""
    Attempt to prove the Birch and Swinnerton-Dyer conjectural
    formula for `E`, returning a list of primes `p` for which this
    function fails to prove BSD(E,p).

    Here, BSD(E,p) is the
    statement: "the Birch and Swinnerton-Dyer formula holds up to a
    rational number coprime to `p`."

    INPUT:

    - ``E`` -- an elliptic curve

    - ``verbosity`` -- integer; how much information about the proof to print

      - 0: print nothing
      - 1: print sketch of proof
      - 2: print information about remaining primes

    - ``two_desc`` -- string (default: ``'mwrank'``); what to use for the
      two-descent. Options are ``'mwrank', 'pari', 'sage'``.

    - ``proof`` -- boolean or ``None`` (default: None, see
      proof.elliptic_curve or sage.structure.proof). If ``False``, this
      function just immediately returns the empty list.

    - ``secs_hi`` -- maximum number of seconds to try to compute the
      Heegner index before switching over to trying to compute the
      Heegner index bound. (Rank 0 only!)

    - ``return_BSD`` -- boolean (default: ``False``); whether to return an object
      which contains information to reconstruct a proof

    .. NOTE::

        When printing verbose output, phrases such as "by Mazur" are referring
        to the following list of papers:

    REFERENCES:

    - [Cha2005]_
    - [Jet2008]_
    - [Kat2004]_
    - [Kol1991]_
    - [LW2015]_
    - [LS]
    - [Maz1978]_
    - [Rub1991]_
    - [SW2013]_
    - [GJPST2009]_

    EXAMPLES::

        sage: EllipticCurve('11a').prove_BSD(verbosity=2)
        p = 2: True by 2-descent
        True for p not in {2, 5} by Kolyvagin.
        Kolyvagin's bound for p = 5 applies by Lawson-Wuthrich
        True for p = 5 by Kolyvagin bound
        []

        sage: EllipticCurve('14a').prove_BSD(verbosity=2)
        p = 2: True by 2-descent
        True for p not in {2, 3} by Kolyvagin.
        Kolyvagin's bound for p = 3 applies by Lawson-Wuthrich
        True for p = 3 by Kolyvagin bound
        []

        sage: E = EllipticCurve("20a1")
        sage: E.prove_BSD(verbosity=2)
        p = 2: True by 2-descent
        True for p not in {2, 3} by Kolyvagin.
        Kato further implies that #Sha[3] is trivial.
        []

        sage: E = EllipticCurve("50b1")
        sage: E.prove_BSD(verbosity=2)
        p = 2: True by 2-descent
        True for p not in {2, 3, 5} by Kolyvagin.
        Kolyvagin's bound for p = 3 applies by Lawson-Wuthrich
        Kolyvagin's bound for p = 5 applies by Lawson-Wuthrich
        True for p = 3 by Kolyvagin bound
        True for p = 5 by Kolyvagin bound
        []
        sage: E.prove_BSD(two_desc='pari')
        []

    A rank two curve::

        sage: E = EllipticCurve('389a')

    We know nothing with proof=True::

        sage: E.prove_BSD()
        Set of all prime numbers: 2, 3, 5, 7, ...

    We (think we) know everything with proof=False::

        sage: E.prove_BSD(proof=False)
        []

    A curve of rank 0 and prime conductor::

        sage: E = EllipticCurve('19a')
        sage: E.prove_BSD(verbosity=2)
        p = 2: True by 2-descent
        True for p not in {2, 3} by Kolyvagin.
        Kolyvagin's bound for p = 3 applies by Lawson-Wuthrich
        True for p = 3 by Kolyvagin bound
        []

        sage: E = EllipticCurve('37a')
        sage: E.rank()
        1
        sage: E._EllipticCurve_rational_field__rank
        (1, True)
        sage: E.analytic_rank = lambda : 0
        sage: E.prove_BSD()
        Traceback (most recent call last):
        ...
        RuntimeError: It seems that the rank conjecture does not hold for this curve
        (Elliptic Curve defined by y^2 + y = x^3 - x over Rational Field)!
        This may be a counterexample to BSD, but is more likely a bug.

    We test the consistency check for the 2-part of Sha::

        sage: E = EllipticCurve('37a')
        sage: S = E.sha(); S
        Tate-Shafarevich group for the Elliptic Curve defined by y^2 + y = x^3 - x
         over Rational Field
        sage: def foo(use_database):
        ....:  return 4
        sage: S.an = foo
        sage: E.prove_BSD()
        Traceback (most recent call last):
        ...
        RuntimeError: Apparent contradiction: 0 <= rank(sha[2]) <= 0, but ord_2(sha_an) = 2

    An example with a Tamagawa number at 5::

        sage: E = EllipticCurve('123a1')
        sage: E.prove_BSD(verbosity=2)
        p = 2: True by 2-descent
        True for p not in {2, 5} by Kolyvagin.
        Kolyvagin's bound for p = 5 applies by Lawson-Wuthrich
        True for p = 5 by Kolyvagin bound
        []

    A curve for which 3 divides the order of the Tate-Shafarevich group::

        sage: E = EllipticCurve('681b')
        sage: E.prove_BSD(verbosity=2)               # long time
        p = 2: True by 2-descent...
        True for p not in {2, 3} by Kolyvagin....
        Remaining primes:
        p = 3: irreducible, surjective, non-split multiplicative
            (0 <= ord_p <= 2)
            ord_p(#Sha_an) = 2
        [3]

    A curve for which we need to use ``heegner_index_bound``::

        sage: E = EllipticCurve('198b')
        sage: E.prove_BSD(verbosity=1, secs_hi=1)
        p = 2: True by 2-descent
        True for p not in {2, 3} by Kolyvagin.
        [3]

    The ``return_BSD`` option gives an object with detailed information
    about the proof::

        sage: E = EllipticCurve('26b')
        sage: B = E.prove_BSD(return_BSD=True)
        sage: B.two_tor_rk
        0
        sage: B.N
        26
        sage: B.gens
        []
        sage: B.primes
        []
        sage: B.heegner_indexes
        {-23: 2}

    TESTS:

    This was fixed by :issue:`8184` and :issue:`7575`::

        sage: EllipticCurve('438e1').prove_BSD(verbosity=1)
        p = 2: True by 2-descent...
        True for p not in {2} by Kolyvagin.
        []

    ::

        sage: E = EllipticCurve('960d1')
        sage: E.prove_BSD(verbosity=1)  # long time (4s on sage.math, 2011)
        p = 2: True by 2-descent
        True for p not in {2} by Kolyvagin.
        []

    ::

        sage: E = EllipticCurve('66b3')
        sage: E.prove_BSD(two_desc="pari",verbosity=1)
        p = 2: True by 2-descent
        True for p not in {2} by Kolyvagin.
        []
    """
    if proof is None:
        from sage.structure.proof.proof import get_flag
        proof = get_flag(proof, "elliptic_curve")
    else:
        proof = bool(proof)
    if not proof:
        return []
    from copy import copy
    BSD = BSD_data()
    # We replace this curve by the optimal curve, which we can do since
    # truth of BSD(E,p) is invariant under isogeny.
    BSD.curve = E.optimal_curve()
    if BSD.curve.has_cm():
        # ensure that CM is by a maximal order
        non_max_j_invs = [-12288000, 54000, 287496, 16581375]
        if BSD.curve.j_invariant() in non_max_j_invs:  # is this possible for optimal curves?
            if verbosity > 0:
                print('CM by non maximal order: switching curves')
            for E in BSD.curve.isogeny_class():
                if E.j_invariant() not in non_max_j_invs:
                    BSD.curve = E
                    break
    BSD.update()
    galrep = BSD.curve.galois_representation()

    if two_desc == 'mwrank':
        M = mwrank_two_descent_work(BSD.curve, BSD.two_tor_rk)
    elif two_desc == 'pari':
        M = pari_two_descent_work(BSD.curve)
    elif two_desc == 'sage':
        M = native_two_isogeny_descent_work(BSD.curve, BSD.two_tor_rk)
    else:
        raise NotImplementedError()
    rank_lower_bd, rank_upper_bd, sha2_lower_bd, sha2_upper_bd, gens = M
    assert sha2_lower_bd <= sha2_upper_bd
    if gens is not None:
        gens = BSD.curve.saturation(gens)[0]
    if rank_lower_bd > rank_upper_bd:
        raise RuntimeError("Apparent contradiction: %d <= rank <= %d." % (rank_lower_bd, rank_upper_bd))
    BSD.two_selmer_rank = rank_upper_bd + sha2_lower_bd + BSD.two_tor_rk
    if sha2_upper_bd == sha2_lower_bd:
        BSD.rank = rank_lower_bd
        BSD.bounds[2] = (sha2_lower_bd, sha2_upper_bd)
    else:
        BSD.rank = BSD.curve.rank(use_database=True)
        sha2_upper_bd -= (BSD.rank - rank_lower_bd)
        BSD.bounds[2] = (sha2_lower_bd, sha2_upper_bd)
        if verbosity > 0:
            print("Unable to compute the rank exactly -- used database.")
    if rank_lower_bd > 1:
        # We do not know BSD(E,p) for even a single p, since it's
        # an open problem to show that L^r(E,1)/(Reg*Omega) is
        # rational for any curve with r >= 2.
        from sage.sets.primes import Primes
        BSD.primes = Primes()
        if return_BSD:
            BSD.rank = rank_lower_bd
            return BSD
        return BSD.primes
    if (BSD.sha_an.ord(2) == 0) != (BSD.bounds[2][1] == 0):
        raise RuntimeError("Apparent contradiction: %d <= rank(sha[2]) <= %d, but ord_2(sha_an) = %d" % (sha2_lower_bd, sha2_upper_bd, BSD.sha_an.ord(2)))
    if BSD.bounds[2][0] == BSD.sha_an.ord(2) and BSD.sha_an.ord(2) == BSD.bounds[2][1]:
        if verbosity > 0:
            print('p = 2: True by 2-descent')
        BSD.primes = []
        BSD.bounds.pop(2)
        BSD.proof[2] = ['2-descent']
    else:
        BSD.primes = [2]
        BSD.proof[2] = [('2-descent',) + BSD.bounds[2]]
    if len(gens) > rank_lower_bd or rank_lower_bd > rank_upper_bd:
        raise RuntimeError("Something went wrong with 2-descent.")
    if BSD.rank != len(gens):
        gens = BSD.curve.gens(proof=True)
        if BSD.rank != len(gens):
            raise RuntimeError("Could not get generators")
    BSD.gens = [BSD.curve.point(x, check=True) for x in gens]

    if BSD.rank != BSD.curve.analytic_rank():
        raise RuntimeError("It seems that the rank conjecture does not hold for this curve (%s)! This may be a counterexample to BSD, but is more likely a bug." % BSD.curve)

    # reduce set of remaining primes to a finite set
    kolyvagin_primes = []
    heegner_index = None
    if BSD.rank == 0:
        for D in BSD.curve.heegner_discriminants_list(10):
            max_height = max(13, BSD.curve.quadratic_twist(D).CPS_height_bound())
            heegner_primes = -1
            while heegner_primes == -1:
                if max_height > 21:
                    break
                heegner_primes, _, exact = BSD.curve.heegner_index_bound(D, max_height=max_height)
                max_height += 1
            if isinstance(heegner_primes, list):
                break
        if not isinstance(heegner_primes, list):
            raise RuntimeError("Tried 10 Heegner discriminants, and heegner_index_bound failed each time.")
        if exact is not False:
            heegner_index = exact
            BSD.heegner_indexes[D] = exact
        else:
            BSD.heegner_index_upper_bound[D] = max(heegner_primes + [1])
        if 2 in heegner_primes:
            heegner_primes.remove(2)
    else:  # rank 1
        for D in BSD.curve.heegner_discriminants_list(10):
            I = BSD.curve.heegner_index(D)
            J = I.is_int()
            if J[0] and J[1] > 0:
                I = J[1]
            else:
                J = (2 * I).is_int()
                if J[0] and J[1] > 0:
                    I = J[1]
                else:
                    continue
            heegner_index = I
            BSD.heegner_indexes[D] = I
            break
        heegner_primes = [p for p in prime_divisors(heegner_index) if p != 2]

    assert BSD.sha_an in ZZ and BSD.sha_an > 0
    if BSD.curve.has_cm():
        if BSD.curve.analytic_rank() == 0:
            if verbosity > 0:
                print('    p >= 5: true by Rubin')
            BSD.primes.append(3)
        else:
            K = QuadraticField(BSD.curve.cm_discriminant(), 'a')
            D_K = K.disc()
            D_E = BSD.curve.discriminant()
            if len(K.factor(3)) == 1:  # 3 does not split in K
                BSD.primes.append(3)
            for p in prime_divisors(D_K):
                if p >= 5:
                    BSD.primes.append(p)
            for p in prime_divisors(D_E):
                if p >= 5 and D_K % p and len(K.factor(p)) == 1:
                    # p is inert in K
                    BSD.primes.append(p)
            for p in heegner_primes:
                if p >= 5 and D_E % p and D_K % p and len(K.factor(p)) == 1:
                    # p is good for E and inert in K
                    kolyvagin_primes.append(p)
            for p in prime_divisors(BSD.sha_an):
                if p >= 5 and D_K % p and len(K.factor(p)) == 1:
                    if BSD.curve.is_good(p):
                        if verbosity > 2 and p in heegner_primes and heegner_index is None:
                            print('ALERT: Prime p (%d) >= 5 dividing sha_an, good for E, inert in K, in heegner_primes, should not divide the actual Heegner index')
                        # Note that the following check is not entirely
                        # exhaustive, in case there is a p not dividing
                        # the Heegner index in heegner_primes,
                        # for which only an outer bound was computed
                        if p not in heegner_primes:
                            raise RuntimeError("p = %d divides sha_an, is of good reduction for E, inert in K, and does not divide the Heegner index. This may be a counterexample to BSD, but is more likely a bug. %s" % (p, BSD.curve))
            if verbosity > 0:
                print('True for p not in {%s} by Kolyvagin (via Stein & Lum -- unpublished) and Rubin.' % str(list(set(BSD.primes).union(set(kolyvagin_primes))))[1:-1])
        BSD.proof['finite'] = copy(BSD.primes)
    else:  # no CM
        # do some tricks to get to a finite set without calling bound_kolyvagin
        BSD.primes += [p for p in galrep.non_surjective() if p != 2]
        for p in heegner_primes:
            if p not in BSD.primes:
                BSD.primes.append(p)
        for p in prime_divisors(BSD.sha_an):
            if p not in BSD.primes and p != 2:
                BSD.primes.append(p)
        if verbosity > 0:
            s = str(BSD.primes)[1:-1]
            if 2 not in BSD.primes:
                if not s:
                    s = '2'
                else:
                    s = '2, ' + s
            print('True for p not in {' + s + '} by Kolyvagin.')
        BSD.proof['finite'] = copy(BSD.primes)
        primes_to_remove = []
        for p in BSD.primes:
            if p == 2:
                continue
            if galrep.is_surjective(p) and not BSD.curve.has_additive_reduction(p):
                if BSD.curve.has_nonsplit_multiplicative_reduction(p):
                    if BSD.rank > 0:
                        continue
                if p == 3:
                    if (not (BSD.curve.is_ordinary(p) and BSD.curve.is_good(p))) and (not BSD.curve.has_split_multiplicative_reduction(p)):
                        continue
                    if BSD.rank > 0:
                        continue
                if verbosity > 1:
                    print('    p = %d: Trying p_primary_bound' % p)
                p_bound = BSD.Sha.p_primary_bound(p)
                if p in BSD.proof:
                    BSD.proof[p].append(('Stein-Wuthrich', p_bound))
                else:
                    BSD.proof[p] = [('Stein-Wuthrich', p_bound)]
                if BSD.sha_an.ord(p) == 0 and p_bound == 0:
                    if verbosity > 0:
                        print('True for p=%d by Stein-Wuthrich.' % p)
                    primes_to_remove.append(p)
                else:
                    if p in BSD.bounds:
                        BSD.bounds[p][1] = min(BSD.bounds[p][1], p_bound)
                    else:
                        BSD.bounds[p] = (0, p_bound)
                    print('Analytic %d-rank is ' % p + str(BSD.sha_an.ord(p)) + ', actual %d-rank is at most %d.' % (p, p_bound))
                    print('    by Stein-Wuthrich.\n')
        for p in primes_to_remove:
            BSD.primes.remove(p)
        kolyvagin_primes = []
        for p in BSD.primes:
            if p == 2:
                continue
            if galrep.is_surjective(p):
                kolyvagin_primes.append(p)
        for p in kolyvagin_primes:
            BSD.primes.remove(p)
    # apply other hypotheses which imply Kolyvagin's bound holds
    D_K = QuadraticField(D, 'a').disc()

    # Cha's hypothesis
    for p in BSD.primes:
        if p == 2:
            continue
        if D_K % p != 0 and BSD.N % (p**2) != 0 and galrep.is_irreducible(p):
            if verbosity > 0:
                print('Kolyvagin\'s bound for p = %d applies by Cha.' % p)
            if p in BSD.proof:
                BSD.proof[p].append('Cha')
            else:
                BSD.proof[p] = ['Cha']
            kolyvagin_primes.append(p)
    # Stein et al replaced
    for p in BSD.primes:
        # the lemma about the vanishing of H^1 is false in Stein et al for p=5 and 11
        # here is the correction from Lawson-Wuthrich. Especially Theorem 14 in
        # [LW2015] above.
        if p in kolyvagin_primes or p == 2 or D_K % p == 0:
            continue
        crit_lw = False
        if p > 11 or p == 7:
            crit_lw = True
        elif p == 11:
            if BSD.N != 121 or BSD.curve.label() != "121c2":
                crit_lw = True
        elif galrep.is_irreducible(p):
            crit_lw = True
        else:
            phis = BSD.curve.isogenies_prime_degree(p)
            if len(phis) != 1:
                crit_lw = True
            else:
                C = phis[0].codomain()
                if p == 3:
                    if BSD.curve.torsion_order() % p != 0 and C.torsion_order() % p != 0:
                        crit_lw = True
                else:  # p == 5
                    Et = BSD.curve.quadratic_twist(5)
                    if Et.torsion_order() % p != 0 and C.torsion_order() % p != 0:
                        crit_lw = True
        if crit_lw:
            if verbosity > 0:
                print('Kolyvagin\'s bound for p = %d applies by Lawson-Wuthrich' % p)
            kolyvagin_primes.append(p)
            if p in BSD.proof:
                BSD.proof[p].append('Lawson-Wuthrich')
            else:
                BSD.proof[p] = ['Lawson-Wuthrich']
    for p in kolyvagin_primes:
        if p in BSD.primes:
            BSD.primes.remove(p)

    # apply Kolyvagin's bound
    primes_to_remove = []
    for p in kolyvagin_primes:
        if p == 2:
            continue
        if p not in heegner_primes:
            ord_p_bound = 0
        elif heegner_index is not None:  # p must divide heegner_index
            ord_p_bound = 2 * heegner_index.ord(p)
            # Here Jetchev's results apply.
            m_max = max([BSD.curve.tamagawa_number(q).ord(p) for q in BSD.N.prime_divisors()])
            if m_max > 0:
                if verbosity > 0:
                    print('Jetchev\'s results apply (at p = %d) with m_max =' % p, m_max)
                if p in BSD.proof:
                    BSD.proof[p].append(('Jetchev', m_max))
                else:
                    BSD.proof[p] = [('Jetchev', m_max)]
            ord_p_bound -= 2 * m_max
        else:  # Heegner index is None
            for D in BSD.heegner_index_upper_bound:
                M = BSD.heegner_index_upper_bound[D]
                ord_p_bound = 0
                while p**(ord_p_bound + 1) <= M**2:
                    ord_p_bound += 1
                # now ord_p_bound is one on I_K!!!
                ord_p_bound *= 2  # by Kolyvagin, now ord_p_bound is one on #Sha
                break
        if p in BSD.proof:
            BSD.proof[p].append(('Kolyvagin', ord_p_bound))
        else:
            BSD.proof[p] = [('Kolyvagin', ord_p_bound)]
        if BSD.sha_an.ord(p) == 0 and ord_p_bound == 0:
            if verbosity > 0:
                print('True for p = %d by Kolyvagin bound' % p)
            primes_to_remove.append(p)
        elif BSD.sha_an.ord(p) > ord_p_bound:
            raise RuntimeError("p = %d: ord_p_bound == %d, but sha_an.ord(p) == %d. This appears to be a counterexample to BSD, but is more likely a bug." % (p, ord_p_bound, BSD.sha_an.ord(p)))
        else:  # BSD.sha_an.ord(p) <= ord_p_bound != 0:
            if p in BSD.bounds:
                low = BSD.bounds[p][0]
                BSD.bounds[p] = (low, min(BSD.bounds[p][1], ord_p_bound))
            else:
                BSD.bounds[p] = (0, ord_p_bound)
    for p in primes_to_remove:
        kolyvagin_primes.remove(p)
    BSD.primes = list(set(BSD.primes).union(set(kolyvagin_primes)))

    # Kato's bound
    if BSD.rank == 0 and not BSD.curve.has_cm():
        L_over_Omega = BSD.curve.lseries().L_ratio()
        kato_primes = BSD.Sha.bound_kato()
        primes_to_remove = []
        for p in BSD.primes:
            if p == 2:
                continue
            if p not in kato_primes:
                if verbosity > 0:
                    print('Kato further implies that #Sha[%d] is trivial.' % p)
                primes_to_remove.append(p)
                if p in BSD.proof:
                    BSD.proof[p].append(('Kato', 0))
                else:
                    BSD.proof[p] = [('Kato', 0)]
            if p not in [2, 3] and BSD.N % p != 0:
                if galrep.is_surjective(p):
                    bd = L_over_Omega.valuation(p)
                    if verbosity > 1:
                        print('Kato implies that ord_p(#Sha[%d]) <= %d ' % (p, bd))
                    if p in BSD.proof:
                        BSD.proof[p].append(('Kato', bd))
                    else:
                        BSD.proof[p] = [('Kato', bd)]
                    if p in BSD.bounds:
                        low = BSD.bounds[p][0]
                        BSD.bounds[p][1] = (low, min(BSD.bounds[p][1], bd))
                    else:
                        BSD.bounds[p] = (0, bd)
        for p in primes_to_remove:
            BSD.primes.remove(p)

    # Mazur
    primes_to_remove = []
    if BSD.N.is_prime():
        for p in BSD.primes:
            if p == 2:
                continue
            if galrep.is_reducible(p):
                primes_to_remove.append(p)
                if verbosity > 0:
                    print('True for p=%s by Mazur' % p)
        for p in primes_to_remove:
            BSD.primes.remove(p)
            if p in BSD.proof:
                BSD.proof[p].append('Mazur')
            else:
                BSD.proof[p] = ['Mazur']

    BSD.primes.sort()

    # Try harder to compute the Heegner index, where it matters
    if heegner_index is None:
        max_height = max(max_height, 18)
        for D in BSD.heegner_index_upper_bound:
            M = BSD.heegner_index_upper_bound[D]
            for p in kolyvagin_primes:
                if p not in BSD.primes or p == 3:
                    continue
                if verbosity > 0:
                    print('    p = %d: Trying harder for Heegner index' % p)
                obt = 0
                while p**(BSD.sha_an.ord(p) / 2 + 1) <= M and max_height < 22:
                    if verbosity > 2:
                        print('    trying max_height =', max_height)
                    old_bound = M
                    M, _, exact = BSD.curve.heegner_index_bound(D, max_height=max_height, secs_dc=secs_hi)
                    if M == -1:
                        max_height += 1
                        continue
                    if exact is not False:
                        heegner_index = exact
                        BSD.heegner_indexes[D] = exact
                        M = exact
                        if verbosity > 2:
                            print('    heegner index =', M)
                    else:
                        M = max(M + [1])
                        if verbosity > 2:
                            print('    bound =', M)
                    if old_bound == M:
                        obt += 1
                        if obt == 2:
                            break
                    max_height += 1
                BSD.heegner_index_upper_bound[D] = min(M, BSD.heegner_index_upper_bound[D])
                low, upp = BSD.bounds[p]
                expn = 0
                while p**(expn + 1) <= M:
                    expn += 1
                if 2 * expn < upp:
                    upp = 2 * expn
                    BSD.bounds[p] = (low, upp)
                    if verbosity > 0:
                        print('    got better bound on ord_p =', upp)
                    if low == upp:
                        if upp != BSD.sha_an.ord(p):
                            raise RuntimeError
                        else:
                            if verbosity > 0:
                                print('    proven!')
                            BSD.primes.remove(p)
            break
        for p in kolyvagin_primes:
            if p not in BSD.primes or p == 3:
                continue
            for D in BSD.curve.heegner_discriminants_list(4):
                if D in BSD.heegner_index_upper_bound:
                    continue
                print('    discriminant', D)
                if verbosity > 0:
                    print('p = %d: Trying discriminant = %d for Heegner index' % (p, D))
                max_height = max(10, BSD.curve.quadratic_twist(D).CPS_height_bound())
                obt = 0
                while True:
                    if verbosity > 2:
                        print('    trying max_height =', max_height)
                    old_bound = M
                    if p**(BSD.sha_an.ord(p) / 2 + 1) > M or max_height >= 22:
                        break
                    M, _, exact = BSD.curve.heegner_index_bound(D, max_height=max_height, secs_dc=secs_hi)
                    if M == -1:
                        max_height += 1
                        continue
                    if exact is not False:
                        heegner_index = exact
                        BSD.heegner_indexes[D] = exact
                        M = exact
                        if verbosity > 2:
                            print('    heegner index =', M)
                    else:
                        M = max(M + [1])
                        if verbosity > 2:
                            print('    bound =', M)
                    if old_bound == M:
                        obt += 1
                        if obt == 2:
                            break
                    max_height += 1
                BSD.heegner_index_upper_bound[D] = M
                low, upp = BSD.bounds[p]
                expn = 0
                while p**(expn + 1) <= M:
                    expn += 1
                if 2 * expn < upp:
                    upp = 2 * expn
                    BSD.bounds[p] = (low, upp)
                    if verbosity > 0:
                        print('    got better bound =', upp)
                    if low == upp:
                        if upp != BSD.sha_an.ord(p):
                            raise RuntimeError
                        else:
                            if verbosity > 0:
                                print('    proven!')
                            BSD.primes.remove(p)
                            break

    # some extra information
    if verbosity > 1:
        if BSD.primes:
            print('Remaining primes:')
        for p in BSD.primes:
            s = 'p = ' + str(p) + ': '
            if galrep.is_irreducible(p):
                s += 'ir'
            s += 'reducible, '
            if not galrep.is_surjective(p):
                s += 'not '
            s += 'surjective, '
            a_p = BSD.curve.an(p)
            if BSD.curve.is_good(p):
                if a_p % p != 0:
                    s += 'good ordinary'
                else:
                    s += 'good, non-ordinary'
            else:
                assert BSD.curve.is_minimal()
                if a_p == 0:
                    s += 'additive'
                elif a_p == 1:
                    s += 'split multiplicative'
                elif a_p == -1:
                    s += 'non-split multiplicative'
            if BSD.curve.tamagawa_product() % p == 0:
                s += ', divides a Tamagawa number'
            if p in BSD.bounds:
                s += '\n    (%d <= ord_p <= %d)' % BSD.bounds[p]
            else:
                s += '\n    (no bounds found)'
            s += '\n    ord_p(#Sha_an) = %d' % BSD.sha_an.ord(p)
            if heegner_index is None:
                may_divide = True
                for D in BSD.heegner_index_upper_bound:
                    if p > BSD.heegner_index_upper_bound[D] or p not in kolyvagin_primes:
                        may_divide = False
                if may_divide:
                    s += '\n    may divide the Heegner index, for which only a bound was computed'
            print(s)

    if BSD.curve.has_cm():
        if BSD.rank == 1:
            BSD.proof['reason_finite'] = 'Rubin&Kolyvagin'
        else:
            BSD.proof['reason_finite'] = 'Rubin'
    else:
        BSD.proof['reason_finite'] = 'Kolyvagin'
    # reduce memory footprint of BSD object:
    BSD.curve = BSD.curve.label()
    BSD.Sha = None
    return BSD if return_BSD else BSD.primes
