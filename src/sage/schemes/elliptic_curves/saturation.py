# sage.doctest: needs sage.rings.finite_rings sage.rings.number_field
r"""
Saturation of Mordell-Weil groups of elliptic curves over number fields

Points `P_1`, `\dots`, `P_r` in `E(K)`, where `E` is an elliptic curve
over a number field `K`, are said to be `p`-saturated if no linear
combination `\sum n_iP_i` is divisible by `p` in `E(K)` except
trivially when all `n_i` are multiples of `p`.  The points are said to
be saturated if they are `p`-saturated at all primes; this is always
true for all but finitely many primes since `E(K)` is a
finitely-generated Abelian group.

The process of `p`-saturating a given set of points is implemented
here.  The naive algorithm simply checks all `(p^r-1)/(p-1)`
projective combinations of the points, testing each to see if it can
be divided by `p`.  If this occurs then we replace one of the points
and continue.  The function :meth:`p_saturation` does one step of
this, while :meth:`full_p_saturation` repeats until the points are
`p`-saturated.  A more sophisticated algorithm for `p`-saturation is
implemented which is much more efficient for large `p` and `r`, and
involves computing the reduction of the points modulo auxiliary primes
to obtain linear conditions modulo `p` which must be satisfied by the
coefficients `a_i` of any nontrivial relation.  When the points are
already `p`-saturated this sieving technique can prove their
saturation quickly.

The method :meth:`saturation` of the class :class:`EllipticCurve_number_field`
applies full `p`-saturation at any given set of primes, or can compute
a bound on the primes `p` at which the given points may not be
`p`-saturated.  This involves computing a lower bound for the
canonical height of points of infinite order, together with estimates
from the geometry of numbers.

AUTHORS:

- Robert Bradshaw

- John Cremona
"""
#*****************************************************************************
#       Copyright (C) 2017 Robert Bradshaw <robertwb@math.washington.edu>
#                          John Cremona <john.cremona@gmail.com>
#                          William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.integer_ring import ZZ
from sage.arith.misc import kronecker as kro
from sage.structure.sage_object import SageObject


def reduce_mod_q(x, amodq):
    r"""The reduction of ``x`` modulo the prime ideal defined by ``amodq``.

    INPUT:

    - ``x`` -- an element of a  number field `K`

    - ``amodq`` -- an element of `GF(q)` which is a root mod `q` of
      the defining polynomial of `K`.  This defines a degree 1 prime
      ideal `Q=(q,\alpha-a)` of `K=\QQ(\alpha)`, where `a \bmod q` =
      ``amodq``.

    OUTPUT: the image of ``x`` in the residue field of `K` at the prime `Q`

    EXAMPLES::

        sage: from sage.schemes.elliptic_curves.saturation import reduce_mod_q
        sage: x = polygen(QQ)
        sage: pol = x^3 - x^2 - 3*x + 1
        sage: K.<a> = NumberField(pol)
        sage: [(q, [(amodq, reduce_mod_q(1 - a + a^4, amodq))
        ....:       for amodq in sorted(pol.roots(GF(q), multiplicities=False))])
        ....:  for q in primes(50, 70)]
        [(53, []),
         (59, [(36, 28)]),
         (61, [(40, 35)]),
         (67, [(10, 8), (62, 28), (63, 60)])]
    """
    Fq = amodq.parent()
    try:
        return x.lift().change_ring(Fq)(amodq)
    except AttributeError: # in case x is in QQ
        return Fq(x)


class EllipticCurveSaturator(SageObject):
    r"""
    Class for saturating points on an elliptic curve over a number field.

    INPUT:

    - ``E`` -- an elliptic curve defined over a number field, or `\QQ`

    - ``verbose`` -- boolean (default: ``False``); verbosity flag

    .. NOTE::

        This function is not normally called directly by users, who
        may access the data via methods of the EllipticCurve
        classes.
    """
    def __init__(self, E, verbose=False):
        r"""
        Initialize the saturator.

        INPUT:

        - ``E`` -- an elliptic curve defined over a number field

        - ``verbose`` -- boolean (default: ``False``); verbosity flag
        """
        self._verbose = verbose
        self._curve = E
        self._N = E.discriminant().norm()
        self._field = K = E.base_field()
        if K.absolute_degree() == 1:
            from sage.rings.rational_field import QQ
            from sage.rings.polynomial.polynomial_ring import polygen
            self._Kpol = polygen(QQ)
        else:
            self._Kpol = K.defining_polynomial()
        self._D = self._Kpol.discriminant()
        self._reductions = {}
        self._lincombs = {}
        self._torsion_gens = [t.element() for t in E.torsion_subgroup().gens()]
        self._reductions = {}
        # This will hold a dictionary with keys (q,aq) with q prime
        # and aq a root of K's defining polynomial mod q, and values
        # (n,gens) where n is the cardinality of the reduction of E
        # and gens are generators of that reduction.

    def add_reductions(self, q):
        r"""
        Add reduction data at primes above ``q`` if not already there.

        INPUT:

        - ``q`` -- a prime number not dividing the defining polynomial
          of ``self.__field``

        OUTPUT:

        Returns nothing, but updates self._reductions dictionary for
        key ``q`` to a dict whose keys are the roots of the defining
        polynomial mod ``q`` and values tuples (``nq``, ``Eq``) where
        ``Eq`` is an elliptic curve over `GF(q)` and ``nq`` its
        cardinality.  If ``q`` divides the conductor norm or order
        discriminant nothing is added.

        EXAMPLES:

        Over `\QQ`::

            sage: from sage.schemes.elliptic_curves.saturation import EllipticCurveSaturator
            sage: E = EllipticCurve('11a1')
            sage: saturator = EllipticCurveSaturator(E)
            sage: saturator._reductions
            {}
            sage: saturator.add_reductions(19)
            sage: saturator._reductions
            {19: {0: (20, Elliptic Curve defined by y^2 + y = x^3 + 18*x^2 + 9*x + 18
                           over Finite Field of size 19)}}

        Over a number field::

            sage: x = polygen(QQ);  K.<a> = NumberField(x^2 + 2)
            sage: E = EllipticCurve(K, [0,1,0,a,a])
            sage: from sage.schemes.elliptic_curves.saturation import EllipticCurveSaturator
            sage: saturator = EllipticCurveSaturator(E)
            sage: for q in primes(20):
            ....:     saturator.add_reductions(q)
            sage: saturator._reductions
            {2:  {},
             3:  {},
             5:  {},
             7:  {},
             11: {3:  (16, Elliptic Curve defined by y^2 = x^3 + x^2 + 3*x + 3
                            over Finite Field of size 11),
                  8:  (8,  Elliptic Curve defined by y^2 = x^3 + x^2 + 8*x + 8
                            over Finite Field of size 11)},
             13: {},
             17: {7:  (20, Elliptic Curve defined by y^2 = x^3 + x^2 + 7*x + 7
                            over Finite Field of size 17),
                  10: (18, Elliptic Curve defined by y^2 = x^3 + x^2 + 10*x + 10
                            over Finite Field of size 17)},
             19: {6:  (16, Elliptic Curve defined by y^2 = x^3 + x^2 + 6*x + 6
                            over Finite Field of size 19),
                  13: (12, Elliptic Curve defined by y^2 = x^3 + x^2 + 13*x + 13
                            over Finite Field of size 19)}}
        """
        if q in self._reductions:
            return
        self._reductions[q] = redmodq = {}
        if q.divides(self._N) or q.divides(self._D):
            return
        from sage.schemes.elliptic_curves.constructor import EllipticCurve
        for amodq in sorted(self._Kpol.roots(GF(q), multiplicities=False)):
            Eq = EllipticCurve([reduce_mod_q(ai, amodq) for ai in self._curve.ainvs()])
            nq = Eq.cardinality()
            redmodq[amodq] = (nq, Eq)

    def full_p_saturation(self, Plist, p):
        r"""Full `p`-saturation of ``Plist``.

        INPUT:

        - ``Plist`` -- list of independent points on one elliptic curve

        - ``p`` -- integer; a prime number

        OUTPUT:

        (``newPlist``, ``exponent``) where ``newPlist`` has the same
        length as ``Plist`` and spans the `p`-saturation of the span
        of ``Plist``, which contains that span with index
        ``p**exponent``.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.saturation import EllipticCurveSaturator
            sage: E = EllipticCurve('389a')
            sage: K.<i> = QuadraticField(-1)
            sage: EK = E.change_ring(K)
            sage: P = EK(1 + i, -1 - 2*i)
            sage: saturator = EllipticCurveSaturator(EK, verbose=True)
            sage: saturator.full_p_saturation([8*P], 2)
             --starting full 2-saturation
            Points were not 2-saturated, exponent was 3
            ([(i + 1 : -2*i - 1 : 1)], 3)

            sage: Q = EK(0, 0)
            sage: R = EK(-1, 1)
            sage: saturator = EllipticCurveSaturator(EK, verbose=False)
            sage: saturator.full_p_saturation([P, Q, R], 3)
            ([(i + 1 : -2*i - 1 : 1), (0 : 0 : 1), (-1 : 1 : 1)], 0)

        An example where the points are not 7-saturated and we gain
        index exponent 1.  Running this example with verbose=True
        would show that it uses the code for when the reduction has
        `p`-rank 2 (which occurs for the reduction modulo `(16-5i)`),
        which uses the Weil pairing::

            sage: saturator.full_p_saturation([P, Q + 3*R, Q - 4*R], 7)
            ([(i + 1 : -2*i - 1 : 1),
              (2869/676 : 154413/17576 : 1),
              (-7095/502681 : -366258864/356400829 : 1)], 1)
        """
        if not Plist:
            return Plist, ZZ.zero()

        exponent = ZZ(0)

        # To handle p-torsion, we must add any torsion generators whose
        # order is divisible by p to the list of points.  Note that it is
        # not correct to add generators of the p-torsion here, we actually
        # need generators of p-cotorsion.  If there are any of these, we
        # cannot use the supplied dict of linear combinations, so we reset
        # this.  The torsion points are removed before returning the
        # saturated list.

        verbose = self._verbose
        if verbose:
            print(" --starting full %s-saturation" % p)

        n = len(Plist)  # number of points supplied & to be returned
        Plist = Plist + [T for T in self._torsion_gens if p.divides(T.order())]
        nx = len(Plist) # number of points including relevant torsion
        extra_torsion = nx-n
        if extra_torsion:
            if verbose:
                print("Adding {} torsion generators before {}-saturation".format(extra_torsion,p))

        res = self.p_saturation(Plist, p)
        while res: # res is either False or (i, newP)
            exponent += 1
            Plist[res[0]] = res[1]
            res = self.p_saturation(Plist, p)

        if extra_torsion:
            # remove the torsion points
            if verbose:
                print("Removing the torsion generators after %s-saturation" % p)
            Plist = Plist[:n]

        if verbose:
            if exponent:
                print("Points were not %s-saturated, exponent was %s" % (p,exponent))
            else:
                print("Points were %s-saturated" % p)

        return (Plist, exponent)

    def p_saturation(self, Plist, p, sieve=True):
        r"""Checks whether the list of points is `p`-saturated.

        INPUT:

        - ``Plist`` -- list of independent points on one elliptic curve

        - ``p`` -- integer; a prime number

        - ``sieve`` -- boolean; if ``True``, use a sieve (when there are at
          least 2 points), otherwise test all combinations

        .. NOTE::

            The sieve is much more efficient when the points are
            saturated and the number of points or the prime are large.

        OUTPUT:

        Either ``False`` if the points are `p`-saturated, or (``i``,
        ``newP``) if they are not `p`-saturated, in which case after
        replacing the i'th point with ``newP``, the subgroup generated
        contains that generated by ``Plist`` with index `p`.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.saturation import EllipticCurveSaturator
            sage: E = EllipticCurve('389a')
            sage: K.<i> = QuadraticField(-1)
            sage: EK = E.change_ring(K)
            sage: P = EK(1 + i, -1 - 2*i)
            sage: saturator = EllipticCurveSaturator(EK)
            sage: saturator.p_saturation([P], 2)
            False
            sage: saturator.p_saturation([2*P], 2)
            (0, (i + 1 : -2*i - 1 : 1))

            sage: Q = EK(0, 0)
            sage: R = EK(-1, 1)
            sage: saturator.p_saturation([P, Q, R], 3)
            False

        Here we see an example where 19-saturation is proved, with the
        verbose flag set to ``True`` so that we can see what is going on::

            sage: saturator = EllipticCurveSaturator(EK, verbose=True)
            sage: saturator.p_saturation([P, Q, R], 19)
            Using sieve method to saturate...
            E has 19-torsion over Finite Field of size 197, projecting points
            --> [(15 : 168 : 1), (0 : 0 : 1), (196 : 1 : 1)]
            --rank is now 1
            E has 19-torsion over Finite Field of size 197, projecting points
            --> [(184 : 27 : 1), (0 : 0 : 1), (196 : 1 : 1)]
            --rank is now 2
            E has 19-torsion over Finite Field of size 293, projecting points
            --> [(139 : 16 : 1), (0 : 0 : 1), (292 : 1 : 1)]
             --rank is now 3
            Reached full rank: points were 19-saturated
            False

        An example where the points are not 11-saturated::

            sage: saturator = EllipticCurveSaturator(EK, verbose=False)
            sage: res = saturator.p_saturation([P + 5*Q, P - 6*Q, R], 11); res
            (0, (-5783311/14600041*i + 1396143/14600041
                 : 37679338314/55786756661*i + 3813624227/55786756661 : 1))

        That means that the 0'th point may be replaced by the displayed
        point to achieve an index gain of 11::

            sage: saturator.p_saturation([res[1], P - 6*Q, R], 11)
            False

        TESTS:

        See :issue:`27387`::

            sage: K.<a> = NumberField(x^2 - x - 26)
            sage: E = EllipticCurve([a, 1 - a, 0, 93 - 16*a, 3150 - 560*a])
            sage: P = E([65 - 35*a/3, (959*a-5377)/9])
            sage: T = E.torsion_points()[0]
            sage: from sage.schemes.elliptic_curves.saturation import EllipticCurveSaturator
            sage: saturator = EllipticCurveSaturator(E, True)
            sage: saturator.p_saturation([P, T], 2)
            Using sieve method to saturate...
            ...
            -- points were not 2-saturated, gaining index 2
            (1, (0 : 1 : 0))

        A CM example where large siecing primes are needed (LMFDB
        label 2.0.3.1-50625.1-CMb2)::

            sage: K.<a> = NumberField(x^2 - x + 1)
            sage: E = EllipticCurve(K, [0, 0, 1, -750, 7906])
            sage: E.has_rational_cm()
            True
            sage: E.cm_discriminant()
            -27
            sage: points = [E([10, -38]), E([15/49*a + 760/49, 675/343*a - 884/343])]
            sage: E.saturation(points, verbose=True) # long time (17s)
            Computing lower height bound..
            ..done: 7.168735020029907e-06
            p-saturating for primes p < 132
            ...
            Saturating at p=131
            --starting full 131-saturation
            Using sieve method to saturate...
            E has 131-torsion over Finite Field of size 617011, projecting points
            --> [(10 : 616973 : 1), (592083 : 192224 : 1)]
            --rank is now 1
            --rank is now 2
            Reached full rank: points were 131-saturated
            Points were 131-saturated
            --already 131-saturated
            ([(10 : -38 : 1), (15/49*a + 760/49 : 675/343*a - 884/343 : 1)],
            1,
            0.123378097374749)
        """
        verbose = self._verbose
        # This code does a lot of elliptic curve group structure
        # computations, and would benefit immensely if those were sped
        # up.
        n = len(Plist)
        if n == 0:
            return False

        if n == 1 and p == 2:
            pts = Plist[0].division_points(p)
            if pts:
                return (0, pts[0])
            else:
                return False

        E = self._curve

        if not sieve:
            from sage.groups.generic import multiples
            from sage.schemes.projective.projective_space import ProjectiveSpace

            mults = [list(multiples(P, p)) for P in Plist[:-1]] + [list(multiples(Plist[-1],2))]
            E0 = E(0)

            for v in ProjectiveSpace(GF(p),n-1): # an iterator
                w = tuple(int(x) for x in v)
                P = sum([m[c] for m,c in zip(mults,w)],E0)
                pts = P.division_points(p)
                if pts:
                    if verbose:
                        print("  points not saturated at {}, increasing index by {}".format(p,p))
                        # w will certainly have a coordinate equal to 1
                    return (w.index(1), pts[0])
            # we only get here if no linear combination is divisible by p,
            # so the points are p-saturated:
            return False

        # Now we use the more sophisticated sieve to either prove
        # p-saturation, or compute a much smaller list of possible points
        # to test for p-divisibility:

        if verbose:
            print("Using sieve method to saturate...")
        from sage.matrix.constructor import matrix
        from sage.sets.primes import Primes

        A = matrix(GF(p), 0, n)
        rankA = 0
        count = 0

        # We reduce the curve modulo primes Q of degree 1 above
        # rational primes q chosen as follows: we avoid primes of bad
        # reduction and primes dividing the discriminant of the
        # defining polynomial of the field, so that we can avoid
        # constructing number field primes entirely and just look for
        # roots amodq of the defining polynomial mod q, then reduction
        # of the curve is easy.

        # We also avoid primes dividing the denominators of any of the
        # points: the points would reduce to 0 modulo such primes
        # anyway, and this way reduction of the points is also easy.

        # Lastly, there is a special case when E has rational CM by
        # some discriminant D.  In the non-split case (kro(D,p)=-1)
        # the image of the mod-p Galois representation is cyclic of
        # order p^2-1 and there will be no p-torsion mod Q unless
        # Frob(Q) has trivial image; a necessary condition for this is
        # q=1 (mod p) so we restrict to such q.  That way the density
        # of q for which there is a point of order p is 1/(p+1)
        # instead of 1/(p^2-1).  This test was put in after running
        # many tests: for example, LMFDB curve 2.0.3.1-50625.1-CMb2
        # has saturation index bound 132 and to saturate at p=131
        # requires q=617011. (In the split case the density is 1/(p-1)
        # and there is no simple test.)

        avoid = [self._N, self._D] + [P[0].denominator_ideal().norm() for P in Plist]
        cm_test = E.has_rational_cm() and kro(E.cm_discriminant(), p) == -1
        for q in Primes():
            if any(q.divides(m) for m in avoid):
                continue
            if cm_test and not p.divides(q-1):
                continue
            self.add_reductions(q) # does nothing if key q is already there
            for amodq in self._reductions[q]:
                (nq, Eq) = self._reductions[q][amodq]
                if not p.divides(nq):
                    continue
                if verbose:
                    print("E has %s-torsion over %s, projecting points" % (p,GF(q)))
                projPlist = [Eq([reduce_mod_q(c, amodq) for c in pt]) for pt in Plist]
                if verbose:
                    print(" --> %s" % projPlist)
                try:
                    vecs = p_projections(Eq, projPlist, p)
                except ValueError:
                    vecs = []
                for v in vecs:
                    A = matrix(A.rows()+[v])
                    newrank = A.rank()
                    if verbose:
                        print(" --rank is now %s" % newrank)
                    if newrank == n:
                        if verbose:
                            print("Reached full rank: points were %s-saturated" % p)
                        return False
                    if newrank == rankA:
                        count += 1
                        if count == 15:
                            if verbose:
                                print("! rank same for 15 steps, checking kernel...")
                            # no increase in rank for the last 15 primes Q
                            # find the points in the kernel and call the no-sieve version
                            vecs = A.right_kernel().basis()
                            if verbose:
                                print("kernel vectors: %s" % vecs)
                            Rlist = [sum([int(vi)*Pi for vi,Pi in zip(v,Plist)],E(0))
                                     for v in vecs]
                            if verbose:
                                print("points generating kernel: %s" % Rlist)

                            # If the nullity of A were 1 (the usual
                            # case) we take any nonzero vector in its
                            # kernel and compute that linear
                            # combination of the points, giving a
                            # point which is certainly a p-multiple
                            # modulo 15 primes Q, and we test if it
                            # actually is a p-multiple:
                            if len(Rlist) == 1:
                                R = Rlist[0]
                                pts = R.division_points(p)
                                if pts:
                                    pt = pts[0]
                                    v = vecs[0]
                                    # replace one of the original
                                    # points with this one; we can
                                    # replace any for which the
                                    # coefficient of v is nonzero
                                    if verbose:
                                        print("-- points were not {}-saturated, gaining index {}".format(p,p))
                                    j = next(i for i,x in enumerate(v) if x)
                                    return (j, pt)
                                else: # R is not a p-multiple so the
                                      # points were p-saturated
                                    return False

                            # Else we call the non-sieve version with
                            # a list of points which are all
                            # p-multiples modulo 15 primes, and we
                            # will just try to divide all linear
                            # combinations of them
                            res = self.p_saturation(Rlist, p, sieve=False)
                            if res:
                                jj, R = res
                                v = vecs[jj]
                                # The v-linear combination of Rlist
                                # really is p*R.  Now to enlarge the
                                # span, we may replce the j'th point
                                # in Plist with R, where v[j] is
                                # nonzero.
                                if verbose:
                                    print("-- points were not {}-saturated, gaining index {}".format(p,p))
                                j = next(i for i,x in enumerate(v) if x)
                                return (j, R)
                            else:
                                # points really were saturated
                                if verbose:
                                    print("-- points were %s-saturated" % p)
                                return False
                    else: # rank went up but is <n; carry on using more Qs
                        rankA = newrank
                        count = 0


def p_projections(Eq, Plist, p, debug=False):
    r"""

    INPUT:

    - ``Eq`` -- an elliptic curve over a finite field

    - ``Plist`` -- list of points on `Eq`

    - ``p`` -- a prime number

    OUTPUT:

    A list of `r\le2` vectors in `\GF{p^n}`, the images of the points in
    `G \otimes \GF{p}`, where `r` is the number of vectors is the
    `p`-rank of `Eq`.

    ALGORITHM:

    First project onto the `p`-primary part of `Eq`.  If that has
    `p`-rank 1 (i.e. is cyclic), use discrete logs there to define a
    map to `\GF{p}`, otherwise use the Weil pairing to define two
    independent maps to `\GF{p}`.

    EXAMPLES:

    This curve has three independent rational points::

        sage: E = EllipticCurve([0,0,1,-7,6])

    We reduce modulo `409` where its order is `3^2\cdot7^2`; the
    `3`-primary part is non-cyclic while the `7`-primary part is
    cyclic of order `49`::

        sage: F = GF(409)
        sage: EF = E.change_ring(F)
        sage: G = EF.abelian_group()
        sage: G
        Additive abelian group isomorphic to Z/147 + Z/3
         embedded in Abelian group of points on Elliptic Curve
          defined by y^2 + y = x^3 + 402*x + 6 over Finite Field of size 409
        sage: G.order().factor()
        3^2 * 7^2

    We construct three points and project them to the `p`-primary
    parts for `p=2,3,5,7`, yielding 0,2,0,1 vectors of length 3 modulo
    `p` respectively.  The exact vectors output depend on the computed
    generators of `G`::

        sage: Plist = [EF([-2,3]), EF([0,2]), EF([1,0])]
        sage: from sage.schemes.elliptic_curves.saturation import p_projections
        sage: [(p, p_projections(EF, Plist, p)) for p in primes(11)]  # random
        [(2, []), (3, [(0, 2, 2), (2, 2, 1)]), (5, []), (7, [(5, 1, 1)])]
        sage: [(p, len(p_projections(EF, Plist, p))) for p in primes(11)]
        [(2, 0), (3, 2), (5, 0), (7, 1)]
    """
    if debug:
        print("In p_projections(Eq,Plist,p) with Eq = {}, Plist = {}, p = {}".format(Eq,Plist,p))
    n = Eq.cardinality()
    m = n.prime_to_m_part(p)      # prime-to-p part of order
    if debug:
        print("m={}, n={}".format(m,n))
    if m == n: # p-primary part trivial, nothing to do
        return []
    G = Eq.abelian_group()
    if debug:
        print("gens = {}".format(G.gens()))

    # project onto p-primary part

    pts = [m*pt for pt in Plist]
    gens = [m*g.element() for g in G.gens()]
    gens = [g for g in gens if g]
    if debug:
        print("gens for {}-primary part of G: {}".format(p, gens))
        print("{}*points: {}".format(m,pts))
    from sage.groups.generic import discrete_log as dlog
    from sage.modules.free_module_element import vector
    Fp = GF(p)

    # If the p-primary part is cyclic we use elliptic discrete logs directly:

    if len(gens) == 1:
        g = gens[0]
        pp = g.order()
        if debug:
            print("Cyclic case, taking dlogs to base {} of order {}".format(g,pp))
        # logs are well-defined mod pp, hence mod p
        v = [dlog(pt, g, ord=pp, operation='+') for pt in pts]
        if debug:
            print("dlogs: {}".format(v))
        return [vector(Fp, v)]

    # We make no assumption about which generator order divides the
    # other, since conventions differ!

    orders = [g.order() for g in gens]
    p1, p2 = min(orders), max(orders)
    g1, g2 = gens
    if debug:
        print("Non-cyclic case, orders = {}, p1={}, p2={}, g1={}, g2={}".format(orders,p1,p2,g1,g2))

    # Now the p-primary part of the reduction is non-cyclic of
    # exponent p2, and we use the Weil pairing, whose values are p1'th
    # roots of unity with p1|p2, together with discrete log in the
    # multiplicative group.

    zeta = g1.weil_pairing(g2, p2)  # a primitive p1'th root of unity
    if debug:
        print("wp of gens = {} with order {}".format(zeta, zeta.multiplicative_order()))
        assert zeta.multiplicative_order() == p1, "Weil pairing error during saturation: p={}, G={}, Plist={}".format(p, G, Plist)

    # logs are well-defined mod p1, hence mod p

    return [vector(Fp, [dlog(pt.weil_pairing(g1, p2), zeta,
                             ord=p1, operation='*') for pt in pts]),
            vector(Fp, [dlog(pt.weil_pairing(g2, p2), zeta,
                             ord=p1, operation='*') for pt in pts])]
