# -*- coding: utf-8 -*-
r"""
Manin Relations for overconvergent modular symbols

Code to create the Manin Relations class, which solves the "Manin relations".
That is, a description of `Div^0(P^1(\QQ))` as a `\ZZ[\Gamma_0(N)]`-module in
terms of generators and relations is found. The method used is geometric,
constructing a nice fundamental domain for `\Gamma_0(N)` and reading the
relevant Manin relations off of that picture. The algorithm follows [PS2011]_.

AUTHORS:

- Robert Pollack, Jonathan Hanke (2012): initial version

"""
#*****************************************************************************
#       Copyright (C) 2012 Robert Pollack <rpollack@math.bu.edu>
#                          Jonathan Hanke <jonhanke@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.matrix.matrix_space import MatrixSpace
from sage.modular.modsym.all import P1List
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.structure.sage_object import SageObject
from sage.misc.cachefunc import cached_method

from .sigma0 import Sigma0

M2ZSpace = MatrixSpace(ZZ,2)


def M2Z(x):
    r"""
    Create an immutable `2 \times 2` integer matrix from ``x``.

    INPUT: anything that can be converted into a `2 \times 2` matrix.

    EXAMPLES::

        sage: from sage.modular.pollack_stevens.fund_domain import M2Z
        sage: M2Z([1,2,3,4])
        [1 2]
        [3 4]
        sage: M2Z(1)
        [1 0]
        [0 1]
    """
    x = M2ZSpace(x)
    x.set_immutable()
    return x

Id = M2Z([1, 0, 0, 1])
sig = M2Z([0, 1, -1, 0])
tau = M2Z([0, -1, 1, -1])
minone_inf_path = M2Z([1, 1, -1, 0])

# We store these so that we do not have to constantly create them.
t00 = (0, 0)
t10 = (1, 0)
t01 = (0, 1)
t11 = (1, 1)


class PollackStevensModularDomain(SageObject):
    r"""
    The domain of a modular symbol.

    INPUT:

    - ``N`` -- a positive integer, the level of the congruence subgroup
      `\Gamma_0(N)`

    - ``reps`` -- a list of `2 \times 2` matrices, the coset
      representatives of `Div^0(P^1(\QQ))`

    - ``indices`` -- a list of integers; indices of elements in
      ``reps`` which are generators

    - ``rels`` -- a list of list of triples ``(d, A, i)``, one for each
      coset representative of ``reps`` which describes how to express the
      elements of ``reps`` in terms of generators specified by ``indices``.
      See :meth:`relations` for a detailed explanations of these triples.

    - ``equiv_ind`` -- a dictionary which maps normalized coordinates on
      `P^1(\ZZ/N\ZZ)` to an integer such that a matrix whose bottom row is
      equivalent to `[a:b]` in `P^1(\ZZ/N\ZZ)` is in the coset of
      ``reps[equiv_ind[(a,b)]]``

    EXAMPLES::

        sage: from sage.modular.pollack_stevens.fund_domain import PollackStevensModularDomain, M2Z
        sage: PollackStevensModularDomain(2 , [M2Z([1,0,0,1]), M2Z([1,1,-1,0]), M2Z([0,-1,1,1])], [0,2], [[(1, M2Z([1,0,0,1]), 0)], [(-1,M2Z([-1,-1,0,-1]),0)], [(1, M2Z([1,0,0,1]), 2)]], {(0,1): 0, (1,0): 1, (1,1): 2})
        Modular Symbol domain of level 2

    TESTS:

    The level ``N`` must be an integer::

        sage: PollackStevensModularDomain(1/2, None, None, None, None)
        Traceback (most recent call last):
        ...
        TypeError: no conversion of this rational to integer
        sage: PollackStevensModularDomain(Gamma0(11), None, None, None, None)
        Traceback (most recent call last):
        ...
        TypeError: unable to coerce <class 'sage.modular.arithgroup.congroup_gamma0.Gamma0_class_with_category'> to an integer

    """
    def __init__(self, N, reps, indices, rels, equiv_ind):
        r"""
        INPUT:

            See :class:`PollackStevensModularDomain`.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import PollackStevensModularDomain, ManinRelations
            sage: isinstance(ManinRelations(11), PollackStevensModularDomain) # indirect doctest
            True
        """
        self._N = ZZ(N)
        self._reps = reps

        self._indices = sorted(indices)
        self._gens = [M2Z(reps[i]) for i in self._indices]
        self._ngens = len(indices)

        if len(rels) != len(reps):
            raise ValueError("length of reps and length of rels must be equal")
        self._rels = rels
        self._rel_dict = {}
        for j, L in enumerate(rels):
            self._rel_dict[reps[j]] = L

        self._equiv_ind = equiv_ind
        self._equiv_rep = {}
        for ky in equiv_ind:
            self._equiv_rep[ky] = reps[equiv_ind[ky]]

    def _repr_(self):
        r"""
        A string representation of this domain.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import PollackStevensModularDomain, M2Z
            sage: PollackStevensModularDomain(2 , [M2Z([1,0,0,1]), M2Z([1,1,-1,0]), M2Z([0,-1,1,1])], [0,2], [[(1, M2Z([1,0,0,1]), 0)], [(-1,M2Z([-1,-1,0,-1]),0)], [(1, M2Z([1,0,0,1]), 2)]], {(0,1): 0, (1,0): 1, (1,1): 2})._repr_()
            'Modular Symbol domain of level 2'
        """
        return "Modular Symbol domain of level %s" % self._N

    def __len__(self):
        r"""
        Return the number of coset representatives.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(11)
            sage: len(A)
            12
        """
        return len(self._reps)

    def __getitem__(self, i):
        r"""
        Return the ``i``-th coset representative.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(11)
            sage: A[4]
            [-1 -2]
            [ 2  3]
        """
        return self._reps[i]

    def __iter__(self):
        r"""
        Return an iterator over all coset representatives.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(11)
            sage: for rep in A:
            ....:     if rep[1,0] == 1:
            ....:         print(rep)
            [ 0 -1]
            [ 1  3]
            [ 0 -1]
            [ 1  2]
            [ 0 -1]
            [ 1  1]
        """
        return iter(self._reps)

    def gens(self):
        r"""
        Return the list of coset representatives chosen as generators.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(11)
            sage: A.gens()
            [
            [1 0]  [ 0 -1]  [-1 -1]
            [0 1], [ 1  3], [ 3  2]
            ]
        """
        return self._gens

    def gen(self, n=0):
        r"""
        Return the ``n``-th generator.

        INPUT:

        - ``n`` -- integer (default: 0), which generator is desired

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(137)
            sage: A.gen(17)
            [-4 -1]
            [ 9  2]
        """
        return self._gens[n]

    def ngens(self):
        r"""
        Return the number of generators.

        OUTPUT:

        The number of coset representatives from which a modular symbol's value
        on any coset can be derived.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(1137)
            sage: A.ngens()
            255
        """
        return len(self._gens)

    def level(self):
        r"""
        Return the level `N` of `\Gamma_0(N)` that we work with.

        OUTPUT:

        The integer `N` of the group `\Gamma_0(N)` for which the Manin
        Relations are being computed.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(11)
            sage: A.level()
            11
        """
        return self._N

    def indices(self, n=None):
        r"""
        Return the `n`-th index of the coset representatives which were
        chosen as our generators.

        In particular, the divisors associated to these coset representatives
        generate all divisors over `\ZZ[\Gamma_0(N)]`, and thus a modular
        symbol is uniquely determined by its values on these divisors.

        INPUT:

        - ``n`` -- integer (default: None)

        OUTPUT:

        The ``n``-th index of the generating set in ``self.reps()`` or all
        indices if ``n`` is ``None``.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(11)
            sage: A.indices()
            [0, 2, 3]

            sage: A.indices(2)
            3

            sage: A = ManinRelations(13)
            sage: A.indices()
            [0, 2, 3, 4, 5]

            sage: A = ManinRelations(101)
            sage: A.indices()
            [0, 2, 3, 4, 5, 6, 8, 9, 11, 13, 14, 16, 17, 19, 20, 23, 24, 26, 28]
        """
        if n is None:
            return self._indices
        else:
            return self._indices[n]

    def reps(self, n=None):
        r"""
        Return the ``n``-th coset representative associated with our
        fundamental domain.

        INPUT:

        - ``n`` -- integer (default: None)

        OUTPUT:

        The ``n``-th coset representative or all coset representatives if ``n``
        is ``None``.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(11)
            sage: A.reps(0)
            [1 0]
            [0 1]
            sage: A.reps(1)
            [ 1  1]
            [-1  0]
            sage: A.reps(2)
            [ 0 -1]
            [ 1  3]
            sage: A.reps()
            [
            [1 0]  [ 1  1]  [ 0 -1]  [-1 -1]  [-1 -2]  [-2 -1]  [ 0 -1]  [ 1  0]
            [0 1], [-1  0], [ 1  3], [ 3  2], [ 2  3], [ 3  1], [ 1  2], [-2  1],
            <BLANKLINE>
            [ 0 -1]  [ 1  0]  [-1 -1]  [ 1 -1]
            [ 1  1], [-1  1], [ 2  1], [-1  2]
            ]
        """
        if n is None:
            return self._reps
        else:
            return self._reps[n]

    def relations(self, A=None):
        r"""
        Express the divisor attached to the coset representative of ``A`` in
        terms of our chosen generators.

        INPUT:

        - ``A`` -- ``None``, an integer, or a coset representative (default:
          ``None``)

        OUTPUT:

        A `\ZZ[\Gamma_0(N)]`-relation expressing the divisor attached to ``A``
        in terms of the generating set. The relation is given as a list of
        triples ``(d, B, i)`` such that the divisor attached to `A`` is the sum
        of ``d`` times the divisor attached to ``B^{-1} * self.reps(i)``.

        If ``A`` is an integer, then return this data for the ``A``-th
        coset representative.

        If ``A`` is ``None``, then return this data in a list for all coset
        representatives.

        .. NOTE::

            These relations allow us to recover the value of a modular symbol
            on any coset representative in terms of its values on our
            generating set.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: MR = ManinRelations(11)
            sage: MR.indices()
            [0, 2, 3]
            sage: MR.relations(0)
            [(1, [1 0]
            [0 1], 0)]
            sage: MR.relations(2)
            [(1, [1 0]
            [0 1], 2)]
            sage: MR.relations(3)
            [(1, [1 0]
            [0 1], 3)]

        The fourth coset representative can be expressed through the
        second coset representative::

            sage: MR.reps(4)
            [-1 -2]
            [ 2  3]
            sage: d, B, i = MR.relations(4)[0]
            sage: P = B.inverse()*MR.reps(i); P
            [ 2 -1]
            [-3  2]
            sage: d # the above corresponds to minus the divisor of A.reps(4) since d is -1
            -1

        The sixth coset representative can be expressed as the sum of
        the second and the third::

            sage: MR.reps(6)
            [ 0 -1]
            [ 1  2]
            sage: MR.relations(6)
            [(1, [1 0]
            [0 1], 2), (1, [1 0]
            [0 1], 3)]
            sage: MR.reps(2), MR.reps(3) # MR.reps(6) is the sum of these divisors
            (
            [ 0 -1]  [-1 -1]
            [ 1  3], [ 3  2]
            )

        TESTS:

        Test that the other ways of calling this method work::

            sage: MR.relations(MR.reps(6))
            [(1, [1 0]
            [0 1], 2), (1, [1 0]
            [0 1], 3)]
            sage: MR.relations(None)
            [[(1, [1 0]
            [0 1], 0)], [(-1, [-1 -1]
            [ 0 -1], 0)], [(1, [1 0]
            [0 1], 2)], [(1, [1 0]
            [0 1], 3)], [(-1, [-3 -2]
            [11  7], 2)], [(-1, [-4 -3]
            [11  8], 3)], [(1, [1 0]
            [0 1], 2), (1, [1 0]
            [0 1], 3)], [(-1, [1 0]
            [0 1], 2), (-1, [1 0]
            [0 1], 3)], [(1, [1 0]
            [0 1], 2), (1, [1 0]
            [0 1], 3), (-1, [-3 -2]
            [11  7], 2), (-1, [-4 -3]
            [11  8], 3)], [(-1, [1 0]
            [0 1], 2), (-1, [1 0]
            [0 1], 3), (1, [-3 -2]
            [11  7], 2), (1, [-4 -3]
            [11  8], 3)], [(-1, [-3 -2]
            [11  7], 2), (-1, [-4 -3]
            [11  8], 3)], [(1, [-3 -2]
            [11  7], 2), (1, [-4 -3]
            [11  8], 3)]]
        """
        if A is None:
            return self._rels
        elif isinstance(A, (int, Integer, slice)):
            return self._rels[A]
        else:
            return self._rel_dict[A]

    def equivalent_index(self, A):
        r"""
        Return the index of the coset representative equivalent to ``A``.

        Here by equivalent we mean the unique coset representative whose bottom
        row is equivalent to the bottom row of ``A`` in `P^1(\ZZ/N\ZZ)`.

        INPUT:

        - ``A`` -- an element of `SL_2(\ZZ)`

        OUTPUT:

        The unique integer ``j`` satisfying that the bottom row of
        ``self.reps(j)`` is equivalent to the bottom row of ``A``.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: MR = ManinRelations(11)
            sage: A = matrix(ZZ,2,2,[1,5,3,16])
            sage: j = MR.equivalent_index(A); j
            11
            sage: MR.reps(11)
            [ 1 -1]
            [-1  2]
            sage: MR.equivalent_rep(A)
            [ 1 -1]
            [-1  2]
            sage: MR.P1().normalize(3,16)
            (1, 9)
        """
        return self._equiv_ind[self._P.normalize(A[t10], A[t11])]

    def equivalent_rep(self, A):
        r"""
        Return a coset representative that is equivalent to ``A`` modulo
        `\Gamma_0(N)`.

        INPUT:

        - ``A`` -- a matrix in `SL_2(\ZZ)`

        OUTPUT:

        The unique generator congruent to ``A`` modulo `\Gamma_0(N)`.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = matrix([[5,3],[38,23]])
            sage: ManinRelations(60).equivalent_rep(A)
            [-7 -3]
            [26 11]
        """
        return self._reps[self.equivalent_index(A)]

    def P1(self):
        r"""
        Return the Sage representation of `P^1(\ZZ/N\ZZ)`.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(11)
            sage: A.P1()
            The projective line over the integers modulo 11

        """
        return self._P


class ManinRelations(PollackStevensModularDomain):
    r"""
    This class gives a description of `Div^0(P^1(\QQ))` as a
    `\ZZ[\Gamma_0(N)]`-module.

    INPUT:

    - ``N`` -- a positive integer, the level of `\Gamma_0(N)` to work with

    EXAMPLES::

        sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
        sage: ManinRelations(1)
        Manin Relations of level 1
        sage: ManinRelations(11)
        Manin Relations of level 11

    Large values of ``N`` are not supported::

        sage: ManinRelations(2^20)
        Traceback (most recent call last):
        ...
        OverflowError: Modulus is too large (must be <= 46340)

    TESTS:

    ``N`` has to be a positive integer::

        sage: ManinRelations(0)
        Traceback (most recent call last):
        ...
        ValueError: N must be a positive integer
        sage: ManinRelations(-5)
        Traceback (most recent call last):
        ...
        ValueError: N must be a positive integer

    """
    def __init__(self, N):
        r"""
        Create an instance of this class.

        INPUT:

        - ``N`` -- a positive integer, the level of `\Gamma_0(N)` to work with

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: type(ManinRelations(30))
            <class 'sage.modular.pollack_stevens.fund_domain.ManinRelations'>
        """
        N = ZZ(N)
        if N <= 0:
            raise ValueError("N must be a positive integer")
        self._N = N
        SN = Sigma0(N)

        ## Creates and stores the Sage representation of P^1(Z/NZ)
        P = P1List(N)
        self._P = P
        IdN = SN([1, 0, 0, 1])

        ## Creates a fundamental domain for Gamma_0(N) whose boundary
        ## is a union of unimodular paths (except in the case of
        ## 3-torsion).  We will call the intersection of this domain
        ## with the real axis the collection of cusps (even if some
        ## are Gamma_0(N) equivalent to one another).
        cusps = self.form_list_of_cusps()

        ## Takes the boundary of this fundamental domain and finds
        ## SL_2(Z) matrices whose associated unimodular path gives
        ## this boundary.  These matrices form the beginning of our
        ## collection of coset reps for Gamma_0(N) / SL_2(Z).
        coset_reps = self.fd_boundary(cusps)

        ## Takes the bottom row of each of our current coset reps,
        ## thinking of them as distinct elements of P^1(Z/NZ)
        p1s = [(coset_reps[j])[1] for j in range(len(coset_reps))]

        ## Initializes relevant Manin data
        gens_index = []
        twotor_index = []
        twotorrels = []
        threetor_index = []
        threetorrels = []
        rels = [0] * len(coset_reps)
        gammas = {}

        ## the list rels (above) will give Z[Gamma_0(N)] relations between
        ## the associated divisor of each coset representatives in terms
        ## of our chosen set of generators.
        ## entries of rels will be lists of elements of the form (c,A,r)
        ## with c a constant, A a Gamma_0(N) matrix, and r the index of a
        ## generator.  The meaning is that the divisor associated to the
        ## j-th coset rep will equal the sum of:
        ##
        ##   c * A^(-1) * (divisor associated to r-th coset rep)
        ##
        ## as one varies over all (c,A,r) in rels[j].
        ## (Here r must be in self.generator_indices().)
        ##
        ## This will be used for modular symbols as then the value of a
        ## modular symbol phi on the (associated divisor) of the j-th
        ## element of coset_reps will be the sum of c * phi (r-th generator) | A
        ## as one varies over the tuples in rels[j]

        boundary_checked = [False] * len(coset_reps)

        ## The list boundary_checked keeps track of which boundary pieces of the
        ## fundamental domain have been already used as we are picking
        ## our generators

        ## The following loop will choose our generators by picking one edge
        ## out of each pair of edges that are glued to each other and picking
        ## each edge glued to itself (arising from two-torsion)
        ## ------------------------------------------------------------------
        for r in range(len(coset_reps)):
            if not boundary_checked[r]:

                ## We now check if this boundary edge is glued to itself by
                ## Gamma_0(N)

                if P.normalize(p1s[r][0], p1s[r][1]) == P.normalize(-p1s[r][1], p1s[r][0]):
                    ## This edge is glued to itself and so coset_reps[r]
                    ## needs to be added to our generator list.

                    ## this relation expresses the fact that
                    ## coset_reps[r] is one of our basic generators
                    rels[r] = [(1, IdN, r)]

                    ## the index r is adding to our list
                    ## of indexes of generators
                    gens_index.append(r)

                    ## the index r is adding to our list of indexes of
                    ## generators which satisfy a 2-torsion relation
                    twotor_index.append(r)

                    # we use the adjugate instead of the inverse for speed
                    gam = SN(coset_reps[r] * sig * coset_reps[r].adjugate())
                    ## gam is 2-torsion matrix and in Gamma_0(N).
                    ## if D is the divisor associated to coset_reps[r]
                    ## then gam * D = - D and so (1+gam)D=0.

                    ## This gives a restriction to the possible values of
                    ## modular symbols on D

                    ## The 2-torsion matrix gam is recorded in our list of
                    ## 2-torsion relations.
                    twotorrels.append(gam)

                    ## We have now finished with this edge.
                    boundary_checked[r] = True

                else:
                    c = coset_reps[r][t10]
                    d = coset_reps[r][t11]

                    ## In the following case the ideal triangle below
                    ## the unimodular path described by coset_reps[r]
                    ## contains a point fixed by a 3-torsion element.
                    if (c ** 2 + d ** 2 + c * d) % N == 0:

                        ## the index r is adding to our list of indexes
                        ## of generators
                        gens_index.append(r)

                        ## this relation expresses the fact that coset_reps[r]
                        ## is one of our basic generators
                        rels[r] = [(1, IdN, r)]

                        ## the index r is adding to our list of indexes of
                        ## generators which satisfy a 3-torsion relation
                        threetor_index.append(r)

                        # Use the adjugate instead of the inverse for speed.
                        gam = SN(coset_reps[r] * tau * coset_reps[r].adjugate())
                        ## gam is 3-torsion matrix and in Gamma_0(N).
                        ## if D is the divisor associated to coset_reps[r]
                        ## then (1+gam+gam^2)D=0.
                        ## This gives a restriction to the possible values of
                        ## modular symbols on D

                        ## The 3-torsion matrix gam is recorded in our list of
                        ## 3-torsion relations.
                        threetorrels.append(gam)

                        ## The reverse of the unimodular path associated to
                        ## coset_reps[r] is not Gamma_0(N) equivalent to it, so
                        ## we need to include it in our list of coset
                        ## representatives and record the relevant relations.

                        a = coset_reps[r][t00]
                        b = coset_reps[r][t01]

                        A = M2Z([-b, a, -d, c])
                        coset_reps.append(A)
                        ## A (representing the reversed edge) is included in
                        ## our list of coset reps

                        rels.append([(-1, IdN, r)])
                        ## This relation means that phi on the reversed edge
                        ## equals -phi on original edge

                        boundary_checked[r] = True
                        ## We have now finished with this edge.

                    else:
                        ## This is the generic case where neither 2 or
                        ## 3-torsion intervenes.
                        ## The below loop searches through the remaining edges
                        ## and finds which one is equivalent to the reverse of
                        ## coset_reps[r]
                        ## ---------------------------------------------------
                        for s in range(r + 1, len(coset_reps)):
                            if boundary_checked[s]:
                                continue
                            if P.normalize(p1s[s][0], p1s[s][1]) == P.normalize(-p1s[r][1], p1s[r][0]):
                                ## the reverse of coset_reps[r] is
                                ## Gamma_0(N)-equivalent to coset_reps[s]
                                ## coset_reps[r] will now be made a generator
                                ## and we need to express phi(coset_reps[s])
                                ## in terms of phi(coset_reps[r])

                                gens_index.append(r)
                                ## the index r is adding to our list of
                                ## indexes of generators

                                rels[r] = [(1, IdN, r)]
                                ## this relation expresses the fact that
                                ## coset_reps[r] is one of our basic generators

                                A = coset_reps[s] * sig
                                ## A corresponds to reversing the orientation
                                ## of the edge corr. to coset_reps[r]
                                # Use adjugate instead of inverse for speed
                                gam = SN(coset_reps[r] * A.adjugate())
                                ## gam is in Gamma_0(N) (by assumption of
                                ## ending up here in this if statement)

                                rels[s] = [(-1, gam, r)]
                                ## this relation means that phi evaluated on
                                ## coset_reps[s] equals -phi(coset_reps[r])|gam
                                ## To see this, let D_r be the divisor
                                ## associated to coset_reps[r] and D_s to
                                ## coset_reps[s]. Then gam D_s = -D_r and so
                                ## phi(gam D_s) = - phi(D_r) and thus
                                ## phi(D_s) = -phi(D_r)|gam
                                ## since gam is in Gamma_0(N)

                                gammas[coset_reps[r]] = gam
                                ## this is a dictionary whose keys are the
                                ## non-torsion generators and whose values
                                ## are the corresponding gamma_i. It is
                                ## eventually stored as self.gammas.

                                boundary_checked[r] = True
                                boundary_checked[s] = True
                                break

        ## We now need to complete our list of coset representatives by
        ## finding all unimodular paths in the interior of the fundamental
        ## domain, as well as express these paths in terms of our chosen set
        ## of generators.
        ## -------------------------------------------------------------------

        for r in range(len(cusps) - 2):
            ## r is the index of the cusp on the left of the path.  We only run
            ## thru to the number of cusps - 2 since you cannot start an
            ## interior path on either of the last two cusps

            for s in range(r + 2, len(cusps)):
            ## s is in the index of the cusp on the right of the path
                cusp1 = cusps[r]
                cusp2 = cusps[s]
                if self.is_unimodular_path(cusp1, cusp2):
                    A, B = self.unimod_to_matrices(cusp1, cusp2)
                    ## A and B are the matrices whose associated paths
                    ## connect cusp1 to cusp2 and cusp2 to cusp1 (respectively)
                    coset_reps.extend([A, B])
                    ## A and B are added to our coset reps
                    vA = []
                    vB = []

                    ## This loop now encodes the relation between the
                    ## unimodular path A and our generators.  This is done
                    ## simply by accounting for all of the edges that lie
                    ## below the path attached to A (as they form a triangle)
                    ## Similarly, this is also done for B.

                    ## Running between the cusps between cusp1 and cusp2
                    for rel in rels[r + 2: s + 2]:
                        ## Add edge relation
                        vA.append(rel[0])
                        ## Add negative of edge relation
                        vB.append((-rel[0][0], rel[0][1], rel[0][2]))
                    ## Add relations for A and B to relations list
                    rels.extend([vA, vB])

        ## Make the translation table between the Sage and Geometric
        ## descriptions of P^1
        equiv_ind = {}
        for i, rep in enumerate(coset_reps):
            ky = P.normalize(rep[t10], rep[t11])
            equiv_ind[ky] = i

        self.gammas = gammas
        PollackStevensModularDomain.__init__(self, N, coset_reps, gens_index,
                                        rels, equiv_ind)

        ## A list of indices of the (geometric) coset representatives whose
        ## paths are identified by some 2-torsion element (which switches the
        ## path orientation)
        self._indices_with_two_torsion = twotor_index
        self._reps_with_two_torsion = [coset_reps[i] for i in twotor_index]

        ## A dictionary of (2-torsion in PSL_2(Z)) matrices in
        ## Gamma_0(N) that give the orientation identification in the
        ## paths listed in twotor_index above!
        self._two_torsion = {}
        for j, tor_elt in zip(twotor_index, twotorrels):
            self._two_torsion[coset_reps[j]] = tor_elt

        ## A list of indices of the (geometric) coset representatives that
        ## form one side of an ideal triangle with an interior fixed point of
        ## a 3-torsion element of Gamma_0(N)
        self._indices_with_three_torsion = threetor_index
        self._reps_with_three_torsion = [coset_reps[i] for i in threetor_index]

        ## A dictionary of (3-torsion in PSL_2(Z)) matrices in
        ## Gamma_0(N) that give the interior fixed point described in
        ## threetor_index above!
        self._three_torsion = {}
        for j, tor_elt in zip(threetor_index, threetorrels):
            self._three_torsion[coset_reps[j]] = tor_elt

    def _repr_(self):
        r"""
        A printable representation of this domain.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: ManinRelations(11)._repr_()
            'Manin Relations of level 11'
        """
        return "Manin Relations of level %s" % self._N

    def indices_with_two_torsion(self):
        r"""
        Return the indices of coset representatives whose associated unimodular path
        contains a point fixed by a `\Gamma_0(N)` element of order 2 (where the
        order is computed in `PSL_2(\ZZ)`).

        OUTPUT:

        A list of integers.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: MR = ManinRelations(11)
            sage: MR.indices_with_two_torsion()
            []
            sage: MR = ManinRelations(13)
            sage: MR.indices_with_two_torsion()
            [3, 4]
            sage: MR.reps(3), MR.reps(4)
            (
            [-1 -1]  [-1 -2]
            [ 3  2], [ 2  3]
            )

        The corresponding matrix of order 2::

            sage: A = MR.two_torsion_matrix(MR.reps(3)); A
            [  5   2]
            [-13  -5]
            sage: A^2
            [-1  0]
            [ 0 -1]

        You can see that multiplication by ``A`` just interchanges the rational
        cusps determined by the columns of the matrix ``MR.reps(3)``::

            sage: MR.reps(3), A*MR.reps(3)
            (
            [-1 -1]  [ 1 -1]
            [ 3  2], [-2  3]
            )
        """
        return self._indices_with_two_torsion

    def reps_with_two_torsion(self):
        r"""
        The coset representatives whose associated unimodular path contains a
        point fixed by a `\Gamma_0(N)` element of order 2 (where the order is
        computed in `PSL_2(\ZZ)`).

        OUTPUT:

        A list of matrices.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: MR = ManinRelations(11)
            sage: MR.reps_with_two_torsion()
            []
            sage: MR = ManinRelations(13)
            sage: MR.reps_with_two_torsion()
            [
            [-1 -1]  [-1 -2]
            [ 3  2], [ 2  3]
            ]
            sage: B = MR.reps_with_two_torsion()[0]

        The corresponding matrix of order 2::

            sage: A = MR.two_torsion_matrix(B); A
            [  5   2]
            [-13  -5]
            sage: A^2
            [-1  0]
            [ 0 -1]

        You can see that multiplication by ``A`` just interchanges the rational
        cusps determined by the columns of the matrix ``MR.reps(3)``::

            sage: B, A*B
            (
            [-1 -1]  [ 1 -1]
            [ 3  2], [-2  3]
            )
        """
        return self._reps_with_two_torsion

    def two_torsion_matrix(self, A):
        r"""
        Return the matrix of order two in `\Gamma_0(N)` which
        corresponds to an ``A`` in ``self.reps_with_two_torsion()``.

        INPUT:

        - ``A`` -- a matrix in ``self.reps_with_two_torsion()``

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: MR = ManinRelations(25)
            sage: B = MR.reps_with_two_torsion()[0]

        The corresponding matrix of order 2::

            sage: A = MR.two_torsion_matrix(B); A
            [  7   2]
            [-25  -7]
            sage: A^2
            [-1  0]
            [ 0 -1]
        """
        return self._two_torsion[A]

    def indices_with_three_torsion(self):
        r"""
        A list of indices of coset representatives whose associated unimodular
        path contains a point fixed by a `\Gamma_0(N)` element of order 3 in
        the ideal triangle directly below that path (the order is computed in
        `PSL_2(\ZZ)`).

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: MR = ManinRelations(11)
            sage: MR.indices_with_three_torsion()
            []
            sage: MR = ManinRelations(13)
            sage: MR.indices_with_three_torsion()
            [2, 5]
            sage: B = MR.reps(2); B
            [ 0 -1]
            [ 1  3]

        The corresponding matrix of order three::

            sage: A = MR.three_torsion_matrix(B); A
            [-4 -1]
            [13  3]
            sage: A^3
            [1 0]
            [0 1]

        The columns of ``B`` and the columns of ``A*B`` and ``A^2*B`` give the
        same rational cusps::

            sage: B
            [ 0 -1]
            [ 1  3]
            sage: A*B, A^2*B
            (
            [-1  1]  [ 1  0]
            [ 3 -4], [-4  1]
            )
        """
        return self._indices_with_three_torsion

    def reps_with_three_torsion(self):
        r"""
        A list of coset representatives whose associated unimodular
        path contains a point fixed by a `\Gamma_0(N)` element of
        order 3 in the ideal triangle directly below that path (the
        order is computed in `PSL_2(\ZZ)`).

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: MR = ManinRelations(13)
            sage: B = MR.reps_with_three_torsion()[0]; B
            [ 0 -1]
            [ 1  3]

        The corresponding matrix of order three::

            sage: A = MR.three_torsion_matrix(B); A
            [-4 -1]
            [13  3]
            sage: A^3
            [1 0]
            [0 1]

        The columns of ``B`` and the columns of ``A*B`` and ``A^2*B``
        give the same rational cusps::

            sage: B
            [ 0 -1]
            [ 1  3]
            sage: A*B, A^2*B
            (
            [-1  1]  [ 1  0]
            [ 3 -4], [-4  1]
            )
        """
        return self._reps_with_three_torsion

    def three_torsion_matrix(self, A):
        r"""
        Return the matrix of order two in `\Gamma_0(N)` which
        corresponds to an ``A`` in ``self.reps_with_two_torsion()``.

        INPUT:

        - ``A`` -- a matrix in ``self.reps_with_two_torsion()``

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: MR = ManinRelations(37)
            sage: B = MR.reps_with_three_torsion()[0]

        The corresponding matrix of order 3::

            sage: A = MR.three_torsion_matrix(B); A
            [-11  -3]
            [ 37  10]
            sage: A^3
            [1 0]
            [0 1]
        """
        return self._three_torsion[A]

    def form_list_of_cusps(self):
        r"""
        Return the intersection of a fundamental domain for `\Gamma_0(N)` with
        the real axis.

        The construction of this fundamental domain follows the arguments of
        [PS2011]_ Section 2.  The boundary of this fundamental domain consists
        entirely of unimodular paths when `\Gamma_0(N)` has no elements of
        order 3.  (See [PS2011]_ Section 2.5 for the case when there are
        elements of order 3.)

        OUTPUT:

        A sorted list of rational numbers marking the intersection of a
        fundamental domain for `\Gamma_0(N)` with the real axis.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(11)
            sage: A.form_list_of_cusps()
            [-1, -2/3, -1/2, -1/3, 0]
            sage: A = ManinRelations(13)
            sage: A.form_list_of_cusps()
            [-1, -2/3, -1/2, -1/3, 0]
            sage: A = ManinRelations(101)
            sage: A.form_list_of_cusps()
            [-1, -6/7, -5/6, -4/5, -7/9, -3/4, -11/15, -8/11, -5/7, -7/10,
            -9/13, -2/3, -5/8, -13/21, -8/13, -3/5, -7/12, -11/19, -4/7, -1/2,
            -4/9, -3/7, -5/12, -7/17, -2/5, -3/8, -4/11, -1/3, -2/7, -3/11,
            -1/4, -2/9, -1/5, -1/6, 0]
        """
        ## Get the level
        N = self.level()

        ## Some convenient shortcuts
        P = self.P1()
        sP = len(P.list())   # Size of P^1(Z/NZ)

        ## Initialize some lists

        C = [QQ(-1), "?", QQ(0)]

        ## Initialize the list of cusps at the bottom of the fund. domain.
        ## The ? denotes that it has not yet been checked if more cusps need
        ## to be added between the surrounding cusps.

        full_domain = False     # Says that we are not done yet!

        v = [False for r in range(sP)]
        ## This initializes a list indexed by P^1(Z/NZ) which keeps track of
        ## which right coset representatives we've found for Gamma_0(N)/SL_2(Z)
        ## thru the construction of a fundamental domain

        ## Includeds the coset repns formed by the original ideal triangle
        ## (with corners at -1, 0, infty)

        v[P.index(0, 1)] = True
        v[P.index(1, -1)] = True
        v[P.index(-1, 0)] = True

        ## Main Loop -- Ideal Triangle Flipping
        ## ====================================
        while (not full_domain):
            full_domain = True

            ## This loop runs through the current set of cusps
            ## and checks to see if more cusps should be added
            ## -----------------------------------------------
            for s in range(1, len(C), 2):  # range over odd indices in the
                                           # final list C
                if C[s] == "?":

                    ## Single out our two cusps (path from cusp2 to cusp1)
                    cusp1 = C[s - 1]
                    cusp2 = C[s + 1]

                    ## Makes the unimodular transform for the path from cusp2
                    ## to cusp1

                    b1 = cusp1.denominator()
                    b2 = cusp2.denominator()

                    ## This is the point where it is determined whether
                    ## or not the adjacent triangle should be added
                    ## ------------------------------------------------
                    pos = P.index(b2, b1)   # The Sage index of the bottom
                                                 ## row of our unimodular
                                           ## transformation gam

                    ## Check if we need to flip (since this P1 element has not
                    ## yet been accounted for!)
                    if not v[pos]:
                        v[pos] = True      # Say this P1 element now occurs
                        v[P.index(b1, -(b1 + b2))] = True
                        ## Say that the other two ideal triangle edges
                        ## also occur!

                        v[P.index(-(b1 + b2), b2)] = True

                        ## Check to see if this triangle contains a fixed
                        ## point by an element of Gamma_0(N).  If such an
                        ## element is present, the fundamental domain can be
                        ## extended no further.

                        if (b1 ** 2 + b2 ** 2 + b1 * b2) % N != 0:

                        ## this congruence is exactly equivalent to
                        ## gam * [0 -1; 1 -1] * gam^(-1) is in Gamma_0(N)
                        ## where gam is the matrix corresponding to the
                        ## unimodular path connecting cusp1 to cusp2

                            C[s] = "i"  # The '?' is changed to an 'i'
                             ## indicating that a new cusp needs to
                                        ##  be inserted here
                            full_domain = False
                        else:
                            C[s] = "x"  # The '?' is changed to an 'x' and no
                                        # more checking below is needed! =)
                    else:
                        C[s] = "x"  # The '?' is changed to an 'x' and no more
                                           ## checking below is needed! =)

            ## Now insert the missing cusps (where there is an 'i' in
            ## the final list)
            ## This will keep the fundamental domain as flat as possible!
            ## ---------------------------------------------------------------
            s = 1
            while s < len(C):   # range over odd indices in the final list C
                if C[s] == "i":
                    C[s] = "?"

                    ## Single out our two cusps (path from cusp2 to cusp1)
                    cusp1 = C[s - 1]
                    cusp2 = C[s + 1]

                    ## Makes the unimodular transform for the path
                    ## from cusp2 to cusp1
                    a1 = cusp1.numerator()
                    b1 = cusp1.denominator()
                    a2 = cusp2.numerator()
                    b2 = cusp2.denominator()

                    ## Inserts the Farey center of these two cusps!
                    a = a1 + a2
                    b = b1 + b2
                    C.insert(s + 1, a / b)
                    C.insert(s + 2, "?")
                    s += 2
                s += 2

        ## Remove the (now superfluous) extra string characters that appear
        ## in the odd list entries
        C = [QQ(C[ss]) for ss in range(0, len(C), 2)]
        return C

    def is_unimodular_path(self, r1, r2):
        r"""
        Determine whether two (non-infinite) cusps are connected by a
        unimodular path.

        INPUT:

        - ``r1, r2`` -- rational numbers

        OUTPUT:

        A boolean expressing whether or not a unimodular path connects ``r1``
        to ``r2``.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(11)
            sage: A.is_unimodular_path(0,1/3)
            True
            sage: A.is_unimodular_path(1/3,0)
            True
            sage: A.is_unimodular_path(0,2/3)
            False
            sage: A.is_unimodular_path(2/3,0)
            False
        """
        a = r1.numerator()
        b = r2.numerator()
        c = r1.denominator()
        d = r2.denominator()
        return (a * d - b * c) ** 2 == 1

    def unimod_to_matrices(self, r1, r2):
        r"""
        Return the two matrices whose associated unimodular paths connect
        ``r1`` and ``r2`` and ``r2`` and ``r1``, respectively.

        INPUT:

        - ``r1, r2`` -- rational numbers (that are assumed to be connected by a
          unimodular path)

        OUTPUT:

        A pair of `2 \times 2` matrices of determinant 1

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(11)
            sage: A.unimod_to_matrices(0,1/3)
            (
            [ 0  1]  [1 0]
            [-1  3], [3 1]
            )
        """
        a = r1.numerator()
        b = r2.numerator()
        c = r1.denominator()
        d = r2.denominator()
        if (a * d - b * c) == 1:
            ans = M2Z([a, b, c, d]), M2Z([-b, a, -d, c])
        else:
            ans = M2Z([-a, b, -c, d]), M2Z([b, a, d, c])
        return ans

    def fd_boundary(self, C):
        r"""
        Find matrices whose associated unimodular paths give the
        boundary of a fundamental domain.

        Here the fundamental domain is for `\Gamma_0(N)`.  (In the
        case when `\Gamma_0(N)` has elements of order three the shape
        cut out by these unimodular matrices is a little smaller than
        a fundamental domain.  See Section 2.5 of [PS2011]_.)

        INPUT:

        - ``C`` -- a list of rational numbers coming from
          ``self.form_list_of_cusps()``

        OUTPUT:

        A list of `2 \times 2` integer matrices of determinant 1 whose associated
        unimodular paths give the boundary of a fundamental domain for
        `\Gamma_0(N)` (or nearly so in the case of 3-torsion).

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.fund_domain import ManinRelations
            sage: A = ManinRelations(11)
            sage: C = A.form_list_of_cusps(); C
            [-1, -2/3, -1/2, -1/3, 0]
            sage: A.fd_boundary(C)
            [
            [1 0]  [ 1  1]  [ 0 -1]  [-1 -1]  [-1 -2]  [-2 -1]
            [0 1], [-1  0], [ 1  3], [ 3  2], [ 2  3], [ 3  1]
            ]
            sage: A = ManinRelations(13)
            sage: C = A.form_list_of_cusps(); C
            [-1, -2/3, -1/2, -1/3, 0]
            sage: A.fd_boundary(C)
            [
            [1 0]  [ 1  1]  [ 0 -1]  [-1 -1]  [-1 -2]  [-2 -1]
            [0 1], [-1  0], [ 1  3], [ 3  2], [ 2  3], [ 3  1]
            ]
            sage: A = ManinRelations(101)
            sage: C = A.form_list_of_cusps(); C
            [-1, -6/7, -5/6, -4/5, -7/9, -3/4, -11/15, -8/11, -5/7, -7/10,
            -9/13, -2/3, -5/8, -13/21, -8/13, -3/5, -7/12, -11/19, -4/7, -1/2,
            -4/9, -3/7, -5/12, -7/17, -2/5, -3/8, -4/11, -1/3, -2/7, -3/11,
            -1/4, -2/9, -1/5, -1/6, 0]
            sage: A.fd_boundary(C)
            [
            [1 0]  [ 1  1]  [ 0 -1]  [-1 -1]  [-1 -2]  [-2 -1]  [-1 -3]  [-3 -2]
            [0 1], [-1  0], [ 1  6], [ 6  5], [ 5  9], [ 9  4], [ 4 11], [11  7],
            <BLANKLINE>
            [-2 -1]  [-1 -4]  [-4 -3]  [-3 -2]  [-2 -7]  [-7 -5]  [-5 -3]  [-3 -4]
            [ 7  3], [ 3 11], [11  8], [ 8  5], [ 5 17], [17 12], [12  7], [ 7  9],
            <BLANKLINE>
            [-4 -1]  [-1 -4]  [ -4 -11]  [-11  -7]  [-7 -3]  [-3 -8]  [ -8 -13]
            [ 9  2], [ 2  7], [  7  19], [ 19  12], [12  5], [ 5 13], [ 13  21],
            <BLANKLINE>
            [-13  -5]  [-5 -2]  [-2 -9]  [-9 -7]  [-7 -5]  [-5 -8]  [ -8 -11]
            [ 21   8], [ 8  3], [ 3 13], [13 10], [10  7], [ 7 11], [ 11  15],
            <BLANKLINE>
            [-11  -3]  [-3 -7]  [-7 -4]  [-4 -5]  [-5 -6]  [-6 -1]
            [ 15   4], [ 4  9], [ 9  5], [ 5  6], [ 6  7], [ 7  1]
            ]
        """
        C.reverse()  # Reverse here to get clockwise orientation of boundary

        ## These matrices correspond to the paths from infty to 0 and
        ## -1 to infty
        mats = [Id, minone_inf_path]

        ## Now find SL_2(Z) matrices whose associated unimodular paths
        ## connect the cusps listed in C.
        for j in range(len(C) - 1):
            a = C[j].numerator()
            b = C[j + 1].numerator()
            c = C[j].denominator()
            d = C[j + 1].denominator()
            new_mat = M2Z([a, b, c, d])
            mats.append(new_mat)

        return mats

    @cached_method
    def prep_hecke_on_gen(self, l, gen, modulus=None):
        r"""
        This function does some precomputations needed to compute `T_l`.

        In particular, if `\phi` is a modular symbol and `D_m` is the divisor
        associated to the generator ``gen``, to compute `(\phi|T_{l})(D_m)` one
        needs to compute `\phi(\gamma_a D_m)|\gamma_a` where `\gamma_a` runs
        through the `l+1` matrices defining `T_l`.  One
        then takes `\gamma_a D_m` and writes it as a sum of unimodular
        divisors.  For each such unimodular divisor, say `[M]` where `M` is a
        `SL_2` matrix, we then write `M=\gamma h` where `\gamma` is in
        `\Gamma_0(N)` and `h` is one of our chosen coset representatives.  Then
        `\phi([M]) = \phi([h]) | \gamma^{-1}`.  Thus, one has

        .. MATH::

            (\phi | \gamma_a)(D_m) = \sum_h \sum_j \phi([h]) | \gamma_{hj}^{-1} \cdot \gamma_a

        as `h` runs over all coset representatives and `j` simply runs over
        however many times `M_h` appears in the above computation.

        Finally, the output of this function is a dictionary ``D``
        whose keys are the coset representatives in ``self.reps()``
        where each value is a list of matrices, and the entries of
        ``D`` satisfy:

        .. MATH::

            D[h][j] = \gamma_{hj} * \gamma_a

        INPUT:

        - ``l`` -- a prime
        - ``gen`` -- a generator

        OUTPUT:

        A list of lists (see above).

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: phi = E.pollack_stevens_modular_symbol()
            sage: phi.values()
            [-1/5, 1, 0]
            sage: M = phi.parent().source()
            sage: w = M.prep_hecke_on_gen(2, M.gens()[0])
            sage: one = Matrix(ZZ,2,2,1)
            sage: one.set_immutable()
            sage: w[one]
            [[1 0]
            [0 2], [1 1]
            [0 2], [2 0]
            [0 1]]
        """
        N = self.level()
        SN = Sigma0(N)

        ans = {}
        # this will be the dictionary D above enumerated by coset reps

        #  This loop will run thru the l+1 (or l) matrices
        #  defining T_l of the form [1, a, 0, l] and carry out the
        #  computation described above.
        #  -------------------------------------
        for a in range(l + 1):
            if ((a < l) or (N % l != 0)) and (modulus is None or a % l == modulus % l):
                # if the level is not prime to l the matrix [l, 0, 0, 1] is avoided.
                gamma = basic_hecke_matrix(a, l)
                t = gamma * gen
                #  In the notation above this is gam_a * D_m
                from .manin_map import unimod_matrices_to_infty, unimod_matrices_from_infty
                v = unimod_matrices_from_infty(t[0, 0], t[1, 0]) + unimod_matrices_to_infty(t[0, 1], t[1, 1])
                #  This expresses t as a sum of unimodular divisors

                # This loop runs over each such unimodular divisor
                # ------------------------------------------------
                for A in v:
                    #  B is the coset rep equivalent to A
                    B = self.equivalent_rep(A)
                    #  gaminv = B*A^(-1), but A is in SL2.
                    gaminv = B * A.adjugate()
                    #  The matrix gaminv * gamma is added to our list in the j-th slot
                    #  (as described above)
                    tmp = SN(gaminv * gamma)
                    try:
                        ans[B].append(tmp)
                    except KeyError:
                        ans[B] = [tmp]

        return ans

    @cached_method
    def prep_hecke_on_gen_list(self, l, gen, modulus=None):
        r"""
        Return the precomputation to compute `T_l` in a way that
        speeds up the Hecke calculation.

        Namely, returns a list of the form [h,A].

        INPUT:

        - ``l`` -- a prime
        - ``gen`` -- a generator

        OUTPUT:

        A list of lists (see above).

        EXAMPLES::

            sage: E = EllipticCurve('11a')
            sage: phi = E.pollack_stevens_modular_symbol()
            sage: phi.values()
            [-1/5, 1, 0]
            sage: M = phi.parent().source()
            sage: len(M.prep_hecke_on_gen_list(2, M.gens()[0]))
            4
        """
        ans = []
        for h, vh in self.prep_hecke_on_gen(l, gen, modulus=modulus).items():
            ans.extend([(h, v) for v in vh])
        return ans


def basic_hecke_matrix(a, l):
    r"""
    Return the `2 \times 2` matrix with entries ``[1, a, 0, l]`` if ``a<l``
    and ``[l, 0, 0, 1]`` if ``a>=l``.

    INPUT:

    - `a` -- an integer or Infinity
    - ``l`` -- a prime

    OUTPUT:

    A `2 \times 2` matrix of determinant l

    EXAMPLES::

        sage: from sage.modular.pollack_stevens.fund_domain import basic_hecke_matrix
        sage: basic_hecke_matrix(0, 7)
        [1 0]
        [0 7]
        sage: basic_hecke_matrix(5, 7)
        [1 5]
        [0 7]
        sage: basic_hecke_matrix(7, 7)
        [7 0]
        [0 1]
        sage: basic_hecke_matrix(19, 7)
        [7 0]
        [0 1]
        sage: basic_hecke_matrix(infinity, 7)
        [7 0]
        [0 1]
    """
    if a < l:
        return M2Z([1, a, 0, l])
    else:
        return M2Z([l, 0, 0, 1])
