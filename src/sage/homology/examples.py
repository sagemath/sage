# -*- coding: utf-8 -*-
"""
Examples of simplicial complexes

AUTHORS:

- John H. Palmieri (2009-04)

This file constructs some examples of simplicial complexes.  There are
two main types: manifolds and examples related to graph theory.

For manifolds, there are functions defining the `n`-sphere for any
`n`, the torus, `n`-dimensional real projective space for any `n`, the
complex projective plane, surfaces of arbitrary genus, and some other
manifolds, all as simplicial complexes.

Aside from surfaces, this file also provides some functions for
constructing some other simplicial complexes: the simplicial complex
of not-`i`-connected graphs on `n` vertices, the matching complex on n
vertices, and the chessboard complex for an `n` by `i` chessboard.
These provide examples of large simplicial complexes; for example,
``simplicial_complexes.NotIConnectedGraphs(7,2)`` has over a million
simplices.

All of these examples are accessible by typing
``simplicial_complexes.NAME``, where ``NAME`` is the name of the example.
You can get a list by typing ``simplicial_complexes.`` and hitting the
TAB key::

   simplicial_complexes.BarnetteSphere
   simplicial_complexes.BrucknerGrunbaumSphere
   simplicial_complexes.ChessboardComplex
   simplicial_complexes.ComplexProjectivePlane
   simplicial_complexes.K3Surface
   simplicial_complexes.KleinBottle
   simplicial_complexes.MatchingComplex
   simplicial_complexes.MooreSpace
   simplicial_complexes.NotIConnectedGraphs
   simplicial_complexes.PoincareHomologyThreeSphere
   simplicial_complexes.PseudoQuaternionicProjectivePlane
   simplicial_complexes.RandomComplex
   simplicial_complexes.RealProjectivePlane
   simplicial_complexes.RealProjectiveSpace
   simplicial_complexes.Simplex
   simplicial_complexes.Sphere
   simplicial_complexes.SumComplex
   simplicial_complexes.SurfaceOfGenus
   simplicial_complexes.Torus

See the documentation for ``simplicial_complexes`` and for each
particular type of example for full details.
"""

from sage.homology.simplicial_complex import SimplicialComplex, Simplex
from sage.sets.set import Set
from sage.misc.functional import is_even
from sage.combinat.subset import Subsets
import sage.misc.prandom as random

def matching(A, B):
    r"""
    List of maximal matchings between the sets ``A`` and ``B``.

    A matching is a set of pairs `(a,b) \in A \times B` where each `a` and
    `b` appears in at most one pair.  A maximal matching is one which is
    maximal with respect to inclusion of subsets of `A \times B`.

    INPUT:

    -  ``A``, ``B`` -- list, tuple, or indeed anything which can be
       converted to a set.

    EXAMPLES::

        sage: from sage.homology.examples import matching
        sage: matching([1,2], [3,4])
        [{(1, 3), (2, 4)}, {(1, 4), (2, 3)}]
        sage: matching([0,2], [0])
        [{(0, 0)}, {(2, 0)}]
    """
    answer = []
    if len(A) == 0 or len(B) == 0:
        return [set([])]
    for v in A:
        for w in B:
            for M in matching(set(A).difference([v]), set(B).difference([w])):
                new = M.union([(v,w)])
                if new not in answer:
                    answer.append(new)
    return answer

def facets_for_RP4():
    """
    Return the list of facets for a minimal triangulation of 4-dimensional
    real projective space.

    We use vertices numbered 1 through 16, define two facets, and define
    a certain subgroup `G` of the symmetric group `S_{16}`. Then the set
    of all facets is the `G`-orbit of the two given facets.

    See the description in Example 3.12 in Datta [Da2007]_.

    EXAMPLES::

        sage: from sage.homology.examples import facets_for_RP4
        sage: A = facets_for_RP4()   # long time (1 or 2 seconds)
        sage: SimplicialComplex(A) == simplicial_complexes.RealProjectiveSpace(4) # long time
        True
    """
    # Define the group:
    from sage.groups.perm_gps.permgroup import PermutationGroup
    g1 = '(2,7)(4,10)(5,6)(11,12)'
    g2 = '(1, 2, 3, 4, 5, 10)(6, 8, 9)(11, 12, 13, 14, 15, 16)'
    G = PermutationGroup([g1, g2])
    # Define the two simplices:
    t1 = (1, 2, 4, 5, 11)
    t2 = (1, 2, 4, 11, 13)
    # Apply the group elements to the simplices:
    facets = []
    for g in G:
        d = g.dict()
        for t in [t1, t2]:
            new = tuple([d[j] for j in t])
            if new not in facets:
                facets.append(new)
    return facets

def facets_for_K3():
    """
    Returns the facets for a minimal triangulation of the K3 surface.

    This is a pure simplicial complex of dimension 4 with 16
    vertices and 288 facets. The facets are obtained by constructing a
    few facets and a permutation group `G`, and then computing the
    `G`-orbit of those facets.

    See Casella and Kühnel in [CK2001]_ and Spreer and Kühnel [SK2011]_;
    the construction here uses the labeling from Spreer and Kühnel.

    EXAMPLES::

        sage: from sage.homology.examples import facets_for_K3
        sage: A = facets_for_K3()   # long time (a few seconds)
        sage: SimplicialComplex(A) == simplicial_complexes.K3Surface()  # long time
        True
    """
    from sage.groups.perm_gps.permgroup import PermutationGroup
    G = PermutationGroup([[(1,3,8,4,9,16,15,2,14,12,6,7,13,5,10)],
                         [(1,11,16),(2,10,14),(3,12,13),(4,9,15),(5,7,8)]])
    return ([tuple([g(i) for i in (1,2,3,8,12)]) for g in G]
            +[tuple([g(i) for i in (1,2,5,8,14)]) for g in G])


# for backwards compatibility:
SimplicialSurface = SimplicialComplex

class SimplicialComplexExamples():
    """
    Some examples of simplicial complexes.

    Here are the available examples; you can also type
    ``simplicial_complexes.``  and hit tab to get a list:

    - :meth:`BarnetteSphere`
    - :meth:`BrucknerGrunbaumSphere`
    - :meth:`ChessboardComplex`
    - :meth:`ComplexProjectivePlane`
    - :meth:`K3Surface`
    - :meth:`KleinBottle`
    - :meth:`MatchingComplex`
    - :meth:`MooreSpace`
    - :meth:`NotIConnectedGraphs`
    - :meth:`PoincareHomologyThreeSphere`
    - :meth:`PseudoQuaternionicProjectivePlane`
    - :meth:`RandomComplex`
    - :meth:`RealProjectivePlane`
    - :meth:`RealProjectiveSpace`
    - :meth:`Simplex`
    - :meth:`Sphere`
    - :meth:`SumComplex`
    - :meth:`SurfaceOfGenus`
    - :meth:`Torus`

    EXAMPLES::

        sage: S = simplicial_complexes.Sphere(2) # the 2-sphere
        sage: S.homology()
        {0: 0, 1: 0, 2: Z}
        sage: simplicial_complexes.SurfaceOfGenus(3)
        Simplicial complex with 15 vertices and 38 facets
        sage: M4 = simplicial_complexes.MooreSpace(4)
        sage: M4.homology()
        {0: 0, 1: C4, 2: 0}
        sage: simplicial_complexes.MatchingComplex(6).homology()
        {0: 0, 1: Z^16, 2: 0}
    """

    def Sphere(self,n):
        """
        A minimal triangulation of the `n`-dimensional sphere.

        INPUT:

        -  ``n`` -- positive integer

        EXAMPLES::

            sage: simplicial_complexes.Sphere(2)
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 2, 3), (0, 1, 2), (1, 2, 3), (0, 1, 3)}
            sage: simplicial_complexes.Sphere(5).homology()
            {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: Z}
            sage: [simplicial_complexes.Sphere(n).euler_characteristic() for n in range(6)]
            [2, 0, 2, 0, 2, 0]
            sage: [simplicial_complexes.Sphere(n).f_vector() for n in range(6)]
            [[1, 2],
             [1, 3, 3],
             [1, 4, 6, 4],
             [1, 5, 10, 10, 5],
             [1, 6, 15, 20, 15, 6],
             [1, 7, 21, 35, 35, 21, 7]]
        """
        S = Simplex(n+1)
        facets = S.faces()
        return SimplicialComplex(facets, is_mutable=False)

    def Simplex(self, n):
        """
        An `n`-dimensional simplex, as a simplicial complex.

        INPUT:

        -  ``n`` -- a non-negative integer

        OUTPUT: the simplicial complex consisting of the `n`-simplex
        on vertices `(0, 1, ..., n)` and all of its faces.

        EXAMPLES::

            sage: simplicial_complexes.Simplex(3)
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 1, 2, 3)}
            sage: simplicial_complexes.Simplex(5).euler_characteristic()
            1
        """
        return SimplicialComplex([Simplex(n)], is_mutable=False)

    def Torus(self):
        r"""
        A minimal triangulation of the torus.

        This is a simplicial complex with 7 vertices, 21 edges and 14
        faces. It is the unique triangulation of the torus with 7
        vertices, and has been found by Möbius in 1861.

        This is also the combinatorial structure of the Császár
        polyhedron (see :wikipedia:`Császár_polyhedron`).

        EXAMPLES::

            sage: T = simplicial_complexes.Torus(); T.homology(1)
            Z x Z
            sage: T.f_vector()
            [1, 7, 21, 14]

        TESTS::

            sage: T.flip_graph().is_isomorphic(graphs.HeawoodGraph())
            True

        REFERENCES:

        .. [LutzCsas] `Császár's Torus <http://www.eg-models.de/models/Classical_Models/2001.02.069/_direct_link.html>`_
        """
        return SimplicialComplex([[0,1,2], [1,2,4], [1,3,4], [1,3,6],
                                  [0,1,5], [1,5,6], [2,3,5], [2,4,5],
                                  [2,3,6], [0,2,6], [0,3,4], [0,3,5],
                                  [4,5,6], [0,4,6]],
                                 is_mutable=False)

    def RealProjectivePlane(self):
        """
        A minimal triangulation of the real projective plane.

        EXAMPLES::

            sage: P = simplicial_complexes.RealProjectivePlane()
            sage: Q = simplicial_complexes.ProjectivePlane()
            sage: P == Q
            True
            sage: P.cohomology(1)
            0
            sage: P.cohomology(2)
            C2
            sage: P.cohomology(1, base_ring=GF(2))
            Vector space of dimension 1 over Finite Field of size 2
            sage: P.cohomology(2, base_ring=GF(2))
            Vector space of dimension 1 over Finite Field of size 2
        """
        return SimplicialComplex([[0,1,2], [0,2,3], [0,1,5], [0,4,5],
                                  [0,3,4], [1,2,4], [1,3,4], [1,3,5],
                                  [2,3,5], [2,4,5]],
                                 is_mutable=False)

    ProjectivePlane = RealProjectivePlane

    def KleinBottle(self):
        """
        A minimal triangulation of the Klein bottle, as presented for example
        in Davide Cervone's thesis [Ce1994]_.

        EXAMPLES::

            sage: simplicial_complexes.KleinBottle()
            Simplicial complex with vertex set (0, 1, 2, 3, 4, 5, 6, 7) and 16 facets

        REFERENCES:

        .. [Ce1994] D. P. Cervone, "Vertex-minimal simplicial immersions of the Klein
           bottle in three-space", Geom. Ded. 50 (1994) 117-141,
           http://www.math.union.edu/~dpvc/papers/1993-03.kb/vmkb.pdf.
        """
        return SimplicialComplex([[2,3,7], [1,2,3], [1,3,5], [1,5,7],
                                  [1,4,7], [2,4,6], [1,2,6], [1,6,0],
                                  [1,4,0], [2,4,0], [3,4,7], [3,4,6],
                                  [3,5,6], [5,6,0], [2,5,0], [2,5,7]],
                                 is_mutable=False)

    def SurfaceOfGenus(self, g, orientable=True):
        """
        A surface of genus `g`.

        INPUT:

        -  ``g`` -- a non-negative integer.  The desired genus

        -  ``orientable`` -- boolean (optional, default ``True``). If
           ``True``, return an orientable surface, and if ``False``,
           return a non-orientable surface.

        In the orientable case, return a sphere if `g` is zero, and
        otherwise return a `g`-fold connected sum of a torus with itself.

        In the non-orientable case, raise an error if `g` is zero.  If
        `g` is positive, return a `g`-fold connected sum of a
        real projective plane with itself.

        EXAMPLES::

            sage: simplicial_complexes.SurfaceOfGenus(2)
            Simplicial complex with 11 vertices and 26 facets
            sage: simplicial_complexes.SurfaceOfGenus(1, orientable=False)
            Simplicial complex with vertex set (0, 1, 2, 3, 4, 5) and 10 facets
        """
        if g == 0:
            if not orientable:
                raise ValueError("No non-orientable surface of genus zero.")
            else:
                return simplicial_complexes.Sphere(2)
        if orientable:
            T = simplicial_complexes.Torus()
        else:
            T = simplicial_complexes.RealProjectivePlane()
        S = T
        for i in range(g-1):
            S = S.connected_sum(T, is_mutable=False)
        return S

    def MooreSpace(self, q):
        """
        Triangulation of the mod `q` Moore space.

        INPUT:

        -  ``q`` -0 integer, at least 2

        This is a simplicial complex with simplices of dimension 0, 1,
        and 2, such that its reduced homology is isomorphic to
        `\\ZZ/q\\ZZ` in dimension 1, zero otherwise.

        If `q=2`, this is the real projective plane.  If `q>2`, then
        construct it as follows: start with a triangle with vertices
        1, 2, 3.  We take a `3q`-gon forming a `q`-fold cover of the
        triangle, and we form the resulting complex as an
        identification space of the `3q`-gon.  To triangulate this
        identification space, put `q` vertices `A_0`, ..., `A_{q-1}`,
        in the interior, each of which is connected to 1, 2, 3 (two
        facets each: `[1, 2, A_i]`, `[2, 3, A_i]`).  Put `q` more
        vertices in the interior: `B_0`, ..., `B_{q-1}`, with facets
        `[3, 1, B_i]`, `[3, B_i, A_i]`, `[1, B_i, A_{i+1}]`, `[B_i,
        A_i, A_{i+1}]`.  Then triangulate the interior polygon with
        vertices `A_0`, `A_1`, ..., `A_{q-1}`.

        EXAMPLES::

            sage: simplicial_complexes.MooreSpace(2)
            Simplicial complex with vertex set (0, 1, 2, 3, 4, 5) and 10 facets
            sage: simplicial_complexes.MooreSpace(3).homology()[1]
            C3
            sage: simplicial_complexes.MooreSpace(4).suspension().homology()[2]
            C4
            sage: simplicial_complexes.MooreSpace(8)
            Simplicial complex with 19 vertices and 54 facets
        """
        if q <= 1:
            raise ValueError("The mod q Moore space is only defined if q is at least 2")
        if q == 2:
            return simplicial_complexes.RealProjectivePlane()
        facets = []
        for i in range(q):
            Ai = "A" + str(i)
            Aiplus = "A" + str((i+1)%q)
            Bi = "B" + str(i)
            facets.append([1, 2, Ai])
            facets.append([2, 3, Ai])
            facets.append([3, 1, Bi])
            facets.append([3, Bi, Ai])
            facets.append([1, Bi, Aiplus])
            facets.append([Bi, Ai, Aiplus])
        for i in range(1, q-1):
            Ai = "A" + str(i)
            Aiplus = "A" + str((i+1)%q)
            facets.append(["A0", Ai, Aiplus])
        return SimplicialComplex(facets, is_immutable=True)

    def ComplexProjectivePlane(self):
        """
        A minimal triangulation of the complex projective plane.

        This was constructed by Kühnel and Banchoff [KB1983]_.

        REFERENCES:

        .. [KB1983] W. Kühnel and T. F. Banchoff, "The 9-vertex complex
           projective plane", Math. Intelligencer 5 (1983), no. 3,
           11-22.

        EXAMPLES::

            sage: C = simplicial_complexes.ComplexProjectivePlane()
            sage: C.f_vector()
            [1, 9, 36, 84, 90, 36]
            sage: C.homology(2)
            Z
            sage: C.homology(4)
            Z
        """
        return SimplicialComplex(
            [[1, 2, 4, 5, 6], [2, 3, 5, 6, 4], [3, 1, 6, 4, 5],
             [1, 2, 4, 5, 9], [2, 3, 5, 6, 7], [3, 1, 6, 4, 8],
             [2, 3, 6, 4, 9], [3, 1, 4, 5, 7], [1, 2, 5, 6, 8],
             [3, 1, 5, 6, 9], [1, 2, 6, 4, 7], [2, 3, 4, 5, 8],
             [4, 5, 7, 8, 9], [5, 6, 8, 9, 7], [6, 4, 9, 7, 8],
             [4, 5, 7, 8, 3], [5, 6, 8, 9, 1], [6, 4, 9, 7, 2],
             [5, 6, 9, 7, 3], [6, 4, 7, 8, 1], [4, 5, 8, 9, 2],
             [6, 4, 8, 9, 3], [4, 5, 9, 7, 1], [5, 6, 7, 8, 2],
             [7, 8, 1, 2, 3], [8, 9, 2, 3, 1], [9, 7, 3, 1, 2],
             [7, 8, 1, 2, 6], [8, 9, 2, 3, 4], [9, 7, 3, 1, 5],
             [8, 9, 3, 1, 6], [9, 7, 1, 2, 4], [7, 8, 2, 3, 5],
             [9, 7, 2, 3, 6], [7, 8, 3, 1, 4], [8, 9, 1, 2, 5]],
            is_mutable=False)

    def PseudoQuaternionicProjectivePlane(self):
        r"""
        Returns a pure simplicial complex of dimension 8 with 490 facets.

        .. WARNING::

            This is expected to be a triangulation of the projective plane
            `HP^2` over the ring of quaternions, but this has not been
            proved yet.

        This simplicial complex has the same homology as `HP^2`. Its
        automorphism group is isomorphic to the alternating group `A_5`
        and acts transitively on vertices.

        This is defined here using the description in [BrK92]_. This
        article deals with three different triangulations. This procedure
        returns the only one which has a transitive group of
        automorphisms.

        EXAMPLES::

            sage: HP2 = simplicial_complexes.PseudoQuaternionicProjectivePlane() ; HP2
            Simplicial complex with 15 vertices and 490 facets
            sage: HP2.f_vector()
            [1, 15, 105, 455, 1365, 3003, 4515, 4230, 2205, 490]

        Checking its automorphism group::

            sage: HP2.automorphism_group().is_isomorphic(AlternatingGroup(5))
            True

        REFERENCES:

        .. [BrK92] Brehm U., Kuhnel W., "15-vertex triangulations of an
                   8-manifold", Math. Annalen 294 (1992), no. 1, 167-193.
        """
        from sage.groups.perm_gps.permgroup import PermutationGroup
        P = [(1,2,3,4,5),(6,7,8,9,10),(11,12,13,14,15)]
        S = [(1,6,11),(2,15,14),(3,13,8),(4,7,5),(9,12,10)]
        start_list = [
            (1,2,3,6,8,11,13,14,15),    # A
            (1,3,6,8,9,10,11,12,13),    # B
            (1,2,6,9,10,11,12,14,15),   # C
            (1,2,3,4,7,9,12,14,15),     # D
            (1,2,4,7,9,10,12,13,14),    # E
            (1,2,6,8,9,10,11,14,15),    # F
            (1,2,3,4,5,6,9,11,13),      # G
            (1,3,5,6,8,9,10,11,12),     # H
            (1,3,5,6,7,8,9,10,11),      # I
            (1,2,3,4,5,7,10,12,15),     # J
            (1,2,3,7,8,10,12,13,14),    # K
            (2,5,6,7,8,9,10,13,14),     # M

            (3,4,6,7,11,12,13,14,15),   # L
            (3,4,6,7,10,12,13,14,15)]   # N
        return SimplicialComplex([ [g(index) for index in tuple]
                for tuple in start_list
                for g in PermutationGroup([P,S]) ])

    def PoincareHomologyThreeSphere(self):
        """
        A triangulation of the Poincare homology 3-sphere.

        This is a manifold whose integral homology is identical to the
        ordinary 3-sphere, but it is not simply connected. In particular,
        its fundamental group is the binary icosahedral group, which has
        order 120. The triangulation given here has 16 vertices and is
        due to Björner and Lutz [BL2000]_.

        REFERENCES:

        .. [BL2000] Anders Björner and Frank H. Lutz, "Simplicial
           manifolds, bistellar flips and a 16-vertex triangulation of
           the Poincaré homology 3-sphere", Experiment. Math. 9
           (2000), no. 2, 275-289.

        EXAMPLES::

            sage: S3 = simplicial_complexes.Sphere(3)
            sage: Sigma3 = simplicial_complexes.PoincareHomologyThreeSphere()
            sage: S3.homology() == Sigma3.homology()
            True
            sage: Sigma3.fundamental_group().cardinality() # long time
            120
        """
        return SimplicialComplex(
            [[1, 2, 4, 9], [1, 2, 4, 15], [1, 2, 6, 14], [1, 2, 6, 15],
             [1, 2, 9, 14], [1, 3, 4, 12], [1, 3, 4, 15], [1, 3, 7, 10],
             [1, 3, 7, 12], [1, 3, 10, 15], [1, 4, 9, 12], [1, 5, 6, 13],
             [1, 5, 6, 14], [1, 5, 8, 11], [1, 5, 8, 13], [1, 5, 11, 14],
             [1, 6, 13, 15], [1, 7, 8, 10], [1, 7, 8, 11], [1, 7, 11, 12],
             [1, 8, 10, 13], [1, 9, 11, 12], [1, 9, 11, 14], [1, 10, 13, 15],
             [2, 3, 5, 10], [2, 3, 5, 11], [2, 3, 7, 10], [2, 3, 7, 13],
             [2, 3, 11, 13], [2, 4, 9, 13], [2, 4, 11, 13], [2, 4, 11, 15],
             [2, 5, 8, 11], [2, 5, 8, 12], [2, 5, 10, 12], [2, 6, 10, 12],
             [2, 6, 10, 14], [2, 6, 12, 15], [2, 7, 9, 13], [2, 7, 9, 14],
             [2, 7, 10, 14], [2, 8, 11, 15], [2, 8, 12, 15], [3, 4, 5, 14],
             [3, 4, 5, 15], [3, 4, 12, 14], [3, 5, 10, 15], [3, 5, 11, 14],
             [3, 7, 12, 13], [3, 11, 13, 14], [3, 12, 13, 14], [4, 5, 6, 7],
             [4, 5, 6, 14], [4, 5, 7, 15], [4, 6, 7, 11], [4, 6, 10, 11],
             [4, 6, 10, 14], [4, 7, 11, 15], [4, 8, 9, 12], [4, 8, 9, 13],
             [4, 8, 10, 13], [4, 8, 10, 14], [4, 8, 12, 14], [4, 10, 11, 13],
             [5, 6, 7, 13], [5, 7, 9, 13], [5, 7, 9, 15], [5, 8, 9, 12],
             [5, 8, 9, 13], [5, 9, 10, 12], [5, 9, 10, 15], [6, 7, 11, 12],
             [6, 7, 12, 13], [6, 10, 11, 12], [6, 12, 13, 15], [7, 8, 10, 14],
             [7, 8, 11, 15], [7, 8, 14, 15], [7, 9, 14, 15], [8, 12, 14, 15],
             [9, 10, 11, 12], [9, 10, 11, 16], [9, 10, 15, 16], [9, 11, 14, 16],
             [9, 14, 15, 16], [10, 11, 13, 16], [10, 13, 15, 16],
             [11, 13, 14, 16], [12, 13, 14, 15], [13, 14, 15, 16]],
            is_mutable=False)

    def RealProjectiveSpace(self, n):
        r"""
        A triangulation of `\Bold{R}P^n` for any `n \geq 0`.

        INPUT:

        - ``n`` -- integer, the dimension of the real projective space
          to construct

        The first few cases are pretty trivial:

        - `\Bold{R}P^0` is a point.

        - `\Bold{R}P^1` is a circle, triangulated as the boundary of a
          single 2-simplex.

        - `\Bold{R}P^2` is the real projective plane, here given its
          minimal triangulation with 6 vertices, 15 edges, and 10
          triangles.

        - `\Bold{R}P^3`: any triangulation has at least 11 vertices by
          a result of Walkup [Wa1970]_; this function returns a
          triangulation with 11 vertices, as given by Lutz [Lu2005]_.

        - `\Bold{R}P^4`: any triangulation has at least 16 vertices by
          a result of Walkup; this function returns a triangulation
          with 16 vertices as given by Lutz; see also Datta [Da2007]_,
          Example 3.12.

        - `\Bold{R}P^n`: Lutz has found a triangulation of
          `\Bold{R}P^5` with 24 vertices, but it does not seem to have
          been published.  Kühnel [Ku1987]_ has described a triangulation of
          `\Bold{R}P^n`, in general, with `2^{n+1}-1` vertices; see
          also Datta, Example 3.21.  This triangulation is presumably
          not minimal, but it seems to be the best in the published
          literature as of this writing.  So this function returns it
          when `n > 4`.

        ALGORITHM: For `n < 4`, these are constructed explicitly by
        listing the facets.  For `n = 4`, this is constructed by
        specifying 16 vertices, two facets, and a certain subgroup `G`
        of the symmetric group `S_{16}`.  Then the set of all facets
        is the `G`-orbit of the two given facets.  This is implemented
        here by explicitly listing all of the facets; the facets
        can be computed by the function :func:`facets_for_RP4`, but
        running the function takes a few seconds.

        For `n > 4`, the construction is as follows: let `S` denote
        the simplicial complex structure on the `n`-sphere given by
        the first barycentric subdivision of the boundary of an
        `(n+1)`-simplex.  This has a simplicial antipodal action: if
        `V` denotes the vertices in the boundary of the simplex, then
        the vertices in its barycentric subdivision `S` correspond to
        nonempty proper subsets `U` of `V`, and the antipodal action
        sends any subset `U` to its complement.  One can show that
        modding out by this action results in a triangulation for
        `\Bold{R}P^n`.  To find the facets in this triangulation, find
        the facets in `S`.  These are indentified in pairs to form
        `\Bold{R}P^n`, so choose a representative from each pair: for
        each facet in `S`, replace any vertex in `S` containing 0 with
        its complement.

        Of course these complexes increase in size pretty quickly as
        `n` increases.

        REFERENCES:

        .. [Da2007] Basudeb Datta, "Minimal triangulations of manifolds",
           J. Indian Inst. Sci. 87 (2007), no. 4, 429-449.

        .. [Ku1987] W. Kühnel, "Minimal triangulations of Kummer varieties",
           Abh. Math. Sem. Univ. Hamburg 57 (1987), 7-20.

        .. [Lu2005] Frank H. Lutz, "Triangulated Manifolds with Few Vertices:
           Combinatorial Manifolds", preprint (2005),
           arXiv:math/0506372.

        .. [Wa1970] David W. Walkup, "The lower bound conjecture for 3- and
           4-manifolds", Acta Math. 125 (1970), 75-107.

        EXAMPLES::

            sage: P3 = simplicial_complexes.RealProjectiveSpace(3)
            sage: P3.f_vector()
            [1, 11, 51, 80, 40]
            sage: P3.homology()
            {0: 0, 1: C2, 2: 0, 3: Z}
            sage: P4 = simplicial_complexes.RealProjectiveSpace(4)
            sage: P4.f_vector()
            [1, 16, 120, 330, 375, 150]
            sage: P4.homology() # long time
            {0: 0, 1: C2, 2: 0, 3: C2, 4: 0}
            sage: P5 = simplicial_complexes.RealProjectiveSpace(5)  # long time (44s on sage.math, 2012)
            sage: P5.f_vector()  # long time
            [1, 63, 903, 4200, 8400, 7560, 2520]

        The following computation can take a long time -- over half an
        hour -- with Sage's default computation of homology groups,
        but if you have CHomP installed, Sage will use that and the
        computation should only take a second or two.  (You can
        download CHomP from http://chomp.rutgers.edu/, or you can
        install it as a Sage package using ``sage -i chomp``). ::

            sage: P5.homology()  # long time # optional - CHomP
            {0: 0, 1: C2, 2: 0, 3: C2, 4: 0, 5: Z}
            sage: simplicial_complexes.RealProjectiveSpace(2).dimension()
            2
            sage: P3.dimension()
            3
            sage: P4.dimension() # long time
            4
            sage: P5.dimension() # long time
            5
        """
        if n == 0:
            return self.Simplex(0)
        if n == 1:
            return self.Sphere(1)
        if n == 2:
            return self.RealProjectivePlane()
        if n == 3:
            # Minimal triangulation found by Walkup and given
            # explicitly by Lutz
            return SimplicialComplex(
                [[1, 2, 3, 7], [1, 4, 7, 9], [2, 3, 4, 8], [2, 5, 8, 10],
                 [3, 6, 7, 10], [1, 2, 3, 11], [1, 4, 7, 10], [2, 3, 4, 11],
                 [2, 5, 9, 10], [3, 6, 8, 9], [1, 2, 6, 9], [1, 4, 8, 9],
                 [2, 3, 7, 8], [2, 6, 9, 10], [3, 6, 9, 10], [1, 2, 6, 11],
                 [1, 4, 8, 10], [2, 4, 6, 10], [3, 4, 5, 9], [4, 5, 6, 7],
                 [1, 2, 7, 9], [1, 5, 6, 8], [2, 4, 6, 11], [3, 4, 5, 11],
                 [4, 5, 6, 11], [1, 3, 5, 10], [1, 5, 6, 11], [2, 4, 8, 10],
                 [3, 4, 8, 9], [4, 5, 7, 9], [1, 3, 5, 11], [1, 5, 8, 10],
                 [2, 5, 7, 8], [3, 5, 9, 10], [4, 6, 7, 10], [1, 3, 7, 10],
                 [1, 6, 8, 9], [2, 5, 7, 9], [3, 6, 7, 8], [5, 6, 7, 8]],
                is_mutable=False)
        if n == 4:
            return SimplicialComplex(
                [(1, 3, 8, 12, 13), (2, 7, 8, 13, 16), (4, 8, 9, 12, 14),
                 (2, 6, 10, 12, 16), (5, 7, 9, 10, 13), (1, 2, 7, 8, 15),
                    (1, 3, 9, 11, 16), (5, 6, 8, 13, 16), (1, 3, 8, 11, 13),
                    (3, 4, 10, 13, 15), (4, 6, 9, 12, 15), (2, 4, 6, 11, 13),
                    (2, 3, 9, 12, 16), (1, 6, 9, 12, 15), (2, 5, 10, 11, 12),
                    (1, 7, 8, 12, 15), (2, 6, 9, 13, 16), (1, 5, 9, 11, 15),
                    (4, 9, 10, 13, 14), (2, 7, 8, 15, 16), (2, 3, 9, 12, 14),
                    (1, 6, 7, 10, 14), (2, 5, 10, 11, 15), (1, 2, 4, 13, 14),
                    (1, 6, 10, 14, 16), (2, 6, 9, 12, 16), (1, 3, 9, 12, 16),
                    (4, 5, 7, 11, 16), (5, 9, 10, 11, 15), (3, 5, 8, 12, 14),
                    (5, 6, 9, 13, 16), (5, 6, 9, 13, 15), (1, 3, 4, 10, 16),
                    (1, 6, 10, 12, 16), (2, 4, 6, 9, 13), (2, 4, 6, 9, 12),
                    (1, 2, 4, 11, 13), (7, 9, 10, 13, 14), (1, 7, 8, 12, 13),
                    (4, 6, 7, 11, 12), (3, 4, 6, 11, 13), (1, 5, 6, 9, 15),
                    (1, 6, 7, 14, 15), (2, 3, 7, 14, 15), (2, 6, 10, 11, 12),
                    (5, 7, 9, 10, 11), (1, 2, 4, 5, 14), (3, 5, 10, 13, 15),
                    (3, 8, 9, 12, 14), (5, 9, 10, 13, 15), (2, 6, 8, 13, 16),
                    (1, 2, 7, 13, 14), (1, 7, 10, 12, 13), (3, 4, 6, 13, 15),
                    (4, 9, 10, 13, 15), (2, 3, 10, 12, 16), (1, 2, 5, 14, 15),
                    (2, 6, 8, 10, 11), (1, 3, 10, 12, 13), (4, 8, 9, 12, 15),
                    (1, 3, 8, 9, 11), (4, 6, 7, 12, 15), (1, 8, 9, 11, 15),
                    (4, 5, 8, 14, 16), (1, 2, 8, 11, 13), (3, 6, 8, 11, 13),
                    (3, 6, 8, 11, 14), (3, 5, 8, 12, 13), (3, 7, 9, 11, 14),
                    (4, 6, 9, 13, 15), (2, 3, 5, 10, 12), (4, 7, 8, 15, 16),
                    (1, 2, 7, 14, 15), (3, 7, 9, 11, 16), (3, 6, 7, 14, 15),
                    (2, 6, 8, 11, 13), (4, 8, 9, 10, 14), (1, 4, 10, 13, 14),
                    (4, 8, 9, 10, 15), (2, 7, 9, 13, 16), (1, 6, 9, 12, 16),
                    (2, 3, 7, 9, 14), (4, 8, 10, 15, 16), (1, 5, 9, 11, 16),
                    (1, 5, 6, 14, 15), (5, 7, 9, 11, 16), (4, 5, 7, 11, 12),
                    (5, 7, 10, 11, 12), (2, 3, 10, 15, 16), (1, 2, 7, 8, 13),
                    (1, 6, 7, 10, 12), (1, 3, 10, 12, 16), (7, 9, 10, 11, 14),
                    (1, 7, 10, 13, 14), (1, 2, 4, 5, 11), (3, 4, 6, 7, 11),
                    (1, 6, 7, 12, 15), (1, 3, 4, 10, 13), (1, 4, 10, 14, 16),
                    (2, 4, 6, 11, 12), (5, 6, 8, 14, 16), (3, 5, 6, 8, 13),
                    (3, 5, 6, 8, 14), (1, 2, 8, 11, 15), (1, 4, 5, 14, 16),
                    (2, 3, 7, 15, 16), (8, 9, 10, 11, 14), (1, 3, 4, 11, 16),
                    (6, 8, 10, 14, 16), (8, 9, 10, 11, 15), (1, 3, 4, 11, 13),
                    (2, 4, 5, 12, 14), (2, 4, 9, 13, 14), (3, 4, 7, 11, 16),
                    (3, 6, 7, 11, 14), (3, 8, 9, 11, 14), (2, 8, 10, 11, 15),
                    (1, 3, 8, 9, 12), (4, 5, 7, 8, 16), (4, 5, 8, 12, 14),
                    (2, 4, 9, 12, 14), (6, 8, 10, 11, 14), (3, 5, 6, 13, 15),
                    (1, 4, 5, 11, 16), (3, 5, 6, 14, 15), (2, 4, 5, 11, 12),
                    (4, 5, 7, 8, 12), (1, 8, 9, 12, 15), (5, 7, 8, 13, 16),
                    (2, 3, 5, 12, 14), (3, 5, 10, 12, 13), (6, 7, 10, 11, 12),
                    (5, 7, 9, 13, 16), (6, 7, 10, 11, 14), (5, 7, 10, 12, 13),
                    (1, 2, 5, 11, 15), (1, 5, 6, 9, 16), (5, 7, 8, 12, 13),
                    (4, 7, 8, 12, 15), (2, 3, 5, 10, 15), (2, 6, 8, 10, 16),
                    (3, 4, 10, 15, 16), (1, 5, 6, 14, 16), (2, 3, 5, 14, 15),
                    (2, 3, 7, 9, 16), (2, 7, 9, 13, 14), (3, 4, 6, 7, 15),
                    (4, 8, 10, 14, 16), (3, 4, 7, 15, 16), (2, 8, 10, 15, 16)],
                 is_mutable=False)
        if n >= 5:
            # Use the construction given by Datta in Example 3.21.
            V = set(range(0, n+2))
            S = simplicial_complexes.Sphere(n).barycentric_subdivision()
            X = S.facets()
            facets = set([])
            for f in X:
                new = []
                for v in f:
                    if 0 in v:
                        new.append(tuple(V.difference(v)))
                    else:
                        new.append(v)
                facets.add(tuple(new))
            return SimplicialComplex(list(facets), is_mutable=False)

    def K3Surface(self):
        """
        Returns a minimal triangulation of the K3 surface.

        This is a pure simplicial complex of dimension 4 with 16 vertices
        and 288 facets. It was constructed by Casella and Kühnel
        in [CK2001]_. The construction here uses the labeling from
        Spreer and Kühnel [SK2011]_.

        REFERENCES:

        .. [CK2001] M. Casella and W. Kühnel, "A triangulated K3 surface
           with the minimum number of vertices", Topology 40 (2001),
           753–772.

        .. [SK2011] J. Spreer and W. Kühnel, "Combinatorial properties
           of the K3 surface: Simplicial blowups and slicings", Experimental
           Mathematics, Volume 20, Issue 2, 2011.

        EXAMPLES::

            sage: K3=simplicial_complexes.K3Surface() ; K3
            Simplicial complex with 16 vertices and 288 facets
            sage: K3.f_vector()
            [1, 16, 120, 560, 720, 288]

        This simplicial complex is implemented just by listing all 288
        facets. The list of facets can be computed by the function
        :func:`facets_for_K3`, but running the function takes a few
        seconds.
        """
        return SimplicialComplex(
            [(2, 10, 13, 15, 16), (2, 8, 11, 15, 16), (2, 5, 7, 8, 10),
             (1, 9, 11, 13, 14), (1, 2, 8, 10, 12), (1, 3, 5, 6, 11),
                (1, 5, 6, 9, 12), (1, 2, 6, 13, 16), (1, 4, 10, 13, 14),
                (1, 9, 10, 14, 15), (2, 4, 7, 8, 12), (3, 4, 6, 10, 12),
                (1, 6, 7, 8, 9), (3, 4, 5, 7, 15), (1, 7, 12, 15, 16),
                (4, 5, 7, 13, 16), (5, 8, 11, 12, 15), (2, 4, 7, 12, 14),
                (1, 4, 5, 14, 16), (2, 5, 6, 10, 11), (1, 6, 8, 12, 14),
                (5, 8, 9, 14, 16), (5, 10, 11, 12, 13), (2, 4, 8, 9, 12),
                (7, 9, 12, 15, 16), (1, 2, 6, 9, 15), (1, 5, 14, 15, 16),
                (2, 3, 4, 5, 9), (6, 8, 10, 11, 15), (1, 5, 8, 10, 12),
                (1, 3, 7, 9, 10), (6, 7, 8, 9, 13), (1, 2, 9, 11, 15),
                (2, 8, 11, 14, 16), (2, 4, 5, 13, 16), (1, 4, 8, 13, 15),
                (4, 7, 8, 10, 11), (2, 3, 9, 11, 14), (2, 3, 4, 9, 13),
                (2, 8, 10, 12, 13), (1, 2, 4, 11, 15), (2, 3, 9, 11, 15),
                (3, 5, 10, 13, 15), (3, 4, 5, 9, 11), (6, 10, 13, 15, 16),
                (8, 10, 11, 15, 16), (6, 7, 11, 13, 15), (1, 5, 7, 15, 16),
                (4, 5, 7, 9, 15), (3, 4, 6, 7, 16), (2, 3, 11, 14, 16),
                (3, 4, 9, 11, 13), (1, 2, 5, 14, 15), (2, 3, 9, 13, 14),
                (1, 2, 5, 13, 16), (2, 3, 7, 8, 12), (2, 9, 11, 12, 14),
                (1, 9, 11, 15, 16), (4, 6, 9, 14, 16), (1, 4, 9, 13, 14),
                (1, 2, 3, 12, 16), (8, 11, 12, 14, 15), (2, 4, 11, 12, 14),
                (1, 4, 10, 12, 13), (1, 2, 6, 7, 13), (1, 3, 6, 10, 11),
                (1, 6, 8, 9, 12), (1, 4, 5, 6, 14), (3, 9, 10, 12, 15),
                (5, 8, 11, 12, 16), (5, 9, 10, 14, 15), (3, 9, 12, 15, 16),
                (3, 6, 8, 14, 15), (2, 4, 9, 10, 16), (5, 8, 9, 13, 15),
                (2, 3, 6, 9, 15), (6, 11, 12, 14, 16), (2, 3, 10, 13, 15),
                (2, 8, 9, 10, 13), (3, 4, 8, 11, 13), (3, 4, 5, 7, 13),
                (5, 7, 8, 10, 14), (4, 12, 13, 14, 15), (6, 7, 10, 14, 16),
                (5, 10, 11, 13, 14), (3, 4, 7, 13, 16), (6, 8, 9, 12, 13),
                (1, 3, 4, 10, 14), (2, 4, 6, 11, 12), (1, 7, 9, 10, 14),
                (4, 6, 8, 13, 14), (4, 9, 10, 11, 16), (3, 7, 8, 10, 16),
                (5, 7, 9, 15, 16), (1, 7, 9, 11, 14), (6, 8, 10, 15, 16),
                (5, 8, 9, 10, 14), (7, 8, 10, 14, 16), (2, 6, 7, 9, 11),
                (7, 9, 10, 13, 15), (3, 6, 7, 10, 12), (2, 4, 6, 10, 11),
                (4, 5, 8, 9, 11), (1, 2, 3, 8, 16), (3, 7, 9, 10, 12),
                (1, 2, 6, 8, 14), (3, 5, 6, 13, 15), (1, 5, 6, 12, 14),
                (2, 5, 7, 14, 15), (1, 5, 10, 11, 12), (3, 7, 8, 10, 11),
                (1, 2, 6, 14, 15), (1, 2, 6, 8, 16), (7, 9, 10, 12, 15),
                (3, 4, 6, 8, 14), (3, 7, 13, 14, 16), (2, 5, 7, 8, 14),
                (6, 7, 9, 10, 14), (2, 3, 7, 12, 14), (4, 10, 12, 13, 14),
                (2, 5, 6, 11, 13), (4, 5, 6, 7, 16), (1, 3, 12, 13, 16),
                (1, 4, 11, 15, 16), (1, 3, 4, 6, 10), (1, 10, 11, 12, 13),
                (6, 9, 11, 12, 14), (1, 4, 7, 8, 15), (5, 8, 9, 10, 13),
                (1, 2, 5, 7, 15), (1, 7, 12, 13, 16), (3, 11, 13, 14, 16),
                (1, 2, 5, 7, 13), (4, 7, 8, 9, 15), (1, 5, 6, 10, 11),
                (6, 7, 10, 13, 15), (3, 4, 7, 14, 15), (7, 11, 13, 14, 16),
                (3, 4, 10, 12, 14), (3, 6, 8, 10, 16), (2, 7, 8, 14, 16),
                (2, 3, 4, 5, 13), (5, 8, 12, 13, 15), (4, 6, 9, 13, 14),
                (2, 4, 5, 6, 12), (1, 3, 7, 8, 9), (8, 11, 12, 14, 16),
                (1, 7, 12, 13, 15), (8, 12, 13, 14, 15), (2, 8, 9, 12, 13),
                (4, 6, 10, 12, 15), (2, 8, 11, 14, 15), (2, 6, 9, 11, 12),
                (8, 9, 10, 11, 16), (2, 3, 6, 13, 15), (2, 3, 12, 15, 16),
                (1, 3, 5, 9, 12), (2, 5, 6, 9, 12), (2, 10, 12, 13, 14),
                (2, 6, 13, 15, 16), (2, 3, 11, 15, 16), (3, 5, 6, 8, 15),
                (2, 4, 5, 9, 12), (5, 6, 8, 11, 15), (6, 8, 12, 13, 14),
                (1, 2, 3, 8, 12), (1, 4, 7, 8, 11), (3, 5, 7, 14, 15),
                (3, 5, 7, 13, 14), (1, 7, 10, 11, 14), (6, 7, 11, 12, 15),
                (3, 4, 6, 7, 12), (1, 2, 4, 7, 11), (6, 9, 10, 14, 16),
                (4, 10, 12, 15, 16), (5, 6, 7, 12, 16), (3, 9, 11, 13, 14),
                (5, 9, 14, 15, 16), (4, 5, 6, 7, 12), (1, 3, 9, 10, 15),
                (4, 7, 8, 9, 12), (5, 9, 10, 13, 15), (1, 3, 8, 13, 16),
                (2, 9, 12, 13, 14), (6, 7, 10, 12, 15), (2, 6, 8, 14, 15),
                (3, 5, 6, 8, 11), (3, 4, 7, 12, 14), (1, 3, 10, 14, 15),
                (7, 11, 12, 13, 16), (3, 11, 12, 13, 16), (3, 4, 5, 8, 15),
                (2, 4, 7, 8, 10), (2, 4, 7, 14, 15), (1, 2, 10, 12, 16),
                (1, 6, 8, 13, 16), (1, 7, 8, 13, 15), (3, 9, 11, 15, 16),
                (4, 6, 10, 11, 15), (2, 4, 11, 14, 15), (1, 3, 8, 9, 12),
                (1, 3, 6, 14, 15), (2, 4, 5, 6, 10), (1, 4, 9, 14, 16),
                (5, 7, 9, 12, 16), (1, 3, 7, 10, 11), (7, 8, 9, 13, 15),
                (3, 5, 10, 14, 15), (1, 4, 10, 12, 16), (3, 4, 5, 8, 11),
                (1, 2, 6, 7, 9), (1, 3, 11, 12, 13), (1, 5, 7, 13, 16),
                (5, 7, 10, 11, 14), (2, 10, 12, 15, 16), (3, 6, 7, 10, 16),
                (1, 2, 5, 8, 10), (4, 10, 11, 15, 16), (5, 8, 10, 12, 13),
                (3, 6, 8, 10, 11), (4, 5, 7, 9, 12), (6, 7, 11, 12, 16),
                (3, 5, 9, 11, 16), (8, 9, 10, 14, 16), (3, 4, 6, 8, 16),
                (1, 10, 11, 13, 14), (2, 9, 10, 13, 16), (1, 2, 5, 8, 14),
                (2, 4, 5, 10, 16), (1, 2, 7, 9, 11), (1, 3, 5, 6, 9),
                (5, 7, 11, 13, 14), (3, 5, 10, 13, 14), (2, 4, 8, 9, 10),
                (4, 11, 12, 14, 15), (2, 3, 7, 14, 16), (3, 4, 8, 13, 16),
                (6, 7, 9, 11, 14), (5, 6, 11, 13, 15), (4, 5, 6, 14, 16),
                (3, 4, 8, 14, 15), (4, 5, 8, 9, 15), (1, 4, 8, 11, 13),
                (5, 6, 12, 14, 16), (2, 3, 10, 12, 14), (1, 2, 5, 10, 16),
                (2, 5, 7, 10, 11), (2, 6, 7, 11, 13), (1, 4, 5, 10, 16),
                (2, 6, 8, 15, 16), (2, 3, 10, 12, 15), (7, 11, 12, 13, 15),
                (1, 3, 8, 11, 13), (4, 8, 9, 10, 11), (1, 9, 14, 15, 16),
                (1, 3, 6, 9, 15), (6, 9, 12, 13, 14), (2, 3, 10, 13, 14),
                (2, 5, 7, 11, 13), (2, 3, 5, 6, 13), (4, 6, 8, 13, 16),
                (6, 7, 9, 10, 13), (5, 8, 12, 14, 16), (4, 6, 9, 13, 16),
                (5, 8, 9, 11, 16), (2, 3, 5, 6, 9), (1, 3, 5, 11, 12),
                (3, 7, 8, 9, 12), (4, 6, 11, 12, 15), (3, 5, 9, 12, 16),
                (5, 11, 12, 13, 15), (1, 3, 4, 6, 14), (3, 5, 11, 12, 16),
                (1, 5, 8, 12, 14), (4, 8, 13, 14, 15), (1, 3, 7, 8, 11),
                (6, 9, 10, 13, 16), (2, 4, 9, 13, 16), (1, 6, 7, 8, 13),
                (1, 4, 12, 13, 15), (2, 4, 7, 10, 11), (1, 4, 9, 11, 13),
                (6, 7, 11, 14, 16), (1, 4, 9, 11, 16), (1, 4, 12, 15, 16),
                (1, 2, 4, 7, 15), (2, 3, 7, 8, 16), (1, 4, 5, 6, 10)],
             is_mutable=False)

    def BarnetteSphere(self):
        r"""
        Returns Barnette's triangulation of the 3-sphere.

        This is a pure simplicial complex of dimension 3 with 8
        vertices and 19 facets, which is a non-polytopal triangulation
        of the 3-sphere. It was constructed by Barnette in
        [B1970]_. The construction here uses the labeling from De
        Loera, Rambau and Santos [DLRS2010]_. Another reference is chapter
        III.4 of Ewald [E1996]_.

        EXAMPLES::

            sage: BS = simplicial_complexes.BarnetteSphere() ; BS
            Simplicial complex with vertex set (1, 2, 3, 4, 5, 6, 7, 8) and 19 facets
            sage: BS.f_vector()
            [1, 8, 27, 38, 19]

        TESTS:

        Checks that this is indeed the same Barnette Sphere as the one
        given on page 87 of [E1996]_.::

            sage: BS2 = SimplicialComplex([[1,2,3,4],[3,4,5,6],[1,2,5,6],
            ...                            [1,2,4,7],[1,3,4,7],[3,4,6,7],
            ...                            [3,5,6,7],[1,2,5,7],[2,5,6,7],
            ...                            [2,4,6,7],[1,2,3,8],[2,3,4,8],
            ...                            [3,4,5,8],[4,5,6,8],[1,2,6,8],
            ...                            [1,5,6,8],[1,3,5,8],[2,4,6,8],
            ...                            [1,3,5,7]])
            sage: BS.is_isomorphic(BS2)
            True

        REFERENCES:

        .. [B1970] Barnette, "Diagrams and Schlegel diagrams", in
           Combinatorial Structures and Their Applications, Proc. Calgary
           Internat. Conference 1969, New York, 1970, Gordon and Breach.

        .. [DLRS2010] De Loera, Rambau and Santos, "Triangulations:
           Structures for Algorithms and Applications", Algorithms and
           Computation in Mathematics, Volume 25, Springer, 2011.

        .. [E1996] Ewald, "Combinatorial Convexity and Algebraic Geometry",
           vol. 168 of Graduate Texts in Mathematics, Springer, 1996

        """
        return SimplicialComplex([
                (1,2,4,5),(2,3,5,6),(1,3,4,6),(1,2,3,7),(4,5,6,7),(1,2,4,7),
                (2,4,5,7),(2,3,5,7),(3,5,6,7),(3,1,6,7),(1,6,4,7),(1,2,3,8),
                (4,5,6,8),(1,2,5,8),(1,4,5,8),(2,3,6,8),(2,5,6,8),(3,1,4,8),
                (3,6,4,8)])

    def BrucknerGrunbaumSphere(self):
        r"""
        Returns Bruckner and Grunbaum's triangulation of the 3-sphere.

        This is a pure simplicial complex of dimension 3 with 8
        vertices and 20 facets, which is a non-polytopal triangulation
        of the 3-sphere. It appeared first in [Br1910]_ and was studied in
        [GrS1967]_.

        It is defined here as the link of any vertex in the unique minimal
        triangulation of the complex projective plane, see chapter 4 of
        [Ku1995]_.

        EXAMPLES::

            sage: BGS = simplicial_complexes.BrucknerGrunbaumSphere() ; BGS
            Simplicial complex with vertex set (1, 2, 3, 4, 5, 6, 7, 8) and 20 facets
            sage: BGS.f_vector()
            [1, 8, 28, 40, 20]

        REFERENCES:

        .. [Br1910] Bruckner, "Uber die Ableitung der allgemeinen
           Polytope und die nach Isomorphismus verschiedenen Typen der
           allgemeinen Achtzelle (Oktatope)", Verhand. Konik. Akad. Wetenschap,
           Erste Sectie, 10 (1910)

        .. [GrS1967] Grunbaum and Sreedharan, "An enumeration of simplicial
           4-polytopes with 8 vertices", J. Comb. Th. 2, 437-465 (1967)

        .. [Ku1995] Kuhnel, "Tight Polyhedral Submanifolds and Tight Triangulations"
           Lecture Notes in Mathematics Volume 1612, 1995
        """
        return simplicial_complexes.ComplexProjectivePlane().link([9])

    ###############################################################
    # examples from graph theory:

    def NotIConnectedGraphs(self, n, i):
        """
        The simplicial complex of all graphs on `n` vertices which are
        not `i`-connected.

        Fix an integer `n>0` and consider the set of graphs on `n`
        vertices.  View each graph as its set of edges, so it is a
        subset of a set of size `n` choose 2.  A graph is
        `i`-connected if, for any `j<i`, if any `j` vertices are
        removed along with the edges emanating from them, then the
        graph remains connected.  Now fix `i`: it is clear that if `G`
        is not `i`-connected, then the same is true for any graph
        obtained from `G` by deleting edges. Thus the set of all
        graphs which are not `i`-connected, viewed as a set of subsets
        of the `n` choose 2 possible edges, is closed under taking
        subsets, and thus forms a simplicial complex.  This function
        produces that simplicial complex.

        INPUT:

        -  ``n``, ``i`` -- non-negative integers with `i` at most `n`

        See Dumas et al. [DHSW2003]_ for information on computing its homology
        by computer, and see Babson et al. [BBLSW1999]_ for theory.  For
        example, Babson et al. show that when `i=2`, the reduced homology of
        this complex is nonzero only in dimension `2n-5`, where it is
        free abelian of rank `(n-2)!`.

        EXAMPLES::

            sage: simplicial_complexes.NotIConnectedGraphs(5,2).f_vector()
            [1, 10, 45, 120, 210, 240, 140, 20]
            sage: simplicial_complexes.NotIConnectedGraphs(5,2).homology(5).ngens()
            6

        REFERENCES:

        .. [BBLSW1999] Babson, Bjorner, Linusson, Shareshian, and Welker,
           "Complexes of not i-connected graphs," Topology 38 (1999),
           271-299

        .. [DHSW2003] Dumas, Heckenbach, Saunders, Welker, "Computing simplicial
           homology based on efficient Smith normal form algorithms,"
           in "Algebra, geometry, and software systems" (2003),
           177-206.
        """
        G_list = range(1,n+1)
        G_vertices = Set(G_list)
        E_list = []
        for w in G_list:
            for v in range(1,w):
                E_list.append((v,w))
        E = Set(E_list)
        facets = []
        i_minus_one_sets = list(G_vertices.subsets(size=i-1))
        for A in i_minus_one_sets:
            G_minus_A = G_vertices.difference(A)
            for B in G_minus_A.subsets():
                if len(B) > 0 and len(B) < len(G_minus_A):
                    C = G_minus_A.difference(B)
                    facet = E
                    for v in B:
                        for w in C:
                            bad_edge = (min(v,w), max(v,w))
                            facet = facet.difference(Set([bad_edge]))
                    facets.append(facet)
        return SimplicialComplex(facets, is_mutable=False)

    def MatchingComplex(self, n):
        """
        The matching complex of graphs on `n` vertices.

        Fix an integer `n>0` and consider a set `V` of `n` vertices.
        A 'partial matching' on `V` is a graph formed by edges so that
        each vertex is in at most one edge.  If `G` is a partial
        matching, then so is any graph obtained by deleting edges from
        `G`.  Thus the set of all partial matchings on `n` vertices,
        viewed as a set of subsets of the `n` choose 2 possible edges,
        is closed under taking subsets, and thus forms a simplicial
        complex called the 'matching complex'.  This function produces
        that simplicial complex.

        INPUT:

        -  ``n`` -- positive integer.

        See Dumas et al. [DHSW2003]_ for information on computing its homology
        by computer, and see Wachs [Wa2003]_ for an expository article about
        the theory.  For example, the homology of these complexes seems to
        have only mod 3 torsion, and this has been proved for the
        bottom non-vanishing homology group for the matching complex `M_n`.

        EXAMPLES::

            sage: M = simplicial_complexes.MatchingComplex(7)
            sage: H = M.homology()
            sage: H
            {0: 0, 1: C3, 2: Z^20}
            sage: H[2].ngens()
            20
            sage: simplicial_complexes.MatchingComplex(8).homology(2)  # long time (6s on sage.math, 2012)
            Z^132

        REFERENCES:

        .. [Wa2003] Wachs, "Topology of Matching, Chessboard and General Bounded
           Degree Graph Complexes" (Algebra Universalis Special Issue
           in Memory of Gian-Carlo Rota, Algebra Universalis, 49 (2003)
           345-385)
        """
        G_vertices = Set(range(1,n+1))
        facets = []
        if is_even(n):
            half = int(n/2)
            half_n_sets = list(G_vertices.subsets(size=half))
        else:
            half = int((n-1)/2)
            half_n_sets = list(G_vertices.subsets(size=half))
        for X in half_n_sets:
            Xcomp = G_vertices.difference(X)
            if is_even(n):
                if 1 in X:
                    A = X
                    B = Xcomp
                else:
                    A = Xcomp
                    B = X
                for M in matching(A, B):
                    facet = []
                    for pair in M:
                        facet.append(tuple(sorted(pair)))
                        facets.append(facet)
            else:
                for w in Xcomp:
                    if 1 in X or (w == 1 and 2 in X):
                        A = X
                        B = Xcomp.difference([w])
                    else:
                        B = X
                        A = Xcomp.difference([w])
                    for M in matching(A, B):
                        facet = []
                        for pair in M:
                            facet.append(tuple(sorted(pair)))
                        facets.append(facet)
        return SimplicialComplex(facets, is_mutable=False)

    def ChessboardComplex(self, n, i):
        r"""
        The chessboard complex for an `n \times i` chessboard.

        Fix integers `n, i > 0` and consider sets `V` of `n` vertices
        and `W` of `i` vertices.  A 'partial matching' between `V` and
        `W` is a graph formed by edges `(v,w)` with `v \in V` and `w
        \in W` so that each vertex is in at most one edge.  If `G` is
        a partial matching, then so is any graph obtained by deleting
        edges from `G`.  Thus the set of all partial matchings on `V`
        and `W`, viewed as a set of subsets of the `n+i` choose 2
        possible edges, is closed under taking subsets, and thus forms
        a simplicial complex called the 'chessboard complex'.  This
        function produces that simplicial complex.  (It is called the
        chessboard complex because such graphs also correspond to ways
        of placing rooks on an `n` by `i` chessboard so that none of
        them are attacking each other.)

        INPUT:

        -  ``n, i`` -- positive integers.

        See Dumas et al. [DHSW2003]_ for information on computing its homology
        by computer, and see Wachs [Wa2003]_ for an expository article about
        the theory.

        EXAMPLES::

            sage: C = simplicial_complexes.ChessboardComplex(5,5)
            sage: C.f_vector()
            [1, 25, 200, 600, 600, 120]
            sage: simplicial_complexes.ChessboardComplex(3,3).homology()
            {0: 0, 1: Z x Z x Z x Z, 2: 0}
        """
        A = range(n)
        B = range(i)
        E_dict = {}
        index = 0
        for v in A:
            for w in B:
                E_dict[(v,w)] = index
                index += 1
        facets = []
        for M in matching(A, B):
            facet = []
            for pair in M:
                facet.append(E_dict[pair])
            facets.append(facet)
        return SimplicialComplex(facets, is_mutable=False)

    def RandomComplex(self, n, d, p=0.5):
        """
        A random ``d``-dimensional simplicial complex on ``n`` vertices.

        INPUT:

        - ``n`` -- number of vertices

        - ``d`` -- dimension of the complex

        -  ``p`` -- floating point number between 0 and 1
           (optional, default 0.5)

        A random `d`-dimensional simplicial complex on `n` vertices,
        as defined for example by Meshulam and Wallach [MW2009]_, is
        constructed as follows: take `n` vertices and include all of
        the simplices of dimension strictly less than `d`, and then for each
        possible simplex of dimension `d`, include it with probability `p`.

        EXAMPLES::

            sage: X = simplicial_complexes.RandomComplex(6, 2); X
            Simplicial complex with vertex set (0, 1, 2, 3, 4, 5) and 10 facets
            sage: len(list(X.vertices()))
            6

        If `d` is too large (if `d+1 > n`, so that there are no
        `d`-dimensional simplices), then return the simplicial complex
        with a single `(n+1)`-dimensional simplex::

            sage: simplicial_complexes.RandomComplex(6, 12)
            Simplicial complex with vertex set (0, 1, 2, 3, 4, 5) and facets {(0, 1, 2, 3, 4, 5)}

        REFERENCES:

        .. [MW2009] Meshulam and Wallach, "Homological connectivity of random
           `k`-dimensional complexes", preprint, math.CO/0609773.
        """
        if d+1 > n:
            return simplicial_complexes.Simplex(n-1)
        else:
            vertices = range(n)
            facets = Subsets(vertices, d).list()
            maybe = Subsets(vertices, d+1)
            facets.extend([f for f in maybe if random.random() <= p])
            return SimplicialComplex(facets, is_mutable=False)

    def SumComplex(self, n, A):
        r"""
        The sum complexes of Linial, Meshulam, and Rosenthal [LMR2010]_.

        If `k+1` is the cardinality of `A`, then this returns a
        `k`-dimensional simplicial complex `X_A` with vertices
        `\ZZ/(n)`, and facets given by all `k+1`-tuples `(x_0, x_1,
        ..., x_k)` such that the sum `\sum x_i` is in `A`. See the
        paper by Linial, Meshulam, and Rosenthal [LMR2010]_, in which
        they prove various results about these complexes; for example,
        if `n` is prime, then `X_A` is rationally acyclic, and if in
        addition `A` forms an arithmetic progression in `\ZZ/(n)`,
        then `X_A` is `\ZZ`-acyclic. Throughout their paper, they
        assume that `n` and `k` are relatively prime, but the
        construction makes sense in general.

        In addition to the results from the cited paper, these
        complexes can have large torsion, given the number of
        vertices; for example, if `n=10`, and `A=\{0,1,2,3,6\}`, then
        `H_3(X_A)` is cyclic of order 2728, and there is a
        4-dimensional complex on 13 vertices with `H_3` having a
        cyclic summand of order

        .. MATH::

            706565607945 = 3 \cdot 5 \cdot 53 \cdot 79 \cdot 131
            \cdot 157 \cdot 547.

        See the examples.

        INPUT:

        - ``n`` -- a positive integer

        - ``A`` -- a subset of `\ZZ/(n)`

        REFERENCES:

        .. [LMR2010] N. Linial, R. Meshulam and M. Rosenthal, "Sum
           complexes -- a new family of hypertrees", Discrete &
           Computational Geometry, 2010, Volume 44, Number 3, Pages
           622-636

        EXAMPLES::

            sage: S = simplicial_complexes.SumComplex(10, [0,1,2,3,6]); S
            Simplicial complex with 10 vertices and 126 facets
            sage: S.homology()
            {0: 0, 1: 0, 2: 0, 3: C2728, 4: 0}
            sage: factor(2728)
            2^3 * 11 * 31

            sage: S = simplicial_complexes.SumComplex(11, [0, 1, 3]); S
            Simplicial complex with 11 vertices and 45 facets
            sage: S.homology(1)
            C23
            sage: S = simplicial_complexes.SumComplex(11, [0,1,2,3,4,7]); S
            Simplicial complex with 11 vertices and 252 facets
            sage: S.homology() # long time
            {0: 0, 1: 0, 2: 0, 3: 0, 4: C645679, 5: 0}
            sage: factor(645679)
            23 * 67 * 419

            sage: S = simplicial_complexes.SumComplex(13, [0, 1, 3]); S
            Simplicial complex with 13 vertices and 66 facets
            sage: S.homology(1)
            C159
            sage: factor(159)
            3 * 53
            sage: S = simplicial_complexes.SumComplex(13, [0,1,2,5]); S
            Simplicial complex with 13 vertices and 220 facets
            sage: S.homology() # long time
            {0: 0, 1: 0, 2: C146989209, 3: 0}
            sage: factor(1648910295)
            3^2 * 5 * 53 * 521 * 1327
            sage: S = simplicial_complexes.SumComplex(13, [0,1,2,3,5]); S
            Simplicial complex with 13 vertices and 495 facets
            sage: S.homology() # long time
            {0: 0, 1: 0, 2: 0, 3: C3 x C237 x C706565607945, 4: 0}
            sage: factor(706565607945)
            3 * 5 * 53 * 79 * 131 * 157 * 547

            sage: S = simplicial_complexes.SumComplex(17, [0, 1, 4]); S
            Simplicial complex with 17 vertices and 120 facets
            sage: S.homology(1)
            C140183
            sage: factor(140183)
            103 * 1361
            sage: S = simplicial_complexes.SumComplex(19, [0, 1, 4]); S
            Simplicial complex with 19 vertices and 153 facets
            sage: S.homology(1)
            C5670599
            sage: factor(5670599)
            11 * 191 * 2699
            sage: S = simplicial_complexes.SumComplex(31, [0, 1, 4]); S
            Simplicial complex with 31 vertices and 435 facets
            sage: S.homology(1) # long time
            C5 x C5 x C5 x C5 x C26951480558170926865
            sage: factor(26951480558170926865)
            5 * 311 * 683 * 1117 * 11657 * 1948909
        """
        from sage.rings.all import Integers
        Zn = Integers(n)
        A = frozenset([Zn(x) for x in A])
        facets = []
        for f in Set(Zn).subsets(len(A)):
            if sum(f) in A:
                facets.append(tuple(f))
        return SimplicialComplex(facets, is_mutable=False)

simplicial_complexes = SimplicialComplexExamples()
