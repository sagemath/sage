# -*- coding: utf-8 -*-
"""
Examples of simplicial complexes

There are two main types: manifolds and examples related to graph
theory.

For manifolds, there are functions defining the `n`-sphere for any
`n`, the torus, `n`-dimensional real projective space for any `n`, the
complex projective plane, surfaces of arbitrary genus, and some other
manifolds, all as simplicial complexes.

Aside from surfaces, this file also provides functions for
constructing some other simplicial complexes: the simplicial complex
of not-`i`-connected graphs on `n` vertices, the matching complex on n
vertices, the chessboard complex for an `n` by `i` chessboard, and
others.  These provide examples of large simplicial complexes; for
example, ``simplicial_complexes.NotIConnectedGraphs(7,2)`` has over a
million simplices.

All of these examples are accessible by typing
``simplicial_complexes.NAME``, where ``NAME`` is the name of the example.

- :func:`BarnetteSphere`
- :func:`BrucknerGrunbaumSphere`
- :func:`ChessboardComplex`
- :func:`ComplexProjectivePlane`
- :func:`DunceHat`
- :func:`K3Surface`
- :func:`KleinBottle`
- :func:`MatchingComplex`
- :func:`MooreSpace`
- :func:`NotIConnectedGraphs`
- :func:`PoincareHomologyThreeSphere`
- :func:`PseudoQuaternionicProjectivePlane`
- :func:`RandomComplex`
- :func:`RandomTwoSphere`
- :func:`RealProjectivePlane`
- :func:`RealProjectiveSpace`
- :func:`RudinBall`
- :func:`ShiftedComplex`
- :func:`Simplex`
- :func:`Sphere`
- :func:`SumComplex`
- :func:`SurfaceOfGenus`
- :func:`Torus`
- :func:`ZieglerBall`

You can also get a list by typing ``simplicial_complexes.`` and hitting the
TAB key.

EXAMPLES::

    sage: S = simplicial_complexes.Sphere(2) # the 2-sphere
    sage: S.homology()
    {0: 0, 1: 0, 2: Z}
    sage: simplicial_complexes.SurfaceOfGenus(3)
    Triangulation of an orientable surface of genus 3
    sage: M4 = simplicial_complexes.MooreSpace(4)
    sage: M4.homology()
    {0: 0, 1: C4, 2: 0}
    sage: simplicial_complexes.MatchingComplex(6).homology()
    {0: 0, 1: Z^16, 2: 0}
"""

from sage.homology.simplicial_complex import SimplicialComplex
from sage.structure.unique_representation import UniqueRepresentation
# Below we define a function Simplex to construct a simplex as a
# simplicial complex. We also need to use actual simplices as
# simplices, hence:
from sage.homology.simplicial_complex import Simplex as TrueSimplex
from sage.sets.set import Set
from sage.misc.functional import is_even
from sage.misc.misc import union
from sage.combinat.subset import Subsets
import sage.misc.prandom as random

# Miscellaneous utility functions.

# The following two functions can be used to generate the facets for
# the corresponding examples in sage.homology.examples. These take a
# few seconds to run, so the actual examples have the facets
# hard-coded. Thus the following functions are not currently used in
# the Sage library.

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

class UniqueSimplicialComplex(SimplicialComplex, UniqueRepresentation):
    """
    This combines :class:`SimplicialComplex` and
    :class:`UniqueRepresentation`. It is intended to be used to make
    standard examples of simplicial complexes unique. See :trac:`13566`.

    INPUTS:

    - the inputs are the same as for a :class:`SimplicialComplex`,
      with one addition and two exceptions. The exceptions are that
      ``is_mutable`` and ``is_immutable`` are ignored: all instances
      of this class are immutable. The addition:

    - ``name`` -- string (optional), the string representation for this complex.

    EXAMPLES::

        sage: from sage.homology.examples import UniqueSimplicialComplex
        sage: SimplicialComplex([[0,1]]) is SimplicialComplex([[0,1]])
        False
        sage: UniqueSimplicialComplex([[0,1]]) is UniqueSimplicialComplex([[0,1]])
        True
        sage: UniqueSimplicialComplex([[0,1]])
        Simplicial complex with vertex set (0, 1) and facets {(0, 1)}
        sage: UniqueSimplicialComplex([[0,1]], name='The 1-simplex')
        The 1-simplex
    """
    @staticmethod
    def __classcall__(self, maximal_faces=None, name=None, **kwds):
        """
        TESTS::

            sage: from sage.homology.examples import UniqueSimplicialComplex
            sage: UniqueSimplicialComplex([[1,2,3], [0,1,3]]) is UniqueSimplicialComplex([(1,2,3), (0,1,3)])
            True
            sage: X = UniqueSimplicialComplex([[1,2,3], [0,1,3]])
            sage: X is UniqueSimplicialComplex(X)
            True

        Testing ``from_characteristic_function``::

            sage: UniqueSimplicialComplex(from_characteristic_function=(lambda x:sum(x)<=4, range(5)))
            Simplicial complex with vertex set (0, 1, 2, 3, 4) and facets {(0, 4), (0, 1, 2), (0, 1, 3)}
        """
        char_fcn = kwds.get('from_characteristic_function', None)
        if char_fcn:
            kwds['from_characteristic_function'] = (char_fcn[0], tuple(char_fcn[1]))
        if maximal_faces:
            # Test to see if maximal_faces is a cell complex or another
            # object which can be converted to a simplicial complex:
            C = None
            if isinstance(maximal_faces, SimplicialComplex):
                C = maximal_faces
            else:
                try:
                    C = maximal_faces._simplicial_()
                except AttributeError:
                    if not isinstance(maximal_faces, (list, tuple, Simplex)):
                        # Convert it into a list (in case it is an iterable)
                        maximal_faces = list(maximal_faces)
                    if len(maximal_faces) != 0:
                        vertex_set = reduce(union, maximal_faces)
            if C is not None:
                maximal_faces = C.facets()
            # Now convert maximal_faces to a tuple of tuples, so that it is hashable.
            maximal_faces = tuple([tuple(_) for _ in maximal_faces])
        return super(UniqueSimplicialComplex, self).__classcall__(self, maximal_faces,
                                                                  name=name,
                                                                  **kwds)

    def __init__(self, maximal_faces=None, name=None, **kwds):
        """
        TESTS::

            sage: from sage.homology.examples import UniqueSimplicialComplex
            sage: UniqueSimplicialComplex([[1,2,3], [0,1,3]], is_mutable=True).is_mutable()
            False
        """
        if 'is_mutable' in kwds:
            del kwds['is_mutable']
        if 'is_immutable' in kwds:
            del kwds['is_immutable']
        self._name = name
        SimplicialComplex.__init__(self, maximal_faces=maximal_faces, is_mutable=False, **kwds)

    def _repr_(self):
        """
        Print representation

        If the argument ``name`` was specified when defining the
        complex, use that. Otherwise, use the print representation
        from the class :class:`SimplicialComplex`.

        TESTS::

            sage: from sage.homology.examples import UniqueSimplicialComplex
            sage: UniqueSimplicialComplex([[0,1]])
            Simplicial complex with vertex set (0, 1) and facets {(0, 1)}
            sage: UniqueSimplicialComplex([[0,1]], name='Joe')
            Joe
        """
        if self._name:
            return self._name
        return SimplicialComplex._repr_(self)

# Now the functions that produce the actual examples...

def Sphere(n):
    """
    A minimal triangulation of the `n`-dimensional sphere.

    INPUT:

    -  ``n`` -- positive integer

    EXAMPLES::

        sage: simplicial_complexes.Sphere(2)
        Minimal triangulation of the 2-sphere
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
    S = TrueSimplex(n+1)
    facets = tuple(S.faces())
    return UniqueSimplicialComplex(facets,
                                   name='Minimal triangulation of the {}-sphere'.format(n))

def Simplex(n):
    """
    An `n`-dimensional simplex, as a simplicial complex.

    INPUT:

    -  ``n`` -- a non-negative integer

    OUTPUT: the simplicial complex consisting of the `n`-simplex
    on vertices `(0, 1, ..., n)` and all of its faces.

    EXAMPLES::

        sage: simplicial_complexes.Simplex(3)
        The 3-simplex
        sage: simplicial_complexes.Simplex(5).euler_characteristic()
        1
    """
    return UniqueSimplicialComplex([TrueSimplex(n)],
                                   name='The {}-simplex'.format(n))

def Torus():
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
    return UniqueSimplicialComplex([[0,1,2], [1,2,4], [1,3,4], [1,3,6],
                              [0,1,5], [1,5,6], [2,3,5], [2,4,5],
                              [2,3,6], [0,2,6], [0,3,4], [0,3,5],
                              [4,5,6], [0,4,6]],
                                   name='Minimal triangulation of the torus')

def RealProjectivePlane():
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
    return UniqueSimplicialComplex([[0,1,2], [0,2,3], [0,1,5], [0,4,5],
                              [0,3,4], [1,2,4], [1,3,4], [1,3,5],
                              [2,3,5], [2,4,5]],
                              name='Minimal triangulation of the real projective plane')

ProjectivePlane = RealProjectivePlane

def KleinBottle():
    """
    A minimal triangulation of the Klein bottle, as presented for example
    in Davide Cervone's thesis [Ce1994]_.

    EXAMPLES::

        sage: simplicial_complexes.KleinBottle()
        Minimal triangulation of the Klein bottle

    REFERENCES:

    .. [Ce1994] \D. P. Cervone, "Vertex-minimal simplicial immersions of the Klein
       bottle in three-space", Geom. Ded. 50 (1994) 117-141,
       http://www.math.union.edu/~dpvc/papers/1993-03.kb/vmkb.pdf.
    """
    return UniqueSimplicialComplex([[2,3,7], [1,2,3], [1,3,5], [1,5,7],
                              [1,4,7], [2,4,6], [1,2,6], [1,6,0],
                              [1,4,0], [2,4,0], [3,4,7], [3,4,6],
                              [3,5,6], [5,6,0], [2,5,0], [2,5,7]],
                            name='Minimal triangulation of the Klein bottle')

def SurfaceOfGenus(g, orientable=True):
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
        Triangulation of an orientable surface of genus 2
        sage: simplicial_complexes.SurfaceOfGenus(1, orientable=False)
        Triangulation of a non-orientable surface of genus 1
    """
    if g == 0:
        if not orientable:
            raise ValueError("No non-orientable surface of genus zero.")
        else:
            return Sphere(2)
    if orientable:
        T = Torus()
    else:
        T = RealProjectivePlane()
    S = T
    for i in range(g-1):
        S = S.connected_sum(T)
    if orientable:
        orient_str = 'n orientable'
    else:
        orient_str = ' non-orientable'
    return UniqueSimplicialComplex(S,
                name='Triangulation of a{} surface of genus {}'.format(orient_str, g))

def MooreSpace(q):
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
        Minimal triangulation of the real projective plane
        sage: simplicial_complexes.MooreSpace(3).homology()[1]
        C3
        sage: simplicial_complexes.MooreSpace(4).suspension().homology()[2]
        C4
        sage: simplicial_complexes.MooreSpace(8)
        Triangulation of the mod 8 Moore space
    """
    if q <= 1:
        raise ValueError("The mod q Moore space is only defined if q is at least 2")
    if q == 2:
        return RealProjectivePlane()
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
    return UniqueSimplicialComplex(facets,
                                   name='Triangulation of the mod {} Moore space'.format(q))

def ComplexProjectivePlane():
    """
    A minimal triangulation of the complex projective plane.

    This was constructed by Kühnel and Banchoff [KB1983]_.

    REFERENCES:

    .. [KB1983] \W. Kühnel and T. F. Banchoff, "The 9-vertex complex
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
    return UniqueSimplicialComplex(
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
        name='Minimal triangulation of the complex projective plane')

def PseudoQuaternionicProjectivePlane():
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
    return UniqueSimplicialComplex([ [g(index) for index in tuple]
            for tuple in start_list
            for g in PermutationGroup([P,S]) ])

def PoincareHomologyThreeSphere():
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
    return UniqueSimplicialComplex(
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
        name='Triangulation of the Poincare homology 3-sphere')

def RealProjectiveSpace(n):
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
    can be computed by the function :func:`~sage.homology.simplicial_complex.facets_for_RP4`, but
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

    .. [Ku1987] \W. Kühnel, "Minimal triangulations of Kummer varieties",
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
        return Simplex(0)
    if n == 1:
        return Sphere(1)
    if n == 2:
        return RealProjectivePlane()
    if n == 3:
        # Minimal triangulation found by Walkup and given
        # explicitly by Lutz
        return UniqueSimplicialComplex(
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
              name='Minimal triangulation of RP^3')
    if n == 4:
        return UniqueSimplicialComplex(
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
            name='Minimal triangulation of RP^4')
    if n >= 5:
        # Use the construction given by Datta in Example 3.21.
        V = set(range(0, n+2))
        S = Sphere(n).barycentric_subdivision()
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
        return UniqueSimplicialComplex(list(facets),
                                       name='Triangulation of RP^{}'.format(n))

def K3Surface():
    """
    Returns a minimal triangulation of the K3 surface.

    This is a pure simplicial complex of dimension 4 with 16 vertices
    and 288 facets. It was constructed by Casella and Kühnel
    in [CK2001]_. The construction here uses the labeling from
    Spreer and Kühnel [SK2011]_.

    REFERENCES:

    .. [CK2001] \M. Casella and W. Kühnel, "A triangulated K3 surface
       with the minimum number of vertices", Topology 40 (2001),
       753--772.

    .. [SK2011] \J. Spreer and W. Kühnel, "Combinatorial properties
       of the K3 surface: Simplicial blowups and slicings", Experimental
       Mathematics, Volume 20, Issue 2, 2011.

    EXAMPLES::

        sage: K3=simplicial_complexes.K3Surface() ; K3
        Minimal triangulation of the K3 surface
        sage: K3.f_vector()
        [1, 16, 120, 560, 720, 288]

    This simplicial complex is implemented just by listing all 288
    facets. The list of facets can be computed by the function
    :func:`~sage.homology.simplicial_complex.facets_for_K3`, but running the function takes a few
    seconds.
    """
    return UniqueSimplicialComplex(
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
        name='Minimal triangulation of the K3 surface')

def BarnetteSphere():
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
        Barnette's triangulation of the 3-sphere
        sage: BS.f_vector()
        [1, 8, 27, 38, 19]

    TESTS:

    Checks that this is indeed the same Barnette Sphere as the one
    given on page 87 of [E1996]_.::

        sage: BS2 = SimplicialComplex([[1,2,3,4],[3,4,5,6],[1,2,5,6],
        ....:                          [1,2,4,7],[1,3,4,7],[3,4,6,7],
        ....:                          [3,5,6,7],[1,2,5,7],[2,5,6,7],
        ....:                          [2,4,6,7],[1,2,3,8],[2,3,4,8],
        ....:                          [3,4,5,8],[4,5,6,8],[1,2,6,8],
        ....:                          [1,5,6,8],[1,3,5,8],[2,4,6,8],
        ....:                          [1,3,5,7]])
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
    return UniqueSimplicialComplex([
            (1,2,4,5),(2,3,5,6),(1,3,4,6),(1,2,3,7),(4,5,6,7),(1,2,4,7),
            (2,4,5,7),(2,3,5,7),(3,5,6,7),(3,1,6,7),(1,6,4,7),(1,2,3,8),
            (4,5,6,8),(1,2,5,8),(1,4,5,8),(2,3,6,8),(2,5,6,8),(3,1,4,8),
            (3,6,4,8)],
          name="Barnette's triangulation of the 3-sphere")

def BrucknerGrunbaumSphere():
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
        Bruckner and Grunbaum's triangulation of the 3-sphere
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
    # X = ComplexProjectivePlane().link([9])
    # return UniqueSimplicialComplex(X.facets(),
    #                                name="Bruckner and Grunbaum's triangulation of the 3-sphere")
    return UniqueSimplicialComplex(ComplexProjectivePlane().link([9]),
                                   name="Bruckner and Grunbaum's triangulation of the 3-sphere")

###############################################################
# examples from graph theory:

def NotIConnectedGraphs(n, i):
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
    return UniqueSimplicialComplex(facets, name='Simplicial complex of not {}-connected graphs on {} vertices'.format(i, n))

def MatchingComplex(n):
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
    return UniqueSimplicialComplex(facets, name='Matching complex on {} vertices'.format(n))

def ChessboardComplex(n, i):
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
    return UniqueSimplicialComplex(facets, name='Chessboard complex for an {}x{} chessboard'.format(n, i))

def RandomComplex(n, d, p=0.5):
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
        Random 2-dimensional simplicial complex on 6 vertices
        sage: len(list(X.vertices()))
        6

    If `d` is too large (if `d+1 > n`, so that there are no
    `d`-dimensional simplices), then return the simplicial complex
    with a single `(n+1)`-dimensional simplex::

        sage: simplicial_complexes.RandomComplex(6, 12)
        The 5-simplex

    REFERENCES:

    .. [MW2009] Meshulam and Wallach, "Homological connectivity of random
       `k`-dimensional complexes", preprint, math.CO/0609773.
    """
    if d+1 > n:
        return Simplex(n-1)
    else:
        vertices = range(n)
        facets = Subsets(vertices, d).list()
        maybe = Subsets(vertices, d+1)
        facets.extend([f for f in maybe if random.random() <= p])
        return UniqueSimplicialComplex(facets,
                                       name='Random {}-dimensional simplicial complex on {} vertices'.format(d,n))

def SumComplex(n, A):
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

    .. [LMR2010] \N. Linial, R. Meshulam and M. Rosenthal, "Sum
       complexes -- a new family of hypertrees", Discrete &
       Computational Geometry, 2010, Volume 44, Number 3, Pages
       622-636

    EXAMPLES::

        sage: S = simplicial_complexes.SumComplex(10, [0,1,2,3,6]); S
        Sum complex on vertices Z/10Z associated to {0, 1, 2, 3, 6}
        sage: S.homology()
        {0: 0, 1: 0, 2: 0, 3: C2728, 4: 0}
        sage: factor(2728)
        2^3 * 11 * 31

        sage: S = simplicial_complexes.SumComplex(11, [0, 1, 3]); S
        Sum complex on vertices Z/11Z associated to {0, 1, 3}
        sage: S.homology(1)
        C23
        sage: S = simplicial_complexes.SumComplex(11, [0,1,2,3,4,7]); S
        Sum complex on vertices Z/11Z associated to {0, 1, 2, 3, 4, 7}
        sage: S.homology(algorithm='no_chomp') # long time
        {0: 0, 1: 0, 2: 0, 3: 0, 4: C645679, 5: 0}
        sage: factor(645679)
        23 * 67 * 419

        sage: S = simplicial_complexes.SumComplex(13, [0, 1, 3]); S
        Sum complex on vertices Z/13Z associated to {0, 1, 3}
        sage: S.homology(1)
        C159
        sage: factor(159)
        3 * 53
        sage: S = simplicial_complexes.SumComplex(13, [0,1,2,5]); S
        Sum complex on vertices Z/13Z associated to {0, 1, 2, 5}
        sage: S.homology(algorithm='no_chomp') # long time
        {0: 0, 1: 0, 2: C146989209, 3: 0}
        sage: factor(1648910295)
        3^2 * 5 * 53 * 521 * 1327
        sage: S = simplicial_complexes.SumComplex(13, [0,1,2,3,5]); S
        Sum complex on vertices Z/13Z associated to {0, 1, 2, 3, 5}
        sage: S.homology(algorithm='no_chomp') # long time
        {0: 0, 1: 0, 2: 0, 3: C3 x C237 x C706565607945, 4: 0}
        sage: factor(706565607945)
        3 * 5 * 53 * 79 * 131 * 157 * 547

        sage: S = simplicial_complexes.SumComplex(17, [0, 1, 4]); S
        Sum complex on vertices Z/17Z associated to {0, 1, 4}
        sage: S.homology(1, algorithm='no_chomp')
        C140183
        sage: factor(140183)
        103 * 1361
        sage: S = simplicial_complexes.SumComplex(19, [0, 1, 4]); S
        Sum complex on vertices Z/19Z associated to {0, 1, 4}
        sage: S.homology(1,algorithm='no_chomp')
        C5670599
        sage: factor(5670599)
        11 * 191 * 2699
        sage: S = simplicial_complexes.SumComplex(31, [0, 1, 4]); S
        Sum complex on vertices Z/31Z associated to {0, 1, 4}
        sage: S.homology(1,algorithm='no_chomp') # long time
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
    return UniqueSimplicialComplex(facets, name='Sum complex on vertices Z/{}Z associated to {}'.format(n, Set(A)))


def RandomTwoSphere(n):
    r"""
    Return a random triangulation of the 2-dimensional sphere with `n`
    vertices.

    INPUT:

    `n` -- an integer

    OUTPUT:

    A random triangulation of the sphere chosen uniformly among
    the *rooted* triangulations on `n` vertices. Because some
    triangulations have nontrivial automorphism groups, this may
    not be equal to the uniform distribution among unrooted
    triangulations.

    ALGORITHM:

    The algorithm is taken from [PS2006]_, section 2.1.

    Starting from a planar tree (represented by its contour as a
    sequence of vertices), one first performs local closures, until no
    one is possible. A local closure amounts to replace in the cyclic
    contour word a sequence ``in1,in2,in3,lf,in3`` by
    ``in1,in3``. After all local closures are done, one has reached
    the partial closure, as in [PS2006]_, figure 5 (a).

    Then one has to perform complete closure by adding two more
    vertices, in order to reach the situation of [PS2006]_, figure 5
    (b). For this, it is necessary to find inside the final contour
    one of the two subsequences ``lf,in,lf``.

    At every step of the algorithm, newly created triangles are added
    in a simplicial complex.

    This algorithm is implemented in
    :meth:`~sage.graphs.generators.random.RandomTriangulation`, which
    creates an embedded graph. The triangles of the simplicial
    complex are recovered from this embedded graph.

    EXAMPLES::

        sage: G = simplicial_complexes.RandomTwoSphere(6); G
        Simplicial complex with vertex set (0, 1, 2, 3, 'a', 'b')
        and 8 facets
        sage: G.homology()
        {0: 0, 1: 0, 2: Z}
        sage: G.is_pure()
        True
        sage: fg = G.flip_graph(); fg
        Graph on 8 vertices
        sage: fg.is_planar() and fg.is_regular(3)
        True
    """
    from sage.graphs.generators.random import RandomTriangulation

    graph = RandomTriangulation(n)

    graph = graph.relabel(inplace=False)
    triangles = [(u, v, w) for u, L in graph._embedding.iteritems()
                 for v, w in zip(L, L[1:] + [L[0]]) if u < v and u < w]

    return SimplicialComplex(triangles, maximality_check=False)

def ShiftedComplex(generators):
    r"""
    Return the smallest shifted simplicial complex containing ``generators``
    as faces.

    Let `V` be a set of vertices equipped with a total order.  The
    'componentwise partial ordering' on k-subsets of `V` is defined as
    follows: if `A = \{a_1 < \cdots < a_k\}` and `B = \{b_1 < \cdots < b_k\}`,
    then `A \leq_C B` iff `a_i \leq b_i` for all `i`.  A simplicial complex
    `X` on vertex set `[n]` is *shifted* if its faces form an order ideal
    under the componentwise partial ordering, i.e., if `B \in X` and
    `A \leq_C B` then `A \in X`.  Shifted complexes of dimension 1 are also
    known as threshold graphs.

    .. NOTE::

        This method assumes that `V` consists of positive integers
        with the natural ordering.

    INPUT:

    - ``generators`` -- a list of generators of the order ideal, which may
      be lists, tuples or simplices

    EXAMPLES::

        sage: X = simplicial_complexes.ShiftedComplex([ Simplex([1,6]), (2,4), [8] ])
        sage: X.facets()
        {(2, 4), (7,), (1, 2), (1, 5), (1, 4), (8,), (2, 3), (1, 6), (1, 3)}
        sage: X = simplicial_complexes.ShiftedComplex([ [2,3,5] ])
        sage: X.facets()
        {(1, 3, 4), (1, 3, 5), (2, 3, 5), (1, 2, 3), (2, 3, 4), (1, 2, 5), (1, 2, 4)}
        sage: X = simplicial_complexes.ShiftedComplex([ [1,3,5], [2,6] ])
        sage: X.facets()
        {(1, 3, 4), (1, 3, 5), (1, 6), (2, 6), (1, 2, 3), (1, 2, 5), (1, 2, 4)}
    """
    from sage.combinat.partition import Partitions
    Facets = []
    for G in generators:
        G = list(reversed(sorted(G)))
        L = len(G)
        for k in range(L * (L+1) // 2, sum(G) + 1):
            for P in Partitions(k, length=L, max_slope=-1, outer=G):
                Facets.append( list(reversed(P)) )
    return SimplicialComplex(Facets)

def RudinBall():
    r"""
    Return the non-shellable ball constructed by Rudin.

    This complex is a non-shellable triangulation of the 3-ball
    with 14 vertices and 41 facets, constructed by Rudin in
    [Ru1958]_.

    EXAMPLES::

        sage: R = simplicial_complexes.RudinBall(); R
        Rudin ball
        sage: R.f_vector()
        [1, 14, 66, 94, 41]
        sage: R.homology()
        {0: 0, 1: 0, 2: 0, 3: 0}
        sage: R.is_cohen_macaulay()
        True

    REFERENCES:

    .. [Ru1958] \M. E. Rudin.
       *An unshellable triangulation of a tetrahedron*.
       Bull. Amer. Math. Soc. 64 (1958), 90-91.
    """
    return UniqueSimplicialComplex(
        [[1,9,2,5], [1,10,2,5], [1,10,5,11], [1,10,7,11], [1,13,5,11],
         [1,13,7,11], [2,10,3,6], [2,11,3,6], [2,11,6,12], [2,11,8,12],
         [2,14,6,12], [2,14,8,12], [3,11,4,7], [3,12,4,7], [3,12,5,9],
         [3,12,7,9], [3,13,5,9], [3,13,7,9], [4,9,1,8], [4,9,6,10],
         [4,9,8,10], [4,12,1,8], [4,14,6,10], [4,14,8,10], [9,10,2,5],
         [9,10,2,6], [9,10,5,11], [9,10,11,12], [9,13,5,11], [10,11,3,6],
         [10,11,3,7], [10,11,6,12], [10,14,6,12], [11,12,4,7], [11,12,4,8],
         [11,12,7,9], [11,13,7,9], [12,9,1,5], [12,9,1,8], [12,9,8,10],
         [12,14,8,10]],
        name="Rudin ball"
    )

def ZieglerBall():
    r"""
    Return the non-shellable ball constructed by Ziegler.

    This complex is a non-shellable triangulation of the 3-ball
    with 10 vertices and 21 facets, constructed by Ziegler in
    [Zi1998]_ and the smallest such complex known.

    EXAMPLES::

        sage: Z = simplicial_complexes.ZieglerBall(); Z
        Ziegler ball
        sage: Z.f_vector()
        [1, 10, 38, 50, 21]
        sage: Z.homology()
        {0: 0, 1: 0, 2: 0, 3: 0}
        sage: Z.is_cohen_macaulay()
        True

    REFERENCES:

    .. [Zi1998] G. Ziegler.
       *Shelling polyhedral 3-balls and 4-polytopes*.
       Discrete Comput. Geom. 19 (1998), 159-174.
    """

    return UniqueSimplicialComplex(
       [[1,2,3,4], [1,2,5,6], [1,5,6,9], [2,5,6,0], [3,6,7,8], [4,5,7,8],
        [2,3,6,7], [1,6,2,9], [2,6,7,0], [3,2,4,8], [4,1,3,7], [3,4,7,8],
        [1,2,4,9], [2,7,3,0], [3,2,6,8], [4,1,5,7], [4,1,8,5], [1,4,8,9],
        [2,3,1,0], [1,8,5,9], [2,1,5,0]],
       name="Ziegler ball"
    )

def DunceHat():
    r"""
    Return the minimal triangulation of the dunce hat given by Hachimori
    [Ha2016]_.

    This is a standard example of a space that is contractible
    but not collapsible.

    EXAMPLES::

        sage: D = simplicial_complexes.DunceHat(); D
        Minimal triangulation of the dunce hat
        sage: D.f_vector()
        [1, 8, 24, 17]
        sage: D.homology()
        {0: 0, 1: 0, 2: 0}
        sage: D.is_cohen_macaulay()
        True

    REFERENCES:

    .. [Ha2016] M. Hachimori.
       http://infoshako.sk.tsukuba.ac.jp/~hachi/math/library/dunce_hat_eng.html
    """

    return UniqueSimplicialComplex(
       [[1,3,5], [2,3,5], [2,4,5], [1,2,4], [1,3,4], [3,4,8],
        [1,2,8], [1,7,8], [1,2,7], [2,3,7], [3,6,7], [1,3,6],
        [1,5,6], [4,5,6], [4,6,8], [6,7,8], [2,3,8]],
       name="Minimal triangulation of the dunce hat"
    )


# For taking care of old pickles
from sage.structure.sage_object import register_unpickle_override
register_unpickle_override('sage.homology.examples', 'SimplicialSurface', SimplicialComplex)
