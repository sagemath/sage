r"""
Group homomorphisms for groups with a GAP backend

EXAMPLES::

    sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
    sage: A = AbelianGroupGap([2, 4])
    sage: F.<a,b> = FreeGroup()
    sage: f = F.hom([g for g in A.gens()])
    sage: K = f.kernel()
    sage: K
    Group(<free, no generators known>)

AUTHORS:

- Simon Brandhorst (2018-02-08): initial version
- Sebastian Oehms  (2018-11-15): have this functionality work for permutation groups (:trac:`26750`)
  and implement :meth:`section` and :meth:`natural_map`
"""

# ****************************************************************************
#       Copyright (C) 2018 Simon Brandhorst <sbrandhorst@web.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.homset import HomsetWithBase
from sage.categories.morphism import Morphism
from sage.rings.integer_ring import ZZ
from sage.groups.libgap_wrapper import ParentLibGAP
from sage.libs.gap.element import GapElement
from sage.misc.latex import latex


class GroupMorphism_libgap(Morphism):
    r"""
    This wraps GAP group homomorphisms.

    Checking if the input defines a group homomorphism can be expensive
    if the group is large.

    INPUT:

    - ``homset`` -- the parent
    - ``gap_hom`` -- a :class:`sage.libs.gap.element.GapElement` consisting of a group homomorphism
    - ``check`` -- (default: ``True``) check if the ``gap_hom`` is a group
      homomorphism; this can be expensive

    EXAMPLES::

        sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
        sage: A = AbelianGroupGap([2, 4])
        sage: A.hom([g^2 for g in A.gens()])
        Group endomorphism of Abelian group with gap, generator orders (2, 4)

    Homomorphisms can be defined between different kinds of GAP groups::

        sage: G = MatrixGroup([Matrix(ZZ, 2, [0,1,1,0])])
        sage: f = A.hom([G.0, G(1)])
        sage: f
        Group morphism:
        From: Abelian group with gap, generator orders (2, 4)
        To:   Matrix group over Integer Ring with 1 generators (
        [0 1]
        [1 0]
        )
        sage: G.<a,b> = FreeGroup()
        sage: H = G / (G([1]), G([2])^3)
        sage: f = G.hom(H.gens())
        sage: f
        Group morphism:
          From: Free Group on generators {a, b}
          To:   Finitely presented group < a, b | a, b^3 >

    Homomorphisms can be defined between GAP groups and permutation groups::

        sage: S = Sp(4,3)
        sage: P = PSp(4,3)
        sage: pr = S.hom(P.gens())
        sage: E = copy(S.one().matrix())
        sage: E[3,0] = 2; e = S(E)
        sage: pr(e)
        (1,16,15)(3,22,18)(4,19,21)(6,34,24)(7,25,33)(9,40,27)(10,28,39)(12,37,30)(13,31,36)

    TESTS:

    Old tests inherited from MatrixGroupMorphisms::

        sage: F = GF(5); MS = MatrixSpace(F,2,2)
        sage: G = MatrixGroup([MS([1,1,0,1])])
        sage: H = MatrixGroup([MS([1,0,1,1])])
        sage: phi = G.hom(H.gens())
        sage: phi
        Group morphism:
        From: Matrix group over Finite Field of size 5 with 1 generators (
        [1 1]
        [0 1]
        )
        To:   Matrix group over Finite Field of size 5 with 1 generators (
        [1 0]
        [1 1]
        )
        sage: phi(MS([1,1,0,1]))
        [1 0]
        [1 1]
        sage: F = GF(7); MS = MatrixSpace(F,2,2)
        sage: F.multiplicative_generator()
        3
        sage: G = MatrixGroup([MS([3,0,0,1])])
        sage: a = G.gens()[0]^2
        sage: phi = G.hom([a])

    Check that :trac:`19406` is fixed::

        sage: G = GL(2, GF(3))
        sage: H = GL(3, GF(2))
        sage: mat1 = H([[-1,0,0],[0,0,-1],[0,-1,0]])
        sage: mat2 = H([[1,1,1],[0,0,-1],[-1,0,0]])
        sage: phi = G.hom([mat1, mat2])
        Traceback (most recent call last):
        ...
        ValueError: images do not define a group homomorphism

        sage: F = GF(5); MS = MatrixSpace(F,2,2)
        sage: G = MatrixGroup([MS([1,1,0,1])])
        sage: H = MatrixGroup([MS([1,0,1,1])])
        sage: phi = G.hom(H.gens())
        sage: phi.gap()
        CompositionMapping( [ (6,7,8,10,9)(11,13,14,12,15)(16,19,20,18,17)(21,25,22,24,23) ]
        -> [ [ [ Z(5)^0, 0*Z(5) ], [ Z(5)^0, Z(5)^0 ] ] ], <action isomorphism> )
        sage: type(_)
        <class 'sage.libs.gap.element.GapElement'>

        sage: F = GF(7); MS = MatrixSpace(F,2,2)
        sage: F.multiplicative_generator()
        3
        sage: G = MatrixGroup([MS([3,0,0,1])])
        sage: a = G.gens()[0]^2
        sage: phi = G.hom([a])
        sage: phi.kernel()
        Subgroup with 1 generators (
        [6 0]
        [0 1]
        ) of Matrix group over Finite Field of size 7 with 1 generators (
        [3 0]
        [0 1]
        )

        sage: F = GF(7); MS = MatrixSpace(F,2,2)
        sage: F.multiplicative_generator()
        3
        sage: G = MatrixGroup([MS([3,0,0,1])])
        sage: a = G.gens()[0]^2
        sage: phi = G.hom([a])
        sage: phi
        Group endomorphism of Matrix group over Finite Field of size 7 with 1 generators (
        [3 0]
        [0 1]
        )
        sage: g = G.gens()[0]
        sage: phi(g)
        [2 0]
        [0 1]
        sage: H = MatrixGroup([MS(a.list())])
        sage: H
        Matrix group over Finite Field of size 7 with 1 generators (
        [2 0]
        [0 1]
        )

    The following tests against :trac:`10659`::

        sage: phi(H)   # indirect doctest
        Subgroup with 1 generators (
        [4 0]
        [0 1]
        ) of Matrix group over Finite Field of size 7 with 1 generators (
        [3 0]
        [0 1]
        )
        sage: F = GF(5); MS = MatrixSpace(F,2,2)
        sage: g = MS([1,1,0,1])
        sage: G = MatrixGroup([g])
        sage: phi = G.hom(G.gens())
        sage: phi(G.0)
        [1 1]
        [0 1]
        sage: phi(G(g^2))
        [1 2]
        [0 1]

        sage: F = GF(5); MS = MatrixSpace(F,2,2)
        sage: gens = [MS([1,2,  -1,1]),MS([1,1,  0,1])]
        sage: G = MatrixGroup(gens)
        sage: phi = G.hom(G.gens())
        sage: phi(G.0)
        [1 2]
        [4 1]
        sage: phi(G.1)
        [1 1]
        [0 1]

    We check that :trac:`19780` is fixed::

        sage: G = groups.matrix.SO(3, 3)
        sage: H = groups.matrix.GL(3, 3)
        sage: phi = G.hom([H(x) for x in G.gens()])
        sage: phi(G.one()).parent()
        General Linear Group of degree 3 over Finite Field of size 3

        sage: MS = MatrixSpace(SR, 2, 2)
        sage: G = MatrixGroup([MS(1), MS([1,2,3,4])])
        sage: G.Hom(G)
        Set of Morphisms from Matrix group over Symbolic Ring with 2 generators (
        [1 0]  [1 2]
        [0 1], [3 4]
        ) to Matrix group over Symbolic Ring with 2 generators (
        [1 0]  [1 2]
        [0 1], [3 4]
        ) in Category of groups

        sage: G = MatrixGroup([matrix(GF(5), [[1,3],[0,1]])])
        sage: H = MatrixGroup([matrix(GF(5), [[1,2],[0,1]])])
        sage: G.hom([H.gen(0)])
        Group morphism:
        From: Matrix group over Finite Field of size 5 with 1 generators (
        [1 3]
        [0 1]
        )
        To:   Matrix group over Finite Field of size 5 with 1 generators (
        [1 2]
        [0 1]
        )

        sage: G = GO(2,2,e=1)
        sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
        sage: A = AbelianGroupGap([2])
        sage: G.hom([A.one(),A.gen(0)])
        Group morphism:
        From: General Orthogonal Group of degree 2 and form parameter 1 over Finite Field of size 2
        To:   Abelian group with gap, generator orders (2,)

    Check that :trac:`19407` is fixed::

        sage: G = GL(2, GF(2))
        sage: H = GL(3, ZZ)
        sage: Hom(G, H)
        Set of Morphisms from General Linear Group of degree 2
         over Finite Field of size 2
         to General Linear Group of degree 3 over Integer Ring in Category of groups
    """
    def __init__(self, homset, gap_hom, check=True):
        r"""
        Constructor method.

        TESTS::

            sage: G = GL(2, ZZ)
            sage: H = GL(2, GF(2))
            sage: P = Hom(G, H)
            sage: gen1 = [g.gap() for g in G.gens()]
            sage: gen2 = [H(g).gap() for g in G.gens()]
            sage: phi = G.gap().GroupHomomorphismByImagesNC(H,gen1, gen2)
            sage: phi = P.element_class(P,phi)
            sage: phi(G.gen(0))
            [0 1]
            [1 0]
        """
        if check:
            if not gap_hom.IsGroupHomomorphism():
                raise ValueError("not a group homomorphism")
            if homset.domain().gap() != gap_hom.Source():
                raise ValueError("domains do not agree")
            if homset.codomain().gap() != gap_hom.Range():
                raise ValueError("ranges do not agree")
        Morphism.__init__(self, homset)
        self._phi = gap_hom


    def __reduce__(self):
        r"""
        Implements pickling.

        We have to work around the fact that GAP does not provide pickling.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: G = AbelianGroupGap([3,2,5])
            sage: f = G.hom(G, G.gens())
            sage: f == loads(dumps(f))
            True
        """
        return self.parent(), (tuple(self(g) for g in self.domain().gens()),)

    def _repr_type(self):
        r"""
        Part of the implementation of :meth:`_repr_`.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([2, 4])
            sage: A.hom([g^2 for g in A.gens()]) # indirect doctest
            Group endomorphism of Abelian group with gap, generator orders (2, 4)
        """
        return "Group"

    def gap(self):
        r"""
        Return the underlying LibGAP group homomorphism.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([2,4])
            sage: f = A.hom([g^2 for g in A.gens()])
            sage: f.gap()
            [ f1, f2 ] -> [ <identity> of ..., f3 ]
        """
        return self._phi

    def _latex_(self):
        r"""
        Return the latex representation of ``self``.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([2,4])
            sage: f = A.hom([g^2 for g in A.gens()])
            sage: latex(f)
            \text{\texttt{Abelian group with gap, generator orders }} \left(2, 4\right)
             \rightarrow \text{\texttt{Abelian group with gap, generator orders }} \left(2, 4\right)
        """
        return r"{} \rightarrow {}".format(latex(self.domain()), latex(self.codomain()))

    def kernel(self):
        r"""
        Return the kernel of ``self``.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A1 = AbelianGroupGap([6, 6])
            sage: A2 = AbelianGroupGap([3, 3])
            sage: f = A1.hom(A2.gens())
            sage: f.kernel()
            Subgroup of Abelian group with gap, generator orders (6, 6)
             generated by (f1*f2, f3*f4)
            sage: f.kernel().order()
            4
            sage: S = Sp(6,3)
            sage: P = PSp(6,3)
            sage: pr = Hom(S, P).natural_map()
            sage: pr.kernel()
            Subgroup with 1 generators (
            [2 0 0 0 0 0]
            [0 2 0 0 0 0]
            [0 0 2 0 0 0]
            [0 0 0 2 0 0]
            [0 0 0 0 2 0]
            [0 0 0 0 0 2]
            ) of Symplectic Group of degree 6 over Finite Field of size 3
        """
        dom = self.domain()
        return dom._subgroup_constructor(self.gap().Kernel())

    def pushforward(self, J, *args, **kwds):
        r"""
        The image of an element or a subgroup.

        INPUT:

        - ``J`` -- a subgroup or an element of the domain of ``self``

        OUTPUT:

        The image of ``J`` under ``self``.

        .. NOTE::

            ``pushforward`` is the method that is used when a map is called
            on anything that is not an element of its domain. For historical
            reasons, we keep the alias ``image()`` for this method.

        EXAMPLES::

            sage: G.<a,b> = FreeGroup()
            sage: H = G / (G([1]), G([2])^3)
            sage: f = G.hom(H.gens())
            sage: S = G.subgroup([a.gap()])
            sage: f.pushforward(S)
            Group([ a ])
            sage: x = f.image(a)
            sage: x
            a
            sage: x.parent()
            Finitely presented group < a, b | a, b^3 >
            sage: G = GU(3,2)
            sage: P = PGU(3,2)
            sage: pr = Hom(G, P).natural_map()
            sage: GS = G.subgroup([G.gen(0)])
            sage: pr.pushforward(GS)
            Subgroup generated by [(3,4,5)(10,18,14)(11,19,15)(12,20,16)(13,21,17)] of (The projective general unitary group of degree 3 over Finite Field of size 2)
        """
        dom = self.domain()
        codom = self.codomain()
        phi = self.gap()
        if isinstance(J, dom.Element) and J in dom:
            return self._call_(dom(J))
        from sage.groups.perm_gps.permgroup import PermutationGroup_generic
        if not isinstance(J, (ParentLibGAP, PermutationGroup_generic)):
            raise TypeError("J (={}) must be a libgap or permutation group".format(J))
        if dom.gap().IsSubgroup(J.gap()).sage():
            return codom._subgroup_constructor(phi.Image(J.gap()))

    image = pushforward

    def _call_(self, g):
        """
        Call syntax for morphisms.

        Some python code for wrapping GAP's ``Images`` function.

        INPUT:

        - ``g`` -- an element of the domain

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([2, 4])
            sage: f = A.hom([g^2 for g in A.gens()])
            sage: a = A.gens()[1]
            sage: f(a)
            f3

        TESTS:

        The following tests we do fall back behind :trac:`10659`::

            sage: O = WeylGroup(['D',6])
            sage: r = prod(O.gens())
            sage: r_ = r^-1
            sage: f = O.hom([r*x*r_ for x in O.gens()])  # long time (19s on sage.math, 2011)
            sage: [f(x) for x in O.gens()]  # long time
            [
            [1 0 0 0 0 0]  [1 0 0 0 0 0]  [1 0 0 0 0 0]  [ 0  0  0  0 -1  0]
            [0 0 1 0 0 0]  [0 1 0 0 0 0]  [0 1 0 0 0 0]  [ 0  1  0  0  0  0]
            [0 1 0 0 0 0]  [0 0 0 1 0 0]  [0 0 1 0 0 0]  [ 0  0  1  0  0  0]
            [0 0 0 1 0 0]  [0 0 1 0 0 0]  [0 0 0 0 1 0]  [ 0  0  0  1  0  0]
            [0 0 0 0 1 0]  [0 0 0 0 1 0]  [0 0 0 1 0 0]  [-1  0  0  0  0  0]
            [0 0 0 0 0 1], [0 0 0 0 0 1], [0 0 0 0 0 1], [ 0  0  0  0  0  1],
            <BLANKLINE>
            [0 0 0 0 0 1]  [ 0  0  0  0  0 -1]
            [0 1 0 0 0 0]  [ 0  1  0  0  0  0]
            [0 0 1 0 0 0]  [ 0  0  1  0  0  0]
            [0 0 0 1 0 0]  [ 0  0  0  1  0  0]
            [0 0 0 0 1 0]  [ 0  0  0  0  1  0]
            [1 0 0 0 0 0], [-1  0  0  0  0  0]
            ]
            sage: f(O)  # long time
            Subgroup with 6 generators of Weyl Group of type ['D', 6] (as a matrix group acting on the ambient space)
            sage: f(O).gens()   # long time
            (
            [1 0 0 0 0 0]  [1 0 0 0 0 0]  [1 0 0 0 0 0]  [ 0  0  0  0 -1  0]
            [0 0 1 0 0 0]  [0 1 0 0 0 0]  [0 1 0 0 0 0]  [ 0  1  0  0  0  0]
            [0 1 0 0 0 0]  [0 0 0 1 0 0]  [0 0 1 0 0 0]  [ 0  0  1  0  0  0]
            [0 0 0 1 0 0]  [0 0 1 0 0 0]  [0 0 0 0 1 0]  [ 0  0  0  1  0  0]
            [0 0 0 0 1 0]  [0 0 0 0 1 0]  [0 0 0 1 0 0]  [-1  0  0  0  0  0]
            [0 0 0 0 0 1], [0 0 0 0 0 1], [0 0 0 0 0 1], [ 0  0  0  0  0  1],
            <BLANKLINE>
            [0 0 0 0 0 1]  [ 0  0  0  0  0 -1]
            [0 1 0 0 0 0]  [ 0  1  0  0  0  0]
            [0 0 1 0 0 0]  [ 0  0  1  0  0  0]
            [0 0 0 1 0 0]  [ 0  0  0  1  0  0]
            [0 0 0 0 1 0]  [ 0  0  0  0  1  0]
            [1 0 0 0 0 0], [-1  0  0  0  0  0]
            )
        """
        img_gap = self.gap().Image(g.gap())
        return self.codomain()(img_gap)

    def lift(self, h):
        r"""
        Return an element of the domain that maps to ``h``.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([2,4])
            sage: f = A.hom([g^2 for g in A.gens()])
            sage: a = A.gens()[1]
            sage: f.lift(a^2)
            f2

        If the element is not in the image, we raise an error::

            sage: f.lift(a)
            Traceback (most recent call last):
            ...
            ValueError: f2 is not an element of the image of Group endomorphism
             of Abelian group with gap, generator orders (2, 4)
        """
        if h not in self.codomain():
            raise TypeError("h (={}) must be an element of the codomain".format(h))
        h = self.codomain()(h)
        phi = self.gap()
        if h.gap() not in phi.Image():
            raise ValueError("{} is not an element of the image of {}".format(h, self))
        return self.domain()( phi.PreImagesRepresentative(h.gap()) )

    def preimage(self, S):
        r"""
        Return the preimage of the subgroup ``S``.

        INPUT:

        - ``S`` -- a subgroup of this group

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([2,4])
            sage: B = AbelianGroupGap([4])
            sage: f = A.hom([B.one(), B.gen(0)^2])
            sage: S = B.subgroup([B.one()])
            sage: f.preimage(S) == f.kernel()
            True
            sage: S = Sp(4,3)
            sage: P = PSp(4,3)
            sage: pr = Hom(S, P).natural_map()
            sage: PS = P.subgroup([P.gen(0)])
            sage: pr.preimage(PS)
            Subgroup with 2 generators (
            [2 0 0 0]  [1 0 0 0]
            [0 2 0 0]  [0 2 0 0]
            [0 0 2 0]  [0 0 2 0]
            [0 0 0 2], [0 0 0 1]
            ) of Symplectic Group of degree 4 over Finite Field of size 3
        """
        phi = self.gap()
        from sage.groups.perm_gps.permgroup import PermutationGroup_generic
        if not isinstance(S, (ParentLibGAP, PermutationGroup_generic)):
            raise TypeError("%s must be a GAP or permutation group of %s"%(S, self))
        if not self.codomain().gap().IsSubgroup(S.gap()).sage():
            raise ValueError("%s must be a subgroup of %s"%(S, self))
        preimage = phi.PreImage(S.gap())
        return self.domain()._subgroup_constructor(preimage)


    def section(self):
        r"""
        This method returns a section map of self by use of :meth:`lift`.
        See :meth:`section` of :class:`sage.categories.map.Map`, as well.

        OUTPUT:

        an instance of :class:`sage.categories.morphism.SetMorphism`
        mapping an element of the codomain of self to one of its preimages

        EXAMPLES::

            sage: G = GU(3,2)
            sage: P = PGU(3,2)
            sage: pr = Hom(G, P).natural_map()
            sage: sect = pr.section()
            sage: sect(P.an_element())
            [a + 1     a     a]
            [    1     1     0]
            [    a     0     0]
        """
        from sage.categories.homset import Hom
        from sage.categories.sets_cat import Sets
        H = Hom(self.codomain(), self.domain(), category=Sets())
        return H(self.lift)


class GroupHomset_libgap(HomsetWithBase):
    r"""
    Homsets of groups with a libgap backend.

    Do not call this directly instead use :func:`Hom`.

    INPUT:

    - ``G`` -- a libgap group
    - ``H`` -- a libgap group
    - ``category`` -- a category

    OUTPUT:

    The homset of two libgap groups.

    EXAMPLES::

        sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
        sage: A = AbelianGroupGap([2,4])
        sage: H = A.Hom(A)
        sage: H
        Set of Morphisms from Abelian group with gap, generator orders (2, 4)
         to Abelian group with gap, generator orders (2, 4)
         in Category of finite enumerated commutative groups
    """
    def __init__(self, G, H, category=None, check=True):
        r"""
        Return the homset of two libgap groups.

        TESTS::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([2,4])
            sage: H = A.Hom(A)
            sage: TestSuite(H).run()
        """
        if check:
            from sage.groups.perm_gps.permgroup import PermutationGroup_generic
            if not isinstance(G, (ParentLibGAP, PermutationGroup_generic)):
                raise TypeError("G (={}) must be a ParentLibGAP or a permutation group".format(G))
            if not isinstance(H, (ParentLibGAP, PermutationGroup_generic)):
                raise TypeError("H (={}) must be a ParentLibGAP or a permutation group".format(H))
        HomsetWithBase.__init__(self, G, H, category, check=check, base=ZZ)

    Element = GroupMorphism_libgap


    def _element_constructor_(self, x, check=True, **options):
        r"""
        Handle conversions and coercions.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([2,4])
            sage: H = A.Hom(A)
            sage: H([g^2 for g in A.gens()])
            Group endomorphism of Abelian group with gap, generator orders (2, 4)
            sage: H(2)
            Traceback (most recent call last):
            ...
            TypeError: unable to convert 2 to an element of ...

        A group homomorphism between a finitely presented group and a subgroup of a permutation group::

            sage: PG = PGU(6,2)
            sage: g, h = PG.gens()
            sage: p1 = h^-3*(h^-1*g^-1)^2*h*g*h^2*g^-1*h^2*g*h^-5*g^-1
            sage: p2 = g*(g*h)^2*g*h^-4*(g*h)^2*(h^2*g*h^-2*g)^2*h^-2*g*h^-2*g^-1*h^-1*g*h*g*h^-1*g
            sage: p3 = h^-3*g^-1*h*g*h^4*g^-1*h^-1*g*h*(h^2*g^-1)^2*h^-4*g*h^2*g^-1*h^-7*g^-2*h^-2*g*h^-2*g^-1*h^-1*(g*h)^2*h^3
            sage: p4 = h*(h^3*g)^2*h*g*h^-1*g*h^2*g^-1*h^-2*g*h^4*g^-1*h^3*g*h^-2*g*h^-1*g^-1*h^2*g*h*g^-1*h^-2*g*h*g^-1*h^2*g*h^2*g^-1
            sage: p5 = h^2*g*h^2*g^-1*h*g*h^-1*g*h*g^-1*h^2*g*h^-2*g*h^2*g*h^-2*(h^-1*g)^2*h^4*(g*h^-1)^2*g^-1
            sage: UPG = PG.subgroup([p1, p2, p3, p4, p5], canonicalize=False)
            sage: B6 = BraidGroup(6)
            sage: reprB6 = B6.hom(UPG.gens())
            sage: b1, b2, b3, b4, b5 = B6.gens()
            sage: reprB6(b1*b2*b3*b4*b5) == p1*p2*p3*p4*p5
            True

        TESTS:

        The following input does not define a valid group homomorphism::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([2])
            sage: G = MatrixGroup([Matrix(ZZ, 2, [0,1,-1,0])])
            sage: A.hom([G.gen(0)])
            Traceback (most recent call last):
            ...
            ValueError: images do not define a group homomorphism
        """
        if isinstance(x, (tuple, list)):
            # there should be a better way
            from sage.libs.gap.libgap import libgap
            dom = self.domain()
            codom = self.codomain()
            gens = dom.gap().GeneratorsOfGroup()
            imgs = [codom(g).gap() for g in x]
            if check:
                if not len(gens) == len(imgs):
                    raise ValueError("provide an image for each generator")
                phi = libgap.GroupHomomorphismByImages(dom.gap(), codom.gap(), gens, imgs)
                # if it is not a group homomorphism, then
                # self._phi is the gap boolean fail
                if phi.is_bool():     # check we did not fail
                    raise ValueError("images do not define a group homomorphism")
            else:
                ByImagesNC = libgap.function_factory("GroupHomomorphismByImagesNC")
                phi = ByImagesNC(dom.gap(), codom.gap(), gens, imgs)
            return self.element_class(self, phi, check=check, **options)
        if isinstance(x, GapElement):
            try:
                return self.element_class(self, x, check=True, **options)
            except ValueError:
                pass
        return super(GroupHomset_libgap, self)._element_constructor_(x, check=check, **options)

    def _an_element_(self):
        r"""
        Return an element of ``self``.

        There is only one group homomorphism which is always defined:
        the trivial one.

        EXAMPLES::

            sage: from sage.groups.abelian_gps.abelian_group_gap import AbelianGroupGap
            sage: A = AbelianGroupGap([2,4])
            sage: G.<a,b> = FreeGroup()
            sage: A.Hom(G).an_element()
            Group morphism:
              From: Abelian group with gap, generator orders (2, 4)
              To:   Free Group on generators {a, b}
        """
        n = len(self.domain().gens())
        one = self.codomain().one()
        return self._element_constructor_((one,) * n, check=False)

    def natural_map(self):
        r"""
        This method from :class:`HomsetWithBase` is overloaded here for cases in which
        both groups have corresponding lists of generators.

        OUTPUT:

        an instance of the element class of self if there exists a group homomorphism
        mapping the generators of the domain of self to the according generators of
        the codomain. Else the method falls back to the default.

        EXAMPLES::

            sage: G = GL(3,2)
            sage: P = PGL(3,2)
            sage: nat = Hom(G, P).natural_map()
            sage: type(nat)
            <class 'sage.groups.libgap_morphism.GroupHomset_libgap_with_category.element_class'>
            sage: g1, g2 = G.gens()
            sage: nat(g1*g2)
            (1,2,4,5,7,3,6)
        """
        if len(self.domain().gens()) == len(self.codomain().gens()):
            dom_gap = self.domain().gap()
            codom_gap = self.codomain().gap()
            from sage.libs.gap.libgap import libgap
            phi = libgap.GroupHomomorphismByImages(dom_gap,
                                                   codom_gap,
                                                   dom_gap.GeneratorsOfGroup(),
                                                   codom_gap.GeneratorsOfGroup()
                                                   )
            if not phi.is_bool():     # phi is indeed a group homomorphism
                return self.element_class(self, phi)
        return super(GroupHomset_libgap, self).natural_map()
