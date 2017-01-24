"""
Mix-in Class for libGAP-based Groups

This class adds access to GAP functionality to groups such that parent
and element have a ``gap()`` method that returns a libGAP object for
the parent/element.

If your group implementation uses libgap, then you should add
:class:`GroupMixinLibGAP` as the first class that you are deriving
from. This ensures that it properly overrides any default methods that
just raise ``NotImplemented``.
"""

from sage.libs.all import libgap
from sage.misc.cachefunc import cached_method

class GroupMixinLibGAP(object):

    @cached_method
    def is_abelian(self):
        r"""
        Test whether the group is Abelian.

        OUTPUT:

        Boolean. ``True`` if this group is an Abelian group.

        EXAMPLES::

            sage: SL(1, 17).is_abelian()
            True
            sage: SL(2, 17).is_abelian()
            False
        """
        return self.gap().IsAbelian().sage()

    @cached_method
    def is_finite(self):
        """
        Test whether the matrix group is finite.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: G = GL(2,GF(3))
            sage: G.is_finite()
            True
            sage: SL(2,ZZ).is_finite()
            False
        """
        return self.gap().IsFinite().sage()

    def cardinality(self):
        """
        Implements :meth:`EnumeratedSets.ParentMethods.cardinality`.

        EXAMPLES::

            sage: G = Sp(4,GF(3))
            sage: G.cardinality()
            51840

            sage: G = SL(4,GF(3))
            sage: G.cardinality()
            12130560

            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,2],[-1,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: G.cardinality()
            480

            sage: G = MatrixGroup([matrix(ZZ,2,[1,1,0,1])])
            sage: G.cardinality()
            +Infinity

            sage: G = Sp(4,GF(3))
            sage: G.cardinality()
            51840

            sage: G = SL(4,GF(3))
            sage: G.cardinality()
            12130560

            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,2],[-1,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: G.cardinality()
            480

            sage: G = MatrixGroup([matrix(ZZ,2,[1,1,0,1])])
            sage: G.cardinality()
            +Infinity
        """
        return self.gap().Size().sage()
        from sage.rings.infinity import Infinity
        return Infinity

    order = cardinality

    @cached_method
    def conjugacy_class_representatives(self):
        """
        Return a set of representatives for each of the conjugacy classes
        of the group.

        EXAMPLES::

            sage: G = SU(3,GF(2))
            sage: len(G.conjugacy_class_representatives())
            16

            sage: G = GL(2,GF(3))
            sage: G.conjugacy_class_representatives()
            (
            [1 0]  [0 2]  [2 0]  [0 2]  [0 2]  [0 1]  [0 1]  [2 0]
            [0 1], [1 1], [0 2], [1 2], [1 0], [1 2], [1 1], [0 1]
            )

            sage: len(GU(2,GF(5)).conjugacy_class_representatives())
            36
        """
        G = self.gap()
        reps = [ cc.Representative() for cc in G.ConjugacyClasses() ]
        return tuple(self(g) for g in reps)

    def conjugacy_classes(self):
        r"""
        Return a list with all the conjugacy classes of ``self``.

        EXAMPLES::

            sage: G = SL(2, GF(2))
            sage: G.conjugacy_classes()
            (Conjugacy class of [1 0]
             [0 1] in Special Linear Group of degree 2 over Finite Field of size 2,
             Conjugacy class of [0 1]
             [1 0] in Special Linear Group of degree 2 over Finite Field of size 2,
             Conjugacy class of [0 1]
             [1 1] in Special Linear Group of degree 2 over Finite Field of size 2)
        """
        from sage.groups.conjugacy_classes import ConjugacyClassGAP
        return tuple(ConjugacyClassGAP(self, self(g)) for g in self.conjugacy_class_representatives())

    def conjugacy_class(self, g):
        r"""
        Return the conjugacy class of ``g``.

        OUTPUT:

        The conjugacy class of ``g`` in the group ``self``. If ``self`` is the
        group denoted by `G`, this method computes the set
        `\{x^{-1}gx\ \vert\ x\in G\}`.

        EXAMPLES::

            sage: G = SL(2, QQ)
            sage: g = G([[1,1],[0,1]])
            sage: G.conjugacy_class(g)
            Conjugacy class of [1 1]
            [0 1] in Special Linear Group of degree 2 over Rational Field
        """
        from sage.groups.conjugacy_classes import ConjugacyClassGAP
        return ConjugacyClassGAP(self, self(g))

    def class_function(self, values):
        """
        Return the class function with given values.

        INPUT:

        - ``values`` -- list/tuple/iterable of numbers. The values of the
          class function on the conjugacy classes, in that order.

        EXAMPLES::

            sage: G = GL(2,GF(3))
            sage: chi = G.class_function(range(8))
            sage: list(chi)
            [0, 1, 2, 3, 4, 5, 6, 7]
        """
        from sage.groups.class_function import ClassFunction_libgap
        return ClassFunction_libgap(self, values)

    @cached_method
    def center(self):
        """
        Return the center of this linear group as a subgroup.

        OUTPUT:

        The center as a subgroup.

        EXAMPLES::

            sage: G = SU(3,GF(2))
            sage: G.center()
            Matrix group over Finite Field in a of size 2^2 with 1 generators (
            [a 0 0]
            [0 a 0]
            [0 0 a]
            )
            sage: GL(2,GF(3)).center()
            Matrix group over Finite Field of size 3 with 1 generators (
            [2 0]
            [0 2]
            )
            sage: GL(3,GF(3)).center()
            Matrix group over Finite Field of size 3 with 1 generators (
            [2 0 0]
            [0 2 0]
            [0 0 2]
            )
            sage: GU(3,GF(2)).center()
            Matrix group over Finite Field in a of size 2^2 with 1 generators (
            [a + 1     0     0]
            [    0 a + 1     0]
            [    0     0 a + 1]
            )

            sage: A = Matrix(FiniteField(5), [[2,0,0], [0,3,0], [0,0,1]])
            sage: B = Matrix(FiniteField(5), [[1,0,0], [0,1,0], [0,1,1]])
            sage: MatrixGroup([A,B]).center()
            Matrix group over Finite Field of size 5 with 1 generators (
            [1 0 0]
            [0 1 0]
            [0 0 1]
            )
        """
        G = self.gap()
        center = list(G.Center().GeneratorsOfGroup())
        if len(center) == 0:
            center = [G.One()]
        return self.subgroup(center)

    def intersection(self, other):
        """
        Return the intersection of two groups (if it makes sense) as a
        subgroup of the first group.

        EXAMPLES::

            sage: A = Matrix([(0, 1/2, 0), (2, 0, 0), (0, 0, 1)])
            sage: B = Matrix([(0, 1/2, 0), (-2, -1, 2), (0, 0, 1)])
            sage: G = MatrixGroup([A,B])
            sage: len(G)  # isomorphic to S_3
            6
            sage: G.intersection(GL(3,ZZ))
            Matrix group over Rational Field with 1 generators (
            [ 1  0  0]
            [-2 -1  2]
            [ 0  0  1]
            )
            sage: GL(3,ZZ).intersection(G)
            Matrix group over Integer Ring with 1 generators (
            [ 1  0  0]
            [-2 -1  2]
            [ 0  0  1]
            )
            sage: G.intersection(SL(3,ZZ))
            Matrix group over Rational Field with 0 generators ()
        """
        G = self.gap()
        H = other.gap()
        C = G.Intersection(H)
        return self.subgroup(C.GeneratorsOfGroup())

    @cached_method
    def irreducible_characters(self):
        """
        Returns the irreducible characters of the group.

        OUTPUT:

        A tuple containing all irreducible characters.

        EXAMPLES::

            sage: G = GL(2,2)
            sage: G.irreducible_characters()
            (Character of General Linear Group of degree 2 over Finite Field of size 2,
             Character of General Linear Group of degree 2 over Finite Field of size 2,
             Character of General Linear Group of degree 2 over Finite Field of size 2)
        """
        Irr = self.gap().Irr()
        L = []
        from sage.groups.class_function import ClassFunction_libgap
        for irr in Irr:
            L.append(ClassFunction_libgap(self, irr))
        return tuple(L)

    def random_element(self):
        """
        Return a random element of this group.

        OUTPUT:

        A group element.

        EXAMPLES::

            sage: G = Sp(4,GF(3))
            sage: G.random_element()  # random
            [2 1 1 1]
            [1 0 2 1]
            [0 1 1 0]
            [1 0 0 1]
            sage: G.random_element() in G
            True

            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,2],[-1,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: G.random_element()  # random
            [1 3]
            [0 3]
            sage: G.random_element() in G
            True
        """
        return self(self.gap().Random())

    def __iter__(self):
        """
        Iterate over the elements of the group.

        EXAMPLES::

            sage: F = GF(3)
            sage: gens = [matrix(F,2, [1,0, -1,1]), matrix(F, 2, [1,1,0,1])]
            sage: G = MatrixGroup(gens)
            sage: next(iter(G))
            [1 0]
            [0 1]
        """
        if self.list.cache is not None:
            for g in self.list():
                yield g
            return
        iterator = self.gap().Iterator()
        while not iterator.IsDoneIterator().sage():
            yield self(iterator.NextIterator(), check=False)

    def __len__(self):
        """
        Return the number of elements in ``self``.

        EXAMPLES::

            sage: F = GF(3)
            sage: gens = [matrix(F,2, [1,-1,0,1]), matrix(F, 2, [1,1,-1,1])]
            sage: G = MatrixGroup(gens)
            sage: len(G)
            48

        An error is raised if the group is not finite::

            sage: len(GL(2,ZZ))
            Traceback (most recent call last):
            ...
            NotImplementedError: group must be finite
        """
        if not self.is_finite():
            raise NotImplementedError('group must be finite')
        return int(self.cardinality())

    @cached_method
    def list(self):
        """
        List all elements of this group.

        OUTPUT:

        A tuple containing all group elements in a random but fixed
        order.

        EXAMPLES::

            sage: F = GF(3)
            sage: gens = [matrix(F,2, [1,0,-1,1]), matrix(F, 2, [1,1,0,1])]
            sage: G = MatrixGroup(gens)
            sage: G.cardinality()
            24
            sage: v = G.list()
            sage: len(v)
            24
            sage: v[:5]
            (
            [1 0]  [2 0]  [0 1]  [0 2]  [1 2]
            [0 1], [0 2], [2 0], [1 0], [2 2]
            )
            sage: all(g in G for g in G.list())
            True

        An example over a ring (see :trac:`5241`)::

            sage: M1 = matrix(ZZ,2,[[-1,0],[0,1]])
            sage: M2 = matrix(ZZ,2,[[1,0],[0,-1]])
            sage: M3 = matrix(ZZ,2,[[-1,0],[0,-1]])
            sage: MG = MatrixGroup([M1, M2, M3])
            sage: MG.list()
            (
            [1 0]  [ 1  0]  [-1  0]  [-1  0]
            [0 1], [ 0 -1], [ 0  1], [ 0 -1]
            )
            sage: MG.list()[1]
            [ 1  0]
            [ 0 -1]
            sage: MG.list()[1].parent()
            Matrix group over Integer Ring with 3 generators (
            [-1  0]  [ 1  0]  [-1  0]
            [ 0  1], [ 0 -1], [ 0 -1]
            )

        An example over a field (see :trac:`10515`)::

            sage: gens = [matrix(QQ,2,[1,0,0,1])]
            sage: MatrixGroup(gens).list()
            (
            [1 0]
            [0 1]
            )

        Another example over a ring (see :trac:`9437`)::

            sage: len(SL(2, Zmod(4)).list())
            48

        An error is raised if the group is not finite::

            sage: GL(2,ZZ).list()
            Traceback (most recent call last):
            ...
            NotImplementedError: group must be finite
        """
        if not self.is_finite():
            raise NotImplementedError('group must be finite')
        elements = self.gap().Elements()
        return tuple(self(x, check=False) for x in elements)

    def is_isomorphic(self, H):
        """
        Test whether ``self`` and ``H`` are isomorphic groups.

        INPUT:

        - ``H`` -- a group.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: m1 = matrix(GF(3), [[1,1],[0,1]])
            sage: m2 = matrix(GF(3), [[1,2],[0,1]])
            sage: F = MatrixGroup(m1)
            sage: G = MatrixGroup(m1, m2)
            sage: H = MatrixGroup(m2)
            sage: F.is_isomorphic(G)
            True
            sage: G.is_isomorphic(H)
            True
            sage: F.is_isomorphic(H)
            True
            sage: F==G, G==H, F==H
            (False, False, False)
        """
        iso = self.gap().IsomorphismGroups(H.gap())
        if iso.is_bool():   # fail means not isomorphic
            try:
                iso.sage()
                assert False
            except ValueError:
                pass
            return False
        return True
