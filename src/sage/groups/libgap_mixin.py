# sage_setup: distribution = sagemath-gap
# sage.doctest: needs sage.modules
"""
Mix-in Class for GAP-based Groups

This class adds access to GAP functionality to groups such that parent
and element have a ``gap()`` method that returns a GAP object for
the parent/element.

If your group implementation uses libgap, then you should add
:class:`GroupMixinLibGAP` as the first class that you are deriving
from. This ensures that it properly overrides any default methods that
just raise :exc:`NotImplementedError`.
"""

from sage.libs.gap.libgap import libgap
from sage.libs.gap.element import GapElement
from sage.structure.element import parent
from sage.misc.cachefunc import cached_method
from sage.groups.class_function import ClassFunction_libgap
from sage.groups.libgap_wrapper import ElementLibGAP


class GroupMixinLibGAP:
    def __contains__(self, elt):
        r"""
        TESTS::

            sage: from sage.groups.libgap_group import GroupLibGAP
            sage: G = GroupLibGAP(libgap.SL(2,3))
            sage: libgap([[1,0],[0,1]]) in G
            False
            sage: o = Mod(1, 3)
            sage: z = Mod(0, 3)
            sage: libgap([[o,z],[z,o]]) in G
            True

            sage: G.an_element() in GroupLibGAP(libgap.GL(2,3))
            True
            sage: G.an_element() in GroupLibGAP(libgap.GL(2,5))
            False
        """
        if parent(elt) is self:
            return True
        elif isinstance(elt, GapElement):
            return elt in self.gap()
        elif isinstance(elt, ElementLibGAP):
            return elt.gap() in self.gap()
        else:
            try:
                elt2 = self(elt)
            except Exception:
                return False
            return elt == elt2

    def is_abelian(self):
        r"""
        Return whether the group is Abelian.

        OUTPUT: boolean; ``True`` if this group is an Abelian group and
        ``False`` otherwise

        EXAMPLES::

            sage: from sage.groups.libgap_group import GroupLibGAP
            sage: GroupLibGAP(libgap.CyclicGroup(12)).is_abelian()
            True
            sage: GroupLibGAP(libgap.SymmetricGroup(12)).is_abelian()
            False

            sage: SL(1, 17).is_abelian()
            True
            sage: SL(2, 17).is_abelian()
            False
        """
        return self.gap().IsAbelian().sage()

    def is_nilpotent(self):
        r"""
        Return whether this group is nilpotent.

        EXAMPLES::

            sage: from sage.groups.libgap_group import GroupLibGAP
            sage: GroupLibGAP(libgap.AlternatingGroup(3)).is_nilpotent()
            True
            sage: GroupLibGAP(libgap.SymmetricGroup(3)).is_nilpotent()
            False
        """
        return self.gap().IsNilpotentGroup().sage()

    def is_solvable(self):
        r"""
        Return whether this group is solvable.

        EXAMPLES::

            sage: from sage.groups.libgap_group import GroupLibGAP
            sage: GroupLibGAP(libgap.SymmetricGroup(4)).is_solvable()
            True
            sage: GroupLibGAP(libgap.SymmetricGroup(5)).is_solvable()
            False
        """
        return self.gap().IsSolvableGroup().sage()

    def is_supersolvable(self):
        r"""
        Return whether this group is supersolvable.

        EXAMPLES::

            sage: from sage.groups.libgap_group import GroupLibGAP
            sage: GroupLibGAP(libgap.SymmetricGroup(3)).is_supersolvable()
            True
            sage: GroupLibGAP(libgap.SymmetricGroup(4)).is_supersolvable()
            False
        """
        return self.gap().IsSupersolvableGroup().sage()

    def is_polycyclic(self):
        r"""
        Return whether this group is polycyclic.

        EXAMPLES::

            sage: from sage.groups.libgap_group import GroupLibGAP
            sage: GroupLibGAP(libgap.AlternatingGroup(4)).is_polycyclic()
            True
            sage: GroupLibGAP(libgap.AlternatingGroup(5)).is_solvable()
            False
        """
        return self.gap().IsPolycyclicGroup().sage()

    def is_perfect(self):
        r"""
        Return whether this group is perfect.

        EXAMPLES::

            sage: from sage.groups.libgap_group import GroupLibGAP
            sage: GroupLibGAP(libgap.SymmetricGroup(5)).is_perfect()
            False
            sage: GroupLibGAP(libgap.AlternatingGroup(5)).is_perfect()
            True

            sage: SL(3,3).is_perfect()
            True
        """
        return self.gap().IsPerfectGroup().sage()

    def is_p_group(self):
        r"""
        Return whether this group is a p-group.

        EXAMPLES::

            sage: from sage.groups.libgap_group import GroupLibGAP
            sage: GroupLibGAP(libgap.CyclicGroup(9)).is_p_group()
            True
            sage: GroupLibGAP(libgap.CyclicGroup(10)).is_p_group()
            False
        """
        return self.gap().IsPGroup().sage()

    def is_simple(self):
        r"""
        Return whether this group is simple.

        EXAMPLES::

            sage: from sage.groups.libgap_group import GroupLibGAP
            sage: GroupLibGAP(libgap.SL(2,3)).is_simple()
            False
            sage: GroupLibGAP(libgap.SL(3,3)).is_simple()
            True

            sage: SL(3,3).is_simple()
            True
        """
        return self.gap().IsSimpleGroup().sage()

    def is_finite(self):
        """
        Test whether the matrix group is finite.

        OUTPUT: boolean

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
        Implement :meth:`EnumeratedSets.ParentMethods.cardinality`.

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

    order = cardinality

    @cached_method
    def conjugacy_classes_representatives(self):
        """
        Return a set of representatives for each of the conjugacy classes
        of the group.

        EXAMPLES::

            sage: G = SU(3,GF(2))                                                       # needs sage.rings.finite_rings
            sage: len(G.conjugacy_classes_representatives())                            # needs sage.rings.finite_rings
            16

            sage: G = GL(2,GF(3))
            sage: G.conjugacy_classes_representatives()
            (
            [1 0]  [0 2]  [2 0]  [0 2]  [0 2]  [0 1]  [0 1]  [2 0]
            [0 1], [1 1], [0 2], [1 2], [1 0], [1 2], [1 1], [0 1]
            )

            sage: len(GU(2,GF(5)).conjugacy_classes_representatives())                  # needs sage.rings.finite_rings
            36

        ::

            sage: GL(2,ZZ).conjugacy_classes_representatives()
            Traceback (most recent call last):
            ...
            NotImplementedError: only implemented for finite groups
        """
        if not self.is_finite():
            raise NotImplementedError("only implemented for finite groups")
        G = self.gap()
        reps = [cc.Representative() for cc in G.ConjugacyClasses()]
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

        ::

            sage: GL(2,ZZ).conjugacy_classes()
            Traceback (most recent call last):
            ...
            NotImplementedError: only implemented for finite groups
        """
        if not self.is_finite():
            raise NotImplementedError("only implemented for finite groups")
        from sage.groups.conjugacy_classes import ConjugacyClassGAP
        return tuple(ConjugacyClassGAP(self, self(g)) for g in self.conjugacy_classes_representatives())

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

        - ``values`` -- list/tuple/iterable of numbers; the values of the
          class function on the conjugacy classes, in that order

        EXAMPLES::

            sage: G = GL(2,GF(3))
            sage: chi = G.class_function(range(8))                                      # needs sage.rings.number_field
            sage: list(chi)                                                             # needs sage.rings.number_field
            [0, 1, 2, 3, 4, 5, 6, 7]
        """
        from sage.groups.class_function import ClassFunction_libgap
        return ClassFunction_libgap(self, values)

    @cached_method
    def center(self):
        """
        Return the center of this group as a subgroup.

        OUTPUT: the center as a subgroup

        EXAMPLES::

            sage: G = SU(3, GF(2))                                                      # needs sage.rings.finite_rings
            sage: G.center()                                                            # needs sage.rings.finite_rings
            Subgroup with 1 generators (
            [a 0 0]
            [0 a 0]
            [0 0 a]
            ) of Special Unitary Group of degree 3 over Finite Field in a of size 2^2
            sage: GL(2, GF(3)).center()
            Subgroup with 1 generators (
            [2 0]
            [0 2]
            ) of General Linear Group of degree 2 over Finite Field of size 3
            sage: GL(3, GF(3)).center()
            Subgroup with 1 generators (
            [2 0 0]
            [0 2 0]
            [0 0 2]
            ) of General Linear Group of degree 3 over Finite Field of size 3
            sage: GU(3, GF(2)).center()                                                 # needs sage.rings.finite_rings
            Subgroup with 1 generators (
            [a + 1     0     0]
            [    0 a + 1     0]
            [    0     0 a + 1]
            ) of General Unitary Group of degree 3 over Finite Field in a of size 2^2

            sage: A = Matrix(FiniteField(5), [[2,0,0], [0,3,0], [0,0,1]])
            sage: B = Matrix(FiniteField(5), [[1,0,0], [0,1,0], [0,1,1]])
            sage: MatrixGroup([A,B]).center()
            Subgroup with 1 generators (
            [1 0 0]
            [0 1 0]
            [0 0 1]
            ) of Matrix group over Finite Field of size 5 with 2 generators (
            [2 0 0]  [1 0 0]
            [0 3 0]  [0 1 0]
            [0 0 1], [0 1 1]
            )

            sage: GL = groups.matrix.GL(3, ZZ)
            sage: GL.center()
            Traceback (most recent call last):
            ...
            NotImplementedError: group must be finite
        """
        if not self.is_finite():
            raise NotImplementedError("group must be finite")
        G = self.gap()
        center = list(G.Center().GeneratorsOfGroup())
        if not center:
            center = [G.One()]
        return self.subgroup(center)

    def centralizer(self, g):
        r"""
        Return the centralizer of ``g`` in ``self``.

        EXAMPLES::

            sage: G = groups.matrix.GL(2, 3)
            sage: g = G([[1,1], [1,0]])
            sage: C = G.centralizer(g); C
            Subgroup with 3 generators (
            [1 1]  [2 0]  [2 1]
            [1 0], [0 2], [1 1]
            ) of General Linear Group of degree 2 over Finite Field of size 3
            sage: C.order()
            8

            sage: S = G.subgroup([G([[2,0],[0,2]]), G([[0,1],[2,0]])]); S
            Subgroup with 2 generators (
            [2 0]  [0 1]
            [0 2], [2 0]
            ) of General Linear Group of degree 2 over Finite Field of size 3
            sage: G.centralizer(S)
            Subgroup with 3 generators (
            [2 0]  [0 1]  [2 2]
            [0 2], [2 0], [1 2]
            ) of General Linear Group of degree 2 over Finite Field of size 3
            sage: G = GL(3,2)
            sage: all(G.order() == G.centralizer(x).order() * G.conjugacy_class(x).cardinality()
            ....:     for x in G)
            True
            sage: H = groups.matrix.Heisenberg(2)
            sage: H.centralizer(H.an_element())
            Traceback (most recent call last):
            ...
            NotImplementedError: group must be finite
        """
        if not self.is_finite():
            raise NotImplementedError("group must be finite")
        G = self.gap()
        centralizer_gens = list(G.Centralizer(g).GeneratorsOfGroup())
        if not centralizer_gens:
            centralizer_gens = [G.One()]
        return self.subgroup(centralizer_gens)

    def subgroups(self):
        r"""
        Return a list of all the subgroups of ``self``.

        OUTPUT:

        Each possible subgroup of ``self`` is contained once in the returned
        list. The list is in order, according to the size of the subgroups,
        from the trivial subgroup with one element on through up to the whole
        group. Conjugacy classes of subgroups are contiguous in the list.

        .. WARNING::

            For even relatively small groups this method can take a very long
            time to execute, or create vast amounts of output. Likely both.
            Its purpose is instructional, as it can be useful for studying
            small groups.

            For faster results, which still exhibit the structure of
            the possible subgroups, use :meth:`conjugacy_classes_subgroups`.

        EXAMPLES::

            sage: G = groups.matrix.GL(2, 2)
            sage: G.subgroups()
            [Subgroup with 0 generators ()
               of General Linear Group of degree 2 over Finite Field of size 2,
             Subgroup with 1 generators (
             [0 1]
             [1 0]
             ) of General Linear Group of degree 2 over Finite Field of size 2,
             Subgroup with 1 generators (
             [1 0]
             [1 1]
             ) of General Linear Group of degree 2 over Finite Field of size 2,
             Subgroup with 1 generators (
             [1 1]
             [0 1]
             ) of General Linear Group of degree 2 over Finite Field of size 2,
             Subgroup with 1 generators (
             [0 1]
             [1 1]
             ) of General Linear Group of degree 2 over Finite Field of size 2,
             Subgroup with 2 generators (
             [0 1]  [1 1]
             [1 1], [0 1]
             ) of General Linear Group of degree 2 over Finite Field of size 2]

            sage: H = groups.matrix.Heisenberg(2)
            sage: H.subgroups()
            Traceback (most recent call last):
            ...
            NotImplementedError: group must be finite
        """
        if not self.is_finite():
            raise NotImplementedError("group must be finite")
        ccs = self.gap().ConjugacyClassesSubgroups()
        return [self.subgroup(h.GeneratorsOfGroup())
                for cc in ccs for h in cc.Elements()]

    def conjugacy_classes_subgroups(self):
        r"""
        Return a complete list of representatives of conjugacy classes of
        subgroups in ``self``.

        The ordering is that given by GAP.

        EXAMPLES::

            sage: G = groups.matrix.GL(2,2)
            sage: G.conjugacy_classes_subgroups()
            [Subgroup with 0 generators ()
               of General Linear Group of degree 2 over Finite Field of size 2,
             Subgroup with 1 generators (
             [1 1]
             [0 1]
             ) of General Linear Group of degree 2 over Finite Field of size 2,
             Subgroup with 1 generators (
             [0 1]
             [1 1]
             ) of General Linear Group of degree 2 over Finite Field of size 2,
             Subgroup with 2 generators (
             [0 1]  [1 1]
             [1 1], [0 1]
             ) of General Linear Group of degree 2 over Finite Field of size 2]

            sage: H = groups.matrix.Heisenberg(2)
            sage: H.conjugacy_classes_subgroups()
            Traceback (most recent call last):
            ...
            NotImplementedError: group must be finite
        """
        if not self.is_finite():
            raise NotImplementedError("group must be finite")
        return [self.subgroup(sub.Representative().GeneratorsOfGroup())
                for sub in self.gap().ConjugacyClassesSubgroups()]

    def group_id(self):
        r"""
        Return the ID code of ``self``, which is a list of two integers.

        It is a unique identified assigned by GAP for groups in the
        ``SmallGroup`` library.

        EXAMPLES::

            sage: PGL(2,3).group_id()
            [24, 12]
            sage: SymmetricGroup(4).group_id()
            [24, 12]

            sage: G = groups.matrix.GL(2, 2)
            sage: G.group_id()
            [6, 1]
            sage: G = groups.matrix.GL(2, 3)
            sage: G.id()
            [48, 29]

            sage: G = groups.matrix.GL(2, ZZ)
            sage: G.group_id()
            Traceback (most recent call last):
            ...
            GAPError: Error, the group identification for groups of size infinity is not available
        """
        from sage.rings.integer import Integer
        return [Integer(n) for n in self.gap().IdGroup()]

    id = group_id

    def exponent(self):
        r"""
        Compute the exponent of the group.

        The exponent `e` of a group `G` is the LCM of the orders of its
        elements, that is, `e` is the smallest integer such that `g^e = 1`
        for all `g \in G`.

        EXAMPLES::

            sage: G = groups.matrix.GL(2, 3)
            sage: G.exponent()
            24

            sage: H = groups.matrix.Heisenberg(2)
            sage: H.exponent()
            Traceback (most recent call last):
            ...
            NotImplementedError: group must be finite
        """
        if not self.is_finite():
            raise NotImplementedError("group must be finite")
        from sage.rings.integer import Integer
        return Integer(self._libgap_().Exponent())

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
            Subgroup with 1 generators (
            [ 1  0  0]
            [-2 -1  2]
            [ 0  0  1]
            ) of Matrix group over Rational Field with 2 generators (
              [  0 1/2   0]  [  0 1/2   0]
              [  2   0   0]  [ -2  -1   2]
              [  0   0   1], [  0   0   1]
              )
            sage: GL(3,ZZ).intersection(G)
            Subgroup with 1 generators (
            [ 1  0  0]
            [-2 -1  2]
            [ 0  0  1]
            ) of General Linear Group of degree 3 over Integer Ring
            sage: G.intersection(SL(3,ZZ))
            Subgroup with 0 generators ()
              of Matrix group over Rational Field with 2 generators (
              [  0 1/2   0]  [  0 1/2   0]
              [  2   0   0]  [ -2  -1   2]
              [  0   0   1], [  0   0   1]
              )
        """
        G = self.gap()
        H = other.gap()
        C = G.Intersection(H)
        return self.subgroup(C.GeneratorsOfGroup())

    @cached_method
    def irreducible_characters(self):
        """
        Return the irreducible characters of the group.

        OUTPUT: tuple containing all irreducible characters

        EXAMPLES::

            sage: G = GL(2,2)
            sage: G.irreducible_characters()                                            # needs sage.rings.number_field
            (Character of General Linear Group of degree 2 over Finite Field of size 2,
             Character of General Linear Group of degree 2 over Finite Field of size 2,
             Character of General Linear Group of degree 2 over Finite Field of size 2)

        ::

            sage: GL(2,ZZ).irreducible_characters()
            Traceback (most recent call last):
            ...
            NotImplementedError: only implemented for finite groups
        """
        if not self.is_finite():
            raise NotImplementedError("only implemented for finite groups")
        Irr = self.gap().Irr()
        L = [ClassFunction_libgap(self, irr) for irr in Irr]
        return tuple(L)

    def character(self, values):
        r"""
        Return a group character from ``values``, where ``values`` is
        a list of the values of the character evaluated on the conjugacy
        classes.

        INPUT:

        - ``values`` -- list of values of the character

        OUTPUT: a group character

        EXAMPLES::

            sage: G = MatrixGroup(AlternatingGroup(4))
            sage: G.character([1]*len(G.conjugacy_classes_representatives()))           # needs sage.rings.number_field
            Character of Matrix group over Integer Ring with 12 generators

        ::

            sage: G = GL(2,ZZ)
            sage: G.character([1,1,1,1])
            Traceback (most recent call last):
            ...
            NotImplementedError: only implemented for finite groups
        """
        if not self.is_finite():
            raise NotImplementedError("only implemented for finite groups")
        return ClassFunction_libgap(self, values)

    def trivial_character(self):
        r"""
        Return the trivial character of this group.

        OUTPUT: a group character

        EXAMPLES::

            sage: MatrixGroup(SymmetricGroup(3)).trivial_character()                    # needs sage.rings.number_field
            Character of Matrix group over Integer Ring with 6 generators

        ::

            sage: GL(2,ZZ).trivial_character()
            Traceback (most recent call last):
            ...
            NotImplementedError: only implemented for finite groups
        """
        if not self.is_finite():
            raise NotImplementedError("only implemented for finite groups")
        values = [1] * self._gap_().NrConjugacyClasses().sage()
        return self.character(values)

    def character_table(self):
        r"""
        Return the matrix of values of the irreducible characters of this
        group `G` at its conjugacy classes.

        The columns represent the conjugacy classes of
        `G` and the rows represent the different irreducible
        characters in the ordering given by GAP.

        OUTPUT: a matrix defined over a cyclotomic field

        EXAMPLES::

            sage: MatrixGroup(SymmetricGroup(2)).character_table()                      # needs sage.rings.number_field
            [ 1 -1]
            [ 1  1]
            sage: MatrixGroup(SymmetricGroup(3)).character_table()                      # needs sage.rings.number_field
            [ 1  1 -1]
            [ 2 -1  0]
            [ 1  1  1]
            sage: MatrixGroup(SymmetricGroup(5)).character_table()  # long time
            [ 1 -1 -1  1 -1  1  1]
            [ 4  0  1 -1 -2  1  0]
            [ 5  1 -1  0 -1 -1  1]
            [ 6  0  0  1  0  0 -2]
            [ 5 -1  1  0  1 -1  1]
            [ 4  0 -1 -1  2  1  0]
            [ 1  1  1  1  1  1  1]
        """
        # code from function in permgroup.py, but modified for
        # how gap handles these groups.
        G = self._gap_()
        cl = self.conjugacy_classes()
        from sage.rings.integer import Integer
        n = Integer(len(cl))
        irrG = G.Irr()
        ct = [[irrG[i][j] for j in range(n)] for i in range(n)]

        from sage.rings.number_field.number_field import CyclotomicField
        e = irrG.Flat().Conductor()
        K = CyclotomicField(e)
        ct = [[K(x) for x in v] for v in ct]

        # Finally return the result as a matrix.
        from sage.matrix.matrix_space import MatrixSpace
        MS = MatrixSpace(K, n)
        return MS(ct)

    def random_element(self):
        """
        Return a random element of this group.

        OUTPUT: a group element

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
            sage: gens = [MS([[1,2],[-1,1]]), MS([[1,1],[0,1]])]
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

            sage: from sage.groups.libgap_group import GroupLibGAP
            sage: G = GroupLibGAP(libgap.AlternatingGroup(5))
            sage: sum(1 for g in G)
            60
        """
        if self.list.cache is not None:
            yield from self.list()
            return
        iterator = self.gap().Iterator()
        while not iterator.IsDoneIterator().sage():
            yield self.element_class(self, iterator.NextIterator())

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
        size = self.gap().Size()
        if size.IsInfinity():
            raise NotImplementedError("group must be finite")
        return int(size)

    @cached_method
    def list(self):
        """
        List all elements of this group.

        OUTPUT: tuple containing all group elements in a random but fixed
        order

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

        An example over a ring (see :issue:`5241`)::

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

        An example over a field (see :issue:`10515`)::

            sage: gens = [matrix(QQ,2,[1,0,0,1])]
            sage: MatrixGroup(gens).list()
            (
            [1 0]
            [0 1]
            )

        Another example over a ring (see :issue:`9437`)::

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
        return tuple(self.element_class(self, g) for g in self.gap().AsList())

    def is_isomorphic(self, H):
        """
        Test whether ``self`` and ``H`` are isomorphic groups.

        INPUT:

        - ``H`` -- a group

        OUTPUT: boolean

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
            sage: F == G, G == H, F == H
            (False, False, False)
        """
        return self.gap().IsomorphismGroups(H.gap()) != libgap.fail
