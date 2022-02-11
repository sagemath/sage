"""
Base classes for Matrix Groups

Loading, saving, ... works::

    sage: G = GL(2,5); G
    General Linear Group of degree 2 over Finite Field of size 5
    sage: TestSuite(G).run()

    sage: g = G.1; g
    [4 1]
    [4 0]
    sage: TestSuite(g).run()

We test that :trac:`9437` is fixed::

    sage: len(list(SL(2, Zmod(4))))
    48

AUTHORS:

- William Stein: initial version

- David Joyner (2006-03-15): degree, base_ring, _contains_, list,
  random, order methods; examples

- William Stein (2006-12): rewrite

- David Joyner (2007-12): Added invariant_generators (with Martin
  Albrecht and Simon King)

- David Joyner (2008-08): Added module_composition_factors (interface
  to GAP's MeatAxe implementation) and as_permutation_group (returns
  isomorphic PermutationGroup).

- Simon King (2010-05): Improve invariant_generators by using GAP
  for the construction of the Reynolds operator in Singular.

- Sebastian Oehms (2018-07): Add :meth:`subgroup` and :meth:`ambient` see :trac:`25894`
"""

# ****************************************************************************
#       Copyright (C) 2006 David Joyner and William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.groups import Groups
from sage.categories.rings import Rings
from sage.rings.integer import is_Integer
from sage.matrix.matrix_space import MatrixSpace
from sage.misc.latex import latex
from sage.structure.richcmp import (richcmp_not_equal, rich_to_bool,
                                    richcmp_method, richcmp)
from sage.misc.cachefunc import cached_method
from sage.groups.group import Group
from sage.groups.libgap_wrapper import ParentLibGAP
from sage.groups.libgap_mixin import GroupMixinLibGAP

from sage.groups.matrix_gps.group_element import (
    MatrixGroupElement_generic, MatrixGroupElement_gap)


def is_MatrixGroup(x):
    """
    Test whether ``x`` is a matrix group.

    EXAMPLES::

        sage: from sage.groups.matrix_gps.matrix_group import is_MatrixGroup
        sage: is_MatrixGroup(MatrixSpace(QQ,3))
        False
        sage: is_MatrixGroup(Mat(QQ,3))
        False
        sage: is_MatrixGroup(GL(2,ZZ))
        True
        sage: is_MatrixGroup(MatrixGroup([matrix(2,[1,1,0,1])]))
        True
    """
    return isinstance(x, MatrixGroup_base)

###################################################################
#
# Base class for all matrix groups
#
###################################################################


class MatrixGroup_base(Group):
    """
    Base class for all matrix groups.

    This base class just holds the base ring, but not the degree. So
    it can be a base for affine groups where the natural matrix is
    larger than the degree of the affine group. Makes no assumption
    about the group except that its elements have a ``matrix()``
    method.

    TESTS::

        sage: G = SO(3, GF(11)); G
        Special Orthogonal Group of degree 3 over Finite Field of size 11
        sage: G.category()
        Category of finite groups
    """
    _ambient = None  # internal attribute to register the ambient group in case this instance is a subgroup

    def _check_matrix(self, x, *args):
        """
        Check whether the matrix ``x`` defines a group element.

        This is used by the element constructor (if you pass
        ``check=True``, the default) that the defining matrix is valid
        for this parent. Derived classes must override this to verify
        that the matrix is, for example, orthogonal or symplectic.

        INPUT:

        - ``x`` -- a Sage matrix in the correct matrix space (degree
          and base ring).

        - ``*args`` -- optional other representations of ``x``,
          depending on the group implementation. Ignored by default.

        OUTPUT:

        A ``TypeError`` must be raised if ``x`` is invalid.

        EXAMPLES::

            sage: G = SU(2,GF(5)); F = G.base_ring() # this is GF(5^2,'a')
            sage: G._check_matrix(identity_matrix(F,2))
            sage: G._check_matrix(matrix(F,[[1,1],[0,1]]))
            Traceback (most recent call last):
            ...
            TypeError: matrix must be unitary with respect to the hermitian form
            [0 1]
            [1 0]
        """
        if not x.is_invertible():
            raise TypeError('matrix is not invertible')

    def as_matrix_group(self):
        """
        Return a new matrix group from the generators.

        This will throw away any extra structure (encoded in a derived
        class) that a group of special matrices has.

        EXAMPLES::

            sage: G = SU(4,GF(5))
            sage: G.as_matrix_group()
            Matrix group over Finite Field in a of size 5^2 with 2 generators (
            [      a       0       0       0]  [      1       0 4*a + 3       0]
            [      0 2*a + 3       0       0]  [      1       0       0       0]
            [      0       0 4*a + 1       0]  [      0 2*a + 4       0       1]
            [      0       0       0     3*a], [      0 3*a + 1       0       0]
            )

            sage: G = GO(3,GF(5))
            sage: G.as_matrix_group()
            Matrix group over Finite Field of size 5 with 2 generators (
            [2 0 0]  [0 1 0]
            [0 3 0]  [1 4 4]
            [0 0 1], [0 2 1]
            )
        """
        from sage.groups.matrix_gps.finitely_generated import MatrixGroup
        return MatrixGroup(self.gens())

    def subgroup(self, generators, check=True):
        """
        Return the subgroup generated by the given generators.

        INPUT:

        - ``generators`` -- a list/tuple/iterable of group elements of self
        - ``check`` -- boolean (optional, default: ``True``). Whether to check that each matrix is invertible.

        OUTPUT: The subgroup generated by ``generators`` as an instance of FinitelyGeneratedMatrixGroup_gap

        EXAMPLES::

            sage: UCF = UniversalCyclotomicField()
            sage: G  = GL(3, UCF)
            sage: e3 = UCF.gen(3); e5 =UCF.gen(5)
            sage: m = matrix(UCF, 3,3, [[e3, 1, 0], [0, e5, 7],[4, 3, 2]])
            sage: S = G.subgroup([m]); S
            Subgroup with 1 generators (
            [E(3)    1    0]
            [   0 E(5)    7]
            [   4    3    2]
            ) of General Linear Group of degree 3 over Universal Cyclotomic Field

            sage: CF3 = CyclotomicField(3)
            sage: G  = GL(3, CF3)
            sage: e3 = CF3.gen()
            sage: m = matrix(CF3, 3,3, [[e3, 1, 0], [0, ~e3, 7],[4, 3, 2]])
            sage: S = G.subgroup([m]); S
            Subgroup with 1 generators (
            [     zeta3          1          0]
            [         0 -zeta3 - 1          7]
            [         4          3          2]
            ) of General Linear Group of degree 3 over Cyclotomic Field of order 3 and degree 2

        TESTS::

            sage: TestSuite(G).run()
            sage: TestSuite(S).run()

            sage: W = CoxeterGroup(['I',7])
            sage: s = W.simple_reflections()
            sage: G = W.subgroup([s[1]])
            sage: G.category()
            Category of finite groups

            sage: W = WeylGroup(['A',2])
            sage: s = W.simple_reflections()
            sage: G = W.subgroup([s[1]])
            sage: G.category()
            Category of finite groups
        """
        try:
            test = self.is_finite()
        except NotImplementedError:
            test = self in Groups().Finite()
        cat = Groups().Finite() if test else Groups()
        # this method enlarges the method with same name of
        # ParentLibGAP to cases where the ambient group is not
        # inherited from ParentLibGAP.
        if isinstance(self, ParentLibGAP):
            return ParentLibGAP.subgroup(self, generators)

        for g in generators:
            if g not in self:
                raise ValueError("generator %s is not in the group" % (g))

        from sage.groups.matrix_gps.finitely_generated import MatrixGroup
        subgroup = MatrixGroup(generators, check=check, category=cat)
        subgroup._ambient = self
        return subgroup

    def ambient(self):
        """
        Return the ambient group of a subgroup.

        OUTPUT:

        A group containing ``self``. If ``self`` has not been defined
        as a subgroup, we just return ``self``.

        EXAMPLES::

            sage: G = GL(2,QQ)
            sage: m = matrix(QQ, 2,2, [[3, 0],[~5,1]])
            sage: S = G.subgroup([m])
            sage: S.ambient() is G
            True
        """
        if self._ambient is None:
            return self
        else:
            return self._ambient

    def _repr_(self):
        """
        Return a string representation.

        OUTPUT:

        String.

        EXAMPLES::

            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,2],[-1,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: G
            Matrix group over Finite Field of size 5 with 2 generators (
            [1 2]  [1 1]
            [4 1], [0 1]
            )

        case of being a subgroup::

            sage: CF3 = CyclotomicField(3)
            sage: G  = GL(2, CF3)
            sage: e3 = CF3.gen()
            sage: m = matrix(CF3, 2,2, [[e3, 1], [0, ~e3]])
            sage: S = G.subgroup([m]); S
            Subgroup with 1 generators (
            [     zeta3          1]
            [         0 -zeta3 - 1]
            ) of General Linear Group of degree 2 over Cyclotomic Field of order 3 and degree 2
        """
        ambient_group = self._ambient

        if ambient_group is None:
            if self.ngens() > 5:
                return 'Matrix group over {0} with {1} generators'.format(
                    self.base_ring(), self.ngens())
            else:
                from sage.repl.display.util import format_list
                return 'Matrix group over {0} with {1} generators {2}'.format(
                    self.base_ring(), self.ngens(), format_list(self.gens()))
        else:
            if self.ngens() > 5:
                return 'Subgroup with {0} generators of {1}'.format(
                    self.ngens(), ambient_group)
            else:
                from sage.repl.display.util import format_list
                return 'Subgroup with {0} generators {1} of {2}'.format(
                    self.ngens(), format_list(self.gens()), ambient_group)

    def _repr_option(self, key):
        """
        Metadata about the :meth:`_repr_` output.

        See :meth:`sage.structure.parent._repr_option` for details.

        EXAMPLES::

            sage: SO3 = groups.matrix.SO(3, QQ)
            sage: SO3._repr_option('element_ascii_art')
            True
        """
        if key == 'element_ascii_art':
            return True
        return super(MatrixGroup_base, self)._repr_option(key)

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: MS = MatrixSpace(GF(5), 2, 2)
            sage: G = MatrixGroup(MS([[1,2],[-1,1]]),MS([[1,1],[0,1]]))
            sage: latex(G)
            \left\langle \left(\begin{array}{rr}
            1 & 2 \\
            4 & 1
            \end{array}\right), \left(\begin{array}{rr}
            1 & 1 \\
            0 & 1
            \end{array}\right) \right\rangle
        """
        gens = ', '.join([latex(x) for x in self.gens()])
        return '\\left\\langle %s \\right\\rangle' % gens

    def sign_representation(self, base_ring=None, side="twosided"):
        r"""
        Return the sign representation of ``self`` over ``base_ring``.

        WARNING: assumes ``self`` is a matrix group over a field which has embedding over real numbers.

        INPUT:

        - ``base_ring`` -- (optional) the base ring; the default is `\ZZ`
        - ``side`` -- ignored

        EXAMPLES::

            sage: G = GL(2, QQ)
            sage: V = G.sign_representation()
            sage: e = G.an_element()
            sage: e
            [1 0]
            [0 1]
            sage: V._default_sign(e)
            1
            sage: m2 = V.an_element()
            sage: m2
            2*B['v']
            sage: m2*e
            2*B['v']
            sage: m2*e*e
            2*B['v']
        """
        if base_ring is None:
            from sage.rings.integer_ring import ZZ
            base_ring = ZZ
        from sage.modules.with_basis.representation import SignRepresentationMatrixGroup
        return SignRepresentationMatrixGroup(self, base_ring)

###################################################################
#
# Matrix group over a generic ring
#
###################################################################


@richcmp_method
class MatrixGroup_generic(MatrixGroup_base):

    Element = MatrixGroupElement_generic

    def __init__(self, degree, base_ring, category=None):
        """
        Base class for matrix groups over generic base rings

        You should not use this class directly. Instead, use one of
        the more specialized derived classes.

        INPUT:

        - ``degree`` -- integer. The degree (matrix size) of the
          matrix group.

        - ``base_ring`` -- ring. The base ring of the matrices.

        TESTS::

            sage: G = GL(2, QQ)
            sage: from sage.groups.matrix_gps.matrix_group import MatrixGroup_generic
            sage: isinstance(G, MatrixGroup_generic)
            True
        """
        assert base_ring in Rings
        assert is_Integer(degree)

        self._deg = degree
        if self._deg <= 0:
            raise ValueError('the degree must be at least 1')

        cat = Groups() if category is None else category
        if base_ring in Rings().Finite():
            cat = cat.Finite()
        super(MatrixGroup_generic, self).__init__(base=base_ring,
                                                  category=cat)

    def degree(self):
        """
        Return the degree of this matrix group.

        OUTPUT:

        Integer. The size (number of rows equals number of columns) of
        the matrices.

        EXAMPLES::

            sage: SU(5,5).degree()
            5
        """
        return self._deg

    @cached_method
    def matrix_space(self):
        """
        Return the matrix space corresponding to this matrix group.

        This is a matrix space over the field of definition of this matrix
        group.

        EXAMPLES::

            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: G = MatrixGroup([MS(1), MS([1,2,3,4])])
            sage: G.matrix_space()
            Full MatrixSpace of 2 by 2 dense matrices over Finite Field of size 5
            sage: G.matrix_space() is MS
            True
        """
        return MatrixSpace(self.base_ring(), self.degree())

    def __richcmp__(self, other, op):
        """
        Implement rich comparison.

        We treat two matrix groups as equal if their generators are
        the same in the same order. Infinitely-generated groups are
        compared by identity.

        INPUT:

        - ``other`` -- anything

        - ``op`` -- comparison operator

        OUTPUT:

        boolean

        EXAMPLES::

            sage: G = GL(2,3)
            sage: H = MatrixGroup(G.gens())
            sage: H == G
            True
            sage: G == H
            True

            sage: MS = MatrixSpace(QQ, 2, 2)
            sage: G = MatrixGroup([MS(1), MS([1,2,3,4])])
            sage: G == G
            True
            sage: G == MatrixGroup(G.gens())
            True

        TESTS::

            sage: G = groups.matrix.GL(4,2)
            sage: H = MatrixGroup(G.gens())
            sage: G == H
            True
            sage: G != H
            False
        """
        if not is_MatrixGroup(other):
            return NotImplemented

        if self is other:
            return rich_to_bool(op, 0)

        lx = self.matrix_space()
        rx = other.matrix_space()
        if lx != rx:
            return richcmp_not_equal(lx, rx, op)

        # compare number of generators
        try:
            n_self = self.ngens()
            n_other = other.ngens()
        except (AttributeError, NotImplementedError):
            return richcmp(id(self), id(other), op)

        if n_self != n_other:
            return richcmp_not_equal(self, other, op)

        from sage.structure.element import InfinityElement as Infinity
        if isinstance(n_self, Infinity) or isinstance(n_other, Infinity):
            return richcmp(id(self), id(other), op)

        # compact generator matrices
        try:
            self_gens = self.gens()
            other_gens = other.gens()
        except (AttributeError, NotImplementedError):
            return richcmp(id(self), id(other), op)

        assert(n_self == n_other)
        for g, h in zip(self_gens, other_gens):
            lx = g.matrix()
            rx = h.matrix()
            if lx != rx:
                return richcmp_not_equal(lx, rx, op)
        return rich_to_bool(op, 0)

###################################################################
#
# Matrix group over a ring that GAP understands
#
###################################################################


class MatrixGroup_gap(GroupMixinLibGAP, MatrixGroup_generic, ParentLibGAP):

    Element = MatrixGroupElement_gap

    def __init__(self, degree, base_ring, libgap_group, ambient=None, category=None):
        """
        Base class for matrix groups that implements GAP interface.

        INPUT:

        - ``degree`` -- integer. The degree (matrix size) of the
          matrix group.

        - ``base_ring`` -- ring. The base ring of the matrices.

        - ``libgap_group`` -- the defining libgap group.

        - ``ambient`` -- A derived class of :class:`ParentLibGAP` or
          ``None`` (default). The ambient class if ``libgap_group``
          has been defined as a subgroup.

        TESTS:

        ::

            sage: from sage.groups.matrix_gps.matrix_group import MatrixGroup_gap
            sage: MatrixGroup_gap(2, ZZ, libgap.eval('GL(2, Integers)'))
            Matrix group over Integer Ring with 3 generators (
            [0 1]  [-1  0]  [1 1]
            [1 0], [ 0  1], [0 1]
            )

        Check that the slowness of GAP iterators and enumerators for matrix groups
        (cf. http://tracker.gap-system.org/issues/369) has been fixed::

            sage: i = iter(GL(6,5))
            sage: [ next(i) for j in range(8) ]
            [
            [1 0 0 0 0 0]  [4 0 0 0 0 1]  [0 4 0 0 0 0]  [0 4 0 0 0 0]
            [0 1 0 0 0 0]  [4 0 0 0 0 0]  [0 0 4 0 0 0]  [0 0 4 0 0 0]
            [0 0 1 0 0 0]  [0 4 0 0 0 0]  [0 0 0 4 0 0]  [0 0 0 4 0 0]
            [0 0 0 1 0 0]  [0 0 4 0 0 0]  [0 0 0 0 4 0]  [0 0 0 0 4 0]
            [0 0 0 0 1 0]  [0 0 0 4 0 0]  [0 0 0 0 0 4]  [0 0 0 0 0 4]
            [0 0 0 0 0 1], [0 0 0 0 4 0], [1 4 0 0 0 0], [2 4 0 0 0 0],
            [3 0 0 0 0 1]  [4 0 0 1 3 3]  [0 0 0 2 0 0]  [1 0 0 0 4 4]
            [3 0 0 0 0 0]  [4 0 0 0 3 3]  [0 0 0 0 4 0]  [1 0 0 0 0 4]
            [0 4 0 0 0 0]  [3 0 0 0 0 1]  [2 2 0 0 0 2]  [1 0 0 0 0 0]
            [0 0 4 0 0 0]  [3 0 0 0 0 0]  [1 4 0 0 0 0]  [0 1 0 0 0 0]
            [0 0 0 4 0 0]  [0 4 0 0 0 0]  [0 2 4 0 0 0]  [0 0 1 0 0 0]
            [4 0 0 0 2 3], [2 0 3 4 4 4], [0 0 1 4 0 0], [0 0 0 1 0 0]
            ]

        And the same for listing the group elements, as well as few other issues::

            sage: F = GF(3)
            sage: gens = [matrix(F,2, [1,0, -1,1]), matrix(F, 2, [1,1,0,1])]
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
        ParentLibGAP.__init__(self, libgap_group, ambient=ambient)
        MatrixGroup_generic.__init__(self, degree, base_ring, category=category)

    def __iter__(self):
        """
        Iterate over the elements of the group.

        This method overrides the matrix group enumerator in GAP which
        does not (and often just cannot) work for infinite groups.

        TESTS:

        infinite groups can be dealt with::

            sage: import itertools
            sage: W = WeylGroup(["A",3,1])
            sage: list(itertools.islice(W, int(4)))
            [
            [1 0 0 0]  [-1  1  0  1]  [ 1  0  0  0]  [ 1  0  0  0]
            [0 1 0 0]  [ 0  1  0  0]  [ 1 -1  1  0]  [ 0  1  0  0]
            [0 0 1 0]  [ 0  0  1  0]  [ 0  0  1  0]  [ 0  1 -1  1]
            [0 0 0 1], [ 0  0  0  1], [ 0  0  0  1], [ 0  0  0  1]
            ]

        and finite groups, too::

            sage: G = GL(6,5)
            sage: list(itertools.islice(G, int(4)))
            [
            [1 0 0 0 0 0]  [4 0 0 0 0 1]  [0 4 0 0 0 0]  [0 4 0 0 0 0]
            [0 1 0 0 0 0]  [4 0 0 0 0 0]  [0 0 4 0 0 0]  [0 0 4 0 0 0]
            [0 0 1 0 0 0]  [0 4 0 0 0 0]  [0 0 0 4 0 0]  [0 0 0 4 0 0]
            [0 0 0 1 0 0]  [0 0 4 0 0 0]  [0 0 0 0 4 0]  [0 0 0 0 4 0]
            [0 0 0 0 1 0]  [0 0 0 4 0 0]  [0 0 0 0 0 4]  [0 0 0 0 0 4]
            [0 0 0 0 0 1], [0 0 0 0 4 0], [1 4 0 0 0 0], [2 4 0 0 0 0]
            ]
        """
        if not self.is_finite():
            # use implementation from category framework
            for g in super(Group, self).__iter__():
                yield g
            return
        # Use the standard GAP iterator for finite groups
        for g in super(MatrixGroup_gap, self).__iter__():
            yield g
        return

    def _check_matrix(self, x_sage, x_gap):
        """
        Check whether the matrix ``x`` defines a group element.

        This is used by the element constructor (if you pass
        ``check=True``, the default) that the defining matrix is valid
        for this parent. Derived classes must override this to verify
        that the matrix is, for example, orthogonal or symplectic.

        INPUT:

        - ``x_sage`` -- a Sage matrix in the correct matrix space (degree
          and base ring).

        - ``x_gap`` -- the corresponding LibGAP matrix.

        OUTPUT:

        A ``TypeError`` must be raised if ``x`` is invalid.

        EXAMPLES::

            sage: m1 = matrix(GF(11), [(0, -1), (1, 0)])
            sage: m2 = matrix(GF(11), [(0, -1), (1, -1)])
            sage: G = MatrixGroup([m1, m2])
            sage: G([1,2,0,1])
            [1 2]
            [0 1]
            sage: G([1,1,1,0])
            Traceback (most recent call last):
            ...
            TypeError: matrix is not in the finitely generated group
        """
        from sage.libs.gap.libgap import libgap
        libgap_contains = libgap.eval(r'\in')
        is_contained = libgap_contains(x_gap, self.gap())
        if not is_contained.sage():
            raise TypeError('matrix is not in the finitely generated group')

    def _subgroup_constructor(self, libgap_subgroup):
        """
        Return a finitely generated subgroup.

        See
        :meth:`sage.groups.libgap_wrapper.ParentLibGAP._subgroup_constructor`
        for details.

        TESTS::

            sage: SL2Z = SL(2,ZZ)
            sage: S, T = SL2Z.gens()
            sage: G = SL2Z.subgroup([T^2]); G   # indirect doctest
            Subgroup with 1 generators (
            [1 2]
            [0 1]
            ) of Special Linear Group of degree 2 over Integer Ring
            sage: G.ambient() is SL2Z
            True
        """
        cat = Groups()
        if self in Groups().Finite():
            cat = cat.Finite()
        from sage.groups.matrix_gps.finitely_generated import FinitelyGeneratedMatrixGroup_gap
        return FinitelyGeneratedMatrixGroup_gap(self.degree(), self.base_ring(),
                                                libgap_subgroup, ambient=self,
                                                category=cat)

    from sage.groups.generic import structure_description
