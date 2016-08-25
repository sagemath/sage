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
"""
from __future__ import absolute_import

##############################################################################
#       Copyright (C) 2006 David Joyner and William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
##############################################################################

import types
from sage.rings.all import ZZ
from sage.rings.integer import is_Integer
from sage.rings.ring import is_Ring
from sage.rings.finite_rings.finite_field_constructor import is_FiniteField
from sage.interfaces.gap import gap
from sage.matrix.matrix import is_Matrix
from sage.matrix.matrix_space import MatrixSpace, is_MatrixSpace
from sage.misc.latex import latex
from sage.structure.sequence import Sequence
from sage.structure.sage_object import SageObject
from sage.misc.decorators import rename_keyword
from sage.misc.cachefunc import cached_method
from sage.groups.generic import structure_description

from sage.groups.group import Group
from sage.groups.libgap_wrapper import ParentLibGAP
from sage.groups.libgap_mixin import GroupMixinLibGAP

from sage.groups.matrix_gps.group_element import (
    is_MatrixGroupElement, MatrixGroupElement_generic, MatrixGroupElement_gap)

#################################################################

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
    """

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

            sage: G = SU(2,GF(5))
            sage: G._check_matrix(identity_matrix(GF(5),2))
            sage: G._check_matrix(matrix(GF(5),[[1,1],[0,1]]))
            Traceback (most recent call last):
            ...
            TypeError: matrix must be unitary
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
        """
        if self.ngens() > 5:
            return 'Matrix group over {0} with {1} generators'.format(
                self.base_ring(), self.ngens())
        else:
            from sage.repl.display.util import format_list
            return 'Matrix group over {0} with {1} generators {2}'.format(
                self.base_ring(), self.ngens(), format_list(self.gens()))

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
        return '\\left\\langle %s \\right\\rangle'%gens



###################################################################
#
# Matrix group over a generic ring
#
###################################################################

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
        assert is_Ring(base_ring)
        assert is_Integer(degree)

        self._deg = degree
        if self._deg <= 0:
            raise ValueError('the degree must be at least 1')

        if (category is None) and is_FiniteField(base_ring):
            from sage.categories.finite_groups import FiniteGroups
            category = FiniteGroups()
        super(MatrixGroup_generic, self).__init__(base=base_ring, category=category)

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

    def _cmp_generators(self, other):
        """
        Implement comparison

        We treat two matrix groups as equal if their generators are
        the same in the same order. Infinitely-generated groups are
        compared by identity.

        INPUT:

        - ``other`` -- anything.

        OUTPUT:

        ``-1``, ``0``, or ``+1``.

        EXAMPLES::

            sage: G = GL(2,3)
            sage: H = MatrixGroup(G.gens())
            sage: G._cmp_generators(H)
            0
            sage: cmp(G,H)
            0
            sage: cmp(H,G)
            0
            sage: G == H
            True
        """
        if not is_MatrixGroup(other):
            return cmp(type(self), type(other))
        c = cmp(self.matrix_space(), other.matrix_space())
        if c != 0:
            return c

        def identity_cmp():
            return cmp(id(self), id(other))

        # compare number of generators
        try:
            n_self = self.ngens()
            n_other = other.ngens()
        except (AttributeError, NotImplementedError):
            return identity_cmp()
        c = cmp(n_self, n_other)
        if c != 0:
            return c
        from sage.structure.element import is_InfinityElement
        if is_InfinityElement(n_self) or is_InfinityElement(n_other):
            return identity_cmp()

        # compacte generator matrices
        try:
            self_gens = self.gens()
            other_gens = other.gens()
        except (AttributeError, NotImplementedError):
            return identity_cmp()
        assert(n_self == n_other)
        for g,h in zip(self_gens, other_gens):
            c = cmp(g.matrix(), h.matrix())
            if c != 0:
                return c
        return c

    def __cmp__(self, other):
        """
        Implement comparison

        EXAMPLES::

            sage: MS = MatrixSpace(SR, 2, 2)
            sage: G = MatrixGroup([MS(1), MS([1,2,3,4])])
            sage: from sage.groups.matrix_gps.matrix_group import MatrixGroup_generic
            sage: MatrixGroup_generic.__cmp__(G, G)
            0
            sage: cmp(G,G)
            0
            sage: cmp(G, MatrixGroup(G.gens()))
            0
            sage: G == MatrixGroup(G.gens())
            True
        """
        return self._cmp_generators(other)

    def _Hom_(self, G, cat=None):
        """
        Construct a homset.

        INPUT:

        - ``G`` -- group; the codomain

        - ``cat`` -- category; must be unset

        OUTPUT:

        The set of homomorphisms from ``self`` to ``G``.

        EXAMPLES::

            sage: MS = MatrixSpace(SR, 2, 2)
            sage: G = MatrixGroup([MS(1), MS([1,2,3,4])])
            sage: G.Hom(G)
            Set of Homomorphisms from Matrix group over Symbolic Ring with 2 generators (
            [1 0]  [1 2]
            [0 1], [3 4]
            ) to Matrix group over Symbolic Ring with 2 generators (
            [1 0]  [1 2]
            [0 1], [3 4]
            )

        TESTS:

        Check that :trac:`19407` is fixed::

            sage: G = GL(2, GF(2))
            sage: H = GL(3, ZZ)
            sage: Hom(G, H)
            Set of Homomorphisms from General Linear Group of degree 2
             over Finite Field of size 2 to General Linear Group of degree 3
             over Integer Ring
        """
        if not is_MatrixGroup(G):
            raise TypeError("G (=%s) must be a matrix group."%G)
        from . import homset
        return homset.MatrixGroupHomset(self, G, cat)

    def hom(self, x):
        """
        Return the group homomorphism defined by ``x``

        INPUT:

        - ``x`` -- a list/tuple/iterable of matrix group elements.

        OUTPUT:

        The group homomorphism defined by ``x``.

        EXAMPLES::

            sage: G = MatrixGroup([matrix(GF(5), [[1,3],[0,1]])])
            sage: H = MatrixGroup([matrix(GF(5), [[1,2],[0,1]])])
            sage: G.hom([H.gen(0)])
            Homomorphism : Matrix group over Finite Field of size 5 with 1 generators (
            [1 3]
            [0 1]
            ) --> Matrix group over Finite Field of size 5 with 1 generators (
            [1 2]
            [0 1]
            )
        """
        v = Sequence(x)
        U = v.universe()
        if not is_MatrixGroup(U):
            raise TypeError("u (=%s) must have universe a matrix group."%U)
        return self.Hom(U)(x)



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

        TESTS::

            sage: from sage.groups.matrix_gps.matrix_group import MatrixGroup_gap
            sage: MatrixGroup_gap(2, ZZ, libgap.eval('GL(2, Integers)'))
            Matrix group over Integer Ring with 3 generators (
            [0 1]  [-1  0]  [1 1]
            [1 0], [ 0  1], [0 1]
            )
        """
        ParentLibGAP.__init__(self, libgap_group, ambient=ambient)
        MatrixGroup_generic.__init__(self, degree, base_ring, category=category)

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
        libgap_contains = libgap.eval('\in')
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
            Matrix group over Integer Ring with 1 generators (
            [1 2]
            [0 1]
            )
            sage: G.ambient() is SL2Z
            True
        """
        from sage.groups.matrix_gps.finitely_generated import FinitelyGeneratedMatrixGroup_gap
        return FinitelyGeneratedMatrixGroup_gap(self.degree(), self.base_ring(),
                                                libgap_subgroup, ambient=self)

    def __iter__(self):
        """
        Iterate over the elements of the group.

        This method overrides the matrix group enumerator in GAP which
        is very slow, see http://tracker.gap-system.org/issues/369.

        EXAMPLES::

            sage: i = iter(GL(6,5))
            sage: [ next(i) for j in range(8) ]
            [
            [1 0 0 0 0 0]  [2 0 0 0 0 0]  [3 0 0 0 0 0]  [3 2 0 0 0 0]
            [0 1 0 0 0 0]  [0 1 0 0 0 0]  [0 1 0 0 0 0]  [0 1 0 0 0 0]
            [0 0 1 0 0 0]  [0 0 1 0 0 0]  [0 0 1 0 0 0]  [0 0 1 0 0 0]
            [0 0 0 1 0 0]  [0 0 0 1 0 0]  [0 0 0 1 0 0]  [0 0 0 1 0 0]
            [0 0 0 0 1 0]  [0 0 0 0 1 0]  [0 0 0 0 1 0]  [0 0 0 0 1 0]
            [0 0 0 0 0 1], [0 0 0 0 0 1], [0 0 0 0 0 1], [0 0 0 0 0 1],
            <BLANKLINE>
            [2 1 0 0 0 0]  [3 3 0 2 3 0]  [2 4 0 1 4 0]  [1 2 4 1 0 3]
            [0 1 0 0 0 0]  [0 1 0 0 0 0]  [0 1 0 0 0 0]  [0 1 0 0 0 0]
            [0 0 1 0 0 0]  [0 0 1 0 0 0]  [0 0 1 0 0 0]  [0 0 1 0 0 0]
            [0 0 0 1 0 0]  [0 0 0 1 0 0]  [0 0 0 1 0 0]  [0 0 0 1 0 0]
            [0 0 0 0 1 0]  [0 0 0 0 1 0]  [0 0 0 0 1 0]  [0 0 0 0 1 0]
            [0 0 0 0 0 1], [0 0 0 0 0 1], [0 0 0 0 0 1], [0 0 0 0 0 1]
            ]

        This is the direct computation in GAP, which will just run
        (almost) forever. If you find that this works then the
        ``MatrixGroup_gap.__iter__`` and ``MatrixGroup_gap.list``
        methods can be removed::

            sage: G = GL(6,5).gap()
            sage: G.Enumerator()   # not tested
            sage: G.Iterator()     # not tested
        """
        if not self.is_finite():
            # use implementation from category framework
            for g in super(Group, self).__iter__():
                yield g
            return
        # Use the standard GAP iterator for small groups
        if self.cardinality() < 1000:
            for g in super(MatrixGroup_gap, self).__iter__():
                yield g
            return
        # Override for large but finite groups
        iso = self.gap().IsomorphismPermGroup()
        P = iso.Image()
        iterator = P.Iterator()
        while not iterator.IsDoneIterator().sage():
            p = iterator.NextIterator()
            g = iso.PreImageElm(p)
            yield self(g, check=False)

    @cached_method
    def list(self):
        """
        List all elements of this group.

        This method overrides the matrix group enumerator in GAP which
        is very slow, see http://tracker.gap-system.org/issues/369.

        OUTPUT:

        A tuple containing all group elements in a random but fixed
        order.

        EXAMPLES::

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
        if not self.is_finite():
            raise NotImplementedError('group must be finite')
        return tuple(iter(self))

MatrixGroup_gap.structure_description = types.MethodType(structure_description, None, MatrixGroup_gap)
