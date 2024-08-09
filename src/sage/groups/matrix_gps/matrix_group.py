"""
Base classes for Matrix Groups

TESTS:

Loading, saving, ... works::

    sage: # needs sage.libs.gap
    sage: G = GL(2,5); G
    General Linear Group of degree 2 over Finite Field of size 5
    sage: TestSuite(G).run()
    sage: g = G.1; g
    [4 1]
    [4 0]
    sage: TestSuite(g).run()

We test that :issue:`9437` is fixed::

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

- Sebastian Oehms (2018-07): Add :meth:`subgroup` and :meth:`ambient` see :issue:`25894`
"""

# ****************************************************************************
#       Copyright (C) 2006 David Joyner and William Stein <wstein@gmail.com>
#                     2009      Mike Hansen
#                     2013      Volker Braun <vbraun.name@gmail.com>
#                     2017-2021 Frédéric Chapoton
#                     2018-2019 Sebastian Oehms
#                     2020      Siddharth Singh
#                     2023      Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.groups import Groups
from sage.categories.rings import Rings
from sage.rings.integer import Integer
from sage.matrix.matrix_space import MatrixSpace
from sage.misc.latex import latex
from sage.structure.richcmp import (richcmp_not_equal, rich_to_bool,
                                    richcmp_method, richcmp)
from sage.misc.cachefunc import cached_method
from sage.groups.group import Group

from sage.groups.matrix_gps.group_element import MatrixGroupElement_generic


def is_MatrixGroup(x):
    """
    Test whether ``x`` is a matrix group.

    EXAMPLES::

        sage: from sage.groups.matrix_gps.matrix_group import is_MatrixGroup
        sage: is_MatrixGroup(MatrixSpace(QQ, 3))
        doctest:warning...
        DeprecationWarning: the function is_MatrixGroup is deprecated;
        use 'isinstance(..., MatrixGroup_base)' instead
        See https://github.com/sagemath/sage/issues/37898 for details.
        False
        sage: is_MatrixGroup(Mat(QQ, 3))
        False
        sage: is_MatrixGroup(GL(2, ZZ))
        True
        sage: is_MatrixGroup(MatrixGroup([matrix(2, [1,1,0,1])]))
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(37898, "the function is_MatrixGroup is deprecated; use 'isinstance(..., MatrixGroup_base)' instead")
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
          and base ring)

        - ``*args`` -- (optional) other representations of ``x``,
          depending on the group implementation. Ignored by default

        OUTPUT: a :exc:`TypeError` must be raised if ``x`` is invalid

        EXAMPLES::

            sage: G = SU(2, GF(5)); F = G.base_ring()  # this is GF(5^2,'a')            # needs sage.rings.finite_rings
            sage: G._check_matrix(identity_matrix(F, 2))                                # needs sage.rings.finite_rings
            sage: G._check_matrix(matrix(F, [[1,1], [0,1]]))                            # needs sage.rings.finite_rings
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

            sage: G = SU(4, GF(5))                                                      # needs sage.rings.finite_rings
            sage: G.as_matrix_group()                                                   # needs sage.libs.gap sage.rings.finite_rings
            Matrix group over Finite Field in a of size 5^2 with 2 generators (
            [      a       0       0       0]  [      1       0 4*a + 3       0]
            [      0 2*a + 3       0       0]  [      1       0       0       0]
            [      0       0 4*a + 1       0]  [      0 2*a + 4       0       1]
            [      0       0       0     3*a], [      0 3*a + 1       0       0]
            )

            sage: # needs sage.libs.gap
            sage: G = GO(3, GF(5))
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

        - ``generators`` -- list/tuple/iterable of group elements of ``self``
        - ``check`` -- boolean (default: ``True``); whether to check that each
          matrix is invertible

        OUTPUT: the subgroup generated by ``generators`` as an instance of
        :class:`FinitelyGeneratedMatrixGroup_gap`

        EXAMPLES::

            sage: # needs sage.libs.gap sage.rings.number_field
            sage: UCF = UniversalCyclotomicField()
            sage: G  = GL(3, UCF)
            sage: e3 = UCF.gen(3); e5 = UCF.gen(5)
            sage: m = matrix(UCF, 3,3, [[e3, 1, 0], [0, e5, 7],[4, 3, 2]])
            sage: S = G.subgroup([m]); S
            Subgroup with 1 generators (
            [E(3)    1    0]
            [   0 E(5)    7]
            [   4    3    2]
            ) of General Linear Group of degree 3 over Universal Cyclotomic Field

            sage: # needs sage.rings.number_field
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

            sage: TestSuite(G).run()                                                    # needs sage.rings.number_field
            sage: TestSuite(S).run()                                                    # needs sage.rings.number_field
        """
        try:
            test = self.is_finite()
        except NotImplementedError:
            test = self in Groups().Finite()
        cat = Groups().Finite() if test else Groups()

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

            sage: G = GL(2, QQ)
            sage: m = matrix(QQ, 2, 2, [[3, 0], [~5,1]])
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

        EXAMPLES::

            sage: F = GF(5); MS = MatrixSpace(F, 2, 2)
            sage: gens = [MS([[1,2], [-1,1]]), MS([[1,1], [0,1]])]
            sage: G = MatrixGroup(gens)
            sage: G
            Matrix group over Finite Field of size 5 with 2 generators (
            [1 2]  [1 1]
            [4 1], [0 1]
            )

        case of being a subgroup::

            sage: # needs sage.rings.number_field
            sage: CF3 = CyclotomicField(3)
            sage: G  = GL(2, CF3)
            sage: e3 = CF3.gen()
            sage: m = matrix(CF3, 2, 2, [[e3, 1], [0, ~e3]])
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
        return super()._repr_option(key)

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: MS = MatrixSpace(GF(5), 2, 2)
            sage: G = MatrixGroup(MS([[1,2], [-1,1]]), MS([[1,1], [0,1]]))
            sage: latex(G)
            \left\langle \left(\begin{array}{rr}
            1 & 2 \\
            4 & 1
            \end{array}\right), \left(\begin{array}{rr}
            1 & 1 \\
            0 & 1
            \end{array}\right) \right\rangle
        """
        gens = ', '.join(latex(x) for x in self.gens())
        return '\\left\\langle %s \\right\\rangle' % gens

    def sign_representation(self, base_ring=None):
        r"""
        Return the sign representation of ``self`` over ``base_ring``.

        INPUT:

        - ``base_ring`` -- (optional) the base ring; the default is the base
          ring of ``self``

        EXAMPLES::

            sage: G = GL(2, QQ)
            sage: V = G.sign_representation()
            sage: e = G.an_element()
            sage: e
            [1 0]
            [0 1]
            sage: m2 = V.an_element()
            sage: m2
            2*B['v']
            sage: m2*e
            2*B['v']
            sage: m2*e*e
            2*B['v']

            sage: W = WeylGroup(["A", 1, 1])
            sage: W.sign_representation()
            Sign representation of
             Weyl Group of type ['A', 1, 1] (as a matrix group acting on the root space)
             over Rational Field

            sage: G = GL(4, 2)
            sage: G.sign_representation() == G.trivial_representation()
            True
        """
        if base_ring is None:
            base_ring = self.base_ring()
        if base_ring.characteristic() == 2:  # characteristic 2
            return self.trivial_representation()
        from sage.modules.with_basis.representation import SignRepresentationMatrixGroup
        return SignRepresentationMatrixGroup(self, base_ring)

    def natural_representation(self, base_ring=None):
        r"""
        Return the natural representation of ``self`` over ``base_ring``.

        INPUT:

        - ``base_ring`` -- (optional) the base ring; the default is the base
          ring of ``self``

        EXAMPLES::

            sage: G = groups.matrix.SL(6, 3)
            sage: V = G.natural_representation()
            sage: V
            Natural representation of Special Linear Group of degree 6
             over Finite Field of size 3
            sage: e = prod(G.gens())
            sage: e
            [2 0 0 0 0 1]
            [2 0 0 0 0 0]
            [0 2 0 0 0 0]
            [0 0 2 0 0 0]
            [0 0 0 2 0 0]
            [0 0 0 0 2 0]
            sage: v = V.an_element()
            sage: v
            2*e[0] + 2*e[1]
            sage: e * v
            e[0] + e[1] + e[2]
        """
        from sage.modules.with_basis.representation import NaturalMatrixRepresentation
        return NaturalMatrixRepresentation(self, base_ring)


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
        Base class for matrix groups over generic base rings.

        You should not use this class directly. Instead, use one of
        the more specialized derived classes.

        INPUT:

        - ``degree`` -- integer; the degree (matrix size) of the
          matrix group

        - ``base_ring`` -- ring; the base ring of the matrices

        TESTS::

            sage: G = GL(2, QQ)
            sage: from sage.groups.matrix_gps.matrix_group import MatrixGroup_generic
            sage: isinstance(G, MatrixGroup_generic)
            True
        """
        assert base_ring in Rings
        assert isinstance(degree, Integer)

        self._deg = degree
        if self._deg <= 0:
            raise ValueError('the degree must be at least 1')

        cat = Groups() if category is None else category
        if base_ring in Rings().Finite():
            cat = cat.Finite()
        super().__init__(base=base_ring, category=cat)

    def degree(self):
        """
        Return the degree of this matrix group.

        OUTPUT: integer; the size (number of rows equals number of columns)
        of the matrices

        EXAMPLES::

            sage: SU(5,5).degree()                                                      # needs sage.rings.finite_rings
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

            sage: F = GF(5); MS = MatrixSpace(F, 2, 2)
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

        OUTPUT: boolean

        EXAMPLES::

            sage: # needs sage.libs.gap
            sage: G = GL(2,3)
            sage: H = MatrixGroup(G.gens())
            sage: H == G
            True
            sage: G == H
            True

            sage: # needs sage.libs.gap
            sage: MS = MatrixSpace(QQ, 2, 2)
            sage: G = MatrixGroup([MS(1), MS([1,2,3,4])])
            sage: G == G
            True
            sage: G == MatrixGroup(G.gens())
            True

        TESTS::

            sage: # needs sage.groups sage.rings.finite_rings
            sage: G = groups.matrix.GL(4,2)
            sage: H = MatrixGroup(G.gens())
            sage: G == H
            True
            sage: G != H
            False
        """
        if not isinstance(other, MatrixGroup_base):
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

        assert n_self == n_other
        for g, h in zip(self_gens, other_gens):
            lx = g.matrix()
            rx = h.matrix()
            if lx != rx:
                return richcmp_not_equal(lx, rx, op)
        return rich_to_bool(op, 0)

    def is_trivial(self):
        r"""
        Return ``True`` if this group is the trivial group.

        A group is trivial, if it consists only of the identity
        element, that is, if all its generators are the identity.

        EXAMPLES::

            sage: MatrixGroup([identity_matrix(3)]).is_trivial()
            True
            sage: SL(2, ZZ).is_trivial()
            False
            sage: CoxeterGroup(['B',3], implementation='matrix').is_trivial()
            False

        TESTS::

            sage: CoxeterGroup(['A',0], implementation='matrix').is_trivial()
            True
            sage: MatrixGroup([matrix(SR, [[1,x], [0,1]])]).is_trivial()
            False
            sage: G = MatrixGroup([identity_matrix(3), identity_matrix(3)])
            sage: G.ngens()
            2
            sage: G.is_trivial()
            True
        """
        return all(g.is_one() for g in self.gens())
