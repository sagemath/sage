r"""
Groups of isometries

Let `M = \ZZ^n` or `\QQ^n`, `b: M \times M \rightarrow \QQ` a bilinear form and
`f: M \rightarrow M` a linear map. We say that `f` is an isometry if for all
elements `x,y` of `M` we have that `b(x,y)=b(f(x),f(y))`.
A group of isometries is a subgroup of `GL(M)` consisting of isometries.

EXAMPLES::

    sage: L = IntegralLattice("D4")                                                     # needs sage.graphs
    sage: O = L.orthogonal_group(); O                                                   # needs sage.graphs
    Group of isometries with 3 generators (
    [0 0 0 1]  [ 1  1  0  0]  [ 1  0  0  0]
    [0 1 0 0]  [ 0  0  1  0]  [-1 -1 -1 -1]
    [0 0 1 0]  [ 0  1  0  1]  [ 0  0  1  0]
    [1 0 0 0], [ 0 -1 -1  0], [ 0  0  0  1]
    )

Basic functionality is provided by GAP::

    sage: O.cardinality()                                                               # needs sage.graphs
    1152
    sage: len(O.conjugacy_classes_representatives())                                    # needs sage.graphs
    25

AUTHORS:

- Simon Brandhorst (2018-02): First created
"""

# ****************************************************************************
#       Copyright (C) 2018 Simon Brandhorst <sbrandhorst@web.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.groups.matrix_gps.finitely_generated_gap import FinitelyGeneratedMatrixGroup_gap
from sage.categories.action import Action


class GroupOfIsometries(FinitelyGeneratedMatrixGroup_gap):
    r"""
    A base class for Orthogonal matrix groups with a gap backend.

    Main difference to :class:`~sage.groups.matrix_gps.orthogonal.OrthogonalMatrixGroup_gap`
    is that we can specify generators and a bilinear form. Following GAP, the group action is
    from the right.

    INPUT:

    - ``degree`` -- integer; the degree (matrix size) of the matrix
    - ``base_ring`` -- ring, the base ring of the matrices
    - ``gens`` -- list of matrices over the base ring
    - ``invariant_bilinear_form`` -- a symmetric matrix
    - ``category`` -- (default: ``None``) a category of groups
    - ``check`` -- boolean (default: ``True``); check if the generators
      preserve the bilinear form
    - ``invariant_submodule`` -- a submodule preserved by the group action
      (default: ``None``); registers an action on this submodule
    - ``invariant_quotient_module`` -- a quotient module preserved by
      the group action (default: ``None``)
      registers an action on this quotient module

    EXAMPLES::

        sage: from sage.groups.matrix_gps.isometries import GroupOfIsometries
        sage: bil = Matrix(ZZ, 2, [3,2,2,3])
        sage: gens = [-Matrix(ZZ, 2, [0,1,1,0])]
        sage: O = GroupOfIsometries(2, ZZ, gens, bil)
        sage: O
        Group of isometries with 1 generator (
        [ 0 -1]
        [-1  0]
        )
        sage: O.order()
        2

    Infinite groups are O.K. too::

        sage: bil = Matrix(ZZ,4,[0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0])
        sage: f = Matrix(ZZ,4,[0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, -1, 1, 1, 1])
        sage: O = GroupOfIsometries(2, ZZ, [f], bil)
        sage: O.cardinality()
        +Infinity
    """

    def __init__(self, degree, base_ring,
                 gens, invariant_bilinear_form,
                 category=None, check=True,
                 invariant_submodule=None,
                 invariant_quotient_module=None):
        r"""
        Create this orthogonal group from the input.

        TESTS::

            sage: from sage.groups.matrix_gps.isometries import GroupOfIsometries
            sage: bil = Matrix(ZZ,2,[3,2,2,3])
            sage: gens = [-Matrix(ZZ,2,[0,1,1,0])]
            sage: cat = Groups().Finite()
            sage: O = GroupOfIsometries(2, ZZ, gens, bil, category=cat)
            sage: TestSuite(O).run()
        """
        from copy import copy
        G = copy(invariant_bilinear_form)
        G.set_immutable()
        self._invariant_bilinear_form = G
        self._invariant_submodule = invariant_submodule
        self._invariant_quotient_module = invariant_quotient_module
        if check:
            I = invariant_submodule
            Q = invariant_quotient_module
            for f in gens:
                self._check_matrix(f)
                if (I is not None) and I * f != I:
                    raise ValueError("the submodule is not preserved")
                if Q is not None and (Q.W() != Q.W()*f or Q.V()*f != Q.V()):
                    raise ValueError("the quotient module is not preserved")
        if len(gens) == 0:    # handle the trivial group
            gens = [G.parent().identity_matrix()]
        from sage.libs.gap.libgap import libgap
        gap_gens = [libgap(matrix_gen) for matrix_gen in gens]
        gap_group = libgap.Group(gap_gens)
        FinitelyGeneratedMatrixGroup_gap.__init__(self,
                                                  degree,
                                                  base_ring,
                                                  gap_group,
                                                  category=category)

    def _repr_(self):
        r"""
        Return the string representation of this matrix group.

        OUTPUT: string

        EXAMPLES::

            sage: from sage.groups.matrix_gps.isometries import GroupOfIsometries
            sage: bil = Matrix(ZZ,2,[3,2,2,3])
            sage: gens = [-Matrix(ZZ,2,[0,1,1,0])]
            sage: O = GroupOfIsometries(2,ZZ,gens,bil)
            sage: O
            Group of isometries with 1 generator (
            [ 0 -1]
            [-1  0]
            )
        """
        n = self.ngens()
        from sage.repl.display.util import format_list
        if n > 5:
            return 'Group of isometries with %s generators ' % n
        elif n == 1:
            return 'Group of isometries with %s generator %s' % (n, format_list(self.gens()))
        else:
            return 'Group of isometries with %s generators %s' % (n, format_list(self.gens()))

    def __reduce__(self):
        r"""
        Implement pickling.

        EXAMPLES::

            sage: from sage.groups.matrix_gps.isometries import GroupOfIsometries
            sage: bil = Matrix(ZZ,2,[3,2,2,3])
            sage: gens = [-Matrix(ZZ,2,[0,1,1,0])]
            sage: cat = Groups().Finite()
            sage: O = GroupOfIsometries(2, ZZ, gens, bil, category=cat)
            sage: loads(dumps(O)) == O
            True
        """
        args = (self.degree(), self.base_ring(),
                tuple(g.matrix() for g in self.gens()), self._invariant_bilinear_form,
                self.category(),
                False,
                self._invariant_submodule,
                self._invariant_quotient_module)
        return (GroupOfIsometries, args)

    def invariant_bilinear_form(self):
        r"""
        Return the symmetric bilinear form preserved by the orthogonal group.

        OUTPUT: the matrix defining the bilinear form

        EXAMPLES::

            sage: from sage.groups.matrix_gps.isometries import GroupOfIsometries
            sage: bil = Matrix(ZZ,2,[3,2,2,3])
            sage: gens = [-Matrix(ZZ,2,[0,1,1,0])]
            sage: O = GroupOfIsometries(2,ZZ,gens,bil)
            sage: O.invariant_bilinear_form()
            [3 2]
            [2 3]
        """
        return self._invariant_bilinear_form

    def _get_action_(self, S, op, self_on_left):
        """
        Provide the coercion system with an action.

        EXAMPLES::

            sage: from sage.groups.matrix_gps.isometries import GroupOfIsometries
            sage: bil = Matrix(ZZ,2,[3,2,2,3])
            sage: gens = [-Matrix(ZZ,2,[0,1,1,0])]
            sage: S = ZZ^2
            sage: T = S/(6*S)
            sage: O = GroupOfIsometries(2, ZZ, gens, bil, invariant_submodule=S, invariant_quotient_module=T)
            sage: O._get_action_(S, operator.mul, False)
            Right action by Group of isometries with 1 generator (
            [ 0 -1]
            [-1  0]
            ) on Ambient free module of rank 2 over the principal ideal domain Integer Ring
            sage: U = T.submodule([2*t for t in T.gens()])
            sage: u = U.an_element()
            sage: f = O.an_element()
            sage: u*f
            (0, 2)
        """
        import operator
        if op == operator.mul and not self_on_left:
            if S is self._invariant_submodule:
                return GroupActionOnSubmodule(self, S)
            if S is self._invariant_quotient_module:
                return GroupActionOnQuotientModule(self, S)
            from sage.modules.fg_pid.fgp_module import FGP_Module_class
            T = self._invariant_quotient_module
            if isinstance(S, FGP_Module_class):
                if S.is_submodule(T):
                    V = S.V()
                    if all(V == V * f.matrix() for f in self.gens()):
                        return GroupActionOnQuotientModule(self, S)
        return None

    def _check_matrix(self, x, *args):
        r"""
        Check whether the matrix ``x`` preserves the bilinear form.

        See :meth:`~sage.groups.matrix_gps.matrix_group._check_matrix`
        for details.

        EXAMPLES::

            sage: from sage.groups.matrix_gps.isometries import GroupOfIsometries
            sage: bil = Matrix(ZZ,2,[3,2,2,3])
            sage: gens = [-Matrix(ZZ,2,[0,1,1,0])]
            sage: O = GroupOfIsometries(2,ZZ,gens,bil)
            sage: g = matrix.identity(2)*2
            sage: O(g)
            Traceback (most recent call last):
            ...
            TypeError: matrix must be orthogonal with respect to the invariant form
        """
        F = self.invariant_bilinear_form()
        if x * F * x.transpose() != F:
            raise TypeError('matrix must be orthogonal '
                'with respect to the invariant form')


class GroupActionOnSubmodule(Action):
    r"""
    Matrix group action on a submodule from the right.

    INPUT:

    - ``MatrixGroup`` -- an instance of :class:`GroupOfIsometries`
    - ``submodule`` -- an invariant submodule
    - ``is_left`` -- boolean (default: ``False``)

    EXAMPLES::

        sage: from sage.groups.matrix_gps.isometries import GroupOfIsometries
        sage: S = span(ZZ, [[0,1]])
        sage: g = Matrix(QQ, 2, [1,0,0,-1])
        sage: G = GroupOfIsometries(2, ZZ, [g],
        ....:                       invariant_bilinear_form=matrix.identity(2),
        ....:                       invariant_submodule=S)
        sage: g = G.an_element()
        sage: x = S.an_element()
        sage: x*g
        (0, -1)
        sage: (x*g).parent()
        Free module of degree 2 and rank 1 over Integer Ring
        Echelon basis matrix:
        [0 1]
    """
    def __init__(self, MatrixGroup,submodule, is_left=False):
        r"""
        Initialize the action.

        TESTS::

            sage: from sage.groups.matrix_gps.isometries import GroupOfIsometries, GroupActionOnSubmodule
            sage: S = span(ZZ, [[0,1]])
            sage: g = Matrix(QQ, 2, [1,0,0,-1])
            sage: e = Matrix.identity(2)
            sage: G = GroupOfIsometries(2, ZZ, [g], e)
            sage: GroupActionOnSubmodule(G, S)
            Right action by Group of isometries with 1 generator (
            [ 1  0]
            [ 0 -1]
            ) on Free module of degree 2 and rank 1 over Integer Ring
            Echelon basis matrix:
            [0 1]
        """
        import operator
        Action.__init__(self, MatrixGroup, submodule, is_left, operator.mul)

    def _act_(self, g, a):
        r"""
        This defines the group action.

        INPUT:

        - ``g`` -- an element of the acting group

        - ``a`` -- an element of the invariant submodule

        OUTPUT: an element of the invariant submodule

        EXAMPLES::

            sage: from sage.groups.matrix_gps.isometries import GroupActionOnSubmodule
            sage: S = span(QQ,[[0,1]])
            sage: g = Matrix(QQ,2,[1,1,0,1/2])
            sage: G = MatrixGroup([g])
            sage: A = GroupActionOnSubmodule(G,S)
            sage: A
            Right action by Matrix group over Rational Field with 1 generators (
            [  1   1]
            [  0 1/2]
            ) on Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [0 1]
            sage: s = S.an_element()
            sage: g = G.an_element()
            sage: A(s,g).parent()
            Vector space of degree 2 and dimension 1 over Rational Field
            Basis matrix:
            [0 1]
        """
        if self.is_left():
            b = g.matrix() * a
        else:
            b = a * g.matrix()
        return a.parent()(b)


class GroupActionOnQuotientModule(Action):
    r"""
    Matrix group action on a quotient module from the right.

    INPUT:

    - ``MatrixGroup`` -- the group acting
      :class:`GroupOfIsometries`
    - ``submodule`` -- an invariant quotient module
    - ``is_left`` -- boolean (default: ``False``)

    EXAMPLES::

        sage: from sage.groups.matrix_gps.isometries import GroupOfIsometries
        sage: S = span(ZZ,[[0,1]])
        sage: Q = S/(6*S)
        sage: g = Matrix(QQ,2,[1,0,0,-1])
        sage: G = GroupOfIsometries(2, ZZ, [g], invariant_bilinear_form=matrix.identity(2), invariant_quotient_module=Q)
        sage: g = G.an_element()
        sage: x = Q.an_element()
        sage: x*g
        (5)
        sage: (x*g).parent()
        Finitely generated module V/W over Integer Ring with invariants (6)
    """
    def __init__(self, MatrixGroup, quotient_module, is_left=False):
        r"""
        Initialize the action.

        TESTS::

            sage: from sage.groups.matrix_gps.isometries import GroupOfIsometries
            sage: S = span(ZZ,[[0,1]])
            sage: Q = S/(6*S)
            sage: g = Matrix(QQ,2,[1,0,0,-1])
            sage: G = GroupOfIsometries(2, ZZ, [g], invariant_bilinear_form=matrix.identity(2), invariant_quotient_module=Q)
            sage: g = G.an_element()
            sage: x = Q.an_element()
            sage: x, x*g
            ((1), (5))
        """
        import operator
        Action.__init__(self, MatrixGroup, quotient_module, is_left, operator.mul)

    def _act_(self, g, a):
        r"""
        This defines the group action.

        INPUT:

        - ``g`` -- an element of the acting group

        - ``a`` -- an element of the invariant submodule

        OUTPUT: an element of the invariant quotient module

        EXAMPLES::

            sage: from sage.groups.matrix_gps.isometries import GroupActionOnQuotientModule
            sage: S = span(ZZ,[[0,1]])
            sage: Q = S/(6*S)
            sage: g = Matrix(QQ,2,[1,1,0,7])
            sage: G = MatrixGroup([g])
            sage: A = GroupActionOnQuotientModule(G,Q)
            sage: A
            Right action by Matrix group over Rational Field with 1 generators (
            [1 1]
            [0 7]
            ) on Finitely generated module V/W over Integer Ring with invariants (6)
            sage: q = Q.an_element()
            sage: g = G.an_element()
            sage: A(q,g).parent()
            Finitely generated module V/W over Integer Ring with invariants (6)
        """
        if self.is_left():
            b = g.matrix() * a.lift()
        else:
            b = a.lift() * g.matrix()
        return a.parent()(b)
