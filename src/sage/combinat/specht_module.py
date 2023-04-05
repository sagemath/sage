r"""
Specht Modules

AUTHORS:

- Travis Scrimshaw (2023-1-22): initial version
"""

# ****************************************************************************
#       Copyright (C) 2023 Travis Scrimshaw <tcscrims (at) gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.combinat.diagram import Diagram
from sage.sets.family import Family
from sage.matrix.constructor import matrix
from sage.rings.rational_field import QQ
from sage.modules.with_basis.subquotient import SubmoduleWithBasis
from sage.categories.modules_with_basis import ModulesWithBasis

class SpechtModule(SubmoduleWithBasis):
    r"""
    A Specht module.

    Let `S_n` be the symmetric group on `n` letters and `R` be a commutative
    ring. The *Specht module* `S^D` for a diagram `D` is an `S_n`-module
    defined as follows. Let

    .. MATH::

        R(D) := \sum_{w \in R_D} w,
        \qquad\qquad
        C(D) := \sum_{w \in C_D} (-1)^w w,

    where `R_D` (resp. `C_D`) is the row (resp. column) stabilizer of `D`.
    Then, we construct the Specht module `S^D` as the left ideal

    .. MATH::

        S^D = R[S_n] C(D) R(D),

    where `R[S_n]` is the group algebra of `S_n` over `R`.

    INPUT:

    - ``SGA`` -- a symmetric group algebra
    - ``D`` -- a diagram

    EXAMPLES:

    We begin by constructing all irreducible Specht modules for the symmetric
    group `S_4` and show that they give a full set of irreducible
    representations both by having distinct Frobenius characters and the
    sum of the square of their dimensions is equal to `4!`::

        sage: SP = [la.specht_module(QQ) for la in Partitions(4)]
        sage: s = SymmetricFunctions(QQ).s()
        sage: [s(S.frobenius_image()) for S in SP]
        [s[4], s[3, 1], s[2, 2], s[2, 1, 1], s[1, 1, 1, 1]]
        sage: sum(S.dimension()^2 for S in SP)
        24

    Next, we compute the Specht module for a more general diagram
    for `S_5` and compute its irreducible decomposition by using
    its Frobenius character::

        sage: D = [(0,0), (0,1), (1,1), (1,2), (0,3)]
        sage: SGA = SymmetricGroupAlgebra(QQ, 5)
        sage: SM = SGA.specht_module(D)
        sage: SM.dimension()
        9
        sage: s(SM.frobenius_image())
        s[3, 2] + s[4, 1]

    This carries a natural (left) action of the symmetric group (algebra)::

        sage: S5 = SGA.group()
        sage: v = SM.an_element(); v
        2*B[0] + 2*B[1] + 3*B[2]
        sage: S5([2,1,5,3,4]) * v
        3*B[0] + 2*B[1] + 2*B[2]
        sage: x = SGA.an_element(); x
        [1, 2, 3, 4, 5] + 2*[1, 2, 3, 5, 4] + 3*[1, 2, 4, 3, 5] + [5, 1, 2, 3, 4]
        sage: x * v
        15*B[0] + 14*B[1] + 16*B[2] - 7*B[5] + 2*B[6] + 2*B[7]

    .. SEEALSO::

        :class:`~sage.combinat.symmetric_group_representations.SpechtRepresentation`
        for an implementation of the representation by matrices.
    """
    @staticmethod
    def __classcall_private__(cls, SGA, D):
        r"""
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: from sage.combinat.specht_module import SpechtModule
            sage: from sage.combinat.diagram import Diagram
            sage: SGA = SymmetricGroupAlgebra(QQ, 3)
            sage: D = [(0,0), (1,1), (1,2)]
            sage: SM1 = SpechtModule(SGA, D)
            sage: SM2 = SpechtModule(SGA, Diagram(D))
            sage: SM1 is SM2
            True
            sage: SM1 is SpechtModule(SGA, [[1,1], [1,2], [0,0]])
            True

            sage: SpechtModule(SGA, [[0,0], [1,1]])
            Traceback (most recent call last):
            ...
            ValueError: the domain size (=3) does not match the number of boxes (=2) of the diagram
        """
        D = _to_diagram(D)
        D = Diagram(D)
        n = len(D)
        if SGA.group().rank() != n - 1:
            rk = SGA.group().rank() + 1
            raise ValueError(f"the domain size (={rk}) does not match the number of boxes (={n}) of the diagram")
        return super().__classcall__(cls, SGA, D)

    def __init__(self, SGA, D):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 4)
            sage: SM = SGA.specht_module([(0,0), (1,1), (1,2), (2,1)])
            sage: TestSuite(SM).run()
        """
        self._diagram = D
        Mod = ModulesWithBasis(SGA.category().base_ring())
        span_set = specht_module_spanning_set(D, SGA)
        support_order = SGA.get_order()
        basis = SGA.echelon_form(span_set, False, order=support_order)
        basis = Family(basis)
        SubmoduleWithBasis.__init__(self, basis, support_order, ambient=SGA,
                                     unitriangular=False, category=Mod.Subobjects())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 4)
            sage: SGA.specht_module([(0,0), (1,1), (1,2), (2,1)])
            Specht module of [(0, 0), (1, 1), (1, 2), (2, 1)] over Rational Field
        """
        return f"Specht module of {self._diagram} over {self.base_ring()}"

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 4)
            sage: SM = SGA.specht_module([(0,0), (1,1), (1,2), (2,1)])
            sage: latex(SM)
                S^{{\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
                \raisebox{-.6ex}{$\begin{array}[b]{*{3}{p{0.6ex}}}\cline{1-1}
                \lr{\phantom{x}}&&\\\cline{1-1}\cline{2-2}\cline{3-3}
                &\lr{\phantom{x}}&\lr{\phantom{x}}\\\cline{2-2}\cline{3-3}\cline{2-2}
                &\lr{\phantom{x}}&\\\cline{2-2}
                \end{array}$}
                }}
        """
        from sage.misc.latex import latex
        return f"S^{{{latex(self._diagram)}}}"

    def _ascii_art_(self):
        r"""
        Return an ascii art representation of ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 4)
            sage: SM = SGA.specht_module([(0,0), (1,1), (1,2), (2,1)])
            sage: ascii_art(SM)
             O . .
             . O O
             . O .
            S
        """
        from sage.typeset.ascii_art import ascii_art
        return ascii_art("S", baseline=0) + ascii_art(self._diagram, baseline=-1)

    def _unicode_art_(self):
        r"""
        Return a unicode art representation of ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 4)
            sage: SM = SGA.specht_module([(0,0), (1,1), (1,2), (2,1)])
            sage: unicode_art(SM)
             ┌─┬─┬─┐
             │X│ │ │
             ├─┼─┼─┤
             │ │X│X│
             ├─┼─┼─┤
             │ │X│ │
             └─┴─┴─┘
            S
        """
        from sage.typeset.unicode_art import unicode_art
        return unicode_art("S", baseline=0) + unicode_art(self._diagram, baseline=-1)

    def representation_matrix(self, elt):
        r"""
        Return the matrix corresponding to the left action of the symmetric
        group (algebra) element ``elt`` on ``self``.

        .. SEEALSO::

            :class:`~sage.combinat.symmetric_group_representations.SpechtRepresentation`

        EXAMPLES::

            sage: SM = Partition([3,1,1]).specht_module(QQ)
            sage: SM.representation_matrix(Permutation([2,1,3,5,4]))
            [-1  0  0  1 -1  0]
            [ 0  0  1  0 -1  1]
            [ 0  1  0 -1  0  1]
            [ 0  0  0  0 -1  0]
            [ 0  0  0 -1  0  0]
            [ 0  0  0  0  0 -1]
            sage: SGA = SymmetricGroupAlgebra(QQ, 5)
            sage: SM.representation_matrix(SGA([3,1,5,2,4]))
            [ 0 -1  0  1  0 -1]
            [ 0  0  0  0  0 -1]
            [ 0  0  0 -1  0  0]
            [ 0  0 -1  0  1 -1]
            [ 1  0  0 -1  1  0]
            [ 0  0  0  0  1  0]
        """
        SGA = self._ambient
        return matrix(self.base_ring(), [self.retract(SGA(elt) * b.lift()).to_vector()
                                         for b in self.basis()])

    @cached_method
    def frobenius_image(self):
        r"""
        Return the Frobenius image of ``self``.

        The Frobenius map is defined as the map to symmetric functions

        .. MATH::

            F(\chi) = \frac{1}{n!} \sum_{w \in S_n} \chi(w) p_{\rho(w)},

        where `\chi` is the character of the `S_n`-module ``self``,
        `p_{\lambda}` is the powersum symmetric function basis element
        indexed by `\lambda`, and `\rho(w)` is partition of the cycle type
        of `w`. Specifically, this map takes irreducible representations
        indexed by `\lambda` to the Schur function `s_{\lambda}`.

        EXAMPLES::

            sage: s = SymmetricFunctions(QQ).s()
            sage: SM = Partition([2,2,1]).specht_module(QQ)
            sage: s(SM.frobenius_image())
            s[2, 2, 1]
            sage: SM = Partition([4,1]).specht_module(CyclotomicField(5))
            sage: s(SM.frobenius_image())
            s[4, 1]

        We verify the regular representation::

            sage: from sage.combinat.diagram import Diagram
            sage: D = Diagram([(0,0), (1,1), (2,2), (3,3), (4,4)])
            sage: F = s(D.specht_module(QQ).frobenius_image()); F
            s[1, 1, 1, 1, 1] + 4*s[2, 1, 1, 1] + 5*s[2, 2, 1]
             + 6*s[3, 1, 1] + 5*s[3, 2] + 4*s[4, 1] + s[5]
            sage: F == sum(StandardTableaux(la).cardinality() * s[la]
            ....:          for la in Partitions(5))
            True
            sage: all(s[la] == s(la.specht_module(QQ).frobenius_image())
            ....:     for n in range(1, 5) for la in Partitions(n))
            True

            sage: D = Diagram([(0,0), (1,1), (1,2), (2,3), (2,4)])
            sage: SM = D.specht_module(QQ)
            sage: s(SM.frobenius_image())
            s[2, 2, 1] + s[3, 1, 1] + 2*s[3, 2] + 2*s[4, 1] + s[5]
        """
        from sage.combinat.sf.sf import SymmetricFunctions
        BR = self._ambient.base_ring()
        p = SymmetricFunctions(BR).p()
        G = self._ambient.group()
        CCR = [(elt, elt.cycle_type()) for elt in G.conjugacy_classes_representatives()]
        return p.sum(self.representation_matrix(elt).trace() / la.centralizer_size() * p[la]
                     for elt, la in CCR)

    class Element(SubmoduleWithBasis.Element):
        def _acted_upon_(self, x, self_on_left=False):
            """
            Return the action of ``x`` on ``self``.

            INPUT:

            - ``x`` -- an element of the base ring or can be converted into
              the defining symmetric group algebra
            - ``self_on_left`` -- boolean (default: ``False``); which side
              ``self`` is on for the action

            EXAMPLES::

                sage: SGA = SymmetricGroupAlgebra(QQ, 4)
                sage: SM = SGA.specht_module([3,1])
                sage: SGA.an_element() * SM.an_element()
                14*B[0] + 18*B[1] + 8*B[2]
                sage: 4 * SM.an_element()
                8*B[0] + 8*B[1] + 12*B[2]
            """
            # Check for a scalar first
            ret = super()._acted_upon_(x, self_on_left)
            if ret is not None:
                return ret
            # Check if it is in the symmetric group algebra
            P = self.parent()
            if x in P._ambient or x in P._ambient.group():
                if self_on_left:  # it is only a left module
                    return None
                else:
                    return P.retract(P._ambient(x) * self.lift())
            return None

def _to_diagram(D):
    r"""
    Convert ``D`` to a list of cells representing a diagram.

    TESTS::

        sage: from sage.combinat.specht_module import _to_diagram
        sage: _to_diagram(Partition([3,1,1]))
        [(0, 0), (0, 1), (0, 2), (1, 0), (2, 0)]
        sage: _to_diagram(SkewPartition([[5,3,1,1],[2,2,1]]))
        [(0, 2), (0, 3), (0, 4), (1, 2), (3, 0)]
        sage: _to_diagram([1,2,0,2])
        [(0, 0), (1, 0), (1, 1), (3, 0), (3, 1)]
        sage: _to_diagram(Composition([2,1,3]))
        [(0, 0), (0, 1), (1, 0), (2, 0), (2, 1), (2, 2)]
        sage: _to_diagram([(1,2), (2,2)])
        [(1, 2), (2, 2)]
    """
    from sage.combinat.integer_vector import IntegerVectors
    from sage.combinat.skew_partition import SkewPartitions
    from sage.combinat.partition import _Partitions
    if isinstance(D, Diagram):
        return D
    if D in _Partitions:
        D = _Partitions(D).cells()
    elif D in SkewPartitions():
        D = SkewPartitions()(D).cells()
    elif D in IntegerVectors():
        cells = []
        for i, row in enumerate(D):
            for j in range(row):
                cells.append((i, j))
        D = cells
    else:
        D = [tuple(cell) for cell in D]
    return D

def specht_module_spanning_set(D, SGA=None):
    r"""
    Return a spanning set of the Specht module of diagram ``D``.

    INPUT:

    - ``D`` -- a list of cells ``(r,c)`` for row ``r`` and column ``c``
    - ``SGA`` -- optional; a symmetric group algebra

    EXAMPLES::

        sage: from sage.combinat.specht_module import specht_module_spanning_set
        sage: specht_module_spanning_set([(0,0), (1,1), (2,2)])
        ([1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1])
        sage: specht_module_spanning_set([(0,0), (1,1), (2,1)])
        ([1, 2, 3] - [1, 3, 2], -[1, 2, 3] + [1, 3, 2], [2, 1, 3] - [3, 1, 2],
         [2, 3, 1] - [3, 2, 1], -[2, 1, 3] + [3, 1, 2], -[2, 3, 1] + [3, 2, 1])

        sage: SGA = SymmetricGroup(3).algebra(QQ)
        sage: specht_module_spanning_set([(0,0), (1,1), (2,1)], SGA)
        (() - (2,3), -(1,2) + (1,3,2), (1,2,3) - (1,3),
         -() + (2,3), -(1,2,3) + (1,3), (1,2) - (1,3,2))
    """
    n = len(D)
    if SGA is None:
        from sage.combinat.symmetric_group_algebra import SymmetricGroupAlgebra
        SGA = SymmetricGroupAlgebra(QQ, n)
    elif SGA.group().rank() != n - 1:
        raise ValueError("the rank does not match the size of the diagram")
    row_diagram = [set() for _ in range(n)]
    col_diagram = [set() for _ in range(n)]
    for i,cell in enumerate(D):
        x,y = cell
        row_diagram[x].add(i)
        col_diagram[y].add(i)
    # Construct the row and column stabilizer elements
    row_stab = SGA.zero()
    col_stab = SGA.zero()
    B = SGA.basis()
    for w in B.keys():
            # Remember that the permutation w is 1-based
            row_perm = [set() for _ in range(n)]
            col_perm = [set() for _ in range(n)]
            for i,cell in enumerate(D):
                    x,y = cell
                    row_perm[x].add(w(i+1)-1)
                    col_perm[y].add(w(i+1)-1)
            if row_diagram == row_perm:
                    row_stab += B[w]
            if col_diagram == col_perm:
                    col_stab += w.sign() * B[w]
    gen = col_stab * row_stab
    return tuple([b * gen for b in B])

def specht_module_rank(D, base_ring=None):
    r"""
    Return the rank of the Specht module of diagram ``D``.

    EXAMPLES::

        sage: from sage.combinat.specht_module import specht_module_rank
        sage: specht_module_rank([(0,0), (1,1), (2,2)])
        6
    """
    D = _to_diagram(D)
    span_set = specht_module_spanning_set(D)
    if base_ring is None:
        base_ring = QQ
    return matrix(base_ring, [v.to_vector() for v in span_set]).rank()
