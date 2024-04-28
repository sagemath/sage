# sage.doctest: needs sage.combinat sage.groups sage.modules
r"""
Specht Modules

AUTHORS:

- Travis Scrimshaw (2023-1-22): initial version
- Travis Scrimshaw (2023-11-23): added simple modules based on code
  from Sacha Goldman

.. TODO::

    Integrate this with the implementations in
    :mod:`sage.modules.with_basis.representation`.
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
from sage.misc.lazy_attribute import lazy_attribute
from sage.combinat.diagram import Diagram
from sage.combinat.partition import _Partitions
from sage.combinat.free_module import CombinatorialFreeModule
from sage.modules.with_basis.representation import Representation_abstract
from sage.sets.family import Family
from sage.matrix.constructor import matrix
from sage.rings.rational_field import QQ
from sage.modules.with_basis.subquotient import SubmoduleWithBasis, QuotientModuleWithBasis
from sage.modules.free_module_element import vector
from sage.categories.modules_with_basis import ModulesWithBasis

class SymmetricGroupRepresentation:
    """
    Mixin class for symmetric group (algebra) representations.
    """
    def __init__(self, SGA):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: SM = Partition([3,1,1]).specht_module(GF(3))
            sage: TestSuite(SM).run()
        """
        self._semigroup = SGA.group()
        self._semigroup_algebra = SGA

    def side(self):
        r"""
        Return the side of the action defining ``self``.

        EXAMPLES::

            sage: SM = Partition([3,1,1]).specht_module(GF(3))
            sage: SM.side()
            'left'
        """
        return "left"

    @cached_method
    def frobenius_image(self):
        r"""
        Return the Frobenius image of ``self``.

        The Frobenius map is defined as the map to symmetric functions

        .. MATH::

            F(\chi) = \frac{1}{n!} \sum_{w \in S_n} \chi(w) p_{\rho(w)},

        where `\chi` is the character of the `S_n`-module ``self``,
        `p_{\lambda}` is the powersum symmetric function basis element
        indexed by `\lambda`, and `\rho(w)` is the cycle type of `w` as a
        partition. Specifically, this map takes irreducible representations
        indexed by `\lambda` to the Schur function `s_{\lambda}`.

        EXAMPLES::

            sage: SM = Partition([2,2,1]).specht_module(QQ)
            sage: SM.frobenius_image()
            s[2, 2, 1]
            sage: SM = Partition([4,1]).specht_module(CyclotomicField(5))
            sage: SM.frobenius_image()
            s[4, 1]

        We verify the regular representation::

            sage: from sage.combinat.diagram import Diagram
            sage: D = Diagram([(0,0), (1,1), (2,2), (3,3), (4,4)])
            sage: F = D.specht_module(QQ).frobenius_image(); F
            s[1, 1, 1, 1, 1] + 4*s[2, 1, 1, 1] + 5*s[2, 2, 1]
             + 6*s[3, 1, 1] + 5*s[3, 2] + 4*s[4, 1] + s[5]
            sage: s = SymmetricFunctions(QQ).s()
            sage: F == sum(StandardTableaux(la).cardinality() * s[la]
            ....:          for la in Partitions(5))
            True
            sage: all(s[la] == la.specht_module(QQ).frobenius_image()
            ....:     for n in range(1, 5) for la in Partitions(n))
            True

            sage: D = Diagram([(0,0), (1,1), (1,2), (2,3), (2,4)])
            sage: SM = D.specht_module(QQ)
            sage: SM.frobenius_image()
            s[2, 2, 1] + s[3, 1, 1] + 2*s[3, 2] + 2*s[4, 1] + s[5]

        An example using the tabloid module::

            sage: SGA = SymmetricGroupAlgebra(QQ, 5)
            sage: TM = SGA.tabloid_module([2, 2, 1])
            sage: TM.frobenius_image()
            s[2, 2, 1] + s[3, 1, 1] + 2*s[3, 2] + 2*s[4, 1] + s[5]
        """
        from sage.combinat.sf.sf import SymmetricFunctions
        p = SymmetricFunctions(QQ).p()
        s = SymmetricFunctions(QQ).s()
        G = self._semigroup
        CCR = [(elt, elt.cycle_type()) for elt in G.conjugacy_classes_representatives()]
        return s(p.sum(QQ(self.representation_matrix(elt).trace()) / la.centralizer_size() * p[la]
                       for elt, la in CCR))

    # TODO: Move these methods up to methods of general representations

    def representation_matrix(self, elt):
        r"""
        Return the matrix corresponding to the left action of the symmetric
        group (algebra) element ``elt`` on ``self``.

        EXAMPLES::

            sage: SM = Partition([3,1,1]).specht_module(QQ)
            sage: SM.representation_matrix(Permutation([2,1,3,5,4]))
            [-1  0  0  0  0  0]
            [ 0  0  0 -1  0  0]
            [ 1  0  0 -1  1  0]
            [ 0 -1  0  0  0  0]
            [ 1 -1  1  0  0  0]
            [ 0 -1  0  1  0 -1]

            sage: SGA = SymmetricGroupAlgebra(QQ, 5)
            sage: SM = SGA.specht_module([(0,0), (0,1), (0,2), (1,0), (2,0)])
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
        return matrix(self.base_ring(), [(elt * b).to_vector() for b in self.basis()])

    @cached_method
    def character(self):
        r"""
        Return the character of ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 5)
            sage: SM = SGA.specht_module([3,2])
            sage: SM.character()
            (5, 1, 1, -1, 1, -1, 0)
            sage: matrix(SGA.specht_module(la).character() for la in Partitions(5))
            [ 1  1  1  1  1  1  1]
            [ 4  2  0  1 -1  0 -1]
            [ 5  1  1 -1  1 -1  0]
            [ 6  0 -2  0  0  0  1]
            [ 5 -1  1 -1 -1  1  0]
            [ 4 -2  0  1  1  0 -1]
            [ 1 -1  1  1 -1 -1  1]

            sage: SGA = SymmetricGroupAlgebra(QQ, SymmetricGroup(5))
            sage: SM = SGA.specht_module([3,2])
            sage: SM.character()
            Character of Symmetric group of order 5! as a permutation group
            sage: SM.character().values()
            [5, 1, 1, -1, 1, -1, 0]
            sage: matrix(SGA.specht_module(la).character().values() for la in reversed(Partitions(5)))
            [ 1 -1  1  1 -1 -1  1]
            [ 4 -2  0  1  1  0 -1]
            [ 5 -1  1 -1 -1  1  0]
            [ 6  0 -2  0  0  0  1]
            [ 5  1  1 -1  1 -1  0]
            [ 4  2  0  1 -1  0 -1]
            [ 1  1  1  1  1  1  1]
            sage: SGA.group().character_table()
            [ 1 -1  1  1 -1 -1  1]
            [ 4 -2  0  1  1  0 -1]
            [ 5 -1  1 -1 -1  1  0]
            [ 6  0 -2  0  0  0  1]
            [ 5  1  1 -1  1 -1  0]
            [ 4  2  0  1 -1  0 -1]
            [ 1  1  1  1  1  1  1]
        """
        G = self._semigroup
        chi = [self.representation_matrix(g).trace()
               for g in G.conjugacy_classes_representatives()]
        try:
            return G.character(chi)
        except AttributeError:
            return vector(chi, immutable=True)

    @cached_method
    def brauer_character(self):
        r"""
        Return the Brauer character of ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(GF(2), 5)
            sage: SM = SGA.specht_module([3,2])
            sage: SM.brauer_character()
            (5, -1, 0)
            sage: SM.simple_module().brauer_character()
            (4, -2, -1)
        """
        from sage.rings.number_field.number_field import CyclotomicField
        from sage.arith.functions import lcm
        G = self._semigroup
        p = self.base_ring().characteristic()
        # We manually compute the order since a Permutation does not implement order()
        chi = []
        for g in G.conjugacy_classes_representatives():
            if p.divides(lcm(g.cycle_type())):
                # ignore the non-p-regular elements
                continue
            evals = self.representation_matrix(g).eigenvalues()
            K = evals[0].parent()
            val = 0
            orders = {la: la.multiplicative_order() for la in evals if la != K.one()}
            zetas = {o: CyclotomicField(o).gen() for o in orders.values()}
            prims = {o: K.zeta(o) for o in orders.values()}
            for la in evals:
                if la == K.one():
                    val += 1
                    continue
                o = la.multiplicative_order()
                zeta = zetas[o]
                prim = prims[o]
                for deg in range(o):
                    if prim ** deg == la:
                        val += zeta ** deg
                        break
            chi.append(val)

        return vector(chi, immutable=True)


class SpechtModule(SubmoduleWithBasis, SymmetricGroupRepresentation, Representation_abstract):
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
        2*S[0] + 2*S[1] + 3*S[2]
        sage: S5([2,1,5,3,4]) * v
        3*S[0] + 2*S[1] + 2*S[2]
        sage: x = SGA.an_element(); x
        [1, 2, 3, 4, 5] + 2*[1, 2, 3, 5, 4] + 3*[1, 2, 4, 3, 5] + [5, 1, 2, 3, 4]
        sage: x * v
        15*S[0] + 14*S[1] + 16*S[2] - 7*S[5] + 2*S[6] + 2*S[7]

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
        if D in _Partitions:
            D = _Partitions(D)
            return TabloidModule(SGA, D).specht_module()
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
        SymmetricGroupRepresentation.__init__(self, SGA)
        self._diagram = D
        Mod = ModulesWithBasis(SGA.category().base_ring())
        span_set = specht_module_spanning_set(D, SGA)
        support_order = SGA.get_order()
        basis = SGA.echelon_form(span_set, False, order=support_order)
        basis = Family(basis)
        SubmoduleWithBasis.__init__(self, basis, support_order, ambient=SGA,
                                    unitriangular=False, category=Mod.Subobjects(),
                                    prefix='S')

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 4)
            sage: SGA.specht_module([(0,0), (1,1), (1,2), (2,1)])
            Specht module of [(0, 0), (1, 1), (1, 2), (2, 1)] over Rational Field

            sage: SGA.specht_module([3, 1])
            Specht module of [3, 1] over Rational Field
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

            sage: SM = SGA.specht_module([3,1])
            sage: latex(SM)
            S^{{\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{3}c}\cline{1-3}
            \lr{\phantom{x}}&\lr{\phantom{x}}&\lr{\phantom{x}}\\\cline{1-3}
            \lr{\phantom{x}}\\\cline{1-1}
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

            sage: SM = SGA.specht_module([3,1])
            sage: ascii_art(SM)
             ***
             *
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

            sage: SM = SGA.specht_module([3,1])
            sage: unicode_art(SM)
             ┌┬┬┐
             ├┼┴┘
             └┘
            S
        """
        from sage.typeset.unicode_art import unicode_art
        return unicode_art("S", baseline=0) + unicode_art(self._diagram, baseline=-1)

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

                sage: SGA = SymmetricGroupAlgebra(QQ, 5)
                sage: from sage.combinat.diagram import Diagram
                sage: D = Diagram([(0,0), (1,1), (2,2), (3,3), (4,4)])
                sage: SM = SGA.specht_module(D)
                sage: SGA.an_element() * SM.an_element()
                15*S[0] + 6*S[1] + 9*S[2] + 6*S[3] + 6*S[4] + 2*S[72] + 2*S[96] + 3*S[97]

                sage: SGA = SymmetricGroupAlgebra(QQ, 4)
                sage: SM = SGA.specht_module([3,1])
                sage: SGA.an_element() * SM.an_element()
                9*S[[1, 2, 3], [4]] + 17*S[[1, 2, 4], [3]] + 14*S[[1, 3, 4], [2]]
                sage: 4 * SM.an_element()
                12*S[[1, 2, 3], [4]] + 8*S[[1, 2, 4], [3]] + 8*S[[1, 3, 4], [2]]

            TESTS::

                sage: SGA = SymmetricGroupAlgebra(QQ, 4)
                sage: SM = SGA.specht_module([3,1])
                sage: SM.an_element() * SGA.an_element()
                Traceback (most recent call last):
                ...
                TypeError: unsupported operand parent(s) for *:
                 'Specht module of [3, 1] over Rational Field'
                 and 'Symmetric group algebra of order 4 over Rational Field'
                sage: groups.permutation.Dihedral(3).an_element() * SM.an_element()
                Traceback (most recent call last):
                ...
                TypeError: unsupported operand parent(s) for *:
                 'Dihedral group of order 6 as a permutation group'
                 and 'Specht module of [3, 1] over Rational Field'
            """
            # Check for a scalar first
            ret = super()._acted_upon_(x, self_on_left)
            if ret is not None:
                return ret
            # Check if it is in the symmetric group algebra
            P = self.parent()
            if x in P._semigroup_algebra or x in P._semigroup_algebra.group():
                if self_on_left:  # it is only a left module
                    return None
                else:
                    return P.retract(P._semigroup_algebra(x) * self.lift())
            return None


class TabloidModule(SymmetricGroupRepresentation, Representation_abstract):
    r"""
    The vector space of all tabloids of a fixed shape with the natural
    symmetric group action.

    A *tabloid* is an :class:`OrderedSetPartition` whose underlying set
    is `\{1, \ldots, n\}`. The symmetric group acts by permuting the
    entries of the set. Hence, this is a representation of the symmetric
    group defined over any field.

    EXAMPLES::

        sage: SGA = SymmetricGroupAlgebra(GF(3), 5)
        sage: TM = SGA.tabloid_module([2, 2, 1])
        sage: TM.dimension()
        30
        sage: TM.brauer_character()
        (30, 6, 2, 0, 0)
        sage: IM = TM.invariant_module()
        sage: IM.dimension()
        1
        sage: IM.basis()[0].lift() == sum(TM.basis())
        True
    """
    @staticmethod
    def __classcall_private__(cls, SGA, shape):
        r"""
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: from sage.combinat.specht_module import TabloidModule
            sage: SGA = SymmetricGroupAlgebra(QQ, 5)
            sage: TM1 = TabloidModule(SGA, [2, 2, 1])
            sage: TM2 = TabloidModule(SGA, Partition([2, 2, 1]))
            sage: TM1 is TM2
            True

            sage: TabloidModule(SGA, [3, 2, 1])
            Traceback (most recent call last):
            ...
            ValueError: the domain size (=5) does not match the number of boxes (=6) of the diagram
        """
        shape = _Partitions(shape)
        if SGA.group().rank() != sum(shape) - 1:
            rk = SGA.group().rank() + 1
            n = sum(shape)
            raise ValueError(f"the domain size (={rk}) does not match the number of boxes (={n}) of the diagram")
        return super().__classcall__(cls, SGA, shape)

    def __init__(self, SGA, shape):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 5)
            sage: TM = SGA.tabloid_module([2,2,1])
            sage: TestSuite(TM).run()
        """
        from sage.combinat.set_partition_ordered import OrderedSetPartitions
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        self._shape = shape
        n = sum(shape)
        self._symgp = SymmetricGroup(n)
        cat = ModulesWithBasis(SGA.base_ring()).FiniteDimensional()
        tabloids = OrderedSetPartitions(n, shape)
        CombinatorialFreeModule.__init__(self, SGA.base_ring(), tabloids,
                                         category=cat, prefix='T', bracket='')
        SymmetricGroupRepresentation.__init__(self, SGA)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 5)
            sage: SGA.tabloid_module([2,2,1])
            Tabloid module of [2, 2, 1] over Rational Field
        """
        return f"Tabloid module of {self._shape} over {self.base_ring()}"

    def _ascii_art_term(self, T):
        r"""
        Return an ascii art representation of the term indexed by ``T``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 5)
            sage: TM = SGA.tabloid_module([2,2,1])
            sage: ascii_art(TM.an_element())  # indirect doctest
            2*T       + 2*T       + 3*T
               {1, 2}      {1, 2}      {1, 2}
               {3, 4}      {3, 5}      {4, 5}
               {5}         {4}         {3}
        """
        # This is basically copied from CombinatorialFreeModule._ascii_art_term
        from sage.typeset.ascii_art import AsciiArt, ascii_art
        pref = AsciiArt([self.prefix()])
        tab = "\n".join("{" + ", ".join(str(val) for val in sorted(row)) + "}" for row in T)
        if not tab:
            tab = '-'
        r = pref * (AsciiArt([" " * len(pref)]) + ascii_art(tab))
        r._baseline = r._h - 1
        return r

    def _unicode_art_term(self, T):
        r"""
        Return a unicode art representation of the term indexed by ``T``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 5)
            sage: TM = SGA.tabloid_module([2,2,1])
            sage: unicode_art(TM.an_element())  # indirect doctest
            2*T       + 2*T       + 3*T
               {1, 2}      {1, 2}      {1, 2}
               {3, 4}      {3, 5}      {4, 5}
               {5}         {4}         {3}
        """
        from sage.typeset.unicode_art import unicode_art
        r = unicode_art(repr(self._ascii_art_term(T)))
        r._baseline = r._h - 1
        return r

    def _latex_term(self, T):
        r"""
        Return a latex representation of the term indexed by ``T``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 5)
            sage: TM = SGA.tabloid_module([2,2,1])
            sage: latex(TM.an_element())  # indirect doctest
            2 T_{{\def\lr#1{\multicolumn{1}{@{\hspace{.6ex}}c@{\hspace{.6ex}}}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{2}c}\cline{1-2}
            \lr{1}&\lr{2}\\\cline{1-2}
            \lr{3}&\lr{4}\\\cline{1-2}
            \lr{5}\\\cline{1-1}
            \end{array}$}
            }} + 2 T_{{\def\lr#1{\multicolumn{1}{@{\hspace{.6ex}}c@{\hspace{.6ex}}}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{2}c}\cline{1-2}
            \lr{1}&\lr{2}\\\cline{1-2}
            \lr{3}&\lr{5}\\\cline{1-2}
            \lr{4}\\\cline{1-1}
            \end{array}$}
            }} + 3 T_{{\def\lr#1{\multicolumn{1}{@{\hspace{.6ex}}c@{\hspace{.6ex}}}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{2}c}\cline{1-2}
            \lr{1}&\lr{2}\\\cline{1-2}
            \lr{4}&\lr{5}\\\cline{1-2}
            \lr{3}\\\cline{1-1}
            \end{array}$}
            }}
        """
        if not T:
            tab = "\\emptyset"
        else:
            from sage.combinat.output import tex_from_array
            A = list(map(sorted, T))
            tab = str(tex_from_array(A))
            tab = tab.replace("|", "")
        return f"{self.prefix()}_{{{tab}}}"

    def _symmetric_group_action(self, osp, g):
        r"""
        Return the action of the symmetric group element ``g`` on the
        ordered set partition ``osp``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 5)
            sage: TM = SGA.tabloid_module([2,2,1])
            sage: osp = TM._indices([[1,4],[3,5],[2]])
            sage: g = SGA.group().an_element(); g
            [5, 1, 2, 3, 4]
            sage: TM._symmetric_group_action(osp, g)
            [{3, 5}, {2, 4}, {1}]
        """
        P = self._indices
        g = self._symgp(g)
        return P.element_class(P, [[g(val) for val in row] for row in osp], check=False)

    def specht_module(self):
        r"""
        Return the Specht submodule of ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 5)
            sage: TM = SGA.tabloid_module([2,2,1])
            sage: TM.specht_module() is SGA.specht_module([2,2,1])
            True
        """
        return SpechtModuleTableauxBasis(self)

    def bilinear_form(self, u, v):
        r"""
        Return the natural bilinear form of ``self`` applied to ``u`` and ``v``.

        The natural bilinear form is given by defining the tabloid basis
        to be orthonormal.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 5)
            sage: TM = SGA.tabloid_module([2,2,1])
            sage: u = TM.an_element(); u
            2*T[{1, 2}, {3, 4}, {5}] + 2*T[{1, 2}, {3, 5}, {4}] + 3*T[{1, 2}, {4, 5}, {3}]
            sage: v = sum(TM.basis())
            sage: TM.bilinear_form(u, v)
            7
            sage: TM.bilinear_form(u, TM.zero())
            0
        """
        if len(v) < len(u):
            u, v = v, u
        R = self.base_ring()
        return R.sum(c * v[T] for T, c in u)

    class Element(CombinatorialFreeModule.Element):
        def _acted_upon_(self, x, self_on_left):
            r"""
            Return the action of ``x`` on ``self``.

            INPUT:

            - ``x`` -- an element of the base ring or can be converted into
              the defining symmetric group algebra
            - ``self_on_left`` -- boolean (default: ``False``); which side
              ``self`` is on for the action

            EXAMPLES::

                sage: SGA = SymmetricGroupAlgebra(QQ, 5)
                sage: SM = SGA.tabloid_module([2,2,1])
                sage: SGA.an_element() * SM.an_element()
                2*T[{1, 5}, {2, 3}, {4}] + 2*T[{1, 5}, {2, 4}, {3}] + 3*T[{1, 5}, {3, 4}, {2}]
                 + 12*T[{1, 2}, {3, 4}, {5}] + 15*T[{1, 2}, {3, 5}, {4}] + 15*T[{1, 2}, {4, 5}, {3}]
                sage: 4 * SM.an_element()
                8*T[{1, 2}, {3, 4}, {5}] + 8*T[{1, 2}, {3, 5}, {4}] + 12*T[{1, 2}, {4, 5}, {3}]
                sage: SM.an_element() * SGA.an_element()
                Traceback (most recent call last):
                ...
                TypeError: unsupported operand parent(s) for *:
                 'Tabloid module of [2, 2, 1] over Rational Field'
                 and 'Symmetric group algebra of order 5 over Rational Field'
            """
            # first check for the base action
            ret = super()._acted_upon_(x, self_on_left)
            if ret is not None:
                return ret

            if self_on_left:
                return None
            P = self.parent()
            if x in P._semigroup_algebra:
                return P.sum(c * (perm * self) for perm, c in x.monomial_coefficients().items())
            if x in P._semigroup_algebra.indices():
                return P.element_class(P, {P._symmetric_group_action(T, x): c
                                           for T, c in self._monomial_coefficients.items()})


class SpechtModuleTableauxBasis(SpechtModule):
    r"""
    A Specht module of a partition in the classical standard
    tableau basis.

    This is constructed as a `S_n`-submodule of the :class:`TabloidModule`
    (also referred to as the standard module).

    .. SEEALSO::

        - :class:`SpechtModule` for the generic diagram implementation
          constructed as a left ideal of the group algebra
        - :class:`~sage.combinat.symmetric_group_representations.SpechtRepresentation`
          for an implementation of the representation by matrices.
    """
    def __init__(self, ambient):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 5)
            sage: SM = SGA.specht_module([2,2,1])
            sage: TestSuite(SM).run()
        """
        self._diagram = ambient._shape
        SymmetricGroupRepresentation.__init__(self, ambient._semigroup_algebra)

        ambient_basis = ambient.basis()
        tabloids = ambient_basis.keys()
        support_order = list(tabloids)

        def elt(T):
            tab = tabloids.element_class(tabloids, list(T), check=False)
            return ambient.sum_of_terms((ambient._symmetric_group_action(tab, sigma), sigma.sign())
                                        for sigma in T.column_stabilizer())

        basis = Family({T: elt(T)
                        for T in self._diagram.standard_tableaux()})
        cat = ambient.category().Subobjects()
        SubmoduleWithBasis.__init__(self, basis, support_order, ambient=ambient,
                                    unitriangular=False, category=cat,
                                    prefix='S', bracket='')

    @lazy_attribute
    def lift(self):
        r"""
        The lift (embedding) map from ``self`` to the ambient space.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 5)
            sage: SM = SGA.specht_module([3, 1, 1])
            sage: SM.lift
            Generic morphism:
              From: Specht module of [3, 1, 1] over Rational Field
              To:   Tabloid module of [3, 1, 1] over Rational Field
        """
        return self.module_morphism(self.lift_on_basis, codomain=self.ambient())

    @lazy_attribute
    def retract(self):
        r"""
        The retract map from the ambient space.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 5)
            sage: X = SGA.tabloid_module([2,2,1])
            sage: Y = X.specht_module()
            sage: Y.retract
            Generic morphism:
              From: Tabloid module of [2, 2, 1] over Rational Field
              To:   Specht module of [2, 2, 1] over Rational Field
            sage: all(Y.retract(u.lift()) == u for u in Y.basis())
            True

            sage: Y.retract(X.zero())
            0
            sage: Y.retract(sum(X.basis()))
            Traceback (most recent call last):
            ...
            ValueError: ... is not in the image
        """
        B = self.basis()
        COB = matrix([b.lift().to_vector() for b in B]).T
        P, L, U = COB.LU()
        # Since U is upper triangular, the nonzero entriesm must be in the
        #   upper square portiion of the matrix
        n = len(B)

        Uinv = U.matrix_from_rows(range(n)).inverse()
        # This is a slight abuse as the codomain should be a module with a different
        #    S_n action, but we only use it internally, so there isn't any problems
        PLinv = (P*L).inverse()

        def retraction(elt):
            vec = PLinv * elt.to_vector(order=self._support_order)
            if not vec:
                return self.zero()
            # vec is now in the image of self under U, which is
            if max(vec.support()) >= n:
                raise ValueError(f"{elt} is not in the image")
            return self._from_dict(dict(zip(B.keys(), Uinv * vec[:n])))

        return self._ambient.module_morphism(function=retraction, codomain=self)

    def bilinear_form(self, u, v):
        r"""
        Return the natural bilinear form of ``self`` applied to ``u`` and ``v``.

        The natural bilinear form is given by the pullback of the natural
        bilinear form on the tabloid module (where the tabloid basis is an
        orthonormal basis).

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 5)
            sage: SM = SGA.specht_module([2,2,1])
            sage: u = SM.an_element(); u
            3*S[[1, 2], [3, 5], [4]] + 2*S[[1, 3], [2, 5], [4]] + 2*S[[1, 4], [2, 5], [3]]
            sage: v = sum(SM.basis())
            sage: SM.bilinear_form(u, v)
            140
        """
        TM = self._ambient
        return TM.bilinear_form(u.lift(), v.lift())

    @cached_method
    def gram_matrix(self):
        r"""
        Return the Gram matrix of the natural bilinear form of ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(QQ, 5)
            sage: SM = SGA.specht_module([2,2,1])
            sage: M = SM.gram_matrix(); M
            [12  4 -4 -4  4]
            [ 4 12  4  4  4]
            [-4  4 12  4  4]
            [-4  4  4 12  4]
            [ 4  4  4  4 12]
            sage: M.det() != 0
            True
        """
        B = self.basis()
        M = matrix([[self.bilinear_form(b, bp) for bp in B] for b in B])
        M.set_immutable()
        return M

    def maximal_submodule(self):
        """
        Return the maximal submodule of ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(GF(3), 5)
            sage: SM = SGA.specht_module([3,2])
            sage: U = SM.maximal_submodule()
            sage: U.dimension()
            4
        """
        return MaximalSpechtSubmodule(self)

    def simple_module(self):
        r"""
        Return the simple (or irreducible) `S_n`-submodule of ``self``.

        .. SEEALSO::

            :class:`~sage.combinat.specht_module.SimpleModule`

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(GF(3), 5)
            sage: SM = SGA.specht_module([3,2])
            sage: L = SM.simple_module()
            sage: L.dimension()
            1

            sage: SGA = SymmetricGroupAlgebra(QQ, 5)
            sage: SM = SGA.specht_module([3,2])
            sage: SM.simple_module() is SM
            True
        """
        if self.base_ring().characteristic() == 0:
            return self
        return SimpleModule(self)


class MaximalSpechtSubmodule(SubmoduleWithBasis, SymmetricGroupRepresentation):
    r"""
    The maximal submodule `U^{\lambda}` of the Specht module `S^{\lambda}`.

    ALGORITHM:

    We construct `U^{\lambda}` as the intersection `S \cap S^{\perp}`,
    where `S^{\perp}` is the orthogonal complement of the Specht module `S`
    inside of the tabloid module `T` (with respect to the natural
    bilinear form on `T`).

    EXAMPLES::

        sage: SGA = SymmetricGroupAlgebra(GF(3), 5)
        sage: SM = SGA.specht_module([3,2])
        sage: U = SM.maximal_submodule()
        sage: u = U.an_element(); u
        2*U[0] + 2*U[1]
        sage: [p * u for p in list(SGA.basis())[:4]]
        [2*U[0] + 2*U[1], 2*U[2] + 2*U[3], 2*U[0] + 2*U[1], U[0] + 2*U[2]]
        sage: sum(SGA.basis()) * u
        0
    """
    def __init__(self, specht_module):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(GF(3), 5)
            sage: SM = SGA.specht_module([3,2])
            sage: U = SM.maximal_submodule()
            sage: TestSuite(U).run()

            sage: SM = SGA.specht_module([2,1,1,1])
            sage: SM.maximal_submodule()
            Traceback (most recent call last):
            ...
            NotImplementedError: only implemented for 3-regular partitions

            sage: SGA = SymmetricGroupAlgebra(QQ, 5)
            sage: SM = SGA.specht_module([3,2])
            sage: U = SM.maximal_submodule()
            sage: TestSuite(U).run()
            sage: U.dimension()
            0
        """
        SymmetricGroupRepresentation.__init__(self, specht_module._semigroup_algebra)

        p = specht_module.base_ring().characteristic()
        if p == 0:
            basis = Family([])
        else:
            TM = specht_module._ambient
            if not TM._shape.is_regular(p):
                raise NotImplementedError(f"only implemented for {p}-regular partitions")
            TV = TM._dense_free_module()
            SV = TV.submodule(specht_module.lift.matrix().columns())
            basis = (SV & SV.complement()).basis()
            basis = [specht_module.retract(TM.from_vector(b)) for b in basis]
            basis = Family(specht_module.echelon_form(basis))

        unitriangular = all(b.leading_support() == 1 for b in basis)
        support_order = list(specht_module.basis().keys())
        cat = specht_module.category().Subobjects()
        SubmoduleWithBasis.__init__(self, basis, support_order, ambient=specht_module,
                                    unitriangular=unitriangular, category=cat,
                                    prefix='U')

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(GF(3), 5)
            sage: SM = SGA.specht_module([3,2])
            sage: SM.maximal_submodule()
            Maximal submodule of Specht module of [3, 2] over Finite Field of size 3
        """
        return f"Maximal submodule of {self._ambient}"

    Element = SpechtModule.Element


class SimpleModule(QuotientModuleWithBasis, SymmetricGroupRepresentation):
    r"""
    The simgle `S_n`-module associated with a partition `\lambda`.

    The simple module `D^{\lambda}` is the quotient of the Specht module
    `S^{\lambda}` by its :class:`maximal submodule <MaximalSpechtSubmodule>`
    `U^{\lambda}`.

    EXAMPLES::

        sage: SGA = SymmetricGroupAlgebra(GF(3), 5)
        sage: SM = SGA.specht_module([3,1,1])
        sage: D = SM.simple_module()
        sage: v = D.an_element(); v
        2*D[[[1, 3, 5], [2], [4]]] + 2*D[[[1, 4, 5], [2], [3]]]
        sage: SGA.an_element() * v
        2*D[[[1, 2, 4], [3], [5]]] + 2*D[[[1, 3, 5], [2], [4]]]

    We give an example on how to construct the decomposition matrix
    (the Specht modules are a complete set of irreducible projective
    modules) and the Cartan matrix of a symmetric group algebra::

        sage: SGA = SymmetricGroupAlgebra(GF(3), 4)
        sage: BM = matrix(SGA.simple_module(la).brauer_character()
        ....:             for la in Partitions(4, regular=3))
        sage: SBT = matrix(SGA.specht_module(la).brauer_character()
        ....:              for la in Partitions(4))
        sage: D = SBT * ~BM; D
        [1 0 0 0]
        [0 1 0 0]
        [1 0 1 0]
        [0 0 0 1]
        [0 0 1 0]
        sage: D.transpose() * D
        [2 0 1 0]
        [0 1 0 0]
        [1 0 2 0]
        [0 0 0 1]

    We verify this against the direct computation (up to reindexing the
    rows and columns)::

        sage: SGA.cartan_invariants_matrix()  # long time
        [1 0 0 0]
        [0 1 0 0]
        [0 0 2 1]
        [0 0 1 2]
    """
    def __init__(self, specht_module):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(GF(3), 5)
            sage: SM = SGA.specht_module([3,1,1])
            sage: D = SM.simple_module()
            sage: TestSuite(D).run()

            sage: SGA = SymmetricGroupAlgebra(GF(3), 5)
            sage: SM = SGA.specht_module([2,1,1,1])
            sage: SM.simple_module()
            Traceback (most recent call last):
            ...
            ValueError: the partition must be 3-regular
        """
        self._diagram = specht_module._diagram
        p = specht_module.base_ring().characteristic()
        if not self._diagram.is_regular(p):
            raise ValueError(f"the partition must be {p}-regular")
        SymmetricGroupRepresentation.__init__(self, specht_module._semigroup_algebra)
        cat = specht_module.category()
        QuotientModuleWithBasis.__init__(self, specht_module.maximal_submodule(), cat, prefix='D')

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: SGA = SymmetricGroupAlgebra(GF(3), 5)
            sage: SM = SGA.specht_module([3,1,1])
            sage: SM.simple_module()
            Simple module of [3, 1, 1] over Finite Field of size 3
        """
        return f"Simple module of {self._diagram} over {self.base_ring()}"

    Element = SpechtModule.Element


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

    TESTS:

    Verify that diagrams bigger than the rank work::

        sage: specht_module_spanning_set([(0,0), (3,5)])
        ([1, 2], [2, 1])
        sage: specht_module_spanning_set([(0,0), (5,3)])
        ([1, 2], [2, 1])
    """
    n = len(D)
    if SGA is None:
        from sage.combinat.symmetric_group_algebra import SymmetricGroupAlgebra
        SGA = SymmetricGroupAlgebra(QQ, n)
    elif SGA.group().rank() != n - 1:
        raise ValueError("the rank does not match the size of the diagram")
    nr = max((c[0] for c in D), default=0) + 1
    nc = max((c[1] for c in D), default=0) + 1
    row_diagram = [set() for _ in range(nr)]
    col_diagram = [set() for _ in range(nc)]
    for i, cell in enumerate(D):
        x, y = cell
        row_diagram[x].add(i)
        col_diagram[y].add(i)
    # Construct the row and column stabilizer elements
    row_stab = SGA.zero()
    col_stab = SGA.zero()
    B = SGA.basis()
    for w in B.keys():
        # Remember that the permutation w is 1-based
        row_perm = [set() for _ in range(nr)]
        col_perm = [set() for _ in range(nc)]
        for i, cell in enumerate(D):
            x, y = cell
            row_perm[x].add(w(i + 1) - 1)
            col_perm[y].add(w(i + 1) - 1)
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


def polytabloid(T):
    r"""
    Compute the polytabloid element associated to a tableau ``T``.

    For a tableau `T`, the polytabloid associated to `T` is

    .. MATH::

        e_T = \sum_{\sigma \in C_T} (-1)^{\sigma} \{\sigma T\},

    where `\{\}` is the row-equivalence class, i.e. a tabloid,
    and `C_T` is the column stabilizer of `T`. The sum takes place in
    the module spanned by tabloids `\{T\}`.

    OUTPUT:

    A ``dict`` whose keys are tabloids represented by tuples of frozensets
    and whose values are the coefficient.

    EXAMPLES::

        sage: from sage.combinat.specht_module import polytabloid
        sage: T = StandardTableau([[1,3,4],[2,5]])
        sage: polytabloid(T)
        {(frozenset({1, 3, 4}), frozenset({2, 5})): 1,
         (frozenset({1, 4, 5}), frozenset({2, 3})): -1,
         (frozenset({2, 3, 4}), frozenset({1, 5})): -1,
         (frozenset({2, 4, 5}), frozenset({1, 3})): 1}
    """
    e_T = {}
    C_T = T.column_stabilizer()
    for perm in C_T:
        TT = tuple([frozenset(perm(val) for val in row) for row in T])
        if TT in e_T:
            e_T[TT] += perm.sign()
        else:
            e_T[TT] = perm.sign()
    return e_T


def tabloid_gram_matrix(la, base_ring):
    r"""
    Compute the Gram matrix of the bilinear form of a Specht module
    pulled back from the tabloid module.

    For the module spanned by all tabloids, we define an bilinear form
    by having the tabloids be an orthonormal basis. We then pull this
    bilinear form back across the natural injection of the Specht module
    into the tabloid module.

    EXAMPLES::

        sage: from sage.combinat.specht_module import tabloid_gram_matrix
        sage: tabloid_gram_matrix([3,2], GF(5))
        [4 2 2 1 4]
        [2 4 1 2 1]
        [2 1 4 2 1]
        [1 2 2 4 2]
        [4 1 1 2 4]
    """
    from sage.combinat.tableau import StandardTableaux
    ST = list(StandardTableaux(la))

    def bilinear_form(p1, p2):
        if len(p2) < len(p1):
            p1, p2 = p2, p1
        return sum(c1 * p2[T1] for T1, c1 in p1.items() if c1 and T1 in p2)

    PT = {T: polytabloid(T) for T in ST}
    gram_matrix = [[bilinear_form(PT[T1], PT[T2]) for T1 in ST] for T2 in ST]
    return matrix(base_ring, gram_matrix)


def simple_module_rank(la, base_ring):
    r"""
    Return the rank of the simple `S_n`-module corresponding to the
    partition ``la`` of size `n` over ``base_ring``.

    EXAMPLES::

        sage: from sage.combinat.specht_module import simple_module_rank
        sage: simple_module_rank([3,2,1,1], GF(3))
        13

    TESTS::

        sage: from sage.combinat.specht_module import simple_module_rank
        sage: simple_module_rank([1,1,1,1], GF(3))
        Traceback (most recent call last):
        ...
        ValueError: the partition [1, 1, 1, 1] is not 3-regular

        sage: from sage.combinat.specht_module import simple_module_rank
        sage: simple_module_rank([2,1], GF(3)['x'])
        Traceback (most recent call last):
        ...
        NotImplementedError: the base must be a field
    """
    from sage.categories.fields import Fields
    from sage.combinat.partition import Partition
    if base_ring not in Fields():
        raise NotImplementedError("the base must be a field")
    p = base_ring.characteristic()
    la = Partition(la)
    if not la.is_regular(p):
        raise ValueError(f"the partition {la} is not {p}-regular")
    return tabloid_gram_matrix(la, base_ring).rank()
