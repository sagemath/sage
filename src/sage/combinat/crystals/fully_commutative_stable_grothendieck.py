r"""
Fully commutative stable Grothendieck crystal

AUTHORS:

- Jianping Pan (2020-08-31): initial version

- Wencin Poh (2020-08-31): initial version

- Anne Schilling (2020-08-31): initial version
"""

# ****************************************************************************
#       Copyright (C) 2020 Jianping Pan <jppan at math dot ucdavis dot edu>
#                          Wencin Poh <wpoh at ucdavis dot edu>
#                          Anne Schilling <anne at math dot ucdavis dot edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.structure.richcmp import richcmp_by_eq_and_lt
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.classical_crystals import ClassicalCrystals
from sage.categories.enumerated_sets import EnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat import permutation
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.rings.integer import Integer
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.lazy_attribute import lazy_attribute

class DecreasingHeckeFactorization(Element, metaclass=InheritComparisonClasscallMetaclass):
    """
    Class of decreasing factorizations in the 0-Hecke monoid.

    INPUT:

    - ``t`` -- decreasing factorization inputted as list of lists

    - ``max_value`` -- maximal value of entries

    EXAMPLES::

        sage: from sage.combinat.crystals.fully_commutative_stable_grothendieck import DecreasingHeckeFactorization
        sage: t = [[3, 2], [], [2, 1]]
        sage: h = DecreasingHeckeFactorization(t, 3); h
        (3, 2)()(2, 1)
        sage: h.excess
        1
        sage: h.factors
        3
        sage: h.max_value
        3
        sage: h.value
        ((3, 2), (), (2, 1))

        sage: u = [[3, 2, 1], [3], [2, 1]]
        sage: h = DecreasingHeckeFactorization(u); h
        (3, 2, 1)(3)(2, 1)
        sage: h.weight()
        (2, 1, 3)
        sage: h.parent()
        Decreasing Hecke factorizations with 3 factors associated to [2, 1, 3, 2, 1] with excess 1
    """
    @staticmethod
    def __classcall_private__(self, t, max_value=None, parent=None):
        """
        Assign the correct parent for ``t`` and call ``t`` as an element of that parent

        EXAMPLES::

            sage: from sage.combinat.crystals.fully_commutative_stable_grothendieck import DecreasingHeckeFactorization
            sage: h1 = DecreasingHeckeFactorization([[3, 1], [], [3, 2]])
            sage: h1.parent()
            Fully commutative stable Grothendieck crystal of type A_2 associated to [1, 3, 2] with excess 1
            sage: h2 = DecreasingHeckeFactorization(h1)
            sage: h1 == h2
            True

            sage: h1 = DecreasingHeckeFactorization([[3, 1], [2, 1], [2, 1]])
            sage: F = h1.parent(); F
            Decreasing Hecke factorizations with 3 factors associated to [1, 3, 2, 1] with excess 2
            sage: h2 = F(h1)
            sage: h1 == h2
            True

        TESTS::

            sage: DecreasingHeckeFactorization([[]])
            ()
        """
        _check_decreasing_hecke_factorization(t)
        if isinstance(t, DecreasingHeckeFactorization):
            u = t.value
            if parent is None:
                parent = t.parent()
        else:
            u = t
            if parent is None:
                if max_value is None:
                    letters = [x for factor in t for x in factor]
                    max_value = max(letters) if letters else 1
                from sage.monoids.hecke_monoid import HeckeMonoid
                S = SymmetricGroup(max_value+1)
                H = HeckeMonoid(S)
                word = H.from_reduced_word(x for factor in t for x in factor).reduced_word()
                factors = len(t)
                excess = sum(len(l) for l in t) - len(word)

                p = permutation.from_reduced_word(word)
                if p.has_pattern([3,2,1]):
                    word = S.from_reduced_word(word)
                    parent = DecreasingHeckeFactorizations(word, factors, excess)
                else:
                    word = S.from_reduced_word(word)
                    parent = FullyCommutativeStableGrothendieckCrystal(word, factors, excess)
        return parent.element_class(parent, u)

    def __init__(self, parent, t):
        """
        Initialize a decreasing factorization for ``self`` given the relevant data.

        EXAMPLES::

            sage: from sage.combinat.crystals.fully_commutative_stable_grothendieck import DecreasingHeckeFactorization
            sage: t = [[2, 1], [2], [], [1]]
            sage: h1 = DecreasingHeckeFactorization(t); h1
            (2, 1)(2)()(1)
            sage: h1.excess
            1
            sage: h2 = DecreasingHeckeFactorization(t, 2)
            sage: h2.value
            ((2, 1), (2,), (), (1,))
            sage: h1 == h2
            True

            sage: t = [[2, 1], [2], [], [3, 1]]
            sage: h = DecreasingHeckeFactorization(t, 5)
            sage: h.max_value
            5
            sage: h.factors
            4
            sage: h.w
            (1, 2, 1, 3)

        TESTS::

            sage: from sage.combinat.crystals.fully_commutative_stable_grothendieck import DecreasingHeckeFactorization
            sage: t1 = [[], [3, 1], [], [2, 1], [1]]
            sage: t2 = [[], [3, 1], [], [2, 1], [2]]
            sage: h1 = DecreasingHeckeFactorization(t1, 3)
            sage: h2 = DecreasingHeckeFactorization(t2, 3)
            sage: h3 = DecreasingHeckeFactorization(t1)
            sage: h1 == h2, h1 == h3
            (False, True)
            sage: h1 != h2, h1 != h3
            (True, False)
        """
        Element.__init__(self, parent)
        self.factors = parent.factors
        self.max_value = parent.max_value
        self.w = parent.w
        self.excess = parent.excess
        self.value = tuple(tuple(factors) for factors in t)

    def _repr_(self):
        """
        Return the representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.crystals.fully_commutative_stable_grothendieck import DecreasingHeckeFactorization
            sage: t = [[], [2, 1], [2], [], [2]]
            sage: h = DecreasingHeckeFactorization(t); h
            ()(2, 1)(2)()(2)
        """
        return "".join("("+repr(list(factor))[1:-1]+")" for factor in self.value)

    def __hash__(self):
        """
        Return hash of ``self``.

        EXAMPLES::

            sage: from sage.combinat.crystals.fully_commutative_stable_grothendieck import DecreasingHeckeFactorization
            sage: t = [[], [2, 1], [2], [], [2]]
            sage: h1 = DecreasingHeckeFactorization(t)
            sage: h2 = DecreasingHeckeFactorization(t, 3)
            sage: h3 = DecreasingHeckeFactorization(t, 2)
            sage: hash(h1) == hash(h2)
            False
            sage: hash(h1) == hash(h3)
            True
        """
        return hash((self.max_value, self.value))

    _richcmp_ = richcmp_by_eq_and_lt("__eq__", "__lt__")

    def __eq__(self, other):
        """
        Return True if ``self`` equals ``other`` and False otherwise.

        EXAMPLES::

            sage: from sage.combinat.crystals.fully_commutative_stable_grothendieck import DecreasingHeckeFactorization
            sage: t = [[], [2, 1], [2], [], [2]]
            sage: h1 = DecreasingHeckeFactorization(t)
            sage: h2 = DecreasingHeckeFactorization(t, 3)
            sage: h1 == h2
            True
        """
        return isinstance(self, type(other)) and self.value == other.value

    def __lt__(self,other):
        """
        Return True if ``self`` comes before ``other`` and False otherwise.

        We say that `h_1` comes before `h_2` if either weight of `h_1 <` weight of `h_2`
        lexicographically, or if both weights of `h_1` and `h_2` are equal,
        but `h_1 < h_2` lexicographically.
        This ordering is mainly used for sorting or comparison.

        EXAMPLES::

            sage: from sage.combinat.crystals.fully_commutative_stable_grothendieck import DecreasingHeckeFactorization
            sage: t1 = [[], [2, 1], [], [2, 1], [1]]
            sage: t2 = [[], [2, 1], [], [2, 1], [2]]
            sage: t3 = [[], [2, 1], [2], [1], [1]]
            sage: h1 = DecreasingHeckeFactorization(t1)
            sage: h2 = DecreasingHeckeFactorization(t2)
            sage: h3 = DecreasingHeckeFactorization(t3)
            sage: h1 < h2
            True
            sage: h1 < h3
            False
        """
        return (self.weight(), self.value) < (other.weight(), other.value)

    def _latex_(self):
        r"""
        Return LaTeX code for ``self``.

        EXAMPLES::

            sage: from sage.combinat.crystals.fully_commutative_stable_grothendieck import DecreasingHeckeFactorization
            sage: t = [[2], [2, 1], [], [4, 3, 1]]
            sage: h = DecreasingHeckeFactorization(t, 6)
            sage: latex(h)
            \left(2\right)\left(2, 1\right)\left(\;\right)\left(4, 3, 1\right)
        """
        s = ""
        for factor in self.value:
            if factor:
                s += r"\left("+repr(list(factor))[1:-1]+r"\right)"
            else:
                s += r"\left(\;\right)"
        return s

    def weight(self):
        """
        Return the weight of ``self``.

        EXAMPLES::

            sage: from sage.combinat.crystals.fully_commutative_stable_grothendieck import DecreasingHeckeFactorization
            sage: t = [[2], [2, 1], [], [4, 3, 1]]
            sage: h = DecreasingHeckeFactorization(t, 6)
            sage: h.weight()
            (3, 0, 2, 1)
        """
        return tuple([len(l) for l in reversed(self.value)])

    def to_word(self):
        """
        Return the word associated to ``self`` in the 0-Hecke monoid.

        EXAMPLES::

            sage: from sage.combinat.crystals.fully_commutative_stable_grothendieck import DecreasingHeckeFactorization
            sage: t = [[2], [], [2, 1], [4, 3, 1]]
            sage: h = DecreasingHeckeFactorization(t)
            sage: h.to_word()
            [2, 2, 1, 4, 3, 1]
        """
        return [j for factors in self.value for j in factors]

    def to_increasing_hecke_biword(self):
        """
        Return the associated increasing Hecke biword of ``self``.

        EXAMPLES::

            sage: from sage.combinat.crystals.fully_commutative_stable_grothendieck import DecreasingHeckeFactorization
            sage: t = [[2], [], [2, 1],[4, 3, 1]]
            sage: h = DecreasingHeckeFactorization(t, 4)
            sage: h.to_increasing_hecke_biword()
            [[1, 1, 1, 2, 2, 4], [1, 3, 4, 1, 2, 2]]
        """
        L = [[],[]]
        for j in range(len(self.value)):
            L[1] += list(self.value[-j-1][::-1])
            L[0] += [j+1]*len(self.value[-j-1])
        return L

class DecreasingHeckeFactorizations(UniqueRepresentation, Parent):
    """
    Set of decreasing factorizations in the 0-Hecke monoid.

    INPUT:

    - ``w`` -- an element in the symmetric group

    - ``factors`` -- the number of factors in the factorization

    - ``excess`` -- the total number of letters in the factorization minus the length of a reduced word for ``w``

    EXAMPLES::

        sage: from sage.combinat.crystals.fully_commutative_stable_grothendieck import DecreasingHeckeFactorizations
        sage: S = SymmetricGroup(3+1)
        sage: w = S.from_reduced_word([1, 3, 2, 1])
        sage: F = DecreasingHeckeFactorizations(w, 3, 3); F
        Decreasing Hecke factorizations with 3 factors associated to [1, 3, 2, 1] with excess 3
        sage: F.list()
        [(3, 1)(3, 1)(3, 2, 1), (3, 1)(3, 2, 1)(2, 1), (3, 2, 1)(2, 1)(2, 1)]
    """
    @staticmethod
    def __classcall_private__(cls, w, factors, excess):
        """
        Classcall to mend the input.

        EXAMPLES::

            sage: from sage.combinat.crystals.fully_commutative_stable_grothendieck import DecreasingHeckeFactorizations
            sage: S = SymmetricGroup(3+1)
            sage: w = S.from_reduced_word([1, 2 ,1])
            sage: F = DecreasingHeckeFactorizations(w, 4, 1); F
            Decreasing Hecke factorizations with 4 factors associated to [1, 2, 1] with excess 1

            sage: from sage.monoids.hecke_monoid import HeckeMonoid
            sage: H = HeckeMonoid(SymmetricGroup(3+1))
            sage: w = H.from_reduced_word([1, 2 ,1])
            sage: G = DecreasingHeckeFactorizations(w, 4, 1); G
            Decreasing Hecke factorizations with 4 factors associated to [1, 2, 1] with excess 1
            sage: F is G
            True
        """
        from sage.monoids.hecke_monoid import HeckeMonoid
        if isinstance(w.parent(), SymmetricGroup):
            H = HeckeMonoid(w.parent())
            w = H.from_reduced_word(w.reduced_word())
        if (not w.reduced_word()) and excess!=0:
            raise ValueError("excess must be 0 for the empty word")
        return super(DecreasingHeckeFactorizations, cls).__classcall__(cls, w, factors, excess)

    def __init__(self, w, factors, excess):
        """
        Initialize a set for ``self`` given reduced word ``w`` in the symmetric group,
        number of factors ``factors`` and``excess`` extra letters.

        EXAMPLES::

            sage: from sage.combinat.crystals.fully_commutative_stable_grothendieck import DecreasingHeckeFactorizations
            sage: S = SymmetricGroup(3+1)
            sage: w = S.from_reduced_word([2, 1, 3, 2, 1])
            sage: F = DecreasingHeckeFactorizations(w, 4, 2)
            sage: F.w
            (2, 1, 3, 2, 1)
            sage: F.factors
            4
            sage: F.excess
            2
            sage: F.H
            0-Hecke monoid of the Symmetric group of order 4! as a permutation group

        TESTS::

            sage: from sage.combinat.crystals.fully_commutative_stable_grothendieck import DecreasingHeckeFactorizations
            sage: S = SymmetricGroup(2+1)
            sage: w = S.from_reduced_word([1, 2, 1])
            sage: F = DecreasingHeckeFactorizations(w, 3, 1)
            sage: TestSuite(F).run()
        """
        Parent.__init__(self, category = FiniteEnumeratedSets())
        self.w = tuple(w.reduced_word())
        self.factors = factors
        self.H = w.parent()
        self.max_value = len(self.H.gens())
        self.excess = excess

    def _repr_(self):
        """
        Return a representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.crystals.fully_commutative_stable_grothendieck import DecreasingHeckeFactorizations
            sage: S = SymmetricGroup(3+1)
            sage: w = S.from_reduced_word([2, 1, 3, 2])
            sage: DecreasingHeckeFactorizations(w, 3, 1)
            Decreasing Hecke factorizations with 3 factors associated to [2, 1, 3, 2] with excess 1
        """
        return "Decreasing Hecke factorizations with {} factors associated to {} with excess {}".format(self.factors, list(self.w), self.excess)

    def list(self):
        """
        Return list of all elements of ``self``.

        EXAMPLES::

            sage: from sage.combinat.crystals.fully_commutative_stable_grothendieck import DecreasingHeckeFactorizations
            sage: S = SymmetricGroup(3+1)
            sage: w = S.from_reduced_word([1, 3, 2, 1])
            sage: F = DecreasingHeckeFactorizations(w, 3, 3)
            sage: F.list()
            [(3, 1)(3, 1)(3, 2, 1), (3, 1)(3, 2, 1)(2, 1), (3, 2, 1)(2, 1)(2, 1)]
        """
        return _generate_decreasing_hecke_factorizations(self.w, self.factors, self.excess, parent=self)

    # temporary workaround while an_element is overridden by Parent
    _an_element_ = EnumeratedSets.ParentMethods._an_element_

    Element = DecreasingHeckeFactorization


class FullyCommutativeStableGrothendieckCrystal(UniqueRepresentation, Parent):
    """
    The crystal on fully commutative decreasing factorizations in the 0-Hecke
    monoid, as introduced by [MPPS2020]_.

    INPUT:

    - ``w`` -- an element in the symmetric group or a (skew) shape

    - ``factors`` -- the number of factors in the factorization

    - ``excess`` -- the total number of letters in the factorization minus the length of a reduced word for ``w``

    - ``shape`` -- (default: ``False``) indicator for input ``w``, True if ``w`` is entered as a (skew) shape and False otherwise.

    EXAMPLES::

        sage: S = SymmetricGroup(3+1)
        sage: w = S.from_reduced_word([1, 3, 2])
        sage: B = crystals.FullyCommutativeStableGrothendieck(w, 3, 2); B
        Fully commutative stable Grothendieck crystal of type A_2 associated to [1, 3, 2] with excess 2
        sage: B.list()
        [(1)(3, 1)(3, 2),
         (3, 1)(1)(3, 2),
         (3, 1)(3, 1)(2),
         (3)(3, 1)(3, 2),
         (3, 1)(3)(3, 2),
         (3, 1)(3, 2)(2)]

    We can also access the crystal by specifying a skew shape::

        sage: crystals.FullyCommutativeStableGrothendieck([[2, 2], [1]], 4, 1, shape=True)
        Fully commutative stable Grothendieck crystal of type A_3 associated to [2, 1, 3] with excess 1

    We can compute the highest weight elements::

        sage: hw = [w for w in B if w.is_highest_weight()]
        sage: hw
        [(1)(3, 1)(3, 2), (3)(3, 1)(3, 2)]
        sage: hw[0].weight()
        (2, 2, 1)

    The crystal operators themselves move elements between adjacent factors::

        sage: b = hw[0]; b
        (1)(3, 1)(3, 2)
        sage: b.f(2)
        (3, 1)(1)(3, 2)
    """
    @staticmethod
    def __classcall_private__(cls, w, factors, excess, shape=False):
        """
        Classcall to mend the input.

        EXAMPLES::

            sage: A = crystals.FullyCommutativeStableGrothendieck([[3, 3], [2, 1]], 4, 1, shape=True); A
            Fully commutative stable Grothendieck crystal of type A_3 associated to [3, 2, 4] with excess 1
            sage: B = crystals.FullyCommutativeStableGrothendieck(SkewPartition([[3, 3], [2, 1]]), 4, 1, shape=True)
            sage: A is B
            True

            sage: C = crystals.FullyCommutativeStableGrothendieck((2, 1), 3, 2, shape=True); C
            Fully commutative stable Grothendieck crystal of type A_2 associated to [1, 3, 2] with excess 2
            sage: D = crystals.FullyCommutativeStableGrothendieck(Partition([2, 1]), 3, 2, shape=True)
            sage: C is D
            True
        """
        from sage.monoids.hecke_monoid import HeckeMonoid
        if shape:
            from sage.combinat.partition import _Partitions
            from sage.combinat.skew_partition import SkewPartition
            cond1 = isinstance(w, (tuple, list)) and len(w)==2 and w[0] in _Partitions and w[1] in _Partitions
            cond2 = isinstance(w, SkewPartition)
            if cond1 or cond2:
                sh = SkewPartition([w[0], w[1]])
            elif w in _Partitions:
                sh = SkewPartition([w, []])
            else:
                raise ValueError("w needs to be a (skew) partition")
            word = _to_reduced_word(sh)
            max_value = max(word) if word else 1
            H = HeckeMonoid(SymmetricGroup(max_value+1))
            w = H.from_reduced_word(word)
        else:
            if isinstance(w.parent(), SymmetricGroup):
                H = HeckeMonoid(w.parent())
                w = H.from_reduced_word(w.reduced_word())
        if (not w.reduced_word()) and excess!=0:
            raise ValueError("excess must be 0 for the empty word")
        return super(FullyCommutativeStableGrothendieckCrystal, cls).__classcall__(cls, w, factors, excess)

    def __init__(self, w, factors, excess):
        """
        Initialize a crystal for ``self`` given reduced word ``w`` in the symmetric group,
        number of factors ``factors`` and``excess`` extra letters.

        EXAMPLES::

            sage: S = SymmetricGroup(3+1)
            sage: w = S.from_reduced_word([1, 3, 2])
            sage: B = crystals.FullyCommutativeStableGrothendieck(w, 3, 2)
            sage: B.w
            (1, 3, 2)
            sage: B.factors
            3
            sage: B.excess
            2
            sage: B.H
            0-Hecke monoid of the Symmetric group of order 4! as a permutation group

        The reduced word ``w`` should be fully commutative, that is, its
        associated permutation should avoid the pattern 321::

            sage: S = SymmetricGroup(3+1)
            sage: w = S.from_reduced_word([1, 2, 1])
            sage: B = crystals.FullyCommutativeStableGrothendieck(w, 4, 2)
            Traceback (most recent call last):
            ...
            ValueError: w should be fully commutative

        TESTS::

            sage: S = SymmetricGroup(3+1)
            sage: w = S.from_reduced_word([2, 3, 1])
            sage: B = crystals.FullyCommutativeStableGrothendieck(w, 4, 2)
            sage: TestSuite(B).run()
        """
        # Check if w is fully commutative
        word = w.reduced_word()
        p = permutation.from_reduced_word(word)
        if p.has_pattern([3,2,1]):
            raise ValueError("w should be fully commutative")

        Parent.__init__(self, category = ClassicalCrystals())
        self.w = tuple(word)
        self.factors = factors
        self.H = w.parent()
        self.max_value = len(self.H.gens())
        self.excess = excess
        self._cartan_type = CartanType(['A', self.factors-1])

    @lazy_attribute
    def module_generators(self):
        """
        Return generators for ``self`` as a crystal.

        EXAMPLES::

            sage: S = SymmetricGroup(3+1)
            sage: w = S.from_reduced_word([1, 3, 2])
            sage: B = crystals.FullyCommutativeStableGrothendieck(w, 3, 2)
            sage: B.module_generators
            ((1)(3, 1)(3, 2), (3)(3, 1)(3, 2))
            sage: C = crystals.FullyCommutativeStableGrothendieck(w, 4, 2)
            sage: C.module_generators
            (()(1)(3, 1)(3, 2),
             ()(3)(3, 1)(3, 2),
             (1)(1)(1)(3, 2),
             (1)(1)(3)(3, 2),
             (1)(3)(3)(3, 2))
        """
        return tuple(self(x).to_highest_weight()[0] for x in _lowest_weights(self.w, self.factors, self.excess, parent=self))

    def _repr_(self):
        """
        Return a representation of ``self``.

        EXAMPLES::

            sage: S = SymmetricGroup(3+1)
            sage: w = S.from_reduced_word([2, 1, 3, 2])
            sage: crystals.FullyCommutativeStableGrothendieck(w, 3, 1)
            Fully commutative stable Grothendieck crystal of type A_2 associated to [2, 1, 3, 2] with excess 1
        """
        return "Fully commutative stable Grothendieck crystal of type A_{} associated to {} with excess {}".format(self.factors-1, list(self.w), self.excess)

    class Element(DecreasingHeckeFactorization):
        def __init__(self, parent, t):
            """
            Create an instance ``self`` of element ``t``.

            This method takes into account the constraints on the word,
            the number of factors, and excess statistic associated to the parent class.

            EXAMPLES::

                sage: S = SymmetricGroup(3+1)
                sage: w = S.from_reduced_word([1, 3, 2])
                sage: B = crystals.FullyCommutativeStableGrothendieck(w, 3, 2)
                sage: from sage.combinat.crystals.fully_commutative_stable_grothendieck import DecreasingHeckeFactorization
                sage: h = DecreasingHeckeFactorization([[3, 1], [3], [3, 2]], 4)
                sage: u = B(h); u.value
                ((3, 1), (3,), (3, 2))
                sage: v = B([[3, 1], [3], [3, 2]]); v.value
                ((3, 1), (3,), (3, 2))
            """
            u = t
            if isinstance(t, DecreasingHeckeFactorization):
                u = t.value
            _check_decreasing_hecke_factorization(t)
            _check_containment(t, parent)
            DecreasingHeckeFactorization.__init__(self, parent, u)

        def e(self, i):
            """
            Return the action of `e_i` on ``self`` using the rules described in [MPPS2020]_.

            EXAMPLES::

                sage: S = SymmetricGroup(4+1)
                sage: w = S.from_reduced_word([2, 1, 4, 3, 2])
                sage: B = crystals.FullyCommutativeStableGrothendieck(w, 4, 3)
                sage: h = B([[4, 2], [4, 2, 1], [3, 2], [2]]); h
                (4, 2)(4, 2, 1)(3, 2)(2)
                sage: h.e(1)
                (4, 2)(4, 2, 1)(3)(3, 2)
                sage: h.e(2)
                (4, 2)(2, 1)(4, 3, 2)(2)
                sage: h.e(3)
            """
            P = self.parent()
            m = P.factors
            L = list(self.value[m-i-1])
            R = list(self.value[m-i])
            b = self.bracketing(i)
            if not b[0]:
                return None
            y = b[0][-1]
            if y-1 in L and y-1 in R:
                # special case: (--x+1--)(--x+1,x--) -->> (--x+1,x--)(--x--)
                L.remove(y-1)
            else:
                L.remove(y)
            R.append(y)
            L.sort(reverse=True)
            R.sort(reverse=True)
            s = [self.value[j] for j in range(m-i-1)]+[L]+[R]+[self.value[j] for j in range(m-i+1, m)]
            return P.element_class(P, s)

        def f(self, i):
            """
            Return the action of `f_i` on ``self`` using the rules described in [MPPS2020]_.

            EXAMPLES::

                sage: S = SymmetricGroup(4+1)
                sage: w = S.from_reduced_word([3, 2, 1, 4, 3])
                sage: B = crystals.FullyCommutativeStableGrothendieck(w, 4, 3)
                sage: h = B([[3, 2], [2, 1], [4, 3], [3, 1]]); h
                (3, 2)(2, 1)(4, 3)(3, 1)
                sage: h.f(1)
                (3, 2)(2, 1)(4, 3, 1)(3)
                sage: h.f(2)
                sage: h.f(3)
                (3, 2, 1)(1)(4, 3)(3, 1)
            """
            P = self.parent()
            m = P.factors
            L = list(self.value[m-i-1])
            R = list(self.value[m-i])
            b = self.bracketing(i)
            if not b[1]:
                return None
            x = b[1][0]
            if x+1 in L and x+1 in R:
                # special case: (--x+1--)(--x+1,x--) -->> (--x+1,x--)(--x--)
                R.remove(x+1)
            else:
                R.remove(x)
            L.append(x)
            L.sort(reverse=True)
            R.sort(reverse=True)
            s = [self.value[j] for j in range(m-i-1)]+[L]+[R]+[self.value[j] for j in range(m-i+1, m)]
            return P.element_class(P, s)

        def bracketing(self, i):
            """
            Remove all bracketed letters between `i`-th and `(i+1)`-th entry.

            EXAMPLES::

                sage: S = SymmetricGroup(4+1)
                sage: w = S.from_reduced_word([3, 2, 1, 4, 3])
                sage: B = crystals.FullyCommutativeStableGrothendieck(w, 3, 2)
                sage: h = B([[3], [4, 2, 1], [4, 3, 1]])
                sage: h.bracketing(1)
                [[], []]
                sage: h.bracketing(2)
                [[], [2, 1]]
            """
            P = self.parent()
            m = P.factors
            L = list(self.value[m-i-1])
            R = list(self.value[m-i])
            right_n = [j for j in R]
            left_n = [j for j in L]
            left_unbracketed = []
            while left_n:
                m = max(left_n)
                left_n.remove(m)
                l = [j for j in right_n if j>=m]
                if l:
                    right_n.remove(min(l))
                else:
                    left_unbracketed += [m]
            return [[j for j in left_unbracketed], [j for j in right_n]]


####################
# Helper functions #
####################

def _check_decreasing_hecke_factorization(t):
    """
    Check if ``t`` is a suitable data type for a decreasing factorization in a 0-Hecke monoid.

    EXAMPLES::

        sage: from sage.combinat.crystals.fully_commutative_stable_grothendieck import _check_decreasing_hecke_factorization
        sage: _check_decreasing_hecke_factorization([[3, 2], [2, 1], [4]])
        sage: _check_decreasing_hecke_factorization([[3, 2, 2], [2, 1], [4]])
        Traceback (most recent call last):
        ...
        ValueError: each nonempty factor should be a strictly decreasing sequence
        sage: _check_decreasing_hecke_factorization([[3, 'a'], [2, 1], [4]])
        Traceback (most recent call last):
        ...
        ValueError: each nonempty factor should contain integers
        sage: _check_decreasing_hecke_factorization([[3, 2], [2, 1], 4])
        Traceback (most recent call last):
        ...
        ValueError: each factor in t should be a list or tuple
    """
    if not isinstance(t, DecreasingHeckeFactorization):
        if not isinstance(t, (tuple, list)):
            raise ValueError("t should be an list or tuple")
        for factor in t:
            if not isinstance(factor, (tuple, list)):
                raise ValueError("each factor in t should be a list or tuple")
            if not all(isinstance(x,(int, Integer)) for x in factor):
                raise ValueError("each nonempty factor should contain integers")
            for i in range(len(factor)-1):
                if factor[i] <= factor[i+1]:
                    raise ValueError("each nonempty factor should be a strictly decreasing sequence")

def _check_containment(t, parent):
    """
    Check if ``t`` is an element of ``parent``.

    EXAMPLES::

        sage: S = SymmetricGroup(3+1)
        sage: w = S.from_reduced_word([1, 3, 2])
        sage: B = crystals.FullyCommutativeStableGrothendieck(w, 3, 2)

        sage: from sage.combinat.crystals.fully_commutative_stable_grothendieck import DecreasingHeckeFactorization, _check_containment
        sage: h1 = DecreasingHeckeFactorization([[3, 1], [3], [3, 2]], 4)
        sage: _check_containment(h1, B)

        sage: h2 = DecreasingHeckeFactorization([[3, 1], [3], [], [3, 2]])
        sage: _check_containment(h2, B)
        Traceback (most recent call last):
        ...
        ValueError: number of factors do not match

        sage: h3 = [[3, 1], [2], [3, 2]]
        sage: _check_containment(h3, B)
        Traceback (most recent call last):
        ...
        ValueError: self and parent must be specified based on equivalent words

        sage: h4 = DecreasingHeckeFactorization([[3, 1], [3, 1], [3, 2]], 3)
        sage: _check_containment(h4, B)
        Traceback (most recent call last):
        ...
        ValueError: number of excess letters do not match
    """
    if isinstance(t, DecreasingHeckeFactorization):
        factors = t.factors
        w = t.w
        excess = t.excess
    else:
        factors = len(t)
        max_value = parent.max_value
        from sage.monoids.hecke_monoid import HeckeMonoid
        H = HeckeMonoid(SymmetricGroup(max_value+1))
        w = tuple(H.from_reduced_word(x for factor in t for x in factor).reduced_word())
        excess = sum(len(l) for l in t) - len(w)

    if factors != parent.factors:
        raise ValueError("number of factors do not match")
    if w != parent.w:
            raise ValueError("self and parent must be specified based on equivalent words")
    if excess != parent.excess:
        raise ValueError("number of excess letters do not match")

def _generate_decreasing_hecke_factorizations(w, factors, ex, weight=None, parent=None):
    """
    Generate all decreasing factorizations of word ``w`` in a 0-Hecke monoid
    with fixed excess and number of factors.

    INPUT:

    - ``w`` -- a reduced word, expressed as an iterable

    - ``factors`` -- number of factors for each decreasing factorization

    - ``ex`` -- number of extra letters in each decreasing factorizations

    - ``weight`` -- (default: None) if None, returns all possible decreasing
                    factorizations, otherwise return all those with the
                    specified weight

    EXAMPLES::

        sage: from sage.combinat.crystals.fully_commutative_stable_grothendieck import _generate_decreasing_hecke_factorizations
        sage: _generate_decreasing_hecke_factorizations([1, 2, 1], 3, 1)
        [()(2, 1)(2, 1),
         (2)(1)(2, 1),
         (1)(2)(2, 1),
         (1)(1)(2, 1),
         (2, 1)()(2, 1),
         (2)(2, 1)(2),
         (1)(2, 1)(2),
         (1)(2, 1)(1),
         (2, 1)(2)(2),
         (2, 1)(2)(1),
         (2, 1)(1)(2),
         (2, 1)(2, 1)()]

        sage: _generate_decreasing_hecke_factorizations([1, 2, 1], 3, 1, weight=[1, 1, 2])
        [(2, 1)(2)(2), (2, 1)(2)(1), (2, 1)(1)(2)]
    """
    if parent is None:
        max_value = max(w) if w else 1
        S = SymmetricGroup(max_value+1)
        v = S.from_reduced_word(w)
        parent = DecreasingHeckeFactorizations(v, factors, ex)

    _canonical_word = lambda w, ex: [list(w)[0]]*ex + list(w)
    wt = lambda t:[len(factor) for factor in reversed(t)]

    L = _list_equivalent_words(_canonical_word(w, ex))
    Factors = []
    for word in L:
        F = _list_all_decreasing_runs(word, factors)
        for f in F:
            t = [[word[j] for j in range(len(word)) if f[j]==i] for i in range(factors, 0, -1)]
            if weight is None or weight == wt(t):
                Factors.append(parent.element_class(parent, t))
    return sorted(Factors, reverse=True)

def _list_all_decreasing_runs(word, m):
    """
    List all possible decreasing runs into ``m`` factors for ``word`` in a
    0-Hecke monoid.

    EXAMPLES::

        sage: from sage.combinat.crystals.fully_commutative_stable_grothendieck import _list_all_decreasing_runs
        sage: _list_all_decreasing_runs([2, 1, 2, 1], 3)
        [[2, 2, 1, 1], [3, 2, 1, 1], [3, 3, 1, 1], [3, 3, 2, 1], [3, 3, 2, 2]]
    """
    from sage.combinat.integer_vector import IntegerVectors
    J = _jumps(word)
    jump_vector = [1]+[int(j in J) for j in range(1, len(word))]
    I = sorted(IntegerVectors(m-1-len(J), len(word)+1), reverse=True)
    P = [[elt[i]+jump_vector[i] for i in range(len(word))] for elt in I]
    V = [[m+1-sum(elt[:i+1]) for i in range(len(elt))] for elt in P]
    return V

def _to_reduced_word(P):
    """
    Return a reduced word associated to skew partition ``P``.

    EXAMPLES::

        sage: from sage.combinat.crystals.fully_commutative_stable_grothendieck import _to_reduced_word
        sage: P = SkewPartition([[2, 2], [1]])
        sage: from sage.combinat.crystals.fully_commutative_stable_grothendieck import _to_reduced_word
        sage: _to_reduced_word(P)
        [2, 1, 3]

        sage: P = SkewPartition([[], []])
        sage: from sage.combinat.crystals.fully_commutative_stable_grothendieck import _to_reduced_word
        sage: _to_reduced_word(P)
        []

        sage: P = SkewPartition([[2, 1], []])
        sage: from sage.combinat.crystals.fully_commutative_stable_grothendieck import _to_reduced_word
        sage: _to_reduced_word(P)
        [1, 3, 2]
    """
    cells = P.cells()
    if not cells:
        return []
    m = max(cell[0] for cell in cells) + 1
    n = max(cell[1] for cell in cells) + 1
    L = []
    for i in range(m, -1, -1):
        for j in range(n, -1, -1):
            if (i, j) in cells:
                L += [j-i+m]
    return L

def _lowest_weights(w, factors, ex, parent=None):
    """
    Generate all decreasing factorizations in the 0-Hecke monoid that correspond
    to some valid semistandard Young tableaux.

    The semistandard Young tableaux should have at most ``factors`` columns and their
    column reading words should be equivalent to ``w`` in a 0-Hecke monoid.

    INPUT:

    - ``w`` -- a fully commutative reduced word, expressed as an iterable

    - ``factors`` -- number of factors for each decreasing factorization

    - ``ex`` -- number of extra letters in each decreasing factorizations

    - ``parent`` -- (default: None) parent of the decreasing factorizations, automatically assigned if it is None

    EXAMPLES::

        sage: from sage.combinat.crystals.fully_commutative_stable_grothendieck import _lowest_weights
        sage: _lowest_weights([1, 2, 1], 3, 1)
        Traceback (most recent call last):
        ...
        ValueError: the word w should be fully commutative

        sage: _lowest_weights([2, 1, 3, 2], 4, 3)
        [(2, 1)(3, 1)(3, 1)(2), (2, 1)(3, 1)(3, 2)(2)]

        sage: _lowest_weights([2, 1, 3, 2], 5, 3)
        [(2, 1)(3, 1)(3, 1)(2)(),
         (2, 1)(3, 1)(3, 2)(2)(),
         (2, 1)(3, 1)(1)(1)(2),
         (2, 1)(3, 1)(1)(2)(2),
         (2, 1)(3, 1)(2)(2)(2),
         (2, 1)(3, 2)(2)(2)(2)]

        sage: _lowest_weights([1, 3], 3, 1)
        [(3, 1)(1)(), (3, 1)(3)(), (1)(1)(3), (1)(3)(3)]

        sage: _lowest_weights([3, 2, 1], 5, 2)
        [(3, 2, 1)(1)(1)()()]
    """
    p = permutation.from_reduced_word(w)
    if p.has_pattern([3,2,1]):
        raise ValueError("the word w should be fully commutative")
    if parent is None:
        k = max(w)
        S = SymmetricGroup(k+1)
        word = S.from_reduced_word(w)
        parent = FullyCommutativeStableGrothendieckCrystal(word, factors, ex)

    _canonical_word = lambda w, ex: [list(w)[0]]*ex + list(w)
    L = _list_equivalent_words(_canonical_word(w, ex))

    M = []
    for v in L:
        if _is_valid_column_word(v, factors):
            J = [0] + _jumps(v) + [len(v)]
            t = [v[J[i]:J[i+1]] for i in range(len(J)-1)]
            if len(J) < factors+1:
                t += [()]*(factors+1-len(J))
            M.append(parent.element_class(parent, t))
    return sorted(M)

def _jumps(w):
    """
    Detect all positions where letters weakly increase in ``w``.

    EXAMPLES::

        sage: from sage.combinat.crystals.fully_commutative_stable_grothendieck import _jumps
        sage: w = [4, 1, 2, 1, 4, 3, 2, 1, 3, 2, 2]
        sage: _jumps(w)
        [2, 4, 8, 10]
    """
    return [i+1 for i in range(len(w)-1) if w[i]<=w[i+1]]

def _is_valid_column_word(w, m=None):
    """
    Determine if ``w`` is actually a valid column reading word of some
    semistandard Young tableau with at most ``m`` columns.

    If ``m`` is None, then we determine if ``w`` is a valid column reading word
    of some semistandard Young tableau.

    EXAMPLES::

        sage: from sage.combinat.crystals.fully_commutative_stable_grothendieck import _is_valid_column_word
        sage: w = [3, 2, 2, 1, 1]
        sage: _is_valid_column_word(w)
        False

        sage: w = [3, 2, 1, 1, 1]
        sage: _is_valid_column_word(w,3)
        True

        sage: w = [3, 2, 1, 1, 1]
        sage: _is_valid_column_word(w,2)
        False

        sage: w = [3, 2, 1, 3, 1]
        sage: _is_valid_column_word(w,2)
        True
    """
    J = [0] + _jumps(w) + [len(w)]
    L = [w[J[i+1]-1:J[i]:-1] for i in range(len(J)-1)]
    if all(len(L[i]) >= len(L[i+1]) for i in range(len(L)-1)):
        if m is None or len(_jumps(w)) <= m-1:
            # By construction the sequences along rows of L are strictly
            # decreasing, so it remains to verify that the sequences along
            # columns of L are weakly increasing
            return all(L[i+1][j] >= L[i][j] for i in range(len(L)-1)
                       for j in range(len(L[i+1])))
    return False

def _list_equivalent_words(w):
    """
    List all words equivalent to ``w`` in a 0-Hecke monoid.

    EXAMPLES::

        sage: from sage.combinat.crystals.fully_commutative_stable_grothendieck import _list_equivalent_words
        sage: _list_equivalent_words([1, 1, 2, 1])
        [(1, 1, 2, 1),
         (1, 2, 1, 1),
         (1, 2, 1, 2),
         (1, 2, 2, 1),
         (2, 1, 1, 2),
         (2, 1, 2, 1),
         (2, 1, 2, 2),
         (2, 2, 1, 2)]

        sage: _list_equivalent_words([2,1,3,1,2])
        [(2, 1, 1, 3, 2),
         (2, 1, 3, 1, 2),
         (2, 1, 3, 2, 2),
         (2, 1, 3, 3, 2),
         (2, 2, 1, 3, 2),
         (2, 2, 3, 1, 2),
         (2, 3, 1, 1, 2),
         (2, 3, 1, 2, 2),
         (2, 3, 1, 3, 2),
         (2, 3, 3, 1, 2)]
    """
    if all(isinstance(i, (int, Integer)) for i in w):
        u = w
    else:
        raise ValueError("w needs to be a tuple of integers")

    def _applicable_relations(word):
        """
        Return all positions where a relation can be applied on ``word``
        along with the type of relation.
        """
        L = []
        for i in range(len(word)-2):
            p, q, r = word[i:(i+2)+1]
            if abs(p-q) > 1:
                L += [[i,"pq=qp"]]
            elif abs(p-q) == 1:
                if p == r:  # p != q by the abs test
                    L += [[i,"pqp=qpq"]]
            elif r != p:  # We must have p == q
                L += [[i,"ppq=pqq"]]
            if q == r and r != p:
                L += [[i,"pqq=ppq"]]
        if len(word) > 1 and abs(word[-2]-word[-1]) > 1:
            L += [[len(word)-2,"pq=qp"]]
        return L

    V = set()
    queue = [tuple(u)]
    while queue:
        v = queue.pop(0)
        if tuple(v) not in V:
            V.add(tuple(v))
            L = _applicable_relations(v)
            for pair in L:
                position, move = pair
                t = _apply_relations(v, position, move)
                queue += [tuple(t)]
    return sorted(v for v in list(V))

def _apply_relations(word, position, move):
    """
    Apply a particular type of ``move`` on ``word`` at the specified
    ``position`` using a relation in a 0-Hecke monoid .

    EXAMPLES::

        sage: from sage.combinat.crystals.fully_commutative_stable_grothendieck import _apply_relations
        sage: w = [2, 1, 3, 4]
        sage: _apply_relations(w, position=1, move="pq=qp")
        [2, 3, 1, 4]

        sage: w = [1, 3, 2, 1, 2, 4]
        sage: _apply_relations(w, position=2, move="pqp=qpq")
        [1, 3, 1, 2, 1, 4]

        sage: w = [2, 3, 1, 2, 2, 3]
        sage: _apply_relations(w, position=3, move="pp=p")
        [2, 3, 1, 2, 3]

        sage: w = [2, 3, 1, 2, 3]
        sage: _apply_relations(w, position=3, move="p=pp")
        [2, 3, 1, 2, 2, 3]

        sage: w = [2, 3, 1, 2, 2, 3]
        sage: _apply_relations(w, position=2, move="pqq=ppq")
        [2, 3, 1, 1, 2, 3]

        sage: w = [2, 3, 1, 1, 2, 3]
        sage: _apply_relations(w, position=2, move="ppq=pqq")
        [2, 3, 1, 2, 2, 3]
    """
    w = list(word)
    # Type 1
    if move == "pq=qp":
        p = w[position]
        q = w[position+1]
        w[position] = q
        w[position+1] = p
    # Type 2
    elif move == "pqp=qpq":
        p = w[position]
        q = w[position+1]
        w[position] = q
        w[position+1] = p
        w[position+2] = q
    # Type 3
    elif move == "pqq=ppq":
        p = w[position]
        q = w[position+2]
        w[position+1] = p
    # Type 4
    elif move == "ppq=pqq":
        p = w[position]
        q = w[position+2]
        w[position+1] = q
    # Type 5
    elif move == "pp=p":
        p = w[position]
        w = w[:position+1] + w[position+2:]
    elif move == "p=pp":
        p = w[position]
        w = w[:position+1] + [p] + w[position+1:]
    return w
