# sage.doctest: needs sage.combinat sage.modules
r"""
`\mathcal{B}(\infty)` Crystal Of PBW Monomials

AUTHORS:

- Dinakar Muthiah (2015-05-11): initial version

.. SEEALSO::

    For information on PBW datum, see
    :ref:`sage.combinat.crystals.pbw_datum`.
"""

# ****************************************************************************
#       Copyright (C) 2015 Dinakar Muthiah <muthiah at ualberta.ca>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.richcmp import richcmp
from sage.categories.highest_weight_crystals import HighestWeightCrystals
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.crystals.pbw_datum import PBWData, PBWDatum


class PBWCrystalElement(Element):
    """
    A crystal element in the PBW model.
    """

    def __init__(self, parent, lusztig_datum, long_word=None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: B = crystals.infinity.PBW(['F', 4])
            sage: u = B.highest_weight_vector()
            sage: b = u.f_string([1,2,3,4,2,3,2,3,4,1,2])
            sage: TestSuite(b).run()
        """
        Element.__init__(self, parent)
        if long_word is None:
            long_word = parent._default_word
        self._pbw_datum = PBWDatum(parent._pbw_datum_parent, long_word, lusztig_datum)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: B = crystals.infinity.PBW(['B', 4])
            sage: u = B.highest_weight_vector()
            sage: u.f_string([1,2,3,4,2,3,2,3,4,1,2])
            PBW monomial with Lusztig datum
            (0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 1, 2)
        """
        pbw_datum = self._pbw_datum.convert_to_new_long_word(self.parent()._default_word)
        return "PBW monomial with Lusztig datum {}".format(pbw_datum.lusztig_datum)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: B = crystals.infinity.PBW(['F', 4])
            sage: u = B.highest_weight_vector()
            sage: b = u.f_string([1,2,3,4,2,3,2,3,4,1,2])
            sage: latex(b)
            f_{\alpha_{4}}^{2}
             f_{\alpha_{3}}
             f_{\alpha_{1} + \alpha_{2} + 2 \alpha_{3}}
             f_{\alpha_{1} + \alpha_{2}}
             f_{\alpha_{2}}^{2}
        """
        pbw_datum = self._pbw_datum.convert_to_new_long_word(self.parent()._default_word)
        lusztig_datum = list(pbw_datum.lusztig_datum)
        al = self.parent()._pbw_datum_parent._root_list_from(self.parent()._default_word)
        from sage.misc.latex import latex
        ret_str = ' '.join("f_{%s}%s" % (latex(al[i]), "^{%s}" % latex(exp) if exp > 1 else "")
                           for i, exp in enumerate(lusztig_datum) if exp)
        if ret_str == '':
            return '1'
        return ret_str

    def lusztig_datum(self, word=None):
        """
        Return the Lusztig datum of ``self`` with respect to the reduced
        expression of the long word ``word``.

        EXAMPLES::

            sage: B = crystals.infinity.PBW(['A', 2])
            sage: u = B.highest_weight_vector()
            sage: b = u.f_string([2,1,2,2,2,2,1,1,2,1,2,1,2,1,2,2])
            sage: b.lusztig_datum()
            (6, 0, 10)
            sage: b.lusztig_datum(word=[2,1,2])
            (4, 6, 0)
        """
        if word is None:
            word = self.parent()._default_word
        else:
            self.parent()._check_is_long_word(word)
        word = tuple(word)
        pbw_datum = self._pbw_datum.convert_to_new_long_word(word)
        return tuple(pbw_datum.lusztig_datum)

    def __eq__(self, other):
        """
        Check equality of ``self`` with ``other``.

        EXAMPLES::

            sage: B = crystals.infinity.PBW(['A', 2])
            sage: u = B.highest_weight_vector()
            sage: b = u.f_string([2,1,2,2,2,2,1,1,2,1,2,1,2,1,2,2])
            sage: bp = u.f_string([2,1,2,2,1,1,2,2,2,1,2,1,2,2,1,2])
            sage: b == bp
            True
        """
        if other not in self.parent():
            return False
        other_long_word = other._pbw_datum.long_word
        other_lusztig_datum = other._pbw_datum.lusztig_datum
        equiv_pbw_datum = self._pbw_datum.convert_to_new_long_word(other_long_word)
        return equiv_pbw_datum.lusztig_datum == other_lusztig_datum

    def __ne__(self, other):
        """
        Check inequality of ``self`` with ``other``.

        EXAMPLES::

            sage: B = crystals.infinity.PBW(['A', 2])
            sage: u = B.highest_weight_vector()
            sage: b = u.f_string([2,1,2,2,2,2,1,1,2,1,2,1,2,1,2,2])
            sage: bp = u.f_string([2,1,2,2,1,1,2,2,2,1,2,1,2,2,1,2])
            sage: b != bp
            False
        """
        return not (self == other)

    # Necessary for displaying subcrystals
    def _richcmp_(self, other, op):
        """
        Return comparison of ``self`` and ``other``.

        EXAMPLES::

            sage: B = crystals.infinity.PBW(['A', 2])
            sage: u = B.highest_weight_vector()
            sage: b = u.f_string([2,1,2,2,2,2,1,1,2,1,2,1,2,1,2,2])
            sage: bp = u.f_string([2,1,2,2,1,1,2,2,2,1,2,1,2])
            sage: w = [1, 2, 1]
            sage: (b < bp) == (b.lusztig_datum(w) < bp.lusztig_datum(w))
            True
            sage: (b > bp) == (b.lusztig_datum(w) > bp.lusztig_datum(w))
            True
        """
        i = self.parent().index_set()[0]
        word = self.parent()._pbw_datum_parent._long_word_begin_with(i)
        lusztig_datum = tuple(self._pbw_datum.convert_to_new_long_word(word).lusztig_datum)
        other_lusztig_datum = tuple(other._pbw_datum.convert_to_new_long_word(word).lusztig_datum)
        return richcmp(lusztig_datum, other_lusztig_datum, op)

    @cached_method
    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: B = crystals.infinity.PBW(['A', 2])
            sage: u = B.highest_weight_vector()
            sage: b = u.f_string([2,1,2,2,2,2,1,1,2,1,2,1,2,1,2,2])
            sage: bp = u.f_string([2,1,2,2,1,1,2,2,2,1,2,1,2,2,1,2])
            sage: hash(b) == hash(bp)
            True
        """
        i = self.parent().index_set()[0]
        word = self.parent()._pbw_datum_parent._long_word_begin_with(i)
        pbw_datum = self._pbw_datum.convert_to_new_long_word(word)
        return hash(tuple(pbw_datum.lusztig_datum))

    def e(self, i):
        """
        Return the action of `e_i` on ``self``.

        EXAMPLES::

            sage: B = crystals.infinity.PBW(['B', 3])
            sage: b = B.highest_weight_vector()
            sage: c = b.f_string([2,1,3,2,1,3,2,2]); c
            PBW monomial with Lusztig datum (0, 1, 0, 1, 0, 0, 0, 1, 2)
            sage: c.e(2)
            PBW monomial with Lusztig datum (0, 1, 0, 1, 0, 0, 0, 1, 1)
            sage: c.e_string([2,2,1,3,2,1,3,2]) == b
            True
        """
        equiv_pbw_datum = self._pbw_datum.convert_to_long_word_with_first_letter(i)
        new_long_word = equiv_pbw_datum.long_word
        new_lusztig_datum = list(equiv_pbw_datum.lusztig_datum)
        if new_lusztig_datum[0] == 0:
            return None
        new_lusztig_datum[0] -= 1
        return type(self)(self.parent(), tuple(new_lusztig_datum), new_long_word)

    def f(self, i):
        """
        Return the action of `f_i` on ``self``.

        EXAMPLES::

            sage: B = crystals.infinity.PBW("D4")
            sage: b = B.highest_weight_vector()
            sage: c = b.f_string([1,2,3,1,2,3,4]); c
            PBW monomial with Lusztig datum (0, 1, 0, 0, 0, 0, 0, 2, 0, 0, 2, 0)
            sage: c == b.f_string([1,2,4,1,2,3,3])
            True
        """
        equiv_PBWDatum = self._pbw_datum.convert_to_long_word_with_first_letter(i)
        new_long_word = equiv_PBWDatum.long_word
        new_lusztig_datum = list(equiv_PBWDatum.lusztig_datum)
        new_lusztig_datum[0] += 1
        return type(self)(self.parent(), tuple(new_lusztig_datum), new_long_word)

    def epsilon(self, i):
        r"""
        Return `\varepsilon_i` of ``self``.

        EXAMPLES::

            sage: B = crystals.infinity.PBW(["A2"])
            sage: s = B((3,0,0), (1,2,1))
            sage: s.epsilon(1)
            3
            sage: s.epsilon(2)
            0
        """
        equiv_pbw_datum = self._pbw_datum.convert_to_long_word_with_first_letter(i)
        return equiv_pbw_datum.lusztig_datum[0]

    def phi(self, i):
        r"""
        Return `\varphi_i` of ``self``.

        EXAMPLES::

            sage: B = crystals.infinity.PBW(['A', 2])
            sage: s = B((3,0,0), (1,2,1))
            sage: s.phi(1)
            -3
            sage: s.phi(2)
            3
        """
        WLR = self.parent().weight_lattice_realization()
        h = WLR.simple_coroots()
        return self.epsilon(i) + self.weight().scalar(h[i])

    def weight(self):
        """
        Return weight of ``self``.

        EXAMPLES::

            sage: B = crystals.infinity.PBW(['A', 2])
            sage: s = B((2,2,2), (1,2,1))
            sage: s.weight()
            (-4, 0, 4)
        """
        WLR = self.parent().weight_lattice_realization()
        al = WLR.simple_roots()
        return WLR.sum(c*al[i] for i,c in self._pbw_datum.weight())

    def star(self):
        r"""
        Return the starred crystal element corresponding
        to ``self``.

        Let `b` be an element of ``self`` with Lusztig datum
        `(b_1, \ldots, b_N)` with respect to `w_0 = s_{i_1} \cdots s_{i_N}`.
        Then `b^*` is the element with Lusztig datum `(b_N, \ldots, b_1)`
        with respect to `w_0 = s_{i_N^*} \cdots s_{i_1^*}`, where
        `i_j^* = \omega(i_j)` with `\omega` being the :meth:`automorphism
        <sage.combinat.root_system.cartan_type.CartanType_standard_finite.opposition_automorphism>`
        given by the action of `w_0` on the simple roots.

        EXAMPLES::

            sage: P = crystals.infinity.PBW(['A', 2])
            sage: P((1,2,3), (1,2,1)).star() == P((3,2,1), (2,1,2))
            True

            sage: B = crystals.infinity.PBW(['E', 6])
            sage: b = B.highest_weight_vector()
            sage: c = b.f_string([1,2,6,3,4,2,5,2,3,4,1,6])
            sage: c == c.star().star()
            True

        TESTS::

            sage: from itertools import product
            sage: def test_star(PBW, depth):
            ....:     S = crystals.infinity.Star(PBW)
            ....:     for f_str in product(*([PBW.index_set()]*depth)):
            ....:         x = PBW.highest_weight_vector().f_string(f_str).star()
            ....:         y = S.highest_weight_vector().f_string(f_str)
            ....:         assert x.lusztig_datum() == y.value.lusztig_datum()
            sage: P = crystals.infinity.PBW(['A', 2])
            sage: test_star(P, 5)
            sage: P = crystals.infinity.PBW(['A', 3])
            sage: test_star(P, 5)
            sage: P = crystals.infinity.PBW(['B', 3])
            sage: test_star(P, 5)
            sage: P = crystals.infinity.PBW(['C', 3])
            sage: test_star(P, 5)
            sage: P = crystals.infinity.PBW(['D', 4])
            sage: test_star(P, 5)  # long time
            sage: P = crystals.infinity.PBW(['D', 5])
            sage: test_star(P, 4)  # long time
            sage: P = crystals.infinity.PBW(['E', 6])
            sage: test_star(P, 4)  # long time
            sage: P = crystals.infinity.PBW(['F', 4])
            sage: test_star(P, 4)  # long time
            sage: P = crystals.infinity.PBW(['G', 2])
            sage: test_star(P, 5)
        """
        starred_pbw_datum = self._pbw_datum.star()
        return type(self)(self.parent(), starred_pbw_datum.lusztig_datum,
                             starred_pbw_datum.long_word)


class PBWCrystal(Parent, UniqueRepresentation):
    r"""
    Crystal of `\mathcal{B}(\infty)` given by PBW monomials.

    A model of the crystal `\mathcal{B}(\infty)` whose elements are
    PBW datum up to equivalence by the tropical Plücker relations.
    The crystal structure on Lusztig data `x = (x_1, \ldots, x_m)`
    for the reduced word `s_{i_1} \cdots s_{i_m} = w_0` is given as
    follows. Suppose `i_1 = j`, then `f_j x = (x_1 + 1, x_2, \ldots, x_m)`.
    If `i_1 \neq j`, then we use the tropical Plücker relations to
    change the reduced expression such that `i_1' = j` and then we
    change back to the original word.

    EXAMPLES::

        sage: PBW = crystals.infinity.PBW(['B', 3])
        sage: hw = PBW.highest_weight_vector()
        sage: x = hw.f_string([1,2,2,3,3,1,3,3,2,3,2,1,3,1,2,3,1,2,1,3,2]); x
        PBW monomial with Lusztig datum (1, 1, 1, 3, 1, 0, 0, 1, 1)

    Elements are expressed in terms of Lusztig datum for a fixed
    reduced expression of `w_0`::

        sage: PBW.default_long_word()
        [1, 3, 2, 3, 1, 2, 3, 1, 2]
        sage: PBW.set_default_long_word([2,1,3,2,1,3,2,3,1])
        sage: x
        PBW monomial with Lusztig datum (3, 1, 1, 0, 1, 0, 1, 3, 4)
        sage: PBW.set_default_long_word([1, 3, 2, 3, 1, 2, 3, 1, 2])

    We can construct elements by giving it Lusztig data (with respect
    to the default long word)::

        sage: PBW([1,1,1,3,1,0,0,1,1])
        PBW monomial with Lusztig datum (1, 1, 1, 3, 1, 0, 0, 1, 1)

    We can also construct elements by passing in a reduced expression
    for a long word::

        sage: x = PBW([1,1,1,3,1,0,0,1,1], [3,2,1,3,2,3,2,1,2]); x
        PBW monomial with Lusztig datum (1, 1, 1, 0, 1, 0, 5, 1, 1)
        sage: x.to_highest_weight()[1]
        [1, 2, 2, 2, 2, 2, 1, 3, 3, 3, 3, 2, 3, 2, 3, 3, 2, 3, 3, 2, 1, 3]
    """
    @staticmethod
    def __classcall__(cls, cartan_type):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: B1 = crystals.infinity.PBW(['A', 2])
            sage: B2 = crystals.infinity.PBW("A2")
            sage: B3 = crystals.infinity.PBW(CartanType("A2"))
            sage: B1 is B2 and B2 is B3
            True
        """
        cartan_type = CartanType(cartan_type)
        if not cartan_type.is_finite():
            raise NotImplementedError("only implemented for finite types")
        return super().__classcall__(cls, cartan_type)

    def __init__(self, cartan_type):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: B = crystals.infinity.PBW(['B', 2])
            sage: TestSuite(B).run()
        """
        self._cartan_type = cartan_type
        self._pbw_datum_parent = PBWData(self._cartan_type)
        category = (HighestWeightCrystals(), InfiniteEnumeratedSets())
        Parent.__init__(self, category=category)

        # There must be a better way to do the following
        i = self._cartan_type.index_set()[0]
        self._default_word = self._pbw_datum_parent._long_word_begin_with(i)
        zero_lusztig_datum = [0]*len(self._default_word)
        self.module_generators = (self.element_class(self,
                                                     zero_lusztig_datum,
                                                     self._default_word),)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: crystals.infinity.PBW(['C', 3])
            Crystal of PBW data of type ['C', 3]
        """
        return "Crystal of PBW data of type {}".format(self._cartan_type)

    def default_long_word(self):
        """
        Return the default long word used to express elements of ``self``.

        EXAMPLES::

            sage: B = crystals.infinity.PBW(['E', 6])
            sage: B.default_long_word()
            [1, 3, 4, 5, 6, 2, 4, 5, 3, 4, 1, 3, 2, 4, 5, 6, 2, 4,
             5, 3, 4, 1, 3, 2, 4, 5, 3, 4, 1, 3, 2, 4, 1, 3, 2, 1]
        """
        return list(self._default_word)

    def _check_is_long_word(self, word):
        """
        Check if ``word`` is a reduced expression of the long of the
        Coxeter group of ``self``.

        EXAMPLES::

            sage: B = crystals.infinity.PBW(['A', 3])
            sage: B._check_is_long_word([1,2,1,3,2,1])
            sage: B._check_is_long_word([1,3,2,3,2,1])
            Traceback (most recent call last):
            ...
            ValueError: not a reduced word of the long element
            sage: B._check_is_long_word([1,2,1,3,2])
            Traceback (most recent call last):
            ...
            ValueError: not a reduced word of the long element
            sage: B._check_is_long_word([1,2,1,3,2,1,2])
            Traceback (most recent call last):
            ...
            ValueError: not a reduced word of the long element
        """
        W = self._pbw_datum_parent.weyl_group
        if (len(word) != len(self._default_word)
            or W.from_reduced_word(word) != W.long_element()):
            raise ValueError("not a reduced word of the long element")

    def set_default_long_word(self, word):
        """
        Set the default long word used to express elements of ``self``.

        EXAMPLES::

            sage: B = crystals.infinity.PBW(['C', 3])
            sage: B.default_long_word()
            [1, 3, 2, 3, 1, 2, 3, 1, 2]
            sage: x = B.highest_weight_vector().f_string([2,1,3,2,3,1,2,3,3,1])
            sage: x
            PBW monomial with Lusztig datum (1, 2, 2, 0, 0, 0, 0, 0, 1)
            sage: B.set_default_long_word([2,1,3,2,1,3,2,3,1])
            sage: B.default_long_word()
            [2, 1, 3, 2, 1, 3, 2, 3, 1]
            sage: x
            PBW monomial with Lusztig datum (2, 0, 0, 0, 0, 0, 1, 3, 2)

        TESTS::

            sage: B = crystals.infinity.PBW(['A', 3])
            sage: B._check_is_long_word([1,2,1,3,2,1,2])
            Traceback (most recent call last):
            ...
            ValueError: not a reduced word of the long element
        """
        self._check_is_long_word(word)
        self._default_word = tuple(word)

    Element = PBWCrystalElement
