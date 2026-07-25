# sage.doctest: needs sage.combinat sage.groups
r"""
Kauffman Tangles

This module implements a semigroup used for the monomials of the
Birman-Murakami-Wenzl algebra, which is considered to be realized as a Kauffman
tangle algebra according to [MW2010]_.

AUTHORS:

- Sebastian Oehms Jan 2026: initial version

REFERENCES:

- [MW2010]_
- [EG2017]_, section 6.2.
"""
#############################################################################
#       Copyright (C) 2026 Sebastian Oehms <seb.oehms@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
#############################################################################

from collections.abc import Callable

from sage.combinat.diagram_algebras import BrauerDiagram
from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.plot.graphics import Graphics
from sage.structure.element_wrapper import ElementWrapper
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

##############################################################################
# Tangle element class
##############################################################################


class KauffmanTangle(ElementWrapper):
    r"""
    Element in the semigroup of Kauffman tangles.

    .. SEEALSO::

        :class:`~sage.algebras.birman_murakami_wenzl_algebra.BirmanMurakamiWenzlAlgebra`

    EXAMPLES::

        sage: BMW = algebras.BirmanMurakamiWenzl(3)
        sage: B = BMW.basis().keys()
        sage: d = dict(B)
        sage: BD = list(d.keys())
        sage: bd0 = BD[0]; bd0
        {{-3, 3}, {-2, -1}, {1, 2}}
        sage: bd1 = BD[1]; bd1
        {{-3, 2}, {-2, -1}, {1, 3}}
        sage: t0 = d[bd0]; t0
        e0
        sage: t1 = d[bd1]; t1
        g1*e0
        sage: t01 = t0 * t1; t01
        e0*g1*e0
        sage: t10 = t1 * t0; t10
        g1*e0^2
        sage: bd01 = bd0.compose(bd1); bd01
        ({{-3, 3}, {-2, -1}, {1, 2}}, 0)
        sage: bd10 = bd1.compose(bd0); bd10
        ({{-3, 2}, {-2, -1}, {1, 3}}, 1)
        sage: t01.connector() == bd01
        True
        sage: t10.connector() == bd10
        True
        sage: d[bd10[0]] == t10
        False
        sage: d[bd01[0]] == t01
        False

    .. PLOT::
        :width: 300 px

        BMW = algebras.BirmanMurakamiWenzl(3)
        V = list(dict(BMW.basis().keys()).values())
        sphinx_plot(V[2].plot())
    """
    _crossing_dict = None
    _crossing_info = None
    _connector = None

    def strands(self):
        """
        Return the number of (unclosed) strands.

        OUTPUT: integer

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KT = KauffmanTangles('g0, g1, e0, e1')
            sage: KT((-1, 2, 3, 4)).strands()
            3
        """
        return self.parent().strands()

    @cached_method
    def _shared_memory(self):
        r"""
        Return another, previously defined mutant instance of the class
        of ``self``, which yields the same calculation result for some of
        the methods. These are all calculations that are independent of
        the sign of the braid generators.

        If there is no such other instance, return ``self``.

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KT = KauffmanTangles('g0, g1, e0, e1')
            sage: t1 = KT((-1, 2, 3, 4))
            sage: t1._shared_memory() == t1
            True
            sage: t2 = KT((1, -2, 3, 4))
            sage: t2._shared_memory() == t1
            True
        """
        P = self.parent()
        positive_word = self.positive_word()
        sh_mem = P._shared_memory
        if positive_word not in sh_mem:
            sh_mem[positive_word] = self
        return sh_mem[positive_word]

    def connector(self) -> tuple:
        r"""
        Return the connector of ``self`` as a Brauer diagramm.

        OUTPUT:

        A pair ``(bd, num_removed_loops)`` of an instance of
        :class:`~sage.algebras.diagram_algebras.BrauerDiagram` ``bd`` and an
        integer ``num_removed_loops`` giving the number of closed loops in ``self``.

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KT = KauffmanTangles('g0, g1, e0, e1')
            sage: KT((-1, 2, 3, 4)).connector()
            ({{-3, -2}, {-1, 1}, {2, 3}}, 0)
            sage: KT((-1, 2, -1, 4, 4, 3, 3)).connector()
            ({{-3, 3}, {-2, -1}, {1, 2}}, 2)
        """
        if self._connector:
            return self._connector
        sh_mem = self._shared_memory()
        if sh_mem != self:
            self._connector = sh_mem.connector()
            return self._connector
        P = self.parent()
        n = P._nstrands
        PA = P.BA.ambient()
        con = P.BA.one_basis()
        num_removed_loop = 0
        for ii in self.defining_word():
            i = abs(ii)
            if i < n:
                bd, = PA.s(i).support()
            else:
                i -= (n - 1)
                bd, = PA.a(i).support()
            con, loops = con.compose(bd)
            num_removed_loop += loops
        self._connector = (con, num_removed_loop)
        return self._connector

    @cached_method
    def defining_word(self) -> tuple:
        r"""
        Return the word defining ``self``.

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KT = KauffmanTangles('g0, g1, e0, e1')
            sage: KT((-1, 2, 3, 4)).defining_word()
            (-1, 2, 3, 4)
        """
        return self.value.Tietze()

    @cached_method
    def positive_word(self) -> tuple:
        r"""
        Return a word for a mutant of ``self``, switching negative
        braid generators to positive.

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KT = KauffmanTangles('g0, g1, e0, e1')
            sage: KT((-1, 2, 3, 4)).positive_word()
            (1, 2, 3, 4)
        """
        return tuple([abs(i) for i in self.defining_word()])

    @cached_method
    def positive_mutant(self):
        r"""
        Return a mutant of ``self`` with negative braid generators
        switched to positive.

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KT = KauffmanTangles('g0, g1, e0, e1')
            sage: KT((-1, 2, 3, 4)).positive_mutant()
            g0*g1*e0*e1
        """
        return self.parent()(self.positive_word())

    @cached_method
    def list_of_strands(self) -> list:
        r"""
        Return a list of instances of ``Strand`` covering all strands of
        ``self`` ordered as indicated below.

        It starts with the propagating strands, ordered by their starting position
        in the top line, followed by inline pairs in the top line and then inline
        pairs in the bottom line, each ordered from left to right. Finally, closed
        loops follow, from top to bottom.

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KT = KauffmanTangles('g0, g1, e0, e1')
            sage: el = KT((-1, 2, -1, 4, 4, 3, 3))
            sage: el.list_of_strands()
            [Propagating strand from position 3 on top to position -3 on bottom,
             Inline strand on top line from position 1 to position 2,
             Inline strand on bottom line from position -1 to position -2,
             The 1-th closed loop on the way from top to bottom,
             The 2-th closed loop on the way from top to bottom]
        """
        conn, loops = self.connector()
        res = [Strand(self, i, j) for i, j in conn]
        res += [Strand(self, i, i) for i in range(1, loops + 1)]
        return sorted(res)

    def crossing_dict(self) -> dict:
        r"""
        Return a (possibly empty) dictionary containing as keys those strands
        of ``self`` that intersect another strand.

        The values are lists of pairs `(st, word_position)``, where ``st`` is
        such another strand and ``word_position`` specifies the position
        of the intersection in the defining word of ``self`` as an integer.

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KT = KauffmanTangles('g0, g1, e0, e1')
            sage: el1 = KT((-1, 2, 4))
            sage: el1.crossing_dict()
            {Propagating strand from position 2 on top to position -1 on bottom: [(Inline strand on top line from position 1 to position 3, 0)],
               Inline strand on top line from position 1 to position 3: [(Propagating strand from position 2 on top to position -1 on bottom, 0),
              (Inline strand on top line from position 1 to position 3, 1)]}
            sage: el2 = KT((-1, 2, 4, 3, 3, 4, 3))
            sage: el2.crossing_dict()
            {Propagating strand from position 2 on top to position -3 on bottom: [(Inline strand on top line from position 1 to position 3, 0)],
               Inline strand on top line from position 1 to position 3: [(Propagating strand from position 2 on top to position -3 on bottom, 0),
              (Inline strand on top line from position 1 to position 3, 1)]}
            sage: k3 = KT((4, 1, 1, 1, 4))
            sage: k3.crossing_dict()
            {Propagating strand from position 1 on top to position -1 on bottom: [(Propagating strand from position 1 on top to position -1 on bottom, 1),
              (Propagating strand from position 1 on top to position -1 on bottom, 2),
              (Propagating strand from position 1 on top to position -1 on bottom, 3)]}
        """
        if self._crossing_dict:
            return self._crossing_dict
        sh_mem = self._shared_memory()
        if sh_mem != self:
            self._crossing_dict = sh_mem.crossing_dict()
            return self._crossing_dict
        P = self.parent()
        n = P._nstrands
        word = self.defining_word()
        lw = len(word)
        if not lw:
            return {}
        g = word[-1]
        los = self.list_of_strands()
        if lw == 1:
            if g >= n:
                # no crossings in ``self``
                return {}
            # one crossing in ``self``
            i = abs(g)
            st1 = los[i - 1]
            st2 = los[i]
            return {st1: [(st2, 0)], st2: [(st1, 0)]}

        left_tangle = P(word[:-1])
        gen = P(word[-1:])
        lcrossing_dict = left_tangle.crossing_dict()
        crossing_dict = {}

        def add_crossings_to_dict(st, crossings):
            r"""
            Accumulate crossings per strand in a set and sort once at the end
            """
            crossing_dict.setdefault(st, set()).update(crossings)

        def sorted_crossing_dict():
            r"""
            Return the sorted crossing_dict at the end
            """
            return {st: sorted(cr) for st, cr in crossing_dict.items()}

        # first add crossings from left_tangle
        for lst1 in lcrossing_dict:
            lcrossings = lcrossing_dict[lst1]
            st1 = lst1.expand_in_product(gen)
            crossings = []
            for lst2, lpos in lcrossings:
                st2 = lst2.expand_in_product(gen)
                crossings.append((st2, lpos))
            add_crossings_to_dict(st1, crossings)

        # now add crossing from gen (for gen a braid generator)
        if g < n:
            i = abs(g)
            bot_transpos = [-i, -(i+1)]
            matches = [st for st in los if [st.start, st.end] == bot_transpos]
            if matches:
                # self crossing
                st1, = matches
                add_crossings_to_dict(st1, [(st1, lw - 1)])
                return sorted_crossing_dict()

            matches = [st for st in los if st.start in bot_transpos or st.end in bot_transpos]
            st1, st2 = matches
            # note that st1 < st2
            add_crossings_to_dict(st1, [(st2, lw - 1)])
            add_crossings_to_dict(st2, [(st1, lw - 1)])

        return sorted_crossing_dict()

    @cached_method
    def _strands_at_position(self, pos):
        r"""
        Return a list of strands of ``self`` that start (or end)
        at position ``pos`` in the top line.

        This is a helper function for :meth:`plot`.

        INPUT:

        - ``pos`` -- an integer that indicates a position in the defining
          word of ``self``

        OUTPUT:

        A list of instances of :class:`Strand` ordered according to their start
        (or end) point in the top line.

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KT = KauffmanTangles('g0, g1, e0, e1')
            sage: el = KT((-1, 2, 4, 3, 3, 4, 3))
            sage: el._strands_at_position(0)
            [Inline strand on top line from position 1 to position 3,
             Propagating strand from position 2 on top to position -3 on bottom,
             Inline strand on top line from position 1 to position 3]
        """
        los = self.list_of_strands()
        n = self.strands()
        res = [None] * n
        for i in range(n):
            for st in los:
                p = st.position_sequence()
                if (i + 1, pos) in p:
                    res[i] = st
                    break
        return res

    def crossing_info(self, pos):
        r"""
        Return a dictionary to describe the crossing given by the braid generator
        at position ``pos`` in the defining word of ``self``.

        INPUT:

        - ``pos`` -- integer giving the position in the defining word of the ``self``
          pointing at the braid generator of the crossing

        OUTPUT:

        A dictionary ``{pos_list1: st1, pos_list2: st2}`` of length two. The two
        keys describe the two strands that meet at the crossing as a pair of
        planar coordinates. The values point the the strands as instances of
        :class:`Strand`.

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KT = KauffmanTangles('g0, g1, e0, e1')
            sage: el = KT((-1, 2, 4, 3, 3, 4, 3))
            sage: el.crossing_info(0)
            {((1, 0), (2, 1)): Inline strand on top line from position 1 to position 3,
             ((2, 0),
              (1, 1)): Propagating strand from position 2 on top to position -3 on bottom}
            sage: el.crossing_info(1)
            {((2, 1), (3, 2)): Inline strand on top line from position 1 to position 3,
             ((2, 2), (3, 1)): Inline strand on top line from position 1 to position 3}
        """
        if self._crossing_info is None:
            self._crossing_info = {}
        if pos in self._crossing_info:
            return self._crossing_info[pos]
        sh_mem = self._shared_memory()
        if sh_mem != self:
            self._crossing_info[pos] = sh_mem.crossing_info(pos)
            return self._crossing_info[pos]
        w = self.defining_word()
        i = w[pos]
        if i >= self.strands():
            # not a braid generator
            return {}
        xi = abs(i)
        # search strands
        crs = self.crossing_dict()
        st2 = None
        for st1 in crs:
            for st, p in crs[st1]:
                if p == pos:
                    st2 = st
                    break
            if st2:
                break

        def positions(st):
            ps = [(x, y) for (x, y) in st.position_sequence() if y in (pos, pos + 1) and x in (xi, xi + 1)]
            if ps[0] == ps[-1]:
                # this may hapen in the case of a closed loop
                ps.pop()
            return ps

        if st1 == st2:
            pos_list = positions(st1)
            assert len(set(pos_list)) == 4
            pos_list1 = (pos_list[0], pos_list[1])
            pos_list2 = (pos_list[2], pos_list[3])
        else:
            pos_list1 = tuple(positions(st1))
            pos_list2 = tuple(positions(st2))
        if st1 > st2:
            return {pos_list2: st2, pos_list1: st1}
        return {pos_list1: st1, pos_list2: st2}

    @cached_method
    def find_unlayered_crossing(self) -> int | None:
        r"""
        Return the first postion in the defining word of ``self`` that
        prevents it to be layered according to the order of its strands given
        by ``<<``; if no such position exists, i.e.  ``self`` is layered,
        ``None`` is returned.

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KT = KauffmanTangles('g0, g1, e0, e1')
            sage: KT((-1, 2, -1, 2)).find_unlayered_crossing()
            0
            sage: KT((1, 2, -1, 2)).find_unlayered_crossing()
            2
            sage: KT((1, 2, 1, 2)).find_unlayered_crossing()
            3
        """
        for pos, gen in enumerate(self.defining_word()):
            d = self.crossing_info(pos)
            if d:
                st1, st2 = d.values()
                if st1 == st2:
                    if st1.closure()[st1]:
                        # st1 has reverse orientation in its closure
                        if st1.cross_over(pos, gen):
                            return pos
                    elif not st1.cross_over(pos, gen):
                        return pos
                elif st1 << st2:
                    if st2.cross_over(pos, gen):
                        return pos
                elif st1.cross_over(pos, gen):
                    return pos
        return None

    @cached_method
    def layered_copy(self):
        r"""
        Return a layered copy of ``self`` switching crossing signs if
        necessary; if ``self`` was already layered it is returned itself.

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KT = KauffmanTangles('g0, g1, g2, e0, e1, e2')
            sage: tang = KT((2, 4, 6))
            sage: tang.layered_copy(), tang
            (g1^-1*e0*e2, g1*e0*e2)
        """
        pos = self.find_unlayered_crossing()
        if pos is None:
            return self
        w = list(self.defining_word())
        w[pos] = -w[pos]
        w = tuple(w)
        return self.parent()(w).layered_copy()

    @cached_method
    def writhe(self, closure: bool = False):
        r"""
        Return the writhe of ``self``. This is the sum of loop-signs
        for all self-crossings of strands of ``self``

        INPUT:

        - ``closure`` -- boolean (default ``False``); if ``True``, the
          writhe is calculated with respect to the closure of ``self``

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KT = KauffmanTangles('g0, g1, e0, e1')
            sage: el = KT((-1, 2, 4, 3, 3, 4, 3))
            sage: el.writhe()
            1
            sage: el.writhe(closure=True)
            0
        """
        res = 0
        w = self.defining_word()
        cr_dict = self.crossing_dict()
        for st1 in cr_dict:
            for st2, pos in cr_dict[st1]:
                if st1 == st2:
                    res += st1.crossing_sign(pos, w[pos])
                elif closure:
                    stc = st1.closure()
                    if stc == st2.closure():
                        cs = st1.crossing_sign(pos, w[pos])
                        if stc[st1] == stc[st2]:
                            res += cs
                        else:
                            res -= cs
        return res

    def plot(self, color='rainbow', orientation='top-bottom', gap=0.05,
             aspect_ratio=1, axes=False, **kwds):
        """
        Plot the tangle.

        The following options are available:

        - ``color`` -- (default: ``'rainbow'``) the color of the
          strands; possible values are:

          * ``'rainbow'``, uses :meth:`~sage.plot.colors.rainbow`
            according to the number of strands.

          * a valid color name for :meth:`~sage.plot.bezier_path`
            and :meth:`~sage.plot.line`. Used for all strands.

          * a list or a tuple of colors for each individual strand.

        - ``orientation`` -- (default: ``'top-bottom'``) determines how
          the braid is printed. The possible values are:

          * ``'bottom-top'``, the braid is printed from bottom to top

          * ``'top-bottom'``, the braid is printed from top to bottom

          * ``'left-right'``, the braid is printed from left to right

          Note that the default doesn't matches the default of the orientation
          in :meth:`~sage.groups.braid.Braid.plot` but is according to
          :meth:`~sage.combinat.diagram_algebras.BrauerDiagram.compose`

        - ``gap`` -- floating point number (default: 0.05); determines
          the size of the gap left when a strand goes under another

        - ``aspect_ratio`` -- floating point number (default:
          ``1``); the aspect ratio

        - ``**kwds`` -- other keyword options that are passed to
          :meth:`~sage.plot.bezier_path` and :meth:`~sage.plot.line`

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KT = KauffmanTangles('g0, g1, e0, e1')
            sage: KT((-1, 2, 4, 3, 3, 4, 3)).plot()
            Graphics object consisting of 23 graphics primitives

        .. PLOT::
            :width: 300 px

            from sage.monoids.tangles import KauffmanTangles
            KT = KauffmanTangles('g0, g1, e0, e1')
            el = KT((-1, 2, 4, 3, 3, 4, 3))
            sphinx_plot(el.plot())
            sphinx_plot(el.plot(orientation='bottom-top', axes=True))
        """
        from sage.plot.colors import rainbow
        from sage.plot.plot import Graphics
        n = self.strands()
        los = self.list_of_strands()
        num_strands = len(los)
        if isinstance(color, (list, tuple)):
            if len(color) < num_strands:
                raise TypeError(f"color (={color}) must contain at least {num_strands} colors")
            col = color
        elif color == "rainbow":
            col = rainbow(num_strands)
        else:
            col = [color] * num_strands
        col_st = {los[i]: col[i] for i in range(num_strands)}
        rotation = 0
        if orientation == 'left-right':
            rotation = 1
        elif orientation == 'top-bottom':
            rotation = 2
        elif orientation != 'bottom-top':
            raise ValueError('unknown value for "orientation"')
        word = self.defining_word()
        a = Graphics()
        for i, m in enumerate(word):
            los_i = self._strands_at_position(i)
            for j in range(n):
                stj = los_i[j]
                if j == m - n + 1 and m >= n:
                    continue
                elif j == m - n:
                    los_i_1 = self._strands_at_position(i + 1)
                    stj_1 = los_i_1[j]
                    cap_cup = _CapCupGenPlot(col_st[stj], col_st[stj_1], pos=(j, i), rotation=rotation)
                    a += cap_cup.plot(**kwds)
                elif j == m:
                    continue
                elif j == m - 1 and m < n:
                    stj_1 = los_i[j + 1]
                    cross = _BraidGenPlot(col_st[stj], col_st[stj_1], pos=(j, i), rotation=rotation)
                    a += cross.plot(**kwds)
                elif j == -m:
                    continue
                elif j == -m - 1:
                    stj_1 = los_i[j + 1]
                    cross = _BraidGenPlot(col_st[stj], col_st[stj_1], positive=False, pos=(j, i), rotation=rotation)
                    a += cross.plot(**kwds)
                else:
                    line = _LinePlot(col_st[stj], pos=(j, i), rotation=rotation)
                    a += line.plot(**kwds)
        a.set_aspect_ratio(aspect_ratio)
        a.axes(axes)
        return a


class KauffmanTangles(UniqueRepresentation, Parent):
    r"""
    The semigroup of Kauffman tangles, which naturally index
    a basis for the Birman-Murakami-Wenzl algebra.

    EXAMPLES::

        sage: from sage.monoids.tangles import KauffmanTangles
        sage: KT = KauffmanTangles('g0, g1, e0, e1'); KT
        Semigroup of tangles with 3 (non closed) strands with generators Family (1, g0, g1, e0, e1, g0^-1, g1^-1)
        sage: KT.gens()
        Family (1, g0, g1, e0, e1, g0^-1, g1^-1)

    Element construction::

        sage: from sage.monoids.tangles import KauffmanTangles
        sage: KT = KauffmanTangles('g0, g1, e0, e1')
        sage: KT((-1, 2, 3, 4))
        g0^-1*g1*e0*e1
        sage: KT(KnotInfo.K6_3.braid())
        g0*(g0*g1^-1)^2*g1^-1
        sage: from sage.combinat.diagram_algebras import BrauerDiagram
        sage: bd = BrauerDiagram(((-3, -2), (-1, 2), (1, 3)))
        sage: KT(bd)
        g1*e0*g1^-1*g0^-1
    """

    Element = KauffmanTangle

    def __init__(self, names):
        r"""
        Constructor of ``self``

        TESTS::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KT = KauffmanTangles('g0, g1, e0, e1')
            sage: TestSuite(KT).run()
        """
        from sage.categories.monoids import Monoids
        from sage.groups.free_group import FreeGroup
        FG = FreeGroup(names)
        n = len(FG.gens()) // 2
        # ambient generators of the tangle monoid: the braid generators with
        # their inverses and the cap-cup generators (which have no inverse)
        gens = FG.semigroup_generators()[:-n]
        self._ambient = FG
        self._semigroup_gens = tuple(gens)
        category = Monoids().FinitelyGenerated().Infinite()
        Parent.__init__(self, category=category)
        self._nstrands = n + 1
        self._mwt_names = {}  # support for the names of the BMW-algebra basis
        self._shared_memory = {}  # shared cache for sign independent results of methods

    def _repr_(self):
        r"""
        Return a string representation of ``self``

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KauffmanTangles('B1, B2, C1, C2')
            Semigroup of tangles with 3 (non closed) strands with generators Family (1, B1, B2, C1, C2, B1^-1, B2^-1)
        """
        return 'Semigroup of tangles with %s (non closed) strands with generators %s' % (self.strands(), self.gens())

    @lazy_attribute
    def BA(self):
        """
        Return the Brauer algebra corresponding to ``self``.

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KauffmanTangles('B1, B2, C1, C2').BA
            Brauer Algebra of rank 3 with parameter x over Univariate Polynomial Ring in x over Integer Ring
        """
        from sage.combinat.diagram_algebras import BrauerAlgebra
        from sage.rings.integer_ring import ZZ
        from sage.rings.polynomial.polynomial_ring import polygen
        return BrauerAlgebra(self._nstrands, polygen(ZZ))

    def ambient(self):
        r"""
        Return the ambient free group of ``self``.

        The elements of ``self`` are retracts of elements of this free group.

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KauffmanTangles('g0, g1, e0, e1').ambient()
            Free Group on generators {g0, g1, e0, e1}
        """
        return self._ambient

    def _retract(self, ambient_element) -> KauffmanTangle:
        r"""
        Wrap an element of the ambient free group into ``self``.

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KT = KauffmanTangles('g0, g1, e0, e1')
            sage: FG = KT.ambient()
            sage: KT._retract(FG.gen(0) * FG.gen(2))
            g0*e0
        """
        return self.element_class(self, ambient_element)

    def product(self, x, y) -> KauffmanTangle:
        r"""
        Return the product of two tangles, computed in the ambient free group.

        The result is cached (via :meth:`_cached_product`): the same tangle
        products recur heavily -- in :meth:`KauffmanTangle.expand_in_product`
        every strand of a factor multiplies by the same generator, and the
        Birman-Murakami-Wenzl recursion reuses products -- so memoizing avoids
        recomputing the ambient free group multiplication.  ``product`` itself
        is kept as a plain method because the ``Magmas`` category introspects
        ``self.product.__func__`` during initialization.

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KT = KauffmanTangles('g0, g1, e0, e1')
            sage: KT((1, 3)) * KT((3,))
            g0*e0^2
        """
        return self._cached_product(x, y)

    @cached_method
    def _cached_product(self, x, y) -> KauffmanTangle:
        r"""
        Cached worker for :meth:`product`.

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KT = KauffmanTangles('g0, g1, e0, e1')
            sage: KT._cached_product(KT((1, 3)), KT((3,)))
            g0*e0^2
        """
        return self._retract(x.value * y.value)

    @cached_method
    def gens(self):
        r"""
        Return the generators of ``self`` (including the identity).

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KauffmanTangles('g0, g1, e0, e1').gens()
            Family (1, g0, g1, e0, e1, g0^-1, g1^-1)
        """
        from sage.sets.family import Family
        return Family([self._retract(g) for g in self._semigroup_gens])

    monoid_generators = gens
    semigroup_generators = gens

    def an_element(self) -> KauffmanTangle:
        r"""
        Return an element of ``self``.

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KauffmanTangles('g0, g1, e0, e1').an_element()
            g0
        """
        return self._retract(self._semigroup_gens[1])

    @cached_method
    def one(self):
        r"""
        Return one as element of ``self``.

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KauffmanTangles('g0, g1, e0, e1').one()
            1
        """
        return self._retract(self._ambient.one())

    @cached_method
    def _element_constructor_(self, x):
        r"""
        Return an element of ``self`` constructed from ``x``.

        INPUT:

        - ``x`` -- can be either

          * a ``tuple`` of integers (which are interpreted as generator indices with
            respect to the ambient group such that ``x`` stands for a word in these
            generators)

          * an instance of :class:`~sage.groups.braid.Braid` interpreted according to
            the ``tuple`` with respect to the ``Tietze`` representation.

          * an instance of :class:`~sage.combinat.diagram_algebras.BrauerDiagram`. In
            this case a Morton-Wasserman tangle is constructed using :meth:`morton_wasserman_tangle`.

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KT = KauffmanTangles('g0, g1, e0, e1')
            sage: KT((-1, 2, 3, -2, 4))
            g0^-1*g1*e0*g1^-1*e1

        More examples are given in the docstring of the class.
        """
        A = self.ambient()
        n = self._nstrands - 1
        if isinstance(x, tuple):
            lx = len(x)
            if not lx:
                return self._retract(self._ambient.one())
            if any(abs(i) > n for i in x if i < 0):
                raise ValueError('inverse generators are only for indices <= %s defined' % n)
            elif any(i > 2 * n for i in x):
                raise ValueError('generators are only for indices <= %s defined' % (2 * n))
            return self._retract(A(x))
        if isinstance(x, BrauerDiagram):
            return self.morton_wasserman_tangle(x)
        from sage.groups.braid import Braid
        if isinstance(x, Braid):
            if x.strands() == self._nstrands:
                return self(x.Tietze())
        if isinstance(x, KauffmanTangle):
            if x.parent() is self:
                return x
            x = x.value
        return self._retract(self._ambient(x))

    @cached_method
    def morton_wasserman_tangle(self, bd: BrauerDiagram, top_bottom: bool = True) -> KauffmanTangle:
        r"""
        Return a an element of ``self`` representing the diagram
        as a connector of a simple layered Morton-Wasserman tangle.

        The defining word of the tangle consists of three parts
        ``wt``, ``we`` and ``wb``. In the case of the first and last
        part the tangle generators correspond to braid generators
        (i.e. signed neighbored transpostions, thus the elements
        returned by :meth:`~sage.combinat.diagram_algebras.PartitionAlgebra.s`
        in the ambient algebra of the Brauer algebra with a sign
        attached). The generators for the word in the middle correspond
        to aligned neighbored inline-pairs (i.e. the elements returned
        by :meth:`~sage.combinat.diagram_algebras.PartitionAlgebra.a` in
        the ambient algebra of the Brauer algebra, also called cap-cup
        generators).

        The order of strands according to which the tangle is layered
        is as follows:

        It starts with the propagating strands and goes from left to
        the right on the top line continuing with the strands connecting
        inline pairs on the top line, again from left to the right. Finally
        the strands connecting the inline pairs on the bottom line follow
        from left to the right. The layers proceed from high to low level.

        INPUT:

        - ``bd`` -- :class:`~sage.combinat.diagram_algebras.BrauerDiagramm`

        - ``top_bottom`` -- boolean (default ``True``); multiplication from left to
          right is interpreted from top to bottom in the tangle diagram; to
          reverse the direction, set this keyword argument to ``False``

        OUTPUT:

        A :class:`~sage.monoids.tangles.KauffmanTangle` of ``self``.

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KT = KauffmanTangles('g0, g1, e0, e1')
            sage: BR = algebras.Brauer(3, x)
            sage: BD = sorted(BR.basis().keys())
            sage: d = {bd: KT.morton_wasserman_tangle(bd) for bd in BD}; d
            {{{-3, -2}, {-1, 1}, {2, 3}}: e1,
             {{-3, -2}, {-1, 2}, {1, 3}}: g1*e0*g1^-1*g0^-1,
             {{-3, -2}, {-1, 3}, {1, 2}}: e0*g1^-1*g0^-1,
             {{-3, -1}, {-2, 1}, {2, 3}}: g0*g1*e0*g1^-1,
             {{-3, -1}, {-2, 2}, {1, 3}}: g1*e0*g1^-1,
             {{-3, -1}, {-2, 3}, {1, 2}}: e0*g1^-1,
             {{-3, 1}, {-2, -1}, {2, 3}}: g0*g1*e0,
             {{-3, 1}, {-2, 2}, {-1, 3}}: g1^-1*g0*g1,
             {{-3, 1}, {-2, 3}, {-1, 2}}: g0*g1,
             {{-3, 2}, {-2, -1}, {1, 3}}: g1*e0,
             {{-3, 2}, {-2, 1}, {-1, 3}}: g1*g0,
             {{-3, 2}, {-2, 3}, {-1, 1}}: g1,
             {{-3, 3}, {-2, -1}, {1, 2}}: e0,
             {{-3, 3}, {-2, 1}, {-1, 2}}: g0,
             {{-3, 3}, {-2, 2}, {-1, 1}}: 1}
            sage: all(d[bd].connector() == (bd, 0) for bd in d)
            True
            sage: KT.morton_wasserman_tangle(BD[1], top_bottom=False)
            g0*g1*e0*g1^-1
        """
        # Idea: adjust horizontal pairs on top and bottom by appropriate
        # permutation of top and bottom endpoints in the middle to give an
        # according list `we` of `e`-generators `(e_1, e_3, ... e_{2r-1})`
        # The permutation on the propagating strands are combined with
        # the adjusting top permutation. Then neighbored transpositions
        # in a reduced word for the permutations are replaced by positive
        # or negative braid generators using :meth:`layered_copy`.
        # Note that the restriction of the construction to tangles with
        # one strand less (on the right) gives according results.
        from sage.misc.flatten import flatten

        n = self.strands()
        nb = max(bd.base_set())
        if n != nb:
            raise ValueError('base set of Brauer diagram is incompatible with the number of strands %s' % n)

        num_prop = bd.propagating_number()
        num_e = (n - num_prop) // 2
        t, b, perm = bd.involution_permutation_triple()
        top_fixed = sorted(flatten(t))
        bottom_fixed = sorted(flatten(b))
        one = list(range(1, n + 1))
        top_prop = sorted(i for i in one if i not in top_fixed)
        bottom_prop = sorted(i for i in one if -i not in bottom_fixed)
        middle_prop = list(one)
        top_perm = list(one)
        bottom_perm = list(one)
        prop_perm = list(one)
        we = []
        for m in range(num_e):
            ti, tj = t[m]
            bi, bj = b[num_e - m - 1]
            k = min(ti, -bj)
            if we:
                k = max(k, max(we) + 2)
            we.append(k)
            middle_prop.remove(k)
            middle_prop.remove(k + 1)
            top_perm[ti - 1] = k
            top_perm[tj - 1] = k + 1
            bottom_perm[-bj - 1] = k
            bottom_perm[-bi - 1] = k + 1

        for m in range(num_prop):
            top_perm[top_prop[m] - 1] = middle_prop[m]
            bottom_perm[bottom_prop[m] - 1] = middle_prop[m]
            prop_perm[middle_prop[m] - 1] = middle_prop[perm[m] - 1]

        from sage.combinat.permutation import Permutation as Perm

        p = Perm(prop_perm)
        p_top = Perm(top_perm)
        p_top_p = p_top * p
        p_bottom = Perm(bottom_perm).inverse()
        wt = p_top_p.reduced_word()
        wb = p_bottom.reduced_word()
        if top_bottom:
            wt.reverse()
            wb.reverse()
        else:
            wt, wb = (wb, wt)

        we_shift = [i + n - 1 for i in we]
        tangle = self(tuple(list(wt) + we_shift + list(wb)))
        res = tangle.layered_copy()
        # supporting names of the basis of the BMW algebra
        if res.value.is_one():
            self._mwt_names[bd] = 'o1'
        else:
            self._mwt_names[bd] = str(res)
        return res

    def strands(self):
        """
        Return the number of strands.

        OUTPUT: integer

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KT = KauffmanTangles('g0, g1, e0, e1')
            sage: KT.strands()
            3
        """
        return self._nstrands


class Strand:
    r"""
    Class to deal with strands of a tangle.

    INPUT:

    - ``start`` -- integer, position on the top or bottom line
      (for bottom inline pairs) of the connector
    - ``end`` -- integer, position on the bottom or top line
      (for top inline pairs) of the connector

    For closed loops both values coincide and give the number
    of the loop according to the top to bottom order.

    EXAMPLES::

        sage: from sage.monoids.tangles import KauffmanTangles
        sage: KT = KauffmanTangles('g0, g1, e0, e1')
        sage: el = KT((-1, 2))
        sage: los = el.list_of_strands(); los
        [Propagating strand from position 1 on top to position -3 on bottom,
         Propagating strand from position 2 on top to position -1 on bottom,
         Propagating strand from position 3 on top to position -2 on bottom]
        sage: type(los[0])
        <class 'sage.monoids.tangles.Strand'>
    """
    def __init__(self, tangle: KauffmanTangle, start: int, end: int):
        r"""
        Constructor

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles, Strand
            sage: KT = KauffmanTangles('g0, g1, e0, e1')
            sage: el = KT((-1, 2))
            sage: Strand(el, 0, 1)
            Traceback (most recent call last):
            ...
            ValueError: 0 and 1 do not describe a strand of g0^-1*g1
        """
        positive_tangle = tangle.positive_mutant()
        conn, loops = positive_tangle.connector()
        if ((start, end) not in conn
            and (end, start) not in conn
                and not (start == end and start in range(1, loops + 1))):
            raise ValueError('%s and %s do not describe a strand of %s' % (start, end, tangle))
        self.tangle = positive_tangle
        self.start = start
        self.end = end
        if start < 0 and end > 0:
            # strand is reversed propagating
            self.start = end
            self.end = start
        elif start * end > 0 and abs(start) > abs(end):
            # inline strand is reversed
            self.start = end
            self.end = start

        if self.start > 0 and self.end < 0:
            # strand is propagating
            self.sort = (1, self.start)
        elif self.start > 0 and self.end > 0 and self.start != self.end:
            # inline-strand on top
            self.sort = (2, self.start)
        elif self.start < 0 and self.end < 0:
            # inline-strand on bottom
            self.sort = (3, self.start)
        elif self.start > 0 and self.start == self.end:
            # closed loop
            self.sort = (4, self.start)
        else:
            raise ValueError('no strand constructible')

    def __repr__(self) -> str:
        r"""
        Return representation of ``self`` as string.

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KT = KauffmanTangles('g0, g1, e0, e1')
            sage: el = KT((-1, 2, 4, 3, 3, 4, 3))
            sage: el.list_of_strands()
            [Propagating strand from position 2 on top to position -3 on bottom,
             Inline strand on top line from position 1 to position 3,
             Inline strand on bottom line from position -1 to position -2,
             The 1-th closed loop on the way from top to bottom]
        """
        if self.propagating():
            return 'Propagating strand from position %s on top to position %s on bottom' % (self.start, self.end)
        if self.inline_top():
            return 'Inline strand on top line from position %s to position %s' % (self.start, self.end)
        if self.inline_bottom():
            return 'Inline strand on bottom line from position %s to position %s' % (self.start, self.end)
        return 'The %s-th closed loop on the way from top to bottom' % self.start

    def propagating(self) -> bool:
        r"""
        Return whether ``self`` is a propagating strand.

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KT = KauffmanTangles('g0, g1, e0, e1')
            sage: el = KT((-1, 2, 4, 3, 3, 4, 3))
            sage: [st.propagating() for st in el.list_of_strands()]
            [True, False, False, False]
        """
        return self.sort[0] == 1

    def inline_top(self) -> bool:
        r"""
        Return whether ``self`` starts and ends on the top line.

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KT = KauffmanTangles('g0, g1, e0, e1')
            sage: el = KT((-1, 2, 4, 3, 3, 4, 3))
            sage: [st.inline_top() for st in el.list_of_strands()]
            [False, True, False, False]
        """
        return self.sort[0] == 2

    def inline_bottom(self) -> bool:
        r"""
        Return whether ``self`` starts and ends on the bottom line.

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KT = KauffmanTangles('g0, g1, e0, e1')
            sage: el = KT((-1, 2, 4, 3, 3, 4, 3))
            sage: [st.inline_bottom() for st in el.list_of_strands()]
            [False, False, True, False]
        """
        return self.sort[0] == 3

    def loop(self) -> bool:
        r"""
        Return whether ``self`` is closed.

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KT = KauffmanTangles('g0, g1, e0, e1')
            sage: el = KT((-1, 2, 4, 3, 3, 4, 3))
            sage: [st.loop() for st in el.list_of_strands()]
            [False, False, False, True]
        """
        return self.sort[0] == 4

    def __lt__(self, other) -> bool:
        r"""
        Return ``True`` if ``self`` comes before ``other`` in the default
        order of the strands. Else return ``False``.

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KT = KauffmanTangles('g0, g1, e0, e1')
            sage: el = KT((-1, 2, 4, 3, 3, 4, 3))
            sage: st1, st2, st3, st4 = el.list_of_strands()
            sage: st1 < st2 < st3 < st4
            True
            sage: st1 < st1
            False
        """
        return self.sort < other.sort

    def __eq__(self, other) -> bool:
        r"""
        Return ``True`` if ``self`` and ``other`` describe the same strand.

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KT = KauffmanTangles('g0, g1, e0, e1')
            sage: el = KT((-1, 2, 4, 3, 3, 4, 3))
            sage: st1, st2, st3, st4 = el.list_of_strands()
            sage: st1.neighbor() == st2
            True
            sage: st1.neighbor() == st3
            False
        """
        return hash(self) == hash(other)

    @cached_method
    def __hash__(self) -> int:
        r"""
        Return a hash for ``self``. Note that we share strands of two tangles
        in memory if their words only differ in signs.

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KT = KauffmanTangles('g0, g1, e0, e1')
            sage: el = KT((-1, 2, 4, 3, 3, 4, 3))
            sage: st1, st2, st3, st4 = el.list_of_strands()
            sage: hash(st1.neighbor()) == hash(st2)
            True
        """
        return hash((self.tangle.defining_word(), self.start, self.end))

    @cached_method
    def __lshift__(self, other) -> bool:
        r"""
        Return ``True`` if ``self`` comes before ``other`` in the *closure*
        order of the strands.

        More precisely ``self`` comes before ``other`` if they are parts of the same strand in the closure
        of the tangle and ``self`` comes before ``other`` in the default
        order (i.e. ``self < other``) or if they belong to different strands
        in the closure and the one of ``self`` has a smaler first item.

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KT = KauffmanTangles('g0, g1, e0, e1')
            sage: el = KT((-1,))
            sage: st1, st2, st3 = el.list_of_strands()
            sage: st1 << st2
            True
            sage: st3 << st2
            False
        """
        scl = self.closure()
        ocl = other.closure()
        scll = list(scl)
        if scl == ocl:
            return scll.index(self) < scll.index(other)
        ocll = list(ocl)
        return scll[0] < ocll[0]

    def overlap(self, other):
        r"""
        Return the number of top line positions of ``other``
        that match bottom line positions of ``self``.

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KT = KauffmanTangles('g0, g1, e0, e1')
            sage: el = KT((-1,))
            sage: st1, st2, st3 = el.list_of_strands()
            sage: st1.overlap(st2)
            1
            sage: st1.overlap(st3)
            0
        """
        if self.inline_top() or self.loop():
            return 0
        if other.inline_bottom() or other.loop():
            return 0
        res = 0
        if self.end < 0:
            if self.end == - other.end:
                res = 1
            elif self.end == - other.start:
                res = 1
        if self.start < 0:
            # inline bottom
            if self.start == - other.end:
                res += 1
            elif self.start == - other.start:
                res += 1
        return res

    @cached_method
    def expand_in_product(self, right_tangle: KauffmanTangle):
        r"""
        Return the strand obtained from ``self`` after ``right_tangle`` is
        multiplied on the tangle of self on the right.

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KT = KauffmanTangles('g0, g1, e0, e1')
            sage: el = KT((-1, 2, 4, 3, 3, 4, 3))
            sage: strands = el.list_of_strands(); strands
            [Propagating strand from position 2 on top to position -3 on bottom,
             Inline strand on top line from position 1 to position 3,
             Inline strand on bottom line from position -1 to position -2,
             The 1-th closed loop on the way from top to bottom]
            sage: for st in strands:
            ....:     st.expand_in_product(KT((-1,)))
            Propagating strand from position 2 on top to position -3 on bottom
            Inline strand on top line from position 1 to position 3
            Inline strand on bottom line from position -1 to position -2
            The 1-th closed loop on the way from top to bottom
            sage: for st in strands:
            ....:     st.expand_in_product(KT((2,)))
            Propagating strand from position 2 on top to position -2 on bottom
            Inline strand on top line from position 1 to position 3
            Inline strand on bottom line from position -1 to position -3
            The 1-th closed loop on the way from top to bottom
            sage: for st in strands:
            ....:     st.expand_in_product(KT((3,)))
            Propagating strand from position 2 on top to position -3 on bottom
            Inline strand on top line from position 1 to position 3
            The 2-th closed loop on the way from top to bottom
            The 1-th closed loop on the way from top to bottom
            sage: for st in strands:
            ....:     st.expand_in_product(KT((4,)))
            Propagating strand from position 2 on top to position -1 on bottom
            Inline strand on top line from position 1 to position 3
            Propagating strand from position 2 on top to position -1 on bottom
            The 1-th closed loop on the way from top to bottom
        """
        right_positive_tangle = right_tangle.positive_mutant()
        pr = self.tangle * right_positive_tangle
        if self.inline_top() or self.loop():
            return Strand(pr, self.start, self.end)

        strands_r = right_positive_tangle.list_of_strands()
        matches = [st for st in strands_r if self.overlap(st) == 2]
        if matches:
            # match of two inline strands gives a new loop
            strands_l = self.tangle.list_of_strands()
            last_st = strands_l[-1]
            new_loop = 1
            if last_st.loop():
                new_loop = last_st.start + 1
            return Strand(pr, new_loop, new_loop)

        matches = [st for st in strands_r if self.overlap(st) == 1]

        if not matches:
            return Strand(pr, self.start, self.end)

        def final_pos(pos):
            r"""
            Return the start position of the strand of ``self`` that ends in
            ``pos``.
            """
            if pos < 0:
                return pos
            strands_l = self.tangle.list_of_strands()
            for st_l in strands_l:
                if st_l.end == -pos:
                    return st_l.start
                if st_l.start == -pos:
                    return st_l.end

        def free_pos(match):
            r"""
            Return the position of the strand ``match`` which is not connected
            to ``self``.
            """
            if match.propagating():
                return match.end
            free = match.end
            if -free in (self.start, self.end):
                free = match.start
            return free

        if len(matches) == 1:
            # if there is only one match self must be propagating
            st = matches[0]
            free_st = free_pos(st)
            return Strand(pr, self.start, final_pos(free_st))
        # two matches are only possible if self is inline at bottom
        m1 = matches[0]
        m2 = matches[1]
        free_m1 = free_pos(m1)
        free_m2 = free_pos(m2)
        return Strand(pr, final_pos(free_m1), final_pos(free_m2))

    @cached_method
    def position_sequence(self):
        r"""
        Return the list of positions ``(x, y)`` along the way of ``self``
        from start to end.

        Here ``x in [1, ..., n]`` is the horizontal position according to
        the connector frame where ``n`` is the number of strands. ``y``
        indicates the position of a generator in the word of the tangle
        of ``self``. ``y = 0`` is on top of the first generator and ``y = 1``
        on the bottom of the first generator.

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KT = KauffmanTangles('g0, g1, e0, e1')
            sage: tang1 = KT((-1, 2, 4, 3, 3, 4, 3))
            sage: los1 = tang1.list_of_strands()
            sage: los1[0].position_sequence()
            [(2, 0), (1, 1), (1, 2), (1, 3), (2, 3), (3, 3), (3, 4),
             (3, 5), (2, 5), (1, 5), (1, 6), (2, 6), (3, 6), (3, 7)]
            sage: los1[1].position_sequence()
            [(1, 0), (2, 1), (3, 2), (2, 2), (3, 1), (3, 0)]
            sage: los1[2].position_sequence()
            [(1, 7), (2, 7)]
            sage: los1[3].position_sequence()
            [(1, 4), (2, 4), (1, 4)]
            sage: tang2 = KT((1, 4, 2, 3, -2, -1))
            sage: los2 = tang2.list_of_strands(); los2
            [Propagating strand from position 2 on top to position -1 on bottom,
             Inline strand on top line from position 1 to position 3,
             Inline strand on bottom line from position -2 to position -3]
            sage: los2[0].position_sequence()
            [(2, 0), (1, 1), (1, 2), (1, 3), (2, 3), (3, 2), (2, 2),
             (3, 3), (3, 4), (2, 5), (1, 6)]
            sage: tang3 = KT((3, 1, 3))
            sage: los3 = tang3.list_of_strands(); los3
            [Propagating strand from position 3 on top to position -3 on bottom,
             Inline strand on top line from position 1 to position 2,
             Inline strand on bottom line from position -1 to position -2,
             The 1-th closed loop on the way from top to bottom]
            sage: los3[3].position_sequence()
            [(1, 2), (2, 1), (1, 1), (2, 2), (1, 2)]
        """
        T = self.tangle
        w = T.defining_word()
        n = T.strands()
        lw = len(w)
        P = T.parent()

        if lw == 1:
            if self.start < 0:
                start_pos = (-self.start, lw)
                end_pos = (-self.end, lw)
            else:
                start_pos = (self.start, 0)
                if self.end < 0:
                    end_pos = (-self.end, lw)
                else:
                    end_pos = (self.end, 0)
            return [start_pos, end_pos]

        left_tangle = P(w[:-1])
        i = w[-1]
        gen = P((abs(i),))
        e_gen = False
        if i >= n:
            i = i - n + 1
            e_gen = True
        else:
            i = abs(i)
        gen_pair = (i, i + 1)

        def final_pos(x):
            return [j for j in gen_pair if j != x][0] if x in gen_pair else x

        def find_join_strand(x):
            assert x in gen_pair
            join = None
            for j in gen_pair:
                if j != x:
                    join = -j
                    break
            assert join is not None
            for lst in left_tangle.list_of_strands():
                if join == lst.start or join == lst.end:
                    return lst
            assert False

        def add_bottom(positions):
            xs, ys = positions[0]
            xe, ye = positions[-1]
            spos = []
            epos = []

            xs = final_pos(xs)
            xe = final_pos(xe)

            start = -self.start
            end = -self.end

            revert_orientation = False
            if start != end:  # not a closed loop
                if ys == lw - 1:
                    if xs == start:
                        spos = [(xs, lw)]
                    elif xs == end:
                        revert_orientation = True
                        epos = [(xs, lw)]
                if ye == lw - 1:
                    if xe == start:
                        revert_orientation = True
                        spos = [(xe, lw)]
                    elif xe == end:
                        epos = [(xe, lw)]

            if revert_orientation:
                positions = list(positions)
                positions.reverse()
            return spos + positions + epos

        def join_positions(st1, st2):
            start_positions = st1.position_sequence()
            end_positions = st2.position_sequence()
            if st1 == st2 and self.loop():
                # just append the first position at the
                return add_bottom(start_positions + [end_positions[0]])
            if -st2.end in gen_pair:
                # revert orientation
                # a copy is needed because the result is cached
                end_positions = list(end_positions)
                end_positions.reverse()
            if -st1.start in gen_pair:
                start_positions = list(start_positions)
                start_positions.reverse()
            return add_bottom(start_positions + end_positions)

        lstrands = [lst for lst in left_tangle.list_of_strands() if lst.expand_in_product(gen) == self]

        if not lstrands:
            assert self.inline_bottom()
            return [(-self.start, lw), (-self.end, lw)]
        if len(lstrands) > 1:
            assert e_gen
            lst1 = lstrands[0]
            lst2 = lstrands[1]
            if lst1.start > 0:
                return join_positions(lst1, lst2)
            if lst2.start > 0:
                return join_positions(lst2, lst1)
            # both are bottom-inline
            if lst1.start > lst2.start:
                return join_positions(lst1, lst2)
            return join_positions(lst2, lst1)
        lst1 = lstrands[0]
        positions = lst1.position_sequence()
        xs, ys = positions[0]
        xe, ye = positions[-1]
        if ye < lw - 1 and ys < lw - 1:
            # the way of this strand is not affected by last generator
            return positions
        if xe in gen_pair:
            if e_gen:
                lst2 = find_join_strand(xe)
                return join_positions(lst1, lst2)
            return add_bottom(positions)
        if xe != xs and xs in gen_pair and ys == lw - 1:
            if e_gen:
                lst2 = find_join_strand(xs)
                return join_positions(lst2, lst1)
            return add_bottom(positions)
        return add_bottom(positions)

    @cached_method
    def cross_over(self, pos: int, gen: int) -> bool:
        r"""
        Return ``True`` if ``self`` over-crosses the other strand at the
        crossing determined by the braid generator at the position ``pos``
        in the defining word of the tangle of ``self``. Return ``None``
        if ``self`` doesn't cross at ``pos``.

        INPUT:

        - ``pos`` -- integer pointing at the position of the crossing in
          the defining word of the tangle of ``self``
        - ``gen`` -- integer, the index of the braid generator; note that
          this may have a different sign as the braid generator at ``pos``
          of the tangle of ``self`` because of usage of ``shared_memory``

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KT = KauffmanTangles('g0, g1, e0, e1')
            sage: tang = KT((-1, 2, 4, 3, 3, 4, 3))
            sage: st = tang.list_of_strands()[0]
            sage: st.cross_over(0, -1)
            True
            sage: st.cross_over(1, 2) is None
            True
        """
        T = self.tangle
        d = T.crossing_info(pos)
        if self not in d.values():
            return None
        pos_list = [p for p in d if d[p] == self][0]
        pos_from, pos_to = pos_list
        fx, fy = pos_from
        tx, ty = pos_to
        if fx < tx and fy < ty:
            # from top left to bottom right
            return gen > 0
        if fx > tx and fy > ty:
            # from bottom right to top left
            return gen > 0
        return gen < 0

    @cached_method
    def crossing_sign(self, pos: int, gen: int) -> int:
        r"""
        Return the sign of the self-crossing given by the braid generator
        at position ``pos`` in the defining word of ``self``.

        INPUT:

        - ``pos`` -- integer pointing at the position of the crossing in
          the defining word of the tangle of ``self``
        - ``gen`` -- integer, the index of the braid generator; note that
          this may have a different sign as the braid generator at ``pos``
          of the tangle of ``self`` because of usage of ``shared_memory``

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KT = KauffmanTangles('g0, g1, e0, e1')
            sage: w = (-1, 2, 4, 3, 3, 4, 3)
            sage: tang = KT(w)
            sage: st = tang.list_of_strands()[0]
            sage: [st.crossing_sign(pos, w[pos]) for pos in range(3)]
            [1, 0, 0]
        """
        d = self.tangle.crossing_info(pos)
        if not d:
            return 0
        pos_list1, pos_list2 = d.keys()
        if not self == d[pos_list1]:
            return 0

        if self.cross_over(pos, gen):
            pos_list_o = pos_list1
            pos_list_u = pos_list2
        else:
            pos_list_o = pos_list2
            pos_list_u = pos_list1

        over_in = pos_list_o[0]
        under_in = pos_list_u[0]
        ox, oy = over_in
        ux, uy = under_in
        xi = abs(gen)
        if ox == ux:
            if ox == xi:
                # both come from the left
                if oy < uy:
                    # over from above -> over goes into a right curve
                    return 1
                # over from below -> over goes into a left curve
                return -1
            # both come from the right
            if oy < uy:
                # over from above -> over goes into a left curve
                return -1
            # over from below -> over goes into a right curve
            return 1
        if oy == uy:
            if oy == pos:
                # both come from above
                if ox < ux:
                    # over from the left -> over goes into a left curve
                    return -1
                # over from the right -> over goes into a right curve
                return 1
            # both come from the below
            if ox < ux:
                # over from the left -> over goes into a right curve
                return 1
            # over from the right -> over goes into a left curve
            return -1

    @cached_method
    def neighbor(self, successor: bool = True):
        r"""
        Return the strand next to ``self`` in the extension of ``self``
        in the closure of the tangle.

        Depending on ``successor`` this is the following or the
        preceeding one.

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KT = KauffmanTangles('g0, g1, e0, e1')
            sage: los = KT((-1, 2, 4, 3, 3, 4, 3)).list_of_strands()
            sage: st1, st2, st3, st4 = los
            sage: st1.neighbor() == st2; st3.neighbor() == st1
            True
            True
            sage: st1.neighbor(successor=False) == st3
            True
            sage: st4.neighbor() == st4
            True
        """
        if self.loop():
            return self
        T = self.tangle
        los = T.list_of_strands()
        if successor:
            pos = -self.end
        else:
            pos = -self.start
        matches = [st for st in los if pos in (st.start, st.end)]
        return [st for st in matches if not st.loop()][0]

    @cached_method
    def closure(self) -> dict:
        r"""
        Return a dictionary ``{strand: reverse}`` where ``strand`` is
        another strand of the tangle of ``self`` which is a part of
        the extension of ``self`` in the closure of the tangle and
        ``reverse`` a boolean indicating whether ``strand`` is oriented
        backwards in the closure of ``self``.

        EXAMPLES::

            sage: from sage.monoids.tangles import KauffmanTangles
            sage: KT = KauffmanTangles('g0, g1, e0, e1')
            sage: los = KT((-1, 2, 4, 3, 3, 4, 3)).list_of_strands()
            sage: st1, st2, st3, st4 = los
            sage: st1.closure()
            {Propagating strand from position 2 on top to position -3 on bottom: False,
             Inline strand on top line from position 1 to position 3: True,
             Inline strand on bottom line from position -1 to position -2: False}
            sage: st2.closure() == st1.closure() == st3.closure()
            True
            sage: st4.closure()
            {The 1-th closed loop on the way from top to bottom: False}
        """
        T = self.tangle
        los = T.list_of_strands()
        for other in los:
            if not other < self:
                break
            closure = other.closure()
            if self in closure:
                return closure
        reverse = False
        prec = self
        res = {prec: reverse}
        suc = self.neighbor()
        while self != suc:
            if reverse:
                if prec.start == -suc.start:
                    reverse = False
            elif prec.end == -suc.end:
                reverse = True
            res[suc] = reverse
            prec = suc
            suc = prec.neighbor(successor=not reverse)
        return res


#############################################################################
# Helper classes for plotting
#############################################################################
class _GeneratorPlot:
    r"""
    Abstract class to keep the data for plotting a tangle generator.

    INPUT:

    - ``pos`` -- pair of integers to fix the position in `\ZZ^2`
      lattice
    - ``rotation`` -- integer, giving the rotation as a factor of `\pi/2`

    EXAMPLES::

        sage: from sage.monoids.tangles import _BraidGenPlot
        sage: from sage.plot.colors import rainbow
        sage: col = rainbow(2)
        sage: cross = _BraidGenPlot(col[0], col[1])
        sage: cross._curves
        {(((0, 0), (0.0, 0.25), (0.5, 0.5)), ((1.0, 0.75), (1.0, 1.0))): '#ff0000',
         (((0.4, 0.55), (0.3, 0.6), (0.0, 0.75), (0.0, 1.0)),): '#00ffff',
         (((1, 0), (1.0, 0.25), (0.7, 0.4), (0.6, 0.45)),): '#00ffff'}
    """
    _curves = {}

    def __init__(self, pos=(0, 0), rotation=0):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.monoids.tangles import _GeneratorPlot
            sage: G = _GeneratorPlot()
            sage: G._curves
            {}
        """
        x, y = pos
        if rotation == 1:
            self.rotate(1)
            self.shift(y, -x)
        elif rotation == 2:
            self.rotate(2)
            self.reflect()
            self.shift(x + 1, -y)
        elif pos != (0, 0):
            self.shift(x, y)

    def _move_points(self, move: Callable[[int, int], tuple]) -> None:
        r"""
        Move all curves according to the given (euclidean) ``move`` function.

        EXAMPLES::

            sage: from sage.monoids.tangles import _BraidGenPlot
            sage: from sage.plot.colors import rainbow
            sage: col = rainbow(2)
            sage: cross = _BraidGenPlot(col[0], col[1])
            sage: cross._move_points(lambda a, b: (a + 3, b - 2))
            sage: cross._curves
            {(((3, -2), (3.0, -1.75), (3.5, -1.5)),
              ((4.0, -1.25), (4.0, -1.0))): '#ff0000',
             (((3.4, -1.45), (3.3, -1.4), (3.0, -1.25), (3.0, -1.0)),): '#00ffff',
             (((4, -2), (4.0, -1.75), (3.7, -1.6), (3.6, -1.55)),): '#00ffff'}
        """
        res = {}
        cv = self._curves
        for c in cv:
            col = cv[c]
            new_seg = []
            for seg in c:
                if type(seg[0]) is tuple:
                    nseg = tuple([move(a, b) for a, b in seg])
                else:
                    a, b = seg
                    nseg = move(a, b)
                new_seg.append(nseg)
            res[tuple(new_seg)] = col
        self._curves = res

    def shift(self, x, y):
        r"""
        Return the generator shifted in ``x`` and ``y`` direction.

        EXAMPLES::

            sage: from sage.monoids.tangles import _CapCupGenPlot
            sage: from sage.plot.colors import rainbow
            sage: col = rainbow(2)
            sage: cap_cup = _CapCupGenPlot(col[0], col[1])
            sage: cap_cup.shift(3, -2)
            sage: cap_cup._curves
            {(((3.0, -2.0), (3.0, -1.6), (3.4, -1.6), (3.5, -1.6)),
              ((3.6, -1.6), (4.0, -1.6), (4.0, -2.0))): '#ff0000',
             (((3.0, -1.0), (3.0, -1.4), (3.4, -1.4), (3.5, -1.4)),
              ((3.6, -1.4), (4.0, -1.4), (4.0, -1.0))): '#00ffff'}
        """
        self._move_points(lambda a, b: (a + x, b + y))

    def rotate(self, half_pi_times):
        r"""
        Return the generator rotated clockwise.

        EXAMPLES::

            sage: from sage.monoids.tangles import _CapCupGenPlot
            sage: from sage.plot.colors import rainbow
            sage: col = rainbow(2)
            sage: cap_cup = _CapCupGenPlot(col[0], col[1])
            sage: cap_cup.rotate(1)
            sage: cap_cup._curves
            {(((0.0, -0.0), (0.4, -0.0), (0.4, -0.4), (0.4, -0.5)),
              ((0.4, -0.6), (0.4, -1.0), (0.0, -1.0))): '#ff0000',
             (((1.0, -0.0), (0.6, -0.0), (0.6, -0.4), (0.6, -0.5)),
              ((0.6, -0.6), (0.6, -1.0), (1.0, -1.0))): '#00ffff'}
        """
        if half_pi_times == 1:
            self._move_points(lambda a, b: (b, -a))
        if half_pi_times == 2:
            self._move_points(lambda a, b: (-a, -b))
        if half_pi_times == 3:
            self._move_points(lambda a, b: (-b, a))

    def reflect(self, x_axis=True):
        r"""
        Return the generator reflected at one of the axises.

        INPUT:

        - ``x_axis`` -- boolean, default ``True``. if set to ``False`` the
          reflextion is at the y-axis.

        EXAMPLES::

            sage: from sage.monoids.tangles import _BraidGenPlot
            sage: from sage.plot.colors import rainbow
            sage: col = rainbow(2)
            sage: cross = _BraidGenPlot(col[0], col[1])
            sage: cross.reflect()
            sage: cross._curves
            {(((-1, 0), (-1.0, 0.25), (-0.7, 0.4), (-0.6, 0.45)),): '#00ffff',
             (((-0.4, 0.55), (-0.3, 0.6), (-0.0, 0.75), (-0.0, 1.0)),): '#00ffff',
             (((0, 0), (-0.0, 0.25), (-0.5, 0.5)), ((-1.0, 0.75), (-1.0, 1.0))): '#ff0000'}
        """
        if x_axis:
            self._move_points(lambda a, b: (-a, b))
        else:
            self._move_points(lambda a, b: (a, -b))

    def plot(self, a: Graphics | None = None, **kwds) -> Graphics:
        r"""
        Return a graphics object in which ``self`` is appended to the
        given object ``a``.

        EXAMPLES::

            sage: from sage.monoids.tangles import _CapCupGenPlot, _BraidGenPlot
            sage: from sage.plot.colors import rainbow
            sage: col = rainbow(2)
            sage: cap_cup = _CapCupGenPlot(col[0], col[1])
            sage: cap_cup.plot()
            Graphics object consisting of 2 graphics primitives
            sage: cross = _BraidGenPlot(col[0], col[1])
            sage: cross.plot()
            Graphics object consisting of 3 graphics primitives
        """
        if isinstance(self, _LinePlot):
            from sage.plot.plot import line as plot
        else:
            from sage.plot.bezier_path import bezier_path as plot
        if not a:
            a = Graphics()
        for curves in self._curves:
            a += plot(curves, color=self._curves[curves], **kwds)
        return a


class _CapCupGenPlot(_GeneratorPlot):
    r"""
    Class for plotting a cap-cup generator with given colours ``col_top`` and
    ``col_bot`` for the both strands on top and bottom.

    INPUT:

    - ``col_top`` -- a valid color name for :meth:`~sage.plot.bezier_path`
      and :meth:`~sage.plot.line` used for the top line strand
    - ``col_bottom`` -- a valid color name for :meth:`~sage.plot.bezier_path`
      and :meth:`~sage.plot.line` used for the bottom line strand
    - ``pos`` -- see :class:`GeneratorPlot`
    - ``rotation`` -- integer, giving the rotation as a factor of ``pi/2``
    """
    def __init__(self, col_top, col_bot, pos=(0, 0), rotation=0):
        r"""
        Set the coordinate data for plotting ``self``.

        EXAMPLES::

            sage: from sage.monoids.tangles import _CapCupGenPlot
            sage: from sage.plot.colors import rainbow
            sage: col = rainbow(2)
            sage: cross = _CapCupGenPlot(col[0], col[1], pos=(1, 1), rotation=1)
            sage: cross._curves
            {(((1.0, -1.0), (1.4, -1.0), (1.4, -1.4), (1.4, -1.5)),
            ((1.4, -1.6), (1.4, -2.0), (1.0, -2.0))): '#ff0000',
            (((2.0, -1.0), (1.6, -1.0), (1.6, -1.4), (1.6, -1.5)),
            ((1.6, -1.6), (1.6, -2.0), (2.0, -2.0))): '#00ffff'}
        """
        # bottom curve left
        s1 = (0.0, 0.0)
        c1 = (0.0, 0.4)
        c2 = (0.4, 0.4)
        e1 = (0.5, 0.4)
        cvtl = (s1, c1, c2, e1)

        # bottom curve right
        c3 = (0.6, 0.4)
        c4 = (1.0, 0.4)
        e2 = (1.0, 0.0)
        cvtr = (c3, c4, e2)

        # top curve left
        t1 = (0.0, 1.0)
        d1 = (0.0, 0.6)
        d2 = (0.4, 0.6)
        f1 = (0.5, 0.6)
        cvbl = (t1, d1, d2, f1)

        # top curve right
        d3 = (0.6, 0.6)
        d4 = (1.0, 0.6)
        f2 = (1.0, 1.0)
        cvbr = (d3, d4, f2)

        self._curves = {(cvtl, cvtr): col_top, (cvbl, cvbr): col_bot}
        super().__init__(pos=pos, rotation=rotation)


class _BraidGenPlot(_GeneratorPlot):
    r"""
    Class for plotting a braid generator with given colours ``col_over`` and
    ``col_under`` for the both strands crossing each other.

    INPUT:

    - ``col_over`` -- a valid color name for :meth:`~sage.plot.bezier_path`
      and :meth:`~sage.plot.line` used for the over-crossing strand
    - ``col_under`` -- a valid color name for :meth:`~sage.plot.bezier_path`
      and :meth:`~sage.plot.line` used for the under-crosing strand
    - ``positive`` -- boolean whether to plot a positive or negative crossing
    - ``gap`` -- floating point number (default: 0.05); see the description
      in :meth:`plot`
    - ``pos`` -- see :class:`_GeneratorPlot`
    - ``rotation`` -- integer, giving the rotation as a factor of ``pi/2``
    """
    def __init__(self, col_over, col_under, positive=True, gap=0.05, pos=(0, 0), rotation=0):
        r"""
        Set the coordinate data for plotting ``self``.

        EXAMPLES::

            sage: from sage.monoids.tangles import _BraidGenPlot
            sage: from sage.plot.colors import rainbow
            sage: col = rainbow(2)
            sage: cross = _BraidGenPlot(col[0], col[1], pos=(1, 1), rotation=1)
            sage: cross._curves
            {(((1, -2), (1.25, -2.0), (1.4, -1.7), (1.45, -1.6)),): '#00ffff',
            (((1, -1), (1.25, -1.0), (1.5, -1.5)),
            ((1.75, -2.0), (2.0, -2.0))): '#ff0000',
            (((1.55, -1.4), (1.6, -1.3), (1.75, -1.0), (2.0, -1.0)),): '#00ffff'}
        """
        if not positive:
            col_over, col_under = (col_under, col_over)
        # over-crossing curve bottom to middle
        s1 = (0, 0)
        c1 = (0.0, 0.25)
        e1 = (0.5, 0.5)
        cvob = (s1, c1, e1)

        # over-crossing curve middle to top
        c2 = (1.0, 0.75)
        e2 = (1.0, 1.0)
        cvot = (c2, e2)

        # under-crossing curve bottom to middle
        t1 = (1, 0)
        d1 = (1.0, 0.25)
        d2 = (0.5 + 4 * gap, 0.5 - 2 * gap)
        f1 = (0.5 + 2 * gap, 0.5 - gap)
        cvub = (t1, d1, d2, f1)

        # under-crossing curve middle to top
        u1 = (0.5 - 2 * gap, 0.5 + gap)
        g1 = (0.5 - 4 * gap, 0.5 + 2 * gap)
        g2 = (0.0, 0.75)
        v1 = (0.0, 1.0)
        cvut = (u1, g1, g2, v1)

        self._curves = {(cvob, cvot): col_over, (cvub,): col_under,  (cvut,): col_under}
        if not positive:
            self.reflect()
            self.shift(1, 0)
        super().__init__(pos=pos, rotation=rotation)


class _LinePlot(_GeneratorPlot):
    r"""
    Class for plotting a sraight line of the braid and cap-cup
    generators with given colour ``col``.

    INPUT:

    - ``col`` -- a valid color name for :meth:`~sage.plot.bezier_path`
      and :meth:`~sage.plot.line` used for the straight strand
    - ``pos`` -- see :class:`_GeneratorPlot`
    - ``rotation`` -- integer, giving the rotation as a factor of ``pi/2``
    """
    def __init__(self, col, pos=(0, 0), rotation=0):
        r"""
        Set the coordinate data for plotting ``self``.

        EXAMPLES::

            sage: from sage.monoids.tangles import _LinePlot
            sage: from sage.plot.colors import rainbow
            sage: col = rainbow(1)
            sage: cross = _LinePlot(col[0], pos=(1, 1), rotation=1)
            sage: cross._curves
            {((1, -1), (2, -1)): '#ff0000'}
        """
        # over-crossing curve bottom to middle
        p1 = (0, 0)
        p2 = (0, 1)

        self._curves = {(p1, p2): col}
        super().__init__(pos=pos, rotation=rotation)
