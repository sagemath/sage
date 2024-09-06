r"""
Knutson-Tao Puzzles

This module implements a generic algorithm to solve Knutson-Tao puzzles. An
instance of this class will be callable: the arguments are the labels of
north-east and north-west sides of the puzzle boundary; the output is the list
of the fillings of the puzzle with the specified pieces.

Acknowledgements
----------------

This code was written during Sage Days 45 at ICERM with Franco Saliola, Anne Schilling, and Avinash Dalal in discussions with Allen Knutson.
The code was tested afterwards by Liz Beazley and Ed Richmond.

.. TODO::

    Functionality to add:

    - plotter will not plot edge labels higher than 2; e.g. in BK puzzles, the labels are
      1,..., n and so in 3-step examples, none of the edge labels with 3 appear

    - we should also have a 3-step puzzle pieces constructor, taken from p22 of
      :arxiv:`math/0610538`

    - implement the bijection from puzzles to tableaux; see for example
      R. Vakil, A geometric Littlewood-Richardson rule, :arxiv:`math/0302294`
      or K. Purbhoo, Puzzles, Tableaux and Mosaics, :arxiv:`0705.1184`.
"""
# ****************************************************************************
#       Copyright (C) 2013 Franco Saliola <saliola@gmail.com>,
#                     2013 Allen Knutson,
#                     2013 Avinash Dalal,
#                     2013 Anne Schilling,
#                     2013 Elizabeth Beazley,
#                     2013 Ed Richmond
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from __future__ import annotations

from sage.misc.lazy_import import lazy_import
lazy_import("sage.plot.graphics", "Graphics")
lazy_import("sage.plot.polygon", "polygon")
lazy_import("sage.plot.line", "line")
lazy_import("sage.plot.text", "text")
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.finite_rings.integer_mod_ring import Integers
lazy_import("sage.plot.plot", "graphics_array")
from sage.misc.cachefunc import cached_method
from sage.structure.unique_representation import UniqueRepresentation


class PuzzlePiece:
    r"""
    Abstract class for puzzle pieces.

    This abstract class contains information on how to test equality of
    puzzle pieces, and sets color and plotting options.
    """

    def __eq__(self, other) -> bool:
        r"""
        TESTS::

            sage: from sage.combinat.knutson_tao_puzzles import DeltaPiece
            sage: delta = DeltaPiece('a','b','c')
            sage: delta1 = DeltaPiece('a','b','c')
            sage: delta == delta1
            True
            sage: delta1 = DeltaPiece('A','b','c')
            sage: delta == delta1
            False
        """
        if isinstance(other, PuzzlePiece):
            return self.border() == other.border()
        else:
            return False

    def __hash__(self):
        r"""
        TESTS::

            sage: from sage.combinat.knutson_tao_puzzles import DeltaPiece
            sage: delta = DeltaPiece('a','b','c')
            sage: hash(delta) == hash(delta)
            True
        """
        return hash((type(self), self.border()))

    def border(self) -> tuple:
        r"""
        Return the border of ``self``.

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import DeltaPiece
            sage: delta = DeltaPiece('a','b','c')
            sage: sorted(delta.border())
            ['a', 'b', 'c']
        """
        return tuple(self.edge_label(edge) for edge in self.edges())

    def color(self) -> str:
        r"""
        Return the color of ``self``.

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import DeltaPiece
            sage: delta = DeltaPiece('a','b','c')
            sage: delta.color()
            'white'
            sage: delta = DeltaPiece('0','0','0')
            sage: delta.color()
            'red'
            sage: delta = DeltaPiece('1','1','1')
            sage: delta.color()
            'blue'
            sage: delta = DeltaPiece('2','2','2')
            sage: delta.color()
            'green'
            sage: delta = DeltaPiece('2','K','2')
            sage: delta.color()
            'orange'
            sage: delta = DeltaPiece('2','T1/2','2')
            sage: delta.color()
            'yellow'
        """
        colors = {('0', '0', '0'): 'red',
                  ('1', '1', '1'): 'blue',
                  ('2', '2', '2'): 'green'}
        border = self.border()
        if border in colors:
            color = colors[border]
        elif 'K' in border:
            color = 'orange'
        elif '10' in border:
            color = 'white'
        elif any(label.startswith('T') for label in border):
            color = 'yellow'
        else:
            color = 'white'
        return color

    def _plot_label(self, label, coords, fontcolor=(0.3, 0.3, 0.3),
                    fontsize=15, rotation=0):
        r"""
        TESTS::

            sage: from sage.combinat.knutson_tao_puzzles import DeltaPiece
            sage: delta = DeltaPiece('2','K','2')
            sage: delta._plot_label('1',(1,1))    # not tested
        """
        if label in ('0', '1', '2'):
            return text(label, coords, color=fontcolor, fontsize=fontsize, rotation=rotation)
        else:
            return Graphics()

    def _plot_piece(self, coords, border_color=(0.5, 0.5, 0.5),
                    border_thickness=1, style='fill'):
        r"""
        TESTS::

            sage: from sage.combinat.knutson_tao_puzzles import DeltaPiece
            sage: delta = DeltaPiece('2','K','2')
            sage: delta._plot_piece([(1,1),(1,2),(2,2)]) # not tested
        """
        if style == 'fill':
            P = polygon(coords, color=self.color())
            P += polygon(coords, fill=False, color=border_color, thickness=border_thickness)
            return P
        elif style == 'edges':
            if isinstance(self, DeltaPiece):
                edges = ('north_west', 'south', 'north_east')
            elif isinstance(self, NablaPiece):
                edges = ('south_west', 'north', 'south_east')
            else:
                edges = self.edges()
            P = Graphics()
            for (i, edge) in enumerate(edges):
                P += line([coords[i], coords[(i + 1) % 3]],
                          color=self.edge_color(edge),
                          thickness=border_thickness)
            return P
        else:
            return NotImplemented

    def edge_color(self, edge) -> str:
        r"""
        Color of the specified edge of ``self`` (to be used when plotting the
        piece).

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import DeltaPiece
            sage: delta = DeltaPiece('1','0','10')
            sage: delta.edge_color('south')
            'blue'
            sage: delta.edge_color('north_west')
            'red'
            sage: delta.edge_color('north_east')
            'white'
        """
        edge_label = self.edge_label(edge)
        colors = {'1': 'blue', '0': 'red'}
        if edge_label in colors:
            color = colors[edge_label]
        elif 'K' in edge_label:
            color = 'orange'
        elif edge_label.startswith('T'):
            color = 'yellow'
        else:
            color = 'white'
        return color

    def edge_label(self, edge) -> str:
        r"""
        Return the edge label of ``edge``.

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import DeltaPiece
            sage: delta = DeltaPiece('2','K','2')
            sage: delta.edge_label('south')
            '2'
            sage: delta.edge_label('north_east')
            '2'
            sage: delta.edge_label('north_west')
            'K'
        """
        return self._edge_labels[edge]

    __getitem__ = edge_label


class NablaPiece(PuzzlePiece):
    r"""
    Nabla Piece takes as input three labels, inputted as strings. They label
    the North, Southeast and Southwest edges, respectively.

    EXAMPLES::

        sage: from sage.combinat.knutson_tao_puzzles import NablaPiece
        sage: NablaPiece('a','b','c')
        c\a/b
    """

    def __init__(self, north, south_east, south_west):
        r"""
        INPUT:

        - ``north``, ``south_east``, ``south_west`` -- strings, which label the edges

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import NablaPiece
            sage: NablaPiece('1','2','3')
            3\1/2
        """
        self._edge_labels = dict(north=north, south_east=south_east, south_west=south_west)

    def __eq__(self, other) -> bool:
        r"""
        TESTS::

            sage: from sage.combinat.knutson_tao_puzzles import NablaPiece
            sage: n = NablaPiece('a','b','c')
            sage: n1 = NablaPiece('a','b','c')
            sage: n == n1
            True
            sage: n1 = NablaPiece('A','b','c')
            sage: n == n1
            False
        """
        if isinstance(other, NablaPiece):
            return (self.border() == other.border() and
                    self._edge_labels == other._edge_labels)
        else:
            return False

    def __hash__(self):
        r"""
        TESTS::

            sage: from sage.combinat.knutson_tao_puzzles import NablaPiece
            sage: n = NablaPiece('a','b','c')
            sage: hash(n) == hash(n)
            True
        """
        return hash((NablaPiece, self.border()))

    def __repr__(self) -> str:
        r"""
        Print the labels of the Nabla piece.

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import NablaPiece
            sage: NablaPiece('1','2','3')
            3\1/2
        """
        return r"%s\%s/%s" % (self['south_west'],
                              self['north'],
                              self['south_east'])

    def clockwise_rotation(self) -> NablaPiece:
        r"""
        Rotate the Nabla piece by 120 degree clockwise.

        OUTPUT: Nabla piece

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import NablaPiece
            sage: nabla = NablaPiece('1','2','3')
            sage: nabla.clockwise_rotation()
            2\3/1
        """
        return NablaPiece(north=self['south_west'],
                          south_east=self['north'],
                          south_west=self['south_east'])

    def half_turn_rotation(self) -> DeltaPiece:
        r"""
        Rotate the Nabla piece by 180 degree.

        OUTPUT: Delta piece

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import NablaPiece
            sage: nabla = NablaPiece('1','2','3')
            sage: nabla.half_turn_rotation()
            2/1\3
        """
        return DeltaPiece(south=self['north'],
                          north_west=self['south_east'],
                          north_east=self['south_west'])

    def edges(self) -> tuple:
        r"""
        Return the tuple of edge names.

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import NablaPiece
            sage: nabla = NablaPiece('1','2','3')
            sage: nabla.edges()
            ('north', 'south_east', 'south_west')
        """
        return ('north', 'south_east', 'south_west')


class DeltaPiece(PuzzlePiece):
    r"""
    Delta Piece takes as input three labels, inputted as strings. They label
    the South, Northwest and Northeast edges, respectively.

    EXAMPLES::

        sage: from sage.combinat.knutson_tao_puzzles import DeltaPiece
        sage: DeltaPiece('a','b','c')
        b/a\c
    """

    def __init__(self, south, north_west, north_east):
        r"""
        INPUT:

        - ``south``, ``north_west``, ``north_east`` -- strings, which label the edges

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import DeltaPiece
            sage: DeltaPiece('1','2','3')
            2/1\3
        """
        self._edge_labels = dict(south=south, north_west=north_west, north_east=north_east)

    def __eq__(self, other) -> bool:
        r"""
        TESTS::

            sage: from sage.combinat.knutson_tao_puzzles import DeltaPiece
            sage: delta = DeltaPiece('a','b','c')
            sage: delta1 = DeltaPiece('a','b','c')
            sage: delta == delta1
            True
            sage: delta1 = DeltaPiece('A','b','c')
            sage: delta == delta1
            False
        """
        if isinstance(other, DeltaPiece):
            return (self.border() == other.border() and
                    self._edge_labels == other._edge_labels)
        else:
            return False

    def __hash__(self):
        r"""
        TESTS::

            sage: from sage.combinat.knutson_tao_puzzles import DeltaPiece
            sage: delta = DeltaPiece('a','b','c')
            sage: hash(delta) == hash(delta)
            True
        """
        return hash((DeltaPiece, self.border()))

    def __repr__(self) -> str:
        r"""
        Print the labels of the Delta piece.

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import DeltaPiece
            sage: DeltaPiece('1','2','3')
            2/1\3
        """
        return r"%s/%s\%s" % (self['north_west'],
                              self['south'],
                              self['north_east'])

    def clockwise_rotation(self) -> DeltaPiece:
        r"""
        Rotate the Delta piece by 120 degree clockwise.

        OUTPUT: Delta piece

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import DeltaPiece
            sage: delta = DeltaPiece('1','2','3')
            sage: delta.clockwise_rotation()
            1/3\2
        """
        return DeltaPiece(south=self['north_east'],
                          north_west=self['south'],
                          north_east=self['north_west'])

    def half_turn_rotation(self) -> NablaPiece:
        r"""
        Rotate the Delta piece by 180 degree.

        OUTPUT: Nabla piece

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import DeltaPiece
            sage: delta = DeltaPiece('1','2','3')
            sage: delta.half_turn_rotation()
            3\1/2
        """
        return NablaPiece(north=self['south'],
                          south_east=self['north_west'],
                          south_west=self['north_east'])

    def edges(self) -> tuple:
        r"""
        Return the tuple of edge names.

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import DeltaPiece
            sage: delta = DeltaPiece('1','2','3')
            sage: delta.edges()
            ('south', 'north_west', 'north_east')
        """
        return ('south', 'north_west', 'north_east')


class RhombusPiece(PuzzlePiece):
    r"""
    Class of rhombi pieces.

    To construct a rhombus piece we input a delta and a nabla piece.
    The delta and nabla pieces are joined along the south and north edge,
    respectively.

    EXAMPLES::

        sage: from sage.combinat.knutson_tao_puzzles import DeltaPiece, NablaPiece, RhombusPiece
        sage: delta = DeltaPiece('1','2','3')
        sage: nabla = NablaPiece('4','5','6')
        sage: RhombusPiece(delta,nabla)
        2/\3  6\/5
    """

    def __init__(self, north_piece, south_piece):
        r"""
        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import DeltaPiece, NablaPiece, RhombusPiece
            sage: delta = DeltaPiece('1','2','3')
            sage: nabla = NablaPiece('4','5','6')
            sage: RhombusPiece(delta,nabla)
            2/\3  6\/5
        """
        self._north_piece = north_piece
        self._south_piece = south_piece
        self._edge_labels = dict(north_west=north_piece['north_west'],
                                 north_east=north_piece['north_east'],
                                 south_east=south_piece['south_east'],
                                 south_west=south_piece['south_west'])

    def __eq__(self, other) -> bool:
        r"""
        TESTS::

            sage: from sage.combinat.knutson_tao_puzzles import DeltaPiece, NablaPiece, RhombusPiece
            sage: delta = DeltaPiece('1','2','3')
            sage: nabla = NablaPiece('4','5','6')
            sage: r = RhombusPiece(delta,nabla)
            sage: r == r
            True
            sage: delta1 = DeltaPiece('A','b','c')
            sage: r == RhombusPiece(delta1,nabla)
            False
        """
        if isinstance(other, RhombusPiece):
            return (self.border() == other.border() and
                    self._north_piece == other._north_piece and
                    self._south_piece == other._south_piece and
                    self._edge_labels == other._edge_labels)
        else:
            return False

    def __hash__(self):
        r"""
        TESTS::

            sage: from sage.combinat.knutson_tao_puzzles import DeltaPiece, NablaPiece, RhombusPiece
            sage: delta = DeltaPiece('1','2','3')
            sage: nabla = NablaPiece('4','5','6')
            sage: r = RhombusPiece(delta,nabla)
            sage: hash(r) == hash(r)
            True
        """
        return hash((RhombusPiece, self.border()))

    def __iter__(self):
        r"""
        Return the list of the north and south piece.

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import DeltaPiece, NablaPiece, RhombusPiece
            sage: delta = DeltaPiece('1','2','3')
            sage: nabla = NablaPiece('4','5','6')
            sage: r = RhombusPiece(delta,nabla)
            sage: list(r)
            [2/1\3, 6\4/5]
        """
        yield self._north_piece
        yield self._south_piece

    def north_piece(self) -> DeltaPiece:
        r"""
        Return the north piece.

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import DeltaPiece, NablaPiece, RhombusPiece
            sage: delta = DeltaPiece('1','2','3')
            sage: nabla = NablaPiece('4','5','6')
            sage: r = RhombusPiece(delta,nabla)
            sage: r.north_piece()
            2/1\3
        """
        return self._north_piece

    def south_piece(self) -> NablaPiece:
        r"""
        Return the south piece.

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import DeltaPiece, NablaPiece, RhombusPiece
            sage: delta = DeltaPiece('1','2','3')
            sage: nabla = NablaPiece('4','5','6')
            sage: r = RhombusPiece(delta,nabla)
            sage: r.south_piece()
            6\4/5
        """
        return self._south_piece

    def __repr__(self) -> str:
        r"""
        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import DeltaPiece, NablaPiece, RhombusPiece
            sage: delta = DeltaPiece('1','2','3')
            sage: nabla = NablaPiece('4','5','6')
            sage: RhombusPiece(delta,nabla)
            2/\3  6\/5
        """
        return r"%s/\%s  %s\/%s" % (self['north_west'], self['north_east'],
                                    self['south_west'], self['south_east'])

    def edges(self) -> tuple:
        r"""
        Return the tuple of edge names.

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import DeltaPiece, NablaPiece, RhombusPiece
            sage: delta = DeltaPiece('1','2','3')
            sage: nabla = NablaPiece('4','5','6')
            sage: RhombusPiece(delta,nabla).edges()
            ('north_west', 'north_east', 'south_east', 'south_west')
        """
        return ('north_west', 'north_east', 'south_east', 'south_west')


class PuzzlePieces:
    r"""
    Construct a valid set of puzzle pieces.

    This class constructs the set of valid puzzle pieces. It can take a list of
    forbidden border labels as input. These labels are forbidden from appearing
    on the south edge of a puzzle filling. The user can add valid nabla or
    delta pieces and specify which rotations of these pieces are legal. For
    example, ``rotations=0`` does not add any additional pieces (only the piece
    itself), ``rotations=60`` adds six pieces (the pieces and its rotations by
    60, 120, 180, 240, 300), etc..

    EXAMPLES::

        sage: from sage.combinat.knutson_tao_puzzles import PuzzlePieces, NablaPiece
        sage: forbidden_border_labels = ['10']
        sage: pieces = PuzzlePieces(forbidden_border_labels)
        sage: pieces.add_piece(NablaPiece('0','0','0'), rotations=60)
        sage: pieces.add_piece(NablaPiece('1','1','1'), rotations=60)
        sage: pieces.add_piece(NablaPiece('1','0','10'), rotations=60)
        sage: pieces
        Nablas : [0\0/0, 0\10/1, 10\1/0, 1\0/10, 1\1/1]
        Deltas : [0/0\0, 0/1\10, 1/10\0, 1/1\1, 10/0\1]

    The user can obtain the list of valid rhombi pieces as follows::

        sage: sorted([p for p in pieces.rhombus_pieces()], key=str)
        [0/\0  0\/0, 0/\0  1\/10, 0/\10  10\/0, 0/\10  1\/1, 1/\0  0\/1,
        1/\1  10\/0, 1/\1  1\/1, 10/\1  0\/0, 10/\1  1\/10]
    """

    def __init__(self, forbidden_border_labels=None):
        r"""
        INPUT:

        - ``forbidden_border_labels`` -- list of forbidden border labels given as strings

        TESTS::

            sage: from sage.combinat.knutson_tao_puzzles import PuzzlePieces
            sage: forbidden_border_labels = ['10']
            sage: pieces = PuzzlePieces(forbidden_border_labels)
            sage: pieces
            Nablas : []
            Deltas : []

            sage: PuzzlePieces('10')
            Traceback (most recent call last):
            ...
            TypeError: Input must be a list
        """
        self._nabla_pieces = set()
        self._delta_pieces = set()
        if forbidden_border_labels is None:
            forbidden_border_labels = []
        if not isinstance(forbidden_border_labels, list):
            raise TypeError("Input must be a list")
        self._forbidden_border_labels = forbidden_border_labels

    def __eq__(self, other) -> bool:
        r"""
        TESTS::

            sage: from sage.combinat.knutson_tao_puzzles import H_grassmannian_pieces
            sage: x = H_grassmannian_pieces()
            sage: y = H_grassmannian_pieces()
            sage: x == y
            True
        """
        if isinstance(other, type(self)):
            return self.__dict__ == other.__dict__
        else:
            return False

    def __hash__(self):
        r"""
        TESTS::

            sage: from sage.combinat.knutson_tao_puzzles import H_grassmannian_pieces
            sage: x = H_grassmannian_pieces()
            sage: hash(x) == hash(x)
            True
        """
        return hash((type(self), repr(self)))

    def add_piece(self, piece, rotations=120):
        r"""
        Add ``piece`` to the list of pieces.

        INPUT:

        - ``piece`` -- a nabla piece or a delta piece
        - ``rotations`` -- (default: 120) 0, 60, 120, 180

        The user can add valid nabla or delta pieces and specify
        which rotations of these pieces are legal. For example, ``rotations=0``
        does not add any additional pieces (only the piece itself), ``rotations=60`` adds
        six pieces (namely three delta and three nabla pieces), while
        ``rotations=120`` adds only delta or nabla (depending on which piece ``self`` is).
        ``rotations=180`` adds the piece and its 180 degree rotation, i.e. one delta and one
        nabla piece.

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import PuzzlePieces, DeltaPiece
            sage: delta = DeltaPiece('a','b','c')
            sage: pieces = PuzzlePieces()
            sage: pieces
            Nablas : []
            Deltas : []
            sage: pieces.add_piece(delta)
            sage: pieces
            Nablas : []
            Deltas : [a/c\b, b/a\c, c/b\a]

            sage: pieces = PuzzlePieces()
            sage: pieces.add_piece(delta,rotations=0)
            sage: pieces
            Nablas : []
            Deltas : [b/a\c]

            sage: pieces = PuzzlePieces()
            sage: pieces.add_piece(delta,rotations=60)
            sage: pieces
            Nablas : [a\b/c, b\c/a, c\a/b]
            Deltas : [a/c\b, b/a\c, c/b\a]
        """
        if isinstance(piece, NablaPiece):
            pieces_list = self._nabla_pieces
        else:
            pieces_list = self._delta_pieces
        pieces_list.add(piece)
        if rotations == 120:
            pieces_list.add(piece.clockwise_rotation())
            pieces_list.add(piece.clockwise_rotation().clockwise_rotation())
        elif rotations == 180:
            self.add_piece(piece.half_turn_rotation(), rotations=0)
        elif rotations == 60:
            self.add_piece(piece, rotations=120)
            self.add_piece(piece.half_turn_rotation(), rotations=120)

    def add_forbidden_label(self, label):
        r"""
        Add forbidden border labels.

        INPUT:

        - ``label`` -- string specifying a new forbidden label

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import PuzzlePieces
            sage: pieces = PuzzlePieces()
            sage: pieces.add_forbidden_label('1')
            sage: pieces._forbidden_border_labels
            ['1']
            sage: pieces.add_forbidden_label('2')
            sage: pieces._forbidden_border_labels
            ['1', '2']
        """
        self._forbidden_border_labels.append(label)

    def add_T_piece(self, label1, label2):
        r"""
        Add a nabla and delta piece with ``label1`` and ``label2``.

        This method adds a nabla piece with edges ``label2``\ T``label1``|``label2`` / ``label1``.
        and a delta piece with edges ``label1``/ T``label1``|``label2`` \ ``label2``.
        It also adds T``label1``|``label2`` to the forbidden list.

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import PuzzlePieces
            sage: pieces = PuzzlePieces()
            sage: pieces.add_T_piece('1','3')
            sage: pieces
            Nablas : [3\T1|3/1]
            Deltas : [1/T1|3\3]
            sage: pieces._forbidden_border_labels
            ['T1|3']
        """
        self.add_forbidden_label('T%s|%s' % (label1, label2))
        self.add_piece(NablaPiece('T%s|%s' % (label1, label2), label1, label2), rotations=180)

    def __repr__(self) -> str:
        r"""
        TESTS::

            sage: from sage.combinat.knutson_tao_puzzles import PuzzlePieces, DeltaPiece
            sage: pieces = PuzzlePieces()
            sage: delta = DeltaPiece('a','b','c')
            sage: pieces.add_piece(delta,rotations=60)
            sage: pieces
            Nablas : [a\b/c, b\c/a, c\a/b]
            Deltas : [a/c\b, b/a\c, c/b\a]
        """
        s = "Nablas : %s\n" % sorted(self._nabla_pieces, key=str)
        s += "Deltas : %s" % sorted(self._delta_pieces, key=str)
        return s

    def delta_pieces(self):
        r"""
        Return the delta pieces as a set.

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import PuzzlePieces, DeltaPiece
            sage: pieces = PuzzlePieces()
            sage: delta = DeltaPiece('a','b','c')
            sage: pieces.add_piece(delta,rotations=60)
            sage: sorted([p for p in pieces.delta_pieces()], key=str)
            [a/c\b, b/a\c, c/b\a]
        """
        return self._delta_pieces

    def nabla_pieces(self):
        r"""
        Return the nabla pieces as a set.

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import PuzzlePieces, DeltaPiece
            sage: pieces = PuzzlePieces()
            sage: delta = DeltaPiece('a','b','c')
            sage: pieces.add_piece(delta,rotations=60)
            sage: sorted([p for p in pieces.nabla_pieces()], key=str)
            [a\b/c, b\c/a, c\a/b]
        """
        return self._nabla_pieces

    def rhombus_pieces(self) -> set:
        r"""
        Return a set of all allowable rhombus pieces.

        Allowable rhombus pieces are those where the south edge of the delta
        piece equals the north edge of the nabla piece.

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import PuzzlePieces, DeltaPiece
            sage: pieces = PuzzlePieces()
            sage: delta = DeltaPiece('a','b','c')
            sage: pieces.add_piece(delta,rotations=60)
            sage: sorted([p for p in pieces.rhombus_pieces()], key=str)
            [a/\b  b\/a, b/\c  c\/b, c/\a  a\/c]
        """
        rhombi = set()
        for nabla in self._nabla_pieces:
            for delta in self._delta_pieces:
                if delta['south'] == nabla['north']:
                    rhombi.add(RhombusPiece(delta, nabla))
        return rhombi

    def boundary_deltas(self) -> tuple:
        r"""
        Return deltas with south edges not in the forbidden list.

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import PuzzlePieces, DeltaPiece
            sage: pieces = PuzzlePieces(['a'])
            sage: delta = DeltaPiece('a','b','c')
            sage: pieces.add_piece(delta,rotations=60)
            sage: sorted([p for p in pieces.boundary_deltas()], key=str)
            [a/c\b, c/b\a]
        """
        return tuple(delta for delta in self.delta_pieces()
                    if delta['south'] not in self._forbidden_border_labels)


def H_grassmannian_pieces():
    r"""
    Define the puzzle pieces used in computing the cohomology of the Grassmannian.

    REFERENCES:

    .. [KTW] Allen Knutson, Terence Tao, Christopher Woodward,
       The honeycomb model of GL(n) tensor products II: Puzzles determine facets of the Littlewood-Richardson cone,
       :arxiv:`math/0107011`

    EXAMPLES::

        sage: from sage.combinat.knutson_tao_puzzles import H_grassmannian_pieces
        sage: H_grassmannian_pieces()
        Nablas : [0\0/0, 0\10/1, 10\1/0, 1\0/10, 1\1/1]
        Deltas : [0/0\0, 0/1\10, 1/10\0, 1/1\1, 10/0\1]
    """
    forbidden_border_labels = ['10']
    pieces = PuzzlePieces(forbidden_border_labels)
    pieces.add_piece(NablaPiece('0', '0', '0'), rotations=60)
    pieces.add_piece(NablaPiece('1', '1', '1'), rotations=60)
    pieces.add_piece(NablaPiece('1', '0', '10'), rotations=60)
    return pieces


def HT_grassmannian_pieces():
    r"""
    Define the puzzle pieces used in computing the torus-equivariant cohomology of the Grassmannian.

    REFERENCES:

    .. [KT2003] Allen Knutson, Terence Tao, Puzzles and (equivariant) cohomology of Grassmannians,
       Duke Math. J. 119 (2003) 221

    EXAMPLES::

        sage: from sage.combinat.knutson_tao_puzzles import HT_grassmannian_pieces
        sage: HT_grassmannian_pieces()
        Nablas : [0\0/0, 0\10/1, 10\1/0, 1\0/10, 1\1/1, 1\T0|1/0]
        Deltas : [0/0\0, 0/1\10, 0/T0|1\1, 1/10\0, 1/1\1, 10/0\1]
    """
    pieces = H_grassmannian_pieces()
    pieces.add_T_piece('0', '1')
    return pieces


def K_grassmannian_pieces():
    r"""
    Define the puzzle pieces used in computing the K-theory of the Grassmannian.

    REFERENCES:

    .. [Buch00] \A. Buch, A Littlewood-Richardson rule for the K-theory of Grassmannians, :arxiv:`math.AG/0004137`

    EXAMPLES::

        sage: from sage.combinat.knutson_tao_puzzles import K_grassmannian_pieces
        sage: K_grassmannian_pieces()
        Nablas : [0\0/0, 0\10/1, 0\K/1, 10\1/0, 1\0/10, 1\0/K, 1\1/1, K\1/0]
        Deltas : [0/0\0, 0/1\10, 1/10\0, 1/1\1, 10/0\1, K/K\K]
    """
    pieces = H_grassmannian_pieces()
    pieces.add_forbidden_label('K')
    pieces.add_piece(NablaPiece('0', 'K', '1'), rotations=120)
    pieces.add_piece(DeltaPiece('K', 'K', 'K'), rotations=0)
    return pieces


def H_two_step_pieces():
    r"""
    Define the puzzle pieces used in two step flags.

    This rule is currently only conjecturally true. See [BuchKreschTamvakis03]_.

    REFERENCES:

    .. [BuchKreschTamvakis03] \A. Buch, A. Kresch, H. Tamvakis, Gromov-Witten invariants on Grassmannian, :arxiv:`math/0306388`

    EXAMPLES::

        sage: from sage.combinat.knutson_tao_puzzles import H_two_step_pieces
        sage: H_two_step_pieces()
        Nablas : [(21)0\21/0, 0\(21)0/21, 0\0/0, 0\10/1, 0\20/2, 10\1/0, 10\2(10)/2, 1\0/10, 1\1/1, 1\21/2,
        2(10)\2/10, 20\2/0, 21\0/(21)0, 21\2/1, 2\0/20, 2\1/21, 2\10/2(10), 2\2/2]
        Deltas : [(21)0/0\21, 0/0\0, 0/1\10, 0/21\(21)0, 0/2\20, 1/10\0, 1/1\1, 1/2\21, 10/0\1, 10/2\2(10),
        2(10)/10\2, 2/2(10)\10, 2/20\0, 2/21\1, 2/2\2, 20/0\2, 21/(21)0\0, 21/1\2]
    """
    forbidden_border_labels = ['10', '20', '21', '(21)0', '2(10)']
    pieces = PuzzlePieces(forbidden_border_labels)
    for i in ('0', '1', '2'):
        pieces.add_piece(DeltaPiece(i, i, i), rotations=60)
    for i, j in (('1', '0'), ('2', '0'), ('2', '1')):
        pieces.add_piece(DeltaPiece(i + j, i, j), rotations=60)
    pieces.add_piece(DeltaPiece('(21)0', '21', '0'), rotations=60)
    pieces.add_piece(DeltaPiece('2(10)', '2', '10'), rotations=60)
    return pieces


def HT_two_step_pieces():
    r"""
    Define the puzzle pieces used in computing the equivariant two step puzzle pieces.

    For the puzzle pieces, see Figure 26 on page 22 of [CoskunVakil06]_.

    REFERENCES:

    .. [CoskunVakil06] \I. Coskun, R. Vakil, Geometric positivity in the cohomology of homogeneous spaces
       and generalized Schubert calculus, :arxiv:`math/0610538`

    EXAMPLES::

       sage: from sage.combinat.knutson_tao_puzzles import HT_two_step_pieces
       sage: HT_two_step_pieces()
       Nablas : [(21)0\21/0, 0\(21)0/21, 0\0/0, 0\10/1, 0\20/2, 10\1/0, 10\2(10)/2,
       1\0/10, 1\1/1, 1\21/2, 1\T0|1/0, 2(10)\2/10, 20\2/0, 21\0/(21)0, 21\2/1, 21\T0|21/0,
       21\T10|21/10, 2\0/20, 2\1/21, 2\10/2(10), 2\2/2, 2\T0|2/0, 2\T10|2/10, 2\T1|2/1]
       Deltas : [(21)0/0\21, 0/0\0, 0/1\10, 0/21\(21)0, 0/2\20, 0/T0|1\1, 0/T0|21\21, 0/T0|2\2,
       1/10\0, 1/1\1, 1/2\21, 1/T1|2\2, 10/0\1, 10/2\2(10), 10/T10|21\21, 10/T10|2\2, 2(10)/10\2,
       2/2(10)\10, 2/20\0, 2/21\1, 2/2\2, 20/0\2, 21/(21)0\0, 21/1\2]
    """
    pieces = H_two_step_pieces()
    for label1, label2 in (('0', '1'), ('0', '2'), ('1', '2'),
                           ('10', '2'), ('0', '21'), ('10', '21')):
        pieces.add_T_piece(label1, label2)
    return pieces


def BK_pieces(max_letter):
    r"""
    The puzzle pieces used in computing the Belkale-Kumar coefficients for any
    partial flag variety in type `A`.

    There are two types of puzzle pieces:

    - a triangle, with each edge labeled with the same letter;
    - a rhombus, with edges labeled `i`, `j`, `i`, `j` in clockwise order with
      `i > j`.

    Each of these is rotated by 60 degrees, but not reflected.

    We model the rhombus pieces as two triangles: a delta piece north-west
    label `i`, north-east label `j` and south label `i(j)`; and a nabla piece
    with south-east label `i`, south-west label `j` and north label `i(j)`.

    INPUT:

    - ``max_letter`` -- positive integer specifying the number of steps in the
      partial flag variety, equivalently, the number of elements in the
      alphabet for the edge labels. The smallest label is `1`.

    REFERENCES:

    .. [KnutsonPurbhoo10] \A. Knutson, K. Purbhoo, Product and puzzle formulae
       for `GL_n` Belkale-Kumar coefficients, :arxiv:`1008.4979`

    EXAMPLES::

        sage: from sage.combinat.knutson_tao_puzzles import BK_pieces
        sage: BK_pieces(3)
        Nablas : [1\1/1, 1\2(1)/2, 1\3(1)/3, 2(1)\2/1, 2\1/2(1), 2\2/2, 2\3(2)/3, 3(1)\3/1, 3(2)\3/2, 3\1/3(1), 3\2/3(2), 3\3/3]
        Deltas : [1/1\1, 1/2\2(1), 1/3\3(1), 2(1)/1\2, 2/2(1)\1, 2/2\2, 2/3\3(2), 3(1)/1\3, 3(2)/2\3, 3/3(1)\1, 3/3(2)\2, 3/3\3]
    """
    forbidden_border_labels = ['%s(%s)' % (i, j)
                               for i in range(1, max_letter + 1)
                               for j in range(1, i)]
    pieces = PuzzlePieces(forbidden_border_labels)
    for i in range(1, max_letter + 1):
        piece = DeltaPiece('%s' % i, '%s' % i, '%s' % i)
        pieces.add_piece(piece, rotations=60)
        for j in range(1, i):
            piece = DeltaPiece(north_west='%s' % i, north_east='%s' % j,
                               south='%s(%s)' % (i, j))
            pieces.add_piece(piece, rotations=60)
    return pieces


class PuzzleFilling:
    r"""
    Create partial puzzles and provides methods to build puzzles from them.
    """

    def __init__(self, north_west_labels, north_east_labels):
        r"""
        TESTS::

            sage: from sage.combinat.knutson_tao_puzzles import PuzzleFilling
            sage: P = PuzzleFilling('0101','0101')
            sage: P
            {}
        """
        self._nw_labels = tuple(north_west_labels)
        self._ne_labels = tuple(north_east_labels)
        self._squares = {}
        self._n = len(self._nw_labels)
        self._kink_coordinates = (1, self._n)

    def __getitem__(self, key):
        r"""
        TESTS::

            sage: from sage.combinat.knutson_tao_puzzles import KnutsonTaoPuzzleSolver
            sage: ps = KnutsonTaoPuzzleSolver("H")
            sage: puzzle = ps('0101','1001')[0]
            sage: puzzle
            {(1, 1): 0/1\10,
             (1, 2): 1/\1  10\/0,
             (1, 3): 0/\10  1\/1,
             (1, 4): 1/\1  10\/0,
             (2, 2): 0/0\0,
             (2, 3): 1/\0  0\/1,
             (2, 4): 0/\0  0\/0,
             (3, 3): 1/1\1,
             (3, 4): 0/\0  1\/10,
             (4, 4): 10/0\1}
            sage: puzzle[(1,2)]  # indirect doctest
            1/\1  10\/0
        """
        return self._squares[key]

    def kink_coordinates(self) -> tuple:
        r"""
        Provide the coordinates of the kinks.

        The kink coordinates are the coordinates up to which the puzzle has already
        been built. The kink starts in the north corner and then moves down the diagonals
        as the puzzles is built.

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import PuzzleFilling
            sage: P = PuzzleFilling('0101','0101')
            sage: P
            {}
            sage: P.kink_coordinates()
            (1, 4)
        """
        return self._kink_coordinates

    def is_in_south_edge(self) -> bool:
        r"""
        Check whether kink coordinates of partial puzzle is in south corner.

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import PuzzleFilling
            sage: P = PuzzleFilling('0101','0101')
            sage: P.is_in_south_edge()
            False
        """
        i, j = self.kink_coordinates()
        return i == j

    def north_west_label_of_kink(self):
        r"""
        Return north-west label of kink.

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import PuzzleFilling
            sage: P = PuzzleFilling('0101','0101')
            sage: P.north_west_label_of_kink()
            '1'
        """
        (i, j) = self.kink_coordinates()
        if i == 1:
            return self._nw_labels[j - 1]
        else:
            return self._squares[i - 1, j]['south_east']

    def north_east_label_of_kink(self):
        r"""
        Return north east label of kink.

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import PuzzleFilling
            sage: P = PuzzleFilling('0101','0101')
            sage: P.north_east_label_of_kink()
            '0'
        """
        (i, j) = self.kink_coordinates()
        if j == self._n:
            return self._ne_labels[i - 1]
        else:
            return self._squares[i, j + 1]['south_west']

    def is_completed(self):
        r"""
        Whether partial puzzle is complete (completely filled) or not.

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import PuzzleFilling
            sage: P = PuzzleFilling('0101','0101')
            sage: P.is_completed()
            False

            sage: from sage.combinat.knutson_tao_puzzles import KnutsonTaoPuzzleSolver
            sage: ps = KnutsonTaoPuzzleSolver("H")
            sage: puzzle = ps('0101','1001')[0]
            sage: puzzle.is_completed()
            True
        """
        i, _ = self.kink_coordinates()
        return i == self._n + 1

    def south_labels(self):
        r"""
        Return south labels for completed puzzle.

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import KnutsonTaoPuzzleSolver
            sage: ps = KnutsonTaoPuzzleSolver("H")
            sage: ps('0101','1001')[0].south_labels()
            ('1', '0', '1', '0')
        """
        # TODO: return ''.join(self[i, i]['south'] for i in range(1, self._n + 1))
        return tuple([self[i, i]['south'] for i in range(1, self._n + 1)])

    def add_piece(self, piece):
        r"""
        Add ``piece`` to partial puzzle.

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import DeltaPiece, PuzzleFilling
            sage: piece = DeltaPiece('0','1','0')
            sage: P = PuzzleFilling('0101','0101'); P
            {}
            sage: P.add_piece(piece); P
            {(1, 4): 1/0\0}
        """
        i, j = self.kink_coordinates()
        self._squares[i, j] = piece
        if isinstance(piece, DeltaPiece):
            i += 1
            j = self._n
        else:
            j -= 1
        self._kink_coordinates = (i, j)

    def add_pieces(self, pieces):
        r"""
        Add ``piece`` to partial puzzle.

        INPUT:

        - ``pieces`` -- tuple of pieces

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import DeltaPiece, PuzzleFilling
            sage: P = PuzzleFilling('0101','0101'); P
            {}
            sage: piece = DeltaPiece('0','1','0')
            sage: pieces = [piece,piece]
            sage: P.add_pieces(pieces)
            sage: P
            {(1, 4): 1/0\0, (2, 4): 1/0\0}
        """
        i, j = self.kink_coordinates()
        for piece in pieces:
            self._squares[i, j] = piece
            if isinstance(piece, DeltaPiece):
                i += 1
                j = self._n
            else:
                j -= 1
        self._kink_coordinates = (i, j)

    def copy(self):
        r"""
        Return copy of ``self``.

        EXAMPLES::


            sage: from sage.combinat.knutson_tao_puzzles import DeltaPiece, PuzzleFilling
            sage: piece = DeltaPiece('0','1','0')
            sage: P = PuzzleFilling('0101','0101'); P
            {}
            sage: PP = P.copy()
            sage: P.add_piece(piece); P
            {(1, 4): 1/0\0}
            sage: PP
            {}
        """
        PP = PuzzleFilling(self._nw_labels, self._ne_labels)
        PP._squares = self._squares.copy()
        PP._kink_coordinates = self._kink_coordinates
        PP._n = self._n
        return PP

    def contribution(self):
        r"""
        Return equivariant contributions from ``self`` in polynomial ring.

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import KnutsonTaoPuzzleSolver
            sage: ps = KnutsonTaoPuzzleSolver("HT")
            sage: puzzles = ps('0101','1001')
            sage: sorted([p.contribution() for p in puzzles], key=str)
            [1, y1 - y3]
        """
        R = PolynomialRing(Integers(), 'y', self._n + 1)
        y = R.gens()
        z = R.one()
        for i in range(1, self._n + 1):
            for j in range(i + 1, self._n + 1):
                if self[i, j].north_piece()['south'].startswith('T'):
                    z *= y[i] - y[j]
                if self[i, j].north_piece()['south'].startswith('K'):
                    z *= -1
        return z

    def __repr__(self):
        r"""
        TESTS::

            sage: from sage.combinat.knutson_tao_puzzles import H_grassmannian_pieces, PuzzleFilling
            sage: P = PuzzleFilling('0101','0101'); P
            {}
            sage: P.__repr__()
            '{}'
        """
        from pprint import pformat
        return pformat(self._squares)

    def __iter__(self):
        r"""
        Iterator.

        TESTS::

            sage: from sage.combinat.knutson_tao_puzzles import KnutsonTaoPuzzleSolver
            sage: ps = KnutsonTaoPuzzleSolver("H")
            sage: puzzle = ps('0101','1001')[0]
            sage: puzzle
            {(1, 1): 0/1\10,
             (1, 2): 1/\1  10\/0,
             (1, 3): 0/\10  1\/1,
             (1, 4): 1/\1  10\/0,
             (2, 2): 0/0\0,
             (2, 3): 1/\0  0\/1,
             (2, 4): 0/\0  0\/0,
             (3, 3): 1/1\1,
             (3, 4): 0/\0  1\/10,
             (4, 4): 10/0\1}
            sage: [p for p in puzzle]
            [1/\1  10\/0,
            0/\10  1\/1,
            0/\0  0\/0,
            1/\1  10\/0,
            1/\0  0\/1,
            0/\0  1\/10,
            0/1\10,
            0/0\0,
            1/1\1,
            10/0\1]
        """
        for d in range(self._n):
            for k in range(d + 1):
                yield self[k + 1, self._n - d + k]

    def plot(self, labels=True, style='fill'):
        r"""
        Plot completed puzzle.

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import KnutsonTaoPuzzleSolver
            sage: ps = KnutsonTaoPuzzleSolver("H")
            sage: puzzle = ps('0101','1001')[0]
            sage: puzzle.plot()  #not tested
            sage: puzzle.plot(style='fill')  #not tested
            sage: puzzle.plot(style='edges')  #not tested
        """
        P = Graphics()
        coords = [(k, -d) for d in range(self._n) for k in range(-d, d + 1, 2)]
        for ((k, d), piece) in zip(coords, self):
            if isinstance(piece, RhombusPiece):
                for (i, triangle) in enumerate(piece):
                    P += triangle._plot_piece([(k, d - 2 * i), (k - 1, d - 1), (k + 1, d - 1)], style=style)
                if labels:
                    P += piece._plot_label(piece['north_west'], (k - 0.5, d - 0.5), rotation=60)
                    P += piece._plot_label(piece['north_east'], (k + 0.5, d - 0.5), rotation=-60)
                    P += piece._plot_label(piece.north_piece()['south'], (k, d - 1))
            else:
                P += piece._plot_piece([(k, d), (k - 1, d - 1), (k + 1, d - 1)], style=style)
                if labels:
                    P += piece._plot_label(piece['north_west'], (k - 0.5, d - 0.5), rotation=60)
                    P += piece._plot_label(piece['north_east'], (k + 0.5, d - 0.5), rotation=-60)
                    P += piece._plot_label(piece['south'], (k, d - 1))
        P.set_aspect_ratio(1.73)
        P.axes(False)
        return P

    def _latex_(self):
        r"""
        Return latex version of ``self``.

        Note that you might need to add tikz to the preamble::

            sage: latex.extra_preamble(r'''\usepackage{tikz}''')
            sage: from sage.combinat.knutson_tao_puzzles import *

            sage: ps = KnutsonTaoPuzzleSolver(H_grassmannian_pieces())
            sage: solns = ps('0101', '0101')
            sage: view(solns[0], viewer='pdf')  # not tested

            sage: ps = KnutsonTaoPuzzleSolver(HT_two_step_pieces())
            sage: solns = ps(list('10212'), list('12012'))
            sage: view(solns[0], viewer='pdf')  # not tested

            sage: ps = KnutsonTaoPuzzleSolver(K_grassmannian_pieces())
            sage: solns = ps('0101', '0101')
            sage: view(solns[0], viewer='pdf')  # not tested
        """
        from collections import defaultdict
        label_colors = defaultdict(lambda: None)
        label_colors.update({'0': 'red', '1': 'blue', '2': 'green'})
        edge_colors = defaultdict(lambda: None)
        edge_colors.update({'0': 'red', '1': 'blue',
                            '2': 'green', 'K': 'orange'})

        s = r"""\begin{tikzpicture}[yscale=1.73]"""
        coords = [(k, -d) for d in range(self._n) for k in range(-d, d + 1, 2)]

        def tikztriangle_fill(color, k, d, i, *args):
            s = r"""\path[color=%s, fill=%s!10]""" % (color, color)
            s += r"""(%s, %s) -- (%s, %s)""" % (k, d - 2 * i, k - 1, d - 1)
            s += r"""-- (%s, %s)""" % (k + 1, d - 1)
            s += r"""-- (%s, %s)""" % (k, d - 2 * i)
            s += ";\n"
            return s

        def tikztriangle_edges(color, k, d, i, label1, label2, label3):
            s = ""
            if i == 1:
                return s
            tikzcmd = r"""\draw[color=%s, fill=none] (%s, %s) -- (%s, %s);""" + "\n"
            if edge_colors[label1]:
                s += tikzcmd % (edge_colors[label1], k - 1, d - 1, k + 1, d - 1)
            if edge_colors[label2]:
                s += tikzcmd % (edge_colors[label2], k, d - 2 * i, k - 1, d - 1)
            if edge_colors[label3]:
                s += tikzcmd % (edge_colors[label3], k + 1, d - 1, k, d - 2 * i)
            return s

        def tikzlabels(color, k, d, i, label1, label2, label3):
            s = r"""\path[] (%s, %s)""" % (k, d - 2 * i)
            s += r"""-- (%s, %s) """ % (k - 1, d - 1)
            if label_colors[label2]:
                s += r"""node[midway, color=%s] {$%s$} """ % (label_colors[label2], label2)
            s += r"""-- (%s, %s) """ % (k + 1, d - 1)
            if label_colors[label1]:
                s += r"""node[midway, color=%s] {$%s$} """ % (label_colors[label1], label1)
            s += r"""-- (%s, %s) """ % (k, d - 2 * i)
            if label_colors[label3]:
                s += r"""node[midway, color=%s] {$%s$} """ % (label_colors[label3], label3)
            s += ";\n"
            return s

        for ((k, d), piece) in zip(coords, self):
            for tikzcmd in (tikztriangle_fill, tikztriangle_edges, tikzlabels):
                if isinstance(piece, RhombusPiece):
                    for (i, triangle) in enumerate([piece.north_piece(), piece.south_piece()]):
                        if i == 0:
                            s += tikzcmd(triangle.color(), k, d, i, *triangle.border())
                        else:
                            s += tikzcmd(triangle.color(), k, d, i, "", "", "")
                else:
                    color = piece.color()
                    s += tikzcmd(color, k, d, 0, *piece.border())

        s += r"""\end{tikzpicture}"""

        return s


class KnutsonTaoPuzzleSolver(UniqueRepresentation):
    r"""
    Return puzzle solver function used to create all puzzles with given boundary conditions.

    This class implements a generic algorithm to solve Knutson-Tao puzzles.
    An instance of this class will be callable: the arguments are the
    labels of north-east and north-west sides of the puzzle boundary; the
    output is the list of the fillings of the puzzle with the specified
    pieces.

    INPUT:

    - ``puzzle_pieces`` -- takes either a collection of puzzle pieces or
      a string indicating a pre-programmed collection of puzzle pieces:

        - ``H`` -- cohomology of the Grassmannian
        - ``HT`` -- equivariant cohomology of the Grassmannian
        - ``K`` -- K-theory
        - ``H2step`` -- cohomology of the *2-step* Grassmannian
        - ``HT2step`` -- equivariant cohomology of the *2-step* Grassmannian
        - ``BK`` -- Belkale-Kumar puzzle pieces

    - ``max_letter`` -- ``None`` or a positive integer(default: ``None``); this
      is only required for Belkale-Kumar puzzles

    EXAMPLES:

    Each puzzle piece is an edge-labelled triangle oriented in such a way
    that it has a south edge (called a *delta* piece) or a north edge
    (called a *nabla* piece). For example, the puzzle pieces corresponding
    to the cohomology of the Grassmannian are the following::

        sage: from sage.combinat.knutson_tao_puzzles import H_grassmannian_pieces
        sage: H_grassmannian_pieces()
        Nablas : [0\0/0, 0\10/1, 10\1/0, 1\0/10, 1\1/1]
        Deltas : [0/0\0, 0/1\10, 1/10\0, 1/1\1, 10/0\1]

    In the string representation, the nabla pieces are depicted as
    ``c\a/b``, where `a` is the label of the north edge, `b` is the label
    of the south-east edge, `c` is the label of the south-west edge.
    A similar string representation exists for the delta pieces.

    To create a puzzle solver, one specifies a collection of puzzle pieces::

        sage: KnutsonTaoPuzzleSolver(H_grassmannian_pieces())
        Knutson-Tao puzzle solver with pieces:
        Nablas : [0\0/0, 0\10/1, 10\1/0, 1\0/10, 1\1/1]
        Deltas : [0/0\0, 0/1\10, 1/10\0, 1/1\1, 10/0\1]

    The following shorthand to create the above puzzle solver is also supported::

        sage: KnutsonTaoPuzzleSolver('H')
        Knutson-Tao puzzle solver with pieces:
        Nablas : [0\0/0, 0\10/1, 10\1/0, 1\0/10, 1\1/1]
        Deltas : [0/0\0, 0/1\10, 1/10\0, 1/1\1, 10/0\1]

    The solver will compute all fillings of the puzzle with the given
    puzzle pieces. The user specifies the labels of north-east and
    north-west sides of the puzzle boundary and the output is a list of the
    fillings of the puzzle with the specified pieces. For example, there is
    one solution to the puzzle whose north-west and north-east edges are
    both labeled '0'::

        sage: ps = KnutsonTaoPuzzleSolver('H')
        sage: ps('0', '0')
        [{(1, 1): 0/0\0}]

    There are two solutions to the puzzle whose north-west and north-east
    edges are both labeled '0101'::

        sage: ps = KnutsonTaoPuzzleSolver('H')
        sage: solns = ps('0101', '0101')
        sage: len(solns)
        2
        sage: solns.sort(key=str)
        sage: solns
        [{(1, 1): 0/0\0,
          (1, 2): 1/\0  0\/1,
          (1, 3): 0/\0  0\/0,
          (1, 4): 1/\0  0\/1,
          (2, 2): 1/1\1,
          (2, 3): 0/\10  1\/1,
          (2, 4): 1/\1  10\/0,
          (3, 3): 1/1\1,
          (3, 4): 0/\0  1\/10,
          (4, 4): 10/0\1}, {(1, 1): 0/1\10,
          (1, 2): 1/\1  10\/0,
          (1, 3): 0/\0  1\/10,
          (1, 4): 1/\0  0\/1,
          (2, 2): 0/0\0,
          (2, 3): 10/\1  0\/0,
          (2, 4): 1/\1  1\/1,
          (3, 3): 0/0\0,
          (3, 4): 1/\0  0\/1,
          (4, 4): 1/1\1}]

    The pieces in a puzzle filling are indexed by pairs of nonnegative
    integers `(i, j)` with `1 \leq i \leq j \leq n`, where `n` is the
    length of the word labelling the triangle edge. The pieces indexed by
    `(i, i)` are the triangles along the south edge of the puzzle. ::

        sage: f = solns[0]
        sage: [f[i, i] for i in range(1,5)]
        [0/0\0, 1/1\1, 1/1\1, 10/0\1]

    The pieces indexed by `(i, j)` for `j > i` are a pair consisting of
    a delta piece and nabla piece glued together along the south edge and
    north edge, respectively (these pairs are called *rhombi*). ::

        sage: f = solns[0]
        sage: f[1, 2]
        1/\0  0\/1

    There are various methods and options to display puzzle solutions.
    A single puzzle can be displayed using the plot method of the puzzle::

        sage: ps = KnutsonTaoPuzzleSolver("H")
        sage: puzzle = ps('0101','1001')[0]
        sage: puzzle.plot()  #not tested
        sage: puzzle.plot(style='fill')  #not tested
        sage: puzzle.plot(style='edges')  #not tested

    To plot several puzzle solutions, use the plot method of the puzzle
    solver::

        sage: ps = KnutsonTaoPuzzleSolver('K')
        sage: solns = ps('0101', '0101')
        sage: ps.plot(solns)        # not tested

    The code can also generate a PDF of a puzzle (using LaTeX and *tikz*)::

        sage: latex.extra_preamble(r'''\usepackage{tikz}''')
        sage: ps = KnutsonTaoPuzzleSolver('H')
        sage: solns = ps('0101', '0101')
        sage: view(solns[0], viewer='pdf')  # not tested


    Below are examples of using each of the currently supported puzzles.

    Cohomology of the Grassmannian::

        sage: ps = KnutsonTaoPuzzleSolver("H")
        sage: solns = ps('0101', '0101')
        sage: sorted(solns, key=str)
        [{(1, 1): 0/0\0,
          (1, 2): 1/\0  0\/1,
          (1, 3): 0/\0  0\/0,
          (1, 4): 1/\0  0\/1,
          (2, 2): 1/1\1,
          (2, 3): 0/\10  1\/1,
          (2, 4): 1/\1  10\/0,
          (3, 3): 1/1\1,
          (3, 4): 0/\0  1\/10,
          (4, 4): 10/0\1}, {(1, 1): 0/1\10,
          (1, 2): 1/\1  10\/0,
          (1, 3): 0/\0  1\/10,
          (1, 4): 1/\0  0\/1,
          (2, 2): 0/0\0,
          (2, 3): 10/\1  0\/0,
          (2, 4): 1/\1  1\/1,
          (3, 3): 0/0\0,
          (3, 4): 1/\0  0\/1,
          (4, 4): 1/1\1}]

    Equivariant puzzles::

        sage: ps = KnutsonTaoPuzzleSolver("HT")
        sage: solns = ps('0101', '0101')
        sage: sorted(solns, key=str)
        [{(1, 1): 0/0\0,
          (1, 2): 1/\0  0\/1,
          (1, 3): 0/\0  0\/0,
          (1, 4): 1/\0  0\/1,
          (2, 2): 1/1\1,
          (2, 3): 0/\1  1\/0,
          (2, 4): 1/\1  1\/1,
          (3, 3): 0/0\0,
          (3, 4): 1/\0  0\/1,
          (4, 4): 1/1\1}, {(1, 1): 0/0\0,
          (1, 2): 1/\0  0\/1,
          (1, 3): 0/\0  0\/0,
          (1, 4): 1/\0  0\/1,
          (2, 2): 1/1\1,
          (2, 3): 0/\10  1\/1,
          (2, 4): 1/\1  10\/0,
          (3, 3): 1/1\1,
          (3, 4): 0/\0  1\/10,
          (4, 4): 10/0\1}, {(1, 1): 0/1\10,
          (1, 2): 1/\1  10\/0,
          (1, 3): 0/\0  1\/10,
          (1, 4): 1/\0  0\/1,
          (2, 2): 0/0\0,
          (2, 3): 10/\1  0\/0,
          (2, 4): 1/\1  1\/1,
          (3, 3): 0/0\0,
          (3, 4): 1/\0  0\/1,
          (4, 4): 1/1\1}]

    K-Theory puzzles::

        sage: ps = KnutsonTaoPuzzleSolver("K")
        sage: solns = ps('0101', '0101')
        sage: sorted(solns, key=str)
        [{(1, 1): 0/0\0,
          (1, 2): 1/\0  0\/1,
          (1, 3): 0/\0  0\/0,
          (1, 4): 1/\0  0\/1,
          (2, 2): 1/1\1,
          (2, 3): 0/\10  1\/1,
          (2, 4): 1/\1  10\/0,
          (3, 3): 1/1\1,
          (3, 4): 0/\0  1\/10,
          (4, 4): 10/0\1}, {(1, 1): 0/1\10,
          (1, 2): 1/\1  10\/0,
          (1, 3): 0/\0  1\/10,
          (1, 4): 1/\0  0\/1,
          (2, 2): 0/0\0,
          (2, 3): 10/\1  0\/0,
          (2, 4): 1/\1  1\/1,
          (3, 3): 0/0\0,
          (3, 4): 1/\0  0\/1,
          (4, 4): 1/1\1}, {(1, 1): 0/1\10,
          (1, 2): 1/\1  10\/0,
          (1, 3): 0/\0  1\/K,
          (1, 4): 1/\0  0\/1,
          (2, 2): 0/0\0,
          (2, 3): K/\K  0\/1,
          (2, 4): 1/\1  K\/0,
          (3, 3): 1/1\1,
          (3, 4): 0/\0  1\/10,
          (4, 4): 10/0\1}]

    Two-step puzzles::

        sage: ps = KnutsonTaoPuzzleSolver("H2step")
        sage: solns = ps('01201', '01021')
        sage: sorted(solns, key=str)
        [{(1, 1): 0/0\0,
          (1, 2): 1/\0  0\/1,
          (1, 3): 2/\0  0\/2,
          (1, 4): 0/\0  0\/0,
          (1, 5): 1/\0  0\/1,
          (2, 2): 1/2\21,
          (2, 3): 2/\2  21\/1,
          (2, 4): 0/\10  2\/21,
          (2, 5): 1/\1  10\/0,
          (3, 3): 1/1\1,
          (3, 4): 21/\2  1\/1,
          (3, 5): 0/\0  2\/20,
          (4, 4): 1/1\1,
          (4, 5): 20/\2  1\/10,
          (5, 5): 10/0\1}, {(1, 1): 0/1\10,
          (1, 2): 1/\1  10\/0,
          (1, 3): 2/\1  1\/2,
          (1, 4): 0/\0  1\/10,
          (1, 5): 1/\0  0\/1,
          (2, 2): 0/2\20,
          (2, 3): 2/\2  20\/0,
          (2, 4): 10/\1  2\/20,
          (2, 5): 1/\1  1\/1,
          (3, 3): 0/0\0,
          (3, 4): 20/\2  0\/0,
          (3, 5): 1/\0  2\/2(10),
          (4, 4): 0/0\0,
          (4, 5): 2(10)/\2  0\/1,
          (5, 5): 1/1\1}, {(1, 1): 0/2\20,
          (1, 2): 1/\21  20\/0,
          (1, 3): 2/\2  21\/1,
          (1, 4): 0/\0  2\/20,
          (1, 5): 1/\0  0\/1,
          (2, 2): 0/0\0,
          (2, 3): 1/\0  0\/1,
          (2, 4): 20/\2  0\/0,
          (2, 5): 1/\1  2\/21,
          (3, 3): 1/1\1,
          (3, 4): 0/\0  1\/10,
          (3, 5): 21/\0  0\/21,
          (4, 4): 10/0\1,
          (4, 5): 21/\2  1\/1,
          (5, 5): 1/1\1}]

    Two-step equivariant puzzles::

        sage: ps = KnutsonTaoPuzzleSolver("HT2step")
        sage: solns = ps('10212', '12012')
        sage: sorted(solns, key=str)
        [{(1, 1): 1/1\1,
          (1, 2): 0/\(21)0  1\/2,
          (1, 3): 2/\1  (21)0\/0,
          (1, 4): 1/\1  1\/1,
          (1, 5): 2/\1  1\/2,
          (2, 2): 2/2\2,
          (2, 3): 0/\2  2\/0,
          (2, 4): 1/\2  2\/1,
          (2, 5): 2/\2  2\/2,
          (3, 3): 0/0\0,
          (3, 4): 1/\0  0\/1,
          (3, 5): 2/\0  0\/2,
          (4, 4): 1/1\1,
          (4, 5): 2/\1  1\/2,
          (5, 5): 2/2\2}, {(1, 1): 1/1\1,
          (1, 2): 0/\(21)0  1\/2,
          (1, 3): 2/\1  (21)0\/0,
          (1, 4): 1/\1  1\/1,
          (1, 5): 2/\1  1\/2,
          (2, 2): 2/2\2,
          (2, 3): 0/\2  2\/0,
          (2, 4): 1/\21  2\/2,
          (2, 5): 2/\2  21\/1,
          (3, 3): 0/0\0,
          (3, 4): 2/\0  0\/2,
          (3, 5): 1/\0  0\/1,
          (4, 4): 2/2\2,
          (4, 5): 1/\1  2\/21,
          (5, 5): 21/1\2}, {(1, 1): 1/1\1,
          (1, 2): 0/\(21)0  1\/2,
          (1, 3): 2/\1  (21)0\/0,
          (1, 4): 1/\1  1\/1,
          (1, 5): 2/\1  1\/2,
          (2, 2): 2/2\2,
          (2, 3): 0/\20  2\/2,
          (2, 4): 1/\21  20\/0,
          (2, 5): 2/\2  21\/1,
          (3, 3): 2/2\2,
          (3, 4): 0/\0  2\/20,
          (3, 5): 1/\0  0\/1,
          (4, 4): 20/0\2,
          (4, 5): 1/\1  2\/21,
          (5, 5): 21/1\2}, {(1, 1): 1/1\1,
          (1, 2): 0/\1  1\/0,
          (1, 3): 2/\1  1\/2,
          (1, 4): 1/\1  1\/1,
          (1, 5): 2/\1  1\/2,
          (2, 2): 0/2\20,
          (2, 3): 2/\2  20\/0,
          (2, 4): 1/\2  2\/1,
          (2, 5): 2/\2  2\/2,
          (3, 3): 0/0\0,
          (3, 4): 1/\0  0\/1,
          (3, 5): 2/\0  0\/2,
          (4, 4): 1/1\1,
          (4, 5): 2/\1  1\/2,
          (5, 5): 2/2\2}, {(1, 1): 1/1\1,
          (1, 2): 0/\1  1\/0,
          (1, 3): 2/\1  1\/2,
          (1, 4): 1/\1  1\/1,
          (1, 5): 2/\1  1\/2,
          (2, 2): 0/2\20,
          (2, 3): 2/\2  20\/0,
          (2, 4): 1/\21  2\/2,
          (2, 5): 2/\2  21\/1,
          (3, 3): 0/0\0,
          (3, 4): 2/\0  0\/2,
          (3, 5): 1/\0  0\/1,
          (4, 4): 2/2\2,
          (4, 5): 1/\1  2\/21,
          (5, 5): 21/1\2}, {(1, 1): 1/1\1,
          (1, 2): 0/\10  1\/1,
          (1, 3): 2/\10  10\/2,
          (1, 4): 1/\1  10\/0,
          (1, 5): 2/\1  1\/2,
          (2, 2): 1/2\21,
          (2, 3): 2/\2  21\/1,
          (2, 4): 0/\2  2\/0,
          (2, 5): 2/\2  2\/2,
          (3, 3): 1/1\1,
          (3, 4): 0/\0  1\/10,
          (3, 5): 2/\0  0\/2,
          (4, 4): 10/0\1,
          (4, 5): 2/\1  1\/2,
          (5, 5): 2/2\2}, {(1, 1): 1/1\1,
          (1, 2): 0/\10  1\/1,
          (1, 3): 2/\10  10\/2,
          (1, 4): 1/\1  10\/0,
          (1, 5): 2/\1  1\/2,
          (2, 2): 1/2\21,
          (2, 3): 2/\2  21\/1,
          (2, 4): 0/\20  2\/2,
          (2, 5): 2/\2  20\/0,
          (3, 3): 1/1\1,
          (3, 4): 2/\1  1\/2,
          (3, 5): 0/\0  1\/10,
          (4, 4): 2/2\2,
          (4, 5): 10/\1  2\/20,
          (5, 5): 20/0\2}, {(1, 1): 1/2\21,
          (1, 2): 0/\20  21\/1,
          (1, 3): 2/\2  20\/0,
          (1, 4): 1/\1  2\/21,
          (1, 5): 2/\1  1\/2,
          (2, 2): 1/1\1,
          (2, 3): 0/\1  1\/0,
          (2, 4): 21/\2  1\/1,
          (2, 5): 2/\2  2\/2,
          (3, 3): 0/0\0,
          (3, 4): 1/\0  0\/1,
          (3, 5): 2/\0  0\/2,
          (4, 4): 1/1\1,
          (4, 5): 2/\1  1\/2,
          (5, 5): 2/2\2}, {(1, 1): 1/2\21,
          (1, 2): 0/\20  21\/1,
          (1, 3): 2/\2  20\/0,
          (1, 4): 1/\1  2\/21,
          (1, 5): 2/\1  1\/2,
          (2, 2): 1/1\1,
          (2, 3): 0/\10  1\/1,
          (2, 4): 21/\2  10\/0,
          (2, 5): 2/\2  2\/2,
          (3, 3): 1/1\1,
          (3, 4): 0/\0  1\/10,
          (3, 5): 2/\0  0\/2,
          (4, 4): 10/0\1,
          (4, 5): 2/\1  1\/2,
          (5, 5): 2/2\2}, {(1, 1): 1/2\21,
          (1, 2): 0/\21  21\/0,
          (1, 3): 2/\2  21\/1,
          (1, 4): 1/\1  2\/21,
          (1, 5): 2/\1  1\/2,
          (2, 2): 0/1\10,
          (2, 3): 1/\1  10\/0,
          (2, 4): 21/\2  1\/1,
          (2, 5): 2/\2  2\/2,
          (3, 3): 0/0\0,
          (3, 4): 1/\0  0\/1,
          (3, 5): 2/\0  0\/2,
          (4, 4): 1/1\1,
          (4, 5): 2/\1  1\/2,
          (5, 5): 2/2\2}]


    Belkale-Kumar puzzles (the following example is Figure 2 of [KnutsonPurbhoo10]_)::

        sage: ps = KnutsonTaoPuzzleSolver('BK', 3)
        sage: solns = ps('12132', '23112')
        sage: len(solns)
        1
        sage: solns[0].south_labels()
        ('3', '2', '1', '2', '1')
        sage: solns
        [{(1, 1): 1/3\3(1),
          (1, 2): 2/\3(2)  3(1)\/1,
          (1, 3): 1/\3(1)  3(2)\/2,
          (1, 4): 3/\3  3(1)\/1,
          (1, 5): 2/\2  3\/3(2),
          (2, 2): 1/2\2(1),
          (2, 3): 2/\2  2(1)\/1,
          (2, 4): 1/\2(1)  2\/2,
          (2, 5): 3(2)/\3  2(1)\/1,
          (3, 3): 1/1\1,
          (3, 4): 2/\1  1\/2,
          (3, 5): 1/\1  1\/1,
          (4, 4): 2/2\2,
          (4, 5): 1/\1  2\/2(1),
          (5, 5): 2(1)/1\2}]
    """

    def __init__(self, puzzle_pieces):
        r"""
        Knutson-Tao puzzle solver.

        TESTS:

        Check that UniqueRepresentation works::

            sage: from sage.combinat.knutson_tao_puzzles import KnutsonTaoPuzzleSolver, H_grassmannian_pieces
            sage: ps = KnutsonTaoPuzzleSolver(H_grassmannian_pieces())
            sage: qs = KnutsonTaoPuzzleSolver("H")
            sage: ps
            Knutson-Tao puzzle solver with pieces:
            Nablas : [0\0/0, 0\10/1, 10\1/0, 1\0/10, 1\1/1]
            Deltas : [0/0\0, 0/1\10, 1/10\0, 1/1\1, 10/0\1]
            sage: qs
            Knutson-Tao puzzle solver with pieces:
            Nablas : [0\0/0, 0\10/1, 10\1/0, 1\0/10, 1\1/1]
            Deltas : [0/0\0, 0/1\10, 1/10\0, 1/1\1, 10/0\1]
            sage: ps == qs
            True
        """
        self._puzzle_pieces = puzzle_pieces
        self._rhombus_pieces = tuple(puzzle_pieces.rhombus_pieces())
        self._bottom_deltas = tuple(puzzle_pieces.boundary_deltas())

    @staticmethod
    def __classcall_private__(cls, puzzle_pieces, max_letter=None):
        r"""
        TESTS::

            sage: from sage.combinat.knutson_tao_puzzles import *
            sage: KnutsonTaoPuzzleSolver(H_grassmannian_pieces()) == KnutsonTaoPuzzleSolver("H") # indirect doctest
            True
            sage: KnutsonTaoPuzzleSolver(HT_grassmannian_pieces()) == KnutsonTaoPuzzleSolver("HT")
            True
            sage: KnutsonTaoPuzzleSolver(K_grassmannian_pieces()) == KnutsonTaoPuzzleSolver("K")
            True
            sage: KnutsonTaoPuzzleSolver(H_two_step_pieces()) == KnutsonTaoPuzzleSolver("H2step")
            True
            sage: KnutsonTaoPuzzleSolver(HT_two_step_pieces()) == KnutsonTaoPuzzleSolver("HT2step")
            True
            sage: KnutsonTaoPuzzleSolver(BK_pieces(3)) == KnutsonTaoPuzzleSolver("BK",3)
            True
        """
        if isinstance(puzzle_pieces, str):
            if puzzle_pieces == "H":
                puzzle_pieces = H_grassmannian_pieces()
            elif puzzle_pieces == "HT":
                puzzle_pieces = HT_grassmannian_pieces()
            elif puzzle_pieces == "K":
                puzzle_pieces = K_grassmannian_pieces()
            elif puzzle_pieces == "H2step":
                puzzle_pieces = H_two_step_pieces()
            elif puzzle_pieces == "HT2step":
                puzzle_pieces = HT_two_step_pieces()
            elif puzzle_pieces == "BK":
                if max_letter is not None:
                    puzzle_pieces = BK_pieces(max_letter)
                else:
                    raise ValueError("max_letter needs to be specified")
        return super().__classcall__(cls, puzzle_pieces)

    def __call__(self, lamda, mu, algorithm='strips'):
        r"""
        TESTS::

            sage: from sage.combinat.knutson_tao_puzzles import KnutsonTaoPuzzleSolver
            sage: ps = KnutsonTaoPuzzleSolver("H")
            sage: ps('0101','1001')
            [{(1, 1): 0/1\10,
              (1, 2): 1/\1  10\/0,
              (1, 3): 0/\10  1\/1,
              (1, 4): 1/\1  10\/0,
              (2, 2): 0/0\0,
              (2, 3): 1/\0  0\/1,
              (2, 4): 0/\0  0\/0,
              (3, 3): 1/1\1,
              (3, 4): 0/\0  1\/10,
              (4, 4): 10/0\1}]
            sage: ps('0101','1001',algorithm='pieces')
            [{(1, 1): 0/1\10,
              (1, 2): 1/\1  10\/0,
              (1, 3): 0/\10  1\/1,
              (1, 4): 1/\1  10\/0,
              (2, 2): 0/0\0,
              (2, 3): 1/\0  0\/1,
              (2, 4): 0/\0  0\/0,
              (3, 3): 1/1\1,
              (3, 4): 0/\0  1\/10,
              (4, 4): 10/0\1}]
        """
        lamda, mu = tuple(lamda), tuple(mu)
        if algorithm == 'pieces':
            return list(self._fill_puzzle_by_pieces(lamda, mu))
        elif algorithm == 'strips':
            return list(self._fill_puzzle_by_strips(lamda, mu))

    solutions = __call__

    def __repr__(self) -> str:
        r"""
        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import KnutsonTaoPuzzleSolver
            sage: KnutsonTaoPuzzleSolver('H')
            Knutson-Tao puzzle solver with pieces:
            Nablas : [0\0/0, 0\10/1, 10\1/0, 1\0/10, 1\1/1]
            Deltas : [0/0\0, 0/1\10, 1/10\0, 1/1\1, 10/0\1]
        """
        return "Knutson-Tao puzzle solver with pieces:\n%s" % self._puzzle_pieces

    def puzzle_pieces(self):
        r"""
        The puzzle pieces used for filling in the puzzles.

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import KnutsonTaoPuzzleSolver
            sage: ps = KnutsonTaoPuzzleSolver('H')
            sage: ps.puzzle_pieces()
            Nablas : [0\0/0, 0\10/1, 10\1/0, 1\0/10, 1\1/1]
            Deltas : [0/0\0, 0/1\10, 1/10\0, 1/1\1, 10/0\1]
        """
        return self._puzzle_pieces

    def _fill_piece(self, nw_label, ne_label, pieces) -> list[PuzzlePiece]:
        r"""
        Fillings of a piece.

        INPUT:

        - ``nw_label``, ``nw_label`` -- label
        - ``pieces`` -- puzzle pieces used for the filling

        OUTPUT: list of the fillings

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import KnutsonTaoPuzzleSolver
            sage: ps = KnutsonTaoPuzzleSolver('H')
            sage: ps._fill_piece('0', '0', ps._bottom_deltas)
            [0/0\0]
        """
        output = []
        for piece in pieces:
            if (piece['north_west'] == nw_label and
                    piece['north_east'] == ne_label):
                output.append(piece)
        return output

    @cached_method
    def _fill_strip(self, nw_labels, ne_label, pieces, final_pieces=None):
        r"""
        Fillings of a strip of height 1.

        INPUT:

        - ``nw_labels`` -- tuple of labels
        - ``nw_label`` -- label
        - ``pieces`` -- puzzle pieces used for the filling
        - ``final_pieces`` -- pieces used for the last piece to be filled in

        OUTPUT: list of lists of the fillings

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import KnutsonTaoPuzzleSolver
            sage: ps = KnutsonTaoPuzzleSolver('H')
            sage: ps._fill_strip(('0',), '0', ps._rhombus_pieces, ps._bottom_deltas)
            [[0/0\0]]
            sage: ps._fill_strip(('0','0'), '0', ps._rhombus_pieces, ps._bottom_deltas)
            [[0/\0  0\/0, 0/0\0]]
            sage: sorted(ps._fill_strip(('0',), '0', ps._rhombus_pieces), key=str)
            [[0/\0  0\/0], [0/\0  1\/10]]
            sage: sorted(ps._fill_strip(('0','1'), '0', ps._rhombus_pieces), key =str)
            [[1/\0  0\/1, 0/\0  0\/0], [1/\0  0\/1, 0/\0  1\/10]]

        TESTS::

            sage: from sage.combinat.knutson_tao_puzzles import KnutsonTaoPuzzleSolver
            sage: ps = KnutsonTaoPuzzleSolver('H')
            sage: ps._fill_strip(('0',), 'goo', ps._rhombus_pieces)
            []
        """
        if final_pieces is None:
            final_pieces = pieces

        output = []
        if len(nw_labels) == 1:
            X = self._fill_piece(nw_labels[0], ne_label, final_pieces)
            if X:
                output = [[x] for x in X]
        else:
            partial_fillings = self._fill_strip(nw_labels[1:], ne_label, pieces)
            for partial_filling in partial_fillings:
                ne_label = partial_filling[-1]['south_west']
                for piece in self._fill_piece(nw_labels[0], ne_label, final_pieces):
                    output.append(partial_filling + [piece])
        return output

    def _fill_puzzle_by_pieces(self, lamda, mu):
        r"""
        Fill puzzle pieces for given outer labels ``lambda`` and ``mu``.

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import KnutsonTaoPuzzleSolver
            sage: ps = KnutsonTaoPuzzleSolver('H')
            sage: list(ps._fill_puzzle_by_pieces('0', '0'))
            [{(1, 1): 0/0\0}]
        """
        queue = [PuzzleFilling(lamda, mu)]
        while queue:
            PP = queue.pop()
            ne_label = PP.north_east_label_of_kink()
            nw_label = PP.north_west_label_of_kink()
            if PP.is_in_south_edge():
                pieces = self._bottom_deltas
            else:
                pieces = self._rhombus_pieces
            for piece in self._fill_piece(nw_label, ne_label, pieces):
                PPcopy = PP.copy()
                PPcopy.add_piece(piece)
                if PPcopy.is_completed():
                    yield PPcopy
                else:
                    queue.append(PPcopy)

    def _fill_puzzle_by_strips(self, lamda, mu):
        r"""
        Fill puzzle pieces by strips for given outer labels ``lambda`` and ``mu``.

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import KnutsonTaoPuzzleSolver
            sage: ps = KnutsonTaoPuzzleSolver('H')
            sage: list(ps._fill_puzzle_by_strips('0', '0'))
            [{(1, 1): 0/0\0}]
            sage: list(ps._fill_puzzle_by_strips('01', '01'))
            [{(1, 1): 0/0\0, (1, 2): 1/\0  0\/1, (2, 2): 1/1\1}]
        """
        queue = [PuzzleFilling(lamda, mu)]
        while queue:
            PP = queue.pop()
            i, _ = PP.kink_coordinates()

            # grab nw labels
            if i == 1:
                nw_labels = PP._nw_labels
            else:
                nw_labels = tuple(PP._squares[i - 1, k]['south_east']
                                  for k in range(i, len(lamda) + 1))

            # grab ne labels
            ne_label = PP._ne_labels[i - 1]

            deltas = self._bottom_deltas
            rhombi = self._rhombus_pieces
            for row in self._fill_strip(nw_labels, ne_label, rhombi, deltas):
                PPcopy = PP.copy()
                PPcopy.add_pieces(row)
                if PPcopy.is_completed():
                    yield PPcopy
                else:
                    queue.append(PPcopy)

    def plot(self, puzzles):
        r"""
        Return plot of puzzles.

        INPUT:

        - ``puzzles`` -- list of puzzles

        EXAMPLES::

            sage: from sage.combinat.knutson_tao_puzzles import KnutsonTaoPuzzleSolver
            sage: ps = KnutsonTaoPuzzleSolver('K')
            sage: solns = ps('0101', '0101')
            sage: ps.plot(solns)        # not tested
        """
        g = [p.plot() for p in puzzles]
        m = len([gg.axes(False) for gg in g])
        return graphics_array(g, (m + 3) / 4, 4)

    def structure_constants(self, lamda, mu, nu=None):
        r"""
        Compute cohomology structure coefficients from puzzles.

        INPUT:

        - ``pieces`` -- puzzle pieces to be used
        - ``lambda``, ``mu`` -- edge labels of puzzle for northwest and north east side
        - ``nu`` -- (default: ``None``) if ``nu`` is not specified a dictionary is returned with
          the structure coefficients corresponding to all south labels; if ``nu`` is given, only
          the coefficients with the specified label is returned

        OUTPUT: dictionary

        EXAMPLES:

        Note: In order to standardize the output of the following examples,
        we output a sorted list of items from the dictionary instead of the
        dictionary itself.

        Grassmannian cohomology::

            sage: ps = KnutsonTaoPuzzleSolver('H')
            sage: cp = ps.structure_constants('0101', '0101')
            sage: sorted(cp.items(), key=str)
            [(('0', '1', '1', '0'), 1), (('1', '0', '0', '1'), 1)]
            sage: ps.structure_constants('001001', '001010', '010100')
            1

        Equivariant cohomology::

            sage: ps = KnutsonTaoPuzzleSolver('HT')
            sage: cp = ps.structure_constants('0101', '0101')
            sage: sorted(cp.items(), key=str)
            [(('0', '1', '0', '1'), y2 - y3),
            (('0', '1', '1', '0'), 1),
            (('1', '0', '0', '1'), 1)]

        K-theory::

            sage: ps = KnutsonTaoPuzzleSolver('K')
            sage: cp = ps.structure_constants('0101', '0101')
            sage: sorted(cp.items(), key=str)
            [(('0', '1', '1', '0'), 1), (('1', '0', '0', '1'), 1), (('1', '0', '1', '0'), -1)]

        Two-step::

            sage: ps = KnutsonTaoPuzzleSolver('H2step')
            sage: cp = ps.structure_constants('01122', '01122')
            sage: sorted(cp.items(), key=str)
            [(('0', '1', '1', '2', '2'), 1)]
            sage: cp = ps.structure_constants('01201', '01021')
            sage: sorted(cp.items(), key=str)
            [(('0', '2', '1', '1', '0'), 1),
             (('1', '2', '0', '0', '1'), 1),
             (('2', '0', '1', '0', '1'), 1)]

        Two-step equivariant::

            sage: ps = KnutsonTaoPuzzleSolver('HT2step')
            sage: cp = ps.structure_constants('10212', '12012')
            sage: sorted(cp.items(), key=str)
            [(('1', '2', '0', '1', '2'), y1*y2 - y2*y3 - y1*y4 + y3*y4),
             (('1', '2', '0', '2', '1'), y1 - y3),
             (('1', '2', '1', '0', '2'), y2 - y4),
             (('1', '2', '1', '2', '0'), 1),
             (('1', '2', '2', '0', '1'), 1),
             (('2', '1', '0', '1', '2'), y1 - y3),
             (('2', '1', '1', '0', '2'), 1)]
        """
        from collections import defaultdict
        R = PolynomialRing(Integers(), 'y', len(lamda) + 1)
        z = defaultdict(R.zero)
        for p in self(lamda, mu):
            z[p.south_labels()] += p.contribution()
        if nu is None:
            return dict(z)
        else:
            return z[tuple(nu)]
