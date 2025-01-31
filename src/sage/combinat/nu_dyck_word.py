# sage.doctest: needs sage.combinat sage.modules
r"""
`\nu`-Dyck words

A class of the `\nu`-Dyck word, see [PRV2017]_ for details.

This file is based off the class
:func:`DyckWords<sage.combinat.dyck_word.DyckWord>` written by Mike Hansen, Dan
Drake, Florent Hivert, Christian Stump, Mike Zabrocki, Jean--Baptiste Priez
and Travis Scrimshaw.

AUTHORS:

- Aram Dermenjian (2020-09-26)
"""
# ****************************************************************************
#       Copyright (C) 2020 Aram Dermenjian <aram.dermenjian@gmail.com>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from __future__ import annotations
from sage.structure.element import Element
from sage.rings.integer import Integer
from sage.combinat.combinat import CombinatorialElement
from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet
from sage.structure.global_options import GlobalOptions
from sage.structure.parent import Parent
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.misc.latex import latex

from sage.combinat.words.paths import WordPaths_north_east
from sage.combinat.words.paths import FiniteWordPath_north_east

ndw_open_symbol = 1
ndw_close_symbol = 0


def update_ndw_symbols(os, cs):
    r"""
    A way to alter the open and close symbols from sage.

    INPUT:

    - ``os`` -- the open symbol
    - ``cs`` -- the close symbol

    EXAMPLES::

        sage: from sage.combinat.nu_dyck_word import update_ndw_symbols
        sage: update_ndw_symbols(0,1)
        sage: dw = NuDyckWord('0101001','0110010'); dw
        [0, 1, 0, 1, 0, 0, 1]

        sage: dw = NuDyckWord('1010110','1001101'); dw
        Traceback (most recent call last):
        ...
        ValueError: invalid nu-Dyck word
        sage: update_ndw_symbols(1,0)
    """
    global ndw_open_symbol
    global ndw_close_symbol
    ndw_open_symbol = os
    ndw_close_symbol = cs


def replace_dyck_char(x):
    r"""
    A map sending an opening character (``'1'``, ``'N'``, and ``'('``) to
    ``ndw_open_symbol``, a closing character (``'0'``, ``'E'``, and ``')'``) to
    ``ndw_close_symbol``, and raising an error on any input other than one of
    the opening or closing characters.

    This is the inverse map of :func:`replace_dyck_symbol`.

    INPUT:

    - ``x`` -- string; a ``'1'``, ``'0'``, ``'N'``, ``'E'``, ``'('`` or ``')'``

    OUTPUT:

    - If ``x`` is an opening character, replace ``x`` with the
      constant ``ndw_open_symbol``.

    - If ``x`` is a closing character, replace ``x`` with the
      constant ``ndw_close_symbol``.

    - Raise a :exc:`ValueError` if ``x`` is neither an opening nor a
      closing character.

    .. SEEALSO:: :func:`replace_dyck_symbol`

    EXAMPLES::

        sage: from sage.combinat.nu_dyck_word import replace_dyck_char
        sage: replace_dyck_char('(')
        1
        sage: replace_dyck_char(')')
        0
        sage: replace_dyck_char(1)
        Traceback (most recent call last):
        ...
        ValueError
    """
    if x == '(' or x == 'N' or x == str(ndw_open_symbol):
        return ndw_open_symbol
    if x == ')' or x == 'E' or x == str(ndw_close_symbol):
        return ndw_close_symbol
    raise ValueError


def replace_dyck_symbol(x, open_char='N', close_char='E') -> str:
    r"""
    A map sending ``ndw_open_symbol`` to ``open_char`` and ``ndw_close_symbol``
    to ``close_char``, and raising an error on any input other than
    ``ndw_open_symbol`` and ``ndw_close_symbol``. The values of the constants
    ``ndw_open_symbol`` and ``ndw_close_symbol`` are subject to change.

    This is the inverse map of :func:`replace_dyck_char`.

    INPUT:

    - ``x`` -- either ``ndw_open_symbol`` or ``ndw_close_symbol``

    - ``open_char`` -- string (optional) default ``'N'``

    - ``close_char`` -- string (optional) default ``'E'``

    OUTPUT: if ``x`` is ``ndw_open_symbol``, replace ``x`` with ``open_char``

    - If ``x`` is ``ndw_close_symbol``, replace ``x`` with ``close_char``.

    - If ``x`` is neither ``ndw_open_symbol`` nor ``ndw_close_symbol``, a
      :exc:`ValueError` is raised.

    .. SEEALSO:: :func:`replace_dyck_char`

    EXAMPLES::

        sage: from sage.combinat.nu_dyck_word import replace_dyck_symbol
        sage: replace_dyck_symbol(1)
        'N'
        sage: replace_dyck_symbol(0)
        'E'
        sage: replace_dyck_symbol(3)
        Traceback (most recent call last):
        ...
        ValueError
    """
    if x == ndw_open_symbol:
        return open_char
    if x == ndw_close_symbol:
        return close_char
    raise ValueError


class NuDyckWord(CombinatorialElement):
    r"""
    A `\nu`-Dyck word.

    Given a lattice path `\nu` in the `\ZZ^2` grid starting at the origin
    `(0,0)` consisting of North `N = (0,1)` and East `E = (1,0)` steps, a
    `\nu`-Dyck path is a lattice path in the `\ZZ^2` grid starting at the
    origin `(0,0)` and ending at the same coordinate as `\nu` such that it is
    weakly above `\nu`. A `\nu`-Dyck word is the representation of a
    `\nu`-Dyck path where a North step is represented by a 1 and an East step
    is represented by a 0.

    INPUT:

    - ``k1`` -- a path for the `\nu`-Dyck word

    - ``k2`` -- a path for `\nu`

    EXAMPLES::

        sage: dw = NuDyckWord([1,0,1,0],[1,0,0,1]); dw
        [1, 0, 1, 0]
        sage: print(dw)
        NENE
        sage: dw.height()
        2

        sage: dw = NuDyckWord('1010',[1,0,0,1]); dw
        [1, 0, 1, 0]

        sage: dw = NuDyckWord('NENE',[1,0,0,1]); dw
        [1, 0, 1, 0]

        sage: NuDyckWord([1,0,1,0],[1,0,0,1]).pretty_print()
           __
         _|x
        | . .

        sage: from sage.combinat.nu_dyck_word import update_ndw_symbols
        sage: update_ndw_symbols(0,1)
        sage: dw = NuDyckWord('0101001','0110010'); dw
        [0, 1, 0, 1, 0, 0, 1]
        sage: dw.pp()
             __
            |x
           _| .
         _|x  .
        | . . .
        sage: update_ndw_symbols(1,0)
    """
    @staticmethod
    def __classcall_private__(cls, dw=None, nu=None, **kwargs):
        """
        Return an element with the appropriate parent.

        EXAMPLES::

            sage: NuDyckWord('110100','101010')
            [1, 1, 0, 1, 0, 0]
            sage: NuDyckWord('010','010')
            [0, 1, 0]
        """
        # if dw is none, then we might have a normal Dyck word
        if dw is None:
            from sage.combinat.dyck_word import DyckWord
            return DyckWord(dw, kwargs)

        if isinstance(dw, NuDyckWord):
            return dw

        if nu is None:
            raise ValueError("nu required")

        dw = to_word_path(dw)
        nu = to_word_path(nu)

        if path_weakly_above_other(dw, nu):
            return NuDyckWords(nu)(dw)

        raise ValueError("invalid nu-Dyck word")

    def __init__(self, parent, dw, latex_options=None):
        """
        Initialize a nu-Dyck word.

        EXAMPLES::

            sage: NuDyckWord('010', '010')
            [0, 1, 0]
            sage: NuDyckWord('110100','101010')
            [1, 1, 0, 1, 0, 0]
        """
        Element.__init__(self, parent)
        self._path = to_word_path(dw)
        self._list = list(self._path)

        if parent is None:
            raise ValueError("need parent object")

        self._nu = parent._nu

        if latex_options is None:
            latex_options = {}
        self._latex_options = dict(latex_options)

    def __eq__(self, other):
        """
        Return if two paths are equal.

        EXAMPLES::

            sage: u = NuDyckWord('010','010')
            sage: w = NuDyckWord('110100','101010')
            sage: w == w
            True
            sage: u == w
            False
            sage: u == 4
            False
        """
        if not isinstance(other, NuDyckWord):
            return False
        return self._path == other._path and self._nu == other._nu

    def __neq__(self, other):
        """
        Return if two paths are not equal.

        EXAMPLES::

            sage: u = NuDyckWord('010','010')
            sage: w = NuDyckWord('110100','101010')
            sage: w != w
            False
            sage: u != w
            True
            sage: u != 4
            True
        """
        return not self.__eq__(other)

    def __le__(self, other):
        """
        Return if one path is included in another.

        EXAMPLES::

            sage: ND1 = NuDyckWord('101', '011')
            sage: ND2 = NuDyckWord('110', '011')
            sage: ND3 = NuDyckWord('011', '011')
            sage: ND1 <= ND1
            True
            sage: ND1 <= ND2
            True
            sage: ND2 <= ND1
            False
            sage: ND3 <= ND2
            True
            sage: ND3 <= ND1
            True
        """
        if self._nu == other._nu:
            return path_weakly_above_other(other._path, self._path)
        return False

    def __lt__(self, other):
        """
        Return if one path is strictly included in another.

        EXAMPLES::

            sage: ND1 = NuDyckWord('101', '011')
            sage: ND2 = NuDyckWord('110', '011')
            sage: ND3 = NuDyckWord('011', '011')
            sage: ND1 < ND1
            False
            sage: ND1 < ND2
            True
            sage: ND2 < ND1
            False
            sage: ND3 < ND2
            True
            sage: ND3 < ND1
            True
        """
        return self.__le__(other) and not self.__eq__(other)

    def __ge__(self, other):
        """
        Return if one path is included in another.

        EXAMPLES::

            sage: ND1 = NuDyckWord('101', '011')
            sage: ND2 = NuDyckWord('110', '011')
            sage: ND3 = NuDyckWord('011', '011')
            sage: ND1 >= ND1
            True
            sage: ND1 >= ND2
            False
            sage: ND2 >= ND1
            True
            sage: ND1 >= ND3
            True
            sage: ND2 >= ND3
            True
        """
        if self._nu == other._nu:
            return path_weakly_above_other(self._path, other._path)
        return False

    def __gt__(self, other):
        """
        Return if one path is strictly included in another.

        EXAMPLES::

            sage: ND1 = NuDyckWord('101', '011')
            sage: ND2 = NuDyckWord('110', '011')
            sage: ND3 = NuDyckWord('011', '011')
            sage: ND1 > ND1
            False
            sage: ND1 > ND2
            False
            sage: ND2 > ND1
            True
            sage: ND1 > ND3
            True
            sage: ND2 > ND3
            True
        """
        return self.__ge__(other) and not self.__eq__(other)

    def _cache_key(self) -> tuple:
        """
        Return a cache key for ``self``.

        EXAMPLES::

            sage: u = NuDyckWord('010','010')
            sage: u._cache_key()
            (0, 1, 0, 0, 1, 0)
        """
        return tuple(self._path) + tuple(self._nu)

    def __hash__(self) -> int:
        """
        Return a hash for ``self``.

        EXAMPLES::

            sage: u = NuDyckWord('010','010')
            sage: hash(u)  # random
            -4577085166836515071
        """
        return hash(''.join(str(i) for i in self._list))

    def set_latex_options(self, D):
        r"""
        Set the latex options for use in the ``_latex_`` function.

        The default values are set in the ``__init__`` function.

        - ``color`` -- (default: black) the line color

        - ``line width`` -- (default: `2 \times` ``tikz_scale``) value
          representing the line width

        - ``nu_options`` -- (default: ``'rounded corners=1, color=red, line
          width=1'``) string to indicate what the tikz options should be for
          path of `\nu`

        - ``points_color`` -- (default: ``'black'``) str to indicate color
          points should be drawn with

        - ``show_grid`` -- boolean (default: ``True``); value to indicate if
          grid should be shown

        - ``show_nu`` -- boolean (default: ``True``); value to indicate if `\nu`
          should be shown

        - ``show_points`` -- boolean (default: ``False``); value to indicate
          if points should be shown on path

        - ``tikz_scale`` -- (default: 1) scale for use with the tikz package

        INPUT:

        - ``D`` -- dictionary with a list of latex parameters to change

        EXAMPLES::

            sage: NDW = NuDyckWord('010','010')
            sage: NDW.set_latex_options({"tikz_scale":2})
            sage: NDW.set_latex_options({"color":"blue", "show_points":True})

        .. TODO::

            This should probably be merged into NuDyckWord.options.
        """
        for opt in D:
            self._latex_options[opt] = D[opt]

    def latex_options(self) -> dict:
        r"""
        Return the latex options for use in the ``_latex_`` function as a
        dictionary.

        The default values are set using the options.

        - ``color`` -- (default: black) the line color

        - ``line width`` -- (default: 2*``tikz_scale``) value representing the
          line width

        - ``nu_options`` -- (default: ``'rounded corners=1, color=red, line
          width=1'``) string to indicate what the tikz options should be for
          path of `\nu`

        - ``points_color`` -- (default: ``'black'``) str to indicate color
          points should be drawn with

        - ``show_grid`` -- boolean (default: ``True``); value to indicate if
          grid should be shown

        - ``show_nu`` -- boolean (default: ``True``); value to indicate if `\nu`
          should be shown

        - ``show_points`` -- boolean (default: ``False``); value to indicate
          if points should be shown on path

        - ``tikz_scale`` -- (default: 1) scale for use with the tikz package

        EXAMPLES::

            sage: NDW = NuDyckWord('010','010')
            sage: NDW.latex_options()
            {'color': black,
             'line width': 2,
             'nu_options': rounded corners=1, color=red, line width=1,
             'points_color': black,
             'show_grid': True,
             'show_nu': True,
             'show_points': False,
             'tikz_scale': 1}

        .. TODO::

            This should probably be merged into NuDyckWord.options.
        """
        d = self._latex_options.copy()
        opts = self.parent().options
        if "tikz_scale" not in d:
            d["tikz_scale"] = opts.latex_tikz_scale
        if "line width" not in d:
            d["line width"] = opts.latex_line_width_scalar * d["tikz_scale"]
        if "color" not in d:
            d["color"] = opts.latex_color
        if "show_points" not in d:
            d["show_points"] = opts.latex_show_points
        if "points_color" not in d:
            d["points_color"] = opts.latex_points_color
        if "show_grid" not in d:
            d["show_grid"] = opts.latex_show_grid
        if "show_nu" not in d:
            d["show_nu"] = opts.latex_show_nu
        if "nu_options" not in d:
            d["nu_options"] = opts.latex_nu_options
        return d

    def _repr_(self) -> str:
        r"""
        Return a string representation of ``self`` depending on
        :meth:`NuDyckWords.options`.

        TESTS::

            sage: NuDyckWord('01010','00011')
            [0, 1, 0, 1, 0]
            sage: NuDyckWord('10010','00011')
            [1, 0, 0, 1, 0]
            sage: NuDyckWords.options.display="lattice"
            sage: NuDyckWord('10010','00011')
                 __
             ___|x
            |x x x

            sage: NuDyckWords.options._reset()
        """
        return self.parent().options._dispatch(self, '_repr_', 'display')

    def _repr_list(self) -> str:
        r"""
        Return a string representation of ``self`` as a list.

        TESTS::

            sage: NuDyckWord([1,1,0],[1,0,1])  # indirect doctest
            [1, 1, 0]
            sage: NuDyckWord('NNEE','NENE')
            [1, 1, 0, 0]
        """
        return str(list(self._path))

    def _repr_lattice(self, style=None, labelling=None):
        r"""
        See :meth:`pretty_print()`.

        TESTS::

            sage: n = NuDyckWord('00011001000100','00011001000100')
            sage: print(n._repr_lattice(style='N-E', labelling=[1,2,3,4]))
                             ____
                       _____| . . 4
                   ___| . . . . . 3
                  | . . . . . . . 2
            ______| . . . . . . . 1


            sage: print(NuDyckWord('100101','010011')._repr_lattice())
                 _|
             ___|x
            |x  . .

            sage: print(NuDyckWord('110010','001011')._repr_lattice())
                 __
             ___|x
            |x x x
            |x x  .
        """
        if style is None:
            style = self.parent().options.diagram_style
            if style == "grid":
                style = "N-E"

        if style == "N-E":
            path_length = self._path.length()
            height = self._path.height()
            width = self._path.width()
            if path_length == 0:
                return ".\n"

            # Handle right-hand side labels
            if labelling is None:
                labels = [" "] * height
            else:
                if len(labelling) != height:
                    raise ValueError(f"the given labelling has the wrong length: {height} needed")
                labels = [str(label) for label in labelling]
                max_length = max(len(label) for label in labels)
                labels = [lbl.rjust(max_length + 1) for lbl in labels]

            rev_path = list(self._path.reversal())
            rev_nu_path = list(self._nu.reversal())
            ts = ""

            # Grab first line
            cur_pos = rev_path.index(ndw_open_symbol)
            cur_nu_pos = rev_nu_path.index(ndw_open_symbol)
            if cur_pos > 0:
                ts += "  " * (width - cur_pos)
                ts += " _" + "__" * (cur_pos - 1)
                ts += "_\n"

            # Middle Lines
            for i in range(height - 1):
                old_pos = cur_pos
                old_nu_pos = cur_nu_pos
                cur_pos = rev_path.index(ndw_open_symbol, cur_pos + 1)
                cur_nu_pos = rev_nu_path.index(ndw_open_symbol, cur_nu_pos + 1)

                ts += "  " * (width - cur_pos + i + 1)
                if cur_pos != old_pos + 1:
                    ts += " _" + "__" * (cur_pos - old_pos - 2)
                ts += "|"
                if old_pos >= 0:
                    ts += "x " * (old_pos - old_nu_pos)
                    ts += " ." * (old_nu_pos - i)
                ts += labels[height - i - 1]
                ts += "\n"

            # Final line
            ts += "__" * (path_length - cur_pos - 1)
            ts += "|"
            ts += "x " * (cur_pos - cur_nu_pos)
            ts += " ." * (cur_nu_pos - i - 1)
            ts += labels[0]
            ts += "\n"
            return ts
        raise ValueError(f"the given style (={style}) is not valid")

    def _ascii_art_(self):
        r"""
        Return an ASCII art representation of ``self``.

        TESTS::

            sage: ascii_art(NuDyckWord('00011001000100','00011001000100'))
                             ____
                       _____| . .
                   ___| . . . . .
                  | . . . . . . .
            ______| . . . . . . .
        """
        from sage.typeset.ascii_art import AsciiArt
        rep = self.parent().options.ascii_art
        if rep == "pretty_output":
            ret = self._repr_lattice()
        return AsciiArt(ret.splitlines(), baseline=0)

    def __str__(self) -> str:
        r"""
        Return a string consisting of N and E steps corresponding to
        the `\nu`-Dyck word.

        EXAMPLES::

            sage: str(NuDyckWord('100101','010011'))
            'NEENEN'
            sage: str(NuDyckWord('101010','100110'))
            'NENENE'
        """
        return "".join(replace_dyck_symbol(let) for let in self._path)

    def pretty_print(self, style=None, labelling=None):
        r"""
        Display a NuDyckWord as a lattice path in the `\ZZ^2` grid.

        If the ``style`` is "N-E", then a cell below the diagonal is
        indicated by a period, whereas a cell below the path but above
        the diagonal is indicated by an x. If a list of labels is
        included, they are displayed along the vertical edges of the
        Dyck path.

        INPUT:

        - ``style`` -- (default: ``None``) can either be:

          - ``None`` to use the option default
          - "N-E" to show ``self`` as a path of north and east steps, or

        - ``labelling`` -- (if style is "N-E") a list of labels assigned to
          the up steps in ``self``

        - ``underpath`` -- (if style is "N-E", default: ``True``) if ``True``,
          an ``x`` to show the boxes between `\nu` and the `\nu`-Dyck Path

        EXAMPLES::

            sage: for ND in NuDyckWords('101010'): ND.pretty_print()
                 __
               _| .
             _| . .
            | . . .
                 __
             ___| .
            |x  . .
            | . . .
               ____
              |x  .
             _| . .
            | . . .
               ____
             _|x  .
            |x  . .
            | . . .
             ______
            |x x  .
            |x  . .
            | . . .

        ::

            sage: nu = [1,0,1,0,1,0,1,0,1,0,1,0]
            sage: ND = NuDyckWord([1,1,1,0,1,0,0,1,1,0,0,0],nu)
            sage: ND.pretty_print()
                   ______
                  |x x  .
               ___|x  . .
             _|x x  . . .
            |x x  . . . .
            |x  . . . . .
            | . . . . . .

        ::

            sage: NuDyckWord([1,1,0,0,1,0],[1,0,1,0,1,0]).pretty_print(
            ....: labelling=[1,3,2])
                 __
             ___| . 2
            |x  . . 3
            | . . . 1

        ::

            sage: NuDyckWord('1101110011010010001101111000110000',
            ....: '1010101010101010101010101010101010').pretty_print(
            ....: labelling=list(range(1,18)))
                                       ________
                                      |x x x  . 17
                                 _____|x x  . . 16
                                |x x x x  . . . 15
                                |x x x  . . . . 14
                                |x x  . . . . . 13
                               _|x  . . . . . . 12
                              |x  . . . . . . . 11
                         _____| . . . . . . . . 10
                     ___|x x  . . . . . . . . .  9
                   _|x x x  . . . . . . . . . .  8
                  |x x x  . . . . . . . . . . .  7
               ___|x x  . . . . . . . . . . . .  6
              |x x x  . . . . . . . . . . . . .  5
              |x x  . . . . . . . . . . . . . .  4
             _|x  . . . . . . . . . . . . . . .  3
            |x  . . . . . . . . . . . . . . . .  2
            | . . . . . . . . . . . . . . . . .  1


        ::

            sage: NuDyckWord().pretty_print()
            .
        """
        print(self._repr_lattice(style, labelling))

    pp = pretty_print

    def _latex_(self):
        r"""
        A latex representation of ``self`` using the tikzpicture package.

        EXAMPLES::

            sage: NDW = NuDyckWord('010','010')
            sage: NDW.set_latex_options({"show_points":True})
            sage: latex(NDW)
            \vcenter{\hbox{$\begin{tikzpicture}[scale=1]
              \draw[dotted] (0, 0) grid (2, 1);
              \draw[line width=2,color=black,fill=black](0, 0) circle (0.21);
              \draw[line width=2,color=black,fill=black](1, 0) circle (0.21);
              \draw[line width=2,color=black,fill=black](1, 1) circle (0.21);
              \draw[line width=2,color=black,fill=black](2, 1) circle (0.21);
              \draw[rounded corners=1, color=red, line width=1] (0, 0) -- (1, 0) -- (1, 1) -- (2, 1);
              \draw[rounded corners=1, color=black, line width=2] (0, 0) -- (1, 0) -- (1, 1) -- (2, 1);
            \end{tikzpicture}$}}
            sage: NuDyckWord('01','01')._latex_()
            '\\vcenter{\\hbox{$\\begin{tikzpicture}[scale=1]\n  \\draw[dotted] (0, 0) grid (1, 1);\n  \\draw[rounded corners=1, color=red, line width=1] (0, 0) -- (1, 0) -- (1, 1);\n  \\draw[rounded corners=1, color=black, line width=2] (0, 0) -- (1, 0) -- (1, 1);\n\\end{tikzpicture}$}}'
            sage: NuDyckWord('101100','101010')._latex_()
            '\\vcenter{\\hbox{$\\begin{tikzpicture}[scale=1]\n  \\draw[dotted] (0, 0) grid (3, 3);\n  \\draw[rounded corners=1, color=red, line width=1] (0, 0) -- (0, 1) -- (1, 1) -- (1, 2) -- (2, 2) -- (2, 3) -- (3, 3);\n  \\draw[rounded corners=1, color=black, line width=2] (0, 0) -- (0, 1) -- (1, 1) -- (1, 2) -- (1, 3) -- (2, 3) -- (3, 3);\n\\end{tikzpicture}$}}'
        """
        latex.add_package_to_preamble_if_available("tikz")
        latex_options = self.latex_options()

        # Start setting up tikz
        res = "\\vcenter{\\hbox{$\\begin{tikzpicture}"
        res += "[scale=" + str(latex_options['tikz_scale']) + "]"
        res += "\n"

        # Setup background grid
        if latex_options['show_grid']:
            grid = [((0, 0), (self.width(), self.height()))]
            for v1, v2 in grid:
                res += f"  \\draw[dotted] {v1} grid {v2};"
                res += "\n"

        # Add points if wanted
        if latex_options['show_points']:
            pt_color = latex_options['points_color']
            radius = 0.15 + .03 * latex_options['line width']
            for v in self.points():
                res += "  \\draw[line width=2,"
                res += f"color={pt_color},fill={pt_color}]"
                res += f"{v} circle ({radius});"
                res += "\n"

        # Add nu if wanted
        if latex_options['show_nu']:
            res += "  \\draw[%s]" % (str(latex_options['nu_options']))
            for k, p in enumerate(self._nu.points()):
                if k == 0:
                    res += " %s" % (str(p))
                else:
                    res += " -- %s" % (str(p))
            res += ";\n"

        # setup Path
        res += "  \\draw[rounded corners=1, color={}, line width={}]".format(
            latex_options['color'],
            str(latex_options['line width'])
        )
        for k, p in enumerate(self._path.points()):
            if k == 0:
                res += " %s" % (str(p))
            else:
                res += " -- %s" % (str(p))
        res += ";\n"
        res += "\\end{tikzpicture}$}}"
        return res

    def plot(self, **kwds):
        r"""
        Plot a `\nu`-Dyck word as a continuous path.

        EXAMPLES::

            sage: NDW = NuDyckWord('010','010')
            sage: NDW.plot()                                                            # needs sage.plot
            Graphics object consisting of 1 graphics primitive
        """
        from sage.plot.plot import list_plot
        return list_plot(list(self.points()), plotjoined=True, **kwds)

    def path(self):
        r"""
        Return the underlying path object.

        EXAMPLES::

            sage: NDW = NuDyckWord('10011001000','00100101001')
            sage: NDW.path()
            Path: 10011001000
        """
        return self._path

    def height(self):
        r"""
        Return the height of ``self``.

        The height is the number of ``north`` steps.

        EXAMPLES::

            sage: NuDyckWord('1101110011010010001101111000110000',
            ....: '1010101010101010101010101010101010').height()
            17
        """
        return self._path.height()

    def width(self):
        r"""
        Return the width of ``self``.

        The width is the number of ``east`` steps.

        EXAMPLES::

            sage: NuDyckWord('110111001101001000110111100011000',
            ....: '101010101010101010101010101010101').width()
            16
        """
        return self._path.width()

    def length(self):
        r"""
        Return the length of ``self``.

        The length is the total number of steps.

        EXAMPLES::

            sage: NDW = NuDyckWord('10011001000','00100101001')
            sage: NDW.length()
            11
        """
        return self._path.length()

    def points(self):
        r"""
        Return an iterator for the points on the `\nu`-Dyck path.

        EXAMPLES::

            sage: list(NuDyckWord('110111001101001000110111100011000',
            ....: '101010101010101010101010101010101')._path.points())
            [(0, 0),
             (0, 1),
             (0, 2),
             (1, 2),
             (1, 3),
             (1, 4),
             (1, 5),
             (2, 5),
             (3, 5),
             (3, 6),
             (3, 7),
             (4, 7),
             (4, 8),
             (5, 8),
             (6, 8),
             (6, 9),
             (7, 9),
             (8, 9),
             (9, 9),
             (9, 10),
             (9, 11),
             (10, 11),
             (10, 12),
             (10, 13),
             (10, 14),
             (10, 15),
             (11, 15),
             (12, 15),
             (13, 15),
             (13, 16),
             (13, 17),
             (14, 17),
             (15, 17),
             (16, 17)]
        """
        return self._path.points()

    def heights(self):
        r"""
        Return the heights of each point on ``self``.

        We view the Dyck word as a Dyck path from `(0,0)` to
        `(x,y)` in the first quadrant by letting ``1``'s represent
        steps in the direction `(0,1)` and ``0``'s represent steps in
        the direction `(1,0)`.

        The heights is the sequence of the `y`-coordinates of all
        `x+y` lattice points along the path.

        EXAMPLES::

            sage: NuDyckWord('010','010').heights()
            [0, 0, 1, 1]
            sage: NuDyckWord('110100','101010').heights()
            [0, 1, 2, 2, 3, 3, 3]
        """
        return self._path.height_vector()

    def widths(self):
        r"""
        Return the widths of each point on ``self``.

        We view the Dyck word as a Dyck path from `(0,0)` to
        `(x,y)` in the first quadrant by letting ``1``'s represent
        steps in the direction `(0,1)` and ``0``'s represent steps in
        the direction `(1,0)`.

        The widths is the sequence of the `x`-coordinates of all
        `x+y` lattice points along the path.

        EXAMPLES::

            sage: NuDyckWord('010','010').widths()
            [0, 1, 1, 2]
            sage: NuDyckWord('110100','101010').widths()
            [0, 0, 0, 1, 1, 2, 3]
        """
        return self._path.width_vector()

    def horizontal_distance(self):
        r"""
        Return a list of how far each point is from `\nu`.

        EXAMPLES::

            sage: NDW = NuDyckWord('10010100','00000111')
            sage: NDW.horizontal_distance()
            [5, 5, 4, 3, 3, 2, 2, 1, 0]
            sage: NDW = NuDyckWord('10010100','00001011')
            sage: NDW.horizontal_distance()
            [4, 5, 4, 3, 3, 2, 2, 1, 0]
            sage: NDW = NuDyckWord('10011001000','00100101001')
            sage: NDW.horizontal_distance()
            [2, 4, 3, 2, 3, 5, 4, 3, 3, 2, 1, 0]
        """
        # Grab furthest east point at each height of nu
        nu_points = list(self._nu.points())
        nu_easts = [max(i for i, j in nu_points if j == k)
                    for k in range(self._nu.height() + 1)]

        points = list(self._path.points())
        return [nu_easts[j] - i for i, j in points]

    def can_mutate(self, i) -> bool | int:
        """
        Return True/False based off if mutable at height `i`.

        Can only mutate if an east step is followed by a north step at height
        `i`.

        OUTPUT: whether we can mutate at height of `i`

        EXAMPLES::

            sage: NDW = NuDyckWord('10010100','00000111')
            sage: NDW.can_mutate(1)
            False
            sage: NDW.can_mutate(3)
            5

        TESTS::

            sage: NDW = NuDyckWord('10010100','00000111')
            sage: NDW.can_mutate(33)
            Traceback (most recent call last):
            ...
            ValueError: cannot mutate above or below path
        """
        if i > self.height() or i <= 0:
            raise ValueError('cannot mutate above or below path')

        # Find the ith north step
        level = 0
        ndw = self._list
        for j, k in enumerate(ndw):
            if k == ndw_open_symbol:
                level += 1
                if level == i:
                    break
        if j > 0 and ndw[j - 1] == ndw_close_symbol:
            return j
        return False

    def mutate(self, i) -> None | NuDyckWord:
        r"""
        Return a new `\nu`-Dyck Word if possible.

        If at height `i` we have an east step E meeting a north step N then we
        calculate all horizontal distances from this point until we find
        the first point that has the same horizontal distance to `\nu`. We let

        - d is everything up until EN (not including EN)

        - f be everything between N and the point with the same horizontal
          distance (including N)

        - g is everything after f

        .. SEEALSO:: :meth:`can_mutate`

        EXAMPLES::

            sage: NDW = NuDyckWord('10010100','00000111')
            sage: NDW.mutate(1)
            sage: NDW.mutate(3)
            [1, 0, 0, 1, 1, 0, 0, 0]
        """
        mutation_index = self.can_mutate(i)
        if not mutation_index:
            return None

        horiz = self.horizontal_distance()
        # Find horiz in between East and North Step
        horiz_num = horiz[mutation_index]
        other_index = len(horiz)
        for i in range(mutation_index + 1, len(horiz)):
            if horiz[i] == horiz_num:
                other_index = i
                break
        ndw = self._list
        d = ndw[0:mutation_index - 1]
        e = ndw[mutation_index:other_index]
        f = ndw[other_index:]
        return NuDyckWord(d + e + [ndw_close_symbol] + f, self._nu)


class NuDyckWords(Parent):
    r"""
    `\nu`-Dyck words.

    Given a lattice path `\nu` in the `\ZZ^2` grid starting at the origin
    `(0,0)` consisting of North `N = (0,1)` and East `E = (1,0)` steps, a
    `\nu`-Dyck path is a lattice path in the`\ZZ^2` grid starting at the
    origin `(0,0)` and ending at the same coordinate as `\nu` such that it is
    weakly above `\nu`. A `\nu`-Dyck word is the representation of a
    `\nu`-Dyck path where a North step is represented by a 1 and an East step
    is represented by a 0.

    INPUT:

    - ``nu`` -- the base lattice path

    EXAMPLES::

        sage: NDW = NuDyckWords('1010'); NDW
        [1, 0, 1, 0] Dyck words
        sage: [1,0,1,0] in NDW
        True
        sage: [1,1,0,0] in NDW
        True
        sage: [1,0,0,1] in NDW
        False
        sage: [0,1,0,1] in NDW
        False
        sage: NDW.cardinality()
        2
    """

    Element = NuDyckWord

    def __init__(self, nu=()):
        """
        Intialize ``self``.

        EXAMPLES::

            sage: TestSuite(NuDyckWords(nu=[1,0,1])).run()
        """
        Parent.__init__(self, category=FiniteEnumeratedSets())

        self._nu = to_word_path(nu)
        if self._nu is None:
            raise ValueError("invalid nu supplied")

    # add options to class
    class options(GlobalOptions):
        r"""
        Set and display the options for `\nu`-Dyck words. If no parameters
        are set, then the function returns a copy of the options dictionary.

        The ``options`` to `\nu`-Dyck words can be accessed as the method
        :meth:`NuDyckWords.options` of :class:`NuDyckWords` and
        related parent classes.

        @OPTIONS

        EXAMPLES::

            sage: ND = NuDyckWords('101')
            sage: ND
            [1, 0, 1] Dyck words
            sage: ND.options
            Current options for NuDyckWords
              - ascii_art:               pretty_output
              - diagram_style:           grid
              - display:                 list
              - latex_color:             black
              - latex_line_width_scalar: 2
              - latex_nu_options:        rounded corners=1, color=red, line width=1
              - latex_points_color:      black
              - latex_show_grid:         True
              - latex_show_nu:           True
              - latex_show_points:       False
              - latex_tikz_scale:        1
        """
        NAME = 'NuDyckWords'
        module = 'sage.combinat.nu_dyck_path'
        display = {'default': "list",
                   'description': 'Specifies how nu Dyck words should be printed',
                   'values': {'list': 'displayed as a list',
                              'lattice': 'displayed on the lattice defined by ``diagram_style``'},
                   'case_sensitive': False}
        ascii_art = {'default': "pretty_output",
                     'description': 'Specifies how the ascii art of nu Dyck words should be printed',
                     'values': {'pretty_output': "Using pretty printing"},
                     'alias': {'pretty_print': "pretty_output"},
                     'case_sensitive': False}
        diagram_style = {'default': "grid",
                         'values': {
                             'grid': 'printing as paths on a grid using N and E steps'},
                         'alias': {'N-E': 'grid'},
                         'case_sensitive': False}
        latex_tikz_scale = {'default': 1,
                            'description': 'The default value for the tikz scale when latexed',
                            'checker': lambda x: True}  # More trouble than it's worth to check
        latex_line_width_scalar = {'default': 2,
                                   'description': 'The default value for the line width as a '
                                   'multiple of the tikz scale when latexed',
                                   'checker': lambda x: True}  # More trouble than it's worth to check
        latex_color = {'default': "black",
                       'description': 'The default value for the color when latexed',
                       'checker': lambda x: isinstance(x, str)}
        latex_show_points = {'default': False,
                             'description': 'The default value for showing points',
                             'checker': lambda x: isinstance(x, bool)}
        latex_points_color = {'default': 'black',
                              'description': 'The default value for path color.',
                              'checker': lambda x: isinstance(x, str)}
        latex_show_grid = {'default': True,
                           'description': 'The default value for showing grid',
                           'checker': lambda x: isinstance(x, bool)}
        latex_show_nu = {'default': True,
                         'description': 'The default value for showing nu',
                         'checker': lambda x: isinstance(x, bool)}
        latex_nu_options = {'default': 'rounded corners=1, color=red, line width=1',
                            'description': 'The default value for options for nu path',
                            'checker': lambda x: isinstance(x, str)}

    def _element_constructor_(self, word):
        """
        Construct an element of ``self``.

        EXAMPLES::

            sage: NDW = NuDyckWords('101')
            sage: elt = NDW('110'); elt
            [1, 1, 0]
            sage: elt.parent() is NDW
            True
        """
        if isinstance(word, NuDyckWord) and word.parent() is self:
            return word
        return self.element_class(self, to_word_path(word))

    def __contains__(self, x) -> bool:
        r"""
        Check for containment.

        TESTS::

            sage: NDW = NuDyckWords([1,0,1,1])
            sage: [1,1,0,1] in NDW
            True
            sage: [1,0,1,1] in NDW
            True
            sage: [0] in NDW
            False
            sage: [1, 0] in NDW
            False
        """
        return path_weakly_above_other(to_word_path(x), self._nu)

    def __eq__(self, other):
        """
        Return equality.

        TESTS::

            sage: A = NuDyckWords([1,0,1,1])
            sage: B = NuDyckWords([1,0,1,1])
            sage: C = NuDyckWords([1,0,1,1,1])
            sage: A == B
            True
            sage: A == C
            False
        """
        if not isinstance(other, NuDyckWords):
            return False
        return self._nu == other._nu

    def __neq__(self, other):
        """
        Return inequality.

        TESTS::

            sage: A = NuDyckWords([1,0,1,1])
            sage: B = NuDyckWords([1,0,1,1])
            sage: C = NuDyckWords([1,0,1,1,1])
            sage: A != B
            False
            sage: A != C
            True
        """
        return not self.__eq__(other)

    def _repr_(self) -> str:
        r"""
        TESTS::

            sage: NuDyckWords([1,0,1,1])
            [1, 0, 1, 1] Dyck words
        """
        return f"{list(self._nu)} Dyck words"

    def _cache_key(self) -> str:
        """
        Return a cache key

        TESTS::

            sage: NuDyckWords([1,0,1,1])._cache_key()
            '1011'
        """
        return str(self._nu)

    def _an_element_(self):
        r"""
        Return an element.

        TESTS::

            sage: NuDyckWords('101').an_element()
            [1, 0, 1]
        """
        return self.element_class(self, self._nu)

    def __iter__(self, N=[], D=[], i=None, X=None):
        """
        Iterate over ``self``.

        The iterator interchanges a 0,1 pair whenever the 0 comes before a 1

        EXAMPLES::

            sage: it = NuDyckWords('101010').__iter__()
            sage: [i for i in it]
            [[1, 0, 1, 0, 1, 0],
             [1, 1, 0, 0, 1, 0],
             [1, 0, 1, 1, 0, 0],
             [1, 1, 0, 1, 0, 0],
             [1, 1, 1, 0, 0, 0]]
        """
        # Define successor function for recursion
        def transpose_close_open(N):
            for k, v in enumerate(N._list):
                if k > 0 and v == ndw_open_symbol:
                    w = N._list[k - 1]
                    if w == ndw_close_symbol:
                        new = N._list[:k - 1] + [v, w] + N._list[k + 1:]
                        yield self.element_class(self, new)

        RES = RecursivelyEnumeratedSet([self.element_class(self, self._nu)],
                                       transpose_close_open)
        return RES.breadth_first_search_iterator()

    def cardinality(self):
        r"""
        Return the number of `\nu`-Dyck words.

        EXAMPLES::

            sage: NDW = NuDyckWords('101010'); NDW.cardinality()
            5
            sage: NDW = NuDyckWords('1010010'); NDW.cardinality()
            7
            sage: NDW = NuDyckWords('100100100'); NDW.cardinality()
            12
        """
        return Integer(len([1 for _ in self.__iter__()]))


def to_word_path(word):
    r"""
    Convert input into a word path over an appropriate alphabet.

    Helper function.

    INPUT:

    - ``word`` -- word to convert to wordpath

    OUTPUT: a ``FiniteWordPath_north_east`` object

    EXAMPLES::

        sage: from sage.combinat.nu_dyck_word import to_word_path
        sage: wp = to_word_path('NEENENEN'); wp
        Path: 10010101
        sage: from sage.combinat.words.paths import FiniteWordPath_north_east
        sage: isinstance(wp,FiniteWordPath_north_east)
        True
        sage: to_word_path('1001')
        Path: 1001
        sage: to_word_path([0,1,0,0,1,0])
        Path: 010010
    """
    # If we already have the object, don't worry
    if isinstance(word, FiniteWordPath_north_east):
        return word

    # If we have a Nu Dyck Path, return the path it contains
    if isinstance(word, NuDyckWord):
        return word.path()

    # if we have a string, convert to list
    if isinstance(word, str):
        word = map(replace_dyck_char, word)

    # By default "north" is first symbol, "east" is second symbol
    P = WordPaths_north_east([ndw_open_symbol, ndw_close_symbol])

    return P(word)


def path_weakly_above_other(path, other) -> bool:
    r"""
    Test if ``path`` is weakly above ``other``.

    A path `P` is wealy above another path `Q` if `P` and `Q` are the same
    length and if any prefix of length `n` of `Q` contains more North steps
    than the prefix of length `n` of `P`.

    INPUT:

    - ``path`` -- the path to verify is weakly above the other path

    - ``other`` -- the other path to verify is weakly below the path

    OUTPUT: boolean

    EXAMPLES::

        sage: from sage.combinat.nu_dyck_word import path_weakly_above_other
        sage: path_weakly_above_other('1001','0110')
        False
        sage: path_weakly_above_other('1001','0101')
        True
        sage: path_weakly_above_other('1111','0101')
        False
        sage: path_weakly_above_other('111100','0101')
        False
    """
    # Ensure we have word paths:
    path = to_word_path(path)
    other = to_word_path(other)

    # Must be same length and must have same height
    if path.length() != other.length() or path.height() != other.height():
        return False

    # path is above other if height is always >= height of other
    p_height = path.height_vector()
    o_height = other.height_vector()
    return all(p_h >= o_h for p_h, o_h in zip(p_height, o_height))
