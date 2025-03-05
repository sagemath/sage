r"""
Dyck Words

A class of an object enumerated by the
:func:`Catalan numbers<sage.combinat.combinat.catalan_number>`,
see [Sta-EC2]_, [StaCat98]_ for details.

AUTHORS:

- Mike Hansen

- Dan Drake (2008-05-30): DyckWordBacktracker support

- Florent Hivert (2009-02-01): Bijections with NonDecreasingParkingFunctions

- Christian Stump (2011-12): added combinatorial maps and statistics

- Mike Zabrocki:

  * (2012-10): added pretty print, characteristic function, more functions
  * (2013-01): added inverse of area/dinv, bounce/area map

- Jean--Baptiste Priez, Travis Scrimshaw (2013-05-17): Added ASCII art

- Travis Scrimshaw (2013-07-09): Removed ``CombinatorialClass`` and added
  global options.

REFERENCES:

.. [Sta-EC2] Richard P. Stanley.
   *Enumerative Combinatorics*, Volume 2.
   Cambridge University Press, 2001.

.. [StaCat98] Richard Stanley. *Exercises on Catalan and Related Numbers
   excerpted from Enumerative Combinatorics, vol. 2 (CUP 1999)*,
   version of 23 June 1998.
   http://www-math.mit.edu/~rstan/ec/catalan.pdf

.. [Hag2008] James Haglund. *The* `q,t` -- *Catalan Numbers and the
   Space of Diagonal Harmonics:
   With an Appendix on the Combinatorics of Macdonald Polynomials*.
   University of Pennsylvania, Philadelphia -- AMS, 2008, 167 pp.

.. [BK2001] \J. Bandlow, K. Killpatrick -- *An area-to_inv bijection
   between Dyck paths and 312-avoiding permutations*, Electronic
   Journal of Combinatorics, Volume 8, Issue 1 (2001).

.. [EP2004] \S. Elizalde, I. Pak. *Bijections for refined restricted
   permutations**. JCTA 105(2) 2004.

.. [CK2008] \A. Claesson, S. Kitaev. *Classification of bijections
   between `321`- and `132`- avoiding permutations*. Séminaire
   Lotharingien de Combinatoire **60** 2008. :arxiv:`0805.1325`.

.. [Knu1973] \D. Knuth. *The Art of Computer Programming, Vol. III*.
   Addison-Wesley. Reading, MA. 1973.

.. [Kra2001] \C. Krattenthaler -- *Permutations with restricted
   patterns and Dyck paths*, Adv. Appl. Math. 27 (2001), 510--530.

.. [DS1992] \A. Denise, R. Simion, *Two combinatorial statistics on
   Dyck paths*, Discrete Math 137 (1992), 155--176.
"""

# ****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
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
from collections.abc import Iterator

from .combinat import CombinatorialElement, catalan_number
from sage.combinat.combinatorial_map import combinatorial_map
from .backtrack import GenericBacktracker

from sage.structure.global_options import GlobalOptions
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.posets import Posets

from sage.rings.integer import Integer
from sage.rings.rational_field import QQ
from sage.combinat.permutation import Permutation, Permutations
from sage.combinat.words.word import Word
from sage.combinat.set_partition import SetPartitions
from sage.misc.latex import latex
from sage.misc.lazy_import import lazy_import

lazy_import('sage.combinat.alternating_sign_matrix', 'AlternatingSignMatrices')


open_symbol = 1
close_symbol = 0


def replace_parens(x):
    r"""
    A map sending ``'('`` to ``open_symbol`` and ``')'`` to
    ``close_symbol``, and raising an error on any input other than
    ``'('``  and ``')'``. The values of the constants ``open_symbol``
    and ``close_symbol`` are subject to change.

    This is the inverse map of :func:`replace_symbols`.

    INPUT:

    - ``x`` -- either an opening or closing parenthesis

    OUTPUT:

    - If ``x`` is an opening parenthesis, replace ``x`` with the
      constant ``open_symbol``.

    - If ``x`` is a closing parenthesis, replace ``x`` with the
      constant ``close_symbol``.

    - Raise a :exc:`ValueError` if ``x`` is neither an opening nor a
      closing parenthesis.

    .. SEEALSO:: :func:`replace_symbols`

    EXAMPLES::

        sage: from sage.combinat.dyck_word import replace_parens
        sage: replace_parens('(')
        1
        sage: replace_parens(')')
        0
        sage: replace_parens(1)
        Traceback (most recent call last):
        ...
        ValueError
    """
    if x == '(':
        return open_symbol
    if x == ')':
        return close_symbol
    raise ValueError


def replace_symbols(x):
    r"""
    A map sending ``open_symbol`` to ``'('`` and ``close_symbol`` to ``')'``,
    and raising an error on any input other than ``open_symbol`` and
    ``close_symbol``. The values of the constants ``open_symbol``
    and ``close_symbol`` are subject to change.

    This is the inverse map of :func:`replace_parens`.

    INPUT:

    - ``x`` -- either ``open_symbol`` or ``close_symbol``

    OUTPUT:

    - If ``x`` is ``open_symbol``, replace ``x`` with ``'('``.

    - If ``x`` is ``close_symbol``, replace ``x`` with ``')'``.

    - If ``x`` is neither ``open_symbol`` nor ``close_symbol``, a
      :exc:`ValueError` is raised.

    .. SEEALSO:: :func:`replace_parens`

    EXAMPLES::

        sage: from sage.combinat.dyck_word import replace_symbols
        sage: replace_symbols(1)
        '('
        sage: replace_symbols(0)
        ')'
        sage: replace_symbols(3)
        Traceback (most recent call last):
        ...
        ValueError
    """
    if x == open_symbol:
        return '('
    if x == close_symbol:
        return ')'
    raise ValueError


class DyckWord(CombinatorialElement):
    r"""
    A Dyck word.

    A Dyck word is a sequence of open and close symbols such that every close
    symbol has a corresponding open symbol preceding it. That is to say, a
    Dyck word of length `n` is a list with `k` entries 1 and `n - k`
    entries 0 such that the first `i` entries always have at least as many 1s
    among them as 0s. (Here, the 1 serves as the open symbol and the 0 as the
    close symbol.)  Alternatively, the alphabet 1 and 0 can be replaced by
    other characters such as '(' and ')'.

    A Dyck word is *complete* if every open symbol moreover has a corresponding
    close symbol.

    A Dyck word may also be specified by either a noncrossing partition or
    by an area sequence or the sequence of heights.

    A Dyck word may also be thought of as a lattice path in the `\ZZ^2`
    grid, starting at the origin `(0,0)`, and with steps in the North
    `N = (0,1)` and east `E = (1,0)` directions such that it does not pass
    below the `x = y` diagonal. The diagonal is referred to as the "main
    diagonal" in the documentation. A North step is represented by a 1 in
    the list and an East step is represented by a 0.

    Equivalently, the path may be represented with steps in
    the `NE = (1,1)` and the `SE = (1,-1)` direction such that it does not
    pass below the horizontal axis.

    .. PLOT::
        :width: 400 px

        d = DyckWord([1,0,1,1,1,1,0,1,0,0,1,0,1,1,0,1,0,1,1,0,0,0,0,0])
        sphinx_plot(d.plot(aspect_ratio=1))

    A path representing a Dyck word (either using `N` and `E` steps, or
    using `NE` and `SE` steps) is called a Dyck path.

    EXAMPLES::

        sage: dw = DyckWord([1, 0, 1, 0]); dw
        [1, 0, 1, 0]
        sage: print(dw)
        ()()
        sage: dw.height()
        1
        sage: dw.to_noncrossing_partition()
        {{1}, {2}}

    ::

        sage: DyckWord('()()')
        [1, 0, 1, 0]
        sage: DyckWord('(())')
        [1, 1, 0, 0]
        sage: DyckWord('((')
        [1, 1]
        sage: DyckWord('')
        []

    ::

        sage: DyckWord(noncrossing_partition=[[1],[2]])
        [1, 0, 1, 0]
        sage: DyckWord(noncrossing_partition=[[1,2]])
        [1, 1, 0, 0]
        sage: DyckWord(noncrossing_partition=[])
        []

    ::

        sage: DyckWord(area_sequence=[0,0])
        [1, 0, 1, 0]
        sage: DyckWord(area_sequence=[0,1])
        [1, 1, 0, 0]
        sage: DyckWord(area_sequence=[0,1,2,2,0,1,1,2])
        [1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0]
        sage: DyckWord(area_sequence=[])
        []

    ::

        sage: DyckWord(heights_sequence=(0,1,0,1,0))
        [1, 0, 1, 0]
        sage: DyckWord(heights_sequence=(0,1,2,1,0))
        [1, 1, 0, 0]
        sage: DyckWord(heights_sequence=(0,))
        []

    ::

        sage: print(DyckWord([1,0,1,1,0,0]).to_path_string())
           /\
        /\/  \
        sage: DyckWord([1,0,1,1,0,0]).pretty_print()
           ___
          | x
         _|  .
        |  . .
    """
    @staticmethod
    def __classcall_private__(cls, dw=None, noncrossing_partition=None,
                              area_sequence=None, heights_sequence=None,
                              catalan_code=None):
        """
        Return an element with the appropriate parent.

        EXAMPLES::

            sage: DyckWord([1,0,1,1,0,0])
            [1, 0, 1, 1, 0, 0]
            sage: DyckWord(heights_sequence=(0,1,2,1,0))
            [1, 1, 0, 0]
            sage: DyckWord(noncrossing_partition=[[1],[2]])
            [1, 0, 1, 0]
        """
        if dw is None:
            if catalan_code is not None:
                return CompleteDyckWords_all().from_Catalan_code(catalan_code)
            if area_sequence is not None:
                return CompleteDyckWords_all().from_area_sequence(area_sequence)
            if noncrossing_partition is not None:
                return CompleteDyckWords_all().from_noncrossing_partition(noncrossing_partition)
            if heights_sequence is not None:
                if heights_sequence[-1] == 0:
                    P = CompleteDyckWords_all()
                else:
                    P = DyckWords_all()
                return P.from_heights(heights_sequence)

            raise ValueError("you have not specified a Dyck word")

        if isinstance(dw, str):
            l = [replace_parens(x) for x in dw]
        else:
            l = dw

        if isinstance(l, DyckWord):
            return l

        if l in CompleteDyckWords_all():
            return CompleteDyckWords_all()(l)

        if is_a(l):
            return DyckWords_all()(l)

        raise ValueError("invalid Dyck word")

    def __init__(self, parent, l, latex_options={}):
        r"""
        TESTS::

            sage: DW = DyckWords(complete=False).from_heights((0,))
            sage: TestSuite(DW).run()
            sage: DW = DyckWords(complete=False).min_from_heights((0,))
            sage: TestSuite(DW).run()
            sage: DW = DyckWords().from_Catalan_code([])
            sage: TestSuite(DW).run()
            sage: DW = DyckWords().from_area_sequence([])
            sage: TestSuite(DW).run()
        """
        CombinatorialElement.__init__(self, parent, l)
        self._latex_options = dict(latex_options)

    def set_latex_options(self, D):
        r"""
        Set the latex options for use in the ``_latex_`` function.

        The default values are set in the ``__init__`` function.

        - ``tikz_scale`` -- (default: 1) scale for use with the tikz package

        - ``diagonal`` -- boolean (default: ``False``); value to draw the
          diagonal or not

        - ``line width`` -- (default: 2*``tikz_scale``) value representing the
          line width

        - ``color`` -- (default: black) the line color

        - ``bounce path`` -- boolean (default: ``False``); value to indicate
          if the bounce path should be drawn

        - ``peaks`` -- boolean (default: ``False``); value to indicate if the
          peaks should be displayed

        - ``valleys`` -- boolean (default: ``False``); value to indicate if the
          valleys should be displayed

        INPUT:

        - ``D`` -- dictionary with a list of latex parameters to change

        EXAMPLES::

            sage: D = DyckWord([1,0,1,0,1,0])
            sage: D.set_latex_options({"tikz_scale":2})
            sage: D.set_latex_options({"valleys":True, "color":"blue"})

        .. TODO::

            This should probably be merged into DyckWord.options.
        """
        for opt in D:
            self._latex_options[opt] = D[opt]

    def latex_options(self):
        r"""
        Return the latex options for use in the ``_latex_`` function as a
        dictionary.

        The default values are set using the options.

        - ``tikz_scale`` -- (default: 1) scale for use with the tikz package

        - ``diagonal`` -- boolean (default: ``False``); value to draw the
          diagonal or not

        - ``line width`` -- (default: ``2*tikz_scale``) value representing the
          line width

        - ``color`` -- (default: black) the line color

        - ``bounce path`` -- boolean (default: ``False``); value to indicate
          if the bounce path should be drawn

        - ``peaks`` -- boolean (default: ``False``); value to indicate if the
          peaks should be displayed

        - ``valleys`` -- boolean (default: ``False``); value to indicate if the
          valleys should be displayed

        EXAMPLES::

            sage: D = DyckWord([1,0,1,0,1,0])
            sage: D.latex_options()
            {'bounce path': False,
             'color': black,
             'diagonal': False,
             'line width': 2,
             'peaks': False,
             'tikz_scale': 1,
             'valleys': False}

        .. TODO::

            This should probably be merged into DyckWord.options.
        """
        d = self._latex_options.copy()
        if "tikz_scale" not in d:
            d["tikz_scale"] = self.parent().options.latex_tikz_scale
        if "diagonal" not in d:
            d["diagonal"] = self.parent().options.latex_diagonal
        if "line width" not in d:
            d["line width"] = self.parent().options.latex_line_width_scalar * d["tikz_scale"]
        if "color" not in d:
            d["color"] = self.parent().options.latex_color
        if "bounce path" not in d:
            d["bounce path"] = self.parent().options.latex_bounce_path
        if "peaks" not in d:
            d["peaks"] = self.parent().options.latex_peaks
        if "valleys" not in d:
            d["valleys"] = self.parent().options.latex_valleys
        return d

    def _repr_(self) -> str:
        r"""
        Return a string representation of ``self`` depending on
        :meth:`DyckWords.options`.

        TESTS::

            sage: DyckWord([1, 0, 1, 0])
            [1, 0, 1, 0]
            sage: DyckWord([1, 1, 0, 0])
            [1, 1, 0, 0]
            sage: DyckWords.options.display="lattice"
            sage: DyckWords.options.diagram_style="line"
            sage: DyckWord([1, 0, 1, 0])
            /\/\
            sage: DyckWord([1, 1, 0, 0])
             /\
            /  \
            sage: DyckWords.options._reset()
        """
        return self.parent().options._dispatch(self, '_repr_', 'display')

    def _repr_list(self) -> str:
        r"""
        Return a string representation of ``self`` as a list.

        TESTS::

            sage: DyckWord([])
            []
            sage: DyckWord([1, 0])
            [1, 0]
            sage: DyckWord('(())')
            [1, 1, 0, 0]
        """
        return super()._repr_()

    def _repr_lattice(self, type=None, labelling=None, underpath=True) -> str:
        r"""
        See :meth:`pretty_print()`.

        TESTS::

            sage: print(DyckWord(area_sequence=[0,1,0])._repr_lattice(type="NE-SE"))
             /\
            /  \/\
            sage: print(DyckWord(area_sequence=[0,1,0])._repr_lattice(labelling=[1,3,2],underpath=False))
                 _
             ___|  2
            | x  . 3
            |  . . 1
        """
        if type is None:
            type = self.parent().options.diagram_style
            if type == "grid":
                type = "N-E"
            elif type == "line":
                type = "NE-SE"

        if type == "NE-SE":
            if labelling is not None or underpath is not True:
                raise ValueError("the labelling cannot be shown with Northeast-Southeast paths")
            return self.to_path_string()
        elif type == "N-E":
            alst = self.to_area_sequence()
            n = len(alst)
            if n == 0:
                return ".\n"
            if labelling is None:
                labels = [" "] * n
            else:
                if len(labelling) != n:
                    raise ValueError("the given labelling has the wrong length")
                labels = [str(label) for label in labelling]
                if not underpath:
                    max_length = max(len(label) for label in labels)
                    labels = [lbl.rjust(max_length + 1) for lbl in labels]

            length_of_final_fall = list(reversed(self)).index(open_symbol)
            if length_of_final_fall == 0:
                final_fall = " "
            else:
                final_fall = " _" + "__" * (length_of_final_fall - 1)
            row = "  " * (n - alst[-1] - 1) + final_fall + "\n"
            for i in range(n - 1):
                c = 0
                row = row + "  " * (n-i-2-alst[-i-2])
                c += n-i-2-alst[-i-2]
                if alst[-i-2]+1 != alst[-i-1]:
                    row += " _"
                c += alst[-i-2] - alst[-i-1]
                if underpath:
                    row += "__" * (alst[-i-2]-alst[-i-1]) + "|" + labels[-1] + "x "*(n-c-2-i) + " ." * i + "\n"
                else:
                    row += "__"*(alst[-i-2]-alst[-i-1])+"| " + "x "*(n-c-2-i) + " ."*i + labels[-1] + "\n"
                labels.pop()
            if underpath:
                row += "|" + labels[-1] + " ." * (n - 1) + "\n"
            else:
                row += "| " + " ." * (n - 1) + labels[-1] + "\n"
            return row
        else:
            raise ValueError("the given type (=%s) is not valid" % type)

    def _ascii_art_(self):
        r"""
        Return an ASCII art representation of ``self``.

        TESTS::

            sage: ascii_art(list(DyckWords(3)))
            [                                   /\   ]
            [            /\    /\      /\/\    /  \  ]
            [ /\/\/\, /\/  \, /  \/\, /    \, /    \ ]
        """
        from sage.typeset.ascii_art import AsciiArt
        rep = self.parent().options.ascii_art
        if rep == "path":
            ret = self.to_path_string()
        elif rep == "pretty_output":
            ret = self._repr_lattice()
        return AsciiArt(ret.splitlines(), baseline=0)

    def _unicode_art_(self):
        r"""
        Return an unicode art representation of this Dyck word.

        EXAMPLES::

            sage: unicode_art(list(DyckWords(3)))
            ⎡                                   ╱╲   ⎤
            ⎢            ╱╲    ╱╲      ╱╲╱╲    ╱  ╲  ⎥
            ⎣ ╱╲╱╲╱╲, ╱╲╱  ╲, ╱  ╲╱╲, ╱    ╲, ╱    ╲ ⎦
        """
        from sage.typeset.unicode_art import UnicodeArt
        return UnicodeArt(self.to_path_string(unicode=True).splitlines())

    def __str__(self) -> str:
        r"""
        Return a string consisting of matched parentheses corresponding to
        the Dyck word.

        EXAMPLES::

            sage: print(DyckWord([1, 0, 1, 0]))
            ()()
            sage: print(DyckWord([1, 1, 0, 0]))
            (())
        """
        return "".join(replace_symbols(x) for x in self)

    def to_path_string(self, unicode=False) -> str:
        r"""
        Return a path representation of the Dyck word consisting of steps
        ``/`` and ``\`` .

        INPUT:

        - ``unicode`` -- boolean (default: ``False``); whether to use unicode

        EXAMPLES::

            sage: print(DyckWord([1, 0, 1, 0]).to_path_string())
            /\/\
            sage: print(DyckWord([1, 1, 0, 0]).to_path_string())
             /\
            /  \
            sage: print(DyckWord([1,1,0,1,1,0,0,1,0,1,0,0]).to_path_string())
                /\
             /\/  \/\/\
            /          \
        """
        if unicode:
            import unicodedata
            space = ' '
            up = unicodedata.lookup('BOX DRAWINGS LIGHT DIAGONAL UPPER RIGHT TO LOWER LEFT')
            down = unicodedata.lookup('BOX DRAWINGS LIGHT DIAGONAL UPPER LEFT TO LOWER RIGHT')
        else:
            space = ' '
            up = '/'
            down = '\\'

        res = [([space] * len(self)) for _ in range(self.height())]
        h = 1
        for i, p in enumerate(self):
            if p == open_symbol:
                res[-h][i] = up
                h += 1
            else:
                h -= 1
                res[-h][i] = down
        return "\n".join("".join(l) for l in res)

    def pretty_print(self, type=None, labelling=None, underpath=True):
        r"""
        Display a DyckWord as a lattice path in the `\ZZ^2` grid.

        If the ``type`` is "N-E", then a cell below the diagonal is
        indicated by a period, whereas a cell below the path but above
        the diagonal is indicated by an x. If a list of labels is
        included, they are displayed along the vertical edges of the
        Dyck path.

        If the ``type`` is "NE-SE", then the path is simply printed
        as up steps and down steps.

        INPUT:

        - ``type`` -- (default: ``None``) can either be:

          - ``None`` to use the option default
          - "N-E" to show ``self`` as a path of north and east steps, or
          - "NE-SE" to show ``self`` as a path of north-east and
            south-east steps.

        - ``labelling`` -- (if type is "N-E") a list of labels assigned to
          the up steps in ``self``

        - ``underpath`` -- (if type is "N-E", default:``True``) if ``True``,
          the labelling is shown under the path; otherwise, it is shown to
          the right of the path

        EXAMPLES::

            sage: for D in DyckWords(3): D.pretty_print()
                 _
               _|
             _|  .
            |  . .
               ___
              | x
             _|  .
            |  . .
                 _
             ___|
            | x  .
            |  . .
               ___
             _| x
            | x  .
            |  . .
             _____
            | x x
            | x  .
            |  . .

        ::

            sage: for D in DyckWords(3): D.pretty_print(type="NE-SE")
            /\/\/\
               /\
            /\/  \
             /\
            /  \/\
             /\/\
            /    \
              /\
             /  \
            /    \

        ::

            sage: D = DyckWord([1,1,1,0,1,0,0,1,1])
            sage: D.pretty_print()
                  | x x
               ___| x  .
             _| x x  . .
            | x x  . . .
            | x  . . . .
            |  . . . . .

            sage: D = DyckWord([1,1,1,0,1,0,0,1,1,0])
            sage: D.pretty_print()
                   _
                  | x x
               ___| x  .
             _| x x  . .
            | x x  . . .
            | x  . . . .
            |  . . . . .

            sage: D = DyckWord([1,1,1,0,1,0,0,1,1,0,0])
            sage: D.pretty_print()
                   ___
                  | x x
               ___| x  .
             _| x x  . .
            | x x  . . .
            | x  . . . .
            |  . . . . .

        ::

            sage: DyckWord(area_sequence=[0,1,0]).pretty_print(labelling=[1,3,2])
                 _
             ___|2
            |3x  .
            |1 . .

            sage: DyckWord(area_sequence=[0,1,0]).pretty_print(labelling=[1,3,2],underpath=False)
                 _
             ___|  2
            | x  . 3
            |  . . 1

        ::

            sage: DyckWord(area_sequence=[0,1,1,2,3,2,3,3,2,0,1,1,2,3,4,2,3]).pretty_print()
                                       _______
                                      | x x x
                                 _____| x x  .
                                | x x x x  . .
                                | x x x  . . .
                                | x x  . . . .
                               _| x  . . . . .
                              | x  . . . . . .
                         _____|  . . . . . . .
                     ___| x x  . . . . . . . .
                   _| x x x  . . . . . . . . .
                  | x x x  . . . . . . . . . .
               ___| x x  . . . . . . . . . . .
              | x x x  . . . . . . . . . . . .
              | x x  . . . . . . . . . . . . .
             _| x  . . . . . . . . . . . . . .
            | x  . . . . . . . . . . . . . . .
            |  . . . . . . . . . . . . . . . .

            sage: DyckWord(area_sequence=[0,1,1,2,3,2,3,3,2,0,1,1,2,3,4,2,3]).pretty_print(labelling=list(range(17)),underpath=False)
                                       _______
                                      | x x x  16
                                 _____| x x  . 15
                                | x x x x  . . 14
                                | x x x  . . . 13
                                | x x  . . . . 12
                               _| x  . . . . . 11
                              | x  . . . . . . 10
                         _____|  . . . . . . .  9
                     ___| x x  . . . . . . . .  8
                   _| x x x  . . . . . . . . .  7
                  | x x x  . . . . . . . . . .  6
               ___| x x  . . . . . . . . . . .  5
              | x x x  . . . . . . . . . . . .  4
              | x x  . . . . . . . . . . . . .  3
             _| x  . . . . . . . . . . . . . .  2
            | x  . . . . . . . . . . . . . . .  1
            |  . . . . . . . . . . . . . . . .  0

        ::

            sage: DyckWord([]).pretty_print()
            .
        """
        print(self._repr_lattice(type, labelling, underpath))

    pp = pretty_print

    def _latex_(self) -> str:
        r"""
        A latex representation of ``self`` using the tikzpicture package.

        EXAMPLES::

            sage: D = DyckWord([1,0,1,1,1,0,1,1,0,0,0,1,0,0])
            sage: D.set_latex_options({"valleys":True, "peaks":True, "bounce path":True})
            sage: latex(D)
            \vcenter{\hbox{$\begin{tikzpicture}[scale=1]
              \draw[line width=2,color=red,fill=red] (2, 0) circle (0.21);
              \draw[line width=2,color=red,fill=red] (6, 2) circle (0.21);
              \draw[line width=2,color=red,fill=red] (11, 1) circle (0.21);
              \draw[line width=2,color=red,fill=red] (1, 1) circle (0.21);
              \draw[line width=2,color=red,fill=red] (5, 3) circle (0.21);
              \draw[line width=2,color=red,fill=red] (8, 4) circle (0.21);
              \draw[line width=2,color=red,fill=red] (12, 2) circle (0.21);
              \draw[rounded corners=1, color=green, line width=4] (0, 0) -- (1, 1) -- (2, 0) -- (3, 1) -- (4, 0) -- (5, 1) -- (6, 2) -- (7, 3) -- (8, 2) -- (9, 1) -- (10, 0) -- (11, 1) -- (12, 2) -- (13, 1) -- (14, 0);
              \draw[dotted] (0, 0) grid (14, 4);
              \draw[rounded corners=1, color=black, line width=2] (0, 0) -- (1, 1) -- (2, 0) -- (3, 1) -- (4, 2) -- (5, 3) -- (6, 2) -- (7, 3) -- (8, 4) -- (9, 3) -- (10, 2) -- (11, 1) -- (12, 2) -- (13, 1) -- (14, 0);
            \end{tikzpicture}$}}
            sage: DyckWord([1,0])._latex_()
            '\\vcenter{\\hbox{$\\begin{tikzpicture}[scale=1]\n  \\draw[dotted] (0, 0) grid (2, 1);\n  \\draw[rounded corners=1, color=black, line width=2] (0, 0) -- (1, 1) -- (2, 0);\n\\end{tikzpicture}$}}'
            sage: DyckWord([1,0,1,1,0,0])._latex_()
            '\\vcenter{\\hbox{$\\begin{tikzpicture}[scale=1]\n  \\draw[dotted] (0, 0) grid (6, 2);\n  \\draw[rounded corners=1, color=black, line width=2] (0, 0) -- (1, 1) -- (2, 0) -- (3, 1) -- (4, 2) -- (5, 1) -- (6, 0);\n\\end{tikzpicture}$}}'
        """
        latex.add_package_to_preamble_if_available("tikz")
        heights = self.heights()
        latex_options = self.latex_options()
        diagonal = latex_options["diagonal"]
        ht = [(0, 0)]
        valleys = []
        peaks = []
        for i in range(1, len(heights)):
            a, b = ht[-1]
            if heights[i] > heights[i - 1]:
                if diagonal:
                    ht.append((a, b + 1))
                else:
                    ht.append((a + 1, b + 1))
                if i < len(heights) - 1 and heights[i + 1] < heights[i]:
                    peaks.append(ht[-1])
            else:
                if diagonal:
                    ht.append((a + 1, b))
                else:
                    ht.append((a + 1, b - 1))
                if i < len(heights) - 1 and heights[i + 1] > heights[i]:
                    valleys.append(ht[-1])
        hti = iter(ht)
        if diagonal:
            grid = [((0, i), (i, i + 1))
                    for i in range(self.number_of_open_symbols())]
        else:
            grid = [((0, 0), (len(self), self.height()))]
        res = "\\vcenter{\\hbox{$\\begin{tikzpicture}[scale=" + str(latex_options['tikz_scale']) + "]\n"
        mark_points = []
        if latex_options['valleys']:
            mark_points.extend(valleys)
        if latex_options['peaks']:
            mark_points.extend(peaks)
        for v in mark_points:
            res += "  \\draw[line width=2,color=red,fill=red] %s circle (%s);\n" % (str(v), 0.15 + .03 * latex_options['line width'])
        if latex_options["bounce path"]:
            D = self.bounce_path()
            D.set_latex_options(latex_options)
            D.set_latex_options({"color": "green",
                                 "line width": 2 * latex_options['line width'],
                                 "bounce path": False,
                                 "peaks": False, "valleys": False})
            res += D._latex_().split("\n")[-2] + "\n"
        for v1, v2 in grid:
            res += "  \\draw[dotted] %s grid %s;\n" % (str(v1), str(v2))
        if diagonal:
            res += "  \\draw (0,0) -- %s;\n" % str((self.number_of_open_symbols(), self.number_of_open_symbols()))
        res += "  \\draw[rounded corners=1, color=%s, line width=%s] (0, 0)" % (latex_options['color'], str(latex_options['line width']))
        next(hti)
        for i, j in hti:
            res += " -- (%s, %s)" % (i, j)
        res += ";\n"
        res += "\\end{tikzpicture}$}}"
        return res

    def _repr_svg_(self) -> str:
        """
        Return the svg picture of ``self``.

        This can be displayed by Jupyter.

        EXAMPLES::

            sage: PP = DyckWords(6).random_element()
            sage: PP._repr_svg_()
            '<?xml...</g></svg>'
        """
        N = self.length()
        width = 0.1 if N < 20 else N / 200
        resu = '<?xml version=\"1.0\" standalone=\"no\"?>'
        resu += '<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" '
        resu += '\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">'
        resu += '<svg xmlns=\"http://www.w3.org/2000/svg\" '
        resu += 'xmlns:xlink=\"http://www.w3.org/1999/xlink\" width=\"500\" viewBox='

        resu1 = '<g style=\"stroke-width:{};stroke-linejoin:bevel; '.format(width)
        resu1 += 'stroke-linecap:butt; stroke:black; fill:snow\">'

        resu3 = '<g style=\"stroke-width:{};stroke-linejoin:bevel;stroke-dasharray:0.25; '.format(width / 2)
        resu3 += 'stroke-linecap:butt; stroke:gray; fill:none\">'

        horizontal = "<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\"/>"
        hori_lines = []
        path = ['<polyline points=\"0,0']
        x, y = 0, 0
        max_y = 0
        last_seen_level = [0]
        for e in self:
            x += 1
            if e == open_symbol:
                y += 1
                last_seen_level.append(x - 1)
                max_y = max(max_y, y)
            else:
                y -= 1
                old_x = last_seen_level.pop()
                hori_lines.append(horizontal.format(old_x, -y, x, -y))
            path.append(f"{x},{-y}")
        path.append('\"/>')
        path.append('</g>')
        resu1 += " ".join(path)
        hori_lines.append('</g></svg>')
        resu3 += "".join(hori_lines)

        margin = 2 * width
        resu += '\"{} {} {} {} \">'.format(-margin, -max_y - margin,
                                           N + 2 * margin, max_y + 2 * margin)

        return resu + resu1 + resu3

    def plot(self, **kwds):
        """
        Plot a Dyck word as a continuous path.

        EXAMPLES::

            sage: w = DyckWords(100).random_element()
            sage: w.plot()                                                              # needs sage.plot
            Graphics object consisting of 1 graphics primitive
        """
        from sage.plot.plot import list_plot
        step = [-1, 1]
        sigma = 0
        list_sigma = [0]
        for l in self:
            sigma += step[l]
            list_sigma.append(sigma)
        return list_plot(list_sigma, plotjoined=True, **kwds)

    def length(self) -> int:
        r"""
        Return the length of ``self``.

        EXAMPLES::

            sage: DyckWord([1, 0, 1, 0]).length()
            4
            sage: DyckWord([1, 0, 1, 1, 0]).length()
            5

        TESTS::

            sage: DyckWord([]).length()
            0
        """
        return len(self)

    def number_of_open_symbols(self) -> int:
        r"""
        Return the number of open symbols in ``self``.

        EXAMPLES::

            sage: DyckWord([1, 0, 1, 0]).number_of_open_symbols()
            2
            sage: DyckWord([1, 0, 1, 1, 0]).number_of_open_symbols()
            3

        TESTS::

            sage: DyckWord([]).number_of_open_symbols()
            0
        """
        return len([x for x in self if x == open_symbol])

    def number_of_close_symbols(self) -> int:
        r"""
        Return the number of close symbols in ``self``.

        EXAMPLES::

            sage: DyckWord([1, 0, 1, 0]).number_of_close_symbols()
            2
            sage: DyckWord([1, 0, 1, 1, 0]).number_of_close_symbols()
            2

        TESTS::

            sage: DyckWord([]).number_of_close_symbols()
            0
        """
        return len([x for x in self if x == close_symbol])

    def is_complete(self) -> bool:
        r"""
        Return ``True`` if ``self`` is complete.

        A Dyck word `d` is complete if `d` contains as many closers as openers.

        EXAMPLES::

            sage: DyckWord([1, 0, 1, 0]).is_complete()
            True
            sage: DyckWord([1, 0, 1, 1, 0]).is_complete()
            False

        TESTS::

            sage: DyckWord([]).is_complete()
            True
        """
        return self.number_of_open_symbols() == self.number_of_close_symbols()

    def height(self) -> int:
        r"""
        Return the height of ``self``.

        We view the Dyck word as a Dyck path from `(0, 0)` to
        `(2n, 0)` in the first quadrant by letting ``1``'s represent
        steps in the direction `(1, 1)` and ``0``'s represent steps in
        the direction `(1, -1)`.

        The height is the maximum `y`-coordinate reached.

        .. SEEALSO:: :meth:`heights`

        EXAMPLES::

            sage: DyckWord([]).height()
            0
            sage: DyckWord([1,0]).height()
            1
            sage: DyckWord([1, 1, 0, 0]).height()
            2
            sage: DyckWord([1, 1, 0, 1, 0]).height()
            2
            sage: DyckWord([1, 1, 0, 0, 1, 0]).height()
            2
            sage: DyckWord([1, 0, 1, 0]).height()
            1
            sage: DyckWord([1, 1, 0, 0, 1, 1, 1, 0, 0, 0]).height()
            3
        """
        # calling max(self.heights()) has a significant overhead (20%)
        height = 0
        height_max = 0
        for letter in self:
            if letter == open_symbol:
                height += 1
                height_max = max(height, height_max)
            elif letter == close_symbol:
                height -= 1
        return height_max

    def heights(self) -> tuple:
        r"""
        Return the heights of ``self``.

        We view the Dyck word as a Dyck path from `(0,0)` to
        `(2n,0)` in the first quadrant by letting ``1``'s represent
        steps in the direction `(1,1)` and ``0``'s represent steps in
        the direction `(1,-1)`.

        The heights is the sequence of the `y`-coordinates of all
        `2n+1` lattice points along the path.

        .. SEEALSO:: :meth:`~DyckWords.from_heights`, :meth:`~DyckWords.min_from_heights`

        EXAMPLES::

            sage: DyckWord([]).heights()
            (0,)
            sage: DyckWord([1,0]).heights()
            (0, 1, 0)
            sage: DyckWord([1, 1, 0, 0]).heights()
            (0, 1, 2, 1, 0)
            sage: DyckWord([1, 1, 0, 1, 0]).heights()
            (0, 1, 2, 1, 2, 1)
            sage: DyckWord([1, 1, 0, 0, 1, 0]).heights()
            (0, 1, 2, 1, 0, 1, 0)
            sage: DyckWord([1, 0, 1, 0]).heights()
            (0, 1, 0, 1, 0)
            sage: DyckWord([1, 1, 0, 0, 1, 1, 1, 0, 0, 0]).heights()
            (0, 1, 2, 1, 0, 1, 2, 3, 2, 1, 0)
        """
        height = 0
        heights = [0] * (len(self) + 1)
        for i, letter in enumerate(self):
            if letter == open_symbol:
                height += 1
            elif letter == close_symbol:
                height -= 1
            heights[i + 1] = height
        return tuple(heights)

    def associated_parenthesis(self, pos) -> int | None:
        r"""
        Report the position for the parenthesis in ``self`` that matches the
        one at position ``pos``.

        The positions in ``self`` are counted from `0`.

        INPUT:

        - ``pos`` -- the index of the parenthesis in the list

        OUTPUT:

        - Integer representing the index of the matching parenthesis.
          If no parenthesis matches, return ``None``.

        EXAMPLES::

            sage: DyckWord([1, 0]).associated_parenthesis(0)
            1
            sage: DyckWord([1, 0, 1, 0]).associated_parenthesis(0)
            1
            sage: DyckWord([1, 0, 1, 0]).associated_parenthesis(1)
            0
            sage: DyckWord([1, 0, 1, 0]).associated_parenthesis(2)
            3
            sage: DyckWord([1, 0, 1, 0]).associated_parenthesis(3)
            2
            sage: DyckWord([1, 1, 0, 0]).associated_parenthesis(0)
            3
            sage: DyckWord([1, 1, 0, 0]).associated_parenthesis(2)
            1
            sage: DyckWord([1, 1, 0]).associated_parenthesis(1)
            2
            sage: DyckWord([1, 1]).associated_parenthesis(0)
        """
        d = 0
        height = 0
        if pos >= len(self):
            raise ValueError("invalid index")

        if self[pos] == open_symbol:
            d += 1
            height += 1
        elif self[pos] == close_symbol:
            d -= 1
            height -= 1
        else:
            raise ValueError("unknown symbol %s" % self[pos])

        while height != 0:
            pos += d
            if pos < 0 or pos >= len(self):
                return None
            if self[pos] == open_symbol:
                height += 1
            elif self[pos] == close_symbol:
                height -= 1
        return pos

    def ascent_prime_decomposition(self) -> list[DyckWord]:
        r"""
        Decompose this Dyck word into a sequence of ascents and prime
        Dyck paths.

        A Dyck word is *prime* if it is complete and has precisely
        one return - the final step.  In particular, the empty Dyck
        path is not prime.  Thus, the factorization is unique.

        This decomposition yields a sequence of odd length: the words
        with even indices consist of up steps only, the words with
        odd indices are prime Dyck paths.  The concatenation of the
        result is the original word.

        EXAMPLES::

            sage: D = DyckWord([1,1,1,0,1,0,1,1,1,1,0,1])
            sage: D.ascent_prime_decomposition()
            [[1, 1], [1, 0], [], [1, 0], [1, 1, 1], [1, 0], [1]]

            sage: DyckWord([]).ascent_prime_decomposition()
            [[]]

            sage: DyckWord([1,1]).ascent_prime_decomposition()
            [[1, 1]]

            sage: DyckWord([1,0,1,0]).ascent_prime_decomposition()
            [[], [1, 0], [], [1, 0], []]
        """
        n = self.length()
        H = self.heights()
        result = []
        i = 0
        height = 0
        up = 0
        while i < n:
            j = i + 1
            while H[j] != height:
                if j == n:
                    i += 1
                    height += 1
                    up += 1
                    break
                j += 1
            else:
                result.extend([DyckWord([open_symbol] * up),  # type:ignore
                               DyckWord(self[i:j])])  # type:ignore
                i = j
                up = 0

        result.append(DyckWord([open_symbol] * up))  # type:ignore
        return result

    def catalan_factorization(self) -> list[DyckWord]:
        r"""
        Decompose this Dyck word into a sequence of complete Dyck
        words.

        Each element of the list returned is a (possibly empty)
        complete Dyck word.  The original word is obtained by placing
        an up step between each of these complete Dyck words.  Thus,
        the number of words returned is one more than the final height.

        See Section 1.2 of [CC1982]_ or Lemma 9.1.1 of [Lot2005]_.

        EXAMPLES::

            sage: D = DyckWord([1,1,1,0,1,0,1,1,1,1,0,1])
            sage: D.catalan_factorization()
            [[], [], [1, 0, 1, 0], [], [], [1, 0], []]

            sage: DyckWord([]).catalan_factorization()
            [[]]

            sage: DyckWord([1,1]).catalan_factorization()
            [[], [], []]

            sage: DyckWord([1,0,1,0]).catalan_factorization()
            [[1, 0, 1, 0]]
        """
        H = self.heights()
        h = 0
        i = 0
        j = n = self.length()
        result = []
        while i <= n:
            if H[j] == h or j == i:
                result.append(DyckWord(self[i:j]))  # type:ignore
                h += 1
                i = j + 1
                j = n
            else:
                j -= 1
        return result

    def number_of_initial_rises(self) -> int:
        r"""
        Return the length of the initial run of ``self``.

        OUTPUT: nonnegative integer indicating the length of the initial rise

        EXAMPLES::

            sage: DyckWord([1, 0, 1, 0]).number_of_initial_rises()
            1
            sage: DyckWord([1, 1, 0, 0]).number_of_initial_rises()
            2
            sage: DyckWord([1, 1, 0, 0, 1, 0]).number_of_initial_rises()
            2
            sage: DyckWord([1, 0, 1, 1, 0, 0]).number_of_initial_rises()
            1

        TESTS::

            sage: DyckWord([]).number_of_initial_rises()
            0
            sage: DyckWord([1, 0]).number_of_initial_rises()
            1
        """
        if not self:
            return 0
        i = 1
        while self[i] == open_symbol:
            i += 1
        return i

    def peaks(self) -> list:
        r"""
        Return a list of the positions of the peaks of a Dyck word.

        A peak is `1` followed by a `0`.  Note that this does not agree with
        the definition given in [Hag2008]_.

        .. SEEALSO:: :meth:`valleys`, :meth:`number_of_peaks`

        EXAMPLES::

            sage: DyckWord([1, 0, 1, 0]).peaks()
            [0, 2]
            sage: DyckWord([1, 1, 0, 0]).peaks()
            [1]
            sage: DyckWord([1,1,0,1,0,1,0,0]).peaks() # Haglund's def gives 2
            [1, 3, 5]
        """
        return [i for i in range(len(self) - 1)
                if self[i] == open_symbol and self[i + 1] == close_symbol]

    def number_of_peaks(self) -> int:
        r"""
        Return the number of peaks of the Dyck path associated to ``self`` .

        .. SEEALSO:: :meth:`peaks`, :meth:`number_of_valleys`

        EXAMPLES::

            sage: DyckWord([1, 0, 1, 0]).number_of_peaks()
            2
            sage: DyckWord([1, 1, 0, 0]).number_of_peaks()
            1
            sage: DyckWord([1,1,0,1,0,1,0,0]).number_of_peaks()
            3
            sage: DyckWord([]).number_of_peaks()
            0
        """
        return len(self.peaks())

    def valleys(self) -> list:
        r"""
        Return a list of the positions of the valleys of a Dyck word.

        A valley is `0` followed by a `1`.

        .. SEEALSO:: :meth:`peaks`, :meth:`number_of_valleys`

        EXAMPLES::

            sage: DyckWord([1, 0, 1, 0]).valleys()
            [1]
            sage: DyckWord([1, 1, 0, 0]).valleys()
            []
            sage: DyckWord([1,1,0,1,0,1,0,0]).valleys()
            [2, 4]
        """
        return [i for i in range(len(self) - 1)
                if self[i] == close_symbol and self[i + 1] == open_symbol]

    def number_of_valleys(self) -> int:
        r"""
        Return the number of valleys of ``self``.

        .. SEEALSO:: :meth:`number_of_peaks`, :meth:`valleys`

        EXAMPLES::

            sage: DyckWord([1, 0, 1, 0]).number_of_valleys()
            1
            sage: DyckWord([1, 1, 0, 0]).number_of_valleys()
            0
            sage: DyckWord([1, 1, 0, 0, 1, 0]).number_of_valleys()
            1
            sage: DyckWord([1, 0, 1, 1, 0, 0]).number_of_valleys()
            1

        TESTS::

            sage: DyckWord([]).number_of_valleys()
            0
            sage: DyckWord([1, 0]).number_of_valleys()
            0
        """
        return len(self.valleys())

    def position_of_first_return(self) -> int:
        r"""
        Return the number of vertical steps before the Dyck path returns to
        the main diagonal.

        EXAMPLES::

            sage: DyckWord([1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0]).position_of_first_return()
            1
            sage: DyckWord([1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0]).position_of_first_return()
            7
            sage: DyckWord([1, 1, 0, 0]).position_of_first_return()
            2
            sage: DyckWord([1, 0, 1, 0]).position_of_first_return()
            1
            sage: DyckWord([]).position_of_first_return()
            0
        """
        touches = self.touch_points()
        if not touches:
            return 0
        else:
            return touches[0]

    def positions_of_double_rises(self) -> list:
        r"""
        Return a list of positions in ``self`` where there are two
        consecutive `1`'s.

        EXAMPLES::

            sage: DyckWord([1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0]).positions_of_double_rises()
            [2, 5]
            sage: DyckWord([1, 1, 0, 0]).positions_of_double_rises()
            [0]
            sage: DyckWord([1, 0, 1, 0]).positions_of_double_rises()
            []
        """
        return [i for i in range(len(self) - 1)
                if self[i] == self[i + 1] == open_symbol]

    def number_of_double_rises(self) -> int:
        r"""
        Return a the number of positions in ``self`` where there are two
        consecutive `1`'s.

        EXAMPLES::

            sage: DyckWord([1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0]).number_of_double_rises()
            2
            sage: DyckWord([1, 1, 0, 0]).number_of_double_rises()
            1
            sage: DyckWord([1, 0, 1, 0]).number_of_double_rises()
            0
        """
        return len(self.positions_of_double_rises())

    def returns_to_zero(self) -> list[int]:
        r"""
        Return a list of positions where ``self`` has height `0`,
        excluding the position `0`.

        EXAMPLES::

            sage: DyckWord([]).returns_to_zero()
            []
            sage: DyckWord([1, 0]).returns_to_zero()
            [2]
            sage: DyckWord([1, 0, 1, 0]).returns_to_zero()
            [2, 4]
            sage: DyckWord([1, 1, 0, 0]).returns_to_zero()
            [4]
        """
        height = 0
        points = []
        for i, letter in enumerate(self):
            if letter == open_symbol:
                height += 1
            elif letter == close_symbol:
                height -= 1
            if not height:
                points.append(i + 1)
        return points

    def touch_points(self) -> list[int]:
        r"""
        Return the abscissae (or, equivalently, ordinates) of the
        points where the Dyck path corresponding to ``self`` (comprising
        `NE` and `SE` steps) touches the main diagonal. This includes
        the last point (if it is on the main diagonal) but excludes the
        beginning point.

        Note that these abscissae are precisely the entries of
        :meth:`returns_to_zero` divided by `2`.

        OUTPUT: list of integers indicating where the path touches the diagonal

        EXAMPLES::

            sage: DyckWord([1, 0, 1, 0]).touch_points()
            [1, 2]
            sage: DyckWord([1, 1, 0, 0]).touch_points()
            [2]
            sage: DyckWord([1, 1, 0, 0, 1, 0]).touch_points()
            [2, 3]
            sage: DyckWord([1, 0, 1, 1, 0, 0]).touch_points()
            [1, 3]
        """
        return [i // 2 for i in self.returns_to_zero()]

    def touch_composition(self):
        r"""
        Return a composition which indicates the positions where ``self``
        returns to the diagonal.

        This assumes ``self`` to be a complete Dyck word.

        OUTPUT: a composition of length equal to the length of the Dyck word

        EXAMPLES::

            sage: DyckWord([1, 0, 1, 0]).touch_composition()
            [1, 1]
            sage: DyckWord([1, 1, 0, 0]).touch_composition()
            [2]
            sage: DyckWord([1, 1, 0, 0, 1, 0]).touch_composition()
            [2, 1]
            sage: DyckWord([1, 0, 1, 1, 0, 0]).touch_composition()
            [1, 2]
            sage: DyckWord([]).touch_composition()
            []
        """
        from sage.combinat.composition import Composition
        if not self:
            return Composition([])
        return Composition(descents=[i - 1 for i in self.touch_points()])

    def number_of_touch_points(self) -> int:
        r"""
        Return the number of touches of ``self`` at the main diagonal.

        OUTPUT: nonnegative integer

        EXAMPLES::

            sage: DyckWord([1, 0, 1, 0]).number_of_touch_points()
            2
            sage: DyckWord([1, 1, 0, 0]).number_of_touch_points()
            1
            sage: DyckWord([1, 1, 0, 0, 1, 0]).number_of_touch_points()
            2
            sage: DyckWord([1, 0, 1, 1, 0, 0]).number_of_touch_points()
            2

        TESTS::

            sage: DyckWord([]).number_of_touch_points()
            0
        """
        return len(self.returns_to_zero())

    def rise_composition(self):
        r"""
        The sequences of lengths of runs of `1`'s in ``self``. Also equal to
        the sequence of lengths of vertical segments in the Dyck path.

        EXAMPLES::

            sage: DyckWord([1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0]).pretty_print()
                       ___
                      | x
               _______|  .
              | x x x  . .
              | x x  . . .
             _| x  . . . .
            | x  . . . . .
            |  . . . . . .

            sage: DyckWord([1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0]).rise_composition()
            [2, 3, 2]
            sage: DyckWord([1,1,0,0]).rise_composition()
            [2]
            sage: DyckWord([1,0,1,0]).rise_composition()
            [1, 1]
        """
        from sage.combinat.composition import Composition
        L = list(self)
        rise_comp = []
        while L:
            i = L.index(0)
            L = L[i + 1:]
            if i:
                rise_comp.append(i)
        return Composition(rise_comp)

    @combinatorial_map(name='to two-row standard tableau')
    def to_standard_tableau(self):
        r"""
        Return a standard tableau of shape `(a,b)` where
        `a` is the number of open symbols and `b` is the number of
        close symbols in ``self``.

        EXAMPLES::

            sage: DyckWord([]).to_standard_tableau()
            []
            sage: DyckWord([1, 0]).to_standard_tableau()
            [[1], [2]]
            sage: DyckWord([1, 1, 0, 0]).to_standard_tableau()
            [[1, 2], [3, 4]]
            sage: DyckWord([1, 0, 1, 0]).to_standard_tableau()
            [[1, 3], [2, 4]]
            sage: DyckWord([1]).to_standard_tableau()
            [[1]]
            sage: DyckWord([1, 0, 1]).to_standard_tableau()
            [[1, 3], [2]]
        """
        open_positions = []
        close_positions = []
        for i in range(len(self)):
            if self[i] == open_symbol:
                open_positions.append(i + 1)
            else:
                close_positions.append(i + 1)
        from sage.combinat.tableau import StandardTableau
        return StandardTableau([x for x in [open_positions, close_positions] if x])

    def to_tamari_sorting_tuple(self) -> list[int]:
        """
        Convert a Dyck word to a Tamari sorting tuple.

        The result is a list of integers, one for every up-step from
        left to right. To each up-step is associated
        the distance to the corresponding down step in the Dyck word.

        This is useful for a faster conversion to binary trees.

        EXAMPLES::

            sage: DyckWord([]).to_tamari_sorting_tuple()
            []
            sage: DyckWord([1, 0]).to_tamari_sorting_tuple()
            [0]
            sage: DyckWord([1, 1, 0, 0]).to_tamari_sorting_tuple()
            [1, 0]
            sage: DyckWord([1, 0, 1, 0]).to_tamari_sorting_tuple()
            [0, 0]
            sage: DyckWord([1, 1, 0, 1, 0, 0]).to_tamari_sorting_tuple()
            [2, 0, 0]

        .. SEEALSO:: :meth:`~DyckWord_complete.to_Catalan_code`
        """
        position = 0
        resu = [-i - 1 for i in range(len(self) // 2)]
        indices_of_active_ups = []
        for letter in self:
            if letter == open_symbol:
                indices_of_active_ups.append(position)
                position += 1
            else:
                previous = indices_of_active_ups.pop()
                resu[previous] += position
        return resu

    @combinatorial_map(name="to binary trees: up step, left tree, down step, right tree")
    def to_binary_tree(self, usemap='1L0R'):
        r"""
        Return a binary tree recursively constructed from the Dyck path
        ``self`` by the map ``usemap``. The default ``usemap`` is ``'1L0R'``
        which means:

        - an empty Dyck word is a leaf,

        - a non empty Dyck word reads `1 L 0 R` where `L` and `R` correspond
          to respectively its left and right subtrees.

        INPUT:

        - ``usemap`` -- string, either ``'1L0R'``, ``'1R0L'``, ``'L1R0'``,
          ``'R1L0'``

        Other valid ``usemap`` are ``'1R0L'``, ``'L1R0'``, and ``'R1L0'``.
        These correspond to different maps from Dyck paths to binary
        trees, whose recursive definitions are hopefully clear from the
        names.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: dw = DyckWord([1,0])
            sage: dw.to_binary_tree()
            [., .]
            sage: dw = DyckWord([])
            sage: dw.to_binary_tree()
            .
            sage: dw = DyckWord([1,0,1,1,0,0])
            sage: dw.to_binary_tree()
            [., [[., .], .]]
            sage: dw.to_binary_tree("L1R0")
            [[., .], [., .]]
            sage: dw = DyckWord([1,0,1,1,0,0,1,1,1,0,1,0,0,0])
            sage: dw.to_binary_tree() == dw.to_binary_tree("1R0L").left_right_symmetry()
            True
            sage: dw.to_binary_tree() == dw.to_binary_tree("L1R0").left_border_symmetry()
            False
            sage: dw.to_binary_tree("1R0L") == dw.to_binary_tree("L1R0").left_border_symmetry()
            True
            sage: dw.to_binary_tree("R1L0") == dw.to_binary_tree("L1R0").left_right_symmetry()
            True
            sage: dw.to_binary_tree("R10L")
            Traceback (most recent call last):
            ...
            ValueError: R10L is not a correct map
        """
        if usemap not in ["1L0R", "1R0L", "L1R0", "R1L0"]:
            raise ValueError("%s is not a correct map" % usemap)
        from sage.combinat.binary_tree import BinaryTree
        if not self:
            return BinaryTree()
        tp = [0]
        tp.extend(self.returns_to_zero())
        l = len(self)
        if usemap[0] == '1':  # we check what kind of reduction we want
            s0 = 1  # start point for first subtree
            e0 = tp[1] - 1  # end point for first subtree
            s1 = e0 + 1  # start point for second subtree
            e1 = l  # end point for second subtree
        else:
            s0 = 0
            e0 = tp[len(tp) - 2]
            s1 = e0 + 1
            e1 = l - 1
        trees = [DyckWord(self[s0:e0]).to_binary_tree(usemap),
                 DyckWord(self[s1:e1]).to_binary_tree(usemap)]
        if usemap[0] == "R" or usemap[1] == "R":
            trees.reverse()
        return BinaryTree(trees)

    @combinatorial_map(name="to the Tamari corresponding Binary tree")
    def to_binary_tree_tamari(self):
        r"""
        Return the binary tree corresponding to ``self`` in a way which
        is consistent with the Tamari orders on the set of Dyck paths and
        on the set of binary trees.

        This is the ``'L1R0'`` map documented in :meth:`to_binary_tree`.

        EXAMPLES::

            sage: DyckWord([1,0]).to_binary_tree_tamari()                               # needs sage.graphs
            [., .]
            sage: DyckWord([1,0,1,1,0,0]).to_binary_tree_tamari()                       # needs sage.graphs
            [[., .], [., .]]
            sage: DyckWord([1,0,1,0,1,0]).to_binary_tree_tamari()                       # needs sage.graphs
            [[[., .], .], .]
        """
        # return self.to_binary_tree("L1R0")  # slower and recursive
        from sage.combinat.binary_tree import from_tamari_sorting_tuple
        tup = self.to_tamari_sorting_tuple()
        return from_tamari_sorting_tuple(tup)

    def tamari_interval(self, other):
        r"""
        Return the Tamari interval between ``self`` and ``other`` as a
        :class:`~sage.combinat.interval_posets.TamariIntervalPoset`.

        A "Tamari interval" means an interval in the Tamari order. The
        Tamari order on the set of Dyck words of size `n` is the
        partial order obtained from the Tamari order on the set of
        binary trees of size `n` (see
        :meth:`~sage.combinat.binary_tree.BinaryTree.tamari_lequal`)
        by means of the Tamari bijection between Dyck words and binary
        trees
        (:meth:`~sage.combinat.binary_tree.BinaryTree.to_dyck_word_tamari`).

        INPUT:

        - ``other`` -- a Dyck word greater or equal to ``self`` in the
          Tamari order

        EXAMPLES::

            sage: # needs sage.graphs
            sage: dw = DyckWord([1, 1, 0, 1, 0, 0, 1, 0])
            sage: ip = dw.tamari_interval(DyckWord([1, 1, 1, 0, 0, 1, 0, 0])); ip
            The Tamari interval of size 4 induced by relations [(2, 4), (3, 4), (3, 1), (2, 1)]
            sage: ip.lower_dyck_word()
            [1, 1, 0, 1, 0, 0, 1, 0]
            sage: ip.upper_dyck_word()
            [1, 1, 1, 0, 0, 1, 0, 0]
            sage: ip.interval_cardinality()
            4
            sage: ip.number_of_tamari_inversions()
            2
            sage: list(ip.dyck_words())
            [[1, 1, 1, 0, 0, 1, 0, 0],
             [1, 1, 1, 0, 0, 0, 1, 0],
             [1, 1, 0, 1, 0, 1, 0, 0],
             [1, 1, 0, 1, 0, 0, 1, 0]]
            sage: dw.tamari_interval(DyckWord([1,1,0,0,1,1,0,0]))
            Traceback (most recent call last):
            ...
            ValueError: the two Dyck words are not comparable on the Tamari lattice
        """
        from sage.combinat.interval_posets import TamariIntervalPosets
        return TamariIntervalPosets.from_dyck_words(self, other)

    def _area_sequence_iter(self) -> Iterator[int]:
        """
        Return an iterator producing the area sequence.

        .. SEEALSO:: :meth:`to_area_sequence`

        EXAMPLES::

            sage: d = DyckWord([1, 0, 1, 0])
            sage: [a for a in d._area_sequence_iter()]
            [0, 0]
        """
        a = 0
        for move in self:
            if move == open_symbol:
                yield a
                a += 1
            else:
                a -= 1

    def to_area_sequence(self) -> list[int]:
        r"""
        Return the area sequence of the Dyck word ``self``.

        The area sequence of a Dyck word `w` is defined as follows:
        Representing the Dyck word `w` as a Dyck path from `(0, 0)` to
        `(n, n)` using `N` and `E` steps (this involves padding `w` by
        `E` steps until `w` reaches the main diagonal if `w` is not
        already a complete Dyck path), the area sequence of `w` is the
        sequence `(a_1, a_2, \ldots, a_n)`, where `a_i` is the number
        of full cells in the `i`-th row of the rectangle
        `[0, n] \times [0, n]` which lie completely above the diagonal.
        (The cells are the regions into which the rectangle is
        subdivided by the lines `x = i` with `i` integer and the lines
        `y = j` with `j` integer. The `i`-th row consists of all the
        cells between the lines `y = i-1` and `y = i`.)

        An alternative definition:
        Representing the Dyck word `w` as a Dyck path consisting of
        `NE` and `SE` steps, the area sequence is the sequence of
        ordinates of all lattice points on the path which are
        starting points of `NE` steps.

        A list of integers `l` is the area sequence of some Dyck path
        if and only if it satisfies `l_0 = 0` and
        `0 \leq l_{i+1} \leq l_i + 1` for `i > 0`.

        EXAMPLES::

            sage: DyckWord([]).to_area_sequence()
            []
            sage: DyckWord([1, 0]).to_area_sequence()
            [0]
            sage: DyckWord([1, 1, 0, 0]).to_area_sequence()
            [0, 1]
            sage: DyckWord([1, 0, 1, 0]).to_area_sequence()
            [0, 0]
            sage: all(dw ==
            ....:     DyckWords().from_area_sequence(dw.to_area_sequence())
            ....:     for i in range(6) for dw in DyckWords(i))
            True
            sage: DyckWord([1,0,1,0,1,0,1,0,1,0]).to_area_sequence()
            [0, 0, 0, 0, 0]
            sage: DyckWord([1,1,1,1,1,0,0,0,0,0]).to_area_sequence()
            [0, 1, 2, 3, 4]
            sage: DyckWord([1,1,1,1,0,1,0,0,0,0]).to_area_sequence()
            [0, 1, 2, 3, 3]
            sage: DyckWord([1,1,0,1,0,0,1,1,0,1,0,1,0,0]).to_area_sequence()
            [0, 1, 1, 0, 1, 1, 1]
        """
        return list(self._area_sequence_iter())


class DyckWord_complete(DyckWord):
    r"""
    The class of complete
    :class:`Dyck words<sage.combinat.dyck_word.DyckWord>`.
    A Dyck word is complete, if it contains as many closers as openers.

    For further information on Dyck words, see
    :class:`DyckWords_class<sage.combinat.dyck_word.DyckWord>`.
    """

    def semilength(self) -> int:
        r"""
        Return the semilength of ``self``.

        The semilength of a complete Dyck word `d` is the number of openers
        and the number of closers.

        EXAMPLES::

            sage: DyckWord([1, 0, 1, 0]).semilength()
            2

        TESTS::

            sage: DyckWord([]).semilength()
            0
        """
        return len(self) // 2

    @combinatorial_map(name='to partition')
    def to_partition(self):
        r"""
        Return the partition associated to ``self`` .

        This partition is determined by thinking of ``self`` as a lattice path
        and considering the cells which are above the path but within the
        `n \times n`  grid and the partition is formed by reading the sequence
        of the number of cells in this collection in each row.

        OUTPUT:

        - a partition representing the rows of cells in the square lattice
          and above the path

        EXAMPLES::

            sage: DyckWord([]).to_partition()
            []
            sage: DyckWord([1,0]).to_partition()
            []
            sage: DyckWord([1,1,0,0]).to_partition()
            []
            sage: DyckWord([1,0,1,0]).to_partition()
            [1]
            sage: DyckWord([1,0,1,0,1,0]).to_partition()
            [2, 1]
            sage: DyckWord([1,1,0,0,1,0]).to_partition()
            [2]
            sage: DyckWord([1,0,1,1,0,0]).to_partition()
            [1, 1]
        """
        from sage.combinat.partition import Partition
        n = len(self) // 2
        res = []
        for c in reversed(self):
            if c == close_symbol:
                n -= 1
            else:
                res.append(n)
        return Partition(res)

    def number_of_parking_functions(self) -> int:
        r"""
        Return the number of parking functions with ``self`` as the supporting
        Dyck path.

        One representation of a parking function is as a pair consisting of a
        Dyck path and a permutation `\pi` such that if
        `[a_0, a_1, \ldots, a_{n-1}]` is the area_sequence of the Dyck path
        (see :meth:`to_area_sequence<DyckWord.to_area_sequence>`) then the
        permutation `\pi` satisfies `\pi_i < \pi_{i+1}` whenever
        `a_{i} < a_{i+1}`. This function counts the number of permutations `\pi`
        which satisfy this condition.

        EXAMPLES::

            sage: DyckWord(area_sequence=[0,1,2]).number_of_parking_functions()
            1
            sage: DyckWord(area_sequence=[0,1,1]).number_of_parking_functions()
            3
            sage: DyckWord(area_sequence=[0,1,0]).number_of_parking_functions()
            3
            sage: DyckWord(area_sequence=[0,0,0]).number_of_parking_functions()
            6
        """
        from sage.arith.misc import multinomial
        return multinomial(self.rise_composition())

    def list_parking_functions(self) -> list:
        r"""
        Return all parking functions whose supporting Dyck path is ``self``.

        EXAMPLES::

            sage: DyckWord([1,1,0,0,1,0]).list_parking_functions()
            [[1, 1, 3], [1, 3, 1], [3, 1, 1]]
        """
        return list(self.parking_functions())

    def parking_functions(self):
        r"""
        Iterate over parking functions whose supporting Dyck path is ``self``.

        EXAMPLES::

            sage: list(DyckWord([1,1,0,1,0,0]).parking_functions())
            [[1, 1, 2], [1, 2, 1], [2, 1, 1]]
        """
        from sage.combinat.parking_functions import ParkingFunction
        alist = self._area_sequence_iter()
        for pi in Permutations([i - ai + 1 for i, ai in enumerate(alist)]):
            yield ParkingFunction(pi)

    def reading_permutation(self) -> Permutation:
        r"""
        Return the reading permutation of ``self``.

        This is the permutation formed by taking the reading word of
        the Dyck path representing ``self`` (with `N` and `E` steps)
        if the vertical edges of the Dyck path are labeled from bottom
        to top with `1` through `n` and the diagonals are read from
        top to bottom starting with the diagonal furthest from the
        main diagonal.

        EXAMPLES::

            sage: DyckWord([1,0,1,0]).reading_permutation()
            [2, 1]
            sage: DyckWord([1,1,0,0]).reading_permutation()
            [2, 1]
            sage: DyckWord([1,1,0,1,0,0]).reading_permutation()
            [3, 2, 1]
            sage: DyckWord([1,1,0,0,1,0]).reading_permutation()
            [2, 3, 1]
            sage: DyckWord([1,0,1,1,0,0,1,0]).reading_permutation()
            [3, 4, 2, 1]
        """
        if not self:
            return Permutation([])  # type:ignore
        alist = self.to_area_sequence()
        m = max(alist)
        p1 = Word([m - alist[-i - 1]
                   for i in range(len(alist))]).standard_permutation()
        return p1.inverse().complement()

    def characteristic_symmetric_function(self, q=None,
                                          R=QQ['q', 't'].fraction_field()):
        r"""
        The characteristic function of ``self`` is the sum of
        `q^{dinv(D,F)} Q_{ides(read(D,F))}` over all permutation
        fillings of the Dyck path representing ``self``, where
        `ides(read(D,F))` is the descent composition of the inverse of the
        reading word of the filling.

        INPUT:

        - ``q`` -- (default: ``q = R('q')``) a parameter for the generating
          function power

        - ``R`` -- (default: ``R = QQ['q','t'].fraction_field()``) the base
          ring to do the calculations over

        OUTPUT:

        - an element of the symmetric functions over the ring ``R``
          (in the Schur basis).

        EXAMPLES::

            sage: R = QQ['q','t'].fraction_field()
            sage: (q,t) = R.gens()
            sage: f = sum(t**D.area() * D.characteristic_symmetric_function()           # needs sage.modules
            ....:         for D in DyckWords(3)); f
            (q^3+q^2*t+q*t^2+t^3+q*t)*s[1, 1, 1] + (q^2+q*t+t^2+q+t)*s[2, 1] + s[3]
            sage: f.nabla(power=-1)                                                     # needs sage.modules
            s[1, 1, 1]
        """
        from sage.combinat.ncsf_qsym.qsym import QuasiSymmetricFunctions
        from sage.combinat.sf.sf import SymmetricFunctions
        if q is None:
            q = R('q')
        else:
            if q not in R:
                raise ValueError("q=%s must be an element of the base ring %s" % (q, R))
        F = QuasiSymmetricFunctions(R).Fundamental()
        p = self.reading_permutation().inverse()
        perms = [Word(perm).standard_permutation()
                 for perm in self.list_parking_functions()]
        QSexpr = sum(q**self.dinv(pv.inverse()) * F(Permutation([p(i) for i in pv]).descents_composition()) for pv in perms)
        s = SymmetricFunctions(R).s()
        return s(QSexpr.to_symmetric_function())

    def to_pair_of_standard_tableaux(self) -> tuple:
        r"""
        Convert ``self`` to a pair of standard tableaux of the same shape and
        of length less than or equal to two.

        EXAMPLES::

            sage: DyckWord([1,0,1,0]).to_pair_of_standard_tableaux()
            ([[1], [2]], [[1], [2]])
            sage: DyckWord([1,1,0,0]).to_pair_of_standard_tableaux()
            ([[1, 2]], [[1, 2]])
            sage: DyckWord([1,1,0,1,0,0,1,1,0,1,0,1,0,0]).to_pair_of_standard_tableaux()
            ([[1, 2, 4, 7], [3, 5, 6]], [[1, 2, 4, 6], [3, 5, 7]])
        """
        from sage.combinat.tableau import Tableau
        n = self.semilength()
        if n == 0:
            return (Tableau([]), Tableau([]))  # type:ignore
        elif self.height() == n:
            T = Tableau([list(range(1, n + 1))])  # type:ignore
            return (T, T)
        else:
            left: list[list[int]] = [[], []]
            right: list[list[int]] = [[], []]
            for pos in range(n):
                if self[pos] == open_symbol:
                    left[0].append(pos + 1)
                else:
                    left[1].append(pos + 1)
                if self[-pos - 1] == close_symbol:
                    right[0].append(pos + 1)
                else:
                    right[1].append(pos + 1)
            return (Tableau(left), Tableau(right))  # type:ignore

    @combinatorial_map(name='to 312 avoiding permutation')
    def to_312_avoiding_permutation(self) -> Permutation:
        r"""
        Convert ``self`` to a `312`-avoiding permutation using the bijection
        by Bandlow and Killpatrick in [BK2001]_.

        This sends the area to the inversion number.

        EXAMPLES::

            sage: DyckWord([1,1,0,0]).to_312_avoiding_permutation()
            [2, 1]
            sage: DyckWord([1,0,1,0]).to_312_avoiding_permutation()
            [1, 2]
            sage: p = DyckWord([1,1,0,1,0,0,1,1,0,1,0,1,0,0]).to_312_avoiding_permutation(); p
            [2, 3, 1, 5, 6, 7, 4]
            sage: DyckWord([1,1,0,1,0,0,1,1,0,1,0,1,0,0]).area()
            5
            sage: p.length()
            5

        TESTS::

            sage: PD = [D.to_312_avoiding_permutation() for D in DyckWords(5)]
            sage: all(pi.avoids([3,1,2]) for pi in PD)
            True
            sage: all(D.area()==D.to_312_avoiding_permutation().length() for D in DyckWords(5))
            True
        """
        n = self.semilength()
        area = self._area_sequence_iter()
        pi = Permutations(n).one()
        for j, aj in enumerate(area):
            for i in range(aj):
                pi = pi.apply_simple_reflection(j - i)
        return pi

    @combinatorial_map(name='to non-crossing permutation')
    def to_noncrossing_permutation(self) -> Permutation:
        r"""
        Use the bijection by C. Stump in [Stu2008]_ to send ``self`` to a
        non-crossing permutation.

        A non-crossing permutation when written in cyclic notation has cycles
        which are strictly increasing.  Sends the area to the inversion number
        and ``self.major_index()`` to `n(n-1) - maj(\sigma) - maj(\sigma^{-1})`.
        Uses the function :func:`~sage.combinat.dyck_word.pealing`

        EXAMPLES::

            sage: DyckWord([1,1,0,0]).to_noncrossing_permutation()
            [2, 1]
            sage: DyckWord([1,0,1,0]).to_noncrossing_permutation()
            [1, 2]
            sage: p = DyckWord([1,1,0,1,0,0,1,1,0,1,0,1,0,0]).to_noncrossing_permutation(); p
            [2, 3, 1, 5, 6, 7, 4]
            sage: DyckWord([1,1,0,1,0,0,1,1,0,1,0,1,0,0]).area()
            5
            sage: p.length()
            5

        TESTS::

            sage: all(D.area()==D.to_noncrossing_permutation().length() for D in DyckWords(5))
            True
            sage: all(20-D.major_index()==D.to_noncrossing_permutation().major_index()
            ....:     +D.to_noncrossing_permutation().imajor_index() for D in DyckWords(5))
            True
        """
        n = self.semilength()
        if n == 0:
            return Permutation([])  # type:ignore
        D, touch_sequence = pealing(self, return_touches=True)
        pi = list(range(1, n + 1))
        while touch_sequence:
            for touches in touch_sequence:
                a = pi[touches[0] - 1]
                for i in range(len(touches) - 1):
                    pi[touches[i] - 1] = pi[touches[i + 1] - 1]
                pi[touches[-1] - 1] = a
            D, touch_sequence = pealing(D, return_touches=True)
        return Permutations()(pi, check=False)

    @combinatorial_map(name='to 321 avoiding permutation')
    def to_321_avoiding_permutation(self) -> Permutation:
        r"""
        Use the bijection (pp. 60-61 of [Knu1973]_ or section 3.1 of [CK2008]_)
        to send ``self`` to a `321`-avoiding permutation.

        It is shown in [EP2004]_ that it sends the number of centered tunnels
        to the number of fixed points, the number of right tunnels to the
        number of excedences, and the semilength plus the height of the middle
        point to 2 times the length of the longest increasing subsequence.

        EXAMPLES::

            sage: DyckWord([1,0,1,0]).to_321_avoiding_permutation()
            [2, 1]
            sage: DyckWord([1,1,0,0]).to_321_avoiding_permutation()
            [1, 2]
            sage: D = DyckWord([1,1,0,1,0,0,1,1,0,1,0,1,0,0])
            sage: p = D.to_321_avoiding_permutation()
            sage: p
            [3, 5, 1, 6, 2, 7, 4]
            sage: D.number_of_tunnels()
            0
            sage: p.number_of_fixed_points()
            0
            sage: D.number_of_tunnels('right')
            4
            sage: len(p.weak_excedences())-p.number_of_fixed_points()
            4
            sage: n = D.semilength()
            sage: D.heights()[n] + n
            8
            sage: 2*p.longest_increasing_subsequence_length()
            8

        TESTS::

            sage: PD = [D.to_321_avoiding_permutation() for D in DyckWords(5)]
            sage: all(pi.avoids([3,2,1]) for pi in PD)
            True
            sage: to_perm = lambda x: x.to_321_avoiding_permutation()
            sage: all(D.number_of_tunnels() == to_perm(D).number_of_fixed_points()
            ....:     for D in DyckWords(5))
            True
            sage: all(D.number_of_tunnels('right') == len(to_perm(D).weak_excedences())
            ....:     -to_perm(D).number_of_fixed_points() for D in DyckWords(5))
            True
            sage: all(D.heights()[5]+5 == 2*to_perm(D).longest_increasing_subsequence_length()
            ....:     for D in DyckWords(5))
            True
        """
        from sage.combinat.rsk import RSK_inverse
        A, B = self.to_pair_of_standard_tableaux()
        return RSK_inverse(A, B, output='permutation')

    @combinatorial_map(name='to 132 avoiding permutation')
    def to_132_avoiding_permutation(self) -> Permutation:
        r"""
        Use the bijection by C. Krattenthaler in [Kra2001]_ to send ``self``
        to a `132`-avoiding permutation.

        EXAMPLES::

            sage: DyckWord([1,1,0,0]).to_132_avoiding_permutation()
            [1, 2]
            sage: DyckWord([1,0,1,0]).to_132_avoiding_permutation()
            [2, 1]
            sage: DyckWord([1,1,0,1,0,0,1,1,0,1,0,1,0,0]).to_132_avoiding_permutation()
            [6, 5, 4, 7, 2, 1, 3]

        TESTS::

            sage: PD = [D.to_132_avoiding_permutation() for D in DyckWords(5)]
            sage: all(pi.avoids([1,3,2]) for pi in PD)
            True
        """
        n = self.semilength()
        area = self.to_area_sequence()
        area.append(0)
        pi = []
        values = list(range(1, n + 1))
        for i in range(n):
            if area[n - i - 1] + 1 > area[n - i]:
                pi.append(n - i - area[n - i - 1])
                values.remove(n - i - area[n - i - 1])
            else:
                v = min(v for v in values if v > n - i - area[n - i - 1])
                pi.append(v)
                values.remove(v)
        return Permutation(pi)  # type:ignore

    def to_permutation(self, map) -> Permutation:
        r"""
        This is simply a method collecting all implemented maps from Dyck
        words to permutations.

        INPUT:

        - ``map`` -- defines the map from Dyck words to permutations.
          These are currently:

          - ``Bandlow-Killpatrick``: :func:`to_312_avoiding_permutation`
          - ``Knuth``: :func:`to_321_avoiding_permutation`
          - ``Krattenthaler``: :func:`to_132_avoiding_permutation`
          - ``Stump``: :func:`to_noncrossing_permutation`

        EXAMPLES::

            sage: D = DyckWord([1,1,1,0,1,0,0,0])
            sage: D.pretty_print()
               _____
             _| x x
            | x x  .
            | x  . .
            |  . . .

            sage: D.to_permutation(map='Bandlow-Killpatrick')
            [3, 4, 2, 1]
            sage: D.to_permutation(map='Stump')
            [4, 2, 3, 1]
            sage: D.to_permutation(map='Knuth')
            [1, 2, 4, 3]
            sage: D.to_permutation(map='Krattenthaler')
            [2, 1, 3, 4]

        TESTS::

            sage: D = DyckWord([1,0,1,0])
            sage: D.to_permutation('Banana')
            Traceback (most recent call last):
            ...
            ValueError: the given map is not valid
        """
        if map == "Bandlow-Killpatrick":
            return self.to_312_avoiding_permutation()
        elif map == "Knuth":
            return self.to_321_avoiding_permutation()
        elif map == "Krattenthaler":
            return self.to_132_avoiding_permutation()
        elif map == "Stump":
            return self.to_noncrossing_permutation()
        else:
            raise ValueError("the given map is not valid")

    def to_noncrossing_partition(self, bijection=None):
        r"""
        Bijection of Biane from ``self`` to a noncrossing partition.

        There is an optional parameter ``bijection`` that indicates if a
        different bijection from Dyck words to non-crossing partitions
        should be used (since there are potentially many).

        If the parameter ``bijection`` is "Stump" then the bijection used is
        from [Stu2008]_, see also the method :meth:`to_noncrossing_permutation`.

        Thanks to Mathieu Dutour for describing the bijection.  See also
        :func:`~CompleteDyckWords.from_noncrossing_partition`.

        EXAMPLES::

            sage: DyckWord([]).to_noncrossing_partition()
            {}
            sage: DyckWord([1, 0]).to_noncrossing_partition()
            {{1}}
            sage: DyckWord([1, 1, 0, 0]).to_noncrossing_partition()
            {{1, 2}}
            sage: DyckWord([1, 1, 1, 0, 0, 0]).to_noncrossing_partition()
            {{1, 2, 3}}
            sage: DyckWord([1, 0, 1, 0, 1, 0]).to_noncrossing_partition()
            {{1}, {2}, {3}}
            sage: DyckWord([1, 1, 0, 1, 0, 0]).to_noncrossing_partition()
            {{1, 3}, {2}}
            sage: DyckWord([]).to_noncrossing_partition("Stump")
            {}
            sage: DyckWord([1, 0]).to_noncrossing_partition("Stump")
            {{1}}
            sage: DyckWord([1, 1, 0, 0]).to_noncrossing_partition("Stump")
            {{1, 2}}
            sage: DyckWord([1, 1, 1, 0, 0, 0]).to_noncrossing_partition("Stump")
            {{1, 3}, {2}}
            sage: DyckWord([1, 0, 1, 0, 1, 0]).to_noncrossing_partition("Stump")
            {{1}, {2}, {3}}
            sage: DyckWord([1, 1, 0, 1, 0, 0]).to_noncrossing_partition("Stump")
            {{1, 2, 3}}
        """
        P = SetPartitions(len(self) // 2)
        if bijection == "Stump":
            return P(self.to_noncrossing_permutation().cycle_tuples(), check=False)
        partition = []
        stack = []
        i = 0
        p = 1

        # Invariants:
        # - self[i] = 0
        # - p is the number of opening parens at position i

        while i < len(self):
            stack.append(p)
            j = i + 1
            while j < len(self) and self[j] == close_symbol:
                j += 1

            # Now j points to the next 1 or past the end of self
            nz = j - (i + 1)  # the number of )'s between i and j
            if nz > 0:
                # Remove the nz last elements of stack and
                # make a new part in partition
                if nz > len(stack):
                    raise ValueError("incorrect Dyck word")

                partition.append(stack[-nz:])

                stack = stack[: -nz]
            i = j
            p += 1

        if stack:
            raise ValueError("incorrect Dyck word")

        return P(partition, check=False)

    def to_Catalan_code(self) -> list:
        r"""
        Return the Catalan code associated to ``self``.

        A Catalan code of length `n` is a sequence
        `(a_1, a_2, \ldots, a_n)` of `n` integers `a_i` such that:

        - `0 \leq a_i \leq n-i` for every `i`;

        - if `i < j` and `a_i > 0` and `a_j > 0` and
          `a_{i+1} = a_{i+2} = \cdots = a_{j-1} = 0`,
          then `a_i - a_j < j-i`.

        It turns out that the Catalan codes of length `n` are in
        bijection with Dyck words.

        The Catalan code of a Dyck word is example (x) in Richard Stanley's
        exercises on combinatorial interpretations for Catalan objects.
        The code in this example is the reverse of the description provided
        there. See [Sta-EC2]_ and [StaCat98]_.

        EXAMPLES::

            sage: DyckWord([]).to_Catalan_code()
            []
            sage: DyckWord([1, 0]).to_Catalan_code()
            [0]
            sage: DyckWord([1, 1, 0, 0]).to_Catalan_code()
            [0, 1]
            sage: DyckWord([1, 0, 1, 0]).to_Catalan_code()
            [0, 0]
            sage: all(dw ==
            ....:     DyckWords().from_Catalan_code(dw.to_Catalan_code())
            ....:     for i in range(6) for dw in DyckWords(i))
            True

        .. SEEALSO:: :meth:`to_tamari_sorting_tuple`
        """
        if not self:
            return []
        cut = self.associated_parenthesis(0)
        if cut is None:
            raise ValueError('not valid for incomplete Dyck words')
        recdw = DyckWord(self[1:cut] + self[cut + 1:])  # type:ignore
        returns = [0] + recdw.returns_to_zero()
        res = recdw.to_Catalan_code()
        res.append(returns.index(cut - 1))
        return res

    @combinatorial_map(name="To Ordered tree")
    def to_ordered_tree(self):
        r"""
        Return the ordered tree corresponding to ``self`` where the depth
        of the tree is the maximal height of ``self``.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: D = DyckWord([1,1,0,0])
            sage: D.to_ordered_tree()
            [[[]]]
            sage: D = DyckWord([1,0,1,0])
            sage: D.to_ordered_tree()
            [[], []]
            sage: D = DyckWord([1, 0, 1, 1, 0, 0])
            sage: D.to_ordered_tree()
            [[], [[]]]
            sage: D = DyckWord([1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0])
            sage: D.to_ordered_tree()
            [[], [[], []], [[], [[]]]]

        TESTS::

            sage: D = DyckWord([1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0])
            sage: D == D.to_ordered_tree().to_dyck_word()                               # needs sage.graphs
            True
        """
        from sage.combinat.ordered_tree import OrderedTree
        levels = [OrderedTree().clone()]
        for u in self:
            if u == 1:
                levels.append(OrderedTree().clone())
            else:
                tree = levels.pop()
                tree.set_immutable()
                root = levels.pop()
                root.append(tree)
                levels.append(root)
        root = levels[0]
        root.set_immutable()
        return root

    def to_triangulation(self) -> list:
        r"""
        Map ``self`` to a triangulation.

        The map from complete Dyck words of length `2n` to
        triangulations of `n+2`-gon given by this function is a
        bijection that can be described as follows.

        Consider the Dyck word as a path from `(0, 0)` to `(n, n)`
        staying above the diagonal, where `1` is an up step and `0` is
        a right step. Then each horizontal step has a co-height (`0`
        at the top and `n-1` at most at the bottom). One reads the
        Dyck word from left to right. At the beginning, all vertices
        from `0` to `n+1` are available. For each horizontal step,
        one creates an edge from the vertex indexed by the co-height
        to the next available vertex. This chops out a triangle from
        the polygon and one removes the middle vertex of this triangle
        from the list of available vertices.

        This bijection has the property that the set of smallest
        vertices of the edges in a triangulation is an encoding of the
        co-heights, from which the Dyck word can be easily recovered.

        OUTPUT:

        a list of pairs `(i, j)` that are the edges of the
        triangulations.

        EXAMPLES::

            sage: DyckWord([1, 1, 0, 0]).to_triangulation()
            [(0, 2)]

            sage: [t.to_triangulation() for t in DyckWords(3)]
            [[(2, 4), (1, 4)],
            [(2, 4), (0, 2)],
            [(1, 3), (1, 4)],
            [(1, 3), (0, 3)],
            [(0, 2), (0, 3)]]

        REFERENCES:

        - [Cha2005]_
        """
        n = self.number_of_open_symbols()
        l = list(range(n + 2))  # from 0 to n + 1
        edges = []
        coheight = n - 1
        for letter in self[1:-1]:
            if letter == 1:
                coheight -= 1
            else:
                edges.append((coheight, l[coheight + 2]))
                l.pop(coheight + 1)

        return edges

    def to_triangulation_as_graph(self):
        r"""
        Map ``self`` to a triangulation and return the result as a graph.

        See :meth:`to_triangulation` for the bijection used to map
        complete Dyck words to triangulations.

        OUTPUT:

        - a graph containing both the perimeter edges and the inner
          edges of a triangulation of a regular polygon.

        EXAMPLES::

            sage: g = DyckWord([1, 1, 0, 0, 1, 0]).to_triangulation_as_graph(); g       # needs sage.graphs
            Graph on 5 vertices
            sage: g.edges(sort=True, labels=False)                                      # needs sage.graphs
            [(0, 1), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (3, 4)]
            sage: g.show()                      # not tested                            # needs sage.graphs
        """
        n = self.number_of_open_symbols()
        edges = self.to_triangulation()
        from sage.graphs.graph import Graph
        peri = [(i, i + 1) for i in range(n + 1)] + [(n + 1, 0)]
        g = Graph(n + 2)
        g.add_edges(peri)
        g.add_edges(edges)
        g.set_pos(g.layout_circular())
        return g

    def to_non_decreasing_parking_function(self):
        r"""
        Bijection to :class:`non-decreasing parking
        functions<sage.combinat.non_decreasing_parking_function.NonDecreasingParkingFunctions>`.

        See there the method
        :meth:`~sage.combinat.non_decreasing_parking_function.NonDecreasingParkingFunction.to_dyck_word`
        for more information.

        EXAMPLES::

            sage: DyckWord([]).to_non_decreasing_parking_function()
            []
            sage: DyckWord([1,0]).to_non_decreasing_parking_function()
            [1]
            sage: DyckWord([1,1,0,0]).to_non_decreasing_parking_function()
            [1, 1]
            sage: DyckWord([1,0,1,0]).to_non_decreasing_parking_function()
            [1, 2]
            sage: DyckWord([1,0,1,1,0,1,0,0,1,0]).to_non_decreasing_parking_function()
            [1, 2, 2, 3, 5]

        TESTS::

            sage: ld = DyckWords(5)
            sage: list(ld) == [dw.to_non_decreasing_parking_function().to_dyck_word() for dw in ld]
            True
        """
        from sage.combinat.non_decreasing_parking_function import NonDecreasingParkingFunction
        return NonDecreasingParkingFunction.from_dyck_word(self)

    def major_index(self) -> int:
        r"""
        Return the major index of ``self`` .

        The major index of a Dyck word `D` is the sum of the positions of
        the valleys of `D` (when started counting at position ``1``).

        EXAMPLES::

            sage: DyckWord([1, 0, 1, 0]).major_index()
            2
            sage: DyckWord([1, 1, 0, 0]).major_index()
            0
            sage: DyckWord([1, 1, 0, 0, 1, 0]).major_index()
            4
            sage: DyckWord([1, 0, 1, 1, 0, 0]).major_index()
            2

        TESTS::

            sage: DyckWord([]).major_index()
            0
            sage: DyckWord([1, 0]).major_index()
            0
        """
        valleys = self.valleys()
        return sum(valleys) + len(valleys)

    def pyramid_weight(self) -> int:
        r"""
        Return the pyramid weight of ``self``.

        A pyramid of ``self`` is a subsequence of the form
        `1^h 0^h`. A pyramid is maximal if it is neither preceded by a `1`
        nor followed by a `0`.

        The pyramid weight of a Dyck path is the sum of the lengths of the
        maximal pyramids and was defined in [DS1992]_.

        EXAMPLES::

            sage: DyckWord([1,1,0,1,1,1,0,0,1,0,0,0,1,1,0,0]).pyramid_weight()
            6
            sage: DyckWord([1,1,1,0,0,0]).pyramid_weight()
            3
            sage: DyckWord([1,0,1,0,1,0]).pyramid_weight()
            3
            sage: DyckWord([1,1,0,1,0,0]).pyramid_weight()
            2
        """
        aseq = self.to_area_sequence() + [0]
        bseq = self.reverse().to_area_sequence() + [0]
        apeak = []
        bpeak = []
        for i in range(len(aseq) - 1):
            if aseq[i + 1] <= aseq[i]:
                apeak.append(i)
            if bseq[i + 1] <= bseq[i]:
                bpeak.append(i)
        out = 0
        for i, apeaki in enumerate(apeak):
            out += min(aseq[apeaki] - aseq[apeaki + 1] + 1,
                       bseq[bpeak[-i - 1]] - bseq[bpeak[-i - 1] + 1] + 1)
        return out

    def tunnels(self):
        r"""
        Return an iterator of ranges of the matching parentheses in the Dyck
        word ``self``.

        That is, if ``(a,b)`` is in ``self.tunnels()``, then the matching
        parenthesis to ``self[a]`` is ``self[b-1]``.

        EXAMPLES::

            sage: list(DyckWord([1, 1, 0, 1, 1, 0, 0, 1, 0, 0]).tunnels())
            [(0, 10), (1, 3), (3, 7), (4, 6), (7, 9)]
        """
        heights = self.heights()
        for i in range(len(heights) - 1):
            height = heights[i]
            if height < heights[i + 1]:
                yield (i, i + 1 + heights[i + 1:].index(height))

    def number_of_tunnels(self, tunnel_type='centered') -> int:
        r"""
        Return the number of tunnels of ``self`` of type ``tunnel_type``.

        A tunnel is a pair `(a,b)` where ``a`` is the position of an open
        parenthesis and ``b`` is the position of the matching close
        parenthesis. If `a + b = n` then the tunnel is called *centered* . If
        `a + b < n` then the tunnel is called *left* and if `a + b > n`, then
        the tunnel is called *right*.

        INPUT:

        - ``tunnel_type`` -- (default: ``'centered'``) can be one of the
          following: ``'left'``, ``'right'``, ``'centered'``, or ``'all'``

        EXAMPLES::

            sage: DyckWord([1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0]).number_of_tunnels()
            0
            sage: DyckWord([1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0]).number_of_tunnels('left')
            5
            sage: DyckWord([1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0]).number_of_tunnels('right')
            2
            sage: DyckWord([1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0]).number_of_tunnels('all')
            7
            sage: DyckWord([1, 1, 0, 0]).number_of_tunnels('centered')
            2
        """
        n = len(self)
        tunnels = self.tunnels()
        if tunnel_type == 'left':
            return len([1 for (i, j) in tunnels if i + j < n])
        elif tunnel_type == 'centered':
            return len([1 for (i, j) in tunnels if i + j == n])
        elif tunnel_type == 'right':
            return len([1 for (i, j) in tunnels if i + j > n])
        elif tunnel_type == 'all':
            return len(list(tunnels))
        else:
            raise ValueError("the given tunnel_type is not valid")

    @combinatorial_map(order=2, name="Reverse path")
    def reverse(self) -> DyckWord:
        r"""
        Return the reverse and complement of ``self``.

        This operation corresponds to flipping the Dyck path across the
        `y=-x` line.

        EXAMPLES::

            sage: DyckWord([1,1,0,0,1,0]).reverse()
            [1, 0, 1, 1, 0, 0]
            sage: DyckWord([1,1,1,0,0,0]).reverse()
            [1, 1, 1, 0, 0, 0]
            sage: len([D for D in DyckWords(5) if D.reverse() == D])
            10

        TESTS::

            sage: DyckWord([]).reverse()
            []
        """
        alist = []
        for i in range(len(self)):
            if self[i] == open_symbol:
                alist.append(close_symbol)
            else:
                alist.append(open_symbol)
        alist.reverse()
        return DyckWord(alist)  # type:ignore

    def first_return_decomposition(self) -> tuple:
        r"""
        Decompose a Dyck word into a pair of Dyck words (potentially empty)
        where the first word consists of the word after the first up step and
        the corresponding matching closing parenthesis.

        EXAMPLES::

            sage: DyckWord([1,1,0,1,0,0,1,1,0,1,0,1,0,0]).first_return_decomposition()
            ([1, 0, 1, 0], [1, 1, 0, 1, 0, 1, 0, 0])
            sage: DyckWord([1,1,0,0]).first_return_decomposition()
            ([1, 0], [])
            sage: DyckWord([1,0,1,0]).first_return_decomposition()
            ([], [1, 0])
        """
        k = self.position_of_first_return() * 2
        return DyckWord(self[1:k - 1]), DyckWord(self[k:])  # type:ignore

    def decomposition_reverse(self) -> DyckWord:
        r"""
        Return the involution of ``self`` with a recursive definition.

        If a Dyck word `D` decomposes as `1 D_1 0 D_2` where `D_1` and
        `D_2` are complete Dyck words then the decomposition reverse is
        `1 \phi(D_2) 0 \phi(D_1)`.

        EXAMPLES::

            sage: DyckWord([1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0]).decomposition_reverse()
            [1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0]
            sage: DyckWord([1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0]).decomposition_reverse()
            [1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0]
            sage: DyckWord([1,1,0,0]).decomposition_reverse()
            [1, 0, 1, 0]
            sage: DyckWord([1,0,1,0]).decomposition_reverse()
            [1, 1, 0, 0]
        """
        if not self:
            return self
        D1, D2 = self.first_return_decomposition()
        D = [1] + list(D2.decomposition_reverse())
        D += [0] + list(D1.decomposition_reverse())
        return DyckWord(D)  # type:ignore

    @combinatorial_map(name="Area-dinv to bounce-area")
    def area_dinv_to_bounce_area_map(self) -> DyckWord:
        r"""
        Return the image of ``self`` under the map which sends a
        Dyck word with ``area`` equal to `r` and ``dinv`` equal to `s` to a
        Dyck word with ``bounce`` equal to `r` and ``area`` equal to `s` .

        The inverse of this map is :meth:`bounce_area_to_area_dinv_map`.

        For a definition of this map, see [Hag2008]_ p. 50 where it is called
        `\zeta`. However, this map differs from Haglund's map by an application
        of :meth:`reverse` (as does the definition of the :meth:`bounce`
        statistic).

        EXAMPLES::

            sage: DyckWord([1,1,0,1,0,0,1,1,0,1,0,1,0,0]).area_dinv_to_bounce_area_map()
            [1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0]
            sage: DyckWord([1,1,0,1,0,0,1,1,0,1,0,1,0,0]).area()
            5
            sage: DyckWord([1,1,0,1,0,0,1,1,0,1,0,1,0,0]).dinv()
            13
            sage: DyckWord([1,1,1,1,1,0,0,0,1,0,0,1,0,0]).area()
            13
            sage: DyckWord([1,1,1,1,1,0,0,0,1,0,0,1,0,0]).bounce()
            5
            sage: DyckWord([1,1,1,1,1,0,0,0,1,0,0,1,0,0]).area_dinv_to_bounce_area_map()
            [1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0]
            sage: DyckWord([1,1,0,0]).area_dinv_to_bounce_area_map()
            [1, 0, 1, 0]
            sage: DyckWord([1,0,1,0]).area_dinv_to_bounce_area_map()
            [1, 1, 0, 0]
        """
        if not self:
            return self
        a = self.to_area_sequence()
        a.reverse()
        image = []
        for i in range(max(a), -2, -1):
            for j in a:
                if j == i:
                    image.append(1)
                elif j == i + 1:
                    image.append(0)
        return DyckWord(image)  # type:ignore

    @combinatorial_map(name="Bounce-area to area-dinv")
    def bounce_area_to_area_dinv_map(self) -> DyckWord:
        r"""
        Return the image of the Dyck word under the map which sends a
        Dyck word with ``bounce`` equal to `r` and ``area`` equal to `s` to a
        Dyck word with ``area`` equal to `r` and ``dinv`` equal to `s` .

        This implementation uses a recursive method by saying that the
        last entry in the area sequence of the Dyck word ``self`` is
        equal to the number of touch points of the Dyck path minus 1
        of the image of this map.

        The inverse of this map is :meth:`area_dinv_to_bounce_area_map`.

        For a definition of this map, see [Hag2008]_ p. 50 where it is called
        `\zeta^{-1}`. However, this map differs from Haglund's map by an
        application of :meth:`reverse` (as does the definition of the
        :meth:`bounce` statistic).

        EXAMPLES::

            sage: DyckWord([1,1,0,1,0,0,1,1,0,1,0,1,0,0]).bounce_area_to_area_dinv_map()
            [1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0]
            sage: DyckWord([1,1,0,1,0,0,1,1,0,1,0,1,0,0]).area()
            5
            sage: DyckWord([1,1,0,1,0,0,1,1,0,1,0,1,0,0]).bounce()
            9
            sage: DyckWord([1,1,0,0,1,1,1,1,0,0,1,0,0,0]).area()
            9
            sage: DyckWord([1,1,0,0,1,1,1,1,0,0,1,0,0,0]).dinv()
            5
            sage: all(D==D.bounce_area_to_area_dinv_map().area_dinv_to_bounce_area_map() for D in DyckWords(6))
            True
            sage: DyckWord([1,1,0,0]).bounce_area_to_area_dinv_map()
            [1, 0, 1, 0]
            sage: DyckWord([1,0,1,0]).bounce_area_to_area_dinv_map()
            [1, 1, 0, 0]
        """
        aseq = self._area_sequence_iter()
        out: list[int] = []
        zeros: list[int] = []
        for ai in aseq:
            p = (zeros + [len(out)])[ai]
            out = [1] + out[p:] + [0] + out[:p]
            zeros = [0] + [j + len(out) - p for j in zeros[:ai]]
        return DyckWord(out)  # type:ignore

    def area(self) -> int:
        r"""
        Return the area for ``self`` corresponding to the area
        of the Dyck path.

        One can view a balanced Dyck word as a lattice path from
        `(0,0)` to `(n,n)` in the first quadrant by letting
        '1's represent steps in the direction `(1,0)` and '0's
        represent steps in the direction `(0,1)`. The resulting
        path will remain weakly above the diagonal `y = x`.

        The area statistic is the number of complete
        squares in the integer lattice which are below the path and above
        the line `y = x`. The 'half-squares' directly above the
        line `y = x` do not contribute to this statistic.

        EXAMPLES::

            sage: dw = DyckWord([1,0,1,0])
            sage: dw.area() # 2 half-squares, 0 complete squares
            0

        ::

            sage: dw = DyckWord([1,1,1,0,1,1,1,0,0,0,1,1,0,0,1,0,0,0])
            sage: dw.area()
            19

        ::

            sage: DyckWord([1,1,1,1,0,0,0,0]).area()
            6
            sage: DyckWord([1,1,1,0,1,0,0,0]).area()
            5
            sage: DyckWord([1,1,1,0,0,1,0,0]).area()
            4
            sage: DyckWord([1,1,1,0,0,0,1,0]).area()
            3
            sage: DyckWord([1,0,1,1,0,1,0,0]).area()
            2
            sage: DyckWord([1,1,0,1,1,0,0,0]).area()
            4
            sage: DyckWord([1,1,0,0,1,1,0,0]).area()
            2
            sage: DyckWord([1,0,1,1,1,0,0,0]).area()
            3
            sage: DyckWord([1,0,1,1,0,0,1,0]).area()
            1
            sage: DyckWord([1,0,1,0,1,1,0,0]).area()
            1
            sage: DyckWord([1,1,0,0,1,0,1,0]).area()
            1
            sage: DyckWord([1,1,0,1,0,0,1,0]).area()
            2
            sage: DyckWord([1,1,0,1,0,1,0,0]).area()
            3
            sage: DyckWord([1,0,1,0,1,0,1,0]).area()
            0
        """
        return sum(self._area_sequence_iter())

    def bounce_path(self) -> DyckWord:
        r"""
        Return the bounce path of ``self`` formed by starting at `(n,n)` and
        traveling West until encountering the first vertical step of ``self``,
        then South until encountering the diagonal, then West again to hit
        the path, etc. until the `(0,0)` point is reached.  The path followed
        by this walk is the bounce path.

        .. SEEALSO:: :meth:`bounce`

        EXAMPLES::

            sage: DyckWord([1,1,0,1,0,0]).bounce_path()
            [1, 0, 1, 1, 0, 0]
            sage: DyckWord([1,1,1,0,0,0]).bounce_path()
            [1, 1, 1, 0, 0, 0]
            sage: DyckWord([1,0,1,0,1,0]).bounce_path()
            [1, 0, 1, 0, 1, 0]
            sage: DyckWord([1,1,1,1,0,0,1,0,0,0]).bounce_path()
            [1, 1, 0, 0, 1, 1, 1, 0, 0, 0]

        TESTS::

            sage: DyckWord([]).bounce_path()
            []
            sage: DyckWord([1,0]).bounce_path()
            [1, 0]
        """
        area_seq = self.to_area_sequence()
        i = len(area_seq) - 1
        n = 5
        while i > 0:
            n -= 1
            a = area_seq[i]
            i_new = i - a
            while i > i_new:
                i -= 1
                area_seq[i] = area_seq[i + 1] - 1
            i -= 1
        return DyckWord(area_sequence=area_seq)  # type:ignore

    def bounce(self) -> int:
        r"""
        Return the bounce statistic of ``self`` due to J. Haglund,
        see [Hag2008]_.

        One can view a balanced Dyck word as a lattice path from `(0,0)` to
        `(n,n)` in the first quadrant by letting '1's represent steps in
        the direction `(0,1)` and '0's represent steps in the direction
        `(1,0)`.  The resulting path will remain weakly above the diagonal
        `y = x`.

        We describe the bounce statistic of such a path in terms of what is
        known as the "bounce path".

        We can think of our bounce path as describing the trail of a billiard
        ball shot West from `(n, n)`, which "bounces" down whenever it
        encounters a vertical step and "bounces" left when it encounters the
        line `y = x`.

        The bouncing ball will strike the diagonal at the places

        .. MATH::

            (0, 0), (j_1, j_1), (j_2, j_2), \ldots, (j_r-1, j_r-1), (j_r, j_r)
            = (n, n).

        We define the bounce to be the sum `\sum_{i=1}^{r-1} j_i`.

        EXAMPLES::

            sage: DyckWord([1,1,1,0,1,1,1,0,0,0,1,1,0,0,1,0,0,0]).bounce()
            7
            sage: DyckWord([1,1,1,1,0,0,0,0]).bounce()
            0
            sage: DyckWord([1,1,1,0,1,0,0,0]).bounce()
            1
            sage: DyckWord([1,1,1,0,0,1,0,0]).bounce()
            2
            sage: DyckWord([1,1,1,0,0,0,1,0]).bounce()
            3
            sage: DyckWord([1,0,1,1,0,1,0,0]).bounce()
            3
            sage: DyckWord([1,1,0,1,1,0,0,0]).bounce()
            1
            sage: DyckWord([1,1,0,0,1,1,0,0]).bounce()
            2
            sage: DyckWord([1,0,1,1,1,0,0,0]).bounce()
            1
            sage: DyckWord([1,0,1,1,0,0,1,0]).bounce()
            4
            sage: DyckWord([1,0,1,0,1,1,0,0]).bounce()
            3
            sage: DyckWord([1,1,0,0,1,0,1,0]).bounce()
            5
            sage: DyckWord([1,1,0,1,0,0,1,0]).bounce()
            4
            sage: DyckWord([1,1,0,1,0,1,0,0]).bounce()
            2
            sage: DyckWord([1,0,1,0,1,0,1,0]).bounce()
            6
        """
        x_pos = len(self) // 2
        y_pos = len(self) // 2

        b = 0

        mode = "left"
        makeup_steps = 0
        l = self._list[:]
        l.reverse()

        for move in l:
            if mode == "left":
                if move == close_symbol:
                    x_pos -= 1
                elif move == open_symbol:
                    y_pos -= 1
                    if x_pos == y_pos:
                        b += x_pos
                    else:
                        mode = "drop"
            elif mode == "drop":
                if move == close_symbol:
                    makeup_steps += 1
                elif move == open_symbol:
                    y_pos -= 1
                    if x_pos == y_pos:
                        b += x_pos
                        mode = "left"
                        x_pos -= makeup_steps
                        makeup_steps = 0

        return b

    def dinv(self, labeling=None) -> int:
        r"""
        Return the dinv statistic of ``self`` due to M. Haiman, see [Hag2008]_.

        If a labeling is provided then this function returns the dinv of the
        labeled Dyck word.

        INPUT:

        - ``labeling`` -- an optional argument to be viewed as the labelings
          of the vertical edges of the Dyck path

        OUTPUT:

        - an integer representing the ``dinv`` statistic of the Dyck path
          or the labelled Dyck path.

        EXAMPLES::

            sage: DyckWord([1,0,1,0,1,0,1,0,1,0]).dinv()
            10
            sage: DyckWord([1,1,1,1,1,0,0,0,0,0]).dinv()
            0
            sage: DyckWord([1,1,1,1,0,1,0,0,0,0]).dinv()
            1
            sage: DyckWord([1,1,0,1,0,0,1,1,0,1,0,1,0,0]).dinv()
            13
            sage: DyckWord([1,1,0,1,0,0,1,1,0,1,0,1,0,0]).dinv([1,2,3,4,5,6,7])
            11
            sage: DyckWord([1,1,0,1,0,0,1,1,0,1,0,1,0,0]).dinv([6,7,5,3,4,2,1])
            2
        """
        alist = self.to_area_sequence()
        cnt = 0
        for j, aj in enumerate(alist):
            if labeling is not None:
                lj = labeling[j]
            for i in range(j):
                if (alist[i] == aj and (labeling is None or labeling[i] < lj)) or (alist[i] - aj == 1 and (labeling is None or labeling[i] > lj)):
                    cnt += 1
        return cnt

    @combinatorial_map(name='to alternating sign matrix')
    def to_alternating_sign_matrix(self):
        r"""
        Return ``self`` as an alternating sign matrix.

        This is an inclusion map from Dyck words of length `2n` to certain
        `n \times n` alternating sign matrices.

        EXAMPLES::

            sage: DyckWord([1,1,1,0,1,0,0,0]).to_alternating_sign_matrix()              # needs sage.modules
            [ 0  0  1  0]
            [ 1  0 -1  1]
            [ 0  1  0  0]
            [ 0  0  1  0]
            sage: DyckWord([1,0,1,0,1,1,0,0]).to_alternating_sign_matrix()              # needs sage.modules
            [1 0 0 0]
            [0 1 0 0]
            [0 0 0 1]
            [0 0 1 0]
        """
        parkfn = self.reverse().to_non_decreasing_parking_function()
        parkfn2 = [len(parkfn) + 1 - parkfn[i] for i in range(len(parkfn))]
        monotone_triangle = [[0] * (len(parkfn2) - j)
                             for j in range(len(parkfn2))]
        for i in range(len(monotone_triangle)):
            for j in range(len(monotone_triangle[i])):
                monotone_triangle[i][j] = len(monotone_triangle[i]) - j
            monotone_triangle[i][0] = parkfn2[i]
        A = AlternatingSignMatrices(len(parkfn))
        return A.from_monotone_triangle(monotone_triangle)


class DyckWords(UniqueRepresentation, Parent):
    r"""
    Dyck words.

    A Dyck word is a sequence `(w_1, \ldots, w_n)` consisting of 1 s and 0 s,
    with the property that for any `i` with `1 \leq i \leq n`, the sequence
    `(w_1, \ldots, w_i)` contains at least as many 1 s as 0 s.

    A Dyck word is balanced if the total number of 1 s is equal to the total
    number of 0 s. The number of balanced Dyck words of length `2k` is given
    by the :func:`Catalan number<sage.combinat.combinat.catalan_number>` `C_k`.

    EXAMPLES:

    This class can be called with three keyword parameters ``k1``, ``k2``
    and ``complete``.

    If neither ``k1`` nor ``k2`` are specified, then :class:`DyckWords`
    returns the combinatorial class of all balanced (=complete) Dyck words,
    unless the keyword ``complete`` is set to False (in which case it
    returns the class of all Dyck words).

    ::

        sage: DW = DyckWords(); DW
        Complete Dyck words
        sage: [] in DW
        True
        sage: [1, 0, 1, 0] in DW
        True
        sage: [1, 1, 0] in DW
        False
        sage: ADW = DyckWords(complete=False); ADW
        Dyck words
        sage: [] in ADW
        True
        sage: [1, 0, 1, 0] in ADW
        True
        sage: [1, 1, 0] in ADW
        True
        sage: [1, 0, 0] in ADW
        False

    If just ``k1`` is specified, then it returns the balanced Dyck words with
    ``k1`` opening parentheses and ``k1`` closing parentheses.

    ::

        sage: DW2 = DyckWords(2); DW2
        Dyck words with 2 opening parentheses and 2 closing parentheses
        sage: DW2.first()
        [1, 0, 1, 0]
        sage: DW2.last()
        [1, 1, 0, 0]
        sage: DW2.cardinality()
        2
        sage: DyckWords(100).cardinality() == catalan_number(100)
        True

    If ``k2`` is specified in addition to ``k1``, then it returns the
    Dyck words with ``k1`` opening parentheses and ``k2`` closing parentheses.

    ::

        sage: DW32 = DyckWords(3,2); DW32
        Dyck words with 3 opening parentheses and 2 closing parentheses
        sage: DW32.list()
        [[1, 0, 1, 0, 1],
         [1, 0, 1, 1, 0],
         [1, 1, 0, 0, 1],
         [1, 1, 0, 1, 0],
         [1, 1, 1, 0, 0]]
    """
    @staticmethod
    def __classcall_private__(cls, k1=None, k2=None, complete=True):
        """
        Choose the correct parent based upon input.

        EXAMPLES::

            sage: DW1 = DyckWords(3,3)
            sage: DW2 = DyckWords(3)
            sage: DW1 is DW2
            True
        """
        if k2 is None:
            if k1 is None:
                if complete:
                    return CompleteDyckWords_all()
                return DyckWords_all()

            k1 = Integer(k1)
            if k1 < 0:
                raise ValueError("k1 (= %s) must be nonnegative" % k1)
            return CompleteDyckWords_size(k1)
        else:
            k1 = Integer(k1)

        k2 = Integer(k2)
        if k1 < 0 or (k2 is not None and k2 < 0):
            raise ValueError("k1 (= %s) and k2 (= %s) must be nonnegative, with k1 >= k2" % (k1, k2))
        if k1 < k2:
            raise ValueError("k1 (= %s) must be >= k2 (= %s)" % (k1, k2))

        if k1 == k2:
            return CompleteDyckWords_size(k1)
        return DyckWords_size(k1, k2)

    Element = DyckWord

    # add options to class
    class options(GlobalOptions):
        r"""
        Set and display the options for Dyck words. If no parameters
        are set, then the function returns a copy of the options dictionary.

        The ``options`` to Dyck words can be accessed as the method
        :meth:`DyckWords.options` of :class:`DyckWords` and
        related parent classes.

        @OPTIONS

        EXAMPLES::

            sage: D = DyckWord([1, 1, 0, 1, 0, 0])
            sage: D
            [1, 1, 0, 1, 0, 0]
            sage: DyckWords.options.display="lattice"
            sage: D
               ___
             _| x
            | x  .
            |  . .
            sage: DyckWords.options(diagram_style='line')
            sage: D
             /\/\
            /    \
            sage: DyckWords.options._reset()
        """
        NAME = 'DyckWords'
        module = 'sage.combinat.dyck_word'
        display = dict(default='list',
            description='Specifies how Dyck words should be printed',
            values=dict(list='displayed as a list',
                        lattice='displayed on the lattice defined by ``diagram_style``'),
            case_sensitive=False)
        ascii_art = dict(default='path',
            description='Specifies how the ascii art of Dyck words should be printed',
            values=dict(path="Using the path string",
                        pretty_output="Using pretty printing"),
            alias=dict(pretty_print='pretty_output', path_string='path'),
            case_sensitive=False)
        diagram_style = dict(default='grid',
            values=dict(grid='printing as paths on a grid using N and E steps',
                        line='printing as paths on a line using NE and SE steps',),
            alias={'N-E': 'grid', 'NE-SE': 'line'},
            case_sensitive=False)
        latex_tikz_scale = dict(default=1,
            description='The default value for the tikz scale when latexed',
            checker=lambda x: True)  # More trouble than it's worth to check
        latex_diagonal = dict(default=False,
            description='The default value for displaying the diagonal when latexed',
            checker=lambda x: isinstance(x, bool))
        latex_line_width_scalar = dict(default=2,
            description='The default value for the line width as a '
                        'multiple of the tikz scale when latexed',
            checker=lambda x: True)  # More trouble than it's worth to check
        latex_color = dict(default='black',
            description='The default value for the color when latexed',
            checker=lambda x: isinstance(x, str))
        latex_bounce_path = dict(default=False,
            description='The default value for displaying the bounce path when latexed',
            checker=lambda x: isinstance(x, bool))
        latex_peaks = dict(default=False,
            description='The default value for displaying the peaks when latexed',
            checker=lambda x: isinstance(x, bool))
        latex_valleys = dict(default=False,
            description='The default value for displaying the valleys when latexed',
            checker=lambda x: isinstance(x, bool))

    def _element_constructor_(self, word):
        """
        Construct an element of ``self``.

        EXAMPLES::

            sage: D = DyckWords()
            sage: elt = D([1, 1, 0, 1, 0, 0]); elt
            [1, 1, 0, 1, 0, 0]
            sage: elt.parent() is D
            True
        """
        if isinstance(word, DyckWord) and word.parent() is self:
            return word
        return self.element_class(self, list(word))

    def __contains__(self, x) -> bool:
        r"""
        TESTS::

            sage: D = DyckWords(complete=False)
            sage: [] in D
            True
            sage: [1] in D
            True
            sage: [0] in D
            False
            sage: [1, 0] in D
            True
        """
        if isinstance(x, DyckWord):
            return True

        if not isinstance(x, list):
            return False

        return is_a(x)

    def from_heights(self, heights) -> DyckWord:
        r"""
        Compute a Dyck word knowing its heights.

        We view the Dyck word as a Dyck path from `(0, 0)` to
        `(2n, 0)` in the first quadrant by letting ``1``'s represent
        steps in the direction `(1, 1)` and ``0``'s represent steps in
        the direction `(1, -1)`.

        The :meth:`~DyckWord.heights` is the sequence of the `y`-coordinates of
        the `2n+1` lattice points along this path.

        EXAMPLES::

            sage: from sage.combinat.dyck_word import DyckWord
            sage: D = DyckWords(complete=False)
            sage: D.from_heights((0,))
            []
            sage: D.from_heights((0, 1, 0))
            [1, 0]
            sage: D.from_heights((0, 1, 2, 1, 0))
            [1, 1, 0, 0]

        This also works for incomplete Dyck words::

            sage: D.from_heights((0, 1, 2, 1, 2, 1))
            [1, 1, 0, 1, 0]
            sage: D.from_heights((0, 1, 2, 1))
            [1, 1, 0]

        .. SEEALSO:: :meth:`~DyckWord.heights`, :meth:`min_from_heights`

        TESTS::

            sage: all(dw == D.from_heights(dw.heights())
            ....:     for i in range(7) for dw in DyckWords(i))
            True

            sage: D.from_heights((1, 2, 1))
            Traceback (most recent call last):
            ...
            ValueError: heights must start with 0: (1, 2, 1)
            sage: D.from_heights((0, 1, 4, 1))
            Traceback (most recent call last):
            ...
            ValueError: consecutive heights must differ by exactly 1: (0, 1, 4, 1)
            sage: D.from_heights(())
            Traceback (most recent call last):
            ...
            ValueError: heights must start with 0: ()
        """
        l1 = len(heights) - 1
        res = [0] * (l1)
        if not heights or heights[0] != 0:
            raise ValueError("heights must start with 0: %s" % (heights,))
        for i in range(l1):
            if heights[i] == heights[i + 1] - 1:
                res[i] = 1
            elif heights[i] != heights[i + 1] + 1:
                raise ValueError("consecutive heights must differ by exactly 1: %s" % (heights,))
        return self.element_class(self, res)

    def min_from_heights(self, heights) -> DyckWord:
        r"""
        Compute the smallest Dyck word which achieves or surpasses
        a given sequence of heights.

        INPUT:

        - ``heights`` -- a nonempty list or iterable consisting of
          nonnegative integers, the first of which is `0`

        OUTPUT:

        - The smallest Dyck word whose sequence of heights is
          componentwise higher-or-equal to ``heights``. Here,
          "smaller" can be understood both in the sense of
          lexicographic order on the Dyck words, and in the sense
          of every vertex of the path having the smallest possible
          height.

        .. SEEALSO::

            - :meth:`~DyckWord.heights`
            - :meth:`from_heights`

        EXAMPLES::

            sage: D = DyckWords(complete=False)
            sage: D.min_from_heights((0,))
            []
            sage: D.min_from_heights((0, 1, 0))
            [1, 0]
            sage: D.min_from_heights((0, 0, 2, 0, 0))
            [1, 1, 0, 0]
            sage: D.min_from_heights((0, 0, 2, 0, 2, 0))
            [1, 1, 0, 1, 0]
            sage: D.min_from_heights((0, 0, 1, 0, 1, 0))
            [1, 1, 0, 1, 0]

        TESTS::

            sage: D.min_from_heights(())
            Traceback (most recent call last):
            ...
            ValueError: heights must start with 0: ()
        """
        if not heights or heights[0] != 0:
            raise ValueError("heights must start with 0: %s" % (heights,))
        # round heights to the smallest even-odd integer
        heights = list(heights)
        for i in range(0, len(heights), 2):
            if heights[i] % 2:
                heights[i] += 1
        for i in range(1, len(heights), 2):
            if heights[i] % 2 == 0:
                heights[i] += 1

        # smooth heights
        for i in range(len(heights) - 1):
            if heights[i + 1] < heights[i]:
                heights[i + 1] = heights[i] - 1
        for i in range(len(heights) - 1, 0, -1):
            if heights[i] > heights[i - 1]:
                heights[i - 1] = heights[i] - 1
        return self.from_heights(heights)


class DyckWords_all(DyckWords):
    """
    All Dyck words.
    """

    def __init__(self):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: TestSuite(DyckWords(complete=False)).run()
        """
        DyckWords.__init__(self, category=InfiniteEnumeratedSets())

    def _repr_(self) -> str:
        r"""
        TESTS::

            sage: DyckWords(complete=False)
            Dyck words
        """
        return "Dyck words"

    def _an_element_(self) -> DyckWord:
        r"""
        TESTS::

            sage: DyckWords(complete=False).an_element()
            [1, 0, 1, 0, 1, 0, 1, 0, 1, 0]
        """
        return DyckWords(5).an_element()

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: it = DyckWords(complete=False).__iter__()
            sage: [next(it) for x in range(10)]
            [[],
             [1],
             [1, 0],
             [1, 1],
             [1, 0, 0],
             [1, 0, 1],
             [1, 1, 0],
             [1, 1, 1],
             [1, 0, 1, 0],
             [1, 1, 0, 0]]
        """
        n = 0
        yield self.element_class(self, [])
        while True:
            for k in range(1, n + 1):
                for x in DyckWords_size(k, n - k):
                    yield self.element_class(self, list(x))
            n += 1


class DyckWordBacktracker(GenericBacktracker):
    r"""
    This class is an iterator for all Dyck words
    with `n` opening parentheses and `n - k` closing parentheses using
    the backtracker class. It is used by the :class:`DyckWords_size` class.

    This is not really meant to be called directly, partially because
    it fails in a couple corner cases: ``DWB(0)`` yields ``[0]``, not the
    empty word, and ``DWB(k, k+1)`` yields something (it shouldn't yield
    anything). This could be fixed with a sanity check in ``_rec()``, but
    then we'd be doing the sanity check *every time* we generate new
    objects; instead, we do *one* sanity check in :class:`DyckWords` and
    assume here that the sanity check has already been made.

    AUTHOR:

    - Dan Drake (2008-05-30)
    """

    def __init__(self, k1, k2):
        r"""
        TESTS::

            sage: from sage.combinat.dyck_word import DyckWordBacktracker
            sage: len(list(DyckWordBacktracker(5, 5)))
            42
            sage: len(list(DyckWordBacktracker(6,4)))
            90
            sage: len(list(DyckWordBacktracker(7,0)))
            1
        """
        GenericBacktracker.__init__(self, [], (0, 0))
        # note that the comments in this class think of our objects as
        # Dyck paths, not words; having k1 opening parens and k2 closing
        # parens corresponds to paths of length k1 + k2 ending at height
        # k1 - k2.
        k1 = Integer(k1)
        k2 = Integer(k2)
        self.n = k1 + k2
        self.endht = k1 - k2

    def _rec(self, path, state):
        r"""
        TESTS::

            sage: from sage.combinat.dyck_word import DyckWordBacktracker
            sage: dwb = DyckWordBacktracker(3, 3)
            sage: list(dwb._rec([1,1,0],(3, 2)))
            [([1, 1, 0, 0], (4, 1), False), ([1, 1, 0, 1], (4, 3), False)]
            sage: list(dwb._rec([1,1,0,0],(4, 0)))
            [([1, 1, 0, 0, 1], (5, 1), False)]
            sage: list(DyckWordBacktracker(4, 4)._rec([1,1,1,1],(4, 4)))
            [([1, 1, 1, 1, 0], (5, 3), False)]
        """
        len, ht = state

        if len < self.n - 1:
            # if length is less than n-1, new path won't have length n, so
            # don't yield it, and keep building paths

            # if the path isn't too low and is not touching the x-axis, we can
            # yield a path with a downstep at the end
            if ht > (self.endht - (self.n - len)) and ht > 0:
                yield path + [0], (len + 1, ht - 1), False

            # if the path isn't too high, it can also take an upstep
            if ht < (self.endht + (self.n - len)):
                yield path + [1], (len + 1, ht + 1), False
        else:
            # length is n - 1, so add a single step (up or down,
            # according to current height and endht), don't try to
            # construct more paths, and yield the path
            if ht < self.endht:
                yield path + [1], None, True
            else:
                yield path + [0], None, True


class DyckWords_size(DyckWords):
    """
    Dyck words with `k_1` openers and `k_2` closers.
    """

    def __init__(self, k1, k2):
        r"""
        TESTS:

        Check that :issue:`18244` is fixed::

            sage: DyckWords(13r, 8r).cardinality()
            87210
            sage: parent(_)
            Integer Ring
            sage: TestSuite(DyckWords(4,2)).run()
        """
        self.k1 = Integer(k1)
        self.k2 = Integer(k2)
        DyckWords.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self) -> str:
        r"""
        TESTS::

            sage: DyckWords(4)
            Dyck words with 4 opening parentheses and 4 closing parentheses
        """
        return "Dyck words with %s opening parentheses and %s closing parentheses" % (self.k1, self.k2)

    def __contains__(self, x) -> bool:
        r"""
        EXAMPLES::

            sage: [1, 0, 0, 1] in DyckWords(2,2)
            False
            sage: [1, 0, 1, 0] in DyckWords(2,2)
            True
            sage: [1, 0, 1, 0, 1] in DyckWords(3,2)
            True
            sage: [1, 0, 1, 1, 0] in DyckWords(3,2)
            True
            sage: [1, 0, 1, 1] in DyckWords(3,1)
            True
        """
        return is_a(x, self.k1, self.k2)

    def __iter__(self):
        r"""
        Return an iterator for Dyck words with ``k1`` opening and ``k2``
        closing parentheses.

        EXAMPLES::

            sage: list(DyckWords(0))
            [[]]
            sage: list(DyckWords(1))
            [[1, 0]]
            sage: list(DyckWords(2))
            [[1, 0, 1, 0], [1, 1, 0, 0]]
            sage: len(DyckWords(5))
            42
            sage: list(DyckWords(3,2))
            [[1, 0, 1, 0, 1],
             [1, 0, 1, 1, 0],
             [1, 1, 0, 0, 1],
             [1, 1, 0, 1, 0],
             [1, 1, 1, 0, 0]]
        """
        if self.k1 == 0:
            yield self.element_class(self, [])
        elif self.k2 == 0:
            yield self.element_class(self, [open_symbol] * self.k1)
        else:
            for w in DyckWordBacktracker(self.k1, self.k2):
                yield self.element_class(self, w)

    def cardinality(self) -> int:
        r"""
        Return the number of Dyck words with `k_1` openers and `k_2` closers.

        This number is

        .. MATH::

            \frac{k_1 - k_2 + 1}{k_1 + 1} \binom{k_1 + k_2}{k_2}
            = \binom{k_1 + k_2}{k_2} - \binom{k_1 + k_2}{k_2 - 1}

        if `k_2 \leq k_1 + 1`, and `0` if `k_2 > k_1` (these numbers are the
        same if `k_2 = k_1 + 1`).

        EXAMPLES::

            sage: DyckWords(3, 2).cardinality()
            5
            sage: all(all(DyckWords(p, q).cardinality()
            ....:           == len(DyckWords(p, q).list()) for q in range(p + 1))
            ....:      for p in range(7))
            True
        """
        return (self.k1 - self.k2 + 1) * (self.k1 + self.k2).binomial(self.k2) // (self.k1 + 1)

################################################################
# Complete Dyck words


class CompleteDyckWords(DyckWords):
    """
    Abstract base class for all complete Dyck words.
    """
    Element = DyckWord_complete

    def __contains__(self, x) -> bool:
        r"""
        TESTS::

            sage: [] in DyckWords()
            True
            sage: [1] in DyckWords()
            False
            sage: [0] in DyckWords()
            False
        """
        if isinstance(x, DyckWord_complete):
            return True

        if not isinstance(x, list):
            return False

        if len(x) % 2:
            return False

        return is_a(x, len(x) // 2)

    def from_Catalan_code(self, code) -> DyckWord:
        r"""
        Return the Dyck word associated to the given Catalan code
        ``code``.

        A Catalan code of length `n` is a sequence
        `(a_1, a_2, \ldots, a_n)` of `n` integers `a_i` such that:

        - `0 \leq a_i \leq n-i` for every `i`;

        - if `i < j` and `a_i > 0` and `a_j > 0` and
          `a_{i+1} = a_{i+2} = \cdots = a_{j-1} = 0`,
          then `a_i - a_j < j-i`.

        It turns out that the Catalan codes of length `n` are in
        bijection with Dyck words.

        The Catalan code of a Dyck word is example (x) in Richard Stanley's
        exercises on combinatorial interpretations for Catalan objects.
        The code in this example is the reverse of the description provided
        there. See [Sta-EC2]_ and [StaCat98]_.

        EXAMPLES::

            sage: DyckWords().from_Catalan_code([])
            []
            sage: DyckWords().from_Catalan_code([0])
            [1, 0]
            sage: DyckWords().from_Catalan_code([0, 1])
            [1, 1, 0, 0]
            sage: DyckWords().from_Catalan_code([0, 0])
            [1, 0, 1, 0]
        """
        code = list(code)
        if not code:
            return self.element_class(self, [])
        res = self.from_Catalan_code(code[:-1])
        cuts = [0] + res.returns_to_zero()
        lst = [1] + res[:cuts[code[-1]]] + [0] + res[cuts[code[-1]]:]
        return self.element_class(self, lst)

    def from_area_sequence(self, code) -> DyckWord:
        r"""
        Return the Dyck word associated to the given area sequence
        ``code``.

        See :meth:`~DyckWord.to_area_sequence` for a definition of the area
        sequence of a Dyck word.

        .. SEEALSO:: :meth:`~DyckWord_complete.area`, :meth:`~DyckWord.to_area_sequence`.

        INPUT:

        - ``code`` -- list of integers satisfying ``code[0] == 0``
          and ``0 <= code[i+1] <= code[i]+1``

        EXAMPLES::

            sage: DyckWords().from_area_sequence([])
            []
            sage: DyckWords().from_area_sequence([0])
            [1, 0]
            sage: DyckWords().from_area_sequence([0, 1])
            [1, 1, 0, 0]
            sage: DyckWords().from_area_sequence([0, 0])
            [1, 0, 1, 0]
        """
        if not is_area_sequence(code):
            raise ValueError("the given sequence is not a sequence giving "
                             "the number of cells between the Dyck path "
                             "and the diagonal")
        dyck_word = []
        for i in range(len(code)):
            if i:
                dyck_word.extend([close_symbol] * (code[i - 1] - code[i] + 1))
            dyck_word.append(open_symbol)
        dyck_word.extend([close_symbol] * (2 * len(code) - len(dyck_word)))
        return self.element_class(self, dyck_word)

    def from_noncrossing_partition(self, ncp) -> DyckWord:
        r"""
        Convert a noncrossing partition ``ncp`` to a Dyck word.

        EXAMPLES::

            sage: DyckWord(noncrossing_partition=[[1,2]]) # indirect doctest
            [1, 1, 0, 0]
            sage: DyckWord(noncrossing_partition=[[1],[2]])
            [1, 0, 1, 0]

            sage: dws = DyckWords(5).list()
            sage: ncps = [x.to_noncrossing_partition() for x in dws]
            sage: dws2 = [DyckWord(noncrossing_partition=x) for x in ncps]
            sage: dws == dws2
            True
        """
        l = [0] * sum(len(v) for v in ncp)
        for v in ncp:
            l[max(v) - 1] = len(v)

        res = []
        for i in l:
            res += [open_symbol] + [close_symbol] * i
        return self.element_class(self, res)

    def from_non_decreasing_parking_function(self, pf) -> DyckWord:
        r"""
        Bijection from :class:`non-decreasing parking
        functions<sage.combinat.non_decreasing_parking_function.NonDecreasingParkingFunctions>`.

        See there the method
        :meth:`~sage.combinat.non_decreasing_parking_function.NonDecreasingParkingFunction.to_dyck_word`
        for more information.

        EXAMPLES::

            sage: D = DyckWords()
            sage: D.from_non_decreasing_parking_function([])
            []
            sage: D.from_non_decreasing_parking_function([1])
            [1, 0]
            sage: D.from_non_decreasing_parking_function([1,1])
            [1, 1, 0, 0]
            sage: D.from_non_decreasing_parking_function([1,2])
            [1, 0, 1, 0]
            sage: D.from_non_decreasing_parking_function([1,1,1])
            [1, 1, 1, 0, 0, 0]
            sage: D.from_non_decreasing_parking_function([1,2,3])
            [1, 0, 1, 0, 1, 0]
            sage: D.from_non_decreasing_parking_function([1,1,3,3,4,6,6])
            [1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0]

        TESTS::

            sage: D.from_non_decreasing_parking_function(NonDecreasingParkingFunction([]))
            []
            sage: D.from_non_decreasing_parking_function(NonDecreasingParkingFunction([1,1,3,3,4,6,6]))
            [1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0]
        """
        return self.from_area_sequence([i - pf[i] + 1 for i in range(len(pf))])


class CompleteDyckWords_all(CompleteDyckWords, DyckWords_all):
    """
    All complete Dyck words.
    """

    def _repr_(self) -> str:
        r"""
        TESTS::

            sage: DyckWords()
            Complete Dyck words
        """
        return "Complete Dyck words"

    def __iter__(self):
        """
        Iterate over ``self``.

        EXAMPLES::

            sage: it = DyckWords().__iter__()
            sage: [next(it) for x in range(10)]
            [[],
             [1, 0],
             [1, 0, 1, 0],
             [1, 1, 0, 0],
             [1, 0, 1, 0, 1, 0],
             [1, 0, 1, 1, 0, 0],
             [1, 1, 0, 0, 1, 0],
             [1, 1, 0, 1, 0, 0],
             [1, 1, 1, 0, 0, 0],
             [1, 0, 1, 0, 1, 0, 1, 0]]
        """
        n = 0
        while True:
            for x in CompleteDyckWords_size(n):
                yield self.element_class(self, list(x))
            n += 1

    class height_poset(UniqueRepresentation, Parent):
        r"""
        The poset of complete Dyck words compared componentwise by ``heights``.

        This is, ``D`` is smaller than or equal to ``D'`` if it is
        weakly below ``D'``.

        This is implemented by comparison of area sequences.
        """

        def __init__(self):
            r"""
            TESTS::

                sage: poset = DyckWords().height_poset()
                sage: TestSuite(poset).run()
            """
            Parent.__init__(self, facade=DyckWords_all(), category=Posets())

        def _an_element_(self):
            r"""
            TESTS::

                sage: DyckWords().height_poset().an_element()   # indirect doctest
                [1, 0, 1, 0, 1, 0, 1, 0, 1, 0]
            """
            return DyckWords_all().an_element()

        def __call__(self, obj):
            r"""
            TESTS::

                sage: poset = DyckWords().height_poset()
                sage: poset([1,0,1,0])
                [1, 0, 1, 0]
            """
            return DyckWords_all()(obj)

        def le(self, dw1, dw2):
            r"""
            Compare two Dyck words of equal size, and return ``True`` if
            all of the heights of ``dw1`` are less than or equal to the
            respective heights of ``dw2`` .

            .. SEEALSO::

                :meth:`~sage.combinat.dyck_word.DyckWord.to_area_sequence`

            EXAMPLES::

                sage: poset = DyckWords().height_poset()
                sage: poset.le(DyckWord([]), DyckWord([]))
                True
                sage: poset.le(DyckWord([1,0]), DyckWord([1,0]))
                True
                sage: poset.le(DyckWord([1,0,1,0]), DyckWord([1,1,0,0]))
                True
                sage: poset.le(DyckWord([1,1,0,0]), DyckWord([1,0,1,0]))
                False
                sage: [poset.le(dw1, dw2)
                ....:     for dw1 in DyckWords(3) for dw2 in DyckWords(3)]
                [True, True, True, True, True, False, True, False, True, True,
                False, False, True, True, True, False, False, False, True,
                True, False, False, False, False, True]
            """
            if len(dw1) != len(dw2):
                raise ValueError(f"length mismatch: {dw1} and {dw2}")
            ar1 = dw1._area_sequence_iter()
            ar2 = dw2._area_sequence_iter()
            return all(a1 <= a2 for a1, a2 in zip(ar1, ar2))


class CompleteDyckWords_size(CompleteDyckWords, DyckWords_size):
    """
    All complete Dyck words of a given size.
    """

    def __init__(self, k):
        """
        Initialize ``self``.

        TESTS::

            sage: TestSuite(DyckWords(4)).run()
        """
        CompleteDyckWords.__init__(self, category=FiniteEnumeratedSets())
        DyckWords_size.__init__(self, k, k)

    def __contains__(self, x) -> bool:
        r"""
        TESTS::

            sage: [1, 0] in DyckWords(1)
            True
            sage: [1, 0] in DyckWords(2)
            False
            sage: [1, 1, 0, 0] in DyckWords(2)
            True
            sage: [1, 0, 0, 1] in DyckWords(2)
            False
        """
        return CompleteDyckWords.__contains__(self, x) and len(x) // 2 == self.k1

    def cardinality(self) -> int:
        r"""
        Return the number of complete Dyck words of semilength `n`, i.e. the
        `n`-th :func:`Catalan number<sage.combinat.combinat.catalan_number>`.

        EXAMPLES::

            sage: DyckWords(4).cardinality()
            14
            sage: ns = list(range(9))
            sage: dws = [DyckWords(n) for n in ns]
            sage: all(dw.cardinality() == len(dw.list()) for dw in dws)
            True
        """
        return catalan_number(self.k1)

    def random_element(self) -> DyckWord:
        """
        Return a random complete Dyck word of semilength `n`.

        The algorithm is based on a classical combinatorial fact. One
        chooses at random a word with `n` 0s and `n+1` 1s. One then
        considers every 1 as an ascending step and every 0 as a
        descending step, and one finds the lowest point of the path
        (with respect to a slightly tilted slope). One then cuts the
        path at this point and builds a Dyck word by exchanging the
        two parts of the word and removing the initial step.

        .. TODO::

            extend this to m-Dyck words

        EXAMPLES::

            sage: dw = DyckWords(8)
            sage: dw.random_element()  # random
            [1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0]

            sage: D = DyckWords(8)
            sage: D.random_element() in D
            True
        """
        from sage.misc.prandom import shuffle
        n = self.k1
        w = [0] * n + [1] * (n + 1)
        shuffle(w)
        idx = 0
        height = 0
        height_min = 0
        for i in range(2 * n):
            if w[i] == 1:
                height += n
            else:
                height -= n + 1
                if height < height_min:
                    height_min = height
                    idx = i + 1
        w = w[idx:] + w[:idx]
        return self.element_class(self, w[1:])

    def _iter_by_recursion(self):
        """
        Iterate over ``self`` by recursively using the position of
        the first return to 0.

        EXAMPLES::

            sage: DW = DyckWords(4)
            sage: L = [d for d in DW._iter_by_recursion()]; L
            [[1, 0, 1, 0, 1, 0, 1, 0],
             [1, 0, 1, 0, 1, 1, 0, 0],
             [1, 0, 1, 1, 0, 0, 1, 0],
             [1, 0, 1, 1, 0, 1, 0, 0],
             [1, 0, 1, 1, 1, 0, 0, 0],
             [1, 1, 0, 0, 1, 0, 1, 0],
             [1, 1, 0, 0, 1, 1, 0, 0],
             [1, 1, 0, 1, 0, 0, 1, 0],
             [1, 1, 1, 0, 0, 0, 1, 0],
             [1, 1, 0, 1, 0, 1, 0, 0],
             [1, 1, 0, 1, 1, 0, 0, 0],
             [1, 1, 1, 0, 0, 1, 0, 0],
             [1, 1, 1, 0, 1, 0, 0, 0],
             [1, 1, 1, 1, 0, 0, 0, 0]]
            sage: len(L) == DW.cardinality()
            True
        """
        # Do a couple of small cases first
        if self.k1 <= 2:
            if self.k1 == 0:
                yield self.element_class(self, [])
            elif self.k1 == 1:
                yield self.element_class(self, [1, 0])
            elif self.k1 == 2:
                yield self.element_class(self, [1, 0, 1, 0])
                yield self.element_class(self, [1, 1, 0, 0])
            return

        # Create all necessary parents
        P = [CompleteDyckWords_size(i) for i in range(self.k1)]

        for i in range(self.k1):
            for p in P[i]._iter_by_recursion():
                list_1p0 = [1] + list(p) + [0]
                for s in P[-i - 1]._iter_by_recursion():
                    yield self.element_class(self, list_1p0 + list(s))


def is_area_sequence(seq) -> bool:
    r"""
    Test if a sequence `l` of integers satisfies `l_0 = 0` and
    `0 \leq l_{i+1} \leq l_i + 1` for `i > 0`.

    EXAMPLES::

        sage: from sage.combinat.dyck_word import is_area_sequence
        sage: is_area_sequence([0,2,0])
        False
        sage: is_area_sequence([1,2,3])
        False
        sage: is_area_sequence([0,1,0])
        True
        sage: is_area_sequence([0,1,2])
        True
        sage: is_area_sequence([])
        True
    """
    if not seq:
        return True
    return seq[0] == 0 and all(0 <= seq[i + 1] <= seq[i] + 1
                               for i in range(len(seq) - 1))


def is_a(obj, k1=None, k2=None) -> bool:
    r"""
    Test if ``obj`` is a Dyck word with exactly ``k1`` open symbols and
    exactly ``k2`` close symbols.

    If ``k1`` is not specified, then the number of open symbols can be
    arbitrary. If ``k1`` is specified but ``k2`` is not, then ``k2`` is
    set to be ``k1``.

    EXAMPLES::

        sage: from sage.combinat.dyck_word import is_a
        sage: is_a([1,1,0,0])
        True
        sage: is_a([1,0,1,0])
        True
        sage: is_a([1,1,0,0], 2)
        True
        sage: is_a([1,1,0,0], 3)
        False
        sage: is_a([1,1,0,0], 3, 2)
        False
        sage: is_a([1,1,0])
        True
        sage: is_a([0,1,0])
        False
        sage: is_a([1,0,0])
        False
        sage: is_a([1,1,0],2,1)
        True
        sage: is_a([1,1,0],2)
        False
        sage: is_a([1,1,0],1,1)
        False
    """
    if k1 is not None:
        if k2 is None:
            k2 = k1
        elif k1 < k2:
            raise ValueError("k1 (= %s) must be >= k2 (= %s)" % (k1, k2))

    n_opens = 0
    n_closes = 0

    for p in obj:
        if p == open_symbol:
            n_opens += 1
        elif p == close_symbol:
            if n_opens == n_closes:
                return False
            n_closes += 1
        else:
            return False

    return (k1 is None and k2 is None) or (n_opens == k1 and n_closes == k2)


def pealing(D, return_touches=False):
    r"""
    A helper function for computing the bijection from a Dyck word to a
    `231`-avoiding permutation using the bijection "Stump".  For details
    see [Stu2008]_.

    .. SEEALSO::

        :meth:`~sage.combinat.dyck_word.DyckWord_complete.to_noncrossing_partition`

    EXAMPLES::

        sage: from sage.combinat.dyck_word import pealing
        sage: pealing(DyckWord([1,1,0,0]))
        [1, 0, 1, 0]
        sage: pealing(DyckWord([1,0,1,0]))
        [1, 0, 1, 0]
        sage: pealing(DyckWord([1, 1, 0, 0, 1, 1, 1, 0, 0, 0]))
        [1, 0, 1, 0, 1, 0, 1, 0, 1, 0]
        sage: pealing(DyckWord([1,1,0,0]),return_touches=True)
        ([1, 0, 1, 0], [[1, 2]])
        sage: pealing(DyckWord([1,0,1,0]),return_touches=True)
        ([1, 0, 1, 0], [])
        sage: pealing(DyckWord([1, 1, 0, 0, 1, 1, 1, 0, 0, 0]),return_touches=True)
        ([1, 0, 1, 0, 1, 0, 1, 0, 1, 0], [[1, 2], [3, 5]])
    """
    n = D.semilength()
    area = D.to_area_sequence()
    new_area = []
    touch_sequences = []
    touches = []
    for i in range(n - 1):
        if area[i + 1] == 0:
            touches.append(i + 1)
            if len(touches) > 1:
                touch_sequences.append(touches)
                touches = []
            elif area[i] == 0:
                touches = []
            new_area.append(0)
        elif area[i + 1] == 1:
            new_area.append(0)
            touches.append(i + 1)
        else:
            new_area.append(area[i + 1] - 2)
    new_area.append(0)
    if area[n - 1]:
        touches.append(n)
        if len(touches) > 1:
            touch_sequences.append(touches)
    D = DyckWords().from_area_sequence(new_area)
    if return_touches:
        return (D, touch_sequences)
    else:
        return D


from sage.misc.persist import register_unpickle_override
register_unpickle_override('sage.combinat.dyck_word', 'DyckWord', DyckWord)
