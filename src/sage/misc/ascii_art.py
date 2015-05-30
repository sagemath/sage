# -*- coding: utf-8 -*-
r"""
ASCII Art

This file contains:

- :class:`AsciiArt` an simple implementation of an ASCII art object,
- :func:`ascii_art` a function to get the ASCII art representation of any
  object in Sage,
- several others functions use to get ASCII art representation of primitive
  python elements (list, tuple, dict, set).

AUTHOR:

- Jean-Baptiste Priez (2013-04): initial version

EXAMPLES::

    sage: n = var('n')
    sage: integrate(n^2/x,x)
    n^2*log(x)
    sage: ascii_art(integrate(n^2/x,x))
     2
    n *log(x)
    sage: ascii_art(integrate(n^2/(pi*x),x))
     2
    n *log(x)
    ---------
        pi
    sage: ascii_art(list(Partitions(6)))
    [                                                       * ]
    [                                                   **  * ]
    [                                      ***      **  *   * ]
    [                      ****       ***  *    **  **  *   * ]
    [         *****  ****  *     ***  **   *    **  *   *   * ]
    [ ******, *    , **  , *   , ***, *  , *  , **, * , * , * ]

This method :meth:`ascii_art` could be automatically use by the display hook
manager activated by the magic function: ``%display ascii_art``::

    sage: from sage.repl.interpreter import get_test_shell
    sage: shell = get_test_shell()
    sage: shell.run_cell('%display ascii_art')
    sage: shell.run_cell("i = var('i')")
    sage: shell.run_cell('sum(factorial(i)*x^i, i, 0, 10)')
             10           9          8         7        6        5       4      3
    3628800*x   + 362880*x  + 40320*x  + 5040*x  + 720*x  + 120*x  + 24*x  + 6*x
    <BLANKLINE>
         2
    + 2*x  + x + 1
    sage: shell.run_cell('3/(7*x)')
     3
    ---
    7*x
    sage: shell.run_cell('list(Compositions(5))')
    [ *
    [ *  **   *        *                   *
    [ *  *   **  ***   *   **    *         *   **    *          *
    [ *  *   *   *    **  **   ***  ****   *   *    **   ***    *    **     *
    [ *, * , * , *  , * , *  , *  , *   , **, ** , ** , **  , ***, *** , ****,
    <BLANKLINE>
          ]
          ]
          ]
          ]
    ***** ]
    sage: shell.run_cell('%display simple')
    sage: shell.quit()

::

                                .      ,    ,              ...
                                .?  .~?$NNO.II7.        ..GOGG
                                 ., ~7NNI7NDG$       ...~~~MGG
                  ..:IG..       ...G7:DDGNNDO..   ...7:~~~~GNN
                  .~~GGGGG...   .O .=$+OD7GI$ ...GG:~+~~~~MG?M.7?
                  .~~~D~~GGGGDGOM$..~+IN=NDG.G::D~~~?~~~7?NDI:???G
                  .~~~~G:~~G~D:~NONGMMOGD$ND~~::I77~~$?++MN7G GI?D.
                   7+I~~DG~~~N:GNGGOD$7OIOOMII::~:I7?G+OGN7G $III?.
                    .7==I~~~++IN?+++?$7ND$==$G:G??????N$OG$  ~III?7
                     .$$+++?GG=+G?GDM?GOG:NGMGN7??????D$D~GODIII???.
                       D++++NG+DD$+=GGONGI$DNGDN??????,IN=DOGODN???.
                     , +~+++N??D$I+I=GNMDGDGNIINN??D+:::O~$GGGDIG?7
                    :DD.:::OG?G=I~~G7GNM7777NIGODDNNG:::I$GGGGG7??O.
                   ~GDDM:::GGMG+?+~$NNN$7IG7NMNMMMN, =::I7GGGGGO???~
                  :GDND7?::OOM.D+O~GGDI77777DMGNM$.  ~:,$IGGGGGO???DO.
                 OGDDNN.D77OO. $?7==GG777$7GNGMGO.     NOIGGGGGO???G$.
               .OODNNN,DIGDM$GGGGG==GGGGIGNDMDMG,      IGIGGGDON???GG:
              .GODNDM.G$I$IOIOMG$G$?GGGIO,7OD7GG.     ,GDIGGG??????GGG.
             .DGDNI7I7MDI+OOODN$$O,$7DMIN,,IOO77O.   G?$DIGGG?ID???OGG
             GGDNNMMO7GD+OOOGIMOG7::NN====:?MMNGIDD,..IINIGGG?I??DIOGG?
           .7ODMMMN.G7IOGOOODIMG,,:::$=~==::7OGG~IGOMDGMNIGGG????.OG$G.
           ?ODDMNNNO,II$OGODMGDMM?:DMG==~MDINNM$.7$IONDGI?GGG????:$GGG.
        .$MDMNNNNN..:?7GDDDGG,GGM?~:GGGDNND.GIM7D+GI$ON.:?GGG????$$OGG.
      .7DNDDDNMDOGG.=IGGND=7II+N??::$GIIO,IIGMG?7I7G$ON?,IGGG????7GIGG.
      ~GGGDNMMOOGGG$MGMMGDGMDGM?,G:GNG,:IIIGDG7IGGGGG$+NMIGGO?????IGGG.
     .GGGDDMM7OGNGMODMNDDDOO.MII?GI$7IIIIING7GGDMM.IGDG.G7GGGG??$?7GGG
    .IGDDDDOINNGMGMDMNDDGDO...$OI+??7OIIIDDGN+$==I=GD ?,NGGGGN7????GGG
     7GGNDM$GONDGDD$MMGNDN. G.:$$G$?$7II$GOO,O=+7O7O~N?OM?GGGGGD+??$GG
     OOGDOI7DGMG..=~DG$DD.  $,,$$$D??7ODOOOODG$777$G OGMM:GGGGGG??GGG.
    .$GDDG?7DMOGDNNMGGDGG ..,=O7$GG$+O$$OG=+O:GI77G$. ...DIGGGGI?$GGG.
    .OGGGGO7OGDGGNGGGDGG.  O:,77$O$$$D     $GN:7777GD   ..$GGGG??GGG$.
    .GGGGDI GD.~NOGDGGG    MG,777G77D7      +$~D77777=   ..7GGG???GGO.
      GGD.   .~NODN...   ..D::77O$7G77      .OG:G???G7.    7$NOM??GG=.
        .     . .      .+~~DD~77$G7OD7      .?GGO?$OM$D     ?????DGG...
                       7I77DN$II7$M.G$G      . .:G?==$G     .?????D??~
                     . .:IIIGO7O$GN..O$         7I=~+?I.    =???7? GGIG
                      $,III7NGGNNNG, .O.         O====7I.  .DG??$?. GG?G
                      .$7III$77777$$~  .          +====D$GG ~$???:  .$GIG.
                     :.+IIII$77$G$$O,+            OI==+7I$7=.:???     .GG$.
                     ,,III$I7?I7$$OG7D            ?+,==+77G?$,???.      GG?
                     ?IIIII=+$I77777             .:++=+++O77:,??:         GG.
                     .IIIIG+III77777.             7,=,++?II7$,?G          .GO .
                      $OII7?OIIIIIO              .I:++??IIIII???.           I$ .
                      .DNGDMOIG777.                +~++?I$IIMN?,             .G,
                       .G.I?IIDGII~                .OG??$$MIMGI.               G,
                        N+?II7OGI7.                  DD?$GGDGID.                 :
                         =?IIOOI$.                    ?~~GGG$II                  .+.
                        ,7=IINGI                       ?=IGGGI?.
                        G,+IGOII                        ?+II7??+
                       .,:?IGOII                         ?+GIOII
                     .:N+=?IMGI7=.                       .~:$I?IO.
                   .$:IGO?$IIIII7+.                       O::I7I?.
                  .:::=IIGIIIIIII.                        .~,$IG?IO
                 .+:$IIIIMIII7G$?.                         I::I7IIG7
                 .$I$+7IIIMMMG7I7                          ?G,DGGNII.
                 .$~:$IIGMNGND77I:                          +I$GOD=?=
                ?=?IIIIGIIIINI7777?                        .7~=~===?.
              ONGNDG??IG?III$N7I777                          D~====I
            .:$??7IIIIIII.....,....                          O::==~$D.
            .,........ ..                                  ..M:I.==$7G.
                                                            I?::IIIII7.
                                                           .~:G:IIIIIG.
                                                           .$:,O7III$$O
                                                             ::~DOGGNNO
                                                             .::,IOODI,
                                                               .7????$.
                                                                 ... .
"""
#*******************************************************************************
#       Copyright (C) 2013 Jean-Baptiste Priez <jbp@kerios.fr>,
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
#                  http://www.gnu.org/licenses/
#*******************************************************************************

import os, sys
from sage.structure.sage_object import SageObject
from sage.rings.integer import Integer


################################################################################
### Global variable use to compute the maximal length allows for ascii art
### object.
MAX_WIDTH = None
################################################################################

class AsciiArt(SageObject):
    r"""
    An Ascii art object is an object with some specific representation for
    *printing*.

    INPUT:

    - ``lines`` -- the list of lines of the representation of the ascii art
      object
    - ``breakpoints`` -- the list of points where the representation can be
      split
    - ``baseline`` -- the reference line (from the bottom)
    - ``atomic`` -- indicate if the ascii art representation is splittable
      (must be coherent with breakpoints)

    EXAMPLES::

        sage: i = var('i')
        sage: ascii_art(sum(pi^i/factorial(i)*x^i, i, 0, oo))
         pi*x
        e
    """
    def __init__(self, lines=[], breakpoints=[], baseline=None, atomic=True):
        r"""
        TESTS::

            sage: from sage.misc.ascii_art import AsciiArt
            sage: aao = AsciiArt()
            sage: aao.is_atomic()
            True
            sage: aao
            <BLANKLINE>
            sage: aa = AsciiArt(["  *  ", " * * ", "*****"]); aa
              *
             * *
            *****
        """
        self._is_uniq = True
        self._matrix = lines
        self._breakpoints = breakpoints
        self._baseline = baseline if baseline is not None else 0
        self._is_atomic = atomic

        self._h = len(lines)
        self._l = max(map(lambda line: len(line), lines) + [0])

    def __getitem__(self, key):
        r"""
        Return the line `key` of the ASCII art object.

        TESTS::

            sage: from sage.misc.ascii_art import AsciiArt
            sage: p5 = AsciiArt(["  *  ", " * * ", "*****"])
            sage: p5[1]
            ' * * '
        """
        return self._matrix[key]

    def __iter__(self):
        r"""
        Iterator on all lines of the ASCII art object.

        TESTS::

            sage: from sage.misc.ascii_art import AsciiArt
            sage: p5 = AsciiArt(["  *  ", " * * ", "*****"])
            sage: for line in p5:
            ....:     print line
              *
             * *
            *****
        """
        for elem in self._matrix:
            yield elem

    def _repr_(self):
        r"""
        TESTS::

            sage: from sage.misc.ascii_art import AsciiArt
            sage: p5 = AsciiArt(["  *  ", " * * ", "*****"])
            sage: repr(p5)
            '  *  \n * * \n*****'
        """
        # Compute the max length of a draw
        global MAX_WIDTH
        if MAX_WIDTH is not None:
            hsize = MAX_WIDTH
        else:
            hsize = self._terminal_width()
        #########
        # if the draw is larger than the max length it try to split...
        if hsize <= self._l and len(self._breakpoints) > 0:
            return self._split_repr_(hsize)
        #########
        output = ""
        if len(self._matrix) > 0:
            for i in range(len(self._matrix) - 1):
                output += self._matrix[i] + "\n"
            return output + self._matrix[len(self._matrix) - 1]
        return output


    def is_atomic(self):
        r"""
        Return ``True`` if the :class:`AsciiArt` object is not splitable and
        ``False`` in otherwise.

        For example, we considere a linear expression::

            sage: a = 14*x^5 + 5*x^4
            sage: ascii_art(a)
                5      4
            14*x  + 5*x

        If ASCII art object is not atomic, it is splittable on the ``+``
        (in fact it is not really true because we use ``sympy`` to make
        ASCII art).

        TESTS::

            sage: from sage.misc.ascii_art import AsciiArt
            sage: aa = AsciiArt(["  *  ", " * * ", "*****"]); aa
              *
             * *
            *****
            sage: aa.is_atomic()
            True
            sage: laa = ascii_art([aa,aa])
            sage: laa.is_atomic()
            False
        """
        return self._is_atomic

    def get_baseline(self):
        r"""
        Return the line where the baseline is, for example::

                5      4
            14*x  + 5*x

        the baseline has at line `0` and ::

            { o       }
            {  \  : 4 }
            {   o     }

        has at line `1`.

        TESTS::

            sage: from sage.misc.ascii_art import AsciiArt
            sage: aa = AsciiArt(["   *   ", "  * *  ", " *   * ", "*******"], baseline=1);aa
               *
              * *
             *   *
            *******
            sage: aa.get_baseline()
            1
            sage: b = AsciiArt(["<-"])
            sage: aa+b
               *
              * *
             *   * <-
            *******
        """
        return self._baseline

    def get_breakpoints(self):
        r"""
        Return an iterator of breakpoints where the object can be split.

        For example the expression::

               5    4
            14x + 5x

        can be split on position 4 (on the ``+``).

        EXAMPLES::

            sage: from sage.misc.ascii_art import AsciiArt
            sage: p3 = AsciiArt([" * ", "***"])
            sage: p5 = AsciiArt(["  *  ", " * * ", "*****"])
            sage: aa = ascii_art([p3, p5])
            sage: aa.get_breakpoints()
            [2, 5, 6, 7, 12]
        """
        return self._breakpoints

    def _isatty(self):
        """
        Test whether stdout is a TTY

        If this test succeeds, you can assume that stdout is directly
        connected to a terminal. Otherwise you should treat stdout as
        being redirected to a file.

        OUTPUT:

        Boolean

        EXAMPLES::

            sage: from sage.misc.ascii_art import empty_ascii_art
            sage: empty_ascii_art._isatty()
            False
        """
        from sage.doctest import DOCTEST_MODE
        if DOCTEST_MODE:
            return False
        try:
            return os.isatty(sys.stdout.fileno())
        except Exception:
            # The IPython zeromq kernel uses a fake stdout that does
            # not support fileno()
            return False
    
    def _terminal_width(self):
        """
        Compute the width size of the terminal.

        EXAMPLES::

            sage: from sage.misc.ascii_art import empty_ascii_art
            sage: empty_ascii_art._terminal_width()
            80
        """
        if not self._isatty():
            return 80
        import fcntl, termios, struct
        rc = fcntl.ioctl(int(0), termios.TIOCGWINSZ,
                         struct.pack('HHHH', sys.stdout.fileno(), 0, 0, 0))
        h, w, hp, wp = struct.unpack('HHHH', rc)
        return w

    def _split_repr_(self, size):
        r"""
        Split the draw and the left part has length ``size``.

        TESTS::

            sage: from sage.misc.ascii_art import AsciiArt
            sage: p3 = AsciiArt([" * ", "***"])
            sage: p5 = AsciiArt(["  *  ", " * * ", "*****"])
            sage: aa = ascii_art([p3, p5])
            sage: print aa._split_repr_(6)
            [
            [  *
            [ ***
            <BLANKLINE>
                *   ]
               * *  ]
            , ***** ]
        """
        import sys
        f_split = self._breakpoints[0]; i = 1
        while i < len(self._breakpoints) and self._breakpoints[i] < size:
            f_split = self._breakpoints[i]
            i += 1
        if size <= f_split:
            import warnings
            warnings.warn("the console size is smaller than the pretty" +
                "representation of the object")
        top, bottom = self.split(f_split)
        return (top * AsciiArt([""])).__repr__() + "\n" + bottom.__repr__()

    def split(self, pos):
        r"""
        Split the representation at the position ``pos``.

        EXMAPLES::

            sage: from sage.misc.ascii_art import AsciiArt
            sage: p3 = AsciiArt([" * ", "***"])
            sage: p5 = AsciiArt(["  *  ", " * * ", "*****"])
            sage: aa = ascii_art([p3, p5])
            sage: a,b= aa.split(6)
            sage: a
            [
            [  *
            [ ***,
            sage: b
               *   ]
              * *  ]
             ***** ]
        """
        left = []; right = []
        for line in self:
            left.append(line[:pos])
            right.append(line[pos:])
        l_bp = []; r_bp = []
        for bp in self._breakpoints:
            if bp < pos:
                l_bp.append(bp)
            elif bp > pos:
                r_bp.append(bp - pos)
        return AsciiArt(left, l_bp), AsciiArt(right, r_bp)

    @staticmethod
    def _compute_new_baseline(obj1, obj2):
        r"""
        TESTS::

            sage: from sage.misc.ascii_art import AsciiArt
            sage: l5 = AsciiArt(lines = ['|' for _ in range(5)], baseline = 2); l5
            |
            |
            |
            |
            |
            sage: l3 = AsciiArt(lines = ['|' for _ in range(3)], baseline = 1); l3
            |
            |
            |
            sage: AsciiArt._compute_new_baseline(l5, l3)
            2
            sage: l5 + l3
            |
            ||
            ||
            ||
            |
            sage: l5._baseline = 0
            sage: AsciiArt._compute_new_baseline(l5, l3)
            1
            sage: l5 + l3
            |
            |
            |
            ||
            ||
             |
            sage: l5._baseline = 4
            sage: AsciiArt._compute_new_baseline(l5, l3)
            4
            sage: l5 + l3
             |
            ||
            ||
            |
            |
            |
            sage: l3._baseline = 0
            sage: AsciiArt._compute_new_baseline(l3, l5)
            4
            sage: l3 + l5
            |
            |
            ||
             |
             |
             |
             |
            sage: l5._baseline = None
            sage: AsciiArt._compute_new_baseline(l3, l5)
            2
            sage: l3._baseline = 2
            sage: AsciiArt._compute_new_baseline(l3, l5)
            4
            sage: l3 + l5
            ||
            ||
            ||
             |
             |

        """
        if obj1.get_baseline() is None:
            if obj2.get_baseline() is None:
                return None
            return obj2.get_baseline() + max(obj1._h - obj2._h, 0)
        if obj2.get_baseline() is None:
            return obj1.get_baseline() + max(obj2._h - obj1._h, 0)
        return max(
            obj1.get_baseline(),
            obj2.get_baseline()
        )

    @staticmethod
    def _compute_new_h(obj1, obj2):
        r"""
        TESTS::

            sage: from sage.misc.ascii_art import AsciiArt
            sage: l5 = AsciiArt(lines=['|' for _ in range(5)], baseline=2); l5
            |
            |
            |
            |
            |
            sage: l3 = AsciiArt(lines=['|' for _ in range(3)], baseline=1); l3
            |
            |
            |
            sage: AsciiArt._compute_new_h(l5, l3)
            5
            sage: l5 + l3
            |
            ||
            ||
            ||
            |
            sage: l5._baseline = 0
            sage: AsciiArt._compute_new_h(l5, l3)
            6
            sage: l5 + l3
            |
            |
            |
            ||
            ||
             |
            sage: l5._baseline = 4
            sage: AsciiArt._compute_new_h(l5, l3)
            6
            sage: l5 + l3
             |
            ||
            ||
            |
            |
            |
            sage: l3._baseline = 0
            sage: AsciiArt._compute_new_h(l3, l5)
            7
            sage: l3 + l5
            |
            |
            ||
             |
             |
             |
             |
        """
        if obj1.get_baseline() is None or obj2.get_baseline() is None:
            return max(obj1._h, obj2._h)
        return max(
            obj1.get_baseline(),
            obj2.get_baseline()
        ) + max(
            obj1._h - obj1.get_baseline(),
            obj2._h - obj2.get_baseline()
        )

    def __len__(self):
        r"""
        Return the length of the ASCII art object.

        TESTS::

            sage: from sage.misc.ascii_art import AsciiArt
            sage: p3 = AsciiArt([" * ", "***"])
            sage: len(p3)
            3
            sage: p5 = AsciiArt(["  *  ", " * * ", "*****"])
            sage: len(p5)
            5
        """
        return self._l

    def __add__(self, Nelt):
        r"""
        Concatenate two ascii art object.

        By default, when two object are concatenated, the new one will be
        splittable between both.

        If the baseline is defined, the concatenation is computed such that the
        new baseline coincidate with the olders.

        For example, let `T` be a tree with it's baseline ascii art
        representation in the middle::

            o
             \
              o
             / \
            o   o

        and let `M` be a matrix with it's baseline ascii art representation at
        the middle two::

            [1 2 3]
            [4 5 6]
            [7 8 9]

        then the concatenation of both will give::

            o
             \   [1 2 3]
              o  [4 5 6]
             / \ [7 8 9]
            o   o

        If one of the objects has not baseline, the concatenation is realized
        from the top::

            o    [1 2 3]
             \   [4 5 6]
              o  [7 8 9]
             / \
            o   o

        TESTS::

            sage: from sage.misc.ascii_art import AsciiArt
            sage: l5 = AsciiArt(lines=['|' for _ in range(5)], baseline=2); l5
            |
            |
            |
            |
            |
            sage: l3 = AsciiArt(lines=['|' for _ in range(3)], baseline=1); l3
            |
            |
            |
            sage: l3 + l5
             |
            ||
            ||
            ||
             |
            sage: l5 + l3
            |
            ||
            ||
            ||
            |
            sage: l5._baseline = 0
            sage: l5 + l3
            |
            |
            |
            ||
            ||
             |
            sage: l5._baseline = 4
            sage: l5 + l3
             |
            ||
            ||
            |
            |
            |
            sage: l3._baseline = 0
            sage: l3 + l5
            |
            |
            ||
             |
             |
             |
             |
        """
        new_matrix = []
        new_h = AsciiArt._compute_new_h(self, Nelt)
        new_baseline = AsciiArt._compute_new_baseline(self, Nelt)

        if self._baseline is not None and Nelt._baseline is not None:
            # left traitement
            for line in self._matrix:
                new_matrix.append(line + " "** Integer(self._l - len(line)))

            if new_h > self._h:
                # |                 new_h > self._h
                # |                 new_baseline > self._baseline
                # ||<-- baseline    number of white lines at the bottom
                #  | }               :: Nelt._baseline - self._baseline
                #  | }
                if new_baseline > self._baseline:
                    for k in range(new_baseline - self._baseline):
                        new_matrix.append(" " ** Integer(self._l))
                #  | }              new_h > self._h
                #  | }              new_h - new_baseline > self._h - self._baseline
                # ||<-- baseline    number of white lines at the top
                # ||                :: new_h - new_baseline - self._h + self._baseline
                # ||
                #  |
                #  |
                if new_h - new_baseline > self._h - self._baseline:
                    for _ in range((new_h - new_baseline) - (self._h - self._baseline)):
                        new_matrix.insert(0, " " ** Integer(self._l))

            # right traitement
            i = 0
            if new_h > Nelt._h:
                # |  }              new_h > Nelt._h
                # |  }              new_h - new_baseline > Nelt._h - self._baseline
                # ||<-- baseline    number of white lines at the top
                # ||                :: new_h - new_baseline - Nelt._h + Nelt._baseline
                # ||
                # ||
                # |

                i = max(new_h - new_baseline - Nelt._h + Nelt._baseline , 0)

            for j in range(Nelt._h):
                new_matrix[i+j] += Nelt._matrix[j]
        else:
            for line in self._matrix:
                new_matrix.append(line + " " ** Integer(self._l - len(line)))
            for i, line_i in enumerate(Nelt._matrix):
                if i == len(new_matrix):
                    new_matrix.append(" "**Integer(self._l) + line_i)
                else: new_matrix[i] += line_i

        # breakpoint
        new_breakpoints = list(self._breakpoints)
        new_breakpoints.append(self._l)
        for bp in Nelt._breakpoints:
            new_breakpoints.append(bp + self._l)
        from sage.misc.misc import uniq
        return AsciiArt(
            lines = new_matrix,
            breakpoints = uniq(new_breakpoints),
            baseline = new_baseline,
            atomic = False
        )

    def __mul__(self, Nelt):
        r"""
        The operator ``*`` is use to the representation ``self`` at
        the top of an other ``Nelt``.

        TESTS::

            sage: from sage.misc.ascii_art import AsciiArt
            sage: cub = AsciiArt(lines=['***' for _ in range(3)]); cub
            ***
            ***
            ***
            sage: pyr = AsciiArt(lines=[' ^ ', '/ \\', '---']); pyr
             ^
            / \
            ---
            sage: cub * pyr
            ***
            ***
            ***
             ^
            / \
            ---
        """
        new_repr = AsciiArt(self._matrix + Nelt._matrix)
        return new_repr

############## PRIMITIVES ################

def ascii_art(obj):
    r"""
    Return an ASCII art reprensentation of ``obj``::

        sage: ascii_art(integral(exp(x+x^2)/(x+1), x))
            /
           |
           |   2
           |  x  + x
           | e
           | ------- dx
           |  x + 1
           |
          /

    TESTS::

        sage: n = var('n')
        sage: ascii_art(sum(binomial(2 * n, n + 1) * x^n, n, 0, oo))
         /        __________    \
        -\2*x + \/ -4*x + 1  - 1/
        --------------------------
                   __________
             2*x*\/ -4*x + 1
        sage: ascii_art(list(DyckWords(3)))
        [                                   /\   ]
        [            /\    /\      /\/\    /  \  ]
        [ /\/\/\, /\/  \, /  \/\, /    \, /    \ ]
        sage: ascii_art(1)
        1
    """
    if isinstance(obj, (tuple, list, dict, set)):
        if obj.__class__ is tuple:
            return ascii_art_tuple(obj)
        elif obj.__class__ is dict:
            return ascii_art_dict(obj)
        elif obj.__class__ is list:
            return ascii_art_list(obj)
        else:
            return ascii_art_set(obj)
    if isinstance(obj, SageObject):
        res = obj._ascii_art_()
    else: res = AsciiArt(str(obj).splitlines())
    return res

def _ascii_art_iter(iter, func=lambda elem, _: ascii_art(elem)):
    r"""
    Intermediate function for ``ascii_art_X`` where ``X`` is ``dict``,
    ``set``, ``list``, or ``tuple``.

    TESTS::

        sage: ascii_art(list(DyckWords(3)))
        [                                   /\   ]
        [            /\    /\      /\/\    /  \  ]
        [ /\/\/\, /\/  \, /  \/\, /    \, /    \ ]
    """
    repr_elems = AsciiArt([""])
    bool = False
    for elem in iter:
        if bool:
            #print repr_elems.get_baseline()
            if repr_elems._baseline is not None:
                repr_elems += AsciiArt(
                    ["" ** Integer(repr_elems._h - 1 - repr_elems._baseline) ] +
                    [", "],
                    baseline=repr_elems._baseline)
            else:
                repr_elems += AsciiArt(
                    (["" ** Integer(repr_elems._h - 1) ] if repr_elems._h > 1 else [])+
                    [", "],
                    baseline=repr_elems._baseline)
            repr_elems._breakpoints.append(repr_elems._l - 1)
        repr_elems += func(elem, iter)
        bool = True
    return repr_elems

def ascii_art_set(set):
    r"""
    Return an ASCII art output of a set.

    TESTS::

        sage: ascii_art(set(DyckWords(3)))
        {                                   /\   }
        {  /\      /\/\              /\    /  \  }
        { /  \/\, /    \, /\/\/\, /\/  \, /    \ }
    """
    repr_elems = _ascii_art_iter(set)
    return _ascii_art_with_borders(repr_elems, "{ ", " }")

def ascii_art_dict(dict):
    r"""
    Return an ASCII art output of a dictionnary.

    TESTS::

        sage: ascii_art({i:dw for i,dw in enumerate(DyckWords(3))})
        {                                             /\   }
        {                /\      /\        /\/\      /  \  }
        { 0:/\/\/\, 1:/\/  \, 2:/  \/\, 3:/    \, 4:/    \ }
    """
    def func(k, dict):
        return ascii_art(k) + AsciiArt([":"], baseline=0) + ascii_art(dict[k])
    repr_elems = _ascii_art_iter(dict, func)
    return _ascii_art_with_borders(repr_elems, "{ ", " }")

def _ascii_art_with_borders(repr_elems, l_border, r_border):
    r"""
    Intermediate function for elements with borders.

    TESTS::

        sage: ascii_art(list(DyckWords(3)))
        [                                   /\   ]
        [            /\    /\      /\/\    /  \  ]
        [ /\/\/\, /\/  \, /  \/\, /    \, /    \ ]
    """
    new_mat = []
    for line in repr_elems:
        new_mat.append(l_border + line + " "**Integer(len(repr_elems) - len(line)) + r_border)
    new_bp = [bp + 2 for bp in repr_elems.get_breakpoints()] + [len(repr_elems) + 2]
    return AsciiArt(new_mat, new_bp, baseline = 0, atomic = False)

def ascii_art_list(list):
    r"""
    Return an ASCII art output of a list.

    TESTS::

        sage: ascii_art(list(DyckWords(3)))
        [                                   /\   ]
        [            /\    /\      /\/\    /  \  ]
        [ /\/\/\, /\/  \, /  \/\, /    \, /    \ ]
    """
    repr_elems = _ascii_art_iter(list)
    return _ascii_art_with_borders(repr_elems, "[ ", " ]")

def ascii_art_tuple(tuple):
    r"""
    Return an ASCII art output of a tuple.

    TESTS::

        sage: ascii_art(tuple(DyckWords(3)))
        (                                   /\   )
        (            /\    /\      /\/\    /  \  )
        ( /\/\/\, /\/  \, /  \/\, /    \, /    \ )
    """
    repr_elems = _ascii_art_iter(tuple)
    return _ascii_art_with_borders(repr_elems, "( ", " )")

empty_ascii_art = AsciiArt([""])

