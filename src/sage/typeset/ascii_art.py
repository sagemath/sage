# -*- coding: utf-8 -*-
r"""
ASCII Art

This file contains:

- :class:`AsciiArt` an simple implementation of an ASCII art object,
- :func:`ascii_art` a function to get the ASCII art representation of any
  object in Sage,


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
# ******************************************************************************
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
#                  https://www.gnu.org/licenses/
# ******************************************************************************

from sage.typeset.character_art import CharacterArt
from sage.typeset.character_art_factory import CharacterArtFactory
import sage.typeset.symbols as symbol


class AsciiArt(CharacterArt):
    r"""
    An Ascii art object is an object with some specific representation for
    *printing*.

    INPUT:

    - ``lines`` -- the list of lines of the representation of the ascii art
      object

    - ``breakpoints`` -- the list of points where the representation can be
      split

    - ``baseline`` -- the reference line (from the bottom)

    EXAMPLES::

        sage: i = var('i')
        sage: ascii_art(sum(pi^i/factorial(i)*x^i, i, 0, oo))
         pi*x
        e
    """
    _string_type = str


_ascii_art_factory = CharacterArtFactory(
    AsciiArt, str, '_ascii_art_',
    (symbol.ascii_left_parenthesis, symbol.ascii_right_parenthesis),
    (symbol.ascii_left_square_bracket, symbol.ascii_right_square_bracket),
    (symbol.ascii_left_curly_brace, symbol.ascii_right_curly_brace),
)


empty_ascii_art = _ascii_art_factory.build_empty()


def ascii_art(*obj, **kwds):
    r"""
    Return an ASCII art representation

    INPUT:

    - ``*obj`` -- any number of positional arguments, of arbitrary
      type. The objects whose ascii art representation we want.

    - ``sep`` -- optional ``'sep=...'`` keyword argument (or ``'separator'``).
      Anything that can be converted to ascii art (default: empty ascii
      art). The separator in-between a list of objects. Only used if
      more than one object given.

    - ``baseline`` -- (default: 0) the baseline for the object

    - ``sep_baseline`` -- (default: 0) the baseline for the separator

    OUTPUT:

    :class:`AsciiArt` instance.

    EXAMPLES::

        sage: result = ascii_art(integral(exp(x+x^2)/(x+1), x))
        ...
        sage: result
            /
           |
           |   2
           |  x  + x
           | e
           | ------- dx
           |  x + 1
           |
          /

    We can specify a separator object::

        sage: ident = lambda n: identity_matrix(ZZ, n)
        sage: ascii_art(ident(1), ident(2), ident(3), sep=' : ')
                      [1 0 0]
              [1 0]   [0 1 0]
        [1] : [0 1] : [0 0 1]

    We can specify the baseline::

        sage: ascii_art(ident(2), baseline=-1) + ascii_art(ident(3))
        [1 0][1 0 0]
        [0 1][0 1 0]
             [0 0 1]

    We can determine the baseline of the separator::

        sage: ascii_art(ident(1), ident(2), ident(3), sep=' -- ', sep_baseline=-1)
                        [1 0 0]
            -- [1 0] -- [0 1 0]
        [1]    [0 1]    [0 0 1]

    If specified, the ``sep_baseline`` overrides the baseline of
    an ascii art separator::

        sage: sep_line = ascii_art('\n'.join(' | ' for _ in range(6)), baseline=6)
        sage: ascii_art(*Partitions(6), separator=sep_line, sep_baseline=0)
               |       |      |      |     |     |     |    |    |    | *
               |       |      |      |     |     |     |    |    | ** | *
               |       |      |      |     |     | *** |    | ** | *  | *
               |       |      | **** |     | *** | *   | ** | ** | *  | *
               | ***** | **** | *    | *** | **  | *   | ** | *  | *  | *
        ****** | *     | **   | *    | *** | *   | *   | ** | *  | *  | *

    TESTS::

        sage: n = var('n')
        sage: ascii_art(sum(binomial(2 * n, n + 1) * x^n, n, 0, oo))
         /        _________    \
        -\2*x + \/ 1 - 4*x  - 1/
        -------------------------
                   _________
             2*x*\/ 1 - 4*x
        sage: ascii_art(list(DyckWords(3)))
        [                                   /\   ]
        [            /\    /\      /\/\    /  \  ]
        [ /\/\/\, /\/  \, /  \/\, /    \, /    \ ]
        sage: ascii_art(1)
        1
    """
    separator, baseline, sep_baseline = _ascii_art_factory.parse_keywords(kwds)
    if kwds:
        raise ValueError('unknown keyword arguments: {0}'.format(list(kwds)))
    if len(obj) == 1:
        return _ascii_art_factory.build(obj[0], baseline=baseline)
    if not isinstance(separator, AsciiArt):
        separator = _ascii_art_factory.build(separator, baseline=sep_baseline)
    elif sep_baseline is not None:
        from copy import copy
        separator = copy(separator)
        separator._baseline = sep_baseline
    return _ascii_art_factory.concatenate(obj, separator, empty_ascii_art,
                                          baseline=baseline)
