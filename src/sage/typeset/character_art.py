# -*- coding: utf-8 -*-
r"""
Base Class for Character-Based Art

This is the common base class for
:class:`sage.typeset.ascii_art.AsciiArt` and
:class:`sage.typeset.ascii_art.UnicodeArt`. They implement simple
graphics by placing characters on a rectangular grid, in other words,
using monospace fonts. The difference is that one is restricted to
7-bit ascii, the other uses all unicode code points.
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
from __future__ import print_function

import os, sys
from sage.structure.sage_object import SageObject


################################################################################
### Global variable use to compute the maximal length allows for ascii art
### object.
MAX_WIDTH = None
################################################################################

class CharacterArt(SageObject):

    def __init__(self, lines=[], breakpoints=[], baseline=None, atomic=None):
        r"""
        Abstract base class for character art

        INPUT:

        - ``lines`` -- the list of lines of the representation of the
          character art object

        - ``breakpoints`` -- the list of points where the representation can be
          split

        - ``baseline`` -- the reference line (from the bottom)

        EXAMPLES::

            sage: i = var('i')
            sage: ascii_art(sum(pi^i/factorial(i)*x^i, i, 0, oo))
             pi*x
            e

        TESTS::

            sage: from sage.typeset.ascii_art import AsciiArt
            sage: aao = AsciiArt()
            sage: aao
            <BLANKLINE>
            sage: aa = AsciiArt(["  *  ", " * * ", "*****"]); aa
              *
             * *
            *****
        """
        if atomic is not None:
            from sage.misc.superseded import deprecation
            deprecation(18357, "the argument atomic is deprecated and will be ignored")
        self._matrix = lines
        self._breakpoints = breakpoints
        self._baseline = baseline if baseline is not None else 0

        self._h = len(lines)
        self._l = 0 if not lines else max([len(line) for line in lines])

    @classmethod
    def empty(cls):
        """
        Return the empty character art object

        EXAMPLES::

            sage: from sage.typeset.ascii_art import AsciiArt
            sage: AsciiArt.empty()
        """
        empty_string = cls._string_type()
        return cls([empty_string])

    def __getitem__(self, key):
        r"""
        Return the line `key` of the ASCII art object.

        TESTS::

            sage: from sage.typeset.ascii_art import AsciiArt
            sage: p5 = AsciiArt(["  *  ", " * * ", "*****"])
            sage: p5[1]
            ' * * '
        """
        return self._matrix[key]

    def __iter__(self):
        r"""
        Iterator on all lines of the ASCII art object.

        TESTS::

            sage: from sage.typeset.ascii_art import AsciiArt
            sage: p5 = AsciiArt(["  *  ", " * * ", "*****"])
            sage: for line in p5:
            ....:     print(line)
              *
             * *
            *****
        """
        for elem in self._matrix:
            yield elem

    def _repr_(self):
        r"""
        TESTS::

            sage: from sage.typeset.ascii_art import AsciiArt
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
        if hsize <= self._l and self._breakpoints:
            return self._split_repr_(hsize)
        #########
        output = ""
        if len(self._matrix) > 0:
            for i in range(len(self._matrix) - 1):
                output += self._matrix[i] + "\n"
            return output + self._matrix[len(self._matrix) - 1]
        return output

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

            sage: from sage.typeset.ascii_art import AsciiArt
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

            sage: from sage.typeset.ascii_art import AsciiArt
            sage: p3 = AsciiArt([" * ", "***"])
            sage: p5 = AsciiArt(["  *  ", " * * ", "*****"])
            sage: aa = ascii_art([p3, p5])
            sage: aa.get_breakpoints()
            [6]
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

            sage: from sage.typeset.ascii_art import empty_ascii_art
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

            sage: from sage.typeset.ascii_art import empty_ascii_art
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

            sage: from sage.typeset.ascii_art import AsciiArt
            sage: p3 = AsciiArt([" * ", "***"])
            sage: p5 = AsciiArt(["  *  ", " * * ", "*****"])
            sage: aa = ascii_art([p3, p5])
            sage: print(aa._split_repr_(10))
            [
            [  *
            [ ***,
            <BLANKLINE>
                *   ]
               * *  ]
              ***** ]
        """
        f_split = self._breakpoints[0]; i = 1
        while i < len(self._breakpoints) and self._breakpoints[i] < size:
            f_split = self._breakpoints[i]
            i += 1
        if size <= f_split:
            import warnings
            warnings.warn("the console size is smaller than the pretty" +
                "representation of the object")
        top, bottom = self.split(f_split)
        return repr(top * self.empty()) + "\n" + repr(bottom)

    def split(self, pos):
        r"""
        Split the representation at the position ``pos``.

        EXAMPLES::

            sage: from sage.typeset.ascii_art import AsciiArt
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
        return self.__class__(left, l_bp), self.__class__(right, r_bp)

    @staticmethod
    def _compute_new_baseline(obj1, obj2):
        r"""
        TESTS::

            sage: from sage.typeset.ascii_art import AsciiArt
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

            sage: from sage.typeset.ascii_art import AsciiArt
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

    def width(self):
        r"""
        Return the length (width) of the ASCII art object.

        OUTPUT:

        Integer. The number of characters in each line.

        TESTS::

            sage: from sage.typeset.ascii_art import AsciiArt
            sage: p3 = AsciiArt([" * ", "***"])
            sage: len(p3), p3.width(), p3.height()
            (3, 3, 2)
            sage: p5 = AsciiArt(["  *  ", " * * ", "*****"])
            sage: len(p5), p5.width(), p5.height()
            (5, 5, 3)
        """
        return self._l

    __len__ = width

    def height(self):
        r"""
        Return the height of the ASCII art object.

        OUTPUT:

        Integer. The number of lines.

        TESTS::

            sage: from sage.typeset.ascii_art import AsciiArt
            sage: p3 = AsciiArt([" * ", "***"])
            sage: p3.height()
            2
            sage: p5 = AsciiArt(["  *  ", " * * ", "*****"])
            sage: p5.height()
            3
        """
        return self._h

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

            sage: from sage.typeset.ascii_art import AsciiArt
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
        new_h = self.__class__._compute_new_h(self, Nelt)
        new_baseline = self.__class__._compute_new_baseline(self, Nelt)

        if self._baseline is not None and Nelt._baseline is not None:
            # left treatement
            for line in self._matrix:
                new_matrix.append(line + " " * (self._l - len(line)))

            if new_h > self._h:
                # |                 new_h > self._h
                # |                 new_baseline > self._baseline
                # ||<-- baseline    number of white lines at the bottom
                #  | }               :: Nelt._baseline - self._baseline
                #  | }
                if new_baseline > self._baseline:
                    for k in range(new_baseline - self._baseline):
                        new_matrix.append(" " * self._l)
                #  | }              new_h > self._h
                #  | }              new_h - new_baseline > self._h - self._baseline
                # ||<-- baseline    number of white lines at the top
                # ||                :: new_h - new_baseline - self._h + self._baseline
                # ||
                #  |
                #  |
                if new_h - new_baseline > self._h - self._baseline:
                    for _ in range((new_h - new_baseline) - (self._h - self._baseline)):
                        new_matrix.insert(0, " "  * self._l)

            # right treatement
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
                new_matrix.append(line + " " * (self._l - len(line)))
            for i, line_i in enumerate(Nelt._matrix):
                if i == len(new_matrix):
                    new_matrix.append(" " * self._l + line_i)
                else: new_matrix[i] += line_i

        # breakpoint
        new_breakpoints = list(self._breakpoints)
        new_breakpoints.append(self._l)
        for bp in Nelt._breakpoints:
            new_breakpoints.append(bp + self._l)
        from sage.misc.misc import uniq
        return self.__class__(
            lines=new_matrix,
            breakpoints=uniq(new_breakpoints),
            baseline=new_baseline,
        )

    def __mul__(self, Nelt):
        r"""
        The operator ``*`` is use to the representation ``self`` at
        the top of an other ``Nelt``.

        TESTS::

            sage: from sage.typeset.ascii_art import AsciiArt
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
        new_repr = self.__class__(self._matrix + Nelt._matrix)
        return new_repr
