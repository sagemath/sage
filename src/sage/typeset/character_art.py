# -*- coding: utf-8 -*-
r"""
Base Class for Character-Based Art

This is the common base class for
:class:`sage.typeset.ascii_art.AsciiArt` and
:class:`sage.typeset.unicode_art.UnicodeArt`. They implement simple
graphics by placing characters on a rectangular grid, in other words,
using monospace fonts. The difference is that one is restricted to
7-bit ascii, the other uses all unicode code points.
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
import sys
from sage.structure.sage_object import SageObject


################################################################################
# Global variable use to compute the maximal length allows for ascii art
# object.
MAX_WIDTH = None
################################################################################


class CharacterArt(SageObject):

    def __init__(self, lines=[], breakpoints=[], baseline=None):
        r"""
        Abstract base class for character art

        INPUT:

        - ``lines`` -- the list of lines of the representation of the
          character art object

        - ``breakpoints`` -- the list of points where the representation can be
          split

        - ``baseline`` -- the reference line (from the bottom)

        Instead of just integers, ``breakpoints`` may also contain tuples
        consisting of an offset and the breakpoints of a nested substring at
        that offset. This is used to prioritize the breakpoints, as line breaks
        inside the substring will be avoided if possible.

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

        If there are nested breakpoints, line breaks are avoided inside the
        nested elements (:trac:`29204`)::

            sage: s = ascii_art([[1..5], [1..17], [1..25]])
            sage: s._breakpoints
            [(2, [4, 7, 10, 13]), 20, (21, [4, 7,..., 56]), 83, (84, [4, 7,..., 88])]
            sage: str(s)
            '[ [ 1, 2, 3, 4, 5 ],\n\n
              [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17 ],\n\n
              [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,\n\n
              22, 23, 24, 25 ] ]'
        """
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
        yield from self._matrix

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
        # if the draw is larger than the max length, try to split it
        if hsize <= self._l and self._breakpoints:
            return self._split_repr_(hsize)
        return '\n'.join(self._matrix)

    def __format__(self, fmt):
        r"""
        Format ``self``.

        EXAMPLES::

            sage: M = matrix([[1,2],[3,4]])                                             # needs sage.modules
            sage: format(ascii_art(M))                                                  # needs sage.modules
            '[1 2]\n[3 4]'
            sage: format(unicode_art(M))                                                # needs sage.modules
            '\u239b1 2\u239e\n\u239d3 4\u23a0'
        """
        return format(self._string_type(self), fmt)

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

        This method is deprecated, as its output is an implementation detail.
        The mere breakpoints of a character art element do not reflect the best
        way to split it if nested structures are involved. For details, see
        :trac:`29204`.

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
            doctest:...: DeprecationWarning: get_breakpoints() is deprecated
            See https://github.com/sagemath/sage/issues/29204 for details.
            [6]
        """
        from sage.misc.superseded import deprecation
        deprecation(29204, "get_breakpoints() is deprecated")
        return self._breakpoints

    def _isatty(self):
        """
        Test whether ``stdout`` is a TTY.

        If this test succeeds, you can assume that ``stdout`` is directly
        connected to a terminal. Otherwise you should treat ``stdout`` as
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
            return sys.stdout.isatty()
        except Exception:
            # for fake ttys that might lead to an error
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
        import fcntl
        import termios
        import struct
        rc = fcntl.ioctl(int(0), termios.TIOCGWINSZ,
                         struct.pack('HHHH', sys.stdout.fileno(), 0, 0, 0))
        h, w, hp, wp = struct.unpack('HHHH', rc)
        return w

    def _splitting_points(self, size, offset=0):
        """
        Iterate over the breakpoints at which the representation can be split
        to obtain chunks of length at most ``size``.

        The final element returned will be the length of this object.

        INPUT:

        - ``size`` -- the maximum width of each chunk

        - ``offset`` -- (default: ``0``); the first chunk has width at most
          ``size - offset``

        TESTS::

            sage: list(ascii_art(*(['a'] * 90))._splitting_points(20, offset=5))
            [15, 35, 55, 75, 90]
        """
        # We implement a custom iterator instead of repeatedly using
        # itertools.chain to prepend elements in order to avoid quadratic time
        # complexity
        class PrependIterator():
            """
            Iterator with support for prepending of elements.
            """
            def __init__(self, stack):
                self._stack = [iter(elems) for elems in stack]

            def prepend(self, elems):
                self._stack.append(iter(elems))

            def __iter__(self):
                return self

            def __next__(self):
                while self._stack:
                    try:
                        return next(self._stack[-1])
                    except StopIteration:
                        self._stack.pop()
                raise StopIteration

        idx = -offset
        breakpoints = PrependIterator([[self._l], self._breakpoints])
        bp = None
        for bp_next in breakpoints:
            if not isinstance(bp_next, tuple):
                if bp_next - idx > size and bp is not None:
                    yield bp
                    idx = bp
                bp = bp_next
            else:
                sub_offset, sub_breakpoints = bp_next
                try:
                    bp_next = next(breakpoints)
                except StopIteration:
                    bp_next = None
                if bp_next is None or isinstance(bp_next, tuple):
                    raise ValueError("nested structure must be followed by a "
                                     "regular breakpoint")
                if bp_next - idx > size:
                    # substructure is too wide for the current line, so force a
                    # line break
                    if bp is not None:
                        yield bp
                        idx = bp
                    breakpoints.prepend([bp_next])
                    breakpoints.prepend(_shifted_breakpoints(sub_breakpoints,
                                                             sub_offset))
                    # at this point, we do not know the next breakpoint yet,
                    # but have already yielded bp, so discard it
                    bp = None
                else:
                    bp = bp_next
        if bp is not None:
            yield bp

    def _split_repr_(self, size):
        r"""
        Split the representation into chunks of length at most ``size``.

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

        ::

            sage: ascii_art(['a' * k for k in (1..10)])._split_repr_(20)
            '[ a, aa, aaa, aaaa,\n\n aaaaa, aaaaaa,\n\n aaaaaaa, aaaaaaaa,\n\n aaaaaaaaa,\n\n aaaaaaaaaa ]'

        Check that wrapping happens exactly at the given size (:trac:`28527`)::

            sage: len(ascii_art(*(['']*90), sep=',')._split_repr_(80).split('\n')[0])
            80
        """
        idx = 0
        parts = []
        for bp in self._splitting_points(size):
            if bp - idx > size:
                import warnings
                warnings.warn("the console size is smaller than the pretty "
                              "representation of the object")
            # Note that this is faster than calling self.split() repeatedly
            parts.append('\n'.join(line[idx:bp] for line in self))
            idx = bp
        return '\n\n'.join(parts)

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
        left = []
        right = []
        for line in self:
            left.append(line[:pos])
            right.append(line[pos:])
        l_bp = []
        r_bp = []
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
        """
        return max(obj1.get_baseline(), obj2.get_baseline())

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
        Concatenate two character art objects.

        By default, when two objects are concatenated, the new one will be
        splittable between both.

        If the baseline is defined, the concatenation is computed such that the
        new baseline coincides with the olders.

        For example, let `T` be a tree with its baseline ascii art
        representation in the middle::

            o
             \
              o
             / \
            o   o

        and let `M` be a matrix with its baseline ascii art representation at
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

        Since :trac:`28527`, the baseline must always be a number.

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
                        new_matrix.insert(0, " " * self._l)

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
                i = max(new_h - new_baseline - Nelt._h + Nelt._baseline, 0)
            for j in range(Nelt._h):
                new_matrix[i + j] += Nelt._matrix[j]
        else:
            for line in self._matrix:
                new_matrix.append(line + " " * (self._l - len(line)))
            for i, line_i in enumerate(Nelt._matrix):
                if i == len(new_matrix):
                    new_matrix.append(" " * self._l + line_i)
                else:
                    new_matrix[i] += line_i

        # breakpoint
        new_breakpoints = list(self._breakpoints)
        if self._l and Nelt._l:
            new_breakpoints.append(self._l)
        new_breakpoints.extend(_shifted_breakpoints(Nelt._breakpoints,
                                                    self._l))
        return self.__class__(
            lines=new_matrix,
            breakpoints=new_breakpoints,
            baseline=new_baseline,
        )

    def __mul__(self, Nelt):
        r"""
        Stack two character art objects vertically.

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


def _shifted_breakpoints(breakpoints, offset):
    """
    Return an iterator that shifts all breakpoints by an offset.

    TESTS::

        sage: b = ascii_art([[1, 2, 3], [4, 5, 6]])._breakpoints; b
        [(2, [4, 7]), 14, (15, [4, 7])]
        sage: from sage.typeset.character_art import _shifted_breakpoints
        sage: list(_shifted_breakpoints(b, 10))
        [(12, [4, 7]), 24, (25, [4, 7])]
    """
    for bp in breakpoints:
        if isinstance(bp, tuple):
            yield bp[0] + offset, bp[1]
        else:
            yield bp + offset
