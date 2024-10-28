"""
Graftal Lace Cellular Automata

AUTHORS:

- Travis Scrimshaw (2020-04-30): Initial version
"""

#*****************************************************************************
#       Copyright (C) 2020 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.sage_object import SageObject
from sage.typeset.ascii_art import AsciiArt
from sage.typeset.unicode_art import UnicodeArt

class GraftalLaceCellularAutomata(SageObject):
    r"""
    Graftal Lace Cellular Automata (GLCA).

    A GLCA is a deterministic cellular automaton whose rule is given
    by an 8-digit octal number `r_7 \cdots r_0`. For a node `s_i`, let
    `b_k`, for `k = -1,0,1` denote if there is an edge from `s_i` to
    `s'_{i+k}`, where `s'_j` is the previous row. We determine the value
    at `t_{i+k}` by considering the value of `r_m`, where the binary
    representation of `m` is `b_{-1} b_0 b_1`. If `r_m` has a binary
    representation of b'_1 b'_0 b'_{-1}`, then we add `b'_k` to `t_{i+k}`.

    INPUT:

    - ``rule`` -- list of length 8 with integer entries `0 \leq x < 8`

    EXAMPLES::

        sage: G = cellular_automata.GraftalLace([0,2,5,4,7,2,3,3])
        sage: G.evolve(3)
        sage: ascii_art(G)
               o
               |
               o
               |
             o o o
            /| |/|
           o o o o o
          /| |/|\|/|
         o o o o o o o

        sage: G = cellular_automata.GraftalLace([3,0,3,4,7,6,3,1])
        sage: G.evolve(3)
        sage: ascii_art(G)
               o
               |
               o
               |\
             o o o
            /  |\ \
           o o o o o
          /|/  |\ \ \
         o o o o o o o

        sage: G = cellular_automata.GraftalLace([2,0,3,3,6,0,2,7])
        sage: G.evolve(20)
        sage: G.plot()                                                                  # needs sage.plot
        Graphics object consisting of 842 graphics primitives

    .. PLOT::

        G = cellular_automata.GraftalLace([2,0,3,3,6,0,2,7])
        G.evolve(20)
        P = G.plot()
        sphinx_plot(P)

    REFERENCES:

    - [Kas2018]_
    """
    def __init__(self, rule):
        """
        Initialize ``self``.

        TESTS::

            sage: G = cellular_automata.GraftalLace([5,1,2,5,4,5,5,0])
            sage: TestSuite(G).run()

            sage: G = cellular_automata.GraftalLace([5,1,2,5,4,5,5])
            Traceback (most recent call last):
            ...
            ValueError: invalid rule
            sage: G = cellular_automata.GraftalLace([5,1,2,5,4,5,5,0,5])
            Traceback (most recent call last):
            ...
            ValueError: invalid rule
            sage: G = cellular_automata.GraftalLace([5,1,2,5,-1,5,5,0])
            Traceback (most recent call last):
            ...
            ValueError: invalid rule
            sage: G = cellular_automata.GraftalLace([8,5,1,2,5,4,5,5])
            Traceback (most recent call last):
            ...
            ValueError: invalid rule
        """
        if len(rule) != 8 or any(x not in range(8) for x in rule):
            raise ValueError("invalid rule")
        # We reverse the rule to make the evolution easier to compute
        self._rule = tuple(reversed(rule))
        self._states = [[2]]

    def __eq__(self, other):
        """
        Check equality.

        Two GLCAs are equal if and only if they have the same rule.

        TESTS::

            sage: G1 = cellular_automata.GraftalLace([5,1,2,5,4,5,5,0])
            sage: G2 = cellular_automata.GraftalLace([0,5,5,4,5,2,1,5])
            sage: G3 = cellular_automata.GraftalLace([5,1,2,5,4,5,5,0])
            sage: G1 == G2
            False
            sage: G1 == G3
            True
            sage: G1 is G3
            False
        """
        return (isinstance(other, GraftalLaceCellularAutomata)
                and self._rule == other._rule)

    def __ne__(self, other):
        """
        Check non equality.

        TESTS::

            sage: G1 = cellular_automata.GraftalLace([5,1,2,5,4,5,5,0])
            sage: G2 = cellular_automata.GraftalLace([0,5,5,4,5,2,1,5])
            sage: G3 = cellular_automata.GraftalLace([5,1,2,5,4,5,5,0])
            sage: G1 != G2
            True
            sage: G1 != G3
            False
            sage: G1 is G3
            False
        """
        return not (self == other)

    # Evolution functions
    # -------------------

    def evolve(self, number=None):
        r"""
        Evolve ``self``.

        INPUT:

        - ``number`` -- (default: 1) the number of times to perform
          the evolution

        EXAMPLES::

            sage: G = cellular_automata.GraftalLace([5,1,2,5,4,5,5,0])
            sage: ascii_art(G)
             o
             |
             o
            sage: G.evolve(2)
            sage: ascii_art(G)
                 o
                 |
                 o
                / \
               o o o
              / \ / \
             o o o o o

            sage: G = cellular_automata.GraftalLace([0,2,1,4,7,2,3,0])
            sage: G.evolve(3)
            sage: ascii_art(G)
                   o
                   |
                   o
                   |
                 o o o
                   |
               o o o o o
                   |
             o o o o o o o
        """
        if number is not None:
            for _ in range(number):
                self.evolve()
            return

        prev_state = self._states[-1]
        next_state = [0] * (len(prev_state) + 2)

        for i, val in enumerate(prev_state):
            next_state[i] += self._rule[val] & 0x1
            next_state[i+1] += self._rule[val] & 0x2
            next_state[i+2] += self._rule[val] & 0x4
        self._states.append(next_state)

    # Output functions
    # ----------------

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: cellular_automata.GraftalLace([5,1,2,5,4,5,5,0])
            Graftal Lace Cellular Automata with rule 51254550
        """
        rule = ''.join(map(str, reversed(self._rule)))
        return "Graftal Lace Cellular Automata with rule {}".format(rule)

    def print_states(self, number=None, use_unicode=False):
        r"""
        Print the first ``num`` states of ``self``.

        .. NOTE::

            If the number of states computed for ``self`` is less than
            ``num``, then this evolves the system using the default
            time evolution.

        INPUT:

        - ``number`` -- the number of states to print

        EXAMPLES::

            sage: G = cellular_automata.GraftalLace([5,1,2,5,4,5,5,0])
            sage: G.evolve(2)
            sage: G.print_states()
                 o
                 |
                 o
                / \
               o o o
              / \ / \
             o o o o o
            sage: G.evolve(20)
            sage: G.print_states(3)
                 o
                 |
                 o
                / \
               o o o
              / \ / \
             o o o o o
        """
        if number is None:
            number = len(self._states)
        if number > len(self._states):
            for dummy in range(number - len(self._states)):
                self.evolve()
        if use_unicode:
            print(self._unicode_art_(number))
        else:
            print(self._ascii_art_(number))

    def _ascii_art_(self, number=None):
        r"""
        Return an ascii art representation of ``self``.

        EXAMPLES::

            sage: G = cellular_automata.GraftalLace([5,1,2,5,4,5,5,0])
            sage: G.evolve(10)
            sage: ascii_art(G)
                                 o
                                 |
                                 o
                                / \
                               o o o
                              / \ / \
                             o o o o o
                            / \  |  / \
                           o o o o o o o
                          / \ / X X \ / \
                         o o o o o o o o o
                        / \  |/ \|/ \|  / \
                       o o o o o o o o o o o
                      / \ / \ \ / \ / / \ / \
                     o o o o o o o o o o o o o
                    / \  |  / \|   |/ \  |  / \
                   o o o o o o o o o o o o o o o
                  / \ / X X \ /     \ / X X \ / \
                 o o o o o o o o o o o o o o o o o
                / \  |/ \|/ \|       |/ \|/ \|  / \
               o o o o o o o o o o o o o o o o o o o
              / \ / \ \ / \ /         \ / \ / / \ / \
             o o o o o o o o o o o o o o o o o o o o o
        """
        if number is None:
            number = len(self._states)

        space = len(self._states[:number]) * 2 - 1
        ret = AsciiArt([' '*space + 'o'])
        space += 1
        for i,state in enumerate(self._states[:number]):
            temp = ' '*(space-2)
            last = ' '
            for x in state:
                if x & 0x4:
                    if last == '/':
                        temp += 'X'
                    else:
                        temp += '\\'
                else:
                    temp += last
                temp += '|' if x & 0x2 else ' '
                last = '/' if x & 0x1 else ' '
            ret *= AsciiArt([temp + last])
            space -= 1
            ret *= AsciiArt([' '*space + ' '.join('o' for dummy in range(2*i+1))])
            space -= 1
        return ret

    def _unicode_art_(self, number=None):
        r"""
        Return a unicode art representation of ``self``.

        EXAMPLES::

            sage: G = cellular_automata.GraftalLace([5,1,2,5,4,5,5,0])
            sage: G.evolve(10)
            sage: unicode_art(G)
                                 ◾
                                 │
                                 ◾
                                ╱ ╲
                               ◾ ◾ ◾
                              ╱ ╲ ╱ ╲
                             ◾ ◾ ◾ ◾ ◾
                            ╱ ╲  │  ╱ ╲
                           ◾ ◾ ◾ ◾ ◾ ◾ ◾
                          ╱ ╲ ╱ ╳ ╳ ╲ ╱ ╲
                         ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾
                        ╱ ╲  │╱ ╲│╱ ╲│  ╱ ╲
                       ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾
                      ╱ ╲ ╱ ╲ ╲ ╱ ╲ ╱ ╱ ╲ ╱ ╲
                     ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾
                    ╱ ╲  │  ╱ ╲│   │╱ ╲  │  ╱ ╲
                   ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾
                  ╱ ╲ ╱ ╳ ╳ ╲ ╱     ╲ ╱ ╳ ╳ ╲ ╱ ╲
                 ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾
                ╱ ╲  │╱ ╲│╱ ╲│       │╱ ╲│╱ ╲│  ╱ ╲
               ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾
              ╱ ╲ ╱ ╲ ╲ ╱ ╲ ╱         ╲ ╱ ╲ ╱ ╱ ╲ ╱ ╲
             ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾ ◾
        """
        if number is None:
            number = len(self._states)

        space = len(self._states[:number]) * 2 - 1
        ret = UnicodeArt([' '*space + '◾'])
        space += 1
        for i,state in enumerate(self._states[:number]):
            temp = ' '*(space-2)
            last = ' '
            for x in state:
                if x & 0x4:
                    if last == '╱':
                        temp += '╳'
                    else:
                        temp += '╲'
                else:
                    temp += last
                temp += '│' if x & 0x2 else ' '
                last = '╱' if x & 0x1 else ' '
            ret *= UnicodeArt([temp + last])
            space -= 1
            ret *= UnicodeArt([' '*space + ' '.join('◾' for dummy in range(2*i+1))])
            space -= 1
        return ret

    def plot(self, number=None):
        r"""
        Return a plot of ``self``.

        INPUT:

        - ``number`` -- the number of states to plot

        EXAMPLES::

            sage: G = cellular_automata.GraftalLace([5,1,2,5,4,5,5,0])
            sage: G.evolve(20)
            sage: G.plot()                                                              # needs sage.plot
            Graphics object consisting of 865 graphics primitives
        """
        if number is None:
            number = len(self._states)
        if number > len(self._states):
            for dummy in range(number - len(self._states)):
                self.evolve()

        from sage.plot.circle import circle
        from sage.plot.line import line

        x = len(self._states[:number])
        rad = 0.1
        ret = circle((x, 1), rad, fill=True)
        for i,state in enumerate(self._states[:number]):
            for j,val in enumerate(state):
                if val & 0x4:
                    ret += line([(x+j,-i), (x+j-1,-i+1)])
                if val & 0x2:
                    ret += line([(x+j,-i), (x+j,-i+1)])
                if val & 0x1:
                    ret += line([(x+j,-i), (x+j+1,-i+1)])
            for j in range(2*i+1):
                ret += circle((x+j, -i), rad, fill=True)
            x -= 1
        ret.set_aspect_ratio(1)
        ret.axes(False)
        return ret

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: G = cellular_automata.GraftalLace([5,1,2,5,4,5,5,0])
            sage: G.evolve(2)
            sage: latex(G)
            \begin{tikzpicture}
            \fill (3,1) circle (2pt);
            \draw[-] (3,0) -- (3,1);
            \fill (3,0) circle (2pt);
            \draw[-] (2,-1) -- (3,0);
            \draw[-] (4,-1) -- (3,0);
            \fill (2,-1) circle (2pt);
            \fill (3,-1) circle (2pt);
            \fill (4,-1) circle (2pt);
            \draw[-] (1,-2) -- (2,-1);
            \draw[-] (3,-2) -- (2,-1);
            \draw[-] (3,-2) -- (4,-1);
            \draw[-] (5,-2) -- (4,-1);
            \fill (1,-2) circle (2pt);
            \fill (2,-2) circle (2pt);
            \fill (3,-2) circle (2pt);
            \fill (4,-2) circle (2pt);
            \fill (5,-2) circle (2pt);
            \end{tikzpicture}
        """
        ret = "\\begin{tikzpicture}\n"

        x = len(self._states)
        rad = 2
        ret += "\\fill ({},{}) circle ({}pt);\n".format(x, 1, rad)
        for i,state in enumerate(self._states):
            for j,val in enumerate(state):
                if val & 0x4:
                    ret += "\\draw[-] ({},{}) -- ({},{});\n".format(x+j,-i, x+j-1,-i+1)
                if val & 0x2:
                    ret += "\\draw[-] ({},{}) -- ({},{});\n".format(x+j,-i, x+j,-i+1)
                if val & 0x1:
                    ret += "\\draw[-] ({},{}) -- ({},{});\n".format(x+j,-i, x+j+1,-i+1)
            for j in range(2*i+1):
                ret += "\\fill ({},{}) circle ({}pt);\n".format(x+j, -i, rad)
            x -= 1
        ret += "\\end{tikzpicture}"
        return ret
