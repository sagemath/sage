"""
Elementary Cellular Automata

AUTHORS:

- Travis Scrimshaw (2018-07-07): Initial version
"""

#*****************************************************************************
#       Copyright (C) 2018 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object import SageObject
from sage.typeset.ascii_art import AsciiArt
from sage.typeset.unicode_art import UnicodeArt
from sage.rings.integer_ring import ZZ
from sage.matrix.constructor import matrix
from sage.misc.lazy_import import lazy_import
lazy_import("sage.plot.matrix_plot", "matrix_plot")
from sage.misc.constant_function import ConstantFunction

class ElementaryCellularAutomata(SageObject):
    r"""
    Elementary cellular automata.

    An *elementary cellular automaton* is a 1-dimensional cellular
    deterministic automaton with two possible values: `X := \{0,1\}`.
    A *state* is therefore a sequence `s \in X^n`, and the *evolution*
    of a state `s \to s'` is given for `s'_i` by looking at the values
    at positions `s_{i-1}, s_i, s_{i+1}` and is determined by the
    *rule* `0 \leq r \leq 255` as follows. Consider the binary
    representation `r = b_7 b_6 b_5 b_4 b_3 b_2 b_1 b_0`. Then, we
    define `s'_i = b_j`, where `j = s_{i-1} s_i s_{i+1}` is the
    corresponding binary representation. In other words, the value
    `s'_i` is given according to the following table:

    .. MATH::

        \begin{array}{cccccccc}
        111 & 110 & 101 & 100 & 011 & 010 & 001 & 000 \\
        b_7 & b_6 & b_5 & b_4 & b_3 & b_2 & b_1 & b_0
        \end{array}

    We consider the boundary values of `s_0 = s_{n+1} = 0`.

    INPUT:

    - ``rule`` -- integer between 0 and 255
    - ``width`` -- (optional) the width of the ECA
    - ``initial_state`` -- (optional) the initial state given
      as a list of ``0`` and ``1``
    - ``boundary`` -- (default: ``(0, 0)``) a tuple of the left and right
      boundary conditions respectively or ``None`` for periodic boundary
      conditions

    Either ``width`` or ``initial_state`` must be given. If ``width``
    is less than the length of ``initial_state``, then ``initial_state``
    has ``0`` prepended so the resulting list has length ``width``.
    If only ``width`` is given, then the initial state is constructed
    randomly.

    The boundary conditions can either be ``0``, ``1``, or a function that
    takes an integer ``n`` corresponding to the state and outputs either
    ``0`` or ``1``.

    EXAMPLES:

    We construct an example with rule `r = 90` using `n = 20`. The
    initial state consists of a single `1` in the rightmost entry::

        sage: ECA = cellular_automata.Elementary(90, width=20, initial_state=[1])
        sage: ECA.evolve(20)
        sage: ascii_art(ECA)
                           X
                          X
                         X X
                        X
                       X X
                      X   X
                     X X X X
                    X
                   X X
                  X   X
                 X X X X
                X       X
               X X     X X
              X   X   X   X
             X X X X X X X X
            X
           X X
          X   X
         X X X X
        X       X
         X     X X

    We now construct it with different boundary conditions. The first is
    with the left boundary being `1` (instead of `0`)::

        sage: ECA = cellular_automata.Elementary(90, width=20, initial_state=[1], boundary=(1,0))
        sage: ECA.evolve(20)
        sage: ascii_art(ECA)
                           X
        X                 X
        XX               X X
         XX             X
         XXX           X X
         X XX         X   X
           XXX       X X X X
        X XX XX     X
        X XX XXX   X X
        X XX X XX X   X
        X XX   XX  X X X
        X XXX XXXXX     X
        X X X X   XX   X X
        X      X XXXX X   X
        XX    X  X  X  X X X
         XX  X XX XX XX
         XXXX  XX XX XXX
         X  XXXXX XX X XX
          XXX   X XX   XXX
        XXX XX X  XXX XX XX
          X XX  XXX X XX XXX

    Now we consider the right boundary as being `1` on every third value::

        sage: def rbdry(n): return 1 if n % 3 == 0 else 0
        sage: ECA = cellular_automata.Elementary(90, width=20, initial_state=[1], boundary=(0,rbdry))
        sage: ECA.evolve(20)
        sage: ascii_art(ECA)
                           X
                          X
                         X X
                        X  X
                       X XX
                      X  XXX
                     X XXX
                    X  X XX
                   X XX  XXX
                  X  XXXXX
                 X XXX   XX
                X  X XX XXXX
               X XX  XX X
              X  XXXXXX  X
             X XXX    XXX X
            X  X XX  XX X
           X XX  XXXXXX  X
          X  XXXXX    XXX X
         X XXX   XX  XX X
        X  X XX XXXXXXX  X
         XX  XX X     XXX X

    Lastly we consider it with periodic boundary condition::

        sage: ECA = cellular_automata.Elementary(90, width=20, initial_state=[1], boundary=None)
        sage: ECA.evolve(20)
        sage: ascii_art(ECA)
                           X
        X                 X
         X               X
        X X             X X
           X           X
          X X         X X
         X   X       X   X
        X X X X     X X X X
               X   X
              X X X X
             X       X
            X X     X X
           X   X   X   X
          X X X X X X X X
         X               X
        X X             X X
           X           X
          X X         X X
         X   X       X   X
        X X X X     X X X X
               X   X

    We show the local evolution rules for rule `110`::

        sage: for t in cartesian_product([[0,1],[0,1],[0,1]]):
        ....:     ECA = cellular_automata.Elementary(110, list(t))
        ....:     ECA.print_states(2)
        ....:     print('#')
        <BLANKLINE>
        <BLANKLINE>
        #
          X
         XX
        #
         X
        XX
        #
         XX
        XXX
        #
        X
        X
        #
        X X
        XXX
        #
        XX
        XX
        #
        XXX
        X X
        #

    We construct an elementary cellular automaton with a random initial
    state with `n = 15` and see the state after `50` evolutions::

        sage: ECA = cellular_automata.Elementary(26, width=25)
        sage: ECA.print_state(50)  # random
            X X   X   X X X

    We construct and plot a larger example with rule `60`::

        sage: ECA = cellular_automata.Elementary(60, width=200)
        sage: ECA.evolve(200)
        sage: ECA.plot()                                                                # needs sage.plot
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        ECA = cellular_automata.Elementary(60, width=200)
        ECA.evolve(200)
        P = ECA.plot()
        sphinx_plot(P)

    With periodic boundary condition for rule `90`::

        sage: ECA = cellular_automata.Elementary(90, initial_state=[1]+[0]*254+[1], boundary=None)
        sage: ECA.evolve(256)
        sage: ECA.plot()                                                                # needs sage.plot
        Graphics object consisting of 1 graphics primitive

    .. PLOT::

        ECA = cellular_automata.Elementary(90, initial_state=[1]+[0]*254+[1], boundary=None)
        ECA.evolve(256)
        P = ECA.plot()
        sphinx_plot(P)

    REFERENCES:

    :wikipedia:`Elementary_cellular_automaton`
    """
    def __init__(self, rule, width=None, initial_state=None, boundary=(0, 0)):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: ECA = cellular_automata.Elementary(110, width=300)
            sage: TestSuite(ECA).run()
        """
        if rule not in ZZ or rule < 0 or rule > 255:
            raise ValueError("invalid rule")
        self._rule = ZZ(rule).binary()
        # We reverse the rule to make it easier to work with
        self._rule = [ZZ(x) for x in reversed('0'*(8-len(self._rule)) + self._rule)]
        if isinstance(width, list):
            initial_state = width
            width = len(initial_state)
        if initial_state is None:
            self._width = width
            initial_state = [ZZ.random_element(0,2) for d in range(width)]
        else:
            if not all(d in [0,1] for d in initial_state):
                raise ValueError("invalid initial state")
            initial_state = list(initial_state) # make sure it is a list and a copy
            if width is None:
                self._width = len(initial_state)
            elif width >= len(initial_state):
                self._width = width
                initial_state = ([0]*(width - len(initial_state))
                                 + initial_state)
            else:
                raise ValueError("the width must be at least the length of"
                                 " the initial state")
        self._states = [initial_state]
        if boundary is not None:
            self._bdry = tuple(boundary)
            self._lbdry = ConstantFunction(boundary[0]) if boundary[0] in [0,1] else boundary[0]
            self._rbdry = ConstantFunction(boundary[1]) if boundary[1] in [0,1] else boundary[1]
        else:
            self._bdry = boundary

    def __eq__(self, other):
        """
        Check equality.

        Two ECAs are equal when they have the same rule, width,
        initial state, and boundary conditions.

        TESTS::

            sage: ECA1 = cellular_automata.Elementary(110, [1,0,0,1,1,0,0,1,0,1])
            sage: ECA2 = cellular_automata.Elementary(101, [1,0,0,1,1,0,0,1,0,1])
            sage: ECA3 = cellular_automata.Elementary(110, [1,0,0,1,1,0,0,1,1,0])
            sage: ECA4 = cellular_automata.Elementary(110, [0,1,0,0,1,1,0,0,1,0,1])
            sage: ECA5 = cellular_automata.Elementary(110, [1,0,0,1,1,0,0,1,0,1])
            sage: ECA6 = cellular_automata.Elementary(110, [1,0,0,1,1,0,0,1,0,1], boundary=(1, 1))
            sage: ECA1 == ECA5
            True
            sage: ECA1 == ECA2
            False
            sage: ECA1 == ECA3
            False
            sage: ECA1 == ECA4
            False
            sage: ECA1 == ECA6
            False
        """
        return (isinstance(other, ElementaryCellularAutomata)
                and self._rule == other._rule
                and self._width == other._width
                and self._states[0] == other._states[0]
                and self._bdry == other._bdry)

    def __ne__(self, other):
        """
        Check non equality.

        TESTS::

            sage: ECA1 = cellular_automata.Elementary(110, [1,0,0,1,1,0,0,1,0,1])
            sage: ECA2 = cellular_automata.Elementary(101, [1,0,0,1,1,0,0,1,0,1])
            sage: ECA3 = cellular_automata.Elementary(110, [1,0,0,1,1,0,0,1,1,0])
            sage: ECA4 = cellular_automata.Elementary(110, [0,1,0,0,1,1,0,0,1,0,1])
            sage: ECA5 = cellular_automata.Elementary(110, [1,0,0,1,1,0,0,1,0,1])
            sage: ECA6 = cellular_automata.Elementary(110, [1,0,0,1,1,0,0,1,0,1], boundary=(1, 1))
            sage: ECA1 != ECA2
            True
            sage: ECA1 != ECA3
            True
            sage: ECA1 != ECA4
            True
            sage: ECA1 != ECA5
            False
            sage: ECA1 != ECA6
            True
        """
        return not (self == other)

    # Evolution functions
    # -------------------

    def evolve(self, number=None):
        r"""
        Evolve ``self``.

        INPUT:

        - ``number`` -- (optional) the number of times to perform
          the evolution

        EXAMPLES::

            sage: ECA = cellular_automata.Elementary(110, [1,0,0,1,1,0,0,1,0,1])
            sage: ascii_art(ECA)
            X  XX  X X
            sage: ECA.evolve()
            sage: ascii_art(ECA)
            X  XX  X X
            X XXX XXXX
            sage: ECA.evolve(10)
            sage: ascii_art(ECA)
            X  XX  X X
            X XXX XXXX
            XXX XXX  X
            X XXX X XX
            XXX XXXXXX
            X XXX    X
            XXX X   XX
            X XXX  XXX
            XXX X XX X
            X XXXXXXXX
            XXX      X
            X X     XX
        """
        if number is not None:
            for _ in range(number):
                self.evolve()
            return

        prev_state = self._states[-1]
        next_state = [None] * self._width

        def to_int(triple):
            return ZZ(list(reversed(triple)), base=2)
        if self._bdry is None:
            next_state[0] = self._rule[to_int([prev_state[-1]] + prev_state[:2])]
            next_state[-1] = self._rule[to_int(prev_state[-2:] + [prev_state[0]])]
        else:
            n = len(self._states)
            next_state[0] = self._rule[to_int([self._lbdry(n)] + prev_state[:2])]
            next_state[-1] = self._rule[to_int(prev_state[-2:] + [self._rbdry(n)])]

        for i in range(1, self._width-1):
            next_state[i] = self._rule[to_int(prev_state[i-1:i+2])]
        self._states.append(next_state)

    # Output functions
    # ----------------

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: ECA = cellular_automata.Elementary(123, width=10, initial_state=[1])
            sage: ECA
            Elementary cellular automata with rule 123 and initial state
             [0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
        """
        return "Elementary cellular automata with rule {} and initial state {}".format(
                        ZZ(self._rule, base=2), self._states[0])

    def print_state(self, number=None):
        r"""
        Print the state ``number``.

        INPUT:

        - ``number`` -- (default: the current state) the state to print

        EXAMPLES::

            sage: ECA = cellular_automata.Elementary(110, width=10,
            ....:                                    initial_state=[1,0,0,1,1,0,1])
            sage: ECA.print_state(15)
            X  X XXXXX
            sage: ECA.print_state(10)
            X    X  XX
            sage: ECA.print_state(20)
            X      XXX
            sage: for i in range(11):
            ....:     ECA.print_state(i)
               X  XX X
              XX XXXXX
             XXXXX   X
            XX   X  XX
            XX  XX XXX
            XX XXXXX X
            XXXX   XXX
            X  X  XX X
            X XX XXXXX
            XXXXXX   X
            X    X  XX
        """
        if number is None:
            number = len(self._states) - 1
        if number + 1 > len(self._states):
            for dummy in range(number + 1 - len(self._states)):
                self.evolve()

        state = self._states[number]
        print(''.join('X' if x else ' ' for x in state))

    def print_states(self, number=None):
        r"""
        Print the first ``num`` states of ``self``.

        .. NOTE::

            If the number of states computed for ``self`` is less than
            ``num``, then this evolves the system using the default
            time evolution.

        INPUT:

        - ``number`` -- the number of states to print

        EXAMPLES::

            sage: ECA = cellular_automata.Elementary(110, width=10,
            ....:                                    initial_state=[1,0,0,1,1,0,1])
            sage: ECA.print_states(10)
               X  XX X
              XX XXXXX
             XXXXX   X
            XX   X  XX
            XX  XX XXX
            XX XXXXX X
            XXXX   XXX
            X  X  XX X
            X XX XXXXX
            XXXXXX   X
        """
        for i in range(number):
            self.print_state(i)

    def _ascii_art_(self):
        r"""
        Return an ascii art representation of ``self``.

        EXAMPLES::

            sage: ECA = cellular_automata.Elementary(22, width=30, initial_state=[1])
            sage: ECA.evolve(30)
            sage: ascii_art(ECA)
                                         X
                                        XX
                                       X
                                      XXX
                                     X   X
                                    XXX XX
                                   X
                                  XXX
                                 X   X
                                XXX XXX
                               X       X
                              XXX     XXX
                             X   X   X   X
                            XXX XXX XXX XX
                           X
                          XXX
                         X   X
                        XXX XXX
                       X       X
                      XXX     XXX
                     X   X   X   X
                    XXX XXX XXX XXX
                   X               X
                  XXX             XXX
                 X   X           X   X
                XXX XXX         XXX XXX
               X       X       X       X
              XXX     XXX     XXX     XXX
             X   X   X   X   X   X   X   X
            XXX XXX XXX XXX XXX XXX XXX XX
        """
        return AsciiArt([''.join('X' if x else ' ' for x in state)
                         for state in self._states])

    def _unicode_art_(self):
        r"""
        Return a unicode art representation of ``self``.

        EXAMPLES::

            sage: ECA = cellular_automata.Elementary(22, width=30, initial_state=[1])
            sage: ECA.evolve(30)
            sage: unicode_art(ECA)
                                         █
                                        ██
                                       █
                                      ███
                                     █   █
                                    ███ ██
                                   █
                                  ███
                                 █   █
                                ███ ███
                               █       █
                              ███     ███
                             █   █   █   █
                            ███ ███ ███ ██
                           █
                          ███
                         █   █
                        ███ ███
                       █       █
                      ███     ███
                     █   █   █   █
                    ███ ███ ███ ███
                   █               █
                  ███             ███
                 █   █           █   █
                ███ ███         ███ ███
               █       █       █       █
              ███     ███     ███     ███
             █   █   █   █   █   █   █   █
            ███ ███ ███ ███ ███ ███ ███ ██
        """
        return UnicodeArt([''.join('█' if x else ' ' for x in state)
                           for state in self._states])

    def plot(self, number=None):
        r"""
        Return a plot of ``self``.

        INPUT:

        - ``number`` -- the number of states to plot

        EXAMPLES::

            sage: ECA = cellular_automata.Elementary(110, width=256)
            sage: ECA.evolve(256)
            sage: ECA.plot()                                                            # needs sage.plot
            Graphics object consisting of 1 graphics primitive
        """
        if number is None:
            number = len(self._states)
        if number > len(self._states):
            for dummy in range(number - len(self._states)):
                self.evolve()
        M = matrix(self._states[:number])
        return matrix_plot(M)
