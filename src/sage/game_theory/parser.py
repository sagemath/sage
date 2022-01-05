"""
Parser For gambit And lrs Nash Equilibria
"""

# ****************************************************************************
#       Copyright (C) 2014 James Campbell james.campbell@tanti.org.uk
#                     2015 Vincent Knight
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

class Parser():
    r"""
    A class for parsing the outputs of different algorithms called in other
    software packages.

    Two parsers are included, one for the ``'lrs'`` algorithm and another for
    the ``'LCP'`` algorithm.
    """

    def __init__(self, raw_string):
        """
        Initialise a Parser instance by storing a ``raw_string``
        (currently only used with H representation of a game).

        EXAMPLES::

            sage: from sage.cpython.string import bytes_to_str
            sage: from sage.game_theory.parser import Parser
            sage: from subprocess import Popen, PIPE
            sage: A = matrix([[1]])
            sage: B = matrix([[5]])
            sage: g = NormalFormGame([A,B])
            sage: game_str = g._lrs_nash_format(A, B)
            sage: game_name = tmp_filename()
            sage: with open(game_name, 'w') as game_file:
            ....:     _ = game_file.write(game_str)
            sage: process = Popen(['lrsnash', game_name], stdout=PIPE, stderr=PIPE)  # optional - lrslib
            sage: lrs_output = [bytes_to_str(row) for row in process.stdout]         # optional - lrslib
            sage: Parser(lrs_output).format_lrs()                                    # optional - lrslib
            [[(1,), (1,)]]

        This class is also used to parse the output of algorithms from
        the gambit python interface using :meth:`format_gambit()`.
        """
        self.raw_string = raw_string

    def format_lrs(self, legacy_format=False):
        r"""
        Parses the output of lrs so as to return vectors
        corresponding to equilibria.

        TESTS::

            sage: from sage.cpython.string import bytes_to_str
            sage: from sage.game_theory.parser import Parser
            sage: from subprocess import Popen, PIPE
            sage: A = matrix([[1, 2], [3, 2]])
            sage: g = NormalFormGame([A])
            sage: game_str = g._lrs_nash_format(A, -A)
            sage: game_name = tmp_filename()
            sage: with open(game_name, 'w') as game_file:
            ....:     _ = game_file.write(game_str)
            sage: process = Popen(['lrsnash', game_name], stdout=PIPE, stderr=PIPE)  # optional - lrslib
            sage: lrs_output = [bytes_to_str(row) for row in process.stdout]         # optional - lrslib

        The above creates a game, writes the H representations to
        temporary files, calls lrs and stores the output in ``lrs_output``
        (here slicing to get rid of some system parameters that get returned)::

            sage: lrs_output[:-2]                                                    # optional - lrslib
            [...,
             '2  0  1  2 \n',
             '1  1/2  1/2 -2 \n',
             '\n',
             '2  0  1  2 \n',
             '1  0  1 -2 \n',
             '\n',
             '*Number of equilibria found: 2\n',
             '*Player 1: vertices=3 bases=3 pivots=5\n',
             '*Player 2: vertices=2 bases=1 pivots=6\n',
             '\n']

        The above is pretty messy, here is the output when we put it through
        the parser::

            sage: nasheq = Parser(lrs_output).format_lrs()                           # optional - lrslib
            sage: nasheq                                                             # optional - lrslib
            [[(1/2, 1/2), (0, 1)], [(0, 1), (0, 1)]]

        Another game::

            sage: A = matrix([[-7, -5, 5],
            ....:             [5, 5, 3],
            ....:             [1, -6, 1]])
            sage: B = matrix([[-9, 7, 9],
            ....:             [6, -2, -3],
            ....:             [-4, 6, -10]])
            sage: g = NormalFormGame([A, B])
            sage: game_str = g._lrs_nash_format(A, B)
            sage: game_name = tmp_filename()
            sage: with open(game_name, 'w') as game_file:
            ....:     _ = game_file.write(game_str)
            sage: process = Popen(['lrsnash', game_name], stdout=PIPE, stderr=PIPE)  # optional - lrslib
            sage: lrs_output = [bytes_to_str(row) for row in process.stdout]         # optional - lrslib
            sage: print(lrs_output[:-2])                                             # optional - lrslib
            [...,
             '2  0  1/6  5/6  10/3 \n',
             '2  1/7  0  6/7  23/7 \n',
             '1  1/3  2/3  0  1 \n',
             '\n',
             '2  0  0  1  5 \n',
             '1  1  0  0  9 \n',
             '\n',
             '2  1  0  0  5 \n',
             '1  0  1  0  6 \n',
             '\n',
             '*Number of equilibria found: 4\n',
             '*Player 1: vertices=6 bases=7 pivots=10\n',
             '*Player 2: vertices=4 bases=2 pivots=14\n',
             '\n']

            sage: nasheq = Parser(lrs_output).format_lrs()                           # optional - lrslib
            sage: sorted(nasheq)                                                     # optional - lrslib
            [[(0, 1, 0), (1, 0, 0)],
             [(1/3, 2/3, 0), (0, 1/6, 5/6)],
             [(1/3, 2/3, 0), (1/7, 0, 6/7)],
             [(1, 0, 0), (0, 0, 1)]]

        TESTS:

        An example with the legacy format::

            sage: from sage.cpython.string import bytes_to_str
            sage: from sage.game_theory.parser import Parser
            sage: from subprocess import Popen, PIPE
            sage: A = matrix([[1, 2], [3, 2]])
            sage: g = NormalFormGame([A])
            sage: game1_str, game2_str = g._Hrepresentation(A, -A)
            doctest:warning...
            DeprecationWarning: NormalFormGame._Hrepresentation is deprecated as it
            creates the legacy input format. Use NormalFormGame._lrs_nash_format instead
            See https://trac.sagemath.org/27745 for details.
            sage: g1_name, g2_name = tmp_filename(), tmp_filename()
            sage: g1_file, g2_file = open(g1_name, 'w'), open(g2_name, 'w')
            sage: _ = g1_file.write(game1_str)
            sage: g1_file.close()
            sage: _ = g2_file.write(game2_str)
            sage: g2_file.close()
            sage: process = Popen(['lrsnash', g1_name, g2_name], stdout=PIPE, stderr=PIPE)  # not tested, optional - lrslib
            sage: lrs_output = [bytes_to_str(row) for row in process.stdout]                # not tested, optional - lrslib
            sage: nasheq = Parser(lrs_output).format_lrs(legacy_format=True)                # not tested, optional - lrslib
            sage: nasheq                                                                    # not tested, optional - lrslib
            [[(1/2, 1/2), (0, 1)], [(0, 1), (0, 1)]]
        """
        equilibria = []
        from sage.misc.sage_eval import sage_eval
        from itertools import groupby, dropwhile
        lines = iter(self.raw_string)
        if legacy_format:
            # Skip until the magic stars announce the beginning of the real output
            while not next(lines).startswith("*****"):
                pass
        else:
            # Skip comment lines starting with a single star
            lines = dropwhile(lambda line: line.startswith('*'), lines)
        for collection in [list(x[1]) for x in groupby(lines, lambda x: x == '\n')]:
            if collection[0].startswith('2'):
                s1 = tuple([sage_eval(k) for k in collection[-1].split()][1:-1])
                for s2 in collection[:-1]:
                    s2 = tuple([sage_eval(k) for k in s2.split()][1:-1])
                    equilibria.append([s1, s2])

        return equilibria

    def format_gambit(self, gambit_game):
        r"""
        Parses the output of gambit so as to return vectors
        corresponding to equilibria obtained using the LCP algorithm.

        TESTS:

        Here we construct a two by two game in gambit::

            sage: import gambit  # optional - gambit
            sage: from sage.game_theory.parser import Parser
            sage: g = gambit.Game.new_table([2,2])  # optional - gambit
            sage: g[int(0), int(0)][int(0)] = int(2)  # optional - gambit
            sage: g[int(0), int(0)][int(1)] = int(1)  # optional - gambit
            sage: g[int(0), int(1)][int(0)] = int(0)  # optional - gambit
            sage: g[int(0), int(1)][int(1)] = int(0)  # optional - gambit
            sage: g[int(1), int(0)][int(0)] = int(0)  # optional - gambit
            sage: g[int(1), int(0)][int(1)] = int(0)  # optional - gambit
            sage: g[int(1), int(1)][int(0)] = int(1)  # optional - gambit
            sage: g[int(1), int(1)][int(1)] = int(2)  # optional - gambit
            sage: solver = gambit.nash.ExternalLCPSolver()  # optional - gambit

        Here is the output of the LCP algorithm::

            sage: LCP_output = solver.solve(g)  # optional - gambit
            sage: LCP_output  # optional - gambit
            [<NashProfile for '': [[1.0, 0.0], [1.0, 0.0]]>,
             <NashProfile for '': [[0.6666666667, 0.3333333333], [0.3333333333, 0.6666666667]]>,
             <NashProfile for '': [[0.0, 1.0], [0.0, 1.0]]>]

        The Parser class outputs the equilibrium::

            sage: nasheq = Parser(LCP_output).format_gambit(g)  # optional - gambit
            sage: nasheq  # optional - gambit
            [[(1.0, 0.0), (1.0, 0.0)], [(0.6666666667, 0.3333333333), (0.3333333333, 0.6666666667)], [(0.0, 1.0), (0.0, 1.0)]]

        Here is another game::

            sage: g = gambit.Game.new_table([2,2])  # optional - gambit
            sage: g[int(0), int(0)][int(0)] = int(4)  # optional - gambit
            sage: g[int(0), int(0)][int(1)] = int(8)  # optional - gambit
            sage: g[int(0), int(1)][int(0)] = int(0)  # optional - gambit
            sage: g[int(0), int(1)][int(1)] = int(1)  # optional - gambit
            sage: g[int(1), int(0)][int(0)] = int(1)  # optional - gambit
            sage: g[int(1), int(0)][int(1)] = int(3)  # optional - gambit
            sage: g[int(1), int(1)][int(0)] = int(1)  # optional - gambit
            sage: g[int(1), int(1)][int(1)] = int(0)  # optional - gambit
            sage: solver = gambit.nash.ExternalLCPSolver()  # optional - gambit

        Here is the LCP output::

            sage: LCP_output = solver.solve(g)  # optional - gambit
            sage: LCP_output  # optional - gambit
            [<NashProfile for '': [[1.0, 0.0], [1.0, 0.0]]>]

        The corresponding parsed equilibrium::

            sage: nasheq = Parser(LCP_output).format_gambit(g)  # optional - gambit
            sage: nasheq  # optional - gambit
            [[(1.0, 0.0), (1.0, 0.0)]]

        Here is a larger degenerate game::

            sage: g = gambit.Game.new_table([3,3])  # optional - gambit
            sage: g[int(0), int(0)][int(0)] = int(-7)  # optional - gambit
            sage: g[int(0), int(0)][int(1)] = int(-9)  # optional - gambit
            sage: g[int(0), int(1)][int(0)] = int(-5)  # optional - gambit
            sage: g[int(0), int(1)][int(1)] = int(7)  # optional - gambit
            sage: g[int(0), int(2)][int(0)] = int(5)  # optional - gambit
            sage: g[int(0), int(2)][int(1)] = int(9)  # optional - gambit
            sage: g[int(1), int(0)][int(0)] = int(5)  # optional - gambit
            sage: g[int(1), int(0)][int(1)] = int(6)  # optional - gambit
            sage: g[int(1), int(1)][int(0)] = int(5)  # optional - gambit
            sage: g[int(1), int(1)][int(1)] = int(-2)  # optional - gambit
            sage: g[int(1), int(2)][int(0)] = int(3)  # optional - gambit
            sage: g[int(1), int(2)][int(1)] = int(-3)  # optional - gambit
            sage: g[int(2), int(0)][int(0)] = int(1)  # optional - gambit
            sage: g[int(2), int(0)][int(1)] = int(-4)  # optional - gambit
            sage: g[int(2), int(1)][int(0)] = int(-6)  # optional - gambit
            sage: g[int(2), int(1)][int(1)] = int(6)  # optional - gambit
            sage: g[int(2), int(2)][int(0)] = int(1)  # optional - gambit
            sage: g[int(2), int(2)][int(1)] = int(-10)  # optional - gambit
            sage: solver = gambit.nash.ExternalLCPSolver()  # optional - gambit

        Here is the LCP output::

            sage: LCP_output = solver.solve(g)  # optional - gambit
            sage: LCP_output  # optional - gambit
            [<NashProfile for '': [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]>,
             <NashProfile for '': [[0.3333333333, 0.6666666667, 0.0], [0.1428571429, 0.0, 0.8571428571]]>,
             <NashProfile for '': [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0]]>]

        The corresponding parsed equilibrium::

            sage: nasheq = Parser(LCP_output).format_gambit(g)  # optional - gambit
            sage: nasheq  # optional - gambit
            [[(1.0, 0.0, 0.0), (0.0, 0.0, 1.0)],
             [(0.3333333333, 0.6666666667, 0.0), (0.1428571429, 0.0, 0.8571428571)],
             [(0.0, 1.0, 0.0), (1.0, 0.0, 0.0)]]

        Note, that this differs from the same output of the lrs algorithm due
        the fact that the game is degenerate.
        """
        nice_stuff = []
        for gambitstrategy in self.raw_string:
            gambitstrategy = list(gambitstrategy)
            profile = [tuple(gambitstrategy[:len(gambit_game.players[int(0)].strategies)])]
            for player in list(gambit_game.players)[1:]:
                previousplayerstrategylength = len(profile[-1])
                profile.append(tuple(gambitstrategy[previousplayerstrategylength: previousplayerstrategylength + len(player.strategies)]))
            nice_stuff.append(profile)

        return nice_stuff

