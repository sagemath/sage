

class Parser():
    r"""
    A class for parsing the outputs of different algorithms called in other
    software packages.

    Two parsers are included, one for the ``'lrs'`` algorithm and another for
    the ``'LCP'`` algorithm.
    """

    def __init__(self, raw_string):
        """
        Initialise a Parser instance by storing a raw_string
        (currently only used with H representation of a game).

        TESTS:

        Simply checking that we have the correct string output
        for the H representation (which is the format required
        for the ``'lrs'`` algorithm)::

            sage: from sage.game_theory.parser import Parser
            sage: A = matrix([[1, 2], [3, 2]])
            sage: g = NormalFormGame([A])
            sage: raw_string = g._Hrepresentation(A, -A)
            sage: P = Parser(raw_string)
            sage: print(P.raw_string[0])
            H-representation
            linearity 1 5
            begin
            5 4 rational
            0 1 0 0
            0 0 1 0
            0 1 3  1
            0 2 2  1
            -1 1 1 0
            end
            <BLANKLINE>

            sage: print(P.raw_string[1])
            H-representation
            linearity 1 5
            begin
            5 4 rational
            0 -1 -2  1
            0 -3 -2  1
            0 1 0 0
            0 0 1 0
            -1 1 1 0
            end
            <BLANKLINE>

        The specific case of a game with 1 strategy for each player::

            sage: A = matrix([[1]])
            sage: B = matrix([[5]])
            sage: g = NormalFormGame([A,B])
            sage: raw_string = g._Hrepresentation(A, B)
            sage: P = Parser(raw_string)
            sage: print(P.raw_string[0])
            H-representation
            linearity 1 3
            begin
            3 3 rational
            0 1 0
            0 -5  1
            -1 1 0
            end
            <BLANKLINE>

            sage: print(P.raw_string[1])
            H-representation
            linearity 1 3
            begin
            3 3 rational
            0 -1  1
            0 1 0
            -1 1 0
            end
            <BLANKLINE>

        Another test::

            sage: from sage.game_theory.parser import Parser
            sage: A = matrix([[-7, -5, 5],
            ....:             [5, 5, 3],
            ....:             [1, -6, 1]])
            sage: B = matrix([[-9, 7, 9],
            ....:             [6, -2, -3],
            ....:             [-4, 6, -10]])
            sage: g = NormalFormGame([A, B])
            sage: raw_string = g._Hrepresentation(A, B)
            sage: P = Parser(raw_string)
            sage: print(P.raw_string[0])
            H-representation
            linearity 1 7
            begin
            7 5 rational
            0 1 0 0 0
            0 0 1 0 0
            0 0 0 1 0
            0  9 -6  4  1
            0 -7  2 -6  1
            0 -9  3  10  1
            -1 1 1 1 0
            end
            <BLANKLINE>

            sage: print(P.raw_string[1])
            H-representation
            linearity 1 7
            begin
            7 5 rational
            0 7 5 -5  1
            0 -5 -5 -3  1
            0 -1 6 -1  1
            0 1 0 0 0
            0 0 1 0 0
            0 0 0 1 0
            -1 1 1 1 0
            end
            <BLANKLINE>

        This class is also used to parse the output of algorithms from the gambit
        python interface using the `format_gambit` function.
        """
        self.raw_string = raw_string

    def format_lrs(self):
        """
        Parses the output of lrs so as to return vectors
        corresponding to equilibria.

        TESTS::

            sage: from sage.cpython.string import bytes_to_str
            sage: from sage.game_theory.parser import Parser
            sage: from subprocess import Popen, PIPE
            sage: A = matrix([[1, 2], [3, 2]])
            sage: g = NormalFormGame([A])
            sage: game1_str, game2_str = g._Hrepresentation(A, -A)
            sage: g1_name = tmp_filename()
            sage: g2_name = tmp_filename()
            sage: g1_file = open(g1_name, 'w')
            sage: g2_file = open(g2_name, 'w')
            sage: _ = g1_file.write(game1_str)
            sage: g1_file.close()
            sage: _ = g2_file.write(game2_str)
            sage: g2_file.close()
            sage: process = Popen(['lrsnash', g1_name, g2_name], stdout=PIPE, stderr=PIPE)  # optional - lrslib
            sage: lrs_output = [bytes_to_str(row) for row in process.stdout]  # optional - lrslib

        The above creates a game, writes the H representation to
        temporary files, calls lrs and stores the output in `lrs_output`
        (here slicing to get rid of some system parameters that get returned)::

            sage: lrs_output[5:16]  # optional - lrslib
            ['\n',
             '***** 4 4 rational\n',
             '2  0  1  2 \n',
             '1  1/2  1/2 -2 \n',
             '\n',
             '2  0  1  2 \n',
             '1  0  1 -2 \n',
             '\n',
             '\n',
             '*Number of equilibria found: 2\n',
             '*Player 1: vertices=3 bases=3 pivots=5\n']

        The above is pretty messy, here is the output when we put it through
        the parser::

            sage: nasheq = Parser(lrs_output).format_lrs()  # optional - lrslib
            sage: nasheq  # optional - lrslib
            [[(1/2, 1/2), (0, 1)], [(0, 1), (0, 1)]]

        Another game::

            sage: A = matrix([[-7, -5, 5],
            ....:             [5, 5, 3],
            ....:             [1, -6, 1]])
            sage: B = matrix([[-9, 7, 9],
            ....:             [6, -2, -3],
            ....:             [-4, 6, -10]])
            sage: g = NormalFormGame([A, B])
            sage: game1_str, game2_str = g._Hrepresentation(A, B)
            sage: g1_name = tmp_filename()
            sage: g2_name = tmp_filename()
            sage: g1_file = open(g1_name, 'w')
            sage: g2_file = open(g2_name, 'w')
            sage: _ = g1_file.write(game1_str)
            sage: g1_file.close()
            sage: _ = g2_file.write(game2_str)
            sage: g2_file.close()
            sage: process = Popen(['lrsnash', g1_name, g2_name], stdout=PIPE, stderr=PIPE)  # optional - lrslib
            sage: lrs_output = [bytes_to_str(row) for row in process.stdout]  # optional - lrslib
            sage: print(lrs_output[5:20])  # optional - lrslib
            ['\n',
             '***** 5 5 rational\n',
             '2  1/7  0  6/7  23/7 \n',
             '2  0  1/6  5/6  10/3 \n',
             '1  1/3  2/3  0  1 \n',
             '\n',
             '2  0  0  1  5 \n',
             '1  1  0  0  9 \n',
             '\n',
             '2  1  0  0  5 \n',
             '1  0  1  0  6 \n',
             '\n',
             '\n',
             '*Number of equilibria found: 4\n',
             '*Player 1: vertices=6 bases=7 pivots=10\n']

            sage: nasheq = Parser(lrs_output).format_lrs()  # optional - lrslib
            sage: sorted(nasheq)  # optional - lrslib
            [[(0, 1, 0), (1, 0, 0)],
             [(1/3, 2/3, 0), (0, 1/6, 5/6)],
             [(1/3, 2/3, 0), (1/7, 0, 6/7)],
             [(1, 0, 0), (0, 0, 1)]]
        """
        equilibria = []
        from sage.misc.sage_eval import sage_eval
        from itertools import groupby
        for collection in [list(x[1]) for x in groupby(self.raw_string[7:], lambda x: x == '\n')]:
            if collection[0].startswith('2'):
                s1 = tuple([sage_eval(k) for k in collection[-1].split()][1:-1])
                for s2 in collection[:-1]:
                    s2 = tuple([sage_eval(k) for k in s2.split()][1:-1])
                    equilibria.append([s1, s2])

        return equilibria

    def format_gambit(self, gambit_game):
        """
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
