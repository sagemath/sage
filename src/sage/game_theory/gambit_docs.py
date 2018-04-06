"""
Using Gambit as a standalone package

This file contains some information and tests for the use of
`Gambit <http://www.gambit-project.org/>`_ as a stand alone package.

To install gambit as an optional package run (from root of Sage)::

    $ sage -i gambit

The `python API documentation for gambit
<http://www.gambit-project.org/gambit14/pyapi.html>`_ shows various examples
that can be run easily in IPython. To run the IPython packaged with Sage run
(from root of Sage)::

    $ ./sage --ipython

Here is an example that constructs the Prisoner's Dilemma::

    In [1]: import gambit
    In [2]: g = gambit.Game.new_table([2,2])
    In [3]: g.title = "A prisoner's dilemma game"
    In [4]: g.players[0].label = "Alphonse"
    In [5]: g.players[1].label = "Gaston"
    In [6]: g
    Out[6]:
    NFG 1 R "A prisoner's dilemma game" { "Alphonse" "Gaston" }

    { { "1" "2" }
    { "1" "2" }
    }
    ""

    {
    { "" 0, 0 }
    { "" 0, 0 }
    { "" 0, 0 }
    { "" 0, 0 }
    }
    1 2 3 4

    In [7]: g.players[0].strategies
    Out[7]: [<Strategy [0] '1' for player 'Alphonse' in game 'A
    prisoner's dilemma game'>,
             <Strategy [1] '2' for player 'Alphonse' in game 'A prisoner's dilemma game'>]
    In [8]: len(g.players[0].strategies)
    Out[8]: 2

    In [9]: g.players[0].strategies[0].label = "Cooperate"
    In [10]: g.players[0].strategies[1].label = "Defect"
    In [11]: g.players[0].strategies
    Out[11]: [<Strategy [0] 'Cooperate' for player 'Alphonse' in game 'A
    prisoner's dilemma game'>,
       <Strategy [1] 'Defect' for player 'Alphonse' in game 'A prisoner's dilemma game'>]

    In [12]: g[0,0][0] = 8
    In [13]: g[0,0][1] = 8
    In [14]: g[0,1][0] = 2
    In [15]: g[0,1][1] = 10
    In [16]: g[1,0][0] = 10
    In [17]: g[1,1][1] = 2
    In [18]: g[1,0][1] = 2
    In [19]: g[1,1][0] = 5
    In [20]: g[1,1][1] = 5

Here is a list of the various solvers available in gambit:

- ExternalEnumPureSolver
- ExternalEnumMixedSolver
- ExternalLPSolver
- ExternalLCPSolver
- ExternalSimpdivSolver
- ExternalGlobalNewtonSolver
- ExternalEnumPolySolver
- ExternalLyapunovSolver
- ExternalIteratedPolymatrixSolver
- ExternalLogitSolver

Here is how to use the ``ExternalEnumPureSolver``::

    In [21]: solver = gambit.nash.ExternalEnumPureSolver()
    In [22]: solver.solve(g)
    Out[22]: [<NashProfile for 'A prisoner's dilemma game': [Fraction(0, 1), Fraction(1, 1), Fraction(0, 1), Fraction(1, 1)]>]

Note that the above finds the equilibria by investigating all potential pure
pure strategy pairs. This will fail to find all Nash equilibria in certain
games.  For example here is an implementation of Matching Pennies::

    In [1]: import gambit
    In [2]: g = gambit.Game.new_table([2,2])
    In [3]: g[0, 0][0] = 1
    In [4]: g[0, 0][1] = -1
    In [5]: g[0, 1][0] = -1
    In [6]: g[0, 1][1] = 1
    In [7]: g[1, 0][0] = -1
    In [8]: g[1, 0][1] = 1
    In [9]: g[1, 1][0] = 1
    In [10]: g[1, 1][1] = -1
    In [11]: solver = gambit.nash.ExternalEnumPureSolver()
    In [12]: solver.solve(g)
    Out[12]: []

If we solve this with the ``LCP`` solver we get the expected Nash equilibrium::

    In [13]: solver = gambit.nash.ExternalLCPSolver()
    In [14]: solver.solve(g)
    Out[14]: [<NashProfile for '': [0.5, 0.5, 0.5, 0.5]>]

Note that the above examples only show how to build and find equilibria for
two player strategic form games. Gambit supports multiple player games as well
as extensive form games: for more details see http://www.gambit-project.org/.

If one really wants to use gambit directly in Sage (without using the
``NormalFormGame`` class as a wrapper) then integers must first be
converted to Python integers (due to the preparser). Here is an example
showing the Battle of the Sexes::

    sage: import gambit  # optional - gambit
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
    sage: solver.solve(g)  # optional - gambit
    [<NashProfile for '': [[1.0, 0.0], [1.0, 0.0]]>,
     <NashProfile for '': [[0.6666666667, 0.3333333333], [0.3333333333, 0.6666666667]]>,
     <NashProfile for '': [[0.0, 1.0], [0.0, 1.0]]>]

AUTHOR:

- Vince Knight (11-2014): Original version

"""
