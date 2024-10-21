# sage_setup: distribution = sagemath-polyhedra
r"""
Normal form games with N players.

This module implements a class for normal form games (strategic form games)
[NN2007]_. At present the following algorithms are implemented to
compute equilibria of these games:

 * ``'enumeration'`` - An implementation of the support enumeration
   algorithm built in Sage.

 * ``'LCP'`` - An interface with the 'gambit' solver's implementation
   of the Lemke-Howson algorithm.

 * ``'lp'`` - A built-in Sage implementation (with a gambit alternative)
   of a zero-sum game solver using linear programming. See
   :class:`MixedIntegerLinearProgram` for more on MILP solvers in Sage.

 * ``'lrs'`` - A solver interfacing with the 'lrslib' library.

The architecture for the class is based on the gambit architecture to
ensure an easy transition between gambit and Sage.  At present the
algorithms for the computation of equilibria only solve 2 player games.

A very simple and well known example of normal form game is referred
to as the 'Battle of the Sexes' in which two players Amy and Bob
are modeled.  Amy prefers to play video games and Bob prefers to
watch a movie.  They both however want to spend their evening together.
This can be modeled using the following two matrices:

.. MATH::

    A = \begin{pmatrix}
        3&1\\
        0&2\\
        \end{pmatrix}


    B = \begin{pmatrix}
        2&1\\
        0&3\\
        \end{pmatrix}

Matrix `A` represents the utilities of Amy and matrix `B` represents the
utility of Bob. The choices of Amy correspond to the rows of the matrices:

* The first row corresponds to video games.

* The second row corresponds to movies.

Similarly Bob's choices are represented by the columns:

* The first column corresponds to video games.

* The second column corresponds to movies.

Thus, if both Amy and Bob choose to play video games: Amy receives a
utility of 3 and Bob a utility of 2. If Amy is indeed going to stick
with video games Bob has no incentive to deviate (and vice versa).

This situation repeats itself if both Amy and Bob choose to watch a movie:
neither has an incentive to deviate.

This loosely described situation is referred to as a Nash Equilibrium.
We can use Sage to find them, and more importantly, see if there is any
other situation where Amy and Bob have no reason to change their choice
of action:

Here is how we create the game in Sage::

    sage: A = matrix([[3, 1], [0, 2]])
    sage: B = matrix([[2, 1], [0, 3]])
    sage: battle_of_the_sexes = NormalFormGame([A, B])
    sage: battle_of_the_sexes
    Normal Form Game with the following utilities: {(0, 0): [3, 2],
     (0, 1): [1, 1], (1, 0): [0, 0], (1, 1): [2, 3]}

To obtain the Nash equilibria we run the ``obtain_nash()`` method. In the
first few examples, we will use the 'support enumeration' algorithm.
A discussion about the different algorithms will be given later::

    sage: battle_of_the_sexes.obtain_nash(algorithm='enumeration')
    [[(0, 1), (0, 1)], [(3/4, 1/4), (1/4, 3/4)], [(1, 0), (1, 0)]]

If we look a bit closer at our output we see that a list of three
pairs of tuples have been returned. Each of these correspond to a
Nash Equilibrium, represented as a probability distribution over the
available strategies:

* `[(1, 0), (1, 0)]` corresponds to the first player only
  playing their first strategy and the second player also only playing
  their first strategy. In other words Amy and Bob both play video games.

* `[(0, 1), (0, 1)]` corresponds to the first player only
  playing their second strategy and the second player also only playing
  their second strategy. In other words Amy and Bob both watch movies.

* `[(3/4, 1/4), (1/4, 3/4)]` corresponds to players `mixing` their
  strategies. Amy plays video games 75% of the time and Bob watches
  movies 75% of the time. At this equilibrium point Amy and Bob will
  only ever do the same activity `3/8` of the time.

We can use Sage to compute the expected utility for any mixed strategy
pair `(\sigma_1, \sigma_2)`. The payoff to player 1 is given by the
vector/matrix multiplication:

.. MATH::

    \sigma_1 A \sigma_2

The payoff to player 2 is given by:

.. MATH::

    \sigma_1 B \sigma_2

To compute this in Sage we have::

    sage: for ne in battle_of_the_sexes.obtain_nash(algorithm='enumeration'):
    ....:     print("Utility for {}: ".format(ne))
    ....:     print("{} {}".format(vector(ne[0]) * A * vector(ne[1]), vector(ne[0]) * B * vector(ne[1])))
    Utility for [(0, 1), (0, 1)]:
    2 3
    Utility for [(3/4, 1/4), (1/4, 3/4)]:
    3/2 3/2
    Utility for [(1, 0), (1, 0)]:
    3 2

Allowing players to play mixed strategies ensures that there will always
be a Nash Equilibrium for a normal form game. This result is called Nash's
Theorem ([Nas1950]_).

Let us consider the game called 'matching pennies' where two players each
present a coin with either HEADS or TAILS showing. If the coins show the
same side then player 1 wins, otherwise player 2 wins:


.. MATH::

    A = \begin{pmatrix}
        1&-1\\
        -1&1\\
        \end{pmatrix}


    B = \begin{pmatrix}
        -1&1\\
        1&-1\\
        \end{pmatrix}

It should be relatively straightforward to observe, that there is no
situation, where both players always do the same thing, and have no
incentive to deviate.

We can plot the utility of player 1 when player 2 is playing a mixed
strategy `\sigma_2 = (y, 1-y)` (so that the utility to player 1 for
playing strategy number `i` is given by the matrix/vector multiplication
`(Ay)_i`, ie element in position `i` of the matrix/vector multiplication
`Ay`) ::

    sage: y = var('y')                                                                  # needs sage.symbolic
    sage: A = matrix([[1, -1], [-1, 1]])
    sage: p = plot((A * vector([y, 1 - y]))[0], y, 0, 1, color='blue',                  # needs sage.symbolic
    ....:          legend_label='$u_1(r_1, (y, 1-y))$', axes_labels=['$y$', ''])
    sage: p += plot((A * vector([y, 1 - y]))[1], y, 0, 1, color='red',                  # needs sage.symbolic
    ....:           legend_label='$u_1(r_2, (y, 1-y))$'); p
    Graphics object consisting of 2 graphics primitives

We see that the only point at which player 1 is indifferent amongst
the available strategies is when `y = 1/2`.

If we compute the Nash equilibria we see that this corresponds to a point
at which both players are indifferent::

    sage: A = matrix([[1, -1], [-1, 1]])
    sage: B = matrix([[-1, 1], [1, -1]])
    sage: matching_pennies = NormalFormGame([A, B])
    sage: matching_pennies.obtain_nash(algorithm='enumeration')
    [[(1/2, 1/2), (1/2, 1/2)]]

The utilities to both players at this Nash equilibrium
is easily computed::

    sage: [vector([1/2, 1/2]) * M * vector([1/2, 1/2])
    ....:  for M in matching_pennies.payoff_matrices()]
    [0, 0]

Note that the above uses the ``payoff_matrices`` method
which returns the payoff matrices for a 2 player game::

    sage: matching_pennies.payoff_matrices()
    (
    [ 1 -1]  [-1  1]
    [-1  1], [ 1 -1]
    )

One can also input a single matrix and then a zero sum game is constructed.
Here is an instance of `Rock-Paper-Scissors-Lizard-Spock
<http://www.samkass.com/theories/RPSSL.html>`_::

    sage: A = matrix([[0, -1, 1, 1, -1],
    ....:             [1, 0, -1, -1, 1],
    ....:             [-1, 1, 0, 1 , -1],
    ....:             [-1, 1, -1, 0, 1],
    ....:             [1, -1, 1, -1, 0]])
    sage: g = NormalFormGame([A])
    sage: g.obtain_nash(algorithm='enumeration')
    [[(1/5, 1/5, 1/5, 1/5, 1/5), (1/5, 1/5, 1/5, 1/5, 1/5)]]

We can also study games where players aim to minimize their utility.
Here is the Prisoner's Dilemma (where players are aiming to reduce
time spent in prison)::

    sage: A = matrix([[2, 5], [0, 4]])
    sage: B = matrix([[2, 0], [5, 4]])
    sage: prisoners_dilemma = NormalFormGame([A, B])
    sage: prisoners_dilemma.obtain_nash(algorithm='enumeration', maximization=False)
    [[(0, 1), (0, 1)]]

When obtaining Nash equilibrium the following algorithms are
currently available:

* ``'lp'``: A solver for constant sum 2 player games using linear
  programming. This constructs a
  :mod:`MixedIntegerLinearProgram <sage.numerical.MILP>` using the
  solver which was passed in with ``solver`` to solve the linear
  programming representation of the game. See
  :class:`MixedIntegerLinearProgram` for more on MILP solvers in Sage.

* ``'lrs'``: Reverse search vertex enumeration for 2 player games. This
  algorithm uses the optional 'lrslib' package. To install it, type
  ``sage -i lrslib`` in the shell. For more information, see [Av2000]_.

* ``'LCP'``: Linear complementarity program algorithm for 2 player games.
  This algorithm uses the open source game theory package:
  `Gambit <http://gambit.sourceforge.net/>`_ [Gambit]_. At present this is
  the only gambit algorithm available in sage but further development will
  hope to implement more algorithms
  (in particular for games with more than 2 players).

* ``'enumeration'``: Support enumeration for 2 player games. This
  algorithm is hard coded in Sage and checks through all potential
  supports of a strategy. Supports of a given size with a conditionally
  dominated strategy are ignored. Note: this is not the preferred
  algorithm. The algorithm implemented is a combination of a basic
  algorithm described in [NN2007]_ and a pruning component described
  in [SLB2008]_.

Below we show how the these algorithms are called::

    sage: matching_pennies.obtain_nash(algorithm='lrs')  # optional - lrslib
    [[(1/2, 1/2), (1/2, 1/2)]]
    sage: matching_pennies.obtain_nash(algorithm='LCP')  # optional - gambit
    [[(0.5, 0.5), (0.5, 0.5)]]
    sage: matching_pennies.obtain_nash(algorithm='lp', solver='PPL')
    [[(1/2, 1/2), (1/2, 1/2)]]
    sage: matching_pennies.obtain_nash(algorithm='lp', solver='gambit') # optional - gambit
    [[(0.5, 0.5), (0.5, 0.5)]]
    sage: matching_pennies.obtain_nash(algorithm='enumeration')
    [[(1/2, 1/2), (1/2, 1/2)]]

Note that if no algorithm argument is passed then the default will be
selected according to the following order (if the corresponding package is
installed):

1. ``'lp'`` (if the game is constant-sum; uses the solver chosen by Sage)
2. ``'lrs'`` (requires 'lrslib')
3. ``'enumeration'``

Here is a game being constructed using gambit syntax (note that a
``NormalFormGame`` object acts like a dictionary with pure strategy tuples as
keys and payoffs as their values)::

    sage: f = NormalFormGame()
    sage: f.add_player(2)  # Adding first player with 2 strategies
    sage: f.add_player(2)  # Adding second player with 2 strategies
    sage: f[0,0][0] = 1
    sage: f[0,0][1] = 3
    sage: f[0,1][0] = 2
    sage: f[0,1][1] = 3
    sage: f[1,0][0] = 3
    sage: f[1,0][1] = 1
    sage: f[1,1][0] = 4
    sage: f[1,1][1] = 4
    sage: f
    Normal Form Game with the following utilities: {(0, 0): [1, 3],
     (0, 1): [2, 3], (1, 0): [3, 1], (1, 1): [4, 4]}

Once this game is constructed we can view the payoff matrices and solve the
game::

    sage: f.payoff_matrices()
    (
    [1 2]  [3 3]
    [3 4], [1 4]
    )
    sage: f.obtain_nash(algorithm='enumeration')
    [[(0, 1), (0, 1)]]

We can add an extra strategy to the first player::

    sage: f.add_strategy(0)
    sage: f
    Normal Form Game with the following utilities: {(0, 0): [1, 3],
     (0, 1): [2, 3],
     (1, 0): [3, 1],
     (1, 1): [4, 4],
     (2, 0): [False, False],
     (2, 1): [False, False]}

If we do this and try and obtain the Nash equilibrium or view the payoff
matrices(without specifying the utilities), an error is returned::

    sage: f.obtain_nash()
    Traceback (most recent call last):
    ...
    ValueError: utilities have not been populated
    sage: f.payoff_matrices()
    Traceback (most recent call last):
    ...
    ValueError: utilities have not been populated

Here we populate the missing utilities::

    sage: f[2, 1] = [5, 3]
    sage: f[2, 0] = [2, 1]
    sage: f.payoff_matrices()
    (
    [1 2]  [3 3]
    [3 4]  [1 4]
    [2 5], [1 3]
    )
    sage: f.obtain_nash()
    [[(0, 0, 1), (0, 1)]]

We can use the same syntax as above to create games with
more than 2 players::

    sage: threegame = NormalFormGame()
    sage: threegame.add_player(2)  # Adding first player with 2 strategies
    sage: threegame.add_player(2)  # Adding second player with 2 strategies
    sage: threegame.add_player(2)  # Adding third player with 2 strategies
    sage: threegame[0, 0, 0][0] = 3
    sage: threegame[0, 0, 0][1] = 1
    sage: threegame[0, 0, 0][2] = 4
    sage: threegame[0, 0, 1][0] = 1
    sage: threegame[0, 0, 1][1] = 5
    sage: threegame[0, 0, 1][2] = 9
    sage: threegame[0, 1, 0][0] = 2
    sage: threegame[0, 1, 0][1] = 6
    sage: threegame[0, 1, 0][2] = 5
    sage: threegame[0, 1, 1][0] = 3
    sage: threegame[0, 1, 1][1] = 5
    sage: threegame[0, 1, 1][2] = 8
    sage: threegame[1, 0, 0][0] = 9
    sage: threegame[1, 0, 0][1] = 7
    sage: threegame[1, 0, 0][2] = 9
    sage: threegame[1, 0, 1][0] = 3
    sage: threegame[1, 0, 1][1] = 2
    sage: threegame[1, 0, 1][2] = 3
    sage: threegame[1, 1, 0][0] = 8
    sage: threegame[1, 1, 0][1] = 4
    sage: threegame[1, 1, 0][2] = 6
    sage: threegame[1, 1, 1][0] = 2
    sage: threegame[1, 1, 1][1] = 6
    sage: threegame[1, 1, 1][2] = 4
    sage: threegame
    Normal Form Game with the following utilities: {(0, 0, 0): [3, 1, 4],
     (0, 0, 1): [1, 5, 9],
     (0, 1, 0): [2, 6, 5],
     (0, 1, 1): [3, 5, 8],
     (1, 0, 0): [9, 7, 9],
     (1, 0, 1): [3, 2, 3],
     (1, 1, 0): [8, 4, 6],
     (1, 1, 1): [2, 6, 4]}

The above requires a lot of input that could be simplified if there is
another data structure with our utilities and/or a structure to the
utilities.  The following example creates a game with a relatively strange
utility function::

    sage: def utility(strategy_triplet, player):
    ....:     return sum(strategy_triplet) * player
    sage: threegame = NormalFormGame()
    sage: threegame.add_player(2)  # Adding first player with 2 strategies
    sage: threegame.add_player(2)  # Adding second player with 2 strategies
    sage: threegame.add_player(2)  # Adding third player with 2 strategies
    sage: for i, j, k in [(i, j, k) for i in [0,1] for j in [0,1] for k in [0,1]]:
    ....:     for p in range(3):
    ....:          threegame[i, j, k][p] = utility([i, j, k], p)
    sage: threegame
    Normal Form Game with the following utilities: {(0, 0, 0): [0, 0, 0],
     (0, 0, 1): [0, 1, 2],
     (0, 1, 0): [0, 1, 2],
     (0, 1, 1): [0, 2, 4],
     (1, 0, 0): [0, 1, 2],
     (1, 0, 1): [0, 2, 4],
     (1, 1, 0): [0, 2, 4],
     (1, 1, 1): [0, 3, 6]}

At present no algorithm has been implemented in Sage for games with
more than 2 players::

    sage: threegame.obtain_nash()
    Traceback (most recent call last):
    ...
    NotImplementedError: Nash equilibrium for games with more than 2 players
     have not been implemented yet. Please see the gambit website
     (http://gambit.sourceforge.net/) that has a variety of available algorithms

There are however a variety of such algorithms available in gambit,
further compatibility between Sage and gambit is actively being developed:
https://github.com/tturocy/gambit/tree/sage_integration.

It can be shown that linear scaling of the payoff matrices conserves the
equilibrium values::

    sage: A = matrix([[2, 1], [1, 2.5]])
    sage: B = matrix([[-1, 3], [2, 1]])
    sage: g = NormalFormGame([A, B])
    sage: g.obtain_nash(algorithm='enumeration')
    [[(1/5, 4/5), (3/5, 2/5)]]
    sage: g.obtain_nash(algorithm='lrs') # optional - lrslib
    [[(1/5, 4/5), (3/5, 2/5)]]
    sage: A = 2 * A
    sage: g = NormalFormGame([A, B])
    sage: g.obtain_nash(algorithm='LCP')  # optional - gambit
    [[(0.2, 0.8), (0.6, 0.4)]]

It is also possible to generate a Normal form game from a gambit Game::

    sage: # optional - gambit
    sage: from gambit import Game
    sage: gambitgame= Game.new_table([2, 2])
    sage: gambitgame[int(0), int(0)][int(0)] = int(8)
    sage: gambitgame[int(0), int(0)][int(1)] = int(8)
    sage: gambitgame[int(0), int(1)][int(0)] = int(2)
    sage: gambitgame[int(0), int(1)][int(1)] = int(10)
    sage: gambitgame[int(1), int(0)][int(0)] = int(10)
    sage: gambitgame[int(1), int(0)][int(1)] = int(2)
    sage: gambitgame[int(1), int(1)][int(0)] = int(5)
    sage: gambitgame[int(1), int(1)][int(1)] = int(5)
    sage: g = NormalFormGame(gambitgame); g
    Normal Form Game with the following utilities: {(0, 0): [8.0, 8.0],
     (0, 1): [2.0, 10.0],
     (1, 0): [10.0, 2.0],
     (1, 1): [5.0, 5.0]}

For more information on using Gambit in Sage see: :mod:`Using Gambit in
Sage<sage.game_theory.gambit_docs>`. This includes how to access Gambit
directly using the version of iPython shipped with Sage and an explanation
as to why the ``int`` calls are needed to handle the Sage preparser.

Here is a slightly longer game that would take too long to solve with
``'enumeration'``. Consider the following:

An airline loses two suitcases belonging to two different travelers. Both
suitcases happen to be identical and contain identical antiques. An
airline manager tasked to settle the claims of both travelers explains
that the airline is liable for a maximum of 10 per suitcase, and in order
to determine an honest appraised value of the antiques the manager
separates both travelers so they can't confer, and asks them to write down
the amount of their value at no less than 2 and no larger than 10. He
also tells them that if both write down the same number, he will treat
that number as the true dollar value of both suitcases and reimburse both
travelers that amount.

However, if one writes down a smaller number than the other, this smaller
number will be taken as the true dollar value, and both travelers will
receive that amount along with a bonus/malus: 2 extra will be paid to the
traveler who wrote down the lower value and a 2 deduction will be taken
from the person who wrote down the higher amount. The challenge is: what
strategy should both travelers follow to decide the value they should
write down?

In the following we create the game (with a max value of 10) and solve it::

    sage: K = 10  # Modifying this value lets us play with games of any size
    sage: A = matrix([[min(i,j) + 2 * sign(j-i)  for j in range(K, 1, -1)]
    ....:             for i in range(K, 1, -1)])
    sage: B = matrix([[min(i,j) + 2 * sign(i-j)  for j in range(K, 1, -1)]
    ....:             for i in range(K, 1, -1)])
    sage: g = NormalFormGame([A, B])
    sage: g.obtain_nash(algorithm='lrs')  # optional - lrslib
    [[(0, 0, 0, 0, 0, 0, 0, 0, 1), (0, 0, 0, 0, 0, 0, 0, 0, 1)]]
    sage: g.obtain_nash(algorithm='LCP')  # optional - gambit
    [[(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0),
      (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0)]]

The output is a pair of vectors (as before) showing the Nash equilibrium.
In particular it here shows that out of the 10 possible strategies both
players should choose the last. Recall that the above considers a reduced
version of the game where individuals can claim integer values from 10
to 2.  The equilibrium strategy is thus for both players to state that
the value of their suitcase is 2.

Several standard Normal Form Games have also been implemented.
For more information on how to access these, see:
:mod:`Game Theory Catalog<sage.game_theory.catalog>`.
Included is information on the situation each Game models.
For example::

    sage: g = game_theory.normal_form_games.PrisonersDilemma()
    sage: g
    Prisoners dilemma - Normal Form Game with the following utilities: ...
    sage: d = {(0, 1): [-5, 0], (1, 0): [0, -5],
    ....:      (0, 0): [-2, -2], (1, 1): [-4, -4]}
    sage: g == d
    True
    sage: g.obtain_nash()
    [[(0, 1), (0, 1)]]

We can easily obtain the best response for a player to a given strategy.  In
this example we obtain the best responses for Player 1, when Player 2 uses two
different strategies::

    sage: A = matrix([[3, 0], [0, 3], [1.5, 1.5]])
    sage: B = matrix([[4, 3], [2, 6], [3, 1]])
    sage: g = NormalFormGame([A, B])
    sage: g.best_responses((1/2, 1/2), player=0)
    [0, 1, 2]
    sage: g.best_responses((3/4, 1/4), player=0)
    [0]

Here we do the same for player 2::

    sage: g.best_responses((4/5, 1/5, 0), player=1)
    [0, 1]

We see that for the game `Rock-Paper-Scissors-Lizard-Spock
<http://www.samkass.com/theories/RPSSL.html>`_ any pure strategy has two best
responses::

    sage: g = game_theory.normal_form_games.RPSLS()
    sage: A, B = g.payoff_matrices()
    sage: A, B
    (
    [ 0 -1  1  1 -1]  [ 0  1 -1 -1  1]
    [ 1  0 -1 -1  1]  [-1  0  1  1 -1]
    [-1  1  0  1 -1]  [ 1 -1  0 -1  1]
    [-1  1 -1  0  1]  [ 1 -1  1  0 -1]
    [ 1 -1  1 -1  0], [-1  1 -1  1  0]
    )
    sage: g.best_responses((1, 0, 0, 0, 0), player=0)
    [1, 4]
    sage: g.best_responses((0, 1, 0, 0, 0), player=0)
    [2, 3]
    sage: g.best_responses((0, 0, 1, 0, 0), player=0)
    [0, 4]
    sage: g.best_responses((0, 0, 0, 1, 0), player=0)
    [0, 2]
    sage: g.best_responses((0, 0, 0, 0, 1), player=0)
    [1, 3]
    sage: g.best_responses((1, 0, 0, 0, 0), player=1)
    [1, 4]
    sage: g.best_responses((0, 1, 0, 0, 0), player=1)
    [2, 3]
    sage: g.best_responses((0, 0, 1, 0, 0), player=1)
    [0, 4]
    sage: g.best_responses((0, 0, 0, 1, 0), player=1)
    [0, 2]
    sage: g.best_responses((0, 0, 0, 0, 1), player=1)
    [1, 3]

Note that degenerate games can cause problems for most algorithms.
The following example in fact has an infinite quantity of equilibria which
is evidenced by the various algorithms returning different solutions::

    sage: A = matrix([[3,3],[2,5],[0,6]])
    sage: B = matrix([[3,3],[2,6],[3,1]])
    sage: degenerate_game = NormalFormGame([A,B])
    sage: degenerate_game.obtain_nash(algorithm='lrs')  # random, optional - lrslib
    [[(0, 1/3, 2/3), (1/3, 2/3)], [(1, 0, 0), (1/2, 3)], [(1, 0, 0), (1, 3)]]
    sage: degenerate_game.obtain_nash(algorithm='LCP')  # optional - gambit
    [[(0.0, 0.3333333333, 0.6666666667), (0.3333333333, 0.6666666667)],
     [(1.0, -0.0, 0.0), (0.6666666667, 0.3333333333)],
     [(1.0, 0.0, 0.0), (1.0, 0.0)]]
    sage: degenerate_game.obtain_nash(algorithm='enumeration')
    [[(0, 1/3, 2/3), (1/3, 2/3)], [(1, 0, 0), (1, 0)]]

We can check the cause of this by using ``is_degenerate()``::

    sage: degenerate_game.is_degenerate()
    True

Note the 'negative' `-0.0` output by gambit. This is due to the numerical
nature of the algorithm used.

Here is an example with the trivial game where all payoffs are 0::

    sage: g = NormalFormGame()
    sage: g.add_player(3)  # Adding first player with 3 strategies
    sage: g.add_player(3)  # Adding second player with 3 strategies
    sage: for key in g:
    ....:     g[key] = [0, 0]
    sage: g.payoff_matrices()
    (
    [0 0 0]  [0 0 0]
    [0 0 0]  [0 0 0]
    [0 0 0], [0 0 0]
    )
    sage: g.obtain_nash(algorithm='enumeration')
    [[(0, 0, 1), (0, 0, 1)], [(0, 0, 1), (0, 1, 0)], [(0, 0, 1), (1, 0, 0)],
     [(0, 1, 0), (0, 0, 1)], [(0, 1, 0), (0, 1, 0)], [(0, 1, 0), (1, 0, 0)],
     [(1, 0, 0), (0, 0, 1)], [(1, 0, 0), (0, 1, 0)], [(1, 0, 0), (1, 0, 0)]]

A good description of degenerate games can be found in [NN2007]_.

REFERENCES:

- [Nas1950]_

- [NN2007]_

- [Av2000]_

- [Gambit]_

- [SLB2008]_

AUTHORS:

- James Campbell and Vince Knight (06-2014): Original version

- Tobenna P. Igwe: Constant-sum game solvers
"""

# ****************************************************************************
#       Copyright (C) 2014 James Campbell james.campbell@tanti.org.uk
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from collections.abc import MutableMapping
from itertools import product
from .parser import Parser
from sage.misc.latex import latex
from sage.combinat.subset import powerset
from sage.rings.rational_field import QQ
from sage.structure.sage_object import SageObject
from sage.matrix.constructor import matrix
from sage.matrix.constructor import vector
from sage.misc.temporary_file import tmp_filename
from sage.numerical.mip import MixedIntegerLinearProgram
from sage.cpython.string import bytes_to_str

try:
    from gambit import Game
    from gambit.nash import ExternalLPSolver, ExternalLCPSolver
except ImportError:
    Game = None
    ExternalLPSolver = None
    ExternalLCPSolver = None


class NormalFormGame(SageObject, MutableMapping):
    r"""
    An object representing a Normal Form Game. Primarily used to compute the
    Nash Equilibria.

    INPUT:

    - ``generator`` -- can be a list of 2 matrices, a single matrix or left
      blank
    """

    def __init__(self, generator=None):
        r"""
        Initialize a Normal Form game and checks the inputs.

        EXAMPLES:

        Can have games with more than 2 players::

            sage: threegame = NormalFormGame()
            sage: threegame.add_player(2)  # Adding first player with 2 strategies
            sage: threegame.add_player(2)  # Adding second player with 2 strategies
            sage: threegame.add_player(2)  # Adding third player with 2 strategies
            sage: threegame[0, 0, 0][0] = 3
            sage: threegame[0, 0, 0][1] = 1
            sage: threegame[0, 0, 0][2] = 4
            sage: threegame[0, 0, 1][0] = 1
            sage: threegame[0, 0, 1][1] = 5
            sage: threegame[0, 0, 1][2] = 9
            sage: threegame[0, 1, 0][0] = 2
            sage: threegame[0, 1, 0][1] = 6
            sage: threegame[0, 1, 0][2] = 5
            sage: threegame[0, 1, 1][0] = 3
            sage: threegame[0, 1, 1][1] = 5
            sage: threegame[0, 1, 1][2] = 8
            sage: threegame[1, 0, 0][0] = 9
            sage: threegame[1, 0, 0][1] = 7
            sage: threegame[1, 0, 0][2] = 9
            sage: threegame[1, 0, 1][0] = 3
            sage: threegame[1, 0, 1][1] = 2
            sage: threegame[1, 0, 1][2] = 3
            sage: threegame[1, 1, 0][0] = 8
            sage: threegame[1, 1, 0][1] = 4
            sage: threegame[1, 1, 0][2] = 6
            sage: threegame[1, 1, 1][0] = 2
            sage: threegame[1, 1, 1][1] = 6
            sage: threegame[1, 1, 1][2] = 4
            sage: threegame.obtain_nash()
            Traceback (most recent call last):
            ...
            NotImplementedError: Nash equilibrium for games with more than
             2 players have not been implemented yet. Please see the gambit
             website (http://gambit.sourceforge.net/) that has a variety of
             available algorithms

        Can initialise a game from a gambit game object::

            sage: # optional - gambit
            sage: from gambit import Game
            sage: gambitgame= Game.new_table([2, 2])
            sage: gambitgame[int(0), int(0)][int(0)] = int(5)
            sage: gambitgame[int(0), int(0)][int(1)] = int(8)
            sage: gambitgame[int(0), int(1)][int(0)] = int(2)
            sage: gambitgame[int(0), int(1)][int(1)] = int(11)
            sage: gambitgame[int(1), int(0)][int(0)] = int(10)
            sage: gambitgame[int(1), int(0)][int(1)] = int(7)
            sage: gambitgame[int(1), int(1)][int(0)] = int(5)
            sage: gambitgame[int(1), int(1)][int(1)] = int(5)
            sage: g = NormalFormGame(gambitgame); g
            Normal Form Game with the following utilities: {(0, 0): [5.0, 8.0],
             (0, 1): [2.0, 11.0],
             (1, 0): [10.0, 7.0],
             (1, 1): [5.0, 5.0]}

        TESTS:

        Raise error if matrices aren't the same size::

            sage: p1 = matrix([[1, 2], [3, 4]])
            sage: p2 = matrix([[3, 3], [1, 4], [6, 6]])
            sage: error = NormalFormGame([p1, p2])
            Traceback (most recent call last):
            ...
            ValueError: matrices must be the same size

        Note that when initializing, a single argument must be passed::

            sage: p1 = matrix([[1, 2], [3, 4]])
            sage: p2 = matrix([[3, 3], [1, 4], [6, 6]])
            sage: error = NormalFormGame(p1, p2)
            Traceback (most recent call last):
            ...
            TypeError: ...__init__() takes from 1 to 2 positional arguments but 3 were given

        When initiating, argument passed must be a list or nothing::

            sage: error = NormalFormGame({4:6, 6:9})
            Traceback (most recent call last):
            ...
            TypeError: Generator function must be a list, gambit game or nothing

        When passing nothing, the utilities then need to be entered manually::

            sage: game = NormalFormGame()
            sage: game
            Normal Form Game with the following utilities: {}
        """
        self.players = []
        self.utilities = {}
        matrices = []
        if generator is not None:
            if type(generator) is not list and type(generator) is not Game:
                raise TypeError("Generator function must be a list, gambit game or nothing")

        if type(generator) is list:
            if len(generator) == 1:
                generator.append(-generator[-1])
            matrices = generator
            if matrices[0].dimensions() != matrices[1].dimensions():
                raise ValueError("matrices must be the same size")
            self._two_matrix_game(matrices)
        elif type(generator) is Game:
            game = generator
            self._gambit_game(game)

    def __delitem__(self, key):
        r"""
        This method is one of a collection that aims to make a game
        instance behave like a dictionary which can be used if a game
        is to be generated without using a matrix.

        Here we set up deleting an element of the utilities dictionary::

            sage: A = matrix([[2, 5], [0, 4]])
            sage: B = matrix([[2, 0], [5, 4]])
            sage: prisoners_dilemma = NormalFormGame([A, B])
            sage: prisoners_dilemma
            Normal Form Game with the following utilities: {(0, 0): [2, 2],
             (0, 1): [5, 0], (1, 0): [0, 5], (1, 1): [4, 4]}
            sage: del(prisoners_dilemma[(0,1)])
            sage: prisoners_dilemma
            Normal Form Game with the following utilities: {(0, 0): [2, 2],
             (1, 0): [0, 5], (1, 1): [4, 4]}
        """
        self.utilities.pop(key, None)

    def __getitem__(self, key):
        r"""
        This method is one of a collection that aims to make a game
        instance behave like a dictionary which can be used if a game
        is to be generated without using a matrix.

        Here we allow for querying a key::

            sage: A = matrix([[2, 5], [0, 4]])
            sage: B = matrix([[2, 0], [5, 4]])
            sage: prisoners_dilemma = NormalFormGame([A, B])
            sage: prisoners_dilemma[(0, 1)]
            [5, 0]
            sage: del(prisoners_dilemma[(0,1)])
            sage: prisoners_dilemma[(0, 1)]
            Traceback (most recent call last):
            ...
            KeyError: (0, 1)
        """

        return self.utilities[key]

    def __iter__(self):
        r"""
        This method is one of a collection that aims to make a game
        instance behave like a dictionary which can be used if a game
        is to be generated without using a matrix.

        Here we allow for iteration over the game to correspond to
        iteration over keys of the utility dictionary::

            sage: A = matrix([[2, 5], [0, 4]])
            sage: B = matrix([[2, 0], [5, 4]])
            sage: prisoners_dilemma = NormalFormGame([A, B])
            sage: for key, value in sorted(prisoners_dilemma.items()):
            ....:     print("The strategy pair {} gives utilities {}".format(key, value))
            The strategy pair (0, 0) gives utilities [2, 2]
            The strategy pair (0, 1) gives utilities [5, 0]
            The strategy pair (1, 0) gives utilities [0, 5]
            The strategy pair (1, 1) gives utilities [4, 4]
        """
        return iter(self.utilities)

    def __setitem__(self, key, value):
        r"""
        This method is one of a collection that aims to make a game
        instance behave like a dictionary which can be used if a game
        is to be generated without using a matrix.

        Here we set up setting the value of a key::

            sage: A = matrix([[2, 5], [0, 4]])
            sage: B = matrix([[2, 0], [5, 4]])
            sage: prisoners_dilemma = NormalFormGame([A, B])
            sage: del(prisoners_dilemma[(0,1)])
            sage: prisoners_dilemma[(0,1)] = [5,6]
            sage: prisoners_dilemma.payoff_matrices()
            (
            [2 5]  [2 6]
            [0 4], [5 4]
            )

        We can use the dictionary-like interface to overwrite a strategy
        profile::

            sage: prisoners_dilemma[(0,1)] = [-3,-30]
            sage: prisoners_dilemma.payoff_matrices()
            (
            [ 2 -3]  [  2 -30]
            [ 0  4], [  5   4]
            )
        """
        self.utilities[key] = value

    def __len__(self):
        r"""
        Return the length of the game to be the length of the utilities.

        EXAMPLES::

            sage: A = matrix([[2, 5], [0, 4]])
            sage: B = matrix([[2, 0], [5, 4]])
            sage: prisoners_dilemma = NormalFormGame([A, B])
            sage: len(prisoners_dilemma)
            4
        """
        return len(self.utilities)

    def _repr_(self):
        r"""
        Return the strategy_profiles of the game.

        EXAMPLES:

        Basic description of the game shown when calling the game instance::

            sage: p1 = matrix([[1, 2], [3, 4]])
            sage: p2 = matrix([[3, 3], [1, 4]])
            sage: g = NormalFormGame([p1, p2])
            sage: g
            Normal Form Game with the following utilities: {(0, 0): [1, 3],
             (0, 1): [2, 3], (1, 0): [3, 1], (1, 1): [4, 4]}
        """
        from pprint import pformat
        base_str = "Normal Form Game with the following utilities: {}"
        return base_str.format(pformat(self.utilities))

    def _latex_(self):
        r"""
        Return the LaTeX code representing the ``NormalFormGame``.

        EXAMPLES:

        LaTeX method shows the two payoff matrices for a two player game::

            sage: A = matrix([[-1, -2], [-12, 2]])
            sage: B = matrix([[1, 0], [1, -1]])
            sage: g = NormalFormGame([A, B])
            sage: latex(g)
            \left(\left(\begin{array}{rr}
            -1 & -2 \\
            -12 & 2
            \end{array}\right), \left(\begin{array}{rr}
            1 & 0 \\
            1 & -1
            \end{array}\right)\right)

        LaTeX method shows nothing interesting for games with more players::

            sage: g = NormalFormGame()
            sage: g.add_player(2)  # Adding first player with 2 strategies
            sage: g.add_player(2)  # Adding second player with 2 strategies
            sage: g.add_player(2)  # Creating a game with three players
            sage: latex(g)
            \begin{array}{l}
            \text{\texttt{Normal{ }Form{ }Game{ }with{ }the{ }...
            ...
            \end{array}
        """
        if len(self.players) == 2:
            M1, M2 = self.payoff_matrices()
            return r"\left(%s, %s\right)" % (M1._latex_(), M2._latex_())
        return latex(str(self))

    def _two_matrix_game(self, matrices):
        r"""
        Populate ``self.utilities`` with the values from 2 matrices.

        EXAMPLES:

        A small example game::

            sage: A = matrix([[1, 0], [-2, 3]])
            sage: B = matrix([[3, 2], [-1, 0]])
            sage: two_game = NormalFormGame()
            sage: two_game._two_matrix_game([A, B])
        """
        self.players = []
        self.utilities = {}
        self.add_player(matrices[0].dimensions()[0])
        self.add_player(matrices[1].dimensions()[1])
        for strategy_profile in self.utilities:
            self.utilities[strategy_profile] = [matrices[0][strategy_profile],
                                                matrices[1][strategy_profile]]

    def _gambit_game(self, game):
        r"""
        Create a ``NormalFormGame`` object from a Gambit game.

        TESTS::

            sage: # optional - gambit
            sage: from gambit import Game
            sage: testgame = Game.new_table([2, 2])
            sage: testgame[int(0), int(0)][int(0)] = int(8)
            sage: testgame[int(0), int(0)][int(1)] = int(8)
            sage: testgame[int(0), int(1)][int(0)] = int(2)
            sage: testgame[int(0), int(1)][int(1)] = int(10)
            sage: testgame[int(1), int(0)][int(0)] = int(10)
            sage: testgame[int(1), int(0)][int(1)] = int(2)
            sage: testgame[int(1), int(1)][int(0)] = int(5)
            sage: testgame[int(1), int(1)][int(1)] = int(5)
            sage: g = NormalFormGame()
            sage: g._gambit_game(testgame); g
            Normal Form Game with the following utilities: {(0, 0): [8.0, 8.0],
            (0, 1): [2.0, 10.0],
            (1, 0): [10.0, 2.0],
            (1, 1): [5.0, 5.0]}
        """
        self.players = []
        self.utilities = {}
        for player in game.players:
            num_strategies = len(player.strategies)
            self.add_player(num_strategies)
        for strategy_profile in self.utilities:
            utility_vector = [float(game[strategy_profile][i]) for i in range(len(self.players))]
            self.utilities[strategy_profile] = utility_vector

    def _gambit_(self, as_integer=False, maximization=True):
        r"""
        Create a Gambit game from a ``NormalFormGame`` object.

        INPUT:

        - ``as_integer`` -- boolean; whether the gambit representation
          should have the payoffs represented as integers or decimals

        - ``maximization`` -- boolean; whether a player is trying to
          maximize their utility or minimize it

        TESTS::

            sage: # optional - gambit
            sage: from gambit import Game
            sage: A = matrix([[2, 1], [1, 2.5]])
            sage: g = NormalFormGame([A])
            sage: gg = g._gambit_(); gg
            NFG 1 R "" { "1" "2" }
            <BLANKLINE>
            { { "1" "2" }
            { "1" "2" }
            }
            ""
            <BLANKLINE>
            {
            { "" 2, -2 }
            { "" 1, -1 }
            { "" 1, -1 }
            { "" 2.5, -2.5 }
            }
            1 2 3 4
            <BLANKLINE>
            sage: gg = g._gambit_(as_integer=True); gg
            NFG 1 R "" { "1" "2" }
            <BLANKLINE>
            { { "1" "2" }
            { "1" "2" }
            }
            ""
            <BLANKLINE>
            {
            { "" 2, -2 }
            { "" 1, -1 }
            { "" 1, -1 }
            { "" 2, -2 }
            }
            1 2 3 4
            <BLANKLINE>

        ::

            sage: # optional - gambit
            sage: A = matrix([[2, 1], [1, 2.5]])
            sage: B = matrix([[3, 2], [5.5, 4]])
            sage: g = NormalFormGame([A, B])
            sage: gg = g._gambit_(); gg
            NFG 1 R "" { "1" "2" }
            <BLANKLINE>
            { { "1" "2" }
            { "1" "2" }
            }
            ""
            <BLANKLINE>
            {
            { "" 2, 3 }
            { "" 1, 5.5 }
            { "" 1, 2 }
            { "" 2.5, 4 }
            }
            1 2 3 4
            <BLANKLINE>
            sage: gg = g._gambit_(as_integer = True); gg
            NFG 1 R "" { "1" "2" }
            <BLANKLINE>
            { { "1" "2" }
            { "1" "2" }
            }
            ""
            <BLANKLINE>
            {
            { "" 2, 3 }
            { "" 1, 5 }
            { "" 1, 2 }
            { "" 2, 4 }
            }
            1 2 3 4
            <BLANKLINE>

        ::

            sage: # optional - gambit
            sage: threegame = NormalFormGame()
            sage: threegame.add_player(2)
            sage: threegame.add_player(2)
            sage: threegame.add_player(2)
            sage: threegame[0, 0, 0][0] = 3
            sage: threegame[0, 0, 0][1] = 1
            sage: threegame[0, 0, 0][2] = 4
            sage: threegame[0, 0, 1][0] = 1
            sage: threegame[0, 0, 1][1] = 5
            sage: threegame[0, 0, 1][2] = 9
            sage: threegame[0, 1, 0][0] = 2
            sage: threegame[0, 1, 0][1] = 6
            sage: threegame[0, 1, 0][2] = 5
            sage: threegame[0, 1, 1][0] = 3
            sage: threegame[0, 1, 1][1] = 5
            sage: threegame[0, 1, 1][2] = 8
            sage: threegame[1, 0, 0][0] = 9
            sage: threegame[1, 0, 0][1] = 7
            sage: threegame[1, 0, 0][2] = 9
            sage: threegame[1, 0, 1][0] = 3
            sage: threegame[1, 0, 1][1] = 2
            sage: threegame[1, 0, 1][2] = 3
            sage: threegame[1, 1, 0][0] = 8
            sage: threegame[1, 1, 0][1] = 4
            sage: threegame[1, 1, 0][2] = 6
            sage: threegame[1, 1, 1][0] = 2
            sage: threegame[1, 1, 1][1] = 6
            sage: threegame[1, 1, 1][2] = 4
            sage: threegame._gambit_(as_integer = True)
            NFG 1 R "" { "1" "2" "3" }
            <BLANKLINE>
            { { "1" "2" }
            { "1" "2" }
            { "1" "2" }
            }
            ""
            <BLANKLINE>
            {
            { "" 3, 1, 4 }
            { "" 9, 7, 9 }
            { "" 2, 6, 5 }
            { "" 8, 4, 6 }
            { "" 1, 5, 9 }
            { "" 3, 2, 3 }
            { "" 3, 5, 8 }
            { "" 2, 6, 4 }
            }
            1 2 3 4 5 6 7 8
            <BLANKLINE>
        """
        from decimal import Decimal
        strategy_sizes = [p.num_strategies for p in self.players]
        g = Game.new_table(strategy_sizes)

        sgn = 1
        if not maximization:
            sgn = -1

        players = len(strategy_sizes)

        for strategy_profile in self.utilities:
            for i in range(players):
                if as_integer:
                    g[strategy_profile][i] = sgn * int(self.utilities[strategy_profile][i])
                else:
                    g[strategy_profile][i] = sgn * Decimal(float(self.utilities[strategy_profile][i]))
        return g

    def is_constant_sum(self):
        r"""
        Check if the game is constant sum.

        EXAMPLES::

            sage: A = matrix([[2, 1], [1, 2.5]])
            sage: g = NormalFormGame([A])
            sage: g.is_constant_sum()
            True
            sage: g = NormalFormGame([A, A])
            sage: g.is_constant_sum()
            False
            sage: A = matrix([[1, 1], [1, 1]])
            sage: g = NormalFormGame([A, A])
            sage: g.is_constant_sum()
            True
            sage: A = matrix([[1, 1, 2], [1, 1, -1], [1, -1, 1]])
            sage: B = matrix([[2, 2, 1], [2, 2, 4], [2, 4, 2]])
            sage: g = NormalFormGame([A, B])
            sage: g.is_constant_sum()
            True
            sage: A = matrix([[1, 1, 2], [1, 1, -1], [1, -1, 1]])
            sage: B = matrix([[2, 2, 1], [2, 2.1, 4], [2, 4, 2]])
            sage: g = NormalFormGame([A, B])
            sage: g.is_constant_sum()
            False
        """
        import sys
        if len(self.players) > 2:
            return False
        m1, m2 = self.payoff_matrices()
        c = m1 + m2
        t = c[0,0]

        for row in c:
            for i in row:
                if abs(t - i) > sys.float_info.epsilon:
                    return False

        return True

    def payoff_matrices(self):
        r"""
        Return 2 matrices representing the payoffs for each player.

        EXAMPLES::

            sage: p1 = matrix([[1, 2], [3, 4]])
            sage: p2 = matrix([[3, 3], [1, 4]])
            sage: g = NormalFormGame([p1, p2])
            sage: g.payoff_matrices()
            (
            [1 2]  [3 3]
            [3 4], [1 4]
            )

        If we create a game with 3 players we will not be able to
        obtain payoff matrices::

            sage: g = NormalFormGame()
            sage: g.add_player(2)  # adding first player with 2 strategies
            sage: g.add_player(2)  # adding second player with 2 strategies
            sage: g.add_player(2)  # adding third player with 2 strategies
            sage: g.payoff_matrices()
            Traceback (most recent call last):
            ...
            ValueError: Only available for 2 player games

        If we do create a two player game but it is not complete
        then an error is also raised::

            sage: g = NormalFormGame()
            sage: g.add_player(1)  # Adding first player with 1 strategy
            sage: g.add_player(1)  # Adding second player with 1 strategy
            sage: g.payoff_matrices()
            Traceback (most recent call last):
            ...
            ValueError: utilities have not been populated

        The above creates a 2 player game where each player has
        a single strategy. Here we populate the strategies and
        can then view the payoff matrices::

            sage: g[0, 0] = [1,2]
            sage: g.payoff_matrices()
            ([1], [2])
        """
        if len(self.players) != 2:
            raise ValueError("Only available for 2 player games")

        if not self._is_complete():
            raise ValueError("utilities have not been populated")

        m1 = matrix(QQ, self.players[0].num_strategies, self.players[1].num_strategies)
        m2 = matrix(QQ, self.players[0].num_strategies, self.players[1].num_strategies)
        for strategy_profile in self.utilities:
            m1[strategy_profile] = self[strategy_profile][0]
            m2[strategy_profile] = self[strategy_profile][1]
        return m1, m2

    def add_player(self, num_strategies):
        r"""
        Add a player to a NormalFormGame.

        INPUT:

        - ``num_strategies`` -- the number of strategies the player should have

        EXAMPLES::

            sage: g = NormalFormGame()
            sage: g.add_player(2)  # Adding first player with 2 strategies
            sage: g.add_player(1)  # Adding second player with 1 strategy
            sage: g.add_player(1)  # Adding third player with 1 strategy
            sage: g
            Normal Form Game with the following utilities:
             {(0, 0, 0): [False, False, False],
              (1, 0, 0): [False, False, False]}
        """
        self.players.append(_Player(num_strategies))
        self._generate_utilities(True)

    def _generate_utilities(self, replacement):
        r"""
        Create all the required keys for ``self.utilities``.

        This is used when generating players and/or adding strategies.

        INPUT:

        - ``replacement`` -- boolean value of whether previously created
          profiles should be replaced or not

        TESTS::

            sage: from sage.game_theory.normal_form_game import _Player
            sage: g = NormalFormGame()
            sage: g.players.append(_Player(2))
            sage: g.players.append(_Player(2))
            sage: g
            Normal Form Game with the following utilities: {}

            sage: g._generate_utilities(True)
            sage: g
            Normal Form Game with the following utilities: {(0, 0): [False, False],
             (0, 1): [False, False],
             (1, 0): [False, False],
             (1, 1): [False, False]}

            sage: g[(0,1)] = [2, 3]
            sage: g.add_strategy(1)
            sage: g._generate_utilities(False)
            sage: g
            Normal Form Game with the following utilities: {(0, 0): [False, False],
             (0, 1): [2, 3],
             (0, 2): [False, False],
             (1, 0): [False, False],
             (1, 1): [False, False],
             (1, 2): [False, False]}

            sage: g._generate_utilities(True)
            sage: g
            Normal Form Game with the following utilities: {(0, 0): [False, False],
             (0, 1): [False, False],
             (0, 2): [False, False],
             (1, 0): [False, False],
             (1, 1): [False, False],
             (1, 2): [False, False]}
        """
        strategy_sizes = [range(p.num_strategies) for p in self.players]
        if replacement is True:
            self.utilities = {}
        for profile in product(*strategy_sizes):
            if profile not in self.utilities.keys():
                self.utilities[profile] = [False] * len(self.players)

    def add_strategy(self, player):
        r"""
        Add a strategy to a player, will not affect already completed
        strategy profiles.

        INPUT:

        - ``player`` -- the index of the player

        EXAMPLES:

        A simple example::

            sage: s = matrix([[1, 0], [-2, 3]])
            sage: t = matrix([[3, 2], [-1, 0]])
            sage: example = NormalFormGame([s, t])
            sage: example
            Normal Form Game with the following utilities: {(0, 0): [1, 3],
             (0, 1): [0, 2], (1, 0): [-2, -1], (1, 1): [3, 0]}
            sage: example.add_strategy(0)
            sage: example
            Normal Form Game with the following utilities: {(0, 0): [1, 3],
             (0, 1): [0, 2],
             (1, 0): [-2, -1],
             (1, 1): [3, 0],
             (2, 0): [False, False],
             (2, 1): [False, False]}
        """
        self.players[player].add_strategy()
        self._generate_utilities(False)

    def _is_complete(self):
        r"""
        Check if ``utilities`` has been completed and return a
        boolean.

        EXAMPLES:

        A simple example::

            sage: s = matrix([[1, 0], [-2, 3]])
            sage: t = matrix([[3, 2], [-1, 0]])
            sage: example = NormalFormGame([s, t])
            sage: example.add_strategy(0)
            sage: example._is_complete()
            False
        """
        results = (all(not isinstance(i, bool) for i in profile)
                   for profile in self.utilities.values())
        return all(results)

    def obtain_nash(self, algorithm=False, maximization=True, solver=None):
        r"""
        A function to return the Nash equilibrium for the game.
        Optional arguments can be used to specify the algorithm used.
        If no algorithm is passed then an attempt is made to use the most
        appropriate algorithm.

        INPUT:

        - ``algorithm`` -- the following algorithms should be available through
          this function:

          * ``'lrs'`` -- this algorithm is only suited for 2 player games.
            See the lrs web site (http://cgm.cs.mcgill.ca/~avis/C/lrs.html).

          * ``'LCP'`` -- this algorithm is only suited for 2 player games.
            See the gambit web site (http://gambit.sourceforge.net/).

          * ``'lp'`` -- this algorithm is only suited for 2 player
            constant sum games. Uses MILP solver determined by the
            ``solver`` argument.

          * ``'enumeration'`` -- this is a very inefficient
            algorithm (in essence a brute force approach).

            1. For each k in 1...min(size of strategy sets)
            2. For each I,J supports of size k
            3. Prune: check if supports are dominated
            4. Solve indifference conditions and check that have Nash Equilibrium.

            Solving the indifference conditions is done by building the
            corresponding linear system.  If  `\rho_1, \rho_2` are the
            supports player 1 and 2 respectively.  Then, indifference implies:

            .. MATH::

                u_1(s_1,\rho_2) = u_1(s_2, \rho_2)

            for all `s_1, s_2` in the support of `\rho_1`. This corresponds to:

            .. MATH::

                \sum_{j\in S(\rho_2)}A_{s_1,j}{\rho_2}_j = \sum_{j\in S(\rho_2)}A_{s_2,j}{\rho_2}_j

            for all `s_1, s_2` in the support of `\rho_1` where `A` is the payoff
            matrix of player 1. Equivalently we can consider consecutive rows of
            `A` (instead of all pairs of strategies). Thus the corresponding
            linear system can be written as:

            .. MATH::

                \left(\sum_{j \in S(\rho_2)}A_{i,j} - A_{i+1,j}\right){\rho_2}_j

            for all `1\leq i \leq |S(\rho_1)|` (where `A` has been modified to only
            contain the rows corresponding to `S(\rho_1)`). We also require all
            elements of `\rho_2` to sum to 1:

            .. MATH::

                \sum_{j\in S(\rho_1)}{\rho_2}_j = 1

        - ``maximization`` -- boolean (default: ``True``); whether a player is
          trying to maximize their utility or minimize it:

          * When set to ``True`` it is assumed that players aim to
            maximise their utility.

          * When set to ``False`` it is assumed that players aim to
            minimise their utility.

        - ``solver`` -- (optional) see :class:`MixedIntegerLinearProgram`
          for more information on the MILP solvers in Sage, may also
          be ``'gambit'`` to use the MILP solver included with the gambit
          library. Note that ``None`` means to use the default Sage LP solver,
          normally GLPK.

        EXAMPLES:

        A game with 1 equilibrium when ``maximization`` is ``True`` and 3 when
        ``maximization`` is ``False``::

            sage: A = matrix([[10, 500, 44],
            ....:       [15, 10, 105],
            ....:       [19, 204, 55],
            ....:       [20, 200, 590]])
            sage: B = matrix([[2, 1, 2],
            ....:             [0, 5, 6],
            ....:             [3, 4, 1],
            ....:             [4, 1, 20]])
            sage: g=NormalFormGame([A, B])
            sage: g.obtain_nash(algorithm='lrs')  # optional - lrslib
            [[(0, 0, 0, 1), (0, 0, 1)]]
            sage: g.obtain_nash(algorithm='lrs', maximization=False)  # optional - lrslib
            [[(2/3, 1/12, 1/4, 0), (6333/8045, 247/8045, 293/1609)],
             [(3/4, 0, 1/4, 0), (0, 11/307, 296/307)],
             [(5/6, 1/6, 0, 0), (98/99, 1/99, 0)]]

        This particular game has 3 Nash equilibria::

            sage: A = matrix([[3,3],
            ....:             [2,5],
            ....:             [0,6]])
            sage: B = matrix([[3,2],
            ....:             [2,6],
            ....:             [3,1]])
            sage: g = NormalFormGame([A, B])
            sage: g.obtain_nash(algorithm='enumeration')
            [[(0, 1/3, 2/3), (1/3, 2/3)],
             [(4/5, 1/5, 0), (2/3, 1/3)],
             [(1, 0, 0), (1, 0)]]

        Here is a slightly larger game::

            sage: A = matrix([[160, 205, 44],
            ....:             [175, 180, 45],
            ....:             [201, 204, 50],
            ....:             [120, 207, 49]])
            sage: B = matrix([[2, 2, 2],
            ....:             [1, 0, 0],
            ....:             [3, 4, 1],
            ....:             [4, 1, 2]])
            sage: g=NormalFormGame([A, B])
            sage: g.obtain_nash(algorithm='enumeration')
            [[(0, 0, 3/4, 1/4), (1/28, 27/28, 0)]]
            sage: g.obtain_nash(algorithm='lrs')  # optional - lrslib
            [[(0, 0, 3/4, 1/4), (1/28, 27/28, 0)]]
            sage: g.obtain_nash(algorithm='LCP')  # optional - gambit
            [[(0.0, 0.0, 0.75, 0.25), (0.0357142857, 0.9642857143, 0.0)]]

        2 random matrices::

            sage: player1 = matrix([[2, 8, -1, 1, 0],
            ....:                   [1, 1, 2, 1, 80],
            ....:                   [0, 2, 15, 0, -12],
            ....:                   [-2, -2, 1, -20, -1],
            ....:                   [1, -2, -1, -2, 1]])
            sage: player2 = matrix([[0, 8, 4, 2, -1],
            ....:                   [6, 14, -5, 1, 0],
            ....:                   [0, -2, -1, 8, -1],
            ....:                   [1, -1, 3, -3, 2],
            ....:                   [8, -4, 1, 1, -17]])
            sage: fivegame = NormalFormGame([player1, player2])
            sage: fivegame.obtain_nash(algorithm='enumeration')
            [[(1, 0, 0, 0, 0), (0, 1, 0, 0, 0)]]
            sage: fivegame.obtain_nash(algorithm='lrs')  # optional - lrslib
            [[(1, 0, 0, 0, 0), (0, 1, 0, 0, 0)]]
            sage: fivegame.obtain_nash(algorithm='LCP')  # optional - gambit
            [[(1.0, 0.0, 0.0, 0.0, 0.0), (0.0, 1.0, 0.0, 0.0, 0.0)]]

        Here are some examples of finding Nash equilibria for constant-sum games::

            sage: A = matrix.identity(2)
            sage: cg = NormalFormGame([A])
            sage: cg.obtain_nash(algorithm='lp')
            [[(0.5, 0.5), (0.5, 0.5)]]
            sage: cg.obtain_nash(algorithm='lp', solver='Coin')                   # optional - sage_numerical_backends_coin
            [[(0.5, 0.5), (0.5, 0.5)]]
            sage: cg.obtain_nash(algorithm='lp', solver='PPL')
            [[(1/2, 1/2), (1/2, 1/2)]]
            sage: cg.obtain_nash(algorithm='lp', solver='gambit')                 # optional - gambit
            [[(0.5, 0.5), (0.5, 0.5)]]
            sage: A = matrix([[2, 1], [1, 3]])
            sage: cg = NormalFormGame([A])
            sage: ne = cg.obtain_nash(algorithm='lp', solver='glpk')
            sage: [[[round(el, 6) for el in v] for v in eq] for eq in ne]
            [[[0.666667, 0.333333], [0.666667, 0.333333]]]
            sage: ne = cg.obtain_nash(algorithm='lp', solver='Coin')              # optional - sage_numerical_backends_coin
            sage: [[[round(el, 6) for el in v] for v in eq] for eq in ne]         # optional - sage_numerical_backends_coin
            [[[0.666667, 0.333333], [0.666667, 0.333333]]]
            sage: cg.obtain_nash(algorithm='lp', solver='PPL')
            [[(2/3, 1/3), (2/3, 1/3)]]
            sage: ne = cg.obtain_nash(algorithm='lp', solver='gambit')            # optional - gambit
            sage: [[[round(el, 6) for el in v] for v in eq] for eq in ne]         # optional - gambit
            [[[0.666667, 0.333333], [0.666667, 0.333333]]]
            sage: A = matrix([[1, 2, 1], [1, 1, 2], [2, 1, 1]])
            sage: B = matrix([[2, 1, 2], [2, 2, 1], [1, 2, 2]])
            sage: cg = NormalFormGame([A, B])
            sage: ne = cg.obtain_nash(algorithm='lp', solver='glpk')
            sage: [[[round(el, 6) for el in v] for v in eq] for eq in ne]
            [[[0.333333, 0.333333, 0.333333], [0.333333, 0.333333, 0.333333]]]
            sage: ne = cg.obtain_nash(algorithm='lp', solver='Coin')              # optional - sage_numerical_backends_coin
            sage: [[[round(el, 6) for el in v] for v in eq] for eq in ne]         # optional - sage_numerical_backends_coin
            [[[0.333333, 0.333333, 0.333333], [0.333333, 0.333333, 0.333333]]]
            sage: cg.obtain_nash(algorithm='lp', solver='PPL')
            [[(1/3, 1/3, 1/3), (1/3, 1/3, 1/3)]]
            sage: ne = cg.obtain_nash(algorithm='lp', solver='gambit')            # optional - gambit
            sage: [[[round(el, 6) for el in v] for v in eq] for eq in ne]         # optional - gambit
            [[[0.333333, 0.333333, 0.333333], [0.333333, 0.333333, 0.333333]]]
            sage: A = matrix([[160, 205, 44],
            ....:             [175, 180, 45],
            ....:             [201, 204, 50],
            ....:             [120, 207, 49]])
            sage: cg = NormalFormGame([A])
            sage: cg.obtain_nash(algorithm='lp', solver='PPL')
            [[(0, 0, 1, 0), (0, 0, 1)]]

        Running the constant-sum solver on a game which is not a constant sum
        game generates a :exc:`ValueError`::

            sage: cg = NormalFormGame([A, A])
            sage: cg.obtain_nash(algorithm='lp', solver='glpk')
            Traceback (most recent call last):
            ...
            ValueError: Input game needs to be a two player constant sum game

        Here is an example of a 3 by 2 game with 3 Nash equilibrium::

            sage: A = matrix([[3,3],
            ....:             [2,5],
            ....:             [0,6]])
            sage: B = matrix([[3,2],
            ....:             [2,6],
            ....:             [3,1]])
            sage: g = NormalFormGame([A, B])
            sage: g.obtain_nash(algorithm='enumeration')
            [[(0, 1/3, 2/3), (1/3, 2/3)], [(4/5, 1/5, 0), (2/3, 1/3)], [(1, 0, 0), (1, 0)]]

        Of the algorithms implemented, only ``'lrs'`` and ``'enumeration'``
        are guaranteed to find all Nash equilibria in a game. The solver for
        constant sum games only ever finds one Nash equilibrium. Although it
        is possible for the ``'LCP'`` solver to find all Nash equilibria
        in some instances, there are instances where it will not be able to
        find all Nash equilibria.::

            sage: A = matrix(2, 2)
            sage: gg = NormalFormGame([A])
            sage: gg.obtain_nash(algorithm='enumeration')
            [[(0, 1), (0, 1)], [(0, 1), (1, 0)], [(1, 0), (0, 1)], [(1, 0), (1, 0)]]
            sage: gg.obtain_nash(algorithm='lrs')  # optional - lrs
            [[(0, 1), (0, 1)], [(0, 1), (1, 0)], [(1, 0), (0, 1)], [(1, 0), (1, 0)]]
            sage: gg.obtain_nash(algorithm='lp', solver='glpk')
            [[(1.0, 0.0), (1.0, 0.0)]]
            sage: gg.obtain_nash(algorithm='LCP')  # optional - gambit
            [[(1.0, 0.0), (1.0, 0.0)]]
            sage: gg.obtain_nash(algorithm='enumeration', maximization=False)
            [[(0, 1), (0, 1)], [(0, 1), (1, 0)], [(1, 0), (0, 1)], [(1, 0), (1, 0)]]
            sage: gg.obtain_nash(algorithm='lrs', maximization=False)  # optional - lrs
            [[(0, 1), (0, 1)], [(0, 1), (1, 0)], [(1, 0), (0, 1)], [(1, 0), (1, 0)]]
            sage: gg.obtain_nash(algorithm='lp', solver='glpk', maximization=False)
            [[(1.0, 0.0), (1.0, 0.0)]]
            sage: gg.obtain_nash(algorithm='LCP', maximization=False)  # optional - gambit
            [[(1.0, 0.0), (1.0, 0.0)]]

        Note that outputs for all algorithms are as lists of lists of
        tuples and the equilibria have been sorted so that all algorithms give
        a comparable output (although ``'LCP'`` returns floats)::

            sage: enumeration_eqs = g.obtain_nash(algorithm='enumeration')
            sage: [[type(s) for s in eq] for eq in enumeration_eqs]
            [[<... 'tuple'>, <... 'tuple'>], [<... 'tuple'>, <... 'tuple'>], [<... 'tuple'>, <... 'tuple'>]]
            sage: lrs_eqs = g.obtain_nash(algorithm='lrs')  # optional - lrslib
            sage: [[type(s) for s in eq] for eq in lrs_eqs]  # optional - lrslib
            [[<... 'tuple'>, <... 'tuple'>], [<... 'tuple'>, <... 'tuple'>], [<... 'tuple'>, <... 'tuple'>]]
            sage: LCP_eqs = g.obtain_nash(algorithm='LCP')  # optional - gambit
            sage: [[type(s) for s in eq] for eq in LCP_eqs]  # optional - gambit
            [[<... 'tuple'>, <... 'tuple'>], [<... 'tuple'>, <... 'tuple'>], [<... 'tuple'>, <... 'tuple'>]]
            sage: enumeration_eqs == sorted(enumeration_eqs)
            True
            sage: lrs_eqs == sorted(lrs_eqs)  # optional - lrslib
            True
            sage: LCP_eqs == sorted(LCP_eqs)  # optional - gambit
            True
            sage: lrs_eqs == enumeration_eqs  # optional - lrslib
            True
            sage: enumeration_eqs == LCP_eqs  # optional - gambit
            False
            sage: [[[round(float(p), 6) for p in str] for str in eq] for eq in enumeration_eqs] == [[[round(float(p), 6) for p in str] for str in eq] for eq in LCP_eqs]  # optional - gambit
            True

        Also, not specifying a valid solver would lead to an error::

            sage: A = matrix.identity(2)
            sage: g = NormalFormGame([A])
            sage: g.obtain_nash(algorithm='invalid')
            Traceback (most recent call last):
            ...
            ValueError: 'algorithm' should be set to 'enumeration', 'LCP', 'lp' or 'lrs'
            sage: g.obtain_nash(algorithm='lp', solver='invalid')
            Traceback (most recent call last):
            ...
            ValueError: 'solver' should be set to 'GLPK', ..., None
             (in which case the default one is used), or a callable.
        """
        if len(self.players) > 2:
            raise NotImplementedError("Nash equilibrium for games with more "
                                      "than 2 players have not been "
                                      "implemented yet. Please see the gambit "
                                      "website (http://gambit.sourceforge.net/) that has a variety of "
                                      "available algorithms")

        if not self._is_complete():
            raise ValueError("utilities have not been populated")

        from sage.features.lrs import LrsNash
        if not algorithm:
            if self.is_constant_sum():
                algorithm = "lp"
            elif LrsNash().is_present():
                algorithm = "lrs"
            else:
                algorithm = "enumeration"

        if algorithm == "lrs":
            LrsNash().require()
            return self._solve_lrs(maximization)

        if algorithm == "LCP":
            if Game is None:
                raise RuntimeError("gambit not found")  # should later become a FeatureNotFoundError
            return self._solve_LCP(maximization)

        if algorithm.startswith('lp'):
            return self._solve_LP(solver=solver, maximization=maximization)

        if algorithm == "enumeration":
            return self._solve_enumeration(maximization)

        raise ValueError("'algorithm' should be set to 'enumeration', 'LCP', 'lp' or 'lrs'")

    def _solve_lrs(self, maximization=True):
        r"""
        EXAMPLES:

        A simple game::

            sage: A = matrix([[1, 2], [3, 4]])
            sage: B = matrix([[3, 3], [1, 4]])
            sage: C = NormalFormGame([A, B])
            sage: C._solve_lrs()  # optional - lrslib
            [[(0, 1), (0, 1)]]

        2 random matrices::

            sage: p1 = matrix([[-1, 4, 0, 2, 0],
            ....:              [-17, 246, -5, 1, -2],
            ....:              [0, 1, 1, -4, -4],
            ....:              [1, -3, 9, 6, -1],
            ....:              [2, 53, 0, -5, 0]])
            sage: p2 = matrix([[0, 1, 1, 3, 1],
            ....:              [3, 9, 44, -1, -1],
            ....:              [1, -4, -1, -3, 1],
            ....:              [1, 0, 0, 0, 0,],
            ....:              [1, -3, 1, 21, -2]])
            sage: biggame = NormalFormGame([p1, p2])
            sage: biggame._solve_lrs()  # optional - lrslib
            [[(0, 0, 0, 20/21, 1/21), (11/12, 0, 0, 1/12, 0)]]

        Another test::

            sage: p1 = matrix([[-7, -5, 5],
            ....:              [5, 5, 3],
            ....:              [1, -6, 1]])
            sage: p2 = matrix([[-9, 7, 9],
            ....:              [6, -2, -3],
            ....:              [-4, 6, -10]])
            sage: biggame = NormalFormGame([p1, p2])
            sage: biggame._solve_lrs()  # optional - lrslib
            [[(0, 1, 0), (1, 0, 0)],
             [(1/3, 2/3, 0), (0, 1/6, 5/6)],
             [(1/3, 2/3, 0), (1/7, 0, 6/7)],
             [(1, 0, 0), (0, 0, 1)]]
        """
        from subprocess import PIPE, Popen
        m1, m2 = self.payoff_matrices()
        if maximization is False:
            m1 = - m1
            m2 = - m2

        game_str = self._lrs_nash_format(m1, m2)
        game_name = tmp_filename()
        with open(game_name, 'w') as game_file:
            game_file.write(game_str)

        from sage.features.lrs import LrsNash
        LrsNash().require()
        process = Popen([LrsNash().absolute_filename(), game_name],
                        stdout=PIPE, stderr=PIPE)

        lrs_output = [bytes_to_str(row) for row in process.stdout]
        process.terminate()

        nasheq = Parser(lrs_output).format_lrs()
        return sorted(nasheq)

    def _solve_LCP(self, maximization):
        r"""
        Solve a :class:`NormalFormGame` using Gambit's LCP algorithm.

        EXAMPLES::

            sage: a = matrix([[1, 0], [1, 4]])
            sage: b = matrix([[2, 3], [2, 4]])
            sage: c = NormalFormGame([a, b])
            sage: c._solve_LCP(maximization=True)  # optional - gambit
            [[(0.0, 1.0), (0.0, 1.0)]]
        """
        g = self._gambit_(maximization)
        output = ExternalLCPSolver().solve(g)
        nasheq = Parser(output).format_gambit(g)
        return sorted(nasheq)

    def _solve_gambit_LP(self, maximization=True):
        r"""
        Solve a constant sum :class:`NormalFormGame` using Gambit's
        LP implementation.

        EXAMPLES::

            sage: A = matrix([[2, 1], [1, 2.5]])
            sage: g = NormalFormGame([A])
            sage: g._solve_gambit_LP()  # optional - gambit
            [[(0.6, 0.4), (0.6, 0.4)]]
            sage: A = matrix.identity(2)
            sage: g = NormalFormGame([A])
            sage: g._solve_gambit_LP()  # optional - gambit
            [[(0.5, 0.5), (0.5, 0.5)]]
            sage: g = NormalFormGame([A,A])
            sage: g._solve_gambit_LP()  # optional - gambit
            Traceback (most recent call last):
            ...
            RuntimeError: Method only valid for constant-sum games.
        """
        if Game is None:
            raise NotImplementedError("gambit is not installed")
        g = self._gambit_(maximization=maximization)
        output = ExternalLPSolver().solve(g)
        nasheq = Parser(output).format_gambit(g)
        return sorted(nasheq)

    def _solve_LP(self, solver='glpk', maximization=True):
        r"""
        Solve a constant sum :class:`NormalFormGame` using
        the specified LP solver.

        INPUT:

        - ``solver`` -- the solver to be used to solve the LP:

          * ``'gambit'`` -- his uses the solver included within the gambit
            library to create and solve the LP

          * for further possible values, see :class:`MixedIntegerLinearProgram`

        EXAMPLES::

            sage: A = matrix.identity(2)
            sage: g = NormalFormGame([A])
            sage: g._solve_LP()
            [[(0.5, 0.5), (0.5, 0.5)]]
            sage: g._solve_LP('gambit')  # optional - gambit
            [[(0.5, 0.5), (0.5, 0.5)]]
            sage: g._solve_LP('Coin')  # optional - sage_numerical_backends_coin
            [[(0.5, 0.5), (0.5, 0.5)]]
            sage: g._solve_LP('PPL')
            [[(1/2, 1/2), (1/2, 1/2)]]
            sage: A = matrix([[2, 1], [1, 3]])
            sage: g = NormalFormGame([A])
            sage: ne = g._solve_LP()
            sage: [[[round(el, 6) for el in v] for v in eq] for eq in ne]
            [[[0.666667, 0.333333], [0.666667, 0.333333]]]
            sage: ne = g._solve_LP('gambit')  # optional - gambit
            sage: [[[round(el, 6) for el in v] for v in eq] for eq in ne]  # optional - gambit
            [[[0.666667, 0.333333], [0.666667, 0.333333]]]
            sage: ne = g._solve_LP('Coin')  # optional - sage_numerical_backends_coin
            sage: [[[round(el, 6) for el in v] for v in eq] for eq in ne]  # optional - sage_numerical_backends_coin
            [[[0.666667, 0.333333], [0.666667, 0.333333]]]
            sage: g._solve_LP('PPL')
            [[(2/3, 1/3), (2/3, 1/3)]]

        An exception is raised if the input game is not constant sum::

            sage: A = matrix.identity(2)
            sage: B = A.transpose()
            sage: g = NormalFormGame([A, B])
            sage: g._solve_LP()
            Traceback (most recent call last):
            ...
            ValueError: Input game needs to be a two player constant sum game
        """
        if not self.is_constant_sum():
            raise ValueError("Input game needs to be a two player constant sum game")
        if solver == 'gambit':
            return self._solve_gambit_LP(maximization)

        sgn = 1
        if not maximization:
            sgn = -1

        strategy_sizes = [p.num_strategies for p in self.players]

        p = MixedIntegerLinearProgram(maximization=False, solver=solver)
        y = p.new_variable(nonnegative=True)
        v = p.new_variable(nonnegative=False)
        p.add_constraint(sgn * self.payoff_matrices()[0] * y - v[0] <= 0)
        p.add_constraint(matrix([[1] * strategy_sizes[1]]) * y == 1)
        p.set_objective(v[0])
        p.solve()
        y = tuple(p.get_values(y).values())

        p = MixedIntegerLinearProgram(maximization=False, solver=solver)
        x = p.new_variable(nonnegative=True)
        u = p.new_variable(nonnegative=False)
        p.add_constraint(sgn * -self.payoff_matrices()[0].T * x - u[0] <= 0)
        p.add_constraint(matrix([[1] * strategy_sizes[0]]) * x == 1)
        p.set_objective(u[0])
        p.solve()
        x = tuple(p.get_values(x).values())
        return [[x, y]]

    def _solve_enumeration(self, maximization=True):
        r"""
        Obtain the Nash equilibria using support enumeration.

        Algorithm implemented here is Algorithm 3.4 of [NN2007]_
        with an aspect of pruning from [SLB2008]_.

        1. For each k in 1...min(size of strategy sets)
        2. For each I,J supports of size k
        3. Prune: check if supports are dominated
        4. Solve indifference conditions and check that have Nash Equilibrium.

        EXAMPLES:

        A Game::

            sage: A = matrix([[160, 205, 44],
            ....:       [175, 180, 45],
            ....:       [201, 204, 50],
            ....:       [120, 207, 49]])
            sage: B = matrix([[2, 2, 2],
            ....:             [1, 0, 0],
            ....:             [3, 4, 1],
            ....:             [4, 1, 2]])
            sage: g=NormalFormGame([A, B])
            sage: g._solve_enumeration()
            [[(0, 0, 3/4, 1/4), (1/28, 27/28, 0)]]

        A game with 3 equilibria::

            sage: A = matrix([[3,3],
            ....:             [2,5],
            ....:             [0,6]])
            sage: B = matrix([[3,2],
            ....:             [2,6],
            ....:             [3,1]])
            sage: g = NormalFormGame([A, B])
            sage: g._solve_enumeration(maximization=False)
            [[(1, 0, 0), (0, 1)]]

        A simple example::

            sage: s = matrix([[1, 0], [-2, 3]])
            sage: t = matrix([[3, 2], [-1, 0]])
            sage: example = NormalFormGame([s, t])
            sage: example._solve_enumeration()
            [[(0, 1), (0, 1)], [(1/2, 1/2), (1/2, 1/2)], [(1, 0), (1, 0)]]

        Another::

            sage: A = matrix([[0, 1, 7, 1],
            ....:             [2, 1, 3, 1],
            ....:             [3, 1, 3, 5],
            ....:             [6, 4, 2, 7]])
            sage: B = matrix([[3, 2, 8, 4],
            ....:             [6, 2, 0, 3],
            ....:             [1, 3, -1, 1],
            ....:             [3, 2, 1, 1]])
            sage: C = NormalFormGame([A, B])
            sage: C._solve_enumeration()
            [[(0, 0, 0, 1), (1, 0, 0, 0)],
             [(2/7, 0, 0, 5/7), (5/11, 0, 6/11, 0)],
             [(1, 0, 0, 0), (0, 0, 1, 0)]]

        Again::

            sage: X = matrix([[1, 4, 2],
            ....:             [4, 0, 3],
            ....:             [2, 3, 5]])
            sage: Y = matrix([[3, 9, 2],
            ....:             [0, 3, 1],
            ....:             [5, 4, 6]])
            sage: Z = NormalFormGame([X, Y])
            sage: Z._solve_enumeration()
            [[(0, 0, 1), (0, 0, 1)], [(2/9, 0, 7/9), (0, 3/4, 1/4)], [(1, 0, 0), (0, 1, 0)]]

        TESTS:

        Due to the nature of the linear equations solved in this algorithm
        some negative vectors can be returned. Here is a test that ensures
        this doesn't happen (the particular payoff matrices chosen give a
        linear system that would have negative valued vectors as solution)::

            sage: a = matrix([[-13, 59],
            ....:             [27, 86]])
            sage: b = matrix([[14, 6],
            ....:             [58, -14]])
            sage: c = NormalFormGame([a, b])
            sage: c._solve_enumeration()
            [[(0, 1), (1, 0)]]

        Testing against an error in `_is_NE`.  Note that 1 equilibrium is
        missing: ``[(2/3, 1/3), (0, 1)]``, however this equilibrium has
        supports of different sizes. This only occurs in degenerate games
        and is not supported in the `enumeration` algorithm::

            sage: N = NormalFormGame([matrix(2,[0,-1,-2,-1]),matrix(2,[1,0,0,2])])
            sage: N._solve_enumeration()
            [[(0, 1), (0, 1)], [(1, 0), (1, 0)]]

        In this instance the `lrs` algorithm is able to find all
        three equilibria::

            sage: N = NormalFormGame([matrix(2,[0,-1,-2,-1]),matrix(2,[1,0,0,2])])
            sage: N.obtain_nash(algorithm='lrs')  # optional - lrslib
            [[(0, 1), (0, 1)], [(2/3, 1/3), (0, 1)], [(1, 0), (1, 0)]]

        Here is another::

            sage: N = NormalFormGame([matrix(2,[7,-8,-4,-8,7,0]),matrix(2,[-9,-1,-8,3,2,3])])
            sage: N._solve_enumeration()
            [[(0, 1), (0, 0, 1)]]
        """

        M1, M2 = self.payoff_matrices()
        if maximization is False:
            M1 = -M1
            M2 = -M2

        potential_supports = [[tuple(support) for support in
                               powerset(range(player.num_strategies))]
                              for player in self.players]

        potential_support_pairs = [pair for pair in product(*potential_supports) if len(pair[0]) == len(pair[1])]

        equilibria = []
        for pair in potential_support_pairs:
            # Check if any supports are dominated for row player
            if (self._row_cond_dominance(pair[0], pair[1], M1)
                # Check if any supports are dominated for col player
               and self._row_cond_dominance(pair[1], pair[0], M2.transpose())):
                a = self._solve_indifference(pair[0], pair[1], M2)
                b = self._solve_indifference(pair[1], pair[0], M1.transpose())
                if a and b and self._is_NE(a, b, pair[0], pair[1], M1, M2):
                    equilibria.append([tuple(a), tuple(b)])

        return sorted(equilibria)

    def _row_cond_dominance(self, p1_sup, p2_sup, matrix):
        r"""
        Check if any row strategies of a sub matrix defined
        by a given pair of supports are conditionally dominated.
        Return ``False`` if a row is conditionally dominated.

        TESTS:

        A matrix that depending on the support for the column player
        has a dominated row::

            sage: g = NormalFormGame()
            sage: A = matrix([[1, 1, 5], [2, 2, 0]])
            sage: g._row_cond_dominance((0, 1), (0, 1), A)
            False

        or does not have a dominated row::

            sage: g._row_cond_dominance((0, 1), (0, 2), A)
            True
        """
        subm = matrix.matrix_from_rows_and_columns(list(p1_sup), list(p2_sup))
        nbr_rows = subm.nrows()
        nbr_cols = subm.ncols()
        for s in range(nbr_rows):
            strategy = subm.rows()[s]
            for r in range(s, nbr_rows):
                row = subm.rows()[r]
                if strategy != row:
                    if all(strategy[i] < row[i] for i in range(nbr_cols)):
                        return False
                    if all(row[i] < strategy[i] for i in range(nbr_cols)):
                        return False
        return True

    def _solve_indifference(self, support1, support2, M):
        r"""
        For support1, returns the strategy with support: support2 that makes the
        column player indifferent for the utilities given by M.

        This is done by building the corresponding linear system.
        If  `\rho_1, \rho_2` are the supports of player 1 and 2 respectively.
        Then, indifference for player 1 implies:

        .. MATH::

            u_1(s_1,\rho_2) = u_1(s_2, \rho_2)

        for all `s_1, s_2` in the support of `\rho_1`. This corresponds to:

        .. MATH::

            \sum_{j\in S(\rho_2)}A_{s_1,j}{\rho_2}_j =
            \sum_{j\in S(\rho_2)}A_{s_2,j}{\rho_2}_j

        for all `s_1, s_2` in the support of `\rho_1` where `A` is the payoff
        matrix of player 1. Equivalently we can consider consecutive rows of
        `A` (instead of all pairs of strategies). Thus the corresponding
        linear system can be written as:

        .. MATH::

            \left(\sum_{j \in S(\rho_2)}^{A_{i,j} - A_{i+1,j}\right){\rho_2}_j

        for all `1\leq i \leq |S(\rho_1)|` (where `A` has been modified to only
        contain the row corresponding to `S(\rho_1)`). We also require all
        elements of `\rho_2` to sum to 1:

        .. MATH::

            \sum_{j\in S(\rho_1)}{\rho_2}_j = 1.

        TESTS:

        Find the indifference vector for a support pair that has
        no dominated strategies::

            sage: A = matrix([[1, 1, 5], [2, 2, 0]])
            sage: g = NormalFormGame([A])
            sage: g._solve_indifference((0, 1), (0, 2), A)
            (1/3, 2/3)
            sage: g._solve_indifference((0, 2), (0, 1), -A.transpose())
            (5/6, 0, 1/6)

        When a support pair has a dominated strategy there is no
        solution to the indifference equation::

            sage: g._solve_indifference((0, 1), (0, 1), -A.transpose())
            <BLANKLINE>

        Particular case of a game with 1 strategy for each for each player::

            sage: A = matrix([[10]])
            sage: g = NormalFormGame([A])
            sage: g._solve_indifference((0,), (0,), -A.transpose())
            (1)
        """
        linearsystem = matrix(QQ, len(support2) + 1, M.nrows())

        # Build linear system for player 1
        for strategy1 in support1:
            # Checking particular case of supports of pure strategies
            if len(support2) == 1:
                for strategy2 in range(M.ncols()):
                    if M[strategy1][support2[0]] < \
                            M[strategy1][strategy2]:
                        return False
            else:
                for strategy_pair2 in range(len(support2)):
                    # Coefficients of linear system that ensure indifference
                    # between two consecutive strategies of the support
                    linearsystem[strategy_pair2, strategy1] = \
                        M[strategy1][support2[strategy_pair2]] -\
                        M[strategy1][support2[strategy_pair2 - 1]]
            # Coefficients of linear system that ensure the vector is
            # a probability vector. ie. sum to 1
            linearsystem[-1, strategy1] = 1
        # Create rhs of linear systems
        linearsystem_rhs = vector([0 for i in range(len(support2))] + [1])

        # Solve both linear systems
        try:
            result = linearsystem.solve_right(linearsystem_rhs)
        except ValueError:
            return None

        return result

    def _is_NE(self, a, b, p1_support, p2_support, M1, M2):
        r"""
        For vectors that obey indifference for a given support pair,
        checks if it corresponds to a Nash equilibria (support is obeyed and
        no negative values, also that no player has incentive to deviate
        out of supports).

        TESTS::

            sage: X = matrix([[1, 4, 2],
            ....:             [4, 0, 3],
            ....:             [2, 3, 5]])
            sage: Y = matrix([[3, 9, 2],
            ....:             [0, 3, 1],
            ....:             [5, 4, 6]])
            sage: Z = NormalFormGame([X, Y])
            sage: Z._is_NE([0, 1/4, 3/4], [3/5, 2/5, 0], (1, 2,), (0, 1,), X, Y)
            False

            sage: Z._is_NE([2/9, 0, 7/9], [0, 3/4, 1/4], (0, 2), (1, 2), X, Y)
            True

        Checking pure strategies are not forgotten::

            sage: A = matrix(2, [0, -1, -2, -1])
            sage: B = matrix(2, [1, 0, 0, 2])
            sage: N = NormalFormGame([A, B])
            sage: N._is_NE([1, 0], [1, 0], (0,), (0,), A, B)
            True
            sage: N._is_NE([0, 1], [0, 1], (1,), (1,), A, B)
            True
            sage: N._is_NE([1, 0], [0, 1], (0,), (1,), A, B)
            False
            sage: N._is_NE([0, 1], [1, 0], (1,), (0,), A, B)
            False

            sage: A = matrix(3, [-7, -5,  5, 5,  5,  3,  1, -6,  1])
            sage: B = matrix(3, [-9, 7, 9, 6, -2, -3, -4, 6, -10])
            sage: N = NormalFormGame([A, B])
            sage: N._is_NE([1, 0, 0], [0, 0, 1], (0,), (2,), A, B)
            True
            sage: N._is_NE([0, 1, 0], [1, 0, 0], (1,), (0,), A, B)
            True
            sage: N._is_NE([0, 1, 0], [0, 1, 0], (1,), (1,), A, B)
            False
            sage: N._is_NE([0, 0, 1], [0, 1, 0], (2,), (1,), A, B)
            False
            sage: N._is_NE([0, 0, 1], [0, 0, 1], (2,), (2,), A, B)
            False
        """
        # Check that supports are obeyed
        if not (all(a[i] > 0 for i in p1_support) and
                all(b[j] > 0 for j in p2_support) and
                all(a[i] == 0 for i in range(len(a))
                    if i not in p1_support) and
                all(b[j] == 0 for j in range(len(b))
                    if j not in p2_support)):
            return False

        # Check that have pair of best responses

        p1_payoffs = [sum(v * row[i] for i, v in enumerate(b))
                      for row in M1.rows()]
        p2_payoffs = [sum(v * col[j] for j, v in enumerate(a))
                      for col in M2.columns()]

        #if p1_payoffs.index(max(p1_payoffs)) not in p1_support:
        if not any(i in p1_support for i, x in enumerate(p1_payoffs)
                   if x == max(p1_payoffs)):
            return False
        if not any(i in p2_support for i, x in enumerate(p2_payoffs)
                   if x == max(p2_payoffs)):
            return False

        return True

    def _Hrepresentation(self, m1, m2):
        r"""
        Create the H-representation strings required to use ``lrsnash``.

        Since lrslib 6.1, this format is referred to as "legacy format".

        This method is deprecated.

        EXAMPLES::

            sage: A = matrix([[1, 2], [3, 4]])
            sage: B = matrix([[3, 3], [1, 4]])
            sage: C = NormalFormGame([A, B])
            sage: print(C._Hrepresentation(A, B)[0])
            doctest:warning...
            DeprecationWarning: NormalFormGame._Hrepresentation is deprecated as it
            creates the legacy input format. Use NormalFormGame._lrs_nash_format instead
            See https://github.com/sagemath/sage/issues/27745 for details.
            H-representation
            linearity 1 5
            begin
            5 4 rational
            0 1 0 0
            0 0 1 0
            0 -3 -1 1
            0 -3 -4 1
            -1 1 1 0
            end
            <BLANKLINE>
            sage: print(C._Hrepresentation(A, B)[1])
            H-representation
            linearity 1 5
            begin
            5 4 rational
            0 -1 -2 1
            0 -3 -4 1
            0 1 0 0
            0 0 1 0
            -1 1 1 0
            end
            <BLANKLINE>
        """
        from sage.misc.superseded import deprecation
        deprecation(27745,
                    "NormalFormGame._Hrepresentation is deprecated as it "
                    "creates the legacy input format. "
                    "Use NormalFormGame._lrs_nash_format instead")

        from sage.geometry.polyhedron.misc import _to_space_separated_string
        m = self.players[0].num_strategies
        n = self.players[1].num_strategies
        midentity = list(matrix.identity(m))
        nidentity = list(matrix.identity(n))

        s = 'H-representation\n'
        s += 'linearity 1 ' + str(m + n + 1) + '\n'
        s += 'begin\n'
        s += str(m + n + 1) + ' ' + str(m + 2) + ' rational\n'
        for f in list(midentity):
            s += '0 ' + _to_space_separated_string(f) + ' 0 \n'
        for e in list(m2.transpose()):
            s += '0 ' + _to_space_separated_string(-e) + '  1 \n'
        s += '-1 '
        for g in range(m):
            s += '1 '
        s += '0 \n'
        s += 'end\n'

        t = 'H-representation\n'
        t += 'linearity 1 ' + str(m + n + 1) + '\n'
        t += 'begin\n'
        t += str(m + n + 1) + ' ' + str(n + 2) + ' rational\n'
        for e in list(m1):
            t += '0 ' + _to_space_separated_string(-e) + '  1 \n'
        for f in list(nidentity):
            t += '0 ' + _to_space_separated_string(f) + ' 0 \n'
        t += '-1 '
        for g in range(n):
            t += '1 '
        t += '0 \n'
        t += 'end\n'
        return s, t

    def _lrs_nash_format(self, m1, m2):
        r"""
        Create the input format for ``lrsnash``, version 6.1 or newer.

        EXAMPLES:

        An example from the ``lrsnash`` manual in the old and the format::

            sage: A = matrix([[0, 6], [2, 5], [3, 3]])
            sage: B = matrix([[1, 0], [0, 2], [4, 3]])
            sage: C = NormalFormGame([A, B])
            sage: print(C._lrs_nash_format(A, B))
            3 2
            <BLANKLINE>
            0 6
            2 5
            3 3
            <BLANKLINE>
            1 0
            0 2
            4 3
            <BLANKLINE>

            sage: legacy_format = C._Hrepresentation(A, B)
            doctest:warning...
            DeprecationWarning: NormalFormGame._Hrepresentation is deprecated as it
            creates the legacy input format. Use NormalFormGame._lrs_nash_format instead
            See https://github.com/sagemath/sage/issues/27745 for details.
            sage: print('*game: player 1\n', legacy_format[0])
            *game: player 1
            H-representation
            linearity 1 6
            begin
            6 5 rational
            0 1 0 0 0
            0 0 1 0 0
            0 0 0 1 0
            0 -1 0 -4  1
            0 0 -2 -3  1
            -1 1 1 1 0
            end
            <BLANKLINE>
            sage: print('*game: player 2\n', legacy_format[1])
            *game: player 2
            H-representation
            linearity 1 6
            begin
            6 4 rational
            0 0 -6  1
            0 -2 -5  1
            0 -3 -3  1
            0 1 0 0
            0 0 1 0
            -1 1 1 0
            end
        """
        from sage.geometry.polyhedron.misc import _to_space_separated_string
        m = self.players[0].num_strategies
        n = self.players[1].num_strategies
        s = f'{m} {n}\n\n'
        for r in m1.rows():
            s += _to_space_separated_string(r) + '\n'
        s += '\n'
        for r in m2.rows():
            s += _to_space_separated_string(r) + '\n'
        return s

    def is_degenerate(self, certificate=False):
        """
        A function to check whether the game is degenerate or not.
        Will return a boolean.

        A two-player game is called nondegenerate if no mixed strategy of
        support size `k` has more than `k` pure best responses [NN2007]_. In a
        degenerate game, this definition is violated, for example if there
        is a pure strategy that has two pure best responses.

        The implementation here transforms the search over mixed strategies to a
        search over supports which is a discrete search. A full explanation of
        this is given in [CK2015]_. This problem is known to be NP-Hard
        [Du2009]_.  Another possible implementation is via best response
        polytopes, see :issue:`18958`.

        The game Rock-Paper-Scissors is an example of a non-degenerate game,::

            sage: g = game_theory.normal_form_games.RPS()
            sage: g.is_degenerate()
            False

        whereas `Rock-Paper-Scissors-Lizard-Spock
        <http://www.samkass.com/theories/RPSSL.html>`_ is degenerate because
        for every pure strategy there are two best responses.::

            sage: g = game_theory.normal_form_games.RPSLS()
            sage: g.is_degenerate()
            True

        EXAMPLES:

        Here is an example of a degenerate game given in [DGRB2010]_::

            sage: A = matrix([[3, 3], [2, 5], [0, 6]])
            sage: B = matrix([[3, 3], [2, 6], [3, 1]])
            sage: degenerate_game = NormalFormGame([A,B])
            sage: degenerate_game.is_degenerate()
            True

        Here is an example of a degenerate game given in [NN2007]_::

            sage: A = matrix([[0, 6], [2, 5], [3, 3]])
            sage: B = matrix([[1, 0], [0, 2], [4, 4]])
            sage: d_game = NormalFormGame([A, B])
            sage: d_game.is_degenerate()
            True

        Here are some other examples of degenerate games::

            sage: M = matrix([[2, 1], [1, 1]])
            sage: N = matrix([[1, 1], [1, 2]])
            sage: game  = NormalFormGame([M, N])
            sage: game.is_degenerate()
            True

        If more information is required, it may be useful to use
        ``certificate=True``. This will return a boolean of whether the game is
        degenerate or not, and if True; a tuple containing the strategy where
        degeneracy was found and the player it belongs to. ``0`` is the row
        player and ``1`` is the column player.::

            sage: M = matrix([[2, 1], [1, 1]])
            sage: N = matrix([[1, 1], [1, 2]])
            sage: g  = NormalFormGame([M, N])
            sage: test, certificate = g.is_degenerate(certificate=True)
            sage: test, certificate
            (True, ((1, 0), 0))

        Using the output, we see that the opponent has more best responses than
        the size of the support of the strategy in question ``(1, 0)``. (We
        specify the player as ``(player + 1) % 2`` to ensure that we have the
        opponent's index.)::

            sage: g.best_responses(certificate[0], (certificate[1] + 1) % 2)
            [0, 1]

        Another example with a mixed strategy causing degeneracy.::

            sage: A = matrix([[3, 0], [0, 3], [1.5, 1.5]])
            sage: B = matrix([[4, 3], [2, 6], [3, 1]])
            sage: g = NormalFormGame([A, B])
            sage: test, certificate = g.is_degenerate(certificate=True)
            sage: test, certificate
            (True, ((1/2, 1/2), 1))

        Again, we see that the opponent has more best responses than the size of
        the support of the strategy in question ``(1/2, 1/2)``.::

            sage: g.best_responses(certificate[0], (certificate[1] + 1) % 2)
            [0, 1, 2]

        Sometimes, the different algorithms for obtaining nash_equilibria don't
        agree with each other. This can happen when games are degenerate::

            sage: a = matrix([[-75, 18, 45, 33],
            ....:            [42, -8, -77, -18],
            ....:            [83, 18, 11, 40],
            ....:            [-10, -38, 76, -9]])
            sage: b = matrix([[62, 64, 87, 51],
            ....:            [-41, -27, -69, 52],
            ....:            [-17, 25, -97, -82],
            ....:            [30, 31, -1, 50]])
            sage: d_game = NormalFormGame([a, b])
            sage: d_game.obtain_nash(algorithm='lrs')  # optional - lrslib
            [[(0, 0, 1, 0), (0, 1, 0, 0)],
             [(17/29, 0, 0, 12/29), (0, 0, 42/73, 31/73)],
             [(122/145, 0, 23/145, 0), (0, 1, 0, 0)]]
            sage: d_game.obtain_nash(algorithm='LCP')  # optional - gambit
            [[(0.5862068966, 0.0, 0.0, 0.4137931034),
              (0.0, 0.0, 0.5753424658, 0.4246575342)]]
            sage: d_game.obtain_nash(algorithm='enumeration')
            [[(0, 0, 1, 0), (0, 1, 0, 0)], [(17/29, 0, 0, 12/29), (0, 0, 42/73, 31/73)]]
            sage: d_game.is_degenerate()
            True

        TESTS::

            sage: g = NormalFormGame()
            sage: g.add_player(3)  # Adding first player with 3 strategies
            sage: g.add_player(3)  # Adding second player with 3 strategies
            sage: for key in g:
            ....:     g[key] = [0, 0]
            sage: g.is_degenerate()
            True

            sage: A = matrix([[3, 0], [0, 3], [1.5, 1.5]])
            sage: B = matrix([[4, 3], [2, 6], [3, 1]])
            sage: g = NormalFormGame([A, B])
            sage: g.is_degenerate()
            True

            sage: A = matrix([[1, -1], [-1, 1]])
            sage: B = matrix([[-1, 1], [1, -1]])
            sage: matching_pennies = NormalFormGame([A, B])
            sage: matching_pennies.is_degenerate()
            False

            sage: A = matrix([[2, 5], [0, 4]])
            sage: B = matrix([[2, 0], [5, 4]])
            sage: prisoners_dilemma = NormalFormGame([A, B])
            sage: prisoners_dilemma.is_degenerate()
            False

            sage: g = NormalFormGame()
            sage: g.add_player(2)
            sage: g.add_player(2)
            sage: g.add_player(2)
            sage: g.is_degenerate()
            Traceback (most recent call last):
            ...
            NotImplementedError: Tests for Degeneracy is not yet implemented for
             games with more than two players.
        """
        if len(self.players) > 2:
            raise NotImplementedError("Tests for Degeneracy is not yet "
                                      "implemented for games with more than "
                                      "two players.")

        d = self._is_degenerate_pure(certificate)
        if d:
            return d

        M1, M2 = self.payoff_matrices()
        potential_supports = [[tuple(support) for support in
                               powerset(range(player.num_strategies))]
                              for player in self.players]

        # filter out all supports that are pure or empty
        potential_supports = [[i for i in k if len(i) > 1]
                              for k in potential_supports]

        potential_support_pairs = [pair for pair in
                                   product(*potential_supports) if
                                   len(pair[0]) != len(pair[1])]

        # Sort so that solve small linear systems first
        potential_support_pairs.sort(key=lambda x: sum([len(k) for k in x]))

        for pair in potential_support_pairs:
            if len(pair[0]) < len(pair[1]):
                strat = self._solve_indifference(pair[0], pair[1], M2)
                if strat and len(self.best_responses(strat, player=0)) > len(pair[0]):
                    if certificate:
                        return True, (strat, 0)
                    else:
                        return True
            elif len(pair[1]) < len(pair[0]):
                strat = self._solve_indifference(pair[1], pair[0], M1.transpose())
                if strat and len(self.best_responses(strat, player=0)) > len(pair[1]):
                    if certificate:
                        return True, (strat, 1)
                    else:
                        return True

        if certificate:
            return False, ()
        else:
            return False

    def best_responses(self, strategy, player):
        """
        For a given strategy for a player and the index of the opponent,
        computes the payoff for the opponent and returns a list of the indices
        of the best responses. Only implemented for two player games

        INPUT:

        - ``strategy`` -- a probability distribution vector

        - ``player`` -- the index of the opponent, ``0`` for the row player,
          ``1`` for the column player

        EXAMPLES::

            sage: A = matrix([[3, 0], [0, 3], [1.5, 1.5]])
            sage: B = matrix([[4, 3], [2, 6], [3, 1]])
            sage: g = NormalFormGame([A, B])

        Now we can obtain the best responses for Player 1, when Player 2 uses
        different strategies::

            sage: g.best_responses((1/2, 1/2), player=0)
            [0, 1, 2]
            sage: g.best_responses((3/4, 1/4), player=0)
            [0]

        To get the best responses for Player 2 we pass the argument :code:`player=1`::

            sage: g.best_responses((4/5, 1/5, 0), player=1)
            [0, 1]

            sage: A = matrix([[1, 0], [0, 1], [0, 0]])
            sage: B = matrix([[1, 0], [0, 1], [0.7, 0.8]])
            sage: g = NormalFormGame([A, B])
            sage: g.best_responses((0, 1, 0), player=1)
            [1]

            sage: A = matrix([[3,3],[2,5],[0,6]])
            sage: B = matrix([[3,3],[2,6],[3,1]])
            sage: degenerate_game = NormalFormGame([A,B])
            sage: degenerate_game.best_responses((1, 0, 0), player=1)
            [0, 1]

            sage: A = matrix([[3, 0], [0, 3], [1.5, 1.5]])
            sage: B = matrix([[4, 3], [2, 6], [3, 1]])
            sage: g = NormalFormGame([A, B])
            sage: g.best_responses((1/3, 1/3, 1/3), player=1)
            [1]

        Note that this has only been implemented for 2 player games::

            sage: g = NormalFormGame()
            sage: g.add_player(2)  # adding first player with 2 strategies
            sage: g.add_player(2)  # adding second player with 2 strategies
            sage: g.add_player(2)  # adding third player with 2 strategies
            sage: g.best_responses((1/2, 1/2), player=2)
            Traceback (most recent call last):
            ...
            ValueError: Only available for 2 player games

        If the strategy is not of the correct dimension for the given player
        then an error is returned::

            sage: A = matrix([[3, 0], [0, 3], [1.5, 1.5]])
            sage: B = matrix([[4, 3], [2, 6], [3, 1]])
            sage: g = NormalFormGame([A, B])
            sage: g.best_responses((1/2, 1/2), player=1)
            Traceback (most recent call last):
            ...
            ValueError: Strategy is not of correct dimension

            sage: g.best_responses((1/3, 1/3, 1/3), player=0)
            Traceback (most recent call last):
            ...
            ValueError: Strategy is not of correct dimension

        If the strategy is not a true probability vector then an error is
        passed::

            sage: A = matrix([[3, 0], [0, 3], [1.5, 1.5]])
            sage: B = matrix([[4, 3], [2, 6], [3, 1]])
            sage: g = NormalFormGame([A, B])
            sage: g.best_responses((1/3, 1/2, 0), player=1)
            Traceback (most recent call last):
            ...
            ValueError: Strategy is not a probability distribution vector

            sage: A = matrix([[3, 0], [0, 3], [1.5, 1.5]])
            sage: B = matrix([[4, 3], [2, 6], [3, 1]])
            sage: g = NormalFormGame([A, B])
            sage: g.best_responses((3/2, -1/2), player=0)
            Traceback (most recent call last):
            ...
            ValueError: Strategy is not a probability distribution vector

        If the player specified is not `0` or `1`, an error is raised::

            sage: A = matrix([[3, 0], [0, 3], [1.5, 1.5]])
            sage: B = matrix([[4, 3], [2, 6], [3, 1]])
            sage: g = NormalFormGame([A, B])
            sage: g.best_responses((1/2, 1/2), player='Player1')
            Traceback (most recent call last):
            ...
            ValueError: Player1 is not an index of the opponent, must be 0 or 1
        """
        if len(self.players) != 2:
            raise ValueError('Only available for 2 player games')

        if player != 0 and player != 1:
            raise ValueError('%s is not an index of the opponent, must be 0 or 1' % player)

        strategy = vector(strategy)

        if sum(strategy) != 1 or min(strategy) < 0:
            raise ValueError('Strategy is not a probability distribution vector')

        if player == 0:
            payoff_matrix = self.payoff_matrices()[0]
        elif player == 1:
            payoff_matrix = self.payoff_matrices()[1].transpose()

        if len(strategy) != payoff_matrix.dimensions()[1]:
            raise ValueError('Strategy is not of correct dimension')

        payoffs = list(payoff_matrix * strategy)
        indices = [i for i, j in enumerate(payoffs) if j == max(payoffs)]

        return indices

    def _is_degenerate_pure(self, certificate=False):
        """
        Check whether a game is degenerate in pure strategies.

        TESTS::

            sage: A = matrix([[3,3],[2,5],[0,6]])
            sage: B = matrix([[3,3],[2,6],[3,1]])
            sage: degenerate_game = NormalFormGame([A,B])
            sage: degenerate_game._is_degenerate_pure()
            True

            sage: A = matrix([[1, 0], [0, 1], [0, 0]])
            sage: B = matrix([[1, 0], [0, 1], [0.7, 0.8]])
            sage: g = NormalFormGame([A, B])
            sage: g._is_degenerate_pure()
            False

            sage: A = matrix([[2, 5], [0, 4]])
            sage: B = matrix([[2, 0], [5, 4]])
            sage: prisoners_dilemma = NormalFormGame([A, B])
            sage: prisoners_dilemma._is_degenerate_pure()
            False

            sage: A = matrix([[0, -1, 1, 1, -1],
            ....:             [1, 0, -1, -1, 1],
            ....:             [-1, 1, 0, 1 , -1],
            ....:             [-1, 1, -1, 0, 1],
            ....:             [1, -1, 1, -1, 0]])
            sage: g = NormalFormGame([A])
            sage: g._is_degenerate_pure()
            True

        Whilst this game is not degenerate in pure strategies, it is
        actually degenerate, but only in mixed strategies::

            sage: A = matrix([[3, 0], [0, 3], [1.5, 1.5]])
            sage: B = matrix([[4, 3], [2, 6], [3, 1]])
            sage: g = NormalFormGame([A, B])
            sage: g._is_degenerate_pure()
            False
        """
        M1, M2 = self.payoff_matrices()
        for i, row in enumerate(M2.rows()):
            if list(row).count(max(row)) > 1:
                if certificate:
                    strat = [0 for k in range(M1.nrows())]
                    strat[i] = 1
                    return True, (tuple(strat), 0)
                else:
                    return True

        for j, col in enumerate(M1.columns()):
            if list(col).count(max(col)) > 1:
                if certificate:
                    strat = [0 for k in range(M1.ncols())]
                    strat[j] = 1
                    return True, (tuple(strat), 1)
                else:
                    return True
        return False


class _Player():
    def __init__(self, num_strategies):
        r"""
        TESTS::

            sage: from sage.game_theory.normal_form_game import _Player
            sage: p = _Player(5)
            sage: p.num_strategies
            5
        """
        self.num_strategies = num_strategies

    def add_strategy(self):
        r"""
        TESTS::

            sage: from sage.game_theory.normal_form_game import _Player
            sage: p = _Player(5)
            sage: p.add_strategy()
            sage: p.num_strategies
            6
        """
        self.num_strategies += 1
