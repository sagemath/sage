r"""
Normal form games with N players.

This module implements a class for normal form games (strategic form games)
[NN2007]_. At present the following algorithms are implemented to
compute equilibria of these games:

 * ``'enumeration'`` - An implementation of the support enumeration
   algorithm built in Sage.

 * An interface with 'gambit', which implements all of its available
   solvers (``'lcp'``, ``'enummixed'``, ``'lp'``, ``'gnm'``,
   ``'enumpure'``, ``'enumpoly'``, ``'liap'``, ``'simpdiv'``, ``'ipa'``
   and ``'logit'``). The ``'lcp'``, ``'enummixed'`` and ``'lp'`` solvers
   are restricted to 2 player games, while the others are able to solve
   games with an arbitrary number of players. See the gambit
   documentation
   (https://gambitproject.readthedocs.io/en/stable/pygambit.api.html).

 * ``'lp'`` - A built-in Sage implementation (with a gambit alternative)
   of a zero-sum game solver using linear programming. See
   :class:`MixedIntegerLinearProgram` for more on MILP solvers in Sage.

 * ``'lrs'`` - A solver interfacing with the 'lrslib' library.

The architecture for the class is based on the gambit architecture to
ensure an easy transition between gambit and Sage.  Most of the
algorithms for the computation of equilibria only solve 2 player games;
the ``'gnm'`` algorithm (which requires gambit) is able to solve games
with an arbitrary number of players.

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

The game can be drawn the way the strategic form of a game is usually
written down, as a table whose rows are Amy's choices and whose columns are
Bob's, with the utilities of a pair of choices written in its cell.  Each
player has a color, Amy red and Bob blue, and it is used for their name,
for their choices and for their utility in every cell::

    sage: battle_of_the_sexes.plot(player_labels=['Amy', 'Bob'],                        # needs sage.plot
    ....:                          strategy_labels=[['video games', 'movie'],
    ....:                                           ['video games', 'movie']],
    ....:                          best_responses=True, pure_nash=True)
    Graphics object consisting of 28 graphics primitives

.. PLOT::

    A = matrix([[3, 1], [0, 2]])
    B = matrix([[2, 1], [0, 3]])
    battle_of_the_sexes = NormalFormGame([A, B])
    choices = ['video games', 'movie']
    payoffs = battle_of_the_sexes.plot(player_labels=['Amy', 'Bob'],
                                       strategy_labels=[choices, choices],
                                       best_responses=True, pure_nash=True)
    sphinx_plot(payoffs)

An underlined utility is one a player gets by doing the best they can
against what the other player is doing, so a shaded cell, in which both
utilities are underlined, is one that neither of them wants to move away
from: exactly the two situations described above.  See
:meth:`~sage.game_theory.normal_form_game.NormalFormGame.plot` for games
with more players, which are drawn as pivot tables.

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
<https://www.samkass.com/theories/RPSSL.html>`_::

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
  :class:`MixedIntegerLinearProgram <sage.numerical.mip.MixedIntegerLinearProgram>` using the
  solver which was passed in with ``solver`` to solve the linear
  programming representation of the game. See
  :class:`~sage.numerical.mip.MixedIntegerLinearProgram` for more on MILP solvers in Sage.

* ``'lrs'``: Reverse search vertex enumeration for 2 player games. This
  algorithm uses the optional 'lrslib' package. To install it, type
  ``sage -i lrslib`` in the shell. For more information, see [Av2000]_.

* Sage also interfaces with all of the Nash equilibrium solvers provided by
  the open source game theory package
  `Gambit <https://www.gambit-project.org/>`_ [Gambit]_, including for games
  with more than 2 players. See :meth:`obtain_nash` for the algorithm names
  accepted here, and the `Gambit API overview
  <https://gambitproject.readthedocs.io/en/stable/pygambit.html>`_ for the
  underlying solvers.

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
    sage: matching_pennies.obtain_nash(algorithm='lcp')  # optional - pygambit
    [[(1/2, 1/2), (1/2, 1/2)]]
    sage: matching_pennies.obtain_nash(algorithm='lp', solver='PPL')
    [[(1/2, 1/2), (1/2, 1/2)]]
    sage: matching_pennies.obtain_nash(algorithm='lp', solver='gambit') # optional - pygambit
    [[(1/2, 1/2), (1/2, 1/2)]]
    sage: matching_pennies.obtain_nash(algorithm='enumeration')
    [[(1/2, 1/2), (1/2, 1/2)]]

Every one of those is exact: the gambit solvers that can compute with rational
arithmetic do so by default, so they agree with Sage's own algorithms on the
nose rather than to within a tolerance.  Pass ``rational=False`` to get the
floating point answer instead::

    sage: matching_pennies.obtain_nash(algorithm='lcp', rational=False)  # abs tol 1e-9 # optional - pygambit
    [[(0.5, 0.5), (0.5, 0.5)]]

The solvers that have no exact mode at all -- ``'gnm'``, ``'ipa'``, ``'liap'``,
``'logit'`` and ``'enumpoly'`` -- are exact here too, because their floating
point answer is rounded back to the equilibrium it approximates and only
accepted once gambit has confirmed, in exact arithmetic, that the rounded
profile really is one::

    sage: matching_pennies.obtain_nash(algorithm='liap')  # optional - pygambit
    [[(1/2, 1/2), (1/2, 1/2)]]

An equilibrium can be irrational once a game has three or more players, and
then there is nothing exact to round to; such an equilibrium is returned in
floating point, as is any equilibrium of a game whose own payoffs are inexact.

Note that if no algorithm argument is passed then the default will be
selected according to the following order (if the corresponding package is
installed):

1. ``'enumpoly'`` (if the game has more than 2 players; requires 'gambit')
2. ``'lp'`` (if the game is constant-sum; uses the solver chosen by Sage)
3. ``'lrs'`` (requires 'lrslib')
4. ``'enumeration'``

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
    ValueError: utilities have not been populated; ...
    sage: f.payoff_matrices()
    Traceback (most recent call last):
    ...
    ValueError: utilities have not been populated; ...

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

Just as a two player game can be created from two payoff matrices, an `N`
player game can be created directly from a list of `N` payoff arrays. Since a
Sage matrix is only two-dimensional, the payoffs of a game with more than two
players are given as `N`-dimensional numpy arrays of the same shape (one per
player)::

    sage: import numpy as np
    sage: A = np.array([[[3, 1], [2, 3]], [[9, 3], [8, 2]]])
    sage: B = np.array([[[1, 5], [6, 5]], [[7, 2], [4, 6]]])
    sage: C = np.array([[[4, 9], [5, 8]], [[9, 3], [6, 4]]])
    sage: NormalFormGame([A, B, C]) == threegame
    True

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

Games with more than 2 players can be solved using the ``'gnm'``
algorithm, which interfaces with gambit's implementation of the global
Newton method::

    sage: threegame.obtain_nash(algorithm='gnm')  # optional - pygambit
    [[(0, 1), (0, 1), (0, 1)]]

When no algorithm is given for a game with more than 2 players, the
``'enumpoly'`` algorithm is selected by default; it interfaces with
gambit's enumeration of equilibria via systems of polynomial equations
and, unlike ``'gnm'``, returns all of the equilibria it finds::

    sage: threegame.obtain_nash()  # optional - pygambit
    [[(0, 1), (0, 1), (0, 1)], [(1, 0), (0, 1), (0, 1)]]

Note that ``'gnm'`` is a numerical algorithm and so returns only a sample of
the equilibria.

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
    sage: g.obtain_nash(algorithm='lcp')  # optional - pygambit
    [[(1/5, 4/5), (3/5, 2/5)]]

It is also possible to generate a Normal form game from a gambit Game::

    sage: # optional - pygambit
    sage: import numpy as np
    sage: from pygambit import Game
    sage: gambitgame = Game.from_arrays(np.array([[8., 2.], [10., 5.]]),
    ....:                               np.array([[8., 10.], [2., 5.]]))
    sage: g = NormalFormGame(gambitgame); g
    Normal Form Game with the following utilities: {(0, 0): [8.0, 8.0],
     (0, 1): [2.0, 10.0],
     (1, 0): [10.0, 2.0],
     (1, 1): [5.0, 5.0]}

For more information on using Gambit in Sage see ``Using Gambit in
Sage``. This includes how to access Gambit
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
    sage: g.obtain_nash(algorithm='lcp')  # optional - pygambit
    [[(0, 0, 0, 0, 0, 0, 0, 0, 1), (0, 0, 0, 0, 0, 0, 0, 0, 1)]]

The output is a pair of vectors (as before) showing the Nash equilibrium.
In particular it here shows that out of the 10 possible strategies both
players should choose the last. Recall that the above considers a reduced
version of the game where individuals can claim integer values from 10
to 2.  The equilibrium strategy is thus for both players to state that
the value of their suitcase is 2.

Several standard Normal Form Games have also been implemented.
For more information on how to access these, see:
:mod:`Game Theory Catalog<sage.game_theory.catalog_normal_form_games>`.
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
<https://www.samkass.com/theories/RPSSL.html>`_ any pure strategy has two best
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
    sage: degenerate_game.obtain_nash(algorithm='lcp')  # optional - pygambit
    [[(0, 1/3, 2/3), (1/3, 2/3)], [(1, 0, 0), (2/3, 1/3)], [(1, 0, 0), (1, 0)]]
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

import os
from collections.abc import MutableMapping
from decimal import Decimal
from itertools import product
from .parser import Parser
from sage.misc.latex import latex
from sage.combinat.subset import powerset
from sage.rings.rational_field import QQ
from sage.rings.real_mpfr import RR
from sage.structure.sage_object import SageObject
from sage.matrix.constructor import matrix
from sage.matrix.constructor import vector
from sage.misc.temporary_file import tmp_filename, atomic_write
from sage.misc.sage_ostools import redirection
from sage.numerical.mip import MixedIntegerLinearProgram
from sage.cpython.string import bytes_to_str

try:
    import numpy as np
except ImportError:
    np = None

from sage.features.gambit import pygambit
from sage.misc.lazy_import import lazy_import
lazy_import('pygambit', ['Game', 'Rational', 'read_nfg', 'catalog'], feature=pygambit())
lazy_import('pygambit', ['MixedBehaviorProfileDouble', 'MixedStrategyProfileDouble'],
            feature=pygambit())
lazy_import('pygambit', 'nash', 'gambit_nash', feature=pygambit())

# The sizes :meth:`NormalFormGame.plot` draws a payoff table with, in the data
# coordinates of the picture, where one row of the table is one unit tall: the
# width one character of a label takes at the reference font size, the blank
# space left around the widest label of a column, how many inches of figure a
# cell is given, and the size in inches past which the figure does not grow any
# further.
_CHAR_WIDTH = 0.183
_CELL_PADDING = 0.7
_INCHES_PER_CELL = 0.5
_MAX_FIGSIZE = 30

# The characters of a label that :func:`_label_width` counts as half a
# character wide, the rest of them being about as wide as a digit.
_NARROW_CHARACTERS = " ,.:;'!|"

# How heavily :meth:`NormalFormGame.plot` rules a payoff table: the weight of
# its border, and the weight the rules between the strategies of a player lose
# with every player nested inside them.
_OUTER_RULE = 2
_RULE_STEP = 0.5

# The color :meth:`NormalFormGame.plot` shades a pure Nash equilibrium with.
# The players own the colors of ``_PLAYER_COLORS``, so it is none of those.
_NASH_COLOR = 'lightyellow'


def _label_width(label):
    r"""
    Return the width :meth:`NormalFormGame.plot` draws ``label`` with, in the
    data coordinates of its picture.

    The font a plot is drawn with is proportional, so this is an estimate: a
    digit, a letter and the like are taken to be ``_CHAR_WIDTH`` wide and the
    punctuation of ``_NARROW_CHARACTERS`` half of that, which is close enough
    for laying out payoffs and centering them in their cell.

    EXAMPLES::

        sage: from sage.game_theory.normal_form_game import _label_width
        sage: _label_width('') == 0
        True
        sage: _label_width('12') == 2 * _label_width('1')
        True

    A comma and a space together are as wide as a single digit::

        sage: _label_width('3, 4') == 3 * _label_width('3')
        True
    """
    return _CHAR_WIDTH * sum(0.5 if character in _NARROW_CHARACTERS else 1
                             for character in label)


def _gambit_payoffs_are_exact(game):
    r"""
    Return whether every payoff of the gambit game ``game`` is exact.

    Gambit stores a payoff either as a rational or as a decimal, and both
    :meth:`NormalFormGame._gambit_` and
    :class:`~sage.game_theory.extensive_form_game.ExtensiveFormGame` arrange for
    an exact Sage payoff to arrive as the former and an inexact one as the
    latter.  So this asks of the gambit game what "the payoffs of this game are
    exact" means for a Sage game.  For a tree game the probabilities of the
    chance moves are payoff data too, and are checked as well.

    This requires the optional gambit package.

    EXAMPLES::

        sage: # optional - pygambit
        sage: from sage.game_theory.normal_form_game import _gambit_payoffs_are_exact
        sage: A = matrix([[1, 2], [3, 4]])
        sage: _gambit_payoffs_are_exact(NormalFormGame([A])._gambit_())
        True
        sage: _gambit_payoffs_are_exact(NormalFormGame([A / 3])._gambit_())
        True
        sage: _gambit_payoffs_are_exact(NormalFormGame([A.change_ring(RDF)])._gambit_())
        False

    A chance move with inexact probabilities makes a tree game inexact::

        sage: # optional - pygambit
        sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
        sage: g = ExtensiveFormGame(players=['Alice'])
        sage: g.append_chance_move(g.root, ['H', 'T'], probs=['1/3', '2/3'])
        sage: for branch in ['H', 'T']:
        ....:     g.append_move(g.root.children[branch], 'Alice', ['L', 'R'])
        ....:     for move in ['L', 'R']:
        ....:         g.set_outcome(g.root.children[branch].children[move],
        ....:                       branch + move, [1])
        sage: _gambit_payoffs_are_exact(g._gambit_())
        True
        sage: g.set_chance_probs(g.root, [0.25, 0.75])
        sage: _gambit_payoffs_are_exact(g._gambit_())
        False
    """
    if any(not isinstance(outcome[player], Rational)
           for outcome in game.outcomes for player in game.players):
        return False
    if not game.is_tree:
        return True
    return all(isinstance(action.prob, Rational)
               for infoset in game.players.chance.infosets
               for action in infoset.actions)


def _rationalize_gambit_profile(profile, tolerance, agent=False):
    r"""
    Return an exact version of the gambit mixed profile ``profile``, or ``None``.

    Several of gambit's solvers -- ``'gnm'``, ``'ipa'``, ``'logit'``, ``'liap'``
    and ``'enumpoly'`` -- have no exact mode and answer in floating point even
    when the game they are given is exact.  Their answer is nonetheless usually
    a floating point rendering of an exact equilibrium, so this rounds each
    probability to the simplest rational within ``tolerance`` of it and then
    asks gambit whether the result is an equilibrium *exactly*: a rational
    profile computes its regret over the rationals, so a maximum regret of zero
    is a proof rather than an estimate.

    ``None`` is returned when the rounded profile fails that test, which is what
    happens when the equilibrium is genuinely irrational -- possible as soon as
    a game has three players.

    A profile that is already exact is returned unchanged.

    INPUT:

    - ``profile`` -- a pygambit mixed strategy or mixed behavior profile

    - ``tolerance`` -- a positive number; how far a probability may be moved to
      reach a rational.  The larger it is the simpler the rationals that are
      tried, so a rounding that is too fine can verify an equilibrium in an
      unenlightening form rather than fail.

    - ``agent`` -- boolean (default: ``False``); whether the profile is an
      *agent* equilibrium of an extensive form game rather than a Nash
      equilibrium, in which case it is verified with ``agent_max_regret()``.
      An agent equilibrium need only be immune to deviations at one information
      set at a time, so testing it for Nash equilibrium -- as the gambit agent
      solvers' answers would otherwise be tested -- can reject a profile that is
      exactly what was asked for.  See
      :meth:`~sage.game_theory.extensive_form_game.ExtensiveFormGame.obtain_nash`.

    OUTPUT: an exact pygambit profile, or ``None``

    This requires the optional gambit package.

    EXAMPLES:

    A game whose only mixed equilibrium is `(3/5, 2/5), (2/5, 3/5)`; ``'liap'``
    finds it to within a few parts in `10^9`, which is plenty to recover it::

        sage: # optional - pygambit
        sage: from pygambit.nash import liap_solve
        sage: from sage.game_theory.normal_form_game import _rationalize_gambit_profile
        sage: battle = NormalFormGame([matrix([[3, 0], [0, 2]]),
        ....:                          matrix([[2, 0], [0, 3]])])
        sage: g = battle._gambit_()
        sage: eq = liap_solve(g.mixed_strategy_profile(rational=False)).equilibria[0]
        sage: _rationalize_gambit_profile(eq, 1e-6)
        [[Rational(3, 5), Rational(2, 5)], [Rational(2, 5), Rational(3, 5)]]

    Round too finely and nothing exact is found::

        sage: _rationalize_gambit_profile(eq, 1e-15) is None    # optional - pygambit
        True

    An exact profile is passed through as it is::

        sage: # optional - pygambit
        sage: exact = g.mixed_strategy_profile(rational=True)
        sage: _rationalize_gambit_profile(exact, 1e-6) is exact
        True

    With ``agent`` the profile is verified as an agent equilibrium instead.  In
    Myerson's game below Player 1 can gain 1 by changing what they do at *both*
    of their information sets, but nothing by changing either one alone, so the
    profile is an agent equilibrium and not a Nash equilibrium::

        sage: # optional - pygambit
        sage: from pygambit import catalog
        sage: myerson = catalog.load('books/myerson1991/fig4_2')
        sage: profile = myerson.mixed_behavior_profile(rational=False)
        sage: for label, prob in [('A1', 0), ('B1', 1), ('Y1', 0),
        ....:                     ('Z1', 1), ('W2', 1), ('X2', 0)]:
        ....:     profile[myerson.actions[label]] = float(prob)
        sage: float(profile.max_regret()), float(profile.agent_max_regret())
        (1.0, 0.0)
        sage: _rationalize_gambit_profile(profile, 1e-6) is None
        True
        sage: _rationalize_gambit_profile(profile, 1e-6, agent=True)
        [[[Rational(0, 1), Rational(1, 1)], [Rational(0, 1), Rational(1, 1)]],
         [[Rational(1, 1), Rational(0, 1)]]]
    """
    if not isinstance(profile, (MixedBehaviorProfileDouble,
                                MixedStrategyProfileDouble)):
        return profile

    game = profile.game
    if isinstance(profile, MixedBehaviorProfileDouble):
        exact = game.mixed_behavior_profile(rational=True)
        # Indexing by the actions of the players' own information sets leaves
        # the chance moves, which are not part of a behavior profile, alone.
        keys = [action for player in game.players
                for infoset in player.infosets for action in infoset.actions]
    else:
        exact = game.mixed_strategy_profile(rational=True)
        keys = [strategy for player in game.players
                for strategy in player.strategies]

    for key in keys:
        # A solver may undershoot zero by a hair; a probability may not.
        nearby = max(RR(profile[key]).nearby_rational(max_error=tolerance),
                     QQ.zero())
        exact[key] = Rational(int(nearby.numerator()), int(nearby.denominator()))

    try:
        # Rounding the probabilities of a player separately need not leave them
        # summing to one; gambit's own normalisation is exact.
        exact = exact.normalize()
    except ValueError:      # all of some player's probabilities rounded to zero
        return None

    # An agent equilibrium is only required to be immune to a deviation at a
    # single information set, so it is the agent regret that must vanish.
    regret = exact.agent_max_regret() if agent else exact.max_regret()
    return exact if regret == 0 else None


class NormalFormGame(SageObject, MutableMapping):
    r"""
    An object representing a Normal Form Game. Primarily used to compute the
    Nash Equilibria.

    INPUT:

    - ``generator`` -- can be any of the following:

      * a list of `N` payoff arrays describing an `N` player game. For a 2
        player game these are the two payoff matrices (Sage matrices); for a
        game with more than 2 players each payoff array is an `N`-dimensional
        numpy array (a Sage matrix being only 2-dimensional cannot describe
        the payoffs of a game with more than 2 players). All arrays must have
        the same shape and the number of arrays must equal their common number
        of dimensions.

      * a single matrix, in which case a 2 player zero-sum game is constructed.

      * a gambit ``Game``.

      * left blank, in which case the utilities are populated manually.
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
            sage: threegame.obtain_nash(algorithm='gnm')  # optional - pygambit
            [[(0, 1), (1, 0), (1, 0)]]

        Rather than populating the utilities by hand, the same game can be
        built directly from a list of payoff arrays, one ``N``-dimensional
        numpy array per player::

            sage: import numpy as np
            sage: A = np.array([[[3, 1], [2, 3]], [[9, 3], [8, 2]]])
            sage: B = np.array([[[1, 5], [6, 5]], [[7, 2], [4, 6]]])
            sage: C = np.array([[[4, 9], [5, 8]], [[9, 3], [6, 4]]])
            sage: arraygame = NormalFormGame([A, B, C])
            sage: arraygame == threegame
            True
            sage: arraygame
            Normal Form Game with the following utilities: {(0, 0, 0): [3, 1, 4],
             (0, 0, 1): [1, 5, 9],
             (0, 1, 0): [2, 6, 5],
             (0, 1, 1): [3, 5, 8],
             (1, 0, 0): [9, 7, 9],
             (1, 0, 1): [3, 2, 3],
             (1, 1, 0): [8, 4, 6],
             (1, 1, 1): [2, 6, 4]}

        Can initialise a game from a gambit game object::

            sage: # optional - pygambit
            sage: import numpy as np
            sage: from pygambit import Game
            sage: gambitgame = Game.from_arrays(np.array([[5., 2.], [10., 5.]]),
            ....:                               np.array([[8., 11.], [7., 5.]]))
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

        Raise an error if the number of payoff arrays does not match the number
        of players, i.e. the number of dimensions of each array (here three
        2-dimensional matrices cannot describe a 3 player game)::

            sage: p1 = matrix([[1, 2], [3, 4]])
            sage: error = NormalFormGame([p1, p1, p1])
            Traceback (most recent call last):
            ...
            ValueError: the number of matrices must match the number of players,
             i.e. the number of dimensions of each matrix

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
        # ``Game`` is a lazy import, so it must be compared with ``isinstance``
        # (which resolves it) rather than with ``is``, and only once pygambit
        # is known to be present -- otherwise resolving it would raise instead
        # of giving the ``TypeError`` below.
        if isinstance(generator, list):
            if len(generator) == 1:
                generator.append(-generator[-1])
            self._n_matrix_game(generator)
        elif pygambit().is_present() and isinstance(generator, Game):
            self._gambit_game(generator)
        elif generator is not None:
            raise TypeError("Generator function must be a list, gambit game or nothing")

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

    def _repr_(self) -> str:
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

    def _latex_(self) -> str:
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

    def plot(self, best_responses=False, pure_nash=False, player_labels=None,
             strategy_labels=None, **kwargs):
        r"""
        Plot the game as a payoff table.

        The table is drawn with Sage's own plotting and needs nothing beyond
        Sage.  The first player owns the rows of the table and the second its
        columns, so a two player game comes out as the usual bimatrix, with the
        payoffs of a strategy profile written in its cell.  Each player has a
        color -- the first red, the second blue, then green, orange, purple and
        brown, cycling for a game with more players than that -- which is used
        for their name, for their strategies and for their payoff in every
        cell, so that a column of the table can be read off at a glance.

        A game with more than two players is drawn as a pivot table: every
        player after the first is nested inside the columns, the second player
        outermost, so that each of their strategies splits every column of the
        table again.  The whole game is therefore visible at once, and the
        rules between the columns say how deep the split is: the boundary
        between two strategies of an outer player is drawn more heavily than
        one between the strategies of a player nested inside it.

        The figure grows with the table, up to a point, so that the cells keep
        their size however large the game is; past that size the labels shrink
        instead.  Both can be overridden with the ``figsize`` and ``fontsize``
        keywords.  The result is a :class:`~sage.plot.graphics.Graphics`
        object, so it composes with the rest of Sage's plotting.

        INPUT:

        - ``best_responses`` -- boolean (default: ``False``); whether to
          underline, in that player's color, every payoff a player receives by
          playing a best response to what the other players play.  All the
          strategies tying for the best payoff are underlined.  A cell whose
          payoffs are all underlined is a pure Nash equilibrium

        - ``pure_nash`` -- boolean (default: ``False``); whether to shade the
          cells that are pure Nash equilibria

        - ``player_labels`` -- list (default: ``None``); a name for each
          player, in the order the game lists them.  By default the players are
          named ``'Player 1'``, ``'Player 2'`` and so on

        - ``strategy_labels`` -- list (default: ``None``); a list of names for
          each player's strategies, in the order the game lists the players.
          By default a strategy is named by its index, so that the labels of a
          cell are the key it is indexed by, as in ``game[0, 1]``

        - ``**kwargs`` -- passed on to the resulting
          :class:`~sage.plot.graphics.Graphics`, so any option of
          :meth:`~sage.plot.graphics.Graphics.show` can be given here; besides
          ``figsize``, ``fontsize`` sets the size the labels are drawn at

        OUTPUT: a :class:`~sage.plot.graphics.Graphics` object

        EXAMPLES:

        The prisoner's dilemma, drawn as the bimatrix it is::

            sage: # needs sage.plot
            sage: from sage.plot.graphics import Graphics
            sage: A = matrix([[2, 5], [0, 4]])
            sage: B = matrix([[2, 0], [5, 4]])
            sage: prisoners_dilemma = NormalFormGame([A, B])
            sage: p = prisoners_dilemma.plot(); p
            Graphics object consisting of 22 graphics primitives
            sage: isinstance(p, Graphics)
            True

        Both players are named beside the table, the first down its side and
        the second above its columns::

            sage: # needs sage.plot
            sage: sorted(t.string for t in p
            ....:        if hasattr(t, 'string') and 'Player' in t.string)
            ['Player 1', 'Player 2']

        Everything written in the table belongs to one of the two players, so
        it is drawn in red or in blue::

            sage: # needs sage.plot
            sage: sorted(set(t.options()['rgbcolor'] for t in p
            ....:            if hasattr(t, 'string')))
            [(0.0, 0.0, 1.0), (1.0, 0.0, 0.0)]

        With ``best_responses`` every payoff a player gets by playing a best
        response is underlined in that player's color.  Here there are four
        such payoffs, two of them the pair in the single pure Nash
        equilibrium::

            sage: # needs sage.plot
            sage: from sage.plot.line import Line
            sage: underlined = [q for q in
            ....:               prisoners_dilemma.plot(best_responses=True)
            ....:               if isinstance(q, Line)
            ....:               and q.options()['rgbcolor'] != 'black']
            sage: len(underlined)
            4
            sage: sorted(set(q.options()['rgbcolor'] for q in underlined))
            ['blue', 'red']

        With ``pure_nash`` that equilibrium is shaded.  It is the cell in the
        top left, where both players play their first strategy::

            sage: # needs sage.plot
            sage: from sage.plot.polygon import Polygon
            sage: shaded = [q for q in prisoners_dilemma.plot(pure_nash=True)
            ....:           if isinstance(q, Polygon)]
            sage: len(shaded)
            1
            sage: shaded[0].ydata   # the top row of the table
            [0.0, 0.0, -1.0, -1.0]
            sage: shaded[0].xdata[0]   # its first column
            0.0

        Since the players and their strategies have no names of their own,
        they are named after their indices, but any names can be given::

            sage: # needs sage.plot
            sage: A = matrix([[3, 1], [0, 2]])
            sage: B = matrix([[2, 1], [0, 3]])
            sage: battle_of_the_sexes = NormalFormGame([A, B])
            sage: named = battle_of_the_sexes.plot(
            ....:     player_labels=['Amy', 'Bob'],
            ....:     strategy_labels=[['video games', 'movie'],
            ....:                      ['video games', 'movie']])
            sage: sorted(set(t.string for t in named
            ....:            if hasattr(t, 'string') and not t.string[0].isdigit()))
            ['Amy', 'Bob', 'movie', 'video games']

        A game with three players is drawn as a pivot table, the third player
        nested inside the columns of the second, and its payoffs are written in
        three colors::

            sage: # needs sage.plot
            sage: import numpy as np
            sage: A = np.array([[[3, 1], [2, 3]], [[9, 3], [8, 2]]])
            sage: B = np.array([[[1, 5], [6, 5]], [[7, 2], [4, 6]]])
            sage: C = np.array([[[4, 9], [5, 8]], [[9, 3], [6, 4]]])
            sage: three = NormalFormGame([A, B, C])
            sage: three.plot()
            Graphics object consisting of 46 graphics primitives
            sage: from sage.plot.colors import Color
            sage: written = set(t.options()['rgbcolor'] for t in three.plot()
            ....:               if hasattr(t, 'string'))
            sage: sorted(Color(c).html_color() for c in written)
            ['#0000ff', '#008000', '#ff0000']

        It has two pure Nash equilibria, so two of its eight cells are
        shaded::

            sage: # needs sage.plot
            sage: len([q for q in three.plot(pure_nash=True)
            ....:      if isinstance(q, Polygon)])
            2

        TESTS:

        The table is drawn to scale with the axes off::

            sage: # needs sage.plot
            sage: p.aspect_ratio()
            1.0
            sage: p.axes()
            False

        The figure grows with the game, and a ``figsize`` given by hand wins::

            sage: # needs sage.plot
            sage: (three.plot()._extra_kwds['figsize'][0]
            ....:  > prisoners_dilemma.plot()._extra_kwds['figsize'][0])
            True
            sage: prisoners_dilemma.plot(figsize=4)._extra_kwds['figsize']
            4

        A game whose utilities have not all been set cannot be drawn::

            sage: # needs sage.plot
            sage: g = NormalFormGame()
            sage: g.add_player(2)
            sage: g.add_player(2)
            sage: g.plot()
            Traceback (most recent call last):
            ...
            ValueError: utilities have not been populated; ...

        Nor can a game with no players at all::

            sage: NormalFormGame().plot()                                       # needs sage.plot
            Traceback (most recent call last):
            ...
            ValueError: a game with no players cannot be plotted

        The labels given by hand must match the game::

            sage: # needs sage.plot
            sage: prisoners_dilemma.plot(player_labels=['Amy'])
            Traceback (most recent call last):
            ...
            ValueError: there must be one player label per player
            sage: prisoners_dilemma.plot(strategy_labels=[['a', 'b', 'c'],
            ....:                                         ['a', 'b']])
            Traceback (most recent call last):
            ...
            ValueError: there must be one strategy label per strategy of every player

        A game with a single player is a column of payoffs::

            sage: # needs sage.plot
            sage: g = NormalFormGame()
            sage: g.add_player(3)
            sage: for s in range(3):
            ....:     g[(s,)] = [s * s]
            sage: g.plot()
            Graphics object consisting of 14 graphics primitives
        """
        from sage.plot.graphics import Graphics
        from sage.plot.line import line2d
        from sage.plot.polygon import polygon2d
        from sage.plot.text import text
        # The two game classes give their players the same colors.
        from sage.game_theory.extensive_form_game import _PLAYER_COLORS

        number_of_players = len(self.players)
        if not number_of_players:
            raise ValueError("a game with no players cannot be plotted")
        self._require_complete()

        shape = [p.num_strategies for p in self.players]
        if player_labels is None:
            player_labels = ["Player {0}".format(i + 1)
                             for i in range(number_of_players)]
        elif len(player_labels) != number_of_players:
            raise ValueError("there must be one player label per player")
        player_labels = [str(label) for label in player_labels]
        if strategy_labels is None:
            strategy_labels = [[str(s) for s in range(k)] for k in shape]
        elif [len(labels) for labels in strategy_labels] != shape:
            raise ValueError("there must be one strategy label per strategy "
                             "of every player")
        strategy_labels = [[str(label) for label in labels]
                           for labels in strategy_labels]

        color_of_player = [_PLAYER_COLORS[i % len(_PLAYER_COLORS)]
                           for i in range(number_of_players)]

        # The first player owns the rows and every other player a band of
        # column headers, the second player's outermost.  A two player game
        # therefore comes out as the usual bimatrix, and every further player
        # splits each of its columns again.
        rows = shape[0]
        column_profiles = list(product(*(range(k) for k in shape[1:])))
        columns = len(column_profiles)
        bands = number_of_players - 1

        # The payoffs of a profile are written out as one string per player, so
        # that each can be drawn in its player's color, the last one carrying no
        # comma.  A cell is made wide enough for the longest of them.
        pieces_of = {}
        for row in range(rows):
            for column_profile in column_profiles:
                profile = (row,) + column_profile
                payoffs = [str(payoff) for payoff in self.utilities[profile]]
                pieces_of[profile] = ([payoff + ", " for payoff in payoffs[:-1]]
                                      + payoffs[-1:])

        def width_of(labels):
            return _CELL_PADDING + max(_label_width(label) for label in labels)

        cell_width = width_of(["".join(pieces) for pieces in pieces_of.values()]
                              + [label for labels in strategy_labels[1:]
                                 for label in labels])
        # The corner of the table holds the name of every column player, beside
        # its band of headers, and the strategies of the row player below them.
        header_width = width_of(strategy_labels[0] + player_labels[1:])

        best = self._pure_best_responses() if best_responses or pure_nash else set()

        plot = Graphics()

        # Shade the cells that are pure Nash equilibria, underneath the rules.
        if pure_nash:
            for row in range(rows):
                for column, column_profile in enumerate(column_profiles):
                    profile = (row,) + column_profile
                    if any((profile, player) not in best
                           for player in range(number_of_players)):
                        continue
                    left, right = column * cell_width, (column + 1) * cell_width
                    plot += polygon2d([(left, -row), (right, -row),
                                       (right, -row - 1), (left, -row - 1)],
                                      color=_NASH_COLOR, zorder=1)

        # Rule the table.  Every horizontal rule runs the whole way across, but
        # a vertical one only reaches up as far as the band of the player whose
        # strategies it separates: above that it would cut a header in two.
        table_left, table_right = -header_width, columns * cell_width
        for y in range(-rows, bands + 1):
            weight = _OUTER_RULE if y in (-rows, 0, bands) else _RULE_STEP
            plot += line2d([(table_left, y), (table_right, y)],
                           color='black', thickness=weight, zorder=3)
        for column in range(columns + 1):
            if column in (0, columns):
                depth, weight = 0, _OUTER_RULE
            else:
                # The rule is as heavy as the outermost player it separates.
                depth = min(band for band in range(bands)
                            if column_profiles[column - 1][band]
                            != column_profiles[column][band])
                weight = max(_RULE_STEP, _OUTER_RULE - depth * _RULE_STEP)
            plot += line2d([(column * cell_width, bands - depth),
                            (column * cell_width, -rows)],
                           color='black', thickness=weight, zorder=3)
        plot += line2d([(table_left, bands), (table_left, -rows)],
                       color='black', thickness=_OUTER_RULE, zorder=3)

        # The figure grows with the table so that the cells keep their size
        # however large the game is; past the cap the labels shrink instead.
        width = _INCHES_PER_CELL * (header_width + columns * cell_width + 1)
        height = _INCHES_PER_CELL * (rows + bands + 1)
        shrink = min(1, _MAX_FIGSIZE / max(width, height))
        kwargs.setdefault('figsize', (width * shrink, height * shrink))
        fontsize = kwargs.pop('fontsize', max(5, 10 * shrink))

        def draw(string, xy, player, **options):
            return text(string, xy, color=color_of_player[player],
                        fontsize=fontsize, zorder=5, **options)

        # The row player is named down the side of the table, and every column
        # player in the corner beside its own band of headers.
        plot += draw(player_labels[0], (table_left - 0.45, -rows / 2), 0,
                     rotation=90)
        for player in range(1, number_of_players):
            depth = player - 1
            plot += draw(player_labels[player],
                         (table_left / 2, bands - depth - 0.5), player)
            # The column profiles run in lexicographic order, so this player's
            # strategy is constant on blocks of as many columns as the players
            # nested inside it have profiles between them.
            block = columns
            for k in shape[1:player + 1]:
                block //= k
            for start in range(0, columns, block):
                strategy = column_profiles[start][depth]
                plot += draw(strategy_labels[player][strategy],
                             ((start + block / 2) * cell_width,
                              bands - depth - 0.5), player)
        for row in range(rows):
            plot += draw(strategy_labels[0][row],
                         (table_left / 2, -row - 0.5), 0)

        # The payoffs of a cell are written out beside one another, each in its
        # player's color, and a best response is underlined in that color too.
        for row in range(rows):
            for column, column_profile in enumerate(column_profiles):
                profile = (row,) + column_profile
                pieces = pieces_of[profile]
                written = "".join(pieces)
                left = (column + 0.5) * cell_width - _label_width(written) / 2
                y = -row - 0.5
                for player, piece in enumerate(pieces):
                    plot += draw(piece, (left, y), player,
                                 horizontal_alignment='left')
                    if best_responses and (profile, player) in best:
                        # The comma of a piece is not part of the payoff.
                        payoff = _label_width(piece.rstrip(", "))
                        plot += line2d([(left, y - 0.25),
                                        (left + payoff, y - 0.25)],
                                       color=color_of_player[player],
                                       thickness=1, zorder=5)
                    left += _label_width(piece)

        # The labels drawn outside the table would otherwise hang over its edge.
        plot.axes(kwargs.get('axes', False))
        plot.set_aspect_ratio(1)
        plot.set_axes_range(table_left - 0.9, table_right + 0.3,
                            -rows - 0.3, bands + 0.3)
        plot._extra_kwds.update(kwargs)
        return plot

    def _pure_best_responses(self):
        r"""
        Return the pairs ``(profile, player)`` of the game at which ``player``
        is playing a best response to what the other players play in
        ``profile``.

        A strategy profile is a pure Nash equilibrium exactly when every player
        of the game is playing a best response at it, so this describes the
        pure equilibria as well.  Every strategy tying for the best payoff is
        returned, so a player may have more than one best response to the same
        choice of their opponents.

        This is used by :meth:`plot` to underline best responses and to shade
        pure Nash equilibria.

        EXAMPLES:

        In the prisoner's dilemma the first strategy is dominant for both
        players, so each of them is playing a best response wherever they play
        it, and the only profile at which both of them are -- hence the only
        pure Nash equilibrium -- is the one where both do::

            sage: A = matrix([[2, 5], [0, 4]])
            sage: B = matrix([[2, 0], [5, 4]])
            sage: prisoners_dilemma = NormalFormGame([A, B])
            sage: sorted(prisoners_dilemma._pure_best_responses())
            [((0, 0), 0), ((0, 0), 1), ((0, 1), 0), ((1, 0), 1)]

        Both players would nevertheless rather be at ``(1, 1)``, which is what
        makes the game a dilemma::

            sage: prisoners_dilemma[0, 0], prisoners_dilemma[1, 1]
            ([2, 2], [4, 4])

        Matching pennies has no pure Nash equilibrium: at every profile exactly
        one of the two players would rather have chosen otherwise::

            sage: A = matrix([[1, -1], [-1, 1]])
            sage: matching_pennies = NormalFormGame([A])
            sage: best = matching_pennies._pure_best_responses()
            sage: [profile for profile in matching_pennies
            ....:  if all((profile, player) in best for player in range(2))]
            []

        Ties are all counted, so in a game where a player is indifferent both
        of their strategies are best responses::

            sage: A = matrix([[1, 1], [0, 0]])
            sage: B = matrix([[3, 3], [3, 3]])
            sage: g = NormalFormGame([A, B])
            sage: sorted(profile for profile, player
            ....:        in g._pure_best_responses() if player == 1)
            [(0, 0), (0, 1), (1, 0), (1, 1)]

        TESTS:

        It works for games with more than two players::

            sage: import numpy as np
            sage: A = np.array([[[3, 1], [2, 3]], [[9, 3], [8, 2]]])
            sage: B = np.array([[[1, 5], [6, 5]], [[7, 2], [4, 6]]])
            sage: C = np.array([[[4, 9], [5, 8]], [[9, 3], [6, 4]]])
            sage: g = NormalFormGame([A, B, C])
            sage: best = g._pure_best_responses()
            sage: [profile for profile in sorted(g)
            ....:  if all((profile, player) in best for player in range(3))]
            [(0, 1, 1), (1, 0, 0)]
        """
        best = set()
        for player in range(len(self.players)):
            # Gather the profiles the other players cannot tell apart, and mark
            # those of them at which this player does as well as they can.
            grouped = {}
            for profile in self.utilities:
                others = profile[:player] + profile[player + 1:]
                grouped.setdefault(others, []).append(profile)
            for profiles in grouped.values():
                highest = max(self.utilities[profile][player]
                              for profile in profiles)
                best.update((profile, player) for profile in profiles
                            if self.utilities[profile][player] == highest)
        return best

    def _n_matrix_game(self, matrices):
        r"""
        Populate ``self.utilities`` from a list of payoff arrays, one per
        player.

        Each entry of ``matrices`` is the payoff array of a player: a Sage
        matrix for a 2 player game, or an `N`-dimensional numpy array for an
        `N` player game. All arrays must have the same shape, and the number
        of arrays must equal their common number of dimensions (one strategy
        axis per player).

        EXAMPLES:

        A small two player game::

            sage: A = matrix([[1, 0], [-2, 3]])
            sage: B = matrix([[3, 2], [-1, 0]])
            sage: two_game = NormalFormGame()
            sage: two_game._n_matrix_game([A, B])
            sage: two_game
            Normal Form Game with the following utilities: {(0, 0): [1, 3],
             (0, 1): [0, 2], (1, 0): [-2, -1], (1, 1): [3, 0]}

        A three player game built from three 3-dimensional numpy arrays::

            sage: import numpy as np
            sage: A = np.array([[[3, 1], [2, 3]], [[9, 3], [8, 2]]])
            sage: B = np.array([[[1, 5], [6, 5]], [[7, 2], [4, 6]]])
            sage: C = np.array([[[4, 9], [5, 8]], [[9, 3], [6, 4]]])
            sage: three_game = NormalFormGame()
            sage: three_game._n_matrix_game([A, B, C])
            sage: three_game
            Normal Form Game with the following utilities: {(0, 0, 0): [3, 1, 4],
             (0, 0, 1): [1, 5, 9],
             (0, 1, 0): [2, 6, 5],
             (0, 1, 1): [3, 5, 8],
             (1, 0, 0): [9, 7, 9],
             (1, 0, 1): [3, 2, 3],
             (1, 1, 0): [8, 4, 6],
             (1, 1, 1): [2, 6, 4]}
        """
        self.players = []
        self.utilities = {}

        def shape(m):
            # Sage matrices expose .dimensions(); numpy arrays expose .shape
            if hasattr(m, "dimensions"):
                return tuple(m.dimensions())
            return tuple(m.shape)

        shapes = [shape(m) for m in matrices]
        if any(s != shapes[0] for s in shapes):
            raise ValueError("matrices must be the same size")
        if len(matrices) != len(shapes[0]):
            raise ValueError("the number of matrices must match the number of "
                             "players, i.e. the number of dimensions of each "
                             "matrix")

        for num_strategies in shapes[0]:
            self.add_player(num_strategies)

        for strategy_profile in self.utilities:
            utility_vector = []
            for m in matrices:
                value = m[strategy_profile]
                if np is not None and isinstance(value, np.generic):
                    value = value.item()
                utility_vector.append(value)
            self.utilities[strategy_profile] = utility_vector

    def _gambit_game(self, game):
        r"""
        Create a ``NormalFormGame`` object from a Gambit game.

        TESTS::

            sage: # optional - pygambit
            sage: import numpy as np
            sage: from pygambit import Game
            sage: testgame = Game.from_arrays(np.array([[8.5, 2.5], [10.1, 5.1]]),
            ....:                             np.array([[8.5, 10.1], [2.5, 5.1]]))
            sage: g = NormalFormGame()
            sage: g._gambit_game(testgame); g
            Normal Form Game with the following utilities: {(0, 0): [8.5, 8.5],
            (0, 1): [2.5, 10.1],
            (1, 0): [10.1, 2.5],
            (1, 1): [5.1, 5.1]}

        ::

            sage: # optional - pygambit
            sage: import numpy as np
            sage: from pygambit import Game
            sage: testgame = Game.from_arrays(np.array([[3, 0], [5, 10]], dtype=int),
            ....:                             np.array([[3, 5], [0, 10]], dtype=int))
            sage: g = NormalFormGame()
            sage: g._gambit_game(testgame); g
            Normal Form Game with the following utilities: {(0, 0): [3, 3],
            (0, 1): [0, 5],
            (1, 0): [5, 0],
            (1, 1): [10, 10]}
        """
        self.players = []
        self.utilities = {}
        for player in game.players:
            num_strategies = len(player.strategies)
            self.add_player(num_strategies)
        # gambit's player collection is indexed by label, not by position
        players = list(game.players)
        for strategy_profile in self.utilities:
            outcome = game[strategy_profile]
            # gambit stores inexact payoffs as decimals and exact ones as rationals
            utility_vector = [float(payoff) if isinstance(payoff, Decimal) else QQ(payoff)
                              for payoff in (outcome[p] for p in players)]
            self.utilities[strategy_profile] = utility_vector

    def _gambit_(self, maximization=True):
        r"""
        Create a Gambit game from a ``NormalFormGame`` object.

        Payoffs are handed to gambit one by one, so each keeps the kind of
        value it is: an exact Sage :class:`~sage.rings.integer.Integer` or
        :class:`~sage.rings.rational.Rational` becomes a gambit ``Rational``
        and an inexact real becomes a gambit ``Decimal``.  Nothing is rounded
        or truncated on the way.

        INPUT:

        - ``maximization`` -- boolean; whether a player is trying to
          maximize their utility or minimize it

        TESTS:

        Gambit games are indexed by strategy labels (and players by their
        labels) rather than by integers; the default labels are the strings
        ``'1'``, ``'2'``, ... ::

            sage: # optional - pygambit
            sage: A = matrix([[2, 1], [1, 2.5]])
            sage: g = NormalFormGame([A])
            sage: gg = g._gambit_()
            sage: print(gg.to_nfg())
            NFG 1 R "Untitled strategic game" { "1" "2" }
            <BLANKLINE>
            { { "1" "2" }
            { "1" "2" }
            }
            ""
            <BLANKLINE>
            {
            { "" 2.0, -2.0 }
            { "" 1.0, -1.0 }
            { "" 1.0, -1.0 }
            { "" 2.5, -2.5 }
            }
            1 2 3 4
            <BLANKLINE>

        ::

            sage: # optional - pygambit
            sage: A = matrix([[2, 1], [1, 2.5]])
            sage: B = matrix([[3, 2], [5.5, 4]])
            sage: g = NormalFormGame([A, B])
            sage: gg = g._gambit_()
            sage: float(gg['1', '1'][gg.players['1']])
            2.0
            sage: float(gg['2', '1'][gg.players['2']])
            5.5

        Exact payoffs stay exact.  A rational payoff reaches gambit as a
        gambit ``Rational``, not as a decimal approximation of one::

            sage: # optional - pygambit
            sage: A = matrix(QQ, [[1/3, 2], [1, 5/7]])
            sage: gg = NormalFormGame([A])._gambit_()
            sage: gg['1', '1'][gg.players['1']]
            Rational(1, 3)
            sage: [line for line in gg.to_nfg().splitlines()
            ....:  if line.startswith('{ ""')]
            ['{ "" 1/3, -1/3 }', '{ "" 1, -1 }', '{ "" 2, -2 }', '{ "" 5/7, -5/7 }']

        Nothing is truncated either; the payoff ``27/10`` is not turned into
        ``2``::

            sage: # optional - pygambit
            sage: gg = NormalFormGame([matrix(QQ, [[27/10, 1], [1, 2]])])._gambit_()
            sage: gg['1', '1'][gg.players['1']]
            Rational(27, 10)

        ::

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
            sage: gg = threegame._gambit_()                      # optional - pygambit
            sage: int(gg['1', '1', '1'][gg.players['1']])         # optional - pygambit
            3
            sage: int(gg['2', '1', '1'][gg.players['2']])         # optional - pygambit
            7
            sage: int(gg['1', '1', '2'][gg.players['3']])         # optional - pygambit
            9

        An unpopulated payoff has no gambit counterpart -- it would be written
        out as a zero, silently turning the game into a different one -- so a
        game whose utilities are incomplete is refused rather than converted::

            sage: # optional - pygambit
            sage: A = matrix([[2, 1], [1, 2.5]])
            sage: g = NormalFormGame([A])
            sage: g.add_strategy(0)
            sage: g._gambit_()
            Traceback (most recent call last):
            ...
            ValueError: utilities have not been populated; ...
        """
        self._require_complete()

        sgn = 1 if maximization else -1
        strategy_sizes = [p.num_strategies for p in self.players]
        n_players = len(strategy_sizes)

        def payoff(value):
            """
            Hand ``value`` to gambit without losing anything.

            An exact value is passed through, so that gambit stores it as a
            rational.  An inexact one is narrowed to a Python float where that
            is lossless, purely so that a written game reads ``2.5`` rather
            than ``2.50000000000000``; a value carrying more precision than a
            float can hold is passed through untouched rather than cut down.
            """
            try:
                ring = value.parent()
            except AttributeError:      # a plain Python number
                return value
            if ring.is_exact():
                return value
            as_float = float(value)
            return as_float if ring(as_float) == value else value

        # The arrays hold the Sage values themselves rather than machine
        # numbers, so that gambit converts each payoff on its own terms and
        # exact ones stay exact; a numeric dtype would force the lot through
        # floating point.  Leaving the cells uninitialised is safe because
        # ``_require_complete`` above has checked that every one is written.
        arrays = [np.empty(strategy_sizes, dtype=object) for _ in range(n_players)]
        for sp in self.utilities:
            for i in range(n_players):
                arrays[i][sp] = payoff(sgn * self.utilities[sp][i])

        return Game.from_arrays(*arrays)

    def save_nfg(self, path):
        r"""
        Save the game to ``path`` in Gambit's strategic-form ``.nfg`` format.

        The game is converted to a Gambit game (see :meth:`_gambit_`),
        serialised with Gambit's writer and written atomically with
        :func:`~sage.misc.temporary_file.atomic_write` so that a partially
        written file is never left behind.

        The game is always written in the ``.nfg`` (strategic form) format.
        A :class:`NormalFormGame` is always a strategic-form game, and Gambit's
        extensive-form (``.efg``) writer is undefined for such games, so that
        format is not available.

        INPUT:

        - ``path`` -- string; the file path to write the game to

        EXAMPLES::

            sage: A = matrix([[2, 1], [1, 2.5]])
            sage: g = NormalFormGame([A])
            sage: path = tmp_filename(ext='.nfg')
            sage: g.save_nfg(path)                           # optional - pygambit
            sage: with open(path) as f:                      # optional - pygambit
            ....:     print(f.read()[:5])
            NFG 1

        A game whose utilities are not fully populated is refused, so the file
        never records payoffs that were never set (see :meth:`_gambit_`).  The
        error comes before anything is written, leaving the file untouched::

            sage: # optional - pygambit
            sage: h = NormalFormGame([A])
            sage: h.add_strategy(0)
            sage: path = tmp_filename(ext='.nfg')
            sage: h.save_nfg(path)
            Traceback (most recent call last):
            ...
            ValueError: utilities have not been populated; ...
            sage: import os; os.path.getsize(path)
            0
        """
        pygambit().require()
        g = self._gambit_()
        with atomic_write(path) as f:   # text mode by default (binary=False)
            f.write(g.to_nfg())

    @classmethod
    def load_nfg(cls, path):
        r"""
        Read a game from a gambit strategic-form ``.nfg`` file.

        The file at ``path`` is read with gambit's ``read_nfg`` reader and the
        resulting gambit game is converted into a new
        :class:`NormalFormGame` (see :meth:`_gambit_game`).  This is the
        inverse of :meth:`save_nfg`.

        INPUT:

        - ``path`` -- string; the path of a ``.nfg`` file to read

        OUTPUT: a new :class:`NormalFormGame` holding the game in the file

        EXAMPLES:

        A game can be saved and then read back in::

            sage: # optional - pygambit
            sage: A = matrix([[2, 1], [1, 2.5]])
            sage: B = matrix([[4, 3], [2, 1]])
            sage: g = NormalFormGame([A, B])
            sage: path = tmp_filename(ext='.nfg')
            sage: g.save_nfg(path)
            sage: h = NormalFormGame.load_nfg(path); h
            Normal Form Game with the following utilities: {(0, 0): [2.0, 4], (0, 1): [1.0, 3], (1, 0): [1.0, 2], (1, 1): [2.5, 1]}

        A file that gambit cannot parse is reported as such::

            sage: # optional - pygambit
            sage: path = tmp_filename(ext='.nfg')
            sage: with open(path, 'w') as f:
            ....:     _ = f.write('not an nfg file')
            sage: NormalFormGame.load_nfg(path)
            Traceback (most recent call last):
            ...
            ValueError: Parse error in game file: ...
        """
        pygambit().require()
        return cls(read_nfg(path))

    @classmethod
    def gambit_catalog_games(cls, **kwargs):
        r"""
        List the games in the gambit catalog.

        The `gambit catalog
        <https://gambitproject.readthedocs.io/en/stable/catalog.html>`_ ships a
        small collection of example games from the literature.  All of them are
        listed here: the strategic form games directly, and the extensive form
        ones because :meth:`load_from_gambit_catalog` converts a tree to its
        reduced strategic form.
        :meth:`~sage.game_theory.extensive_form_game.ExtensiveFormGame.gambit_catalog_games`
        lists just the trees, which are the ones that keep their tree structure.

        INPUT:

        - ``**kwargs`` -- passed on to ``pygambit.catalog.games``, whose
          keywords filter the listing, for instance ``n_players``, ``is_tree``,
          ``is_const_sum`` or ``include_descriptions``

        OUTPUT: a :class:`pandas.DataFrame` with a ``Game`` column holding the
        slugs and a ``Title`` column holding the titles.  Slugs are full paths,
        such as ``'journals/geb/bagwell1995'`` or ``'books/myerson1991/fig2_1'``.

        EXAMPLES::

            sage: # optional - pygambit
            sage: games = NormalFormGame.gambit_catalog_games()
            sage: 'journals/geb/bagwell1995' in list(games['Game'])
            True
            sage: list(games.columns)
            ['Game', 'Title']

        The keywords of gambit's own catalog listing filter the table::

            sage: # optional - pygambit
            sage: trees = NormalFormGame.gambit_catalog_games(is_tree=True)
            sage: len(trees) < len(games)
            True
        """
        pygambit().require()
        return catalog.games(**kwargs)

    @classmethod
    def load_from_gambit_catalog(cls, slug):
        r"""
        Load a game from the gambit catalog.

        The available games are listed by :meth:`gambit_catalog_games`.  A
        catalog game given in extensive form is converted to its reduced
        strategic form, via :meth:`_gambit_game`.

        INPUT:

        - ``slug`` -- string; the slug of a catalog game, a full path such as
          ``'journals/geb/bagwell1995'``

        OUTPUT: a new :class:`NormalFormGame` holding the catalog game

        EXAMPLES::

            sage: # optional - pygambit
            sage: g = NormalFormGame.load_from_gambit_catalog(
            ....:     'journals/geb/bagwell1995')
            sage: len(g.players)
            2

        ``'journals/geb/bagwell1995'`` is a tree game in the catalog; loading
        it here gives its reduced strategic form, so each player's strategies
        are complete contingent plans rather than single moves::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: t = ExtensiveFormGame.load_from_gambit_catalog(
            ....:     'journals/geb/bagwell1995')
            sage: [p.num_strategies for p in g.players]
            [2, 4]
            sage: len(t.infosets)
            3

        A slug that is not in the catalog is rejected::

            sage: # optional - pygambit
            sage: NormalFormGame.load_from_gambit_catalog('not_a_real_game')
            Traceback (most recent call last):
            ...
            ValueError: 'not_a_real_game' is not a game in the gambit catalog; ...
        """
        pygambit().require()
        try:
            loaded = catalog.load(slug)
        except FileNotFoundError:
            raise ValueError(
                f"{slug!r} is not a game in the gambit catalog; call "
                "gambit_catalog_games() to see the available games"
            )
        return cls(loaded)

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
        t = c[0, 0]

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
            ValueError: utilities have not been populated; ...

        The above creates a 2 player game where each player has
        a single strategy. Here we populate the strategies and
        can then view the payoff matrices::

            sage: g[0, 0] = [1,2]
            sage: g.payoff_matrices()
            ([1], [2])
        """
        if len(self.players) != 2:
            raise ValueError("Only available for 2 player games")

        self._require_complete()

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

    def _require_complete(self):
        r"""
        Raise a :class:`ValueError` unless every strategy profile has a payoff.

        Adding a player or a strategy leaves the new profiles unpopulated, and
        everything that reads the payoffs off the game -- converting it to a
        gambit game, forming the payoff matrices, solving it or plotting it --
        needs all of them.  This is the one place that check is phrased.

        .. SEEALSO:: :meth:`_is_complete`

        EXAMPLES::

            sage: s = matrix([[1, 0], [-2, 3]])
            sage: t = matrix([[3, 2], [-1, 0]])
            sage: example = NormalFormGame([s, t])
            sage: example._require_complete()
            sage: example.add_strategy(0)
            sage: example._require_complete()
            Traceback (most recent call last):
            ...
            ValueError: utilities have not been populated; ...
        """
        if not self._is_complete():
            raise ValueError(
                "utilities have not been populated; set a payoff value for "
                "every strategy profile, e.g. game[0, 1][0] = 3, or "
                "construct the game from payoff matrices via "
                "NormalFormGame([A, B])"
            )

    def obtain_nash(self, algorithm=False, maximization=True, solver=None,
                    phc_path=None, rational=True, tolerance=1e-4):
        r"""
        A function to return the Nash equilibrium for the game.
        Optional arguments can be used to specify the algorithm used.
        If no algorithm is passed then an attempt is made to use the most
        appropriate algorithm.

        INPUT:

        - ``algorithm`` -- the following algorithms should be available through
          this function:

          * ``'lrs'`` -- this algorithm is only suited for 2 player games.
            See the lrs web site (https://cgm.cs.mcgill.ca/~avis/C/lrs.html).

          * ``'enummixed'`` -- this algorithm is only suited for 2 player games. It
            computes all mixed strategy Nash equilibria based on the gambit implementation,
            see the gambit web site (https://gambitproject.readthedocs.io/en/stable/api/pygambit.nash.enummixed_solve.html).

          * ``'lcp'`` -- this algorithm is only suited for 2 player games.
            See the gambit web site (https://gambitproject.readthedocs.io/en/stable/api/pygambit.nash.lcp_solve.html).

          * ``'gnm'``, ``'enumpure'``, ``'enumpoly'``, ``'liap'``,
            ``'simpdiv'``, ``'ipa'``, ``'logit'`` -- these algorithms are
            suited for games with an arbitrary number of players (in
            particular for games with more than 2 players) and use the
            corresponding solver implemented in gambit (the global Newton
            method, enumeration of pure equilibria, enumeration of equilibria
            via polynomial systems, Lyapunov function minimisation,
            simplicial subdivision, iterated polymatrix approximation and
            logit quantal response tracing respectively). They are numerical
            algorithms and so in general return floating point approximations
            of a sample of the equilibria. See the gambit web site
            (https://gambitproject.readthedocs.io/en/stable/pygambit.api.html). When no ``algorithm`` is given
            for a game with more than 2 players, ``'enumpoly'`` is used.

          * ``'lp'`` -- this algorithm is only suited for 2 player
            constant sum games. Uses MILP solver or the gambit solver, determined by the
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

          The name of the algorithm is matched without regard to case, so
          ``'LCP'`` and ``'lcp'`` name the same solver.

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

        - ``rational`` -- boolean (default: ``True``); whether to answer with
          rational probabilities rather than floating point ones, which is done
          in whichever of two ways the algorithm allows.  ``'lcp'``, ``'lp'``
          and ``'enummixed'`` are asked to compute in exact arithmetic
          throughout; ``'gnm'``, ``'ipa'``, ``'logit'``, ``'liap'`` and
          ``'enumpoly'`` have no exact mode, so they compute in floating point
          and their answer is then rounded to the exact equilibrium it
          approximates, as described under ``tolerance``.  ``'enumpure'``,
          ``'simpdiv'`` and Sage's own ``'enumeration'`` and ``'lrs'`` always
          work exactly and ignore the flag.

          Note that exact arithmetic does not by itself make an answer exact:
          ``'simpdiv'`` refines a grid, so it returns a rational that
          approximates an equilibrium.

        - ``tolerance`` -- a positive number (default: ``1e-4``); how far a
          probability computed by one of the numerical algorithms may be moved
          in order to round it to a rational.  A rounded profile is returned
          only once gambit has confirmed, in exact arithmetic, that it is an
          equilibrium; when it is not -- as happens when the equilibrium is
          genuinely irrational, which a game of three or more players may have
          -- the floating point answer is returned instead.  The larger the
          tolerance the simpler the rationals that are tried, so a value that is
          too small can round an equilibrium to an unenlightening exact form
          rather than fail to round it; the default matches the accuracy the
          gambit solvers themselves aim for.  Has no effect when ``rational``
          is ``False``, nor on a game whose own payoffs are inexact.

        - ``phc_path`` -- (optional) a path (a string or
          :class:`~pathlib.Path`) to the PHCpack ``phc`` executable. When
          given, the ``'enumpoly'`` algorithm solves the underlying systems of
          polynomial equations with PHCpack instead of gambit's built-in
          solver. This is only supported by the ``'enumpoly'`` algorithm (the
          default for games with more than 2 players); passing it for any other
          algorithm raises a :class:`ValueError`. PHCpack is available from
          https://homepages.math.uic.edu/~jan/download.html.

        OUTPUT:

        A sorted list holding one *mixed strategy profile* per computed
        equilibrium.  Each profile is a list with one entry per player, in the
        order the players were added to the game, and each entry is a tuple
        giving the probability that player assigns to each of their pure
        strategies, in the order of that player's strategies.  So for a
        2-player game ``[[(1/2, 1/2), (1/3, 2/3)]]`` is a single equilibrium in
        which the first player mixes uniformly over two strategies and the
        second plays their strategies with probabilities `1/3` and `2/3`.
        The probabilities of each player sum to one.

        A probability that was found exactly is an element of
        :class:`~sage.rings.rational.Rational` and one that was not is a Python
        float; see ``rational`` and ``tolerance`` above for when each occurs.
        The numerical gambit solvers may in addition find only a sample of the
        equilibria rather than all of them.

        .. SEEALSO:: :meth:`_extract_gambit_equilibria`, which converts gambit's
            own profile objects into this format.

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
            sage: g.obtain_nash(algorithm='lcp')  # optional - pygambit
            [[(0, 0, 3/4, 1/4), (1/28, 27/28, 0)]]

        The ``'enummixed'`` algorithm (2 player games only) enumerates all the
        extreme mixed strategy Nash equilibria; for a coordination game it
        returns the two pure equilibria together with the mixed one::

            sage: A = matrix([[1, 0], [0, 1]])
            sage: coordination = NormalFormGame([A, A])
            sage: coordination.obtain_nash(algorithm='enummixed')  # optional - pygambit
            [[(0, 1), (0, 1)], [(1/2, 1/2), (1/2, 1/2)], [(1, 0), (1, 0)]]

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
            sage: fivegame.obtain_nash(algorithm='lcp')  # optional - pygambit
            [[(1, 0, 0, 0, 0), (0, 1, 0, 0, 0)]]

        Here are some examples of finding Nash equilibria for constant-sum games::

            sage: A = matrix.identity(2)
            sage: cg = NormalFormGame([A])
            sage: cg.obtain_nash(algorithm='lp')
            [[(0.5, 0.5), (0.5, 0.5)]]
            sage: cg.obtain_nash(algorithm='lp', solver='Coin')                   # optional - sage_numerical_backends_coin
            [[(0.5, 0.5), (0.5, 0.5)]]
            sage: cg.obtain_nash(algorithm='lp', solver='PPL')
            [[(1/2, 1/2), (1/2, 1/2)]]
            sage: cg.obtain_nash(algorithm='lp', solver='gambit')                 # optional - pygambit
            [[(1/2, 1/2), (1/2, 1/2)]]
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
            sage: ne = cg.obtain_nash(algorithm='lp', solver='gambit')            # optional - pygambit
            sage: [[[round(el, 6) for el in v] for v in eq] for eq in ne]         # optional - pygambit
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
            sage: ne = cg.obtain_nash(algorithm='lp', solver='gambit')            # optional - pygambit
            sage: [[[round(el, 6) for el in v] for v in eq] for eq in ne]         # optional - pygambit
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
            ValueError: the 'lp' algorithm only works for two-player constant-sum games, ...

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
        is possible for the ``'lcp'`` solver to find all Nash equilibria
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
            sage: gg.obtain_nash(algorithm='lcp')  # optional - pygambit
            [[(1, 0), (1, 0)]]
            sage: gg.obtain_nash(algorithm='enumeration', maximization=False)
            [[(0, 1), (0, 1)], [(0, 1), (1, 0)], [(1, 0), (0, 1)], [(1, 0), (1, 0)]]
            sage: gg.obtain_nash(algorithm='lrs', maximization=False)  # optional - lrs
            [[(0, 1), (0, 1)], [(0, 1), (1, 0)], [(1, 0), (0, 1)], [(1, 0), (1, 0)]]
            sage: gg.obtain_nash(algorithm='lp', solver='glpk', maximization=False)
            [[(1.0, 0.0), (1.0, 0.0)]]
            sage: gg.obtain_nash(algorithm='lcp', maximization=False)  # optional - pygambit
            [[(1, 0), (1, 0)]]

        Note that outputs for all algorithms are as lists of lists of
        tuples and the equilibria have been sorted so that all algorithms give
        a comparable output (although ``'lcp'`` returns floats)::

            sage: enumeration_eqs = g.obtain_nash(algorithm='enumeration')
            sage: [[type(s) for s in eq] for eq in enumeration_eqs]
            [[<... 'tuple'>, <... 'tuple'>], [<... 'tuple'>, <... 'tuple'>], [<... 'tuple'>, <... 'tuple'>]]
            sage: lrs_eqs = g.obtain_nash(algorithm='lrs')  # optional - lrslib
            sage: [[type(s) for s in eq] for eq in lrs_eqs]  # optional - lrslib
            [[<... 'tuple'>, <... 'tuple'>], [<... 'tuple'>, <... 'tuple'>], [<... 'tuple'>, <... 'tuple'>]]
            sage: lcp_eqs = g.obtain_nash(algorithm='lcp')  # optional - pygambit
            sage: [[type(s) for s in eq] for eq in lcp_eqs]  # optional - pygambit
            [[<... 'tuple'>, <... 'tuple'>], [<... 'tuple'>, <... 'tuple'>], [<... 'tuple'>, <... 'tuple'>]]
            sage: enumeration_eqs == sorted(enumeration_eqs)
            True
            sage: lrs_eqs == sorted(lrs_eqs)  # optional - lrslib
            True
            sage: lcp_eqs == sorted(lcp_eqs)  # optional - pygambit
            True
            sage: lrs_eqs == enumeration_eqs  # optional - lrslib
            True
            sage: enumeration_eqs == lcp_eqs  # optional - pygambit
            True
            sage: [[[round(float(p), 6) for p in str] for str in eq] for eq in enumeration_eqs] == [[[round(float(p), 6) for p in str] for str in eq] for eq in lcp_eqs]  # optional - pygambit
            True

        The :math:`3\times 3` coordination game (both payoff matrices the
        identity) has :math:`2^3 - 1 = 7` equilibria, one uniform mix over
        each nonempty set of matching strategies.  ``'enumpoly'`` finds all
        seven, ``'enumpure'`` finds the three pure ones and ``'simpdiv'``
        returns the fully mixed one::

            sage: I3 = matrix.identity(3)
            sage: coordination = NormalFormGame([I3, I3])
            sage: coordination.obtain_nash(algorithm='enumpure')  # optional - pygambit
            [[(0, 0, 1), (0, 0, 1)], [(0, 1, 0), (0, 1, 0)], [(1, 0, 0), (1, 0, 0)]]
            sage: coordination.obtain_nash(algorithm='enumpoly')  # optional - pygambit
            [[(0, 0, 1), (0, 0, 1)],
             [(0, 1/2, 1/2), (0, 1/2, 1/2)],
             [(0, 1, 0), (0, 1, 0)],
             [(1/3, 1/3, 1/3), (1/3, 1/3, 1/3)],
             [(1/2, 0, 1/2), (1/2, 0, 1/2)],
             [(1/2, 1/2, 0), (1/2, 1/2, 0)],
             [(1, 0, 0), (1, 0, 0)]]
            sage: coordination.obtain_nash(algorithm='simpdiv')  # optional - pygambit
            [[(1/3, 1/3, 1/3), (1/3, 1/3, 1/3)]]

        A 3 player game, coming from a local max cut instance, with two pure equilibria where
        either player is on one side of the cut and one mixed equilibrium where all players mix uniformly::

            sage: import numpy as np
            sage: A = np.array([[[0, -1], [2, 1]], [[1, 2],[-1, 0]]])
            sage: B = np.array([[[0, 2], [4, 2]], [[2, 4], [2, 0]]])
            sage: C = np.array([[[0, 1], [2, -1]], [[-1, 2], [1, 0]]])
            sage: max_cut_game = NormalFormGame([A, B, C])
            sage: max_cut_game.obtain_nash(algorithm='enumpoly')  # optional - pygambit
            [[(0, 1), (1, 0), (0, 1)],
             [(1/2, 1/2), (1/2, 1/2), (1/2, 1/2)],
             [(1, 0), (0, 1), (1, 0)]]
            sage: max_cut_game.obtain_nash(algorithm='enumpure')  # optional - pygambit
            [[(0, 1), (1, 0), (0, 1)], [(1, 0), (0, 1), (1, 0)]]

        The solvers with no exact mode answer exactly as well, their floating
        point answer being rounded back to the equilibrium it approximates.
        Here ``'liap'`` minimises the Lyapunov function of the same game and
        lands on the mixed equilibrium::

            sage: max_cut_game.obtain_nash(algorithm='liap')  # optional - pygambit
            [[(1/2, 1/2), (1/2, 1/2), (1/2, 1/2)]]

        Rounding is only accepted once gambit has confirmed, in exact
        arithmetic, that the rounded profile is an equilibrium, so it never
        passes off an approximation as exact.  Ask for a rounding so fine that
        nothing verifies and the floating point answer comes back unchanged::

            sage: # optional - pygambit
            sage: battle = NormalFormGame([matrix([[3, 0], [0, 2]]),
            ....:                          matrix([[2, 0], [0, 3]])])
            sage: battle.obtain_nash(algorithm='liap')
            [[(3/5, 2/5), (2/5, 3/5)]]
            sage: battle.obtain_nash(algorithm='liap', tolerance=1e-15)  # abs tol 1e-6
            [[(0.6, 0.4), (0.4, 0.6)]]

        That is what happens of its own accord when an equilibrium is
        irrational, which a game of three or more players may well have.  The
        one below is `((101 - \sqrt{3129})/136, \ldots)`, so no rounding of it
        can be exact and ``'liap'`` reports what it found::

            sage: A = np.array([[[-5, -3], [2, -1]], [[1, -4], [1, -3]]])
            sage: B = np.array([[[2, 4], [-5, -3]], [[1, -3], [3, 1]]])
            sage: C = np.array([[[0, -4], [-2, -2]], [[2, 2], [-1, 3]]])
            sage: irrational = NormalFormGame([A, B, C])
            sage: irrational.obtain_nash(algorithm='liap')  # abs tol 1e-5 # optional - pygambit
            [[(0.331342, 0.668658), (0.668658, 0.331342), (0.265634, 0.734366)]]
            sage: float((101 - sqrt(3129)) / 136)                                   # needs sage.symbolic
            0.33134216850155895

        Also, not specifying a valid solver would lead to an error::

            sage: A = matrix.identity(2)
            sage: g = NormalFormGame([A])
            sage: g.obtain_nash(algorithm='invalid')
            Traceback (most recent call last):
            ...
            ValueError: unknown algorithm 'invalid' for a 2-player game; ...

        The name of the algorithm is not case-sensitive, and an unknown one is
        quoted back the way it was spelled::

            sage: g.obtain_nash(algorithm='Enumeration')
            [[(1/2, 1/2), (1/2, 1/2)]]
            sage: g.obtain_nash(algorithm='Invalid')
            Traceback (most recent call last):
            ...
            ValueError: unknown algorithm 'Invalid' for a 2-player game; ...
            sage: g.obtain_nash(algorithm='lp', solver='invalid')
            Traceback (most recent call last):
            ...
            ValueError: 'solver' should be set to 'GLPK', ..., None
             (in which case the default one is used), or a callable.

        Passing ``phc_path`` with an algorithm other than ``'enumpoly'``
        raises an error::

            sage: g.obtain_nash(algorithm='lrs', phc_path='/usr/bin/phc')
            Traceback (most recent call last):
            ...
            ValueError: 'phc_path' is only supported by the 'enumpoly' algorithm; got algorithm 'lrs'
        """
        self._require_complete()

        from sage.features.lrs import LrsNash
        # The algorithm is matched without regard to case, but the error
        # messages below quote the name back the way the caller spelled it.
        requested = algorithm
        if not algorithm:
            if len(self.players) > 2:
                # Only the gambit solvers handle games with more than two
                # players; enumpoly enumerates the equilibria by solving the
                # corresponding systems of polynomial equations.
                algorithm = "enumpoly"
            elif self.is_constant_sum():
                algorithm = "lp"
            elif LrsNash().is_present():
                algorithm = "lrs"
            else:
                algorithm = "enumeration"
            requested = algorithm
        else:
            algorithm = algorithm.lower()

        if phc_path is not None and algorithm != "enumpoly":
            raise ValueError(
                "'phc_path' is only supported by the 'enumpoly' algorithm; "
                f"got algorithm {requested!r}"
            )

        if len(self.players) < 3:
            if algorithm == "lrs":
                LrsNash().require()
                return self._solve_lrs(maximization)

            if algorithm == "lcp":
                pygambit().require()
                return self._use_gambit_solver('lcp', maximization,
                                               rational=rational,
                                               tolerance=tolerance)

            if algorithm == 'lp':
                return self._solve_LP(solver=solver, maximization=maximization,
                                      rational=rational)

            if algorithm == "enumeration":
                return self._solve_enumeration(maximization)

            if algorithm == "enummixed":
                pygambit().require()
                return self._use_gambit_solver('enummixed', maximization,
                                               rational=rational,
                                               tolerance=tolerance)

        # The remaining gambit solvers all handle an arbitrary number of
        # players, so they are routed here (outside the two player branch
        # above) and are available for two player games as well.
        gambit_algorithms = {"gnm", "enumpure", "enumpoly", "liap",
                             "simpdiv", "ipa", "logit"}
        if algorithm in gambit_algorithms:
            pygambit().require()
            return self._use_gambit_solver(algorithm, maximization,
                                           phc_path=phc_path,
                                           rational=rational,
                                           tolerance=tolerance)

        n = len(self.players)
        raise ValueError(
            f"unknown algorithm {requested!r} for a {n}-player game; "
            "for 2-player games use 'enumeration', 'lrs', 'lcp', 'lp', or "
            "'enummixed', and for any number of players use one of the gambit "
            "solvers: 'gnm', 'enumpure', 'enumpoly', 'liap', 'simpdiv', 'ipa', 'logit'"
        )

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

    def _extract_gambit_equilibria(self, equilibria):
        r"""
        Convert a pygambit equilibria collection to a list of strategy profiles.

        The Gambit solvers (see :meth:`_use_gambit_solver`) return their
        results as pygambit objects.  Each
        equilibrium ``eq`` is a mixed strategy profile that behaves like a
        mapping from a player's strategies to the probability with which that
        strategy is played.  This helper flattens that representation into the
        format used throughout the rest of :class:`NormalFormGame`: a list of
        equilibria, where each equilibrium is a list with one entry per player
        and each entry is a tuple giving the probability assigned to each of
        that player's pure strategies.

        Concretely, the returned object has the shape::

            [
              [ (p_0^0, p_0^1, ...),   # mixed strategy of player 0
                (p_1^0, p_1^1, ...),   # mixed strategy of player 1
                ... ],                 # one tuple per player
              ...                      # one such list per equilibrium
            ]

        where `p_i^j` is the probability that player `i` plays their pure
        strategy `j`.  A probability gambit computed exactly is returned as a
        Sage :class:`~sage.rings.rational.Rational` and an inexact one as a
        Python ``float``, mirroring the conversion :meth:`_gambit_game` makes
        for payoffs.

        INPUT:

        - ``equilibria`` -- an iterable of pygambit mixed strategy profiles,
          such as the ``equilibria`` attribute of the result returned by a
          Gambit solver

        OUTPUT: a list of equilibria, each represented as a list of tuples of
        floats (one tuple per player)

        EXAMPLES:

        Solving a two player game with Gambit's GNM solver and then extracting
        the equilibria into Sage's native format::

            sage: from pygambit.nash import gnm_solve            # optional - pygambit
            sage: A = matrix([[1, 2], [3, 4]])
            sage: B = matrix([[3, 3], [1, 4]])
            sage: C = NormalFormGame([A, B])
            sage: g = C._gambit_()                               # optional - pygambit
            sage: result = gnm_solve(g)                          # optional - pygambit
            sage: C._extract_gambit_equilibria(result.equilibria)  # optional - pygambit
            [[(0.0, 1.0), (0.0, 1.0)]]

        The method works for any number of players; here each of the three
        players is given a single tuple in every equilibrium::

            sage: threegame = NormalFormGame()
            sage: threegame.add_player(2)
            sage: threegame.add_player(2)
            sage: threegame.add_player(2)
            sage: for i in range(2):
            ....:     for j in range(2):
            ....:         for k in range(2):
            ....:             threegame[i, j, k][0] = i + 1
            ....:             threegame[i, j, k][1] = j + 1
            ....:             threegame[i, j, k][2] = k + 1
            sage: g = threegame._gambit_()                       # optional - pygambit
            sage: eqs = threegame._extract_gambit_equilibria(gnm_solve(g).equilibria)  # optional - pygambit
            sage: all(len(eq) == 3 for eq in eqs)                # optional - pygambit
            True
        """
        # gambit hands back rationals from the exact solvers and floats from
        # the numerical ones; keep whichever it produced rather than flattening
        # an exact answer down to a machine number.
        def entry(value):
            return QQ(value) if isinstance(value, Rational) else float(value)

        return [[tuple(entry(eq[s]) for s in player.strategies)
                 for player in eq.game.players]
                for eq in equilibria]

    def _use_gambit_solver(self, algorithm, maximization=True, phc_path=None,
                           rational=True, tolerance=1e-4):
        r"""
        Solve a :class:`NormalFormGame` using one of Gambit's solvers.

        This is a single entry point for the Gambit ([Gambit]_) Nash
        solvers.  The desired solver is selected through the ``algorithm``
        argument so that the conversion to a Gambit game, the call to the
        solver and the extraction of the equilibria need only be written
        once.

        The ``'lcp'`` and ``'lp'`` solvers are restricted to two player
        games; all of the other solvers can compute a Nash equilibrium for a
        game with an arbitrary number of players.  These are numerical
        algorithms and so, in general, return floating point approximations
        of (and not necessarily all of) the equilibria.

        INPUT:

        - ``algorithm`` -- string; the Gambit solver to use (matched
          case-insensitively), one of

          * ``'lcp'`` -- the Linear Complementarity solver (two player games)
          * ``'lp'`` -- the Linear Programming solver (two player constant
            sum games)
          * ``'enummixed'`` -- enumeration of extreme points of convex sets of
            all Nash equilibria (two player games)
          * ``'gnm'`` -- the global Newton method (any number of players)
          * ``'enumpure'`` -- enumeration of the pure strategy equilibria
            (any number of players)
          * ``'enumpoly'`` -- enumeration of equilibria by solving systems of
            polynomial equations (any number of players)
          * ``'liap'`` -- minimisation of the Lyapunov function starting from
            the centroid (any number of players)
          * ``'simpdiv'`` -- simplicial subdivision starting from the centroid
            (any number of players)
          * ``'ipa'`` -- iterated polymatrix approximation (any number of
            players)
          * ``'logit'`` -- tracing of the logit quantal response equilibrium
            correspondence (any number of players)

        - ``maximization`` -- boolean (default: ``True``); whether the
          players maximize (``True``) or minimize (``False``) their utility

        - ``rational`` -- boolean (default: ``True``); whether to answer with
          rational probabilities rather than floating point ones.  ``'lcp'``,
          ``'lp'`` and ``'enummixed'`` are asked to use exact arithmetic
          throughout, while ``'gnm'``, ``'ipa'``, ``'logit'``, ``'liap'`` and
          ``'enumpoly'`` have no exact mode and instead have their floating
          point answer rounded by :func:`_rationalize_gambit_profile`, which
          accepts the result only if it is exactly an equilibrium.
          ``'enumpure'`` and ``'simpdiv'`` always work exactly and ignore the
          flag.  Exact arithmetic does not by itself make an answer exact --
          ``'simpdiv'`` refines a grid, so its rational output approximates an
          equilibrium.

        - ``tolerance`` -- a positive number (default: ``1e-4``); how far a
          probability may be moved to round it to a rational, passed to
          :func:`_rationalize_gambit_profile`.  Has no effect when ``rational``
          is ``False``, nor on a game whose payoffs are inexact -- such a game
          keeps its inexact equilibria, as decided by
          :func:`_gambit_payoffs_are_exact`.

        - ``phc_path`` -- (optional) a path (a string or
          :class:`~pathlib.Path`) to the PHCpack ``phc`` executable. Only
          supported by the ``'enumpoly'`` algorithm, for which it makes the
          underlying systems of polynomial equations be solved with PHCpack
          (using enumeration on the strategic game). Passing it for any other
          algorithm raises a :class:`ValueError`.

        OUTPUT: a sorted list of Nash equilibria, each a list with one tuple of
        probabilities per player, rational where they were found exactly and
        floats where they were not (see :meth:`_extract_gambit_equilibria`)

        EXAMPLES:

        The ``'lcp'`` solver on a two player game::

            sage: a = matrix([[1, 0], [1, 4]])
            sage: b = matrix([[2, 3], [2, 4]])
            sage: c = NormalFormGame([a, b])
            sage: c._use_gambit_solver('lcp')  # optional - pygambit
            [[(0, 1), (0, 1)]]

        The GNM solver on a two player game::

            sage: A = matrix([[1, 2], [3, 4]])
            sage: B = matrix([[3, 3], [1, 4]])
            sage: C = NormalFormGame([A, B])
            sage: C._use_gambit_solver('gnm')  # optional - pygambit
            [[(0, 1), (0, 1)]]

        GNM (like every solver other than ``'lcp'``, ``'lp'``, and ``'enummixed'``) can also
        solve games with more than two players.  Here is a three player game::

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
            sage: threegame._use_gambit_solver('gnm')  # optional - pygambit
            [[(0, 1), (1, 0), (1, 0)]]

        The LP solver on a constant sum game::

            sage: A = matrix([[2, 1], [1, 2.5]])
            sage: g = NormalFormGame([A])
            sage: g._use_gambit_solver('lp')  # optional - pygambit
            [[(3/5, 2/5), (3/5, 2/5)]]

        Players can also be set to minimize their utility, here using the
        Prisoner's Dilemma where players minimize time spent in prison::

            sage: A = matrix([[2, 5], [0, 4]])
            sage: B = matrix([[2, 0], [5, 4]])
            sage: prisoners_dilemma = NormalFormGame([A, B])
            sage: prisoners_dilemma._use_gambit_solver('gnm', maximization=False)  # optional - pygambit
            [[(0, 1), (0, 1)]]

        The remaining solvers are used in the same way::

            sage: c._use_gambit_solver('enumpure')  # optional - pygambit
            [[(0, 1), (0, 1)]]
            sage: c._use_gambit_solver('enumpoly')  # optional - pygambit
            [[(0, 1), (0, 1)]]
            sage: c._use_gambit_solver('simpdiv')  # optional - pygambit
            [[(0, 1), (0, 1)]]

        When the PHCpack ``phc`` executable is available, ``'enumpoly'`` can
        solve the underlying systems of polynomial equations with PHCpack
        instead of gambit's built-in solver by passing ``phc_path``.  PHCpack
        is a numerical homotopy continuation solver, and it finds the same
        equilibria as the built-in solver::

            sage: from shutil import which                          # optional - phc pygambit
            sage: c._use_gambit_solver('enumpoly', phc_path=which('phc'))  # optional - phc pygambit
            [[(0, 1), (0, 1)]]

        The ``'ipa'``, ``'liap'`` and ``'logit'`` solvers are iterative and
        compute in floating point, but their answers are rounded back to the
        exact equilibria they approximate (see ``tolerance`` above)::

            sage: # optional - pygambit
            sage: c._use_gambit_solver('ipa')
            [[(0, 1), (0, 1)]]
            sage: c._use_gambit_solver('liap')
            [[(0, 1), (0, 1)]]
            sage: c._use_gambit_solver('logit')
            [[(0, 1), (0, 1)]]

        The following examples cross-check the solvers against the equilibria
        recorded in Gambit's own test suite (``gambit/tests/test_nash.py``).
        The :math:`2\times 2` zero-sum game whose payoff matrices are
        :math:`I` and :math:`-I` is matching-pennies-like: it has the single
        equilibrium in which both players mix uniformly.  The ``'lcp'``,
        ``'enumpoly'`` and ``'simpdiv'`` solvers all recover it::

            sage: A = matrix.identity(2)
            sage: zero_sum = NormalFormGame([A, -A])
            sage: zero_sum._use_gambit_solver('lcp')  # optional - pygambit
            [[(1/2, 1/2), (1/2, 1/2)]]
            sage: zero_sum._use_gambit_solver('lcp', rational=False)  # abs tol 1e-9 # optional - pygambit
            [[(0.5, 0.5), (0.5, 0.5)]]
            sage: zero_sum._use_gambit_solver('enumpoly')  # optional - pygambit
            [[(1/2, 1/2), (1/2, 1/2)]]

        The game is constant sum, so the ``'lp'`` solver applies to its
        single-matrix form and finds the same equilibrium, while
        ``'enumpure'`` correctly reports that there is no equilibrium in pure
        strategies::

            sage: NormalFormGame([A])._use_gambit_solver('lp')  # optional - pygambit
            [[(1/2, 1/2), (1/2, 1/2)]]
            sage: zero_sum._use_gambit_solver('enumpure')  # optional - pygambit
            []

        The :math:`3\times 3` coordination game (both payoff matrices the
        identity) has :math:`2^3 - 1 = 7` equilibria, one uniform mix over
        each nonempty set of matching strategies.  ``'enumpoly'`` finds all
        seven, ``'enumpure'`` finds the three pure ones and ``'simpdiv'``
        returns the fully mixed one::

            sage: I3 = matrix.identity(3)
            sage: coordination = NormalFormGame([I3, I3])
            sage: coordination._use_gambit_solver('enumpure')  # optional - pygambit
            [[(0, 0, 1), (0, 0, 1)], [(0, 1, 0), (0, 1, 0)], [(1, 0, 0), (1, 0, 0)]]
            sage: coordination._use_gambit_solver('enumpoly')  # optional - pygambit
            [[(0, 0, 1), (0, 0, 1)],
             [(0, 1/2, 1/2), (0, 1/2, 1/2)],
             [(0, 1, 0), (0, 1, 0)],
             [(1/3, 1/3, 1/3), (1/3, 1/3, 1/3)],
             [(1/2, 0, 1/2), (1/2, 0, 1/2)],
             [(1/2, 1/2, 0), (1/2, 1/2, 0)],
             [(1, 0, 0), (1, 0, 0)]]
            sage: coordination._use_gambit_solver('simpdiv')  # optional - pygambit
            [[(1/3, 1/3, 1/3), (1/3, 1/3, 1/3)]]

        A 3 player game, coming from a local max cut instance, with two pure equilibria where
        either player is on one side of the cut and one mixed equilibrium where all players mix uniformly::

            sage: import numpy as np
            sage: A = np.array([[[0, -1], [2, 1]], [[1, 2],[-1, 0]]])
            sage: B = np.array([[[0, 2], [4, 2]], [[2, 4], [2, 0]]])
            sage: C = np.array([[[0, 1], [2, -1]], [[-1, 2], [1, 0]]])
            sage: max_cut_game = NormalFormGame([A, B, C])
            sage: max_cut_game._use_gambit_solver('enumpoly')  # optional - pygambit
            [[(0, 1), (1, 0), (0, 1)],
             [(1/2, 1/2), (1/2, 1/2), (1/2, 1/2)],
             [(1, 0), (0, 1), (1, 0)]]
            sage: max_cut_game._use_gambit_solver('enumpure')  # optional - pygambit
            [[(0, 1), (1, 0), (0, 1)], [(1, 0), (0, 1), (1, 0)]]

        Finally, a :math:`6\times 6` game with long Lemke-Howson paths and a
        unique equilibrium; the ``'gnm'``, ``'ipa'`` and ``'lcp'`` solvers all
        recover it::

            sage: A = matrix([[-180, 72, -333, 297, -153, 270],
            ....:             [-30, 17, -33, 42, -3, 20],
            ....:             [-81, 36, -126, 126, -36, 90],
            ....:             [90, -36, 126, -126, 36, -81],
            ....:             [20, -3, 42, -33, 17, -30],
            ....:             [270, -153, 297, -333, 72, -180]])
            sage: B = matrix([[72, 36, 17, -3, -36, -153],
            ....:             [-180, -81, -30, 20, 90, 270],
            ....:             [297, 126, 42, -33, -126, -333],
            ....:             [-333, -126, -33, 42, 126, 297],
            ....:             [270, 90, 20, -30, -81, -180],
            ....:             [-153, -36, -3, 17, 36, 72]])
            sage: long_lh = NormalFormGame([A, B])
            sage: long_lh._use_gambit_solver('gnm')  # optional - pygambit
            [[(1/30, 1/6, 3/10, 3/10, 1/6, 1/30), (1/6, 1/30, 3/10, 3/10, 1/30, 1/6)]]
            sage: long_lh._use_gambit_solver('ipa')  # optional - pygambit
            [[(1/30, 1/6, 3/10, 3/10, 1/6, 1/30), (1/6, 1/30, 3/10, 3/10, 1/30, 1/6)]]
            sage: long_lh._use_gambit_solver('lcp')  # optional - pygambit
            [[(1/30, 1/6, 3/10, 3/10, 1/6, 1/30), (1/6, 1/30, 3/10, 3/10, 1/30, 1/6)]]

        An unknown algorithm raises an error::

            sage: c._use_gambit_solver('invalid')
            Traceback (most recent call last):
            ...
            ValueError: unknown gambit algorithm 'invalid'; ...

        ``phc_path`` may only be given for the ``'enumpoly'`` algorithm::

            sage: c._use_gambit_solver('gnm', phc_path='/usr/bin/phc')
            Traceback (most recent call last):
            ...
            ValueError: 'phc_path' is only supported by the 'enumpoly' algorithm; got algorithm 'gnm'

        There is nothing to round to within a tolerance of zero::

            sage: c._use_gambit_solver('gnm', tolerance=0)
            Traceback (most recent call last):
            ...
            ValueError: 'tolerance' must be positive; got 0

        A game whose payoffs are inexact keeps its inexact equilibria, whatever
        ``rational`` is set to; there is nothing exact about them to recover::

            sage: # optional - pygambit
            sage: battle = NormalFormGame([matrix([[3, 0], [0, 2]]),
            ....:                          matrix([[2, 0], [0, 3]])])
            sage: battle._use_gambit_solver('liap')
            [[(3/5, 2/5), (2/5, 3/5)]]
            sage: battle._use_gambit_solver('liap', rational=False)  # abs tol 1e-6
            [[(0.6, 0.4), (0.4, 0.6)]]
            sage: inexact = NormalFormGame([matrix(RDF, [[3, 0], [0, 2]]),
            ....:                           matrix(RDF, [[2, 0], [0, 3]])])
            sage: inexact._use_gambit_solver('liap')  # abs tol 1e-6
            [[(0.6, 0.4), (0.4, 0.6)]]
        """

        if phc_path is not None and algorithm != "enumpoly":
            raise ValueError(
                "'phc_path' is only supported by the 'enumpoly' algorithm; "
                f"got algorithm {algorithm!r}"
            )

        if tolerance <= 0:
            raise ValueError(f"'tolerance' must be positive; got {tolerance!r}")

        pygambit().require()

        g = self._gambit_(maximization=maximization)

        def enumpoly_solve():
            # When ``phc_path`` is given the polynomial systems are solved with
            # PHCpack, which only supports enumeration on the strategic game,
            # hence ``use_strategic``.
            if phc_path is None:
                return gambit_nash.enumpoly_solve(g).equilibria
            # phc is an interactive program, driven here in batch mode, so it
            # writes its prompts, and the Ada exception raised when it reads
            # the end of its input, to the terminal.  That happens on
            # successful runs too, so it carries no information; gambit runs
            # phc without capturing either stream, hence the redirection at
            # file descriptor level, which is what a subprocess inherits.
            with open(os.devnull, 'w') as devnull:
                with redirection(1, devnull, close=False), \
                     redirection(2, devnull, close=False):
                    equilibria = gambit_nash.enumpoly_solve(
                        g, use_strategic=True,
                        phcpack_path=phc_path).equilibria
            # When PHCpack returns a solution that does not name every
            # variable of a support, gambit fills the missing probabilities in
            # with zeros; the all zero profile that results has zero regret
            # and so passes gambit's own check for being an equilibrium.  Such
            # a profile is not a mixed strategy profile at all, so drop
            # anything whose probabilities do not sum to one for every player.
            return [eq for eq in equilibria
                    if all(abs(sum(eq[s] for s in player.strategies) - 1) < 1e-6
                           for player in eq.game.players)]

        # Each solver has its own calling convention: ``lcp``/``lp`` take a
        # ``rational`` flag while ``liap``/``simpdiv`` need a starting mixed
        # strategy profile rather than the game itself.
        solvers = {
            'lcp': lambda: gambit_nash.lcp_solve(g, rational=rational).equilibria,
            'lp': lambda: gambit_nash.lp_solve(g, rational=rational).equilibria,
            'gnm': lambda: gambit_nash.gnm_solve(g).equilibria,
            'enumpure': lambda: gambit_nash.enumpure_solve(g).equilibria,
            'enummixed': lambda: gambit_nash.enummixed_solve(g, rational=rational).equilibria,
            'enumpoly': enumpoly_solve,
            'ipa': lambda: gambit_nash.ipa_solve(g).equilibria,
            'logit': lambda: gambit_nash.logit_solve(g).equilibria,
            'liap': lambda: gambit_nash.liap_solve(g.mixed_strategy_profile(rational=False)).equilibria,
            'simpdiv': lambda: gambit_nash.simpdiv_solve(g.mixed_strategy_profile(rational=True)).equilibria,
        }

        if algorithm not in solvers:
            raise ValueError(
                f"unknown gambit algorithm {algorithm!r}; "
                "supported values are 'lcp', 'lp', and 'enummixed' (2-player games only) "
                "and 'gnm', 'enumpure', 'enumpoly', 'liap', 'simpdiv', "
                "'ipa', 'logit' (any number of players)"
            )

        equilibria = solvers[algorithm]()
        if rational and _gambit_payoffs_are_exact(g):
            # The numerical solvers answer in floating point whatever the game
            # is made of, so ask for the exact equilibrium their answer renders;
            # ``None`` comes back when there is none, and the approximation is
            # then all there is to return.
            equilibria = [_rationalize_gambit_profile(eq, tolerance) or eq
                          for eq in equilibria]
        nasheq = self._extract_gambit_equilibria(equilibria)
        return sorted(nasheq)

    def _solve_LP(self, solver='glpk', maximization=True, rational=True):
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
            sage: g._solve_LP('gambit')  # optional - pygambit
            [[(1/2, 1/2), (1/2, 1/2)]]
            sage: g._solve_LP('Coin')  # optional - sage_numerical_backends_coin
            [[(0.5, 0.5), (0.5, 0.5)]]
            sage: g._solve_LP('PPL')
            [[(1/2, 1/2), (1/2, 1/2)]]
            sage: A = matrix([[2, 1], [1, 3]])
            sage: g = NormalFormGame([A])
            sage: ne = g._solve_LP()
            sage: [[[round(el, 6) for el in v] for v in eq] for eq in ne]
            [[[0.666667, 0.333333], [0.666667, 0.333333]]]
            sage: ne = g._solve_LP('gambit')  # optional - pygambit
            sage: [[[round(el, 6) for el in v] for v in eq] for eq in ne]  # optional - pygambit
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
            ValueError: the 'lp' algorithm only works for two-player constant-sum games, ...
        """
        if not self.is_constant_sum():
            raise ValueError(
                "the 'lp' algorithm only works for two-player constant-sum "
                "games, but this game is not constant-sum; use "
                "algorithm='enumeration', 'lrs', or 'lcp' instead"
            )
        if solver == 'gambit':
            return self._use_gambit_solver('lp', maximization,
                                           rational=rational)

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

        Testing against an error in ``_is_NE``.  Note that 1 equilibrium is
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

        potential_support_pairs = (pair for pair in product(*potential_supports) if len(pair[0]) == len(pair[1]))

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

        # if p1_payoffs.index(max(p1_payoffs)) not in p1_support:
        if not any(i in p1_support for i, x in enumerate(p1_payoffs)
                   if x == max(p1_payoffs)):
            return False
        return any(i in p2_support for i, x in enumerate(p2_payoffs)
                   if x == max(p2_payoffs))

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

        .. NOTE::

            The former legacy format has been removed in :issue:`39464`.
        """
        from sage.geometry.polyhedron.misc import _to_space_separated_string
        m = self.players[0].num_strategies
        n = self.players[1].num_strategies
        s = f'{m} {n}\n\n'
        s += '\n'.join(_to_space_separated_string(r) for r in m1.rows())
        s += '\n\n'
        s += '\n'.join(_to_space_separated_string(r) for r in m2.rows())
        s += '\n'
        return s

    def is_degenerate(self, certificate=False) -> bool:
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
        <https://www.samkass.com/theories/RPSSL.html>`_ is degenerate because
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
            sage: d_game.obtain_nash(algorithm='lcp')  # optional - pygambit
            [[(17/29, 0, 0, 12/29), (0, 0, 42/73, 31/73)]]
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
                    return True
            elif len(pair[1]) < len(pair[0]):
                strat = self._solve_indifference(pair[1], pair[0], M1.transpose())
                if strat and len(self.best_responses(strat, player=0)) > len(pair[1]):
                    if certificate:
                        return True, (strat, 1)
                    return True

        if certificate:
            return False, ()
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
                return True

        for j, col in enumerate(M1.columns()):
            if list(col).count(max(col)) > 1:
                if certificate:
                    strat = [0 for k in range(M1.ncols())]
                    strat[j] = 1
                    return True, (tuple(strat), 1)
                return True
        return False


class _Player:
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
