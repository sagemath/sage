"""
Normal Form games with N players.

This module implements a class for normal form games (strategic form games)
[NN2007]_. At present 3 algorithms are implemented to compute equilibria
of these games (``'lrs'`` - interfaced with the 'lrs' library, ``'LCP'`` interfaced
with the 'gambit' library and support enumeration built in Sage). The architecture
for the class is based on the gambit architecture to ensure an easy transition
between gambit and Sage.  At present the algorithms for the computation of equilibria
only solve 2 player games.

A very simple and well known example of normal form game is referred
to as the 'Battle of the Sexes' in which two players Amy and Bob
are modeled.  Amy prefers to play video games and Bob prefers to
watch a movie.  They both however want to spend their evening together.
This can be modeled using the following two matrices:

.. MATH::

    A = \\begin{pmatrix}
        3&1\\\\
        0&2\\\\
        \end{pmatrix}


    B = \\begin{pmatrix}
        2&1\\\\
        0&3\\\\
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
    Normal Form Game with the following utilities: {(0, 1): [1, 1], (1, 0): [0, 0], (0, 0): [3, 2], (1, 1): [2, 3]}

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
    ....:     print "Utility for {}: ".format(ne)
    ....:     print vector(ne[0]) * A * vector(ne[1]), vector(ne[0]) * B * vector(ne[1])
    Utility for [(0, 1), (0, 1)]:
    2 3
    Utility for [(3/4, 1/4), (1/4, 3/4)]:
    3/2 3/2
    Utility for [(1, 0), (1, 0)]:
    3 2

Allowing players to play mixed strategies ensures that there will always
be a Nash Equilibrium for a normal form game. This result is called Nash's
Theorem ([N1950]_).

Let us consider the game called 'matching pennies' where two players each
present a coin with either HEADS or TAILS showing. If the coins show the
same side then player 1 wins, otherwise player 2 wins:


.. MATH::

    A = \\begin{pmatrix}
        1&-1\\\\
        -1&1\\\\
        \end{pmatrix}


    B = \\begin{pmatrix}
        -1&1\\\\
        1&-1\\\\
        \end{pmatrix}

It should be relatively straightforward to observe, that there is no
situation, where both players always do the same thing, and have no
incentive to deviate.

We can plot the utility of player 1 when player 2 is playing a mixed
strategy `\sigma_2=(y,1-y)` (so that the utility to player 1 for
playing strategy number `i` is given by the matrix/vector multiplication
`(Ay)_i`, ie element in position `i` of the matrix/vector multiplication
`Ay`) ::

    sage: y = var('y')
    sage: A = matrix([[1, -1], [-1, 1]])
    sage: p = plot((A * vector([y, 1 - y]))[0], y, 0, 1, color='blue', legend_label='$u_1(r_1, (y, 1-y))$', axes_labels=['$y$', ''])
    sage: p += plot((A * vector([y, 1 - y]))[1], y, 0, 1, color='red', legend_label='$u_1(r_2, (y, 1-y))$'); p
    Graphics object consisting of 2 graphics primitives

We see that the only point at which player 1 is indifferent amongst
the available strategies is when `y=1/2`.

If we compute the Nash equilibria we see that this corresponds to a point
at which both players are indifferent::

    sage: A = matrix([[1, -1], [-1, 1]])
    sage: B = matrix([[-1, 1], [1, -1]])
    sage: matching_pennies = NormalFormGame([A, B])
    sage: matching_pennies.obtain_nash(algorithm='enumeration')
    [[(1/2, 1/2), (1/2, 1/2)]]

The utilities to both players at this Nash equilibrium
is easily computed::

    sage: [vector([1/2, 1/2]) * M * vector([1/2, 1/2]) for M in matching_pennies.payoff_matrices()]
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

When obtaining Nash equilibrium there are 3 algorithms currently available:

* ``'lrs'``: Reverse search vertex enumeration for 2 player games. This
  algorithm uses the optional 'lrs' package. To install it type ``sage -i
  lrs`` at the command line. For more information see [A2000]_.

* ``'LCP'``: Linear complementarity program algorithm for 2 player games.
  This algorithm uses the open source game theory package:
  `Gambit <http://gambit.sourceforge.net/>`_ [MMAT2014]_. At present this is
  the only gambit algorithm available in sage but further development will
  hope to implement more algorithms
  (in particular for games with more than 2 players). To install it
  type ``sage -i gambit`` at the command line.

* ``'enumeration'``: Support enumeration for 2 player games. This
  algorithm is hard coded in Sage and checks through all potential
  supports of a strategy. Supports of a given size with a conditionally
  dominated strategy are ignored. Note: this is not the preferred
  algorithm. The algorithm implemented is a combination of a basic
  algorithm described in [NN2007]_ and a pruning component described
  in [SLB2008]_.

Below we show how the three algorithms are called::

    sage: matching_pennies.obtain_nash(algorithm='lrs')  # optional - lrs
    [[(1/2, 1/2), (1/2, 1/2)]]
    sage: matching_pennies.obtain_nash(algorithm='LCP')  # optional - gambit
    [[(0.5, 0.5), (0.5, 0.5)]]
    sage: matching_pennies.obtain_nash(algorithm='enumeration')
    [[(1/2, 1/2), (1/2, 1/2)]]

Note that if no algorithm argument is passed then the default will be
selected according to the following order (if the corresponding package is
installed):

    1. ``'lrs'`` (requires 'lrs')
    2. ``'enumeration'``

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
    Normal Form Game with the following utilities: {(0, 1): [2, 3], (1, 0): [3, 1], (0, 0): [1, 3], (1, 1): [4, 4]}

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
    Normal Form Game with the following utilities: {(0, 1): [2, 3], (0, 0): [1, 3], (2, 1): [False, False], (2, 0): [False, False], (1, 0): [3, 1], (1, 1): [4, 4]}

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
    Normal Form Game with the following utilities: {(0, 1, 1): [3, 5, 8], (1, 1, 0): [8, 4, 6], (1, 0, 0): [9, 7, 9], (0, 0, 1): [1, 5, 9], (1, 0, 1): [3, 2, 3], (0, 0, 0): [3, 1, 4], (0, 1, 0): [2, 6, 5], (1, 1, 1): [2, 6, 4]}

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
    Normal Form Game with the following utilities: {(0, 1, 1): [0, 2, 4], (1, 1, 0): [0, 2, 4], (1, 0, 0): [0, 1, 2], (0, 0, 1): [0, 1, 2], (1, 0, 1): [0, 2, 4], (0, 0, 0): [0, 0, 0], (0, 1, 0): [0, 1, 2], (1, 1, 1): [0, 3, 6]}

At present no algorithm has been implemented in Sage for games with
more than 2 players::

    sage: threegame.obtain_nash()
    Traceback (most recent call last):
    ...
    NotImplementedError: Nash equilibrium for games with more than 2 players have not been implemented yet. Please see the gambit website (http://gambit.sourceforge.net/) that has a variety of available algorithms

There are however a variety of such algorithms available in gambit,
further compatibility between Sage and gambit is actively being developed:
https://github.com/tturocy/gambit/tree/sage_integration.

Note that the Gambit implementation of ``LCP`` can only handle integer
payoffs. If a non integer payoff is used an error will be raised::

    sage: A = matrix([[2, 1], [1, 2.5]])
    sage: B = matrix([[-1, 3], [2, 1]])
    sage: g = NormalFormGame([A, B])
    sage: g.obtain_nash(algorithm='LCP')  # optional - gambit
    Traceback (most recent call last):
    ...
    ValueError: The Gambit implementation of LCP only allows for integer valued payoffs. Please scale your payoff matrices.

Other algorithms can handle these payoffs::

    sage: g.obtain_nash(algorithm='enumeration')
    [[(1/5, 4/5), (3/5, 2/5)]]
    sage: g.obtain_nash(algorithm='lrs') # optional - lrs
    [[(1/5, 4/5), (3/5, 2/5)]]

It can be shown that linear scaling of the payoff matrices conserves the
equilibrium values::

    sage: A = 2 * A
    sage: g = NormalFormGame([A, B])
    sage: g.obtain_nash(algorithm='LCP')  # optional - gambit
    [[(0.2, 0.8), (0.6, 0.4)]]

It is also possible to generate a Normal Form Game from a gambit Game::

    sage: from gambit import Game  # optional - gambit
    sage: gambitgame= Game.new_table([2, 2])  # optional - gambit
    sage: gambitgame[int(0), int(0)][int(0)] = int(8)  # optional - gambit
    sage: gambitgame[int(0), int(0)][int(1)] = int(8)  # optional - gambit
    sage: gambitgame[int(0), int(1)][int(0)] = int(2)  # optional - gambit
    sage: gambitgame[int(0), int(1)][int(1)] = int(10)  # optional - gambit
    sage: gambitgame[int(1), int(0)][int(0)] = int(10)  # optional - gambit
    sage: gambitgame[int(1), int(0)][int(1)] = int(2)  # optional - gambit
    sage: gambitgame[int(1), int(1)][int(0)] = int(5)  # optional - gambit
    sage: gambitgame[int(1), int(1)][int(1)] = int(5)  # optional - gambit
    sage: g = NormalFormGame(gambitgame)  # optional - gambit
    sage: g  # optional - gambit
    Normal Form Game with the following utilities: {(0, 1): [2.0, 10.0], (1, 0): [10.0, 2.0], (0, 0): [8.0, 8.0], (1, 1): [5.0, 5.0]}

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
the amount of their value at no less than 2 and no larger than 100. He
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
    sage: A = matrix([[min(i,j) + 2 * sign(j-i)  for j in range(2, K+1)]  for i in range(2, K+1)])
    sage: B = matrix([[min(i,j) + 2 * sign(i-j)  for j in range(2, K+1)]  for i in range(2, K+1)])
    sage: g = NormalFormGame([A, B])
    sage: g.obtain_nash(algorithm='lrs') # optional - lrs
    [[(1, 0, 0, 0, 0, 0, 0, 0, 0), (1, 0, 0, 0, 0, 0, 0, 0, 0)]]
    sage: g.obtain_nash(algorithm='LCP') # optional - gambit
    [[(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0), (1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)]]

The output is a pair of vectors (as before) showing the Nash equilibrium.
In particular it here shows that out of the 10 possible strategies both
players should choose the first. Recall that the above considers a reduced
version of the game where individuals can claim integer values between 2
and 10.  The equilibrium strategy is thus for both players to state that
the value of their suitcase is 2.

Note that degenerate games can cause problems for most algorithms.
The following example in fact has an infinite quantity of equilibria which
is evidenced by the various algorithms returning different solutions::

    sage: A = matrix([[3,3],[2,5],[0,6]])
    sage: B = matrix([[3,3],[2,6],[3,1]])
    sage: degenerate_game = NormalFormGame([A,B])
    sage: degenerate_game.obtain_nash(algorithm='lrs') # optional - lrs
    [[(0, 1/3, 2/3), (1/3, 2/3)], [(1, 0, 0), (2/3, 1/3)], [(1, 0, 0), (1, 0)]]
    sage: degenerate_game.obtain_nash(algorithm='LCP') # optional - gambit
    [[(0.0, 0.3333333333, 0.6666666667), (0.3333333333, 0.6666666667)],
     [(1.0, -0.0, 0.0), (0.6666666667, 0.3333333333)],
     [(1.0, 0.0, 0.0), (1.0, 0.0)]]
    sage: degenerate_game.obtain_nash(algorithm='enumeration')
    [[(0, 1/3, 2/3), (1/3, 2/3)], [(1, 0, 0), (1, 0)]]

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
    [[(0, 0, 1), (0, 0, 1)], [(0, 0, 1), (0, 1, 0)], [(0, 0, 1), (1, 0, 0)], [(0, 1, 0), (0, 0, 1)], [(0, 1, 0), (0, 1, 0)], [(0, 1, 0), (1, 0, 0)], [(1, 0, 0), (0, 0, 1)], [(1, 0, 0), (0, 1, 0)], [(1, 0, 0), (1, 0, 0)]]

A good description of degenerate games can be found in [NN2007]_.

REFERENCES:

.. [N1950] John Nash.
   *Equilibrium points in n-person games.*
   Proceedings of the National Academy of Sciences 36.1 (1950): 48-49.

.. [NN2007] Nisan, Noam, et al., eds.
   *Algorithmic game theory.*
   Cambridge University Press, 2007.

.. [A2000] Avis, David.
   *A revised implementation of the reverse search vertex enumeration algorithm.*
   Polytopes-combinatorics and computation
   Birkhauser Basel, 2000.

.. [MMAT2014] McKelvey, Richard D., McLennan, Andrew M., and Turocy, Theodore L.
   *Gambit: Software Tools for Game Theory, Version 13.1.2.*
   http://www.gambit-project.org (2014).

.. [SLB2008] Shoham, Yoav, and Kevin Leyton-Brown.
   *Multiagent systems: Algorithmic, game-theoretic, and logical foundations.*
   Cambridge University Press, 2008.

AUTHOR:

- James Campbell and Vince Knight (06-2014): Original version

"""

#*****************************************************************************
#       Copyright (C) 2014 James Campbell james.campbell@tanti.org.uk
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from collections import MutableMapping
from itertools import product
from parser import Parser
from sage.combinat.cartesian_product import CartesianProduct
from sage.misc.latex import latex
from sage.misc.misc import powerset
from sage.rings.all import QQ
from sage.structure.sage_object import SageObject
from sage.matrix.constructor import matrix
from sage.matrix.constructor import vector
from sage.misc.package import is_package_installed
from sage.misc.temporary_file import tmp_filename

try:
    from gambit import Game
except ImportError:
    Game = None

class NormalFormGame(SageObject, MutableMapping):
    r"""
    An object representing a Normal Form Game. Primarily used to compute the
    Nash Equilibria.

    INPUT:

    - ``generator`` - Can be a list of 2 matrices, a single matrix or left
                      blank.

    """

    def __init__(self, generator=None):
        r"""
        Initializes a Normal Form game and checks the inputs.

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
            NotImplementedError: Nash equilibrium for games with more than 2 players have not been implemented yet. Please see the gambit website (http://gambit.sourceforge.net/) that has a variety of available algorithms

        Can initialise a game from a gambit game object::

            sage: from gambit import Game  # optional - gambit
            sage: gambitgame= Game.new_table([2, 2])  # optional - gambit
            sage: gambitgame[int(0), int(0)][int(0)] = int(5)  # optional - gambit
            sage: gambitgame[int(0), int(0)][int(1)] = int(8)  # optional - gambit
            sage: gambitgame[int(0), int(1)][int(0)] = int(2)  # optional - gambit
            sage: gambitgame[int(0), int(1)][int(1)] = int(11)  # optional - gambit
            sage: gambitgame[int(1), int(0)][int(0)] = int(10)  # optional - gambit
            sage: gambitgame[int(1), int(0)][int(1)] = int(7)  # optional - gambit
            sage: gambitgame[int(1), int(1)][int(0)] = int(5)  # optional - gambit
            sage: gambitgame[int(1), int(1)][int(1)] = int(5)  # optional - gambit
            sage: g = NormalFormGame(gambitgame)  # optional - gambit
            sage: g  # optional - gambit
            Normal Form Game with the following utilities: {(0, 1): [2.0, 11.0], (1, 0): [10.0, 7.0], (0, 0): [5.0, 8.0], (1, 1): [5.0, 5.0]}

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
            TypeError: __init__() takes at most 2 arguments (3 given)

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
            Normal Form Game with the following utilities: {(0, 1): [5, 0], (1, 0): [0, 5], (0, 0): [2, 2], (1, 1): [4, 4]}
            sage: del(prisoners_dilemma[(0,1)])
            sage: prisoners_dilemma
            Normal Form Game with the following utilities: {(1, 0): [0, 5], (0, 0): [2, 2], (1, 1): [4, 4]}
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
            sage: for key in prisoners_dilemma:
            ....:     print "The strategy pair {} gives utilities {}".format(key, prisoners_dilemma[key])
            The strategy pair (0, 1) gives utilities [5, 0]
            The strategy pair (1, 0) gives utilities [0, 5]
            The strategy pair (0, 0) gives utilities [2, 2]
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
        Here we return the length of the
        game to be the length of the utilities::

            sage: A = matrix([[2, 5], [0, 4]])
            sage: B = matrix([[2, 0], [5, 4]])
            sage: prisoners_dilemma = NormalFormGame([A, B])
            sage: len(prisoners_dilemma)
            4
        """
        return len(self.utilities)

    def _repr_(self):
        r"""
        Returns the strategy_profiles of the game.

        EXAMPLES:

        Basic description of the game shown when calling the game instance::

            sage: p1 = matrix([[1, 2], [3, 4]])
            sage: p2 = matrix([[3, 3], [1, 4]])
            sage: g = NormalFormGame([p1, p2])
            sage: g
            Normal Form Game with the following utilities: {(0, 1): [2, 3], (1, 0): [3, 1], (0, 0): [1, 3], (1, 1): [4, 4]}
        """
        base_str = "Normal Form Game with the following utilities: {}"
        return base_str.format(self.utilities)

    def _latex_(self):
        r"""
        Returns the LaTeX code representing the ``NormalFormGame``.

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
            \text{\texttt{Normal{ }Form{ }Game{ }...[False,{ }False,{ }False]{\char`\}}}}
        """
        if len(self.players) == 2:
            M1, M2 = self.payoff_matrices()
            return "\left(%s, %s\\right)" % (M1._latex_(), M2._latex_())
        return latex(self.__str__())

    def _two_matrix_game(self, matrices):
        r"""
        Populates ``self.utilities`` with the values from 2 matrices.

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
        Creates a ``NormalFormGame`` object from a Gambit game.

        TESTS::

            sage: from gambit import Game  # optional - gambit
            sage: testgame = Game.new_table([2, 2])  # optional - gambit
            sage: testgame[int(0), int(0)][int(0)] = int(8)  # optional - gambit
            sage: testgame[int(0), int(0)][int(1)] = int(8)  # optional - gambit
            sage: testgame[int(0), int(1)][int(0)] = int(2)  # optional - gambit
            sage: testgame[int(0), int(1)][int(1)] = int(10)  # optional - gambit
            sage: testgame[int(1), int(0)][int(0)] = int(10)  # optional - gambit
            sage: testgame[int(1), int(0)][int(1)] = int(2)  # optional - gambit
            sage: testgame[int(1), int(1)][int(0)] = int(5)  # optional - gambit
            sage: testgame[int(1), int(1)][int(1)] = int(5)  # optional - gambit
            sage: g = NormalFormGame()  # optional - gambit
            sage: g._gambit_game(testgame)  # optional - gambit
            sage: g  # optional - gambit
            Normal Form Game with the following utilities: {(0, 1): [2.0, 10.0], (1, 0): [10.0, 2.0], (0, 0): [8.0, 8.0], (1, 1): [5.0, 5.0]}
        """
        self.players = []
        self.utilities = {}
        for player in game.players:
            num_strategies = len(player.strategies)
            self.add_player(num_strategies)
        for strategy_profile in self.utilities:
            utility_vector = [float(game[strategy_profile][i]) for i in range(len(self.players))]
            self.utilities[strategy_profile] = utility_vector

    def payoff_matrices(self):
        r"""
        Returns 2 matrices representing the payoffs for each player.

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
            sage: g.add_player(2)  # Adding first player with 2 strategies
            sage: g.add_player(2)  # Adding second player with 2 strategies
            sage: g.add_player(2)  # Adding third player with 2 strategies
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
        Adds a player to a NormalFormGame.

        INPUT:

        - ``num_strategies`` - the number of strategies the player should have.

        EXAMPLES::

            sage: g = NormalFormGame()
            sage: g.add_player(2)  # Adding first player with 2 strategies
            sage: g.add_player(1)  # Adding second player with 1 strategy
            sage: g.add_player(1)  # Adding third player with 1 strategy
            sage: g
            Normal Form Game with the following utilities: {(1, 0, 0): [False, False, False], (0, 0, 0): [False, False, False]}
        """
        self.players.append(_Player(num_strategies))
        self._generate_utilities(True)

    def _generate_utilities(self, replacement):
        r"""
        Creates all the required keys for ``self.utilities``.
        This is used when generating players and/or adding strategies.

        INPUT:

        - ``replacement`` -- Boolean value of whether previously created
                            profiles should be replaced or not.

        TESTS::

            sage: from sage.game_theory.normal_form_game import _Player
            sage: g = NormalFormGame()
            sage: g.players.append(_Player(2))
            sage: g.players.append(_Player(2))
            sage: g
            Normal Form Game with the following utilities: {}

            sage: g._generate_utilities(True)
            sage: g
            Normal Form Game with the following utilities: {(0, 1): [False, False], (1, 0): [False, False], (0, 0): [False, False], (1, 1): [False, False]}

            sage: g[(0,1)] = [2, 3]
            sage: g.add_strategy(1)
            sage: g._generate_utilities(False)
            sage: g
            Normal Form Game with the following utilities: {(0, 1): [2, 3], (1, 2): [False, False], (0, 0): [False, False], (0, 2): [False, False], (1, 0): [False, False], (1, 1): [False, False]}

            sage: g._generate_utilities(True)
            sage: g
            Normal Form Game with the following utilities: {(0, 1): [False, False], (1, 2): [False, False], (0, 0): [False, False], (1, 1): [False, False], (1, 0): [False, False], (0, 2): [False, False]}
        """
        strategy_sizes = [range(p.num_strategies) for p in self.players]
        if replacement is True:
            self.utilities = {}
        for profile in product(*strategy_sizes):
            if profile not in self.utilities.keys():
                self.utilities[profile] = [False]*len(self.players)

    def add_strategy(self, player):
        r"""
        Adds a strategy to a player, will not affect already completed
        strategy profiles.

        INPUT:

        - ``player`` - the index of the player.

        EXAMPLES:

        A simple example::

            sage: s = matrix([[1, 0], [-2, 3]])
            sage: t = matrix([[3, 2], [-1, 0]])
            sage: example = NormalFormGame([s, t])
            sage: example
            Normal Form Game with the following utilities: {(0, 1): [0, 2], (1, 0): [-2, -1], (0, 0): [1, 3], (1, 1): [3, 0]}
            sage: example.add_strategy(0)
            sage: example
            Normal Form Game with the following utilities: {(0, 1): [0, 2], (0, 0): [1, 3], (2, 1): [False, False], (2, 0): [False, False], (1, 0): [-2, -1], (1, 1): [3, 0]}

        """
        self.players[player].add_strategy()
        self._generate_utilities(False)

    def _is_complete(self):
        r"""
        Checks if ``utilities`` has been completed and returns a
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
        results = []
        for profile in self.utilities.values():
            results.append(all(type(i) is not bool for i in profile))
        return all(results)

    def obtain_nash(self, algorithm=False, maximization=True):
        r"""
        A function to return the Nash equilibrium for the game.
        Optional arguments can be used to specify the algorithm used.
        If no algorithm is passed then an attempt is made to use the most
        appropriate algorithm.

        INPUT:

        - ``algorithm`` - the following algorithms should be available through
                          this function:

          * ``'lrs'`` - This algorithm is only suited for 2 player games.
            See the lrs web site (http://cgm.cs.mcgill.ca/~avis/C/lrs.html).

          * ``'LCP'`` - This algorithm is only suited for 2 player games.
            See the gambit web site (http://gambit.sourceforge.net/). Note
            that the output differs from the other algorithms: floats are
            returned.

          * ``'enumeration'`` - This is a very inefficient
            algorithm (in essence a brute force approach).

            1. For each k in 1...min(size of strategy sets)
            2. For each I,J supports of size k
            3. Prune: check if supports are dominated
            4. Solve indifference conditions and check that have Nash Equilibrium.

            Solving the indifference conditions is done by building the
            corresponding linear system.  If  `\rho_1, \rho_2` are the
            supports player 1 and 2 respectively.  Then, indifference implies:

            .. MATH::

                u_1(s_1,\rho_2) = u_2(s_2, \rho_2)

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

        - ``maximization`` - Whether a player is trying to maximize their
                             utility or minimize it.

          * When set to ``True`` (default) it is assumed that players
            aim to maximise their utility.

          * When set to ``False`` it is assumed that players aim to
            minimise their utility.

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
            sage: g.obtain_nash(algorithm='lrs') # optional - lrs
            [[(0, 0, 0, 1), (0, 0, 1)]]
            sage: g.obtain_nash(algorithm='lrs', maximization=False) # optional - lrs
            [[(2/3, 1/12, 1/4, 0), (6333/8045, 247/8045, 293/1609)], [(3/4, 0, 1/4, 0), (0, 11/307, 296/307)], [(5/6, 1/6, 0, 0), (98/99, 1/99, 0)]]

        This particular game has 3 Nash equilibria::

            sage: A = matrix([[3,3],
            ....:             [2,5],
            ....:             [0,6]])
            sage: B = matrix([[3,2],
            ....:             [2,6],
            ....:             [3,1]])
            sage: g = NormalFormGame([A, B])
            sage: g.obtain_nash(algorithm='enumeration')
            [[(0, 1/3, 2/3), (1/3, 2/3)], [(4/5, 1/5, 0), (2/3, 1/3)], [(1, 0, 0), (1, 0)]]

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
            sage: g.obtain_nash(algorithm='lrs')  # optional - lrs
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
            sage: fivegame.obtain_nash(algorithm='lrs') # optional - lrs
            [[(1, 0, 0, 0, 0), (0, 1, 0, 0, 0)]]
            sage: fivegame.obtain_nash(algorithm='LCP') # optional - gambit
            [[(1.0, 0.0, 0.0, 0.0, 0.0), (0.0, 1.0, 0.0, 0.0, 0.0)]]

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

        Note that outputs for all algorithms are as lists of lists of
        tuples and the equilibria have been sorted so that all algorithms give
        a comparable output (although ``'LCP'`` returns floats)::

            sage: enumeration_eqs = g.obtain_nash(algorithm='enumeration')
            sage: [[type(s) for s in eq] for eq in enumeration_eqs]
            [[<type 'tuple'>, <type 'tuple'>], [<type 'tuple'>, <type 'tuple'>], [<type 'tuple'>, <type 'tuple'>]]
            sage: lrs_eqs = g.obtain_nash(algorithm='lrs')  # optional - lrs
            sage: [[type(s) for s in eq] for eq in lrs_eqs]  # optional - lrs
            [[<type 'tuple'>, <type 'tuple'>], [<type 'tuple'>, <type 'tuple'>], [<type 'tuple'>, <type 'tuple'>]]
            sage: LCP_eqs = g.obtain_nash(algorithm='LCP')  # optional - gambit
            sage: [[type(s) for s in eq] for eq in LCP_eqs]  # optional - gambit
            [[<type 'tuple'>, <type 'tuple'>], [<type 'tuple'>, <type 'tuple'>], [<type 'tuple'>, <type 'tuple'>]]
            sage: enumeration_eqs == sorted(enumeration_eqs)
            True
            sage: lrs_eqs == sorted(lrs_eqs)  # optional - lrs
            True
            sage: LCP_eqs == sorted(LCP_eqs)  # optional - gambit
            True
            sage: lrs_eqs == enumeration_eqs  # optional - lrs
            True
            sage: enumeration_eqs == LCP_eqs  # optional - gambit
            False
            sage: [[[round(float(p), 6) for p in str] for str in eq] for eq in enumeration_eqs] == [[[round(float(p), 6) for p in str] for str in eq] for eq in LCP_eqs]  # optional - gambit
            True

        """
        if len(self.players) > 2:
            raise NotImplementedError("Nash equilibrium for games with more "
                                      "than 2 players have not been "
                                      "implemented yet. Please see the gambit "
                                      "website (http://gambit.sourceforge.net/) that has a variety of "
                                      "available algorithms")

        if not self._is_complete():
            raise ValueError("utilities have not been populated")

        if not algorithm:
            if is_package_installed('lrs'):
                algorithm = "lrs"
            else:
                algorithm = "enumeration"

        if algorithm == "lrs":
            if not is_package_installed('lrs'):
                raise NotImplementedError("lrs is not installed")

            return self._solve_lrs(maximization)

        if algorithm == "LCP":
            if Game is None:
                raise NotImplementedError("gambit is not installed")
            for strategy_profile in self.utilities:
                payoffs = self.utilities[strategy_profile]
                if payoffs != [int(payoffs[0]), int(payoffs[1])]:
                    raise ValueError("""The Gambit implementation of LCP only
                                     allows for integer valued payoffs.
                                     Please scale your payoff matrices.""")
            return self._solve_LCP(maximization)

        if algorithm == "enumeration":
            return self._solve_enumeration(maximization)

    def _solve_lrs(self, maximization=True):
        r"""
        EXAMPLES:

        A simple game::

            sage: A = matrix([[1, 2], [3, 4]])
            sage: B = matrix([[3, 3], [1, 4]])
            sage: C = NormalFormGame([A, B])
            sage: C._solve_lrs() # optional - lrs
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
            sage: biggame._solve_lrs() # optional - lrs
            [[(0, 0, 0, 20/21, 1/21), (11/12, 0, 0, 1/12, 0)]]

        Another test::

            sage: p1 = matrix([[-7, -5, 5],
            ....:              [5, 5, 3],
            ....:              [1, -6, 1]])
            sage: p2 = matrix([[-9, 7, 9],
            ....:              [6, -2, -3],
            ....:              [-4, 6, -10]])
            sage: biggame = NormalFormGame([p1, p2])
            sage: biggame._solve_lrs() # optional - lrs
            [[(0, 1, 0), (1, 0, 0)], [(1/3, 2/3, 0), (0, 1/6, 5/6)], [(1/3, 2/3, 0), (1/7, 0, 6/7)], [(1, 0, 0), (0, 0, 1)]]
        """
        from subprocess import PIPE, Popen
        m1, m2 = self.payoff_matrices()
        if maximization is False:
            m1 = - m1
            m2 = - m2
        game1_str, game2_str = self._Hrepresentation(m1, m2)

        g1_name = tmp_filename()
        g2_name = tmp_filename()
        g1_file = file(g1_name, 'w')
        g2_file = file(g2_name, 'w')
        g1_file.write(game1_str)
        g1_file.close()
        g2_file.write(game2_str)
        g2_file.close()

        process = Popen(['nash', g1_name, g2_name], stdout=PIPE)
        lrs_output = [row for row in process.stdout]
        nasheq = Parser(lrs_output).format_lrs()
        return sorted(nasheq)

    def _solve_LCP(self, maximization):
        r"""
        Solves a NormalFormGame using Gambit's LCP algorithm.

        EXAMPLES::

            sage: a = matrix([[1, 0], [1, 4]])
            sage: b = matrix([[2, 3], [2, 4]])
            sage: c = NormalFormGame([a, b])
            sage: c._solve_LCP(maximization=True) # optional - gambit
            [[(0.0, 1.0), (0.0, 1.0)]]
        """
        from gambit.nash import ExternalLCPSolver
        strategy_sizes = [p.num_strategies for p in self.players]
        g = Game.new_table(strategy_sizes)
        scalar = 1
        if maximization is False:
            scalar *= -1
        for strategy_profile in self.utilities:
            g[strategy_profile][0] = int(scalar *
                                            self.utilities[strategy_profile][0])
            g[strategy_profile][1] = int(scalar *
                                            self.utilities[strategy_profile][1])
        output = ExternalLCPSolver().solve(g)
        nasheq = Parser(output).format_LCP(g)
        return sorted(nasheq)

    def _solve_enumeration(self, maximization=True):
        r"""
        Obtains the Nash equilibria using support enumeration.
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
            [[(0, 0, 0, 1), (1, 0, 0, 0)], [(2/7, 0, 0, 5/7), (5/11, 0, 6/11, 0)], [(1, 0, 0, 0), (0, 0, 1, 0)]]

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
        missing: `[(2/3, 1/3), (0, 1)]` however this equilibrium has
        supports of different sizes. This only occurs in degenerate games
        and is not supported in the `enumeration` algorithm::

            sage: N = NormalFormGame([matrix(2,[0,-1,-2,-1]),matrix(2,[1,0,0,2])])
            sage: N._solve_enumeration()
            [[(0, 1), (0, 1)], [(1, 0), (1, 0)]]

        In this instance the `lrs` algorithm is able to find all three equilibria::

            sage: N = NormalFormGame([matrix(2,[0,-1,-2,-1]),matrix(2,[1,0,0,2])])
            sage: N.obtain_nash(algorithm='lrs')  # optional - lrs
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

        potential_support_pairs = [pair for pair in CartesianProduct(*potential_supports) if len(pair[0]) == len(pair[1])]

        equilibria = []
        for pair in potential_support_pairs:
            # Check if any supports are dominated for row player
            if (self._row_cond_dominance(pair[0], pair[1], M1)
                # Check if any supports are dominated for col player
               and self._row_cond_dominance(pair[1], pair[0], M2.transpose())):
                    result = self._solve_indifference(pair[0], pair[1], M1, M2)
                    if result:
                        equilibria.append([tuple(result[0]), tuple(result[1])])
        return sorted(equilibria)

    def _row_cond_dominance(self, p1_sup, p2_sup, matrix):
        r"""
        Checks if any row strategies of a sub matrix defined
        by a given pair of supports are conditionally dominated.
        Returns False if a row is conditionally dominated.

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

    def _solve_indifference(self, p1_support, p2_support, M1, M2):
        r"""
        For a support pair obtains vector pair that ensures indifference
        amongst support strategies.

        This is done by building the corresponding linear system.
        If  `\rho_1, \rho_2` are the supports player 1 and 2 respectively.
        Then, indifference implies:

        .. MATH::

            u_1(s_1,\rho_2) = u_2(s_2, \rho_2)

        for all `s_1, s_2` in the support of `\rho_1`. This corresponds to:

        .. MATH::

            \sum_{j\in S(\rho_2)}A_{s_1,j}{\rho_2}_j = \sum_{j\in S(\rho_2)}A_{s_2,j}{\rho_2}_j

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

            \sum_{j\in S(\rho_1)}{\rho_2}_j = 1

        TESTS:

        Find the indifference vector for a support pair that has
        no dominated strategies::

            sage: A = matrix([[1, 1, 5], [2, 2, 0]])
            sage: g = NormalFormGame([A])
            sage: g._solve_indifference((0, 1), (0, 2), A, -A)
            [(1/3, 2/3), (5/6, 0, 1/6)]

        When a support pair has a dominated strategy there is no
        solution to the indifference equation::

            sage: g._solve_indifference((0, 1), (0, 1), A, -A)
            <BLANKLINE>

        Particular case of a game with 1 strategy for each for each player::

            sage: A = matrix([[10]])
            sage: g = NormalFormGame([A])
            sage: g._solve_indifference((0,), (0,), A, -A)
            [(1), (1)]
        """
        linearsystem1 = matrix(QQ, len(p2_support)+1, self.players[0].num_strategies)
        linearsystem2 = matrix(QQ, len(p1_support)+1, self.players[1].num_strategies)

        # Build linear system for player 1
        for p1_strategy in p1_support:
            # Checking particular case of supports of pure strategies
            if len(p2_support) == 1:
                for p2_strategy in range(self.players[1].num_strategies):
                    if M2[p1_strategy][p2_support[0]] < \
                            M2[p1_strategy][p2_strategy]:
                        return False
            else:
                for p2_strategy_pair in range(len(p2_support)):
                    # Coefficients of linear system that ensure indifference between two consecutive strategies of the support of p1
                    linearsystem1[p2_strategy_pair, p1_strategy] = \
                            M2[p1_strategy][p2_support[p2_strategy_pair]] -\
                              M2[p1_strategy][p2_support[p2_strategy_pair-1]]
            linearsystem1[-1, p1_strategy] = 1  # Coefficients of linear system to ensure that vector is probability

        # Build linear system for player 2
        for p2_strategy in p2_support:
            # Checking particular case of supports of pure strategies
            if len(p1_support) == 1:
                for p1_strategy in range(self.players[0].num_strategies):
                    if M1[p1_support[0]][p2_strategy] < \
                            M1[p1_strategy][p2_strategy]:
                        return False
            else:
                for p1_strategy_pair in range(len(p1_support)):
                    # Coefficients of linear system that ensure indifference between two consecutive strategies of the support of p1
                    linearsystem2[p1_strategy_pair, p2_strategy] = \
                            M1[p1_support[p1_strategy_pair]][p2_strategy] -\
                              M1[p1_support[p1_strategy_pair-1]][p2_strategy]
            linearsystem2[-1, p2_strategy] = 1  # Coefficients of linear system that ensure that vector is probability

        # Create rhs of linear systems
        linearsystemrhs1 = vector([0 for i in range(len(p2_support))] + [1])
        linearsystemrhs2 = vector([0 for i in range(len(p1_support))] + [1])

        # Solve both linear systems
        try:
            a = linearsystem1.solve_right(linearsystemrhs1)
            b = linearsystem2.solve_right(linearsystemrhs2)
        except ValueError:
            return None

        if self._is_NE(a, b, p1_support, p2_support, M1, M2):
            return [a, b]
        return None

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
        if not (all([a[i] > 0 for i in p1_support]) and
            all([b[j] > 0 for j in p2_support]) and
            all([a[i] == 0 for i in range(len(a)) if i not in p1_support]) and
            all([b[j] == 0 for j in range(len(b)) if j not in p2_support])):
            return False

        # Check that have pair of best responses

        p1_payoffs = [sum(v * row[i] for i, v in enumerate(b)) for row
                                                                  in M1.rows()]
        p2_payoffs = [sum(v * col[j] for j, v in enumerate(a)) for col
                                                               in M2.columns()]

        #if p1_payoffs.index(max(p1_payoffs)) not in p1_support:
        if not any(i in p1_support for i, x in enumerate(p1_payoffs) if x == max(p1_payoffs)):
            return False
        if not any(i in p2_support for i, x in enumerate(p2_payoffs) if x == max(p2_payoffs)):
            return False

        return True

    def _Hrepresentation(self, m1, m2):
        r"""
        Creates the H-representation strings required to use lrs nash.

        EXAMPLES::

            sage: A = matrix([[1, 2], [3, 4]])
            sage: B = matrix([[3, 3], [1, 4]])
            sage: C = NormalFormGame([A, B])
            sage: print C._Hrepresentation(A, B)[0]
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
            sage: print C._Hrepresentation(A, B)[1]
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
